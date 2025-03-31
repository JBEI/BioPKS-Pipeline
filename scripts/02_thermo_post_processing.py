"""
This post-processing script computes reaction thermodynamics of post-PKS pathways using eQuilibrator.
eQuilibrator requires a local compound cache which can be downloaded from zenodo.org/records/4128543
This compound cache allows eQuilibrator to lookup thermodynamic contributions of groups in existing compounds
ChemAxon is needed only to decompose new compounds into their constituent groups.
If reactions involve known compounds only that have been decomposed, their thermodynamics can still be calculated.
"""
import cvxpy
import sqlalchemy
import numpy as np
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_api import ComponentContribution, Q_
from equilibrator_cache.compound_cache import CompoundCache

import warnings
warnings.simplefilter('ignore')

# --------------------- Helper functions for calculating thermodynamics ---------------------
def parse_rxn(rxn_str: str):
    """
    Separates a reaction string into a list of reactants and of products
    Reaction str of the form substrate_smiles + cofactor_smiles = product_smiles + cofactor_smiles
    This will return [substrate, cofactor] and [product, cofactor]

    :param rxn_str: reaction string of the form above
    :return reactants_list (list): list of reactants involved, eg:[substrate, cofactor]
    :return products_list (list): list of products involved, eg: [product, cofactor]
    """
    reactants, products = rxn_str.split(" = ")
    reactants_list = reactants.split(" + ")
    products_list = products.split(" + ")
    return reactants_list, products_list

def rxn_constructor_in_accession_IDs(reactants_list: list, products_list: list, lc: any):
    """
    Construct a reaction string of the form substrate + cofactor = product + cofactor
    Except instead of SMILES, eQuilibrator accession IDs are used instead for all species
    This output can directly be fed into eQuilibrator to calculate reaction dG

    :param reactants_list: list of reactants on the LHS (includes cofactors)
    :param products_list: list of products on the RHS (includes cofactors)
    :param lc: eQuilibrator compounds' cache object
    :return rxn_str_in_db_accessions: reaction string in eQuilibrator accession IDs
    """

    # note that compound objects will be created even if compound cannot be decomposed

    # create eQuilibrator compound objects for all LHS species
    reactant_cpd_objects = lc.get_compounds(reactants_list)

    # create eQuilibrator compound objects for all RHS species
    product_cpd_objects = lc.get_compounds(products_list)

    # initialize empty string
    rxn_str_in_db_accessions = ""

    try:
        for reactant in reactant_cpd_objects:
            # get eQuilibrator accession IDs for all LHS species (if decomposable)
            # if compound is not decomposable into groups, code will error out
            reactant_accession = reactant.get_accession()
            rxn_str_in_db_accessions += reactant_accession
            rxn_str_in_db_accessions += " + "

        # remove extra ' + ' sign on the right
        rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(" + ")

        # add ' = ' sign before switching over to products side
        rxn_str_in_db_accessions += " = "

        for product in product_cpd_objects:
            # get eQuilibrator accession IDs for all RHS species (if decomposable)
            product_accession = product.get_accession()
            rxn_str_in_db_accessions += product_accession
            rxn_str_in_db_accessions += " + "

        # remove extra ' + ' sign on the right
        rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(" + ")

        # return reaction string in terms of eQuilibrator accession IDs
        return rxn_str_in_db_accessions

    except:
        # if even one species on either LHS or RHS cannot be decomposed, return None
        return None

def calc_dG_frm_rxn_str(new_rxn_str: str, pH: float, pMg: float, ionic_strength: str, temp: str, cc:any):
    rxn_object = cc.parse_reaction_formula(new_rxn_str)

    cc.p_h = Q_(pH)
    cc.p_mg = Q_(pMg)
    cc.ionic_strength = Q_(ionic_strength)
    cc.temperature = Q_(temp)

    phys_dG_value = float(
        str(cc.physiological_dg_prime(rxn_object).value).rstrip(" kilojoule/mole")
    )
    phys_dG_error = float(
        str(cc.physiological_dg_prime(rxn_object).error).rstrip(" kilojoule/mole")
    )
    std_dG_value = float(
        str(cc.standard_dg_prime(rxn_object).value).rstrip(" kilojoule/mole")
    )
    std_dG_error = float(
        str(cc.standard_dg_prime(rxn_object).error).rstrip(" kilojoule/mole")
    )

    return rxn_object, phys_dG_value, phys_dG_error, std_dG_value, std_dG_error

def pick_MDF_constraints(S: any, Nc: any, Nr: any, ln_conc: any, dg_prime: any, B: any, lb: float, ub: float):
    """
     Determines the appropriate set of constraints for calculating the Maximal Driving Force (MDF)
     based on the presence of specific cofactor pairs within a reaction system.

     This function checks for the presence of certain cofactor pairs (e.g., ATP/ADP, NADH/NAD+, etc.)
     in the reaction's stoichiometric matrix. Depending on the cofactors present, it applies specific
     constraints related to their concentration ratios. If none of the specific cofactor pairs are present,
     regular concentration and thermodynamic constraints are applied.

     Parameters:
     - S (any): Stoichiometric matrix of the reaction system.
     - Nc (any): Number of compounds in the system.
     - Nr (any): Number of reactions in the system.
     - ln_conc (any): Natural logarithm of the concentrations of compounds.
     - dg_prime (any): Standard Gibbs free energy changes of reactions.
     - B (any): Predefined biochemical parameter for constraints.
     - lb (float): Lower bound for the concentration of compounds.
     - ub (float): Upper bound for the concentration of compounds.

     Returns:
     - List: A list of constraints that should be applied to the reaction system to calculate the MDF.
             The specific constraints depend on which cofactor pairs are present in the system.

     The function supports multiple cofactor pairs, and it will return the constraints for the first
     matched pair found. If no specific cofactor pairs are found, it returns general regulatory constraints.
     """
    ### Cofactors to track for specific constraints

    # [ATP]/[ADP] = 10
    ATP_inchi_key = "ZKHQWZAMYRWXGA-UHFFFAOYSA-N" # checked
    ADP_inchi_key = "XTWYTFMLZFPYCI-UHFFFAOYSA-N" # checked

    # [ADP]/ [AMP] = 1
    AMP_inchi_key = "UDMBCSSLTHHNCD-UHFFFAOYSA-N" # checked

    # [NADH]/[NAD+] = 0.1
    NADH_inchi_key = "BOPGDPNILDQYTO-UHFFFAOYSA-N" # checked
    NAD_plus_inchi_key = "BAWFJGJZGIEFAR-UHFFFAOYSA-O" # checked

    # [NADPH]/[NADP+] = 10
    NADP_plus_inchi_key = "XJLXINKUBYWONI-UHFFFAOYSA-O" # checked
    NADPH_inchi_key = "ACFIXJIJDZMPPO-UHFFFAOYSA-N" # checked

    # Get all compounds' inchi keys from the stoichiometric matrix, S
    compound_inchi_keys_list = [compound.inchi_key for compound in list(S.index)]

    ATP_ADP_present = [ATP_inchi_key in compound_inchi_keys_list,
                         ADP_inchi_key in compound_inchi_keys_list]

    ADP_AMP_present = [ADP_inchi_key in compound_inchi_keys_list,
                           AMP_inchi_key in compound_inchi_keys_list]

    NADH_NADplus_present = [NADH_inchi_key in compound_inchi_keys_list,
                                NAD_plus_inchi_key in compound_inchi_keys_list]

    NADP_NADPplus_present = [NADPH_inchi_key in compound_inchi_keys_list,
                                 NADP_plus_inchi_key in compound_inchi_keys_list]

    # initialize a list of basic constraints
    constraints_list = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime]

    if all(ATP_ADP_present):
        for i, inchi_key in enumerate(compound_inchi_keys_list):
            if inchi_key == ATP_inchi_key:
                ATP_index = i
            if inchi_key == ADP_inchi_key:
                ADP_index = i

        constraints_list.append( ln_conc[ATP_index] - ln_conc[ADP_index] == 2.303 )

    if all(ADP_AMP_present):
        for i,inchi_key in enumerate(compound_inchi_keys_list):
            if inchi_key == ADP_inchi_key:
                ADP_index = i
            if inchi_key == AMP_inchi_key:
                AMP_index = i

        constraints_list.append( ln_conc[ADP_index] == ln_conc[AMP_index] )

    if all(NADH_NADplus_present):
        for i,inchi_key in enumerate(compound_inchi_keys_list):
            if NAD_plus_inchi_key == inchi_key:
                NAD_plus_index = i
            if NADH_inchi_key == inchi_key:
                NADH_index = i

        constraints_list.append( ln_conc[NAD_plus_index] - ln_conc[NADH_index] == 2.303 )

    if all(NADP_NADPplus_present):
        for i,inchi_key in enumerate(compound_inchi_keys_list):
            if inchi_key == NADP_plus_inchi_key:
                NADP_plus_index = i

            if inchi_key == NADPH_inchi_key:
                NADPH_index = i

        constraints_list.append( ln_conc[NADPH_index] - ln_conc[NADP_plus_index] == 2.303 )

    return constraints_list

def calc_MDF(rxn_object, lb: float, ub: float, cc:any):
    # Setting up optimization problem with variables
    standard_dgr_prime_mean, standard_dgr_Q = cc.standard_dg_prime_multi(
        rxn_object, uncertainty_representation="fullrank"
    )

    S = cc.create_stoichiometric_matrix_from_reaction_objects(rxn_object)
    Nc, Nr = S.shape

    ln_conc = cvxpy.Variable(shape=Nc, name="metabolite log concentration")  # vector
    B = cvxpy.Variable()  # scalar
    dg_prime = -(standard_dgr_prime_mean.m_as("kJ/mol") + cc.RT.m_as("kJ/mol") * S.values.T @ ln_conc)

    constraints = pick_MDF_constraints(S, Nc, Nr, ln_conc, dg_prime, B, lb, ub)

    prob_max = cvxpy.Problem(cvxpy.Maximize(B), constraints)
    prob_max.solve()
    max_df = prob_max.value
    return max_df

# -------------------------------------------------------------------------------------------

# read in filepath to SQlite compounds database (taken from Zenodo) for eQuilibrator
with open('../equilibrator_filepath.txt','r') as f:
    eq_compounds_db_filepath = f.read()

URI_EQ = eq_compounds_db_filepath

lc = LocalCompoundCache()
lc.ccache = CompoundCache(sqlalchemy.create_engine(f'sqlite:////{URI_EQ}'))
cc = ComponentContribution(ccache=lc.ccache)

target_name = '1_Amino_2_propanol'

# calculating thermodynamics of PKS reactions is currently not supported so sequence has to be ['pks','bio'
pathway_sequence = ['pks','bio']
num_bio_steps = 2

# enter parameters to calculate reaction thermodynamics
pH = 7.4
pMg = 3.0
ionic_strength = "0.25M"
temp = "298.15K"
lb = 1e-4 # 0.1 millimolar (lower bound for MDF calculation)
ub = 1e-1 # 100 millimolar (upper bound for MDF calculation)

# calculating thermodynamics is final step so enzymes have to be queried first using 01_find_enzymes_for_bio.py
input_filepath = f'../data/results_logs/{target_name}_PKS_BIO{num_bio_steps}_with_enzymes.txt'
output_filepath = f"../data/results_logs/{target_name}_PKS_BIO{num_bio_steps}_with_enzymes_and_thermo.txt"

rxns_count = 0

with open(input_filepath, 'r') as file:
    lines = file.readlines()
    for i,line in enumerate(lines):
        if ' = ' in line:
            rxns_count += 1
            rxn_index = i
            rxn = line.strip('\n').strip('    ')
            reactants_list, products_list = parse_rxn(rxn)
            rxn_str_in_eq_accessions = rxn_constructor_in_accession_IDs(reactants_list,products_list,lc)

            rxn_object, phys_dG_value, \
            phys_dG_error, std_dG_value, std_dG_error = calc_dG_frm_rxn_str(new_rxn_str = rxn_str_in_eq_accessions,
                                                                            pH = pH,
                                                                            pMg = pMg,
                                                                            ionic_strength = ionic_strength,
                                                                            temp = temp,
                                                                            cc = cc)

            # reaction object for a single reaction has to be passed in as a list
            MDF = calc_MDF(rxn_object = [rxn_object], lb = lb, ub = ub, cc = cc)

            lines[i] = f'    {rxn}, Rxn MDF: {MDF:.2f} kJ/mol, std dG value: {std_dG_value:.2f} kJ/mol, rxn std dG error: {std_dG_error:.2f} kJ/mol, rxn phys dG value: {phys_dG_value:.2f} kJ/mol, rxn phys dG error: {phys_dG_error:.2f} kJ/mol\n'

if rxns_count != 0:
    with open(output_filepath,'w') as new_file:
        new_file.writelines(lines)
else:
    print('No pathways and reactions found. As such, no reaction thermodynamic values were computed.')