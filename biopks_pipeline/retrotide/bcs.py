# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 21:25:37 2022

@author: Tyler Backman, Vincent Blay
"""

import cobra
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import OrderedDict
from copy import copy
from json import dumps
import pkg_resources


# %% PREPARE STARTERS AND EXTENDERS

def set_starters_extenders(path_starters=None, path_extenders=None):
    """
    Function loading the lists of starters and extenders used in RetroTide.
    List of starters and extenders are included by default with the package
    and are loaded whenever `import retrotide` is called.

    :param path_starters: String specifying the path to a text file containing
    the list of starters to be used.
    :type type_fp: str, optional
    :param path_extenders: String specifying the path to a text file containing
    the list of extenders to be used.
    :type type_fp: str, optional

    .. note::
        `starters` and `extenders` are reserved variables when using RetroTide,
        so they should not be manually modified.

    """
    # DATABASE OF PKS STARTERS
    if path_starters is None:
        path_starters = 'data/starters.smi'  # always use slash
    filepath = pkg_resources.resource_filename(__name__, path_starters)
    suppl = Chem.rdmolfiles.SmilesMolSupplier(
        filepath,
        delimiter='\t',
        titleLine=True,
        sanitize=True,
    )

    # now let's process these into a dict, and 'virtually attach to the ACP' by removal of the CoA
    global starters
    starters = {}

    for m in suppl:
        Chem.SanitizeMol(m)
        starters[m.GetProp('_Name')] = m

    # DATABASE OF PKS EXTENDER SUBSTRATES
    if path_extenders is None:
        path_extenders = 'data/extenders.smi'  # always use slash
    filepath = pkg_resources.resource_filename(__name__, path_extenders)
    extenderSuppl = Chem.rdmolfiles.SmilesMolSupplier(
        filepath,
        delimiter='\t',
        titleLine=True,
        sanitize=True,
    )

    # now let's process these into a dict, and 'virtually attach to the ACP' by removal of the CoA
    global extenders
    extenders = {}

    for m in extenderSuppl:
        rxn = AllChem.ReactionFromSmarts(
            # This reaction removes the CoA
            '[O:1]=[C:2]([O:3])[C:4][C:5](=[O:6])[S:7]'
            '>>'
            '[*:10]-[C:2](=[O:1])[C:4][C:5](=[O:6])[S:7].[O:3]')
        prod = rxn.RunReactants([m])[0][0]
        Chem.SanitizeMol(prod)
        extenders[m.GetProp('_Name')] = prod

    return


# Loads the defaul starters and extenders lists.
set_starters_extenders()


# Users can overwrite the default list of starters and extenders
# by calling this function with adequate inputs.


# %% DEFINE OBJECTS FOR REPRESENTING PKS

class Cluster:
    # Class representing a PKS design, which is just a list of modules

    def __init__(self, modules=None):
        if modules:
            self.modules = modules
        else:
            self.modules = []

    def computeProduct(self, structureDB, chain=False):
        '''
        This function computes the chemical product of this PKS design.
        If a mol object is passed as 'chain' only the final module operation is performed
        on this chain, and returned. This final module feature is to accelerate retrobiosynthesis.
        '''
        if chain:
            prod = chain
            modulesToExecute = [self.modules[-1]]  # last module only
        else:
            prod = False
            modulesToExecute = self.modules

        for module in modulesToExecute:
            if TE in module.domains:
                return module.domains[TE].operation(prod)

            if prod:
                moduleStructure = structureDB[module]

                # perform condensation if this isn't in the starter
                rxn = AllChem.ReactionFromSmarts(('[#6:10][C:1](=[O:2])[S:3].'
                                                  '[*:4][C:5]~[C:6]>>'
                                                  '[#6:10][C:5]~[C:6]'
                                                  '.[*:4].[C:1](=[O:2])[S:3]'))

                prod = rxn.RunReactants((prod, moduleStructure))[0][0]
                Chem.SanitizeMol(prod)

            else:  # starter module
                # rxn = AllChem.ReactionFromSmarts(
                #    # KSQ operation - this reaction removes the * and replaces it with C
                #    '[#0:1][C:2][C:3][C:4](=[O:5])-[S:6]'
                #    '>>'
                #    '[C:3][C:4](=[O:5])-[S:6].[#0:1][C:2]')
                # assert len(moduleStructure.GetSubstructMatches(Chem.MolFromSmiles('CC(=O)S'),
                #     useChirality=True)) == 1, Chem.MolToSmiles(moduleStructure)
                # prod = rxn.RunReactants([moduleStructure])[0][0]
                # prod = moduleStructure
                prod = starters[module.domains[AT].substrate]

        return prod


class Module:
    # Class representing a PKS module
    # self.domains is an OrderedDict where keys are domain classes,
    # and values are domain objects

    def __init__(self, product='', iterations=1, domains=None, loading=False):
        self.product = product
        self.iterations = iterations
        self.loading = loading
        if domains:
            self.domains = domains
        else:
            self.domains = OrderedDict()

    @staticmethod
    def domainTypes():
        '''
        Returns all domain types that can occur in this PKS in the catalytic order
        in which they operate.
        '''
        # return Domain.__subclasses__()
        return [AT, KR, DH, ER, TE]

    def computeProduct(self):
        '''
        computes the chemical product of this module
        '''
        chain = False
        for domaintype, domain in self.domains.items():
            chain = domain.operation(chain)

        return chain

    def __repr__(self):
        return repr(
            [cls.__name__ + repr(domain) for cls, domain in self.domains.items()] + ['loading: ' + repr(self.loading)])

    def __hash__(self):
        # produce a unique hash key for each domain configuration
        return hash(tuple(self.domains.values()) + (self.loading,))

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return hash(self) == hash(other)
        else:
            return False

    def __ne__(self, other):
        return (not self.__eq__(other))


class Domain:
    # Abstract base class used to build PKS catalytic domains

    def __init__(self, active):
        '''
        Initiates a new domain with a design as reported by designSpace
        '''

        self.active = active

    def design(self):
        '''
        Reports the design of this object
        '''

        return vars(self)

    @classmethod
    def designSpace(cls, module=False):
        '''
        Returns a set of objects representing the full design space
        of this domain.
        Can optionally take a PKS module to report only the compatible
        configurations of this domain with that design. Domains of this type
        in the design are ignored. If incompatible domains are included in the
        design, it just returns an empty list.
        '''

        raise NotImplementedError

    def operation(self, chain):
        '''
        Executes this domains operation on top of an existing PKS chain as an
        RDKit mol object, and returns
        the chemical product, as well as a cobrapy reacton representing the
        stoichiometry.
        '''

        raise NotImplementedError

    def reactants(self):
        '''
        Returns all reactants of this domain, excluding the substrate (polyketide chain).
        The format is a list of Cobrapy metabolites.
        '''

        raise NotImplementedError

    def products(self):
        '''
        Returns all products of this domain, excluding the polyketide chain.
        The format is a list of Cobrapy metabolites.
        '''

        raise NotImplementedError

    def __repr__(self):
        '''
        Returns a string representing this domain type for text based storage of
        PKS designs, or reporting to the user.
        Only prints activity if active=False to keep things concise.
        '''
        if self.active:
            designCopy = copy(self.design())
            del designCopy['active']
            return (repr(designCopy))
        else:
            return (repr(self.design()))

    def __hash__(self):
        # produce a unique hash key for each domain configuration
        return hash(dumps(self.design(), sort_keys=True))

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return hash(self) == hash(other)
        else:
            return False

    def __ne__(self, other):
        return (not self.__eq__(other))

    # Acyltransferase - selects and transfers Acyl-CoA units to the growing polyketide chain


class AT(Domain):

    def __init__(self, active, substrate):
        '''
        Initiates a new domain with a design as reported by designSpace
        '''
        self.active = active
        self.substrate = substrate

    @classmethod
    def designSpace(cls, module=False):
        if module:
            if module.loading:
                return [cls(active=True, substrate=s) for s in starters.keys()]

        # return only extension ATs unless passed a loading module for context
        return [cls(active=True, substrate=s) for s in extenders.keys()]

    def operation(self, chain, loading=False):
        if not chain:
            if loading:
                return starters[self.substrate]
            else:
                return extenders[self.substrate]
        else:
            # ATs here don't perform condensation, so need to operate first
            # the condensation is performed afterwards
            raise NotImplementedError

    def reactants(self):
        '''
        Returns all reactants of this domain, excluding the substrate (polyketide chain).
        The format is a list of Cobrapy metabolites.
        '''

        if self.substrate == 'Malonyl-CoA':
            substratecobrapy = cobra.Metabolite('malcoa_c', compartment='c')
        elif self.substrate == 'Methylmalonyl-CoA':
            substratecobrapy = cobra.Metabolite('mmcoa__S_c', compartment='c')
        else:
            substratecobrapy = cobra.Metabolite(self.substrate + '_c', compartment='c')

        return [
            cobra.Metabolite('h_c', compartment='c'),
            substratecobrapy
        ]

    def products(self):
        '''
        Returns all products of this domain, excluding the polyketide chain.
        The format is a list of Cobrapy metabolites.
        '''

        return [
            cobra.Metabolite('coa_c', compartment='c'),
            cobra.Metabolite('co2_c', compartment='c')
        ]

        raise NotImplementedError

    # Ketoreductase - reduces ketone groups to hydroxyl groups


# Ketoreductase - reduces ketone group to hydroxyl group
class KR(Domain):
    TYPE_CHOICES = {'B1', 'B', 'C1'}  # 2D change

    # TYPE_CHOICES = {'A1', 'A2', 'A', 'B1', 'B2', 'B', 'C1', 'C2'}
    # TYPE_CHOICES = {'A1', 'A2', 'A', 'B1', 'B2', 'B', 'C1', 'C2', 'U'}

    def __init__(self, active, type):
        '''
        Initiates a new domain with a design as reported by designSpace
        '''
        assert type in self.TYPE_CHOICES
        self.active = active
        self.type = type

    @classmethod
    def designSpace(cls, module=False):
        updatedTypeChoices = copy(cls.TYPE_CHOICES)

        if module and module.domains[AT].substrate != 'Malonyl-CoA':
            # if the domain occurs in a module WITHOUT a Malonyl-CoA AT, remove the A/B type
            updatedTypeChoices.difference_update({'A', 'B'})
        elif module:
            # if the domain occurs in a module WITH a Malonyl-CoA AT, keep only the A/B type
            updatedTypeChoices.difference_update({'A1', 'A2', 'B1', 'B2', 'C1', 'C2'})

        return [cls(active=True, type=type) for type in updatedTypeChoices] + [cls(active=False, type='B1')]

    def operation(self, chain):
        if self.type == 'A1':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C@:2]([O:3])[C@:4]'
                                              '[C:5](=[O:6])[S:7]'))
        elif self.type == 'A2':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C@:2]([O:3])[C@@:4]'
                                              '[C:5](=[O:6])[S:7]'))
        elif self.type == 'A':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C@:2]([O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B1':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C@@:2]([O:3])[C@:4]'
                                              '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B2':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C@@:2]([O:3])[C@@:4]'
                                              '[C:5](=[O:6])[S:7]'))
        elif self.type == 'B':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C@@:2]([O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]'))
        elif self.type == 'C1':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C:2](=[O:3])[C@:4]'
                                              '[C:5](=[O:6])[S:7]'))
        elif self.type == 'C2':
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C:2](=[O:3])[C@@:4]'
                                              '[C:5](=[O:6])[S:7]'))
        else:
            # By first specifying some stereochemistry in the reactants
            # and then explicitly "losing" the stereochemistry in the products
            # we can forget the stereochemistry in our molecule
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2](=[O:3])[C@:4]'
                                              '[C:5](=[O:6])[S:7]>>'
                                              '[#0:1][C:2]([O:3])[C:4]'
                                              '[C:5](=[O:6])[S:7]'))

        assert len(chain.GetSubstructMatches(Chem.MolFromSmiles('C(=O)CC(=O)S'),
                                             useChirality=True)) == 1, Chem.MolToSmiles(chain)
        prod = rxn.RunReactants((chain,))[0][0]
        Chem.SanitizeMol(prod)
        return prod

    def reactants(self):
        '''
        Returns all reactants of this domain, excluding the substrate (polyketide chain).
        The format is a list of Cobrapy metabolites.

        Stoich: ketone_pks_product + NADPH + H+ -> hydroxyl_pks_product + NADP+
        '''

        return [
            cobra.Metabolite('nadph_c', compartment='c'),
            cobra.Metabolite('h_c', compartment='c'),
        ]

    def products(self):
        '''
        Returns all products of this domain, excluding the polyketide chain.
        The format is a list of Cobrapy metabolites.

        Stoich: ketone_pks_product + NADPH + H+ -> hydroxyl_pks_product + NADP+
        '''

        return [
            cobra.Metabolite('nadp_c', compartment='c'),
        ]

    # Dehydratase - remove water molecule from a beta-hydroxy group, resulting in double bond formation


class DH(Domain):

    @classmethod
    def designSpace(cls, module=False):
        # adding False as a design type specifies that this domain is optional,
        # e.g. a PKS can exist without it
        if not module:
            return [cls(active=True), cls(active=False)]

        if not KR in module.domains:
            return [cls(active=False)]

        # require that we have an active B/B1 KR type
        if module.domains[KR].active and (module.domains[KR].type in {'B', 'B1', 'U'}):
            return [cls(active=True), cls(active=False)]
        else:
            return [cls(active=False)]

    def operation(self, chain):
        # try setting CH unchanged
        rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2]([O:3])[C:4][C:6](=[O:7])[S:8]>>'
                                          '[#0:1][CH1:2]=[CH1:4][C:6](=[O:7])[S:8].[O:3]'))
        # '[#0:1][CH1:2]=[CH1:4][C:6](=[O:7])[S:8].[O:3]'))
        assert len(chain.GetSubstructMatches(Chem.MolFromSmiles('C(O)CC(=O)S'),
                                             useChirality=True)) == 1, Chem.MolToSmiles(chain)
        prod = rxn.RunReactants((chain,))[0][0]
        try:
            Chem.SanitizeMol(prod)
        except ValueError:
            # if this has a methyl attached on the alpha carbon, we'll set CH0
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2]([O:3])[C:4][C:6](=[O:7])[S:8]>>'
                                              '[#0:1][CH1:2]=[CH0:4][C:6](=[O:7])[S:8].[O:3]'))
            assert len(chain.GetSubstructMatches(Chem.MolFromSmiles('C(O)CC(=O)S'),
                                                 useChirality=True)) == 1, Chem.MolToSmiles(chain)
            prod = rxn.RunReactants((chain,))[0][0]
            Chem.SanitizeMol(prod)
        return prod

    def reactants(self):
        '''
        Returns all reactants of this domain, excluding the substrate (polyketide chain).
        The format is a list of Cobrapy metabolites.

        Stoich: hydroxyl_pks_product -> alkene_pks_product + H2O
        '''

        return []

    def products(self):
        '''
        Returns all products of this domain, excluding the polyketide chain.
        The format is a list of Cobrapy metabolites.
        '''

        return [
            cobra.Metabolite('h2o_c', compartment='c'),
        ]

    # Enoylreductase - reduces the double bond introduced by a dehydratase, forming a saturated bond


class ER(Domain):

    @classmethod
    def designSpace(cls, module=False):
        # adding False as a design type specifies that this domain is optional,
        # e.g. a PKS can exist without it
        if not module:
            return [cls(active=True), cls(active=False)]

        if not DH in module.domains:
            return [cls(active=False)]

        # require that we have an active DH type
        if module.domains[DH].active:
            return [cls(active=True), cls(active=False)]
        else:
            return [cls(active=False)]

    def operation(self, chain):
        # try setting CH unchanged
        rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2]=[C:3][C:4](=[O:5])[S:6]>>'
                                          '[#0:1][CH2:2][CH2:3][C:4](=[O:5])[S:6]'))
        # '[#0:1][CH2:2][CH2:3][C:4](=[O:5])[S:6]'))
        assert len(chain.GetSubstructMatches(Chem.MolFromSmiles('C=CC(=O)S'),
                                             useChirality=True)) == 1, Chem.MolToSmiles(chain)
        prod = rxn.RunReactants((chain,))[0][0]
        try:
            Chem.SanitizeMol(prod)
        except ValueError:
            rxn = AllChem.ReactionFromSmarts(('[#0:1][C:2]=[C:3][C:4](=[O:5])[S:6]>>'
                                              '[#0:1][CH2:2][C@@H1:3][C:4](=[O:5])[S:6]'))
            # '[#0:1][CH2:2][CH:3][C:4](=[O:5])[S:6]'))
            assert len(chain.GetSubstructMatches(Chem.MolFromSmiles('C=CC(=O)S'),
                                                 useChirality=True)) == 1, Chem.MolToSmiles(chain)
            prod = rxn.RunReactants((chain,))[0][0]
        return prod

    def reactants(self):
        '''
        Returns all reactants of this domain, excluding the substrate (polyketide chain).
        The format is a list of Cobrapy metabolites.

        Stoich: alkene_pks_product + NADPH + H+ -> alkane_pks_product + NADP+
        '''

        return [
            cobra.Metabolite('nadph_c', compartment='c'),
            cobra.Metabolite('h_c', compartment='c')
        ]

    def products(self):
        '''
        Returns all products of this domain, excluding the polyketide chain.
        The format is a list of Cobrapy metabolites.
        '''

        return [
            cobra.Metabolite('nadp_c', compartment='c')
        ]

    # Thioesterase - releases the final polyketide product from the PKS enzymes


class TE(Domain):

    def __init__(self, active, cyclic, ring):
        '''
        Initiates a new domain with a design as reported by designSpace
        '''
        self.active = active
        self.cyclic = cyclic
        self.ring = ring

    @classmethod
    def designSpace(cls, module=False):
        # adding False as a design type specifies that this domain is optional,
        # e.g. a PKS can exist without it

        # For now this returns false so it doesn't get included in designs
        # later we need to deal with terminal domains in a better way
        return [cls(active=False, cyclic=False, ring=0)]

    def operation(self, chain):
        assert len(chain.GetSubstructMatches(Chem.MolFromSmiles('C(=O)S'),
                                             useChirality=True)) == 1, Chem.MolToSmiles(chain)

        index = -1
        if self.cyclic:
            rxn = AllChem.ReactionFromSmarts('([C:1](=[O:2])[S:3].[O,N:4][C:5][C:6])>>'
                                             '[C:1](=[O:2])[*:4][C:5][C:6].[S:3]')
            index -= self.ring
        else:
            rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[S:3]>>[C:1](=[O:2])[O].[S:3]')

        # Using index -1 will yield the largest ring
        prod = rxn.RunReactants((chain,))[index][0]
        Chem.SanitizeMol(prod)

        return prod

    def reactants(self):
        '''
        Returns all reactants of this domain, excluding the substrate (polyketide chain).
        The format is a list of Cobrapy metabolites.

        Stoich: alkane_pks_product + H2O -> free_product + H+
        # NOTE: this is for non-cyclic version, we should add support for the cyclic version
        '''

        return [cobra.Metabolite('h2o_c', compartment='c')]

    def products(self):
        '''
        Returns all products of this domain, excluding the polyketide chain.
        The format is a list of Cobrapy metabolites.
        '''

        return [cobra.Metabolite('h_c', compartment='c')]