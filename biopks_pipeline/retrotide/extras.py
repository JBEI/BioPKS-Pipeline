# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 21:38:16 2022

@author: Vincent Blay
"""

from .bcs import Module, Cluster, TE, AT
from itertools import product
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from queue import PriorityQueue
import pkg_resources
import re

from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')


# %% PREPARE LISTS OF PARTS

def _get_allModuleTypes_allStarterTypes():
    """
    Function loading the lists of starter and extender Modules available.
    These lists are used by the random PK library evolver and by Retrotide.
    This function is run automatically when the library is imported.

    .. note::
        `allStarterTypes` and `allModuleTypes` are reserved variables when
        using RetroTide.

    """

    # create DB of possible extension designs

    # get a list of lists showing each possible domain design
    domainTypes = Module.domainTypes()
    allDomainDesigns = [domaintype.designSpace() for domaintype in domainTypes]

    # the cartesian product of all designs specifies the full design space
    fullDesignSpace = product(*allDomainDesigns)

    # create a set of all possible modules
    designs = set()
    for design in fullDesignSpace:
        newModule = Module()
        for domain in design:
            if not domain.active:
                continue
            newModule.domains[type(domain)] = domain

        # check against every domain to confirm that this is a valid design
        # if a domain reports no valid configurations, this design is skipped
        try:
            for domain in newModule.domains.values():
                assert domain in domain.designSpace(newModule)
        except AssertionError:
            continue

        designs.add(newModule)

    global allModuleTypes
    allModuleTypes = list(designs)

    # Create DB of possible Starters

    global allStarterTypes
    allStarterTypes = []
    for mol1 in AT.designSpace(module=Module(loading=True)):
        starter = Module(domains={AT: mol1}, loading=True)
        allStarterTypes.append(starter)


_get_allModuleTypes_allStarterTypes()


def _get_structureDB():
    structureDB = {}
    mols = []
    legends = []
    for module in allModuleTypes:
        try:
            structureDB[module] = module.computeProduct()
            mols.append(structureDB[module])
            legends.append(repr(module))
        except IndexError:
            continue
    return structureDB


structureDB = _get_structureDB()


# %% PK library maker

def random_PK_library_maker(n_gen, max_modules=8, return_modules=False,
                            p=2, random_state=None, remove_duplicates=True):
    """
    Function generating a random library of polyketide molecules. It does this
    by randomly combining PKS modules. Returns the molecules in SMILES format.

    :param n_gen: Target number of polyketides to generate. The number actually
        obtained may be smaller after removing duplicates.
    :type n_gen: int
    :param max_modules: Maximum number of modules when simulating the synthases
        that produce the polyketides. Defaults to 8.
    :type max_modules: int, optional
    :param return_modules: If true, the output will include a list of strings
        specifying the PKS that produces each molecule. Defaults to False.
    :type max_modules: bool, optional
    :param p: Specifies the probability distribution used to sample the
        number of modules used to simulate each PKS (up to `max_modules`).
        Defaults to 2 (2nd order polynomial).
    :type p: str or int, optional
    :param random_state: Seed to ensure repeatable results.
    :type random_state: int, optional
    :param remove_duplicates: If True (default), ensures that all output SMILES are
        unique.
    :type remove_duplicates: bool, optional

    :return: Returns a list containing the polyketide molecules generated at
        random in SMILES format.
        If `return_modules` is set to True, the return will be a tuple with two
        lists. The first list contains the SMILES and
        the second list contains the strings representing the PKS that generate
        each of the polyketides.
    :rtype: list or tuple

    """
    if random_state is not None:
        np.random.seed(random_state)

    modules_lib = []
    smiles_lib = []

    # Let us decide the lengths of PKS that we will generate
    lengths, max_modules = __lengths_generator(max_modules, n_gen, p)

    for i in range(n_gen):
        if ((i + 1) % 5000) == 0:
            print(i + 1)

        starter = np.random.choice(allStarterTypes, 1)[0]
        modules = [starter]

        # We next randomly append N extender modules
        N = lengths[i]
        extender_modules = np.random.choice(allModuleTypes, size=N, replace=True)

        modules.extend(extender_modules)

        smiles = __modules_to_smiles(modules)

        if smiles is not None:
            modules_lib.append(modules)
            smiles_lib.append(smiles)

    if remove_duplicates:
        smiles_lib, idx = np.unique(smiles_lib, return_index=True)
        smiles_lib = list(smiles_lib)
        modules_lib = [modules_lib[i] for i in idx]

    print(f'{len(smiles_lib)} unique polyketides generated.')

    if return_modules:
        return smiles_lib, modules_lib
    else:
        return smiles_lib


def __modules_to_smiles(modules, remove_duplicates=True):
    # Let us compute the PK products that the PKS will produce
    try:
        cluster = Cluster(modules)
        attached_product = __cluster_to_attached_product(cluster)
        product = __attached_product_to_product(attached_product)
        smiles = __product_to_smiles(product)
    except:
        print(f'Skipping invalid set of modules: {modules}.')
        smiles = None
    return smiles


def __cluster_to_attached_product(cluster):
    attached_product = cluster.computeProduct(structureDB)
    return attached_product


def __attached_product_to_product(attached_product, fast=True):
    cyclic = np.random.choice([True, False], size=1)[0]
    ring = 0
    if cyclic:

        validring = []

        if fast:
            rango = [0]
        else:
            rango = range(0, 12)

        for ring in rango:
            try:
                debsTE = TE(active=True, cyclic=True, ring=ring)
                debsProduct = debsTE.operation(attached_product)
                validring.append(ring)
            except:
                pass

        if len(validring) == 0:
            cyclic = False
        else:
            ring = np.random.choice(validring, size=1)[0]

    debsTE = TE(active=True, cyclic=cyclic, ring=ring)

    debsProduct = debsTE.operation(attached_product)
    return debsProduct


def __product_to_smiles(product):
    return Chem.MolToSmiles(product)


def __lengths_generator(max_len, n_gen, p, lengths=None):
    """Returns a list of n_gen molecule lengths up to max_len.
    If max_len is 0, then max_len will be obtained from p or lengths.
    p describes the probability of generating different length values."""
    if max_len <= 0:
        if isinstance(p, np.ndarray):
            max_len = len(p)
        else:
            max_len = max(lengths)
        print(f'max_len set to {max_len}.')

    if isinstance(p, (int, float)):
        p = np.power(range(max_len), p)
    elif isinstance(p, str):
        p = p.lower()
        if p == 'exp':
            p = np.exp(range(max_len))
        elif p == 'empirical':
            p = np.zeros(max_len)
            for i in lengths:
                p[i - 1] += 1
        elif p == 'cumsum':
            # This makes drawing longer molecules at least as likely as shorter ones
            p = np.zeros(max_len)
            for i in lengths:
                p[i - 1] += 1
            p = p.cumsum()

        else:
            raise IOError("Invalid p input value: {p}.")

    p = p / sum(p)

    lengths = np.random.choice(range(1, max_len + 1), size=n_gen,
                               replace=True, p=p)
    return lengths, max_len


# %% PKSABLE FRAGMENT-BASED CLASSIFIER

def smiles_clean(smiles, return_idx=False):
    smiles = list(smiles)
    ind = []
    ind_bad = []
    clean_smiles = []
    for i, smi in enumerate(smiles):
        try:
            smi = max(smi.split("."), key=len)
            smi = smi.replace("/C", "C")
            smi = smi.replace("\\C", "C")
            smi = smi.replace("/c", "c")
            smi = smi.replace("\\c", "c")
            m = Chem.MolFromSmiles(smi, sanitize=True)
            if m is not None:
                ind.append(i)
                clean_smiles.append(smi)
            else:
                ind_bad.append(i)
                print(f"Warning: invalid SMILES in position {i}: {smiles[i]}")
        except:
            ind_bad.append(i)
            print(f"Warning: invalid SMILES in position {i}: {smiles[i]}")

    if return_idx:
        return clean_smiles, ind, ind_bad
    else:
        return clean_smiles


def remove_chirality(smi):
    for m in re.findall("\[C@.*?\]", smi):
        smi = smi.replace(m, "C")
    return smi


def _make_refkeys_file(positive_train_smiles, filename=None, radius=2):
    """
    Generates a file containing a list of keys used by `PKSable_fragment_classifier`.

    """

    refkeys = set()
    for smi in positive_train_smiles:
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprint(mol, radius=radius, useChirality=True)
        newpos = list(fp.GetNonzeroElements().keys())
        refkeys.update(newpos)

        smi_nochirality = remove_chirality(smi)
        mol = Chem.MolFromSmiles(smi_nochirality)
        fp = AllChem.GetMorganFingerprint(mol, radius=radius, useChirality=True)
        newpos = list(fp.GetNonzeroElements().keys())
        refkeys.update(newpos)

    if filename is None:
        filename = f'./data/refkeys_2_{len(positive_train_smiles)}_dual.pickle'
    with open(filename, 'w') as f:
        for item in refkeys:
            f.write("%s\n" % item)


def _get_refkeys_file(path=None):
    if path is None:
        path = 'data/refkeys_2_5199932_dual.txt'  # always use slash
    filepath = pkg_resources.resource_filename(__name__, path)
    with open(filepath, "r") as f:
        refkeys = f.read().splitlines()
    refkeys = set(map(int, refkeys))
    return refkeys


def PKSable_fragment_classifier(qsmiles, path_refkeys=None):
    """
    Estimates whether input molecules are PKSable (values
    close to 1) or not PKSable (values close to 0) according to the chemical
    space accessible to Retrotide. Predictions are made by
    evaluating what fraction of on fingeprint keys in the input molecule
    have also been observed in a reference set of PKSable molecules. A
    reference set of keys is used by default, which was obtained
    from a large set of molecules generated using `random_PK_library_maker`.

    :param qsmiles: List of query molecules in SMILES format.
    :type qsmiles: list
    :param path_refkeys: Allows specifying a custom list of keys to be used
        in assessing PKSability. If not provided, will use the list by default.
    :type path_refkeys: str, optional

    :return: Returns a numpy array of floats between 0 and 1 corresponding to
        the PKSability estimate for each input molecule.
    :rtype: np.array

    """

    refkeys = _get_refkeys_file(path_refkeys)

    radius = 2
    qsmiles = list(qsmiles)

    qsmiles, idx_good, idx_bad = smiles_clean(qsmiles, return_idx=True)

    i = 0
    probs = np.zeros(len(qsmiles))

    for smi in qsmiles:
        if ((i + 1) % 5000) == 0:
            print(i + 1)

        qmol = Chem.MolFromSmiles(smi)

        fp = AllChem.GetMorganFingerprint(qmol, radius=radius, useChirality=True)
        keys = set(fp.GetNonzeroElements().keys())
        intersect = keys.intersection(refkeys)
        prob = len(intersect) / len(keys)
        probs[i] = prob

        i += 1

    # let us insert nan's for any invalid SMILES
    ix = [idx_bad[i] - i for i in range(len(idx_bad))]
    probs = np.insert(probs, ix, np.nan)

    return probs


# %% POLYKETIDE EVOLVER

def PKS_library_evolver(
        mcw,
        model,
        spec,
        k1=800,
        k2=80,
        n_rounds=6,
        n_hits=10,
        max_modules=8,
        noise_factor=0.2,
        random_state=None,
        **kwargs):
    """
    Generates polyketide molecules predicted to meet a desired property specification.

    :param mcw: Embedder that takes a SMILES string as input and produces
        an input to model.
    :type mcw: function or MACAW object
    :param model: Function that predicts a property of interest for the `mcw`
        encoding of a molecule.
    :type model: function
    :param spec: Desired property specification value (desired output of `model`).
    :type spec: float
    :param k1: Number of polyketides to generate in each diversification
        step of the algorithm. Defaults to 800.
    :type k1: int, optional
    :param k2: Number of polyketides to choose in each selection step
        of the algorithm. Defaults to 80.
    :type k2: int, optional
    :param n_rounds: Number of rounds for the algorithm. Defaults to 6.
    :type n_rounds: int, optional
    :param n_hits: Number of molecules to return in the output. Defaults to 10.
    :type n_hits: int, optional
    :param max_modules: Maximum number of modules composing a PKS used in
        the first round. Defaults to 8.
    :type max_modules: int, optional
    :param noise_factor: Controls the level of uniform noise in generating the
        PKS. Should be set between 0 (no noise) and 1 (completely random). Defaults
        to 0.2.
    :type noise_factor: float, optional
    :param random_state: Seed to ensure repeatable results.
    :type random_state: int, optional

    :return: Returns a tuple containing three lists. The first list contains the
        SMILES of the hit molecules generated. The second list contains a numpy
        array with the predicted property values for each of the molecules. The
        third list contains the strings representing the PKS that make the molecules.
    :rtype: tuple


    """

    if random_state is not None:
        np.random.seed(random_state)

    # Generate PKS molecules
    smiles_lib, modules_lib = random_PK_library_maker(k1, max_modules=max_modules, return_modules=True)

    # Evaluate properties

    if not callable(mcw):
        try:
            mcw = mcw.transform
        except AttributeError:
            raise IOError("mcw input is not callable.")

    if not callable(model):
        try:
            model = model.predict
        except AttributeError:
            raise IOError("model input is not callable.")

    X = mcw(smiles_lib)
    Y = model(X)

    # Select best k2 subset
    idx = find_Knearest_idx(spec, Y, k=k2)
    smiles_lib_old = [smiles_lib[i] for i in idx]
    modules_lib_old = [modules_lib[i] for i in idx]
    Y_old = Y[idx]

    # Start FOR loop
    for k in range(n_rounds):
        print(f'\nRound {k + 1}\n')
        # Diversify from the subset
        print(f'len(smiles_lib_old)={len(smiles_lib_old)}')
        print(f'len(modules_lib_old)={len(modules_lib_old)}')
        smiles_lib, modules_lib = biased_PK_library_maker(modules_lib_old, n_gen=k1, max_modules=0)

        # Evaluate properties
        X = mcw(smiles_lib)
        Y = model(X)

        # Append best old molecules and remove duplicates
        smiles_lib += smiles_lib_old  # concatenates lists
        smiles_lib, idx = np.unique(smiles_lib, return_index=True)
        smiles_lib = list(smiles_lib)

        Y = np.concatenate((Y, Y_old))
        Y = Y[idx]

        modules_lib += modules_lib_old
        modules_lib = [modules_lib[i] for i in idx]

        # Select best k2 subset
        idx = find_Knearest_idx(spec, Y, k=k2)
        smiles_lib_old = [smiles_lib[i] for i in idx]
        Y_old = Y[idx]
        modules_lib_old = [modules_lib[i] for i in idx]

    # Return subset
    idx = find_Knearest_idx(spec, Y_old, k=n_hits)
    smiles_lib = [smiles_lib_old[i] for i in idx]  # Access multiple elements of a list
    Y = Y_old[idx]
    modules_lib = [modules_lib_old[i] for i in idx]

    return smiles_lib, Y, modules_lib


def biased_PK_library_maker(modules_lib, n_gen, max_modules=0,
                            p=2, random_state=None, remove_duplicates=True):
    """
    This function generates a library of SMILES and PKS designs around an input list of
    PKS. It is used to create diversity in the `PKS_library_evolver`
    function.

    :param modules_lib: List of strings specifying PKS.
    :type modules_lib: list
    :param n_gen: Target number of PKS to generate. The actual number in the
        output may be smaller after removing duplicates.
    :type n_gen: int
    :param max_modules: Maximum number of modules composing a PKS. If 0 (default),
        the maximum number observed in the input `modules_list` is used.
    :type max_modules: int, optional
    :param p: Specifies the probability distribution used to sample the
        number of modules used to simulate each PKS (up to `max_modules`).
        Defaults to 2 (2nd order polynomial).
    :type p: str or int, optional
    :param random_state: Seed to ensure repeatable results.
    :type random_state: int, optional
    :param remove_duplicates: If True (default), ensures that all output SMILES
        are unique.
    :type remove_duplicates: bool, optional

    :return: Returns a tuple containing two lists. The first list contains the
        SMILES strings produced by the library of PKS generated. The second list
        contains the strings representing the corresponding PKS.
    :rtype: tuple

    """

    if random_state is not None:
        np.random.seed(random_state)

    # Let us retrieve the alphabet of modules we will be using

    starter_modules = set()
    extender_modules = set()
    lengths = []

    for modules in modules_lib:
        starter_modules.add(modules[0])
        extender_modules.update(modules[1:])
        lengths.append(len(modules[1:]))

    starter_modules = list(starter_modules)
    extender_modules = list(extender_modules)

    # Let us decide the lengths of PKS that we will generate
    lengths, max_modules = __lengths_generator(max_modules, n_gen, p, lengths)

    # Maps between modules and indices
    starter_to_idx = dict(zip(starter_modules, range(len(starter_modules))))
    # idx_to_starter = dict(enumerate(starter_modules))

    emodule_to_idx = dict(zip(extender_modules, range(len(extender_modules))))
    # idx_to_emodule = dict(enumerate(extender_modules))

    # Let us build the probability matrix for each module
    prob_starter = np.zeros((1, len(starter_modules)))
    prob_emodule = np.zeros((max_modules, len(extender_modules)))

    for modules in modules_lib:
        for i in range(len(modules)):
            module = modules[i]
            if i == 0:
                idx = starter_to_idx[module]
                prob_starter[0, idx] += 1
            else:
                idx = emodule_to_idx[module]
                prob_emodule[i - 1, idx] += 1

    # Normalize the probability matrices
    prob_starter = __noise_adder(prob_starter)
    prob_emodule = __noise_adder(prob_emodule)

    # Let us sample from the probability matrix to generate new PKS designs
    modules_lib = []
    smiles_lib = []
    for i in range(n_gen):

        starter = np.random.choice(starter_modules, 1, p=prob_starter[0, :])[0]
        modules = [starter]

        N = lengths[i]
        for i in range(N):
            emodules = np.random.choice(extender_modules, size=1, p=prob_emodule[i, :])
            modules.extend(emodules)

        # Let us compute the PK products that the PKS will produce
        smiles = __modules_to_smiles(modules)

        if smiles is not None:
            modules_lib.append(modules)
            smiles_lib.append(smiles)

    if remove_duplicates:
        smiles_lib, idx = np.unique(smiles_lib, return_index=True)
        smiles_lib = list(smiles_lib)
        modules_lib = [modules_lib[i] for i in idx]

    print(f'{len(smiles_lib)} unique polyketides generated.')

    return smiles_lib, modules_lib


def __noise_adder(matrix, noise_factor=0.3):
    """Normalizes the input matrix row-wise and mixes it linearly with a
    uniform matrix"""
    # Here we will add some noise to the prob matrix and normalize it
    row_sums = matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # prevent division by zero
    prob_matrix = matrix / row_sums

    B = np.ones(prob_matrix.shape) / prob_matrix.shape[1]  # uniform matrix
    prob_matrix = (1 - noise_factor) * prob_matrix + noise_factor * B

    # normalize prob_matrix row-wise again, although should not be necessary
    row_sums = prob_matrix.sum(axis=1, keepdims=True)
    prob_matrix = prob_matrix / row_sums

    return prob_matrix


def find_Knearest_idx(x, arr, k=1):
    """
    Finds the `k` nearest values to number `x` in unsorted array `arr` using a
    heap data structue.

    Adapted from
    https://www.geeksforgeeks.org/find-k-closest-numbers-in-an-unsorted-array/

    """

    n = len(arr)
    k = min(n, k)
    # Make a max heap of difference with
    # first k elements.
    pq = PriorityQueue()

    idx = []
    for i in range(k):
        pq.put((-abs(arr[i] - x), i))

    # Now process remaining elements
    for i in range(k, n):
        diff = abs(arr[i] - x)
        p, pi = pq.get()
        curr = -p

        # If difference with current
        # element is more than root,
        # then put it back.
        if diff > curr:
            pq.put((-curr, pi))
            continue
        else:

            # Else remove root and insert
            pq.put((-diff, i))

    # Print contents of heap.
    while not pq.empty():
        p, q = pq.get()
        idx.append(q)

    idx = np.array(idx)[np.argsort(arr[idx])]  # sort idx by arr value
    return idx