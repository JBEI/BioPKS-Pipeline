# -*- coding: utf-8 -*-
"""
@author: Tyler Backman, Vincent Blay
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit import DataStructs
from rdkit.Chem.AtomPairs import Pairs
import itertools

from .extras import allStarterTypes, allModuleTypes, structureDB
from .bcs import Cluster
from .AtomAtomPathSimilarity import getpathintegers, AtomAtomPathSimilarity


def compareToTarget(structure, target, similarity='atompairs'):
    # convert to smiles and back to fix stereochem- this shouldn't be needed!
    # testProduct = Chem.MolFromSmiles(Chem.MolToSmiles(testProduct, isomericSmiles=True))

    # remove C(=O)S from testProduct before comparison
    testProduct = Chem.rdmolops.ReplaceSubstructs(structure, Chem.MolFromSmiles('C(=O)S'), Chem.MolFromSmiles('C'))[0]

    if similarity == 'mcs_with_stereo':
        # MCS
        # result=rdFMCS.FindMCS([target, testProduct], timeout=1, matchChiralTag=True, ringMatchesRingOnly=True) # search for 1 second max
        result = rdFMCS.FindMCS([target, testProduct], timeout=1, matchValences=True, matchChiralTag=True,
                                bondCompare=Chem.rdFMCS.BondCompare.CompareOrderExact)  # search for 1 second max

        if result.canceled:
            print('MCS timeout')
        score = result.numAtoms / (len(testProduct.GetAtoms()) + len(target.GetAtoms()) - result.numAtoms)

    elif similarity == 'mcs_without_stereo':
        result = rdFMCS.FindMCS([target, testProduct], timeout=1, matchValences=True, matchChiralTag=False,
                                bondCompare=Chem.rdFMCS.BondCompare.CompareOrderExact)  # search for 1 second max

        if result.canceled:
            print('MCS timeout')
        score = result.numAtoms / (len(testProduct.GetAtoms()) + len(target.GetAtoms()) - result.numAtoms)

    elif similarity == 'atompairs':
        # atom pair
        ms = [target, testProduct]
        pairFps = [Pairs.GetAtomPairFingerprint(x) for x in ms]
        score = DataStructs.TanimotoSimilarity(pairFps[0], pairFps[1])  # can also try DiceSimilarity

    elif similarity == 'atomatompath':
        score = AtomAtomPathSimilarity(target, testProduct, m1pathintegers=targetpathintegers)

    elif callable(similarity):
        score = similarity(target, testProduct)

    else:
        raise IOError('Invalid similarity input')

    return score

def designPKS(targetMol, previousDesigns=False, maxDesignsPerRound=25, similarity='atompairs'):
    # This function recursively designs PKSs

    if previousDesigns:
        print('computing module ' + str(len(previousDesigns)))
    else:
        print('computing module 1')
        if similarity == 'atomatompath':
            global targetpathintegers
            targetpathintegers = getpathintegers(targetMol)

    # get all starter and extender modules
    # allModuleTypes, allStarterTypes = _get_allModuleTypes_allStarterTypes() # these variables are global

    # if previousDesigns=False this is the loading module iteration so we will
    # use that as our starting point
    if not previousDesigns:
        previousDesigns = [
            [[Cluster(modules=[i]), 0.0, Cluster(modules=[i]).computeProduct(structureDB)] for i in allStarterTypes]]

    # create cartesian product of all previous designs with all possible new modules
    extendedSets = list(itertools.product([design for design in previousDesigns[-1]], allModuleTypes))

    # perform each extension
    designs = [Cluster(modules=x[0][0].modules + [x[1]]) for x in extendedSets]

    print('   testing ' + str(len(designs)) + ' designs')

    # compute structures
    prevStructures = [x[0][2] for x in extendedSets]
    structures = [design.computeProduct(structureDB, chain=prevStructure) for design, prevStructure in
                  zip(designs, prevStructures)]

    # compare modules to target
    scores = [compareToTarget(structure, targetMol, similarity) for structure in structures]

    # assemble scores
    assembledScores = list(zip(designs, scores, structures))

    # sort designs by score
    assembledScores.sort(reverse=True, key=lambda x: x[1])

    # find best score from previous design round
    bestPreviousScore = previousDesigns[-1][0][1]

    # get the score of the first (best) design from this round
    bestCurrentScore = assembledScores[0][1]

    print('   best score is ' + str(bestCurrentScore))

    if bestCurrentScore > bestPreviousScore:
        # run another round if the scores are still improving

        # keep just top designs for the next round
        if len(assembledScores) > maxDesignsPerRound:
            assembledScores = assembledScores[0:maxDesignsPerRound]

        # recursively call self
        return designPKS(targetMol, previousDesigns=previousDesigns + [assembledScores],
                         maxDesignsPerRound=maxDesignsPerRound, similarity=similarity)

    else:
        # if these designs are no better than before, just return the last round
        return previousDesigns
