import numpy as np
from prune_checks import isFitterSubset
import itertools

def fittestMinSubnets(satRxnVec, rxnMat, prodMat, sumRxnVec, coreProdRxns, 
                      coreTBP, nutrientSet, Currency, stoich_matrix, Energy):
    satRxns = np.nonzero(satRxnVec)[0]

    # Figuring out which single reactions can be removed.
    canRemoveVec = np.array([isFitterSubset(remRxn, satRxnVec, rxnMat, prodMat, sumRxnVec, 
                                            coreProdRxns, nutrientSet, Currency, coreTBP, stoich_matrix, Energy)
                             for remRxn in satRxns]) * 1
    
    # Filtering them out from the satisfied set of reactions.
    removableMets = satRxns[ np.where( canRemoveVec ) ]

    # Iterating over subsets of removable metabolites to find the largest coremovable set(s).
    for k in range(len(removableMets), 0, -1):
        removableSubsets = np.array(list(itertools.combinations(removableMets, k)))
        canRemoveSubset = np.array([isFitterSubset(remRxns, satRxnVec, rxnMat, prodMat, sumRxnVec,
                                                   coreProdRxns, nutrientSet, Currency, coreTBP, stoich_matrix, Energy)
                                    for remRxns in removableSubsets]) * 1

        # Checking if any subset is removable; if even one exists, it has to be the largest
        # possible coremovable one.
        if canRemoveSubset.any():
            return [set(satRxns).difference(set(remSet)) 
                    for remSet in removableSubsets[np.where(canRemoveSubset)]]
                    