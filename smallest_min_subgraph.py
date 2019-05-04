import numpy as np
from prune_checks import isCoreProduced
import itertools

def smallestMinSubnets(satRxnVec, rxnMat, prodMat, sumRxnVec, coreProdRxns, 
                       coreTBP, nutrientSet, Currency):
    satRxns = np.nonzero(satRxnVec)[0]

    # Figuring out which single reactions can be removed.
    canRemoveVec = np.array([isCoreProduced(remRxn, satRxnVec, rxnMat, prodMat, sumRxnVec, 
                                            coreProdRxns, nutrientSet, Currency, coreTBP)
                             for remRxn in satRxns]) * 1
    
    # Filtering them out from the satisfied set of reactions.
    removableMets = satRxns[ np.where( canRemoveVec ) ]

    # Iterating over subsets of removable metabolites to find the largest coremovable set(s).
    for k in range(len(removableMets), 0, -1):
        removableSubsets, canRemoveSubset = [], []
        for removableSubset in itertools.combinations(removableMets, k):
            canRemoveSubset.append(isCoreProduced(list(removableSubset), satRxnVec, rxnMat, prodMat, sumRxnVec, 
                                                  coreProdRxns, nutrientSet, Currency, coreTBP))
            removableSubsets.append(removableSubset)

        # Checking if any subset is removable; if even one exists, it has to be the largest
        # possible coremovable one.
        canRemoveSubset = np.array(canRemoveSubset)
        if canRemoveSubset.any():
            removableSubsets = np.array(removableSubsets)
            return [set(satRxns).difference(set(remSet)) 
                    for remSet in removableSubsets[np.where(canRemoveSubset)]]
