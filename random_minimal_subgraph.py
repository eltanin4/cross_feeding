import numpy as np
from sat_marker import markSatMetsRxns
from prune_checks import isCoreProduced

def randMinSubnet(satRxnVec, rxnMat, prodMat, sumRxnVec, 
                  coreProdRxns, coreTBP, nutrientSet, Currency):
    currSatRxnVec = np.copy(satRxnVec)

    while True:
        # Keeping track of what the graph currently looks like.
        currSatRxns = np.nonzero(currSatRxnVec)[0]

        # Marking out which reactions are singly removable.
        canRemoveVec = np.array([isCoreProduced(remRxn, currSatRxnVec, rxnMat, prodMat, sumRxnVec, 
                                                coreProdRxns, nutrientSet, Currency, coreTBP)
                                 for remRxn in currSatRxns]) * 1

        # Keeping a list of what is removable.
        removableMets = currSatRxns[ np.where( canRemoveVec ) ]

        # If nothing can be removed, minimality condition satisfied. 
        # Quitting with what we have now.
        if not removableMets.any():
            return np.nonzero(currSatRxnVec)[0]

        # Randomly permuting the singly removable reactions.
        removalOrder = np.random.permutation(removableMets)

        # Calling a vector of reaction that can be removed.
        for remRxn in removalOrder:
            if isCoreProduced(remRxn, currSatRxnVec, rxnMat, prodMat, 
                              sumRxnVec, coreProdRxns, nutrientSet, 
                              Currency, coreTBP):
                currSatRxnVec[remRxn] = 0
            else:
                break
