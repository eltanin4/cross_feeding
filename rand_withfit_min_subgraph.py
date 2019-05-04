import numpy as np
from sat_marker import markSatMetsRxns
from prune_checks import isFitterSubset

def randWFitnessMinSubnet(satRxnVec, rxnMat, prodMat, sumRxnVec, 
                  coreProdRxns, coreTBP, nutrientSet, Currency,
                  stoich_matrix, Energy):   
    currSatRxnVec = np.copy(satRxnVec)
    mulFac = [1, 3, 3]

    while True:
        # Keeping track of what the graph currently looks like.
        currSatRxns = np.nonzero(currSatRxnVec)[0]

        # Marking out which reactions are singly removable.
        canRemoveVec = np.array([isFitterSubset(remRxn, currSatRxnVec, rxnMat, prodMat, sumRxnVec, 
                                                coreProdRxns, nutrientSet, Currency, coreTBP, stoich_matrix, Energy)
                                 for remRxn in currSatRxns]) * 1
        removableMets = currSatRxns[ np.where( canRemoveVec ) ]

        # Randomly permuting the singly removable reactions.
        removalOrder = np.random.permutation(removableMets)

        # Calling a vector of reaction that can be removed.
        for remRxn in removalOrder:
            tempSatRxnVec = np.copy(currSatRxnVec)
            tempSatRxnVec[remRxn] = 0
            
            # Calculating the marked set of reactions from the temporary set.
            tempSatMetVec, tempSatRxnVec = markSatMetsRxns(tempSatRxnVec, rxnMat, prodMat, sumRxnVec, 
                                                           coreProdRxns, nutrientSet, Currency)

            # Checking if this still produces the core molecule.
            if coreTBP in np.nonzero(tempSatMetVec)[0]:
                newCost= sum([sum(stoich_matrix[rxn][Energy] * mulFac) 
                              for rxn in np.nonzero(tempSatRxnVec)[0]])
                oldCost = sum([sum(stoich_matrix[rxn][Energy] * mulFac) 
                              for rxn in np.nonzero(currSatRxnVec)[0]])
                if newCost >= oldCost:
                    currSatRxnVec[remRxn] = 0
            else:
                isEverythingRemovable = False
                break

        if not removableMets.any():
            return np.nonzero(currSatRxnVec)[0]
