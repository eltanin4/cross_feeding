import numpy as np
from sat_marker import markSatMetsRxns

def prunedSatsMets(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                   coreProdRxns, nutrientSet, Currency, coreTBP):
    # Creating a temporary set of reactions with
    # some reactions pruned.
    tempSatRxns = np.copy(satRxns)
    tempSatRxns[remRxns] = 0
    
    # Calculating the marked set of reactions from the temporary set.
    return markSatMetsRxns(tempSatRxns, rxnMat, prodMat, sumRxnVec, 
                           coreProdRxns, nutrientSet, Currency)

#-------------------------------------------------------------------------

def isCoreProduced(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                   coreProdRxns, nutrientSet, Currency, coreTBP):

    tempSatMets, tempSatRxns = prunedSatsMets(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                                              coreProdRxns, nutrientSet, Currency, coreTBP)
    
    # Checking if this still produces the core molecule.
    if coreTBP in np.nonzero(tempSatMets)[0]:
        return True
    return False

#-------------------------------------------------------------------------

def isFitterSubset(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                   coreProdRxns, nutrientSet, Currency, coreTBP, stoich_matrix, Energy):
    
    tempSatMets, tempSatRxns = prunedSatsMets(remRxns, satRxns, rxnMat, prodMat, sumRxnVec, 
                                              coreProdRxns, nutrientSet, Currency, coreTBP)

    mulFac = [1, 3, 3]
    # Checking if this still produces the core molecule, and 
    # generates at least as much ATP/NADH/NADPH.
    if coreTBP in np.nonzero(tempSatMets)[0]:
        newCost= sum([sum(stoich_matrix[rxn][Energy] * mulFac) 
                      for rxn in np.nonzero(tempSatRxns)[0]])
        oldCost = sum([sum(stoich_matrix[rxn][Energy] * mulFac) 
                      for rxn in np.nonzero(satRxns)[0]])
        if newCost >= oldCost:
            return True
    return False
