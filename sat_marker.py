import numpy as np

def markSatMetsRxns(rxnProc, rxnMat, prodMat, sumRxnVec, 
                        coreProdRxns, nutrientSet, Currency):
    """
    Takes in a bunch of reactions, core producing reactions and 
    nutrients and using the KEGG provided chemistry, markes all
    reactions and metabolites in the list of reactions that are 
    'satisfied', i.e. which can be reached via simply a seed set 
    of the given nutrients and currency metabolites.

    RETURNS:

        satMetVec, satRxnVec: the current sets of metabolites and reactions
                          that are said to be satisfied.
    """
    satMetVec = np.zeros(len(np.transpose(rxnMat)))
    satMetVec[nutrientSet + Currency], satRxnVec = 1, np.zeros(len(rxnMat))

    while True:
        oldSatRxnVec = np.copy(satRxnVec)

        # Marking first reactions, then metabolites, iteratively.
        satRxnVec = np.logical_and((np.dot(rxnMat, satMetVec) - sumRxnVec == 0) * 1, 
                                  rxnProc) * 1
        satMetVec = (np.dot(np.transpose(prodMat), satRxnVec) + satMetVec > 0) * 1

        # Checking if all satisfied nodes have been marked.
        if np.array_equal(oldSatRxnVec, satRxnVec):
            break
    
    return satMetVec, satRxnVec
