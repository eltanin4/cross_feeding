import numpy as np
from smallest_min_subgraph import smallestMinSubnets

def givePathSubgraph(rxnMat, prodMat, sumRxnVec, nutrientSet, Currency, coreTBP):
    """
    Takes in the stoichiometric matrix in the form of a reaction 
    and product matrix, along with a sum vector mentioning the 
    number of reactants in each reaction in the matrix. Also takes
    in a set of nutrients available in the medium and a set of currency
    metabolites. Receives information about which core molecule must
    be produced.

    By reverse scope-expansion, discovers a subgraph that uses only the 
    nutrient set and currency metabolites to successfully generate the 
    core molecule.

    RETURNS:
    Returns the satisfied metabolites and reactions in the subgraph
    (corresponding to a pathway.)

    Returns empty vectors if no such wholly satisfied subgraph is found.

    satMets, satRxns are sets of metabolites and reactions, with their custom IDs.
    """

    # Initializing all the vectors to propagate the satisfied subgraph search.
    seedVec, rxnProc = np.zeros(len(np.transpose(rxnMat))), np.zeros(len(rxnMat))
    seedVec[coreTBP] = 1
    coreProdRxns = (prodMat[:, coreTBP] == 1) * 1
    currScopeMets, prevScopeMets, deltaMetVec = np.copy(seedVec), np.copy(seedVec), np.copy(seedVec)

    while True:
        prevScopeMets = np.logical_or(currScopeMets, prevScopeMets)

        # Propagating reverse scope.
        rxnProc = (np.dot(prodMat, deltaMetVec) + rxnProc > 0) * 1
        currScopeMets = (np.dot(np.transpose(rxnMat), 
                                np.dot(prodMat, deltaMetVec)) > 0) * 1
        currScopeMets = np.logical_xor(currScopeMets, 
                                       np.logical_and(currScopeMets, prevScopeMets))

        # Marking the satisfied metabolites and reactions.
        satMets, satRxns = markSatMetsRxns(rxnProc, rxnMat, prodMat, sumRxnVec, 
                                           coreProdRxns, nutrientSet, Currency)

        # If core has been reached, pruning the satisfied sets, then returning.
        if coreTBP in np.nonzero(satMets)[0]:
            satMets, satRxns = smallestMinSubnets(satMets, satRxns, rxnMat, prodMat, 
                                                  sumRxnVec, coreProdRxns, coreTBP, 
                                                  nutrientSet, Currency)
            return satMets, satRxns

        # Calculating the new metabolites that need to be produced
        # due to the new steps added to the reverse scope.
        deltaMetVec = np.logical_xor(currScopeMets, 
                                     np.logical_and(currScopeMets, satMets)) * 1

        # Quitting if fail condition is reached.
        if set(np.nonzero(currScopeMets)[0]).issubset(set(np.nonzero(prevScopeMets)[0])):
            return np.zeros(len(np.transpose(rxnMat))), np.zeros(len(rxnMat))
