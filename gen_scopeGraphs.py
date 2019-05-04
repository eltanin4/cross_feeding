import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools

#-------------------------------------------------------------------------
# Generating the reaction graph.
#-------------------------------------------------------------------------
def genRxnGraph(scopeRxns, rxnMat, prodMat):
    # Initializing the reaction graph.
    rxnGraph = nx.DiGraph()
    # rxnGraph.add_nodes_from(['R' + str(sr) for sr in scopeRxns])

    # Generating the reaction pairs to get links over.
    rxnPairs = list(zip(* [ iter(scopeRxns) ] * 2))
    rxnPairs += [i[::-1] for i in rxnPairs]

    for currPair in rxnPairs:
        if np.logical_and(prodMat[currPair[0]], 
                         rxnMat[currPair[1]]).any():
            rxnGraph.add_edge('R' + str(currPair[0]), 'R' + str(currPair[1]))
        if np.logical_and(rxnMat[currPair[0]], 
                         prodMat[currPair[1]]).any():
            rxnGraph.add_edge('R' + str(currPair[1]), 'R' + str(currPair[0]))

    return rxnGraph

#-------------------------------------------------------------------------
# Generating the bipartite graph.
#-------------------------------------------------------------------------
def genScopeGraph(scopeRxns, rxnMat, prodMat):
    # Calculating the scopeMets.
    scopeRxnVec = np.zeros(len(rxnMat))
    scopeRxnVec[scopeRxns] = 1
    scopeMetVec = (np.dot(np.transpose(rxnMat), scopeRxnVec) > 0) * 1
    scopeMetVec = (np.dot(np.transpose(prodMat), scopeRxnVec) + scopeMetVec > 0) * 1
    scopeMets = np.nonzero(scopeMetVec)[0]

    # Trying to write a bipartite graph of the expanded scope.
    scopeGraph = nx.DiGraph()
    scopeGraph.add_nodes_from(['C' + str(sm) for sm in scopeMets], bipartite = 0)
    scopeGraph.add_nodes_from(['R' + str(sr) for sr in scopeRxns], bipartite = 1)

    # Connect all the metabolites to reactions.
    for rxnInd in scopeRxns:
        scopeGraph.add_edges_from([('C' + str(currReactant), 'R' + str(rxnInd)) 
                                   for currReactant in np.nonzero(rxnMat[rxnInd])[0]])
        scopeGraph.add_edges_from([('R' + str(rxnInd), 'C' + str(currProd)) 
                                   for currProd in np.nonzero(prodMat[rxnInd])[0]])

    return scopeGraph

#-------------------------------------------------------------------------
# Generating the metabolite graph.
#-------------------------------------------------------------------------
def genMetGraph(scopeRxnVec, rxnMat, prodMat):
    # Calculating which metabolites are in the scope of the reaction set.
    scopeMetVec = (np.dot(np.transpose(rxnMat), scopeRxnVec) > 0) * 1
    scopeMetVec = (np.dot(np.transpose(prodMat), scopeRxnVec) + scopeMetVec > 0) * 1
    scopeMets = np.nonzero(scopeMetVec)[0]

    metGraph = nx.DiGraph()
    for fromMet, toMet in itertools.product(scopeMets, scopeMets):
        if fromMet != toMet:
            if np.logical_and(
                              np.logical_and((rxnMat[:, fromMet] > 0) * 1, 
                                             (prodMat[:, toMet] > 0) * 1), 
                              scopeRxnVec).any():
                metGraph.add_edge('C' + str(fromMet), 'C' + str(toMet))
                
    return metGraph
