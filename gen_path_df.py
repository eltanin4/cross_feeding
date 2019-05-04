import pandas as pd
from expand_seed import *
from gen_scopeGraphs import *
import numpy as np
from give_path_subgraph import givePathSubgraph
from org_generators import randOrg
from sec_byproducts import secByproducts

#-------------------------------------------------------------------------

def pathDF(rxnMat, prodMat, nutrientSet, Currency, sumRxnVec):
    pathDict, metDict = {}, {}
    for coreTBP in Core:
        satMetVec, satRxnVec = givePathSubgraph(rxnMat, prodMat, sumRxnVec, 
                                                nutrientSet, Currency, coreTBP)
        satMets, satRxns = np.nonzero(satMetVec)[0], np.nonzero(satRxnVec)[0]
        print(coreTBP, '\t', cpd_string_dict[id_to_kegg[coreTBP]], '\t', len(satRxns))
        if not satMetVec.any():
            print(coreTBP, '\t', cpd_string_dict[id_to_kegg[coreTBP]])
        else:
            pathDict[coreTBP] = satRxns
            metDict[coreTBP] = satMets

    currencyNodes = ['C' + str(n) for n in Currency]
    mulFac = [1, 3, 3]

    pathList = []
    filtPaths = list(pathDict.values())

    # Need to now calculate the costs, and organize it into
    # a neat, clean data frame.
    pathCosts = [sum([sum(stoich_matrix[rxn][Energy] * mulFac) 
                         for rxn in path]) for path in filtPaths]
    coresProduced = [set(np.nonzero(
                         np.logical_or.reduce(prodMat[path]) * 1)[0]) 
                     & set(Core) for path in filtPaths]
    numCores = [len(cP) for cP in coresProduced]
    nutrientsUsed = [set(np.nonzero(
                         np.logical_or.reduce(rxnMat[path]) * 1)[0]) 
                     & set(nutrientSet) for path in filtPaths]

    # Now I need to add these entries to the data frame if they are unique
    # paths
    for i, p in enumerate(filtPaths):
        pathList.append([nutrientsUsed[i], coresProduced[i], numCores[i],
                         p, pathCosts[i]])

    pathArr = np.array(pathList)
    pathDF = pd.DataFrame(pathArr, columns=['nutrients', 'cores_prod', 'nCores', 'path_rxns', 'cost'])

    return pathDF

#-------------------------------------------------------------------------
# Pulling the data frame and fetching the best metabolic generalist 
# from it.
#-------------------------------------------------------------------------
currPaths = pathDF(rxnMat, prodMat, nutrientSet, Currency, sumRxnVec)
currPaths = currPaths.sort_values(by = ['nCores', 'cost'], ascending=False)
# orgRxns = randOrg(currPaths, prodMat, Core)
# bypVec = secByproducts(orgRxns, rxnMat, prodMat)
print(currPaths)
metGraph = genMetGraph(scopeMets, scopeRxnVec, rxnMat, prodMat)

mGMapping = {n: cpd_string_dict[id_to_kegg[int(n[1:])]] for n in o1g.nodes()}
currencyNodes = ['C' + str(n) for n in Currency]    
o1g.remove_nodes_from(currencyNodes)
o1g = nx.relabel_nodes(o1g, mGMapping)
