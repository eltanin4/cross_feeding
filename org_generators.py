import numpy as np
import random
from uniqify import uniqify
from gen_scopeGraphs import genMetGraph
from load_data import *
from fit_cost import fitCost, nutrientsConsumed

def genPathDict( saveDir, prefix, Core ):
    pathDict = {}

    # Iterating over all core molecules.
    for coreTBP in Core:
        pathDict[ coreTBP ] = []

        for i in range( 10 ):
            # Loading the relevant files.
            currRxns = np.genfromtxt( saveDir + prefix + '_core' + str(coreTBP) + '_' + str(i) + '.txt' ).astype(int)

            # Failsafe for when only one reaction is in the list.
            currRxns = np.append( np.array( [] ), currRxns )

            # Checking if the pathway consumes a nutrient in the set.
            # if not nutrientsConsumed( currRxns ):
            #     continue

            # Adding it to the path dictionary.
            pathDict[ coreTBP ].append( currRxns )

    return pathDict

#-------------------------------------------------------------------------

# def secGenPathDict( saveDir, prefix, matchLib, Core ):
#     pathDict = { ID : {} for ID in matchLib.values() }
#     exceptions  = set()

#     # Iterating over all core molecules.
#     for ID in list( matchLib.values() ):
#         for coreTBP in Core:
#             pathDict[ ID ][ coreTBP ] = []

#             for i in range( 10 ):
#                 # Loading the relevant files.
#                 try:
#                     currRxns = np.genfromtxt( saveDir +  str(ID) + '_' + prefix + '_core' + str(coreTBP) + '_' + str(i) + '.txt' ).astype(int)

#                     # Failsafe for when only one reaction is in the list.
#                     currRxns = np.append( np.array( [] ), currRxns )

#                     # Adding it to the path dictionary.
#                     pathDict[ ID ][ coreTBP ].append( currRxns )
#                 except:
#                     exceptions.add( ID )

#     for thisID in exceptions:
#         pathDict.pop( thisID )

#     return pathDict

#-------------------------------------------------------------------------

def secGenPathDict( saveDir, prefix, matchLib, Core ):
    pathDict = {}
    exceptions  = set()

    # Iterating over all core molecules.
    for ID in list( matchLib.values() ):
        for coreTBP in Core:
            pathDict[ coreTBP ] = []

            for i in range( 10 ):
                # Loading the relevant files.
                try:
                    currRxns = np.genfromtxt( saveDir +  str(ID) + '_' + prefix + '_core' + str(coreTBP) + '_' + str(i) + '.txt' ).astype(int)

                    # Failsafe for when only one reaction is in the list.
                    currRxns = np.append( np.array( [] ), currRxns )

                    # Adding it to the path dictionary.
                    if currRxns not in pathDict[ coreTBP ]:
                        pathDict[ coreTBP ].append( currRxns )
                except:
                    exceptions.add( ID )

    return pathDict

#-------------------------------------------------------------------------
def genRandOrg( pathDict ):
    orgRxns = np.array( [ ] )
    for coreTBP in Core:

        # Picking a path at random from the dictionary that generates
        # the current core molecule.
        orgRxns = np.append( orgRxns, random.choice( pathDict[ coreTBP ] ) )
    
    # Returning the unique bunch of reactions that correspond 
    # to the individual.
    return np.array( uniqify( orgRxns ) ).astype(int)

#-------------------------------------------------------------------------

def sortedGenOrg( pathDict, sortFunc = fitCost, optFunc = max ):
    orgRxns = np.array( [ ] )
    for coreTBP in Core:

        # Creating a list of size of pathways that produce the
        # current core molecule.
        sortList = [ sortFunc( path ) for path in pathDict[ coreTBP ] ]

        # Picking that path which has smallest size.
        orgRxns = np.append( orgRxns, pathDict[ coreTBP ][ sortList.index( optFunc( sortList ) ) ] )
    
    # Returning the unique bunch of reactions that correspond 
    # to the individual.
    return np.array( uniqify( orgRxns ) ).astype(int)

#-------------------------------------------------------------------------

def jitterer(inp):
    return inp * (random.random() * 0.1 + 0.95)

def mjt(inp):
    return [e * (random.random() * 0.4 + .8) for e in inp]

#-------------------------------------------------------------------------

def saveOrgMetGraph( orgRxns, orgFileName ):
    import networkx as nx
    orgRxnVec = np.zeros( len( rxnMat ) )
    orgRxnVec[ orgRxns ] = 1

    # Generating the corresponding metabolite graph for the organism.
    orgMetGraph = genMetGraph( orgRxnVec, rxnMat, prodMat )

    # Mapping the metabolite IDs in the graph to names.
    mGMapping = { n: cpd_string_dict[ id_to_kegg[ int( n[ 1 : ] ) ] ] 
                  for n in orgMetGraph.nodes() }

    # Removing the currency nodes so that they don't muck 
    # up the visualization.
    currencyNodes = ['C' + str(n) for n in Currency]
    orgMetGraph.remove_nodes_from(currencyNodes)

    # Relabeling the nodes to names now.
    orgMetGraph = nx.relabel_nodes(orgMetGraph, mGMapping)

    # Saving the metabolite graph corresponding to the organism.
    nx.write_graphml( orgMetGraph, orgFileName )

#-------------------------------------------------------------------------

def giveOrgMetGraph( orgRxns ):
    import networkx as nx
    orgRxnVec = np.zeros( len( rxnMat ) )
    orgRxnVec[ orgRxns ] = 1

    # Generating the corresponding metabolite graph for the organism.
    orgMetGraph = genMetGraph( orgRxnVec, rxnMat, prodMat )

    # Mapping the metabolite IDs in the graph to names.
    mGMapping = { n: cpd_string_dict[ id_to_kegg[ int( n[ 1 : ] ) ] ] 
                  for n in orgMetGraph.nodes() }

    # Removing the currency nodes so that they don't muck 
    # up the visualization.
    currencyNodes = ['C' + str(n) for n in Currency]
    orgMetGraph.remove_nodes_from(currencyNodes)

    # Relabeling the nodes to names now.
    orgMetGraph = nx.relabel_nodes(orgMetGraph, mGMapping)

    # Returning the metabolite graph corresponding to the organism.
    return orgMetGraph

#-------------------------------------------------------------------------
