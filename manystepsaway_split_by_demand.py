from unlistify import unlistify
from uniqify import uniqify
import numpy as np
from load_data import *
from print_progress_bar import print_progress_bar
from sec_byproducts import *
import random
from org_generators import *
from sec_byproducts import secByproducts
from fit_cost import fitCost, biomassCost, compFitness
from tqdm import tqdm
from split_by_demand import splitByDemand

saveDir = 'prebSet_networks/'
secSaveDir = 'sec_' + saveDir
prefix = 'prebSeed_randRxnSets'

medType = 'preb'
mediumFile = medType + '_medium.txt'
nutrientSet = [kegg_to_id[i] for i in np.genfromtxt(mediumFile)]

# Pulling out the pathway dictionary from the pregenerated folder.
pathDict = genPathDict( saveDir, prefix, Core )
secPathDict = genPathDict( secSaveDir, prefix, Core )

# Defining the new fitness measurements
def newFitness( orgRxns ):
    return splitByDemand( stoich_matrix, rxnMat, prodMat, 
                          sumRxnVec, rho, pi, nutrientSet, 
                          Energy, Currency, Core, orgRxns )[0]

def splitDemandByps( orgRxns, givenByps ):
    return splitByDemand( stoich_matrix, rxnMat, prodMat, 
                          sumRxnVec, rho, pi, nutrientSet, 
                          Energy, Currency, Core, orgRxns, givenByps )

def newBiomass( orgRxns ):
    return splitByDemand( stoich_matrix, rxnMat, prodMat, 
                          sumRxnVec, rho, pi, nutrientSet, 
                          Energy, Currency, Core, orgRxns )[1]

fitterSecDB = []
allDB = []
neutralCFNDB, consumeCFNDB, produceCFNDB = [], [], []
fitterDB, pathDB, secDB = [], [], []
delEDB, delSDB = [], []
tempdepKinds = []
depKinds = []
secsUsed = np.array([])
tTried = 0
while len( fitterDB ) < 100:
    print_progress_bar( len(secDB), 1000, 'Building mutualisms database')
    while True:
        tFlag = False
        tTried += 1
        # Generate a random organism.
        while True:
            orgPathDict = {}
            for coreTBP in Core:
                orgPathDict[ coreTBP ] = list( random.choice( pathDict[ coreTBP ] ).astype(int) )
            
            # Storing the list of reactions.
            orgRxns = np.array( uniqify( unlistify( orgPathDict.values() ) ) ).astype( int )

            # Creating a dictionary of secretions.
            orgSecDict = {}
            for coreTBP in Core:
                orgSecDict[ coreTBP ] = list( np.nonzero( pathwaySecByproducts( 
                                              orgPathDict[ coreTBP ], orgRxns, 
                                              rxnMat, prodMat, Core ) )[0] )
            
            # Purging duplicates to have only unique byproducts.
            tempSet = set()
            fullSet = unlistify( orgSecDict.values() )
            duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
            orgSecDict = { coreTBP: list( set( orgSecDict[ coreTBP ] ).difference( duplicates ) ) 
                           for coreTBP in Core }

            # Updating the full list of byproducts in the network.
            bypArr = uniqify( unlistify( orgSecDict.values() ) )

            if newFitness( orgRxns ) > 0.0:
                # print('Successfuly found a working host.')
                ogCost = newFitness( orgRxns )
                break

    #-------------------------------------------------------------------------

        # Generate a random organism.
        while True:
            org2PathDict = {}
            for coreTBP in Core:
                org2PathDict[ coreTBP ] = list( random.choice( pathDict[ coreTBP ] ).astype(int) )
            
            # Storing the list of reactions.
            org2Rxns = np.array( uniqify( unlistify( org2PathDict.values() ) ) ).astype( int )

            # Creating a dictionary of secretions.
            org2SecDict = {}
            for coreTBP in Core:
                org2SecDict[ coreTBP ] = list( np.nonzero( pathwaySecByproducts( 
                                              org2PathDict[ coreTBP ], org2Rxns, 
                                              rxnMat, prodMat, Core ) )[0] )
            
            # Purging duplicates to have only unique byproducts.
            tempSet = set()
            fullSet = unlistify( org2SecDict.values() )
            duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
            org2SecDict = { coreTBP: list( set( org2SecDict[ coreTBP ] ).difference( duplicates ) ) 
                           for coreTBP in Core }

            # Updating the full list of byproducts in the network.
            bypArr = uniqify( unlistify( org2SecDict.values() ) )

            if newFitness( org2Rxns ) > 0.0:
                # print('Successfuly found a working host.')
                ogCost = newFitness( org2Rxns )
                break

    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------

        # Creating a list of pathways that work on the respective secretions.
        usablePaths = {}
        for coreTBP in Core:
            usablePaths[ coreTBP ] = []
            tempBypList = orgSecDict[ coreTBP ][:]
            if not tempBypList:
                continue
            for currPath in secPathDict[ coreTBP ]:
                if ((np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)
                    [ tempBypList ].any()):
                    usablePaths[ coreTBP ].append( currPath )
                    firstSecUsed = tempBypList[:]

        #-------------------------------------------------------------------------

        for coreTBP in Core:
            tFlag = False
            if not usablePaths[ coreTBP ]:
                continue

        #-------------------------------------------------------------------------

            # Randomly picking a usable path.
            # print('First replacing a pathway for core ' + str( coreTBP ) )
            tempdepKinds.append([coreTBP])
            PATH_ID = random.choice( range( len( usablePaths[ coreTBP ] ) ) )
            tempOrgPathDict = org2PathDict.copy()
            tempOrgPathDict[ coreTBP ] = list( usablePaths[ coreTBP ][ PATH_ID ].astype( int ) )
            tempOrgRxns = np.array( uniqify( unlistify( tempOrgPathDict.values() ) ) ).astype( int )
            replacedPathSec = list( np.nonzero( pathwaySecByproducts( 
                                    tempOrgPathDict[ coreTBP ], tempOrgRxns, 
                                    rxnMat, prodMat, Core ) )[0] )

            # Creating a dictionary of secretions.
            tempOrgSecDict = {}
            for i in Core:
                tempOrgSecDict[ i ] = list( np.nonzero( pathwaySecByproducts( 
                                              tempOrgPathDict[ i ], tempOrgRxns, 
                                              rxnMat, prodMat, Core ) )[0] )
            
            # Purging duplicates to have only unique byproducts.
            tempSet = set()
            fullSet = unlistify( tempOrgSecDict.values() )
            duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
            tempOrgSecDict = { i: list( set( tempOrgSecDict[ i ] ).difference( duplicates ) ) 
                           for i in Core }

        #-------------------------------------------------------------------------

            if replacedPathSec:
                newUsablePaths = {}
                for i in Core:
                    if i != coreTBP:
                        newUsablePaths[ i ] = []
                        for currPath in secPathDict[ i ]:
                            if ((np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)
                                [ replacedPathSec ].any()):
                                newUsablePaths[ i ].append( currPath )
                                secondSecUsed = replacedPathSec[:]
            else:
                continue
        #-------------------------------------------------------------------------
            
            flag = False
            for i in np.random.permutation(Core):
                if i != coreTBP:
                    if newUsablePaths[ i ]:
                        # print('Finally replacing pathway for core ' + str( i ) )
                        tempdepKinds[-1].append(i)
                        flag = True
                        PATH_ID = random.choice( range ( len( newUsablePaths[ i ] ) ) )
                        newOrgPathDict = orgPathDict.copy()
                        newOrgPathDict[ i ] = list( newUsablePaths[ i ][ PATH_ID ].astype( int ) )
                        newOrgRxns = np.array( uniqify( unlistify( newOrgPathDict.values() ) ) ).astype( int )
                        break
            if not flag:
                continue

            # Creating a dictionary of secretions.
            newOrgSecDict = {}
            for i in Core:
                newOrgSecDict[ i ] = list( np.nonzero( pathwaySecByproducts( 
                                              newOrgPathDict[ i ], newOrgRxns, 
                                              rxnMat, prodMat, Core ) )[0] )
            
            # Purging duplicates to have only unique byproducts.
            tempSet = set()
            fullSet = unlistify( newOrgSecDict.values() )
            duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
            newOrgSecDict = { i: list( set( newOrgSecDict[ i ] ).difference( duplicates ) ) 
                           for i in Core }
            tFlag = True
            break

        #-------------------------------------------------------------------------
        
        if tFlag:
            allDB.append( ( newOrgRxns, tempOrgRxns ) ) 
            # print('At least made a pair.')
            tempOrgSecs = uniqify(unlistify(tempOrgSecDict.values()))
            newOrgSecs = uniqify(unlistify(newOrgSecDict.values()))

            tempOrgSecPathDict = {}
            for coreTBP in Core:
                currPath = tempOrgPathDict[ coreTBP ]
                tempOrgSecPathDict[ coreTBP ] = [ sec for sec in newOrgSecs 
                                                  if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ sec ] ]

            newOrgSecPathDict = {}
            for coreTBP in Core:
                currPath = newOrgPathDict[ coreTBP ]
                newOrgSecPathDict[ coreTBP ] = [ sec for sec in tempOrgSecs 
                                                  if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ sec ] ]

            secsUsed = np.append(secsUsed, unlistify(tempOrgSecPathDict.values()) + unlistify(newOrgSecPathDict.values()))
            
            secDB.append( ( firstSecUsed, secondSecUsed ) )
            delEDB.append( ( newFitness(newOrgRxns) - fittestOrgFit, 
                             newFitness(tempOrgRxns) - fittestOrgFit ) )
            # delSDB.append( ( len(newOrgRxns) - fittestOrgLen, 
            #                  len(tempOrgRxns) - fittestOrgLen ) )

            newProf = np.array( [ sum( stoich_matrix[ rxn ][ Energy ] * [ 1, 3, 3 ] ) 
                                  for rxn in newOrgRxns ] )
            oldProf = np.array( [ sum( stoich_matrix[ rxn ][ Energy ] * [ 1, 3, 3 ] ) 
                                  for rxn in tempOrgRxns ] )
            fracNeutral1, fracNeutral2 = np.count_nonzero( newProf == 0 ) / len( newOrgRxns ), np.count_nonzero( oldProf == 0 ) / len( tempOrgRxns )
            fracConsuming1, fracConsuming2 = np.count_nonzero( newProf < 0 ) / len( newOrgRxns ), np.count_nonzero( oldProf < 0 ) / len( tempOrgRxns )
            fracProducing1, fracProducing2 = np.count_nonzero( newProf > 0 ) / len( newOrgRxns ), np.count_nonzero( oldProf > 0 ) / len( tempOrgRxns )

            neutralCFNDB.append( ( fracNeutral1, fracNeutral2 ) )
            consumeCFNDB.append( ( fracConsuming1, fracConsuming2 ) )
            produceCFNDB.append( ( fracProducing1, fracProducing2 ) )

            # neutralAUTDB, consumeAUTDB, produceAUTDB = [], [], []
            # for currOrgRxns in tqdm( randOrgList[ : 20000] ):
            #     currProf = np.array( [ sum( stoich_matrix[ rxn ][ Energy ] * [ 1, 3, 3 ] ) 
            #                       for rxn in currOrgRxns ] )
            #     neutralAUTDB.append( np.count_nonzero( currProf == 0 ) / len( currOrgRxns) )
            #     consumeAUTDB.append( np.count_nonzero( currProf < 0 ) / len( currOrgRxns) )
            #     produceAUTDB.append( np.count_nonzero( currProf > 0 ) / len( currOrgRxns) )

            # if smallestOrgLen * 1.5 >= len(newOrgRxns) and smallestOrgLen * 1.5 >= len(tempOrgRxns):
            if fittestOrgFit <= newFitness( newOrgRxns ) and fittestOrgFit <= newFitness( tempOrgRxns ):
                fitterDB.append((newOrgRxns, tempOrgRxns))
                fitterSecDB.append( ( firstSecUsed, secondSecUsed ) )
                pathDB.append( ( newOrgPathDict, tempOrgPathDict ) )
                depKinds.append(tempdepKinds[-1])

            break
                # tempOrgSecs = uniqify(unlistify(tempOrgSecDict.values()))
                # newOrgSecs = uniqify(unlistify(newOrgSecDict.values()))

                # tempOrgSecPathDict = {}
                # for coreTBP in Core:
                #     currPath = tempOrgPathDict[ coreTBP ]
                #     tempOrgSecPathDict[ coreTBP ] = [ sec for sec in newOrgSecs 
                #                                       if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ sec ] ]

                # newOrgSecPathDict = {}
                # for coreTBP in Core:
                #     currPath = newOrgPathDict[ coreTBP ]
                #     newOrgSecPathDict[ coreTBP ] = [ sec for sec in tempOrgSecs 
                #                                       if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ sec ] ]

                # secsUsed = np.append(secsUsed, unlistify(tempOrgSecPathDict.values()) + unlistify(newOrgSecPathDict.values()))
                
                # secDB.append( ( firstSecUsed, secondSecUsed ) )

                # print('SUCCESS')
                # print('Energy yield: ', fittestOrgFit, fitCost(newOrgRxns), fitCost(tempOrgRxns))
                # # print('Biomass yield: ', maxBmsOrgBms, biomassCost(newOrgRxns), biomassCost(tempOrgRxns))
                # print('Lengths: ', fittestOrgLen, len(newOrgRxns), len(tempOrgRxns))
            # elif maxBmsOrgBms <= biomassCost(newOrgRxns) or maxBmsOrgBms <= biomassCost(tempOrgRxns):
            #     print('Biomass yield: ', maxBmsOrgBms, biomassCost(newOrgRxns), biomassCost(tempOrgRxns))
            #     print('Lengths: ', fittestOrgLen, len(newOrgRxns), len(tempOrgRxns))


# import seaborn as sns
# sns.set( style = 'ticks' )
# sns.jointplot( np.array(unlistify(delSDB)), np.array(unlistify(delEDB)), kind='kde', color = 'dodgerblue' ).set_axis_labels('S', 'E')
# plt.tight_layout()
# plt.show()

# ovFitter, ovRandom = [], []
# for thisPair in fitterDB:
#     sec1 = set( np.nonzero( secByproducts( thisPair[0], rxnMat, prodMat, Core ) )[0] )
#     sec2 = set( np.nonzero( secByproducts( thisPair[1], rxnMat, prodMat, Core ) )[0] )

#     ovFitter.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )

# for thisPair in allDB:
#     sec1 = set( np.nonzero( secByproducts( thisPair[0], rxnMat, prodMat, Core ) )[0] )
#     sec2 = set( np.nonzero( secByproducts( thisPair[1], rxnMat, prodMat, Core ) )[0] )

#     ovRandom.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )

# ovSecs = []
# for thisSecPair in secDB:
#     sec1, sec2 = set( thisSecPair[0][:] ), set( thisSecPair[1][:] )
#     ovSecs.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )    

# ovFitSecs = []
# for thisSecPair in fitterSecDB:
#     sec1, sec2 = set( thisSecPair[0][:] ), set( thisSecPair[1][:] )
#     ovFitSecs.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )    

allFitnesses = unlistify( [ [ e[0] + fittestOrgFit, e[1] + fittestOrgFit ] 
                              for e in delEDB ] )
allFitDB = [ [ e[0] + fittestOrgFit, e[1] + fittestOrgFit ] 
                              for e in delEDB ]
meanFitDB = [ np.mean([ e[0] + fittestOrgFit, e[1] + fittestOrgFit ] )
                              for e in delEDB ]
unionFitnesses = [ newFitness(e) for e in tqdm(unionDB) ]
