from print_progress_bar import print_progress_bar
from unlistify import unlistify
from uniqify import uniqify
import numpy as np
from org_generators import *
from sec_byproducts import secByproducts, pathwaySecByproducts
from load_data import *
from tqdm import tqdm

nutrientSet = [25, 12]

def nutrientsConsumed( rxns, currNutrientSet=nutrientSet ):
    return abs( sum( [ sum( stoich_matrix[ rxn ][ nutrientSet ] ) for rxn in rxns.astype(int)
                       if sum( stoich_matrix[ rxn ][ currNutrientSet ] ) < 0.0 ] ) )

def redoxFitCost( rxns, currExchangeRatio=3, currNutrientSet=nutrientSet ):
    mulFac = [1, currExchangeRatio, currExchangeRatio]
    return int( sum( [ sum( stoich_matrix[ rxn ][ Energy ] * mulFac ) 
                for rxn in rxns.astype(int) ] ) ) / nutrientsConsumed( rxns, currNutrientSet )

saveDir = 'prebSet_networks/'
secSaveDir = 'sec_' + saveDir
prefix = 'prebSeed_randRxnSets'

medType = 'preb'
mediumFile = medType + '_medium.txt'
nutrientSet = [kegg_to_id[i] for i in np.genfromtxt(mediumFile)]

# Pulling out the pathway dictionary from the pregenerated folder.
pathDict = genPathDict( saveDir, prefix, Core )
secPathDict = genPathDict( secSaveDir, prefix, Core )

# Generating random organism generalist reaction sets.
pNeeded = []
oNeeded = []
for currExchangeRatio in np.arange(0.0, 3.1, 0.5):
    oNeeded.append( [] )
    NUM_RAND_ORGS = 50000

    randOrgLenList, randOrgFitList, randOrgBmsList, randOrgComList = np.array( [] ), np.array( [] ), np.array( [] ), np.array( [] )
    randOrgList = []
    numByps = np.array( [] )
    bypArr = np.array( [] )
    for i in tqdm( range( NUM_RAND_ORGS ) ):
        # Generating a random organism from a pregenerated 
        # path dictionary.
        orgRxns = genRandOrg( pathDict )

        # Calculating the size of the organism and its fitness.
        randOrgLenList = np.append( randOrgLenList, len( orgRxns ) )
        randOrgFitList = np.append( randOrgFitList, redoxFitCost( orgRxns, currExchangeRatio ) )

        # Keeping a list of generated organisms.
        randOrgList.append( orgRxns )
        
        # Updating the list of byproducts in the network.
        bypArr = np.append( bypArr, np.nonzero( secByproducts( orgRxns, rxnMat, prodMat, Core ) )[0] )
        numByps = np.append( numByps, len( np.nonzero( secByproducts( orgRxns, rxnMat, prodMat, Core ) )[0] ) )

    # Generating the fittest generalist characteristics.
    fittestOrgLen = min( randOrgLenList[ np.where( randOrgFitList == max( randOrgFitList ) )[0] ] )
    fittestOrgFit = max( randOrgFitList[ np.where( randOrgFitList == max( randOrgFitList ) )[0] ] )

#-------------------------------------------------------------------------

    if fittestOrgFit < 0.0:
        pNeeded.append(0.0)
        continue

#-------------------------------------------------------------------------

    fitterSecDB = []
    allDB = []
    neutralCFNDB, consumeCFNDB, produceCFNDB = [], [], []
    fitterDB, pathDB, secDB = [], [], []
    delEDB, delSDB = [], []
    tempdepKinds = []
    depKinds = []
    secsUsed = np.array([])
    tTried = 0
    while len( secDB ) < 10000:
        print_progress_bar( len(secDB), 10000, 'Building mutualisms database')
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

                if redoxFitCost( orgRxns, currExchangeRatio ) > 0.0:
                    # print('Successfuly found a working host.')
                    ogCost = redoxFitCost( orgRxns, currExchangeRatio )
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

                if redoxFitCost( org2Rxns, currExchangeRatio ) > 0.0:
                    # print('Successfuly found a working host.')
                    ogCost = redoxFitCost( org2Rxns, currExchangeRatio )
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
                delEDB.append( ( redoxFitCost(newOrgRxns, currExchangeRatio) - fittestOrgFit, 
                                 redoxFitCost(tempOrgRxns, currExchangeRatio) - fittestOrgFit ) )
                delSDB.append( ( len(newOrgRxns) - fittestOrgLen, 
                                 len(tempOrgRxns) - fittestOrgLen ) )

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
                if fittestOrgFit <= redoxFitCost( newOrgRxns, currExchangeRatio ) and fittestOrgFit <= redoxFitCost( tempOrgRxns, currExchangeRatio ):
                    fitterDB.append((newOrgRxns, tempOrgRxns))
                    fitterSecDB.append( ( firstSecUsed, secondSecUsed ) )
                    pathDB.append( ( newOrgPathDict, tempOrgPathDict ) )
                    depKinds.append(tempdepKinds[-1])
                    oNeeded[ -1 ] += [ float( np.mean( [ redoxFitCost( newOrgRxns,                                 currExchangeRatio ) - fittestOrgFit,
                                                         redoxFitCost( tempOrgRxns, currExchangeRatio) - fittestOrgFit ] )               ) ]

                break

    pNeeded.append( len( fitterDB ) / len( secDB ) )
