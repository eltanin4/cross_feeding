from unlistify import unlistify
from uniqify import uniqify
from sec_byproducts import *
from tqdm import tqdm
from print_progress_bar import print_progress_bar
# from init_generalist_ground import *

othList = np.array( [ 0.0, 1e-2, 5e-2, 1e-1 ])
decList = np.arange(0.0, 1.1, 0.1)
autHM = np.zeros( ( len( othList ), len( decList ) ) )
buildHM = np.zeros( ( len( othList ), len( decList ) ) )

for gamma in tqdm( othList ):
    for alpha, beta in zip( decList, 1 - decList ):
        print( '(A, B, G) = (' + str(alpha) + ', ' + str(beta) + ', ' + str(gamma) + ')')
        currComp = max( [ compFitness( thisOrg, alpha, beta, gamma ) 
                          for thisOrg in randOrgList ] )
        fitterDB = []
        tempdepKinds = []
        depKinds = []
        secsUsed = np.array([])
        tTried = 0
        quitFlag = False
        while len( fitterDB ) < 100:
            print_progress_bar( len(fitterDB), 100, 'Testing scenario')
            if quitFlag:
                break
            while True:
                if tTried > 1e5:
                    quitFlag = True
                    break
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

                    if fitCost( orgRxns ) > 0.0:
                        # print('Successfuly found a working host.')
                        ogCost = fitCost( orgRxns )
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

                    if fitCost( org2Rxns ) > 0.0:
                        # print('Successfuly found a working host.')
                        ogCost = fitCost( org2Rxns )
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
                    # print('At least made a pair.')
                    # if smallestOrgLen * 1.5 >= len(newOrgRxns) and smallestOrgLen * 1.5 >= len(tempOrgRxns):
                    if currComp <= compFitness( newOrgRxns, alpha, beta, gamma ) and currComp <= compFitness( tempOrgRxns, alpha, beta, gamma ):
                        fitterDB.append((newOrgRxns, tempOrgRxns))
                        depKinds.append(tempdepKinds[-1])

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

                        # print('SUCCESS')
                        # print('Fitness: ', currComp, compFitness( newOrgRxns, alpha, beta, gamma ), compFitness( tempOrgRxns, alpha, beta, gamma ))
                        # # print('Biomass yield: ', maxBmsOrgBms, biomassCost(newOrgRxns), biomassCost(tempOrgRxns))
                        # print('Lengths: ', fittestOrgLen, len(newOrgRxns), len(tempOrgRxns))
                        break

        # Received fitterDB.
        if not quitFlag:
            autHM[ np.where( othList == gamma )[0][0], np.where( decList == alpha )[0][0] ] = currComp
            buildHM[ np.where( othList == gamma )[0][0], np.where( decList == alpha )[0][0] ] = np.mean( [ np.mean( [ compFitness( p1, alpha, beta, gamma ) - currComp, compFitness( p2, alpha, beta, gamma ) ] ) for p1, p2 in fitterDB ] )
