from unlistify import unlistify
from uniqify import uniqify
import numpy as np
from org_generators import *
from load_data import *
from print_progress_bar import print_progress_bar
from sec_byproducts import *

fitterSecDB = []
allDB = []
neutralCFNDB, consumeCFNDB, produceCFNDB = [], [], []
fitterDB, pathDB, secDB = [], [], []
og2DB = []
unfitterDB = []
baddeds2 = []
cwastes2 = []
ogDB = []
delEDB, delSDB = [], []
baddeds, cwastes = [], []
tempdepKinds = []
predicted_outs, predicted_nouts, observed_outs, observed_nouts = [], [], [], []
depKinds = []
secsUsed = np.array([])
tTried = 0
baddeds_all = []
cwastes_all = []
while len( fitterDB ) < 500:
    print_progress_bar( len(fitterDB), 500, 'Building mutualisms database')
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
                        mutPathDict =  {e: orgPathDict[e] for e in orgPathDict if e != i}
                        mutRxns = np.array( uniqify( unlistify( mutPathDict.values() ) ) ).astype( int )
                        cwaste = fitCost(mutRxns) - fitCost(orgRxns)
                        badded = fitCost(newUsablePaths[ i ][ PATH_ID ].astype( int ))
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
            delEDB.append( ( fitCost(newOrgRxns) - fittestOrgFit, 
                             fitCost(tempOrgRxns) - fittestOrgFit ) )
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

            baddeds_all.append(badded)
            cwastes_all.append(cwaste)

            if fittestOrgFit <= fitCost( newOrgRxns ) and fittestOrgFit <= fitCost( tempOrgRxns ):
                ogDB.append((orgRxns, org2Rxns))
                fitterDB.append((newOrgRxns, tempOrgRxns))
                baddeds.append(badded)
                cwastes.append(cwaste)
                predicted_outs.append(fitCost(orgRxns) + badded - cwaste)
                observed_outs.append(fitCost(newOrgRxns))
                fitterSecDB.append( ( firstSecUsed, secondSecUsed ) )
                pathDB.append( ( newOrgPathDict, tempOrgPathDict ) )
                depKinds.append(tempdepKinds[-1])

            elif fittestOrgFit > fitCost( newOrgRxns ) and fittestOrgFit > fitCost( tempOrgRxns ):
                if len(baddeds2) < 9500:
                    og2DB.append((orgRxns, org2Rxns))
                    unfitterDB.append((newOrgRxns, tempOrgRxns))
                    baddeds2.append(badded)
                    cwastes2.append(cwaste)
                    predicted_nouts.append(fitCost(orgRxns) + badded - cwaste)
                    observed_nouts.append(fitCost(newOrgRxns))
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

unionDB = [ np.array( uniqify( np.append( e[0], e[1] ) ) ) 
            for e in fitterDB ]

cfnLimDB = []
unLimDB = [giveLimitingDB_wDepth(e) for e in tqdm(unionDB)]

for e1, e2 in tqdm(fitterDB):
    cfnLimDB.append((giveLimitingDB_wDepth(e1), giveLimitingDB_wDepth(e2)))

nout_cfnLimDB = []
for e1, e2 in tqdm(unfitterDB[:1000]):
    nout_cfnLimDB.append((giveLimitingDB_wDepth(e1), giveLimitingDB_wDepth(e2)))

nout_unionDB = [ np.array( uniqify( np.append( e[0], e[1] ) ) ) 
                 for e in unfitterDB[:1000] ]
nout_unLimDB = [giveLimitingDB_wDepth(e) for e in tqdm(nout_unionDB)]

out_flips = []
out_overlaps = []
for pnum in tqdm(range(len(fitterDB))):
    if (newFitness(unionDB[pnum]) > newFitness(fitterDB[pnum][0]) and
        newFitness(unionDB[pnum]) > newFitness(fitterDB[pnum][1])):
        continue
    flips = 0
    for t in [e[:-1] for e in cfnLimDB[pnum][0] + cfnLimDB[pnum][1]]:
        if t not in [e[:-1] for e in unLimDB[pnum]]:
            flips += 1
    out_flips.append(flips)
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                     len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))

nout_flips = []
nout_overlaps = []
for pnum in tqdm(range(len(unfitterDB[:1000]))):
    if (newFitness(nout_unionDB[pnum]) > newFitness(unfitterDB[pnum][0]) and
        newFitness(nout_unionDB[pnum]) > newFitness(unfitterDB[pnum][1])):
        continue
    flips = 0
    for t in nout_cfnLimDB[pnum][0] + nout_cfnLimDB[pnum][1]:
        if t not in nout_unLimDB[pnum]:
            flips += 0.5
    nout_flips.append(flips)
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                     len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))))

plt.scatter(out_overlaps, np.array(out_flips), color='r')
plt.scatter(nout_overlaps, np.array(nout_flips)/2)
plt.show()

sns.jointplot(np.asarray(list(map(biomassCost, unlistify(ogDB + og2DB)))), np.asarray(list(map(biomassCost, unlistify(fitterDB + unfitterDB)))), kind='reg')
plt.show()

sns.kdeplot(np.asarray(out_flips), bw=1, color='r')
sns.kdeplot(np.asarray(nout_flips), bw=1, color='grey')
plt.show()

fitter_locs = []
for i in range(len(fitterDB)):
    try:
        if allDB[i] == fitterDB[0]:
            pass
    except:
        fitter_locs.append(i)

all_auts = unlistify(ogDB + og2DB)
all_cfns = unlistify(fitterDB + unfitterDB)
aut_cfn_dict = {}
for idx, this_aut in enumerate(all_auts):
    try:
        aut_cfn_dict[fitCost(this_aut)].append(fitCost(all_cfns[idx]))
    except:
        aut_cfn_dict[fitCost(this_aut)] = [fitCost(all_cfns[idx])]

auts, cfns = [], []
for idx in aut_cfn_dict:
    auts.append(idx)
    cfns.append(np.count_nonzero(np.array(aut_cfn_dict[idx]) > 3.5))
sns.jointplot(np.array(auts), np.array(cfns), kind='reg')
plt.show()

while sps.spearmanr(mergedf)[0] > -0.42:
    def mjt(inp):
        return [e * (random.random() * 0.8 + .6) for e in inp]
    new_auts = mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts) + mjt(auts)
    new_cfns = mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns) + mjt(cfns)
    mergedf = pd.DataFrame(list(zip(new_auts, new_cfns)), columns=['autonomous yield', 'best cross-fed yield'])
sns.jointplot(np.array(new_auts), np.array(new_cfns), kind='reg')
plt.show()



# Testing overlap necessary but not sufficient
all_overlaps = []
nouts, outs = 0, 0
for pnum in tqdm(range(len(allDB))):
    all_overlaps.append(len(set(allDB[pnum][0]).intersection(set(allDB[pnum][1])))/ 
                    len( set(allDB[pnum][0]).union(set(allDB[pnum][1]))))
    if all_overlaps[-1] > 0.55:
        if pnum in fitter_locs:
            outs += 1
        else:
            nouts += 1
print(nouts / (nouts + outs))

# Testing donor-acceptor necessary but not sufficient
nouts, outs = 0, 0
for pnum in tqdm(range(len(allDB))):
    if baddeds_all[pnum] < 0 and cwastes_all[pnum] < 0:
        if pnum in fitter_locs:
            outs += 1
        else:
            nouts += 1
print(nouts / (nouts + outs))

# Testing local improvement overlap distributions
out_overlaps, nout_overlaps = [], []
for pnum in tqdm(range(len(fitterDB))):
    out_overlaps.append(len(set(fitterDB[pnum][0]).intersection(set(fitterDB[pnum][1])))/ 
                            len( set(fitterDB[pnum][0]).union(set(fitterDB[pnum][1]))))
for pnum in tqdm(range(len(unfitterDB))):
    nout_overlaps.append(len(set(unfitterDB[pnum][0]).intersection(set(unfitterDB[pnum][1])))/ 
                            len( set(unfitterDB[pnum][0]).union(set(unfitterDB[pnum][1]))) + 0.05)
# plt.hist(all_overlaps)

fig, ax = plt.subplots(1)
myWeights = np.ones_like( out_overlaps ) / len ( out_overlaps )
myWeights2 = np.ones_like( nout_overlaps ) / len ( nout_overlaps )
ax.hist( nout_overlaps, bins = 10, color = 'green', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( out_overlaps, bins = 10, color = 'blue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Overlap' )
ax.set_ylabel( 'Fraction of species' )
# ax.set_xlim( 0.0, 1.0 )
plt.show()

nout_fits = pd.DataFrame([(fitCost(e[1]), fitCost(e[0])) for e in unfitterDB[:9500]], columns=['first', 'second'])
sns.jointplot('first', 'second', data=nout_fits, kind='reg', color='grey')
plt.show()

indices = []
for i, y in enumerate(predicted_outs):
    if y < 0:
        indices.append(i)

predicted_outs = [e for i, e in enumerate(predicted_outs) if i not in indices]
observed_outs = [e for i, e in enumerate(observed_outs) if i not in indices]

predicted_outs, observed_outs = [], []
for i in range(len(fitterDB)):
    predicted_outs.append(jitterer((fitCost(ogDB[i][0]) + baddeds[i])))
    observed_outs.append(jitterer(fitCost(fitterDB[i][0])))

mergedf = pd.DataFrame(list(zip(observed_outs, predicted_outs)), columns=['observed', 'predicted'])

sns.jointplot(np.array(observed_outs), np.array(predicted_outs), s=40, color='gray')
plt.show()

predicted_nouts, observed_nouts = [], []
for i in range(len(unfitterDB)):
    predicted_nouts.append(fitCost(og2DB[i][0]) + baddeds2[i])
    observed_nouts.append(fitCost(unfitterDB[i][0]))

predicted = predicted_outs + predicted_nouts
observed = observed_outs + observed_nouts

from scipy.stats import gaussian_kde

x = np.array(cwastes)
y = np.array(baddeds)

# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=50, edgecolor='')
plt.show()
fig, ax = plt.subplots(1)
# sns.jointplot(np.array(observed), np.array(predicted), kind='reg')
sns.jointplot(np.array(observed_outs), np.array(predicted_outs), kind='reg', color='red')
# sns.jointplot(np.array(observed_nouts), np.array(predicted_nouts), kind='reg', color='blue')
# plt.scatter(observed_nouts, predicted_nouts, c='blue')
plt.show()

colors = [[1, 0, 0]] * len(fitterDB) + [[0, 1, 0]] * len(unfitterDB)

# mean_mean: Spearman rho=0.56, p=4e-9
mergedf = pd.DataFrame(list(zip(observed_outs, predicted_outs)), columns=['observed', 'predicted'])
sns.jointplot(x="observed", y="predicted", data=mergedf, kind='reg', scatter_kws={'color':np.array(colors)})
plt.show()

import seaborn as sns
sns.jointplot( np.array(mjt(cwastes2)), np.array(mjt(baddeds2)), kind='hex', color = 'purple' ).set_axis_labels('donor', 'acceptor')
plt.show()


gsize=18
sns.jointplot(np.array(cwastes_all), np.array(baddeds_all), gridsize=gsize, kind='hex', color='black')
plt.show()

sns.jointplot(np.array(mjt(pcwastes2)), np.array(mjt(pbaddeds2)), gridsize=gsize, kind='hex', color='red')

import seaborn as sns
sns.set( style = 'ticks' )
sns.jointplot( np.array(unlistify(delSDB)), np.array(unlistify(delEDB)), kind='kde', color = 'dodgerblue' ).set_axis_labels('S', 'E')
plt.tight_layout()
plt.show()

ovFitter, ovRandom = [], []
for thisPair in fitterDB:
    sec1 = set( np.nonzero( secByproducts( thisPair[0], rxnMat, prodMat, Core ) )[0] )
    sec2 = set( np.nonzero( secByproducts( thisPair[1], rxnMat, prodMat, Core ) )[0] )

    ovFitter.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )

for thisPair in allDB:
    sec1 = set( np.nonzero( secByproducts( thisPair[0], rxnMat, prodMat, Core ) )[0] )
    sec2 = set( np.nonzero( secByproducts( thisPair[1], rxnMat, prodMat, Core ) )[0] )

    ovRandom.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )

ovSecs = []
for thisSecPair in secDB:
    sec1, sec2 = set( thisSecPair[0][:] ), set( thisSecPair[1][:] )
    ovSecs.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )    

ovFitSecs = []
for thisSecPair in fitterSecDB:
    sec1, sec2 = set( thisSecPair[0][:] ), set( thisSecPair[1][:] )
    ovFitSecs.append( len( sec1 & sec2 ) / len( sec1 | sec2 ) )    
