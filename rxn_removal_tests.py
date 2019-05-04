import numpy as np
from load_data import *
from uniqify import uniqify
from unlistify import unlistify

coreProdRxns = { coreTBP: (prodMat[:, coreTBP] == 1) * 1 for coreTBP in Core }
frglList = [ ] 
NUM_RAND_ORGS = 1000

for thisIter in tqdm( range( NUM_RAND_ORGS ) ):
    orgPathDict = {}
    for coreTBP in Core:
        orgPathDict[ coreTBP ] = list( random.choice( pathDict[ coreTBP ] ).astype(int) )
    
    # Storing the list of reactions.
    orgRxns = np.array( uniqify( unlistify( orgPathDict.values() ) ) ).astype( int )

    # Moving through all reactions in the organism.
    remRxn = np.random.choice( orgRxns )
    thisRxnList = [ ]
    for coreTBP in Core:
        if remRxn in orgPathDict[ coreTBP ]:
            thisRxnList.append( coreTBP )

    frglList.append( len( thisRxnList ) )

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1)
myWeights = np.ones_like( frglList ) / len ( frglList )
ax.hist( frglList, bins = 7, color = 'gray', weights = myWeights)
ax.set_xlabel( 'Number of precursors lost' )
ax.set_ylabel( 'Normalised frequency' )
ax.set_xlim(0, len(Core))
fig.tight_layout()
plt.show()

# Now plotting this for the cross-feeding networks.
cfnList = [ ]
thisIter = -1
for thisPair in tqdm( pathDB ):
    thisIter += 1
    o1PathDict, o2PathDict = thisPair[0], thisPair[1]
    o1Rxns = set( uniqify( unlistify( o1PathDict.values() ) ) )
    o2Rxns = set( uniqify( unlistify( o2PathDict.values() ) ) )
    allRxns = np.array( list( o1Rxns | o2Rxns ) )
    o1Rxns, o2Rxns = np.array( list( o1Rxns ) ), np.array( list( o2Rxns ) )
    
    # Checking for main module loss.
    remRxn1, remRxn2 = np.random.choice( o1Rxns ), np.random.choice( o2Rxns )
    o1RxnList, o2RxnList = [ ], [ ]
    for coreTBP in Core:
        if remRxn1 in o1PathDict[ coreTBP ]:
            o1RxnList.append( coreTBP )
        if remRxn2 in o2PathDict[ coreTBP ]:
            o2RxnList.append( coreTBP )

    # Checking for secreted goods now
    o1SatRxnVec, o2SatRxnVec = np.zeros( len( rxnMat ) ), np.zeros( len( rxnMat ) )
    o1SatRxnVec[ o1Rxns ], o2SatRxnVec[ o2Rxns ] = 1, 1
    o1SatRxnVec[ remRxn1 ], o2SatRxnVec[ remRxn2 ] = 0, 0
    o1SatMetVec, o1SatRxnVec = markSatMetsRxns( o1SatRxnVec, rxnMat, prodMat, 
                                                sumRxnVec, [], nutrientSet, 
                                                Currency )

    o2SatMetVec, o2SatRxnVec = markSatMetsRxns( o2SatRxnVec, rxnMat, prodMat, 
                                                sumRxnVec, [], nutrientSet, 
                                                Currency )

    if not o1SatMetVec[ secDB[ thisIter ][ 0 ] ].all():
        if coreTBP not in o2RxnList:
            o2RxnList.append( coreTBP )

    if not o2SatMetVec[ secDB[ thisIter ][ 1 ] ].all():
        if coreTBP not in o1RxnList:
            o1RxnList.append( coreTBP )

    cfnList.append( len( o1RxnList ) )
    cfnList.append( len( o2RxnList ) )

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1)
myWeights = np.ones_like( cfnList ) / len ( cfnList )
ax.hist( cfnList, bins = 5, color = 'gray', weights = myWeights)
ax.set_xlabel( 'Number of precursors lost' )
ax.set_ylabel( 'Normalised frequency' )
ax.set_xlim(0, len(Core))
fig.tight_layout()
plt.show()
