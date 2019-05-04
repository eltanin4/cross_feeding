import numpy as np
from sat_marker import markSatMetsRxns

def pruneOrgNet(satRxns, rxnMat, prodMat, 
                sumRxnVec, Core, nutrientSet, Currency):
    currSatRxnVec = np.zeros( len( rxnMat ) )
    currSatRxnVec[ satRxns ] = 1

    while True:
        # Keeping track of what the graph currently looks like.
        currSatRxns = np.nonzero( currSatRxnVec )[ 0 ]

        # Marking out which reactions are singly removable.
        removableMets = np.array( [ ] )
        for remRxn in currSatRxns:
            tempSatRxnVec = np.copy( currSatRxnVec )
            tempSatRxnVec[ remRxn ] = 0
            tempSatMetVec, tempSatRxnVec = markSatMetsRxns( tempSatRxnVec, rxnMat, prodMat, 
                                                            sumRxnVec, [], nutrientSet, Currency )
            if tempSatMetVec[ Core ].all():
                removableMets = np.append( removableMets, remRxn )

        # If no more reactions can be pruned, exit with the current vector state.
        if not removableMets.any():
            return np.nonzero( currSatRxnVec )[ 0 ]

        # Randomly permuting the singly removable reactions.
        removalOrder = np.random.permutation( removableMets )

        # Calling a vector of reaction that can be removed.
        for remRxn in removalOrder.astype(int):
            # Mark what is achievable when the reaction is removed.
            tempSatRxnVec = np.copy( currSatRxnVec )
            tempSatRxnVec[ remRxn ] = 0
            tempSatMetVec, tempSatRxnVec = markSatMetsRxns( tempSatRxnVec, rxnMat, prodMat, 
                                                            sumRxnVec, [], nutrientSet, Currency )

            # If all core molecules can still be produced post removal, remove the reaction.
            if tempSatMetVec[ Core ].all():
                currSatRxnVec[ remRxn ] = 0
            else:
                break
