import sys

thisNum = int(sys.argv[1])
medType = 'preb'

from give_reverse_scope import giveRevScope
from load_data import *
from sec_byproducts import secByproducts
from random_minimal_subgraph import randMinSubnet
from bool2int import bool2int
from org_generators import *
from unlistify import unlistify
import os
import numpy as np

mediumFile = medType + '_medium.txt'
pullDir = medType + 'Set_networks/'
prefix = medType + 'Seed_randRxnSets'
saveDir = 'new_sec_' + medType + 'Set_networks/'

# Define the seed set, this would the nutrients in the medium (see line 1)
nutrientSet = [kegg_to_id[i] for i in np.genfromtxt(mediumFile)]

NUM_PATHS = 10

def savePathsWithSecVec( currSecVec ):
    currSecMets = list( np.nonzero( currSecVec )[0] )
    augNutrientSet = currSecMets[:]

    if not currSecVec.any():
        return
    elif not os.path.exists( saveDir + str( bool2int(currSecVec) ) + '_prebSeed_randRxnSets_core' + str( Core[-1] ) + '_' + str( NUM_PATHS - 1 ) + '.txt' ):
        for coreTBP in Core:
            print('Finding paths to ' + cpd_string_dict[id_to_kegg[coreTBP]] + '.')
            coreProdRxns = (prodMat[:, coreTBP] == 1) * 1
            satMetVec, satRxnVec = giveRevScope(rxnMat, prodMat, sumRxnVec, 
                                                augNutrientSet, Currency, coreTBP)
            satMets, satRxns = np.nonzero(satMetVec)[0], np.nonzero(satRxnVec)[0]

            # Randomizing pathways to coreTBP.
            rMinSubnetsRxns = []
            for netNum in range(NUM_PATHS):
                rMinSubnetsRxns.append(randMinSubnet( satRxnVec, rxnMat, prodMat, sumRxnVec, 
                                                      coreProdRxns, coreTBP, augNutrientSet, 
                                                      Currency) )

                # Saving the made pathway.
                np.savetxt( saveDir + str( bool2int( currSecVec ) ) + '_prebSeed_randRxnSets_core' + str( coreTBP ) + '_' + str( netNum ) + '.txt', rMinSubnetsRxns[-1] )

vecsAdd = np.genfromtxt( 'vecsAdd.txt' ).astype( int )
savePathsWithSecVec( vecsAdd[ thisNum ] )
