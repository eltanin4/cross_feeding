medType = 'preb'

import numpy as np
import matplotlib.pyplot as plt
from org_generators import *
from sec_byproducts import secByproducts
from load_data import *
from fit_cost import *
from unlistify import unlistify
from bool2int import bool2int
from print_progress_bar import print_progress_bar
from give_reverse_scope import giveRevScope
from random_minimal_subgraph import randMinSubnet
import os

mediumFile = medType + '_medium.txt'
pullDir = medType + 'Set_networks/'
prefix = medType + 'Seed_randRxnSets'
saveDir = 'new_sec_' + medType + 'Set_networks/'

# Define the seed set, this would the nutrients in the medium (see line 1)
nutrientSet = [kegg_to_id[i] for i in np.genfromtxt(mediumFile)]

# Pulling out the pathway dictionary from the pregenerated folder.
pathDict = genPathDict( pullDir, prefix, Core )

# Pulling out a list of all paths.
allPaths = unlistify( pathDict.values() )

# Getting only the unique pathways from this list.
allPaths = np.array( list( { a.tostring(): a for a in allPaths}.values() ) )

NUM_PATHS = 10
# Finding paths for all sets of secreted goods.
for pI, currPath in enumerate(allPaths):
    print_progress_bar(pI, len(allPaths), 'Creating secreted path library')
    currSecVec = secByproducts( currPath.astype(int), rxnMat, prodMat, Core ) * 1
    currSecMets = list( np.nonzero( currSecVec )[0] )

    augNutrientSet = currSecMets

    if not currSecVec.any():
        continue
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
                rMinSubnetsRxns.append(randMinSubnet(satRxnVec, rxnMat, prodMat, sumRxnVec, 
                                                     coreProdRxns, coreTBP, augNutrientSet, Currency))

                # Saving the made pathway.
                np.savetxt(saveDir + str( bool2int(currSecVec) ) + '_prebSeed_randRxnSets_core' + str(coreTBP) + '_' + str(netNum) + '.txt', rMinSubnetsRxns[-1])

