medType = 'preb'
mediumFile = medType + '_medium.txt'

import numpy as np
import pickle
from display_community import drxn
from give_reverse_scope import giveRevScope
# from smallest_min_subgraph import smallestMinSubnets
# from fittest_minimal_subgraph import fittestMinSubnets
from random_minimal_subgraph import randMinSubnet
# from rand_withfit_min_subgraph import randWFitnessMinSubnet
from new_load_data import *
import os

#-------------------------------------------------------------------------
# Beginning to load things to the network expansion function.
#-------------------------------------------------------------------------
# Define the seed set, this would the nutrients in the medium (see line 1)
nutrientSet = [kegg_to_id[i] for i in np.genfromtxt(mediumFile)]

saveDir = 'newCore_' + medType + 'Set_networks/'

if not os.path.exists( saveDir ):
    os.mkdir( saveDir )

print( Core )

for coreTBP in Core:
    print('Finding paths to ' + cpd_string_dict[id_to_kegg[coreTBP]] + '.')
    coreProdRxns = (prodMat[:, coreTBP] == 1) * 1
    satMetVec, satRxnVec = giveRevScope(rxnMat, prodMat, sumRxnVec, 
                                        nutrientSet, Currency, coreTBP)
    satMets, satRxns = np.nonzero(satMetVec)[0], np.nonzero(satRxnVec)[0]

    rMinSubnetsRxns = []
    for netNum in range(10):
        print('Constructing pathway number ' + str( netNum ) )
        rMinSubnetsRxns.append(randMinSubnet(satRxnVec, rxnMat, prodMat, sumRxnVec, 
                                             coreProdRxns, coreTBP, nutrientSet, Currency))
        np.savetxt(saveDir + 'prebSeed_randRxnSets_core' + str(coreTBP) + '_' + str(netNum) + '.txt', np.array(list(rMinSubnetsRxns[-1])))

    print(coreTBP, '\t', cpd_string_dict[id_to_kegg[coreTBP]], '\t')

#-------------------------------------------------------------------------
# Calculating paths using only the secreted nutrients.
#-------------------------------------------------------------------------
augNutrientSet = list(np.array(secID).astype(int))
saveDir = 'newCore_sec_' + medType + 'Set_networks/'
for coreTBP in Core:
 print('Finding paths to ' + cpd_string_dict[id_to_kegg[coreTBP]] + '.')
 coreProdRxns = (prodMat[:, coreTBP] == 1) * 1
 satMetVec, satRxnVec = giveRevScope(rxnMat, prodMat, sumRxnVec, 
                                     augNutrientSet, Currency, coreTBP)
 satMets, satRxns = np.nonzero(satMetVec)[0], np.nonzero(satRxnVec)[0]

 rMinSubnetsRxns = []
 for netNum in range(10):
     print('Constructing pathway number ' + str( netNum ) )
     rMinSubnetsRxns.append(randMinSubnet(satRxnVec, rxnMat, prodMat, sumRxnVec, 
                                          coreProdRxns, coreTBP, augNutrientSet, Currency))
     np.savetxt(saveDir + 'prebSeed_randRxnSets_core' + str(coreTBP) + '_' + str(netNum) + '.txt', np.array(list(rMinSubnetsRxns[-1])))

 print(coreTBP, '\t', cpd_string_dict[id_to_kegg[coreTBP]], '\t')

