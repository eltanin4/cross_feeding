import numpy as np
from org_generators import *
from sec_byproducts import secByproducts
from load_data import *
from fit_cost import fitCost, biomassCost, compFitness
from display_community import drxn
from tqdm import tqdm
from split_by_demand import splitByDemand

saveDir = 'prebSet_networks/'
secSaveDir = 'sec_prebSet_networks/' # + saveDir
prefix = 'prebSeed_randRxnSets'

medType = 'preb'
mediumFile = medType + '_medium.txt'
nutrientSet = [kegg_to_id[i] for i in np.genfromtxt(mediumFile)]

# Pulling out the pathway dictionary from the pregenerated folder.
pathDict = genPathDict( saveDir, prefix, Core )
secPathDict = genPathDict( secSaveDir, prefix, Core )

# Generating random organism generalist reaction sets.
NUM_RAND_ORGS = 100000

randOrgLenList, randOrgFitList, randOrgBmsList, randOrgComList = np.array( [] ), np.array( [] ), np.array( [] ), np.array( [] )
sandOrgFitList, sandOrgBmsList = np.array( [] ), np.array( [] )
randOrgList = []
numByps = np.array( [] )
bypArr = np.array( [] )
for i in tqdm( range( NUM_RAND_ORGS ) ):
    # Generating a random organism from a pregenerated 
    # path dictionary.
    orgRxns = genRandOrg( pathDict )

    # Calculating the size of the organism and its fitness.
    randOrgLenList = np.append( randOrgLenList, len( orgRxns ) )
    randOrgFitList = np.append( randOrgFitList, fitCost( orgRxns ) )
    # sandE, sandB, isValid = splitByDemand( stoich_matrix, rxnMat, prodMat, 
    #                                        sumRxnVec, rho, pi, nutrientSet, 
    #                                        Energy, Currency, Core, orgRxns )
    
    # if not isValid:
    #     continue

    # sandOrgFitList = np.append( sandOrgFitList, sandE )
    # sandOrgBmsList = np.append( sandOrgBmsList, sandB )
    randOrgBmsList = np.append( randOrgBmsList, biomassCost( orgRxns ) )
    randOrgComList = np.append( randOrgComList, compFitness( orgRxns ) )

    # Keeping a list of generated organisms.
    randOrgList.append( orgRxns )
    
    # Updating the list of byproducts in the network.
    bypArr = np.append( bypArr, np.nonzero( secByproducts( orgRxns, rxnMat, prodMat, Core ) )[0] )
    numByps = np.append( numByps, len( np.nonzero( secByproducts( orgRxns, rxnMat, prodMat, Core ) )[0] ) )

# Generating the smallest generalist characteristics.
smallestOrgLen = min( randOrgLenList[ np.where( randOrgLenList == min( randOrgLenList ) )[0] ] )
fittestOrgLen = min( randOrgLenList[ np.where( randOrgFitList == max( randOrgFitList ) )[0] ] )
maxBmsOrgLen = min( randOrgLenList[ np.where( randOrgBmsList == max( randOrgBmsList ) )[0] ] )

# Generating the fittest generalist characteristics.
smallestOrgFit = max( randOrgFitList[ np.where( randOrgLenList == min( randOrgLenList ) )[0] ] )
fittestOrgFit = max( randOrgFitList[ np.where( randOrgFitList == max( randOrgFitList ) )[0] ] )
maxBmsOrgFit = max( randOrgFitList[ np.where( randOrgBmsList == max( randOrgBmsList ) )[0] ] )

# Generating the most biomass-producing generalist characteristics.
smallestOrgBms = max( randOrgBmsList[ np.where( randOrgLenList == min( randOrgLenList ) )[0] ] )
fittestOrgBms = max( randOrgBmsList[ np.where( randOrgFitList == max( randOrgFitList ) )[0] ] )
maxBmsOrgBms = max( randOrgBmsList[ np.where( randOrgBmsList == max( randOrgBmsList ) )[0] ] )

#-------------------------------------------------------------------------

import matplotlib.pyplot as plt
# Generating the size distribution.
fig, ax = plt.subplots(1)

myWeights = np.ones_like( randOrgLenList ) / len ( randOrgLenList )
ax.hist( randOrgLenList, bins = 25, color = 'darkorange', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Size of metabolic network' )
ax.set_ylabel( 'Fraction of species' )

ax.set_xlim( min( randOrgLenList ) - 10, max( randOrgLenList ) + 10 )
ax.axvline( smallestOrgLen, color='darkorange', linestyle='dashed', linewidth=2.5 )
ax.axvline( fittestOrgLen, color='dodgerblue', linestyle='dashed', linewidth=2.5 )
ax.axvline( maxBmsOrgLen, color='mediumseagreen', linestyle='dashed', linewidth=2.5)
plt.show()

#-------------------------------------------------------------------------


# Generating the energy yield distribution.
fig, ax = plt.subplots(1)

myWeights = np.ones_like( sandOrgFitList ) / len ( sandOrgFitList )
ax.hist( sandOrgFitList, bins = 10, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
# ax.hist( cfnOrgFitList, bins = 25, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Energy yield of metabolic network' )
ax.set_ylabel( 'Fraction of species' )

ax.set_xlim( min( sandOrgFitList ) - 6, max( sandOrgFitList ) + 6 )
ax.axvline( smallestOrgFit, color='darkorange', linestyle='dashed', linewidth=2.5)
ax.axvline( fittestOrgFit, color='dodgerblue', linestyle='dashed', linewidth=2.5)
ax.axvline( maxBmsOrgFit, color='mediumseagreen', linestyle='dashed', linewidth=2.5)

plt.show()

#-------------------------------------------------------------------------

# Generating the biomass yield distribution.
fig, ax = plt.subplots(1)

myWeights = np.ones_like( sandOrgBmsList ) / len ( sandOrgBmsList )
ax.hist( sandOrgBmsList, bins = 5, color = 'mediumseagreen', weights = myWeights, histtype = 'stepfilled' )

ax.set_xlabel( 'Biomass yield of metabolic network' )
ax.set_ylabel( 'Fraction of species' )
ax.set_xlim( min( sandOrgBmsList ) - 2, max( sandOrgBmsList ) + 2 )

# ax.axvline( smallestOrgBms, color='darkorange', linestyle='dashed', linewidth=2.5)
# ax.axvline( maxBmsOrgBms, color='mediumseagreen', linestyle='dashed', linewidth=2.5)
# ax.axvline( fittestOrgBms, color='dodgerblue', linestyle='dashed', linewidth=2.5)

plt.show()

#-------------------------------------------------------------------------
orgNames = []
with open('endo_removed_prok_abbr_kegg.txt', 'r') as f:
    for thisLine in f.readlines():
        orgNames.append( thisLine.strip() )

# Generating the byproduct distribution.
secID, counts = np.unique( bypArr, return_counts = True )

# Removing AMP from this list, have to add it to currency metabolites.
secID, counts = np.delete( secID, np.where( secID == 15 ) ), np.delete( counts, np.where( secID == 15 ) )
secID, counts = np.delete( secID, np.where( secID == 15 ) ), np.delete( counts, np.where( secID == 457 ) )
secID = [ secID for ( counts, secID ) in sorted( zip( counts, secID ), key=lambda pair: pair[0], reverse = True ) ]
counts = [ counts for ( counts, secID ) in sorted( zip( counts, secID ), key=lambda pair: pair[0], reverse = True  ) ]

# Normalizing counts to number of organisms generated.
counts = np.divide( counts, NUM_RAND_ORGS )

# Generating a bar chart.
ind = np.arange( len( secID ) )
width = 0.25
fig, ax = plt.subplots( 1 )
rects1 = ax.bar( ind, counts, width, color='darkorange' )
ax.set_xticks( ind + width / 2 )
ax.set_xticklabels( [ cpd_string_dict[ id_to_kegg [ i ] ] for i in secID ], rotation = 'vertical' )
ax.set_ylabel( 'Frequency of secretion' )
fig.tight_layout()
plt.show()

#-------------------------------------------------------------------------

# Generating the size-yield scatter.
import seaborn as sns
sns.set( style = 'ticks' )
sns.jointplot( randOrgLenList, sandOrgFitList, kind="hex", gridsize=14, color="darkseagreen" ).set_axis_labels('Size of metabolic network', 'Energy yield')
plt.tight_layout()
plt.show()

#-------------------------------------------------------------------------

# Generating the byproduct-yield scatter.
import seaborn as sns
sns.set( style = 'ticks' )
sns.jointplot( numByps, sandOrgFitList, kind = 'hex', gridsize=9, color = 'gray' ).set_axis_labels('Number of byproducts', 'Energy yield of network')
plt.tight_layout()
plt.show()

# Generating the byproduct-yield scatter.
import seaborn as sns
sns.set( style = 'ticks' )
sns.jointplot( numByps, sandOrgBmsList, kind = 'hex', gridsize=9, color = 'gray' ).set_axis_labels('Number of byproducts', 'Biomass yield of network')
plt.tight_layout()
plt.show()

import seaborn as sns
sns.set( style = 'ticks' )
sns.jointplot( numByps, randOrgLenList, kind = 'hex', gridsize=9, color = 'gray' ).set_axis_labels('Number of byproducts', 'Size of network')
plt.tight_layout()
plt.show()
