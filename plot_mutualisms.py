import pickle
import matplotlib.pyplot as plt
from unlistify import unlistify

fitterDB = pickle.load( open( 'fitterDB_prebMed.dat', 'rb' ) )

unionDB = [ np.array( uniqify( np.append( e[0], e[1] ) ) ) 
            for e in allDB ]

# Generating the biomass yield distribution of the fitter organisms.
mutualBiomasses = unlistify( [ [ biomassCost( e[0] ), biomassCost( e[1] ) ] 
                              for e in fitterDB ] )

fig, ax = plt.subplots(1)

myWeights = np.ones_like( mutualBiomasses ) / len ( mutualBiomasses )
ax.hist( mutualBiomasses, bins = 10, color = 'mediumseagreen', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Biomass yield of metabolic network' )
ax.set_ylabel( 'Fraction of species' )

ax.set_xlim( -2, max( mutualBiomasses ) + 1 )
ax.axvline( smallestOrgBms, color='darkorange', linestyle='dashed', linewidth=2.5)
ax.axvline( fittestOrgBms, color='dodgerblue', linestyle='dashed', linewidth=2.5)
ax.axvline( maxBmsOrgBms, color='mediumseagreen', linestyle='dashed', linewidth=2.5)

plt.savefig('plots_170422/mutual_met_biom_prebMed.svg')
plt.savefig('plots_170422/mutual_met_biom_prebMed.png')
plt.show()

#-------------------------------------------------------------------------

fig, ax = plt.subplots(1)

myWeights = np.ones_like( unionBiomasses ) / len ( unionBiomasses )
ax.hist( unionBiomasses, bins = 10, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
myWeights = np.ones_like( randOrgBmsList ) / len ( randOrgBmsList )
ax.hist( randOrgBmsList, bins = 10, color = 'mediumseagreen', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Biomass yield of metabolic network' )
ax.set_ylabel( 'Fraction of species' )

ax.set_xlim( -2, max( unionBiomasses ) + 1 )
# ax.axvline( smallestOrgBms, color='darkorange', linestyle='dashed', linewidth=2.5)
# ax.axvline( fittestOrgBms, color='dodgerblue', linestyle='dashed', linewidth=2.5)
ax.axvline( maxBmsOrgBms, color='mediumseagreen', linestyle='dashed', linewidth=2.5)

plt.savefig('plots_180226/mutual_met_biom_prebMed.svg')
plt.savefig('plots_180226/mutual_met_biom_prebMed.png')
plt.show()

#-------------------------------------------------------------------------

# Generating the energy yield distribution of the fitter organisms.
mutualFitnesses = unlistify( [ [ fitCost( e[0] ), fitCost( e[1] ) ] 
                              for e in fitterDB ] )

fig, ax = plt.subplots(1)

myWeights = np.ones_like( mutualFitnesses ) / len ( mutualFitnesses )
ax.hist( mutualFitnesses, bins = 10, color = 'dodgerblue', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Energy yield of metabolic network' )
ax.set_ylabel( 'Fraction of species' )

ax.set_xlim( -8, max( mutualFitnesses ) + 2 )
ax.axvline( smallestOrgFit, color='darkorange', linestyle='dashed', linewidth=2.5)
ax.axvline( fittestOrgFit, color='dodgerblue', linestyle='dashed', linewidth=2.5)
ax.axvline( maxBmsOrgFit, color='mediumseagreen', linestyle='dashed', linewidth=2.5)

plt.savefig('plots_170422/mutual_met_fitn_prebMed.svg')
plt.savefig('plots_170422/mutual_met_fitn_prebMed.png')
plt.show()

#-------------------------------------------------------------------------

fig, ax = plt.subplots(1)

myWeights = np.ones_like( unionFitnesses ) / len ( unionFitnesses )
myWeights2 = np.ones_like( allFitnesses ) / len ( allFitnesses )
ax.hist( allFitnesses, bins = 10, color = 'dodgerblue', weights = myWeights2, histtype = 'stepfilled' )
ax.hist( unionFitnesses, bins = 7, color = 'mediumseagreen', weights = myWeights, histtype = 'stepfilled' )
plt.show()

myWeights = np.ones_like( randOrgFitList ) / len ( randOrgFitList )
ax.hist( randOrgFitList, bins = 10, color = 'mediumseagreen', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Energy yield of metabolic network' )
ax.set_ylabel( 'Fraction of species' )

ax.set_xlim( -8, max( unionFitnesses ) + 2 )
# ax.axvline( smallestOrgFit, color='darkorange', linestyle='dashed', linewidth=2.5)
ax.axvline( fittestOrgFit, color='dodgerblue', linestyle='dashed', linewidth=2.5)
# ax.axvline( maxBmsOrgFit, color='mediumseagreen', linestyle='dashed', linewidth=2.5)

plt.savefig('plots_180226/mutual_met_fitn_prebMed.svg')
plt.savefig('plots_180226/mutual_met_fitn_prebMed.png')
plt.show()

#-------------------------------------------------------------------------

# Generating the size distribution of the fitter organisms.
mutualSizes = unlistify([ [ len( e[0] ), len( e[1] ) ]
                              for e in fitterDB  ])

fig, ax = plt.subplots(1)

myWeights = np.ones_like( mutualSizes ) / len ( mutualSizes )
ax.hist( mutualSizes, bins = 20, color = 'darkorange', weights = myWeights, histtype = 'stepfilled' )
ax.set_xlabel( 'Size of metabolic network' )
ax.set_ylabel( 'Fraction of species' )

ax.set_xlim( min( mutualSizes ) - 30, max( mutualSizes ) + 10 )
ax.axvline( smallestOrgLen, color='darkorange', linestyle='dashed', linewidth=2.5)
ax.axvline( fittestOrgLen, color='dodgerblue', linestyle='dashed', linewidth=2.5)
ax.axvline( maxBmsOrgLen, color='mediumseagreen', linestyle='dashed', linewidth=2.5)

plt.savefig('plots_170422/mutual_met_size_prebMed.svg')
plt.savefig('plots_170422/mutual_met_size_prebMed.png')
plt.show()

#-------------------------------------------------------------------------

# Getting a histogram of shared reactions.
numRxnsShared = [ len( set( e[0] ).intersection( set(e[1]) ) ) / 
                  len( set( e[0] ).union( set (e[1]) ) )
                  for e in fitterDB ] 

fig, ax = plt.subplots(1)
myWeights = np.ones_like( numRxnsShared ) / len ( numRxnsShared )
ax.hist( numRxnsShared, bins = 7, color = 'gray', weights = myWeights)
ax.set_xlabel( 'Fraction of reactions shared in mutualistic pairs' )
ax.set_ylabel( 'Frequency of cases' )

ax.set_xlim( 0, 1 )
fig.tight_layout()
plt.savefig('plots_170422/mutual_frac_shared_prebMed.svg')
plt.savefig('plots_170422/mutual_frac_shared_prebMed.png')
plt.show()

#-------------------------------------------------------------------------

# Control with all pairs in the set.
allPairs = [ (allOrgs[2 * i], allOrgs[2 * i + 1]) for i in range( int( len( allOrgs ) / 2 ) ) ]

pnumRxnsShared = [ len( set( e[0] ).intersection( set(e[1]) ) ) / 
                  len( set( e[0] ).union( set (e[1]) ) )
                  for e in allPairs ] 

fig, ax = plt.subplots(1)
myWeights = np.ones_like( pnumRxnsShared ) / len ( pnumRxnsShared )
ax.hist( pnumRxnsShared, bins = 15, color = 'gray', weights = myWeights)
ax.set_xlabel( 'Fraction of reactions shared in random pairs' )
ax.set_ylabel( 'Frequency of cases' )

ax.set_xlim( 0, 1 )
fig.tight_layout()
plt.savefig('plots_170422/control_mutual_frac_shared_prebMed.svg')
plt.savefig('plots_170422/control_mutual_frac_shared_prebMed.png')
plt.show()


#-------------------------------------------------------------------------

# Getting a scatter for reactions shared versus increase in fitness.
numRxnsShared = np.array([ len( set( e[0] ).intersection( set(e[1]) ) ) / 
                           len( set( e[0] ).union( set (e[1]) ) )
                           for e in fitterDB ])
allFits = np.array( [ np.mean([mutualFitnesses[2 * i], mutualFitnesses[2 * i + 1]]) for i in range( int( len( mutualFitnesses ) / 2 ) ) ] ) - fittestOrgFit

import seaborn as sns
sns.set( style = 'ticks' )
sns.jointplot( np.array(numRxnsShared), allFits, kind = 'reg', color = 'gray' ).set_axis_labels('Fraction of reactions shared in mutualistic pair', 'Average energy yield benefit')
plt.tight_layout()
plt.savefig('plots_170401/corr_overlap_benefit_prebMed.svg')
plt.savefig('plots_170401/corr_overlap_benefit_prebMed.png')
plt.show()


#-------------------------------------------------------------------------------
# Creating pandas dataframe and organizing genome size data.
#-------------------------------------------------------------------------------
nut_removal_dataset = []
nut_removal_dataset += [ [ e, 'Autonomous networks' ] for e in nAUT ]
nut_removal_dataset += [ [ e, 'Cross-feeding networks' ] for e in nCFN ]

nut_removal_dataset = pd.DataFrame(nut_removal_dataset, index=range(len(nut_removal_dataset)),
                                columns=['Biomass precursors not producible', 'Network types'])

#-------------------------------------------------------------------------------
# Plotting violin plot from the genome data.
#-------------------------------------------------------------------------------
import seaborn as sns
sns.set(style='white', palette='Set2', font='Helvetica')
fig, ax = plt.subplots(1)
ax.set_title('Comparison of genome sizes between core and peripheral species')
sns.violinplot(x='Network types', y='Biomass precursors not producible', data=nut_removal_dataset,
               scale="width", lw=1.0, inner='quartile')
# sns.stripplot(x='Network types', y='Biomass precursors not producible', data=rxn_removal_dataset,
#               jitter=False, color="black", edgecolor="white")
ax.set_ylim(0.0)
plt.savefig('plots_170430/' + 'distfig_nut.svg')
plt.show()
