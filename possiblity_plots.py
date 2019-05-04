import numpy as np
import matplotlib.pyplot as plt
from org_generators import *
from sec_byproducts import secByproducts
import itertools

saveDir = 'richSet_networks/'
secSaveDir = 'sec_' + saveDir
prefix = 'prebSeed_randRxnSets'

# Pulling out the pathway dictionary from the pregenerated folder.
pathDict = genPathDict( saveDir, prefix, Core )
secPathDict = genPathDict( secSaveDir, prefix, Core )

ratioDict = { }
for coreTBP in Core:
    ratioDict[ coreTBP ] = []
    for currPair in itertools.product( 
                    pathDict[ coreTBP ], secPathDict[ coreTBP ] ):
        ratioDict[ coreTBP ].append( fitCost( currPair[ 1 ] ) - fitCost( currPair[ 0 ] ) )


counts = [ np.mean( ratioDict[ coreTBP ] ) for coreTBP in Core ]
yerrs = [ np.std( ratioDict[ coreTBP ] ) for coreTBP in Core ]
# yerrs = [ counts[ i ] - min( ratioDict[ coreTBP ] ) for i, coreTBP in enumerate( Core ) ]
coresOrd = [ coresOrd for ( counts, coresOrd ) in sorted( zip( counts, Core ), key=lambda pair: pair[0], reverse = True ) ]
yerrs = [ yerrd for ( counts, yerrd ) in sorted( zip( counts, yerrs ), key=lambda pair: pair[0], reverse = True ) ]
counts = [ counts for ( counts, coresOrd ) in sorted( zip( counts, Core ), key=lambda pair: pair[0], reverse = True  ) ]


ind = np.arange( len( Core ) )
width = 0.25
fig, ax = plt.subplots( 1 )
rects1 = ax.bar( ind, counts, width, color='darkorange', yerr=yerrs, 
                 error_kw=dict(ecolor='black', lw=1.2, capsize=3, capthick=1.2) )
ax.set_xticks( ind + width / 2 )
ax.set_xticklabels( [ cpd_string_dict[ id_to_kegg [ i ] ] for i in coresOrd ], rotation = 'vertical' )
ax.set_ylabel( 'Energy yield benefit with secretions ' )
ax.axhline( 0.0, color='black', linestyle='dashed', linewidth=1.5)
# ax.set_ylim(0.0, 9.0)
fig.tight_layout()
plt.savefig('plots_170328/fitnComp_byp_richMed.svg')
plt.savefig('plots_170328/fitnComp_byp_richMed.png')
plt.show()
