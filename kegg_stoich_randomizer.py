import numpy as np
import random
from load_data import *
from tqdm import tqdm
from copy import deepcopy

NUM_SWAPS = 20000

def give_random_reaction_id( stoich_matrix ):
    return random.choice( range( np.shape( stoich_matrix )[0] ) )

def swap_once( curr_stoich_matrix, Currency, prob_reactant=0.5 ):
    rID1, rID2 = ( give_random_reaction_id( curr_stoich_matrix ), 
                   give_random_reaction_id( curr_stoich_matrix ) )
    
    if random.random() <  prob_reactant:
        r1Rs = np.where( curr_stoich_matrix[ rID1 ] < 0 )[0]
        r2Rs = np.where( curr_stoich_matrix[ rID2 ] < 0 )[0]
    else:
        r1Rs = np.where( curr_stoich_matrix[ rID1 ] > 0 )[0]
        r2Rs = np.where( curr_stoich_matrix[ rID2 ] > 0 )[0]

    # Removing currency metabolites from available set.
    r1Rs_currRemoved = [ e for e in r1Rs ]
    r2Rs_currRemoved = [ e for e in r2Rs ]

    # Checking if only currency metabolites are involved.
    if not r1Rs_currRemoved or not r2Rs_currRemoved:
        return curr_stoich_matrix

    # Choosing swapping reactants/products from available set.
    toSwapR1, toSwapR2 = ( random.choice( r1Rs_currRemoved ), 
                           random.choice( r2Rs_currRemoved ) )

    # Saving stoichiometry of the chosen reactant/products.
    stoichOf_r1R, stoichOf_r2R = ( curr_stoich_matrix[ rID1, toSwapR1 ], 
                                   curr_stoich_matrix[ rID2, toSwapR2 ] )

    # Swapping in the current version of the stoichiometric matrix.
    if not ( curr_stoich_matrix[ rID1, toSwapR2 ] 
             or curr_stoich_matrix[ rID2, toSwapR1 ] ):
        curr_stoich_matrix[ rID1, toSwapR1 ] = 0
        curr_stoich_matrix[ rID2, toSwapR2 ] = 0
        curr_stoich_matrix[ rID1, toSwapR2 ] = stoichOf_r2R
        curr_stoich_matrix[ rID2, toSwapR1 ] = stoichOf_r1R

    return curr_stoich_matrix

curr_stoich_matrix = deepcopy( stoich_matrix[ :, : ] )
for this_swap_num in tqdm( range( NUM_SWAPS ) ):
    curr_stoich_matrix = swap_once( curr_stoich_matrix, Currency, prob_reactant=0.5 )
