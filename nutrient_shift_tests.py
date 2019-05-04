import numpy as np
from tqdm import tqdm
from sat_marker import markSatMetsRxns
from fit_cost import fitCost

def aut_getNutrientShiftResults( randOrgList, nutrientSet, rxnMat, prodMat, 
                                 sumRxnVec, Core, Currency ):
    """
    Gets a set of organisms and a set of nutrients and performs random nutrient shifts on each. Marks all achievable molecules on this shifted nutrient set, and tests if the organism survives, and when it does, does its yield increase or decrease.

    ARGUMENTS:
    randOrgList (list): set of reaction indices for each organism's metabolic network.
    nutrientSet (list): custom IDs of nutrients originally in the environment.
    rxnMat (np.array): matrix of shape (R x M) that lists reactants in reactions.
    prodMat (np.array): matrix of shape(R x M) that lists products in reactions.
    sumRxnVec (np.array): vector of size R that lists the number of reactants in each reaction.
    Core (list): custom IDs of core molecules.
    Currency (list): custom IDs of currency molecules.

    RETURNS:
    notSurvive (int): number of shift cases where the organism cannot make all core molecules.
    yieldInc (int): number of shift cases where the organism's yield increases.
    yieldDec (int): number of shift cases where the organism's yield decreases.
    numCoresActive (list): number of core molecules that could still be made after shift if dead.
    """
    yieldDec, yieldInc, notSurvive = 0, 0, 0
    numCoresActive = []
    for orgRxns in tqdm( randOrgList[ 6000 : 7000] ):
        # Get a vector of reactants in the current network.
        reactVec = np.logical_or.reduce( rxnMat[ orgRxns ] ) * 1
        allowedNutrients = np.where( reactVec == 1 )[0]

        # Nutrient shifting pipeline.
        randShiftedNutrient = np.random.choice( [ n for n in allowedNutrients 
                                                  if n not in nutrientSet ] )
        shiftedNutrientSet = list( np.append( np.random.choice(nutrientSet),
                                              randShiftedNutrient ) )

        # Testing for survival on this nutrient set.
        satRxns = orgRxns[:]
        currSatRxnVec = np.zeros( len( rxnMat ) )
        currSatRxnVec[ satRxns ] = 1
        tempSatRxnVec = np.copy( currSatRxnVec )
        tempSatMetVec, tempSatRxnVec = markSatMetsRxns( tempSatRxnVec, rxnMat, prodMat, 
                                                        sumRxnVec, [], shiftedNutrientSet, 
                                                        Currency )

        # Survival guaranteed.
        if tempSatMetVec[Core].all():
            # Checking if yield increases or decreases.
            if ( fitCost( np.nonzero( tempSatRxnVec )[0], 
                          shiftedNutrientSet ) < fitCost( orgRxns, nutrientSet ) ):
                yieldDec += 1
            else:
                yieldInc += 1
        else:
            notSurvive += 1
            numCoresActive.append( np.count_nonzero( tempSatMetVec[ Core ] ) )

    return notSurvive, yieldInc, yieldDec, numCoresActive

#--------------------------------------------------------------------------------------------------#--------------------------------------------------------------------------------------------------

def cfn_getNutrientShiftResults( randFitterSet, secList, nutrientSet, 
                                 rxnMat, prodMat, sumRxnVec, Core, Currency ):
    """
    Gets a set of pairs of cross-feeding organisms and a set of nutrients and performs random nutrient shifts on each. Marks all achievable molecules on this shifted nutrient set, and tests if both organisms survive, and when they do, does their yield increase or decrease.

    ARGUMENTS:
    randFitterSet (list): set of pair of reaction indices for each pair of organisms' metabolic networks.
    secList (list): set of pairs of custom IDs of metabolites secreted by first to second and vice versa.
    nutrientSet (list): custom IDs of nutrients originally in the environment.
    rxnMat (np.array): matrix of shape (R x M) that lists reactants in reactions.
    prodMat (np.array): matrix of shape(R x M) that lists products in reactions.
    sumRxnVec (np.array): vector of size R that lists the number of reactants in each reaction.
    Core (list): custom IDs of core molecules.
    Currency (list): custom IDs of currency molecules.

    RETURNS:
    notSurvive (int): number of shift cases where both organisms cannot make all core molecules.
    notSec (int): number of shift cases where the secretions where lost.
    yieldInc (int): number of shift cases where both organisms yield increases.
    yieldDec (int): number of shift cases where both organisms yield decreases.
    numCoresActive (list): number of core molecules that could still be made after shift if dead.
    """
    yieldDec, yieldInc, notSurvive, notSec = 0, 0, 0, 0
    numCoresActive = []
    thisIter = -1
    for ( o1Rxns, o2Rxns ) in tqdm( randFitterSet ):
        thisIter += 1
        # Get a vector of reactants present in both networks.
        reactVec = np.logical_or( np.logical_or.reduce( rxnMat[ o1Rxns ] ) * 1,
                                   np.logical_or.reduce( rxnMat[ o2Rxns ] ) * 1 )
        allowedNutrients = np.where( reactVec == 1 )[0]

        # Nutrient shifting pipeline.
        randShiftedNutrient = np.random.choice( [ n for n in allowedNutrients 
                                                  if n not in nutrientSet ] )
        shiftedNutrientSet = list( np.append( np.random.choice(nutrientSet),
                                              randShiftedNutrient ) )

        # Testing for survival on this nutrient set.
        sat1Rxns, sat2Rxns = o1Rxns[:], o2Rxns[:]
        currSatRxnVec1, currSatRxnVec2 = np.zeros( len( rxnMat ) ), np.zeros( len( rxnMat ) )
        currSatRxnVec1[ sat1Rxns ], currSatRxnVec2[ sat2Rxns ] = 1, 1
        tempSatRxnVec1, tempSatRxnVec2 = np.copy( currSatRxnVec1 ), np.copy( currSatRxnVec2 )
        tempSatMetVec1, tempSatRxnVec1 = markSatMetsRxns( tempSatRxnVec1, rxnMat, prodMat, 
                                                          sumRxnVec, [], shiftedNutrientSet, 
                                                          Currency )
        tempSatMetVec2, tempSatRxnVec2 = markSatMetsRxns( tempSatRxnVec2, rxnMat, prodMat, 
                                                          sumRxnVec, [], shiftedNutrientSet, 
                                                          Currency )
        # Survival guaranteed.
        if tempSatMetVec1[ Core ].all() and tempSatMetVec2[ Core ].all():

            if ( tempSatMetVec1[ secList[ thisIter ][ 0 ] ].all() and 
                 tempSatMetVec2[ secList[ thisIter ][ 1 ] ].all() ):

                # Checking if yield increases or decreases.
                if ( fitCost( np.nonzero( tempSatRxnVec1 )[0], 
                              shiftedNutrientSet ) < fitCost( o1Rxns, nutrientSet )  
                    and
                   ( fitCost( np.nonzero( tempSatRxnVec2 )[0], 
                              shiftedNutrientSet ) < fitCost( o2Rxns, nutrientSet ) ) ):
                    yieldDec += 1
                else:
                    yieldInc += 1
            else:
                notSec += 1
                numCoresActive.append( np.count_nonzero( tempSatMetVec1[ Core ] ) - 1 )
                numCoresActive.append( np.count_nonzero( tempSatMetVec2[ Core ] ) - 1)

        else:
            notSurvive += 1
            numCoresActive.append( np.count_nonzero( tempSatMetVec1[ Core ] ) )
            numCoresActive.append( np.count_nonzero( tempSatMetVec2[ Core ] ) )

    return notSurvive, notSec, yieldInc, yieldDec, numCoresActive

