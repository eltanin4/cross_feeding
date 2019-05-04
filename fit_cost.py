from load_data import *

medType = 'preb'
mediumFile = medType + '_medium.txt'
nutrientSet = [kegg_to_id[i] for i in np.genfromtxt(mediumFile)]
# nutrientSet = [25, 12]

def fitCost( rxns, currNutrientSet=nutrientSet ):
    mulFac = [1, 3, 3]
    if nutrientsConsumed( rxns, currNutrientSet ) == 0.0:
        return 0.0
    return int( sum( [ sum( stoich_matrix[ rxn ][ Energy ] * mulFac ) 
                for rxn in rxns.astype(int) ] ) ) / nutrientsConsumed( rxns, currNutrientSet )

def nutrientsConsumed( rxns, currNutrientSet=nutrientSet ):
    return abs( sum( [ sum( stoich_matrix[ rxn ][ nutrientSet ] ) for rxn in rxns.astype(int)
                       if sum( stoich_matrix[ rxn ][ currNutrientSet ] ) < 0.0 ] ) )

def biomassCost( rxns, currNutrientSet=nutrientSet ):
    if nutrientsConsumed( rxns, currNutrientSet ):
        return int( sum( [ sum( stoich_matrix[ rxn ][ Core ] ) 
                    for rxn in rxns.astype(int) ] ) ) / nutrientsConsumed( rxns, currNutrientSet )
    else:
        return 0.0

def compFitness( orgRxns, currNutrientSet=nutrientSet, alpha=1.0, beta=0.0, gamma=0.10 ):
    return ( alpha * fitCost( orgRxns, currNutrientSet ) 
             + beta * biomassCost( orgRxns, currNutrientSet )
             - gamma * len( orgRxns ) )

def redoxFitCost( rxns, currExchangeRatio=3, currNutrientSet=nutrientSet ):
    mulFac = [1, currExchangeRatio, currExchangeRatio]
    return int( sum( [ sum( stoich_matrix[ rxn ][ Energy ] * mulFac ) 
                for rxn in rxns.astype(int) ] ) ) / nutrientsConsumed( rxns, currNutrientSet )
