from sat_marker import markSatMetsRxns

allowedNutrients = np.where(np.logical_or.reduce(rxnMat[necoliRxns]))[0]
viableNutrientSets = []
SIZE_OF_NUTRIENT_SET = 2

while not viableNutrientSets:
    SIZE_OF_NUTRIENT_SET += 1
    for thisNutSet in tqdm(list(itertools.combinations(allowedNutrients, SIZE_OF_NUTRIENT_SET))):
        # Testing for survival on this nutrient set.
        satRxns = necoliRxns[:]
        currSatRxnVec = np.zeros( len( rxnMat ) )
        currSatRxnVec[ satRxns ] = 1
        tempSatRxnVec = np.copy( currSatRxnVec )
        tempSatMetVec, tempSatRxnVec = markSatMetsRxns( tempSatRxnVec, rxnMat, prodMat, 
                                                        sumRxnVec, [], list(thisNutSet), 
                                                        Currency )

        # Survival guaranteed.
        if np.sum(tempSatMetVec[Core]) == len(Core) - 1:
            # Checking if yield increases or decreases.
            viableNutrientSets.append( thisNutSet )

#---------------------------------------------------------------------------
# Results
#---------------------------------------------------------------------------
# [17, 19, 42, 69]
# [17, 19, 69, 228]
# [17, 20, 69, 228]
# [17, 69, 133, 228]
# [17, 69, 228, 458]
# [19, 42, 69, 81]
# [19, 69, 81, 228]
# [20, 69, 81, 228]
# [69, 81, 133, 228]
# [69, 81, 228, 45]

prunedRxns, prunedMets = [], []
for thisNutSet in viableNutrientSets:
    # Testing for survival on this nutrient set.
    satRxns = necoliRxns[:]
    currSatRxnVec = np.zeros( len( rxnMat ) )
    currSatRxnVec[ satRxns ] = 1
    tempSatRxnVec = np.copy( currSatRxnVec )
    tempSatMetVec, tempSatRxnVec = markSatMetsRxns( tempSatRxnVec, rxnMat, prodMat, 
                                                    sumRxnVec, [], list(thisNutSet), 
                                                    Currency )
    prunedRxns.append( tempSatRxnVec )
    prunedMets.append( tempSatMetVec )
