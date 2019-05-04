yDiffs = []
NUM_TO_DO = 5000
num_done = 0
while num_done < NUM_TO_DO:
    print_progress_bar(num_done, NUM_TO_DO, 'Building mutant dictionary')
    while True:
        orgPathDict = {}
        for coreTBP in Core:
            orgPathDict[ coreTBP ] = list( random.choice( pathDict[ coreTBP ] ).astype(int) )
        
        # Storing the list of reactions.
        orgRxns = np.array( uniqify( unlistify( orgPathDict.values() ) ) ).astype( int )

        # Creating a dictionary of secretions.
        orgSecDict = {}
        for coreTBP in Core:
            orgSecDict[ coreTBP ] = list( np.nonzero( pathwaySecByproducts( 
                                          orgPathDict[ coreTBP ], orgRxns, 
                                          rxnMat, prodMat, Core ) )[0] )
        
        # Purging duplicates to have only unique byproducts.
        tempSet = set()
        fullSet = unlistify( orgSecDict.values() )
        duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
        orgSecDict = { coreTBP: list( set( orgSecDict[ coreTBP ] ).difference( duplicates ) ) 
                       for coreTBP in Core }

        # Updating the full list of byproducts in the network.
        bypArr = uniqify( unlistify( orgSecDict.values() ) )

        if newFitness( orgRxns ) > 0.0 and bypArr:
            # print('Successfuly found a working host.')
            ogCost = newFitness( orgRxns )
            break
    
    # Calculating byproduct yields.    
    E, B, isFull, bypYields = splitDemandByps(orgRxns, bypArr)
    
    # Creating mutant knock-out networks.
    muts = list(itertools.combinations(list(orgPathDict.values()), len(Core) - 1))
    for thisMut in muts:
        mutRxns = np.array( uniqify( unlistify( thisMut ) ) ).astype( int )
        mE, mB, misFull, mbypYields = splitDemandByps(mutRxns, bypArr)

        # Calculating yield differences.
        yDiffs.append([mbypYields[tb] - bypYields[tb] 
                       for tb in bypArr 
                       if bypYields[tb]])

    num_done += 1
