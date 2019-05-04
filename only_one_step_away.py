import timeit
import pandas as pd

saveDir = 'randKEGG_prebSet_networks/'
secSaveDir = 'sec_' + saveDir
prefix = 'prebSeed_randRxnSets'

# Pulling out the pathway dictionary from the pregenerated folder.
pathDict = genPathDict( saveDir, prefix, Core )
secPathDict = genPathDict( secSaveDir, prefix, Core )

# Generate a random organism.
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

    if fitCost( orgRxns ) > 0.0:
        print('Successfuly found a working host.')
        ogCost = fitCost( orgRxns )
        break

#-------------------------------------------------------------------------

# Generate a random organism.
while True:
    org2PathDict = {}
    for coreTBP in Core:
        org2PathDict[ coreTBP ] = list( random.choice( pathDict[ coreTBP ] ).astype(int) )
    
    # Storing the list of reactions.
    org2Rxns = np.array( uniqify( unlistify( org2PathDict.values() ) ) ).astype( int )

    # Creating a dictionary of secretions.
    org2SecDict = {}
    for coreTBP in Core:
        org2SecDict[ coreTBP ] = list( np.nonzero( pathwaySecByproducts( 
                                      org2PathDict[ coreTBP ], org2Rxns, 
                                      rxnMat, prodMat, Core ) )[0] )
    
    # Purging duplicates to have only unique byproducts.
    tempSet = set()
    fullSet = unlistify( org2SecDict.values() )
    duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
    org2SecDict = { coreTBP: list( set( org2SecDict[ coreTBP ] ).difference( duplicates ) ) 
                   for coreTBP in Core }

    # Updating the full list of byproducts in the network.
    bypArr = uniqify( unlistify( org2SecDict.values() ) )

    if fitCost( org2Rxns ) > 0.0:
        print('Successfuly found a working host.')
        ogCost = fitCost( org2Rxns )
        break

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------


# Creating a list of pathways that work on the respective secretions.
usablePaths = {}
for coreTBP in Core:
    usablePaths[ coreTBP ] = []
    tempBypList = orgSecDict[ coreTBP ][:]
    if not tempBypList:
        continue
    for currPath in secPathDict[ coreTBP ]:
        if ((np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)
            [ tempBypList ].any()):
            usablePaths[ coreTBP ].append( currPath )

#-------------------------------------------------------------------------

for coreTBP in Core:
    if not usablePaths[ coreTBP ]:
        continue

#-------------------------------------------------------------------------

    # Randomly picking a usable path.
    PATH_ID = random.choice( range( len( usablePaths[ coreTBP ] ) ) )
    tempOrgPathDict = org2PathDict.copy()
    tempOrgPathDict[ coreTBP ] = list( usablePaths[ coreTBP ][ PATH_ID ].astype( int ) )
    tempOrgRxns = np.array( uniqify( unlistify( tempOrgPathDict.values() ) ) ).astype( int )
    replacedPathSec = list( np.nonzero( pathwaySecByproducts( 
                            tempOrgPathDict[ coreTBP ], tempOrgRxns, 
                            rxnMat, prodMat, Core ) )[0] )

    # Creating a dictionary of secretions.
    tempOrgSecDict = {}
    for i in Core:
        tempOrgSecDict[ i ] = list( np.nonzero( pathwaySecByproducts( 
                                      tempOrgPathDict[ i ], tempOrgRxns, 
                                      rxnMat, prodMat, Core ) )[0] )
    
    # Purging duplicates to have only unique byproducts.
    tempSet = set()
    fullSet = unlistify( tempOrgSecDict.values() )
    duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
    tempOrgSecDict = { i: list( set( tempOrgSecDict[ i ] ).difference( duplicates ) ) 
                   for i in Core }

#-------------------------------------------------------------------------

    if replacedPathSec:
        newUsablePaths = {}
        for i in Core:
            if i != coreTBP:
                newUsablePaths[ i ] = []
                for currPath in secPathDict[ i ]:
                    if ((np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)
                        [ replacedPathSec ].any()):
                        newUsablePaths[ i ].append( currPath )

#-------------------------------------------------------------------------

    for i in Core:
        if newUsablePaths[ i ]:
            PATH_ID = random.choice( range ( len( newUsablePaths[ i ] ) ) )
            newOrgPathDict = orgPathDict.copy()
            newOrgPathDict[ i ] = list( newUsablePaths[ i ][ PATH_ID ].astype( int ) )
            newOrgRxns = np.array( uniqify( unlistify( newOrgPathDict.values() ) ) ).astype( int )
            break

    # Creating a dictionary of secretions.
    newOrgSecDict = {}
    for i in Core:
        newOrgSecDict[ i ] = list( np.nonzero( pathwaySecByproducts( 
                                      newOrgPathDict[ i ], newOrgRxns, 
                                      rxnMat, prodMat, Core ) )[0] )
    
    # Purging duplicates to have only unique byproducts.
    tempSet = set()
    fullSet = unlistify( newOrgSecDict.values() )
    duplicates = set(x for x in fullSet if x in tempSet or tempSet.add(x))
    newOrgSecDict = { i: list( set( newOrgSecDict[ i ] ).difference( duplicates ) ) 
                   for i in Core }

#-------------------------------------------------------------------------

while True:
    randPath = random.choice( pathDict[ coreTBP ] )
    if not (np.logical_or.reduce(prodMat[np.array( randPath ).astype(int)]) * 1)[68]:
        break

tempOrgNutPathDict = {}
for coreTBP in Core:
    currPath = tempOrgPathDict[ coreTBP ]
    tempOrgNutPathDict[ coreTBP ] = [ nutrient for nutrient in nutrientSet 
                                      if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ nutrient ] ]

tempOrgSecPathDict = {}
for coreTBP in Core:
    currPath = tempOrgPathDict[ coreTBP ]
    tempOrgSecPathDict[ coreTBP ] = [ sec for sec in newOrgSecs 
                                      if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ sec ] ]

newOrgNutPathDict = {}
for coreTBP in Core:
    currPath = newOrgPathDict[ coreTBP ]
    newOrgNutPathDict[ coreTBP ] = [ nutrient for nutrient in nutrientSet 
                                      if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ nutrient ] ]

newOrgSecPathDict = {}
for coreTBP in Core:
    currPath = newOrgPathDict[ coreTBP ]
    newOrgSecPathDict[ coreTBP ] = [ sec for sec in tempOrgSecs 
                                      if ( np.logical_or.reduce(rxnMat[np.array( currPath ).astype(int)]) * 1)[ sec ] ]

#-------------------------------------------------------------------------

while True:
    num_tries = 0
    good_inds, margins, removed_paths = [], [], []
    while True and num_tries < 1000:
        num_tries += 1
        try:
            rand_path = random.choice( secPathDict )
            new_cores = set(return_core(rand_path, Core))
            for n, i in enumerate(a):
                fluxes = dict(og_fluxes)
                if set(return_core(i, Core)) & new_cores:
                    new_ind = Individual([p for m, p in enumerate(a) if n != m] + [rand_path], stoich_matrix)
                    if cost_path(rand_path) > cost_path(i) and cost_path(i) > 0.0:
                        print('Success. New cost margin is ' + str(cost_path(rand_path) - cost_path(i)))
                        good_inds.append(new_ind)
                        margins.append(cost_path(rand_path) - cost_path(i))
                        removed_paths.append(i)
        except:
            pass

    break

# #-------------------------------------------------------------------------

# for new_ind in set(good_inds):
#     for other_ind in set(good_inds):
#         if other_ind != new_ind:
            new_ind = deepcopy(m)
            other_ind = deepcopy(n)
            N = 20
            fitness_lib = {}
            f, f1 = {}, {}
            flux = {}
            for q in range(N + 1):
                fluxes = {i: 1.0 for i in Env}
                init_pop = curr_pop = [other_ind] * q + [new_ind] * (N - q)
                Secreted_dict = return_secretion_dict(curr_pop)
                fitness_lib = update_fitness_lib(curr_pop, Env, ALPHA)
                if new_ind.__hash__() in list(fitness_lib.keys()):
                    f[q] = fitness_lib[new_ind.__hash__()]
                    if q > 0:
                        f1[q] = fitness_lib[other_ind.__hash__()]
                    else:
                        f1[q] = calc_fitness(other_ind, ALPHA, Env, curr_pop, fluxes, [])
                else:
                    f[q] = calc_fitness(new_ind, ALPHA, Env, curr_pop, fluxes, [])
                    f1[q] = calc_fitness(other_ind, ALPHA, Env, curr_pop, fluxes, [])

                print_progress_bar(q, N, 'Plotting')

            fig, ax = plt.subplots(1)
            ax.plot(list(range(N + 1)), list(f.values()))
            ax.plot(list(range(N + 1)), list(f1.values()))
            plt.show()

# #-------------------------------------------------------------------------

# NUM_SAMPLES = 100
# fixation_times = []
# for sample_num in range(NUM_SAMPLES):
#     print_progress_bar(sample_num + 1, NUM_SAMPLES, 'Calculating fixation times')
#     N = 100
#     q = 50
#     fluxes = {i: 1.0 for i in Env}
#     init_pop = curr_pop = [other_ind] * q + [new_ind] * (N - q)
#     Secreted_dict = return_secretion_dict(curr_pop)
#     fitness_lib = update_fitness_lib(curr_pop, Env, alpha)
#     NUM_ITERS = 1000
#     pop_comps = [tuple(return_pop_composition(curr_pop, fitness_lib).values())]
#     for iterNum in range(NUM_ITERS):
#         curr_pop = adv_gen(curr_pop, Env, 0.0)
#         pop_comps.append(tuple(return_pop_composition(curr_pop, fitness_lib).values()))
#         if 0 in pop_comps[-1]:
#             # print('Reached fixation')
#             break
#     fixation_times.append(len(pop_comps))

# fig, ax = plt.subplots(1)
# X_plot = np.linspace(0, 500, 1000)[:, np.newaxis]
# bwdith = 10

# kde1 = KernelDensity(kernel='gaussian', bandwidth=bwdith).fit(D)
# log_dens1 = kde1.score_samples(X_plot)
# ax.plot(X_plot[:, 0], np.exp(log_dens1), color='green', lw=1.3, label='KEGG')
# plt.show()

# fig, ax = plt.subplots(1)
# ax.plot(list(range(len(pop_comps))), [x[0] for x in pop_comps])
# ax.plot(list(range(len(pop_comps))), [x[1] for x in pop_comps])
# plt.show()

#-------------------------------------------------------------------------

