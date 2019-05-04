a = des_BMG(Env)
search_times, gms = [], []
for ns in range(1):
    print_progress_bar(ns + 1, 1000, 'Finding search times to first obligate mutualism')
    cpDict = {deepcopy(a): 1.0}
    secGoods = uniqify(
    unlistify(return_secretion_dict(list(cpDict.keys())).values()))

    tests = 0
    while len(search_times) < ns + 1:
        if len(cpDict.keys()) > 1.0:
            chosen_parent = np.random.choice(
                list(cpDict.keys()), 1, p=list(cpDict.values()))[0]
        else:
            chosen_parent = deepcopy(list(cpDict.keys())[0])
        secGoods = uniqify(unlistify(return_secretion_dict(list(cpDict.keys())).values()))

        chosen_parent.mutate(secGoods, path_array)

        if calc_fitness(chosen_parent, ALPHA, Env, list(cpDict.keys()) + [chosen_parent], fluxes, []):
            cpDict[chosen_parent] = 0.01

        cpDict, updFitness = returnPopSS(cpDict)
        tests += 1
        # print(cpDict.values())
        if len(cpDict.values()) == 2:
            species = list(cpDict.keys())
            if np.array_equal(get_intMatrix(species, Env, Core), [[0,1],[1,0]]):
                if not (calc_fitness(species[0], ALPHA, Env, species[0], [], []) and
                        calc_fitness(species[1], ALPHA, Env, species[1], [], [])):
                    search_times.append(tests)
                    gms.append(chosen_parent)
                    break
