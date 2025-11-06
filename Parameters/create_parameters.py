def create_parameters(filename):
    
    popsize = 50000
    rec1 = 0.5
    rec2 = 0.5
    maxrepro = 300
    overdom = 0
    numgens = 10000
    parthreduction = 0.1

    # 5k: [int(popsize*0.1), int(popsize*0.02), int(popsize*0.005)]
    # 50k: [int(popsize*0.1), int(popsize*0.02), int(popsize*0.004), int(popsize*0.0005)]
    # 500k: [int(popsize*0.1), int(popsize*0.02), int(popsize*0.004), int(popsize*0.0008), int(popsize*0.00005)] 

    maxsubsize_combos = [int(popsize*0.1), int(popsize*0.02), int(popsize*0.004), int(popsize*0.0005)]
    migration_combos = [popsize*10]
    mutation_combos = [10**4, 10**5, 10**7, 10**8]
    encounter_combos = [(1, 1), (1, 2), (2, 1), (2, 2), (1, 10), (10, 1), (12, 6), (24, 12), (50, 25)]
    parthpenalty_combos = [0, 1/300, 2/300, 9/300, 12/300]

    header = ("popsize, maxsubsize, migration, mutation, rec1, rec2, numinds, numindssub, maxrepro, overdom, numgens, parthreduction, parthpenalty, plotnumber")
    parameters = [header, ""]
    total_combos = len(maxsubsize_combos) * len(migration_combos) * len(mutation_combos) * len(encounter_combos) * len(parthpenalty_combos)

    counter = 1
    def create_plotnumber(i):
        plotnumber = f"{i:0{len(str(total_combos))}d}"
        return plotnumber

    for maxsubsize in maxsubsize_combos:
        for migration in migration_combos:
            for mutation in mutation_combos:
                for encounter in encounter_combos:
                    for parthpenalty in parthpenalty_combos:

                        numinds, numindssub = encounter
                        plotnumber = create_plotnumber(counter)
                        counter += 1

                        parameter = [popsize, maxsubsize, migration, mutation, rec1, rec2, numinds, numindssub, maxrepro, overdom, numgens, parthreduction, parthpenalty, plotnumber]
                        
                        parameter_string = " ".join(map(str, parameter))
                        parameters.append(parameter_string)
                    
    with open(filename, 'w') as file:
        file.write("\n".join(parameters))

create_parameters("popsize_50k_gens_10k_parthreduction_0.1.txt")