# Parthenogenesis Metapopulation Model v9

# changes:
# allele values back to continuous [0, 0.5) (from discrete)
# combined encounter + mating phase
# all female parthenogenetic offspring 
# initialized all arrays with -1 for safety
# plotting sexual + parthenogenetic offspring
# plotting sex ratio

# performance:
# using numba for jit compilation (speed)
# removed meeting matrix (memory)
# removed originalindex array (memory)
# slicing arrays instead of looping (speed)
# using numpy arrays instead of lists (speed)

import multiprocessing
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
from numba import jit
import os
import csv
import sys
inputs = sys.argv[1:]

# parameters
popsize = int(inputs[0])
maxsubsize = int(inputs[1])
migrecip = int(inputs[2])
migration = 0
if (migrecip != 0): migration = 1/migrecip
mutrecip = int(inputs[3])
mutation = 1/mutrecip
rec1 = float(inputs[4])
rec2 = float(inputs[5])
numinds = int(inputs[6])
numindssub = int(inputs[7])
maxrepro = int(inputs[8])
overdom = int(inputs[9]) # 0 = fitness dominance, 1 = fitness overdominance
numgens = int(inputs[10])
parthreduction = float(inputs[11])
parthrepro = int(parthreduction*maxrepro + 0.5)
parthpenalty = float(inputs[12])
plot_path = f"/gpfs/scratch/blukacsy/graphs/sim{inputs[13]}.png"

@jit(nopython=True)
def countMaleFemale(sex_list):
    counter_males = 0
    counter_females = 0
    for i in range(1, len(sex_list)):
        if (sex_list[i] == 0): counter_males += 1
        if (sex_list[i] == 1): counter_females += 1
    return counter_males, counter_females

@jit(nopython=True)
def run_phases(current_popsize, sex, location, loc1allele1, loc1allele2, loc2allele1, loc2allele2, loc3allele1, loc3allele2, offspringalive, offspringsex, offspringlocation, offspringloc1allele1, offspringloc1allele2, offspringloc2allele1, offspringloc2allele2, offspringloc3allele1, offspringloc3allele2, score, newoffspringalive, newoffspringsex, newoffspringlocation, newoffspringloc1allele1, newoffspringloc1allele2, newoffspringloc2allele1, newoffspringloc2allele2, newoffspringloc3allele1, newoffspringloc3allele2, mutation, rec1, rec2, numinds_for_popsize, maxrepro, overdom, parthrepro, parthpenalty):

    offspringalive[:] = -1
    offspringsex[:] = -1
    offspringlocation[:] = -1
    offspringloc1allele1[:] = -1
    offspringloc1allele2[:] = -1
    offspringloc2allele1[:] = -1.0
    offspringloc2allele2[:] = -1.0
    offspringloc3allele1[:] = -1.0
    offspringloc3allele2[:] = -1.0
    score[:] = -1.0
    newoffspringalive[:] = -1
    newoffspringsex[:] = -1
    newoffspringlocation[:] = -1
    newoffspringloc1allele1[:] = -1
    newoffspringloc1allele2[:] = -1
    newoffspringloc2allele1[:] = -1.0
    newoffspringloc2allele2[:] = -1.0
    newoffspringloc3allele1[:] = -1.0
    newoffspringloc3allele2[:] = -1.0

    # mutation phase
    for i in range(1, current_popsize+1):

        # locus 1 allele 1
        z = np.random.random()
        if (z < mutation): loc1allele1[i] = abs(loc1allele1[i] - 1)

        # locus 1 allele 2
        z = np.random.random()
        if (z < mutation): loc1allele2[i] = abs(loc1allele2[i] - 1)

        # locus 2 allele 1
        z = np.random.random()
        if (z < mutation): loc2allele1[i] = np.random.random()/2
                
        # locus 2 allele 2
        z = np.random.random()
        if (z < mutation): loc2allele2[i] = np.random.random()/2
            
        # locus 3 allele 1
        z = np.random.random()
        if (z < mutation): loc3allele1[i] = np.random.random()/2

        # locus 3 allele 2
        z = np.random.random()
        if (z < mutation): loc3allele2[i] = np.random.random()/2

    sexmothers = np.full(current_popsize+1, -1, np.int32)
    sexfathers = np.full(current_popsize+1, -1, np.int32)
    parthmothers = np.full(current_popsize+1, -1, np.int32)

    sexualmatingscount = 1
    parthmatingscount = 0

    # encounter + mating phase
    for a in range(1, current_popsize+1):
        matingflag = 0
        if (sex[a] == 1):
            x = numinds_for_popsize

            for b in range(1, x+1):

                zz = 0
                temp_list = []

                while (True):
                    if (len(temp_list) >= current_popsize): break
                    zz = np.random.randint(1, current_popsize+1)
                    if (a != zz): break
                    if (zz not in temp_list): temp_list.append(zz)

                if (len(temp_list) < current_popsize):
                    if (sex[zz] == 0):
                        matingflag = 1
                        sexmothers[sexualmatingscount] = a
                        sexfathers[sexualmatingscount] = zz
                        break # female will mate with first male she encounters

            if (matingflag == 1): 
                sexualmatingscount += 1

            else:
                parthmatingscount += 1
                parthmothers[parthmatingscount] = a

    # reproduction phase
    offspringcount = 0
    max_offspring = 0
    mhaplotype = 0
    phaplotype = 0

    # maternal alleles
    m_alleles = [
        (loc1allele1, loc2allele1, loc3allele1),
        (loc1allele1, loc2allele1, loc3allele2),
        (loc1allele1, loc2allele2, loc3allele2),
        (loc1allele1, loc2allele2, loc3allele1),
        (loc1allele2, loc2allele2, loc3allele2),
        (loc1allele2, loc2allele2, loc3allele1),
        (loc1allele2, loc2allele1, loc3allele1),
        (loc1allele2, loc2allele1, loc3allele2)
    ]

    # sexual reproduction
    for a in range(1, sexualmatingscount):
        mother = sexmothers[a]
        father = sexfathers[a]
        
        if (loc2allele1[mother] + loc2allele2[mother] == 0):
            max_offspring = maxrepro
        else:
            max_offspring = int(maxrepro * (1 - parthpenalty) + 0.5)

        for b in range(1, max_offspring+1):
            offspringcount += 1

            # determine maternal haplotype
            x = np.random.randint(0,2)
            y = np.random.random()
            z = np.random.random()

            # assign mhaplotype based on recombination
            mhaplotype = (x << 2) | ((y < rec1) << 1) | (z < rec2)

            alleles = m_alleles[mhaplotype]
            offspringloc1allele1[offspringcount] = alleles[0][mother]
            offspringloc2allele1[offspringcount] = alleles[1][mother]
            offspringloc3allele1[offspringcount] = alleles[2][mother]

            # determine paternal haplotype
            x = np.random.randint(0,2)
            y = np.random.random()
            z = np.random.random()

            # assign phaplotype based on recombination
            phaplotype = (x << 2) | ((y < rec1) << 1) | (z < rec2)

            # paternal alleles
            alleles = m_alleles[phaplotype]
            offspringloc1allele2[offspringcount] = alleles[0][father]
            offspringloc2allele2[offspringcount] = alleles[1][father]
            offspringloc3allele2[offspringcount] = alleles[2][father]

            offspringalive[offspringcount] = 1
            offspringsex[offspringcount] = 1 if np.random.random() < 0.5 else 0
            offspringlocation[offspringcount] = location[mother]

    sexual_offspring = offspringcount

    # parthenogenetic reproduction
    for a in range(1, parthmatingscount+1):
        mother = parthmothers[a]
        if (np.random.random() < (loc2allele1[mother] + loc2allele2[mother])):
            for b in range(1, parthrepro+1):
                offspringcount += 1

                # determine maternal haplotype
                x = np.random.randint(0,2)
                y = np.random.random()
                z = np.random.random()

                # assign mhaplotype based on recombination
                mhaplotype = (x << 2) | ((y < rec1) << 1) | (z < rec2)

                # parthenogenetic offspring represent duplicated maternal meiotic haplotypes
                alleles = m_alleles[mhaplotype]
                offspringloc1allele1[offspringcount] = alleles[0][mother]
                offspringloc1allele2[offspringcount] = alleles[0][mother]
                offspringloc2allele1[offspringcount] = alleles[1][mother]
                offspringloc2allele2[offspringcount] = alleles[1][mother]
                offspringloc3allele1[offspringcount] = alleles[2][mother]
                offspringloc3allele2[offspringcount] = alleles[2][mother]

                offspringalive[offspringcount] = 1
                offspringsex[offspringcount] = 1 # all female parthenogenetic offspring 
                offspringlocation[offspringcount] = location[mother]

    parth_offspring = offspringcount - sexual_offspring

    for a in range(1, offspringcount+1):
        score[a] = np.random.random()

    # randomly sort offspring pool based on score
    sorted_indices = np.argsort(score[1:offspringcount+1]) + 1

    for i in range(len(sorted_indices)):
        idx = i + 1
        a = sorted_indices[i]

        newoffspringalive[idx] = offspringalive[a]
        newoffspringsex[idx] = offspringsex[a]
        newoffspringlocation[idx] = offspringlocation[a]
        newoffspringloc1allele1[idx] = offspringloc1allele1[a]
        newoffspringloc1allele2[idx] = offspringloc1allele2[a]
        newoffspringloc2allele1[idx] = offspringloc2allele1[a]
        newoffspringloc2allele2[idx] = offspringloc2allele2[a]
        newoffspringloc3allele1[idx] = offspringloc3allele1[a]
        newoffspringloc3allele2[idx] = offspringloc3allele2[a]

    # selection phase
    for a in range(1, offspringcount+1):
        sum_alleles = newoffspringloc1allele1[a] + newoffspringloc1allele2[a]

        # fitness dominance
        if (overdom == 0):
            if (sum_alleles == 0): newoffspringalive[a] = 0

        # fitness overdominance
        elif (overdom == 1):
            if ((sum_alleles == 0) or (sum_alleles == 2)): newoffspringalive[a] = 0

    return offspringcount, sexual_offspring, parth_offspring

@jit(nopython=True)
def run_simulation(run):

    randomseed = 5002 + run
    np.random.seed(randomseed)

    # main population primary arrays
    sex_main = np.full(popsize+1, -1, dtype=np.int8)
    location_main = np.full(popsize+1, -1, dtype=np.int8)
    loc1allele1_main = np.full(popsize+1, -1, dtype=np.int8)
    loc1allele2_main = np.full(popsize+1, -1, dtype=np.int8)
    loc2allele1_main = np.full(popsize+1, -1.0, dtype=np.float32)
    loc2allele2_main = np.full(popsize+1, -1.0, dtype=np.float32)
    loc3allele1_main = np.full(popsize+1, -1.0, dtype=np.float32)
    loc3allele2_main = np.full(popsize+1, -1.0, dtype=np.float32)

    # subpopulation primary arrays
    sex_sub = np.full(maxsubsize+1, -1, dtype=np.int8)
    location_sub = np.full(maxsubsize+1, -1, dtype=np.int8)
    loc1allele1_sub = np.full(maxsubsize+1, -1, dtype=np.int8)
    loc1allele2_sub = np.full(maxsubsize+1, -1, dtype=np.int8)
    loc2allele1_sub = np.full(maxsubsize+1, -1.0, dtype=np.float32)
    loc2allele2_sub = np.full(maxsubsize+1, -1.0, dtype=np.float32)
    loc3allele1_sub = np.full(maxsubsize+1, -1.0, dtype=np.float32)
    loc3allele2_sub = np.full(maxsubsize+1, -1.0, dtype=np.float32)

    max_offspring_main = popsize*maxrepro
    max_offspring_sub = maxsubsize*maxrepro

    # main population offspring arrays
    offspringalive_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    offspringsex_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    offspringlocation_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    offspringloc1allele1_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    offspringloc1allele2_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    offspringloc2allele1_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)
    offspringloc2allele2_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)
    offspringloc3allele1_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)
    offspringloc3allele2_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)

    # subpopulation offspring arrays
    offspringalive_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    offspringsex_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    offspringlocation_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    offspringloc1allele1_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    offspringloc1allele2_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    offspringloc2allele1_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)
    offspringloc2allele2_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)
    offspringloc3allele1_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)
    offspringloc3allele2_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)

    # main population sorting + newoffspring arrays
    score_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)
    newoffspringalive_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    newoffspringsex_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    newoffspringlocation_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    newoffspringloc1allele1_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    newoffspringloc1allele2_main = np.full(max_offspring_main+1, -1, dtype=np.int8)
    newoffspringloc2allele1_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)
    newoffspringloc2allele2_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)
    newoffspringloc3allele1_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)
    newoffspringloc3allele2_main = np.full(max_offspring_main+1, -1.0, dtype=np.float32)

    # subpopulation sorting + newoffspring arrays
    score_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)
    newoffspringalive_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    newoffspringsex_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    newoffspringlocation_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    newoffspringloc1allele1_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    newoffspringloc1allele2_sub = np.full(max_offspring_sub+1, -1, dtype=np.int8)
    newoffspringloc2allele1_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)
    newoffspringloc2allele2_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)
    newoffspringloc3allele1_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)
    newoffspringloc3allele2_sub = np.full(max_offspring_sub+1, -1.0, dtype=np.float32)

    # tracking results arrays
    loc2freq_main = np.full(numgens+1, -1.0, dtype=np.float32)
    loc2freq_sub = np.full(numgens+1, -1.0, dtype=np.float32)
    loc3freq_main = np.full(numgens+1, -1.0, dtype=np.float32)
    loc3freq_sub = np.full(numgens+1, -1.0, dtype=np.float32)
    num_in_main = np.full(numgens+1, -1, dtype=np.int32)
    num_in_sub = np.full(numgens+1, -1, dtype=np.int32)
    num_sexual_offspring = np.full(numgens+1, -1, dtype=np.int32)
    num_parth_offspring = np.full(numgens+1, -1, dtype=np.int32)
    num_males = np.full(numgens+1, -1, dtype=np.int32)
    num_females = np.full(numgens+1, -1, dtype=np.int32)
    female_male_ratio = np.full(numgens+1, -1.0, dtype=np.float32)

    num_in_main[0] = popsize
    num_in_sub[0] = 0

    # initialize main population
    for i in range(1, popsize+1):
        aa = np.random.random()
        if (aa < 0.5): sex_main[i] = 1
        else: sex_main[i] = 0
        location_main[i] = 0
        loc1allele1_main[i] = 1
        if (overdom == 1): loc1allele2_main[i] = 0
        else: loc1allele2_main[i] = 1
        loc2allele1_main[i] = 0
        loc2allele2_main[i] = 0
        loc3allele1_main[i] = 0
        loc3allele2_main[i] = 0

    male_female_count = countMaleFemale(sex_main)
    num_males[0] = male_female_count[0]
    num_females[0] = male_female_count[1]
    female_male_ratio[0] = num_females[0] / max(1, num_males[0])

    # start main loop
    for gen in range(1, numgens+1):

        offspringcount_main = 0
        offspringcount_sub = 0

        sexual_offspringcount_main = 0
        sexual_offspringcount_sub = 0
        parth_offspringcount_main = 0
        parth_offspringcount_sub = 0
        
        current_num_in_main = num_in_main[gen-1]
        current_num_in_sub = num_in_sub[gen-1]

        # run main population
        if (current_num_in_main > 0):
            results = run_phases(current_num_in_main, sex_main, location_main, loc1allele1_main, loc1allele2_main, loc2allele1_main, loc2allele2_main, loc3allele1_main, loc3allele2_main, offspringalive_main, offspringsex_main, offspringlocation_main, offspringloc1allele1_main, offspringloc1allele2_main, offspringloc2allele1_main, offspringloc2allele2_main, offspringloc3allele1_main, offspringloc3allele2_main, score_main, newoffspringalive_main, newoffspringsex_main, newoffspringlocation_main, newoffspringloc1allele1_main, newoffspringloc1allele2_main, newoffspringloc2allele1_main, newoffspringloc2allele2_main, newoffspringloc3allele1_main, newoffspringloc3allele2_main, mutation, rec1, rec2, numinds, maxrepro, overdom, parthrepro, parthpenalty)
            offspringcount_main = results[0]
            sexual_offspringcount_main = results[1]
            parth_offspringcount_main = results[2]

        # run subpopulation
        if (current_num_in_sub > 0):
            results = run_phases(current_num_in_sub, sex_sub, location_sub, loc1allele1_sub, loc1allele2_sub, loc2allele1_sub, loc2allele2_sub, loc3allele1_sub, loc3allele2_sub, offspringalive_sub, offspringsex_sub, offspringlocation_sub, offspringloc1allele1_sub, offspringloc1allele2_sub, offspringloc2allele1_sub, offspringloc2allele2_sub, offspringloc3allele1_sub, offspringloc3allele2_sub, score_sub, newoffspringalive_sub, newoffspringsex_sub, newoffspringlocation_sub, newoffspringloc1allele1_sub, newoffspringloc1allele2_sub, newoffspringloc2allele1_sub, newoffspringloc2allele2_sub, newoffspringloc3allele1_sub, newoffspringloc3allele2_sub, mutation, rec1, rec2, numindssub, maxrepro, overdom, parthrepro, parthpenalty)
            offspringcount_sub = results[0]
            sexual_offspringcount_sub = results[1]
            parth_offspringcount_sub = results[2]

        # select main to sub migrants
        go_sub = []
        stay_main = []
        offspringcounter = 0

        for i in range(1, popsize+1):
            offspringcounter += 1
            while (offspringcounter <= offspringcount_main and newoffspringalive_main[offspringcounter] == 0):
                offspringcounter += 1

            if (offspringcounter > offspringcount_main):
                break

            # migration chance
            if (np.random.random() < migration): go_sub.append(offspringcounter)
            else: stay_main.append(offspringcounter)

        # select sub to main migrants
        go_main = []
        stay_sub = []
        offspringcounter = 0

        for i in range(1, maxsubsize+1):
            offspringcounter += 1
            while (offspringcounter <= offspringcount_sub and newoffspringalive_sub[offspringcounter] == 0):
                offspringcounter += 1

            if (offspringcounter > offspringcount_sub):
                break

            # migration chance
            if (np.random.random() < migration): go_main.append(offspringcounter)
            else: stay_sub.append(offspringcounter)

        # rebuild main
        counter = 1
        for i in stay_main:
            sex_main[counter] = newoffspringsex_main[i]
            location_main[counter] = 0
            loc1allele1_main[counter] = newoffspringloc1allele1_main[i]
            loc1allele2_main[counter] = newoffspringloc1allele2_main[i]
            loc2allele1_main[counter] = newoffspringloc2allele1_main[i]
            loc2allele2_main[counter] = newoffspringloc2allele2_main[i]
            loc3allele1_main[counter] = newoffspringloc3allele1_main[i]
            loc3allele2_main[counter] = newoffspringloc3allele2_main[i]

            counter += 1

        # rebuild sub
        counter = 1
        for i in stay_sub:
            sex_sub[counter] = newoffspringsex_sub[i]
            location_sub[counter] = 1
            loc1allele1_sub[counter] = newoffspringloc1allele1_sub[i]
            loc1allele2_sub[counter] = newoffspringloc1allele2_sub[i]
            loc2allele1_sub[counter] = newoffspringloc2allele1_sub[i]
            loc2allele2_sub[counter] = newoffspringloc2allele2_sub[i]
            loc3allele1_sub[counter] = newoffspringloc3allele1_sub[i]
            loc3allele2_sub[counter] = newoffspringloc3allele2_sub[i]

            counter += 1

        # append main migrants to sub
        current_num_in_sub = len(stay_sub)
        for i in go_sub:
            if (current_num_in_sub >= maxsubsize): break

            current_num_in_sub += 1

            sex_sub[current_num_in_sub] = newoffspringsex_main[i]
            location_sub[current_num_in_sub] = 1
            loc1allele1_sub[current_num_in_sub] = newoffspringloc1allele1_main[i]
            loc1allele2_sub[current_num_in_sub] = newoffspringloc1allele2_main[i]
            loc2allele1_sub[current_num_in_sub] = newoffspringloc2allele1_main[i]
            loc2allele2_sub[current_num_in_sub] = newoffspringloc2allele2_main[i]
            loc3allele1_sub[current_num_in_sub] = newoffspringloc3allele1_main[i]
            loc3allele2_sub[current_num_in_sub] = newoffspringloc3allele2_main[i]

        # append sub migrants to main
        current_num_in_main = len(stay_main)
        for i in go_main:
            if (current_num_in_main >= popsize): break

            current_num_in_main += 1

            sex_main[current_num_in_main] = newoffspringsex_sub[i]
            location_main[current_num_in_main] = 0
            loc1allele1_main[current_num_in_main] = newoffspringloc1allele1_sub[i]
            loc1allele2_main[current_num_in_main] = newoffspringloc1allele2_sub[i]
            loc2allele1_main[current_num_in_main] = newoffspringloc2allele1_sub[i]
            loc2allele2_main[current_num_in_main] = newoffspringloc2allele2_sub[i]
            loc3allele1_main[current_num_in_main] = newoffspringloc3allele1_sub[i]
            loc3allele2_main[current_num_in_main] = newoffspringloc3allele2_sub[i]

        current_num_in_main = min((len(stay_main) + len(go_main)), popsize)
        current_num_in_sub = min((len(stay_sub) + len(go_sub)), maxsubsize)

        # clear zombies
        sex_main[current_num_in_main+1:] = -1
        location_main[current_num_in_main+1:] = -1
        loc1allele1_main[current_num_in_main+1:] = -1
        loc1allele2_main[current_num_in_main+1:] = -1
        loc2allele1_main[current_num_in_main+1:] = -1.0
        loc2allele2_main[current_num_in_main+1:] = -1.0
        loc3allele1_main[current_num_in_main+1:] = -1.0
        loc3allele2_main[current_num_in_main+1:] = -1.0
        sex_sub[current_num_in_sub+1:] = -1
        location_sub[current_num_in_sub+1:] = -1
        loc1allele1_sub[current_num_in_sub+1:] = -1
        loc1allele2_sub[current_num_in_sub+1:] = -1
        loc2allele1_sub[current_num_in_sub+1:] = -1.0
        loc2allele2_sub[current_num_in_sub+1:] = -1.0
        loc3allele1_sub[current_num_in_sub+1:] = -1.0
        loc3allele2_sub[current_num_in_sub+1:] = -1.0

        # compute and store generation data
        loc2count_main = sum(loc2allele1_main[1:current_num_in_main+1]) + sum(loc2allele2_main[1:current_num_in_main+1])
        loc3count_main = sum(loc3allele1_main[1:current_num_in_main+1]) + sum(loc3allele2_main[1:current_num_in_main+1])

        if (current_num_in_main > 0):
            loc2freq_main[gen] = loc2count_main / (2 * current_num_in_main)
            loc3freq_main[gen] = loc3count_main / (2 * current_num_in_main)
        else:
            loc2freq_main[gen] = 0
            loc3freq_main[gen] = 0

        loc2count_sub = sum(loc2allele1_sub[1:current_num_in_sub+1]) + sum(loc2allele2_sub[1:current_num_in_sub+1])
        loc3count_sub = sum(loc3allele1_sub[1:current_num_in_sub+1]) + sum(loc3allele2_sub[1:current_num_in_sub+1])

        if (current_num_in_sub > 0):
            loc2freq_sub[gen] = loc2count_sub / (2 * current_num_in_sub)
            loc3freq_sub[gen] = loc3count_sub / (2 * current_num_in_sub)
        else:
            loc2freq_sub[gen] = 0
            loc3freq_sub[gen] = 0
        
        num_in_main[gen] = current_num_in_main
        num_in_sub[gen] = current_num_in_sub

        num_sexual_offspring[gen] = sexual_offspringcount_main + sexual_offspringcount_sub
        num_parth_offspring[gen] = parth_offspringcount_main + parth_offspringcount_sub

        male_female_count_main = countMaleFemale(sex_main)
        male_female_count_sub = countMaleFemale(sex_sub)
        num_males[gen] = male_female_count_main[0] + male_female_count_sub[0]
        num_females[gen] = male_female_count_main[1] + male_female_count_sub[1]
        female_male_ratio[gen] = num_females[gen] / max(1, num_males[gen])

    y_axis_loc2_main = loc2freq_main[0:numgens+1]
    y_axis_loc2_sub = loc2freq_sub[0:numgens+1]
    y_axis_loc3_main = loc3freq_main[0:numgens+1]
    y_axis_loc3_sub = loc3freq_sub[0:numgens+1]
    y_axis_num_main = num_in_main[0:numgens+1]
    y_axis_num_sub = num_in_sub[0:numgens+1]
    y_axis_num_sexual_offspring = num_sexual_offspring[0:numgens+1]
    y_axis_num_parth_offspring = num_parth_offspring[0:numgens+1]
    y_axis_female_male_ratio = female_male_ratio[0:numgens+1]

    return y_axis_loc2_main, y_axis_loc2_sub, y_axis_loc3_main, y_axis_loc3_sub, y_axis_num_main, y_axis_num_sub, y_axis_num_sexual_offspring, y_axis_num_parth_offspring, y_axis_female_male_ratio

if __name__ == '__main__':
    total_runs = 12
    num_processes = 12

    total_loc2_allele_freq_main = []
    total_loc2_allele_freq_sub = []
    total_loc3_allele_freq_main = []
    total_loc3_allele_freq_sub = []
    total_num_main = []
    total_num_sub = []
    total_num_sexual_offspring = []
    total_num_parth_offspring = []
    total_female_male_ratio = []

    with multiprocessing.Pool(processes=num_processes) as pool:
        for result in tqdm(pool.imap_unordered(run_simulation, range(total_runs)), total=total_runs):
            total_loc2_allele_freq_main.append(result[0])
            total_loc2_allele_freq_sub.append(result[1])
            total_loc3_allele_freq_main.append(result[2])
            total_loc3_allele_freq_sub.append(result[3])
            total_num_main.append(result[4])
            total_num_sub.append(result[5])
            total_num_sexual_offspring.append(result[6])
            total_num_parth_offspring.append(result[7])
            total_female_male_ratio.append(result[8])

    param_text = (f"runs={total_runs}, numgens={numgens}, popsize={popsize}, maxsubsize={maxsubsize}, migration={migration}, mutation={mutation}, numinds={numinds}, numindssub={numindssub}, maxrepro={maxrepro}, overdom={overdom}, parthreduction={parthreduction}, parthpenalty={parthpenalty}")

    # plotting
    figure, axis = plt.subplots(3, 3, figsize=(23, 20))

    # plot 1
    axis[0][0].set_xlim(0, numgens)
    axis[0][0].set_ylim(0, 0.55)
    axis[0][0].set_xlabel("Generations")
    axis[0][0].set_ylabel("Frequency")
    axis[0][0].set_title("Average Locus 2 (Parthenogenetic) Allele Frequency Main Population")
    for run in range(total_runs):
        axis[0][0].plot(total_loc2_allele_freq_main[run], color=plt.cm.rainbow(run / total_runs))

    # plot 2
    axis[1][0].set_xlim(0, numgens)
    axis[1][0].set_ylim(0, 0.55)
    axis[1][0].set_xlabel("Generations")
    axis[1][0].set_ylabel("Frequency")
    axis[1][0].set_title("Average Locus 2 (Parthenogenetic) Allele Frequency Subpopulation")
    for run in range(total_runs):
        axis[1][0].plot(total_loc2_allele_freq_sub[run], color=plt.cm.rainbow(run / total_runs))

    # plot 3
    axis[0][1].set_xlim(0, numgens)
    axis[0][1].set_ylim(0, 0.55)
    axis[0][1].set_xlabel("Generations")
    axis[0][1].set_ylabel("Frequency")
    axis[0][1].set_title("Average Locus 3 (Neutral) Allele Frequency Main Population")
    for run in range(total_runs):
        axis[0][1].plot(total_loc3_allele_freq_main[run], color=plt.cm.rainbow(run / total_runs))

    # plot 4
    axis[1][1].set_xlim(0, numgens)
    axis[1][1].set_ylim(0, 0.55)
    axis[1][1].set_xlabel("Generations")
    axis[1][1].set_ylabel("Frequency")
    axis[1][1].set_title("Average Locus 3 (Neutral) Allele Frequency Subpopulation")
    for run in range(total_runs):
        axis[1][1].plot(total_loc3_allele_freq_sub[run], color=plt.cm.rainbow(run / total_runs))

    # plot 5
    axis[0][2].set_xlim(0, numgens)
    axis[0][2].set_ylim(0, popsize)
    axis[0][2].set_xlabel("Generations")
    axis[0][2].set_ylabel("Quantity")
    axis[0][2].set_title("Number of Individuals Main Population")
    for run in range(total_runs):
        axis[0][2].plot(total_num_main[run], color=plt.cm.rainbow(run / total_runs))

    # plot 6
    axis[1][2].set_xlim(0, numgens)
    axis[1][2].set_ylim(0, maxsubsize)
    axis[1][2].set_xlabel("Generations")
    axis[1][2].set_ylabel("Quantity")
    axis[1][2].set_title("Number of Individuals Subpopulation")
    for run in range(total_runs):
        axis[1][2].plot(total_num_sub[run], color=plt.cm.rainbow(run / total_runs))

    # plot 7
    axis[2][0].set_xlim(0, numgens)
    axis[2][0].set_ylim(0, popsize*maxrepro)
    axis[2][0].set_xlabel("Generations")
    axis[2][0].set_ylabel("Quantity")
    axis[2][0].set_title("Number of Sexual Offspring")
    for run in range(total_runs):
        axis[2][0].plot(total_num_sexual_offspring[run], color=plt.cm.rainbow(run / total_runs), marker = 'o', linestyle='none', mfc='none', markersize=5)

    # plot 8
    axis[2][1].set_xlim(0, numgens)
    axis[2][1].set_ylim(0, popsize*parthrepro)
    axis[2][1].set_xlabel("Generations")
    axis[2][1].set_ylabel("Quantity")
    axis[2][1].set_title("Number of Parthenogenetic Offspring")
    for run in range(total_runs):
        axis[2][1].plot(total_num_parth_offspring[run], color=plt.cm.rainbow(run / total_runs), marker = 'o', linestyle='none', mfc='none', markersize=5)

    # plot 9
    axis[2][2].set_xlim(0, numgens)
    axis[2][2].set_ylim(0, 3)
    axis[2][2].set_xlabel("Generations")
    axis[2][2].set_ylabel("Ratio")
    axis[2][2].set_title("Ratio of Females to Males")
    for run in range(total_runs):
        axis[2][2].plot(total_female_male_ratio[run], color=plt.cm.rainbow(run / total_runs), marker = 'o', linestyle='none', mfc='none', markersize=5)

    figure.suptitle(param_text, fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.97])
    plt.savefig(plot_path)

    start = int(0.9*numgens)
    def count_fixed_alleles(total_allele_freq):
        epsilon = 0.05
        count_fixed = 0
        for run in total_allele_freq:
            truncated_run = run[start:]
            saw_max = False
            attempts = 0
            for freq in truncated_run:
                if (freq >= (0.5 - epsilon)): saw_max = True
                if (freq < (0.5 - epsilon)):
                    if (attempts >= 3): 
                        saw_max = False
                        break
                    attempts += 1
            if (saw_max): count_fixed += 1
        return count_fixed

    count_locus2_main = count_fixed_alleles(total_loc2_allele_freq_main)
    count_locus2_sub = count_fixed_alleles(total_loc2_allele_freq_sub)
    count_locus3_main = count_fixed_alleles(total_loc3_allele_freq_main)
    count_locus3_sub = count_fixed_alleles(total_loc3_allele_freq_sub)

    output = "/gpfs/scratch/blukacsy/graphs/results.csv"
    file_not_exist = ((not os.path.exists(output)) or (os.path.getsize(output) == 0))
    with open(output, "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        if (file_not_exist): writer.writerow(["ID", "Locus_2_Main", "Locus_2_Sub", "Locus_3_Main", "Locus_3_Sub"])
        writer.writerow([inputs[13], count_locus2_main, count_locus2_sub, count_locus3_main, count_locus3_sub])