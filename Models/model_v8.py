# Parthenogenesis Metapopulation Model v8
# bug fix: was double counting migrants
# added explicit clearing of previous generation arrays (redundancy)
# reordered migration structure: fill main, fill sub, collect migrants, append migrants

import random
import multiprocessing
from tqdm import tqdm
import matplotlib.pyplot as plt
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
numindssub = min(int(0.04*maxsubsize), numinds//2)
maxrepro = int(inputs[7])
overdom = int(inputs[8]) # 0 = fitness dominance, 1 = fitness overdominance
numgens = int(inputs[9])
parthreduction = float(inputs[10])
parthrepro = int(parthreduction*maxrepro + 0.5)
parthpenalty = float(inputs[11])
plot_path = f"../graphs/sim{inputs[12]}.png"

def createMatrix(rows, cols):
    rows += 1
    cols += 1
    return [[0 for x in range(cols)] for y in range(rows)]

def clearArray(array):
    for i in range(len(array)):
        array[i] = 0

def clearZombies(array, index):
    for i in range(index+1, len(array)):
        array[i] = 0

def run_phases(current_popsize, sex, location, loc1allele1, loc1allele2, loc2allele1, loc2allele2, loc3allele1, loc3allele2, offspringalive, offspringsex, offspringlocation, offspringloc1allele1, offspringloc1allele2, offspringloc2allele1, offspringloc2allele2, offspringloc3allele1, offspringloc3allele2, score, originalindex, newoffspringalive, newoffspringsex, newoffspringlocation, newoffspringloc1allele1, newoffspringloc1allele2, newoffspringloc2allele1, newoffspringloc2allele2, newoffspringloc3allele1, newoffspringloc3allele2, meeting, mutation, rec1, rec2, numinds_for_popsize, maxrepro, overdom, parthrepro, parthpenalty):

    clearArray(offspringalive)
    clearArray(offspringsex)
    clearArray(offspringlocation)
    clearArray(offspringloc1allele1)
    clearArray(offspringloc1allele2)
    clearArray(offspringloc2allele1)
    clearArray(offspringloc2allele2)
    clearArray(offspringloc3allele1)
    clearArray(offspringloc3allele2)
    clearArray(score)
    clearArray(originalindex)
    clearArray(newoffspringalive)
    clearArray(newoffspringsex)
    clearArray(newoffspringlocation)
    clearArray(newoffspringloc1allele1)
    clearArray(newoffspringloc1allele2)
    clearArray(newoffspringloc2allele1)
    clearArray(newoffspringloc2allele2)
    clearArray(newoffspringloc3allele1)
    clearArray(newoffspringloc3allele2)

    # mutation phase
    for i in range(1, current_popsize+1):

        # locus 1 allele 1
        z = random.random()
        if z < mutation: loc1allele1[i] = abs(loc1allele1[i] - 1)

        # locus 1 allele 2
        z = random.random()
        if z < mutation: loc1allele2[i] = abs(loc1allele2[i] - 1)

        # locus 2 allele 1
        z = random.random()
        if z < mutation: loc2allele1[i] = 0.5
        
        # locus 2 allele 2
        z = random.random()
        if z < mutation: loc2allele2[i] = 0.5
        
        # locus 3 allele 1
        z = random.random()
        if z < mutation: loc3allele1[i] = 0.5

        # locus 3 allele 2
        z = random.random()
        if z < mutation: loc3allele2[i] = 0.5

    # encounter phase
    sexmothers = [0]*(current_popsize + 1)
    sexfathers = [0]*(current_popsize + 1)
    parthmothers = [0]*(current_popsize + 1)

    for a in range(1, current_popsize+1):
        x = numinds_for_popsize
        
        for b in range(1, x+1):
            temp_list = []

            while True:
                if len(temp_list) >= current_popsize: break
                zz = random.randint(1, current_popsize)
                if a != zz and location[a] == location[zz]: break
                if zz not in temp_list: temp_list.append(zz)

            if len(temp_list) < current_popsize:
                if b <= len(meeting[a]) - 1:
                    meeting[a][b] = zz

    sexualmatingscount = 1
    parthmatingscount = 0

    # find male that each female mates with
    for a in range(1, current_popsize+1):
        matingflag = 0
        if sex[a] == 1:
            x = numinds_for_popsize

            for b in range(1, x+1):
                if sex[meeting[a][b]] == 0:
                    matingflag = 1
                    sexmothers[sexualmatingscount] = a
                    sexfathers[sexualmatingscount] = meeting[a][b]
                    break # female will mate with first male she encounters instead of last

            if matingflag == 1: 
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
        
        if loc2allele1[mother] + loc2allele2[mother] == 0:
            max_offspring = maxrepro
        else:
            max_offspring = int(maxrepro * (1 - parthpenalty) + 0.5)

        for b in range(1, max_offspring+1):
            offspringcount += 1

            # determine maternal haplotype
            x = random.randint(0,1)
            y = random.random()
            z = random.random()

            # assign mhaplotype based on recombination
            mhaplotype = (x << 2) | ((y < rec1) << 1) | (z < rec2)

            alleles = m_alleles[mhaplotype]
            offspringloc1allele1[offspringcount] = alleles[0][mother]
            offspringloc2allele1[offspringcount] = alleles[1][mother]
            offspringloc3allele1[offspringcount] = alleles[2][mother]

            # determine paternal haplotype
            x = random.randint(0,1)
            y = random.random()
            z = random.random()

            # assign phaplotype based on recombination
            phaplotype = (x << 2) | ((y < rec1) << 1) | (z < rec2)

            # paternal alleles
            alleles = m_alleles[phaplotype]
            offspringloc1allele2[offspringcount] = alleles[0][father]
            offspringloc2allele2[offspringcount] = alleles[1][father]
            offspringloc3allele2[offspringcount] = alleles[2][father]

            # set other offspring parameters
            offspringalive[offspringcount] = 1
            offspringsex[offspringcount] = 1 if random.random() < 0.5 else 0
            offspringlocation[offspringcount] = location[mother]

    # record number of sexual offspring
    sexual_offspring = offspringcount

    # parthenogenetic reproduction
    for a in range(1, parthmatingscount+1):
        mother = parthmothers[a]
        if random.random() < (loc2allele1[mother] + loc2allele2[mother]):
            for b in range(1, parthrepro+1):
                offspringcount += 1

                # determine maternal haplotype
                x = random.randint(0,1)
                y = random.random()
                z = random.random()

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

                # set other offspring parameters
                offspringalive[offspringcount] = 1
                offspringsex[offspringcount] = 1 if random.random() < 0.5 else 0
                offspringlocation[offspringcount] = location[mother]

    # record number of parthenogenetic offspring
    parth_offspring = offspringcount - sexual_offspring

    # randomly sort offspring pool
    for a in range(1, offspringcount+1):
        score[a] = random.random()
        originalindex[a] = a

    # sort based on score
    combined = list(zip(score[1:offspringcount+1], originalindex[1:offspringcount+1]))
    combined.sort()

    sorted_indices = [idx for _, idx in combined]
    for idx, a in enumerate(sorted_indices, start=1):
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
        if overdom == 0:
            # fitness dominance
            if newoffspringloc1allele1[a] + newoffspringloc1allele2[a] == 0:
                newoffspringalive[a] = 0
        else:
            # fitness overdominance
            total_alleles = newoffspringloc1allele1[a] + newoffspringloc1allele2[a]
            if total_alleles == 0 or total_alleles == 2:
                newoffspringalive[a] = 0

    return offspringcount

def run_simulation(run):

    randomseed = random.random() + run  # ensure different seeds
    random.seed(randomseed)

    # main population primary arrays
    sex_main = [0]*(popsize+1)
    sex_main[0] = -1 # prevents unintended mating
    location_main = [0]*(popsize+1)
    loc1allele1_main = [0]*(popsize+1)
    loc1allele2_main = [0]*(popsize+1)
    loc2allele1_main = [0]*(popsize+1)
    loc2allele2_main = [0]*(popsize+1)
    loc3allele1_main = [0]*(popsize+1)
    loc3allele2_main = [0]*(popsize+1)

    # subpopulation primary arrays
    sex_sub = [0]*(maxsubsize+1)
    sex_sub[0] = -1 # prevents unintended mating
    location_sub = [1]*(maxsubsize+1)
    loc1allele1_sub = [0]*(maxsubsize+1)
    loc1allele2_sub = [0]*(maxsubsize+1)
    loc2allele1_sub = [0]*(maxsubsize+1)
    loc2allele2_sub = [0]*(maxsubsize+1)
    loc3allele1_sub = [0]*(maxsubsize+1)
    loc3allele2_sub = [0]*(maxsubsize+1)

    # main population offspring arrays
    offspringalive_main = [0]*(popsize*maxrepro + 1)
    offspringsex_main = [0]*(popsize*maxrepro + 1)
    offspringlocation_main = [0]*(popsize*maxrepro + 1)
    offspringloc1allele1_main = [0]*(popsize*maxrepro + 1)
    offspringloc1allele2_main = [0]*(popsize*maxrepro + 1)
    offspringloc2allele1_main = [0]*(popsize*maxrepro + 1)
    offspringloc2allele2_main = [0]*(popsize*maxrepro + 1)
    offspringloc3allele1_main = [0]*(popsize*maxrepro + 1)
    offspringloc3allele2_main = [0]*(popsize*maxrepro + 1)

    # subpopulation offspring arrays
    offspringalive_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringsex_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringlocation_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringloc1allele1_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringloc1allele2_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringloc2allele1_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringloc2allele2_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringloc3allele1_sub = [0]*(maxsubsize*maxrepro + 1)
    offspringloc3allele2_sub = [0]*(maxsubsize*maxrepro + 1)

    # main population sorting + newoffspring arrays
    score_main = [0]*(popsize*maxrepro + 1)
    originalindex_main = [0]*(popsize*maxrepro + 1)
    newoffspringalive_main = [0]*(popsize*maxrepro + 1)
    newoffspringsex_main = [0]*(popsize*maxrepro + 1)
    newoffspringlocation_main = [0]*(popsize*maxrepro + 1)
    newoffspringloc1allele1_main = [0]*(popsize*maxrepro + 1)
    newoffspringloc1allele2_main = [0]*(popsize*maxrepro + 1)
    newoffspringloc2allele1_main = [0]*(popsize*maxrepro + 1)
    newoffspringloc2allele2_main = [0]*(popsize*maxrepro + 1)
    newoffspringloc3allele1_main = [0]*(popsize*maxrepro + 1)
    newoffspringloc3allele2_main = [0]*(popsize*maxrepro + 1)

    # subpopulation sorting + newoffspring arrays
    score_sub = [0]*(maxsubsize*maxrepro + 1)
    originalindex_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringalive_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringsex_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringlocation_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringloc1allele1_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringloc1allele2_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringloc2allele1_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringloc2allele2_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringloc3allele1_sub = [0]*(maxsubsize*maxrepro + 1)
    newoffspringloc3allele2_sub = [0]*(maxsubsize*maxrepro + 1)

    # meeting arrays
    meeting_main = createMatrix(popsize, numinds)
    meeting_sub = createMatrix(maxsubsize, numindssub)

    # arrays for tracking main and subpopulation frequencies separately
    loc2freq_main = [0] * (2*numgens + 1)
    loc2freq_sub = [0] * (2*numgens + 1)
    loc3freq_main = [0] * (2*numgens + 1)
    loc3freq_sub = [0] * (2*numgens + 1)

    # arrays for tracking number of individuals in main and subpopulation
    num_in_main = [0] * (2*numgens + 1)
    num_in_sub = [0] * (2*numgens + 1)
    num_in_main[0] = popsize
    num_in_sub[0] = 0

    # initialize main population
    for i in range(1, popsize+1):
        aa = random.random()
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

    # start main loop
    for gen in range(1, numgens+1):

        offspringcount_main = 0
        offspringcount_sub = 0
        
        current_num_in_main = num_in_main[gen-1]
        current_num_in_sub = num_in_sub[gen-1]

        # run main population
        if current_num_in_main > 0:
            offspringcount_main = run_phases(current_num_in_main, sex_main, location_main, loc1allele1_main, loc1allele2_main, loc2allele1_main, loc2allele2_main, loc3allele1_main, loc3allele2_main, offspringalive_main, offspringsex_main, offspringlocation_main, offspringloc1allele1_main, offspringloc1allele2_main, offspringloc2allele1_main, offspringloc2allele2_main, offspringloc3allele1_main, offspringloc3allele2_main, score_main, originalindex_main, newoffspringalive_main, newoffspringsex_main, newoffspringlocation_main, newoffspringloc1allele1_main, newoffspringloc1allele2_main, newoffspringloc2allele1_main, newoffspringloc2allele2_main, newoffspringloc3allele1_main, newoffspringloc3allele2_main, meeting_main, mutation, rec1, rec2, numinds, maxrepro, overdom, parthrepro, parthpenalty)

        # run subpopulation
        if current_num_in_sub > 0:
            offspringcount_sub = run_phases(current_num_in_sub, sex_sub, location_sub, loc1allele1_sub, loc1allele2_sub, loc2allele1_sub, loc2allele2_sub, loc3allele1_sub, loc3allele2_sub, offspringalive_sub, offspringsex_sub, offspringlocation_sub, offspringloc1allele1_sub, offspringloc1allele2_sub, offspringloc2allele1_sub, offspringloc2allele2_sub, offspringloc3allele1_sub, offspringloc3allele2_sub, score_sub, originalindex_sub, newoffspringalive_sub, newoffspringsex_sub, newoffspringlocation_sub, newoffspringloc1allele1_sub, newoffspringloc1allele2_sub, newoffspringloc2allele1_sub, newoffspringloc2allele2_sub, newoffspringloc3allele1_sub, newoffspringloc3allele2_sub, meeting_sub, mutation, rec1, rec2, numindssub, maxrepro, overdom, parthrepro, parthpenalty)

        # select main to sub migrants
        go_sub = []
        stay_main = []
        offspringcounter = 0

        for i in range(1, popsize+1):
            offspringcounter += 1
            while offspringcounter <= offspringcount_main and newoffspringalive_main[offspringcounter] == 0:
                offspringcounter += 1

            if offspringcounter > offspringcount_main:
                break

            # migration chance
            if random.random() < migration: go_sub.append(offspringcounter)
            else: stay_main.append(offspringcounter)

        # select sub to main migrants
        go_main = []
        stay_sub = []
        offspringcounter = 0

        for i in range(1, maxsubsize+1):
            offspringcounter += 1
            while offspringcounter <= offspringcount_sub and newoffspringalive_sub[offspringcounter] == 0:
                offspringcounter += 1

            if offspringcounter > offspringcount_sub:
                break

            # migration chance
            if random.random() < migration: go_main.append(offspringcounter)
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
            if current_num_in_sub >= maxsubsize: break

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
            if current_num_in_main >= popsize: break

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

        clearZombies(sex_main, current_num_in_main)
        clearZombies(loc1allele1_main, current_num_in_main)
        clearZombies(loc1allele2_main, current_num_in_main)
        clearZombies(loc2allele1_main, current_num_in_main)
        clearZombies(loc2allele2_main, current_num_in_main)
        clearZombies(loc3allele1_main, current_num_in_main)
        clearZombies(loc3allele2_main, current_num_in_main)

        clearZombies(sex_sub, current_num_in_sub)
        clearZombies(loc1allele1_sub, current_num_in_sub)
        clearZombies(loc1allele2_sub, current_num_in_sub)
        clearZombies(loc2allele1_sub, current_num_in_sub)
        clearZombies(loc2allele2_sub, current_num_in_sub)
        clearZombies(loc3allele1_sub, current_num_in_sub)
        clearZombies(loc3allele2_sub, current_num_in_sub)

        sex_main[0] = -1
        sex_sub[0] = -1

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

    # prepare data for plotting
    y_axis_loc2_main = loc2freq_main[1:numgens+1]
    y_axis_loc2_sub = loc2freq_sub[1:numgens+1]
    y_axis_loc3_main = loc3freq_main[1:numgens+1]
    y_axis_loc3_sub = loc3freq_sub[1:numgens+1]
    y_axis_num_main = num_in_main[0:numgens+1]
    y_axis_num_sub = num_in_sub[0:numgens+1]

    return y_axis_loc2_main, y_axis_loc2_sub, y_axis_loc3_main, y_axis_loc3_sub, y_axis_num_main, y_axis_num_sub

if __name__ == '__main__':
    total_runs = 96
    num_processes = 96 # number of CPU logical processors

    total_loc2_allele_freq_main = []
    total_loc2_allele_freq_sub = []
    total_loc3_allele_freq_main = []
    total_loc3_allele_freq_sub = []
    total_num_main = []
    total_num_sub = []

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for result in tqdm(pool.imap_unordered(run_simulation, range(total_runs)), total=total_runs):
            total_loc2_allele_freq_main.append(result[0])
            total_loc2_allele_freq_sub.append(result[1])
            total_loc3_allele_freq_main.append(result[2])
            total_loc3_allele_freq_sub.append(result[3])
            total_num_main.append(result[4])
            total_num_sub.append(result[5])

    param_text = (f"runs={total_runs}, numgens={numgens}, popsize={popsize}, maxsubsize={maxsubsize}, migration={migration}, mutation={mutation}, numinds={numinds}, numindssub={numindssub}, maxrepro={maxrepro}, overdom={overdom}, parthreduction={parthreduction}, parthpenalty={parthpenalty}")

    # plotting
    figure, axis = plt.subplots(2, 3, figsize=(20, 13))

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

    figure.suptitle(param_text, fontsize=12, fontweight="bold")
    plt.subplots_adjust()
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(plot_path)
    plt.show()