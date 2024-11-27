# Parthenogenesis Metapopulation Model v5
# tracking main and subpopulation separately
# made bug fix in mating phase

import random
import multiprocessing
from tqdm import tqdm
import matplotlib.pyplot as plt

# helper function to create a matrix
def createMatrix(rows, cols):
    rows += 1
    cols += 1
    return [[0 for x in range(cols)] for y in range(rows)]

# simulation function for each run
def run_simulation(run):
    
    # parameters
    popsize=5000
    maxsubsize=100
    migration = (1/popsize)
    mutrecip=10000000
    mutation = 1/mutrecip
    rec1 = 0.5
    rec2 = 0.5
    numinds=10
    numindssub=min(int(0.04*maxsubsize), numinds//2)
    maxrepro=10
    f=1
    overdom = 0
    if (f != 1): overdom = 1
    numgens=2000
    parthreduction=0.2
    parthrepro = int(parthreduction*maxrepro+0.5)
    parthpenalty = 0.25

    randomseed = random.randint(0, 4999) + run  # ensure different seeds
    random.seed(randomseed)

    # dimensioning
    sex = [0] * (popsize*2 + 1)
    sex[0] = -1  # prevents unintended mating
    location = [0] * (popsize + 1)
    loc1allele1 = [0] * (popsize + 1)
    loc1allele2 = [0] * (popsize + 1)
    loc2allele1 = [0] * (popsize + 1)
    loc2allele2 = [0] * (popsize + 1)
    loc3allele1 = [0] * (popsize + 1)
    loc3allele2 = [0] * (popsize + 1)
    offspringalive = [0] * (maxrepro*popsize + 1)
    offspringsex = [0] * (maxrepro*popsize + 1)
    offspringlocation = [0] * (maxrepro*popsize + 1)
    offspringloc1allele1 = [0] * (maxrepro*popsize + 1)
    offspringloc1allele2 = [0] * (maxrepro*popsize + 1)
    offspringloc2allele1 = [0] * (maxrepro*popsize + 1)
    offspringloc2allele2 = [0] * (maxrepro*popsize + 1)
    offspringloc3allele1 = [0] * (maxrepro*popsize + 1)
    offspringloc3allele2 = [0] * (maxrepro*popsize + 1)
    sexualoffspring = [0] * (2*numgens + 1)
    parthoffspring = [0] * (2*numgens + 1)
    sexuals = [0] * (2*numgens + 1)
    parths = [0] * (2*numgens + 1)
    loc1freq = [0] * (2*numgens + 1)
    loc2freq = [0] * (2*numgens + 1)
    loc3freq = [0] * (2*numgens + 1)

    # lists for tracking main and subpopulation separately
    loc2freq_main = [0] * (2*numgens + 1)
    loc2freq_sub = [0] * (2*numgens + 1)

    meeting = createMatrix(popsize, numinds)
    sexmothers = [0] * (popsize + 1)
    sexfathers = [0] * (popsize + 1)
    parthmothers = [0] * (popsize + 1)
    score = [0] * (maxrepro*popsize + 1)
    originalindex = [0] * (maxrepro*popsize + 1)
    newoffspringalive = [0] * (maxrepro*popsize + 1)
    newoffspringsex = [0] * (maxrepro*popsize + 1)
    newoffspringlocation = [0] * (maxrepro*popsize + 1)
    newoffspringloc1allele1 = [0] * (maxrepro*popsize + 1)
    newoffspringloc1allele2 = [0] * (maxrepro*popsize + 1)
    newoffspringloc2allele1 = [0] * (maxrepro*popsize + 1)
    newoffspringloc2allele2 = [0] * (maxrepro*popsize + 1)
    newoffspringloc3allele1 = [0] * (maxrepro*popsize + 1)
    newoffspringloc3allele2 = [0] * (maxrepro*popsize + 1)

    # keep track of the maximum value of locus 2 allele that has appeared
    maxparth = 0

    # set initial arrays
    for i in range(1, popsize+1):
        aa = random.random()
        if (aa < 0.5): sex[i] = 1
        else: sex[i] = 0
        location[i] = 0
        loc1allele1[i] = 1
        if (overdom == 1): loc1allele2[i] = 0
        else: loc1allele2[i] = 1
        loc2allele1[i] = 0
        loc2allele2[i] = 0
        loc3allele1[i] = 0
        loc3allele2[i] = 0

    # start main loop
    for i in range(1, numgens+1):
        # mutation phase
        for j in range(1, popsize+1):
            # locus 1 allele 1
            z = random.random()
            if (z < mutation): loc1allele1[j] = abs(loc1allele1[j] - 1)

            # locus 1 allele 2
            z = random.random()
            if (z < mutation): loc1allele2[j] = abs(loc1allele2[j] - 1)

            # locus 2 allele 1
            z = random.random()
            if (z < mutation):
                loc2allele1[j] = 0.5
            if (loc2allele1[j] > maxparth): maxparth = loc2allele1[j]

            # locus 2 allele 2
            z = random.random()
            if (z < mutation):
                loc2allele2[j] = 0.5
            if (loc2allele2[j] > maxparth): maxparth = loc2allele2[j]

            # locus 3 allele 1
            z = random.random()
            if (z < mutation): loc3allele1[j] = 0.5

            # locus 3 allele 2
            z = random.random()
            if (z < mutation): loc3allele2[j] = 0.5

        # encounter phase
        for a in range(1, popsize+1):
            x = numindssub if location[a] == 1 else numinds

            for b in range(1, x+1):
                temp_list = []

                while True:
                    if len(temp_list) >= popsize:
                        break
                    zz = random.randint(1, popsize)
                    if a != zz and location[a] == location[zz]:
                        break
                    if zz not in temp_list:
                        temp_list.append(zz)

                if len(temp_list) < popsize:
                    if b <= len(meeting[a]) - 1:
                        meeting[a][b] = zz

        sexualmatingscount = 1
        parthmatingscount = 0

        # find male that each female mates with
        for a in range(1, popsize):
            matingflag = 0
            if sex[a] == 1:
                x = numindssub if location[a] == 1 else numinds # redundant because sex[0] = -1

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

        # record number of sexual and parthenogenetic reproducing females this generation
        sexuals[i] = sexualmatingscount - 1
        parths[i] = parthmatingscount

        # reproduction phase
        offspringcount = 0
        max_offspring = 0
        mhaplotype = 0
        phaplotype = 0

        # sexual reproduction
        for a in range(1, sexualmatingscount):
            if loc2allele1[sexmothers[a]] + loc2allele2[sexmothers[a]] == 0:
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

                alleles = m_alleles[mhaplotype]
                offspringloc1allele1[offspringcount] = alleles[0][sexmothers[a]]
                offspringloc2allele1[offspringcount] = alleles[1][sexmothers[a]]
                offspringloc3allele1[offspringcount] = alleles[2][sexmothers[a]]

                # determine paternal haplotype
                x = random.randint(0,1)
                y = random.random()
                z = random.random()

                # assign phaplotype based on recombination
                phaplotype = (x << 2) | ((y < rec1) << 1) | (z < rec2)

                # paternal alleles
                alleles = m_alleles[phaplotype]
                offspringloc1allele2[offspringcount] = alleles[0][sexfathers[a]]
                offspringloc2allele2[offspringcount] = alleles[1][sexfathers[a]]
                offspringloc3allele2[offspringcount] = alleles[2][sexfathers[a]]

                # set other offspring parameters
                offspringalive[offspringcount] = 1
                offspringsex[offspringcount] = 1 if random.random() < 0.5 else 0
                offspringlocation[offspringcount] = location[sexmothers[a]]

        # record number of sexual offspring
        sexualoffspring[i] = offspringcount

        # parthenogenetic reproduction
        for a in range(1, parthmatingscount+1):
            if random.random() < (loc2allele1[parthmothers[a]] + loc2allele2[parthmothers[a]]):
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
                    offspringloc1allele1[offspringcount] = alleles[0][parthmothers[a]]
                    offspringloc1allele2[offspringcount] = alleles[0][parthmothers[a]]
                    offspringloc2allele1[offspringcount] = alleles[1][parthmothers[a]]
                    offspringloc2allele2[offspringcount] = alleles[1][parthmothers[a]]
                    offspringloc3allele1[offspringcount] = alleles[2][parthmothers[a]]
                    offspringloc3allele2[offspringcount] = alleles[2][parthmothers[a]]

                    # set other offspring parameters
                    offspringalive[offspringcount] = 1
                    offspringsex[offspringcount] = 1 if random.random() < 0.5 else 0
                    offspringlocation[offspringcount] = location[parthmothers[a]]

        # record number of parthenogenetic offspring
        parthoffspring[i] = offspringcount - sexualoffspring[i]

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

        # migration and adult pool for next generation
        maincount = 0
        subcount = 0
        offspringcounter = 0

        for a in range(1, popsize+1):
            offspringcounter += 1
            while offspringcounter <= offspringcount and newoffspringalive[offspringcounter] == 0:
                offspringcounter += 1

            if offspringcounter > offspringcount:
                break

            # migration
            if random.random() < migration:
                newoffspringlocation[offspringcounter] = 1 - newoffspringlocation[offspringcounter]

            if newoffspringlocation[offspringcounter] == 0:
                maincount += 1
            else:
                subcount += 1

            if subcount > maxsubsize:
                newoffspringlocation[offspringcounter] = 0

            # update adult pool
            sex[a] = newoffspringsex[offspringcounter]
            location[a] = newoffspringlocation[offspringcounter]
            loc1allele1[a] = newoffspringloc1allele1[offspringcounter]
            loc1allele2[a] = newoffspringloc1allele2[offspringcounter]
            loc2allele1[a] = newoffspringloc2allele1[offspringcounter]
            loc2allele2[a] = newoffspringloc2allele2[offspringcounter]
            loc3allele1[a] = newoffspringloc3allele1[offspringcounter]
            loc3allele2[a] = newoffspringloc3allele2[offspringcounter]

        # compute and store generation data
        loc1count = sum(loc1allele1[1:popsize+1]) + sum(loc1allele2[1:popsize+1])
        loc2count = sum(loc2allele1[1:popsize+1]) + sum(loc2allele2[1:popsize+1])
        loc3count = sum(loc3allele1[1:popsize+1]) + sum(loc3allele2[1:popsize+1])

        loc1freq[i] = loc1count / (2 * popsize)
        loc2freq[i] = loc2count / (2 * popsize)
        loc3freq[i] = loc3count / (2 * popsize)

        num_main = 0
        num_sub = 0
        loc2count_main = 0
        loc2count_sub = 0

        for a in range(1, popsize+1):
            if (location[a] == 0):
                num_main += 1
                loc2count_main += loc2allele1[a] + loc2allele2[a]
            elif (location[a] == 1):
                num_sub += 1
                loc2count_sub += loc2allele1[a] + loc2allele2[a]

        if num_main > 0: loc2freq_main[i] = loc2count_main / (2 * num_main)   
        else: loc2freq_main[i] = 0
        if num_sub > 0: loc2freq_sub[i] = loc2count_sub / (2 * num_sub)
        else: loc2freq_sub[i] = 0
            
    # prepare data for plotting
    y_axis_loc2 = loc2freq[1:numgens+1]
    y_axis_loc3 = loc3freq[1:numgens+1]
    y_axis_loc2_main = loc2freq_main[1:numgens+1]
    y_axis_loc2_sub = loc2freq_sub[1:numgens+1]

    return y_axis_loc2, y_axis_loc3, y_axis_loc2_main, y_axis_loc2_sub

if __name__ == '__main__':
    total_runs = 32
    num_processes = 32  # number of CPU logical processors
    numgens = 2000

    total_loc2_allele_freq = []
    total_loc3_allele_freq = []
    total_loc2_allele_freq_main = []
    total_loc2_allele_freq_sub = []

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for result in tqdm(pool.imap_unordered(run_simulation, range(total_runs)), total=total_runs):
            total_loc2_allele_freq.append(result[0])
            total_loc3_allele_freq.append(result[1])
            total_loc2_allele_freq_main.append(result[2])
            total_loc2_allele_freq_sub.append(result[3])

    param_text = f"parameters: popsize=5000, maxsubsize=100, numinds=10, maxrepro=10, numgens=2000, parthreduction=0.2, parthpenalty=0.25, migration=(1/popsize)"

    # plotting
    figure, axis = plt.subplots(2, 2, figsize=(16, 10))

    # plot 1
    axis[0][0].set_xlim(0, numgens)
    axis[0][0].set_ylim(0, 0.6)
    axis[0][0].set_xlabel("Generations")
    axis[0][0].set_ylabel("Frequency")
    axis[0][0].set_title("Average Locus 2 (Parthenogenetic) Allele Frequency")
    for run in range(total_runs):
        axis[0][0].plot(total_loc2_allele_freq[run], color=plt.cm.rainbow(run / total_runs))

    # plot 2
    axis[0][1].set_xlim(0, numgens)
    axis[0][1].set_ylim(0, 0.6)
    axis[0][1].set_xlabel("Generations")
    axis[0][1].set_ylabel("Frequency")
    axis[0][1].set_title("Average Locus 3 (Neutral) Allele Frequency")
    for run in range(total_runs):
        axis[0][1].plot(total_loc3_allele_freq[run], color=plt.cm.rainbow(run / total_runs))

    # plot 3
    axis[1][0].set_xlim(0, numgens)
    axis[1][0].set_ylim(0, 0.6)
    axis[1][0].set_xlabel("Generations")
    axis[1][0].set_ylabel("Frequency")
    axis[1][0].set_title("Average Locus 2 (Parthenogenetic) Allele Frequency Main Population")
    for run in range(total_runs):
        axis[1][0].plot(total_loc2_allele_freq_main[run], color=plt.cm.rainbow(run / total_runs))

    # plot 4
    axis[1][1].set_xlim(0, numgens)
    axis[1][1].set_ylim(0, 0.6)
    axis[1][1].set_xlabel("Generations")
    axis[1][1].set_ylabel("Frequency")
    axis[1][1].set_title("Average Locus 2 (Parthenogenetic) Allele Frequency Subpopulation")
    for run in range(total_runs):
        axis[1][1].plot(total_loc2_allele_freq_sub[run], color=plt.cm.rainbow(run / total_runs))

    figure.suptitle(param_text, fontsize=12)
    plt.subplots_adjust()  
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()