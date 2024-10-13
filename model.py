# Facultative Parthenogenesis Metapopulation model v1
# parth 19 adds option to have no migration, i.e. no subpopulation
# get run parameters from user

popsize=100
# popsize = int(input("Enter size of metapopulation: "))

maxsubsize=50
# maxsubsize = int(input("Enter maximum size of subpopulation: "))

migrecip=1000
# migrecip = int(input("Enter reciprocal of migration rate [enter 0 for no migration; i.e. no subpopulation]: "))
migration = 0
if (migrecip != 0): migration = 1/migrecip

mutrecip=1000
# mutrecip = int(input("Enter reciprocal of mutation rate: "))
mutation = 1/mutrecip

'''
1300 rem ? "enter recombination rate between Locus 1 and Locus 2"
1400 rem input Rec1
1500 rem ? "enter recombination rate between Locus 2 and Locus 3”
1600 rem input Rec2
'''

rec1 = 0.5
rec2 = 0.5

numinds=4
# numinds = int(input("Enter number of individuals encountered in main population: "))

numindssub=10
# numindssub = int(input("Enter number of individuals encountered in sub population: "))

maxrepro=10
# maxrepro = int(input("Enter number of offspring per female: "))

f=1
# f = int(input("Choose fitness dominance [1] or fitness overdominance [2]: "))

overdom = 0
if (f != 1): overdom = 1

numgens=5
# numgens = int(input("Enter number of generations in the run: "))

outfile="test.txt"
# outfile = input("Enter name of output file: ")

'''
2210 rem ? "enter random seed"
2220 rem input RandomSeed
'''

import random
randomseed = random.randint(0, 4999)
random.seed(randomseed)

parthreduction=0.9
# parthreduction = float(input("Enter % reproduction of Parthenogenesis vs. sexual reproduction: "))

parthrepro = int(parthreduction*maxrepro+0.5)

print("Sexual reproduction produces "+str(maxrepro)+" offspring. Parthenogenetic reproduction produces "+str(parthrepro)+" offspring.")

x = 0.9
# x = float(input("Enter % reproduction for parthenogenetic-capable females reproducing sexually: "))

parthpenalty = round(1-x, 3)

# input("Hit any key followed by <return> to continue: ")

# print all the variables
print(f"PopSize: {popsize}")
print(f"MaxSubSize: {maxsubsize}")
print(f"migRecip: {migrecip}")
print(f"Migration: {migration}")
print(f"mutRecip: {mutrecip}")
print(f"Mutation: {mutation}")
print(f"Rec1: {rec1}")
print(f"Rec2: {rec2}")
print(f"NumInds: {numinds}")
print(f"NumIndsSub: {numindssub}")
print(f"MaxRepro: {maxrepro}")
print(f"f: {f}")
print(f"OverDom: {overdom}")
print(f"NumGens: {numgens}")
print(f"Outfile: {outfile}")
print(f"RandomSeed: {randomseed}")
print(f"ParthReduction: {parthreduction}")
print(f"ParthRepro: {parthrepro}")
print(f"ParthPenalty: {parthpenalty}")

# helper function to create matrix
def createMatrix(rows, cols):
    rows += 1
    cols += 1
    return [[0 for x in range(cols)] for y in range(rows)]

# dimensioning
sex = [0] * (popsize*2 + 1)
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
popinds = [0] * (popsize + 1)
subinds = [0] * (popsize + 1)
meeting = createMatrix(popsize,numinds)
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
loc1freq = [0] * (2*numgens + 1)
loc2freq = [0] * (2*numgens + 1)
loc3freq = [0] * (2*numgens + 1)
mainpopcount = [0] * (2*numgens + 1)
subpopcount = [0] * (2*numgens + 1)
sexualoffspring = [0] * (2*numgens + 1)
parthoffspring = [0] * (2*numgens + 1)
sexuals = [0] * (2*numgens + 1)
parths = [0] * (2*numgens + 1)

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
        if (z < mutation): loc2allele1[j] = random.random()/2
        if (loc2allele1[j] > maxparth): maxparth = loc2allele1[j]

        # locus 2 allele 2
        z = random.random()
        if (z < mutation): loc2allele2[j] = random.random()/2
        if (loc2allele2[j] > maxparth): maxparth = loc2allele2[j]

        # locus 3 allele 1
        z = random.random()
        if (z < mutation): loc3allele1[j] = random.random()/2

        # locus 3 allele 2
        # no new random number z? 
        if (z < mutation): loc3allele2[j] = random.random()/2
    
    # encounter phase
    # each individual encounters NumInds or NumIndsSub individuals depending on its environment
    # for programming simplicity, self-encounters and duplicate encounters count (bad luck)
    # values of the Meeting variable are the subscripts of the encountered individual

    for a in range(1, popsize+1):
        x = numindssub
        if (location[a] == 0): x = numinds

        for b in range(1, x+1):
            while (True):
                zz = random.randint(0, popsize)
                if (location[a] == location[zz]): break
            meeting[a][b] = zz
    
    '''
    11000 rem encounters check
    11100 rem for a = 1 to PopSize
    11200 rem if Sex(a)=1 then goto 11300 else goto 11800
    11300 rem ? "female ";a; "meets ";
    11400 rem for b = 1 to NumInds
    11500 rem if Sex(Meeting(a,b))=1 then ? "female ";Meeting(a,b);" ";
    11600 rem if Sex(Meeting(a,b))=0 then ? "male ";Meeting(a,b);" ";
    11700 rem next b
    11750 rem ?
    11800 rem next a
    '''

    sexualmatingscount = 1
    parthmatingscount = 0

    # find male that each female mates with
    # for computational simplicity, each female will mate with the last male she encountered

    for a in range(1, popsize):
        matingflag = 0
        if (sex[a] == 1):
            # focal individual is female
            for b in range(1, numinds+1):
                if (sex[meeting[a][b]] == 0):
                    matingflag = 1
                    sexmothers[sexualmatingscount] = a
                    sexfathers[sexualmatingscount] = meeting[a][b]
            if (matingflag == 1): sexualmatingscount = sexualmatingscount+1
            else:
                parthmatingscount = parthmatingscount+1
                parthmothers[parthmatingscount] = a
        else: continue

    # record number of sexual and parthenogenetic reproducing females this generation
    sexuals[i] = sexualmatingscount-1
    parths[i] = parthmatingscount

    '''
    15200 rem check matings
    15300 rem for a = 1 to SexualMatingsCount-1
    15400 rem ? "female ";SexMothers(a); " mates with male ";SexFathers(a)
    15500 rem next a
    15600 rem for a = 1 to ParthMatingsCount
    15700 rem ? "female ";ParthMothers(a); " attempts parthenogenesis"
    15800 rem next a
    16000 rem reproduction phase
    16100 rem sexual reproduction
    '''

    offspringcount = 0
    max = 0
    mhaplotype = 0
    phaplotype = 0
    for a in range(1, sexualmatingscount):
        if (loc2allele1[sexmothers[a]] + loc2allele2[sexmothers[a]] == 0):
            max = maxrepro
        else:
            max = int(maxrepro*(1-parthpenalty)+0.5)
        
        for b in range(1, max+1):
            offspringcount = offspringcount+1

            # determine maternal haplotype
            x = random.randint(0,1)
            y = random.random()
            z = random.random()

            if (x == 0 and y > rec1 and z > rec2): mhaplotype = 1
            if (x == 0 and y > rec1 and z < rec2): mhaplotype = 2
            if (x == 0 and y < rec1 and z > rec2): mhaplotype = 3
            if (x == 0 and y < rec1 and z < rec2): mhaplotype = 4
            if (x == 1 and y > rec1 and z > rec2): mhaplotype = 5
            if (x == 1 and y > rec1 and z < rec2): mhaplotype = 6
            if (x == 1 and y < rec1 and z > rec2): mhaplotype = 7
            if (x == 1 and y < rec1 and z < rec2): mhaplotype = 8
            
            if (mhaplotype == 1):
                offspringloc1allele1[offspringcount] = loc1allele1[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele1[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele1[sexmothers[a]]

            if (mhaplotype == 2):
                offspringloc1allele1[offspringcount] = loc1allele1[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele1[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele2[sexmothers[a]]

            if (mhaplotype == 3):
                offspringloc1allele1[offspringcount] = loc1allele1[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele2[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele2[sexmothers[a]]

            if (mhaplotype == 4):
                offspringloc1allele1[offspringcount] = loc1allele1[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele2[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele1[sexmothers[a]]
            
            if (mhaplotype == 5):
                offspringloc1allele1[offspringcount] = loc1allele2[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele2[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele2[sexmothers[a]]
                
            if (mhaplotype == 6):
                offspringloc1allele1[offspringcount] = loc1allele2[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele2[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele1[sexmothers[a]]

            if (mhaplotype == 7):
                offspringloc1allele1[offspringcount] = loc1allele2[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele1[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele1[sexmothers[a]]

            if (mhaplotype == 8):
                offspringloc1allele1[offspringcount] = loc1allele2[sexmothers[a]]
                offspringloc2allele1[offspringcount] = loc2allele1[sexmothers[a]]
                offspringloc3allele1[offspringcount] = loc3allele2[sexmothers[a]]

            # determine paternal haplotype
            x = random.randint(0,1)
            y = random.random()
            z = random.random()

            if (x == 0 and y > rec1 and z > rec2): phaplotype = 1
            if (x == 0 and y > rec1 and z < rec2): phaplotype = 2
            if (x == 0 and y < rec1 and z > rec2): phaplotype = 3
            if (x == 0 and y < rec1 and z < rec2): phaplotype = 4
            if (x == 1 and y > rec1 and z > rec2): phaplotype = 5
            if (x == 1 and y > rec1 and z < rec2): phaplotype = 6
            if (x == 1 and y < rec1 and z > rec2): phaplotype = 7
            if (x == 1 and y < rec1 and z < rec2): phaplotype = 8

            if (phaplotype == 1):
                offspringloc1allele2[offspringcount] = loc1allele1[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele1[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele1[sexfathers[a]]

            if (phaplotype == 2):
                offspringloc1allele2[offspringcount] = loc1allele1[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele1[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele2[sexfathers[a]]

            if (phaplotype == 3):
                offspringloc1allele2[offspringcount] = loc1allele1[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele2[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele2[sexfathers[a]]

            if (phaplotype == 4):
                offspringloc1allele2[offspringcount] = loc1allele1[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele2[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele1[sexfathers[a]]
            
            if (phaplotype == 5):
                offspringloc1allele2[offspringcount] = loc1allele2[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele2[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele2[sexfathers[a]]
                
            if (phaplotype == 6):
                offspringloc1allele2[offspringcount] = loc1allele2[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele2[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele1[sexfathers[a]]

            if (phaplotype == 7):
                offspringloc1allele2[offspringcount] = loc1allele2[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele1[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele1[sexfathers[a]]

            if (phaplotype == 8):
                offspringloc1allele2[offspringcount] = loc1allele2[sexfathers[a]]
                offspringloc2allele2[offspringcount] = loc2allele1[sexfathers[a]]
                offspringloc3allele2[offspringcount] = loc3allele2[sexfathers[a]]

            # set other offspring parameters
            offspringalive[offspringcount] = 1
            zz = random.random() 
            if (zz < 0.5): offspringsex[offspringcount] = 1
            else: offspringsex[offspringcount] = 0
            if (location[sexmothers[a]] == 1): offspringlocation[offspringcount] = 1
            else: offspringlocation[offspringcount] = 0

    # save current OffspringCount number as the number of offspring produced by sexual reproduction this generation
    sexualoffspring[i] = offspringcount
    # parthenogenic reproduction
    # produce parthogenetic offspring by gamete cell duplication

    mhaplotype = 0 # ensuring reset of mhaplotype
    for a in range(1, parthmatingscount+1):
        zz = random.random()
        if (zz < (loc2allele1[parthmothers[a]] + loc2allele2[parthmothers[a]])):
            # female will reproduce parthenogenetic offspring
            for b in range(1, parthrepro+1):
                offspringcount = offspringcount+1

                # determine maternal haplotype
                x = random.randint(0,1)
                y = random.random()
                z = random.random()

                if (x == 0 and y > rec1 and z > rec2): mhaplotype = 1
                if (x == 0 and y > rec1 and z < rec2): mhaplotype = 2
                if (x == 0 and y < rec1 and z > rec2): mhaplotype = 3
                if (x == 0 and y < rec1 and z < rec2): mhaplotype = 4
                if (x == 1 and y > rec1 and z > rec2): mhaplotype = 5
                if (x == 1 and y > rec1 and z < rec2): mhaplotype = 6
                if (x == 1 and y < rec1 and z > rec2): mhaplotype = 7
                if (x == 1 and y < rec1 and z < rec2): mhaplotype = 8

                # parthenogenetic offspring represent duplicated maternal meiotic haplotypes

                if (mhaplotype == 1):
                    offspringloc1allele1[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele1[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele1[parthmothers[a]]

                if (mhaplotype == 2):
                    offspringloc1allele1[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele2[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele2[parthmothers[a]]

                if (mhaplotype == 3):
                    offspringloc1allele1[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele2[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele2[parthmothers[a]]

                if (mhaplotype == 4):
                    offspringloc1allele1[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele1[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele1[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele1[parthmothers[a]]

                if (mhaplotype == 5):
                    offspringloc1allele1[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele2[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele2[parthmothers[a]]

                if (mhaplotype == 6):
                    offspringloc1allele1[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele1[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele2[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele1[parthmothers[a]]

                if (mhaplotype == 7):
                    offspringloc1allele1[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele1[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele1[parthmothers[a]]

                if (mhaplotype == 8):
                    offspringloc1allele1[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele1[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele1[offspringcount] = loc3allele2[parthmothers[a]]
                    offspringloc1allele2[offspringcount] = loc1allele2[parthmothers[a]]
                    offspringloc2allele2[offspringcount] = loc2allele1[parthmothers[a]]
                    offspringloc3allele2[offspringcount] = loc3allele2[parthmothers[a]]
                
                # set other offspring parameters
                offspringalive[offspringcount] = 1
                zz = random.random()
                if (zz < 0.5): offspringsex[offspringcount] = 1
                else: offspringsex[offspringcount] = 0
                if (location[parthmothers[a]] == 1): offspringlocation[offspringcount] = 1
                else: offspringlocation[offspringcount] = 0
            
        else: continue
    
    # set number of offspring produced in parthenogenetic loop as the number of parthenogenetic offspring produced this generation
    parthoffspring[i] = offspringcount - sexualoffspring[i]
    # randomly sort offspring pool
    # assign RND (0,1) to each member of OffspringCount

    for a in range(1, offspringcount+1):
        score[a] = random.random()
        originalindex[a] = a

    # sort OffspringCount subscripts by the random number and get new subscript for each member
    
    # changed bubble sort to built-in sort (nlogn)
    # custom code to maintain scores and original index
    combined = []
    for a in range(1, offspringcount+1):
        combined.append((score[a], originalindex[a]))
    
    combined.sort(key=lambda x: x[0])

    for a in range(0, offspringcount):
        curr = combined[a]
        score[a+1] = curr[0]
        originalindex[a+1] = curr[1]
    # custom code end

    '''
    43000 rem Print the sorted array and original indices
    43100 rem PRINT "Sorted Scores and Original Indices:"
    43200 rem FOR i = 1 TO OffspringCount
    43300 rem   PRINT "Score: "; score(i); " Original Index: "; originalIndex(i)
    43400 rem NEXT i
    43500 rem set new reordered offspring array
    '''

    for a in range(1, offspringcount+1):
        newoffspringalive[a] = offspringalive[originalindex[a]]
        newoffspringsex[a] = offspringsex[originalindex[a]]
        newoffspringlocation[a] = offspringlocation[originalindex[a]]
        newoffspringloc1allele1[a] = offspringloc1allele1[originalindex[a]]
        newoffspringloc1allele2[a] = offspringloc1allele2[originalindex[a]]
        newoffspringloc2allele1[a] = offspringloc2allele1[originalindex[a]]
        newoffspringloc2allele2[a] = offspringloc2allele2[originalindex[a]]
        newoffspringloc3allele1[a] = offspringloc3allele1[originalindex[a]]
        newoffspringloc3allele2[a] = offspringloc3allele2[originalindex[a]]

    # selection phase
    for a in range(1, offspringcount+1):

        if (overdom == 0):
            # fitness dominance
            if (newoffspringloc1allele1[a] + newoffspringloc1allele2[a] == 0): newoffspringalive[a] = 0
    
        else:
            # fitness overdominance
            if (newoffspringloc1allele1[a] + newoffspringloc1allele2[a] == 0): newoffspringalive[a] = 0
            if (newoffspringloc1allele1[a] + newoffspringloc1allele2[a] == 2): newoffspringalive[a] = 0

    # simultaneously determine migration and send living offspring into adult pool for next generation
    # need to keep count of number of individuals in the main and sub populations

    maincount = 0
    subcount = 0
    offspringcounter = 0

    for a in range(1, popsize+1):
        offspringcounter = offspringcounter+1
        while True:
            if (newoffspringalive[offspringcounter] == 0): offspringcounter = offspringcounter+1
            else: break

        z = random.random()
        if (z < migration): newoffspringlocation[offspringcounter] = abs(newoffspringlocation[offspringcounter]-1)

        if (newoffspringlocation[offspringcounter] == 0): maincount = maincount+1
        else: subcount = subcount+1

        if (subcount > maxsubsize): newoffspringlocation[offspringcounter] = 0 

        # Migration determination is complete; now the individual goes into the adult pool for the next generation
        sex[a] = newoffspringsex[offspringcounter]
        location[a] = newoffspringlocation[offspringcounter]
        loc1allele1[a] = newoffspringloc1allele1[offspringcounter]
        loc1allele2[a] = newoffspringloc1allele2[offspringcounter]
        loc2allele1[a] = newoffspringloc2allele1[offspringcounter]
        loc2allele2[a] = newoffspringloc2allele2[offspringcounter]
        loc3allele1[a] = newoffspringloc3allele1[offspringcounter]
        loc3allele2[a] = newoffspringloc3allele2[offspringcounter]
    
    '''
    49910 rem check genotypes
    49920 rem for a = 1 to PopSize
    49930 rem ? "ind# ";a;" ";
    49940 rem ? "sex ";Sex(a);" ";
    49950 rem ? "loc ";Location(a);" ";
    49960 rem ? "genotype ";Loc1Allele1(a);Loc1Allele2(a);" ";Loc2Allele1(a);Loc2Allele2(a);" ";Loc3Allele1(a);Loc3Allele2(a)
    49965 rem next a
    49970 rem ? "hit any key then [return] to continue"
    49980 rem input x$
    50000 rem compute and store generation data
    50100 rem calculate allele frequencies at each locus for this generation
    '''

    print()
    maincount = 0
    subcount = 0
    print("generation "+str(i))
    loc1count = 0
    loc2count = 0
    loc3count = 0
    for a in range(1, popsize+1):
        loc1count = loc1count+loc1allele1[a]+loc1allele2[a]
        loc2count = loc2count+loc2allele1[a]+loc2allele2[a]
        loc3count = loc3count+loc3allele1[a]+loc3allele2[a]
        if (location[a] == 0): maincount = maincount+1
        if (location[a] == 1): subcount = subcount+1
    
    loc1freq[i] = loc1count/(2*popsize)
    print("mean fitness (determined by Locus 1) = "+str(loc1freq[i]))
    loc2freq[i] = loc2count/popsize
    print("Mean Parthenogenetic Capability (determined by Locus 2) = "+str(loc2freq[i]))
    loc3freq[i] = loc3count/popsize
    print("Locus 3 mean (neutral locus) = "+str(loc3freq[i]))
    print("main population: "+str(maincount)+" individuals")
    print("subpopulation: "+str(subcount)+" individuals")
    print(str(sexuals[i])+" females reproduced sexually")
    print(str(sexualoffspring[i])+" sexual offspring produced")
    print(str(parths[i])+" females attempted to reproduce parthenogenetically")
    print(str(parthoffspring[i])+" parthenogenetic offspring produced")
    
    if (overdom == 1 and parthoffspring[i] > 0): print("(but these die due to homoz. Lethality)")
    else: print("overdom==0 OR parthoffspring[i] >= 0")
    print("Maximum value of parthenogenetic allele that has appeared: "+str(maxparth))

# print all data to output file
# 
# 
#