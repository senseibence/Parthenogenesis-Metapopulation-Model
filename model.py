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
rem ? "enter recombination rate between Locus 1 and Locus 2"
rem input Rec1
rem ? "enter recombination rate between Locus 2 and Locus 3”
rem input Rec2
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
rem ? "enter random seed"
rem input RandomSeed
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
Meeting = createMatrix(popsize,numinds)
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
    pass