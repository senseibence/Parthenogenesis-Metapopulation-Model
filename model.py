# Facultative Parthenogenesis Metapopulation model v1
# parth 19 adds option to have no migration, i.e. no subpopulation
# get run parameters from user

# PopSize=100
popsize = int(input("Enter size of metapopulation: "))

# MaxSubSize=50
maxsubsize = int(input("Enter maximum size of subpopulation: "))

# migRecip=1000
migrecip = int(input("Enter reciprocal of migration rate [enter 0 for no migration; i.e. no subpopulation]: "))
migration = 0
if (migrecip != 0): migration = 1/migrecip

# mutRecip=1000
mutrecip = int(input("Enter reciprocal of mutation rate: "))
mutation = 1/mutrecip

'''
rem ? "enter recombination rate between Locus 1 and Locus 2"
rem input Rec1
rem ? "enter recombination rate between Locus 2 and Locus 3”
rem input Rec2
'''

rec1 = 0.5
rec2 = 0.5

# NumInds=4
numinds = int(input("Enter number of individuals encountered in main population: "))

numindssub = int(input("Enter number of individuals encountered in sub population: "))

# MaxRepro=10
maxrepro = int(input("Enter number of offspring per female: "))

# f=1
f = int(input("Choose fitness dominance [1] or fitness overdominance [2]: "))

overdom = 0
if (f != 1): overdom = 1

# NumGens=5
numgens = int(input("Enter number of generations in the run: "))

# Outfile$="outfile"
outfile = input("Enter name of output file: ")

'''
rem ? "enter random seed"
rem input RandomSeed
'''

import random
randomseed = random.randint(0, 4999)
random.seed(randomseed)

parthreduction = float(input("Enter % reproduction of Parthenogenesis vs. sexual reproduction: "))

parthrepro = int(parthreduction*maxrepro+0.5)

print("Sexual reproduction produces "+str(maxrepro)+" offspring. Parthenogenetic reproduction produces "+str(parthrepro)+" offspring.")

x = float(input("Enter % reproduction for parthenogenetic-capable females reproducing sexually: "))

parthpenalty = 1-x

input("Hit any key followed by <return> to continue: ")

# dimensioning

print(popsize, maxsubsize, migrecip, mutrecip, rec1, rec2, numinds, numindssub, maxrepro, f, numgens, outfile, randomseed, parthreduction, parthrepro, parthpenalty)