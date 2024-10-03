# Facultative Parthenogenesis Metapopulation model v1
# parth 19 adds option to have no migration, i.e. no subpopulation
# get run parameters from user

PopSize=100
# PopSize = int(input("Enter size of metapopulation: "))

MaxSubSize=50
# MaxSubSize = int(input("Enter maximum size of subpopulation: "))

migRecip=1000
# migRecip = int(input("Enter reciprocal of migration rate [enter 0 for no migration; i.e. no subpopulation]: "))
Migration = 0
if (migRecip != 0): Migration = 1/migRecip

mutRecip=1000
# mutRecip = int(input("Enter reciprocal of mutation rate: "))
Mutation = 1/mutRecip

'''
rem ? "enter recombination rate between Locus 1 and Locus 2"
rem input Rec1
rem ? "enter recombination rate between Locus 2 and Locus 3”
rem input Rec2
'''

Rec1 = 0.5
Rec2 = 0.5

NumInds=4
# NumInds = int(input("Enter number of individuals encountered in main population: "))

NumIndsSub=10
# NumIndsSub = int(input("Enter number of individuals encountered in sub population: "))

MaxRepro=10
# MaxRepro = int(input("Enter number of offspring per female: "))

f=1
# f = int(input("Choose fitness dominance [1] or fitness overdominance [2]: "))

OverDom = 0
if (f != 1): OverDom = 1

NumGens=5
# NumGens = int(input("Enter number of generations in the run: "))

Outfile="test.txt"
# Outfile = input("Enter name of output file: ")

'''
rem ? "enter random seed"
rem input RandomSeed
'''

import random
RandomSeed = random.randint(0, 4999)
random.seed(RandomSeed)

ParthReduction=0.9
# ParthReduction = float(input("Enter % reproduction of Parthenogenesis vs. sexual reproduction: "))

ParthRepro = int(ParthReduction*MaxRepro+0.5)

print("Sexual reproduction produces "+str(MaxRepro)+" offspring. Parthenogenetic reproduction produces "+str(ParthRepro)+" offspring.")

x = 0.9
# x = float(input("Enter % reproduction for parthenogenetic-capable females reproducing sexually: "))

ParthPenalty = round(1-x, 3)

# input("Hit any key followed by <return> to continue: ")

# print all the variables
print(f"PopSize: {PopSize}")
print(f"MaxSubSize: {MaxSubSize}")
print(f"migRecip: {migRecip}")
print(f"Migration: {Migration}")
print(f"mutRecip: {mutRecip}")
print(f"Mutation: {Mutation}")
print(f"Rec1: {Rec1}")
print(f"Rec2: {Rec2}")
print(f"NumInds: {NumInds}")
print(f"NumIndsSub: {NumIndsSub}")
print(f"MaxRepro: {MaxRepro}")
print(f"f: {f}")
print(f"OverDom: {OverDom}")
print(f"NumGens: {NumGens}")
print(f"Outfile: {Outfile}")
print(f"RandomSeed: {RandomSeed}")
print(f"ParthReduction: {ParthReduction}")
print(f"ParthRepro: {ParthRepro}")
print(f"ParthPenalty: {ParthPenalty}")

# helper function to create matrix
def createMatrix(rows, cols):
    return [[0 for x in range(cols)] for y in range(rows)]

# dimensioning
Sex = [0] * (PopSize*2 + 1)
Location = [0] * (PopSize + 1)
Loc1Allele1 = [0] * (PopSize + 1)
Loc1Allele2 = [0] * (PopSize + 1)
Loc2Allele1 = [0] * (PopSize + 1)
Loc2Allele2 = [0] * (PopSize + 1)
Loc3Allele1 = [0] * (PopSize + 1)
Loc3Allele2 = [0] * (PopSize + 1)
OffspringAlive = [0] * (MaxRepro*PopSize + 1)
OffspringSex = [0] * (MaxRepro*PopSize + 1)
OffspringLocation = [0] * (MaxRepro*PopSize + 1)
OffspringLoc1Allele1 = [0] * (MaxRepro*PopSize + 1)
OffspringLoc1Allele2 = [0] * (MaxRepro*PopSize + 1)
OffspringLoc2Allele1 = [0] * (MaxRepro*PopSize + 1)
OffspringLoc2Allele2 = [0] * (MaxRepro*PopSize + 1)
OffspringLoc3Allele1 = [0] * (MaxRepro*PopSize + 1)
OffspringLoc3Allele2 = [0] * (MaxRepro*PopSize + 1)
PopInds = [0] * (PopSize + 1)
SubInds = [0] * (PopSize + 1)
Meeting = createMatrix(PopSize,NumInds)
SexMothers = [0] * (PopSize + 1)
SexFathers = [0] * (PopSize + 1)
ParthMothers = [0] * (PopSize + 1)
score = [0] * (MaxRepro*PopSize + 1)
OriginalIndex = [0] * (MaxRepro*PopSize + 1)
NewOffspringAlive = [0] * (MaxRepro*PopSize + 1)
NewOffspringSex = [0] * (MaxRepro*PopSize + 1)
NewOffspringLocation = [0] * (MaxRepro*PopSize + 1)
NewOffspringLoc1Allele1 = [0] * (MaxRepro*PopSize + 1)
NewOffspringLoc1Allele2 = [0] * (MaxRepro*PopSize + 1)
NewOffspringLoc2Allele1 = [0] * (MaxRepro*PopSize + 1)
NewOffspringLoc2Allele2 = [0] * (MaxRepro*PopSize + 1)
NewOffspringLoc3Allele1 = [0] * (MaxRepro*PopSize + 1)
NewOffspringLoc3Allele2 = [0] * (MaxRepro*PopSize + 1)
Loc1Freq = [0] * (2*NumGens + 1)
Loc2Freq = [0] * (2*NumGens + 1)
Loc3Freq = [0] * (2*NumGens + 1)
MainPopCount = [0] * (2*NumGens + 1)
SubPopCount = [0] * (2*NumGens + 1)
SexualOffspring = [0] * (2*NumGens + 1)
ParthOffspring = [0] * (2*NumGens + 1)
Sexuals = [0] * (2*NumGens + 1)
Parths = [0] * (2*NumGens + 1)
ParthMode = [0] * (2*NumGens + 1)
frequency = [0] * (2*PopSize + 1)
ModeFreq = [0] * (2*NumGens + 1)

# keep track of the maximum value of locus 2 allele that has appeared
MaxParth = 0

# set initial arrays
for i in range(1, PopSize+1):
    aa = random.random()
    if (aa < 0.5): Sex[i] = 1
    else: Sex[i] = 0
    Location[i] = 0
    Loc1Allele1[i] = 1
    if (OverDom == 1): Loc1Allele2[i] = 0
    else: Loc1Allele2[i] = 1
    Loc2Allele1[i] = 0
    Loc2Allele2[i] = 0
    Loc3Allele1[i] = 0
    Loc3Allele2[i] = 0

# start main loop
for i in range(1, NumGens+1):
    pass