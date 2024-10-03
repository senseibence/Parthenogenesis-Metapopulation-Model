# Facultative Parthenogenesis Metapopulation model v1
# parth 19 adds option to have no migration, i.e. no subpopulation
# get run parameters from user

# PopSize=100
PopSize = int(input("Enter size of metapopulation: "))

# MaxSubSize=50
MaxSubSize = int(input("Enter maximum size of subpopulation: "))

# migRecip=1000
migRecip = int(input("Enter reciprocal of migration rate [enter 0 for no migration; i.e. no subpopulation]: "))
Migration = 0
if (migRecip != 0): Migration = 1/migRecip

# mutRecip=1000
mutRecip = int(input("Enter reciprocal of mutation rate: "))
Mutation = 1/mutRecip

'''
rem ? "enter recombination rate between Locus 1 and Locus 2"
rem input Rec1
rem ? "enter recombination rate between Locus 2 and Locus 3”
rem input Rec2
'''

Rec1 = 0.5
Rec2 = 0.5

# NumInds=4
NumInds = int(input("Enter number of individuals encountered in main population: "))

# NumIndsSub=10
NumIndsSub = int(input("Enter number of individuals encountered in sub population: "))

# MaxRepro=10
MaxRepro = int(input("Enter number of offspring per female: "))

# f=1
f = int(input("Choose fitness dominance [1] or fitness overdominance [2]: "))

OverDom = 0
if (f != 1): OverDom = 1

# NumGens=5
NumGens = int(input("Enter number of generations in the run: "))

# Outfile$="outfile"
Outfile = input("Enter name of output file: ")

'''
rem ? "enter random seed"
rem input RandomSeed
'''

import random
RandomSeed = random.randint(0, 4999)
random.seed(RandomSeed)

ParthReduction = float(input("Enter % reproduction of Parthenogenesis vs. sexual reproduction: "))

ParthRepro = int(ParthReduction*MaxRepro+0.5)

print("Sexual reproduction produces "+str(MaxRepro)+" offspring. Parthenogenetic reproduction produces "+str(ParthRepro)+" offspring.")

x = float(input("Enter % reproduction for parthenogenetic-capable females reproducing sexually: "))

ParthPenalty = round(1-x, 3)

input("Hit any key followed by <return> to continue: ")

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

# dimensioning
