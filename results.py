import os
import re
import numpy as np
import pandas as pd
from collections import defaultdict

popsize = "50k"
root = f"Outputs/popsize_{popsize}_gens_10k"
arr_path_main_global = "loc2_freq_main.npy"
arr_path_sub_global = "loc2_freq_sub.npy"
results_path = f"Results/popsize_{popsize}_gens_10k.xlsx"

groups = defaultdict(list)
pattern = re.compile(r"sim(\d+) ([01]{5})$")
for directory in os.listdir(root):
    match = pattern.match(directory)
    if (not match): continue
    num, bits = match.groups()
    groups[bits].append(num)

def process(run):
    res = "no_parth"
    for freq in run:
        if (freq >= 0.475): 
            res = "max_parth"
            break
        elif (freq >= 0.1): res = "some_parth"
    return res

def process_all_allele_freq(total_allele_freq):
    no_parth = 0
    some_parth = 0
    max_parth = 0
    start = int(0.9*len(total_allele_freq[0]))
    for run in total_allele_freq:
        truncated_run = run[start:]
        res = process(truncated_run)
        if (res == "no_parth"): no_parth += 1
        if (res == "some_parth"): some_parth += 1
        if (res == "max_parth"): max_parth += 1
    return (no_parth, some_parth, max_parth)

results = {}
for bits, nums in groups.items():
    total_counts_main = np.zeros(3, dtype=int)
    total_counts_sub = np.zeros(3, dtype=int)

    for num in nums:
        directory = os.path.join(root, f"sim{num} {bits}")
        arr_path_main = os.path.join(directory, arr_path_main_global)
        arr_path_sub = os.path.join(directory, arr_path_sub_global)
        arr_main = np.load(arr_path_main)
        arr_sub = np.load(arr_path_sub)
        main_counts = process_all_allele_freq(arr_main)
        sub_counts = process_all_allele_freq(arr_sub)
        total_counts_main += main_counts
        total_counts_sub += sub_counts

    main_sum = total_counts_main.sum()
    sub_sum = total_counts_sub.sum()
    proportions_main = tuple(total_counts_main/main_sum)
    proportions_sub = tuple(total_counts_sub/sub_sum)
    results[bits] = (proportions_main, proportions_sub)

rows = []
for bits, (prop_main, prop_sub) in sorted(results.items()):
    high_low = ["high" if bit == "1" else "low" for bit in bits]
    rows.append({
        "subpop size": high_low[0],
        "mutation": high_low[1],
        "encounter main": high_low[2],
        "encounter sub": high_low[3],
        "penalty": high_low[4],
        "main": f"({prop_main[0]:.3f},{prop_main[1]:.3f},{prop_main[2]:.3f})",
        "sub": f"({prop_sub[0]:.3f},{prop_sub[1]:.3f},{prop_sub[2]:.3f})"
    })

df = pd.DataFrame(rows)
df.to_excel(results_path, index=False)