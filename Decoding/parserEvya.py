import pandas as pd
import edlib
from tqdm import tqdm

# Open output file
file_as = open("evyaBinning2Lev_T100M7_sup.txt", "w")

# Input files
filename = "matched_results_2.csv"
filename_des = "design.csv"

# Read data
seq_df = pd.read_csv(filename, low_memory=False)
des_df = pd.read_csv(filename_des)
var_list = des_df['sequence'].tolist()

# Group by variant_id
grouped = seq_df.groupby(['variant_id'], sort=False)

# Iterate over clusters with progress bar
for name_of_the_group, gr in tqdm(grouped, desc="Processing clusters"):
    if name_of_the_group == -1:
        continue

    list_of_sequence = gr['sequence'].to_numpy()
    list_of_counts = gr['count'].to_numpy()
    reference = var_list[name_of_the_group]

    reads = []

    for item, count in zip(list_of_sequence, list_of_counts):
        count_of_read = int(count)
        for _ in range(count_of_read):
            if edlib.align(item, reference)['editDistance'] < 200:
                reads.append(item)

    if reads:
        # Sort reads: first those of length 140, preserving original order within group
        reads.sort(key=lambda x: 0 if len(x) == 140 else 1)

        file_as.write(reference + "\n")
        file_as.write("*****************************\n")
        for read in reads:
            file_as.write(read + "\n")
        file_as.write("\n\n")

file_as.close()
