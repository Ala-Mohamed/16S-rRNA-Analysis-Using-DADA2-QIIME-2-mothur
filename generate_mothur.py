import re

input_file = "mock3.count_table"
count_out = "mock3.dada2.count_table"
shared_out = "mock3.shared"

with open(input_file, 'r') as f:
    lines = [line.strip() for line in f if line.strip()]

header = re.split(r'\s+', lines[0])
samples = header[2:]
num_samples = len(samples)

asvs = []
totals = []
counts_matrix = []

for line in lines[1:]:
    parts = re.split(r'\s+', line)
    asvs.append(parts[0])
    totals.append(parts[1])
    counts_matrix.append([int(x) for x in parts[2:]])

# Write .count_table
with open(count_out, 'w') as f:
    f.write("#Compressed Format: groupIndex,abundance.\n")
    for idx, sample in enumerate(samples, 1):
        f.write(f"#{idx},{sample}\n")
    f.write("Representative_Sequence\ttotal\t" + "\t".join(samples) + "\n")
    
    for asv, total, row in zip(asvs, totals, counts_matrix):
        group_pairs = [f"{i+1},{cnt}" for i, cnt in enumerate(row) if cnt > 0]
        f.write(f"{asv}\t{total}\t" + "\t".join(group_pairs) + "\n")

# Write .shared
with open(shared_out, 'w') as f:
    f.write("label\tGroup\tnumOtus\t" + "\t".join(asvs) + "\n")
    for s_idx, sample in enumerate(samples):
        sample_counts = [str(counts_matrix[a_idx][s_idx]) for a_idx in range(len(asvs))]
        f.write(f"0.03\t{sample}\t{len(asvs)}\t" + "\t".join(sample_counts) + "\n")

print("Files generated cleanly!")
