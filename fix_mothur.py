import re

input_path = "/home/alaa/Downloads/mock3.count_table"
count_out = "/home/alaa/Downloads/mock3.dada2.count_table"
shared_out = "/home/alaa/Downloads/mock3.shared"

with open(input_path, "r") as f:
    lines = [line.strip() for line in f if line.strip()]

# Extract header sample names
header = re.split(r"\s+", lines[0])
samples = header[2:]

asvs, totals, matrix = [], [], []

for line in lines[1:]:
    parts = re.split(r"\s+", line)
    asvs.append(parts[0])
    totals.append(parts[1])
    matrix.append([int(x) for x in parts[2:]])

# Write compressed .count_table
with open(count_out, "w") as f:
    f.write("#Compressed Format: groupIndex,abundance.\n")
    for idx, sample in enumerate(samples, 1):
        f.write(f"#{idx},{sample}\n")
    f.write("Representative_Sequence\ttotal\t" + "\t".join(samples) + "\n")
    
    for asv, total, row in zip(asvs, totals, matrix):
        pairs = [f"{i+1},{cnt}" for i, cnt in enumerate(row) if cnt > 0]
        f.write(f"{asv}\t{total}\t" + "\t".join(pairs) + "\n")

# Write .shared table
with open(shared_out, "w") as f:
    f.write("label\tGroup\tnumOtus\t" + "\t".join(asvs) + "\n")
    for s_idx, sample in enumerate(samples):
        counts = [str(matrix[a_idx][s_idx]) for a_idx in range(len(asvs))]
        f.write(f"0.03\t{sample}\t{len(asvs)}\t" + "\t".join(counts) + "\n")

print("Files generated cleanly!")
