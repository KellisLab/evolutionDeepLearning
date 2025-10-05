import os
import pandas as pd

output_dir = "accessibilityMatrices"
os.makedirs(output_dir, exist_ok=True)

biosampleFiles = set()
for biosampleFile in os.listdir("bigWigAverageOverBed/altai.Neanderthal.A"):
	biosampleFiles.add(biosampleFile.replace(".txt", ""))

for biosample in biosampleFiles:
	print(biosample)
	biosample_matrix = {}
	for haplotype in os.listdir("bigWigAverageOverBed"):
		with open(os.path.join("bigWigAverageOverBed", haplotype, f"{biosample}.txt")) as f:
			for line in f:
				fields = line.strip().split("\t")
				peak_id = fields[0]
				accessibility = float(fields[4])
				if peak_id not in biosample_matrix:
					biosample_matrix[peak_id] = {}
				biosample_matrix[peak_id][haplotype] = accessibility
	df = pd.DataFrame.from_dict(biosample_matrix, orient="index")
	df = df.sort_index(axis=1)
	df.to_csv(os.path.join(output_dir, f"{biosample}.txt"), sep="\t", header=True, index=True)

print("DONE")
