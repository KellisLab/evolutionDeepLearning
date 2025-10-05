#!/bin/csh -ef

import os 
import pandas as pd 

folder_path = "regionQuantificationOutput" 

data = []

for filename in os.listdir(folder_path):
	sample_part, remainder = filename.split('.A', 1)
	sample = sample_part + '.A'
	
	track = remainder.split('.')[1]
	
	filepath = os.path.join(folder_path, filename)
	with open(filepath, 'r') as f:
		lines = f.readlines()
		for line in lines:
			fields = line.split()
			region = fields[0]
			proportion_covered = float(fields[2])/float(fields[1])
			mean0_score = float(fields[4])
			if proportion_covered > 0.0:
				data.append({'sample': sample, 'track': track, 'region': region, 'mean0_score': mean0_score})
df = pd.DataFrame(data)

output_file = 'regionQuantification.A.txt'
df.to_csv(output_file, sep='\t', index=False)

print(f"Data saved to {output_file}")
