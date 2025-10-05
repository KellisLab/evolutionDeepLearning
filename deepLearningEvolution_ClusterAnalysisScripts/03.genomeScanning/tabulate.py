#!/bin/csh -ef

import os 
import pandas as pd 

parent_folder = "regionQuantification"

data = []

for track_dir in os.listdir(parent_folder):
	track_path = os.path.join(parent_folder, track_dir)
	track = track_dir
	for filename in os.listdir(track_path):
		sample_part, remainder = filename.split('.A', 1)
		sample = sample_part + '.A'
		filepath = os.path.join(track_path, filename)
		with open(filepath, 'r') as f:
			lines = f.readlines()
			for line in lines:
				fields = line.split()
				region = fields[0]
				mean0_scpre = fields[4]
				data.append({'sample': sample, 'track': track, 'region': region, 'mean0_score': mean0_score})

df = pd.DataFrame(data)

output_file = 'regionQuantification.txt'
df.to_csv(output_file, sep='\t', index=False)

print(f"Data saved to {output_file}")

