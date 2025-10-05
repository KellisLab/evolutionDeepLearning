import os
import pyBigWig

liftedRegionsAFilepath = "/home/tnikitha/data/ancoraUpdated/liftedRegionsA"
bigWigOutput = "/home/tnikitha/data/ancoraUpdated/bigWigOutput"

def visualize(regionName, trackName, padding=1000):
	output_matrix = []
	for sampleCoordsFile in os.listdir(liftedRegionsAFilepath):
		filename_fields = sampleCoordsFile.split('.')
		sampleName = ".".join(filename_fields[1:-1])
		with open(liftedRegionsAFilepath + f"/{sampleCoordsFile}", 'r') as coordConversionFile:
			lines = coordConversionFile.readlines()
			for line in lines:
				if regionName in line:
					regionLine = line.split('\t')
					regionChrom, regionStart, regionEnd = regionLine[0], regionLine[1], regionLine[2]
					regionCenter = int(regionStart) + (int(regionEnd) - int(regionStart))//2
					regionToVis = (int(regionCenter) - padding, int(regionCenter) + padding)
					bigWigFilepath = f"/home/tnikitha/data/ancoraUpdated/bigWigOutput/{sampleName}/{sampleName}.{trackName}.bw"
					bw = pyBigWig.open(bigWigFilepath)
					accessibility = bw.values(regionChrom, regionToVis[0], regionToVis[1])
					output_matrix.append([sampleName] + accessibility)
	return output_matrix

nested_output = visualize("HAQER0710", "bipolar_neuron_from_GM23338_doxycycline_4_days", padding=2500)
with open('visOutput.txt', 'w') as file:
	for row in nested_output:
		file.write("\t".join(map(str, row)) + "\n")
