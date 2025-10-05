import os
import pyBigWig

# File paths
liftedRegionsAFilePath = "/net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/liftedRegionsA"
liftedRegionsBFilePath = "/net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/liftedRegionsB"
bigWigOutput = "/net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/bigWigOutput"

def visualize(regionName, biosample, padding=10000):
    output_matrix = []
    def process_directory(directoryPath):
        for sampleCoordsFile in os.listdir(directoryPath):
            print(sampleCoordsFile)
            filename_fields = sampleCoordsFile.split('.')
            sampleName = ".".join(filename_fields[1:-1])
            with open(os.path.join(directoryPath, sampleCoordsFile), 'r') as coordConversionFile:
                lines = coordConversionFile.readlines()
                for line in lines:
                    if regionName in line:
                        regionLine = line.split('\t')
                        regionChrom, regionStart, regionEnd = regionLine[0], regionLine[1], regionLine[2]
                        regionCenter = int(regionStart) + (int(regionEnd) - int(regionStart)) // 2
                        regionToVis = (int(regionCenter) - padding, int(regionCenter) + padding)
                        bigWigFilePath = f"{bigWigOutput}/{sampleName}/{sampleName}.{biosample}.bw"
                        if os.path.exists(bigWigFilePath):  # Ensure BigWig file exists
                            with pyBigWig.open(bigWigFilePath) as bw:
                                accessibility = bw.values(regionChrom, regionToVis[0], regionToVis[1])
                                output_matrix.append([sampleName] + accessibility)
    process_directory(liftedRegionsAFilePath)
    process_directory(liftedRegionsBFilePath)
    return output_matrix


nested_output = visualize("chr20.45095781.45096191", "brain_male_embryo_105_days", padding = 2500)
print(nested_output)
with open('visOutput.txt', 'w') as file:
	for row in nested_output:
		file.write("\t".join(map(str, row)) + "\n")

