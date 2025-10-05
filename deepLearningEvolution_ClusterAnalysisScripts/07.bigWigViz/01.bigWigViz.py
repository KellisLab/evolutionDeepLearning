import os
import pyBigWig
import sys
import csv

# Paths
assemblyCoordsDir = "/home/rimangan/data/enformer/macs3/finalReproduciblePeaks/assemblyCoords"
bigWigOutputDir = "/net/bmc-lab4/data/kellis/users/tnikitha/ancoraUpdated/bigWigOutput"

def visualize(regionName, biosample, padding=10000):
    output_matrix = []
    for haplotypeDir in os.listdir(assemblyCoordsDir):
        haplotypePath = os.path.join(assemblyCoordsDir, haplotypeDir)
        if os.path.isdir(haplotypePath):
            bedFilePath = os.path.join(haplotypePath, f"{biosample}.bed")
            print(bedFilePath)
            if os.path.exists(bedFilePath):
                print("Good so far")
                with open(bedFilePath, 'r') as bedFile:
                    for line in bedFile:
                        fields = line.strip().split('\t')
                        if len(fields) < 4:
                            continue  # Skip lines without a name field
                        chrom, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
                        # Check if the regionName matches the name field
                        if regionName in name:
                            regionCenter = start + (end - start) // 2
                            regionToVis = (regionCenter - padding, regionCenter + padding)
                            # Fetch accessibility values from the BigWig file
                            bigWigFilePath = os.path.join(bigWigOutputDir, haplotypeDir, f"{haplotypeDir}.{biosample}.bw")
                            if os.path.exists(bigWigFilePath):
                                with pyBigWig.open(bigWigFilePath) as bw:
                                    accessibility = bw.values(chrom, regionToVis[0], regionToVis[1])
                                    output_matrix.append([haplotypeDir, biosample] + list(accessibility))
    return output_matrix

def main(tsv_file):
    os.makedirs("output", exist_ok=True)  # ensure output/ exists
    print("started")
    with open(tsv_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            biosample, region, padding_str = row[:3]
            padding = int(padding_str)
            output_matrix = visualize(region, biosample, padding)
            
            output_file = os.path.join("output", f"{biosample}.{region}.txt")
            with open(output_file, 'w') as f:
                for row in output_matrix:
                    formatted_row = [row[0]] + [f"{x:.4f}" for x in row[2:]]
                    f.write("\t".join(map(str, formatted_row)) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 01.bigWigViz.py <input.tsv>")
        print("Expected TSV format: biosample<TAB>region<TAB>padding")
        sys.exit(1)
    main(sys.argv[1])
