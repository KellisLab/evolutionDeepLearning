import sys, getopt, os, kipoiseq, pyfaidx
import numpy as np
import tensorflow as tf
import tensorflow_hub as hub
import pandas as pd
from kipoiseq import Interval
import copy

### Input file paths ###
def usage():
    print('\nsaturatedMutagenesisGenerator outputs a txt file for saturated mutagenesis of on cell type specified by trackIndex. ')
    print('Usage:\n')
    print('saturatedMutagenesisGenerator.py -e <enformerPath> -f <fastaFilePath> -b <bedFilePath> -o <outputPath>')
    print('-e, --enformerPath\t\tSpecify a path to an Enformer model (model_code + variables).')
    print('-f, --fastaFilePath\t\tFasta format file (end in .fa; ex. hg38.fa for human genome) for genome of interest.')
    print('-b, --bedFilePath\t\tBED format file (end in .bed) listing genes regions to analyze.')
    print('-o, --outputPath\t\tSpecity full path of the file to write the tab delimited result table.')
    print('-t, --trackIndex\t\t(Optional) Track index for cell type, default=485 (bipolar neuron).')
    print('-w, --numWindows\t\t(Optional) number of windows left & right of a region specfied to make the amplitude measurement, default=3.')

def main(argv):
    enformerPath = ''
    fastaFilePath = ''
    bedFilePath = ''
    outputPath = ''
    trackIndex = 485 # indicates of cell type for Enformer: https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt
    numWindows = 3
    try:
      opts, args = getopt.getopt(argv, "he:f:b:o:t:w:", ["enformerPath", "fastaFilePath", "bedFilePath", "outputPath", "trackIndex", "numWindows"]) # : means look for long names
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:   # entered as opt arg pairs in command
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-e", "--enformerPath"):
            enformerPath = arg
        elif opt in ("-f", "--fastaFilePath"):
            fastaFilePath = arg
        elif opt in ("-b", "--bedFilePath"):
            bedFilePath = arg
        elif opt in ("-o", "--outputPath"):
            outputPath = arg
        elif opt in ("-t", "--trackIndex"):
            trackIndex = int(arg)
        elif opt in ("-w", "--numWindows"):
            numWindows = arg
    if enformerPath == '' or fastaFilePath == '' or bedFilePath == '' or outputPath == '': # ensure required options entered
        usage()
        sys.exit("Required options not found.")
    return enformerPath, fastaFilePath, bedFilePath, outputPath, trackIndex, numWindows

if __name__ == "__main__":
    enformerPath, fastaFilePath, bedFilePath, outputPath, trackIndex, numWindows = main(sys.argv[1:])

### Code from the Enformer repository ###
SEQUENCE_LENGTH = 393216

class Enformer:

  def __init__(self, tfhub_path):
    self._model = hub.KerasLayer(tfhub_path)
    #self._model = hub.load(tfhub_url).model# hub is a library for transfer learning, handles loading models from urls

  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}

  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]

    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:# records differentiation operations so the tape can be "watched" for gradient-based interpretability
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass

    input_grad = tape.gradient(prediction, input_sequence) * input_sequence
    input_grad = tf.squeeze(input_grad, axis=0)# squeeze removes dimensions of size one from a tensor.
    return tf.reduce_sum(input_grad, axis=-1)

### FastaStringExtractor - enables fasta file I/O to extract subsequences from a fasta file as strings based on a target interval ###
class FastaStringExtractor:

    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()

def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

# resize intervals to 114_688 bp, the Enformer prediction window size
def ResizeRegion(inter):
      midpoint = int((inter.start + inter.end) / 2)
      newStart = midpoint - 57344
      newEnd = midpoint + 57344 - 1
      return kipoiseq.Interval(inter.chrom, newStart, newEnd)

# func for parsing the regionsBed argument into a pandas data frame. Catches some common errors.
def ParseRegionFile(regionsBed):
    headerNames = ["chrom", "start", "end", "name"]
    try:
        bedData = pd.read_csv(regionsBed, sep="\t", header=None, names=headerNames)
        if len(bedData.columns) != 4:
            raise ValueError("Invalid number of columns in the regions TSS BED file. Expected 4.")
    except pd.errors.EmptyDataError:
        print("Empty Bed File.")
        return None
    except pd.errors.ParserError:
        print("BED parsing error.")
        return None
    except ValueError as ve:
        print(ve)
        return None
    return bedData

# base is an index (A = 0; C = 1; G = 2; T = 3)
def mutateBase(matrix, position, base):
   matrix[position] = 0
   matrix[position, base] = 1
   return matrix

def print_one_hot_region(matrix, region_start, region_end, print_width):
   region_start = max(0, region_start)
   region_end = min(len(matrix), region_end)
   for i in range(region_start, region_end, print_width):
      subset_matrix = matrix[i:i+print_width]
      for base_vector in subset_matrix:
         base = ['A', 'C', 'G', 'T'][np.argmax(base_vector)]
         print(base, end='')
      print()

def print_middle_region_of_one_hot(matrix, region_size, print_width=50):
   region_start = max(0, len(matrix) // 2 - region_size // 2)
   region_end = min(len(matrix), region_start + region_size)
   print_one_hot_region(matrix, region_start, region_end, print_width)

def maxAmplitudeFromPredictions(predictions, trackIndex, windows):
    currWindow = int((len(predictions[:, int(trackIndex)]) / 2) - windows)
    rightWindow = currWindow + 2*windows
    maxAmplitude = 0.0
    while currWindow <= rightWindow:
        if predictions[:, int(trackIndex)][currWindow] > maxAmplitude:
           maxAmplitude = predictions[:, int(trackIndex)][currWindow]
        currWindow += 1
    return maxAmplitude

### Main code ###
# parse input files
# Code to specify one region; now read BED file instead - region = kipoiseq.Interval("chr20", 63102094, 63102278) # kipoiseq.Interval: Object from package conceptually similar to a BED storing genomic regions; defined by which chromosome, start, stop point. Currently have coordinates for HAR1 (famous fast-evolving region in the human genome)
fastaExtractor = FastaStringExtractor(fastaFilePath) # converted into Obj that makes easy to pull specific subsequence out of a fasta file
model = hub.load(enformerPath).model # Loads the actual enformer model (https://www.kaggle.com/models/deepmind/enformer/frameworks/tensorFlow2/variations/enformer/versions/1?tfhub-redirect=true)
regions = ParseRegionFile(bedFilePath)
print(regions)

# clear previous output file, if exist, and append values for each region
if os.path.exists(outputPath):
    os.remove(outputPath)


with open(outputPath, "w") as file:
    file.write("Position\tMutation\tAmplitude\tName\n")
    # for each region, extract fasta & encode matrix
    for index, region in regions.iterrows():
        print(f"chromosome: {region['chrom']}", f"start: {region['start']}")

        interval = ResizeRegion(kipoiseq.Interval(region['chrom'], int(region['start']), int(region['end']))) # resizes input region to size of the set prediction window (200kb); region of interest will be at the midpoint of a much larger sequence window
        refOneHot = one_hot_encode(fastaExtractor.extract(interval.resize(SEQUENCE_LENGTH))) # extracts 200kb seq from fasta file & converts into One Hot Encode matrix
        # print(refOneHot)
        # print_middle_region_of_one_hot(refOneHot, region.end - region.start, 50)
        # print()

        # calls tensorflow function to get prediction back from sequence
        refPredictions = model.predict_on_batch(refOneHot[np.newaxis])['human'][0] # Enformer has human & mouse predictions, selected human here
        """
        Predictions will be a vector of values, where the index corresponds to position and the value corresponds to the amplitude of the signal.
        we have a 200kb input sequence, but the prediction vector is much shorter than that because it is not at single base-pair resolution.
        Enformer reports an amplitude for 128bp windows. Enhancers and promoters are often much larger than this (500-1000bp) so this isn't an issue.

        To measure the 'accessibility' of an enhancer given this, record the amplitude at the center of prediction window (region of interest) ± a few windows
        to ensure we get highest value (spread controlled by 'numWindows')
        """

        refAmplitude = maxAmplitudeFromPredictions(refPredictions, trackIndex,  numWindows) # ref for reference as 'wild-type' sequence
        # print("Reference\t{}\n\n".format(refAmplitude))

        region_length = int(region['end']) - int(region['start'])
        region_start_in_one_hot = max(0, len(refOneHot) // 2 - region_length // 2) # this finds the position of the start of our region in the prediction interval
        region_end = region_start_in_one_hot + region_length
        # range through all positions in sequence to test every possible variant
        for currPos in range(region_start_in_one_hot, region_end, 1):
            # for every possible mutation (standard alphabetic ordering 0-3, where 0:A, 1:C, 2:G, 3:T)
            for currMut in range(0, 4, 1):
                currOneHot = copy.deepcopy(refOneHot) # make copy of original sequence to mutate
                mutateBase(currOneHot, currPos, currMut)
                altPredictions = model.predict_on_batch(currOneHot[np.newaxis])['human'][0]
                #print_middle_region_of_one_hot(currOneHot, region.end - region.start, 50)
                altAmplitude = maxAmplitudeFromPredictions(altPredictions, trackIndex, numWindows)
                file.write("{}\t{}\t{}\t{}\n".format(currPos - region_start_in_one_hot+int(region['start']), currMut,  altAmplitude - refAmplitude, region['name']))

