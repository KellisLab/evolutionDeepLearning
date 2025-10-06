import sys, getopt, kipoiseq, pyfaidx
import numpy as np
import tensorflow as tf
import tensorflow_hub as hub
import pandas as pd
from kipoiseq import Interval

def usage():
    print('\n amplitudeAnalysis.py is a script to report peak amplitudes for select genomic regions using Enformer.\n')
    print('Usage:')
    print('amplitudeAnalysis -f <fastaFile> -r <regionsBed> -t <trackFile> -o <outFile> [-c <category>] [-w <windows>] [-m <modelDir>]\n')
    print('Options and Arguments:')
    print('-f, --fastaFile\t\tFasta format file for genome sequence. Note that this program expects a fasta index .fai file to be found in the same directory as the fasta file.')
    print('-r, --regionsBed\tBED file containing peak locations. Expects 4 columns, with the region name in the fourth field.')
    print('-t, --trackFile\t\tSpecify which tracks from which to report peak amplitude information.')
    print('-o, --outFile\t\tOutput filename. Will be a tabular file.')
    print('-c, --category\t\t(optional) Sample group for differential expression/accessibility analysis.')
    print('-w, --windows\t\t(optional) Specify the number of peak adjacent windows to consider for maximum amplitude measurement. (default 3)')
    print('-m, --modelDir\t\t(optional) Specify the path to a directory containing the pre-trained Enformer model. (default "/net/bmc-lab4/data/kellis/users/rimangan/enformer/enformer_1")')

def main(argv):
    fastaFile = ''
    category = ''
    regionsBed = ''
    trackFile = ''
    outFile = ''
    modelDir = '/net/bmc-lab4/data/kellis/users/rimangan/enformer/enformer_1' # default location on luria
    windows = 3 # default value for windows
    try:
        opts, args = getopt.getopt(argv, "hm:f:c:r:t:o:w:", ["modelDir", "fastaFile", "category", "regionsBed", "trackFile", "outFile", "windows"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-m", "--modelDir"):
            modelDir = arg    
        elif opt in ("-f", "--fastaFile"):
            fastaFile = arg
        elif opt in ("-c", "--category"):
            category = arg
        elif opt in ("-r", "--regionsBed"):
            regionsBed = arg
        elif opt in ("-t", "--trackFile"):
            trackFile = arg    
        elif opt in ("-o", "--outFile"):
            outFile = arg
        elif opt in ("-w", "--windows"):
            windows = int(arg)
    if fastaFile == '' or regionsBed == '' or trackFile == '' or outFile == '':
        usage()
        sys.exit("Required options not found.")
    return fastaFile, category, regionsBed, trackFile, outFile, windows, modelDir

if __name__ == "__main__":
    fastaFile, category, regionsBed, trackFile, outFile, windows, modelDir = main(sys.argv[1:])

#Now I've copied code from the Enformer repository
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

# FastaStringExtractor is a class that enables fasta file I/O and allows us 
  # to extract subsequences from a fasta file as strings based on a target 
  # interval.
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

# one_hot_encode converts a string into a one-hot matrix representation of a DNA sequence.
def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

# resize intervals to 114_688 bp, the Enformer prediction window size
def RegionToInterval(inter):
      midpoint = int((inter.start + inter.end) / 2)
      newStart = midpoint - 57344
      newEnd = midpoint + 57344 - 1
      return kipoiseq.Interval(inter.chrom, newStart, newEnd)

# func for parsing the regionsBed argument into a pandas data frame. Catches some common errors.
def ParseRegionFile(regionsBed):
    headerNames = ["Chrom", "ChromStart", "ChromEnd", "Name"]
    try:
        bedData = pd.read_csv(regionsBed, sep="\t", header=None, names=headerNames)
        if len(bedData.columns) != len(headerNames):
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


# on to the main event, first we parse the input files
# note that the fasta file is converted into a FastaStringExtractor object
fastaExtractor = FastaStringExtractor(fastaFile)
regions = ParseRegionFile(regionsBed)
tracks = pd.read_csv(trackFile, sep="\t")

# now we create the output file and write the header
out = open(outFile, "w")
out.write("Region\tSample\tTrack\tGroup\tAmplitude\n")

# now we load the Enformer model
model = hub.load(modelDir).model

# for every region from which we want a prediction, as specified by the regionsBed argument
for regionIndex, region in regions.iterrows():
    # first, convert the current region into an 'Interval'
    currInterval = RegionToInterval(kipoiseq.Interval(region['Chrom'], int(region['ChromStart']), int(region['ChromEnd'])))
    # next, we convert the Interval into a one-hot encoded sequence matrix, extracting sequence info from the fasta.
    currOneHot = one_hot_encode(fastaExtractor.extract(currInterval.resize(SEQUENCE_LENGTH)))
    currPredictions = model.predict_on_batch(currOneHot[np.newaxis])['human'][0]
    # for each track we'd like to report, as specified by the trackFile.
    for trackIndex, track in tracks.iterrows():
        currWindow = int((len(currPredictions[:, int(track['Index'])]) / 2) - windows)
        rightWindow = currWindow + 2*windows
        maxAmplitude = 0.0
        while currWindow <= rightWindow:
            if currPredictions[:, int(track['Index'])][currWindow] > maxAmplitude:
                maxAmplitude = currPredictions[:, int(track['Index'])][currWindow]
            currWindow = currWindow + 1
        outputLine = "{}\t{}\t{}\t{}\t{}\n".format(region['Name'], fastaFile, track['Name'], category, maxAmplitude)    
        out.write(outputLine)

# all done, time to close the open files
fastaExtractor.close()
out.close()