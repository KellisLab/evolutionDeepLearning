import sys, getopt, kipoiseq, pyfaidx
import numpy as np
import tensorflow as tf
import tensorflow_hub as hub
import pandas as pd
from kipoiseq import Interval

def usage():
    print('\nexpressionAnalysis.py is a script to run gene expression analysis using Enformer.\n')
    print('Usage:')
    print('expressionAnalysis -f <fastaFile> -g <geneTssBed> -t <trackFile> -o <outFile> [-c <category>] [-w <windows>] [-m <modelDir>]\n')
    print('Options and Arguments:')
    print('-f, --fastaFile\t\tFasta format file for gene sequence. Note that this program expects a fasta index .fai file to be found in the same directory as the fasta file.')
    print('-g, --geneTssBed\tBED file containing gene TSS locations. Expects 4 columns, with the gene name in the fourth field.')
    print('-t, --trackFile\t\tSpecify which expression tracks from which to report gene expression information.')
    print('-o, --outFile\t\tOutput filename. Will be a tabular file.')
    print('-c, --category\t\t(optional) Sample group for differential expression analysis.')
    print('-w, --windows\t\t(optional) Specify the number of TSS adjacent windows to consider for expression analysis. (default 5)')
    print('-m, --modelDir\t\t(optional) Specify the path to a directory containing the pre-trained Enformer model. (default "/net/bmc-lab4/data/kellis/users/rimangan/enformer/enformer_1")')

def main(argv):
    fastaFile = ''
    category = ''
    geneTssBed = ''
    trackFile = ''
    outFile = ''
    modelDir = '/net/bmc-lab4/data/kellis/users/rimangan/enformer/enformer_1' # default location on luria
    windows = 3 # default value for windows
    try:
        opts, args = getopt.getopt(argv, "hm:f:c:g:t:o:w:", ["modelDir", "fastaFile", "category", "geneTssBed", "trackFile", "outFile", "windows"])
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
        elif opt in ("-g", "--geneTssBed"):
            geneTssBed = arg
        elif opt in ("-t", "--trackFile"):
            trackFile = arg    
        elif opt in ("-o", "--outFile"):
            outFile = arg
        elif opt in ("-w", "--windows"):
            windows = int(arg)
    if fastaFile == '' or geneTssBed == '' or trackFile == '' or outFile == '':
        usage()
        sys.exit("Required options not found.")
    return fastaFile, category, geneTssBed, trackFile, outFile, windows, modelDir

if __name__ == "__main__":
    fastaFile, category, geneTssBed, trackFile, outFile, windows, modelDir = main(sys.argv[1:])

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
def TssToInterval(inter):
      newStart = inter.start - 57344
      newEnd = inter.end + 57344 - 1
      return kipoiseq.Interval(inter.chrom, newStart, newEnd)

# func for parsing the gene TSS bed argument into a pandas data frame. Catches some common errors.
def ParseTssFile(geneTssBed):
    headerNames = ["Chrom", "ChromStart", "ChromEnd", "Name"]
    try:
        bedData = pd.read_csv(geneTssBed, sep="\t", header=None, names=headerNames)
        if len(bedData.columns) != len(headerNames):
            raise ValueError("Invalid number of columns in the gene TSS BED file. Expected 4.")    
    except pd.errors.EmptyDataError:
        print("Empty TSS Bed File.")
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
tssBed = ParseTssFile(geneTssBed)
tracks = pd.read_csv(trackFile, sep="\t")

# now we create the output file and write the header
out = open(outFile, "w")
out.write("Gene\tSample\tTrack\tGroup\tAmplitude\n")

# now we load the Enformer model
model = hub.load(modelDir).model

# for every gene from which we want a prediction, as specified by the geneTssBed argument
for geneIndex, gene in tssBed.iterrows():
    # first, convert the current TSS into an 'Interval'
    currGeneInterval = TssToInterval(kipoiseq.Interval(gene['Chrom'], int(gene['ChromStart']), int(gene['ChromEnd'])))
    # next, we convert the Interval into a one-hot encoded sequence matrix, extracting sequence info from the fasta.
    currOneHot = one_hot_encode(fastaExtractor.extract(currGeneInterval.resize(SEQUENCE_LENGTH)))
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
        outputLine = "{}\t{}\t{}\t{}\t{}\n".format(gene['Name'], fastaFile, track['Name'], category, maxAmplitude)    
        out.write(outputLine)

# all done, time to close the open files
fastaExtractor.close()
out.close()