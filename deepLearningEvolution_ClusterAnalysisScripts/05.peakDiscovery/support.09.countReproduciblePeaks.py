import sys
from collections import defaultdict

def read_bed(file_path):
    peaks = []
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split()
            peaks.append(tuple(parts[:3]))  # Only store chrom, start, end
    return peaks

def count_overlaps(tmp_overlap_file):
    overlap_counts = defaultdict(int)
    with open(tmp_overlap_file) as f:
        for line in f:
            parts = line.strip().split()
            peak = tuple(parts[:3])  # Only store chrom, start, end
            overlap_counts[peak] += 1
    return overlap_counts

def filter_reproducible_peaks(overlap_counts, union_peaks, threshold=6):
    reproducible_peaks = []
    for peak in union_peaks:
        if overlap_counts[peak] >= threshold:
            reproducible_peaks.append(peak)
    return reproducible_peaks

def write_bed(file_path, peaks):
    with open(file_path, 'w') as f:
        for chrom, start, end in peaks:
            f.write(f"{chrom}\t{start}\t{end}\n")

def main(tmp_overlap_file, union_file, output_file):
    union_peaks = read_bed(union_file)
    overlap_counts = count_overlaps(tmp_overlap_file)
    reproducible_peaks = filter_reproducible_peaks(overlap_counts, union_peaks, threshold=3)
    write_bed(output_file, reproducible_peaks)

if __name__ == "__main__":
    tmp_overlap_file = sys.argv[1]
    union_file = sys.argv[2]
    output_file = sys.argv[3]
    main(tmp_overlap_file, union_file, output_file)

