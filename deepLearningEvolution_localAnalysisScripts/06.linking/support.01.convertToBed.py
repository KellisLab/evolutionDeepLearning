import csv
import gzip
import argparse
import re
from pathlib import Path

def parse_gtf_for_mapping(gtf_path):
    """Parses GTF, creates dict of gene/id to gene symbol."""
    mapping = {}
    with gzip.open(gtf_path, 'rt') as gtf_file:
        for line in gtf_file:
            if line.startswith('#') or '\tgene\t' not in line:
                continue
            fields = line.strip().split('\t')
            attr = fields[8]
            gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attr)
            if gene_id_match and gene_name_match:
                gene_id = gene_id_match.group(1).split('.')[0]
                gene_name = gene_name_match.group(1)
                mapping[gene_id] = gene_name
    return mapping

def convert_tsv_to_bed(tsv_gz_path, mapping, output_path):
    """Reads a gzipped TSV file and writes a BED file with gene symbols."""
    with gzip.open(tsv_gz_path, 'rt') as infile, open(output_path, 'w') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            chrom = row['chr']
            start = row['start']
            end = row['end']
            gene_id = row['gene'].split('.')[0]
            gene_symbol = mapping.get(gene_id, gene_id)
            outfile.write(f"{chrom}\t{start}\t{end}\t{gene_symbol}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert a directory of enhancer–gene .tsv.gz files to BED with gene symbols.")
    parser.add_argument('--input_dir', required=True, type=Path, help='Directory containing .tsv.gz enhancer–gene link files')
    parser.add_argument('--gtf', required=True, type=Path, help='Path to gzipped GTF file')
    parser.add_argument('--output_dir', required=True, type=Path, help='Directory to save BED output files')
    args = parser.parse_args()

    input_dir = args.input_dir
    gtf_path = args.gtf
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    mapping = parse_gtf_for_mapping(gtf_path)
    tsv_files = sorted(input_dir.glob("*.tsv.gz"))

    for tsv_path in tsv_files:
        out_name = tsv_path.name.replace(".tsv.gz", ".bed")
        output_path = output_dir / out_name
        convert_tsv_to_bed(tsv_path, mapping, output_path)

if __name__ == '__main__':
    main()
