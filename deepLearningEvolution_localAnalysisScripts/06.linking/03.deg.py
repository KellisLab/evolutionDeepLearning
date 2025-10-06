import csv


def get_genes_from_links(filepath):
    genes = set()
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 11:
                gene = fields[10].strip()
                if gene:
                    genes.add(gene)
    return genes

def get_human_chimp_degs(filepath):
    deg_genes = set()
    with open(filepath, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        # Fix BOM issue by renaming first field
        if reader.fieldnames and reader.fieldnames[0].startswith('\ufeff'):
            reader.fieldnames[0] = reader.fieldnames[0].replace('\ufeff', '')
        for row in reader:
            if row['group_1'] == 'human' and row['group_2'] == 'chimp':
                gene = row['genes'].strip()
                if gene:
                    deg_genes.add(gene)
    return deg_genes

def main():
    #brainLinks = "mergedOutput/upInHominins.merged.links_by_group.brain.links.txt"
    brainLinks = "mergedOutput/brain.upInHominin.links_by_group.brain.hg38.links.txt"
    degFile = "science.ade9516_tables_s1_to_s8/science.ade9516_table_s3.csv"
    degs = get_human_chimp_degs(degFile)
    linkedGenes = get_genes_from_links(brainLinks)
    shared = linkedGenes & degs
    if shared:
        for gene in sorted(shared):
            print(gene)


if __name__ == '__main__':
    main()