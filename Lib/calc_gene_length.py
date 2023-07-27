from Bio import SeqIO
import argparse

def extract_gene_info(fasta_files):
    gene_info = []
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.description
            gene_length = len(record.seq)
            aro_term = header.split("|")[0].replace(">", "")
            gene_name = header.split("|")[2].split(":")[1]
            gene_info.append([aro_term, gene_name, gene_length])
    return gene_info

def write_gene_info_to_file(gene_info, output_file):
    with open(output_file, "w") as f:
        f.write("ARO_Term\tGene_Name\tGene_Length\n")
        for info in gene_info:
            f.write(f"{info[0]}\t{info[1]}\t{info[2]}\n")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Take CARD database fasta files and calculate gene lengths for GCPM calculations')
    parser.add_argument('-i', '--fasta', nargs='+', help='Path for input fasta file(s)', required=True)
    parser.add_argument('-o', '--outfile', help='Path for output file', required=True)
    args = parser.parse_args()

    gene_info = extract_gene_info(args.fasta)
    write_gene_info_to_file(gene_info, args.outfile)
