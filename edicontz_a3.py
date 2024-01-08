"""
Author : Edicon Chan
Date   : 2023-12-13
Assignment3 : Binding Motif Search
Linting: Pylint, flake8, Mypy
"""

from collections import namedtuple
import re
import os
import argparse
import random
# pip install biopython
from Bio.SeqIO.FastaIO import SimpleFastaParser


def get_args():
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description=' Transcription Factor Binding Site Detection ',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta',
                        metavar='fasta_dir',
                        type=str,
                        help='Input the directory containing the fasta files')

    parser.add_argument('gff3',
                        metavar='gff3_dir',
                        type=str,
                        help='Input the directory containing the gff3 files')

    parser.add_argument('pro',
                        metavar='promoters_file',
                        type=str,
                        help='Input file containing motifs of interest')

    parser.add_argument('gene',
                        metavar='genes_file',
                        type=str,
                        help='Input file containing list of selected genes')

    parser.add_argument('-o',
                        '--out',
                        help='Name of the output file',
                        metavar='file',
                        type=str,
                        default='output.txt')

    parser.add_argument('-n',
                        '--num',
                        help='Number of random samples',
                        metavar='int',
                        type=int,
                        default=5)

    args = parser.parse_args()
    Args = namedtuple("Args", "fasta gff3 pro gene out sample")

    return Args(args.fasta, args.gff3, args.pro, args.gene, args.out, args.num)


def reverse_complement(dna_string: str) -> str:
    """ Return the reverse complementary sequence of a DNA strand
    containing both uppercase and lowercase bases """

    comp_dna = ""
    comp_dict = {"A": "T",
                 "T": "A",
                 "C": "G",
                 "G": "C",
                 "a": "t",
                 "t": "a",
                 "c": "g",
                 "g": "c"}

    for base in dna_string:
        if base not in comp_dict:
            comp_dna += base
        else:
            comp_dna += comp_dict[base]

    return comp_dna[::-1]


def dir_handler(fasta_dir: str, gff3_dir: str) -> tuple:
    """ Search directories for all fasta and gff3 files """

    fasta_list = []
    gff3_list = []

    for file in os.listdir(fasta_dir):
        if file.endswith(".fa"):
            fasta_list.append(file)
    for file in os.listdir(gff3_dir):
        if file.endswith(".gff3"):
            gff3_list.append(file)

    return fasta_list, gff3_list


def file_handler(file: str) -> list:
    """ Returns a list containing all of the lines in a file"""

    with open(file, 'r', encoding='utf8') as file_lines:
        file_contents = []
        for line in file_lines.readlines():
            file_contents.append(line.rstrip())

    return file_contents


def gff3_parser(gff3_list: list, gff3_dir: str) -> list:
    """ Return a list containing object-like data structures
    of genes extracted from the provided gff3 files """

    print("GFF3 files parsing...")
    gff3_genes = []
    dir_name = gff3_dir + "/"

    # For each gff3 file provided to the command line
    for name in gff3_list:
        with open(dir_name + name, 'r', encoding='utf8') as file:
            # Search and keep all lines with the word gene
            for line in file.readlines():
                if line[0] != '#' and line.find("gene") != -1:
                    gff3_genes.append(line)

    # Initialize list to contain the namedtuple data structures for the genes
    gene_list = []
    GFF3 = namedtuple('GFF3', 'name chrom start end strand')

    # For each gff3 gene line, split and assign attributes for name,
    # chromosome num, start and end indices, and positive/negative strand
    for line in gff3_genes:
        gene_info = line.split("\t")
        gene_name = gene_info[8].split(":")[1].split(";")[0]
        gff3_line = GFF3(gene_name, int(gene_info[0]),
                         int(gene_info[3]), int(gene_info[4]),
                         gene_info[6])
        gene_list.append(gff3_line)

    print("GFF3 files finished parsing!")
    return gene_list


def fasta_parser(fasta_list: list, fasta_dir: str) -> dict:
    """Return a dict containing the sequence for each chromosome"""

    print("FASTA files parsing...")
    total_seq = {}
    dir_name = fasta_dir + "/"

    for name in fasta_list:
        with open(dir_name + name, 'r', encoding='utf8') as file:
            for title, seq in SimpleFastaParser(file):
                # Extract and save chromosome number
                chrom = title.split(":")[3]
                total_seq[chrom] = seq

    print("FASTA files finished parsing!")
    return total_seq


def promoters(selected_genes: list, total_seq: dict) -> list:
    """Returns a list containing the 500 bp promoter regions
    upstream of the transcription start site for every selected
    gene provided. Regions may be smaller if there is a presence
    of gaps symbolized by 100 Ns"""

    promoter_regions = []

    for gff3 in selected_genes:
        # If the strand is positive, search upstream
        if gff3.strand == "+":
            start = gff3.start
            region = total_seq[str(gff3.chrom)][start-501:start-1]
        # If the strand is negative, search downstream of the 5'-3'
        elif gff3.strand == "-":
            start = gff3.end
            region = total_seq[str(gff3.chrom)][start:start+500]
            region = reverse_complement(region)
        # If there are gaps in the promoter, take the region between
        # the gene and the appearance of the closest gap
        if region.find("N"*100) != -1:
            rev_region = region[::-1]
            new_start = rev_region.find("N"*100)
            rev_region = rev_region[:new_start]
            region = rev_region[::-1]

        promoter_regions.append(region)

    return promoter_regions


def get_motif_count(motif_list: list, promoter_list: list) -> dict:
    """ Returns a dictionary containing the motifs and the
    number of counts detected in the promoters provided"""

    motif_count = {}

    for motif in motif_list:
        if motif not in motif_count:
            motif_count[motif] = 0
        for promoter_seq in promoter_list:
            # ?= overlooking for counting overlapping motifs
            motif_count[motif] += len(re.findall(rf"(?={motif})",
                                                 promoter_seq,
                                                 re.IGNORECASE))

    return motif_count


def random_promoter_sampling(gene_list: list, promoter_txt: list,
                             total_seq: dict, num_of_genes: int) -> dict:
    """ Randomly sample genes in the provided gff3 file and
    return a motif count for those promoters"""

    sampled_promoters = promoters(random.sample(gene_list, num_of_genes),
                                  total_seq)

    return get_motif_count(promoter_txt, sampled_promoters)


def motif_output(output_path: str, motif_count: dict,
                 gene_list: list, promoter_txt: list,
                 total_seq: dict, num_of_genes: int,
                 sample_number: int) -> None:
    """ Output a file containing the motifs and the number of
    counts found in the selected promoters. Subsequent random
    sampling of promoters and their motif counts are included
    Default random sampling set to 5"""

    print("Writing motif counts out...")

    with open(output_path, 'w', encoding='utf8') as output:
        count = 0
        for key in motif_count:
            count += 1
            row_output = key + "\t" + str(motif_count[key])
            loop_index = 0
            while loop_index < sample_number:
                sample = random_promoter_sampling(gene_list, promoter_txt,
                                                  total_seq, num_of_genes)
                row_output += "\t" + str(sample[key])
                loop_index += 1

            # Output percentage of completed motif counts to console
            percentage = int(count/len(motif_count) * 100)
            print("Completed: " + str(percentage) + "% " + row_output)
            # Write output to provided path
            output.write(row_output + "\n")

    print("Finished writing to output file. " +
          f"Output file can be found at {output_path}")


def main():
    """Main function"""

    args = get_args()

    # Loading the input files
    file_lists = dir_handler(args.fasta, args.gff3)
    promoter_txt = file_handler(args.pro)
    gene_names = file_handler(args.gene)
    output_path = args.out
    sample_number = args.sample

    # Parsing the sequences and gene information
    total_seq = fasta_parser(file_lists[0], args.fasta)
    gene_list = gff3_parser(file_lists[1], args.gff3)

    # Selecting the genes given by the input file
    selected_genes = []
    for gene in gene_list:
        if gene.name in gene_names:
            selected_genes.append(gene)

    # Obtaining the motif counts for the selected promoters
    selected_promoters = promoters(selected_genes, total_seq)
    motif_count = get_motif_count(promoter_txt, selected_promoters)

    # Outputting the motif counts, as well as the counts
    # for a random sampling of other promoters
    num_of_genes = len(selected_genes)
    motif_output(output_path, motif_count,
                 gene_list, promoter_txt,
                 total_seq, num_of_genes, sample_number)


if __name__ == "__main__":
    main()
