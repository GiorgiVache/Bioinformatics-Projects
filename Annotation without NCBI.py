"""The program takes a gene marks file as an input and writes a fasta file where the genes from gene marks will be
translated to amino acid sequences in fasta"""

TRANSLATION_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

import re
import textwrap


def embl_generator(genemarks_file, fasta_file, locus_tag):
    """Generates an embl file with information about its genes and the proteins the genes synthesize."""
    with open(genemarks_file, "r") as source_file:
        content = source_file.read()
        split_content = content.split("\n\n")
        if split_content[-1] == '\n':
            split_content = split_content[:-1]
        with open("annotation_file_without_NCBI.embl", "w") as embl_file:
            pattern = r"\d+\|\d+"
            for gene in split_content:
                info_tuple = ()
                if "+" in gene:
                    match = re.findall(pattern, gene)
                    new_match = match[0].replace("|", "..")
                    # print(new_match)
                    info_tuple += (f"{new_match}",)
                    # print(info_tuple[0])
                elif "-" in gene:
                    match = re.findall(pattern, gene)
                    new_match = match[0].replace("|", "..")
                    info_tuple += (f"complement({new_match})",)
                    # print(info_tuple)
                lines = gene.split('\n')
                nucleotide_sequence = ''.join(lines[1:])
                amino_acid_sequence = ''.join(
                    TRANSLATION_TABLE[nucleotide_sequence[i:i + 3]] for i in
                    range(0, len(nucleotide_sequence), 3)).replace("*", "")
                info_tuple += (amino_acid_sequence,)
                embl_content = f"""FT   gene            {info_tuple[0]}
FT                   /locus_tag="{locus_tag + "-" + str(split_content.index(gene) + 1).zfill(3)}"
FT   CDS             {info_tuple[0]}
FT                   /locus_tag="{locus_tag + "-" + str(split_content.index(gene) + 1).zfill(3)}"
FT                   /codon_start=1
FT                   /transl_table=11
FT                   /product=""
{textwrap.fill(info_tuple[1], width=79, initial_indent='FT                   /translation="', subsequent_indent=
                'FT                   ')}"
"""
                embl_file.write(embl_content)
                print(embl_content)
            with open(fasta_file, 'r') as fast_file:
                embl_file.write(fast_file.read())


x = input("Please write the 'Genemarks' file name with txt extension: ")
y = input("Please write the Fasta file name with txt extension, containing whole genome: ")
z = input("Please write the locus tag: ")
embl_generator(x, y, z)
