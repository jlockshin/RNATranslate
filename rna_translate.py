"""
Jane Lockshin
CSCI 598A - Bioinformatics
Homework 1 Problem 1
September 20, 2018
Colorado School of Mines

"""

rna_codons = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y',  # Stop codons are here
    'UGU': 'C', 'UGC': 'C', 'UGG': 'W',  # and here

    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'N', 'AAG': 'N',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',

    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def rna_translation(input_file):
    """Returns encoded amino acids from RNA bases """

    # Open the input file (new lines separate different RNA sequences)
    with open(input_file, "r"):
        sequences = [line.rstrip('\n') for line in open(input_file)]

    # Translation: RNA to Proteins
    translations = []
    for sequence in sequences:

        # If the length of the sequence isn't divisible by 3, ignore remaining bases
        if len(sequence) % 3 == 1:
            sequence = sequence[:- 1]
        elif len(sequence) % 3 == 2:
            sequence = sequence[:- 2]

        # Split up the RNA sequence into the appropriate codons
        bases = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        amino_acids = ''

        # Loop through the codons
        for i in bases:
            # Stop the translation if stop codons are detected
            if i == 'UAA' or i == 'UAG' or i == 'UGA':
                break
            else:
                # Check that the sequence contains valid codons
                if i not in rna_codons.keys():
                    print("Value not found in RNA codon table: " + i)
                else:
                    amino_acids += rna_codons[i]

        translations += [amino_acids]

    return translations


print(rna_translation('rna_input.txt'))