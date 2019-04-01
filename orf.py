"""
Jane Lockshin
CSCI 598A - Bioinformatics
Homework 1 Problem 2
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

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def reverse_complement(sequence):
    """Returns the reverse-complement of a DNA sequence"""

    # Retrieve the complement for each base, then reverse the list of bases
    bases = [complements.get(base) for base in sequence[0]][::-1]
    reverse = ''

    # Convert the reversed list of bases to a string
    for i in bases:
        reverse += i

    return reverse


def start_codons(sequence):
    """Returns indices of start codons"""

    # Find the index of first start codons
    start = sequence.find('AUG')
    starts = []

    # Keep adding indices to list of starts as they are found
    while start >= 0:
        starts.append(start)
        start = sequence.find('AUG', start + 1)

    return starts


def find_frame_length(sequence):
    """Returns indices of stop codons"""

    # Find length of frame (until a stop codon is detected)
    frame_len = 0
    while frame_len < len(sequence) - 2 and sequence[frame_len:frame_len + 3] not in ('UAG', 'UGA', 'UAA'):
        frame_len += 3

    # Return frame length if a stop codon is found
    if frame_len < len(sequence) - 2:
        return frame_len
    else:
        return -1


def find_orfs(sequence):
    """Return indices of open reading frames of RNA sequences"""

    # Find indices of start codons
    starts = start_codons(sequence)

    # Find stop codon following the start
    frames = []
    for start_position in starts:

        # Find the length of the frame
        frame_length = find_frame_length(sequence[start_position:])

        # Find the position of the stop codon
        if frame_length != -1:
            stop_position = start_position + frame_length

            # Add start/stop positions to frames
            frames.append((start_position, stop_position))

    return frames


def generate_orfs(input_file):
    """Returns candidate protein strings from open reading frames of DNA bases"""

    # Open the input file (new lines separate different DNA sequences)
    with open(input_file, "r"):
        dna_sequences = [line.rstrip('\n') for line in open(input_file)]

    # Add the reverse complements to the DNA sequences
    dna_sequences += [reverse_complement(dna_sequences)]

    # Transcription: DNA to RNA
    rna_sequences = []
    rna_sequences += [seq.replace('T', 'U') for seq in dna_sequences]

    # Translation: RNA to Proteins
    translations = []
    for bases in rna_sequences:

        # Find the indices of all open reading frames
        orf_indices = find_orfs(bases)
        proteins = ''

        # Find the candidate protein strings for all frames
        for frames in orf_indices:

            # Convert the frame into RNA codons
            frame = bases[frames[0]:frames[1]]
            translate = [frame[i:i + 3] for i in range(0, len(frame), 3)]

            # Translate the RNA codons to proteins
            for i in translate:
                proteins += rna_codons[i]

            # Avoid including duplicates in output
            if proteins not in translations:
                translations += [proteins]

            # Re-initiate proteins for the next frame
            proteins = ''

    return translations


print(generate_orfs('orf_input.txt'))
