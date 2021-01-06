#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import re
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Alexis Gouthey"
__copyright__ = "Universite CyTech"
__credits__ = ["Alexis Gouthey"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Alexis Gouthey"
__email__ = "goutheyale@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """
    Check if path is an existing file.

    Parameters:
        path: (String) Path to the file

    Returns: The path of the file
    """

    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """
    Retrieve the arguments of the program.

    Returns: An object that contains the arguments
    """

    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
    Read the file and return every line that contain nucleic sequences.

    Parameters:
        amplicon_file: (String) Path to the fastq_file
        minseqlen: (Int) Minimum size for a nucleic sequence

    Returns: A generator of nucleic sequences
    """

    # Open the file in "read" mode
    file = open(amplicon_file, "r")

    # String which will contains nucleic sequences.
    # A single nucleic sequence can be written on multiple lines
    string = ""

    for line in file:
        # Match the line with regular expression containing
        # nucleic sequences characters
        if re.match(r'[ATCG]+', line):
            # Add the line to the string and remove the "\n" at the end
            string += line.rstrip("\n")
        else:
            # Only yield sequences which have a length greater than minseqlen
            if len(string) >= minseqlen:
                yield string

            # Reset the variable
            string = ""

    # Yield the last sequence if the file ends with a nucleic sequence
    if len(string) >= minseqlen:
        yield string

    file.close()

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    Filter nucleic sequences from read_fasta and return it with number of
    occurence in descending order.

    Parameters:
        amplicon_file: (String) Path to the fastq_file
        minseqlen: (Int) Minimum size for a nucleic sequence
        mincount: (Int) Minimum number of occurence

    Returns: A generator of nucleic sequences with their number of occurence
    """

    liste = []

    for line in read_fasta(amplicon_file, minseqlen):
        liste.append(line)

    # Filter sequences that appear less than mincount times
    list_mincount = list(filter(lambda x: liste.count(x) >= mincount, liste))

    # Remove duplicates
    list_mincount = list(set(list_mincount))

    # Transform the list by adding number of occurence
    list_mincount = list(map(lambda x: [x, liste.count(x)], list_mincount))

    # Sort the list by occurence in descending order before yielding
    list_mincount.sort(key = lambda x : x[1], reverse = True)

    for elem in list_mincount:
        yield elem

def get_chunks(sequence, chunk_size):
    """
    Divide sequence in 4 or more sub sequences.

    Parameters:
        sequence: (String) Nucleic sequence
        chunk_size: (Int) Size of cut sequences

    Returns: A list containg sub sequences
    """

    if chunk_size * 4 > len(sequence):
        raise ValueError
    else:
        list = []

        for i in range(int(len(sequence) / chunk_size)):
            list.append(sequence[i * chunk_size : (i + 1) * chunk_size])

        return list

def get_unique(ids):
    """
    Get unique ids in a dictionnary

    Parameters:
        ids: (List) List of ids

    Returns: A dictionnary which contains unique ids
    """

    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    """
    Get shared elements between 2 lists

    Parameters:
        lst1: (List) List of elements
        lst2: (List) List of elements

    Returns: A list containg shared elements between the two lists
    """

    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """
    Cut a nucleic sequence with a kmer_size and return every cut sequences.

    Parameters:
        sequence: (String) A nucleic sequence
        kmer_size: (Int) Size of returned kmer

    Returns: A generator of nucleic sequences with a size of kmer_size
    """

    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i : i + kmer_size]

def get_identity(alignment_list):
    """
    Compute the similarity rate between 2 sequences.

    Parameters:
        alignment_list: (List/Tuple) List or Tuple of sequences

    Returns: Similarity rate between sequences
    """

    shared_nucleotids = 0

    if len(alignment_list[0]) == len(alignment_list[1]):
        for i in range(len(alignment_list[0])):
            # Compare every caracter between 2 sequences
            if alignment_list[0][i] == alignment_list[1][i]:
                shared_nucleotids += 1

    # Divide the result by the length of one of the sequence
    return shared_nucleotids / len(alignment_list[0]) * 100

def compute_similarity_matrix(chunks, parents_sequence, non_chimeric_list, chunk_size):
    """
    Compute the similarity matrix.

    Parameters:
        chunks: (List) Chunks of the current sequence
        parents_sequence: (List) Parent sequences from the sequence
        non_chimeric_list: (List) List of non_chimeric sequences
        chunk_size: (Int) Size of a chunk

    Returns: The similarity matrix
    """

    # Initialize the matrix
    perc_identity_matrix = [[] for _ in range(len(chunks))]

    for parent_sequence in parents_sequence:
        # Get sub sequences from non_chimeric_sequence list with a size of chunk_size
        non_chimeric_chunk = get_chunks(non_chimeric_list[parent_sequence], chunk_size)

        for index, chunk in enumerate(chunks):
            # Compute alignement list between current chunk and chunk from the non_chimeric list
            global_alignement = nw.global_align(chunk, non_chimeric_chunk[index], gap_open=-1,
            gap_extend=-1, matrix = os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))

            # Compute the similarity rate from the alignement list
            similarity_rate = get_identity(global_alignement)

            # Add the result to the matrix
            perc_identity_matrix[index].append(similarity_rate)

    return perc_identity_matrix

def search_mates(kmer_dict, sequence, kmer_size):
    """
    Get the 8 most common nucleic sequence from the kmer_dict.

    Parameters:
        kmer_dict: (Dict) A kmer dictionnary
        sequence: (String) A nucleic sequence
        kmer_size: (Int) Size for cutting kmer

    Returns: A list
    """

    liste = []

    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            for i in Counter(kmer_dict[kmer]).most_common(8):
                liste.append(i[0])

    return liste

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    Add every unique kmer contained in sequence to kmer_dict with id_seq.

    Parameters:
        kmer_dict: (Dict) Dictionnary of kmer
        sequence: (String) A nucleic sequence
        id_seq: (Int) Id of the nucleic sequence
        kmer_size: (Int) Size for cutting nucleic sequences

    Returns: A dictionnary containing new unique sequences
    """

    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict.keys():
            kmer_dict[kmer] = []

        kmer_dict[kmer].append(id_seq)

    return kmer_dict

def detect_chimera(perc_identity_matrix):
    """
    Check if a candidate nucleic sequence is chimeric or not.

    Parameters:
        perc_identity_matrix: (Matrix) Matrix which indicates similarity between
                              one sequence and its parent sequences

    Returns: A boolean
    """

    standard_deviation = 0
    sequence1_similarity, sequence2_similarity = set(), set()

    for percentage in perc_identity_matrix:
        standard_deviation += statistics.stdev(percentage)
        sequence1_similarity.add(percentage[0])
        sequence2_similarity.add(percentage[1])

    if len(sequence1_similarity) >= 2 or len(sequence2_similarity) >= 2:
        mean_standard_deviation = standard_deviation / len(perc_identity_matrix)
        if mean_standard_deviation > 5:
            return True

    return False

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    For each sequence, we divide it in 4 segments. For each segment, we take
    the 8 most common sequences which are similar to the segment.
    If at least 2 sequences are among the 8 most common, we'll call it "parent sequences"
    We compute similarity scores between parents sequences and current sequence.
    Finally, current sequence will be considered as chimeric if the standard deviation
    is greater than 5 and if at least 2 segments from the sequence have a different
    similarity from 1 out of 2 parents

    Parameters:
        amplicon_file: (String) Path to the fastq_file
        minseqlen: (Int) Minimum size for a nucleic sequence
        mincount: (Int) Minimum number of occurence
        chunk_size: (Int) Size of a chunk
        kmer_size: (Int) Size for cutting kmer

    Returns: A generator of non chimeric nucleic sequences with their number of occurence
    """

    kmer_dict = {}

    non_chimeric_sequence, perc_identity_matrix = [], []

    # Start with an id from 0 and increase it during the loop
    id_seq = 0

    for sequence, occurence in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        # Divide sequence in 4 parts
        chunks_list = get_chunks(sequence, chunk_size)[:4]

        # Get chunk mates
        chunk_mates = []
        for sub_sequence in chunks_list:
            chunk_mates.append(search_mates(kmer_dict, sub_sequence, kmer_size))

        parents_sequence = []
        for chunk in parents_sequence:
            parents_sequence = common(parents_sequence, chunk)

        # If at least 2 sub_sequences from the current sequence have differents
        # similarities from 1 out of 2 parents
        if len(parents_sequence) >= 2:
            perc_identity_matrix = compute_similarity_matrix(chunks_list,
                parents_sequence[:2], non_chimeric_sequence, chunk_size)

        if not detect_chimera(perc_identity_matrix):
            # Get fresh kmer_dict with sub_sequences from sequence
            kmer_dict = get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size)

            non_chimeric_sequence.append(sequence)

            # Incrementing id for next loops
            id_seq += 1

            yield [sequence, occurence]

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    Parameters:
        amplicon_file: (String) Path to the fastq_file
        minseqlen: (Int) Minimum size for a nucleic sequence
        mincount: (Int) Minimum number of occurence
        chunk_size: (Int) Size of a chunk
        kmer_size: (Int) Size for cutting kmer

    Returns: A list of OTU (sequence) with their respective number of occurence
    """

    liste = []

    for test in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
        liste.append(test)

    return liste

def fill(text, width = 80):
    """
    Split text with a line return to respect fasta format
    """

    return os.linesep.join(text[i : i + width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """
    Save OTUs in a file.

    Parameters:
        OTU_list: A list of OTU
        output_file: (String) Path to the output file
    """

    # Open the file in "write" mode
    file = open(output_file, "w")

    for index, tuple in enumerate(OTU_list):
        file.write(">OTU_" + str(index + 1) + " occurrence:" + str(tuple[1]) + "\n")
        file.write(fill(tuple[0]) + "\n")

    file.close()

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """

    # Get arguments
    args = get_arguments()

    fasta_file = args.amplicon_file
    output_file = args.output_file

    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = args.chunk_size
    kmer_size = args.kmer_size

    # Get OTU list
    OTU_list = abundance_greedy_clustering(fasta_file, minseqlen,
               mincount, chunk_size, kmer_size)

    # Write OTU in a file
    write_OTU(OTU_list, output_file)

if __name__ == '__main__':
    main()
