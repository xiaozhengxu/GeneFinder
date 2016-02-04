# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Xiaozheng Xu
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """
    Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way 
    """
    return ''.join(random.sample(s, len(s)))

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide=='A':
        return 'T'
    elif nucleotide=='T':
        return 'A'
    elif nucleotide=='G':
        return 'C'
    elif nucleotide=='C':
        return 'G'
    else:
        return None

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_dna=dna[::-1]    #First, reverse the string
    #this uses the slicing method: [start:end:step] to reverse string
    # can also use ''.join(reversed(a))
    reverse_complement=''
    for char in reverse_dna:
        reverse_complement=reverse_complement+get_complement(char) #get compliment for each char in string
    return reverse_complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        Stop codons are 'TAG','TAA','TGA'
        dna: a DNA sequence
        returns: the open reading frame represented as a string

        added Test case of string with no stop codon and 
        Test case of string with a different stop codon than TAG
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGAGGGGGGG")
    'ATGAGAGGGGGGG'
    >>> rest_of_ORF("ATGAGAGCTTGA")
    'ATGAGAGCT'
    """
    stop_codons=['TAG','TAA','TGA']
    for i in range((len(dna))/3): # break it down into segments of 3
        for codon in stop_codons:
            segment=dna[i*3:i*3+3] 
            #print(segment)
            if segment==codon:   # check if segment is a stop codon
                return dna[:i*3]   
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I added test cases for string with the first 'ATG' in a different frame, string with no 'ATG' in frame, 
        and string with no 'ATG' at all 
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    >>> find_all_ORFs_oneframe('TGCATGAATGTAG')
    ['ATGAATGTAG']

    >>> find_all_ORFs_oneframe('TGCCTGAATGTAG')
    []

    >>> find_all_ORFs_oneframe('GCGTGACTGTAGAT')
    []

    """
    stop_codons=['TAG','TAA','TGA']
    start_codons=['ATG']
    result_list=[]
    i=0;
    while i<len(dna)/3: # break it down into segments of 3
            segment=dna[i*3:i*3+3] 
            if segment=='ATG':   
                ORF=rest_of_ORF(dna[i*3:])
                result_list.append(ORF)
                i+=len(ORF)/3           
            else:
                i+=1

    return result_list

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        The test cases 
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # result_list=[]
    # for i in range(3):
    #     ORF_in_frame=find_all_ORFs_oneframe(dna[i:])
    #     result_list.extend(ORF_in_frame)   #extend add elements of a list to a new list
    result_list=[find_all_ORFs_oneframe(dna[i:]) for i in range(3)] # a appended nested list 
    result_list=[s for l in result_list for s in l] # A way to flatten a nexted list: for lists in result_list, for string in each sublists
    return result_list


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    ORF_in_original=find_all_ORFs(dna)
    result_list=ORF_in_original
    ORF_in_complement=find_all_ORFs(get_reverse_complement(dna))
    result_list.extend(ORF_in_complement)

    return result_list


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    ORFs=find_all_ORFs_both_strands(dna)
    # maxlength=max(len(s) for s in ORFs)  #not using [] makes max operate on iterative instead of creating a full list before hand
    # longest=[s for s in ORFs if len(s)==maxlength]
    # return longest[0]
    length=0
    for s in ORFs:
        if len(s)>length:
            length=len(s)
            longest=s
    return longest

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    length=0
    for i in range(num_trials):
        shuffled_dna=shuffle_string(dna)
        longest=longest_ORF(shuffled_dna)
        if len(longest)>length:
            length=len(longest)
    return length

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    result_string=''
    for i in range((len(dna))/3): # break it down into segments of 3
            codon=dna[i*3:i*3+3] 
            amino_acid=aa_table[codon]
            result_string=result_string+amino_acid 
    return result_string


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    print threshold
    ORFs=find_all_ORFs_both_strands(dna)
    return [coding_strand_to_AA(ORF) for ORF in ORFs if len(ORF)>=threshold]


print gene_finder(dna)

# if __name__ == "__main__"
import doctest
#     doctest.run_docstring_examples(rest_of_ORF, globals())

doctest.testmod()
