"""
Created on Sun Feb  2 10:10:24 2015
@author: Ziyi Lan
"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table
import random
from load import load_seq
dna = load_seq("./data/X73525.fa")

def shuffle_string(s):
    """ Shuffles the characters in the input string
        NOTE: this is a helper function, you do not have to modify this in any way """
    return ''.join(random.sample(s,len(s)))

### YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = nucleotide 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)





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
    # TODO: implement this
    return get_complement(dna)[::-1]

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    for i in range(0,len(dna),3):
        if dna[i:i+3] in ["TAG", "TAA", "TGA"]:
            return dna[:i]

    return dna

# print rest_of_ORF("ATGAGATAGG")

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    ORFs=[]
    i=0
    while i<len(dna):
        if dna[i:i+3]=="ATG":
            j=i
            dna2=rest_of_ORF(dna[j:])
            ORFs.append(dna2)
            i=i+len(dna2)
        i=i+3
    return ORFs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATG', 'ATGAATGTAG', 'ATGCATGAATGTAG']
    """
    # TODO: implement this
    all_ORFs=[]
    for i in range(3):
        all_ORFs+=find_all_ORFs_oneframe(dna[i:])

    all_ORFs=list(set(all_ORFs))
    return all_ORFs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    strand2=get_reverse_complement(dna)
    strands=find_all_ORFs(dna)
    strands+=find_all_ORFs(strand2)
    return strands


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    # return max(find_all_ORFs_both_strands(dna), key=len)
    list=find_all_ORFs_both_strands(dna)
    output=""
    if len(list)>=1:
        output=max(find_all_ORFs_both_strands(dna), key=len)
    return output

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    a=[]
    for i in range(num_trials):
        a.append(longest_ORF(shuffle_string(dna)))
        # print a
        # print a
    return len(max(a, key=len))

        



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
    # TODO: implement this
    proteinsequence=''
    for n in range(0,len(dna),3):
        if dna[n:n+3] in aa_table.keys():
            proteinsequence += aa_table[dna[n:n+3]]
    return proteinsequence

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna, 1500)
   
    a=[]
    for i in find_all_ORFs_both_strands(dna):
        if len(i)>=threshold:
            a.append(coding_strand_to_AA(i))

    return a

if __name__ == "__main__":
    import doctest
    doctest.testmod()

print gene_finder(dna)
    