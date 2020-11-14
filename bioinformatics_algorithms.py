#This module collects some useful bioinformatics algorithms

'''
Read a (reference) genome file
******************************
Problem: Load a locally saved genome file and read it line by line, stripping any whitespace, tabulators etc.
'''

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

'''
Read short reads in Fastq format
********************************
Problem: Load a fastq file with n short reads line by line, extract and store in two lists the extracted sequences and base qualities
'''

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

'''
Convert phred33 quality scores extracted in the 'readFastq' function to decimal quality score
*********************************************************************************************
Problem: A fastq file, output after sequencing, contains a line per read with the respective base qualities (probability of the base being
correctly sequenced). These qualities are each exactly one byte long, so as to match with the one-byte size of each sequenced nucleotide,
'A', 'T'. 'G', or 'C'. This is done by converting the actual integer-based quality score to an ASCII character.
This function converts this ASCII character back to an integer, i.e. an analyzable quality score. Since the numbers from 0 to 32 have no
character as ASCII equivalent, but [SPACE], [ESC], etc., only ASCII equivalents of 33 and above can be used as byte-long (phred33) quality
scores in a fastq file. Thus we need to subtract 33 from the identified decimal quality score in order to get the actual base quality.
Otherwise we wouldn't be able to output very low quality scores (i.e. quality scores up to 32).
As an example, let's say we identify a '#' hashtag character in the quality line of the fastq file. This corresponds to the decimal number
35 when converting from ASCII to decimal. In reality, we have a quality score of 2 here, since we need to subtract 33 from 35.
'''

def phred33toQ(qual):
    return ord(qual) - 33 # the built-in function 'ord' converts an ASCII character to a decimal number.

'''
Return the reverse complement of an input sequence
**************************************************
Problem: Find the reverse complement (e.g. 'GTAAG' -> 'CTTAC') of a given sequence. **Here 'N' nucleotides due to sequencing errors are 
accounted for as well, so as to avoid dictionary key errors in later use of this function.**
'''

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

'''
Naive exact matching
********************
Problem: Match short reads to a reference genome, e.g. where in the reference genome is the sequence 'ATAAGCC' located.
Naive exact matching iterates through the reference genome character by character in an outer loop and through the read
character by character in an inner loop and checks if the characters (nucleotides) in both sequences match.
'''

def naive(p, t):
    occurrences = [] # initialize an empty list in which the offsets of all matches are going to be stored
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break      
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

'''
Naive exact matching with reverse complement
********************************************
Problem: Same as for the naive exact matching algorithm above, but taking into account the reverse strand (e.g. 'AGC' -> 'GCT') as well,
i.e. matches are accounted for both the input sequence and its reverse complement. It is made use of above 'reverseComplement' code.
'''

def naive_with_rc(p, t):
    occurrences = [] # initialize an empty list in which the offsets of all matches are going to be stored
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break       
        if match:
            occurrences.append(i)  # all chars matched; record
    r = reverseComplement(p) # call the reverseComplement function

    for m in range(len(t) - len(p) + 1): # repeat above loop for the reverse complement
        match = True
        for k in range(len(r)):
          if t[m+k] != r[k]:
                  match = False
                  break
        if match:
          occurrences.append(m)
    # return the 'occurrences' list after converting to a set in order to eliminate doubles which occur if the forward strand equals the reverse
    # strand (e.g. ATAT -> ATAT)
    return sorted(list(set(occurrences))) 

'''
Naive exact matching with n mismatches allowed
**********************************************
Problem: Since DNA sequencers are not perfect, exact matching might not always be the most appropriate way to map sequencing reads to a reference
genome. Here, n mismatches per comparison are allowed, e.g. for 'ATCGGA' and 2 allowed mismatches, 'TTCAGA' would still qualify as
a valid match. **In this algorithm, the reverse complement is not accounted for, but can of course be implemented.**
'''

def naive_with_mismatches(p, t, n):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mm = 0 #initialize the number of mismatches
        match = True
        for j in range(len(p)):  # loop over characters
           if t[i+j] != p[j]:  # compare characters
                mm += 1 #increment the no. of mismatches
                if mm > n: #only n mismatches are allowed
                  match = False
                  break        
        if match:
            occurrences.append(i)  # all chars matched; record

    return occurrences





