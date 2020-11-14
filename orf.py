import Bio
from Bio import SeqIO
from textwrap import wrap

'''
def orf(fasta_file, reading_frame):
    #creates a list of all open reading frames (ORF), i.e. triplets in a DNA sequence that contain a start codon and a stop codon and therefore
    #have the potential for protein encoding
    #The 'reading_frame' argument goes from 0 to 2, i.e. 0 represents reading frame 1, 1 is r.f. 2, and 2 is r.f. 3
    parsed_fasta = SeqIO.parse(fasta_file, "fasta")
    seq_list = [seq_record.seq for seq_record in parsed_fasta]
    orf_list = []
    for seq in seq_list:
        orf_in_seq = []
        wrapped_seq = wrap(str(seq[reading_frame:]),3) 
        i = 0
        while i < len(wrapped_seq):
            if wrapped_seq[i] == 'ATG':
                s = wrapped_seq[i]
                i+=1
                s = s + wrapped_seq[i]
                while wrapped_seq[i] not in ['TAG', 'TAA', 'TGA']:
                    i+=1
                    if i < len(wrapped_seq):
                        s = s + wrapped_seq[i]
                    else:
                        s = ''
                        break
                orf_in_seq.append((len(s),(i+1)*3-len(s)+1+reading_frame,len(wrapped_seq)*3))  
                #orf_in_seq.append(len(s))
            i+=1
        orf_list.append(orf_in_seq)
    #The record ids need to be parsed separately
    parsed2 = SeqIO.parse(fasta_file, "fasta")
    id_list = [seq_record.id for seq_record in parsed2]
    #append the record id to each sublist of ORFs
    for j in range(len(id_list)):
       orf_list[j].append(id_list[j])
    #for elem in orf_list:
        #if len(elem) > 1:
            #print (max(elem[:-1]))
    return orf_list

print(orf("dna2.fasta",1))
'''

#Alternate version with dict

def orf_alt(fasta_file, reading_frame):
    #creates a list of all open reading frames (ORF), i.e. triplets in a DNA sequence that contain a start codon and a stop codon and therefore
    #have the potential for protein encoding
    #The 'reading_frame' argument goes from 0 to 2, i.e. 0 represents reading frame 1, 1 is r.f. 2, and 2 is r.f. 3
    parsed_sequences = SeqIO.parse(fasta_file, "fasta")
    seq_list = [seq_record.seq for seq_record in parsed_sequences]
    parsed_ids = SeqIO.parse(fasta_file, "fasta")
    id_list = [seq_record.id for seq_record in parsed_ids]
    orf_list = {}
    for j in range(len(seq_list)):
        orf_in_seq = {}
        wrapped_seq = wrap(str(seq_list[j][reading_frame:]),3) 
        i = 0
        while i < len(wrapped_seq):
            if wrapped_seq[i] == 'ATG':
                s = wrapped_seq[i]
                i+=1
                s = s + wrapped_seq[i]
                while wrapped_seq[i] not in ['TAG', 'TAA', 'TGA']:
                    i+=1
                    if i < len(wrapped_seq):
                        s = s + wrapped_seq[i]
                    else:
                        s = ''
                        break
                if len(s) > 0:
                    orf_in_seq[len(s)] = (i+1)*3-len(s)+1+reading_frame
            i+=1
        orf_list[id_list[j]] = orf_in_seq
    return orf_list

res=orf_alt("dna2.fasta",2)

#print (res)

#find longest ORF in the respective reading frames (1, 2 or 3)
longest = 0
for id in res:
    if len(res[id]) > 0:
        if max(res[id]) > longest:
            longest = max(res[id])
            id_of_longest = id
print ((id_of_longest,longest))
#identify starting position of ORF
print (res[id_of_longest][longest])




#The following is a test function without outer loop:

'''
def orf(fasta_file):
    #creates a list of all open reading frames (ORF), i.e. triplets in a DNA sequence that contain a start codon and a stop codon and therefore
    #have the potential for protein encoding
    parsed_fasta = SeqIO.parse(fasta_file,"fasta")
    seq_list = [seq_record.seq for seq_record in parsed_fasta]
    id_list = [seq_record.id for seq_record in parsed_fasta]
    orf_list = []
    seq = seq_list[2]
    wrapped_seq = wrap(str(seq), 3)
    #wrapped_seq.append('ATG')
    wrapped_seq.append('TAG')
    print (wrapped_seq)
    i = 0
    while i < len(wrapped_seq):
        if wrapped_seq[i] == 'ATG':
            s = wrapped_seq[i]
            i+=1
            s = s + wrapped_seq[i]
            while wrapped_seq[i] not in ['TAG', 'TAA', 'TGA']:
                i+=1
                s = s + wrapped_seq[i]
                #print (s)
            orf_list.append(s)
            #print (orf_list) 
        i+=1
    return orf_list

#seq = 'AAACCCAAATTTATGAAAAAAAAACCCTTTTAGAAAGGGTTTATGCCCGGGTGACCCCCCCCCCCCATGGGGGGGGGGTAGCCCATGCCCTAAAAAAAAAAAAACCCCCCCCAAAAAAAAAAAAAAAAAAATGAAACCCTAGAAA'
print(orf("dna2.fasta"))
'''

