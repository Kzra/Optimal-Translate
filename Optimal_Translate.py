#!usr/bin/env python3

def transcribe(sequence,rev): #turns sequence into a list, cycles through letters and produces reverse complement, rejoins list
 if sequence.find("U") == -1:  #if you have dna 
    if rev == 1: #if the sequence is reversed we need to generate the reverse complement, otherwise we should just replace T for U
     seq = list(sequence) 
     rna_seq = ['Blank']*len(seq)
     rc = range(0,len(seq))
     for c,dna in zip(rc,seq):
        if dna == 'A':
           rna_seq[c] = 'U'
        if dna == 'T':
           rna_seq[c] = 'A'
        if dna == 'G':
           rna_seq[c] = 'C'
        if dna == 'C':
           rna_seq[c] = 'G'
     rna_seq = ("".join(rna_seq))
    else: 
     rna_seq = sequence.replace('T','U') 
 else:
    if rev == 1: #also the reverse complement for RNA
     seq = list(sequence) 
     rna_seq = ['Blank']*len(seq)
     rc = range(0,len(seq))
     for c,rna in zip(rc,seq):
        if rna == 'A':
           rna_seq[c] = 'U'
        if rna == 'U':
           rna_seq[c] = 'A'
        if rna == 'G':
           rna_seq[c] = 'C'
        if rna == 'C':
           rna_seq[c] = 'G'
     rna_seq = ("".join(rna_seq))
    else: 
        rna_seq = sequence
 return rna_seq

def translate(sequence): #this creates a dictionary {} which is essentially a hash table linking codon to amino acid.
     
    codon2aa = {"AAA":"K", "AAC":"N", "AAG":"K", "AAU":"N", 
                "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T", 
                "AGA":"R", "AGC":"S", "AGG":"R", "AGU":"S", 
                "AUA":"I", "AUC":"I", "AUG":"M", "AUU":"I", 

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAU":"H", 
                "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P", 
                "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R", 
                "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L", 

                "GAA":"E", "GAC":"D", "GAG":"E", "GAU":"D", 
                "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A", 
                "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G", 
                "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V", 

                "UAA":"*", "UAC":"Y", "UAG":"*", "UAU":"T", 
                "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S", 
                "UGA":"*", "UGC":"C", "UGG":"W", "UGU":"C", 
                "UUA":"L", "UUC":"F", "UUG":"L", "UUU":"F"}
    
    r = int(len(sequence)/3)
    codon = ['None']*r
    prim_seq = ['None']*r
    codon[0] = sequence[0:3]
    
    for i in range(1, r): #here range generates integers from 1 up to r (but not including r)
        n = i*3
        codon[i] = sequence[n:n+3] #in python between (:) excludes the last integer
    

    for i in range(0, r):
        prim_seq[i] = codon2aa[codon[i]]
    
    
    translation = ("".join(prim_seq))
    return translation

##This function transcribes DNA and then counts the number of stopcodons in each reading frame (1 - 6) and finally translates the reading frame that produces the  least number of stop codons. 
def count_stop_codons(sequence):
    
    nsequence = transcribe(sequence,0)
    rsequence = transcribe(sequence,1)
                
   
    frame_1 = nsequence 
    frame_2 = nsequence[1:len(sequence)]
    frame_3 = nsequence[2:len(sequence)]
    
    reverse_sequence = rsequence[::-1] #this is a function of [] indexing called 'slicing'
    
    
    rframe_1 =reverse_sequence
    rframe_2 =reverse_sequence[1:len(reverse_sequence)]
    rframe_3 =reverse_sequence[2:len(reverse_sequence)]
    frames = [frame_1,frame_2,frame_3,rframe_1,rframe_2,rframe_3]
    reading_frame =['Frame_1','Frame_2','Frame_3','rFrame_1','rFrame_2','rFrame_3']
    
    stop_count = ['None']*6
    
    count = [0,1,2,3,4,5]
    
    for f,a in zip(frames, count): #Zip is a nice command that allows simultaneous iteration of two same sized variables
        stopcodon = ['UAA','UGA','UAG']
        b = int(len(f)/3)
        codon = ['None']*b
        codon[0]=f[0:3]
        
        for i in range(1, b):
            c =i*3
            codon[i] = f[c:c+3]
        
        
        stopcodon[2]=  codon.count('UAG')
        stopcodon[1] = codon.count('UGA')
        stopcodon[0] = codon.count('UAA')
        
        stop_count[a] = sum(stopcodon)
    
    
    min_value = min(stop_count)
    min_index = ([i for i, x in enumerate(stop_count) if x == min_value]) #.index returns the first index value from a find-type function, as is standard output in MATLAB. Here I want to see if there are multiple minimum values.
    # a note on enumerate, the syntax is for counter,value. So in the above x takes on the value of the stop_counts and i the index position. The command creates a list by expanding the [] printing i for each iteration that x = minvalue.
   
    
    optimal_translation = ['Blank']*len(min_index)
    optimal_frame = ['Blank']*len(min_index)
    tc = range(0, len(min_index))
    
    for v,index in zip(tc,min_index):
        optimal_translation[v] = translate(frames[index])
        optimal_frame[v] = reading_frame[index]  
    return optimal_translation, optimal_frame

#Here I work on combining them all 
#Now I want to work on importing .fasta files and translating them 
# sequence contains the sequences and sequence_name the names

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("DNA_fasta_file", help="[Input]: list of DNA sequences to be translated")
parser.add_argument("output_fasta_file", help ="[Output]: list of translated sequences")  
args = parser.parse_args()
ln = 1*10**7
with open(args.DNA_fasta_file, "rt") as data: #with open automatically calls file close
    

    text = data.read()
    sp = ([i for i, x in enumerate(text) if x == ">"]) #> demarks new sequences
    sequence_name = ['Blank']*len(sp)
    sequence = ['Blank']*len(sp)
    if sp == [0]: 
        query = text
        #print(query)
        sequence_name = query.splitlines()[0]
        #print(sequence_name)
        sequence = query.splitlines()[1:ln] #obviously there are fewer lines than the ln but this works
        #sequence = sequence[0] #removes double list
        #print(sequence)

   ##if there is more than one
    else:
        for ind in range(0,len(sp)-1):
            query = (text[sp[ind]:sp[ind+1]])
            #print(query)
            sequence_name[ind] = query.splitlines()[0]
            #print(sequence_name)
            sequence[ind] = query.splitlines()[1:ln] #obviously there are fewer lines than the len(sequence) but this works
            #print(sequence)
            sequence[ind] = sequence[ind][0] #removes double list
   
        
        query = (text[sp[len(sp)-1]:len(text)]) #deal with the last sequence outside of the for loop
        sequence_name[len(sp)-1] = query.splitlines()[0]
        sequence[len(sp)-1] = query.splitlines()[1:len(sequence)]
        sequence[len(sp)-1] = sequence[len(sp)-1][0]
optimal_translation = ['Blank']*len(sp)
reading_frame = ['Blank']*len(sp)
count = range(0,len(optimal_translation))


for ct,sq in zip(count,sequence):
    optimal_translation[ct] = count_stop_codons(sq)[0]
    reading_frame[ct] = count_stop_codons(sq)[1]

    

    
with open(args.output_fasta_file,'w') as output:

     for name,frame,trans in zip(sequence_name,reading_frame,optimal_translation):
         for f in range(0, len(trans)):
             temp = [name,frame[f]]
             temp = "_".join(temp)
             output.write(temp)
             output.write('\n')   
             output.write(trans[f])
             output.write('\n')   
         output.write('\n')     
    
        
   

