#/usr/bin/python3

# Concatenate the sequences of multiple-sequence alignments with the same taxa present in a single directory.

import os
from Bio import AlignIO

my_files = os.listdir('.')

my_seqs = {}

for file in my_files:

    if file.endswith(".fasta"):

        print(file)

        alignment = AlignIO.read(open(file), "fasta")

        for record in alignment:

            my_seqs[record.id] = ""


for file in my_files:

    if file.endswith(".fasta"):

        alignment = AlignIO.read(open(file), "fasta")

        unknown = "-"*alignment.get_alignment_length()

        my_records = {}

        for record in alignment:

            my_records[record.id] = record.seq

        for seq in my_seqs.keys():

            if seq in my_records.keys():

                new_seq = my_seqs[seq] + my_records[seq]

                my_seqs[seq] = new_seq

            elif seq not in my_records:
                
                new_seq = my_seqs[seq] + unknown

                my_seqs[seq] = new_seq

for key in my_seqs.keys():

    print(">"+key)
    print(my_seqs[key],"\n")


        
        
