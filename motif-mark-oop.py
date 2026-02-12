#!/usr/bin/env python

import argparse
import cairo
import re

motifs_file:str = 'Fig_1_motifs.txt'
fasta_file:str = 'Figure_1.fasta'

#initialize variables
motif_list:list = []
longest_gene:int = 0
num_genes:int = 0


gene_dict:dict = {}
header:str = ""
sequence:str = ""
seq_length:int = 0


with open(motifs_file, "r") as m_file: #write all motifs into a list
    for line in (m_file):
        line = line.strip()
        motif_list.append(line)

with open(fasta_file, "r") as f_file: #grab important info about each sequence ans store it in a dictionary
    for line in (f_file):
        line = line.strip()
        if line.startswith(">") == True and num_genes != 0: #at the header lines
            if len(sequence) > longest_gene: #update longest sequence if this is the longest sequence
                longest_gene = len(sequence)
                
            #this gene is complete, update the dictionary
            gene_dict[header] = (sequence)

            #reset the values for the next gene
            header = line #store the header for this gene
            sequence = "" #reset sequence for this gene
            num_genes += 1
            print(num_genes)
            continue
        
        if line.startswith(">") == True and num_genes == 0: #at the first header line
            header = line #store the header for this gene
            num_genes += 1
            print(num_genes)
            print("start")
            continue

        else:
            sequence = sequence + line

    #at the end of the file
    if len(sequence) > longest_gene: #update longest sequence if the last sequence is the longest sequence
                longest_gene = len(sequence)
                
    #this gene is complete, update the dictionary
    gene_dict[header] = (sequence)

for value in gene_dict.values():
     print(value)


class Motif:
    def __init__(self, the_name, the_sex, the_language):
        '''This is how a human is made.'''

        # Data
        self.name = the_name
        self.sex = the_sex
        self.language = the_language

    # Methods
    def introduce(self):
        print(f'My name is {self.name} and I speak {self.language}.')