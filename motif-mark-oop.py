#!/usr/bin/env python

import argparse
import cairo
import re

motifs_file:str = 'Fig_1_motifs.txt'
fasta_file:str = 'Figure_1.fasta'

#initialize variables
motif_list:list = []
longest_gene:int = 0


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
        if line.startswith(">") == True: #at the header lines
            if len(sequence) > longest_gene: #update longest sequence if this is the longest sequence
                longest_gene = len(sequence)
                print(longest_gene)
                break

            header = line #store the header for this gene
            sequence = "" #reset sequence for this gene

        else:
            sequence = sequence + line



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