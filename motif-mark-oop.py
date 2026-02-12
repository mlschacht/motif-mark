#!/usr/bin/env python

import argparse
import cairo
import re

motifs_file:str = 'Fig_1_motifs.txt'
fasta_file:str = 'Figure_1.fasta'

motif_list:list = []

with open(motifs_file, "r") as m_file: #write all motifs into a list
    for line in (m_file):
        line = line.strip()
        motif_list.append(line)




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