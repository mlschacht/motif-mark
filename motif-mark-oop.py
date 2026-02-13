#!/usr/bin/env python

import argparse
import cairo
import re
import bioinfo

#update with argparse when ready
motifs_file:str = 'Fig_1_motifs.txt'
fasta_file:str = 'Figure_1.fasta'
prefix:str = fasta_file.split(".")[0] #save the prefix of the input fasta file to use as the output png file name


#initialize variables
motif_dict:dict = {} #keys are regex motifs and values are the original sequence lengths
longest_gene:int = 0
num_genes:int = 0


gene_dict:dict = {}
header:str = ""
sequence:str = ""
seq_length:int = 0

is_DNA_file = bioinfo.validate_base_file(fasta_file) #True if fasta is a DNA file. False if RNA

class Motif:
    def __init__(self, start, length, gene_position):
        '''Take in the regex motif, the motif start on the gene, the length of the .'''

        # Data
        self.start = start
        self.length = length
        self.gene_position = gene_position

    # Methods
    def draw_motif(self):
        #draw the motif based on the start and the length
        context.rectangle(self.start, self.gene_position, self.length, 20)        #(x0, y0, length, height)
        context.set_source_rgba(0,1,0, 0.5)
        context.fill()


with open(motifs_file, "r") as m_file: #write all motifs into a list
    for line in (m_file):
        seq = line.strip()

        #convert the sequence to DNA or RNA to match what is seen in the sequencing file
        is_DNA_motif = "u" not in line.lower()
        if is_DNA_file != is_DNA_motif:
            seq = bioinfo.convert_DNA_RNA(seq, DNAflag=is_DNA_motif)
        
        # motif to regex conversion
        print(seq)
        regex_motif = bioinfo.motif_to_regex(seq, DNAflag = is_DNA_file)

        motif_dict[regex_motif.lower()] = len(seq) #add the regex motif to the motif list

print(motif_dict)

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
            continue
        
        if line.startswith(">") == True and num_genes == 0: #at the first header line
            header = line #store the header for this gene
            num_genes += 1
            continue

        else:
            sequence = sequence + line

    #at the end of the file
    if len(sequence) > longest_gene: #update longest sequence if the last sequence is the longest sequence
                longest_gene = len(sequence)
                
    #this gene is complete, update the dictionary
    gene_dict[header] = (sequence)


# draw the image canvas based on the number of genes and the longest gene length
height:int = num_genes *100
width:int = longest_gene + 20
image_file_name:str = prefix + ".png"
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
context = cairo.Context(surface)

x_start:int = 50 #initialize the starting point for each gene
y1:int = 40     #10 pixels above the line
text_start:int = 20 #initialize the starting point for each header

#go through the dictionary
for header in gene_dict.keys():
    #write out the header for each gene
    context.set_font_size(12)
    context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.move_to(10, text_start) #header start
    context.show_text(header)
    context.stroke() #write out the header

    sequence = gene_dict[header]
    #draw the length of the gene
    seq_length = len(sequence)
    context.set_line_width(2)
    context.move_to(10,x_start)
    context.line_to(seq_length, x_start)
    context.set_source_rgb(0,0,0)
    context.stroke()
    x_start += 100

    #find the exon position and length
    exon_find = re.search(r'[A-Z]+', sequence)
    exon_start:int = exon_find.start() + 1 # add 1 because python uses index 0 so the start of the sequence would be 0 instead of 1
    exon_length:int = len(exon_find.group())
    exon_end = exon_start + exon_length

    #draw the exon
    context.rectangle(exon_start, y1, exon_length, 20)        #(x0, y0, length, height)
    context.set_source_rgba(0,0,0, 1)
    context.fill()


#---------------START HERE AFTER MAKING THE MOTIF TO REGEX CONVERSION FUNCTION IN BIOINFO-----------------------------
    #find the motif in the sequence
    for motif in motif_dict.keys():
        motif_length = motif_dict[motif]
        for i, char in enumerate(sequence):
             check = sequence[i:i+motif_length]
             if re.match(motif, check):
                motif_obj = Motif(i+1, motif_length, y1)
                motif_obj.draw_motif()
                print("motif drawn")


    #move the y positions to the next exon
    y1 += 100
    text_start += 100





surface.write_to_png(image_file_name)

