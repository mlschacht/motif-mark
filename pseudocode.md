# Pseudocode for motif mark

## Save the motifs in a list
- go through the motif files
- save them in a list

initialize a longest length holder

## Store important info per sequence
- go through the file and store the info below in a dictionary with the keys as the header lines
    - sequence, length, and header
- also store longest length

## Create a canvas to draw shapes
The total number of length of the keys in the dictionaries is how many sequences in the file so
- multiply the number of keys by 100 (that's the height)
- add longest length + 20 pixels (thats the width of the canvas)


go through the file again
## Draw the lines and exons for each 
- grab the header
- add a line at pixel 0+50 with the start line at 10 and the end at the length of that sequence
- go through the sequence
    - identify the start of the exon and the length of the exon
    - add a box that spans the start (plus 10) to the end with 10 pixels above and below the line
- go through the sequence per motif (5 in this case)
    - each time a motif is found, make a motif class object
    - draw the object
- proceed to the next sequence

- close the image object
-export the image object



motif class holds
location