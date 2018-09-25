import os
import sys
import getopt
import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


if __name__ == "__main__":

    filename = ''
    outfilename = ''

    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "i:o:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-i":
            filename = arg
        elif opt == "-o":
            outfilename = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)

    with open(filename, 'r') as handle, open(outfilename, 'w') as outfile:

        #read in data without id line
        lines = handle.readlines()[1:]
        #create empty list for complements
        complements = []
        #get complements and fill in list
        for line in lines:
            new_line = line.strip('\n')
            for element in new_line:
                if element == "A":
                    complements.append("T")
                elif element == "T":
                    complements.append("A")
                elif element == "C":
                    complements.append("G")
                elif element == "G":
                    complements.append("C")

        #reverse complements list, add fasta id, print to fasta file
        reverse_comp_list = complements[::-1]
        reverse_comp = ''.join(reverse_comp_list)
        seq = SeqRecord(Seq(reverse_comp), id = "reversed", description = '')
        SeqIO.write(seq, outfile, "fasta")
