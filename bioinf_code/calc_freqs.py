import os
import sys
import getopt
import numpy


if __name__ == "__main__":

    filename = ''

    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "i:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-i":
            filename = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)

    with open(filename, 'r') as handle:

        #read in data without id line and last return line, create variables
        lines = handle.readlines()[1:-1]
        lines_list = []
        A_count, T_count, C_count, G_count, total_count = 0, 0, 0, 0, 0

        #loop over lines and count total and nucleotide counts
        for line in lines:
            new_line = line.replace('\n','')
            new_line2 = new_line.replace("N", "A")
            lines_list.append(new_line2)
            for element in new_line:
                total_count +=1
                if element == "A":
                    A_count += 1
                elif element == "T":
                    T_count += 1
                elif element == "C":
                    C_count += 1
                elif element == "G":
                    G_count += 1

        #calcculate frequencies and print
        A_freq = float(A_count) / total_count
        T_freq = float(T_count) / total_count
        C_freq = float(C_count) / total_count
        G_freq = float(G_count) / total_count
        print "A:", A_freq
        print "C:", C_freq
        print "G:", G_freq
        print "T:", T_freq
        print '\n'


        #dinucleotide frequencies, use dictionary instead much shorter
        seq = ''.join(lines_list)
        dinucleotide_dict = {}
        #loop over length -1 so no loner, add new to dict, count existing
        for i in range(len(seq)-1):
            dinucleotide = seq[i:i+2]
            if dinucleotide not in dinucleotide_dict:
                dinucleotide_dict[dinucleotide] = 1
            elif dinucleotide in dinucleotide_dict:
                dinucleotide_dict[dinucleotide] +=1

        #sum over dict values to get total number of dinucleotides
        total_dinuc = sum(dinucleotide_dict.values())

        #calculate frequencies and print
        for key in dinucleotide_dict:
            freq = float(dinucleotide_dict[key]) / total_dinuc
            print key, ":", freq
