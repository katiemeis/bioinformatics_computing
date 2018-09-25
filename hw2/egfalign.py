import os
import sys
import getopt
import numpy

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

    with open(filename, 'r') as handle, open(outfilename, 'r') as outfile:
        #read in data without id line, append lines to get seq
        human_lines = handle.readlines()[1:]
        human_lines_list = []
        human_lines_list.append('-')
        fly_lines = outfile.readlines()[1:]
        fly_lines_list = []
        fly_lines_list.append('-')
        for line in human_lines:
            new_line = line.strip('\n')
            new_line2 = new_line.replace("N", "A")
            human_lines_list.append(new_line2)
        for line in fly_lines:
            new_line = line.strip('\n')
            new_line2 = new_line.replace("N", "A")
            fly_lines_list.append(new_line2)
        human_seq = ''.join(human_lines_list)
        fly_seq = ''.join(fly_lines_list)

        align_table = numpy.zeros((len(fly_seq), len(human_seq)))
        #path_table = numpy.zeros((len(fly_seq), len(human_seq)))
        #base case already all zeros
        #for traceback 1 is up, 2 is from left, 3 is diagonal
        max_val = 0
        gap_val = -2
        mis_val = -1
        match_val = 2
        for i in range(1,len(fly_seq)):
            for j in range(1,len(human_seq)):
                val_choices = {}
                val_choices[1] = align_table[i,j-1] + gap_val
                val_choices[2] = align_table[i-1,j] + gap_val
                if fly_seq[i] == human_seq[j]:
                    val_choices[3] = align_table[i-1,j-1] + match_val
                elif fly_seq[i] != human_seq[j]:
                    val_choices[3] = align_table[i-1,j-1] + mis_val
                align_table[i,j] = max(val_choices.values())
                #path_table[i,j] = val_choices.keys()[val_choices.values().index(max(val_choices.values()))]

        max_last_row = max(align_table[-1,])
        max_last_col = max(align_table[-1])
        max_val = max(max_last_row, max_last_col)
        print "Optimal alignment score is:", max_val
