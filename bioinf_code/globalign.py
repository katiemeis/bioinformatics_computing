import os
import sys
import getopt
import numpy

if __name__ == "__main__":
    filename = ''
    filename2 = ''
    outfilename = ''
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "i:f:o:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-i":
            filename = arg
        elif opt == "-f":
            filename2 = arg
        elif opt == "-o":
            outfilename = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)

    with open(filename, 'r') as handle, open(filename2, 'r') as handle2, open(outfilename, 'w') as outfile:
        #read in data without id line, append lines to get seq
        human_lines = handle.readlines()[1:]
        human_lines_list = []
        human_lines_list.append('-')
        fly_lines = handle2.readlines()[1:]
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
        path_table = numpy.zeros((len(fly_seq), len(human_seq)))
        #base case
        for i in range(1,len(fly_seq)):
            align_table[i, 0] = i*(-2)
        for j in range(1,len(human_seq)):
            align_table[0, j] = j*(-2)
        #for traceback 1 is up, 2 is from left, 3 is diagonal
        max_val = 0
        gap_val = -2
        mis_val = -1
        match_val = 2
        for i in range(1,len(fly_seq)):
            for j in range(1,len(human_seq)):
                val_choices = {}
                val_choices[1] = align_table[i-1,j] + gap_val #same column, up one row
                val_choices[2] = align_table[i,j-1] + gap_val #same row, left one column
                if fly_seq[i] == human_seq[j]:
                    val_choices[3] = align_table[i-1,j-1] + match_val
                elif fly_seq[i] != human_seq[j]:
                    val_choices[3] = align_table[i-1,j-1] + mis_val
                align_table[i,j] = max(val_choices.values())
                #choose key with largest value of dictionary
                #if tie, will always choose the first in the dictionary - same order each time
                path_table[i,j] = val_choices.keys()[val_choices.values().index(max(val_choices.values()))]

        print "Optimal alignment score is:", align_table[-1, -1]

        alignment1 = []
        alignment2 = []
        #start from bottom right corner
        i = len(fly_seq)-1
        j = len(human_seq)-1
        #if value in i,j cell is 1, go back to cell above, if 2 go to left, if 3 go to diagonal
        #continue until get to zero
        while i > 0 and j > 0:
            if path_table[i,j] == 1:
                alignment1.append(fly_seq[i])
                alignment2.append('-')
                i -= 1
            elif path_table[i,j] == 2:
                alignment1.append('-')
                alignment2.append(human_seq[j])
                j -= 1
            elif path_table[i,j] == 3:
                alignment1.append(fly_seq[i])
                alignment2.append(human_seq[j])
                i -= 1
                j -= 1

        reverse_align1 = alignment1[::-1]
        reverse_align2 = alignment2[::-1]
        reverse_al1_seq = ''.join(reverse_align1)
        reverse_al2_seq = ''.join(reverse_align2)
        sep_list = []
        for i in range(0,len(reverse_al1_seq)):
            if reverse_al1_seq[i] == reverse_al2_seq[i]:
                sep_list.append('|')
            else:
                sep_list.append(' ')
        sep_seq = ''.join(sep_list)

        outfile.write(">fly (top) vs human (bottom)")
        outfile.write('\n')
        for i in range(0, len(reverse_al1_seq), 70):
            outfile.write(reverse_al1_seq[i:i+70])
            outfile.write('\n')
            outfile.write(sep_seq[i:i+70])
            outfile.write('\n')
            outfile.write(reverse_al2_seq[i:i+70])
            outfile.write('\n')
            outfile.write('\n')
