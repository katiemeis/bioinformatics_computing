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

    with open(filename, 'r') as handle, open(outfilename, 'w') as outfile:
        line = handle.readlines()
        #get string without return, add - for t=0 spot
        seq = line[0].replace("\n", "")
        roll_seq = "-" + seq
        #make empty matrices and fill in t=0, fair is first row loaded is second row
        #if came from fair its 0, if came from loaded its 1 for traceback
        prob_table = numpy.zeros((2, len(roll_seq)))
        path_table = numpy.zeros((2, len(roll_seq)))
        prob_table[0,0] = 1
        #probabilities
        p_fair_fair, p_fair_loaded, p_loaded_loaded, p_loaded_fair = 0.95, 0.05, 0.9, 0.1
        pf_any_fair = 1/float(6)
        #dictionary for loaded
        loaded_probs = {1:0.1, 2:0.1, 3:0.1, 4:0.1, 5:0.1, 6:0.5}

        for j in range(1,len(roll_seq)):
            for i in range(2):
                #if in jth column of true row (first row)
                if i == 0:
                    #dictionary for choices
                    cell_choices = {}
                    #add values to choices, key is state it comes from (0 fair, 1 loaded)
                    cell_choices[0] = prob_table[i, j-1] * p_fair_fair * pf_any_fair
                    cell_choices[1] = prob_table[i+1, j-1] * p_loaded_fair * pf_any_fair
                    #add max of values to prob table
                    prob_table[i,j] = max(cell_choices.values())
                    #add the key of the max value to the path table
                    path_table[i,j] = cell_choices.keys()[cell_choices.values().index(max(cell_choices.values()))]
                #do the same steps as above but for loaded cell, jth column of 2nd row
                elif i == 1:
                    cell_choices = {}
                    cell_choices[0] = prob_table[i-1, j-1] * p_fair_loaded * loaded_probs[int(roll_seq[j])]
                    cell_choices[1] = prob_table[i, j-1] * p_loaded_loaded * loaded_probs[int(roll_seq[j])]
                    prob_table[i,j] = max(cell_choices.values())
                    path_table[i,j] = cell_choices.keys()[cell_choices.values().index(max(cell_choices.values()))]

        #dictionary for F/L and last column vals
        ending_vals = {0:prob_table[0,-1], 1:prob_table[1,-1]}
        state_seq = []
        #start cell of largest value of final column
        i = ending_vals.keys()[ending_vals.values().index(max(ending_vals.values()))] #which row has max value in last column
        j = len(roll_seq)-1

        while j > 0:
            if path_table[i,j] == 0:
                state_seq.append("F")
                j -= 1
                if i == 1:
                    i -= 1

            elif path_table[i,j] == 1:
                state_seq.append("L")
                j -= 1
                if i == 0:
                    i += 1

        reverse_seq = state_seq[::-1]
        rev_seq_str = ''.join(reverse_seq)
        outfile.write(rev_seq_str)
