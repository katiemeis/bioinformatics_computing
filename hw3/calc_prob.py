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
        line = handle.readlines()
        #get string without return, add - for t=0 spot
        seq = line[0].replace("\n", "")
        roll_seq = "-" + seq
        #make empty matrices and fill in t=0, fair is first row loaded is second row
        prob_table = numpy.zeros((2, len(roll_seq)))
        prob_table[0,0] = 1
        #probabilities
        p_fair_fair, p_fair_loaded, p_loaded_loaded, p_loaded_fair = 0.95, 0.05, 0.9, 0.1
        pf_any_fair = 1/float(6)
        #dictionary for loaded
        loaded_probs = {1:0.1, 2:0.1, 3:0.1, 4:0.1, 5:0.1, 6:0.5}

        for j in range(1,len(roll_seq)):
            for i in range(2):
                if i == 0:
                    cell_sum = (prob_table[i, j-1] * p_fair_fair * pf_any_fair) + (prob_table[i+1, j-1] * p_loaded_fair * pf_any_fair)
                    prob_table[i,j] = cell_sum
                elif i == 1:
                    cell_sum = (prob_table[i, j-1] * p_loaded_loaded * loaded_probs[int(roll_seq[j])]) + (prob_table[i-1, j-1] * p_fair_loaded * loaded_probs[int(roll_seq[j])])
                    prob_table[i,j] = cell_sum

        print prob_table[0,-1] + prob_table[1,-1]
