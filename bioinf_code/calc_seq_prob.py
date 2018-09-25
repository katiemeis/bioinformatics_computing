import os
import sys
import getopt
import numpy
import math

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
            filename2 = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)
    with open(filename, 'r') as handle, open(filename2, 'r') as handle2:

        #read in data without id line
        human_lines = handle.readlines()[1:]
        neand_lines = handle2.readlines()[1:]
        human_lines_list = []
        neand_lines_list = []
        #multinomial model
        #make counts 0
        A_count, T_count, C_count, G_count = 0, 0, 0, 0
        #get counts for probabilities
        for line in human_lines:
            new_line = line.strip('\n')
            new_line2 = new_line.replace("N", "A")
            human_lines_list.append(new_line2)
            for element in new_line:
                if element == "A":
                    A_count += 1
                elif element == "T":
                    T_count += 1
                elif element == "C":
                    C_count += 1
                elif element == "G":
                    G_count += 1

        total_count = float(A_count + T_count + C_count + G_count)
        #dictionary base, probability
        prob_dict = {"A":A_count/total_count, "T":T_count/total_count, "C":C_count/total_count, "G":G_count/total_count}
        multinom_log_prob = 0
        #add up log probabilities for neader sequence based on human sequence
        for line in neand_lines:
            new_line = line.strip('\n')
            neand_lines_list.append(new_line)
            for element in new_line:
                if element == "A" or element == "T" or element == "C" or element == "G":
                    multinom_log_prob += math.log(prob_dict[element])

        print "Log probability (multinomial, ln):", multinom_log_prob

        #markov model
        human_seq = ''.join(human_lines_list)
        fourmer_dict = {}
        #fill in dictionary with fourmer, count
        for i in range(len(human_seq)-3):
            fourmer = human_seq[i:i+4]
            if fourmer not in fourmer_dict:
                fourmer_dict[fourmer] = 1
            elif fourmer in fourmer_dict:
                fourmer_dict[fourmer] +=1

        neand_seq = ''.join(neand_lines_list)
        markov_log_prob = 0
        #look up fourmers in neander seq, calulate their probs, add log prob to seq prob
        for i in range(len(neand_seq)-3):
            nfourmer = neand_seq[i:i+4]
            first3 = nfourmer[0:3]
            freq_total = 0
            for key, value in fourmer_dict.iteritems():
                if first3 == key[0:3]:
                    freq_total += value
            four_prob = fourmer_dict[nfourmer] / float(freq_total)
            markov_log_prob += math.log(four_prob)
        #since loop misses the first, second and third lets add those
        #correction Scott wants probs based on multinomial for all thre in beginning instead of markov as done below
        markov_log_prob += math.log(prob_dict[neand_seq[0]])
        markov_log_prob += math.log(prob_dict[neand_seq[1]])
        markov_log_prob += math.log(prob_dict[neand_seq[2]])
        '''#di
        di_dict = {}
        di_total = 0
        for i in range(len(human_seq)-1):
            di = human_seq[i:i+2]
            if di not in di_dict:
                di_dict[di] = 1
            elif di in di_dict:
                di_dict[di] +=1
        for key, value in di_dict.iteritems():
            if neand_seq[0] == key[0]:
                di_total += value
        di_freq = di_dict[neand_seq[0:2]]/float(di_total)
        markov_log_prob += math.log(di_freq)
        seq_prob_count +=1
        #tri
        tri_dict = {}
        tri_total = 0
        for i in range(len(human_seq)-2):
            tri = human_seq[i:i+3]
            if tri not in tri_dict:
                tri_dict[tri] = 1
            elif tri in tri_dict:
                tri_dict[tri] +=1
        for key, value in tri_dict.iteritems():
            if neand_seq[0:2] == key[0:2]:
                tri_total += value
        tri_freq = tri_dict[neand_seq[0:3]]/float(tri_total)
        markov_log_prob += math.log(tri_freq)
        seq_prob_count +=1'''

        print "Log probability (Markov, ln):", markov_log_prob
