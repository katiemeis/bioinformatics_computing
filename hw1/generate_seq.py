import os
import sys
import getopt
import numpy
import random
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

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
        human_lines = handle.readlines()[1:]
        human_lines_list = []
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
        human_seq = ''.join(human_lines_list)
        #dictionary base, probability
        prob_dict = {"A":A_count/total_count, "T":T_count/total_count, "C":C_count/total_count, "G":G_count/total_count}
        #make dictionary for pairs
        #ignore commented section
        '''di_dict = {}
        for i in range(len(human_seq)-1):
            di = human_seq[i:i+2]
            if di not in di_dict:
                di_dict[di] = 1
            elif di in di_dict:
                di_dict[di] +=1
        di_freq_dict = {}
        di_total = sum(di_dict.values())
        for key, value in di_dict.iteritems():
            di_freq_dict[key] = di_dict[key]/float(di_total)
        #make dictionary for trios
        tri_dict = {}
        for i in range(len(human_seq)-2):
            tri = human_seq[i:i+3]
            if tri not in tri_dict:
                tri_dict[tri] = 1
            elif tri in tri_dict:
                tri_dict[tri] +=1
        tri_freq_dict = {}
        tri_total = sum(tri_dict.values())
        for key, value in tri_dict.iteritems():
            tri_freq_dict[key] = tri_dict[key]/float(tri_total)'''
        #dictionary for fourmers
        fourmer_dict = {}
        for i in range(len(human_seq)-3):
            fourmer = human_seq[i:i+4]
            if fourmer not in fourmer_dict:
                fourmer_dict[fourmer] = 1
            elif fourmer in fourmer_dict:
                fourmer_dict[fourmer] +=1
        fourmer_freq_dict = {}
        fourmer_total = sum(fourmer_dict.values())
        for key, value in fourmer_dict.iteritems():
            fourmer_freq_dict[key] = fourmer_dict[key]/float(fourmer_total)

        #make sequence
        new_seq_list = []
        #add first base based on individual frequency
        #correction to code, Scott wants all three in beginning added this way instead of with first and second order markov
        for i in range(3):
            #choose a random number between 0 and 1
            rand_num1 = random.uniform(0,1)
            #order dictionary by increasing value (frequency)
            ranked_base1 = sorted(prob_dict.items(), key=lambda prob_dict:prob_dict[1])
            #choose a base based on their frequencies and append to sequence list
            if rand_num1 < ranked_base1[0][1]:
                new_seq_list.append(ranked_base1[0][0])
            elif ranked_base1[0][1] < rand_num1 < ranked_base1[0][1]+ranked_base1[1][1]:
                new_seq_list.append(ranked_base1[1][0])
            elif ranked_base1[0][1]+ranked_base1[1][1] < rand_num1 < ranked_base1[0][1]+ranked_base1[1][1]+ranked_base1[2][1]:
                new_seq_list.append(ranked_base1[2][0])
            elif ranked_base1[0][1]+ranked_base1[1][1]+ranked_base1[2][1] < rand_num1 < 1.0:
                new_seq_list.append(ranked_base1[3][0])

        #ignore commented section
        #do first order to get second base in sequence
        '''rand_num2 = random.uniform(0,1)
        #grab the previous base for conditional
        prev_base = new_seq_list[0]
        #set of pair options
        new_di_opts = [prev_base+"A", prev_base+"T", prev_base+"C", prev_base+"G"]
        new_di_count = {}
        total_di_count = 0
        #calculate counts and totals for options in human sequence, if doesn't appear in human add with freq 0
        for option in new_di_opts:
            if option not in di_dict:
                new_di_count[option] = 0
            elif option in di_dict:
                new_di_count[option] = di_dict[option]
                total_di_count += di_dict[option]
        new_di_prob = {}
        #use above to calculate frequencies that add up to 1
        for key, value in new_di_count.iteritems():
            new_di_prob[key] = new_di_count[key]/float(total_di_count)
        #need to do last two loops because need freq based on previous, not out of total
        #order dictionary by increasing value (frequency)
        ranked_base2 = sorted(new_di_prob.items(), key=lambda new_di_prob:new_di_prob[1])
        #choose a base based on their frequencies and append to sequence list
        if rand_num2 < ranked_base2[0][1]:
            new_seq_list.append(ranked_base2[0][0][-1])
        elif ranked_base2[0][1] < rand_num2 < ranked_base2[0][1]+ranked_base2[1][1]:
            new_seq_list.append(ranked_base2[1][0][-1])
        elif ranked_base2[0][1]+ranked_base2[1][1] < rand_num2 < ranked_base2[0][1]+ranked_base2[1][1]+ranked_base2[2][1]:
            new_seq_list.append(ranked_base2[2][0][-1])
        elif ranked_base2[0][1]+ranked_base2[1][1]+ranked_base2[2][1] < rand_num2 < 1.0:
            new_seq_list.append(ranked_base2[3][0][-1])
        #do second order to get third base in sequence
        #same steps as above but conditioned on the two previous bases
        rand_num3 = random.uniform(0,1)
        prev_base3 = ''.join(new_seq_list[0:2])
        new_tri_opts = [prev_base3+"A", prev_base3+"T", prev_base3+"C", prev_base3+"G"]
        new_tri_count = {}
        total_tri_count = 0
        for option in new_tri_opts:
            if option not in tri_dict:
                new_tri_count[option] = 0
            elif option in tri_dict:
                new_tri_count[option] = tri_dict[option]
                total_tri_count += tri_dict[option]
        new_tri_prob = {}
        for key, value in new_tri_count.iteritems():
            new_tri_prob[key] = new_tri_count[key]/float(total_tri_count)
        ranked_base3 = sorted(new_tri_prob.items(), key=lambda new_tri_prob:new_tri_prob[1])
        if rand_num3 < ranked_base3[0][1]:
            new_seq_list.append(ranked_base3[0][0][-1])
        elif ranked_base3[0][1] < rand_num3 < ranked_base3[0][1]+ranked_base3[1][1]:
            new_seq_list.append(ranked_base3[1][0][-1])
        elif ranked_base3[0][1]+ranked_base3[1][1] < rand_num3 < ranked_base3[0][1]+ranked_base3[1][1]+ranked_base3[2][1]:
            new_seq_list.append(ranked_base3[2][0][-1])
        elif ranked_base3[0][1]+ranked_base3[1][1]+ranked_base3[2][1] < rand_num3 < 1.0:
            new_seq_list.append(ranked_base3[3][0][-1])'''

        #do the rest of the bases with in 3rd order
        #same steps as above but conditional on previous three bases
        for i in range(3,20000):
            #choose random number
            rand_num = random.uniform(0,1)
            #get previous three bases
            tri_base = ''.join(new_seq_list[i-3:i])
            #make strings of new options
            new_base_options = [tri_base+"A", tri_base+"T", tri_base+"C", tri_base+"G"]
            new_base_count = {}
            total_count_four = 0
            #calculate counts and totals for options in human sequence, if doesn't appear in human add with freq 0
            for option in new_base_options:
                if option not in fourmer_dict:
                    new_base_count[option] = 0
                elif option in fourmer_dict:
                    new_base_count[option] = fourmer_dict[option]
                    total_count_four += fourmer_dict[option]
            new_base_prob = {}
            #get frequencies from above loop counts
            for key, value in new_base_count.iteritems():
                new_base_prob[key] = new_base_count[key]/float(total_count_four)
            #order by increasing frequency
            ranked_base = sorted(new_base_prob.items(), key=lambda new_base_prob:new_base_prob[1])
            #choose new base based on how it compares to random number
            if rand_num < ranked_base[0][1]:
                new_seq_list.append(ranked_base[0][0][-1])
            elif ranked_base[0][1] < rand_num < ranked_base[0][1]+ranked_base[1][1]:
                new_seq_list.append(ranked_base[1][0][-1])
            elif ranked_base[0][1]+ranked_base[1][1] < rand_num < ranked_base[0][1]+ranked_base[1][1]+ranked_base[2][1]:
                new_seq_list.append(ranked_base[2][0][-1])
            elif ranked_base[0][1]+ranked_base[1][1]+ranked_base[2][1] < rand_num < 1.0:
                new_seq_list.append(ranked_base[3][0][-1])
        # join bases in list and print to file
        gened_sequence = ''.join(new_seq_list)
        outfile.write(gened_sequence)
        #can use following two lines instead of line 187 if want fasta format - also uncomment imports lines 6-8
        #seq = SeqRecord(Seq(gened_sequence), id = "generated sequence by Markov model", description = '')
        #SeqIO.write(seq, outfile, "fasta")
