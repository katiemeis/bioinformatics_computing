import os
import sys
import getopt
import random

if __name__ == "__main__":
    outfilename = ''
    try:
        (opts, args) = getopt.getopt(sys.argv[1:], "o:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-o":
            outfilename = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)


    with open(outfilename, 'w') as outfile:
        #define probabilities
        p_fair_fair = 0.95
        p_fair_loaded = 0.05
        p_loaded_loaded = 0.9
        p_loaded_fair = 0.1

        pf1, pf2, pf3, pf4, pf5, pf6 = 1/float(6), 1/float(6), 1/float(6), 1/float(6), 1/float(6), 1/float(6)
        pl1, pl2, pl3, pl4, pl5, pl6 = 0.1, 0.1, 0.1, 0.1, 0.1, 0.5

        #empty list to add rolls to, add first roll knoing it has to be fair
        rolls = []
        fair_choices = [1,2,3,4,5,6]
        rolls.append(random.choice(fair_choices))
        loaded_choices = [1,2,3,4,5]
        #empty list for states, fair 0 loaded 1, add first state knowing it's fair
        states = []
        states.append(0)
        #loop through the rest adding to both lists
        for i in range(1,300):
            #generate random number
            rand_num = random.uniform(0,1)
            #if the previous state was fair
            if states[i-1] == 0:
                #switch to loaded
                if rand_num < p_fair_loaded:
                    states.append(1)
                    rand_num_loaded = random.uniform(0,1)
                    if rand_num_loaded < 0.5:
                        rolls.append(6)
                    elif rand_num_loaded > 0.5:
                        rolls.append(random.choice(loaded_choices))
                #stay fair
                elif rand_num > p_fair_loaded:
                    states.append(0)
                    rolls.append(random.choice(fair_choices))
            #if the previous state was loaded
            elif states[i-1] == 1:
                #switch to fair
                if rand_num < p_loaded_fair:
                    states.append(0)
                    rolls.append(random.choice(fair_choices))
                #stay loaded
                elif rand_num > p_loaded_fair:
                    states.append(1)
                    rand_num_loaded = random.uniform(0,1)
                    if rand_num_loaded < 0.5:
                        rolls.append(6)
                    elif rand_num_loaded > 0.5:
                        rolls.append(random.choice(loaded_choices))
        print rolls
        print len(rolls)
        for item in rolls:
            outfile.write("%s" % item)
