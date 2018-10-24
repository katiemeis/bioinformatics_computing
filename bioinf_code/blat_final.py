""" PWM related functions. """

import os
import sys
import getopt
import numpy as np
from collections import deque
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class axtFile:
    " class to handle axt file produced by Blastz"
    def __init__(self, filename):
        """ reads axt file and builds a dictionary of start/end positions
          and sequence matches """

        try:
            filehandle = open(filename, "r")
            try:
                stream = filehandle.read()
            finally:
                filehandle.close()
        except IOError:
            print "File %s doesn't exist\n" % filename
            sys.exit(2)
        
	records = stream.split("\n\n")
	self.records = records
	
        self.matches = []

        "do the first record separately"
        firstrec = records[0]
        lines = firstrec.splitlines()
        fields = lines[0].split()
        bs = int(fields[2])

        seq1 = lines[1].upper()
        seq2 = lines[2].upper()

	
        "Add matches within the aligned segment"
        gaps_removed = []
        
        "Remove gaps in the first sequence"
        for alignpair in zip(seq1,seq2):
            if alignpair[0] == "-":
                continue
            gaps_removed.append(alignpair)    

        self.matches.extend([x == y for (x,y) in gaps_removed])

        bend = int(fields[3])
        

        self.startpos = bs

        "do the rest of the records, adding false for intergap regions"
        for record in records[1:]:
            lines = record.splitlines()
            fields = lines[0].split()
            bs = int(fields[2])

            seq1 = lines[1].upper()
            seq2 = lines[2].upper()

            "Add false for gaps between aligned segments"
            self.matches.extend([False for x in range(bs-bend-1)])

            "Add matches within the aligned segment"
            gaps_removed = []

            "Remove gaps in the first sequence"
            for alignpair in zip(seq1,seq2):
                if alignpair[0] == "-":
                    continue
                gaps_removed.append(alignpair)    

            self.matches.extend([x == y for (x,y) in gaps_removed])

            bend = int(fields[3])
        
        self.endpos = bend

        print "Finished loading data...\n"

	self.computeMovingAverage(width)
	self.makeOutlist(threshold)
        

    def computeMovingAverage(self, width):
        " Compute moving average of sequence matches for a window width "

        "because the window includes the central nucleotide at center"
        self.width = width
        wd = width - 1
        "faster container datatype"
        m = deque(self.matches)

        matches = np.array(m, dtype=float)
        nummatches = len(matches)
        movavg = np.array(m, dtype=float)

        "initialize movavg"
        for j in range(width):
            movavg[j] = 0.0
            movavg[-(j+1)] = 0.0

        """ sum the first window. The center nucleotide is at width+1 
        nucleotides. We sum the matches width nucleotides before and after.
        """
        sum = matches[width]
        for j in range(width):
            sum = sum + matches[j] + matches[width+1+j]

        "This is the consensus sum at position width"
        movavg[width] = sum

        for j in range(width+1, nummatches-width):
            movavg[j] = movavg[j-1] - matches[j-width-1] + movavg[j+width]

        self.movavg = movavg/(2*width+1) 
        
        print "Finished computing averages...\n"


    def makeOutlist(self, threshold=0.0):
	" Makes a list of position and thresholded conservation score "

        t = threshold
        positions = range(self.startpos, self.endpos+1) 
	outlist = ["%d  %f" % (pos, sco > t and sco or 0.0) \
                            for  (pos, sco) in \
                            zip(positions, self.movavg)]
	self.outlist = outlist	


    def printWig(self, outfilename, skippos = 0):
        """ prints out the moving average in wig format """

	outlist = self.outlist
	
        if skippos:
            outlist = outlist[::skippos]

        outstr = "\n".join(outlist[self.width/2-1:-1*(self.width/2)])
        try:
            filehandle = open(outfilename, "w")
            try:
                stream = filehandle.write(outstr)
            finally:
                filehandle.close()
        except IOError:
            print "Error writing to file %s\n" % outfilename
            sys.exit(2)

        print "Finished writing wig file...\n"


    def printFasta(self, fastafile):
	"prints out crms in fasta format"
	
	outlist = self.outlist
	crms = []
	start = False
	sequences = []
	prevend = -1

	for item in outlist: 

		part = item.split()
		pos = int(part[0])
		sco = float(part[1])

		if sco >= threshold:
			if start == False:
				start = True
				crmstart = pos
			prevend = pos

		elif start:
				start = False
				crmend = prevend
				crms.append((crmstart, crmend))

	for crm in crms:

		crmstart = crm[0]
		crmend = crm[1]

		for record in self.records:

			lines = record.splitlines()
	                fields = lines[0].split()
			chrom = fields[1]
	                bs = int(fields[2])
			bend = int(fields[3])
			sequence = lines[1].replace('-', '')

			if bs < crmstart < bend:
				crmseq = sequence[crmstart-bs:crmend-bs+1]
				seq = SeqRecord(Seq(crmseq),
					id = str(chrom),
					description = str(crm))
				sequences.append(seq)
			
	SeqIO.write(sequences, fastafile, "fasta")

	print "Finished writing fasta file...\n"


    def blatSearch(self, fastafile, pslfile):
	"blat search for sequences in fasta file"

	try:
            fasta = open(fastafile, "r")
        except IOError:
            print "File %s doesn't exist\n" % filename
            sys.exit(2)

	psl = open(pslfile, "w")	

	os.system('blat /labcommon/Genome/Mouse/mm9/Mus_musculus.NCBIM37.64.dna.toplevel.fa {} {}'.format(fastafile, pslfile))

	try:
            blathandle = open(pslfile, "r")
	    try:
                stream = blathandle.read()
            finally:
                blathandle.close()
        except IOError:
            print "File %s doesn't exist\n" % filename
            sys.exit(2)

	line = stream.splitlines()

	for line in line[5:]:
		part = line.split("\t")
		print part[9]
		print (part[15], part[16])



if __name__ == "__main__":


    axtfilename = ''
    outfilename = ''
    fastafile = ''
    pslfile = ''
    width = 10
    skippos = 0
    threshold = 0.0

    try:    
        (opts, args) = getopt.getopt(sys.argv[1:], "o:s:t:w:f:b:")
    except getopt.GetoptError:
        print "Problem getting options\n"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-o":
            outfilename = arg
        elif opt == "-s":
            skippos = int(arg)
        elif opt == "-t":
            threshold = float(arg)
        elif opt == "-w":
            width = int(arg)
	elif opt == "-f":
	    fastafile = arg
	elif opt == "-b":
	    pslfile = arg
        else:
            print "Don't know %s \n" % opt
            sys.exit(2)

    if len(args) == 1:
        axtfilename = args[0]
    else:    
        print """Need one argument: the name of the file
                containing the axt Blastz alignment .\n"""
        sys.exit(2)

    if not outfilename:
        print """Specify an ouput filename with -o\n"""
        sys.exit(2)

    axt = axtFile(axtfilename)
    axt.printWig(outfilename, skippos)
    axt.printFasta(fastafile)
    axt.blatSearch(fastafile, pslfile)
