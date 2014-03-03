#!/usr/bin/python

def enum(collection,st):
    #just like "enumerate" but you define your own starting position.
    #this returns indices RELATIVE TO ORIGINAL LIST
    i = st
    while i < len(collection):
        yield (i,collection[i])
        i += 1

def getfirstindex(L,st,value,K):
    for pos,t in enum(L,st):
        if t[1] > value-K:
            return pos #returns first read in the list that CAN contain given bp
    return 0

def getlastindex(L,st,value):
    for pos,t in enum(L,st):
        if t[1] > value:
            return pos #returns first read in the list that CANNOT contain given bp
    return len(L)

# ^these 2 were inspired by code here: http://stackoverflow.com/questions/946860/using-pythons-list-index-method-on-a-list-of-tuples-or-objects (see answer labeled "10", superperformant)

# function for converting SAM "flags" to bits:
def to_bin(n):
    return bin(n)[2:].zfill(11)

def countReadlets(fname, outfname, k, chromosome, stranded, rev):
    import pysam
    #from datetime import datetime #for debugging
    samfile = pysam.Samfile(fname,"rb")
    id_start_end = []
    cigar = []
    maxpos = 0
    minpos = 3000000000
    for read in samfile.fetch(chromosome):
        if to_bin(read.flag)[6] == '1' and not rev:
            strand = "-"
        elif to_bin(read.flag)[6] == '1' and rev:
            strand = "+"
        elif to_bin(read.flag)[6] == '0' and not rev:
            strand = "+"
        else:
            strand = "-"
        readstart = read.pos+1
        readend = read.aend
        if len(read.cigar)>1:
            currentloc = readstart
            for a in range(len(read.cigar)):
                if read.cigar[a][0] == 0: #match: count this segment as covered, obviously
                    id_start_end.append([read.qname, currentloc, currentloc+read.cigar[a][1]-1, strand])
                    currentloc = currentloc + read.cigar[a][1]
                if read.cigar[a][0] == 1: #insertion - not in reference, so our location doesn't move
                    continue
                if read.cigar[a][0] == 2: #deletion - count this segment as covered, since a read overlaps it.
                    id_start_end.append([read.qname, currentloc, currentloc+read.cigar[a][1]-1, strand])
                    currentloc = currentloc + read.cigar[a][1]
                if read.cigar[a][0] == 3: #skipped region: don't count and move the location
                    currentloc = currentloc + read.cigar[a][1]
                if read.cigar[a][0] > 3: #I don't think the rest of these flags exist in these alignments, but just make sure.
                    print("WARNING: you have a cigar flag in here that you haven't accounted for.")
        else:
            id_start_end.append([read.qname, readstart, readend, strand])
        if read.aend > maxpos:
            maxpos = read.aend
        if read.pos+1 < minpos:
            minpos = read.pos+1
    # sort list by starting points, since we've now introduced extra complication...:
    import operator
    id_start_end.sort(key=operator.itemgetter(1))
    # credit: http://stackoverflow.com/questions/5201191/sorting-a-list-of-lists-in-python

    if stranded:
        fplus = open(outfname+'_plus', 'w')
        fminus = open(outfname+'_minus', 'w')
        first = 0
        last = 0
        for z in xrange(1,minpos):
            fplus.write("%s\t%s\n" % (z,0))
            fminus.write("%s\t%s\n" % (z,0))
        for i in xrange(minpos, maxpos+1):
            last = getlastindex(id_start_end,first,i)
            first = getfirstindex(id_start_end,first,i,k)
            if first == last:
                fplus.write("%s\t%s\n" % (i,0))
                fminus.write("%s\t%s\n" % (i,0))
                continue
            overlaps = id_start_end[first:last]
            readnames_plus = set()
            readnames_minus = set()
            for j in xrange(len(overlaps)):
                if i>=overlaps[j][1] and i<=overlaps[j][2]:
                    if overlaps[j][3] == "+":
                        readnames_plus.add(overlaps[j][0])
                    else:
                        readnames_minus.add(overlaps[j][0])
            numreads_plus = len(readnames_plus) #count readlets only once per read
            numreads_minus = len(readnames_minus) #count readlets only once per read
            fplus.write("%s\t%s\n" % (i,numreads_plus))
            fminus.write("%s\t%s\n" % (i,numreads_minus))
        fplus.close()
        fminus.close()
        
    else:
        f = open(outfname, 'w')
        #g = open('keeptrackofiterations','w') #for debugging
        first = 0
        last = 0
        #npos = 0 #for debugging
        for z in xrange(1,minpos):
            f.write("%s\t%s\n" % (z,0))
        for i in xrange(minpos, maxpos+1):
        #for i in xrange(maxpos-2999,maxpos+1): #for debugging
            #npos = npos+1 #for debugging
            #if npos % 10000 == 0: g.write("Did 10000 "+str(datetime.now())+"\n") #for debugging
            last = getlastindex(id_start_end,first,i)
            first = getfirstindex(id_start_end,first,i,k)
            if first == last:
                f.write("%s\t%s\n" % (i,0))
                continue
            overlaps = id_start_end[first:last]
            readnames = set()
            for j in xrange(len(overlaps)):
                if i>=overlaps[j][1] and i<=overlaps[j][2]:
                    readnames.add(overlaps[j][0])
            numreads = len(readnames) #count readlets only once per read
            f.write("%s\t%s\n" % (i,numreads))
        f.close()

    return None

# get arguments from command line
from optparse import OptionParser
opts = OptionParser()
opts.add_option("--file","-f",type="string",help="input file name (must be .bam)")
opts.add_option("--output","-o",type="string",help="output file name")
opts.add_option("--kmer","-k",type="int",help="kmer length")
opts.add_option("--chrom","-c",type="string",help="chromosome to parse")
opts.add_option("--stranded","-s",type="string",help="stranded or reverse-stranded protocol?",default="FALSE")
options,arguments = opts.parse_args()

if options.stranded == "TRUE":
    stranded = True
    rev = False
elif options.stranded == "REVERSE":
    stranded = True
    rev = True
elif options.stranded != "FALSE":
    ValueError('stranded must either be TRUE, REVERSE, or FALSE')
else:
    stranded = False
    rev = False

countReadlets(options.file, options.output, options.kmer, options.chrom, stranded, rev)
