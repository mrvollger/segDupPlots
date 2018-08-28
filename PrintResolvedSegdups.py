#!/usr/bin/env python

import sys
import argparse
ap = argparse.ArgumentParser(description="Determine which segdups are assembled and which split contigs.")
ap.add_argument("tab", help="Table of contig/segdup intersection.  bedtools intersect -a contigs.bed -b segdups.bed -wa -wb")
ap.add_argument("segdups", help="Original segdups file.")
ap.add_argument("--extra", help="Required amount on the side of a segdup to count as resolved.", type=int, default=20000)
ap.add_argument("--resolved", help="Print resolved segdups here.", default=None)
ap.add_argument("--unresolved", help="Print unresolved segdups here.", default=None)
ap.add_argument("--allseg", help="Print all segdups here.", default=None)

args=ap.parse_args()
tabFile = open(args.tab)

#chr13   53057721        53093613        utg7180000000010/0_35898        254     chr13   53051402        53076388
#chr13   53057721        53093613        utg7180000000010/0_35898        254     chr13   53076392        53170864
#chr13   52867234        52912667        utg7180000000010/0_35898        0       chr13   52889195        52892882
#chr13   52867234        52912667        utg7180000000010/0_35898        0       chr13   52885495        52889187
#chr13   52867234        52912667        utg7180000000010/0_35898        0       chr13   52892898        52918144
#chr13   52867234        52912667        utg7180000000010/0_35898        0       chr13   52771866        52884436
if (args.resolved is None):
	args.resolved = args.tab + ".resolved"
if (args.unresolved is None):
	args.unresolved = args.tab + ".unresolved"

segdups = {}
segdupFile = open(args.segdups)
segdupLines = segdupFile.readlines()
for line in segdupLines:
    vals = line.split()
    segdups["_".join(vals[0:3])] = (False, -1)
    
r = open(args.resolved, 'w')
u = open(args.unresolved, 'w')
a = open(args.allseg, 'w')


for line in tabFile:
	vals = line.split()
	mapq = int(vals[4])
	if (mapq < 30):
		continue
	(ctgChrom, ctgStart, ctgEnd, ctgName) = (vals[0], int(vals[1]), int(vals[2]), vals[3])
	(segChrom, segStart, segEnd) = (vals[5], int(vals[6]), int(vals[7]))
	pref = segStart - ctgStart
	suff = ctgEnd - segEnd
	segdup = "_".join(vals[5:8])
	if (pref > args.extra and suff > args.extra):
		segdups[segdup] = (True, min(pref, suff))
	else:
		segdups[segdup] = (False, min(pref, suff))
	#		r.write("\t".join(vals[5:8]) + "\n")
	#	else:
	#		u.write("\t".join(vals[5:8]) + "\n")

for line in segdupLines:
    vals = line.split()
    segdup = "_".join(vals[0:3])
    if (segdups[segdup][0] == True):
        r.write(line)
    else:
        u.write(line)
    # wirte all entreis
    vals.append( str(segdups[segdup][1]) )
    line = "\t".join(vals)  + "\n"
    a.write(line)

r.close()
u.close()
a.close()

		

