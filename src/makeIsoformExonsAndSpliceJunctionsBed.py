#!/usr/bin/env python
#
# NOTE: Exons in the input gtf file are expected to be specified in order (i.e. exon1, exon2, ...)
#

import gzip
import re
import sys

def getIsoformExons(input_gff):
    reGtfLine = re.compile("transcript_id .(\S+)..\s+")
    isoform_exons = {}

    if (input_gff.endswith(".gz")):
        ip = gzip.open(input_gff, 'rb')
    else:
        ip = open(input_gff, 'r')

    line = ip.readline()
    while (line != ""):
        if (line[0] != '#'):
            elems = line.split('\t')
            if (elems[2] == "exon"):
                mo = reGtfLine.search(line)
                transcript_id = mo.group(1)
                if (not isoform_exons.has_key(transcript_id)):
                    isoform_exons[transcript_id] = []
                tup = (elems[0], int(elems[3])-1, int(elems[4]), transcript_id, elems[6])
                isoform_exons[transcript_id].append( tup )
        line = ip.readline()
    ip.close()
    
    return isoform_exons

    
def writeExonsAndSpliceJunctions(isoform_exons, output_bed):
    op = open(output_bed, 'w')

    # Form the isoform exons and splice junction "parts"
    for transcript_id, transcript_exons in isoform_exons.items():
        transcript_exons.sort(lambda x,y:cmp(x[1],y[1]))
        num_exons = len(transcript_exons)
        for i in xrange(num_exons):
            tup = transcript_exons[i]
            strand = tup[-1]
            op.write("%s\t%d\t%d\t%s\t1\t%s\t%d\t%d\t-\t1\t%d\t0\n" % (tup[0], tup[1], tup[2], tup[3], strand, tup[1], tup[2], tup[2]-tup[1]))
            if (i < num_exons-1):
                region_start = tup[2]-1
                region_end = transcript_exons[i+1][1] + 1
                second_block_start = region_end - region_start - 1
                op.write("%s\t%d\t%d\t%s\t1\t%s\t%d\t%d\t-\t2\t1,1\t0,%d\n" % \
                         (tup[0], region_start, region_end, tup[3], strand, region_start, region_end, second_block_start))
                         
                #op.write("%s\t%d\t%d\t%s 0\t+\t0\t0\t-\t%2\t1,1\t%0,%d\n" % (tup[0], tup[2]-1, tup[2], tup[-1]))
                #next_exon_start = transcript_exons[i+1][1] + 1
                #op.write("%s\t%d\t%d\t%s 0\t+\t0\t0\t-\t%1\t%d\t%d\n" % (tup[0], next_exon_start, next_exon_start+1, tup[-1]))

    op.close()

    
if (__name__ == "__main__"):
    input_gtf, output_bed = sys.argv[1:]

    isoform_exons = getIsoformExons(input_gtf)
    writeExonsAndSpliceJunctions(isoform_exons, output_bed)

    sys.exit(0)
