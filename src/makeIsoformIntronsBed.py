#!/usr/bin/env python

import sys

def writeIsoformIntrons(isoforms_bed, output_bed):
    ip = open(isoforms_bed, 'r')
    op = open(output_bed, 'w')
    
    for line in ip:
        fields = line.strip().split("\t")
        start = int(fields[1])
        strand = fields[5]
        num_exons = int(fields[9])
        block_sizes = map(int, fields[10].split(','))
        block_offsets = map(int, fields[11].split(','))

        intron_num = 1 if (strand=='+') else num_exons-1
        for i in xrange(num_exons-1):
            intron_start = start + block_offsets[i] + block_sizes[i]
            intron_stop = start + block_offsets[i+1]
            
            try:
                assert (intron_start < intron_stop + 4), "Illegal intron for intron %d of %s" % (intron_num, fields[3])
            except AssertionError, ae:
                # Known problem for CG_67465.4.2.4, CG_79973.6.18.13, CG_88860.1.1.1
                print >> sys.stderr, "ERROR: %s" % ae.message
                continue
            
            op.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (fields[0], intron_start, intron_stop, fields[3], intron_num, strand))
            if (strand == '+'):
                intron_num += 1
            else:
                intron_num -= 1

    ip.close()
    op.close()


if (__name__ == "__main__"):
    isoforms_bed, output_bed = sys.argv[1:]
    writeIsoformIntrons(isoforms_bed, output_bed)
    sys.exit(0)
