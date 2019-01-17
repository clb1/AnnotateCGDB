#!/usr/bin/env python
# Output of positions follows BED format (as opposed to GTF), where start is zero-based and end is 1-based.

import sys

ip = sys.stdin
op = sys.stdout

line = ip.readline()
curr_exonic_region = None # Format: [chr, region_start, region_end, transcript_id, is_unique_region]
last_genome_pos = None
while (line != ""):
    chrom, chromStart, chromEnd, transcript_id, index, curr_coverage = line.split("\t")
    curr_genome_pos = int(chromStart) + int(index) - 1
    curr_coverage = int(curr_coverage)

    #is_unique_region = 1 if (curr_coverage==1) else 0
    #if (curr_exonic_region == None):
    #    curr_exonic_region = [chrom, curr_genome_pos, None, transcript_id, is_unique_region]
    #elif ((last_genome_pos != curr_genome_pos-1) or (is_unique_region != curr_exonic_region[-1]) or (transcript_id != curr_exonic_region[3])):
    #    curr_exonic_region[2] = last_genome_pos+1
    #   op.write("%s\t%d\t%d\t%s\t%d\n" % tuple(curr_exonic_region))
    #    curr_exonic_region = [chrom, curr_genome_pos, None, transcript_id, is_unique_region]

    if (curr_exonic_region == None):
        curr_exonic_region = [chrom, curr_genome_pos, None, transcript_id, curr_coverage]
    elif ((last_genome_pos != curr_genome_pos-1) or (curr_coverage != curr_exonic_region[-1]) or (transcript_id != curr_exonic_region[3])):
        curr_exonic_region[2] = last_genome_pos+1
        op.write("%s\t%d\t%d\t%s\t%d\n" % tuple(curr_exonic_region))
        curr_exonic_region = [chrom, curr_genome_pos, None, transcript_id, curr_coverage]

    last_genome_pos = curr_genome_pos
    line = ip.readline()

# Finish the last region
curr_exonic_region[2] = last_genome_pos+1
op.write("%s\t%d\t%d\t%s\t%d\n" % tuple(curr_exonic_region))
    
sys.exit(0)
