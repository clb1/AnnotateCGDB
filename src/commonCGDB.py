from collections import defaultdict
import gzip
import re
import sys


def partsSort(A,B):
    ret_val = cmp(A[0][1], B[0][1])
    if (ret_val == 0):
        ret_val = cmp(A[0][2], B[0][2])
    return ret_val


def associateIsoformsToTargetGenes(input_gtf, target_isoforms):
    target_genes = set(map(lambda x: x.split('.')[0], target_isoforms))
    all_target_loci_isoforms = set()
    
    isoforms_per_gene = defaultdict(set)

    reLine = re.compile("gene_id .(\S+).. transcript_id .(\S+).. ") # .+ oId .(\S+)..
    if (input_gtf.endswith(".gz")):
        ip = gzip.open(input_gtf, 'rb')
    else:
        ip = open(input_gtf, 'r')
    for line in ip:
        if (line[0] != '#'):
            mo = reLine.search(line)
            try:
                gene, isoform = mo.groups()
            except AttributeError:
                print >> sys.stderr, "ERROR: processing gtf file to associate isoforms to genes. Line is:"
                print >> sys.stderr, line
                sys.exit(1)

            if (gene in target_genes):
                isoforms_per_gene[gene].add( isoform )
                all_target_loci_isoforms.add( isoform )
    ip.close()
    return isoforms_per_gene, all_target_loci_isoforms


def associateIsoformsToTargetGenesPlusOldID(input_gtf, target_isoforms):
    target_genes = set(map(lambda x: x.split('.')[0], target_isoforms))
    all_target_loci_isoforms = set()
    isoformOID = {}
    
    isoforms_per_gene = defaultdict(set)

    reLine = re.compile("gene_id .(\S+).. transcript_id .(\S+).. .*oID .(\S+)\|.+..")
    reLineSimpler = re.compile("gene_id .(\S+).. transcript_id .(\S+).. .*oID .(\S+)..$")
    if (input_gtf.endswith(".gz")):
        ip = gzip.open(input_gtf, 'rb')
    else:
        ip = open(input_gtf, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        if (line[0] != '#' and fields[2] == "exon"):
            mo = reLine.search(fields[-1])
            try:
                gene, isoform, oID = mo.groups()
            except AttributeError:
                mo = reLineSimpler.search(fields[-1])
                gene, isoform, oID = mo.groups()
            except AttributeError:
                print >> sys.stderr, "ERROR: processing gtf file to associate isoforms to genes. Line is:"
                print >> sys.stderr, line
                sys.exit(1)

            if (gene in target_genes):
                isoforms_per_gene[gene].add( isoform )
                all_target_loci_isoforms.add( isoform )
                isoformOID[isoform] = oID
    ip.close()

    return isoforms_per_gene, all_target_loci_isoforms, isoformOID


def associateIsoformsToGenes(input_gtf):
    isoforms_per_gene = defaultdict(set)
    reLine = re.compile("gene_id .(\S+).. .*transcript_id .(\S+).. ")
    if (input_gtf.endswith(".gz")):
        ip = gzip.open(input_gtf, 'rb')
    else:
        ip = open(input_gtf, 'r')
    for line in ip:
        fields = line.strip().split("\t")
        if (line[0] != '#' and fields[2] == "exon"):
            mo = reLine.search(fields[-1])
            try:
                gene, isoform = mo.groups()
            except AttributeError:
                print >> sys.stderr, "ERROR: processing gtf file to associate isoforms to genes. Line is:"
                print >> sys.stderr, line
                sys.exit(1)

            isoforms_per_gene[gene].add( isoform )
    ip.close()
    return isoforms_per_gene


def associateIsoformsToGenesPlusOldID(input_gtf):
    isoforms_per_gene = defaultdict(set)
    isoformOID = {}
    reLine = re.compile("gene_id .(\S+).. transcript_id .(\S+).. .*oID .(\S+)\|.+..")
    if (input_gtf.endswith(".gz")):
        ip = gzip.open(input_gtf, 'rb')
    else:
        ip = open(input_gtf, 'r')
    for line in ip:
        if (line[0] != '#'):
            mo = reLine.search(line)
            try:
                gene, isoform, oID = mo.groups()
            except AttributeError:
                print >> sys.stderr, "ERROR: processing gtf file to associate isoforms to genes. Line is:"
                print >> sys.stderr, line
                sys.exit(1)

            isoforms_per_gene[gene].add( isoform )
            isoformOID[isoform] = oID
            
    ip.close()
    return isoforms_per_gene, isoformOID


def buildTargetIsoformPartsList(input_bed, target_isoforms):
    """Coordinates of parts are in BED format, which is 0-based, start
    position is inclusive and end position is exclusive (ie 1 position
    passed end of elemement."""
    reBedLine = re.compile("(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+1\s+(.)\s+\d+\s+\d+\s+-\s+(\d)\s+\S+\s+\S+")
    isoform_parts_list = {}
    for isoform in target_isoforms:
        isoform_parts_list[isoform] = {}

    print >> sys.stderr, "Building parts list...",
    ip = open(input_bed, 'r')
    for line in ip:
        mo = reBedLine.match(line)
        chrom, chromStart, chromEnd, transcript_id, strand, num_blocks = mo.groups()
        num_blocks = int(num_blocks)
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
            
        if (transcript_id in target_isoforms):
            # Format: ["e"|"ue"|"se"|"sj", part number placeholder, is_supported_by_an_isoform-consistent_read]
            if (num_blocks == 1): # Is exon
                isoform_parts_list[transcript_id][(chrom, chromStart, chromEnd)] = ["e", -1, False] 
            elif (num_blocks == 2): # Is splice junction
                isoform_parts_list[transcript_id][(chrom, chromStart, chromEnd)] = ["sj", -1, False]
            else:
                print >> sys.stderr, "ERROR: number of blocks > 2.", line
                sys.exit(1)

            if (isoform_parts_list[transcript_id].has_key("strand")):
                assert (isoform_parts_list[transcript_id]["strand"] == strand)
            else:
                isoform_parts_list[transcript_id]["strand"] = strand
    ip.close()
    print >> sys.stderr, "done"
    
    return isoform_parts_list


def buildIsoformPartsList2(input_bed, isoform_parts_list):
    """Coordinates of parts are in BED format, which is 0-based, start
    position is inclusive and end position is exclusive (ie 1 position
    passed end of elemement."""
    reBedLine = re.compile("(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+1\s+(.)\s+\d+\s+\d+\s+-\s+(\d)\s+\S+\s+\S+")

    print >> sys.stderr, "Building parts list...",
    ip = open(input_bed, 'r')
    for line in ip:
        mo = reBedLine.match(line)
        chrom, chromStart, chromEnd, transcript_id, strand, num_blocks = mo.groups()
        num_blocks = int(num_blocks)
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
        if (not isoform_parts_list.has_key(transcript_id)):
            isoform_parts_list[transcript_id] = {}
            
        # Format: ["e"|"ue"|"se"|"sj", part number placeholder, is_supported_by_an_isoform-consistent_read]
        if (num_blocks == 1): # Is exon
            isoform_parts_list[transcript_id][(chrom, chromStart, chromEnd)] = ["e", -1, False] 
        elif (num_blocks == 2): # Is splice junction
            isoform_parts_list[transcript_id][(chrom, chromStart, chromEnd)] = ["sj", -1, False]
        else:
            print >> sys.stderr, "ERROR: number of blocks > 2.", line
            sys.exit(1)
            
        if (isoform_parts_list[transcript_id].has_key("strand")):
            assert (isoform_parts_list[transcript_id]["strand"] == strand)
        else:
            isoform_parts_list[transcript_id]["strand"] = strand

    ip.close()
    print >> sys.stderr, "done"


def buildIsoformPartsList(input_bed):
    """Coordinates of parts are in BED format, which is 0-based, start
    position is inclusive and end position is exclusive (ie 1 position
    passed end of elemement."""
    reBedLine = re.compile("(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+1\s+(.)\s+\d+\s+\d+\s+-\s+(\d)\s+\S+\s+\S+")
    isoform_parts_list = {}

    ip = open(input_bed, 'r')
    for line in ip:
        mo = reBedLine.match(line)
        chrom, chromStart, chromEnd, transcript_id, strand, num_blocks = mo.groups()
        num_blocks = int(num_blocks)
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
        if (not isoform_parts_list.has_key(transcript_id)):
            isoform_parts_list[transcript_id] = {}
            
        # Format: ["e"|"ue"|"se"|"sj", part number placeholder, is_supported_by_an_isoform-consistent_read]
        if (num_blocks == 1): # Is exon
            isoform_parts_list[transcript_id][(chrom, chromStart, chromEnd)] = ["e", -1, False] 
        elif (num_blocks == 2): # Is splice junction
            isoform_parts_list[transcript_id][(chrom, chromStart, chromEnd)] = ["sj", -1, False]
        else:
            print >> sys.stderr, "ERROR: number of blocks > 2.", line
            sys.exit(1)
            
        if (isoform_parts_list[transcript_id].has_key("strand")):
            assert (isoform_parts_list[transcript_id]["strand"] == strand)
        else:
            isoform_parts_list[transcript_id]["strand"] = strand

    ip.close()
    
    return isoform_parts_list


def annotateAndOrderTargetIsoformParts(all_exonic_parts_tbl, isoform_parts_list, target_isoforms, longform_otp_annotation=False):
    """Coordinates of parts are in BED format, which is 0-based, start
    position is inclusive and end position is exclusive (ie 1 position
    passed end of elemement."""
    # Add all of the exonic parts
    rePartsLine = re.compile("^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)")
    ip = open(all_exonic_parts_tbl, 'r')
    for line in ip:
        mo = rePartsLine.match(line)
        partChrom, partChromStart, partChromEnd, transcript_id, is_unique_region = mo.groups()
        if (transcript_id in target_isoforms):
            is_unique_region = int(is_unique_region) == 1
            partChromStart = int(partChromStart)
            partChromEnd = int(partChromEnd)
            part_tup = (partChrom, partChromStart, partChromEnd)
            if (isoform_parts_list[transcript_id].has_key(part_tup)):
                isoform_parts_list[transcript_id][part_tup][0] = "ue" if (is_unique_region) else "se"
            else:
                isoform_parts_list[transcript_id][part_tup] = ["ue", -1, False] if (is_unique_region) else ["se", -1, False]
    ip.close()

    # Order the parts. Only order "ue" or "se" exonic parts (in addition to "sj"), as "e" parts are superflous to them.
    for transcript_id in isoform_parts_list.keys(): 
        ordered_transcript_parts = []
        part_def_and_data = []
        for tup in isoform_parts_list[transcript_id].items():
            if (tup[0] != "strand"):
                part_def_and_data.append(tup)
        part_def_and_data.sort(partsSort)
        counter = 0
        for def_data in part_def_and_data:
            if (def_data[1][0] != "e"): 
                def_data[1][1] = counter
                counter += 1
                if (longform_otp_annotation):
                    ordered_transcript_parts.append( (def_data[0][0], def_data[0][1], def_data[0][2], def_data[1][0]) ) # chrom, start, end, part_type (e.g. "ue", "se", "sj")
                else:
                    ordered_transcript_parts.append( (def_data[0][1], def_data[0][2], def_data[1][0] == "ue" or def_data[1][0] == "se") ) # start, end, is_exon

        isoform_parts_list[transcript_id]["otp"] = ordered_transcript_parts


def annotateAndOrderIsoformParts(all_exonic_parts_tbl, isoform_parts_list, longform_otp_annotation=False):
    """Coordinates of parts are in BED format, which is 0-based, start
    position is inclusive and end position is exclusive (ie 1 position
    passed end of elemement."""
    # Add all of the exonic parts
    rePartsLine = re.compile("^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)")
    ip = open(all_exonic_parts_tbl, 'r')
    for line in ip:
        mo = rePartsLine.match(line)
        partChrom, partChromStart, partChromEnd, transcript_id, is_unique_region = mo.groups()
        is_unique_region = int(is_unique_region) == 1
        partChromStart = int(partChromStart)
        partChromEnd = int(partChromEnd)
        part_tup = (partChrom, partChromStart, partChromEnd)
        if (isoform_parts_list[transcript_id].has_key(part_tup)):
            isoform_parts_list[transcript_id][part_tup][0] = "ue" if (is_unique_region) else "se"
        else:
            isoform_parts_list[transcript_id][part_tup] = ["ue", -1, False] if (is_unique_region) else ["se", -1, False]
    ip.close()

    # Order the parts. Only order "ue" or "se" exonic parts (in addition to "sj"), as "e" parts are superflous to them.
    for transcript_id in isoform_parts_list.keys(): 
        ordered_transcript_parts = []
        part_def_and_data = []
        for tup in isoform_parts_list[transcript_id].items(): # tup is ((partChrom, partChromStart, partChromEnd), ["ue"|"se"|"sj", -1, False])
            if (tup[0] != "strand"):
                part_def_and_data.append(tup)
        part_def_and_data.sort(partsSort)
        counter = 0
        for def_data in part_def_and_data:
            if (def_data[1][0] != "e"): 
                def_data[1][1] = counter
                counter += 1
                if (longform_otp_annotation):
                    ordered_transcript_parts.append( (def_data[0][0], def_data[0][1], def_data[0][2], def_data[1][0]) ) # chrom, start, end, part_type (e.g. "ue", "se", "sj")
                else:
                    ordered_transcript_parts.append( (def_data[0][1], def_data[0][2], def_data[1][0] == "ue" or def_data[1][0] == "se") ) # start, end, is_exon

        isoform_parts_list[transcript_id]["otp"] = ordered_transcript_parts


# Creates a 0-based mRNA index to 1-based genomic index tuple list
def correspondmRNACoordsToGenomicCoords(mRNA_sequence, otp, strand):
    mRNA_indices = range(len(mRNA_sequence))
    genomic_indices = []
    for tup in otp:
        if (tup[3] != "sj"):
            genomic_indices.extend( range(tup[1]+1,tup[2]+1) )
     assert (len(mRNA_indices) == len(genomic_indices))
    if (strand == "-"):
        genomic_indices.reverse()

    return zip(mRNA_indices, genomic_indices)


def annotateSpliceJunctionUniqueness(all_splice_junction_parts_tbl, isoform_parts_list):
    """All input is assumed to be for unique splice junctions."""
    ip = open(all_splice_junction_parts_tbl, 'r')
    for line in ip:
        elems = line.strip().split("\t")
        if (int(elems[-1]) == 1):
            isoform_parts_list[elems[3]][(elems[0], int(elems[1]), int(elems[2]))][0] = "usj" 

    ip.close()
