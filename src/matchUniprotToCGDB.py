#!/usr/bin/env python

from operator import itemgetter
from itertools import combinations
from collections import defaultdict
import bz2
import copy
import networkx as nx
import gzip
import os
import pandas as pd
import re
import sys

import pdb
    
def linkUniprotToCGDB(uniprot_to_refseq, uniprot_to_gencode, CGDB_absorption):
    reAbsorption1 = re.compile("^(\S+) <- (\S+)$")
    reAbsorption2 = re.compile("^(CG_\S+) \((\S+)\) <- (\S+)$")

    transcript_to_uniprot = defaultdict(set)
    with open(uniprot_to_refseq, 'r') as ip:
        header = ip.readline()
        for line in ip:
            uniprot_id, refseq_transcript_id = line.strip().split("\t")
            refseq_transcript_id = refseq_transcript_id.split('.')[0]
            transcript_to_uniprot[refseq_transcript_id].add(uniprot_id)
            
    with open(uniprot_to_gencode, 'r') as ip:
        header = ip.readline()
        for line in ip:
            uniprot_id, gencode_transcript_id = line.strip().split("\t")
            gencode_transcript_id = gencode_transcript_id.split('.')[0]
            transcript_to_uniprot[gencode_transcript_id].add(uniprot_id)

    # Example lines:
    #   CG_102132.1.1.1 <- klorfeyby.aAug10
    #   None (ENST00000544132.2) <- uc058yim
    #   CG_99547.1.6.14 (XM_017020835) <- HIT000220173
    uniprot_to_CGDB_tuples = []
    ip = open(CGDB_absorption,'r')
    for line in ip:
        mo1 = reAbsorption1.match(line)
        mo2 = reAbsorption2.match(line)

        if (mo1 != None):
            absorber_isoform, absorbee_isoform = mo1.groups()
            assert (absorber_isoform != "None" and absorber_isoform.startswith("CG_"))
            absorbee_isoform = absorbee_isoform.rsplit('.',1)[0]
            if (absorbee_isoform in transcript_to_uniprot):
                uniprot_to_CGDB_tuples.append( (transcript_to_uniprot[absorbee_isoform], absorber_isoform) )
        elif (mo2 != None):
            CG_isoform, absorber_isoform, absorbee_isoform = mo2.groups()
            assert (CG_isoform.startswith("CG_"))

            absorber_isoform = absorber_isoform.rsplit('.',1)[0]
            if (absorber_isoform in transcript_to_uniprot):
                uniprot_to_CGDB_tuples.append( (transcript_to_uniprot[absorber_isoform], CG_isoform) )

            absorbee_isoform = absorbee_isoform.rsplit('.',1)[0]
            if (absorbee_isoform in transcript_to_uniprot):
                uniprot_to_CGDB_tuples.append( (transcript_to_uniprot[absorbee_isoform], CG_isoform) )
        elif (not line.startswith("None")):
            print >> sys.stderr, "ERROR: did not parse line %s" % line
            sys.exit(1)
            
    ip.close()

    uniprot_to_CGDB_link = {}

    uniprot_w_transcript_link = set.union(*transcript_to_uniprot.values())

    # NOTE: UniProt transcripts that don't have an associated CGDB isoform reside on non-canonical chromosomes/scaffolds that are not included in CGDB.
    uniprot_linked_to_CGDB_transcript = set.union(*map(itemgetter(0), uniprot_to_CGDB_tuples))
    for uniprot_id in uniprot_w_transcript_link - uniprot_linked_to_CGDB_transcript:
        uniprot_to_CGDB_link[uniprot_id] = None

    for uniprot_ids, CGDB_ID in uniprot_to_CGDB_tuples:
        for uniprot_id in uniprot_ids:
            if (not uniprot_to_CGDB_link.has_key(uniprot_id)):
                uniprot_to_CGDB_link[uniprot_id] = set()
            uniprot_to_CGDB_link[uniprot_id].add(CGDB_ID)

    return (uniprot_to_CGDB_link, uniprot_w_transcript_link)


def mapUniprotToCGDB(uniprot_w_transcript_link, annotated_PM_proteins, blast_directory):
    uniprot_to_CGDB_map = {}

    blast_file_columns = ["query", "subject", "perc_ident", "query_len", "subject_len", "align_len", "num_mismatches", "num_gaps", \
                          "query_start", "query_end", "subject_start", "subject_end", "E_value", "bit_score"]

    ip = gzip.open(annotated_PM_proteins, 'rb')
    header = ip.readline()

    # Collect the isoforms that match each UniProt
    for line in ip:
        fields = line.strip().split("\t")
        uniprot_id = fields[0]
        if (uniprot_id not in uniprot_w_transcript_link):
            print >> sys.stderr, "\rINFO: processing Blast results for %s" % uniprot_id
            blast_file = "%s/%s.tblastn.bz2" % (blast_directory, uniprot_id)
            assert(os.path.exists(blast_file))
            
            assert (not uniprot_to_CGDB_map.has_key(uniprot_id))
            uniprot_to_CGDB_map[uniprot_id] = set()

            blast_df = pd.read_csv(blast_file, sep="\t", header=None, skiprows=5, names=blast_file_columns, compression="bz2", skipfooter=1)
            most_identities_match = 0
            for tup in blast_df.itertuples(index=False):
                perc_ident = float(tup[2])
                query_len = float(tup[3])
                align_len = float(tup[5])
                num_identities = perc_ident * align_len / 100.0
                if (num_identities >= most_identities_match or (perc_ident >= 98.0 and align_len/query_len >= 0.98)):
                    most_identities_match = max(num_identities,most_identities_match)
                    uniprot_to_CGDB_map[uniprot_id].add(tup[1])

            print >> sys.stderr, "INFO: added %d for %s" % (len(uniprot_to_CGDB_map[uniprot_id]), uniprot_id)

    ip.close()

    return uniprot_to_CGDB_map


def annotateCGDB(uniprot_to_CGDB_link, uniprot_to_CGDB_map, annotated_PM_proteins):
    CGDB_from_uniprot = []

    assert(len(set(uniprot_to_CGDB_link.keys()) & set(uniprot_to_CGDB_map.keys())) == 0), "UniProt proteins used for linking and mapping to CGDB isoforms should be disjoint sets"
    uniprot_to_CGDB = dict(uniprot_to_CGDB_link.items() +  uniprot_to_CGDB_map.items())

    ip = gzip.open(annotated_PM_proteins, 'rb')
    header = ip.readline()

    # Collect the isoforms that match each UniProt
    for line in ip:
        fields = line.strip().split("\t")
        uniprot_id = fields[0]
        if (uniprot_to_CGDB[uniprot_id]==None): # See Note above
            continue
        has_TM_segments = fields[2] != "None"
        has_GPI = fields[3] == "True"
                
        CGDB_from_uniprot.append( (uniprot_to_CGDB[uniprot_id], uniprot_id, has_TM_segments, has_GPI) )
    ip.close()

    return CGDB_from_uniprot


########## Below unused ######
def mapIsoformGroupToCGDBIsoforms(isoform_group_definitions):
    IG_to_CGDB_isoforms = {}
    ip = bz2.BZ2File(isoform_group_definitions, 'rb')
    header = ip.readline()
    header_fields = header.strip().split('\t')
    assert(header_fields == ['GroupID', 'IsoformIDs', 'GeneSymbol(s)', 'pORFs', 'ProteinSequence']), "Unexpected columns in %s" % isoform_group_definitions

    for line in ip:
        fields = line.split('\t')
        IG_to_CGDB_isoforms[fields[0]] = set(fields[1].split('|'))
    ip.close()

    return IG_to_CGDB_isoforms


def matchUniprotProteinsToIsoformGroup(IG_to_CGDB_isoforms, annotated_PM_proteins, blast_directory):
    CGDB_from_uniprot = []

    blast_file_columns = ["query", "subject", "perc_ident", "query_len", "subject_len", "align_len", "num_mismatches", "num_gaps", \
                          "query_start", "query_end", "subject_start", "subject_end", "E_value", "bit_score"]

    ip = gzip.open(annotated_PM_proteins, 'rb')
    header = ip.readline()

    # Collect the isoforms that match each UniProt
    for line in ip:
        fields = line.strip().split("\t")
        uniprot_id = fields[0]
        has_TM_segments = fields[2] != "None"
        has_GPI = fields[3] == "True"
        print >> sys.stderr, "\rINFO: processing Blast results for %s" % uniprot_id
        blast_file = "%s/%s.tblastn.bz2" % (blast_directory, uniprot_id)
        assert(os.path.exists(blast_file))

        blast_df = pd.read_csv(blast_file, sep="\t", header=None, skiprows=5, names=blast_file_columns, compression="bz2", skipfooter=1)
        most_identities_match = 0
        CGDB_isoforms = set()
        for tup in blast_df.itertuples(index=False):
            num_identities = float(tup[2]) * float(tup[5])
            if (num_identities >= most_identities_match or tup[3]==tup[5]):
                most_identities_match = max(num_identities,most_identities_match)
                CGDB_isoforms.update( IG_to_CGDB_isoforms[tup[1]] )
                
        CGDB_from_uniprot.append( (CGDB_isoforms, uniprot_id, has_TM_segments, has_GPI) )

    ip.close()
    print >> sys.stderr, "INFO: done\n"

    return CGDB_from_uniprot


def expandIsoformMatches(CGDB_from_uniprot, CGDB_exons):
    '''Expand the group of isoforms associated with each uniprot to those
    isoforms that share an exon boundary with, or are from the same
    locus as, the core set of isoforms already identified.
    '''
    expanded_CGDB_from_uniprot = []
    locus_isoforms = defaultdict(set)

    core_isoforms = set.union(*map(itemgetter(0), CGDB_from_uniprot))
    core_loci = set(map(lambda x:x.split('.')[0], core_isoforms))

    G = nx.Graph()
    G.add_nodes_from(core_isoforms)
    for isoforms_set in map(itemgetter(0), CGDB_from_uniprot):
        for isoform1, isoform2 in combinations(isoforms_set,2):
            G.add_edge(isoform1, isoform2)

    # Record the exons boundaries for each core isoform
    core_isoforms_w_exon_boundary = defaultdict(set)
    ip_exons = open(CGDB_exons, 'r')
    for line in ip_exons:
        chrom, start, stop, isoform_ID, score, strand = line.strip().split("\t")
        locus = isoform_ID.split('.',1)[0]
        locus_isoforms[locus].add(isoform_ID)
        if (isoform_ID in core_isoforms):
            core_isoforms_w_exon_boundary[(chrom, start, strand)].add(isoform_ID)
            core_isoforms_w_exon_boundary[(chrom, stop, strand)].add(isoform_ID)

    # Connect new isoforms that share an exon
    ip_exons.seek(0)
    

    # For each core isoform, link it to every isoform with which it shares an exon boundary
    ## and to the rest of the isoforms from those isoforms' loci
    for line in ip_exons:
        chrom, start, stop, isoform_ID, score, strand = line.strip().split("\t")
        exon_boundary_descr1 = (chrom, start, strand)
        exon_boundary_descr2 = (chrom, stop, strand)
        if (core_isoforms_w_exon_boundary.has_key(exon_boundary_descr1) and isoform_ID not in core_isoforms_w_exon_boundary[exon_boundary_descr1]):
            #noncore_isoform_locus = isoform_ID.split('.',1)[0]
            for core_isoform_ID in core_isoforms_w_exon_boundary[exon_boundary_descr1]:
                #for locus_isoform in locus_isoforms[noncore_isoform_locus]:
                #    G.add_edge(locus_isoform, core_isoform_ID)
                G.add_edge(isoform_ID, core_isoform_ID)
                
        if (core_isoforms_w_exon_boundary.has_key(exon_boundary_descr2) and isoform_ID not in core_isoforms_w_exon_boundary[exon_boundary_descr2]):
            #noncore_isoform_locus = isoform_ID.split('.',1)[0]
            for core_isoform_ID in core_isoforms_w_exon_boundary[exon_boundary_descr2]:
                #for locus_isoform in locus_isoforms[noncore_isoform_locus]:
                #    G.add_edge(locus_isoform, core_isoform_ID)
                G.add_edge(isoform_ID, core_isoform_ID)

    ip_exons.close()


    for core_CGDB_isoforms, uniprot_id, has_TM_segments, has_GPI in CGDB_from_uniprot:
        expanded_set = copy.deepcopy(core_CGDB_isoforms)
        for core_isoform in core_CGDB_isoforms:
            expanded_set.update(G.neighbors(core_isoform))
        expanded_CGDB_from_uniprot.append( (expanded_set, uniprot_id, has_TM_segments, has_GPI) )

    return (expanded_CGDB_from_uniprot, core_loci)


def writeOutput(expanded_CGDB_from_uniprot, core_loci, output_tsv):
    print >> sys.stderr, "INFO: writing output to %s" % output_tsv
    op = gzip.open(output_tsv, 'wb')
    op.write("CGDB_ISOFORM\tTM_PROTEIN_ID\tGPI_PROTEIN_ID\tUNIPROT_ASSOC\n")
    for CGDB_isoforms, uniprot_id, has_TM_segments, has_GPI in expanded_CGDB_from_uniprot:
        tm_val = uniprot_id if (has_TM_segments) else '-'
        gpi_val = uniprot_id if (has_GPI) else '-'
        for CGDB_isoform in CGDB_isoforms:
            isoform_locus = CGDB_isoform.split('.',1)[0]
            if (isoform_locus in core_loci):
                op.write("%s\t%s\t%s\tdirect\n" % (CGDB_isoform, tm_val, gpi_val))
            else:
                op.write("%s\t%s\t%s\tindirect\n" % (CGDB_isoform, tm_val, gpi_val))
    op.close()


if (__name__ == "__main__"):
    uniprot_to_refseq, uniprot_to_gencode, CGDB_absorption, isoform_group_definitions, annotated_PM_proteins, CGDB_exons, blast_directory, output_tsv = sys.argv[1:]

    print >> sys.stderr, "INFO: Linking...",
    uniprot_to_CGDB_link, uniprot_w_transcript_link = linkUniprotToCGDB(uniprot_to_refseq, uniprot_to_gencode, CGDB_absorption)
    print >> sys.stderr, "done"

    print >> sys.stderr, "INFO: Mapping...",
    uniprot_to_CGDB_map = mapUniprotToCGDB(uniprot_w_transcript_link, annotated_PM_proteins, blast_directory)
    print >> sys.stderr, "done"

    #print >> sys.stderr, "INFO: Mapping isoform group to isoforms...",
    #IG_to_CGDB_isoforms = mapIsoformGroupToCGDBIsoforms(isoform_group_definitions)
    #print >> sys.stderr, "done"

    CGDB_from_uniprot = annotateCGDB(uniprot_to_CGDB_link, uniprot_to_CGDB_map, annotated_PM_proteins)
    # Associate each UniProt to it's source isoform group
    #CGDB_from_uniprot = matchUniprotProteinsToIsoformGroup(IG_to_CGDB_isoforms, annotated_PM_proteins, blast_directory)

    print >> sys.stderr, "INFO: Expanding UniProt matches...",
    expanded_CGDB_from_uniprot, core_loci = expandIsoformMatches(CGDB_from_uniprot, CGDB_exons)
    print >> sys.stderr, "done"

    writeOutput(expanded_CGDB_from_uniprot, core_loci, output_tsv)

    sys.exit(0)
