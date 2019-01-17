#!/usr/bin/env python

from collections import defaultdict
import cPickle
import itertools
import gzip
import re
from operator import itemgetter
import sys

import pdb

class HMPAS:
    def __init__(self, hmpas_ID, membrane_interaction_type, subcellular_location="", protein_info=[]):
        self.hmpas_ID = hmpas_ID
        self.membrane_interaction_type = membrane_interaction_type
        self.subcellular_locations = []
        if (subcellular_location != ""):
            self.subcellular_locations.append(subcellular_location)
        self.sequence = None
        self.sequence_incomplete = None
        self.membrane_topology = {"HMMTOP":[], "PHOBIUS":[], "S-TMHMM":[], "SCAMPI":[], "TMHMM":[], "PDBTM":[], "UniProt Sequence annotation":[], "Fallback":[]}
        self.gpi = {"FragAnchor":None}
        
        # Names
        self.gene_name = "None"
        self.gene_symbol = "None"
        for name_type, name in filter(lambda x:len(x)==2, protein_info):
            if (name_type == "Protein Name"):
                self.gene_name = name
            elif (name_type == "Gene Name"):
                self.gene_symbol = name

                
    def evalCompleteness(self, target_type):
        evaluation = ""
        if (self.gene_name == "None" and self.gene_symbol == "None"):
            evaluation += "No associated gene name/symbol."
        if (self.sequence == None):
            evaluation += " No amino acid sequence."

        if (target_type == "Integral membrane protein"):
            if (len(self.getTransmembraneSegments("all")) == 0):
                if ("No amino acid sequence" not in evaluation):
                    evaluation += " No membrane topology specified. Setting Fallback to entire sequence."
                    #self.membrane_topology["Fallback"].append( ("transmembrane", 1, len(self.sequence)) )
        elif (target_type == "Lipid-anchored protein"):
            if (all(map(lambda x:x==None, self.gpi.values()))):
                evaluation += " No GPI sites specified."
        else:
            print >> sys.stderr, "ERROR: unrecognized type. Exiting."
            sys.exit(1)
            
        return evaluation

    
    def incompleteCDS(self):
        return self.sequence_incomplete


    def addAminoAcidSequence(self, sequence):
        self.sequence = sequence.strip()
        self.sequence_incomplete = self.sequence[0]=='X' or self.sequence[-1]=='X'


    def addSubcellularLocation(self, subcellular_location):
        self.subcellular_locations.append(subcellular_location)        


    def getSubcellularLocation(self):
        if (len(self.subcellular_locations)==0):
            self.subcellular_locations = [('',"Unknown")]
        return self.subcellular_locations
    

    def hasPlasmaMembraneLocalization(self):
        return any(map(lambda x: x[0].startswith("SL:1"), self.subcellular_locations))

        
    def addMembraneTopologyInfo(self, prediction_method, topology_segment):
        """topology_segment is a tuple (["inside"|"outside"|"transmembrane"], start, end)"""
        self.membrane_topology[prediction_method].append( topology_segment )


    def addGPIInfo(self, prediction_method, cleavage_position):
        self.gpi[prediction_method] = cleavage_position


    def normalizeTopologyData(self):
        for method, topology_segments in self.membrane_topology.items():
            self.membrane_topology[method] = sorted(topology_segments, key=itemgetter(1))


    def getGeneSymbolAndName(self):
        return (self.gene_symbol, self.gene_name)


    def getAminoAcidSequence(self):
        return self.sequence
    
    
    def getTransmembraneSegments(self, method):
        """Returns just the tuples (start, end) for the transmembrane segments"""
        segments = set()

        if (method == "all"):
            for topology_segments in self.membrane_topology.values():
                segments.update( map(lambda y:y[1:], filter(lambda x:x[0]=="transmembrane", topology_segments)) )
        else:
            topology_segments = self.membrane_topology[method]
            segments.update( map(lambda y:y[1:], filter(lambda x:x[0]=="transmembrane", topology_segments)) )

        return segments


    def getGPICleavageSites(self, method):
        """Returns a list of the cleavage site basepair positions."""
        sites = []
        if (method == "all"):
            sites.extend( filter(lambda x:x!=None, self.gpi.values()) )
        else:
            if (self.gpi[method] != None):
                sites.append( self.gpi[method] )
        return sites


def readProteinInfo(HMPAS_ProteinInfo):
    protein_info = defaultdict(list)
    ip = open(HMPAS_ProteinInfo, 'r')
    for line in ip:
        elems = line.strip().split("\t")
        protein_info[elems[0]].append( tuple(elems[1:]) )
    ip.close()
    return protein_info


def getAllHMPASofType(HMPAS_ClassMember_file, target_type, HMPAS_ProteinInfo_file=None, HMPAS_Sequence_file=None, HMPAS_SequenceFeature_file=None):
    hmpas_of_type = {}

    protein_info = readProteinInfo(HMPAS_ProteinInfo_file)
        
    acceptable_evidence = set(["Experiment Evidence", "Literature Evidence", "Manual Curation Evidence", "Ortholog based Prediction"])
    unacceptable_evidence = set(["Computational Prediction Evidence", "IMP Classifier", "LAP Classifier", "Non-traceable Evidence", "PMP Classifier"])
    
    # Collect the HMPAS entries of the target type
    ip = open(HMPAS_ClassMember_file, 'r')
    for line in ip:
        elems = line.strip().split("\t")
        if (elems[-1] in acceptable_evidence):
            if (elems[1] == "Interaction type with membrane" and elems[3]==target_type):
                hmpas_of_type[elems[0]] = HMPAS(elems[0], elems[3], "", protein_info[elems[0]])
            elif (elems[1] == "Subcellualr Localization" and hmpas_of_type.has_key(elems[0])):
                hmpas_of_type[elems[0]].addSubcellularLocation( (elems[2],elems[3]) )
        else:
            try:
                assert (elems[-1] in unacceptable_evidence)
            except AssertionError:
                print elems
                pdb.set_trace()
    ip.close()

    # Add the amino acid sequences
    ip = open(HMPAS_Sequence_file, 'r')
    for line in ip:
        hmpas_ID, seqlen, aa_seq = line.split("\t")
        if (hmpas_of_type.has_key(hmpas_ID)):
            hmpas_of_type[hmpas_ID].addAminoAcidSequence(aa_seq)
    ip.close()
    
    # Add sequence feature information
    ip = open(HMPAS_SequenceFeature_file, 'r')
    if (target_type == "Integral membrane protein"):
        for line in ip:
            elems = line.strip().split("\t")
            if (elems[1] == "Membrane Topology" and hmpas_of_type.has_key(elems[0])):
                hmpas_of_type[elems[0]].addMembraneTopologyInfo(elems[2], (elems[3],int(elems[5]),int(elems[6])))

    elif (target_type == "Lipid-anchored protein"):
        for line in ip:
            elems = line.strip().split("\t")
            if (elems[1] == "Lipid-Anchor" and elems[3] == "GPI-anchor" and hmpas_of_type.has_key(elems[0])):
                hmpas_of_type[elems[0]].addGPIInfo(elems[2], int(elems[-1]))
    ip.close()

    # Validate completeness
    to_remove = set()
    for hmpas_ID, hmpas_object in hmpas_of_type.items():
        evaluation = hmpas_object.evalCompleteness(target_type)
        if (evaluation != ""):
            #print >> sys.stderr, "Removing", hmpas_ID, hmpas_object.gene_name, hmpas_object.gene_symbol
            #print >> sys.stderr, "WARNING: Validation issues for %s: %s" % (hmpas_ID, evaluation)
            to_remove.add(hmpas_ID)
            #if ("Fallback" not in evaluation):
            #    pdb.set_trace()
        
    for hmpas_ID in to_remove:
        del hmpas_of_type[hmpas_ID]

    return hmpas_of_type
    

def identifyTMpORFs(match_data, hmpas_TM, min_TM_match_len=15):
    pORFs_aligned_to_TM_segments = defaultdict(set)

    for hmpas_ID, hmpas_seq_info in hmpas_TM.items():
        pORF_matches = match_data[hmpas_ID]
        TM_segments = hmpas_seq_info.getTransmembraneSegments("all")
        hmpas_seq = hmpas_seq_info.getAminoAcidSequence()
        gene_symbol, gene_name = hmpas_seq_info.getGeneSymbolAndName()
        incomplete_CDS = hmpas_seq_info.incompleteCDS()

        TM_positions = []
        for TM_start, TM_stop in TM_segments:
            TM_positions.append( set(range(TM_start,TM_stop+1)) )

        # If the pORF sufficiently aligned to a known HMPAS protein TM segment, record
        for HMPAS_align_start, HMPAS_align_stop, pORF_ID in pORF_matches:
            match_positions = set(range(HMPAS_align_start, HMPAS_align_stop+1))
            if (any(map(lambda x: len(match_positions & x) >= min_TM_match_len, TM_positions))):
                pORFs_aligned_to_TM_segments[pORF_ID].add(hmpas_ID)

    return pORFs_aligned_to_TM_segments


def getAllHPA(HPA_PM_proteins):
    all_hpa = defaultdict(set)
    ip = gzip.open(HPA_PM_proteins, 'rb')
    header = ip.readline().strip().split("\t")
    assert (header[2] == "UniProt/SwissProt Accession" and header[3] == "UniProt/TrEMBL Accession" and
            header[6] == "Main location" and header[7] == "Other location" and header[9] == "Reliability"), "Incorrectly formatted file -> %s" % HPA_PM_proteins
        
    for line in ip:
        fields = line.strip().split("\t")
        #try:
        assert ("Plasma membrane" in fields[6] or "Plasma membrane" in fields[7])
        assert (fields[9] == "Supportive")
        #except AssertionError:
        #    pdb.set_trace()

        ensembl_gene = fields[0]
        swiss_acc = fields[2].strip()
        if (swiss_acc != ''):
            all_hpa[ensembl_gene].add(swiss_acc)
        else:
            trembl_acc = fields[3].strip()
            all_hpa[ensembl_gene].add(trembl_acc)
    ip.close()
    return all_hpa


def getAllFromHMPAS(HMPAS_ClassMember_file, HMPAS_ProteinInfo_file, HMPAS_Sequence_file, HMPAS_SequenceFeature_file, has_GPI, TM_segments, protein_seqs):
    # Collect annotation information for HMPAS TM proteins and GPI proteins
    hmpas_TM = getAllHMPASofType(HMPAS_ClassMember_file, "Integral membrane protein", HMPAS_ProteinInfo_file, HMPAS_Sequence_file, HMPAS_SequenceFeature_file)
    hmpas_gpi = getAllHMPASofType(HMPAS_ClassMember_file, "Lipid-anchored protein", HMPAS_ProteinInfo_file, HMPAS_Sequence_file, HMPAS_SequenceFeature_file)

    for uniprot_ID in set(hmpas_TM.keys()) | set(hmpas_gpi.keys()):
        if (hmpas_TM.has_key(uniprot_ID)):
            hmpas_object = hmpas_TM[uniprot_ID]
            if (hmpas_object.hasPlasmaMembraneLocalization()):
                TM_ranges = hmpas_object.getTransmembraneSegments("UniProt Sequence annotation") # "all"
                if (len(TM_ranges) == 0):
                    print >> sys.stderr, "WARNING: HMPAS integral membrane UniProt protein %s has no UniProt transmembrane segments. Skipping." % uniprot_ID
                else:
                    TM_segments[uniprot_ID].update( TM_ranges )
                    protein_seqs[uniprot_ID] = hmpas_object.getAminoAcidSequence()
        if (hmpas_gpi.has_key(uniprot_ID)):
            has_GPI[uniprot_ID] = True
            if (not protein_seqs.has_key(uniprot_ID)):
                hmpas_object = hmpas_gpi[uniprot_ID]
                protein_seqs[uniprot_ID] = hmpas_object.getAminoAcidSequence()

    print >> sys.stderr, "INFO: added %d Plasma Membrane from HMPAS" % len(protein_seqs.keys())


def getAllFromUniProt(Uniprot_PM_transmembrane_proteins, Uniprot_PM_GPI_proteins, has_GPI, TM_segments, protein_seqs):
    re_TM_segment = re.compile(r"\s+(\d+)\s+(\d+)\s+([A-Za-z ]+)")

    # Get all Uniprot that are annotated as "cell membrane" and have a transmembrane segment
    for uniprot_file in [Uniprot_PM_transmembrane_proteins, Uniprot_PM_GPI_proteins]:
        num_new_added = 0
        ip = gzip.open(uniprot_file, 'rb')
        header = ip.readline().strip().split("\t")
        assert (header == ['Entry', 'Entry name', 'Protein names', 'Gene names', 'Transmembrane', 'Intramembrane', 'Topological domain', 'Subcellular location [CC]', 'Sequence']), \
            "Uniprot file doesn't have expected columns/column order"

        for line in ip:
            fields = line.strip().split("\t")
            if (len(fields) != 9):
                if (fields[2].startswith("Deleted")):
                    continue
                else:
                    pdb.set_trace()
            uniprot_acc = fields[0].strip()
            TM_field = fields[4].strip()
            CC_field = fields[7].strip()
            aa_seq = fields[8].strip()
            if (TM_field != ''):
                for mo in map(lambda x: re_TM_segment.search(x), TM_field.split('TRANSMEM')[1:]):
                    start, stop, tm_descriptor = mo.groups()
                    assert (tm_descriptor in ["Helical", "Discontinuously helical", "Beta stranded"]), "Unrecognized transmembrane segment descriptor"
                    TM_segments[uniprot_acc].add( (int(start), int(stop)) )
            if (protein_seqs.has_key(uniprot_acc)):
                if (protein_seqs[uniprot_acc] != aa_seq):
                    protein_seqs[uniprot_acc] = aa_seq # Use the presumably more up-to-date data directly from Uniprot over that from the dated HMPAS database
            else:
                if (uniprot_file == HPA_PM_proteins):
                    print >> sys.stderr, uniprot_acc
                protein_seqs[uniprot_acc] = aa_seq
                num_new_added += 1

            if ("GPI-anchor" in CC_field):
                has_GPI[uniprot_acc] = True
                if (protein_seqs.has_key(uniprot_acc)):
                    assert (protein_seqs[uniprot_acc] == aa_seq)
                else:
                    protein_seqs[uniprot_acc] = aa_seq
                    num_new_added += 1

                    
        ip.close()
        print >> sys.stderr, "INFO: added %d from %s" % (num_new_added, uniprot_file)


def getAllFromHPA(HPA_PM_proteins, has_GPI, TM_segments, protein_seqs):
    re_TM_segment = re.compile(r"\s+(\d+)\s+(\d+)\s+([A-Za-z ]+)")

    # Get all Uniprot that are annotated as "cell membrane" and have a transmembrane segment
    for uniprot_file in [Uniprot_PM_transmembrane_proteins, Uniprot_PM_GPI_proteins, HPA_PM_proteins]:
        num_new_added = 0
        ip = gzip.open(uniprot_file, 'rb')
        header = ip.readline().strip().split("\t")
        assert (header == ['Entry', 'Entry name', 'Protein names', 'Gene names', 'Transmembrane', 'Intramembrane', 'Topological domain', 'Subcellular location [CC]', 'Sequence']), \
            "Uniprot file doesn't have expected columns/column order"

        for line in ip:
            fields = line.strip().split("\t")
            if (len(fields) != 9):
                if (fields[2].startswith("Deleted")):
                    continue
                else:
                    pdb.set_trace()
            uniprot_acc = fields[0].strip()
            TM_field = fields[4].strip()
            CC_field = fields[7].strip()
            aa_seq = fields[8].strip()
            if (TM_field != ''):
                for mo in map(lambda x: re_TM_segment.search(x), TM_field.split('TRANSMEM')[1:]):
                    start, stop, tm_descriptor = mo.groups()
                    assert (tm_descriptor in ["Helical", "Discontinuously helical", "Beta stranded"]), "Unrecognized transmembrane segment descriptor"
                    TM_segments[uniprot_acc].add( (int(start), int(stop)) )
                if (protein_seqs.has_key(uniprot_acc)):
                    if (protein_seqs[uniprot_acc] != aa_seq):
                        protein_seqs[uniprot_acc] = aa_seq # Use the presumably more up-to-date data directly from Uniprot over that from the dated HMPAS database
                else:
                    print >> sys.stderr, uniprot_acc
                    protein_seqs[uniprot_acc] = aa_seq
                    num_new_added += 1

            if ("GPI-anchor" in CC_field):
                has_GPI[uniprot_acc] = True
                if (protein_seqs.has_key(uniprot_acc)):
                    assert (protein_seqs[uniprot_acc] == aa_seq)
                else:
                    protein_seqs[uniprot_acc] = aa_seq
                    num_new_added += 1

                    
        ip.close()
        print >> sys.stderr, "INFO: added %d from %s" % (num_new_added, uniprot_file)



def getAllUniProt(Uniprot_PM_proteins, Uniprot_all_proteins, all_hpa, uniprot_from_hmpas):
    re_TM_segment = re.compile(r"\s+(\d+)\s+(\d+)\s+([A-Za-z ]+)")

    new_uniprot = {}

    # Get all Uniprot that are annotated as "cell membrane" and have a transmembrane segment
    num_new_uniprot_cell_membrane = 0
    ip = gzip.open(Uniprot_PM_proteins, 'rb')
    header = ip.readline().strip().split("\t")
    assert (header[0] == "Entry" and header[4] == "Transmembrane" and header[5] == "Topological domain" and header[-1] == "Sequence")
        
    for line in ip:
        fields = line.strip().split("\t")
        uniprot_acc = fields[0].strip()
        TM_field = fields[4].strip()
        if (TM_field != '' and uniprot_acc not in uniprot_from_hmpas):
            TM_segments = set()
            for mo in map(lambda x: re_TM_segment.search(x), TM_field.split('TRANSMEM')[1:]):
                start, stop, tm_descriptor = mo.groups()
                assert (tm_descriptor in ["Helical", "Discontinuously helical", "Beta stranded"]), "Unrecognized transmembrane segment descriptor"
                TM_segments.add( (int(start), int(stop)) )
            assert (len(TM_segments) > 0)
            new_uniprot[fields[0]] = (fields[-1], TM_segments)
            num_new_uniprot_cell_membrane += 1
    ip.close()

    print >> sys.stderr, "INFO: added %d annotated 'cell membrane' in UniProt that were not in HMPAS" % num_new_uniprot_cell_membrane

    # Get all Uniprot that
    #    1) are not annotated as "cell membrane" but are plasma membrane proteins (because they have an "extracellular" topological
    #       component) that have a transmembrane segment, or
    #    2) HPA says are plasma membrane associated, are not annotated in Uniprot as cell membrane, but that have transmembrane segments.
    all_HPA_PM = set(itertools.chain(*all_hpa.values()))

    ip = gzip.open(Uniprot_all_proteins, 'rb')
    header = ip.readline().strip().split("\t")
    assert (header[0] == "Entry" and header[4] == "Transmembrane" and header[-1] == "Sequence")
        
    num_new_by_extracellular_annotation = 0
    num_new_by_HPA = 0
    num_new_by_both_extracellular_and_HPA = 0
    for line in ip:
        fields = line.strip().split("\t")
        uniprot_acc = fields[0].strip()
        if (not new_uniprot.has_key(uniprot_acc) and not uniprot_acc in uniprot_from_hmpas):
            TM_field = fields[4].strip()
            TOPO_field = fields[5].strip()

            has_TM_and_extracellular = TM_field != '' and "xtracellular" in TOPO_field
            is_PM_by_HPA_and_has_TM = fields[0] in all_HPA_PM and TM_field != ''

            if (has_TM_and_extracellular or is_PM_by_HPA_and_has_TM):
                TM_segments = []
                for mo in map(lambda x: re_TM_segment.search(x), TM_field.split('TRANSMEM')[1:]):
                    start, stop, tm_descriptor = mo.groups()
                    assert (tm_descriptor in ["Helical", "Discontinuously helical", "Beta stranded"]), "Unrecognized transmembrane segment descriptor"
                    TM_segments.append( (int(start), int(stop)) )
                assert (len(TM_segments) > 0)
                new_uniprot[uniprot_acc] = (fields[-1], TM_segments)
                if (has_TM_and_extracellular):
                    num_new_by_extracellular_annotation += 1
                    if (is_PM_by_HPA_and_has_TM):
                        num_new_by_both_extracellular_and_HPA += 1
                    #print "Found by annotation"
                else:
                    num_new_by_HPA +=1 
                    #print line
                    #print "Found by HPA"
                    #pdb.set_trace()
    ip.close()

    print >> sys.stderr, "INFO: added additional %d annotated 'extracellular' in Uniprot (%d also PM in HPA) and %d only annotated 'plasma membrane' in HPA" % \
      (num_new_by_extracellular_annotation, num_new_by_both_extracellular_and_HPA, num_new_by_HPA)

    return new_uniprot


def assessHPACoverageByUniprot(all_hpa, all_uniprot):
    '''How many genes from HPA have a uniprot entry?'''
    num_gene_w_uniprot, num_gene_wo_uniprot = 0,0

    for gene, uniprot_accessions in all_hpa.items():
        if (any(map(lambda x: all_uniprot.has_key(x), uniprot_accessions))):
            num_gene_w_uniprot += 1
        else:
            print >> sys.stderr, gene
            num_gene_wo_uniprot += 1

    print >> sys.stderr, "HPA w/ Uniprot = %d\tHPA w/o Uniprot = %d" % (num_gene_w_uniprot, num_gene_wo_uniprot)


def writeOutput(has_GPI, TM_segments, protein_seqs, output_tsv):
    #assert (set(protein_seqs.keys()) == set(has_GPI.keys() + TM_segments.keys()))

    op = gzip.open(output_tsv, 'wb')
    op.write("UNIPROT_ID\tAA_SEQ\tPM_TM_SEGMENTS\tHAS_GPI\n")

    for uniprot_ID, aa_seq in protein_seqs.items():
        if (TM_segments.has_key(uniprot_ID)):
            tm_segments = TM_segments[uniprot_ID]
            tm_segments = sorted(tm_segments)
            tm_segments_str = ",".join(map(lambda x: "%d-%d" % tuple(x), tm_segments))
        else:
            tm_segments_str = "None"

        gpi_str = "True" if has_GPI[uniprot_ID] else "False"

        op.write("%s\t%s\t%s\t%s\n" % (uniprot_ID, aa_seq, tm_segments_str, gpi_str))

    op.close()


if (__name__ == "__main__"):
    HPA_PM_proteins, Uniprot_PM_transmembrane_proteins, Uniprot_PM_GPI_proteins, \
        HMPAS_Sequence_file, HMPAS_ClassMember_file, HMPAS_SequenceFeature_file, HMPAS_ProteinInfo_file, output_tsv = sys.argv[1:]

    # Using HMPAS-, HPA-, and Uniprot-identified plasma membrane-localized proteins, compile for each protein is GPI and transmembrane segment annotations
    has_GPI = defaultdict(lambda: False)
    TM_segments = defaultdict(set)
    protein_seqs = {}

    getAllFromHMPAS(HMPAS_ClassMember_file, HMPAS_ProteinInfo_file, HMPAS_Sequence_file, HMPAS_SequenceFeature_file, has_GPI, TM_segments, protein_seqs)
    getAllFromUniProt(Uniprot_PM_transmembrane_proteins, Uniprot_PM_GPI_proteins, has_GPI, TM_segments, protein_seqs)
    print >> sys.stderr, "WARNING: skipping HPA for now. HPA plasma membrane proteins not necessarily exposed extracellularly."
    getAllFromHPA(HPA_PM_proteins, has_GPI, TM_segments, protein_seqs)

    uniprot_from_hmpas = set(new_hmpas.keys()) | set(hmpas_gpi.keys())
    all_hpa = getAllHPA(HPA_PM_proteins)
    new_uniprot = getAllUniProt(Uniprot_PM_proteins, Uniprot_all_proteins, all_hpa, uniprot_from_hmpas):

    writeOutput(has_GPI, TM_segments, protein_seqs, output_tsv)

    assessHPACoverageByUniprot(all_hpa, all_uniprot)

    # Record which regions of each HMPAS protein aligned to which pORFs
    match_data = collectDataForMatchesToHMPAS(CG_v_HMPAS_align_gz)
    #cPickle.dump(match_data, open("CG_v_HMPAS_match_data.pkl", 'wb'))
    #match_data = cPickle.load(open("CG_v_HMPAS_match_data.pkl", 'rb'))
    
    # Determine which pORFs aligned to TM regions of HMPAS proteins
    pORFs_aligned_to_TM_segments = identifyTMpORFs(match_data, hmpas_TM)

    sys.exit(0)
