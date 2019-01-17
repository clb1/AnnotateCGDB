#!/usr/bin/env python

import os
import sys
import tempfile

from Bio import SeqIO

from subprocess import call
from multiprocessing import Process
from multiprocessing.queues import JoinableQueue


def blastQueries(blast_output_dir, query_queue, search_database):
    probable_protein = query_queue.get()

    #userfields = "+".join( ["query", "target", "id", "alnlen", "mism", "opens", "qlor", "qhir", "tlor", "thir", "evalue", "bits"] )
    blastp_fmt = "\"7 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore\""
    while (probable_protein != "STOP"):
        seqid, aa_seq = probable_protein.split()

        blast_output_file = "%s/%s.tblastn" % (blast_output_dir, seqid.split("||")[0])
        name_final_file = "%s.bz2" % blast_output_file
        if (not os.path.exists(name_final_file)):

            # Put the sequence in its own file
            tmp_fasta_file = tempfile.NamedTemporaryFile(mode='w+t',suffix=".fa",bufsize=0)
            tmp_fasta_file.write(">%s\n%s\n" % (seqid, aa_seq))

            # Run BLAST
            blast_cmd = ["tblastn", "-num_threads", "2", "-query", tmp_fasta_file.name, "-db", search_database,
                         "-matrix", "BLOSUM62", "-comp_based_stats", "D", "-outfmt", blastp_fmt, "-out", blast_output_file, "&> /dev/null"]

            # Run USEARCH UBLAST. Need 64-bit USEARCH
            # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
            #blast_cmd = ["usearch", "-ublast", "-threads", "4",
            #             "-query", tmp_fasta_file.name, "-db", search_database,  "-mid", "50.0",
            #             "-userout", blast_output_file, "-userfields", userfields, "&> /dev/null"]
        
            blast_cmd_str = " ".join(blast_cmd)
            #print >> sys.stderr, blast_cmd_str
            err_code = call(blast_cmd_str,shell=True)
            err_code = call("bzip2 %s" % blast_output_file, shell=True)
            tmp_fasta_file.close()
        
        query_queue.task_done()
        probable_protein = query_queue.get()

    query_queue.task_done()

    
if (__name__ == "__main__"):
    num_worker_processes, query_proteins_fasta, search_database, blast_output_dir = sys.argv[1:]
    num_worker_processes = int(num_worker_processes)

    print >> sys.stderr, "INFO: main isoform parts server has PID %d" % os.getpid()

    query_queue = JoinableQueue(200)

    worker_pids = []
    for i in xrange(num_worker_processes):
        p = Process(target=blastQueries, args=(blast_output_dir, query_queue, search_database))
        p.start()
        worker_pids.append(str(p.pid))
    print >> sys.stderr, "INFO: worker PIDs are %s" % ",".join(worker_pids)        

    for seq_record in SeqIO.parse(query_proteins_fasta, "fasta"):
        porf_ID = seq_record.id
        aa_seq = seq_record.seq.tostring().upper()
        query_queue.put( "%s %s" % (porf_ID, aa_seq) )
        
    for j in xrange(num_worker_processes):
        query_queue.put( "STOP" )

    query_queue.join()
    
    sys.exit(0)
