#!/usr/bin/env python

import cPickle
from itertools import combinations, product
import networkx as nx
import sys
import sqlite3

sys.path.append("/usr/local/src/MFEprimer")
sys.path.append("/usr/local/src/MFEprimer/chilli")
import chilli
from chilli.chilli import int2DNA

import pdb


def mismatch(word, letters, num_mismatches):
    for locs in combinations(range(len(word)), num_mismatches):
        this_word = [[char] for char in word]
        for loc in locs:
            orig_char = word[loc]
            this_word[loc] = [l for l in letters if l != orig_char]
        for poss in product(*this_word):
            yield ''.join(poss)
            

def readKmersAndIndexes(kmer_db, k):
    kmer_to_index = {}

    conn = sqlite3.connect(kmer_db)
    cursor = conn.cursor()
    cursor.execute("SELECT mer_id FROM pos WHERE plus != '' AND minus != ''")
    for kmer_id in cursor.fetchall():
        kmer = int2DNA(kmer_id[0], k).upper()
        kmer_to_index[kmer] = kmer_id

    return kmer_to_index


def connectKmers(kmers_to_index, mm_limit):
    G = nx.Graph()

    for kmer, kmer_index in kmers_to_index.items():
        for mm_kmer in mismatch(kmer, "ACGT", mm_limit):
            G.add_edge(kmer_index, kmers_to_index[mm_kmer])

    return G


if (__name__ == "__main__"):
    mm_limit, k, kmer_db, output_graph_pkl = sys.argv[1:]
    mm_limit = int(mm_limit)
    k = int(k)
    
    kmers_to_index = readKmersAndIndexes(kmer_db, k)

    pdb.set_trace()

    G = connectKmers(kmers_to_index, mm_limit)

    cPickle.dump(G, open(output_graph_pkl, 'wb'))
    
    sys.exit(0)
