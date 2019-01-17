SHELL = /bin/bash

THIS_DIR = /raid1/projects/CGDB/viz
#DB_MERGING_DIR = ${THIS_DIR}/database_merging

HUMAN_GENOME_FASTA = /raid1/references_and_indexes/hg38/hg38.fa.gz

JBROWSE_ROOT = /srv/www/htdocs/jbrowse/JBrowse-1.11.6
JBROWSE_BIN = ${JBROWSE_ROOT}/bin
JBROWSE_CGDB_DATA_DIR = ${JBROWSE_ROOT}/CGDB


.PHONY: prepare_ref_seq

# Add --compress
prepare_ref_seq:
	${JBROWSE_BIN}/prepare-refseqs.pl --fasta ${HUMAN_GENOME_FASTA} --out ${JBROWSE_CGDB_DATA_DIR}

# Add --compress
prepare_CGDB: ${DB_MERGING_DIR}/CGDB.gff3
	${JBROWSE_BIN}/flatfile-to-json.pl --out ${JBROWSE_CGDB_DATA_DIR} \
	--gff $< \
	--trackType CanvasFeatures \
	--type mRNA,exon \
	--trackLabel CGDB \
	--key CGDB \
	--getSubfeatures \
	--autocomplete all
