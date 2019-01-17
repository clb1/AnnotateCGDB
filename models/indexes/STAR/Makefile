SHELL = /bin/bash

TRANSCRIPTOME_FASTA = /raid1/projects/CGDB/models/CG.fa

STAR = /usr/local/src/STAR/bin/Linux_x86_64_static/STAR

.PHONY: transcriptome_index

transcriptome_index: ${TRANSCRIPTOME_FASTA}
	${STAR} --runThreadN 5 --genomeChrBinNbits 13 --runMode genomeGenerate --genomeFastaFiles ${TRANSCRIPTOME_FASTA}
