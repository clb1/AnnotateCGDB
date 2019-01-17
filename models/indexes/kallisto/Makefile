SHELL = /bin/bash

TRANSCRIPTOME_FASTA = /raid1/projects/CGDB/models/CGDBv1.0.fa

KALLISTO = /usr/local/src/kallisto_linux-v0.43.0/kallisto

.PHONY: kallisto_index

kallisto_index: ${TRANSCRIPTOME_FASTA}
	${KALLISTO} index -i CGDBv1.0.index ${TRANSCRIPTOME_FASTA}
