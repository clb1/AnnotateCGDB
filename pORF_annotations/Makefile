SHELL=/bin/bash

CGDB = /raid1/projects/CGDB
THIS_DIR = ${CGDB}/pORF_annotations
PORFS_DIR = ${CGDB}/pORFs
SCRIPTS = ${CGDB}/scripts
#CG_FASTA=${CGDB}/models/CG.fa
#CG_GTF=${CGDB}/models/CG.gtf
#UNIREF=/raid/databases/UniRef/uniref80.fa
#BLAST_BIN=/raid/software/src/ncbi-blast-2.2.26+/bin

HMPAS = /raid1/databases/HMPAS

#START_CODONS = atg ctg gtg ttg acg # As cDNA. Implied initiator amino acids: M L V L T
#ALL_PORF_FILES = $(shell ls ${THIS_DIR}/CG-{atg,ctg,gtg,ttg,acg,xxx}-pORFs-??50.fa.bz2)
#GE_PORF_FILES = $(shell ls ${THIS_DIR}/CG-{atg,ctg,gtg,ttg,acg,xxx}-pORFs-ge50.fa.bz2)
#PREDGPI_FILES = $(shell ls ${THIS_DIR}/CG-{atg,ctg,gtg,ttg,acg}-pORFs-ge50.predGPI.bz2)
#SIGNALP_FILES = $(shell ls ${THIS_DIR}/CG-{atg,ctg,gtg,ttg,acg}-pORFs-ge50.signalp.bz2)

.PHONY: 
.SECONDARY:


CG-longest-pORFs-ge25_v_HMPAS.aln.gz: ${PORFS_DIR}/CG-longest-pORFs-ge25.fa ${HMPAS}/HMPAS_Sequence.fa
	usearch -threads 16 -usearch_local $< -db ${HMPAS}/HMPAS_Sequence.fa -id 0.5 -alnout CG-longest-pORFs-ge25_v_HMPAS.aln
	gzip CG-longest-pORFs-ge25_v_HMPAS.aln


pORF_HMPAS_annotations.tsv: CG-longest-pORFs-ge25_v_HMPAS.aln.gz ${HMPAS}/HMPAS_Sequence.txt \
			${HMPAS}/HMPAS_ClassMember.txt ${HMPAS}/HMPAS_SequenceFeature.txt ${HMPAS}/HMPAS_ProteinInfo.txt
	${SCRIPTS}/annotatepORFsWithHMPAS.py $^ $@


PM_pORFs.txt: pORF_HMPAS_annotations.tsv
	tail -n +2 $< | cut -f 1 -d, | sort | uniq > $@


HMPAS-v-longest_similarity: ${HMPAS}/HMPAS_Sequence.txt
	-mkdir HMPAS-v-longest_similarity
	${SCRIPTS}/sequenceSimilaritySearch.py 8 <(awk '{printf ">%s||%saa\n%s\n",$$1,$$2,$$3}' $<) HMPAS-v-longest_similarity


%.predGPI.bz2: %.fa.bz2
	${SCRIPTS}/predGPI.py <(bzcat $<) $(basename $@)
	bzip2 $(basename $@)


pORF_TM_and_GPI_annotations.tsv: pORF_HMPAS_annotations.tsv CG-longest-pORFs-ge50.fa
	egrep "^>" CG-longest-pORFs-ge50.fa | cut -b 2- | sort -k1,1 > longest.tsv
	bzcat $(PREDGPI_FILES) | egrep -v "^#" | awk '{if ($$2<=0.001) printf "%s\t%s\n",$$1,$$3}' | sort -k1,1 > intermed_tmp.tsv
	cat <(echo "pORF_ID     PREDGPI") <(join -t $$'\t' longest.tsv intermed_tmp.tsv)  > longest_predgpi_intermed.tsv
	${SCRIPTS}/mergeTMandGPIAnnotations.py $< longest_predgpi_intermed.tsv $@
	rm longest.tsv intermed_tmp.tsv longest_predgpi_intermed.tsv


all-signalp-pORFs-%:
	for x in ${START_CODONS}; do \
		${MAKE} CG-$$x-pORFs-ge$*.signalp.bz2 ; \
	done


%.signalp.bz2: %.fa.bz2
	bzcat $< | sed s/^[LVT]/M/ | /raid/software/src/signalp-4.1/signalp -t euk -u sensitive -s best -f short | bzip2 --stdout > $@


proteins_w_signal_peptides.txt: ${SIGNALP_FILES}
	:> $@
	for f in $^; do \
		bzcat $$f | egrep "^CG" | awk '($$9>=0.45)' >> $@ ; \
	done

