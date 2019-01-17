SHELL=/bin/bash

VERSION = $(shell basename `pwd`)

CGDB = /raid1/projects/CGDB
THIS_DIR = ${CGDB}/pORFs/${VERSION}
SCRIPTS = ${CGDB}/scripts
CG_FASTA = ${CGDB}/models/CGDB${VERSION}.fa
CG_GTF = ${CGDB}/models/CGDB${VERSION}.gtf

HPA = /raid1/databases/HPA
HMPAS = /raid1/databases/HMPAS
UNIPROT = /raid1/databases/UniProt

START_CODONS = atg ctg gtg ttg acg # As cDNA. Implied initiator amino acids: M L V L T
#ALL_PORF_FILES = $(shell ls ${THIS_DIR}/CG-{atg,ctg,gtg,ttg,acg,xxx}-pORFs-??50.fa.bz2)
GE_PORF_FILES = $(shell ls CG-{atg,ctg,gtg,ttg,acg,xxx}-pORFs-ge25.fa.bz2)
#PREDGPI_FILES = $(shell ls ${THIS_DIR}/CG-{atg,ctg,gtg,ttg,acg}-pORFs-ge50.predGPI.bz2)
#SIGNALP_FILES = $(shell ls ${THIS_DIR}/CG-{atg,ctg,gtg,ttg,acg}-pORFs-ge50.signalp.bz2)


.PHONY: all_pORF_files annotated-PM-v-longest_%_similarity
.SECONDARY:


all_pORF_files: ${CG_FASTA}
	${SCRIPTS}/mRNAAllORFs.py 25 $<


same-pORF_%_isoform_groups.tsv.bz2:
	${SCRIPTS}/compileIsoformGroups <(bzcat $(GE_PORF_FILES) | awk '{ l=$$0; getline; printf "%s\t%s\n",$$0,substr(l,2) }') \
		${CGDB}/models/CGDB${VERSION}_loci_names.tsv same-pORF_$*_isoform_groups.tsv
	bzip2 same-pORF_$*_isoform_groups.tsv


same-pORF_%_isoform_groups-blastdb: same-pORF_%_isoform_groups.tsv.bz2
	-mkdir blast_index
	bzcat $< | tail -n +2 | awk '{printf ">%s\n%s\n",$$1,$$5}' > same-pORF_$*_isoform_groups.fa
	makeblastdb -in same-pORF_$*_isoform_groups.fa -dbtype prot -out blast_index/same-pORF_$*_isoform_groups -title same-pORF_$*_isoform_groups
	rm same-pORF_$*_isoform_groups.fa


annotated_PM_proteins.tsv.gz: ${HPA}/cell_membrane_proteins.tsv.gz ${UNIPROT}/cell_membrane_proteins.tsv.gz ${UNIPROT}/all_membrane_proteins.tsv.gz
	${SCRIPTS}/compileMasterAnnotationSource.py $^ ${HMPAS}/HMPAS_Sequence.txt ${HMPAS}/HMPAS_ClassMember.txt ${HMPAS}/HMPAS_SequenceFeature.txt \
	${HMPAS}/HMPAS_ProteinInfo.txt $@


# This should use the blast db made from the same-pORFs
# ||%saa length($$2),
annotated-PM-v-%_isoform_groups: annotated_PM_proteins.tsv.gz
	-mkdir $@
	${SCRIPTS}/sequenceSimilaritySearch.py 2 <(zcat $< | tail -n +2 | awk '{printf ">%s\n%s\n",$$1,$$2}') blast_index/same-pORF_$*_isoform_groups $@


pORF_%_PM_annotations.tsv.gz: annotated_PM_proteins.tsv.gz
	${SCRIPTS}/annotatepORFs.py $< annotated-PM-v-$*_isoform_groups $@


# Next to do: 
# % is, for example, ge25_PM
same-pORF_ge25_PM_isoform_groups.tsv.bz2: same-pORF_ge25_isoform_groups.tsv.bz2 pORF_ge25_PM_annotations.tsv.gz
	#t=`echo $* | cut -f 2 -d_`; \#
	${SCRIPTS}/selectIsoformGroupsOfType.py PM $^ $@


# Needed?
CG-longest-pORFs-ge25.fa.bz2: $(addsuffix -pORFs-ge25.fa.bz2, $(addprefix CG-, $(START_CODONS)))
	bzcat $^ | ${SCRIPTS}/selectLongestORFForEachStopCodon.py $@



