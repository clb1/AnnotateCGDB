SHELL = /bin/bash


all_udb:
	for x in `ls CGDB_*.fa`; do \
		t=`echo $$x | sed 's/fa$$/udb/'` ; \
		${MAKE} $$t ; \
	done


CGDB_%.udb: CGDB_%.fa
	usearch -makeudb_usearch $< -output $@ -wordlength 8 -dbstep 1 -alpha nt
