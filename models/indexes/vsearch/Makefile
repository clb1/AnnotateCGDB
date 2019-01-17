SHELL = /bin/bash

CGDBv2.0.udb: CGDBv2.0.fa
	vsearch --makeudb_usearch $< --wordlength 7 --dbmask none --qmask none --output $@
