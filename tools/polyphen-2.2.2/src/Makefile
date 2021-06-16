# $LastChangedDate: 2012-08-31 14:45:35 -0400 (Fri, 31 Aug 2012) $
# $LastChangedRevision: 403 $
# $LastChangedBy: ivan $

pph2 = $(shell pwd | sed 's/\/src//')

CC = gcc
# static flag may be required on mixed architecture clusters
#CFLAGS = -O3 -static -s
#LFLAGS = -static -s
CFLAGS = -O3 -s
LFLAGS = -s
#CLUOPT = $(LFLAGS) -static -O3
CLUOPT = $(LFLAGS) -O3

all: 	download fetchseq dist bfact psic cluspack leon coils cluspack hbplus dsspcmbi mafft

.PHONY:	fetchseq dist bfact psic leon cluspack coils cluspack hbplus dsspcmbi mafft

fetchseq:	fetchseq/fetchseq.c
	make -C fetchseq CC='$(CC)' CFLAGS='$(CFLAGS)'

dist: 	dist/dist.c
	make -C dist CC='$(CC)' CFLAGS='$(CFLAGS)'

bfact: 	bfact/bfact.c
	make -C bfact CC='$(CC)' CFLAGS='$(CFLAGS)'

psic: 	psic/psic.c
	make -C psic CC='$(CC)' CFLAGS='$(CFLAGS)'

leon:
	(test -d leon) && make -C leon CC='$(CC)' CFLAGS='$(CFLAGS) -c' LFLAGS='$(LFLAGS) -lm'

cluspack:
	(test -d cluspack) && make OPT='$(CLUOPT)' -C cluspack

coils:
	(test -d coils) && cd coils; $(CC) $(CFLAGS) -o ncoils -I. ncoils.c read_matrix.c -lm

hbplus:
	(test -d hbplus) && make -C hbplus CC='$(CC)' CFLAGS='$(CFLAGS)'

dsspcmbi:
	(test -d dsspcmbi) && make -C dsspcmbi CC='$(CC)' CFLAGS='$(CFLAGS) -Wall -DGCC'

mafft:
	(test -d mafft) && make PREFIX=$(pph2) -C mafft/core

install:
	cp -f fetchseq/fetchseq dist/dist bfact/bfact psic/psic \
	cluspack/cluspack \
	coils/ncoils \
	hbplus/hbplus \
	dsspcmbi/dsspcmbi \
	../bin/
	cp -f coils/*.mat ../matr
	cd leon; cp -f resbias convseq bias2relacs query_seqerrs comp_query \
	coils2relacs query2relacs chain_blocks \
	../../bin/; cd ../
	make PREFIX=$(pph2) -C mafft/core install

clean:
	rm -f fetchseq/fetchseq dist/dist bfact/bfact psic/psic
	rm -f leon/*.o
	(test ! -d cluspack) || make -C cluspack clean
	rm -f coils/ncoils
	rm -f hbplus/*.o
	rm -f dsspcmbi/*.o
	(test ! -d mafft/core) || make -C mafft/core clean

download:
	@echo "PLEASE READ LICENSE AGREEMENTS FOR ALL PACKAGES CAREFULLY!"
	(test -f coils.tar.gz || test -d coils) || \
		wget -c http://www.russell.embl-heidelberg.de/coils/coils.tar.gz
	(test -d coils) || tar vxzf coils.tar.gz
# 	(test -f hbplus.tar.Z || test -d hbplus) || \
# 		wget -c ftp://ftp.biochem.ucl.ac.uk/pub/hbplus/hbplus.tar.Z
	(test -d hbplus) || tar vxzf hbplus.tar.Z
	(test -f dsspcmbi.zip || test -d dsspcmbi) || \
 		wget -c ftp://ftp.cmbi.ru.nl/pub/molbio/software/old-dssp/dsspcmbi.zip
	(test -d dsspcmbi) || (unzip -j dsspcmbi.zip -d dsspcmbi && \
		cp -f Makefile.dsspcmbi dsspcmbi/Makefile)
	(test -f mafft-7.221-without-extensions-src.tgz || test -d mafft) || \
		wget -c http://mafft.cbrc.jp/alignment/software/mafft-7.221-without-extensions-src.tgz
	(test -d mafft) || \
		(tar vxzf mafft-7.221-without-extensions-src.tgz && \
		mv mafft-7.221-without-extensions/ mafft)
	(test -f leon_1.0.tar.gz || test -d cluspack) || \
		wget -c ftp://ftp-igbmc.u-strasbg.fr/pub/Leon/leon_1.0.tar.gz
	(test -d cluspack) || \
		(tar vxzf leon_1.0.tar.gz && \
		mv leon/cluspack_src/ cluspack && \
		rm -rf leon && \
		tar vxzf cluspack_patched.tar.gz && \
		make -C cluspack clean)
	(test -f leon_1.1.tar.gz || test -d leon) || \
		wget -c ftp://ftp-igbmc.u-strasbg.fr/pub/Leon/leon_1.1.tar.gz
	(test -d leon) || \
		(mkdir leon && tar vxzf leon_1.1.tar.gz -C leon && tar vxzf leon_1.2.tar.gz -C leon)
	(test -f weka-3-6-12.zip || test -d ../weka) || \
		wget -c http://prdownloads.sourceforge.net/weka/weka-3-6-12.zip
	(test -d ../weka) || \
		(unzip weka-3-6-12.zip && mv weka-3-6-12/ ../weka)
