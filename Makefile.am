#
# Makefile.am
#
# Copyright (C) 2007 Uday Bondhugula 
# 
#

RM = rm -f
LN = ln -s

if EXTERNAL_ISL
MAY_ISL =
# ISL_INCLUDE and ISL_LIBADD were already defined
else
MAY_ISL = isl
ISL_INCLUDE += -I$(top_srcdir)/isl/include
ISL_LIBADD += $(top_srcdir)/isl/libisl.la
endif


#SUBDIRS = piplib openscop clan candl polylib $(MAY_ISL) cloog-isl src
SUBDIRS = $(MAY_ISL) pet src

ACLOCAL_AMFLAGS = -I m4

dist-hook:
	rm -rf `find $(distdir)/doc -name CVS`


#bin_SCRIPTS = polycc vloog ploog plorc plutune getversion.sh
bin_SCRIPTS = getversion.sh
CLEANFILES: $(bin_SCRIPTS) parsetab.py
EXTRA_DIST = polycc.sh.in  examples test


pkginclude_HEADERS = include/smo/libsmo.h

#polycc: polycc.sh
#	rm -f polycc
#	echo "#! " $(BASH) > polycc
#	cat $(srcdir)/polycc.sh >> polycc
#	chmod ugo+x polycc

.PHONY: bin binit

bindir:
	@if [ ! -d bin ];then \
		echo mkdir bin; \
		mkdir bin; \
		fi

binit:
	cp -f src/smo bin/
#	cp -r polycc bin/

pclean:
	$(MAKE) -C src clean

test_libsmo: src/test_libsmo.c
	libtool --mode=link gcc -g src/test_libsmo.c -Lsrc/.libs/ $(ISL_INCLUDE) $(ISL_LIBADD) -I include -lsmo -o test_libsmo

test: .PHO

.PHO: 
	./test.sh
	./test.sh --tile --parallel --unroll
	./test.sh --maxfuse
	./test.sh --nofuse
