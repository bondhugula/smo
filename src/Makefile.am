# SMO
#
# Copyright (C) Somashekaracharya G B 2014
#

SUBDIRS = 

MAINTAINERCLEANFILES = Makefile.in


if EXTERNAL_ISL
# ISL_INCLUDE, ISL_LIBADD were already defined
else
ISL_INCLUDE += -I$(top_srcdir)/isl/include
ISL_LIBADD += $(top_srcdir)/isl/libisl.la
endif

bin_PROGRAMS = smo 

smo_SOURCES = main.cpp ConflictSpec.cpp ConflictSpec.h ConflictSpecBuilder.cpp ConflictSpecBuilder.h ISLUtils.cpp ISLUtils.h IterativeStoragePartition.cpp IterativeStoragePartition.h Namer.cpp Namer.h StoragePartition.cpp StoragePartition.h SMO.h
smo_CPPFLAGS = -O3 -funroll-loops -fopenmp \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../pet/include \
   -I../include \
   $(ISL_INCLUDE)

smo_LDFLAGS = -static -fopenmp
smo_LDADD =  \
   ../pet/libpet.la \
   $(ISL_LIBADD) \
   -lm \
   -lglpk

lib_LTLIBRARIES = libsmo.la

# libsmo doesn't use polylib-based functions
libsmo_la_SOURCES = main.cpp ConflictSpec.cpp ConflictSpec.h ConflictSpecBuilder.cpp ConflictSpecBuilder.h ISLUtils.cpp ISLUtils.h IterativeStoragePartition.cpp IterativeStoragePartition.h Namer.cpp Namer.h StoragePartition.cpp StoragePartition.h SMO.h
libsmo_la_CPPFLAGS = -O3 -funroll-loops -fopenmp \
   -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../pet/include \
   $(ISL_INCLUDE)

# libsmo doesn't use polylib
libsmo_la_LIBADD = \
   ../pet/libpet.la \
   $(ISL_LIBADD) \
   -lm -lgomp
