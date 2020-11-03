bin_PROGRAMS = sortsim hausdorff chunking

AM_CXXFLAGS = $(JELLYFISH_CFLAGS) -I$(top_srcdir)/gzstream -fopenmp
AM_LDFLAGS = $(JELLYFISH_LDFLAGS) -fopenmp
LDADD = $(JELLYFISH_LIBS) -lz


sortsim_SOURCES = sortsim.cc gzstream/gzstream.C gzstream/gzstream.h

hausdorff_SOURCES = hausdorff.cc gzstream/gzstream.C gzstream/gzstream.h

chunking_SOURCES = chunking.cc gzstream/gzstream.C gzstream/gzstream.h
