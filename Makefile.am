bin_PROGRAMS = jellyfishsim

AM_CXXFLAGS = $(JELLYFISH_CFLAGS) -I$(top_srcdir)/gzstream -fopenmp
AM_LDFLAGS = $(JELLYFISH_LDFLAGS) -fopenmp
LDADD = $(JELLYFISH_LIBS) -lz

jellyfishsim_SOURCES = jellyfishsim.cc gzstream/gzstream.C	\
                       gzstream/gzstream.h
