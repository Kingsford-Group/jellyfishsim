bin_PROGRAMS = jellyfishsim

AM_CXXFLAGS = $(JELLYFISH_CFLAGS) -I$(top_srcdir)/gzstream
AM_LDFLAGS = $(JELLYFISH_LDFLAGS)
LDADD = $(JELLYFISH_LIBS) -lz

jellyfishsim_SOURCES = jellyfishsim.cc gzstream/gzstream.C	\
                       gzstream/gzstream.h
