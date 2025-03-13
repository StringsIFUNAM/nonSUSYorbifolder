SUBDIRS = \
    src/io \
    src/linalg \
    src/groupTheory \
    src/orbifolder

ACLOCAL_AMFLAGS = -I m4

dist_noinst_DATA = GaugeInvariance/* Geometry/* src/examples/* doc/*

# Add readline library to relevant targets
src_io_LDADD = -lreadline
src_orbifolder_LDADD = -lreadline

# Here include compilation rules
AM_CXXFLAGS = $(CXXFLAGS)
