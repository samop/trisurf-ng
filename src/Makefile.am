bin_PROGRAMS = trisurf tsmeasure
trisurf_SOURCES = general.c vertex.c bond.c triangle.c cell.c vesicle.c initial_distribution.c io.c dumpstate.c frame.c energy.c timestep.c vertexmove.c bondflip.c main.c poly.c stats.c sh.c shcomplex.c constvol.c snapshot.c restore.c
GITVERSION:=$(shell git --no-pager describe --tags --always --dirty)
AM_CFLAGS = -Wall -Werror -DTS_VERSION=\"$(GITVERSION)\" -fgnu89-inline
AM_CPPFLAGS = ${libxml2_CFLAGS} -fgnu89-inline
trisurf_LDADD = ${libcurl_LIBS} ${libxml2_LIBS}
#shdiscoverdir=../
#shdiscover_PROGRAMS= shdiscover
#shdiscover_SOURCES= general.c vertex.c bond.c triangle.c cell.c vesicle.c initial_distribution.c io.c energy.c sh.c shdiscover.c poly.c stats.c shcomplex.c
#co_testdir=../
#co_test_PROGRAMS=co_test
#co_test_SOURCES= general.c vertex.c bond.c triangle.c cell.c vesicle.c initial_distribution.c io.c energy.c sh.c co_test.c frame.c bondflip.c poly.c stats.c shcomplex.c
#spherical_trisurfdir=../
#spherical_trisurf_PROGRAMS = spherical_trisurf
#spherical_trisurf_SOURCES = general.c vertex.c bond.c triangle.c cell.c vesicle.c initial_distribution.c io.c frame.c energy.c timestep.c vertexmove.c spherical_trisurf.c sh.c bondflip.c poly.c stats.c shcomplex.c
#spherical_trisurf_ffdir=../
#spherical_trisurf_ff_PROGRAMS = spherical_trisurf_ff
#spherical_trisurf_ff_SOURCES = general.c vertex.c bond.c triangle.c cell.c vesicle.c initial_distribution.c io.c frame.c energy.c timestep.c vertexmove.c spherical_trisurf_ff.c sh.c bondflip.c poly.c stats.c shcomplex.c


tsmeasure_SOURCES = general.c vertex.c bond.c triangle.c cell.c vesicle.c initial_distribution.c io.c dumpstate.c frame.c energy.c timestep.c vertexmove.c bondflip.c tsmeasure.c poly.c stats.c sh.c shcomplex.c constvol.c snapshot.c restore.c
tsmeasure_LDADD = ${libcurl_LIBS} ${libxml2_LIBS}
#gitversion.c: .git/HEAD .git/index
#    echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" > $@
