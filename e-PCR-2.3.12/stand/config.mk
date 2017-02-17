## $Id: config.mk,v 1.7 2007/07/05 16:06:04 rotmistr Exp $
########################################################################
##
##                            PUBLIC DOMAIN NOTICE
##               National Center for Biotechnology Information
##
##  This software/database is a "United States Government Work" under the
##  terms of the United States Copyright Act.  It was written as part of
##  the author's official duties as a United States Government employee and
##  thus cannot be copyrighted.  This software/database is freely available
##  to the public for use. The National Library of Medicine and the U.S.
##  Government have not placed any restriction on its use or reproduction.
##
##  Although all reasonable efforts have been taken to ensure the accuracy
##  and reliability of the software and data, the NLM and the U.S.
##  Government do not and cannot warrant the performance or results that
##  may be obtained by using this software or data. The NLM and the U.S.
##  Government disclaim all warranties, express or implied, including
##  warranties of performance, merchantability or fitness for any particular
##  purpose.
##
##  Please cite the author in any work or product based on this material.
##
########################################################################

tgtdir = .
objdir = .
srcdir = .
COMMON_CC_FLAGS = 
ifdef OPTIMIZE
CC_FLAGS = $(COMMON_CC_FLAGS) -O$(OPTIMIZE)
LD_FLAGS = 
else
ifdef PROFILE
CC_FLAGS = $(COMMON_CC_FLAGS) -g -pg
LD_FLAGS = -g -pg
else
CC_FLAGS = $(COMMON_CC_FLAGS) -g2 
LD_FLAGS = -g2 
endif
endif

#arch   = $(shell echo `uname -s`-`uname -m`)
prefix = /usr/local/
#$(arch)
BINDIR = $(prefix)/bin/
INCDIR = $(prefix)/include/
LIBDIR = $(prefix)/lib/

FIXCMD = perl -ne's/^([^\s\#]+)/\$$(objdir)\/$$1/;print'
include $(srcdir)/stand/version.mk

#########################################################################
# GNU compiler flags
CC = gcc
CXX = g++
CXXFLAGS = -I$(srcdir) -I$(INCDIR) $(CC_FLAGS) $(PART_CXXFLAGS) \
	-DDEALLOCATE=0 $(LF64CCFLAGS) $(VERSION_FLAGS) -DSTANDALONE=1
LDFLAGS = $(LD_FLAGS) $(LF64LDFLAGS) -L$(tgtdir) -L$(LIBDIR) $(PART_LDFLAGS) 
#	$(PART_PRELIBS) $(LIBS:%=-l%) $(PART_POSTLIBS)

LF64CCFLAGS = `getconf LFS_CFLAGS` 
LF64LDFLAGS = `getconf LFS_LDFLAGS` `getconf LFS_LIBS` 

## Use following lines if you don't have getconf but need to 
## explicitely turn on largefile support 
# LF64CCFLAGS = -D_LARGEFILE64_SOURCE -DFILE_OFFSET_BITS=64 
# LF64LDFLAGS =

## Use following lines for Mac OS X and other systems that lack *64 functions
# LF64CCFLAGS = -DNATIVE_LARGEFILE 
# LF64LDFLAGS =

VERSION_FLAGS = -DVERSION=\"$(VERSION)\" \
	-DVER_MAJOR=$(VER_MAJOR) \
	-DVER_MINOR=$(VER_MINOR) \
	-DVER_BUILD=$(VER_BUILD) 

LIBS = seq epcr

src = $(SRC:%=$(srcdir)/%)
hdr = $(HDR:%=$(srcdir)/%)

all: links target

links:
	if test -n "$(LIBNAME)" ; then \
		test -L $(LIBNAME) || ln -s $(srcdir) $(LIBNAME) ; \
	fi

dirs: 
	for i in $(INCDIR)/$(LIBNAME) $(BINDIR) $(LIBDIR) ; do \
		test -d $$i || mkdir -p $$i ; \
	done

clean: 
	-rm $(OBJ) $(HDR:%=%~) $(SRC:%=%~)

clean-all: clean
	-rm $(TARGET) 

dist-clean: clean-all
	-rm *~
	-test -L $(LIBNAME) && rm $(LIBNAME)

$(objdir)/%.o: $(srcdir)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

########################################################################
## $Log: config.mk,v $
## Revision 1.7  2007/07/05 16:06:04  rotmistr
## Made things compileable by MS Visual C++ 8.0
##
## Revision 1.6  2004/09/08 18:30:59  rotmistr
## Fixed typo
##
## Revision 1.5  2004/09/03 15:54:48  rotmistr
## Compilation for Mac OS/X
##
## Revision 1.4  2004/06/03 23:37:29  rotmistr
## New aligner added.
##
## Revision 1.3  2004/03/31 05:04:11  rotmistr
## Search range fix
##
## Revision 1.2  2004/02/04 21:23:46  rotmistr
## - gcc-3.3.2 compatible
## - better postfiltering for reverse-e-PCR for discontiguos words
## - cgi added, that supports:
##  -- contig to chromosome mapping
##  -- simple mapviewer links
##  -- unists links
##  -- discontiguos words
##
## Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
## Package that includes e-PCR, reverse e-PCR, and sequence data preparation
## program for reverse e-PCR looks ready
##
## Revision 1.8  2003/12/10 19:55:48  rotmistr
## Plain fasta interface is about to be substituted to blastdb interface
##
## Revision 1.7  2003/12/04 21:29:34  rotmistr
## Looks like faindex branch works better!
##
## Revision 1.6  2003/11/24 19:33:40  rotmistr
## Optimised. Added OneTimeRun flag.
##
## Revision 1.5  2003/11/23 03:40:53  rotmistr
## Looks like working, requires optimisation.
##
## Revision 1.4  2003/11/20 23:05:58  rotmistr
## Contiguos words work.
## Discontiguos need to be modified.
##
## Revision 1.3  2003/11/20 18:27:32  rotmistr
## Sample files updated
## Program does not crush
##
## Revision 1.2  2003/11/20 02:12:28  rotmistr
## Fixed id, log tags and copyright notice
##
########################################################################
