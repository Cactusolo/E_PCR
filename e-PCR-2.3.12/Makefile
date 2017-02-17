## $Id: Makefile,v 1.7 2007/07/05 16:06:04 rotmistr Exp $
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

LIBS = epcr 
BINS = cmd_epcr cmd_famap cmd_fahash cmd_rpcr
srcdir = .

all links dirs clean dist clean-all install install-lib dist-clean depend: 	
	for i in $(LIBS:%=lib%) $(BINS) ; do \
		$(MAKE) -ef $(srcdir)/stand/Makefile.$$i $@ ; \
	done

#

########################################################################
## $Log: Makefile,v $
## Revision 1.7  2007/07/05 16:06:04  rotmistr
## Made things compileable by MS Visual C++ 8.0
##
## Revision 1.6  2004/05/27 20:36:03  rotmistr
## Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
##
## Revision 1.5  2004/04/06 16:44:57  rotmistr
## *** empty log message ***
##
## Revision 1.4  2004/03/29 03:16:47  rotmistr
## *** empty log message ***
##
## Revision 1.3  2004/02/04 21:23:46  rotmistr
## - gcc-3.3.2 compatible
## - better postfiltering for reverse-e-PCR for discontiguos words
## - cgi added, that supports:
##  -- contig to chromosome mapping
##  -- simple mapviewer links
##  -- unists links
##  -- discontiguos words
##
## Revision 1.2  2003/12/23 21:30:57  rotmistr
## - gaps/mismatches reporting
## - lo/hi fixup
## - reverse sts in re-PCR_main
##
## Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
## Package that includes e-PCR, reverse e-PCR, and sequence data preparation
## program for reverse e-PCR looks ready
##
## Revision 1.3  2003/11/20 05:56:02  rotmistr
## Loading looks working
##
## Revision 1.2  2003/11/20 02:12:28  rotmistr
## Fixed id, log tags and copyright notice
##
########################################################################
