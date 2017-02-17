## $Id: version.mk,v 1.26 2008/06/18 14:48:20 rotmistr Exp $
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
#set -x
#echo "HELLO!!!"
VER_MAJOR=2
VER_MINOR=3
VER_BUILD=12
VERSION = $(VER_MAJOR).$(VER_MINOR).$(VER_BUILD)
########################################################################
## $Log: version.mk,v $
## Revision 1.26  2008/06/18 14:48:20  rotmistr
## *** empty log message ***
##
## Revision 1.25  2008/04/28 16:39:19  rotmistr
## Applied patch to build with gcc-4.3
##
## Revision 1.24  2008/03/26 16:04:35  rotmistr
## Added support for blastdb files
##
## Revision 1.23  2007/07/11 20:49:33  rotmistr
## Made 64bit-compatible
##
## Revision 1.22  2007/07/05 16:06:04  rotmistr
## Made things compileable by MS Visual C++ 8.0
##
## Revision 1.21  2005/06/14 16:46:51  rotmistr
## Changed report format for floppy tails
##
## Revision 1.20  2005/02/11 20:42:59  rotmistr
## Fixed "margin" bug, added primer search from file
##
## Revision 1.19  2004/10/26 17:16:41  rotmistr
## Added 5'-end masking for primers
##
## Revision 1.18  2004/09/03 15:54:48  rotmistr
## Compilation for Mac OS/X
##
## Revision 1.17  2004/06/08 16:14:59  rotmistr
## *** empty log message ***
##
## Revision 1.16  2004/06/07 16:25:03  rotmistr
## Bug fixes to previos version.
##
## Revision 1.15  2004/06/03 23:37:29  rotmistr
## New aligner added.
##
## Revision 1.14  2004/04/27 00:01:55  rotmistr
## Second version of reverse hash file started
##
## Revision 1.13  2004/04/06 04:53:18  rotmistr
## All is compileable with BCC5.5 and runnable on WIndows
##
## Revision 1.12  2004/04/01 17:24:46  rotmistr
## *** empty log message ***
##
## Revision 1.11  2004/03/30 19:08:08  rotmistr
## default STS size is tunnable now
##
## Revision 1.10  2004/03/26 17:02:18  rotmistr
## Compat-options are now allowed everywhere, and multiple fasta files can be used.
##
## Revision 1.9  2004/03/25 19:36:52  rotmistr
## API: separate left and right primers mism/gaps in forward API
##
## Revision 1.8  2004/03/23 22:36:02  rotmistr
## 2.0 release
##
## Revision 1.7  2004/02/18 05:44:40  rotmistr
## Changes in CGI: sort order, separate misalignments for l and r primers, reload button
##
## Revision 1.6  2004/02/12 21:39:29  rotmistr
## New version
##
## Revision 1.5  2004/01/28 23:27:09  rotmistr
## "Best of overlapping" hit selection postprocessor added.
##
## Revision 1.4  2004/01/08 23:22:47  rotmistr
## Fixed init error in faread,
## Adjusted output to standard,
## Added output format style and output file to parameters.
##
## Revision 1.3  2004/01/07 16:57:48  rotmistr
## Fragment size is now configurable.
##
## Revision 1.2  2004/01/06 21:54:28  rotmistr
## Statistics for word repetitions API added
##
## Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
## Package that includes e-PCR, reverse e-PCR, and sequence data preparation
## program for reverse e-PCR looks ready
##
## Revision 1.2  2003/11/20 02:12:32  rotmistr
## Fixed id, log tags and copyright notice
##
########################################################################
