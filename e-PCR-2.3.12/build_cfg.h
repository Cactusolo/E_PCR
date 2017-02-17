/* $Id: build_cfg.h,v 1.13 2008/03/26 16:04:29 rotmistr Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * =========================================================================
 *
 * Author: Kirill Rotmistrovsky
 *
 * ========================================================================= */

#ifndef EPCR_BUILD__HPP
#define EPCR_BUILD__HPP

#ifdef STANDALONE

#include <sys/types.h>

namespace std {}

#define BEGIN_SCOPE(a) namespace a {
#define END_SCOPE(a)   }

#define BEGIN_NCBI_SCOPE BEGIN_SCOPE(ncbi) USING_SCOPE(std);
#define END_NCBI_SCOPE   END_SCOPE(ncbi)

#define USING_SCOPE(a) using namespace a
#define USING_NCBI_SCOPE USING_SCOPE(ncbi)

#ifdef _WIN32
#include <epcr/mswin.h>
#define FILE_BINARY "b"
#define FILE_TEXT   "t"
//#warning "Using Borland C/C++ Builder config"
#define madvise(a,b,c) // no madvise
#define MADV_SEQUENTIAL 0
#define MADV_DONTNEED 0

#ifdef __cplusplus
BEGIN_NCBI_SCOPE
#endif

typedef char Int1;
typedef short Int2;
typedef int   Int4;
typedef long long Int8;

typedef unsigned char Uint1;
typedef unsigned short Uint2;
typedef unsigned int   Uint4;
typedef unsigned long long Uint8;

#ifdef __cplusplus
END_NCBI_SCOPE
#endif

#else // _WIN32

#include <unistd.h>
#include <sys/mman.h>
#include <sys/fcntl.h>

#ifdef NATIVE_LARGEFILES
#include <epcr/native64.h>
#endif // NATIVE_LARGEFILES

#include <stdint.h>

#ifdef __cplusplus
BEGIN_NCBI_SCOPE
#endif

typedef   int8_t  Int1;
typedef  int16_t  Int2;
typedef  int32_t  Int4;
typedef  int64_t  Int8;

typedef  uint8_t Uint1;
typedef uint16_t Uint2;
typedef uint32_t Uint4;
typedef uint64_t Uint8;

#ifdef __cplusplus
END_NCBI_SCOPE
#endif

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
    int madvise(void* addr, size_t len, int advice);
#ifdef __cplusplus
}
#endif // __cplusplus


#endif // _WIN32

#else // STANDALONE

#include <corelib/ncbistd.hpp>

#include <unistd.h>
#include <sys/mman.h>
#include <sys/fcntl.h>

#endif // STANDALONE

#define EPCR_SCOPE pcr_tools

#ifndef FILE_BINARY 
#define FILE_BINARY 
#endif

#ifndef FILE_TEXT
#define FILE_TEXT
#endif   

#endif

/*
 * $Log: build_cfg.h,v $
 * Revision 1.13  2008/03/26 16:04:29  rotmistr
 * Added support for blastdb files
 *
 * Revision 1.12  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.11  2007/07/05 16:23:08  rotmistr
 * Forgot two changes
 *
 * Revision 1.10  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.9  2004/09/03 21:28:49  rotmistr
 * Fixes to compile with Borland C++ 5.5
 *
 * Revision 1.8  2004/09/03 15:54:43  rotmistr
 * Compilation for Mac OS/X
 *
 * Revision 1.7  2004/05/27 20:35:46  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.6  2004/04/15 14:18:22  rotmistr
 * Fix to compile with NCBI toolkit (CGI)
 *
 * Revision 1.5  2004/04/06 04:53:17  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.4  2004/04/01 16:37:40  rotmistr
 * Cleaned after adding windows capabilities
 *
 * Revision 1.3  2004/04/01 05:57:52  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.2  2004/02/04 21:23:21  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
