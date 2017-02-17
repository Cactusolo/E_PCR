/* $Id: mswin.h,v 1.6 2007/07/05 16:05:58 rotmistr Exp $
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
#ifndef EPCR_MSWIN__HPP
#define EPCR_MSWIN__HPP

#include <windows.h>
#include <io.h>
#include <stdio.h>

#define NO_POPEN 1
#define USE_WIN  1
#define O_LARGEFILE 0

typedef char * caddr_t;
typedef long int32_t;
typedef __int64 off64_t;
typedef __int64 huge;
typedef unsigned __int64 uhuge;

#define fopen64 fopen
#define snprintf _snprintf

inline off64_t lseek64(int fd, off64_t off, int dir) 
{
	long lo=off;
	long hi=off>>32;
	
	switch(dir) {
	case SEEK_SET: dir=FILE_BEGIN; break;
	case SEEK_CUR: dir=FILE_CURRENT; break;
	case SEEK_END: dir=FILE_END; break;
	}
	HANDLE h=(HANDLE)_get_osfhandle(fd);
	lo=SetFilePointer(h, lo, &hi, dir);
	if(lo==INVALID_SET_FILE_POINTER && GetLastError()!=NO_ERROR) return -1;
	return lo+((off64_t(hi)<<32)&0xffffffff);
}

inline int fseeko64(FILE* f, off64_t off, int dir)
{
	fflush(f);
	off64_t rc = lseek64(fileno(f), off, dir);
	return rc == (off64_t)-1 ? -1 : 0;
} 
inline off64_t ftello64(FILE* f)
{
	fflush(f);
	return lseek64(fileno(f), 0, SEEK_CUR);
} 

#ifdef __cplusplus
extern "C" {
#endif

    extern int optind;
    extern int optopt;
    extern int opterr;
    extern const char* optarg;

    int getopt(int argc, char ** argv, const char* optstring);
    
#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: mswin.h,v $
 * Revision 1.6  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.5  2004/09/03 19:10:21  rotmistr
 * Public domain notice added.
 *
 */
