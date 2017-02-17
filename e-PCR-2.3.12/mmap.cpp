/* $Id: mmap.cpp,v 1.7 2007/07/05 16:05:58 rotmistr Exp $
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

#include <epcr/mmap.hpp>

#include <stdexcept>
//#include <sstream>

#include <string.h>
#include <errno.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

#ifndef USE_WIN

CMMap::CMMap(unsigned size, unsigned prot, unsigned flags,
             TFileHandle fd, TOffset offset) 
{
    static unsigned page=sysconf(_SC_PAGESIZE);
    
    offset-=(m_Delta=(offset%page));
    m_Data=(char*)mmap64(0,(m_Size=size)+m_Delta,prot,flags,fd,offset);
    if(m_Data==0 || m_Data==(char*)-1) {
        throw runtime_error("mmap failed: "+string(strerror(errno)));
    }
    else m_Data+=m_Delta;
}

CMMap::~CMMap() throw()
{
    if(m_Data != (char *)-1) munmap(m_Data-m_Delta,m_Size+m_Delta);
}

#else

static string errmsg() 
{
	char * x;
	if(!FormatMessage(
		   FORMAT_MESSAGE_ALLOCATE_BUFFER|
		   FORMAT_MESSAGE_FROM_SYSTEM|
		   FORMAT_MESSAGE_IGNORE_INSERTS,
		   NULL,
		   GetLastError(),
		   MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		   (char*)&x,0,NULL)) return "UNKNOWN";
	return x?x:"null";
}

CMMap::CMMap(unsigned size, unsigned prot, unsigned flags,
             TFileHandle fd, TOffset offset) 
{
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	
	static unsigned page=si.dwAllocationGranularity;
    
    offset-=(m_Delta=(offset%page));

	m_Size=size;
	m_Map=CreateFileMapping((HANDLE)_get_osfhandle(fd),0,
							prot&fProtWrite?PAGE_READWRITE:PAGE_READONLY,
							0,0,0);
	if(m_Map) {
		m_Data=(char*)MapViewOfFile(
			m_Map,
			flags&fMapPrivate?FILE_MAP_COPY:FILE_MAP_ALL_ACCESS,
			offset>>32,offset,size+m_Delta+1);
	} else m_Data=0;
	
	if(m_Map==0 || m_Data==0 || m_Data==(char*)-1) {
		string err( "mmap failed: " );
		err.append( errmsg() );
        throw runtime_error( err );
    }
    else m_Data+=m_Delta;
}

CMMap::~CMMap() throw()
{
    if(m_Data != (char *)-1 && m_Data != 0) {
		UnmapViewOfFile(m_Data);
		CloseHandle(m_Map);
	}
}

#endif

/*
 * $Log: mmap.cpp,v $
 * Revision 1.7  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.6  2005/01/27 19:09:13  rotmistr
 * Fixed mmap for win32
 *
 * Revision 1.5  2004/04/06 04:53:18  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.4  2004/04/01 16:37:41  rotmistr
 * Cleaned after adding windows capabilities
 *
 * Revision 1.3  2004/04/01 05:57:53  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.2  2004/01/06 21:54:19  rotmistr
 * Statistics for word repetitions API added
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
