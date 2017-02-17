/* $Id: fast_seqio_read.cpp,v 1.6 2007/07/11 20:49:29 rotmistr Exp $
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

#include <epcr/fast_seqio.hpp>
#include <epcr/bin-io.hpp>

#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);


CFastaMap::~CFastaMap() throw()
{
    Close(); 
}

CFastaMap::CFastaMap():m_Fd(-1)
{}

CFastaMap::CFastaMap(const string& path):m_Fd(-1)
{
	Open(path);
}
	
void CFastaMap::Open(const string& path)
{
	if(m_Fd!=-1) Close();
	int fd=open(path.c_str(),O_RDONLY|O_LARGEFILE);
	if(fd==-1) throw runtime_error(path+": "+strerror(errno));

	m_Fd=fd;
	
	ReadPrologue();
	ReadEpilogue();
}

void CFastaMap::Close()
{
	if(m_Fd!=-1) {
        m_Size.clear();
        m_Offset.clear();
        m_Defline.clear();
        m_Ident.clear();

		close(m_Fd);
	}
}

void CFastaMap::ReadPrologue()
{
	if( sizeof( Int2 ) != 2 ) throw logic_error( "Bad compile options: sizeof( Int2 ) != 2" );
	if( sizeof( Int4 ) != 4 ) throw logic_error( "Bad compile options: sizeof( Int4 ) != 4" );
	if( sizeof( Int8 ) != 8 ) throw logic_error( "Bad compile options: sizeof( Int8 ) != 8" );

	char buff[128];
	if(read(m_Fd,buff,8) != 8) throw runtime_error(strerror(errno));
	if(memcmp(buff,"FASTAMAP",8)) 
		throw runtime_error("FastaMap Signature is not found");
	
	switch(Read<unsigned>(m_Fd)){
	case eHiEndian: m_SwapBytes=true; break;
	case eLoEndian: m_SwapBytes=false; break;
	default: throw runtime_error("Bad format: wrong byteorder signature");
	}
//    if(m_SwapBytes) fprintf(stderr,"* FaMap: Swapping bytes\n");

	short ver=BoCvt(Read<short>(m_Fd),m_SwapBytes);
	short rev=BoCvt(Read<short>(m_Fd),m_SwapBytes);

	if(ver!=1 && rev!=0) throw runtime_error("Wrong FastaMap file version");
}

void CFastaMap::ReadEpilogue()
{
	/*TOffset filesize=*/lseek64(m_Fd,0,SEEK_END);
	TOffset directory_pos=lseek64(m_Fd,-(TOffset)sizeof(TOffset),SEEK_CUR);	
	TOffset directory=BoCvt(Read<TOffset>(m_Fd),m_SwapBytes);

	if(TOffset(directory+3*sizeof(TOffset)) != directory_pos) 
		throw runtime_error("Wrong directory record");
    
	lseek64(m_Fd,directory,SEEK_SET);
	TOffset epilogue=BoCvt(Read<TOffset>(m_Fd),m_SwapBytes);
	TOffset ident   =BoCvt(Read<TOffset>(m_Fd),m_SwapBytes);
	TOffset defline =BoCvt(Read<TOffset>(m_Fd),m_SwapBytes);

	lseek64(m_Fd,epilogue,SEEK_SET);
	char buff[128];
	if(read(m_Fd,buff,8)!=8) throw runtime_error("read failed");
	if(memcmp(buff,"EPILOGUE",8)) 
		throw runtime_error("epilogue is not found!");
	unsigned sz=BoCvt(Read<unsigned>(m_Fd),m_SwapBytes);
	
	m_Size.resize(sz);
	m_Offset.resize(sz);
	m_Defline.resize(sz);
	m_Ident.resize(sz);
	for(unsigned i=0; i<sz; ++i) {
		m_Offset[i]=BoCvt(Read<TOffset>(m_Fd),m_SwapBytes);
	}
	for(unsigned i=0; i<sz; ++i) {
		m_Size[i]=BoCvt(Read<TOffset>(m_Fd),m_SwapBytes);
	}
	lseek64(m_Fd,ident,SEEK_SET);
	for(unsigned i=0; i<sz; ++i)  m_Ident[i]=Read<string>(m_Fd);
	lseek64(m_Fd,defline,SEEK_SET);
	for(unsigned i=0; i<sz; ++i)  m_Defline[i]=Read<string>(m_Fd);
}

CMmSequence::CMmSequence(CFastaMap& fm, unsigned i):
    CIndexAssert(i,fm.GetCount()),
    CMMap(fm.m_Size[i]+1,CMMap::fProtRead,CMMap::fMapPrivate,fm.m_Fd,fm.m_Offset[i])
{ if(m_Size>1) m_Size--; }

CMmSequence::~CMmSequence() throw ()
{ if(m_Size) ++m_Size; }

/*
 * $Log: fast_seqio_read.cpp,v $
 * Revision 1.6  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.5  2004/04/01 16:37:41  rotmistr
 * Cleaned after adding windows capabilities
 *
 * Revision 1.4  2004/04/01 05:57:53  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.3  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.2  2003/12/30 15:27:22  rotmistr
 * Fixed bug with sequence end
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
