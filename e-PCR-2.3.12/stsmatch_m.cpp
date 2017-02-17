/* $Id: stsmatch_m.cpp,v 1.16 2008/06/18 14:45:33 rotmistr Exp $
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

#include <epcr/defaults.h>
#include <epcr/align.hpp>
#include <epcr/hashset.hpp>
#include <epcr/stsmatch_m.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <ctype.h>
#include <stdexcept>


USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

CMmFileSts::CMmFileSts(const char * ref, const CStrRef& p1, const CStrRef& p2,
					   unsigned lo, unsigned hi, unsigned flags,
                       unsigned char ovhg1, unsigned char ovhg2):
	m_SizeLo(lo),
	m_SizeHi(hi),
	m_Reference(ref),
	m_Flags(flags)
{
	m_Primer[0]=p1.data();
	m_Length[0]=p1.length();
	m_Primer[1]=p2.data();
	m_Length[1]=p2.length();
    m_OvhgChars[0]=ovhg1;
    m_OvhgChars[1]=ovhg2;
}

void CMmFileSts::ParseRange(const CStrRef& range, int& lo, int& hi, 
                            int defLo, int defHi)
{
	char * c=const_cast<char*>(range.data());
	if(isdigit(*c)) {
		lo=strtol(c,&c,10);
		hi=(*c=='-')?strtol(c+1,0,10):lo;
	}
	else {
		lo=defLo; //ePCR_DEFAULT_size_lo;
		hi=defHi; //ePCR_DEFAULT_size_hi;
	}
}

int CMmFileSts::Parse(const char * ref, CStrRef * dest, int maxf) 
{
	int i=0;
	while(i<maxf) {
		while(*ref && isspace(*ref)) ++ref;
		const char * p=strpbrk(ref,"\t\n\r");
		if(p==0 || *p==0) dest[i++]=CStrRef(ref);
		else {
			const char * pp = p;
			while( pp > ref && isspace(pp[-1]) ) --pp;
			dest[i++]=CStrRef(ref,pp-ref); 
		}
		if(p==0 || *p!='\t') return i;
		ref=p+1;
	}
	return i;
}

CStrRef CMmFileSts::GetName() const 
{
	CStrRef fld("?");
	Parse(m_Reference,&fld,1);
	return fld;
}

CStrRef CMmFileSts::GetDescription() const 
{
	CStrRef fld[5];
	if(Parse(m_Reference,fld,5)>=5) {
		const char * c=fld[4].data();
		const char * cc=strpbrk(c,"\n\r");
		return CStrRef(c,cc?cc-c:strlen(c));
	}
	else return "";
}

////////////////////////////////////////////////////////////////////////

void CStsFileHash::AttachFile(const string& fname)
{
	int fd=open(fname.c_str(),O_RDONLY);
	if(fd==-1) throw runtime_error("opening "+fname+": "+strerror(errno));
	struct stat st;
	if(fstat(fd,&st)) {
		close(fd);
		throw runtime_error("getting status "+fname+": "+strerror(errno));
	}
	m_MemorySize=st.st_size;
#ifdef USE_WIN
	m_MemoryBase=new char[m_MemorySize+1];
	read(fd,m_MemoryBase,m_MemorySize);
#else
	m_MemoryBase=(char*)mmap(0,m_MemorySize,PROT_READ|PROT_WRITE,
                             MAP_PRIVATE|MAP_NORESERVE,fd,0);
#endif
	close(fd);
    
	if(m_MemoryBase==0 || m_MemoryBase==caddr_t(-1))
		throw runtime_error("memory mapping "+fname+": "+strerror(errno));
	if(strchr("\r\n",m_MemoryBase[m_MemorySize-1])==0) {
		DetachFile();
		throw runtime_error("format error "+fname+": need newline before EOF");
	}
}

void CStsFileHash::DetachFile()
{
#ifdef USE_WIN
	delete[] m_MemoryBase;
#else
	munmap(const_cast<char*>(m_MemoryBase),m_MemorySize);
#endif
	m_MemoryBase=0;
	m_MemorySize=0;
}

void CStsFileHash::ReadStsFile(const string& fname, IStsFileCallback * cbk)
{
    // TODO::   dsjkghdgifdjgd not complete reset!!!
	Clear();

	AttachFile(fname);

#ifndef USE_WIN
	madvise(const_cast<char*>(m_MemoryBase),m_MemorySize,
            MADV_SEQUENTIAL|MADV_WILLNEED);
#endif

	if(cbk==0 || cbk->Start()) {
		const char * pos=m_MemoryBase;
		const char * end=m_MemoryBase+m_MemorySize;

		for(int to_parse=m_MemorySize; to_parse>1; to_parse=end-pos) {
			const char * nl=strpbrk(pos,"\n\r");

			if(!ParseLine(cbk,pos,nl-pos)) break;

			if(nl[0]!=nl[1] && strchr("\r\n",nl[1])) pos=nl+2;
			else pos=nl+1;
		}
		if(cbk) cbk->Done();
	}

#ifndef USE_WIN
	madvise(const_cast<char*>(m_MemoryBase),m_MemorySize,
            MADV_RANDOM|MADV_WILLNEED);
#endif
}

// lowercase characters are masked. 
inline bool IsMasked(char c) { return islower(c); }

static int Clip5PrimeLowercase(CStrRef& fld) 
{
    int len=fld.length();
    int oldlen=len;
    const char * c=fld.data();
    while( len>0 && IsMasked(*c) ) --len, ++c;
    fld=CStrRef(c,len);
    return oldlen-len;
}

bool CStsFileHash::ParseLine(IStsFileCallback * cbk, 
							 const char * pos, unsigned len) 
{
	if(cbk && !cbk->NextLine(pos,len)) return false;
	if(*pos!='#') {
		CStrRef fld[5];

		int cnt=CMmFileSts::Parse(pos,fld,4);
		if(cnt>3) {
            int oh1=0, oh2=0;
            
            if(AllowOverhang()) {
                oh1=Clip5PrimeLowercase(fld[1]);
                oh2=Clip5PrimeLowercase(fld[2]);
            }

			if(fld[1].length()<GetWordSize() ||
			   fld[2].length()<GetWordSize()) {
				return (cbk&&cbk->Error(IStsFileCallback::eErrShortPrimer));
			}
			const char * rev1=FlipSequence(fld[1].data(),fld[1].length());
			const char * rev2=FlipSequence(fld[2].data(),fld[2].length());
					
			int lo,hi;
			CMmFileSts::ParseRange(fld[3],lo,hi,
								   m_DefaultSizeLo,m_DefaultSizeHi);
			//cerr << "\e[31m" << __PRETTY_FUNCTION__ << "\e[032m: def-lo = " << m_DefaultSizeLo << ", def-hi = " << m_DefaultSizeHi << ", lo = " << lo << ", hi = " << hi << "\e[0m\n";

            const char * fwd1=(UnmaskPrimers()?
                               UCaseSequence(fld[1].data(),fld[1].length()):
                               fld[1].data());
            const char * fwd2=(UnmaskPrimers()?
                               UCaseSequence(fld[2].data(),fld[2].length()):
                               fld[2].data());
            int flags=(UnmaskPrimers()?CMmFileSts::fAllocLeft:0)|
                CMmFileSts::fAllocRight;

			CMmFileSts * fsts = new CMmFileSts(
				pos,
                CStrRef(fwd1,fld[1].length()),
                CStrRef(rev2,fld[2].length()),
                lo,hi,flags,oh1,oh2);

			CMmFileSts * rsts = new CMmFileSts(
				pos,
                CStrRef(fwd2,fld[2].length()),
                CStrRef(rev1,fld[1].length()),
                lo,hi,flags|CMmFileSts::fReverse,oh2,oh1);
						
			bool fok=AddStsEntry(fsts);
			bool rok=AddStsEntry(rsts);

			if(!sm_OneTimeRun) {
				if(!fok) delete fsts;
				if(!rok) delete rsts;
			}
			
			if(!fok && !rok)
				return (cbk&&cbk->Error(IStsFileCallback::eErrAmbiquosPrimer));
		} else if(cbk && !cbk->Error(IStsFileCallback::eErrBadLine))
			return false;
	}
	return true;
}

/*
 * $Log: stsmatch_m.cpp,v $
 * Revision 1.16  2008/06/18 14:45:33  rotmistr
 * Fixed problem with -d x-X parameter being reset if -w or some others are used after it.
 *
 * Revision 1.15  2007/07/11 20:49:30  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.14  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.13  2004/10/26 17:16:35  rotmistr
 * Added 5'-end masking for primers
 *
 * Revision 1.12  2004/09/03 21:28:50  rotmistr
 * Fixes to compile with Borland C++ 5.5
 *
 * Revision 1.11  2004/06/03 23:37:23  rotmistr
 * New aligner added.
 *
 * Revision 1.10  2004/05/27 20:35:49  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.9  2004/04/01 17:23:20  rotmistr
 * *** empty log message ***
 *
 * Revision 1.8  2004/04/01 16:37:42  rotmistr
 * Cleaned after adding windows capabilities
 *
 * Revision 1.7  2004/04/01 05:57:53  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.6  2004/03/30 21:06:53  rotmistr
 * Fixes for setting default STS size range.
 *
 * Revision 1.5  2004/03/30 19:08:03  rotmistr
 * default STS size is tunnable now
 *
 * Revision 1.4  2004/03/23 22:35:26  rotmistr
 * Fixed processing of -mid flag in cmdline
 * Fixed destructor for fasta reader
 * Removed cgi
 *
 * Revision 1.3  2004/03/07 06:36:00  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.2  2004/02/04 21:23:22  rotmistr
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
