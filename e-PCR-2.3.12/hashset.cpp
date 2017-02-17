/* $Id: hashset.cpp,v 1.2 2008/06/16 16:02:40 rotmistr Exp $
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

#include <epcr/hashset.hpp>
#include <string.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

unsigned char CHashSet::sm_HashId[]=
"********************************"
"********************************"
"*\0*\1***\2************\3***********"
//=Ab=Cdef=Ghijklmnopqrs=Tuvwxyz*****//
"*\0*\1***\2************\3***********"
"********************************"
"********************************"
"********************************"
"********************************";
inline bool s_IsAmbiq(unsigned char a) { return a>3; }

typedef CHashSet::hash_type hash_type;
typedef CHashSet::size_type size_type;
typedef CHashSet::bits_type bits_type;

/*
static char * s_BitsRepr(unsigned u, int w=32) 
{
	static char ret[65];
	
	if(w>32) w=31;
	if(w<1) w=1;
	
	for(unsigned m=(1<<(w-1)), i=0; m; ++i,(m>>=1)) {
		ret[i]=(m&u)?'1':'0';
	}
	ret[w]=0;
	return ret;
}
*/

/*
************************************************************************
m_mask[*] 0000111111111
m_word[0] 0000110110110
m_word[1] 0000101101101
m_word[2] 0000011011011
m_ambiq   0000001000000
str(rev)  acgtACNGTCATGcgatagat
                       ^ ptr, offset
************************************************************************
*/

bool CHashSet::Begin(const char * ptr)
{
	m_Ptr=ptr;
	m_Offset=0;
	m_AmbMask= ~0U; // all are ambiquos in the very beginning
//	for(size_type t=0; t<m_Period; ++t) m_Hash[t]=0; 
	memset(m_Hash,0,m_Period*sizeof(*m_Hash));

	if(ptr==0) return false;
	
	if(m_Period==1) {
		for(;m_Offset < m_WdSize && *m_Ptr; ++m_Ptr, ++m_Offset) {
			m_AmbMask<<=1;
					
			unsigned char nt=sm_HashId[(int)*m_Ptr];
			if(nt>3) { m_AmbMask|=1; nt=0; }
			
					// Circular 'rotation' with update of values
			m_Hash[0]=m_Mask[0]&((m_Hash[0]<<2)|nt);
		}
	}
	else {
		for(;m_Offset < m_WdSize && *m_Ptr; ++m_Ptr, ++m_Offset) {
			m_AmbMask<<=1;
			
			unsigned char nt=sm_HashId[(int)*m_Ptr];
			if(nt>3) { m_AmbMask|=1; nt=0; }
			hash_type x=m_Hash[m_Period-1];
			for(size_type t=m_Period-1; t>0; --t) {
				m_Hash[t]=m_Mask[t]&((m_Hash[t-1]<<2)|nt);
			}
			m_Hash[0]=x&m_Mask[0];
		}
	}

	return !End();
}

bool CHashSet::Good() const 
{
	for(size_type t=0; t< m_Period; ++t) if(!Good(t)) return false;
	return true;
}

CHashSet::~CHashSet() throw ()
{
	_intl_free();
}

void CHashSet::_intl_free() 
{
	
	delete[] m_Word;
	delete[] m_Mask;
	delete[] m_Hash;
	delete[] m_TbSize;
}

CHashSet::CHashSet(TSize wdsize, TSize period)
	:m_Ptr(0),m_Offset(0),m_AmbMask(0)
{
	if(period==0 || period>wdsize) period=1;
#if 0
	//  Assure that hash ranged are equal
	if(wdsize%period) throw logic_error("Uneven hash masks");
#endif
	m_TbSize = new TSize[period];
	m_Word   = new TBitMask[period];
	m_Mask   = new THashValue[period];
	m_Hash   = new THashValue[period];
	m_WdSize = wdsize;
	m_Period = period;
	if(period<=1) {
		m_TbSize[0] = 1<<(2*wdsize);
		m_Mask[0]=m_TbSize[0]-1;
		m_Word[0]=(1<<wdsize)-1;
		m_Hash[0]=0;
	}
	else {
		for(size_type i=0; i<period; ++i) m_Word[i]=m_Mask[i]=0;

		for(size_type i=0, bit=1; i<m_WdSize; ++i, (bit<<=1)) {
			m_Word[i%period]|=bit;
			m_Mask[i%period]++; // count bits off
		}
		bits_type wdmask=((1<<m_WdSize)-1);
		for(size_type i=0; i<period; ++i) {
			m_TbSize[i]=1<<((wdsize-m_Mask[i])*2);
			m_Mask[i]=m_TbSize[i]-1;
			m_Word[i]= (~m_Word[i]) & wdmask;//&m_Mask[i];
			m_Hash[i]=0;
		}
	}
}


CHashSet::CHashSet(const CHashSet& s)
{
	_intl_copy(s);
}

void CHashSet::_intl_copy(const CHashSet& s)
{
	m_Ptr=s.m_Ptr;
	m_Offset=s.m_Offset;
	m_AmbMask=s.m_AmbMask;
	m_WdSize=s.m_WdSize;
	m_Period=s.m_Period;

	m_TbSize = new TSize[m_Period];
	m_Word   = new TBitMask[m_Period];
	m_Mask   = new THashValue[m_Period];
	m_Hash   = new THashValue[m_Period];

	for(size_type t=0; t<m_Period; ++t) {
		m_TbSize[t]=s.m_TbSize[t];
		m_Word  [t]=s.m_Word  [t];
		m_Mask  [t]=s.m_Mask  [t];
		m_Hash  [t]=s.m_Hash  [t];
	}
}

CHashSet& CHashSet::operator = (const CHashSet& s)
{
	if(&s!=this) {
		_intl_free();
		_intl_copy(s);
	}
	return *this;
}

////////////////////////////////////////////////////////////////////////
#if defined TEST_MAIN
#include <string>
#include <iostream>
#include <iomanip>

#include <unistd.h>

using namespace std;

int main(int argc, char ** argv) 
{
	try {
		int optchar;
		int period=0;
		int wdsize=7;
		
		while((optchar=getopt(argc,argv,"hVp:w:"))!=-1) {
			switch(optchar){
			case 'V':
			case 'h': 
				cout<< "usage: [-hV] [-p period] [-w wordsize] [string ...]\n";
				return 0;
			case 'p': period=atoi(optarg); break;
			case 'w': wdsize=atoi(optarg); break;
			}
		}
		
		CHashSet hset(wdsize,period);
		cout 
			<< "WordSize = " << hset.GetWordSize() << endl
			<< "TabCount = " << hset.GetWordCount() << endl
			<< "Id\t" 
			<< setw(16) << "Word" << "\t" 
			<< setw(16) << "Mask" << "\t" 
			<< "Table size\n";
		
		for(size_type i=0; i<hset.GetWordCount(); ++i) {
			cout 
				<< i << "\t"
				<< s_BitsRepr(hset.GetWord(i),16) << "\t"
				<< s_BitsRepr(hset.GetMask(i),16) << "\t"
				<< hset.GetTableSize(i) << "\n";
		}
		
		for(int i=optind; i<argc; ++i) {
			const char * ptr=argv[i];
			if(strlen(ptr)<wdsize) continue;
			
			hset.Begin(ptr);

 			for(size_type p=0; p<wdsize; ++p)
 				cout << "\n" << setw(5) << p << " " << ptr[p];
			
 			for(size_type t=0; t<hset.GetWordCount(); ++t) {
 				cout << " " << (hset.Good(t)?'+':'-') 
 					 << s_BitsRepr(hset.GetValue(t),16);
 			}
 			cout << " [" << s_BitsRepr(hset.GetAmbiquityMask(),16) << "]\n";
			
			for(; !hset.End();) {
				cout << setw(5) << hset.GetPosition() 
					 << " " << hset.GetPtr()[0];
				if(!hset.Next()) break;
				
				for(size_type t=0; t<hset.GetWordCount(); ++t) {
					cout << " " << (hset.Good(t)?'+':'-') 
						 << s_BitsRepr(hset.GetValue(t),16);
				}
				cout << " [" 
                     << s_BitsRepr(hset.GetAmbiquityMask(),16) << "]\n";
			}
		}
		return 0;
	} 
	catch(exception& e) {
		cerr << "! Error: " << e.what() << endl;
	}
	return 100;
}

#endif			



/*
 * $Log: hashset.cpp,v $
 * Revision 1.2  2008/06/16 16:02:40  rotmistr
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
