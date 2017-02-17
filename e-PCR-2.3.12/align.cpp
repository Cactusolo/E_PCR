/* $Id: align.cpp,v 1.5 2004/10/26 17:16:32 rotmistr Exp $
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

//#include <epcr/hashset.hpp>
//#include <epcr/stsmatch_i.hpp>

#include <epcr/align.hpp>
#include <stdio.h>
#include <string.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

static char s_Compl[]=
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NTBGHNNCDNNMNKNNNNYNAABWXRNNNNNN"
//ABCDEFGHIJKLMNOPQRSTUVWXYZNNNNN
"NTBGHNNCDNNMNKNNNNYNAABWXRNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
;
static char s_Lc2Uc[]=
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NABCDNNGHNNKNMNNNNRNTTVWXYNNNNNN"
//ABCDEFGHIJKLMNOPQRSTUVWXYZNNNNN
"NABCDNNGHNNKNMNNNNRNTTVWXYNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
;

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

inline char UCase(char c) { return s_Lc2Uc[int((unsigned char)c)]; }
char * UCaseSequence(const char * seq, unsigned len) 
{
    if(len == ~0U) len=strlen(seq);
    char * ret=new char[len+1];
    ret[len]=0;
    for(char * x=ret; len!=0; --len) *x++=s_Lc2Uc[int(*seq++)];
    return ret;
}

char * FlipSequence(const char * seq, unsigned len) 
{
    if(len == ~0U) len=strlen(seq);
    char * ret=new char[len+1];
    ret[len]=0;
    for(char * x=ret+len; len!=0; --len) *--x=s_Compl[int(*seq++)];
    return ret;
}

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

bool CAlignExact::Forward(const char * seq_ptr,
                          const char * seq_end,
                          const char * primer,
                          int length) 
{ 
    if(strncmp(seq_ptr,primer,length)==0) {
        m_Identities=length;
        return true;
    } else {
        m_Identities=0;
        return false;
    }
}

bool CAlignExact::Reverse(const char * seq_start,
                          const char * seq_ptr,
                          const char * primer,
                          int length) 
{
    seq_ptr-=length;
    if(seq_ptr>=seq_start &&
       strncmp(seq_ptr,primer,length)==0) {
        m_Identities=length;
        return true;
    } else {
        m_Identities=0;
        return false;
    }
}

bool CAlignNoGaps::Forward(const char * seq_ptr,
                           const char * seq_end,
                           const char * primer,
                           int length) 
{ 
    return x_Compare(seq_ptr,primer,length);
}

bool CAlignNoGaps::Reverse(const char * seq_start,
                           const char * seq_ptr,
                           const char * primer,
                           int length) 
{
    seq_ptr-=length;
    if(seq_ptr>=seq_start) {
        return x_Compare(seq_ptr,primer,length);
    }
    m_Identities=0;
    return false;
}


bool CAlignNoGaps::x_Compare(const char * a, 
                             const char * b, int l)
{
    m_Identities=0;
    if(l<=0) return false;

    int mism=m_Mism=m_MaxMismatch;
 	for(;l; --l) {
 		if(*a++ != *b++) {
            if(a[-1]==0 || --mism<0) return false;
        } 
    }

    m_Mism-=mism;
    m_Identities=l-m_Mism;

    return l==0;
}

bool CAlignFast::Forward(const char * a, const char * A,
                         const char * b, int l) 
{
	const char * B=b+l;

    int gaps=m_Gaps=m_MaxGaps;
    int mism=m_Mism=m_MaxMismatch;
    m_Identities=0;
    
 	for(;b<B && *a && *b && a<A; ++a, ++b) {
 		if(*a == *b) { m_Identities++; continue; }
 		if(a[1] == b[1]) {
 			if(--mism<0) return false;
 			++a,++b;
 		} else if(gaps && a[1] == b[0]) {
 			// mismatches are preferred to gaps
 			if(a[1] == b[1]) {
 				if(mism >=2 ) { mism-=2; a+=2; b+=2; continue; }
 			}
 			if(--gaps<0) return false;
 			++a;
 		} else if(gaps && a[0] == b[1]) {
 			if(--gaps<0) return false;
 			++b;
 		} else {
            if(--mism<0) return false;
        }
 	}
 	for(;*b && b<B; ++b) if(--gaps<0) return false;

    m_Gaps-=gaps;
    m_Mism-=mism;
	
 	return true;
}

bool CAlignFast::Reverse(const char * A, const char * a, 
                         const char * b, int l) 
{
	const char * B=b;
	b+=l-1;
	--a;

    int gaps=m_Gaps=m_MaxGaps;
    int mism=m_Mism=m_MaxMismatch;
    m_Identities=0;
    
 	for(;a>=A && b>=B; --a, --b) {
 		if(*a == *b) { m_Identities++; continue; }

		if(a>A && b>B && a[-1] == b[-1]) {
			if(--mism<0)  return false;
			--a,--b;
		} else if(gaps && a>A && a[-1] == b[0]) {
			// mismatches are preferred to gaps
			if(b>B && a[0] == b[-1]) {
				if(mism >=2 ) { mism-=2; a-=2; b-=2; continue; }
			}
			if(--gaps<0) return false;
			--a;
		} else if(gaps && b>B && a[0] == b[-1]) {
 			if(--gaps<0) return false;
 			--b;
 		} else {
            if(--mism<0) return false;
        }
 	}
 	for(;b>=B; --b) if(--gaps<0) return false;

    m_Gaps-=gaps;
    m_Mism-=mism;
	
 	return true;
}

bool CAlignLCS::Forward(const char * seq_ptr,
                        const char * seq_end,
                        const char * primer,
                        int length) 
{
    m_Matrix.Build<const char*>(seq_ptr,seq_end,primer,length);
    m_Matrix.Stat<const char*>(seq_ptr,seq_end,primer,length);
    return 
        m_Matrix.GetMismatches()<=m_MaxMismatch &&
        m_Matrix.GetGaps()<=m_MaxGaps;
}
        
bool CAlignLCS::Reverse(const char * seq_start,
                        const char * seq_ptr,
                        const char * primer,
                        int length)
{
    m_Matrix.Build<CReverseConstSeqIterator<const char> >(
        seq_ptr-1,seq_start-1,primer+length-1,length);
    m_Matrix.Stat<CReverseConstSeqIterator<const char> >(
        seq_ptr-1,seq_start-1,primer+length-1,length);
    return 
        m_Matrix.GetMismatches()<=m_MaxMismatch &&
        m_Matrix.GetGaps()<=m_MaxGaps;
}

CAlignLCS::CAlignLCS(int mm, int gg):
    m_MaxMismatch(mm),m_MaxGaps(gg),m_Matrix(256,gg)
{}

CAlignCompromise::CAlignCompromise(int mm, int gg):
    CAlignLCS(mm,gg), CAlignNoGaps(mm) 
{}

bool CAlignCompromise::Forward(const char * seq_ptr,
                               const char * seq_end,
                               const char * primer,
                               int length)
{
    if(CAlignNoGaps::Forward(seq_ptr,seq_end,primer,length)) {
        m_Gaps=CAlignNoGaps::GetGaps();
        m_Mismatches=CAlignNoGaps::GetMismatches();
        m_Identities=CAlignNoGaps::GetIdentities();
        return true;
    } 
    if(CAlignLCS::Forward(seq_ptr,seq_end,primer,length)) {
        m_Gaps=CAlignLCS::GetGaps();
        m_Mismatches=CAlignLCS::GetMismatches();
        m_Identities=CAlignLCS::GetIdentities();
        return true;
    }
    return false;
}

bool CAlignCompromise::Reverse(const char * seq_start,
                               const char * seq_ptr,
                               const char * primer,
                               int length)
{
    if(CAlignNoGaps::Reverse(seq_start,seq_ptr,primer,length)) {
        m_Gaps=CAlignNoGaps::GetGaps();
        m_Mismatches=CAlignNoGaps::GetMismatches();
        m_Identities=CAlignNoGaps::GetIdentities();
        return true;
    }
    if(CAlignLCS::Reverse(seq_start,seq_ptr,primer,length)) {
        m_Gaps=CAlignLCS::GetGaps();
        m_Mismatches=CAlignLCS::GetMismatches();
        m_Identities=CAlignLCS::GetIdentities();
        return true;
    }
    m_Identities=0;
    return false;
}

    
/*
 * $Log: align.cpp,v $
 * Revision 1.5  2004/10/26 17:16:32  rotmistr
 * Added 5'-end masking for primers
 *
 * Revision 1.4  2004/06/08 20:32:50  rotmistr
 * Fixup for gap+insert special case
 *
 * Revision 1.3  2004/06/08 16:14:55  rotmistr
 * *** empty log message ***
 *
 * Revision 1.2  2004/06/07 16:24:56  rotmistr
 * Bug fixes to previos version.
 *
 * Revision 1.1  2004/06/03 23:37:19  rotmistr
 * New aligner added.
 *
 *
 */
