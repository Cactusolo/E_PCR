/* $Id: stsmatch_i.cpp,v 1.13 2005/06/14 16:46:44 rotmistr Exp $
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
#include <epcr/stsmatch_i.hpp>

#include <stdexcept>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

#define DBG fprintf(stderr,__FILE__":%d\n",__LINE__)
#define INT(i) fprintf(stderr,__FILE__":%d - %s = %d\n",__LINE__,#i,i)
#define STR(s,i) fprintf(stderr,__FILE__":%d - %s = %.*s\n",__LINE__,#s,i,s)

bool CStsHash::sm_OneTimeRun(false);

void CStsHash::Reset()
{
    if(m_Table) {
        // NB: It DOES NOT delete STSes themselve!!!
        for(unsigned i=0; i<m_Hash.GetWordCount(); ++i)
            delete[] m_Table[i];
        
        delete[] m_Table;
    }
    m_Table=0;
    if(!sm_OneTimeRun) {
        for(TStsList::iterator i=m_All.begin(); i!=m_All.end(); ++i) {
            delete *i;
        }
        m_All.clear();
    }
}

void CStsHash::Clear()
{
    if(m_Table) {
        for(unsigned i=0; i<m_Hash.GetWordCount(); ++i)
            m_Table[i]->clear();
    }
    
    for(TStsList::iterator i=m_All.begin(); i!=m_All.end(); ++i) {
        delete *i;
    }
    m_All.clear();
}


void CStsHash::SetHash(const CHashSet& hs)
{
    Reset();
    
    m_Hash=hs;

	m_Table=new THashTable[m_Hash.GetWordCount()];
	for(size_t i=0; i<m_Hash.GetWordCount(); ++i) {
		m_Table[i]=new TStsList[m_Hash.GetTableSize(i)];
	}
}

bool CStsHash::AddStsEntry(ISts * sts) 
{
    CStrRef primer=sts->GetPrimer(ISts::eLeft);
    if(primer.length()<GetWordSize()) return false;

    bool ok=false;
    
    m_Hash.Begin(primer.data()+(primer.length()-GetWordSize()));
    for(unsigned it=0; it<GetWordCount(); ++it) {
        if(m_Hash.Good(it)) {
            m_Table[it][m_Hash[it]].push_back(sts);
            ok=true;
        }
    }    
    if(ok && !sm_OneTimeRun) m_All.push_back(sts);
    return ok;
}

////////////////////////////////////////////////////////////////////////
void CPcrMachine::ProcessSequence(const char * label,
                                  const char * seq_data,
                                  unsigned seq_len) 
{
    if(!m_Callback || !m_StsHash) 
        throw logic_error("CPcrMachine needs SetHash and SetCallback calls");
    if(m_StsHash->GetWordCount()==0)
        throw logic_error("Word count should be greater then zero");
    if(m_StsHash->m_Table==0)
        throw logic_error("Sts table should be initialized");
    if(seq_len==0) seq_len=strlen(seq_data);
    
    m_Callback->CbkSequence(label);
    m_Callback->CbkSequenceData(seq_data,seq_len);

    CHashSet& hash=m_HashL;

    if(m_Progress) 
        m_Progress->PgsSequenceStart(label,seq_data,seq_len,
                                     hash.GetWordSize());

    if(seq_len<hash.GetWordSize()) {
        m_Callback->CbkWarning("too short sequence");
        m_Callback->CbkSequenceEnd();
        if(m_Progress) m_Progress->PgsSequenceEnd();
        return;
    }
    
    if(hash.GetWordCount()==1) {
        CStsHash::THashTable tab=m_StsHash->m_Table[0];
        for(hash.Begin(seq_data); !hash.End(); hash.Next()) {
            if(m_Progress) m_Progress->PgsSequenceAt(hash.GetPosition());
            if(!hash.Good(0)) continue;
            if(tab[hash[0]].size()) Match(0, tab[hash[0]], seq_data, seq_len);
        }
    } else {
        CStsHash::TData tab=m_StsHash->m_Table;
    
        for(hash.Begin(seq_data); !hash.End(); hash.Next()) {
            if(m_Progress) m_Progress->PgsSequenceAt(hash.GetPosition());
            for(unsigned word=0; word<hash.GetWordCount(); ++word) {
                if(!hash.Good(word)) continue;
                if(tab[word][hash[word]].size())
                    Match(word, tab[word][hash[word]], seq_data, seq_len);
            }

        }
    }
    if(m_Progress) m_Progress->PgsSequenceEnd();
    m_Callback->CbkSequenceEnd();
}

void CPcrMachine::Match(unsigned word, const TStsList& lst,
                        const char * seq_data, int seq_len) 
{
    int seq_pos = m_HashL.GetPosition();
    const char * seq_cur = seq_data+seq_pos;//+pcr_p1_len;
    int cur_len = seq_len - seq_pos;

    for(TStsList::const_iterator ists=lst.begin(); ists!=lst.end(); ++ists) {

        const ISts * sts=*ists;
    
    
        int pcr_p1_len  = sts->GetPrimerLength(ISts::eLeft);
        
        if(m_HashL.GetPosition() < (unsigned)pcr_p1_len) continue;
        const char * lprimer = sts->GetPrimerData(ISts::eLeft);
                
        bool match=m_AlignL->Reverse(seq_data,seq_cur,lprimer,pcr_p1_len);
        
        if(match) {
            int pcr_p2_len = sts->GetPrimerLength(ISts::eRight);

            int pcr_len = pcr_p1_len + pcr_p2_len;
            
            const char * p = seq_cur + sts->GetSizeLo() - pcr_len - m_Margin;
            const char * P = seq_cur + sts->GetSizeHi() - pcr_len + m_Margin;
            
            if( p < seq_cur - pcr_p1_len ) p = seq_cur - pcr_p1_len;
            if( p < seq_data ) p = seq_data;
            if( P > seq_cur + cur_len - pcr_p2_len )
                P = seq_cur + cur_len - pcr_p2_len;

            if( p > P ) continue;

            const char * rprimer = sts->GetPrimerData(ISts::eRight);
                
#ifndef USE_HASH_FOR_RIGHT_PRIMER
#define USE_HASH_FOR_RIGHT_PRIMER 1
#endif

            int ovhg1 = sts->GetOverhangChars(ISts::eLeft);
            int ovhg2 = sts->GetOverhangChars(ISts::eRight);
            int ovhgall = ovhg1+ovhg2;

#if !USE_HASH_FOR_RIGHT_PRIMER
            for(; p<=P; --P) {

                if(m_AlignR->Forward(P,seq_data+seq_len,rprimer,pcr_p2_len)) {
                    
                    IPcrMachineCallback::SScore score(
                        P-p+pcr_len+ovhgall,
                        m_AlignL->GetMismatches(),
                        m_AlignR->GetMismatches(),
                        m_AlignL->GetGaps(),
                        m_AlignR->GetGaps());
                    
                    m_Callback->CbkMatch(
                        sts,
                        seq_pos-pcr_p1_len-ovhg1,
                        P-seq_data+pcr_p2_len+ovhg2,
                        &score);
                }
            }
#else
            CHashSet& rhash=m_HashR;
            CHashSet& mhash=m_HashS;
            rhash.Begin(sts->GetPrimerData(ISts::eRight));

            for (mhash.Begin(p); p<=P; p++, mhash.Next()) {
                for(unsigned wd=0; wd<rhash.GetWordSize(); ++wd) {
                    if(!rhash.Good(wd) || !mhash.Good(wd)) continue;
			
                    if(mhash.GetValue(wd)!=rhash.GetValue(wd)) continue;
                
                    if(m_AlignR->Forward(p,seq_data+seq_len,
                                         rprimer,pcr_p2_len)) {
                        
                        IPcrMachineCallback::SScore score(
                            P-p+pcr_len+ovhgall,
                            m_AlignL->GetMismatches(),
                            m_AlignR->GetMismatches(),
                            m_AlignL->GetGaps(),
                            m_AlignR->GetGaps());
                    
                        m_Callback->CbkMatch(sts,
                                             seq_pos-pcr_p1_len-ovhg1, 
                                             p-seq_data+pcr_p2_len+ovhg2,
                                             &score);
                        
                        break; // no other reports for this positions 
                    }
                }
            }
#endif
        }
	
    }
}

/*
 * $Log: stsmatch_i.cpp,v $
 * Revision 1.13  2005/06/14 16:46:44  rotmistr
 * Changed report format for floppy tails
 *
 * Revision 1.12  2004/10/26 17:16:35  rotmistr
 * Added 5'-end masking for primers
 *
 * Revision 1.11  2004/06/03 23:37:22  rotmistr
 * New aligner added.
 *
 * Revision 1.10  2004/05/27 21:12:46  rotmistr
 * Some warnings fixed.
 *
 * Revision 1.9  2004/05/27 20:35:49  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.8  2004/03/31 05:04:00  rotmistr
 * Search range fix
 *
 * Revision 1.7  2004/03/25 19:36:52  rotmistr
 * API: separate left and right primers mism/gaps in forward API
 *
 * Revision 1.6  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.5  2004/02/18 05:43:22  rotmistr
 * Bug fix with search range
 *
 * Revision 1.4  2004/02/11 04:34:57  rotmistr
 * Optimised lookup speed and memory usage
 * Fixed bug with end of sequence in stsmatch
 * Changing CGI look
 *
 * Revision 1.3  2004/01/28 23:27:02  rotmistr
 * "Best of overlapping" hit selection postprocessor added.
 *
 * Revision 1.2  2004/01/08 23:22:41  rotmistr
 * Fixed init error in faread,
 * Adjusted output to standard,
 * Added output format style and output file to parameters.
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
