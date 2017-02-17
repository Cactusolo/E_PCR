/* $Id: stsmatch_i.hpp,v 1.10 2005/06/14 16:46:44 rotmistr Exp $
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

#ifndef EPCR_STSMATCH_I__HPP
#define EPCR_STSMATCH_I__HPP

#include <epcr/build_cfg.h>
#include <epcr/hashset.hpp>
#include <epcr/sts_i.hpp>
#include <epcr/align.hpp>

#include <stdexcept>
#include <vector>
#include <list>
#include <map>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class CStsHash;
class CPcrMachine;
class IPcrMachineCallback;
class IPcrProgressCallback;

// This class implements STS hash table, but lacks table loading methods
class CStsHash
{
public:
    typedef vector<ISts*> TStsList;
    typedef TStsList *  THashTable;
    typedef THashTable* TData;

    CStsHash():m_Table(0) {}

    void SetHash(const CHashSet& hs);
    const CHashSet& GetHash() const { return m_Hash; }
    
    virtual void Reset(); // deletes everything, if sm_OneTimeRun is not set
    virtual void Clear(); // clears hash lists and deletes data

    virtual ~CStsHash() throw () { if(sm_OneTimeRun) Reset(); }

    unsigned GetWordSize()  const { return m_Hash.GetWordSize(); }
    unsigned GetWordCount() const { return m_Hash.GetWordCount(); }

    static void SetOneTimeRun(bool v) { sm_OneTimeRun=v; }
protected:
    virtual bool AddStsEntry(ISts * sts);

protected:
    static bool sm_OneTimeRun;
    friend class CPcrMachine;

    TData    m_Table;
    CHashSet m_Hash;
    TStsList m_All;
};

// This class implements STS match algorithm
class CPcrMachine
{
public:
    typedef CStsHash::TStsList   TStsList;
    typedef CStsHash::THashTable THashTable;
    typedef CStsHash::TData      TData;

    virtual void SetMargin(unsigned margin) { m_Margin=margin; }
    unsigned GetMargin() const { return m_Margin; }

    virtual void SetAligner(IAlign * alignl, IAlign * alignr) 
    { m_AlignL=alignl; m_AlignR=alignr; }
        
    virtual void SetProgressCallback(IPcrProgressCallback * callback) 
        { m_Progress=callback; }
    virtual void SetCallback(IPcrMachineCallback * callback) 
        { m_Callback=callback; }
    virtual void SetStsHash(const CStsHash* stshash) { 
        m_StsHash=stshash; 
        if(m_StsHash) { m_HashL=m_HashR=m_HashS=stshash->m_Hash; } 
    }

    virtual void ProcessSequence(const char * label,
                                 const char * seqdata,
                                 unsigned     seqlen);
    CPcrMachine():m_AlignL(0),m_AlignR(0),m_Callback(0),m_StsHash(0),m_Progress(0) {}
    virtual ~CPcrMachine() throw () {}

protected:
    virtual void Match(unsigned word, const TStsList& lst,
                       const char * seq_data, int seq_len);
protected:
    unsigned m_Margin;
    
    IAlign * m_AlignL, * m_AlignR;
    CHashSet m_HashL, m_HashR, m_HashS;

    IPcrMachineCallback  * m_Callback;
    const CStsHash       * m_StsHash;
    IPcrProgressCallback * m_Progress;
};


// Callback that is used in CPcrMachine::ProcessSequence
class IPcrMachineCallback
{
public:
    struct SScore
    {
//         int n_mism;
//         int n_gaps;
        int actlen;
		char mism_l, mism_r, gaps_l, gaps_r;
        SScore(int a=0, 
			   char m_l=0, char m_r=0, 
			   char g_l=0, char g_r=0):
            actlen(a),
            mism_l(m_l),mism_r(m_r),
            gaps_l(g_l),gaps_r(g_r){}
    };
    
    virtual void CbkStart() {};
    virtual void CbkEnd() {};
    virtual void CbkSequenceEnd() {};
    virtual void CbkSequence(const char * label) =0;
    virtual void CbkSequenceData(const char * data, unsigned size) {};
    virtual void CbkMatch(const ISts * sts, unsigned pos1, unsigned pos2,
                          const SScore *score) =0;
    virtual void CbkWarning(const char * message) =0;

	virtual ~IPcrMachineCallback() {}
	
};

// Separate to optimize CPcrMachine::ProcessSequence scanning loop 
// (callback only if it is not null)
class IPcrProgressCallback
{
public:
    virtual ~IPcrProgressCallback() {}
    virtual void PgsSequenceStart(const char * label, const char * data, 
                                  unsigned length, unsigned wsize) = 0;
    virtual void PgsSequenceEnd() = 0;
    virtual void PgsSequenceAt(unsigned pos) = 0;
};

// Accumulates hits, removes redundant hits as well as suboptimal hits
class CPcrMachinePostprocess:public IPcrMachineCallback
{
public:
	~CPcrMachinePostprocess() {}
	
    CPcrMachinePostprocess(IPcrMachineCallback * out):m_Callback(out) {}

    virtual void CbkStart() { m_Callback->CbkStart(); }
    virtual void CbkEnd() { m_Callback->CbkEnd(); }
    virtual void CbkSequenceEnd() { Flush(); m_Callback->CbkSequenceEnd(); }
    virtual void CbkSequence(const char * label) 
        { Flush(); m_Callback->CbkSequence(label); }
    virtual void CbkMatch(const ISts * sts, unsigned pos1, unsigned pos2,
                          const SScore *score);
    virtual void CbkWarning(const char * message) {
        m_Callback->CbkWarning(message); 
    }
    virtual void CbkSequenceData(const char * data, unsigned size) {
        m_Callback->CbkSequenceData(data,size);
    };
    void Flush();
    
protected:
    IPcrMachineCallback * m_Callback;
    
    struct SOutput 
    {
        int pos1, pos2;
        char mism_l, mism_r, gaps_l, gaps_r;
        int length() const { return pos2-pos1; }

		int mism() const { return mism_l+mism_r; }
		int gaps() const { return gaps_l+gaps_r; }
	
		bool operator == (const SOutput& o) const {
			return pos2==o.pos2 && pos1==o.pos1 &&
				mism_l==o.mism_l && gaps_l==o.gaps_l &&
				mism_r==o.mism_r && gaps_r==o.gaps_r;
		}
        SOutput(int p1, int p2, char m_l, char m_r, char g_l, char g_r):
            pos1(p1),pos2(p2),mism_l(m_l),mism_r(m_r),gaps_l(g_l),gaps_r(g_r){}
        SOutput(){}
    };
    static bool OrderByPos2Pos1(const SOutput&, const SOutput&);
    static int  Compare(const SOutput&, const SOutput&, int, int);
    static bool Overlap(int, int, const SOutput&, const SOutput&);

    typedef vector<SOutput> TStsHits;
    typedef map<const ISts*,TStsHits> TAllHits;

    typedef TStsHits::iterator TStsHits_I;
    typedef TAllHits::iterator TAllHits_I;

    typedef TStsHits::const_iterator TStsHits_CI;
    typedef TAllHits::const_iterator TAllHits_CI;
    
    TAllHits m_OutQueues;
};


END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: stsmatch_i.hpp,v $
 * Revision 1.10  2005/06/14 16:46:44  rotmistr
 * Changed report format for floppy tails
 *
 * Revision 1.9  2004/06/07 16:24:57  rotmistr
 * Bug fixes to previos version.
 *
 * Revision 1.8  2004/06/03 23:37:22  rotmistr
 * New aligner added.
 *
 * Revision 1.7  2004/03/26 17:02:13  rotmistr
 * Compat-options are now allowed everywhere, and multiple fasta files can be used.
 *
 * Revision 1.6  2004/03/25 19:36:52  rotmistr
 * API: separate left and right primers mism/gaps in forward API
 *
 * Revision 1.5  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.4  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
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
