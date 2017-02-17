/* $Id: fahash.hpp,v 1.15 2007/07/11 20:49:29 rotmistr Exp $
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

#ifndef EPCR_HASH__HPP
#define EPCR_HASH__HPP

#include <epcr/fast_seqio.hpp>
#include <epcr/hashset.hpp>
#include <epcr/sts_i.hpp>
#include <epcr/align.hpp>

#include <list>
#include <map>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class CFaData;
class CFaIndexer;
class CFaLookup;
class IFaLookupCallback;

class CFaData
{
public:
    struct  SSeqDesc
    {
        Int2    fid; 
        off64_t data;
        SSeqDesc( Int2 f, off64_t o ) : fid( f ), data( o ) {}
        SSeqDesc() : fid( -1 ), data( 0 ) {}
    };
    typedef std::string string;
    struct SFileDesc 
    {
        Int2 flags;
        string path;
        SFileDesc() : flags( 0 ) {}
        SFileDesc( Int2 fl, const string& s ) : flags( fl ), path( s ) {}
    };
				
    typedef std::vector<SFileDesc> TPathLst;
    typedef std::vector<SSeqDesc> TSeqLst;
    typedef Int4 THashElement;
//	typedef util::CChunkArray<4092,1,THashElement> THashList;
    typedef std::vector<THashElement> THashList;
//  (-)seq-id, offset, offset, (-)seq-id, offset, ..., -1
    typedef std::vector<THashList> THashTable;
    typedef std::vector<THashTable> TData;

    typedef Uint8 TBigCount;
    typedef vector<TBigCount> TBranchStat;
    typedef vector<vector<TBigCount> > TStatVector;

    TPathLst m_Fapath; // is accumulated and appended to file
    TSeqLst  m_Seqlst; // is accumulated and appended to file
    CHashSet m_Hash;   // hash itself
    std::vector< std::pair<off64_t,off64_t> > m_Tabloc;
    const CHashSet& GetHash() const { return m_Hash; }   // hash itself
};

class IFaIndexerCallback
{
public:
    virtual ~IFaIndexerCallback() {}
    virtual void CbkSequence( const char * id ) = 0;
    virtual void CbkProgress( unsigned pos, unsigned length ) = 0;
    virtual void CbkFile( const char * name ) =0;
    virtual void CbkDumpProgress( unsigned pos, unsigned length ) = 0;
    virtual void CbkDumpStart() {};
    virtual void CbkDumpEnd() {};
    virtual void CbkResetStart() {};
    virtual void CbkResetEnd() {};
};

class AFaIndexerBase : public CFaData
{
public:
    virtual void AttachFile( const string& path ) = 0;
				
    virtual void SetHash( const CHashSet& hs ) = 0;
    virtual void AddFile( const string& path );
    virtual void Finish();

    virtual Uint4 GetFlags() const { return m_Flags; }
    virtual void SetFlags( unsigned f ) { m_Flags = f; }
				
    virtual void SetCallback( IFaIndexerCallback * cbk ) { m_Cbk = cbk; }
    
    virtual ~AFaIndexerBase() throw() {}
    AFaIndexerBase() : m_File( 0 ), m_Count( 0 ), m_Cbk( 0 ), m_Flags( 0 ) {}
protected:
    // writes default header
    virtual void AttachHeaderDefault( const string& path, Uint4 version );

    virtual void AddSequence( const char * seq, unsigned len, off64_t off ) = 0;
    virtual void DumpTables() = 0;
    
protected:
    FILE   * m_File;
    Uint4 m_Count;
    IFaIndexerCallback * m_Cbk;
    Uint4 m_Flags;
};

class CFaIndexer1 : public AFaIndexerBase
{
public:
    enum EFlags {
        fSkipRepeatitive = 0x01,
        fStatTable       = 0x02
    };
    
    void AttachFile( const string& path );
    void SetHash( const CHashSet& hs );

    void SetFragmentSizeRange( unsigned lo, unsigned hi ) { m_Lo = lo; m_Hi = hi; }
    unsigned GetFragmentLo() const { return m_Lo; }
    unsigned GetFragmentHi() const { return m_Hi; }

    ~CFaIndexer1() throw() {}
    CFaIndexer1() : m_Lo( 512 * 1024 * 1024 ), m_Hi( 1536 * 1024 * 1024 ) {}
protected:
    void AddSequence( const char * seq, unsigned len, off64_t off );
    void DumpTables();

protected:
    TData    m_Data;   // vector of hash tables (one per word)
    unsigned m_Lo, m_Hi;
};

class CFaIndexer2 : public AFaIndexerBase
{
public:
    void AttachFile( const string& path );
				
    void SetHash( const CHashSet& hs );

	void SetCacheSize( unsigned cs ) { m_CacheSize = cs; }

    ~CFaIndexer2() throw() {}
    CFaIndexer2() : m_CacheSize( 200000000 ) {}

protected:
    void AddSequence( const char * seq, unsigned len, off64_t off );
    void DumpTables();

	void StoreSequence( unsigned sid, const char * seq, unsigned len );
    void WriteCache();
protected:
    TStatVector m_Data;    // vector of hash lists
    TStatVector m_Seqs;    // vector of seq counts
    TStatVector m_LastSid; // vector of last sids
	TStatVector m_Cursor;  // vectpr of "current positions" for 2nd pass
    TData       m_Cache;
	unsigned    m_CacheSize; // cache size
};

class CFaLookup : public CFaData
{
public:
    CFaLookup() : m_Fd( -1 ), m_AlignL( 0 ), m_AlignR( 0 ) {}

    typedef list<ISts*> TStsList;
    
    enum EPrimerEnd { eLeft = 'l', eRight = 'r' };

    void AttachFile( const string& path );
    void Find( IFaLookupCallback * cbk,
               const string& label,
               char report_strand, // which strand to report
               const string& primer );

    void Find( IFaLookupCallback * cbk, 
               ISts* sts, int window = 0 );
    void Stat();
    TBigCount CalcStat( const char * primer, unsigned wd = 0 ) {
        m_Hash.Begin( primer ); 
        return m_Hash.Good( wd ) ? CalcStat( m_Hash.GetValue( wd ), wd ) : 0 ;
    }
    TBigCount CalcStat( THashElement hashval, unsigned wd = 0 );
    TBigCount GetStat( const char * primer, unsigned wd = 0 ) {
        m_Hash.Begin( primer ); 
        return m_Hash.Good( wd ) ? GetStat( m_Hash.GetValue( wd ), wd ) : 0 ;
    }
    TBigCount GetStat( THashElement hashval, unsigned wd = 0 ) {
        return m_Counts.size() ? m_Counts[ wd ][ hashval ] : 0 ;
    }
    const TStatVector& GetStat() const { return m_Counts; }
    
    void Find( IFaLookupCallback * cbk, 
               const TStsList& sts,
               bool syscall_optimize = true,
               int window = 0);

    void SetAligner( IAlign * left, IAlign * right ) {
        m_AlignL = left;
        m_AlignR = right;
    }

protected:
    off64_t GetHashEntries( const char * table, off64_t tab_off,
                            unsigned word, THashElement value,
                            unsigned& size ) const;
    void InitTableOffsets();
protected:
    int m_Fd;
    Uint4    m_Version;
    TStatVector m_Counts;
    vector<Uint4> m_TableOffset; // offsets of tables for words
    Uint4  m_ElSize;
    IAlign * m_AlignL, * m_AlignR;
};

struct SFaMatchBlock
{
    enum EStrand { eUnkn='0', ePos = '+', eNeg = '-' };
    enum EType { eSTS, ePrimer };
    typedef std::string string;

    EType    type;
    EStrand  strand;
    string   seq_label;
    string   sts_label;
    unsigned from, to;
    unsigned char mism, gaps;
    unsigned char mism_l, mism_r;
    unsigned char gaps_l, gaps_r;
    const char * sequence;
    unsigned     seqlen;
};

class IFaLookupCallback
{
public:
    virtual ~IFaLookupCallback() {}
    virtual bool Start() = 0;
    virtual bool Done() = 0;
    virtual bool Fail( const std::string& msg ) = 0;
    virtual bool Warn( const std::string& msg, 
                       const ISts * ) = 0;
    virtual bool Match( const SFaMatchBlock * info) = 0;
    virtual bool Match( const SFaMatchBlock * info, 
                        const ISts * sts) = 0;
    virtual void Fragment( unsigned i, unsigned total ) {};
    virtual void Progress( unsigned i, unsigned total ) {};
};
 
END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: fahash.hpp,v $
 * Revision 1.15  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.14  2004/06/07 16:24:56  rotmistr
 * Bug fixes to previos version.
 *
 * Revision 1.13  2004/06/03 23:37:19  rotmistr
 * New aligner added.
 *
 * Revision 1.12  2004/05/27 20:35:46  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.11  2004/04/28 14:35:35  rotmistr
 * hashfile ver2 build/search works now
 *
 * Revision 1.10  2004/04/27 00:01:54  rotmistr
 * Second version of reverse hash file started
 *
 * Revision 1.9  2004/04/01 05:57:52  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.8  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.7  2004/02/18 05:44:24  rotmistr
 * Changes in CGI: sort order, separate misalignments for l and r primers, reload button
 *
 * Revision 1.6  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.5  2004/01/07 16:57:42  rotmistr
 * Fragment size is now configurable.
 *
 * Revision 1.4  2004/01/06 21:54:19  rotmistr
 * Statistics for word repetitions API added
 *
 * Revision 1.3  2003/12/30 21:36:31  rotmistr
 * Syscall optimisation mode added.
 *
 * Revision 1.2  2003/12/23 21:30:50  rotmistr
 * - gaps/mismatches reporting
 * - lo/hi fixup
 * - reverse sts in re-PCR_main
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
