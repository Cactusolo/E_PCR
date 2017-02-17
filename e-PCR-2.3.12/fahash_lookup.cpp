/* $Id: fahash_lookup.cpp,v 1.24 2008/03/26 16:03:16 rotmistr Exp $
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

#include <epcr/fahash.hpp>
#include <epcr/bin-io.hpp>
#include "fahash_defines.h"

#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <set>
#include <map>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

#include "fahash_internal.hpp"

////////////////////////////////////////////////////////////////////////

#define DBG fprintf(stderr,__FILE__":%d\n",__LINE__)
#define HEX(a) fprintf(stderr,__FILE__":%d\t"#a"\t= 0x%08x\n",__LINE__,a)
#define INT(a) fprintf(stderr,__FILE__":%d\t"#a"\t= %d\n",__LINE__,a)
#define STR(a) fprintf(stderr,__FILE__":%d\t"#a"\t= %s\n",__LINE__,a.c_str())

class CHashListIterator
{
protected:
    Uint4 * cur, * end;
    Uint4 sid;
public:
    enum { 
        eHighBit = ( Uint4(1) << ( 8 * sizeof( Uint4 ) - 1 ) ), 
        eNullSid = ( 0 ^ eHighBit ), 
        eSidMask = ~eHighBit 
    };
	    
    CHashListIterator( void * ptr = 0, unsigned cnt = 0 ) {
        cur = (Uint4 *) ptr;
        end = cur + cnt;
        sid = eNullSid;
        if( cur && cur < end ) {
            while( cur < end && *cur & eHighBit ) sid = *cur++ & eSidMask;
            if( sid == eNullSid ) 
                throw runtime_error( "hash list format error: "
                                     "does not start with sid" );
        } 
    }
    unsigned GetSid() const { return sid; }
    operator bool () const { return cur && cur < end; }
    int operator * () const { return int(*cur); } // !!!
    CHashListIterator& operator ++ () {
        if( cur ) {
            if( cur < end ) ++cur;
            while( cur < end && *cur & eHighBit ) sid = *cur++ & eSidMask;
        }
        return *this;
    }
    CHashListIterator operator ++ (int) { 
        CHashListIterator i( *this ); ++*this;
        return i;
    }
    bool operator < ( const CHashListIterator& i ) const 
        { return sid < i.sid || ( sid == i.sid && cur < i.cur ); }
    bool operator > ( const CHashListIterator& i ) const 
        { return sid > i.sid || ( sid == i.sid && cur > i.cur ); }
    bool operator == ( const CHashListIterator& i ) const 
        { return sid == i.sid && cur == i.cur; }
    bool operator != ( const CHashListIterator& i ) const 
        { return sid != i.sid || cur != i.cur; }
    bool operator <= ( const CHashListIterator& i ) const 
        { return sid <= i.sid || ( sid == i.sid && cur <= i.cur ); }
    bool operator >= ( const CHashListIterator& i ) const 
        { return sid >= i.sid || ( sid == i.sid && cur >= i.cur ); }
};

class CHashCoIterator 
{
protected:
    CHashListIterator x, y, z;
    int lo, hi;
public:
    CHashCoIterator( CHashListIterator& _x, CHashListIterator& _y,
                     int l, int h ) : x( _x ), y( _y ), z( _x ), lo( l ), hi( h ) { fit(); }
    void report( const string& state, FILE * f = stdout ) {
        fprintf( stdout ,"#%s# x=%d:%d; y=%d:%d; z=%d:%d; "
                 "lo=%d; len=%d; hi=%d;\n",state.c_str(),
                 x.GetSid(), *x, y.GetSid(), *y, z.GetSid(), *z,
                 lo, *y - *x, hi );
    }
    void fit() {
        if( ! *this ) return;
        
        while( true ) {
            do {
                while( x.GetSid() < y.GetSid() ) if( ! ++x ) return;
                while( x.GetSid() > y.GetSid() ) if( ! ++y ) return;
            } while( x.GetSid() != y.GetSid() );
            
            while( (*y - *x) > hi) {
                if( ! ++x ) return;
                if( x.GetSid() != y.GetSid() ) goto nextsid;
            }
            if( ( *y - *x ) >= lo ) { z = x; return; } 
        nextsid: if( ! ++y ) return;
        }
    }
    operator bool () const { return x && y && z; }
    const CHashListIterator& GetX() const { return z; }
    const CHashListIterator& GetY() const { return y; }
    CHashCoIterator& IncX() { return ++*this; }
    CHashCoIterator& IncY() { if(++y) fit();  return *this; }
            
    CHashCoIterator& operator ++ () {
        if( ++z && z.GetSid() == y.GetSid() && ( *y - *z ) >= lo ) return *this; 
        if( ++y ) fit(); 
        return *this;
    }
};

////////////////////////////////////////////////////////////////////////

void CFaLookup::AttachFile( const string& path )
{
	m_Fd = open( path.c_str(), O_RDONLY | O_LARGEFILE );
	if( m_Fd == -1 ) SYSERROR( path );
	// TODO: here signature, version, byte order and similar 
	// stuff will be written
	char buff[16];
	read( m_Fd, buff, 16 ); // signature
	if( memcmp( buff, SIGNATURE, 16 ) ) 
		throw runtime_error( "Wrong signature in " + path );
    m_Version = Read<Uint4>( m_Fd );
	if( Read<Uint4>( m_Fd ) != BYTE_ORDER_WORD )
		throw runtime_error( "Wrong byte order in " + path );
	if( m_Version != FILE_VERSION && m_Version != FILE_VERSION2 )
		throw runtime_error( "Wrong version in " + path );
	if( Read<Uint4>( m_Fd ) != HEADER_SIZE )
		throw runtime_error( "Wrong header size order in " + path );
	Uint4 wdsize = Read<Uint4>( m_Fd );
	Uint4 period = Read<Uint4>( m_Fd );
	m_Hash = CHashSet( wdsize, period );
	off64_t last_off = lseek64( m_Fd, 0, SEEK_END );
	if( last_off == -1 ) SYSERROR( string( "seek-end(0) failed" ) );
	last_off = lseek64( m_Fd, -1 * int(sizeof(off64_t)), SEEK_CUR );
	if( last_off == -1 ) SYSERROR( string( "seek-cur(-8) failed" ) );
//	sizeof(off64_t)
	off64_t tail_off = Read<off64_t>( m_Fd );
// 	fprintf(stderr,
// 			"ATTACH: last_off  = %lld\n"
// 			"ATTACH: tail_off  = %lld\n",
// 			last_off,
// 			tail_off);
	
	if( last_off != off64_t(tail_off + 3 * sizeof( off64_t )) ) {
		throw runtime_error( "Wrong file trailer in " + path );
	}
	
	lseek64( m_Fd, tail_off, SEEK_SET );
	off64_t seq_tbl = Read<off64_t>( m_Fd );
	off64_t fil_tbl = Read<off64_t>( m_Fd );
	off64_t off_tbl = Read<off64_t>( m_Fd );

	lseek64( m_Fd, seq_tbl, SEEK_SET );
	m_Seqlst.resize( Read<Uint4>( m_Fd ) );
	for( unsigned i = 0; i < m_Seqlst.size(); ++i ) {
		Int2 file = Read<Int2>(m_Fd);
		m_Seqlst[i] = SSeqDesc( file, Read<off64_t>( m_Fd ) );
	}

	lseek64( m_Fd, fil_tbl, SEEK_SET );
	m_Fapath.resize( Read<Uint4>( m_Fd ) );
	for( unsigned i = 0; i < m_Fapath.size(); ++i ) {
		Int2 flags = Read<Int2>( m_Fd );
		m_Fapath[i] = SFileDesc( flags, Read<string>( m_Fd ) );
	}
	
	lseek64( m_Fd, off_tbl, SEEK_SET );
	m_Tabloc.resize( Read<Uint4>( m_Fd ) );
	for( unsigned i = 0; i < m_Tabloc.size(); ++i ) {
		m_Tabloc[i].first = Read<off64_t>( m_Fd );
		m_Tabloc[i].second = Read<off64_t>( m_Fd );
	}

    InitTableOffsets();
}
	
void CFaLookup::InitTableOffsets()
{
    m_TableOffset.resize( m_Hash.GetWordCount() + 1 );
    m_ElSize=
        m_Version == FILE_VERSION ? 2 * sizeof(Uint4) :
        m_Version == FILE_VERSION2 ? sizeof(off64_t) + sizeof(Uint4):
        0;
    unsigned i = 0, pos = 0;
    for( ; i < m_Hash.GetWordCount(); ++i ) {
        m_TableOffset[i] = pos;
        pos += m_Hash.GetTableSize( i );
    }
	m_TableOffset[i] = pos;
}

off64_t CFaLookup::GetHashEntries( const char * table, off64_t tab_off, 
                                   unsigned wd, THashElement val,
                                   unsigned& size ) const
{
    const char * entry = table + ( m_TableOffset[wd] + val ) * m_ElSize;
    
    switch( m_Version ) {
    case FILE_VERSION2: 
        size = ( (const Uint4*) entry )[2];
        return ( (const off64_t*) entry )[0];
    case FILE_VERSION:
        size = ( (const Uint4*) entry )[1];
        return ( (const Uint4*) entry )[0] + tab_off;
    default: throw logic_error( "Unknown file version" );
    }
}

void CFaLookup::Find( IFaLookupCallback * cbk,
					  const string& label,
					  const char report_strand,
					  const string& primer)
{
	vector<CFastaMap> famap( m_Fapath.size() );
	
	SFaMatchBlock info;
	info.type = SFaMatchBlock::ePrimer;
	info.strand =
		report_strand == '+' ? SFaMatchBlock::ePos :
		report_strand == '-' ? SFaMatchBlock::eNeg :
		SFaMatchBlock::eUnkn;
	info.sts_label = label;

	if( cbk && !cbk->Start() ) return;
	if( primer.length() < m_Hash.GetWordSize() ) {
		cbk && cbk->Fail( "primer " + label + " is too short" );
		return;
	}

    unsigned htbl_size = m_TableOffset[ m_Hash.GetWordCount() ] * m_ElSize;
	
	m_Hash.Begin( primer.c_str() ); 

    bool warned = false;
	for( unsigned frag = 0; frag < m_Tabloc.size(); ++frag ) {
		CMMap table( htbl_size,
                     CMMap::fProtRead, CMMap::fMapPrivate, m_Fd,
                     m_Tabloc[frag].first );
		
		for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
            if( cbk ) cbk->Progress( wd, m_Hash.GetWordCount() );
			if( ! m_Hash.Good( wd ) ) continue;
            
            unsigned count;
            off64_t start = GetHashEntries( table.data(), m_Tabloc[frag].first,
                                             wd, m_Hash.GetValue( wd ), count );
            if( count == ~0U ) {
                if( cbk && ! warned ) 
                    cbk->Warn( "Repeated word in primer " + label, 0 );
                warned = true;
                continue;
            } 

            CMMap xmap( count * sizeof(Uint4),
                        CMMap::fProtRead, CMMap::fMapPrivate, m_Fd,
                        start );
            madvise( xmap.data(), xmap.size(), MADV_SEQUENTIAL | MADV_DONTNEED );

			for( CHashListIterator u( xmap.data(), count ); u; ++u ) {
                Int4 pos1 = *u - m_Hash.GetWordSize();
                Int4 sid = u.GetSid();
                Int2 fil = m_Seqlst[sid].fid;
                
                if( ! famap[fil].IsOpen() ) 
                    famap[fil].Open( m_Fapath[fil].path );
            
                CMmSequence seq( famap[fil], m_Seqlst[sid].data );
                if( m_AlignR->Forward( seq.GetData() + pos1,
                                       seq.GetData() + seq.Length(),
                                       primer.c_str(), primer.length() ) ) {
                    info.seq_label = famap[fil].GetIdent( m_Seqlst[sid].data );
                    info.mism = m_AlignR->GetMismatches();
                    info.gaps = m_AlignR->GetGaps();
                    info.mism_l = info.mism_r = 0;
                    info.gaps_l = info.gaps_r = 0;
                    info.from = pos1;
                    info.to = info.from + primer.length();
                    info.sequence = seq.GetData();
                    info.seqlen = seq.Length();
                    if( cbk && ! ( cbk->Match( &info ) ) ) return;
				}
			}
		}
	}
	cbk && cbk->Done();
}

void CFaLookup::Find( IFaLookupCallback * cbk, ISts* sts, int window )
{
    TStsList lst;
    lst.push_back( sts );
    Find( cbk, lst, false, window );
}

struct SMatchPos
{
	unsigned sid;
	unsigned pos;
	SMatchPos( unsigned s, unsigned p ) : sid( s ), pos( p ) {}
	bool operator < ( const SMatchPos& s ) const {
		return sid < s.sid || sid == s.sid && pos < s.pos;
	}
};

struct SStsPreMatch 
{
    ISts * sts;
    int le, rb;
    SStsPreMatch( ISts * s, int p1, int p2 ) : sts( s ), le( p1 ), rb( p2 ) {}
};


class CMatchSts
{
public:    
    struct SHit 
    {
        int pos1, pos2;
        unsigned char mism_l, mism_r, gaps_l, gaps_r;
        int length() const { return pos2 - pos1; }

		bool operator == (const SHit& o) const {
			return 
                pos2   == o.pos2   && pos1   == o.pos1 && 
                mism_l == o.mism_l && gaps_l == o.gaps_l &&
                mism_r == o.mism_r && gaps_r == o.gaps_r;
            
		}
        SHit( int p1, int p2, char ml, char mr, char gl, char gr ) :
            pos1( p1 ), pos2( p2 ), mism_l( ml ), mism_r( mr ), gaps_l( gl ), gaps_r( gr ) {}
        SHit() {}
        unsigned char mism() const { return mism_l + mism_r; }
        unsigned char gaps() const { return gaps_l + gaps_r; }
        
        static bool OrderByPos2Pos1( const SHit&, const SHit& );
        static int  Compare( const SHit&, const SHit&, 
                             unsigned, unsigned);
        static bool Overlap( unsigned, unsigned, 
                             const SHit&, const SHit&);
        bool operator < ( const SHit& h )  const {
            return pos1 < h.pos1 || pos1 == h.pos1 && 
                 ( pos2 < h.pos2 || pos2 == h.pos2 &&
                 ( gaps() < h.gaps() || gaps() == h.gaps() &&
                   mism() < h.mism() ) );
        }
    };
    enum {
        fLeftNoMatch = 2,
        fRightNoMatch = 1,
    };
    
    typedef vector<SHit> TStsHits;
    typedef map<const ISts*,TStsHits> TAllHits;
    typedef map<int,TAllHits> TSidHits;
    
    typedef TStsHits::iterator TStsHits_I;
    typedef TAllHits::iterator TAllHits_I;
    typedef TSidHits::iterator TSidHits_I;

    typedef TStsHits::const_iterator TStsHits_CI;
    typedef TAllHits::const_iterator TAllHits_CI;
    typedef TSidHits::const_iterator TSidHits_CI;
public:
    ~CMatchSts() { delete m_Seq; }

    void Flush();
    
    CMatchSts( const CFaData::TPathLst& fapath, 
               const CFaData::TSeqLst&  seqlst,
               IFaLookupCallback * cbk, 
               IAlign * alignl, IAlign * alignr):
        m_Fapath( fapath ), m_Seqlst( seqlst ), m_Famap( fapath.size() ),
        m_OldSid( -1 ), m_Seq( 0 ), m_Cbk( cbk ), m_AlignL( alignl ), m_AlignR( alignr ) 
        { m_Info.type = SFaMatchBlock::eSTS; }
    // returns false if no processing should continue
    bool Match( ISts * sts, int sid, int posle, int posrb );
    int  MatchEx( ISts * sts, int sid, int posle, int posrb );
protected:
    TSidHits m_OutQueues;
    const CFaData::TPathLst&   m_Fapath;
    const CFaData::TSeqLst&    m_Seqlst;
	vector<CFastaMap> m_Famap;
    int m_OldSid;
    CMmSequence * m_Seq;
    IFaLookupCallback *m_Cbk;
//    int m_Mism, m_Gaps;
    IAlign * m_AlignL, * m_AlignR;
	SFaMatchBlock m_Info;
};

bool CMatchSts::SHit::OrderByPos2Pos1( const SHit& o1, const SHit& o2 )
{
    return ( o1.pos2 < o2.pos2 || o1.pos2 == o2.pos2 && o1.pos1 < o2.pos1 );
}


int CMatchSts::SHit::Compare( const SHit& a, const SHit& b, 
                              unsigned min_l, unsigned max_l ) 
{
    if( a.gaps() > b.gaps() ) return -1;
    if( a.gaps() < b.gaps() ) return +1;
    if( a.mism() > b.mism() ) return -1;
    if( a.mism() < b.mism() ) return +1;
    int lb = ( b.pos2 - b.pos1 );
    int la = ( a.pos2 - a.pos1 );
    int min_len = min_l;
    int max_len = max_l;
    if( la >= min_len && la <= max_len ) {
        if( lb >= min_len && lb <= max_len ) {
            return lb - la;
        } else return 1;
    } else {
        if( lb >= min_len && lb <= max_len ) {
            return -1;
        } else {
            int da = ( la < min_len ) ? min_len - la : la - max_len;
            int db = ( lb < min_len ) ? min_len - lb : lb - max_len;
            return db - da;
        }
    }
}

bool CMatchSts::SHit::Overlap( unsigned l1, unsigned l2, 
                               const SHit& a, const SHit& b )
{
    return ( (unsigned)abs( a.pos2 - b.pos2 ) <= l2 && 
             (unsigned)abs( a.pos1 - b.pos1 ) <= l1 );
}

void CMatchSts::Flush()
{
    for( TSidHits_I s = m_OutQueues.begin(); s != m_OutQueues.end(); ++s ) {
        int sid = s->first;
        int fil = m_Seqlst[sid].fid;
        m_Info.seq_label = m_Famap[fil].GetIdent( m_Seqlst[sid].data );
        CMmSequence seq( m_Famap[fil], m_Seqlst[sid].data );

        TAllHits& all = s->second;

        for( TAllHits_I i = all.begin(); i != all.end(); ++ i ) {
            // 1st: cluster 
            // 2nd: select best per cluster
            if( i->second.size() == 0 ) continue;
            
            int l1 = i->first->GetPrimerLength( 0 );
            int l2 = i->first->GetPrimerLength( 1 );
            int lo = i->first->GetSizeLo();
            int hi = i->first->GetSizeHi();
            
            TStsHits& hits = i->second;
            sort( hits.begin(), hits.end() );
        
            
            while( hits.size() ) {
                
                TStsHits todo;
                TStsHits cluster;
                cluster.push_back( hits.back() );
                SHit best = cluster.front();
                
                hits.pop_back();
                
                for( TStsHits_CI j = hits.begin(); j != hits.end(); ++j ) {
                    bool overlaps = false;
                    for( TStsHits_CI h = cluster.begin(); h != cluster.end(); ++h )
                        if( SHit::Overlap( l1, l2, *j, *h ) ) 
                            { overlaps = true; break; }
                    if( overlaps ) {
                        if( SHit::Compare( *j, best, lo, hi ) >= 0 ) best = *j;
                        cluster.push_back( *j );
                    } else {
                        todo.push_back( *j );
                    }
                }
                
                m_Info.strand =
                    i->first->GetDirection() == '+' ? SFaMatchBlock::ePos :
                    i->first->GetDirection() == '-' ? SFaMatchBlock::eNeg :
                    SFaMatchBlock::eUnkn;

                m_Info.sts_label = i->first->GetName();
                m_Info.from = best.pos1;
                m_Info.to   = best.pos2;
                m_Info.mism = best.mism();
                m_Info.gaps = best.gaps();
                m_Info.mism_l = best.mism_l;
                m_Info.mism_r = best.mism_r;
                m_Info.gaps_l = best.gaps_l;
                m_Info.gaps_r = best.gaps_r;
                m_Info.sequence = seq.GetData();
                m_Info.seqlen = seq.Length();
                m_Cbk->Match( &m_Info, i->first );
                
                hits = todo;
            }
        }
        all.clear();
    }
    m_OutQueues.clear();
}

bool CMatchSts::Match( ISts * sts, int sid, int posle, int posrb ) 
{
    CMatchSts::MatchEx( sts, sid, posle, posrb ) ;
    return true;
}

int CMatchSts::MatchEx( ISts * sts, int sid, int posle, int posrb ) 
{
    if( sid != m_OldSid ) {
//        Flush();
        int fil = m_Seqlst[sid].fid;
        if( ! m_Famap[fil].IsOpen() ) {
            m_Famap[fil].Open( m_Fapath[fil].path );
        }
        delete m_Seq;
        m_Seq = new CMmSequence( m_Famap[fil], m_Seqlst[sid].data );
        if( m_Seq->GetSize() < (unsigned)posrb || 
            m_Seq->GetSize() < (unsigned)posle ) {
            ostringstream o;
            o << "range error for hash word, sid=" << sid 
              << ": 0 < " << posle << " < " << posrb << " < " 
              << m_Seq->GetSize();
            
            throw range_error( o.str() );
        }
        m_OldSid = sid;
    }

    CStrRef lp( sts->GetPrimer( ISts::eLeft ) );
    CStrRef rp( sts->GetPrimer( ISts::eRight ) );
                    
    int poslb = posle-lp.length();
    int posre = posrb+rp.length();

    int flags = 0;
    
    if( ! m_AlignR->Forward( m_Seq->GetData() + posrb,
                             m_Seq->GetData() + m_Seq->Length(),
                             rp.data(), rp.length() ) ) flags |= fRightNoMatch;
    if( ! m_AlignL->Reverse( m_Seq->GetData(),
                             m_Seq->GetData() + posle,
                             lp.data(), lp.length())) flags |= fLeftNoMatch;
        
    if( ! flags ) {
        m_OutQueues[m_OldSid][sts].push_back(
            (sts->GetDirection() == ISts::ePlus)?
            SHit( poslb, posre,
                  m_AlignL->GetMismatches(),
                  m_AlignR->GetMismatches(),
                  m_AlignL->GetGaps(),
                  m_AlignR->GetGaps()) :
            SHit( poslb, posre,
                  m_AlignR->GetMismatches(),
                  m_AlignL->GetMismatches(),
                  m_AlignR->GetGaps(),
                  m_AlignL->GetGaps() )
            );
    }
    
    return flags;
} 


void CFaLookup::Stat()
{
    m_Counts.clear();
    
    m_Counts.resize( m_Hash.GetWordCount() );
    for( unsigned i = 0; i < m_Hash.GetWordCount(); ++i ) {
        m_Counts[i].resize( m_Hash.GetTableSize(i) );
    }
    unsigned htbl_size = m_TableOffset[m_Hash.GetWordCount()] * m_ElSize;
	
	for( unsigned frag = 0; frag < m_Tabloc.size(); ++frag ) {
        CMMap loc( htbl_size,
                   CMMap::fProtRead, CMMap::fMapPrivate, m_Fd,
                   m_Tabloc[frag].first );


        for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
            
            for( unsigned e = 0; e < m_Hash.GetTableSize(wd); ++e ) {
                unsigned count;
                GetHashEntries( loc.data(), m_Tabloc[frag].first, wd, e, count );

                if( count != ~0U ) m_Counts[wd][e] += count;
            }
        }
    }
}

CFaLookup::TBigCount CFaLookup::CalcStat( THashElement hash, unsigned wd )
{
    TBigCount cnt = 0;
    
    unsigned base = 0;
        
    for( unsigned i = 0; i < wd; ++i ) 
        base += 2 * sizeof(THashElement) * m_Hash.GetTableSize( wd );

    base += m_ElSize * hash;

	for( unsigned frag = 0; frag < m_Tabloc.size(); ++frag ) {
        CMMap loc( m_Tabloc[frag].second,
                   CMMap::fProtRead, CMMap::fMapPrivate, m_Fd,
                   m_Tabloc[frag].first );

        Uint4 * entry = (Uint4*)( loc.data() + base );
//        unsigned start=entry[0];
        Uint4 size = entry[1];
        if( size != ~0U ) cnt += size;
    }

    return cnt;
}

void CFaLookup::Find( IFaLookupCallback * cbk, const TStsList& lsts,
                      bool syscall_optimize, int window )
{
	if( cbk && ! cbk->Start() ) return;

    CMatchSts stsmatch( m_Fapath, m_Seqlst, cbk, m_AlignL, m_AlignR );
	
    for( TStsList::const_iterator i = lsts.begin(); i != lsts.end(); ++i ) {
        ISts * const& sts = *i;
        
        if( sts->GetPrimerLength( ISts::eLeft )  < m_Hash.GetWordSize() ||
            sts->GetPrimerLength( ISts::eRight ) < m_Hash.GetWordSize() ) {
            cbk && cbk->Fail( "sts " + string( sts->GetName() ) + " is too short\n" );
            return;
        }
    }
    
    map<ISts*,bool> warned;

	CHashSet r_Hash( m_Hash );
    unsigned htbl_size = m_TableOffset[m_Hash.GetWordCount()] * m_ElSize;
	
	unsigned overall = m_Tabloc.size() * lsts.size() * m_Hash.GetWordCount();
    unsigned progress = 0;
    
	for( unsigned frag = 0; frag < m_Tabloc.size(); ++frag ) {
        if(cbk) cbk->Fragment( frag, m_Tabloc.size() );
        
        CMMap table( htbl_size, //m_Tabloc[frag].second,
                     CMMap::fProtRead, CMMap::fMapPrivate, m_Fd,
                     m_Tabloc[frag].first );

        typedef multimap<int,SStsPreMatch> TPrematch;
		typedef TPrematch::const_iterator TPrematch_CI;
		TPrematch prematch;
		// walk around stupidestest MS Visual C++ 8.x (.NET) which can't convert const_iterator to iterator :-/
		const TPrematch& cprematch( prematch ); 

        for( TStsList::const_iterator i = lsts.begin(); i != lsts.end(); ++i ) {
             ISts * const& sts = *i;
            
            m_Hash.Begin( sts->GetPrimerData( ISts::eLeft ) +
                          sts->GetPrimerLength( ISts::eLeft ) -
                          m_Hash.GetWordSize() ); 
            r_Hash.Begin( sts->GetPrimerData( ISts::eRight ) ); 

            int fixup=
                sts->GetPrimerLength( ISts::eLeft ) +
                sts->GetPrimerLength( ISts::eRight ) -
                m_Hash.GetWordSize();
            
            int limitsize = max( sts->GetPrimerLength( ISts::eLeft ),
                                 sts->GetPrimerLength( ISts::eRight ) );
	    
            int lo = sts->GetSizeLo() - fixup - window;
            int hi = sts->GetSizeHi() - fixup + window;

            int minlo = limitsize - fixup;

            if( lo < minlo ) lo = minlo;
            if( hi < lo ) hi = lo;

            for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
                if( cbk ) cbk->Progress( progress++, overall );
                if( ! m_Hash.Good(wd) || ! r_Hash.Good(wd) ) continue;

                unsigned ucount, vcount;
                off64_t  ustart, vstart;
                {{ 
                    ustart = GetHashEntries( table.data(), m_Tabloc[frag].first,
                                             wd, m_Hash.GetValue(wd), ucount );
                    
                    vstart = GetHashEntries( table.data(), m_Tabloc[frag].first,
                                             wd, r_Hash.GetValue(wd), vcount );
                    

                    if( cbk && ! warned[sts] ) {
                        if( ucount == ~0U )
                            cbk->Warn( "STS " + string( sts->GetName() ) +
                                      " has repeatable word for left primer ",
                                      sts );
                        if( vcount == ~0U )
                            cbk->Warn( "STS " + string( sts->GetName() ) +
                                      " has repeatable word for right primer ",
                                      sts );
                        if( ucount == ~0U || vcount == ~0U ) {
                            warned[sts] = true;
                            continue;
                        }
                    }
                    if( ucount == ~0U || vcount == ~0U) continue;
                    if( ucount == 0   || vcount == 0) continue;
                }}
                
                CMMap umap( ucount * sizeof( Uint4 ),
                            CMMap::fProtRead, CMMap::fMapPrivate,
                            m_Fd, ustart );
                
                CMMap vmap( vcount * sizeof( Uint4 ),
                            CMMap::fProtRead, CMMap::fMapPrivate,
                            m_Fd, vstart );
                
                madvise( umap.data(), umap.size(), MADV_SEQUENTIAL | MADV_DONTNEED );
                madvise( vmap.data(), vmap.size(), MADV_SEQUENTIAL | MADV_DONTNEED );

                CHashListIterator u( umap.data(), ucount );
                CHashListIterator v( vmap.data(), vcount );
                
                for( CHashCoIterator i( u, v, lo, hi ); i; ) {
                    int posrb = *i.GetY() - m_Hash.GetWordSize();
                    int posle = *i.GetX();
                    int sid = i.GetX().GetSid();
                    if( (Uint4)sid != i.GetY().GetSid()) {
                        throw runtime_error( "CoIterator error!" );
                    }
                    
                    if( syscall_optimize ) {
                        
                        prematch.insert(
                            make_pair( sid, SStsPreMatch( sts, posle, posrb ) ) );
                            
                        if( prematch.size() > 10000 ) {
                            
                            for(TPrematch_CI im = cprematch.begin(); im != cprematch.end(); ++im) {
//                                    int sid=im->first;
                                if( ! stsmatch.Match(im->second.sts,
                                                     im->first,
                                                     im->second.le,
                                                     im->second.rb ) )
                                    return;
                            }
                            prematch.clear();
                        }
                    } else {
                        int rc = stsmatch.MatchEx( sts, sid, posle, posrb );
                        if(rc) {
                            if( rc & CMatchSts::fRightNoMatch )  i.IncY();
                            else  if( rc & CMatchSts::fLeftNoMatch ) i.IncX();
                            continue;
                        }
                    }
                    
                    ++i;
                }
            }
        }

        for( TPrematch_CI im = cprematch.begin(); im != cprematch.end(); ++im ) {

            if( ! stsmatch.Match( im->second.sts, im->first,
                               im->second.le, im->second.rb ) )
                return;

        }
        stsmatch.Flush();
	}
    
	cbk && cbk->Done();
}
	
////////////////////////////////////////////////////////////////////////



/*
 * $Log: fahash_lookup.cpp,v $
 * Revision 1.24  2008/03/26 16:03:16  rotmistr
 * Fixed bug with false negatives in unoptimized mode
 *
 * Revision 1.23  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.22  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.21  2005/02/11 21:53:50  rotmistr
 * Fixed one more "margin" bug
 *
 * Revision 1.20  2005/02/11 20:42:54  rotmistr
 * Fixed "margin" bug, added primer search from file
 *
 * Revision 1.19  2004/09/03 21:28:49  rotmistr
 * Fixes to compile with Borland C++ 5.5
 *
 * Revision 1.18  2004/06/07 16:24:56  rotmistr
 * Bug fixes to previos version.
 *
 * Revision 1.17  2004/06/03 23:37:20  rotmistr
 * New aligner added.
 *
 * Revision 1.16  2004/05/27 21:12:46  rotmistr
 * Some warnings fixed.
 *
 * Revision 1.15  2004/05/27 20:35:47  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.14  2004/04/28 14:35:36  rotmistr
 * hashfile ver2 build/search works now
 *
 * Revision 1.13  2004/04/27 00:01:55  rotmistr
 * Second version of reverse hash file started
 *
 * Revision 1.12  2004/04/01 05:57:52  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.11  2004/02/18 05:44:25  rotmistr
 * Changes in CGI: sort order, separate misalignments for l and r primers, reload button
 *
 * Revision 1.10  2004/02/12 21:38:20  rotmistr
 * Fixed typo in seqcmp
 * Optimized and fixed lookup
 * Better look for reverse.cgi
 *
 * Revision 1.9  2004/02/11 04:34:55  rotmistr
 * Optimised lookup speed and memory usage
 * Fixed bug with end of sequence in stsmatch
 * Changing CGI look
 *
 * Revision 1.8  2004/02/04 21:37:40  rotmistr
 * Optimized filter for labeling sequences
 *
 * Revision 1.7  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.6  2004/01/28 23:27:02  rotmistr
 * "Best of overlapping" hit selection postprocessor added.
 *
 * Revision 1.5  2004/01/07 16:57:42  rotmistr
 * Fragment size is now configurable.
 *
 * Revision 1.4  2004/01/06 21:54:19  rotmistr
 * Statistics for word repetitions API added
 *
 * Revision 1.3  2003/12/30 21:36:32  rotmistr
 * Syscall optimisation mode added.
 *
 * Revision 1.2  2003/12/23 21:30:50  rotmistr
 * - gaps/mismatches reporting
 * - lo/hi fixup
 * - reverse sts in re-PCR_main
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
