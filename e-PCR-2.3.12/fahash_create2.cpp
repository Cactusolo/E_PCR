/* $Id: fahash_create2.cpp,v 1.6 2008/04/28 16:38:45 rotmistr Exp $
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

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cstring>
#include <memory>

#include <errno.h>
#include <assert.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

#include "fahash_internal.hpp"

////////////////////////////////////////////////////////////////////////

void CFaIndexer2::AddSequence( const char * sequence, unsigned len, off64_t pos )
{
    m_Seqlst.push_back( SSeqDesc( m_Fapath.size() - 1, pos ) );

	m_Hash.Begin( sequence );
		
	THashElement sid = ( m_Seqlst.size() - 1 ) | kHighBit;

	assert( m_Hash.GetWordCount() == m_Data.size() );
	for( ; m_Hash.GetPosition() < len/*!m_Hash.End()*/; m_Hash.Next() ) {
        if( m_Cbk && m_Hash.GetPosition() % 100 == 0 ) m_Cbk->CbkProgress( m_Hash.GetPosition(), len );
		for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
			assert( m_Hash.GetTableSize( wd ) == m_Data[wd].size() );
			if( m_Hash.Good( wd ) ) {
				THashElement val = m_Hash.GetValue( wd );
				if( m_LastSid[wd][val] != sid ) {
					m_LastSid[wd][val] = sid;
					m_Seqs[wd][val]++;
				}
				m_Data[wd][val]++;
			}
		}
	}
}

void CFaIndexer2::AttachFile( const string& path )
{
    AttachHeaderDefault( path, FILE_VERSION2 );
}

void CFaIndexer2::DumpTables()
{
    if( m_Cbk ) m_Cbk->CbkResetStart();

	// Allocate Space
	unsigned entries = 0;
	for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		entries += m_Hash.GetTableSize(wd);
	}

	off64_t start_off = SeekAlign( m_File );
	off64_t curr = start_off;

	TStatVector& pos = m_Cursor;
	pos.resize( m_Hash.GetWordCount() );
    m_Cache.resize( m_Hash.GetWordCount() );
	for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		pos[wd].resize( m_Hash.GetTableSize(wd) );
        m_Cache[wd].resize( m_Hash.GetTableSize(wd) );
		curr += m_Hash.GetTableSize( wd ) * ( sizeof(off64_t) + sizeof(THashElement) );
	}

	m_Tabloc.push_back( make_pair( start_off, curr - start_off ) );

	off64_t totalseq = 0, maxseq = 0, nullseq = 0, totalwd = 0, maxwd = 0;
	for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		for( unsigned el = 0; el < m_Hash.GetTableSize( wd ); ++el ) {
			// Set starting position
            Write<off64_t>( m_File, curr );
            Write<Uint4>( m_File, m_Seqs[wd][el] + m_Data[wd][el] );
            
			pos[wd][el] = curr;
			curr += ( m_Seqs[wd][el] + m_Data[wd][el] ) * sizeof(Uint4);
			// Reset lastsid
			m_LastSid[wd][el] = ~0U; 
			// Statistics
			totalseq += m_Seqs[wd][el];
			if( m_Seqs[wd][el] > maxseq ) maxseq = m_Seqs[wd][el];
			if( m_Seqs[wd][el] == 0 ) nullseq++;
			totalwd += m_Data[wd][el];
			if( m_Data[wd][el] > maxwd ) maxwd = m_Data[wd][el];
		}
	}

    if( m_Cbk ) m_Cbk->CbkResetEnd();

	// Print stat
	cerr << "Total number of sequences: " << m_Seqlst.size() << endl;

	cerr << "Number of hash entries: " 
		 << entries << endl;
	cerr << "Average number of sequences per hash entry: " 
		 << double(totalseq) / entries << endl;
	cerr << "Maximal number of sequences per hash entry: " 
		 << maxseq << endl;
	cerr << "Null hash entries: " 
		 << nullseq << endl;
	cerr << "Total number of hits: " 
		 << totalwd << endl;
	cerr << "Average number of hits per hash entry: " 
		 << double(totalwd) / entries << endl;
	cerr << "Maximal number of hits per hash entry: " 
		 << maxwd << endl;

	// Second pass
	auto_ptr<CFastaMap> fmap( 0 );
	Uint2 oldfid = ~Uint2(0);
    m_Count = 0;
    
	for( int sid = 0; sid < m_Seqlst.size(); ++sid ) {
		const SSeqDesc& sd = m_Seqlst[sid];
		if( oldfid != sd.fid ) {
			if( fmap.get() != 0 ) delete fmap.get();
			fmap.reset( new CFastaMap( m_Fapath[ oldfid = sd.fid ].path ) );
		}
        if( m_Cbk ) m_Cbk->CbkSequence( fmap->GetIdent( sd.data ).c_str() );
		CMmSequence seq( *fmap, sd.data );
		StoreSequence( sid, seq.data(), seq.size() );
	}
	WriteCache();
}

void CFaIndexer2::StoreSequence( unsigned sid, const char * seq, unsigned len )
{
	m_Hash.Begin( seq );
		
//	THashElement maxsid=(m_Seqlst.size()-1)|kHighBit;

	assert( m_Hash.GetWordCount() == m_Data.size() );

	for( ; m_Hash.GetPosition() < len/*!m_Hash.End()*/; m_Hash.Next() ) {
        if( m_Cbk && m_Hash.GetPosition() % 100 == 0 ) m_Cbk->CbkProgress( m_Hash.GetPosition(), len );
		for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
			if( m_Hash.Good(wd) ) {
				THashElement val = m_Hash.GetValue( wd );
				if( m_LastSid[wd][val] != sid ) {
					m_LastSid[wd][val] = sid;
                    m_Cache[wd][val].push_back( sid | kHighBit );
// 					Write<unsigned>(m_File,sid|kHighBit);
// 					m_Cursor[wd][val]+=sizeof(unsigned);
                    ++m_Count;
				}
                m_Cache[wd][val].push_back( m_Hash.GetPosition() );
// 				Write<unsigned>(m_File,m_Hash.GetPosition());
// 				m_Cursor[wd][val]+=sizeof(unsigned);
                ++m_Count;
			}
		}
        if( m_Count > m_CacheSize ) {
            // this value should be calculated based on statistics
            WriteCache();
            m_Count = 0;
        }
	}
}

void CFaIndexer2::WriteCache()
{
    if( m_Cbk ) m_Cbk->CbkDumpStart();
//    fprintf(stderr,"\nWriting cache...");
    for( unsigned wd = 0, cnt = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		TData::value_type & wdCache( m_Cache[wd] );
        for( unsigned val = 0; val < wdCache.size(); ++val ) {
			THashList & cache( wdCache[val] );
			
            if( cache.size() ) {
				Uint8 & cursor = m_Cursor[wd][val];
                if( fseeko64( m_File, cursor, SEEK_SET ) ) {
					ostringstream err;
					err << "Failed to seek to pos " << cursor 
						<< "(uint" << ( 8 * sizeof( cursor ) ) << "): "
						<< strerror(errno) << "\n";
					throw runtime_error( err.str() );
				}

                Write( m_File, cache[0], cache.size() );
                cursor += sizeof( cache[0] ) * cache.size();
				
				cnt += cache.size();
                cache.resize(0);
				
				if(m_Cbk) m_Cbk->CbkDumpProgress( cnt, m_Count );
            }
        }
    }
//    fprintf(stderr,"...Done\n");
    if( m_Cbk ) m_Cbk->CbkDumpEnd();

}

void CFaIndexer2::SetHash( const CHashSet& hs )
{
// 	for(unsigned i=0; i<m_data.size();++i) {
// 		for(unsigned j=0; j<m_data[i].size();++j) {
// 			m_data[i][j].clear();
// 		}
// 		m_data[i].clear();
// 	}
 	if( m_Data.size() ) m_Data.clear();
 	if( m_Seqs.size() ) m_Seqs.clear();
 	if( m_LastSid.size() ) m_LastSid.clear();
	m_Data.resize( hs.GetWordCount() );
	m_Seqs.resize( hs.GetWordCount() );
	m_LastSid.resize( hs.GetWordCount() );
	for( unsigned i = 0; i < hs.GetWordCount(); ++i ) {
//		m_data[i].clear();
		m_Data[i].resize( hs.GetTableSize(i) );
		m_Seqs[i].resize( hs.GetTableSize(i) );
		m_LastSid[i].resize( hs.GetTableSize(i) );
// 		for(unsigned j=0; j<m_data[i].size();++j) {
// 			m_data[i][j].clear();
// 		}
	}
	m_Hash = hs;
	m_Count = 0;
}


/*
 * $Log: fahash_create2.cpp,v $
 * Revision 1.6  2008/04/28 16:38:45  rotmistr
 * Applied patch to build with gcc-4.3
 *
 * Revision 1.5  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.4  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.3  2004/05/27 20:35:46  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.2  2004/04/28 14:35:36  rotmistr
 * hashfile ver2 build/search works now
 *
 * Revision 1.1  2004/04/27 00:03:59  rotmistr
 * Test files for second version of re-PCR
 *
 * Revision 1.11  2004/04/06 04:53:17  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.10  2004/04/01 05:57:52  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.9  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.8  2004/02/18 05:44:25  rotmistr
 * Changes in CGI: sort order, separate misalignments for l and r primers, reload button
 *
 * Revision 1.7  2004/02/12 21:38:20  rotmistr
 * Fixed typo in seqcmp
 * Optimized and fixed lookup
 * Better look for reverse.cgi
 *
 * Revision 1.6  2004/02/11 04:34:55  rotmistr
 * Optimised lookup speed and memory usage
 * Fixed bug with end of sequence in stsmatch
 * Changing CGI look
 *
 * Revision 1.5  2004/02/05 23:41:21  rotmistr
 * Better reload, fixed margin report in commandline, unists tab in CGI form.
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
 * Revision 1.3  2004/01/07 16:57:42  rotmistr
 * Fragment size is now configurable.
 *
 * Revision 1.2  2003/12/30 15:27:22  rotmistr
 * Fixed bug with sequence end
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
