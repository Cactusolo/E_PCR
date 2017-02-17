/* $Id: fahash_create1.cpp,v 1.2 2007/07/11 20:49:29 rotmistr Exp $
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
#include <string.h>
#include <errno.h>
#include <assert.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

#include "fahash_internal.hpp"

////////////////////////////////////////////////////////////////////////

void CFaIndexer1::AddSequence( const char * sequence, unsigned len, off64_t pos )
{
    Uint8 exp_count = len * m_Hash.GetWordCount();
    
    if( exp_count > m_Hi ) throw runtime_error( "Too large sequence" );
    if( ( m_Count + exp_count ) > m_Hi ) DumpTables();
    
    m_Seqlst.push_back( SSeqDesc( m_Fapath.size() - 1, pos ) );

	m_Hash.Begin( sequence );
		
	assert( m_Hash.GetWordCount() == m_Data.size() );
	for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		assert( m_Hash.GetTableSize( wd ) == m_Data[ wd ].size() );
		for( unsigned u = 0; u < m_Hash.GetTableSize( wd ); ++u ) {
					
            THashElement sid = ( m_Seqlst.size() - 1 ) | kHighBit;
            if( m_Data[ wd ][ u ].size() && 
                m_Data[ wd ][ u ][ m_Data[ wd ][ u ].size() - 1 ] & kHighBit ) 
                // subst last sid
                m_Data[ wd ][ u ][ m_Data[ wd ][ u ].size() - 1 ] = sid;
            else {
                m_Data[ wd ][ u ].push_back( sid );
                m_Count++;
            }
		}
	}
		
	for( ; m_Hash.GetPosition() < len/*!m_Hash.End()*/; m_Hash.Next() ) {
        if( m_Hash.GetPosition() % 100 == 0 && m_Cbk ) m_Cbk->CbkProgress( m_Hash.GetPosition(), len );
    
		for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
			if( m_Hash.Good(wd) ) {
							
                m_Data[ wd ][ m_Hash.GetValue(wd) ].push_back( THashElement( m_Hash.GetPosition() ) );
                m_Count++;

			}
		}
	}
	if( m_Count > m_Lo ) DumpTables();
}

void CFaIndexer1::AttachFile( const string& path )
{
    AttachHeaderDefault( path, FILE_VERSION );
}

void CFaIndexer1::DumpTables()
{
//	INDICATE("DUMPING...");
    if( m_Cbk ) m_Cbk->CbkDumpStart();
	off64_t start_off = SeekAlign( m_File );
	Uint4 pos = 0;
	for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		pos += m_Hash.GetTableSize(wd) * 2 * sizeof(Uint4);
	}
	Uint4 cur=pos;

    long double avg=0, dev=0;
    Uint8 cnt=0;

	if( m_Flags & fSkipRepeatitive ) {
        for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
            for( unsigned el = 0; el < m_Hash.GetTableSize( wd ); ++el, ++cnt ) {
                avg += m_Data[ wd ][ el ].size();
                dev += m_Data[ wd ][ el ].size() * m_Data[ wd ][ el ].size();
            }
        }
        if( cnt ) {
            avg /= cnt;
            dev /= cnt;
        }
        dev -= avg;
    }
    
	for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		for(unsigned el = 0; el < m_Hash.GetTableSize( wd ); ++el ) {
			Write<Uint4>( m_File, pos );
            if( m_Flags & fSkipRepeatitive && m_Data[wd][el].size() > avg + 3 * dev ) {
                Write<Uint4>( m_File, ~0 );
            } else {
                Write<Uint4>( m_File, m_Data[wd][el].size() );
            }
			pos += m_Data[wd][el].size() * sizeof( THashElement );
		}
	}
	for( unsigned wd = 0; wd < m_Hash.GetWordCount(); ++wd ) {
		for( unsigned el = 0; el < m_Hash.GetTableSize( wd ); ++el ) {
            if( m_Cbk ) m_Cbk->CbkDumpProgress( cur, pos );
            if( m_Flags & fSkipRepeatitive && m_Data[wd][el].size() > avg + 3 * dev ) {
            } else {
                if( fwrite( &m_Data[wd][el][0], sizeof(THashElement),
                             m_Data[wd][el].size(), m_File ) !=
                    m_Data[wd][el].size() ) SYSERROR( string( "Block write" ) );
                cur += m_Data[wd][el].size() * sizeof(THashElement);
            }
		}
	}
	off64_t end_off = ftello64(m_File);
	
    if( m_Cbk ) m_Cbk->CbkDumpEnd();
// 	fprintf(stderr,
// 			"DUMP: start offset    = %lld\n"
// 			"DUMP: end   offset    = %lld\n"
// 			"DUMP: distance        = %lld\n"
// 			"DUMP: pos calculated  = %u\n"
// 			"DUMP: pos accumulated = %u\n",
// 			start_off, end_off, end_off-start_off,pos, cur);
	
	if( end_off - start_off != pos || pos != cur ) 
		throw runtime_error( "File format error!" );

    if( m_Cbk ) m_Cbk->CbkResetStart();
	m_Tabloc.push_back( make_pair( start_off, end_off - start_off ) );
//	INDICATE("RESETTING...");
	for( unsigned i = 0; i < m_Hash.GetWordCount(); ++i ) {
        for( unsigned j = 0; j < m_Hash.GetTableSize(i); ++j ) {
            m_Data[i][j].clear();
        }
    }
	m_Count = 0;
//	SetHash(m_Hash);
    if( m_Cbk ) m_Cbk->CbkResetEnd();
//	INDICATE("DUMPING...DONE");
}

void CFaIndexer1::SetHash( const CHashSet& hs )
{
// 	for(unsigned i=0; i<m_data.size();++i) {
// 		for(unsigned j=0; j<m_data[i].size();++j) {
// 			m_data[i][j].clear();
// 		}
// 		m_data[i].clear();
// 	}
 	m_Data.clear();
	m_Data.resize(hs.GetWordCount());
	for( unsigned i = 0; i < hs.GetWordCount(); ++i ) {
//		m_data[i].clear();
		m_Data[i].resize( hs.GetTableSize( i ) );
// 		for(unsigned j=0; j<m_data[i].size();++j) {
// 			m_data[i][j].clear();
// 		}
	}
	m_Hash = hs;
	m_Count = 0;
}

/*
 * $Log: fahash_create1.cpp,v $
 * Revision 1.2  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.1  2004/05/27 20:35:46  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 */


