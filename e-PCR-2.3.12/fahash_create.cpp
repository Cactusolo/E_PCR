/* $Id: fahash_create.cpp,v 1.14 2007/07/11 20:49:29 rotmistr Exp $
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

void AFaIndexerBase::AddFile( const string& path )
{
	CFastaMap fmap( path );
	m_Fapath.push_back( SFileDesc( 0, path ) );

    if( m_Cbk ) m_Cbk->CbkFile( path.c_str() ); 

	for( unsigned i = 0; i < fmap.SequenceCount(); ++i ) {
        if( m_Cbk ) m_Cbk->CbkSequence( fmap.GetIdent(i).c_str() );

		CMmSequence seq( fmap, i );
		AddSequence( seq.GetData(), seq.GetSize(), i );
	}
}

void AFaIndexerBase::AttachHeaderDefault( const string& path, Uint4 version )
{
	m_File = fopen64( path.c_str(), "w" FILE_BINARY );
	if( m_File == 0 ) SYSERROR( path );

	setvbuf( m_File, 0, _IOFBF, 1024 * 1024 );
	
	if( fwrite( SIGNATURE, 1, 16, m_File ) != 16 ) SYSERROR( path ); // signature
	Write( m_File, Uint4( version ) );                // version
	Write( m_File, Uint4( BYTE_ORDER_WORD ) );        // architecture 
	Write( m_File, Uint4( HEADER_SIZE ) ) ;           // header size
	Write( m_File, Uint4( m_Hash.GetWordSize() ) );   // word size
	Write( m_File, Uint4( m_Hash.GetWordCount() ) );  // period 
	fseeko64( m_File, HEADER_SIZE, SEEK_SET );
	m_Count = 0;
}

void AFaIndexerBase::Finish()
{
	DumpTables();

	off64_t seq_tbl = SeekAlign( m_File );
	Write( m_File, Uint4( m_Seqlst.size() ) );
	for( unsigned i = 0; i < m_Seqlst.size(); ++i ) {
		Write( m_File, m_Seqlst[i].fid );
		Write( m_File, m_Seqlst[i].data );
	}
	off64_t fil_tbl = SeekAlign( m_File );
	Write( m_File, Uint4( m_Fapath.size() ) );
	for( unsigned i = 0; i < m_Fapath.size(); ++i ) {
		Write( m_File, m_Fapath[i].flags );
		Write( m_File, m_Fapath[i].path );
	}
	off64_t off_tbl = SeekAlign( m_File );
	Write( m_File, Uint4( m_Tabloc.size() ) );
	for( unsigned i = 0; i < m_Tabloc.size(); ++i ) {
		Write( m_File, m_Tabloc[i].first );
		Write( m_File, m_Tabloc[i].second );
	}
	off64_t endoffile = SeekAlign( m_File );
	Write( m_File, seq_tbl );
	Write( m_File, fil_tbl );
	Write( m_File, off_tbl );
	Write( m_File, endoffile );
	fclose( m_File );
	m_File = 0;
}


/*
 * $Log: fahash_create.cpp,v $
 * Revision 1.14  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.13  2004/05/27 20:35:46  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.12  2004/04/27 00:01:54  rotmistr
 * Second version of reverse hash file started
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
