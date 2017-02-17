/* $Id: fast_seqio_write.cpp,v 1.5 2007/07/11 20:49:29 rotmistr Exp $
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

#include <epcr/fast_seqio.hpp>
#include <epcr/bin-io.hpp>

#include <stdexcept>
#include <string.h>
#include <errno.h>


USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

void CFastaMapPrepare::Open( const string& fname ) 
{
    if( IsOpen() ) Close();

    FILE * f = fopen64( fname.c_str(), "w"FILE_BINARY );
    if( f == 0 ) throw runtime_error( fname + ": " + strerror( errno ) );
    setvbuf( f, 0, _IOFBF, 16192 );
    m_Fptr = f;

    WritePrologue();
}

void CFastaMapPrepare::Close() 
{
    if( m_Fptr ) {
        WriteEpilogue();

        m_Offset.clear();
        m_Size.clear();
        m_Defline.clear();
        m_Ident.clear();

        fclose( (FILE*)m_Fptr );
        m_Fptr = 0;
    }
}

void CFastaMapPrepare::AddFile( const string& fname, const char * cvtTable )
{
    if( IsOpen() ) {
    
        CFastaReader reader( fname );
		if( cvtTable ) reader.SetCvtTable( cvtTable );
        reader.ReadFile( this );
    }
}

void CFastaMapPrepare::WritePrologue()
{
	if( sizeof( Int2 ) != 2 ) throw logic_error( "Bad compile options: sizeof( Int2 ) != 2" );
	if( sizeof( Int4 ) != 4 ) throw logic_error( "Bad compile options: sizeof( Int4 ) != 4" );
	if( sizeof( Int8 ) != 8 ) throw logic_error( "Bad compile options: sizeof( Int8 ) != 8" );
    
	if( fwrite( "FASTAMAP", 1, 8, (FILE*)m_Fptr ) != 8 ) 
		throw runtime_error( strerror( errno ) +
							 string(" while writing signature") );
	Write( (FILE*)m_Fptr, Uint4( eLoEndian ) );
	Write( (FILE*)m_Fptr, Uint2( 1 ));
	Write( (FILE*)m_Fptr, Uint2( 0 ));
	Write( (FILE*)m_Fptr, Uint4( 20 )); // header size
}

void CFastaMapPrepare::WriteEpilogue()
{
	TOffset epilogue = SeekAlign( (FILE*)m_Fptr );
	fwrite( "EPILOGUE", 1, 8, (FILE*)m_Fptr );
	Write( (FILE*)m_Fptr, Uint4( m_Offset.size() ) );
	Write( (FILE*)m_Fptr, m_Offset[0], m_Offset.size() );
	Write( (FILE*)m_Fptr, m_Size[0], m_Size.size() );

	TOffset ident = ftello64( (FILE*)m_Fptr );
	for( unsigned i = 0; i < m_Ident.size(); ++i ) 
        Write( (FILE*)m_Fptr, m_Ident[i] );
	
	TOffset defline = ftello64( (FILE*)m_Fptr );
	for( unsigned i = 0; i < m_Defline.size(); ++i ) 
        Write( (FILE*)m_Fptr, m_Defline[i] );
	
	TOffset directory = SeekAlign( (FILE*)m_Fptr );
	Write( (FILE*)m_Fptr, epilogue );
	Write( (FILE*)m_Fptr, ident );
	Write( (FILE*)m_Fptr, defline );
	Write( (FILE*)m_Fptr, directory );
}

void CFastaMapPrepare::CbkDefline( const char * defline, unsigned length )
{
    m_Defline.back().assign( defline, length );
}

void CFastaMapPrepare::CbkIdent( const char * ident, unsigned length )
{
    m_Ident.back().assign( ident, length );
}
                                
void CFastaMapPrepare::CbkSeqline( const char * seqline, unsigned length )
{
    Write( (FILE*)m_Fptr, *seqline, length );
    m_Size.back() += length;
}

void CFastaMapPrepare::CbkEntryBegin()
{
    m_Defline.push_back( string("") );
    m_Ident.push_back( string("") );
    m_Size.push_back( 0 );
    m_Offset.push_back( SeekAlign( (FILE*)m_Fptr, 2 ) );
}

void CFastaMapPrepare::CbkEntryEnd()
{
    Write( (FILE*)m_Fptr, char(0) );
}

void CFastaMapPrepare::CbkFileBegin()
{}

void CFastaMapPrepare::CbkFileEnd()
{}

/*
 * $Log: fast_seqio_write.cpp,v $
 * Revision 1.5  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.4  2004/04/06 04:53:17  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.3  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.2  2003/12/30 15:27:22  rotmistr
 * Fixed bug with sequence end
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
