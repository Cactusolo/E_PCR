/* $Id: bin-io.hpp,v 1.5 2008/04/28 16:38:45 rotmistr Exp $
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

#ifndef EPCR_BIN_IO__HPP
#define EPCR_BIN_IO__HPP

#include <epcr/build_cfg.h>

#include <stdio.h>
#include <errno.h>

#include <stdexcept>
#include <cstring>
#include <string>
#include <vector>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

// This code is being used INTERNALLY by e-PCR library

enum EByteOrder { eHiEndian = 0x78563412, eLoEndian = 0x12345678 };

template <class T>
inline T BoCvt( T x, bool do_swap ) 
{
    if( do_swap ) {
        unsigned char * a = (unsigned char*) &x;
        unsigned char * b = a + sizeof( x ) - 1;
					
        for( ; a < b ; ++a, --b ) swap( *a, *b );
    }
    return x;
}

template<class T>
inline void Write( int fd, const T& t, unsigned sz = 1 )
{
    if( (size_t) write( fd, &t, sz * sizeof(T) ) != sz * sizeof( T ) )
        throw runtime_error( "write failed: " + string( std::strerror( errno ) ) );
}

template<>
inline void Write<string>( int fd, const string& t, unsigned )
{
    Write( fd, Uint4( t.length() ) );
    if( (size_t)write( fd, t.data(), t.length() ) != t.length() )
        throw runtime_error( "write failed: " + string( strerror( errno ) ) );
}

template<class T>
inline void Write( FILE* f, const T& t, unsigned sz = 1 )
{
    if( (size_t)fwrite( &t, sizeof( T ), sz, f ) != sz )
        throw runtime_error( "write failed: " + string( strerror( errno ) ) );
}

template<>
inline void Write<string>(FILE* f, const string& t, unsigned )
{
    Write( f, Uint4( t.length() ) );
    if((size_t)fwrite( t.data(), 1, t.length(), f ) != t.length() )
        throw runtime_error( "write failed: " + string( strerror( errno ) ) );
}

template<class T>
inline T Read( int fd )
{
    T t( -1 );
    if( (size_t)read( fd, &t, sizeof( T ) ) != sizeof( T ) && errno )
        throw runtime_error( "read failed: " + string( strerror( errno ) ) );
    return t;
}

template<>
inline string Read<string>( int fd )
{
    vector<char> t( Read<Uint4>( fd ) );
    if( t.size() && (size_t)read( fd, &t[0], t.size() ) != t.size() && errno )
        throw runtime_error( "read failed: " + string( strerror( errno ) ) );
    return string( &t[0], t.size() );
}
			
inline off64_t SeekAlign( int fd, unsigned page = 4096 ) 
{
    off64_t        c = lseek64( fd, 0, SEEK_CUR );
    if( c % page ) c = lseek64( fd, page - c % page, SEEK_CUR );
    return c;
}

inline off64_t SeekAlign( FILE* f, unsigned page = 4096 ) 
{
    off64_t    c = ftello64( f );
    if( c % page ) {
        fseeko64( f, page - c % page, SEEK_CUR );
        c = ftello64( f );
    }
    return c;
}

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: bin-io.hpp,v $
 * Revision 1.5  2008/04/28 16:38:45  rotmistr
 * Applied patch to build with gcc-4.3
 *
 * Revision 1.4  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.3  2004/04/01 05:57:52  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.2  2004/02/04 21:23:21  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
