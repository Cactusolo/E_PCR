/* $Id: hashset.hpp,v 1.3 2007/07/11 20:49:29 rotmistr Exp $
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

#ifndef EPCR_HASHSET__HPP
#define EPCR_HASHSET__HPP

#include <epcr/build_cfg.h>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

    class CHashSet;


class CHashSet
{
public:
    // stl-style typedefs
    typedef Uint4 hash_type;
    typedef Uint4 size_type;
    typedef Uint4 bits_type;
	
    // ncbi-style typedefs
    typedef hash_type THashValue;
    typedef size_type TSize;
    typedef bits_type TBitMask;

    TSize GetWordSize() const { return m_WdSize; }
    TSize GetWordCount() const { return m_Period; }
    TBitMask GetWord( TSize id ) const { return m_Word[id]; }
    TBitMask GetMask( TSize id ) const { return m_Mask[id]; }
    THashValue GetValue( TSize id ) const { return m_Hash[id]; }
    THashValue GetTableSize( TSize id ) const 
        { return m_TbSize[id]; }
    THashValue GetAmbiquityMask() const { return m_AmbMask; }
	
    THashValue operator [] ( TSize id ) const 
        { return GetValue( id ); }

    bool Begin( const char * ptr );	
    bool Next() { // inline
        if( End() ) return false;
					
        m_AmbMask <<= 1;
					
        unsigned char nt = sm_HashId[ int(*m_Ptr) ];
        if( nt > 3 ) /* ambiquos */ { m_AmbMask |= 1; nt = 0; }
					
        // Circular 'rotation' with update of values
        if( m_Period == 1 ) {
            m_Hash[0] = m_Mask[0] & ( ( m_Hash[0] << 2 ) | nt );
        } else {
            hash_type x = m_Hash[ m_Period - 1 ];
            for( size_type t = m_Period - 1; t > 0; --t ) {
                m_Hash[t] = m_Mask[t] & ( ( m_Hash[ t - 1 ] << 2 ) | nt );
            }
            m_Hash[0] = x & m_Mask[0];
        }
					
        ++m_Offset;
        ++m_Ptr;
					
        return true;
    };
    bool Good() const;
    bool Good( TSize id ) const { return ( m_Word[id] & m_AmbMask ) == 0; }
    bool End() const { return m_Ptr == 0 || *m_Ptr == 0; }
	
    TSize GetPosition() const { return m_Offset; }
    const char * GetPtr() const { return m_Ptr; }

    virtual ~CHashSet() throw ();

    CHashSet( TSize wdsize = 0, TSize period = 0 );
    CHashSet( const CHashSet& );

    CHashSet& operator = ( const CHashSet& s );
protected:

    void _intl_free();
    void _intl_copy( const CHashSet& s );
	
protected:
    const char * m_Ptr;     // Pointer to current position 
    // of string being processed
    TSize        m_Offset;  // Current offset in original string

    TSize        m_WdSize;  // Wordsize
    TSize        m_Period;  // Number of words to be used
    TSize      * m_TbSize;  // Hash table sizes for each word; 
    // used by users

    TBitMask     m_AmbMask; // Mask showing ambiquity positions
    TBitMask   * m_Word;    // Words; needed for ambiquity test
	
    THashValue * m_Hash;    // Current hash values for each word
    THashValue * m_Mask;    // Masks for each word; 
    // used for value update
protected:
    static unsigned char sm_HashId[];
};

		
END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: hashset.hpp,v $
 * Revision 1.3  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.2  2004/02/04 21:23:22  rotmistr
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
