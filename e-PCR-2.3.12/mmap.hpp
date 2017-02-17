/* $Id: mmap.hpp,v 1.4 2004/04/01 16:37:41 rotmistr Exp $
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

#ifndef EPCR_MMAP__HPP
#define EPCR_MMAP__HPP

#include <epcr/build_cfg.h>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

    class CMMap;

class CMMap
{
public:
    typedef off64_t TOffset;
    typedef int TFileHandle;

#ifndef USE_WIN     
    enum EProtFlags { 
        fProtRead  = PROT_READ,
        fProtWrite = PROT_WRITE,
        fProtExec  = PROT_EXEC,
        fProtNone  = PROT_NONE,
        fProtFlagsNONE = 0
    };
    
    enum EAccessFlags {
        fMapPrivate = MAP_PRIVATE,
        fMapShared  = MAP_SHARED,
        fMapNoReserve = MAP_NORESERVE,
        fMapAnon      = MAP_ANON,
        fAccessFlagsNONE = 0
    };
#else
    enum EProtFlags { 
        fProtRead  = 0x001,
        fProtWrite = 0x002,
        fProtExec  = 0x000,
        fProtNone  = 0x000,
        fProtFlagsNONE = 0
    };
    
    enum EAccessFlags {
        fMapPrivate = 0x001,
        fMapShared  = 0x002,
        fMapNoReserve = 0x004,
        fMapAnon      = 0x008,
        fAccessFlagsNONE = 0
    };

#endif
    explicit CMMap(unsigned size, unsigned prot, unsigned flags, 
                   TFileHandle fd, off64_t offset=0);
    ~CMMap() throw ();
    // stl-style methods
    char *   data() const { return m_Data; }
    unsigned size() const { return m_Size; }
    // ncbi-style methods
    char *   GetData() const { return m_Data; }
    unsigned GetSize() const { return m_Size; }
    // operators
    operator char *   () const { return m_Data; }
    operator unsigned () const { return m_Size; }

protected:
    char *   m_Data;
    unsigned m_Size;
    unsigned m_Delta;
#ifdef USE_WIN
	HANDLE m_Map;
#endif
};

		
END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: mmap.hpp,v $
 * Revision 1.4  2004/04/01 16:37:41  rotmistr
 * Cleaned after adding windows capabilities
 *
 * Revision 1.3  2004/04/01 05:57:53  rotmistr
 * Compilable with borland C++
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
