/* $Id: faread.hpp,v 1.5 2007/07/11 20:49:29 rotmistr Exp $
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

#ifndef EPCR_FAREAD__HPP
#define EPCR_FAREAD__HPP

#include <epcr/build_cfg.h>
#include <stdexcept>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class CFastaReader;
class IFastaReaderCallback;

class CFastaReader
{
public:
    virtual ~CFastaReader() throw () { Close(); }
    CFastaReader() : m_Fptr(0), m_CvtTable(0) {}
    explicit CFastaReader( const string& fname ):
        m_Fptr( 0 ), m_CvtTable( 0 ) { Open( fname ); }
    bool IsOpen() const { return m_Fptr; }
    void Open( const string& fname );
    void PipeIn( const string& cmd );
    void Close();
    virtual void ReadFile( IFastaReaderCallback* cbk );
//    virtual bool ReadEntry(IFastaReaderCallback* cbk);
    void SetCvtTable( const char * table = 0 ) { m_CvtTable = table; }
    
    // These tables should have '\0' in all positions to be "eaten" 
    // (f.e. for space and tab), everything else will be xlated
    static char sm_Nucleotides[];
    static char sm_NucleotidesUc[];
    static char sm_NucleotidesExt[];
    static char sm_NucleotidesExtUc[];
protected:
    void SeqlineCbk( IFastaReaderCallback * cbk, char * data, unsigned length );
    
protected:
    void * m_Fptr;
    const char * m_CvtTable;
};

class IFastaReaderCallback
{
public:
    virtual ~IFastaReaderCallback() throw () {}
    
    virtual void CbkDefline( const char * defline, unsigned length ) = 0;
    virtual void CbkIdent( const char * ident, unsigned length ) = 0;
    virtual void CbkSeqline( const char * data, unsigned length ) = 0;
    virtual void CbkEntryBegin() {}
    virtual void CbkEntryEnd() = 0;
    virtual void CbkFileBegin() {}
    virtual void CbkFileEnd() {}
};


END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: faread.hpp,v $
 * Revision 1.5  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.4  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
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
 * Revision 1.2  2004/01/08 23:22:41  rotmistr
 * Fixed init error in faread,
 * Adjusted output to standard,
 * Added output format style and output file to parameters.
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
