/* $Id: fast_seqio.hpp,v 1.3 2004/04/06 04:53:17 rotmistr Exp $
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

#ifndef EPCR_FAST_SEQIO__HPP
#define EPCR_FAST_SEQIO__HPP

#include <epcr/build_cfg.h>
#include <epcr/faread.hpp>
#include <epcr/mmap.hpp>
#include <string>
#include <vector>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class CFastaMapData;
class CFastaMapPrepare;
class CFastaMap;
class CMmSequence;

class CFastaMapData
{
public:
    typedef off64_t TOffset;
    typedef off64_t TSize;

protected:
    vector<string>  m_Ident;
    vector<string>  m_Defline;
    vector<TOffset> m_Offset;
    vector<TSize>   m_Size;
};

class CFastaMapPrepare:public IFastaReaderCallback,protected CFastaMapData
{
public:
    virtual void CbkDefline(const char * defline, unsigned length);
    virtual void CbkIdent(const char * ident, unsigned length);
    virtual void CbkSeqline(const char * data, unsigned length);
    virtual void CbkEntryBegin();
    virtual void CbkEntryEnd();
    virtual void CbkFileBegin();
    virtual void CbkFileEnd();

    void Open(const string& fname);
    void Close();
    bool IsOpen() const { return m_Fptr!=0; }

    void AddFile(const string& fname, const char * cvtTable=0);

    virtual ~CFastaMapPrepare() throw () { Close(); }
    CFastaMapPrepare():m_Fptr(0) {}
    explicit CFastaMapPrepare(const string& fname):m_Fptr(0) { Open(fname); }
            
protected:
    void WritePrologue();
    void WriteEpilogue();

protected:
    
    void * m_Fptr;
};

class CFastaMap:public CFastaMapData
{
public:
    friend class CMmSequence;
				
    ~CFastaMap() throw ();
    CFastaMap();
    CFastaMap(const string& file);
				
    bool IsOpen() const { return m_Fd!=-1; }
				
    void Open(const string& file);
    void Close();
				
    unsigned SequenceCount() const { return m_Size.size(); }
				
    const string& GetDefline(unsigned i) const 
        { return m_Defline[i]; }
    const string& GetIdent(unsigned i) const 
        { return m_Ident[i]; } 
    unsigned      GetSize(unsigned i) const
        { return m_Size[i]; } 
    unsigned      GetCount() const 
        { return m_Size.size(); }
    
protected:
    void ReadPrologue();
    void ReadEpilogue();
				
protected:
    int m_Fd;
    bool m_SwapBytes;
};

class CIndexAssert
{
public:
    CIndexAssert(unsigned index, unsigned size, const string& msg="") {
        if(index >= size) throw range_error(msg);
    }
};
    
class CMmSequence:public CIndexAssert, public CMMap
{
public:
    explicit CMmSequence(CFastaMap&, unsigned);
    ~CMmSequence() throw ();
    operator const char * () const { return GetData(); }
//      const char * Data()  const { return GetData(); }
     unsigned     Length() const { return GetSize(); }
protected:
//     const char * m_Data;
//     unsigned m_Size;
};

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: fast_seqio.hpp,v $
 * Revision 1.3  2004/04/06 04:53:17  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.2  2003/12/30 15:27:22  rotmistr
 * Fixed bug with sequence end
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
