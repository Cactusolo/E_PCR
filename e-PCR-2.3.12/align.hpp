/* $Id: align.hpp,v 1.5 2008/03/26 16:04:29 rotmistr Exp $
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

#ifndef EPCR_ALIGN__HPP
#define EPCR_ALIGN__HPP

#include <epcr/build_cfg.h>
#include <epcr/minilcs.hpp>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class IAlign;
class CAlignExact;
class CAlignNoGaps;
class CAlignFast;
class CAlignLCS;
class CAlignCompromise;

char * FlipSequence(const char * seq, unsigned len= ~0U);
char * UCaseSequence(const char * seq, unsigned len= ~0U);

class IAlign
{
public:        
    virtual ~IAlign() throw () {};
    
    virtual bool Forward(const char * seq_ptr,
                         const char * seq_end,
                         const char * primer,
                         int length) = 0;
    virtual bool Reverse(const char * seq_start,
                         const char * seq_ptr,
                         const char * primer,
                         int length) = 0;
    virtual int GetIdentities() const { return 0; }
    virtual int GetMismatches() const { return 0; }
    virtual int GetGaps() const { return 0; }
	virtual int GetFloppy() const { return 0; }
};

class CAlignExact:public virtual IAlign
{
public:
    CAlignExact():m_Identities(0) {}
    
    virtual ~CAlignExact() throw () {};
    
    bool Forward(const char * seq_ptr,
                 const char * seq_end,
                 const char * primer,
                 int length);

    bool Reverse(const char * seq_start,
                 const char * seq_ptr,
                 const char * primer,
                 int length);
    int GetIdentities() const { return m_Identities; }
protected:
    int m_Identities;
};
 
class CAlignNoGaps:public CAlignExact
{
public:
    virtual ~CAlignNoGaps() throw () {};
    
    bool Forward(const char * seq_ptr,
                 const char * seq_end,
                 const char * primer,
                 int length);

    bool Reverse(const char * seq_start,
                 const char * seq_ptr,
                 const char * primer,
                 int length);

    CAlignNoGaps(int mm=1):m_MaxMismatch(mm),m_Mism(0) {}

    int GetMismatches() const { return m_Mism; }
    
protected:
    bool x_Compare(const char * a, const char * b, int l);
    
protected:
    int m_MaxMismatch;
    int m_Mism;
};
 
class CAlignFast:public CAlignNoGaps
{
public:
    virtual ~CAlignFast() throw () {};
    
    bool Forward(const char * seq_ptr,
                 const char * seq_end,
                 const char * primer,
                 int length);
    bool Reverse(const char * seq_start,
                 const char * seq_ptr,
                 const char * primer,
                 int length);
    CAlignFast(int mm=1, int gg=1):CAlignNoGaps(mm),m_MaxGaps(gg),m_Gaps(0) {}

    int GetGaps() const { return m_Gaps; }
    
protected:
    int m_MaxGaps;
    int m_Gaps;
};
 
class CAlignLCS:public virtual IAlign
{
public:
    virtual ~CAlignLCS() throw () {};
    
    bool Forward(const char * seq_ptr,
                 const char * seq_end,
                 const char * primer,
                 int length);
    bool Reverse(const char * seq_start,
                 const char * seq_ptr,
                 const char * primer,
                 int length);
    CAlignLCS(int mm=1, int gg=1);

    int GetIdentities() const { return m_Matrix.GetMatches(); }
    int GetMismatches() const { return m_Matrix.GetMismatches(); }
    int GetGaps() const { return m_Matrix.GetGaps(); }
    
protected:
    int m_MaxMismatch;
    int m_MaxGaps;
    CLcsMatrix<char> m_Matrix;
};

class CAlignCompromise:public virtual CAlignLCS, public virtual CAlignNoGaps
{
public:
    virtual ~CAlignCompromise() throw () {};
    
    bool Forward(const char * seq_ptr,
                 const char * seq_end,
                 const char * primer,
                 int length);
    bool Reverse(const char * seq_start,
                 const char * seq_ptr,
                 const char * primer,
                 int length);
    CAlignCompromise(int mm=1, int gg=1);

    int GetIdentities() const { return m_Identities; }
    int GetMismatches() const { return m_Mismatches; }
    int GetGaps() const { return m_Gaps; }
    
protected:
    int m_Identities, m_Mismatches, m_Gaps;
};

//class CAlignFloppy : public virtual CAlignLCS

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: align.hpp,v $
 * Revision 1.5  2008/03/26 16:04:29  rotmistr
 * Added support for blastdb files
 *
 * Revision 1.4  2004/10/26 17:16:33  rotmistr
 * Added 5'-end masking for primers
 *
 * Revision 1.3  2004/06/08 20:32:50  rotmistr
 * Fixup for gap+insert special case
 *
 * Revision 1.2  2004/06/08 16:14:55  rotmistr
 * *** empty log message ***
 *
 * Revision 1.1  2004/06/03 23:37:19  rotmistr
 * New aligner added.
 *
 */

