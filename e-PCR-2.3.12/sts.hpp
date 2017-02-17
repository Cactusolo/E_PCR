/* $Id: sts.hpp,v 1.2 2003/12/23 21:30:50 rotmistr Exp $
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

#ifndef EPCR_STS__HPP
#define EPCR_STS__HPP

#include <epcr/sts_i.hpp>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class CSimpleSTS;

class CSimpleSTS:public ISts
{
public:
    ~CSimpleSTS() throw () {}
    CSimpleSTS(const string& l, const string& r, 
               int a, int b=0, ISts::EDirect s=ISts::ePlus,
               const string& n="*", const string& d="...")
        :m_name(n),m_desc(d),m_left(l),m_right(r),
         m_lo(a),m_hi(b>a?b:a),m_strand(s) {}
				
    virtual const char * GetPrimerData(int s) const 
        { return s==eLeft?m_left.c_str():m_right.c_str(); }
    virtual unsigned     GetPrimerLength(int s) const 
        { return s==eLeft?m_left.length():m_right.length(); }
    virtual EDirect      GetDirection() const 
        { return m_strand; }
    virtual unsigned     GetSizeLo() const 
        { return m_lo; }
    virtual unsigned     GetSizeHi() const 
        { return m_hi; }
    virtual CStrRef      GetName() const 
        { return CStrRef(m_name.c_str(),m_name.length()); }
    virtual CStrRef      GetDescription() const 
        { return CStrRef(m_desc.c_str(),m_desc.length()); }
protected:
    string m_name;
    string m_desc;
    string m_left;
    string m_right;
    int m_lo, m_hi;
    EDirect m_strand;
};
			
END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: sts.hpp,v $
 * Revision 1.2  2003/12/23 21:30:50  rotmistr
 * - gaps/mismatches reporting
 * - lo/hi fixup
 * - reverse sts in re-PCR_main
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
