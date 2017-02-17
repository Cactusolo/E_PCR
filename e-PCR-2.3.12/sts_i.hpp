/* $Id: sts_i.hpp,v 1.2 2004/10/26 17:16:34 rotmistr Exp $
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

#ifndef EPCR_STS_I__HPP
#define EPCR_STS_I__HPP

#include <epcr/strref.hpp>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

	class ISts;

class ISts
{
public:
    enum ESide { eLeft=0, eRight=1 };
    enum EDirect { ePlus='+', eMinus='-' };
			
    CStrRef GetPrimer(ESide s) const 
        { return CStrRef(GetPrimerData(s),GetPrimerLength(s)); }
						
    virtual const char * GetPrimerData(int) const =0;
    virtual unsigned     GetPrimerLength(int) const =0;
    virtual EDirect      GetDirection() const =0;
    virtual unsigned     GetSizeLo() const =0;
    virtual unsigned     GetSizeHi() const =0;
    virtual CStrRef      GetName() const =0;
    virtual CStrRef      GetDescription() const =0;

    virtual int GetOverhangChars(int) const { return 0; }
    
    virtual ~ISts() throw () {}
};
		
END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: sts_i.hpp,v $
 * Revision 1.2  2004/10/26 17:16:34  rotmistr
 * Added 5'-end masking for primers
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
