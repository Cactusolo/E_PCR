/* $Id: stsmatch.hpp,v 1.3 2004/05/27 20:35:48 rotmistr Exp $
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

#ifndef EPCR_STSMATCH__HPP
#define EPCR_STSMATCH__HPP

#include <epcr/build_cfg.h>
#include <epcr/stsmatch_m.hpp>
#include <epcr/defaults.h>


BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class CPcrMachineCompat;

class CStsFileCallbackDefault:public IStsFileCallback
{
public:
    CStsFileCallbackDefault();
    virtual ~CStsFileCallbackDefault() throw () {}
    virtual bool NextLine(const char * , int);
    virtual bool Error(EError);
    virtual bool Start();
    virtual bool Done();
protected:
    const char * m_String;
    unsigned     m_Length;
    int m_Line;
    int m_Bad[3];
};
			

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: stsmatch.hpp,v $
 * Revision 1.3  2004/05/27 20:35:48  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.2  2004/03/29 21:25:40  rotmistr
 * Dist files are prepared
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
