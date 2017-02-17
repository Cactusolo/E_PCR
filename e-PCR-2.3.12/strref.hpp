/* $Id: strref.hpp,v 1.3 2008/04/28 16:38:45 rotmistr Exp $
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

#ifndef EPCR_STRREF__HPP
#define EPCR_STRREF__HPP

#include <epcr/build_cfg.h>
#include <cstring>
#include <string>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

	class CStrRef;

// This class has string interface, but does NOT create copy of data
// Not safe, but VERY useful
class CStrRef
{
public:
	CStrRef():m_Data(0),m_Length(0) {}
	CStrRef(const char * s):m_Data(s),m_Length(strlen(s)) {}
	CStrRef(const char * s, unsigned l):m_Data(s),m_Length(s?l:0) {}
	CStrRef(const string& s):m_Data(s.c_str()),m_Length(s.length()) {}
	unsigned     length() const { return m_Length; }
    void resize(unsigned u) { m_Length=u; }
	const char * data() const { return m_Data; }
	operator string () const { return string(data(),length()); }
protected:
	const char * m_Data;
	unsigned     m_Length;
};

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: strref.hpp,v $
 * Revision 1.3  2008/04/28 16:38:45  rotmistr
 * Applied patch to build with gcc-4.3
 *
 * Revision 1.2  2004/10/26 14:25:29  rotmistr
 * Added resize()
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
