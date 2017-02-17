/* $Id: stsmatch.cpp,v 1.5 2004/06/03 23:37:22 rotmistr Exp $
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

#include <epcr/defaults.h>
#include <epcr/hashset.hpp>
#include <epcr/stsmatch.hpp>

#include <errno.h>
#include <ctype.h>
#include <cstdio>
#include <stdexcept>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

////////////////////////////////////////////////////////////////////////

CStsFileCallbackDefault::CStsFileCallbackDefault() 
{
	Start();
}

bool CStsFileCallbackDefault::Start()
{
	m_String="";
	m_Length=0;
	m_Line=0;
	memset(m_Bad,0,sizeof(m_Bad));
	return true;
}

bool CStsFileCallbackDefault::NextLine(const char * s, int l)
{
	m_String=s;
	m_Length=l;
	++m_Line;
	return true;
}

const char * errorMessage(IStsFileCallback::EError e)
{
	switch(e){
	case IStsFileCallback::eErrOK:return "OK";
	case IStsFileCallback::eErrShortPrimer:return "short primer";
	case IStsFileCallback::eErrAmbiquosPrimer:return "ambiquities in primer";
	case IStsFileCallback::eErrBadLine:return "bad line";
	default: return "unknown error";
	}
}

bool CStsFileCallbackDefault::Error(EError err)
{
	const char * msg="OK";
	
	switch(err) {
    case eErrTOTAL: throw logic_error("err can't be eErrTOTAL here");
	case eErrOK: return true;
	case eErrBadLine:        m_Bad[0]++; msg="bad line"; break;
	case eErrShortPrimer:    m_Bad[1]++; msg="short primer"; break;
	case eErrAmbiquosPrimer: m_Bad[2]++; msg="ambiquity in primer"; break;
	case eErrSystem: 
		throw runtime_error(string("STS file error: ")+strerror(errno));
	}
	
	return true;
}

bool CStsFileCallbackDefault::Done()
{
	if(m_Bad[0]) 
		fprintf(stderr,"WARNING: %d STSs have incomplete description line\n",
				m_Bad[0]);
	if(m_Bad[1]) 
		fprintf(stderr,"WARNING: %d STSs have primer shorter than W\n",
				m_Bad[1]);
	if(m_Bad[2]) 
		fprintf(stderr,
				"WARNING: %d STSs have ambiguities within W of 3\' end\n",
				m_Bad[2]);
	return true;
}

////////////////////////////////////////////////////////////////////////

/*
 * $Log: stsmatch.cpp,v $
 * Revision 1.5  2004/06/03 23:37:22  rotmistr
 * New aligner added.
 *
 * Revision 1.4  2004/05/27 20:35:48  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.3  2004/04/01 05:57:53  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.2  2004/03/29 21:25:40  rotmistr
 * Dist files are prepared
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
