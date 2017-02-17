/* $Id: fahash_internal.hpp,v 1.3 2007/07/11 20:49:29 rotmistr Exp $
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

#ifndef FAHASH_INTERNAL__HPP
#define FAHASH_INTERNAL__HPP

#define SYSERROR(a) throw runtime_error(a+": "+strerror(errno))

typedef AFaIndexerBase::THashElement THashElement;
static const THashElement kHighBit = THashElement( 1 ) << ( 8 * sizeof(THashElement) - 1 );

static const unsigned kKilobyte = 1024;
static const unsigned kMegabyte = 1024 * kKilobyte;
static const unsigned kGigabyte = 1024 * kMegabyte;

static const unsigned kMinFragSize = kGigabyte / 4 / 4 * 3;
static const unsigned kMaxFragSize = kGigabyte / 4 / 4 * 6;
//static const unsigned min_frag_size=100*Megabyte/4*3;
//static const unsigned max_frag_size=100*Megabyte/4*6;

#endif

/*
 * $Log: fahash_internal.hpp,v $
 * Revision 1.3  2007/07/11 20:49:29  rotmistr
 * Made 64bit-compatible
 *
 * Revision 1.2  2004/09/03 19:06:41  rotmistr
 * Code formatting changes
 *
 */
