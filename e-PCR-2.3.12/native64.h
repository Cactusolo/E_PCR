/* $Id: native64.h,v 1.3 2004/09/03 19:10:21 rotmistr Exp $
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

#ifndef EPCR_NATIVE64__HPP
#define EPCR_NATIVE64__HPP

#include <stdio.h>

#define O_LARGEFILE 0

typedef off_t off64_t;

#define fopen64  fopen
#define fseeko64 fseeko
#define ftello64 ftello
#define lseek64  lseek
#define mmap64   mmap

#endif

/*
 * $Log: native64.h,v $
 * Revision 1.3  2004/09/03 19:10:21  rotmistr
 * Public domain notice added.
 *
 */
