// $Id: defaults.h,v 1.2 2004/03/30 18:52:19 rotmistr Exp $
/* ===========================================================================
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
 * ========================================================================= */

#ifndef EPCR_DEFAULTS__H
#define EPCR_DEFAULTS__H

#define ePCR_WDSIZE_DEFAULT   7
#define ePCR_WDSIZE_MIN       3
#define ePCR_WDSIZE_MAX       8
//// Number of mismatches allowed
#define ePCR_MMATCH_DEFAULT   0
#define ePCR_MMATCH_MIN       0
#define ePCR_MMATCH_MAX       10
#define ePCR_GAPS_DEFAULT   0
#define ePCR_GAPS_MIN       0
#define ePCR_GAPS_MAX       5
//// Margin (allowed deviation in product size)
#define ePCR_MARGIN_DEFAULT   50
#define ePCR_MARGIN_MIN       0
#define ePCR_MARGIN_MAX       10000

#define ePCR_DEFAULT_size_lo 100
#define ePCR_DEFAULT_size_hi 350

#endif

/*
 * $Log: defaults.h,v $
 * Revision 1.2  2004/03/30 18:52:19  rotmistr
 * Updated default STS size
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 * Revision 1.2  2003/11/20 02:12:29  rotmistr
 * Fixed id, log tags and copyright notice
 *
 */
