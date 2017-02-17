/* $Id: getopt.c,v 1.2 2004/09/03 19:59:25 rotmistr Exp $
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

#include <string.h>
#include <stdio.h>

int optind=1;
int optopt=-1;
int opterr=0;
const char* optarg=0;
static char * nextarg=0;

int getopt(int argc, char ** argv, const char* optstring) 
{
    while(optind<argc) {

        if(nextarg==0) {
            if(argv[optind][0]=='-' || argv[optind][0]=='/') {
                nextarg=argv[optind]+1;
            } else {
                // first non-flag
                return optopt=-1;
            }
        } 
        switch(optopt=*nextarg) {
        case '-': ++optind; // -- explicitely stops processing, ignored
        case 0:   return optopt=-1; // - stops processing, passed
        default: 
            do {
                const char * x=strchr(optstring,optopt);
                if(x) {
                    if(x[1]==':') {
                        if(nextarg[1]) {
                            optarg=nextarg+1;
                        }
                        else {
                            if(++optind>=argc) {
                                fprintf(stderr,
                                        "getopt: need argument for %c\n",
                                        optopt);
                                ++opterr;
                                return -1;
                            }
                            optarg=argv[optind];
                        }
                        nextarg=0;
                        ++optind;
                    }
                    else {
                        optarg=0;
                        if(nextarg[1]) ++nextarg;
                        else { nextarg=0; ++optind; }
                    }
                    return optopt;
                } else {
                    fprintf(stderr,
                            "getopt: invalid option %c\n",
                            optopt);
                    if(nextarg[1]) ++nextarg;
                    else { nextarg=0; ++optind; }
                    ++opterr;
                }
            } while(0);
            break;
        case ':': // ':' should not be used as option
            fprintf(stderr,"getopt: bad option %c\n",optopt);
            if(nextarg[1]) ++nextarg;
            else { nextarg=0; ++optind; }
            ++opterr;
            break;
        }
    }
    
    return optopt=-1;
}
    

/*
 * $Log: getopt.c,v $
 * Revision 1.2  2004/09/03 19:59:25  rotmistr
 * *** empty log message ***
 *
 * Revision 1.1  2004/04/02 15:43:55  rotmistr
 * *** empty log message ***
 *
 */
