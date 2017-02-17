/* $Id: faread.cpp,v 1.6 2004/04/06 04:53:17 rotmistr Exp $
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

#include <epcr/faread.hpp>
#include <epcr/bin-io.hpp>

#include <stdexcept>
#include <string.h>
#include <ctype.h>
#include <errno.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

struct SFile
{
	bool m_pipe;
	FILE * f;
	SFile(bool pipe, const char * name) {
#ifdef NO_POPEN
            if(pipe) 
				throw runtime_error("Pipes are not supported in windows!");
            f=fopen(name,"r");
#else
            f=pipe?popen(name,"r"):fopen64(name,"r");
#endif
            if(f==0)  throw runtime_error(string(name)+": "+strerror(errno));
		setvbuf(f,0,_IOFBF,16192);
		m_pipe=pipe;
	}
        ~SFile() {
#ifndef NO_POPEN
            if(m_pipe) pclose(f); else
#endif
                fclose(f);
        }
};

void CFastaReader::Open(const string& fname) 
{
    if(IsOpen()) Close();
    
//     FILE * f=fopen64(fname.c_str(),"r");
//     if(f==0) throw runtime_error(fname+": "+strerror(errno));
//     setvbuf(f,0,_IOFBF,16192);
    m_Fptr=new SFile(false,fname.c_str());
}

void CFastaReader::PipeIn(const string& command)
{
	if(IsOpen()) Close();
	m_Fptr=new SFile(true,command.c_str());
}

void CFastaReader::Close() 
{
    if(m_Fptr) {
//        fclose((FILE*)m_Fptr);
		delete (SFile*)m_Fptr;
        m_Fptr=0;
    }
}

void CFastaReader::ReadFile(IFastaReaderCallback * cbk)
{
    if(!cbk) return;
    if(!IsOpen()) return;
    
    FILE * f=((SFile*)m_Fptr)->f;

    char buffer[16192];
    bool is_defline=false;
    bool is_newline=true;
    bool beginning_of_file=true;
    
    cbk->CbkFileBegin();

    string defline;
    while(!feof(f)) {
        if(!fgets(buffer,sizeof(buffer),f)) break;

        char * eol=strlen(buffer)+buffer;
        bool complete_line=eol>buffer && (eol[-1]=='\r' || eol[-1]=='\n');
        while(eol>buffer && isspace(eol[-1])) --eol;
        *eol=0;

        if(!is_newline) {
            if(is_defline) { defline.append(buffer); }
            else SeqlineCbk(cbk,buffer,eol-buffer);
        }
        else {
            if(*buffer=='>') {
                if(!beginning_of_file)  
                    cbk->CbkEntryEnd();
                cbk->CbkEntryBegin();

                beginning_of_file=false;

                is_defline=true;
                defline.assign(buffer);
            }
            else {
                if(beginning_of_file)
                    throw runtime_error("Leading garbage in fasta file!");
                
                if(defline.length()) { // should have at least '>'
                    cbk->CbkDefline(defline.c_str(),defline.length());
                    const char * s=defline.c_str()+1;
                    for(;*s &&  isspace(*s);++s);
                    const char * ss=s;
                    for(;*s && !isspace(*s);++s);
                    cbk->CbkIdent(ss,s-ss);
                    defline.clear();
                }
                
                is_defline=false;
                SeqlineCbk(cbk,buffer,eol-buffer);
            }
        }
        is_newline=complete_line;
    }
    if(!beginning_of_file) cbk->CbkEntryEnd();
    cbk->CbkFileEnd();
}

void CFastaReader::SeqlineCbk(IFastaReaderCallback* cbk, 
                              char * buffer, unsigned length) 
{
    char * eol=buffer+length;
    
    if(m_CvtTable) {
        
        for(char * x=buffer; x<eol; ++x) {
            char c=m_CvtTable[*x];
            if(c==0) {
                char * xx=x+1;
                while(xx<eol && m_CvtTable[*xx]==0) ++xx;
                memmove(x,xx,eol-xx);
                eol -= xx-x;
                if(x>=eol) break;
                c=m_CvtTable[*x];
            }
            *x=c;
        }
    }
    cbk->CbkSeqline(buffer,eol-buffer);
}


char CFastaReader::sm_NucleotidesExt[]=
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0ABCDNNGHNNKNMNNNNRNTTNWNYN\0\0\0\0\0"
//.abcdefghijklmnopqrstuvwxyz.....
"\0abcdnnghnnknmnnnnrnttnwnyn\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
;
char CFastaReader::sm_NucleotidesExtUc[]=
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0ABCDNNGHNNKNMNNNNRNTTNWNYN\0\0\0\0\0"
//.abcdefghijklmnopqrstuvwxyz.....
"\0ABCDNNGHNNKNMNNNNRNTTNWNYN\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
;
char CFastaReader::sm_Nucleotides[]=
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0ANCNNNGNNNNNNNNNNNNTNNNNNN\0\0\0\0\0"
//.abcdefghijklmnopqrstuvwxyz.....
"\0ancnnngnnnnnnnnnnnntnnnnnn\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
;
char CFastaReader::sm_NucleotidesUc[]=
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0ANCNNNGNNNNNNNNNNNNTNNNNNN\0\0\0\0\0"
//.abcdefghijklmnopqrstuvwxyz.....
"\0ANCNNNGNNNNNNNNNNNNTNNNNNN\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
;

/*
 * $Log: faread.cpp,v $
 * Revision 1.6  2004/04/06 04:53:17  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.5  2004/04/01 05:57:52  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.4  2004/03/23 22:35:25  rotmistr
 * Fixed processing of -mid flag in cmdline
 * Fixed destructor for fasta reader
 * Removed cgi
 *
 * Revision 1.3  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.2  2004/01/08 23:22:41  rotmistr
 * Fixed init error in faread,
 * Adjusted output to standard,
 * Added output format style and output file to parameters.
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
