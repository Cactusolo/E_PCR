/* $Id: seqcmp_main.cpp,v 1.2 2004/06/03 23:37:21 rotmistr Exp $
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

#include <epcr/align.hpp>
#include <stdexcept>
#include <string.h>
#include <errno.h>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);


#ifndef VERSION
#define VERSION "(devel:" __DATE__ ")"
#endif

class CMain
{
public:
	CMain(int c, char ** v):
		argc(c),argv(v),done(false),mism(0),gaps(0),fwd(true),
        aligner(eAlignerLcs) 
        { rev[0]=rev[1]=false; } 
	int Run();
    enum EAligner {
        eAlignerExact,
        eAlignerNoGaps,
        eAlignerFast,
        eAlignerLcs
    }
    aligner;
protected:
	int Execute();
	int ParseCmdline();
	int Help(FILE* = stdout);
	int Version();

protected:
	int argc;
	char ** argv;
	bool done;
protected:
    int mism, gaps;
    bool fwd, rev[2];
};

int CMain::Help(FILE* out)
{
	fprintf(out,"usage: [-hV] [-n mism] [-g gaps] [-s +|-] [-r 1|2] "
            "[-a aligner] "
            "{sequence} {primer}\n");
    fprintf(out,"where:\n"
            "\t-s +|-     - side (begin/end)\n"
            "\t-r 1|2     - flip sequence 1|2\n"
            "\t-a aligner - aligner to use\n"
            "aligeners available are:\n"
            "\texact  - exact match\n"
            "\tnogaps - allow n mismatches\n"
            "\tfast   - allow n mismatches and g gaps,\n"
            "\t         some patterns are incorrecly recognised\n"
            "\tlcs    - use lcs aligner\n");
    
	done=true;
}

struct 
{
    const char * label;
    CMain::EAligner aligner;
}
aligners[]={
    {"exact",CMain::eAlignerExact},
    {"ex",CMain::eAlignerExact},
    {"e",CMain::eAlignerExact},
    {"nogaps",CMain::eAlignerNoGaps},
    {"ngaps",CMain::eAlignerNoGaps},
    {"ng",CMain::eAlignerNoGaps},
    {"n",CMain::eAlignerNoGaps},
    {"fast",CMain::eAlignerFast},
    {"f",CMain::eAlignerFast},
    {"lcs",CMain::eAlignerLcs},
    {"l",CMain::eAlignerLcs},
    {0,CMain::eAlignerLcs}
};


int CMain::ParseCmdline()
{
	int optchar;
	while((optchar=getopt(argc,argv,"hVn:g:s:r:a:"))!=-1) {
		switch(optchar) {
		case 'h': Help(); break;
		case 'V': Version(); break;
		case 'n': mism=atoi(optarg); break;
		case 'g': gaps=atoi(optarg); break;
		case 's': fwd=*optarg=='+'?true:*optarg=='-'?false:fwd; break;
        case 'r': optarg[0]=='1'?rev[0]=true:optarg[0]=='2'?rev[1]=true:false;
            break;
        case 'a': 
            for(int i=0; aligners[i].label; ++i) {
                if(strcmp(optarg,aligners[i].label)==0) {
                    aligner=aligners[i].aligner;
                    goto next;
                }
            }
            throw runtime_error("Unknown aligner: "+string(optarg));
		}
    next:;
	}
	if(done) return 0;
	if(optind >= argc+2) { Help(stderr); return 1; }
	return 0;
}

int CMain::Execute()
{
    IAlign * align;
    
    switch(aligner) {
    case eAlignerExact:  align=new CAlignExact(); break;
    case eAlignerNoGaps: align=new CAlignNoGaps(mism); break;
    case eAlignerFast:   align=new CAlignFast(mism,gaps); break;
    case eAlignerLcs:    align=new CAlignLCS(mism,gaps); break;
    }
    
    char * seq=argv[optind];
    if(rev[0]) seq=FlipSequence(seq);
    char * end=seq+strlen(seq);
    char * primer=argv[optind+1];
    if(rev[1]) primer=FlipSequence(primer);
    int len=strlen(primer);
    bool rc=fwd?
        align->Forward(seq,end,primer,len):
        align->Reverse(seq,end,primer,len);
    printf("seq: %s\npri: %s\n%s [%d mism, %d gaps]\n",
           seq,primer,rc?"OK":"fail",
           align->GetMismatches(),align->GetGaps());

    if(gaps) {
        CLcsMatrix<char> matrix(256,gaps);
        vector<string> graph;
        if(fwd) {
            matrix.Build<const char *>(seq,end,primer,len);
            matrix.Graph<const char *>(seq,end,primer,len,graph);
            matrix.Stat<const char *>(seq,end,primer,len);
        } else {
            matrix.Build<CReverseConstSeqIterator<const char> >(
                end-1,seq-1,primer+len-1,len);
            matrix.Graph<CReverseConstSeqIterator<const char> >(
                end-1,seq-1,primer+len-1,len,graph);
            matrix.Stat<CReverseConstSeqIterator<const char> >(
                end-1,seq-1,primer+len-1,len);
        }
        printf("Primer: %s\n",graph[0].c_str());
        printf("        %s\n",graph[2].c_str());
        printf("Seqnce: %s\n",graph[1].c_str());
    }
    delete align;
    if(rev[0]) delete[] seq;
    if(rev[1]) delete[] primer;
}



int CMain::Version()
{
	done=true;
	puts("Fasta converter for e-PCR version " VERSION);
}

int CMain::Run() 
{
	if(int rc=ParseCmdline() ) return rc;
	if(done) return 0;
	return Execute();
}

int main(int argc, char ** argv) 
{
	try {
		CMain app(argc,argv);
		return app.Run();
	}
	catch(logic_error& e) {
		fprintf(stderr,"! Fatal: Internal error %s\n",e.what());
	}
	catch(exception& e) {
		fprintf(stderr,"! Fatal: %s\n",e.what());
	}
	catch(...) {
		fprintf(stderr,"! Fatal: Unknown error\n");
	}
	return 100;
}
	


/*
 * $Log: seqcmp_main.cpp,v $
 * Revision 1.2  2004/06/03 23:37:21  rotmistr
 * New aligner added.
 *
 * Revision 1.1  2004/02/12 21:38:20  rotmistr
 * Fixed typo in seqcmp
 * Optimized and fixed lookup
 * Better look for reverse.cgi
 *
 * Revision 1.2  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
