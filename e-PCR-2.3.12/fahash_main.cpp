/* $Id: fahash_main.cpp,v 1.5 2008/04/28 16:38:45 rotmistr Exp $
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
#include <epcr/stsmatch.hpp>
#include <epcr/fahash.hpp>
#include <epcr/sts.hpp>

//#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include <cstdio>
#include <memory>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

class CMain
{
public:
	CMain(int c, char ** v):
		argc(c),argv(v),done(false),
        wdsize(8),period(0),
        lo(0),hi(0),skip_repetative(false),cachesize(0),
        command(eNone),version(eVer2) {} 
	int Run();
protected:
	int Execute();
	int ParseCmdline();
	int Help(FILE* = stdout);
	int Version();

	int Build();
	int Stat();
	
protected:
	int argc;
	char ** argv;
	bool done;
protected:
	string findex;
	int    wdsize, period;
    int lo, hi;
	bool skip_repetative;
    unsigned cachesize;
	enum ECommand { eNone, eBuild, eStat } command;
    enum EVersion { eVer1=1, eVer2=2 } version;
    string outfile;
};

int CMain::Help(FILE* out)
{
	fprintf(out,"usage: [-hV] -b hash-file [-w wdsize] [-f period] "
            "[-F fragment_min,fragment_max] [-k] [-c cachesize] [-v 1|2] "
			"famap-file ...\n");
	fprintf(out,"   or: [-hV] -T hash-file [-o outfile]\n");
    fprintf(out,"where:\n"
            "\t-T hash-file\tPrint word usage statistics for hash-file\n"
            "\t-b hash-file\tBuild hash tables (hash-file) "
            "from sequence files,\n"
            "\t-w wordsize \tSet word size when building hash tables\n"
            "\t-f period   \tSet discontiguity when building hash tables\n"
            "\t-k          \tSkip repeats when building hash-file\n"
            "\t-F min,max  \tSet watermarks for fragment size (in Mb) "
            "(version 1 only)\n"
            "\t-c cachesize\tSet cache size (version 2 only)\n"
            "\t-v ver      \tUse format version (1|2, 2 is default)\n"
            "\t-o outfile  \tWrite output to file `outfile'\n");
            
	done=true;
	return 0;
}

int CMain::ParseCmdline()
{
	int optchar;
    while((optchar=getopt(argc,argv,"hVw:F:f:b:kT:v:o:c:"))!=-1) {
        switch(optchar) {
        case 'h': Help(); break;
        case 'V': Version(); break;
        case 'w': wdsize=atoi(optarg); break;
        case 'f': period=atoi(optarg); break;
        case 'c': cachesize=atoi(optarg); break;
        case 'k': skip_repetative=true; break;
        case 'b': findex=optarg; command=eBuild; break;
        case 'T': findex=optarg; command=eStat; break;
        case 'o': outfile=optarg; break;
		case 'F': 
            do { 
                char * x=0; 
                lo=strtol(optarg,&x,10);  
                if(x==0 || strchr(":,-;/",*x)==0) 
                    throw runtime_error("fragment size should be RANGE");
                hi=strtol(x+1,0,10);
            } while(0);
            break;
        case 'v': 
            do {
                int v=atoi(optarg);
                switch(v) {
                case 1: version=eVer1; break;
                case 2: version=eVer2; break;
                default: 
                    throw runtime_error("only file versions 1 and 2 "
                                        "are supported");
                }
            } while(0);
            break;
		}
	}
	if(done) return 0;
	if(command !=eStat &&
       (command==eNone || optind >= argc)) { Help(stderr); return 1; }
	return 0;
}

#ifndef USE_WIN
#define CLREOL "\r\x1b[K"
#else
#define CLREOL "\r"
#endif

class CFaIndexerCallback:public IFaIndexerCallback
{
    time_t t0;
public:
    CFaIndexerCallback():t0(time(0)) {}
    
    virtual void CbkSequence(const char * id) { 
        fprintf(stderr,CLREOL " - Adding Sequence %s\n",id);
        fflush(stderr);
    }
    virtual void CbkProgress(unsigned pos, unsigned length) {
		time_t t1=time(0);
        if(t0==t1) return;
		if( length > 10000 ) length/=100;
		else pos*=100;
		if(length) {
	        fprintf(stderr,CLREOL "        %50s\r   %3d%% %.*s",
    	            "..................................................",
        	        pos/length,pos/length/2,
            	    "##################################################");
        	fflush(stderr);
		}
        t0=t1;
    }
    virtual void CbkFile(const char * name) {
        fprintf(stderr,CLREOL "* Adding File %s\n",name);
        fflush(stderr);
    }
    virtual void CbkDumpProgress(unsigned pos, unsigned length) {
		time_t t1=time(0);
        if(t0==t1) return;
		if( length > 10000 ) length/=100;
		else pos*=100;
		if(length) {
	        fprintf(stderr,CLREOL "        %50s\r   %3d%% %.*s",
    	            "oooooooooooooooooooooooooooooooooooooooooooooooooo",
        	        pos/length,pos/length/2,
            	    "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO");
	        fflush(stderr);
		}
        t0=t1;
    }
    virtual void CbkDumpStart() { 
        fprintf(stderr,CLREOL "= Dumping \n"); 
        fflush(stderr);
    }
    virtual void CbkDumpEnd() { 
        fprintf(stderr,CLREOL "= Dumping OK\n"); 
        fflush(stderr);
    }
    virtual void CbkResetStart() {
        fprintf(stderr,CLREOL "= Resetting..."); 
        fflush(stderr);
    };
    virtual void CbkResetEnd() {
        fprintf(stderr,"\r\r\r OK"); 
    };
};

int CMain::Build()
{
    auto_ptr<AFaIndexerBase> indexer(
        version==eVer1?new CFaIndexer1():
        version==eVer2?new CFaIndexer2():
        (AFaIndexerBase*)0);

    switch(version) {
    case eVer1: 
        if(CFaIndexer1* ind=dynamic_cast<CFaIndexer1*>(indexer.get())) {
            if(lo && hi) 
                ind->SetFragmentSizeRange(unsigned(lo)*1024U*1024U,
                                          unsigned(hi)*1024U*1024U);
            if(skip_repetative) 
                ind->SetFlags(ind->GetFlags()&CFaIndexer1::fSkipRepeatitive);
        }
        break;
        
    case eVer2: 
        if(CFaIndexer2* ind=dynamic_cast<CFaIndexer2*>(indexer.get())) {
            if(cachesize)
                ind->SetCacheSize(cachesize);
        }
        break;

    default:
        throw logic_error("Unknown index version");
    }

    CFaIndexerCallback cbk;
    indexer->SetCallback(&cbk);
	indexer->SetHash(CHashSet(wdsize,period));
	indexer->AttachFile(findex+"~");

	for(int i=optind; i<argc; ++i) {
		indexer->AddFile(argv[i]);
	}
	indexer->Finish();
	if(rename((findex+"~").c_str(),findex.c_str())) {
		throw runtime_error("Rename failed for "+findex+": "+strerror(errno));
	}
	return 0;
}


int CMain::Execute()
{
	switch(command) {
	case eBuild: return Build(); 
	case eStat: return Stat();
	}
}

int CMain::Stat()
{
	CFaLookup lookup;
	lookup.AttachFile(findex);
    fputs("* Calculating statistics...\n",stderr);
    lookup.Stat();
    fputs("* Done\n",stderr);
    FILE * out=outfile.length()?fopen64(outfile.c_str(),"w"):stdout;
    if(out==0) 
        throw runtime_error("Failed to open "+outfile+": "+strerror(errno));
    
    const CFaLookup::TStatVector& st=lookup.GetStat();
    int count=0;
    double average=0;
    
    for(unsigned wd=0; wd<st.size(); ++wd) {
        const CFaLookup::TBranchStat& x=st[wd];
        for(unsigned i=0; i!=x.size(); ++i, ++count) {
            average+=x[i];
        }
    }
    average/=count;
    
    for(unsigned wd=0; wd<st.size(); ++wd) {
        const CFaLookup::TBranchStat& x=st[wd];
        for(unsigned i=0; i!=x.size(); ++i) {
            char wtype[128];
            memset(wtype,0,sizeof(wtype));
            char * ptr=wtype+lookup.GetHash().GetWordSize();
            if(lookup.GetHash().GetWordCount()==1) {
                for(int k=0; k<lookup.GetHash().GetWordSize(); ++k) {
                    *--ptr="ACGT"[(i>>(k*2)) & 3];
                }
            } else {
                for(int k=0; ptr>wtype; ++k) {
                    if(((ptr-wtype+wd)%lookup.GetHash().GetWordCount())==0) 
                        *--ptr='*';
                    *--ptr="ACGT"[(i>>(k*2)) & 3];
                }
            }
            fprintf(out,"%d\t%d\t%s\t%lld\t%lf\n",wd,i,wtype,
                    x[i],x[i]/average);
        }
    }
	if(outfile.length()) fclose(out);
    
	return 0;
}

int CMain::Version()
{
	done=true;
	puts("Reverse e-PCR: sequence hash builder version " VERSION);
	return 0;
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
 * $Log: fahash_main.cpp,v $
 * Revision 1.5  2008/04/28 16:38:45  rotmistr
 * Applied patch to build with gcc-4.3
 *
 * Revision 1.4  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.3  2004/09/03 21:28:49  rotmistr
 * Fixes to compile with Borland C++ 5.5
 *
 * Revision 1.2  2004/06/03 23:37:20  rotmistr
 * New aligner added.
 *
 * Revision 1.1  2004/05/27 20:35:47  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 *
 */
