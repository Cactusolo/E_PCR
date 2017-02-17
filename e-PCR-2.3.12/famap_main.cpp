/* $Id: famap_main.cpp,v 1.7 2008/04/28 16:38:45 rotmistr Exp $
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

#include <epcr/fast_seqio.hpp>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cstring>

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
		argc(c),argv(v),done(false),command(eStat),cvt(0) {} 
	int Run();
protected:
	int Execute();
	int ParseCmdline();
	int Help(FILE* = stdout);
	int Version();

	int Build();
	int Stat();
	int Dump();
	int List();
	
protected:
	int argc;
	char ** argv;
	bool done;
protected:
	string famap;
	enum ECommand { eBuild, eStat, eDump, eList } command;
	const char * cvt;
};

int CMain::Help(FILE* out)
{
	fprintf(out,"usage: [-hV] -b mmapped-file [-t cvt] [fafile ...]\n");
	fprintf(out,"   or: [-hV] -d mmapped-file [ord ...]\n");
	fprintf(out,"   or: [-hV] -l mmapped-file [ord ...]\n");
	fprintf(out,"where cvt (convertion table) is one of:\n"
			"\toff - as is (default)\n"
			"\tn   - nucleotide [acgtnACGTN] allowed,\n"
			"\tN   - nucleotide uppercase allowed [ACGTN]\n"
			"\tnx  - nucleotide with ambiquity codes allowed\n"
			"\tNX  - nucleotide with ambiquity codes uppercase\n");
	done=true;
	return 0;
}

int CMain::ParseCmdline()
{
	int optchar;
	while((optchar=getopt(argc,argv,"hVb:d:l:t:"))!=-1) {
		switch(optchar) {
		case 'h': Help(); break;
		case 'V': Version(); break;
		case 'b': command=eBuild; famap=optarg; break;
		case 'd': command=eDump; famap=optarg; break;
		case 'l': command=eList; famap=optarg; break;
		case 't':
			if(strcmp(optarg,"off")==0||strcmp(optarg,"-")==0) cvt=0;
			else if(strcmp(optarg,"n")==0) cvt=CFastaReader::sm_Nucleotides;
			else if(strcmp(optarg,"N")==0) cvt=CFastaReader::sm_NucleotidesUc;
			else if(strcmp(optarg,"nx")==0) 
				cvt=CFastaReader::sm_NucleotidesExt;
			else if(strcmp(optarg,"NX")==0) 
				cvt=CFastaReader::sm_NucleotidesExtUc;
			else fprintf(stderr,"Unknown table %s\n",optarg);
			break;
		}
	}
	if(done) return 0;
	if(command == eBuild && optind >= argc) { Help(stderr); return 1; }
	return 0;
}

int CMain::Build()
{
	do {
		
		CFastaMapPrepare prepare(famap+"~");
		
		for(int i=optind; i<argc; ++i) {
			fprintf(stderr,"%s* Adding %s\n",
					isatty(fileno(stderr))?"\r\x1b[K":"",
					argv[i]);
			prepare.AddFile(argv[i],cvt);
		}
	} while(0);
	
	if(rename((famap+"~").c_str(),famap.c_str())) {
		throw runtime_error("Rename failed for "+famap+": "+strerror(errno));
	}
	
	return 0;
}


int CMain::Execute()
{
	switch(command) {
	case eBuild: return Build(); 
	case eStat: return Stat();
	case eDump: return Dump();
	case eList: return List();
	}
}

int CMain::Stat()
{
	CFastaMap fmap(famap);
	fprintf(stdout,"Sequence count: %d",fmap.SequenceCount());
	
	return 0;
}

int CMain::Dump()
{
	CFastaMap fmap(famap);
    if(optind >= argc) {
        
        for(unsigned i=0; i<fmap.SequenceCount(); ++i) {
            printf("%s\n",fmap.GetDefline(i).c_str());
            CMmSequence seq(fmap,i);
            for(unsigned j=0; j<seq.GetSize(); j+=75) {
                fwrite(seq.GetData()+j,1,min(75U,seq.GetSize()-j),stdout);
                putc('\n',stdout);
            }
        }
        
	} else {
        
        for(unsigned a=optind; a<argc; ++a) {
            int x=atoi(argv[a]);
            if(x>=fmap.SequenceCount()) 
                fprintf(stderr,"? Id %d is out of range\n",x);
            else {
                printf("%s\n",fmap.GetDefline(x).c_str());
                CMmSequence seq(fmap,x);
                for(unsigned j=0; j<seq.GetSize(); j+=75) {
                    fwrite(seq.GetData()+j,1,min(75U,seq.GetSize()-j),stdout);
                    putc('\n',stdout);
                }
            }
        }
    }
    

    return 0;
}

int CMain::List()
{
	CFastaMap fmap(famap);
    if(optind >= argc) {
        
        for(unsigned i=0; i<fmap.SequenceCount(); ++i) {
            printf("%d\t%s\t%d\n",i,fmap.GetIdent(i).c_str(),fmap.GetSize(i));
//             CMmSequence seq(fmap,i);
//             for(unsigned j=0; j<seq.GetSize(); j+=75) {
//                 fwrite(seq.GetData()+j,1,min(75U,seq.GetSize()-j),stdout);
//                 putc('\n',stdout);
//             }
        }
        
	} else {
        
        for(unsigned a=optind; a<argc; ++a) {
            int x=atoi(argv[a]);
            if(x>=fmap.SequenceCount()) 
                fprintf(stderr,"? Id %d is out of range\n",x);
            else {
                printf("%d\t%s\t%d\n",x,fmap.GetIdent(x).c_str(),
                       fmap.GetSize(x));
            }
        }
    }
    

    return 0;
}

int CMain::Version()
{
	done=true;
	puts("Fasta converter for e-PCR version " VERSION);
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
 * $Log: famap_main.cpp,v $
 * Revision 1.7  2008/04/28 16:38:45  rotmistr
 * Applied patch to build with gcc-4.3
 *
 * Revision 1.6  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.5  2004/04/06 04:53:17  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.4  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.3  2004/02/12 21:38:20  rotmistr
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
