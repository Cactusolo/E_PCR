/* $Id: re-PCR_main.cpp,v 1.21 2008/04/28 16:38:45 rotmistr Exp $
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
#include <epcr/align.hpp>
#include <epcr/sts.hpp>

//#include <unistd.h>
#include <errno.h>

#include <stdexcept>
#include <cstdio>
#include <string>
#include <memory>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

class CMain
{
public:
	CMain(int c, char ** v):
		argc(c),argv(v),done(false),
        margin(50),gaps(0),mism(0),rev_lookup(true),presize(false),
        command(eNone),mode(eCmdline),optimize(false),batchcnt(1000),
        defLo(ePCR_DEFAULT_size_lo),defHi(ePCR_DEFAULT_size_hi),
        m_PrintAlignments(false),m_Quiet(false) {} 
	int Run();
protected:
	int Execute();
	int ParseCmdline();
	int Help(FILE* = stdout);
	int Version();

	int PrimerLookup();
	int STSLookup();
	
protected:
	int argc;
	char ** argv;
	bool done;
protected:
	string findex;
	int margin;
    int gaps, mism;
    bool rev_lookup;
    bool presize;
	enum ECommand { eNone, ePrimerLookup, eSTSLookup } command;
    enum EMode { eCmdline, eFile } mode;
    bool optimize;
    string outfile;
    int batchcnt;
	int defLo, defHi;
    IAlign * m_AlignL, * m_AlignR;
    bool m_PrintAlignments;
    bool m_Quiet;
};

int CMain::Help(FILE* out)
{
	fprintf(out,"usage: [-hV] -p hash-file [-g gaps] [-n mism] [-lq] "
            "[primer ...]\n");
	fprintf(out,"   or: [-hV] -P hash-file [-g gaps] [-n mism] [-l] "
	    "[-m margin] [-O+|-] [-C batchcnt] [-o outfile] [-r+|-] "
	    "[primers-file ...]\n");
	fprintf(out,"   or: [-hV] -s hash-file [-g gaps] [-n mism] [-lq] "
            "[-m margin] [-o outfile] [-r+|-] "
            "[left right lo[-hi] [...]]\n");
	fprintf(out,"   or: [-hV] -S hash-file [-g gaps] [-n mism] [-lq] "
            "[-m margin] [-O+|-] [-C batchcnt] [-o outfile] [-r+|-] "
            "[stsfile ...]\n");
    fprintf(out,"where:\n"
            "\t-p hash-file\tPerform primer lookup using hash-file\n"
            "\t-P hash-file\tPerform primer lookup using hash-file\n"
            "\t-s hash-file\tPerform STS lookup using hash-file\n"
            "\t-S hash-file\tPerform STS lookup using hash-file\n"
            "\t-n mism      \tSet max allowed mismatches per primer "
            "for lookup\n"
            "\t-g gaps      \tSet max allowed indels per primer for lookup\n"
            "\t-m margin    \tSet variability for STS size for lookup\n"
            "\t-l           \tUse presize alignments (only if gaps>0)\n"
            "\t-G           \tPrint alignments in comments\n"
			"\t-d min-max   \tSet default STS size\n"
            "\t-r +|-       \tEnable/disable reverse STS lookup\n"
            "\t-O +|-       \tEnable/disable syscall optimisation\n"
            "\t-C batchcnt  \tSet number of STSes per batch\n"
            "\t-o outfile   \tSet output file name\n"
	    "\t-q           \tQuiet (no progress indicator)\n"
            "Use famap and fahash to generate hash files\n");
    
	done=true;

	return 0;
}

int CMain::ParseCmdline()
{
	int optchar;
	while((optchar=getopt(argc,argv,"hVp:P:s:S:d:m:r:g:n:O:C:o:lGq"))!=-1) {
		switch(optchar) {
		case 'h': Help(); break;
		case 'V': Version(); break;
		case 'p': findex=optarg; command=ePrimerLookup; break;
		case 's': findex=optarg; command=eSTSLookup; break;
		case 'm': margin=strtol(optarg,0,10); break;
        case 'r': rev_lookup=*optarg=='+'?true:*optarg=='-'?false:rev_lookup;
            break;
        case 'g': gaps=strtol(optarg,0,10); break;
        case 'o': outfile=optarg; break;
        case 'n': mism=strtol(optarg,0,10); break;
        case 'l': presize=true; break;
        case 'q': m_Quiet=true; break;
        case 'G': m_PrintAlignments=true; break;
		case 'd': 
			do { 
				char * x=const_cast<char*>(optarg); 
				int lo=strtol(x,&x,10); 
				int hi=(x && *x=='-')?strtol(x+1,&x,10):lo;
				if(lo>0 && lo<hi) {
					defLo=lo;
					defHi=hi;
				}
				else throw runtime_error("bad range: "+string(optarg));
			} while(0);
			break;
		case 'S': findex=optarg; mode=eFile; command=eSTSLookup; break;
		case 'P': findex=optarg; mode=eFile; command=ePrimerLookup; break;
        case 'O': optimize=*optarg=='+'?true:*optarg=='-'?false:optimize; 
            break;
        case 'C': batchcnt=strtol(optarg,0,10); break;
		}
	}
	if(done) return 0;
    if(command==eNone || optind >= argc) { Help(stderr); return 1; }
    if(!gaps) presize=false;
	return 0;
}

#ifndef USE_WIN
#define CLREOL "\r\x1b[K"
#else
#define CLREOL "\r"
#endif

class CFaLookupCallbackBase:public IFaLookupCallback
{
protected:
    FILE * out;
    bool quiet;
	int old;
public:
    CFaLookupCallbackBase(FILE * o=stdout, bool q=false):out(o),quiet(q),old(-1) {}
    virtual bool Fail(const std::string& msg) { 
        fprintf(out,"#- Error: %s\n",msg.c_str()); 
		return true;
    }
    virtual bool Warn(const std::string& msg, const ISts *) { 
        fprintf(out,"#- Warning: %s\n",msg.c_str()); 
		return true;
    }
    virtual bool Done(){ 
        if((!quiet) && isatty(fileno(stderr))) {
            fprintf(stderr,CLREOL "");
            fflush(out);
        }
        fprintf(out,"#- Done\n"); 
		return true;
    }
    virtual void Progress(unsigned i, unsigned total) {
        if((!quiet) && isatty(fileno(stderr))) {
            int I=(+i*50)/total;
			if( I == old ) return;
			old = I;
            fprintf(stderr,CLREOL "- Progress: %3d%% %.*sO%.*s\r",
                    int(0.5+100*(i+0.5)/total),I,
                    "=================================================="
                    "==================================================",
                    50-I-1,
                    "--------------------------------------------------"
                    "--------------------------------------------------");
            fflush(out);
        }
    };
//     virtual void Fragment(unsigned i, unsigned total) {
//         if(isatty(fileno(stderr))) {
//             fprintf(stderr,CLREOL "- Fragment %5d of %-5d %.*sO%.*s\r",
//                     i+1,total,i,
//                     "=================================================="
//                     "==================================================",
//                     total-i-1,
//                     "--------------------------------------------------"
//                     "--------------------------------------------------");
//             fflush(out);
//         }
//     };
};
				
class CFaLookupPrimerCallback:public CFaLookupCallbackBase
{
public:
    CFaLookupPrimerCallback(FILE * o=stdout, bool q=false):CFaLookupCallbackBase(o,q) {}
    virtual bool Start() { 
        fprintf(out,"#- sts\tseq\tstrand\tfrom\tto\tmism\tgaps\n");
		return true;
    }
    
    virtual bool Match(const SFaMatchBlock * info, const ISts *) {
        return Match(info);
		return true;
    }
    virtual bool Match(const SFaMatchBlock * info) {
        fprintf(out,"%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d/%d\n",
                info->sts_label.c_str(),
                info->seq_label.c_str(),
                char(info->strand),
                info->from+1,info->to,
                info->mism,info->gaps,
                info->to-info->from,0);
        return true;
    }
};


class CFaLookupStsCallback:public CFaLookupCallbackBase
{
public:
    CFaLookupStsCallback(FILE * o=stdout, bool pa=false, bool q=false):
        CFaLookupCallbackBase(o,q),m_PrintAlignments(pa),m_Matrix(127,2) {}
    virtual bool Start() { 
        fprintf(out,"#- sts\tseq\tstrand\tfrom\tto\tmism\tgaps\t"
                "act_len/exp_len\n");
		return true;
    }
    virtual bool Match(const SFaMatchBlock * info) {
        return Match(info,0);
    }
    virtual bool Match(const SFaMatchBlock * info, const ISts * sts) {
        fprintf(out,"%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d/%d-%d\n",
                info->sts_label.c_str(),
                info->seq_label.c_str(),
                char(info->strand),
                info->from+1,info->to,
                info->mism,info->gaps,
                info->to-info->from,
                sts?sts->GetSizeLo():0,
                sts?sts->GetSizeHi():0);
        if(m_PrintAlignments && info->sequence && info->seqlen>=info->to) {
            vector<string> left, right;
            m_Matrix.Build<const char*>(
                info->sequence+info->to-sts->GetPrimerLength(ISts::eRight),
                info->sequence+info->seqlen,
                sts->GetPrimerData(ISts::eRight),
                sts->GetPrimerLength(ISts::eRight));
            m_Matrix.Graph<const char*>(
                info->sequence+info->to-sts->GetPrimerLength(ISts::eRight),
                info->sequence+info->seqlen,
                sts->GetPrimerData(ISts::eRight),
                sts->GetPrimerLength(ISts::eRight),
                right);
            m_Matrix.Build<CReverseConstSeqIterator<const char> >(
                info->sequence+info->from+sts->GetPrimerLength(ISts::eLeft)-1,
                info->sequence,
                sts->GetPrimerData(ISts::eLeft)+
                sts->GetPrimerLength(ISts::eLeft)-1,
                sts->GetPrimerLength(ISts::eLeft));
            m_Matrix.Graph<CReverseConstSeqIterator<const char> >(
                info->sequence+info->from+sts->GetPrimerLength(ISts::eLeft)-1,
                info->sequence,
                sts->GetPrimerData(ISts::eLeft)+
                sts->GetPrimerLength(ISts::eLeft)-1,
                sts->GetPrimerLength(ISts::eLeft),
                left);
            int l=max(int(info->seq_label.length()),int(sts->GetName().length()));
            int d=info->to-info->from-
                sts->GetPrimerLength(ISts::eLeft)-
                sts->GetPrimerLength(ISts::eRight);
            
            string stsname(sts->GetName().data(),sts->GetName().length());
            fprintf(out,
                    "#   STS %*s %s...%d...%s\n"
                    "#       %.*s %s   %d   %s\n"
                    "#   Seq %*s %s...%d...%s\n"
                    "#############################################"
                    "#############################\n",
                    l,stsname.c_str(),
                    left[0].c_str(),d,right[0].c_str(),
                    l,"                                                   ",
                    left[2].c_str(),d,right[2].c_str(),
                    l,info->seq_label.c_str(),
                    left[1].c_str(),d,right[1].c_str());
            
        }
        
        return true;
    }
protected:
    bool m_PrintAlignments;
    CLcsMatrix<char> m_Matrix;
};

int CMain::Execute()
{
    if(presize) {
        m_AlignL = new CAlignLCS(mism,gaps);
        m_AlignR = new CAlignLCS(mism,gaps);
    } else if(gaps) {
        m_AlignL = new CAlignFast(mism,gaps);
        m_AlignR = new CAlignFast(mism,gaps);
    } else if(mism) {
        m_AlignL = new CAlignNoGaps(mism);
        m_AlignR = new CAlignNoGaps(mism);
    } else {
        m_AlignL = new CAlignExact();
        m_AlignR = new CAlignExact();
    }
    
	switch(command) {
	case ePrimerLookup: return PrimerLookup();
	case eSTSLookup: return STSLookup();
	}
}

int CMain::PrimerLookup()
{
	CFaLookup lookup;
    lookup.SetAligner(m_AlignL, m_AlignR);
	lookup.AttachFile(findex);
	CFaLookupPrimerCallback cbk(stdout,m_Quiet);
	
    if(mode == eCmdline) {
        
        for(int i=optind; i<argc; ++i) {
            char buffer[1024];
            snprintf(buffer,sizeof(buffer),"PRIMER-%d",(i-optind+1));
            lookup.Find(&cbk,buffer,'+',string(argv[i]));
            if(rev_lookup) {
                auto_ptr<char> rev(FlipSequence(argv[i]));
                lookup.Find(&cbk,buffer,'-',string(rev.get()));
            }
            
        }
        
    } else {
        
        for(int i=optind; i<argc; ++i) {
			if(FILE * in = fopen(argv[i], "r")) {
                
				char buffer[4096];
				while(fgets(buffer,sizeof(buffer),in)) {
					if(feof(in)) break;
					buffer[sizeof(buffer)-1]=0;
					if(*buffer == '#') continue;
					const char * c = buffer;
					while(*c && isspace(*c)) ++c;
					if(*c == 0) continue;
					const char * cc = c;
					while(*cc && !isspace(*cc)) ++cc;
					string pname(c,cc);
					c = cc;
					while(*c && isspace(*c)) ++c;
					cc = c;
					while(*cc && !isspace(*cc)) ++cc;
					if(cc == c) continue;
					string primer(c,cc);
					lookup.Find(&cbk,pname,'+',primer);
                    if(rev_lookup) {
                        auto_ptr<char> rev(FlipSequence(primer.c_str()));
                        lookup.Find(&cbk,pname,'-',string(rev.get()));
                    }
				}	

				fclose(in);
			} else {
				fprintf(stderr,"! Could not open file %s -- ignored\n",
                        argv[i]);
			}
		}	
	}
	
	return 0;
}

int CMain::STSLookup()
{
	CFaLookup lookup;
    lookup.SetAligner(m_AlignL, m_AlignR);
	lookup.AttachFile(findex);
    FILE * out=outfile.length()?fopen64(outfile.c_str(),"w"):stdout;
	CFaLookupStsCallback cbk(out,m_PrintAlignments,m_Quiet);

    if(mode==eCmdline) {
        
        list<ISts*> stslist;
        for(int i=optind; i<=argc-3; i+=3) {
            char buffer[1024];
            snprintf(buffer,sizeof(buffer),"STS-%d",(i-optind+1));

            char * x=0;
            int lo=strtol(argv[i+2],&x,10);
            int hi=(x&&(*x=='-'))?max(atoi(x+1),lo):lo;
            if(lo==hi && lo==0) { lo=defLo; hi=defHi; }

            char* rev2(FlipSequence(argv[i+1],strlen(argv[i+1])));
            stslist.push_back(new CSimpleSTS(argv[i],rev2,lo,hi,
                                             ISts::ePlus,buffer));
            delete[] rev2;
            if(rev_lookup) {
                char* rev1(FlipSequence(argv[i+0],strlen(argv[i+0])));
                stslist.push_back(new CSimpleSTS(argv[i+1],rev1,lo,hi,
                                                 ISts::eMinus,buffer));
            }
        }
        lookup.Find(&cbk,stslist,optimize,margin);
    } else {
        unsigned tcount=0;
        for(int i=optind; i<argc; i++) {
            FILE * f=fopen(argv[i],"r");
            if(!f) throw runtime_error(string(argv[i])+": "+strerror(errno));
            char buffer[4096];
            list<ISts*> stslist;
//            fprintf(stderr,"* File %s\n",argv[i]);
            unsigned fcount=0;
            
            while(fgets(buffer,sizeof(buffer),f)) {
                if(feof(f)) break;
                if(*buffer=='#') continue;

                CStrRef fld[5];
                int cnt=CMmFileSts::Parse(buffer,fld);
                if(cnt<3) 
                    throw runtime_error("format error in file "+
                                        string(argv[i]));
                int lo=defLo, hi=defHi;
				
                if(cnt>3)
                    CMmFileSts::ParseRange(fld[3],lo,hi);

                char* rev2(FlipSequence(fld[2].data(),
                                            fld[2].length()));

                string name(fld[0]);
                string lprimer(fld[1]);
                string rprimer(rev2);
                delete[] rev2;
                
                stslist.push_back(new CSimpleSTS(lprimer,rprimer,
                                                 lo,hi,ISts::ePlus,name));

                if(rev_lookup) {
                    char* rev1(FlipSequence(fld[1].data(),
                                                fld[1].length()));
                    lprimer.assign(fld[2]);
                    rprimer.assign(rev1);
                    delete[] rev1;
                    stslist.push_back(new CSimpleSTS(lprimer,rprimer,
                                                     lo,hi,ISts::eMinus,name));
                }
                if(stslist.size()>=batchcnt) {
                    fprintf(stderr,"= File %s, STSes %u-%u\n",
                            argv[i],fcount+1,fcount+stslist.size());
                    
                    lookup.Find(&cbk,stslist,optimize,margin);
                    for(list<ISts*>::const_iterator i=stslist.begin(); 
                        i!=stslist.end(); ++i) {
                        delete *i;
                    }
                    fcount+=stslist.size();
                    
                    stslist.clear();
                }
            }
            fclose(f);
            fprintf(stderr,"= File %s, STSes %u-%u\n",
                    argv[i],fcount+1,fcount+stslist.size());
            
            lookup.Find(&cbk,stslist,optimize,margin);
            for(list<ISts*>::const_iterator i=stslist.begin(); 
                i!=stslist.end(); ++i) {
                        delete *i;
            }
            stslist.clear();
        }
    }
	return 0;
}

int CMain::Version()
{
	done=true;
	puts("Reverse e-PCR search cmdline tool version " VERSION);
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
 * $Log: re-PCR_main.cpp,v $
 * Revision 1.21  2008/04/28 16:38:45  rotmistr
 * Applied patch to build with gcc-4.3
 *
 * Revision 1.20  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.19  2005/04/11 18:17:52  rotmistr
 * Fixed range-parsing code for commandline
 *
 * Revision 1.18  2005/02/11 20:42:54  rotmistr
 * Fixed "margin" bug, added primer search from file
 *
 * Revision 1.17  2004/09/03 21:28:50  rotmistr
 * Fixes to compile with Borland C++ 5.5
 *
 * Revision 1.16  2004/06/07 16:24:57  rotmistr
 * Bug fixes to previos version.
 *
 * Revision 1.15  2004/06/03 23:37:21  rotmistr
 * New aligner added.
 *
 * Revision 1.14  2004/05/27 20:35:47  rotmistr
 * Version 2.1.0 with appropriate changes (see Changes) is ready for tests.
 *
 * Revision 1.13  2004/04/06 04:53:18  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.12  2004/03/30 21:06:53  rotmistr
 * Fixes for setting default STS size range.
 *
 * Revision 1.11  2004/03/29 21:25:40  rotmistr
 * Dist files are prepared
 *
 * Revision 1.10  2004/02/05 23:41:21  rotmistr
 * Better reload, fixed margin report in commandline, unists tab in CGI form.
 *
 * Revision 1.9  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.8  2004/01/28 23:27:02  rotmistr
 * "Best of overlapping" hit selection postprocessor added.
 *
 * Revision 1.7  2004/01/08 23:22:41  rotmistr
 * Fixed init error in faread,
 * Adjusted output to standard,
 * Added output format style and output file to parameters.
 *
 * Revision 1.6  2004/01/07 16:57:42  rotmistr
 * Fragment size is now configurable.
 *
 * Revision 1.5  2004/01/06 21:54:19  rotmistr
 * Statistics for word repetitions API added
 *
 * Revision 1.4  2003/12/30 21:36:32  rotmistr
 * Syscall optimisation mode added.
 *
 * Revision 1.3  2003/12/30 15:27:22  rotmistr
 * Fixed bug with sequence end
 *
 * Revision 1.2  2003/12/23 21:30:50  rotmistr
 * - gaps/mismatches reporting
 * - lo/hi fixup
 * - reverse sts in re-PCR_main
 *
 * Revision 1.1.1.1  2003/12/23 18:17:28  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
