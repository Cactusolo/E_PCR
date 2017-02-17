/* $Id: e-PCR_main.cpp,v 1.25 2008/06/18 14:45:33 rotmistr Exp $
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
#include <epcr/faread.hpp>

#include <epcr/stsmatch.hpp>

#include <assert.h>
#include <errno.h>
#include <time.h>

#include <cstdio>
#include <memory>
#include <sstream>

#ifndef STANDALONE
//#include <objtools/readers/seqdb/seqdb.hpp>
//#include <objtools/readers/seqdb/seqdbcommon.hpp>
#include <objtools/blast/seqdb_reader/seqdb.hpp>
#include <objtools/blast/seqdb_reader/seqdbcommon.hpp>
#endif

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

class CMain
{
public:
    enum EAlignMode { eNever, eAlways, eFallback };
	CMain(int c, char ** v):
		argc(c),argv(v),done(false),ofmt(1),
        postprocess(true),have_postprocess(false),verbose(0),
        m_MaxMismatch(ePCR_MMATCH_DEFAULT),
        m_MaxGaps(ePCR_GAPS_DEFAULT),
        m_AlignL(0), m_AlignR(0), m_AlignMode(eNever)
#ifndef STANDALONE
		, m_blastdbs( false ) 
#endif
	{
        stsFileHash.SetHash(CHashSet(ePCR_WDSIZE_DEFAULT,0));
        pcrMachine.SetMargin(ePCR_MARGIN_DEFAULT);
    } 
	int Run();
protected:
	int Execute();
	int ParseCmdline();
	int Help(FILE* = stdout);
	int Version();
    void ParseVerbose(const char * opt);
protected:
	int argc;
	char ** argv;
	bool done;
protected:
//	CPcrMachineCompat pcrMachine;
    CPcrMachine pcrMachine;
    CStsFileHash stsFileHash;
    
	string stsfile;
	list<string> fafiles;

    int ofmt;
    string ofile;

    bool postprocess, have_postprocess;
    int verbose;
    int m_MaxMismatch, m_MaxGaps;
    IAlign * m_AlignL, * m_AlignR;
    EAlignMode m_AlignMode;
#ifndef STANDALONE
	bool m_blastdbs;
	string m_gilist;
#endif
};

int CMain::Help(FILE* out)
{
	done=true;
	fprintf(out,
			"usage: [-hV] [posix-options] stsfile [fasta ...] "
			"[compat-options]\n"
			"where posix-options are:\n");
	fprintf(out,"\t-m ##\tMargin (default %d)\n",ePCR_MARGIN_DEFAULT);
	fprintf(out,"\t-w ##\tWordsize  (default %d)\n",ePCR_WDSIZE_DEFAULT);
	fprintf(out,"\t-n ##\tMax mismatches allowed (default %d)\n",
			ePCR_MMATCH_DEFAULT);
	fprintf(out,"\t-g ##\tMax indels allowed (default %d)\n",
			ePCR_GAPS_DEFAULT);
	fprintf(out,"\t-f ##\tUse ## discontiguos words, slow if ##>1\n");
	fprintf(out,"\t-o ##\tSet output file\n");
	fprintf(out,"\t-t ##\tSet output format:\n"
            "\t\t1 - classic, range (pos1..pos2)\n"
            "\t\t2 - classic, midpoint\n"
            "\t\t3 - tabular\n"
            "\t\t4 - tabular with alignment in comments (slow)\n"
        );
	fprintf(out,"\t-d##-##\tSet default size range (default %d-%d)\n",
			ePCR_DEFAULT_size_lo,ePCR_DEFAULT_size_hi);
	fprintf(out,"\t-p +-\tTurn hits postprocess on/off\n");
    fprintf(out,"\t-v ##\tVerbosity flags\n");
    fprintf(out,"\t-a a|f\tUse presize alignmens (only if gaps>0), slow\n"
            "\t\t a - Allways or f - as Fallback\n");
    fprintf(out,"\t-x +-\tUse 5'-end lowercase masking of primers "
            "(default %s)\n",stsFileHash.AllowOverhang()?"+":"-");
    fprintf(out,"\t-u +-\tUppercase all primers "
            "(default %s)\n",stsFileHash.UnmaskPrimers()?"+":"-");
#ifndef STANDALONE
	fprintf(out,"\t-b +-\tInput sequences are in blastdb\n");
	fprintf(out,"\t-l file\tLimit blastdb sequences to list of gis from the file\n");
#endif
	fprintf(out,"and compat-options (duplicate posix-options) are:\n");
	fprintf(out,"\tM=##\tMargin (default %d)\n",ePCR_MARGIN_DEFAULT);
	fprintf(out,"\tW=##\tWordsize  (default %d)\n",ePCR_WDSIZE_DEFAULT);
	fprintf(out,"\tN=##\tMax mismatches allowed (default %d)\n",
			ePCR_MMATCH_DEFAULT);
	fprintf(out,"\tG=##\tMax indels allowed (default %d)\n",
			ePCR_GAPS_DEFAULT);
	fprintf(out,"\tF=##\tUse ## discontinuos words\n");
	fprintf(out,"\tO=##\tSet output file to ##\n");
	fprintf(out,"\tT=##\tSet output format (1..3)\n");
	fprintf(out,"\tD=##-##\tSet default size range\n");
	fprintf(out,"\tP=+-\tPostprocess hits on/off\n");
	fprintf(out,"\tV=##\tVerbosity flags\n");
    fprintf(out,"\tA=a|f\tUse presize alignmens (only if gaps>0), slow\n"
            "\t\t a - Allways or f - as Fallback\n");
    fprintf(out,"\tX=+-\tUse 5'-end lowercase masking of primers "
            "(default %s)\n",stsFileHash.AllowOverhang()?"+":"-");
    fprintf(out,"\tU=+-\tUppercase all primers "
            "(default %s)\n",stsFileHash.UnmaskPrimers()?"+":"-");
#ifndef STANDALONE
	fprintf(out,"\tB=+-\tInput sequences are in blastdb\n");
	fprintf(out,"\tL=file\tLimit blastdb sequences to list of gis from the file\n");
#endif
	fprintf(out,"\t-mid\tSame as T=2\n");
    fprintf(out,"verbosity flags are (flags may be changed in future):\n"
            "\t-  set all progress reporting off (default)\n"
            "\t+  switch error reporting to basic (same as Sl)\n"
            "\tt  display time\n"
            "\tl  display fasta identifiers\n"
            "\to  display sequence offset (currently: 3' position of first primer)\n"
            "\tp  display percent of sequence processed\n"
            "\ts  report every sequence start\n"
            "\te  report every sequence end\n"
            "\tS  newline after sequence start report\n"
            "\tE  newline after sequence end report\n"
            "\tP  newline after sequence progress report\n");

	return 0;    
}

void SetDefaultSize(CStsFileHash& stsFileHash, const char * str) 
{

	char * x=0;
	int hi=0, lo=strtol(str,&x,10);
	if(x!=0 && *x=='-') {
		hi=atoi(x+1);
	}
	else hi=lo;
	if(lo>0 && hi>=lo) stsFileHash.SetDefaultSize(lo,hi);
	else throw runtime_error("bad range: "+string(str));
}

int CMain::ParseCmdline()
{
	int optchar;
	while((optchar=getopt(argc,argv,"+hVf:m:n:w:g:o:t:p:v:d:a:x:u:b:l:"))!=-1) {
		switch(optchar) {
		case 'h': Help(); break;
		case 'V': Version(); break;
        case 'a': 
            m_AlignMode=(optarg[0]=='a'?eAlways:
                         optarg[0]=='f'?eFallback:eNever); 
            if(m_AlignMode==eNever) {
                fprintf(stderr,"? Unknown alignment mode `%s' ignored\n", 
                        optarg);
            }
            break;
		case 'm': 
			if(strcmp(optarg,"id")==0) ofmt=2;
			else pcrMachine.SetMargin(atoi(optarg)); 
			break;
 		case 'w': stsFileHash.SetHash(CHashSet(atoi(optarg), stsFileHash.GetWordCount())); break;
 		case 'n': m_MaxMismatch=atoi(optarg); break;
 		case 'g': m_MaxGaps=atoi(optarg); break;
 		case 'f': stsFileHash.SetHash(CHashSet(stsFileHash.GetWordSize(), atoi(optarg))); break;
 		case 'd': SetDefaultSize(stsFileHash,optarg); break;
        case 'o': ofile=optarg; break;
        case 't': ofmt=atoi(optarg); break;
        case 'p': have_postprocess=true;
            postprocess=*optarg=='+'?true:*optarg=='-'?false:postprocess;
            break;
#ifndef STANDALONE
        case 'b': m_blastdbs=*optarg=='+'?true:*optarg=='-'?false:m_blastdbs;
            break;
		case 'l': m_gilist = optarg; break;
#endif
        case 'v': 
            ParseVerbose(optarg);
            break;
        case 'x': 
            if(*optarg=='+') 
                stsFileHash.SetFlags(CStsFileHash::fAllowOverhang,true);
            else if(*optarg=='-') 
                stsFileHash.SetFlags(CStsFileHash::fAllowOverhang,false);
            break;
        case 'u': 
            if(*optarg=='+') 
                stsFileHash.SetFlags(CStsFileHash::fUnmaskPrimers,true);
            else if(*optarg=='-') 
                stsFileHash.SetFlags(CStsFileHash::fUnmaskPrimers,false);
            break;
		}
	}
	if(done) return 0;
	if(optind >= argc) { Help(stderr); return 1; }
	// Parse compat options
 	for(; optind<argc; ++optind) {
 		if(argv[optind][1]!='=') { 
            if(strcmp(argv[optind],"-mid")==0) { ofmt=2; continue; }
			if(stsfile.length()==0) { 
				stsfile=argv[optind]; 
				if(stsfile == "-") { Help(stderr); return 1; }
				continue; 
			} else {
				fafiles.push_back(argv[optind]);
				continue;
			}
		
// 			Help(stderr); 
// 			return 2; 
		}
 		switch(argv[optind][0]) {
 		case 'M': pcrMachine.SetMargin(atoi(argv[optind]+2)); break;
  		case 'W': stsFileHash.SetHash(CHashSet(atoi(argv[optind]+2),
                                               stsFileHash.GetWordCount())); 
        break;
        case 'N': m_MaxMismatch=atoi(argv[optind]+2); break;
  		case 'G': m_MaxGaps=atoi(argv[optind]+2); break;
  		case 'F': stsFileHash.SetHash(CHashSet(stsFileHash.GetWordSize(),
                                               atoi(argv[optind]+2))); break;
        case 'O': ofile=argv[optind]+2; break;
        case 'T': ofmt=atoi(argv[optind]+2); break;
        case 'D': SetDefaultSize(stsFileHash,argv[optind]+2); break;
        case 'P': 
            postprocess=
                argv[optind][2]=='+'?true:
                argv[optind][2]=='-'?false:postprocess;
            break;
        case 'V': ParseVerbose(argv[optind]+2); break;
        case 'A': 
            m_AlignMode=(tolower(argv[optind][2])=='a'?eAlways:
                         tolower(argv[optind][2])=='f'?eFallback:eNever); 
            if(m_AlignMode==eNever) {
                fprintf(stderr,"? Unknown alignment mode `%s' ignored\n", 
                        argv[optind]+2);
            }
            break;
#ifndef STANDALONE
        case 'B': 
            if(argv[optind][2]=='+') 
                m_blastdbs = true;
            else if(argv[optind][2]=='-') 
                m_blastdbs = false;
            break;
		case 'L': m_gilist = argv[optind]+2; break;
#endif
        case 'X': 
            if(argv[optind][2]=='+') 
                stsFileHash.SetFlags(CStsFileHash::fAllowOverhang,true);
            else if(argv[optind][2]=='-') 
                stsFileHash.SetFlags(CStsFileHash::fAllowOverhang,false);
            break;
        case 'U': 
            if(argv[optind][2]=='+') 
                stsFileHash.SetFlags(CStsFileHash::fUnmaskPrimers,true);
            else if(argv[optind][2]=='-') 
                stsFileHash.SetFlags(CStsFileHash::fUnmaskPrimers,false);
            break;
		default:
			fprintf(stderr,"No option \"%c=##\" allowed\n",argv[optind][0]);
			break;
 		}
 	}
	if(stsfile.length()==0) {
		Help(stderr);
		return 2;
	}
	return 0;
}

class CPcrMachineCallback:public IPcrMachineCallback
{
public:
    CPcrMachineCallback(FILE * out):
        m_Out(out) {}
    
    virtual ~CPcrMachineCallback() throw () {}
    
    virtual void CbkSequence(const char * label) { m_SeqId=label; }
    virtual void CbkWarning(const char * ) {}
protected:
    string m_SeqId;
    FILE * m_Out;
};

class CPcrProgressCallback: public IPcrProgressCallback
{
public:
    enum EFlags {
        fSequenceLabel    = 1,
        fSequencePosition = 2,
        fSequencePercent  = 4,
        fCurrentTime      = 8,
        fSequenceStart    = 0x10,
        fSequenceEnd      = 0x20,
        fSequenceProgress = ( fSequencePosition | fSequencePercent ),
        fNewlineAtSequenceStart = 0x100,
        fNewlineAtSequenceEnd   = 0x200,
        fNewlineAtProgress      = 0x400,
        fNoBeol                 = 0x800,
        fNoTTY = ( fNoBeol |
                   fNewlineAtSequenceStart |
                   fNewlineAtSequenceEnd |
                   fNewlineAtProgress ),
        fBasic = ( fSequenceLabel | fSequenceStart | fNewlineAtSequenceStart | fNoBeol ),
        fNone = 0
    };
    CPcrProgressCallback(int flags):m_Flags(flags),m_Timer(time(0)) {
        if(! isatty(fileno(stderr))) m_Flags |= fNoTTY;
    }
    void PgsSequenceStart(const char * label, const char * data, unsigned len, unsigned wsize) {
        m_Label = label; m_SeqData = data; m_SeqLen = len; m_WdSize = wsize; 
        if(m_Flags & fSequenceStart) Report(0,true);
    }
    void PgsSequenceEnd() { if(m_Flags & fSequenceEnd) Report(m_SeqLen,true); }
    void PgsSequenceAt(unsigned pos) { if(m_Flags & fSequenceProgress) Report(pos,false); }
    void Report(unsigned pos, bool force) {
        time_t t = time(0);
        if(t != m_Timer) { m_Timer = t; } else if (!force) return;
        if(!(m_Flags & fNoBeol)) fputs("\r\033[K* ",stderr);
        
        if(m_Flags&fCurrentTime) {
            char buffer[1024];
            strftime(buffer,sizeof buffer,"%Y/%m/%d %H:%M:%S ",localtime(&m_Timer));
            fputs(buffer,stderr);
        }
        if(m_Flags&fSequenceLabel) fprintf(stderr,"%s ",m_Label);
        if(pos) {
            if(pos != m_SeqLen) {
                if(m_Flags&fSequenceProgress ) {
                    if(m_Flags&fSequencePosition) fprintf(stderr,"%u/%u ",pos,m_SeqLen);
                    if(m_Flags&fSequencePercent) fprintf(stderr,"%5.1lf%% ",(pos*100.0)/m_SeqLen);
                    if(m_Flags&fNewlineAtProgress) putc('\n',stderr);
                }
            } else {
                fputs("Done",stderr);
                if(m_Flags&fNewlineAtSequenceEnd) putc('\n',stderr);
            }
        } else 
            if(m_Flags&fNewlineAtSequenceStart) putc('\n',stderr);
        fflush(stderr);
    }
            
protected:
    int m_Flags;
    time_t m_Timer;
    
    const char * m_Label;
    const char * m_SeqData;
    unsigned m_SeqLen;
    unsigned m_WdSize;
};

void CMain::ParseVerbose(const char * opt)
{
    while(*opt) {
        switch(*opt++) {
        case '-': verbose = CPcrProgressCallback::fNone; break;
        case '+': verbose = CPcrProgressCallback::fBasic; break;
        case 'l': verbose |= CPcrProgressCallback::fSequenceLabel;    break;
        case 'o': verbose |= CPcrProgressCallback::fSequencePosition; break;
        case 'p': verbose |= CPcrProgressCallback::fSequencePercent;  break;
        case 't': verbose |= CPcrProgressCallback::fCurrentTime;      break;
        case 'S': verbose |= CPcrProgressCallback::fNewlineAtSequenceStart;
        case 's': verbose |= CPcrProgressCallback::fSequenceStart;    break;
        case 'E': verbose |= CPcrProgressCallback::fNewlineAtSequenceEnd; 
        case 'e': verbose |= CPcrProgressCallback::fSequenceEnd;    break;
        case 'P': verbose |= CPcrProgressCallback::fNewlineAtProgress; break;
        default: fprintf(stderr,"? Unknown verbosity flag `%c'\n",opt[-1]);
            break;
        }
    }
}

class CPcrMachineCallbackClassic:public CPcrMachineCallback
{
public:
    CPcrMachineCallbackClassic(FILE * out, bool midpnt=false):
        CPcrMachineCallback(out),show_midpt(midpnt) {}
    virtual ~CPcrMachineCallbackClassic() throw () {}
    virtual void CbkMatch(const ISts * sts, unsigned pos1, unsigned pos2,
                          const SScore*) 
        { ReportHit(m_SeqId.c_str(),pos1,pos2,sts); }
protected:
    virtual int  ReportHit(const char *seq_label, int pos1, int pos2,
                           const ISts *sts);
protected:
    bool show_midpt;
};

int CPcrMachineCallbackClassic::ReportHit (
    const char *seq_label, 
    int pos1, int pos2,
    const ISts *sts) 
{
    char position[32];
    
    //int ovhg1 = sts->GetOverhangChars(ISts::eLeft);
    //int ovhg2 = sts->GetOverhangChars(ISts::eRight);

//     pos1 += ovhg1;
//     pos2 -= ovhg2;

    if (show_midpt)
        sprintf(position,"%d", 1 + (pos1+pos2-1)/2);
    else
        sprintf(position,"%d..%d",pos1+1,pos2);
    
    fprintf(m_Out,"%-10s %-16s %-14.*s  %.*s\n",
            seq_label,position,
            sts->GetName().length(),sts->GetName().data(),
            sts->GetDescription().length(),
            sts->GetDescription().data());

	return 1;
}

class CPcrMachineCallbackTabular:public CPcrMachineCallback
{
public:
    CPcrMachineCallbackTabular(FILE * out, bool showalign, int gaps):
        CPcrMachineCallback(out),
        m_ShowAlign(showalign),m_Matrix(127,gaps),m_SeqData(0),m_SeqLength(0) {}
    virtual ~CPcrMachineCallbackTabular() throw () {}
    virtual void CbkMatch(const ISts * sts, unsigned pos1, unsigned pos2,
                          const SScore* score) ;
    virtual void CbkSequenceData(const char * data, unsigned size) {
//         delete[] m_SeqData;
//         m_SeqData=new char[(m_SeqLength=size)+1];
//         memcpy(m_SeqData,data,size);
//         m_SeqData[size]=0;
        m_SeqData=data;
        m_SeqLength=size;
    }
protected:
    bool m_ShowAlign;
    CLcsMatrix<char> m_Matrix;
    const char * m_SeqData;
    unsigned m_SeqLength;
};

void CPcrMachineCallbackTabular::CbkMatch (
    const ISts * sts,
    unsigned pos1, unsigned pos2,
    const SScore* score) 
{
    int mism=score->mism_l+score->mism_r;
    int gaps=score->gaps_l+score->gaps_r;
    
    int len1=sts->GetPrimerLength(ISts::eLeft);
    int len2=sts->GetPrimerLength(ISts::eRight);
    int ovhg1=sts->GetOverhangChars(ISts::eLeft);
    int ovhg2=sts->GetOverhangChars(ISts::eRight);

    const char * data1=sts->GetPrimerData(ISts::eLeft);
    const char * data2=sts->GetPrimerData(ISts::eRight);

    if(m_ShowAlign && m_SeqData && pos2<=m_SeqLength) {
        
        vector<string> left, right;

        m_Matrix.Build<const char*>(
            m_SeqData+pos2-len2-ovhg2,m_SeqData+m_SeqLength-ovhg2,data2,len2);
        m_Matrix.Graph<const char*>(
            m_SeqData+pos2-len2-ovhg2,m_SeqData+m_SeqLength-ovhg2,data2,len2,
            right);
        m_Matrix.Stat<const char*>(
            m_SeqData+pos2-len2-ovhg2,m_SeqData+m_SeqLength-ovhg2,data2,len2);

        mism = m_Matrix.GetMismatches();
        gaps = m_Matrix.GetGaps();

        m_Matrix.Build<CReverseConstSeqIterator<const char> >(
            m_SeqData+pos1+len1+ovhg1-1,m_SeqData,data1+len1-1,len1);
        m_Matrix.Graph<CReverseConstSeqIterator<const char> >(
            m_SeqData+pos1+len1+ovhg1-1,m_SeqData,data1+len1-1,len1,left);
        m_Matrix.Stat<CReverseConstSeqIterator<const char> >(
            m_SeqData+pos1+len1+ovhg1-1,m_SeqData,data1+len1-1,len1);

        mism += m_Matrix.GetMismatches();
        gaps += m_Matrix.GetGaps();
        
        int l=max(int(m_SeqId.length()),int(sts->GetName().length()));
        int d=score->actlen-len1-len2;
        
        string stsname(sts->GetName().data(),sts->GetName().length());
        fprintf(m_Out,
                "#####################################"
                "#####################################\n"
                "#   STS %*s %s...%d...%s\n"
                "#       %.*s %s   %d   %s\n"
                "#   Seq %*s %s...%d...%s\n",
                l,stsname.c_str(),
                left[0].c_str(),d,right[0].c_str(),
                l,"                                                   ",
                left[2].c_str(),d,right[2].c_str(),
                l,m_SeqId.c_str(),
                left[1].c_str(),d,right[1].c_str());
    }

    pos1 += ovhg1;
    pos2 -= ovhg2;
    
	fprintf(m_Out,"%s\t%.*s\t%c\t%d\t%d\t%d/%d-%d\t%d\t%d\t%.*s\n",
            m_SeqId.c_str(),
            sts->GetName().length(),sts->GetName().data(),
            sts->GetDirection(),
            pos1+1,pos2,
            score->actlen,
            sts->GetSizeLo(),sts->GetSizeHi(),
            mism, gaps,
            sts->GetDescription().length(),
            sts->GetDescription().data());
}

class CPcrFastaProcessor:public IFastaReaderCallback
{
public:
    virtual ~CPcrFastaProcessor() throw () { if( !m_NoCopySeq ) free(m_Sequence); }
    CPcrFastaProcessor(CPcrMachine* pmachine, bool noCopySeq = false ):
        m_Sequence(0),m_Size(0),m_Capacity(0), m_NoCopySeq( noCopySeq )
        { m_PcrMachine=pmachine; }
    virtual void CbkDefline(const char * , unsigned ) {}
    virtual void CbkIdent(const char * ident, unsigned length);
    virtual void CbkSeqline(const char * data, unsigned length);
    virtual void CbkEntryEnd();
protected:
    CPcrMachine * m_PcrMachine;
    string m_Ident;
    char * m_Sequence;
    unsigned m_Size;
    unsigned m_Capacity;
	bool m_NoCopySeq;
};

void CPcrFastaProcessor::CbkIdent(const char * ident, unsigned length)
{
    m_Ident.assign(ident,length);
}

void CPcrFastaProcessor::CbkSeqline(const char * seq, unsigned length)
{
	if( m_NoCopySeq ) {
		assert( m_Size == 0 && m_Sequence == 0 );
		m_Size = length;
		m_Sequence = const_cast<char*>(seq);
	} else {
	    while(m_Size+length>=m_Capacity) 
    	    m_Sequence=(char*)realloc(m_Sequence,m_Capacity+=16192);
	    memcpy(m_Sequence+m_Size,seq,length);
    	m_Size+=length;
	}
}

void CPcrFastaProcessor::CbkEntryEnd()
{
	if( !m_NoCopySeq ) { if(m_Sequence) m_Sequence[m_Size]=0; }
   	m_PcrMachine->ProcessSequence(m_Ident.c_str(),m_Sequence,m_Size);
	m_Size=0;
	if( m_NoCopySeq ) { m_Sequence = 0; }
    m_Ident.clear();
}

int CMain::Execute()
{
	stsFileHash.SetOneTimeRun(true);
	
	do {
		CStsFileCallbackDefault cbk;
		stsFileHash.ReadStsFile(stsfile, &cbk);
	} while(0);

    do {
        FILE * out=ofile.length()?fopen64(ofile.c_str(),"w"):stdout;
        if(out==0) 
            throw runtime_error(ofile+": "+strerror(errno));
            
        auto_ptr<CPcrMachineCallback> cbk(0);
        switch(ofmt) {
        case 4:
        case 3: 
            cbk.reset(new CPcrMachineCallbackTabular(out,ofmt==4,m_MaxGaps));
            break;
        case 2: 
        case 1: 
        default: 
            cbk.reset(new CPcrMachineCallbackClassic(out,ofmt==2));
            break;
        }
        
//         if(!have_postprocess) {
//             if(pcrMachine.GetMaxIndels() || 
//                stsFileHash.GetHash().GetWordCount()>1) {
//                 postprocess=true;
//             } else {
//                 postprocess=false;
//             }
//         }

        CPcrMachinePostprocess post(cbk.get());
        if(postprocess)
            pcrMachine.SetCallback(&post);
        else
            pcrMachine.SetCallback(cbk.get());
        
        CPcrProgressCallback pgscbk(verbose);
        if(verbose) pcrMachine.SetProgressCallback(&pgscbk);
            
        pcrMachine.SetStsHash(&stsFileHash);

        if(m_MaxGaps) {
            switch(m_AlignMode) {
            case eNever:
                m_AlignL=new CAlignFast(m_MaxMismatch,m_MaxGaps);
                m_AlignR=new CAlignFast(m_MaxMismatch,m_MaxGaps);
                break;
            case eAlways:
                m_AlignL=new CAlignLCS(m_MaxMismatch,m_MaxGaps);
                m_AlignR=new CAlignLCS(m_MaxMismatch,m_MaxGaps);
                break;
            case eFallback:
                m_AlignL=new CAlignCompromise(m_MaxMismatch,m_MaxGaps);
                m_AlignR=new CAlignCompromise(m_MaxMismatch,m_MaxGaps);
                break;
            default:
                throw logic_error("Invalig align mode");
            }
        }
        else if(m_MaxMismatch) {
            m_AlignL=new CAlignNoGaps(m_MaxMismatch);
            m_AlignR=new CAlignNoGaps(m_MaxMismatch);
        }
        else {
            m_AlignL=new CAlignExact();
            m_AlignR=new CAlignExact();
        }

        pcrMachine.SetAligner(m_AlignL,m_AlignR);
        
        CPcrFastaProcessor processor(&pcrMachine
#ifndef STANDALONE 
				, m_blastdbs
#endif
				);
		if(fafiles.size()) {
			for(list<string>::const_iterator f=fafiles.begin(); 
				f!=fafiles.end(); ++f) {
#ifndef STANDALONE
				if( m_blastdbs ) {
                    vector<string> volumes;
                    try {
                        CSeqDB::FindVolumePaths( *f, CSeqDB::eNucleotide, volumes, 0, true );
                    } catch(exception& e) {
                        cerr << "? Warning: CSeqDB::FindVolumePaths( \"" << *f << "\", CSeqDB::eNucleotide, volumes, 0, true ); failed with error: " << e.what() << "\n";
                        volumes.clear();
                        volumes.push_back( *f );
                    } catch(...) {
                        cerr << "? Warning: CSeqDB::FindVolumePaths( \"" << *f << "\", CSeqDB::eNucleotide, volumes, 0, true ); failed with unknown exception\n";
                        volumes.clear();
                        volumes.push_back( *f );
                    }
                    for( vector<string>::const_iterator v = volumes.begin(); v != volumes.end(); ++v ) {
						auto_ptr<CSeqDB> seqDB( 0 );
                        try {
					        CSeqDBGiList * lst = (m_gilist.length() ? new CSeqDBFileGiList( m_gilist ) : 0 );
                            seqDB.reset( new CSeqDB( *v, CSeqDB::eNucleotide, lst ) );
                        } catch(exception& e) {
                            throw runtime_error( "Failed to open blastdb volume " + *v + ": " + e.what() );
                        } catch(...) {
                            throw runtime_error( "Failed to open blastdb volume " + *v + ": unknown error" );
                        }
						processor.CbkFileBegin();
						for( CSeqDBIter i = seqDB->Begin(); i; ++i ) {
							list<CRef<CSeq_id> > ids = seqDB->GetSeqIDs( i.GetOID() );
                            if( ids.size() == 0 ) {
                                ostringstream err;
                                err << "Bad entry in " << *f << " (" << *v << ") " << " oid " << i.GetOID() << ": no seqids\n";
                                throw runtime_error( err.str() );
                                //cerr << "? Warning: " << err.str();
                                //continue;
                            }
							string ident;
							for( list<CRef<CSeq_id> >::const_iterator x = ids.begin(); x != ids.end(); ++x ) {
								if( x != ids.begin() ) ident += "|";
								ident += (*x)->AsFastaString();
							}
							processor.CbkEntryBegin();
							processor.CbkIdent( ident.c_str(), ident.length() );
							string seq;
							seqDB->GetSequenceAsString( i.GetOID(), seq );
							processor.CbkSeqline( seq.c_str(), seq.length() );
							processor.CbkEntryEnd();
						}
                    }
					processor.CbkFileEnd();
				} else {
#endif
					if(*f=="-") {
						CFastaReader reader("/dev/stdin");
						reader.SetCvtTable(CFastaReader::sm_NucleotidesUc);
						reader.ReadFile(&processor);
					} else {
						CFastaReader reader(*f);
						reader.SetCvtTable(CFastaReader::sm_NucleotidesUc);
						reader.ReadFile(&processor);
					}
#ifndef STANDALONE
				}
#endif
			}
		}
		else {
			CFastaReader reader("/dev/stdin");
			reader.SetCvtTable(CFastaReader::sm_NucleotidesUc);
			reader.ReadFile(&processor);
		}

        delete m_AlignR;
        delete m_AlignL;

        fclose(out);
    } while(0);

	return 0;
}

int CMain::Version()
{
	done=true;
	puts("e-PCR cmdline tool version " VERSION);
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
//	try {
		CMain app(argc,argv);
		return app.Run();
        /*
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
    */
}


/*
 * $Log: e-PCR_main.cpp,v $
 * Revision 1.25  2008/06/18 14:45:33  rotmistr
 * Fixed problem with -d x-X parameter being reset if -w or some others are used after it.
 *
 * Revision 1.24  2008/06/16 16:02:40  rotmistr
 * *** empty log message ***
 *
 * Revision 1.23  2008/04/28 16:38:45  rotmistr
 * Applied patch to build with gcc-4.3
 *
 * Revision 1.22  2008/03/27 14:36:58  rotmistr
 * Added assert.h to make it compiling with VC8
 *
 * Revision 1.21  2008/03/26 16:04:29  rotmistr
 * Added support for blastdb files
 *
 * Revision 1.20  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.19  2005/06/14 16:46:41  rotmistr
 * Changed report format for floppy tails
 *
 * Revision 1.18  2004/10/26 17:16:33  rotmistr
 * Added 5'-end masking for primers
 *
 * Revision 1.17  2004/06/08 20:32:51  rotmistr
 * Fixup for gap+insert special case
 *
 * Revision 1.16  2004/06/08 16:14:55  rotmistr
 * *** empty log message ***
 *
 * Revision 1.15  2004/06/03 23:37:19  rotmistr
 * New aligner added.
 *
 * Revision 1.14  2004/04/06 04:53:17  rotmistr
 * All is compileable with BCC5.5 and runnable on WIndows
 *
 * Revision 1.13  2004/04/01 16:37:41  rotmistr
 * Cleaned after adding windows capabilities
 *
 * Revision 1.12  2004/04/01 05:57:52  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.11  2004/03/30 21:06:53  rotmistr
 * Fixes for setting default STS size range.
 *
 * Revision 1.10  2004/03/30 19:11:18  rotmistr
 * STS default size
 *
 * Revision 1.9  2004/03/30 19:08:03  rotmistr
 * default STS size is tunnable now
 *
 * Revision 1.8  2004/03/26 17:02:13  rotmistr
 * Compat-options are now allowed everywhere, and multiple fasta files can be used.
 *
 * Revision 1.7  2004/03/25 19:36:52  rotmistr
 * API: separate left and right primers mism/gaps in forward API
 *
 * Revision 1.6  2004/03/23 22:35:25  rotmistr
 * Fixed processing of -mid flag in cmdline
 * Fixed destructor for fasta reader
 * Removed cgi
 *
 * Revision 1.5  2004/03/07 06:35:59  rotmistr
 * Many bugfixes and optimisations -- cgi is to go to production
 *
 * Revision 1.4  2004/02/04 21:23:22  rotmistr
 * - gcc-3.3.2 compatible
 * - better postfiltering for reverse-e-PCR for discontiguos words
 * - cgi added, that supports:
 *  -- contig to chromosome mapping
 *  -- simple mapviewer links
 *  -- unists links
 *  -- discontiguos words
 *
 * Revision 1.3  2004/01/28 23:27:02  rotmistr
 * "Best of overlapping" hit selection postprocessor added.
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
