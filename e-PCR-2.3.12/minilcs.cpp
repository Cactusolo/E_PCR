/* $Id: minilcs.cpp,v 1.5 2004/09/03 19:06:41 rotmistr Exp $
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

#include <iostream>
#include <stdexcept>
#include <vector>

#include "minilcs.hpp"

#include <unistd.h>

using namespace std;
using namespace ncbi;
using namespace EPCR_SCOPE;

int main(int argc, char ** argv) 
{
    try {
        
        bool done=false;
        bool reverse=false;
        
        int optopt;

        int maxgaps=2;
        
        while((optopt=getopt(argc,argv,"hVg:r"))!=-1) {
            switch(optopt) {
            case 'h': 
                cout << "Usage: [-hV] [-r] [-g gaps] sequence primer\n"; 
                done=true; 
                break;
            case 'V': 
                cout << "Version " __DATE__ "\n"; 
                done=true; 
                break;
            case 'g': 
                maxgaps=atoi(optarg);
                break;
            case 'r': reverse=true;
                break;
            }
        }
    
        if(done) return 0;
        if(argc-optind<2) throw runtime_error("Need primer and sequence");
        CLcsMatrix<char> matrix(256,maxgaps);
        vector<string> rc;
        double pcti;
        
        string sequence=argv[optind];
        string primer=argv[optind+1];

        if(!reverse) {
            
            pcti=matrix.Build<const char*>(
                sequence.c_str(),
                sequence.c_str()+sequence.length(),
                primer.c_str(),primer.length());
        
            matrix.Graph<const char*>(
                sequence.c_str(),
                sequence.c_str()+sequence.length(),
                primer.c_str(),primer.length(),rc);

            matrix.Stat<const char*>(
                sequence.c_str(),
                sequence.c_str()+sequence.length(),
                primer.c_str(),primer.length());

        } else {
        
            pcti=matrix.Build<CReverseConstSeqIterator<const char> >(
                sequence.c_str()+sequence.length()-1,
                sequence.c_str()-1,
                primer.c_str()+primer.length()-1,
                primer.length());
        
            matrix.Graph<CReverseConstSeqIterator<const char> >(
                sequence.c_str()+sequence.length()-1,
                sequence.c_str()-1,
                primer.c_str()+primer.length()-1,
                primer.length(),rc);
            matrix.Stat<CReverseConstSeqIterator<const char> >(
                sequence.c_str()+sequence.length()-1,
                sequence.c_str()-1,
                primer.c_str()+primer.length()-1,
                primer.length());
        }
        
        cout << "PCTI:   " << pcti << "\n";
        
        cout << "BestX:  " << matrix.GetBestX() << endl
             << "BestY:  " << matrix.GetBestY() << endl
             << "BestVal:" << matrix.GetBestVal() << endl;
        
        cout << "Mism:   " << matrix.GetMismatches() << endl
             << "Gaps:   " << matrix.GetGaps() << endl
             << "Match:  " << matrix.GetMatches() << endl;

        cout << "primer: " << rc[0] << endl
             << "        " << rc[2] << endl
             << "seqnce: " << rc[1] << endl;

        cout << "        ";
        for(int j=0; j<sequence.length(); ++j) {
            if(!reverse)
                cout << "  " << sequence[j] << "  |";
            else
                cout << "  " << sequence[sequence.length()-j-1] << "  |";
        }
        cout << "\n";
        for(int i=0; i<=matrix.GetBestX(); ++i) {
            if(i) {
                if(!reverse) {
                    cout << primer[i-1] << ": ";
                }
                else {
                    cout << primer[primer.length()-i-1] << ": ";
                }
            }
            else cout << " : ";
            for(int j=0; j<=sequence.length(); ++j) {
                cout.width(2);
                
                cout << matrix.Get(i,j);
                cout << matrix.Sym(i,j);
                cout << " | ";
            }
            cout << endl;
        }

        return 0;
    }
    catch(exception& e) {
        cerr << "! " << e.what() << endl;
    }
    catch(...) {
        cerr << "! Unknown error\n";
    }
    return 100;
}

#if 0
template<class T>
class CReverseConstSeqIterator
{
public:
    typedef CReverseConstSeqIterator<T> TClass;
    typedef CReverseConstSeqIterator<T> class_type;
    typedef T data_type;
    typedef T TDataType;
    
    CReverseConstSeqIterator(const T* ptr=0):m_Ptr(ptr) {}
    CReverseConstSeqIterator(const TClass& i):m_Ptr(i.m_Ptr) {}

    T* Get() const { return m_Ptr; }
    TClass& Set(const T* ptr=0) { m_Ptr=ptr; return *this; }
        
    T  operator * () const { return *m_Ptr; }
    T  operator [] (int i) const { return m_Ptr[-i]; }

    TClass& operator ++ () { --m_Ptr; return *this; }
    TClass& operator -- () { ++m_Ptr; return *this; }

    TClass  operator ++ (int) { TClass i(this); ++*this; return i; }
    TClass  operator -- (int) { TClass i(this); --*this; return i; }

    TClass& operator += (int i) { m_Ptr-=i; return *this; }
    TClass& operator -= (int i) { m_Ptr+=i; return *this; }

    TClass& operator  = (const TClass& i) { m_Ptr=i.m_Ptr; return *this; }

    bool    operator == (const TClass& i) const { return m_Ptr==i.m_Ptr; }
    bool    operator != (const TClass& i) const { return m_Ptr!=i.m_Ptr; }
    bool    operator >= (const TClass& i) const { return m_Ptr<=i.m_Ptr; }
    bool    operator <= (const TClass& i) const { return m_Ptr>=i.m_Ptr; }
    bool    operator <  (const TClass& i) const { return m_Ptr> i.m_Ptr; }
    bool    operator >  (const TClass& i) const { return m_Ptr< i.m_Ptr; }
    
    friend TClass operator + (const TClass& c, int i) { 
        return TClass(c.m_Ptr-i);
    }
    friend TClass operator - (const TClass& c, int i) { 
        return TClass(c.m_Ptr+i);
    }
protected:
    TDataType * m_Ptr;
};
        

template <class T>
class CLcsMatrix
{
public:
    typedef T data_type;
    typedef data_type TDataType;
    
    CLcsMatrix(int length, int maxgaps);
    ~CLcsMatrix() throw ();

    typedef const char * TForwardStr;
    typedef const char * TReverseStr;

    template<class CSeqIterator>
    int Build(CSeqIterator genome, CSeqIterator gend,
              CSeqIterator primer, int length);
    
protected:
    int m_Size;
    int m_MaxGaps;
    TDataType ** m_Data;
};

template<class T>
template<class CSeqIterator>
inline int CLcsMatrix<T>::Build(CSeqIterator genome, CSeqIterator gend,
                                CSeqIterator primer, int len)
{
    double best_pcti=0;
    for(int x=1; x <= len; ++x, ++genome, ++primer) {
        cout << *primer << ":\t";
        for(int y=-min(x-1,m_MaxGaps); y<=m_MaxGaps&&genome+y<gend; ++y) {
            if(genome[y]==*primer) m_Data[y][x]=m_Data[y][x-1]+1;
            else m_Data[y][x]=max(m_Data[y-1][x],m_Data[y+1][x-1]);
            if(x>len-m_MaxGaps) {
                double pcti=double(m_Data[y][x])/max(x,x+y);
                if(pcti>best_pcti) best_pcti=pcti;
            }
            cout.width(5);
            cout << (int)m_Data[y][x] << ", ";
        }
            
        cout << "\n";
    }
    cout << "pcti: " << 100*best_pcti << endl;
    
    return 0;
}

template<class T> 
inline CLcsMatrix<T>::CLcsMatrix(int length, int maxgaps):
    m_Size(length),
    m_MaxGaps(maxgaps),
    m_Data((new (TDataType*)[(1+m_MaxGaps)*2+1])+m_MaxGaps+1)
{
    for(int i=-m_MaxGaps-1; i<=m_MaxGaps+1; ++i) {
        m_Data[i]=new TDataType[m_Size+m_MaxGaps+1];
        memset(m_Data[i],0,(m_Size+m_MaxGaps+1)*sizeof(TDataType));
    }
}

template <class T>
inline CLcsMatrix<T>::~CLcsMatrix() throw () 
{
    for(int i=-m_MaxGaps-1; i<=m_MaxGaps+1; ++i) {
        delete[] m_Data[i];
    }
    m_Data-=m_MaxGaps+1;
    delete[] m_Data;
}

#endif

/*
 * $Log: minilcs.cpp,v $
 * Revision 1.5  2004/09/03 19:06:41  rotmistr
 * Code formatting changes
 *
 *.
