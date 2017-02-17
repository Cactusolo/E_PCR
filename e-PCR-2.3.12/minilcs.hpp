/* $Id: minilcs.hpp,v 1.9 2007/07/05 16:23:08 rotmistr Exp $
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

#ifndef EPCR_MINILCS__HPP
#define EPCR_MINILCS__HPP

#include <epcr/build_cfg.h>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cstring>

#include <ctype.h>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

template<class T>
class CReverseConstSeqIterator
{
public:
    typedef CReverseConstSeqIterator<T> TClass;
    typedef CReverseConstSeqIterator<T> class_type;
    typedef T data_type;
    typedef T TDataType;
protected:
    TDataType * m_Ptr;
public:
    
    CReverseConstSeqIterator(const T* ptr=0):m_Ptr(ptr) {}
    CReverseConstSeqIterator(const TClass& i):m_Ptr(i.m_Ptr) {}

    T* Get() const { return m_Ptr; }
    TClass& Set(const T* ptr=0) { m_Ptr=ptr; return *this; }
        
    const T& operator * () const { return *m_Ptr; }
    const T& operator [] (int i) const { return m_Ptr[-i]; }

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
};

template <class T>
class CLcsMatrix
{
public:
    typedef T data_type;
    typedef data_type TDataType;
    typedef TDataType * TDataTypePtr;
    
    CLcsMatrix(int length, int maxgaps);
    ~CLcsMatrix() throw ();

    typedef const char * TForwardStr;
    typedef CReverseConstSeqIterator<const char> TReverseStr;

    template<class CSeqIterator>
    double Build(CSeqIterator genome, CSeqIterator gend,
                 CSeqIterator primer, int length);

    template<class CSeqIterator>
    void Graph(CSeqIterator genome, CSeqIterator gend,
               CSeqIterator primer, int length, vector<string>& dest, 
               int extra=2, int intra=0);

    template<class CSeqIterator>
    void Stat(CSeqIterator genome, CSeqIterator gend,
              CSeqIterator primer, int length);
    
    int GetMatches() const { return m_Matches; }
    int GetMismatches() const { return m_Mismatches; }
    int GetSeqInsertions() const { return m_Insertions; }
    int GetSeqDeletions() const { return m_Deletions; }
    int GetGaps() const { return m_Insertions+m_Deletions; }
    int GetResultLength() const { return m_ResultLength; }
    int GetIdentities() const { return m_Identities; }
    int GetBestX() const { return m_BestX; }
    int GetBestY() const { return m_BestY; }
    int GetBestVal() const { return m_BestVal; }
    
    int Get(int i, int j) 
    {
        int y=j-i;
        return abs(y)>m_MaxGaps?0:m_Data[y][i];
    }
    char Sym(int i, int j) 
    {
        int y=j-i;
        char c=abs(y)>m_MaxGaps?0:m_Path[y][i];
        return c?c:'.';
    }

protected:
    int m_Size;
    int m_MaxGaps;
    TDataType ** m_Data;
    TDataType ** m_Path;
    int m_BestX, m_BestY, m_BestVal;
    int m_PrimerLength;
    int m_Matches, m_Mismatches, m_Insertions, m_Deletions;
    int m_Identities, m_ResultLength;
};

template<class T>
template<class CSeqIterator>
inline double CLcsMatrix<T>::Build(CSeqIterator genome, CSeqIterator gend,
                                   CSeqIterator primer, int len)
{
    m_PrimerLength=len;
    m_ResultLength=m_Identities=m_BestY=m_BestX=0;
    CSeqIterator gg=genome;
    
    for(int x=1; x <= len; ++x, ++gg, ++primer) {
        int start=-min(x-1,m_MaxGaps);
        CSeqIterator g=gg+start;
        for(int y=start; y<=m_MaxGaps && g<gend; ++y, ++g) {
            if(*g==*primer) {
                m_Data[y][x]=m_Data[y][x-1]+1;
                m_Path[y][x]='=';
            }
            else {
                char cy=m_Data[y-1][x];
                char cx=m_Data[y+1][x-1];
                if(cx>cy) {
                    m_Data[y][x]=cx;
                    m_Path[y][x]='i';
                } else if(cy>cx) {
                    m_Data[y][x]=cy;
                    m_Path[y][x]='d';
                } else {
//                    char cz=m_Data[y][x-1];
                    m_Data[y][x]=cx;//max(cx,cz);
                    m_Path[y][x]='!';
                }
            }
        }
    }


    int z = m_BestX = len;
    CSeqIterator g = genome + len - m_MaxGaps;
    m_BestVal = 0;
    
    for(int y = -m_MaxGaps; y <= m_MaxGaps && g < gend; ++y, ++g) {
        switch( m_Path[y][z] ) {
        case 'i': case 'd': continue;
        }
        char val = m_Data[y][z];
        if( val > m_BestVal ) { // || (val==m_BestVal && y<=0)) {
            m_BestVal = val;
            m_BestY = y;
        }
    } 
    
    return m_BestVal / max( len, m_BestY + m_BestX - 1);
}

template<class T>
template<class CSeqIterator>
inline void CLcsMatrix<T>::Stat(CSeqIterator genome, CSeqIterator gend,
                                CSeqIterator primer, int length)
{
    int x=m_BestX;
    int y=m_BestY;

    m_Mismatches=length>x?length-x:0;
    m_Matches=m_Insertions=m_Deletions=0;
    
    while(x>0) {
        switch(m_Path[y][x]) {
        default:
        case '=': m_Matches++; x--; break;
        case '!': 
        mismatch: m_Mismatches++; x--; break;
        case 'd': 
            if(m_Path[y-1][x]=='i') goto mismatch; 
            m_Insertions++; y--; 
            break;
        case 'i': 
            if(m_Path[y+1][x-1]=='d') goto mismatch; 
            m_Deletions++; y++; x--; 
            break;
        }
    }
}

template<class T>
template<class CSeqIterator>
inline void CLcsMatrix<T>::Graph(CSeqIterator genome, CSeqIterator gend,
                                 CSeqIterator primer, int length, 
                                 vector<string>& dest, int extra, int intra)
{
    dest.clear();
    dest.resize(3);
    int x=m_BestX;
    int y=m_BestY;

    m_Mismatches=length>x?length-x:0;
    m_Matches=m_Insertions=m_Deletions=0;
    
    for(int pos=y+x-1+extra; extra>0; --pos, --extra) {
        dest[0].push_back(' ');
        dest[1].push_back(genome+pos>=gend?' ':tolower(genome[pos]));
        dest[2].push_back(' ');
    }

    while(x>0) {
        switch(m_Path[y][x]) {
        default:
        case '=': 
            x--; 
            dest[0].push_back(primer[x]);
            dest[1].push_back(genome[x+y]);
            dest[2].push_back('|');
            break;
        case '!': 
        mismatch:
            x--; 
            dest[0].push_back(primer[x]);
            dest[1].push_back(genome[x+y]);
            dest[2].push_back(' ');
            break;
        case 'd': 
            if(m_Path[y-1][x]=='i') goto mismatch; 
            y--; 
            dest[0].push_back('-');
            dest[1].push_back(genome[x+y]);
            dest[2].push_back(' ');
            break;
        case 'i': 
            if(m_Path[y+1][x-1]=='d') goto mismatch; 
            y++; 
            x--; 
            dest[0].push_back(primer[x]);
            dest[1].push_back('-');
            dest[2].push_back(' ');
            break;
        }
    }
    for(int pos=1; pos<=intra; ++pos) {
        dest[0].push_back(' ');
        dest[1].push_back(tolower(genome[-pos]));
        dest[2].push_back(' ');
    }
    if(&genome[0]<&genome[1]) {
        reverse(dest[0].begin(),dest[0].end());
        reverse(dest[1].begin(),dest[1].end());
        reverse(dest[2].begin(),dest[2].end());
    }
}

template<class T> 
inline CLcsMatrix<T>::CLcsMatrix(int length, int maxgaps):
    m_Size(length),
    m_MaxGaps(maxgaps),
    m_Data((new TDataTypePtr[(1+m_MaxGaps)*2+1])+m_MaxGaps+1),
    m_Path((new TDataTypePtr[(1+m_MaxGaps)*2+1])+m_MaxGaps+1)
{
    for(int i=-m_MaxGaps-1; i<=m_MaxGaps+1; ++i) {
        m_Data[i]=new TDataType[m_Size+m_MaxGaps+1];
        memset(m_Data[i],0,(m_Size+m_MaxGaps+1)*sizeof(TDataType));
        m_Path[i]=new TDataType[m_Size+m_MaxGaps+1];
        memset(m_Path[i],0,(m_Size+m_MaxGaps+1)*sizeof(TDataType));
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

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: minilcs.hpp,v $
 * Revision 1.9  2007/07/05 16:23:08  rotmistr
 * Forgot two changes
 *
 * Revision 1.8  2007/07/05 16:05:58  rotmistr
 * Made things compileable by MS Visual C++ 8.0
 *
 * Revision 1.7  2004/09/03 19:06:41  rotmistr
 * Code formatting changes
 *
 * Revision 1.6  2004/07/22 20:40:08  rotmistr
 * Fixed to work with gcc-3.4.0
 *
 * Revision 1.5  2004/06/08 20:32:51  rotmistr
 * Fixup for gap+insert special case
 *
 * Revision 1.4  2004/06/08 16:14:55  rotmistr
 * *** empty log message ***
 *
 * Revision 1.3  2004/06/07 16:24:56  rotmistr
 * Bug fixes to previos version.
 *
 * Revision 1.2  2004/06/03 23:37:20  rotmistr
 * New aligner added.
 *
 * Revision 1.1  2004/06/02 21:37:54  rotmistr
 * Added minilcs function
 *
 */
