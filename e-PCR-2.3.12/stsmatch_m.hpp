/* $Id: stsmatch_m.hpp,v 1.7 2008/06/18 14:45:33 rotmistr Exp $
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

#ifndef EPCR_STSMATCH_M__HPP
#define EPCR_STSMATCH_M__HPP

#include <epcr/build_cfg.h>
#include <epcr/stsmatch_i.hpp>
#include <stdexcept>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(EPCR_SCOPE)

class CMmFileSts;
class CStsFileHash;
class IStsFileCallback;

class CMmFileSts:public ISts
{
public:
    enum EFlags {
        fNone = 0,
        fReverse = 0x01,
        fAllocLeft = 0x10,
        fAllocRight = 0x20
    };

    virtual const char * GetPrimerData(int s) const 
        { return m_Primer[s]; }
    virtual unsigned     GetPrimerLength(int s) const 
        { return m_Length[s]; }
    virtual EDirect      GetDirection() const 
        { return m_Flags&fReverse?eMinus:ePlus; }
    virtual unsigned     GetSizeLo() const { return m_SizeLo; }
    virtual unsigned     GetSizeHi() const { return m_SizeHi; }
    virtual CStrRef      GetName() const;
    virtual CStrRef      GetDescription() const;
    virtual int GetOverhangChars(int s) const { return m_OvhgChars[s]; }
			
    virtual ~CMmFileSts() throw () { 
        if(m_Flags&fAllocLeft)  delete[] m_Primer[0]; 
        if(m_Flags&fAllocRight) delete[] m_Primer[1]; 
    }

    CMmFileSts(const char * ref, 
               const CStrRef& p1, 
               const CStrRef& p2,
               unsigned lo, unsigned hi, unsigned flags,
               unsigned char ovhg1 = 0, unsigned char ovhg2 = 0);

    static int Parse(const char * ref, CStrRef * dest, int maxf=5);
    static void ParseRange(const CStrRef& range, int& lo, int& hi, 
						   int defLo=ePCR_DEFAULT_size_lo,
						   int defHi=ePCR_DEFAULT_size_hi);

    bool Valid() const { return m_Reference!=0; }

protected:
    const char *   m_Primer[2];
    char           m_Length[2];
    unsigned       m_SizeLo;
    unsigned       m_SizeHi;
    const char *   m_Reference;
    char           m_Flags;
    char           m_OvhgChars[2];
};

// This class implements loading of STS hash table from mmapped file
class CStsFileHash:public CStsHash
{
public:
    enum EFlags {
        fAllowOverhang = 1,
        fUnmaskPrimers = 2,
        fNONE = 0
    };
    virtual void Reset() { 
		CStsHash::Reset(); DetachFile(); 
		//SetDefaultSize(ePCR_DEFAULT_size_lo,ePCR_DEFAULT_size_hi); 
	}
    virtual ~CStsFileHash() throw () {}
	CStsFileHash():m_MemoryBase(0),m_MemorySize(0),
				   m_DefaultSizeLo(ePCR_DEFAULT_size_lo),
                   m_DefaultSizeHi(ePCR_DEFAULT_size_hi),m_Flags(0) {}
	
	unsigned GetDefaultSizeLo() const { return m_DefaultSizeLo; }
	unsigned GetDefaultSizeHi() const { return m_DefaultSizeHi; }
	
	void SetDefaultSize(unsigned lo, unsigned hi) {
		//cerr << "\e[31m" << __PRETTY_FUNCTION__ << "\e[32m: lo = " << lo << ", hi = " << hi << "\e[0m\n";
		m_DefaultSizeLo=lo;
		m_DefaultSizeHi=hi;
	}
	
    virtual void ReadStsFile(const string& fname, IStsFileCallback* cbk = 0);
    
    bool AllowOverhang() const { return m_Flags & fAllowOverhang; }
    bool UnmaskPrimers() const { return m_Flags & fUnmaskPrimers; }

    void SetFlags(int flags, bool on=true) { 
        if(on) m_Flags|=flags; else m_Flags&=~flags;
    }

protected:
    virtual void AttachFile(const string& fname);
    virtual void DetachFile();
    bool ParseLine(IStsFileCallback * cbk, const char * pos, unsigned len);
    
protected:
    char * m_MemoryBase;
    unsigned     m_MemorySize;
	unsigned     m_DefaultSizeLo;
	unsigned     m_DefaultSizeHi;
	unsigned     m_Flags;
};

// Callback for warnings and errors
class IStsFileCallback
{
public:
    enum EError {
        eErrOK,
        eErrBadLine,
        eErrShortPrimer,
        eErrAmbiquosPrimer,
        eErrSystem,
        eErrTOTAL
    };
				
    virtual ~IStsFileCallback() throw () {}
    virtual bool NextLine(const char * line, int sz) =0;
    virtual bool Error(EError) =0;
    virtual bool Start() =0;
    virtual bool Done() =0;
};

END_SCOPE(EPCR_SCOPE)
END_NCBI_SCOPE

#endif

/*
 * $Log: stsmatch_m.hpp,v $
 * Revision 1.7  2008/06/18 14:45:33  rotmistr
 * Fixed problem with -d x-X parameter being reset if -w or some others are used after it.
 *
 * Revision 1.6  2004/10/26 17:16:36  rotmistr
 * Added 5'-end masking for primers
 *
 * Revision 1.5  2004/04/01 05:57:53  rotmistr
 * Compilable with borland C++
 *
 * Revision 1.4  2004/03/30 21:06:53  rotmistr
 * Fixes for setting default STS size range.
 *
 * Revision 1.3  2004/03/30 19:08:03  rotmistr
 * default STS size is tunnable now
 *
 * Revision 1.2  2004/01/28 23:27:02  rotmistr
 * "Best of overlapping" hit selection postprocessor added.
 *
 * Revision 1.1.1.1  2003/12/23 18:17:27  rotmistr
 * Package that includes e-PCR, reverse e-PCR, and sequence data preparation
 * program for reverse e-PCR looks ready
 *
 */
