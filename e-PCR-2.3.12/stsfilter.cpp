/* $Id: stsfilter.cpp,v 1.4 2004/06/03 23:37:22 rotmistr Exp $
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

#include <epcr/hashset.hpp>
#include <epcr/stsmatch_i.hpp>

#include <stdexcept>

USING_NCBI_SCOPE;
USING_SCOPE(EPCR_SCOPE);

#define DBG      fprintf(stdout,__FILE__"[%d]\n",__LINE__);
#define SHOWI(a) fprintf(stdout,__FILE__"[%d] "#a"=%d\n",__LINE__,a);
#define SHOWS(a) \
fprintf(stdout,__FILE__"[%d] "#a"=%.*s\n",__LINE__,a.length(),a.data());


bool CPcrMachinePostprocess::OrderByPos2Pos1(const SOutput& o1,
                                             const SOutput& o2) 
{
    return (o1.pos2<o2.pos2 || o1.pos2==o2.pos2 && o1.pos1<o2.pos1);
}

int CPcrMachinePostprocess::Compare(const SOutput& a, const SOutput& b, 
                                    int min_len, int max_len) 
{
    if(a.gaps()>b.gaps()) return -1;
    if(a.gaps()<b.gaps()) return +1;
    if(a.mism()>b.mism()) return -1;
    if(a.mism()<b.mism()) return +1;
    int lb=(b.pos2-b.pos1);
    int la=(a.pos2-a.pos1);
    if(la >= min_len && la <= max_len) {
        if(lb >= min_len && lb <= max_len) {
            return lb-la;
        } else return 1;
    } else {
        if(lb >= min_len && lb <= max_len) {
            return -1;
        } else {
            int da=(la < min_len)?min_len-la:la-max_len;
            int db=(lb < min_len)?min_len-lb:lb-max_len;
            return db-da;
        }
    }
}

bool CPcrMachinePostprocess::Overlap(int l1, int l2, 
                                     const SOutput& a, const SOutput& b)
{
    return (abs(a.pos2-b.pos2)<l2 && abs(a.pos1-b.pos1)<l1);
}

void CPcrMachinePostprocess::Flush()
{
    for(TAllHits_I i=m_OutQueues.begin(); i!=m_OutQueues.end(); ++i) {
        // 1st: cluster 
        // 2nd: select best per cluster
        if(i->second.size()==0) continue;

        int l1=i->first->GetPrimerLength(0);
        int l2=i->first->GetPrimerLength(1);
        int lo=i->first->GetSizeLo();
        int hi=i->first->GetSizeHi();

        TStsHits& hits=i->second;

        while(hits.size()) {
            
            TStsHits todo;
            TStsHits cluster;
            cluster.push_back(hits.back());
            SOutput best=cluster.front();

            hits.pop_back();

            for(TStsHits_CI j=hits.begin(); j!=hits.end(); ++j) {
                bool overlaps=false;
                for(TStsHits_CI h=cluster.begin(); h!=cluster.end(); ++h)
                    if(Overlap(l1,l2,*j,*h)) { overlaps=true; break; }
                if(overlaps) {
                    if(Compare(*j,best,lo,hi)>=0) best=*j;
                    cluster.push_back(*j);
                } else {
                    todo.push_back(*j);
                }
            }

            SScore sc(best.length(),best.mism_l, best.mism_r, 
					  best.gaps_l, best.gaps_r);
            m_Callback->CbkMatch(i->first, best.pos1, best.pos2, &sc);

            hits=todo;
        }
    }
    m_OutQueues.clear();
}

void CPcrMachinePostprocess::CbkMatch(const ISts * sts,
                                      unsigned pos1, unsigned pos2,
                                      const SScore* score)
{
    SOutput o(pos1,pos2,
			  score->mism_l,score->mism_r,
			  score->gaps_l,score->gaps_r);
    TStsHits& hits=m_OutQueues[sts];
	if(hits.size() && hits.back()==o) return;
    hits.push_back(o);
}


/*
 * $Log: stsfilter.cpp,v $
 * Revision 1.4  2004/06/03 23:37:22  rotmistr
 * New aligner added.
 *
 * Revision 1.3  2004/03/25 19:36:52  rotmistr
 * API: separate left and right primers mism/gaps in forward API
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
 * Revision 1.1  2004/01/28 23:27:02  rotmistr
 * "Best of overlapping" hit selection postprocessor added.
 *
 */
