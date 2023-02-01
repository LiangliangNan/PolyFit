/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxleastsqsc.cpp
 * @brief LP least squares scaling.
 */
#include <cmath>
#include <assert.h>
#include "spxleastsqsc.h"
#include "spxout.h"
#include "basevectors.h"
#include "svsetbase.h"
#include "svectorbase.h"
#include "ssvectorbase.h"
#include <array>

namespace soplex
{

/* update scaling vector */
static void updateScale(
   const SSVector vecnnzeroes,
   const SSVector resnvec,
   SSVector& tmpvec,
   SSVector*& psccurr,
   SSVector*& pscprev,
   Real qcurr,
   Real qprev,
   Real eprev1,
   Real eprev2)
{
   assert(psccurr != NULL);
   assert(pscprev != NULL);
   assert(qcurr * qprev != 0.0);

   Real fac = -(eprev1 * eprev2);

   SSVector* pssv;

   *pscprev -= *psccurr;

   if( isZero(fac) )
      (*pscprev).clear();
   else
      *pscprev *= fac;

   *pscprev += tmpvec.assignPWproduct4setup(resnvec, vecnnzeroes);

   *pscprev *= 1.0 / (qcurr * qprev);
   *pscprev += *psccurr;

   /* swap pointers */
   pssv = psccurr;
   psccurr = pscprev;
   pscprev = pssv;
}

/* update scaling vector after main loop */
static void updateScaleFinal(
   const SSVector vecnnzeroes,
   const SSVector resnvec,
   SSVector& tmpvec,
   SSVector*& psccurr,
   SSVector*& pscprev,
   Real q,
   Real eprev1,
   Real eprev2)
{
   assert(q != 0);
   assert(psccurr != NULL);
   assert(pscprev != NULL);

   Real fac = -(eprev1 * eprev2);

   *pscprev -= *psccurr;

   if( isZero(fac) )
      (*pscprev).clear();
   else
      *pscprev *= fac;

   *pscprev += tmpvec.assignPWproduct4setup(resnvec, vecnnzeroes);
   *pscprev *= 1.0 / q;
   *pscprev += *psccurr;

   psccurr = pscprev;
}

/* update residual vector */
static inline void updateRes(
   const SVSet facset,
   const SSVector resvecprev,
   SSVector& resvec,
   SSVector& tmpvec,
   Real eprev,
   Real qcurr)
{
   assert(qcurr != 0.0);

   if( isZero(eprev) )
      resvec.clear();
   else
      resvec *= eprev;

   tmpvec.assign2product4setup(facset, resvecprev);
   tmpvec.setup();
   resvec += tmpvec;

   resvec *= (-1.0 / qcurr);
   resvec.setup();
}


/* initialize constant vectors and matrices */
static void initConstVecs(
   const SVSet* vecset,
   SVSet& facset,
   SSVector& veclogs,
   SSVector& vecnnzinv)
{
   assert(vecset != NULL);

   const int nvec = vecset->num();

   for( int k = 0; k < nvec; ++k )
   {
      Real logsum = 0.0;
      int nnz = 0;
      // get kth row or column of LP
      const SVector& lpvec = (*vecset)[k];
      const int size = lpvec.size();

      for( int i = 0; i < size; ++i )
      {
         const Real a = lpvec.value(i);

         if( !isZero(a) )
         {
            logsum += log2(double(spxAbs(a))); // todo spxLog2?
            nnz++;
         }
      }

      Real nnzinv;
      if( nnz > 0)
      {
         nnzinv = 1.0 / nnz;
      }
      else
      {
         /* all-0 entries */
         logsum = 1.0;
         nnzinv = 1.0;
      }

      veclogs.add(k, logsum);
      vecnnzinv.add(k, nnzinv);

      /* create new vector for facset */
      SVector& vecnew = (*(facset.create(nnz)));

      for( int i = 0; i < size; ++i )
      {
         if( !isZero(lpvec.value(i)) )
            vecnew.add(lpvec.index(i), nnzinv);
      }
      vecnew.sort();
   }

   assert(veclogs.isSetup());
   assert(vecnnzinv.isSetup());
}

/* return name of scaler */
static const char* makename()
{
   return "Least squares";
}

SPxLeastSqSC::SPxLeastSqSC()
   : SPxScaler(makename(), false, false)
{}

SPxLeastSqSC::SPxLeastSqSC(const SPxLeastSqSC& old)
   : SPxScaler(old), acrcydivisor(old.acrcydivisor), maxrounds(old.maxrounds)
{}

SPxLeastSqSC& SPxLeastSqSC::operator=(const SPxLeastSqSC& rhs)
{
   if(this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}


void SPxLeastSqSC::setRealParam(Real param, const char* name)
{
   assert(param >= 1.0);
   acrcydivisor = param;
}

void SPxLeastSqSC::setIntParam(int param, const char* name)
{
   assert(param >= 0);
   maxrounds = param;
}

void SPxLeastSqSC::scale(SPxLP& lp,  bool persistent)
{
   MSG_INFO1( (*spxout), (*spxout) << "Least squares LP scaling" << (persistent ? " (persistent)" : "") << std::endl; )

   setup(lp);

   const int nrows = lp.nRows();
   const int ncols = lp.nCols();
   const int lpnnz = lp.nNzos();

   /* constant factor matrices;
    * in Curtis-Reid article
    * facnrows equals E^T M^(-1)
    * facncols equals E N^(-1)
    * */
   SVSet facnrows(nrows, nrows, 1.1, 1.2);
   SVSet facncols(ncols, ncols, 1.1, 1.2);

   /* column scaling factor vectors */
   SSVector colscale1(ncols);
   SSVector colscale2(ncols);

   /* row scaling factor vectors */
   SSVector rowscale1(nrows);
   SSVector rowscale2(nrows);

   /* residual vectors */
   SSVector resnrows(nrows);
   SSVector resncols(ncols);

   /* vectors to store temporary values */
   SSVector tmprows(nrows);
   SSVector tmpcols(ncols);

   /* vectors storing the row and column sums (respectively) of logarithms of
    *(absolute values of) non-zero elements of left hand matrix of LP
    */
   SSVector rowlogs(nrows);
   SSVector collogs(ncols);

   /* vectors storing the inverted number of non-zeros in each row and column
    *(respectively) of left hand matrix of LP
    */
   SSVector rownnzinv(nrows);
   SSVector colnnzinv(ncols);

   /* vector pointers */
   SSVector* csccurr = &colscale1;
   SSVector* cscprev = &colscale2;
   SSVector* rsccurr = &rowscale1;
   SSVector* rscprev = &rowscale2;

   MSG_INFO2( (*spxout), (*spxout) << "before scaling:"
      << " min= " << lp.minAbsNzo()
      << " max= " << lp.maxAbsNzo()
      << " col-ratio= " << maxColRatio(lp)
      << " row-ratio= " << maxRowRatio(lp)
      << std::endl; )

   /* initialize scalars, vectors and matrices */

   assert(acrcydivisor > 0.0);

   const Real smax = lpnnz / acrcydivisor;
   Real qcurr = 1.0;
   Real qprev = 0.0;

   std::array<Real, 3> eprev;
   eprev.fill(0.0);

   initConstVecs(lp.rowSet(), facnrows, rowlogs, rownnzinv);
   initConstVecs(lp.colSet(), facncols, collogs, colnnzinv);

   assert(tmprows.isSetup());
   assert(tmpcols.isSetup());
   assert(rowscale1.isSetup());
   assert(rowscale2.isSetup());
   assert(colscale1.isSetup());
   assert(colscale2.isSetup());

   // compute first residual vector r0
   resncols = collogs - tmpcols.assign2product4setup(facnrows, rowlogs);

   resncols.setup();
   resnrows.setup();

   rowscale1.assignPWproduct4setup(rownnzinv, rowlogs);
   rowscale2 = rowscale1;

   Real scurr = resncols * tmpcols.assignPWproduct4setup(colnnzinv, resncols);

   int k;

   /* conjugate gradient loop */
   for( k = 0; k < maxrounds; ++k )
   {
      const Real sprev = scurr;

      // is k even?
      if( (k % 2) == 0 )
      {
         // not in first iteration?
         if( k != 0 ) // true, then update row scaling factor vector
            updateScale(rownnzinv, resnrows, tmprows, rsccurr, rscprev, qcurr, qprev, eprev[1], eprev[2]);

         updateRes(facncols, resncols, resnrows, tmprows, eprev[0], qcurr);
         scurr = resnrows * tmprows.assignPWproduct4setup(resnrows, rownnzinv);
      }
      else // k is odd
      {
         // update column scaling factor vector
         updateScale(colnnzinv, resncols, tmpcols, csccurr, cscprev, qcurr, qprev, eprev[1], eprev[2]);

         updateRes(facnrows, resnrows, resncols, tmpcols, eprev[0], qcurr);
         scurr = resncols * tmpcols.assignPWproduct4setup(resncols, colnnzinv);
      }

      // shift eprev entries one to the right
      for( unsigned l = 2; l > 0; --l)
         eprev[l] = eprev[l - 1];

      eprev[0] = (qcurr * scurr) / sprev;

      const Real tmp = qcurr;
      qcurr = 1.0 - eprev[0];
      qprev = tmp;

      // termination criterion met?
      if( scurr < smax )
         break;
   }

   // is k even?
   if( (k % 2) == 0 )
   {
      // update column scaling factor vector
      updateScaleFinal(colnnzinv, resncols, tmpcols, csccurr, cscprev, qprev, eprev[1], eprev[2]);
   }
   else // k is odd
   {
      // update row scaling factor vector
      updateScaleFinal(rownnzinv, resnrows, tmprows, rsccurr, rscprev, qprev, eprev[1], eprev[2]);
   }

   /* compute actual scaling factors */

   const SSVector& rowscale = *rsccurr;
   const SSVector& colscale = *csccurr;

   DataArray<int>& colscaleExp = *m_activeColscaleExp;
   DataArray<int>& rowscaleExp = *m_activeRowscaleExp;

   for( k = 0; k < nrows; ++k )
      rowscaleExp[k] = -int( rowscale[k] + ((rowscale[k] >= 0.0)? (+0.5) : (-0.5)) );

   for( k = 0; k < ncols; ++k )
      colscaleExp[k] = -int( colscale[k] + ((colscale[k] >= 0.0)? (+0.5) : (-0.5)) );

   // scale
   applyScaling(lp);

   MSG_INFO3( (*spxout), (*spxout) << "Row scaling min= " << minAbsRowscale()
      << " max= " << maxAbsRowscale()
      << std::endl
      << "Col scaling min= " << minAbsColscale()
      << " max= " << maxAbsColscale()
      << std::endl; )

   MSG_INFO2( (*spxout), (*spxout) << "after scaling: "
      << " min= " << lp.minAbsNzo(false)
      << " max= " << lp.maxAbsNzo(false)
      << " col-ratio= " << maxColRatio(lp)
      << " row-ratio= " << maxRowRatio(lp)
      << std::endl; )
}

} // namespace soplex
