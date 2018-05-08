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

/**@file  spxgeometsc.cpp
 * @brief Geometric mean row/column scaling.
 */
#include <assert.h>

#include "spxgeometsc.h"
#include "spxout.h"
#include "spxlpbase.h"
#include "spxequilisc.h"

namespace soplex
{

static Real computeScalingVec(
      const SVSet*             vecset,
      const std::vector<Real>& coScaleval,
      std::vector<Real>&       scaleval)
{
   Real pmax = 0.0;

   assert(scaleval.size() == unsigned(vecset->num()));

   for( int i = 0; i < vecset->num(); ++i )
   {
      const SVector& vec = (*vecset)[i];

      Real maxi = 0.0;
      Real mini = infinity;

      for( int j = 0; j < vec.size(); ++j )
      {
         const Real x = spxAbs(vec.value(j) * coScaleval[unsigned(vec.index(j))]);

         if (!isZero(x))
         {
            if (x > maxi)
               maxi = x;
            if (x < mini)
               mini = x;
         }
      }
      // empty rows/cols are possible
      if (mini == infinity || maxi == 0.0)
      {
         mini = 1.0;
         maxi = 1.0;
      }
      assert(mini < infinity);
      assert(maxi > 0.0);

      scaleval[unsigned(i)] = 1.0 / spxSqrt(mini * maxi);

      const Real p = maxi / mini;

      if (p > pmax)
         pmax = p;
   }
   return pmax;
}


SPxGeometSC::SPxGeometSC(bool equilibrate, int maxIters, Real minImpr, Real goodEnough)
   : SPxScaler("Geometric")
   , postequilibration(equilibrate)
   , m_maxIterations(maxIters)
   , m_minImprovement(minImpr)
   , m_goodEnoughRatio(goodEnough)
{
   assert(maxIters > 0);
   assert(minImpr > 0.0 && minImpr <= 1.0);
   assert(goodEnough >= 0.0);
}

SPxGeometSC::SPxGeometSC(const SPxGeometSC& old)
   : SPxScaler(old)
   , postequilibration(old.postequilibration)
   , m_maxIterations(old.m_maxIterations)
   , m_minImprovement(old.m_minImprovement)
   , m_goodEnoughRatio(old.m_goodEnoughRatio)
{
   assert(m_maxIterations > 0);
   assert(m_minImprovement > 0.0 && m_minImprovement <= 1.0);
   assert(m_goodEnoughRatio >= 0.0);
}

SPxGeometSC& SPxGeometSC::operator=(const SPxGeometSC& rhs)
{
   if (this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}

void SPxGeometSC::scale(SPxLPBase<Real>& lp, bool persistent)
{

   MSG_INFO1( (*spxout), (*spxout) << "Geometric scaling LP" << (persistent ? " (persistent)" : "") << (postequilibration ? " with post-equilibration" : "") << std::endl; )

   setup(lp);

   /* We want to do that direction first, with the lower ratio.
    * See SPxEquiliSC::scale() for a reasoning.
    */
   const Real colratio = maxColRatio(lp);
   const Real rowratio = maxRowRatio(lp);

   const bool colFirst = colratio < rowratio;

   Real p0start;
   Real p1start;

   if( colFirst )
   {
     p0start = colratio;
     p1start = rowratio;
   }
   else
   {
     p0start = rowratio;
     p1start = colratio;
   }

   MSG_INFO2( (*spxout), (*spxout) << "before scaling:"
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio
                        << " row-ratio= " << rowratio
                        << std::endl; )

   // perform geometric scaling only if maximum ratio is above threshold
   bool geoscale = p1start > m_goodEnoughRatio;

   if( !geoscale )
   {
      MSG_INFO2( (*spxout), (*spxout) << "No geometric scaling done, ratio good enough" << std::endl; )

      if( !postequilibration )
      {
         lp.setScalingInfo(true);
         return;
      }

      MSG_INFO2( (*spxout), (*spxout) << " ... but will still perform equilibrium scaling" << std::endl; )
   }

   std::vector<Real> rowscale(unsigned(lp.nRows()), 1.0);
   std::vector<Real> colscale(unsigned(lp.nCols()), 1.0);

   Real p0 = 0.0;
   Real p1 = 0.0;

   if( geoscale )
   {
      Real p0prev = p0start;
      Real p1prev = p1start;

      // we make at most maxIterations.
      for( int count = 0; count < m_maxIterations; count++ )
      {
         if( colFirst )
         {
            p0 = computeScalingVec(lp.colSet(), rowscale, colscale);
            p1 = computeScalingVec(lp.rowSet(), colscale, rowscale);
         }
         else
         {
            p0 = computeScalingVec(lp.rowSet(), colscale, rowscale);
            p1 = computeScalingVec(lp.colSet(), rowscale, colscale);
         }

         MSG_INFO3( (*spxout), (*spxout) << "Geometric scaling round " << count
                              << " col-ratio= " << (colFirst ? p0 : p1)
                              << " row-ratio= " << (colFirst ? p1 : p0)
                              << std::endl; )

         if( p0 > m_minImprovement * p0prev && p1 > m_minImprovement * p1prev )
            break;

         p0prev = p0;
         p1prev = p1;
      }

      // perform geometric scaling only if there is enough (default 15%) improvement.
      geoscale = (p0 <= m_minImprovement * p0start || p1 <= m_minImprovement * p1start);
   }

   if( !geoscale && !postequilibration )
   {
      MSG_INFO2( (*spxout), (*spxout) << "No geometric scaling done." << std::endl; )
      lp.setScalingInfo(true);
   }
   else
   {
      DataArray<int>& colscaleExp = *m_activeColscaleExp;
      DataArray<int>& rowscaleExp = *m_activeRowscaleExp;

      if( postequilibration )
      {
         if( !geoscale )
         {
            std::fill(rowscale.begin(), rowscale.end(), 1.0);
            std::fill(colscale.begin(), colscale.end(), 1.0);
         }
         SPxEquiliSC::computePostequiExpVecs(lp, rowscale, colscale, rowscaleExp, colscaleExp);
      }
      else
      {
         computeExpVec(colscale, colscaleExp);
         computeExpVec(rowscale, rowscaleExp);
      }

      applyScaling(lp);

      MSG_INFO3( (*spxout), (*spxout) << "Row scaling min= " << minAbsRowscale()
                           << " max= " << maxAbsRowscale()
                           << std::endl
                           << "IGEOSC06 Col scaling min= " << minAbsColscale()
                           << " max= " << maxAbsColscale()
                           << std::endl; )

      MSG_INFO2( (*spxout), (*spxout) << "after scaling: "
                           << " min= " << lp.minAbsNzo(false)
                           << " max= " << lp.maxAbsNzo(false)
                           << " col-ratio= " << maxColRatio(lp) 
                           << " row-ratio= " << maxRowRatio(lp) 
                           << std::endl; )
   }
}


} // namespace soplex
