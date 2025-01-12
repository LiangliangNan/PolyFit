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

/**@file  spxequilisc.cpp
 * @brief Equilibrium row/column scaling.
 */
#include <assert.h>

#include "spxequilisc.h"
#include "spxout.h"
#include "spxlpbase.h"
#include "spxlp.h"
#include "soplex.h"

namespace soplex
{
static const char* makename(bool doBoth)
{
   return doBoth ? "bi-Equilibrium" : "uni-Equilibrium";
}

/// maximum ratio between absolute biggest and smallest element in any scaled row/column.
static Real maxPrescaledRatio(const SPxLP& lp, const std::vector<Real>& coScaleval, bool rowRatio)
{
   Real pmax = 0.0;
   const int n = rowRatio ? lp.nRows() : lp.nCols();

   for( int i = 0; i < n; ++i )
   {
      const SVector& vec = rowRatio ? lp.rowVector(i) : lp.colVector(i);
      Real mini = infinity;
      Real maxi = 0.0;

      for( int j = 0; j < vec.size(); ++j )
      {
         assert(vec.index(j) >= 0);
         const Real x = spxAbs(vec.value(j)) * coScaleval[unsigned(vec.index(j))];

         if( isZero(x) )
            continue;
         if( x < mini )
            mini = x;
         if( x > maxi )
            maxi = x;
      }

      if( mini == infinity )
         continue;

      const Real p = maxi / mini;

      if( p > pmax )
         pmax = p;
   }
   return pmax;
}


void SPxEquiliSC::computeEquiExpVec(const SVSet* vecset, const DataArray<int>& coScaleExp, DataArray<int>& scaleExp)
{
   assert(vecset != nullptr);

   for( int i = 0; i < vecset->num(); ++i )
   {
      const SVector& vec = (*vecset)[i];

      Real maxi = 0.0;

      for( int j = 0; j < vec.size(); ++j )
      {
         const Real x = spxAbs(spxLdexp(vec.value(j), coScaleExp[vec.index(j)]));

         if( GT(x, maxi) )
            maxi = x;
      }
      // empty rows/cols are possible
      if( maxi == 0.0 )
         maxi = 1.0;

      assert(maxi > 0.0);

      spxFrexp(1.0 / maxi, &(scaleExp[i]));

      scaleExp[i] -= 1;
   }
}

void SPxEquiliSC::computeEquiExpVec(const SVSet* vecset, const std::vector<Real>& coScaleVal, DataArray<int>& scaleExp)
{
   assert(vecset != nullptr);

   for( int i = 0; i < vecset->num(); ++i )
   {
      const SVector& vec = (*vecset)[i];

      Real maxi = 0.0;

      for( int j = 0; j < vec.size(); ++j )
      {
         assert(vec.index(j) >= 0);
         const Real x = spxAbs(vec.value(j) * coScaleVal[unsigned(vec.index(j))]);

         if( GT(x, maxi) )
            maxi = x;
      }
      // empty rows/cols are possible
      if( maxi == 0.0 )
         maxi = 1.0;

      assert(maxi > 0.0);

      spxFrexp(1.0 / maxi, &(scaleExp[i]));

      scaleExp[i] -= 1;
   }
}

void SPxEquiliSC::computePostequiExpVecs(const SPxLPBase<Real>& lp, const std::vector<Real>& preRowscale, const std::vector<Real>& preColscale,
      DataArray<int>& rowscaleExp, DataArray<int>& colscaleExp)
{
   const Real colratio = maxPrescaledRatio(lp, preRowscale, false);
   const Real rowratio = maxPrescaledRatio(lp, preColscale, true);

   const bool colFirst = colratio < rowratio;

   // see SPxEquiliSC::scale for reason behind this branch
   if( colFirst )
   {
      computeEquiExpVec(lp.colSet(), preRowscale, colscaleExp);
      computeEquiExpVec(lp.rowSet(), colscaleExp, rowscaleExp);
   }
   else
   {
      computeEquiExpVec(lp.rowSet(), preColscale, rowscaleExp);
      computeEquiExpVec(lp.colSet(), rowscaleExp, colscaleExp);
   }
}


SPxEquiliSC::SPxEquiliSC(bool doBoth)
   : SPxScaler(makename(doBoth), false, doBoth)
{}

SPxEquiliSC::SPxEquiliSC(const SPxEquiliSC& old)
   : SPxScaler(old)
{}

SPxEquiliSC& SPxEquiliSC::operator=(const SPxEquiliSC& rhs)
{
   if(this != &rhs)
   {
      SPxScaler::operator=(rhs);
   }

   return *this;
}


void SPxEquiliSC::scale(SPxLP& lp, bool persistent)
{

   MSG_INFO1( (*spxout), (*spxout) << "Equilibrium scaling LP" << (persistent ? " (persistent)" : "") << std::endl; )

   setup(lp);

   /* We want to do the direction first, which has a lower maximal ratio,
    * since the lowest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first
    * Example:
    *
    *                     Rowratio
    *            0.1  1   10
    *            10   1   10
    *
    * Colratio   100  1
    *
    * Row first =>         Col next =>
    *            0.1  1          0.1  1
    *            1    0.1        1    0.1
    *
    * Col first =>         Row next =>
    *            0.01 1          0.01 1
    *            1    1          1    1
    *
    */
   Real colratio = maxColRatio(lp);
   Real rowratio = maxRowRatio(lp);

   bool colFirst = colratio < rowratio;

   MSG_INFO2( (*spxout), (*spxout) << "before scaling:"
                        << " min= " << lp.minAbsNzo()
                        << " max= " << lp.maxAbsNzo()
                        << " col-ratio= " << colratio
                        << " row-ratio= " << rowratio
                        << std::endl; )

   if (colFirst)
   {
      computeEquiExpVec(lp.colSet(), *m_activeRowscaleExp, *m_activeColscaleExp);

      if (m_doBoth)
         computeEquiExpVec(lp.rowSet(), *m_activeColscaleExp, *m_activeRowscaleExp);
   }
   else
   {
      computeEquiExpVec(lp.rowSet(), *m_activeColscaleExp, *m_activeRowscaleExp);

      if (m_doBoth)
         computeEquiExpVec(lp.colSet(), *m_activeRowscaleExp, *m_activeColscaleExp);
   }

   /* scale */
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
