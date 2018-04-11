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

#include <iostream>

#include "spxdefines.h"
#include "spxsumst.h"
#include "vector.h"

namespace soplex
{

void SPxSumST::setupWeights(SPxSolver& base)
{
   int count;
   int i;
   Real x;
   DVector work, delta, rowLen;

   assert(base.nRows() > 0);
   assert(base.nCols() > 0);

   rowLen.reDim(base.nRows(), true);
   work.reDim(base.nCols(), true);
   delta.reDim(base.nCols(), true);

   Real* wrk = work.get_ptr();
   const Real* lhs = base.lhs().get_const_ptr();
   const Real* rhs = base.rhs().get_const_ptr();
   const Real* up = base.upper().get_const_ptr();
   const Real* low = base.lower().get_const_ptr();

   for (i = base.nRows(); --i >= 0;)
   {
      rowLen[i] = base.rowVector(i).length2();
      if (lhs[i] > 0)
         delta.multAdd(lhs[i] / rowLen[i], base.rowVector(i));
      else if (rhs[i] < 0)
         delta.multAdd(rhs[i] / rowLen[i], base.rowVector(i));
   }

   for (count = 0;; count++)
   {
      work += delta;
      for (i = base.nCols(); --i >= 0;)
      {
         if (wrk[i] > up[i])
            wrk[i] = up[i];
         if (wrk[i] < low[i])
            wrk[i] = low[i];
      }

      //      std::cout << -(work * base.maxObj()) << std::endl;
      if (count >= 12)
         break;

      delta.clear();
      for (i = base.nRows(); --i >= 0;)
      {
         x = base.rowVector(i) * work;
         if (lhs[i] > x)
            delta.multAdd((lhs[i] - x) / rowLen[i], base.rowVector(i));
         else if (rhs[i] < x)
            delta.multAdd((rhs[i] - x) / rowLen[i], base.rowVector(i));
      }
   }

   primal(work);
   SPxVectorST::setupWeights(base);
}
} // namespace soplex
