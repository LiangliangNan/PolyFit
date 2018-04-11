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

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "spxsolver.h"

namespace soplex
{

void SPxSolver::qualConstraintViolation(Real& maxviol, Real& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   DVector solu( nCols() );

   getPrimal( solu );

   for( int row = 0; row < nRows(); ++row )
   {
      const SVector& rowvec = rowVector( row );

      Real val = 0.0;         

      for( int col = 0; col < rowvec.size(); ++col )
         val += rowvec.value( col ) * solu[rowvec.index( col )];

      Real viol = 0.0;

      assert(lhs( row ) <= rhs( row ) + (100 * epsilon()));

      if (val < lhs( row )) 
         viol = spxAbs(val - lhs( row ));
      else
         if (val > rhs( row ))
            viol = spxAbs(val - rhs( row ));

      if (viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

void SPxSolver::qualBoundViolation(
   Real& maxviol, Real& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   DVector solu( nCols() );

   getPrimal( solu );

   for( int col = 0; col < nCols(); ++col )
   {
      assert(lower( col ) <= upper( col ) + (100 * epsilon()));

      Real viol = 0.0;

      if (solu[col] < lower( col ))
         viol = spxAbs( solu[col] - lower( col ));
      else
         if (solu[col] > upper( col ))
            viol = spxAbs( solu[col] - upper( col ));
         
      if (viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

void SPxSolver::qualSlackViolation(Real& maxviol, Real& sumviol) const
{
   maxviol = 0.0;
   sumviol = 0.0;

   DVector solu( nCols() );
   DVector slacks( nRows() );

   getPrimal( solu );
   getSlacks( slacks );

   for( int row = 0; row < nRows(); ++row )
   {
      const SVector& rowvec = rowVector( row );

      Real val = 0.0;         

      for( int col = 0; col < rowvec.size(); ++col )
         val += rowvec.value( col ) * solu[rowvec.index( col )];

      Real viol = spxAbs(val - slacks[row]);

      if (viol > maxviol)
         maxviol = viol;

      sumviol += viol;
   }
}

void SPxSolver::qualRedCostViolation(Real& maxviol, Real& sumviol) const
{   
   maxviol = 0.0;
   sumviol = 0.0;

   int i;
   // TODO:   y = c_B * B^-1  => coSolve(y, c_B)
   //         redcost = c_N - yA_N 
   // solve system "x = e_i^T * B^-1" to get i'th row of B^-1
   // DVector y( nRows() );
   // basis().coSolve( x, spx->unitVector( i ) );
   // DVector rdcost( nCols() );
#if 0 // un-const
   if (lastUpdate() > 0)
      factorize();

   computePvec();

   if (type() == ENTER)
      computeTest();
#endif
   if (type() == ENTER)
   {
      for(i = 0; i < dim(); ++i)
      {
         Real x = coTest()[i];
         
         if (x < 0.0)
         {
            sumviol -= x;
            
            if (x < maxviol)
               maxviol = x;
         }
      }
      for(i = 0; i < coDim(); ++i)
      {
         Real x = test()[i];
         
         if (x < 0.0)
         {
            sumviol -= x;
            
            if (x < maxviol)
               maxviol = x;
         }
      } 
   }
   else
   {
      assert(type() == LEAVE);

      for(i = 0; i < dim(); ++i)
      {
         Real x = fTest()[i];
         
         if (x < 0.0)
         {
            sumviol -= x;
            
            if (x < maxviol)
               maxviol = x;
         }
      }
   }
   maxviol *= -1;
}

} // namespace soplex
