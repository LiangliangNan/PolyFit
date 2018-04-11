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
#include "spxout.h"

namespace soplex
{
void SPxSolver::shiftFvec()
{

   /* the allowed tolerance is (rep() == COLUMN) ? feastol() : opttol() because theFvec is the primal vector in COLUMN
    * and the dual vector in ROW representation; this is equivalent to entertol()
    */
   Real minrandom = 10.0 * entertol();
   Real maxrandom = 100.0 * entertol();
   Real allow = entertol() - epsilon();

   assert(type() == ENTER);
   assert(allow > 0);

   for (int i = dim() - 1; i >= 0; --i)
   {
      if (theUBbound[i] + allow < (*theFvec)[i])
      {
         MSG_DEBUG( std::cout << "DSHIFT08 theUBbound[" << i << "] violated by " << (*theFvec)[i] - theUBbound[i] - allow << std::endl );

         if (theUBbound[i] != theLBbound[i])
            shiftUBbound(i, (*theFvec)[i] + random.next(minrandom, maxrandom));
         else
         {
            shiftUBbound(i, (*theFvec)[i]);
            theLBbound[i] = theUBbound[i];
         }
      }
      else if ((*theFvec)[i] < theLBbound[i] - allow)
      {
         MSG_DEBUG( std::cout << "DSHIFT08 theLBbound[" << i << "] violated by " << theLBbound[i] - (*theFvec)[i] - allow << std::endl );

         if (theUBbound[i] != theLBbound[i])
            shiftLBbound(i, (*theFvec)[i] - random.next(minrandom, maxrandom));
         else
         {
            shiftLBbound(i, (*theFvec)[i]);
            theUBbound[i] = theLBbound[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   MSG_DEBUG( std::cout << "DSHIFT01 shiftFvec: OK" << std::endl; )
#endif
}

// -----------------------------------------------------------------

/*
    This methods assumes correctly setup vectors |pVec| and |coPvec| and bound
    vectors for leaving simplex. Then it checks all values of |pVec| and
    |coPvec| to obey these bounds and enlarges them if neccessary.
 */
void SPxSolver::shiftPvec()
{

   /* the allowed tolerance is (rep() == ROW) ? feastol() : opttol() because thePvec is the primal vector in ROW and the
    * dual vector in COLUMN representation; this is equivalent to leavetol()
    */
   Real minrandom = 10.0 * leavetol();
   Real maxrandom = 100.0 * leavetol();
   Real allow = leavetol() - epsilon();
   bool tmp;
   int i;

   assert(type() == LEAVE);
   assert(allow > 0.0);

   for (i = dim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(coId(i));
      if ((*theCoUbound)[i] + allow <= (*theCoPvec)[i] && tmp)
      {
         if ((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftUCbound(i, (*theCoPvec)[i] + random.next(minrandom,maxrandom));
         else
         {
            shiftUCbound(i, (*theCoPvec)[i]);
            (*theCoLbound)[i] = (*theCoUbound)[i];
         }
      }
      else if ((*theCoLbound)[i] - allow >= (*theCoPvec)[i] && tmp)
      {
         if ((*theCoUbound)[i] != (*theCoLbound)[i])
            shiftLCbound(i, (*theCoPvec)[i] - random.next(minrandom,maxrandom));
         else
         {
            shiftLCbound(i, (*theCoPvec)[i]);
            (*theCoUbound)[i] = (*theCoLbound)[i];
         }
      }
   }

   for (i = coDim() - 1; i >= 0; --i)
   {
      tmp = !isBasic(id(i));
      if ((*theUbound)[i] + allow <= (*thePvec)[i] && tmp)
      {
         if ((*theUbound)[i] != (*theLbound)[i])
            shiftUPbound(i, (*thePvec)[i] + random.next(minrandom,maxrandom));
         else
         {
            shiftUPbound(i, (*thePvec)[i]);
            (*theLbound)[i] = (*theUbound)[i];
         }
      }
      else if ((*theLbound)[i] - allow >= (*thePvec)[i] && tmp)
      {
         if ((*theUbound)[i] != (*theLbound)[i])
            shiftLPbound(i, (*thePvec)[i] - random.next(minrandom,maxrandom));
         else
         {
            shiftLPbound(i, (*thePvec)[i]);
            (*theUbound)[i] = (*theLbound)[i];
         }
      }
   }

#ifndef NDEBUG
   testBounds();
   MSG_DEBUG( std::cout << "DSHIFT02 shiftPvec: OK" << std::endl; )
#endif
}
// -----------------------------------------------------------------
void SPxSolver::perturbMin(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real p_delta,
   int start,
   int incr)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   Real minrandom = 10.0 * p_delta;
   Real maxrandom = 100.0 * p_delta;
   Real x, l, u;
   int i;

   if( fullPerturbation )
   {
      eps = p_delta;

      for( i = uvec.dim() - start - 1; i >= 0; i -= incr )
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if( LT(u, infinity) && NE(l, u) && u <= x + eps )
         {
            p_up[i] = x + random.next(minrandom,maxrandom);
            theShift += p_up[i] - u;
         }
         if( GT(l, -infinity) && NE(l, u) && l >= x - eps )
         {
            p_low[i] = x - random.next(minrandom,maxrandom);
            theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const Real* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for( int j = uvec.delta().size() - start - 1; j >= 0; j -= incr )
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];

         // do not permute these bounds! c.f. with computeFrhs2() in spxvecs.cpp
         if( dualStatus(baseId(i)) == SPxBasis::Desc::D_ON_BOTH )
         {
            continue;
         }

         if (x < -eps)
         {
            if( LT(u, infinity) && NE(l, u) && vec[i] >= u - eps )
            {
               p_up[i] = vec[i] + random.next(minrandom,maxrandom);
               theShift += p_up[i] - u;
            }
         }
         else if (x > eps)
         {
            if( GT(l, -infinity) && NE(l, u) && vec[i] <= l + eps )
            {
               p_low[i] = vec[i] - random.next(minrandom,maxrandom);
               theShift -= p_low[i] - l;
            }
         }
      }
   }
}
// -----------------------------------------------------------------
void SPxSolver::perturbMax(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real p_delta,
   int start,
   int incr) 
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   Real minrandom = 10.0 * p_delta;
   Real maxrandom = 100.0 * p_delta;
   Real x, l, u;
   int i;

   if( fullPerturbation )
   {
      eps = p_delta;
      for (i = uvec.dim() - start - 1; i >= 0; i -= incr)
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if( LT(u, infinity) && NE(l, u) && u <= x + eps )
         {
            p_up[i] = x + random.next(minrandom,maxrandom);
            theShift += p_up[i] - u;
         }
         if( GT(l, -infinity) && NE(l, u) && l >= x - eps )
         {
            p_low[i] = x - random.next(minrandom,maxrandom);
            theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const Real* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for( int j = uvec.delta().size() - start - 1; j >= 0; j -= incr )
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];

         // do not permute these bounds! c.f. computeFrhs2() in spxvecs.cpp
         if( dualStatus(baseId(i)) == SPxBasis::Desc::D_ON_BOTH )
         {
            continue;
         }

         if (x > eps)
         {
            if( LT(u, infinity) && NE(l, u) && vec[i] >= u - eps )
            {
               p_up[i] = vec[i] + random.next(minrandom,maxrandom);
               theShift += p_up[i] - u;
            }
         }
         else if (x < -eps)
         {
            if( GT(l, -infinity) && NE(l, u) && vec[i] <= l + eps )
            {
               p_low[i] = vec[i] - random.next(minrandom,maxrandom);
               theShift -= p_low[i] - l;
            }
         }
      }
   }
}

void SPxSolver::perturbMinEnter(void)
{
   MSG_DEBUG( std::cout << "DSHIFT03 iteration= " << iteration() << ": perturbing " << shift(); )
   fVec().delta().setup();
   perturbMin(fVec(), lbBound(), ubBound(), epsilon(), entertol());
   MSG_DEBUG( std::cout << "\t->" << shift() << std::endl; )
}


void SPxSolver::perturbMaxEnter(void)
{
   MSG_DEBUG( std::cout << "DSHIFT04 iteration= " << iteration() << ": perturbing " << shift(); )
   fVec().delta().setup();
   perturbMax(fVec(), lbBound(), ubBound(), epsilon(), entertol());
   MSG_DEBUG( std::cout << "\t->" << shift() << std::endl; )
}


Real SPxSolver::perturbMin(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real p_delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   Real minrandom = 10.0 * p_delta;
   Real maxrandom = 100.0 * p_delta;
   Real x, l, u;
   int i;
   Real l_theShift = 0;

   if( fullPerturbation )
   {
      eps = p_delta;
      for( i = uvec.dim() - start - 1; i >= 0; i -= incr )
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if( LT(u, infinity) && NE(l, u) && u <= x + eps && rep() * stat[i] < 0 )
         {
            p_up[i] = vec[i] + random.next(minrandom,maxrandom);
            l_theShift += p_up[i] - u;
         }
         if( GT(l, -infinity) && NE(l, u) && l >= x - eps && rep() * stat[i] < 0 )
         {
            p_low[i] = vec[i] - random.next(minrandom,maxrandom);
            l_theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const Real* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for( int j = uvec.delta().size() - start - 1; j >= 0; j -= incr )
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];
         if (x < -eps)
         {
            if( LT(u, infinity) && NE(l, u) && vec[i] >= u - eps && rep() * stat[i] < 0 )
            {
               p_up[i] = vec[i] + random.next(minrandom,maxrandom);
               l_theShift += p_up[i] - u;
            }
         }
         else if (x > eps)
         {
            if( GT(l, -infinity) && NE(l, u) && vec[i] <= l + eps && rep() * stat[i] < 0 )
            {
               p_low[i] = vec[i] - random.next(minrandom,maxrandom);
               l_theShift -= p_low[i] - l;
            }
         }
      }
   }
   return l_theShift;
}

Real SPxSolver::perturbMax(
   const UpdateVector& uvec,
   Vector& p_low,
   Vector& p_up,
   Real eps,
   Real p_delta,
   const SPxBasis::Desc::Status* stat,
   int start,
   int incr)
{
   assert(uvec.dim() == p_low.dim());
   assert(uvec.dim() == p_up.dim());

   const Real* vec = uvec.get_const_ptr();
   Real minrandom = 10.0 * p_delta;
   Real maxrandom = 100.0 * p_delta;
   Real x, l, u;
   int i;
   Real l_theShift = 0;

   if( fullPerturbation )
   {
      eps = p_delta;
      for( i = uvec.dim() - start - 1; i >= 0; i -= incr )
      {
         u = p_up[i];
         l = p_low[i];
         x = vec[i];

         if( LT(u, infinity) && NE(l, u) && u <= x + eps && rep() * stat[i] < 0 )
         {
            p_up[i] = vec[i] + random.next(minrandom,maxrandom);
            l_theShift += p_up[i] - u;
         }
         if( GT(l, -infinity) && NE(l, u) && l >= x - eps && rep() * stat[i] < 0 )
         {
            p_low[i] = vec[i] - random.next(minrandom,maxrandom);
            l_theShift -= p_low[i] - l;
         }
      }
   }
   else
   {
      const Real* upd = uvec.delta().values();
      const IdxSet& idx = uvec.delta().indices();

      for( int j = uvec.delta().size() - start - 1; j >= 0; j -= incr )
      {
         i = idx.index(j);
         x = upd[i];
         u = p_up[i];
         l = p_low[i];
         if( x > eps )
         {
            if( LT(u, infinity) && NE(l, u) && vec[i] >= u - eps && rep() * stat[i] < 0 )
            {
               p_up[i] = vec[i] + random.next(minrandom,maxrandom);
               l_theShift += p_up[i] - u;
            }
         }
         else if( x < -eps )
         {
            if( GT(l, -infinity) && NE(l, u) && vec[i] <= l + eps && rep() * stat[i] < 0 )
            {
               p_low[i] = vec[i] - random.next(minrandom,maxrandom);
               l_theShift -= p_low[i] - l;
            }
         }
      }
   }
   return l_theShift;
}


void SPxSolver::perturbMinLeave(void)
{
   MSG_DEBUG( std::cout << "DSHIFT05 iteration= " << iteration() << ": perturbing " << shift(); )
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMin(pVec(), lpBound(), upBound(), epsilon(), leavetol(),
      desc().status(), 0, 1);
   theShift += perturbMin(coPvec(), lcBound(), ucBound(), epsilon(), leavetol(),
      desc().coStatus(), 0, 1);
   MSG_DEBUG( std::cout << "\t->" << shift() << std::endl; )
}


void SPxSolver::perturbMaxLeave(void)
{
   MSG_DEBUG( std::cout << "DSHIFT06 iteration= " << iteration() << ": perturbing " << shift(); )
   pVec().delta().setup();
   coPvec().delta().setup();
   theShift += perturbMax(pVec(), lpBound(), upBound(), epsilon(), leavetol(),
      desc().status(), 0, 1);
   theShift += perturbMax(coPvec(), lcBound(), ucBound(), epsilon(), leavetol(),
      desc().coStatus(), 0, 1);
   MSG_DEBUG( std::cout << "\t->" << shift() << std::endl; )
}


void SPxSolver::unShift(void)
{
   MSG_INFO3( (*spxout), (*spxout) << "DSHIFT07 = " << "unshifting ..." << std::endl; );

   if (isInitialized())
   {
      int i;
      Real t_up, t_low;
      const SPxBasis::Desc& ds = desc();

      theShift = 0;
      if (type() == ENTER)
      {
         Real eps = entertol();

         if (rep() == COLUMN)
         {
            for (i = dim(); i-- > 0;)
            {
               SPxId l_id = baseId(i);
               int l_num = number(l_id);
               if (l_id.type() == SPxId::ROW_ID)
               {
                  t_up = -lhs(l_num);
                  t_low = -rhs(l_num);
               }
               else
               {
                  assert(l_id.type() == SPxId::COL_ID);
                  t_up = upper(l_num);
                  t_low = lower(l_num);
               }
               if (t_up != t_low)
               {
                  if ((*theFvec)[i] < t_up + eps) // check allowed violation
                     theUBbound[i] = t_up; // reset shifted bound to original
                  else if ((*theFvec)[i] > t_up) // shifted bound is required for feasibility
                     theShift += theUBbound[i] - t_up;
                  if ((*theFvec)[i] > t_low - eps) // check allowed violation
                     theLBbound[i] = t_low; // reset shifted bound to original
                  else if ((*theFvec)[i] < t_low) // shifted bound is required for feasibility
                     theShift -= theLBbound[i] - t_low;
               }
               else
               {
                  if (theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else if (theLBbound[i] < t_low)
                     theShift += t_low - theLBbound[i];
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  t_up = -lhs(i);
                  t_low = -rhs(i);
                  if (theURbound[i] > t_up) // what about t_up == t_low ?
                     theShift += theURbound[i] - t_up;
                  if (t_low > theLRbound[i]) // what about t_up == t_low ?
                     theShift += t_low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  t_up = upper(i);
                  t_low = lower(i);
                  if (theUCbound[i] > t_up) // what about t_up == t_low ?
                     theShift += theUCbound[i] - t_up;
                  if (t_low > theLCbound[i]) // what about t_up == t_low ?
                     theShift += t_low - theLCbound[i];
               }
            }
         }
         else
         {
            assert(rep() == ROW);
            for (i = dim(); i-- > 0;)
            {
               SPxId l_id = baseId(i);
               int l_num = number(l_id);
               t_up = t_low = 0;
               if (l_id.type() == SPxId::ROW_ID)
                  clearDualBounds(ds.rowStatus(l_num), t_up, t_low);
               else
                  clearDualBounds(ds.colStatus(l_num), t_up, t_low);
               if (theUBbound[i] != theLBbound[i])
               {
                  if (theUBbound[i] > t_up)
                     theShift += theUBbound[i] - t_up;
                  else
                     theShift -= theUBbound[i] - t_up;
               }
               else
               {
                  /* if the basic (primal or dual) variable is fixed (e.g., basis status P_FREE in row representation)
                   * then shiftFvec() and shiftPvec() do not relax the bounds, but shift both, hence they may be outside
                   * of [t_low,t_up] */
                  assert(theLBbound[i] == theUBbound[i] || theUBbound[i] >= t_up);
                  assert(theLBbound[i] == theUBbound[i] || theLBbound[i] <= t_low);

                  if ((*theFvec)[i] < t_up - eps)
                     theUBbound[i] = t_up;
                  else if ((*theFvec)[i] > t_up)
                     theShift += theUBbound[i] - t_up;

                  if ((*theFvec)[i] > t_low + eps)
                     theLBbound[i] = t_low;
                  else if ((*theFvec)[i] < t_low)
                     theShift -= theLBbound[i] - t_low;
               }
            }
            for (i = nRows(); i-- > 0;)
            {
               if (!isBasic(ds.rowStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.rowStatus(i), t_up, t_low);
                  if (theURbound[i] > t_up) // what about t_up == t_low ?
                     theShift += theURbound[i] - t_up;
                  if (t_low > theLRbound[i]) // what about t_up == t_low ?
                     theShift += t_low - theLRbound[i];
               }
            }
            for (i = nCols(); i-- > 0;)
            {
               if (!isBasic(ds.colStatus(i)))
               {
                  t_up = t_low = 0;
                  clearDualBounds(ds.colStatus(i), t_up, t_low);
                  if (theUCbound[i] > t_up) // what about t_up == t_low ?
                     theShift += theUCbound[i] - t_up;
                  if (t_low > theLCbound[i]) // what about t_up == t_low ?
                     theShift += t_low - theLCbound[i];
               }
            }
         }
      }
      else
      {
         assert(type() == LEAVE);

         Real eps = leavetol();

         if (rep() == COLUMN)
         {
            for (i = nRows(); i-- > 0;)
            {
               t_up = t_low = maxRowObj(i);
               clearDualBounds(ds.rowStatus(i), t_up, t_low);
               if (!isBasic(ds.rowStatus(i)))
               {
                  if ((*theCoPvec)[i] < t_up + eps)
                  {
                     theURbound[i] = t_up; // reset bound to original value
                     if( t_up == t_low )
                        theLRbound[i] = t_low; // for fixed rows we change both bounds
                  }
                  else
                     theShift += theURbound[i] - t_up;
                  if ((*theCoPvec)[i] > t_low - eps)
                  {
                     theLRbound[i] = t_low; // reset bound to original value
                     if( t_up == t_low )
                        theURbound[i] = t_up; // for fixed rows we change both bounds
                  }
                  else
                     theShift += t_low - theLRbound[i];
               }
               else if (theURbound[i] > t_up)
                  theShift += theURbound[i] - t_up;
               else if (theLRbound[i] < t_low)
                  theShift += t_low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               t_up = t_low = -maxObj(i);
               clearDualBounds(ds.colStatus(i), t_low, t_up);
               if (!isBasic(ds.colStatus(i)))
               {
                  if ((*thePvec)[i] < -t_up + eps)
                  {
                     theUCbound[i] = -t_up; // reset bound to original value
                     if( t_up == t_low )
                        theLCbound[i] = -t_low; // for fixed variables we change both bounds
                  }
                  else
                     theShift += theUCbound[i] - (-t_up);
                  if ((*thePvec)[i] > -t_low - eps)
                  {
                     theLCbound[i] = -t_low; // reset bound to original value
                     if( t_up == t_low )
                        theUCbound[i] = -t_up; // for fixed variables we change both bounds
                  }
                  else
                     theShift += (-t_low) - theLCbound[i];
               }
               else if (theUCbound[i] > -t_up)
                  theShift += theUCbound[i] - (-t_up);
               else if (theLCbound[i] < -t_low)
                  theShift += (-t_low) - theLCbound[i];
            }
         }
         else
         {
            assert(rep() == ROW);
            for (i = nRows(); i-- > 0;)
            {
               t_up = rhs(i);
               t_low = lhs(i);
               if (t_up == t_low)
               {
                  if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  else
                     theShift += t_low - theLRbound[i];
               }
               else
                  if (!isBasic(ds.rowStatus(i)))
                  {
                     if ((*thePvec)[i] < t_up + eps)
                        theURbound[i] = t_up; // reset bound to original value
                     else
                        theShift += theURbound[i] - t_up;
                     if ((*thePvec)[i] > t_low - eps)
                        theLRbound[i] = t_low; // reset bound to original value
                     else
                        theShift += t_low - theLRbound[i];
                  }
                  else if (theURbound[i] > t_up)
                     theShift += theURbound[i] - t_up;
                  else if (theLRbound[i] < t_low)
                     theShift += t_low - theLRbound[i];
            }
            for (i = nCols(); i-- > 0;)
            {
               t_up = upper(i);
               t_low = lower(i);
               if (t_up == t_low)
               {
                  if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  else
                     theShift += t_low - theLCbound[i];
               }
               else
                  if (!isBasic(ds.colStatus(i)))
                  {
                     if ((*theCoPvec)[i] < t_up + eps)
                        theUCbound[i] = t_up; // reset bound to original value
                     else
                        theShift += theUCbound[i] - t_up;
                     if ((*theCoPvec)[i] > t_low - eps)
                        theLCbound[i] = t_low; // reset bound to original value
                     else
                        theShift += t_low - theLCbound[i];
                  }
                  else if (theUCbound[i] > t_up)
                     theShift += theUCbound[i] - t_up;
                  else if (theLCbound[i] < t_low)
                     theShift += t_low - theLCbound[i];
            }
         }
      }
   }
}
} // namespace soplex
