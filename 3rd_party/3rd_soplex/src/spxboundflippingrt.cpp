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
#include "spxdefines.h"
#include "spxboundflippingrt.h"
#include "sorter.h"
#include "spxsolver.h"
#include "spxout.h"
#include "spxid.h"

namespace soplex
{

#define LOWSTAB          1e-10
#define MAX_RELAX_COUNT  2
#define LONGSTEP_FREQ    100


/** perform necessary bound flips to restore dual feasibility */
void SPxBoundFlippingRT::flipAndUpdate(
   int&                  nflips              /**< number of bounds that should be flipped */
   )
{
   assert(nflips > 0);

   // number of bound flips that are not performed
   int skipped;

   updPrimRhs.setup();
   updPrimRhs.reDim(thesolver->dim());
   updPrimVec.reDim(thesolver->dim());
   updPrimRhs.clear();
   updPrimVec.clear();

   skipped = 0;
   for( int i = 0; i < nflips; ++i )
   {
      int idx;
      idx = breakpoints[i].idx;
      if( idx < 0 )
      {
         ++skipped;
         continue;
      }
      Real range;
      Real upper;
      Real lower;
      Real objChange = 0.0;
      SPxBasis::Desc::Status stat;
      SPxBasis::Desc& ds = thesolver->basis().desc();

      range = 0;
      if( breakpoints[i].src == PVEC )
      {
         assert(thesolver->rep() == SPxSolver::COLUMN);
         stat = ds.status(idx);
         upper = thesolver->upper(idx);
         lower = thesolver->lower(idx);
         switch( stat )
         {
            case SPxBasis::Desc::P_ON_UPPER :
               ds.status(idx) = SPxBasis::Desc::P_ON_LOWER;
               range = lower - upper;
               assert((*thesolver->theLbound)[idx] == -infinity);
               (*thesolver->theLbound)[idx] = (*thesolver->theUbound)[idx];
               (*thesolver->theUbound)[idx] = infinity;
               objChange = range * (*thesolver->theLbound)[idx];
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               ds.status(idx) = SPxBasis::Desc::P_ON_UPPER;
               range = upper - lower;
               assert((*thesolver->theUbound)[idx] == infinity);
               (*thesolver->theUbound)[idx] = (*thesolver->theLbound)[idx];
               (*thesolver->theLbound)[idx] = -infinity;
               objChange = range * (*thesolver->theUbound)[idx];
               break;
            default :
               ++skipped;
               MSG_WARNING( (*thesolver->spxout), (*thesolver->spxout) << "PVEC unexpected status: " << stat
                                   << " index: " << idx
                                   << " val: " << thesolver->pVec()[idx]
                                   << " upd: " << thesolver->pVec().delta()[idx]
                                   << " lower: " << lower
                                   << " upper: " << upper
                                   << " bp.val: " << breakpoints[i].val
                                   << std::endl; )
         }
         MSG_DEBUG( std::cout << "PVEC flipped from: " << stat
                           << " index: " << idx
                           << " val: " << thesolver->pVec()[idx]
                           << " upd: " << thesolver->pVec().delta()[idx]
                           << " lower: " << lower
                           << " upper: " << upper
                           << " bp.val: " << breakpoints[i].val
                           << " UCbound: " << thesolver->theUCbound[idx]
                           << " LCbound: " << thesolver->theLCbound[idx]
                           << std::endl; )
         assert(spxAbs(range) < 1e20);
         updPrimRhs.multAdd(range, thesolver->vector(idx));
         if( objChange != 0.0 )
            thesolver->updateNonbasicValue(objChange);
      }
      else if( breakpoints[i].src == COPVEC )
      {
         assert(thesolver->rep() == SPxSolver::COLUMN);
         stat = ds.coStatus(idx);
         upper = thesolver->rhs(idx);
         lower = thesolver->lhs(idx);
         switch( stat )
         {
            case SPxBasis::Desc::P_ON_UPPER :
               ds.coStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
               range = lower - upper;
               assert((*thesolver->theCoUbound)[idx] == infinity);
               (*thesolver->theCoUbound)[idx] = -(*thesolver->theCoLbound)[idx];
               (*thesolver->theCoLbound)[idx] = -infinity;
               objChange = range * (*thesolver->theCoUbound)[idx];
               break;
            case SPxBasis::Desc::P_ON_LOWER :
               ds.coStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
               range = upper - lower;
               assert((*thesolver->theCoLbound)[idx] == -infinity);
               (*thesolver->theCoLbound)[idx] = -(*thesolver->theCoUbound)[idx];
               (*thesolver->theCoUbound)[idx] = infinity;
               objChange = range * (*thesolver->theCoLbound)[idx];
               break;
            default :
               ++skipped;
               MSG_WARNING( (*thesolver->spxout), (*thesolver->spxout) << "COPVEC unexpected status: " << stat
                                   << " index: " << idx
                                   << " val: " << thesolver->coPvec()[idx]
                                   << " upd: " << thesolver->coPvec().delta()[idx]
                                   << " lower: " << lower
                                   << " upper: " << upper
                                   << " bp.val: " << breakpoints[i].val
                                   << std::endl; )
         }
         MSG_DEBUG( std::cout << "COPVEC flipped from: " << stat
                           << " index: " << idx
                           << " val: " << thesolver->coPvec()[idx]
                           << " upd: " << thesolver->coPvec().delta()[idx]
                           << " lower: " << lower
                           << " upper: " << upper
                           << " bp.val: " << breakpoints[i].val
                           << " URbound: " << thesolver->theURbound[idx]
                           << " LRbound: " << thesolver->theLRbound[idx]
                           << std::endl; )
         assert(spxAbs(range) < 1e20);
         updPrimRhs.setValue(idx, updPrimRhs[idx] - range);
         if( objChange != 0.0 )
            thesolver->updateNonbasicValue(objChange);
      }
      else if( breakpoints[i].src == FVEC )
      {
         assert(thesolver->rep() == SPxSolver::ROW);
         SPxId baseId = thesolver->basis().baseId(idx);
         int IdNumber;

         if( baseId.isSPxRowId() )
         {
            IdNumber = thesolver->number(SPxRowId(baseId));
            stat = ds.rowStatus(IdNumber);
            upper = thesolver->rhs(IdNumber);
            lower = thesolver->lhs(IdNumber);
            switch( stat )
            {
               case SPxBasis::Desc::P_ON_UPPER :
                  ds.rowStatus(IdNumber) = SPxBasis::Desc::P_ON_LOWER;
                  range = upper - lower;
                  assert(thesolver->theUBbound[idx] == infinity);
                  thesolver->theUBbound[idx] = -thesolver->theLBbound[idx];
                  thesolver->theLBbound[idx] = -infinity;
                  break;
               case SPxBasis::Desc::P_ON_LOWER :
                  ds.rowStatus(IdNumber) = SPxBasis::Desc::P_ON_UPPER;
                  range = lower - upper;
                  assert(thesolver->theLBbound[idx] == -infinity);
                  thesolver->theLBbound[idx] = -thesolver->theUBbound[idx];
                  thesolver->theUBbound[idx] = infinity;
                  break;
               default :
                  ++skipped;
                  MSG_WARNING( (*thesolver->spxout), (*thesolver->spxout) << "unexpected basis status: " << stat
                                    << " index: " << idx
                                    << " val: " << thesolver->fVec()[idx]
                                    << " upd: " << thesolver->fVec().delta()[idx]
                                    << " lower: " << lower
                                    << " upper: " << upper
                                    << " bp.val: " << breakpoints[i].val
                                    << std::endl; )
            }
         }
         else
         {
            assert(baseId.isSPxColId());
            IdNumber = thesolver->number(SPxColId(baseId));
            stat = ds.colStatus(IdNumber);
            upper = thesolver->upper(IdNumber);
            lower = thesolver->lower(IdNumber);

            switch( stat )
            {
               case SPxBasis::Desc::P_ON_UPPER :
                  ds.colStatus(IdNumber) = SPxBasis::Desc::P_ON_LOWER;
                  range = upper - lower;
                  assert(thesolver->theUBbound[idx] == infinity);
                  thesolver->theUBbound[idx] = -thesolver->theLBbound[idx];
                  thesolver->theLBbound[idx] = -infinity;
                  break;
               case SPxBasis::Desc::P_ON_LOWER :
                  ds.colStatus(IdNumber) = SPxBasis::Desc::P_ON_UPPER;
                  range = lower - upper;
                  assert(thesolver->theLBbound[idx] == -infinity);
                  thesolver->theLBbound[idx] = -thesolver->theUBbound[idx];
                  thesolver->theUBbound[idx] = infinity;
                  break;
               default :
                  ++skipped;
                  MSG_WARNING( (*thesolver->spxout), (*thesolver->spxout) << "FVEC unexpected status: " << stat
                                    << " index: " << idx
                                    << " val: " << thesolver->fVec()[idx]
                                    << " upd: " << thesolver->fVec().delta()[idx]
                                    << " lower: " << lower
                                    << " upper: " << upper
                                    << " bp.val: " << breakpoints[i].val
                                    << std::endl; )
            }
         }
         MSG_DEBUG( std::cout << "basic row/col flipped from: " << stat
                           << " index: " << idx
                           << " val: " << thesolver->fVec()[idx]
                           << " upd: " << thesolver->fVec().delta()[idx]
                           << " lower: " << lower
                           << " upper: " << upper
                           << " bp.val: " << breakpoints[i].val
                           << std::endl; )
         assert(spxAbs(range) < 1e20);
         assert(updPrimRhs[idx] == 0);
         updPrimRhs.add(idx, range);
      }
   }
   nflips -= skipped;
   if( nflips > 0 )
   {
      if(thesolver->rep() == SPxSolver::ROW)
      {
         assert(m_type == SPxSolver::ENTER);
         (*thesolver->theCoPrhs) -= updPrimRhs;
         thesolver->setup4coSolve2(&updPrimVec, &updPrimRhs);
      }
      else
      {
         assert(thesolver->rep() == SPxSolver::COLUMN);
         assert(m_type == SPxSolver::LEAVE);
         (*thesolver->theFrhs) -= updPrimRhs;
         thesolver->setup4solve2(&updPrimVec, &updPrimRhs);
      }
   }

   return;
}

/** store all available pivots/breakpoints in an array (positive pivot search direction) */
void SPxBoundFlippingRT::collectBreakpointsMax(
   int&                  nBp,                /**< number of found breakpoints so far */
   int&                  minIdx,             /**< index to current minimal breakpoint */
   const int*            idx,                /**< pointer to indices of current vector */
   int                   nnz,                /**< number of nonzeros in current vector */
   const Real*           upd,                /**< pointer to update values of current vector */
   const Real*           vec,                /**< pointer to values of current vector */
   const Real*           upp,                /**< pointer to upper bound/rhs of current vector */
   const Real*           low,                /**< pointer to lower bound/lhs of current vector */
   BreakpointSource      src                 /**< type of vector (pVec, coPvec or fVec)*/
   )
{
   Real minVal;
   Real curVal;
   const int* last;

   minVal = ( nBp == 0 ) ? infinity : breakpoints[minIdx].val;

   last = idx + nnz;
   for( ; idx < last; ++idx )
   {
      int i = *idx;
      Real x = upd[i];
      if( x > epsilon )
      {
         if( upp[i] < infinity )
         {
            Real y = upp[i] - vec[i];
            curVal = (y <= 0) ? fastDelta / x : (y + fastDelta) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      else if( x < -epsilon )
      {
         if (low[i] > -infinity)
         {
            Real y = low[i] - vec[i];
            curVal = (y >= 0) ? -fastDelta / x : (y - fastDelta) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      if( nBp >= breakpoints.size() )
         breakpoints.reSize(nBp * 2);
   }

   return;
}

/** store all available pivots/breakpoints in an array (negative pivot search direction) */
void SPxBoundFlippingRT::collectBreakpointsMin(
   int&                  nBp,                /**< number of found breakpoints so far */
   int&                  minIdx,             /**< index to current minimal breakpoint */
   const int*            idx,                /**< pointer to indices of current vector */
   int                   nnz,                /**< number of nonzeros in current vector */
   const Real*           upd,                /**< pointer to update values of current vector */
   const Real*           vec,                /**< pointer to values of current vector */
   const Real*           upp,                /**< pointer to upper bound/rhs of current vector */
   const Real*           low,                /**< pointer to lower bound/lhs of current vector */
   BreakpointSource      src                 /**< type of vector (pVec, coPvec or fVec)*/
   )
{
   Real minVal;
   Real curVal;
   const int* last;

   minVal = ( nBp == 0 ) ? infinity : breakpoints[minIdx].val;

   last = idx + nnz;

   for( ; idx < last; ++idx )
   {
      int i = *idx;
      Real x = upd[i];
      if( x > epsilon )
      {
         if( low[i] > -infinity )
         {
            Real y = low[i] - vec[i];

            curVal = (y >= 0) ? fastDelta / x : (fastDelta - y) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      else if( x < -epsilon )
      {
         if (upp[i] < infinity)
         {
            Real y = upp[i] - vec[i];
            curVal = (y <= 0) ? -fastDelta / x : -(y + fastDelta) / x;
            assert(curVal > 0);

            breakpoints[nBp].idx = i;
            breakpoints[nBp].src = src;
            breakpoints[nBp].val = curVal;

            if( curVal < minVal )
            {
               minVal = curVal;
               minIdx = nBp;
            }

            nBp++;
         }
      }
      if( nBp >= breakpoints.size() )
         breakpoints.reSize(nBp * 2);
   }
   return;
}

/** get values for entering index and perform shifts if necessary */
bool SPxBoundFlippingRT::getData(
   Real&                 val,
   SPxId&                enterId,
   int                   idx,
   Real                  stab,
   Real                  degeneps,
   const Real*           upd,
   const Real*           vec,
   const Real*           low,
   const Real*           upp,
   BreakpointSource      src,
   Real                  max
   )
{
   if( src == PVEC )
   {
      thesolver->pVec()[idx] = thesolver->vector(idx) * thesolver->coPvec();
      Real x = upd[idx];
      // skip breakpoint if it is too small
      if( spxAbs(x) < stab )
      {
         return false;
      }
      enterId = thesolver->id(idx);
      val = (max * x > 0) ? upp[idx] : low[idx];
      val = (val - vec[idx]) / x;
      if( upp[idx] == low[idx] )
      {
         val = 0.0;
         if( vec[idx] > upp[idx] )
            thesolver->theShift += vec[idx] - upp[idx];
         else
            thesolver->theShift += low[idx] - vec[idx];
         thesolver->upBound()[idx] = thesolver->lpBound()[idx] = vec[idx];
      }
      else if( (max > 0 && val < -degeneps) || (max < 0 && val > degeneps) )
      {
         val = 0.0;
         if( max * x > 0 )
            thesolver->shiftUPbound(idx, vec[idx]);
         else
            thesolver->shiftLPbound(idx, vec[idx]);
      }
   }
   else // src == COPVEC
   {
      Real x = upd[idx];
      if( spxAbs(x) < stab )
      {
         return false;
      }
      enterId = thesolver->coId(idx);
      val = (max * x > 0.0) ? upp[idx] : low[idx];
      val = (val - vec[idx]) / x;
      if( upp[idx] == low[idx] )
      {
         val = 0.0;
         if( vec[idx] > upp[idx] )
            thesolver->theShift += vec[idx] - upp[idx];
         else
            thesolver->theShift += low[idx] - vec[idx];
         thesolver->ucBound()[idx] = thesolver->lcBound()[idx] = vec[idx];
      }
      else if( (max > 0 && val < -degeneps) || (max < 0 && val > degeneps) )
      {
         val = 0.0;
         if( max * x > 0 )
            thesolver->shiftUCbound(idx, vec[idx]);
         else
            thesolver->shiftLCbound(idx, vec[idx]);
      }
   }
   return true;
}

/** get values for leaving index and perform shifts if necessary */
bool SPxBoundFlippingRT::getData(
   Real&                 val,
   int&                  leaveIdx,
   int                   idx,
   Real                  stab,
   Real                  degeneps,
   const Real*           upd,
   const Real*           vec,
   const Real*           low,
   const Real*           upp,
   BreakpointSource      src,
   Real                  max
   )
{
   assert( src == FVEC );

   Real x = upd[idx];
   // skip breakpoint if it is too small
   if( spxAbs(x) < stab )
   {
      return false;
   }
   leaveIdx = idx;
   val = (max * x > 0) ? upp[idx] : low[idx];
   val = (val - vec[idx]) / x;
   if( upp[idx] == low[idx] )
   {
      val = 0.0;
      thesolver->shiftLBbound(idx, vec[idx]);
      thesolver->shiftUBbound(idx, vec[idx]);
   }
   else if( (max > 0 && val < -degeneps) || (max < 0 && val > degeneps) )
   {
      val = 0.0;
      if( thesolver->dualStatus(thesolver->baseId(idx)) != SPxBasis::Desc::D_ON_BOTH )
      {
         if( max * x > 0 )
            thesolver->shiftUBbound(idx, vec[idx]);
         else
            thesolver->shiftLBbound(idx, vec[idx]);
      }
   }
   return true;
}

/** determine entering row/column */
SPxId SPxBoundFlippingRT::selectEnter(
   Real&                 val,
   int                   leaveIdx,
   bool                  polish
   )
{
   assert( m_type == SPxSolver::LEAVE );
   assert(thesolver->boundflips == 0);

   // reset the history and try again to do some long steps
   if( thesolver->leaveCount % LONGSTEP_FREQ == 0 )
   {
      MSG_DEBUG( std::cout << "DLBFRT06 resetting long step history" << std::endl; )
      flipPotential = 1;
   }
   if( !enableBoundFlips || polish || thesolver->rep() == SPxSolver::ROW || flipPotential <= 0 )
   {
      MSG_DEBUG( std::cout << "DLBFRT07 switching to fast ratio test" << std::endl; )
      return SPxFastRT::selectEnter(val, leaveIdx, polish);
   }
   const Real*  pvec = thesolver->pVec().get_const_ptr();
   const Real*  pupd = thesolver->pVec().delta().values();
   const int*   pidx = thesolver->pVec().delta().indexMem();
   int          pupdnnz = thesolver->pVec().delta().size();
   const Real*  lpb  = thesolver->lpBound().get_const_ptr();
   const Real*  upb  = thesolver->upBound().get_const_ptr();

   const Real*  cvec = thesolver->coPvec().get_const_ptr();
   const Real*  cupd = thesolver->coPvec().delta().values();
   const int*   cidx = thesolver->coPvec().delta().indexMem();
   int          cupdnnz = thesolver->coPvec().delta().size();
   const Real*  lcb  = thesolver->lcBound().get_const_ptr();
   const Real*  ucb  = thesolver->ucBound().get_const_ptr();

   resetTols();

   Real max;

   // index in breakpoint array of minimal value (i.e. choice of normal RT)
   int minIdx;

   // temporary breakpoint data structure to make swaps possible
   Breakpoint tmp;

   // most stable pivot value in candidate set
   Real moststable;

   // initialize invalid enterId
   SPxId enterId;

   // slope of objective function improvement
   Real slope;

   // number of found breakpoints
   int nBp;

   // number of passed breakpoints
   int npassedBp;

   Real degeneps;
   Real stab;
   bool instable;

   max = val;
   val = 0.0;
   moststable = 0.0;
   nBp = 0;
   minIdx = -1;

   // get breakpoints and and determine the index of the minimal value
   if( max > 0 )
   {
      collectBreakpointsMax(nBp, minIdx, pidx, pupdnnz, pupd, pvec, upb, lpb, PVEC);
      collectBreakpointsMax(nBp, minIdx, cidx, cupdnnz, cupd, cvec, ucb, lcb, COPVEC);
   }
   else
   {
      collectBreakpointsMin(nBp, minIdx, pidx, pupdnnz, pupd, pvec, upb, lpb, PVEC);
      collectBreakpointsMin(nBp, minIdx, cidx, cupdnnz, cupd, cvec, ucb, lcb, COPVEC);
   }

   if( nBp == 0 )
   {
      val = max;
      return enterId;
   }

   assert(minIdx >= 0);

   // swap smallest breakpoint to the front to skip the sorting phase if no bound flip is possible
   tmp = breakpoints[minIdx];
   breakpoints[minIdx] = breakpoints[0];
   breakpoints[0] = tmp;

   // get initial slope
   slope = spxAbs(thesolver->fTest()[leaveIdx]);
   if( slope == 0 )
   {
      // this may only happen if SoPlex decides to make an instable pivot
      assert(thesolver->instableLeaveNum >= 0);
      // restore original slope
      slope = spxAbs(thesolver->instableLeaveVal);
   }

   // set up structures for the quicksort implementation
   BreakpointCompare compare;
   compare.entry = breakpoints.get_const_ptr();

   // pointer to end of sorted part of breakpoints
   int sorted = 0;
   // minimum number of entries that are supposed to be sorted by partial sort
   int sortsize = 4;

   // get all skipable breakpoints
   for( npassedBp = 0; npassedBp < nBp && slope > 0; ++npassedBp)
   {
      // sort breakpoints only partially to save time
      if( npassedBp > sorted )
      {
         sorted = SPxQuicksortPart(breakpoints.get_ptr(), compare, sorted + 1, nBp, sortsize);
      }
      int i = breakpoints[npassedBp].idx;
      // compute new slope
      if( breakpoints[npassedBp].src == PVEC )
      {
         if( thesolver->isBasic(i) )
         {
            // mark basic indices
            breakpoints[npassedBp].idx = -1;
            thesolver->pVec().delta().clearIdx(i);
         }
         else
         {
            Real absupd = spxAbs(pupd[i]);
            slope -= (thesolver->upper(i) * absupd) - (thesolver->lower(i) * absupd);
            // get most stable pivot
            if( absupd > moststable )
               moststable = absupd;
         }
      }
      else
      {
         assert(breakpoints[npassedBp].src == COPVEC);
         if( thesolver->isCoBasic(i) )
         {
            // mark basic indices
            breakpoints[npassedBp].idx = -1;
            thesolver->coPvec().delta().clearIdx(i);
         }
         else
         {
            Real absupd = spxAbs(cupd[i]);
            slope -= (thesolver->rhs(i) * absupd) - (thesolver->lhs(i) * absupd);
            if( absupd > moststable )
               moststable = absupd;
         }
      }
   }
   --npassedBp;
   assert(npassedBp >= 0);

   // check for unboundedness/infeasibility
   if( slope > delta && npassedBp >= nBp - 1 )
   {
      MSG_DEBUG( std::cout << "DLBFRT02 " << thesolver->basis().iteration()
                        << ": unboundedness in ratio test" << std::endl; )
      flipPotential -= 0.5;
      val = max;
      return SPxFastRT::selectEnter(val, leaveIdx);
   }

   MSG_DEBUG( std::cout << "DLBFRT01 "
                     << thesolver->basis().iteration()
                     << ": number of flip candidates: "
                     << npassedBp
                     << std::endl; )

   // try to get a more stable pivot by looking at those with similar step length
   int stableBp;              // index to walk over additional breakpoints (after slope change)
   int bestBp = -1;           // breakpoints index with best possible stability
   Real bestDelta = breakpoints[npassedBp].val;  // best step length (after bound flips)

   for( stableBp = npassedBp + 1; stableBp < nBp; ++stableBp )
   {
      Real stableDelta = 0;
      // get next breakpoints in increasing order
      if( stableBp > sorted )
      {
         sorted = SPxQuicksortPart(breakpoints.get_ptr(), compare, sorted + 1, nBp, sortsize);
      }
      int idx = breakpoints[stableBp].idx;
      if( breakpoints[stableBp].src == PVEC )
      {
         if( thesolver->isBasic(idx) )
         {
            // mark basic indices
            breakpoints[stableBp].idx = -1;
            thesolver->pVec().delta().clearIdx(idx);
            continue;
         }
         Real x = pupd[idx];
         if( spxAbs(x) > moststable )
         {
            thesolver->pVec()[idx] = thesolver->vector(idx) * thesolver->coPvec();
            stableDelta = (x > 0.0) ? upb[idx] : lpb[idx];
            stableDelta = (stableDelta - pvec[idx]) / x;

            if( stableDelta <= bestDelta)
            {
               moststable = spxAbs(x);
               bestBp = stableBp;
            }
         }
      }
      else
      {
         if( thesolver->isCoBasic(idx) )
         {
            // mark basic indices
            breakpoints[stableBp].idx = -1;
            thesolver->coPvec().delta().clearIdx(idx);
            continue;
         }
         Real x = cupd[idx];
         if( spxAbs(x) > moststable )
         {
            stableDelta = (x > 0.0) ? ucb[idx] : lcb[idx];
            stableDelta = (stableDelta - cvec[idx]) / x;

            if( stableDelta <= bestDelta )
            {
               moststable = spxAbs(x);
               bestBp = stableBp;
            }
         }
      }

      // stop searching if the step length is too big
      if( stableDelta > delta + bestDelta )
         break;
   }

   degeneps = fastDelta / moststable;  /* as in SPxFastRT */
   // get stability requirements
   instable = thesolver->instableLeave;
   assert(!instable || thesolver->instableLeaveNum >= 0);
   stab = instable ? LOWSTAB : SPxFastRT::minStability(moststable);

   bool foundStable = false;

   if( bestBp >= 0 )
   {
      // found a more stable pivot
      if( moststable > stab )
      {
         // stability requirements are satisfied
         int idx = breakpoints[bestBp].idx;
         assert(idx >= 0);
         if( breakpoints[bestBp].src == PVEC )
            foundStable = getData(val, enterId, idx, stab, degeneps, pupd, pvec, lpb, upb, PVEC, max);
         else
            foundStable = getData(val, enterId, idx, stab, degeneps, cupd, cvec, lcb, ucb, COPVEC, max);
      }
   }

   else
   {
      // scan passed breakpoints from back to front and stop as soon as a good one is found
      while( !foundStable && npassedBp >= 0 )
      {
         int idx = breakpoints[npassedBp].idx;

         // only look for non-basic variables
         if( idx >= 0 )
         {
            if( breakpoints[npassedBp].src == PVEC )
               foundStable = getData(val, enterId, idx, stab, degeneps, pupd, pvec, lpb, upb, PVEC, max);
            else
               foundStable = getData(val, enterId, idx, stab, degeneps, cupd, cvec, lcb, ucb, COPVEC, max);
         }
         --npassedBp;
      }
      ++npassedBp;
   }

   if( !foundStable )
   {
      assert(!enterId.isValid());
      if( relax_count < MAX_RELAX_COUNT )
      {
         MSG_DEBUG( std::cout << "DLBFRT04 "
                           << thesolver->basis().iteration()
                           << ": no valid enterId found - relaxing..."
                           << std::endl; )
         relax();
         ++relax_count;
         // restore original value
         val = max;
         // try again with relaxed delta
         return SPxBoundFlippingRT::selectEnter(val, leaveIdx);
      }
      else
      {
         MSG_DEBUG( std::cout << "DLBFRT05 "
                           << thesolver->basis().iteration()
                           << " no valid enterId found - breaking..."
                           << std::endl; )
         return enterId;
      }
   }
   else
   {
      relax_count = 0;
      tighten();
   }

   // flip bounds of skipped breakpoints only if a nondegenerate step is to be performed
   if( npassedBp > 0 && spxAbs(breakpoints[npassedBp].val) > fastDelta )
   {
      flipAndUpdate(npassedBp);
      thesolver->boundflips = npassedBp;
      if( npassedBp >= 10 )
         flipPotential = 1;
      else
         flipPotential -= 0.05;
   }
   else
   {
      thesolver->boundflips = 0;
      flipPotential -= 0.1;
   }

   MSG_DEBUG( std::cout << "DLBFRT06 "
                     << thesolver->basis().iteration()
                     << ": selected Id: "
                     << enterId
                     << " number of candidates: "
                     << nBp
                     << std::endl; )
   return enterId;
}

/** determine leaving row/column */
int SPxBoundFlippingRT::selectLeave(
   Real&                 val,
   Real                  enterTest,
   bool                  polish
   )
{
   assert( m_type == SPxSolver::ENTER );
   assert(thesolver->boundflips == 0);

   // reset the history and try again to do some long steps
   if( thesolver->enterCount % LONGSTEP_FREQ == 0 )
   {
      MSG_DEBUG( std::cout << "DEBFRT06 resetting long step history" << std::endl; )
      flipPotential = 1;
   }

   if( polish || !enableBoundFlips || !enableRowBoundFlips || thesolver->rep() == SPxSolver::COLUMN || flipPotential <= 0 )
   {
      MSG_DEBUG( std::cout << "DEBFRT07 switching to fast ratio test" << std::endl; )
      return SPxFastRT::selectLeave(val, enterTest, polish);
   }

   const Real*  vec = thesolver->fVec().get_const_ptr();         /**< pointer to values of current vector */
   const Real*  upd = thesolver->fVec().delta().values();        /**< pointer to update values of current vector */
   const int*   idx = thesolver->fVec().delta().indexMem();      /**< pointer to indices of current vector */
   int          updnnz = thesolver->fVec().delta().size();       /**< number of nonzeros in update vector */
   const Real*  lb  = thesolver->lbBound().get_const_ptr();      /**< pointer to lower bound/lhs of current vector */
   const Real*  ub  = thesolver->ubBound().get_const_ptr();      /**< pointer to upper bound/rhs of current vector */

   resetTols();

   Real max;

   // index in breakpoint array of minimal value (i.e. choice of normal RT)
   int minIdx;

   // temporary breakpoint data structure to make swaps possible
   Breakpoint tmp;

   // most stable pivot value in candidate set
   Real moststable;

   // initialize invalid leaving index
   int leaveIdx = -1;

   // slope of objective function improvement
   Real slope;

   // number of found breakpoints
   int nBp;

   // number of passed breakpoints
   int npassedBp;

   Real degeneps;
   Real stab;
   bool instable;

   max = val;
   val = 0.0;
   moststable = 0.0;
   nBp = 0;
   minIdx = -1;

   assert(thesolver->fVec().delta().isSetup());

   // get breakpoints and and determine the index of the minimal value
   if( max > 0 )
   {
      collectBreakpointsMax(nBp, minIdx, idx, updnnz, upd, vec, ub, lb, FVEC);
   }
   else
   {
      collectBreakpointsMin(nBp, minIdx, idx, updnnz, upd, vec, ub, lb, FVEC);
   }

   // return -1 if no BP was found
   if( nBp == 0 )
   {
      val = max;
      return leaveIdx;
   }

   assert(minIdx >= 0);

   // swap smallest breakpoint to the front to skip the sorting phase if no bound flip is possible
   tmp = breakpoints[minIdx];
   breakpoints[minIdx] = breakpoints[0];
   breakpoints[0] = tmp;

   // get initial slope
   slope = spxAbs(enterTest);
   if( slope == 0 )
   {
      // this may only happen if SoPlex decides to make an instable pivot
      assert(thesolver->instableEnterId.isValid());
      // restore original slope
      slope = thesolver->instableEnterVal;
   }

   // set up structures for the quicksort implementation
   BreakpointCompare compare;
   compare.entry = breakpoints.get_const_ptr();

   // pointer to end of sorted part of breakpoints
   int sorted = 0;
   // minimum number of entries that are supposed to be sorted by partial sort
   int sortsize = 4;

   // get all skipable breakpoints
   for( npassedBp = 0; npassedBp < nBp && slope > 0; ++npassedBp)
   {
      // sort breakpoints only partially to save time
      if( npassedBp > sorted )
      {
         sorted = SPxQuicksortPart(breakpoints.get_ptr(), compare, sorted + 1, nBp, sortsize);
      }
      assert( breakpoints[npassedBp].src == FVEC );
      int breakpointidx = breakpoints[npassedBp].idx;
      // compute new slope
      Real upper;
      Real lower;
      Real absupd = spxAbs(upd[breakpointidx]);
      SPxId baseId = thesolver->baseId(breakpointidx);
      int i = thesolver->number(baseId);
      if( baseId.isSPxColId() )
      {
         upper = thesolver->upper(i);
         lower = thesolver->lower(i);
      }
      else
      {
         assert(baseId.isSPxRowId());
         upper = thesolver->rhs(i);
         lower = thesolver->lhs(i);
      }

      slope -= (upper * absupd) - (lower * absupd);
      // get most stable pivot
      if( absupd > moststable )
         moststable = absupd;
   }
   --npassedBp;
   assert(npassedBp >= 0);

   // check for unboundedness/infeasibility
   if( slope > delta && npassedBp >= nBp - 1 )
   {
      MSG_DEBUG( std::cout << "DEBFRT02 " << thesolver->basis().iteration()
                        << ": unboundedness in ratio test" << std::endl; )
      flipPotential -= 0.5;
      val = max;
      return SPxFastRT::selectLeave(val, enterTest);
   }

   MSG_DEBUG( std::cout << "DEBFRT01 "
                     << thesolver->basis().iteration()
                     << ": number of flip candidates: "
                     << npassedBp
                     << std::endl; )

   // try to get a more stable pivot by looking at those with similar step length
   int stableBp;              // index to walk over additional breakpoints (after slope change)
   int bestBp = -1;           // breakpoints index with best possible stability
   Real bestDelta = breakpoints[npassedBp].val;  // best step length (after bound flips)

   for( stableBp = npassedBp + 1; stableBp < nBp; ++stableBp )
   {
      Real stableDelta = 0;
      // get next breakpoints in increasing order
      if( stableBp > sorted )
      {
         sorted = SPxQuicksortPart(breakpoints.get_ptr(), compare, sorted + 1, nBp, sortsize);
      }
      int breakpointidx = breakpoints[stableBp].idx;
      assert( breakpoints[stableBp].src == FVEC );
      Real x = upd[breakpointidx];
      if( spxAbs(x) > moststable )
      {
         stableDelta = (x > 0.0) ? ub[breakpointidx] : lb[breakpointidx];
         stableDelta = (stableDelta - vec[breakpointidx]) / x;

         if( stableDelta <= bestDelta)
         {
            moststable = spxAbs(x);
            bestBp = stableBp;
         }
      }
      // stop searching if the step length is too big
      else if( stableDelta > delta + bestDelta )
         break;
   }

   degeneps = fastDelta / moststable;  /* as in SPxFastRT */
   // get stability requirements
   instable = thesolver->instableEnter;
   assert(!instable || thesolver->instableEnterId.isValid());
   stab = instable ? LOWSTAB : SPxFastRT::minStability(moststable);

   bool foundStable = false;

   if( bestBp >= 0 )
   {
      // found a more stable pivot
      if( moststable > stab )
      {
         // stability requirements are satisfied
         int breakpointidx = breakpoints[bestBp].idx;
         assert(breakpointidx >= 0);
         foundStable = getData(val, leaveIdx, breakpointidx, moststable, degeneps, upd, vec, lb, ub, FVEC, max);
      }
   }

   else
   {
      // scan passed breakpoints from back to front and stop as soon as a good one is found
      while( !foundStable && npassedBp >= 0 )
      {
         int breakpointidx = breakpoints[npassedBp].idx;

         // only look for non-basic variables
         if( breakpointidx >= 0 )
         {
            foundStable = getData(val, leaveIdx, breakpointidx, moststable, degeneps, upd, vec, lb, ub, FVEC, max);
         }
         --npassedBp;
      }
      ++npassedBp;
   }

   if( !foundStable )
   {
      assert(leaveIdx < 0);
      if( relax_count < MAX_RELAX_COUNT )
      {
         MSG_DEBUG( std::cout << "DEBFRT04 "
                           << thesolver->basis().iteration()
                           << ": no valid leaveIdx found - relaxing..."
                           << std::endl; )
         relax();
         ++relax_count;
         // restore original value
         val = max;
         // try again with relaxed delta
         return SPxBoundFlippingRT::selectLeave(val, enterTest);
      }
      else
      {
         MSG_DEBUG( std::cout << "DEBFRT05 "
                           << thesolver->basis().iteration()
                           << " no valid leaveIdx found - breaking..."
                           << std::endl; )
         return leaveIdx;
      }
   }
   else
   {
      relax_count = 0;
      tighten();
   }

   // flip bounds of skipped breakpoints only if a nondegenerate step is to be performed
   if( npassedBp > 0 && spxAbs(breakpoints[npassedBp].val) > fastDelta )
   {
      flipAndUpdate(npassedBp);
      thesolver->boundflips = npassedBp;
      if( npassedBp >= 10 )
         flipPotential = 1;
      else
         flipPotential -= 0.05;
   }
   else
   {
      thesolver->boundflips = 0;
      flipPotential -= 0.1;
   }

   MSG_DEBUG( std::cout << "DEBFRT06 "
                     << thesolver->basis().iteration()
                     << ": selected Index: "
                     << leaveIdx
                     << " number of candidates: "
                     << nBp
                     << std::endl; )

   return leaveIdx;
}


} // namespace soplex
