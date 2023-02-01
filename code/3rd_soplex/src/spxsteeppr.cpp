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

//TODO may be faster to have a greater zero tolerance for sparse pricing vectors
//     to reduce the number of nonzero entries, e.g. for workVec

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "spxsteeppr.h"
#include "random.h"

#define STEEP_REFINETOL 2.0

namespace soplex
{

// #define EQ_PREF 1000

void SPxSteepPR::clear()
{
   thesolver = 0;
}

void SPxSteepPR::load(SPxSolver* base)
{
   thesolver = base;

   if (base)
   {
      workVec.clear();
      workVec.reDim(base->dim());
      workRhs.clear();
      workRhs.reDim(base->dim());
   }
}

void SPxSteepPR::setType(SPxSolver::Type type)
{
   workRhs.setEpsilon(thesolver->epsilon());

   setupWeights(type);
   workVec.clear();
   workRhs.clear();
   refined = false;

   bestPrices.clear();
   bestPrices.setMax(thesolver->dim());
   prices.reMax(thesolver->dim());

   if( type == SPxSolver::ENTER )
   {
      bestPricesCo.clear();
      bestPricesCo.setMax(thesolver->coDim());
      pricesCo.reMax(thesolver->coDim());
   }
}

void SPxSteepPR::setupWeights(SPxSolver::Type type)
{
   int i;
   int endDim = 0;
   int endCoDim = 0;
   DVector& weights = thesolver->weights;
   DVector& coWeights = thesolver->coWeights;

   if( setup == DEFAULT )
   {
      if( type == SPxSolver::ENTER )
      {
         if( thesolver->weightsAreSetup )
         {
            // check for added/removed rows and adapt norms accordingly
            if (coWeights.dim() < thesolver->dim())
               endDim = coWeights.dim();
            else
               endDim = thesolver->dim();
            if (weights.dim() < thesolver->coDim())
               endCoDim = weights.dim();
            else
               endCoDim = thesolver->coDim();
         }

         coWeights.reDim(thesolver->dim(), false);
         for (i = thesolver->dim() - 1; i >= endDim; --i)
            coWeights[i] = 2.0;
         weights.reDim(thesolver->coDim(), false);
         for (i = thesolver->coDim() - 1; i >= endCoDim; --i)
            weights[i] = 1.0;
      }
      else
      {
         assert(type == SPxSolver::LEAVE);

         if( thesolver->weightsAreSetup )
         {
            // check for added/removed rows and adapt norms accordingly
            if (coWeights.dim() < thesolver->dim())
               endDim = coWeights.dim();
            else
               endDim = thesolver->dim();
         }

         coWeights.reDim(thesolver->dim(), false);
         for (i = thesolver->dim() - 1; i >= endDim; --i)
            coWeights[i]   = 1.0;
      }
   }
   else
   {
      MSG_INFO1( (*thesolver->spxout), (*thesolver->spxout) << " --- initializing steepest edge multipliers" << std::endl; )

      if (type == SPxSolver::ENTER)
      {
         coWeights.reDim(thesolver->dim(), false);
         for (i = thesolver->dim() - 1; i >= endDim; --i)
            coWeights[i] = 1.0;
         weights.reDim(thesolver->coDim(), false);
         for (i = thesolver->coDim() - 1; i >= endCoDim; --i)
            weights[i] = 1.0 + thesolver->vector(i).length2();
      }
      else
      {
         assert(type == SPxSolver::LEAVE);
         coWeights.reDim(thesolver->dim(), false);
         SSVector tmp(thesolver->dim(), thesolver->epsilon());
         for( i = thesolver->dim() - 1; i >= endDim && !thesolver->isTimeLimitReached(); --i )
         {
            thesolver->basis().coSolve(tmp, thesolver->unitVector(i));
            coWeights[i] = tmp.length2();
         }
      }
   }
   thesolver->weightsAreSetup = true;
}

void SPxSteepPR::setRep(SPxSolver::Representation)
{
   if (workVec.dim() != thesolver->dim())
   {
      DVector tmp = thesolver->weights;
      thesolver->weights = thesolver->coWeights;
      thesolver->coWeights = tmp;

      workVec.clear();
      workVec.reDim(thesolver->dim());
   }
}

void SPxSteepPR::left4(int n, SPxId id)
{
   assert(thesolver->type() == SPxSolver::LEAVE);

   if (id.isValid())
   {
      Real        delta         = 0.1 + 1.0 / thesolver->basis().iteration();
      Real*       coWeights_ptr = thesolver->coWeights.get_ptr();
      const Real* workVec_ptr   = workVec.get_const_ptr();
      const Real* rhoVec        = thesolver->fVec().delta().values();
      Real        rhov_1        = 1.0 / rhoVec[n];
      Real        beta_q        = thesolver->coPvec().delta().length2() * rhov_1 * rhov_1;

#ifndef NDEBUG
      if (spxAbs(rhoVec[n]) < theeps * 0.5)
      {
         MSG_INFO3( (*thesolver->spxout), (*thesolver->spxout) << "WSTEEP04: rhoVec = "
                           << rhoVec[n] << " with smaller absolute value than 0.5*theeps = " << 0.5*theeps << std::endl; )
      }
#endif

      const IdxSet& rhoIdx = thesolver->fVec().idx();
      int           len    = thesolver->fVec().idx().size();

      for(int i = 0; i < len; ++i)
      {
         int  j = rhoIdx.index(i);
         coWeights_ptr[j] += rhoVec[j] * (beta_q * rhoVec[j] - 2.0 * rhov_1 * workVec_ptr[j]);

         if (coWeights_ptr[j] < delta)
            coWeights_ptr[j] = delta; // coWeights_ptr[j] = delta / (1+delta-x);
         else if (coWeights_ptr[j] >= infinity)
            coWeights_ptr[j] = 1.0 / theeps;
      }
      coWeights_ptr[n] = beta_q;
   }
}

Real inline computePrice(Real viol, Real weight, Real tol)
{
   if( weight < tol )
      return viol * viol / tol;
   else
      return viol * viol / weight;
}

int SPxSteepPR::buildBestPriceVectorLeave( Real feastol )
{
   int idx;
   int nsorted;
   Real x;
   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   IdxElement price;
   prices.clear();
   bestPrices.clear();

   // construct vector of all prices
   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = fTest[idx];
      if (x < -feastol)
      {
         // it might happen that we call the pricer with a tighter tolerance than what was used when computing the violations
         thesolver->isInfeasible[idx] = VIOLATED;
         price.val = computePrice(x, cpen[idx], feastol);
         price.idx = idx;
         prices.append(price);
      }
   }
   // set up structures for the quicksort implementation
   compare.elements = prices.get_const_ptr();
   // do a partial sort to move the best ones to the front
   // TODO this can be done more efficiently, since we only need the indices
   nsorted = SPxQuicksortPart(prices.get_ptr(), compare, 0, prices.size(), HYPERPRICINGSIZE);
   // copy indices of best values to bestPrices
   for( int i = 0; i < nsorted; ++i )
   {
      bestPrices.addIdx(prices[i].idx);
      thesolver->isInfeasible[prices[i].idx] = VIOLATED_AND_CHECKED;
   }

   if( nsorted > 0 )
      return prices[0].idx;
   else
      return -1;
}


int SPxSteepPR::selectLeave()
{
   assert(isConsistent());

   int retid;

   if (thesolver->hyperPricingLeave && thesolver->sparsePricingLeave)
   {
      if ( bestPrices.size() < 2 || thesolver->basis().lastUpdate() == 0 )
      {
         // call init method to build up price-vector and return index of largest price
         retid = buildBestPriceVectorLeave(theeps);
      }
      else
         retid = selectLeaveHyper(theeps);
   }
   else if (thesolver->sparsePricingLeave)
      retid = selectLeaveSparse(theeps);
   else
      retid = selectLeaveX(theeps);

   if( retid < 0 && !refined )
   {
      refined = true;
      MSG_INFO3( (*thesolver->spxout), (*thesolver->spxout) << "WSTEEP03 trying refinement step..\n"; )
      retid = selectLeaveX(theeps/STEEP_REFINETOL);
   }

   if( retid >= 0 )
   {
      assert( thesolver->coPvec().delta().isConsistent() );
      // coPvec().delta() might be not setup after the solve when it contains too many nonzeros.
      // This is intended and forcing to keep the sparsity information leads to a slowdown
      // TODO implement a dedicated solve method for unitvectors
      thesolver->basis().coSolve(thesolver->coPvec().delta(),
                                 thesolver->unitVector(retid));
      assert( thesolver->coPvec().delta().isConsistent() );
      workRhs.setup_and_assign(thesolver->coPvec().delta());
      thesolver->setup4solve(&workVec, &workRhs);
   }

   return retid;
}

int SPxSteepPR::selectLeaveX(Real tol)
{
   const Real* coWeights_ptr = thesolver->coWeights.get_const_ptr();
   const Real* fTest         = thesolver->fTest().get_const_ptr();
   Real best = -infinity;
   Real x;
   int lastIdx = -1;

   for (int i = thesolver->dim() - 1; i >= 0; --i)
   {
      x = fTest[i];

      if (x < -tol)
      {
         x = computePrice(x, coWeights_ptr[i], tol);

         if (x > best)
         {
            best = x;
            lastIdx = i;
         }
      }
   }

   return lastIdx;
}

int SPxSteepPR::selectLeaveSparse(Real tol)
{
   const Real* coWeights_ptr = thesolver->coWeights.get_const_ptr();
   const Real* fTest         = thesolver->fTest().get_const_ptr();
   Real best = -infinity;
   Real x;
   int lastIdx = -1;
   int idx;

   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = fTest[idx];

      if (x < -tol)
      {
         x = computePrice(x, coWeights_ptr[idx], tol);

         if (x > best)
         {
            best = x;
            lastIdx = idx;
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         assert(thesolver->isInfeasible[idx] == VIOLATED || thesolver->isInfeasible[idx] == VIOLATED_AND_CHECKED);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }

   return lastIdx;
}

int SPxSteepPR::selectLeaveHyper(Real tol)
{
   const Real* coPen = thesolver->coWeights.get_const_ptr();
   const Real* fTest = thesolver->fTest().get_const_ptr();
   Real leastBest = infinity;
   Real best = -infinity;
   Real x;
   int bestIdx = -1;
   int idx = 0;

   // find the best price from the short candidate list
   for( int i = bestPrices.size() - 1; i >= 0; --i )
   {
      idx = bestPrices.index(i);
      x = fTest[idx];
      if( x < -tol )
      {
         assert(thesolver->isInfeasible[idx] == VIOLATED || thesolver->isInfeasible[idx] == VIOLATED_AND_CHECKED);
         x = computePrice(x, coPen[idx], tol);

         if( x > best )
         {
            best = x;
            bestIdx = idx;
         }
         if( x < leastBest )
            leastBest = x;
      }
      else
      {
         bestPrices.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }

   // make sure we do not skip potential candidates due to a high leastBest value
   if( leastBest == infinity )
   {
      assert(bestPrices.size() == 0);
      leastBest = 0;
   }

   // scan the updated indices for a better price
   for( int i = thesolver->updateViols.size() - 1; i >= 0; --i )
   {
      idx = thesolver->updateViols.index(i);
      // is this index a candidate for bestPrices?
      if( thesolver->isInfeasible[idx] == VIOLATED )
      {
         x = fTest[idx];
         assert(x < -tol);
         x = computePrice(x, coPen[idx], tol);

         if( x > leastBest )
         {
            if( x > best )
            {
               best = x;
               bestIdx = idx;
            }
            thesolver->isInfeasible[idx] = VIOLATED_AND_CHECKED;
            bestPrices.addIdx(idx);
         }
      }
   }

   return bestIdx;
}

/* Entering Simplex
 */
void SPxSteepPR::entered4(SPxId /* id */, int n)
{
   assert(thesolver->type() == SPxSolver::ENTER);

   if (n >= 0 && n < thesolver->dim())
   {
      Real delta = 2 + 1.0 / thesolver->basis().iteration();
      Real* coWeights_ptr = thesolver->coWeights.get_ptr();
      Real* weights_ptr = thesolver->weights.get_ptr();
      const Real* workVec_ptr = workVec.get_const_ptr();
      const Real* pVec = thesolver->pVec().delta().values();
      const IdxSet& pIdx = thesolver->pVec().idx();
      const Real* coPvec = thesolver->coPvec().delta().values();
      const IdxSet& coPidx = thesolver->coPvec().idx();
      Real xi_p = 1 / thesolver->fVec().delta()[n];
      int i, j;
      Real xi_ip;

      assert(thesolver->fVec().delta()[n] > thesolver->epsilon()
              || thesolver->fVec().delta()[n] < -thesolver->epsilon());

      for (j = coPidx.size() - 1; j >= 0; --j)
      {
         i = coPidx.index(j);
         xi_ip = xi_p * coPvec[i];
         coWeights_ptr[i] += xi_ip * (xi_ip * pi_p - 2.0 * workVec_ptr[i]);
         /*
         if(coWeights_ptr[i] < 1)
             coWeights_ptr[i] = 1 / (2-x);
         */
         if (coWeights_ptr[i] < delta)
            coWeights_ptr[i] = delta;
         // coWeights_ptr[i] = 1;
         else if (coWeights_ptr[i] > infinity)
            coWeights_ptr[i] = 1 / thesolver->epsilon();
      }

      for (j = pIdx.size() - 1; j >= 0; --j)
      {
         i = pIdx.index(j);
         xi_ip = xi_p * pVec[i];
         weights_ptr[i] += xi_ip * (xi_ip * pi_p - 2.0 * (thesolver->vector(i) * workVec));
         /*
         if(weights_ptr[i] < 1)
             weights_ptr[i] = 1 / (2-x);
         */
         if (weights_ptr[i] < delta)
            weights_ptr[i] = delta;
         // weights_ptr[i] = 1;
         else if (weights_ptr[i] > infinity)
            weights_ptr[i] = 1.0 / thesolver->epsilon();
      }
   }

   /*@
       if(thesolver->isId(id))
           weights[   thesolver->number(id) ] *= 1.0001;
       else if(thesolver->isCoId(id))
           coWeights[ thesolver->number(id) ] *= 1.0001;
   */

}


SPxId SPxSteepPR::buildBestPriceVectorEnterDim( Real& best, Real feastol )
{
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coWeights_ptr = thesolver->coWeights.get_const_ptr();
   int idx;
   int nsorted;
   Real x;
   IdxElement price;

   prices.clear();
   bestPrices.clear();

   // construct vector of all prices
   for( int i = thesolver->infeasibilities.size() - 1; i >= 0; --i )
   {
      idx = thesolver->infeasibilities.index(i);
      x = coTest[idx];
      if ( x < -feastol)
      {
         // it might happen that we call the pricer with a tighter tolerance than what was used when computing the violations
         thesolver->isInfeasible[idx] = VIOLATED;
         price.val = computePrice(x, coWeights_ptr[idx], feastol);
         price.idx = idx;
         prices.append(price);
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }
   // set up structures for the quicksort implementation
   compare.elements = prices.get_const_ptr();
   // do a partial sort to move the best ones to the front
   // TODO this can be done more efficiently, since we only need the indices
   nsorted = SPxQuicksortPart(prices.get_ptr(), compare, 0, prices.size(), HYPERPRICINGSIZE);
   // copy indices of best values to bestPrices
   for( int i = 0; i < nsorted; ++i )
   {
      bestPrices.addIdx(prices[i].idx);
      thesolver->isInfeasible[prices[i].idx] = VIOLATED_AND_CHECKED;
   }

   if( nsorted > 0 )
   {
      best = prices[0].val;
      return thesolver->coId(prices[0].idx);
   }
   else
      return SPxId();
}


SPxId SPxSteepPR::buildBestPriceVectorEnterCoDim( Real& best, Real feastol )
{
   const Real* test        = thesolver->test().get_const_ptr();
   const Real* weights_ptr = thesolver->weights.get_const_ptr();
   int idx;
   int nsorted;
   Real x;
   IdxElement price;

   pricesCo.clear();
   bestPricesCo.clear();

   // construct vector of all prices
   for( int i = thesolver->infeasibilitiesCo.size() - 1; i >= 0; --i )
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      x = test[idx];
      if ( x < -feastol)
      {
         // it might happen that we call the pricer with a tighter tolerance than what was used when computing the violations
         thesolver->isInfeasibleCo[idx] = VIOLATED;
         price.val = computePrice(x, weights_ptr[idx], feastol);
         price.idx = idx;
         pricesCo.append(price);
      }
      else
      {
         thesolver->infeasibilitiesCo.remove(i);
         thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
      }
   }
   // set up structures for the quicksort implementation
   compare.elements = pricesCo.get_const_ptr();
   // do a partial sort to move the best ones to the front
   // TODO this can be done more efficiently, since we only need the indices
   nsorted = SPxQuicksortPart(pricesCo.get_ptr(), compare, 0, pricesCo.size(), HYPERPRICINGSIZE);
   // copy indices of best values to bestPrices
   for( int i = 0; i < nsorted; ++i )
   {
      bestPricesCo.addIdx(pricesCo[i].idx);
      thesolver->isInfeasibleCo[pricesCo[i].idx] = VIOLATED_AND_CHECKED;
   }

   if( nsorted > 0 )
   {
      best = pricesCo[0].val;
      return thesolver->id(pricesCo[0].idx);
   }
   else
      return SPxId();
}


SPxId SPxSteepPR::selectEnter()
{
   assert(thesolver != 0);
   SPxId enterId;

   enterId = selectEnterX(theeps);

   if( !enterId.isValid() && !refined )
   {
      refined = true;
      MSG_INFO3( (*thesolver->spxout), (*thesolver->spxout) << "WSTEEP05 trying refinement step..\n"; )
      enterId = selectEnterX(theeps/STEEP_REFINETOL);
   }

   assert(isConsistent());

   if (enterId.isValid())
   {
      SSVector& delta = thesolver->fVec().delta();

      thesolver->basis().solve4update(delta, thesolver->vector(enterId));

      workRhs.setup_and_assign(delta);
      pi_p = 1 + delta.length2();

      thesolver->setup4coSolve(&workVec, &workRhs);
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterX(Real tol)
{
   SPxId enterId;
   SPxId enterCoId;
   Real best;
   Real bestCo;

   best = -infinity;
   bestCo = -infinity;

   if( thesolver->hyperPricingEnter && !refined )
   {
      if( bestPrices.size() < 2 || thesolver->basis().lastUpdate() == 0 )
         enterCoId = (thesolver->sparsePricingEnter) ? buildBestPriceVectorEnterDim(best, tol) : selectEnterDenseDim(best, tol);
      else
         enterCoId = (thesolver->sparsePricingEnter) ? selectEnterHyperDim(best, tol) : selectEnterDenseDim(best, tol);

      if( bestPricesCo.size() < 2 || thesolver->basis().lastUpdate() == 0 )
         enterId = (thesolver->sparsePricingEnterCo) ? buildBestPriceVectorEnterCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
      else
         enterId = (thesolver->sparsePricingEnterCo) ? selectEnterHyperCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
   }
   else
   {
      enterCoId = (thesolver->sparsePricingEnter && !refined) ? selectEnterSparseDim(best, tol) : selectEnterDenseDim(best, tol);
      enterId = (thesolver->sparsePricingEnterCo && !refined) ? selectEnterSparseCoDim(bestCo, tol) : selectEnterDenseCoDim(bestCo, tol);
   }

   // prefer slack indices to reduce nonzeros in basis matrix
   if( enterCoId.isValid() && (best > SPARSITY_TRADEOFF * bestCo || !enterId.isValid()) )
      return enterCoId;
   else
      return enterId;
}


SPxId SPxSteepPR::selectEnterHyperDim(Real& best, Real tol)
{
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coWeights_ptr = thesolver->coWeights.get_const_ptr();

   Real leastBest = infinity;
   Real x;
   int enterIdx = -1;
   int idx;

   // find the best price from short candidate list
   for( int i = bestPrices.size() - 1; i >= 0; --i )
   {
      idx = bestPrices.index(i);
      x = coTest[idx];
      if( x < -tol )
      {
         x = computePrice(x, coWeights_ptr[idx], tol);
         if( x > best )
         {
            best = x;
            enterIdx = idx;
         }
         if( x < leastBest )
            leastBest = x;
      }
      else
      {
         bestPrices.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }

   // make sure we do not skip potential candidates due to a high leastBest value
   if( leastBest == infinity )
   {
      assert(bestPrices.size() == 0);
      leastBest = 0;
   }

   // scan the updated indices for a better price
   for( int i = thesolver->updateViols.size() -1; i >= 0; --i )
   {
      idx = thesolver->updateViols.index(i);
      // only look at indices that were not checked already
      if( thesolver->isInfeasible[idx] == VIOLATED )
      {
         x = coTest[idx];
         if( x < -tol )
         {
            x = computePrice(x, coWeights_ptr[idx], tol);
            if( x > leastBest )
            {
               if (x > best)
               {
                  best = x;
                  enterIdx = idx;
               }
               // put index into candidate list
               thesolver->isInfeasible[idx] = VIOLATED_AND_CHECKED;
               bestPrices.addIdx(idx);
            }
         }
         else
         {
            thesolver->isInfeasible[idx] = NOT_VIOLATED;
         }
      }
   }

   if( enterIdx >= 0 )
      return thesolver->coId(enterIdx);
   else
      return SPxId();
}


SPxId SPxSteepPR::selectEnterHyperCoDim(Real& best, Real tol)
{
   const Real* test        = thesolver->test().get_const_ptr();
   const Real* weights_ptr = thesolver->weights.get_const_ptr();

   Real leastBest = infinity;
   Real x;
   int enterIdx = -1;
   int idx;

   // find the best price from short candidate list
   for( int i = bestPricesCo.size() - 1; i >= 0; --i )
   {
      idx = bestPricesCo.index(i);
      x = test[idx];
      if( x < -tol )
      {
         x = computePrice(x, weights_ptr[idx], tol);
         if( x > best )
         {
            best = x;
            enterIdx = idx;
         }
         if( x < leastBest )
            leastBest = x;
      }
      else
      {
         bestPricesCo.remove(i);
         thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
      }
   }

   // make sure we do not skip potential candidates due to a high leastBest value
   if( leastBest == infinity )
   {
      assert(bestPricesCo.size() == 0);
      leastBest = 0;
   }

   // scan the updated indices for a better price
   for( int i = thesolver->updateViolsCo.size() -1; i >= 0; --i )
   {
      idx = thesolver->updateViolsCo.index(i);
      // only look at indices that were not checked already
      if( thesolver->isInfeasibleCo[idx] == VIOLATED )
      {
         x = test[idx];
         if( x < -tol )
         {
            x = computePrice(x, weights_ptr[idx], tol);
            if( x > leastBest )
            {
               if (x > best)
               {
                  best = x;
                  enterIdx = idx;
               }
               // put index into candidate list
               thesolver->isInfeasibleCo[idx] = VIOLATED_AND_CHECKED;
               bestPricesCo.addIdx(idx);
            }
         }
         else
         {
            thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
         }
      }
   }

   if( enterIdx >= 0 )
      return thesolver->id(enterIdx);
   else
      return SPxId();
}


SPxId SPxSteepPR::selectEnterSparseDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coWeights_ptr = thesolver->coWeights.get_const_ptr();

   int idx;
   Real x;

   for (int i = thesolver->infeasibilities.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = coTest[idx];

      if (x < -tol)
      {
         x = computePrice(x, coWeights_ptr[idx], tol);
         if (x > best)
         {
            best = x;
            enterId = thesolver->coId(idx);
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterSparseCoDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* test          = thesolver->test().get_const_ptr();
   const Real* weights_ptr   = thesolver->weights.get_const_ptr();

   int idx;
   Real x;

   for (int i = thesolver->infeasibilitiesCo.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      x = test[idx];

      if (x < -tol)
      {
         x = computePrice(x, weights_ptr[idx], tol);
         if (x > best)
         {
            best   = x;
            enterId = thesolver->id(idx);
         }
      }
      else
      {
         thesolver->infeasibilitiesCo.remove(i);
         thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
      }
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterDenseDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* coTest        = thesolver->coTest().get_const_ptr();
   const Real* coWeights_ptr = thesolver->coWeights.get_const_ptr();

   Real x;

   for (int i = 0, end = thesolver->dim(); i < end; ++i)
   {
      x = coTest[i];
      if (x < -tol)
      {
         x = computePrice(x, coWeights_ptr[i], tol);
         if (x > best)
         {
            best = x;
            enterId = thesolver->coId(i);
         }
      }
   }
   return enterId;
}

SPxId SPxSteepPR::selectEnterDenseCoDim(Real& best, Real tol)
{
   SPxId enterId;
   const Real* test          = thesolver->test().get_const_ptr();
   const Real* weights_ptr   = thesolver->weights.get_const_ptr();

   Real x;

   for(int i = 0, end = thesolver->coDim(); i < end; ++i)
   {
      x = test[i];
      if (x < -tol)
      {
         x = computePrice(x, weights_ptr[i], tol);
         if (x > best)
         {
            best   = x;
            enterId = thesolver->id(i);
         }
      }
   }
   return enterId;
}


void SPxSteepPR::addedVecs(int n)
{
   DVector& weights = thesolver->weights;
   n = weights.dim();
   weights.reDim(thesolver->coDim());

   if (thesolver->type() == SPxSolver::ENTER)
   {
      for (; n < weights.dim(); ++n)
         weights[n] = 2;
   }
}

void SPxSteepPR::addedCoVecs(int n)
{
   DVector& coWeights = thesolver->coWeights;
   n = coWeights.dim();
   workVec.reDim (thesolver->dim());
   coWeights.reDim (thesolver->dim());
   for (; n < coWeights.dim(); ++n)
      coWeights[n] = 1;
}

void SPxSteepPR::removedVec(int i)
{
   assert(thesolver != 0);
   DVector& weights = thesolver->weights;
   weights[i] = weights[weights.dim()];
   weights.reDim(thesolver->coDim());
}

void SPxSteepPR::removedVecs(const int perm[])
{
   assert(thesolver != 0);
   DVector& weights = thesolver->weights;
   if (thesolver->type() == SPxSolver::ENTER)
   {
      int i;
      int j = weights.dim();
      for (i = 0; i < j; ++i)
      {
         if (perm[i] >= 0)
            weights[perm[i]] = weights[i];
      }
   }
   weights.reDim(thesolver->coDim());
}

void SPxSteepPR::removedCoVec(int i)
{
   assert(thesolver != 0);
   DVector& coWeights = thesolver->coWeights;
   coWeights[i] = coWeights[coWeights.dim()];
   coWeights.reDim(thesolver->dim());
}

void SPxSteepPR::removedCoVecs(const int perm[])
{
   assert(thesolver != 0);
   DVector& coWeights = thesolver->coWeights;
   int i;
   int j = coWeights.dim();
   for (i = 0; i < j; ++i)
   {
      if (perm[i] >= 0)
         coWeights[perm[i]] = coWeights[i];
   }
   coWeights.reDim(thesolver->dim());
}

bool SPxSteepPR::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   DVector& w = thesolver->weights;
   DVector& coW = thesolver->coWeights;
   if (thesolver != 0 && thesolver->type() == SPxSolver::LEAVE && setup == EXACT)
   {
      int i;
      SSVector tmp(thesolver->dim(), thesolver->epsilon());
      Real x;
      for (i = thesolver->dim() - 1; i >= 0; --i)
      {
         thesolver->basis().coSolve(tmp, thesolver->unitVector(i));
         x = coW[i] - tmp.length2();
         if (x > thesolver->leavetol() || -x > thesolver->leavetol())
         {
            MSG_ERROR( std::cerr << "ESTEEP03 x[" << i << "] = " << x << std::endl; )
         }
      }
   }

   if (thesolver != 0 && thesolver->type() == SPxSolver::ENTER)
   {
      int i;
      for (i = thesolver->dim() - 1; i >= 0; --i)
      {
         if (coW[i] < thesolver->epsilon())
            return MSGinconsistent("SPxSteepPR");
      }

      for (i = thesolver->coDim() - 1; i >= 0; --i)
      {
         if (w[i] < thesolver->epsilon())
            return MSGinconsistent("SPxSteepPR");
      }
   }
#endif

   return true;
}
} // namespace soplex
