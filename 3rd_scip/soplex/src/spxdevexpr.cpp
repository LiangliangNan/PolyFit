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

#include "spxdefines.h"
#include "spxdevexpr.h"

#define DEVEX_REFINETOL 2.0

namespace soplex
{

void SPxDevexPR::load(SPxSolver* base)
{
   thesolver = base;
   setRep(base->rep());
   assert(isConsistent());
}

bool SPxDevexPR::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (thesolver != 0)
      if (weights.dim() != thesolver->coDim()
           || coWeights.dim() != thesolver->dim())
         return MSGinconsistent("SPxDevexPR");
#endif

   return true;
}

void SPxDevexPR::setupWeights(SPxSolver::Type tp)
{
   int i;
   int coWeightSize = 0;
   int weightSize = 0;

   DVector& weights = thesolver->weights;
   DVector& coWeights = thesolver->coWeights;

   if( tp == SPxSolver::ENTER )
   {
      coWeights.reDim(thesolver->dim(), false);
      for( i = thesolver->dim() - 1; i >= coWeightSize; --i )
         coWeights[i] = 2.0;

      weights.reDim(thesolver->coDim(), false);
      for( i = thesolver->coDim() - 1; i >= weightSize; --i )
         weights[i] = 2.0;
   }
   else
   {
      coWeights.reDim(thesolver->dim(), false);
      for( i = thesolver->dim() - 1; i >= coWeightSize; --i )
         coWeights[i] = 1.0;
   }
   thesolver->weightsAreSetup = true;
}

void SPxDevexPR::setType(SPxSolver::Type tp)
{
   setupWeights(tp);
   refined = false;

   bestPrices.clear();
   bestPrices.setMax(thesolver->dim());
   prices.reMax(thesolver->dim());

   if( tp == SPxSolver::ENTER )
   {
      bestPricesCo.clear();
      bestPricesCo.setMax(thesolver->coDim());
      pricesCo.reMax(thesolver->coDim());
   }

   assert(isConsistent());
}

/**@todo suspicious: Shouldn't the relation between dim, coDim, Vecs, 
 *       and CoVecs be influenced by the representation ?
 */
void SPxDevexPR::setRep(SPxSolver::Representation)
{
   if (thesolver != 0)
   {
      // resize weights and initialize new entries
      addedVecs(thesolver->coDim());
      addedCoVecs(thesolver->dim());
      assert(isConsistent());
   }
}

Real inline computePrice(Real viol, Real weight, Real tol)
{
   if( weight < tol )
      return viol * viol / tol;
   else
      return viol * viol / weight;
}


int SPxDevexPR::buildBestPriceVectorLeave( Real feastol )
{
   int idx;
   int nsorted;
   Real fTesti;
   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   IdxElement price;
   prices.clear();
   bestPrices.clear();

   // TODO we should check infeasiblities for duplicates or loop over dimension
   //      bestPrices may then also contain duplicates!
   // construct vector of all prices
   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      fTesti = fTest[idx];
      if (fTesti < -feastol)
      {
         thesolver->isInfeasible[idx] = VIOLATED;
         price.idx = idx;
         price.val = computePrice(fTesti, cpen[idx], feastol);
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

int SPxDevexPR::selectLeave()
{
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

   if ( retid < 0 && !refined )
   {
      refined = true;
      MSG_INFO3( (*thesolver->spxout), (*thesolver->spxout) << "WDEVEX02 trying refinement step..\n"; )
      retid = selectLeaveX(theeps/DEVEX_REFINETOL);
   }

   assert(retid < thesolver->dim());

   return retid;
}

int SPxDevexPR::selectLeaveX(Real feastol, int start, int incr)
{
   Real x;

   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   Real best = 0;
   int bstI = -1;
   int end = thesolver->coWeights.dim();

   for (; start < end; start += incr)
   {
      if (fTest[start] < -feastol)
      {
         x = computePrice(fTest[start], cpen[start], feastol);
         if (x > best)
         {
            best = x;
            bstI = start;
            last = cpen[start];
         }
      }
   }
   return bstI;
}

int SPxDevexPR::selectLeaveSparse(Real feastol)
{
   Real x;

   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   Real best = 0;
   int bstI = -1;
   int idx = -1;

   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = fTest[idx];
      if (x < -feastol)
      {
         x = computePrice(x, cpen[idx], feastol);
         if (x > best)
         {
            best = x;
            bstI = idx;
            last = cpen[idx];
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         assert(thesolver->isInfeasible[idx] == VIOLATED || thesolver->isInfeasible[idx] == VIOLATED_AND_CHECKED);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }
   return bstI;
}

int SPxDevexPR::selectLeaveHyper(Real feastol)
{
   Real x;

   const Real* fTest = thesolver->fTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   Real best = 0;
   Real leastBest = infinity;
   int bstI = -1;
   int idx = -1;

   // find the best price from the short candidate list
   for( int i = bestPrices.size() - 1; i >= 0; --i )
   {
      idx = bestPrices.index(i);
      x = fTest[idx];
      if( x < -feastol )
      {
         x = computePrice(x, cpen[idx], feastol);
         if( x > best )
         {
            best = x;
            bstI = idx;
            last = cpen[idx];
         }
         // get the smallest price of candidate list
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
      // only look at indeces that were not checked already
      if( thesolver->isInfeasible[idx] == VIOLATED )
      {
         x = fTest[idx];
         assert(x < -feastol);
         x = computePrice(x, cpen[idx], feastol);
         if( x > leastBest )
         {
            if( x > best )
            {
               best = x;
               bstI = idx;
               last = cpen[idx];
            }
            // put index into candidate list
            thesolver->isInfeasible[idx] = VIOLATED_AND_CHECKED;
            bestPrices.addIdx(idx);
         }
      }
   }

   return bstI;
}

void SPxDevexPR::left4(int n, SPxId id)
{
   DVector& coWeights = thesolver->coWeights;
   if (id.isValid())
   {
      int i, j;
      Real x;
      const Real* rhoVec = thesolver->fVec().delta().values();
      Real rhov_1 = 1 / rhoVec[n];
      Real beta_q = thesolver->coPvec().delta().length2() * rhov_1 * rhov_1;

#ifndef NDEBUG
      if (spxAbs(rhoVec[n]) < theeps)
      {
         MSG_INFO3( (*thesolver->spxout), (*thesolver->spxout) << "WDEVEX01: rhoVec = "
                           << rhoVec[n] << " with smaller absolute value than theeps = " << theeps << std::endl; )
      }
#endif  // NDEBUG

      //  Update #coPenalty# vector
      const IdxSet& rhoIdx = thesolver->fVec().idx();
      int len = thesolver->fVec().idx().size();
      for (i = len - 1; i >= 0; --i)
      {
         j = rhoIdx.index(i);
         x = rhoVec[j] * rhoVec[j] * beta_q;
         // if(x > coPenalty[j])
         coWeights[j] += x;
      }

      coWeights[n] = beta_q;
   }
}

SPxId SPxDevexPR::buildBestPriceVectorEnterDim( Real& best, Real feastol )
{
   int idx;
   int nsorted;
   Real x;
   const Real* coTest = thesolver->coTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   IdxElement price;
   prices.clear();
   bestPrices.clear();

   // construct vector of all prices
   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = coTest[idx];
      if ( x < -feastol)
      {
         thesolver->isInfeasible[idx] = VIOLATED;
         price.idx = idx;
         price.val = computePrice(x, cpen[idx], feastol);
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

SPxId SPxDevexPR::buildBestPriceVectorEnterCoDim( Real& best, Real feastol )
{
   int idx;
   int nsorted;
   Real x;
   const Real* test = thesolver->test().get_const_ptr();
   const Real* pen = thesolver->weights.get_const_ptr();
   IdxElement price;
   pricesCo.clear();
   bestPricesCo.clear();

   // construct vector of all prices
   for (int i = thesolver->infeasibilitiesCo.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      x = test[idx];
      if ( x < -feastol)
      {
         thesolver->isInfeasibleCo[idx] = VIOLATED;
         price.idx = idx;
         price.val = computePrice(x, pen[idx], feastol);
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

SPxId SPxDevexPR::selectEnter()
{
   assert(thesolver != 0);

   SPxId enterId;

   enterId = selectEnterX(theeps);

   if( !enterId.isValid() && !refined )
   {
      refined = true;
      MSG_INFO3( (*thesolver->spxout), (*thesolver->spxout) << "WDEVEX02 trying refinement step..\n"; )
      enterId = selectEnterX(theeps/DEVEX_REFINETOL);
   }

   return enterId;
}

// choose the best entering index among columns and rows but prefer sparsity
SPxId SPxDevexPR::selectEnterX(Real tol)
{
   SPxId enterId;
   SPxId enterCoId;
   Real best;
   Real bestCo;

   best = 0;
   bestCo = 0;
   last = 1.0;

   // avoid uninitialized value later on in entered4X()
   last = 1.0;

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

   // prefer coIds to increase the number of unit vectors in the basis matrix, i.e., rows in colrep and cols in rowrep
   if( enterCoId.isValid() && (best > SPARSITY_TRADEOFF * bestCo || !enterId.isValid()) )
      return enterCoId;
   else
      return enterId;
}

SPxId SPxDevexPR::selectEnterHyperDim(Real& best, Real feastol)
{
   const Real* cTest = thesolver->coTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   Real leastBest = infinity;
   Real x;
   int enterIdx = -1;
   int idx;

   // find the best price from short candidate list
   for( int i = bestPrices.size() - 1; i >= 0; --i )
   {
      idx = bestPrices.index(i);
      x = cTest[idx];
      if( x < -feastol )
      {
         x = computePrice(x, cpen[idx], feastol);
         if( x > best )
         {
            best = x;
            enterIdx = idx;
            last = cpen[idx];
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

   // scan the updated indeces for a better price
   for( int i = thesolver->updateViols.size() -1; i >= 0; --i )
   {
      idx = thesolver->updateViols.index(i);
      // only look at indeces that were not checked already
      if( thesolver->isInfeasible[idx] == VIOLATED )
      {
         x = cTest[idx];
         if( x < -feastol )
         {
            x = computePrice(x, cpen[idx], feastol);
            if(x > leastBest)
            {
               if( x > best )
               {
                  best = x;
                  enterIdx = idx;
                  last = cpen[idx];
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

   if (enterIdx >= 0)
      return thesolver->coId(enterIdx);
   else
      return SPxId();
}


SPxId SPxDevexPR::selectEnterHyperCoDim(Real& best, Real feastol)
{
   const Real* test = thesolver->test().get_const_ptr();
   const Real* pen = thesolver->weights.get_const_ptr();
   Real leastBest = infinity;
   Real x;
   int enterIdx = -1;
   int idx;

   // find the best price from short candidate list
   for( int i = bestPricesCo.size() - 1; i >= 0; --i )
   {
      idx = bestPricesCo.index(i);
      x = test[idx];
      if( x < -feastol )
      {
         x = computePrice(x, pen[idx], feastol);
         if( x > best )
         {
            best = x;
            enterIdx = idx;
            last = pen[idx];
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

   //scan the updated indeces for a better price
   for( int i = thesolver->updateViolsCo.size() -1; i >= 0; --i )
   {
      idx = thesolver->updateViolsCo.index(i);
      // only look at indeces that were not checked already
      if( thesolver->isInfeasibleCo[idx] == VIOLATED )
      {
         x = test[idx];
         if( x < -feastol )
         {
            x = computePrice(x, pen[idx], feastol);
            if( x > leastBest )
            {
               if( x > best )
               {
                  best = x;
                  enterIdx = idx;
                  last = pen[idx];
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


SPxId SPxDevexPR::selectEnterSparseDim(Real& best, Real feastol)
{
   const Real* cTest = thesolver->coTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   int enterIdx = -1;
   int idx;
   Real x;

   assert(thesolver->coWeights.dim() == thesolver->coTest().dim());
   for(int i = thesolver->infeasibilities.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = cTest[idx];
      if( x < -feastol )
      {
         x = computePrice(x, cpen[idx], feastol);
         if (x > best)
         {
            best = x;
            enterIdx = idx;
            last = cpen[idx];
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         thesolver->isInfeasible[idx] = NOT_VIOLATED;
      }
   }
   if (enterIdx >= 0)
      return thesolver->coId(enterIdx);

   return SPxId();
}


SPxId SPxDevexPR::selectEnterSparseCoDim(Real& best, Real feastol)
{
   const Real* test = thesolver->test().get_const_ptr();
   const Real* pen = thesolver->weights.get_const_ptr();
   int enterIdx = -1;
   int idx;
   Real x;

   assert(thesolver->weights.dim() == thesolver->test().dim());
   for (int i = thesolver->infeasibilitiesCo.size() -1; i >= 0; --i)
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      x = test[idx];
      if (x < -feastol)
      {
         x = computePrice(x, pen[idx], feastol);
         if (x > best)
         {
            best = x;
            enterIdx = idx;
            last = pen[idx];
         }
      }
      else
      {
         thesolver->infeasibilitiesCo.remove(i);
         thesolver->isInfeasibleCo[idx] = NOT_VIOLATED;
      }
   }

   if (enterIdx >= 0)
      return thesolver->id(enterIdx);

   return SPxId();
}


SPxId SPxDevexPR::selectEnterDenseDim(Real& best, Real feastol, int start, int incr)
{
   const Real* cTest = thesolver->coTest().get_const_ptr();
   const Real* cpen = thesolver->coWeights.get_const_ptr();
   int end = thesolver->coWeights.dim();
   int enterIdx = -1;
   Real x;

   assert(end == thesolver->coTest().dim());
   for (; start < end; start += incr)
   {
      x = cTest[start];
      if( x < -feastol )
      {
         x = computePrice(x, cpen[start], feastol);
         if (x > best)
         {
            best = x;
            enterIdx = start;
            last = cpen[start];
         }
      }
   }

   if (enterIdx >= 0)
      return thesolver->coId(enterIdx);

   return SPxId();
}


SPxId SPxDevexPR::selectEnterDenseCoDim(Real& best, Real feastol, int start, int incr)
{
   const Real* test = thesolver->test().get_const_ptr();
   const Real* pen = thesolver->weights.get_const_ptr();
   int end = thesolver->weights.dim();
   int enterIdx = -1;
   Real x;

   assert(end == thesolver->test().dim());
   for (; start < end; start += incr)
   {
      x = test[start];
      if (test[start] < -feastol)
      {
         x = computePrice(x, pen[start], feastol);
         if (x > best)
         {
            best = x;
            enterIdx = start;
            last = pen[start];
         }
      }
   }

   if (enterIdx >= 0)
      return thesolver->id(enterIdx);

   return SPxId();
}


/**@todo suspicious: the pricer should be informed, that variable id 
    has entered the basis at position n, but the id is not used here 
    (this is true for all pricers)
*/
void SPxDevexPR::entered4(SPxId /*id*/, int n)
{
   DVector& weights = thesolver->weights;
   DVector& coWeights = thesolver->coWeights;

   if (n >= 0 && n < thesolver->dim())
   {
      const Real* pVec = thesolver->pVec().delta().values();
      const IdxSet& pIdx = thesolver->pVec().idx();
      const Real* coPvec = thesolver->coPvec().delta().values();
      const IdxSet& coPidx = thesolver->coPvec().idx();
      Real xi_p = 1 / thesolver->fVec().delta()[n];
      int i, j;

      assert(thesolver->fVec().delta()[n] > thesolver->epsilon()
              || thesolver->fVec().delta()[n] < -thesolver->epsilon());

      xi_p = xi_p * xi_p * last;

      for (j = coPidx.size() - 1; j >= 0; --j)
      {
         i = coPidx.index(j);
         coWeights[i] += xi_p * coPvec[i] * coPvec[i];
         if (coWeights[i] <= 1 || coWeights[i] > 1e+6)
         {
            setupWeights(SPxSolver::ENTER);
            return;
         }
      }

      for (j = pIdx.size() - 1; j >= 0; --j)
      {
         i = pIdx.index(j);
         weights[i] += xi_p * pVec[i] * pVec[i];
         if (weights[i] <= 1 || weights[i] > 1e+6)
         {
            setupWeights(SPxSolver::ENTER);
            return;
         }
      }
   }
}

void SPxDevexPR::addedVecs (int n)
{
   int initval = (thesolver->type() == SPxSolver::ENTER) ? 2 : 1;
   DVector& weights = thesolver->weights;
   n = weights.dim();
   weights.reDim (thesolver->coDim());
   for( int i = weights.dim()-1; i >= n; --i )
      weights[i] = initval;
}

void SPxDevexPR::addedCoVecs(int n)
{
   int initval = (thesolver->type() == SPxSolver::ENTER) ? 2 : 1;
   DVector& coWeights = thesolver->coWeights;
   n = coWeights.dim();
   coWeights.reDim(thesolver->dim());
   for( int i = coWeights.dim()-1; i >= n; --i )
      coWeights[i] = initval;
}

} // namespace soplex
