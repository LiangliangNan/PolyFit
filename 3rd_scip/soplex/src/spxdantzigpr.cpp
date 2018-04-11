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

// #define EQ_PREF 1000 

#include "spxdefines.h"
#include "spxdantzigpr.h"

namespace soplex
{

int SPxDantzigPR::selectLeave()
{
   assert(thesolver != 0);

   if( thesolver->sparsePricingLeave )
      return selectLeaveSparse();

   //    const Real* up  = thesolver->ubBound();
   //    const Real* low = thesolver->lbBound();

   Real best = -theeps;
   int  n    = -1;

   for(int i = thesolver->dim() - 1; i >= 0; --i)
   {
      Real x = thesolver->fTest()[i];

      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (up[i] == low[i]));
         if (x < best)
         {
            n    = i;
            best = x;
         }
      }
   }
   return n;
}

int SPxDantzigPR::selectLeaveSparse()
{
   assert(thesolver != 0);

   Real best   = -theeps;
   Real x;
   int  n      = -1;
   int  index;

   for(int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      index = thesolver->infeasibilities.index(i);
      x = thesolver->fTest()[index];
      if (x < -theeps)
      {
         if (x < best)
         {
            n    = index;
            best = x;
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);
         assert(thesolver->isInfeasible[index] > 0);
         thesolver->isInfeasible[index] = 0;
      }
   }
   return n;
}


SPxId SPxDantzigPR::selectEnter()
{
   assert(thesolver != 0);

   // const SPxBasis::Desc&    ds   = thesolver->basis().desc();

   SPxId enterId;
   enterId = selectEnterX();

   return enterId;
}

SPxId SPxDantzigPR::selectEnterX()
{
   SPxId enterId;
   SPxId enterIdCo;
   Real best;
   Real bestCo;

   best = -theeps;
   bestCo = -theeps;
   enterId = (thesolver->sparsePricingEnter) ? selectEnterSparseDim(best,enterId) : selectEnterDenseDim(best,enterId);
   enterIdCo = (thesolver->sparsePricingEnterCo) ? selectEnterSparseCoDim(bestCo,enterId) : selectEnterDenseCoDim(bestCo,enterId);

   // prefer slack indices to reduce nonzeros in basis matrix
   if( enterId.isValid() && (best > SPARSITY_TRADEOFF * bestCo || !enterIdCo.isValid()) )
      return enterId;
   else
      return enterIdCo;
}


SPxId SPxDantzigPR::selectEnterSparseDim(Real& best,SPxId& enterId)
{
   assert(thesolver != 0);

   int idx;
   Real x;

   for (int i = thesolver->infeasibilities.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilities.index(i);
      x = thesolver->coTest()[idx];

      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasis::Desc::P_FREE
         //                || ds.coStatus(i) == SPxBasis::Desc::D_FREE));
         if (x < best)
         {
            enterId = thesolver->coId(idx);
            best = x;
         }
      }
      else
      {
         thesolver->infeasibilities.remove(i);

         assert(thesolver->isInfeasible[idx]);
         thesolver->isInfeasible[idx] = 0;
      }
   }
   return enterId;
}

SPxId SPxDantzigPR::selectEnterSparseCoDim(Real& best,SPxId& enterId)
{
   assert(thesolver != 0);

   int idx;
   Real x;

   for (int i = thesolver->infeasibilitiesCo.size() - 1; i >= 0; --i)
   {
      idx = thesolver->infeasibilitiesCo.index(i);
      x = thesolver->test()[idx];

      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasis::Desc::P_FREE
         //                || ds.coStatus(i) == SPxBasis::Desc::D_FREE));
         if (x < best)
         {
            enterId = thesolver->id(idx);
            best = x;
         }
      }
      else
      {
         thesolver->infeasibilitiesCo.remove(i);
         assert(thesolver->isInfeasibleCo[idx] > 0);
         thesolver->isInfeasibleCo[idx] = 0;
      }
   }
   return enterId;
}

SPxId SPxDantzigPR::selectEnterDenseDim(Real& best,SPxId& enterId)
{
   assert(thesolver != 0);

   Real x;

   for (int i = thesolver->dim() - 1; i >= 0; --i)
   {
      x = thesolver->coTest()[i];

      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasis::Desc::P_FREE
         //                || ds.coStatus(i) == SPxBasis::Desc::D_FREE));
         if (x < best)
         {
            enterId   = thesolver->coId(i);
            best = x;
         }
      }
   }
   return enterId;
}

SPxId SPxDantzigPR::selectEnterDenseCoDim(Real& best,SPxId& enterId)
{
   assert(thesolver != 0);

   Real x;

   for (int i = thesolver->coDim() - 1; i >= 0; --i)
   {
      x = thesolver->test()[i];

      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (ds.status(i) == SPxBasis::Desc::P_FREE
         //                || ds.status(i) == SPxBasis::Desc::D_FREE));
         if (x < best)
         {
            enterId   = thesolver->id(i);
            best = x;
         }
      }
   }
   return enterId;
}

} // namespace soplex
