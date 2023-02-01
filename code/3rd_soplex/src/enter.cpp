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

/*      \SubSection{Updating the Basis for Entering Variables}
 */
#include <assert.h>

#include "spxdefines.h"
#include "spxratiotester.h"
#include "spxpricer.h"
#include "spxout.h"
#include "exceptions.h"

namespace soplex
{

/*
In the entering simplex algorithms (i.e. iteratively a vector is selected to
\em enter the simplex basis as in the dual rowwise and primal columnwise case)
let $A$ denote the current basis, $x$ and entering vector and $f$ the
feasibility vector. For a feasible basis $l \le f \le u$ holds.
For the rowwise case $f$ is obtained by solving $f^T = c^T A^{-1}$,
wherease in columnwisecase $f = A^{-1} b$.

Let us further consider the rowwise case. Exchanging $x$ with the $i$-th
vector of $A$ yields

\begin{equation}\label{update.eq}
    A^{(i)} = E_i A \hbox{, with } E_i = I + e_i (x^T A^{-1} - e_i^T).
\end{equation}

With $E_i^{-1} = I + e_i \frac{e_i^T - \delta^T}{\delta}$,
$\delta^T = x^T A^{-1}$ one gets the new feasibility vector

\begin{eqnarray*}
        (f^{(i)})^T
    &=& c^T (A^{(i)})^{-1}      \\
    &=& c^T A^{-1} + c^T A^{-1} e_i \frac{e_i^T - \delta^T}{\delta_i} \\
    &=& f^T + \frac{f_i}{\delta_i} e_i^T - \frac{f_i}{\delta_i} \delta^T. \\
\end{eqnarray*}

The selection of the leaving vector $i^*$ for the basis must ensure, that for
all $j \ne i^*$ $f^{(i^*)}_j$ remains within its bounds $l_j$ and $u_j$.
 */


/*
    Testing all values of |pVec| against its bounds. If $i$, say, is violated
    the violation is saved as negative value in |theTest[i]|.
 */
Real SPxSolver::test(int i, SPxBasis::Desc::Status stat) const
{
   assert(type() == ENTER);
   assert(!isBasic(stat));

   Real x;

   switch (stat)
   {
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_BOTH:
      assert(rep() == ROW);
      x = (*thePvec)[i] - lhs(i);
      if (x < 0)
         return x;
      // no break: next is else case
      //lint -fallthrough
   case SPxBasis::Desc::D_ON_LOWER:
      assert(rep() == ROW);
      return rhs(i) - (*thePvec)[i];
   case SPxBasis::Desc::D_ON_UPPER:
      assert(rep() == ROW);
      return (*thePvec)[i] - lhs(i);

   case SPxBasis::Desc::P_ON_UPPER:
      assert(rep() == COLUMN);
      return maxObj(i) - (*thePvec)[i];
   case SPxBasis::Desc::P_ON_LOWER:
      assert(rep() == COLUMN);
      return (*thePvec)[i] - maxObj(i);
   case SPxBasis::Desc::P_FREE :
      x = maxObj(i) - (*thePvec)[i];
      return (x < 0) ? x : -x;

   default:
      return 0;
   }
}

void SPxSolver::computeTest()
{

   const SPxBasis::Desc& ds = desc();
   Real pricingTol = leavetol();
   m_pricingViolCoUpToDate = true;
   m_pricingViolCo = 0;
   infeasibilitiesCo.clear();
   int ninfeasibilities = 0;
   int sparsitythreshold = (int) (sparsePricingFactor * coDim());

   for(int i = 0; i < coDim(); ++i)
   {
      SPxBasis::Desc::Status stat = ds.status(i);

      if(isBasic(stat))
      {
         theTest[i] = 0.0;
         if( remainingRoundsEnterCo == 0 )
            isInfeasibleCo[i] = SPxPricer::NOT_VIOLATED;
      }
      else
      {
         assert(!isBasic(stat));
         theTest[i] = test(i, stat);

         if( remainingRoundsEnterCo == 0 )
         {
            if( theTest[i] < -pricingTol )
            {
               assert(infeasibilitiesCo.size() < infeasibilitiesCo.max());
               m_pricingViolCo -= theTest[i];
               infeasibilitiesCo.addIdx(i);
               isInfeasibleCo[i] = SPxPricer::VIOLATED;
               ++ninfeasibilities;
            }
            else
               isInfeasibleCo[i] = SPxPricer::NOT_VIOLATED;
            if( ninfeasibilities > sparsitythreshold)
            {
               MSG_INFO2( (*spxout), (*spxout) << " --- using dense pricing"
                                 << std::endl; )
               remainingRoundsEnterCo = DENSEROUNDS;
               sparsePricingEnterCo = false;
               ninfeasibilities = 0;
            }
         }
         else if( theTest[i] < -pricingTol )
            m_pricingViolCo -= theTest[i];
      }
   }
   if( ninfeasibilities == 0 && !sparsePricingEnterCo )
      --remainingRoundsEnterCo;
   else if( ninfeasibilities <= sparsitythreshold && !sparsePricingEnterCo )
   {
      MSG_INFO2( (*spxout),
         std::streamsize prec = spxout->precision();
         if( hyperPricingEnter )
            (*spxout) << " --- using hypersparse pricing, ";
         else
            (*spxout) << " --- using sparse pricing, ";
         (*spxout) << "sparsity: "
                << std::setw(6) << std::fixed << std::setprecision(4)
                << (Real) ninfeasibilities/coDim()
                << std::scientific << std::setprecision(int(prec))
                << std::endl;
      )
      sparsePricingEnterCo = true;
   }
}

Real SPxSolver::computePvec(int i)
{

   return (*thePvec)[i] = vector(i) * (*theCoPvec);
}

Real SPxSolver::computeTest(int i)
{
   SPxBasis::Desc::Status stat = desc().status(i);
   if (isBasic(stat))
      return theTest[i] = 0;
   else
      return theTest[i] = test(i, stat);
}

/*
    Testing all values of #coPvec# against its bounds. If $i$, say, is violated
    the violation is saved as negative value in |theCoTest[i]|.
 */
Real SPxSolver::coTest(int i, SPxBasis::Desc::Status stat) const
{
   assert(type() == ENTER);
   assert(!isBasic(stat));

   Real x;

   switch (stat)
   {
   case SPxBasis::Desc::D_FREE:
   case SPxBasis::Desc::D_ON_BOTH :
      assert(rep() == ROW);
      x = (*theCoPvec)[i] - SPxLP::lower(i);
      if (x < 0)
         return x;
      // no break: next is else case
      //lint -fallthrough
   case SPxBasis::Desc::D_ON_LOWER:
      assert(rep() == ROW);
      return SPxLP::upper(i) - (*theCoPvec)[i];
   case SPxBasis::Desc::D_ON_UPPER:
      assert(rep() == ROW);
      return (*theCoPvec)[i] - SPxLP::lower(i);

   case SPxBasis::Desc::P_ON_UPPER:
      assert(rep() == COLUMN);
      return (*theCoPvec)[i] - maxRowObj(i);             // slacks !
   case SPxBasis::Desc::P_ON_LOWER:
      assert(rep() == COLUMN);
      return maxRowObj(i) - (*theCoPvec)[i];             // slacks !

   default:
      return 0;
   }
}

void SPxSolver::computeCoTest()
{
   int i;
   Real pricingTol = leavetol();
   m_pricingViolUpToDate = true;
   m_pricingViol = 0;
   infeasibilities.clear();
   int ninfeasibilities = 0;
   int sparsitythreshold = (int) (sparsePricingFactor * dim());
   const SPxBasis::Desc& ds = desc();

   for (i = dim() - 1; i >= 0; --i)
   {
      SPxBasis::Desc::Status stat = ds.coStatus(i);
      if (isBasic(stat))
      {
         theCoTest[i] = 0;
         if( remainingRoundsEnter == 0 )
            isInfeasible[i] = SPxPricer::NOT_VIOLATED;
      }
      else
      {
         theCoTest[i] = coTest(i, stat);
         if( remainingRoundsEnter == 0 )
         {
            if( theCoTest[i] < -pricingTol )
            {
               assert(infeasibilities.size() < infeasibilities.max());
               m_pricingViol -= theCoTest[i];
               infeasibilities.addIdx(i);
               isInfeasible[i] = SPxPricer::VIOLATED;
               ++ninfeasibilities;
            }
            else
               isInfeasible[i] = SPxPricer::NOT_VIOLATED;
            if( ninfeasibilities > sparsitythreshold )
            {
               MSG_INFO2( (*spxout), (*spxout) << " --- using dense pricing"
                                 << std::endl; )
               remainingRoundsEnter = DENSEROUNDS;
               sparsePricingEnter = false;
               ninfeasibilities = 0;
            }
         }
         else if( theCoTest[i] < -pricingTol )
            m_pricingViol -= theCoTest[i];
      }
   }
   if( ninfeasibilities == 0 && !sparsePricingEnter )
      --remainingRoundsEnter;
   else if( ninfeasibilities <= sparsitythreshold && !sparsePricingEnter )
   {
      MSG_INFO2( (*spxout),
         std::streamsize prec = spxout->precision();
         if( hyperPricingEnter )
            (*spxout) << " --- using hypersparse pricing, ";
         else
            (*spxout) << " --- using sparse pricing, ";
         (*spxout) << "sparsity: "
                << std::setw(6) << std::fixed << std::setprecision(4)
                << (Real) ninfeasibilities/dim()
                << std::scientific << std::setprecision(int(prec))
                << std::endl;
      )
      sparsePricingEnter = true;
   }
}


/*
    The following methods require propersy initialized vectors |fVec| and
    #coPvec#.
 */
void SPxSolver::updateTest()
{
   thePvec->delta().setup();

   const IdxSet& idx = thePvec->idx();
   const SPxBasis::Desc& ds = desc();
   Real pricingTol = leavetol();

   int i;
   updateViolsCo.clear();
   for (i = idx.size() - 1; i >= 0; --i)
   {
      int j = idx.index(i);
      SPxBasis::Desc::Status stat = ds.status(j);
      if (!isBasic(stat))
      {
         if( m_pricingViolCoUpToDate && theTest[j] < -pricingTol )
            m_pricingViolCo += theTest[j];

         theTest[j] = test(j, stat);

         if( sparsePricingEnterCo )
         {
            if( theTest[j] < -pricingTol )
            {
               assert(remainingRoundsEnterCo == 0);
               m_pricingViolCo -= theTest[j];
               if( isInfeasibleCo[j] == SPxPricer::NOT_VIOLATED )
               {
                  infeasibilitiesCo.addIdx(j);
                  isInfeasibleCo[j] = SPxPricer::VIOLATED;
               }
               if( hyperPricingEnter )
                  updateViolsCo.addIdx(j);
            }
            else
            {
               isInfeasibleCo[j] = SPxPricer::NOT_VIOLATED;
            }
         }
         else if( theTest[j] < -pricingTol )
            m_pricingViolCo -= theTest[j];
      }
      else
      {
         isInfeasibleCo[j] = SPxPricer::NOT_VIOLATED;
         theTest[j] = 0;
      }
   }
}

void SPxSolver::updateCoTest()
{
   theCoPvec->delta().setup();

   const IdxSet& idx = theCoPvec->idx();
   const SPxBasis::Desc& ds = desc();
   Real pricingTol = leavetol();

   int i;
   updateViols.clear();
   for (i = idx.size() - 1; i >= 0; --i)
   {
      int j = idx.index(i);
      SPxBasis::Desc::Status stat = ds.coStatus(j);
      if (!isBasic(stat))
      {
         if( m_pricingViolUpToDate && theCoTest[j] < -pricingTol )
            m_pricingViol += theCoTest[j];

         theCoTest[j] = coTest(j, stat);

         if( sparsePricingEnter )
         {
            if( theCoTest[j] < -pricingTol )
            {
               assert(remainingRoundsEnter == 0);
               m_pricingViol -= theCoTest[j];
               if( isInfeasible[j] == SPxPricer::NOT_VIOLATED )
               {
                  //                if( !hyperPricingEnter )
                  infeasibilities.addIdx(j);
                  isInfeasible[j] = SPxPricer::VIOLATED;
               }
               if( hyperPricingEnter )
                  updateViols.addIdx(j);
            }
            else
            {
               // @todo do we need to remove index j from infeasibilitiesCo?
               isInfeasible[j] = SPxPricer::NOT_VIOLATED;
            }
         }
         else if( theCoTest[j] < -pricingTol )
            m_pricingViol -= theCoTest[j];
      }
      else
      {
         isInfeasible[j] = SPxPricer::NOT_VIOLATED;
         theCoTest[j] = 0;
      }
   }
}



/*  \Section{Compute statistics on entering variable}
    Here is a list of variables relevant when including |Id| to the basis.
    They are computed by |computeEnterStats()|.
 */
void SPxSolver::getEnterVals
(
   SPxId enterId,
   Real& enterTest,
   Real& enterUB,
   Real& enterLB,
   Real& enterVal,
   Real& enterMax,
   Real& enterPric,
   SPxBasis::Desc::Status& enterStat,
   Real& enterRO,
   Real& objChange
)
{
   int enterIdx;
   SPxBasis::Desc& ds = desc();

   if (enterId.isSPxColId())
   {
      enterIdx = number(SPxColId(enterId));
      enterStat = ds.colStatus(enterIdx);
      assert(!isBasic(enterStat));

      /*      For an #Id# to enter the basis we better recompute the Test value.
       */
      if (rep() == COLUMN)
      {
         computePvec(enterIdx);
         enterTest = computeTest(enterIdx);
         theTest[enterIdx] = 0;
      }
      else
      {
         enterTest = coTest()[enterIdx];
         theCoTest[enterIdx] = 0;
      }

      switch (enterStat)
      {
         // primal/columnwise cases:
      case SPxBasis::Desc::P_ON_UPPER :
         assert( rep() == COLUMN );
         enterUB = theUCbound[enterIdx];
         enterLB = theLCbound[enterIdx];
         enterVal = enterUB;
         enterMax = enterLB - enterUB;
         enterPric = (*thePvec)[enterIdx];
         enterRO = maxObj(enterIdx);
         objChange -= enterVal * enterRO;
         if( enterLB <= -infinity )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_LOWER;
         else if( EQ( enterLB, enterUB ) )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert( rep() == COLUMN );
         enterUB = theUCbound[enterIdx];
         enterLB = theLCbound[enterIdx];
         enterVal = enterLB;
         enterMax = enterUB - enterLB;
         enterPric = (*thePvec)[enterIdx];
         enterRO = maxObj(enterIdx);
         objChange -= enterVal * enterRO;
         if( enterUB >= infinity )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_UPPER;
         else if( EQ( enterLB, enterUB ) )
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.colStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_FREE :
         assert( rep() == COLUMN );
         enterUB = theUCbound[enterIdx];
         enterLB = theLCbound[enterIdx];
         enterVal = 0;
         enterPric = (*thePvec)[enterIdx];
         enterRO = maxObj(enterIdx);
         ds.colStatus(enterIdx) = SPxBasis::Desc::D_UNDEFINED;
         enterMax = (enterRO - enterPric > 0) ? infinity : -infinity;
         break;

         // dual/rowwise cases:
      case SPxBasis::Desc::D_ON_UPPER :
         assert( rep() == ROW );
         assert(theUCbound[enterIdx] < infinity);
         enterUB = theUCbound[enterIdx];
         enterLB = -infinity;
         enterMax = -infinity;
         enterVal = enterUB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = SPxLP::lower(enterIdx);
         objChange -= enterRO * enterVal;
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert( rep() == ROW );
         assert(theLCbound[enterIdx] > -infinity);
         enterLB = theLCbound[enterIdx];
         enterUB = infinity;
         enterMax = infinity;
         enterVal = enterLB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = SPxLP::upper(enterIdx);
         objChange -= enterRO * enterVal;
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_FREE:
         assert( rep() == ROW );
         assert(SPxLP::lower(enterIdx) == SPxLP::upper(enterIdx));
         enterUB = infinity;
         enterLB = -infinity;
         enterVal = 0;
         enterRO = SPxLP::upper(enterIdx);
         enterPric = (*theCoPvec)[enterIdx];
         if (enterPric > enterRO)
            enterMax = infinity;
         else
            enterMax = -infinity;
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert( rep() == ROW );
         enterPric = (*theCoPvec)[enterIdx];
         if (enterPric > SPxLP::upper(enterIdx))
         {
            enterLB = theLCbound[enterIdx];
            enterUB = infinity;
            enterMax = infinity;
            enterVal = enterLB;
            enterRO = SPxLP::upper(enterIdx);
            ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            enterUB = theUCbound[enterIdx];
            enterVal = enterUB;
            enterRO = SPxLP::lower(enterIdx);
            enterLB = -infinity;
            enterMax = -infinity;
            ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         }
         objChange -= theLCbound[enterIdx] * SPxLP::upper(enterIdx);
         objChange -= theUCbound[enterIdx] * SPxLP::lower(enterIdx);
         break;
      default:
         throw SPxInternalCodeException("XENTER01 This should never happen.");
      }
      MSG_DEBUG( std::cout << "DENTER03 SPxSolver::getEnterVals() : col " << enterIdx
                        << ": " << enterStat
                        << " -> " << ds.colStatus(enterIdx)
                        << " objChange: " << objChange
                        << std::endl; )
   }

   else
   {
      assert(enterId.isSPxRowId());
      enterIdx = number(SPxRowId(enterId));
      enterStat = ds.rowStatus(enterIdx);
      assert(!isBasic(enterStat));

      /*      For an #Id# to enter the basis we better recompute the Test value.
       */
      if (rep() == ROW)
      {
         computePvec(enterIdx);
         enterTest = computeTest(enterIdx);
         theTest[enterIdx] = 0;
      }
      else
      {
         enterTest = coTest()[enterIdx];
         theCoTest[enterIdx] = 0;
      }

      switch (enterStat)
      {
         // primal/columnwise cases:
      case SPxBasis::Desc::P_ON_UPPER :
         assert( rep() == COLUMN );
         enterUB = theURbound[enterIdx];
         enterLB = theLRbound[enterIdx];
         enterVal = enterLB;
         enterMax = enterUB - enterLB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = maxRowObj(enterIdx);
         objChange -= enterRO * enterVal;
         if( enterUB >= infinity )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_LOWER;
         else if( EQ( enterLB, enterUB ) )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert( rep() == COLUMN );
         enterUB = theURbound[enterIdx];
         enterLB = theLRbound[enterIdx];
         enterVal = enterUB;
         enterMax = enterLB - enterUB;
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = maxRowObj(enterIdx);
         objChange -= enterRO * enterVal;
         if( enterLB <= -infinity )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_UPPER;
         else if( EQ( enterLB, enterUB ) )
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_FREE;
         else
            ds.rowStatus(enterIdx) = SPxBasis::Desc::D_ON_BOTH;
         break;
      case SPxBasis::Desc::P_FREE :
         assert( rep() == COLUMN );
#if 1
         throw SPxInternalCodeException("XENTER02 This should never happen.");
#else
         MSG_ERROR( std::cerr << "EENTER99 ERROR: not yet debugged!" << std::endl; )
         enterPric = (*theCoPvec)[enterIdx];
         enterRO = maxRowObj(enterIdx);
         ds.rowStatus(enterIdx) = SPxBasis::Desc::D_UNDEFINED;
#endif
         break;

         // dual/rowwise cases:
      case SPxBasis::Desc::D_ON_UPPER :
         assert( rep() == ROW );
         assert(theURbound[enterIdx] < infinity);
         enterUB = theURbound[enterIdx];
         enterLB = -infinity;
         enterVal = enterUB;
         enterMax = -infinity;
         enterPric = (*thePvec)[enterIdx];
         enterRO = lhs(enterIdx);
         objChange -= enterRO * enterVal;
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert( rep() == ROW );
         assert(theLRbound[enterIdx] > -infinity);
         enterLB = theLRbound[enterIdx];
         enterUB = infinity;
         enterVal = enterLB;
         enterMax = infinity;
         enterPric = (*thePvec)[enterIdx];
         enterRO = rhs(enterIdx);
         objChange -= enterRO * enterVal;
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_FREE:
         assert( rep() == ROW );
         assert(rhs(enterIdx) == lhs(enterIdx));
         enterUB = infinity;
         enterLB = -infinity;
         enterVal = 0;
         enterPric = (*thePvec)[enterIdx];
         enterRO = rhs(enterIdx);
         enterMax = (enterPric > enterRO) ? infinity : -infinity;
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert( rep() == ROW );
         enterPric = (*thePvec)[enterIdx];
         if (enterPric > rhs(enterIdx))
         {
            enterLB = theLRbound[enterIdx];
            enterVal = enterLB;
            enterUB = infinity;
            enterMax = infinity;
            enterRO = rhs(enterIdx);
            ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            enterUB = theURbound[enterIdx];
            enterVal = enterUB;
            enterLB = -infinity;
            enterMax = -infinity;
            enterRO = lhs(enterIdx);
            ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         }
         objChange -= theLRbound[enterIdx] * rhs(enterIdx);
         objChange -= theURbound[enterIdx] * lhs(enterIdx);
         break;

      default:
         throw SPxInternalCodeException("XENTER03 This should never happen.");
      }
      MSG_DEBUG( std::cout << "DENTER05 SPxSolver::getEnterVals() : row "
                        << enterIdx << ": " << enterStat
                        << " -> " << ds.rowStatus(enterIdx)
                        << " objChange: " << objChange
                        << std::endl; )
   }
}

/*      process leaving variable
 */
void SPxSolver::getEnterVals2
(
   int leaveIdx,
   Real enterMax,
   Real& leavebound,
   Real& objChange
)
{
   int idx;
   SPxBasis::Desc& ds = desc();
   SPxId leftId = baseId(leaveIdx);

   if (leftId.isSPxRowId())
   {
      idx = number(SPxRowId(leftId));
      SPxBasis::Desc::Status leaveStat = ds.rowStatus(idx);

      switch (leaveStat)
      {
      case SPxBasis::Desc::P_FIXED :
         assert(rep() == ROW);
         throw SPxInternalCodeException("XENTER04 This should never happen.");
         break;
      case SPxBasis::Desc::P_ON_UPPER :
         assert(rep() == ROW);
         leavebound = theLBbound[leaveIdx];
         theLRbound[idx] = leavebound;
         ds.rowStatus(idx) = dualRowStatus(idx);
         switch (ds.rowStatus(idx))
         {
         case SPxBasis::Desc::D_ON_UPPER :
            objChange += theURbound[idx] * lhs(idx);
            break;
         case SPxBasis::Desc::D_ON_LOWER :
            objChange += theLRbound[idx] * rhs(idx);
            break;
         case SPxBasis::Desc::D_ON_BOTH :
            objChange += theURbound[idx] * lhs(idx);
            objChange += theLRbound[idx] * rhs(idx);
            break;
         default:
            break;
         }
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert(rep() == ROW);
         leavebound = theUBbound[leaveIdx];
         theURbound[idx] = leavebound;
         ds.rowStatus(idx) = dualRowStatus(idx);
         switch (ds.rowStatus(idx))
         {
         case SPxBasis::Desc::D_ON_UPPER :
            objChange += theURbound[idx] * lhs(idx);
            break;
         case SPxBasis::Desc::D_ON_LOWER :
            objChange += theLRbound[idx] * rhs(idx);
            break;
         case SPxBasis::Desc::D_ON_BOTH :
            objChange += theURbound[idx] * lhs(idx);
            objChange += theLRbound[idx] * rhs(idx);
            break;
         default:
            break;
         }
         break;
      case SPxBasis::Desc::P_FREE :
         assert(rep() == ROW);
#if 1
         throw SPxInternalCodeException("XENTER05 This should never happen.");
#else
         MSG_ERROR( std::cerr << "EENTER98 ERROR: not yet debugged!" << std::endl; )
         if ((*theCoPvec)[leaveIdx] - theLBbound[leaveIdx] <
              theUBbound[leaveIdx] - (*theCoPvec)[leaveIdx])
         {
            leavebound = theLBbound[leaveIdx];
            theLRbound[idx] = leavebound;
         }
         else
         {
            leavebound = theUBbound[leaveIdx];
            theURbound[idx] = leavebound;
         }
         ds.rowStatus(idx) = SPxBasis::Desc::D_UNDEFINED;
#endif
         break;
         // primal/columnwise cases:
      case SPxBasis::Desc::D_UNDEFINED :
         assert(rep() == COLUMN);
         throw SPxInternalCodeException("XENTER06 This should never happen.");
         break;
      case SPxBasis::Desc::D_FREE :
         assert(rep() == COLUMN);
         if (theFvec->delta()[leaveIdx] * enterMax < 0)
            leavebound = theUBbound[leaveIdx];
         else
            leavebound = theLBbound[leaveIdx];
         theLRbound[idx] = leavebound;
         theURbound[idx] = leavebound;
         objChange += leavebound * maxRowObj(leaveIdx);
         ds.rowStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == COLUMN);
         leavebound = theUBbound[leaveIdx];
         theURbound[idx] = leavebound;
         objChange += leavebound * maxRowObj(leaveIdx);
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == COLUMN);
         leavebound = theLBbound[leaveIdx];
         theLRbound[idx] = leavebound;
         objChange += leavebound * maxRowObj(leaveIdx);
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert(rep() == COLUMN);
         if (enterMax * theFvec->delta()[leaveIdx] < 0)
         {
            leavebound = theUBbound[leaveIdx];
            theURbound[idx] = leavebound;
            objChange += leavebound * maxRowObj(leaveIdx);
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         }
         else
         {
            leavebound = theLBbound[leaveIdx];
            theLRbound[idx] = leavebound;
            objChange += leavebound * maxRowObj(leaveIdx);
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         }
         break;

      default:
         throw SPxInternalCodeException("XENTER07 This should never happen.");
      }
      MSG_DEBUG( std::cout << "DENTER06 SPxSolver::getEnterVals2(): row "
                        << idx << ": " << leaveStat
                        << " -> " << ds.rowStatus(idx)
                        << " objChange: " << objChange
                        << std::endl; )
   }

   else
   {
      assert(leftId.isSPxColId());
      idx = number(SPxColId(leftId));
      SPxBasis::Desc::Status leaveStat = ds.colStatus(idx);

      switch (leaveStat)
      {
      case SPxBasis::Desc::P_ON_UPPER :
         assert(rep() == ROW);
         leavebound = theLBbound[leaveIdx];
         theLCbound[idx] = leavebound;
         ds.colStatus(idx) = dualColStatus(idx);
         switch (ds.colStatus(idx))
         {
         case SPxBasis::Desc::D_ON_UPPER :
            objChange += theUCbound[idx] * lower(idx);
            break;
         case SPxBasis::Desc::D_ON_LOWER :
            objChange += theLCbound[idx] * upper(idx);
            break;
         case SPxBasis::Desc::D_ON_BOTH :
            objChange += theLCbound[idx] * upper(idx);
            objChange += theUCbound[idx] * lower(idx);
            break;
         default:
            break;
         }
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert(rep() == ROW);
         leavebound = theUBbound[leaveIdx];
         theUCbound[idx] = leavebound;
         ds.colStatus(idx) = dualColStatus(idx);
         switch (ds.colStatus(idx))
         {
         case SPxBasis::Desc::D_ON_UPPER :
            objChange += theUCbound[idx] * lower(idx);
            break;
         case SPxBasis::Desc::D_ON_LOWER :
            objChange += theLCbound[idx] * upper(idx);
            break;
         case SPxBasis::Desc::D_ON_BOTH :
            objChange += theLCbound[idx] * upper(idx);
            objChange += theUCbound[idx] * lower(idx);
            break;
         default:
            break;
         }
         break;
      case SPxBasis::Desc::P_FREE :
         assert(rep() == ROW);
         if (theFvec->delta()[leaveIdx] * enterMax > 0)
         {
            leavebound = theLBbound[leaveIdx];
            theLCbound[idx] = leavebound;
         }
         else
         {
            leavebound = theUBbound[leaveIdx];
            theUCbound[idx] = leavebound;
         }
         ds.colStatus(idx) = SPxBasis::Desc::D_UNDEFINED;
         break;
      case SPxBasis::Desc::P_FIXED:
         assert(rep() == ROW);
         throw SPxInternalCodeException("XENTER08 This should never happen.");
         break;
         // primal/columnwise cases:
      case SPxBasis::Desc::D_FREE :
         assert(rep() == COLUMN);
         if (theFvec->delta()[leaveIdx] * enterMax > 0)
            leavebound = theLBbound[leaveIdx];
         else
            leavebound = theUBbound[leaveIdx];
         theUCbound[idx] =
            theLCbound[idx] = leavebound;
         objChange += maxObj(idx) * leavebound;
         ds.colStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == COLUMN);
         leavebound = theLBbound[leaveIdx];
         theLCbound[idx] = leavebound;
         objChange += maxObj(idx) * leavebound;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == COLUMN);
         leavebound = theUBbound[leaveIdx];
         theUCbound[idx] = leavebound;
         objChange += maxObj(idx) * leavebound;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
      case SPxBasis::Desc::D_UNDEFINED :
         assert(rep() == COLUMN);
         if (enterMax * theFvec->delta()[leaveIdx] < 0)
         {
            leavebound = theUBbound[leaveIdx];
            theUCbound[idx] = leavebound;
            objChange += maxObj(idx) * leavebound;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            leavebound = theLBbound[leaveIdx];
            theLCbound[idx] = leavebound;
            objChange += maxObj(idx) * leavebound;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         }
         break;

      default:
         throw SPxInternalCodeException("XENTER09 This should never happen.");
      }
      MSG_DEBUG( std::cout << "DENTER07 SPxSolver::getEnterVals2(): col "
                        << idx << ": " << leaveStat
                        << " -> " << ds.colStatus(idx)
                        << " objChange: " << objChange
                        << std::endl; )
   }
}


void
SPxSolver::ungetEnterVal(
   SPxId enterId,
   SPxBasis::Desc::Status enterStat,
   Real leaveVal,
   const SVector& vec,
   Real& objChange
)
{
   assert(rep() == COLUMN);
   int enterIdx;
   SPxBasis::Desc& ds = desc();

   if (enterId.isSPxColId())
   {
      enterIdx = number(SPxColId(enterId));
      if (enterStat == SPxBasis::Desc::P_ON_UPPER)
      {
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         objChange += theLCbound[enterIdx] * maxObj(enterIdx);
      }
      else
      {
         ds.colStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         objChange += theUCbound[enterIdx] * maxObj(enterIdx);
      }
      theFrhs->multAdd(leaveVal, vec);
   }
   else
   {
      enterIdx = number(SPxRowId(enterId));
      assert(enterId.isSPxRowId());
      if (enterStat == SPxBasis::Desc::P_ON_UPPER)
      {
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_LOWER;
         objChange += (theURbound[enterIdx]) * maxRowObj(enterIdx);
      }
      else
      {
         ds.rowStatus(enterIdx) = SPxBasis::Desc::P_ON_UPPER;
         objChange += (theLRbound[enterIdx]) * maxRowObj(enterIdx);
      }
      (*theFrhs)[enterIdx] += leaveVal;
   }
   if (isId(enterId))
   {
      theTest[enterIdx] = 0;
      isInfeasibleCo[enterIdx] = SPxPricer::NOT_VIOLATED;
   }
   else
   {
      theCoTest[enterIdx] = 0;
      isInfeasible[enterIdx] = SPxPricer::NOT_VIOLATED;
   }
}

void SPxSolver::rejectEnter(
   SPxId enterId,
   Real enterTest,
   SPxBasis::Desc::Status enterStat
)
{
   int enterIdx = number(enterId);
   if (isId(enterId))
   {
      theTest[enterIdx] = enterTest;
      desc().status(enterIdx) = enterStat;
   }
   else
   {
      theCoTest[enterIdx] = enterTest;
      desc().coStatus(enterIdx) = enterStat;
   }
}


void SPxSolver::computePrimalray4Col(Real direction, SPxId enterId)
{
   Real sign = (direction > 0 ? 1.0 : -1.0);

   primalRay.clear();
   primalRay.setMax(fVec().delta().size() + 1);

   for( int j = 0; j < fVec().delta().size(); ++j )
   {
      SPxId i = baseId(fVec().idx().index(j));
      if( i.isSPxColId() )
         primalRay.add(number(SPxColId(i)), sign*fVec().delta().value(j));
   }

   if( enterId.isSPxColId() )
      primalRay.add(number(SPxColId(enterId)), -sign);
}


void SPxSolver::computeDualfarkas4Row(Real direction, SPxId enterId)
{
   Real sign = (direction > 0 ? -1.0 : 1.0);

   dualFarkas.clear();
   dualFarkas.setMax(fVec().delta().size() + 1);

   for( int j = 0; j < fVec().delta().size(); ++j )
   {
      SPxId spxid = baseId(fVec().idx().index(j));
      if( spxid.isSPxRowId() )
         dualFarkas.add(number(SPxRowId(spxid)), sign * fVec().delta().value(j));
   }

   if( enterId.isSPxRowId() )
      dualFarkas.add(number(SPxRowId(enterId)), -sign);
}


bool SPxSolver::enter(SPxId& enterId, bool polish)
{
   assert(enterId.isValid());
   assert(type() == ENTER);
   assert(initialized);

   SPxId none;          // invalid id used when enter fails
   Real enterTest;      // correct test value of entering var
   Real enterUB;        // upper bound of entering variable
   Real enterLB;        // lower bound of entering variable
   Real enterVal;       // current value of entering variable
   Real enterMax;       // maximum value for entering shift
   Real enterPric;      // priced value of entering variable
   SPxBasis::Desc::Status enterStat;      // status of entering variable
   Real enterRO;        // rhs/obj of entering variable
   Real objChange = 0.0;
   const SVector* enterVec = enterVector(enterId);

   bool instable = instableEnter;
   assert(!instable || instableEnterId.isValid());

   getEnterVals(enterId, enterTest, enterUB, enterLB,
      enterVal, enterMax, enterPric, enterStat, enterRO, objChange);

   if (!polish && enterTest > -epsilon())
   {
      rejectEnter(enterId, enterTest, enterStat);
      change(-1, none, 0);

      MSG_DEBUG( std::cout << "DENTER08 rejecting false enter pivot" << std::endl; )

      return false;
   }

   /*  Before performing the actual basis update, we must determine, how this
       is to be accomplished.
    */
   // BH 2005-11-15: Obviously solve4update() is only called if theFvec.delta()
   // is setup (i.e. the indices of the NZEs are stored within it) and there are
   // 0 NZEs (???).
   // In that case theFvec->delta() is set such that
   //   Base * theFvec->delta() = enterVec
   if (theFvec->delta().isSetup() && theFvec->delta().size() == 0)
      SPxBasis::solve4update(theFvec->delta(), *enterVec);
#ifdef ENABLE_ADDITIONAL_CHECKS
   else
   {
      // BH 2005-11-29: This code block seems to check the assertion
      //   || Base * theFvec->delta() - enterVec ||_2 <= entertol()
      DVector tmp(dim());
      // BH 2005-11-15: This cast is necessary since SSVector inherits protected from DVector.
      tmp = reinterpret_cast<DVector&>(theFvec->delta());
      multBaseWith(tmp);
      tmp -= *enterVec;
      if (tmp.length() > entertol()) {
         // This happens frequently and does usually not hurt, so print these
         // warnings only with verbose level INFO2 and higher.
         MSG_INFO2( (*spxout), (*spxout) << "WENTER09 fVec updated error = "
                              << tmp.length() << std::endl; )
      }
   }
#endif  // ENABLE_ADDITIONAL_CHECKS

   if (!polish && m_numCycle > m_maxCycle)
   {
      if (-enterMax > 0)
         perturbMaxEnter();
      else
         perturbMinEnter();
   }

   Real leaveVal = -enterMax;

   boundflips = 0;
   int leaveIdx = theratiotester->selectLeave(leaveVal, enterTest, polish);

   /* in row representation, fixed columns and rows should not leave the basis */
   assert(leaveIdx < 0 || !baseId(leaveIdx).isSPxColId() || desc().colStatus(number(SPxColId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);
   assert(leaveIdx < 0 || !baseId(leaveIdx).isSPxRowId() || desc().rowStatus(number(SPxRowId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);

   instableEnterVal = 0;
   instableEnterId = SPxId();
   instableEnter = false;

   /*
       We now tried to find a variable to leave the basis. If one has been
       found, a regular basis update is to be performed.
    */
   if (leaveIdx >= 0)
   {
      if (spxAbs(leaveVal) < entertol())
      {
         if (NE(theUBbound[leaveIdx], theLBbound[leaveIdx])
            && enterStat != Desc::P_FREE && enterStat != Desc::D_FREE)
         {
            m_numCycle++;
            enterCycles++;
         }
      }
      else
         m_numCycle /= 2;

      // setup for updating the copricing vector
      if (coSolveVector3 && coSolveVector2 )
      {
         assert(coSolveVector2->isConsistent());
         assert(coSolveVector2rhs->isSetup());
         assert(coSolveVector3->isConsistent());
         assert(coSolveVector3rhs->isSetup());
         assert(boundflips > 0);
         SPxBasis::coSolve(theCoPvec->delta(), *coSolveVector2, *coSolveVector3
                        , unitVecs[leaveIdx], *coSolveVector2rhs, *coSolveVector3rhs);
         (*theCoPvec) -= (*coSolveVector3);
      }
      else if (coSolveVector3)
      {
         assert(coSolveVector3->isConsistent());
         assert(coSolveVector3rhs->isSetup());
         assert(boundflips > 0);
         SPxBasis::coSolve(theCoPvec->delta(), *coSolveVector3, unitVecs[leaveIdx], *coSolveVector3rhs);
         (*theCoPvec) -= (*coSolveVector3);
      }
      else if (coSolveVector2)
         SPxBasis::coSolve(theCoPvec->delta(), *coSolveVector2, unitVecs[leaveIdx], *coSolveVector2rhs);
      else
         SPxBasis::coSolve(theCoPvec->delta(), unitVecs[leaveIdx]);

      if( boundflips > 0 )
      {
         for( int i = coSolveVector3->dim()-1; i >= 0; --i)
         {
            if( spxAbs((*coSolveVector3)[i]) > epsilon() )
               (*thePvec).multAdd(-(*coSolveVector3)[i],(*thecovectors)[i]);
         }
         // we need to update enterPric in case it was changed by bound flips
         if( enterId.isSPxColId() )
            enterPric = (*theCoPvec)[number(SPxColId(enterId))];
         else
            enterPric = (*thePvec)[number(SPxRowId(enterId))];
         MSG_DEBUG( std::cout << "IEBFRT02 breakpoints passed / bounds flipped = " << boundflips << std::endl; )
         totalboundflips += boundflips;
      }

      (*theCoPrhs)[leaveIdx] = enterRO;
      theCoPvec->value() = (enterRO - enterPric) / theFvec->delta()[leaveIdx];

      if (theCoPvec->value() > epsilon() || theCoPvec->value() < -epsilon())
      {
         if (pricing() == FULL)
         {
            thePvec->value() = theCoPvec->value();
            setupPupdate();
         }
         doPupdate();
      }
      assert(thePvec->isConsistent());
      assert(theCoPvec->isConsistent());

      assert(!baseId(leaveIdx).isSPxRowId() || desc().rowStatus(number(SPxRowId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);
      assert(!baseId(leaveIdx).isSPxColId() || desc().colStatus(number(SPxColId(baseId(leaveIdx)))) != SPxBasis::Desc::P_FIXED);

      Real leavebound;             // bound on which leaving variable moves
      try
      {
         getEnterVals2(leaveIdx, enterMax, leavebound, objChange);
      }
      catch( const SPxException& F )
      {
         rejectEnter(enterId, enterTest, enterStat);
         change(-1, none, 0);
         throw F;
      }

      //  process entering variable
      theUBbound[leaveIdx] = enterUB;
      theLBbound[leaveIdx] = enterLB;

      //  compute tests:
      updateCoTest();
      if (pricing() == FULL)
         updateTest();

      // update feasibility vectors
      theFvec->value() = leaveVal;
      theFvec->update();
      (*theFvec)[leaveIdx] = enterVal - leaveVal;

      if (leavebound > epsilon() || leavebound < -epsilon())
         theFrhs->multAdd(-leavebound, baseVec(leaveIdx));

      if (enterVal > epsilon() || enterVal < -epsilon())
         theFrhs->multAdd(enterVal, *enterVec);

      // update objective funtion value
      updateNonbasicValue(objChange);

      //  change basis matrix
      change(leaveIdx, enterId, enterVec, &(theFvec->delta()));

      return true;
   }
   /*  No leaving vector could be found that would yield a stable pivot step.
    */
   else if (NE(leaveVal, -enterMax))
   {
      /* In the ENTER algorithm, when for a selected entering variable we find only
         an instable leaving variable, then the basis change is not conducted.
         Instead, we save the entering variable's id in instableEnterId and set
         the test value to zero, hoping to find a different leaving
         variable with a stable leavingvariable.
         If this fails, however, and no more entering variable is found, we have to
         perform the instable basis change using instableEnterId. In this (and only
         in this) case, the flag instableEnter is set to true.

         leaveVal != enterMax is the case that selectLeave has found only an instable leaving
         variable. We store this leaving variable for later if we are not already in the
         instable case */

      if (!instable)
      {
         instableEnterId = enterId;
         instableEnterVal = enterTest;

         MSG_DEBUG( std::cout << "DENTER09 rejecting enter pivot and looking for others" << std::endl; )

         rejectEnter(enterId, enterTest / 10.0, enterStat);
         change(-1, none, 0);
      }
      else
      {
         MSG_DEBUG( std::cout << "DENTER10 rejecting enter pivot in instable state, resetting values" << std::endl; )
         rejectEnter(enterId, enterTest, enterStat);
         change(-1, none, 0);
      }

      return false;
   }
   /*  No leaving vector has been selected from the basis. However, if the
       shift amount for |fVec| is bounded, we are in the case, that the
       entering variable is moved from one bound to its other, before any of
       the basis feasibility variables reaches their bound. This may only
       happen in primal/columnwise case with upper and lower bounds on
       variables.
    */
   else if (!polish && leaveVal < infinity && leaveVal > -infinity)
   {
      assert(rep() == COLUMN);
      assert(leaveVal == -enterMax);

      change(-1, enterId, enterVec);

      theFvec->value() = leaveVal;
      theFvec->update();

      ungetEnterVal(enterId, enterStat, leaveVal, *enterVec, objChange);

      // update objective funtion value
      updateNonbasicValue(objChange);

      MSG_DEBUG( std::cout << "DENTER11 moving entering variable from one bound to the other" << std::endl; )

      return false;
   }
   /*  No variable could be selected to leave the basis and even the entering
       variable is unbounded --- this is a failure.
    */
   else
   {
      /* The following line originally was in the "lastUpdate() > 1" case;
         we need it in the INFEASIBLE/UNBOUNDED case, too, to have the
         basis descriptor at the correct size.
       */
      rejectEnter(enterId, enterTest, enterStat);
      change(-1, none, 0);

      if (polish)
         return false;

      else if (lastUpdate() > 1)
      {
         MSG_INFO3( (*spxout), (*spxout) << "IENTER01 factorization triggered in "
                              << "enter() for feasibility test" << std::endl; )
         try
         {
            factorize();
         }
         catch( const SPxStatusException& E )
         {
            // don't exit immediately but handle the singularity correctly
            assert(SPxBasis::status() == SPxBasis::SINGULAR);
            MSG_INFO3( (*spxout), (*spxout) << "Caught exception in factorization: " << E.what() << std::endl; )
         }

         /* after a factorization, the entering column/row might not be infeasible or suboptimal anymore, hence we do
          * not try to call leave(leaveIdx), but rather return to the main solving loop and call the pricer again
          */
         return false;
      }

      /* do not exit with status infeasible or unbounded if there is only a very small violation
       * ROW: recompute the primal variables and activities for another, more precise, round of pricing
       */
      else if (spxAbs(enterTest) < entertol())
      {
         MSG_INFO3( (*spxout), (*spxout) << "IENTER11 clean up step to reduce numerical errors" << std::endl; )

         SPxBasis::coSolve(*theCoPvec, *theCoPrhs);
         computePvec();
         computeCoTest();
         computeTest();

         return false;
      }

      MSG_INFO3( (*spxout), (*spxout) << "IENTER02 unboundedness/infeasibility found in "
                           << "enter()" << std::endl; )

      if (rep() == ROW)
      {
         computeDualfarkas4Row(leaveVal, enterId);
         setBasisStatus(SPxBasis::INFEASIBLE);
      }
      /**@todo if shift() is not zero, we must not conclude primal unboundedness */
      else
      {
         computePrimalray4Col(leaveVal, enterId);
         setBasisStatus(SPxBasis::UNBOUNDED);
      }

      return false;
   }
}
} // namespace soplex
