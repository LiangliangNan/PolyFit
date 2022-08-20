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

/* Updating the Basis for Leaving Variables
 */
#include <assert.h>
#include <stdio.h>

#include "spxdefines.h"
#include "spxpricer.h"
#include "spxsolver.h"
#include "spxratiotester.h"
#include "spxout.h"
#include "exceptions.h"

namespace soplex
{
static const Real reject_leave_tol = 1e-10; // = LOWSTAB as defined in spxfastrt.cpp

/*
    Vector |fTest| gives the feasibility test of all basic variables. For its
    computation |fVec|, |theUBbound| and |theLBbound| must be setup correctly.
    Values of |fTest| $<0$ represent infeasible variables, which are eligible
    for leaving the basis in the simplex loop.
 */
void SPxSolver::computeFtest()
{

   assert(type() == LEAVE);

   Real theeps = entertol();
   m_pricingViolUpToDate = true;
   m_pricingViolCoUpToDate = true;
   m_pricingViol = 0;
   m_pricingViolCo = 0;
   infeasibilities.clear();
   int ninfeasibilities = 0;
   int sparsitythreshold = (int) (sparsePricingFactor * dim());

   for( int i = 0; i < dim(); ++i )
   {
      theCoTest[i] = ((*theFvec)[i] > theUBbound[i])
         ? theUBbound[i] - (*theFvec)[i]
         : (*theFvec)[i] - theLBbound[i];

      if( remainingRoundsLeave == 0 )
      {
         if( theCoTest[i] < -theeps )
         {
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
            remainingRoundsLeave = DENSEROUNDS;
            sparsePricingLeave = false;
            ninfeasibilities = 0;
         }
      }
      else if( theCoTest[i] < -theeps )
            m_pricingViol -= theCoTest[i];
   }

   if( ninfeasibilities == 0 && !sparsePricingLeave )
   {
      --remainingRoundsLeave;
   }
   else if( ninfeasibilities <= sparsitythreshold && !sparsePricingLeave )
   {
      MSG_INFO2( (*spxout),
         std::streamsize prec = spxout->precision();
         if( hyperPricingLeave )
            (*spxout) << " --- using hypersparse pricing, ";
         else
            (*spxout) << " --- using sparse pricing, ";
         (*spxout) << "sparsity: "
                << std::setw(6) << std::fixed << std::setprecision(4)
                << (Real) ninfeasibilities/dim()
                << std::scientific << std::setprecision(int(prec))
                << std::endl;
      )
      sparsePricingLeave = true;
   }
}

void SPxSolver::updateFtest()
{
   const IdxSet& idx = theFvec->idx();
   Vector& ftest = theCoTest;      // |== fTest()|
   assert(&ftest == &fTest());

   assert(type() == LEAVE);

   updateViols.clear();
   Real theeps = entertol();
   for (int j = idx.size() - 1; j >= 0; --j)
   {
      int i = idx.index(j);

      if( m_pricingViolUpToDate && ftest[i] < -theeps )
         m_pricingViol += ftest[i];

      ftest[i] = ((*theFvec)[i] > theUBbound[i])
         ? theUBbound[i] - (*theFvec)[i]
         : (*theFvec)[i] - theLBbound[i];


      if( sparsePricingLeave && ftest[i] < -theeps )
      {
         assert(remainingRoundsLeave == 0);
         if( m_pricingViolUpToDate )
            m_pricingViol -= ftest[i];
         if( isInfeasible[i] == SPxPricer::NOT_VIOLATED )
         {
            // this can cause problems - we cannot keep on adding indeces to infeasibilities,
            // because they are not deleted in hyper mode...
//             if( !hyperPricingLeave )
            infeasibilities.addIdx(i);
            isInfeasible[i] = SPxPricer::VIOLATED;
         }
         if( hyperPricingLeave )
            updateViols.addIdx(i);
      }
      else if( m_pricingViolUpToDate && ftest[i] < -theeps )
         m_pricingViol -= ftest[i];

   }
   // if boundflips were performed, we need to update these indices as well
   if( boundflips > 0 )
   {
      Real eps = epsilon();
      for( int j = 0; j < solveVector3->size(); ++j )
      {
         if( spxAbs(solveVector3->value(j)) > eps )
         {
            int i = solveVector3->index(j);

            if( m_pricingViolUpToDate && ftest[i] < -theeps )
               m_pricingViol += ftest[i];

            ftest[i] = ((*theFvec)[i] > theUBbound[i]) ? theUBbound[i] - (*theFvec)[i] : (*theFvec)[i] - theLBbound[i];

            if( sparsePricingLeave && ftest[i] < -theeps )
            {
               assert(remainingRoundsLeave == 0);
               if( m_pricingViolUpToDate )
                  m_pricingViol -= ftest[i];
               if( !isInfeasible[i] )
               {
                  infeasibilities.addIdx(i);
                  isInfeasible[i] = true;
               }
            }
            else if( m_pricingViolUpToDate && ftest[i] < -theeps )
               m_pricingViol -= ftest[i];
         }
      }
   }
}


/* compute statistics on leaving variable 
   Compute a set of statistical values on the variable selected for leaving the
   basis.
 */
void SPxSolver::getLeaveVals(
   int leaveIdx,
   SPxBasis::Desc::Status& leaveStat,
   SPxId& leaveId,
   Real& leaveMax,
   Real& leavebound,
   int& leaveNum,
   Real& objChange)
{
   SPxBasis::Desc& ds = desc();
   leaveId = baseId(leaveIdx);

   if (leaveId.isSPxRowId())
   {
      leaveNum = number(SPxRowId(leaveId));
      leaveStat = ds.rowStatus(leaveNum);

      assert(isBasic(leaveStat));
      switch (leaveStat)
      {
      case SPxBasis::Desc::P_ON_UPPER :
         assert( rep() == ROW );
         ds.rowStatus(leaveNum) = dualRowStatus(leaveNum);
         leavebound = 0;
         leaveMax = -infinity;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert( rep() == ROW );
         ds.rowStatus(leaveNum) = dualRowStatus(leaveNum);
         leavebound = 0;
         leaveMax = infinity;
         break;
      case SPxBasis::Desc::P_FREE :
         assert( rep() == ROW );
         throw SPxInternalCodeException("XLEAVE01 This should never happen.");
      case SPxBasis::Desc::D_FREE :
         assert( rep() == COLUMN );
         ds.rowStatus(leaveNum) = SPxBasis::Desc::P_FIXED;
         assert(lhs(leaveNum) == rhs(leaveNum));
         leavebound = -rhs(leaveNum);
         if ((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
            leaveMax = infinity;
         else
            leaveMax = -infinity;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert( rep() == COLUMN );
         ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
         leavebound = -rhs(leaveNum);                // slack !!
         leaveMax = infinity;
         objChange += theLRbound[leaveNum] * rhs(leaveNum);
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert( rep() == COLUMN );
         ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
         leavebound = -lhs(leaveNum);                // slack !!
         leaveMax = -infinity;
         objChange += theURbound[leaveNum] * lhs(leaveNum);
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert( rep() == COLUMN );
         if ((*theFvec)[leaveIdx] > theLBbound[leaveIdx])
         {
            ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
            theLRbound[leaveNum] = -infinity;
            leavebound = -lhs(leaveNum);            // slack !!
            leaveMax = -infinity;
            objChange += theURbound[leaveNum] * lhs(leaveNum);
         }
         else
         {
            ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
            theURbound[leaveNum] = infinity;
            leavebound = -rhs(leaveNum);            // slack !!
            leaveMax = infinity;
            objChange += theLRbound[leaveNum] * rhs(leaveNum);
         }
         break;

      default:
         throw SPxInternalCodeException("XLEAVE02 This should never happen.");
      }
      MSG_DEBUG( std::cout << "DLEAVE51 SPxSolver::getLeaveVals() : row " << leaveNum
                        << ": " << leaveStat
                        << " -> " << ds.rowStatus(leaveNum)
                        << " objChange: " << objChange
                        << std::endl; )
   }

   else
   {
      assert(leaveId.isSPxColId());
      leaveNum = number(SPxColId(leaveId));
      leaveStat = ds.colStatus(leaveNum);

      assert(isBasic(leaveStat));
      switch (leaveStat)
      {
      case SPxBasis::Desc::P_ON_UPPER :
         assert( rep() == ROW );
         ds.colStatus(leaveNum) = dualColStatus(leaveNum);
         leavebound = 0;
         leaveMax = -infinity;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert( rep() == ROW );
         ds.colStatus(leaveNum) = dualColStatus(leaveNum);
         leavebound = 0;
         leaveMax = infinity;
         break;
      case SPxBasis::Desc::P_FREE :
         assert( rep() == ROW );
         ds.colStatus(leaveNum) = dualColStatus(leaveNum);
         if ((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
         {
            leavebound = theLBbound[leaveIdx];
            leaveMax = -infinity;
         }
         else
         {
            leavebound = theUBbound[leaveIdx];
            leaveMax = infinity;
         }
         break;

      case SPxBasis::Desc::D_FREE :
         assert( rep() == COLUMN );
         assert(SPxLP::upper(leaveNum) == SPxLP::lower(leaveNum));
         ds.colStatus(leaveNum) = SPxBasis::Desc::P_FIXED;
         leavebound = SPxLP::upper(leaveNum);
         objChange += maxObj(leaveNum) * leavebound;
         if ((*theFvec)[leaveIdx] < theLBbound[leaveIdx])
            leaveMax = infinity;
         else
            leaveMax = -infinity;
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert( rep() == COLUMN );
         ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
         leavebound = SPxLP::upper(leaveNum);
         objChange += theUCbound[leaveNum] * leavebound;
         leaveMax = -infinity;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert( rep() == COLUMN );
         ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
         leavebound = SPxLP::lower(leaveNum);
         objChange += theLCbound[leaveNum] * leavebound;
         leaveMax = infinity;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert( rep() == COLUMN );
         if ((*theFvec)[leaveIdx] > theUBbound[leaveIdx])
         {
            leaveMax = -infinity;
            leavebound = SPxLP::upper(leaveNum);
            objChange += theUCbound[leaveNum] * leavebound;
            theLCbound[leaveNum] = -infinity;
            ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
         }
         else
         {
            leaveMax = infinity;
            leavebound = SPxLP::lower(leaveNum);
            objChange += theLCbound[leaveNum] * leavebound;
            theUCbound[leaveNum] = infinity;
            ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
         }
         break;
      default:
         throw SPxInternalCodeException("XLEAVE03 This should never happen.");
      }
      MSG_DEBUG( std::cout << "DLEAVE52 SPxSolver::getLeaveVals() : col " << leaveNum
                        << ": " << leaveStat
                        << " -> " << ds.colStatus(leaveNum)
                        << " objChange: " << objChange
                        << std::endl; )
   }
}

void SPxSolver::getLeaveVals2(
   Real leaveMax,
   SPxId enterId,
   Real& enterBound,
   Real& newUBbound,
   Real& newLBbound,
   Real& newCoPrhs,
   Real& objChange
)
{
   SPxBasis::Desc& ds = desc();

   enterBound = 0;
   if (enterId.isSPxRowId())
   {
      int idx = number(SPxRowId(enterId));
      SPxBasis::Desc::Status enterStat = ds.rowStatus(idx);

      switch (enterStat)
      {
      case SPxBasis::Desc::D_FREE :
         assert(rep() == ROW);
         if (thePvec->delta()[idx] * leaveMax < 0)
            newCoPrhs = theLRbound[idx];
         else
            newCoPrhs = theURbound[idx];
         newUBbound = infinity;
         newLBbound = -infinity;
         ds.rowStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == ROW);
         newUBbound = 0;
         newLBbound = -infinity;
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         newCoPrhs = theLRbound[idx];
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == ROW);
         newUBbound = infinity;
         newLBbound = 0;
         ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         newCoPrhs = theURbound[idx];
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert(rep() == ROW);
         if (leaveMax * thePvec->delta()[idx] < 0)
         {
            newUBbound = 0;
            newLBbound = -infinity;
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
            newCoPrhs = theLRbound[idx];
         }
         else
         {
            newUBbound = infinity;
            newLBbound = 0;
            ds.rowStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
            newCoPrhs = theURbound[idx];
         }
         break;

      case SPxBasis::Desc::P_ON_UPPER :
         assert(rep() == COLUMN);
         ds.rowStatus(idx) = dualRowStatus(idx);
         if (lhs(idx) > -infinity)
            theURbound[idx] = theLRbound[idx];
         newCoPrhs = theLRbound[idx];        // slack !!
         newUBbound = -lhs(idx);
         newLBbound = -rhs(idx);
         enterBound = -rhs(idx);
         objChange -= newCoPrhs * rhs(idx);
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert(rep() == COLUMN);
         ds.rowStatus(idx) = dualRowStatus(idx);
         if (rhs(idx) < infinity)
            theLRbound[idx] = theURbound[idx];
         newCoPrhs = theURbound[idx];        // slack !!
         newLBbound = -rhs(idx);
         newUBbound = -lhs(idx);
         enterBound = -lhs(idx);
         objChange -= newCoPrhs * lhs(idx);
         break;
      case SPxBasis::Desc::P_FREE :
         assert(rep() == COLUMN);
#if 1
         throw SPxInternalCodeException("XLEAVE04 This should never happen.");
#else
         MSG_ERROR( std::cerr << "ELEAVE53 ERROR: not yet debugged!" << std::endl; )
         ds.rowStatus(idx) = dualRowStatus(idx);
         newCoPrhs = theURbound[idx];        // slack !!
         newUBbound = infinity;
         newLBbound = -infinity;
         enterBound = 0;
#endif
         break;
      case SPxBasis::Desc::P_FIXED :
         assert(rep() == COLUMN);
         MSG_ERROR( std::cerr << "ELEAVE54 "
                           << "ERROR! Tried to put a fixed row variable into the basis: "
                           << "idx="   << idx
                           << ", lhs=" << lhs(idx)
                           << ", rhs=" << rhs(idx) << std::endl; )
         throw SPxInternalCodeException("XLEAVE05 This should never happen.");

      default:
         throw SPxInternalCodeException("XLEAVE06 This should never happen.");
      }
      MSG_DEBUG( std::cout << "DLEAVE55 SPxSolver::getLeaveVals2(): row " << idx
                        << ": " << enterStat
                        << " -> " << ds.rowStatus(idx)
                        << " objChange: " << objChange
                        << std::endl; )
   }

   else
   {
      assert(enterId.isSPxColId());
      int idx = number(SPxColId(enterId));
      SPxBasis::Desc::Status enterStat = ds.colStatus(idx);

      switch (enterStat)
      {
      case SPxBasis::Desc::D_ON_UPPER :
         assert(rep() == ROW);
         newUBbound = 0;
         newLBbound = -infinity;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
         newCoPrhs = theLCbound[idx];
         break;
      case SPxBasis::Desc::D_ON_LOWER :
         assert(rep() == ROW);
         newUBbound = infinity;
         newLBbound = 0;
         ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
         newCoPrhs = theUCbound[idx];
         break;
      case SPxBasis::Desc::D_FREE :
         assert(rep() == ROW);
         newUBbound = infinity;
         newLBbound = -infinity;
         newCoPrhs = theLCbound[idx];
         ds.colStatus(idx) = SPxBasis::Desc::P_FIXED;
         break;
      case SPxBasis::Desc::D_ON_BOTH :
         assert(rep() == ROW);
         if (leaveMax * theCoPvec->delta()[idx] < 0)
         {
            newUBbound = 0;
            newLBbound = -infinity;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_LOWER;
            newCoPrhs = theLCbound[idx];
         }
         else
         {
            newUBbound = infinity;
            newLBbound = 0;
            ds.colStatus(idx) = SPxBasis::Desc::P_ON_UPPER;
            newCoPrhs = theUCbound[idx];
         }
         break;

      case SPxBasis::Desc::P_ON_UPPER :
         assert(rep() == COLUMN);
         ds.colStatus(idx) = dualColStatus(idx);
         if (SPxLP::lower(idx) > -infinity)
            theLCbound[idx] = theUCbound[idx];
         newCoPrhs = theUCbound[idx];
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = SPxLP::upper(idx);
         objChange -= newCoPrhs * enterBound;
         break;
      case SPxBasis::Desc::P_ON_LOWER :
         assert(rep() == COLUMN);
         ds.colStatus(idx) = dualColStatus(idx);
         if (SPxLP::upper(idx) < infinity)
            theUCbound[idx] = theLCbound[idx];
         newCoPrhs = theLCbound[idx];
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = SPxLP::lower(idx);
         objChange -= newCoPrhs * enterBound;
         break;
      case SPxBasis::Desc::P_FREE :
         assert(rep() == COLUMN);
         ds.colStatus(idx) = dualColStatus(idx);
         if (thePvec->delta()[idx] * leaveMax > 0)
            newCoPrhs = theUCbound[idx];
         else
            newCoPrhs = theLCbound[idx];
         newUBbound = SPxLP::upper(idx);
         newLBbound = SPxLP::lower(idx);
         enterBound = 0;
         break;
      case SPxBasis::Desc::P_FIXED :
         assert(rep() == COLUMN);
         MSG_ERROR( std::cerr << "ELEAVE56 "
                           << "ERROR! Tried to put a fixed column variable into the basis. "
                           << "idx="     << idx
                           << ", lower=" << lower(idx)
                           << ", upper=" << upper(idx) << std::endl; )
         throw SPxInternalCodeException("XLEAVE07 This should never happen.");
      default:
         throw SPxInternalCodeException("XLEAVE08 This should never happen.");
      }

      MSG_DEBUG( std::cout << "DLEAVE57 SPxSolver::getLeaveVals2(): col " << idx
                        << ": " << enterStat
                        << " -> " << ds.colStatus(idx)
                        << " objChange: " << objChange
                        << std::endl; )
   }

}

void SPxSolver::rejectLeave(
   int leaveNum,
   SPxId leaveId,
   SPxBasis::Desc::Status leaveStat,
   const SVector* //newVec
)
{
   SPxBasis::Desc& ds = desc();
   if (leaveId.isSPxRowId())
   {
      MSG_DEBUG( std::cout << "DLEAVE58 rejectLeave()  : row " << leaveNum
                        << ": " << ds.rowStatus(leaveNum)
                        << " -> " << leaveStat << std::endl; )

      if (leaveStat == SPxBasis::Desc::D_ON_BOTH)
      {
         if (ds.rowStatus(leaveNum) == SPxBasis::Desc::P_ON_LOWER)
            theLRbound[leaveNum] = theURbound[leaveNum];
         else
            theURbound[leaveNum] = theLRbound[leaveNum];
      }
      ds.rowStatus(leaveNum) = leaveStat;
   }
   else
   {
      MSG_DEBUG( std::cout << "DLEAVE59 rejectLeave()  : col " << leaveNum
                        << ": " << ds.colStatus(leaveNum)
                        << " -> " << leaveStat << std::endl; )

      if (leaveStat == SPxBasis::Desc::D_ON_BOTH)
      {
         if (ds.colStatus(leaveNum) == SPxBasis::Desc::P_ON_UPPER)
            theLCbound[leaveNum] = theUCbound[leaveNum];
         else
            theUCbound[leaveNum] = theLCbound[leaveNum];
      }
      ds.colStatus(leaveNum) = leaveStat;
   }
}


void SPxSolver::computePrimalray4Row(Real direction)
{
   Real sign = (direction > 0 ? 1.0 : -1.0);

   primalRay.clear();
   primalRay.setMax(coPvec().delta().size());

   for( int i = 0; i < coPvec().delta().size(); ++i )
      primalRay.add(coPvec().delta().index(i), sign * coPvec().delta().value(i));
}

void SPxSolver::computeDualfarkas4Col(Real direction)
{
   Real sign = (direction > 0 ? -1.0 : 1.0);

   dualFarkas.clear();
   dualFarkas.setMax(coPvec().delta().size());

   for( int i = 0; i < coPvec().delta().size(); ++i )
      dualFarkas.add(coPvec().delta().index(i), sign * coPvec().delta().value(i));
}

bool SPxSolver::leave(int leaveIdx, bool polish)
{
   assert(leaveIdx < dim() && leaveIdx >= 0);
   assert(type() == LEAVE);
   assert(initialized);

   bool instable = instableLeave;
   assert(!instable || instableLeaveNum >= 0);

   /*
       Before performing the actual basis update, we must determine, how this
       is to be accomplished.
       When using steepest edge pricing this solve is already performed by the pricer
    */
   if (theCoPvec->delta().isSetup() && theCoPvec->delta().size() == 0)
   {
      coSolve(theCoPvec->delta(), unitVecs[leaveIdx]);
   }
#ifdef ENABLE_ADDITIONAL_CHECKS
   else
   {
      SSVector tmp(dim(), epsilon());
      tmp.clear();
      coSolve(tmp, unitVecs[leaveIdx]);
      tmp -= theCoPvec->delta();
      if (tmp.length() > leavetol()) {
         // This happens very frequently and does usually not hurt, so print
         // these warnings only with verbose level INFO2 and higher.
         MSG_INFO2( (*spxout), (*spxout) << "WLEAVE60 iteration=" << basis().iteration()
                              << ": coPvec.delta error = " << tmp.length() 
                              << std::endl; )
      }
   }
#endif  // ENABLE_ADDITIONAL_CHECKS

   setupPupdate();

   assert(thePvec->isConsistent());
   assert(theCoPvec->isConsistent());

   SPxBasis::Desc::Status leaveStat;      // status of leaving var
   SPxId leaveId;        // id of leaving var
   SPxId none;           // invalid id used if leave fails
   Real leaveMax;       // maximium lambda of leaving var
   Real leavebound;     // current fVec value of leaving var
   int  leaveNum;       // number of leaveId in bounds
   Real objChange = 0.0; // amount of change in the objective function

   getLeaveVals(leaveIdx, leaveStat, leaveId, leaveMax, leavebound, leaveNum, objChange);

   if (!polish && m_numCycle > m_maxCycle)
   {
      if (leaveMax > 0)
         perturbMaxLeave();
      else
         perturbMinLeave();
      //@ m_numCycle /= 2;
      // perturbation invalidates the currently stored nonbasic value
      forceRecompNonbasicValue();
   }
   //@ testBounds();

   Real enterVal = leaveMax;
   boundflips = 0;
   Real oldShift = theShift;
   SPxId enterId = theratiotester->selectEnter(enterVal, leaveIdx, polish);
   if (NE(theShift, oldShift))
   {
      MSG_DEBUG( std::cout << "DLEAVE71 trigger recomputation of nonbasic value due to shifts in ratiotest" << std::endl; )
      forceRecompNonbasicValue();
   }

   assert(!enterId.isValid() || !isBasic(enterId));

   instableLeaveNum = -1;
   instableLeave = false;

   /*
       No variable could be selected to enter the basis and even the leaving
       variable is unbounded.
    */
   if (!enterId.isValid())
   {
      /* the following line originally was below in "rejecting leave" case;
         we need it in the unbounded/infeasible case, too, to have the
         correct basis size */
      rejectLeave(leaveNum, leaveId, leaveStat);
      change(-1, none, 0);
      objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case

      if (polish)
         return false;

      if (NE(enterVal, leaveMax))
      {
         MSG_DEBUG( std::cout << "DLEAVE61 rejecting leave A (leaveIdx=" << leaveIdx
                           << ", theCoTest=" << theCoTest[leaveIdx] << ")"
                           << std::endl; )

         /* In the LEAVE algorithm, when for a selected leaving variable we find only
            an instable entering variable, then the basis change is not conducted.
            Instead, we save the leaving variable's index in instableLeaveNum and scale
            theCoTest[leaveIdx] down by some factor, hoping to find a different leaving
            variable with a stable entering variable.
            If this fails, however, and no more leaving variable is found, we have to
            perform the instable basis change using instableLeaveNum. In this (and only
            in this) case, the flag instableLeave is set to true.

            enterVal != leaveMax is the case that selectEnter has found only an instable entering
            variable. We store this leaving variable for later -- if we are not already in the
            instable case: then we continue and conclude unboundedness/infeasibility */
         if (!instable)
         {
            instableLeaveNum = leaveIdx;

            // Note: These changes do not survive a refactorization
            instableLeaveVal = theCoTest[leaveIdx];
            theCoTest[leaveIdx] = instableLeaveVal / 10.0;

            return true;
         }
      }

      if (lastUpdate() > 1)
      {
         MSG_INFO3( (*spxout), (*spxout) << "ILEAVE01 factorization triggered in "
                              << "leave() for feasibility test" << std::endl; )
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

         /* after a factorization, the leaving column/row might not be infeasible or suboptimal anymore, hence we do
          * not try to call leave(leaveIdx), but rather return to the main solving loop and call the pricer again
          */
         return true;
      }

      /* do not exit with status infeasible or unbounded if there is only a very small violation */
      if (spxAbs(enterVal) < leavetol())
      {
         MSG_INFO3( (*spxout), (*spxout) << "ILEAVE11 clean up step to reduce numerical errors" << std::endl; )

         computeFrhs();
         SPxBasis::solve(*theFvec, *theFrhs);
         computeFtest();

         return true;
      }
      MSG_INFO3( (*spxout), (*spxout) << "ILEAVE02 unboundedness/infeasibility found "
                           << "in leave()" << std::endl; )

      if (rep() != COLUMN)
      {
         computePrimalray4Row(enterVal);
         setBasisStatus(SPxBasis::UNBOUNDED);
      }
      else
      {
         computeDualfarkas4Col(enterVal);
         setBasisStatus(SPxBasis::INFEASIBLE);
      }
      return false;
   }
   else
   {
      /*
        If an entering variable has been found, a regular basis update is to
        be performed.
      */
      if (enterId != baseId(leaveIdx))
      {
         const SVector& newVector = *enterVector(enterId);
         // update feasibility vectors
         if( solveVector2 != NULL && solveVector3 != NULL )
         {
            assert(solveVector2->isConsistent());
            assert(solveVector2rhs->isSetup());
            assert(solveVector3->isConsistent());
            assert(solveVector3rhs->isSetup());
            assert(boundflips > 0);
            SPxBasis::solve4update(theFvec->delta(),
                                   *solveVector2,
                                   *solveVector3,
                                   newVector,
                                   *solveVector2rhs,
                                   *solveVector3rhs);

            // perform update of basic solution
            primVec -= (*solveVector3);
            MSG_DEBUG( std::cout << "ILBFRT02 breakpoints passed / bounds flipped = " << boundflips << std::endl; )
            totalboundflips += boundflips;
         }
         else if( solveVector2 != NULL )
         {
            assert(solveVector2->isConsistent());
            assert(solveVector2rhs->isSetup());

            SPxBasis::solve4update(theFvec->delta(),
                                   *solveVector2,
                                   newVector,
                                   *solveVector2rhs);
         }
         else if( solveVector3 != NULL )
         {
            assert(solveVector3->isConsistent());
            assert(solveVector3rhs->isSetup());
            assert(boundflips > 0);
            SPxBasis::solve4update(theFvec->delta(),
                                   *solveVector3,
                                   newVector,
                                   *solveVector3rhs);

            // perform update of basic solution
            primVec -= (*solveVector3);
            MSG_DEBUG( std::cout << "ILBFRT02 breakpoints passed / bounds flipped = " << boundflips << std::endl; )
            totalboundflips += boundflips;
         }
         else
            SPxBasis::solve4update (theFvec->delta(), newVector);

#ifdef ENABLE_ADDITIONAL_CHECKS
         {
            SSVector tmp(dim(), epsilon());
            SPxBasis::solve(tmp, newVector);
            tmp -= fVec().delta();
            if (tmp.length() > entertol()) {
               // This happens very frequently and does usually not hurt, so print
               // these warnings only with verbose level INFO2 and higher.
               MSG_INFO2( (*spxout), (*spxout) << "WLEAVE62\t(" << tmp.length() << ")\n"; )
                  }
         }
#endif  // ENABLE_ADDITIONAL_CHECKS


         if (spxAbs(theFvec->delta()[leaveIdx]) < reject_leave_tol)
         {
            if (instable)
            {
               /* We are in the case that for all leaving variables only instable entering
                  variables were found: Thus, above we already accepted such an instable
                  entering variable. Now even this seems to be impossible, thus we conclude
                  unboundedness/infeasibility. */
               MSG_INFO3( (*spxout), (*spxout) << "ILEAVE03 unboundedness/infeasibility found "
                  << "in leave()" << std::endl; )

               rejectLeave(leaveNum, leaveId, leaveStat);
               change(-1, none, 0);
               objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case

               /**@todo if shift() is not zero we must not conclude unboundedness */
               if (rep() == ROW)
               {
                  computePrimalray4Row(enterVal);
                  setBasisStatus(SPxBasis::UNBOUNDED);
               }
               else
               {
                  computeDualfarkas4Col(enterVal);
                  setBasisStatus(SPxBasis::INFEASIBLE);
               }

               return false;
            }
            else
            {
               theFvec->delta().clear();
               rejectLeave(leaveNum, leaveId, leaveStat, &newVector);
               change(-1, none, 0);
               objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case

               MSG_DEBUG( std::cout << "DLEAVE63 rejecting leave B (leaveIdx=" << leaveIdx
                  << ", theCoTest=" << theCoTest[leaveIdx]
                  << ")" << std::endl; )

               // Note: These changes do not survive a refactorization
               theCoTest[leaveIdx] *= 0.01;

               return true;
            }
         }

         //      process leaving variable
         if (leavebound > epsilon() || leavebound < -epsilon())
            theFrhs->multAdd(-leavebound, baseVec(leaveIdx));

         //      process entering variable
         Real enterBound;
         Real newUBbound;
         Real newLBbound;
         Real newCoPrhs;

         try
         {
            getLeaveVals2(leaveMax, enterId, enterBound, newUBbound, newLBbound, newCoPrhs, objChange);
         }
         catch( const SPxException& F )
         {
            rejectLeave(leaveNum, leaveId, leaveStat);
            change(-1, none, 0);
            objChange = 0.0; // the nonbasicValue is not supposed to be updated in this case
            throw F;
         }

         theUBbound[leaveIdx] = newUBbound;
         theLBbound[leaveIdx] = newLBbound;
         (*theCoPrhs)[leaveIdx] = newCoPrhs;

         if (enterBound > epsilon() || enterBound < -epsilon())
            theFrhs->multAdd(enterBound, newVector);

         // update pricing vectors
         theCoPvec->value() = enterVal;
         thePvec->value() = enterVal;
         if (enterVal > epsilon() || enterVal < -epsilon())
            doPupdate();

         // update feasibility vector
         theFvec->value() = -((*theFvec)[leaveIdx] - leavebound)
            / theFvec->delta()[leaveIdx];
         theFvec->update();
         (*theFvec)[leaveIdx] = enterBound - theFvec->value();
         updateFtest();

         // update objective funtion value
         updateNonbasicValue(objChange);

         //  change basis matrix
         change(leaveIdx, enterId, &newVector, &(theFvec->delta()));
      }

      /*
        No entering vector has been selected from the basis. However, if the
        shift amount for |coPvec| is bounded, we are in the case, that the
        entering variable is moved from one bound to its other, before any of
        the basis feasibility variables reaches their bound. This may only
        happen in primal/columnwise case with upper and lower bounds on
        variables.
      */
      else
      {
         // @todo update obj function value here!!!
         assert(rep() == ROW);
         SPxBasis::Desc& ds = desc();

         change(leaveIdx, none, 0);

         if (leaveStat == SPxBasis::Desc::P_ON_UPPER)
         {
            if (leaveId.isSPxRowId())
            {
               ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
               (*theCoPrhs)[leaveIdx] = theLRbound[leaveNum];
            }
            else
            {
               ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_LOWER;
               (*theCoPrhs)[leaveIdx] = theLCbound[leaveNum];
            }
            theUBbound[leaveIdx] = 0;
            theLBbound[leaveIdx] = -infinity;
         }
         else
         {
            assert( leaveStat == SPxBasis::Desc::P_ON_LOWER );
            if (leaveId.isSPxRowId())
            {
               ds.rowStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
               (*theCoPrhs)[leaveIdx] = theURbound[leaveNum];
            }
            else
            {
               ds.colStatus(leaveNum) = SPxBasis::Desc::P_ON_UPPER;
               (*theCoPrhs)[leaveIdx] = theUCbound[leaveNum];
            }
            theUBbound[leaveIdx] = infinity;
            theLBbound[leaveIdx] = 0;
         }

         // update copricing vector
         theCoPvec->value() = enterVal;
         thePvec->value() = enterVal;
         if (enterVal > epsilon() || enterVal < -epsilon())
            doPupdate();

         // update feasibility vectors
         theFvec->value() = 0;
         assert(theCoTest[leaveIdx] < 0.0);
         m_pricingViol += theCoTest[leaveIdx];
         theCoTest[leaveIdx] *= -1;
      }

      if ((leaveMax > entertol() && enterVal <= entertol()) || (leaveMax < -entertol() && enterVal >= -entertol()))
      {
         if ((theUBbound[leaveIdx] < infinity || theLBbound[leaveIdx] > -infinity)
            && leaveStat != SPxBasis::Desc::P_FREE
            && leaveStat != SPxBasis::Desc::D_FREE)
            {
            m_numCycle++;
               leaveCycles++;
            }
      }
      else
         m_numCycle /= 2;

#ifdef ENABLE_ADDITIONAL_CHECKS
      {
         DVector tmp = fVec();
         multBaseWith(tmp);
         tmp -= fRhs();
         if (tmp.length() > entertol())
         {
            // This happens very frequently and does usually not hurt, so print
            // these warnings only with verbose level INFO2 and higher.
            MSG_INFO2( (*spxout), (*spxout) << "WLEAVE64\t" << basis().iteration()
               << ": fVec error = " << tmp.length() << std::endl; )
               SPxBasis::solve(tmp, fRhs());
            tmp -= fVec();
            MSG_INFO2( (*spxout), (*spxout) << "WLEAVE65\t(" << tmp.length() << ")\n"; )
               }
      }
#endif  // ENABLE_ADDITIONAL_CHECKS

      return true;
   }
}
} // namespace soplex
