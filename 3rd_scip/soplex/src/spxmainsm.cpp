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

#include "spxmainsm.h"
#include "array.h"
#include "dataarray.h"
#include "sorter.h"
#include "spxout.h"
#include <sstream>
#include <iostream>
#include <fstream>


//rows
#define FREE_LHS_RHS            1
#define FREE_CONSTRAINT         1
#define EMPTY_CONSTRAINT        1
#define ROW_SINGLETON           1
#define FORCE_CONSTRAINT        1
//cols
#define FREE_BOUNDS             1
#define EMPTY_COLUMN            1
#define FIX_VARIABLE            1
#define FREE_ZERO_OBJ_VARIABLE  1
#define ZERO_OBJ_COL_SINGLETON  1
#define DOUBLETON_EQUATION      1
#define FREE_COL_SINGLETON      1
//dual
#define DOMINATED_COLUMN        1
#define WEAKLY_DOMINATED_COLUMN 1
#define MULTI_AGGREGATE         1
//other
#define TRIVIAL_HEURISTICS      1
#define PSEUDOOBJ               1


#define EXTREMES                1
#define ROWS                    1
#define COLS                    1
#define DUAL                    1
///@todo check: with this simplification step, the unsimplified basis seems to be slightly suboptimal for some instances
#define DUPLICATE_ROWS          1
#define DUPLICATE_COLS          1


#ifndef NDEBUG
#define CHECK_BASIC_DIM
#endif  // NDEBUG

namespace soplex
{
bool SPxMainSM::PostStep::checkBasisDim(DataArray<SPxSolver::VarStatus> rows,
                                        DataArray<SPxSolver::VarStatus> cols) const
{
   int numBasis = 0;
   for(int rs = 0; rs < nRows; ++rs)
   {
      if(rows[rs] == SPxSolver::BASIC)
         numBasis++;
   }

   for(int cs = 0; cs < nCols; ++cs)
   {
      if(cols[cs] == SPxSolver::BASIC)
         numBasis++;
   }
   assert(numBasis == nRows);
   return numBasis == nRows;
}

void SPxMainSM::RowObjPS::execute(DVector& x, DVector& y, DVector& s, DVector&,
                                          DataArray<SPxSolver::VarStatus>& cStatus,
                                          DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   assert(isZero(s[m_i], 1e-9));
   s[m_i] = -x[m_j];

   assert(rStatus[m_i] != SPxSolver::UNDEFINED);
   assert(cStatus[m_j] != SPxSolver::UNDEFINED);
   assert(rStatus[m_i] != SPxSolver::BASIC || cStatus[m_j] != SPxSolver::BASIC);

   MSG_DEBUG( std::cout << "RowObjPS: removing slack column " << m_j << " (" << cStatus[m_j] << ") for row " << m_i << " (" << rStatus[m_i] << ").\n" );

   if( rStatus[m_i] != SPxSolver::BASIC )
   {
      switch( cStatus[m_j] )
      {
      case SPxSolver::ON_UPPER:
         rStatus[m_i] = SPxSolver::ON_LOWER;
         break;
      case SPxSolver::ON_LOWER:
         rStatus[m_i] = SPxSolver::ON_UPPER;
         break;
      default:
         rStatus[m_i] = cStatus[m_j];
      }

      // otherwise checkBasisDim() may fail
      cStatus[m_j] = SPxSolver::ZERO;
   }

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      assert(false);
      throw SPxInternalCodeException("XMAISM15 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::FreeConstraintPS::execute(DVector& x, DVector& y, DVector& s, DVector&,
                                          DataArray<SPxSolver::VarStatus>& cStatus,
                                          DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // correcting the change of idx by deletion of the row:
   s[m_old_i] = s[m_i];
   y[m_old_i] = y[m_i];
   rStatus[m_old_i] = rStatus[m_i];

   // primal:
   Real slack = 0.0;

   for (int k = 0; k < m_row.size(); ++k)
      slack += m_row.value(k) * x[m_row.index(k)];

   s[m_i] = slack;

   // dual:
   y[m_i] = m_row_obj;

   // basis:
   rStatus[m_i] = SPxSolver::BASIC;

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
       throw SPxInternalCodeException("XMAISM15 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::EmptyConstraintPS::execute(DVector&, DVector& y, DVector& s, DVector&,
                                           DataArray<SPxSolver::VarStatus>& cStatus,
                                           DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // correcting the change of idx by deletion of the row:
   s[m_old_i] = s[m_i];
   y[m_old_i] = y[m_i];
   rStatus[m_old_i] = rStatus[m_i];

   // primal:
   s[m_i] = 0.0;

   // dual:
   y[m_i] = m_row_obj;

   // basis:
   rStatus[m_i] = SPxSolver::BASIC;

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
       throw SPxInternalCodeException("XMAISM16 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::RowSingletonPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                        DataArray<SPxSolver::VarStatus>& cStatus,
                                        DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // correcting the change of idx by deletion of the row:
   s[m_old_i] = s[m_i];
   y[m_old_i] = y[m_i];
   rStatus[m_old_i] = rStatus[m_i];

   Real aij = m_col[m_i];

   // primal:
   s[m_i] = aij * x[m_j];

   // dual & basis:
   Real val = m_obj;

   for(int k = 0; k < m_col.size(); ++k)
   {
      if (m_col.index(k) != m_i)
         val -= m_col.value(k) * y[m_col.index(k)];
   }

   Real newLo = (aij > 0) ? m_lhs/aij : m_rhs/aij;  // implicit lhs
   Real newUp = (aij > 0) ? m_rhs/aij : m_lhs/aij;  // implicit rhs

   switch(cStatus[m_j])
   {
   case SPxSolver::FIXED:
      if(newLo <= m_oldLo && newUp >= m_oldUp)
      {
         // this row is totally redundant, has not changed bound of xj
         rStatus[m_i] = SPxSolver::BASIC;
         y[m_i] = m_row_obj;
      }
      else if(EQrel(newLo, newUp, eps()))
      {
         // row is in the type  aij * xj = b
         assert(EQrel(newLo, x[m_j], eps()));

         if(EQrel(m_oldLo, m_oldUp, eps()))
         {
            // xj has been fixed in other row
            rStatus[m_i] = SPxSolver::BASIC;
            y[m_i] = m_row_obj;
         }
         else if((EQrel(m_oldLo, x[m_j], eps()) && r[m_j] <= -eps())
                 || (EQrel(m_oldUp, x[m_j], eps()) && r[m_j] >= eps())
                 || (!EQrel(m_oldLo, x[m_j], eps()) && !(EQrel(m_oldUp, x[m_j], eps()))))
         {
            // if x_j on lower but reduced cost is negative, or x_j on upper but reduced cost is positive, or x_j not on bound: basic
            rStatus[m_i] = (EQrel(m_lhs, x[m_j]*aij, eps())) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
            cStatus[m_j] = SPxSolver::BASIC;
            y[m_i] = val / aij;
            r[m_j] = 0.0;
         }
         else
         {
            // set x_j on one of the bound
            assert(EQrel(m_oldLo, x[m_j], eps()) || EQrel(m_oldUp, x[m_j], eps()));

            cStatus[m_j] = EQrel(m_oldLo, x[m_j], eps()) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
            rStatus[m_i] = SPxSolver::BASIC;
            y[m_i] = m_row_obj;
            r[m_j] = val;
         }
      }
      else if(EQrel(newLo, m_oldUp, eps()))
      {
         // row is in the type  xj >= b/aij, try to set xj on upper
         if(r[m_j] >= eps())
         {
            // the reduced cost is positive, xj should in the basic
            assert(EQrel(m_rhs, x[m_j]*aij, eps()) || EQrel(m_lhs, x[m_j]*aij, eps()));

            rStatus[m_i] = (EQrel(m_lhs, x[m_j]*aij, eps())) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
            cStatus[m_j] = SPxSolver::BASIC;
            y[m_i] = val / aij;
            r[m_j] = 0.0;
         }
         else
         {
            assert(EQrel(m_oldUp, x[m_j], eps()));

            cStatus[m_j] = SPxSolver::ON_UPPER;
            rStatus[m_i] = SPxSolver::BASIC;
            y[m_i] = m_row_obj;
            r[m_j] = val;
         }
      }
      else if(EQrel(newUp, m_oldLo, eps()))
      {
         // row is in the type  xj <= b/aij, try to set xj on lower
         if(r[m_j] <= -eps())
         {
            // the reduced cost is negative, xj should in the basic
            assert(EQrel(m_rhs, x[m_j]*aij, eps()) || EQrel(m_lhs, x[m_j]*aij, eps()));

            rStatus[m_i] = (EQrel(m_lhs, x[m_j]*aij, eps())) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
            cStatus[m_j] = SPxSolver::BASIC;
            y[m_i] = val / aij;
            r[m_j] = 0.0;
         }
         else
         {
            assert(EQrel(m_oldLo, x[m_j], eps()));

            cStatus[m_j] = SPxSolver::ON_LOWER;
            rStatus[m_i] = SPxSolver::BASIC;
            y[m_i] = m_row_obj;
            r[m_j] = val;
         }
      }
      else
      {
         // the variable is set to FIXED by other constraints, i.e., this singleton row is redundant
         rStatus[m_i] = SPxSolver::BASIC;
         y[m_i] = m_row_obj;
      }
      break;
   case SPxSolver::BASIC:
      rStatus[m_i] = SPxSolver::BASIC;
      y[m_i] = m_row_obj;
      r[m_j] = 0.0;
      break;
   case SPxSolver::ON_LOWER:
      if(EQrel(m_oldLo, x[m_j], eps())) // xj may stay on lower
      {
         rStatus[m_i] = SPxSolver::BASIC;
         y[m_i] = m_row_obj;
         r[m_j] = val;
      }
      else // if reduced costs are negative or old lower bound not equal to xj, we need to change xj into the basis
      {
         assert(EQrel(m_rhs, x[m_j]*aij, eps()) || EQrel(m_lhs, x[m_j]*aij, eps()));

         cStatus[m_j] = SPxSolver::BASIC;
         rStatus[m_i] = (EQrel(m_lhs, x[m_j]*aij, eps())) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
         y[m_i] = val / aij;
         r[m_j] = 0.0;
      }
      break;
   case SPxSolver::ON_UPPER:
      if(EQrel(m_oldUp, x[m_j], eps())) // xj may stay on upper
      {
         rStatus[m_i] = SPxSolver::BASIC;
         y[m_i] = m_row_obj;
         r[m_j] = val;
      }
      else // if reduced costs are positive or old upper bound not equal to xj, we need to change xj into the basis
      {
         assert(EQrel(m_rhs, x[m_j]*aij, eps()) || EQrel(m_lhs, x[m_j]*aij, eps()));

         cStatus[m_j] = SPxSolver::BASIC;
         rStatus[m_i] = (EQrel(m_lhs, x[m_j]*aij, eps())) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
         y[m_i] = val / aij;
         r[m_j] = 0.0;
      }
      break;
   case SPxSolver::ZERO:
      rStatus[m_i] = SPxSolver::BASIC;
      y[m_i] = m_row_obj;
      r[m_j] = val;
      break;
   default:
      break;
   }

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
       throw SPxInternalCodeException("XMAISM17 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::ForceConstraintPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                           DataArray<SPxSolver::VarStatus>& cStatus,
                                           DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // correcting the change of idx by deletion of the row:
   s[m_old_i] = s[m_i];
   y[m_old_i] = y[m_i];
   rStatus[m_old_i] = rStatus[m_i];

   // primal:
   s[m_i] = m_lRhs;

   // basis:
   int cBasisCandidate = -1;
   Real maxViolation = -1.0;
   int bas_k = -1;

   for(int k = 0; k < m_row.size(); ++k)
   {
      int  cIdx  = m_row.index(k);
      Real aij   = m_row.value(k);
      Real oldLo = m_oldLowers[k];
      Real oldUp = m_oldUppers[k];

      switch(cStatus[cIdx])
      {
      case SPxSolver::FIXED:
         if(m_fixed[k])
         {
            assert(EQrel(oldLo, x[cIdx], eps()) || EQrel(oldUp, x[cIdx], eps()));

            Real violation = spxAbs(r[cIdx]/aij);

            cStatus[cIdx] = EQrel(oldLo, x[cIdx], eps()) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;

            if( violation > maxViolation && ( (EQrel(oldLo, x[cIdx], eps()) && r[cIdx] < -eps()) || (EQrel(oldUp, x[cIdx], eps()) && r[cIdx] > eps()) ) )
            {
               maxViolation = violation;
               cBasisCandidate = cIdx;
               bas_k = k;
            }
         } // do nothing, if the old bounds are equal, i.e. variable has been not fixed in this row
         break;
      case SPxSolver::ON_LOWER:
      case SPxSolver::ON_UPPER:
      case SPxSolver::BASIC:
         break;
      default:
         break;
      }
   }

   // dual and basis :
   if(cBasisCandidate >= 0)  // one of the variable in the row should in the basis
   {
      assert(EQrel(m_lRhs, m_rhs, eps()) || EQrel(m_lRhs, m_lhs, eps()));
      assert(bas_k >= 0);
      assert(cBasisCandidate == m_row.index(bas_k));

      cStatus[cBasisCandidate] = SPxSolver::BASIC;
      rStatus[m_i] = (EQrel(m_lRhs, m_lhs, eps())) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;

      Real aij = m_row.value(bas_k);
      Real multiplier = r[cBasisCandidate]/aij;
      r[cBasisCandidate] = 0.0;

      for(int k = 0; k < m_row.size(); ++k)  // update the reduced cost
      {
         if(k == bas_k)
         {
            continue;
         }
         r[m_row.index(k)] -= m_row.value(k) * multiplier;
      }

      // compute the value of new dual variable (because we have a new row)
      Real val = m_objs[bas_k];
      DSVector basis_col = m_cols[bas_k];

      for(int k = 0; k < basis_col.size(); ++k)
      {
         if (basis_col.index(k) != m_i)
            val -= basis_col.value(k) * y[basis_col.index(k)];
      }

      y[m_i] = val/aij;
   }
   else // slack in the basis
   {
      rStatus[m_i] = SPxSolver::BASIC;
      y[m_i] = m_rowobj;
   }

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM18 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::FixVariablePS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                       DataArray<SPxSolver::VarStatus>& cStatus,
                                       DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // update the index mapping; if m_correctIdx is false, we assume that this has happened already
   if(m_correctIdx)
   {
      x[m_old_j] = x[m_j];
      r[m_old_j] = r[m_j];
      cStatus[m_old_j] = cStatus[m_j];
   }

   // primal:
   x[m_j] = m_val;

   for(int k = 0; k < m_col.size(); ++k)
      s[m_col.index(k)] += m_col.value(k) * x[m_j];

   // dual:
   Real val = m_obj;

   for(int k = 0; k < m_col.size(); ++k)
      val -= m_col.value(k) * y[m_col.index(k)];

   r[m_j] = val;

   // basis:
   if( m_lower == m_upper )
   {
      assert(EQrel(m_lower, m_val));

      cStatus[m_j] = SPxSolver::FIXED;
   }
   else
   {
      assert(EQrel(m_val, m_lower) || EQrel(m_val, m_upper) || m_val == 0.0);

      cStatus[m_j] = EQrel(m_val, m_lower) ? SPxSolver::ON_LOWER : (EQrel(m_val, m_upper) ? SPxSolver::ON_UPPER : SPxSolver::ZERO);
   }

#ifdef CHECK_BASIC_DIM
   if(m_correctIdx)
   {
      if (!checkBasisDim(rStatus, cStatus))
      {
         throw SPxInternalCodeException("XMAISM19 Dimension doesn't match after this step.");
      }
   }
#endif
}

void SPxMainSM::FixBoundsPS::execute(DVector&, DVector&, DVector&, DVector&,
                                     DataArray<SPxSolver::VarStatus>& cStatus,
                                     DataArray<SPxSolver::VarStatus>&, bool isOptimal) const
{
   // basis:
   cStatus[m_j] = m_status;
}

void SPxMainSM::FreeZeroObjVariablePS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                               DataArray<SPxSolver::VarStatus>& cStatus,
                                               DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // correcting the change of idx by deletion of the column and corresponding rows:
   x[m_old_j] = x[m_j];
   r[m_old_j] = r[m_j];
   cStatus[m_old_j] = cStatus[m_j];

   int rIdx = m_old_i - m_col.size() + 1;

   for(int k = 0; k < m_col.size(); ++k)
   {
      int rIdx_new = m_col.index(k);
      s[rIdx] = s[rIdx_new];
      y[rIdx] = y[rIdx_new];
      rStatus[rIdx] = rStatus[rIdx_new];
      rIdx++;
   }

   // primal:
   int      domIdx = -1;
   DSVector slack(m_col.size());

   if (m_loFree)
   {
      Real minRowUp = infinity;

      for(int k = 0; k < m_rows.size(); ++k)
      {
         Real           val = 0.0;
         const SVector& row = m_rows[k];

         for(int l = 0; l < row.size(); ++l)
         {
            if (row.index(l) != m_j)
               val += row.value(l) * x[row.index(l)];
         }

         Real scale = maxAbs(m_lRhs[k], val);

         if (scale < 1.0)
            scale = 1.0;

         Real z = (m_lRhs[k] / scale) - (val / scale);

         if (isZero(z))
            z = 0.0;

         Real up = z * scale / row[m_j];
         slack.add(k, val);

         if (up < minRowUp)
         {
            minRowUp = up;
            domIdx   = k;
         }
      }

      if (m_bnd < minRowUp)
      {
         x[m_j] = m_bnd;
         domIdx = -1;
      }
      else
         x[m_j] = minRowUp;
   }
   else
   {
      Real maxRowLo = -infinity;

      for(int k = 0; k < m_rows.size(); ++k)
      {
         Real val = 0.0;
         const SVector& row = m_rows[k];

         for(int l = 0; l < row.size(); ++l)
         {
            if (row.index(l) != m_j)
               val += row.value(l) * x[row.index(l)];
         }

         Real scale = maxAbs(m_lRhs[k], val);

         if (scale < 1.0)
            scale = 1.0;

         Real z = (m_lRhs[k] / scale) - (val / scale);

         if (isZero(z))
            z = 0.0;

         Real lo = z * scale / row[m_j];
         slack.add(k, val);

         if (lo > maxRowLo)
         {
            maxRowLo = lo;
            domIdx   = k;
         }
      }

      if (m_bnd > maxRowLo)
      {
         x[m_j] = m_bnd;
         domIdx = -1;
      }
      else
         x[m_j] = maxRowLo;
   }

   for(int k = 0; k < m_col.size(); ++k)
      s[m_col.index(k)] = slack[k] + m_col.value(k) * x[m_j];

   // dual:
   r[m_j] = 0.0;

   for(int k = 0; k < m_col.size(); ++k)
   {
      int idx = m_col.index(k);
      y[idx] = m_rowObj[idx];
   }

   // basis:
   for(int k = 0; k < m_col.size(); ++k)
   {
      if (k != domIdx)
         rStatus[m_col.index(k)] = SPxSolver::BASIC;

      else
      {
         cStatus[m_j] = SPxSolver::BASIC;
         if (m_loFree)
            rStatus[m_col.index(k)] = (m_col.value(k) > 0) ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
         else
            rStatus[m_col.index(k)] = (m_col.value(k) > 0) ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
      }
   }
   if (domIdx == -1)
   {
      if (m_loFree)
         cStatus[m_j] = SPxSolver::ON_UPPER;
      else
         cStatus[m_j] = SPxSolver::ON_LOWER;
   }

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM20 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::ZeroObjColSingletonPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                               DataArray<SPxSolver::VarStatus>& cStatus,
                                               DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // correcting the change of idx by deletion of the column and corresponding rows:
   x[m_old_j] = x[m_j];
   r[m_old_j] = r[m_j];
   cStatus[m_old_j] = cStatus[m_j];

   // primal & basis:
   Real aij = m_row[m_j];

   if (isZero(s[m_i], 1e-6))
      s[m_i] = 0.0;
   else if (s[m_i] >= soplex::infinity)
      // this is a fix for a highly ill conditioned instance that is "solved" in presolving (ilaser0 from MINLP, mittelmann)
      throw SPxException("Simplifier: infinite activities - aborting unsimplification");

   Real scale1 = maxAbs(m_lhs, s[m_i]);
   Real scale2 = maxAbs(m_rhs, s[m_i]);

   if (scale1 < 1.0)
      scale1 = 1.0;
   if (scale2 < 1.0)
      scale2 = 1.0;

   Real z1 = (m_lhs / scale1) - (s[m_i] / scale1);
   Real z2 = (m_rhs / scale2) - (s[m_i] / scale2);

   if (isZero(z1))
      z1 = 0.0;
   if (isZero(z2))
      z2 = 0.0;

   Real lo = (aij > 0) ? z1 * scale1 / aij : z2 * scale2 / aij;
   Real up = (aij > 0) ? z2 * scale2 / aij : z1 * scale1 / aij;

   if (isZero(lo, eps()))
      lo = 0.0;
   if (isZero(up, eps()))
      up = 0.0;

   assert(LErel(lo, up));
   ASSERT_WARN( "WMAISM01", isNotZero(aij, 1.0 / infinity) );

   if (rStatus[m_i] == SPxSolver::ON_LOWER)
   {
      if ( m_lower <= -infinity && m_upper >= infinity )
      {
         x[m_j] = 0.0;
         cStatus[m_j] = SPxSolver::ZERO;
      }
      else if ( m_lower == m_upper )
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::FIXED;
      }
      else if (aij > 0)
      {
         x[m_j]       = m_upper;
         cStatus[m_j] = SPxSolver::ON_UPPER;
      }
      else if (aij < 0)
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::ON_LOWER;
      }
      else
         throw SPxInternalCodeException("XMAISM01 This should never happen.");
   }
   else if (rStatus[m_i] == SPxSolver::ON_UPPER)
   {
      if ( m_lower <= -infinity && m_upper >= infinity )
      {
         x[m_j] = 0.0;
         cStatus[m_j] = SPxSolver::ZERO;
      }
      else if ( m_lower == m_upper )
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::FIXED;
      }
      else if (aij > 0)
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = SPxSolver::ON_LOWER;
      }
      else if (aij < 0)
      {
         x[m_j]       = m_upper;
         cStatus[m_j] = SPxSolver::ON_UPPER;
      }
      else
         throw SPxInternalCodeException("XMAISM02 This should never happen.");
   }
   else if (rStatus[m_i] == SPxSolver::FIXED)
   {
      if ( m_lower <= -infinity && m_upper >= infinity )
      {
         x[m_j] = 0.0;
         cStatus[m_j] = SPxSolver::ZERO;
      }
      else
      {
         assert(EQrel(m_lower, m_upper));

         x[m_j]        = m_lower;
         cStatus[m_j]  = SPxSolver::FIXED;
      }
   }
   else if (rStatus[m_i] == SPxSolver::BASIC)
   {
      if (GErel(m_lower, lo, eps()) && m_lower > -infinity)
      {
         x[m_j]       = m_lower;
         cStatus[m_j] = (m_lower == m_upper) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
      else if (LErel(m_upper, up, eps()) && m_upper < infinity)
      {
         x[m_j]       = m_upper;
         cStatus[m_j] = (m_lower == m_upper) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
      else if (lo > -infinity)
      {
         // make m_i non-basic and m_j basic
         x[m_j]       = lo;
         cStatus[m_j] = SPxSolver::BASIC;
         rStatus[m_i] = (aij > 0 ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER);
      }
      else if (up < infinity)
      {
         // make m_i non-basic and m_j basic
         x[m_j]       = up;
         cStatus[m_j] = SPxSolver::BASIC;
         rStatus[m_i] = (aij > 0 ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER);
      }
      else
         throw SPxInternalCodeException("XMAISM03 This should never happen.");
   }
   else
      throw SPxInternalCodeException("XMAISM04 This should never happen.");

   s[m_i] += aij * x[m_j];

   // dual:
   r[m_j] = -1.0 * aij * y[m_i];

   assert(!isOptimal || (cStatus[m_j] != SPxSolver::BASIC || isZero(r[m_j], eps())));

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM21 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::FreeColSingletonPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                            DataArray<SPxSolver::VarStatus>& cStatus,
                                            DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{

   // correcting the change of idx by deletion of the row:
   s[m_old_i] = s[m_i];
   y[m_old_i] = y[m_i];
   rStatus[m_old_i] = rStatus[m_i];

   // correcting the change of idx by deletion of the column:
   x[m_old_j] = x[m_j];
   r[m_old_j] = r[m_j];
   cStatus[m_old_j] = cStatus[m_j];

   // primal:
   Real val = 0.0;
   Real aij = m_row[m_j];

   for(int k = 0; k < m_row.size(); ++k)
   {
      if (m_row.index(k) != m_j)
         val += m_row.value(k) * x[m_row.index(k)];
   }

   Real scale = maxAbs(m_lRhs, val);

   if (scale < 1.0)
      scale = 1.0;

   Real z = (m_lRhs / scale) - (val / scale);

   if (isZero(z))
      z = 0.0;

   x[m_j] = z * scale / aij;
   s[m_i] = m_lRhs;

   // dual:
   y[m_i] = m_obj / aij;
   r[m_j] = 0.0;

   // basis:
   cStatus[m_j] = SPxSolver::BASIC;

   if (m_eqCons)
      rStatus[m_i] = SPxSolver::FIXED;
   else if (m_onLhs)
      rStatus[m_i] = SPxSolver::ON_LOWER;
   else
      rStatus[m_i] = SPxSolver::ON_UPPER;

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM22 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::DoubletonEquationPS::execute(DVector& x, DVector& y, DVector&, DVector& r,
                                             DataArray<SPxSolver::VarStatus>& cStatus,
                                             DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // dual:
   if ((cStatus[m_k]  != SPxSolver::BASIC) &&
       ((cStatus[m_k] == SPxSolver::ON_LOWER && m_strictLo) ||
        (cStatus[m_k] == SPxSolver::ON_UPPER && m_strictUp) ||
        (cStatus[m_k] == SPxSolver::FIXED    &&
         (( m_maxSense && ((r[m_j] > 0 && m_strictUp) || (r[m_j] < 0 && m_strictLo))) ||
          (!m_maxSense && ((r[m_j] > 0 && m_strictLo) || (r[m_j] < 0 && m_strictUp)))))))
   {
      Real val  = m_kObj;
      Real aik  = m_col[m_i];

      for(int _k = 0; _k < m_col.size(); ++_k)
      {
         if (m_col.index(_k) != m_i)
            val -= m_col.value(_k) * y[m_col.index(_k)];
      }

      y[m_i] = val / aik;
      r[m_k] = 0.0;

      r[m_j] = m_jObj - val * m_aij / aik;

      ASSERT_WARN( "WMAISM73", isNotZero(m_aij * aik) );

      // basis:
      if( m_jFixed)
         cStatus[m_j] = SPxSolver::FIXED;
      else
      {
         if( GT(r[m_j], 0) || (isZero(r[m_j]) && EQ(x[m_j], m_Lo_j)) )
            cStatus[m_j] = SPxSolver::ON_LOWER;
         else
            cStatus[m_j] = SPxSolver::ON_UPPER;
      }

      cStatus[m_k] = SPxSolver::BASIC;
   }

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM23 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::DuplicateRowsPS::execute(DVector&, DVector& y, DVector& s, DVector&,
                                         DataArray<SPxSolver::VarStatus>& cStatus,
                                         DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // correcting the change of idx by deletion of the duplicated rows:
   if(m_isLast)
   {
      for(int i = m_perm.size() - 1; i >= 0; --i)
      {
         if (m_perm[i] >= 0)
         {
            int rIdx_new = m_perm[i];
            int rIdx = i;
            s[rIdx] = s[rIdx_new];
            y[rIdx] = y[rIdx_new];
            rStatus[rIdx] = rStatus[rIdx_new];
         }
      }
   }

   // primal:
   for(int k = 0; k < m_scale.size(); ++k)
   {
      if (m_scale.index(k) != m_i)
         s[m_scale.index(k)] = s[m_i] / m_scale.value(k);
   }

   // dual & basis:
   bool haveSetBasis = false;

   for(int k = 0; k < m_scale.size(); ++k)
   {
      int i = m_scale.index(k);

      if (rStatus[m_i] == SPxSolver::BASIC || (haveSetBasis && i!=m_i))
         // if the row with tightest lower and upper bound in the basic, every duplicate row should in basic
         // or basis status of row m_i has been set, this row should be in basis
      {
         y[i]       = m_rowObj.value(k);
         rStatus[i] = SPxSolver::BASIC;
         continue;
      }

      ASSERT_WARN( "WMAISM02", isNotZero(m_scale.value(k)) );

      if (rStatus[m_i] == SPxSolver::FIXED && (i == m_maxLhsIdx || i == m_minRhsIdx))
      {
         // this row leads to the tightest lower or upper bound, slack should not be in the basis
         y[i]   = y[m_i] * m_scale.value(k);
         y[m_i] = m_i_rowObj;

         if(m_isLhsEqualRhs[k])
         {
            rStatus[i] = SPxSolver::FIXED;
         }
         else if(i == m_maxLhsIdx)
         {
            rStatus[i] = m_scale.value(k)*m_scale.value(0) > 0 ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
         }
         else
         {
            assert(i == m_minRhsIdx);

            rStatus[i] = m_scale.value(k)*m_scale.value(0) > 0 ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
         }
         if (i != m_i)
            rStatus[m_i] = SPxSolver::BASIC;
         haveSetBasis = true;
      }
      else if (i == m_maxLhsIdx && rStatus[m_i] == SPxSolver::ON_LOWER)
      {
         // this row leads to the tightest lower bound, slack should not be in the basis
         y[i]   = y[m_i] * m_scale.value(k);
         y[m_i] = m_i_rowObj;

         rStatus[i] = m_scale.value(k)*m_scale.value(0) > 0 ? SPxSolver::ON_LOWER : SPxSolver::ON_UPPER;
         if (i != m_i)
            rStatus[m_i] = SPxSolver::BASIC;
         haveSetBasis = true;
      }
      else if (i == m_minRhsIdx && rStatus[m_i] == SPxSolver::ON_UPPER)
      {
         // this row leads to the tightest upper bound, slack should not be in the basis
         y[i]   = y[m_i] * m_scale.value(k);
         y[m_i] = m_i_rowObj;

         rStatus[i] = m_scale.value(k)*m_scale.value(0) > 0 ? SPxSolver::ON_UPPER : SPxSolver::ON_LOWER;
         if (i != m_i)
            rStatus[m_i] = SPxSolver::BASIC;
         haveSetBasis = true;
      }
      else if (i != m_i)
      {
         // this row does not lead to the tightest lower or upper bound, slack should be in the basis
         y[i]       = m_rowObj.value(k);
         rStatus[i] = SPxSolver::BASIC;
      }
   }

#ifdef CHECK_BASIC_DIM
   if(m_isFirst && !checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM24 Dimension doesn't match after this step.");
   }
#endif

   // nothing to do for the reduced cost values
}

void SPxMainSM::DuplicateColsPS::execute(DVector& x,
                                         DVector&,
                                         DVector&,
                                         DVector& r,
                                         DataArray<SPxSolver::VarStatus>& cStatus,
                                         DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{

   if(m_isFirst)
   {
#ifdef CHECK_BASIC_DIM
      if (!checkBasisDim(rStatus, cStatus))
      {
         throw SPxInternalCodeException("XMAISM25 Dimension doesn't match after this step.");
      }
#endif
      return;
   }


   // correcting the change of idx by deletion of the columns:
   if(m_isLast)
   {
      for(int i = m_perm.size() - 1; i >= 0; --i)
      {
         if (m_perm[i] >= 0)
         {
            int cIdx_new = m_perm[i];
            int cIdx = i;
            x[cIdx] = x[cIdx_new];
            r[cIdx] = r[cIdx_new];
            cStatus[cIdx] = cStatus[cIdx_new];
         }
      }
      return;
   }

   // primal & basis:
   ASSERT_WARN( "WMAISM03", isNotZero(m_scale) );

   if (cStatus[m_k] == SPxSolver::ON_LOWER)
   {
      x[m_k] = m_loK;

      if (m_scale > 0)
      {
         x[m_j]       = m_loJ;
         cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
      else
      {
         x[m_j]       = m_upJ;
         cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
   }
   else if (cStatus[m_k] == SPxSolver::ON_UPPER)
   {
      x[m_k] = m_upK;

      if (m_scale > 0)
      {
         x[m_j]       = m_upJ;
         cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
      }
      else
      {
         x[m_j]       = m_loJ;
         cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
      }
   }
   else if (cStatus[m_k] == SPxSolver::FIXED)
   {
      // => x[m_k] and x[m_j] are also fixed before the corresponding preprocessing step
      x[m_j]       = m_loJ;
      cStatus[m_j] = SPxSolver::FIXED;
   }
   else if (cStatus[m_k] == SPxSolver::ZERO)
   {
      /* we only aggregate duplicate columns if 0 is contained in their bounds, so we can handle this case properly */
      assert(isZero(x[m_k]));
      assert(LErel(m_loJ, 0.0));
      assert(GErel(m_upJ, 0.0));
      assert(LErel(m_loK, 0.0));
      assert(GErel(m_upK, 0.0));

      if (isZero(m_loK) && isZero(m_upK) && m_loK == m_upK)
         cStatus[m_k] = SPxSolver::FIXED;
      else if (isZero(m_loK))
         cStatus[m_k] = SPxSolver::ON_LOWER;
      else if (isZero(m_upK))
         cStatus[m_k] = SPxSolver::ON_UPPER;
      else if (LErel(m_loK, 0.0) && GErel(m_upK, 0.0))
         cStatus[m_k] = SPxSolver::ZERO;
      else
         throw SPxInternalCodeException("XMAISM05 This should never happen.");

      x[m_j] = 0.0;
      if (isZero(m_loJ) && isZero(m_upJ) && m_loJ == m_upJ)
         cStatus[m_j] = SPxSolver::FIXED;
      else if (isZero(m_loJ))
         cStatus[m_j] = SPxSolver::ON_LOWER;
      else if (isZero(m_upJ))
         cStatus[m_j] = SPxSolver::ON_UPPER;
      else if (LErel(m_loJ, 0.0) && GErel(m_upJ, 0.0))
            cStatus[m_j] = SPxSolver::ZERO;
      else
         throw SPxInternalCodeException("XMAISM06 This should never happen.");
   }
   else if (cStatus[m_k] == SPxSolver::BASIC)
   {
      Real scale1 = maxAbs(x[m_k], m_loK);
      Real scale2 = maxAbs(x[m_k], m_upK);

      if (scale1 < 1.0)
         scale1 = 1.0;
      if (scale2 < 1.0)
         scale2 = 1.0;

      Real z1 = (x[m_k] / scale1) - (m_loK / scale1);
      Real z2 = (x[m_k] / scale2) - (m_upK / scale2);

      if (isZero(z1))
         z1 = 0.0;
      if (isZero(z2))
         z2 = 0.0;

      if( m_loJ <= -infinity && m_upJ >= infinity && m_loK <= -infinity && m_upK >= infinity )
      {
         cStatus[m_j] = SPxSolver::ZERO;
         x[m_j] = 0.0;
      }
      else if( m_scale > 0.0 )
      {
         if( GErel(x[m_k], m_upK + m_scale * m_upJ) )
         {
            assert(m_upJ < infinity);
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
            x[m_j] = m_upJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else if( GErel(x[m_k], m_loK + m_scale * m_upJ) && m_upJ < infinity )
         {
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
            x[m_j] = m_upJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else if( GErel(x[m_k], m_upK + m_scale * m_loJ) && m_upK < infinity )
         {
            cStatus[m_k] = (m_loK == m_upK) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
            x[m_k] = m_upK;
            cStatus[m_j] = SPxSolver::BASIC;
            x[m_j] = z2 * scale2 / m_scale;
         }
         else if( GErel(x[m_k], m_loK + m_scale * m_loJ) && m_loJ > -infinity )
         {
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
            x[m_j] = m_loJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else if( GErel(x[m_k], m_loK + m_scale * m_loJ) && m_loK > -infinity )
         {
            cStatus[m_k] = (m_loK == m_upK) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
            x[m_k] = m_loK;
            cStatus[m_j] = SPxSolver::BASIC;
            x[m_j] = z1 * scale1 / m_scale;
         }
         else if( LTrel(x[m_k], m_loK + m_scale * m_loJ) )
         {
            assert(m_loJ > -infinity);
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
            x[m_j] = m_loJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else
         {
            throw SPxInternalCodeException("XMAISM08 This should never happen.");
         }
      }
      else
      {
         assert(m_scale < 0.0);

         if( GErel(x[m_k], m_upK + m_scale * m_loJ) )
         {
            assert(m_loJ > -infinity);
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
            x[m_j] = m_loJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else if( GErel(x[m_k], m_loK + m_scale * m_loJ) && m_loJ > -infinity )
         {
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
            x[m_j] = m_loJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else if( GErel(x[m_k], m_upK + m_scale * m_upJ) && m_upK < infinity )
         {
            cStatus[m_k] = (m_loK == m_upK) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
            x[m_k] = m_upK;
            cStatus[m_j] = SPxSolver::BASIC;
            x[m_j] = z2 * scale2 / m_scale;
         }
         else if( GErel(x[m_k], m_loK + m_scale * m_upJ) && m_upJ < infinity )
         {
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
            x[m_j] = m_upJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else if( GErel(x[m_k], m_loK + m_scale * m_upJ) && m_loK > -infinity )
         {
            cStatus[m_k] = (m_loK == m_upK) ? SPxSolver::FIXED : SPxSolver::ON_LOWER;
            x[m_k] = m_loK;
            cStatus[m_j] = SPxSolver::BASIC;
            x[m_j] = z1 * scale1 / m_scale;
         }
         else if( LTrel(x[m_k], m_loK + m_scale * m_upJ) )
         {
            assert(m_upJ < infinity);
            cStatus[m_j] = (m_loJ == m_upJ) ? SPxSolver::FIXED : SPxSolver::ON_UPPER;
            x[m_j] = m_upJ;
            x[m_k] -= m_scale * x[m_j];
         }
         else
         {
            throw SPxInternalCodeException("XMAISM09 This should never happen.");
         }
      }
   }

   // dual:
   r[m_j] = m_scale * r[m_k];
}

void SPxMainSM::MultiAggregationPS::execute(DVector& x, DVector& y, DVector& s, DVector& r,
                                            DataArray<SPxSolver::VarStatus>& cStatus,
                                            DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{

   // correcting the change of idx by deletion of the row:
   s[m_old_i] = s[m_i];
   y[m_old_i] = y[m_i];
   rStatus[m_old_i] = rStatus[m_i];

   // correcting the change of idx by deletion of the column:
   x[m_old_j] = x[m_j];
   r[m_old_j] = r[m_j];
   cStatus[m_old_j] = cStatus[m_j];

   // primal:
   Real val = 0.0;
   Real aij = m_row[m_j];

   for(int k = 0; k < m_row.size(); ++k)
   {
      if(m_row.index(k) != m_j)
         val += m_row.value(k) * x[m_row.index(k)];
   }

   Real scale = maxAbs(m_const, val);

   if (scale < 1.0)
      scale = 1.0;

   Real z = (m_const / scale) - (val / scale);

   if (isZero(z))
      z = 0.0;

   x[m_j] = z * scale / aij;
   s[m_i] = 0.0;

   assert(!isOptimal || (GE(x[m_j], m_lower) && LE(x[m_j], m_upper)));

   // dual:
   Real dualVal = 0.0;

   for(int k = 0; k < m_col.size(); ++k)
   {
      if(m_col.index(k) != m_i)
         dualVal += m_col.value(k) * y[m_col.index(k)];
   }

   z = m_obj - dualVal;

   y[m_i] = z / aij;
   r[m_j] = 0.0;

   // basis:
   cStatus[m_j] = SPxSolver::BASIC;

   if (m_eqCons)
      rStatus[m_i] = SPxSolver::FIXED;
   else if (m_onLhs)
      rStatus[m_i] = SPxSolver::ON_LOWER;
   else
      rStatus[m_i] = SPxSolver::ON_UPPER;

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM22 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::TightenBoundsPS::execute(DVector& x, DVector&, DVector&, DVector&,
                                     DataArray<SPxSolver::VarStatus>& cStatus,
                                     DataArray<SPxSolver::VarStatus>& rStatus, bool isOptimal) const
{
   // basis:
   switch(cStatus[m_j])
   {
      case SPxSolver::FIXED:
         if(LT(x[m_j], m_origupper, eps()) && GT(x[m_j], m_origlower, eps()))
            cStatus[m_j] = SPxSolver::BASIC;
         else if(LT(x[m_j], m_origupper, eps()))
            cStatus[m_j] = SPxSolver::ON_LOWER;
         else if(GT(x[m_j], m_origlower, eps()))
            cStatus[m_j] = SPxSolver::ON_UPPER;

         break;
      case SPxSolver::ON_LOWER:
         if(GT(x[m_j], m_origlower, eps()))
            cStatus[m_j] = SPxSolver::BASIC;

         break;
      case SPxSolver::ON_UPPER:
         if(LT(x[m_j], m_origupper, eps()))
            cStatus[m_j] = SPxSolver::BASIC;

         break;
      default:
         break;
   }

#ifdef CHECK_BASIC_DIM
   if (!checkBasisDim(rStatus, cStatus))
   {
      throw SPxInternalCodeException("XMAISM22 Dimension doesn't match after this step.");
   }
#endif
}

void SPxMainSM::handleRowObjectives(SPxLP& lp)
{
   for( int i = lp.nRows() - 1; i >= 0; --i )
   {
      if( lp.maxRowObj(i) != 0.0 )
      {
         RowObjPS* RowObjPSptr = 0;
         spx_alloc(RowObjPSptr);
         m_hist.append(new (RowObjPSptr) RowObjPS(lp, i, lp.nCols()));
         lp.addCol(lp.rowObj(i), -lp.rhs(i), UnitVector(i), -lp.lhs(i));
         lp.changeRange(i, 0.0, 0.0);
         lp.changeRowObj(i, 0.0);
         m_addedcols++;
      }
   }
}

void SPxMainSM::handleExtremes(SPxLP& lp)
{

   // This method handles extreme value of the given LP by
   //
   // 1. setting numbers of very small absolute values to zero and
   // 2. setting numbers of very large absolute values to -infinity or +infinity, respectively.

   Real maxVal  = infinity / 5.0;
   Real tol = feastol() * 1e-2;
   tol = (tol < epsZero()) ? epsZero() : tol;
   int  remRows = 0;
   int  remNzos = 0;
   int  chgBnds = 0;
   int  chgLRhs = 0;
   int  objCnt  = 0;

   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      // lhs
      Real lhs = lp.lhs(i);

      if (lhs != 0.0 && isZero(lhs, epsZero()))
      {
         lp.changeLhs(i, 0.0);
         ++chgLRhs;
      }
      else if (lhs > -infinity && lhs < -maxVal)
      {
         lp.changeLhs(i, -infinity);
         ++chgLRhs;
      }
      else if (lhs <  infinity && lhs >  maxVal)
      {
         lp.changeLhs(i,  infinity);
         ++chgLRhs;
      }

      // rhs
      Real rhs = lp.rhs(i);

      if (rhs != 0.0 && isZero(rhs, epsZero()))
      {
         lp.changeRhs(i, 0.0);
         ++chgLRhs;
      }
      else if (rhs > -infinity && rhs < -maxVal)
      {
         lp.changeRhs(i, -infinity);
         ++chgLRhs;
      }
      else if (rhs <  infinity && rhs >  maxVal)
      {
         lp.changeRhs(i,  infinity);
         ++chgLRhs;
      }

      if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
      {
         FreeConstraintPS* FreeConstraintPSptr = 0;
         spx_alloc(FreeConstraintPSptr);
         m_hist.append(new (FreeConstraintPSptr) FreeConstraintPS(lp, i));

         removeRow(lp, i);
         ++remRows;

         ++m_stat[FREE_ROW];
      }
   }

   for(int j = 0; j < lp.nCols(); ++j)
   {
      // lower bound
      Real lo = lp.lower(j);

      if (lo != 0.0 && isZero(lo, epsZero()))
      {
         lp.changeLower(j, 0.0);
         ++chgBnds;
      }
      else if (lo > -infinity && lo < -maxVal)
      {
         lp.changeLower(j, -infinity);
         ++chgBnds;
      }
      else if (lo <  infinity && lo >  maxVal)
      {
         lp.changeLower(j,  infinity);
         ++chgBnds;
      }

      // upper bound
      Real up = lp.upper(j);

      if (up != 0.0 && isZero(up, epsZero()))
      {
         lp.changeUpper(j, 0.0);
         ++chgBnds;
      }
      else if (up > -infinity && up < -maxVal)
      {
         lp.changeUpper(j, -infinity);
         ++chgBnds;
      }
      else if (up <  infinity && up >  maxVal)
      {
         lp.changeUpper(j,  infinity);
         ++chgBnds;
      }

      // fixed columns will be eliminated later
      if (NE(lo, up))
      {
         lo = spxAbs(lo);
         up = spxAbs(up);

         Real absBnd = (lo > up) ? lo : up;

         if (absBnd < 1.0)
            absBnd = 1.0;

         // non-zeros
         SVector& col = lp.colVector_w(j);
         int        i = 0;

         while(i < col.size())
         {
            Real aij = spxAbs(col.value(i));

            if (isZero(aij * absBnd, tol))
            {
               SVector& row = lp.rowVector_w(col.index(i));
               int row_j = row.pos(j);

               // this changes col.size()
               if( row_j >= 0 )
                  row.remove(row_j);
               col.remove(i);

               MSG_DEBUG( (*spxout) << "IMAISM04 aij=" << aij
                                 << " removed, absBnd=" << absBnd
                                 << std::endl; )
               ++remNzos;
            }
            else
            {
               if (aij > maxVal)
               {
                  MSG_WARNING( (*spxout), (*spxout) << "WMAISM05 Warning! Big matrix coefficient " << aij << std::endl );
               }
               else if( isZero(aij, tol) )
               {
                  MSG_WARNING( (*spxout), (*spxout) << "WMAISM06 Warning! Tiny matrix coefficient " << aij << std::endl );
               }

               ++i;
            }
         }
      }

      // objective
      Real obj = lp.obj(j);

      if (obj != 0.0 && isZero(obj, epsZero()))
      {
         lp.changeObj(j, 0.0);
         ++objCnt;
      }
      else if (obj > -infinity && obj < -maxVal)
      {
         lp.changeObj(j, -infinity);
         ++objCnt;
      }
      else if (obj <  infinity && obj >  maxVal)
      {
         lp.changeObj(j,  infinity);
         ++objCnt;
      }
   }

   if (remRows + remNzos + chgLRhs + chgBnds + objCnt > 0)
   {
      m_remRows += remRows;
      m_remNzos += remNzos;
      m_chgLRhs += chgLRhs;
      m_chgBnds += chgBnds;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (extremes) removed "
                        << remRows << " rows, "
                        << remNzos << " non-zeros, "
                        << chgBnds << " col bounds, "
                        << chgLRhs << " row bounds, "
                        << objCnt  << " objective coefficients" << std::endl; )
   }
   assert(lp.isConsistent());
}

/// computes the minimum and maximum residual activity for a given variable
void SPxMainSM::computeMinMaxResidualActivity(SPxLP& lp, int rowNumber, int colNumber, Real& minAct, Real& maxAct)
{
   const SVector& row = lp.rowVector(rowNumber);
   bool minNegInfinite = false;
   bool maxInfinite = false;

   minAct = 0;   // this is the minimum value that the aggregation can attain
   maxAct = 0;   // this is the maximum value that the aggregation can attain

   for (int l = 0; l < row.size(); ++l)
   {
      if (colNumber < 0 || row.index(l) != colNumber)
      {
         // computing the minimum activity of the aggregated variables
         if (GT(row.value(l), 0.0))
         {
            if (LE(lp.lower(row.index(l)), -infinity))
               minNegInfinite = true;
            else
               minAct += row.value(l)*lp.lower(row.index(l));
         }
         else if (LT(row.value(l), 0.0))
         {
            if (GE(lp.upper(row.index(l)), infinity))
               minNegInfinite = true;
            else
               minAct += row.value(l)*lp.upper(row.index(l));
         }

         // computing the maximum activity of the aggregated variables
         if (GT(row.value(l), 0.0))
         {
            if (GE(lp.upper(row.index(l)), infinity))
               maxInfinite = true;
            else
               maxAct += row.value(l)*lp.upper(row.index(l));
         }
         else if (LT(row.value(l), 0.0))
         {
            if (LE(lp.lower(row.index(l)), -infinity))
               maxInfinite = true;
            else
               maxAct += row.value(l)*lp.lower(row.index(l));
         }
      }
   }

   // if an infinite value exists for the minimum activity, then that it taken
   if (minNegInfinite)
      minAct = -infinity;

   // if an -infinite value exists for the maximum activity, then that value is taken
   if (maxInfinite)
      maxAct = infinity;
}


/// calculate min/max value for the multi aggregated variables
void SPxMainSM::computeMinMaxValues(SPxLP& lp, Real side, Real val, Real minRes, Real maxRes, Real& minVal, Real& maxVal)
{
   minVal = 0;
   maxVal = 0;

   if(LT(val, 0.0))
   {
      if(LE(minRes, -infinity))
         minVal = -infinity;
      else
         minVal = (side - minRes)/val;

      if(GE(maxRes, infinity))
         maxVal = infinity;
      else
         maxVal = (side - maxRes)/val;
   }
   else if(GT(val, 0.0))
   {
      if(GE(maxRes, infinity))
         minVal = -infinity;
      else
         minVal = (side - maxRes)/val;

      if(LE(minRes, -infinity))
         maxVal = infinity;
      else
         maxVal = (side - minRes)/val;
   }
}


/// tries to find good lower bound solutions by applying some trivial heuristics
void SPxMainSM::trivialHeuristic(SPxLP& lp)
{
   DVector         zerosol(lp.nCols());  // the zero solution vector
   DVector         lowersol(lp.nCols()); // the lower bound solution vector
   DVector         uppersol(lp.nCols()); // the upper bound solution vector
   DVector         locksol(lp.nCols());  // the locks solution vector

   DVector         upLocks(lp.nCols());
   DVector         downLocks(lp.nCols());

   Real            zeroObj = m_objoffset;
   Real            lowerObj = m_objoffset;
   Real            upperObj = m_objoffset;
   Real            lockObj = m_objoffset;

   bool            zerovalid = true;

   Real largeValue = infinity;
   if(LT(1.0/feastol(), infinity))
      largeValue = 1.0/feastol();



   for(int j = lp.nCols()-1; j >= 0; --j)
   {
      upLocks[j] = 0;
      downLocks[j] = 0;

      // computing the locks on the variables
      const SVector& col = lp.colVector(j);
      for(int k = 0; k < col.size(); ++k)
      {
         Real aij = col.value(k);

         ASSERT_WARN( "WMAISM45", isNotZero(aij, 1.0 / infinity) );

         if(GT(lp.lhs(col.index(k)), -infinity) && LT(lp.rhs(col.index(k)), infinity))
         {
            upLocks[j]++;
            downLocks[j]++;
         }
         else if (GT(lp.lhs(col.index(k)), -infinity))
         {
            if(aij > 0)
               downLocks[j]++;
            else if (aij < 0)
               upLocks[j]++;
         }
         else if (LT(lp.rhs(col.index(k)), infinity))
         {
            if(aij > 0)
               upLocks[j]++;
            else if (aij < 0)
               downLocks[j]++;
         }
      }

      Real lower = lp.lower(j);
      Real upper = lp.upper(j);

      if(LE(lower, -infinity))
         lower = MINIMUM(-largeValue, upper);
      if(GE(upper, infinity))
         upper = MAXIMUM(lp.lower(j), largeValue);

      if(zerovalid)
      {
         if(LE(lower, 0.0, feastol()) && GE(upper, 0.0, feastol()))
            zerosol[j] = 0.0;
         else
            zerovalid = false;
      }

      lowersol[j] = lower;
      uppersol[j] = upper;

      if(downLocks[j] > upLocks[j])
         locksol[j] = upper;
      else if(downLocks[j] < upLocks[j])
         locksol[j] = lower;
      else
         locksol[j] = (lower + upper)/2.0;

      lowerObj += lp.maxObj(j)*lowersol[j];
      upperObj += lp.maxObj(j)*uppersol[j];
      lockObj += lp.maxObj(j)*locksol[j];
   }

   // trying the lower bound solution
   if(checkSolution(lp, lowersol))
   {
      if(lowerObj > m_cutoffbound)
         m_cutoffbound = lowerObj;
   }

   // trying the upper bound solution
   if(checkSolution(lp, uppersol))
   {
      if(upperObj > m_cutoffbound)
         m_cutoffbound = upperObj;
   }

   // trying the zero solution
   if(zerovalid && checkSolution(lp, zerosol))
   {
      if(zeroObj > m_cutoffbound)
         m_cutoffbound = zeroObj;
   }

   // trying the lock solution
   if(checkSolution(lp, locksol))
   {
      if(lockObj > m_cutoffbound)
         m_cutoffbound = lockObj;
   }
}


/// checks a solution for feasibility
bool SPxMainSM::checkSolution(SPxLP& lp, DVector sol)
{
   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      const SVector& row = lp.rowVector(i);
      Real activity = 0;

      for(int k = 0; k < row.size(); k++)
         activity += row.value(k)*sol[row.index(k)];

      if(!GE(activity, lp.lhs(i), feastol()) || !LE(activity, lp.rhs(i), feastol()))
         return false;
   }

   return true;
}



/// tightens variable bounds by propagating the pseudo objective function value.
void SPxMainSM::propagatePseudoobj(SPxLP& lp)
{
   Real pseudoObj = m_objoffset;

   for( int j = lp.nCols()-1; j >= 0; --j )
   {
      Real val = lp.maxObj(j);
      if( val < 0 )
      {
         if( lp.lower(j) <= -infinity )
            return;
         pseudoObj += val*lp.lower(j);
      }
      else if( val > 0 )
      {
         if( lp.upper(j) >= -infinity )
            return;
         pseudoObj += val*lp.upper(j);
      }
   }

   if( GT(m_cutoffbound, -infinity) && LT(m_cutoffbound, infinity) )
   {
      if(pseudoObj > m_pseudoobj)
         m_pseudoobj = pseudoObj;

      for(int j = lp.nCols()-1; j >= 0; --j)
      {
         Real objval = lp.maxObj(j);

         if(EQ(objval, 0.0))
            continue;

         if(objval < 0.0)
         {
            Real newbound = lp.lower(j) + (m_cutoffbound - m_pseudoobj) / objval;

            if(LT(newbound, lp.upper(j)))
            {
               TightenBoundsPS* TightenBoundsPSptr = 0;
               spx_alloc(TightenBoundsPSptr);
               m_hist.append(new (TightenBoundsPSptr) TightenBoundsPS(lp, j, lp.upper(j), lp.lower(j)));
               lp.changeUpper(j, newbound);
            }
         }
         else if(objval > 0.0)
         {
            Real newbound = lp.upper(j) + (m_cutoffbound - m_pseudoobj) / objval;
            if(GT(newbound, lp.lower(j)))
            {
               TightenBoundsPS* TightenBoundsPSptr = 0;
               spx_alloc(TightenBoundsPSptr);
               m_hist.append(new (TightenBoundsPSptr) TightenBoundsPS(lp, j, lp.upper(j), lp.lower(j)));
               lp.changeLower(j, newbound);
            }
         }
      }
   }
}



SPxSimplifier::Result SPxMainSM::removeEmpty(SPxLP& lp)
{

   // This method removes empty rows and columns from the LP.

   int remRows = 0;
   int remCols = 0;

   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      const SVector& row = lp.rowVector(i);

      if (row.size() == 0)
      {
         MSG_DEBUG( (*spxout) << "IMAISM07 row " << i
                           << ": empty ->"; )

         if (LT(lp.rhs(i), 0.0, feastol()) || GT(lp.lhs(i), 0.0, feastol()))
         {
            MSG_DEBUG( (*spxout) << " infeasible lhs=" << lp.lhs(i)
                              << " rhs=" << lp.rhs(i) << std::endl; )
            return INFEASIBLE;
         }
         MSG_DEBUG( (*spxout) << " removed" << std::endl; )

         EmptyConstraintPS* EmptyConstraintPSptr = 0;
         spx_alloc(EmptyConstraintPSptr);
         m_hist.append(new (EmptyConstraintPSptr) EmptyConstraintPS(lp, i));

         ++remRows;
         removeRow(lp, i);

         ++m_stat[EMPTY_ROW];
      }
   }

   for(int j = lp.nCols()-1; j >= 0; --j)
   {
      const SVector& col = lp.colVector(j);

      if (col.size() == 0)
      {
         MSG_DEBUG( (*spxout) << "IMAISM08 col " << j
                           << ": empty -> maxObj=" << lp.maxObj(j)
                           << " lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j); )

         Real val;

         if (GT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.upper(j) >= infinity)
            {
               MSG_DEBUG( (*spxout) << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.upper(j);
         }
         else if (LT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.lower(j) <= -infinity)
            {
               MSG_DEBUG( (*spxout) << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.lower(j);
         }
         else
         {
            ASSERT_WARN( "WMAISM09", isZero(lp.maxObj(j), epsZero()) );
            // any value within the bounds is ok
            if (lp.lower(j) > -infinity)
               val = lp.lower(j);
            else if (lp.upper(j) < infinity)
               val = lp.upper(j);
            else
               val = 0.0;
         }
         MSG_DEBUG( (*spxout) << " removed" << std::endl; )

         FixBoundsPS* FixBoundsPSptr = 0;
         FixVariablePS* FixVariablePSptr = 0;
         spx_alloc(FixBoundsPSptr);
         spx_alloc(FixVariablePSptr);
         m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, val));
         m_hist.append(new (FixVariablePSptr) FixVariablePS(lp, *this, j, val));

         ++remCols;
         removeCol(lp, j);

         ++m_stat[EMPTY_COL];
      }
   }

   if (remRows + remCols > 0)
   {
      m_remRows += remRows;
      m_remCols += remCols;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (empty rows/colums) removed "
                        << remRows << " rows, "
                        << remCols << " cols"
                        << std::endl; )

   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::removeRowSingleton(SPxLP& lp, const SVector& row, int& i)
{
   assert(row.size() == 1);

   Real aij = row.value(0);
   int  j   = row.index(0);
   Real lo  = -infinity;
   Real up  =  infinity;

   MSG_DEBUG( (*spxout) << "IMAISM22 row " << i
                     << ": singleton -> val=" << aij
                     << " lhs=" << lp.lhs(i)
                     << " rhs=" << lp.rhs(i); )

   if (GT(aij, 0.0, epsZero()))           // aij > 0
   {
      lo = (lp.lhs(i) <= -infinity) ? -infinity : lp.lhs(i) / aij;
      up = (lp.rhs(i) >=  infinity) ?  infinity : lp.rhs(i) / aij;
   }
   else if (LT(aij, 0.0, epsZero()))      // aij < 0
   {
      lo = (lp.rhs(i) >=  infinity) ? -infinity : lp.rhs(i) / aij;
      up = (lp.lhs(i) <= -infinity) ?  infinity : lp.lhs(i) / aij;
   }
   else if (LT(lp.rhs(i), 0.0, feastol()) || GT(lp.lhs(i), 0.0, feastol()))
   {
      // aij == 0, rhs < 0 or lhs > 0
      MSG_DEBUG( (*spxout) << " infeasible" << std::endl; )
      return INFEASIBLE;
   }

   if (isZero(lo, epsZero()))
      lo = 0.0;

   if (isZero(up, epsZero()))
      up = 0.0;

   MSG_DEBUG( (*spxout) << " removed, lower=" << lo
                     << " (" << lp.lower(j)
                     << ") upper=" << up
                     << " (" << lp.upper(j)
                     << ")" << std::endl; )

   bool stricterUp = false;
   bool stricterLo = false;

   Real oldLo = lp.lower(j);
   Real oldUp = lp.upper(j);

   if (LTrel(up, lp.upper(j), feastol()))
   {
      lp.changeUpper(j, up);
      stricterUp = true;
   }
   if (GTrel(lo, lp.lower(j), feastol()))
   {
      lp.changeLower(j, lo);
      stricterLo = true;
   }

   RowSingletonPS* RowSingletonPSptr = 0;
   spx_alloc(RowSingletonPSptr);
   m_hist.append(new (RowSingletonPSptr) RowSingletonPS(lp, i, j, stricterLo, stricterUp, lp.lower(j), lp.upper(j), oldLo, oldUp));

   removeRow(lp, i);

   m_remRows++;
   m_remNzos++;
   ++m_stat[SINGLETON_ROW];

   return OKAY;
}

SPxSimplifier::Result SPxMainSM::simplifyRows(SPxLP& lp, bool& again)
{

   // This method simplifies the rows of the LP.
   //
   // The following operations are done:
   // 1. detect implied free variables
   // 2. detect implied free constraints
   // 3. detect infeasible constraints
   // 4. remove unconstrained constraints
   // 5. remove empty constraints
   // 6. remove row singletons and tighten the corresponding variable bounds if necessary
   // 7. detect forcing rows and fix the corresponding variables

   int remRows = 0;
   int remNzos = 0;
   int chgLRhs = 0;
   int chgBnds = 0;
   int keptBnds = 0;
   int keptLRhs = 0;

   int oldRows = lp.nRows();

   bool redundantLower;
   bool redundantUpper;
   bool redundantLhs;
   bool redundantRhs;

   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      const SVector& row = lp.rowVector(i);


      // compute bounds on constraint value
      Real lhsBnd = 0.0; // minimal activity (finite summands)
      Real rhsBnd = 0.0; // maximal activity (finite summands)
      int  lhsCnt = 0; // number of -infinity summands in minimal activity
      int  rhsCnt = 0; // number of +infinity summands in maximal activity

      for(int k = 0; k < row.size(); ++k)
      {
         Real aij = row.value(k);
         int  j   = row.index(k);

         if( !isNotZero(aij, 1.0 / infinity) )
         {
            MSG_WARNING( (*spxout), (*spxout) << "Warning: tiny nonzero coefficient " << aij << " in row " << i << "\n" );
         }

         if (aij > 0.0)
         {
            if (lp.lower(j) <= -infinity)
               ++lhsCnt;
            else
               lhsBnd += aij * lp.lower(j);

            if (lp.upper(j) >= infinity)
               ++rhsCnt;
            else
               rhsBnd += aij * lp.upper(j);
         }
         else if (aij < 0.0)
         {
            if (lp.lower(j) <= -infinity)
               ++rhsCnt;
            else
               rhsBnd += aij * lp.lower(j);

            if (lp.upper(j) >= infinity)
               ++lhsCnt;
            else
               lhsBnd += aij * lp.upper(j);
         }
      }

#if FREE_BOUNDS
      // 1. detect implied free variables
      if (rhsCnt <= 1 || lhsCnt <= 1)
      {
         for(int k = 0; k < row.size(); ++k)
         {
            Real aij = row.value(k);
            int  j   = row.index(k);

            redundantLower = false;
            redundantUpper = false;

            ASSERT_WARN( "WMAISM12", isNotZero(aij, 1.0 / infinity) );

            if (aij > 0.0)
            {
               if (lp.lhs(i) > -infinity && lp.lower(j) > -infinity && rhsCnt <= 1 && NErel(lp.lhs(i), rhsBnd, feastol())
                  // do not perform if strongly different orders of magnitude occur
                  && spxAbs(lp.lhs(i) / maxAbs(rhsBnd, 1.0)) > Param::epsilon())
               {
                  Real lo    = -infinity;
                  Real scale = maxAbs(lp.lhs(i), rhsBnd);

                  if (scale < 1.0)
                     scale = 1.0;

                  Real z = (lp.lhs(i) / scale) - (rhsBnd / scale);

                  if (isZero(z, epsZero()))
                     z = 0.0;

                  assert(rhsCnt > 0 || lp.upper(j) < infinity);

                  if (rhsCnt == 0)
                     lo = lp.upper(j) + z * scale / aij;
                  else if (lp.upper(j) >= infinity)
                     lo = z * scale / aij;

                  if (isZero(lo, epsZero()))
                     lo = 0.0;

                  if (GErel(lo, lp.lower(j), feastol()))
                  {
                     MSG_DEBUG( (*spxout) << "IMAISM13 row " << i
                                       << ": redundant lower bound on x" << j
                                       << " -> lower=" << lo
                                       << " (" << lp.lower(j)
                                       << ")" << std::endl; )

                     redundantLower = true;
                  }

               }
               if (lp.rhs(i) < infinity && lp.upper(j) < infinity && lhsCnt <= 1 && NErel(lp.rhs(i), lhsBnd, feastol())
                  // do not perform if strongly different orders of magnitude occur
                  && spxAbs(lp.rhs(i) / maxAbs(lhsBnd, 1.0)) > Param::epsilon())
               {
                  Real up    = infinity;
                  Real scale = maxAbs(lp.rhs(i), lhsBnd);

                  if (scale < 1.0)
                     scale = 1.0;

                  Real z = (lp.rhs(i) / scale) - (lhsBnd / scale);

                  if (isZero(z, epsZero()))
                     z = 0.0;

                  assert(lhsCnt > 0 || lp.lower(j) > -infinity);

                  if (lhsCnt == 0)
                     up = lp.lower(j) + z * scale / aij;
                  else if (lp.lower(j) <= -infinity)
                     up = z * scale / aij;

                  if (isZero(up, epsZero()))
                     up = 0.0;

                  if (LErel(up, lp.upper(j), feastol()))
                  {
                     MSG_DEBUG( (*spxout) << "IMAISM14 row " << i
                                       << ": redundant upper bound on x" << j
                                       << " -> upper=" << up
                                       << " (" << lp.upper(j)
                                       << ")" << std::endl; )

                     redundantUpper = true;
                  }
               }
               if (redundantLower)
               {
                  // no upper bound on x_j OR redundant upper bound
                  if ((lp.upper(j) >= infinity) || redundantUpper || (!m_keepbounds))
                  {
                     ++lhsCnt;
                     lhsBnd -= aij * lp.lower(j);

                     lp.changeLower(j, -infinity);
                     ++chgBnds;
                  }
                  else
                     ++keptBnds;
               }
               if (redundantUpper)
               {
                  // no lower bound on x_j OR redundant lower bound
                  if ((lp.lower(j) <= -infinity) || redundantLower || (!m_keepbounds))
                  {
                     ++rhsCnt;
                     rhsBnd -= aij * lp.upper(j);

                     lp.changeUpper(j, infinity);
                     ++chgBnds;
                  }
                  else
                     ++keptBnds;
               }
            }
            else if (aij < 0.0)
            {
               if (lp.lhs(i) > -infinity && lp.upper(j) < infinity && rhsCnt <= 1 && NErel(lp.lhs(i), rhsBnd, feastol())
                  // do not perform if strongly different orders of magnitude occur
                  && spxAbs(lp.lhs(i) / maxAbs(rhsBnd, 1.0)) > Param::epsilon())
               {
                  Real up    = infinity;
                  Real scale = maxAbs(lp.lhs(i), rhsBnd);

                  if (scale < 1.0)
                     scale = 1.0;

                  Real z = (lp.lhs(i) / scale) - (rhsBnd / scale);

                  if (isZero(z, epsZero()))
                     z = 0.0;

                  assert(rhsCnt > 0 || lp.lower(j) > -infinity);

                  if (rhsCnt == 0)
                     up = lp.lower(j) + z * scale / aij;
                  else if (lp.lower(j) <= -infinity)
                     up = z * scale / aij;

                  if (isZero(up, epsZero()))
                     up = 0.0;

                  if (LErel(up, lp.upper(j), feastol()))
                  {
                     MSG_DEBUG( (*spxout) << "IMAISM15 row " << i
                                       << ": redundant upper bound on x" << j
                                       << " -> upper=" << up
                                       << " (" << lp.upper(j)
                                       << ")" << std::endl; )

                     redundantUpper = true;
                  }
               }
               if (lp.rhs(i) < infinity && lp.lower(j) > -infinity && lhsCnt <= 1 && NErel(lp.rhs(i), lhsBnd, feastol())
                  // do not perform if strongly different orders of magnitude occur
                  && spxAbs(lp.rhs(i) / maxAbs(lhsBnd, 1.0)) > Param::epsilon())
               {
                  Real lo    = -infinity;
                  Real scale = maxAbs(lp.rhs(i), lhsBnd);

                  if (scale < 1.0)
                     scale = 1.0;

                  Real z = (lp.rhs(i) / scale) - (lhsBnd / scale);

                  if (isZero(z, epsZero()))
                     z = 0.0;

                  assert(lhsCnt > 0 || lp.upper(j) < infinity);

                  if (lhsCnt == 0)
                     lo = lp.upper(j) + z * scale / aij;
                  else if (lp.upper(j) >= infinity)
                     lo = z * scale / aij;

                  if (isZero(lo, epsZero()))
                     lo = 0.0;

                  if (GErel(lo, lp.lower(j)))
                  {
                     MSG_DEBUG( (*spxout) << "IMAISM16 row " << i
                                       << ": redundant lower bound on x" << j
                                       << " -> lower=" << lo
                                       << " (" << lp.lower(j)
                                       << ")" << std::endl; )

                     redundantLower = true;
                  }
               }
               if (redundantUpper)
               {
                  // no lower bound on x_j OR redundant lower bound
                  if ((lp.lower(j) <= -infinity) || redundantLower || (!m_keepbounds))
                  {
                     ++lhsCnt;
                     lhsBnd -= aij * lp.upper(j);

                     lp.changeUpper(j, infinity);
                     ++chgBnds;
                  }
                  else
                     ++keptBnds;
               }
               if (redundantLower)
               {
                  // no upper bound on x_j OR redundant upper bound
                  if ((lp.upper(j) >= infinity) || redundantUpper || (!m_keepbounds))
                  {
                     ++rhsCnt;
                     rhsBnd -= aij * lp.lower(j);

                     lp.changeLower(j, -infinity);
                     ++chgBnds;
                  }
                  else
                     ++keptBnds;
               }
            }
         }
      }
#endif

#if FREE_LHS_RHS

      redundantLhs = false;
      redundantRhs = false;

      // 2. detect implied free constraints
      if (lp.lhs(i) > -infinity && lhsCnt == 0 && GErel(lhsBnd, lp.lhs(i), feastol()))
      {
         MSG_DEBUG( (*spxout) << "IMAISM17 row " << i
                           << ": redundant lhs -> lhsBnd=" << lhsBnd
                           << " lhs=" << lp.lhs(i)
                           << std::endl; )

         redundantLhs = true;
      }
      if (lp.rhs(i) <  infinity && rhsCnt == 0 && LErel(rhsBnd, lp.rhs(i), feastol()))
      {
         MSG_DEBUG( (*spxout) << "IMAISM18 row " << i
                           << ": redundant rhs -> rhsBnd=" << rhsBnd
                           << " rhs=" << lp.rhs(i)
                           << std::endl; )

         redundantRhs = true;
      }
      if (redundantLhs)
      {
         // no rhs for constraint i OR redundant rhs
         if ((lp.rhs(i) >= infinity) || redundantRhs || (!m_keepbounds))
         {
            lp.changeLhs(i, -infinity);
            ++chgLRhs;
         }
         else
            ++keptLRhs;
      }
      if (redundantRhs)
      {
         // no lhs for constraint i OR redundant lhs
         if ((lp.lhs(i) <= -infinity) || redundantLhs || (!m_keepbounds))
         {
            lp.changeRhs(i, infinity);
            ++chgLRhs;
         }
         else
            ++keptLRhs;
      }
#endif

      // 3. infeasible constraint
      if (LTrel(lp.rhs(i), lp.lhs(i), feastol())                 ||
          (LTrel(rhsBnd,   lp.lhs(i), feastol()) && rhsCnt == 0) ||
          (GTrel(lhsBnd,   lp.rhs(i), feastol()) && lhsCnt == 0))
      {
         MSG_DEBUG( (*spxout) << "IMAISM19 row " << std::setprecision(20) << i
                           << ": infeasible -> lhs=" << lp.lhs(i)
                           << " rhs=" << lp.rhs(i)
                           << " lhsBnd=" << lhsBnd
                           << " rhsBnd=" << rhsBnd
                           << std::endl; )
         return INFEASIBLE;
      }

#if FREE_CONSTRAINT
      // 4. unconstrained constraint
      if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
      {
         MSG_DEBUG( (*spxout) << "IMAISM20 row " << i
                           << ": unconstrained -> removed" << std::endl; )

         FreeConstraintPS* FreeConstraintPSptr = 0;
         spx_alloc(FreeConstraintPSptr);
         m_hist.append(new (FreeConstraintPSptr) FreeConstraintPS(lp, i));

         ++remRows;
         remNzos += row.size();
         removeRow(lp, i);

         ++m_stat[FREE_ROW];

         continue;
      }
#endif

#if EMPTY_CONSTRAINT
      // 5. empty constraint
      if (row.size() == 0)
      {
         MSG_DEBUG( (*spxout) << "IMAISM21 row " << i
                           << ": empty ->"; )

         if (LT(lp.rhs(i), 0.0, feastol()) || GT(lp.lhs(i), 0.0, feastol()))
         {
            MSG_DEBUG( (*spxout) << " infeasible lhs=" << lp.lhs(i)
                              << " rhs=" << lp.rhs(i) << std::endl; )
            return INFEASIBLE;
         }
         MSG_DEBUG( (*spxout) << " removed" << std::endl; )

         EmptyConstraintPS* EmptyConstraintPSptr = 0;
         spx_alloc(EmptyConstraintPSptr);
         m_hist.append(new (EmptyConstraintPSptr) EmptyConstraintPS(lp, i));

         ++remRows;
         removeRow(lp, i);

         ++m_stat[EMPTY_ROW];

         continue;
      }
#endif

#if ROW_SINGLETON
      // 6. row singleton
      if (row.size() == 1)
      {
         removeRowSingleton(lp, row, i);
         continue;
      }
#endif

#if FORCE_CONSTRAINT
      // 7. forcing constraint (postsolving)
      // fix variables to obtain the upper bound on constraint value
      if (rhsCnt == 0 && EQrel(rhsBnd, lp.lhs(i), feastol()))
      {
         MSG_DEBUG( (*spxout) << "IMAISM24 row " << i
                           << ": forcing constraint fix on lhs ->"
                           << " lhs=" << lp.lhs(i)
                           << " rhsBnd=" << rhsBnd
                           << std::endl; )

         DataArray<bool> fixedCol(row.size());
         DataArray<Real> lowers(row.size());
         DataArray<Real> uppers(row.size());

         for(int k = 0; k < row.size(); ++k)
         {
            Real aij = row.value(k);
            int  j   = row.index(k);

            fixedCol[k] = !(EQrel(lp.upper(j), lp.lower(j), m_epsilon));

            lowers[k] = lp.lower(j);
            uppers[k] = lp.upper(j);

            ASSERT_WARN( "WMAISM25", isNotZero(aij, 1.0 /infinity) );

            if (aij > 0.0)
               lp.changeLower(j, lp.upper(j));
            else
               lp.changeUpper(j, lp.lower(j));
         }

         ForceConstraintPS* ForceConstraintPSptr = 0;
         spx_alloc(ForceConstraintPSptr);
         m_hist.append(new (ForceConstraintPSptr) ForceConstraintPS(lp, i, true, fixedCol, lowers, uppers));

         ++remRows;
         remNzos += row.size();
         removeRow(lp, i);

         ++m_stat[FORCE_ROW];

         continue;
      }
      // fix variables to obtain the lower bound on constraint value
      if (lhsCnt == 0 && EQrel(lhsBnd, lp.rhs(i), feastol()))
      {
         MSG_DEBUG( (*spxout) << "IMAISM26 row " << i
                           << ": forcing constraint fix on rhs ->"
                           << " rhs=" << lp.rhs(i)
                           << " lhsBnd=" << lhsBnd
                           << std::endl; )

         DataArray<bool> fixedCol(row.size());
         DataArray<Real> lowers(row.size());
         DataArray<Real> uppers(row.size());

         for(int k = 0; k < row.size(); ++k)
         {
            Real aij   = row.value(k);
            int  j     = row.index(k);

            fixedCol[k] = !(EQrel(lp.upper(j), lp.lower(j), m_epsilon));

            lowers[k] = lp.lower(j);
            uppers[k] = lp.upper(j);

            ASSERT_WARN( "WMAISM27", isNotZero(aij, 1.0 / infinity) );

            if (aij > 0.0)
               lp.changeUpper(j, lp.lower(j));
            else
               lp.changeLower(j, lp.upper(j));
         }

         ForceConstraintPS* ForceConstraintPSptr = 0;
         spx_alloc(ForceConstraintPSptr);
         m_hist.append(new (ForceConstraintPSptr) ForceConstraintPS(lp, i, false, fixedCol, lowers, uppers));

         ++remRows;
         remNzos += row.size();
         removeRow(lp, i);

         ++m_stat[FORCE_ROW];

         continue;
      }
#endif
   }

   assert(remRows > 0 || remNzos == 0);

   if (remRows + chgLRhs + chgBnds > 0)
   {
      m_remRows += remRows;
      m_remNzos += remNzos;
      m_chgLRhs += chgLRhs;
      m_chgBnds += chgBnds;
      m_keptBnds += keptBnds;
      m_keptLRhs += keptLRhs;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (rows) removed "
                        << remRows << " rows, "
                        << remNzos << " non-zeros, "
                        << chgBnds << " col bounds, "
                        << chgLRhs << " row bounds; kept "
                        << keptBnds << " column bounds, "
                        << keptLRhs << " row bounds"
                        << std::endl; )
      if( remRows > m_minReduction * oldRows )
         again = true;
   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::simplifyCols(SPxLP& lp, bool& again)
{

   // This method simplifies the columns of the LP.
   //
   // The following operations are done:
   // 1. detect empty columns and fix corresponding variables
   // 2. detect variables that are unconstrained from below or above
   //    and fix corresponding variables or remove involved constraints
   // 3. fix variables
   // 4. use column singleton variables with zero objective to adjust constraint bounds
   // 5. free column singleton combined with doubleton equation are
   //    used to make the column singleton variable free
   // 6. substitute (implied) free column singletons

   int remRows = 0;
   int remCols = 0;
   int remNzos = 0;
   int chgBnds = 0;

   int oldCols = lp.nCols();
   int oldRows = lp.nRows();

   for(int j = lp.nCols()-1; j >= 0; --j)
   {
       const SVector& col = lp.colVector(j);

      // infeasible bounds
      if (GTrel(lp.lower(j), lp.upper(j), feastol()))
      {
         MSG_DEBUG( (*spxout) << "IMAISM29 col " << j
                           << ": infeasible bounds on x" << j
                           << " -> lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j)
                           << std::endl; )
         return INFEASIBLE;
      }

      // 1. empty column
      if (col.size() == 0)
      {
#if EMPTY_COLUMN
	 MSG_DEBUG( (*spxout) << "IMAISM30 col " << j
                           << ": empty -> maxObj=" << lp.maxObj(j)
                           << " lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j); )

         Real val;

         if (GT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.upper(j) >= infinity)
            {
               MSG_DEBUG( (*spxout) << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.upper(j);
         }
         else if (LT(lp.maxObj(j), 0.0, epsZero()))
         {
            if (lp.lower(j) <= -infinity)
            {
               MSG_DEBUG( (*spxout) << " unbounded" << std::endl; )
               return UNBOUNDED;
            }
            val = lp.lower(j);
         }
         else
         {
            assert(isZero(lp.maxObj(j), epsZero()));
            // any value within the bounds is ok
            if (lp.lower(j) > -infinity)
               val = lp.lower(j);
            else if (lp.upper(j) < infinity)
               val = lp.upper(j);
            else
               val = 0.0;
         }
         MSG_DEBUG( (*spxout) << " removed" << std::endl; )

         FixBoundsPS* FixBoundsPSptr = 0;
         FixVariablePS* FixVariablePSptr = 0;
         spx_alloc(FixBoundsPSptr);
         spx_alloc(FixVariablePSptr);
         m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, val));
         m_hist.append(new (FixVariablePSptr) FixVariablePS(lp, *this, j, val));

         ++remCols;
         removeCol(lp, j);

         ++m_stat[EMPTY_COL];

         continue;
#endif
      }

      if (NErel(lp.lower(j), lp.upper(j), feastol()))
      {
         // will be set to false if any constraint implies a bound on the variable
         bool loFree = true;
         bool upFree = true;

         // 1. fix and remove variables
         for(int k = 0; k < col.size(); ++k)
         {
            if (!loFree && !upFree)
               break;

            int i = col.index(k);

            // warn since this unhandled case may slip through unnoticed otherwise
            ASSERT_WARN( "WMAISM31", isNotZero(col.value(k), 1.0 / infinity) );

            if (col.value(k) > 0.0)
            {
               if (lp.rhs(i) <  infinity)
                  upFree = false;

               if (lp.lhs(i) > -infinity)
                  loFree = false;
            }
            else if (col.value(k) < 0.0)
            {
               if (lp.rhs(i) <  infinity)
                  loFree = false;

               if (lp.lhs(i) > -infinity)
                  upFree = false;
            }
         }

         // 2. detect variables that are unconstrained from below or above
         // max  3 x
         // s.t. 5 x >= 8
         if (GT(lp.maxObj(j), 0.0, epsZero()) && upFree)
         {
#if FIX_VARIABLE
            MSG_DEBUG( (*spxout) << "IMAISM32 col " << j
                              << ": x" << j
                              << " unconstrained above ->"; )

            if (lp.upper(j) >= infinity)
            {
               MSG_DEBUG( (*spxout) << " unbounded" << std::endl; )

               return UNBOUNDED;
            }
            MSG_DEBUG( (*spxout) << " fixed at upper=" << lp.upper(j) << std::endl; )

            FixBoundsPS* FixBoundsPSptr = 0;
            spx_alloc(FixBoundsPSptr);
            m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, lp.upper(j)));
            lp.changeLower(j, lp.upper(j));
         }
         // max -3 x
         // s.t. 5 x <= 8
         else if (LT(lp.maxObj(j), 0.0, epsZero()) && loFree)
         {
            MSG_DEBUG( (*spxout) << "IMAISM33 col " << j
                              << ": x" << j
                              << " unconstrained below ->"; )

            if (lp.lower(j) <= -infinity)
            {
               MSG_DEBUG( (*spxout) << " unbounded" << std::endl; )

               return UNBOUNDED;
            }
            MSG_DEBUG( (*spxout) << " fixed at lower=" << lp.lower(j) << std::endl; )

            FixBoundsPS* FixBoundsPSptr = 0;
            spx_alloc(FixBoundsPSptr);
            m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, lp.lower(j)));
            lp.changeUpper(j, lp.lower(j));
#endif
         }
         else if (isZero(lp.maxObj(j), epsZero()))
         {
#if FREE_ZERO_OBJ_VARIABLE
            bool unconstrained_below = loFree && lp.lower(j) <= -infinity;
            bool unconstrained_above = upFree && lp.upper(j) >= infinity;

            if (unconstrained_below || unconstrained_above)
            {
               MSG_DEBUG( (*spxout) << "IMAISM34 col " << j
                                 << ": x" << j
                                 << " unconstrained "
                                 << (unconstrained_below ? "below" : "above")
                                 << " with zero objective (" << lp.maxObj(j)
                                 << ")" << std::endl; )

               SVector col_idx_sorted(col);

               // sort col elements by increasing idx
               IdxCompare compare;
               SPxQuicksort(col_idx_sorted.mem(), col_idx_sorted.size(), compare);

               FreeZeroObjVariablePS* FreeZeroObjVariablePSptr = 0;
               spx_alloc(FreeZeroObjVariablePSptr);
               m_hist.append(new (FreeZeroObjVariablePSptr) FreeZeroObjVariablePS(lp, j, unconstrained_below, col_idx_sorted));

               // we have to remove the rows with larger idx first, because otherwise the rows are reorder and indices
               // are out-of-date
               remRows += col.size();
               for(int k = col_idx_sorted.size()-1; k >= 0; --k)
                  removeRow(lp, col_idx_sorted.index(k));

               // remove column
               removeCol(lp, j);

               // statistics
               for(int k = 0; k < col.size(); ++k)
               {
                  int l   =  col.index(k);
                  remNzos += lp.rowVector(l).size();
               }

               ++m_stat[FREE_ZOBJ_COL];
               ++remCols;

               continue;
            }
#endif
         }
      }

#if FIX_VARIABLE
      // 3. fix variable
      if (EQrel(lp.lower(j), lp.upper(j), feastol()))
      {
         MSG_DEBUG( (*spxout) << "IMAISM36 col " << j
                           << ": x" << j
                           << " fixed -> lower=" << lp.lower(j)
                           << " upper=" << lp.upper(j) << std::endl; )

         fixColumn(lp, j);

         ++remCols;
         remNzos += col.size();
         removeCol(lp, j);

         ++m_stat[FIX_COL];

         continue;
      }
#endif

      // handle column singletons
      if (col.size() == 1)
      {
         Real aij = col.value(0);
         int  i   = col.index(0);

         // 4. column singleton with zero objective
         if (isZero(lp.maxObj(j), epsZero()))
         {
#if ZERO_OBJ_COL_SINGLETON
            MSG_DEBUG( (*spxout) << "IMAISM37 col " << j
                              << ": singleton in row " << i
                              << " with zero objective"; )

            Real lhs = -infinity;
            Real rhs = +infinity;

            if (GT(aij, 0.0, epsZero()))
            {
               if (lp.lhs(i) > -infinity && lp.upper(j) <  infinity)
                  lhs = lp.lhs(i) - aij * lp.upper(j);
               if (lp.rhs(i) <  infinity && lp.lower(j) > -infinity)
                  rhs = lp.rhs(i) - aij * lp.lower(j);
            }
            else if (LT(aij, 0.0, epsZero()))
            {
               if (lp.lhs(i) > -infinity && lp.lower(j) > -infinity)
                  lhs = lp.lhs(i) - aij * lp.lower(j);
               if (lp.rhs(i) <  infinity && lp.upper(j) <  infinity)
                  rhs = lp.rhs(i) - aij * lp.upper(j);
            }
            else
            {
               lhs = lp.lhs(i);
               rhs = lp.rhs(i);
            }

            if (isZero(lhs, epsZero()))
               lhs = 0.0;
            if (isZero(rhs, epsZero()))
               rhs = 0.0;

            MSG_DEBUG( (*spxout) << " removed -> lhs=" << lhs
                              << " (" << lp.lhs(i)
                              << ") rhs=" << rhs
                              << " (" << lp.rhs(i)
                              << ")" << std::endl; )

            ZeroObjColSingletonPS* ZeroObjColSingletonPSptr = 0;
            spx_alloc(ZeroObjColSingletonPSptr);
            m_hist.append(new (ZeroObjColSingletonPSptr) ZeroObjColSingletonPS(lp, *this, j, i));

            lp.changeRange(i, lhs, rhs);

            ++remCols;
            ++remNzos;
            removeCol(lp, j);

            ++m_stat[ZOBJ_SINGLETON_COL];

            if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
            {
               FreeConstraintPS* FreeConstraintPSptr = 0;
               spx_alloc(FreeConstraintPSptr);
               m_hist.append(new (FreeConstraintPSptr) FreeConstraintPS(lp, i));

               ++remRows;
               removeRow(lp, i);

               ++m_stat[FREE_ROW];
            }

            continue;
#endif
         }

         // 5. not free column singleton combined with doubleton equation
         else if (EQrel(lp.lhs(i), lp.rhs(i), feastol())             &&
                  lp.rowVector(i).size() == 2                         &&
                  (lp.lower(j) > -infinity || lp.upper(j) < infinity))
         {
#if DOUBLETON_EQUATION
            MSG_DEBUG( (*spxout) << "IMAISM38 col " << j
                              << ": singleton in row " << i
                              << " with doubleton equation ->"; )

            Real lhs = lp.lhs(i);

            const SVector& row = lp.rowVector(i);

            Real aik;
            int  k;

            if (row.index(0) == j)
            {
               aik = row.value(1);
               k   = row.index(1);
            }
            else if (row.index(1) == j)
            {
               aik = row.value(0);
               k   = row.index(0);
            }
            else
               throw SPxInternalCodeException("XMAISM11 This should never happen.");

            ASSERT_WARN( "WMAISM39", isNotZero(aik, 1.0 / infinity) );

            Real lo, up;
            Real oldLower = lp.lower(k);
            Real oldUpper = lp.upper(k);

            Real scale1 = maxAbs(lhs, aij * lp.upper(j));
            Real scale2 = maxAbs(lhs, aij * lp.lower(j));

            if (scale1 < 1.0)
               scale1 = 1.0;
            if (scale2 < 1.0)
               scale2 = 1.0;

            Real z1 = (lhs / scale1) - (aij * lp.upper(j) / scale1);
            Real z2 = (lhs / scale2) - (aij * lp.lower(j) / scale2);

            if (isZero(z1, epsZero()))
               z1 = 0.0;
            if (isZero(z2, epsZero()))
               z2 = 0.0;

            if (GT(aij * aik, 0.0, epsZero()))
            {
               lo = (lp.upper(j) >=  infinity) ? -infinity : z1 * scale1 / aik;
               up = (lp.lower(j) <= -infinity) ?  infinity : z2 * scale2 / aik;
            }
            else if (LT(aij * aik, 0.0, epsZero()))
            {
               lo = (lp.lower(j) <= -infinity) ? -infinity : z2 * scale2 / aik;
               up = (lp.upper(j) >=  infinity) ?  infinity : z1 * scale1 / aik;
            }
            else
               throw SPxInternalCodeException("XMAISM12 This should never happen.");

            if (GTrel(lo, lp.lower(k), epsZero()))
               lp.changeLower(k, lo);

            if (LTrel(up, lp.upper(k), epsZero()))
               lp.changeUpper(k, up);

            MSG_DEBUG( (*spxout) << " made free, bounds on x" << k
                              << ": lower=" << lp.lower(k)
                              << " (" << oldLower
                              << ") upper=" << lp.upper(k)
                              << " (" << oldUpper
                              << ")" << std::endl; )

            // infeasible bounds
            if (GTrel(lp.lower(k), lp.upper(k), feastol()))
            {
               MSG_DEBUG( (*spxout) << "new bounds are infeasible"
                                 << std::endl; )
               return INFEASIBLE;
            }

            DoubletonEquationPS* DoubletonEquationPSptr = 0;
            spx_alloc(DoubletonEquationPSptr);
            m_hist.append(new (DoubletonEquationPSptr) DoubletonEquationPS(lp, j, k, i, oldLower, oldUpper));

            if (lp.lower(j) > -infinity && lp.upper(j) < infinity)
               chgBnds += 2;
            else
               ++chgBnds;

            lp.changeBounds(j, -infinity, infinity);

            ++m_stat[DOUBLETON_ROW];
#endif
         }

         // 6. (implied) free column singleton
         if (lp.lower(j) <= -infinity && lp.upper(j) >= infinity)
         {
#if FREE_COL_SINGLETON
            Real slackVal = lp.lhs(i);

            // constraint i is an inequality constraint -> transform into equation type
            if (NErel(lp.lhs(i), lp.rhs(i), feastol()))
            {
               MSG_DEBUG( (*spxout) << "IMAISM40 col " << j
                                 << ": free singleton in inequality constraint" << std::endl; )

               // do nothing if constraint i is unconstrained
               if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
                  continue;

               // introduce slack variable to obtain equality constraint
               Real sMaxObj = lp.maxObj(j) / aij; // after substituting variable j in objective
               Real sLo     = lp.lhs(i);
               Real sUp     = lp.rhs(i);

               if (GT(sMaxObj, 0.0, epsZero()))
               {
                  if (sUp >= infinity)
                  {
                     MSG_DEBUG( (*spxout) << " -> problem unbounded" << std::endl; )
                     return UNBOUNDED;
                  }
                  slackVal = sUp;
               }
               else if (LT(sMaxObj, 0.0, epsZero()))
               {
                  if (sLo <= -infinity)
                  {
                     MSG_DEBUG( (*spxout) << " -> problem unbounded" << std::endl; )
                     return UNBOUNDED;
                  }
                  slackVal = sLo;
               }
               else
               {
                  assert(isZero(sMaxObj, epsZero()));
                  // any value within the bounds is ok
                  if (sLo > -infinity)
                     slackVal = sLo;
                  else if (sUp < infinity)
                     slackVal = sUp;
                  else
                     throw SPxInternalCodeException("XMAISM13 This should never happen.");
               }
            }

            FreeColSingletonPS* FreeColSingletonPSptr = 0;
            spx_alloc(FreeColSingletonPSptr);
            m_hist.append(new (FreeColSingletonPSptr) FreeColSingletonPS(lp, *this, j, i, slackVal));

            MSG_DEBUG( (*spxout) << "IMAISM41 col " << j
                              << ": free singleton removed" << std::endl; )

            const SVector& row = lp.rowVector(i);

            for (int h = 0; h < row.size(); ++h)
            {
               int k = row.index(h);

               if (k != j)
               {
                  Real new_obj = lp.obj(k) - (lp.obj(j) * row.value(h) / aij);
                  lp.changeObj(k, new_obj);
               }
            }

            ++remRows;
            ++remCols;
            remNzos += row.size();
            removeRow(lp, i);
            removeCol(lp, j);

            ++m_stat[FREE_SINGLETON_COL];

            continue;
#endif
         }
      }
   }

   if (remCols + remRows > 0)
   {
      m_remRows += remRows;
      m_remCols += remCols;
      m_remNzos += remNzos;
      m_chgBnds += chgBnds;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (columns) removed "
                        << remRows << " rows, "
                        << remCols << " cols, "
                        << remNzos << " non-zeros, "
                        << chgBnds << " col bounds"
                        << std::endl; )
      if( remCols + remRows > m_minReduction * (oldCols + oldRows) )
         again = true;
   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::simplifyDual(SPxLP& lp, bool& again)
{

   // This method simplifies LP using the following dual structures:
   //
   // 1. dominated columns
   // 2. weakly dominated columns
   //
   // For constructing the dual variables, it is assumed that the objective sense is max

   int remRows = 0;
   int remCols = 0;
   int remNzos = 0;

   int oldRows = lp.nRows();
   int oldCols = lp.nCols();

   DataArray<bool> colSingleton(lp.nCols());
   DVector         dualVarLo(lp.nRows());
   DVector         dualVarUp(lp.nRows());
   DVector         dualConsLo(lp.nCols());
   DVector         dualConsUp(lp.nCols());

   // init
   for(int i = lp.nRows()-1; i >= 0; --i)
   {
      // check for unconstrained constraints
      if (lp.lhs(i) <= -infinity && lp.rhs(i) >= infinity)
      {
         MSG_DEBUG( (*spxout) << "IMAISM43 row " << i
                           << ": unconstrained" << std::endl; )

         FreeConstraintPS* FreeConstraintPSptr = 0;
         spx_alloc(FreeConstraintPSptr);
         m_hist.append(new (FreeConstraintPSptr) FreeConstraintPS(lp, i));

         ++remRows;
         remNzos += lp.rowVector(i).size();
         removeRow(lp, i);

         ++m_stat[FREE_ROW];

         continue;
      }

      // corresponds to maximization sense
      dualVarLo[i] = (lp.lhs(i) <= -infinity) ? 0.0 : -infinity;
      dualVarUp[i] = (lp.rhs(i) >=  infinity) ? 0.0 :  infinity;
   }

   // compute bounds on the dual variables using column singletons
   for(int j = 0; j < lp.nCols(); ++j)
   {
      if (lp.colVector(j).size() == 1)
      {
         int  i   = lp.colVector(j).index(0);
         Real aij = lp.colVector(j).value(0);

         ASSERT_WARN( "WMAISM44", isNotZero(aij, 1.0 / infinity) );

         Real bound = lp.maxObj(j) / aij;

         if (aij > 0)
         {
            if (lp.lower(j) <= -infinity && bound < dualVarUp[i])
               dualVarUp[i] = bound;
            if (lp.upper(j) >=  infinity && bound > dualVarLo[i])
               dualVarLo[i] = bound;
         }
         else if (aij < 0)
         {
            if (lp.lower(j) <= -infinity && bound > dualVarLo[i])
               dualVarLo[i] = bound;
            if (lp.upper(j) >=  infinity && bound < dualVarUp[i])
               dualVarUp[i] = bound;
         }
      }

   }

   // compute bounds on the dual constraints
   for(int j = 0; j < lp.nCols(); ++j)
   {
      dualConsLo[j] = dualConsUp[j] = 0.0;

      const SVector& col = lp.colVector(j);

      for(int k = 0; k < col.size(); ++k)
      {
         if (dualConsLo[j] <= -infinity && dualConsUp[j] >= infinity)
            break;

         Real aij = col.value(k);
         int  i   = col.index(k);

         ASSERT_WARN( "WMAISM45", isNotZero(aij, 1.0 / infinity) );

         if (aij > 0)
         {
            if (dualVarLo[i] <= -infinity)
               dualConsLo[j] = -infinity;
            else
               dualConsLo[j] += aij * dualVarLo[i];

            if (dualVarUp[i] >= infinity)
               dualConsUp[j] = infinity;
            else
               dualConsUp[j] += aij * dualVarUp[i];
         }
         else if (aij < 0)
         {
            if (dualVarLo[i] <= -infinity)
               dualConsUp[j] = infinity;
            else
               dualConsUp[j] += aij * dualVarLo[i];

            if (dualVarUp[i] >= infinity)
               dualConsLo[j] = -infinity;
            else
               dualConsLo[j] += aij * dualVarUp[i];
         }
      }
   }

   for(int j = lp.nCols()-1; j >= 0; --j)
   {
      if (lp.colVector(j).size() <= 1)
         continue;

      // dual infeasibility checks
      if (LTrel(dualConsUp[j], dualConsLo[j], opttol()))
      {
         MSG_DEBUG( (*spxout) << "IMAISM46 col " << j
                           << ": dual infeasible -> dual lhs bound=" << dualConsLo[j]
                           << " dual rhs bound=" << dualConsUp[j] << std::endl; )
         return DUAL_INFEASIBLE;
      }

      Real obj = lp.maxObj(j);

      // 1. dominated column
      // Is the problem really unbounded in the cases below ??? Or is only dual infeasibility be shown
      if (GTrel(obj, dualConsUp[j], opttol()))
      {
#if DOMINATED_COLUMN
         MSG_DEBUG( (*spxout) << "IMAISM47 col " << j
                           << ": dominated -> maxObj=" << obj
                           << " dual rhs bound=" << dualConsUp[j] << std::endl; )

         if (lp.upper(j) >= infinity)
         {
            MSG_INFO2( (*spxout), (*spxout) << " unbounded" << std::endl; )
            return UNBOUNDED;
         }

         MSG_DEBUG( (*spxout) << " fixed at upper=" << lp.upper(j) << std::endl; )

         FixBoundsPS* FixBoundsPSptr = 0;
         spx_alloc(FixBoundsPSptr);
         m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, lp.upper(j)));
         lp.changeLower(j, lp.upper(j));

         ++m_stat[DOMINATED_COL];
#endif
      }
      else if (LTrel(obj, dualConsLo[j], opttol()))
      {
#if DOMINATED_COLUMN
         MSG_DEBUG( (*spxout) << "IMAISM48 col " << j
                           << ": dominated -> maxObj=" << obj
                           << " dual lhs bound=" << dualConsLo[j] << std::endl; )

         if (lp.lower(j) <= -infinity)
         {
            MSG_INFO2( (*spxout), (*spxout) << " unbounded" << std::endl; )
            return UNBOUNDED;
         }

         MSG_DEBUG( (*spxout) << " fixed at lower=" << lp.lower(j) << std::endl; )

         FixBoundsPS* FixBoundsPSptr = 0;
         spx_alloc(FixBoundsPSptr);
         m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, lp.lower(j)));
         lp.changeUpper(j, lp.lower(j));

         ++m_stat[DOMINATED_COL];
#endif
      }

      // 2. weakly dominated column (no postsolving)
      else if (lp.upper(j) < infinity && EQrel(obj, dualConsUp[j], opttol()))
      {
#if WEAKLY_DOMINATED_COLUMN
         MSG_DEBUG( (*spxout) << "IMAISM49 col " << j
                           << ": weakly dominated -> maxObj=" << obj
                           << " dual rhs bound=" << dualConsUp[j] << std::endl; )

         FixBoundsPS* FixBoundsPSptr = 0;
         spx_alloc(FixBoundsPSptr);
         m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, lp.upper(j)));
         lp.changeLower(j, lp.upper(j));

         ++m_stat[WEAKLY_DOMINATED_COL];
#endif
      }
      else if (lp.lower(j) > -infinity && EQrel(obj, dualConsLo[j], opttol()))
      {
#if WEAKLY_DOMINATED_COLUMN
         MSG_DEBUG( (*spxout) << "IMAISM50 col " << j
                           << ": weakly dominated -> maxObj=" << obj
                           << " dual lhs bound=" << dualConsLo[j] << std::endl; )

         FixBoundsPS* FixBoundsPSptr = 0;
         spx_alloc(FixBoundsPSptr);
         m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j, lp.lower(j)));
         lp.changeUpper(j, lp.lower(j));

         ++m_stat[WEAKLY_DOMINATED_COL];
#endif
      }

      // fix column
      if (EQrel(lp.lower(j), lp.upper(j), feastol()))
      {
#if FIX_VARIABLE
         fixColumn(lp, j);

         ++remCols;
         remNzos += lp.colVector(j).size();
         removeCol(lp, j);

         ++m_stat[FIX_COL];
#endif
      }
   }


   assert(remRows > 0 || remCols > 0 || remNzos == 0);

   if (remCols + remRows > 0)
   {
      m_remRows += remRows;
      m_remCols += remCols;
      m_remNzos += remNzos;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (dual) removed "
                        << remRows << " rows, "
                        << remCols << " cols, "
                        << remNzos << " non-zeros"
                        << std::endl; )
      if( remCols + remRows > m_minReduction * (oldCols + oldRows) )
         again = true;
   }
   return OKAY;
}



SPxSimplifier::Result SPxMainSM::multiaggregation(SPxLP& lp, bool& again)
{
   // this simplifier eliminates rows and columns by performing multi aggregations as identified by the constraint
   // activities.
   int remRows = 0;
   int remCols = 0;
   int remNzos = 0;

   int oldRows = lp.nRows();
   int oldCols = lp.nCols();

   DVector upLocks(lp.nCols());
   DVector downLocks(lp.nCols());

   for(int j = lp.nCols()-1; j >= 0; --j)
   {
      // setting the locks on the variables
      upLocks[j] = 0;
      downLocks[j] = 0;

      if (lp.colVector(j).size() <= 1)
         continue;

      const SVector& col = lp.colVector(j);
      for(int k = 0; k < col.size(); ++k)
      {
         Real aij = col.value(k);

         ASSERT_WARN( "WMAISM45", isNotZero(aij, 1.0 / infinity) );

         if(GT(lp.lhs(col.index(k)), -infinity) && LT(lp.rhs(col.index(k)), infinity))
         {
            upLocks[j]++;
            downLocks[j]++;
         }
         else if (GT(lp.lhs(col.index(k)), -infinity))
         {
            if(aij > 0)
               downLocks[j]++;
            else if (aij < 0)
               upLocks[j]++;
         }
         else if (LT(lp.rhs(col.index(k)), infinity))
         {
            if(aij > 0)
               upLocks[j]++;
            else if (aij < 0)
               downLocks[j]++;
         }
      }

      // multi-aggregate column
      if (upLocks[j] == 1 || downLocks[j] == 1)
      {
         Real lower = lp.lower(j);
         Real upper = lp.upper(j);
         int maxOtherLocks;
         int bestpos = -1;
         bool bestislhs = true;



         for(int k = 0; k < col.size(); ++k)
         {
            int rowNumber;
            Real lhs;
            Real rhs;
            bool lhsExists;
            bool rhsExists;
            bool aggLhs;
            bool aggRhs;

            Real val = col.value(k);

            rowNumber = col.index(k);
            lhs = lp.lhs(rowNumber);
            rhs = lp.rhs(rowNumber);

            if( EQ(lhs, rhs, feastol()) )
               continue;

            lhsExists = GT(lhs, -infinity);
            rhsExists = LT(rhs, infinity);

            if (lp.rowVector(rowNumber).size() <= 2)
               maxOtherLocks = INT_MAX;
            else if (lp.rowVector(rowNumber).size() == 3)
               maxOtherLocks = 3;
            else if (lp.rowVector(rowNumber).size() == 4)
               maxOtherLocks = 2;
            else
               maxOtherLocks = 1;

            aggLhs = lhsExists
               && ((col.value(k) > 0.0 && lp.maxObj(j) <= 0.0 && downLocks[j] == 1 && upLocks[j] <= maxOtherLocks)
               || (col.value(k) < 0.0 && lp.maxObj(j) >= 0.0 && upLocks[j] == 1 && downLocks[j] <= maxOtherLocks));
            aggRhs = rhsExists
               && ((col.value(k) > 0.0 && lp.maxObj(j) >= 0.0 && upLocks[j] == 1 && downLocks[j] <= maxOtherLocks)
               || (col.value(k) < 0.0 && lp.maxObj(j) <= 0.0 && downLocks[j] == 1 && upLocks[j] <= maxOtherLocks));

            if (aggLhs || aggRhs)
            {
               Real minRes = 0;   // this is the minimum value that the aggregation can attain
               Real maxRes = 0;   // this is the maximum value that the aggregation can attain

               // computing the minimum and maximum residuals if variable j is set to zero.
               computeMinMaxResidualActivity(lp, rowNumber, j, minRes, maxRes);

               // we will try to aggregate to the lhs
               if (aggLhs)
               {
                  Real minVal;
                  Real maxVal;

                  // computing the values of the upper and lower bounds for the aggregated variables
                  computeMinMaxValues(lp, lhs, val, minRes, maxRes, minVal, maxVal);

                  assert(LE(minVal, maxVal));

                  // if the bounds of the aggregation and the original variable are equivalent, then we can reduce
                  if ((minVal > -infinity && GT(minVal, lower, feastol()))
                     && (maxVal < infinity && LT(maxVal, upper, feastol())))
                  {
                     bestpos = col.index(k);
                     bestislhs = true;
                     break;
                  }
               }

               // we will try to aggregate to the rhs
               if (aggRhs)
               {
                  Real minVal;
                  Real maxVal;

                  // computing the values of the upper and lower bounds for the aggregated variables
                  computeMinMaxValues(lp, rhs, val, minRes, maxRes, minVal, maxVal);

                  assert(LE(minVal, maxVal));

                  if ((minVal > -infinity && GT(minVal, lower, feastol()))
                     && (maxVal < infinity && LT(maxVal, upper, feastol())))
                  {
                     bestpos = col.index(k);
                     bestislhs = false;
                     break;
                  }
               }
            }
         }

         // it is only possible to aggregate if a best position has been found
         if( bestpos >= 0 )
         {
            const SVector& bestRow = lp.rowVector(bestpos);
            // aggregating the variable and applying the fixings to the all other constraints
            Real aggConstant = (bestislhs ? lp.lhs(bestpos) : lp.rhs(bestpos));   // this is the lhs or rhs of the aggregated row
            Real aggAij = bestRow[j];                                   // this is the coefficient of the deleted col

            MSG_DEBUG(
               (*spxout) << "IMAISM51 col " << j
                                            << ": Aggregating row: " << bestpos
                                            << " Aggregation Constant=" << aggConstant
                                            << " Coefficient of aggregated col=" << aggAij << std::endl;
               )

            MultiAggregationPS* MultiAggregationPSptr = 0;
            spx_alloc(MultiAggregationPSptr);
            m_hist.append(new (MultiAggregationPSptr) MultiAggregationPS(lp, *this, bestpos, j, aggConstant));

            for(int k = 0; k < col.size(); ++k)
            {
               if(col.index(k) != bestpos)
               {
                  int rowNumber = col.index(k);
                  DVector updateRow(lp.nCols());
                  Real updateRhs = lp.rhs(col.index(k));
                  Real updateLhs = lp.lhs(col.index(k));

                  updateRow = lp.rowVector(col.index(k));

                  // updating the row with the best row
                  for(int l = 0; l < bestRow.size(); l++)
                  {
                     if(bestRow.index(l) != j)
                     {
                        if(lp.rowVector(rowNumber).pos(bestRow.index(l)) >= 0)
                           lp.changeElement(rowNumber, bestRow.index(l), updateRow[bestRow.index(l)]
                              - updateRow[j]*bestRow.value(l)/aggAij);
                        else
                           lp.changeElement(rowNumber, bestRow.index(l), -1.0*updateRow[j]*bestRow.value(l)/aggAij);
                     }
                  }

                  // NOTE: I don't know whether we should change the LHS and RHS if they are currently at infinity
                  if(LT(lp.rhs(rowNumber), infinity))
                     lp.changeRhs(rowNumber, updateRhs - updateRow[j]*aggConstant/aggAij);
                  if(GT(lp.lhs(rowNumber), -infinity))
                     lp.changeLhs(rowNumber, updateLhs - updateRow[j]*aggConstant/aggAij);

                  assert(LE(lp.lhs(rowNumber), lp.rhs(rowNumber)));
               }
            }

            for(int l = 0; l < bestRow.size(); l++)
            {
               if(bestRow.index(l) != j)
                  lp.changeMaxObj(bestRow.index(l), lp.maxObj(bestRow.index(l)) - lp.maxObj(j)*bestRow.value(l)/aggAij);
            }

            ++remCols;
            remNzos += lp.colVector(j).size();
            removeCol(lp, j);
            ++remRows;
            remNzos += lp.rowVector(bestpos).size();
            removeRow(lp, bestpos);

            ++m_stat[MULTI_AGG];
         }
      }
   }


   assert(remRows > 0 || remCols > 0 || remNzos == 0);

   if (remCols + remRows > 0)
   {
      m_remRows += remRows;
      m_remCols += remCols;
      m_remNzos += remNzos;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (multi-aggregation) removed "
                        << remRows << " rows, "
                        << remCols << " cols, "
                        << remNzos << " non-zeros"
                        << std::endl; )
      if( remCols + remRows > m_minReduction * (oldCols + oldRows) )
         again = true;
   }
   return OKAY;
}



SPxSimplifier::Result SPxMainSM::duplicateRows(SPxLP& lp, bool& again)
{

   // This method simplifies the LP by removing duplicate rows
   // Duplicates are detected using the algorithm of Bixby and Wagner [1987]

   // Possible extension: use generalized definition of duplicate rows according to Andersen and Andersen
   // However: the resulting sparsification is often very small since the involved rows are usually very sparse

   int remRows = 0;
   int remNzos = 0;

   int oldRows = lp.nRows();

   // remove empty rows and columns
   SPxSimplifier::Result ret = removeEmpty(lp);
   if (ret != OKAY)
      return ret;

#if ROW_SINGLETON
   int rs_remRows = 0;
   for (int i = 0; i < lp.nRows(); ++i)
   {
      const SVector& row = lp.rowVector(i);

      if (row.size() == 1)
      {
         removeRowSingleton(lp, row, i);
         rs_remRows++;
      }
   }

   if (rs_remRows > 0)
   {
      MSG_INFO2( (*spxout), (*spxout) << "Simplifier duplicate rows (row singleton stage) removed "
                        << rs_remRows << " rows, "
                        << rs_remRows << " non-zeros"
                        << std::endl; )
   }
#endif

   if (lp.nRows() < 2)
      return OKAY;

   DataArray<int>    pClass(lp.nRows());           // class of parallel rows
   DataArray<int>    classSize(lp.nRows());        // size of each class
   DataArray<Real>   scale(lp.nRows());            // scaling factor for each row
   int*              idxMem = 0;

   try
   {
      spx_alloc(idxMem, lp.nRows());
   }
   catch( const SPxMemoryException& x )
   {
      spx_free(idxMem);
      throw x;
   }

   IdxSet idxSet(lp.nRows(), idxMem);           // set of feasible indices for new pClass

   // init
   pClass[0]    = 0;
   scale[0]     = 0.0;
   classSize[0] = lp.nRows();

   for(int i = 1; i < lp.nRows(); ++i)
   {
      pClass[i] = 0;
      scale[i]  = 0.0;
      classSize[i] = 0;
      idxSet.addIdx(i);
   }

   Real oldVal = 0.0;

   // main loop
   for(int j = 0; j < lp.nCols(); ++j)
   {
      const SVector& col = lp.colVector(j);

      for(int k = 0; k < col.size(); ++k)
      {
         Real aij = col.value(k);
         int  i   = col.index(k);

         if (scale[i] == 0.0)
            scale[i] = aij;

         m_classSetRows[pClass[i]].add(i, aij / scale[i]);
         if (--classSize[pClass[i]] == 0)
            idxSet.addIdx(pClass[i]);
      }

      // update each parallel class with non-zero column entry
      for(int m = 0; m < col.size(); ++m)
      {
         int k = pClass[col.index(m)];

         if (m_classSetRows[k].size() > 0)
         {
            // sort classSet[k] w.r.t. scaled column values
            ElementCompare compare;

            if (m_classSetRows[k].size() > 1)
               SPxQuicksort(m_classSetRows[k].mem(), m_classSetRows[k].size(), compare);

            // use new index first
            int classIdx = idxSet.index(0);
            idxSet.remove(0);

            for(int l = 0; l < m_classSetRows[k].size(); ++l)
            {
               if (l != 0 && NErel(m_classSetRows[k].value(l), oldVal, epsZero()))
               {
                  classIdx = idxSet.index(0);
                  idxSet.remove(0);
               }

               pClass[m_classSetRows[k].index(l)] = classIdx;
               ++classSize[classIdx];

               oldVal = m_classSetRows[k].value(l);
            }

            m_classSetRows[k].clear();
         }
      }
   }

   spx_free(idxMem);

   DataArray<bool> remRow(lp.nRows());

   for(int k = 0; k < lp.nRows(); ++k )
      m_dupRows[k].clear();

   for(int k = 0; k < lp.nRows(); ++k)
   {
      remRow[k] = false;
      m_dupRows[pClass[k]].add(k, 0.0);
   }

   const int nRowsOld_tmp = lp.nRows();
   int* perm_tmp = 0;
   spx_alloc(perm_tmp, nRowsOld_tmp);

   for(int j = 0; j < nRowsOld_tmp; ++j)
   {
      perm_tmp[j] = 0;
   }

   int idxFirstDupRows = -1;
   int idxLastDupRows = -1;
   int numDelRows = 0;

   for(int k = 0; k < lp.nRows(); ++k)
   {
      if (m_dupRows[k].size() > 1 && !(lp.rowVector(m_dupRows[k].index(0)).size() == 1))
      {
         idxLastDupRows = k;

         if(idxFirstDupRows < 0)
         {
            idxFirstDupRows = k;
         }

         for(int l = 1; l < m_dupRows[k].size(); ++l)
         {
            int i = m_dupRows[k].index(l);
            perm_tmp[i] = -1;
         }

         numDelRows += (m_dupRows[k].size()-1);
      }
   }

   {
      int k_tmp, j_tmp = -1;
      for (k_tmp = j_tmp = 0; k_tmp < nRowsOld_tmp; ++k_tmp)
      {
         if (perm_tmp[k_tmp] >= 0)
            perm_tmp[k_tmp] = j_tmp++;
      }
   }

   // store rhs and lhs changes for combined update
   bool doChangeRanges = false;
   DVector newLhsVec(lp.lhs());
   DVector newRhsVec(lp.rhs());

   for(int k = 0; k < lp.nRows(); ++k)
   {
      if (m_dupRows[k].size() > 1 && !(lp.rowVector(m_dupRows[k].index(0)).size() == 1))
      {
         MSG_DEBUG( (*spxout) << "IMAISM53 " << m_dupRows[k].size()
                           << " duplicate rows found" << std::endl; )

         m_stat[DUPLICATE_ROW] += m_dupRows[k].size()-1;

         // index of one non-column singleton row in dupRows[k]
         int  rowIdx    = -1;
         int  maxLhsIdx = -1;
         int  minRhsIdx = -1;
         Real maxLhs    = -infinity;
         Real minRhs    = +infinity;

         DataArray<bool> isLhsEqualRhs(m_dupRows[k].size());

         // determine strictest bounds on constraint
         for(int l = 0; l < m_dupRows[k].size(); ++l)
         {
            int i = m_dupRows[k].index(l);
            isLhsEqualRhs[l] = (lp.lhs(i) == lp.rhs(i));

            ASSERT_WARN( "WMAISM54", isNotZero(scale[i], 1.0 / infinity) );

            if (rowIdx == -1)
            {
               rowIdx = i;
               maxLhs = lp.lhs(rowIdx);
               minRhs = lp.rhs(rowIdx);
            }
            else
            {
               Real scaledLhs, scaledRhs;
               Real factor = scale[rowIdx] / scale[i];

               if (factor > 0)
               {
                  scaledLhs = (lp.lhs(i) <= -infinity) ? -infinity : lp.lhs(i) * factor;
                  scaledRhs = (lp.rhs(i) >=  infinity) ?  infinity : lp.rhs(i) * factor;
               }
               else
               {
                  scaledLhs = (lp.rhs(i) >=  infinity) ? -infinity : lp.rhs(i) * factor;
                  scaledRhs = (lp.lhs(i) <= -infinity) ?  infinity : lp.lhs(i) * factor;
               }
               if (scaledLhs > maxLhs)
               {
                  maxLhs    = scaledLhs;
                  maxLhsIdx = i;
               }
               if (scaledRhs < minRhs)
               {
                  minRhs    = scaledRhs;
                  minRhsIdx = i;
               }

               remRow[i] = true;
            }
         }

         if (rowIdx != -1)
         {
            Real newLhs = (maxLhs > lp.lhs(rowIdx)) ? maxLhs : lp.lhs(rowIdx);
            Real newRhs = (minRhs < lp.rhs(rowIdx)) ? minRhs : lp.rhs(rowIdx);

            if(k == idxLastDupRows)
            {
               DataArray<int> da_perm(nRowsOld_tmp);
               for(int j = 0; j < nRowsOld_tmp; ++j)
               {
                  da_perm[j] = perm_tmp[j];
               }
               DuplicateRowsPS* DuplicateRowsPSptr = 0;
               spx_alloc(DuplicateRowsPSptr);
               m_hist.append(new (DuplicateRowsPSptr) DuplicateRowsPS(lp, rowIdx, maxLhsIdx, minRhsIdx, m_dupRows[k], scale, da_perm, isLhsEqualRhs, true, EQrel(newLhs, newRhs), k==idxFirstDupRows));
            }
            else
            {
               DataArray<int> da_perm_empty(0);
               DuplicateRowsPS* DuplicateRowsPSptr = 0;
               spx_alloc(DuplicateRowsPSptr);
               m_hist.append(new (DuplicateRowsPSptr) DuplicateRowsPS(lp, rowIdx, maxLhsIdx, minRhsIdx, m_dupRows[k], scale, da_perm_empty, isLhsEqualRhs, false, EQrel(newLhs, newRhs), k == idxFirstDupRows));
            }

            if (maxLhs > lp.lhs(rowIdx) || minRhs < lp.rhs(rowIdx))
            {
               // modify lhs and rhs of constraint rowIdx
               doChangeRanges = true;

               if (LTrel(newRhs, newLhs, feastol()))
               {
                  MSG_DEBUG( (*spxout) << "IMAISM55 duplicate rows yield infeasible bounds:"
                                    << " lhs=" << newLhs
                                    << " rhs=" << newRhs << std::endl; )
                  spx_free(perm_tmp);
                  return INFEASIBLE;
               }
               // if we accept the infeasibility we should clean up the values to avoid problems later
               if( newRhs < newLhs )
                  newRhs = newLhs;

               newLhsVec[rowIdx] = newLhs;
               newRhsVec[rowIdx] = newRhs;
            }
         }
      }
   }

   // change ranges for all modified constraints by one single call (more efficient)
   if (doChangeRanges)
   {
      lp.changeRange(newLhsVec, newRhsVec);
   }

   // remove all rows by one single method call (more efficient)
   const int nRowsOld = lp.nRows();
   int* perm = 0;
   spx_alloc(perm, nRowsOld);

   for(int i = 0; i < nRowsOld; ++i)
   {
      if (remRow[i])
      {
         perm[i] = -1;
         ++remRows;
         remNzos += lp.rowVector(i).size();
      }
      else
         perm[i] = 0;
   }
   lp.removeRows(perm);

   for(int i = 0; i < nRowsOld; ++i)
   {
      // assert that the pre-computed permutation was correct
      assert(perm[i] == perm_tmp[i]);

      // update the global index mapping
      if (perm[i] >= 0)
         m_rIdx[perm[i]] = m_rIdx[i];
   }

   spx_free(perm);
   spx_free(perm_tmp);

   if (remRows + remNzos > 0)
   {
      m_remRows += remRows;
      m_remNzos += remNzos;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (duplicate rows) removed "
                        << remRows << " rows, "
                        << remNzos << " non-zeros"
                        << std::endl; )
      if( remRows > m_minReduction * oldRows )
         again = true;

   }
   return OKAY;
}

SPxSimplifier::Result SPxMainSM::duplicateCols(SPxLP& lp, bool& again)
{

   // This method simplifies the LP by removing duplicate columns
   // Duplicates are detected using the algorithm of Bixby and Wagner [1987]

   int remCols = 0;
   int remNzos = 0;

   // remove empty rows and columns
   SPxSimplifier::Result ret = removeEmpty(lp);
   if (ret != OKAY)
      return ret;

   if (lp.nCols() < 2)
      return OKAY;

   DataArray<int>    pClass(lp.nCols());          // class of parallel columns
   DataArray<int>    classSize(lp.nCols());       // size of each class
   DataArray<Real>   scale(lp.nCols());           // scaling factor for each column
   int*              idxMem = 0;

   try
   {
      spx_alloc(idxMem, lp.nCols());
   }
   catch( const SPxMemoryException& x )
   {
      spx_free(idxMem);
      throw x;
   }

   IdxSet idxSet(lp.nCols(), idxMem);  // set of feasible indices for new pClass

   // init
   pClass[0]    = 0;
   scale[0]     = 0.0;
   classSize[0] = lp.nCols();

   for(int j = 1; j < lp.nCols(); ++j)
   {
      pClass[j] = 0;
      scale[j]  = 0.0;
      classSize[j] = 0;
      idxSet.addIdx(j);
   }

   Real oldVal = 0.0;

   // main loop
   for(int i = 0; i < lp.nRows(); ++i)
   {
      const SVector& row = lp.rowVector(i);

      for(int k = 0; k < row.size(); ++k)
      {
         Real aij = row.value(k);
         int  j   = row.index(k);

         if (scale[j] == 0.0)
            scale[j] = aij;

         m_classSetCols[pClass[j]].add(j, aij / scale[j]);
         if (--classSize[pClass[j]] == 0)
            idxSet.addIdx(pClass[j]);
      }

      // update each parallel class with non-zero row entry
      for(int m = 0; m < row.size(); ++m)
      {
         int k = pClass[row.index(m)];

         if (m_classSetCols[k].size() > 0)
         {
            // sort classSet[k] w.r.t. scaled row values
            ElementCompare compare;

            if (m_classSetCols[k].size() > 1)
               SPxQuicksort(m_classSetCols[k].mem(), m_classSetCols[k].size(), compare);

            // use new index first
            int classIdx = idxSet.index(0);
            idxSet.remove(0);

            for(int l = 0; l < m_classSetCols[k].size(); ++l)
            {
               if (l != 0 && NErel(m_classSetCols[k].value(l), oldVal, epsZero()))
               {
                  // start new parallel class
                  classIdx = idxSet.index(0);
                  idxSet.remove(0);
               }

               pClass[m_classSetCols[k].index(l)] = classIdx;
               ++classSize[classIdx];

               oldVal = m_classSetCols[k].value(l);
            }

            m_classSetCols[k].clear();
         }
      }
   }

   spx_free(idxMem);

   DataArray<bool> remCol(lp.nCols());
   DataArray<bool> fixAndRemCol(lp.nCols());

   for( int k = 0; k < lp.nCols(); ++k )
      m_dupCols[k].clear();

   for(int k = 0; k < lp.nCols(); ++k)
   {
      remCol[k] = false;
      fixAndRemCol[k] = false;
      m_dupCols[pClass[k]].add(k, 0.0);
   }

   bool hasDuplicateCol = false;
   DataArray<int>  m_perm_empty(0);

   for(int k = 0; k < lp.nCols(); ++k)
   {
      if (m_dupCols[k].size() > 1 && !(lp.colVector(m_dupCols[k].index(0)).size() == 1))
      {
         MSG_DEBUG( (*spxout) << "IMAISM58 " << m_dupCols[k].size()
                           << " duplicate columns found" << std::endl; )

         if (!hasDuplicateCol)
         {
            DuplicateColsPS* DuplicateColsPSptr = 0;
            spx_alloc(DuplicateColsPSptr);
            m_hist.append(new (DuplicateColsPSptr) DuplicateColsPS(lp, 0, 0, 1.0, m_perm_empty, true));
            hasDuplicateCol = true;
         }

         for(int l = 0; l < m_dupCols[k].size(); ++l)
         {
            for(int m = 0; m < m_dupCols[k].size(); ++m)
            {
               int j1  = m_dupCols[k].index(l);
               int j2  = m_dupCols[k].index(m);

               if (l != m && !remCol[j1] && !remCol[j2])
               {
                  Real cj1 = lp.maxObj(j1);
                  Real cj2 = lp.maxObj(j2);

                  // A.j1 = factor * A.j2
                  Real factor = scale[j1] / scale[j2];
                  Real objDif = cj1 - cj2 * scale[j1] / scale[j2];

                  ASSERT_WARN( "WMAISM59", isNotZero(factor, epsZero()) );

                  if (isZero(objDif, epsZero()))
                  {
                     // case 1: objectives also duplicate

                     // if 0 is not within the column bounds, we are not able to postsolve if the aggregated column has
                     // status ZERO, hence we skip this case
                     if (LErel(lp.lower(j1), 0.0) && GErel(lp.upper(j1), 0.0)
                        && LErel(lp.lower(j2), 0.0) && GErel(lp.upper(j2), 0.0))
                     {
                        DuplicateColsPS* DuplicateColsPSptr = 0;
                        spx_alloc(DuplicateColsPSptr);
                        // variable substitution xj2' := xj2 + factor * xj1 <=> xj2 = -factor * xj1 + xj2'
                        m_hist.append(new (DuplicateColsPSptr) DuplicateColsPS(lp, j1, j2, factor, m_perm_empty));

                        // update bounds of remaining column j2 (new column j2')
                        if (factor > 0)
                        {
                           if (lp.lower(j2) <= -infinity || lp.lower(j1) <= -infinity)
                              lp.changeLower(j2, -infinity);
                           else
                              lp.changeLower(j2, lp.lower(j2) + factor * lp.lower(j1));

                           if (lp.upper(j2) >= infinity || lp.upper(j1) >= infinity)
                              lp.changeUpper(j2, infinity);
                           else
                              lp.changeUpper(j2, lp.upper(j2) + factor * lp.upper(j1));
                        }
                        else if (factor < 0)
                        {
                           if (lp.lower(j2) <= -infinity || lp.upper(j1) >= infinity)
                              lp.changeLower(j2, -infinity);
                           else
                              lp.changeLower(j2, lp.lower(j2) + factor * lp.upper(j1));

                           if (lp.upper(j2) >= infinity || lp.lower(j1) <= -infinity)
                              lp.changeUpper(j2, infinity);
                           else
                              lp.changeUpper(j2, lp.upper(j2) + factor * lp.lower(j1));
                        }

                        MSG_DEBUG( (*spxout) << "IMAISM60 two duplicate columns " << j1
                           << ", " << j2
                           << " replaced by one" << std::endl; )

                           remCol[j1] = true;

                        ++m_stat[SUB_DUPLICATE_COL];
                     }
                     else
                     {
                        MSG_DEBUG( (*spxout) << "IMAISM80 not removing two duplicate columns " << j1
                           << ", " << j2
                           << " because zero not contained in their bounds" << std::endl; )
                     }
                  }
                  else
                  {
                     // case 2: objectives not duplicate
                     // considered for maximization sense
                     if (lp.lower(j2) <= -infinity)
                     {
                        if (factor > 0 && objDif > 0)
                        {
                           if (lp.upper(j1) >= infinity)
                           {
                              MSG_DEBUG( (*spxout) << "IMAISM75 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }

                           // fix j1 at upper bound
                           MSG_DEBUG( (*spxout) << "IMAISM61 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at upper bound=" << lp.upper(j1) << std::endl; )

                           FixBoundsPS* FixBoundsPSptr = 0;
                           spx_alloc(FixBoundsPSptr);
                           m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j1, lp.upper(j1)));
                           lp.changeLower(j1, lp.upper(j1));
                        }
                        else if (factor < 0 && objDif < 0)
                        {
                           if (lp.lower(j1) <= -infinity)
                           {
                              MSG_DEBUG( (*spxout) << "IMAISM76 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }

                           // fix j1 at lower bound
                           MSG_DEBUG( (*spxout) << "IMAISM62 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at lower bound=" << lp.lower(j1) << std::endl; )

                           FixBoundsPS* FixBoundsPSptr = 0;
                           spx_alloc(FixBoundsPSptr);
                           m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j1, lp.lower(j1)));
                           lp.changeUpper(j1, lp.lower(j1));
                        }
                     }
                     else if (lp.upper(j2) >= infinity)
                     {
                        // fix j1 at upper bound
                        if (factor < 0 && objDif > 0)
                        {
                           if (lp.upper(j1) >= infinity)
                           {
                              MSG_DEBUG( (*spxout) << "IMAISM77 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }

                           // fix j1 at upper bound
                           MSG_DEBUG( (*spxout) << "IMAISM63 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at upper bound=" << lp.upper(j1) << std::endl; )

                           FixBoundsPS* FixBoundsPSptr = 0;
                           spx_alloc(FixBoundsPSptr);
                           m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j1, lp.upper(j1)));
                           lp.changeLower(j1, lp.upper(j1));
                        }

                        // fix j1 at lower bound
                        else if (factor > 0 && objDif < 0)
                        {
                           if (lp.lower(j1) <= -infinity)
                           {
                              MSG_DEBUG( (*spxout) << "IMAISM78 LP unbounded" << std::endl; )
                              return UNBOUNDED;
                           }

                           // fix j1 at lower bound
                           MSG_DEBUG( (*spxout) << "IMAISM64 two duplicate columns " << j1
                                             << ", " << j2
                                             << " first one fixed at lower bound=" << lp.lower(j1) << std::endl; )

                           FixBoundsPS* FixBoundsPSptr = 0;
                           spx_alloc(FixBoundsPSptr);
                           m_hist.append(new (FixBoundsPSptr) FixBoundsPS(lp, j1, lp.lower(j1)));
                           lp.changeUpper(j1, lp.lower(j1));
                        }
                     }
                     if (EQrel(lp.lower(j1), lp.upper(j1), feastol()))
                     {
                        remCol[j1] = true;
                        fixAndRemCol[j1] = true;

                        ++m_stat[FIX_DUPLICATE_COL];
                     }
                  }
               }
            }
         }
      }
   }

   for(int j = 0; j < lp.nCols(); ++j)
   {
      if(fixAndRemCol[j])
      {
         assert(remCol[j]);

         // correctIdx == false, because the index mapping will be handled by the postsolving in DuplicateColsPS
         fixColumn(lp, j, false);
      }
   }

   // remove all columns by one single method call (more efficient)
   const int nColsOld = lp.nCols();
   int* perm = 0;
   spx_alloc(perm, nColsOld);

   for(int j = 0; j < nColsOld; ++j)
   {
      if (remCol[j])
      {
         perm[j] = -1;
         ++remCols;
         remNzos += lp.colVector(j).size();
      }
      else
         perm[j] = 0;
   }
   lp.removeCols(perm);

   for(int j = 0; j < nColsOld; ++j)
   {
      if (perm[j] >= 0)
         m_cIdx[perm[j]] = m_cIdx[j];
   }

   DataArray<int> da_perm(nColsOld);
   for(int j = 0; j < nColsOld; ++j)
   {
       da_perm[j] = perm[j];
   }

   if (hasDuplicateCol)
   {
      DuplicateColsPS* DuplicateColsPSptr = 0;
      spx_alloc(DuplicateColsPSptr);
      m_hist.append(new (DuplicateColsPSptr) DuplicateColsPS(lp, 0, 0, 1.0, da_perm, false, true));
   }

   spx_free(perm);

   assert(remCols > 0 || remNzos == 0);

   if (remCols > 0)
   {
      m_remCols += remCols;
      m_remNzos += remNzos;

      MSG_INFO2( (*spxout), (*spxout) << "Simplifier (duplicate columns) removed "
                        << remCols << " cols, "
                        << remNzos << " non-zeros"
                        << std::endl; )
      if( remCols > m_minReduction * nColsOld )
         again = true;
   }
   return OKAY;
}

void SPxMainSM::fixColumn(SPxLP& lp, int j, bool correctIdx)
{

   assert(EQrel(lp.lower(j), lp.upper(j), feastol()));

   Real lo            = lp.lower(j);
   const SVector& col = lp.colVector(j);

   assert(NE(lo, infinity) && NE(lo, -infinity));

   MSG_DEBUG( (*spxout) << "IMAISM66 fix variable x" << j
                     << ": lower=" << lp.lower(j)
                     << " upper=" << lp.upper(j)
                     << std::endl; )

   if (isNotZero(lo, epsZero()))
   {
      for(int k = 0; k < col.size(); ++k)
      {
         int i = col.index(k);

         if (lp.rhs(i) < infinity)
         {
            Real y     = lo * col.value(k);
            Real scale = maxAbs(lp.rhs(i), y);

            if (scale < 1.0)
               scale = 1.0;

            Real rhs = (lp.rhs(i) / scale) - (y / scale);

            if (isZero(rhs, epsZero()))
               rhs = 0.0;
            else
               rhs *= scale;

            MSG_DEBUG( (*spxout) << "IMAISM67 row " << i
                              << ": rhs=" << rhs
                              << " (" << lp.rhs(i)
                              << ") aij=" << col.value(k)
                              << std::endl; )

            lp.changeRhs(i, rhs);
         }
         if (lp.lhs(i) > -infinity)
         {
            Real y     = lo * col.value(k);
            Real scale = maxAbs(lp.lhs(i), y);

            if (scale < 1.0)
               scale = 1.0;

            Real lhs = (lp.lhs(i) / scale) - (y / scale);

            if (isZero(lhs, epsZero()))
               lhs = 0.0;
            else
               lhs *= scale;

            MSG_DEBUG( (*spxout) << "IMAISM68 row " << i
                              << ": lhs=" << lhs
                              << " (" << lp.lhs(i)
                              << ") aij=" << col.value(k)
                              << std::endl; )

            lp.changeLhs(i, lhs);
         }
         assert(lp.lhs(i) <= lp.rhs(i));
      }
   }

   FixVariablePS* FixVariablePSptr = 0;
   spx_alloc(FixVariablePSptr);
   m_hist.append(new (FixVariablePSptr) FixVariablePS(lp, *this, j, lp.lower(j), correctIdx));
}

SPxSimplifier::Result SPxMainSM::simplify(SPxLP& lp, Real eps, Real ftol, Real otol, bool keepbounds)
{
   // transfer message handler
   spxout = lp.spxout;
   assert(spxout != 0);

   m_thesense = lp.spxSense();
   m_timeUsed->reset();
   m_timeUsed->start();

   m_objoffset = 0.0;
   m_cutoffbound = -infinity;
   m_pseudoobj = -infinity;

   m_remRows = 0;
   m_remCols = 0;
   m_remNzos = 0;
   m_chgBnds = 0;
   m_chgLRhs = 0;
   m_keptBnds = 0;
   m_keptLRhs = 0;

   m_result     = OKAY;
   bool   again = true;

   if(m_hist.size() > 0)
   {
      // delete pointers in old m_hist
      for(int k = 0; k < m_hist.size(); ++k)
      {
         m_hist[k]->~PostStep();
         spx_free(m_hist[k]);
      }
      m_hist.clear();
   }

   m_hist.reSize(0);
   m_postsolved = false;

   if( eps < 0.0 )
      throw SPxInterfaceException("XMAISM30 Cannot use negative epsilon in simplify().");

   if( ftol < 0.0 )
      throw SPxInterfaceException("XMAISM31 Cannot use negative feastol in simplify().");

   if( otol < 0.0 )
      throw SPxInterfaceException("XMAISM32 Cannot use negative opttol in simplify().");

   m_epsilon = eps;
   m_feastol = ftol;
   m_opttol = otol;


   MSG_INFO2( (*spxout),
      int numRangedRows = 0;
      int numBoxedCols = 0;

      for(int i = 0; i < lp.nRows(); ++i)
         if (lp.lhs(i) > -infinity && lp.rhs(i) < infinity)
            ++numRangedRows;
      for(int j = 0; j < lp.nCols(); ++j)
         if (lp.lower(j) > -infinity && lp.upper(j) < infinity)
            ++numBoxedCols;

      (*spxout) << "LP has "
                << numRangedRows << " ranged rows, "
                << numBoxedCols << " boxed columns"
                << std::endl;
   )

   m_stat.reSize(16);

   for(int k = 0; k < m_stat.size(); ++k)
      m_stat[k] = 0;

   m_addedcols = 0;
   handleRowObjectives(lp);

   m_prim.reDim(lp.nCols());
   m_slack.reDim(lp.nRows());
   m_dual.reDim(lp.nRows());
   m_redCost.reDim(lp.nCols());
   m_cBasisStat.reSize(lp.nCols());
   m_rBasisStat.reSize(lp.nRows());
   m_cIdx.reSize(lp.nCols());
   m_rIdx.reSize(lp.nRows());

   m_classSetRows.reSize(lp.nRows());
   m_classSetCols.reSize(lp.nCols());
   m_dupRows.reSize(lp.nRows());
   m_dupCols.reSize(lp.nCols());

   m_keepbounds = keepbounds;

   for(int i = 0; i < lp.nRows(); ++i)
      m_rIdx[i] = i;

   for(int j = 0; j < lp.nCols(); ++j)
      m_cIdx[j] = j;

   // round extreme values (set all values smaller than eps to zero and all values bigger than infinity/5 to infinity)
#if EXTREMES
   handleExtremes(lp);
#endif

   // main presolving loop
   while(again && m_result == OKAY)
   {
      again = false;

#if ROWS
      if (m_result == OKAY)
         m_result = simplifyRows(lp, again);
#endif

#if COLS
      if (m_result == OKAY)
         m_result = simplifyCols(lp, again);
#endif

#if DUAL
      if (m_result == OKAY)
         m_result = simplifyDual(lp, again);
#endif

#if DUPLICATE_ROWS
      if (m_result == OKAY)
         m_result = duplicateRows(lp, again);
#endif

#if DUPLICATE_COLS
      if (m_result == OKAY)
         m_result = duplicateCols(lp, again);
#endif

      if( !again )
      {
#if TRIVIAL_HEURISTICS
         trivialHeuristic(lp);
#endif

#if PSEUDOOBJ
         propagatePseudoobj(lp);
#endif

#if MULTI_AGGREGATE
      if (m_result == OKAY)
         m_result = multiaggregation(lp, again);
#endif
      }

   }

   // preprocessing detected infeasibility or unboundedness
   if (m_result != OKAY)
   {
      MSG_INFO1( (*spxout), (*spxout) << "Simplifier result: " << m_result << std::endl; )
      return m_result;
   }

   m_remCols -= m_addedcols;
   m_remNzos -= m_addedcols;
   MSG_INFO1( (*spxout), (*spxout) << "Simplifier removed "
                     << m_remRows << " rows, "
                     << m_remCols << " columns, "
                     << m_remNzos << " nonzeros, "
                     << m_chgBnds << " col bounds, "
                     << m_chgLRhs << " row bounds"
                     << std::endl; )

   if (keepbounds)
      MSG_INFO2( (*spxout), (*spxout) << "Simplifier kept "
                        << m_keptBnds << " column bounds, "
                        << m_keptLRhs << " row bounds"
                        << std::endl; )

   MSG_INFO1( (*spxout), (*spxout) << "Reduced LP has "
                     << lp.nRows() << " rows "
                     << lp.nCols() << " columns "
                     << lp.nNzos() << " nonzeros"
                     << std::endl; )

   MSG_INFO2( (*spxout),
      int numRangedRows = 0;
      int numBoxedCols  = 0;

      for(int i = 0; i < lp.nRows(); ++i)
         if (lp.lhs(i) > -infinity && lp.rhs(i) < infinity)
            ++numRangedRows;

      for(int j = 0; j < lp.nCols(); ++j)
         if (lp.lower(j) > -infinity && lp.upper(j) < infinity)
            ++numBoxedCols;

      (*spxout) << "Reduced LP has "
                << numRangedRows << " ranged rows, "
                << numBoxedCols << " boxed columns"
                << std::endl;
   )

   if (lp.nCols() == 0 && lp.nRows() == 0)
   {
      MSG_INFO1( (*spxout), (*spxout) << "Simplifier removed all rows and columns" << std::endl; )
      m_result = VANISHED;
   }

   MSG_INFO2( (*spxout), (*spxout) << "\nSimplifier performed:\n"
                     << m_stat[EMPTY_ROW]            << " empty rows\n"
                     << m_stat[FREE_ROW]             << " free rows\n"
                     << m_stat[SINGLETON_ROW]        << " singleton rows\n"
                     << m_stat[FORCE_ROW]            << " forcing rows\n"
                     << m_stat[EMPTY_COL]            << " empty columns\n"
                     << m_stat[FIX_COL]              << " fixed columns\n"
                     << m_stat[FREE_ZOBJ_COL]        << " free columns with zero objective\n"
                     << m_stat[ZOBJ_SINGLETON_COL]   << " singleton columns with zero objective\n"
                     << m_stat[DOUBLETON_ROW]        << " singleton columns combined with a doubleton equation\n"
                     << m_stat[FREE_SINGLETON_COL]   << " free singleton columns\n"
                     << m_stat[DOMINATED_COL]        << " dominated columns\n"
                     << m_stat[WEAKLY_DOMINATED_COL] << " weakly dominated columns\n"
                     << m_stat[DUPLICATE_ROW]        << " duplicate rows\n"
                     << m_stat[FIX_DUPLICATE_COL]    << " duplicate columns (fixed)\n"
                     << m_stat[SUB_DUPLICATE_COL]    << " duplicate columns (substituted)\n"
                     << m_stat[MULTI_AGG]            << " multi aggregation of variables\n"
                     << std::endl; );

   m_timeUsed->stop();

   return m_result;
}

void SPxMainSM::unsimplify(const Vector& x, const Vector& y, const Vector& s, const Vector& r,
                           const SPxSolver::VarStatus rows[], const SPxSolver::VarStatus cols[], bool isOptimal)
{
   MSG_INFO1( (*spxout), (*spxout) << " --- unsimplifying solution and basis" << std::endl; )
   assert(x.dim() <= m_prim.dim());
   assert(y.dim() <= m_dual.dim());
   assert(x.dim() == r.dim());
   assert(y.dim() == s.dim());

   // assign values of variables in reduced LP
   // NOTE: for maximization problems, we have to switch signs of dual and reduced cost values,
   // since simplifier assumes minimization problem
   for(int j = 0; j < x.dim(); ++j)
   {
      m_prim[j] = isZero(x[j], epsZero()) ? 0.0 : x[j];
      m_redCost[j] = isZero(r[j], epsZero()) ? 0.0 : (m_thesense == SPxLP::MAXIMIZE ? -r[j] : r[j]);
      m_cBasisStat[j] = cols[j];
   }
   for(int i = 0; i < y.dim(); ++i)
   {
      m_dual[i] = isZero(y[i], epsZero()) ? 0.0 : (m_thesense == SPxLP::MAXIMIZE ? -y[i] : y[i]);
      m_slack[i] = isZero(s[i], epsZero()) ? 0.0 : s[i];
      m_rBasisStat[i] = rows[i];
   }

   // undo preprocessing
   for(int k = m_hist.size()-1; k >= 0; --k)
   {
      MSG_DEBUG( std::cout << "unsimplifying " << m_hist[k]->getName() << "\n" );

      try
      {
         m_hist[k]->execute(m_prim, m_dual, m_slack, m_redCost, m_cBasisStat, m_rBasisStat, isOptimal);
      }
      catch( const SPxException& ex )
      {
         MSG_INFO1( (*spxout), (*spxout) << "Exception thrown while unsimplifying " << m_hist[k]->getName() << ":\n" << ex.what() << "\n" );
         throw SPxInternalCodeException("XMAISM00 Exception thrown during unsimply().");
      }

      m_hist[k]->~PostStep();
      spx_free(m_hist[k]);
      m_hist.reSize(k);
   }

   // for maximization problems, we have to switch signs of dual and reduced cost values back
   if(m_thesense == SPxLP::MAXIMIZE)
   {
      for(int j = 0; j < m_redCost.dim(); ++j)
         m_redCost[j] = -m_redCost[j];

      for(int i = 0; i < m_dual.dim(); ++i)
         m_dual[i] = -m_dual[i];
   }

   if( m_addedcols > 0 )
   {
      assert(m_prim.dim() >= m_addedcols);
      m_prim.reDim(m_prim.dim() - m_addedcols);
      m_redCost.reDim(m_redCost.dim() - m_addedcols);
      m_cBasisStat.reSize(m_cBasisStat.size() - m_addedcols);
      m_cIdx.reSize(m_cIdx.size() - m_addedcols);
   }

#ifdef CHECK_BASIC_DIM
   int numBasis = 0;
   for(int rs = 0; rs < m_rBasisStat.size(); ++rs)
   {
      if(m_rBasisStat[rs] == SPxSolver::BASIC)
         numBasis ++;
   }

   for(int cs = 0; cs < m_cBasisStat.size(); ++cs)
   {
      if(m_cBasisStat[cs] == SPxSolver::BASIC)
         numBasis ++;
   }

   if( numBasis != m_rBasisStat.size())
   {
      throw SPxInternalCodeException("XMAISM26 Dimension doesn't match after this step.");
   }
#endif

   m_hist.clear();
   m_postsolved = true;
}

// Pretty-printing of solver status.
std::ostream& operator<<(std::ostream& os, const SPxSimplifier::Result& status)
{
   switch( status )
   {
   case SPxSimplifier::OKAY:
      os << "SUCCESS";
      break;
   case SPxSimplifier::INFEASIBLE:
      os << "INFEASIBLE";
      break;
   case SPxSimplifier::DUAL_INFEASIBLE:
      os << "DUAL_INFEASIBLE";
      break;
   case SPxSimplifier::UNBOUNDED:
      os << "UNBOUNDED";
      break;
   case SPxSimplifier::VANISHED:
      os << "VANISHED";
      break;
   default:
      os << "UNKNOWN";
      break;
   }
   return os;
}

} //namespace soplex
