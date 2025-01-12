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

/**@file  spxmainsm.h
 * @brief General methods in LP preprocessing.
 */
#ifndef _SPXMAINSM_H_
#define _SPXMAINSM_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxsimplifier.h"
#include "array.h"
#include "exceptions.h"

namespace soplex
{
//---------------------------------------------------------------------
//  class SPxMainSM
//---------------------------------------------------------------------

/**@brief   LP simplifier for removing uneccessary row/columns.
   @ingroup Algo

   This #SPxSimplifier is mainly based on the paper "Presolving in
   linear programming" by E. Andersen and K. Andersen (Mathematical
   Programming, 1995).  It implements all proposed methods and some
   other preprocessing techniques for removing redundant rows and
   columns and bounds.  Also infeasibility and unboundedness may be
   detected.

   Removed are:
   - empty rows / columns
   - unconstraint rows
   - row singletons
   - forcing rows
   - zero objective column singletons
   - (implied) free column singletons
   - doubleton equations combined with a column singleton
   - (implicitly) fixed columns
   - redundant lhs / rhs
   - redundant variable bounds
   - variables that are free in one direction
   - (weakly) dominated columns
   - duplicate rows / columns
*/
class SPxMainSM : public SPxSimplifier
{
private:
   //---------------------------------------------------------------------
   //  class PostsolveStep
   //---------------------------------------------------------------------

   /**@brief   Base class for postsolving operations.
      @ingroup Algo

      Class #PostStep is an abstract base class providing the
      interface for operations in the postsolving process.
   */
   class PostStep
   {
   private:
      /// name of the simplifier
      const char* m_name;
      /// number of cols
      int nCols;
      /// number of rows
      int nRows;

   public:
      /// constructor.
      PostStep(const char* p_name, int nR = 0, int nC = 0)
         : m_name(p_name)
         , nCols(nC)
         , nRows(nR)
      {}
      /// copy constructor.
      PostStep(const PostStep& old)
         : m_name(old.m_name)
         , nCols(old.nCols)
         , nRows(old.nRows)
      {}
      /// assignment operator
      PostStep& operator=(const PostStep& /*rhs*/)
      {
         return *this;
      }
      /// destructor.
      virtual ~PostStep()
      {
         m_name = 0;
      }
      /// get name of simplifying step.
      virtual const char* getName() const
      {
         return m_name;
      }
      /// clone function for polymorphism
      virtual PostStep* clone() const = 0;
      /// executes the postsolving.
      virtual void execute(
        DVector& x,                                 //*< Primal solution vector */
        DVector& y,                                 //*< Dual solution vector */
        DVector& s,                                 //*< Vector of slacks */
        DVector& r,                                 //*< Reduced cost vector */
        DataArray<SPxSolver::VarStatus>& cBasis,    //*< Basis status of column basis */
        DataArray<SPxSolver::VarStatus>& rBasis,    //*< Basis status of row basis */
        bool isOptimal
        ) const = 0;

      virtual bool checkBasisDim(DataArray<SPxSolver::VarStatus> rows,  DataArray<SPxSolver::VarStatus> cols) const;

      static Real eps()
      {
         return 1e-6;
      }
   };

   /**@brief   Postsolves row objectives.
      @ingroup Algo
   */
   class RowObjPS : public PostStep
   {
   private:
      int m_i; ///< row index
      int m_j; ///< slack column index

   public:
      ///
   RowObjPS(const SPxLP& lp, int _i, int _j)
         : PostStep("RowObj", lp.nRows(), lp.nCols())
         , m_i(_i)
         , m_j(_j)
      {}
      /// copy constructor
      RowObjPS(const RowObjPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_j(old.m_j)
      {}
      /// assignment operator
      RowObjPS& operator=( const RowObjPS& rhs)
      {
         if(this != &rhs)
         {
            m_i = rhs.m_i;
            m_j = rhs.m_j;
         }

         return *this;
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new RowObjPS(*this);
      }
   };

   /**@brief   Postsolves unconstraint constraints.
      @ingroup Algo
   */
   class FreeConstraintPS : public PostStep
   {
   private:
      int m_i;
      int m_old_i;
      DSVector  m_row;
      Real m_row_obj;

   public:
      ///
   FreeConstraintPS(const SPxLP& lp, int _i)
         : PostStep("FreeConstraint", lp.nRows(), lp.nCols())
         , m_i(_i)
         , m_old_i(lp.nRows()-1)
         , m_row(lp.rowVector(_i))
         , m_row_obj(lp.rowObj(_i))
      {}
      /// copy constructor
      FreeConstraintPS(const FreeConstraintPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_old_i(old.m_old_i)
         , m_row(old.m_row)
         , m_row_obj(old.m_row_obj)
      {}
      /// assignment operator
      FreeConstraintPS& operator=( const FreeConstraintPS& rhs)
      {
         if(this != &rhs)
         {
            m_i = rhs.m_i;
            m_old_i = rhs.m_old_i;
            m_row = rhs.m_row;
            m_row_obj = rhs.m_row_obj;
         }

         return *this;
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new FreeConstraintPS(*this);
      }
   };

   /**@brief   Postsolves empty constraints.
      @ingroup Algo
   */
   class EmptyConstraintPS : public PostStep
   {
   private:
      int m_i;
      int m_old_i;
      Real m_row_obj;

   public:
      ///
      EmptyConstraintPS(const SPxLP& lp, int _i)
         : PostStep("EmptyConstraint", lp.nRows(), lp.nCols())
         , m_i(_i)
         , m_old_i(lp.nRows()-1)
         , m_row_obj(lp.rowObj(_i))
      {}
      /// copy constructor
      EmptyConstraintPS(const EmptyConstraintPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_old_i(old.m_old_i)
         , m_row_obj(old.m_row_obj)
      {}
      /// assignment operator
      EmptyConstraintPS& operator=( const EmptyConstraintPS& rhs)
      {
         if(this != &rhs)
         {
            m_i = rhs.m_i;
            m_old_i = rhs.m_old_i;
            m_row_obj = rhs.m_row_obj;
         }

         return *this;
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new EmptyConstraintPS(*this);
      }
   };

   /**@brief   Postsolves row singletons.
      @ingroup Algo
   */
   class RowSingletonPS : public PostStep
   {
   private:
      const int  m_i;
      const int  m_old_i;
      const int  m_j;
      const Real m_lhs;
      const Real m_rhs;
      const bool m_strictLo;
      const bool m_strictUp;
      const bool m_maxSense;
      const Real m_obj;
      DSVector   m_col;
      const Real m_newLo;
      const Real m_newUp;
      const Real m_oldLo;
      const Real m_oldUp;
      const Real m_row_obj;

   public:
      ///
      RowSingletonPS(const SPxLP& lp, int _i, int _j, bool strictLo, bool strictUp,
                     Real newLo, Real newUp, Real oldLo, Real oldUp)
         : PostStep("RowSingleton", lp.nRows(), lp.nCols())
         , m_i(_i)
         , m_old_i(lp.nRows()-1)
         , m_j(_j)
         , m_lhs(lp.lhs(_i))
         , m_rhs(lp.rhs(_i))
         , m_strictLo(strictLo)
         , m_strictUp(strictUp)
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_obj(lp.spxSense() == SPxLP::MINIMIZE ? lp.obj(_j) : -lp.obj(_j))
         , m_col(lp.colVector(_j))
         , m_newLo(newLo)
         , m_newUp(newUp)
         , m_oldLo(oldLo)
         , m_oldUp(oldUp)
         , m_row_obj(lp.rowObj(_i))
      {}
      /// copy constructor
      RowSingletonPS(const RowSingletonPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_old_i(old.m_old_i)
         , m_j(old.m_j)
         , m_lhs(old.m_lhs)
         , m_rhs(old.m_rhs)
         , m_strictLo(old.m_strictLo)
         , m_strictUp(old.m_strictUp)
         , m_maxSense(old.m_maxSense)
         , m_obj(old.m_obj)
         , m_col(old.m_col)
         , m_newLo(old.m_newLo)
         , m_newUp(old.m_newUp)
         , m_oldLo(old.m_oldLo)
         , m_oldUp(old.m_oldUp)
         , m_row_obj(old.m_row_obj)
      {}
      /// assignment operator
      RowSingletonPS& operator=( const RowSingletonPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new RowSingletonPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves forcing constraints.
      @ingroup Algo
   */
   class ForceConstraintPS : public PostStep
   {
   private:
      const int       m_i;
      const int       m_old_i;
      const Real      m_lRhs;
      DSVector        m_row;
      DataArray<Real> m_objs;
      DataArray<bool> m_fixed;
      Array<DSVector> m_cols;
      const bool      m_lhsFixed;
      const bool      m_maxSense;
      DataArray<Real> m_oldLowers;
      DataArray<Real> m_oldUppers;
      const Real      m_lhs;
      const Real      m_rhs;
      const Real      m_rowobj;

   public:
      ///
      ForceConstraintPS(const SPxLP& lp, int _i, bool lhsFixed, DataArray<bool>& fixCols, DataArray<Real>& lo, DataArray<Real>& up)
         : PostStep("ForceConstraint", lp.nRows(), lp.nCols())
         , m_i(_i)
         , m_old_i(lp.nRows()-1)
         , m_lRhs(lhsFixed ? lp.lhs(_i) : lp.rhs(_i))
         , m_row(lp.rowVector(_i))
         , m_objs(lp.rowVector(_i).size())
         , m_fixed(fixCols)
         , m_cols(lp.rowVector(_i).size())
         , m_lhsFixed(lhsFixed)
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_oldLowers(lo)
         , m_oldUppers(up)
         , m_lhs(lp.lhs(_i))
         , m_rhs(lp.rhs(_i))
         , m_rowobj(lp.rowObj(_i))
      {
         for(int k = 0; k < m_row.size(); ++k)
         {
            m_objs[k] = (lp.spxSense() == SPxLP::MINIMIZE ? lp.obj(m_row.index(k)) : -lp.obj(m_row.index(k)));
            m_cols[k] = lp.colVector(m_row.index(k));
         }
      }
      /// copy constructor
      ForceConstraintPS(const ForceConstraintPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_old_i(old.m_old_i)
         , m_lRhs(old.m_lRhs)
         , m_row(old.m_row)
         , m_objs(old.m_objs)
         , m_fixed(old.m_fixed)
         , m_cols(old.m_cols)
         , m_lhsFixed(old.m_lhsFixed)
         , m_maxSense(old.m_maxSense)
         , m_oldLowers(old.m_oldLowers)
         , m_oldUppers(old.m_oldUppers)
         , m_lhs(old.m_lhs)
         , m_rhs(old.m_rhs)
         , m_rowobj(old.m_rowobj)
      {}
      /// assignment operator
      ForceConstraintPS& operator=( const ForceConstraintPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
            m_objs = rhs.m_objs;
            m_fixed = rhs.m_fixed;
            m_cols = rhs.m_cols;
            m_oldLowers = rhs.m_oldLowers;
            m_oldUppers = rhs.m_oldUppers;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new ForceConstraintPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves variable fixing.
      @ingroup Algo
   */
   class FixVariablePS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_old_j;
      const Real m_val;
      const Real m_obj;
      const Real m_lower;
      const Real m_upper;
      const bool m_correctIdx; /// does the index mapping have to be updated in postsolving?
      DSVector   m_col;

   public:
      ///
      FixVariablePS(const SPxLP& lp, SPxMainSM& simplifier, int _j, const Real val, bool correctIdx = true)
         : PostStep("FixVariable", lp.nRows(), lp.nCols())
         , m_j(_j)
         , m_old_j(lp.nCols()-1)
         , m_val(val)
         , m_obj(lp.spxSense() == SPxLP::MINIMIZE ? lp.obj(_j) : -lp.obj(_j))
         , m_lower(lp.lower(_j))
         , m_upper(lp.upper(_j))
         , m_correctIdx(correctIdx)
         , m_col(lp.colVector(_j))
      {
         simplifier.addObjoffset(m_val*lp.obj(m_j));
      }
      /// copy constructor
      FixVariablePS(const FixVariablePS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_old_j(old.m_old_j)
         , m_val(old.m_val)
         , m_obj(old.m_obj)
         , m_lower(old.m_lower)
         , m_upper(old.m_upper)
         , m_correctIdx(old.m_correctIdx)
         , m_col(old.m_col)
      {}
      /// assignment operator
      FixVariablePS& operator=( const FixVariablePS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         return new FixVariablePS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves variable bound fixing.
      @ingroup Algo
   */
   class FixBoundsPS : public PostStep
   {
   private:
      const int            m_j;
      SPxSolver::VarStatus m_status;

   public:
      ///
      FixBoundsPS(const SPxLP& lp, int j, Real val)
         : PostStep("FixBounds", lp.nRows(), lp.nCols())
         , m_j(j)
      {
         if (EQrel(lp.lower(j), lp.upper(j), eps()))
            m_status = SPxSolver::FIXED;
         else if (EQrel(val, lp.lower(j), eps()))
            m_status = SPxSolver::ON_LOWER;
         else if (EQrel(val, lp.upper(j), eps()))
            m_status = SPxSolver::ON_UPPER;
         else if (lp.lower(j) <= -infinity && lp.upper(j) >= infinity)
            m_status = SPxSolver::ZERO;
         else
         {
            throw SPxInternalCodeException("XMAISM14 This should never happen.");
         }
      }
      /// copy constructor
      FixBoundsPS(const FixBoundsPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_status(old.m_status)
      {}
      /// assignment operator
      FixBoundsPS& operator=( const FixBoundsPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_status = rhs.m_status;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         FixBoundsPS* FixBoundsPSptr = 0;
         spx_alloc(FixBoundsPSptr);
         return new (FixBoundsPSptr) FixBoundsPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief Postsolves the case when constraints are removed due to a
             variable with zero objective that is free in one direction.
      @ingroup Algo
   */
   class FreeZeroObjVariablePS : public PostStep
   {
   private:
      const int       m_j;
      const int       m_old_j;
      const int       m_old_i;
      const Real      m_bnd;
      DSVector        m_col;
      DSVector        m_lRhs;
      DSVector        m_rowObj;
      Array<DSVector> m_rows;
      const bool      m_loFree;

   public:
      ///
      FreeZeroObjVariablePS(const SPxLP& lp, int _j, bool loFree, SVector col_idx_sorted)
         : PostStep("FreeZeroObjVariable", lp.nRows(), lp.nCols())
         , m_j(_j)
         , m_old_j(lp.nCols()-1)
         , m_old_i(lp.nRows()-1)
         , m_bnd(loFree ? lp.upper(_j) : lp.lower(_j))
         , m_col(col_idx_sorted)
         , m_lRhs(lp.colVector(_j).size())
         , m_rowObj(lp.colVector(_j).size())
         , m_rows(lp.colVector(_j).size())
         , m_loFree(loFree)
      {
         for(int k = 0; k < m_col.size(); ++k)
         {
            assert(isNotZero(m_col.value(k)));

            int r = m_col.index(k);
            if ((m_loFree  && m_col.value(k) > 0) ||
                (!m_loFree && m_col.value(k) < 0))
               m_lRhs.add(k, lp.rhs(r));
            else
               m_lRhs.add(k, lp.lhs(r));

            m_rows[k] = lp.rowVector(r);
            m_rowObj.add(k, lp.rowObj(r));
         }
      }
      /// copy constructor
      FreeZeroObjVariablePS(const FreeZeroObjVariablePS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_old_j(old.m_old_j)
         , m_old_i(old.m_old_i)
         , m_bnd(old.m_bnd)
         , m_col(old.m_col)
         , m_lRhs(old.m_lRhs)
         , m_rowObj(old.m_rowObj)
         , m_rows(old.m_rows)
         , m_loFree(old.m_loFree)
      {}
      /// assignment operator
      FreeZeroObjVariablePS& operator=( const FreeZeroObjVariablePS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
            m_lRhs = rhs.m_lRhs;
            m_rowObj = rhs.m_rowObj;
            m_rows = rhs.m_rows;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         FreeZeroObjVariablePS* FreeZeroObjVariablePSptr = 0;
         spx_alloc(FreeZeroObjVariablePSptr);
         return new (FreeZeroObjVariablePSptr) FreeZeroObjVariablePS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves column singletons with zero objective.
      @ingroup Algo
   */
   class ZeroObjColSingletonPS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_i;
      const int  m_old_j;
      const Real m_lhs;
      const Real m_rhs;
      const Real m_lower;
      const Real m_upper;
      DSVector   m_row;

    public:
      ///
      ZeroObjColSingletonPS(const SPxLP& lp, const SPxMainSM& , int _j, int _i)
         : PostStep("ZeroObjColSingleton", lp.nRows(), lp.nCols())
         , m_j(_j)
         , m_i(_i)
         , m_old_j(lp.nCols()-1)
         , m_lhs(lp.lhs(_i))
         , m_rhs(lp.rhs(_i))
         , m_lower(lp.lower(_j))
         , m_upper(lp.upper(_j))
         , m_row(lp.rowVector(_i))
      {}
      /// copy constructor
      ZeroObjColSingletonPS(const ZeroObjColSingletonPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_i(old.m_i)
         , m_old_j(old.m_old_j)
         , m_lhs(old.m_lhs)
         , m_rhs(old.m_rhs)
         , m_lower(old.m_lower)
         , m_upper(old.m_upper)
         , m_row(old.m_row)
      {}
      /// assignment operator
      ZeroObjColSingletonPS& operator=( const ZeroObjColSingletonPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         ZeroObjColSingletonPS* ZeroObjColSingletonPSptr = 0;
         spx_alloc(ZeroObjColSingletonPSptr);
         return new (ZeroObjColSingletonPSptr) ZeroObjColSingletonPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves free column singletons.
      @ingroup Algo
   */
   class FreeColSingletonPS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_i;
      const int  m_old_j;
      const int  m_old_i;
      const Real m_obj;
      const Real m_lRhs;
      const bool m_onLhs;
      const bool m_eqCons;
      DSVector   m_row;

   public:
      ///
      FreeColSingletonPS(const SPxLP& lp, SPxMainSM& simplifier, int _j, int _i, Real slackVal)
         : PostStep("FreeColSingleton", lp.nRows(), lp.nCols())
         , m_j(_j)
         , m_i(_i)
         , m_old_j(lp.nCols()-1)
         , m_old_i(lp.nRows()-1)
         , m_obj(lp.spxSense() == SPxLP::MINIMIZE ? lp.obj(_j) : -lp.obj(_j))
         , m_lRhs(slackVal)
         , m_onLhs(EQ(slackVal, lp.lhs(_i)))
         , m_eqCons(EQ(lp.lhs(_i), lp.rhs(_i)))
         , m_row(lp.rowVector(_i))
      {
         assert(m_row[m_j] != 0.0);
         simplifier.addObjoffset(m_lRhs*(lp.obj(m_j)/m_row[m_j]));
      }
      /// copy constructor
      FreeColSingletonPS(const FreeColSingletonPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_i(old.m_i)
         , m_old_j(old.m_old_j)
         , m_old_i(old.m_old_i)
         , m_obj(old.m_obj)
         , m_lRhs(old.m_lRhs)
         , m_onLhs(old.m_onLhs)
         , m_eqCons(old.m_eqCons)
         , m_row(old.m_row)
      {}
      /// assignment operator
      FreeColSingletonPS& operator=( const FreeColSingletonPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         FreeColSingletonPS* FreeColSingletonPSptr = 0;
         spx_alloc(FreeColSingletonPSptr);
         return new (FreeColSingletonPSptr) FreeColSingletonPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves doubleton equations combined with a column singleton.
      @ingroup Algo
   */
   class DoubletonEquationPS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_k;
      const int  m_i;
      const bool m_maxSense;
      const bool m_jFixed;
      const Real m_jObj;
      const Real m_kObj;
      const Real m_aij;
      const bool m_strictLo;
      const bool m_strictUp;
      const Real m_newLo;
      const Real m_newUp;
      const Real m_oldLo;
      const Real m_oldUp;
      const Real m_Lo_j;
      const Real m_Up_j;
      const Real m_lhs;
      const Real m_rhs;
      DSVector   m_col;

   public:
      ///
      DoubletonEquationPS(const SPxLP& lp, int _j, int _k, int _i, Real oldLo, Real oldUp)
         : PostStep("DoubletonEquation", lp.nRows(), lp.nCols())
         , m_j(_j)
         , m_k(_k)
         , m_i(_i)
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_jFixed(EQ(lp.lower(_j), lp.upper(_j)))
         , m_jObj(lp.spxSense() == SPxLP::MINIMIZE ? lp.obj(_j) : -lp.obj(_j))
         , m_kObj(lp.spxSense() == SPxLP::MINIMIZE ? lp.obj(_k) : -lp.obj(_k))
         , m_aij(lp.colVector(_j).value(0))
         , m_strictLo(lp.lower(_k) > oldLo)
         , m_strictUp(lp.upper(_k) < oldUp)
         , m_newLo(lp.lower(_k))
         , m_newUp(lp.upper(_k))
         , m_oldLo(oldLo)
         , m_oldUp(oldUp)
         , m_Lo_j(lp.lower(_j))
         , m_Up_j(lp.upper(_j))
         , m_lhs(lp.lhs(_i))
         , m_rhs(lp.rhs(_i))
         , m_col(lp.colVector(_k))
      {}
      /// copy constructor
      DoubletonEquationPS(const DoubletonEquationPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_k(old.m_k)
         , m_i(old.m_i)
         , m_maxSense(old.m_maxSense)
         , m_jFixed(old.m_jFixed)
         , m_jObj(old.m_jObj)
         , m_kObj(old.m_kObj)
         , m_aij(old.m_aij)
         , m_strictLo(old.m_strictLo)
         , m_strictUp(old.m_strictUp)
         , m_newLo(old.m_newLo)
         , m_newUp(old.m_newUp)
         , m_oldLo(old.m_oldLo)
         , m_oldUp(old.m_oldUp)
         , m_Lo_j(old.m_Lo_j)
         , m_Up_j(old.m_Up_j)
         , m_lhs(old.m_lhs)
         , m_rhs(old.m_rhs)
         , m_col(old.m_col)
      {}
      /// assignment operator
      DoubletonEquationPS& operator=( const DoubletonEquationPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_col = rhs.m_col;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         DoubletonEquationPS* DoubletonEquationPSptr = 0;
         spx_alloc(DoubletonEquationPSptr);
         return new (DoubletonEquationPSptr) DoubletonEquationPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves duplicate rows.
      @ingroup Algo
   */
   class DuplicateRowsPS : public PostStep
   {
   private:
      const int       m_i;
      const Real      m_i_rowObj;
      const int       m_maxLhsIdx;
      const int       m_minRhsIdx;
      const bool      m_maxSense;
      const bool      m_isFirst;
      const bool      m_isLast;
      const bool      m_fixed;
      const int       m_nCols;
      DSVector        m_scale;
      DSVector        m_rowObj;
      DataArray<int>  m_rIdxLocalOld;
      DataArray<int>  m_perm;
      DataArray<bool> m_isLhsEqualRhs;

   public:
      DuplicateRowsPS(const SPxLP& lp, int _i,
                      int maxLhsIdx, int minRhsIdx, const DSVector& dupRows,
                      const DataArray<Real> scale, const DataArray<int> perm, const DataArray<bool> isLhsEqualRhs,
                      bool isTheLast, bool isFixedRow, bool isFirst = false)
         : PostStep("DuplicateRows", lp.nRows(), lp.nCols())
         , m_i(_i)
         , m_i_rowObj(lp.rowObj(_i))
         , m_maxLhsIdx((maxLhsIdx == -1) ? -1 : maxLhsIdx)
         , m_minRhsIdx((minRhsIdx == -1) ? -1 : minRhsIdx)
         , m_maxSense(lp.spxSense() == SPxLP::MAXIMIZE)
         , m_isFirst(isFirst)
         , m_isLast(isTheLast)
         , m_fixed(isFixedRow)
         , m_nCols(lp.nCols())
         , m_scale(dupRows.size())
         , m_rowObj(dupRows.size())
         , m_rIdxLocalOld(dupRows.size())
         , m_perm(perm)
         , m_isLhsEqualRhs(isLhsEqualRhs)
      {
         Real rowScale = scale[_i];

         for(int k = 0; k < dupRows.size(); ++k)
         {
            m_scale.add(dupRows.index(k), rowScale / scale[dupRows.index(k)]);
            m_rowObj.add(dupRows.index(k), lp.rowObj(dupRows.index(k)));
            m_rIdxLocalOld[k] = dupRows.index(k);
         }
      }
      /// copy constructor
      DuplicateRowsPS(const DuplicateRowsPS& old)
         : PostStep(old)
         , m_i(old.m_i)
         , m_i_rowObj(old.m_i_rowObj)
         , m_maxLhsIdx(old.m_maxLhsIdx)
         , m_minRhsIdx(old.m_minRhsIdx)
         , m_maxSense(old.m_maxSense)
         , m_isFirst(old.m_isFirst)
         , m_isLast(old.m_isLast)
         , m_fixed(old.m_fixed)
         , m_nCols(old.m_nCols)
         , m_scale(old.m_scale)
         , m_rowObj(old.m_rowObj)
         , m_rIdxLocalOld(old.m_rIdxLocalOld)
         , m_perm(old.m_perm)
         , m_isLhsEqualRhs(old.m_isLhsEqualRhs)
      {}
      /// assignment operator
      DuplicateRowsPS& operator=( const DuplicateRowsPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_scale = rhs.m_scale;
            m_rowObj = rhs.m_rowObj;
            m_rIdxLocalOld = rhs.m_rIdxLocalOld;
            m_perm = rhs.m_perm;
            m_isLhsEqualRhs = rhs.m_isLhsEqualRhs;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         DuplicateRowsPS* DuplicateRowsPSptr = 0;
         spx_alloc(DuplicateRowsPSptr);
         return new (DuplicateRowsPSptr) DuplicateRowsPS(*this);
      }
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves duplicate columns.
      @ingroup Algo
   */
   class DuplicateColsPS : public PostStep
   {
   private:
      const int            m_j;
      const int            m_k;
      const Real           m_loJ;
      const Real           m_upJ;
      const Real           m_loK;
      const Real           m_upK;
      const Real           m_scale;
      const bool           m_isFirst;
      const bool           m_isLast;
      DataArray<int>       m_perm;

   public:
      DuplicateColsPS(const SPxLP& lp, int _j, int _k, Real scale, DataArray<int>  perm, bool isFirst = false, bool isTheLast = false)
         : PostStep("DuplicateCols", lp.nRows(), lp.nCols())
         , m_j(_j)
         , m_k(_k)
         , m_loJ(lp.lower(_j))
         , m_upJ(lp.upper(_j))
         , m_loK(lp.lower(_k))
         , m_upK(lp.upper(_k))
         , m_scale(scale)
         , m_isFirst(isFirst)
         , m_isLast(isTheLast)
         , m_perm(perm)
      {}
      /// copy constructor
      DuplicateColsPS(const DuplicateColsPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_k(old.m_k)
         , m_loJ(old.m_loJ)
         , m_upJ(old.m_upJ)
         , m_loK(old.m_loK)
         , m_upK (old.m_upK)
         , m_scale (old.m_scale)
         , m_isFirst(old.m_isFirst)
         , m_isLast(old.m_isLast)
         , m_perm(old.m_perm)
      {}
      /// assignment operator
      DuplicateColsPS& operator=( const DuplicateColsPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         DuplicateColsPS* DuplicateColsPSptr = 0;
         spx_alloc(DuplicateColsPSptr);
         return new (DuplicateColsPSptr) DuplicateColsPS(*this);
      }
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves multi aggregation.
      @ingroup Algo
   */
   class MultiAggregationPS : public PostStep
   {
   private:
      const int  m_j;
      const int  m_i;
      const int  m_old_j;
      const int  m_old_i;
      const Real m_upper;
      const Real m_lower;
      const Real m_obj;
      const Real m_const;
      const bool m_onLhs;
      const bool m_eqCons;
      DSVector   m_row;
      DSVector   m_col;

   public:
      ///
      MultiAggregationPS(const SPxLP& lp, SPxMainSM& simplifier, int _i, int _j, Real constant)
         : PostStep("MultiAggregation", lp.nRows(), lp.nCols())
         , m_j(_j)
         , m_i(_i)
         , m_old_j(lp.nCols()-1)
         , m_old_i(lp.nRows()-1)
         , m_upper(lp.upper(_j))
         , m_lower(lp.lower(_j))
         , m_obj(lp.spxSense() == SPxLP::MINIMIZE ? lp.obj(_j) : -lp.obj(_j))
         , m_const(constant)
         , m_onLhs(EQ(constant, lp.lhs(_i)))
         , m_eqCons(EQ(lp.lhs(_i), lp.rhs(_i)))
         , m_row(lp.rowVector(_i))
         , m_col(lp.colVector(_j))
      {
         assert(m_row[m_j] != 0.0);
         simplifier.addObjoffset(m_obj*m_const/m_row[m_j]);
      }
      /// copy constructor
      MultiAggregationPS(const MultiAggregationPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_i(old.m_i)
         , m_old_j(old.m_old_j)
         , m_old_i(old.m_old_i)
         , m_upper(old.m_upper)
         , m_lower(old.m_lower)
         , m_obj(old.m_obj)
         , m_const(old.m_const)
         , m_onLhs(old.m_onLhs)
         , m_eqCons(old.m_eqCons)
         , m_row(old.m_row)
         , m_col(old.m_col)
      {}
      /// assignment operator
      MultiAggregationPS& operator=( const MultiAggregationPS& rhs)
      {
         if(this != &rhs)
         {
            PostStep::operator=(rhs);
            m_row = rhs.m_row;
            m_col = rhs.m_col;
         }

         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         MultiAggregationPS* MultiAggregationPSptr = 0;
         spx_alloc(MultiAggregationPSptr);
         return new (MultiAggregationPSptr) MultiAggregationPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };

   /**@brief   Postsolves variable bound tightening from pseudo objective propagation.
      @ingroup Algo
   */
   class TightenBoundsPS : public PostStep
   {
   private:
      const int            m_j;
      const Real           m_origupper;
      const Real           m_origlower;

   public:
      ///
      TightenBoundsPS(const SPxLP& lp, int j, Real origupper, Real origlower)
         : PostStep("TightenBounds", lp.nRows(), lp.nCols())
         , m_j(j)
         , m_origupper(origupper)
         , m_origlower(origlower)
      {
      }
      /// copy constructor
      TightenBoundsPS(const TightenBoundsPS& old)
         : PostStep(old)
         , m_j(old.m_j)
         , m_origupper(old.m_origupper)
         , m_origlower(old.m_origlower)
      {}
      /// assignment operator
      TightenBoundsPS& operator=( const TightenBoundsPS& rhs)
      {
         return *this;
      }
      /// clone function for polymorphism
      inline virtual PostStep* clone() const
      {
         TightenBoundsPS* TightenBoundsPSptr = 0;
         spx_alloc(TightenBoundsPSptr);
         return new (TightenBoundsPSptr) TightenBoundsPS(*this);
      }
      ///
      virtual void execute(DVector& x, DVector& y, DVector& s, DVector& r,
                           DataArray<SPxSolver::VarStatus>& cBasis, DataArray<SPxSolver::VarStatus>& rBasis, bool isOptimal) const;
   };
   // friends
   friend class FreeConstraintPS;
   friend class EmptyConstraintPS;
   friend class RowSingletonPS;
   friend class ForceConstraintPS;
   friend class FixVariablePS;
   friend class FixBoundsPS;
   friend class FreeZeroObjVariablePS;
   friend class ZeroObjColSingletonPS;
   friend class FreeColSingletonPS;
   friend class DoubletonEquationPS;
   friend class DuplicateRowsPS;
   friend class DuplicateColsPS;

private:
   //------------------------------------
   //**@name Types */
   //@{
   /// Different simplification steps.
   enum SimpleStep
   {
      EMPTY_ROW            =  0,
      FREE_ROW             =  1,
      SINGLETON_ROW        =  2,
      FORCE_ROW            =  3,
      EMPTY_COL            =  4,
      FIX_COL              =  5,
      FREE_ZOBJ_COL        =  6,
      ZOBJ_SINGLETON_COL   =  7,
      DOUBLETON_ROW        =  8,
      FREE_SINGLETON_COL   =  9,
      DOMINATED_COL        = 10,
      WEAKLY_DOMINATED_COL = 11,
      DUPLICATE_ROW        = 12,
      FIX_DUPLICATE_COL    = 13,
      SUB_DUPLICATE_COL    = 14,
      MULTI_AGG            = 15
   };
   //@}

   //------------------------------------
   //**@name Data */
   //@{
   ///
   DVector                         m_prim;       ///< unsimplified primal solution vector.
   DVector                         m_slack;      ///< unsimplified slack vector.
   DVector                         m_dual;       ///< unsimplified dual solution vector.
   DVector                         m_redCost;    ///< unsimplified reduced cost vector.
   DataArray<SPxSolver::VarStatus> m_cBasisStat; ///< basis status of columns.
   DataArray<SPxSolver::VarStatus> m_rBasisStat; ///< basis status of rows.
   DataArray<int>                  m_cIdx;       ///< column index vector in original LP.
   DataArray<int>                  m_rIdx;       ///< row index vector in original LP.
   DataArray<PostStep*>            m_hist;       ///< vector of presolve history.
   Array<DSVector>                 m_classSetRows; ///< stores parallel classes with non-zero colum entry
   Array<DSVector>                 m_classSetCols; ///< stores parallel classes with non-zero row entry
   Array<DSVector>                 m_dupRows;    ///< arrange duplicate rows using bucket sort w.r.t. their pClass values
   Array<DSVector>                 m_dupCols;    ///< arrange duplicate columns w.r.t. their pClass values
   bool                            m_postsolved; ///< status of postsolving.
   Real                            m_epsilon;    ///< epsilon zero.
   Real                            m_feastol;    ///< primal feasibility tolerance.
   Real                            m_opttol;     ///< dual feasibility tolerance.
   DataArray<int>                  m_stat;       ///< preprocessing history.
   SPxLP::SPxSense                 m_thesense;   ///< optimization sense.
   bool                            m_keepbounds;  ///< keep some bounds (for boundflipping)
   int                             m_addedcols;  ///< columns added by handleRowObjectives()
   Result                          m_result;     ///< result of the simplification.
   Real                            m_cutoffbound;  ///< the cutoff bound that is found by heuristics
   Real                            m_pseudoobj;    ///< the pseudo objective function value
   //@}

private:
   //------------------------------------
   //**@name Private helpers */
   //@{
   /// handle row objectives
   void handleRowObjectives(SPxLP& lp);

   /// handles extreme values by setting them to zero or infinity.
   void handleExtremes(SPxLP& lp);

   /// computes the minimum and maximum residual activity for a given row and column. If colNumber is set to -1, then
   //  the activity of the row is returned.
   void computeMinMaxResidualActivity(SPxLP& lp, int rowNumber, int colNumber, Real& minAct, Real& maxAct);

   /// calculate min/max value for the multi aggregated variables
   void computeMinMaxValues(SPxLP& lp, Real side, Real val, Real minRes, Real maxRes, Real& minVal, Real& maxVal);

   /// tries to find good lower bound solutions by applying some trivial heuristics
   void trivialHeuristic(SPxLP& lp);

   /// checks a solution for feasibility
   bool checkSolution(SPxLP& lp, DVector sol);

   /// tightens variable bounds by propagating the pseudo objective function value.
   void propagatePseudoobj(SPxLP& lp);

   /// removed empty rows and empty columns.
   Result removeEmpty(SPxLP& lp);

   /// remove row singletons.
   Result removeRowSingleton(SPxLP& lp, const SVector& row, int& i);

   /// performs simplification steps on the rows of the LP.
   Result simplifyRows(SPxLP& lp, bool& again);

   /// performs simplification steps on the columns of the LP.
   Result simplifyCols(SPxLP& lp, bool& again);

   /// performs simplification steps on the LP based on dual concepts.
   Result simplifyDual(SPxLP& lp, bool& again);

   /// performs multi-aggregations of variable based upon constraint activitu.
   Result multiaggregation(SPxLP& lp, bool& again);

   /// removes duplicate rows.
   Result duplicateRows(SPxLP& lp, bool& again);

   /// removes duplicate columns
   Result duplicateCols(SPxLP& lp, bool& again);

   /// handles the fixing of a variable. correctIdx is true iff the index mapping has to be updated.
   void fixColumn(SPxLP& lp, int i, bool correctIdx = true);

   /// removes a row in the LP.
   void removeRow(SPxLP& lp, int i)
   {
      m_rIdx[i] = m_rIdx[lp.nRows()-1];
      lp.removeRow(i);
   }
   /// removes a column in the LP.
   void removeCol(SPxLP& lp, int j)
   {
      m_cIdx[j] = m_cIdx[lp.nCols()-1];
      lp.removeCol(j);
   }
   /// returns for a given row index of the (reduced) LP the corresponding row index in the unsimplified LP.
   int rIdx(int i) const
   {
      return m_rIdx[i];
   }
   /// returns for a given column index of the (reduced) LP the corresponding column index in the unsimplified LP.
   int cIdx(int j) const
   {
      return m_cIdx[j];
   }
   ///
   Real epsZero() const
   {
      return m_epsilon;
   }
   ///
   Real feastol() const
   {
      return m_feastol;
   }
   ///
   Real opttol() const
   {
      return m_opttol;
   }
   //@}

public:

   //------------------------------------
   //**@name Constructors / destructors */
   //@{
   /// default constructor.
   SPxMainSM(Timer::TYPE ttype = Timer::USER_TIME)
      : SPxSimplifier("MainSM", ttype)
      , m_postsolved(0)
      , m_epsilon(DEFAULT_EPS_ZERO)
      , m_feastol(DEFAULT_BND_VIOL)
      , m_opttol(DEFAULT_BND_VIOL)
      , m_stat(16)
      , m_thesense(SPxLP::MAXIMIZE)
      , m_keepbounds(false)
      , m_addedcols(0)
      , m_result(OKAY)
      , m_cutoffbound(-infinity)
      , m_pseudoobj(-infinity)
   {}
   /// copy constructor.
   SPxMainSM(const SPxMainSM& old)
      : SPxSimplifier(old)
      , m_prim(old.m_prim)
      , m_slack(old.m_slack)
      , m_dual(old.m_dual)
      , m_redCost(old.m_redCost)
      , m_cBasisStat(old.m_cBasisStat)
      , m_rBasisStat(old.m_rBasisStat)
      , m_cIdx(old.m_cIdx)
      , m_rIdx(old.m_rIdx)
      , m_postsolved(old.m_postsolved)
      , m_epsilon(old.m_epsilon)
      , m_feastol(old.m_feastol)
      , m_opttol(old.m_opttol)
      , m_stat(old.m_stat)
      , m_thesense(old.m_thesense)
      , m_keepbounds(old.m_keepbounds)
      , m_addedcols(old.m_addedcols)
      , m_result(old.m_result)
      , m_cutoffbound(old.m_cutoffbound)
      , m_pseudoobj(old.m_pseudoobj)
   {
      // copy pointers in m_hist
      m_hist.reSize(0);
      for(int k = 0; k < old.m_hist.size(); ++k)
      {
         if(old.m_hist[k] != 0)
            m_hist.append(old.m_hist[k]->clone());
         else
            m_hist.append(0);
      }
   }
   /// assignment operator
   SPxMainSM& operator=( const SPxMainSM& rhs)
   {
      if(this != &rhs)
      {
         SPxSimplifier::operator=(rhs);
         m_prim = rhs.m_prim;
         m_slack = rhs.m_slack;
         m_dual = rhs.m_dual;
         m_redCost = rhs.m_redCost;
         m_cBasisStat = rhs.m_cBasisStat;
         m_rBasisStat = rhs.m_rBasisStat;
         m_cIdx = rhs.m_cIdx;
         m_rIdx = rhs.m_rIdx;
         m_postsolved = rhs.m_postsolved;
         m_epsilon = rhs.m_epsilon;
         m_feastol = rhs.m_feastol;
         m_opttol = rhs.m_opttol;
         m_stat = rhs.m_stat;
         m_thesense = rhs.m_thesense;
         m_keepbounds = rhs.m_keepbounds;
         m_addedcols = rhs.m_addedcols;
         m_result = rhs.m_result;
         m_cutoffbound = rhs.m_cutoffbound;
         m_pseudoobj = rhs.m_pseudoobj;

         // delete pointers in m_hist
         for(int k = 0; k < m_hist.size(); ++k)
         {
            m_hist[k]->~PostStep();
            spx_free(m_hist[k]);
         }

         m_hist.clear();

         // copy pointers in m_hist
         for(int k = 0; k < rhs.m_hist.size(); ++k)
         {
            if(rhs.m_hist[k] != 0)
               m_hist.append(rhs.m_hist[k]->clone());
            else
               m_hist.append(0);
         }
      }

      return *this;
   }
   /// destructor.
   virtual ~SPxMainSM()
   {
      // delete pointers in m_hist
      for(int k = 0; k < m_hist.size(); ++k)
      {
         if( m_hist[k] != 0 )
         {
            m_hist[k]->~PostStep();
            spx_free(m_hist[k]);
         }
      }
   }
   /// clone function for polymorphism
   inline virtual SPxSimplifier* clone() const
   {
      return new SPxMainSM(*this);
   }
   //@}

   //------------------------------------
   //**@name LP simplification */
   //@{
   /// simplify SPxLP \p lp with identical primal and dual feasibility tolerance.
   virtual Result simplify(SPxLP& lp, Real eps, Real delta)
   {
      return simplify(lp, eps, delta, delta);
   }
   /// simplify SPxLP \p lp with independent primal and dual feasibility tolerance.
   virtual Result simplify(SPxLP& lp, Real eps, Real ftol, Real otol, bool keepbounds = false);

   /// reconstructs an optimal solution for the unsimplified LP.
   virtual void unsimplify(const Vector& x, const Vector& y, const Vector& s, const Vector& r,
                           const SPxSolver::VarStatus rows[], const SPxSolver::VarStatus cols[], bool isOptimal = true);

   /// returns result status of the simplification
   virtual Result result() const
   {
      return m_result;
   }

   /// specifies whether an optimal solution has already been unsimplified.
   virtual bool isUnsimplified() const
   {
      return m_postsolved;
   }
   /// returns a reference to the unsimplified primal solution.
   virtual const Vector& unsimplifiedPrimal()
   {
      assert(m_postsolved);
      return m_prim;
   }
   /// returns a reference to the unsimplified dual solution.
   virtual const Vector& unsimplifiedDual()
   {
      assert(m_postsolved);
      return m_dual;
   }
   /// returns a reference to the unsimplified slack values.
   virtual const Vector& unsimplifiedSlacks()
   {
      assert(m_postsolved);
      return m_slack;
   }
   /// returns a reference to the unsimplified reduced costs.
   virtual const Vector& unsimplifiedRedCost()
   {
      assert(m_postsolved);
      return m_redCost;
   }
   /// gets basis status for a single row.
   virtual SPxSolver::VarStatus getBasisRowStatus(int i) const
   {
      assert(m_postsolved);
      return m_rBasisStat[i];
   }
   /// gets basis status for a single column.
   virtual SPxSolver::VarStatus getBasisColStatus(int j) const
   {
      assert(m_postsolved);
      return m_cBasisStat[j];
   }
   /// get optimal basis.
   virtual void getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[], const int rowsSize = -1, const int colsSize = -1) const
   {
      assert(m_postsolved);
      assert(rowsSize < 0 || rowsSize >= m_rBasisStat.size());
      assert(colsSize < 0 || colsSize >= m_cBasisStat.size());

      for(int i = 0; i < m_rBasisStat.size(); ++i)
         rows[i] = m_rBasisStat[i];

      for(int j = 0; j < m_cBasisStat.size(); ++j)
         cols[j] = m_cBasisStat[j];
   }
   //@}

private:
   //------------------------------------
   //**@name Types */
   //@{
   /// comparator for class SVector::Element: compare nonzeros according to value
   struct ElementCompare
   {
   public:
      ElementCompare() {}

      int operator()(const SVector::Element& e1, const SVector::Element& e2) const
      {
	 if (EQ(e1.val, e2.val))
            return 0;
         if (e1.val < e2.val)
	    return -1;
	 else // (e1.val > e2.val)
	    return 1;
      }
   };
   /// comparator for class SVector::Element: compare nonzeros according to index
   struct IdxCompare
   {
   public:
      IdxCompare() {}

      int operator()(const SVector::Element& e1, const SVector::Element& e2) const
      {
	 if (EQ(e1.idx, e2.idx))
            return 0;
         if (e1.idx < e2.idx)
	    return -1;
	 else // (e1.idx > e2.idx)
	    return 1;
      }
   };
   //@}
};

} // namespace soplex
#endif // _SPXMAINSM_H_
