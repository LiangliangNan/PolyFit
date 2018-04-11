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

/**@file  spxscaler.h
 * @brief LP scaling base class.
 */
#ifndef _SPXSCALER_H_
#define _SPXSCALER_H_

#include <assert.h>

#include "spxdefines.h"
#include "dataarray.h"
#include "vector.h"
#include "svector.h"
#include "svset.h"
#include "dsvector.h"
#include "dvector.h"
#include <vector>

namespace soplex
{

template < class R >
class SPxLPBase;
/**@brief   LP scaler abstract base class.
   @ingroup Algo

   Instances of classes derived from SPxScaler may be loaded to SoPlex in
   order to scale LPs before solving them. SoPlex will load() itself to
   the SPxScaler and then call #scale(). Generally any SPxLP can be
   loaded to a SPxScaler for #scale()%ing it. The scaling can
   be undone by calling unscale().

   Mathematically, the scaling of a constraint matrix A can be written
   as \f$ A' = R A C \f$, with \f$ R \f$ and \f$ C \f$, being diagonal matrices
   corresponding to the row and column scale factors, respectively. Besides the
   constraints matrix, also the upper and lower bounds of both columns and rows
   need to be scaled.

   Note that by default scaling is performed both before and after presolving and
   the former scaling factors are retained during branch-and-bound (persistent scaling).
   However, while within SoPlex the scaled problem is used, data accessed through
   the soplex.cpp interface is provided w.r.t. the original problem (i.e., in unscaled form).
   For instance, consider a scaled constraints matrix A' that is extended by artificial slack
   variables to the matrix (A',I).
   A basis \f$ B' = [(A',I)P]_{[1:m][1:m] }\f$ (with P being a permutation matrix)
   for the scaled problem corresponds to the basis
   \f$ B = R^{-1} [(A',I)P]_{[1:m][1:m]} [P^{T} \tilde{C}^{-1} P]_{[1:m][1:m] } \f$. In
   this equation, \f$ \tilde{C} \f$is of the form

   \f[
    \begin{array}{cc}
         C & 0 \\
         O & R^{-1}
   \end{array}$
    \f]

   Note that in SoPlex only scaling factors \f$ 2^k, k \in \mathbb{Z} \f$ are used.


*/

class SPxScaler
{
protected:

   //-------------------------------------
   /**@name Data */
   //@{
   const char*        m_name;      ///< Name of the scaler
   DataArray < int >* m_activeColscaleExp; ///< pointer to currently active column scaling factors
   DataArray < int >* m_activeRowscaleExp; ///< pointer to currently active row scaling factors
   bool               m_colFirst;  ///< do column scaling first 
   bool               m_doBoth;    ///< do columns and rows
   SPxOut*            spxout;      ///< message handler
   //@}

   //-------------------------------------
   /**@name Protected helpers */
   //@{

   /// clear and setup scaling arrays in the LP
   virtual void setup(SPxLPBase<Real>& lp);
   //@}

public:

   /// compute a single scaling vector , e.g. of a newly added row
   virtual int computeScaleExp(const SVector& vec, const DataArray<int>& oldScaleExp) const;

#ifndef SOPLEX_LEGACY
   virtual int computeScaleExp(const SVectorBase<Rational>& vec, const DataArray<int>& oldScaleExp) const;
#endif

   /// applies m_colscale and m_rowscale to the \p lp.
   virtual void applyScaling(SPxLPBase<Real>& lp);

   friend std::ostream& operator<<(std::ostream& s, const SPxScaler& sc);

   //-------------------------------------
   /**@name Construction / destruction */
   //@{
   /// constructor
   explicit SPxScaler(const char* name, bool colFirst = false, bool doBoth = true, SPxOut* spxout = NULL);
   /// copy constructor
   SPxScaler(const SPxScaler& );
   /// assignment operator
   SPxScaler& operator=(const SPxScaler& );
   /// destructor.
   virtual ~SPxScaler();
   /// clone function for polymorphism
   virtual SPxScaler* clone() const = 0;
   //@}

   //-------------------------------------
   /**@name Access / modification */
   //@{
   /// get name of scaler
   virtual const char* getName() const;
   /// set scaling order
   virtual void setOrder(bool colFirst); 
   /// set wether column and row scaling should be performed
   virtual void setBoth(bool both); 
   /// set message handler
   virtual void setOutstream(SPxOut& newOutstream)
   {
      spxout = &newOutstream;
   }
   /// set real parameter
   virtual void setRealParam(Real param, const char* name = "realparam");
   /// set int parameter
   virtual void setIntParam(int param, const char* name = "intparam");
   //@}

   //-------------------------------------
   /**@name Scaling */
   //@{
   /// scale SPxLP.
   virtual void scale(SPxLPBase<Real>& lp, bool persistent = true) = 0;
   /// unscale SPxLP
   virtual void unscale(SPxLPBase<Real>& lp);
   /// returns scaling factor for column \p i
   virtual int getColScaleExp(int i) const;
   /// returns scaling factor for row \p i
   virtual int getRowScaleExp(int i) const;
   /// gets unscaled column \p i
   virtual void getColUnscaled(const SPxLPBase<Real>& lp, int i, DSVector& vec) const;
   /// returns maximum absolute value of unscaled column \p i
   virtual Real getColMaxAbsUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// returns minumum absolute value of unscaled column \p i
   virtual Real getColMinAbsUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// returns unscaled upper bound \p i
   virtual Real upperUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// returns unscaled upper bound vector of \p lp
   virtual void getUpperUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const;
   /// returns unscaled lower bound \p i
   virtual Real lowerUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// gets unscaled lower bound vector
   virtual void getLowerUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const;
   /// returns unscaled objective function coefficient of \p i
   virtual Real maxObjUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// gets unscaled objective function
   virtual void getMaxObjUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const;
   /// returns unscaled row \p i
   virtual void getRowUnscaled(const SPxLPBase<Real>& lp, int i, DSVector& vec) const;
   /// returns maximum absolute value of unscaled row \p i
   virtual Real getRowMaxAbsUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// returns minimum absolute value of unscaled row \p i
   virtual Real getRowMinAbsUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// returns unscaled right hand side \p i
   virtual Real rhsUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// gets unscaled right hand side vector
   virtual void getRhsUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const;
   /// returns unscaled left hand side \p i of \p lp
   virtual Real lhsUnscaled(const SPxLPBase<Real>& lp, int i) const;
   /// returns unscaled left hand side vector of \p lp
   virtual void getLhsUnscaled(const SPxLPBase<Real>& lp, Vector& vec) const;
   /// returns unscaled coefficient of \p lp
   virtual Real getCoefUnscaled(const SPxLPBase<Real>& lp, int row, int col) const;
   /// unscale dense primal solution vector given in \p x. 
   virtual void unscalePrimal(const SPxLPBase<Real>& lp, Vector& x) const;
   /// unscale dense slack vector given in \p s.
   virtual void unscaleSlacks(const SPxLPBase<Real>& lp, Vector& s) const;
   /// unscale dense dual solution vector given in \p pi. 
   virtual void unscaleDual(const SPxLPBase<Real>& lp, Vector& pi) const;
   /// unscale dense reduced cost vector given in \p r.
   virtual void unscaleRedCost(const SPxLPBase<Real>& lp, Vector& r) const;
   /// unscale primal ray given in \p ray.
   virtual void unscalePrimalray(const SPxLPBase<Real>& lp, Vector& ray) const;
   /// unscale dual ray given in \p ray.
   virtual void unscaleDualray(const SPxLPBase<Real>& lp, Vector& ray) const;
   /// apply scaling to objective function vector \p origObj.
   virtual void scaleObj(const SPxLPBase<Real>& lp, VectorReal& origObj) const;
   /// returns scaled objective function coefficient \p origObj.
   virtual Real scaleObj(const SPxLPBase<Real>& lp, int i, Real origObj) const;
   /// returns scaled LP element in \p row and \p col.
   virtual Real scaleElement(const SPxLPBase<Real>& lp, int row, int col, Real val) const;
   /// returns scaled lower bound of column \p col.
   virtual Real scaleLower(const SPxLPBase<Real>& lp, int col, Real lower) const;
   /// returns scaled upper bound of column \p col.
   virtual Real scaleUpper(const SPxLPBase<Real>& lp, int col, Real upper) const;
   /// returns scaled left hand side of row \p row.
   virtual Real scaleLhs(const SPxLPBase<Real>& lp, int row, Real lhs) const;
   /// returns scaled right hand side of row \p row.
   virtual Real scaleRhs(const SPxLPBase<Real>& lp, int row, Real rhs) const;
   /// absolute smallest column scaling factor
   virtual Real minAbsColscale() const;
   /// absolute biggest column scaling factor
   virtual Real maxAbsColscale() const;
   /// absolute smallest row scaling factor
   virtual Real minAbsRowscale() const;
   /// absolute biggest row scaling factor
   virtual Real maxAbsRowscale() const;
   /// maximum ratio between absolute biggest and smallest element in any column.
   virtual Real maxColRatio(const SPxLPBase<Real>& lp) const;
   /// maximum ratio between absolute biggest and smallest element in any row.
   virtual Real maxRowRatio(const SPxLPBase<Real>& lp) const;
   /// round vector entries to power of 2
   void computeExpVec(const std::vector<Real>& vec, DataArray<int>& vecExp);
   //@}

   //-------------------------------------
   /**@name Debugging */
   //@{
   /// consistency check
   virtual bool isConsistent() const;
   //@}
};
} // namespace soplex
#endif // _SPXSCALER_H_
