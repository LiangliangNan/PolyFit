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

/**@file  lprowbase.h
 * @brief (In)equality for LPs.
 */
#ifndef _LPROWBASE_H_
#define _LPROWBASE_H_

#include <assert.h>

#include "spxdefines.h"
#include "basevectors.h"

namespace soplex
{
/**@brief   (In)equality for LPs.
 * @ingroup Algo
 *
 *  Class LPRowBase provides constraints for linear programs in the form \f[ l \le a^Tx \le r, \f] where \em a is a
 *  DSVector. \em l is referred to as %left hand side, \em r as %right hand side and \em a as \em row \em vector or the
 *  constraint vector. \em l and \em r may also take values \f$\pm\f$ #infinity.  This static member is predefined, but
 *  may be overridden to meet the needs of the LP solver to be used.
 *
 *  LPRowBases allow to specify regular inequalities of the form \f[ a^Tx \sim \alpha, \f] where \f$\sim\f$ can take any
 *  value of \f$\le, =, \ge\f$, by setting rhs and lhs to the same value or setting one of them to \f$\infty\f$.
 *
 *  Since constraints in the regular form occur often, LPRowBases offers methods type() and value() for retreiving
 *  \f$\sim\f$ and \f$\alpha\f$ of an LPRowBase in this form, respectively. Also, a constructor for LPRowBases given in
 *  regular form is provided.
 */
template < class R >
class LPRowBase
{
   template < class S > friend class LPRowBase;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   R left;                 ///< left-hand side of the constraint
   R right;                ///< right-hand side of the constraint
   R object;               ///< objective coefficient of corresponding slack variable s = vec times primal
   DSVectorBase<R> vec;    ///< the row vector

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Types */
   //@{

   /// (In)Equality type of an LP row.
   /** #LPRowBase%s may be of one of the following Types. This datatype may be used for constructing new #LPRowBase%s in the
    *  regular form.
    */
   enum Type
   {
      LESS_EQUAL,          ///< \f$a^Tx \le \alpha\f$.
      EQUAL,               ///< \f$a^Tx = \alpha\f$.
      GREATER_EQUAL,       ///< \f$a^Tx \ge \alpha\f$.
      RANGE                ///< \f$\lambda \le a^Tx \le \rho\f$.
   };

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction / destruction */
   //@{

   /// Constructs LPRowBase with a vector ready to hold \p defDim nonzeros.
   explicit LPRowBase<R>(int defDim = 0)
      : left(0), right(infinity), object(0), vec(defDim)
   {
      assert(isConsistent());
   }

   /// Copy constructor.
   LPRowBase<R>(const LPRowBase<R>& row)
      : left(row.left), right(row.right), object(row.object), vec(row.vec)
   {
      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   LPRowBase<R>(const LPRowBase<S>& row)
      : left(row.left), right(row.right), object(row.object), vec(row.vec)
   {
      assert(isConsistent());
   }

   /// Constructs LPRowBase with the given left-hand side, right-hand side and rowVector.
   LPRowBase<R>(const R& p_lhs, const SVectorBase<R>& p_rowVector, const R& p_rhs, const R& p_obj = 0)
      : left(p_lhs), right(p_rhs), object(p_obj), vec(p_rowVector)
   {
      assert(isConsistent());
   }

   /// Constructs LPRowBase from passed \p rowVector, \p type and \p value.
   LPRowBase<R>(const SVectorBase<R>& p_rowVector, Type p_type, const R& p_value, const R& p_obj = 0)
      : object(p_obj), vec(p_rowVector)
   {
      switch( p_type )
      {
      case LESS_EQUAL:
         left = -infinity;
         right = p_value;
         break;
      case EQUAL:
         left = p_value;
         right = p_value;
         break;
      case GREATER_EQUAL:
         left = p_value;
         right = infinity;
         break;
      default:
         throw SPxInternalCodeException("XLPROW03 This should never happen.");
      }

      assert(isConsistent());
   }

   /// Destructor.
   ~LPRowBase()
   {}

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Access / modification */
   //@{

   /// Gets type of row.
   Type type() const
   {
      if( rhs() >= infinity )
         return GREATER_EQUAL;

      if( lhs() <= -infinity )
         return LESS_EQUAL;

      if( lhs() == rhs() )
         return EQUAL;

      return RANGE;
   }

   /// Sets type of (in)equality.
   void setType(Type p_type)
   {
      switch( p_type )
      {
      case LESS_EQUAL:
         left = -infinity;
         break;
      case EQUAL:
         if( lhs() > -infinity )
            right = lhs();
         else
            left = rhs();
         break;
      case GREATER_EQUAL:
         right = infinity;
         break;
      case RANGE:
         MSG_ERROR( std::cerr << "ELPROW01 RANGE not supported in LPRow::setType()"
            << std::endl; )
            throw SPxInternalCodeException("XLPROW01 This should never happen.");
      default:
         throw SPxInternalCodeException("XLPROW02 This should never happen.");
      }
   }

   /// Right hand side value of (in)equality.
   /** This method returns \f$\alpha\f$ for a LPRowBase in regular form.  However, value() may only be called for
    *  LPRowBase%s with type() != \c RANGE.
    */
   R value() const
   {
      assert(type() != RANGE);

      return (rhs() < infinity) ? rhs() : lhs();
   }

   /// Left-hand side value.
   R lhs() const
   {
      return left;
   }

   /// Sets left-hand side value.
   void setLhs(const R& p_left)
   {
      left = p_left;
   }

   /// Right-hand side value.
   R rhs() const
   {
      return right;
   }

   /// Sets right-hand side value.
   void setRhs(const R& p_right)
   {
      right = p_right;
   }

   /// Objective coefficient value.
   R obj() const
   {
      return object;
   }

   /// Sets objective coefficient value.
   void setObj(const R& p_obj)
   {
      object = p_obj;
   }

   /// Constraint row vector.
   const SVectorBase<R>& rowVector() const
   {
      return vec;
   }

   /// access constraint row vector.
   void setRowVector(const DSVectorBase<R>& p_vec)
   {
      vec = p_vec;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Consistency check */
   //@{

   /// Checks consistency.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      return vec.isConsistent();
#else
      return true;
#endif
   }

   //@}
};

} // namespace soplex
#endif // _LPROWBASE_H_
