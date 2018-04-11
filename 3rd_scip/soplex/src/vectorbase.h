/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996      Roland Wunderling                              */
/*                  1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  vectorbase.h
 * @brief Dense vector.
 */
#ifndef _VECTORBASE_H_
#define _VECTORBASE_H_

#include <assert.h>
#include <string.h>
#include <math.h>
#include <iostream>

namespace soplex
{
template < class R > class SVectorBase;
template < class R > class SSVectorBase;

/**@brief   Dense vector.
 * @ingroup Algebra
 *
 *  Class VectorBase provides dense linear algebra vectors.  It does not provide memory management for the %array of
 *  values. Instead, the constructor requires a pointer to a memory block large enough to fit the desired dimension of
 *  Real or Rational values.
 *
 *  After construction, the values of a VectorBase can be accessed with the subscript operator[]().  Safety is provided by
 *  qchecking of array bound when accessing elements with the subscript operator[]() (only when compiled without \c
 *  -DNDEBUG).
 *
 *  A VectorBase is distinguished from a simple array of %Reals or %Rationals by providing a set of mathematical
 *  operations.  Since VectorBase does not provide any memory management features, no operations are available that would
 *  require allocation of temporary memory space.
 *
 *  The following mathematical operations are provided by class VectorBase (VectorBase \p a, \p b; R \p x):
 *
 *  <TABLE>
 *  <TR><TD>Operation</TD><TD>Description   </TD><TD></TD>&nbsp;</TR>
 *  <TR><TD>\c -=    </TD><TD>subtraction   </TD><TD>\c a \c -= \c b </TD></TR>
 *  <TR><TD>\c +=    </TD><TD>addition      </TD><TD>\c a \c += \c b </TD></TR>
 *  <TR><TD>\c *     </TD><TD>scalar product</TD>
 *      <TD>\c x = \c a \c * \c b </TD></TR>
 *  <TR><TD>\c *=    </TD><TD>scaling       </TD><TD>\c a \c *= \c x </TD></TR>
 *  <TR><TD>maxAbs() </TD><TD>infinity norm </TD>
 *      <TD>\c a.maxAbs() == \f$\|a\|_{\infty}\f$ </TD></TR>
 *  <TR><TD>minAbs() </TD><TD> </TD>
 *      <TD>\c a.minAbs() == \f$\min |a_i|\f$ </TD></TR>
 *
 *  <TR><TD>length() </TD><TD>euclidian norm</TD>
 *      <TD>\c a.length() == \f$\sqrt{a^2}\f$ </TD></TR>
 *  <TR><TD>length2()</TD><TD>square norm   </TD>
 *      <TD>\c a.length2() == \f$a^2\f$ </TD></TR>
 *  <TR><TD>multAdd(\c x,\c b)</TD><TD>add scaled vector</TD>
 *      <TD> \c a +=  \c x * \c b </TD></TR>
 *  </TABLE>
 *
 *  When using any of these operations, the vectors involved must be of the same dimension.  Also an SVectorBase \c b is
 *  allowed if it does not contain nonzeros with index greater than the dimension of \c a.q
 */
template < class R >
class VectorBase
{
protected:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   /// Dimension of vector.
   int dimen;

   /// Values of vector.
   /** The memory block pointed to by val must at least have size dimen * sizeof(R). */
   R* val;

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction and assignment */
   //@{

   /// Constructor.
   /** There is no default constructor since the storage for a VectorBase must be provided externally.  Storage must be
    *  passed as a memory block val at construction. It must be large enough to fit at least dimen values.
    */
   VectorBase<R>(int p_dimen, R* p_val)
      : dimen(p_dimen)
      , val(p_val)
   {
      assert(dimen >= 0);
      assert(isConsistent());
   }

   /// Assignment operator.
   template < class S >
   VectorBase<R>& operator=(const VectorBase<S>& vec)
   {
      if( (VectorBase<S>*)this != &vec )
      {
         assert(dim() == vec.dim());

         for( int i = 0; i < dimen; i++ )
            val[i] = vec[i];

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   VectorBase<R>& operator=(const VectorBase<R>& vec)
   {
      if( this != &vec )
      {
         assert(dim() == vec.dim());

         for( int i = 0; i < dimen; i++ )
            val[i] = vec[i];

         assert(isConsistent());
      }

      return *this;
   }

   /// scale and assign
   VectorBase<Real>& scaleAssign(int scaleExp, const VectorBase<Real>& vec)
   {
      if( this != &vec )
      {
         assert(dim() == vec.dim());

         for( int i = 0; i < dimen; i++ )
            val[i] = spxLdexp(vec[i], scaleExp);

         assert(isConsistent());
      }

      return *this;
   }

   /// scale and assign
   VectorBase<Real>& scaleAssign(const int* scaleExp, const VectorBase<Real>& vec, bool negateExp = false)
   {
      if( this != &vec )
      {
         assert(dim() == vec.dim());

         if( negateExp)
         {
            for( int i = 0; i < dimen; i++ )
               val[i] = spxLdexp(vec[i], -scaleExp[i]);
         }
         else
         {
            for( int i = 0; i < dimen; i++ )
               val[i] = spxLdexp(vec[i], scaleExp[i]);
         }

         assert(isConsistent());
      }

      return *this;
   }


   /// Assignment operator.
   /** Assigning an SVectorBase to a VectorBase using operator=() will set all values to 0 except the nonzeros of \p vec.
    *  This is diffent in method assign().
    */
   template < class S >
   VectorBase<R>& operator=(const SVectorBase<S>& vec);

   /// Assignment operator.
   /** Assigning an SSVectorBase to a VectorBase using operator=() will set all values to 0 except the nonzeros of \p
    *  vec.  This is diffent in method assign().
    */
   /**@todo do we need this also in non-template version, because SSVectorBase can be automatically cast to VectorBase? */
   template < class S >
   VectorBase<R>& operator=(const SSVectorBase<S>& vec);

   /// Assign values of \p vec.
   /** Assigns all nonzeros of \p vec to the vector.  All other values remain unchanged. */
   template < class S >
   VectorBase<R>& assign(const SVectorBase<S>& vec);

   /// Assign values of \p vec.
   /** Assigns all nonzeros of \p vec to the vector.  All other values remain unchanged. */
   template < class S >
   VectorBase<R>& assign(const SSVectorBase<S>& vec);

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Access */
   //@{

   /// Dimension of vector.
   int dim() const
   {
      return dimen;
   }

   /// Return \p n 'th value by reference.
   R& operator[](int n)
   {
      assert(n >= 0 && n < dimen);
      return val[n];
   }

   /// Return \p n 'th value.
   const R& operator[](int n) const
   {
      assert(n >= 0 && n < dimen);
      return val[n];
   }

   /// Equality operator.
   friend bool operator==(const VectorBase<R>& vec1, const VectorBase<R>& vec2)
   {
      if( &vec1 == &vec2 )
         return true;
      else if( vec1.dim() != vec2.dim() )
         return false;
      else
      {
         for( int i = 0; i < vec1.dim(); i++ )
         {
            if( vec1[i] != vec2[i] )
               return false;
         }
      }

      return true;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Arithmetic operations */
   //@{

   /// Set vector to 0.
   void clear()
   {
      if( dimen > 0 )
      {
         for( int i = 0; i < dimen; i++ )
            val[i] = 0;
      }
   }

   /// Addition.
   template < class S >
   VectorBase<R>& operator+=(const VectorBase<S>& vec)
   {
      assert(dim() == vec.dim());
      assert(dim() == dimen);

      for( int i = 0; i < dimen; i++ )
         val[i] += vec[i];

      return *this;
   }

   /// Addition.
   template < class S >
   VectorBase<R>& operator+=(const SVectorBase<S>& vec);

   /// Addition.
   template < class S >
   VectorBase<R>& operator+=(const SSVectorBase<S>& vec);

   /// Subtraction.
   template < class S >
   VectorBase<R>& operator-=(const VectorBase<S>& vec)
   {
      assert(dim() == vec.dim());
      assert(dim() == dimen);

      for( int i = 0; i < dimen; i++ )
         val[i] -= vec[i];

      return *this;
   }

   /// Subtraction.
   template < class S >
   VectorBase<R>& operator-=(const SVectorBase<S>& vec);

   /// Subtraction.
   template < class S >
   VectorBase<R>& operator-=(const SSVectorBase<S>& vec);

   /// Scaling.
   template < class S >
   VectorBase<R>& operator*=(const S& x)
   {
      assert(dim() == dimen);

      for( int i = 0; i < dimen; i++ )
         val[i] *= x;

      return *this;
   }

   /// Division.
   template < class S >
   VectorBase<R>& operator/=(const S& x)
   {
      assert(x != 0);

      for( int i = 0; i < dim(); i++ )
         val[i] /= x;

      return *this;
   }

   /// Inner product.
   R operator*(const VectorBase<R>& vec) const
   {
      assert(vec.dim() == dimen);

      R x = 0.0;

      for( int i = 0; i < dimen; i++ )
         x += val[i] * vec.val[i];

      return x;
   }

   /// Inner product.
   R operator*(const SVectorBase<R>& vec) const;

   /// Inner product.
   R operator*(const SSVectorBase<R>& vec) const;

   /// Maximum absolute value, i.e., infinity norm.
   R maxAbs() const
   {
      assert(dim() > 0);
      assert(dim() == dimen);

      R maxi = 0.0;

      for( int i = 0; i < dimen; i++ )
      {
         R x = spxAbs(val[i]);

         if( x > maxi )
            maxi = x;
      }

      assert(maxi >= 0.0);

      return maxi;
   }

   /// Minimum absolute value.
   R minAbs() const
   {
      assert(dim() > 0);
      assert(dim() == dimen);

      R mini = spxAbs(val[0]);

      for( int i = 1; i < dimen; i++ )
      {
         R x = spxAbs(val[i]);

         if( x < mini )
            mini = x;
      }

      assert(mini >= 0.0);

      return mini;
   }

   /// Floating point approximation of euclidian norm (without any approximation guarantee).
   Real length() const
   {
      return spxSqrt((Real)length2());
   }

   /// Squared norm.
   R length2() const
   {
      return (*this) * (*this);
   }

   /// Addition of scaled vector.
   template < class S, class T >
   VectorBase<R>& multAdd(const S& x, const VectorBase<T>& vec)
   {
      assert(vec.dim() == dimen);

      for( int i = 0; i < dimen; i++ )
         val[i] += x * vec.val[i];

      return *this;
   }

   /// Addition of scaled vector.
   template < class S, class T >
   VectorBase<R>& multAdd(const S& x, const SVectorBase<T>& vec);

   /// Subtraction of scaled vector.
   template < class S, class T >
   VectorBase<R>& multSub(const S& x, const SVectorBase<T>& vec);

   /// Addition of scaled vector.
   template < class S, class T >
   VectorBase<R>& multAdd(const S& x, const SSVectorBase<T>& vec);

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Utilities */
   //@{

   /// Conversion to C-style pointer.
   /** This function serves for using a VectorBase in an C-style function. It returns a pointer to the first value of
    *  the array.
    *
    *  @todo check whether this non-const c-style access should indeed be public
    */
   R* get_ptr()
   {
      return val;
   }

   /// Conversion to C-style pointer.
   /** This function serves for using a VectorBase in an C-style function. It returns a pointer to the first value of
    *  the array.
    */
   const R* get_const_ptr() const
   {
      return val;
   }

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( dim() > 0 && val == 0 )
         return MSGinconsistent("VectorBase");
#endif

      return true;
   }

   //@}

};



/// Assignment operator (specialization for Real).
template <>
inline
VectorBase<Real>& VectorBase<Real>::operator=(const VectorBase<Real>& vec)
{
   if( this != &vec )
   {
      assert(dim() == vec.dim());

      memcpy(val, vec.val, (unsigned int)dimen*sizeof(Real));

      assert(isConsistent());
   }
   return *this;
}



/// Assignment operator (specialization for Real).
template <>
template <>
inline
VectorBase<Real>& VectorBase<Real>::operator=(const VectorBase<Rational>& vec)
{
   if( (VectorBase<Rational>*)this != &vec )
   {
      assert(dim() == vec.dim());

      for( int i = 0; i < dimen; i++ )
         val[i] = Real(vec[i]);

      assert(isConsistent());
   }

   return *this;
}



/// Set vector to 0 (specialization for Real).
template<>
inline
void VectorBase<Real>::clear()
{
   if( dimen > 0 )
      memset(val, 0, (unsigned int)dimen * sizeof(Real));
}



#ifndef SOPLEX_LEGACY
/// Inner product.
template<>
inline
Rational VectorBase<Rational>::operator*(const VectorBase<Rational>& vec) const
{
   assert(vec.dim() == dimen);

   if( dimen <= 0 || vec.dim() <= 0 )
      return 0;

   Rational x = val[0];
   x *= vec.val[0];

   for( int i = 1; i < dimen && i < vec.dim(); i++ )
      x.addProduct(val[i], vec.val[i]);

   return x;
}
#endif

} // namespace soplex
#endif // _VECTORBASE_H_
