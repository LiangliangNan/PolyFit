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

/**@file  basevectors.h
 * @brief Collection of dense, sparse, and semi-sparse vectors.
 */
#ifndef _BASEVECTORS_H_
#define _BASEVECTORS_H_

/* undefine SOPLEX_DEBUG flag from including files; if SOPLEX_DEBUG should be defined in this file, do so below */
#ifdef SOPLEX_DEBUG
#define SOPLEX_DEBUG_BASEVECTORS
#undef SOPLEX_DEBUG
#endif

#include "spxdefines.h"
#include "rational.h"
#include "vectorbase.h"
#include "dvectorbase.h"
#include "ssvectorbase.h"
#include "svectorbase.h"
#include "dsvectorbase.h"
#include "unitvectorbase.h"
#include "svsetbase.h"

#define SOPLEX_VECTOR_MARKER   1e-100

namespace soplex
{

// ---------------------------------------------------------------------------------------------------------------------
//  Methods of VectorBase
// ---------------------------------------------------------------------------------------------------------------------

/// Assignment operator.
/** Assigning an SVectorBase to a VectorBase using operator=() will set all values to 0 except the nonzeros of \p vec.
 *  This is different in method assign().
 */



template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::operator=(const SVectorBase<S>& vec)
{
   clear();

   for( int i = 0; i < vec.size(); ++i )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)] = vec.value(i);
   }

   assert(isConsistent());

   return *this;
}



/// Assign values of \p vec.
/** Assigns all nonzeros of \p vec to the vector.  All other values remain unchanged. */
template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::assign(const SVectorBase<S>& vec)
{
   for( int i = vec.size() - 1; i >= 0; --i )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)] = vec.value(i);
   }

   assert(isConsistent());

   return *this;
}



/// Assignment operator.
/** Assigning an SSVectorBase to a VectorBase using operator=() will set all values to 0 except the nonzeros of \p vec.
 *  This is different in method assign().
 */
template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::operator=(const SSVectorBase<S>& vec)
{
   if( vec.isSetup() )
   {
      clear();
      assign(vec);
   }
   else
      operator=(static_cast<const VectorBase<R>&>(vec));

   return *this;
}



/// Assign values of \p vec.
/** Assigns all nonzeros of \p vec to the vector.  All other values remain unchanged. */
template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::assign(const SSVectorBase<S>& vec)
{
   assert(vec.dim() <= dim());

   if (vec.isSetup())
   {
      const int* idx = vec.indexMem();

      for(int i = vec.size() - 1; i >= 0; --i)
      {
         val[*idx] = vec.val[*idx];
         idx++;
      }
   }
   else
      operator=(static_cast<const VectorBase<R>&>(vec));

   return *this;
}



/// Addition.
template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::operator+=(const SVectorBase<S>& vec)
{
   for( int i = vec.size() - 1; i >= 0; --i )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] += vec.value(i);
   }

   return *this;
}



/// Addition.
template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::operator+=(const SSVectorBase<S>& vec)
{
   assert(dim() == vec.dim());

   if ( vec.isSetup() )
   {
      for( int i = vec.size() - 1; i >= 0 ; --i )
         val[vec.index(i)] += vec.value(i);
   }
   else
   {
      for( int i = dim() - 1; i >= 0; --i )
         val[i] += vec[i];
   }

   return *this;
}



/// Subtraction.
template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::operator-=(const SVectorBase<S>& vec)
{
   for( int i = vec.size() - 1; i >= 0; --i )
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }

   return *this;
}



/// Subtraction.
template < class R >
template < class S >
inline
VectorBase<R>& VectorBase<R>::operator-=(const SSVectorBase<S>& vec)
{
   assert(dim() == vec.dim());

   if ( vec.isSetup() )
   {
      for( int i = vec.size() - 1; i >= 0; --i )
         val[vec.index(i)] -= vec.value(i);
   }
   else
   {
      for( int i = dim() - 1; i >= 0; --i )
         val[i] -= vec[i];
   }

   return *this;
}



/// Inner product.
template < class R >
inline
R VectorBase<R>::operator*(const SVectorBase<R>& vec) const
{
   assert(dim() >= vec.dim());

   R x(0);

   for( int i = vec.size() - 1; i >= 0; --i )
      x += val[vec.index(i)] * vec.value(i);

   return x;
}



/// Inner product.
template < class R >
inline
R VectorBase<R>::operator*(const SSVectorBase<R>& vec) const
{
   assert(dim() == vec.dim());

   if( vec.isSetup() )
   {
      const int* idx = vec.indexMem();

      R x(0);

      for( int i = vec.size() - 1; i >= 0; --i )
      {
         x += val[*idx] * vec.val[*idx];
         idx++;
      }

      return x;
   }
   else
      return operator*(static_cast<const VectorBase<R>&>(vec));
}



/// Addition of scaled vector.
template < class R >
template < class S, class T >
inline
VectorBase<R>& VectorBase<R>::multAdd(const S& x, const SVectorBase<T>& vec)
{
   for( int i = vec.size() - 1; i >= 0; --i )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)] += x * vec.value(i);
   }

   return *this;
}



/// Subtraction of scaled vector.
template < class R >
template < class S, class T >
inline
VectorBase<R>& VectorBase<R>::multSub(const S& x, const SVectorBase<T>& vec)
{
   for( int i = vec.size() - 1; i >= 0; --i )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= x * vec.value(i);
   }

   return *this;
}



#ifndef SOPLEX_LEGACY
/// Addition of scaled vector, specialization for rationals
template <>
template <>
inline
VectorBase<Rational>& VectorBase<Rational>::multAdd(const Rational& x, const SVectorBase<Rational>& vec)
{
   for( int i = vec.size() - 1; i >= 0; --i )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)].addProduct(x, vec.value(i));
   }

   return *this;
}



/// Subtraction of scaled vector, specialization for rationals
template <>
template <>
inline
VectorBase<Rational>& VectorBase<Rational>::multSub(const Rational& x, const SVectorBase<Rational>& vec)
{
   for( int i = vec.size() - 1; i >= 0; --i )
   {
      assert(vec.index(i) < dim());
      val[vec.index(i)].subProduct(x, vec.value(i));
   }

   return *this;
}
#endif



/// Addition of scaled vector.
template < class R >
template < class S, class T >
inline
VectorBase<R>& VectorBase<R>::multAdd(const S& x, const SSVectorBase<T>& vec)
{
   assert(vec.dim() <= dim());

   if( vec.isSetup() )
   {
      const int* idx = vec.indexMem();

      for( int i = vec.size() - 1; i>= 0; --i )
         val[idx[i]] += x * vec[idx[i]];
   }
   else
   {
      assert(vec.dim() == dim());

      for( int i = dim() - 1; i >= 0; --i )
         val[i] += x * vec.val[i];
   }

   return *this;
}



// ---------------------------------------------------------------------------------------------------------------------
//  Methods of DVectorBase
// ---------------------------------------------------------------------------------------------------------------------



/// Assignment operator.
template < class R >
template < class S >
inline
DVectorBase<R>& DVectorBase<R>::operator=(const SVectorBase<S>& vec)
{
   // dim() of SVector is not the actual dimension, rather the largest index plus 1
   // avoiding the reDim() saves unnecessary clearing of values
   if( vec.dim() > VectorBase<R>::dim() )
      reDim(vec.dim());

   VectorBase<R>::operator=(vec);

   assert(isConsistent());

   return *this;
}



// ---------------------------------------------------------------------------------------------------------------------
// Methods of SSVectorBase
// ---------------------------------------------------------------------------------------------------------------------



/// Addition.
template < class R >
template < class S >
inline
SSVectorBase<R>& SSVectorBase<R>::operator+=(const SVectorBase<S>& vec)
{
   VectorBase<R>::operator+=(vec);

   if( isSetup() )
   {
      setupStatus = false;
      setup();
   }

   return *this;
}



/// Subtraction.
template < class R >
template < class S >
inline
SSVectorBase<R>& SSVectorBase<R>::operator-=(const SVectorBase<S>& vec)
{
   VectorBase<R>::operator-=(vec);

   if( isSetup() )
   {
      setupStatus = false;
      setup();
   }

   return *this;
}



/// Addition of a scaled vector.
///@todo SSVectorBase::multAdd() should be rewritten without pointer arithmetic.
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::multAdd(S xx, const SVectorBase<T>& vec)
{
   if( isSetup() )
   {
      R* v = VectorBase<R>::val;
      R x;
      bool adjust = false;
      int j;

      for( int i = vec.size() - 1; i >= 0; --i )
      {
         j = vec.index(i);

         if( v[j] != 0 )
         {
            x = v[j] + xx * vec.value(i);
            if( isNotZero(x, epsilon) )
               v[j] = x;
            else
            {
               adjust = true;
               v[j] = SOPLEX_VECTOR_MARKER;
            }
         }
         else
         {
            x = xx * vec.value(i);
            if( isNotZero(x, epsilon) )
            {
               v[j] = x;
               addIdx(j);
            }
         }
      }

      if( adjust )
      {
         int* iptr = idx;
         int* iiptr = idx;
         int* endptr = idx + num;

         for( ; iptr < endptr; ++iptr )
         {
            x = v[*iptr];
            if( isNotZero(x, epsilon) )
               *iiptr++ = *iptr;
            else
               v[*iptr] = 0;
         }

         num = int(iiptr - idx);
      }
   }
   else
      VectorBase<R>::multAdd(xx, vec);

   assert(isConsistent());

   return *this;
}


/// Assigns pair wise vector product of setup x and setup y to SSVectorBase.
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::assignPWproduct4setup(const SSVectorBase<S>& x, const SSVectorBase<T>& y)
{
   assert(dim() == x.dim());
   assert(x.dim() == y.dim());
   assert(x.isSetup());
   assert(y.isSetup());

   clear();
   setupStatus = false;

   int i = 0;
   int j = 0;
   int n = x.size() - 1;
   int m = y.size() - 1;

   /* both x and y non-zero vectors? */
   if( m >= 0 && n >= 0 )
   {
      int xi = x.index(i);
      int yj = y.index(j);

      while( i < n && j < m )
      {
         if( xi == yj )
         {
            VectorBase<R>::val[xi] = R(x.val[xi]) * R(y.val[xi]);
            xi = x.index(++i);
            yj = y.index(++j);
         }
         else if( xi < yj )
            xi = x.index(++i);
         else
            yj = y.index(++j);
      }

      /* check (possible) remaining indices */

      while( i < n && xi != yj )
         xi = x.index(++i);

      while( j < m && xi != yj )
         yj = y.index(++j);

      if( xi == yj )
         VectorBase<R>::val[xi] = R(x.val[xi]) * R(y.val[xi]);
   }

   setup();

   assert(isConsistent());

   return *this;
}



/// Assigns \f$x^T \cdot A\f$ to SSVectorBase.
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::assign2product(const SSVectorBase<S>& x, const SVSetBase<T>& A)
{
   assert(A.num() == dim());

   R y;

   clear();

   for( int i = dim() - 1; i >= 0; --i )
   {
      y = A[i] * x;

      if( isNotZero(y, epsilon) )
      {
         VectorBase<R>::val[i] = y;
         IdxSet::addIdx(i);
      }
   }

   assert(isConsistent());

   return *this;
}



/// Assigns SSVectorBase to \f$A \cdot x\f$ for a setup \p x.
#define shortProductFactor 0.5
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::assign2product4setup(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(A.num() == x.dim());
   assert(x.isSetup());
   clear();

   if( x.size() == 1 )
   {
      assign2product1(A, x);
      setupStatus = true;
   }
   else if( isSetup() && (double(x.size()) * A.memSize() <= shortProductFactor * dim() * A.num()) )
   {
      assign2productShort(A, x);
      setupStatus = true;
   }
   else
   {
      assign2productFull(A, x);
      setupStatus = false;
   }

   assert(isConsistent());

   return *this;
}



/// Assignment helper.
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::assign2product1(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(x.isSetup());
   assert(x.size() == 1);

   // get the nonzero value of x and the corresponding vector in A:
   const int nzidx = x.idx[0];
   const T nzval = x.val[nzidx];
   const SVectorBase<S>& Ai = A[nzidx];

   // compute A[nzidx] * nzval:
   if( isZero(nzval, epsilon) || Ai.size() == 0 )
      clear();    // this := zero vector
   else
   {
      num = Ai.size();
      for( int j = num - 1; j >= 0; --j )
      {
         const Nonzero<S>& Aij = Ai.element(j);
         idx[j] = Aij.idx;
         VectorBase<R>::val[Aij.idx] = nzval * Aij.val;
      }
   }

   assert(isConsistent());

   return *this;
}



/// Assignment helper.
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::assign2productShort(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(x.isSetup());

   if( x.size() == 0 ) // x can be setup but have size 0 => this := zero vector
   {
      clear();
      return *this;
   }

   // compute x[0] * A[0]
   int curidx = x.idx[0];
   const T x0 = x.val[curidx];
   const SVectorBase<S>& A0 = A[curidx];
   int nonzero_idx = 0;
   int xsize = x.size();
   int Aisize;

   num = A0.size();
   if( isZero(x0, epsilon) || num == 0 )
   {
      // A[0] == 0 or x[0] == 0 => this := zero vector
      clear();
   }
   else
   {
      for( int j = 0; j < num; ++j )
      {
         const Nonzero<S>& elt = A0.element(j);
         const R product = x0 * elt.val;

         // store the value in any case
         idx[nonzero_idx] = elt.idx;
         VectorBase<R>::val[elt.idx] = product;

         // count only non-zero values; not 'isNotZero(product, epsilon)'
         if( product != 0 )
            ++nonzero_idx;
      }
   }

   // Compute the other x[i] * A[i] and add them to the existing vector.
   for( int i = 1; i < xsize; ++i )
   {
      curidx = x.idx[i];
      const T xi     = x.val[curidx];
      const SVectorBase<S>& Ai = A[curidx];

      // If A[i] == 0 or x[i] == 0, do nothing.
      Aisize = Ai.size();
      if ( isNotZero(xi, epsilon) || Aisize == 0 )
      {
         // Compute x[i] * A[i] and add it to the existing vector.
         for( int j = 0; j < Aisize; ++j )
         {
            const Nonzero<S>& elt = Ai.element(j);
            idx[nonzero_idx] = elt.idx;
            R oldval  = VectorBase<R>::val[elt.idx];

            // An old value of exactly 0 means the position is still unused.
            // It will be used now (either by a new nonzero or by a SOPLEX_VECTOR_MARKER),
            // so increase the counter. If oldval != 0, we just
            // change an existing NZ-element, so don't increase the counter.
            if( oldval == 0 )
               ++nonzero_idx;

            // Add the current product x[i] * A[i][j]; if oldval was
            // SOPLEX_VECTOR_MARKER before, it does not hurt because SOPLEX_VECTOR_MARKER is really small.
            oldval += xi * elt.val;

            // If the new value is exactly 0, mark the index as used
            // by setting a value which is nearly 0; otherwise, store
            // the value. Values below epsilon will be removed later.
            if( oldval == 0 )
               VectorBase<R>::val[elt.idx] = SOPLEX_VECTOR_MARKER;
            else
               VectorBase<R>::val[elt.idx] = oldval;
         }
      }
   }

   // Clean up by shifting all nonzeros (w.r.t. epsilon) to the front of idx,
   // zeroing all values which are nearly 0, and setting #num# appropriately.
   int nz_counter = 0;

   for( int i = 0; i < nonzero_idx; ++i )
   {
      curidx = idx[i];

      if( isZero( VectorBase<R>::val[curidx], epsilon ) )
         VectorBase<R>::val[curidx] = 0;
      else
      {
         idx[nz_counter] = curidx;
         ++nz_counter;
      }

      num = nz_counter;
   }

   assert(isConsistent());

   return *this;
}



/// Assignment helper.
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::assign2productFull(const SVSetBase<S>& A, const SSVectorBase<T>& x)
{
   assert(x.isSetup());

   if( x.size() == 0 ) // x can be setup but have size 0 => this := zero vector
   {
      clear();
      return *this;
   }

   bool A_is_zero = true;
   int xsize = x.size();
   int Aisize;

   for( int i = 0; i < xsize; ++i )
   {
      const int curidx = x.idx[i];
      const T xi = x.val[curidx];
      const SVectorBase<S>& Ai = A[curidx];
      Aisize = Ai.size();

      if( A_is_zero && Aisize > 0 )
         A_is_zero = false;

      for( int j = 0; j < Aisize; ++j )
      {
         const Nonzero<S>& elt = Ai.element(j);
         VectorBase<R>::val[elt.idx] += xi * elt.val;
      }
   }

   if( A_is_zero )
      clear(); // case x != 0 but A == 0

   return *this;
}



/// Assigns SSVectorBase to \f$A \cdot x\f$ thereby setting up \p x.
template < class R >
template < class S, class T >
inline
SSVectorBase<R>& SSVectorBase<R>::assign2productAndSetup(const SVSetBase<S>& A, SSVectorBase<T>& x)
{
   if( x.isSetup() )
      return assign2product4setup(A, x);

   if( x.dim() == 0 )
   { // x == 0 => this := zero vector
      clear();
      x.num = 0;
   }
   else
   {
      // x is not setup, so walk through its value vector
      int nzcount = 0;
      int end = x.dim();

      for( int i = 0; i < end; ++i )
      {
         // advance to the next element != 0
         T& xval = x.val[i];

         if( xval != 0 )
         {
            // If x[i] is really nonzero, compute A[i] * x[i] and adapt x.idx,
            // otherwise set x[i] to 0.
            if( isNotZero(xval, epsilon) )
            {
               const SVectorBase<S>& Ai = A[i];
               x.idx[ nzcount++ ] = i;

               for( int j = Ai.size() - 1; j >= 0; --j )
               {
                  const Nonzero<S>& elt = Ai.element(j);
                  VectorBase<R>::val[elt.idx] += xval * elt.val;
               }
            }
            else
               xval = 0;
         }
      }

      x.num = nzcount;
      setupStatus = false;
   }

   x.setupStatus = true;

   assert(isConsistent());

   return *this;
}



/// Assigns only the elements of \p rhs.
template < class R >
template < class S >
inline
SSVectorBase<R>& SSVectorBase<R>::assign(const SVectorBase<S>& rhs)
{
   assert(rhs.dim() <= VectorBase<R>::dim());

   int s = rhs.size();
   num = 0;

   for( int i = 0; i < s; ++i )
   {
      int k = rhs.index(i);
      S v = rhs.value(i);

      if( isZero(v, epsilon) )
         VectorBase<R>::val[k] = 0;
      else
      {
         VectorBase<R>::val[k] = v;
         idx[num++] = k;
      }
   }

   setupStatus = true;

   assert(isConsistent());

   return *this;
}



/// Assigns only the elements of \p rhs.
template <  >
template <  >
inline
SSVectorBase<Rational>& SSVectorBase<Rational>::assign(const SVectorBase<Rational>& rhs)
{
   assert(rhs.dim() <= VectorBase<Rational>::dim());

   int s = rhs.size();
   num = 0;

   for( int i = 0; i < s; ++i )
   {
      int k = rhs.index(i);
      const Rational& v = rhs.value(i);

      if( v == 0 )
         VectorBase<Rational>::val[k] = 0;
      else
      {
         VectorBase<Rational>::val[k] = v;
         idx[num++] = k;
      }
   }

   setupStatus = true;

   assert(isConsistent());

   return *this;
}



/// Assignment operator.
template < class R >
template < class S >
inline
SSVectorBase<R>& SSVectorBase<R>::operator=(const SVectorBase<S>& rhs)
{
   clear();

   return assign(rhs);
}



// ---------------------------------------------------------------------------------------------------------------------
//  Methods of SVectorBase
// ---------------------------------------------------------------------------------------------------------------------



/// Assignment operator.
template < class R >
template < class S >
inline
SVectorBase<R>& SVectorBase<R>::operator=(const VectorBase<S>& vec)
{
   int n = 0;
   Nonzero<R> *e = m_elem;

   clear();

   for( int i = vec.dim() - 1; i >= 0; --i )
   {
      if( vec[i] != 0 )
      {
         assert(n < max());

         e->idx = i;
         e->val = vec[i];
         ++e;
         ++n;
      }
   }

   set_size(n);

   return *this;
}



/// Assignment operator (specialization for Real).
template <>
template < class S >
inline
SVectorBase<Real>& SVectorBase<Real>::operator=(const VectorBase<S>& vec)
{
   int n = 0;
   Nonzero<Real> *e = m_elem;

   clear();

   for( int i = vec.dim() - 1; i >= 0; --i )
   {
      if( vec[i] != 0 )
      {
         assert(n < max());

         e->idx = i;
         e->val = Real(vec[i]);
         ++e;
         ++n;
      }
   }

   set_size(n);

   return *this;
}



/// Assignment operator.
template < class R >
template < class S >
inline
SVectorBase<R>& SVectorBase<R>::operator=(const SSVectorBase<S>& sv)
{
   assert(sv.isSetup());
   assert(max() >= sv.size());

   int nnz = 0;
   int idx;

   Nonzero<R> *e = m_elem;

   for( int i = 0; i < nnz; ++i )
   {
      idx = sv.index(i);
      if( sv.value(idx) != 0.0 )
      {
         e->idx = idx;
         e->val = sv[idx];
         ++e;
         ++nnz;
      }
   }
   set_size(nnz);

   return *this;
}



/// Inner product.
template < class R >
inline
R SVectorBase<R>::operator*(const VectorBase<R>& w) const
{
   R x = 0;
   Nonzero<R>* e = m_elem;

   for( int i = size() - 1; i >= 0; --i )
   {
      x += e->val * w[e->idx];
      e++;
   }

   return x;
}



// ---------------------------------------------------------------------------------------------------------------------
//  Methods of DSVectorBase
// ---------------------------------------------------------------------------------------------------------------------



/// Copy constructor.
template < class R >
template < class S >
inline
DSVectorBase<R>::DSVectorBase(const VectorBase<S>& vec)
   : theelem(0)
{
   allocMem((vec.dim() < 1) ? 2 : vec.dim());
   *this = vec;

   assert(isConsistent());
}



/// Copy constructor.
template < class R >
template < class S >
inline
DSVectorBase<R>::DSVectorBase(const SSVectorBase<S>& old)
   : theelem(0)
{
   allocMem(old.size() < 1 ? 2 : old.size());
   SVectorBase<R>::operator=(old);

   assert(isConsistent());
}



/// Assignment operator.
template < class R >
template < class S >
inline
DSVectorBase<R>& DSVectorBase<R>::operator=(const VectorBase<S>& vec)
{
   assert(this != (const DSVectorBase<R>*)(&vec));

   SVectorBase<R>::clear();
   setMax(vec.dim());
   SVectorBase<R>::operator=(vec);

   assert(isConsistent());

   return *this;
}



/// Assignment operator.
template < class R >
template < class S >
inline
DSVectorBase<R>& DSVectorBase<R>::operator=(const SSVectorBase<S>& vec)
{
   assert(this != &vec);

   SVectorBase<R>::clear();
   makeMem(vec.size());
   SVectorBase<R>::operator=(vec);

   return *this;
}



// ---------------------------------------------------------------------------------------------------------------------
//  Operators
// ---------------------------------------------------------------------------------------------------------------------



/// Output operator.
template < class R >
inline
std::ostream& operator<<(std::ostream& s, const VectorBase<R>& vec)
{
   int i;

   s << '(';

   for( i = 0; i < vec.dim() - 1; ++i )
      s << vec[i] << ", ";

   s << vec[i] << ')';

   return s;
}



/// Negation.
template < class R >
inline
DVectorBase<R> operator-(const VectorBase<R>& vec)
{
   DVectorBase<R> res(vec.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = -vec[i];

   return res;
}



/// Addition.
template < class R >
inline
DVectorBase<R> operator+(const VectorBase<R>& v, const VectorBase<R>& w)
{
   assert(v.dim() == w.dim());

   DVectorBase<R> res(v.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = v[i] + w[i];

   return res;
}



/// Addition.
template < class R >
inline
DVectorBase<R> operator+(const VectorBase<R>& v, const SVectorBase<R>& w)
{
   DVectorBase<R> res(v);

   res += w;

   return res;
}



/// Addition.
template < class R >
inline
DVectorBase<R> operator+(const SVectorBase<R>& v, const VectorBase<R>& w)
{
   return w + v;
}



/// Subtraction.
template < class R >
inline
DVectorBase<R> operator-(const VectorBase<R>& v, const VectorBase<R>& w)
{
   assert(v.dim() == w.dim());

   DVectorBase<R> res(v.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = v[i] - w[i];

   return res;
}



/// Subtraction.
template < class R >
inline
DVectorBase<R> operator-(const VectorBase<R>& v, const SVectorBase<R>& w)
{
   DVectorBase<R> res(v);

   res -= w;

   return res;
}



/// Subtraction.
template < class R >
inline
DVectorBase<R> operator-(const SVectorBase<R>& v, const VectorBase<R>& w)
{
   DVectorBase<R> res(w.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = -w[i];

   res += v;

   return res;
}



/// Scaling.
template < class R >
inline
DVectorBase<R> operator*(const VectorBase<R>& v, R x)
{
   DVectorBase<R> res(v.dim());

   for( int i = 0; i < res.dim(); ++i )
      res[i] = x * v[i];

   return res;
}



/// Scaling.
template < class R >
inline
DVectorBase<R> operator*(R x, const VectorBase<R>& v)
{
   return v * x;
}



/// Scaling.
template < class R >
inline
DSVectorBase<R> operator*(const SVectorBase<R>& v, R x)
{
   DSVectorBase<R> res(v.size());

   for( int i = 0; i < v.size(); ++i )
      res.add(v.index(i), v.value(i) * x);

   return res;
}



/// Scaling.
template < class R >
inline
DSVectorBase<R> operator*(R x, const SVectorBase<R>& v)
{
   return v * x;
}



/// Input operator.
template < class R >
inline
std::istream& operator>>(std::istream& s, DVectorBase<R>& vec)
{
   char c;
   R val;
   int i = 0;

   while( s.get(c).good() )
   {
      if( c != ' ' && c != '\t' && c != '\n' )
         break;
   }

   if( c != '(' )
      s.putback(c);
   else
   {
      do
      {
         s >> val;

         if( i >= vec.dim() - 1 )
            vec.reDim(i + 16);
         vec[i++] = val;

         while( s.get(c).good() )
         {
            if( c != ' ' && c != '\t' && c != '\n' )
               break;
         }

         if( c != ',' )
         {
            if (c != ')')
               s.putback(c);
            break;
         }
      }
      while( s.good() );
   }

   vec.reDim(i);

   return s;
}



/// Output operator.
template < class R >
inline
std::ostream& operator<<(std::ostream& os, const SVectorBase<R>& v)
{
   for( int i = 0, j = 0; i < v.size(); ++i )
   {
      if( j )
      {
         if( v.value(i) < 0 )
            os << " - " << -v.value(i);
         else
            os << " + " << v.value(i);
      }
      else
         os << v.value(i);

      os << " x" << v.index(i);
      j = 1;

      if( (i + 1) % 4 == 0 )
         os << "\n\t";
   }

   return os;
}



/// Output operator.
template < class R >
inline
std::ostream& operator<<(std::ostream& os, const SVSetBase<R>& s)
{
   for( int i = 0; i < s.num(); ++i )
      os << s[i] << "\n";

   return os;
}
}

/* reset the SOPLEX_DEBUG flag to its original value */
#undef SOPLEX_DEBUG
#ifdef SOPLEX_DEBUG_BASEVECTORS
#define SOPLEX_DEBUG
#undef SOPLEX_DEBUG_BASEVECTORS
#endif

#endif // _BASEVECTORS_H_
