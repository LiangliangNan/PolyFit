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


/**@file  ssvectorbase.h
 * @brief Semi sparse vector.
 */
#ifndef _SSVECTORBASE_H_
#define _SSVECTORBASE_H_

#include <assert.h>

#include "spxdefines.h"
#include "dvectorbase.h"
#include "idxset.h"
#include "spxalloc.h"

namespace soplex
{
template < class R > class SVectorBase;
template < class R > class SVSetBase;

/**@brief   Semi sparse vector.
 * @ingroup Algebra
 *
 *  This class implements semi-sparse vectors.  Such are #DVectorBase%s where the indices of its nonzero elements can be
 *  stored in an extra IdxSet.  Only elements with absolute value > #epsilon are considered to be nonzero.  Since really
 *  storing the nonzeros is not always convenient, an SSVectorBase provides two different stati: setup and not setup.
 *  An SSVectorBase being setup means that the nonzero indices are available, otherwise an SSVectorBase is just an
 *  ordinary DVectorBase with an empty IdxSet.  Note that due to arithmetic operation, zeros can slip in, i.e., it is
 *  only guaranteed that at least every non-zero is in the IdxSet.
 */
template < class R >
class SSVectorBase : public DVectorBase<R>, protected IdxSet
{
private:

   friend class DVectorBase<R>;
   friend class VectorBase<R>;
   template < class S > friend class DSVectorBase;

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   /// Is the SSVectorBase set up?
   bool setupStatus;

   /// A value x with |x| < epsilon is considered zero.
   R epsilon;

   /// Allocates enough space to accommodate \p newmax values.
   void setMax(int newmax)
   {
      assert(idx != 0);
      assert(newmax != 0);
      assert(newmax >= IdxSet::size());

      len = newmax;
      spx_realloc(idx, len);
   }

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Status of an SSVectorBase
    *
    *  An SSVectorBase can be set up or not. In case it is set up, its IdxSet correctly contains all indices of nonzero
    *  elements of the SSVectorBase.  Otherwise, it does not contain any useful data. Whether or not an SSVectorBase is
    *  setup can be determined with the method \ref soplex::SSVectorBase::isSetup() "isSetup()".
    *
    *  There are three methods for directly affecting the setup status of an SSVectorBase:
    *
    *  - unSetup():    This method sets the status to ``not setup''.
    *
    *  - setup():      This method initializes the IdxSet to the SSVectorBase's nonzero indices and sets the status to
    *                  ``setup''.
    *
    *  - forceSetup(): This method sets the status to ``setup'' without verifying that the IdxSet correctly contains all
    *                  nonzero indices. It may be used when the nonzero indices have been computed externally.
    */
   //@{

   /// Only used in slufactor.cpp.
   R* get_ptr()
   {
      return DVectorBase<R>::get_ptr();
   }
   /// Returns the non-zero epsilon used.
   R getEpsilon() const
   {
      return epsilon;
   }

   /// Changes the non-zero epsilon, invalidating the setup. */
   void setEpsilon(R eps)
   {
      if( eps != epsilon )
      {
         epsilon = eps;
         setupStatus = false;
      }
   }

   /// Returns setup status.
   bool isSetup() const
   {
      return setupStatus;
   }

   /// Makes SSVectorBase not setup.
   void unSetup()
   {
      setupStatus = false;
   }

   /// Initializes nonzero indices for elements with absolute values above #epsilon and sets all other elements to 0.
   void setup()
   {
      if( !isSetup() )
      {
         IdxSet::clear();

         int d = dim();
         num = 0;

         for( int i = 0; i < d; ++i )
         {
            if( VectorBase<R>::val[i] != R(0) )
            {
               if( spxAbs(VectorBase<R>::val[i]) <= epsilon )
                  VectorBase<R>::val[i] = R(0);
               else
               {
                  idx[num] = i;
                  num++;
               }
            }
         }

         setupStatus = true;

         assert(isConsistent());
      }
   }

   /// Forces setup status.
   void forceSetup()
   {
      setupStatus = true;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Methods for setup SSVectorBases */
   //@{

   /// Returns index of the \p n 'th nonzero element.
   int index(int n) const
   {
      assert(isSetup());

      return IdxSet::index(n);
   }

   /// Returns value of the \p n 'th nonzero element.
   R value(int n) const
   {
      assert(isSetup());
      assert(n >= 0 && n < size());

      return VectorBase<R>::val[idx[n]];
   }

   /// Finds the position of index \p i in the #IdxSet, or -1 if \p i doesn't exist.
   int pos(int i) const
   {
      assert(isSetup());

      return IdxSet::pos(i);
   }

   /// Returns the number of nonzeros.
   int size() const
   {
      assert(isSetup());

      return IdxSet::size();
   }

   /// Adds nonzero (\p i, \p x) to SSVectorBase.
   /** No nonzero with index \p i must exist in the SSVectorBase. */
   void add(int i, R x)
   {
      assert(VectorBase<R>::val[i] == R(0));
      assert(pos(i) < 0);

      addIdx(i);
      VectorBase<R>::val[i] = x;
   }

   /// Sets \p i 'th element to \p x.
   void setValue(int i, R x)
   {
      assert(i >= 0);
      assert(i < DVectorBase<R>::dim());

      if( isSetup() )
      {
         int n = pos(i);

         if( n < 0 )
         {
            if( spxAbs(x) > epsilon )
               IdxSet::add(1, &i);
         }
         else if( x == R(0) )
            clearNum(n);
      }

      VectorBase<R>::val[i] = x;

      assert(isConsistent());
   }

   /// Scale \p i 'th element by a
   void scaleValue(int i, int scaleExp)
   {
      assert(i >= 0);
      assert(i < DVectorBase<R>::dim());

      VectorBase<R>::val[i] = spxLdexp(VectorBase<R>::val[i], scaleExp);

      assert(isConsistent());
   }

   /// Clears element \p i.
   void clearIdx(int i)
   {
      if( isSetup() )
      {
         int n = pos(i);

         if( n >= 0 )
            remove(n);
      }

      VectorBase<R>::val[i] = 0;

      assert(isConsistent());
   }

   /// Sets \p n 'th nonzero element to 0 (index \p n must exist).
   void clearNum(int n)
   {
      assert(isSetup());
      assert(index(n) >= 0);

      VectorBase<R>::val[index(n)] = 0;
      remove(n);

      assert(isConsistent());
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Methods independent of the Status */
   //@{

   /// Returns \p i 'th value.
   R operator[](int i) const
   {
      return VectorBase<R>::val[i];
   }

   /// Returns array indices.
   const int* indexMem() const
   {
      return idx;
   }

   /// Returns array values.
   const R* values() const
   {
      return VectorBase<R>::val;
   }

   /// Returns indices.
   const IdxSet& indices() const
   {
      return *this;
   }

   /// Returns array indices.
   int* altIndexMem()
   {
      unSetup();
      return idx;
   }

   /// Returns array values.
   R* altValues()
   {
      unSetup();
      return VectorBase<R>::val;
   }

   /// Returns indices.
   IdxSet& altIndices()
   {
      unSetup();
      return *this;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Arithmetic operations */
   //@{

   /// Addition.
   template < class S >
   SSVectorBase<R>& operator+=(const VectorBase<S>& vec)
   {
      VectorBase<S>::operator+=(vec);

      if (isSetup())
      {
         setupStatus = false;
         setup();
      }
      return *this;
   }

   /// Addition.
   template < class S >
   SSVectorBase<R>& operator+=(const SVectorBase<S>& vec);

   /// Addition.
   template < class S >
   SSVectorBase<R>& operator+=(const SSVectorBase<S>& vec)
   {
      assert(vec.isSetup());
      for( int i = vec.size() - 1; i >= 0; --i )
         VectorBase<R>::val[vec.index(i)] += vec.value(i);

      if( isSetup() )
      {
         setupStatus = false;
         setup();
      }

      return *this;
   }

   /// Subtraction.
   template < class S >
   SSVectorBase<R>& operator-=(const VectorBase<S>& vec)
   {
      VectorBase<R>::operator-=(vec);

      if( isSetup() )
      {
         setupStatus = false;
         setup();
      }

      return *this;
   }

   /// Subtraction.
   template < class S >
   SSVectorBase<R>& operator-=(const SVectorBase<S>& vec);

   /// Subtraction.
   template < class S >
   SSVectorBase<R>& operator-=(const SSVectorBase<S>& vec)
   {
      if( vec.isSetup() )
      {
         for( int i = vec.size() - 1; i >= 0; --i )
            VectorBase<R>::val[vec.index(i)] -= vec.value(i);
      }
      else
         VectorBase<R>::operator-=(VectorBase<S>(vec));

      if( isSetup() )
      {
         setupStatus = false;
         setup();
      }

      return *this;
   }

   /// Scaling.
   template < class S >
   SSVectorBase<R>& operator*=(S x)
   {
      assert(isSetup());
      assert(x != S(0));

      for( int i = size() - 1; i >= 0; --i )
         VectorBase<R>::val[index(i)] *= x;

      assert(isConsistent());

      return *this;
   }

   // Inner product.
   template < class S >
   R operator*(const SSVectorBase<S>& w)
   {
      setup();

      R x = R(0);
      int i = size() - 1;
      int j = w.size() - 1;

      // both *this and w non-zero vectors?
      if( i >= 0 && j >= 0 )
      {
         int vi = index(i);
         int wj = w.index(j);

         while( i != 0 && j != 0 )
         {
            if( vi == wj )
            {
               x += VectorBase<R>::val[vi] * R(w.val[wj]);
               vi = index(--i);
               wj = w.index(--j);
            }
            else if( vi > wj )
               vi = index(--i);
            else
               wj = w.index(--j);
         }

         /* check remaining indices */

         while( i != 0 && vi != wj )
            vi = index(--i);

         while( j != 0 && vi != wj )
            wj = w.index(--j);

         if( vi == wj )
            x += VectorBase<R>::val[vi] * R(w.val[wj]);
      }

      return x;
   }

   /// Addition of a scaled vector.
   ///@todo SSVectorBase::multAdd() should be rewritten without pointer arithmetic.
   template < class S, class T >
   SSVectorBase<R>& multAdd(S xx, const SVectorBase<T>& vec);

   /// Addition of a scaled vector.
   template < class S, class T >
   SSVectorBase<R>& multAdd(S x, const VectorBase<T>& vec)
   {
      VectorBase<R>::multAdd(x, vec);

      if( isSetup() )
      {
         setupStatus = false;
         setup();
      }

      return *this;
   }

   /// Assigns pair wise vector product to SSVectorBase.
   template < class S, class T >
   SSVectorBase<R>& assignPWproduct4setup(const SSVectorBase<S>& x, const SSVectorBase<T>& y);

   /// Assigns \f$x^T \cdot A\f$ to SSVectorBase.
   template < class S, class T >
   SSVectorBase<R>& assign2product(const SSVectorBase<S>& x, const SVSetBase<T>& A);

   /// Assigns SSVectorBase to \f$A \cdot x\f$ for a setup \p x.
   template < class S, class T >
   SSVectorBase<R>& assign2product4setup(const SVSetBase<S>& A, const SSVectorBase<T>& x);

public:

   /// Assigns SSVectorBase to \f$A \cdot x\f$ thereby setting up \p x.
   template < class S, class T >
   SSVectorBase<R>& assign2productAndSetup(const SVSetBase<S>& A, SSVectorBase<T>& x);

   /// Maximum absolute value, i.e., infinity norm.
   R maxAbs() const
   {
      if( isSetup() )
      {
         R maxabs = 0;

         for( int i = 0; i < num; ++i )
         {
            R x = spxAbs(VectorBase<R>::val[idx[i]]);

            if( x > maxabs )
               maxabs = x;
         }

         return maxabs;
      }
      else
         return VectorBase<R>::maxAbs();
   }

   /// Squared euclidian norm.
   R length2() const
   {
      R x = 0;

      if( isSetup() )
      {
         for( int i = 0; i < num; ++i )
            x += VectorBase<R>::val[idx[i]] * VectorBase<R>::val[idx[i]];
      }
      else
         x = VectorBase<R>::length2();

      return x;
   }

   /// Floating point approximation of euclidian norm (without any approximation guarantee).
   Real length() const
   {
      return spxSqrt((Real)length2());
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Miscellaneous */
   //@{

   /// Dimension of VectorBase.
   int dim() const
   {
      return VectorBase<R>::dimen;
   }

   /// Resets dimension to \p newdim.
   void reDim(int newdim)
   {
      for( int i = IdxSet::size() - 1; i >= 0; --i )
      {
         if( index(i) >= newdim )
            remove(i);
      }

      DVectorBase<R>::reDim(newdim);
      setMax(DVectorBase<R>::memSize() + 1);

      assert(isConsistent());
   }

   /// Sets number of nonzeros (thereby unSetup SSVectorBase).
   void setSize(int n)
   {
      assert(n >= 0);
      assert(n <= IdxSet::max());

      unSetup();
      num = n;
   }

   /// Resets memory consumption to \p newsize.
   void reMem(int newsize)
   {
      DVectorBase<R>::reSize(newsize);
      assert(isConsistent());

      setMax(DVectorBase<R>::memSize() + 1);
   }

   /// Clears vector.
   void clear()
   {
      if( isSetup() )
      {
         for( int i = 0; i < num; ++i )
            VectorBase<R>::val[idx[i]] = 0;
      }
      else
         VectorBase<R>::clear();

      IdxSet::clear();
      setupStatus = true;

      assert(isConsistent());
   }

   /// consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( VectorBase<R>::dim() > IdxSet::max() )
         return MSGinconsistent("SSVectorBase");

      if( VectorBase<R>::dim() < IdxSet::dim() )
         return MSGinconsistent("SSVectorBase");

      if( isSetup() )
      {
         for( int i = 0; i < VectorBase<R>::dim(); ++i )
         {
            int j = pos(i);

            if( j < 0 && spxAbs(VectorBase<R>::val[i]) > 0 )
            {
               MSG_ERROR( std::cerr << "ESSVEC01 i = " << i
                  << "\tidx = " << j
                  << "\tval = " << std::setprecision(16) << VectorBase<R>::val[i]
                  << std::endl; )

               return MSGinconsistent("SSVectorBase");
            }
         }
      }

      return DVectorBase<R>::isConsistent() && IdxSet::isConsistent();
#else
      return true;
#endif
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Constructors / Destructors */
   //@{

   /// Default constructor.
   explicit SSVectorBase<R>(int p_dim, R p_eps = Param::epsilon())
      : DVectorBase<R>(p_dim)
      , IdxSet()
      , setupStatus(true)
      , epsilon(p_eps)
   {
      len = (p_dim < 1) ? 1 : p_dim;
      spx_alloc(idx, len);
      VectorBase<R>::clear();

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   SSVectorBase<R>(const SSVectorBase<S>& vec)
      : DVectorBase<R>(vec)
      , IdxSet()
      , setupStatus(vec.setupStatus)
      , epsilon(vec.epsilon)
   {
      len = (vec.dim() < 1) ? 1 : vec.dim();
      spx_alloc(idx, len);
      IdxSet::operator=(vec);

      assert(isConsistent());
   }

   /// Copy constructor.
   /** The redundancy with the copy constructor below is necessary since otherwise the compiler doesn't realize that it
    *  could use the more general one with S = R and generates a shallow copy constructor.
    */
   SSVectorBase<R>(const SSVectorBase<R>& vec)
      : DVectorBase<R>(vec)
      , IdxSet()
      , setupStatus(vec.setupStatus)
      , epsilon(vec.epsilon)
   {
      len = (vec.dim() < 1) ? 1 : vec.dim();
      spx_alloc(idx, len);
      IdxSet::operator=(vec);

      assert(isConsistent());
   }

   /// Constructs nonsetup copy of \p vec.
   template < class S >
   explicit SSVectorBase<R>(const VectorBase<S>& vec, R eps = Param::epsilon())
      : DVectorBase<R>(vec)
      , IdxSet()
      , setupStatus(false)
      , epsilon(eps)
   {
      len = (vec.dim() < 1) ? 1 : vec.dim();
      spx_alloc(idx, len);

      assert(isConsistent());
   }

   /// Sets up \p rhs vector, and assigns it.
   template < class S >
   void setup_and_assign(SSVectorBase<S>& rhs)
   {
      clear();
      epsilon = rhs.epsilon;
      setMax(rhs.max());
      DVectorBase<R>::reDim(rhs.dim());

      if( rhs.isSetup() )
      {
         IdxSet::operator=(rhs);

         for( int i = size() - 1; i >= 0; --i )
         {
            int j  = index(i);
            VectorBase<R>::val[j] = rhs.val[j];
         }
      }
      else
      {
         int d = rhs.dim();
         num = 0;

         for( int i = 0; i < d; ++i )
         {
            if( rhs.val[i] != 0 )
            {
               if( spxAbs(rhs.val[i]) > epsilon )
               {
                  rhs.idx[num] = i;
                  idx[num] = i;
                  VectorBase<R>::val[i] = rhs.val[i];
                  num++;
               }
               else
                  rhs.val[i] = 0;
            }
         }

         rhs.num = num;
         rhs.setupStatus = true;
      }

      setupStatus = true;

      assert(rhs.isConsistent());
      assert(isConsistent());
   }

   /// Assigns only the elements of \p rhs.
   template < class S >
   SSVectorBase<R>& assign(const SVectorBase<S>& rhs);

   /// Assignment operator.
   template < class S >
   SSVectorBase<R>& operator=(const SSVectorBase<S>& rhs)
   {
      assert(rhs.isConsistent());

      if( this != &rhs )
      {
         clear();
         epsilon = rhs.epsilon;
         setMax(rhs.max());
         DVectorBase<R>::reDim(rhs.dim());

         if( rhs.isSetup() )
         {
            IdxSet::operator=(rhs);

            for( int i = size() - 1; i >= 0; --i )
            {
               int j = index(i);
               VectorBase<R>::val[j] = rhs.val[j];
            }
         }
         else
         {
            int d = rhs.dim();
            num = 0;

            for( int i = 0; i < d; ++i )
            {
               if( spxAbs(rhs.val[i]) > epsilon )
               {
                  VectorBase<R>::val[i] = rhs.val[i];
                  idx[num] = i;
                  num++;
               }
            }
         }

         setupStatus = true;
      }

      assert(isConsistent());

      return *this;
   }

   /// Assignment operator.
   SSVectorBase<R>& operator=(const SSVectorBase<R>& rhs)
   {
      assert(rhs.isConsistent());

      if( this != &rhs )
      {
         clear();
         epsilon = rhs.epsilon;
         setMax(rhs.max());
         DVectorBase<R>::reDim(rhs.dim());

         if( rhs.isSetup() )
         {
            IdxSet::operator=(rhs);

            for( int i = size() - 1; i >= 0; --i )
            {
               int j = index(i);
               VectorBase<R>::val[j] = rhs.val[j];
            }
         }
         else
         {
            num = 0;

            for( int i = 0; i < rhs.dim(); ++i )
            {
               if( spxAbs(rhs.val[i]) > epsilon )
               {
                  VectorBase<R>::val[i] = rhs.val[i];
                  idx[num] = i;
                  num++;
               }
            }
         }

         setupStatus = true;
      }

      assert(isConsistent());

      return *this;
   }

   /// Assignment operator.
   template < class S >
   SSVectorBase<R>& operator=(const SVectorBase<S>& rhs);

   /// Assignment operator.
   template < class S >
   SSVectorBase<R>& operator=(const VectorBase<S>& rhs)
   {
      unSetup();
      VectorBase<R>::operator=(rhs);

      assert(isConsistent());

      return *this;
   }

   /// destructor
   ~SSVectorBase<R>()
   {
      if( idx )
         spx_free(idx);
   }

   //@}

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Private helpers */
   //@{

   /// Assignment helper.
   template < class S, class T >
   SSVectorBase<R>& assign2product1(const SVSetBase<S>& A, const SSVectorBase<T>& x);

   /// Assignment helper.
   template < class S, class T >
   SSVectorBase<R>& assign2productShort(const SVSetBase<S>& A, const SSVectorBase<T>& x);

   /// Assignment helper.
   template < class S, class T >
   SSVectorBase<R>& assign2productFull(const SVSetBase<S>& A, const SSVectorBase<T>& x);

   //@}
};

} // namespace soplex
#endif // _SSVECTORBASE_H_
