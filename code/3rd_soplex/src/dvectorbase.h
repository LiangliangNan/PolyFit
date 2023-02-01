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

/**@file  dvectorbase.h
 * @brief Dynamic dense vectors.
 */
#ifndef _DVECTORBASE_H_
#define _DVECTORBASE_H_

#include <iostream>
#include <assert.h>

#include "spxdefines.h"
#include "spxalloc.h"
#include "vectorbase.h"

namespace soplex
{
template < class R > class SVectorBase;

/**@brief   Dynamic dense vectors.
 * @ingroup Algebra
 *
 *  Class DVectorBase is a derived class of VectorBase adding automatic memory management to such objects.  This allows
 *  to implement maths operations operator+() and operator-().  Further, it is possible to reset the dimension of a
 *  DVectorBase via method reDim().  However, this may render all references to values of a #reDim()%ed DVectorBase
 *  invalid.
 *
 *  For vectors that are often subject to reDim() it may be unconvenient to reallocate the required memory every time.
 *  Instead, an array of values of length memSize() is kept, where only the first dim() elements are used.  Initially,
 *  memSize() == dim().  However, if the dimension is increased, memSize() will be increased at least by a factor of 1.2
 *  to be prepared for future (small) #reDim()%s.  Finally, one can explicitly set memSize() with method reSize(), but
 *  not lower than dim().
 */
template < class R >
class DVectorBase : public VectorBase<R>
{
   template < class S > friend class DVectorBase;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   /// Length of array of values \ref soplex::DVectorBase::mem "mem"
   int memsize;

   /// Array of values.
   R* mem;

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction, destruction, and assignment */
   //@{

   /// Default constructor. \p d is the initial dimension.
   explicit DVectorBase<R>(int d = 0)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      memsize = (d > 0) ? d : 4;

      spx_alloc(mem, memsize);
      for( int i = 0; i < memsize; i++ )
         new (&(mem[i])) R();

      VectorBase<R>::val = mem;
      VectorBase<R>::dimen = d;

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   explicit DVectorBase<R>(const VectorBase<S>& old)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      VectorBase<R>::dimen = old.dim();
      memsize = VectorBase<R>::dimen;

      spx_alloc(mem, memsize);

      for( int i = 0; i < VectorBase<R>::dimen; i++ )
         new (&(mem[i])) R(old[i]);
      for( int i = VectorBase<R>::dimen; i < memsize; i++ )
         new (&(mem[i])) R();

      VectorBase<R>::val = mem;

      assert(isConsistent());
   }

   /// Copy constructor.
   /** The redundancy with the copy constructor below is necessary since otherwise the compiler doesn't realize that it
    *  could use the more general one with S = R and generates a shallow copy constructor.
    */
   DVectorBase<R>(const DVectorBase<R>& old)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      VectorBase<R>::dimen = old.dim();
      memsize = old.memsize;

      spx_alloc(mem, memsize);

      for( int i = 0; i < VectorBase<R>::dimen; i++ )
         new (&(mem[i])) R(old[i]);
      for( int i = VectorBase<R>::dimen; i < memsize; i++ )
         new (&(mem[i])) R();

      VectorBase<R>::val = mem;

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   DVectorBase<R>(const DVectorBase<S>& old)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      VectorBase<R>::dimen = old.dim();
      memsize = old.memsize;

      spx_alloc(mem, memsize);

      for( int i = 0; i < VectorBase<R>::dimen; i++ )
         new (&(mem[i])) R(old[i]);
      for( int i = VectorBase<R>::dimen; i < memsize; i++ )
         new (&(mem[i])) R();

      VectorBase<R>::val = mem;

      assert(isConsistent());
   }

   /// Assignment operator.
   DVectorBase<R>& operator=(const VectorBase<R>& vec)
   {
      if( (VectorBase<R>*)this != &vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DVectorBase<R>& operator=(const VectorBase<S>& vec)
   {
      if( (VectorBase<S>*)this != &vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   DVectorBase<R>& operator=(const DVectorBase<R>& vec)
   {
      if( this != &vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DVectorBase<R>& operator=(const DVectorBase<S>& vec)
   {
      if( this != (const DVectorBase<R>*)&vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DVectorBase<R>& operator=(const SVectorBase<S>& vec);

   /// Destructor.
   virtual ~DVectorBase<R>()
   {
      if( mem != 0 )
      {
         for( int i = memsize-1; i >= 0; i-- )
            mem[i].~R();

         spx_free(mem);
      }
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Access and modification */
   //@{

   /// Returns \ref soplex::DVectorBase "DVectorBase"'s memory size.
   int memSize() const
   {
      return memsize;
   }

   /// Resets \ref soplex::DVectorBase "DVectorBase"'s dimension to \p newdim.
   void reDim(int newdim, const bool setZero = true )
   {
      assert(memsize >= 0);

      if( newdim > memsize )
         reSize(int(newdim + 0.2 * memsize));

      if( setZero )
      {
         for( int i = VectorBase<R>::dimen; i < newdim; i++ )
            mem[i] = 0;
      }

      VectorBase<R>::dimen = newdim;
   }

   /// Resets \ref soplex::DVectorBase "DVectorBase"'s memory size to \p newsize.
   void reSize(int newsize)
   {
      assert(newsize > VectorBase<R>::dim());

      /* allocate new memory */
      R* newmem = 0;
      spx_alloc(newmem, newsize);

      /* call copy constructor for first elements */
      int i;
      for( i = 0; i < VectorBase<R>::dim(); i++ )
         new (&(newmem[i])) R(mem[i]);

      /* call default constructor for remaining elements */
      for( ; i < newsize; i++ )
         new (&(newmem[i])) R();

      /* free old memory */
      for( i = memsize-1; i >= 0; i-- )
         mem[i].~R();

      spx_free(mem);

      /* assign new memory */
      mem = newmem;
      memsize = newsize;
      VectorBase<R>::val = mem;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Utilities */
   //@{

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( VectorBase<R>::val != mem || VectorBase<R>::dimen > memsize || VectorBase<R>::dimen < 0 )
         return MSGinconsistent("DVectorBase");

      return VectorBase<R>::isConsistent();
#else
      return true;
#endif
   }

   //@}
};



/// Resets \ref soplex::DVectorBase "DVectorBase"'s memory size to \p newsize (specialization for Real).
template<>
inline
void DVectorBase<Real>::reSize(int newsize)
{
   assert(newsize >= dim());

   spx_realloc(mem, newsize);

   val = mem;
   memsize = newsize;
}



/// Default constructor with \p d as the initial dimension (specialization for Real).
template<>
inline
DVectorBase<Real>::DVectorBase(int d)
   : VectorBase<Real>(0, 0)
   , mem(0)
{
   memsize = (d > 0) ? d : 4;
   spx_alloc(mem, memsize);
   val = mem;
   dimen = d;

   assert(DVectorBase<Real>::isConsistent());
}



/// Copy constructor (specialization for Real).
template<>
template<>
inline
DVectorBase<Real>::DVectorBase(const VectorBase<Real>& old)
   : VectorBase<Real>(0, 0)
   , mem( 0 )
{
   dimen = old.dim();
   memsize = dimen;
   spx_alloc(mem, memsize);
   val = mem;
   *this = old;

   assert(DVectorBase<Real>::isConsistent());
}



/// Copy constructor (specialization for Real).
template<>
inline
DVectorBase<Real>::DVectorBase(const DVectorBase<Real>& old)
   : VectorBase<Real>(0, 0)
   , mem(0)
{
   dimen = old.dim();
   memsize = old.memsize;
   spx_alloc(mem, memsize);
   val = mem;
   *this = old;

   assert(DVectorBase<Real>::isConsistent());
}
} // namespace soplex
#endif // _DVECTORBASE_H_
