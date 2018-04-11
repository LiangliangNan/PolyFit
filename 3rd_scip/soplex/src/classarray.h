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

/**@file  classarray.h
 * @brief Save arrays of data objects.
 */
#ifndef _CLASSARRAY_H_
#define _CLASSARRAY_H_

#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <iostream>

#include "spxdefines.h"
#include "spxalloc.h"

namespace soplex
{
/**@brief   Safe arrays of class objects.
 * @ingroup Elementary
 *
 *  Class ClassArray provides safe arrays of general C++ objects (in contrast to data objects).  The elements of an
 *  instance of ClassArray can be accessed just like ordinary C++ array elements by means of the index
 *  operator[](). Safety is provided by
 *
 *  - automatic memory management in constructor and destructor preventing memory leaks
 *  - checking of array bounds when accessing elements with the indexing operator[]() when compiled without \c -DNDEBUG
 *
 *  Moreover, #ClassArray%s may easily be extended by #insert%ing or #append%ing elements to the ClassArray or shrunken
 *  by \ref remove() "removing" elements.  Method reSize(int n) resets the ClassArray%s length to \p n thereby possibly
 *  appending elements or truncating the ClassArray to the required size.
 *
 *  A ClassArray may be used as arguments for standard C functions requiring pointers through the use of get_ptr() and
 *  get_const_ptr().
 *
 *  Internally, a ClassArray object allocates a block of memory that fits up to max() elements, only size() of them are
 *  used. This makes extension and shrinking methods perform better.
 *
 *  @see Array, \ref DataObjects "Data Objects"
 */
template < class T >
class ClassArray
{
protected:
   int thesize;           ///< number of used elements in array data
   int themax;            ///< the length of array data and
   T*  data;              ///< the array of elements

protected:
   /** When a ClassArray is reSize()%d to more than max() elements, the new value for max() is not just set to the new
    *  size but rather to \p memFactor * \p size. This makes #reSize%ing perform better in codes where a ClassArray is
    *  extended often by a small number of elements only.
    */
   double memFactor;      ///< memory extension factor.

public:

   /// Reference to \p n 'th element.
   T& operator[](int n)
   {
      assert(n >= 0);
      assert(n < thesize);
      return data[n];
   }

   /// Reference to \p n 'th const element.
   const T& operator[](int n) const
   {
      assert(n >= 0);
      assert(n < thesize);
      return data[n];
   }

   /// Reference to last element.
   T& last()
   {
      assert(thesize > 0);
      return data[thesize-1];
   }

   /// Reference to last const element.
   const T& last() const
   {
      assert(thesize > 0);
      return data[thesize-1];
   }

   /// Gets a C pointer to the data.
   T* get_ptr()
   {
      return data;
   }

   /// Gets a const C pointer to the data.
   const T* get_const_ptr() const
   {
      return data;
   }

   /// Appends element \p t.
   void append(const T& t)
   {
      insert(thesize, 1, &t);
   }

   /// Appends \p n elements from \p t.
   void append(int n, const T t[])
   {
      insert(thesize, n, t);
   }

   /// Appends all elements from \p t.
   void append(const ClassArray<T>& t)
   {
      insert(thesize, t);
   }

   /// Inserts \p n uninitialized elements before \p i 'th element.
   void insert(int i, int n)
   {
      assert(n >= 0);
      assert(i >= 0);
      assert(i <= thesize);

      if( n > 0 )
      {
         int j = thesize;

         reSize(thesize + n);
         assert(thesize == j + n);

         /// move \p n elements in memory from insert position \p i to the back
         while( j > i )
         {
            j--;
            data[j + n] = data[j];
         }
      }
   }

   /// Inserts \p n elements from \p t before \p i 'the element.
   void insert(int i, int n, const T t[])
   {
      if( n > 0 )
      {
         insert(i, n);

         for( int j = 0; j < n; j++ )
            data[i + j] = t[j];
      }
   }

   /// Inserts all elements from \p t before \p i 'th element.
   void insert(int i, const ClassArray<T>& t)
   {
      if( t.size() )
      {
         insert(i, t.size());

         for( int j = 0; j < t.size(); j++ )
            data[i + j] = t[j];
      }
   }

   /// Removes \p m elements starting at \p n.
   void remove(int n = 0, int m = 1)
   {
      assert(n >= 0);
      assert(n < size());
      assert(m >= 0);
      assert(n + m <= size());

      for( int j = n + m; j < size(); j++ )
         data[j - m] = data[j];

      thesize -= m;
   }

   /// Removes \p m last elements.
   void removeLast(int m = 1)
   {
      assert(m >= 0);
      assert(m <= size());

      thesize -= m;
   }

   /// Removes all elements.
   void clear()
   {
      thesize = 0;
   }

   /// Returns number of elements.
   int size() const
   {
      return thesize;
   }

   /// Resets size to \p newsize.
   /** Resizing a ClassArray to less than the previous size, involves discarding its last elements. Resizing to a larger
    *  value involves adding uninitialized elements (similar to append()). If neccessary, also memory will be
    *  reallocated.
    */
   void reSize(int newsize)
   {
      assert(memFactor >= 1);

      if( newsize > themax )
         reMax(int(memFactor * newsize), newsize);
      else if( newsize < 0 )
         thesize = 0;
      else
         thesize = newsize;
   }

   /// Returns maximum number of elements.
   /** Even though the ClassArray currently holds no more than size() elements, up to max() elements could be added
    *  without need to reallocated free store.
    */
   int max() const
   {
      return themax;
   }

   /// Resets maximum number of elements.
   /** The value of max() is reset to \p newMax thereby setting size() to \p newSize. However, if \p newSize has a value
    *  \c < \c 0 (as the default argument does) size() remains unchanged and max() is set to MIN(size(), newMax). Hence,
    *  calling reMax() without the default arguments, will reduce the memory consumption to a minimum.  In no instance
    *  max() will be set to a value less than 1 (even if specified).
    *
    *  @return reMax returns the difference in bytes of the new and the old memory block, which can be used to update
    *          pointers pointing to elements of the memory block.
    */
   ptrdiff_t reMax(int newMax = 1, int newSize = -1)
   {
      /* process input */
      if( newSize < 0 )
         newSize = size();

      if( newMax < 1 )
         newMax = 1;

      if( newMax < newSize )
         newMax = newSize;

      /* nothing to reallocate */
      if( newMax == themax )
      {
         thesize = newSize;
         return 0;
      }

      /* allocate new memory */
      T* newMem = 0;
      spx_alloc(newMem, newMax);

      /* call copy constructor for first elements */
      int i;
      for( i = 0; i < size() && i < newSize; i++ )
         new (&(newMem[i])) T(data[i]);

      /* call default constructor for remaining elements */
      for( ; i < newMax; i++ )
         new (&(newMem[i])) T();

      /* compute pointer difference */
      ptrdiff_t pshift = reinterpret_cast<char*>(newMem) - reinterpret_cast<char*>(data);

      /* free old memory */
      for( i = themax-1; i >= 0; i-- )
         data[i].~T();

      spx_free(data);

      /* assign new memory */
      data = newMem;
      themax = newMax;
      thesize = newSize;

      return pshift;
   }

   /// Assignment operator.
   ClassArray& operator=(const ClassArray& rhs)
   {
      if( this != &rhs )
      {
         reSize(rhs.size());

         for( int i = 0; i < size(); i++ )
            data[i] = rhs.data[i];

         assert(isConsistent());
      }

      return *this;
   }

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( (data == 0)
         || (themax < 1)
         || (themax < thesize)
         || (thesize < 0)
         || (memFactor < 1.0) )
         return MSGinconsistent("ClassArray");
#endif
      return true;
   }

   /// Copy constructor.
   ClassArray(const ClassArray& old)
      : thesize(old.thesize)
      , themax(old.themax)
      , data(0)
      , memFactor(old.memFactor)
   {
      /* allocate memory */
      spx_alloc(data, max());

      /* call copy constructor for first elements */
      int i;
      for( i = 0; i < size(); i++ )
         new (&(data[i])) T(old.data[i]);

      /* call default constructor for remaining elements */
      for( ; i < max(); i++ )
         new (&(data[i])) T();

      assert(isConsistent());
   }

   /// Default constructor.
   /** The constructor allocates a ClassArray containing \p size uninitialized elements. The internal array is allocated
    *  to have \p max nonzeros, and the memory extension factor is set to \p fac.
    *
    *  @param p_size number of unitialised elements.
    *  @param p_max  maximum number of elements the array can hold.
    *  @param p_fac  value for memFactor.
    */
   explicit ClassArray(int p_size = 0, int p_max = 0, double p_fac = 1.2)
      : data(0)
      , memFactor(p_fac)
   {
      thesize = (p_size < 0) ? 0 : p_size;

      if( p_max > thesize )
         themax = p_max;
      else
         themax = (thesize == 0) ? 1 : thesize;

      spx_alloc(data, max());

      /* call default constructor for each element */
      for( int i = 0; i < max(); i++ )
         new (&(data[i])) T();

      assert(isConsistent());
   }

   /// Destructor.
   virtual ~ClassArray()
   {
      if( data )
      {
         for( int i = themax-1; i >= 0; i-- )
            data[i].~T();

         spx_free(data);
      }
   }
};

} // namespace soplex
#endif // _CLASSARRAY_H_
