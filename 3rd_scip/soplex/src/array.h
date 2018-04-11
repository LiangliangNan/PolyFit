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

/**@file  array.h
 * @brief Save arrays of arbitrary types.
 */
#ifndef _ARRAY_H_
#define _ARRAY_H_

#include <assert.h>
#include <string.h>
#include "spxalloc.h"

namespace soplex
{
/**@brief   Safe arrays of arbitrary types.
   @ingroup Elementary

   Class Array provides safe arrays of arbitrary type. Array elements are
   accessed just like ordinary C++ array elements by means of the index
   operator[](). Safety is provided by

    - automatic memory management in constructor and destructure
      preventing memory leaks
    - checking of array bound when accessing elements with the
      indexing operator[]() (only when compiled without \c -DNDEBUG).
 
    Moreover, #Array%s may easily be extended by #insert%ing or
    #append%ing elements to the Array or shrunken by 
    \ref remove() "removing"
    elements. Method reSize(int n) resets the Array's length to \p n,
    thereby appending elements or truncating the Array to the
    required size.
 
    An Array is implemented in a C++-compliant way with respect to
    how memory is managed: Only operators new and delete are
    used for allocating memory. This involves some overhead for all
    methods effecting the length of an Array, i.e., all methods
    insert(), append(), remove() and reSize(). This involves
    allocating a new C++ array of the new size and copying all
    elements with the template parameters operator=().
 
    For this reason, it is not convenient to use class Array if its elements
    are \ref DataObjects "Data Objects". In this case use class DataArray 
    instead.
 
    @see DataArray, \ref DataObjects "Data Objects" 
*/
template < class T >
class Array
{
protected:

   //----------------------------------------
   /**@name Data */
   //@{
   int num;     ///< the length of array data 
   T*  data;    ///< the array of elements
   //@}

public:

   //----------------------------------------
   /**@name Access / modification */
   //@{
   /// reference \p n 'th element.
   T& operator[](int n)
   {
      assert(n >= 0 && n < num);
      return data[n];
   }
   /// reference \p n 'th element.
   const T& operator[](int n) const
   {
      assert(n >= 0 && n < num);
      return data[n];
   }

   /** This function serves for using a Vector in an C-style
    *  function. It returns a pointer to the first value of the array.
    */
   T* get_ptr()
   {
      return data;
   }

   /// append \p n uninitialized elements.
   void append(int n)
   {
      insert(num, n);
   }
   /// append \p n elements from \p p_array.
   void append(int n, const T* p_array)
   {
      insert(num, n, p_array);
   }
   /// append all elements from \p p_array.
   void append(const Array<T>& p_array)
   {
      insert(num, p_array);
   }

   /** insert \p n uninitialized elements before \p i 'th element.
    *
    * You must not! use realloc, memcpy or memmove, because some data element points inside itself, and therefore you
    * always need to copy all elements by hand.
    */
   void insert(int i, int n)
   {
      assert(i <= num);
      assert(num >= 0);
      if (n > 0)
      {
         T* newdata = 0;
         int k;

         spx_alloc(newdata, num + n);
         assert(newdata != 0);

         // copy front segment to new array
         for( k = 0; k < i; ++k )
         {
            new (&(newdata[k])) T();
            newdata[k] = data[k];
            data[k].~T();
         }

         // call constructor for new elements
         for( ; k < i+n; ++k )
            new (&(newdata[k])) T();

         // copy rear segment to new array
         for( k = i; k < num; ++k )
         {
            new (&(newdata[n + k])) T();
            newdata[n + k] = data[k];
            data[k].~T();
         }

         if( data )
            spx_free(data);
         data = newdata;
         num += n;
      }
   }

   /// insert \p n elements from \p p_array before \p i 'th element.
   void insert(int i, int n, const T* p_array)
   {
      insert(i, n);
      for (n--; n >= 0; --n)
         data[n + i] = p_array[n];
   }

   /// insert all elements from \p p_array before \p i 'th element.
   void insert(int i, const Array<T>& p_array)
   {
      int n = p_array.size();
      insert(i, n);
      for (n--; n >= 0; --n)
         data[n + i] = p_array.data[n];
   }

   /** remove \p m elements starting at \p n.
    *
    * You must not! use realloc, memcpy or memmove, because some data element points inside itself, and therefore you
    * always need to copy all elements by hand.
    */
   void remove(int n = 0, int m = 1)
   {
      assert(n >= 0 && m >= 0);
      if (m > 0 && n < num)
      {
         T* newdata = 0;
         int k;

         assert(num == size());
         m -= (n + m <= num) ? 0 : n + m - num;

         if( num > m )
         {
            spx_alloc(newdata, num - m);
            assert(newdata != 0);

            // copy rear segment to new array
            for( k = num - 1; k >= n + m; --k )
            {
               new (&(newdata[k - m])) T();
               newdata[k - m] = data[k];
               data[k].~T();
            }

            // call destructor for old elements
            for( ; k >= n; --k )
               data[k].~T();

            // copy front segment to new array
            for( ; k >= 0; --k )
            {
               new (&(newdata[k])) T();
               newdata[k] = data[k];
               data[k].~T();
            }
         }
         else
         {
            assert(num == m);
            assert(n == 0);

            // call destructor for old elements
            for( k = m - 1; k >= 0; --k )
               data[k].~T();
         }

         assert(data != 0);
         spx_free(data);
         data = newdata;

         num -= m;
         assert((data == 0) == (num == 0));
      }
   }

   /// remove all elements.
   void clear()
   {
      // call destructors of all elements
      while( num > 0 )
      {
         --num;
         data[num].~T();
      }
      if( data )
         spx_free(data);
      assert(num == 0);
   }

   /// return the number of elements.
   int size() const
   {
      return num;
   }

   /// reset the number of elements.
   void reSize(int newsize)
   {
      if (newsize < num)
         remove(newsize, num - newsize);
      else if (newsize > num)
         append(newsize - num);
   }
   //@}

   //----------------------------------------
   /**@name Construction / destruction */
   //@{
   /// assignment operator.
   /** Assigning an rvalue Array to an lvalue Array involves resizing
    *  the lvalue to the rvalues size() and copying all elements via
    *  the Array element's assignment operator=().
    */
   Array<T>& operator=(const Array<T>& rhs)
   {
      if (this != &rhs)
      {
         reSize(rhs.size());
         for (int i = 0; i < num; ++i)
            data[i] = rhs.data[i];
         assert(Array::isConsistent());
      }
      return *this;
   }

   /// default constructor.
   /** The constructor allocates an Array of \p n uninitialized elements.
    */
   explicit
   Array(int n = 0) 
      : data(0)
   {
      assert(n >= 0);
      num = n;
      if (num > 0)
      {
         int k;

         spx_alloc(data, num);
         assert(data != 0);

         for( k = 0; k < num; ++k )
            new (&(data[k])) T();
      }
      assert(Array::isConsistent());
   }

   /// copy constructor
   Array(const Array<T>& old) 
      : num(old.num)
   {
      if (num > 0)
      {
         int k;

         data = 0;
         spx_alloc(data, num);
         assert(data != 0);

         for( k = 0; k < num; ++k )
         {
            new (&(data[k])) T();
            data[k] = old.data[k];
         }
      }
      else
         data = 0;
      assert(Array::isConsistent());
   }

   /// destructor
   ~Array()
   {
      while( num > 0 )
      {
         --num;
         data[num].~T();
      }
      if( data )
         spx_free(data);
   }

   /// consistency check
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if (num < 0 || (num > 0 && data == 0))
         return MSGinconsistent("Array");
#endif

      return true;
   }
   //@}
};
} // namespace soplex
#endif // _ARRAY_H_
