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

/**@file  sorter.h
 * @brief Generic QuickSort implementation.
 */
#ifndef _SORTER_H_
#define _SORTER_H_

#include <assert.h>

namespace soplex
{

#define SHELLSORTMAX 25

/** shell-sort an array of data elements; use it only for arrays smaller than 25 entries */
template < class T, class COMPARATOR >
void SPxShellsort(T* keys, int end, COMPARATOR& compare, int start = 0)
{
   static const int incs[3] = {1, 5, 19}; /* sequence of increments */
   int k;

   assert(start <= end);

   for( k = 2; k >= 0; --k )
   {
      int h = incs[k];
      int first = h + start;
      int i;

      for( i = first; i <= end; ++i )
      {
         int j;
         T tempkey = keys[i];

         j = i;
         while( j >= first && compare(tempkey, keys[j-h]) < 0 )
         {
            keys[j] = keys[j-h];
            j -= h;
         }

         keys[j] = tempkey;
      }
   }
}



/// Generic QuickSort implementation.
/** This template function sorts an array \p t holding \p n elements of
    type T using \p compare for comparisons. Class COMPARATOR must provide
    an overloaded operator()(const T& t1,const T& t2) which returns
    - < 0, if \p t1 is to appear before \p t2,
    - = 0, if \p t1 and \p t2 can appear in any order, or
    - > 0, if \p t1 is to appear after \p t2.
*/

template < class T, class COMPARATOR >
void SPxQuicksort(T* keys, int end, COMPARATOR& compare, int start = 0, bool type = true)
{
   assert(start >= 0);

   /* nothing to sort for less than two elements */
   if( end <= start + 1 )
      return;

   /* reduce end position to last element index */
   --end;

   /* use quick-sort for long lists */
   while( end - start >= SHELLSORTMAX )
   {
      T pivotkey;
      T tmp;
      int lo;
      int hi;
      int mid;

      /* select pivot element */
      mid = (start+end)/2;
      pivotkey = keys[mid];

      /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
      lo = start;
      hi = end;
      for( ;; )
      {
         if( type )
         {
            while( lo < end && compare(keys[lo], pivotkey) < 0 )
               lo++;
            while( hi > start && compare(keys[hi], pivotkey) >= 0 )
               hi--;
         }
         else
         {
            while( lo < end && compare(keys[lo], pivotkey) <= 0 )
               lo++;
            while( hi > start && compare(keys[hi], pivotkey) > 0 )
               hi--;
         }

         if( lo >= hi )
            break;

         tmp = keys[lo];
         keys[lo] = keys[hi];
         keys[hi] = tmp;

         lo++;
         hi--;
      }
      assert((hi == lo-1) || (type && hi == start) || (!type && lo == end));

      /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
      if( type )
      {
         while( lo < end && compare(pivotkey, keys[lo]) >= 0 )
            lo++;

         /* make sure that we have at least one element in the smaller partition */
         if( lo == start )
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

            tmp = keys[lo];
            keys[lo] = keys[mid];
            keys[mid] = tmp;

            lo++;
         }
      }
      else
      {
         while( hi > start && compare(pivotkey, keys[hi]) <= 0 )
            hi--;

         /* make sure that we have at least one element in the smaller partition */
         if( hi == end )
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

            tmp = keys[hi];
            keys[hi] = keys[mid];
            keys[mid] = tmp;

            hi--;
         }
      }

      /* sort the smaller partition by a recursive call, sort the larger part without recursion */
      if( hi - start <= end - lo )
      {
         /* sort [start,hi] with a recursive call */
         if( start < hi )
         {
            SPxQuicksort(keys, hi+1, compare, start, !type);
         }

         /* now focus on the larger part [lo,end] */
         start = lo;
      }
      else
      {
         if( lo < end )
         {
            SPxQuicksort(keys, end+1, compare, lo, !type);
         }

         /* now focus on the larger part [start,hi] */
         end = hi;
      }
      type = !type;
   }

   /* use shell sort on the remaining small list */
   if( end - start >= 1 )
   {
      SPxShellsort(keys, end, compare, start);
   }


#ifdef CHECK_SORTING
   for( int i = start; i < end; ++i )
      assert(compare(keys[i], keys[i+1]) <= 0);
#endif
}


/**@brief  Generic implementation of Partial QuickSort.
 *
 * This template function sorts an array \p t holding \p n elements of
 * type T partially using \p compare for comparisons, i.e. ensures that
 * the \p size smallest elements are sorted to the front.
 *
 * Class COMPARATOR must provide an overloaded
 * operator()(const T& t1,const T& t2) which returns
 * - < 0, if \p t1 is to appear before \p t2,
 * - = 0, if \p t1 and \p t2 can appear in any order, or
 * - > 0, if \p t1 is to appear after \p t2.
 */
template < class T, class COMPARATOR >
int SPxQuicksortPart(
   T*                    keys,               /**< array of elements to be sorted between index start and end */
   COMPARATOR&           compare,            /**< comparator */
   int                   start,              /**< index of first element in range to be sorted */
   int                   end,                /**< index of last element in range to be sorted, plus 1 */
   int                   size,               /**< guaranteed number of sorted elements */
   int                   start2 = 0,         /**< auxiliary start index of sub range used for recursive call */
   int                   end2 = 0,           /**< auxiliary end index of sub range used for recursive call */
   bool                  type = true         /**< type of sorting, to be more flexable on degenerated cases */
   )
{
   assert(start >= 0);

   /* nothing to sort for less than two elements */
   if( end < start + 1 )
      return 0;
   else if( end == start + 1 )
      return 1;

   /* we assume that range {start, ..., start2-1} already contains the start2-start smallest elements in sorted order;
    * start2 has to lie in {start, ..., end-1} */
   if( start2 < start )
      start2 = start;

#ifdef CHECK_SORTING
   assert(start2 < end);

   for( int i = start; i < start2 - 1; ++i )
      assert(compare(keys[i], keys[i+1]) <= 0);
#endif
   assert(end2 <= end);

   /* if all remaining elements should be sorted, we simply call standard quicksort */
   if( start2 + size >= end - 1 )
   {
      SPxQuicksort(keys, end, compare, start2, type);
      return end-1;
   }

   T pivotkey;
   T tmp;
   int lo;
   int hi;
   int mid;

   /* reduce end position to last element index */
   --end;

   /* select pivot element */
   mid = (start2 + end)/2;
   pivotkey = keys[mid];

   /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
   lo = start2;
   hi = end;
   for( ;; )
   {
      if( type )
      {
         while( lo < end && compare(keys[lo], pivotkey) < 0 )
            lo++;
         while( hi > start2 && compare(keys[hi], pivotkey) >= 0 )
            hi--;
      }
      else
      {
         while( lo < end && compare(keys[lo], pivotkey) <= 0 )
            lo++;
         while( hi > start2 && compare(keys[hi], pivotkey) > 0 )
            hi--;
      }

      if( lo >= hi )
         break;

      tmp = keys[lo];
      keys[lo] = keys[hi];
      keys[hi] = tmp;

      lo++;
      hi--;
   }
   assert((hi == lo-1) || (type && hi == start2) || (!type && lo == end));

   /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
   if( type )
   {
      while( lo < end && compare(pivotkey, keys[lo]) >= 0 )
         lo++;

      /* make sure that we have at least one element in the smaller partition */
      if( lo == start2 )
      {
         /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
         assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

         tmp = keys[lo];
         keys[lo] = keys[mid];
         keys[mid] = tmp;

         lo++;
      }
   }
   else
   {
      while( hi > start2 && compare(pivotkey, keys[hi]) <= 0 )
         hi--;

      /* make sure that we have at least one element in the smaller partition */
      if( hi == end )
      {
         /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
         assert(compare(keys[mid], pivotkey) == 0); /* the pivot element did not change its position */

         tmp = keys[hi];
         keys[hi] = keys[mid];
         keys[mid] = tmp;

         hi--;
      }
   }

#ifdef CHECK_SORTING
   for( int i = start2; i < lo; ++i )
      assert(compare(keys[i], pivotkey) <= 0);
#endif

   /* if we only need to sort less than half of the "<" part, use partial sort again */
   if( 2*size <= hi - start2 )
   {
      return SPxQuicksortPart(keys, compare, start, hi+1, size, start2, end2, !type);
   }
   /* otherwise, and if we do not need to sort the ">" part, use standard quicksort on the "<" part */
   else if( size <= lo - start2 )
   {
      SPxQuicksort(keys, hi+1, compare, start2, !type);
      return lo-1;
   }
   /* otherwise we have to sort the "<" part fully (use standard quicksort) and the ">" part partially */
   else
   {
      SPxQuicksort(keys, hi+1, compare, start2, !type);
      return SPxQuicksortPart(keys, compare, start, end, size+start2-lo, lo, end2, !type);
   }
}

} // namespace soplex
#endif // _SORTER_H_
