/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sorttpl.c
 * @brief  template functions for sorting
 * @author Michael Winkler
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* template parameters that have to be passed in as #define's:
 * #define SORTTPL_NAMEEXT      <ext>      extension to be used for SCIP method names, for example DownIntRealPtr
 * #define SORTTPL_KEYTYPE      <type>     data type of the key array
 * #define SORTTPL_FIELD1TYPE   <type>     data type of first additional array which should be sorted in the same way (optional)
 * #define SORTTPL_FIELD2TYPE   <type>     data type of second additional array which should be sorted in the same way (optional)
 * #define SORTTPL_FIELD3TYPE   <type>     data type of third additional array which should be sorted in the same way (optional)
 * #define SORTTPL_FIELD4TYPE   <type>     data type of fourth additional array which should be sorted in the same way (optional)
 * #define SORTTPL_FIELD5TYPE   <type>     data type of fifth additional array which should be sorted in the same way (optional)
 * #define SORTTPL_FIELD6TYPE   <type>     data type of fifth additional array which should be sorted in the same way (optional)
 * #define SORTTPL_PTRCOMP                 ptrcomp method should be used for comparisons (optional)
 * #define SORTTPL_INDCOMP                 indcomp method should be used for comparisons (optional)
 * #define SORTTPL_BACKWARDS               should the array be sorted other way around
 */
#include "scip/def.h"
#define SORTTPL_SHELLSORTMAX 25
#define SORTTPL_MINSIZENINTHER 729 /* minimum input size to use ninther (median of nine) for pivot selection */

#ifndef SORTTPL_NAMEEXT
#error You need to define SORTTPL_NAMEEXT.
#endif
#ifndef SORTTPL_KEYTYPE
#error You need to define SORTTPL_KEYTYPE.
#endif

#ifdef SORTTPL_EXPANDNAME
#undef SORTTPL_EXPANDNAME
#endif
#ifdef SORTTPL_NAME
#undef SORTTPL_NAME
#endif

/* enabling and disabling additional lines in the code */
#ifdef SORTTPL_FIELD1TYPE
#define SORTTPL_HASFIELD1(x)    x
#define SORTTPL_HASFIELD1PAR(x) x,
#else
#define SORTTPL_HASFIELD1(x)    /**/
#define SORTTPL_HASFIELD1PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD2TYPE
#define SORTTPL_HASFIELD2(x)    x
#define SORTTPL_HASFIELD2PAR(x) x,
#else
#define SORTTPL_HASFIELD2(x)    /**/
#define SORTTPL_HASFIELD2PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD3TYPE
#define SORTTPL_HASFIELD3(x)    x
#define SORTTPL_HASFIELD3PAR(x) x,
#else
#define SORTTPL_HASFIELD3(x)    /**/
#define SORTTPL_HASFIELD3PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD4TYPE
#define SORTTPL_HASFIELD4(x)    x
#define SORTTPL_HASFIELD4PAR(x) x,
#else
#define SORTTPL_HASFIELD4(x)    /**/
#define SORTTPL_HASFIELD4PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD5TYPE
#define SORTTPL_HASFIELD5(x)    x
#define SORTTPL_HASFIELD5PAR(x) x,
#else
#define SORTTPL_HASFIELD5(x)    /**/
#define SORTTPL_HASFIELD5PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD6TYPE
#define SORTTPL_HASFIELD6(x)    x
#define SORTTPL_HASFIELD6PAR(x) x,
#else
#define SORTTPL_HASFIELD6(x)    /**/
#define SORTTPL_HASFIELD6PAR(x) /**/
#endif
#ifdef SORTTPL_PTRCOMP
#define SORTTPL_HASPTRCOMP(x)    x
#define SORTTPL_HASPTRCOMPPAR(x) x,
#else
#define SORTTPL_HASPTRCOMP(x)    /**/
#define SORTTPL_HASPTRCOMPPAR(x) /**/
#endif
#ifdef SORTTPL_INDCOMP
#define SORTTPL_HASINDCOMP(x)    x
#define SORTTPL_HASINDCOMPPAR(x) x,
#else
#define SORTTPL_HASINDCOMP(x)    /**/
#define SORTTPL_HASINDCOMPPAR(x) /**/
#endif


/* the two-step macro definition is needed, such that macro arguments
 * get expanded by prescan of the C preprocessor (see "info cpp",
 * chapter 3.10.6: Argument Prescan)
 */
#define SORTTPL_EXPANDNAME(method, methodname) \
   method ## methodname
#define SORTTPL_NAME(method, methodname) \
  SORTTPL_EXPANDNAME(method, methodname)

/* comparator method */
#ifdef SORTTPL_PTRCOMP
#ifdef SORTTPL_BACKWARDS
#define SORTTPL_CMP(x,y) (-ptrcomp((x), (y)))
#else
#define SORTTPL_CMP(x,y) (ptrcomp((x), (y)))
#endif
#else
#ifdef SORTTPL_INDCOMP
#ifdef SORTTPL_BACKWARDS
#define SORTTPL_CMP(x,y) (-indcomp(dataptr, (x), (y)))
#else
#define SORTTPL_CMP(x,y) (indcomp(dataptr, (x), (y)))
#endif
#else
#ifdef SORTTPL_BACKWARDS
#define SORTTPL_CMP(x,y) ((y) - (x))
#else
#define SORTTPL_CMP(x,y) ((x) - (y))
#endif
#endif
#endif

#define SORTTPL_ISBETTER(x,y) (SORTTPL_CMP(x,y) < 0)
#define SORTTPL_ISWORSE(x,y) (SORTTPL_CMP(x,y) > 0)

/* swapping two variables */
#define SORTTPL_SWAP(T,x,y) \
   {                \
      T temp = x;   \
      x = y;        \
      y = temp;     \
   }


/** shell-sort an array of data elements; use it only for arrays smaller than 25 entries */
static
void SORTTPL_NAME(sorttpl_shellSort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SCIP_Real*            weights,            /**< real, nonnegative weights that should be permuted like key, or NULL */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   start,              /**< starting index */
   int                   end                 /**< ending index */
   )
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
         SORTTPL_KEYTYPE tempkey = key[i];

         SCIP_Real tmpweight = weights != NULL ? weights[i] : 1;

         SORTTPL_HASFIELD1( SORTTPL_FIELD1TYPE tempfield1 = field1[i]; )
         SORTTPL_HASFIELD2( SORTTPL_FIELD2TYPE tempfield2 = field2[i]; )
         SORTTPL_HASFIELD3( SORTTPL_FIELD3TYPE tempfield3 = field3[i]; )
         SORTTPL_HASFIELD4( SORTTPL_FIELD4TYPE tempfield4 = field4[i]; )
         SORTTPL_HASFIELD5( SORTTPL_FIELD5TYPE tempfield5 = field5[i]; )
         SORTTPL_HASFIELD6( SORTTPL_FIELD6TYPE tempfield6 = field6[i]; )

         j = i;
         while( j >= first && SORTTPL_ISBETTER(tempkey, key[j-h]) )
         {
            key[j] = key[j-h];

            if( weights != NULL )
               weights[j] = weights[j - h];

            SORTTPL_HASFIELD1( field1[j] = field1[j-h]; )
            SORTTPL_HASFIELD2( field2[j] = field2[j-h]; )
            SORTTPL_HASFIELD3( field3[j] = field3[j-h]; )
            SORTTPL_HASFIELD4( field4[j] = field4[j-h]; )
            SORTTPL_HASFIELD5( field5[j] = field5[j-h]; )
            SORTTPL_HASFIELD6( field6[j] = field6[j-h]; )
            j -= h;
         }

         key[j] = tempkey;

         if( weights != NULL )
            weights[j] = tmpweight;

         SORTTPL_HASFIELD1( field1[j] = tempfield1; )
         SORTTPL_HASFIELD2( field2[j] = tempfield2; )
         SORTTPL_HASFIELD3( field3[j] = tempfield3; )
         SORTTPL_HASFIELD4( field4[j] = tempfield4; )
         SORTTPL_HASFIELD5( field5[j] = tempfield5; )
         SORTTPL_HASFIELD6( field6[j] = tempfield6; )
      }
   }
}

/** returns the index a,b, or c of the median element among key[a], key[b], and key[c] */
static
int SORTTPL_NAME(sorttpl_medianThree, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   a,                  /**< first index of the key array to consider */
   int                   b,                  /**< second index of the key array to consider */
   int                   c                   /**< third index of the array to consider */
   )
{
   assert(a >= 0);
   assert(b >= 0);
   assert(c >= 0);
   assert(a != b);
   assert(b != c);
   assert(c != a);
   /* let the elements in the unsorted order be a,b,c at positions start, mid, and end */
   if( SORTTPL_ISBETTER( key[a], key[b]) ) /* a <= b */
   {
      if( SORTTPL_ISBETTER( key[b], key[c]) ) /* b <= c */
         /* resulting permutation: a b c */
         return b;
      else /* b > c */
      {
         if( SORTTPL_ISBETTER( key[a], key[c]) ) /* a <= c */
            /* resulting permutation: a c b */
            return c;
         else
            /* resulting permutation: c a b */
            return a;
      }
   }
   else /* a > b */
   {
      if( SORTTPL_ISBETTER( key[b], key[c] ) )
      {
         if( SORTTPL_ISBETTER( key[a], key[c]) )
            /* resulting permutation: b a c */
            return a;
         else
            /* resulting permutation: b c a */
            return c;
      }
      else
         /* resulting permutation: c b a */
         return b;
   }
}
/* guess a median for the key array [start, ..., end] by using the median of the first, last, and middle element */
static
int SORTTPL_NAME(sorttpl_selectPivotIndex, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   start,              /**< first index of the key array to consider */
   int                   end                 /**< last index of the key array to consider */
   )
{
   int pivotindex;

   /* use the middle index on small arrays */
   if( end - start + 1 <= SORTTPL_SHELLSORTMAX )
      pivotindex = (start + end) / 2;
   else if( end - start + 1 < SORTTPL_MINSIZENINTHER )
   {
      /* select the median of the first, last, and middle element as pivot element */
      int mid = (start + end) / 2;
      pivotindex = SORTTPL_NAME(sorttpl_medianThree, SORTTPL_NAMEEXT)
            (key,
                  SORTTPL_HASPTRCOMPPAR(ptrcomp)
                  SORTTPL_HASINDCOMPPAR(indcomp)
                  SORTTPL_HASINDCOMPPAR(dataptr)
                  start, mid, end);
   }
   else
   {
      /* use the median of medians of nine evenly distributed elements of the key array */
      int gap = (end - start + 1) / 9;
      int median1;
      int median2;
      int median3;

      /* this should always hold */
      assert(start + 8 * gap <= end);

      /* collect 3 medians evenly distributed over the array */
      median1 = SORTTPL_NAME(sorttpl_medianThree, SORTTPL_NAMEEXT)
            (key,
               SORTTPL_HASPTRCOMPPAR(ptrcomp)
               SORTTPL_HASINDCOMPPAR(indcomp)
               SORTTPL_HASINDCOMPPAR(dataptr)
               start, start + gap, start + 2 * gap);

      median2 = SORTTPL_NAME(sorttpl_medianThree, SORTTPL_NAMEEXT)
            (key,
               SORTTPL_HASPTRCOMPPAR(ptrcomp)
               SORTTPL_HASINDCOMPPAR(indcomp)
               SORTTPL_HASINDCOMPPAR(dataptr)
               start + 3 * gap, start + 4 * gap, start + 5 * gap);
      median3 = SORTTPL_NAME(sorttpl_medianThree, SORTTPL_NAMEEXT)
            (key,
               SORTTPL_HASPTRCOMPPAR(ptrcomp)
               SORTTPL_HASINDCOMPPAR(indcomp)
               SORTTPL_HASINDCOMPPAR(dataptr)
               start + 6 * gap, start + 7 * gap, start + 8 * gap);

      /* compute and return the median of the medians */
      pivotindex = SORTTPL_NAME(sorttpl_medianThree, SORTTPL_NAMEEXT)
            (key,
               SORTTPL_HASPTRCOMPPAR(ptrcomp)
               SORTTPL_HASINDCOMPPAR(indcomp)
               SORTTPL_HASINDCOMPPAR(dataptr)
               median1, median2, median3);
   }

   return pivotindex;
}


/** quick-sort an array of pointers; pivot is the medial element */
static
void SORTTPL_NAME(sorttpl_qSort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   start,              /**< starting index */
   int                   end,                /**< ending index */
   SCIP_Bool             type                /**< TRUE, if quick-sort should start with with key[lo] < pivot <= key[hi], key[lo] <= pivot < key[hi] otherwise */
   )
{
   assert(start <= end);

   /* use quick-sort for long lists */
   while( end - start >= SORTTPL_SHELLSORTMAX )
   {
      SORTTPL_KEYTYPE pivotkey;
      int lo;
      int hi;
      int mid;

      /* select pivot element */
      mid = SORTTPL_NAME(sorttpl_selectPivotIndex, SORTTPL_NAMEEXT)
            (key,
                  SORTTPL_HASPTRCOMPPAR(ptrcomp)
                  SORTTPL_HASINDCOMPPAR(indcomp)
                  SORTTPL_HASINDCOMPPAR(dataptr)
                  start, end);
      pivotkey = key[mid];

      /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
      lo = start;
      hi = end;
      for( ;; )
      {
         if( type )
         {
            while( lo < end && SORTTPL_ISBETTER(key[lo], pivotkey) )
               lo++;
            while( hi > start && !SORTTPL_ISBETTER(key[hi], pivotkey) )
               hi--;
         }
         else
         {
            while( lo < end && !SORTTPL_ISWORSE(key[lo], pivotkey) )
               lo++;
            while( hi > start && SORTTPL_ISWORSE(key[hi], pivotkey) )
               hi--;
         }

         if( lo >= hi )
            break;

         SORTTPL_SWAP(SORTTPL_KEYTYPE, key[lo], key[hi]);
         SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[lo], field1[hi]); )
         SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[lo], field2[hi]); )
         SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[lo], field3[hi]); )
         SORTTPL_HASFIELD4( SORTTPL_SWAP(SORTTPL_FIELD4TYPE, field4[lo], field4[hi]); )
         SORTTPL_HASFIELD5( SORTTPL_SWAP(SORTTPL_FIELD5TYPE, field5[lo], field5[hi]); )
         SORTTPL_HASFIELD6( SORTTPL_SWAP(SORTTPL_FIELD6TYPE, field6[lo], field6[hi]); )

         lo++;
         hi--;
      }
      assert((hi == lo-1) || (type && hi == start) || (!type && lo == end));

      /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
      if( type )
      {
         while( lo < end && !SORTTPL_ISBETTER(pivotkey, key[lo]) )
            lo++;

         /* make sure that we have at least one element in the smaller partition */
         if( lo == start )
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert(!SORTTPL_ISBETTER(key[mid], pivotkey)); /* the pivot element did not change its position */
            assert(!SORTTPL_ISBETTER(pivotkey, key[mid]));
            SORTTPL_SWAP(SORTTPL_KEYTYPE, key[lo], key[mid]);
            SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[lo], field1[mid]); )
            SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[lo], field2[mid]); )
            SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[lo], field3[mid]); )
            SORTTPL_HASFIELD4( SORTTPL_SWAP(SORTTPL_FIELD4TYPE, field4[lo], field4[mid]); )
            SORTTPL_HASFIELD5( SORTTPL_SWAP(SORTTPL_FIELD5TYPE, field5[lo], field5[mid]); )
            SORTTPL_HASFIELD6( SORTTPL_SWAP(SORTTPL_FIELD6TYPE, field6[lo], field6[mid]); )
            lo++;
         }
      }
      else
      {
         while( hi > start && !SORTTPL_ISWORSE(pivotkey, key[hi]) )
            hi--;

         /* make sure that we have at least one element in the smaller partition */
         if( hi == end )
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert(!SORTTPL_ISBETTER(key[mid], pivotkey)); /* the pivot element did not change its position */
            assert(!SORTTPL_ISBETTER(pivotkey, key[mid]));
            SORTTPL_SWAP(SORTTPL_KEYTYPE, key[hi], key[mid]);
            SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[hi], field1[mid]); )
            SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[hi], field2[mid]); )
            SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[hi], field3[mid]); )
            SORTTPL_HASFIELD4( SORTTPL_SWAP(SORTTPL_FIELD4TYPE, field4[hi], field4[mid]); )
            SORTTPL_HASFIELD5( SORTTPL_SWAP(SORTTPL_FIELD5TYPE, field5[hi], field5[mid]); )
            SORTTPL_HASFIELD6( SORTTPL_SWAP(SORTTPL_FIELD6TYPE, field6[hi], field6[mid]); )
            hi--;
         }
      }

      /* sort the smaller partition by a recursive call, sort the larger part without recursion */
      if( hi - start <= end - lo )
      {
         /* sort [start,hi] with a recursive call */
         if( start < hi )
         {
            SORTTPL_NAME(sorttpl_qSort, SORTTPL_NAMEEXT)
               (key,
                SORTTPL_HASFIELD1PAR(field1)
                SORTTPL_HASFIELD2PAR(field2)
                SORTTPL_HASFIELD3PAR(field3)
                SORTTPL_HASFIELD4PAR(field4)
                SORTTPL_HASFIELD5PAR(field5)
                SORTTPL_HASFIELD6PAR(field6)
                SORTTPL_HASPTRCOMPPAR(ptrcomp)
                SORTTPL_HASINDCOMPPAR(indcomp)
                SORTTPL_HASINDCOMPPAR(dataptr)
                  start, hi, !type);
         }

         /* now focus on the larger part [lo,end] */
         start = lo;
      }
      else
      {
         if( lo < end )
         {
            /* sort [lo,end] with a recursive call */
            SORTTPL_NAME(sorttpl_qSort, SORTTPL_NAMEEXT)
               (key,
                SORTTPL_HASFIELD1PAR(field1)
                SORTTPL_HASFIELD2PAR(field2)
                SORTTPL_HASFIELD3PAR(field3)
                SORTTPL_HASFIELD4PAR(field4)
                SORTTPL_HASFIELD5PAR(field5)
                SORTTPL_HASFIELD6PAR(field6)
                SORTTPL_HASPTRCOMPPAR(ptrcomp)
                SORTTPL_HASINDCOMPPAR(indcomp)
                SORTTPL_HASINDCOMPPAR(dataptr)
                  lo, end, !type);
         }

         /* now focus on the larger part [start,hi] */
         end = hi;
      }
      type = !type;
   }

   /* use shell sort on the remaining small list */
   if( end - start >= 1 )
   {
      SORTTPL_NAME(sorttpl_shellSort, SORTTPL_NAMEEXT)
         (key, NULL,
            SORTTPL_HASFIELD1PAR(field1)
            SORTTPL_HASFIELD2PAR(field2)
            SORTTPL_HASFIELD3PAR(field3)
            SORTTPL_HASFIELD4PAR(field4)
            SORTTPL_HASFIELD5PAR(field5)
            SORTTPL_HASFIELD6PAR(field6)
            SORTTPL_HASPTRCOMPPAR(ptrcomp)
            SORTTPL_HASINDCOMPPAR(indcomp)
            SORTTPL_HASINDCOMPPAR(dataptr)
            start, end);
   }
}

#ifndef NDEBUG
/** verifies that an array is indeed sorted */
static
void SORTTPL_NAME(sorttpl_checkSort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of the array */
   )
{
   int i;

   for( i = 0; i < len-1; i++ )
   {
      assert(!SORTTPL_ISBETTER(key[i+1], key[i]));
   }
}
#endif

/** SCIPsort...(): sorts array 'key' and performs the same permutations on the additional 'field' arrays */
void SORTTPL_NAME(SCIPsort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of arrays */
   )
{
   /* ignore the trivial cases */
   if( len <= 1 )
      return;

   /* use shell sort on the remaining small list */
   if( len <= SORTTPL_SHELLSORTMAX )
   {
      SORTTPL_NAME(sorttpl_shellSort, SORTTPL_NAMEEXT)
         (key, NULL,
            SORTTPL_HASFIELD1PAR(field1)
            SORTTPL_HASFIELD2PAR(field2)
            SORTTPL_HASFIELD3PAR(field3)
            SORTTPL_HASFIELD4PAR(field4)
            SORTTPL_HASFIELD5PAR(field5)
            SORTTPL_HASFIELD6PAR(field6)
            SORTTPL_HASPTRCOMPPAR(ptrcomp)
            SORTTPL_HASINDCOMPPAR(indcomp)
            SORTTPL_HASINDCOMPPAR(dataptr)
            0, len-1);
   }
   else
   {
      SORTTPL_NAME(sorttpl_qSort, SORTTPL_NAMEEXT)
         (key,
            SORTTPL_HASFIELD1PAR(field1)
            SORTTPL_HASFIELD2PAR(field2)
            SORTTPL_HASFIELD3PAR(field3)
            SORTTPL_HASFIELD4PAR(field4)
            SORTTPL_HASFIELD5PAR(field5)
            SORTTPL_HASFIELD6PAR(field6)
            SORTTPL_HASPTRCOMPPAR(ptrcomp)
            SORTTPL_HASINDCOMPPAR(indcomp)
            SORTTPL_HASINDCOMPPAR(dataptr)
            0, len-1, TRUE);
   }
#ifndef NDEBUG
   SORTTPL_NAME(sorttpl_checkSort, SORTTPL_NAMEEXT)
      (key,
       SORTTPL_HASPTRCOMPPAR(ptrcomp)
       SORTTPL_HASINDCOMPPAR(indcomp)
       SORTTPL_HASINDCOMPPAR(dataptr)
       len);
#endif
}


/** SCIPsortedvecInsert...(): adds an element to a sorted multi-vector;
 *  This method does not do any memory allocation! It assumes that the arrays are large enough
 *  to store the additional values.
 */
void SORTTPL_NAME(SCIPsortedvecInsert, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   SORTTPL_KEYTYPE       keyval,             /**< key value of new element */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE     field1val  )  /**< field1 value of new element */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE     field2val  )  /**< field1 value of new element */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE     field3val  )  /**< field1 value of new element */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE     field4val  )  /**< field1 value of new element */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE     field5val  )  /**< field1 value of new element */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE     field6val  )  /**< field1 value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insert position, or NULL */
   )
{
   int j;

   for( j = *len; j > 0 && SORTTPL_ISBETTER(keyval, key[j-1]); j-- )
   {
      key[j] = key[j-1];
      SORTTPL_HASFIELD1( field1[j] = field1[j-1]; )
      SORTTPL_HASFIELD2( field2[j] = field2[j-1]; )
      SORTTPL_HASFIELD3( field3[j] = field3[j-1]; )
      SORTTPL_HASFIELD4( field4[j] = field4[j-1]; )
      SORTTPL_HASFIELD5( field5[j] = field5[j-1]; )
      SORTTPL_HASFIELD6( field6[j] = field6[j-1]; )
   }

   key[j] = keyval;
   SORTTPL_HASFIELD1( field1[j] = field1val; )
   SORTTPL_HASFIELD2( field2[j] = field2val; )
   SORTTPL_HASFIELD3( field3[j] = field3val; )
   SORTTPL_HASFIELD4( field4[j] = field4val; )
   SORTTPL_HASFIELD5( field5[j] = field5val; )
   SORTTPL_HASFIELD6( field6[j] = field6val; )

   (*len)++;

   if( pos != NULL )
      (*pos) = j;
}

/** SCIPsortedvecDelPos...(): deletes an element at a given position from a sorted multi-vector */
void SORTTPL_NAME(SCIPsortedvecDelPos, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   )
{
   int j;

   assert(0 <= pos && pos < *len);

   (*len)--;

   for( j = pos; j < *len; j++ )
   {
      key[j] = key[j+1];
      SORTTPL_HASFIELD1( field1[j] = field1[j+1]; )
      SORTTPL_HASFIELD2( field2[j] = field2[j+1]; )
      SORTTPL_HASFIELD3( field3[j] = field3[j+1]; )
      SORTTPL_HASFIELD4( field4[j] = field4[j+1]; )
      SORTTPL_HASFIELD5( field5[j] = field5[j+1]; )
      SORTTPL_HASFIELD6( field6[j] = field6[j+1]; )
   }
}


/* The SCIPsortedvecFind...() method only has needs the key array but not the other field arrays. In order to
 * avoid defining the same method multiple times, only include this method if we do not have any additional fields.
 */
#ifndef SORTTPL_FIELD1TYPE

/** SCIPsortedvecFind...(): Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
SCIP_Bool SORTTPL_NAME(SCIPsortedvecFind, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   SORTTPL_KEYTYPE       val,                /**< data field to find position for */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store the insert position */
   )
{
   int left;
   int right;

   assert(key != NULL);
   assert(pos != NULL);

   left = 0;
   right = len-1;
   while( left <= right )
   {
      int middle;

      middle = (left+right)/2;
      assert(0 <= middle && middle < len);

      if( SORTTPL_ISBETTER(val, key[middle]) )
         right = middle-1;
      else if( SORTTPL_ISBETTER(key[middle], val) )
         left = middle+1;
      else
      {
         *pos = middle;
         return TRUE;
      }
   }
   assert(left == right+1);

   *pos = left;
   return FALSE;
}

#endif



/** macro that performs an exchange in the weighted selection algorithm, including weights */
#define EXCH(x,y)                                                                              \
   do                                                                                          \
   {                                                                                           \
      SORTTPL_SWAP(SORTTPL_KEYTYPE, key[x], key[y]);                                           \
                                                                                               \
      if( weights != NULL )                                                                    \
         SORTTPL_SWAP(SCIP_Real, weights[x], weights[y]);                                      \
                                                                                               \
      SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[x], field1[y]); )             \
      SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[x], field2[y]); )             \
      SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[x], field3[y]); )             \
      SORTTPL_HASFIELD4( SORTTPL_SWAP(SORTTPL_FIELD4TYPE, field4[x], field4[y]); )             \
      SORTTPL_HASFIELD5( SORTTPL_SWAP(SORTTPL_FIELD5TYPE, field5[x], field5[y]); )             \
      SORTTPL_HASFIELD6( SORTTPL_SWAP(SORTTPL_FIELD6TYPE, field6[x], field6[y]); )             \
   }                                                                                           \
   while( FALSE )

/** partially sorts a given keys array around the weighted median w.r.t. the \p capacity and permutes the additional 'field' arrays
 *  in the same way.
 *
 *  If no weights-array is passed, the algorithm assumes weights equal to 1.
 */
void SORTTPL_NAME(SCIPselectWeighted, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   )
{
   int hi;
   int lo;
   int j;
   int recursiondepth;
   SCIP_Real residualcapacity;

   lo = 0;
   hi = len - 1;
   residualcapacity = capacity;
   recursiondepth = 0;

   while( hi - lo + 1 > SORTTPL_SHELLSORTMAX )
   {
      int i;
      int bt;
      int wt;
      int p;
      int pivotindex;
      SCIP_Real betterweightsum;
      SCIP_Real pivotweight;
      SORTTPL_KEYTYPE pivot;

      ++recursiondepth;

      /* guess a median as pivot */
      pivotindex = SORTTPL_NAME(sorttpl_selectPivotIndex, SORTTPL_NAMEEXT)
            (key,
                  SORTTPL_HASPTRCOMPPAR(ptrcomp)
                  SORTTPL_HASINDCOMPPAR(indcomp)
                  SORTTPL_HASINDCOMPPAR(dataptr)
                  lo, hi);

      pivot = key[pivotindex];

      /* swap pivot element to the end of the array */
      if( pivotindex != lo )
      {
         EXCH(lo, pivotindex);
      }

      /* initialize array indices for the current element, the better elements, and the worse elements */
      i = lo;
      bt = lo;
      wt = hi;

      /* iterate through elements once to establish three partition into better elements, equal elements, and worse elements
       *
       * at every iteration, i denotes the current, previously unseen element, starting from the position lo
       * all elements [lo,...bt - 1] are better than the pivot
       * all elements [wt + 1,... hi] are worse than the pivot
       *
       * at termination, all elements [bt,...wt] are equal to the pivot element
       * */
      while( i <= wt )
      {
         /* element i is better than pivot; exchange elements i and bt, increase both */
         if( SORTTPL_ISBETTER(key[i], pivot) )
         {
            EXCH(i, bt);
            i++;
            bt++;
         }
         /* element i is worse than pivot: exchange it with the element at position wt; no increment of i
          * because an unseen element is waiting at index i after the swap
          */
         else if( SORTTPL_ISWORSE(key[i], pivot) )
         {
           EXCH(i, wt);
           wt--;
         }
         else
            i++;
      }

      assert(wt >= bt);

      if( weights != NULL )
      {
         betterweightsum = 0.0;
         /* collect weights of elements larger than the pivot  */
         for( i = lo; i < bt; ++i )
         {
            assert(SORTTPL_ISBETTER(key[i], pivot));
            betterweightsum += weights[i];
         }
      }
      else
      {
         /* if all weights are equal to one, we directly know the larger and the equal weight sum */
         betterweightsum = bt - lo;
      }

      /* the weight in the better half of the array exceeds the capacity. Continue the search there */
      if( betterweightsum >= residualcapacity )
      {
         hi = bt - 1;
      }
      else
      {
         SCIP_Real weightsum = betterweightsum;
         /* loop through duplicates of pivot element and check if one is the weighted median */
         for( p = bt; p <= wt; ++p )
         {
            assert(SORTTPL_CMP(key[p], pivot) == 0);
            pivotweight = weights != NULL ? weights[p] : 1;
            weightsum += pivotweight;

            /* the element at index p is exactly the weighted median */
            if( weightsum >= residualcapacity )
            {
               if( medianpos != NULL )
                  *medianpos = p;

               return;
            }
         }

         /* continue loop by searching the remaining elements [wt+1,...,hi] */
         residualcapacity -= weightsum;
         lo = wt + 1;
      }
   }

   assert(hi - lo + 1 <= SORTTPL_SHELLSORTMAX);

   /* use shell sort to solve the remaining elements completely */
   if( hi - lo + 1 > 1 )
   {
      SORTTPL_NAME(sorttpl_shellSort, SORTTPL_NAMEEXT)
      (key, weights,
         SORTTPL_HASFIELD1PAR(field1)
         SORTTPL_HASFIELD2PAR(field2)
         SORTTPL_HASFIELD3PAR(field3)
         SORTTPL_HASFIELD4PAR(field4)
         SORTTPL_HASFIELD5PAR(field5)
         SORTTPL_HASFIELD6PAR(field6)
         SORTTPL_HASPTRCOMPPAR(ptrcomp)
         SORTTPL_HASINDCOMPPAR(indcomp)
         SORTTPL_HASINDCOMPPAR(dataptr)
         lo, hi);
   }
   /* determine the median position among the remaining elements */
   for( j = lo; j <= MAX(lo, hi); ++j )
   {
      SCIP_Real weight = (weights != NULL ? weights[j] : 1);
      /* we finally found the median element */
      if( weight >= residualcapacity )
      {
         if( medianpos != NULL )
            *medianpos = j;

         break;
      }
      else
         residualcapacity -= weight;
   }

   /* the capacity is not exceeded by the elements in the array */
   if( j == len && medianpos != NULL )
   {
      assert(residualcapacity > 0);
      *medianpos = len;
   }
}

/** partially sorts a given keys array around the given index \p k and permutes the additional 'field' arrays are in the same way */
void SORTTPL_NAME(SCIPselect, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key,                /**< pointer to data array that defines the order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */
   SORTTPL_HASPTRCOMPPAR( SCIP_DECL_SORTPTRCOMP((*ptrcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( SCIP_DECL_SORTINDCOMP((*indcomp)) )  /**< data element comparator */
   SORTTPL_HASINDCOMPPAR( void*                  dataptr    )  /**< pointer to data field that is given to the external compare method */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   )
{
   SCIP_Real capacity;
   int pos;

   /* return directly in cases that make no sense at all */
   if( k < 0 || k >= len )
      return;

   /* The summand 0.5 is necessary because the elements are zero-indexed. */
   capacity = k + 0.5;

   pos = -1;

   /* call the general algorithm for the weighted median with weights equal to -1 (by passing NULL) */
   SORTTPL_NAME(SCIPselectWeighted, SORTTPL_NAMEEXT)
   (key,
      SORTTPL_HASFIELD1PAR(field1)
      SORTTPL_HASFIELD2PAR(field2)
      SORTTPL_HASFIELD3PAR(field3)
      SORTTPL_HASFIELD4PAR(field4)
      SORTTPL_HASFIELD5PAR(field5)
      SORTTPL_HASFIELD6PAR(field6)
      SORTTPL_HASPTRCOMPPAR(ptrcomp)
      SORTTPL_HASINDCOMPPAR(indcomp)
      SORTTPL_HASINDCOMPPAR(dataptr)
      NULL, capacity, len, &pos);

   /* the weighted median position should be exactly at position k */
   assert(pos == k);
}

/* undefine template parameters and local defines */
#undef SORTTPL_NAMEEXT
#undef SORTTPL_KEYTYPE
#undef SORTTPL_FIELD1TYPE
#undef SORTTPL_FIELD2TYPE
#undef SORTTPL_FIELD3TYPE
#undef SORTTPL_FIELD4TYPE
#undef SORTTPL_FIELD5TYPE
#undef SORTTPL_FIELD6TYPE
#undef SORTTPL_PTRCOMP
#undef SORTTPL_INDCOMP
#undef SORTTPL_HASFIELD1
#undef SORTTPL_HASFIELD2
#undef SORTTPL_HASFIELD3
#undef SORTTPL_HASFIELD4
#undef SORTTPL_HASFIELD5
#undef SORTTPL_HASFIELD6
#undef SORTTPL_HASPTRCOMP
#undef SORTTPL_HASINDCOMP
#undef SORTTPL_HASFIELD1PAR
#undef SORTTPL_HASFIELD2PAR
#undef SORTTPL_HASFIELD3PAR
#undef SORTTPL_HASFIELD4PAR
#undef SORTTPL_HASFIELD5PAR
#undef SORTTPL_HASFIELD6PAR
#undef SORTTPL_HASPTRCOMPPAR
#undef SORTTPL_HASINDCOMPPAR
#undef SORTTPL_ISBETTER
#undef SORTTPL_ISWORSE
#undef SORTTPL_CMP
#undef EXCH
#undef SORTTPL_SWAP
#undef SORTTPL_SHELLSORTMAX
#undef SORTTPL_MINSIZENINTHER
#undef SORTTPL_BACKWARDS
