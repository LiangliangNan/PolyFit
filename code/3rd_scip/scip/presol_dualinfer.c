/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file    presol_dualinfer.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief   dual inference presolver
 * @author  Dieter Weninger
 * @author  Patrick Gemander
 *
 * This presolver does bound strengthening on continuous variables (columns) for getting bounds on dual variables y.
 * The bounds of the dual variables are then used to fix primal variables or change the side of constraints.
 * For ranged rows one needs to decide which side (rhs or lhs) determines the equality.
 *
 * We distinguish two cases concerning complementary slackness:
 * i)  reduced cost fixing:       c_j - sup_y(y^T A_{.j}) > 0 => x_j = l_j
 *                                c_j - inf_y(y^T A_{.j}) < 0 => x_j = u_j
 * ii) positive dual lower bound: y_i > 0 =>  A_{i.}x = b_i
 *
 * Further information on this presolving approach are given in
 * Achterberg et al. "Presolve reductions in mixed integer programming"
 * and for a two-column extension in
 * Chen et al. "Two-row and two-column mixed-integer presolve using hasing-based pairing methods".
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"
#include "scip/pub_matrix.h"
#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/presol_dualinfer.h"
#include "scip/pub_cons.h"
#include "scip/pub_matrix.h"
#include "scip/pub_message.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_presol.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_var.h"

#define PRESOL_NAME             "dualinfer"
#define PRESOL_DESC             "exploit dual information for fixings and side changes"
#define PRESOL_PRIORITY                -3000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS                0    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING                   SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_TWOCOLUMN_COMBINE       TRUE /**< should two column convex combination be used per default */
#define DEFAULT_MAXLOOPS_DUALBNDSTR     12   /**< default maximal number of loops for dual bound strengthening */
#define DEFAULT_MAXCONSIDEREDNONZEROS   100  /**< default maximal number of considered non-zeros within one row */
#define DEFAULT_MAXRETRIEVEFAILS        1000 /**< default maximal number of consecutive useless hashtable retrieves */
#define DEFAULT_MAXCOMBINEFAILS         1000 /**< default maximal number of consecutive useless row combines */
#define DEFAULT_MAXHASHFAC              10   /**< default maximal number of hashlist entries as multiple of number of rows in the problem */
#define DEFAULT_MAXPAIRFAC              1    /**< default maximal number of processed row pairs as multiple of the number of rows in the problem */
#define DEFAULT_MAXROWSUPPORT           3    /**< default maximal number of non-zeros in one row for turning an inequality into an equality */


/*
 * Data structures
 */

/** control parameters */
struct SCIP_PresolData
{
   SCIP_Bool usetwocolcombine;               /**< use convex combination of two columns */
   int maxdualbndloops;                      /**< default number of dual bound strengthening loops */
   int maxpairfac;                           /**< maximal number of processed row pairs as multiple of the number of rows in the problem (-1: no limit) */
   int maxhashfac;                           /**< maximal number of hashlist entries as multiple of number of rows in the problem (-1: no limit) */
   int maxretrievefails;                     /**< maximal number of consecutive useless hashtable retrieves */
   int maxcombinefails;                      /**< maximal number of consecutive useless row combines */
   int maxconsiderednonzeros;                /**< maximal number of considered non-zeros within one row (-1: no limit) */
   int maxrowsupport;                        /**< maximal number of non-zeros in one row for turning an inequality into an equality */
};

/** type of variable fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,                             /** fix variable at its lower bound */
   NOFIX   =  0,                             /** no fixing */
   FIXATUB =  1                              /** fix variable at its upper bound */
};
typedef enum Fixingdirection FIXINGDIRECTION;

/** type of constraint side change */
enum SideChange
{
   RHSTOLHS = -1,                            /** set rhs to value of lhs */
   NOCHANGE =  0,                            /** no side change */
   LHSTORHS =  1                             /** set lhs to value of rhs */
};
typedef enum SideChange SIDECHANGE;

/** Signum for convex-combined variable coefficients \f$(\lambda * A_{ri} + (1 - \lambda) * A_{si})\f$
 *  UP  - Coefficient changes from negative to positive for increasing lambda
 *  DN  - Coefficient changes from positive to negative for increasing lambda
 *  POS - Coefficient is positive for all lambda in (0,1)
 *  NEG - Coefficient is negative for all lambda in (0,1)
 */
enum signum {UP, DN, POS, NEG};

/** structure representing a pair of column indices; used for lookup in a hashtable */
struct ColPair
{
   int col1idx;                              /**< first row index */
   int col2idx;                              /**< second row index */
};
typedef struct ColPair COLPAIR;

/*
 * Local methods
 */

/** encode contents of a colpair as void* pointer */
static
void*
encodeColPair(
   COLPAIR*              colpair             /**< pointer to colpair */
   )
{
   uint64_t a;
   uint64_t b;

   assert(colpair->col1idx >= 0);
   assert(colpair->col2idx >= 0);

   a = (uint64_t)(long)colpair->col1idx;
   b = (uint64_t)(long)colpair->col2idx;
   return (void*)((a << 32) | b);
}

/** compute single positive int hashvalue for two ints */
static
int
hashIndexPair(
   int                   idx1,               /**< first integer index */
   int                   idx2                /**< second integer index */
   )
{
   uint32_t hash = SCIPhashTwo(idx1, idx2);
   return (int)(hash>>1);
}

/** add hash/rowidx pair to hashlist/rowidxlist **/
static
SCIP_RETCODE addEntry(
   SCIP*                 scip,               /**< SCIP datastructure */
   int*                  pos,                /**< position of last entry added */
   int*                  listsize,           /**< size of hashlist and rowidxlist */
   int**                 hashlist,           /**< block memory array containing hashes */
   int**                 colidxlist,         /**< block memory array containing column indices */
   int                   hash,               /**< hash to be inserted */
   int                   colidx              /**< column index to be inserted */
   )
{
   if( (*pos) >= (*listsize) )
   {
      int newsize  = SCIPcalcMemGrowSize(scip, (*pos) + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, hashlist, (*listsize), newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, colidxlist, (*listsize), newsize) );
      (*listsize) = newsize;
   }

   (*hashlist)[(*pos)] = hash;
   (*colidxlist)[(*pos)] = colidx;
   (*pos)++;

   return SCIP_OKAY;
}

/** Within a sorted list, get next block with same value
 *  E.g. for [h1, h1, h1, h2, h2, h2, h2, h3,...] and end = 0
 *  returns start = 0, end = 3
 *  and on a second call with end = 3 on the same list
 *  returns start = 3, end = 7.
 */
static
void findNextBlock(
   int*                  list,               /**< list of integers */
   int                   len,                /**< length of list */
   int*                  start,              /**< variable to contain start index of found block */
   int*                  end                 /**< variable to contain end index of found block */
   )
{
   int i;
   (*start) = (*end);
   i = (*end) + 1;
   while( i < len && list[i] == list[i - 1] )
      i++;

   (*end) = i;
}

/**
 * The algorithm described in Belotti P. "Bound reduction using pairs of linear inequalities"
 * tries to derive upper and lower bounds for all variables via convex combinations of linear inequalities
 * We apply Belotti's algorithm to pairs of columns of continuous variables.
 */
static
SCIP_RETCODE combineCols(
   SCIP*                 scip,               /**< SCIP datastructure */
   int*                  row1idxptr,         /**< indices specifying bound positions in lbs and ubs for first row */
   int*                  row2idxptr,         /**< indices specifying bound positions in lbs und ubs for second row */
   SCIP_Real*            row1valptr,         /**< first row coefficients */
   SCIP_Real*            row2valptr,         /**< second row coefficients */
   SCIP_Real             b1,                 /**< rhs of first row */
   SCIP_Real             b2,                 /**< rhs of second row*/
   int                   row1len,            /**< length of first row (e.g. row1idxptr and row1valptr)*/
   int                   row2len,            /**< length of second row (e.g. row2idxptr and row2valptr)*/
   int                   ncols,              /**< length of bound arrays lbs and ubs */
   SCIP_Bool             swaprow1,           /**< should the sense of the first row be swapped to <= ? */
   SCIP_Bool             swaprow2,           /**< should the sense of the second row be swapped to <= ? */
   SCIP_Real*            lbs,                /**< lower bound array */
   SCIP_Real*            ubs,                /**< upper bound array */
   SCIP_Bool*            success             /**< we return (success || found better bounds") */
   )
{
   int i;
   int j;
   int nvars;
   int* varinds;
   int nbreakpoints;
   SCIP_Real* breakpoints;
   int idx;
   int idx1;
   int idx2;
   SCIP_Real* row1coefs;
   SCIP_Real* row2coefs;
   enum signum* signs;
   int ninfs;
   int l1infs;
   SCIP_Real  l1;
   SCIP_Real  l2;
   SCIP_Real* newlbs;
   SCIP_Real* newubs;
   SCIP_Real coef;
   int sign;
   int shift;

   SCIP_CALL( SCIPallocBufferArray(scip, &row1coefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &row2coefs, ncols) );

   SCIPsortIntReal(row1idxptr, row1valptr, row1len);
   SCIPsortIntReal(row2idxptr, row2valptr, row2len);

   /* swap rows if necessary */
   if( swaprow1 )
   {
      for( i = 0; i < row1len; i++ )
         row1coefs[row1idxptr[i]] = -row1valptr[i];
      b1 = -b1;
   }
   else
   {
      for( i = 0; i < row1len; i++ )
         row1coefs[row1idxptr[i]] = row1valptr[i];
   }

   if( swaprow2 )
   {
      for( i = 0; i < row2len; i++ )
         row2coefs[row2idxptr[i]] = -row2valptr[i];
      b2 = -b2;
   }
   else
   {
      for( i = 0; i < row2len; i++ )
         row2coefs[row2idxptr[i]] = row2valptr[i];
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varinds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &signs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &breakpoints, ncols) );

   /* calculate cancellation breakpoints and sign behaviour */
   i = 0;
   j = 0;
   nvars = 0;
   nbreakpoints = 0;
   while( i < row1len && j < row2len )
   {
      assert(i + 1 == row1len || row1idxptr[i] < row1idxptr[i + 1]);
      assert(j + 1 == row2len || row2idxptr[j] < row2idxptr[j + 1]);

      idx1 = row1idxptr[i];
      idx2 = row2idxptr[j];

      /* We use 2.0 as default value for "no cancellation". For cancellations, this will be replaced by values in (0,1).
       * A value larger than 1.0 is used because we sort the array and want to put non-cancellations to the end. */
      breakpoints[nvars] = 2.0;

      if( idx1 == idx2 )
      {
         if( (SCIPisNegative(scip, row1coefs[idx1]) && SCIPisPositive(scip, row2coefs[idx2])) ||
            (SCIPisPositive(scip, row1coefs[idx1]) && SCIPisNegative(scip, row2coefs[idx2])) )
         {
            if( SCIPisNegative(scip, row2coefs[idx2]) )
               signs[idx1] = UP;
            else
               signs[idx1] = DN;

            breakpoints[nvars] = row2coefs[idx2] / (row2coefs[idx2] - row1coefs[idx1]);
            nbreakpoints++;
         }
         else if( SCIPisPositive(scip, row1coefs[idx1]) )
            signs[idx1] = POS;
         else
            signs[idx1] = NEG;

         varinds[nvars] = idx1;
         i++;
         j++;
      }
      else if( idx1 < idx2 )
      {
         if( SCIPisPositive(scip, row1coefs[idx1]) )
            signs[idx1] = POS;
         else
            signs[idx1] = NEG;

         /* We will access this entry later on, so we explicitly write a zero here */
         row2coefs[idx1] = 0.0;

         varinds[nvars] = idx1;
         i++;
      }
      else
      {
         assert(idx1 > idx2);
         if( SCIPisPositive(scip, row2coefs[idx2]) )
            signs[idx2] = POS;
         else
            signs[idx2] = NEG;

         /* We will access this entry later on, so we explicitly write a zero here */
         row1coefs[idx2] = 0.0;

         varinds[nvars] = idx2;
         j++;
      }
      nvars++;
   }

   while( i < row1len )
   {
      idx1 = row1idxptr[i];

      if( SCIPisPositive(scip, row1coefs[idx1]) )
         signs[idx1] = POS;
      else
         signs[idx1] = NEG;

      /* We will access this entry later on, so we explicitly write a zero here */
      row2coefs[idx1] = 0.0;

      varinds[nvars] = idx1;
      breakpoints[nvars] = 2.0;
      nvars++;
      i++;
   }

   while( j < row2len )
   {
      idx2 = row2idxptr[j];

      if( SCIPisPositive(scip, row2coefs[idx2]) )
         signs[idx2] = POS;
      else
         signs[idx2] = NEG;

      /* We will access this entry later on, so we explicitly write a zero here */
      row1coefs[idx2] = 0.0;

      varinds[nvars] = idx2;
      breakpoints[nvars] = 2.0;
      nvars++;
      j++;
   }

   SCIPsortRealInt(breakpoints, varinds, nvars);

   /* The obvious preconditions for bound tightenings are met, so we try to calculate new bounds. */
   if( nbreakpoints >= 1 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &newlbs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newubs, nvars) );

      for( i = 0; i < nvars; i++)
      {
         idx = varinds[i];
         newlbs[i] = lbs[idx];
         newubs[i] = ubs[idx];
      }

      /* calculate activity contributions of each row */
      l1 = b1;
      l2 = b2;
      l1infs = 0;
      ninfs = 0;
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];
         if( !SCIPisZero(scip, row2coefs[idx]) )
         {
            if( SCIPisNegative(scip, row2coefs[idx]) )
            {
               if( !SCIPisInfinity(scip, -lbs[idx]) )
               {
                  l1 -= row1coefs[idx] * lbs[idx];
                  l2 -= row2coefs[idx] * lbs[idx];
               }
               else
                  ninfs++;
            }
            else
            {
               /* coefficient of second row is positive */
               if( !SCIPisInfinity(scip, ubs[idx]) )
               {
                  l1 -= row1coefs[idx] * ubs[idx];
                  l2 -= row2coefs[idx] * ubs[idx];
               }
               else
                  ninfs++;
            }
         }
         else
         {
            /* since row2coefs[idx] is zero, we have to choose the bound using row1coefs[idx] */
            assert(!SCIPisZero(scip, row1coefs[idx]) && SCIPisZero(scip, row2coefs[idx]));
            if( SCIPisNegative(scip, row1coefs[idx]) )
            {
               if( !SCIPisInfinity(scip, -lbs[idx]) )
                  l1 -= row1coefs[idx] * lbs[idx];
               else
                  l1infs++;
            }
            else
            {
               /* coefficient of first row is positive */
               if( !SCIPisInfinity(scip, ubs[idx]) )
                  l1 -= row1coefs[idx] * ubs[idx];
               else
                  l1infs++;
            }
         }
      }

      /* Calculate bounds for lambda = 0 */
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(scip, "lambda = 0, l1 = %g, l2 = %g, ninfs = %d\n", i, breakpoints[i], l1, l2, ninfs);
#endif

      if( ninfs <= 1 )
      {
#ifdef SCIP_MORE_DEBUG
         SCIP_Real oldlb;
         SCIP_Real oldub;
#endif
         for( i = 0; i < nvars; i++ )
         {
#ifdef SCIP_MORE_DEBUG
            oldlb = newlbs[i];
            oldub = newubs[i];
#endif
            idx = varinds[i];
            if( SCIPisPositive(scip, row2coefs[idx]) )
            {
               if( ninfs == 0 )
                  newlbs[i] = MAX(newlbs[i], (l2 + row2coefs[idx] * ubs[idx]) / row2coefs[idx]);
               else if( SCIPisInfinity(scip, ubs[idx]) )
                  newlbs[i] = MAX(newlbs[i], l2 / row2coefs[idx]);
            }
            else if ( SCIPisNegative(scip, row2coefs[idx]) )
            {
               if( ninfs == 0 )
                  newubs[i] = MIN(newubs[i], (l2 + row2coefs[idx] * lbs[idx]) / row2coefs[idx]);
               else if( SCIPisInfinity(scip, -lbs[idx]) )
                  newubs[i] = MIN(newubs[i], l2 / row2coefs[idx]);
            }
#ifdef SCIP_MORE_DEBUG
            if( !SCIPisEQ(scip, oldlb, newlbs[i]) || !SCIPisEQ(scip, oldub, newubs[i]) )
               SCIPdebugMsg(scip, "%g <= %g <= var_%d <= %g <= %g\n", oldlb, newlbs[i], i, newubs[i], oldub);
#endif
         }
      }

      ninfs += l1infs;

      i = 0;
      while( i < nbreakpoints )
      {
         int nnewinfs;
         SCIP_Real l1update;
         SCIP_Real l2update;
         SCIP_Bool updated;

         /* determine number of infinities and compute update for l1 and l2 */
         shift = 0;
         nnewinfs = 0;
         l1update = 0.0;
         l2update = 0.0;
         updated = FALSE;
         j = i;
         while( !updated )
         {
            idx = varinds[j];
            assert(signs[idx] == UP || signs[idx] == DN);
            if( signs[idx] == UP )
               sign = 1;
            else
               sign = -1;

            if( !SCIPisInfinity(scip, -lbs[idx]) )
            {
               l1update += sign * row1coefs[idx] * lbs[idx];
               l2update += sign * row2coefs[idx] * lbs[idx];
            }
            else
            {
               if( signs[idx] == UP  )
                  ninfs--;
               else
                  nnewinfs++;
            }

            if( !SCIPisInfinity(scip, ubs[idx]) )
            {
               l1update -= sign * row1coefs[idx] * ubs[idx];
               l2update -= sign * row2coefs[idx] * ubs[idx];
            }
            else
            {
               if( signs[idx] == UP  )
                  nnewinfs++;
               else
                  ninfs--;
            }

            if( signs[idx] == UP )
               signs[idx] = POS;
            else
               signs[idx] = NEG;

            if( j + 1 >= nbreakpoints || !SCIPisEQ(scip, breakpoints[j], breakpoints[j + 1]) )
               updated = TRUE;

            shift++;
            j++;
         }

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "lambda_%d = %g, l1 = %g, l2 = %g, ninfs = %d\n", i, breakpoints[i], l1, l2, ninfs);
#endif

         assert(ninfs >= 0);

         /* if more than one infinity destroys our bounds we cannot tighten anything */
         if( ninfs <= 1 )
         {
            /* check for bounds to be tightened */
            for( j = 0; j < nvars; j++ )
            {
#ifdef SCIP_MORE_DEBUG
               SCIP_Real oldlb;
               SCIP_Real oldub;
#endif

               /* catch the special case where the entire remaining constraint is cancelled */
               if( j >= nvars )
                  break;

#ifdef SCIP_MORE_DEBUG
               oldlb = newlbs[j];
               oldub = newubs[j];
#endif

               idx = varinds[j];
               coef = breakpoints[i] * row1coefs[idx] + (1 - breakpoints[i]) * row2coefs[idx];
               assert(!SCIPisEQ(scip, breakpoints[i], 2.0));

               /* skip if the coefficient is too close to zero as it becomes numerically unstable */
               if( SCIPisZero(scip, coef) )
                  continue;

               if( signs[idx] == POS || signs[idx] == DN )
               {
                  if( ninfs == 0 )
                     newlbs[j] = MAX(newlbs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2 + coef * ubs[idx]) / coef);
                  else if( SCIPisInfinity(scip, ubs[idx]) )
                     newlbs[j] = MAX(newlbs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2) / coef);
               }
               else if ( signs[idx] == NEG || signs[idx] == UP )
               {
                  if( ninfs == 0 )
                     newubs[j] = MIN(newubs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2 + coef * lbs[idx]) / coef);
                  else if( SCIPisInfinity(scip, -lbs[idx]) )
                     newubs[j] = MIN(newubs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2) / coef);
               }
#ifdef SCIP_MORE_DEBUG
               if( !SCIPisEQ(scip, oldlb, newlbs[j]) || !SCIPisEQ(scip, oldub, newubs[j]) )
                  SCIPdebugMsg(scip, "%g <= %g <= var_%d <= %g <= %g\n", oldlb, newlbs[j], j, newubs[j], oldub);
#endif
            }
         }

         i += shift;
         ninfs += nnewinfs;
         l1 += l1update;
         l2 += l2update;
      }

      /* check infinities in first row */
      ninfs = 0;
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];
         if( (SCIPisPositive(scip, row1coefs[idx]) && SCIPisInfinity(scip, ubs[idx]))
            || (SCIPisNegative(scip, row1coefs[idx]) && SCIPisInfinity(scip, -lbs[idx])) )
            ninfs++;
      }

      /* calculate bounds for lambda = 1 */
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(scip, "lambda = 1, l1 = %g, l2 = %g, ninfs = %d\n", i, breakpoints[i], l1, l2, ninfs);
#endif
      if( ninfs <= 1 )
      {
#ifdef SCIP_MORE_DEBUG
         SCIP_Real oldlb;
         SCIP_Real oldub;
#endif
         for( i = 0; i < nvars; i++ )
         {
#ifdef SCIP_MORE_DEBUG
            oldlb = newlbs[i];
            oldub = newubs[i];
#endif
            idx = varinds[i];
            if( SCIPisPositive(scip, row1coefs[idx]) )
            {
               if( ninfs == 0 )
                  newlbs[i] = MAX(newlbs[i], (l1 + row1coefs[idx] * ubs[idx]) / row1coefs[idx]);
               else if( SCIPisInfinity(scip, ubs[idx]) )
                  newlbs[i] = MAX(newlbs[i], l1 / row1coefs[idx]);
            }
            else if ( SCIPisNegative(scip, row1coefs[idx]) )
            {
               if( ninfs == 0 )
                  newubs[i] = MIN(newubs[i], (l1 + row1coefs[idx] * lbs[idx]) / row1coefs[idx]);
               else if( SCIPisInfinity(scip, -lbs[idx]) )
                  newubs[i] = MIN(newubs[i], l1 / row1coefs[idx]);
            }
#ifdef SCIP_MORE_DEBUG
            if( !SCIPisEQ(scip, oldlb, newlbs[i]) || !SCIPisEQ(scip, oldub, newubs[i]) )
               SCIPdebugMsg(scip, "%g <= %g <= var_%i <= %g <= %g\n", oldlb, newlbs[i], i, newubs[i], oldub);
#endif
         }
      }

      /* update bound arrays and determine success */
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];

         assert(SCIPisLE(scip, lbs[idx], newlbs[i]));
         assert(SCIPisGE(scip, ubs[idx], newubs[i]));

         if( SCIPisGT(scip, newlbs[i], lbs[idx]) || SCIPisLT(scip, newubs[i], ubs[idx]) )
         {
            (*success) = TRUE;

            lbs[idx] = newlbs[i];
            ubs[idx] = newubs[i];
         }
      }
      SCIPfreeBufferArray(scip, &newubs);
      SCIPfreeBufferArray(scip, &newlbs);
   }

   SCIPfreeBufferArray(scip, &breakpoints);
   SCIPfreeBufferArray(scip, &signs);
   SCIPfreeBufferArray(scip, &varinds);
   SCIPfreeBufferArray(scip, &row2coefs);
   SCIPfreeBufferArray(scip, &row1coefs);

   return SCIP_OKAY;
}

/** get minimal and maximal residual activities without one specific column */
static
void getMinMaxActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real*            lbs,                /**< lower bounds */
   SCIP_Real*            ubs,                /**< upper bounds */
   SCIP_Real*            minresactivity,     /**< minimum residual activity of this row */
   SCIP_Real*            maxresactivity,     /**< maximum residual activity of this row */
   SCIP_Bool*            isminsettoinfinity, /**< flag indicating if minresactiviy is set to infinity */
   SCIP_Bool*            ismaxsettoinfinity  /**< flag indicating if maxresactiviy is set to infinity */
   )
{
   SCIP_Real coef;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int nmaxactneginf;
   int nmaxactposinf;
   int nminactneginf;
   int nminactposinf;
   SCIP_Real maxresact;
   SCIP_Real minresact;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   *isminsettoinfinity = FALSE;
   *ismaxsettoinfinity = FALSE;

   nmaxactneginf = 0;
   nmaxactposinf = 0;
   nminactneginf = 0;
   nminactposinf = 0;
   maxresact = 0;
   minresact = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for( ; rowpnt < rowend; rowpnt++, valpnt++ )
   {
      if(*rowpnt == col)
         continue;

      coef = *valpnt;

      /* positive coefficient */
      if( coef > 0.0 )
      {
         if( SCIPisInfinity(scip, ubs[col]) )
            nmaxactposinf++;
         else
            maxresact += coef * ubs[col];

         if( SCIPisInfinity(scip, -lbs[col]) )
            nminactneginf++;
         else
            minresact += coef * lbs[col];
      }
      else /* negative coefficient */
      {
         if( SCIPisInfinity(scip, -lbs[col]) )
            nmaxactneginf++;
         else
            maxresact += coef * lbs[col];

         if( SCIPisInfinity(scip, ubs[col]) )
            nminactposinf++;
         else
            minresact += coef * ubs[col];
      }
   }

   if( (nmaxactneginf + nmaxactposinf) > 0 )
      *ismaxsettoinfinity = TRUE;
   else
      *maxresactivity = maxresact;

   if( (nminactneginf + nminactposinf) > 0 )
      *isminsettoinfinity = TRUE;
   else
      *minresactivity = minresact;
}

/** calculate the upper and lower bound of one variable from one row */
static
void getVarBoundsOfRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index of variable */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this column in this row */
   SCIP_Real*            lbs,                /**< lower bounds */
   SCIP_Real*            ubs,                /**< upper bounds */
   SCIP_Real*            rowub,              /**< upper bound of row */
   SCIP_Bool*            ubfound,            /**< flag indicating that an upper bound was calculated */
   SCIP_Real*            rowlb,              /**< lower bound of row */
   SCIP_Bool*            lbfound             /**< flag indicating that a lower bound was caluclated */
   )
{
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(rowub != NULL);
   assert(ubfound != NULL);
   assert(rowlb != NULL);
   assert(lbfound != NULL);

   *rowub = SCIPinfinity(scip);
   *ubfound = FALSE;
   *rowlb = -SCIPinfinity(scip);
   *lbfound = FALSE;

   getMinMaxActivityResiduals(scip, matrix, col, row, lbs, ubs,
      &minresactivity, &maxresactivity,
      &isminsettoinfinity, &ismaxsettoinfinity);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);

   if( val > 0.0 )
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) )
      {
         *rowub = (rhs - minresactivity) / val; // maybe one wants some kind of numerical guard of check that values is not too small for all these
         *ubfound = TRUE;
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) )
      {
         *rowlb = (lhs - maxresactivity) / val;
         *lbfound = TRUE;
      }
   }
   else
   {
      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) )
      {
         *rowub = (lhs - maxresactivity) / val;
         *ubfound = TRUE;
      }

      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) )
      {
         *rowlb = (rhs - minresactivity) / val;
         *lbfound = TRUE;
      }
   }
}


/** detect implied variable bounds */
static
void getImpliedBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index for implied free test */
   SCIP_Real*            lbs,                /**< lower bounds */
   SCIP_Real*            ubs,                /**< upper bounds */
   SCIP_Bool*            ubimplied,          /**< flag indicating an implied upper bound */
   SCIP_Bool*            lbimplied           /**< flag indicating an implied lower bound */
   )
{
   SCIP_Real impliedub;
   SCIP_Real impliedlb;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbs != NULL);
   assert(ubs != NULL);
   assert(ubimplied != NULL);
   assert(lbimplied != NULL);

   *ubimplied = FALSE;
   impliedub = SCIPinfinity(scip);

   *lbimplied = FALSE;
   impliedlb = -SCIPinfinity(scip);

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);
   for( ; (colpnt < colend); colpnt++, valpnt++ )
   {
      SCIP_Real rowub;
      SCIP_Bool ubfound;
      SCIP_Real rowlb;
      SCIP_Bool lbfound;

      getVarBoundsOfRow(scip, matrix, col, *colpnt, *valpnt, lbs, ubs,
         &rowub, &ubfound, &rowlb, &lbfound);

      if( ubfound && (rowub < impliedub) )
         impliedub = rowub;

      if( lbfound && (rowlb > impliedlb) )
         impliedlb = rowlb;
   }

   /* we consider +/-inf bounds as implied bounds */
   if( SCIPisInfinity(scip, ubs[col]) ||
      (!SCIPisInfinity(scip, ubs[col]) && SCIPisLE(scip, impliedub, ubs[col])) )
      *ubimplied = TRUE;

   if( SCIPisInfinity(scip, -lbs[col]) ||
      (!SCIPisInfinity(scip, -lbs[col]) && SCIPisGE(scip, impliedlb, lbs[col])) )
      *lbimplied = TRUE;
}


/** calculate minimal column activity from one variable without one row */
static
SCIP_Real getMinColActWithoutRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   withoutrow,         /**< exclude this row index */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual              /**< upper bounds of dual variables */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   SCIP_Real mincolactivity;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual != NULL);
   assert(ubdual != NULL);

   mincolactivity = 0;

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; colpnt < colend; colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      if( row == withoutrow )
         continue;

      if( val > 0.0 )
      {
         assert(!SCIPisInfinity(scip, -lbdual[row]));
         mincolactivity += val * lbdual[row];
      }
      else if( val < 0.0 )
      {
         assert(!SCIPisInfinity(scip, ubdual[row]));
         mincolactivity += val * ubdual[row];
      }
   }

   return mincolactivity;
}


/** In the primal the residual activity of a constraint w.r.t. a variable is the activity of the constraint without the variable.
 *  This function does the same but in the dual.
 *  It computes the residual activity of column 'col' w.r.t. variable 'row'
 */
static
void calcMinColActResidual(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< matrix coefficient */
   SCIP_Real*            lbdual,             /**< lower bounds of the dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of the dual variables */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   int*                  mincolactinf,       /**< number of infinite contributions to minimal column activity */
   SCIP_Real*            mincolresact        /**< minimal residual column activity */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(mincolact != NULL);
   assert(mincolactinf != NULL);
   assert(mincolresact != NULL);

   *mincolresact = -SCIPinfinity(scip);

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, -lbdual[row]) )
      {
         assert(mincolactinf[col] >= 1);
         if( mincolactinf[col] == 1 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row, lbdual, ubdual);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( mincolactinf[col] > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * lbdual[row];
      }
   }
   else if( val < 0.0 )
   {
      if( SCIPisInfinity(scip, ubdual[row]) )
      {
         assert(mincolactinf[col] >= 1);
         if( mincolactinf[col] == 1 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row, lbdual, ubdual);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( mincolactinf[col] > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * ubdual[row];
      }
   }
}

/** calculate minimal column activity of one column */
static
void calcMinColActivity(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column for activity calculations */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   int*                  mincolactinf        /**< number of -inf contributions to minimal column activity */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(mincolact != NULL);
   assert(mincolactinf != NULL);

   /* init activities */
   mincolact[col] = 0.0;
   mincolactinf[col] = 0;

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   /* calculate column activities */
   for( ; colpnt < colend; colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      if( val > 0.0 )
      {
         if(SCIPisInfinity(scip, -lbdual[row]))
            mincolactinf[col]++;
         else
            mincolact[col] += val * lbdual[row];
      }
      else if( val < 0.0 )
      {
         if(SCIPisInfinity(scip, ubdual[row]))
            mincolactinf[col]++;
         else
            mincolact[col] += val * ubdual[row];
      }
   }

   /* update column activities if infinity counters are greater 0 */
   if( mincolactinf[col] > 0 )
      mincolact[col] = -SCIPinfinity(scip);
}

/** calculate maximal column activity of one column */
static
void calcMaxColActivity(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column for activity calculations */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   SCIP_Real*            maxcolact,          /**< minimal column activities */
   int*                  maxcolactinf        /**< number of -inf contributions to minimal column activity */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactinf != NULL);

   /* init activities */
   maxcolact[col] = 0.0;
   maxcolactinf[col] = 0;

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   /* calculate column activities */
   for( ; colpnt < colend; colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      if( val > 0.0 )
      {
         if(SCIPisInfinity(scip, ubdual[row]))
            maxcolactinf[col]++;
         else
            maxcolact[col] += val * ubdual[row];
      }
      else if( val < 0.0 )
      {
         if(SCIPisInfinity(scip, -lbdual[row]))
            maxcolactinf[col]++;
         else
            maxcolact[col] += val * lbdual[row];
      }
   }

   /* update column activities if infinity counters are greater 0 */
   if( maxcolactinf[col] > 0 )
      maxcolact[col] = SCIPinfinity(scip);
}


/** update minimal/maximal column activity infinity counters */
static
void infinityCountUpdate(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   SCIP_Bool*            isubimplied,        /**< flags indicating of the upper bound is implied */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   int*                  mincolactinf,       /**< number of infinity contributions to minimal column activity */
   SCIP_Bool             ubinfchange,        /**< flag indicating if the upper bound has changed from infinity to a finite value */
   SCIP_Bool             lbinfchange         /**< flag indicating if the lower bound has changed from -infinity to a finite value */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   SCIP_Real val;
   int col;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   /* look at all column entries present within row and update the
    * corresponding infinity counters. if one counter gets to zero,
    * then calculate this column activity new.
    */

   for(; (rowpnt < rowend); rowpnt++, valpnt++ )
   {
      col = *rowpnt;
      val = *valpnt;

      if( isubimplied[col] )
      {
         if( val < 0 )
         {
            if( ubinfchange )
            {
               assert(mincolactinf[col] > 0);
               mincolactinf[col]--;
            }
         }
         else if( val > 0 )
         {
            if( lbinfchange )
            {
               assert(mincolactinf[col] > 0);
               mincolactinf[col]--;
            }
         }

         if( mincolactinf[col] == 0 )
            calcMinColActivity(scip, matrix, col, lbdual, ubdual, mincolact, mincolactinf);
      }
   }
}

#ifdef SCIP_DEBUG
/** use LP calculations for determining the best dual variable bounds from a specific row index */
static
SCIP_RETCODE determineBestBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index for dual bound calculations */
   SCIP_Bool             solveLP,            /**< flag indicating to solve subscip LP */
   SCIP_Real*            lowerbnddual,       /**< lower bound of dual variable */
   SCIP_Real*            upperbnddual        /**< upper bound of dual variable */
   )
{
   int i;
   int nrows;
   int ncols;
   int numberconvars;
   SCIP_VAR* var;
   SCIP_VAR** variables;
   SCIP_VAR** tmpvars;
   SCIP_Real* tmpcoef;
   SCIP_CONS** constraints;
   int numDualVars;
   SCIP* subscip;
   SCIP_RETCODE retcode;
   char name[SCIP_MAXSTRLEN+3];
   int fillcnt;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;
   int* colmap;

   *lowerbnddual = -SCIPinfinity(scip);
   *upperbnddual = SCIPinfinity(scip);

   nrows = SCIPmatrixGetNRows(matrix);
   assert(0 <= row && row < nrows);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPcreateProbBasic(subscip, "subscip") );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* avoid recursive calls */
   SCIP_CALL( SCIPsetIntParam(subscip, "presolving/dualinfer/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", TRUE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &colmap, ncols) );
   numberconvars = 0;
   for(i = 0; i < ncols; i++)
   {
      var = SCIPmatrixGetVar(matrix, i);
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
      {
         colmap[i] = numberconvars; /* start numbering with 0 */
         numberconvars++;
      }
      else
         colmap[i] = -1;
   }
   numDualVars = nrows + 2 * numberconvars;

   /* create dual variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &variables, numDualVars) );
   for( i = 0; i < nrows; i++ )
   {
      variables[i] = NULL;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y%d", i);
      if( !SCIPmatrixIsRowRhsInfinity(matrix, i ) )
      {
         /* dual variable for equation or ranged row */
         SCIP_CALL( SCIPcreateVarBasic(subscip, &variables[i], name,
               -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      }
      else
      {
         /* dual variable for >= inequality */
         SCIP_CALL( SCIPcreateVarBasic(subscip, &variables[i], name,
               0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      }
      SCIP_CALL( SCIPaddVar(subscip, variables[i]) );
      assert( variables[i] != NULL );
   }

   /* in addition, we introduce dual variables for the bounds,
      because we treat each continuous variable as a free variable */
   fillcnt = nrows;
   for( i = 0; i < numberconvars; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ylb%d", fillcnt);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &variables[fillcnt], name,
            0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(subscip, variables[fillcnt]) );
      assert( variables[fillcnt] != NULL );
      fillcnt++;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yub%d", fillcnt);
      SCIP_CALL( SCIPcreateVarBasic(subscip, &variables[fillcnt], name,
            0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(subscip, variables[fillcnt]) );
      assert( variables[fillcnt] != NULL );
      fillcnt++;
   }
   assert(numDualVars == fillcnt);

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, numDualVars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcoef, numDualVars) );

   SCIP_CALL( SCIPallocBufferArray(scip, &constraints, numberconvars) );
   for( i = 0; i <numberconvars; i++)
      constraints[i] = NULL;

   for(i = 0; i < ncols; i++)
   {
      var = SCIPmatrixGetVar(matrix, i);
      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
      {
         SCIP_Real objval = SCIPvarGetObj(var);
         int cidx = colmap[i];
         assert(0 <= cidx && cidx < numberconvars);

         colpnt = SCIPmatrixGetColIdxPtr(matrix, i);
         colend = colpnt + SCIPmatrixGetColNNonzs(matrix, i);
         valpnt = SCIPmatrixGetColValPtr(matrix, i);
         fillcnt = 0;
         for( ; colpnt < colend; colpnt++, valpnt++ )
         {
            assert(0 <= *colpnt && *colpnt < nrows);
            assert(variables[*colpnt] != NULL);
            tmpvars[fillcnt] = variables[*colpnt];
            tmpcoef[fillcnt] = *valpnt;
            fillcnt++;
         }

         /* consider dual variable for a lower bound */
         if(SCIPisGT(scip, SCIPvarGetLbGlobal(var), -SCIPinfinity(scip)))
         {
            assert(variables[nrows + 2 * cidx] != NULL);
            tmpvars[fillcnt] = variables[nrows + 2 * cidx];
            tmpcoef[fillcnt] = 1.0;
            fillcnt++;
         }

         /* consider dual variable for an upper bound */
         if(SCIPisLT(scip, SCIPvarGetUbGlobal(var), SCIPinfinity(scip)))
         {
            assert(variables[nrows + 2 * cidx + 1] != NULL);
            tmpvars[fillcnt] = variables[nrows + 2 * cidx + 1];
            tmpcoef[fillcnt] = -1.0;
            fillcnt++;
         }

         /* because we treat the continuous columns as free variable,
            we need here an equality */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "c%d", cidx);
         SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &constraints[cidx], name,
               fillcnt, tmpvars, tmpcoef, objval, objval) );
         SCIP_CALL( SCIPaddCons(subscip, constraints[cidx]) );
      }
   }

   /* determine lower dual bound via a minimization problem */
   SCIP_CALL( SCIPsetObjsense(subscip,SCIP_OBJSENSE_MINIMIZE) );
   SCIP_CALL( SCIPchgVarObj(subscip, variables[row], 1.0) );
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dbg_min_%s.lp", SCIPvarGetName(variables[row]));
   SCIP_CALL( SCIPwriteOrigProblem(subscip, name, "lp", FALSE) );
   if( solveLP )
   {
      retcode = SCIPsolve(subscip);
      if( retcode != SCIP_OKAY )
         SCIPwarningMessage(scip, "Error subscip: <%d>\n", retcode);
      else
      {
         if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
         {
            SCIP_SOL* sol;
            SCIP_Bool feasible;
            sol = SCIPgetBestSol(subscip);
            SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, TRUE, TRUE) );

            if(feasible)
               *lowerbnddual = SCIPgetSolOrigObj(subscip, sol);
         }
      }
      SCIP_CALL( SCIPfreeTransform(subscip) );
   }

   /* determine upper dual bound via a maximization problem */
   SCIP_CALL( SCIPsetObjsense(subscip,SCIP_OBJSENSE_MAXIMIZE) );
   SCIP_CALL( SCIPchgVarObj(subscip, variables[row], 1.0) );
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dbg_max_%s.lp", SCIPvarGetName(variables[row]));
   SCIP_CALL( SCIPwriteOrigProblem(subscip, name, "lp", FALSE) );
   if( solveLP )
   {
      retcode = SCIPsolve(subscip);
      if( retcode != SCIP_OKAY )
         SCIPwarningMessage(scip, "Error subscip: <%d>\n", retcode);
      else
      {
         if( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
         {
            SCIP_SOL* sol;
            SCIP_Bool feasible;
            sol = SCIPgetBestSol(subscip);
            SCIP_CALL( SCIPcheckSolOrig(subscip, sol, &feasible, TRUE, TRUE) );

            if(feasible)
               *upperbnddual = SCIPgetSolOrigObj(subscip, sol);
         }
      }
      SCIP_CALL( SCIPfreeTransform(subscip) );
   }

   /* release variables and constraints */
   for( i = 0; i < numDualVars; i++ )
   {
      if(variables[i] != NULL)
         SCIP_CALL( SCIPreleaseVar(subscip, &variables[i]) );
   }
   for( i = 0; i < numberconvars; i++ )
   {
      if(constraints[i] != NULL)
         SCIP_CALL( SCIPreleaseCons(subscip, &constraints[i]) );
   }

   SCIPfreeBufferArray(scip, &constraints);
   SCIPfreeBufferArray(scip, &tmpcoef);
   SCIPfreeBufferArray(scip, &tmpvars);
   SCIPfreeBufferArray(scip, &variables);
   SCIP_CALL( SCIPfree(&subscip) );
   SCIPfreeBufferArray(scip, &colmap);

   return SCIP_OKAY;
}
#endif

/** update bounds of the dual variables */
static
void updateDualBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real             objval,             /**< objective function value */
   SCIP_Real             val,                /**< matrix coefficient */
   int                   row,                /**< row index */
   SCIP_Real             mincolresact,       /**< minimal column residual activity */
   SCIP_Real*            lbdual,             /**< dual lower bounds */
   SCIP_Real*            ubdual,             /**< dual upper bounds */
   int*                  boundchanges,       /**< counter for the number of bound changes */
   SCIP_Bool*            ubinfchange,        /**< flag indicating an upper bound change from infinite to finite */
   SCIP_Bool*            lbinfchange         /**< flag indicating a lower bound change from infinite to finite */
   )
{
   SCIP_Real newlbdual;
   SCIP_Real newubdual;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(boundchanges != NULL);
   assert(ubinfchange != NULL);
   assert(lbinfchange != NULL);

   *ubinfchange = FALSE;
   *lbinfchange = FALSE;

   if( !SCIPisInfinity(scip, -mincolresact) )
   {
      if( val > 0 )
      {
         newubdual = (objval - mincolresact) / val;

         if( newubdual < ubdual[row] )
         {
            /* accept the new upper bound only if the numerics are reliable */
            if( SCIPisLE(scip,lbdual[row],newubdual) )
            {
               if( SCIPisInfinity(scip, ubdual[row]) )
                  *ubinfchange = TRUE;

               ubdual[row] = newubdual;
               (*boundchanges)++;
            }
         }
      }
      else if( val < 0 )
      {
         newlbdual = (objval - mincolresact) / val;

         if( newlbdual > lbdual[row] )
         {
            /* accept the new lower bound only if the numerics are reliable */
            if( SCIPisLE(scip,newlbdual,ubdual[row]) )
            {
               if( SCIPisInfinity(scip, -lbdual[row]) )
                  *lbinfchange = TRUE;

               lbdual[row] = newlbdual;
               (*boundchanges)++;
            }
         }
      }
   }
}

/** dual bound strengthening */
static
SCIP_RETCODE dualBoundStrengthening(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data structure */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  npossiblefixings,   /**< number of possible fixings */
   SIDECHANGE*           sidestochange,      /**< array holding if this is an implied equality */
   int*                  npossiblesidechanges/**< number of possible equality changes */
   )
{
   SCIP_Real* lbdual;
   SCIP_Real* ubdual;
   SCIP_Real* mincolact;
   int* mincolactinf;
   SCIP_Real* maxcolact;
   int* maxcolactinf;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;
   int boundchanges;
   int loops;
   int i;
   int j;
   int k;
   int nrows;
   int ncols;
   SCIP_Bool* isubimplied;
   SCIP_Bool* islbimplied;
   SCIP_Real* tmplbs;
   SCIP_Real* tmpubs;
   SCIP_VAR* var;
   int* implubvars;
   int nimplubvars;

   SCIP_Longint maxhashes;
   int maxlen;
   int pospp;
   int listsizepp;
   int posmm;
   int listsizemm;
   int pospm;
   int listsizepm;
   int posmp;
   int listsizemp;

   int* hashlistpp;
   int* hashlistmm;
   int* hashlistpm;
   int* hashlistmp;

   int* colidxlistpp;
   int* colidxlistmm;
   int* colidxlistpm;
   int* colidxlistmp;

   int block1start;
   int block1end;
   int block2start;
   int block2end;

   SCIP_HASHSET* pairhashset;
   SCIP_Real* colvalptr;
   int* colidxptr;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(varstofix != NULL);
   assert(npossiblefixings != NULL);
   assert(sidestochange != NULL);
   assert(npossiblesidechanges != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &tmplbs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpubs, ncols) );
   for( i = 0; i < ncols; i++ )
   {
      var = SCIPmatrixGetVar(matrix, i);
      tmplbs[i] = SCIPvarGetLbLocal(var);
      tmpubs[i] = SCIPvarGetUbLocal(var);
   }

   /* verify which bounds of continuous variables are implied */
   SCIP_CALL( SCIPallocBufferArray(scip, &isubimplied, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &islbimplied, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &implubvars, ncols) );
   nimplubvars = 0;
   for( i = 0; i < ncols; i++ )
   {
      var = SCIPmatrixGetVar(matrix, i);

      if( SCIPmatrixUplockConflict(matrix, i) || SCIPmatrixDownlockConflict(matrix, i) ||
         (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(var) != SCIP_VARTYPE_IMPLINT) )
      {
         /* we don't care about integral variables or variables that have conflicting locks */
         isubimplied[i] = FALSE;
         islbimplied[i] = FALSE;
      }
      else
      {
         getImpliedBounds(scip, matrix, i, tmplbs, tmpubs, &(isubimplied[i]), &(islbimplied[i]));

         /* if a continuous variable has a not implied upper bound we can
          * not use this variable (column) for propagating dual bounds.
          * not implied lowers bound can usually be treated.
          */

         /* collect continuous variables with implied upper bound */
         if( isubimplied[i] )
         {
            implubvars[nimplubvars] = i;
            nimplubvars++;
         }

         /* reset implied bounds for further detections of other implied bounds */
         if( isubimplied[i] )
            tmpubs[i] = SCIPinfinity(scip);

         if( islbimplied[i] )
            tmplbs[i] = -SCIPinfinity(scip);
      }
   }

   /* initialize bounds of the dual variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &lbdual, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubdual, nrows) );
   for( i = 0; i < nrows; i++ )
   {
      if( !SCIPmatrixIsRowRhsInfinity(matrix, i) )
      {
         /* dual free variable for equation or ranged row */
         lbdual[i] = -SCIPinfinity(scip);
         ubdual[i] = SCIPinfinity(scip);
      }
      else
      {
         /* dual variable for >= inequality */
         lbdual[i] = 0.0;
         ubdual[i] = SCIPinfinity(scip);
      }
   }

   /* run convex combination on pairs of continuous variables (columns) using Belotti's algorithm */
   if( nimplubvars >= 2 && presoldata->usetwocolcombine )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistpp, ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistmm, ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistpm, ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistmp, ncols) );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colidxlistpp, ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colidxlistmm, ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colidxlistpm, ncols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colidxlistmp, ncols) );

      pospp = 0;
      posmm = 0;
      pospm = 0;
      posmp = 0;
      listsizepp = ncols;
      listsizemm = ncols;
      listsizepm = ncols;
      listsizemp = ncols;
      maxhashes = presoldata->maxhashfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)ncols) * presoldata->maxhashfac);

      for( i = 0; i < nimplubvars; i++)
      {
         if( ((SCIP_Longint)pospp) + posmm + pospm + posmp > maxhashes )
            break;

         colvalptr = SCIPmatrixGetColValPtr(matrix, implubvars[i]);
         colidxptr = SCIPmatrixGetColIdxPtr(matrix, implubvars[i]);
         maxlen = MIN(presoldata->maxconsiderednonzeros, SCIPmatrixGetColNNonzs(matrix, implubvars[i])); /*lint !e666*/
         for( j = 0; j < maxlen; j++)
         {
            for( k = j + 1; k < maxlen; k++)
            {
               if( SCIPisPositive(scip, colvalptr[j]) )
               {
                  if(SCIPisPositive(scip, colvalptr[k]) )
                  {
                     SCIP_CALL( addEntry(scip, &pospp, &listsizepp, &hashlistpp, &colidxlistpp,
                        hashIndexPair(colidxptr[j], colidxptr[k]), implubvars[i]) );
                  }
                  else
                  {
                     SCIP_CALL( addEntry(scip, &pospm, &listsizepm, &hashlistpm, &colidxlistpm,
                        hashIndexPair(colidxptr[j], colidxptr[k]), implubvars[i]) );
                  }
               }
               else
               {
                  if(SCIPisPositive(scip, colvalptr[k]) )
                  {
                     SCIP_CALL( addEntry(scip, &posmp, &listsizemp, &hashlistmp, &colidxlistmp,
                        hashIndexPair(colidxptr[j], colidxptr[k]), implubvars[i]) );
                  }
                  else
                  {
                     SCIP_CALL( addEntry(scip, &posmm, &listsizemm, &hashlistmm, &colidxlistmm,
                        hashIndexPair(colidxptr[j], colidxptr[k]), implubvars[i]) );
                  }
               }
            }
         }
      }
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(scip, "hashlist sizes: pp %d, mm %d, pm %d, mp %d \n", pospp, posmm, pospm, posmp);
#endif
      SCIPsortIntInt(hashlistpp, colidxlistpp, pospp);
      SCIPsortIntInt(hashlistmm, colidxlistmm, posmm);
      SCIPsortIntInt(hashlistpm, colidxlistpm, pospm);
      SCIPsortIntInt(hashlistmp, colidxlistmp, posmp);

      SCIP_CALL( SCIPhashsetCreate(&pairhashset, SCIPblkmem(scip), 1) );

      /* Process pp and mm lists */
      if( pospp > 0 && posmm > 0 )
      {
         SCIP_Longint ncombines;
         SCIP_Longint maxcombines;
         SCIP_Bool finished;
         SCIP_Bool success;
         int combinefails;
         int retrievefails;
         COLPAIR colpair;

         finished = FALSE;
         block1start = 0;
         block1end = 0;
         block2start = 0;
         block2end = 0;

         maxcombines = presoldata->maxpairfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)ncols) * presoldata->maxpairfac);

         ncombines = 0;
         combinefails = 0;
         retrievefails = 0;
         findNextBlock(hashlistpp, pospp, &block1start, &block1end);
         findNextBlock(hashlistmm, posmm, &block2start, &block2end);
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "processing pp and mm\n");
#endif

         // same as in the rworowbnd presolver - both while loops to basically the same with one using pp and mm and the other pm and mp
         // I would write an additional function and remove the code duplication
         while( !finished )
         {
            if( hashlistpp[block1start] == hashlistmm[block2start] )
            {
               for( i = block1start; i < block1end; i++ )
               {
                  for( j = block2start; j < block2end; j++ )
                  {
                     if( colidxlistpp[i] != colidxlistmm[j] )
                     {
                        colpair.col1idx = MIN(colidxlistpp[i], colidxlistmm[j]);
                        colpair.col2idx = MAX(colidxlistpp[i], colidxlistmm[j]);

                        if( !SCIPhashsetExists(pairhashset, encodeColPair(&colpair)) )
                        {
                           int* colpnt1       = SCIPmatrixGetColIdxPtr(matrix, colpair.col1idx);
                           SCIP_Real* valpnt1 = SCIPmatrixGetColValPtr(matrix, colpair.col1idx);
                           int* colpnt2       = SCIPmatrixGetColIdxPtr(matrix, colpair.col2idx);
                           SCIP_Real* valpnt2 = SCIPmatrixGetColValPtr(matrix, colpair.col2idx);
                           SCIP_Real obj1     = SCIPvarGetObj(SCIPmatrixGetVar(matrix, colpair.col1idx));
                           SCIP_Real obj2     = SCIPvarGetObj(SCIPmatrixGetVar(matrix, colpair.col2idx));
                           int collen1        = SCIPmatrixGetColNNonzs(matrix, colpair.col1idx);
                           int collen2        = SCIPmatrixGetColNNonzs(matrix, colpair.col2idx);

                           success = FALSE;

                           SCIP_CALL( combineCols(scip, colpnt1, colpnt2, valpnt1, valpnt2, obj1, obj2, collen1,
                                    collen2, nrows, TRUE, TRUE, lbdual, ubdual, &success) );

                           if( success )
                              combinefails = 0;
                           else
                              combinefails++;

                           SCIP_CALL( SCIPhashsetInsert(pairhashset, SCIPblkmem(scip), encodeColPair(&colpair)) );
                           ncombines++;
                           if( ncombines >= maxcombines || combinefails >= presoldata->maxcombinefails )
                              finished = TRUE;
#ifdef SCIP_MORE_DEBUG
                           SCIPdebugMsg(scip, "pm/mp: %d retrievefails before reset, %d combines\n", retrievefails, ncombines);
#endif
                           retrievefails = 0;
                        }
                        else if( retrievefails < presoldata->maxretrievefails )
                           retrievefails++;
                        else
                           finished = TRUE;
                     }
                     if( finished )
                        break;
                  }
                  if( finished )
                     break;
               }

               if( block1end < pospp && block2end < posmm )
               {
                  findNextBlock(hashlistpp, pospp, &block1start, &block1end);
                  findNextBlock(hashlistmm, posmm, &block2start, &block2end);
               }
               else
                  finished = TRUE;
            }
            else if( hashlistpp[block1start] < hashlistmm[block2start] && block1end < pospp )
               findNextBlock(hashlistpp, pospp, &block1start, &block1end);
            else if( hashlistpp[block1start] > hashlistmm[block2start] && block2end < posmm )
               findNextBlock(hashlistmm, posmm, &block2start, &block2end);
            else
               finished = TRUE;
         }
      }

      /* Process pm and mp lists */
      if( pospm > 0 && posmp > 0 )
      {
         SCIP_Longint maxcombines;
         SCIP_Longint ncombines;
         SCIP_Bool finished;
         SCIP_Bool success;
         int combinefails;
         int retrievefails;
         COLPAIR colpair;

         finished = FALSE;
         block1start = 0;
         block1end = 0;
         block2start = 0;
         block2end = 0;

         maxcombines = presoldata->maxpairfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)ncols) * presoldata->maxpairfac);

         ncombines = 0;
         combinefails = 0;
         retrievefails = 0;
         findNextBlock(hashlistpm, pospm, &block1start, &block1end);
         findNextBlock(hashlistmp, posmp, &block2start, &block2end);
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "processing pm and mp\n");
#endif

         while( !finished )
         {
            if( hashlistpm[block1start] == hashlistmp[block2start] )
            {
               for( i = block1start; i < block1end; i++ )
               {
                  for( j = block2start; j < block2end; j++ )
                  {
                     if( colidxlistpm[i] != colidxlistmp[j] )
                     {
                        colpair.col1idx = MIN(colidxlistpm[i], colidxlistmp[j]);
                        colpair.col2idx = MAX(colidxlistpm[i], colidxlistmp[j]);

                        if( !SCIPhashsetExists(pairhashset, encodeColPair(&colpair)) )
                        {
                           int* colpnt1       = SCIPmatrixGetColIdxPtr(matrix, colpair.col1idx);
                           SCIP_Real* valpnt1 = SCIPmatrixGetColValPtr(matrix, colpair.col1idx);
                           int* colpnt2       = SCIPmatrixGetColIdxPtr(matrix, colpair.col2idx);
                           SCIP_Real* valpnt2 = SCIPmatrixGetColValPtr(matrix, colpair.col2idx);
                           SCIP_Real obj1     = SCIPvarGetObj(SCIPmatrixGetVar(matrix, colpair.col1idx));
                           SCIP_Real obj2     = SCIPvarGetObj(SCIPmatrixGetVar(matrix, colpair.col2idx));
                           int collen1        = SCIPmatrixGetColNNonzs(matrix, colpair.col1idx);
                           int collen2        = SCIPmatrixGetColNNonzs(matrix, colpair.col2idx);

                           success = FALSE;

                           SCIP_CALL( combineCols(scip, colpnt1, colpnt2, valpnt1, valpnt2, obj1, obj2, collen1,
                                    collen2, nrows, TRUE, TRUE, lbdual, ubdual, &success) );

                           if( success )
                              combinefails = 0;
                           else
                              combinefails++;

                           SCIP_CALL( SCIPhashsetInsert(pairhashset, SCIPblkmem(scip), encodeColPair(&colpair)) );
                           ncombines++;
                           if( ncombines >= maxcombines || combinefails >= presoldata->maxcombinefails )
                              finished = TRUE;

                           retrievefails = 0;
                        }
                        else if( retrievefails < presoldata->maxretrievefails )
                           retrievefails++;
                        else
                           finished = TRUE;
                     }
                     if( finished )
                        break;
                  }
                  if( finished )
                     break;
               }

               if( block1end < pospm && block2end < posmp )
               {
                  findNextBlock(hashlistpm, pospm, &block1start, &block1end);
                  findNextBlock(hashlistmp, posmp, &block2start, &block2end);
               }
               else
                  finished = TRUE;
            }
            else if( hashlistpm[block1start] < hashlistmp[block2start] && block1end < pospm )
               findNextBlock(hashlistpm, pospm, &block1start, &block1end);
            else if( hashlistpm[block1start] > hashlistmp[block2start] && block2end < posmp )
               findNextBlock(hashlistmp, posmp, &block2start, &block2end);
            else
               finished = TRUE;
         }
      }

      SCIPhashsetFree(&pairhashset, SCIPblkmem(scip));
      SCIPfreeBlockMemoryArray(scip, &colidxlistmp, listsizemp);
      SCIPfreeBlockMemoryArray(scip, &colidxlistpm, listsizepm);
      SCIPfreeBlockMemoryArray(scip, &colidxlistmm, listsizemm);
      SCIPfreeBlockMemoryArray(scip, &colidxlistpp, listsizepp);
      SCIPfreeBlockMemoryArray(scip, &hashlistmp, listsizemp);
      SCIPfreeBlockMemoryArray(scip, &hashlistpm, listsizepm);
      SCIPfreeBlockMemoryArray(scip, &hashlistmm, listsizemm);
      SCIPfreeBlockMemoryArray(scip, &hashlistpp, listsizepp);

#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(scip, "CombCols:\n");
      for( i = 0; i < nrows; i++ )
      {
         assert(SCIPisLE(scip,lbdual[i],ubdual[i]));
         SCIPdebugMsg(scip, "y%d=[%g,%g]\n",i,lbdual[i],ubdual[i]);
      }
      SCIPdebugMsg(scip,"\n");
#endif
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &mincolact, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincolactinf, ncols) );

   /* apply dual bound strengthening */
   loops = 0;
   boundchanges = 1;
   while( 0 < boundchanges && loops < presoldata->maxdualbndloops )
   {
      loops++;
      boundchanges = 0;

      for( i = 0; i < nimplubvars; i++ )
      {
         assert(SCIPvarGetType(SCIPmatrixGetVar(matrix, implubvars[i])) == SCIP_VARTYPE_CONTINUOUS ||
            SCIPvarGetType(SCIPmatrixGetVar(matrix, implubvars[i])) == SCIP_VARTYPE_IMPLINT);
         calcMinColActivity(scip, matrix, implubvars[i], lbdual, ubdual, mincolact, mincolactinf);
      }

      for( i = 0; i < nimplubvars; i++ )
      {
         SCIP_Real objval;
         SCIP_Bool ubinfchange;
         SCIP_Bool lbinfchange;
         int col;

         col = implubvars[i];
         var = SCIPmatrixGetVar(matrix, col);

         objval = SCIPvarGetObj(var);
         colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
         colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
         valpnt = SCIPmatrixGetColValPtr(matrix, col);

         for( ; colpnt < colend; colpnt++, valpnt++ )
         {
            int row;
            SCIP_Real val;
            SCIP_Real mincolresact;

            row = *colpnt;
            val = *valpnt;

            calcMinColActResidual(scip, matrix, col, row, val, lbdual, ubdual,
               mincolact, mincolactinf, &mincolresact);

            updateDualBounds(scip, matrix, objval, val, row, mincolresact,
               lbdual, ubdual, &boundchanges, &ubinfchange, &lbinfchange);

            if( ubinfchange || lbinfchange )
               infinityCountUpdate(scip, matrix, row, lbdual, ubdual, isubimplied,
                  mincolact, mincolactinf, ubinfchange, lbinfchange);
         }
      }
   }

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(scip, "BndStr:\n");
   for( i = 0; i < nrows; i++ )
   {
      assert(SCIPisLE(scip,lbdual[i],ubdual[i]));
      SCIPdebugMsg(scip, "y%d=[%g,%g]\n",i,lbdual[i],ubdual[i]);
   }
   SCIPdebugMsg(scip,"\n");
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolact, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolactinf, ncols) );

   /* calculate final minimal and maximal column activities */
   for( i = 0; i < ncols; i++ )
   {
      calcMinColActivity(scip, matrix, i, lbdual, ubdual, mincolact, mincolactinf);
      calcMaxColActivity(scip, matrix, i, lbdual, ubdual, maxcolact, maxcolactinf);
   }

   for( i = 0; i < ncols; i++ )
   {
      SCIP_Real objval;

      var = SCIPmatrixGetVar(matrix, i);

      /* do not fix variables if the locks do not match */
      if( SCIPmatrixUplockConflict(matrix, i) || SCIPmatrixDownlockConflict(matrix, i) )
         continue;

      objval = SCIPvarGetObj(var);

      /* c_j - sup(y^T A_{.j}) > 0 => fix x_j to its lower bound */
      if( SCIPisGT(scip, objval, maxcolact[i]) && varstofix[i] == NOFIX )
      {
         if( SCIPisGT(scip, SCIPvarGetLbGlobal(var), -SCIPinfinity(scip)) )
         {
            varstofix[i] = FIXATLB;
            (*npossiblefixings)++;
         }
      }

      /* c_j - inf(y^T A_{.j}) < 0 => fix x_j to its upper bound */
      if( SCIPisLT(scip, objval, mincolact[i]) && varstofix[i] == NOFIX )
      {
         if( SCIPisLT(scip, SCIPvarGetUbGlobal(var), SCIPinfinity(scip)) )
         {
            varstofix[i] = FIXATUB;
            (*npossiblefixings)++;
         }
      }
   }

   for( i = 0; i < nrows; i++ )
   {
      /* implied equality: y_i > 0 =>  A_{i.}x - b_i = 0 */
      if( SCIPmatrixIsRowRhsInfinity(matrix, i) )
      {
         if( SCIPisGT(scip, lbdual[i], 0.0) && (sidestochange[i] == NOCHANGE) )
         {
            /* change >= inequality to equality */
            sidestochange[i] = RHSTOLHS;
            (*npossiblesidechanges)++;
         }
      }
      else
      {
         if( !SCIPmatrixIsRowRhsInfinity(matrix, i) &&
            !SCIPisEQ(scip,SCIPmatrixGetRowLhs(matrix, i),SCIPmatrixGetRowRhs(matrix, i)) )
         {
            /* for ranged rows we have to decide which side (lhs or rhs) determines the equality */
            if( SCIPisGT(scip, lbdual[i], 0.0) && sidestochange[i]==NOCHANGE )
            {
               sidestochange[i] = RHSTOLHS;
               (*npossiblesidechanges)++;
            }

            if( SCIPisLT(scip, ubdual[i], 0.0) && sidestochange[i]==NOCHANGE)
            {
               sidestochange[i] = LHSTORHS;
               (*npossiblesidechanges)++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &maxcolactinf);
   SCIPfreeBufferArray(scip, &maxcolact);
   SCIPfreeBufferArray(scip, &mincolactinf);
   SCIPfreeBufferArray(scip, &mincolact);

   SCIPfreeBufferArray(scip, &ubdual);
   SCIPfreeBufferArray(scip, &lbdual);
   SCIPfreeBufferArray(scip, &implubvars);
   SCIPfreeBufferArray(scip, &islbimplied);
   SCIPfreeBufferArray(scip, &isubimplied);
   SCIPfreeBufferArray(scip, &tmpubs);
   SCIPfreeBufferArray(scip, &tmplbs);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualinfer)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDualinfer(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeDualinfer)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualinfer)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_Bool infeasible;
   SCIP_PRESOLDATA* presoldata;
   FIXINGDIRECTION* varstofix;
   int npossiblefixings;
   int nconvarsfixed;
   int nintvarsfixed;
   int nbinvarsfixed;
   SIDECHANGE* sidestochange;
   int npossiblesidechanges;
   int nsideschanged;
   int i;
   int nrows;
   int ncols;
   SCIP_VAR* var;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNContVars(scip) + SCIPgetNImplVars(scip) == 0 )
      return SCIP_OKAY;

   /* the reductions made in this presolver apply to all optimal solutions because of complementary slackness */
   if( !SCIPallowWeakDualReds(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, TRUE, &initialized, &complete, &infeasible,
      naddconss, ndelconss, nchgcoefs, nchgbds, nfixedvars) );

   /* if infeasibility was detected during matrix creation, return here */
   if( infeasible )
   {
      if( initialized )
         SCIPmatrixFree(scip, &matrix);

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( !initialized )
      return SCIP_OKAY;

   npossiblefixings = 0;
   nconvarsfixed = 0;
   nintvarsfixed = 0;
   nbinvarsfixed = 0;
   npossiblesidechanges = 0;
   nsideschanged = 0;

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sidestochange, nrows) );

   BMSclearMemoryArray(varstofix, ncols);
   BMSclearMemoryArray(sidestochange, nrows);

   SCIP_CALL( dualBoundStrengthening(scip, matrix, presoldata,
         varstofix, &npossiblefixings, sidestochange, &npossiblesidechanges) );

   if( npossiblefixings > 0 )
   {
      for( i = ncols - 1; i >= 0; --i )
      {
         SCIP_Bool fixed;

         var = SCIPmatrixGetVar(matrix, i);

         /* there should be no fixings for variables with inconsistent locks */
         assert(varstofix[i] == NOFIX || (!SCIPmatrixUplockConflict(matrix, i) && !SCIPmatrixDownlockConflict(matrix, i)));

         fixed = FALSE;

         if( varstofix[i] == FIXATLB )
         {
            SCIP_Real lb;
            lb = SCIPvarGetLbLocal(var);

            /* fix at lower bound */
            SCIP_CALL( SCIPfixVar(scip, var, lb, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, " -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               break;
            }
            assert(fixed);
            (*nfixedvars)++;
            *result = SCIP_SUCCESS;
         }
         else if( varstofix[i] == FIXATUB )
         {
            SCIP_Real ub;
            ub = SCIPvarGetUbLocal(var);

            /* fix at upper bound */
            SCIP_CALL( SCIPfixVar(scip, var, ub, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, " -> infeasible fixing\n");
               *result = SCIP_CUTOFF;
               break;
            }
            assert(fixed);
            (*nfixedvars)++;
            *result = SCIP_SUCCESS;
         }

         /* keep a small statistic which types of variables are fixed */
         if( fixed )
         {
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
               nconvarsfixed++;
            else if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
               nbinvarsfixed++;
            else
               nintvarsfixed++;
         }
      }
   }

   if( npossiblesidechanges > 0 )
   {
      for( i = 0; i < nrows; i++ )
      {
         SCIP_CONS* cons;
         SCIP_CONSHDLR* conshdlr;
         const char* conshdlrname;

         if( sidestochange[i] == NOCHANGE )
            continue;

         if( presoldata->maxrowsupport < SCIPmatrixGetRowNNonzs(matrix, i) )
            continue;

         cons = SCIPmatrixGetCons(matrix,i);
         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrname = SCIPconshdlrGetName(conshdlr);

         if( strcmp(conshdlrname, "linear") == 0 )
         {
            SCIP_Real lhs;
            SCIP_Real rhs;
            SCIP_Real matrixlhs;
            SCIP_Real matrixrhs;

            lhs = SCIPgetLhsLinear(scip, cons);
            rhs = SCIPgetRhsLinear(scip, cons);
            matrixlhs = SCIPmatrixGetRowLhs(matrix, i);
            matrixrhs = SCIPmatrixGetRowRhs(matrix, i);

            assert(!SCIPisEQ(scip, matrixlhs, matrixrhs));

            /* when creating the matrix, contraints are multiplied if necessary by (-1)
               * to ensure that the following representation is obtained:
               * infty >= a x >= b
               * or
               * c >= ax >= b (ranged rows)
               */

            /* for ranged constraints we have to distinguish between both sides */
            if( sidestochange[i] == RHSTOLHS )
            {
               if( SCIPisEQ(scip, matrixlhs, lhs) )
               {
                  /* change rhs to lhs */
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, matrixlhs) );
               }
               else
               {
                  /* consider multiplication by (-1) in the matrix */
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, -matrixlhs) );
               }

               nsideschanged++;
               (*nchgsides)++;
            }
            else if( sidestochange[i] == LHSTORHS )
            {
               if( SCIPisEQ(scip, matrixrhs, rhs) )
               {
                  /* change lhs to rhs */
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, matrixrhs) );
               }
               else
               {
                  /* consider multiplication by (-1) in the matrix */
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, -matrixrhs) );
               }

               nsideschanged++;
               (*nchgsides)++;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &sidestochange);
   SCIPfreeBufferArray(scip, &varstofix);

   if( (nconvarsfixed + nintvarsfixed + nbinvarsfixed) > 0 || npossiblesidechanges > 0 )
   {
      SCIPdebugMsg(scip, "### fixed vars [cont: %d, int: %d, bin: %d], changed sides [%d]\n",
         nconvarsfixed, nintvarsfixed, nbinvarsfixed, nsideschanged);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the dual inference presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualinfer(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;
   SCIP_PRESOLDATA* presoldata;

   /* create presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDualinfer, presoldata) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDualinfer) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeDualinfer) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/dualinfer/twocolcombine",
         "use convex combination of columns for determining dual bounds",
         &presoldata->usetwocolcombine, FALSE, DEFAULT_TWOCOLUMN_COMBINE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualinfer/maxdualbndloops",
         "maximal number of dual bound strengthening loops",
         &presoldata->maxdualbndloops, FALSE, DEFAULT_MAXLOOPS_DUALBNDSTR, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualinfer/maxconsiderednonzeros",
         "maximal number of considered non-zeros within one column (-1: no limit)",
         &presoldata->maxconsiderednonzeros, TRUE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualinfer/maxretrievefails",
         "maximal number of consecutive useless hashtable retrieves",
         &presoldata->maxretrievefails, TRUE, DEFAULT_MAXRETRIEVEFAILS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualinfer/maxcombinefails",
         "maximal number of consecutive useless column combines",
         &presoldata->maxcombinefails, TRUE, DEFAULT_MAXCOMBINEFAILS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualinfer/maxhashfac",
         "Maximum number of hashlist entries as multiple of number of columns in the problem (-1: no limit)",
         &presoldata->maxhashfac, TRUE, DEFAULT_MAXHASHFAC, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualinfer/maxpairfac",
         "Maximum number of processed column pairs as multiple of the number of columns in the problem (-1: no limit)",
         &presoldata->maxpairfac, TRUE, DEFAULT_MAXPAIRFAC, -1, INT_MAX, NULL, NULL) );

    SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualinfer/maxrowsupport",
         "Maximum number of row's non-zeros for changing inequality to equality",
         &presoldata->maxrowsupport, FALSE, DEFAULT_MAXROWSUPPORT, 2, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
