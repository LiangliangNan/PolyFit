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

/**@file   presol_tworowbnd.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  do bound tightening by using two rows
 * @author Dieter Weninger
 * @author Patrick Gemander
 *
 * Perform bound tightening on two inequalities with some common variables.
 * Two possible methods are being used.
 *
 * 1. LP-bound
 * Let two constraints be given:
 * \f{eqnarray*}{
 *   A_{iR} x_R + A_{iS} x_S              \geq b_i\\
 *   A_{kR} x_R              + A_{kT} x_T \geq b_k
 * \f}
 * with \f$N\f$ the set of variable indexes, \f$R \subseteq N\f$, \f$S \subseteq N\f$, \f$T \subseteq N\f$,
 * \f$R \cap S = \emptyset\f$, \f$R \cap T = \emptyset\f$, \f$S \cap T = \emptyset\f$ and row indices \f$i \not= k\f$.
 *
 * Let \f$\ell\f$ and \f$u\f$ be bound vectors for x and solve the following two LPs
 * \f{eqnarray*}{
 *   L = \min \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i, \ell \leq x \leq u \}\\
 *   U = \max \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i, \ell \leq x \leq u \}
 * \f}
 * and use \f$L\f$ and \f$U\f$ for getting bounds on \f$x_T\f$.
 *
 * If \f$L + \mbox{infimum}(A_{kT}x_T) \geq b_k\f$, then the second constraint above is redundant.
 *
 * More details can be found in
 * - Chen W. et. al "Two-row and two-column mixed-integer presolve using hashing-based pairing methods"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*
 * Additional debug defines in this presolver
 * SCIP_DEBUG_HASHING
 * SCIP_DEBUG_BOUNDS
 * SCIP_DEBUG_SINGLEROWLP
 */

#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_matrix.h"
#include "scip/presol_tworowbnd.h"
#include <string.h>

#define PRESOL_NAME                    "tworowbnd"
#define PRESOL_DESC                    "do bound tigthening by using two rows"
#define PRESOL_PRIORITY                -2000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS               0        /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING                  SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_ENABLECOPY             TRUE     /**< should tworowbnd presolver be copied to sub-SCIPs? */
#define DEFAULT_MAXCONSIDEREDNONZEROS  100      /**< maximal number of considered non-zeros within one row (-1: no limit) */
#define DEFAULT_MAXRETRIEVEFAILS       1000     /**< maximal number of consecutive useless hashtable retrieves */
#define DEFAULT_MAXCOMBINEFAILS        1000     /**< maximal number of consecutive useless row combines */
#define DEFAULT_MAXHASHFAC             10       /**< maximal number of hashlist entries as multiple of number of rows in the problem (-1: no limit) */
#define DEFAULT_MAXPAIRFAC             1        /**< maximal number of processed row pairs as multiple of the number of rows in the problem (-1: no limit) */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int maxpairfac;            /**< maximal number of processed row pairs as multiple of the number of rows in the problem (-1: no limit) */
   int maxhashfac;            /**< maximal number of hashlist entries as multiple of number of rows in the problem (-1: no limit) */
   int maxretrievefails;      /**< maximal number of consecutive useless hashtable retrieves */
   int maxcombinefails;       /**< maximal number of consecutive useless row combines */
   int maxconsiderednonzeros; /**< maximal number of considered non-zeros within one row (-1: no limit) */
   int nchgbnds;              /**< number of variable bounds changed by this presolver */
   int nuselessruns;          /**< number of runs where this presolver did not apply any changes */
   SCIP_Bool enablecopy;      /**< should tworowbnd presolver be copied to sub-SCIPs? */
};

/** structure representing a pair of row indices; used for lookup in a hashtable */
struct RowPair
{
   int row1idx;               /**< first row index */
   int row2idx;               /**< second row index */
};

typedef struct RowPair ROWPAIR;


/*
 * Local methods
 */

/** encode contents of a rowpair as void* pointer */
static
void* encodeRowPair(
   ROWPAIR*              rowpair             /**< pointer to rowpair */
   )
{
  uint64_t a = (uint64_t)(long)rowpair->row1idx;
  uint64_t b = (uint64_t)(long)rowpair->row2idx;
   return (void*)((a << 32) | b);
}

/** compute single positive int hashvalue for two ints */
static
int hashIndexPair(
   int                   idx1,               /**< first integer index */
   int                   idx2                /**< second integer index */
   )
{
   uint32_t hash = SCIPhashTwo(idx1, idx2);
   return (int)(hash >> 1);
}

/** add hash/rowidx pair to hashlist/rowidxlist */
static
SCIP_RETCODE addEntry(
   SCIP*                 scip,               /**< SCIP datastructure */
   int*                  pos,                /**< position of last entry added */
   int*                  listsize,           /**< size of hashlist and rowidxlist */
   int**                 hashlist,           /**< block memory array containing hashes */
   int**                 rowidxlist,         /**< block memory array containing row indices */
   int                   hash,               /**< hash to be inserted */
   int                   rowidx              /**< row index to be inserted */
   )
{
   if( (*pos) >= (*listsize) )
   {
      int newsize  = SCIPcalcMemGrowSize(scip, (*pos) + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, hashlist, (*listsize), newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, rowidxlist, (*listsize), newsize) );
      (*listsize) = newsize;
   }

   (*hashlist)[(*pos)] = hash;
   (*rowidxlist)[(*pos)] = rowidx;
   (*pos)++;

   return SCIP_OKAY;
}

/*  Within a sorted list, get next block with same value
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

/*  Solve single-row LP of the form
 *  min c^T x
 *  s.t. a^T x >= b
 *  lbs <= x <= ubs
 *
 *  First, the problem is transformed such that
 *  SCIPselectWeightedReal() can be applied, which
 *  then solves the problem as a continuous knapsack
 *  in linear time.
 */
static
SCIP_RETCODE solveSingleRowLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            a,                  /**< constraint coefficients */
   SCIP_Real             b,                  /**< right hand side */
   SCIP_Real*            c,                  /**< objective coefficients */
   SCIP_Real*            lbs,                /**< lower variable bounds */
   SCIP_Real*            ubs,                /**< upper variable bounds */
   int                   len,                /**< length of arrays */
   SCIP_Real*            obj,                /**< objective value of solution */
   SCIP_Bool*            solvable            /**< status whether LP was solvable */
   )
{
   int i;
   int k;
   int nvars;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real mincost;
   SCIP_Real maxgain;

#ifdef SCIP_DEBUG_SINGLEROWLP
   SCIPdebugMsg(scip, "solving single row LP with %d variables\n", len);
#endif

   nvars = 0;
   (*obj) = 0;
   (*solvable) = TRUE;
   mincost = SCIPinfinity(scip);
   maxgain = 0;
   for( i = 0; i < len; i++)
   {
      /* Handle variables with zero weight */
      if( SCIPisZero(scip, a[i]) )
      {
         /* a[i] = 0, c[i] > 0 */
         if( SCIPisPositive(scip, c[i]) )
         {
            if( SCIPisInfinity(scip, -lbs[i]) )
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
            else
               (*obj) += c[i] * lbs[i];
         }
         /* a[i] = 0, c[i] < 0 */
         else if( SCIPisNegative(scip, c[i]) )
         {
            if( SCIPisInfinity(scip, ubs[i]) )
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
            else
               (*obj) += c[i] * ubs[i];
         }
         /* Note that variables with a[i] = 0, c[i] = 0 can be ignored */
         continue;
      }

      /* Handle free variables */
      if( SCIPisInfinity(scip, -lbs[i]) && SCIPisInfinity(scip, ubs[i]) )
      {
         /* The problem is unbounded */
         if( (SCIPisPositive(scip, c[i]) && SCIPisNegative(scip, a[i])) ||
             (SCIPisNegative(scip, c[i]) && SCIPisPositive(scip, a[i])) )
         {
            (*solvable) = FALSE;
            return SCIP_OKAY;
         }
         else
         {
            mincost = MIN(mincost, c[i] / a[i]);
            maxgain = MAX(maxgain, c[i] / a[i]);
         }
         continue;
      }

      /* Swap variable orientation if lower bound is infinite */
      if( SCIPisInfinity(scip, -lbs[i]) )
      {
         c[i] = -c[i];
         a[i] = -a[i];
         lb = -ubs[i];
         ub = -lbs[i];
      }
      else
      {
         lb = lbs[i];
         ub = ubs[i];
      }

      /* Handle variables with infinite upper bound */
      if( SCIPisInfinity(scip, ub) )
      {
         if( SCIPisPositive(scip, a[i]) )
         {
            /* a[i] > 0, c[i] >= 0 */
            if( !SCIPisNegative(scip, c[i]) )
            {
               mincost = MIN(mincost, c[i]/a[i]);
            }
            /* a[i] > 0, c[i] < 0 */
            else
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
         }
         /* a[i] < 0, c[i] < 0 */
         else if( SCIPisNegative(scip, c[i]) )
         {
            maxgain = MAX(maxgain, c[i] / a[i]);
         }
         /* a[i] < 0, c[i] >= 0 results in dual fixing of this variable, which is included in the bound shift below */

         /* Shift lower bound to zero */
         if( !SCIPisZero(scip, lb) )
         {
            (*obj) += c[i] * lb;
            b -= a[i] * lb;
         }
         continue;
      }

      /* Handle fixed variables */
      if( SCIPisEQ(scip, lb, ub) )
      {
         (*obj) += c[i] * lb;
         b -= a[i] * lb;
         continue;
      }

      /* Dual fixing for variables with finite bounds */
      if( !SCIPisNegative(scip, c[i]) && SCIPisNegative(scip, a[i]) )
      {
         (*obj) += c[i] * lb;
         b -= a[i] * lb;
         continue;
      }
      else if( !SCIPisPositive(scip, c[i]) && SCIPisPositive(scip, a[i]) )
      {
         (*obj) += c[i] * ub;
         b -= a[i] * ub;
         continue;
      }

      assert(!SCIPisInfinity(scip, -lb));
      assert(!SCIPisInfinity(scip, ub));

      /* At this point the variable has finite bounds and a[i],c[i] are both positive or both negative.
       * Normalize variable such that
       *  1. x_i \in [0,1]
       *  2. a[i] > 0
       *  3. c[i] >= 0
       * and calculate its "unit price" c[i]/a[i].
       */
      if( SCIPisNegative(scip, a[i]) )
      {
         c[i] = -c[i];
         a[i] = -a[i];
         lb = -ubs[i];
         ub = -lbs[i];
      }

      /* All variables with a <= 0 have been handled and variables with a[i] = 0, c[i] = 0 ignored */
      assert(SCIPisPositive(scip, a[i]) && SCIPisPositive(scip, c[i]));

      /* Adjust objective offset and b to shift lower bound to zero */
      (*obj) += c[i] * lb;
      b -= a[i] * lb;

      /* Calculate unit price */
      c[nvars] = c[i] / a[i];

      /* Normalize bound [0, ub] to [0,1] */
      a[nvars] = (ub - lb) * a[i];
      nvars++;
   }

#ifdef SCIP_DEBUG_SINGLEROWLP
   SCIPdebugMsg(scip, "After preprocessing: obj = %g, b = %g, nvars = %d, mincost = %g, maxgain = %g\n", (*obj), b, nvars, mincost, maxgain);
#endif

   /* Actual solving starts here.
    * If maxgain > 0 holds, we have a variable that can relax the constraint to an arbitrary degree while yielding
    * a certain profit per unit. This will be called downslack. If mincost < inf holds, we have a variable that can
    * always satisfy the constraint at a certain unit price. This will be called upslack.
    */

   /* Problem is unbounded since the downslack variable yields higher gains than the upslack variable costs */
   if( SCIPisLT(scip, mincost, maxgain) )
   {
      (*solvable) = FALSE;
      return SCIP_OKAY;
   }
   /* Solution is trivial as we have slack variables of equal price for both directions */
   else if( SCIPisEQ(scip, mincost, maxgain) )
   {
      /* Use all elements with cost smaller than maxgain */
      for( i = 0; i < nvars; i++ )
      {
         if( SCIPisLT(scip, c[i], maxgain) )
         {
            (*obj) += c[i] * a[i];
            b -= a[i];
         }
      }
      /* Use slack variable to satisfy constraint */
      (*obj) += mincost * b;
      return SCIP_OKAY;
   }
   /* mincost > maxgain
    * In this case we need to solve the problem for the remaining variables with mincost > c[i] > maxgain.
    */
   else
   {
      /* Only keep variables that are cheaper than the upslack variable */
      if( !SCIPisInfinity(scip, mincost) )
      {
         k = 0;
         for( i = 0; i < nvars; i++ )
         {
            if( SCIPisLT(scip, c[i], mincost) )
            {
               c[k] = c[i];
               a[k] = a[i];
               k++;
            }
         }
         nvars = k;
      }

      /* Exploit all variables that are cheaper than the downslack variable */
      if( !SCIPisZero(scip, maxgain) )
      {
         k = 0;
         for( i = 0; i < nvars; i++ )
         {
            if( SCIPisLE(scip, c[i], maxgain) )
            {
               (*obj) += c[i] * a[i];
               b -= a[i];
            }
            else
            {
               c[k] = c[i];
               a[k] = a[i];
               k++;
            }
         }
         if( !SCIPisPositive(scip, b) )
         {
            (*obj) += maxgain * b;
            return SCIP_OKAY;
         }
         nvars = k;
      }

#ifdef SCIP_DEBUG_SINGLEROWLP
      SCIPdebugMsg(scip, "After exploiting slacks: obj = %g, nvars = %d\n", (*obj), nvars);
#endif

      /* If there are no variables left we can trivially put together a solution or determine infeasibility */
      if( nvars == 0 )
      {
         if( !SCIPisInfinity(scip, mincost) )
         {
            (*obj) += mincost * b;
            return SCIP_OKAY;
         }
         else
         {
            (*solvable) = FALSE;
            return SCIP_OKAY;
         }
      }
      /* Solve the remaining part of the problem */
      else
      {
         assert(nvars > 0);
#ifdef SCIP_DEBUG_SINGLEROWLP
         for( i = 0; i < nvars; i++ )
            SCIPdebugMsg(scip, "c[%d] = %g, a[%d] = %g\n", i, c[i], i, a[i]);
#endif

         SCIPselectWeightedReal(c, a, b, nvars, &k);

#ifdef SCIP_DEBUG_SINGLEROWLP
         SCIPdebugMsg(scip, "k-mean = %g at index %d\n", c[k], k, b);
         for( i = 0; i < nvars; i++ )
            SCIPdebugMsg(scip, "c[%d] = %g, a[%d] = %g\n", i, c[i], i, a[i]);
#endif

         /* Finalize objective value of solution. First we use all elements cheaper than the k-median */
         for( i = 0; i < k; i++ )
         {
            (*obj) += c[i] * a[i];
            b -= a[i];
         }

#ifdef SCIP_DEBUG_SINGLEROWLP
         SCIPdebugMsg(scip, "LP is solved: b = %g\n", b);
#endif

         /* If the constraint is not yet satisfied, we have to fix that */
         if( SCIPisPositive(scip, b) )
         {
            /* There exists an element to satisfy the constraint */
            if( k < nvars )
            {
               (*obj) += c[k] * b;
               return SCIP_OKAY;
            }
            /* There is an upslack variable to satisfy the constraint */
            else if( !SCIPisInfinity(scip, mincost) )
            {
#ifdef SCIP_DEBUG_SINGLEROWLP
               SCIPdebugMsg(scip, "We use %g units of upslack to satisfy the constraint\n", b);
#endif
               (*obj) += mincost * b;
               return SCIP_OKAY;
            }
            /* We cannot satisfy the constraint so the problem is infeasible */
            else
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
         }
         /* The constraint is already satisfied, i.e. b <= 0 */
         else
         {
            return SCIP_OKAY;
         }
      }
   }
}

/** Transform rows into single row LPs, solve them and and tighten bounds
 *
 *  During transformation, create coefficient arrays where variables with a zero coefficient in both rows are ignored
 *  and bring the LP in the form min c^T x, s.t. a^T x >= b, lbs <= x <= ubs.
 *  These LPs are then solved and bounds tightened as described in LP-bound (see above).
 */
static
SCIP_RETCODE transformAndSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object, rows specified by row1idx/row2idx must be sorted */
   int                   row1idx,            /**< index of first row */
   int                   row2idx,            /**< index of second row */
   SCIP_Bool             swaprow1,           /**< should row1 <= rhs be used in addition to lhs <= row1 */
   SCIP_Bool             swaprow2,           /**< should row2 <= rhs be used in addition to lhs <= row2 */
   SCIP_Real*            aoriginal,          /**< buffer array for original constraint coefficients */
   SCIP_Real*            acopy,              /**< buffer array for coefficients adjusted to single-row LP to be solved */
   SCIP_Real*            coriginal,          /**< buffer array for original objective coefficients */
   SCIP_Real*            ccopy,              /**< buffer array for coefficients adjusted to single-row LP to be solved */
   SCIP_Bool*            cangetbnd,          /**< buffer array for flags of which variables a bound can be generated */
   SCIP_Real*            lbs,                /**< buffer array for lower bounds for single-row LP */
   SCIP_Real*            ubs,                /**< buffer array for upper bounds for single-row LP */
   SCIP_Real*            newlbsoriginal,     /**< buffer array for new lower bounds not adjusted to individual single-row LPs */
   SCIP_Real*            newlbscopy,         /**< buffer array for adjusted lower bounds */
   SCIP_Real*            newubsoriginal,     /**< buffer array for new upper bounds not adjusted to individual single-row LPs */
   SCIP_Real*            newubscopy,         /**< buffer array for adjusted upper bounds */
   SCIP_Bool*            success,            /**< return (success || "found better bounds") */
   SCIP_Bool*            infeasible          /**< we return (infeasible || "detected infeasibility") */
   )
{
   int i;
   int j;
   int idx1;
   int idx2;
   int row1len;
   int row2len;
   int* row1idxptr;
   int* row2idxptr;
   SCIP_Real* row1valptr;
   SCIP_Real* row2valptr;
   int nvars;
   SCIP_Real minact;
   SCIP_Real maxact;
   int maxinfs;
   int mininfs;

   SCIP_Bool minsolvable;
   SCIP_Real minobj = SCIP_INVALID;
   SCIP_Bool maxsolvable;
   SCIP_Real maxobj;
   SCIP_Bool minswapsolvable;
   SCIP_Real minswapobj = 0.0;
   SCIP_Bool maxswapsolvable;
   SCIP_Real maxswapobj;

   SCIP_Real newbnd;

   assert(!swaprow1 || (swaprow1 && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, row1idx))));
   assert(!swaprow2 || (swaprow2 && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, row2idx))));

   row1len = SCIPmatrixGetRowNNonzs(matrix, row1idx);
   row2len = SCIPmatrixGetRowNNonzs(matrix, row2idx);
   row1idxptr = SCIPmatrixGetRowIdxPtr(matrix, row1idx);
   row2idxptr = SCIPmatrixGetRowIdxPtr(matrix, row2idx);
   row1valptr = SCIPmatrixGetRowValPtr(matrix, row1idx);
   row2valptr = SCIPmatrixGetRowValPtr(matrix, row2idx);

   /*  Preprocess rows:
    *  1. Calculate minimal and maximal activity of variables not appearing in both rows,
    *     as this represents the right-hand sides of the single-row LPs to be solved.
    *  2. Transform rows into format required by solveSingleRowLP where
    *     first row represents the objective vector c and second row represents the constraint vector a.
    *  3. Determine for which variables new bounds can be calculated.
    */
   i = 0;
   j = 0;
   nvars = 0;
   mininfs = 0;
   maxinfs = 0;
   minact = 0;
   maxact = 0;
   while( i < row1len && j < row2len )
   {
      idx1 = row1idxptr[i];
      idx2 = row2idxptr[j];

      if( idx1 == idx2 )
      {
         coriginal[nvars] = row1valptr[i];
         aoriginal[nvars] = row2valptr[j];
         newlbsoriginal[nvars] = lbs[idx1];
         newubsoriginal[nvars] = ubs[idx1];
         cangetbnd[idx1] = FALSE;
         nvars++;
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "%g <= (%s) <= %g  has coefs %g and %g, %d LP vars\n",
                      lbs[idx1], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx1)),
                      ubs[idx1], row1valptr[i], row2valptr[j], nvars);
#endif
         i++;
         j++;
      }
      else if( idx1 < idx2 )
      {
         if( SCIPisPositive(scip, row1valptr[i]) )
         {
            if( SCIPisInfinity(scip, ubs[idx1]) )
               maxinfs++;
            else
               maxact -= row1valptr[i] * ubs[idx1];

            if( SCIPisInfinity(scip, -lbs[idx1]) )
               mininfs++;
            else
               minact -= row1valptr[i] * lbs[idx1];
         }
         else
         {
            if( SCIPisInfinity(scip, -lbs[idx1]) )
               maxinfs++;
            else
               maxact -= row1valptr[i] * lbs[idx1];

            if( SCIPisInfinity(scip, ubs[idx1]) )
               mininfs++;
            else
               minact -= row1valptr[i] * ubs[idx1];

            cangetbnd[idx1] = TRUE;
         }
         if( maxinfs > 1 && mininfs > 1 )
         {
            (*success) = FALSE;
            return SCIP_OKAY;
         }
         i++;
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "%g <= (%s) <= %g has coefs %g and 0.0, minact = %g, maxact = %g\n",
                      lbs[idx1], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx1)),
                      ubs[idx1], row1valptr[i], minact, maxact);
#endif
      }
      else
      {
         coriginal[nvars] = 0.0;
         aoriginal[nvars] = row2valptr[j];
         newlbsoriginal[nvars] = lbs[idx2];
         newubsoriginal[nvars] = ubs[idx2];
         cangetbnd[idx2] = FALSE;
         nvars++;
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "%g <= (%s) <= %g has coefs 0.0 and %g, %d LP vars\n",
                      lbs[idx2], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx2)),
                      ubs[idx2], row2valptr[j], nvars);
#endif
         j++;
      }
   }
   while( i < row1len )
   {
      idx1 = row1idxptr[i];
      if( SCIPisPositive(scip, row1valptr[i]) )
      {
         if( SCIPisInfinity(scip, ubs[idx1]) )
            maxinfs++;
         else
            maxact -= row1valptr[i] * ubs[idx1];

         if( SCIPisInfinity(scip, -lbs[idx1]) )
            mininfs++;
         else
            minact -= row1valptr[i] * lbs[idx1];
      }
      else
      {
         if( SCIPisInfinity(scip, -lbs[idx1]) )
            maxinfs++;
         else
            maxact -= row1valptr[i] * lbs[idx1];

         if( SCIPisInfinity(scip, ubs[idx1]) )
            mininfs++;
         else
            minact -= row1valptr[i] * ubs[idx1];
      }
      cangetbnd[idx1] = TRUE;
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "%g <= (%s) <= %g  has coefs %g and 0.0, minact = %g, maxact = %g\n",
                   lbs[idx1], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx1)),
                   ubs[idx1], row1valptr[i], minact, maxact);
#endif
      i++;
   }
   while( j < row2len )
   {
      idx2 = row2idxptr[j];
      coriginal[nvars] = 0.0;
      aoriginal[nvars] = row2valptr[j];
      newlbsoriginal[nvars] = lbs[idx2];
      newubsoriginal[nvars] = ubs[idx2];
      nvars++;
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "%g <= (%s) <= %g has coefs 0.0 and %g, %d LP vars\n",
                   lbs[idx2], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx2)),
                   ubs[idx2], row2valptr[j], nvars);
#endif
      j++;
   }

#ifdef SCIP_DEBUG_2RB
   SCIPdebugMsg(scip, "right hand sides: %g and %g\n",
                SCIPmatrixGetRowLhs(matrix, row1idx), SCIPmatrixGetRowLhs(matrix, row2idx));
#endif

   /* solve single-row LPs */
   maxsolvable = FALSE;
   minsolvable = FALSE;
   maxswapsolvable = FALSE;
   minswapsolvable = FALSE;
   /* maximize overlap in first row with lhs <= row2 as constraint */
   if( maxinfs <= 1 )
   {
      for( i = 0; i < nvars; i++ )
      {
         acopy[i] = aoriginal[i];
         ccopy[i] = -coriginal[i];
         newlbscopy[i] = newlbsoriginal[i];
         newubscopy[i] = newubsoriginal[i];
      }
      SCIP_CALL( solveSingleRowLP(scip, acopy, SCIPmatrixGetRowLhs(matrix, row2idx),
                                  ccopy, newlbscopy, newubscopy, nvars, &maxobj, &maxsolvable) );
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "max-LP solved: obj = %g\n", maxobj);
#endif
   }

   /* minimize overlap in first row with lhs <= row2 as constraint */
   if( mininfs == 0 || (mininfs == 1 && swaprow1) )
   {
      /* copy coefficients */
      for( i = 0; i < nvars; i++ )
      {
         acopy[i] = aoriginal[i];
         ccopy[i] = coriginal[i];
         newlbscopy[i] = newlbsoriginal[i];
         newubscopy[i] = newubsoriginal[i];
      }
      SCIP_CALL( solveSingleRowLP(scip, acopy, SCIPmatrixGetRowLhs(matrix, row2idx),
                                  ccopy, newlbscopy, newubscopy, nvars, &minobj, &minsolvable) );
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "min-LP solved: obj = %g\n", minobj);
#endif
   }

   if( swaprow2 )
   {
     /* maximize overlap in first row with row2 <= rhs as constraint */
      if( maxinfs <= 1 )
      {
         /* copy coefficients */
         for( i = 0; i < nvars; i++ )
         {
            acopy[i] = -aoriginal[i];
            ccopy[i] = -coriginal[i];
            newlbscopy[i] = newlbsoriginal[i];
            newubscopy[i] = newubsoriginal[i];
         }
         SCIP_CALL( solveSingleRowLP(scip, acopy, -SCIPmatrixGetRowRhs(matrix, row2idx),
                                     ccopy, newlbscopy, newubscopy, nvars, &maxswapobj, &maxswapsolvable) );
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "maxswap-LP solved: obj = %g\n", maxswapobj);
#endif
      }

      /* minimize overlap in first row with row2 <= rhs as constraint */
      if( mininfs == 0 || (mininfs == 1 && swaprow1) )
      {
         /* copy coefficients */
         for( i = 0; i < nvars; i++ )
         {
            acopy[i] = -aoriginal[i];
            ccopy[i] = coriginal[i];
            newlbscopy[i] = newlbsoriginal[i];
            newubscopy[i] = newubsoriginal[i];
         }
         SCIP_CALL( solveSingleRowLP(scip, acopy, -SCIPmatrixGetRowRhs(matrix, row2idx),
                                     ccopy, newlbscopy, newubscopy, nvars, &minswapobj, &minswapsolvable) );
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "minswap-LP solved: obj = %g\n", minswapobj);
#endif
      }
   }

   /* perform bound tightening, infeasibility checks and redundancy checks */
   if( maxinfs <= 1 && (maxsolvable || maxswapsolvable) )
   {
      SCIP_Real activity;

      if( maxsolvable && maxswapsolvable )
         activity = MAX(maxobj, maxswapobj) + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/
      else if( maxsolvable )
         activity = maxobj + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/
      else
         activity = maxswapobj + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/

      /* infeasibility check */
      if( maxinfs == 0 && SCIPisPositive(scip, activity) )
      {
         (*infeasible) = TRUE;
         (*success) = TRUE;
         return SCIP_OKAY;
      }
      /* strengthen bounds of all variables outside overlap */
      else if( maxinfs == 0 )
      {
         for( i = 0; i < row1len; i++ )
         {
            idx1 = row1idxptr[i];
            if( cangetbnd[idx1] )
            {
               if( SCIPisPositive(scip, row1valptr[i]) )
               {
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPceil(scip, (activity + row1valptr[i] * ubs[idx1]) / row1valptr[i]);
                  else
                     newbnd = (activity + row1valptr[i] * ubs[idx1]) / row1valptr[i];

                  if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n",
                                  lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                     lbs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
               else
               {
                  assert(SCIPisNegative(scip, row1valptr[i]));
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPfloor(scip, (activity + row1valptr[i] * lbs[idx1]) / row1valptr[i]);
                  else
                     newbnd = (activity + row1valptr[i] * lbs[idx1]) / row1valptr[i];

                  if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n",
                                  lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                     ubs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
            }
         }
      }
      /* strengthen bound of the single variable contributing the infinity */
      else
      {
         assert(maxinfs == 1);
         for( i = 0; i < row1len; i++ )
         {
            idx1 = row1idxptr[i];
            if( cangetbnd[idx1] )
            {
               if( SCIPisPositive(scip, row1valptr[i]) && SCIPisInfinity(scip, ubs[idx1]) )
               {
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPceil(scip, activity / row1valptr[i]);
                  else
                     newbnd = activity / row1valptr[i];

                  if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n",
                                  lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                     lbs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
               else if( SCIPisInfinity(scip, -lbs[idx1]) )
               {
                  assert(SCIPisNegative(scip, row1valptr[i]));
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPfloor(scip, activity / row1valptr[i]);
                  else
                     newbnd = activity / row1valptr[i];

                  if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n",
                                  lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                     ubs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
            }
         }
      }
   }

   /* in this case the objective is swapped. therefore the minimum and the maximum of the support switch roles */
   if( swaprow1 )
   {
      /* perform bound tightening, infeasibility checks and redundancy checks */
      if( mininfs <= 1 && (minsolvable || minswapsolvable) )
      {
         SCIP_Real activity;

         assert(minobj != SCIP_INVALID); /*lint !e777*/
         if( minsolvable && minswapsolvable )
            activity = MAX(minobj, minswapobj) - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;
         else if( minsolvable )
            activity = minobj - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;
         else
            activity = minswapobj - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;

         /* infeasibility check */
         if( mininfs == 0 && SCIPisPositive(scip, activity) )
         {
            (*infeasible) = TRUE;
            (*success) = TRUE;
            return SCIP_OKAY;
         }
         /* strengthen bounds of all variables outside overlap */
         else if( mininfs == 0 )
         {
            for( i = 0; i < row1len; i++ )
            {
               idx1 = row1idxptr[i];
               if( cangetbnd[idx1] )
               {
                  if( SCIPisNegative(scip, row1valptr[i]) ) /* since we look at the swapped case, this represents a positive coefficient */
                  {
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPceil(scip, (activity - row1valptr[i] * ubs[idx1]) / (-row1valptr[i]));
                     else
                        newbnd = (activity - row1valptr[i] * ubs[idx1]) / (-row1valptr[i]);

                     if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n",
                                     lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                        lbs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
                  else
                  {
                     /* since we look at the swapped case, this represents a negative coefficient */
                     assert(SCIPisPositive(scip, row1valptr[i]));
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPfloor(scip, (activity - row1valptr[i] * lbs[idx1]) / (-row1valptr[i]));
                     else
                        newbnd = (activity - row1valptr[i] * lbs[idx1]) / (-row1valptr[i]);

                     if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n",
                                     lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                        ubs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
               }
            }
         }
         /* strengthen bound of the single variable contributing the infinity */
         else
         {
            assert(mininfs == 1);
            for( i = 0; i < row1len; i++ )
            {
               idx1 = row1idxptr[i];
               if( cangetbnd[idx1] )
               {
                  /* since we look at the swapped case, this represents a positive coefficient */
                  if( SCIPisNegative(scip, row1valptr[i]) && SCIPisInfinity(scip, ubs[idx1]) )
                  {
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPceil(scip, activity / (-row1valptr[i]));
                     else
                        newbnd = activity / (-row1valptr[i]);

                     if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n", lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                        lbs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
                  else if( SCIPisInfinity(scip, -lbs[idx1]) )
                  {
                     /* since we look at the swapped case, this represents a negative coefficient */
                     assert(SCIPisPositive(scip, row1valptr[i]));
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPfloor(scip, activity / (-row1valptr[i]));
                     else
                        newbnd = activity / (-row1valptr[i]);

                     if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n",
                                     lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                        ubs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** create required buffer arrays and apply LP-based bound tightening in both directions */
static
SCIP_RETCODE applyLPboundTightening(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row1,               /**< index of first row */
   int                   row2,               /**< index of seond row */
   SCIP_Bool             swaprow1,           /**< should row1 <= rhs be used in addition to lhs <= row1 */
   SCIP_Bool             swaprow2,           /**< should row2 <= rhs be used in addition to lhs <= row2 */
   SCIP_Real*            lbs,                /**< lower variable bounds */
   SCIP_Real*            ubs,                /**< upper variable bounds */
   SCIP_Bool*            success             /**< return (success || "found better bounds") */
   )
{
   SCIP_Real* aoriginal;
   SCIP_Real* acopy;
   SCIP_Real* coriginal;
   SCIP_Real* ccopy;
   SCIP_Real* newlbsoriginal;
   SCIP_Real* newlbscopy;
   SCIP_Real* newubsoriginal;
   SCIP_Real* newubscopy;
   SCIP_Bool* cangetbnd;
   SCIP_Bool infeasible;

#ifdef SCIP_DEBUG_2RB
   SCIPdebugMsg(scip, "combining rows %d (%s) and %d (%s)\n",
                row1, SCIPmatrixGetRowName(matrix, row1), row2, SCIPmatrixGetRowName(matrix, row2));
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &aoriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &acopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ccopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlbsoriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlbscopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newubsoriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newubscopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cangetbnd, SCIPmatrixGetNColumns(matrix)) );

   /* Sort matrix rows */
   SCIPsortIntReal(SCIPmatrixGetRowIdxPtr(matrix, row1), SCIPmatrixGetRowValPtr(matrix, row1),
                   SCIPmatrixGetRowNNonzs(matrix, row1));
   SCIPsortIntReal(SCIPmatrixGetRowIdxPtr(matrix, row2), SCIPmatrixGetRowValPtr(matrix, row2),
                   SCIPmatrixGetRowNNonzs(matrix, row2));

   /* Use row2 to strengthen row1 */
   infeasible = FALSE;
   SCIP_CALL( transformAndSolve(scip, matrix, row1, row2, swaprow1, swaprow2, aoriginal, acopy,
                                coriginal, ccopy, cangetbnd, lbs, ubs, newlbsoriginal, newlbscopy,
                                newubsoriginal, newubscopy, success, &infeasible) );

   /* Switch roles and use row1 to strengthen row2 */
   SCIP_CALL( transformAndSolve(scip, matrix, row2, row1, swaprow2, swaprow1, aoriginal, acopy,
                                coriginal, ccopy, cangetbnd, lbs, ubs, newlbsoriginal, newlbscopy,
                                newubsoriginal, newubscopy, success, &infeasible) );

   SCIPfreeBufferArray(scip, &cangetbnd);
   SCIPfreeBufferArray(scip, &newubscopy);
   SCIPfreeBufferArray(scip, &newubsoriginal);
   SCIPfreeBufferArray(scip, &newlbscopy);
   SCIPfreeBufferArray(scip, &newlbsoriginal);
   SCIPfreeBufferArray(scip, &ccopy);
   SCIPfreeBufferArray(scip, &coriginal);
   SCIPfreeBufferArray(scip, &acopy);
   SCIPfreeBufferArray(scip, &aoriginal);

   return SCIP_OKAY;
}

/* Find hashes contained in both hashlists, and apply LP-bound
 * on their corresponding rows. Both hashlists must be sorted.
 */
static
SCIP_RETCODE processHashlists(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int*                  hashlist1,          /**< first list of hashes */
   int*                  hashlist2,          /**< second list of hashes */
   int                   lenhashlist1,       /**< length of first hashlist */
   int                   lenhashlist2,       /**< length of second hashlist */
   int*                  rowidxlist1,        /**< list of row indices corresponding to hashes in hashlist1 */
   int*                  rowidxlist2,        /**< list of row indices corresponding to hashes in hashlist2 */
   SCIP_Real*            newlbs,             /**< lower variable bounds, new bounds will be written here */
   SCIP_Real*            newubs              /**< upper variable bounds, new bound will be written here */
   )
{
   int i;
   int j;
   int block1start;
   int block1end;
   int block2start;
   int block2end;
   SCIP_Longint maxcombines;
   SCIP_Bool finished;
   SCIP_Bool success;
   SCIP_Bool swaprow1;
   SCIP_Bool swaprow2;
   int ncombines;
   int combinefails;
   int retrievefails;
   ROWPAIR rowpair;
   SCIP_HASHSET* pairhashset;

   SCIP_CALL( SCIPhashsetCreate(&pairhashset, SCIPblkmem(scip), 1) );

   finished = FALSE;
   block1start = 0;
   block1end = 0;
   block2start = 0;
   block2end = 0;
   maxcombines = presoldata->maxpairfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)SCIPmatrixGetNRows(matrix)) * presoldata->maxpairfac);

   ncombines = 0;
   combinefails = 0;
   retrievefails = 0;
   findNextBlock(hashlist1, lenhashlist1, &block1start, &block1end);
   findNextBlock(hashlist2, lenhashlist2, &block2start, &block2end);
   while( !finished )
   {
      if( hashlist1[block1start] == hashlist2[block2start] )
      {
         for( i = block1start; i < block1end; i++ )
         {
            for( j = block2start; j < block2end; j++ )
            {
               if( rowidxlist1[i] != rowidxlist2[j] )
               {
                  rowpair.row1idx = MIN(rowidxlist1[i], rowidxlist2[j]);
                  rowpair.row2idx = MAX(rowidxlist1[i], rowidxlist2[j]);
                  if( !SCIPhashsetExists(pairhashset, encodeRowPair(&rowpair)) )
                  {
                     assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row1idx)));
                     assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row2idx)));

                     success = FALSE;

                     /* apply lp-based bound tightening */
                     swaprow1 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx));
                     swaprow2 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx));

                     SCIP_CALL( applyLPboundTightening(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                           swaprow1, swaprow2, newlbs, newubs, &success) );

                     if( success )
                        combinefails = 0;
                     else
                        combinefails++;

                     SCIP_CALL( SCIPhashsetInsert(pairhashset, SCIPblkmem(scip), encodeRowPair(&rowpair)) );
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
               /* check if SCIP ran into a time limit already */
               if( j % 10 == 0 && SCIPisStopped(scip) )
                  finished = TRUE;
               if( finished )
                  break;
            }
            /* check if SCIP ran into a time limit already */
            if( SCIPisStopped(scip) )
               finished = TRUE;
            if( finished )
               break;
         }

         if( block1end < lenhashlist1 && block2end < lenhashlist2 )
         {
            findNextBlock(hashlist1, lenhashlist1, &block1start, &block1end);
            findNextBlock(hashlist2, lenhashlist2, &block2start, &block2end);
         }
         else
            finished = TRUE;
      }
      else if( hashlist1[block1start] < hashlist2[block2start] && block1end < lenhashlist1 )
         findNextBlock(hashlist1, lenhashlist1, &block1start, &block1end);
      else if( hashlist1[block1start] > hashlist2[block2start] && block2end < lenhashlist2 )
         findNextBlock(hashlist2, lenhashlist2, &block2start, &block2end);
      else
         finished = TRUE;
   }

   SCIPhashsetFree(&pairhashset, SCIPblkmem(scip));

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyTworowbnd)
{
   SCIP_PRESOLDATA* presoldata;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver if copying is enabled */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   if( presoldata->enablecopy )
   {
      SCIP_CALL( SCIPincludePresolTworowbnd(scip) );
   }

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeTworowbnd)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitTworowbnd)
{
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   presoldata->nchgbnds = 0;
   presoldata->nuselessruns = 0;

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecTworowbnd)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_Bool infeasible;
   SCIP_PRESOLDATA* presoldata;
   int oldnchgbds;
   int oldnfixedvars;
   int nrows;
   int ncols;
   SCIP_Real* oldlbs;
   SCIP_Real* oldubs;
   SCIP_Real* newlbs;
   SCIP_Real* newubs;
   int* rowidxptr;
   SCIP_Real* rowvalptr;
   SCIP_VAR* var;

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

   int* rowidxlistpp;
   int* rowidxlistmm;
   int* rowidxlistpm;
   int* rowidxlistmp;

   SCIP_Bool finiterhs;

   int i;
   int j;
   int k;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   infeasible = FALSE;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   if( presoldata->nuselessruns >= 5 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

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

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   if( nrows <= 1 )
   {
      SCIPmatrixFree(scip, &matrix);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistpp, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistmm, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistpm, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistmp, nrows) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistpp, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistmm, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistpm, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistmp, nrows) );

   pospp = 0;
   posmm = 0;
   pospm = 0;
   posmp = 0;
   listsizepp = nrows;
   listsizemm = nrows;
   listsizepm = nrows;
   listsizemp = nrows;
   maxhashes = presoldata->maxhashfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)nrows) * presoldata->maxhashfac);

   /* skim through the problem and create hashlists for combination candidates */
   for( i = 0; i < nrows; i++)
   {
      if( ((SCIP_Longint)pospp) + posmm + pospm + posmp > maxhashes )
         break;

      rowvalptr = SCIPmatrixGetRowValPtr(matrix, i);
      rowidxptr = SCIPmatrixGetRowIdxPtr(matrix, i);
      finiterhs = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, i));
      maxlen = MIN(presoldata->maxconsiderednonzeros, SCIPmatrixGetRowNNonzs(matrix, i)); /*lint !e666*/
      for( j = 0; j < maxlen; j++)
      {
         for( k = j+1; k < maxlen; k++)
         {
            if( SCIPisPositive(scip, rowvalptr[j]) )
            {
               if(SCIPisPositive(scip, rowvalptr[k]) )
               {
                  SCIP_CALL( addEntry(scip, &pospp, &listsizepp, &hashlistpp, &rowidxlistpp,
                     hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
                  if( finiterhs )
                     SCIP_CALL( addEntry(scip, &posmm, &listsizemm, &hashlistmm, &rowidxlistmm,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
               }
               else
               {
                  SCIP_CALL( addEntry(scip, &pospm, &listsizepm, &hashlistpm, &rowidxlistpm,
                     hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
                  if( finiterhs )
                     SCIP_CALL( addEntry(scip, &posmp, &listsizemp, &hashlistmp, &rowidxlistmp,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
               }
            }
            else
            {
               if(SCIPisPositive(scip, rowvalptr[k]) )
               {
                  SCIP_CALL( addEntry(scip, &posmp, &listsizemp, &hashlistmp, &rowidxlistmp,
                     hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
                  if( finiterhs )
                     SCIP_CALL( addEntry(scip, &pospm, &listsizepm, &hashlistpm, &rowidxlistpm,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
               }
               else
               {
                  SCIP_CALL( addEntry(scip, &posmm, &listsizemm, &hashlistmm, &rowidxlistmm,
                     hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
                  if( finiterhs )
                     SCIP_CALL( addEntry(scip, &pospp, &listsizepp, &hashlistpp, &rowidxlistpp,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i) );
               }
            }
         }
      }
   }

#ifdef SCIP_DEBUG_HASHING
   SCIPdebugMsg(scip, "pp\n");
   for( i = 0; i < pospp; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpp[i], rowidxlistpp[i]);
   SCIPdebugMsg(scip, "mm\n");
   for( i = 0; i < posmm; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmm[i], rowidxlistmm[i]);
   SCIPdebugMsg(scip, "pm\n");
   for( i = 0; i < pospm; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpm[i], rowidxlistpm[i]);
   SCIPdebugMsg(scip, "mp\n");
   for( i = 0; i < posmp; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmp[i], rowidxlistmp[i]);
#endif
   SCIPdebugMsg(scip, "hashlist sizes: pp %d, mm %d, pm %d, mp %d \n", pospp, posmm, pospm, posmp);

   SCIPsortIntInt(hashlistpp, rowidxlistpp, pospp);
   SCIPsortIntInt(hashlistmm, rowidxlistmm, posmm);
   SCIPsortIntInt(hashlistpm, rowidxlistpm, pospm);
   SCIPsortIntInt(hashlistmp, rowidxlistmp, posmp);

#ifdef SCIP_DEBUG_HASHING
   SCIPdebugMsg(scip, "sorted pp\n");
   for( i = 0; i < pospp; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpp[i], rowidxlistpp[i]);
   SCIPdebugMsg(scip, "sorted mm\n");
   for( i = 0; i < posmm; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmm[i], rowidxlistmm[i]);
   SCIPdebugMsg(scip, "sorted pm\n");
   for( i = 0; i < pospm; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpm[i], rowidxlistpm[i]);
   SCIPdebugMsg(scip, "sorted mp\n");
   for( i = 0; i < posmp; i++)
      SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmp[i], rowidxlistmp[i]);
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &oldlbs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldubs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlbs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newubs, ncols) );

   for( i = 0; i < SCIPmatrixGetNColumns(matrix); i++ )
   {
      var = SCIPmatrixGetVar(matrix, i);
      oldlbs[i] = SCIPvarGetLbLocal(var);
      oldubs[i] = SCIPvarGetUbLocal(var);
      newlbs[i] = oldlbs[i];
      newubs[i] = oldubs[i];
   }

   /* Process pp and mm hashlists */
   if( pospp > 0 && posmm > 0 )
   {
     SCIPdebugMsg(scip, "processing pp and mm\n");
     SCIP_CALL( processHashlists(scip, presoldata, matrix, hashlistpp, hashlistmm, pospp, posmm, rowidxlistpp,
                                 rowidxlistmm, newlbs, newubs) );
   }

   /* Process pm and mp hashlists */
   if( pospm > 0 && posmp > 0 )
   {
     SCIPdebugMsg(scip, "processing pm and mp\n");
     SCIP_CALL( processHashlists(scip, presoldata, matrix, hashlistpm, hashlistmp, pospm, posmp, rowidxlistpm,
                                 rowidxlistmp, newlbs, newubs) );
   }

   /* Apply reductions */
   oldnchgbds = *nchgbds;
   oldnfixedvars = *nfixedvars;
   for( i = 0; i < SCIPmatrixGetNColumns(matrix); i++ )
   {
      SCIP_Bool bndwastightened;
      SCIP_Bool fixed;

      var = SCIPmatrixGetVar(matrix, i);

      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT
         || (SCIPisEQ(scip, newlbs[i], SCIPceil(scip, newlbs[i])) && (SCIPisEQ(scip, newubs[i], SCIPfloor(scip, newubs[i])))));

      if( SCIPisEQ(scip, newlbs[i], newubs[i]) )
      {
         SCIP_CALL( SCIPfixVar(scip, var, newlbs[i], &infeasible, &fixed) );

         if( infeasible )
         {
            SCIPdebugMessage(" -> infeasible fixing of variable %s\n", SCIPvarGetName(var));
            break;
         }

         if( fixed )
         {
            SCIPdebugMessage("variable %s fixed to %g\n", SCIPvarGetName(var), newlbs[i]);
            (*nfixedvars)++;
         }
      }

      if( SCIPisLT(scip, oldlbs[i], newlbs[i]) )
      {
         SCIP_CALL( SCIPtightenVarLb(scip, var, newlbs[i], FALSE, &infeasible, &bndwastightened) );

         if( infeasible )
         {
            SCIPdebugMessage(" -> infeasible lower bound tightening: %s >= %g\n", SCIPvarGetName(var), newlbs[i]);
            break;
         }

         if( bndwastightened )
         {
            SCIPdebugMessage("lower bound of %s changed from %g to %g\n", SCIPvarGetName(var), oldlbs[i], newlbs[i]);
            (*nchgbds)++;
         }
      }

      if( SCIPisGT(scip, oldubs[i], newubs[i]) )
      {
         SCIP_CALL( SCIPtightenVarUb(scip, var, newubs[i], FALSE, &infeasible, &bndwastightened) );

         if( infeasible )
         {
            SCIPdebugMessage(" -> infeasible upper bound tightening: %s <= %g\n", SCIPvarGetName(var), newubs[i]);
            break;
         }

         if( bndwastightened )
         {
            SCIPdebugMessage("upper bound of %s changed from %g to %g\n", SCIPvarGetName(var), oldubs[i], newubs[i]);
            (*nchgbds)++;
         }
      }
   }

   /* set result */
   if( *nchgbds > oldnchgbds || *nfixedvars > oldnfixedvars )
   {
      *result = SCIP_SUCCESS;
      presoldata->nuselessruns = 0;
   }
   else if( infeasible )
   {
      *result = SCIP_CUTOFF;
   }
   else
   {
      presoldata->nuselessruns++;
   }

   SCIPfreeBufferArray(scip, &newubs);
   SCIPfreeBufferArray(scip, &newlbs);
   SCIPfreeBufferArray(scip, &oldubs);
   SCIPfreeBufferArray(scip, &oldlbs);
   SCIPfreeBlockMemoryArray(scip, &rowidxlistmp, listsizemp);
   SCIPfreeBlockMemoryArray(scip, &rowidxlistpm, listsizepm);
   SCIPfreeBlockMemoryArray(scip, &rowidxlistmm, listsizemm);
   SCIPfreeBlockMemoryArray(scip, &rowidxlistpp, listsizepp);
   SCIPfreeBlockMemoryArray(scip, &hashlistmp, listsizemp);
   SCIPfreeBlockMemoryArray(scip, &hashlistpm, listsizepm);
   SCIPfreeBlockMemoryArray(scip, &hashlistmm, listsizemm);
   SCIPfreeBlockMemoryArray(scip, &hashlistpp, listsizepp);

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the tworowbndb presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolTworowbnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create tworowbnd presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   presol = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecTworowbnd,
         presoldata) );

   assert(presol != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyTworowbnd) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeTworowbnd) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitTworowbnd) );

   /* add tworowbnd presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/tworowbnd/enablecopy",
         "should tworowbnd presolver be copied to sub-SCIPs?",
         &presoldata->enablecopy, TRUE, DEFAULT_ENABLECOPY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxconsiderednonzeros",
         "maximal number of considered non-zeros within one row (-1: no limit)",
         &presoldata->maxconsiderednonzeros, FALSE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxretrievefails",
         "maximal number of consecutive useless hashtable retrieves",
         &presoldata->maxretrievefails, FALSE, DEFAULT_MAXRETRIEVEFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxcombinefails",
         "maximal number of consecutive useless row combines",
         &presoldata->maxcombinefails, FALSE, DEFAULT_MAXCOMBINEFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxhashfac",
         "Maximum number of hashlist entries as multiple of number of rows in the problem (-1: no limit)",
         &presoldata->maxhashfac, FALSE, DEFAULT_MAXHASHFAC, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxpairfac",
         "Maximum number of processed row pairs as multiple of the number of rows in the problem (-1: no limit)",
         &presoldata->maxpairfac, FALSE, DEFAULT_MAXPAIRFAC, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
