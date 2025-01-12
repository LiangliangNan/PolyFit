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

/**@file   presol_dualsparsify.c
 * @brief  cancel nonzeros of the constraint matrix based on the columns
 * @author Dieter Weninger
 * @author Leona Gottwald
 * @author Ambros Gleixner
 * @author Weikun Chen
 *
 * This presolver attempts to cancel non-zero entries of the constraint
 * matrix by adding scaled columns to other columns.
 *
 * In more detail, for two columns A_{j.} and A_{k.}, suppose for a given value s, we have
 *
 *                  | A_{j.} | - | A_{j.} - s*A_{k.} | > eta,
 *
 * where eta is an nonnegative integer. Then we introduce a new variable y := s*x_j + x_k
 * and aggregate the variable x_k = y - s*x_j. After aggregation, the column of the variable
 * x_j is A_{j.} + s*A_{j.} which is sparser than A_{j.}. In the case that x_k is nonimplied
 * free variable, we need to add a new constraint l_k <= y - weight*x_j <= u_k into the problem
 * to keep the bounds constraints of variable x_k.
 *
 * Further information can be found in
 * Chen et al. "Two-row and two-column mixed-integer presolve using hasing-based pairing methods".
 *
 * @todo add infrastructure to SCIP for handling aggregated binary variables
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/debug.h"
#include "scip/presol_dualsparsify.h"
#include "scip/pub_cons.h"
#include "scip/pub_matrix.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_presol.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"
#include <string.h>

#define PRESOL_NAME            "dualsparsify"
#define PRESOL_DESC            "eliminate non-zero coefficients"

#define PRESOL_PRIORITY           -240000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_ENABLECOPY           TRUE    /**< should dualsparsify presolver be copied to sub-SCIPs? */
#define DEFAULT_PRESERVEINTCOEFS    FALSE    /**< should we forbid cancellations that destroy integer coefficients? */
#define DEFAULT_PRESERVEGOODLOCKS   FALSE    /**< should we preserve good locked properties of variables (at most one lock in one direction)? */
#define DEFAULT_MAX_CONT_FILLIN         1    /**< default value for the maximal fillin for continuous variables */
#define DEFAULT_MAX_BIN_FILLIN          1    /**< default value for the maximal fillin for binary variables */
#define DEFAULT_MAX_INT_FILLIN          1    /**< default value for the maximal fillin for integer variables (including binary) */
#define DEFAULT_MAXCONSIDEREDNONZEROS  70    /**< maximal number of considered nonzeros within one column (-1: no limit) */
#define DEFAULT_MINELIMINATEDNONZEROS 100    /**< minimal eleminated nonzeros within one column if we need to add a constraint to the problem */
#define DEFAULT_MAXRETRIEVEFAC      100.0    /**< limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints */
#define DEFAULT_WAITINGFAC            2.0    /**< number of calls to wait until next execution as a multiple of the number of useless calls */

#define MAXSCALE                   1000.0    /**< maximal allowed scale for cancelling nonzeros */


/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int                   ncancels;           /**< total number of canceled nonzeros (net value, i.e., removed minus added nonzeros) */
   int                   nfillin;            /**< total number of added nonzeros */
   int                   nfailures;          /**< number of calls to presolver without success */
   int                   nwaitingcalls;      /**< number of presolver calls until next real execution */
   int                   naggregated;        /**< number of aggregated variables */
   int                   maxcontfillin;      /**< maximal fillin for continuous variables */
   int                   maxintfillin;       /**< maximal fillin for integer variables*/
   int                   maxbinfillin;       /**< maximal fillin for binary variables */
   int                   maxconsiderednonzeros;/**< maximal number of considered nonzeros within one column (-1: no limit) */
   int                   mineliminatednonzeros;/**< minimal eliminated nonzeros within one column if we need to add a constraint to the problem */
   SCIP_Real             maxretrievefac;     /**< limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints */
   SCIP_Real             waitingfac;         /**< number of calls to wait until next execution as a multiple of the number of useless calls */
   SCIP_Bool             enablecopy;         /**< should dualsparsify presolver be copied to sub-SCIPs? */
   SCIP_Bool             preserveintcoefs;   /**< should we forbid cancellations that destroy integer coefficients? */
   SCIP_Bool             preservegoodlocks;  /**< should we preserve good locked properties of variables (at most one lock in one direction)? */
};

/** structure representing a pair of constraints in a column; used for lookup in a hashtable */
struct ColConsPair
{
   int colindex;                             /**< index of the column */
   int consindex1;                           /**< index of the first constraint */
   int consindex2;                           /**< index of the second constraint */
   SCIP_Real conscoef1;                      /**< coefficient of the first constraint */
   SCIP_Real conscoef2;                      /**< coefficient of the second constriant */
};

typedef struct ColConsPair COLCONSPAIR;

/*
 * Local methods
 */

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(consPairsEqual)
{  /*lint --e{715}*/
   SCIP* scip;
   COLCONSPAIR* conspair1;
   COLCONSPAIR* conspair2;
   SCIP_Real ratio1;
   SCIP_Real ratio2;

   scip = (SCIP*) userptr;
   conspair1 = (COLCONSPAIR*) key1;
   conspair2 = (COLCONSPAIR*) key2;

   if( conspair1->consindex1 != conspair2->consindex1 )
      return FALSE;

   if( conspair1->consindex2 != conspair2->consindex2 )
      return FALSE;

   ratio1 = conspair1->conscoef2 / conspair1->conscoef1;
   ratio2 = conspair2->conscoef2 / conspair2->conscoef1;

   return SCIPisEQ(scip, ratio1, ratio2);
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(consPairHashval)
{  /*lint --e{715}*/
   COLCONSPAIR* conspair;

   conspair = (COLCONSPAIR*) key;

   return SCIPhashThree(conspair->consindex1, conspair->consindex2, SCIPrealHashCode(conspair->conscoef2 / conspair->conscoef1));
}

/** calculate maximal activity of one row without one specific column */
static
SCIP_Real getMaxActivitySingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col                 /**< column index */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   SCIP_Real maxactivity;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);

   maxactivity = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      if( val > 0.0 )
      {
         ub = SCIPmatrixGetColUb(matrix, c);
         assert(!SCIPisInfinity(scip, ub));

         maxactivity += val * ub;
      }
      else if( val < 0.0 )
      {
         lb = SCIPmatrixGetColLb(matrix, c);
         assert(!SCIPisInfinity(scip, -lb));

         maxactivity += val * lb;
      }
   }

   return maxactivity;
}

/** calculate minimal activity of one row without one specific column */
static
SCIP_Real getMinActivitySingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col                 /**< column index */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   SCIP_Real minactivity;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);

   minactivity = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      if( val > 0.0 )
      {
         lb = SCIPmatrixGetColLb(matrix, c);
         assert(!SCIPisInfinity(scip, -lb));

         minactivity += val * lb;
      }
      else if( val < 0.0 )
      {
         ub = SCIPmatrixGetColUb(matrix, c);
         assert(!SCIPisInfinity(scip, ub));

         minactivity += val * ub;
      }
   }

   return minactivity;
}

/** get minimal and maximal residual activity without one specified column */
static
void getMinMaxActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this variable in this row */
   SCIP_Real*            minresactivity,     /**< minimum residual activity of this row */
   SCIP_Real*            maxresactivity,     /**< maximum residual activity of this row */
   SCIP_Bool*            isminsettoinfinity, /**< flag indicating if minresactiviy is set to infinity */
   SCIP_Bool*            ismaxsettoinfinity  /**< flag indicating if maxresactiviy is set to infinity */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   int nmaxactneginf;
   int nmaxactposinf;
   int nminactneginf;
   int nminactposinf;
   SCIP_Real maxactivity;
   SCIP_Real minactivity;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   lb = SCIPmatrixGetColLb(matrix, col);
   ub = SCIPmatrixGetColUb(matrix, col);

   *isminsettoinfinity = FALSE;
   *ismaxsettoinfinity = FALSE;

   nmaxactneginf = SCIPmatrixGetRowNMaxActNegInf(matrix, row);
   nmaxactposinf = SCIPmatrixGetRowNMaxActPosInf(matrix, row);
   nminactneginf = SCIPmatrixGetRowNMinActNegInf(matrix, row);
   nminactposinf = SCIPmatrixGetRowNMinActPosInf(matrix, row);

   maxactivity = SCIPmatrixGetRowMaxActivity(matrix, row);
   minactivity = SCIPmatrixGetRowMinActivity(matrix, row);

   if( val >= 0.0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(nmaxactposinf >= 1);
         if( nmaxactposinf == 1 && nmaxactneginf == 0 )
            *maxresactivity = getMaxActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nmaxactneginf + nmaxactposinf) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = maxactivity - val * ub;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nminactneginf >= 1);
         if( nminactneginf == 1 && nminactposinf == 0 )
            *minresactivity = getMinActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nminactneginf + nminactposinf) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = minactivity - val * lb;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nmaxactneginf >= 1);
         if( nmaxactneginf == 1 && nmaxactposinf == 0 )
            *maxresactivity = getMaxActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nmaxactneginf + nmaxactposinf) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = maxactivity - val * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(nminactposinf >= 1);
         if( nminactposinf == 1 && nminactneginf == 0 )
            *minresactivity = getMinActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nminactneginf + nminactposinf) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = minactivity - val * ub;
      }
   }
}

/** calculate the upper and lower bound of one variable from one row */
static
void getVarBoundsOfRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index of variable */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this column in this row */
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

   getMinMaxActivityResiduals(scip, matrix, col, row, val,
      &minresactivity, &maxresactivity,
      &isminsettoinfinity, &ismaxsettoinfinity);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);

   if( val > 0.0 )
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) )
      {
         *rowub = (rhs - minresactivity) / val;
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
   SCIP_Bool*            ubimplied,          /**< flag indicating an implied upper bound */
   SCIP_Bool*            lbimplied           /**< flag indicating an implied lower bound */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real impliedub;
   SCIP_Real impliedlb;
   SCIP_Real ub;
   SCIP_Real lb;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(ubimplied != NULL);
   assert(lbimplied != NULL);

   *ubimplied = FALSE;
   impliedub = SCIPinfinity(scip);

   *lbimplied = FALSE;
   impliedlb = -SCIPinfinity(scip);

   ub =  SCIPmatrixGetColUb(matrix, col);
   lb =  SCIPmatrixGetColLb(matrix, col);

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);
   for( ; (colpnt < colend); colpnt++, valpnt++ )
   {
      SCIP_Real rowub;
      SCIP_Bool ubfound;
      SCIP_Real rowlb;
      SCIP_Bool lbfound;

      getVarBoundsOfRow(scip, matrix, col, *colpnt, *valpnt, &rowub, &ubfound, &rowlb, &lbfound);

      if( ubfound && (rowub < impliedub) )
         impliedub = rowub;

      if( lbfound && (rowlb > impliedlb) )
         impliedlb = rowlb;
   }

   /* we consider +/-inf bounds as implied bounds */
   if( SCIPisInfinity(scip, ub) ||
      (!SCIPisInfinity(scip, ub) && SCIPisLE(scip, impliedub, ub)) )
      *ubimplied = TRUE;

   if( SCIPisInfinity(scip, -lb) ||
      (!SCIPisInfinity(scip, -lb) && SCIPisGE(scip, impliedlb, lb)) )
      *lbimplied = TRUE;
}

/** y = weight1 * var[colidx1] + var[colidx2] */
static
SCIP_RETCODE aggregation(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_MATRIX*          matrix,             /**< matrix datastructure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_VAR**            vars,               /**< the current variables */
   int                   colidx1,            /**< one of the indexes of column to try nonzero cancellation for */
   int                   colidx2,            /**< one of the indexes of column to try nonzero cancellation for */
   SCIP_Bool             isimpliedfree,      /**< is the aggregated variable implied free? */
   SCIP_Real             weight1             /**< weight variable one in the aggregated expression */
   )
{
   SCIP_VAR* tmpvars[2];
   SCIP_Real coefs[2];
   char newvarname[SCIP_MAXSTRLEN];
   char newconsname[SCIP_MAXSTRLEN];
   SCIP_CONS* newcons;
   SCIP_VAR* aggregatedvar;
   SCIP_VAR* newvar;
   SCIP_VARTYPE newvartype;
   SCIP_Real constant;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool infeasible;
   SCIP_Bool aggregated;
#ifndef NDEBUG
   if( isimpliedfree )
   {
      SCIP_Bool lbimplied;
      SCIP_Bool ubimplied;

      getImpliedBounds(scip, matrix, colidx2, &ubimplied, &lbimplied);
      assert(lbimplied && ubimplied);
   }
#endif

   assert( !SCIPisZero(scip, weight1) );

   presoldata->naggregated += 1;
   aggregatedvar = vars[colidx2];

   /* if the variable is implied free, we make sure that the columns bounds are removed,
    * so that subsequent checks for implied bounds do not interfere with the exploitation
    * of this variables implied bounds
    */
   if( isimpliedfree )
   {
      SCIPdebugMsg(scip, "remove column bounds of column %d\n", colidx2);
      SCIPmatrixRemoveColumnBounds(scip, matrix, colidx2);
   }

   assert(!SCIPdoNotMultaggrVar(scip, aggregatedvar));

   (void) SCIPsnprintf(newvarname, SCIP_MAXSTRLEN, "dualsparsifyvar_%d", presoldata->naggregated);

   constant = 0.0;

   if( weight1 > 0.0 )
   {
      if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(vars[colidx1])) ||
         SCIPisInfinity(scip, -SCIPvarGetLbGlobal(vars[colidx2])) )
         newlb = -SCIPinfinity(scip);
      else
         newlb = weight1 * SCIPvarGetLbGlobal(vars[colidx1]) + SCIPvarGetLbGlobal(vars[colidx2]);

      if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx1])) ||
         SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx2])) )
         newub = SCIPinfinity(scip);
      else
         newub = weight1 * SCIPvarGetUbGlobal(vars[colidx1]) + SCIPvarGetUbGlobal(vars[colidx2]);
   }
   else
   {
      if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx1])) ||
         SCIPisInfinity(scip, -SCIPvarGetLbGlobal(vars[colidx2])) )
         newlb = -SCIPinfinity(scip);
      else
         newlb = weight1 * SCIPvarGetUbGlobal(vars[colidx1]) + SCIPvarGetLbGlobal(vars[colidx2]);

      if( SCIPisInfinity(scip, SCIPvarGetLbGlobal(vars[colidx1])) ||
         SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx2])) )
         newub = SCIPinfinity(scip);
      else
         newub = weight1 * SCIPvarGetLbGlobal(vars[colidx1]) + SCIPvarGetUbGlobal(vars[colidx2]);
   }

   if( SCIPvarIsIntegral(aggregatedvar) )
      newvartype = (SCIPvarGetType(aggregatedvar) == SCIP_VARTYPE_IMPLINT) ?
         SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_INTEGER;
   else
      newvartype = SCIP_VARTYPE_CONTINUOUS;

   lhs = SCIPvarGetLbGlobal(vars[colidx2]);
   rhs = SCIPvarGetUbGlobal(vars[colidx2]);

   SCIP_CALL( SCIPcreateVar(scip, &newvar, newvarname, newlb, newub, 0.0, newvartype,
         SCIPvarIsInitial(aggregatedvar), SCIPvarIsRemovable(aggregatedvar), NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, newvar) );

   /* set the debug solution value for the new variable */
#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      SCIP_Real val1;
      SCIP_Real val2;

      SCIP_CALL( SCIPdebugGetSolVal(scip, vars[colidx1], &val1) );
      SCIP_CALL( SCIPdebugGetSolVal(scip, vars[colidx2], &val2) );
      SCIP_CALL( SCIPdebugAddSolVal(scip, newvar, weight1 * val1 + val2) );

      SCIPdebugMsg(scip, "set debug solution value of %s to %g\n", SCIPvarGetName(newvar), weight1 * val1 + val2);
   }
#endif

   tmpvars[0] = vars[colidx1];
   tmpvars[1] = newvar;
   coefs[0] = -weight1;
   coefs[1] = 1.0;

   SCIP_CALL( SCIPmultiaggregateVar(scip, aggregatedvar, 2, tmpvars, coefs, constant, &infeasible, &aggregated) );

   assert(!infeasible);
   assert(aggregated);

   vars[colidx2] = newvar;

   /* create a linear constraint that ensures that var[colidx2].lb <= y - weight1 * var[colidx1] <= var[colidx2].ub;
    * note that it might happen that vars[colidx2] is not implied free even though it has infinite bounds because
    * getImpliedBounds() considers infinite bounds to be implied
    */
   if( !isimpliedfree && (!SCIPisInfinity(scip, rhs) || !SCIPisInfinity(scip, -lhs)) )
   {
      SCIPdebugMsg(scip, "create a linear constraint to ensure %g <= %g %s + %g %s <= %g\n", lhs, coefs[0], SCIPvarGetName(tmpvars[0]),
         coefs[1], SCIPvarGetName(tmpvars[1]), rhs);
      (void) SCIPsnprintf(newconsname, SCIP_MAXSTRLEN, "dualsparsifycons_%d", presoldata->naggregated);

      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, newconsname, 2, tmpvars, coefs,
            lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );

      SCIPdebugPrintCons(scip, newcons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   }

   SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

   return SCIP_OKAY;
}

/** try nonzero cancellation for given column */
static
SCIP_RETCODE cancelCol(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_MATRIX*          matrix,             /**< the constraint matrix */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_HASHTABLE*       pairtable,          /**< the hashtable containing COLCONSPAIR's of equations */
   SCIP_Bool*            ishashingcols,      /**< array to indicates whether it is impliedfree or not */
   SCIP_VAR**            vars,               /**< array to store the current variables */
   SCIP_Bool*            isblockedvar,       /**< array to indicates whether it is blocked or not */
   int                   colidx,             /**< index of column to try nonzero cancellation for */
   int                   maxcontfillin,      /**< maximal fill-in allowed for continuous variables */
   int                   maxintfillin,       /**< maximal fill-in allowed for integral variables */
   int                   maxbinfillin,       /**< maximal fill-in allowed for binary variables */
   int                   maxconsiderednonzeros, /**< maximal number of nonzeros to consider for cancellation */
   SCIP_Bool             preserveintcoefs,   /**< only perform nonzero cancellation if integrality of coefficients is preserved? */
   SCIP_Longint*         nuseless,           /**< pointer to update number of useless hashtable retrieves */
   int*                  nchgcoefs,          /**< pointer to update number of changed coefficients */
   int*                  ncanceled,          /**< pointer to update number of canceled nonzeros */
   int*                  nfillin,            /**< pointer to update the produced fill-in */
   SCIP_Bool             isaddedcons         /**< whether a linear constraint required to added to keep the validity */
   )
{
   SCIP_VAR* cancelvar;
   SCIP_Real* cancelcolvals;
   SCIP_Real* colvalptr;
   SCIP_Real* tmpvals;
   SCIP_Real* scores;
   int* cancelcolinds;
   int* colidxptr;
   int* tmpinds;
   SCIP_Real bestcancelrate;
   SCIP_Real bestscale;
   SCIP_Real ncols;
   SCIP_Bool colishashing;
   SCIP_Bool swapped = FALSE;
   int cancelcollen;
   int bestnfillin;
   int nretrieves;
   int maxfillin;
   int bestcand;
   int nchgcoef;

   ncols = SCIPmatrixGetNColumns(matrix);
   colishashing = ishashingcols[colidx];
   cancelcollen = SCIPmatrixGetColNNonzs(matrix, colidx);
   colidxptr = SCIPmatrixGetColIdxPtr(matrix, colidx);
   colvalptr = SCIPmatrixGetColValPtr(matrix, colidx);
   cancelvar = vars[colidx];

   if( SCIPvarIsIntegral(cancelvar) )
   {
      if( SCIPvarIsBinary(cancelvar) )
         maxfillin = maxbinfillin;
      else
         maxfillin = maxintfillin;
   }
   else
      maxfillin = maxcontfillin;

   SCIP_CALL( SCIPduplicateBufferArray(scip, &cancelcolinds, colidxptr, cancelcollen) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &cancelcolvals, colvalptr, cancelcollen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpinds, cancelcollen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvals, cancelcollen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, cancelcollen) );

   nchgcoef = 0;
   nretrieves = 0;
   while( TRUE ) /*lint !e716 */
   {
      COLCONSPAIR colconspair;
      int maxlen;
      int i;
      int j;

      bestcand = -1;
      bestnfillin = 0;
      bestscale = 1.0;
      bestcancelrate = 0.0;

      /* sort the rows non-decreasingly by number of nonzeros
       * if the number of nonzeros, we use the colindex as tie-breaker
       */
      for( i = 0; i < cancelcollen; ++i )
      {
         tmpinds[i] = i;
         scores[i] = -SCIPmatrixGetRowNNonzs(matrix, cancelcolinds[i]) - 1.0 * cancelcolinds[i] / (ncols);
      }
      SCIPsortRealInt(scores, tmpinds, cancelcollen);

      maxlen = cancelcollen;
      if( maxconsiderednonzeros >= 0 )
         maxlen = MIN(cancelcollen, maxconsiderednonzeros);

      for( i = 0; i < maxlen; ++i )
      {
         for( j = i + 1; j < maxlen; ++j )
         {
            COLCONSPAIR* hashingcolconspair;
            SCIP_VAR* hashingcolvar;
            SCIP_Real* hashingcolvals;
            int* hashingcolinds;
            SCIP_Real hashingcollb;
            SCIP_Real hashingcolub;
            SCIP_Real cancelrate;
            SCIP_Real rowlhs;
            SCIP_Real rowrhs;
            SCIP_Real scale;
            SCIP_Bool hashingcolisbin;
            SCIP_Bool abortpair;
            int hashingcollen;
            int ntotfillin;
            int ncancel;
            int a,b;
            int i1,i2;

            i1 = tmpinds[i];
            i2 = tmpinds[j];

            assert(cancelcolinds[i] < cancelcolinds[j]);

            if( cancelcolinds[i1] < cancelcolinds[i2] )
            {
               colconspair.consindex1 = cancelcolinds[i1];
               colconspair.consindex2 = cancelcolinds[i2];
               colconspair.conscoef1 = cancelcolvals[i1];
               colconspair.conscoef2 = cancelcolvals[i2];
            }
            else
            {
               colconspair.consindex1 = cancelcolinds[i2];
               colconspair.consindex2 = cancelcolinds[i1];
               colconspair.conscoef1 = cancelcolvals[i2];
               colconspair.conscoef2 = cancelcolvals[i1];
            }

            hashingcolconspair = (COLCONSPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &colconspair);
            nretrieves++;

            if( hashingcolconspair == NULL ||
               hashingcolconspair->colindex == colidx || isblockedvar[hashingcolconspair->colindex] )
               continue;

            /* if the column we want to cancel is a hashing column (which we stored for canceling other columns),
             * we will only use the hashing columns for canceling with less nonzeros and if the number of nonzeros
             * is equal we use the colindex as tie-breaker to avoid cyclic nonzero cancellation
             */
            hashingcollen = SCIPmatrixGetColNNonzs(matrix, hashingcolconspair->colindex);
            if( colishashing && (cancelcollen < hashingcollen ||
                  (cancelcollen == hashingcollen && colidx < hashingcolconspair->colindex)) )
               continue;

            hashingcolvals = SCIPmatrixGetColValPtr(matrix, hashingcolconspair->colindex);
            hashingcolinds = SCIPmatrixGetColIdxPtr(matrix, hashingcolconspair->colindex);
            hashingcolvar = vars[hashingcolconspair->colindex];
            hashingcollb = SCIPvarGetLbGlobal(hashingcolvar);
            hashingcolub = SCIPvarGetUbGlobal(hashingcolvar);
            hashingcolisbin = (SCIPvarGetType(hashingcolvar) == SCIP_VARTYPE_BINARY) ||
               (SCIPvarIsIntegral(hashingcolvar) && SCIPisZero(scip, hashingcollb) &&
                  SCIPisEQ(scip, hashingcolub, 1.0));
            scale = -colconspair.conscoef1 / hashingcolconspair->conscoef1;

            if( SCIPisZero(scip, scale) )
               continue;

            if( REALABS(scale) > MAXSCALE )
               continue;

            /* @todo do more reduction if knspsack constraint handler supports downgrading constraint,
             * i.e., converting into a linear constraint
             */
            if( hashingcolisbin )
               continue;
            else if( SCIPvarIsIntegral(hashingcolvar) )
            {
               if( SCIPvarIsIntegral(cancelvar) )
               {
                  /* skip if the hashing variable is an integer variable and
                   * the canceled variable is an implied integer variable
                   */
                  if( (SCIPvarGetType(hashingcolvar) != SCIP_VARTYPE_IMPLINT) &&
                     (SCIPvarGetType(cancelvar) == SCIP_VARTYPE_IMPLINT) )
                     continue;

                  /* skip if the scale is non-integral */
                  if( !SCIPisIntegral(scip, scale) )
                     continue;

                  /* round scale to be exactly integral */
                  scale = floor(scale + 0.5);
               }
               /* skip if the canceled variable is a continuous variable */
               else
                  continue;
            }

            a = 0;
            b = 0;
            ncancel = 0;
            ntotfillin = 0;
            abortpair = FALSE;

            while( a < cancelcollen && b < hashingcollen )
            {
               if( cancelcolinds[a] == hashingcolinds[b] )
               {
                  SCIP_Real newcoef;

                  newcoef = cancelcolvals[a] + scale * hashingcolvals[b];

                  /* check if coefficient is canceled */
                  if( SCIPisZero(scip, newcoef) )
                  {
                     ++ncancel;
                  }
                  /* otherwise, check if integral coefficients are preserved if the column is integral */
                  else if( (preserveintcoefs && SCIPvarIsIntegral(cancelvar) &&
                        SCIPisIntegral(scip, cancelcolvals[a]) && !SCIPisIntegral(scip, newcoef)) )
                  {
                     abortpair = TRUE;
                     break;
                  }
                  /* finally, check if locks could be modified in a bad way due to flipped signs */
                  else if( COPYSIGN(1.0, newcoef) != COPYSIGN(1.0, cancelcolvals[a]) ) /*lint !e777*/
                  {
                     /* do not flip signs for non-canceled coefficients if this adds a lock to a variable that
                      * had at most one lock in that direction before, except if the other direction gets unlocked
                      */
                     rowrhs = SCIPmatrixGetRowRhs(matrix, cancelcolinds[a]);
                     rowlhs = SCIPmatrixGetRowLhs(matrix, cancelcolinds[a]);
                     if( (cancelcolvals[a] > 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                        (cancelcolvals[a] < 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                     {
                        /* if we get into this case the variable had a positive coefficient in a <= constraint or
                         * a negative coefficient in a >= constraint, e.g. an uplock. If this was the only uplock
                         * we do not abort their cancelling, otherwise we abort if we had a single or no downlock
                         * and add one
                         */
                        if( presoldata->preservegoodlocks && (SCIPmatrixGetColNUplocks(matrix, colidx) > 1 &&
                              SCIPmatrixGetColNDownlocks(matrix, colidx) <= 1) )
                        {
                           abortpair = TRUE;
                           break;
                        }
                     }

                     if( (cancelcolvals[a] < 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                        (cancelcolvals[a] > 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                     {
                        /* symmetric case where the variable had a downlock */
                        if( presoldata->preservegoodlocks && (SCIPmatrixGetColNDownlocks(matrix, colidx) > 1 &&
                              SCIPmatrixGetColNUplocks(matrix, colidx) <= 1) )
                        {
                           abortpair = TRUE;
                           break;
                        }
                     }
                  }

                  ++a;
                  ++b;
               }
               else if( cancelcolinds[a] < hashingcolinds[b] )
               {
                  ++a;
               }
               else
               {
                  SCIP_Real newcoef;

                  newcoef = scale * hashingcolvals[b];
                  rowrhs = SCIPmatrixGetRowRhs(matrix, hashingcolinds[b]);
                  rowlhs = SCIPmatrixGetRowLhs(matrix, hashingcolinds[b]);

                  if( (newcoef > 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                     (newcoef < 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                  {
                     if( presoldata->preservegoodlocks && SCIPmatrixGetColNUplocks(matrix, colidx) <= 1 )
                     {
                        abortpair = TRUE;
                        ++b;
                        break;
                     }
                  }

                  if( (newcoef < 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                     (newcoef > 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                  {
                     if( presoldata->preservegoodlocks && SCIPmatrixGetColNDownlocks(matrix, colidx) <= 1 )
                     {
                        abortpair = TRUE;
                        ++b;
                        break;
                     }
                  }

                  ++b;

                  if( ++ntotfillin > maxfillin )
                  {
                     abortpair = TRUE;
                     break;
                  }
               }
            }

            if( abortpair )
               continue;

            while( b < hashingcollen )
            {
               ++b;

               if( ++ntotfillin > maxfillin )
                  break;
            }
CHECKFILLINAGAIN:
            if( ntotfillin > maxfillin || ntotfillin >= ncancel )
               continue;

            cancelrate = (ncancel - ntotfillin) / (SCIP_Real) cancelcollen;

            /* if a linear constraint is needed to keep the validity, we require a large nonzero cancellation */
            if( isaddedcons && (ncancel - ntotfillin < presoldata->mineliminatednonzeros) )
               continue;

            if( cancelrate > bestcancelrate )
            {
               if( ishashingcols[hashingcolconspair->colindex] )
               {
                  SCIP_Bool lbimplied;
                  SCIP_Bool ubimplied;

                  /* recompute whether a variable is still implied free; after some previous multi-aggregations of
                   * some variables, it might be that other variables that are contained in the same linear rows of the
                   * matrix are not implied free anymore (see #2971)
                   */
                  getImpliedBounds(scip, matrix, hashingcolconspair->colindex, &ubimplied, &lbimplied);

                  if( !lbimplied || !ubimplied )
                  {
                     ishashingcols[hashingcolconspair->colindex] = FALSE;
                     ntotfillin += 2;
                     goto CHECKFILLINAGAIN;
                  }
               }

               bestnfillin = ntotfillin;
               bestcand = hashingcolconspair->colindex;
               bestscale = scale;
               bestcancelrate = cancelrate;

               /* stop looking if the current candidate does not create any fill-in or alter coefficients */
               if( cancelrate == 1.0 )
                  break;
            }

            /* we accept the best candidate immediately if it does not create any fill-in or alter coefficients */
            if( bestcand != -1 && bestcancelrate == 1.0 )
               break;
         }
      }

      if( bestcand != -1 )
      {
         SCIP_Real* hashingcolvals;
         int* hashingcolinds;
         int hashingcollen;
         int tmpcollen;
         int a;
         int b;

         SCIPdebugMsg(scip, "cancelcol %d (%s) candidate column %d (%s) (bestcancelrate = %g, bestscale = %g)\n",
            colidx, SCIPvarGetName(cancelvar), bestcand, SCIPvarGetName(vars[bestcand]), bestcancelrate, bestscale);

         hashingcolvals = SCIPmatrixGetColValPtr(matrix, bestcand);
         hashingcolinds = SCIPmatrixGetColIdxPtr(matrix, bestcand);
         hashingcollen = SCIPmatrixGetColNNonzs(matrix, bestcand);

         a = 0;
         b = 0;
         tmpcollen = 0;

         while( a < cancelcollen && b < hashingcollen )
         {
            if( cancelcolinds[a] == hashingcolinds[b] )
            {
               SCIP_Real val = cancelcolvals[a] + bestscale * hashingcolvals[b];

               if( !SCIPisZero(scip, val) )
               {
                  tmpinds[tmpcollen] = cancelcolinds[a];
                  tmpvals[tmpcollen] = val;
                  ++tmpcollen;
               }
               ++nchgcoef;

               ++a;
               ++b;
            }
            else if( cancelcolinds[a] < hashingcolinds[b] )
            {
               tmpinds[tmpcollen] = cancelcolinds[a];
               tmpvals[tmpcollen] = cancelcolvals[a];
               ++tmpcollen;
               ++a;
            }
            else
            {
               tmpinds[tmpcollen] = hashingcolinds[b];
               tmpvals[tmpcollen] = hashingcolvals[b] * bestscale;
               ++nchgcoef;
               ++tmpcollen;
               ++b;
            }
         }

         while( a < cancelcollen )
         {
            tmpinds[tmpcollen] = cancelcolinds[a];
            tmpvals[tmpcollen] = cancelcolvals[a];
            ++tmpcollen;
            ++a;
         }

         while( b < hashingcollen )
         {
            tmpinds[tmpcollen] = hashingcolinds[b];
            tmpvals[tmpcollen] = hashingcolvals[b] * bestscale;
            ++nchgcoef;
            ++tmpcollen;
            ++b;
         }

         /* update fill-in counter */
         *nfillin += bestnfillin;

         /* swap the temporary arrays so that the cancelcolinds and cancelcolvals arrays, contain the new
          * changed column, and the tmpinds and tmpvals arrays can be overwritten in the next iteration
          */
         SCIPswapPointers((void**) &tmpinds, (void**) &cancelcolinds);
         SCIPswapPointers((void**) &tmpvals, (void**) &cancelcolvals);
         swapped = ! swapped;
         cancelcollen = tmpcollen;
         SCIP_CALL( aggregation(scip, matrix, presoldata, vars, colidx, bestcand, ishashingcols[bestcand], -bestscale) );

         /* the newly created variable is now at the position bestcand and is assumed to have the same coefficients.
          * this is not the case if the variable is not implied free since then a new constraint was added and the
          * nonzero fillin would not be counted correctly if we do not block this variable
          */
         if( !ishashingcols[bestcand] )
            isblockedvar[bestcand] = TRUE;
      }
      else
         break;
   }

   if( nchgcoef != 0 )
   {
      /* update counters */
      *nchgcoefs += nchgcoef;
      *ncanceled += SCIPmatrixGetColNNonzs(matrix, colidx) - cancelcollen;

      isblockedvar[colidx] = TRUE;

      /* if successful, decrease the useless hashtable retrieves counter; the rationale here is that we want to keep
       * going if, after many useless calls that almost exceeded the budget, we finally reach a useful section; but we
       * don't allow a negative build-up for the case that the useful section is all at the beginning and we just want
       * to quit quickly afterwards
       */
      *nuseless -= nretrieves;
      *nuseless = MAX(*nuseless, 0);
   }
   else
   {
      /* if not successful, increase useless hashtable retrieves counter */
      *nuseless += nretrieves;
   }

   SCIPfreeBufferArray(scip, &scores);
   if( swapped )
   {
      SCIPfreeBufferArray(scip, &cancelcolvals);
      SCIPfreeBufferArray(scip, &cancelcolinds);
      SCIPfreeBufferArray(scip, &tmpvals);
      SCIPfreeBufferArray(scip, &tmpinds);
   }
   else
   {
      SCIPfreeBufferArray(scip, &tmpvals);
      SCIPfreeBufferArray(scip, &tmpinds);
      SCIPfreeBufferArray(scip, &cancelcolvals);
      SCIPfreeBufferArray(scip, &cancelcolinds);
   }

   return SCIP_OKAY;
}

/** updates failure counter after one execution */
static
void updateFailureStatistic(
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_Bool             success             /**< was this execution successful? */
   )
{
   assert(presoldata != NULL);

   if( success )
   {
      presoldata->nfailures = 0;
      presoldata->nwaitingcalls = 0;
   }
   else
   {
      presoldata->nfailures++;
      presoldata->nwaitingcalls = (int)(presoldata->waitingfac*(SCIP_Real)presoldata->nfailures);
   }
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualsparsify)
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
      SCIP_CALL( SCIPincludePresolDualsparsify(scip) );
   }

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualsparsify)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   int* perm;
   int* colidxsorted;
   int* colsparsity;
   SCIP_Real* scores;
   COLCONSPAIR* conspairs;
   SCIP_HASHTABLE* pairtable;
   SCIP_PRESOLDATA* presoldata;
   SCIP_Bool* ishashingcols;
   SCIP_Bool* isblockedvar;
   SCIP_VAR** vars;
   SCIP_Longint maxuseless;
   SCIP_Longint nuseless;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_Bool infeasible;
   int ncols;
   int c;
   int i;
   int j;
   int conspairssize;
   int nconspairs;
   int numcancel;
   int nfillin;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPdoNotAggr(scip) )
      return SCIP_OKAY;

   /* If restart is performed, some cuts will be tranformed into linear constraints.
    * However, SCIPmatrixCreate() only collects the original constraints (not the constraints transformed from cuts)
    * For this reason, we only perform this method in the first run of branch-and-cut.
    * */
   if( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);

   if( presoldata->nwaitingcalls > 0 )
   {
      presoldata->nwaitingcalls--;
      SCIPdebugMsg(scip, "skipping dualsparsify: nfailures=%d, nwaitingcalls=%d\n", presoldata->nfailures,
         presoldata->nwaitingcalls);
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "starting dualsparsify. . .\n");
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

   ncols = SCIPmatrixGetNColumns(matrix);

   /* sort columns by row indices */
   for( i = 0; i < ncols; i++ )
   {
      int* colpnt = SCIPmatrixGetColIdxPtr(matrix, i);
      SCIP_Real* valpnt = SCIPmatrixGetColValPtr(matrix, i);
      SCIPsortIntReal(colpnt, valpnt, SCIPmatrixGetColNNonzs(matrix, i));
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &scores, SCIPmatrixGetNRows(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, SCIPmatrixGetNRows(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ishashingcols, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isblockedvar, SCIPmatrixGetNColumns(matrix)) );

   /* loop over all columns and create cons pairs */
   conspairssize = 0;
   nconspairs = 0;
   conspairs = NULL;
   SCIP_CALL( SCIPhashtableCreate(&pairtable, SCIPblkmem(scip), 1,
         SCIPhashGetKeyStandard, consPairsEqual, consPairHashval, (void*) scip) );

   /* collect implied free variables and their number of nonzeros */
   for( c = 0; c < ncols; c++ )
   {
      SCIP_Bool lbimplied;
      SCIP_Bool ubimplied;
      int nnonz;

      vars[c] = SCIPmatrixGetVar(matrix, c);

      /* if the locks do not match do not consider the column for sparsification */
      if( SCIPmatrixDownlockConflict(matrix, c) || SCIPmatrixUplockConflict(matrix, c) )
      {
         isblockedvar[c] = TRUE;
         ishashingcols[c] = FALSE;
         continue;
      }

      /* skip if the variable is not allowed to be multi-aggregated */
      if( SCIPdoNotMultaggrVar(scip, vars[c]) )
      {
         isblockedvar[c] = TRUE;
         ishashingcols[c] = FALSE;
         continue;
      }

      nnonz = SCIPmatrixGetColNNonzs(matrix, c);

      getImpliedBounds(scip, matrix, c, &ubimplied, &lbimplied);

      ishashingcols[c] = FALSE;

      if( lbimplied && ubimplied )
         ishashingcols[c] = TRUE;

      isblockedvar[c] = FALSE;

      /* only consider implied free variables
       * skip singleton variables, because either the constraint is redundant
       * or the variables can be canceled by variable substitution
       */
      if( nnonz >= 2 && (lbimplied && ubimplied) )
      {
         SCIP_Real* colvals;
         int* colinds;
         int failshift;
         int npairs;

         colinds = SCIPmatrixGetColIdxPtr(matrix, c);
         colvals = SCIPmatrixGetColValPtr(matrix, c);

         /* sort the rows non-decreasingly by number of nonzeros
          * if the number of nonzeros is equal, we use the colindex as tie-breaker
          */
         for( i = 0; i < nnonz; ++i )
         {
            perm[i] = i;
            scores[i] = -SCIPmatrixGetRowNNonzs(matrix, colinds[i]) - 1.0  *colinds[i] / ncols;
         }
         SCIPsortRealInt(scores, perm, nnonz);

         if( presoldata->maxconsiderednonzeros >= 0 )
            nnonz = MIN(nnonz, presoldata->maxconsiderednonzeros);

         npairs = (nnonz * (nnonz - 1)) / 2;
         if( nconspairs + npairs > conspairssize )
         {
            int newsize = SCIPcalcMemGrowSize(scip, nconspairs + npairs);
            SCIP_CALL( SCIPreallocBufferArray(scip, &conspairs, newsize) );
            conspairssize = newsize;
         }

         /* if we are called after one or more failures, i.e., executions without finding cancellations, then we
          * shift the section of nonzeros considered; in the case that the maxconsiderednonzeros limit is hit, this
          * results in different constraint pairs being tried and avoids trying the same useless cancellations
          * repeatedly
          */
         failshift = presoldata->nfailures*presoldata->maxconsiderednonzeros;

         for( i = 0; i < nnonz; ++i )
         {
            for( j = i + 1; j < nnonz; ++j )
            {
               int i1;
               int i2;

               assert(nconspairs < conspairssize);
               assert(conspairs != NULL);

               i1 = perm[(i + failshift) % nnonz];
               i2 = perm[(j + failshift) % nnonz];
               /* coverity[var_deref_op] */
               conspairs[nconspairs].colindex = c;

               if( colinds[i1] < colinds[i2])
               {
                  conspairs[nconspairs].consindex1 = colinds[i1];
                  conspairs[nconspairs].consindex2 = colinds[i2];
                  conspairs[nconspairs].conscoef1 = colvals[i1];
                  conspairs[nconspairs].conscoef2 = colvals[i2];
               }
               else
               {
                  conspairs[nconspairs].consindex1 = colinds[i2];
                  conspairs[nconspairs].consindex2 = colinds[i1];
                  conspairs[nconspairs].conscoef1 = colvals[i2];
                  conspairs[nconspairs].conscoef2 = colvals[i1];
               }
               ++nconspairs;
            }
         }
      }
   }

   /* insert conspairs into hash table */
   for( c = 0; c < nconspairs; ++c )
   {
      COLCONSPAIR* otherconspair;
      SCIP_Bool insert;

      assert(conspairs != NULL);

      insert = TRUE;

      /* check if this pair is already contained in the hash table;
       * The loop is required due to the non-transitivity of the hash functions
       */
      while( (otherconspair = (COLCONSPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &conspairs[c])) != NULL )
      {
         /* if the previous constraint pair has fewer or the same number of nonzeros in the attached column
          * we keep that pair and skip this one
          */
         if( SCIPmatrixGetColNNonzs(matrix, otherconspair->colindex) <=
            SCIPmatrixGetColNNonzs(matrix, conspairs[c].colindex) )
         {
            insert = FALSE;
            break;
         }

         /* this pairs column has fewer nonzeros, so remove the other pair from the hash table and loop */
         SCIP_CALL( SCIPhashtableRemove(pairtable, (void*) otherconspair) );
      }

      if( insert )
      {
         SCIP_CALL( SCIPhashtableInsert(pairtable, (void*) &conspairs[c]) );
      }
   }

   /* sort cols according to decreasing sparsity */
   SCIP_CALL( SCIPallocBufferArray(scip, &colidxsorted, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colsparsity, ncols) );
   for( c = 0; c < ncols; ++c )
   {
      colidxsorted[c] = c;
      colsparsity[c] = -SCIPmatrixGetColNNonzs(matrix, c);
   }
   SCIPsortIntInt(colsparsity, colidxsorted, ncols);

   /* loop over the columns and cancel nonzeros until maximum number of retrieves is reached */
   maxuseless = (SCIP_Longint)(presoldata->maxretrievefac * (SCIP_Real)ncols);
   nuseless = 0;
   numcancel = 0;
   for( c = 0; c < ncols && nuseless <= maxuseless && !SCIPisStopped(scip); c++ )
   {
      int colidx;

      colidx = colidxsorted[c];

      if( isblockedvar[colidx] )
         continue;

      /* since the function parameters for the max fillin are unsigned we do not need to handle the
       * unlimited (-1) case due to implicit conversion rules */
      SCIP_CALL( cancelCol(scip, matrix, presoldata, pairtable, ishashingcols, vars, isblockedvar, colidx, \
            presoldata->maxcontfillin == -1 ? INT_MAX : presoldata->maxcontfillin, \
            presoldata->maxintfillin == -1 ? INT_MAX : presoldata->maxintfillin, \
            presoldata->maxbinfillin == -1 ? INT_MAX : presoldata->maxbinfillin, \
            presoldata->maxconsiderednonzeros, presoldata->preserveintcoefs, \
            &nuseless, nchgcoefs, &numcancel, &nfillin, FALSE) );
   }

   if( numcancel > 0 )
      *result = SCIP_SUCCESS;
   else /* do reductions on variables that contain larger nonzero entries */
   {
      SCIPhashtableRemoveAll(pairtable);
      nconspairs = 0;

      /* collect large nonzero entries variables and their number of nonzeros */
      for( c = 0; c < ncols; c++ )
      {
         int nnonz;

         nnonz = SCIPmatrixGetColNNonzs(matrix, c);
         vars[c] = SCIPmatrixGetVar(matrix, c);

         /* if the locks do not match do not consider the column for sparsification */
         if( SCIPmatrixDownlockConflict(matrix, c) || SCIPmatrixUplockConflict(matrix, c) )
         {
            isblockedvar[c] = TRUE;
            ishashingcols[c] = FALSE;
            continue;
         }

         isblockedvar[c] = FALSE;

         /* only consider nonimplied free variables, i.e., non-hashing columns in the previous step,
          * with large nonzero entries
          * skip singleton variables, because either the constraint is redundant
          * or the variables can be canceled by variables substitution
          */
         if( nnonz >= presoldata->mineliminatednonzeros && !ishashingcols[c] )
         {
            int* colinds;
            SCIP_Real* colvals;
            int npairs;
            int failshift;

            ishashingcols[c] = TRUE;
            colinds = SCIPmatrixGetColIdxPtr(matrix, c);
            colvals = SCIPmatrixGetColValPtr(matrix, c);

            /* sort the rows non-decreasingly by number of nonzeros
             * if the number of nonzeros, we use the colindex as tie-breaker
             */
            for( i = 0; i < nnonz; ++i )
            {
               perm[i] = i;
               scores[i] = -SCIPmatrixGetRowNNonzs(matrix, colinds[i]) - 1.0 * colinds[i] / ncols;
            }
            SCIPsortRealInt(scores, perm, nnonz);

            if( presoldata->maxconsiderednonzeros >= 0 )
               nnonz = MIN(nnonz, presoldata->maxconsiderednonzeros);

            npairs = (nnonz * (nnonz - 1)) / 2;
            if( nconspairs + npairs > conspairssize )
            {
               int newsize = SCIPcalcMemGrowSize(scip, nconspairs + npairs);
               SCIP_CALL( SCIPreallocBufferArray(scip, &conspairs, newsize) );
               conspairssize = newsize;
            }

            /* if we are called after one or more failures, i.e., executions without finding cancellations, then we
             * shift the section of nonzeros considered; in the case that the maxconsiderednonzeros limit is hit,
             * this results in different constraint pairs being tried and avoids trying the same useless
             * cancellations repeatedly
             */
            failshift = presoldata->nfailures*presoldata->maxconsiderednonzeros;

            for( i = 0; i < nnonz; ++i )
            {
               for( j = i + 1; j < nnonz; ++j )
               {
                  int i1;
                  int i2;

                  assert(nconspairs < conspairssize);
                  assert(conspairs != NULL);

                  i1 = perm[(i + failshift) % nnonz];
                  i2 = perm[(j + failshift) % nnonz];
                  conspairs[nconspairs].colindex = c;

                  if( colinds[i1] < colinds[i2])
                  {
                     conspairs[nconspairs].consindex1 = colinds[i1];
                     conspairs[nconspairs].consindex2 = colinds[i2];
                     conspairs[nconspairs].conscoef1 = colvals[i1];
                     conspairs[nconspairs].conscoef2 = colvals[i2];
                  }
                  else
                  {
                     conspairs[nconspairs].consindex1 = colinds[i2];
                     conspairs[nconspairs].consindex2 = colinds[i1];
                     conspairs[nconspairs].conscoef1 = colvals[i2];
                     conspairs[nconspairs].conscoef2 = colvals[i1];
                  }
                  ++nconspairs;
               }
            }
         }
         else
         {
            ishashingcols[c] = FALSE;
         }
      }

      /* insert conspairs into hash table */
      for( c = 0; c < nconspairs; ++c )
      {
         SCIP_Bool insert;
         COLCONSPAIR* otherconspair;

         assert(conspairs != NULL);

         insert = TRUE;

         /* check if this pair is already contained in the hash table;
          * The loop is required due to the non-transitivity of the hash functions
          */
         while( (otherconspair = (COLCONSPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &conspairs[c])) != NULL )
         {
            /* if the previous constraint pair has fewer or the same number of nonzeros in the attached column
             * we keep that pair and skip this one
             */
            if( SCIPmatrixGetColNNonzs(matrix, otherconspair->colindex) <=
               SCIPmatrixGetColNNonzs(matrix, conspairs[c].colindex) )
            {
               insert = FALSE;
               break;
            }

            /* this pairs column has fewer nonzeros, so remove the other pair from the hash table and loop */
            SCIP_CALL( SCIPhashtableRemove(pairtable, (void*) otherconspair) );
         }

         if( insert )
         {
            SCIP_CALL( SCIPhashtableInsert(pairtable, (void*) &conspairs[c]) );
         }
      }

      /* sort rows according to decreasingly sparsity */
      assert(colidxsorted != NULL);
      assert(colsparsity != NULL);
      for( c = 0; c < ncols; ++c )
      {
         colidxsorted[c] = c;
         colsparsity[c] = -SCIPmatrixGetColNNonzs(matrix, c);
      }
      SCIPsortIntInt(colsparsity, colidxsorted, ncols);

      /* loop over the columns and cancel nonzeros until maximum number of retrieves is reached */
      maxuseless = (SCIP_Longint)(presoldata->maxretrievefac * (SCIP_Real)ncols);
      nuseless = 0;
      for( c = 0; c < ncols && nuseless <= maxuseless; c++ )
      {
         int colidx;
         int nnonz;

         colidx = colidxsorted[c];
         nnonz = SCIPmatrixGetColNNonzs(matrix, colidx);

         if( isblockedvar[colidx] || nnonz < presoldata->mineliminatednonzeros )
            continue;

         /* since the function parameters for the max fillin are unsigned we do not need to handle the
          * unlimited (-1) case due to implicit conversion rules */
         SCIP_CALL( cancelCol(scip, matrix, presoldata, pairtable, ishashingcols, vars, isblockedvar, colidx, \
               presoldata->maxcontfillin == -1 ? INT_MAX : presoldata->maxcontfillin, \
               presoldata->maxintfillin == -1 ? INT_MAX : presoldata->maxintfillin, \
               presoldata->maxbinfillin == -1 ? INT_MAX : presoldata->maxbinfillin, \
               presoldata->maxconsiderednonzeros, presoldata->preserveintcoefs, \
               &nuseless, nchgcoefs, &numcancel, &nfillin, TRUE) );
      }

      if( numcancel > 0 )
      {
         *result = SCIP_SUCCESS;
      }
   }

   updateFailureStatistic(presoldata, numcancel > 0);

   SCIPfreeBufferArray(scip, &colsparsity);
   SCIPfreeBufferArray(scip, &colidxsorted);

   SCIPhashtableFree(&pairtable);
   SCIPfreeBufferArrayNull(scip, &conspairs);

   SCIPfreeBufferArray(scip, &isblockedvar);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &ishashingcols);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &scores);

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeDualsparsify)
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
SCIP_DECL_PRESOLINIT(presolInitDualsparsify)
{
   SCIP_PRESOLDATA* presoldata;

   /* set the counters in the init (and not in the initpre) callback such that they persist across restarts */
   presoldata = SCIPpresolGetData(presol);
   presoldata->ncancels = 0;
   presoldata->nfillin = 0;
   presoldata->nfailures = 0;
   presoldata->nwaitingcalls = 0;
   presoldata->naggregated = 0;

   return SCIP_OKAY;
}

/** creates the dualsparsify presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualsparsify(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create dualsparsify presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDualsparsify, presoldata) );

   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDualsparsify) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeDualsparsify) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitDualsparsify) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/dualsparsify/enablecopy",
         "should dualsparsify presolver be copied to sub-SCIPs?",
         &presoldata->enablecopy, TRUE, DEFAULT_ENABLECOPY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/dualsparsify/preserveintcoefs",
         "should we forbid cancellations that destroy integer coefficients?",
         &presoldata->preserveintcoefs, TRUE, DEFAULT_PRESERVEINTCOEFS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/dualsparsify/preservegoodlocks",
         "should we preserve good locked properties of variables (at most one lock in one direction)?",
         &presoldata->preservegoodlocks, TRUE, DEFAULT_PRESERVEGOODLOCKS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxcontfillin",
         "maximal fillin for continuous variables (-1: unlimited)",
         &presoldata->maxcontfillin, FALSE, DEFAULT_MAX_CONT_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxbinfillin",
         "maximal fillin for binary variables (-1: unlimited)",
         &presoldata->maxbinfillin, FALSE, DEFAULT_MAX_BIN_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxintfillin",
         "maximal fillin for integer variables including binaries (-1: unlimited)",
         &presoldata->maxintfillin, FALSE, DEFAULT_MAX_INT_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxconsiderednonzeros",
         "maximal number of considered nonzeros within one column (-1: no limit)",
         &presoldata->maxconsiderednonzeros, TRUE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/mineliminatednonzeros",
         "minimal eliminated nonzeros within one column if we need to add a constraint to the problem",
         &presoldata->mineliminatednonzeros, FALSE, DEFAULT_MINELIMINATEDNONZEROS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/maxretrievefac",
         "limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints",
         &presoldata->maxretrievefac, TRUE, DEFAULT_MAXRETRIEVEFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/waitingfac",
         "number of calls to wait until next execution as a multiple of the number of useless calls",
         &presoldata->waitingfac, TRUE, DEFAULT_WAITINGFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
