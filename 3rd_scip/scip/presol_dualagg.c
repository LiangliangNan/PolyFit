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

/**@file   presol_dualagg.c
 * @brief  aggregate variables by dual arguments
 * @author Dieter Weninger
 *
 * This presolver looks for variables which could not be handled by
 * duality fixing because of one up-/downlock.
 * If the constraint which delivers the up-/downlock has
 * a specific structure, we can aggregate the corresponding variable.
 *
 * In more detail (for a minimization problem and the case of only one uplock):
 *
 * Given a variable \f$x_i\f$ with \f$c_i \leq 0\f$ and only one up lock (originating from a constraint c),
 * we are looking for a binary variable \f$x_j\f$ such that:
 * 1. if \f$x_j = 0\f$, constraint c can only be fulfilled for \f$x_i = lb_i\f$, and
 * 2. if \f$x_j = 1\f$, constraint c becomes redundant and \f$x_i\f$ can be dual-fixed to its upper bound \f$ub_i\f$
 * (or vice versa). Then we can perform the following aggregation: \f$x_i = lb_i + x_j (ub_i - lb_i)\f$.
 *
 * Similar arguments apply for the case of only one down lock and \f$c_i \geq 0\f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/pub_matrix.h"
#include "presol_dualagg.h"

#define PRESOL_NAME            "dualagg"
#define PRESOL_DESC            "aggregate variables by dual arguments"
#define PRESOL_PRIORITY           -12000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               0     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

/** type of aggregation */
enum AggrType
{
   BIN0UBOUND = -1,         /**< x_j = u_j + (l_j-u_j)x_i with x_i binary and x_j aggregation variable */
   NOAGG      =  0,         /**< do not aggregate */
   BIN0LBOUND =  1          /**< x_j = l_j + (u_j-l_j)x_i with x_i binary and x_j aggregation variable */
};
typedef enum AggrType AGGRTYPE;

/*
 * Local methods
 */

/** find row which leads to the uplock of the given variable */
static
void getUplockRowIdx(
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int                   aggvaridx,          /**< index of variable which should be aggregated */
   int*                  rowidx,             /**< pointer to store row index of uplock */
   SCIP_Real*            coef                /**< pointer to store coefficient of variable */
   )
{
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;

   assert(rowidx != NULL);
   assert(coef != NULL);
   assert(SCIPmatrixGetColNUplocks(matrix, aggvaridx) == 1);

   /* get nonzero entries of the variable in the matrix */
   colpnt = SCIPmatrixGetColIdxPtr(matrix, aggvaridx);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, aggvaridx);
   valpnt = SCIPmatrixGetColValPtr(matrix, aggvaridx);

   /* iterate over all non-zero coefficients of the column */
   *rowidx = -1;
   for(; (colpnt < colend); colpnt++, valpnt++)
   {
      /* currently we support only >= relation */
      if( !SCIPmatrixIsRowRhsInfinity(matrix, *colpnt) )
         break;

      /* coef < 0 for >= relation: this row provides an uplock for the variable */
      if( *valpnt < 0.0 )
      {
         *rowidx = *colpnt;
         *coef = *valpnt;
         break;
      }
   }
#ifndef NDEBUG
   /* in debug mode, we check that the lock number is correct */
   assert(colpnt < colend);
   for(colpnt++, valpnt++; (colpnt < colend); colpnt++, valpnt++)
   {
      assert(*valpnt > 0.0);
   }
#endif

}

/** find row which leads to the downlock of the given variable */
static
void getDownlockRowIdx(
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int                   aggvaridx,          /**< index of variable which should be aggregated */
   int*                  rowidx,             /**< pointer to store row index of downlock */
   SCIP_Real*            coef                /**< pointer to store coefficient of variable */
   )
{
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;

   assert(rowidx != NULL);
   assert(coef != NULL);
   assert(SCIPmatrixGetColNDownlocks(matrix, aggvaridx) == 1);

   /* get nonzero entries of the variable in the matrix */
   colpnt = SCIPmatrixGetColIdxPtr(matrix, aggvaridx);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, aggvaridx);
   valpnt = SCIPmatrixGetColValPtr(matrix, aggvaridx);

   /* iterate over all non-zero coefficients of the column */
   *rowidx = -1;
   for(; (colpnt < colend); colpnt++, valpnt++)
   {
      /* currently we support only >= relation */
      if( !SCIPmatrixIsRowRhsInfinity(matrix, *colpnt) )
         break;

      /* coef > 0 for >= relation: this row provides a downlock for the variable */
      if( *valpnt > 0.0 )
      {
         *rowidx = *colpnt;
         *coef = *valpnt;
         break;
      }
   }
#ifndef NDEBUG
   /* in debug mode, we check that the lock number is correct */
   assert(colpnt < colend);
   for(colpnt++, valpnt++; (colpnt < colend); colpnt++, valpnt++)
   {
      assert(*valpnt < 0.0);
   }
#endif
}

/** find fitting binary variable aggregation for uplock case */
static
void getBinVarIdxInUplockRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int                   aggvaridx,          /**< index of variable which should be aggregated */
   int*                  binvaridx,          /**< pointer to store index of binary variable */
   AGGRTYPE*             aggtype             /**< pointer to store type of aggregation */
   )
{
   int rowidx;
   SCIP_Real coef;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   SCIP_Real minact;
   SCIP_Real maxact;
   SCIP_Real lhs;
   SCIP_Real lb;

   assert(binvaridx != NULL);
   assert(aggtype != NULL);

   *binvaridx = -1;
   *aggtype = NOAGG;

   getUplockRowIdx(matrix, aggvaridx, &rowidx, &coef);

   if( rowidx < 0 )
      return;

   assert(coef < 0);
   minact = SCIPmatrixGetRowMinActivity(matrix, rowidx);
   maxact = SCIPmatrixGetRowMaxActivity(matrix, rowidx);

   if( SCIPisInfinity(scip, -minact) || SCIPisInfinity(scip, maxact) )
      return;

   lhs = SCIPmatrixGetRowLhs(matrix, rowidx);
   lb = SCIPmatrixGetColLb(matrix, aggvaridx);

   /* search for appropriate binary variables */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, rowidx);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, rowidx);
   valpnt = SCIPmatrixGetRowValPtr(matrix, rowidx);
   for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
   {
      SCIP_VAR* var;

      if( *rowpnt == aggvaridx )
         continue;

      var = SCIPmatrixGetVar(matrix, *rowpnt);

      /* avoid cases where the binary variable has lb=ub=1 or lb=ub=0 */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY &&
          SCIPmatrixGetColLb(matrix, *rowpnt) < 0.5 &&
          SCIPmatrixGetColUb(matrix, *rowpnt) > 0.5 )
      {
         SCIP_Real bincoef;
         bincoef = *valpnt;

         if( bincoef < 0 )
         {
            /* binvar = 0 implies that the constraint is redundant */
            if( SCIPisGE(scip, minact-bincoef, lhs) )
            {
               /* binvar = 1 implies that aggvar = lb */
               SCIP_Real bnd;
               bnd = (lhs - maxact + coef*lb - bincoef) / coef;
               if( SCIPisGE(scip, lb, bnd) )
               {
                  *binvaridx = *rowpnt;
                  *aggtype = BIN0UBOUND;
                  break;
               }
            }
         }

         if( bincoef > 0 )
         {
            /* binvar = 1 implies that the constraint is redundant */
            if( SCIPisGE(scip, minact+bincoef, lhs) )
            {
               /* binvar = 0 implies that aggvar = lb */
               SCIP_Real bnd;
               bnd = (lhs - maxact + coef*lb + bincoef) / coef;
               if( SCIPisGE(scip, lb, bnd) )
               {
                  *binvaridx = *rowpnt;
                  *aggtype = BIN0LBOUND;
               }
            }
         }
      }
   }
}

/** find fitting binary variable aggregation for downlock case */
static
void getBinVarIdxInDownlockRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int                   aggvaridx,          /**< index of variable which should be aggregated */
   int*                  binvaridx,          /**< pointer to store index of binary variable */
   AGGRTYPE*             aggtype             /**< pointer to store type of aggregation */
   )
{
   int rowidx;
   SCIP_Real coef;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   SCIP_Real minact;
   SCIP_Real maxact;
   SCIP_Real lhs;
   SCIP_Real ub;

   assert(binvaridx != NULL);
   assert(aggtype != NULL);

   *binvaridx = -1;
   *aggtype = NOAGG;

   getDownlockRowIdx(matrix, aggvaridx, &rowidx, &coef);

   if( rowidx < 0 )
      return;

   assert(coef > 0);
   minact = SCIPmatrixGetRowMinActivity(matrix, rowidx);
   maxact = SCIPmatrixGetRowMaxActivity(matrix, rowidx);

   if( SCIPisInfinity(scip, -minact) || SCIPisInfinity(scip, maxact) )
      return;

   lhs = SCIPmatrixGetRowLhs(matrix, rowidx);
   ub = SCIPmatrixGetColUb(matrix, aggvaridx);

   /* search for appropriate binary variables */
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, rowidx);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, rowidx);
   valpnt = SCIPmatrixGetRowValPtr(matrix, rowidx);
   for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
   {
      SCIP_VAR* var;

      if( *rowpnt == aggvaridx )
         continue;

      var = SCIPmatrixGetVar(matrix, *rowpnt);

      /* avoid cases where the binary variable has lb=ub=1 or lb=ub=0 */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY &&
          SCIPmatrixGetColLb(matrix, *rowpnt) < 0.5 &&
          SCIPmatrixGetColUb(matrix, *rowpnt) > 0.5 )
      {
         SCIP_Real bincoef;

         bincoef = *valpnt;

         if( bincoef < 0 )
         {
            /* binvar = 0 implies that the constraint is redundant */
            if( SCIPisGE(scip, minact-bincoef, lhs) )
            {
               /* binvar = 1 implies that aggvar = ub */
               SCIP_Real bnd;
               bnd = (lhs - maxact + coef*ub - bincoef) / coef;
               if( SCIPisGE(scip, bnd, ub) )
               {
                  *binvaridx = *rowpnt;
                  *aggtype = BIN0LBOUND;
                  break;
               }
            }
         }

         if( bincoef > 0 )
         {
            /* binvar = 1 implies that the constraint is redundant */
            if( SCIPisGE(scip, minact+bincoef, lhs) )
            {
               /* binvar = 0 implies that aggvar = ub */
               SCIP_Real bnd;
               bnd = (lhs - maxact + coef*ub + bincoef) / coef;
               if( SCIPisGE(scip, bnd, ub) )
               {
                  *binvaridx = *rowpnt;
                  *aggtype = BIN0UBOUND;
                  break;
               }
            }
         }
      }
   }
}

/** find variable aggregations for uplock case */
static
SCIP_RETCODE findUplockAggregations(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int*                  nvaragg,            /**< number of redundant variables */
   AGGRTYPE*             aggtypes,           /**< type of aggregations (in same order as variables in matrix) */
   SCIP_VAR**            binvars             /**< pointers to the binary variables (in same order as variables in matrix) */
   )
{
   int nvars;
   int i;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(nvaragg != NULL);
   assert(aggtypes != NULL);
   assert(binvars != NULL);

   nvars = SCIPmatrixGetNColumns(matrix);

   for( i = 0; i < nvars; i++ )
   {
      /* column has only one uplock which keeps it from being fixed by duality fixing */
      if( SCIPmatrixGetColNUplocks(matrix, i) == 1 &&
         SCIPisLE(scip, SCIPvarGetObj(SCIPmatrixGetVar(matrix, i)), 0.0) )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPmatrixGetColLb(matrix, i);
         ub = SCIPmatrixGetColUb(matrix, i);
         assert(lb == SCIPvarGetLbGlobal(SCIPmatrixGetVar(matrix, i))); /*lint !e777*/
         assert(ub == SCIPvarGetUbGlobal(SCIPmatrixGetVar(matrix, i))); /*lint !e777*/

         /* the variable needs to have finite bounds to allow an agregation */
         if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
         {
            int binvaridx;
            AGGRTYPE aggtype;

            getBinVarIdxInUplockRow(scip, matrix, i, &binvaridx, &aggtype);

            if( binvaridx >= 0 )
            {
               aggtypes[i] = aggtype;
               binvars[i] = SCIPmatrixGetVar(matrix, binvaridx);
               (*nvaragg)++;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** find variable aggregations for downlock case */
static
SCIP_RETCODE findDownlockAggregations(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int*                  nvaragg,            /**< number of redundant variables */
   AGGRTYPE*             aggtypes,           /**< type of aggregations (in same order as variables in matrix) */
   SCIP_VAR**            binvars             /**< pointers to the binary variables (in same order as variables in matrix) */
   )
{
   int nvars;
   int i;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(nvaragg != NULL);
   assert(aggtypes != NULL);
   assert(binvars != NULL);

   nvars = SCIPmatrixGetNColumns(matrix);

   for( i = 0; i < nvars; i++ )
   {
      /* column has only one downlock which keeps it from being fixed by duality fixing;
       * only handle variable if it was not yet aggregated due to a single uplock
       */
      if( SCIPmatrixGetColNDownlocks(matrix, i) == 1 &&
         SCIPisGE(scip, SCIPvarGetObj(SCIPmatrixGetVar(matrix, i)), 0.0) &&
         aggtypes[i] == NOAGG )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPmatrixGetColLb(matrix, i);
         ub = SCIPmatrixGetColUb(matrix, i);
         assert(lb == SCIPvarGetLbGlobal(SCIPmatrixGetVar(matrix, i))); /*lint !e777*/
         assert(ub == SCIPvarGetUbGlobal(SCIPmatrixGetVar(matrix, i))); /*lint !e777*/

         /* the variable needs to have finite bounds to allow an agregation */
         if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
         {
            int binvaridx;
            AGGRTYPE aggtype;
            getBinVarIdxInDownlockRow(scip, matrix, i, &binvaridx, &aggtype);

            if( binvaridx >= 0 )
            {
               aggtypes[i] = aggtype;
               binvars[i] = SCIPmatrixGetVar(matrix, binvaridx);
               (*nvaragg)++;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualagg)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   if( SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   /* we only work on pure MIPs currently */
   if( initialized && complete )
   {
      AGGRTYPE* aggtypes;
      SCIP_VAR** binvars;
      int nvaragg;
      int ncols;

      ncols = SCIPmatrixGetNColumns(matrix);
      nvaragg = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &aggtypes, ncols) );
      BMSclearMemoryArray(aggtypes, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &binvars, ncols) );
      SCIPdebug( BMSclearMemoryArray(binvars, ncols) );

      /* search for aggregations */
      SCIP_CALL( findUplockAggregations(scip, matrix, &nvaragg, aggtypes, binvars) );
      SCIP_CALL( findDownlockAggregations(scip, matrix, &nvaragg, aggtypes, binvars) );

      /* apply aggregations, if we found any */
      if( nvaragg > 0 )
      {
         int v;

         for( v = 0; v < ncols; v++ )
         {
            if( aggtypes[v] != NOAGG )
            {
               SCIP_Bool infeasible;
               SCIP_Bool redundant;
               SCIP_Bool aggregated;
               SCIP_Real ub;
               SCIP_Real lb;

               ub = SCIPmatrixGetColUb(matrix, v);
               lb = SCIPmatrixGetColLb(matrix, v);

               /* aggregate variable */
               assert(binvars[v] != NULL);
               if( aggtypes[v] == BIN0UBOUND )
               {
                  SCIP_CALL( SCIPaggregateVars(scip, SCIPmatrixGetVar(matrix, v), binvars[v], 1.0, ub-lb,
                        ub, &infeasible, &redundant, &aggregated) );
               }
               else
               {
                  assert(aggtypes[v] == BIN0LBOUND);
                  SCIP_CALL( SCIPaggregateVars(scip, SCIPmatrixGetVar(matrix, v), binvars[v], 1.0, lb-ub,
                        lb, &infeasible, &redundant, &aggregated) );
               }

               /* infeasible aggregation */
               if( infeasible )
               {
                  SCIPdebugMsg(scip, " -> infeasible aggregation\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }

               if( aggregated )
                  (*naggrvars)++;
            }
         }

         /* set result pointer */
         if( (*naggrvars) > 0 )
            *result = SCIP_SUCCESS;
      }

      SCIPfreeBufferArray(scip, &binvars);
      SCIPfreeBufferArray(scip, &aggtypes);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the dualagg presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualagg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDualagg, NULL) );

   return SCIP_OKAY;
}
