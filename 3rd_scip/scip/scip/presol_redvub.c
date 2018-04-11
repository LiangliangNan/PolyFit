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

/**@file   presol_redvub.c
 * @brief  remove redundant variable upper bound constraints
 * @author Dieter Weninger
 *
 * This presolver looks for dominating variable bound constraints
 * on the same continuous variable and discards them. For example let x be a
 * continuous variable and y, y' are binary variables. In addition, let two variable
 * upper bound constraints ax - by <= e and cx - dy' <= f are given. If
 * ax - by <= e implies cx - dy' <= f, then we can remove the second constraint
 * and substitute/aggregate y' := y. The same can be done with variable lower
 * bound constraints.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>


#include "scip/pub_matrix.h"
#include "presol_redvub.h"

#define PRESOL_NAME            "redvub"
#define PRESOL_DESC            "detect redundant variable bound constraints"
#define PRESOL_PRIORITY        -9000000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               0     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define MAXPAIRCOMP                 1000     /**< maximal number of pairwise comparisons */

/*
 * Local methods
 */

/** verify if the constraint is a variable upper bound constraint */
static
SCIP_Bool isVub(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row,                /**< row index */
   SCIP_Real*            lowthreshold,       /**< low switching threshold */
   SCIP_Real*            highthreshold,      /**< high switching threshold */
   int*                  conidx,             /**< variable index of continuous variable */
   int*                  binidx              /**< variable index of binary variable */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   SCIP_Bool isvub;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < SCIPmatrixGetNRows(matrix));
   assert(lowthreshold != NULL);
   assert(highthreshold != NULL);
   assert(conidx != NULL);
   assert(binidx != NULL);

   isvub = FALSE;

   if( SCIPmatrixGetRowNNonzs(matrix, row) == 2 && SCIPmatrixIsRowRhsInfinity(matrix, row) )
   {
      SCIP_VARTYPE type1;
      SCIP_VARTYPE type2;
      int idx1;
      int idx2;
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real val1;
      SCIP_Real val2;

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
      valpnt = SCIPmatrixGetRowValPtr(matrix, row);

      idx1 = *rowpnt;
      val1 = *valpnt;
      var1 = SCIPmatrixGetVar(matrix, idx1);
      type1 = SCIPvarGetType(var1);

      rowpnt++;
      valpnt++;

      idx2 = *rowpnt;
      val2 = *valpnt;
      var2 = SCIPmatrixGetVar(matrix, idx2);
      type2 = SCIPvarGetType(var2);

      /* we claim that the vub has the structure ax + cy >= b
       * with a<0, c>0, x continuous, x>=0, y binary and obj(y)>=0
       */
      if( (type1 == SCIP_VARTYPE_CONTINUOUS && type2 == SCIP_VARTYPE_BINARY)
         && val1 < 0.0 && val2 > 0.0 && SCIPisGE(scip, SCIPvarGetLbGlobal(var1), 0.0)
         && SCIPisGE(scip, SCIPvarGetObj(var2), 0.0) )
      {
         *lowthreshold = SCIPmatrixGetRowLhs(matrix, row) / val1;
         *highthreshold = (SCIPmatrixGetRowLhs(matrix, row) - val2) / val1;
         *conidx = idx1;
         *binidx = idx2;
         isvub = TRUE;
      }
      else if( (type1 == SCIP_VARTYPE_BINARY && type2 == SCIP_VARTYPE_CONTINUOUS)
         && val1 > 0.0 && val2 < 0.0 && SCIPisGE(scip, SCIPvarGetLbGlobal(var2), 0.0)
         && SCIPisGE(scip, SCIPvarGetObj(var1), 0.0) )
      {
         *lowthreshold = SCIPmatrixGetRowLhs(matrix, row) / val2;
         *highthreshold = (SCIPmatrixGetRowLhs(matrix, row) - val1) / val2;
         *conidx = idx2;
         *binidx = idx1;
         isvub = TRUE;
      }
   }

   return isvub;
}

/** verify if the constraint is a variable lower bound constraint */
static
SCIP_Bool isVlb(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row,                /**< row index */
   SCIP_Real*            lowthreshold,       /**< low switching threshold */
   SCIP_Real*            highthreshold,      /**< high switching threshold */
   int*                  conidx,             /**< variable index of continuous variable */
   int*                  binidx              /**< variable index of binary variable */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   SCIP_Bool isvlb;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < SCIPmatrixGetNRows(matrix));
   assert(lowthreshold != NULL);
   assert(highthreshold != NULL);
   assert(conidx != NULL);
   assert(binidx != NULL);

   isvlb = FALSE;

   if( SCIPmatrixGetRowNNonzs(matrix, row) == 2 && SCIPmatrixIsRowRhsInfinity(matrix, row) )
   {
      SCIP_VARTYPE type1;
      SCIP_VARTYPE type2;
      int idx1;
      int idx2;
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real val1;
      SCIP_Real val2;

      rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
      valpnt = SCIPmatrixGetRowValPtr(matrix, row);

      idx1 = *rowpnt;
      val1 = *valpnt;
      var1 = SCIPmatrixGetVar(matrix, idx1);
      type1 = SCIPvarGetType(var1);

      rowpnt++;
      valpnt++;

      idx2 = *rowpnt;
      val2 = *valpnt;
      var2 = SCIPmatrixGetVar(matrix, idx2);
      type2 = SCIPvarGetType(var2);

      /* we claim that the vlb has the structure ax + cy >= b
       * with a>0, c<0, x continuous, x>=0, y binary and obj(y)>=0
       */
      if( (type1 == SCIP_VARTYPE_CONTINUOUS && type2 == SCIP_VARTYPE_BINARY)
         && val1 > 0.0 && val2 < 0.0 && SCIPisGE(scip, SCIPvarGetLbGlobal(var1), 0.0)
         && SCIPisGE(scip, SCIPvarGetObj(var2), 0.0) )
      {
         *lowthreshold = SCIPmatrixGetRowLhs(matrix, row) / val1;
         *highthreshold = (SCIPmatrixGetRowLhs(matrix, row) - val2) / val1;
         *conidx = idx1;
         *binidx = idx2;
         isvlb = TRUE;
      }
      else if( (type1 == SCIP_VARTYPE_BINARY && type2 == SCIP_VARTYPE_CONTINUOUS)
         && val1 < 0.0 && val2 > 0.0 && SCIPisGE(scip, SCIPvarGetLbGlobal(var2), 0.0)
         && SCIPisGE(scip, SCIPvarGetObj(var1), 0.0) )
      {
         *lowthreshold = SCIPmatrixGetRowLhs(matrix, row) / val2;
         *highthreshold = (SCIPmatrixGetRowLhs(matrix, row) - val1) / val2;
         *conidx = idx2;
         *binidx = idx1;
         isvlb = TRUE;
      }
   }

   return isvlb;
}


/** searches for variable upper bound constraints on the same continuous variable with a dominance relation */
static
SCIP_RETCODE detectDominatingVubs(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   nvubs,              /**< number of vubs */
   int*                  vubs,               /**< row indices of the vubs */
   SCIP_Real*            lowthresholds,      /**< low switching thresholds */
   SCIP_Real*            highthresholds,     /**< high switching thresholds */
   int*                  conidxs,            /**< variable indexes of continuous variable */
   int*                  binidxs,            /**< variable indexes of binary variable */
   int*                  nvaragg,            /**< number of variables for aggregation */
   SCIP_Bool*            isvartoagg,         /**< flags indicating if variable could be aggregated */
   SCIP_VAR**            aggvars,            /**< pointers to the variables by which the aggregation should be done */
   int*                  ndeletecons,        /**< number of deleteable constraints */
   SCIP_Bool*            deletecons          /**< flags which constraints could be deleted */
   )
{
   int i;
   int j;
   SCIP_Bool uselinearscan;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(vubs != NULL);
   assert(nvubs >= 2);
   assert(lowthresholds != NULL);
   assert(highthresholds != NULL);
   assert(conidxs != NULL);
   assert(binidxs != NULL);
   assert(nvaragg != NULL);
   assert(isvartoagg != NULL);
   assert(aggvars != NULL);
   assert(ndeletecons != NULL);
   assert(deletecons != NULL);

   /* use pairwise comparison only for a small number of vub constraints */
   if( nvubs >= MAXPAIRCOMP )
      uselinearscan = TRUE;
   else
      uselinearscan = FALSE;

   for( i = 0; i < nvubs; i++ )
   {
      for( j = i+1; j < nvubs; j++ )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;

         if( !SCIPisEQ(scip, lowthresholds[i], lowthresholds[j]) )
            continue;

         var1 = SCIPmatrixGetVar(matrix, binidxs[i]);
         var2 = SCIPmatrixGetVar(matrix, binidxs[j]);

         /* make sure that the binary variables switch synchronously */
         if((SCIPmatrixGetColNDownlocks(matrix, binidxs[j]) == 1 &&
               SCIPmatrixGetColNDownlocks(matrix, binidxs[i]) == 1 &&
               SCIPmatrixGetColNUplocks(matrix, binidxs[j]) == 0 &&
               SCIPmatrixGetColNUplocks(matrix, binidxs[i]) == 0) ||
            (SCIPvarsHaveCommonClique(var1, FALSE, var2, TRUE, TRUE) &&
               SCIPvarsHaveCommonClique(var1, TRUE, var2, FALSE, TRUE)) )
         {

            if( SCIPisLE(scip, highthresholds[i], highthresholds[j]) )
            {
#ifdef SCIP_DEBUG
               SCIPdebugMsg(scip, "Aggregate variable %s by %s\n", SCIPvarGetName(var2), SCIPvarGetName(var1));
               SCIPdebugMsg(scip, "Delete variable upper bound constraint:\n");
               SCIP_CALL( SCIPprintCons(scip, SCIPmatrixGetCons(matrix, vubs[j]), NULL));
               SCIPinfoMessage(scip, NULL, "\n");
#endif

               isvartoagg[binidxs[j]] = TRUE;
               aggvars[binidxs[j]] = SCIPmatrixGetVar(matrix, binidxs[i]);
               (*nvaragg)++;

               deletecons[vubs[j]] = TRUE;
               (*ndeletecons)++;
            }
            else
            {
               assert(SCIPisGT(scip, highthresholds[i], highthresholds[j]));
#ifdef SCIP_DEBUG
               SCIPdebugMsg(scip, "Aggregate variable %s by %s\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
               SCIPdebugMsg(scip, "Delete variable upper bound constraint:\n");
               SCIP_CALL( SCIPprintCons(scip, SCIPmatrixGetCons(matrix, vubs[i]), NULL));
               SCIPinfoMessage(scip, NULL, "\n");
#endif

               isvartoagg[binidxs[i]] = TRUE;
               aggvars[binidxs[i]] = SCIPmatrixGetVar(matrix, binidxs[j]);
               (*nvaragg)++;

               deletecons[vubs[i]] = TRUE;
               (*ndeletecons)++;
            }
         }

         if( uselinearscan )
            break;
      }
   }

   return SCIP_OKAY;
}

/** searches for variable lower bound constraints on the same continuous variable with a dominance relation */
static
SCIP_RETCODE detectDominatingVlbs(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   nvlbs,              /**< number of vlbs */
   int*                  vlbs,               /**< row indices of the vlbs */
   SCIP_Real*            lowthresholds,      /**< low switching thresholds */
   SCIP_Real*            highthresholds,     /**< high switching thresholds */
   int*                  conidxs,            /**< variable indexes of continuous variable */
   int*                  binidxs,            /**< variable indexes of binary variable */
   int*                  nvaragg,            /**< number of variables for aggregation */
   SCIP_Bool*            isvartoagg,         /**< flags indicating if variable could be aggregated */
   SCIP_VAR**            aggvars,            /**< pointers to the variables by which the aggregation should be done */
   int*                  ndeletecons,        /**< number of deleteable constraints */
   SCIP_Bool*            deletecons          /**< flags which constraints could be deleted */

   )
{
   int i;
   int j;
   SCIP_Bool uselinearscan;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(vlbs != NULL);
   assert(nvlbs >= 2);
   assert(lowthresholds != NULL);
   assert(highthresholds != NULL);
   assert(conidxs != NULL);
   assert(binidxs != NULL);
   assert(nvaragg != NULL);
   assert(isvartoagg != NULL);
   assert(aggvars != NULL);
   assert(ndeletecons != NULL);
   assert(deletecons != NULL);

   /* use pairwise comparison only for a small number of vlb constraints */
   if( nvlbs >= MAXPAIRCOMP )
      uselinearscan = TRUE;
   else
      uselinearscan = FALSE;

   for( i = 0; i < nvlbs; i++ )
   {
      for( j = i+1; j < nvlbs; j++ )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;

         if( !SCIPisEQ(scip, lowthresholds[i], lowthresholds[j]) )
            continue;

         var1 = SCIPmatrixGetVar(matrix, binidxs[i]);
         var2 = SCIPmatrixGetVar(matrix, binidxs[j]);

         /* make sure that the binary variables switch synchronously */
         if((SCIPmatrixGetColNUplocks(matrix, binidxs[j]) == 1 &&
               SCIPmatrixGetColNUplocks(matrix, binidxs[i]) == 1 &&
               SCIPmatrixGetColNDownlocks(matrix, binidxs[j]) == 0 &&
               SCIPmatrixGetColNDownlocks(matrix, binidxs[i]) == 0) ||
            (SCIPvarsHaveCommonClique(var1, FALSE, var2, TRUE, TRUE) &&
               SCIPvarsHaveCommonClique(var1, TRUE, var2, FALSE, TRUE)) )
         {

            if( SCIPisGE(scip, highthresholds[i], highthresholds[j]) )
            {
#ifdef SCIP_DEBUG
               SCIPdebugMsg(scip, "Aggregate variable %s by %s\n", SCIPvarGetName(var2), SCIPvarGetName(var1));
               SCIPdebugMsg(scip, "Delete variable lower bound constraint:\n");
               SCIP_CALL( SCIPprintCons(scip, SCIPmatrixGetCons(matrix, vlbs[j]), NULL));
               SCIPinfoMessage(scip, NULL, "\n");
#endif

               isvartoagg[binidxs[j]] = TRUE;
               aggvars[binidxs[j]] = SCIPmatrixGetVar(matrix, binidxs[i]);
               (*nvaragg)++;

               deletecons[vlbs[j]] = TRUE;
               (*ndeletecons)++;
            }
            else
            {
               assert(SCIPisLT(scip, highthresholds[i], highthresholds[j]));
#ifdef SCIP_DEBUG
               SCIPdebugMsg(scip, "Aggregate variable %s by %s\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
               SCIPdebugMsg(scip, "Delete variable lower bound constraint:\n");
               SCIP_CALL( SCIPprintCons(scip, SCIPmatrixGetCons(matrix, vlbs[i]), NULL));
               SCIPinfoMessage(scip, NULL, "\n");
#endif

               isvartoagg[binidxs[i]] = TRUE;
               aggvars[binidxs[i]] = SCIPmatrixGetVar(matrix, binidxs[j]);
               (*nvaragg)++;

               deletecons[vlbs[i]] = TRUE;
               (*ndeletecons)++;
            }
         }

         if( uselinearscan )
            break;
      }
   }

   return SCIP_OKAY;
}

/** find variable aggregations and redundant variable bound constraints */
static
SCIP_RETCODE findVarAggrRedVbcons(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int*                  nvaragg,            /**< number of redundant variables */
   SCIP_Bool*            isvartoagg,         /**< flags indicating which variables could be substituted/aggregated */
   SCIP_VAR**            aggvars,            /**< pointers to the variables by which the aggregation should be done */
   int*                  ndeletecons,        /**< number of redundant constraints */
   SCIP_Bool*            deletecons          /**< flags indicating which constraints could be deleted */
   )
{
   int c;
   int* colpnt;
   int* colend;
   int* vbcons;
   int nvbcons;
   int ncols;
   int nrows;
   SCIP_Real* lowthresholds;
   SCIP_Real* highthresholds;
   int* conidxs;
   int* binidxs;

   ncols = SCIPmatrixGetNColumns(matrix);
   nrows = SCIPmatrixGetNRows(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &binidxs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conidxs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lowthresholds, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &highthresholds, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbcons, nrows) );

   for( c = 0; c < ncols; c++ )
   {
      SCIP_VAR* var;

      var = SCIPmatrixGetVar(matrix, c);

      if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         continue;

      /* search vubs per variable */
      nvbcons = 0;
      colpnt = SCIPmatrixGetColIdxPtr(matrix, c);
      colend = colpnt + SCIPmatrixGetColNNonzs(matrix, c);
      for( ; (colpnt < colend); colpnt++ )
      {
         SCIP_Real lowthreshold;
         SCIP_Real highthreshold;
         int conidx;
         int binidx;

         if( isVub(scip, matrix, *colpnt, &lowthreshold, &highthreshold, &conidx, &binidx) )
         {
            vbcons[nvbcons] = *colpnt;
            lowthresholds[nvbcons] = lowthreshold;
            highthresholds[nvbcons] = highthreshold;
            conidxs[nvbcons] = conidx;
            binidxs[nvbcons] = binidx;
            nvbcons++;
         }
      }
      if( nvbcons >= 2 )
      {
         SCIP_CALL( detectDominatingVubs(scip, matrix, nvbcons, vbcons,
               lowthresholds, highthresholds, conidxs, binidxs,
               nvaragg, isvartoagg, aggvars, ndeletecons, deletecons) );
      }

      /* search vlbs per variable */
      nvbcons = 0;
      colpnt = SCIPmatrixGetColIdxPtr(matrix, c);
      colend = colpnt + SCIPmatrixGetColNNonzs(matrix, c);
      for( ; (colpnt < colend); colpnt++ )
      {
         SCIP_Real lowthreshold;
         SCIP_Real highthreshold;
         int conidx;
         int binidx;

         if( isVlb(scip, matrix, *colpnt, &lowthreshold, &highthreshold, &conidx, &binidx) )
         {
            vbcons[nvbcons] = *colpnt;
            lowthresholds[nvbcons] = lowthreshold;
            highthresholds[nvbcons] = highthreshold;
            conidxs[nvbcons] = conidx;
            binidxs[nvbcons] = binidx;
            nvbcons++;
         }
      }
      if( nvbcons >= 2 )
      {
         SCIP_CALL( detectDominatingVlbs(scip, matrix, nvbcons, vbcons,
               lowthresholds, highthresholds, conidxs, binidxs,
               nvaragg, isvartoagg, aggvars, ndeletecons, deletecons) );
      }
   }

   SCIPfreeBufferArray(scip, &vbcons);
   SCIPfreeBufferArray(scip, &highthresholds);
   SCIPfreeBufferArray(scip, &lowthresholds);
   SCIPfreeBufferArray(scip, &conidxs);
   SCIPfreeBufferArray(scip, &binidxs);

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecRedvub)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPgetNContVars(scip) == 0 || SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      int nvaragg;
      SCIP_Bool* isvartoagg;
      int ndeletecons;
      SCIP_Bool* deletecons;
      SCIP_VAR** aggvars;
      int ncols;
      int nrows;

      ncols = SCIPmatrixGetNColumns(matrix);
      nrows = SCIPmatrixGetNRows(matrix);

      nvaragg = 0;
      ndeletecons = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &isvartoagg, ncols) );
      BMSclearMemoryArray(isvartoagg, ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &deletecons, nrows) );
      BMSclearMemoryArray(deletecons, nrows);

      SCIP_CALL( SCIPallocBufferArray(scip, &aggvars, ncols) );
      BMSclearMemoryArray(aggvars, ncols);

      SCIP_CALL( findVarAggrRedVbcons(scip, matrix, &nvaragg, isvartoagg, aggvars, &ndeletecons, deletecons) );

      if( nvaragg > 0 )
      {
         int v;
         for( v = 0; v < ncols; v++ )
         {
            if( isvartoagg[v] )
            {
               SCIP_Bool infeasible;
               SCIP_Bool redundant;
               SCIP_Bool aggregated;

               /* substitute/aggregate binary variable */
               assert(aggvars[v] != NULL);
               SCIP_CALL( SCIPaggregateVars(scip, SCIPmatrixGetVar(matrix,v), aggvars[v], 1.0, -1.0,
                     0.0, &infeasible, &redundant, &aggregated) );

               if( infeasible )
               {
                  SCIPdebugMsg(scip, " -> infeasible aggregation\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }

               if( aggregated )
               {
                  (*naggrvars)++;

                  /* set result pointer */
                  if( (*result) == SCIP_DIDNOTFIND )
                     *result = SCIP_SUCCESS;
               }
            }
         }
      }

      if( ndeletecons > 0 )
      {
         int r;
         for( r = 0; r < nrows; r++ )
         {
            if( deletecons[r] )
            {
               SCIP_CONS* cons;

               /* remove redundant variable bound constraint */
               cons = SCIPmatrixGetCons(matrix, r);
               SCIP_CALL( SCIPdelCons(scip, cons) );

               (*ndelconss)++;

               /* set result pointer */
               if( (*result) == SCIP_DIDNOTFIND )
                  *result = SCIP_SUCCESS;
            }
         }
      }

      SCIPfreeBufferArray(scip, &aggvars);
      SCIPfreeBufferArray(scip, &deletecons);
      SCIPfreeBufferArray(scip, &isvartoagg);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the redvub presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolRedvub(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecRedvub, NULL) );

   return SCIP_OKAY;
}
