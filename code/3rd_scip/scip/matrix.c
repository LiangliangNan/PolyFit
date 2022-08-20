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

/**@file   matrix.c
 * @brief  methods for MIP matrix data structure
 * @author Dieter Weninger
 * @author Gerald Gamrath
 *
 * The MIP matrix is organized as sparse data structure in row and
 * and column major format.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/struct_matrix.h"
#include "scip/pub_matrix.h"

#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"

/*
 * private functions
 */

/** transforms given variables, scalars and constant to the corresponding active variables, scalars and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR***           vars,               /**< vars array to get active variables for */
   SCIP_Real**           scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant            /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
   )
{
   int requiredsize;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(scalars != NULL);
   assert(*vars != NULL);
   assert(*scalars != NULL);
   assert(nvars != NULL);
   assert(constant != NULL);

   SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

   if( requiredsize > *nvars )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );

      /* call function a second time with enough memory */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
      assert(requiredsize <= *nvars);
   }

   return SCIP_OKAY;
}

/** add one row to the constraint matrix */
static
SCIP_RETCODE addRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   SCIP_VAR**            vars,               /**< variables of this row */
   SCIP_Real*            vals,               /**< coefficients of this row */
   int                   nvars,              /**< number of variables of this row */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   int                   maxnnonzsmem,       /**< maximal number of fillable elements */
   SCIP_Bool*            rowadded            /**< flag indicating if constraint was added to matrix */
   )
{
   int j;
   int probindex;
   int rowidx;
   SCIP_Real factor;
   SCIP_Bool rangedorequality;

   assert(vars != NULL);
   assert(vals != NULL);

   rowidx = matrix->nrows;
   rangedorequality = FALSE;

   if( SCIPisInfinity(scip, -lhs) )
   {
      factor = -1.0;
      matrix->lhs[rowidx] = -rhs;
      matrix->rhs[rowidx] = SCIPinfinity(scip);
      matrix->isrhsinfinite[rowidx] = TRUE;
   }
   else
   {
      factor = 1.0;
      matrix->lhs[rowidx] = lhs;
      matrix->rhs[rowidx] = rhs;
      matrix->isrhsinfinite[rowidx] = SCIPisInfinity(scip, matrix->rhs[rowidx]);

      if( !SCIPisInfinity(scip, rhs) )
         rangedorequality = TRUE;
   }

   if(SCIPisInfinity(scip, -matrix->lhs[rowidx]))
   {
      /* ignore redundant constraint */
      *rowadded = FALSE;
      return SCIP_OKAY;
   }

   matrix->rowmatbeg[rowidx] = matrix->nnonzs;

   /* = or ranged */
   if( rangedorequality )
   {
      assert(factor > 0);

      for( j = 0; j < nvars; j++ )
      {
         assert(maxnnonzsmem > matrix->nnonzs);

         /* ignore variables with very small coefficients */
         if( SCIPisZero(scip, vals[j]) )
            continue;

         matrix->rowmatval[matrix->nnonzs] = factor * vals[j];
         probindex = SCIPvarGetProbindex(vars[j]);
         assert(matrix->vars[probindex] == vars[j]);

         matrix->nuplocks[probindex]++;
         matrix->ndownlocks[probindex]++;

         assert(0 <= probindex && probindex < matrix->ncols);
         matrix->rowmatind[matrix->nnonzs] = probindex;

         (matrix->nnonzs)++;
      }
   }
   /* >= or <= */
   else
   {
      for( j = 0; j < nvars; j++ )
      {
         assert(maxnnonzsmem > matrix->nnonzs);

         /* ignore variables with very small coefficients */
         if( SCIPisZero(scip, vals[j]) )
            continue;

         /* due to the factor, <= constraints will be transfered to >= */
         matrix->rowmatval[matrix->nnonzs] = factor * vals[j];
         probindex = SCIPvarGetProbindex(vars[j]);
         assert(matrix->vars[probindex] == vars[j]);

         if( matrix->rowmatval[matrix->nnonzs] > 0 )
            matrix->ndownlocks[probindex]++;
         else
         {
            assert(matrix->rowmatval[matrix->nnonzs] < 0);
            matrix->nuplocks[probindex]++;
         }

         assert(0 <= probindex && probindex < matrix->ncols);
         matrix->rowmatind[matrix->nnonzs] = probindex;

         (matrix->nnonzs)++;
      }
   }

   matrix->rowmatcnt[rowidx] = matrix->nnonzs - matrix->rowmatbeg[rowidx];

   ++(matrix->nrows);
   *rowadded = TRUE;

   return SCIP_OKAY;
}

/** add one constraint to matrix */
static
SCIP_RETCODE addConstraint(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   SCIP_VAR**            vars,               /**< variables of this constraint */
   SCIP_Real*            vals,               /**< variable coefficients of this constraint */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   int                   maxnnonzsmem,       /**< maximal number of fillable elements */
   SCIP_Bool*            rowadded            /**< flag indicating of row was added to matrix */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   SCIP_Real activeconstant;
   int nactivevars;
   int v;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(vars != NULL || nvars == 0);
   assert(SCIPisLE(scip, lhs, rhs));
   assert(rowadded != NULL);

   *rowadded = FALSE;

   /* constraint is redundant */
   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   /* we do not add empty constraints to the matrix */
   if( nvars == 0 )
      return SCIP_OKAY;

   activevars = NULL;
   activevals = NULL;
   nactivevars = nvars;
   activeconstant = 0.0;

   /* duplicate variable and value array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   if( vals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; v++ )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, &activevars, &activevals, &nactivevars, &activeconstant) );

   /* adapt left and right hand side */
   if( !SCIPisInfinity(scip, -lhs) )
      lhs -= activeconstant;
   if( !SCIPisInfinity(scip, rhs) )
      rhs -= activeconstant;

   /* add single row to matrix */
   if( nactivevars > 0 )
   {
      SCIP_CALL( addRow(scip, matrix, activevars, activevals, nactivevars, lhs, rhs, maxnnonzsmem, rowadded) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevals);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}

/** transform row major format into column major format */
static
SCIP_RETCODE setColumnMajorFormat(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX*          matrix              /**< constraint matrix */
   )
{
   int colidx;
   int i;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int* fillidx;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(matrix->colmatval != NULL);
   assert(matrix->colmatind != NULL);
   assert(matrix->colmatbeg != NULL);
   assert(matrix->colmatcnt != NULL);
   assert(matrix->rowmatval != NULL);
   assert(matrix->rowmatind != NULL);
   assert(matrix->rowmatbeg != NULL);
   assert(matrix->rowmatcnt != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &fillidx, matrix->ncols) );
   BMSclearMemoryArray(fillidx, matrix->ncols);
   BMSclearMemoryArray(matrix->colmatcnt, matrix->ncols);

   for( i = 0; i < matrix->nrows; i++ )
   {
      rowpnt = matrix->rowmatind + matrix->rowmatbeg[i];
      rowend = rowpnt + matrix->rowmatcnt[i];
      for( ; rowpnt < rowend; rowpnt++ )
      {
         colidx = *rowpnt;
         (matrix->colmatcnt[colidx])++;
      }
   }

   matrix->colmatbeg[0] = 0;
   for( i = 0; i < matrix->ncols-1; i++ )
   {
      matrix->colmatbeg[i+1] = matrix->colmatbeg[i] + matrix->colmatcnt[i];
   }

   for( i = 0; i < matrix->nrows; i++ )
   {
      rowpnt = matrix->rowmatind + matrix->rowmatbeg[i];
      rowend = rowpnt + matrix->rowmatcnt[i];
      valpnt = matrix->rowmatval + matrix->rowmatbeg[i];

      for( ; rowpnt < rowend; rowpnt++, valpnt++ )
      {
         assert(*rowpnt < matrix->ncols);
         colidx = *rowpnt;
         matrix->colmatval[matrix->colmatbeg[colidx] + fillidx[colidx]] = *valpnt;
         matrix->colmatind[matrix->colmatbeg[colidx] + fillidx[colidx]] = i;
         fillidx[colidx]++;
      }
   }

   SCIPfreeBufferArray(scip, &fillidx);

   return SCIP_OKAY;
}

/** calculate min/max activity per row */
static
SCIP_RETCODE calcActivityBounds(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX*          matrix              /**< constraint matrix */
   )
{
   SCIP_Real val;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int col;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);

   for( row = 0; row < matrix->nrows; row++ )
   {
      matrix->minactivity[row] = 0;
      matrix->maxactivity[row] = 0;
      matrix->minactivityneginf[row] = 0;
      matrix->minactivityposinf[row] = 0;
      matrix->maxactivityneginf[row] = 0;
      matrix->maxactivityposinf[row] = 0;

      rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
      rowend = rowpnt + matrix->rowmatcnt[row];
      valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

      for( ; rowpnt < rowend; rowpnt++, valpnt++ )
      {
         /* get column index */
         col = *rowpnt;

         /* get variable coefficient */
         val = *valpnt;
         assert(!SCIPisZero(scip, val));

         assert(matrix->ncols > col);

         assert(!SCIPisInfinity(scip, matrix->lb[col]));
         assert(!SCIPisInfinity(scip, -matrix->ub[col]));

         /* positive coefficient */
         if( val > 0.0 )
         {
            if( SCIPisInfinity(scip, matrix->ub[col]) )
               matrix->maxactivityposinf[row]++;
            else
               matrix->maxactivity[row] += val * matrix->ub[col];

            if( SCIPisInfinity(scip, -matrix->lb[col]) )
               matrix->minactivityneginf[row]++;
            else
               matrix->minactivity[row] += val * matrix->lb[col];
         }
         /* negative coefficient */
         else
         {
            if( SCIPisInfinity(scip, -matrix->lb[col]) )
               matrix->maxactivityneginf[row]++;
            else
               matrix->maxactivity[row] += val * matrix->lb[col];

            if( SCIPisInfinity(scip, matrix->ub[col]) )
               matrix->minactivityposinf[row]++;
            else
               matrix->minactivity[row] += val * matrix->ub[col];
         }
      }

      /* consider infinite bound contributions for the activities */
      if( matrix->maxactivityneginf[row] + matrix->maxactivityposinf[row] > 0 )
         matrix->maxactivity[row] = SCIPinfinity(scip);

      if( matrix->minactivityneginf[row] + matrix->minactivityposinf[row] > 0 )
         matrix->minactivity[row] = -SCIPinfinity(scip);

   }

   return SCIP_OKAY;
}

/*
 * public functions
 */

/** initialize matrix */
SCIP_RETCODE SCIPmatrixCreate(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX**         matrixptr,          /**< pointer to constraint matrix object to be initialized */
   SCIP_Bool*            initialized,        /**< was the initialization successful? */
   SCIP_Bool*            complete            /**< are all constraint represented within the matrix? */
   )
{
   SCIP_MATRIX* matrix;
   SCIP_CONSHDLR** conshdlrs;
   const char* conshdlrname;
   SCIP_Bool stopped;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   int nconshdlrs;
   int nconss;
   int nconssall;
   int nnonzstmp;
   int nvars;
   int c;
   int i;
   int v;
   int cnt;

   nnonzstmp = 0;

   assert(scip != NULL);
   assert(matrixptr != NULL);
   assert(initialized != NULL);
   assert(complete != NULL);

   *initialized = FALSE;
   *complete = FALSE;

   /* return if no variables or constraints are present */
   if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
      return SCIP_OKAY;

   /* loop over all constraint handlers and collect the number of checked constraints */
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   nconss = 0;
   nconssall = 0;

   for( i = 0; i < nconshdlrs; ++i )
   {
      int nconshdlrconss;

      nconshdlrconss = SCIPconshdlrGetNCheckConss(conshdlrs[i]);

      if( nconshdlrconss > 0 )
      {
         conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);

         if( (strcmp(conshdlrname, "linear") == 0) || (strcmp(conshdlrname, "setppc") == 0)
            || (strcmp(conshdlrname, "logicor") == 0) || (strcmp(conshdlrname, "knapsack") == 0)
            || (strcmp(conshdlrname, "varbound") == 0) )
         {
            /* increment number of supported constraints */
            nconss += nconshdlrconss;
         }

         /* increment number of supported and unsupported constraints */
         nconssall += nconshdlrconss;
      }
   }

   /* print warning if we have unsupported constraint types.
    * we do not abort the matrix creation process here, because
    * it makes sometimes sense to work on an incomplete
    * matrix as long as the number of interesting variable
    * uplocks or downlocks of the matrix and scip
    * are the same.
    */
   if( nconss < nconssall )
   {
      SCIPdebugMsg(scip, "Warning: milp matrix not complete!\n");
   }
   else
   {
      /* all constraints represented within the matrix */
      *complete = TRUE;
   }

   /* do nothing if we have no checked constraints */
   if( nconss == 0 )
      return SCIP_OKAY;

   stopped = FALSE;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* approximate number of nonzeros by taking for each variable the number of up- and downlocks;
    * this counts nonzeros in equalities twice, but can be at most two times as high as the exact number
    */
   for( i = nvars - 1; i >= 0; --i )
   {
      nnonzstmp += SCIPvarGetNLocksDown(vars[i]);
      nnonzstmp += SCIPvarGetNLocksUp(vars[i]);
   }

   /* do nothing if we have no entries */
   if( nnonzstmp == 0 )
      return SCIP_OKAY;

   /* build the matrix structure */
   SCIP_CALL( SCIPallocBuffer(scip, matrixptr) );
   matrix = *matrixptr;

   /* copy vars array and set number of variables */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &matrix->vars, vars, nvars) );
   matrix->ncols = nvars;

   matrix->nrows = 0;
   matrix->nnonzs = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatbeg, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatcnt, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lb, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->ub, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->nuplocks, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->ndownlocks, matrix->ncols) );

   BMSclearMemoryArray(matrix->nuplocks, matrix->ncols);
   BMSclearMemoryArray(matrix->ndownlocks, matrix->ncols);

   /* init bounds */
   for( v = 0; v < matrix->ncols; v++ )
   {
      var = matrix->vars[v];
      assert(var != NULL);

      matrix->lb[v] = SCIPvarGetLbGlobal(var);
      matrix->ub[v] = SCIPvarGetUbGlobal(var);
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatbeg, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatcnt, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lhs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rhs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->cons, nconss) );
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &matrix->isrhsinfinite, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivity, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivity, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityneginf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityposinf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityneginf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityposinf, nconss) );

   cnt = 0;

   /* loop a second time over constraints handlers and add supported constraints to the matrix */
   for( i = 0; i < nconshdlrs; ++i )
   {
      SCIP_CONS** conshdlrconss;
      int nconshdlrconss;
      SCIP_Bool rowadded;

      if( SCIPisStopped(scip) )
      {
         stopped = TRUE;
         break;
      }

      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);
      conshdlrconss = SCIPconshdlrGetCheckConss(conshdlrs[i]);
      nconshdlrconss = SCIPconshdlrGetNCheckConss(conshdlrs[i]);

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
         {
            cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsLinear(scip, cons),
                  SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons), nnonzstmp, &rowadded) );

            if(rowadded)
            {
               assert(cnt < nconss);
               matrix->cons[cnt] = cons;
               cnt++;
            }
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
         {
            SCIP_Real lhs;
            SCIP_Real rhs;

            cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            switch( SCIPgetTypeSetppc(scip, cons) )
            {
            case SCIP_SETPPCTYPE_PARTITIONING :
               lhs = 1.0;
               rhs = 1.0;
               break;
            case SCIP_SETPPCTYPE_PACKING :
               lhs = -SCIPinfinity(scip);
               rhs = 1.0;
               break;
            case SCIP_SETPPCTYPE_COVERING :
               lhs = 1.0;
               rhs = SCIPinfinity(scip);
               break;
            default:
               return SCIP_ERROR;
            }

            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsSetppc(scip, cons), NULL,
                  SCIPgetNVarsSetppc(scip, cons), lhs, rhs, nnonzstmp, &rowadded) );

            if(rowadded)
            {
               assert(cnt < nconss);
               matrix->cons[cnt] = cons;
               cnt++;
            }
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
         {
            cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsLogicor(scip, cons),
                  NULL, SCIPgetNVarsLogicor(scip, cons), 1.0, SCIPinfinity(scip), nnonzstmp, &rowadded) );

            if(rowadded)
            {
               assert(cnt < nconss);
               matrix->cons[cnt] = cons;
               cnt++;
            }
         }
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         if( nconshdlrconss > 0 )
         {
            SCIP_Real* consvals;
            int valssize;

            valssize = 100;
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, valssize) );

            for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
            {
               SCIP_Longint* weights;

               cons = conshdlrconss[c];
               assert(SCIPconsIsTransformed(cons));

               weights = SCIPgetWeightsKnapsack(scip, cons);
               nvars = SCIPgetNVarsKnapsack(scip, cons);

               if( nvars > valssize )
               {
                  valssize = (int) (1.5 * nvars);
                  SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, valssize) );
               }

               for( v = 0; v < nvars; v++ )
                  consvals[v] = (SCIP_Real)weights[v];

               SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsKnapsack(scip, cons), consvals,
                     SCIPgetNVarsKnapsack(scip, cons), -SCIPinfinity(scip),
                     (SCIP_Real)SCIPgetCapacityKnapsack(scip, cons), nnonzstmp, &rowadded) );

               if(rowadded)
               {
                  assert(cnt < nconss);
                  matrix->cons[cnt] = cons;
                  cnt++;
               }
            }

            SCIPfreeBufferArray(scip, &consvals);
         }
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         if( nconshdlrconss > 0 )
         {
            SCIP_VAR** consvars;
            SCIP_Real* consvals;

            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );
            consvals[0] = 1.0;

            for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
            {
               cons = conshdlrconss[c];
               assert(SCIPconsIsTransformed(cons));

               consvars[0] = SCIPgetVarVarbound(scip, cons);
               consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

               consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

               SCIP_CALL( addConstraint(scip, matrix, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
                     SCIPgetRhsVarbound(scip, cons), nnonzstmp, &rowadded) );

               if(rowadded)
               {
                  assert(cnt < nconss);
                  matrix->cons[cnt] = cons;
                  cnt++;
               }
            }

            SCIPfreeBufferArray(scip, &consvals);
            SCIPfreeBufferArray(scip, &consvars);
         }
      }
   }
   assert(matrix->nrows == cnt);
   assert(matrix->nrows <= nconss);
   assert(matrix->nnonzs <= nnonzstmp);

   if( !stopped )
   {
      /* calculate row activity bounds */
      SCIP_CALL( calcActivityBounds(scip, matrix) );

      /* transform row major format into column major format */
      SCIP_CALL( setColumnMajorFormat(scip, matrix) );

      *initialized = TRUE;
   }

   return SCIP_OKAY;
}


/** frees the constraint matrix */
void SCIPmatrixFree(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX**         matrix              /**< constraint matrix object */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);

   if( (*matrix) != NULL )
   {
      assert((*matrix)->colmatval != NULL);
      assert((*matrix)->colmatind != NULL);
      assert((*matrix)->colmatbeg != NULL);
      assert((*matrix)->colmatcnt != NULL);
      assert((*matrix)->lb != NULL);
      assert((*matrix)->ub != NULL);
      assert((*matrix)->nuplocks != NULL);
      assert((*matrix)->ndownlocks != NULL);

      assert((*matrix)->rowmatval != NULL);
      assert((*matrix)->rowmatind != NULL);
      assert((*matrix)->rowmatbeg != NULL);
      assert((*matrix)->rowmatcnt != NULL);
      assert((*matrix)->lhs != NULL);
      assert((*matrix)->rhs != NULL);

      SCIPfreeBufferArray(scip, &((*matrix)->maxactivityposinf));
      SCIPfreeBufferArray(scip, &((*matrix)->maxactivityneginf));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivityposinf));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivityneginf));
      SCIPfreeBufferArray(scip, &((*matrix)->maxactivity));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivity));

      SCIPfreeMemoryArray(scip, &((*matrix)->isrhsinfinite));
      SCIPfreeBufferArray(scip, &((*matrix)->cons));

      SCIPfreeBufferArray(scip, &((*matrix)->rhs));
      SCIPfreeBufferArray(scip, &((*matrix)->lhs));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatcnt));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatbeg));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatval));

      SCIPfreeBufferArray(scip, &((*matrix)->ndownlocks));
      SCIPfreeBufferArray(scip, &((*matrix)->nuplocks));
      SCIPfreeBufferArray(scip, &((*matrix)->ub));
      SCIPfreeBufferArray(scip, &((*matrix)->lb));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatcnt));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatbeg));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatval));

      (*matrix)->nrows = 0;
      (*matrix)->ncols = 0;
      (*matrix)->nnonzs = 0;

      SCIPfreeBufferArrayNull(scip, &((*matrix)->vars));

      SCIPfreeBuffer(scip, matrix);
   }
}

/** print one row of the matrix */
void SCIPmatrixPrintRow(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row                 /**< row index */
   )
{
   int* rowpnt;
   int* rowend;
   int col;
   SCIP_Real val;
   SCIP_Real* valpnt;

   rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
   rowend = rowpnt + matrix->rowmatcnt[row];
   valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

   printf("### %s: %.15g <=", SCIPconsGetName(matrix->cons[row]), matrix->lhs[row]);
   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      col = *rowpnt;
      val = *valpnt;
      if( val < 0 )
         printf(" %.15g %s [%.15g,%.15g]", val, SCIPvarGetName(matrix->vars[col]),
            SCIPvarGetLbGlobal(matrix->vars[col]), SCIPvarGetUbGlobal(matrix->vars[col]));
      else
         printf(" +%.15g %s [%.15g,%.15g]", val, SCIPvarGetName(matrix->vars[col]),
            SCIPvarGetLbGlobal(matrix->vars[col]), SCIPvarGetUbGlobal(matrix->vars[col]));
   }
   printf(" <= %.15g ###\n", matrix->rhs[row]);
}

/** detect parallel rows of matrix. rhs/lhs are ignored. */
SCIP_RETCODE SCIPmatrixGetParallelRows(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real*            scale,              /**< scale factors of rows */
   int*                  pclass              /**< parallel row classes */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* values;
   int* classsizes;
   int* pcset;
   int* colpnt;
   int* colend;
   int* rowindices;
   int* pcs;
   SCIP_Real startval;
   SCIP_Real aij;
   int startpc;
   int startk;
   int startt;
   int pcsetfill;
   int rowidx;
   int k;
   int t;
   int m;
   int i;
   int c;
   int newpclass;
   int pc;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(pclass != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &classsizes, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcset, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &values, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowindices, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcs, matrix->nrows) );

   /* init */
   BMSclearMemoryArray(scale, matrix->nrows);
   BMSclearMemoryArray(pclass, matrix->nrows);
   BMSclearMemoryArray(classsizes, matrix->nrows);
   classsizes[0] = matrix->nrows;
   pcsetfill = 0;
   for( t = 1; t < matrix->nrows; ++t )
      pcset[pcsetfill++] = t;

   /* loop over all columns */
   for( c = 0; c < matrix->ncols; ++c )
   {
      if( matrix->colmatcnt[c] == 0 )
         continue;

      colpnt = matrix->colmatind + matrix->colmatbeg[c];
      colend = colpnt + matrix->colmatcnt[c];
      valpnt = matrix->colmatval + matrix->colmatbeg[c];

      i = 0;
      for( ; (colpnt < colend); colpnt++, valpnt++ )
      {
         aij = *valpnt;
         rowidx = *colpnt;

         if( scale[rowidx] == 0.0 )
            scale[rowidx] = aij;
         assert(scale[rowidx] != 0.0);

         rowindices[i] = rowidx;
         values[i] = aij / scale[rowidx];
         pc = pclass[rowidx];
         assert(pc < matrix->nrows);

         /* update class sizes and pclass set */
         assert(classsizes[pc] > 0);
         classsizes[pc]--;
         if( classsizes[pc] == 0 )
         {
            assert(pcsetfill < matrix->nrows);
            pcset[pcsetfill++] = pc;
         }
         pcs[i] = pc;

         i++;
      }

      /* sort on the pclass values */
      if( i > 1 )
      {
         SCIPsortIntIntReal(pcs, rowindices, values, i);
      }

      k = 0;
      while( TRUE ) /*lint !e716*/
      {
         assert(k < i);
         startpc = pcs[k];
         startk = k;

         /* find pclass-sets */
         while( k < i && pcs[k] == startpc )
            k++;

         /* sort on the A values which have equal pclass values */
         if( k - startk > 1 )
            SCIPsortRealInt(&(values[startk]), &(rowindices[startk]), k - startk);

         t = 0;
         while( TRUE ) /*lint !e716*/
         {
            assert(startk + t < i);
            startval = values[startk + t];
            startt = t;

            /* find A-sets */
            while( t < k - startk && SCIPisEQ(scip, startval, values[startk + t]) )
               t++;

            /* get new pclass */
            newpclass = pcset[0];
            assert(pcsetfill > 0);
            pcset[0] = pcset[--pcsetfill];

            /* renumbering */
            for( m = startk + startt; m < startk + t; m++ )
            {
               assert(m < i);
               assert(rowindices[m] < matrix->nrows);
               assert(newpclass < matrix->nrows);

               pclass[rowindices[m]] = newpclass;
               classsizes[newpclass]++;
            }

            if( t == k - startk )
               break;
         }

         if( k == matrix->colmatcnt[c] )
            break;
      }
   }

   SCIPfreeBufferArray(scip, &pcs);
   SCIPfreeBufferArray(scip, &rowindices);
   SCIPfreeBufferArray(scip, &values);
   SCIPfreeBufferArray(scip, &pcset);
   SCIPfreeBufferArray(scip, &classsizes);

   return SCIP_OKAY;
}

/** detect parallel rows of matrix.
 * obj coefficients are ignored.
 */
SCIP_RETCODE SCIPmatrixGetParallelCols(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real*            scale,              /**< scale factors of cols */
   int*                  pclass,             /**< parallel column classes */
   SCIP_Bool*            varineq             /**< indicating if variable is within an equation */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* values;
   int* classsizes;
   int* pcset;
   int* rowpnt;
   int* rowend;
   int* colindices;
   int* pcs;
   SCIP_Real startval;
   SCIP_Real aij;
   int startpc;
   int startk;
   int startt;
   int pcsetfill;
   int colidx;
   int k;
   int t;
   int m;
   int i;
   int r;
   int newpclass;
   int pc;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(pclass != NULL);
   assert(varineq != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &classsizes, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcset, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &values, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colindices, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcs, matrix->ncols) );

   /* init */
   BMSclearMemoryArray(scale, matrix->ncols);
   BMSclearMemoryArray(pclass, matrix->ncols);
   BMSclearMemoryArray(classsizes, matrix->ncols);
   classsizes[0] = matrix->ncols;
   pcsetfill = 0;
   for( t = 1; t < matrix->ncols; ++t )
      pcset[pcsetfill++] = t;

   /* loop over all rows */
   for( r = 0; r < matrix->nrows; ++r )
   {
      /* we consider only equations or ranged rows */
      if( !matrix->isrhsinfinite[r] )
      {
         rowpnt = matrix->rowmatind + matrix->rowmatbeg[r];
         rowend = rowpnt + matrix->rowmatcnt[r];
         valpnt = matrix->rowmatval + matrix->rowmatbeg[r];

         i = 0;
         for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
         {
            aij = *valpnt;
            colidx = *rowpnt;

            /* remember variable was part of an equation or ranged row */
            varineq[colidx] = TRUE;

            if( scale[colidx] == 0.0 )
               scale[colidx] = aij;
            assert(scale[colidx] != 0.0);

            colindices[i] = colidx;
            values[i] = aij / scale[colidx];
            pc = pclass[colidx];
            assert(pc < matrix->ncols);

            /* update class sizes and pclass set */
            assert(classsizes[pc] > 0);
            classsizes[pc]--;
            if( classsizes[pc] == 0 )
            {
               assert(pcsetfill < matrix->ncols);
               pcset[pcsetfill++] = pc;
            }
            pcs[i] = pc;

            i++;
         }

         /* sort on the pclass values */
         if( i > 1 )
         {
            SCIPsortIntIntReal(pcs, colindices, values, i);
         }

         k = 0;
         while( TRUE ) /*lint !e716*/
         {
            assert(k < i);
            startpc = pcs[k];
            startk = k;

            /* find pclass-sets */
            while( k < i && pcs[k] == startpc )
               k++;

            /* sort on the A values which have equal pclass values */
            if( k - startk > 1 )
               SCIPsortRealInt(&(values[startk]), &(colindices[startk]), k - startk);

            t = 0;
            while( TRUE ) /*lint !e716*/
            {
               assert(startk + t < i);
               startval = values[startk + t];
               startt = t;

               /* find A-sets */
               while( t < k - startk && SCIPisEQ(scip, startval, values[startk + t]) )
                  t++;

               /* get new pclass */
               newpclass = pcset[0];
               assert(pcsetfill > 0);
               pcset[0] = pcset[--pcsetfill];

               /* renumbering */
               for( m = startk + startt; m < startk + t; m++ )
               {
                  assert(m < i);
                  assert(colindices[m] < matrix->ncols);
                  assert(newpclass < matrix->ncols);

                  pclass[colindices[m]] = newpclass;
                  classsizes[newpclass]++;
               }

               if( t == k - startk )
                  break;
            }

            if( k == matrix->rowmatcnt[r] )
               break;
         }
      }
   }

   SCIPfreeBufferArray(scip, &pcs);
   SCIPfreeBufferArray(scip, &colindices);
   SCIPfreeBufferArray(scip, &values);
   SCIPfreeBufferArray(scip, &pcset);
   SCIPfreeBufferArray(scip, &classsizes);

   return SCIP_OKAY;
}


/*
 * access functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPmatrixGetColValPtr
#undef SCIPmatrixGetColIdxPtr
#undef SCIPmatrixGetColNNonzs
#undef SCIPmatrixGetNColumns
#undef SCIPmatrixGetColUb
#undef SCIPmatrixGetColLb
#undef SCIPmatrixGetColNUplocks
#undef SCIPmatrixGetColNDownlocks
#undef SCIPmatrixGetVar
#undef SCIPmatrixGetColName
#undef SCIPmatrixGetRowValPtr
#undef SCIPmatrixGetRowIdxPtr
#undef SCIPmatrixGetRowNNonzs
#undef SCIPmatrixGetRowName
#undef SCIPmatrixGetNRows
#undef SCIPmatrixGetRowLhs
#undef SCIPmatrixGetRowRhs
#undef SCIPmatrixIsRowRhsInfinity
#undef SCIPmatrixGetNNonzs
#undef SCIPmatrixGetRowMinActivity
#undef SCIPmatrixGetRowMaxActivity
#undef SCIPmatrixGetRowNMinActNegInf
#undef SCIPmatrixGetRowNMinActPosInf
#undef SCIPmatrixGetRowNMaxActNegInf
#undef SCIPmatrixGetRowNMaxActPosInf
#undef SCIPmatrixGetCons
#undef SCIPmatrixUplockConflict
#undef SCIPmatrixDownlockConflict

/** get column based start pointer of values */
SCIP_Real* SCIPmatrixGetColValPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return matrix->colmatval + matrix->colmatbeg[col];
}

/** get column based start pointer of row indices */
int* SCIPmatrixGetColIdxPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return matrix->colmatind + matrix->colmatbeg[col];
}

/** get the number of non-zero entries of this column */
int SCIPmatrixGetColNNonzs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return matrix->colmatcnt[col];
}

/** get number of columns of the matrix */
int SCIPmatrixGetNColumns(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   )
{
   assert(matrix != NULL);

   return matrix->ncols;
}

/** get upper bound of column */
SCIP_Real SCIPmatrixGetColUb(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);

   return matrix->ub[col];
}

/** get lower bound of column */
SCIP_Real SCIPmatrixGetColLb(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);

   return matrix->lb[col];
}

/** get number of uplocks of column */
int SCIPmatrixGetColNUplocks(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return matrix->nuplocks[col];
}

/** get number of downlocks of column */
int SCIPmatrixGetColNDownlocks(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return matrix->ndownlocks[col];
}

/** get variable pointer of column */
SCIP_VAR* SCIPmatrixGetVar(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return matrix->vars[col];
}

/** get name of column/variable */
const char* SCIPmatrixGetColName(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return SCIPvarGetName(matrix->vars[col]);
}

/** get row based start pointer of values */
SCIP_Real* SCIPmatrixGetRowValPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->rowmatval + matrix->rowmatbeg[row];
}

/** get row based start pointer of column indices */
int* SCIPmatrixGetRowIdxPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->rowmatind + matrix->rowmatbeg[row];
}

/** get number of non-zeros of this row */
int SCIPmatrixGetRowNNonzs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->rowmatcnt[row];
}

/** get name of row */
const char* SCIPmatrixGetRowName(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return SCIPconsGetName(matrix->cons[row]);
}

/** get number of rows of the matrix */
int SCIPmatrixGetNRows(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   )
{
   assert(matrix != NULL);

   return matrix->nrows;
}

/** get left-hand-side of row */
SCIP_Real SCIPmatrixGetRowLhs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->lhs[row];
}

/** get right-hand-side of row */
SCIP_Real SCIPmatrixGetRowRhs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->rhs[row];
}

/** flag indicating if right-hand-side of row is infinity */
SCIP_Bool SCIPmatrixIsRowRhsInfinity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->isrhsinfinite[row];
}

/** get number of non-zeros of matrix */
int SCIPmatrixGetNNonzs(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   )
{
   assert(matrix != NULL);

   return matrix->nnonzs;
}

/** get minimal activity of row */
SCIP_Real SCIPmatrixGetRowMinActivity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->minactivity[row];
}

/** get maximal activity of row */
SCIP_Real SCIPmatrixGetRowMaxActivity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->maxactivity[row];
}

/** get number of negative infinities present within minimal activity */
int SCIPmatrixGetRowNMinActNegInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->minactivityneginf[row];
}

/** get number of positive infinities present within minimal activity */
int SCIPmatrixGetRowNMinActPosInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->minactivityposinf[row];
}

/** get number of negative infinities present within maximal activity */
int SCIPmatrixGetRowNMaxActNegInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->maxactivityneginf[row];
}

/** get number of positive infinities present within maximal activity */
int SCIPmatrixGetRowNMaxActPosInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->maxactivityposinf[row];
}

/** get constraint pointer for constraint representing row */
SCIP_CONS* SCIPmatrixGetCons(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   )
{
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   return matrix->cons[row];
}

/** get if conflicting uplocks of a specific variable present */
SCIP_Bool SCIPmatrixUplockConflict(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return (SCIPvarGetNLocksUp(matrix->vars[col]) == matrix->nuplocks[col]);
}

/** get if conflicting downlocks of a specific variable present */
SCIP_Bool SCIPmatrixDownlockConflict(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   )
{
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);

   return (SCIPvarGetNLocksDown(matrix->vars[col]) == matrix->ndownlocks[col]);
}
