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

/**@file   presol_stuffing.c
 * @brief  fix singleton continuous variables
 * @author Dieter Weninger
 *
 * Investigate singleton continuous variables if one can be fixed at a bound.
 *
 * @todo enhancement from singleton continuous variables to continuous
 *       variables with only one lock in a common row
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>


#include "scip/pub_matrix.h"
#include "presol_stuffing.h"

#define PRESOL_NAME            "stuffing"
#define PRESOL_DESC            "fix redundant singleton continuous variables"
#define PRESOL_PRIORITY             -100     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               0     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,         /**< fix variable at lower bound */
   NOFIX   =  0,         /**< do not fix variable */
   FIXATUB =  1          /**< fix variable at upper bound */
};
typedef enum Fixingdirection FIXINGDIRECTION;

/*
 * Local methods
 */

/** try to fix singleton continuous variables */
static
SCIP_RETCODE singletonColumnStuffing(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   int*                  nfixings            /**< number of possible fixings */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* colratios;
   SCIP_Real* colcoeffs;
   SCIP_Bool* rowprocessed;
   int* rowpnt;
   int* rowend;
   int* colindices;
   int* dummy;
   SCIP_Bool* swapped;
   SCIP_Real upperconst;
   SCIP_Real lowerconst;
   SCIP_Real rhs;
   SCIP_Bool tryfixing;
   int idx;
   int col;
   int row;
   int fillcnt;
   int k;
   int nrows;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(varstofix != NULL);
   assert(nfixings != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &colindices, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colratios, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colcoeffs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dummy, ncols) );

   SCIP_CALL( SCIPallocBufferArray(scip, &swapped, ncols) );
   BMSclearMemoryArray(swapped, ncols);

   SCIP_CALL( SCIPallocBufferArray(scip, &rowprocessed, nrows) );
   BMSclearMemoryArray(rowprocessed, nrows);

   for( col = 0; col < ncols; col++ )
   {
      /* consider only rows with minimal one continuous singleton column */
      if( SCIPmatrixGetColNNonzs(matrix, col) == 1
         && SCIPvarGetType(SCIPmatrixGetVar(matrix, col)) == SCIP_VARTYPE_CONTINUOUS
         && SCIPvarGetNLocksUp(SCIPmatrixGetVar(matrix, col)) == SCIPmatrixGetColNUplocks(matrix, col)
         && SCIPvarGetNLocksDown(SCIPmatrixGetVar(matrix, col)) == SCIPmatrixGetColNDownlocks(matrix, col) )
      {
         row = *(SCIPmatrixGetColIdxPtr(matrix, col));
         if( rowprocessed[row] )
            continue;

         rowprocessed[row] = TRUE;

         /* treat all >= rows from matrix, but internally we transform to <= relation */
         if( SCIPmatrixIsRowRhsInfinity(matrix, row) )
         {
            fillcnt = 0;
            tryfixing = TRUE;
            upperconst = 0.0;
            lowerconst = 0.0;
            rhs = -SCIPmatrixGetRowLhs(matrix, row);

            rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
            rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
            valpnt = SCIPmatrixGetRowValPtr(matrix, row);

            for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
            {
               SCIP_Real coef;
               SCIP_VAR* var;
               SCIP_Real lb;
               SCIP_Real ub;

               coef = -(*valpnt);
               idx = *rowpnt;
               var = SCIPmatrixGetVar(matrix, idx);
               lb = SCIPvarGetLbGlobal(var);
               ub = SCIPvarGetUbGlobal(var);

               /* we need to check if this is a singleton continuous variable and
                * all constraints containing this variable are present inside
                * the mixed integer linear matrix
                */
               if( SCIPmatrixGetColNNonzs(matrix, idx) == 1 &&
                   (SCIPvarGetNLocksUp(var) + SCIPvarGetNLocksDown(var)) == 1 &&
                   SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  if( SCIPisLT(scip, SCIPvarGetObj(var), 0.0) && SCIPisGT(scip, coef, 0.0) )
                  {
                     /* case 1: obj < 0 and coef > 0 */
                     if( SCIPisInfinity(scip, -lb) )
                     {
                        tryfixing = FALSE;
                        break;
                     }

                     upperconst += coef * lb;
                     lowerconst += coef * lb;
                     colratios[fillcnt] = SCIPvarGetObj(var) / coef;
                     colindices[fillcnt] = idx;
                     colcoeffs[fillcnt] = coef;
                     fillcnt++;
                  }
                  else if( SCIPisGT(scip, SCIPvarGetObj(var), 0.0) && SCIPisLT(scip, coef, 0.0) )
                  {
                     /* case 2: obj > 0 and coef < 0 */
                     if( SCIPisInfinity(scip, ub) )
                     {
                        tryfixing = FALSE;
                        break;
                     }

                     /* multiply column by (-1) to become case 1.
                      * now bounds are swapped: ub := -lb, lb := -ub
                      */
                     swapped[idx] = TRUE;
                     upperconst += coef * ub;
                     lowerconst += coef * ub;
                     colratios[fillcnt] = SCIPvarGetObj(var) / coef;
                     colindices[fillcnt] = idx;
                     colcoeffs[fillcnt] = -coef;
                     fillcnt++;
                  }
                  else if( SCIPisGE(scip, SCIPvarGetObj(var), 0.0) && SCIPisGE(scip, coef, 0.0) )
                  {
                     /* case 3: obj >= 0 and coef >= 0 is handled by duality fixing.
                      *  we only consider the lower bound for the constants
                      */
                     if( SCIPisInfinity(scip, -lb) )
                     {
                        /* maybe unbounded */
                        tryfixing = FALSE;
                        break;
                     }

                     upperconst += coef * lb;
                     lowerconst += coef * lb;
                  }
                  else
                  {
                     /* case 4: obj <= 0 and coef <= 0 is also handled by duality fixing.
                      * we only consider the upper bound for the constants
                      */
                     assert(SCIPisLE(scip, SCIPvarGetObj(var), 0.0));
                     assert(SCIPisLE(scip, coef, 0.0));

                     if( SCIPisInfinity(scip, ub) )
                     {
                        /* maybe unbounded */
                        tryfixing = FALSE;
                        break;
                     }

                     upperconst += coef * ub;
                     lowerconst += coef * ub;
                  }
               }
               else
               {
                  /* consider contribution of discrete variables, non-singleton
                   * continuous variables and variables with more than one lock
                   */
                  if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
                  {
                     tryfixing = FALSE;
                     break;
                  }

                  if( coef > 0 )
                  {
                     upperconst += coef * ub;
                     lowerconst += coef * lb;
                  }
                  else
                  {
                     upperconst += coef * lb;
                     lowerconst += coef * ub;
                  }
               }
            }

            if( tryfixing )
            {
               SCIPsortRealRealIntInt(colratios, colcoeffs, colindices, dummy, fillcnt);

               /* verify which singleton continuous variable can be fixed */
               for( k = 0; k < fillcnt; k++ )
               {
                  SCIP_VAR* var;
                  SCIP_Real lb;
                  SCIP_Real ub;
                  SCIP_Real delta;

                  idx = colindices[k];
                  var = SCIPmatrixGetVar(matrix, idx);
                  lb = SCIPvarGetLbGlobal(var);
                  ub = SCIPvarGetUbGlobal(var);

                  /* stop fixing if variable bounds are not finite */
                  if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
                     break;

                  assert(SCIPmatrixGetColNNonzs(matrix, idx) == 1);
                  assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
                  assert((SCIPvarGetNLocksUp(var) + SCIPvarGetNLocksDown(var)) == 1);
                  assert(colcoeffs[k] >= 0);

                  /* calculate the change in the row activities if this variable changes
                   * its value from to its other bound
                   */
                  if( swapped[idx] )
                     delta = -(lb - ub) * colcoeffs[k];
                  else
                     delta =  (ub - lb) * colcoeffs[k];

                  assert(delta >= 0);

                  if( SCIPisLE(scip, delta, rhs - upperconst) )
                  {
                     if( swapped[idx] )
                        varstofix[idx] = FIXATLB;
                     else
                        varstofix[idx] = FIXATUB;

                     (*nfixings)++;
                  }
                  else if( SCIPisLE(scip, rhs, lowerconst) )
                  {
                     if( swapped[idx] )
                        varstofix[idx] = FIXATUB;
                     else
                        varstofix[idx] = FIXATLB;

                     (*nfixings)++;
                  }

                  upperconst += delta;
                  lowerconst += delta;
               }
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &rowprocessed);
   SCIPfreeBufferArray(scip, &swapped);
   SCIPfreeBufferArray(scip, &dummy);
   SCIPfreeBufferArray(scip, &colcoeffs);
   SCIPfreeBufferArray(scip, &colratios);
   SCIPfreeBufferArray(scip, &colindices);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyStuffing)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolStuffing(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecStuffing)
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

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized )
   {
      FIXINGDIRECTION* varstofix;
      int nfixings;
      int ncols;

      nfixings = 0;
      ncols = SCIPmatrixGetNColumns(matrix);

      SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, ncols) );
      BMSclearMemoryArray(varstofix, ncols);

      SCIP_CALL( singletonColumnStuffing(scip, matrix, varstofix, &nfixings) );

      if( nfixings > 0 )
      {
         int v;
         int oldnfixedvars;

         oldnfixedvars = *nfixedvars;

         /* look for fixable variables */
         for( v = ncols - 1; v >= 0; --v )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIP_VAR* var;

            var = SCIPmatrixGetVar(matrix, v);

            if( varstofix[v] == FIXATLB )
            {
               SCIP_Real lb;

               assert(SCIPvarGetNLocksUp(var) == SCIPmatrixGetColNUplocks(matrix, v) &&
                  SCIPvarGetNLocksDown(var) == SCIPmatrixGetColNDownlocks(matrix, v));

               lb = SCIPvarGetLbGlobal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

               /* avoid fixings to infinite values */
               assert(!SCIPisInfinity(scip, -lb));

               SCIPdebugMsg(scip, "Fix variable %s at lower bound %.15g\n", SCIPvarGetName(var), lb);

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
            }
            else if( varstofix[v] == FIXATUB )
            {
               SCIP_Real ub;

               assert(SCIPvarGetNLocksUp(var) == SCIPmatrixGetColNUplocks(matrix, v) &&
                  SCIPvarGetNLocksDown(var) == SCIPmatrixGetColNDownlocks(matrix, v));

               ub = SCIPvarGetUbGlobal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

               /* avoid fixings to infinite values */
               assert(!SCIPisInfinity(scip, ub));

               SCIPdebugMsg(scip, "Fix variable %s at upper bound %.15g\n", SCIPvarGetName(var), ub);

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
            }
         }

         if( *result != SCIP_CUTOFF && *nfixedvars > oldnfixedvars )
            *result = SCIP_SUCCESS;
      }

      SCIPfreeBufferArray(scip, &varstofix);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the stuffing presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolStuffing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecStuffing, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyStuffing) );

   return SCIP_OKAY;
}

/*lint --e{749}*/
