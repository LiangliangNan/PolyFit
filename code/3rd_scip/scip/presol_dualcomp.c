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

/**@file   presol_dualcomp.c
 * @brief  dual compensation presolver
 * @author Dieter Weninger
 *
 * This presolver looks for variables with
 *         i) objcoef >= 0 and exactly one downlock
 *        ii) objcoef <= 0 and exactly one uplock
 * and fixes the variable in case i) at the lower bound and in case ii) at the
 * upper bound if a combination of singleton continuous variables can compensate
 * the downlock in case i) and the uplock in case ii).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/pub_matrix.h"
#include "presol_dualcomp.h"

#define PRESOL_NAME            "dualcomp"
#define PRESOL_DESC            "compensate single up-/downlocks by singleton continuous variables"

/* we need singleton continuous variables for the lock compensation,
 * thus it is presumably a good idea to call this presolver before stuffing, which
 * fixes singleton continuous variables
 */
#define PRESOL_PRIORITY              -50     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_COMP_ONLY_DIS_VARS FALSE     /**< should only discrete variables be compensated? */

/*
 * Data structures
 */

/** control parameters */
struct SCIP_PresolData
{
   SCIP_Bool             componlydisvars;    /**< flag indicating if only discrete variables should be compensated */
};

/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,         /**< fix variable at lower bound */
   NOFIX   =  0,         /**< do not fix variable */
   FIXATUB =  1          /**< fix variable at upper bound */
};
typedef enum Fixingdirection FIXINGDIRECTION;

/** type of variable lock compensation */
enum Lockcompensation
{
   COMPENSATE_DOWNLOCK = 0,
   COMPENSATE_UPLOCK   = 1
};
typedef enum Lockcompensation LOCKCOMPENSATION;

/*
 * Local methods
 */

/** try to compensate a variable with a single opposite lock
    by using singleton continuous variables */
static
SCIP_RETCODE compensateVarLock(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< variable fixing candidate */
   int                   row,                /**< row index with opposite lock */
   SCIP_Real             val,                /**< value of fixing candidate in the opposite lock constraint */
   SCIP_Bool             twosides,           /**< flag indicating that two sides are present */
   LOCKCOMPENSATION      compensation,       /**< type of lock compensation */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   int*                  nfixings            /**< number of possible fixings */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   SCIP_VAR* var;
   int colidx;
   SCIP_Real coef;
   SCIP_Real lhs;
   SCIP_Real delta;
   SCIP_Bool trytofix;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool deltaisinf;
   SCIP_Real ratio;
   SCIP_Bool multrowbyminusone;
   SCIP_Bool singleton;
   SCIP_Real offset;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= col && col < SCIPmatrixGetNColumns(matrix));
   assert(0 <= row && row < SCIPmatrixGetNRows(matrix));
   assert(compensation == COMPENSATE_DOWNLOCK || compensation == COMPENSATE_UPLOCK);
   assert(varstofix != NULL);
   assert(nfixings != NULL);

   /* the variable for compensation should not be a compensation variable itself */
   assert(!(SCIPmatrixGetColNNonzs(matrix,col) == 1 && SCIPvarGetType(SCIPmatrixGetVar(matrix,col)) == SCIP_VARTYPE_CONTINUOUS));

   /* try lock compensation only if minimum one singleton continuous variable is present */
   singleton = FALSE;
   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   for( ; rowpnt < rowend; rowpnt++ )
   {
      var = SCIPmatrixGetVar(matrix, *rowpnt);

      if( SCIPmatrixGetColNNonzs(matrix, *rowpnt) == 1 &&
         SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS &&
         SCIPvarGetNLocksUp(var) == SCIPmatrixGetColNUplocks(matrix, *rowpnt) &&
         SCIPvarGetNLocksDown(var) == SCIPmatrixGetColNDownlocks(matrix, *rowpnt)
         )
      {
         /* minimal one valid compensation variable is present in this row */
         singleton = TRUE;
         break;
      }
   }

   /* return if no compensation variable is available */
   if( !singleton )
      return SCIP_OKAY;


   /* we perform the following transformations afterwards:
    *
    * lhs <= a1 x1 + a2 x2 + ... an xn <= rhs
    * with a1, a2, ..., an >= 0.
    *
    * for the downlock case we multiply the constraint in thought by (-1)
    * if the corresponding coefficient is negative.
    *
    * we attribute the uplock case to the downlock case by multiplying
    * in thought the corresponding column by (-1).
    */
   multrowbyminusone = FALSE;
   if( compensation == COMPENSATE_DOWNLOCK )
   {
      if( SCIPisLT(scip,val,0.0) )
         multrowbyminusone = TRUE;
   }
   else
   {
      assert(compensation == COMPENSATE_UPLOCK);

      /* in the uplock case we multiply the column in thought by (-1) and
       * thus we need to multiply the constraint by (-1) to get a positive coefficient
       */
      if( SCIPisGT(scip,val,0.0) )
         multrowbyminusone = TRUE;
   }

   /* we need the objective coefficient and constraint coefficient ratio
    * to later preserve optimality.
    * further we need to consider multiplications of the constraint by (-1).
    * for ranged rows and equalities we switch to the rhs.
    */
   lhs = SCIPmatrixGetRowLhs(matrix, row);
   ratio = SCIPvarGetObj( SCIPmatrixGetVar(matrix,col) ) / val;
   if( multrowbyminusone )
   {
      if( twosides )
         lhs = -SCIPmatrixGetRowRhs(matrix, row);
      else
         lhs = -lhs;

      ratio = -ratio;
   }

   offset = 0.0;
   trytofix = TRUE;
   delta = 0;
   deltaisinf = FALSE;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for( ; rowpnt < rowend; rowpnt++, valpnt++ )
   {
      colidx = *rowpnt;
      coef = *valpnt;
      var = SCIPmatrixGetVar(matrix, colidx);
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      if( colidx == col )
      {
         /* this is the variable which we want to compensate */

         if( compensation == COMPENSATE_DOWNLOCK )
         {
            if( SCIPisInfinity(scip, -lb) )
            {
               trytofix = FALSE;
               break;
            }
            else
            {
               if( multrowbyminusone )
                  offset += (-coef) * lb;
               else
                  offset += coef * lb;
            }
         }
         else
         {
            if( SCIPisInfinity(scip, ub) )
            {
               trytofix = FALSE;
               break;
            }
            else
            {
               /* for the uplock case we have opposed sign for the coefficient as
                * in the downlock case.
                * the multiplication of the column results in swapping the negative bounds.
                */
               if( multrowbyminusone )
                  offset += coef * (-ub);
               else
                  offset += (-coef) * (-ub);
            }
         }
      }
      else if( SCIPmatrixGetColNNonzs(matrix, colidx) == 1 &&
         SCIPvarGetNLocksUp(var) == SCIPmatrixGetColNUplocks(matrix, colidx) &&
         SCIPvarGetNLocksDown(var) == SCIPmatrixGetColNDownlocks(matrix, colidx) &&
         SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         /* this is singleton continuous variable and
          * thus a valid compensation candidate
          */

         if( SCIPisLT(scip,coef,0.0) )
         {
            /* coef < 0 */

            if( multrowbyminusone )
            {
               if( SCIPisInfinity(scip, -lb) )
               {
                  trytofix = FALSE;
                  break;
               }

               /* we have a negative coefficient and the row is multiplied by (-1)
                * thus actually we have a positive coefficient
                */
               offset += (-coef) * lb;

               /* only consider singleton continuous variables with a better or the same
                * obj/coef ratio for preserving optimality
                */
               if( SCIPisLE(scip,SCIPvarGetObj(SCIPmatrixGetVar(matrix, colidx))/(-coef), ratio) )
               {
                  if( SCIPisInfinity(scip, ub) )
                  {
                     deltaisinf = TRUE;
                     break;
                  }

                  /* calculate the contribution to the compensation value */
                  delta += (-coef) * (ub - lb);
               }
            }
            else
            {
               if( SCIPisInfinity(scip, ub) )
               {
                  trytofix = FALSE;
                  break;
               }

               /* we have a negative coefficient and hence need to multiply the column by (-1).
                * this means the bounds swap and change the sign
                */
               offset += (-coef) * (-ub);

               /* only consider singleton continuous variables with a better or the same
                * obj/coef ratio for preserving optimality
                */
               if( SCIPisLE(scip,SCIPvarGetObj(SCIPmatrixGetVar(matrix, colidx))/coef, ratio) )
               {
                  if( SCIPisInfinity(scip, -lb) )
                  {
                     deltaisinf = TRUE;
                     break;
                  }

                  /* calculate the contribution to the compensation value */
                  delta += (-coef) * (ub - lb);
               }
            }
         }
         else
         {
            /* coef >= 0 */

            if( multrowbyminusone )
            {
               /* we have a positive or zero coefficient and the row is multiplied by (-1) */
               if( SCIPisInfinity(scip, ub) )
               {
                  trytofix = FALSE;
                  break;
               }

               /* we have a positive or zero coefficient and multiply in thought the constraint
                * by (-1) thus we have actually a negative coefficient and multiply the column by (-1).
                * therefore the sign of the coefficient does not change but the bounds swap and change
                * the sign.
                */
               offset += coef * (-ub);

               /* we have a positive or zero coefficient and multiply in thought the constraint
                * by (-1) which delivers the ratio.
                * a further multiplication of the column does not change anything.
                */
               if( SCIPisLE(scip,SCIPvarGetObj(SCIPmatrixGetVar(matrix, colidx))/(-coef), ratio) )
               {
                  if( SCIPisInfinity(scip, -lb) )
                  {
                     deltaisinf = TRUE;
                     break;
                  }

                  /* calculate the contribution to the compensation value */
                  delta += coef * (ub - lb);
               }
            }
            else
            {
               if( SCIPisInfinity(scip, -lb) )
               {
                  trytofix = FALSE;
                  break;
               }

               /* we have positive coefficient and do not need to multiply anything by (-1) */
               offset += coef * lb;

               if( SCIPisLE(scip,SCIPvarGetObj(SCIPmatrixGetVar(matrix, colidx))/coef, ratio) )
               {
                  if( SCIPisInfinity(scip, ub) )
                  {
                     deltaisinf = TRUE;
                     break;
                  }

                  /* calculate the contribution to the compensation value */
                  delta += coef * (ub - lb);
               }
            }
         }
      }
      else
      {
         /* remaining variables */

         /* the reasons for the following signs are the same as for the singleton
          * continuous variables
          */
         if( SCIPisLT(scip,coef,0.0) )
         {
            if( multrowbyminusone )
            {
               if( SCIPisInfinity(scip, -lb) )
               {
                  trytofix = FALSE;
                  break;
               }

               offset += (-coef) * lb;
            }
            else
            {
               if( SCIPisInfinity(scip, ub) )
               {
                  trytofix = FALSE;
                  break;
               }

               offset += (-coef) * (-ub);
            }
         }
         else
         {
            if( multrowbyminusone )
            {
               if( SCIPisInfinity(scip, ub) )
               {
                  trytofix = FALSE;
                  break;
               }

               offset += coef * (-ub);
            }
            else
            {
               if( SCIPisInfinity(scip, -lb) )
               {
                  trytofix = FALSE;
                  break;
               }

               offset += coef * lb;
            }
         }
      }
   }

   /* avoid fixings to infinite values or fixings of already fixed variables */
   if( trytofix && varstofix[col] == NOFIX)
   {
      /* feasibility is secured if the compensation value delta
       * is large enough to compensate the value lhs-offset
       */
      if( deltaisinf || SCIPisLE(scip, lhs-offset, delta) )
      {
         if( compensation == COMPENSATE_UPLOCK )
         {
            if( !SCIPisInfinity(scip,SCIPvarGetUbGlobal(SCIPmatrixGetVar(matrix, col))) )
            {
               varstofix[col] = FIXATUB;
               (*nfixings)++;

#ifdef SCIP_MORE_DEBUG
               SCIPmatrixPrintRow(scip, matrix, row);
               SCIPdebugMsg(scip, "%s, bds=[%.2f,%.2f], obj=%.2f, nnonzs=%d, type=%s, fix=ub, %.1f <= %.1f\n",
                  SCIPvarGetName(SCIPmatrixGetVar(matrix, col)),SCIPvarGetLbGlobal(SCIPmatrixGetVar(matrix, col)),
                  SCIPvarGetUbGlobal(SCIPmatrixGetVar(matrix, col)), SCIPvarGetObj(SCIPmatrixGetVar(matrix, col)),
                  SCIPmatrixGetColNNonzs(matrix, col),
                  SCIPvarGetType(SCIPmatrixGetVar(matrix, col))==SCIP_VARTYPE_CONTINUOUS ? "con" : "dis",
                  lhs-offset, delta);
#endif
            }
         }
         else
         {
            if( !SCIPisInfinity(scip,-SCIPvarGetLbGlobal(SCIPmatrixGetVar(matrix, col))) )
            {
               varstofix[col] = FIXATLB;
               (*nfixings)++;

#ifdef SCIP_MORE_DEBUG
               SCIPmatrixPrintRow(scip, matrix, row);
               SCIPdebugMsg(scip, "%s, bds=[%.2f,%.2f], obj=%.2f, nnonzs=%d, type=%s, fix=lb, %.1f <= %.1f\n",
                  SCIPvarGetName(SCIPmatrixGetVar(matrix, col)),SCIPvarGetLbGlobal(SCIPmatrixGetVar(matrix, col)),
                  SCIPvarGetUbGlobal(SCIPmatrixGetVar(matrix, col)), SCIPvarGetObj(SCIPmatrixGetVar(matrix, col)),
                  SCIPmatrixGetColNNonzs(matrix, col),
                  SCIPvarGetType(SCIPmatrixGetVar(matrix, col))==SCIP_VARTYPE_CONTINUOUS ? "con" : "dis",
                  lhs-offset, delta);
#endif
            }
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualcomp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDualcomp(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualcomp)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* don't run if no compensation variables are present */
   if( SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

   if( !SCIPallowDualReds(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   /* we only work on pure MIPs currently */
   if( initialized && complete )
   {
      int ncols;
      int i;
      SCIP_Real* valpnt;
      int* colpnt;
      int* colend;
      int row;
      SCIP_VAR* var;
      SCIP_Bool inspect;
      SCIP_Real val;
      FIXINGDIRECTION* varstofix;
      int nfixings;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Bool twosides;

      ncols = SCIPmatrixGetNColumns(matrix);
      nfixings = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, ncols) );
      BMSclearMemoryArray(varstofix, ncols);

      for(i = 0; i < ncols; i++)
      {
         var = SCIPmatrixGetVar(matrix, i);

         /* exclude compensation variables itself for compensation */
         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS &&
            SCIPmatrixGetColNNonzs(matrix, i) == 1 )
            continue;

         /* if requested exclude continuous variables for compensation */
         if( presoldata->componlydisvars && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            continue;

         /* verifiy that this variable has one uplock and that the uplocks are consistent */
         if( SCIPvarGetNLocksUp(var) == 1 &&
            SCIPmatrixGetColNUplocks(matrix, i) == 1 &&
            SCIPisLE(scip, SCIPvarGetObj(var), 0.0) )
         {
            row = -1;
            val = 0.0;
            inspect = FALSE;
            twosides = FALSE;
            colpnt = SCIPmatrixGetColIdxPtr(matrix, i);
            colend = colpnt + SCIPmatrixGetColNNonzs(matrix, i);
            valpnt = SCIPmatrixGetColValPtr(matrix, i);

            /* search row which causes the uplock */
            for( ; (colpnt < colend); colpnt++, valpnt++ )
            {
               row = *colpnt;
               val = *valpnt;
               lhs = SCIPmatrixGetRowLhs(matrix, row);
               rhs = SCIPmatrixGetRowRhs(matrix, row);

               if( SCIPisEQ(scip, lhs, rhs) )
               {
                  /* equation */
                  inspect = TRUE;
                  twosides = TRUE;
                  break;
               }
               else if( SCIPmatrixIsRowRhsInfinity(matrix, row) )
               {
                  /* >= */
                  if( SCIPisLT(scip, val, 0.0) )
                  {
                     inspect = TRUE;
                     break;
                  }
               }
               else if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) )
               {
                  /* ranged row */
                  inspect = TRUE;
                  twosides = TRUE;
                  break;
               }
            }

            assert(inspect);

            if( inspect ) /*lint !e774*/
            {
               assert(row >= 0);
               assert(!SCIPisZero(scip, val));

               /* try to fix variable i at the upper bound */
               SCIP_CALL( compensateVarLock(scip, matrix, i, row, val,
                     twosides, COMPENSATE_UPLOCK, varstofix, &nfixings) );
            }
         }
         /* verifiy that this variable has one downlock and that the downlocks are consistent */
         else if( SCIPvarGetNLocksDown(var) == 1 &&
            SCIPmatrixGetColNDownlocks(matrix, i) == 1 &&
            SCIPisGE(scip, SCIPvarGetObj(var), 0.0) )
         {
            row = -1;
            val = 0.0;
            inspect = FALSE;
            twosides = FALSE;
            colpnt = SCIPmatrixGetColIdxPtr(matrix, i);
            colend = colpnt + SCIPmatrixGetColNNonzs(matrix, i);
            valpnt = SCIPmatrixGetColValPtr(matrix, i);

            /* search row which causes the downlock */
            for( ; (colpnt < colend); colpnt++, valpnt++ )
            {
               row = *colpnt;
               val = *valpnt;
               lhs = SCIPmatrixGetRowLhs(matrix, row);
               rhs = SCIPmatrixGetRowRhs(matrix, row);

               if( SCIPisEQ(scip, lhs, rhs) )
               {
                  /* equation */
                  inspect = TRUE;
                  twosides = TRUE;
                  break;
               }
               else if( SCIPmatrixIsRowRhsInfinity(matrix, row) )
               {
                  /* >= */
                  if( SCIPisGT(scip, val, 0.0) )
                  {
                     inspect = TRUE;
                     break;
                  }
               }
               else if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) )
               {
                  /* ranged row */
                  inspect = TRUE;
                  twosides = TRUE;
                  break;
               }
            }

            assert(inspect);

            if( inspect ) /*lint !e774*/
            {
               assert(row >= 0);
               assert(!SCIPisZero(scip, val));

               /* try to fix variable i at the lower bound */
               SCIP_CALL( compensateVarLock(scip, matrix, i, row, val,
                     twosides, COMPENSATE_DOWNLOCK, varstofix, &nfixings) );
            }
         }
      }

      if( nfixings > 0 )
      {
         int v;
         int oldnfixedvars;
         int numupperboundfixings;
         int numlowerboundfixings;
         int numcontinuousfixings;
         int numdiscretefixings;

         oldnfixedvars = *nfixedvars;
         numupperboundfixings = 0;
         numlowerboundfixings = 0;
         numcontinuousfixings = 0;
         numdiscretefixings = 0;

         /* look for fixable variables */
         for( v = ncols - 1; v >= 0; --v )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            var = SCIPmatrixGetVar(matrix, v);

            if( varstofix[v] == FIXATLB )
            {
               SCIP_Real lb;

               lb = SCIPvarGetLbGlobal(var);

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
               numlowerboundfixings++;

               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  numcontinuousfixings++;
               else
                  numdiscretefixings++;
            }
            else if( varstofix[v] == FIXATUB )
            {
               SCIP_Real ub;

               ub = SCIPvarGetUbGlobal(var);

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
               numupperboundfixings++;

               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  numcontinuousfixings++;
               else
                  numdiscretefixings++;
            }
         }

         if( *result != SCIP_CUTOFF && *nfixedvars > oldnfixedvars )
            *result = SCIP_SUCCESS;

         SCIPdebugMsg(scip, "### lbfixes: %d, ubfixes: %d, con: %d, dis: %d\n",
            numlowerboundfixings, numupperboundfixings,
            numcontinuousfixings, numdiscretefixings);
      }

      SCIPfreeBufferArray(scip, &varstofix);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeDualcomp)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** creates the dualcomp presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualcomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create dualcomp presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDualcomp, presoldata) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDualcomp) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeDualcomp) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/dualcomp/componlydisvars",
         "should only discrete variables be compensated?",
         &presoldata->componlydisvars, FALSE, DEFAULT_COMP_ONLY_DIS_VARS, NULL, NULL) );

   return SCIP_OKAY;
}
