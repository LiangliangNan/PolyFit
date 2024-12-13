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

/**@file   nlhdlr_default.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  default nonlinear handler that calls expression handler methods
 * @author Stefan Vigerske
 */

#include <string.h>

#include "scip/nlhdlr_default.h"
#include "scip/pub_nlhdlr.h"
#include "scip/cons_nonlinear.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME            "default"
#define NLHDLR_DESC            "default handler for expressions"
#define NLHDLR_DETECTPRIORITY  0
#define NLHDLR_ENFOPRIORITY    0

/** translate from one value of infinity to another
 *
 *  if val is &ge; infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

#define UNDERESTIMATEUSESACTIVITY 0x1u  /**< whether underestimation uses activity */
#define OVERESTIMATEUSESACTIVITY  0x2u  /**< whether overestimation uses activity */

/*lint -e666*/
/*lint -e850*/

/** evaluates an expression w.r.t. the values in the auxiliary variables */
static
SCIP_RETCODE evalExprInAux(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            val,                /**< buffer to store value of expression */
   SCIP_SOL*             sol                 /**< solution to be evaluated */
   )
{
   SCIP_Real* childvals;
   SCIP_VAR* childvar;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(val != NULL);
   assert(SCIPexprGetNChildren(expr) > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &childvals, SCIPexprGetNChildren(expr)) );

   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      childvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[c]);
      /* there should be an auxiliary variable, because we created them in detect for every child if we said that we will separate;
       * at the moment, EVALAUX should only be called for nlhdlrs that said they will separate
       * if that changes, then we should handle this here, e.g., via *val = SCIPexprGetEvalValue(expr); break;
       */
      assert(childvar != NULL);

      childvals[c] = SCIPgetSolVal(scip, sol, childvar);
   }

   SCIP_CALL( SCIPcallExprEval(scip, expr, childvals, val) );

   SCIPfreeBufferArray(scip, &childvals);

   return SCIP_OKAY;
}

/** check whether expression should be handled by the default nlhdlr
 *
 * if no nlhdlr so far provides enforcement or boundtightening for expr, then the default nlhdlr takes over
 */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectDefault)
{ /*lint --e{715}*/
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_Bool estimatebelowusesactivity = FALSE;
   SCIP_Bool estimateaboveusesactivity = FALSE;
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);
   assert(nlhdlrexprdata != NULL);

   exprhdlr = SCIPexprGetHdlr(expr);
   assert(exprhdlr != NULL);

   if( (*enforcing & SCIP_NLHDLR_METHOD_ACTIVITY) == 0 )
   {
      /* expr handlers having reverseprop but no inteval is something that we don't support at the moment for simplicity */
      assert(!SCIPexprhdlrHasReverseProp(exprhdlr) || SCIPexprhdlrHasIntEval(exprhdlr));

      /* participate in inteval and/or reverseprop if that is not yet provided in enforcing and we have inteval */
      if( SCIPexprhdlrHasIntEval(exprhdlr) )
         *participating = SCIP_NLHDLR_METHOD_ACTIVITY;
   }

   /* participate in sepa if exprhdlr for expr has an estimate callback and sepa below or above is still missing */
   if( ((*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) != SCIP_NLHDLR_METHOD_SEPABOTH) && SCIPexprhdlrHasEstimate(exprhdlr) )
   {
      /* communicate back that the nlhdlr will provide the separation on the currently missing sides */
      if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABELOW) == 0 )
         *participating |= SCIP_NLHDLR_METHOD_SEPABELOW;

      if( (*enforcing & SCIP_NLHDLR_METHOD_SEPAABOVE) == 0 )
         *participating |= SCIP_NLHDLR_METHOD_SEPAABOVE;
   }

   if( !*participating )
      return SCIP_OKAY;

   /* since this is the default handler, we enforce where we participate */
   *enforcing |= *participating;

   /* increment activity usage counter and create auxiliary variables if necessary
    * if separating, first guess whether we will use activities in estimate (distinguish under- and overestimation)
    * we assume that the exprhdlr will use activity on all children iff we are estimating on a nonconvex side
    * TODO it would be better to request this information directly from the exprhdlr than inferring it from curvature,
    * but with the currently available exprhdlr that wouldn't make a difference
    */
   if( *participating & SCIP_NLHDLR_METHOD_SEPABOTH )
   {
      SCIP_EXPRCURV* childcurv;

      /* allocate memory to store the required curvature of the children (though we don't use it) */
      SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, SCIPexprGetNChildren(expr)) );

      if( *participating & SCIP_NLHDLR_METHOD_SEPABELOW )
      {
         /* check whether the expression is convex */
         SCIP_Bool isconvex;
         SCIP_CALL( SCIPcallExprCurvature(scip, expr, SCIP_EXPRCURV_CONVEX, &isconvex, childcurv) );
         estimatebelowusesactivity = !isconvex;
      }

      if( *participating & SCIP_NLHDLR_METHOD_SEPAABOVE )
      {
         /* check whether the expression is concave */
         SCIP_Bool isconcave;
         SCIP_CALL( SCIPcallExprCurvature(scip, expr, SCIP_EXPRCURV_CONCAVE, &isconcave, childcurv) );
         estimateaboveusesactivity = !isconcave;
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &childcurv);
   }

   /* indicate enforcement methods required in children:
    * - if separating, make sure that (auxiliary) variable will exist
    * - if activity computation, then register activity usage
    * - if estimating on a non-convex side, then indicate activity usage for separation for that side
    */
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      /* todo skip auxvarusage for value-expressions? would then need update in evalExprInAux, too */
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(expr)[c],
         *participating & SCIP_NLHDLR_METHOD_SEPABOTH,
         *participating & SCIP_NLHDLR_METHOD_ACTIVITY, estimatebelowusesactivity, estimateaboveusesactivity) );
   }

   /* remember estimatebelowusesactivity and estimateaboveusesactivity in nlhdlrexprdata */
   *nlhdlrexprdata = (SCIP_NLHDLREXPRDATA*)(size_t)((estimatebelowusesactivity ? UNDERESTIMATEUSESACTIVITY : 0x0u)
      | (estimateaboveusesactivity ? OVERESTIMATEUSESACTIVITY : 0x0u));

   return SCIP_OKAY;
}

/** evaluate expression w.r.t. values of auxiliary variables in children */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalAuxDefault)
{ /*lint --e{715}*/
   assert(expr != NULL);
   assert(auxvalue != NULL);

   SCIP_CALL( evalExprInAux(scip, expr, auxvalue, sol) );

   return SCIP_OKAY;
}

/** initialize LP relaxation by initial estimators */
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaDefault)
{ /*lint --e{715}*/
   SCIP_INTERVAL* childrenbounds;
   SCIP_Real* coefs[SCIP_EXPR_MAXINITESTIMATES];
   SCIP_Real constant[SCIP_EXPR_MAXINITESTIMATES];
   SCIP_VAR* auxvar;
   SCIP_ROWPREP* rowprep;
   int nreturned;
   int i, j;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   if( !SCIPexprhdlrHasInitEstimates(SCIPexprGetHdlr(expr)) )
      return SCIP_OKAY;

   SCIPdebug( SCIPinfoMessage(scip, NULL, "initsepa exprhdlr %s for expr ", SCIPexprhdlrGetName(SCIPexprGetHdlr(expr))) );
   SCIPdebug( SCIPprintExpr(scip, expr, NULL) );
   SCIPdebug( SCIPinfoMessage(scip, NULL, "\n") );

   /* use global bounds of auxvar as global valid bounds for children
    * if at root node (thus local=global) and estimate actually uses bounds, then intersect with (local) activity of expression
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &childrenbounds, SCIPexprGetNChildren(expr)) );
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      auxvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[i]);
      assert(auxvar != NULL);

      SCIPintervalSetBounds(&childrenbounds[i],
         -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -SCIPvarGetLbGlobal(auxvar)),
          infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  SCIPvarGetUbGlobal(auxvar)));

      if( SCIPgetDepth(scip) == 0 &&
          ((underestimate && ((size_t)nlhdlrexprdata & UNDERESTIMATEUSESACTIVITY)) ||
           (overestimate  && ((size_t)nlhdlrexprdata & OVERESTIMATEUSESACTIVITY ))) )
      {
         SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(expr)[i]) );
         SCIPintervalIntersect(&childrenbounds[i], childrenbounds[i], SCIPexprGetActivity(SCIPexprGetChildren(expr)[i]));
      }

      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childrenbounds[i]) )
      {
         SCIPdebugMsg(scip, "activity for expression %d (unexpectedly) empty in initsepa\n", i);
         *infeasible = TRUE;
         SCIPfreeBufferArray(scip, &childrenbounds);
         return SCIP_OKAY;
      }
   }

   /* allocate each coefficients array */
   for( i = 0; i < SCIP_EXPR_MAXINITESTIMATES; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &coefs[i], SCIPexprGetNChildren(expr)) );
   }

   /* create rowprep */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, SCIPexprGetNChildren(expr)+1) );

   /* call the separation initialization callback of the expression handler and turn estimates into SCIP rows */
   for( i = 0; i < 2 && !*infeasible; ++i )
   {
      nreturned = 0;
      if( i == 0 && underestimate )
      {
         SCIP_CALL( SCIPcallExprInitestimates(scip, expr, childrenbounds, FALSE, coefs, constant, &nreturned) );
         assert(SCIProwprepGetSidetype(rowprep) == SCIP_SIDETYPE_RIGHT);
      }
      if( i == 1 && overestimate )
      {
         SCIP_CALL( SCIPcallExprInitestimates(scip, expr, childrenbounds, TRUE, coefs, constant, &nreturned) );
         SCIProwprepSetSidetype(rowprep, SCIP_SIDETYPE_LEFT);
      }

      for( j = 0; j < nreturned && !*infeasible; ++j )
      {
         SCIP_Bool success;
         int v;

         SCIProwprepReset(rowprep);

         for( v = 0; v < SCIPexprGetNChildren(expr); ++v )
         {
            SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[v]), coefs[j][v]) );
         }
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetExprAuxVarNonlinear(expr), -1.0) );
         SCIProwprepAddConstant(rowprep, constant[j]);  /*lint !e644*/

         /* special treatment for sums to get equality rows */
         if( j == 0 && SCIPisExprSum(scip, expr) )
         {
            SCIP_Real scalefactor;
            SCIP_ROW* row;

            /* improve numerics by scaling only (does not relax inequality) */
            scalefactor = SCIPscaleupRowprep(scip, rowprep, 1.0, &success);
            if( success && scalefactor == 1.0 && underestimate && overestimate )
            {
               /* if the rowprep didn't have to be changed, then turn it into a row, change this to an equality, and add it to the LP */
               /* TODO do this also if not actually needing both under- and overestimator (should still be valid, but also stronger?) */
               (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "initestimate_sum%d", j);

               SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

               /* since we did not relax the estimator, we can turn the row into an equality */
               if( SCIPisInfinity(scip, SCIProwGetRhs(row)) )
               {
                  SCIP_CALL( SCIPchgRowRhs(scip, row, SCIProwGetLhs(row)) );
               }
               else
               {
                  SCIP_CALL( SCIPchgRowLhs(scip, row, SCIProwGetRhs(row)) );
               }
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );

               SCIPdebug( SCIPinfoMessage(scip, NULL, "  added %scut ", *infeasible ? "infeasible " : "") );
               SCIPdebug( SCIPprintRow(scip, row, NULL) );

               SCIP_CALL( SCIPreleaseRow(scip, &row) );

               i = 2;  /* to break outside loop on i, too */
               break;
            }
         }

         /* straighten out numerics */
         SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, SCIPgetHugeValue(scip), &success) );

         /* if cleanup removed all but one variable, then the cut is essentially a bound; we can skip this and rely on boundtightening */
         if( success && SCIProwprepGetNVars(rowprep) > 1 )
         {
            /* add the cut */
            SCIP_ROW* row;

            (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "init%sestimate%d_%s",
                  i == 0 ? "under" : "over", j, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));

            SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );

            SCIPdebug( SCIPinfoMessage(scip, NULL, "  added %scut ", *infeasible ? "infeasible " : "") );
            SCIPdebug( SCIPprintRow(scip, row, NULL) );

            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }
      }
   }

   SCIPfreeRowprep(scip, &rowprep);

   for( i = SCIP_EXPR_MAXINITESTIMATES-1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &coefs[i]);
   }

   SCIPfreeBufferArray(scip, &childrenbounds);

   return SCIP_OKAY;
}

/** compute linear estimator */
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateDefault)
{ /*lint --e{715}*/
   SCIP_Real constant;
   SCIP_Bool local;
   SCIP_Bool* branchcand = NULL;
   int nchildren;
   int c;
   SCIP_INTERVAL* localbounds;
   SCIP_INTERVAL* globalbounds;
   SCIP_Real* refpoint;
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* auxvar;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(rowpreps != NULL);
   assert(success != NULL);

   *addedbranchscores = FALSE;

   nchildren = SCIPexprGetNChildren(expr);

   SCIP_CALL( SCIPallocBufferArray(scip, &localbounds, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &globalbounds, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &refpoint, nchildren) );
   /* we need to pass a branchcand array to exprhdlr's estimate also if not asked to add branching scores */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchcand, nchildren) );

   SCIPdebug( SCIPinfoMessage(scip, NULL, "estimate exprhdlr %s for expr ", SCIPexprhdlrGetName(SCIPexprGetHdlr(expr))) );
   SCIPdebug( SCIPprintExpr(scip, expr, NULL) );
   SCIPdebug( SCIPinfoMessage(scip, NULL, "\n") );

   for( c = 0; c < nchildren; ++c )
   {
      auxvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[c]);
      assert(auxvar != NULL);

      SCIPintervalSetBounds(&localbounds[c],
         -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -SCIPvarGetLbLocal(auxvar)),
          infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  SCIPvarGetUbLocal(auxvar)));

      if( ((size_t)nlhdlrexprdata & (overestimate ? OVERESTIMATEUSESACTIVITY : UNDERESTIMATEUSESACTIVITY)) )
      {
         /* if expr estimate uses bounds, then intersect the auxvar bounds with the current activity, in case the latter is a bit tighter */
         SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(expr)[c]) );
         SCIPintervalIntersectEps(&localbounds[c], SCIPepsilon(scip), localbounds[c], SCIPexprGetActivity(SCIPexprGetChildren(expr)[c]));

         if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, localbounds[c]) )
         {
            *success = FALSE;
            goto TERMINATE;
         }
      }
      else
      {
         /* if we think that expr estimate wouldn't use bounds, then just set something valid */
      }

      SCIPintervalSetBounds(&globalbounds[c],
         -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -SCIPvarGetLbGlobal(auxvar)),
          infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  SCIPvarGetUbGlobal(auxvar)));

      refpoint[c] = SCIPgetSolVal(scip, sol, auxvar);

      branchcand[c] = TRUE;
   }

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );

   /* make sure enough space is available in rowprep arrays */
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nchildren) );

   /* call the estimation callback of the expression handler */
   SCIP_CALL( SCIPcallExprEstimate(scip, expr, localbounds, globalbounds, refpoint, overestimate, targetvalue,
         SCIProwprepGetCoefs(rowprep), &constant, &local, success, branchcand) );

   if( *success )
   {
      int i;

      SCIProwprepSetLocal(rowprep, local);

      /* add variables to rowprep (coefs were already added by SCIPexprhdlrEstimateExpr) */
      for( i = 0; i < nchildren; ++i )
      {
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(expr)[i]),
               SCIProwprepGetCoefs(rowprep)[i]) );
      }

      SCIProwprepAddConstant(rowprep, constant);

      SCIPdebug( SCIPinfoMessage(scip, NULL, "  found rowprep ") );
      SCIPdebug( SCIPprintRowprepSol(scip, rowprep, sol, NULL) );

      SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );

      (void) SCIPsnprintf(SCIProwprepGetName(rowprep), SCIP_MAXSTRLEN, "%sestimate_%s%p_%s%" SCIP_LONGINT_FORMAT,
         overestimate ? "over" : "under",
         SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)),
         (void*)expr,
         sol != NULL ? "sol" : "lp",
         sol != NULL ? (SCIP_Longint) SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));
   }
   else
   {
      SCIPfreeRowprep(scip, &rowprep);
   }

   if( addbranchscores )
   {
      SCIP_Real violation;

#ifndef BRSCORE_ABSVIOL
      SCIP_CALL( SCIPgetExprRelAuxViolationNonlinear(scip, expr, auxvalue, sol, &violation, NULL, NULL) );
#else
      SCIP_CALL( SCIPgetExprAbsAuxViolationNonlinear(scip, expr, auxvalue, sol, &violation, NULL, NULL) );
#endif
      assert(violation > 0.0);  /* there should be a violation if we were called to enforce */

      if( nchildren == 1 )
      {
         if( branchcand[0] )
         {
            SCIP_CALL( SCIPaddExprsViolScoreNonlinear(scip, SCIPexprGetChildren(expr), 1, violation, sol, addedbranchscores) );
         }
      }
      else
      {
         SCIP_EXPR** exprs;
         int nexprs = 0;

         /* get list of those children that have the branchcand-flag set */
         SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nchildren) );

         for( c = 0; c < nchildren; ++c )
            if( branchcand[c] )
               exprs[nexprs++] = SCIPexprGetChildren(expr)[c];

         SCIP_CALL( SCIPaddExprsViolScoreNonlinear(scip, exprs, nexprs, violation, sol, addedbranchscores) );

         SCIPfreeBufferArray(scip, &exprs);
      }

      if( *addedbranchscores )
      {
         /* count this branchscore as belonging to the exprhdlr, too
          * thus, it will be counted for the default nlhdlr, but also for this exprhdlr
          */
         SCIPexprhdlrIncrementNBranchings(SCIPexprGetHdlr(expr));
      }
   }

TERMINATE:
   SCIPfreeBufferArray(scip, &branchcand);
   SCIPfreeBufferArray(scip, &refpoint);
   SCIPfreeBufferArray(scip, &globalbounds);
   SCIPfreeBufferArray(scip, &localbounds);

   return SCIP_OKAY;
}

/** interval-evaluate expression w.r.t. activity of children */
static
SCIP_DECL_NLHDLRINTEVAL(nlhdlrIntevalDefault)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* call the interval evaluation callback of the expression handler */
   SCIP_CALL( SCIPcallExprInteval(scip, expr, interval, intevalvar, intevalvardata) );

   return SCIP_OKAY;
}

/** tighten bounds on children from bounds on expression and bounds on children */
static
SCIP_DECL_NLHDLRREVERSEPROP(nlhdlrReversepropDefault)
{ /*lint --e{715}*/
   SCIP_INTERVAL* childrenbounds;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   *nreductions = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &childrenbounds, SCIPexprGetNChildren(expr)) );
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      childrenbounds[c] = SCIPgetExprBoundsNonlinear(scip, SCIPexprGetChildren(expr)[c]);

   /* call the reverse propagation callback of the expression handler */
   SCIP_CALL( SCIPcallExprReverseprop(scip, expr, bounds, childrenbounds, infeasible) );

   if( !*infeasible )
   {
      for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      {
         SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, SCIPexprGetChildren(expr)[c], childrenbounds[c],
               infeasible, nreductions) );
      }
      SCIPexprhdlrIncrementNDomainReductions(SCIPexprGetHdlr(expr), *nreductions);
   }

   SCIPfreeBufferArray(scip, &childrenbounds);

   return SCIP_OKAY;
}

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrDefault)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrDefault(targetscip) );

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataDefault)
{  /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   *nlhdlrexprdata = NULL;

   return SCIP_OKAY;
}

/** includes default nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrDefault(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
         NLHDLR_ENFOPRIORITY, nlhdlrDetectDefault, nlhdlrEvalAuxDefault, NULL) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrDefault);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataDefault);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaDefault, NULL, nlhdlrEstimateDefault, NULL);
   SCIPnlhdlrSetProp(nlhdlr, nlhdlrIntevalDefault, nlhdlrReversepropDefault);

   return SCIP_OKAY;
}
