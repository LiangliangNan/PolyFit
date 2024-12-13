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

/**@file   expr_exp.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  exponential expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_E on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>

#include "scip/expr_exp.h"
#include "scip/expr_value.h"

#define EXPRHDLR_NAME         "exp"
#define EXPRHDLR_DESC         "exponential expression"
#define EXPRHDLR_PRECEDENCE   85000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(10181.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/** computes coefficients of secant of an exponential term */
static
void addExpSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of secant */
   SCIP_Real*            linconstant,        /**< buffer to add constant of secant */
   SCIP_Bool*            success             /**< buffer to set to FALSE if secant has failed due to large numbers or unboundedness */
   )
{
   SCIP_Real coef;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(!SCIPisInfinity(scip,  lb));
   assert(!SCIPisInfinity(scip, -ub));
   assert(SCIPisLE(scip, lb, ub));
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
   {
      /* unboundedness */
      *success = FALSE;
      return;
   }

   /* if lb and ub are too close use a safe secant */
   if( SCIPisEQ(scip, lb, ub) )
   {
      coef = 0.0;
      constant = exp(ub);
   }
   else
   {
      coef = (exp(ub) - exp(lb)) / (ub - lb);
      constant = exp(ub) - coef * ub;
   }

   if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoef     += coef;
   *linconstant += constant;
}

/** computes coefficients of linearization of an exponential term in a reference point */
static
void addExpLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             refpoint,           /**< point for which to compute value of linearization */
   SCIP_Bool             isint,              /**< whether corresponding variable is a discrete variable, and thus linearization could be moved */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of secant */
   SCIP_Real*            linconstant,        /**< buffer to add constant of secant */
   SCIP_Bool*            success             /**< buffer to set to FALSE if secant has failed due to large numbers or unboundedness */
   )
{
   SCIP_Real constant;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(lincoef != NULL);
   assert(linconstant != NULL);
   assert(success != NULL);

   if( SCIPisInfinity(scip, REALABS(refpoint)) )
   {
      *success = FALSE;
      return;
   }

   if( !isint || SCIPisIntegral(scip, refpoint) )
   {
      coef = exp(refpoint);
      constant = exp(refpoint) * (1.0 - refpoint);
   }
   else
   {
      /* exp(x) -> secant between f=floor(refpoint) and f+1 = ((e-1)*e^f) * x + e^f - f * ((e-1)*e^f) */
      SCIP_Real f;

      f = SCIPfloor(scip, refpoint);

      coef = (M_E - 1.0) * exp(f);
      constant = exp(f) - f * coef;
   }

   if( SCIPisInfinity(scip, REALABS(coef)) || SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      return;
   }

   *lincoef     += coef;
   *linconstant += constant;
}

/*
 * Callback methods of expression handler
 */

/** simplifies an exp expression
 *
 * Evaluates the exponential function when its child is a value expression.
 *
 * TODO: exp(log(*)) = *
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyExp)
{
   SCIP_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /**! [SnippetExprSimplifyExp] */
   /* check for value expression */
   if( SCIPisExprValue(scip, child) )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, exp(SCIPgetValueExprValue(child)), ownercreate, ownercreatedata) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }
   /**! [SnippetExprSimplifyExp] */

   return SCIP_OKAY;
}

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrExp)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrExp(scip) );

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataExp)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPexprGetData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataExp)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseExp)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;

   assert(expr != NULL);

   /**! [SnippetExprParseExp] */
   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   /* create exponential expression */
   SCIP_CALL( SCIPcreateExprExp(scip, expr, childexpr, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the exponential expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;
   /**! [SnippetExprParseExp] */

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalExp)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) == NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = exp(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffExp)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(childidx == 0);
   assert(!SCIPisExprValue(scip, SCIPexprGetChildren(expr)[0]));
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/

   *val = SCIPexprGetEvalValue(expr);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalExp)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) == NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalExp(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimateExp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   *success = TRUE;
   *coefs = 0.0;
   *constant = 0.0;

   if( overestimate )
   {
      addExpSecant(scip, localbounds[0].inf, localbounds[0].sup, coefs, constant, success);
      *islocal = TRUE; /* secants are only valid locally */
   }
   else
   {
      addExpLinearization(scip, refpoint[0], SCIPexprIsIntegral(SCIPexprGetChildren(expr)[0]), coefs, constant, success);
      *islocal = FALSE; /* linearization are globally valid */
      *branchcand = FALSE;
   }

   return SCIP_OKAY;
}

/** initital estimates callback for an exponential expression */
static
SCIP_DECL_EXPRINITESTIMATES(initestimatesExp)
{
   SCIP_Real refpointsunder[3] = {SCIP_INVALID, SCIP_INVALID, SCIP_INVALID};
   SCIP_Bool overest[4] = {FALSE, FALSE, FALSE, TRUE};
   SCIP_EXPR* child;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool success;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   /* get expression data */
   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   lb = bounds[0].inf;
   ub = bounds[0].sup;

   if( !overestimate )
   {
      SCIP_Real lbfinite;
      SCIP_Real ubfinite;

      /* make bounds finite */
      lbfinite = SCIPisInfinity(scip, -lb) ? MIN(-5.0, ub - 0.1 * REALABS(ub)) : lb; /*lint !e666*/
      ubfinite = SCIPisInfinity(scip, ub) ? MAX( 3.0, lb + 0.1 * REALABS(lb)) : ub; /*lint !e666*/

      refpointsunder[0] = (7.0 * lbfinite + ubfinite) / 8.0;
      refpointsunder[1] = (lbfinite + ubfinite) / 2.0;
      refpointsunder[2] = (lbfinite + 7.0 * ubfinite) / 8.0;
   }

   *nreturned = 0;
   for( i = 0; i < 4; ++i )
   {
      if( !overest[i] && overestimate )
         continue;

      if( overest[i] && (!overestimate || SCIPisInfinity(scip, ub) || SCIPisInfinity(scip, -lb)) )
         continue;

      assert(overest[i] || (SCIPisLE(scip, refpointsunder[i], ub) && SCIPisGE(scip, refpointsunder[i], lb))); /*lint !e661*/

      coefs[*nreturned][0] = 0.0;
      constant[*nreturned] = 0.0;

      success = TRUE;

      if( !overest[i] )
      {
         assert(i < 3);
         /* coverity[overrun] */
         addExpLinearization(scip, refpointsunder[i], SCIPexprIsIntegral(child), coefs[*nreturned], &constant[*nreturned], &success); /*lint !e661*/
      }
      else
         addExpSecant(scip, lb, ub, coefs[*nreturned], &constant[*nreturned], &success);

      if( success )
         ++*nreturned;
   }

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropExp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPintervalGetInf(bounds) >= 0.0);

   if( SCIPintervalGetSup(bounds) <= 0.0 )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* f = exp(c0) -> c0 = log(f) */
   SCIPintervalLog(SCIP_INTERVAL_INFINITY, childrenbounds, bounds);

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_EXPRHASH(hashExp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= childrenhashes[0];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureExp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   /* expression is convex if child is convex; expression cannot be concave or linear */
   if( exprcurvature == SCIP_EXPRCURV_CONVEX )
   {
      *success = TRUE;
      *childcurv = SCIP_EXPRCURV_CONVEX;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityExp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   *result = SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** creates the handler for exponential expressions and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrExp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalExp, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrExp, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataExp, freedataExp);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyExp);
   SCIPexprhdlrSetParse(exprhdlr, parseExp);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalExp);
   SCIPexprhdlrSetEstimate(exprhdlr, initestimatesExp, estimateExp);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropExp);
   SCIPexprhdlrSetHash(exprhdlr, hashExp);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffExp, NULL, NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureExp);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityExp);

   return SCIP_OKAY;
}

/** creates an exponential expression */
SCIP_RETCODE SCIPcreateExprExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindExprhdlr(scip, EXPRHDLR_NAME) != NULL);

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPfindExprhdlr(scip, EXPRHDLR_NAME), NULL, 1, &child, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** indicates whether expression is of exp-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0;
}
