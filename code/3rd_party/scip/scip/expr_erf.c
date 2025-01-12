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

/**@file   expr_erf.c
 * @brief  handler for Gaussian error function expressions
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_erf.h"
#include "scip/expr_value.h"

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "erf"
#define EXPRHDLR_DESC         "Gaussian error function"
#define EXPRHDLR_PRECEDENCE   79000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(131071.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/** evaluates the Gaussian error function at a given point */
static
SCIP_Real errorf(
   SCIP_Real             x                   /**< point to evaluate */
   )
{
   SCIP_Real a1 = +0.254829592;
   SCIP_Real a2 = -0.284496736;
   SCIP_Real a3 = +1.421413741;
   SCIP_Real a4 = -1.453152027;
   SCIP_Real a5 = +1.061405429;
   SCIP_Real p  = +0.3275911;
   int sign  = (x >= 0.0) ? 1 : -1;
   SCIP_Real t = 1.0 / (1.0 + p * REALABS(x));
   SCIP_Real y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x*x);

   return sign*y;
}

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrErf)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrErf(scip) );

   return SCIP_OKAY;
}

/** simplifies an erf expression */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyErf)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPisExprValue(scip, child) )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, errorf(SCIPgetValueExprValue(child)), ownercreate, ownercreatedata) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseErf)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   /* create gaussian error function expression */
   SCIP_CALL( SCIPcreateExprErf(scip, expr, childexpr, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the gaussian error function expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_EXPREVAL(evalErf)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) == NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = errorf(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffErf)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of erf expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalErf)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) == NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
   {
      SCIP_Real childinf = SCIPintervalGetInf(childinterval);
      SCIP_Real childsup = SCIPintervalGetSup(childinterval);
      SCIP_Real inf = childinf <= -SCIP_INTERVAL_INFINITY ? -1.0 : errorf(childinf);
      SCIP_Real sup = childsup >= +SCIP_INTERVAL_INFINITY ? +1.0 : errorf(childsup);
      assert(inf <= sup);
      SCIPintervalSetBounds(interval, inf, sup);
   }

   return SCIP_OKAY;
}

/** erf hash callback */
static
SCIP_DECL_EXPRHASH(hashErf)
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
SCIP_DECL_EXPRCURVATURE(curvatureErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* expression is
    *  - convex if child is convex and child <= 0
    *  - concave if child is concave and child >= 0
    */
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
SCIP_DECL_EXPRMONOTONICITY(monotonicityErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);

   *result = SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityErf)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   *isintegral = FALSE;

   return SCIP_OKAY;
}

/** creates an erf expression
 *
 * @attention The implementation of `erf` expressions is incomplete.
 * They are not usable for most use cases so far.
 */
SCIP_RETCODE SCIPcreateExprErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< child expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   assert(expr != NULL);
   assert(child != NULL);

   exprhdlr = SCIPfindExprhdlr(scip, EXPRHDLR_NAME);
   if( exprhdlr == NULL )
   {
      SCIPerrorMessage("could not find %s expression handler -> abort\n", EXPRHDLR_NAME);
      SCIPABORT();
      return SCIP_PLUGINNOTFOUND;
   }

   /* create expression */
   SCIP_CALL( SCIPcreateExpr(scip, expr, exprhdlr, NULL, 1, &child, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** indicates whether expression is of erf-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprErf(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0;
}

/** creates the handler for erf expressions and includes it into SCIP
 *
 * @attention The implementation of this expression handler is incomplete.
 * It is not usable for most use cases so far.
 */
SCIP_RETCODE SCIPincludeExprhdlrErf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalErf, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrErf, NULL);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyErf);
   SCIPexprhdlrSetParse(exprhdlr, parseErf);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalErf);
   SCIPexprhdlrSetHash(exprhdlr, hashErf);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffErf, NULL, NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureErf);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityErf);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityErf);

   return SCIP_OKAY;
}
