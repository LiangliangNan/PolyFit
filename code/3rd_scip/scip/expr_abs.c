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

/**@file   expr_abs.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  absolute expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/expr_value.h"
#include "scip/expr_abs.h"
#include "scip/expr.h"

#define EXPRHDLR_NAME         "abs"
#define EXPRHDLR_DESC         "absolute value expression"
#define EXPRHDLR_PRECEDENCE   70000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(7187.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/** computes both tangent underestimates and secant */
static
SCIP_RETCODE computeCutsAbs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         bounds,             /**< bounds of child */
   SCIP_Bool             overestimate,       /**< whether the expression shall be overestimated or underestimated */
   SCIP_Real**           coefs,              /**< buffer to store coefficients of computed estimators */
   SCIP_Real*            constant,           /**< buffer to store constant of computed estimators */
   int*                  nreturned           /**< buffer to store number of estimators that have been computed */
   )
{
   assert(scip != NULL);

   *nreturned = 0;

   /**! [SnippetExprInitestimatesAbs] */
   if( !overestimate )
   {
      /* compute left tangent -x <= z */
      coefs[*nreturned][0] = -1.0;
      constant[*nreturned] = 0.0;
      (*nreturned)++;

      /* compute right tangent x <= z */
      coefs[*nreturned][0] = 1.0;
      constant[*nreturned] = 0.0;
      (*nreturned)++;
   }

   /* compute secant */
   if( overestimate )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = bounds.inf;
      ub = bounds.sup;

      /* it does not make sense to add a cut if child variable is unbounded or fixed */
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) && !SCIPisEQ(scip, lb, ub) )
      {
         if( !SCIPisPositive(scip, ub) )
         {
            /* z = -x, so add z <= -x here (-x <= z is the underestimator that is added above) */
            coefs[*nreturned][0] = -1.0;
            constant[*nreturned] = 0.0;
            (*nreturned)++;
         }
         else if( !SCIPisNegative(scip, lb) )
         {
            /* z =  x, so add z <= x here (x <= z is the underestimator that is added above) */
            coefs[*nreturned][0] = 1.0;
            constant[*nreturned] = 0.0;
            (*nreturned)++;
         }
         else
         {
            /* z = abs(x), x still has mixed sign */
            SCIP_Real alpha;

            /* let alpha = (|ub|-|lb|) / (ub-lb) then the resulting secant looks like
             *
             * z - |ub| <= alpha * (x - ub)  <=> z <= alpha * x + |ub| - alpha * ub
             */
            alpha = (REALABS(ub) - REALABS(lb)) / (ub - lb);

            coefs[*nreturned][0] = alpha;
            constant[*nreturned] = REALABS(ub) - alpha * ub;
            (*nreturned)++;
         }
      }
   }
   /**! [SnippetExprInitestimatesAbs] */

   return SCIP_OKAY;
}


/*
 * Callback methods of expression handler
 */

/** simplifies an abs expression
 *
 * Evaluates the absolute value function when its child is a value expression.
 *
 * TODO: abs(*) = * if * >= 0 or - * if * < 0
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyAbs)
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
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, REALABS(SCIPgetValueExprValue(child)), ownercreate, ownercreatedata) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrAbs)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrAbs(scip) );

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRPARSE(parseAbs)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   /* create absolute expression */
   SCIP_CALL( SCIPcreateExprAbs(scip, expr, childexpr, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the absolute expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalAbs)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = REALABS(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]));

   return SCIP_OKAY;
}


/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffAbs)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(child)), "val") != 0);

   *val = (SCIPexprGetEvalValue(child) >= 0.0) ? 1.0 : -1.0;

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalAbs)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalAbs(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimateAbs)
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

   SCIPdebugMsg(scip, "%sestimate |child| over locdom=[%g,%g] glbdom=[%g,%g]\n", overestimate ? "over" : "under",
      localbounds[0].inf, localbounds[0].sup, globalbounds[0].inf, globalbounds[0].sup);

   /**! [SnippetExprEstimateAbs] */
   if( !overestimate )
   {
      *constant = 0.0;

      if( refpoint[0] <= 0.0 )
         *coefs = -1.0;
      else
         *coefs = 1.0;

      *islocal = FALSE;
      *branchcand = FALSE;
   }
   else
   {
      /* overestimator */
      SCIP_Real lb;
      SCIP_Real ub;

      lb = localbounds[0].inf;
      ub = localbounds[0].sup;

      if( !SCIPisPositive(scip, ub) )
      {
         /* |x| = -x */
         *coefs = -1.0;
         *constant = 0.0;
         *islocal = SCIPisPositive(scip, globalbounds[0].sup);
         *branchcand = FALSE;
      }
      else if( !SCIPisNegative(scip, lb) )
      {
         /* |x| = x */
         *coefs =  1.0;
         *constant = 0.0;
         *islocal = SCIPisNegative(scip, globalbounds[0].inf);
         *branchcand = FALSE;
      }
      else if( !SCIPisRelEQ(scip, lb, -ub) )
      {
         /* |x| with x having mixed sign and ub+lb does not cancel out -> secant */
         SCIP_Real alpha;

         assert(lb < 0.0);
         assert(ub > 0.0);

         /* let alpha = (|ub|-|lb|) / (ub-lb) = (ub+lb)/(ub-lb)
          * then the resulting secant is -lb + alpha * (x - lb) = -lb - alpha*lb + alpha*x
          */
         alpha = (ub + lb) / (ub - lb);

         *coefs = alpha;
         *constant = -lb - alpha * lb;
         *islocal = TRUE;
      }
      else if( lb == -ub ) /*lint !e777*/
      {
         /* alpha = 0 */
         *coefs = 0.0;
         *constant = -lb;
         *islocal = TRUE;
      }
      else
      {
         *success = FALSE;
         return SCIP_OKAY;
      }
   }
   /**! [SnippetExprEstimateAbs] */

   SCIPdebugMsg(scip, "-> %g * <child> %+g, local=%u branchcand=%u\n", *coefs, *constant, *islocal, *branchcand);

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression estimate initialization callback */
static
SCIP_DECL_EXPRINITESTIMATES(initEstimatesAbs)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   /* compute initial cuts */
   SCIP_CALL( computeCutsAbs(scip, bounds[0], overestimate, coefs, constant, nreturned) );

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropAbs)
{  /*lint --e{715}*/
   SCIP_INTERVAL childbounds;
   SCIP_INTERVAL left;
   SCIP_INTERVAL right;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(bounds.inf >= 0.0);  /* bounds should have been intersected with activity, which is >= 0 */

   /**! [SnippetExprReversepropAbs] */
   /* abs(x) in I -> x \in (-I \cup I) \cap bounds(x) */
   right = bounds;  /* I */
   SCIPintervalSetBounds(&left, -right.sup, -right.inf); /* -I */

   childbounds = childrenbounds[0];
   /* childbounds can be empty here already, but that should work fine here */

   SCIPintervalIntersect(&left, left, childbounds);    /* -I \cap bounds(x), could become empty */
   SCIPintervalIntersect(&right, right, childbounds);  /*  I \cap bounds(x), could become empty */

   /* compute smallest interval containing (-I \cap bounds(x)) \cup (I \cap bounds(x)) = (-I \cup I) \cap bounds(x)
    * this works also if left or right is empty
    */
   SCIPintervalUnify(&childrenbounds[0], left, right);
   /**! [SnippetExprReversepropAbs] */

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_EXPRHASH(hashAbs)
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
SCIP_DECL_EXPRCURVATURE(curvatureAbs)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_INTERVAL childbounds;
   SCIP_Real childinf;
   SCIP_Real childsup;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprcurvature != SCIP_EXPRCURV_UNKNOWN);
   assert(success != NULL);
   assert(childcurv != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /**! [SnippetExprCurvatureAbs] */
   /* expression is |child|, get domain of child */
   SCIP_CALL( SCIPevalExprActivity(scip, child) );
   childbounds = SCIPexprGetActivity(child);
   childinf = SCIPintervalGetInf(childbounds);
   childsup = SCIPintervalGetSup(childbounds);

   *success = TRUE;
   if( childinf >= 0.0 )  /* |f(x)| = f(x) */
      childcurv[0] = exprcurvature;
   else if( childsup <= 0.0 ) /* |f(x)| = -f(x) */
      childcurv[0] = SCIPexprcurvNegate(exprcurvature);
   else if( exprcurvature == SCIP_EXPRCURV_CONVEX )   /* |f(x)|, f mixed sign, is convex if f is linear */
      childcurv[0] = SCIP_EXPRCURV_LINEAR;
   else /* |f(x)|, f mixed sign, is never concave nor linear */
      *success = FALSE;
   /**! [SnippetExprCurvatureAbs] */

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityAbs)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_INTERVAL childbounds;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /**! [SnippetExprMonotonicityAbs] */
   SCIP_CALL( SCIPevalExprActivity(scip, child) );
   childbounds = SCIPexprGetActivity(child);

   if( childbounds.sup <= 0.0 )
      *result = SCIP_MONOTONE_DEC;
   else if( childbounds.inf >= 0.0 )
      *result = SCIP_MONOTONE_INC;
   else
      *result = SCIP_MONOTONE_UNKNOWN;
   /**! [SnippetExprMonotonicityAbs] */

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityAbs)
{  /*lint --e{715}*/
   SCIP_EXPR* child;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   *isintegral = SCIPexprIsIntegral(child);

   return SCIP_OKAY;
}


/** creates the handler for absolute expression and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrAbs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalAbs, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrAbs, NULL);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyAbs);
   SCIPexprhdlrSetParse(exprhdlr, parseAbs);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalAbs);
   SCIPexprhdlrSetEstimate(exprhdlr, initEstimatesAbs, estimateAbs);
   SCIPexprhdlrSetHash(exprhdlr, hashAbs);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropAbs);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffAbs, NULL, NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureAbs);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityAbs);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityAbs);

   return SCIP_OKAY;
}

/** creates an absolute value expression */
SCIP_RETCODE SCIPcreateExprAbs(
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

/** indicates whether expression is of abs-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprAbs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0;
}
