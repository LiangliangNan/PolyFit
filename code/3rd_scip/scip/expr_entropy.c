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

/**@file   expr_entropy.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  handler for -x*log(x) expressions
 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 * @author Ksenia Bestuzheva
 *
 * @todo replace exp(-1.0) by 1.0/M_E
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_entropy.h"
#include "scip/expr_value.h"
#include "scip/expr.h"

#include <string.h>

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "entropy"
#define EXPRHDLR_DESC         "entropy expression (-x*log(x))"
#define EXPRHDLR_PRECEDENCE   81000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(7477.0)

/*
 * Data structures
 */

/*
 * Local methods
 */

/** helper function for reverseProp() which returns an x* in [xmin,xmax] s.t. the distance -x*log(x) and a given target
 *  value is minimized; the function assumes that -x*log(x) is monotone on [xmin,xmax];
 */
static
SCIP_Real reversePropBinarySearch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             xmin,               /**< smallest possible x */
   SCIP_Real             xmax,               /**< largest possible x */
   SCIP_Bool             increasing,         /**< -x*log(x) is increasing or decreasing on [xmin,xmax] */
   SCIP_Real             targetval           /**< target value */
   )
{
   SCIP_Real xminval = (xmin == 0.0) ? 0.0 : -xmin * log(xmin);
   SCIP_Real xmaxval = (xmax == 0.0) ? 0.0 : -xmax * log(xmax);
   int i;

   assert(xmin <= xmax);
   assert(increasing ? xminval <= xmaxval : xminval >= xmaxval);

   /* function can not achieve -x*log(x) -> return xmin or xmax */
   if( SCIPisGE(scip, xminval, targetval) && SCIPisGE(scip, xmaxval, targetval) )
      return increasing ? xmin : xmax;
   else if( SCIPisLE(scip, xminval, targetval) && SCIPisLE(scip, xmaxval, targetval) )
      return increasing ? xmax : xmin;

   /* binary search */
   for( i = 0; i < 1000; ++i )
   {
      SCIP_Real x = (xmin + xmax) / 2.0;
      SCIP_Real xval = (x == 0.0) ? 0.0 : -x * log(x);

      /* found the corresponding point -> skip */
      if( SCIPisEQ(scip, xval, targetval) )
         return x;
      else if( SCIPisLT(scip, xval, targetval) )
      {
         if( increasing )
            xmin = x;
         else
            xmax = x;
      }
      else
      {
         if( increasing )
            xmax = x;
         else
            xmin = x;
      }
   }

   return SCIP_INVALID;
}

/** helper function for reverse propagation; needed for proper unittest */
static
SCIP_RETCODE reverseProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         exprinterval,       /**< bounds on the expression */
   SCIP_INTERVAL         childinterval,      /**< bounds on the interval of the child */
   SCIP_INTERVAL*        interval            /**< resulting interval */
   )
{
   SCIP_INTERVAL childentropy;
   SCIP_INTERVAL intersection;
   SCIP_INTERVAL tmp;
   SCIP_Real childinf;
   SCIP_Real childsup;
   SCIP_Real extremum;
   SCIP_Real boundinf;
   SCIP_Real boundsup;

   assert(scip != NULL);
   assert(interval != NULL);

   /* check whether domain is empty, i.e., bounds on -x*log(x) > 1/e */
   if( SCIPisGT(scip, SCIPintervalGetInf(exprinterval), exp(-1.0))
      || SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
   {
      SCIPintervalSetEmpty(interval);
      return SCIP_OKAY;
   }

   /* compute the intersection between entropy([childinf,childsup]) and [expr.inf, expr.sup] */
   SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, &childentropy, childinterval);
   SCIPintervalIntersect(&intersection, childentropy, exprinterval);

   /* intersection empty -> infeasible */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, intersection) )
   {
      SCIPintervalSetEmpty(interval);
      return SCIP_OKAY;
   }

   /* intersection = childentropy -> nothing can be learned */
   if( SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, childentropy, intersection) )
   {
      SCIPintervalSetBounds(interval, 0.0, SCIP_INTERVAL_INFINITY);
      SCIPintervalIntersect(interval, *interval, childinterval);
      return SCIP_OKAY;
   }

   childinf = MAX(0.0, SCIPintervalGetInf(childinterval)); /*lint !e666*/
   childsup = SCIPintervalGetSup(childinterval);
   extremum = exp(-1.0);
   boundinf = SCIP_INVALID;
   boundsup = SCIP_INVALID;

   /*
    * check whether lower bound of child can be improved
    */
   SCIPintervalSet(&tmp, childinf);
   SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, &tmp, tmp);

   /* entropy(childinf) < intersection.inf -> consider [childinf, MIN(childsup, extremum)] */
   if( SCIPintervalGetInf(intersection) > -SCIP_INTERVAL_INFINITY && SCIPintervalGetSup(tmp) - SCIPintervalGetInf(intersection) < -SCIPepsilon(scip) )
   {
      boundinf = reversePropBinarySearch(scip, childinf, MIN(extremum, childsup), TRUE,
         SCIPintervalGetInf(intersection));
   }
   /* entropy(childinf) > intersection.sup -> consider [MAX(childinf,extremum), childsup] */
   else if( SCIPintervalGetSup(intersection) < SCIP_INTERVAL_INFINITY && SCIPintervalGetInf(tmp) - SCIPintervalGetSup(intersection) > SCIPepsilon(scip) )
   {
      boundinf = reversePropBinarySearch(scip, MAX(childinf, extremum), childsup, FALSE,
         SCIPintervalGetSup(intersection));
   }
   /* using a strict greater-than here because we expect a tightening because we saw an at-least-epsilon-potential above */
   assert(boundinf == SCIP_INVALID || boundinf > childinf); /*lint !e777*/

   /*
    * check whether upper bound of child can be improved
    */
   if( childsup < SCIP_INTERVAL_INFINITY )
   {
      SCIPintervalSet(&tmp, childsup);
      SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, &tmp, tmp);
   }
   else
      SCIPintervalSetBounds(&tmp, -SCIP_INTERVAL_INFINITY, -SCIP_INTERVAL_INFINITY);  /* entropy(inf) = -inf */

   /* entropy(childsup) < intersection.inf -> consider [MAX(childinf,extremum), childsup] */
   if( SCIPintervalGetInf(intersection) > -SCIP_INTERVAL_INFINITY && SCIPintervalGetSup(tmp) - SCIPintervalGetInf(intersection) < -SCIPepsilon(scip) )
   {
      boundsup = reversePropBinarySearch(scip, MAX(childinf, extremum), childsup, FALSE,
         SCIPintervalGetInf(intersection));
   }
   /* entropy(childsup) > intersection.sup -> consider [childinf, MIN(childsup,extremum)] */
   else if( SCIPintervalGetSup(intersection) < SCIP_INTERVAL_INFINITY && SCIPintervalGetInf(tmp) - SCIPintervalGetSup(intersection) > SCIPepsilon(scip) )
   {
      boundsup = reversePropBinarySearch(scip, childinf, MIN(childsup, extremum), TRUE,
         SCIPintervalGetSup(intersection));
   }
   /* using a strict smaller-than here because we expect a tightening because we saw an at-least-epsilon-potential above */
   assert(boundsup == SCIP_INVALID || boundsup < childsup); /*lint !e777*/

   if( boundinf != SCIP_INVALID ) /*lint !e777*/
   {
      childinf = MAX(childinf, boundinf);
   }
   if( boundsup != SCIP_INVALID ) /*lint !e777*/
   {
      childsup = boundsup;
   }
   assert(childinf <= childsup); /* infeasible case has been handled already */

   /* set the resulting bounds */
   SCIPintervalSetBounds(interval, childinf, childsup);

   return SCIP_OKAY;
}

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrEntropy)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrEntropy(scip) );

   return SCIP_OKAY;
}

/** simplifies an entropy expression */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyEntropy)
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
      SCIP_Real childvalue = SCIPgetValueExprValue(child);

      /* TODO how to handle a negative value? */
      assert(childvalue >= 0.0);

      if( childvalue == 0.0 || childvalue == 1.0 )
      {
         SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, 0.0, ownercreate, ownercreatedata) );
      }
      else
      {
         SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, -childvalue * log(childvalue), ownercreate,
               ownercreatedata) );
      }
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
   }

   /* TODO handle -x*log(x) = 0 if x in {0,1} */

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataEntropy)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPexprGetData(sourceexpr) == NULL);

   *targetexprdata = NULL;
   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataEntropy)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPexprSetData(expr, NULL);
   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseEntropy)
{  /*lint --e{715}*/
   SCIP_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseExpr(scip, &childexpr, string, endstring, ownercreate, ownercreatedata) );
   assert(childexpr != NULL);

   /* create entropy expression */
   SCIP_CALL( SCIPcreateExprEntropy(scip, expr, childexpr, ownercreate, ownercreatedata) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the entropy expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}


/** expression (point-) evaluation callback */
static
SCIP_DECL_EXPREVAL(evalEntropy)
{  /*lint --e{715}*/
   SCIP_Real childvalue;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) == NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   childvalue = SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[0]);

   if( childvalue < 0.0 )
   {
      SCIPdebugMsg(scip, "invalid evaluation of entropy expression\n");
      *val = SCIP_INVALID;
   }
   else if( childvalue == 0.0 || childvalue == 1.0 )
   {
      /* -x*log(x) = 0 iff x in {0,1} */
      *val = 0.0;
   }
   else
   {
      *val = -childvalue * log(childvalue);
   }

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffEntropy)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_Real childvalue;

   assert(expr != NULL);
   assert(childidx == 0);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);
   assert(!SCIPisExprValue(scip, child));

   childvalue = SCIPexprGetEvalValue(child);

   /* derivative is not defined for x = 0 */
   if( childvalue <= 0.0 )
      *val = SCIP_INVALID;
   else
      *val = -1.0 - log(childvalue);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalEntropy)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPexprGetData(expr) == NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[0]);

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      SCIPintervalSetEmpty(interval);
   else
      SCIPintervalEntropy(SCIP_INTERVAL_INFINITY, interval, childinterval);

   return SCIP_OKAY;
}

/** expression estimator callback */
static
SCIP_DECL_EXPRESTIMATE(estimateEntropy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(localbounds != NULL);
   assert(globalbounds != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(refpoint != NULL);
   assert(coefs != NULL);
   assert(constant != NULL);
   assert(islocal != NULL);
   assert(branchcand != NULL);
   assert(*branchcand == TRUE);
   assert(success != NULL);

   *success = FALSE;

   /* use secant for underestimate (locally valid) */
   if( !overestimate )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real vallb;
      SCIP_Real valub;

      lb = localbounds[0].inf;
      ub = localbounds[0].sup;

      if( lb < 0.0 || SCIPisInfinity(scip, ub) || SCIPisEQ(scip, lb, ub) )
         return SCIP_OKAY;

      assert(lb >= 0.0 && ub >= 0.0);
      assert(ub - lb != 0.0);

      vallb = (lb == 0.0) ? 0.0 : -lb * log(lb);
      valub = (ub == 0.0) ? 0.0 : -ub * log(ub);

      coefs[0] = (valub - vallb) / (ub - lb);
      *constant = valub - coefs[0] * ub;
      assert(SCIPisEQ(scip, *constant, vallb - coefs[0] * lb));

      *islocal = TRUE;
   }
   /* use gradient cut for overestimate (globally valid) */
   else
   {
      if( !SCIPisPositive(scip, refpoint[0]) )
      {
         /* if refpoint is 0 (then lb=0 probably) or negative, then slope is infinite (or not defined), then try to move away from 0 */
         if( SCIPisZero(scip, localbounds[0].sup) )
            return SCIP_OKAY;

         refpoint[0] = SCIPepsilon(scip);
      }

      /* -x*(1+log(x*)) + x* <= -x*log(x) */
      coefs[0] = -(1.0 + log(refpoint[0]));
      *constant = refpoint[0];

      *islocal = FALSE;
      *branchcand = FALSE;
   }

   /* give up if the constant or coefficient is too large */
   if( SCIPisInfinity(scip, REALABS(*constant)) || SCIPisInfinity(scip, REALABS(coefs[0])) )
      return SCIP_OKAY;

   *success = TRUE;

   return SCIP_OKAY;
}

/** initial estimates callback */
static
SCIP_DECL_EXPRINITESTIMATES(initestimatesEntropy)
{  /*lint --e{715}*/
   SCIP_Real refpointsover[3] = {SCIP_INVALID, SCIP_INVALID, SCIP_INVALID};
   SCIP_Bool overest[4] = {TRUE, TRUE, TRUE, FALSE};
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   lb = bounds[0].inf;
   ub = bounds[0].sup;

   if( SCIPisEQ(scip, lb, ub) )
      return SCIP_OKAY;

   if( overestimate )
   {
      /* adjust lb */
      lb = MAX(lb, SCIPepsilon(scip)); /*lint !e666*/

      refpointsover[0] = lb;
      refpointsover[1] = SCIPisInfinity(scip, ub) ? lb + 2.0 : (lb + ub) / 2;
      refpointsover[2] = SCIPisInfinity(scip, ub) ? lb + 20.0 : ub;
   }

   *nreturned = 0;

   for( i = 0; i < 4; ++i )
   {
      if( (overest[i] && !overestimate) || (!overest[i] && (overestimate || SCIPisInfinity(scip, ub))) )
         continue;

      assert(!overest[i] || (SCIPisLE(scip, refpointsover[i], ub) && SCIPisGE(scip, refpointsover[i], lb))); /*lint !e661*/

      if( overest[i] )
      { /*lint !e661*/
         /* -x*(1+log(x*)) + x* <= -x*log(x) */
         assert(i < 3);
         /* coverity[overrun] */
         coefs[*nreturned][0] = -(1.0 + log(refpointsover[i]));
         /* coverity[overrun] */
         constant[*nreturned] = refpointsover[i];
      }
      else
      {
         assert(lb > 0.0 && ub >= 0.0);
         assert(ub - lb != 0.0);

         coefs[*nreturned][0] = (-ub * log(ub) + lb * log(lb)) / (ub - lb);
         constant[*nreturned] = -ub * log(ub) - coefs[*nreturned][0] * ub;
         assert(SCIPisEQ(scip, constant[*nreturned], -lb * log(lb) - coefs[*nreturned][0] * lb));
      }

      ++(*nreturned);
   }

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropEntropy)
{  /*lint --e{715}*/
   SCIP_INTERVAL newinterval;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);
   assert(childrenbounds != NULL);
   assert(infeasible != NULL);

   /* compute resulting intervals (reverseProp handles childinterval being empty) */
   SCIP_CALL( reverseProp(scip, bounds, childrenbounds[0], &newinterval) );
   assert(SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newinterval) || newinterval.inf >= 0.0);

   childrenbounds[0] = newinterval;

   return SCIP_OKAY;
}

/** entropy hash callback */
static
SCIP_DECL_EXPRHASH(hashEntropy)
{  /*lint --e{715}*/
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
SCIP_DECL_EXPRCURVATURE(curvatureEntropy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(childcurv != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 1);

   /* to be concave, the child needs to be concave, too; we cannot be convex or linear */
   if( exprcurvature == SCIP_EXPRCURV_CONCAVE )
   {
      *childcurv = SCIP_EXPRCURV_CONCAVE;
      *success = TRUE;
   }
   else
      *success = FALSE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityEntropy)
{  /*lint --e{715}*/
   SCIP_EXPR* child;
   SCIP_INTERVAL childbounds;
   SCIP_Real brpoint = exp(-1.0);

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx == 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   SCIP_CALL( SCIPevalExprActivity(scip, child) );
   childbounds = SCIPexprGetActivity(child);

   if( childbounds.sup <= brpoint )
      *result = SCIP_MONOTONE_INC;
   else if( childbounds.inf >= brpoint )
      *result = SCIP_MONOTONE_DEC;
   else
      *result = SCIP_MONOTONE_UNKNOWN;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityEntropy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   /* TODO it is possible to check for the special case that the child is integral and its bounds are [0,1]; in
    * this case the entropy expression can only achieve 0 and is thus integral
    */
   *isintegral = FALSE;

   return SCIP_OKAY;
}

/** creates the handler for entropy expressions and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrEntropy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLRDATA* exprhdlrdata;
   SCIP_EXPRHDLR* exprhdlr;

   /* create expression handler data */
   exprhdlrdata = NULL;

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE,
         evalEntropy, exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrEntropy, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataEntropy, freedataEntropy);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyEntropy);
   SCIPexprhdlrSetParse(exprhdlr, parseEntropy);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalEntropy);
   SCIPexprhdlrSetEstimate(exprhdlr, initestimatesEntropy, estimateEntropy);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropEntropy);
   SCIPexprhdlrSetHash(exprhdlr, hashEntropy);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffEntropy, NULL ,NULL);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureEntropy);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityEntropy);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityEntropy);

   return SCIP_OKAY;
}

/** creates an entropy expression */
SCIP_RETCODE SCIPcreateExprEntropy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< child expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(child != NULL);

   exprhdlr = SCIPfindExprhdlr(scip, EXPRHDLR_NAME);
   assert(exprhdlr != NULL);

   /* create expression data */
   exprdata = NULL;

   /* create expression */
   SCIP_CALL( SCIPcreateExpr(scip, expr, exprhdlr, exprdata, 1, &child, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** indicates whether expression is of entropy-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprEntropy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0;
}
