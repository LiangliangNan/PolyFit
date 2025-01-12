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

/**@file   expr_value.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  constant value expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/expr_value.h"

#define EXPRHDLR_NAME            "val"
#define EXPRHDLR_DESC            "constant value"
#define EXPRHDLR_PRECEDENCE      10000
#define EXPRHDLR_HASHKEY         SCIPcalcFibHash(36787.0)

/*
 * Data structures
 */

/** expression data */
struct SCIP_ExprData
{
   SCIP_Real             value;              /**< value that expression represents */
};

/*
 * Callback methods of expression handler
 */

/** the order of two values is the real order */
static
SCIP_DECL_EXPRCOMPARE(compareValue)
{  /*lint --e{715}*/
   SCIP_Real val1;
   SCIP_Real val2;

   assert(SCIPexprGetData(expr1) != NULL);
   assert(SCIPexprGetData(expr2) != NULL);

   val1 = SCIPexprGetData(expr1)->value;
   val2 = SCIPexprGetData(expr2)->value;

   return val1 < val2 ? -1 : val1 == val2 ? 0 : 1; /*lint !e777*/
}

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrValue)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrValue(scip) );

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataValue)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(targetscip, targetexprdata) );
   (*targetexprdata)->value = SCIPexprGetData(sourceexpr)->value;

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataValue)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemory(scip, &exprdata);
   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_EXPRPRINT(printValue)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   if( stage == SCIP_EXPRITER_ENTEREXPR )
   {
      SCIP_Real v = SCIPexprGetData(expr)->value;
      if( v < 0.0 && EXPRHDLR_PRECEDENCE <= parentprecedence )
      {
         SCIPinfoMessage(scip, file, "(%g)", v);
      }
      else
      {
         SCIPinfoMessage(scip, file, "%g", v);
      }
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalValue)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   *val = SCIPexprGetData(expr)->value;

   return SCIP_OKAY;
}

/** expression backward derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffValue)
{  /*lint --e{715}*/
   /* should never be called since value expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression forward derivative evaluation callback */
static
SCIP_DECL_EXPRFWDIFF(fwdiffValue)
{  /*lint --e{715}*/
   assert(expr != NULL);

   *dot = 0.0;

   return SCIP_OKAY;
}

/** derivative evaluation callback for Hessian directions (backward over forward) */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffValue)
{  /*lint --e{715}*/
   /* should never be called since value expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalValue)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   SCIPintervalSet(interval, SCIPexprGetData(expr)->value);

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_EXPRHASH(hashValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);
   assert(hashkey != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(SCIPexprGetData(expr)->value);

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   *result = SCIP_MONOTONE_CONST;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityValue)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   *isintegral = EPSISINT(SCIPexprGetData(expr)->value, 0.0); /*lint !e835 !e666*/

   return SCIP_OKAY;
}

/** creates the handler for constant value expression and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrValue(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE,
         evalValue, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrValue, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataValue, freedataValue);
   SCIPexprhdlrSetCompare(exprhdlr, compareValue);
   SCIPexprhdlrSetPrint(exprhdlr, printValue);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalValue);
   SCIPexprhdlrSetHash(exprhdlr, hashValue);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffValue, fwdiffValue, bwfwdiffValue);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureValue);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityValue);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityValue);

   return SCIP_OKAY;
}

/** creates constant value expression */
SCIP_RETCODE SCIPcreateExprValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_Real             value,              /**< value to be stored */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(SCIPisFinite(value));

   SCIP_CALL( SCIPallocBlockMemory(scip, &exprdata) );
   exprdata->value = value;

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPgetExprhdlrValue(scip), exprdata, 0, NULL, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/* from pub_expr.h */

/** gets the value of a constant value expression */
SCIP_Real SCIPgetValueExprValue(
   SCIP_EXPR*            expr                /**< sum expression */
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   return SCIPexprGetData(expr)->value;
}
