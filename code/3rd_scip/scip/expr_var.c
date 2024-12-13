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

/**@file   expr_var.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  variable expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "scip/expr_var.h"
#include "scip/expr_sum.h"

#ifdef NDEBUG
#undef SCIPgetVarExprVar
#endif

#define EXPRHDLR_NAME         "var"
#define EXPRHDLR_DESC         "SCIP variable expression"
#define EXPRHDLR_PRECEDENCE   0
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(22153.0)

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))


/** simplifies a variable expression
 *
 * We replace the variable when fixed by its value.
 * If a variable is fixed, (multi)aggregated or more generally, inactive, we replace it with its active counterpart
 *
 * Implementation note:
 * - we follow the general approach of the simplify, where we replace the var expression for its
 *   simplified expression only in the current parent. So if we see that there is any performance issue in the simplify
 *   we might have to revisit this decision.
 * - we build the sum expression by appending variable expressions one at a time. This may be
 *   speed-up if we allocate memory for all the variable expressions and build the sum directly.
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real constant;
   int nvars;
   int varssize;
   int i;
   SCIP_EXPR* sumexpr;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);

   var = SCIPgetVarExprVar(expr);
   assert(var != NULL);

   /* if var is active then there is nothing to simplify */
   if( SCIPvarIsActive(var) )
   {
      *simplifiedexpr = expr;
      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
      return SCIP_OKAY;
   }

   /* var is not active; obtain active representation var = constant + sum coefs_i vars_i */
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, varssize) );

   vars[0]  = var;
   coefs[0] = 1.0;
   constant = 0.0;
   nvars = 1;
   if( !SCIPvarIsOriginal(var) )
   {
      int requsize;

      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, coefs, &nvars, varssize, &constant, &requsize, TRUE) );

      if( requsize > varssize )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  requsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, requsize) );
         varssize = requsize;
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, coefs, &nvars, varssize, &constant, &requsize, TRUE) );
         assert(requsize <= nvars);
      }
   }

   /* create expression for constant + sum coefs_i vars_i */
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 0, NULL, NULL, constant, ownercreate, ownercreatedata) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_EXPR* child;

      SCIP_CALL( SCIPcreateExprVar(scip, &child, vars[i], ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPappendExprSumExpr(scip, sumexpr, child, coefs[i]) );
      SCIP_CALL( SCIPreleaseExpr(scip, &child) );
   }

   /* simplify since it might not really be a sum */
   SCIP_CALL( SCIPcallExprSimplify(scip, sumexpr, simplifiedexpr, ownercreate, ownercreatedata) );

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "expr_var simplify: <%s> := ", SCIPvarGetName(var));
   SCIPprintExpr(scip, *simplifiedexpr, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   /* we cannot handle fixings to infinity at the moment (TODO we should) */
   assert(!SCIPisInfinity(scip, REALABS(constant)));

   /* release no longer used sumexpr */
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** the order of two variable is given by their indices
 *
 * @note this is affected by permutations in the problem
 */
static
SCIP_DECL_EXPRCOMPARE(compareVar)
{  /*lint --e{715}*/
   int index1;
   int index2;

   index1 = SCIPvarGetIndex(SCIPgetVarExprVar(expr1));
   index2 = SCIPvarGetIndex(SCIPgetVarExprVar(expr2));

   return index1 < index2 ? -1 : index1 == index2 ? 0 : 1;
}

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrVar)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrVar(scip) );

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   /* copying into a different SCIP should be handled on the SCIPexprCopy() level (via mapexpr) */
   assert(targetscip == sourcescip);

   var = SCIPgetVarExprVar(sourceexpr);
   assert(var != NULL);

   *targetexprdata = (SCIP_EXPRDATA*)var;

   SCIP_CALL( SCIPcaptureVar(targetscip, var) );

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataVar)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIP_CALL( SCIPreleaseVar(scip, (SCIP_VAR**)&exprdata) );

   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_EXPRPRINT(printVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetVarExprVar(expr) != NULL);

   if( stage == SCIP_EXPRITER_ENTEREXPR )
   {
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(SCIPgetVarExprVar(expr)));
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetVarExprVar(expr) != NULL);

   *val = SCIPgetSolVal(scip, sol, SCIPgetVarExprVar(expr));

   return SCIP_OKAY;
}

/** expression backward derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffVar)
{  /*lint --e{715}*/
   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression forward derivative evaluation callback */
static
SCIP_DECL_EXPRFWDIFF(fwdiffVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetVarExprVar(expr) != NULL);

   *dot = SCIPgetSolVal(scip, direction, SCIPgetVarExprVar(expr));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffVar)
{  /*lint --e{715}*/
   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(expr != NULL);

   var = SCIPgetVarExprVar(expr);
   assert(var != NULL);

   if( intevalvar != NULL )
      *interval = intevalvar(scip, var, intevalvardata);
   else
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      SCIPintervalSetBounds(interval,  /*lint !e666*/
         -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb),    /*lint !e666*/
          infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  ub));   /*lint !e666*/
   }

   return SCIP_OKAY;
}

/** variable hash callback */
static
SCIP_DECL_EXPRHASH(hashVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);
   assert(hashkey != NULL);

   var = SCIPgetVarExprVar(expr);
   assert(var != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash((SCIP_Real)SCIPvarGetIndex(var));

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   /* x -> x is linear, convex, and concave */
   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   *result = SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   *isintegral = SCIPvarIsIntegral(SCIPgetVarExprVar(expr));

   return SCIP_OKAY;
}

/** creates the handler for variable expression and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrVar(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalVar, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrVar, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataVar, freedataVar);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyVar);
   SCIPexprhdlrSetCompare(exprhdlr, compareVar);
   SCIPexprhdlrSetPrint(exprhdlr, printVar);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalVar);
   SCIPexprhdlrSetHash(exprhdlr, hashVar);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffVar, fwdiffVar, bwfwdiffVar);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureVar);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityVar);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityVar);

   return SCIP_OKAY;
}

/** creates a variable expression */
SCIP_RETCODE SCIPcreateExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_VAR*             var,                /**< variable to be stored */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(var != NULL);

   /* capture the variable so that it doesn't disappear while the expr still points to it */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   exprdata = (SCIP_EXPRDATA*)var;

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPgetExprhdlrVar(scip), exprdata, 0, NULL, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/* from pub_expr.h */

/** gets the variable of a variable expression */
SCIP_VAR* SCIPgetVarExprVar(
   SCIP_EXPR*            expr                /**< variable expression */
   )
{
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(SCIPexprGetData(expr) != NULL);

   return (SCIP_VAR*)SCIPexprGetData(expr);
}
