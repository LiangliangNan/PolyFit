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

/**@file   expr_varidx.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  handler for variable index expressions
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_varidx.h"
#include "scip/intervalarith.h"
#include "scip/pub_expr.h"
#include "scip/scip_expr.h"
#include "scip/scip_message.h"
#include "scip/pub_misc.h"

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "varidx"
#define EXPRHDLR_DESC         "expression that represents a variable index (typically used for NLPI)"
#define EXPRHDLR_PRECEDENCE   0
#define EXPRHDLR_HASHKEY      20201210

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrVaridx)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrVaridx(scip) );

   return SCIP_OKAY;
}

/** expression compare callback */
static
SCIP_DECL_EXPRCOMPARE(compareVaridx)
{  /*lint --e{715}*/
   int idx1, idx2;

   assert(expr1 != NULL);
   assert(expr2 != NULL);

   idx1 = SCIPgetIndexExprVaridx(expr1);
   idx2 = SCIPgetIndexExprVaridx(expr2);

   if( idx1 < idx2 )
      return -1;
   if( idx1 > idx2 )
      return 1;
   return 0;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataVaridx)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);

   *targetexprdata = SCIPexprGetData(sourceexpr);

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataVaridx)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_EXPRPRINT(printVaridx)
{  /*lint --e{715}*/
   assert(expr != NULL);

   if( stage == SCIP_EXPRITER_ENTEREXPR )
   {
      SCIPinfoMessage(scip, file, "x%d", SCIPgetIndexExprVaridx(expr));
   }

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_EXPREVAL(evalVaridx)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of varidx expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression backward derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffVaridx)
{  /*lint --e{715}*/
   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression forward derivative evaluation callback */
static
SCIP_DECL_EXPRFWDIFF(fwdiffVaridx)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of varidx expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression backward-forward derivative evaluation callback */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffVaridx)
{  /*lint --e{715}*/
   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** varidx hash callback */
static
SCIP_DECL_EXPRHASH(hashVaridx)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(hashkey != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash((SCIP_Real)SCIPgetIndexExprVaridx(expr));

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureVaridx)
{  /*lint --e{715}*/
   assert(success != NULL);

   /* x -> x is linear, convex, and concave */
   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityVaridx)
{  /*lint --e{715}*/
   assert(result != NULL);

   *result = SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** creates the handler for variable index expressions and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrVaridx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLRDATA* exprhdlrdata;
   SCIP_EXPRHDLR* exprhdlr;

   /* create expression handler data */
   exprhdlrdata = NULL;

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalVaridx,
         exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrVaridx, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataVaridx, freedataVaridx);
   SCIPexprhdlrSetCompare(exprhdlr, compareVaridx);
   SCIPexprhdlrSetPrint(exprhdlr, printVaridx);
   SCIPexprhdlrSetHash(exprhdlr, hashVaridx);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffVaridx, fwdiffVaridx, bwfwdiffVaridx);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureVaridx);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityVaridx);

   return SCIP_OKAY;
}

/** creates a variable index expression */
SCIP_RETCODE SCIPcreateExprVaridx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   varidx,             /**< variable index to represent */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprhdlr = SCIPfindExprhdlr(scip, EXPRHDLR_NAME);

   if( exprhdlr == NULL )
   {
      SCIPerrorMessage("could not find %s expression handler -> abort\n", EXPRHDLR_NAME);
      SCIPABORT();
      return SCIP_ERROR;
   }

   /* create expression data */
   exprdata = (SCIP_EXPRDATA*)(size_t)varidx;  /*lint !e571*/

   /* create expression */
   SCIP_CALL( SCIPcreateExpr(scip, expr, exprhdlr, exprdata, 0, NULL, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** indicates whether expression is varidx expression */  /*lint -e{715}*/
SCIP_Bool SCIPisExprVaridx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{ /*lint --e{715}*/
   assert(expr != NULL);

   /* quick inconclusive check first */
   if( SCIPexprGetNChildren(expr) > 0 )
      return FALSE;

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0;
}

/** gives the index stored in a varidx expression */
int SCIPgetIndexExprVaridx(
   SCIP_EXPR*            expr                /**< varindex expression */
   )
{
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   return (int)(size_t)SCIPexprGetData(expr);
}

/** sets the index stored in a varidx expression */
void SCIPsetIndexExprVaridx(
   SCIP_EXPR*            expr,               /**< varindex expression */
   int                   newindex            /**< new index */
   )
{
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(newindex >= 0);

   SCIPexprSetData(expr, (SCIP_EXPRDATA*)(size_t)newindex);  /*lint !e571*/
}
