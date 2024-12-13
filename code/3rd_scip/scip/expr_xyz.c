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

/**@file   expr_xyz.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  handler for xyz expressions
 * @author Jane Doe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_xyz.h"
#include "scip/pub_expr.h"

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "xyz"
#define EXPRHDLR_DESC         "expression handler template"
#define EXPRHDLR_PRECEDENCE   0

/*
 * Data structures
 */

/* TODO: fill in the necessary data */

/** expression handler data */
struct SCIP_ExprhdlrData
{
};

/** expression data */
struct SCIP_ExprData
{
};

/*
 * Local methods
 */

/* TODO: put your local methods here, and declare them static */

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression handler free callback */
static
SCIP_DECL_EXPRFREEHDLR(freehdlrXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** simplifies a xyz expression */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression compare callback */
static
SCIP_DECL_EXPRCOMPARE(compareXyz)
{  /*lint --e{715}*/
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return 0;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_EXPRPRINT(printXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression parse callback */
static
SCIP_DECL_EXPRPARSE(parseXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_EXPREVAL(evalXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression backward derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression forward derivative evaluation callback */
static
SCIP_DECL_EXPRFWDIFF(fwdiffXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression backward-forward derivative evaluation callback */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** initial estimates callback */
static
SCIP_DECL_EXPRINITESTIMATES(initEstimatesXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression under/overestimation callback */
static
SCIP_DECL_EXPRESTIMATE(estimateXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** xyz hash callback */
static
SCIP_DECL_EXPRHASH(hashXyz)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureXyz)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityXyz)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityXyz)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   SCIPerrorMessage("method of xyz expression handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** creates the handler for xyz expressions and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLRDATA* exprhdlrdata;
   SCIP_EXPRHDLR* exprhdlr;

   /* create expression handler data */
   exprhdlrdata = NULL;

   /* TODO: create and store expression handler specific data here */

   /* include expression handler */
   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalXyz,
         exprhdlrdata) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrXyz, freehdlrXyz);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataXyz, freedataXyz);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyXyz);
   SCIPexprhdlrSetCompare(exprhdlr, compareXyz);
   SCIPexprhdlrSetPrint(exprhdlr, printXyz);
   SCIPexprhdlrSetParse(exprhdlr, parseXyz);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalXyz);
   SCIPexprhdlrSetEstimate(exprhdlr, initEstimatesXyz, estimateXyz);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropXyz);
   SCIPexprhdlrSetHash(exprhdlr, hashXyz);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffXyz, fwdiffXyz, bwfwdiffXyz);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureXyz);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityXyz);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityXyz);

   return SCIP_OKAY;
}

/** creates a xyz expression */
SCIP_RETCODE SCIPcreateExprXyz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children (can be NULL if nchildren is 0) */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprhdlr = SCIPfindExprhdlr(scip, EXPRHDLR_NAME);

   if( exprhdlr != NULL )
   {
      SCIPerrorMessage("could not find %s expression handler -> abort\n", EXPRHDLR_NAME);
      SCIPABORT();
      return SCIP_ERROR;
   }

   /* create expression data */
   exprdata = NULL;

   /* TODO: create and store expression specific data here */

   /* create expression */
   SCIP_CALL( SCIPcreateExpr(scip, expr, exprhdlr, exprdata, nchildren, children, ownercreate,
         ownercreatedata) );

   return SCIP_OKAY;
}

/** indicates whether expression is of xyz-type */  /*lint -e{715}*/
SCIP_Bool SCIPisExprXyz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);

   return strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0;
}
