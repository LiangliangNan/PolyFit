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

/**@file   sepa_convexproj.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  convexproj separator
 * @author Felipe Serrano
 *
 * @todo should separator only be run when SCIPallColsInLP is true?
 * @todo check if it makes sense to implement the copy callback
 * @todo add SCIPisStopped(scip) to the condition of time consuming loops
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "blockmemshell/memory.h"
#include "scip/scip_expr.h"
#include "scip/scip_nlpi.h"
#include "scip/expr_varidx.h"
#include "scip/expr_pow.h"
#include "scip/expr_sum.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nlp.h"
#include "scip/pub_sepa.h"
#include "scip/pub_var.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sepa.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/sepa_convexproj.h"
#include <string.h>


#define SEPA_NAME              "convexproj"
#define SEPA_DESC              "separate at projection of point onto convex region"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE      /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                 TRUE      /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXDEPTH             -1      /**< maximum depth at which the separator is applied; -1 means no limit */
#define DEFAULT_NLPITERLIM          250      /**< default NLP iteration limit */

#define VIOLATIONFAC                100      /**< points regarded violated if max violation > VIOLATIONFAC*SCIPfeastol() */

/*
 * Data structures
 */

/** side that makes an nlrow convex */
enum ConvexSide
{
   LHS = 0,                                  /**< left hand side */
   RHS = 1                                   /**< right hand side */
};
typedef enum ConvexSide CONVEXSIDE;

/** separator data
 * it keeps the nlpi which represents the projection problem (see sepa_convexproj.h); it also keeps the convex nlrows
 * and the side which actually makes them convex; when separating, we use the nlpi to compute the projection and then
 * the convex nlrows to compute the actual gradient cuts */
struct SCIP_SepaData
{
   SCIP_NLPI*            nlpi;               /**< nlpi used to create the nlpi problem */
   SCIP_NLPIPROBLEM*     nlpiprob;           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_VAR**            nlpivars;           /**< array containing all variables of the nlpi */
   SCIP_HASHMAP*         var2nlpiidx;        /**< mapping between variables and nlpi indices */
   int                   nlpinvars;          /**< total number of nlpi variables */

   SCIP_Bool             skipsepa;           /**< should separator be skipped? */

   SCIP_NLROW**          nlrows;             /**< convex nlrows */
   CONVEXSIDE*           convexsides;        /**< which sides make the nlrows convex */
   SCIP_Real*            constraintviolation;/**< array storing the violation of constraint by current solution; 0.0 if it is not violated */
   int                   nnlrows;            /**< total number of nlrows */
   int                   nlrowssize;         /**< memory allocated for nlrows, convexsides and nlrowsidx */

   /* parameter */
   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for no limit */
   int                   maxdepth;           /**< maximal depth at which the separator is applied */

   int                   ncuts;              /**< number of cuts generated */
};


/*
 * Local methods
 */

/** clears the sepadata data */
static
SCIP_RETCODE sepadataClear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   assert(sepadata != NULL);

   /* nlrowssize gets allocated first and then its decided whether to create the nlpiprob */
   if( sepadata->nlrowssize > 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &sepadata->constraintviolation, sepadata->nlrowssize);
      SCIPfreeBlockMemoryArray(scip, &sepadata->convexsides, sepadata->nlrowssize);
      SCIPfreeBlockMemoryArray(scip, &sepadata->nlrows, sepadata->nlrowssize);
      sepadata->nlrowssize = 0;
   }

   if( sepadata->nlpiprob != NULL )
   {
      assert(sepadata->nlpi != NULL);

      SCIPfreeBlockMemoryArray(scip, &sepadata->nlpivars, sepadata->nlpinvars);

      SCIPhashmapFree(&sepadata->var2nlpiidx);
      SCIP_CALL( SCIPfreeNlpiProblem(scip, sepadata->nlpi, &sepadata->nlpiprob) );

      sepadata->nlpinvars = 0;
      sepadata->nnlrows = 0;
   }
   assert(sepadata->nlpinvars == 0);
   assert(sepadata->nnlrows == 0);
   assert(sepadata->nlrowssize == 0);

   sepadata->skipsepa = FALSE;

   return SCIP_OKAY;
}

/** computes gradient cut (linearization) of nlrow at projection */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             projection,         /**< point where we compute gradient cut */
   SCIP_NLROW*           nlrow,              /**< constraint for which we generate gradient cut */
   CONVEXSIDE            convexside,         /**< which side makes the nlrow convex */
   SCIP_Real             activity,           /**< activity of constraint at projection */
   SCIP_EXPRITER*        exprit,             /**< expression iterator that can be used */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   char rowname[SCIP_MAXSTRLEN];
   SCIP_SEPADATA* sepadata;
   SCIP_Real gradx0; /* <grad f(x_0), x_0> */
   SCIP_EXPR* expr;
   int i;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(nlrow != NULL);
   assert(row != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   gradx0 = 0.0;

   /* an nlrow has a linear part and expression; ideally one would just build the gradient but we
    * do not know if the different parts share variables or not, so we can't just build the gradient; for this reason
    * we create the row right away and compute the gradients of each part independently and add them to the row; the
    * row takes care to add coeffs corresponding to the same variable when they appear in different parts of the nlrow
    * NOTE: a gradient cut is globally valid whenever the constraint from which it is deduced is globally valid; since
    *       we build the convex relaxation using only globally valid constraints, the cuts are globally valid
    */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "proj_cut_%s_%u", SCIPnlrowGetName(nlrow), ++(sepadata->ncuts));
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, row, sepa, rowname, -SCIPinfinity(scip), SCIPinfinity(scip), TRUE, FALSE ,
            TRUE) );

   SCIP_CALL( SCIPcacheRowExtensions(scip, *row) );

   /* linear part */
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
   {
      gradx0 += SCIPgetSolVal(scip, projection, SCIPnlrowGetLinearVars(nlrow)[i]) * SCIPnlrowGetLinearCoefs(nlrow)[i];
      SCIP_CALL( SCIPaddVarToRow(scip, *row, SCIPnlrowGetLinearVars(nlrow)[i], SCIPnlrowGetLinearCoefs(nlrow)[i]) );
   }

   expr = SCIPnlrowGetExpr(nlrow);
   assert(expr != NULL);

   SCIP_CALL( SCIPevalExprGradient(scip, expr, projection, 0L) );

   SCIP_CALL( SCIPexpriterInit(exprit, expr, SCIP_EXPRITER_DFS, FALSE) );
   for( ; !SCIPexpriterIsEnd(exprit); expr = SCIPexpriterGetNext(exprit) )  /*lint !e441*/ /*lint !e440*/
   {
      SCIP_Real grad;
      SCIP_VAR* var;

      if( !SCIPisExprVar(scip, expr) )
         continue;

      grad = SCIPexprGetDerivative(expr);
      var = SCIPgetVarExprVar(expr);
      assert(var != NULL);

      gradx0 += grad * SCIPgetSolVal(scip, projection, var);
      SCIP_CALL( SCIPaddVarToRow(scip, *row, var, grad) );
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, *row) );

   SCIPdebugPrintf("gradient: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *row, NULL) ) );
   SCIPdebugPrintf("gradient dot x_0: %g\n", gradx0);

   /* gradient cut is f(x_0) - <grad f(x_0), x_0> + <grad f(x_0), x> <= rhs or >= lhs */
   if( convexside == RHS )
   {
      assert(!SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)));
      SCIP_CALL( SCIPchgRowRhs(scip, *row, SCIPnlrowGetRhs(nlrow) - activity + gradx0) );
   }
   else
   {
      assert(convexside == LHS);
      assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));
      SCIP_CALL( SCIPchgRowLhs(scip, *row, SCIPnlrowGetLhs(nlrow) - activity + gradx0) );
   }

   SCIPdebugPrintf("gradient cut: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *row, NULL) ) );

   return SCIP_OKAY;
}

/** set quadratic part of objective function: \f$ \sum_i x_i^2 \f$
 *
 * the objective function is \f$ ||x - x_0||^2 \f$,
 * where \f$ x_0 \f$ is the point to separate; the only part that changes is the term \f$ -2 \langle x_0, x \rangle \f$
 * which is linear and is set every time we want to separate a point, see separateCuts()
 */
static
SCIP_RETCODE setQuadraticObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< the cut separator data */
   )
{
   SCIP_EXPR* exprsum;
   SCIP_EXPR** exprspow;
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->nlpi != NULL);
   assert(sepadata->nlpiprob != NULL);
   assert(sepadata->var2nlpiidx != NULL);
   assert(sepadata->nlpinvars > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &exprspow, sepadata->nlpinvars) );
   for( i = 0; i < sepadata->nlpinvars; i++ )
   {
      SCIP_VAR* var;
      SCIP_EXPR* varexpr;

      var = sepadata->nlpivars[i];
      assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

      SCIP_CALL( SCIPcreateExprVaridx(scip, &varexpr, SCIPhashmapGetImageInt(sepadata->var2nlpiidx, (void*)var), NULL, NULL) );
      SCIP_CALL( SCIPcreateExprPow(scip, &exprspow[i], varexpr, 2.0, NULL, NULL) );
      SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );
   }

   SCIP_CALL( SCIPcreateExprSum(scip, &exprsum, sepadata->nlpinvars, exprspow, NULL, 0.0, NULL, NULL) );

   /* set quadratic part of objective function */
   SCIP_CALL( SCIPsetNlpiObjective(scip, sepadata->nlpi, sepadata->nlpiprob, 0, NULL, NULL, exprsum, 0.0) );

   /* free memory */
   SCIP_CALL( SCIPreleaseExpr(scip, &exprsum) );
   for( i = sepadata->nlpinvars-1; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &exprspow[i]) );
   }
   SCIPfreeBufferArray(scip, &exprspow);

   return SCIP_OKAY;
}

/** projects sol onto convex relaxation (stored in sepadata) and tries to generate gradient cuts at the projection
 *
 * it generates cuts only for the constraints that were violated by the LP solution and are now active or still
 * violated (in case we don't solve to optimality).
 * @todo: store a feasible solution if one is found to use as warmstart
 */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             sol,                /**< solution that should be separated */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SOL*      projection;
   SCIP_Real*     linvals;
   SCIP_Real*     nlpisol;
   int            nlpinvars;
   int            i;
   int*           lininds;
   SCIP_Bool      nlpunstable;
   SCIP_EXPRITER* exprit;

   nlpunstable = FALSE;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(result != NULL);
   assert(sepadata != NULL);
   assert(sepadata->nnlrows > 0);
   assert(sepadata->nlpi != NULL);
   assert(sepadata->nlpinvars > 0);
   assert(sepadata->nlrows != NULL);
   assert(sepadata->nlpiprob != NULL);
   assert(sepadata->var2nlpiidx != NULL);
   assert(sepadata->convexsides != NULL);
   assert(sepadata->constraintviolation != NULL);

   nlpinvars = sepadata->nlpinvars;
   /* set linear part of objective function: \norm(x - x^0)^2 = \norm(x)^2 - \sum 2 * x_i * x^0_i + const
    * we ignore the constant; x0 is `sol`
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nlpinvars) );
   for( i = 0; i < nlpinvars; i++ )
   {
      SCIP_VAR* var;

      var = sepadata->nlpivars[i];
      assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

      lininds[i] = SCIPhashmapGetImageInt(sepadata->var2nlpiidx, (void*)var);
      linvals[i] = - 2.0 * SCIPgetSolVal(scip, sol, var);

      /* if coefficient is too large, don't separate */
      if( SCIPisHugeValue(scip, REALABS(linvals[i])) )
      {
         SCIPdebugMsg(scip, "Don't separate points too close to infinity\n");
         goto CLEANUP;
      }
   }

   /* set linear part of objective function */
   SCIP_CALL( SCIPchgNlpiLinearCoefs(scip, sepadata->nlpi, sepadata->nlpiprob, -1, nlpinvars, lininds, linvals) );

   /* compute the projection onto the convex NLP relaxation */
   SCIP_CALL( SCIPsolveNlpi(scip, sepadata->nlpi, sepadata->nlpiprob,
      .iterlimit = sepadata->nlpiterlimit > 0 ? sepadata->nlpiterlimit : INT_MAX,
      .feastol = SCIPfeastol(scip) / 10.0, /* use tighter tolerances for the NLP solver */
      .opttol = MAX(SCIPfeastol(scip), SCIPdualfeastol(scip))) );  /*lint !e666*/
   SCIPdebugMsg(scip, "NLP solstat = %d\n", SCIPgetNlpiSolstat(scip, sepadata->nlpi, sepadata->nlpiprob));

   /* if solution is feasible, add cuts */
   switch( SCIPgetNlpiSolstat(scip, sepadata->nlpi, sepadata->nlpiprob) )
   {
      case SCIP_NLPSOLSTAT_GLOBOPT:
      case SCIP_NLPSOLSTAT_LOCOPT:
         /* @todo: if solution is optimal, we might as well add the cut <x - P(x_0), x_0 - P(x_0)> <= 0
          * even though this cut is implied by all the gradient cuts of the rows active at the projection,
          * we do not add them all (only the gradient cuts of constraints that violated the LP solution */
      case SCIP_NLPSOLSTAT_FEASIBLE:
         /* generate cuts for violated constraints (at sol) that are active or still violated at the projection, since
          * a suboptimal solution or numerical issues could give a solution of the projection problem where constraints
          * are not active; if the solution of the projection problem is in the interior of the region, we do nothing
          */

         /* get solution: build SCIP_SOL out of nlpi sol */
         SCIP_CALL( SCIPgetNlpiSolution(scip, sepadata->nlpi, sepadata->nlpiprob, &nlpisol, NULL, NULL, NULL, NULL) );
         assert(nlpisol != NULL);

         SCIP_CALL( SCIPcreateSol(scip, &projection, NULL) );
         for( i = 0; i < nlpinvars; i++ )
         {
            SCIP_VAR* var;

            var = sepadata->nlpivars[i];
            assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

            SCIP_CALL( SCIPsetSolVal(scip, projection, var,
                     nlpisol[SCIPhashmapGetImageInt(sepadata->var2nlpiidx, (void *)var)]) );
         }
         SCIPdebug( SCIPprintSol(scip, projection, NULL, TRUE) );

         /** @todo this could just be created inside generateCut and the extra argument removed */
         SCIP_CALL( SCIPcreateExpriter(scip, &exprit) );

         /* check for active or violated constraints */
         for( i = 0; i < sepadata->nnlrows; ++i )
         {
            SCIP_NLROW* nlrow;
            CONVEXSIDE convexside;
            SCIP_Real activity;

            /* ignore constraints that are not violated by `sol` */
            if( SCIPisFeasZero(scip, sepadata->constraintviolation[i]) )
               continue;

            convexside = sepadata->convexsides[i];
            nlrow = sepadata->nlrows[i];
            assert(nlrow != NULL);

            /* check for currently active constraints at projected point */
            SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, projection, &activity) );

            SCIPdebugMsg(scip, "NlRow activity at nlpi solution: %g <= %g <= %g\n", SCIPnlrowGetLhs(nlrow), activity,
                  SCIPnlrowGetRhs(nlrow) );

            /* if nlrow is active or violates the projection, build gradient cut at projection */
            if( (convexside == RHS && SCIPisFeasGE(scip, activity, SCIPnlrowGetRhs(nlrow)))
               || (convexside == LHS && SCIPisFeasLE(scip, activity, SCIPnlrowGetLhs(nlrow))) )
            {
               SCIP_ROW* row;

               SCIP_CALL( generateCut(scip, sepa, projection, nlrow, convexside, activity, exprit,
                        &row) );

               SCIPdebugMsg(scip, "active or violated nlrow: (sols vio: %e)\n", sepadata->constraintviolation[i]);
               SCIPdebug( SCIP_CALL( SCIPprintNlRow(scip, nlrow, NULL) ) );
               SCIPdebugMsg(scip, "cut with efficacy %g generated\n", SCIPgetCutEfficacy(scip, sol, row));
               SCIPdebug( SCIPprintRow(scip, row, NULL) );

               /* add cut if it is efficacious for the point we want to separate (sol) */
               if( SCIPisCutEfficacious(scip, sol, row) )
               {
                  SCIP_Bool infeasible;

                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

                  if( infeasible )
                  {
                     *result = SCIP_CUTOFF;
                     SCIP_CALL( SCIPreleaseRow(scip, &row) );
                     break;
                  }
                  else
                  {
                     *result = SCIP_SEPARATED;
                  }
               }

               /* release the row */
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
         }

         SCIPfreeExpriter(&exprit);

#ifdef SCIP_DEBUG
         {
            SCIP_Real distance;

            /* compute distance between LP sol and its projection (only makes sense when it is optimal) */
            distance = 0.0;
            for( i = 0; i < SCIPgetNNLPVars(scip); ++i )
            {
               SCIP_VAR* var;

               var = SCIPgetNLPVars(scip)[i];
               assert(var != NULL);

               /* assert NLP solution is within the bounds of the variable (only make sense when sol is optimal) */
               if( !SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
                  assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(var), SCIPvarGetNLPSol(var)));
               if( !SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
                  assert(SCIPisFeasLE(scip, SCIPvarGetNLPSol(var), SCIPvarGetUbLocal(var)));

               /*SCIPdebugMsg(scip, "NLP sol (LP sol): %s = %f (%g)\n", SCIPvarGetName(var),
                *     SCIPvarGetNLPSol(var), SCIPgetSolVal(scip, sol, var));
                */

               distance += SQR( SCIPvarGetNLPSol(var) - SCIPgetSolVal(scip, sol, var) );
            }

            SCIPdebugMsg(scip, "NLP objval: %e, distance: %e\n", SCIPgetNLPObjval(scip), distance);
         }
#endif

         /* free solution */
         SCIP_CALL( SCIPfreeSol(scip, &projection) );
         break;

      case SCIP_NLPSOLSTAT_GLOBINFEASIBLE:
      case SCIP_NLPSOLSTAT_LOCINFEASIBLE:
         /* fallthrough;
          * @todo: write what it means to be locinfeasible and why it can't be used to cutoff the node */
      case SCIP_NLPSOLSTAT_UNKNOWN:
         /* unknown... assume numerical issues */
         nlpunstable = TRUE;
         break;

      case SCIP_NLPSOLSTAT_UNBOUNDED:
      default:
         SCIPerrorMessage("Projection NLP is not unbounded by construction, should not get here!\n");
         SCIPABORT();
         nlpunstable = TRUE;
   }

   /* if nlp is detected to be unstable, don't try to separate again */
   if( nlpunstable )
   {
      /* @todo: maybe change objective function to \sum [(x_i - x_i^*)/max(|x_i^*|, 1)]^2
       * or some other scaling when unstable and try again.
       *       maybe free it here */
      sepadata->skipsepa = TRUE;
   }

   /* reset objective */
   BMSclearMemoryArray(linvals, nlpinvars);
   SCIP_CALL( SCIPchgNlpiLinearCoefs(scip, sepadata->nlpi, sepadata->nlpiprob, -1, nlpinvars, lininds, linvals) );

CLEANUP:
   /* free memory */
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);

   return SCIP_OKAY;
}

/** computes the violation and maximum violation of the convex nlrows stored in sepadata wrt sol */
static
SCIP_RETCODE computeMaxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< solution that should be separated */
   SCIP_Real*            maxviolation        /**< buffer to store maximum violation */
   )
{
   SCIP_NLROW*    nlrow;
   int            i;

   assert(sepadata != NULL);
   assert(sepadata->nnlrows > 0);
   assert(sepadata->nlrows != NULL);
   assert(sepadata->convexsides != NULL);
   assert(sepadata->constraintviolation != NULL);

   *maxviolation = 0.0;
   for( i = 0; i < sepadata->nnlrows; i++ )
   {
      SCIP_Real activity;
      SCIP_Real violation;

      nlrow = sepadata->nlrows[i];

      /* get activity of nlrow */
      SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, sol, &activity) );

      /* violation = max{activity - rhs, 0.0} when convex and max{lhs - activity, 0.0} when concave */
      if( sepadata->convexsides[i] == RHS )
      {
         assert(SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONVEX);
         assert(!SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)));

         violation = activity - SCIPnlrowGetRhs(nlrow);
         sepadata->constraintviolation[i] = MAX(violation, 0.0);
      }
      if( sepadata->convexsides[i] == LHS )
      {
         assert(SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONCAVE);
         assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));

         violation = SCIPnlrowGetLhs(nlrow) - activity;
         sepadata->constraintviolation[i] = MAX(violation, 0.0);
      }

      /* compute maximum */
      if( *maxviolation < sepadata->constraintviolation[i] )
         *maxviolation = sepadata->constraintviolation[i];
   }

   SCIPdebugMsg(scip, "Maximum violation %g\n", *maxviolation);

   return SCIP_OKAY;
}


/** stores, from the constraints represented by nlrows, the nonlinear convex ones in sepadata */
static
SCIP_RETCODE storeNonlinearConvexNlrows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW**          nlrows,             /**< nlrows from which to store convex ones */
   int                   nnlrows             /**< number of nlrows */
   )
{
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);

   SCIPdebugMsg(scip, "storing convex nlrows\n");

   sepadata->nlrowssize = nnlrows;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->nlrows), nnlrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->convexsides), nnlrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->constraintviolation), nnlrows) );

   /* count the number of nonlinear convex rows and store them */
   sepadata->nnlrows = 0;
   for( i = 0; i < nnlrows; ++i )
   {
      SCIP_NLROW* nlrow;

      nlrow = nlrows[i];
      assert(nlrow != NULL);

      /* linear case */
      if( SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_LINEAR || SCIPnlrowGetExpr(nlrow) == NULL )
         continue;

      /* nonlinear case */
      if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)) && SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONVEX )
      {
         sepadata->convexsides[sepadata->nnlrows] = RHS;
         sepadata->nlrows[sepadata->nnlrows] = nlrow;
         ++(sepadata->nnlrows);
      }
      else if( !SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)) && SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONCAVE )
      {
         sepadata->convexsides[sepadata->nnlrows] = LHS;
         sepadata->nlrows[sepadata->nnlrows] = nlrow;
         ++(sepadata->nnlrows);
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeConvexproj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( sepadataClear(scip, sepadata) );

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolConvexproj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   SCIP_CALL( sepadataClear(scip, sepadata) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpConvexproj)
{  /*lint --e{715}*/
   SCIP_Real maxviolation;
   SCIP_SOL* lpsol;
   SCIP_SEPADATA* sepadata;

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* do not run if there is no interesting convex relaxation (with at least one nonlinear convex constraint),
    * or if we have found it to be numerically unstable
    * @todo: should it be with at least 2 nonlinear convex constraints?
    */
   if( sepadata->skipsepa )
   {
      SCIPdebugMsg(scip, "not running because convex relaxation is uninteresting or numerically unstable\n");
      return SCIP_OKAY;
   }

   /* the separator needs an NLP solver */
   if( SCIPgetNNlpis(scip) == 0 )
      return SCIP_OKAY;

   /* only call separator up to a maximum depth */
   if( sepadata->maxdepth >= 0 && depth > sepadata->maxdepth )
      return SCIP_OKAY;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* do not run if SCIP does not have constructed an NLP */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "NLP not constructed, skipping convex projection separator\n");
      return SCIP_OKAY;
   }

   /* recompute convex NLP relaxation if the variable set changed and we are still at the root node */
   if( sepadata->nlpiprob != NULL && SCIPgetNVars(scip) != sepadata->nlpinvars  && SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( sepadataClear(scip, sepadata) );
      assert(sepadata->nlpiprob == NULL);
   }

   /* create or update convex NLP relaxation */
   if( sepadata->nlpiprob == NULL )
   {
      /* store convex nonlinear constraints */
      SCIP_CALL( storeNonlinearConvexNlrows(scip, sepadata, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip)) );

      /* check that convex NLP relaxation is interesting (more than one nonlinear constraint) */
      if( sepadata->nnlrows < 1 )
      {
         SCIPdebugMsg(scip, "convex relaxation uninteresting, don't run\n");
         sepadata->skipsepa = TRUE;
         return SCIP_OKAY;
      }

      sepadata->nlpinvars = SCIPgetNVars(scip);
      sepadata->nlpi = SCIPgetNlpis(scip)[0];
      assert(sepadata->nlpi != NULL);

      SCIP_CALL( SCIPhashmapCreate(&sepadata->var2nlpiidx, SCIPblkmem(scip), sepadata->nlpinvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &sepadata->nlpivars, SCIPgetVars(scip), sepadata->nlpinvars) ); /*lint !e666*/

      SCIP_CALL( SCIPcreateNlpiProblemFromNlRows(scip, sepadata->nlpi, &sepadata->nlpiprob, "convexproj-nlp", SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip),
            sepadata->var2nlpiidx, NULL, NULL, SCIPgetCutoffbound(scip), FALSE, TRUE) );

      /* add rows of the LP
       * we do not sue the depth argument of the callback because we want to build a globally valid initia lrelaxation
       */
      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( SCIPaddNlpiProblemRows(scip, sepadata->nlpi, sepadata->nlpiprob, sepadata->var2nlpiidx,
                  SCIPgetLPRows(scip), SCIPgetNLPRows(scip)) );
      }

      /* set quadratic part of objective function */
      SCIP_CALL( setQuadraticObj(scip, sepadata) );
   }
   else
   {
      SCIP_CALL( SCIPupdateNlpiProblem(scip, sepadata->nlpi, sepadata->nlpiprob, sepadata->var2nlpiidx,
            sepadata->nlpivars, sepadata->nlpinvars, SCIPgetCutoffbound(scip)) );
   }

   /* assert that the lp solution satisfies the cutoff bound; if this fails then we shouldn't have a cutoff bound in the
    * nlpi, since then the projection could be in the interior of the actual convex relaxation */
   assert(SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL ||
         SCIPisLE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)));

   /* get current sol: LP or pseudo solution if LP sol is not available */
   SCIP_CALL( SCIPcreateCurrentSol(scip, &lpsol, NULL) );

   /* do not run if current solution's violation is small */
   SCIP_CALL( computeMaxViolation(scip, sepadata, lpsol, &maxviolation) );
   if( maxviolation < VIOLATIONFAC * SCIPfeastol(scip) )
   {
      SCIPdebugMsg(scip, "solution doesn't violate constraints enough, do not separate\n");
      SCIP_CALL( SCIPfreeSol(scip, &lpsol) );
      return SCIP_OKAY;
   }

   /* run the separator */
   *result = SCIP_DIDNOTFIND;

   /* separateCuts computes the projection and then gradient cuts on each constraint that was originally violated */
   SCIP_CALL( separateCuts(scip, sepa, lpsol, result) );

   /* free memory */
   SCIP_CALL( SCIPfreeSol(scip, &lpsol) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the convexproj separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaConvexproj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create convexproj separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   /* this sets all data in sepadata to 0 */
   BMSclearMemory(sepadata);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpConvexproj, NULL,
         sepadata) );
   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeConvexproj) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolConvexproj) );

   /* add convexproj separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxdepth",
         "maximal depth at which the separator is applied (-1: unlimited)",
         &sepadata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nlpiterlimit",
         "iteration limit of NLP solver; 0 for no limit",
         &sepadata->nlpiterlimit, TRUE, DEFAULT_NLPITERLIM, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
