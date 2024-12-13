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

/**@file   scip_nlpi.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for NLP interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 *
 * @todo check SCIP_STAGE_* switches
 * @todo allow for optional callbacks
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip_nlp.h"
#include "blockmemshell/memory.h"
#include "scip/scip_expr.h"
#include "scip/scip_lp.h"
#include "scip/scip_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_nlpi.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/pub_expr.h"
#include "scip/pub_lp.h"
#include "scip/pub_var.h"
#include "scip/expr_varidx.h"
#include "scip/debug.h"
#include "scip/nlpi.h"
#include "scip/paramset.h"
#include "scip/set.h"
#include "scip/struct_scip.h"


/** method to call, when the priority of an NLPI was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNlpiPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSetPriorityNlpi() to mark the nlpis unsorted */
   SCIP_CALL( SCIPsetNlpiPriority(scip, (SCIP_NLPI*)paramdata, SCIPparamGetInt(param)) );

   return SCIP_OKAY;
}

/** create varidx expression for var expression
 *
 * called when expr is duplicated for addition to NLPI
 */
static
SCIP_DECL_EXPR_MAPEXPR(mapvar2varidx)
{
   SCIP_HASHMAP* var2idx;
   int varidx;

   assert(sourcescip != NULL);
   assert(sourcescip == targetscip);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);
   assert(*targetexpr == NULL);
   assert(mapexprdata != NULL);

   /* do not provide map if not variable */
   if( !SCIPisExprVar(sourcescip, sourceexpr) )
      return SCIP_OKAY;

   assert(SCIPvarIsActive(SCIPgetVarExprVar(sourceexpr)));

   var2idx = (SCIP_HASHMAP*)mapexprdata;
   assert(SCIPhashmapExists(var2idx, SCIPgetVarExprVar(sourceexpr)));

   varidx = SCIPhashmapGetImageInt(var2idx, SCIPgetVarExprVar(sourceexpr));

   SCIP_CALL( SCIPcreateExprVaridx(targetscip, targetexpr, varidx, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** creates an NLPI and includes it into SCIP */
SCIP_RETCODE SCIPincludeNlpi(
   SCIP*                           scip,                        /**< SCIP data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying an NLPI, can be NULL */
   SCIP_DECL_NLPIFREE              ((*nlpifree)),               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer)),   /**< get solver pointer, can be NULL */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem)),      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem)),        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer)),  /**< get problem pointer, can be NULL */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars)),            /**< add variables */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints)),     /**< add constraints */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective)),       /**< set objective */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds)),       /**< change variable bounds */
   SCIP_DECL_NLPICHGCONSSIDES      ((*nlpichgconssides)),       /**< change constraint sides */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset)),          /**< delete a set of constraints */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset)),         /**< delete a set of constraints */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs)),     /**< change coefficients in linear part of a constraint or objective */
   SCIP_DECL_NLPICHGEXPR           ((*nlpichgexpr)),            /**< change nonlinear expression a constraint or objective */
   SCIP_DECL_NLPICHGOBJCONSTANT    ((*nlpichgobjconstant)),     /**< change the constant offset in the objective */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess)),    /**< set initial guess, can be NULL */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve)),              /**< solve NLP */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat)),         /**< get solution status */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat)),        /**< get termination status */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution)),        /**< get solution */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics)),      /**< get solve statistics */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
   )
{
   SCIP_NLPI* nlpi = NULL;
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeNlpi", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether NLPI of given name is already present */
   if( SCIPfindNlpi(scip, name) != NULL )
   {
      SCIPerrorMessage("NLPI <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPnlpiCreate(&nlpi, name, description, priority,
      nlpicopy, nlpifree, nlpigetsolverpointer,
      nlpicreateproblem, nlpifreeproblem, nlpigetproblempointer,
      nlpiaddvars, nlpiaddconstraints, nlpisetobjective, nlpichgvarbounds, nlpichgconssides, nlpidelvarset, nlpidelconsset, nlpichglinearcoefs, nlpichgexpr, nlpichgobjconstant,
      nlpisetinitialguess, nlpisolve, nlpigetsolstat, nlpigettermstat, nlpigetsolution, nlpigetstatistics,
      nlpidata) );
   assert(nlpi != NULL);

   SCIP_CALL( SCIPsetIncludeNlpi(scip->set, nlpi) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "nlpi/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of NLPI <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, paramname, paramdesc,
         NULL, FALSE, SCIPnlpiGetPriority(nlpi), INT_MIN/4, INT_MAX/4,
         paramChgdNlpiPriority, (SCIP_PARAMDATA*)nlpi) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_NLPI* SCIPfindNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of NLPI */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindNlpi(scip->set, name);
}

/** returns the array of currently available NLPIs (sorted by priority) */
SCIP_NLPI** SCIPgetNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortNlpis(scip->set);

   return scip->set->nlpis;
}

/** returns the number of currently available NLPIs */
int SCIPgetNNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nnlpis;
}

/** sets the priority of an NLPI */
SCIP_RETCODE SCIPsetNlpiPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of the NLPI */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSetPriorityNlpi(scip->set, nlpi, priority);

   return SCIP_OKAY;
}

/** gets internal pointer to NLP solver */
SCIP_DECL_NLPIGETSOLVERPOINTER(SCIPgetNlpiSolverPointer)
{
   assert(scip != NULL);

   return SCIPnlpiGetSolverPointer(scip->set, nlpi, problem);
}

/** creates an empty problem instance */
SCIP_DECL_NLPICREATEPROBLEM(SCIPcreateNlpiProblem)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiCreateProblem(scip->set, nlpi, problem, name) );

   return SCIP_OKAY;
}

/** frees a problem instance */
SCIP_DECL_NLPIFREEPROBLEM(SCIPfreeNlpiProblem)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiFreeProblem(scip->set, nlpi, problem) );

   return SCIP_OKAY;
}

/** gets internal pointer to solver-internal problem instance */
SCIP_DECL_NLPIGETPROBLEMPOINTER(SCIPgetNlpiProblemPointer)
{
   assert(scip != NULL);

   return SCIPnlpiGetProblemPointer(scip->set, nlpi, problem);
}

/** add variables to nlpi */
SCIP_DECL_NLPIADDVARS(SCIPaddNlpiVars)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiAddVars(scip->set, nlpi, problem, nvars, lbs, ubs, varnames) );

   return SCIP_OKAY;
}

/** add constraints to nlpi */
SCIP_DECL_NLPIADDCONSTRAINTS(SCIPaddNlpiConstraints)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiAddConstraints(scip->set, nlpi, problem, nconss, lhss, rhss, nlininds, lininds, linvals, exprs, names) );

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected */
SCIP_DECL_NLPISETOBJECTIVE(SCIPsetNlpiObjective)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetObjective(scip->set, nlpi, problem, nlins, lininds, linvals, expr, constant) );

   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_DECL_NLPICHGVARBOUNDS(SCIPchgNlpiVarBounds)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgVarBounds(scip->set, nlpi, problem, nvars, indices, lbs, ubs) );

   return SCIP_OKAY;
}

/** change constraint sides */
SCIP_DECL_NLPICHGCONSSIDES(SCIPchgNlpiConsSides)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgConsSides(scip->set, nlpi, problem, nconss, indices, lhss, rhss) );

   return SCIP_OKAY;
}

/** delete a set of variables */
SCIP_DECL_NLPIDELVARSET(SCIPdelNlpiVarSet)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiDelVarSet(scip->set, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** delete a set of constraints */
SCIP_DECL_NLPIDELCONSSET(SCIPdelNlpiConsSet)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiDelConsSet(scip->set, nlpi, problem, dstats, dstatssize) );

   return SCIP_OKAY;
}

/** changes or adds linear coefficients in a constraint or objective */
SCIP_DECL_NLPICHGLINEARCOEFS(SCIPchgNlpiLinearCoefs)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgLinearCoefs(scip->set, nlpi, problem, idx, nvals, varidxs, vals) );

   return SCIP_OKAY;
}

/** change the expression in the nonlinear part */
SCIP_DECL_NLPICHGEXPR(SCIPchgNlpiExpr)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgExpr(scip->set, nlpi, problem, idxcons, expr) );

   return SCIP_OKAY;
}

/** change the constant offset in the objective */
SCIP_DECL_NLPICHGOBJCONSTANT(SCIPchgNlpiObjConstant)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiChgObjConstant(scip->set, nlpi, problem, objconstant) );

   return SCIP_OKAY;
}

/** sets initial guess */
SCIP_DECL_NLPISETINITIALGUESS(SCIPsetNlpiInitialGuess)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSetInitialGuess(scip->set, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues) );

   return SCIP_OKAY;
}

/** try to solve NLP with all parameters given as SCIP_NLPPARAM struct
 *
 * Typical use is
 *
 *     SCIP_NLPPARAM nlparam = { SCIP_NLPPARAM_DEFAULT(scip); }
 *     nlpparam.iterlim = 42;
 *     SCIP_CALL( SCIPsolveNlpiParam(scip, nlpi, nlpiproblem, nlpparam) );
 *
 * or, in "one" line:
 *
 *     SCIP_CALL( SCIPsolveNlpiParam(scip, nlpi, nlpiproblem,
 *        (SCIP_NLPPARAM){ SCIP_NLPPARAM_DEFAULT(scip), .iterlimit = 42 }) );
 *
 * To get the latter, also \ref SCIPsolveNlpi can be used.
 */
SCIP_DECL_NLPISOLVE(SCIPsolveNlpiParam)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiSolve(scip->set, scip->stat, nlpi, problem, &param) );

   return SCIP_OKAY;
}

#if defined(_MSC_VER) && _MSC_VER < 1800
/* warn that SCIPsolveNlpi() macro isn't perfect with ancient MSVC */
#pragma message ( "Warning: designated initializers not supported by this version of MSVC. Parameters given to NLP solves will be ignored." )
#endif

/** gives solution status */
SCIP_DECL_NLPIGETSOLSTAT(SCIPgetNlpiSolstat)
{
   assert(scip != NULL);

   return SCIPnlpiGetSolstat(scip->set, nlpi, problem);
}

/** gives termination reason */
SCIP_DECL_NLPIGETTERMSTAT(SCIPgetNlpiTermstat)
{
   assert(scip != NULL);

   return SCIPnlpiGetTermstat(scip->set, nlpi, problem);
}

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_DECL_NLPIGETSOLUTION(SCIPgetNlpiSolution)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetSolution(scip->set, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues, objval) );

   return SCIP_OKAY;
}

/** gives solve statistics */
SCIP_DECL_NLPIGETSTATISTICS(SCIPgetNlpiStatistics)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPnlpiGetStatistics(scip->set, nlpi, problem, statistics) );

   return SCIP_OKAY;
}

/** creates a NLPI problem from given nonlinear rows
 *
 * The function computes for each variable the number of non-linear occurrences and stores it in the nlscore array.
 *
 * @note the first row corresponds always to the cutoff row (even if cutoffbound is SCIPinfinity(scip))
 **/
SCIP_RETCODE SCIPcreateNlpiProblemFromNlRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM**    nlpiprob,           /**< buffer to store pointer to created nlpi problem */
   const char*           name,               /**< name to give to problem */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< number of nonlinear rows */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_HASHMAP*         nlrow2idx,          /**< empty hash map to store mapping between variables and indices in nlpiprob, can be NULL */
   SCIP_Real*            nlscore,            /**< array to store the score of each nonlinear variable (NULL if not needed) */
   SCIP_Real             cutoffbound,        /**< cutoff bound */
   SCIP_Bool             setobj,             /**< whether the objective function should be set to one of the SCIP problem */
   SCIP_Bool             onlyconvex          /**< filter only for convex constraints */
   )
{
   SCIP_EXPR** exprs;
   SCIP_Real** linvals;
   int** lininds;
   int* nlininds;
   SCIP_Real* lhss;
   SCIP_Real* rhss;
   const char** names;
   SCIP_VAR** vars;
   int nvars;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real* objvals = NULL;
   int* objinds = NULL;
   const char** varnames;
   int nobjinds;
   int nconss;
   SCIP_EXPRITER* it = NULL;
   int i;

   assert(nlpiprob != NULL);
   assert(name != NULL);
   assert(var2idx != NULL);
   assert(nlrows != NULL);
   assert(nnlrows > 0);
   assert(nlpi != NULL);

   SCIPdebugMsg(scip, "SCIPcreateNlpiProblemFromNlRows() called with cutoffbound %g\n", cutoffbound);

   SCIP_CALL( SCIPnlpiCreateProblem(scip->set, nlpi, nlpiprob, name) );

   if( nlscore != NULL )
   {
      BMSclearMemoryArray(nlscore, SCIPgetNVars(scip));
   }
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nconss = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &names, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nnlrows + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nnlrows + 1) );

   if( setobj )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &objinds, nvars) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varnames, nvars) );

   /* create a unique mapping between variables and {0,..,nvars-1} */
   nobjinds = 0;
   for( i = 0; i < nvars; ++i )
   {
      assert(vars[i] != NULL);
      SCIP_CALL( SCIPhashmapInsertInt(var2idx, (void*)vars[i], i) );

      lbs[i] = SCIPvarGetLbLocal(vars[i]);
      ubs[i] = SCIPvarGetUbLocal(vars[i]);
      varnames[i] = SCIPvarGetName(vars[i]);

      /* collect non-zero objective coefficients */
      if( setobj && !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
      {
         assert(objvals != NULL);
         assert(objinds != NULL);

         objvals[nobjinds] = SCIPvarGetObj(vars[i]);
         objinds[nobjinds] = i;
         ++nobjinds;
      }
   }

   /* add variables */
   SCIP_CALL( SCIPaddNlpiVars(scip, nlpi, *nlpiprob, nvars, lbs, ubs, varnames) );
   SCIPfreeBufferArray(scip, &varnames);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   /* set the objective function */
   if( setobj )
   {
      if( nobjinds > 0 )
      {
         SCIP_CALL( SCIPsetNlpiObjective(scip, nlpi, *nlpiprob, nobjinds, objinds, objvals, NULL, 0.0) );
      }

      SCIPfreeBufferArray(scip, &objinds);
      SCIPfreeBufferArray(scip, &objvals);
   }

   /* add row for cutoff bound even if cutoffbound == SCIPinfinity() */
   lhss[nconss] = -SCIPinfinity(scip);
   rhss[nconss] = cutoffbound;
   names[nconss] = "objcutoff";
   lininds[nconss] = NULL;
   linvals[nconss] = NULL;
   nlininds[nconss] = 0;
   exprs[nconss] = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &lininds[nconss], nvars) ); /*lint !e866*/
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals[nconss], nvars) ); /*lint !e866*/

   for( i = 0; i < nvars; ++i )
   {
      if( !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
      {
         linvals[nconss][nlininds[nconss]] = SCIPvarGetObj(vars[i]);
         lininds[nconss][nlininds[nconss]] = i;
         ++nlininds[nconss];
      }
   }
   ++nconss;

   if( nlscore != NULL )
   {
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   }

   /* add convex nonlinear rows to NLPI problem */
   for( i = 0; i < nnlrows; ++i )
   {
      SCIP_Bool userhs;
      SCIP_Bool uselhs;
      int k;
      SCIP_NLROW* nlrow;

      nlrow = nlrows[i];
      assert(nlrow != NULL);

      uselhs = FALSE;
      userhs = FALSE;

      /* check curvature together with constraint sides of a nonlinear row */
      if( SCIPnlrowGetExpr(nlrow) == NULL )
      {
         uselhs = TRUE;
         userhs = TRUE;
      }
      else
      {
         if( (!onlyconvex || SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONVEX)
            && !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)) )
            userhs = TRUE;
         if( (!onlyconvex || SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONCAVE)
            && !SCIPisInfinity(scip, SCIPnlrowGetLhs(nlrow)) )
            uselhs = TRUE;
      }

      if( !uselhs && !userhs )
         continue;

      lhss[nconss] = uselhs ? SCIPnlrowGetLhs(nlrow) - SCIPnlrowGetConstant(nlrow) : -SCIPinfinity(scip);
      rhss[nconss] = userhs ? SCIPnlrowGetRhs(nlrow) - SCIPnlrowGetConstant(nlrow) :  SCIPinfinity(scip);
      names[nconss] = SCIPnlrowGetName(nlrow);
      nlininds[nconss] = 0;
      lininds[nconss] = NULL;
      linvals[nconss] = NULL;

      /* copy linear part */
      if( SCIPnlrowGetNLinearVars(nlrow) > 0 )
      {
         SCIP_VAR* var;

         nlininds[nconss] = SCIPnlrowGetNLinearVars(nlrow);

         SCIP_CALL( SCIPallocBufferArray(scip, &lininds[nconss], nlininds[nconss]) ); /*lint !e866*/
         SCIP_CALL( SCIPallocBufferArray(scip, &linvals[nconss], nlininds[nconss]) ); /*lint !e866*/

         for( k = 0; k < nlininds[nconss]; ++k )
         {
            var = SCIPnlrowGetLinearVars(nlrow)[k];
            assert(var != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var));

            lininds[nconss][k] = SCIPhashmapGetImageInt(var2idx, (void*)var);
            assert(var == vars[lininds[nconss][k]]);
            linvals[nconss][k] = SCIPnlrowGetLinearCoefs(nlrow)[k];
         }
      }

      if( SCIPnlrowGetExpr(nlrow) != NULL )
      {
         /* create copy of expr that uses varidx expressions corresponding to variables indices in NLPI */
         SCIP_CALL( SCIPduplicateExpr(scip, SCIPnlrowGetExpr(nlrow), &exprs[nconss], mapvar2varidx, var2idx, NULL, NULL) );
      }
      else
      {
         exprs[nconss] = NULL;
      }

      /* update nlscore */
      if( nlscore != NULL && exprs[nconss] != NULL )
      {
         SCIP_EXPR* expr;
         int varidx;

         SCIP_CALL( SCIPexpriterInit(it, exprs[nconss], SCIP_EXPRITER_DFS, FALSE) );
         for( expr = exprs[nconss]; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )  /*lint !e441*/ /*lint !e440*/
         {
            if( !SCIPisExprVaridx(scip, expr) )
               continue;

            varidx = SCIPgetIndexExprVaridx(expr);
            assert(varidx >= 0);
            assert(varidx < nvars);

            /* update nlscore */
            nlscore[varidx] += 1.0;
         }
      }

      /* if the row to index hash map is provided, we need to store the row index */
      if( nlrow2idx != NULL )
      {
         SCIP_CALL( SCIPhashmapInsertInt(nlrow2idx, nlrow, nconss) );
      }

      ++nconss;
   }
   assert(nconss > 0);

   /* pass all constraint information to nlpi */
   SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpi, *nlpiprob, nconss, lhss, rhss, nlininds, lininds, linvals,
         exprs, names) );

   if( it != NULL )
   {
      SCIPfreeExpriter(&it);
   }

   /* free memory */
   for( i = nconss - 1; i > 0; --i )
   {
      if( nlininds[i] > 0 )
      {
         assert(linvals[i] != NULL);
         assert(lininds[i] != NULL);
         SCIPfreeBufferArray(scip, &linvals[i]);
         SCIPfreeBufferArray(scip, &lininds[i]);
      }
      if( exprs[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &exprs[i]) );
      }
   }
   /* free row for cutoff bound even if objective is 0 */
   SCIPfreeBufferArray(scip, &linvals[i]);
   SCIPfreeBufferArray(scip, &lininds[i]);

   SCIPfreeBufferArray(scip, &rhss);
   SCIPfreeBufferArray(scip, &lhss);
   SCIPfreeBufferArray(scip, &names);
   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &exprs);

   return SCIP_OKAY;
}

/** updates variable bounds and the cutoff row in a NLPI problem
 *
 * The NLPI problem must have been setup by SCIPcreateNlpiProblemFromNlRows().
 */
SCIP_RETCODE SCIPupdateNlpiProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_HASHMAP*         var2nlpiidx,        /**< mapping between variables and nlpi indices */
   SCIP_VAR**            nlpivars,           /**< array containing all variables of the nlpi */
   int                   nlpinvars,          /**< total number of nlpi variables */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   )
{
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int* inds;
   int i;

   SCIPdebugMsg(scip, "SCIPupdateNlpiProblem() called\n");

   /* update variable bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nlpinvars) );

   for( i = 0; i < nlpinvars; ++i )
   {
      assert(nlpivars[i] != NULL);
      assert(SCIPhashmapExists(var2nlpiidx, (void*)nlpivars[i]));

      lbs[i] = SCIPvarGetLbLocal(nlpivars[i]);
      ubs[i] = SCIPvarGetUbLocal(nlpivars[i]);
      inds[i] = SCIPhashmapGetImageInt(var2nlpiidx, (void*)nlpivars[i]);
      assert(inds[i] >= 0 && inds[i] < nlpinvars);
   }

   SCIP_CALL( SCIPchgNlpiVarBounds(scip, nlpi, nlpiprob, nlpinvars, inds, lbs, ubs) );

   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ubs);
   SCIPfreeBufferArray(scip, &lbs);

   /* update cutoff row */
   lhs = -SCIPinfinity(scip);
   rhs = cutoffbound;
   i = 0;

   SCIP_CALL( SCIPchgNlpiConsSides(scip, nlpi, nlpiprob, 1, &i, &lhs, &rhs) );

   return SCIP_OKAY;
}

/** adds SCIP_ROWs to a NLPI problem */
SCIP_RETCODE SCIPaddNlpiProblemRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_ROW**            rows,               /**< rows to add */
   int                   nrows               /**< number of rows to add */
   )
{
   const char** names;
   SCIP_Real* lhss;
   SCIP_Real* rhss;
   SCIP_Real** linvals;
   int** lininds;
   int* nlininds;
   int i;

   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(var2idx != NULL);
   assert(nrows == 0 || rows != NULL);

   SCIPdebugMsg(scip, "SCIPaddNlpiProblemRows() called with %d rows\n", nrows);

   if( nrows <= 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &names, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nrows) );

   for( i = 0; i < nrows; ++i )
   {
      int k;

      assert(rows[i] != NULL);
      assert(SCIProwGetNNonz(rows[i]) <= SCIPgetNVars(scip));

      names[i] = SCIProwGetName(rows[i]);
      lhss[i] = SCIProwGetLhs(rows[i]) - SCIProwGetConstant(rows[i]);
      rhss[i] = SCIProwGetRhs(rows[i]) - SCIProwGetConstant(rows[i]);
      nlininds[i] = SCIProwGetNNonz(rows[i]);
      linvals[i] = SCIProwGetVals(rows[i]);
      lininds[i] = NULL;

      SCIP_CALL( SCIPallocBufferArray(scip, &lininds[i], SCIProwGetNNonz(rows[i])) ); /*lint !e866*/

      for( k = 0; k < SCIProwGetNNonz(rows[i]); ++k )
      {
         SCIP_VAR* var;

         var = SCIPcolGetVar(SCIProwGetCols(rows[i])[k]);
         assert(var != NULL);
         assert(SCIPhashmapExists(var2idx, (void*)var));

         lininds[i][k] = SCIPhashmapGetImageInt(var2idx, (void*)var);
         assert(lininds[i][k] >= 0 && lininds[i][k] < SCIPgetNVars(scip));
      }
   }

   /* pass all linear rows to the nlpi */
   SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpi, nlpiprob, nrows, lhss, rhss, nlininds, lininds, linvals,
         NULL, names) );

   /* free memory */
   for( i = nrows - 1; i >= 0; --i )
   {
      SCIPfreeBufferArray(scip, &lininds[i]);
   }
   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &rhss);
   SCIPfreeBufferArray(scip, &lhss);
   SCIPfreeBufferArray(scip, &names);

   return SCIP_OKAY;
}

/** adds SCIP_NLROWs to a NLPI problem */
SCIP_RETCODE SCIPaddNlpiProblemNlRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_NLROW**          nlrows,             /**< rows to add */
   int                   nnlrows             /**< number of rows to add */
   )
{
   const char** names;
   SCIP_Real* lhss;
   SCIP_Real* rhss;
   SCIP_Real** linvals;
   int** lininds;
   int* nlininds;
   SCIP_EXPR** exprs;
   int i;

   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(var2idx != NULL);
   assert(nnlrows == 0 || nlrows != NULL);

   SCIPdebugMsg(scip, "SCIPaddNlpiProblemNlRows() called with %d rows\n", nnlrows);

   if( nnlrows <= 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &names, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhss, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlininds, nnlrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nnlrows) );

   for( i = 0; i < nnlrows; ++i )
   {
      SCIP_NLROW* nlrow;

      nlrow = nlrows[i];
      assert(nlrow != NULL);

      lhss[i] = !SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)) ? SCIPnlrowGetLhs(nlrow) - SCIPnlrowGetConstant(nlrow) : -SCIPinfinity(scip);
      rhss[i] = !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)) ? SCIPnlrowGetRhs(nlrow) - SCIPnlrowGetConstant(nlrow) :  SCIPinfinity(scip);
      names[i] = SCIPnlrowGetName(nlrow);
      nlininds[i] = 0;
      lininds[i] = NULL;
      linvals[i] = NULL;

      /* copy linear part */
      if( SCIPnlrowGetNLinearVars(nlrow) > 0 )
      {
         SCIP_VAR* var;
         int k;

         nlininds[i] = SCIPnlrowGetNLinearVars(nlrow);

         SCIP_CALL( SCIPallocBufferArray(scip, &lininds[i], nlininds[i]) ); /*lint !e866*/
         SCIP_CALL( SCIPallocBufferArray(scip, &linvals[i], nlininds[i]) ); /*lint !e866*/

         for( k = 0; k < nlininds[i]; ++k )
         {
            var = SCIPnlrowGetLinearVars(nlrow)[k];
            assert(var != NULL);
            assert(SCIPhashmapExists(var2idx, (void*)var));

            lininds[i][k] = SCIPhashmapGetImageInt(var2idx, (void*)var);
            linvals[i][k] = SCIPnlrowGetLinearCoefs(nlrow)[k];
         }
      }

      if( SCIPnlrowGetExpr(nlrow) != NULL )
      {
         /* create copy of expr that uses varidx expressions corresponding to variables indices in NLPI */
         SCIP_CALL( SCIPduplicateExpr(scip, SCIPnlrowGetExpr(nlrow), &exprs[i], mapvar2varidx, var2idx, NULL, NULL) );
      }
      else
      {
         exprs[i] = NULL;
      }
   }

   /* pass all rows to the nlpi */
   SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpi, nlpiprob, nnlrows, lhss, rhss, nlininds, lininds, linvals, exprs, names) );

   /* free memory */
   for( i = nnlrows - 1; i >= 0; --i )
   {
      SCIPfreeBufferArrayNull(scip, &linvals[i]);
      SCIPfreeBufferArrayNull(scip, &lininds[i]);
      if( exprs[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &exprs[i]) );
      }
   }
   SCIPfreeBufferArray(scip, &exprs);
   SCIPfreeBufferArray(scip, &nlininds);
   SCIPfreeBufferArray(scip, &lininds);
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &rhss);
   SCIPfreeBufferArray(scip, &lhss);
   SCIPfreeBufferArray(scip, &names);

   return SCIP_OKAY;
}
