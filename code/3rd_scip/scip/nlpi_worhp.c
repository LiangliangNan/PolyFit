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

/**@file    nlpi_worhp.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   Worhp NLP interface
 * @author  Benjamin Mueller
 * @author  Renke Kuhlmann
 *
 * @todo So far, Worhp can not handle the case that variables have been fixed before warm-starting. Remove the code in
 * nlpiChgVarBoundsWorhp when this has changed.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_worhp.h"
#include "scip/nlpioracle.h"
#include "scip/exprinterpret.h"
#include "scip/interrupt.h"
#include "scip/scip_nlpi.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_solve.h"
#include "scip/pub_misc.h"

#include <stdio.h>
#include <stdlib.h>

#include "worhp/worhp.h"

#if WORHP_MAJOR < 2 && WORHP_MINOR < 10
#error "Require at least Worhp 1.10"
#endif

#define NLPI_DESC              "Worhp interface"            /**< description of solver */
#define NLPI_PRIORITY_IP       0                            /**< priority of NLP solver (Interior Point) */
#define NLPI_PRIORITY_SQP      -2000                        /**< priority of NLP solver (SQP) */

#define DEFAULT_VERBLEVEL      0                            /**< default verbosity level (0: normal 1: full 2: debug >2: more debug) */
#define DEFAULT_SCALEDKKT      TRUE                         /**< default whether KKT conditions are allowed to be scaled in the solver */
#define DEFAULT_RANDSEED       107                          /**< initial random seed */

#define MAXPERTURB             0.01                         /**< maximal perturbation of bounds in starting point heuristic */

/*
 * Data structures
 */

struct SCIP_NlpiData
{
   SCIP_Bool                   useip;        /**< should the Interior Point solver of Worhp be used? */
};

struct SCIP_NlpiProblem
{
   SCIP_NLPIORACLE*            oracle;       /**< Oracle-helper to store and evaluate NLP */
   SCIP_RANDNUMGEN*            randnumgen;   /**< random number generator */

   SCIP_NLPTERMSTAT            lasttermstat; /**< termination status from last run */
   SCIP_NLPSOLSTAT             lastsolstat;  /**< solution status from last run */
   SCIP_Real                   lasttime;     /**< time spend in last run */
   int                         lastniter;    /**< number of iterations in last run */

   SCIP_Real*                  lastprimal;   /**< primal solution from last run, if available */
   SCIP_Real*                  lastdualcons; /**< dual solution from last run, if available */
   SCIP_Real*                  lastduallb;   /**< dual solution for lower bounds from last run, if available */
   SCIP_Real*                  lastdualub;   /**< dual solution for upper bounds from last run, if available */
   int                         lastprimalsize; /**< size of lastprimal array */
   int                         lastdualconssize; /**< size of lastdualcons array */
   int                         lastduallbsize; /**< size of lastduallb array */
   int                         lastdualubsize; /**< size of lastdualub array */

   SCIP_Bool                   firstrun;     /**< whether the next NLP solve will be the first one (with the current problem structure) */
   SCIP_Real*                  initguess;    /**< initial values for primal variables, or NULL if not known */

   /* Worhp data structures */
   OptVar*                     opt;          /**< Worhp variables  */
   Workspace*                  wsp;          /**< Worhp working space */
   Params*                     par;          /**< Worhp parameters */
   Control*                    cnt;          /**< Worhp control */
};

/*
 * Local methods
 */

/** clears the last solution information */
static
void invalidateSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(problem->lastprimal),   problem->lastprimalsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(problem->lastdualcons), problem->lastdualconssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(problem->lastduallb),   problem->lastduallbsize);
   SCIPfreeBlockMemoryArrayNull(scip, &(problem->lastdualub),   problem->lastdualubsize);

   problem->lastprimalsize = 0;
   problem->lastdualconssize = 0;
   problem->lastduallbsize = 0;
   problem->lastdualubsize = 0;
   problem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
}

/** evaluate last Worhp run */
static
SCIP_RETCODE evaluateWorhpRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   int i;

   assert(problem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->cnt != NULL);

   switch( problem->cnt->status )
   {
   case InitError:
   {
      /* initialization error */
      SCIPdebugMsg(scip, "Worhp failed because of initialization error!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OUTOFMEMORY;
      break;
   }

   case DataError:
   {
      /* data error */
      SCIPdebugMsg(scip, "Worhp failed because of data error!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
      break;
   }

   case LicenseError:
   {
      /* license error */
      SCIPerrorMessage("Worhp failed because of license error!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_LICENSEERROR;
      break;
   }

   case evalsNaN:
   {
      /* evaluation errors */
      SCIPdebugMsg(scip, "Worhp failed because of a NaN value in an evaluation!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERROR;
      break;
   }

   case QPerror:
   case MinimumStepsize:
   case TooBig:
   case LinearSolverFailed:
   {
      /* numerical errors during solution of NLP */
      SCIPdebugMsg(scip, "Worhp failed because of a numerical error during optimization!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERICERROR;
      break;
   }

   case MaxCalls:
   case MaxIter:
   {
      /* maximal number of calls or iteration */
      SCIPdebugMsg(scip, "Worhp failed because maximal number of calls or iterations is reached!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_ITERLIMIT;
      break;
   }

   case Timeout:
   {
      /* time limit reached */
      SCIPdebugMsg(scip, "Worhp failed because time limit is reached!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_TIMELIMIT;
      break;
   }

   case DivergingPrimal:
   case DivergingDual:
   {
      /* iterates diverge */
      SCIPdebugMsg(scip, "Worhp failed because of diverging iterates!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNBOUNDED;
      problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERICERROR;
      break;
   }

   case LocalInfeas:
   case LocalInfeasOptimal:
   {
      /* infeasible stationary point found */
      SCIPdebugMsg(scip, "Worhp failed because of convergence against infeasible stationary point!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case GlobalInfeas:
   {
      /* infeasible stationary point found */
      SCIPdebugMsg(scip, "Worhp failed because of convergence against infeasible stationary point!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_GLOBINFEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case RegularizationFailed:
   {
      /* regularization of Hessian matrix failed */
      SCIPdebugMsg(scip, "Worhp failed because of regularization of Hessian matrix failed!\n");
      invalidateSolution(scip, problem);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERICERROR;
      break;
   }

   case OptimalSolution:
   {
      /* everything went fine */
      SCIPdebugMsg(scip, "Worhp terminated successfully at a local optimum!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_LOCOPT;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case OptimalSolutionConstantF:
   {
      /* feasible point, KKT conditions are not satisfied, and Worhp thinks that there is no objective function */
      SCIPdebugMsg(scip, "Worhp terminated successfully with a feasible point but KKT are not met!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case AcceptableSolutionSKKT:
   case AcceptableSolutionScaled:
   case AcceptablePreviousScaled:
   {
      /* feasible point but KKT conditions are violated in unscaled space */
      SCIPdebugMsg(scip, "Worhp terminated successfully with a feasible point but KKT are violated in unscaled space!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case LowPassFilterOptimal:
   {
      /* feasible and no further progress */
      SCIPdebugMsg(scip, "Worhp terminated at feasible solution without further progress!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case FeasibleSolution:
   {
      /* feasible and in feasibility mode, i.e., optimality not required */
      SCIPdebugMsg(scip, "Worhp terminated at feasible solution, optimality was not required!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case AcceptableSolution:
   case AcceptableSolutionConstantF:
   {
      /* acceptable solution found, but stopped due to limit or error */
      SCIPdebugMsg(scip, "Worhp terminated at acceptable solution due to limit or error!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case AcceptablePrevious:
   case AcceptablePreviousConstantF:
   {
      /* previously acceptable solution was found, but stopped due to limit or error */
      SCIPdebugMsg(scip, "Worhp previously found acceptable solution but terminated due to limit or error!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case LowPassFilterAcceptable:
   {
      /* acceptable solution found, and no further progress */
      SCIPdebugMsg(scip, "Worhp found acceptable solution but terminated due no further progress!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case SearchDirectionZero:
   case SearchDirectionSmall:
   {
      /* acceptable solution found, but search direction is small or zero */
      SCIPdebugMsg(scip, "Worhp found acceptable solution but search direction is small or zero!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   case FritzJohn:
   case NotDiffable:
   case Unbounded:
   {
      /* acceptable solution found, but not optimal */
      SCIPdebugMsg(scip, "Worhp found acceptable solution but terminated perhaps due to nondifferentiability, unboundedness or at Fritz John point!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
      break;
   }

   default:
   {
      SCIPerrorMessage("Worhp returned with unknown solution status %d\n", problem->cnt->status);
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
      return SCIP_OKAY;
   }
   }

  /* store solution */
  if( problem->lastprimal == NULL )
  {
     if( problem->opt->m > 0 )
     {
        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &problem->lastdualcons, problem->opt->Mu, problem->opt->m) );
        problem->lastdualconssize = problem->opt->m;
     }

     if( problem->opt->n > 0 )
     {
        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &problem->lastprimal, problem->opt->X, problem->opt->n) );
        problem->lastprimalsize = problem->opt->n;

        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->lastduallb, problem->opt->n) );
        problem->lastduallbsize = problem->opt->n;

        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &problem->lastdualub, problem->opt->n) );
        problem->lastdualubsize = problem->opt->n;
     }
  }
  else
  {
     BMScopyMemoryArray(problem->lastprimal, problem->opt->X, problem->opt->n);
     BMScopyMemoryArray(problem->lastdualcons, problem->opt->Mu, problem->opt->m);
  }

  for( i = 0; i < problem->opt->n; ++i )
  {
     if( problem->opt->Lambda[i] <= 0.0 )
     {
        problem->lastduallb[i] = -problem->opt->Lambda[i];
        problem->lastdualub[i] = 0.0;
     }
     else
     {
        problem->lastduallb[i] = 0.0;
        problem->lastdualub[i] = problem->opt->Lambda[i];
     }
  }

  assert(problem->lastprimal != NULL || problem->opt->n == 0);
  assert(problem->lastdualcons != NULL || problem->opt->m == 0);
  assert(problem->lastduallb != NULL || problem->opt->n == 0);
  assert(problem->lastdualub != NULL || problem->opt->n == 0);

  return SCIP_OKAY;
}

/** evaluates objective function and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userF(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_Real objval;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPnlpiOracleEvalObjectiveValue(scip, problem->oracle, problem->opt->X, &objval) );
   problem->opt->F = problem->wsp->ScaleObj * objval;

#ifdef SCIP_DEBUG_USERF
   {
      int i;

      printf("userF()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);
      printf("  obj = %g\n", problem->opt->F);
   }
#endif

   return SCIP_OKAY;
}

/** evaluates constraints and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPnlpiOracleEvalConstraintValues(scip, problem->oracle, problem->opt->X, problem->opt->G) );

#ifdef SCIP_DEBUG_USERG
   {
      int i;

      printf("userG()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);

      for( i = 0; i < problem->opt->m; ++i )
         printf("  cons[%d] = %g\n", i, problem->opt->G[i]);
   }
#endif

   return SCIP_OKAY;
}

/** computes objective gradient and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userDF(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_Real objval;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   /* TODO this needs to be changed if we store the gradient of the objective function in a sparse format */
   SCIP_CALL( SCIPnlpiOracleEvalObjectiveGradient(scip, problem->oracle, problem->opt->X, TRUE, &objval,
         problem->wsp->DF.val) );

   /* scale gradient if necessary */
   if( problem->wsp->ScaleObj != 1.0 )
   {
      int i;
      for( i = 0; i < problem->opt->n; ++i )
         problem->wsp->DF.val[i] *= problem->wsp->ScaleObj;
   }

#ifdef SCIP_DEBUG_USERDF
   {
      int i;

      printf("userDF()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);

      for( i = 0; i < problem->opt->n; ++i )
         printf("  DF[%d] = %g\n", i, problem->wsp->DF.val[i]);
   }
#endif

   return SCIP_OKAY;
}

/** computes jacobian matrix and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userDG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real* jacvals;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &jacvals, problem->wsp->DG.nnz) );
   retcode = SCIPnlpiOracleEvalJacobian(scip, problem->oracle, problem->opt->X, TRUE, NULL, jacvals);

   if( retcode == SCIP_OKAY )
   {
      int i;

      /* map values with DG indices */
      for( i = 0; i < problem->wsp->DG.nnz; ++i )
      {
         problem->wsp->DG.val[i] = jacvals[ problem->wsp->DG.perm[i]-1 ];
      }

#ifdef SCIP_DEBUG_USERDG
      printf("userDG()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);
      for( i = 0; i < problem->wsp->DG.nnz; ++i )
         printf("  DG[%d] = %g\n", i, problem->wsp->DG.val[i]);
#endif
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &jacvals, problem->wsp->DG.nnz);

   return retcode;
}

/** computes hessian matrix and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userHM(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   const int* offset;
   SCIP_Real* hessianvals;
   SCIP_RETCODE retcode;
   int nnonz;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   /* get nonzero entries in HM of SCIP (excludes unused diagonal entries) */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(scip, problem->oracle, &offset, NULL) );
   nnonz = offset[problem->opt->n];

   /* evaluate hessian */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hessianvals, problem->wsp->HM.nnz) );
   retcode = SCIPnlpiOracleEvalHessianLag(scip, problem->oracle, problem->opt->X, TRUE, TRUE, problem->wsp->ScaleObj,
         problem->opt->Mu, hessianvals);

   if( retcode == SCIP_OKAY )
   {
      int i;

      assert(problem->wsp->HM.nnz >= nnonz);
      for( i = 0; i < problem->wsp->HM.nnz; ++i )
      {
         /* an entry i with HM.perm[i] - 1 >= nnonz corresponds to an in SCIP non-existing diagonal element */
         if( problem->wsp->HM.perm[i] - 1 >= nnonz )
            problem->wsp->HM.val[i] = 0.0;
         else
            problem->wsp->HM.val[i] = hessianvals[ problem->wsp->HM.perm[i] - 1 ];
      }

#ifdef SCIP_DEBUG_HM
      printf("userHM()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);
      for( i = 0; i < problem->wsp->HM.nnz; ++i )
         printf("  HM[%d] = %g\n", i, problem->wsp->HM.val[i]);
#endif
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &hessianvals, problem->wsp->HM.nnz);

   return retcode;
}

/** Worhp print callback function that does nothing */ /*lint -e{715}*/
static void noprint(
   int                   mode,               /**< the mode */
   const char            s[]                 /**< a string */
   )
{ /*lint --e{715}*/
}

/** initialize Worhp data */
static
SCIP_RETCODE initWorhp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   Workspace* wsp;
   Control* cnt;
   OptVar* opt;
   Params* par;
   const SCIP_Real* lbs;
   const SCIP_Real* ubs;
   const int* offset;
   const int* cols;
   int i;
   int j;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->cnt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->opt != NULL);
   assert(problem->firstrun);

   wsp = problem->wsp;
   cnt = problem->cnt;
   opt = problem->opt;
   par = problem->par;

   /* properly zeros everything */
   WorhpPreInit(opt, wsp, par, cnt);

   /* set problem dimensions */
   opt->n = SCIPnlpiOracleGetNVars(problem->oracle);
   opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle);
   SCIPdebugMsg(scip, "nvars %d nconss %d\n", opt->n, opt->m);

   /* assume that objective function is dense; TODO use sparse representation */
   wsp->DF.nnz = opt->n;

   /* get number of non-zero entries in Jacobian */
   SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(scip, problem->oracle, &offset, NULL) );
   wsp->DG.nnz = offset[opt->m];
   SCIPdebugMsg(scip, "nnonz jacobian %d\n", wsp->DG.nnz);

   /* get number of non-zero entries in hessian
    *
    * note that Worhp wants to have the full diagonal in ANY case
    */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(scip, problem->oracle, &offset, &cols) );
   wsp->HM.nnz = 0;

   j = offset[0];
   for( i = 0; i < opt->n; ++i )
   {
      /* diagonal element */
      ++(wsp->HM.nnz);

      /* strict lower triangle elements */
      for( ; j < offset[i+1]; ++j )
      {
         if( i != cols[j] )
         {
            assert(i > cols[j]);
            ++(wsp->HM.nnz);
         }
      }
   }
   assert(offset[opt->n] <= wsp->HM.nnz);
   SCIPdebugMsg(scip, "nnonz hessian %d\n", wsp->HM.nnz);

   /* initialize data in Worhp */
   WorhpInit(opt, wsp, par, cnt);
   if (cnt->status != FirstCall)
   {
      SCIPerrorMessage("Initialisation failed.\n");
      return SCIP_ERROR;
   }

   /* set variable bounds */
   lbs = SCIPnlpiOracleGetVarLbs(problem->oracle);
   ubs = SCIPnlpiOracleGetVarUbs(problem->oracle);

   BMScopyMemoryArray(opt->XL, lbs, opt->n);
   BMScopyMemoryArray(opt->XU, ubs, opt->n);

#ifdef SCIP_DEBUG
   for( i = 0; i < opt->n; ++i )
   {
      SCIPdebugMsg(scip, "bounds %d [%g,%g]\n", i, opt->XL[i], opt->XU[i]);
   }
#endif

   /* set constraint sides */
   for( i = 0; i < opt->m; ++i )
   {
      opt->GL[i] = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
      opt->GU[i] = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);

      /* adjust constraint sides when both are infinite */
      if( SCIPisInfinity(scip, -opt->GL[i]) && SCIPisInfinity(scip, opt->GU[i]) )
      {
         SCIPwarningMessage(scip, "Lhs and rhs of constraint %d are infinite.\n", i);
         opt->GL[i] = -SCIPinfinity(scip) / 10.0;
         opt->GU[i] = SCIPinfinity(scip) / 10.0;
      }

      SCIPdebugMsg(scip, "sides %d [%g,%g]\n", i, opt->GL[i], opt->GU[i]);
   }

   /* set column indices of objective function; note that indices go from 1 to n */
   /* if( wsp->DF.NeedStructure ) evaluates to FALSE if DF is dense */
   {
      SCIPdebugPrintf("column indices of objective function:");
      for( i = 0; i < opt->n; ++i )
      {
         wsp->DF.row[i] = i + 1;
         SCIPdebugPrintf(" %d", wsp->DF.row[i]);
      }
      SCIPdebugPrintf("\n");
   }

   /* set column and row indices of non-zero entries in Jacobian matrix */
   /* if( wsp->DG.NeedStructure ) evaluates to FALSE if DG is dense */
   {
      int nnonz;

      SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(scip, problem->oracle, &offset, &cols) );
      assert(offset[opt->m] == wsp->DG.nnz);

      nnonz = 0;
      j = offset[0];
      for( i = 0; i < opt->m; ++i )
      {
         for( ; j < offset[i+1]; ++j )
         {
            wsp->DG.row[nnonz] = i + 1;
            wsp->DG.col[nnonz] = cols[j] + 1;
            ++nnonz;
         }
      }
      assert(nnonz == wsp->DG.nnz);

      /* sort arrays w.r.t the column-major order */
      SortWorhpMatrix(&wsp->DG);
   }

   /* set column and row indices of non-zero entries in hessian matrix */
   if( problem->par->UserHM || problem->par->FidifHM || problem->par->BFGSmethod > 1 )
   {
      int nnonz;
      int k;

      SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(scip, problem->oracle, &offset, &cols) );
      assert(offset[opt->n] <= wsp->HM.nnz);

      k = offset[opt->n];
      nnonz = 0;
      j = offset[0];
      for( i = 0; i < opt->n; ++i )
      {
         SCIP_Bool adddiag = TRUE;

         for( ; j < offset[i+1]; ++j )
         {
            problem->wsp->HM.row[nnonz] = i + 1;
            problem->wsp->HM.col[nnonz] = cols[j] + 1;
            ++nnonz;

            if( i == cols[j] )
               adddiag = FALSE;
         }

         /* Worhp wants to have each diagonal element */
         if( adddiag )
         {
            problem->wsp->HM.row[k] = i + 1;
            problem->wsp->HM.col[k] = i + 1;
            ++k;
         }
      }
      assert(nnonz == offset[opt->n]);
      assert(k == wsp->HM.nnz);

      /* sort arrays w.r.t the LT column-major order */
      SortWorhpMatrix(&wsp->HM);

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "column and row indices of hessian:\n");
      for( i = 0; i < wsp->HM.nnz; ++i )
      {
         SCIPdebugMsg(scip, "  entry %d: (row,col) = (%d,%d)\n", i, wsp->HM.row[i], wsp->HM.col[i]);
      }
#endif
   }

   return SCIP_OKAY;
}

/** update Worhp data */
static
SCIP_RETCODE updateWorhp(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   const SCIP_Real* lbs;
   const SCIP_Real* ubs;
   int i;

   assert(problem != NULL);
   assert(problem->cnt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->opt != NULL);
   assert(problem->oracle != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   WorhpRestart(problem->opt, problem->wsp, problem->par, problem->cnt);

   /* update variable bounds */
   lbs = SCIPnlpiOracleGetVarLbs(problem->oracle);
   ubs = SCIPnlpiOracleGetVarUbs(problem->oracle);
   for( i = 0; i < problem->opt->n; ++i )
   {
      problem->opt->XL[i] = lbs[i];
      problem->opt->XU[i] = ubs[i];
   }

   /* update constraint sides */
   for( i = 0; i < problem->opt->m; ++i )
   {
      problem->opt->GL[i] = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
      problem->opt->GU[i] = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);
   }

   return SCIP_OKAY;
}

/** frees Worhp data */
static
SCIP_RETCODE freeWorhp(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->cnt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->opt != NULL);

   if( problem->opt->initialised )
      WorhpFree(problem->opt, problem->wsp, problem->par, problem->cnt);

   return SCIP_OKAY;
}

/** pass NLP solve parameters to Ipopt */
static
SCIP_RETCODE handleNlpParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< Worhp interface */
   Params*               par,                /**< Worhp parameters */
   const SCIP_NLPPARAM   nlpparam            /**< NLP solve parameters */
   )
{
   SCIP_NLPIDATA* nlpidata;

   assert(par != NULL);
   assert(nlpi != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   par->Algorithm = nlpidata->useip ? 2 : 1;
   par->ScaledKKT = DEFAULT_SCALEDKKT;
   par->sKKTOnlyAcceptable = DEFAULT_SCALEDKKT;
   par->Infty = SCIPinfinity(scip);

   if( nlpparam.warmstart )
   {
      SCIPdebugMsg(scip, "warmstart parameter not supported by Worhp interface yet. Ignored.\n");
   }

   if( nlpparam.lobjlimit > -SCIP_REAL_MAX )
   {
      SCIPwarningMessage(scip, "lobjlimit parameter not supported by Worhp interface yet. Ignored.\n");
   }

   if( nlpparam.fastfail )
   {
      SCIPdebugMsg(scip, "fastfail parameter not supported by Worhp interface yet. Ignored.\n");
   }

   par->TolFeas = nlpparam.feastol;
   par->TolOpti = nlpparam.opttol;
   par->TolComp = nlpparam.opttol;
   par->Timeout = nlpparam.timelimit;
   par->MaxIter = nlpparam.iterlimit;
   par->NLPprint = nlpparam.verblevel - 1; /* Worhp verbosity levels: -1 = off, 0 = normal, 1 = debug, >1 = more debug */

#ifdef CHECKFUNVALUES
   /* activate gradient and hessian check */
   par->CheckValuesDF = TRUE;
   par->CheckValuesDG = TRUE;
   par->CheckValuesHM = TRUE;
#endif

   return SCIP_OKAY;
}

/*
 * Callback methods of NLP solver interface
 */

/** copy method of NLP interface (called when SCIP copies plugins) */
static
SCIP_DECL_NLPICOPY(nlpiCopyWorhp)
{
   SCIP_NLPIDATA* sourcedata;

   assert(sourcenlpi != NULL);

   sourcedata = SCIPnlpiGetData(sourcenlpi);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPincludeNlpSolverWorhp(scip, sourcedata->useip) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data */
static
SCIP_DECL_NLPIFREE(nlpiFreeWorhp)
{
   assert(nlpi != NULL);
   assert(nlpidata != NULL);
   assert(*nlpidata != NULL);

   SCIPfreeBlockMemory(scip, nlpidata);
   assert(*nlpidata == NULL);

   return SCIP_OKAY;
}  /*lint !e715*/

/** creates a problem instance */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemWorhp)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, problem) );
   assert( *problem != NULL );

   /* initialize problem */
   (*problem)->firstrun = TRUE;
   SCIP_CALL( SCIPnlpiOracleCreate(scip, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName(scip, (*problem)->oracle, name) );

   /* allocate memory for Worhp data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*problem)->opt) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*problem)->wsp) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*problem)->par) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(*problem)->cnt) );
   WorhpPreInit((*problem)->opt, (*problem)->wsp, (*problem)->par, (*problem)->cnt);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &(*problem)->randnumgen, DEFAULT_RANDSEED, TRUE) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** free a problem instance */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemWorhp)
{
   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   if( (*problem)->opt != NULL )
   {
      /* free memory for last solution information */
      invalidateSolution(scip, *problem);

      assert((*problem)->wsp != NULL);
      assert((*problem)->par != NULL);
      assert((*problem)->cnt != NULL);

      /* free Worhp data */
      SCIP_CALL( freeWorhp(*problem) );
      SCIPfreeBlockMemory(scip, &(*problem)->cnt);
      SCIPfreeBlockMemory(scip, &(*problem)->par);
      SCIPfreeBlockMemory(scip, &(*problem)->wsp);
      SCIPfreeBlockMemory(scip, &(*problem)->opt);
   }

   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(scip, &(*problem)->oracle) );
   }

   SCIPfreeRandom(scip, &(*problem)->randnumgen);
   SCIPfreeMemoryArrayNull(scip, &(*problem)->initguess);
   SCIPfreeBlockMemory(scip, problem);
   *problem = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

/** add variables */
static
SCIP_DECL_NLPIADDVARS( nlpiAddVarsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddVars(scip, problem->oracle, nvars, lbs, ubs, varnames) );

   SCIPfreeMemoryArrayNull(scip, &problem->initguess);
   invalidateSolution(scip, problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/


/** add constraints */
static
SCIP_DECL_NLPIADDCONSTRAINTS(nlpiAddConstraintsWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddConstraints(scip, problem->oracle,
         nconss, lhss, rhss,
         nlininds, lininds, linvals,
         exprs, names) );

   invalidateSolution(scip, problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected */
static
SCIP_DECL_NLPISETOBJECTIVE(nlpiSetObjectiveWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   /* We pass the objective gradient in dense form to WORHP, so if the sparsity of that gradient changes, we do not need
    * to reset WORHP (firstrun=TRUE).  However, if the sparsity of the Hessian matrix of the objective changes, then the
    * sparsity pattern of the Hessian of the Lagrangian may change.  Thus, reset Worhp if the objective was and/or
    * becomes nonlinear, but leave firstrun untouched if it was and stays linear.
    */
   if( expr != NULL || SCIPnlpiOracleIsConstraintNonlinear(problem->oracle, -1) )
      problem->firstrun = TRUE;

   SCIP_CALL( SCIPnlpiOracleSetObjective(scip, problem->oracle,
         constant, nlins, lininds, linvals, expr) );

   invalidateSolution(scip, problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change variable bounds */
static
SCIP_DECL_NLPICHGVARBOUNDS(nlpiChgVarBoundsWorhp)
{
#ifdef SCIP_DISABLED_CODE
   const SCIP_Real* oldlbs = SCIPnlpiOracleGetVarLbs(problem->oracle);
   const SCIP_Real* oldubs = SCIPnlpiOracleGetVarUbs(problem->oracle);
   int i;
#endif

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

#ifdef SCIP_DISABLED_CODE
   /* TODO check WORHP version here */
   /* So far, Worhp can not handle fixed variables (and fixed variables that have been unfixed) when applying a
    * restart. The following code needs to be removed when this has changed.
    */
   for( i = 0; i < nvars; ++i )
   {
      int index = indices[i];
      SCIPdebugMsg(scip, "change bounds of %d from [%g,%g] -> [%g,%g]\n", index, oldlbs[index], oldubs[index],
         lbs[i], ubs[i]);

      if( REALABS(lbs[i] - ubs[i]) <= problem->feastol )
         problem->firstrun = TRUE;
      else
         if( REALABS(oldlbs[index] - oldubs[index]) <= problem->feastol )
            problem->firstrun = TRUE;
   }
#endif

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(scip, problem->oracle, nvars, indices, lbs, ubs) );

   invalidateSolution(scip, problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change constraint bounds */
static
SCIP_DECL_NLPICHGCONSSIDES(nlpiChgConsSidesWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

#ifdef SCIP_DEBUG
   {
      SCIP_Real oldlhs;
      SCIP_Real oldrhs;
      int i;

      for( i = 0; i < nconss; ++i )
      {
         int index = indices[i];
         oldlhs = SCIPnlpiOracleGetConstraintLhs(problem->oracle, index);
         oldrhs = SCIPnlpiOracleGetConstraintRhs(problem->oracle, index);
         SCIPdebugMsg(scip, "change constraint side of %d from [%g,%g] -> [%g,%g]\n", index, oldlhs, oldrhs, lhss[i], rhss[i]);
      }
   }
#endif

   SCIP_CALL( SCIPnlpiOracleChgConsSides(scip, problem->oracle, nconss, indices, lhss, rhss) );

   invalidateSolution(scip, problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of variables */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(scip, problem->oracle, dstats) );

   SCIPfreeMemoryArrayNull(scip, &problem->initguess); /* @TODO keep initguess for remaining variables */

   invalidateSolution(scip, problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of constraints */
static
SCIP_DECL_NLPIDELCONSSET(nlpiDelConstraintSetWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(scip, problem->oracle, dstats) );

   invalidateSolution(scip, problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(scip, problem->oracle, idx, nvals, varidxs, vals) );

   invalidateSolution(scip, problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** replaces the expression of a constraint or objective */
static
SCIP_DECL_NLPICHGEXPR(nlpiChgExprWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExpr(scip, problem->oracle, idxcons, expr) );

   invalidateSolution(scip, problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change the constant offset in the objective */
static
SCIP_DECL_NLPICHGOBJCONSTANT(nlpiChgObjConstantWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(scip, problem->oracle, objconstant) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets initial guess for primal variables */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   if( primalvalues != NULL )
   {
      if( !problem->initguess )
      {
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle)) );
      }
      else
      {
         BMScopyMemoryArray(problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle));
      }
   }
   else
   {
      SCIPfreeMemoryArrayNull(scip, &problem->initguess);
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** tries to solve NLP */
static
SCIP_DECL_NLPISOLVE(nlpiSolveWorhp)
{
   Workspace* wsp = problem->wsp;
   Control* cnt = problem->cnt;
   OptVar* opt = problem->opt;
   Params* par = problem->par;
   int status;
   int i;

   SCIPdebugMsg(scip, "solve with parameters " SCIP_NLPPARAM_PRINT(param));

   SCIP_CALL( SCIPnlpiOracleResetEvalTime(scip, problem->oracle) );

   if( param.timelimit == 0.0 )
   {
      /* there is nothing we can do if we are not given any time */
      problem->lastniter = 0;
      problem->lasttime = 0.0;
      problem->lasttermstat = SCIP_NLPTERMSTAT_TIMELIMIT;
      problem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;

      return SCIP_OKAY;
   }

   problem->lastniter = -1;
   problem->lasttime  = -1.0;

   if( param.verblevel == 0 )
   {
      SetWorhpPrint(noprint);
   }
   else
   {
      /* TODO this should go to a function that prints to the SCIP message handler
       * all this doesn't seem threadsafe at all!
       */
      SetWorhpPrint(WorhpDefaultPrintFunction);
   }

   /* initialize Worhp data if necessary */
   if( problem->firstrun )
   {
      SCIP_CALL( freeWorhp(problem) );
      SCIP_CALL( initWorhp(scip, nlpi, problem) );
      problem->firstrun = FALSE;
   }
   else
   {
      SCIP_CALL( updateWorhp(problem) );
   }

   /* set parameters */
   InitParams(&status, par);

   if( status != OK )
      return SCIP_INVALIDCALL;

   SCIP_CALL( handleNlpParam(scip, nlpi, par, param) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPnlpiOraclePrintProblem(problem->oracle, nlpidata->messagehdlr, NULL) );
#endif

   /* set initial guess (if available) */
   if( problem->initguess != NULL )
   {
      BMScopyMemoryArray(problem->opt->X, problem->initguess, problem->opt->n);
   }
   else
   {
      SCIP_Real lb, ub;

      assert(problem->randnumgen != NULL);

      SCIPdebugMsg(scip, "Worhp started without initial primal values; make up starting guess by projecting 0 onto variable bounds\n");

      for( i = 0; i < problem->opt->n; ++i )
      {
         lb = SCIPnlpiOracleGetVarLbs(problem->oracle)[i];
         ub = SCIPnlpiOracleGetVarUbs(problem->oracle)[i];

         if( lb > 0.0 )
            problem->opt->X[i] = SCIPrandomGetReal(problem->randnumgen, lb, lb + MAXPERTURB*MIN(1.0, ub-lb));
         else if( ub < 0.0 )
            problem->opt->X[i] = SCIPrandomGetReal(problem->randnumgen, ub - MAXPERTURB*MIN(1.0, ub-lb), ub);
         else
            problem->opt->X[i] = SCIPrandomGetReal(problem->randnumgen,
               MAX(lb, -MAXPERTURB*MIN(1.0, ub-lb)), MIN(ub, MAXPERTURB*MIN(1.0, ub-lb)));
      }
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "start point:\n");
   for( i = 0; i < problem->opt->n; ++i )
   {
      SCIPdebugMsg(scip, "x[%d] = %f\n", i, problem->opt->X[i]);
   }
#endif

   /*
    * Worhp Reverse Communication loop.
    * In every iteration poll GetUserAction for the requested action, i.e. one
    * of {callWorhp, iterOutput, evalF, evalG, evalDF, evalDG, evalHM, fidif}.
    *
    * Make sure to reset the requested user action afterwards by calling
    * DoneUserAction, except for 'callWorhp' and 'fidif'.
    */
   while( cnt->status < TerminateSuccess && cnt->status > TerminateError && !SCIPisSolveInterrupted(scip) )
   {
      /*
       * Worhp's main routine.
       * Do not manually reset callWorhp, this is only done by the FD routines.
       */
      if( GetUserAction(cnt, callWorhp) )
      {
         Worhp(opt, wsp, par, cnt);
         /* No DoneUserAction! */
      }

      /*
       * Show iteration output.
       * The call to IterationOutput() may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, iterOutput) )
      {
         IterationOutput(opt, wsp, par, cnt);
         DoneUserAction(cnt, iterOutput);
      }

      /*
       * Evaluate the objective function.
       * The call to UserF may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalF) )
      {
         if( userF(scip, problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalF);
      }

      /*
       * Evaluate the constraints.
       * The call to UserG may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalG) )
      {
         if( userG(scip, problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalG);
      }

      /*
       * Evaluate the gradient of the objective function.
       * The call to UserDF may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalDF) )
      {
         if( userDF(scip, problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalDF);
      }

      /*
       * Evaluate the Jacobian of the constraints.
       * The call to UserDG may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalDG) )
      {
         if( userDG(scip, problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalDG);
      }

      /*
       * Evaluate the Hessian matrix of the Lagrange function (L = f + mu*g)
       * The call to UserHM may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalHM) )
      {
         if( userHM(scip, problem) != SCIP_OKAY)
            break;
         DoneUserAction(cnt, evalHM);
      }

      /*
       * Use finite differences with RC to determine derivatives
       * Do not reset fidif, this is done by the FD routine.
       */
      if( GetUserAction(cnt, fidif) )
      {
         WorhpFidif(opt, wsp, par, cnt);
         /* No DoneUserAction! */
      }
   }

   /* interpret Worhp result */
   if( SCIPisSolveInterrupted(scip) )
   {
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_INTERRUPT;
   }
   else if( cnt->status < TerminateSuccess && cnt->status > TerminateError )
   {
      SCIPwarningMessage(scip, "Worhp failed because of an invalid function evaluation!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERICERROR;
   }
   else
   {
      SCIP_CALL( evaluateWorhpRun(scip, problem) );
   }

   /* prints a status message with information about the current solver status */
   StatusMsg(opt, wsp, par, cnt);

   /* store statistics */
   problem->lastniter = wsp->MajorIter;
   problem->lasttime = GetTimerCont(&cnt->Timer);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives solution status */
static
SCIP_DECL_NLPIGETSOLSTAT(nlpiGetSolstatWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->lastsolstat;
}  /*lint !e715*/

/** gives termination reason */
static
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetTermstatWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->lasttermstat;
}  /*lint !e715*/

/** gives primal and dual solution values */
static
SCIP_DECL_NLPIGETSOLUTION( nlpiGetSolutionWorhp )
{
   assert(problem != NULL);

   if( primalvalues != NULL )
      *primalvalues = problem->lastprimal;

   if( consdualvalues != NULL )
      *consdualvalues = problem->lastdualcons;

   if( varlbdualvalues != NULL )
      *varlbdualvalues = problem->lastduallb;

   if( varubdualvalues != NULL )
      *varubdualvalues = problem->lastdualub;

   if( objval != NULL )
   {
      if( problem->lastprimal != NULL )
      {
         /* TODO store last solution value instead of reevaluating the objective function */
         SCIP_CALL( SCIPnlpiOracleEvalObjectiveValue(scip, problem->oracle, problem->lastprimal, objval) );
      }
      else
         *objval = SCIP_INVALID;
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives solve statistics */
static
SCIP_DECL_NLPIGETSTATISTICS(nlpiGetStatisticsWorhp)
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(statistics != NULL);

   statistics->niterations = problem->lastniter;
   statistics->totaltime = problem->lasttime;
   statistics->evaltime = SCIPnlpiOracleGetEvalTime(scip, problem->oracle);
   statistics->consviol = problem->wsp->FeasOrigMax;
   statistics->boundviol = 0.0;

   return SCIP_OKAY;
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for Worhp solver and includes it into SCIP, if Worhp is available */
SCIP_RETCODE SCIPincludeNlpSolverWorhp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             useip               /**< TRUE for using Interior Point, FALSE for SQP */
   )
{
   SCIP_NLPIDATA* nlpidata;
   char name[SCIP_MAXSTRLEN];
   int priority;

   /* create Worhp solver interface data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlpidata) );
   nlpidata->useip = useip;

   /* disable Worhp's keyboard handler, not useful here and not threadsafe */
   (void) setenv("WORHP_DISABLE_KEYBOARD_HANDLER", "1", 0);

#if DEFAULT_VERBLEVEL == 0
   /* disable Worhp output by default */
   SetWorhpPrint(noprint);
#endif

   /* checks the version of the library and header files */
   CHECK_WORHP_VERSION

   /* create solver interface */
   if( useip )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "worhp-ip");
      priority = NLPI_PRIORITY_IP;
   }
   else
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "worhp-sqp");
      priority = NLPI_PRIORITY_SQP;
   }

   SCIP_CALL( SCIPincludeNlpi(scip,
         name, NLPI_DESC, priority,
         nlpiCopyWorhp, nlpiFreeWorhp, NULL,
         nlpiCreateProblemWorhp, nlpiFreeProblemWorhp, NULL,
         nlpiAddVarsWorhp, nlpiAddConstraintsWorhp, nlpiSetObjectiveWorhp,
         nlpiChgVarBoundsWorhp, nlpiChgConsSidesWorhp, nlpiDelVarSetWorhp, nlpiDelConstraintSetWorhp,
         nlpiChgLinearCoefsWorhp, nlpiChgExprWorhp,
         nlpiChgObjConstantWorhp, nlpiSetInitialGuessWorhp, nlpiSolveWorhp, nlpiGetSolstatWorhp, nlpiGetTermstatWorhp,
         nlpiGetSolutionWorhp, nlpiGetStatisticsWorhp,
         nlpidata) );

   if( useip )  /* TODO lookup whether Worhp info has already been included instead of assuming that worhp-up will be included */
   {
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameWorhp(), SCIPgetSolverDescWorhp()) );
   }

   return SCIP_OKAY;
}

/** gets string that identifies Worhp (version number) */
const char* SCIPgetSolverNameWorhp(
   void
   )
{
#ifdef WORHP_VERSION
   return "WORHP " WORHP_VERSION;
#else
   static char solvername[20];
   sprintf(solvername, "WORHP %d.%d." WORHP_PATCH, WORHP_MAJOR, WORHP_MINOR);
   return solvername;
#endif
}

/** gets string that describes Worhp (version number) */
const char* SCIPgetSolverDescWorhp(
   void
   )
{
   return "Nonlinear programming solver developed at Research Institute Steinbeis (www.worhp.de)";
}

/** returns whether Worhp is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisWorhpAvailableWorhp(
   void
   )
{
   return TRUE;
}
