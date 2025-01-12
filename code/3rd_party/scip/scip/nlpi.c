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

/**@file   nlpi.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling NLP solver interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/nlpi.h"
#include "scip/pub_message.h"
#include "scip/pub_nlpi.h"
#include "scip/clock.h"
#include "scip/struct_nlpi.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"

/** compares two NLPIs w.r.t. their priority */
SCIP_DECL_SORTPTRCOMP(SCIPnlpiComp)
{  /*lint --e{715}*/
   return ((SCIP_NLPI*)elem2)->priority - ((SCIP_NLPI*)elem1)->priority;
}

/** creates an NLP solver interface */
SCIP_RETCODE SCIPnlpiCreate(
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying of NLPI, can be NULL */
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
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;

   assert(nlpi != NULL);
   assert(name != NULL);
   assert(description != NULL);
   assert(nlpifree != NULL);
   assert(nlpicreateproblem != NULL);
   assert(nlpifreeproblem != NULL);
   assert(nlpiaddvars != NULL);
   assert(nlpiaddconstraints != NULL);
   assert(nlpisetobjective != NULL);
   assert(nlpichgvarbounds != NULL);
   assert(nlpichgconssides != NULL);
   assert(nlpidelconsset != NULL);
   assert(nlpichglinearcoefs != NULL);
   assert(nlpichgobjconstant != NULL);
   assert(nlpisolve != NULL);
   assert(nlpigetsolstat != NULL);
   assert(nlpigettermstat != NULL);
   assert(nlpigetsolution != NULL);
   assert(nlpigetstatistics != NULL);

   SCIP_ALLOC( BMSallocClearMemory(nlpi) );

   if( BMSduplicateMemoryArray(&(*nlpi)->name, name, strlen(name)+1) == NULL )
   {
      BMSfreeMemory(nlpi);
      return SCIP_NOMEMORY;
   }

   if( BMSduplicateMemoryArray(&(*nlpi)->description, description, strlen(description)+1) == NULL )
   {
      BMSfreeMemoryArray(&(*nlpi)->name);
      BMSfreeMemory(nlpi);
      return SCIP_NOMEMORY;
   }

   (*nlpi)->priority = priority;
   (*nlpi)->nlpicopy = nlpicopy;
   (*nlpi)->nlpifree = nlpifree;
   (*nlpi)->nlpigetsolverpointer = nlpigetsolverpointer;
   (*nlpi)->nlpicreateproblem = nlpicreateproblem;
   (*nlpi)->nlpifreeproblem = nlpifreeproblem;
   (*nlpi)->nlpigetproblempointer = nlpigetproblempointer;
   (*nlpi)->nlpiaddvars = nlpiaddvars;
   (*nlpi)->nlpiaddconstraints = nlpiaddconstraints;
   (*nlpi)->nlpisetobjective = nlpisetobjective;
   (*nlpi)->nlpichgvarbounds = nlpichgvarbounds;
   (*nlpi)->nlpichgconssides = nlpichgconssides;
   (*nlpi)->nlpidelvarset = nlpidelvarset;
   (*nlpi)->nlpidelconsset = nlpidelconsset;
   (*nlpi)->nlpichglinearcoefs = nlpichglinearcoefs;
   (*nlpi)->nlpichgobjconstant = nlpichgobjconstant;
   (*nlpi)->nlpisetinitialguess = nlpisetinitialguess;
   (*nlpi)->nlpisolve = nlpisolve;
   (*nlpi)->nlpigetsolstat = nlpigetsolstat;
   (*nlpi)->nlpigettermstat = nlpigettermstat;
   (*nlpi)->nlpigetsolution = nlpigetsolution;
   (*nlpi)->nlpigetstatistics = nlpigetstatistics;
   (*nlpi)->nlpidata = nlpidata;

   retcode = SCIPclockCreate(&(*nlpi)->problemtime, SCIP_CLOCKTYPE_DEFAULT);
   if( retcode != SCIP_OKAY )
   {
      BMSfreeMemoryArray(&(*nlpi)->description);
      BMSfreeMemoryArray(&(*nlpi)->name);
      BMSfreeMemory(nlpi);
   }

   return retcode;
}

/** sets NLP solver priority */
void SCIPnlpiSetPriority(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   int                   priority            /**< new priority of NLPI */
   )
{
   assert(nlpi != NULL);

   nlpi->priority = priority;
}

/** copies an NLPI and includes it into another SCIP instance */
SCIP_RETCODE SCIPnlpiCopyInclude(
   SCIP_NLPI*            sourcenlpi,         /**< the NLP interface to copy */
   SCIP_SET*             targetset           /**< global SCIP settings where to include copy */
   )
{
   assert(sourcenlpi != NULL);
   assert(targetset != NULL);

   if( sourcenlpi->nlpicopy != NULL )
   {
      SCIP_CALL( sourcenlpi->nlpicopy(targetset->scip, sourcenlpi) );
   }

   return SCIP_OKAY;
}

/** frees NLPI */
SCIP_RETCODE SCIPnlpiFree(
   SCIP_NLPI**           nlpi,               /**< pointer to NLPI data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nlpi  != NULL);
   assert(*nlpi != NULL);
   assert(set   != NULL);

   if( (*nlpi)->nlpifree != NULL )
   {
      SCIP_CALL( (*nlpi)->nlpifree(set->scip, *nlpi, &(*nlpi)->nlpidata) );
      assert((*nlpi)->nlpidata == NULL);
   }
   BMSfreeMemoryArray(&(*nlpi)->name);
   BMSfreeMemoryArray(&(*nlpi)->description);

   SCIPclockFree(&(*nlpi)->problemtime);

   BMSfreeMemory(nlpi);

   assert(*nlpi == NULL);

   return SCIP_OKAY;
}

/** initializes NLPI */
void SCIPnlpiInit(
   SCIP_NLPI*            nlpi                /**< solver interface */
   )
{
   assert(nlpi != NULL);

   nlpi->nproblems = 0;
   nlpi->nsolves = 0;
   SCIPclockReset(nlpi->problemtime);
   nlpi->solvetime = 0.0;
   nlpi->evaltime = 0.0;
   nlpi->niter = 0L;
   BMSclearMemoryArray(nlpi->ntermstat, (int)SCIP_NLPTERMSTAT_OTHER+1);
   BMSclearMemoryArray(nlpi->nsolstat, (int)SCIP_NLPSOLSTAT_UNKNOWN+1);
}

/** gets pointer for NLP solver */
void* SCIPnlpiGetSolverPointer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance, or NULL */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);

   if( nlpi->nlpigetsolverpointer != NULL )
      return nlpi->nlpigetsolverpointer(set->scip, nlpi, problem);
   else
      return NULL;
}

/** creates a problem instance */
SCIP_RETCODE SCIPnlpiCreateProblem(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM**    problem,            /**< problem pointer to store the problem data */
   const char*           name                /**< name of problem, can be NULL */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpicreateproblem != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpicreateproblem(set->scip, nlpi, problem, name) );
   SCIPclockStop(nlpi->problemtime, set);

   ++nlpi->nproblems;

   return SCIP_OKAY;
}

/** frees a problem instance */
SCIP_RETCODE SCIPnlpiFreeProblem(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM**    problem             /**< pointer where problem instance is stored */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpifreeproblem != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpifreeproblem(set->scip, nlpi, problem) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** gets pointer to solver-internal problem instance */
void* SCIPnlpiGetProblemPointer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(problem != NULL);

   if( nlpi->nlpigetproblempointer != NULL )
      return nlpi->nlpigetproblempointer(set->scip, nlpi, problem);
   else
      return NULL;
}

/** add variables to nlpi */
SCIP_RETCODE SCIPnlpiAddVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      lbs,                /**< lower bounds of variables, can be NULL if -infinity */
   const SCIP_Real*      ubs,                /**< upper bounds of variables, can be NULL if +infinity */
   const char**          varnames            /**< names of variables, can be NULL */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpiaddvars != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpiaddvars(set->scip, nlpi, problem, nvars, lbs, ubs, varnames) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** add constraints to nlpi */
SCIP_RETCODE SCIPnlpiAddConstraints(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nconss,             /**< number of constraints */
   const SCIP_Real*      lhss,               /**< left hand sides of constraints, can be NULL if -infinity */
   const SCIP_Real*      rhss,               /**< right hand sides of constraints, can be NULL if +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   SCIP_EXPR**           exprs,              /**< expressions for nonlinear part of constraints, entry of array may be NULL in case of no nonlinear part, may be NULL in case of no nonlinear part in any constraint */
   const char**          names               /**< names of constraints, may be NULL or entries may be NULL */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpiaddconstraints != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpiaddconstraints(set->scip, nlpi, problem, nconss, lhss, rhss, nlininds, lininds, linvals, exprs, names) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** sets or overwrites objective, a minimization problem is expected */
SCIP_RETCODE SCIPnlpiSetObjective(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nlins,              /**< number of linear variables */
   const int*            lininds,            /**< variable indices, may be NULL in case of no linear part */
   const SCIP_Real*      linvals,            /**< coefficient values, may be NULL in case of no linear part */
   SCIP_EXPR*            expr,               /**< expression for nonlinear part of objective function, may be NULL in case of no nonlinear part */
   const SCIP_Real       constant            /**< objective value offset */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpisetobjective != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpisetobjective(set->scip, nlpi, problem, nlins, lininds, linvals, expr, constant) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** change variable bounds */
SCIP_RETCODE SCIPnlpiChgVarBounds(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgvarbounds != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpichgvarbounds(set->scip, nlpi, problem, nvars, indices, lbs, ubs) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** change constraint sides */
SCIP_RETCODE SCIPnlpiChgConsSides(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nconss,             /**< number of constraints to change sides */
   const int*            indices,            /**< indices of constraints to change sides */
   const SCIP_Real*      lhss,               /**< new left hand sides */
   const SCIP_Real*      rhss                /**< new right hand sides */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgconssides != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpichgconssides(set->scip, nlpi, problem, nconss, indices, lhss, rhss) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** delete a set of variables */
SCIP_RETCODE SCIPnlpiDelVarSet(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int*                  dstats,             /**< deletion status of vars; 1 if var should be deleted, 0 if not */
   int                   dstatssize          /**< size of the dstats array */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpidelvarset != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpidelvarset(set->scip, nlpi, problem, dstats, dstatssize) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** delete a set of constraints */
SCIP_RETCODE SCIPnlpiDelConsSet(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int*                  dstats,             /**< deletion status of constraints; 1 if constraint should be deleted, 0 if not */
   int                   dstatssize          /**< size of the dstats array */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpidelconsset != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpidelconsset(set->scip, nlpi, problem, dstats, dstatssize) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** changes or adds linear coefficients in a constraint or objective */
SCIP_RETCODE SCIPnlpiChgLinearCoefs(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   idx,                /**< index of constraint or -1 for objective */
   int                   nvals,              /**< number of values in linear constraint to change */
   const int*            varidxs,            /**< indices of variables which coefficient to change */
   const SCIP_Real*      vals                /**< new values for coefficients */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichglinearcoefs != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpichglinearcoefs(set->scip, nlpi, problem, idx, nvals, varidxs, vals) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** change the expression in the nonlinear part */
SCIP_RETCODE SCIPnlpiChgExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   idxcons,            /**< index of constraint or -1 for objective */
   SCIP_EXPR*            expr                /**< new expression for constraint or objective, or NULL to only remove previous tree */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgexpr != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpichgexpr(set->scip, nlpi, problem, idxcons, expr) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** change the constant offset in the objective */
SCIP_RETCODE SCIPnlpiChgObjConstant(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real             objconstant         /**< new value for objective constant */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpichgobjconstant != NULL);
   assert(problem != NULL);

   SCIPclockStart(nlpi->problemtime, set);
   SCIP_CALL( nlpi->nlpichgobjconstant(set->scip, nlpi, problem, objconstant) );
   SCIPclockStop(nlpi->problemtime, set);

   return SCIP_OKAY;
}

/** sets initial guess for primal variables */
SCIP_RETCODE SCIPnlpiSetInitialGuess(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real*            primalvalues,       /**< initial primal values for variables, or NULL to clear previous values */
   SCIP_Real*            consdualvalues,     /**< initial dual values for constraints, or NULL to clear previous values */
   SCIP_Real*            varlbdualvalues,    /**< initial dual values for variable lower bounds, or NULL to clear previous values */
   SCIP_Real*            varubdualvalues     /**< initial dual values for variable upper bounds, or NULL to clear previous values */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(problem != NULL);

   if( nlpi->nlpisetinitialguess != NULL )
   {
      SCIP_CALL( nlpi->nlpisetinitialguess(set->scip, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues) );
   }

   return SCIP_OKAY;
}

/** tries to solve NLP */
SCIP_RETCODE SCIPnlpiSolve(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM*        param               /**< solve parameters */
   )
{
   SCIP_NLPSTATISTICS stats;

   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpisolve != NULL);
   assert(problem != NULL);
   assert(param != NULL);

   /* check that parameter values are in accepted range (if type allows more than we would accept) */
   if( param->iterlimit < 0 )
   {
      SCIPerrorMessage("Value %d for parameter iteration limit must be non-negative.\n", param->iterlimit);
      return SCIP_PARAMETERWRONGVAL;
   }
   if( param->feastol < 0.0 )
   {
      SCIPerrorMessage("Value %g for parameter feasibility tolerance cannot be negative\n", param->feastol);
      return SCIP_PARAMETERWRONGVAL;
   }
   if( param->opttol < 0.0 )
   {
      SCIPerrorMessage("Value %g for parameter optimality tolerance cannot be negative\n", param->opttol);
      return SCIP_PARAMETERWRONGVAL;
   }
   if( param->solvertol < 0.0 )
   {
      SCIPerrorMessage("Value %g for parameter solver tolerance cannot be negative\n", param->solvertol);
      return SCIP_PARAMETERWRONGVAL;
   }
   if( param->timelimit < 0.0 )
   {
      SCIPerrorMessage("Value %g for parameter time limit cannot be negative\n", param->timelimit);
      return SCIP_PARAMETERWRONGVAL;
   }

   if( param->timelimit == SCIP_REAL_MAX && set->istimelimitfinite )  /*lint !e777*/
   {
      /* set timelimit to time remaining if limits/time has been set */
      param->timelimit = set->limit_time - SCIPclockGetTime(stat->solvingtime);
      if( param->timelimit < 0.0 )
         param->timelimit = 0.0;
      /* still call NLP solver if no time left to ensure proper termination codes */
   }

   ++nlpi->nsolves;

   SCIP_CALL( nlpi->nlpisolve(set->scip, nlpi, problem, *param) );

   /* coverity[overrun] */
   ++nlpi->ntermstat[nlpi->nlpigettermstat(set->scip, nlpi, problem)];
   /* coverity[overrun] */
   ++nlpi->nsolstat[nlpi->nlpigetsolstat(set->scip, nlpi, problem)];

   SCIP_CALL( nlpi->nlpigetstatistics(set->scip, nlpi, problem, &stats) );
   nlpi->solvetime += stats.totaltime;
   nlpi->evaltime += stats.evaltime;
   nlpi->niter += stats.niterations;

   return SCIP_OKAY;
}

/** gives solution status */
SCIP_NLPSOLSTAT SCIPnlpiGetSolstat(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetsolstat != NULL);
   assert(problem != NULL);

   return nlpi->nlpigetsolstat(set->scip, nlpi, problem);
}

/** gives termination reason */
SCIP_NLPTERMSTAT SCIPnlpiGetTermstat(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigettermstat != NULL);
   assert(problem != NULL);

   return nlpi->nlpigettermstat(set->scip, nlpi, problem);
}

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_RETCODE SCIPnlpiGetSolution(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real**           primalvalues,       /**< buffer to store pointer to array to primal values, or NULL if not needed */
   SCIP_Real**           consdualvalues,     /**< buffer to store pointer to array to dual values of constraints, or NULL if not needed */
   SCIP_Real**           varlbdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real**           varubdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real*            objval              /**< pointer to store the objective value, or NULL if not needed */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetsolution != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetsolution(set->scip, nlpi, problem, primalvalues, consdualvalues, varlbdualvalues, varubdualvalues, objval) );

   return SCIP_OKAY;
}

/** gives solve statistics */
SCIP_RETCODE SCIPnlpiGetStatistics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);
   assert(nlpi->nlpigetstatistics != NULL);
   assert(problem != NULL);

   SCIP_CALL( nlpi->nlpigetstatistics(set->scip, nlpi, problem, statistics) );

   return SCIP_OKAY;
}

/* from pub_nlpi.h */

#ifdef NDEBUG
/* Undo the defines from pub_nlpi.h, which exist if NDEBUG is defined. */
#undef SCIPnlpiGetData
#undef SCIPnlpiGetName
#undef SCIPnlpiGetDesc
#undef SCIPnlpiGetPriority
#undef SCIPnlpiSetPriority
#undef SCIPnlpiGetNProblems
#undef SCIPnlpiGetProblemTime
#undef SCIPnlpiGetNSolves
#undef SCIPnlpiGetSolveTime
#undef SCIPnlpiGetEvalTime
#undef SCIPnlpiGetNIterations
#undef SCIPnlpiGetNTermStat
#undef SCIPnlpiGetNSolStat
#endif

/** gets data of an NLPI */
SCIP_NLPIDATA* SCIPnlpiGetData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->nlpidata;
}

/** gets NLP solver name */
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->name;
}

/** gets NLP solver description */
const char* SCIPnlpiGetDesc(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->description;
}

/** gets NLP solver priority */
int SCIPnlpiGetPriority(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);

   return nlpi->priority;
}


/**@name Statistics */
/**@{ */

/** gives number of problems created for NLP solver so far */
int SCIPnlpiGetNProblems(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);
   return nlpi->nproblems;
}

/** gives total time spend in problem creation/modification/freeing */
SCIP_Real SCIPnlpiGetProblemTime(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);
   return SCIPclockGetTime(nlpi->problemtime);
}

/** total number of NLP solves so far */
int SCIPnlpiGetNSolves(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);
   return nlpi->nsolves;
}

/** gives total time spend in NLP solves (as reported by solver) */
SCIP_Real SCIPnlpiGetSolveTime(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);
   return nlpi->solvetime;
}

/** gives total time spend in function evaluation during NLP solves
 *
 * If parameter `timing/nlpieval` is off (the default), depending on the NLP solver, this may just return 0.
 */
SCIP_Real SCIPnlpiGetEvalTime(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);
   return nlpi->evaltime;
}

/** gives total number of iterations spend by NLP solver so far */
SCIP_Longint SCIPnlpiGetNIterations(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   )
{
   assert(nlpi != NULL);
   return nlpi->niter;
}

/** gives number of times a solve ended with a specific termination status */
int SCIPnlpiGetNTermStat(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   SCIP_NLPTERMSTAT      termstatus          /**< the termination status to query for */
   )
{
   assert(nlpi != NULL);
   return nlpi->ntermstat[termstatus];
}

/** gives number of times a solve ended with a specific solution status */
int SCIPnlpiGetNSolStat(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   SCIP_NLPSOLSTAT       solstatus           /**< the solution status to query for */
   )
{
   assert(nlpi != NULL);
   return nlpi->nsolstat[solstatus];
}

/** adds statistics from one NLPI to another */
void SCIPnlpiMergeStatistics(
   SCIP_NLPI*            targetnlpi,         /**< NLP interface where to add statistics */
   SCIP_NLPI*            sourcenlpi,         /**< NLP interface from which to add statistics */
   SCIP_Bool             reset               /**< whether to reset statistics in sourcescip */
   )
{
   int i;

   assert(targetnlpi != NULL);
   assert(sourcenlpi != NULL);

   targetnlpi->nproblems += sourcenlpi->nproblems;
   targetnlpi->nsolves += sourcenlpi->nsolves;
   SCIPclockSetTime(targetnlpi->problemtime, SCIPclockGetTime(targetnlpi->problemtime) + SCIPclockGetTime(sourcenlpi->problemtime));
   targetnlpi->solvetime += sourcenlpi->solvetime;
   targetnlpi->evaltime += sourcenlpi->evaltime;
   targetnlpi->niter += sourcenlpi->niter;

   for( i = (int)SCIP_NLPTERMSTAT_OKAY; i <= (int)SCIP_NLPTERMSTAT_OTHER; ++i )
      targetnlpi->ntermstat[i] += sourcenlpi->ntermstat[i];
   for( i = (int)SCIP_NLPSOLSTAT_GLOBOPT; i <= (int)SCIP_NLPSOLSTAT_UNKNOWN; ++i )
      targetnlpi->nsolstat[i] += sourcenlpi->nsolstat[i];

   if( reset )
   {
      sourcenlpi->nproblems = 0;
      sourcenlpi->nsolves = 0;
      SCIPclockReset(sourcenlpi->problemtime);
      sourcenlpi->solvetime = 0.0;
      sourcenlpi->evaltime = 0.0;
      sourcenlpi->niter = 0;

      for( i = (int)SCIP_NLPTERMSTAT_OKAY; i <= (int)SCIP_NLPTERMSTAT_OTHER; ++i )
         sourcenlpi->ntermstat[i] = 0;
      for( i = (int)SCIP_NLPSOLSTAT_GLOBOPT; i <= (int)SCIP_NLPSOLSTAT_UNKNOWN; ++i )
         sourcenlpi->nsolstat[i] = 0;
   }
}

/**@} */ /* Statistics */
