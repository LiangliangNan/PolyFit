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

/**@file    nlpi_all.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   NLP interface that uses all available NLP interfaces
 * @author  Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nlpi_all.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_nlpi.h"

#include <string.h>

#define NLPI_NAME              "all"                       /**< short concise name of solver */
#define NLPI_DESC              "NLP interface that uses all available NLP interfaces" /**< description of solver */
#define NLPI_PRIORITY          -3000                       /**< priority of NLP solver */

/*
 * Data structures
 */

struct SCIP_NlpiData
{
   SCIP_NLPI**           nlpis;              /**< array containing all nlpis */
   int                   nnlpis;             /**< total number of nlpis */
};

struct SCIP_NlpiProblem
{
   SCIP_NLPIPROBLEM**    nlpiproblems;       /**< array containing all nlpi problems */
   int                   nnlpiproblems;      /**< total number of nlpi problems */
   int                   bestidx;            /**< index of NLP solver with the best solution */
};

#ifdef SCIP_STATISTIC
static int _nnlps = 0;                       /**< number of NLPs that have been solved */
#endif

/*
 * Local methods
 */

/*
 * Callback methods of NLP solver interface
 */

/** copy method of NLP interface (called when SCIP copies plugins) */
static
SCIP_DECL_NLPICOPY(nlpiCopyAll)
{
   /* include NLPI */
   SCIP_CALL( SCIPincludeNlpSolverAll(scip) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data */
static
SCIP_DECL_NLPIFREE(nlpiFreeAll)
{
   assert(nlpi != NULL);
   assert(nlpidata != NULL);
   assert(*nlpidata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlpidata)->nlpis, (*nlpidata)->nnlpis);
   SCIPfreeBlockMemory(scip, nlpidata);
   assert(*nlpidata == NULL);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** creates a problem instance */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemAll)
{
   SCIP_NLPIDATA* data;
   int i;

   assert(nlpi    != NULL);
   assert(problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, problem) );

   /* initialize problem */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*problem)->nlpiproblems, data->nnlpis) );
   (*problem)->nnlpiproblems = data->nnlpis;

   for( i = 0; i < data->nnlpis; ++i )
   {
      assert(data->nlpis[i] != NULL);
      SCIP_CALL( SCIPcreateNlpiProblem(scip, data->nlpis[i], &((*problem)->nlpiproblems[i]), name) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** free a problem instance */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemAll)
{
   SCIP_NLPIDATA* data;
   int i;

   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   for( i = 0; i < data->nnlpis; ++i )
   {
      assert(data->nlpis[i] != NULL);
      SCIP_CALL( SCIPfreeNlpiProblem(scip, data->nlpis[i], &(*problem)->nlpiproblems[i]) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*problem)->nlpiproblems, data->nnlpis);
   SCIPfreeBlockMemory(scip, problem);

   return SCIP_OKAY;
}  /*lint !e715*/

/** add variables */
static
SCIP_DECL_NLPIADDVARS(nlpiAddVarsAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPaddNlpiVars(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], nvars, lbs, ubs, varnames) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/


/** add constraints */
static
SCIP_DECL_NLPIADDCONSTRAINTS(nlpiAddConstraintsAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPaddNlpiConstraints(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], nconss, lhss, rhss,
         nlininds, lininds, linvals, exprs, names) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected */
static
SCIP_DECL_NLPISETOBJECTIVE(nlpiSetObjectiveAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPsetNlpiObjective(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], nlins, lininds, linvals, expr, constant) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change variable bounds */
static
SCIP_DECL_NLPICHGVARBOUNDS(nlpiChgVarBoundsAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPchgNlpiVarBounds(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], nvars, indices, lbs, ubs) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change constraint bounds */
static
SCIP_DECL_NLPICHGCONSSIDES(nlpiChgConsSidesAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPchgNlpiConsSides(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], nconss, indices, lhss, rhss) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** delete a set of variables */
static
SCIP_DECL_NLPIDELVARSET(nlpiDelVarSetAll)
{
   SCIP_NLPIDATA* nlpidata;
   int* tmpdstats;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tmpdstats, dstatssize) );

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      if( i < nlpidata->nnlpis -1 )
      {
         /* restore dstats entries */
         BMScopyMemoryArray(tmpdstats, dstats, dstatssize);

         SCIP_CALL( SCIPdelNlpiVarSet(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], tmpdstats, dstatssize) );
      }
      else
      {
         /* NOTE this works only when all dstats array are the same after calling the nlpidelvarset callback
          * As long as all solvers use the SCIP NLPI oracle to store the NLP problem data, this is the case.
          * @TODO Assert that the returned dstats are all the same?
          */
         SCIP_CALL( SCIPdelNlpiVarSet(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], dstats, dstatssize) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &tmpdstats, dstatssize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** delete a set of constraints */
static
SCIP_DECL_NLPIDELCONSSET( nlpiDelConstraintSetAll )
{
   SCIP_NLPIDATA* nlpidata;
   int* tmpdstats;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tmpdstats, dstatssize) );

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      if( i < nlpidata->nnlpis - 1 )
      {
         /* restore dstats entries */
         BMScopyMemoryArray(tmpdstats, dstats, dstatssize);

         SCIP_CALL( SCIPdelNlpiConsSet(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], tmpdstats, dstatssize) );
      }
      else
      {
         /* NOTE this works only when all dstats array are the same after calling the nlpidelconsset callback
          * As long as all solvers use the SCIP NLPI oracle to store the NLP problem data, this is the case.
          * @TODO Assert that the returned dstats are all the same?
          */
         SCIP_CALL( SCIPdelNlpiConsSet(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], dstats, dstatssize) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &tmpdstats, dstatssize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective */
static
SCIP_DECL_NLPICHGLINEARCOEFS(nlpiChgLinearCoefsAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPchgNlpiLinearCoefs(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], idx, nvals, varidxs, vals) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** replaces the expression of a constraint or objective */
static
SCIP_DECL_NLPICHGEXPR(nlpiChgExprAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPchgNlpiExpr(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], idxcons, expr) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change the constant offset in the objective */
static
SCIP_DECL_NLPICHGOBJCONSTANT(nlpiChgObjConstantAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPchgNlpiObjConstant(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], objconstant) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets initial guess for primal variables */
static
SCIP_DECL_NLPISETINITIALGUESS(nlpiSetInitialGuessAll)
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPsetNlpiInitialGuess(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], primalvalues, consdualvalues,
            varlbdualvalues, varubdualvalues) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** tries to solve NLP */
static
SCIP_DECL_NLPISOLVE(nlpiSolveAll)
{
   SCIP_NLPIDATA* nlpidata;
   SCIP_NLPTERMSTAT besttermstat;
   SCIP_NLPSOLSTAT bestsolstat;
   SCIP_Real bestsolval;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   /* use first solver per default */
   problem->bestidx = 0;

   /* initialize best solution values */
   besttermstat = SCIP_NLPTERMSTAT_OTHER;
   bestsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   bestsolval = SCIPinfinity(scip);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      SCIP_NLPTERMSTAT termstat;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Real solval;
      SCIP_Bool update;

      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      /* solve NLP */
      SCIP_CALL( SCIPsolveNlpiParam(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], param) );

      termstat = SCIPgetNlpiTermstat(scip, nlpidata->nlpis[i], problem->nlpiproblems[i]);
      solstat = SCIPgetNlpiSolstat(scip, nlpidata->nlpis[i], problem->nlpiproblems[i]);
      solval = SCIPinfinity(scip);
      update = FALSE;

      /* collect solution value */
      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         SCIP_CALL( SCIPgetNlpiSolution(scip, nlpidata->nlpis[i], problem->nlpiproblems[i],
               NULL, NULL, NULL, NULL, &solval) );
         assert(!SCIPisInfinity(scip, solval));
      }

      /* better termination status -> update best solver */
      if( termstat < besttermstat )
         update = TRUE;

      /* no feasible solutions have been found so far -> update best solver */
      else if( bestsolstat >= SCIP_NLPSOLSTAT_LOCINFEASIBLE && solstat <= SCIP_NLPSOLSTAT_LOCINFEASIBLE )
         update = TRUE;

      /* use solver with the better solution value */
      else if( solval < bestsolval )
         update = TRUE;

      /* update best solver */
      if( update )
      {
         besttermstat = termstat;
         bestsolstat = solstat;
         bestsolval = solval;
         problem->bestidx = i;
      }

#ifdef SCIP_STATISTIC
      {
         SCIP_NLPSTATISTICS stats;

         SCIP_CALL( SCIPgetNlpiStatistics(scip, nlpidata->nlpis[i], problem->nlpiproblems[i], &stats) );

         SCIPstatisticMessage("%d solver %s termstat %d solstat %d solval %e iters %d time %g\n",
            _nnlps, SCIPnlpiGetName(nlpidata->nlpis[i]), termstat, solstat, solval,
            stats.niterations, stats.totaltime);
      }
#endif

      /* don't try more NLP solvers if allowed time is exceeded or SCIP is asked to interrupt */
      if( termstat == SCIP_NLPTERMSTAT_TIMELIMIT || termstat == SCIP_NLPTERMSTAT_INTERRUPT )
         break;
   }

#ifdef SCIP_STATISTIC
   ++_nnlps;
#endif

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives solution status */
static
SCIP_DECL_NLPIGETSOLSTAT(nlpiGetSolstatAll)
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* return the solution status of the first nlpi */
   return SCIPgetNlpiSolstat(scip, nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx]);
}

/** gives termination reason */
static
SCIP_DECL_NLPIGETTERMSTAT(nlpiGetTermstatAll)
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* return the solution status of the first nlpi */
   return SCIPgetNlpiTermstat(scip, nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx]);
}

/** gives primal and dual solution values */
static
SCIP_DECL_NLPIGETSOLUTION(nlpiGetSolutionAll)
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* return the solution status of the first nlpi */
   SCIP_CALL( SCIPgetNlpiSolution(scip, nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx],
         primalvalues, consdualvalues, varlbdualvalues, varubdualvalues, objval) );

   return SCIP_OKAY;
}

/** gives solve statistics */
static
SCIP_DECL_NLPIGETSTATISTICS(nlpiGetStatisticsAll)
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* collect statistics of the first solver */
   SCIP_CALL( SCIPgetNlpiStatistics(scip, nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx],
         statistics) );

   return SCIP_OKAY;
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for the solver "All" and includes it into SCIP, if at least 2 NLPIs have already been included
 *
 * This method should be called after all other NLP solver interfaces have been included.
 */
SCIP_RETCODE SCIPincludeNlpSolverAll(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   assert(scip != NULL);

   /* the number of NLPIs so far must be >= 2 */
   if( SCIPgetNNlpis(scip) < 2 )
      return SCIP_OKAY;

   /* create all solver interface data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &nlpidata) );

   nlpidata->nnlpis = SCIPgetNNlpis(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlpidata->nlpis, nlpidata->nnlpis) );

   /* copy nlpi pointers TODO should not need that */
   for( i = 0; i < nlpidata->nnlpis; ++i )
      nlpidata->nlpis[i] = SCIPgetNlpis(scip)[i];

   /* create solver interface */
   SCIP_CALL( SCIPincludeNlpi(scip,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyAll, nlpiFreeAll, NULL,
         nlpiCreateProblemAll, nlpiFreeProblemAll, NULL,
         nlpiAddVarsAll, nlpiAddConstraintsAll, nlpiSetObjectiveAll,
         nlpiChgVarBoundsAll, nlpiChgConsSidesAll, nlpiDelVarSetAll, nlpiDelConstraintSetAll,
         nlpiChgLinearCoefsAll, nlpiChgExprAll,
         nlpiChgObjConstantAll, nlpiSetInitialGuessAll, nlpiSolveAll, nlpiGetSolstatAll, nlpiGetTermstatAll,
         nlpiGetSolutionAll, nlpiGetStatisticsAll,
         nlpidata) );

   return SCIP_OKAY;
}
