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

/**@file   heur_reoptsols.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  reoptsols primal heuristic
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/heur_reoptsols.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_reopt.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include <string.h>


#define HEUR_NAME             "reoptsols"
#define HEUR_DESC             "primal heuristic updating solutions found in a previous optimization round"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_PROP
#define HEUR_PRIORITY         40000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
   int maxsols;                     /**< maximal number of solution to update per run */
   int maxruns;                     /**< check solutions of the last k runs only */

   /* statistics */
   int ncheckedsols;                /**< number of updated solutions */
   int nimprovingsols;              /**< number of improving solutions */
};


/*
 * Local methods
 */


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< the current heuristic */
   SCIP_SOL*             sol,                /**< solution of the subproblem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int        nvars;                         /* the original problem's number of variables */
   SCIP_Real* solvals;                       /* solution values of the subproblem */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, solvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, solvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyReoptsols)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurReoptsols(scip) );

   return SCIP_OKAY;
}

/* free data of the heuristic */
static
SCIP_DECL_HEURFREE(heurFreeReoptsols)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL );
   assert(heur != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL );

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/* initialize the heuristic */
static SCIP_DECL_HEURINIT(heurInitReoptsols)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL );
   assert(heur != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL );

   heurdata->ncheckedsols = 0;
   heurdata->nimprovingsols = 0;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecReoptsols)
{/*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   SCIP_SOL** sols;
   SCIP_Real objsimsol;
   SCIP_Bool sepabestsol;
   int allocmem;
   int nchecksols;
   int nsolsadded;
#ifdef SCIP_MORE_DEBUG
   int nsolsaddedrun;
#endif
   int run;
   int max_run;

   assert(heur != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   if( !SCIPisReoptEnabled(scip) )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   max_run = heurdata->maxruns == -1 ? 0 : MAX(0, SCIPgetNReoptRuns(scip)-1-heurdata->maxruns); /*lint !e666*/
   nchecksols = heurdata->maxsols == -1 ? INT_MAX : heurdata->maxsols;

   SCIP_CALL( SCIPgetRealParam(scip, "reoptimization/objsimsol", &objsimsol) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/sepabestsol", &sepabestsol) );

   /* allocate a buffer array to store the solutions */
   allocmem = 20;
   SCIP_CALL( SCIPallocBufferArray(scip, &sols, allocmem) );

   nsolsadded = 0;

   for( run = SCIPgetNReoptRuns(scip); run > max_run && nchecksols > 0; run-- )
   {
      SCIP_Real sim;
      int nsols;

#ifdef SCIP_MORE_DEBUG
      nsolsaddedrun = 0;
#endif
      nsols = 0;

      if( objsimsol == -1 )
         sim = 1;
      else
         sim = SCIPgetReoptSimilarity(scip, run, SCIPgetNReoptRuns(scip)-1);

      if( sim == SCIP_INVALID ) /*lint !e777*/
         return SCIP_INVALIDRESULT;

      if( sim >= objsimsol )
      {
         int s;

         /* get solutions of a specific run */
         SCIP_CALL( SCIPgetReoptSolsRun(scip, run, sols, allocmem, &nsols) );

         /* check memory and reallocate */
         if( nsols >= allocmem )
         {
            allocmem = nsols;
            SCIP_CALL( SCIPreallocBufferArray(scip, &sols, allocmem) );

            SCIP_CALL( SCIPgetReoptSolsRun(scip, run, sols, allocmem, &nsols) );
         }
         assert(nsols <= allocmem);

         /* update the solutions
          * stop, if the maximal number of solutions to be checked is reached
          */
         for( s = 0; s < nsols && nchecksols > 0; s++ )
         {
            SCIP_SOL* sol;
            SCIP_Real objsol;

            sol = sols[s];

            SCIP_CALL( SCIPrecomputeSolObj(scip, sol) );
            objsol = SCIPgetSolTransObj(scip, sol);

            /* we do not want to add solutions with objective value +infinity.
             * moreover, the solution should improve the current primal bound
             */
            if( !SCIPisInfinity(scip, objsol) && !SCIPisInfinity(scip, -objsol)
              && SCIPisFeasLT(scip, objsol, SCIPgetCutoffbound(scip)) )
            {
               SCIP_Bool stored;

               /* create a new solution */
               SCIP_CALL( createNewSol(scip, heur, sol, &stored) );

               if( stored )
               {
                  nsolsadded++;
#ifdef SCIP_MORE_DEBUG
                  nsolsaddedrun++;
#endif
                  heurdata->nimprovingsols++;
               }
            }

            nchecksols--;
            heurdata->ncheckedsols++;
         }
      }
#ifdef SCIP_MORE_DEBUG
         printf(">> heuristic <%s> found %d of %d improving solutions from run %d.\n", HEUR_NAME, nsolsaddedrun, nsols, run);
#endif
      }

   SCIPdebugMsg(scip, ">> heuristic <%s> found %d improving solutions.\n", HEUR_NAME, nsolsadded);

   if( nsolsadded > 0 )
      *result = SCIP_FOUNDSOL;
   else
      *result = SCIP_DIDNOTFIND;

   /* reset the marks of the checked solutions */
   SCIPresetReoptSolMarks(scip);

   /* free the buffer array */
   SCIPfreeBufferArray(scip, &sols);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/* returns the number of checked solutions */
int SCIPreoptsolsGetNCheckedsols(
   SCIP*                 scip
   )
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);

   heur = SCIPfindHeur(scip, HEUR_NAME);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->ncheckedsols;
}

/* returns the number of found improving solutions */
int SCIPreoptsolsGetNImprovingsols(
   SCIP*                 scip
   )
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);

   heur = SCIPfindHeur(scip, HEUR_NAME);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->nimprovingsols;
}

/** creates the reoptsols primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurReoptsols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create reoptsols primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecReoptsols, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyReoptsols) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeReoptsols) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitReoptsols) );

   /* parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxsols", "maximal number solutions which should be checked. (-1: all)",
         &heurdata->maxsols, TRUE, 1000, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxruns", "check solutions of the last k runs. (-1: all)",
         &heurdata->maxruns, TRUE, -1, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
