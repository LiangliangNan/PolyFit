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

/**@file   heur_adaptivediving.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  diving heuristic that selects adaptively between the existing, public dive sets
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_adaptivediving.h"
#include "scip/heuristics.h"
#include "scip/scipdefplugins.h"

#define HEUR_NAME             "adaptivediving"
#define HEUR_DESC             "diving heuristic that selects adaptively between the existing, public divesets"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_DIVING
#define HEUR_PRIORITY         -70000
#define HEUR_FREQ             5
#define HEUR_FREQOFS          3
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DIVESETS_INITIALSIZE 10
#define DEFAULT_INITIALSEED 13

/*
 * Default parameter settings
 */
#define DEFAULT_SELTYPE 'w'
#define DEFAULT_SCORETYPE 'c'                /**< score parameter for selection: minimize either average 'n'odes, LP 'i'terations,
                                              *  backtrack/'c'onflict ratio, 'd'epth, 1 / 's'olutions, or
                                              *  1 / solutions'u'ccess */
#define DEFAULT_USEADAPTIVECONTEXT FALSE
#define DEFAULT_SELCONFIDENCECOEFF 10.0      /**< coefficient c to decrease initial confidence (calls + 1.0) / (calls + c) in scores */
#define DEFAULT_EPSILON             1.0      /**< parameter that increases probability of exploration among divesets (only active if seltype is 'e') */
#define DEFAULT_MAXLPITERQUOT       0.1      /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS      1500L      /**< additional number of allowed LP iterations */
#define DEFAULT_BESTSOLWEIGHT      10.0      /**< weight of incumbent solutions compared to other solutions in computation of LP iteration limit */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   /* data structures used internally */
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator for selection */
   SCIP_DIVESET**        divesets;           /**< publicly available divesets from diving heuristics */
   int                   ndivesets;          /**< number of publicly available divesets from diving heuristics */
   int                   divesetssize;       /**< array size for divesets array */
   int                   lastselection;      /**< stores the last selected diveset when the heuristics was run */
   /* user parameters */
   SCIP_Real             epsilon;            /**< parameter that increases probability of exploration among divesets (only active if seltype is 'e') */
   SCIP_Real             selconfidencecoeff; /**< coefficient c to decrease initial confidence (calls + 1.0) / (calls + c) in scores */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   SCIP_Longint          maxlpiterofs;       /**< additional number of allowed LP iterations */
   SCIP_Real             bestsolweight;      /**< weight of incumbent solutions compared to other solutions in computation of LP iteration limit */
   char                  seltype;            /**< selection strategy: (e)psilon-greedy, (w)eighted distribution, (n)ext diving */
   char                  scoretype;          /**< score parameter for selection: minimize either average 'n'odes, LP 'i'terations,
                                               *  backtrack/'c'onflict ratio, 'd'epth, 1 / 's'olutions, or
                                               *  1 / solutions'u'ccess */
   SCIP_Bool             useadaptivecontext; /**< should the heuristic use its own statistics, or shared statistics? */
};

/*
 * local methods
 */


/** get the selection score for this dive set */
static
SCIP_RETCODE divesetGetSelectionScore(
   SCIP_DIVESET*         diveset,            /**< diving settings data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_DIVECONTEXT      divecontext,        /**< context for diving statistics */
   SCIP_Real*            scoreptr            /**< pointer to store the score */
   )
{
   SCIP_Real confidence;

   assert(scoreptr != NULL);

   /* compute confidence scalar (converges towards 1 with increasing number of calls) */
   confidence = (SCIPdivesetGetNCalls(diveset, divecontext) + 1.0) /
            (SCIPdivesetGetNCalls(diveset, divecontext) + heurdata->selconfidencecoeff);

   switch (heurdata->scoretype) {
      case 'n': /* min average nodes */
         *scoreptr = confidence * SCIPdivesetGetNProbingNodes(diveset, divecontext) / (SCIPdivesetGetNCalls(diveset, divecontext) + 1.0);
         break;
      case 'i': /* min avg LP iterations */
         *scoreptr = confidence * SCIPdivesetGetNLPIterations(diveset, divecontext) / (SCIPdivesetGetNCalls(diveset, divecontext) + 1.0);
         break;
      case 'c': /* min backtrack / conflict ratio (the current default) */
         *scoreptr = confidence * (SCIPdivesetGetNBacktracks(diveset, divecontext)) / (SCIPdivesetGetNConflicts(diveset, divecontext) + 10.0);
         break;
      case 'd': /* minimum average depth */
         *scoreptr = SCIPdivesetGetAvgDepth(diveset, divecontext) * confidence;
         break;
      case 's': /* maximum number of solutions */
         *scoreptr = confidence / (SCIPdivesetGetNSols(diveset, divecontext) + 1.0);
         break;
      case 'u': /* maximum solution success (which weighs best solutions higher) */
         *scoreptr = confidence / (SCIPdivesetGetSolSuccess(diveset, divecontext) + 1.0);
         break;
      default:
         SCIPerrorMessage("Unsupported scoring parameter '%c'\n", heurdata->scoretype);
         SCIPABORT();
         *scoreptr = SCIP_INVALID;
         return SCIP_PARAMETERWRONGVAL;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyAdaptivediving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurAdaptivediving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeAdaptivediving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( heurdata->divesets != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->divesets, heurdata->divesetssize);
   }

   SCIPfreeRandom(scip, &heurdata->randnumgen);

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** find publicly available divesets and store them */
static
SCIP_RETCODE findAndStoreDivesets(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   int h;
   SCIP_HEUR** heurs;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);

   heurs = SCIPgetHeurs(scip);

   heurdata->divesetssize = DIVESETS_INITIALSIZE;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->divesets, heurdata->divesetssize) );
   heurdata->ndivesets = 0;

   for( h = 0; h < SCIPgetNHeurs(scip); ++h )
   {
      int d;
      assert(heurs[h] != NULL);

      /* loop over divesets of this heuristic and check whether they are public */
      for( d = 0; d < SCIPheurGetNDivesets(heurs[h]); ++d )
      {
         SCIP_DIVESET* diveset = SCIPheurGetDivesets(heurs[h])[d];
         if( SCIPdivesetIsPublic(diveset) )
         {
            SCIPdebugMsg(scip, "Found publicly available diveset %s\n", SCIPdivesetGetName(diveset));

            if( heurdata->ndivesets == heurdata->divesetssize )
            {
               int newsize = 2 * heurdata->divesetssize;
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &heurdata->divesets, heurdata->divesetssize, newsize) );
               heurdata->divesetssize = newsize;
            }
            heurdata->divesets[heurdata->ndivesets++] = diveset;
         }
         else
         {
            SCIPdebugMsg(scip, "Skipping private diveset %s\n", SCIPdivesetGetName(diveset));
         }
      }
   }
   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitAdaptivediving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get and reset heuristic data */
   heurdata = SCIPheurGetData(heur);
   heurdata->lastselection = -1;
   if( heurdata->divesets != NULL )
   {
      /* we clear the list of collected divesets to ensure reproducability and consistent state across multiple runs
       * within the same SCIP data structure */
      SCIPfreeBlockMemoryArray(scip, &heurdata->divesets, heurdata->divesetssize);
      assert(heurdata->divesets == NULL);
   }

   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize random seed; use problem dimensions to vary initial order between different instances */
   SCIPsetRandomSeed(scip, heurdata->randnumgen,
         (unsigned int)(DEFAULT_INITIALSEED + SCIPgetNOrigVars(scip) + SCIPgetNOrigConss(scip)));

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitAdaptivediving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** get LP iteration limit for diving */
static
SCIP_Longint getLPIterlimit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_Real nsolsfound = SCIPheurGetNSolsFound(heur) + heurdata->bestsolweight * SCIPheurGetNBestSolsFound(heur);
   SCIP_Longint nlpiterations = SCIPgetNNodeLPIterations(scip);
   SCIP_Longint ncalls = SCIPheurGetNCalls(heur);

   SCIP_Longint nlpiterationsdive = 0;
   SCIP_Longint lpiterlimit;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(heurdata != NULL);

   /* loop over the divesets and collect their individual iterations */
   for( i = 0; i < heurdata->ndivesets; ++i )
   {
      nlpiterationsdive += SCIPdivesetGetNLPIterations(heurdata->divesets[i], SCIP_DIVECONTEXT_ADAPTIVE);
   }

   /* compute the iteration limit */
   lpiterlimit = (SCIP_Longint)(heurdata->maxlpiterquot * (nsolsfound+1.0)/(ncalls+1.0) * nlpiterations);
   lpiterlimit += heurdata->maxlpiterofs;
   lpiterlimit -= nlpiterationsdive;

   return lpiterlimit;
}

#ifdef SCIP_DEBUG
/** print array for debug purpose */
static
char* printRealArray(
   char*                 strbuf,             /**< string buffer array */
   SCIP_Real*            elems,              /**< array elements */
   int                   nelems              /**< number of elements */
   )
{
   int c;
   char* pos = strbuf;

   for( c = 0; c < nelems; ++c )
   {
      pos += sprintf(pos, "%.4f ", elems[c]);
   }

   return strbuf;
}
#endif

/** sample from a distribution defined by weights */ /*lint -e715*/
static
int sampleWeighted(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      rng,                /**< random number generator */
   SCIP_Real*            weights,            /**< weights of a ground set that define the sampling distribution */
   int                   nweights            /**< number of elements in the ground set */
   )
{
   SCIP_Real weightsum;
   SCIP_Real randomnr;
   int w;
#ifdef SCIP_DEBUG
   char strbuf[SCIP_MAXSTRLEN];
   SCIPdebugMsg(scip, "Weights: %s\n", printRealArray(strbuf, weights, nweights));
#endif

   weightsum = 0.0;
   /* collect sum of weights */
   for( w = 0; w < nweights; ++w )
   {
      weightsum += weights[w];
   }
   assert(weightsum > 0);

   randomnr = SCIPrandomGetReal(rng, 0.0, weightsum);

   weightsum = 0.0;
   /* choose first element i such that the weight sum exceeds the random number */
   for( w = 0; w < nweights - 1; ++w )
   {
      weightsum += weights[w];

      if( weightsum >= randomnr )
         break;
   }
   assert(w < nweights);
   assert(weights[w] > 0.0);

   return w;
}

/** select the diving method to apply */
static
SCIP_RETCODE selectDiving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  selection           /**< selection made */
   )
{
   SCIP_Bool* methodunavailable;
   SCIP_DIVESET** divesets;
   int ndivesets;
   int d;
   SCIP_RANDNUMGEN* rng;
   SCIP_DIVECONTEXT divecontext;
   SCIP_Real* weights;
   SCIP_Real epsilon_t;

   divesets = heurdata->divesets;
   ndivesets = heurdata->ndivesets;
   assert(ndivesets > 0);
   assert(divesets != NULL);

   SCIP_CALL( SCIPallocClearBufferArray(scip, &methodunavailable, ndivesets) );

   divecontext = heurdata->useadaptivecontext ? SCIP_DIVECONTEXT_ADAPTIVE : SCIP_DIVECONTEXT_TOTAL;

   /* check availability of divesets */
   for( d = 0; d < heurdata->ndivesets; ++d )
   {
      SCIP_Bool available;
      SCIP_CALL( SCIPisDivesetAvailable(scip, heurdata->divesets[d], &available) );
      methodunavailable[d] = ! available;
   }

   *selection = -1;

   rng = heurdata->randnumgen;
   assert(rng != NULL);

   switch (heurdata->seltype) {
   case 'e':
      epsilon_t = heurdata->epsilon * sqrt(ndivesets / (SCIPheurGetNCalls(heur) + 1.0));
      epsilon_t = MAX(epsilon_t, 0.05);

      /* select one of the available methods at random */
      if( epsilon_t >= 1.0 || SCIPrandomGetReal(rng, 0.0, 1.0) < epsilon_t )
      {
         do
         {
            *selection = SCIPrandomGetInt(rng, 0, ndivesets - 1);
         }
         while( methodunavailable[*selection] );
      }
      else
      {
         SCIP_Real bestscore = SCIP_REAL_MAX;
         for( d = 0; d < heurdata->ndivesets; ++d )
         {
            SCIP_Real score;

            if( methodunavailable[d] )
               continue;

            SCIP_CALL( divesetGetSelectionScore(divesets[d], heurdata, divecontext, &score) );

            if( score < bestscore )
            {
               bestscore = score;
               *selection = d;
            }
         }
      }
      break;
   case 'w':
      SCIP_CALL( SCIPallocBufferArray(scip, &weights, ndivesets) );

      /* initialize weights as inverse of the score + a small positive epsilon */
      for( d = 0; d < ndivesets; ++d )
      {
         SCIP_Real score;

         SCIP_CALL( divesetGetSelectionScore(divesets[d], heurdata, divecontext, &score) );

         weights[d] = methodunavailable[d] ? 0.0 : 1 / (score + 1e-4);
      }

      *selection = sampleWeighted(scip, rng, weights, ndivesets);

      SCIPfreeBufferArray(scip, &weights);
      break;
   case 'n':
         /* continue from last selection and stop at the next available method */
         *selection = heurdata->lastselection;

         do
         {
            *selection = (*selection + 1) % ndivesets;
         }
         while( methodunavailable[*selection] );
         heurdata->lastselection = *selection;
      break;
   default:
      SCIPerrorMessage("Error: Unknown selection method %c\n", heurdata->seltype);

      return SCIP_INVALIDDATA;
   }

   assert(*selection >= 0 && *selection < ndivesets);
   SCIPfreeBufferArray(scip, &methodunavailable);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecAdaptivediving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;
   SCIP_DIVESET** divesets;
   SCIP_Longint lpiterlimit;
   int selection;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   heurdata = SCIPheurGetData(heur);
   if( heurdata->divesets == NULL )
   {
      SCIP_CALL( findAndStoreDivesets(scip, heur, heurdata) );
   }

   divesets = heurdata->divesets;
   assert(divesets != NULL);
   assert(heurdata->ndivesets > 0);

   SCIPdebugMsg(scip, "heurExecAdaptivediving: depth %d sols %d inf %u node %lld (last dive at %lld)\n",
         SCIPgetDepth(scip),
         SCIPgetNSols(scip),
         nodeinfeasible,
         SCIPgetNNodes(scip),
         SCIPgetLastDivenode(scip)
         );

   *result = SCIP_DELAYED;

   /* do not call heuristic in node that was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
   {
      SCIPdebugMsg(scip, "already dived at node here\n");

      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTRUN;

   lpiterlimit = getLPIterlimit(scip, heur, heurdata);

   if( lpiterlimit <= 0 )
      return SCIP_OKAY;

   /* select the next diving strategy based on previous success */
   SCIP_CALL( selectDiving(scip, heur, heurdata, &selection) );
   assert(selection >= 0 && selection < heurdata->ndivesets);

   diveset = divesets[selection];
   assert(diveset != NULL);

   SCIPdebugMsg(scip, "Selected diveset %s\n", SCIPdivesetGetName(diveset));

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible,
         lpiterlimit, SCIP_DIVECONTEXT_ADAPTIVE) );

   if( *result == SCIP_FOUNDSOL )
   {
      SCIPdebugMsg(scip, "Solution found by diveset %s\n", SCIPdivesetGetName(diveset));
   }

   return SCIP_OKAY;
}

/** creates the adaptivediving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurAdaptivediving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RETCODE retcode;
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create adaptivediving data */
   heurdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   heurdata->divesets = NULL;
   heurdata->ndivesets = 0;
   heurdata->divesetssize = -1;

   SCIP_CALL_TERMINATE( retcode, SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_INITIALSEED, TRUE), TERMINATE );

   /* include adaptive diving primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecAdaptivediving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyAdaptivediving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeAdaptivediving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitAdaptivediving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitAdaptivediving) );

   /* add parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/epsilon",
         "parameter that increases probability of exploration among divesets (only active if seltype is 'e')",
         &heurdata->epsilon, FALSE, DEFAULT_EPSILON, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/scoretype",
         "score parameter for selection: minimize either average 'n'odes, LP 'i'terations,"
         "backtrack/'c'onflict ratio, 'd'epth, 1 / 's'olutions, or 1 / solutions'u'ccess",
         &heurdata->scoretype, FALSE, DEFAULT_SCORETYPE, "cdinsu", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/seltype",
         "selection strategy: (e)psilon-greedy, (w)eighted distribution, (n)ext diving",
         &heurdata->seltype, FALSE, DEFAULT_SELTYPE, "enw", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useadaptivecontext",
         "should the heuristic use its own statistics, or shared statistics?", &heurdata->useadaptivecontext, TRUE,
         DEFAULT_USEADAPTIVECONTEXT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/selconfidencecoeff",
         "coefficient c to decrease initial confidence (calls + 1.0) / (calls + c) in scores",
         &heurdata->selconfidencecoeff, FALSE, DEFAULT_SELCONFIDENCECOEFF, 1.0, (SCIP_Real)INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0L, (SCIP_Longint)INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/bestsolweight",
         "weight of incumbent solutions compared to other solutions in computation of LP iteration limit",
         &heurdata->bestsolweight, FALSE, DEFAULT_BESTSOLWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

/* cppcheck-suppress unusedLabel */
TERMINATE:
   if( retcode != SCIP_OKAY )
   {
      SCIPfreeMemory(scip, &heurdata);
      return retcode;
   }

   return SCIP_OKAY;
}
