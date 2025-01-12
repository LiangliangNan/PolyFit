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

/**@file   cutsel_hybrid.c
 * @ingroup DEFPLUGINS_CUTSEL
 * @brief  hybrid cut selector
 * @author Leona Gottwald
 * @author Felipe Serrano
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip_cutsel.h"
#include "scip/scip_cut.h"
#include "scip/scip_lp.h"
#include "scip/scip_randnumgen.h"
#include "scip/cutsel_hybrid.h"


#define CUTSEL_NAME              "hybrid"
#define CUTSEL_DESC              "weighted sum of efficacy, dircutoffdist, objparal, and intsupport"
#define CUTSEL_PRIORITY           8000

#define RANDSEED                  0x5EED
#define GOODSCORE                 0.9
#define BADSCORE                  0.0

#define DEFAULT_EFFICACYWEIGHT          1.0  /**< weight of efficacy in score calculation */
#define DEFAULT_DIRCUTOFFDISTWEIGHT     0.0  /**< weight of directed cutoff distance in score calculation */
#define DEFAULT_OBJPARALWEIGHT          0.1  /**< weight of objective parallelism in score calculation */
#define DEFAULT_INTSUPPORTWEIGHT        0.1  /**< weight of integral support in cut score calculation */
#define DEFAULT_MINORTHO                0.90 /**< minimal orthogonality for a cut to enter the LP */
#define DEFAULT_MINORTHOROOT            0.90 /**< minimal orthogonality for a cut to enter the LP in the root node */

/*
 * Data structures
 */

/** cut selector data */
struct SCIP_CutselData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random generator for tiebreaking */
   SCIP_Real             goodscore;          /**< threshold for score of cut relative to best score to be considered good,
                                              *   so that less strict filtering is applied */
   SCIP_Real             badscore;           /**< threshold for score of cut relative to best score to be discarded */
   SCIP_Real             objparalweight;     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             efficacyweight;     /**< weight of efficacy in cut score calculation */
   SCIP_Real             dircutoffdistweight;/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             intsupportweight;   /**< weight of integral support in cut score calculation */
   SCIP_Real             minortho;           /**< minimal orthogonality for a cut to enter the LP */
   SCIP_Real             minorthoroot;       /**< minimal orthogonality for a cut to enter the LP in the root node */
};


/*
 * Local methods
 */

/** returns the maximum score of cuts; if scores is not NULL, then stores the individual score of each cut in scores */
static
SCIP_Real scoring(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to score */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for tie-breaking, or NULL */
   SCIP_Real             dircutoffdistweight,/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             efficacyweight,     /**< weight of efficacy in cut score calculation */
   SCIP_Real             objparalweight,     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             intsupportweight,   /**< weight of integral support in cut score calculation */
   int                   ncuts,              /**< number of cuts in cuts array */
   SCIP_Real*            scores              /**< array to store the score of cuts or NULL */
   )
{
   SCIP_Real maxscore = 0.0;
   SCIP_SOL* sol;
   int i;

   sol = SCIPgetBestSol(scip);

   /* if there is an incumbent and the factor is not 0.0, compute directed cutoff distances for the incumbent */
   if( sol != NULL && dircutoffdistweight > 0.0 )
   {
      for( i = 0; i < ncuts; ++i )
      {
         SCIP_Real score;
         SCIP_Real objparallelism;
         SCIP_Real intsupport;
         SCIP_Real efficacy;

         if( intsupportweight > 0.0 )
            intsupport = intsupportweight * SCIPgetRowNumIntCols(scip, cuts[i]) / (SCIP_Real) SCIProwGetNNonz(cuts[i]);
         else
            intsupport = 0.0;

         if( objparalweight > 0.0 )
            objparallelism = objparalweight * SCIPgetRowObjParallelism(scip, cuts[i]);
         else
            objparallelism = 0.0;

         efficacy = SCIPgetCutEfficacy(scip, NULL, cuts[i]);

         if( SCIProwIsLocal(cuts[i]) )
         {
            score = dircutoffdistweight * efficacy;
         }
         else
         {
            score = SCIPgetCutLPSolCutoffDistance(scip, sol, cuts[i]);
            score = dircutoffdistweight * MAX(score, efficacy);
         }

         efficacy *= efficacyweight;
         score += objparallelism + intsupport + efficacy;

         /* add small term to prefer global pool cuts */
         if( SCIProwIsInGlobalCutpool(cuts[i]) )
            score += 1e-4;

         if( randnumgen != NULL )
         {
            score += SCIPrandomGetReal(randnumgen, 0.0, 1e-6);
         }

         maxscore = MAX(maxscore, score);

         if( scores != NULL )
            scores[i] = score;
      }
   }
   else
   {
      /* in case there is no solution add the directed cutoff distance weight to the efficacy weight
       * since the efficacy underestimates the directed cuttoff distance
       */
      efficacyweight += dircutoffdistweight;
      for( i = 0; i < ncuts; ++i )
      {
         SCIP_Real score;
         SCIP_Real objparallelism;
         SCIP_Real intsupport;
         SCIP_Real efficacy;

         if( intsupportweight > 0.0 )
            intsupport = intsupportweight * SCIPgetRowNumIntCols(scip, cuts[i]) / (SCIP_Real) SCIProwGetNNonz(cuts[i]);
         else
            intsupport = 0.0;

         if( objparalweight > 0.0 )
            objparallelism = objparalweight * SCIPgetRowObjParallelism(scip, cuts[i]);
         else
            objparallelism = 0.0;

         efficacy = efficacyweight > 0.0 ?  efficacyweight * SCIPgetCutEfficacy(scip, NULL, cuts[i]) : 0.0;

         score = objparallelism + intsupport + efficacy;

         /* add small term to prefer global pool cuts */
         if( SCIProwIsInGlobalCutpool(cuts[i]) )
            score += 1e-4;

         if( randnumgen != NULL )
         {
            score += SCIPrandomGetReal(randnumgen, 0.0, 1e-6);
         }

         maxscore = MAX(maxscore, score);

         if( scores != NULL )
            scores[i] = score;
      }
   }
   return maxscore;
}


/** move the cut with the highest score to the first position in the array; there must be at least one cut */
static
void selectBestCut(
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_Real*            scores,             /**< array with scores of cuts to perform selection algorithm */
   int                   ncuts               /**< number of cuts in given array */
   )
{
   int i;
   int bestpos;
   SCIP_Real bestscore;

   assert(ncuts > 0);
   assert(cuts != NULL);
   assert(scores != NULL);

   bestscore = scores[0];
   bestpos = 0;

   for( i = 1; i < ncuts; ++i )
   {
      if( scores[i] > bestscore )
      {
         bestpos = i;
         bestscore = scores[i];
      }
   }

   SCIPswapPointers((void**) &cuts[bestpos], (void**) &cuts[0]);
   SCIPswapReals(&scores[bestpos], &scores[0]);
}

/** filters the given array of cuts to enforce a maximum parallelism constraint
 *  w.r.t the given cut; moves filtered cuts to the end of the array and returns number of selected cuts */
static
int filterWithParallelism(
   SCIP_ROW*             cut,                /**< cut to filter orthogonality with */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_Real*            scores,             /**< array with scores of cuts to perform selection algorithm */
   int                   ncuts,              /**< number of cuts in given array */
   SCIP_Real             goodscore,          /**< threshold for the score to be considered a good cut */
   SCIP_Real             goodmaxparall,      /**< maximal parallelism for good cuts */
   SCIP_Real             maxparall           /**< maximal parallelism for all cuts that are not good */
   )
{
   int i;

   assert( cut != NULL );
   assert( ncuts == 0 || cuts != NULL );
   assert( ncuts == 0 || scores != NULL );

   for( i = ncuts - 1; i >= 0; --i )
   {
      SCIP_Real thisparall;
      SCIP_Real thismaxparall;

      thisparall = SCIProwGetParallelism(cut, cuts[i], 'e');
      thismaxparall = scores[i] >= goodscore ? goodmaxparall : maxparall;

      if( thisparall > thismaxparall )
      {
         --ncuts;
         SCIPswapPointers((void**) &cuts[i], (void**) &cuts[ncuts]);
         SCIPswapReals(&scores[i], &scores[ncuts]);
      }
   }

   return ncuts;
}


/*
 * Callback methods of cut selector
 */


/** copy method for cut selector plugin (called when SCIP copies plugins) */
static
SCIP_DECL_CUTSELCOPY(cutselCopyHybrid)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(cutsel != NULL);
   assert(strcmp(SCIPcutselGetName(cutsel), CUTSEL_NAME) == 0);

   /* call inclusion method of cut selector */
   SCIP_CALL( SCIPincludeCutselHybrid(scip) );

   return SCIP_OKAY;
}

/** destructor of cut selector to free user data (called when SCIP is exiting) */
/**! [SnippetCutselFreeHybrid] */
static
SCIP_DECL_CUTSELFREE(cutselFreeHybrid)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);

   SCIPfreeBlockMemory(scip, &cutseldata);

   SCIPcutselSetData(cutsel, NULL);

   return SCIP_OKAY;
}
/**! [SnippetCutselFreeHybrid] */

/** initialization method of cut selector (called after problem was transformed) */
static
SCIP_DECL_CUTSELINIT(cutselInitHybrid)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &(cutseldata)->randnumgen, RANDSEED, TRUE) );

   return SCIP_OKAY;
}

/** deinitialization method of cut selector (called before transformed problem is freed) */
static
SCIP_DECL_CUTSELEXIT(cutselExitHybrid)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->randnumgen != NULL);

   SCIPfreeRandom(scip, &cutseldata->randnumgen);

   return SCIP_OKAY;
}

/** cut selection method of cut selector */
static
SCIP_DECL_CUTSELSELECT(cutselSelectHybrid)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;
   SCIP_Real goodmaxparall;
   SCIP_Real maxparall;

   assert(cutsel != NULL);
   assert(result != NULL);

   *result = SCIP_SUCCESS;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);

   if( root )
   {
      maxparall = 1.0 - cutseldata->minorthoroot;
      goodmaxparall = MAX(0.5, 1.0 - cutseldata->minorthoroot);
   }
   else
   {
      maxparall = 1.0 - cutseldata->minortho;
      goodmaxparall = MAX(0.5, 1.0 - cutseldata->minortho);
   }

   SCIP_CALL( SCIPselectCutsHybrid(scip, cuts, forcedcuts, cutseldata->randnumgen, cutseldata->goodscore, cutseldata->badscore,
         goodmaxparall, maxparall, cutseldata->dircutoffdistweight, cutseldata->efficacyweight,
         cutseldata->objparalweight, cutseldata->intsupportweight, ncuts, nforcedcuts, maxnselectedcuts, nselectedcuts) );

   return SCIP_OKAY;
}


/*
 * cut selector specific interface methods
 */

/** creates the hybrid cut selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeCutselHybrid(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CUTSELDATA* cutseldata;
   SCIP_CUTSEL* cutsel;

   /* create hybrid cut selector data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &cutseldata) );
   BMSclearMemory(cutseldata);
   cutseldata->goodscore = GOODSCORE;
   cutseldata->badscore  = BADSCORE;

   SCIP_CALL( SCIPincludeCutselBasic(scip, &cutsel, CUTSEL_NAME, CUTSEL_DESC, CUTSEL_PRIORITY, cutselSelectHybrid,
         cutseldata) );

   assert(cutsel != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetCutselCopy(scip, cutsel, cutselCopyHybrid) );

   SCIP_CALL( SCIPsetCutselFree(scip, cutsel, cutselFreeHybrid) );
   SCIP_CALL( SCIPsetCutselInit(scip, cutsel, cutselInitHybrid) );
   SCIP_CALL( SCIPsetCutselExit(scip, cutsel, cutselExitHybrid) );

   /* add hybrid cut selector parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "cutselection/" CUTSEL_NAME "/efficacyweight",
         "weight of efficacy in cut score calculation",
         &cutseldata->efficacyweight, FALSE, DEFAULT_EFFICACYWEIGHT, 0.0, SCIP_INVALID/10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "cutselection/" CUTSEL_NAME "/dircutoffdistweight",
         "weight of directed cutoff distance in cut score calculation",
         &cutseldata->dircutoffdistweight, FALSE, DEFAULT_DIRCUTOFFDISTWEIGHT, 0.0, SCIP_INVALID/10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "cutselection/" CUTSEL_NAME "/objparalweight",
         "weight of objective parallelism in cut score calculation",
         &cutseldata->objparalweight, FALSE, DEFAULT_OBJPARALWEIGHT, 0.0, SCIP_INVALID/10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "cutselection/" CUTSEL_NAME "/intsupportweight",
         "weight of integral support in cut score calculation",
         &cutseldata->intsupportweight, FALSE, DEFAULT_INTSUPPORTWEIGHT, 0.0, SCIP_INVALID/10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "cutselection/" CUTSEL_NAME "/minortho",
         "minimal orthogonality for a cut to enter the LP",
         &cutseldata->minortho, FALSE, DEFAULT_MINORTHO, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "cutselection/" CUTSEL_NAME "/minorthoroot",
         "minimal orthogonality for a cut to enter the LP in the root node",
         &cutseldata->minorthoroot, FALSE, DEFAULT_MINORTHOROOT, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}


/** perform a cut selection algorithm for the given array of cuts
 *
 *  This is the selection method of the hybrid cut selector which uses a weighted sum of the
 *  efficacy, parallelism, directed cutoff distance, and the integral support.
 *  The input cuts array gets resorted s.t the selected cuts come first and the remaining
 *  ones are the end.
 */
SCIP_RETCODE SCIPselectCutsHybrid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_ROW**            forcedcuts,         /**< array with forced cuts */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for tie-breaking, or NULL */
   SCIP_Real             goodscorefac,       /**< factor of best score among the given cuts to consider a cut good
                                              *   and filter with less strict settings of the maximum parallelism */
   SCIP_Real             badscorefac,        /**< factor of best score among the given cuts to consider a cut bad
                                              *   and discard it regardless of its parallelism to other cuts */
   SCIP_Real             goodmaxparall,      /**< maximum parallelism for good cuts */
   SCIP_Real             maxparall,          /**< maximum parallelism for non-good cuts */
   SCIP_Real             dircutoffdistweight,/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             efficacyweight,     /**< weight of efficacy in cut score calculation */
   SCIP_Real             objparalweight,     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             intsupportweight,   /**< weight of integral support in cut score calculation */
   int                   ncuts,              /**< number of cuts in cuts array */
   int                   nforcedcuts,        /**< number of forced cuts */
   int                   maxselectedcuts,    /**< maximal number of cuts from cuts array to select */
   int*                  nselectedcuts       /**< pointer to return number of selected cuts from cuts array */
   )
{
   SCIP_Real* scores;
   SCIP_Real* scoresptr;
   SCIP_Real maxforcedscores;
   SCIP_Real maxnonforcedscores;
   SCIP_Real goodscore;
   SCIP_Real badscore;
   int i;

   assert(cuts != NULL && ncuts > 0);
   assert(forcedcuts != NULL || nforcedcuts == 0);
   assert(nselectedcuts != NULL);

   *nselectedcuts = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &scores, ncuts) );

   /* compute scores of cuts and max score of cuts and forced cuts (used to define goodscore) */
   maxforcedscores = scoring(scip, forcedcuts, randnumgen, dircutoffdistweight, efficacyweight, objparalweight, intsupportweight, nforcedcuts, NULL);
   maxnonforcedscores = scoring(scip, cuts, randnumgen, dircutoffdistweight, efficacyweight, objparalweight, intsupportweight, ncuts, scores);

   goodscore = MAX(maxforcedscores, maxnonforcedscores);

   /* compute values for filtering cuts */
   badscore = goodscore * badscorefac;
   goodscore *= goodscorefac;

   /* perform cut selection algorithm for the cuts */

   /* forced cuts are going to be selected so use them to filter cuts */
   for( i = 0; i < nforcedcuts && ncuts > 0; ++i )
   {
      ncuts = filterWithParallelism(forcedcuts[i], cuts, scores, ncuts, goodscore, goodmaxparall, maxparall);
   }

   /* now greedily select the remaining cuts */
   scoresptr = scores;
   while( ncuts > 0 )
   {
      SCIP_ROW* selectedcut;

      selectBestCut(cuts, scores, ncuts);
      selectedcut = cuts[0];

      /* if the best cut of the remaining cuts is considered bad, we discard it and all remaining cuts */
      if( scores[0] < badscore )
         goto TERMINATE;

      ++(*nselectedcuts);

      /* if the maximal number of cuts was selected, we can stop here */
      if( *nselectedcuts == maxselectedcuts )
         goto TERMINATE;

      /* move the pointers to the next position and filter the remaining cuts to enforce the maximum parallelism constraint */
      ++cuts;
      ++scores;
      --ncuts;

      ncuts = filterWithParallelism(selectedcut, cuts, scores, ncuts, goodscore, goodmaxparall, maxparall);
   }

TERMINATE:
   SCIPfreeBufferArray(scip, &scoresptr);

   return SCIP_OKAY;
}
