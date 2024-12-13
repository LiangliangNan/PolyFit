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

/**@file   heur_farkasdiving.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  LP diving heuristic that tries to construct a Farkas-proof
 * @author Jakob Witzig
 *
 * The heuristic dives into the direction of the pseudosolution, i.e., variables get rounded
 * towards their best bound w.r.t there objective coefficient. This strategy is twofold, if
 * a feasible solution is found the solution has potentially a very good objective value; on the other
 * hand, the left-hand side of a potential Farkas-proof \f$y^Tb - y^TA{l',u'} > 0\f$ (i.e., infeasibility proof)
 * gets increased, where \f$l',u'\f$ are the local bounds. The contribution of each variable \f$x_i\f$ to the
 * Farkas-proof can be approximated by \f$c_i = y^TA_i\f$ because we only dive on basic variables with
 * reduced costs \f$c_i - y^TA_i = 0\f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/heur_farkasdiving.h"
#include "scip/heuristics.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_tree.h"
#include <string.h>

#define HEUR_NAME             "farkasdiving"
#define HEUR_DESC             "LP diving heuristic that tries to construct a Farkas-proof"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_DIVING
#define HEUR_PRIORITY         -900000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES     SCIP_DIVETYPE_INTEGRALITY | SCIP_DIVETYPE_SOS1VARIABLE /**< bit mask that represents all supported dive types */
#define DIVESET_ISPUBLIC      FALSE  /**< is this dive set publicly available (ie., can be used by other primal heuristics?) */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15 /**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           1 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_RANDSEED            151 /**< initial seed for random number generation */

#define DEFAULT_MAXOBJOCC           1.0 /**< maximal occurance factor of an objective coefficient */
#define DEFAULT_OBJDYN           0.0001 /**< minimal objective dynamism (log) */
#define DEFAULT_CHECKCANDS        FALSE /**< should diving candidates be checked before running? */
#define DEFAULT_SCALESCORE         TRUE /**< should the score be scaled? */
#define DEFAULT_ROOTSUCCESS        TRUE /**< should the heuristic only run within the tree if at least one solution
                                         *   was found at the root node? */
#define DEFAULT_SCALETYPE           'i' /**< scale score by [f]ractionality or [i]mpact on farkasproof */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             maxobjocc;          /**< maximal occurance factor of an objective coefficient */
   SCIP_Real             objdynamism;        /**< minimal objective dynamism (log) */
   SCIP_Bool             disabled;           /**< remember if the heuristic should not run at all */
   SCIP_Bool             glbchecked;         /**< remember whether one global check was performed */
   SCIP_Bool             checkcands;         /**< should diving candidates be checked before running? */
   SCIP_Bool             scalescore;         /**< should score be scaled by fractionality */
   SCIP_Bool             rootsuccess;        /**< run if successfull at root */
   SCIP_Bool             foundrootsol;       /**< was a solution found at the root node? */
   char                  scaletype;          /**< scale score by [f]ractionality or [i]mpact on farkasproof */
};


/** check whether the diving candidates fulfill requirements */
static
SCIP_RETCODE checkDivingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VAR**            divecandvars,       /**< diving candidates to check */
   int                   ndivecands,         /**< number of diving candidates */
   SCIP_Bool*            success             /**< pointer to store whether the check was successfull */
   )
{
   SCIP_Real* objcoefs;
   SCIP_Real lastobjcoef;
   SCIP_Real objdynamism;
   int maxfreq;
   int nnzobjcoefs;
#ifdef SCIP_DEBUG
   int ndiffnnzobjs;
#endif
   int i;

   *success = TRUE;

   assert(heurdata != NULL);
   assert(divecandvars != NULL);
   assert(ndivecands >= 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &objcoefs, ndivecands) );

   /* collect all absolute values of objective coefficients and count binary, integer,
    * and implicit integer in the objective function
    */
   nnzobjcoefs = 0;

   if( SCIPgetNObjVars(scip) > 0 )
   {
      for( i = 0; i < ndivecands; i++ )
      {
         SCIP_Real obj = SCIPvarGetObj(divecandvars[i]);

         if( SCIPisZero(scip, obj) )
            continue;

         objcoefs[nnzobjcoefs] = REALABS(obj);
         ++nnzobjcoefs;
      }
   }

   if( nnzobjcoefs == 0 )
   {
      *success = FALSE;
      goto TERMINATE;
   }
   assert(nnzobjcoefs > 0);

   /* skip here if we are cheching the global properties and want to check the local candidates, too */
   if( !heurdata->glbchecked && heurdata->checkcands )
      goto TERMINATE;

   /* sort in increasing order */
   SCIPsortReal(objcoefs, nnzobjcoefs);
   assert(!SCIPisZero(scip, objcoefs[0]));

   lastobjcoef = objcoefs[0];
#ifdef SCIP_DEBUG
   ndiffnnzobjs = 1;
#endif

   objdynamism = log10(objcoefs[nnzobjcoefs-1] / objcoefs[0]);

   SCIPdebugMsg(scip, "minabscoef: %g, maxabscoef: %g, objdym: %g\n", objcoefs[0], objcoefs[nnzobjcoefs-1], objdynamism);

   /* CHECK#2: check objective dynamism */
   if( objdynamism < heurdata->objdynamism )
   {
      SCIPdebugMsg(scip, " ---> disable farkasdiving at node %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

      *success = FALSE;
      goto TERMINATE;
   }

   /* CHECK#4: the check would be always fulfilled for heurdata->maxobjocc = 1.0 */
   if( heurdata->maxobjocc < 1.0 )
   {
      int tmpmaxfreq;

      tmpmaxfreq = 0;
      maxfreq = 0;

      /* count number of different absolute objective values */
      for( i = 1; i < nnzobjcoefs; i++ )
      {
         if( SCIPisGT(scip, objcoefs[i], lastobjcoef) )
         {
            if( tmpmaxfreq > maxfreq )
               maxfreq = tmpmaxfreq;
            tmpmaxfreq = 0;

            lastobjcoef = objcoefs[i];
#ifdef SCIP_DEBUG
            ++ndiffnnzobjs;
#endif
         }
         else
         {
            ++tmpmaxfreq;
         }
      }

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "%d divecands; %d nnzobjs; %d diffnnzobjs; %d maxfreq\n", ndivecands, nnzobjcoefs, ndiffnnzobjs,
         maxfreq);
#endif

      if( maxfreq > heurdata->maxobjocc * nnzobjcoefs )
      {
         SCIPdebugMsg(scip, " ---> disable farkasdiving at node %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

         *success = FALSE;
      }
   }

  TERMINATE:
   SCIPfreeBufferArray(scip, &objcoefs);

   return SCIP_OKAY;
}


/** check whether the objective functions has nonzero coefficients corresponding to binary and integer variables */
static
SCIP_RETCODE checkGlobalProperties(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_VAR** vars;
   SCIP_Bool success;
   int nbinvars;
   int nintvars;

   assert(scip != NULL);
   assert(heurdata != NULL);

   vars = SCIPgetVars(scip);
   nbinvars = SCIPgetNBinVars(scip);
   nintvars = SCIPgetNIntVars(scip);

   SCIP_CALL( checkDivingCandidates(scip, heurdata, vars, nbinvars+nintvars, &success) );

   if( !success )
   {
      SCIPdebugMsg(scip, " ---> disable farkasdiving (at least one global property is violated)\n");
      heurdata->disabled = TRUE;
   }

   heurdata->glbchecked = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyFarkasdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurFarkasdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeFarkasdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitFarkasdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   heurdata->disabled = FALSE;
   heurdata->glbchecked = FALSE;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitFarkasdiving)
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

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolFarkasdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->glbchecked = FALSE;
   heurdata->disabled = FALSE;
   heurdata->foundrootsol = FALSE;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecFarkasdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;
   SCIP_Bool success;

   heurdata = SCIPheurGetData(heur);
   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);

   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check some simple global properties that are needed to run this heuristic */
   if( !heurdata->glbchecked )
   {
      SCIP_CALL( checkGlobalProperties(scip, heurdata) );
   }

   /* terminate if the heuristic has been disabled */
   if( heurdata->disabled )
      return SCIP_OKAY;

   if( heurdata->rootsuccess && !heurdata->foundrootsol && SCIPgetDepth(scip) > 0 )
   {
      heurdata->disabled = TRUE;
      return SCIP_OKAY;
   }

   success = TRUE;

   /* check diving candidates in detail */
   if( heurdata->checkcands )
   {
      SCIP_VAR** divecandvars;
      int ndivecands;

      /* we can only access the branching candidates if the LP is solved to optimality */
      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIP_CALL( SCIPgetLPBranchCands(scip, &divecandvars, NULL, NULL, &ndivecands, NULL, NULL) );

         SCIP_CALL( checkDivingCandidates(scip, heurdata, divecandvars, ndivecands, &success) );
      }
      else
      {
         success = FALSE;
      }
   }

   if( success )
   {
      SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible, -1L, SCIP_DIVECONTEXT_SINGLE) );

      if( heurdata->rootsuccess && SCIPgetDepth(scip) == 0 && SCIPdivesetGetNSols(diveset, SCIP_DIVECONTEXT_SINGLE) > 0 )
         heurdata->foundrootsol = TRUE;
   }

   return SCIP_OKAY;
}

#define MIN_RAND 1e-06
#define MAX_RAND 1e-05

/** calculate score and preferred rounding direction for the candidate variable */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreFarkasdiving)
{  /*lint --e{715}*/
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real obj;

   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   randnumgen = SCIPdivesetGetRandnumgen(diveset);
   assert(randnumgen != NULL);

   obj = SCIPvarGetObj(cand);

   /* dive towards the pseudosolution, at the same time approximate the contribution to
    * a potentially Farkas-proof (infeasibility proof) by y^TA_i = c_i.
    */
   if( SCIPisNegative(scip, obj) )
   {
      *roundup = TRUE;
   }
   else if( SCIPisPositive(scip, obj) )
   {
      *roundup = FALSE;
   }
   else
   {
      if( SCIPisEQ(scip, candsfrac, 0.5) )
         *roundup = !SCIPrandomGetInt(randnumgen, 0, 1);
      else
         *roundup = (candsfrac > 0.5);
   }

   /* larger score is better */
   *score = REALABS(obj) + SCIPrandomGetReal(randnumgen, MIN_RAND, MAX_RAND);

   if( heurdata->scalescore )
   {
      if( heurdata->scaletype == 'f' )
      {
         if( *roundup )
            *score *= (1.0 - candsfrac);
         else
            *score *= candsfrac;
      }
      else
      {
         assert(heurdata->scaletype == 'i');

         if( *roundup )
            *score *= (SCIPceil(scip, candsol) - SCIPvarGetLbLocal(cand));
         else
            *score *= (SCIPvarGetUbLocal(cand) - SCIPfloor(scip, candsol));
      }
   }

   /* prefer decisions on binary variables */
   if( SCIPvarGetType(cand) != SCIP_VARTYPE_BINARY )
      *score = -1.0 / *score;

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */
#define divesetAvailableFarkasdiving NULL

/** creates the farkasdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurFarkasdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Farkasdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecFarkasdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyFarkasdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeFarkasdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitFarkasdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitFarkasdiving) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolFarkasdiving) );

   /* farkasdiving heuristic parameters */
   /* create a diveset (this will automatically install some additional parameters for the heuristic) */
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED, DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS, DIVESET_ISPUBLIC, DIVESET_DIVETYPES,
         divesetGetScoreFarkasdiving, divesetAvailableFarkasdiving) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/checkcands",
         "should diving candidates be checked before running?",
         &heurdata->checkcands, TRUE, DEFAULT_CHECKCANDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/scalescore",
         "should the score be scaled?",
         &heurdata->scalescore, TRUE, DEFAULT_SCALESCORE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/rootsuccess",
         "should the heuristic only run within the tree if at least one solution was found at the root node?",
         &heurdata->rootsuccess, TRUE, DEFAULT_ROOTSUCCESS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxobjocc",
         "maximal occurance factor of an objective coefficient",
         &heurdata->maxobjocc, TRUE, DEFAULT_MAXOBJOCC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/objdynamism",
         "minimal objective dynamism (log) to run",
         &heurdata->objdynamism, TRUE, DEFAULT_OBJDYN, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/scaletype",
         "scale score by [f]ractionality or [i]mpact on farkasproof",
         &heurdata->scaletype, TRUE, DEFAULT_SCALETYPE, "fi", NULL, NULL) );

   return SCIP_OKAY;
}
