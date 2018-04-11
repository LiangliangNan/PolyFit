/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heuristics.c
 * @brief  methods commonly used by primal heuristics
 * @author Gregor Hendel
 */

#include "scip/heuristics.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"

#include "pub_heur.h"

/* the indicator and SOS1 constraint handlers are included for the diving algorithm SCIPperformGenericDivingAlgorithm() */
#include "scip/cons_indicator.h"
#include "scip/cons_sos1.h"

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */


/** solve probing LP */
static
SCIP_RETCODE solveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_Longint          maxnlpiterations,   /**< maximum number of allowed LP iterations */
   SCIP_Bool*            lperror,            /**< pointer to store if an internal LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the LP was infeasible */
   )
{
   int lpiterationlimit;
   SCIP_RETCODE retstat;
   SCIP_Longint nlpiterations;

   assert(lperror != NULL);
   assert(cutoff != NULL);

   nlpiterations = SCIPgetNLPIterations(scip);

   /* allow at least MINLPITER more iterations */
   lpiterationlimit = (int)(maxnlpiterations - SCIPdivesetGetNLPIterations(diveset));
   lpiterationlimit = MAX(lpiterationlimit, MINLPITER);

   retstat = SCIPsolveProbingLP(scip, lpiterationlimit, lperror, cutoff);

   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
    * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
#ifdef NDEBUG
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while solving LP in %s diving heuristic; LP solve terminated with code <%d>.\n", SCIPdivesetGetName(diveset), retstat);
   }
#else
   SCIP_CALL( retstat );
#endif

   /* update iteration count */
   SCIPupdateDivesetLPStats(scip, diveset, SCIPgetNLPIterations(scip) - nlpiterations);

   return SCIP_OKAY;
}

/** select the next variable and type of diving */
static
SCIP_RETCODE selectNextDiving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< dive set */
   SCIP_SOL*             worksol,            /**< current working solution */
   SCIP_Bool             onlylpbranchcands,  /**< should only LP branching candidates be considered? */
   SCIP_Bool             storelpcandscores,  /**< should the scores of the LP candidates be updated? */
   SCIP_VAR**            lpcands,            /**< LP branching candidates, or NULL if not needed */
   SCIP_Real *           lpcandssol,         /**< solution values LP branching candidates, or NULL if not needed */
   SCIP_Real*            lpcandsfrac,        /**< fractionalities of LP branching candidates, or NULL if not needed*/
   SCIP_Real*            lpcandsscores,      /**< array with LP branching candidate scores, or NULL */
   SCIP_Bool*            lpcandroundup,      /**< array to remember whether the preferred branching direction is upwards */
   int*                  nviollpcands,       /**< pointer to store the number of LP candidates whose solution value already violates local bounds */
   int                   nlpcands,           /**< number of current LP cands */
   SCIP_Bool*            enfosuccess,        /**< pointer to store whether a candidate was sucessfully found */
   SCIP_Bool*            infeasible          /**< pointer to store whether the diving can be immediately aborted because it is infeasible */
   )
{
   assert(scip != NULL);
   assert(worksol != NULL);
   assert(!onlylpbranchcands || lpcandsscores != NULL);
   assert(!onlylpbranchcands || lpcandroundup != NULL);
   assert(enfosuccess != NULL);
   assert(infeasible != NULL);
   assert(nviollpcands != NULL);

   *nviollpcands = 0;
   /* we use diving solution enforcement provided by the constraint handlers */
   if( !onlylpbranchcands )
   {
      SCIP_CALL( SCIPgetDiveBoundChanges(scip, diveset, worksol, enfosuccess, infeasible) );
   }
   else
   {
      int c;
      int bestcandidx;
      SCIP_Real bestscore;
      SCIP_Real score;

      bestscore = SCIP_REAL_MIN;
      bestcandidx = -1;

      SCIPclearDiveBoundChanges(scip);

      /* search for the candidate that maximizes the dive set score function and whose solution value is still feasible */
      for( c = 0; c < nlpcands; ++c )
      {
         assert(SCIPgetSolVal(scip, worksol, lpcands[c]) == lpcandssol[c]); /*lint !e777 doesn't like comparing floats for equality */

         /* scores are kept in arrays for faster reuse */
         if( storelpcandscores )
         {
            SCIP_CALL( SCIPgetDivesetScore(scip, diveset, SCIP_DIVETYPE_INTEGRALITY, lpcands[c], lpcandssol[c],
                  lpcandsfrac[c], &lpcandsscores[c], &lpcandroundup[c]) );
         }

         score = lpcandsscores[c];
         /* update the best candidate if it has a higher score and a solution value which does not violate one of the local bounds */
         if( SCIPisLE(scip, SCIPvarGetLbLocal(lpcands[c]), lpcandssol[c]) && SCIPisGE(scip, SCIPvarGetUbLocal(lpcands[c]), lpcandssol[c]) )
         {
            if( score > bestscore )
            {
               bestcandidx = c;
               bestscore = score;
            }
         }
         else
            ++(*nviollpcands);
      }

      /* there is no guarantee that a candidate is found since local bounds might render all solution values infeasible */
      *enfosuccess = (bestcandidx >= 0);
      if( *enfosuccess )
      {
         /* if we want to round up the best candidate, it is added as the preferred bound change */
         SCIP_CALL( SCIPaddDiveBoundChange(scip, lpcands[bestcandidx], SCIP_BRANCHDIR_UPWARDS,
               SCIPceil(scip, lpcandssol[bestcandidx]), lpcandroundup[bestcandidx]) );
         SCIP_CALL( SCIPaddDiveBoundChange(scip, lpcands[bestcandidx], SCIP_BRANCHDIR_DOWNWARDS,
               SCIPfloor(scip, lpcandssol[bestcandidx]), ! lpcandroundup[bestcandidx]) );
      }
   }
   return SCIP_OKAY;
}



/** performs a diving within the limits of the diveset parameters
 *
 *  This method performs a diving according to the settings defined by the diving settings @p diveset; Contrary to the
 *  name, SCIP enters probing mode (not diving mode) and dives along a path into the tree. Domain propagation
 *  is applied at every node in the tree, whereas probing LPs might be solved less frequently.
 *
 *  Starting from the current LP solution, the algorithm selects candidates which maximize the
 *  score defined by the @p diveset and whose solution value has not yet been rendered infeasible by propagation,
 *  and propagates the bound change on this candidate.
 *
 *  The algorithm iteratively selects the the next (unfixed) candidate in the list, until either enough domain changes
 *  or the resolve frequency of the LP trigger an LP resolve (and hence, the set of potential candidates changes),
 *  or the last node is proven to be infeasible. It optionally backtracks and tries the
 *  other branching direction.
 *
 *  After the set of remaining candidates is empty or the targeted depth is reached, the node LP is
 *  solved, and the old candidates are replaced by the new LP candidates.
 *
 *  @see heur_guideddiving.c for an example implementation of a dive set controlling the diving algorithm.
 *
 *  @note the node from where the algorithm is called is checked for a basic LP solution. If the solution
 *        is non-basic, e.g., when barrier without crossover is used, the method returns without performing a dive.
 *
 *  @note currently, when multiple diving heuristics call this method and solve an LP at the same node, only the first
 *        call will be executed, @see SCIPgetLastDiveNode()
 *
 *  @todo generalize method to work correctly with pseudo or external branching/diving candidates
 */
SCIP_RETCODE SCIPperformGenericDivingAlgorithm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< settings for diving */
   SCIP_SOL*             worksol,            /**< non-NULL working solution */
   SCIP_HEUR*            heur,               /**< the calling primal heuristic */
   SCIP_RESULT*          result,             /**< SCIP result pointer */
   SCIP_Bool             nodeinfeasible      /**< is the current node known to be infeasible? */
   )
{
   SCIP_CONSHDLR* indconshdlr;               /* constraint handler for indicator constraints */
   SCIP_CONSHDLR* sos1conshdlr;              /* constraint handler for SOS1 constraints */
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;

   SCIP_VAR** previouscands;
   SCIP_Real* lpcandsscores;
   SCIP_Real* previousvals;
   SCIP_Real* lpcandsfrac;
   SCIP_Bool* lpcandroundup;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real ubquot;
   SCIP_Real avgquot;
   SCIP_Longint ncalls;
   SCIP_Longint oldsolsuccess;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   SCIP_Longint domreds;
   int startndivecands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int totalnbacktracks;
   int totalnprobingnodes;
   int lastlpdepth;
   int previouscandssize;
   int lpcandsscoressize;
   int nviollpcands;
   SCIP_Longint oldnsolsfound;
   SCIP_Longint oldnbestsolsfound;

   SCIP_Bool success;
   SCIP_Bool leafsol;
   SCIP_Bool enfosuccess;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Bool backtrack;
   SCIP_Bool onlylpbranchcands;

   int nlpcands;
   int lpsolvefreq;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);
   assert(worksol != NULL);

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
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < SCIPdivesetGetMinRelDepth(diveset) * maxdepth || depth > SCIPdivesetGetMaxRelDepth(diveset) * maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPdivesetGetNCalls(diveset);
   oldsolsuccess = SCIPdivesetGetSolSuccess(diveset);

   /*todo another factor of 10, REALLY? */
   maxnlpiterations = (SCIP_Longint)((1.0 + 10*(oldsolsuccess+1.0)/(ncalls+1.0)) * SCIPdivesetGetMaxLPIterQuot(diveset) * nlpiterations);
   maxnlpiterations += SCIPdivesetGetMaxLPIterOffset(diveset);

   /* don't try to dive, if we took too many LP iterations during diving */
   if( SCIPdivesetGetNLPIterations(diveset) >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   if( SCIPdivesetGetNLPIterations(diveset) + MINLPITER > maxnlpiterations )
      maxnlpiterations = SCIPdivesetGetNLPIterations(diveset) + MINLPITER;

   /* these constraint handlers are required for polishing an LP relaxation solution beyond rounding */
   indconshdlr = SCIPfindConshdlr(scip, "indicator");
   sos1conshdlr = SCIPfindConshdlr(scip, "SOS1");

   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );

   onlylpbranchcands = SCIPdivesetUseOnlyLPBranchcands(diveset);
   /* don't try to dive, if there are no diving candidates */
   if( onlylpbranchcands && nlpcands == 0 )
      return SCIP_OKAY;

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      ubquot = SCIPdivesetGetUbQuotNoSol(diveset);
      avgquot = SCIPdivesetGetAvgQuotNoSol(diveset);
   }
   else
   {
      ubquot = SCIPdivesetGetUbQuot(diveset);
      avgquot = SCIPdivesetGetAvgQuot(diveset);
   }

   if( ubquot > 0.0 )
      searchubbound = SCIPgetLowerbound(scip) + ubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
   else
      searchubbound = SCIPinfinity(scip);

   if( avgquot > 0.0 )
      searchavgbound = SCIPgetLowerbound(scip) + avgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   else
      searchavgbound = SCIPinfinity(scip);

   searchbound = MIN(searchubbound, searchavgbound);

   if( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if ( sos1conshdlr != NULL )
      maxdivedepth += SCIPgetNSOS1Vars(sos1conshdlr);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;

   *result = SCIP_DIDNOTFIND;

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   SCIPdebugMsg(scip, "(node %" SCIP_LONGINT_FORMAT ") executing %s heuristic: depth=%d, %d fractionals, dualbound=%g, avgbound=%g, cutoffbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPheurGetName(heur), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPgetAvgDualbound(scip),
      SCIPretransformObj(scip, SCIPgetCutoffbound(scip)), SCIPretransformObj(scip, searchbound));


   /* allocate buffer storage for previous candidates and their branching values for pseudo cost updates */
   lpsolvefreq = SCIPdivesetGetLPSolveFreq(diveset);
   previouscandssize = MAX(1, lpsolvefreq);
   SCIP_CALL( SCIPallocBufferArray(scip, &previouscands, previouscandssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &previousvals, previouscandssize) );

   /* keep some statistics */
   lperror = FALSE;
   cutoff = FALSE;
   lastlpdepth = -1;
   startndivecands = nlpcands;
   totalnbacktracks = 0;
   totalnprobingnodes = 0;
   oldnsolsfound = SCIPgetNSolsFound(scip);
   oldnbestsolsfound = SCIPgetNBestSolsFound(scip);

   /* link the working solution to the dive set */
   SCIPdivesetSetWorkSolution(diveset, worksol);

   if( onlylpbranchcands )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lpcandsscores, nlpcands) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lpcandroundup, nlpcands) );

      lpcandsscoressize = nlpcands;
   }
   else
   {
      lpcandroundup = NULL;
      lpcandsscores = NULL;
      lpcandsscoressize = 0;
   }

   enfosuccess = TRUE;
   leafsol = FALSE;

   /* LP loop; every time a new LP was solved, conditions are checked
    * dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   while( !lperror && !cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL && enfosuccess
      && (SCIPgetProbingDepth(scip) < 10
         || nlpcands <= startndivecands - SCIPgetProbingDepth(scip) / 2
         || (SCIPgetProbingDepth(scip) < maxdivedepth && SCIPdivesetGetNLPIterations(diveset) < maxnlpiterations && SCIPgetLPObjval(scip) < searchbound))
         && !SCIPisStopped(scip) )
   {
      SCIP_Real lastlpobjval;
      int c;
      SCIP_Bool allroundable;
      SCIP_Bool infeasible;
      int prevcandsinsertidx;

      /* remember the last LP depth  */
      assert(lastlpdepth < SCIPgetProbingDepth(scip));
      lastlpdepth = SCIPgetProbingDepth(scip);
      domreds = 0;

      SCIPdebugMsg(scip, "%s heuristic continues diving at depth %d, %d candidates left\n",
         SCIPdivesetGetName(diveset), lastlpdepth, nlpcands);


      /* loop over candidates and determine if they are roundable */
      allroundable = TRUE;
      c = 0;
      while( allroundable && c < nlpcands )
      {
         if( SCIPvarMayRoundDown(lpcands[c]) || SCIPvarMayRoundUp(lpcands[c]) )
            allroundable = TRUE;
         else
            allroundable = FALSE;
         ++c;
      }

      /* if all candidates are roundable, try to round the solution */
      if( allroundable )
      {
         success = FALSE;

         /* working solution must be linked to LP solution */
         SCIP_CALL( SCIPlinkLPSol(scip, worksol) );
         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIProundSol(scip, worksol, &success) );

         /* adjust SOS1 constraints */
         if( success && sos1conshdlr != NULL )
         {
            SCIP_Bool changed;
            SCIP_CALL( SCIPmakeSOS1sFeasible(scip, sos1conshdlr, worksol, &changed, &success) );
         }

         /* succesfully rounded solutions are tried for primal feasibility */
         if( success )
         {
            SCIP_Bool changed = FALSE;
            SCIPdebugMsg(scip, "%s found roundable primal solution: obj=%g\n", SCIPdivesetGetName(diveset), SCIPgetSolOrigObj(scip, worksol));

            /* adjust indicator constraints */
            if( indconshdlr != NULL )
            {
               SCIP_CALL( SCIPmakeIndicatorsFeasible(scip, indconshdlr, worksol, &changed) );
            }

            success = FALSE;
            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, changed, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMsg(scip, " -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
               leafsol = (nlpcands == 0);

               /* the rounded solution found above led to a cutoff of the node LP solution */
               if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OBJLIMIT )
               {
                  cutoff = TRUE;
                  break;
               }
            }
         }
      }

      /* working solution must be linked to LP solution */
      assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);
      lastlpobjval = SCIPgetLPObjval(scip);
      SCIP_CALL( SCIPlinkLPSol(scip, worksol) );

      /* in case we do not solve LP's at every probing node, we must explicitly store the solution values by unlinking the
       * solution, otherwise, the working solution may contain wrong entries, if, e.g., a backtrack occurred after an
       * intermediate LP resolve
       */
      if( lpsolvefreq != 1 )
      {
         SCIP_CALL( SCIPunlinkSol(scip, worksol) );
      }


      /* ensure array sizes for the diving on the fractional variables */
      if( onlylpbranchcands && nlpcands > lpcandsscoressize )
      {
         assert(nlpcands > 0);
         assert(lpcandsscores != NULL);
         assert(lpcandroundup != NULL);

         SCIP_CALL( SCIPreallocBufferArray(scip, &lpcandsscores, nlpcands) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &lpcandroundup, nlpcands) );

         lpcandsscoressize = nlpcands;
      }


      enfosuccess = FALSE;
      /* select the next diving action by selecting appropriate dive bound changes for the preferred and alternative child */
      SCIP_CALL( selectNextDiving(scip, diveset, worksol, onlylpbranchcands, SCIPgetProbingDepth(scip) == lastlpdepth,
             lpcands, lpcandssol, lpcandsfrac, lpcandsscores, lpcandroundup, &nviollpcands, nlpcands,
             &enfosuccess, &infeasible) );

      /* if we did not succeed finding an enforcement, the solution is potentially feasible and we break immediately */
      if( ! enfosuccess )
         break;

      prevcandsinsertidx = -1;

      /* start propagating candidate variables
       *   - until the desired targetdepth is reached,
       *   - or there is no further candidate variable left because of intermediate bound changes,
       *   - or a cutoff is detected
       */
      do
      {
         SCIP_VAR* bdchgvar;
         SCIP_Real bdchgvalue;
         SCIP_Longint localdomreds;
         SCIP_BRANCHDIR bdchgdir;
         int nbdchanges;

         /* ensure that a new candidate was successfully determined (usually at the end of the previous loop iteration) */
         assert(enfosuccess);
         bdchgvar = NULL;
         bdchgvalue = SCIP_INVALID;
         bdchgdir = SCIP_BRANCHDIR_AUTO;

         backtracked = FALSE;
         do
         {
            int d;
            SCIP_VAR** bdchgvars;
            SCIP_BRANCHDIR* bdchgdirs;
            SCIP_Real* bdchgvals;

            nbdchanges = 0;
            /* get the bound change information stored in the dive set */
            SCIPgetDiveBoundChangeData(scip, &bdchgvars, &bdchgdirs, &bdchgvals, &nbdchanges, !backtracked);

            assert(nbdchanges > 0);
            assert(bdchgvars != NULL);
            assert(bdchgdirs != NULL);
            assert(bdchgvals != NULL);

            /* return if we reached the depth limit */
            if( SCIP_MAXTREEDEPTH <= SCIPgetDepth(scip) )
            {
               SCIPdebugMsg(scip, "depth limit reached, we stop diving immediately.\n");
               goto TERMINATE;
            }

            /* dive deeper into the tree */
            SCIP_CALL( SCIPnewProbingNode(scip) );
            ++totalnprobingnodes;

            /* apply all suggested domain changes of the variables */
            for( d = 0; d < nbdchanges; ++d )
            {
               SCIP_Real lblocal;
               SCIP_Real ublocal;
               SCIP_Bool infeasbdchange;

               /* get the next bound change data: variable, direction, and value */
               bdchgvar = bdchgvars[d];
               bdchgvalue = bdchgvals[d];
               bdchgdir = bdchgdirs[d];

               assert(bdchgvar != NULL);
               lblocal = SCIPvarGetLbLocal(bdchgvar);
               ublocal = SCIPvarGetUbLocal(bdchgvar);


               SCIPdebugMsg(scip, "  dive %d/%d, LP iter %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ": var <%s>, oldbounds=[%g,%g],",
                     SCIPgetProbingDepth(scip), maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
                     SCIPvarGetName(bdchgvar), lblocal, ublocal);

               infeasbdchange = FALSE;
               /* tighten the lower and/or upper bound depending on the bound change type */
               switch( bdchgdir )
               {
                  case SCIP_BRANCHDIR_UPWARDS:
                     /* test if bound change is possible, otherwise propagation might have deduced the same
                      * bound already or numerical troubles might have occurred */
                     if( SCIPisFeasLE(scip, bdchgvalue, lblocal) || SCIPisFeasGT(scip, bdchgvalue, ublocal) )
                        infeasbdchange = TRUE;
                     else
                     {
                        /* round variable up */
                        SCIP_CALL( SCIPchgVarLbProbing(scip, bdchgvar, bdchgvalue) );
                     }
                     break;
                  case SCIP_BRANCHDIR_DOWNWARDS:
                     /* test if bound change is possible, otherwise propagation might have deduced the same
                      * bound already or numerical troubles might have occurred */
                     if( SCIPisFeasGE(scip, bdchgvalue, ublocal) || SCIPisFeasLT(scip, bdchgvalue, lblocal) )
                        infeasbdchange = TRUE;
                     else
                     {
                        /* round variable down */
                        SCIP_CALL( SCIPchgVarUbProbing(scip, bdchgvar, bdchgvalue) );
                     }
                     break;
                  case SCIP_BRANCHDIR_FIXED:
                     /* test if bound change is possible, otherwise propagation might have deduced the same
                      * bound already or numerical troubles might have occurred */
                     if( SCIPisFeasLT(scip, bdchgvalue, lblocal) || SCIPisFeasGT(scip, bdchgvalue, ublocal) ||
                           (SCIPisFeasEQ(scip, lblocal, ublocal) && nbdchanges < 2) )
                        infeasbdchange = TRUE;
                     else
                     {
                        /* fix variable to the bound change value */
                        if( SCIPisFeasLT(scip, lblocal, bdchgvalue) )
                        {
                           SCIP_CALL( SCIPchgVarLbProbing(scip, bdchgvar, bdchgvalue) );
                        }
                        if( SCIPisFeasGT(scip, ublocal, bdchgvalue) )
                        {
                           SCIP_CALL( SCIPchgVarUbProbing(scip, bdchgvar, bdchgvalue) );
                        }
                     }
                     break;
                  case SCIP_BRANCHDIR_AUTO:
                  default:
                     SCIPerrorMessage("Error: Unsupported bound change direction <%d> specified for diving, aborting\n",bdchgdirs[d]);
                     SCIPABORT();
                     return SCIP_INVALIDDATA; /*lint !e527*/
               }
               /* if the variable domain has been shrunk in the meantime, numerical troubles may have
                * occured or variable was fixed by propagation while backtracking => Abort diving!
                */
               if( infeasbdchange )
               {
                  SCIPdebugMsg(scip, "\nSelected variable <%s> domain already [%g,%g] as least as tight as desired bound change, diving aborted \n",
                     SCIPvarGetName(bdchgvar), SCIPvarGetLbLocal(bdchgvar), SCIPvarGetUbLocal(bdchgvar));
                  cutoff = TRUE;
                  break;
               }

               SCIPdebugMsg(scip, "newbounds=[%g,%g]\n", SCIPvarGetLbLocal(bdchgvar), SCIPvarGetUbLocal(bdchgvar));
            }
            /* break loop immediately if we detected a cutoff */
            if( cutoff )
               break;

            /* apply domain propagation */
            localdomreds = 0;
            SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, &localdomreds) );

            /* add the number of bound changes we applied by ourselves after propagation, otherwise the counter would have been reset */
            localdomreds += nbdchanges;

            /* resolve the diving LP if the diving resolve frequency is reached or a sufficient number of intermediate bound changes
             * was reached
             */
            if( ! cutoff
                  && ((lpsolvefreq > 0 && ((SCIPgetProbingDepth(scip) - lastlpdepth) % lpsolvefreq) == 0)
                  || (domreds + localdomreds > SCIPdivesetGetLPResolveDomChgQuot(diveset) * SCIPgetNVars(scip))
                  || (onlylpbranchcands && nviollpcands > (int)(SCIPdivesetGetLPResolveDomChgQuot(diveset) * nlpcands))) )
            {
               SCIP_CALL( solveLP(scip, diveset, maxnlpiterations, &lperror, &cutoff) );

               /* lp errors lead to early termination */
               if( lperror )
               {
                  cutoff = TRUE;
                  break;
               }
            }

            /* perform backtracking if a cutoff was detected */
            if( cutoff && !backtracked && SCIPdivesetUseBacktrack(diveset) )
            {
               SCIPdebugMsg(scip, "  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
               SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );
               ++totalnbacktracks;
               backtracked = TRUE;
               backtrack = TRUE;
               cutoff = FALSE;
            }
            else
               backtrack = FALSE;
         }
         while( backtrack );

         /* we add the domain reductions from the last evaluated node */
         domreds += localdomreds; /*lint !e771 lint thinks localdomreds has not been initialized */

         /* store candidate for pseudo cost update and choose next candidate only if no cutoff was detected */
         if( ! cutoff )
         {
            if( nbdchanges == 1 && (bdchgdir == SCIP_BRANCHDIR_UPWARDS || bdchgdir == SCIP_BRANCHDIR_DOWNWARDS) )
            {
               ++prevcandsinsertidx;
               assert(prevcandsinsertidx <= SCIPgetProbingDepth(scip) - lastlpdepth - 1);
               assert(SCIPgetProbingDepth(scip) > 0);
               assert(bdchgvar != NULL);
               assert(bdchgvalue != SCIP_INVALID); /*lint !e777 doesn't like comparing floats for equality */

               /* extend array in case of a dynamic, domain change based LP resolve strategy */
               if( prevcandsinsertidx >= previouscandssize )
               {
                  previouscandssize *= 2;
                  SCIP_CALL( SCIPreallocBufferArray(scip, &previouscands, previouscandssize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &previousvals, previouscandssize) );
               }
               assert(previouscandssize > prevcandsinsertidx);

               /* store candidate for pseudo cost update */
               previouscands[prevcandsinsertidx] = bdchgvar;
               previousvals[prevcandsinsertidx] = bdchgvalue;
            }

            /* choose next candidate variable and resolve the LP if none is found. */
            if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED )
            {
               assert(SCIPgetProbingDepth(scip) > lastlpdepth);
               enfosuccess = FALSE;

               /* select the next diving action */
               SCIP_CALL( selectNextDiving(scip, diveset, worksol, onlylpbranchcands, SCIPgetProbingDepth(scip) == lastlpdepth,
                      lpcands, lpcandssol, lpcandsfrac, lpcandsscores, lpcandroundup, &nviollpcands, nlpcands,
                      &enfosuccess, &infeasible) );


               /* in case of an unsuccesful candidate search, we solve the node LP */
               if( !enfosuccess )
               {
                  SCIP_CALL( solveLP(scip, diveset, maxnlpiterations, &lperror, &cutoff) );

                  /* check for an LP error and terminate in this case, cutoffs lead to termination anyway */
                  if( lperror )
                     cutoff = TRUE;

                  /* enfosuccess must be set to TRUE for entering the main LP loop again */
                  enfosuccess = TRUE;
               }
            }
         }
      }
      while( !cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED );

      assert(cutoff || lperror || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_NOTSOLVED);

      assert(cutoff || (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OBJLIMIT && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_INFEASIBLE &&
            (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL || SCIPisLT(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))));


      /* check new LP candidates and use the LP Objective gain to update pseudo cost information */
      if( ! cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         int v;
         SCIP_Real gain;

         SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );

         SCIPdebugMsg(scip, "   -> lpsolstat=%d, objval=%g/%g, nfrac=%d\n", SCIPgetLPSolstat(scip), SCIPgetLPObjval(scip), searchbound, nlpcands);
         /* distribute the gain equally over all variables that we rounded since the last LP */
         gain = SCIPgetLPObjval(scip) - lastlpobjval;
         gain = MAX(gain, 0.0);
         gain /= (1.0 * (SCIPgetProbingDepth(scip) - lastlpdepth));

         /* loop over previously fixed candidates and share gain improvement */
         for( v = 0; v <= prevcandsinsertidx; ++v )
         {
            SCIP_VAR* cand = previouscands[v];
            SCIP_Real val = previousvals[v];
            SCIP_Real solval = SCIPgetSolVal(scip, worksol, cand);

            /* todo: should the gain be shared with a smaller weight, instead of dividing the gain itself? */
            /* it may happen that a variable had an integral solution value beforehand, e.g., for indicator variables */
            if( ! SCIPisZero(scip, val - solval) )
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, cand, val - solval, gain, 1.0) );
            }
         }
      }
      else
         nlpcands = 0;
   }

   success = FALSE;
   /* check if a solution has been found */
   if( !enfosuccess && !lperror && !cutoff && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, worksol) );
      SCIPdebugMsg(scip, "%s found primal solution: obj=%g\n", SCIPdivesetGetName(diveset), SCIPgetSolOrigObj(scip, worksol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, worksol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMsg(scip, " -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
         leafsol = TRUE;
      }
   }

   SCIPupdateDivesetStats(scip, diveset, totalnprobingnodes, totalnbacktracks, SCIPgetNSolsFound(scip) - oldnsolsfound,
         SCIPgetNBestSolsFound(scip) - oldnbestsolsfound, leafsol);

   SCIPdebugMsg(scip, "(node %" SCIP_LONGINT_FORMAT ") finished %s heuristic: %d fractionals, dive %d/%d, LP iter %" SCIP_LONGINT_FORMAT "/%" SCIP_LONGINT_FORMAT ", objval=%g/%g, lpsolstat=%d, cutoff=%u\n",
      SCIPgetNNodes(scip), SCIPdivesetGetName(diveset), nlpcands, SCIPgetProbingDepth(scip), maxdivedepth, SCIPdivesetGetNLPIterations(diveset), maxnlpiterations,
      SCIPretransformObj(scip, SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL ? SCIPgetLPObjval(scip) : SCIPinfinity(scip)), SCIPretransformObj(scip, searchbound), SCIPgetLPSolstat(scip), cutoff);

  TERMINATE:
   /* end probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   SCIPdivesetSetWorkSolution(diveset, NULL);

   if( onlylpbranchcands )
   {
      SCIPfreeBufferArray(scip, &lpcandroundup);
      SCIPfreeBufferArray(scip, &lpcandsscores);
   }
   SCIPfreeBufferArray(scip, &previousvals);
   SCIPfreeBufferArray(scip, &previouscands);

   return SCIP_OKAY;
}


/** creates the rows of the subproblem */
static
SCIP_RETCODE createRows(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem */
   SCIP_HASHMAP*         varmap              /**< a hashmap to store the mapping of source variables to the corresponding
                                               *   target variables */
   )
{
   SCIP_ROW** rows;                          /* original scip rows                       */
   SCIP_CONS* cons;                          /* new constraint                           */
   SCIP_VAR** consvars;                      /* new constraint's variables               */
   SCIP_COL** cols;                          /* original row's columns                   */

   SCIP_Real constant;                       /* constant added to the row                */
   SCIP_Real lhs;                            /* left hand side of the row                */
   SCIP_Real rhs;                            /* left right side of the row               */
   SCIP_Real* vals;                          /* variables' coefficient values of the row */

   int nrows;
   int nnonz;
   int i;
   int j;

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      /* ignore rows that are only locally valid */
      if( SCIProwIsLocal(rows[i]) )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, (SCIPcolGetVar(cols[j])));

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}


/** get a sub-SCIP copy of the transformed problem */
SCIP_RETCODE SCIPcopyLargeNeighborhoodSearch(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP used by the heuristic */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables */
   const char*           suffix,             /**< suffix for the problem name */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             uselprows,          /**< should the linear relaxation of the problem defined by LP rows be copied? */
   SCIP_Bool             copycuts,           /**< should cuts be copied (only if uselprows == FALSE) */
   SCIP_Bool*            success,            /**< was the copying successful? */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid, or NULL */
   )
{
   assert(sourcescip != NULL);
   assert(suffix != NULL);
   assert(subscip != NULL);
   assert(varmap != NULL);
   assert(success != NULL);

   if( uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_crossoversub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_%s", SCIPgetProbName(sourcescip), suffix);

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(sourcescip, subscip, varmap, NULL, fixedvars, fixedvals, nfixedvars, TRUE) );

      /* copy parameter settings */
      SCIP_CALL( SCIPcopyParamSettings(sourcescip, subscip) );

      /* create linear constraints from LP rows of the source problem */
      SCIP_CALL( createRows(sourcescip, subscip, varmap) );
   }
   else
   {
      SCIP_CALL( SCIPcopyConsCompression(sourcescip, subscip, varmap, NULL, suffix, fixedvars, fixedvals, nfixedvars,
            TRUE, FALSE, TRUE, valid) );

      if( copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(sourcescip, subscip, varmap, NULL, TRUE, NULL) );
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}
