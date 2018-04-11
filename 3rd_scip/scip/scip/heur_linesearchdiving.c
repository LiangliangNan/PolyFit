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

/**@file   heur_linesearchdiving.c
 * @brief  LP diving heuristic that fixes variables with a large difference to their root solution
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_linesearchdiving.h"

#define HEUR_NAME             "linesearchdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings following the line from root solution to current solution"
#define HEUR_DISPCHAR         'l'
#define HEUR_PRIORITY         -1006000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          6
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES     SCIP_DIVETYPE_INTEGRALITY | SCIP_DIVETYPE_SOS1VARIABLE /**< bit mask that represents all supported dive types */


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
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_RANDSEED            137 /**< default initialization for random seed number generation */
/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
};


/*
 * Local methods
 */


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLinesearchdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurLinesearchdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLinesearchdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
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
SCIP_DECL_HEURINIT(heurInitLinesearchdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitLinesearchdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLinesearchdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);
   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}

/* diving setting callbacks */

/** returns a score for the given candidate -- the best candidate maximizes the diving score */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreLinesearchdiving)
{  /*lint --e{715}*/
   SCIP_Real rootsolval;
   SCIP_Real distquot;

   rootsolval = SCIPvarGetRootSol(cand);

   /* preferred branching direction is further away from the root LP solution */
   if( SCIPisLT(scip, candsol, rootsolval) )
   {
      /* round down*/
      *roundup = FALSE;

      switch( divetype )
      {
         case SCIP_DIVETYPE_INTEGRALITY:
            distquot = (candsfrac + SCIPsumepsilon(scip)) / (rootsolval - candsol);
            break;
         case SCIP_DIVETYPE_SOS1VARIABLE:
            if ( SCIPisFeasPositive(scip, candsol) )
               distquot = (candsfrac + SCIPsumepsilon(scip)) / (rootsolval - candsol);
            else
               distquot = (1.0 - candsfrac + SCIPsumepsilon(scip)) / (candsol - rootsolval);
            break;
         default:
            SCIPerrorMessage("Error: Unsupported diving type\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
      } /*lint !e788*/

      /* avoid roundable candidates */
      if( SCIPvarMayRoundDown(cand) )
         distquot *= 1000.0;
   }
   else if( SCIPisGT(scip, candsol, rootsolval) )
   {
      /* round up */
      switch( divetype )
      {
         case SCIP_DIVETYPE_INTEGRALITY:
            distquot = (1.0 - candsfrac + SCIPsumepsilon(scip)) / (candsol - rootsolval);
            break;
         case SCIP_DIVETYPE_SOS1VARIABLE:
            if ( SCIPisFeasPositive(scip, candsol) )
               distquot = (1.0 - candsfrac + SCIPsumepsilon(scip)) / (candsol - rootsolval);
            else
               distquot = (candsfrac + SCIPsumepsilon(scip)) / (rootsolval - candsol);
            break;
         default:
            SCIPerrorMessage("Error: Unsupported diving type\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
      } /*lint !e788*/

      /* avoid roundable candidates */
      if( SCIPvarMayRoundUp(cand) )
         distquot *= 1000.0;
      *roundup = TRUE;
   }
   else
   {
      /* if the solution values are equal, we arbitrarily select branching downwards;
       * candidates with equal LP solution values are penalized with an infinite score
       */
      *roundup = FALSE;
      distquot = SCIPinfinity(scip);
   }

   *score = -distquot;

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the linesearchdiving primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLinesearchdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Linesearchdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecLinesearchdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyLinesearchdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeLinesearchdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitLinesearchdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitLinesearchdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL,
         DEFAULT_LPRESOLVEDOMCHGQUOT, DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED,
         DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS, DIVESET_DIVETYPES, divesetGetScoreLinesearchdiving) );

   return SCIP_OKAY;
}
