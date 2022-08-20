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

/**@file   heur_locks.c
 * @brief  rounding locks primal heuristic
 * @author Michael Winkler
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_locks.h"

#define HEUR_NAME             "locks"
#define HEUR_DESC             "heuristic that fixes variables based on their rounding locks"
#define HEUR_DISPCHAR         'k'
#define HEUR_PRIORITY         3000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE                       /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXNODES      5000LL                     /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_ROUNDUPPROBABILITY 0.67                  /**< probability for rounding a variable up in case of ties */
#define DEFAULT_MINFIXINGRATE 0.65                       /**< minimum percentage of variables that have to be fixed */
#define DEFAULT_MINIMPROVE    0.01                       /**< factor by which locks heuristic should at least improve the
                                                          *   incumbent
                                                          */
#define DEFAULT_MINNODES      500LL                      /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_NODESOFS      500LL                      /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT     0.1                        /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_MAXPROPROUNDS 2                          /**< maximum number of propagation rounds during probing */
#define DEFAULT_UPDATELOCKS   TRUE                       /**< should the locks be updated based on LP rows? */
#define DEFAULT_COPYCUTS      TRUE                       /**< should all active cuts from the cutpool of the
                                                          *   original scip be copied to constraints of the subscip? */
#define DEFAULT_USEFINALSUBMIP TRUE                      /**< should a final sub-MIP be solved to construct a feasible
                                                          *   solution if the LP was not roundable? */
#define DEFAULT_RANDSEED      73                         /**< initial random seed */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generation */
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Longint          usednodes;          /**< nodes already used by locks heuristic in earlier calls */
   SCIP_Real             roundupprobability; /**< probability for rounding a variable up in case of ties */
   SCIP_Real             minfixingrate;      /**< minimum percentage of variables that have to be fixed */
   SCIP_Real             minimprove;         /**< factor by which locks heuristic should at least improve the incumbent */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   int                   maxproprounds;      /**< maximum number of propagation rounds during probing */
   SCIP_Bool             updatelocks;        /**< should the locks be updated based on LP rows? */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   the subproblem? */
   SCIP_Bool             usefinalsubmip;     /**< should a final sub-MIP be solved to costruct a feasible solution if
                                              *   the LP was not roundable? */
};

/*
 * Local methods
 */

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             newsol,             /**< working solution */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);
   assert(success != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyLocks)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurLocks(scip) );

   return SCIP_OKAY;
}

/** free method for primal heuristic plugins (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeLocks)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);

    /* free primal heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitLocks) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* initialize data */
   heurdata->usednodes = 0;

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitLocks) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   return SCIP_OKAY;
}

#define heurInitsolLocks NULL
#define heurExitsolLocks NULL

/** apply fix-and-propagate scheme based on variable locks
 *
 *  @note probing mode of SCIP needs to be enabled before
 */
SCIP_RETCODE SCIPapplyLockFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< primal heuristic data */
   SCIP_Bool*            cutoff,             /**< pointer to store if a cutoff was detected */
   SCIP_Bool*            allrowsfulfilled    /**< pointer to store if all rows became redundant */
   )
{
   SCIP_ROW** lprows;
   SCIP_VAR** vars;
   SCIP_VAR** sortvars;
   SCIP_Real* minact;
   SCIP_Real* maxact;
   SCIP_Bool* fulfilled;
   SCIP_VAR* var;
   SCIP_ROW* row;
   SCIP_COL* col;
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   int ncolrows;
   int* ndownlocks;
   int* nuplocks;
   int* varpos = NULL;
   SCIP_Real lastfixval;
   SCIP_Real randnumber;
   SCIP_Real roundupprobability;
   SCIP_Bool propagate;
   SCIP_Bool propagated;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   SCIP_Bool updatelocks;
   int lastfixlocks;
   int maxproprounds;
   int nglbfulfilledrows;
   int rowpos;
   int nbinvars;
   int nvars;
   int nlprows;
   int nfulfilledrows;
   int bestpos;
   int lastbestscore;
   int bestscore;
   int score;
   int v;
   int r;
   int i;

   assert(scip != NULL);
   assert(cutoff != NULL);
   assert(allrowsfulfilled != NULL);
   assert(SCIPinProbing(scip));

   if( heurdata == NULL )
   {
      SCIP_HEUR* heur = SCIPfindHeur(scip, HEUR_NAME);
      heurdata = SCIPheurGetData(heur);
   }
   assert(heurdata != NULL);

   *cutoff = FALSE;
   *allrowsfulfilled = FALSE;

   propagate = (heurdata->maxproprounds != 0);

   if( heurdata->maxproprounds == -2 )
      maxproprounds = 0;
   else
      maxproprounds = heurdata->maxproprounds;

   roundupprobability = heurdata->roundupprobability;


   updatelocks = heurdata->updatelocks && (SCIPgetNCheckConss(scip) == SCIPgetNLPRows(scip));

   SCIPdebugMsg(scip, "%d constraints: %d logicor, updatelocks=%d\n", SCIPgetNConss(scip), SCIPconshdlrGetNCheckConss(SCIPfindConshdlr(scip, "logicor")), updatelocks);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   assert(vars != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortvars, vars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nuplocks, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ndownlocks, nbinvars) );

   /* get LP data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minact, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxact, nlprows) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &fulfilled, nlprows) );

   /* @todo add objective value as second sorting criteria */

   nglbfulfilledrows = 0;

   /* get locks of variables */
   for( v = 0; v < nbinvars; ++v )
   {
      var = sortvars[v];
      nuplocks[v] = SCIPvarGetNLocksUp(var);
      ndownlocks[v] = SCIPvarGetNLocksDown(var);
   }

   /* get activities of rows */
   for( r = 0; r < nlprows; ++r )
   {
      row = lprows[r];
      assert(SCIProwGetLPPos(row) == r);

      /* no trivial rows */
      assert(!SCIPisInfinity(scip, -SCIProwGetLhs(row)) || !SCIPisInfinity(scip, SCIProwGetRhs(row)));

      minact[r] = SCIPgetRowMinActivity(scip, row);
      maxact[r] = SCIPgetRowMaxActivity(scip, row);
   }

   propagated = TRUE;
   lastbestscore = INT_MAX;

   /* fix variables */
   for( v = 0; v < nbinvars; v++ )
   {
      if( SCIPisStopped(scip) )
         break;

      assert(!(*cutoff));

      nfulfilledrows = 0;

      while( v < nbinvars && (SCIPvarGetLbLocal(sortvars[v]) > 0.5 || SCIPvarGetUbLocal(sortvars[v]) < 0.5) )
      {
         ++v;
      }
      if( v == nbinvars )
         break;

      bestpos = v;
      bestscore = nuplocks[v] + ndownlocks[v];

      /* get best variable */
      if( bestscore < lastbestscore )
      {
         for( i = v + 1; i < nbinvars; ++i )
         {
            var = sortvars[i];

            /* variable is already fixed; move it to the front and increment v to ignore it */
            if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5 )
            {
               int locks;

               sortvars[i] = sortvars[v];
               sortvars[v] = var;

               locks = nuplocks[i];
               nuplocks[i] = nuplocks[v];
               nuplocks[v] = locks;

               locks = ndownlocks[i];
               ndownlocks[i] = ndownlocks[v];
               ndownlocks[v] = locks;

               if( varpos != NULL )
               {
                  varpos[SCIPvarGetProbindex(sortvars[i])] = i;
                  varpos[SCIPvarGetProbindex(sortvars[v])] = v;
               }

               if( bestpos == v )
                  bestpos = i;

               ++v;

               continue;
            }

            score = nuplocks[i] + ndownlocks[i];
            assert(score <= lastbestscore);

            if( score > bestscore )
            {
               bestscore = score;
               bestpos = i;

               if( bestscore == lastbestscore )
                  break;
            }
         }
         if( v == nbinvars )
            break;
      }
      lastbestscore = bestscore;

      /* move best variable to the front (at position v) */
      if( bestpos != v )
      {
         int locks;

         var = sortvars[bestpos];
         sortvars[bestpos] = sortvars[v];
         sortvars[v] = var;

         locks = nuplocks[bestpos];
         nuplocks[bestpos] = nuplocks[v];
         nuplocks[v] = locks;

         locks = ndownlocks[bestpos];
         ndownlocks[bestpos] = ndownlocks[v];
         ndownlocks[v] = locks;

         if( varpos != NULL )
         {
            varpos[SCIPvarGetProbindex(sortvars[bestpos])] = bestpos;
            varpos[SCIPvarGetProbindex(sortvars[v])] = v;
         }
      }

      var = sortvars[v];

      /* all remaining variables are fixed, we can break the fix-and-propagate loop */
      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5 )
      {
         assert(v == nbinvars);

         break;
      }

      /* stop if we reached the depth limit */
      if( SCIP_MAXTREEDEPTH <= SCIPgetDepth(scip) )
         break;

      if( propagated )
      {
         SCIP_CALL( SCIPnewProbingNode(scip) );
         propagated = FALSE;
      }

      /* set variables to the bound with fewer locks, if tie choose an average value */
      if( ndownlocks[v] > nuplocks[v] )
         lastfixval = 1.0;
      else if( ndownlocks[v] < nuplocks[v] )
         lastfixval = 0.0;
      else
      {
         /* prefer one-fixing if objective value is not positive */
         if( !SCIPisPositive(scip, SCIPvarGetObj(var)) )
            lastfixval = 1.0;
         else
         {
            randnumber = SCIPrandomGetReal(heurdata->randnumgen, 0.0, 1.0);

            /* if a tie occurs, we randomly round the variable based on the parameter 'roundupprobability' */
            if( randnumber < roundupprobability )
               lastfixval = 1.0;
            else
               lastfixval = 0.0;
         }
      }

      lastfixlocks = lastfixval > 0.5 ? nuplocks[v] : ndownlocks[v];

      SCIP_CALL( SCIPfixVarProbing(scip, var, lastfixval) );

      SCIPdebugMsg(scip, "iteration %d: fixing variable <%s> to %d with locks (%d, %d)\n", v, SCIPvarGetName(var), lastfixval > 0.5 ? 1 : 0, ndownlocks[v], nuplocks[v]);

      if( propagate && lastfixlocks > 0 )
      {
         /* apply propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, maxproprounds, cutoff, NULL) );
         propagated = TRUE;

         if( *cutoff )
         {
            SCIPdebugMsg(scip, "last fixing led to infeasibility trying other bound\n");

            /* fix cutoff variable in other direction */
            SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );
            *cutoff = FALSE;

            if( lastfixval < 0.5 )
            {
               lastfixval = 1.0;

               if( SCIPvarGetUbLocal(var) > 0.5 )
               {
                  SCIP_CALL( SCIPfixVarProbing(scip, var, 1.0) );
               }
               /* because of the limited number of propagation rounds, it may happen that conflict analysis finds a
                * valid global fixing for the last fixed variable that conflicts with applying the reverse fixing
                * after backtracking; in that case, we ran into a deadend and stop
                */
               else
                  *cutoff = TRUE;
            }
            else
            {
               lastfixval = 0.0;

               if( SCIPvarGetLbLocal(var) < 0.5 )
               {
                  SCIP_CALL( SCIPfixVarProbing(scip, var, 0.0) );
               }
               /* because of the limited number of propagation rounds, it may happen that conflict analysis finds a
                * valid global fixing for the last fixed variable that conflicts with applying the reverse fixing
                * after backtracking; in that case, we ran into a deadend and stop
                */
               else
                  *cutoff = TRUE;
            }

            if( !(*cutoff) )
            {
               SCIP_CALL( SCIPpropagateProbing(scip, maxproprounds, cutoff, NULL) );
            }
            if( *cutoff )
            {
               SCIPdebugMsg(scip, "probing was infeasible\n");

               break;
            }
         }
         /* @todo collect propagated bounds and use them to update row activities as well */
      }

      if( updatelocks )
      {
         if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
            continue;

         col = SCIPvarGetCol(var);
         assert(col != NULL);

         colrows = SCIPcolGetRows(col);
         colvals = SCIPcolGetVals(col);
         ncolrows = SCIPcolGetNNonz(col);

         /* update activities */
         for( r = 0; r < ncolrows; ++r )
         {
            row = colrows[r];
            rowpos = SCIProwGetLPPos(row);

            /* the row is not in the LP */
            if( rowpos == -1 )
               continue;

            assert(lprows[rowpos] == row);

            /* we disregard cuts */
            if( SCIProwGetRank(row) > 0 )
               continue;

            /* the row is already fulfilled */
            if( fulfilled[rowpos] )
               continue;

            haslhs = !SCIPisInfinity(scip, -SCIProwGetLhs(row));
            hasrhs = !SCIPisInfinity(scip, SCIProwGetRhs(row));

            /* no trivial rows */
            assert(hasrhs || haslhs);

            if( ((colvals[r] > 0) == (lastfixval < 0.5)) )
            {
               maxact[rowpos] -= REALABS(colvals[r]);
            }
            if( ((colvals[r] > 0) == (lastfixval > 0.5)) )
            {
               minact[rowpos] += REALABS(colvals[r]);
            }

            /* check if the row cannot be violated anymore */
            if( (!haslhs || SCIPisFeasGE(scip, minact[rowpos], SCIProwGetLhs(row)))
               && (!hasrhs || SCIPisFeasLE(scip, maxact[rowpos], SCIProwGetRhs(row))) )
            {
               SCIP_COL** cols;
               SCIP_VAR* colvar;
               SCIP_Real* vals;
               int ncols;
               int pos;
               int w;

               SCIPdebugMsg(scip, "Row <%s> has activity [%g, %g], lhs=%g, rhs=%g\n", SCIProwGetName(row), minact[rowpos], maxact[rowpos], SCIProwGetLhs(row), SCIProwGetRhs(row));
               SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

               if( varpos == NULL )
               {
                  SCIP_CALL( SCIPallocBufferArray(scip, &varpos, nbinvars) );

                  for( pos = 0; pos < nbinvars; ++pos )
                     varpos[SCIPvarGetProbindex(sortvars[pos])] = pos;
               }

               ++nfulfilledrows;
               fulfilled[rowpos] = TRUE;
               cols = SCIProwGetCols(row);
               vals = SCIProwGetVals(row);
               ncols = SCIProwGetNNonz(row);

               /* we assume that all rows are locking the variables */
               for( w = ncols - 1; w >= 0; --w  )
               {
                  colvar = SCIPcolGetVar(cols[w]);
                  if( SCIPvarGetType(colvar) == SCIP_VARTYPE_BINARY && colvar != var )
                  {
                     assert(sortvars[varpos[SCIPvarGetProbindex(colvar)]] == colvar);
                     pos = varpos[SCIPvarGetProbindex(colvar)];

                     if( haslhs )
                     {
                        if( vals[w] > 0.0 )
                           --(ndownlocks[pos]);
                        else
                           --(nuplocks[pos]);
                     }
                     if( hasrhs )
                     {
                        if( vals[w] > 0.0 )
                           --(nuplocks[pos]);
                        else
                           --(ndownlocks[pos]);
                     }
                  }
               }

               continue;
            }
            else if( SCIPisFeasLT(scip, maxact[rowpos], SCIProwGetLhs(row)) || SCIPisFeasGT(scip, minact[rowpos], SCIProwGetRhs(row)) )
            {
               *cutoff = TRUE;
               break;
            }
         }

         if( *cutoff )
         {
            SCIPdebugMsg(scip, "found infeasible row, stopping heur\n");
            break;
         }

         nglbfulfilledrows += nfulfilledrows;
         SCIPdebugMsg(scip, "last fixing led to %d fulfilled rows, now %d of %d rows are fulfilled\n", nfulfilledrows, nglbfulfilledrows, nlprows);

         if( nglbfulfilledrows == nlprows )
         {
            *allrowsfulfilled = TRUE;
            break;
         }
      }
   } /*lint --e{850}*/

   SCIPfreeBufferArrayNull(scip, &varpos);
   SCIPfreeBufferArray(scip, &fulfilled);
   SCIPfreeBufferArray(scip, &maxact);
   SCIPfreeBufferArray(scip, &minact);
   SCIPfreeBufferArray(scip, &ndownlocks);
   SCIPfreeBufferArray(scip, &nuplocks);
   SCIPfreeBufferArray(scip, &sortvars);

   return SCIP_OKAY;
}




/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecLocks)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_LPSOLSTAT lpstatus = SCIP_LPSOLSTAT_ERROR;
   SCIP_Real lowerbound;
   SCIP_Bool cutoff;
   SCIP_Bool lperror;
   SCIP_Bool allrowsfulfilled = FALSE;
#ifdef NOCONFLICT
   SCIP_Bool enabledconflicts;
#endif
   int oldnpscands;
   int npscands;

   int nvars;
   int i;

   *result = SCIP_DIDNOTRUN;

   /* only run once */
   if( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   if( SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;

   /* only run if we are allowed to solve an LP at the current node in the tree */
   if( !SCIPhasCurrentNodeLP(scip) )
      return SCIP_OKAY;

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

      /* manually cut off the node if the LP construction detected infeasibility (heuristics cannot return such a result) */
      if( cutoff )
      {
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetCurrentNode(scip)) );
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPflushLP(scip) );

      /* we need an LP */
      if( SCIPgetNLPRows(scip) == 0 )
         return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

#ifdef NOCONFLICT
   /* disable conflict analysis */
   SCIP_CALL( SCIPgetBoolParam(scip, "conflict/enable", &enabledconflicts) );

   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", FALSE) );
   }
#endif

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   lowerbound = SCIPgetLowerbound(scip);
   oldnpscands = SCIPgetNPseudoBranchCands(scip);

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

#ifdef COLLECTSTATISTICS
   SCIPenableVarHistory(scip);
#endif

   cutoff = FALSE;
   lperror = FALSE;

   SCIP_CALL( SCIPapplyLockFixings(scip, heurdata, &cutoff, &allrowsfulfilled) );

   if( cutoff || SCIPisStopped(scip) )
      goto TERMINATE;

   /* check that we had enough fixings */
   npscands = SCIPgetNPseudoBranchCands(scip);

   SCIPdebugMsg(scip, "npscands=%d, oldnpscands=%d, allrowsfulfilled=%u heurdata->minfixingrate=%g\n",
      npscands, oldnpscands, allrowsfulfilled, heurdata->minfixingrate);

   if( !allrowsfulfilled && npscands > oldnpscands * (1 - heurdata->minfixingrate) )
   {
      SCIPdebugMsg(scip, "--> too few fixings\n");

      goto TERMINATE;
   }
   else
   {
      SCIPdebugMsg(scip, "starting solving locks-lp at time %g\n", SCIPgetSolvingTime(scip));

      /* solve LP;
       * errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      {
         SCIP_Bool retstat;
         retstat = SCIPsolveProbingLP(scip, -1, &lperror, &cutoff);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while solving LP in LOCKS heuristic; LP solve terminated with code <%d>\n",
               retstat);
         }
      }
#else
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
#endif
      SCIPdebugMsg(scip, "ending solving locks-lp at time %g\n", SCIPgetSolvingTime(scip));

      lpstatus = SCIPgetLPSolstat(scip);

      SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
      SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));

      /* check if this is a feasible solution */
      if( !lperror && lpstatus == SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIP_Bool success;

         lowerbound = SCIPgetLPObjval(scip);

         /* copy the current LP solution to the working solution */
         SCIP_CALL( SCIPlinkLPSol(scip, sol) );

         SCIP_CALL( SCIProundSol(scip, sol, &success) );

         if( success )
         {
            SCIP_Bool stored;

            /* check solution for feasibility, and add it to solution store if possible.
             * Neither integrality nor feasibility of LP rows have to be checked, because they
             * are guaranteed by the heuristic at this stage.
             */
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, FALSE, &stored) );

            if( stored )
            {
#ifdef SCIP_MORE_DEBUG
               SCIP_Bool feasible;
               SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, TRUE, TRUE, TRUE, &feasible) );
               assert(feasible);
#endif

               SCIPdebugMsg(scip, "found feasible solution:\n");
               SCIPdebug(SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE)) );
               *result = SCIP_FOUNDSOL;
            }

            /* we found a solution, so we are done */
            goto TERMINATE;
         }
      }
   }

   if( heurdata->usefinalsubmip && !cutoff && !lperror && lpstatus != SCIP_LPSOLSTAT_INFEASIBLE && lpstatus != SCIP_LPSOLSTAT_OBJLIMIT )
   {
      SCIP* subscip;
      SCIP_VAR** subvars;
      SCIP_HASHMAP* varmap;
      SCIP_Longint nstallnodes;
      SCIP_Bool valid;

      /* calculate the maximal number of branching nodes until heuristic is aborted */
      nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

      /* reward locks heuristic if it succeeded often */
      nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
      nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-MIP as 100 nodes */
      nstallnodes += heurdata->nodesofs;

      /* determine the node limit for the current process */
      nstallnodes -= heurdata->usednodes;
      nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

      /* check whether we have enough nodes left to call subproblem solving */
      if( nstallnodes < heurdata->minnodes )
      {
         SCIPdebugMsg(scip, "skipping " HEUR_NAME ": nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
         goto TERMINATE;
      }

      /* check whether there is enough time and memory left */
      SCIP_CALL( SCIPcheckCopyLimits(scip, &valid) );

      if( !valid )
         goto TERMINATE;

      /* get all variables */
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

      /* create subproblem */
      SCIP_CALL( SCIPcreate(&subscip) );

      /* allocate temporary memory for subscip variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), nvars) );

      SCIP_CALL( SCIPcopy(scip, subscip, varmap, NULL, "_locks", FALSE, FALSE, TRUE, &valid) );

      if( heurdata->copycuts )
      {
         /* copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, FALSE, NULL) );
      }

      for( i = 0; i < nvars; i++ )
         subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);

      /* free hash map */
      SCIPhashmapFree(&varmap);

      /* do not abort subproblem on CTRL-C */
      SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifdef SCIP_DEBUG
      /* for debugging, enable full output */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
      SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
      /* disable statistic timing inside sub SCIP and output to console */
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

      /* set limits for the subproblem */
      SCIP_CALL( SCIPcopyLimits(scip, subscip) );
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );

      /* forbid call of heuristics and separators solving sub-CIPs */
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

      /* disable cutting plane separation */
      SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

      /* disable expensive presolving */
      SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

      /* use inference branching */
      if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
      }

      /* speed up sub-SCIP by not checking dual LP feasibility */
      SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

      /* employ a limit on the number of enforcement rounds in the quadratic constraint handler; this fixes the issue that
       * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
       * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
       * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no deductions shall be
       * made for the original SCIP
       */
      if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
      }

      /* if there is already a solution, add an objective cutoff */
      if( SCIPgetNSols(scip) > 0 )
      {
         SCIP_Real upperbound;
         SCIP_Real minimprove;
         SCIP_Real cutoffbound;

         minimprove = heurdata->minimprove;
         assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

         upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

         if( !SCIPisInfinity(scip, -1.0 * lowerbound) )
         {
            cutoffbound = (1-minimprove) * SCIPgetUpperbound(scip) + minimprove * lowerbound;
         }
         else
         {
            if( SCIPgetUpperbound ( scip ) >= 0 )
               cutoffbound = (1 - minimprove) * SCIPgetUpperbound(scip);
            else
               cutoffbound = (1 + minimprove) * SCIPgetUpperbound(scip);
         }
         cutoffbound = MIN(upperbound, cutoffbound);
         SCIP_CALL( SCIPsetObjlimit(subscip, cutoffbound) );
         SCIPdebugMsg(scip, "setting objlimit for subscip to %g\n", cutoffbound);
      }

      SCIPdebugMsg(scip, "starting solving locks-submip at time %g\n", SCIPgetSolvingTime(scip));

      /* solve the subproblem */
      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      {
         SCIP_RETCODE retstat;
         retstat = SCIPpresolve(subscip);
         if( retstat != SCIP_OKAY )
         {
            SCIPwarningMessage(scip, "Error while presolving subMIP in locks heuristic; sub-SCIP terminated with code <%d>\n", retstat);

            goto FREESCIPANDTERMINATE;
         }
      }
#else
      SCIP_CALL_ABORT( SCIPpresolve(subscip) );
#endif

      SCIPdebugMsg(scip, "locks heuristic presolved subproblem at time %g : %d vars, %d cons; fixing value = %g\n", SCIPgetSolvingTime(scip), SCIPgetNVars(subscip), SCIPgetNConss(subscip), ((nvars - SCIPgetNVars(subscip)) / (SCIP_Real)nvars));

      /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
       * to ensure that not only the MIP but also the LP relaxation is easy enough
       */
      if( ((nvars - SCIPgetNVars(subscip)) / (SCIP_Real)nvars) >= heurdata->minfixingrate )
      {
         SCIP_SOL** subsols;
         SCIP_Bool success;
         int nsubsols;

         SCIPdebugMsg(scip, "solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);

#ifdef NDEBUG
         {
            SCIP_RETCODE retstat;
            retstat = SCIPsolve(subscip);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving subMIP in locks heuristic; sub-SCIP terminated with code <%d>\n",retstat);

               goto FREESCIPANDTERMINATE;
            }
         }
#else
         SCIP_CALL_ABORT( SCIPsolve(subscip) );
#endif
         SCIPdebugMsg(scip, "ending solving locks-submip at time %g, status = %d\n", SCIPgetSolvingTime(scip), SCIPgetStatus(subscip));

         /* check, whether a solution was found; due to numerics, it might happen that not all solutions are feasible ->
          * try all solutions until one was accepted
          */
         nsubsols = SCIPgetNSols(subscip);
         subsols = SCIPgetSols(subscip);
         success = FALSE;

         for( i = 0; i < nsubsols && !success; ++i )
         {
            SCIP_CALL( createNewSol(scip, subscip, subvars, sol, subsols[i], &success) );
         }
         if( success )
            *result = SCIP_FOUNDSOL;
      }

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      heurdata->usednodes += SCIPgetNNodes(subscip);
#ifdef NDEBUG
   FREESCIPANDTERMINATE:
#endif
      /* free subproblem */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }


 TERMINATE:
   /* exit probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

#ifdef NOCONFLICT
   /* reset the conflict analysis */
   if( !SCIPisParamFixed(scip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", enabledconflicts) );
   }
#endif

   /* free all allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the locks primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurLocks(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyLocks,
         heurFreeLocks, heurInitLocks, heurExitLocks,
         heurInitsolLocks, heurExitsolLocks, heurExecLocks,
         heurdata) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "maximum number of propagation rounds to be performed in each propagation call (-1: no limit, -2: parameter settings)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -2, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "minimum percentage of integer variables that have to be fixable",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/roundupprobability",
         "probability for rounding a variable up in case of ties",
         &heurdata->roundupprobability, FALSE, DEFAULT_ROUNDUPPROBABILITY, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usefinalsubmip",
         "should a final sub-MIP be solved to costruct a feasible solution if the LP was not roundable?",
         &heurdata->usefinalsubmip, TRUE, DEFAULT_USEFINALSUBMIP, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which " HEUR_NAME " heuristic should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/updatelocks",
         "should the locks be updated based on LP rows?",
         &heurdata->updatelocks, TRUE, DEFAULT_UPDATELOCKS, NULL, NULL) );

   return SCIP_OKAY;
}
