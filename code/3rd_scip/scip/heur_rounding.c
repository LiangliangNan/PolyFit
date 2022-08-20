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

/**@file   heur_rounding.c
 * @brief  LP rounding heuristic that tries to recover from intermediate infeasibilities
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_rounding.h"


#define HEUR_NAME             "rounding"
#define HEUR_DESC             "LP rounding heuristic with infeasibility recovering"
#define HEUR_DISPCHAR         'R'
#define HEUR_PRIORITY         -1000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_SUCCESSFACTOR 100            /**< number of calls per found solution that are considered as standard success, 
                                              * a higher factor causes the heuristic to be called more often 
                                              */
#define DEFAULT_ONCEPERNODE   FALSE          /**< should the heuristic only be called once per node? */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Longint          lastlp;             /**< last LP number where the heuristic was applied */
   int                   successfactor;      /**< number of calls per found solution that are considered as standard success, 
                                              * a higher factor causes the heuristic to be called more often 
                                              */
   SCIP_Bool             oncepernode;        /**< should the heuristic only be called once per node? */
};


/*
 * local methods
 */

/** update row violation arrays after a row's activity value changed */
static
void updateViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_ROW**            violrows,           /**< array with currently violated rows */
   int*                  violrowpos,         /**< position of LP rows in violrows array */
   int*                  nviolrows,          /**< pointer to the number of currently violated rows */
   SCIP_Real             oldactivity,        /**< old activity value of LP row */
   SCIP_Real             newactivity         /**< new activity value of LP row */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool oldviol;
   SCIP_Bool newviol;

   assert(violrows != NULL);
   assert(violrowpos != NULL);
   assert(nviolrows != NULL);

   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);
   oldviol = (SCIPisFeasLT(scip, oldactivity, lhs) || SCIPisFeasGT(scip, oldactivity, rhs));
   newviol = (SCIPisFeasLT(scip, newactivity, lhs) || SCIPisFeasGT(scip, newactivity, rhs));
   if( oldviol != newviol )
   {
      int rowpos;
      int violpos;

      rowpos = SCIProwGetLPPos(row);
      assert(rowpos >= 0);

      if( oldviol )
      {
         /* the row violation was repaired: remove row from violrows array, decrease violation count */
         violpos = violrowpos[rowpos];
         assert(0 <= violpos && violpos < *nviolrows);
         assert(violrows[violpos] == row);
         violrowpos[rowpos] = -1;
         if( violpos != *nviolrows-1 )
         {
            violrows[violpos] = violrows[*nviolrows-1];
            violrowpos[SCIProwGetLPPos(violrows[violpos])] = violpos;
         }
         (*nviolrows)--;
      }
      else
      {
         /* the row is now violated: add row to violrows array, increase violation count */
         assert(violrowpos[rowpos] == -1);
         violrows[*nviolrows] = row;
         violrowpos[rowpos] = *nviolrows;
         (*nviolrows)++;
      }
   }
}

/** update row activities after a variable's solution value changed */
static
SCIP_RETCODE updateActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            activities,         /**< LP row activities */
   SCIP_ROW**            violrows,           /**< array with currently violated rows */
   int*                  violrowpos,         /**< position of LP rows in violrows array */
   int*                  nviolrows,          /**< pointer to the number of currently violated rows */
   int                   nlprows,            /**< number of rows in current LP */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldsolval,          /**< old solution value of variable */
   SCIP_Real             newsolval           /**< new solution value of variable */
   )
{
   SCIP_COL* col;
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   SCIP_Real delta;
   int ncolrows;
   int r;

   assert(activities != NULL);
   assert(nviolrows != NULL);
   assert(0 <= *nviolrows && *nviolrows <= nlprows);

   delta = newsolval - oldsolval;
   col = SCIPvarGetCol(var);
   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   ncolrows = SCIPcolGetNLPNonz(col);
   assert(ncolrows == 0 || (colrows != NULL && colvals != NULL));

   for( r = 0; r < ncolrows; ++r )
   {
      SCIP_ROW* row;
      int rowpos;

      row = colrows[r];
      rowpos = SCIProwGetLPPos(row);
      assert(-1 <= rowpos && rowpos < nlprows);

      if( rowpos >= 0 && !SCIProwIsLocal(row) )
      {
         SCIP_Real oldactivity;
         SCIP_Real newactivity;

         assert(SCIProwIsInLP(row));

         /* update row activity */
         oldactivity = activities[rowpos];
         if( !SCIPisInfinity(scip, -oldactivity) && !SCIPisInfinity(scip, oldactivity) )
         {
            newactivity = oldactivity + delta * colvals[r];
            if( SCIPisInfinity(scip, newactivity) )
               newactivity = SCIPinfinity(scip);
            else if( SCIPisInfinity(scip, -newactivity) )
               newactivity = -SCIPinfinity(scip);
            activities[rowpos] = newactivity;

            /* update row violation arrays */
            updateViolations(scip, row, violrows, violrowpos, nviolrows, oldactivity, newactivity);
         }
      }
   }

   return SCIP_OKAY;
}

/** returns a variable, that pushes activity of the row in the given direction with minimal negative impact on other rows;
 *  if variables have equal impact, chooses the one with best objective value improvement in corresponding direction;
 *  rounding in a direction is forbidden, if this forces the objective value over the upper bound
 */
static
SCIP_RETCODE selectRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   SCIP_ROW*             row,                /**< LP row */
   int                   direction,          /**< should the activity be increased (+1) or decreased (-1)? */
   SCIP_VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   SCIP_Real*            oldsolval,          /**< pointer to store old (fractional) solution value of rounding variable */
   SCIP_Real*            newsolval           /**< pointer to store new (rounded) solution value of rounding variable */
   )
{
   SCIP_COL* col;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   SCIP_Real solval;
   SCIP_Real roundval;
   SCIP_Real obj;
   SCIP_Real deltaobj;
   SCIP_Real bestdeltaobj;
   SCIP_VARTYPE vartype;
   int nrowcols;
   int nlocks;
   int minnlocks;
   int c;

   assert(direction == +1 || direction == -1);
   assert(roundvar != NULL);
   assert(oldsolval != NULL);
   assert(newsolval != NULL);

   /* get row entries */
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   nrowcols = SCIProwGetNLPNonz(row);

   /* select rounding variable */
   minnlocks = INT_MAX;
   bestdeltaobj = SCIPinfinity(scip);
   *roundvar = NULL;
   for( c = 0; c < nrowcols; ++c )
   {
      col = rowcols[c];
      var = SCIPcolGetVar(col);

      vartype = SCIPvarGetType(var);
      if( vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER )
      {
         solval = SCIPgetSolVal(scip, sol, var);

         if( !SCIPisFeasIntegral(scip, solval) )
         {
            val = rowvals[c];
            obj = SCIPvarGetObj(var);

            if( direction * val < 0.0 )
            {
               /* rounding down */
               nlocks = SCIPvarGetNLocksDown(var);
               if( nlocks <= minnlocks )
               {
                  roundval = SCIPfeasFloor(scip, solval);
                  deltaobj = obj * (roundval - solval);
                  if( (nlocks < minnlocks || deltaobj < bestdeltaobj) && minobj - obj < SCIPgetCutoffbound(scip) )
                  {
                     minnlocks = nlocks;
                     bestdeltaobj = deltaobj;
                     *roundvar = var;
                     *oldsolval = solval;
                     *newsolval = roundval;
                  }
               }
            }
            else
            {
               /* rounding up */
               assert(direction * val > 0.0);
               nlocks = SCIPvarGetNLocksUp(var);
               if( nlocks <= minnlocks )
               {
                  roundval = SCIPfeasCeil(scip, solval);
                  deltaobj = obj * (roundval - solval);
                  if( (nlocks < minnlocks || deltaobj < bestdeltaobj) && minobj + obj < SCIPgetCutoffbound(scip) )
                  {
                     minnlocks = nlocks;
                     bestdeltaobj = deltaobj;
                     *roundvar = var;
                     *oldsolval = solval;
                     *newsolval = roundval;
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** returns a variable, that increases the activity of the row */
static
SCIP_RETCODE selectIncreaseRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   SCIP_Real*            oldsolval,          /**< old (fractional) solution value of rounding variable */
   SCIP_Real*            newsolval           /**< new (rounded) solution value of rounding variable */
   )
{
   return selectRounding(scip, sol, minobj, row, +1, roundvar, oldsolval, newsolval);
}

/** returns a variable, that decreases the activity of the row */
static
SCIP_RETCODE selectDecreaseRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   SCIP_Real*            oldsolval,          /**< old (fractional) solution value of rounding variable */
   SCIP_Real*            newsolval           /**< new (rounded) solution value of rounding variable */
   )
{
   return selectRounding(scip, sol, minobj, row, -1, roundvar, oldsolval, newsolval);
}

/** returns a fractional variable, that has most impact on rows in opposite direction, i.e. that is most crucial to
 *  fix in the other direction;
 *  if variables have equal impact, chooses the one with best objective value improvement in corresponding direction;
 *  rounding in a direction is forbidden, if this forces the objective value over the upper bound
 */
static
SCIP_RETCODE selectEssentialRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Real             minobj,             /**< minimal objective value possible after rounding remaining fractional vars */
   SCIP_VAR**            lpcands,            /**< fractional variables in LP */
   int                   nlpcands,           /**< number of fractional variables in LP */
   SCIP_VAR**            roundvar,           /**< pointer to store the rounding variable, returns NULL if impossible */
   SCIP_Real*            oldsolval,          /**< old (fractional) solution value of rounding variable */
   SCIP_Real*            newsolval           /**< new (rounded) solution value of rounding variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real solval;
   SCIP_Real roundval;
   SCIP_Real obj;
   SCIP_Real deltaobj;
   SCIP_Real bestdeltaobj;
   int maxnlocks;
   int nlocks;
   int v;

   assert(roundvar != NULL);
   assert(oldsolval != NULL);
   assert(newsolval != NULL);

   /* select rounding variable */
   maxnlocks = -1;
   bestdeltaobj = SCIPinfinity(scip);
   *roundvar = NULL;
   for( v = 0; v < nlpcands; ++v )
   {
      var = lpcands[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      solval = SCIPgetSolVal(scip, sol, var);
      if( !SCIPisFeasIntegral(scip, solval) )
      {
         obj = SCIPvarGetObj(var);

         /* rounding down */
         nlocks = SCIPvarGetNLocksUp(var);
         if( nlocks >= maxnlocks )
         {
            roundval = SCIPfeasFloor(scip, solval);
            deltaobj = obj * (roundval - solval);
            if( (nlocks > maxnlocks || deltaobj < bestdeltaobj) && minobj - obj < SCIPgetCutoffbound(scip) )
            {
               maxnlocks = nlocks;
               bestdeltaobj = deltaobj;
               *roundvar = var;
               *oldsolval = solval;
               *newsolval = roundval;
            }
         }

         /* rounding up */
         nlocks = SCIPvarGetNLocksDown(var);
         if( nlocks >= maxnlocks )
         {
            roundval = SCIPfeasCeil(scip, solval);
            deltaobj = obj * (roundval - solval);
            if( (nlocks > maxnlocks || deltaobj < bestdeltaobj) && minobj + obj < SCIPgetCutoffbound(scip) )
            {
               maxnlocks = nlocks;
               bestdeltaobj = deltaobj;
               *roundvar = var;
               *oldsolval = solval;
               *newsolval = roundval;
            }
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRounding)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRounding(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRounding) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitRounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create heuristic data */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   heurdata->lastlp = -1;

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolRounding)
{
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->lastlp = -1;

   /* change the heuristic's timingmask, if nit should be called only once per node */
   if( heurdata->oncepernode )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_AFTERLPNODE);

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolRounding)
{
   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_ROW** lprows;
   SCIP_Real* activities;
   SCIP_ROW** violrows;
   int* violrowpos;
   SCIP_Real obj;
   SCIP_Real bestroundval;
   SCIP_Real minobj;
   int nlpcands;
   int nlprows;
   int nfrac;
   int nviolrows;
   int c;
   int r;
   SCIP_Longint nlps;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   SCIP_Longint nnodes;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't call heuristic, if we have already processed the current LP solution */
   nlps = SCIPgetNLPs(scip);
   if( nlps == heurdata->lastlp )
      return SCIP_OKAY;
   heurdata->lastlp = nlps;

   /* don't call heuristic, if it was not successful enough in the past */
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + SCIPheurGetNSolsFound(heur);
   nnodes = SCIPgetNNodes(scip);
   if( nnodes % ((ncalls/heurdata->successfactor)/(nsolsfound+1)+1) != 0 )
      return SCIP_OKAY;

   /* get fractional variables, that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );
   nfrac = nlpcands;

   /* only call heuristic, if LP solution is fractional */
   if( nfrac == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get LP rows */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );

   SCIPdebugMsg(scip, "executing rounding heuristic: %d LP rows, %d fractionals\n", nlprows, nfrac);

   /* get memory for activities, violated rows, and row violation positions */
   SCIP_CALL( SCIPallocBufferArray(scip, &activities, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violrows, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violrowpos, nlprows) );

   /* get the activities for all globally valid rows;
    * the rows should be feasible, but due to numerical inaccuracies in the LP solver, they can be violated
    */
   nviolrows = 0;
   for( r = 0; r < nlprows; ++r )
   {
      SCIP_ROW* row;

      row = lprows[r];
      assert(SCIProwGetLPPos(row) == r);

      if( !SCIProwIsLocal(row) )
      {
         activities[r] = SCIPgetRowActivity(scip, row);
         if( SCIPisFeasLT(scip, activities[r], SCIProwGetLhs(row))
            || SCIPisFeasGT(scip, activities[r], SCIProwGetRhs(row)) )
         {
            violrows[nviolrows] = row;
            violrowpos[r] = nviolrows;
            nviolrows++;
         }
         else
            violrowpos[r] = -1;
      }
   }

   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert(sol != NULL);

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   /* calculate the minimal objective value possible after rounding fractional variables */
   minobj = SCIPgetSolTransObj(scip, sol);
   assert(minobj < SCIPgetCutoffbound(scip));
   for( c = 0; c < nlpcands; ++c )
   {
      obj = SCIPvarGetObj(lpcands[c]);
      bestroundval = obj > 0.0 ? SCIPfeasFloor(scip, lpcandssol[c]) : SCIPfeasCeil(scip, lpcandssol[c]);
      minobj += obj * (bestroundval - lpcandssol[c]);
   }

   /* try to round remaining variables in order to become/stay feasible */
   while( nfrac > 0 )
   {
      SCIP_VAR* roundvar;
      SCIP_Real oldsolval;
      SCIP_Real newsolval;

      SCIPdebugMsg(scip, "rounding heuristic: nfrac=%d, nviolrows=%d, obj=%g (best possible obj: %g)\n",
         nfrac, nviolrows, SCIPgetSolOrigObj(scip, sol), SCIPretransformObj(scip, minobj));

      /* minobj < SCIPgetCutoffbound(scip) should be true, otherwise the rounding variable selection
       * should have returned NULL. Due to possible cancellation we use SCIPisLE. */
      assert( SCIPisLE(scip, minobj, SCIPgetCutoffbound(scip)) );

      /* choose next variable to process:
       *  - if a violated row exists, round a variable decreasing the violation, that has least impact on other rows
       *  - otherwise, round a variable, that has strongest devastating impact on rows in opposite direction
       */
      if( nviolrows > 0 )
      {
         SCIP_ROW* row;
         int rowpos;

         row = violrows[nviolrows-1];
         rowpos = SCIProwGetLPPos(row);
         assert(0 <= rowpos && rowpos < nlprows);
         assert(violrowpos[rowpos] == nviolrows-1);

         SCIPdebugMsg(scip, "rounding heuristic: try to fix violated row <%s>: %g <= %g <= %g\n",
            SCIProwGetName(row), SCIProwGetLhs(row), activities[rowpos], SCIProwGetRhs(row));
         if( SCIPisFeasLT(scip, activities[rowpos], SCIProwGetLhs(row)) )
         {
            /* lhs is violated: select a variable rounding, that increases the activity */
            SCIP_CALL( selectIncreaseRounding(scip, sol, minobj, row, &roundvar, &oldsolval, &newsolval) );
         }
         else
         {
            assert(SCIPisFeasGT(scip, activities[rowpos], SCIProwGetRhs(row)));
            /* rhs is violated: select a variable rounding, that decreases the activity */
            SCIP_CALL( selectDecreaseRounding(scip, sol, minobj, row, &roundvar, &oldsolval, &newsolval) );
         }
      }
      else
      {
         SCIPdebugMsg(scip, "rounding heuristic: search rounding variable and try to stay feasible\n");
         SCIP_CALL( selectEssentialRounding(scip, sol, minobj, lpcands, nlpcands, &roundvar, &oldsolval, &newsolval) );
      }

      /* check, whether rounding was possible */
      if( roundvar == NULL )
      {
         SCIPdebugMsg(scip, "rounding heuristic:  -> didn't find a rounding variable\n");
         break;
      }

      SCIPdebugMsg(scip, "rounding heuristic:  -> round var <%s>, oldval=%g, newval=%g, obj=%g\n",
         SCIPvarGetName(roundvar), oldsolval, newsolval, SCIPvarGetObj(roundvar));

      /* update row activities of globally valid rows */
      SCIP_CALL( updateActivities(scip, activities, violrows, violrowpos, &nviolrows, nlprows, 
            roundvar, oldsolval, newsolval) );

      /* store new solution value and decrease fractionality counter */
      SCIP_CALL( SCIPsetSolVal(scip, sol, roundvar, newsolval) );
      nfrac--;

      /* update minimal objective value possible after rounding remaining variables */
      obj = SCIPvarGetObj(roundvar);
      if( obj > 0.0 && newsolval > oldsolval )
         minobj += obj;
      else if( obj < 0.0 && newsolval < oldsolval )
         minobj -= obj;

      SCIPdebugMsg(scip, "rounding heuristic:  -> nfrac=%d, nviolrows=%d, obj=%g (best possible obj: %g)\n",
         nfrac, nviolrows, SCIPgetSolOrigObj(scip, sol), SCIPretransformObj(scip, minobj));
   }

   /* check, if the new solution is feasible */
   if( nfrac == 0 && nviolrows == 0 )
   {
      SCIP_Bool stored;

      /* check solution for feasibility, and add it to solution store if possible
       * neither integrality nor feasibility of LP rows has to be checked, because this is already
       * done in the rounding heuristic itself; however, be better check feasibility of LP rows,
       * because of numerical problems with activity updating
       */
      SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, TRUE, &stored) );

      if( stored )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "found feasible rounded solution:\n");
         SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
         *result = SCIP_FOUNDSOL;
      }
   }

   /* free memory buffers */
   SCIPfreeBufferArray(scip, &violrowpos);
   SCIPfreeBufferArray(scip, &violrows);
   SCIPfreeBufferArray(scip, &activities);

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the rounding heuristic with infeasibility recovering and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Rounding primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRounding, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRounding) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRounding) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRounding) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRounding) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRounding) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRounding) );

   /* add rounding primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/successfactor",
         "number of calls per found solution that are considered as standard success, a higher factor causes the heuristic to be called more often",
         &heurdata->successfactor, TRUE, DEFAULT_SUCCESSFACTOR, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/oncepernode",
         "should the heuristic only be called once per node?",
         &heurdata->oncepernode, TRUE, DEFAULT_ONCEPERNODE, NULL, NULL) );

   return SCIP_OKAY;
}
