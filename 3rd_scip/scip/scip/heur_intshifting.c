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

/**@file   heur_intshifting.c
 * @brief  LP rounding heuristic that tries to recover from intermediate infeasibilities, shifts integer variables, and
 *         solves a final LP to calculate feasible values for continuous variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_intshifting.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "intshifting"
#define HEUR_DESC             "LP rounding heuristic with infeasibility recovering and final LP solving"
#define HEUR_DISPCHAR         'i'
#define HEUR_PRIORITY         -10000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define MAXSHIFTINGS          50        /**< maximal number of non improving shiftings */
#define WEIGHTFACTOR          1.1
#define DEFAULT_RANDSEED      17

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Longint          lastlp;             /**< last LP number where the heuristic was applied */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
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
   int*                  nviolfracrows,      /**< pointer to the number of violated rows with fractional candidates */
   int*                  nfracsinrow,        /**< array with number of fractional variables for every row */
   SCIP_Real             oldminactivity,     /**< old minimal activity value of LP row */
   SCIP_Real             oldmaxactivity,     /**< old maximal activity value of LP row */
   SCIP_Real             newminactivity,     /**< new minimal activity value of LP row */
   SCIP_Real             newmaxactivity      /**< new maximal activity value of LP row */
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
   oldviol = (SCIPisFeasLT(scip, oldmaxactivity, lhs) || SCIPisFeasGT(scip, oldminactivity, rhs));
   newviol = (SCIPisFeasLT(scip, newmaxactivity, lhs) || SCIPisFeasGT(scip, newminactivity, rhs));
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
         /* first, move the row to the end of the subset of violated rows with fractional variables */
         if( nfracsinrow[rowpos] > 0 )
         {
            violrows[violpos] = violrows[*nviolrows-1];
            assert(violpos < *nviolfracrows);

            /* replace with last violated row containing fractional variables */
            if( violpos != *nviolfracrows - 1 )
            {
               violrows[violpos] = violrows[*nviolfracrows - 1];
               violrowpos[SCIProwGetLPPos(violrows[violpos])] = violpos;
               violpos = *nviolfracrows - 1;
            }
            (*nviolfracrows)--;
         }

         assert(violpos >= *nviolfracrows);

         /* swap row at the end of the violated array to the position of this row and decrease the counter */
         if( violpos != *nviolrows - 1 )
         {
            violrows[violpos] = violrows[*nviolrows - 1];
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

         /* if the row contains fractional variables, swap with the last violated row that has no fractional variables
          * at position *nviolfracrows
          */
         if( nfracsinrow[rowpos] > 0 )
         {
            if( *nviolfracrows < *nviolrows - 1 )
            {
               assert(nfracsinrow[SCIProwGetLPPos(violrows[*nviolfracrows])] == 0);

               violrows[*nviolrows - 1] = violrows[*nviolfracrows];
               violrowpos[SCIProwGetLPPos(violrows[*nviolrows - 1])] = *nviolrows - 1;

               violrows[*nviolfracrows] = row;
               violrowpos[rowpos] = *nviolfracrows;
            }
            (*nviolfracrows)++;
         }
      }
   }
}

/** update row activities after a variable's solution value changed */
static
SCIP_RETCODE updateActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            minactivities,      /**< LP row minimal activities */
   SCIP_Real*            maxactivities,      /**< LP row maximal activities */
   SCIP_ROW**            violrows,           /**< array with currently violated rows */
   int*                  violrowpos,         /**< position of LP rows in violrows array */
   int*                  nviolrows,          /**< pointer to the number of currently violated rows */
   int*                  nviolfracrows,      /**< pointer to the number of violated rows with fractional candidates */
   int*                  nfracsinrow,        /**< array with number of fractional variables for every row */
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

   assert(minactivities != NULL);
   assert(maxactivities != NULL);
   assert(nviolrows != NULL);
   assert(0 <= *nviolrows && *nviolrows <= nlprows);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

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
         SCIP_Real oldminactivity;
         SCIP_Real oldmaxactivity;
         SCIP_Real newminactivity;
         SCIP_Real newmaxactivity;

         assert(SCIProwIsInLP(row));

         /* update row activities */
         oldminactivity = minactivities[rowpos];
         oldmaxactivity = maxactivities[rowpos];

         if( !SCIPisInfinity(scip, -oldminactivity) )
         {
            newminactivity = oldminactivity + delta * colvals[r];
            minactivities[rowpos] = newminactivity;
         }
         else
            newminactivity = -SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, oldmaxactivity) )
         {
            newmaxactivity = oldmaxactivity + delta * colvals[r];
            maxactivities[rowpos] = newmaxactivity;
         }
         else
            newmaxactivity = SCIPinfinity(scip);

         /* update row violation arrays */
         updateViolations(scip, row, violrows, violrowpos, nviolrows, nviolfracrows, nfracsinrow, oldminactivity, oldmaxactivity,
            newminactivity, newmaxactivity);
      }
   }

   return SCIP_OKAY;
}

/** returns an integer variable, that pushes activity of the row in the given direction with minimal negative impact on
 *  other rows;
 *  if variables have equal impact, chooses the one with best objective value improvement in corresponding direction;
 *  prefer fractional integers over other variables in order to become integral during the process;
 *  shifting in a direction is forbidden, if this forces the objective value over the upper bound, or if the variable
 *  was already shifted in the opposite direction
 */
static
SCIP_RETCODE selectShifting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             rowactivity,        /**< activity of LP row */
   int                   direction,          /**< should the activity be increased (+1) or decreased (-1)? */
   SCIP_Real*            nincreases,         /**< array with weighted number of increasings per variables */
   SCIP_Real*            ndecreases,         /**< array with weighted number of decreasings per variables */
   SCIP_Real             increaseweight,     /**< current weight of increase/decrease updates */
   SCIP_VAR**            shiftvar,           /**< pointer to store the shifting variable, returns NULL if impossible */
   SCIP_Real*            oldsolval,          /**< pointer to store old solution value of shifting variable */
   SCIP_Real*            newsolval           /**< pointer to store new (shifted) solution value of shifting variable */
   )
{
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int nrowcols;
   SCIP_Real activitydelta;
   SCIP_Real bestshiftscore;
   SCIP_Real bestdeltaobj;
   int c;

   assert(direction == +1 || direction == -1);
   assert(nincreases != NULL);
   assert(ndecreases != NULL);
   assert(shiftvar != NULL);
   assert(oldsolval != NULL);
   assert(newsolval != NULL);

   /* get row entries */
   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   nrowcols = SCIProwGetNLPNonz(row);

   /* calculate how much the activity must be shifted in order to become feasible */
   activitydelta = (direction == +1 ? SCIProwGetLhs(row) - rowactivity : SCIProwGetRhs(row) - rowactivity);
   assert((direction == +1 && SCIPisPositive(scip, activitydelta))
      || (direction == -1 && SCIPisNegative(scip, activitydelta)));

   /* select shifting variable */
   bestshiftscore = SCIP_REAL_MAX;
   bestdeltaobj = SCIPinfinity(scip);
   *shiftvar = NULL;
   *newsolval = 0.0;
   *oldsolval = 0.0;
   for( c = 0; c < nrowcols; ++c )
   {
      SCIP_COL* col;
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real solval;
      SCIP_Real shiftval;
      SCIP_Real shiftscore;
      SCIP_Bool isfrac;
      SCIP_Bool increase;
      int probindex;

      col = rowcols[c];
      var = SCIPcolGetVar(col);
      val = rowvals[c];
      assert(!SCIPisZero(scip, val));
      solval = SCIPgetSolVal(scip, sol, var);

      /* only accept integer variables */
      if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY && SCIPvarGetType(var) != SCIP_VARTYPE_INTEGER )
         continue;

      isfrac = !SCIPisFeasIntegral(scip, solval);
      increase = (direction * val > 0.0);
      probindex = SCIPvarGetProbindex(var);

      /* calculate the score of the shifting (prefer smaller values) */
      if( isfrac )
         shiftscore = increase ? -1.0 / (SCIPvarGetNLocksUp(var) + 1.0) :
            -1.0 / (SCIPvarGetNLocksDown(var) + 1.0);
      else
      {
         if( increase )
            shiftscore = ndecreases[probindex]/increaseweight;
         else
            shiftscore = nincreases[probindex]/increaseweight;
         shiftscore += 1.0;
      }

      if( shiftscore <= bestshiftscore )
      {
         SCIP_Real deltaobj;

         if( !increase )
         {
            /* shifting down */
            assert(direction * val < 0.0);
            if( isfrac )
               shiftval = SCIPfeasFloor(scip, solval);
            else
            {
               SCIP_Real lb;

               assert(activitydelta/val < 0.0);
               shiftval = solval + activitydelta/val;
               assert(shiftval <= solval); /* may be equal due to numerical digit erasement in the subtraction */
               shiftval = SCIPfeasFloor(scip, shiftval);
               lb = SCIPvarGetLbGlobal(var);
               shiftval = MAX(shiftval, lb);
            }
         }
         else
         {
            /* shifting up */
            assert(direction * val > 0.0);
            if( isfrac )
               shiftval = SCIPfeasCeil(scip, solval);
            else
            {
               SCIP_Real ub;

               assert(activitydelta/val > 0.0);
               shiftval = solval + activitydelta/val;
               assert(shiftval >= solval); /* may be equal due to numerical digit erasement in the subtraction */
               shiftval = SCIPfeasCeil(scip, shiftval);
               ub = SCIPvarGetUbGlobal(var);
               shiftval = MIN(shiftval, ub);
            }
         }

         if( SCIPisEQ(scip, shiftval, solval) )
            continue;

         deltaobj = SCIPvarGetObj(var) * (shiftval - solval);
         if( (shiftscore < bestshiftscore || deltaobj < bestdeltaobj)
            && !SCIPisHugeValue(scip, REALABS(shiftval)) ) /* ignore candidates for which shiftval is too large */
         {
            bestshiftscore = shiftscore;
            bestdeltaobj = deltaobj;
            *shiftvar = var;
            *oldsolval = solval;
            *newsolval = shiftval;
         }
      }
   }

   return SCIP_OKAY;
}

/** returns a fractional variable, that has most impact on rows in opposite direction, i.e. that is most crucial to
 *  fix in the other direction;
 *  if variables have equal impact, chooses the one with best objective value improvement in corresponding direction;
 *  shifting in a direction is forbidden, if this forces the objective value over the upper bound
 */
static
SCIP_RETCODE selectEssentialRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Real             minobj,             /**< minimal objective value possible after shifting remaining fractional vars */
   SCIP_VAR**            lpcands,            /**< fractional variables in LP */
   int                   nlpcands,           /**< number of fractional variables in LP */
   SCIP_VAR**            shiftvar,           /**< pointer to store the shifting variable, returns NULL if impossible */
   SCIP_Real*            oldsolval,          /**< old (fractional) solution value of shifting variable */
   SCIP_Real*            newsolval           /**< new (shifted) solution value of shifting variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real solval;
   SCIP_Real shiftval;
   SCIP_Real obj;
   SCIP_Real deltaobj;
   SCIP_Real bestdeltaobj;
   int maxnlocks;
   int nlocks;
   int v;

   assert(shiftvar != NULL);
   assert(oldsolval != NULL);
   assert(newsolval != NULL);

   /* select shifting variable */
   maxnlocks = -1;
   bestdeltaobj = SCIPinfinity(scip);
   *shiftvar = NULL;
   for( v = 0; v < nlpcands; ++v )
   {
      var = lpcands[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      solval = SCIPgetSolVal(scip, sol, var);
      if( !SCIPisFeasIntegral(scip, solval) )
      {
         obj = SCIPvarGetObj(var);

         /* shifting down */
         nlocks = SCIPvarGetNLocksUp(var);
         if( nlocks >= maxnlocks )
         {
            shiftval = SCIPfeasFloor(scip, solval);
            deltaobj = obj * (shiftval - solval);
            if( (nlocks > maxnlocks || deltaobj < bestdeltaobj) && minobj - obj < SCIPgetCutoffbound(scip) )
            {
               maxnlocks = nlocks;
               bestdeltaobj = deltaobj;
               *shiftvar = var;
               *oldsolval = solval;
               *newsolval = shiftval;
            }
         }

         /* shifting up */
         nlocks = SCIPvarGetNLocksDown(var);
         if( nlocks >= maxnlocks )
         {
            shiftval = SCIPfeasCeil(scip, solval);
            deltaobj = obj * (shiftval - solval);
            if( (nlocks > maxnlocks || deltaobj < bestdeltaobj) && minobj + obj < SCIPgetCutoffbound(scip) )
            {
               maxnlocks = nlocks;
               bestdeltaobj = deltaobj;
               *shiftvar = var;
               *oldsolval = solval;
               *newsolval = shiftval;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** adds a given value to the fractionality counters of the rows in which the given variable appears */
static
void addFracCounter(
   int*                  nfracsinrow,        /**< array to store number of fractional variables per row */
   SCIP_ROW**            violrows,           /**< array with currently violated rows */
   int*                  violrowpos,         /**< position of LP rows in violrows array */
   int*                  nviolfracrows,      /**< pointer to store the number of violated rows with fractional variables */
   int                   nviolrows,          /**< the number of currently violated rows (stays unchanged in this method) */
   int                   nlprows,            /**< number of rows in LP */
   SCIP_VAR*             var,                /**< variable for which the counting should be updated */
   int                   incval              /**< value that should be added to the corresponding array entries */
   )
{
   SCIP_COL* col;
   SCIP_ROW** rows;
   int nrows;
   int r;

   assert(incval != 0);
   assert(nviolrows >= *nviolfracrows);

   col = SCIPvarGetCol(var);
   rows = SCIPcolGetRows(col);
   nrows = SCIPcolGetNLPNonz(col);
   for( r = 0; r < nrows; ++r )
   {
      int rowlppos;
      int theviolrowpos;
      SCIP_ROW* row;

      row = rows[r];
      assert(NULL != row);
      rowlppos = SCIProwGetLPPos(row);
      assert(0 <= rowlppos && rowlppos < nlprows);
      assert(!SCIProwIsLocal(row) || violrowpos[rowlppos] == -1);

      if( SCIProwIsLocal(row) )
         continue;

      nfracsinrow[rowlppos] += incval;
      assert(nfracsinrow[rowlppos] >= 0);

      theviolrowpos = violrowpos[rowlppos];

      /* swap positions in violrows array if fractionality has changed to 0 */
      if( theviolrowpos >= 0 )
      {
         /* first case: the number of fractional variables has become zero: swap row in violrows array to the
          * second part, containing only violated rows without fractional variables
          */
         if( nfracsinrow[rowlppos] == 0 )
         {
            assert(theviolrowpos <= *nviolfracrows - 1);

            /* nothing to do if row is already at the end of the first part, otherwise, swap it to the last position
             * and decrease the counter */
            if( theviolrowpos < *nviolfracrows - 1 )
            {
               violrows[theviolrowpos] = violrows[*nviolfracrows - 1];
               violrows[*nviolfracrows - 1] = row;


               violrowpos[SCIProwGetLPPos(violrows[theviolrowpos])] = theviolrowpos;
               violrowpos[rowlppos] = *nviolfracrows - 1;
            }
            (*nviolfracrows)--;
         }
         /* second interesting case: the number of fractional variables was zero before, swap it with the first
          * violated row without fractional variables
          */
         else if( nfracsinrow[rowlppos] == incval )
         {
            assert(theviolrowpos >= *nviolfracrows);

            /* nothing to do if the row is exactly located at index *nviolfracrows */
            if( theviolrowpos > *nviolfracrows )
            {
               violrows[theviolrowpos] = violrows[*nviolfracrows];
               violrows[*nviolfracrows] = row;


               violrowpos[SCIProwGetLPPos(violrows[theviolrowpos])] = theviolrowpos;
               violrowpos[rowlppos] = *nviolfracrows;
            }
            (*nviolfracrows)++;
         }
      }
   }
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyIntshifting)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurIntshifting(scip) );

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitIntshifting) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(SCIPheurGetData(heur) == NULL);

   /* create heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   heurdata->lastlp = -1;
   SCIPheurSetData(heur, heurdata);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen,
         DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitIntshifting) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   /* free random number generator */
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolIntshifting)
{
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->lastlp = -1;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIntshifting) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_ROW** lprows;
   SCIP_Real* minactivities;
   SCIP_Real* maxactivities;
   SCIP_ROW** violrows;
   SCIP_Real* nincreases;
   SCIP_Real* ndecreases;
   int* violrowpos;
   int* nfracsinrow;
   SCIP_Real increaseweight;
   SCIP_Real obj;
   SCIP_Real bestshiftval;
   SCIP_Real minobj;
   int nlpcands;
   int nlprows;
   int nvars;
   int nfrac;
   int nviolrows;
   int nviolfracrows;
   int nprevviolrows;
   int minnviolrows;
   int nnonimprovingshifts;
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

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* don't call heuristic, if no continuous variables are present
    *  -> in this case, it is equivalent to shifting heuristic
    */
   if( SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

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
   if( nnodes % (ncalls/(nsolsfound+1)+1) != 0 )  /*?????????? ncalls/100 */
      return SCIP_OKAY;

   /* get fractional variables, that should be integral */
   /* todo check if heuristic should include implicit integer variables for its calculations */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, NULL) );
   nfrac = nlpcands;

   /* only call heuristic, if LP solution is fractional */
   if( nfrac == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get LP rows */
   SCIP_CALL( SCIPgetLPRowsData(scip, &lprows, &nlprows) );

   SCIPdebugMsg(scip, "executing intshifting heuristic: %d LP rows, %d fractionals\n", nlprows, nfrac);

   /* get memory for activities, violated rows, and row violation positions */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &minactivities, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxactivities, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violrows, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &violrowpos, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nfracsinrow, nlprows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nincreases, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ndecreases, nvars) );
   BMSclearMemoryArray(nfracsinrow, nlprows);
   BMSclearMemoryArray(nincreases, nvars);
   BMSclearMemoryArray(ndecreases, nvars);

   /* get the minimal and maximal activity for all globally valid rows for continuous variables in their full range;
    * these are the values of a*x' with x' being the LP solution for integer variables and the lower or upper bound
    * for the continuous variables
    */
   nviolrows = 0;
   for( r = 0; r < nlprows; ++r )
   {
      SCIP_ROW* row;

      row = lprows[r];
      assert(SCIProwGetLPPos(row) == r);

      if( !SCIProwIsLocal(row) )
      {
         SCIP_COL** cols;
         SCIP_Real* vals;
         int nnonz;
         SCIP_Bool mininf;
         SCIP_Bool maxinf;

         mininf = FALSE;
         maxinf = FALSE;
         minactivities[r] = 0.0;
         maxactivities[r] = 0.0;
         cols = SCIProwGetCols(row);
         vals = SCIProwGetVals(row);
         nnonz = SCIProwGetNNonz(row);
         for( c = 0; c < nnonz && !(mininf && maxinf); ++c )
         {
            SCIP_VAR* var;

            var = SCIPcolGetVar(cols[c]);
            if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
            {
               SCIP_Real act;

               act = vals[c] * SCIPcolGetPrimsol(cols[c]);
               minactivities[r] += act;
               maxactivities[r] += act;
            }
            else if( vals[c] > 0.0 )
            {
               SCIP_Real lb;
               SCIP_Real ub;

               lb = SCIPvarGetLbGlobal(var);
               ub = SCIPvarGetUbGlobal(var);
               if( SCIPisInfinity(scip, -lb) )
                  mininf = TRUE;
               else
                  minactivities[r] += vals[c] * lb;
               if( SCIPisInfinity(scip, ub) )
                  maxinf = TRUE;
               else
                  maxactivities[r] += vals[c] * ub;
            }
            else if( vals[c] < 0.0 )
            {
               SCIP_Real lb;
               SCIP_Real ub;

               lb = SCIPvarGetLbGlobal(var);
               ub = SCIPvarGetUbGlobal(var);
               if( SCIPisInfinity(scip, ub) )
                  mininf = TRUE;
               else
                  minactivities[r] += vals[c] * ub;
               if( SCIPisInfinity(scip, -lb) )
                  maxinf = TRUE;
               else
                  maxactivities[r] += vals[c] * lb;
            }

            if( mininf )
               minactivities[r] = -SCIPinfinity(scip);
            if( maxinf )
               maxactivities[r] = SCIPinfinity(scip);
         }

         if( SCIPisFeasLT(scip, maxactivities[r], SCIProwGetLhs(row))
            || SCIPisFeasGT(scip, minactivities[r], SCIProwGetRhs(row)) )
         {
            violrows[nviolrows] = row;
            violrowpos[r] = nviolrows;
            nviolrows++;
         }
         else
            violrowpos[r] = -1;
      }
      else
         /* if row is a local row */
         violrowpos[r] = -1;
   }

   nviolfracrows = 0;
   /* calc the current number of fractional variables in rows */
   for( c = 0; c < nlpcands; ++c )
      addFracCounter(nfracsinrow, violrows, violrowpos, &nviolfracrows, nviolrows, nlprows, lpcands[c], +1);

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
      bestshiftval = obj > 0.0 ? SCIPfeasFloor(scip, lpcandssol[c]) : SCIPfeasCeil(scip, lpcandssol[c]);
      minobj += obj * (bestshiftval - lpcandssol[c]);
   }

   /* try to shift remaining variables in order to become/stay feasible */
   nnonimprovingshifts = 0;
   minnviolrows = INT_MAX;
   increaseweight = 1.0;
   while( (nfrac > 0 || nviolrows > 0) && nnonimprovingshifts < MAXSHIFTINGS && !SCIPisStopped(scip) )
   {
      SCIP_VAR* shiftvar;
      SCIP_Real oldsolval;
      SCIP_Real newsolval;
      SCIP_Bool oldsolvalisfrac;
      int probindex;

      SCIPdebugMsg(scip, "intshifting heuristic: nfrac=%d, nviolrows=%d, obj=%g (best possible obj: %g), cutoff=%g\n",
         nfrac, nviolrows, SCIPgetSolOrigObj(scip, sol), SCIPretransformObj(scip, minobj),
         SCIPretransformObj(scip, SCIPgetCutoffbound(scip)));

      nprevviolrows = nviolrows;

      /* choose next variable to process:
       *  - if a violated row exists, shift a variable decreasing the violation, that has least impact on other rows
       *  - otherwise, shift a variable, that has strongest devastating impact on rows in opposite direction
       */
      shiftvar = NULL;
      oldsolval = 0.0;
      newsolval = 0.0;
      if( nviolrows > 0 && (nfrac == 0 || nnonimprovingshifts < MAXSHIFTINGS-1) )
      {
         SCIP_ROW* row;
         int rowidx;
         int rowpos;
         int direction;

         assert(nviolfracrows == 0 || nfrac > 0);
         /* violated rows containing fractional variables are preferred; if such a row exists, choose the last one from the list
          * (at position nviolfracrows - 1) because removing this row will cause one swapping operation less than other rows
          */
         if( nviolfracrows > 0 )
            rowidx = nviolfracrows - 1;
         else
            rowidx = SCIPrandomGetInt(heurdata->randnumgen, 0, nviolrows-1);

         assert(0 <= rowidx && rowidx < nviolrows);
         row = violrows[rowidx];
         rowpos = SCIProwGetLPPos(row);
         assert(0 <= rowpos && rowpos < nlprows);
         assert(violrowpos[rowpos] == rowidx);
         assert(nfracsinrow[rowpos] == 0 || rowidx == nviolfracrows - 1);

         SCIPdebugMsg(scip, "intshifting heuristic: try to fix violated row <%s>: %g <= [%g,%g] <= %g\n",
            SCIProwGetName(row), SCIProwGetLhs(row), minactivities[rowpos], maxactivities[rowpos], SCIProwGetRhs(row));
         SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

         /* get direction in which activity must be shifted */
         assert(SCIPisFeasLT(scip, maxactivities[rowpos], SCIProwGetLhs(row))
            || SCIPisFeasGT(scip, minactivities[rowpos], SCIProwGetRhs(row)));
         direction = SCIPisFeasLT(scip, maxactivities[rowpos], SCIProwGetLhs(row)) ? +1 : -1;

         /* search an integer variable that can shift the activity in the necessary direction */
         SCIP_CALL( selectShifting(scip, sol, row, direction == +1 ? maxactivities[rowpos] : minactivities[rowpos],
               direction, nincreases, ndecreases, increaseweight, &shiftvar, &oldsolval, &newsolval) );
      }

      if( shiftvar == NULL && nfrac > 0 )
      {
         SCIPdebugMsg(scip, "intshifting heuristic: search rounding variable and try to stay feasible\n");
         SCIP_CALL( selectEssentialRounding(scip, sol, minobj, lpcands, nlpcands, &shiftvar, &oldsolval, &newsolval) );
      }

      /* check, whether shifting was possible */
      if( shiftvar == NULL || SCIPisEQ(scip, oldsolval, newsolval) )
      {
         SCIPdebugMsg(scip, "intshifting heuristic:  -> didn't find a shifting variable\n");
         break;
      }

      assert(SCIPvarGetType(shiftvar) == SCIP_VARTYPE_BINARY || SCIPvarGetType(shiftvar) == SCIP_VARTYPE_INTEGER);

      SCIPdebugMsg(scip, "intshifting heuristic:  -> shift var <%s>[%g,%g], type=%d, oldval=%g, newval=%g, obj=%g\n",
         SCIPvarGetName(shiftvar), SCIPvarGetLbGlobal(shiftvar), SCIPvarGetUbGlobal(shiftvar), SCIPvarGetType(shiftvar),
         oldsolval, newsolval, SCIPvarGetObj(shiftvar));

      /* update row activities of globally valid rows */
      SCIP_CALL( updateActivities(scip, minactivities, maxactivities, violrows, violrowpos, &nviolrows, &nviolfracrows,
            nfracsinrow, nlprows, shiftvar, oldsolval, newsolval) );

      if( nviolrows >= nprevviolrows )
         nnonimprovingshifts++;
      else if( nviolrows < minnviolrows )
      {
         minnviolrows = nviolrows;
         nnonimprovingshifts = 0;
      }

      /* store new solution value and decrease fractionality counter */
      SCIP_CALL( SCIPsetSolVal(scip, sol, shiftvar, newsolval) );

      /* update fractionality counter and minimal objective value possible after shifting remaining variables */
      oldsolvalisfrac = !SCIPisFeasIntegral(scip, oldsolval);
      obj = SCIPvarGetObj(shiftvar);
      if( oldsolvalisfrac )
      {
         assert(SCIPisFeasIntegral(scip, newsolval));
         nfrac--;
         nnonimprovingshifts = 0;
         minnviolrows = INT_MAX;
         addFracCounter(nfracsinrow, violrows, violrowpos, &nviolfracrows, nviolrows, nlprows, shiftvar, -1);

         /* the rounding was already calculated into the minobj -> update only if rounding in "wrong" direction */
         if( obj > 0.0 && newsolval > oldsolval )
            minobj += obj;
         else if( obj < 0.0 && newsolval < oldsolval )
            minobj -= obj;
      }
      else
      {
         /* update minimal possible objective value */
         minobj += obj * (newsolval - oldsolval);
      }

      /* update increase/decrease arrays */
      if( !oldsolvalisfrac )
      {
         probindex = SCIPvarGetProbindex(shiftvar);
         assert(0 <= probindex && probindex < nvars);
         increaseweight *= WEIGHTFACTOR;
         if( newsolval < oldsolval )
            ndecreases[probindex] += increaseweight;
         else
            nincreases[probindex] += increaseweight;
         if( increaseweight >= 1e+09 )
         {
            int i;

            for( i = 0; i < nvars; ++i )
            {
               nincreases[i] /= increaseweight;
               ndecreases[i] /= increaseweight;
            }
            increaseweight = 1.0;
         }
      }

      SCIPdebugMsg(scip, "intshifting heuristic:  -> nfrac=%d, nviolrows=%d, obj=%g (best possible obj: %g)\n",
         nfrac, nviolrows, SCIPgetSolOrigObj(scip, sol), SCIPretransformObj(scip, minobj));
   }

   /* check, if the new solution is potentially feasible and solve the LP to calculate values for the continuous
    * variables
    */
   if( nfrac == 0 && nviolrows == 0 )
   {
      SCIP_VAR** vars;
      SCIP_Bool lperror;
      int nintvars;
      int v;
#ifdef NDEBUG
      SCIP_RETCODE retstat;
#endif

      SCIPdebugMsg(scip, "shifted solution is potentially feasible -> solve LP to fix continuous variables\n");

      /* start diving to calculate the LP relaxation */
      SCIP_CALL( SCIPstartDive(scip) );

      /* set the bounds of the variables: fixed for integers, global bounds for continuous */
      vars = SCIPgetVars(scip);
      nintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPchgVarLbDive(scip, vars[v], SCIPvarGetLbGlobal(vars[v])) );
            SCIP_CALL( SCIPchgVarUbDive(scip, vars[v], SCIPvarGetUbGlobal(vars[v])) );
         }
      }
      for( v = 0; v < nintvars; ++v ) /* apply this after global bounds to not cause an error with intermediate empty domains */
      {
         if( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_Real solval;
            solval = SCIPgetSolVal(scip, sol, vars[v]);
            SCIP_CALL( SCIPchgVarLbDive(scip, vars[v], solval) );
            SCIP_CALL( SCIPchgVarUbDive(scip, vars[v], solval) );
         }
      }

      /* solve LP */
      SCIPdebugMsg(scip, " -> old LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));

      /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
       * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
#ifdef NDEBUG
      retstat = SCIPsolveDiveLP(scip, -1, &lperror, NULL);
      if( retstat != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving LP in Intshifting heuristic; LP solve terminated with code <%d>\n",retstat);
      }
#else
      SCIP_CALL( SCIPsolveDiveLP(scip, -1, &lperror, NULL) );
#endif

      SCIPdebugMsg(scip, " -> new LP iterations: %" SCIP_LONGINT_FORMAT "\n", SCIPgetNLPIterations(scip));
      SCIPdebugMsg(scip, " -> error=%u, status=%d\n", lperror, SCIPgetLPSolstat(scip));

      /* check if this is a feasible solution */
      if( !lperror && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         SCIP_Bool stored;

         /* copy the current LP solution to the working solution */
         SCIP_CALL( SCIPlinkLPSol(scip, sol) );

         /* check solution for feasibility, and add it to solution store if possible
          * neither integrality nor feasibility of LP rows has to be checked, because this is already
          * done in the intshifting heuristic itself and due to the LP resolve
          */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, FALSE, &stored) );

         if( stored )
         {
            SCIPdebugMsg(scip, "found feasible shifted solution:\n");
            SCIPdebug( SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) ) );
            *result = SCIP_FOUNDSOL;
         }
      }

      /* terminate the diving */
      SCIP_CALL( SCIPendDive(scip) );
   }

   /* free memory buffers */
   SCIPfreeBufferArray(scip, &ndecreases);
   SCIPfreeBufferArray(scip, &nincreases);
   SCIPfreeBufferArray(scip, &nfracsinrow);
   SCIPfreeBufferArray(scip, &violrowpos);
   SCIPfreeBufferArray(scip, &violrows);
   SCIPfreeBufferArray(scip, &maxactivities);
   SCIPfreeBufferArray(scip, &minactivities);

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the intshifting heuristic with infeasibility recovering and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIntshifting(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEUR* heur;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecIntshifting, NULL) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyIntshifting) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitIntshifting) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitIntshifting) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolIntshifting) );


   return SCIP_OKAY;
}

