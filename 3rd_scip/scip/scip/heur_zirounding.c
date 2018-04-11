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

/**@file   heur_zirounding.c
 * @brief  zirounding primal heuristic
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_zirounding.h"
#include "scip.h"

#define HEUR_NAME             "zirounding"
#define HEUR_DESC             "LP rounding heuristic as suggested by C. Wallace taking row slacks and bounds into account"
#define HEUR_DISPCHAR         'z'
#define HEUR_PRIORITY         -500
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_MAXROUNDINGLOOPS   2      /**< delimits the number of main loops, can be set to -1 for no limit */
#define DEFAULT_STOPZIROUND        TRUE   /**< deactivation check is enabled by default */
#define DEFAULT_STOPPERCENTAGE     0.02   /**< the tolerance percentage after which zirounding will not be executed anymore */
#define DEFAULT_MINSTOPNCALLS      1000   /**< number of heuristic calls before deactivation check */

#define MINSHIFT                   1e-4   /**< minimum shift value for every single step */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Longint          lastlp;             /**< the number of the last LP for which ZIRounding was called */
   int                   maxroundingloops;   /**< limits rounding loops in execution */
   SCIP_Bool             stopziround;        /**< sets deactivation check */
   SCIP_Real             stoppercentage;     /**< threshold for deactivation check */
   int                   minstopncalls;      /**< number of heuristic calls before deactivation check */
};

enum Direction
{
   DIRECTION_UP          =  1,
   DIRECTION_DOWN        = -1
};
typedef enum Direction DIRECTION;

/*
 * Local methods
 */

/** returns the fractionality of a value x, which is calculated as zivalue(x) = min(x-floor(x), ceil(x)-x) */
static
SCIP_Real getZiValue(
   SCIP*                 scip,               /**< pointer to current SCIP data structure */
   SCIP_Real             val                 /**< the value for which the fractionality should be computed */
   )
{
   SCIP_Real upgap;     /* the gap between val and ceil(val) */
   SCIP_Real downgap;   /* the gap between val and floor(val) */

   assert(scip != NULL);

   upgap   = SCIPfeasCeil(scip, val) - val;
   downgap = val - SCIPfeasFloor(scip, val);

   return MIN(upgap, downgap);
}

/** determines shifting bounds for variable */
static
void calculateBounds(
   SCIP*                 scip,               /**< pointer to current SCIP data structure */
   SCIP_VAR*             var,                /**< the variable for which lb and ub have to be calculated */
   SCIP_Real             currentvalue,       /**< the current value of var in the working solution */
   SCIP_Real*            upperbound,         /**< pointer to store the calculated upper bound on the variable shift */
   SCIP_Real*            lowerbound,         /**< pointer to store the calculated lower bound on the variable shift */
   SCIP_Real*            upslacks,           /**< array that contains the slacks between row activities and the right hand sides of the rows */
   SCIP_Real*            downslacks,         /**< array that contains lhs slacks */
   int                   nslacks,            /**< current number of slacks */
   SCIP_Bool*            numericalerror      /**< flag to determine whether a numerical error occurred */
   )
{
   SCIP_COL*      col;
   SCIP_ROW**     colrows;
   SCIP_Real*     colvals;
   int            ncolvals;
   int i;

   assert(scip != NULL);
   assert(var != NULL);
   assert(upslacks != NULL);
   assert(downslacks != NULL);
   assert(upperbound != NULL);
   assert(lowerbound != NULL);

   /* get the column associated to the variable, the nonzero rows and the nonzero coefficients */
   col       = SCIPvarGetCol(var);
   colrows   = SCIPcolGetRows(col);
   colvals   = SCIPcolGetVals(col);
   ncolvals  = SCIPcolGetNLPNonz(col);

   /* only proceed, when variable has nonzero coefficients */
   if( ncolvals == 0 )
      return;

   assert(colvals != NULL);
   assert(colrows != NULL);

   /* initialize the bounds on the shift to be the gap of the current solution value to the bounds of the variable */
   if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
      *upperbound = SCIPinfinity(scip);
   else
      *upperbound = SCIPvarGetUbGlobal(var) - currentvalue;

   if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
      *lowerbound = SCIPinfinity(scip);
   else
      *lowerbound = currentvalue - SCIPvarGetLbGlobal(var);

   /* go through every nonzero row coefficient corresponding to var to determine bounds for shifting
    * in such a way that shifting maintains feasibility in every LP row.
    * a lower or upper bound as it is calculated in zirounding always has to be >= 0.0.
    * if one of these values is significantly < 0.0, this will cause the abort of execution of the heuristic so that
    * infeasible solutions are avoided
    */
   for( i = 0; i < ncolvals && MAX(*lowerbound, *upperbound) >= MINSHIFT; ++i )
   {
      SCIP_ROW* row;
      int       rowpos;

      row = colrows[i];
      rowpos = SCIProwGetLPPos(row);

      /* the row might currently not be in the LP, ignore it! */
      if( rowpos == -1 )
         continue;

      assert(0 <= rowpos && rowpos < nslacks);

      /* all bounds and slacks as they are calculated in zirounding always have to be greater equal zero.
       * It might however be due to numerical issues, e.g. with scaling, that they are not. Better abort in this case.
       */
      if( SCIPisFeasLT(scip, *lowerbound, 0.0) || SCIPisFeasLT(scip, *upperbound, 0.0)
         || SCIPisFeasLT(scip, upslacks[rowpos], 0.0) || SCIPisFeasLT(scip, downslacks[rowpos] , 0.0) )
      {
         *numericalerror = TRUE;
         return;
      }

      SCIPdebugMsg(scip, "colval: %15.8g, downslack: %15.8g, upslack: %5.2g, lb: %5.2g, ub: %5.2g\n", colvals[i], downslacks[rowpos], upslacks[rowpos],
         *lowerbound, *upperbound);

      /* if coefficient > 0, rounding up might violate up slack and rounding down might violate down slack
       * thus search for the minimum so that no constraint is violated; vice versa for coefficient < 0
       */
      if( colvals[i] > 0 )
      {
         if( !SCIPisInfinity(scip, upslacks[rowpos]) )
         {
            SCIP_Real upslack;
            upslack = MAX(upslacks[rowpos], 0.0); /* avoid errors due to numerically slightly infeasible rows */
            *upperbound = MIN(*upperbound, upslack/colvals[i]);
         }

         if( !SCIPisInfinity(scip, downslacks[rowpos]) )
         {
            SCIP_Real downslack;
            downslack = MAX(downslacks[rowpos], 0.0); /* avoid errors due to numerically slightly infeasible rows */
            *lowerbound = MIN(*lowerbound, downslack/colvals[i]);
         }
      }
      else
      {
         assert(colvals[i] != 0.0);

         if( !SCIPisInfinity(scip, upslacks[rowpos]) )
         {
            SCIP_Real upslack;
            upslack = MAX(upslacks[rowpos], 0.0); /* avoid errors due to numerically slightly infeasible rows */
            *lowerbound = MIN(*lowerbound, -upslack/colvals[i]);
         }

         if( !SCIPisInfinity(scip, downslacks[rowpos]) )
         {
            SCIP_Real downslack;
            downslack = MAX(downslacks[rowpos], 0.0); /* avoid errors due to numerically slightly infeasible rows */
            *upperbound = MIN(*upperbound, -downslack/colvals[i]);
         }
      }
   }
}

/**  when a variable is shifted, the activities and slacks of all rows it appears in have to be updated */
static
SCIP_RETCODE updateSlacks(
   SCIP*                 scip,               /**< pointer to current SCIP data structure */
   SCIP_SOL*             sol,                /**< working solution */
   SCIP_VAR*             var,                /**< pointer to variable to be modified */
   SCIP_Real             shiftvalue,         /**< the value by which the variable is shifted */
   SCIP_Real*            upslacks,           /**< upslacks of all rows the variable appears in */
   SCIP_Real*            downslacks,         /**< downslacks of all rows the variable appears in */
   SCIP_Real*            activities,         /**< activities of the LP rows */
   SCIP_VAR**            slackvars,          /**< the slack variables for equality rows */
   SCIP_Real*            slackcoeffs,        /**< the slack variable coefficients */
   int                   nslacks             /**< size of the arrays */
   )
{
   SCIP_COL*    col;        /* the corresponding column of variable var */
   SCIP_ROW**   rows;       /* pointer to the nonzero coefficient rows for variable var */
   int          nrows;      /* the number of nonzeros */
   SCIP_Real*   colvals;    /* array to store the nonzero coefficients */
   int i;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(var != NULL);
   assert(upslacks != NULL);
   assert(downslacks != NULL);
   assert(activities != NULL);
   assert(nslacks >= 0);

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   rows     = SCIPcolGetRows(col);
   nrows    = SCIPcolGetNLPNonz(col);
   colvals  = SCIPcolGetVals(col);
   assert(nrows == 0 || (rows != NULL && colvals != NULL));

   /* go through all rows the shifted variable appears in */
   for( i = 0; i < nrows; ++i )
   {
      int rowpos;

      rowpos = SCIProwGetLPPos(rows[i]);
      assert(-1 <= rowpos && rowpos < nslacks);

      /* if the row is in the LP, update its activity, up and down slack */
      if( rowpos >= 0 )
      {
         SCIP_Real val;

         val = colvals[i] * shiftvalue;

         /* if the row is an equation, we update its slack variable instead of its activities */
         if( SCIPisFeasEQ(scip, SCIProwGetLhs(rows[i]), SCIProwGetRhs(rows[i])) )
         {
            SCIP_Real slackvarshiftval;
            SCIP_Real slackvarsolval;

            assert(slackvars[rowpos] != NULL);
            assert(!SCIPisZero(scip, slackcoeffs[rowpos]));

            slackvarsolval = SCIPgetSolVal(scip, sol, slackvars[rowpos]);
            slackvarshiftval = -val / slackcoeffs[rowpos];

            assert(SCIPisFeasGE(scip, slackvarsolval + slackvarshiftval, SCIPvarGetLbGlobal(slackvars[rowpos])));
            assert(SCIPisFeasLE(scip, slackvarsolval + slackvarshiftval, SCIPvarGetUbGlobal(slackvars[rowpos])));

            SCIP_CALL( SCIPsetSolVal(scip, sol, slackvars[rowpos], slackvarsolval + slackvarshiftval) );
         }
         else if( !SCIPisInfinity(scip, -activities[rowpos]) && !SCIPisInfinity(scip, activities[rowpos]) )
            activities[rowpos] += val;

         /* the slacks of the row now can be updated independently of its type */
         if( !SCIPisInfinity(scip, upslacks[rowpos]) )
            upslacks[rowpos] -= val;
         if( !SCIPisInfinity(scip, -downslacks[rowpos]) )
            downslacks[rowpos] += val;

         assert(!SCIPisFeasNegative(scip, upslacks[rowpos]));
         assert(!SCIPisFeasNegative(scip, downslacks[rowpos]));
      }
   }
   return SCIP_OKAY;
}

/** finds a continuous slack variable for an equation row, NULL if none exists */
static
void rowFindSlackVar(
   SCIP*                 scip,               /**< pointer to current SCIP data structure */
   SCIP_ROW*             row,                /**< the row for which a slack variable is searched */
   SCIP_VAR**            varpointer,         /**< pointer to store the slack variable */
   SCIP_Real*            coeffpointer        /**< pointer to store the coefficient of the slack variable */
   )
{
   int v;
   SCIP_COL** rowcols;
   SCIP_Real* rowvals;
   int nrowvals;

   assert(row != NULL);
   assert(varpointer != NULL);
   assert(coeffpointer != NULL);

   rowcols = SCIProwGetCols(row);
   rowvals = SCIProwGetVals(row);
   nrowvals = SCIProwGetNNonz(row);

   assert(nrowvals == 0 || rowvals != NULL);
   assert(nrowvals == 0 || rowcols != NULL);

   /* iterate over the row variables. Stop after the first unfixed continuous variable was found. */
   for( v = nrowvals - 1; v >= 0; --v )
   {
      SCIP_VAR* colvar;

      assert(rowcols[v] != NULL);
      if( SCIPcolGetLPPos(rowcols[v]) == -1 )
         continue;

      colvar = SCIPcolGetVar(rowcols[v]);

      if( SCIPvarGetType(colvar) == SCIP_VARTYPE_CONTINUOUS
         && !SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(colvar), SCIPvarGetUbGlobal(colvar))
         && SCIPcolGetNLPNonz(rowcols[v]) == 1 )
      {
         SCIPdebugMsg(scip, "  slack variable for row %s found: %s\n", SCIProwGetName(row), SCIPvarGetName(colvar));

         *coeffpointer = rowvals[v];
         *varpointer = colvar;

         return;
      }
   }

   *varpointer = NULL;
   *coeffpointer = 0.0;

   SCIPdebugMsg(scip, "No slack variable for row %s found. \n", SCIProwGetName(row));
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyZirounding)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurZirounding(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeZirounding)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitZirounding)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitZirounding)  /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolZirounding)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->lastlp = -1;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecZirounding)
{  /*lint --e{715}*/
   SCIP_HEURDATA*     heurdata;
   SCIP_SOL*          sol;
   SCIP_VAR**         lpcands;
   SCIP_VAR**         zilpcands;

   SCIP_VAR**         slackvars;
   SCIP_Real*         upslacks;
   SCIP_Real*         downslacks;
   SCIP_Real*         activities;
   SCIP_Real*         slackvarcoeffs;
   SCIP_Bool*         rowneedsslackvar;

   SCIP_ROW**         rows;
   SCIP_Real*         lpcandssol;
   SCIP_Real*         solarray;

   SCIP_Longint       nlps;
   int                currentlpcands;
   int                nlpcands;
   int                nimplfracs;
   int                i;
   int                c;
   int                nslacks;
   int                nroundings;

   SCIP_Bool          improvementfound;
   SCIP_Bool          numericalerror;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic if an optimal LP-solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* Do not call heuristic if deactivation check is enabled and percentage of found solutions in relation
    * to number of calls falls below heurdata->stoppercentage */
   if( heurdata->stopziround && SCIPheurGetNCalls(heur) >= heurdata->minstopncalls
      && SCIPheurGetNSolsFound(heur)/(SCIP_Real)SCIPheurGetNCalls(heur) < heurdata->stoppercentage )
      return SCIP_OKAY;

   /* assure that heuristic has not already been called after the last LP had been solved */
   nlps = SCIPgetNLPs(scip);
   if( nlps == heurdata->lastlp )
      return SCIP_OKAY;

   heurdata->lastlp = nlps;

   /* get fractional variables */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, &nimplfracs) );
   nlpcands = nlpcands + nimplfracs;
   /* make sure that there is at least one fractional variable that should be integral */
   if( nlpcands == 0 )
      return SCIP_OKAY;

   assert(nlpcands > 0);
   assert(lpcands != NULL);
   assert(lpcandssol != NULL);

   /* get LP rows data */
   rows    = SCIPgetLPRows(scip);
   nslacks = SCIPgetNLPRows(scip);

   /* cannot do anything if LP is empty */
   if( nslacks == 0 )
      return SCIP_OKAY;

   assert(rows != NULL);
   assert(nslacks > 0);

   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert(sol != NULL);

   *result = SCIP_DIDNOTFIND;

   solarray = NULL;
   zilpcands = NULL;

   /* copy the current LP solution to the working solution and allocate memory for local data */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solarray, nlpcands) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zilpcands, nlpcands) );

   /* copy necessary data to local arrays */
   BMScopyMemoryArray(solarray, lpcandssol, nlpcands);
   BMScopyMemoryArray(zilpcands, lpcands, nlpcands);

   /* allocate buffer data arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &slackvars, nslacks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &upslacks, nslacks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &downslacks, nslacks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &slackvarcoeffs, nslacks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowneedsslackvar, nslacks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activities, nslacks) );

   BMSclearMemoryArray(slackvars, nslacks);
   BMSclearMemoryArray(slackvarcoeffs, nslacks);
   BMSclearMemoryArray(rowneedsslackvar, nslacks);

   numericalerror = FALSE;
   nroundings = 0;

   /* loop over fractional variables and involved LP rows to find all rows which require a slack variable */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* cand;
      SCIP_ROW** candrows;
      int r;
      int ncandrows;

      cand = zilpcands[c];
      assert(cand != NULL);
      assert(SCIPcolGetLPPos(SCIPvarGetCol(cand)) >= 0);

      candrows = SCIPcolGetRows(SCIPvarGetCol(cand));
      ncandrows = SCIPcolGetNLPNonz(SCIPvarGetCol(cand));

      assert(candrows == NULL || ncandrows > 0);

      for( r = 0; r < ncandrows; ++r )
      {
         int rowpos;

         assert(candrows != NULL); /* to please flexelint */
         assert(candrows[r] != NULL);
         rowpos = SCIProwGetLPPos(candrows[r]);

         if( rowpos >= 0 && SCIPisFeasEQ(scip, SCIProwGetLhs(candrows[r]), SCIProwGetRhs(candrows[r])) )
         {
            rowneedsslackvar[rowpos] = TRUE;
            SCIPdebugMsg(scip, "  Row %s needs slack variable for variable %s\n", SCIProwGetName(candrows[r]), SCIPvarGetName(cand));
         }
      }
   }

   /* calculate row slacks for every every row that belongs to the current LP and ensure, that the current solution
    * has no violated constraint -- if any constraint is violated, i.e. a slack is significantly smaller than zero,
    * this will cause the termination of the heuristic because Zirounding does not provide feasibility recovering
    */
   for( i = 0; i < nslacks; ++i )
   {
      SCIP_ROW*          row;
      SCIP_Real          lhs;
      SCIP_Real          rhs;

      row = rows[i];

      assert(row != NULL);

      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);

      /* get row activity */
      activities[i] = SCIPgetRowActivity(scip, row);
      assert(SCIPisFeasLE(scip, lhs, activities[i]) && SCIPisFeasLE(scip, activities[i], rhs));

      /* in special case if LHS or RHS is (-)infinity slacks have to be initialized as infinity */
      if( SCIPisInfinity(scip, -lhs) )
         downslacks[i] = SCIPinfinity(scip);
      else
         downslacks[i] = activities[i] - lhs;

      if( SCIPisInfinity(scip, rhs) )
         upslacks[i] = SCIPinfinity(scip);
      else
         upslacks[i] = rhs - activities[i];

      SCIPdebugMsg(scip, "lhs:%5.2f <= act:%5.2g <= rhs:%5.2g --> down: %5.2g, up:%5.2g\n", lhs, activities[i], rhs, downslacks[i], upslacks[i]);

      /* row is an equation. Try to find a slack variable in the row, i.e.,
       * a continuous variable which occurs only in this row. If no such variable exists,
       * there is no hope for an IP-feasible solution in this round
       */
      if( SCIPisFeasEQ(scip, lhs, rhs) && rowneedsslackvar[i] )
      {
         /* @todo: This is only necessary for rows containing fractional variables. */
         rowFindSlackVar(scip, row, &(slackvars[i]), &(slackvarcoeffs[i]));

         if( slackvars[i] == NULL )
         {
            SCIPdebugMsg(scip, "No slack variable found for equation %s, terminating ZI Round heuristic\n", SCIProwGetName(row));
            goto TERMINATE;
         }
         else
         {
            SCIP_Real ubslackvar;
            SCIP_Real lbslackvar;
            SCIP_Real solvalslackvar;
            SCIP_Real coeffslackvar;
            SCIP_Real ubgap;
            SCIP_Real lbgap;

            assert(SCIPvarGetType(slackvars[i]) == SCIP_VARTYPE_CONTINUOUS);
            solvalslackvar = SCIPgetSolVal(scip, sol, slackvars[i]);
            ubslackvar = SCIPvarGetUbGlobal(slackvars[i]);
            lbslackvar = SCIPvarGetLbGlobal(slackvars[i]);

            coeffslackvar = slackvarcoeffs[i];
            assert(!SCIPisFeasZero(scip, coeffslackvar));

            ubgap = MAX(0.0, ubslackvar - solvalslackvar);
            lbgap = MAX(0.0, solvalslackvar - lbslackvar);

            if( SCIPisPositive(scip, coeffslackvar) )
            {
              if( !SCIPisInfinity(scip, lbslackvar) )
                upslacks[i] += coeffslackvar * lbgap;
              else
                upslacks[i] = SCIPinfinity(scip);
              if( !SCIPisInfinity(scip, ubslackvar) )
                downslacks[i] += coeffslackvar * ubgap;
              else
                downslacks[i] = SCIPinfinity(scip);
            }
            else
            {
               if( !SCIPisInfinity(scip, ubslackvar) )
                  upslacks[i] -= coeffslackvar * ubgap;
               else
                  upslacks[i] = SCIPinfinity(scip);
               if( !SCIPisInfinity(scip, lbslackvar) )
                  downslacks[i] -= coeffslackvar * lbgap;
               else
                  downslacks[i] = SCIPinfinity(scip);
            }
            SCIPdebugMsg(scip, "  Slack variable for row %s at pos %d: %g <= %s = %g <= %g; Coeff %g, upslack = %g, downslack = %g  \n",
               SCIProwGetName(row), SCIProwGetLPPos(row), lbslackvar, SCIPvarGetName(slackvars[i]), solvalslackvar, ubslackvar, coeffslackvar,
               upslacks[i], downslacks[i]);
         }
      }
      /* due to numerical inaccuracies, the rows might be feasible, even if the slacks are
       * significantly smaller than zero -> terminate
       */
      if( SCIPisFeasLT(scip, upslacks[i], 0.0) || SCIPisFeasLT(scip, downslacks[i], 0.0) )
         goto TERMINATE;
   }

   assert(nslacks == 0 || (upslacks != NULL && downslacks != NULL && activities != NULL));

   /* initialize number of remaining variables and flag to enter the main loop */
   currentlpcands = nlpcands;
   improvementfound = TRUE;

   /* iterate over variables as long as there are fractional variables left */
   while( currentlpcands > 0 && improvementfound && (heurdata->maxroundingloops == -1 || nroundings < heurdata->maxroundingloops) )
   {  /*lint --e{850}*/
      improvementfound = FALSE;
      nroundings++;
      SCIPdebugMsg(scip, "zirounding enters while loop for %d time with %d candidates left. \n", nroundings, currentlpcands);

      /* check for every remaining fractional variable if a shifting decreases ZI-value of the variable */
      for( c = 0; c < currentlpcands; ++c )
      {
         SCIP_VAR* var;
         SCIP_Real oldsolval;
         SCIP_Real upperbound;
         SCIP_Real lowerbound;
         SCIP_Real up;
         SCIP_Real down;
         SCIP_Real ziup;
         SCIP_Real zidown;
         SCIP_Real zicurrent;
         SCIP_Real shiftval;

         DIRECTION direction;

         /* get values from local data */
         oldsolval = solarray[c];
         var = zilpcands[c];

         assert(!SCIPisFeasIntegral(scip, oldsolval));
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

         /* calculate bounds for variable and make sure that there are no numerical inconsistencies */
         upperbound = SCIPinfinity(scip);
         lowerbound = SCIPinfinity(scip);
         calculateBounds(scip, var, oldsolval, &upperbound, &lowerbound, upslacks, downslacks, nslacks, &numericalerror);

         if( numericalerror )
            goto TERMINATE;

         /* continue if only marginal shifts are possible */
         if( MAX(upperbound, lowerbound) < MINSHIFT )
         {
            /* stop immediately if a variable has not been rounded during the last round */
            if( nroundings == heurdata->maxroundingloops )
               break;
            else
               continue;
         }

         /* calculate the possible values after shifting */
         up   = oldsolval + upperbound;
         down = oldsolval - lowerbound;

         /* if the variable is integer or implicit binary, do not shift further than the nearest integer */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY)
         {
            SCIP_Real ceilx;
            SCIP_Real floorx;

            ceilx = SCIPfeasCeil(scip, oldsolval);
            floorx = SCIPfeasFloor(scip, oldsolval);
            up   = MIN(up, ceilx);
            down = MAX(down, floorx);
         }

         /* calculate necessary values */
         ziup      = getZiValue(scip, up);
         zidown    = getZiValue(scip, down);
         zicurrent = getZiValue(scip, oldsolval);

         /* calculate the shifting direction that reduces ZI-value the most,
          * if both directions improve ZI-value equally, take the direction which improves the objective
          */
         if( SCIPisFeasLT(scip, zidown, zicurrent) || SCIPisFeasLT(scip, ziup, zicurrent) )
         {
            if( SCIPisFeasEQ(scip,ziup, zidown) )
               direction  = SCIPisFeasGE(scip, SCIPvarGetObj(var), 0.0) ? DIRECTION_DOWN : DIRECTION_UP;
            else if( SCIPisFeasLT(scip, zidown, ziup) )
               direction = DIRECTION_DOWN;
            else
               direction = DIRECTION_UP;

            /* once a possible shifting direction and value have been found, variable value is updated */
            shiftval = (direction == DIRECTION_UP ? up - oldsolval : down - oldsolval);

            /* this improves numerical stability in some cases */
            if( direction == DIRECTION_UP )
               shiftval = MIN(shiftval, upperbound);
            else
               shiftval = MIN(shiftval, lowerbound);
            /* update the solution */
            solarray[c] = direction == DIRECTION_UP ? up : down;
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, solarray[c]) );

            /* update the rows activities and slacks */
            SCIP_CALL( updateSlacks(scip, sol, var, shiftval, upslacks,
                  downslacks, activities, slackvars, slackvarcoeffs, nslacks) );

            SCIPdebugMsg(scip, "zirounding update step : %d var index, oldsolval=%g, shiftval=%g\n",
               SCIPvarGetIndex(var), oldsolval, shiftval);
            /* since at least one improvement has been found, heuristic will enter main loop for another time because the improvement
             * might affect many LP rows and their current slacks and thus make further rounding steps possible */
            improvementfound = TRUE;
         }

         /* if solution value of variable has become feasibly integral due to rounding step,
          * variable is put at the end of remaining candidates array so as not to be considered in future loops
          */
         if( SCIPisFeasIntegral(scip, solarray[c]) )
         {
            zilpcands[c] = zilpcands[currentlpcands - 1];
            solarray[c] = solarray[currentlpcands - 1];
            currentlpcands--;

            /* counter is decreased if end of candidates array has not been reached yet */
            if( c < currentlpcands )
               c--;
         }
         else if( nroundings == heurdata->maxroundingloops )
            goto TERMINATE;
      }
   }

   /* in case that no candidate is left for rounding after the final main loop
    * the found solution has to be checked for feasibility in the original problem
    */
   if( currentlpcands == 0 )
   {
      SCIP_Bool stored;
      SCIP_CALL(SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, TRUE, FALSE, &stored));
      if( stored )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "found feasible rounded solution:\n");
         SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
         SCIPstatisticMessage("  ZI Round solution value: %g \n", SCIPgetSolOrigObj(scip, sol));

         *result = SCIP_FOUNDSOL;
      }
   }

   /* free memory for all locally allocated data */
 TERMINATE:
   SCIPfreeBufferArrayNull(scip, &activities);
   SCIPfreeBufferArrayNull(scip, &rowneedsslackvar);
   SCIPfreeBufferArrayNull(scip, &slackvarcoeffs);
   SCIPfreeBufferArrayNull(scip, &downslacks);
   SCIPfreeBufferArrayNull(scip, &upslacks);
   SCIPfreeBufferArrayNull(scip, &slackvars);
   SCIPfreeBufferArrayNull(scip, &zilpcands);
   SCIPfreeBufferArrayNull(scip, &solarray);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the zirounding primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurZirounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create zirounding primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecZirounding, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyZirounding) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeZirounding) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitZirounding) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitZirounding) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolZirounding) );

   /* add zirounding primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/zirounding/maxroundingloops",
         "determines maximum number of rounding loops",
         &heurdata->maxroundingloops, TRUE, DEFAULT_MAXROUNDINGLOOPS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/zirounding/stopziround",
         "flag to determine if Zirounding is deactivated after a certain percentage of unsuccessful calls",
         &heurdata->stopziround, TRUE, DEFAULT_STOPZIROUND, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,"heuristics/zirounding/stoppercentage",
         "if percentage of found solutions falls below this parameter, Zirounding will be deactivated",
         &heurdata->stoppercentage, TRUE, DEFAULT_STOPPERCENTAGE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/zirounding/minstopncalls",
         "determines the minimum number of calls before percentage-based deactivation of Zirounding is applied",
         &heurdata->minstopncalls, TRUE, DEFAULT_MINSTOPNCALLS, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
