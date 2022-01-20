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

/**@file   heur_simplerounding.c
 * @brief  simple and fast LP rounding heuristic
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 *
 * The heuristic also tries to round relaxation solutions if available.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_simplerounding.h"


#define HEUR_NAME             "simplerounding"
#define HEUR_DESC             "simple and fast LP rounding heuristic"
#define HEUR_DISPCHAR         'r'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_DURINGPRICINGLOOP
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_ONCEPERNODE   FALSE          /**< should the heuristic only be called once per node? */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Longint          lastlp;             /**< last LP number where the heuristic was applied */
   int                   nroundablevars;     /**< number of variables that can be rounded (-1 if not yet calculated) */
   SCIP_Bool             oncepernode;        /**< should the heuristic only be called once per node? */
};


/*
 * Local methods
 */

/** perform rounding */
static
SCIP_RETCODE performSimpleRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_SOL*             sol,                /**< solution to round */
   SCIP_VAR**            cands,              /**< candidate variables */
   SCIP_Real*            candssol,           /**< solutions of candidate variables */
   int                   ncands,             /**< number of candidates */
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   int c;
   int nunroundableimplints = 0;

   /* round all roundable fractional columns in the corresponding direction as long as no unroundable column was found */
   for (c = 0; c < ncands; ++c)
   {
      SCIP_VAR* var;
      SCIP_Real oldsolval;
      SCIP_Real newsolval;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;

      oldsolval = candssol[c];
      assert( ! SCIPisFeasIntegral(scip, oldsolval) );
      var = cands[c];
      assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      SCIPdebugMsg(scip, "simple rounding heuristic: var <%s>, val=%g, rounddown=%u, roundup=%u\n",
         SCIPvarGetName(var), oldsolval, mayrounddown, mayroundup);

      /* choose rounding direction */
      if ( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if ( SCIPvarGetObj(var) >= 0.0 )
            newsolval = SCIPfeasFloor(scip, oldsolval);
         else
            newsolval = SCIPfeasCeil(scip, oldsolval);
      }
      else if ( mayrounddown )
         newsolval = SCIPfeasFloor(scip, oldsolval);
      else if ( mayroundup )
         newsolval = SCIPfeasCeil(scip, oldsolval);
      else if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
      {
         ++nunroundableimplints;
         continue;
      }
      else
         break;

      /* store new solution value */
      SCIP_CALL( SCIPsetSolVal(scip, sol, var, newsolval) );
   }

   /* check, if rounding was successful */
   if( c == ncands )
   {
      SCIP_Bool stored;
      SCIP_Bool checklprows;

      /* unroundable implicit integers are adjusted. LP rows must be checked afterwards */
      if( nunroundableimplints > 0 )
      {
         SCIP_CALL( SCIPadjustImplicitSolVals(scip, sol, TRUE) );
         checklprows = TRUE;
      }
      else
         checklprows = FALSE;

      if( SCIPallColsInLP(scip) )
      {
         /* check solution for feasibility, and add it to solution store if possible
          * integrality need not be checked, because all fractional
          * variables were already moved in feasible direction to the next integer
          *
          * feasibility of LP rows must be checked again at the presence of
          * unroundable, implicit integer variables with fractional LP solution
          * value
          */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, checklprows, &stored) );
      }
      else
      {
         /* if there are variables which are not present in the LP, e.g., for 
          * column generation, we need to check their bounds
          */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, FALSE, checklprows, &stored) );
      }

      if( stored )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "found feasible rounded solution:\n");
         SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
         *result = SCIP_FOUNDSOL;
      }
   }
   return SCIP_OKAY;
}

/** perform LP-rounding */
static
SCIP_RETCODE performLPSimpleRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Longint nlps;
   int nlpcands;
   int nfracimplvars;

   /* only call heuristic, if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get fractional variables, that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL, &nfracimplvars) );

   /* only call heuristic, if LP solution is fractional; except we are called during pricing, in this case we
    * want to detect a (mixed) integer (LP) solution which is primal feasible
    */
   if ( nlpcands == 0  && heurtiming != SCIP_HEURTIMING_DURINGPRICINGLOOP )
      return SCIP_OKAY;

   /* don't call heuristic, if there are more fractional variables than roundable ones. We do not consider
    * fractional implicit integer variables here, because simple rounding may adjust those separately,
    * even if they aren't roundable
    */
   if ( nlpcands > heurdata->nroundablevars )
      return SCIP_OKAY;

   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert( sol != NULL );

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   /* don't call heuristic, if we have already processed the current LP solution */
   nlps = SCIPgetNLPs(scip);
   if( nlps == heurdata->lastlp )
      return SCIP_OKAY;
   heurdata->lastlp = nlps;

   /* perform simple rounding */
   SCIPdebugMsg(scip, "executing simple LP-rounding heuristic, fractionals: %d + %d\n", nlpcands, nfracimplvars);
   SCIP_CALL( performSimpleRounding(scip, sol, lpcands, lpcandssol, nlpcands + nfracimplvars, result) );

   return SCIP_OKAY;
}

/** perform relaxation solution rounding */
static
SCIP_RETCODE performRelaxSimpleRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_VAR** relaxcands;
   SCIP_Real* relaxcandssol;
   int nrelaxcands = 0;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int ndiscretevars;
   int v;

   /* do not call heuristic if no relaxation solution is available */
   if ( ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   /* get variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, &nimplvars, NULL) );
   ndiscretevars = nbinvars + nintvars + nimplvars; /* consider binary, integral, and implicit integer variables */

   /* get storage */
   SCIP_CALL( SCIPallocBufferArray(scip, &relaxcands, ndiscretevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &relaxcandssol, ndiscretevars) );

   /* get fractional variables, that should be integral */
   for (v = 0; v < nbinvars + nintvars; ++v)
   {
      SCIP_Real val;

      val = SCIPgetRelaxSolVal(scip, vars[v]);
      if ( ! SCIPisFeasIntegral(scip, val) )
      {
         relaxcands[nrelaxcands] = vars[v];
         relaxcandssol[nrelaxcands++] = val;
      }
   }

   /* don't call heuristic, if there are more fractional variables than roundable ones. We explicitly
    * do not consider implicit integer variables with fractional relaxation solution here
    * because they may be feasibly adjusted, although they are not roundable
    */
   if ( nrelaxcands > heurdata->nroundablevars )
   {
      SCIPfreeBufferArray(scip, &relaxcands);
      SCIPfreeBufferArray(scip, &relaxcandssol);
      return SCIP_OKAY;
   }

   /* collect implicit integer variables with fractional solution value */
   for( v = nbinvars + nintvars; v < ndiscretevars; ++v )
   {
      SCIP_Real val;

      val = SCIPgetRelaxSolVal(scip, vars[v]);
      if ( ! SCIPisFeasIntegral(scip, val) )
      {
         relaxcands[nrelaxcands] = vars[v];
         relaxcandssol[nrelaxcands++] = val;
      }
   }
   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert( sol != NULL );

   /* copy the current relaxation solution to the working solution */
   SCIP_CALL( SCIPlinkRelaxSol(scip, sol) );

   /* perform simple rounding */
   SCIPdebugMsg(scip, "executing simple rounding heuristic on relaxation solution: %d fractionals\n", nrelaxcands);
   SCIP_CALL( performSimpleRounding(scip, sol, relaxcands, relaxcandssol, nrelaxcands, result) );

   /* free storage */
   SCIPfreeBufferArray(scip, &relaxcands);
   SCIPfreeBufferArray(scip, &relaxcandssol);

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySimplerounding)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSimplerounding(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSimplerounding) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create heuristic data */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   heurdata->lastlp = -1;
   heurdata->nroundablevars = -1;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitSimplerounding) /*lint --e{715}*/
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
SCIP_DECL_HEURINITSOL(heurInitsolSimplerounding)
{
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->lastlp = -1;

   /* change the heuristic's timingmask, if it should be called only once per node */
   if( heurdata->oncepernode )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_AFTERLPNODE);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolSimplerounding)
{
   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand or if relaxation solution is available */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't call heuristic, if we have already processed the current LP solution but no relaxation solution is available */
   if ( SCIPgetNLPs(scip) == heurdata->lastlp && ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   /* on our first call or after each pricing round, calculate the number of roundable variables */
   if( heurdata->nroundablevars == -1  || heurtiming == SCIP_HEURTIMING_DURINGPRICINGLOOP )
   {
      SCIP_VAR** vars;
      int nbinintvars;
      int nroundablevars;
      int i;

      vars = SCIPgetVars(scip);
      nbinintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
      nroundablevars = 0;
      for( i = 0; i < nbinintvars; ++i )
      {
         if( SCIPvarMayRoundDown(vars[i]) || SCIPvarMayRoundUp(vars[i]) )
            nroundablevars++;
      }
      heurdata->nroundablevars = nroundablevars;
   }

   /* don't call heuristic if there are no roundable variables; except we are called during pricing, in this case we
    * want to detect a (mixed) integer (LP) solution which is primal feasible */
   if( heurdata->nroundablevars == 0 && heurtiming != SCIP_HEURTIMING_DURINGPRICINGLOOP )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* try to round LP solution */
   SCIP_CALL( performLPSimpleRounding(scip, heurdata, heurtiming, result) );

   /* try to round relaxation solution */
   SCIP_CALL( performRelaxSimpleRounding(scip, heurdata, result) );

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** creates the simple rounding heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSimplerounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSimplerounding, heurdata) );
   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySimplerounding) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSimplerounding) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitSimplerounding) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolSimplerounding) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolSimplerounding) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSimplerounding) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/oncepernode",
         "should the heuristic only be called once per node?",
         &heurdata->oncepernode, TRUE, DEFAULT_ONCEPERNODE, NULL, NULL) );

   return SCIP_OKAY;
}
