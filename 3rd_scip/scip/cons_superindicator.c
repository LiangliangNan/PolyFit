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

/**@file   cons_superindicator.c
 * @brief  constraint handler for indicator constraints over arbitrary constraint types
 * @author Ambros Gleixner
 * @author Frederic Pythoud
 */

/**@todo allow more types for slack constraint */
/**@todo implement more upgrades, e.g., for nonlinear, quadratic, logicor slack constraints; upgrades could also help to
 *       handle difficult slack constraints such as pseudoboolean or indicator
 */
/**@todo unify enfolp and enfops, sepalp and sepaps callbacks */
/**@todo enforce by branching on binary variable if slack constraint only returns SCIP_INFEASIBLE */
/**@todo consider enforcing by adding slack constraint (or copy of it) locally if binary variable is fixed to 1
 *       (some constraint handler cannot enforce constraints that are not active)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_superindicator.h"
#include "scip/dialog_default.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME                        "superindicator"
#define CONSHDLR_DESC                        "constraint handler for indicator constraints over arbitrary constraint types"
#define CONSHDLR_SEPAPRIORITY              0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -5000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -5000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ                 -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ                  1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ               100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS             -1 /**< maximal number of presolving rounds the constraint handler
                                              *   participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA             FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP             FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS              TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define DEFAULT_CHECKSLACKTYPE          TRUE /**< should type of slack constraint be checked when creating superindicator constraint? */
#define DEFAULT_UPGDPRIOINDICATOR          1 /**< priority for upgrading to an indicator constraint (-1: never) */
#define DEFAULT_UPGDPRIOLINEAR             2 /**< priority for upgrading to a linear constraint (-1: never) */
#define DEFAULT_MAXUPGDCOEFLINEAR        1e4 /**< maximum big-M coefficient of binary variable in upgrade to a linear constraint
                                              *   (relative to smallest coefficient) */


/*
 * Data structures
 */

/** constraint data for superindicator constraints */
struct SCIP_ConsData
{
   SCIP_CONS*            slackcons;          /**< constraint corresponding to the handled constraint */
   SCIP_VAR*             binvar;             /**< binary variable for indicator constraint */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             checkslacktype;     /**< should type of slack constraint be checked when creating superindicator constraint? */
   SCIP_Real             maxupgdcoeflinear;  /**< maximum big-M coefficient of binary variable in upgrade to a linear constraint
                                              *   (relative to smallest coefficient) */
   int                   upgdprioindicator;  /**< priority for upgrading to an indicator constraint (-1: never) */
   int                   upgdpriolinear;     /**< priority for upgrading to a linear constraint (-1: never) */
   int                   nrejects;           /**< number of rejected calls to create method */
};

/*
 * Local methods
 */

/** creates superindicator constraint data */
static
SCIP_RETCODE consdataCreateSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   SCIP_VAR*             binvar,             /**< binary variable */
   SCIP_CONS*            slackcons           /**< slack constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->binvar = binvar;
   (*consdata)->slackcons = slackcons;

   if( SCIPisTransformed(scip) )
   {
      SCIPdebugMsg(scip, "creating the transformed data\n");

      /* do not capture the slack constraint when scip is in transformed mode; this automatically happens in
       * SCIPtransformCons() if necessary
       */
      SCIP_CALL( SCIPtransformCons(scip, (*consdata)->slackcons, &(*consdata)->slackcons) );

      /* get transformed binary variable */
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->binvar, &(*consdata)->binvar) );
   }
   else
   {
      /* we need to capture the constraint to avoid that SCIP deletes them since they are not (yet) added to the problem */
      SCIP_CALL( SCIPcaptureCons(scip, slackcons) );
   }

   assert((*consdata)->slackcons != NULL);

   return SCIP_OKAY;
}

/** checks the feasibility of a superindicator constraint */
static
SCIP_RETCODE consdataCheckSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< pointer to superindicator constraint data */
   SCIP_SOL*             sol,                /**< pointer to the solution to be checked */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the test */
   )
{
   SCIP_Real binval;

   /* not to be called if infeasibility is already detected */
   assert(*result == SCIP_FEASIBLE || *result == SCIP_DIDNOTRUN);

   binval = SCIPgetSolVal(scip, sol, consdata->binvar);

   /* check integrality of binary variable */
   if( checkintegrality && !SCIPisIntegral(scip, binval) )
   {
      if( printreason )
      {
         SCIPinfoMessage(scip, NULL, "violation: binvar takes fractional value %.15g\n", binval);
      }

      *result = SCIP_INFEASIBLE;
   }
   /* if binvar is one, call SCIPcheckCons() for the slack constraint */
   else if( binval > 0.5 )
   {
      assert(SCIPisFeasEQ(scip, binval, 1.0));

      SCIP_CALL( SCIPcheckCons(scip, consdata->slackcons, sol, checkintegrality, checklprows, printreason, result) );

      if( printreason && *result != SCIP_FEASIBLE )
      {
         SCIPinfoMessage(scip, NULL, "violation: SCIPcheckCons() for slack constraint <%s> returns infeasible while binvar <%s> == 1\n",
            SCIPconsGetName(consdata->slackcons), SCIPvarGetName(consdata->binvar));
      }

#ifdef SCIP_DEBUG
      {
         /* checking in debug mode that different flags don't give us different results */
         SCIP_RESULT testresultnotintegrality;
         SCIP_RESULT testresultnotlprows;

         SCIP_CALL( SCIPcheckCons(scip, consdata->slackcons, sol, checkintegrality, TRUE, TRUE, &testresultnotintegrality) );
         SCIP_CALL( SCIPcheckCons(scip, consdata->slackcons, sol, TRUE, checklprows, TRUE, &testresultnotlprows) );

         assert(*result == testresultnotintegrality);
         assert(*result == testresultnotlprows);
      }
#endif

      SCIPdebugMsg(scip, "binvar <%s> == 1, sol=%p --> SCIPcheckCons() on constraint <%s> --> %s\n",
         SCIPvarGetName(consdata->binvar), (void*)sol, SCIPconsGetName(consdata->slackcons),
         *result == SCIP_FEASIBLE ? "satisfied" : "violated");
   }
   /* if binval is zero, the superindicator constraint is feasible */
   else
   {
      *result = SCIP_FEASIBLE;
   }

   return SCIP_OKAY;
}

/** computes the minactivity, maxactivity, and minimal absolute value of nonzero coefficients of a linear constraint
 *  with respect to its global bounds
 */
static
SCIP_RETCODE extractLinearValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< pointer to linear constraint */
   SCIP_Real*            minactivity,        /**< pointer to return the minimal activity */
   SCIP_Real*            maxactivity,        /**< pointer to return the maximal activity */
   SCIP_Real*            minabscoef          /**< pointer to return the minimal absolute value of the coefficients */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool ismininfinity;
   SCIP_Bool ismaxinfinity;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(minactivity != NULL);
   assert(maxactivity != NULL);
   assert(minabscoef != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0);

   /* get nonzero elements */
   vars = SCIPgetVarsLinear(scip, cons);
   vals = SCIPgetValsLinear(scip, cons);
   nvars = SCIPgetNVarsLinear(scip, cons);

   /* initialize values */
   *minactivity = 0.0;
   *maxactivity = 0.0;
   *minabscoef = SCIPinfinity(scip);
   ismininfinity = FALSE;
   ismaxinfinity = FALSE;

   /* we loop over all the coefficients of the constraint and we cannot end if the minactivity is infinite as we
    * still need to compute the minimum absolute coefficient value
    */
   for( i = nvars-1; i >= 0; i-- )
   {
      SCIP_Real val;
      SCIP_Real lb;
      SCIP_Real ub;

      val = vals[i];
      lb = SCIPvarGetLbGlobal(vars[i]);
      ub = SCIPvarGetUbGlobal(vars[i]);

      /* update flags for infinite bounds */
      ismininfinity = ismininfinity || (val > 0.0 && (SCIPisInfinity(scip, lb) || SCIPisInfinity(scip, -lb)))
         || (val < 0.0 && (SCIPisInfinity(scip, ub) || SCIPisInfinity(scip, -ub)));

      ismaxinfinity = ismaxinfinity || (val > 0.0 && (SCIPisInfinity(scip, ub) || SCIPisInfinity(scip, -ub)))
         || (val < 0.0 && (SCIPisInfinity(scip, lb) || SCIPisInfinity(scip, -lb)));

      /* update activities if not infinite */
      if( !ismininfinity )
         *minactivity += (val > 0.0) ? val * lb : val * ub;

      if( !ismaxinfinity )
         *maxactivity += (val > 0.0) ? val * ub : val * lb;

      /* update minimal absolute coefficient value */
      if( val > 0.0 && val < *minabscoef )
         *minabscoef = val;
      else if( val < 0.0 && -val < *minabscoef )
         *minabscoef = -vals[i];
   }

   if( ismininfinity )
      *minactivity = -SCIPinfinity(scip);

   if( ismaxinfinity )
      *maxactivity = SCIPinfinity(scip);

   return SCIP_OKAY;
}

/** tries to upgrade superindicator constraint to an indicator constraint */
static
SCIP_RETCODE upgradeIndicatorSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< superindicator constraint to be upgraded */
   SCIP_Bool*            success,            /**< pointer to store if the upgrading was successful */
   SCIP_Bool*            deleted             /**< pointer to store if the constraint was deleted */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* indcons;

   SCIP_Real lhs;
   SCIP_Real rhs;
   char name[SCIP_MAXSTRLEN];
   int i;

#ifdef SCIP_DEBUG
   int nnewconss;
#endif

   assert(scip != NULL);
   assert(cons != NULL);
   assert(success != NULL);
   assert(deleted != NULL);

   *success = FALSE;
   *deleted = FALSE;

   SCIPdebug( nnewconss = 0 );

   /* get data of superindicator constraint */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* upgrade only for linear slack constraint */
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(consdata->slackcons)), "linear") != 0 )
      return SCIP_OKAY;

   /* upgrade only if indicator constraint handler found */
   conshdlr = SCIPfindConshdlr(scip, "indicator");
   if( conshdlr == NULL )
      return SCIP_OKAY;

   /* if linear slack constraint is free we can delete the superindicator constraint */
   lhs = SCIPgetLhsLinear(scip, consdata->slackcons);
   rhs = SCIPgetRhsLinear(scip, consdata->slackcons);
   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );
      *deleted = TRUE;

      SCIPdebugMsg(scip, "constraint <%s> deleted because of free slack constraint\n", SCIPconsGetName(cons));

      return SCIP_OKAY;
   }

   /* upgrade rhs inequality */
   if( !SCIPisInfinity(scip, rhs) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_upgd_indrhs", SCIPconsGetName(cons));

      SCIP_CALL( SCIPcreateConsIndicator(scip, &indcons, name, consdata->binvar, SCIPgetNVarsLinear(scip, consdata->slackcons),
            SCIPgetVarsLinear(scip, consdata->slackcons), SCIPgetValsLinear(scip, consdata->slackcons), rhs,
            SCIPconsIsInitial(cons),  SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );

      SCIP_CALL( SCIPaddCons(scip, indcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &indcons) );

      SCIPdebug( nnewconss++ );
   }

   /* upgrade lhs inequality */
   if( !SCIPisInfinity(scip, -lhs) )
   {
      SCIP_Real* negvals;
      SCIP_Real* vals;
      int nvars;

      vals = SCIPgetValsLinear(scip, consdata->slackcons);
      nvars = SCIPgetNVarsLinear(scip, consdata->slackcons);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_upgd_indlhs", SCIPconsGetName(cons));

      /* create array of negated coefficient values */
      SCIP_CALL( SCIPallocBufferArray(scip, &negvals, nvars) );
      for( i = nvars-1; i >= 0; i-- )
         negvals[i] = -vals[i];

      SCIP_CALL( SCIPcreateConsIndicator(scip, &indcons, name, consdata->binvar, nvars,
            SCIPgetVarsLinear(scip, consdata->slackcons), negvals, -lhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );

      SCIP_CALL( SCIPaddCons(scip, indcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &indcons) );

      SCIPfreeBufferArray(scip, &negvals);

      SCIPdebug( nnewconss++ );
   }

   SCIPdebug( SCIPdebugMsg(scip, "constraint <%s> upgraded to %d indicator constraint%s\n",
         SCIPconsGetName(cons), nnewconss, nnewconss == 1 ? "" : "s") );

   /* delete the superindicator constraint */
   SCIP_CALL( SCIPdelCons(scip, cons) );
   *success = TRUE;

   return SCIP_OKAY;
}

/** upgrades a superindicator constraint to a linear constraint if possible */
static
SCIP_RETCODE upgradeLinearSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< superindicator constraint to be upgraded */
   SCIP_Bool*            success,            /**< pointer to store if the upgrading was successful */
   SCIP_Bool*            deleted             /**< pointer to store if the constraint was deleted */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* slackcons;
   SCIP_VAR** slackvars;
   SCIP_VAR** newvars;
   SCIP_Real* slackvals;
   SCIP_Real* newvals;

   SCIP_Real maxcoef;
   SCIP_Real minabscoef;
   SCIP_Real minact;
   SCIP_Real maxact;
   SCIP_Real lhs;
   SCIP_Real rhs;

   int nvars;
   int i;

#ifdef SCIP_DEBUG
   int nnewconss;
#endif

   assert(scip != NULL);
   assert(cons != NULL);
   assert(success != NULL);
   assert(deleted != NULL);

   *success = FALSE;
   *deleted = FALSE;

   SCIPdebug( nnewconss = 0 );

   /* get data of superindicator constraint */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   slackcons = consdata->slackcons;
   assert(slackcons != NULL);

   /* upgrade only for linear slack constraint */
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "linear") != 0 )
      return SCIP_OKAY;

   /**@todo store in conshdlrdata */

   /* upgrade only if linear constraint handler found */
   conshdlr = SCIPfindConshdlr(scip, "linear");
   if( conshdlr == NULL )
      return SCIP_OKAY;

   /* if linear slack constraint is free we can delete the superindicator constraint */
   rhs = SCIPgetRhsLinear(scip, slackcons);
   lhs = SCIPgetLhsLinear(scip, slackcons);

   if( SCIPisInfinity(scip, rhs) && SCIPisInfinity(scip, -lhs) )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );
      *deleted = TRUE;

      SCIPdebugMsg(scip, "constraint <%s> deleted because of free slack constraint\n", SCIPconsGetName(cons));

      return SCIP_OKAY;
   }

   /* if linear slack constraint is redundant due to bounded activities we can delete the superindicator constraint */
   SCIP_CALL( extractLinearValues(scip, slackcons, &minact, &maxact, &minabscoef) );
   assert(!SCIPisInfinity(scip, minact));
   assert(!SCIPisInfinity(scip, -maxact));

   if( (SCIPisInfinity(scip, -lhs) || SCIPisLE(scip, lhs, minact)) && (SCIPisInfinity(scip, rhs) || SCIPisGE(scip, rhs, maxact)) )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );
      *deleted = TRUE;

      SCIPdebugMsg(scip, "constraint <%s> deleted because of redundant slack constraint\n", SCIPconsGetName(cons));

      return SCIP_OKAY;
   }

   /* if the big-M coefficient is too large compared to the coefficients of the slack constraint, we do not upgrade to
    * avoid numerical problems
    */
   maxcoef = minabscoef * SCIPconshdlrGetData(SCIPconsGetHdlr(cons))->maxupgdcoeflinear;

   if( (!SCIPisInfinity(scip, rhs) && (SCIPisInfinity(scip, maxact) || SCIPisInfinity(scip, maxact - rhs) ||
            maxact - rhs > maxcoef)) ||
      (!SCIPisInfinity(scip, -lhs) && (SCIPisInfinity(scip, -minact) || SCIPisInfinity(scip, lhs - minact) ||
         lhs - minact > maxcoef)) )
   {
      SCIPdebugMsg(scip, "constraint <%s> not upgraded to a linear constraint due to large big-M coefficient\n",
         SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   /* allocating memory for new constraint */
   nvars = SCIPgetNVarsLinear(scip, slackcons);
   SCIP_CALL( SCIPallocBufferArray(scip, &newvars, nvars+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newvals, nvars+1) );

   /* copy the vars and the vals array */
   slackvars = SCIPgetVarsLinear(scip, slackcons);
   slackvals = SCIPgetValsLinear(scip, slackcons);

   assert(slackvars != NULL);
   assert(slackvals != NULL);

   for( i = nvars-1; i >= 0; i-- )
   {
      newvars[i] = slackvars[i];
      newvals[i] = slackvals[i];
   }

   /* add binary variable */
   newvars[nvars] = consdata->binvar;
   assert(newvars[nvars] != NULL);

   assert(!SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs));

   /* create the upgraded constraint for rhs inequality */
   if( !SCIPisInfinity(scip, rhs) )
   {
      SCIP_CONS* newcons;
      char name[SCIP_MAXSTRLEN];

      assert(!SCIPisInfinity(scip, -maxact) );
      assert(!SCIPisInfinity(scip, maxact));

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_upgd_linrhs", SCIPconsGetName(cons));

      /* compute big-M */
      newvals[nvars] = maxact - rhs;
      assert(!SCIPisInfinity(scip, newvals[nvars]));
      assert(!SCIPisInfinity(scip, -newvals[nvars]));

      /* rhs inequality is redundant if maxact is less equal rhs */
      if( SCIPisPositive(scip, newvals[nvars]) )
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, nvars+1, newvars, newvals, -SCIPinfinity(scip), maxact,
               SCIPconsIsInitial(cons),  SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

         SCIPdebug( nnewconss++ );
      }
   }

   /* create the upgraded constraint for rhs inequality */
   if( !SCIPisInfinity(scip, -lhs) )
   {
      SCIP_CONS* newcons;
      char name[SCIP_MAXSTRLEN];

      assert(!SCIPisInfinity(scip, minact));
      assert(!SCIPisInfinity(scip, -minact));

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_upgd_linlhs", SCIPconsGetName(cons));

      /* compute big-M */
      newvals[nvars] = minact - lhs;
      assert(!SCIPisInfinity(scip, newvals[nvars]));
      assert(!SCIPisInfinity(scip, -newvals[nvars]));

      /* lhs inequality is redundant if minact is greater equal lhs */
      if( SCIPisNegative(scip, newvals[nvars]) )
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, nvars+1, newvars, newvals, minact, SCIPinfinity(scip),
               SCIPconsIsInitial(cons),  SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
               SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
               SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

         SCIPdebug( nnewconss++ );
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &newvals);
   SCIPfreeBufferArray(scip, &newvars);

   SCIPdebug( SCIPdebugMsg(scip, "constraint <%s> upgraded to %d indicator constraint%s\n",
         SCIPconsGetName(cons), nnewconss, nnewconss == 1 ? "" : "s") );

   /* delete the superindicator constraint */
   SCIP_CALL( SCIPdelCons(scip, cons) );
   *success = TRUE;

   return SCIP_OKAY;
}

/** tries to upgrade a superindicator constraint in order of the upgrade priority parameters */
static
SCIP_RETCODE upgradeSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< superindicator constraint to be updated */
   SCIP_Bool*            success,            /**< pointer to store if the constraint was upgraded */
   SCIP_Bool*            deleted             /**< pointer to store if the constraint was deleted */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(success != NULL);
   assert(deleted != NULL);

   *success = FALSE;
   *deleted = FALSE;

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));

   /* indicator upgrade before linear upgrade */
   if( conshdlrdata->upgdprioindicator > conshdlrdata->upgdpriolinear )
   {
      assert(conshdlrdata->upgdprioindicator >= 0);

      SCIP_CALL( upgradeIndicatorSuperindicator(scip, cons, success, deleted) );

      if( !*deleted && !*success && conshdlrdata->upgdpriolinear >= 0 )
      {
         SCIP_CALL( upgradeLinearSuperindicator(scip, cons, success, deleted) );
      }
   }
   /* linear upgrade before indicator upgrade */
   else if( conshdlrdata->upgdpriolinear >= 0 )
   {
      SCIP_CALL( upgradeLinearSuperindicator(scip, cons, success, deleted) );

      if( !*deleted && !*success && conshdlrdata->upgdprioindicator >= 0 )
      {
         SCIP_CALL( upgradeIndicatorSuperindicator(scip, cons, success, deleted) );
      }
   }

   return SCIP_OKAY;
}

/** helper function to enforce constraints */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{  /*lint --e{715}*/
   SCIP_Bool cont;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   /* if the solution is infeasible anyway, skip the enforcement */
   if( solinfeasible )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "executing enforcement callback for %s solution\n", sol == NULL ? "LP" : "relaxation");

   cont = TRUE;
   *result = SCIP_FEASIBLE;

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif

   /* check all constraints */
   for( i = nconss-1; i >= 0 && cont; i-- )
   {
      SCIP_CONSDATA* consdata;
      SCIP_RESULT locresult;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      locresult = SCIP_FEASIBLE;

      /* enforce only if binvar is fixed to one */
      if( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->binvar), 1.0));

         if( sol == NULL )
         {
            SCIPdebugMsg(scip, "binvar <%s> == 1 locally --> SCIPenfolpCons() on constraint <%s>\n",
               SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

            SCIP_CALL( SCIPenfolpCons(scip, consdata->slackcons, solinfeasible, &locresult) );
         }
         else
         {
            SCIPdebugMsg(scip, "binvar <%s> == 1 locally --> SCIPenforelaxCons() on constraint <%s>\n",
               SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

            SCIP_CALL( SCIPenforelaxCons(scip, consdata->slackcons, sol, solinfeasible, &locresult) );
         }

         SCIPdebugPrintf(" --> %slocresult=%d\n", locresult == SCIP_FEASIBLE ? "satisfied, " : "", locresult);
      }
      /* otherwise check if we have not yet detected infeasibility */
      else if( *result == SCIP_FEASIBLE )
      {
         SCIP_CALL( consdataCheckSuperindicator(scip, consdata, sol, TRUE, FALSE, FALSE, &locresult) );
      }

      /* evaluate result */
      switch( locresult )
      {
      case SCIP_CUTOFF:
      case SCIP_BRANCHED:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         *result = locresult;
         cont = FALSE;
         break;
      case SCIP_CONSADDED:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF )
            *result = locresult;
         break;
      case SCIP_REDUCEDDOM:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED )
            *result = locresult;
         break;
      case SCIP_SEPARATED:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM )
            *result = locresult;
         break;
      case SCIP_INFEASIBLE:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_BRANCHED )
            *result = locresult;
         break;
      case SCIP_FEASIBLE:
         break;
      default:
         SCIPerrorMessage("invalid SCIP result %d\n", locresult);
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/
   }

   SCIPdebugMsg(scip, "enforcement result=%d\n", *result);

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySuperindicator)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSuperindicator(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   SCIPdebugMsg(scip, "freeing superindicator constraint handler data\n");

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   SCIPdebugMsg(scip, "initializing presolving\n");

   for( i = nconss-1; i >= 0; i-- )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* make the constraint local to avoid wrong propagation */
      SCIP_CALL( SCIPsetConsLocal(scip, consdata->slackcons, TRUE) );
   }

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSuperindicator)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->slackcons != NULL);

   SCIPdebugMsg(scip, "deleting constraint <%s>\n", SCIPconsGetName(cons));

   /* we have to release the slack constraint also in case we transformed it manually since it is captured automatically
    * in SCIPtransformCons()
    */
   SCIP_CALL( SCIPreleaseCons(scip, &((*consdata)->slackcons)) );

   /* free memory */
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   char newname[SCIP_MAXSTRLEN];

   SCIPdebugMsg(scip, "transforming superindicator constraint <%s>\n", SCIPconsGetName(sourcecons));

   /* get constraint data of source constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   (void) SCIPsnprintf(newname, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons) );
   SCIP_CALL( consdataCreateSuperindicator(scip, &targetdata, sourcedata->binvar, sourcedata->slackcons) );

   /* create target constraint and capture it at the same time */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, newname, conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpSuperindicator)
{
   int c;

   assert(scip != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   SCIPdebugMsg(scip, "executing initlp callback\n");

   for( c = nconss-1; c >= 0 && !(*infeasible); c-- )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);

      assert(consdata != NULL);
      assert(SCIPconsIsInitial(conss[c]));

      if( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->binvar), 1.0));

         SCIPdebugMsg(scip, "binvar <%s> == 1 --> SCIPinitlpCons() on constraint <%s>\n",
            SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

         SCIP_CALL( SCIPinitlpCons(scip, consdata->slackcons, infeasible) );
      }
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSuperindicator)
{  /*lint --e{715}*/
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;

   SCIPdebugMsg(scip, "executing sepalp callback\n");

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPprintSol(scip, NULL, NULL, FALSE) );
#endif

   /* check all useful constraints */
   for( c = nusefulconss-1; c >= 0 && *result != SCIP_CUTOFF; c-- )
   {
      SCIP_CONSDATA* consdata;
      SCIP_RESULT locresult;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      locresult = SCIP_DELAYED;

      /* separate only if binvar is fixed to one */
      if( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->binvar), 1.0));

         SCIPdebugMsg(scip, "binvar <%s> == 1 --> SCIPsepalpCons() on constraint <%s>\n",
            SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

         SCIP_CALL( SCIPsepalpCons(scip, consdata->slackcons, &locresult) );

         SCIPdebugMsgPrint(scip, " --> locresult=%d\n", locresult);
      }

      /* evaluate result value */
      switch( locresult )
      {
      case SCIP_CUTOFF:
      case SCIP_CONSADDED:
         assert(*result != SCIP_CUTOFF);
         *result = locresult;
         break;
      case SCIP_REDUCEDDOM:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED )
            *result = locresult;
         break;
      case SCIP_SEPARATED:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM )
            *result = locresult;
         break;
      case SCIP_NEWROUND:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED )
            *result = locresult;
         break;
      case SCIP_DIDNOTFIND:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_NEWROUND
            && *result != SCIP_SEPARATED )
            *result = locresult;
         break;
      case SCIP_DIDNOTRUN:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_NEWROUND
            && *result != SCIP_SEPARATED
            && *result != SCIP_DIDNOTFIND )
            *result = locresult;
         break;
      case SCIP_INFEASIBLE:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_NEWROUND )
            *result = locresult;
         break;
      case SCIP_DELAYED:
         break;
      default:
         SCIPerrorMessage("invalid SCIP result %d\n", locresult);
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/
   }

   SCIPdebugMsg(scip, "sepalp result=%d\n", *result);

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolSuperindicator)
{  /*lint --e{715}*/
   int c;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;

   SCIPdebugMsg(scip, "executing sepasol callback\n");

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPprintSol(scip, NULL, NULL, FALSE) );
#endif


   /* check all the useful constraint */
   for( c = 0; c < nusefulconss && *result != SCIP_CUTOFF; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_RESULT locresult;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      locresult = SCIP_DELAYED;

      /* separate only if binvar is fixed to one */
      if( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->binvar), 1.0));

         SCIPdebugMsg(scip, "binvar <%s> == 0 --> SCIPsepasolCons() on constraint <%s>\n",
            SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

         SCIP_CALL( SCIPsepasolCons(scip, consdata->slackcons, sol, &locresult) );

         SCIPdebugMsgPrint(scip, " --> result=%d\n", locresult);
      }

      /* evaluate result value */
      switch( locresult )
      {
      case SCIP_CUTOFF:
      case SCIP_CONSADDED:
         assert(*result != SCIP_CUTOFF);
         *result = locresult;
         break;
      case SCIP_REDUCEDDOM:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED )
            *result = locresult;
         break;
      case SCIP_SEPARATED:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM )
            *result = locresult;
         break;
      case SCIP_NEWROUND:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED )
            *result = locresult;
         break;
      case SCIP_DIDNOTFIND:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_NEWROUND
            && *result != SCIP_SEPARATED )
            *result = locresult;
         break;
      case SCIP_DIDNOTRUN:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_NEWROUND
            && *result != SCIP_SEPARATED
            && *result != SCIP_DIDNOTFIND )
            *result = locresult;
         break;
      case SCIP_INFEASIBLE:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_NEWROUND )
            *result = locresult;
         break;
      case SCIP_DELAYED:
         break;
      default:
         SCIPerrorMessage("invalid SCIP result %d\n", locresult);
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/
   }

   SCIPdebugMsg(scip, "sepa sol result=%d\n", *result);

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSuperindicator)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, solinfeasible, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSuperindicator)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, solinfeasible, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSuperindicator)
{  /*lint --e{715}*/
   SCIP_Bool cont;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   /* if the solution is infeasible anyway, skip the enforcement */
   if( solinfeasible )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }
   else if( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "executing enfops callback\n");

   *result = SCIP_FEASIBLE;
   cont = TRUE;

   /* check all contraints */
   for( i = nconss-1; i >= 0 && cont; i-- )
   {
      SCIP_CONSDATA* consdata;
      SCIP_RESULT locresult;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      locresult = SCIP_DIDNOTRUN;

      /* enforce only if binvar is fixed to one */
      if( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->binvar), 1.0));

         SCIPdebugMsg(scip, "binvar <%s> == 1 locally --> SCIPenfopsCons() on constraint <%s>\n",
            SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

         SCIP_CALL( SCIPenfopsCons(scip, consdata->slackcons, solinfeasible, objinfeasible, &locresult) );

         SCIPdebugMsgPrint(scip, " --> %slocresult=%d\n", locresult == SCIP_FEASIBLE ? "satisfied, " : "", locresult);
      }
      /* otherwise check if we have not yet detected infeasibility */
      else if( *result == SCIP_FEASIBLE || *result == SCIP_DIDNOTRUN )
      {
         SCIP_CALL( consdataCheckSuperindicator(scip, consdata, NULL, TRUE, FALSE, FALSE, &locresult) );

      }

      /* evaluate result value */
      switch( locresult )
      {
      case SCIP_CUTOFF:
      case SCIP_BRANCHED:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         *result = locresult;
         cont = FALSE;
         break;
      case SCIP_CONSADDED:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF )
            *result = locresult;
         break;
      case SCIP_REDUCEDDOM:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED )
            *result = locresult;
         break;
      case SCIP_SOLVELP:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED )
            *result = locresult;
         break;
      case SCIP_INFEASIBLE:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED
            && *result != SCIP_SOLVELP )
            *result = locresult;
         break;
      case SCIP_DIDNOTRUN:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED
            && *result != SCIP_SOLVELP
            && *result != SCIP_INFEASIBLE )
            *result = locresult;
         break;
      case SCIP_FEASIBLE:
         assert(*result != SCIP_CUTOFF);
         assert(*result != SCIP_BRANCHED);
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED
            && *result != SCIP_SOLVELP
            && *result != SCIP_INFEASIBLE
            && *result != SCIP_DIDNOTRUN )
            *result = locresult;
         break;
      default:
         SCIPerrorMessage("invalid SCIP result %d\n", locresult);
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/
   }

   SCIPdebugMsg(scip, "enfops result=%d\n", *result);

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSuperindicator)
{  /*lint --e{715}*/
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);
   assert(sol != NULL);

   *result = SCIP_FEASIBLE;

   for( i = nconss-1; i >= 0 && (*result == SCIP_FEASIBLE || completely); i-- )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[i]);
      SCIP_CALL( consdataCheckSuperindicator(scip, consdata, sol, checkintegrality, checklprows, printreason, result) );
   }

   SCIPdebugMsg(scip, "checked solution from <%s> (checkintegrality=%u, checklprows=%u) --> result=%d (%sfeasible)\n",
      SCIPsolGetHeur(sol) == NULL ? "NULL" : SCIPheurGetName(SCIPsolGetHeur(sol)), checkintegrality, checklprows,
      *result, *result == SCIP_INFEASIBLE ? "in" : "");

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSuperindicator)
{  /*lint --e{715}*/
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "executing prop callback\n");

   /* loop over all useful contraints */
   for( i = nusefulconss-1; i >= 0 && *result != SCIP_CUTOFF; i-- )
   {
      SCIP_CONSDATA* consdata;
      SCIP_RESULT locresult;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      locresult = SCIP_DIDNOTRUN;

      /* propagate only if binvar is fixed to one */
      if( SCIPvarGetLbGlobal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->binvar), 1.0));

         SCIPdebugMsg(scip, "binvar <%s> == 1 globally --> deleting superindicator and adding slack constraint <%s>\n",
            SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

         SCIP_CALL( SCIPsetConsLocal(scip, consdata->slackcons, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, consdata->slackcons) );
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );

         locresult = SCIP_DIDNOTFIND;
      }
      else if( SCIPvarGetLbLocal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->binvar), 1.0));

         SCIPdebugMsg(scip, "binvar <%s> == 1 locally --> propagating slack constraint <%s>\n",
            SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

         SCIP_CALL( SCIPpropCons(scip, consdata->slackcons, proptiming, &locresult) );

         SCIPdebugMsgPrint(scip, " --> locresult=%d\n", locresult);
      }
      /**@todo else propagate the domain of the binvar as well: start probing mode, fix binvar to one, propagate
       *       constraint, and see whether we become infeasible; if this is implemented, the resprop callback must be
       *       updated
       */

      /* evaluate result value */
      switch( locresult )
      {
      case SCIP_CUTOFF:
      case SCIP_DELAYED:
         /* if propagation of one constraint is delayed, we want to propagate again unless the node is cut off */
         assert(*result != SCIP_CUTOFF);
         *result = locresult;
         break;
      case SCIP_REDUCEDDOM:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_DELAYED )
            *result = locresult;
         break;
      case SCIP_DIDNOTFIND:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_REDUCEDDOM
            && *result != SCIP_DELAYED )
            *result = locresult;
         break;
      case SCIP_DIDNOTRUN:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_REDUCEDDOM
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DELAYED )
            *result = locresult;
         break;
      default:
         SCIPerrorMessage("invalid SCIP result %d\n", locresult);
         return SCIP_INVALIDRESULT;
      }  /*lint !e788*/
   }

   SCIPdebugMsg(scip, "prop result=%d\n", *result);

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSuperindicator)
{  /*lint --e{715}*/
   int i;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(conshdlr != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "executing presol callback\n");

   for( i = nconss-1; i >= 0 && *result != SCIP_CUTOFF; i-- )
   {
      SCIP_CONSDATA* consdata;
      SCIP_RESULT locresult;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      locresult = SCIP_DIDNOTFIND;

      /**@todo check whether the slack constraint is added to SCIP; in this case the superindicator can be deleted */

      /**@todo check whether the slack constraint is a superindicator constraint and presolve */

      /* if binvar is globally fixed to 1, we add the slack constraint and remove the superindicator */
      if( SCIPvarGetLbGlobal(consdata->binvar) > 0.5 )
      {
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->binvar), 1.0));

         SCIPdebugMsg(scip, "binvar <%s> == 1 globally --> deleting superindicator and adding slack constraint <%s>\n",
            SCIPvarGetName(consdata->binvar), SCIPconsGetName(consdata->slackcons));

         SCIP_CALL( SCIPsetConsLocal(scip, consdata->slackcons, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, consdata->slackcons) );
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );

         locresult = SCIP_SUCCESS;
      }
      /* otherwise try upgrading */
      else
      {
         SCIP_Bool success;
         SCIP_Bool deleted;

         SCIP_CALL( upgradeSuperindicator(scip, conss[i], &success, &deleted) );

         /* update statistics */
         if( deleted )
            (*ndelconss)++;
         else if( success )
            (*nupgdconss)++;

         /**@todo mark if upgrading failed to avoid trying too often; however, since upgrading might fail only due to
          *       large domains, we may want to try again later, e.g., if SCIPisPresolveFinished() is TRUE
          */

         if( deleted || success )
            locresult = SCIP_SUCCESS;
      }
      /**@todo else propagate the domain of the binvar as well: start probing mode, fix binvar to one, propagate
       *       constraint, and see whether we become infeasible
       */

      /* evaluate result value */
      switch( locresult )
      {
      case SCIP_SUCCESS:
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_DELAYED )
            *result = locresult;
         break;
      default:
         assert(locresult == SCIP_DIDNOTFIND);
         assert(*result != SCIP_CUTOFF);
         if( *result != SCIP_UNBOUNDED && *result != SCIP_DELAYED && *result != SCIP_SUCCESS )
            *result = locresult;
         break;
      } /*lint !e788*/
   }

   SCIPdebugMsg(scip, "presol result=%d\n", *result);

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(infervar != NULL);
   assert(bdchgidx != NULL);
   assert(result != NULL);

   SCIPdebugMsg(scip, "executing resprop callback for constraint <%s>\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTFIND;

   /* check that we only propagated if the binvar is fixed to one */
   assert(SCIPisFeasEQ(scip, SCIPgetVarUbAtIndex(scip, consdata->binvar, bdchgidx, TRUE), 1.0));

   /* add tightened lower bound on binvar to conflict set */
   SCIP_CALL( SCIPaddConflictLb(scip, consdata->binvar, bdchgidx) );

   /* call propagation conflict resolving method for the slack constraint */
   SCIP_CALL( SCIPrespropCons(scip, consdata->slackcons, infervar, inferinfo, boundtype, bdchgidx, relaxedbd, result) );

   SCIPdebugMsgPrint(scip, " --> result=%d\n", *result);

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);

   SCIPdebugMsg(scip, "locking variables for constraint <%s>\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* lock binvar up */
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->binvar, nlocksneg, nlockspos) );

   /* call lock method for the slack constraint */
   SCIP_CALL( SCIPaddConsLocks(scip, consdata->slackcons, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR* binvar;
   int zeroone;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get binary variable */
   binvar = consdata->binvar;
   assert(binvar != NULL);

   /* resolve negation if necessary */
   zeroone = 1;
   if ( SCIPvarGetStatus(binvar) == SCIP_VARSTATUS_NEGATED )
   {
      zeroone = 0;
      binvar = SCIPvarGetNegatedVar(binvar);
      assert(binvar != NULL);
   }

   /* print name of the binary variable */
   SCIP_CALL( SCIPwriteVarName(scip, file, binvar, TRUE) );

   /* print implication */
   SCIPinfoMessage(scip, file, " = %d ->", zeroone);

   /* print slack constraint */
   assert(consdata->slackcons != NULL);
   SCIP_CALL( SCIPprintCons(scip, consdata->slackcons, file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlrslack;
   SCIP_CONSDATA* sourceconsdata;
   SCIP_CONS* sourceslackcons;
   SCIP_CONS* targetslackcons;
   SCIP_VAR* targetbinvar;
   const char* consname;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0);

   *valid = TRUE;

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIPdebugMsg(scip, "copying superindicator constraint <%s> to <%s>\n", SCIPconsGetName(sourcecons), consname);

   if( modifiable )
   {
      SCIPwarningMessage(scip, "cannot create modifiable superindicator constraint when trying to copy constraint <%s>\n",
         SCIPconsGetName(sourcecons));
      *valid = FALSE;
      return SCIP_OKAY;
   }

   sourceconsdata = SCIPconsGetData(sourcecons);
   assert(sourceconsdata != NULL);

   /* get slack constraint */
   sourceslackcons = sourceconsdata->slackcons;
   assert(sourceslackcons != NULL);

   /* if the slack constraint has been deleted, create an empty linear constraint */
   if( SCIPconsIsDeleted(sourceslackcons) )
   {
      SCIPdebugMsg(scip, "slack constraint <%s> deleted; creating empty linear constraint\n",
         SCIPconsGetName(sourceslackcons));

      SCIP_CALL( SCIPcreateConsLinear(scip, &targetslackcons, "dummy", 0, NULL, NULL, 0.0, SCIPinfinity(scip),
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, targetslackcons) );
   }
   else
   {
      /* get copied version of slack constraint */
      conshdlrslack = SCIPconsGetHdlr(sourceslackcons);
      assert(conshdlrslack != NULL);

      /* if copying scip after transforming the original instance before presolving, we need to correct the slack
       * constraint pointer
       */
      assert(!SCIPisTransformed(sourcescip) || SCIPconsIsTransformed(sourceslackcons));
      if( SCIPisTransformed(sourcescip) && !SCIPconsIsTransformed(sourceslackcons) )
      {
	 SCIP_CONS* transslackcons;

         SCIP_CALL( SCIPgetTransformedCons(sourcescip, sourceslackcons, &transslackcons) );
         assert(transslackcons != NULL);
         SCIP_CALL( SCIPreleaseCons(sourcescip, &sourceconsdata->slackcons) );
         SCIP_CALL( SCIPcaptureCons(sourcescip, transslackcons) );

         sourceconsdata->slackcons = transslackcons;
         sourceslackcons = transslackcons;
      }

      SCIP_CALL( SCIPgetConsCopy(sourcescip, scip, sourceslackcons, &targetslackcons, conshdlrslack, varmap, consmap,
            SCIPconsGetName(sourceslackcons), SCIPconsIsInitial(sourceslackcons), SCIPconsIsSeparated(sourceslackcons),
            SCIPconsIsEnforced(sourceslackcons), SCIPconsIsChecked(sourceslackcons), SCIPconsIsPropagated(sourceslackcons),
            SCIPconsIsLocal(sourceslackcons), SCIPconsIsModifiable(sourceslackcons), SCIPconsIsDynamic(sourceslackcons),
            SCIPconsIsRemovable(sourceslackcons), SCIPconsIsStickingAtNode(sourceslackcons), global, valid) );
   }

   /* find copied variable corresponding to binvar */
   if( *valid )
   {
      SCIP_VAR* sourcebinvar;

      sourcebinvar = sourceconsdata->binvar;
      assert(sourcebinvar != NULL);

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcebinvar, &targetbinvar, varmap, consmap, global, valid) );
   }
   else
      targetbinvar = NULL;

   /* create superindicator constraint */
   if( *valid )
   {
      assert(targetslackcons != NULL);
      assert(targetbinvar != NULL);
      assert(!modifiable);

      SCIP_CALL( SCIPcreateConsSuperindicator(scip, cons, consname, targetbinvar, targetslackcons,
            initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );
   }

   /* relase slack constraint */
   if( targetslackcons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &targetslackcons) );
   }

   if( !(*valid) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "could not copy superindicator constraint <%s>\n", SCIPconsGetName(sourcecons));
   }

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSuperindicator)
{  /*lint --e{715}*/
   SCIP_VAR* binvar;
   SCIP_CONS* slackcons;
   char binvarname[1024];
   const char* slackstr;
   int zeroone;
   int nargs;

   assert(cons != NULL);
   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);

   *success = FALSE;

   /* extract binary variable name and value which triggers slack constraint */
   nargs = sscanf(str, " <%1023[^>]>[B] = %d", binvarname, &zeroone);

   if( nargs != 2 || (zeroone != 0 && zeroone != 1) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected the following form: <var> = [0|1] ->  <cons>\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "got: %s\n", str);
      return SCIP_OKAY;
   }

   /* extract string describing slack constraint */
   slackstr = strstr(str, "->");

   if( slackstr == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected the following form: <var> = [0|1] ->  <cons>\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "got: %s\n", str);
      return SCIP_OKAY;
   }

   slackstr = strstr(slackstr, "[");

   if( slackstr == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: expected the following form: <var> = [0|1] ->  <cons>\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "got: %s\n", str);
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "binvarname=%s, zeroone=%d, slackstr=%s\n", binvarname, zeroone, slackstr);

   /* get binary variable */
   binvar = SCIPfindVar(scip, binvarname);
   if( binvar == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable <%s>\n", binvarname);
      return SCIP_OKAY;
   }

   /* resolve negation if necessary */
   if( zeroone == 0 )
   {
      SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &binvar) );
   }

   /**@todo get slack constraint name and check whether constraint already exists; however, using only SCIPfindCons() is
    *       not sufficient since slack constraints are not added to the problem; do we need something like
    *       SCIPfindConsInConshdlr()?; currently, if there are two superindicator constraints with same slack constraint
    *       (binvars may be different), then after writing and reading, the slack constraint will be created twice with
    *       identical constraint name; this is not incorrect, but might consume more memory or time
    */

   /* parse slack constraint string */
   SCIP_CALL( SCIPparseCons(scip, &slackcons, slackstr, initial, separate, enforce, check, propagate, local, modifiable,
         dynamic, removable, stickingatnode, success) );

   if( *success )
   {
      assert(binvar != NULL);
      assert(slackcons != NULL);

      /* create the superindicator constraint */
      SCIP_CALL( SCIPcreateConsSuperindicator(scip, cons, name, binvar, slackcons,
            initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );

      /* the new superindicator constraint captured the slack constraint, so we can release it now */
      SCIP_CALL( SCIPreleaseCons(scip, &slackcons) );
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* must be ready to hold at least the binary variable */
   if( varssize <= 0 )
      *success = FALSE;
   else
   {
      /* add binary variable */
      vars[0] = consdata->binvar;

      /* add variables of slack constraint */
      SCIP_CALL( SCIPgetConsVars(scip, consdata->slackcons, &(vars[1]), varssize-1, success) );
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSuperindicator)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get number of variables in slack constraint */
   SCIP_CALL( SCIPgetConsNVars(scip, consdata->slackcons, nvars, success) );

   /* add binary variable */
   if( *success )
      (*nvars)++;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for superindicator constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSuperindicator(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_DIALOG* root;
   SCIP_DIALOG* changemenu;
   SCIP_DIALOG* dialog;

   /* create superindicator constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   conshdlrdata->nrejects = 0;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSuperindicator, consEnfopsSuperindicator, consCheckSuperindicator, consLockSuperindicator,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySuperindicator, consCopySuperindicator) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSuperindicator, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSuperindicator, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSuperindicator, consSepasolSuperindicator, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSuperindicator) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSuperindicator) );

   /* includes or updates the default dialog menus in SCIP */
   SCIP_CALL( SCIPincludeDialogDefault(scip) );

   root = SCIPgetRootDialog(scip);
   assert(root != NULL);

   /* find change menu */
   if( !SCIPdialogHasEntry(root, "change") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &changemenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "change", "change the problem", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, changemenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &changemenu) );
   }

   if( SCIPdialogFindEntry(root, "change", &changemenu) != 1 )
   {
      SCIPerrorMessage("change sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* add minuc dialog */
   if( !SCIPdialogHasEntry(changemenu, "minuc") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecChangeMinUC, NULL, NULL,
            "minuc", "transforms the current problem into a MinUC problem minimizing the number of unsatisfied constraints",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, changemenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* add constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/checkslacktype",
         "should type of slack constraint be checked when creating superindicator constraint?",
         &conshdlrdata->checkslacktype, TRUE, DEFAULT_CHECKSLACKTYPE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/maxupgdcoeflinear",
         "maximum big-M coefficient of binary variable in upgrade to a linear constraint (relative to smallest coefficient)",
         &conshdlrdata->maxupgdcoeflinear, TRUE, DEFAULT_MAXUPGDCOEFLINEAR, 0.0, 1e15, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/upgdprioindicator",
         "priority for upgrading to an indicator constraint (-1: never)",
         &conshdlrdata->upgdprioindicator, TRUE, DEFAULT_UPGDPRIOINDICATOR, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/upgdpriolinear",
         "priority for upgrading to an indicator constraint (-1: never)",
         &conshdlrdata->upgdpriolinear, TRUE, DEFAULT_UPGDPRIOLINEAR, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a superindicator constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< pointer to the indicator constraint  */
   SCIP_CONS*            slackcons,          /**< constraint corresponding to the handled constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_Bool modifiable;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(name != NULL);
   assert(binvar != NULL);
   assert(slackcons != NULL);

   modifiable = FALSE;

   /* find the superindicator constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("superindicator constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* only allow types of slack constraints that can be handled */
   if( conshdlrdata->checkslacktype &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "abspower") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "and") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "bivariate") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "bounddisjunction") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "conjunction") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "disjunction") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "knapsack") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "linear") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "linking") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "logicor") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "nonlinear") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "or") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "quadratic") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "soc") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "SOS1") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "SOS2") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "cumulative") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "varbound") != 0 &&
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)), "superindicator") != 0
      )
   {
      if( conshdlrdata->nrejects < 5 )
      {
         SCIPwarningMessage(scip, "rejected creation of superindicator with slack constraint <%s> of type <%s> "
            "(use parameter <checkslacktype> to disable check)\n",
            SCIPconsGetName(slackcons), SCIPconshdlrGetName(SCIPconsGetHdlr(slackcons)));
         conshdlrdata->nrejects++;
      }

      if( conshdlrdata->nrejects == 5 )
      {
         SCIPwarningMessage(scip, "suppressing further warning messages of this type\n");
         conshdlrdata->nrejects++;
      }

      return SCIP_INVALIDCALL;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreateSuperindicator(scip, &consdata, binvar, slackcons) );
   assert(consdata != NULL);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a superindicator constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsSuperindicator(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsSuperindicator() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSuperindicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             binvar,             /**< pointer to the indicator constraint  */
   SCIP_CONS*            slackcons           /**< constraint corresponding to the handled constraint */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(name != NULL);
   assert(binvar != NULL);
   assert(slackcons != NULL);

   SCIP_CALL( SCIPcreateConsSuperindicator(scip, cons, name, binvar, slackcons,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** gets binary variable corresponding to the general indicator constraint */
SCIP_VAR* SCIPgetBinaryVarSuperindicator(
   SCIP_CONS*            cons                /**< superindicator constraint */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->binvar;
}

/** gets the slack constraint corresponding to the general indicator constraint */
SCIP_CONS* SCIPgetSlackConsSuperindicator(
   SCIP_CONS*            cons                /**< superindicator constraint */
   )
{
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->slackcons;
}


/*
 *  constraint-dependent SCIP methods
 */

/** transforms the current problem into a MinUC problem (minimizing the number of unsatisfied constraints),
 *  a CIP generalization of the MinULR (min. unsatisfied linear relations) problem
 */
SCIP_RETCODE SCIPtransformMinUC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            success             /**< pointer to store whether all constraints could be transformed */
   )
{
   SCIP_CONS** conss;
   SCIP_CONS** probconss;
   SCIP_VAR** vars;
   char consname[SCIP_MAXSTRLEN];
   char varname[SCIP_MAXSTRLEN];
   int maxbranchprio;
   int ntransconss;
   int nconss;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( SCIPgetStage(scip) !=  SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("method <SCIPtransformMinUC> can only be called in problem stage\n");
      return SCIP_INVALIDCALL;
   }

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* copy the conss array because it changes when adding and deleting constraints */
   nconss = SCIPgetNConss(scip);
   probconss = SCIPgetConss(scip);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &conss, probconss, nconss) );

   /* clear objective function and compute maximal branching priority */
   maxbranchprio = 0;
   for( i = nvars-1; i >= 0; i-- )
   {
      SCIP_CALL( SCIPchgVarObj(scip, vars[i], 0.0) );

      if( SCIPvarGetBranchPriority(vars[i]) > maxbranchprio )
         maxbranchprio = SCIPvarGetBranchPriority(vars[i]);
   }

   maxbranchprio++;

   /* transform each constraint to slack constraint in a newly created superindicator constraint; note that we also need
    * to transform superindicator constraints, since their binary variable might have down-locks
    */
   ntransconss = 0;
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONS* cons;
      SCIP_CONS* supindcons;
      SCIP_VAR* binvar;
      SCIP_VAR* negbinvar;
      SCIP_RETCODE retcode;

      cons = conss[i];
      assert(cons != NULL);

      /* create a new binary variable with objective coefficient one */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_master", SCIPconsGetName(cons));

      SCIP_CALL( SCIPcreateVar(scip, &binvar, varname, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY,
            TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

      /* get negated variable, since we want to minimize the number of violated constraints */
      SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &negbinvar) );

      /* create superindicator constraint */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_super", SCIPconsGetName(cons));

      retcode = SCIPcreateConsSuperindicator(scip, &supindcons, consname, negbinvar, cons,
         SCIPconsIsInitial(cons),  SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons));

      if( retcode == SCIP_OKAY )
      {
         /* add binary variable and increase its branching priority */
         SCIP_CALL( SCIPaddVar(scip, binvar) );
         SCIP_CALL( SCIPchgVarBranchPriority(scip, binvar, maxbranchprio) );

         /* add superindicator constraint */
         SCIP_CALL( SCIPaddCons(scip, supindcons) );

         /* release binary variable and superindicator constraint */
         SCIP_CALL( SCIPreleaseVar(scip, &binvar) );
         SCIP_CALL( SCIPreleaseCons(scip, &supindcons) );

         /* delete slack constraint; it is still captured by the superindicator constraint */
         SCIP_CALL( SCIPdelCons(scip, cons) );

         ntransconss++;
      }
      else if( retcode == SCIP_INVALIDCALL )
      {
         SCIPdebugMsg(scip, "constraint <%s> of type <%s> could not be transformed to superindicator and was removed\n",
            SCIPconsGetName(cons), SCIPconshdlrGetName(SCIPconsGetHdlr(cons)));

         /* release binary variable */
         SCIP_CALL( SCIPreleaseVar(scip, &binvar) );

         /* delete slack constraint; this is necessary, because, e.g., the indicator expects its linear slack constraint
          * present in the problem, but this has just be transformed; hence, it cannot function any more and we have to
          * remove it
          */
         SCIP_CALL( SCIPdelCons(scip, cons) );
      }
      else
      {
         /* return all other error codes */
         SCIP_CALL( retcode );
      }
   }

   if( ntransconss == nconss )
      *success = TRUE;

   /* minimize the number of violated constraints */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* free the allocated memory for the copied constraint array */
   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}


/*
 *  constraint-dependent dialog entries
 */

/** dialog execution method for the SCIPtransformMinUC() method */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeMinUC)
{  /*lint --e{715}*/
   SCIP_Bool success;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );
   SCIPdialogMessage(scip, NULL, "\n");

   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;
   case SCIP_STAGE_PROBLEM:
      SCIPdialogMessage(scip, NULL, "change problem to MinUC\n");
      SCIPdialogMessage(scip, NULL, "==============\n");

      SCIP_CALL( SCIPtransformMinUC(scip, &success) );

      if( !success )
      {
         SCIPdialogMessage(scip, NULL, "some constraints could not be transformed to superindicator constraints and were removed\n");
      }

      SCIPdialogMessage(scip, NULL, "\n");
      SCIPdialogMessage(scip, NULL, "changed problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
         SCIPgetNVars(scip), SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNImplVars(scip), SCIPgetNContVars(scip),
         SCIPgetNConss(scip));

      SCIPdialogMessage(scip, NULL, "increased branching priority of new binary variables");

      break;
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
      SCIPdialogMessage(scip, NULL, "problem has to be in problem stage to create MinUC problem\n");
      break;
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }  /*lint --e{616}*/

   SCIPdialogMessage(scip, NULL, "\n");
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}
