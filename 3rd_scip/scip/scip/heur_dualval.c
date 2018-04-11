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

/**@file   heur_dualval.c
 * @brief  dualval primal heuristic
 * @author Tobias Buchwald
 *
 * This heuristic tries to find solutions by taking the LP or NLP, rounding solution values, fixing the variables to the
 * rounded values and then changing some of the values.To determine which variable is changed we give each variable a
 * ranking dependent on its dualvalue.  We work with a transformed problem that is always feasible and has objective = 0
 * iff the original problem is also feasible. Thus we cannot expect to find really good solutions.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/heur_dualval.h"
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/cons_indicator.h"
#include "scip/cons_varbound.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_knapsack.h"

#include "nlpi/nlpi.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/nlpi_ipopt.h"
#include "nlpi/exprinterpret.h"

#define HEUR_NAME             "dualval"
#define HEUR_DESC             "primal heuristic using dual values"
#define HEUR_DISPCHAR         'Y'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define EVENTHDLR_NAME        "lpsol_dualval"
#define EVENTHDLR_DESC        "event handler for lp solution found"

/* default values for user parameters */
/* boolean parameters */
#define DEFAULT_FORCEIMPROVEMENTS   FALSE    /**< exit if objective doesn't improve */
#define DEFAULT_ONLYCHEAPER         TRUE     /**< add constraint to ensure that discrete vars are improving */
#define DEFAULT_ONLYLEAVES          FALSE    /**< disable the heuristic if it was not called at a leaf of the B&B tree */
#define DEFAULT_RELAXINDICATORS     FALSE    /**< relax the indicator variables by introducing continuous copies */
#define DEFAULT_RELAXCONTVARS       FALSE    /**< enable relaxation of continous variables */

/* integer parameters */
#define DEFAULT_HEURVERBLEVEL       0        /**< verblevel of the heuristic, default is 0 to display nothing */
#define DEFAULT_NLPVERBLEVEL        0        /**< verblevel of the nlp solver, can be 0 or 1 */
#define DEFAULT_RANKVALUE           10       /**< number of ranks that should be displayed when the heuristic is called */
#define DEFAULT_MAXCALLS            25       /**< maximal number of recursive calls of the heuristic (if dynamicdepth is off) */
#define DEFAULT_DYNAMICDEPTH        0        /**< says if and how the recursion depth is computed at runtime */
#define DEFAULT_MAXEQUALRANKS       50       /**< maximal number of variables that may have maximal rank, quit if there are more, turn off by setting -1 */

/* real value parameters */
#define DEFAULT_MINGAP              5.0      /**< minimal gap for which we still run the heuristic, if gap is less we return without doing anything */
#define DEFAULT_LAMBDASLACK         1.0      /**< value added to objective of slack variables, must not be zero */
#define DEFAULT_LAMBDAOBJ           0.0      /**< scaling factor for the objective function */


/**primal heuristic data */
struct SCIP_HeurData
{
   SCIP*                 subscip;            /**< copy of CIP */
   SCIP_VAR**            integervars;        /**< array of all binary and integer variables of the original scip */
   SCIP_HASHMAP*         varsciptosubscip;   /**< mapping variables in SCIP to sub-SCIP variables */
   SCIP_HASHMAP*         varsubsciptoscip;   /**< mapping variables in sub-SCIP to SCIP variables */
   SCIP_HASHMAP*         origsubscipConsMap; /**< maps constraints from the transformed problem to corresponding constraints in subproblem */
   SCIP_HASHMAP*         switchedvars;       /**< stores the last value of switched var to avoid cycling */
   SCIP_HASHMAP*         switchedvars2;      /**< stores the second last value of switched vars to avoid cycling */
   SCIP_HASHMAP*         relaxcons;          /**< maps subscip variables to their relaxation constraints */
   SCIP_HASHMAP*         relaxconsindi;      /**< maps indicator variables and their copies to relaxation constraint */
   SCIP_HASHMAP*         slacktoindivarsmap; /**< maps slack variables of indicator constraint to indicator variable */
   SCIP_HASHMAP*         indicators;         /**< maps indicator variables to their indicator constraint */
   SCIP_HASHMAP*         conss2nlrow;        /**< maps constraint to the corresponding nlrow */
   SCIP_HASHMAP*         dualvalues;         /**< maps constraints of the subscip to their dual values */
   SCIP_HASHMAP*         slack2var;          /**< maps slack variables to the variable they actually relax */
   SCIP_HASHMAP*         indicopymap;        /**< maps indicator variables to their copy variables */
   SCIP_HASHMAP*         indicopymapback;    /**< maps copy variables to their indicator variables */
   SCIP_HASHMAP*         slackvarlbMap;      /**< mapping used indicators to slack variables lower bound*/
   SCIP_HASHMAP*         slackvarubMap;      /**< mapping used indicators to slack variables upper bound*/
   SCIP_CONS*            objbound;           /**< contraint for upper bound of the objective function */
   SCIP_Real             prevobjective;      /**< stores objective value (of the original) so we know if it improved */
   SCIP_Real             mingap;             /**< don't run the heuristic if the gap is less than mingap */
   SCIP_Real             lambdaslack;        /**< the value added to the objective function */
   SCIP_Real             lambdaobj;          /**< the value the original objective function is scaled with */
   int                   integervarssize;    /**< size of integervars array */
   int                   nintegervars;       /**< number of integer variables in the original problem */
   int                   heurverblevel;      /**< verblevel, range is 0 to 4 */
   int                   nlpverblevel;       /**< sets verblevel of the included nlp */
   int                   rankvalue;          /**< print out the 'rankvalue' highest ranks during iterations */
   int                   maxcalls;           /**< maximum number of allowed iterations */
   int                   nonimprovingRounds; /**< nr of rounds, where the algorithm has not improved */
   int                   dynamicdepth;       /**< how should the number of calls be computed? */
   int                   maxequalranks;      /**< maximum number of variables that may have maximal (absolute) rank */
   int                   nvars;              /**< number of active transformed variables in SCIP */
   int                   nsubvars;           /**< number of original variables in sub-SCIP */
   int                   usedcalls;          /**< number of currently used iterations */
   SCIP_Bool             isnlp;              /**< tells us, whether we have nonlinearities in our program or not */
   SCIP_Bool             forceimprovements;  /**< whether we exit on nonimproving objective in the relaxation or not */
   SCIP_Bool             prevInfeasible;     /**< will tell us if the previous call led to an infeasible fixing */
   SCIP_Bool             solfound;           /**< parameter says, if we already found a solution and have to go back */
   SCIP_Bool             subscipisvalid;     /**< whether all constraints have been copied */
   SCIP_Bool             switchdifferent;    /**< tells us that we want to go up one level and switch another variable */
   SCIP_Bool             triedsetupsubscip;  /**< whether we have tried to setup a sub-SCIP */
   SCIP_Bool             onlycheaper;        /**< add constraint to ensure that discrete vars are improving */
   SCIP_Bool             onlyleaves;         /**< don't use heuristic if we are not in a leaf of the B&B tree */
   SCIP_Bool             relaxindicators;    /**< additionally relax indicator variables */
   SCIP_Bool             relaxcontvars;      /**< additionally relax continous variables */
};

/*
 * event handler method
 */

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitLPsol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* notify SCIP that your event handler wants to react on the event type best solution found */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_FIRSTLPSOLVED | SCIP_EVENTTYPE_LPSOLVED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitLPsol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);

   /* notify SCIP that your event handler wants to drop the event type best solution found */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_FIRSTLPSOLVED | SCIP_EVENTTYPE_LPSOLVED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecLPsol)
{  /*lint --e{715}*/
   int i;
   int nsubconss;
   SCIP_HEURDATA* heurdata;
   SCIP_CONS** subconss;
   SCIP_Real* dualval;

   assert(eventhdlr != NULL);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   heurdata = (SCIP_HEURDATA* )SCIPeventhdlrGetData(eventhdlr);
   nsubconss = SCIPgetNOrigConss(heurdata->subscip);
   subconss = SCIPgetOrigConss(heurdata->subscip);

   /* free memory of all entries and clear the hashmap before filling it */
   for( i = 0; i < nsubconss; i++ )
   {
      dualval = (SCIP_Real*)SCIPhashmapGetImage(heurdata->dualvalues, subconss[i]);
      if( dualval != NULL )
         SCIPfreeBlockMemoryArray(heurdata->subscip, &dualval, 1);
   }
   SCIP_CALL( SCIPhashmapRemoveAll(heurdata->dualvalues) );

   /* insert dualvalues from LP into a hashmap */
   for( i = 0; i < nsubconss; i++ )
   {
      SCIP_CONS* transcons = NULL;
      SCIP_CALL( SCIPgetTransformedCons(heurdata->subscip, subconss[i], &transcons) );

      if( transcons == NULL )
         continue;

      if( SCIPconsGetHdlr(transcons) != SCIPfindConshdlr(heurdata->subscip, "linear") )
         continue;

      SCIP_CALL( SCIPallocBlockMemoryArray(heurdata->subscip, &dualval, 1) ); /*lint !e506*/
      *dualval = -SCIPgetDualsolLinear(heurdata->subscip, transcons );
      SCIP_CALL( SCIPhashmapInsert(heurdata->dualvalues, subconss[i], dualval) );
   }
   if( heurdata->heurverblevel > 2 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "LP solved event!\n");

   return SCIP_OKAY;
}

/** includes event handler for best solution found */
static
SCIP_RETCODE SCIPincludeEventHdlrLPsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr = NULL;

   eventhdlrdata = (SCIP_EVENTHDLRDATA*)heurdata;

   /* create event handler */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecLPsol, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitLPsol) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitLPsol) );

   return SCIP_OKAY;
}

/*
 * Local methods
 */

/** releases all variables or constraints from given hash map */
static
SCIP_RETCODE releaseHashmapEntries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         hashmap,            /**< hashmap */
   SCIP_Bool             isvarmap            /**< are the entries variables or constraints? */
   )
{
   int nentries;
   int i;

   assert(scip != NULL);
   assert(hashmap != NULL);

   nentries = SCIPhashmapGetNEntries(hashmap);

   for( i = 0; i < nentries; ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      entry = SCIPhashmapGetEntry(hashmap, i);

      if( entry != NULL )
      {
         if( isvarmap )
         {
            SCIP_VAR* var;
            var = (SCIP_VAR*) SCIPhashmapEntryGetImage(entry);

            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
         else
         {
            SCIP_CONS* cons;
            cons = (SCIP_CONS*) SCIPhashmapEntryGetImage(entry);

            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
   }

   return SCIP_OKAY;
}

/** releases all NLP rows from given hash map */
static
SCIP_RETCODE releaseHashmapNLPRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         hashmap             /**< hashmap */
   )
{
   int nentries;
   int i;

   assert(scip != NULL);
   assert(hashmap != NULL);

   nentries = SCIPhashmapGetNEntries(hashmap);

   for( i = 0; i < nentries; ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      entry = SCIPhashmapGetEntry(hashmap, i);
      if( entry != NULL )
      {
         SCIP_NLROW* nlrow;
         nlrow = (SCIP_NLROW*) SCIPhashmapEntryGetImage(entry);

         SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
      }
   }

   return SCIP_OKAY;
}


/** adds linear constraints from a SCIP instance to its NLP */
static
SCIP_RETCODE addLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for linear constraints */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints to NLP */
   SCIP_Bool             addcontconss,       /**< whether to add continuous    linear constraints to NLP */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_CONS**   conss;
   SCIP_VAR**    vars;
   SCIP_NLROW*   nlrow;
   int           nconss;
   int           i;
   int           j;
   int           nvars;
   SCIP_Bool     iscombinatorial;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNActiveConss(conshdlr);
   conss  = SCIPconshdlrGetConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      /* under some circumstances, this method may be called even though the problem has been shown to be
       * infeasible in presolve already.
       * this infeasibility may come from a linear constraint with lhs > rhs
       * the NLP does not allow such constraints, so we skip them here
       */
      if( !SCIPisRelLE(scip, SCIPgetLhsLinear(scip, conss[i]), SCIPgetRhsLinear(scip, conss[i])) )
         continue;

      nvars = SCIPgetNVarsLinear(scip, conss[i]);
      vars  = SCIPgetVarsLinear(scip, conss[i]);

      /* check if constraint should be added, only need this check if we do not wanna any constraint anyway */
      if( !addcombconss || !addcontconss )
      {
         iscombinatorial = TRUE;

         for( j = 0; j < nvars; ++j )
         {
            if( SCIPvarGetType(vars[j]) >= SCIP_VARTYPE_CONTINUOUS )
            {
               iscombinatorial = FALSE;
               break;
            }
         }

         /* skip constraint, if not of interest */
         if( (iscombinatorial && !addcombconss) || (!iscombinatorial && !addcontconss) )
            continue;
      }

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            SCIPgetNVarsLinear(scip, conss[i]), SCIPgetVarsLinear(scip, conss[i]), SCIPgetValsLinear(scip, conss[i]),
            0, NULL, 0, NULL, NULL,
            SCIPgetLhsLinear(scip, conss[i]), SCIPgetRhsLinear(scip, conss[i]),
            SCIP_EXPRCURV_LINEAR) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->conss2nlrow, conss[i], nlrow) );
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }

   return SCIP_OKAY;
}

/** adds variable bound constraints from a SCIP instance to its NLP */
static
SCIP_RETCODE addVarboundConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for linear constraints */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints to NLP */
   SCIP_Bool             addcontconss,       /**< whether to add continuous    linear constraints to NLP */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   SCIP_VAR*     vars[2];
   SCIP_Real     coefs[2];
   SCIP_Bool     iscombinatorial;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNActiveConss(conshdlr);
   conss  = SCIPconshdlrGetConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      vars[0] = SCIPgetVarVarbound(scip, conss[i]);
      vars[1] = SCIPgetVbdvarVarbound(scip, conss[i]);

      iscombinatorial = SCIPvarGetType(vars[0]) < SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(vars[1]) < SCIP_VARTYPE_CONTINUOUS;

      /* skip constraint, if not of interest */
      if( (iscombinatorial && !addcombconss) || (!iscombinatorial && !addcontconss) )
         continue;

      coefs[0] = 1.0;
      coefs[1] = SCIPgetVbdcoefVarbound(scip, conss[i]);

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            2, vars, coefs,
            0, NULL, 0, NULL, NULL,
            SCIPgetLhsVarbound(scip, conss[i]), SCIPgetRhsVarbound(scip, conss[i]),
            SCIP_EXPRCURV_LINEAR) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->conss2nlrow, conss[i], nlrow) );
   }

   return SCIP_OKAY;
}


/** adds logic-or constraints to NLP */
static
SCIP_RETCODE addLogicOrConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for linear constraints */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   int           j;
   SCIP_Real*    coefs;
   int           coefssize;
   int           nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNActiveConss(conshdlr);
   if( !nconss )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);

   coefs = NULL;
   coefssize = 0;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      nvars = SCIPgetNVarsLogicor(scip, conss[i]);

      if( coefssize < nvars )
      {
         if( coefs == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, nvars) );
         }
         for( j = coefssize; j < nvars; ++j )
            coefs[j] = 1.0;
         coefssize = nvars;
      }

      /* logic or constraints: 1 == sum_j x_j */
      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            nvars, SCIPgetVarsLogicor(scip, conss[i]), coefs,
            0, NULL, 0, NULL, NULL,
            1.0, SCIPinfinity(scip),
            SCIP_EXPRCURV_LINEAR) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->conss2nlrow, conss[i], nlrow) );
   }

   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** adds setppc constraints to NLP */
static
SCIP_RETCODE addSetppcConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for linear constraints */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   int           j;
   SCIP_Real*    coefs;
   int           coefssize;
   int           nvars;
   SCIP_Real     lhs;
   SCIP_Real     rhs;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNActiveConss(conshdlr);
   if( nconss == 0 )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);

   coefs = NULL;
   coefssize = 0;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      nvars = SCIPgetNVarsSetppc(scip, conss[i]);

      if( coefssize < nvars )
      {
         if( coefs == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, nvars) );
         }
         for( j = coefssize; j < nvars; ++j )
            coefs[j] = 1.0;
         coefssize = nvars;
      }

      /* setppc constraint: 1 ~ sum_j x_j */

      switch( SCIPgetTypeSetppc(scip, conss[i]) )
      {
      case SCIP_SETPPCTYPE_PARTITIONING:
         lhs = 1.0;
         rhs = 1.0;
         break;

      case SCIP_SETPPCTYPE_PACKING:
         lhs = -SCIPinfinity(scip);
         rhs = 1.0;
         break;

      case SCIP_SETPPCTYPE_COVERING:
         lhs = 1.0;
         rhs = SCIPinfinity(scip);
         break;

      default:
         SCIPerrorMessage("unexpected setppc type\n");
         return SCIP_ERROR;
      }

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            nvars, SCIPgetVarsSetppc(scip, conss[i]), coefs,
            0, NULL, 0, NULL, NULL,
            lhs, rhs,
            SCIP_EXPRCURV_LINEAR) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->conss2nlrow, conss[i], nlrow) );
   }

   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** adds knapsack constraints to NLP */
static
SCIP_RETCODE addKnapsackConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for linear constraints */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   int           j;
   SCIP_Real*    coefs;
   int           coefssize;
   int           nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNActiveConss(conshdlr);
   if( nconss == 0 )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   coefs = NULL;
   coefssize = 0;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_Longint* weights;

      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      nvars = SCIPgetNVarsKnapsack(scip, conss[i]);

      if( coefssize < nvars )
      {
         if( coefs == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, nvars) );
         }
         coefssize = nvars;
      }

      weights = SCIPgetWeightsKnapsack(scip, conss[i]);
      for( j = 0; j < nvars; ++j )
         coefs[j] = (SCIP_Real)weights[j];  /*lint !e613*/

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            nvars, SCIPgetVarsKnapsack(scip, conss[i]), coefs,
            0, NULL, 0, NULL, NULL,
            -SCIPinfinity(scip), (SCIP_Real)SCIPgetCapacityKnapsack(scip, conss[i]),
            SCIP_EXPRCURV_LINEAR) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->conss2nlrow, conss[i], nlrow) );
   }

   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}


/** adds combinatorial and/or continuous variants of linear constraints from a SCIP instance to its NLP */
static
SCIP_RETCODE addLinearConstraintsToNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints to NLP */
   SCIP_Bool             addcontconss,       /**< whether to add continuous    linear constraints to NLP */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* add linear constraints */
   conshdlr = SCIPfindConshdlr(scip, "linear");
   if( conshdlr != NULL )
   {
      SCIP_CALL( addLinearConstraints(scip, conshdlr, addcombconss, addcontconss, heurdata) );
   }

   /* add varbound constraints */
   conshdlr = SCIPfindConshdlr(scip, "varbound");
   if( conshdlr != NULL )
   {
      SCIP_CALL( addVarboundConstraints(scip, conshdlr, addcombconss, addcontconss, heurdata) );
   }

   if( addcombconss )
   {
      /* add logic-or constraints */
      conshdlr = SCIPfindConshdlr(scip, "logicor");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addLogicOrConstraints(scip, conshdlr, heurdata) );
      }

      /* add setppc constraints */
      conshdlr = SCIPfindConshdlr(scip, "setppc");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addSetppcConstraints(scip, conshdlr, heurdata) );
      }

      /* add knapsack constraints */
      conshdlr = SCIPfindConshdlr(scip, "knapsack");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addKnapsackConstraints(scip, conshdlr, heurdata) );
      }
   }

   return SCIP_OKAY;
}



/** creates a SCIP_SOL in our SCIP space out of the SCIP_SOL from a sub-SCIP */
static
SCIP_RETCODE createSolFromSubScipSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL**            sol,                /**< buffer to store solution value; if pointing to NULL, a new solution
                                              *   is created, otherwise values in the given one are overwritten */
   SCIP_SOL*             subsol              /**< solution of sub-SCIP */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR**     vars;
   SCIP_VAR**     subvars;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   SCIP_Real      scalar;
   SCIP_Real      constant;
   SCIP_Real      val;
   int            i;
   int            nvars;

   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );
   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nvars, NULL, NULL, NULL, NULL) );

   if( *sol == NULL )
   {
      SCIP_CALL( SCIPcreateOrigSol(scip, sol, heur) );
   }

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);

   for( i = 0; i < nvars; ++i )
   {
      var = vars[i];

      constant = 0;
      scalar   = 1.0;
      var = SCIPvarGetTransVar(var);
      val = 0;

      if( REALABS(scalar) > 0 )
      {
         SCIP_Real transval = 0.0;

         subvar = (SCIP_VAR*) SCIPhashmapGetImage(heurdata->varsciptosubscip, (void*)var);
         if( subvar == NULL )
         {
            SCIPdebugMsg(scip, "return14 : abort building solution since a variable was not in our list\n");

            SCIP_CALL( SCIPfreeSol(scip, sol) );
            return SCIP_OKAY;
         }

         if( SCIPvarIsBinary(subvar) )
            transval = SCIPvarGetLbGlobal(subvar);
         else
         {
            SCIP_Real tconstant = 0.0;
            SCIP_Real tscalar   = 1.0;
            SCIP_CALL( SCIPgetProbvarSum(heurdata->subscip, &subvar, &tscalar, &tconstant) );

            transval = 0.0;

            if( REALABS(tscalar) > 0.0 )
            {
               assert(subvar != NULL);
               transval = SCIPgetSolVal(heurdata->subscip, subsol, subvar);
            }

            /* recompute aggregations */
            transval = tscalar * transval + tconstant;
         }
         val = scalar * transval + constant;
      }
      else
      {
         /* recompute aggregations */
         val = scalar * val + constant;
      }

      assert( val != SCIP_INVALID ); /*lint !e777*/
      SCIP_CALL( SCIPsetSolVal(scip, *sol, vars[i], val) );
   }

   return SCIP_OKAY;
}



/** creates copy of CIP from problem in SCIP */
static
SCIP_RETCODE createSubSCIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_HASHMAP* varsmap;
   SCIP_HASHMAP* conssmap;
   SCIP_CONSHDLR*  conshdlrindicator;
   SCIP_CONSHDLR*  conshdlrindi;
   SCIP_CONSHDLR*  conshdlrlin;
   SCIP_CONSHDLR*  conshdlrabspow;
   SCIP_CONSHDLR*  conshdlrquad;
   SCIP_CONSHDLR*  conshdlrnonlin;
   SCIP_CONSHDLR*  conshdlrvarbound;
   SCIP_CONSHDLR*  conshdlrknapsack;
   SCIP_CONSHDLR*  conshdlrlogicor;
   SCIP_CONSHDLR*  conshdlrsetppc;
   SCIP_CONSHDLR*  currentconshdlr;
   SCIP_CONSHDLR*  conshdlrsignpower;
   SCIP_CONS**  conss;
   SCIP_CONS*   subcons;
   SCIP_CONS*   transcons;
   SCIP_CONS*   linindicons;
   SCIP_CONS*   indicons;
   SCIP_CONS*   cons = NULL;
   SCIP_VAR**   vars;
   SCIP_VAR**   subvars;
   SCIP_VAR*    var;
   SCIP_VAR*    tmpvar;
   SCIP_VAR*    subvar;
   SCIP_VAR*    slackvarpos;
   SCIP_VAR*    slackvarneg;
   SCIP_VAR*    indislackvarpos;
   SCIP_VAR*    indislackvarneg;
   SCIP_VAR*    indicatorcopy;
   char         probname[SCIP_MAXSTRLEN];
   char         varname[SCIP_MAXSTRLEN];
   char         consname[SCIP_MAXSTRLEN];
   SCIP_Real    varobjective;
   int          nconss;
   int          nconsindicator;
   int          i;
   int          j;
   int          k;
   int          nvars;
   int          ncontvars;
   SCIP_Bool    feasible;
   SCIP_Bool    success;

   assert( heurdata != NULL );
   assert( heurdata->subscip == NULL );

   heurdata->usedcalls = 0;
   heurdata->solfound = FALSE;
   heurdata->nonimprovingRounds = 0;

   /* we can't change the vartype in some constraints, so we have to check that only the right constraints are present*/
   conshdlrindi = SCIPfindConshdlr(scip, "indicator");
   conshdlrlin = SCIPfindConshdlr(scip, "linear");
   conshdlrabspow = SCIPfindConshdlr(scip, "abspower");
   conshdlrquad = SCIPfindConshdlr(scip, "quadratic");
   conshdlrnonlin = SCIPfindConshdlr(scip, "nonlinear");
   conshdlrvarbound = SCIPfindConshdlr(scip, "varbound");
   conshdlrknapsack = SCIPfindConshdlr(scip, "knapsack");
   conshdlrlogicor = SCIPfindConshdlr(scip, "logicor");
   conshdlrsetppc = SCIPfindConshdlr(scip, "setppc");
   conshdlrsignpower = SCIPfindConshdlr(scip, "signpower");

   nconss = SCIPgetNOrigConss(scip);
   conss = SCIPgetOrigConss(scip);

   /* for each constraint ask if it has an allowed type */
   for (i = 0; i < nconss; i++ )
   {
      cons = conss[i];
      currentconshdlr = SCIPconsGetHdlr(cons);

      if( currentconshdlr == conshdlrindi ||
         currentconshdlr == conshdlrabspow ||
         currentconshdlr == conshdlrquad ||
         currentconshdlr == conshdlrnonlin ||
         currentconshdlr == conshdlrvarbound ||
         currentconshdlr == conshdlrknapsack ||
         currentconshdlr == conshdlrlogicor ||
         currentconshdlr == conshdlrsetppc ||
         currentconshdlr == conshdlrlin ||
         currentconshdlr == conshdlrsignpower)
      {
         continue;
      }
      else
      {
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, &ncontvars) );

   if( heurdata->dynamicdepth == 1 )
   {
      heurdata->maxcalls = (int)SCIPfloor(scip, sqrt((double)(nvars - ncontvars)));
   }

   heurdata->triedsetupsubscip = TRUE;

   /* initializing the subproblem */
   SCIP_CALL( SCIPcreate(&heurdata->subscip) );

   /* create variable hash mapping scip -> subscip */
   SCIP_CALL( SCIPhashmapCreate(&varsmap, SCIPblkmem(scip), nvars) );

   SCIP_CALL( SCIPhashmapCreate(&heurdata->switchedvars, SCIPblkmem(scip), heurdata->maxcalls) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->switchedvars2, SCIPblkmem(scip), heurdata->maxcalls) );

   /* create sub-SCIP copy of CIP, copy interesting plugins */
   success = TRUE;
   SCIP_CALL( SCIPcopyPlugins(scip, heurdata->subscip, TRUE, FALSE, TRUE, FALSE, TRUE,
         FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, &success) );

   if( success == FALSE )
   {
      SCIPdebugMsg(scip, "In heur_dualval: failed to copy some plugins to sub-SCIP, continue anyway\n");
   }

   /* copy parameter settings */
   SCIP_CALL( SCIPcopyParamSettings(scip, heurdata->subscip) );

   /* create problem in sub-SCIP */

   /* get name of the original problem and add "dualval" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_dualval", SCIPgetProbName(scip));
   SCIP_CALL( SCIPcreateProb(heurdata->subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPincludeEventHdlrLPsol(heurdata->subscip, heurdata) );

   /* copy all variables */
   SCIP_CALL( SCIPcopyVars(scip, heurdata->subscip, varsmap, NULL, NULL, NULL, 0, TRUE) );

   /* copy as many constraints as possible */
   SCIP_CALL( SCIPhashmapCreate(&conssmap, SCIPblkmem(scip), SCIPgetNConss(scip)) );
   SCIP_CALL( SCIPcopyConss(scip, heurdata->subscip, varsmap, conssmap, TRUE, FALSE, &heurdata->subscipisvalid) );

   SCIP_CALL( SCIPhashmapCreate(&heurdata->origsubscipConsMap, SCIPblkmem(scip), SCIPgetNConss(scip)) );

   nconss = SCIPgetNOrigConss(scip);
   conss = SCIPgetOrigConss(scip);

   /* fill constraint mapping from original scip to the subscip */
   for( i = 0; i < nconss; ++i )
   {
      transcons = NULL;
      SCIP_CALL( SCIPgetTransformedCons(scip, conss[i], &transcons) );

      subcons = (SCIP_CONS*)SCIPhashmapGetImage(conssmap, transcons);
      assert( subcons != NULL );

      SCIP_CALL( SCIPcaptureCons(heurdata->subscip, subcons) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->origsubscipConsMap, transcons, subcons) );
   }

   SCIP_CALL( SCIPhashmapCreate(&heurdata->conss2nlrow, SCIPblkmem(scip), SCIPgetNConss(scip)) );

   if( !heurdata->subscipisvalid )
   {
      SCIPdebugMsg(scip, "In heur_dualval: failed to copy some constraints to sub-SCIP, continue anyway\n");
   }

   SCIP_CALL( SCIPgetVarsData(heurdata->subscip, &subvars, &heurdata->nsubvars, NULL, NULL, NULL, NULL) );
   heurdata->nvars = nvars;

   /* create hashmaps from scip transformed vars to subscip original vars, and vice versa
    * capture variables in SCIP and sub-SCIP
    * catch global bound change events */
   SCIP_CALL( SCIPhashmapCreate(&heurdata->varsubsciptoscip, SCIPblkmem(scip), SCIPgetNOrigVars(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->varsciptosubscip, SCIPblkmem(scip), SCIPgetNOrigVars(scip)) );

   /* we need to get all subscip variables, also those which are copies of fixed variables from the main scip,
    * therefore we iterate over the hashmap */
   for( i = 0; i < SCIPhashmapGetNEntries(varsmap); ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      entry = SCIPhashmapGetEntry(varsmap, i);
      if( entry != NULL )
      {
         var    = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
         subvar = (SCIP_VAR*) SCIPhashmapEntryGetImage(entry);

         assert( SCIPvarGetProbindex(subvar) >= 0 );
         assert( SCIPvarGetProbindex(subvar) <= heurdata->nsubvars );

         if( SCIPvarIsActive(var) )
         {
            assert( SCIPvarGetProbindex(var) <= heurdata->nvars );
            /* assert that we have no mapping for this var yet */
            assert( SCIPhashmapGetImage(heurdata->varsciptosubscip,var) == NULL );
            SCIP_CALL( SCIPhashmapInsert(heurdata->varsciptosubscip, var, subvar) );
         }

         assert( SCIPhashmapGetImage(heurdata->varsubsciptoscip, subvar) == NULL );
         SCIP_CALL( SCIPhashmapInsert(heurdata->varsubsciptoscip, subvar, var) );

         SCIP_CALL( SCIPcaptureVar(scip, var) );
         SCIP_CALL( SCIPcaptureVar(heurdata->subscip, subvar) );

         assert( SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetLbGlobal(subvar)) );
         assert( SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), SCIPvarGetUbGlobal(subvar)) );
      }
   }

   /* we map all slack variables of indicator constraints to their indicator variables */
   conshdlrindicator  = SCIPfindConshdlr(scip, "indicator");
   nconsindicator = SCIPconshdlrGetNConss(conshdlrindicator);

   SCIP_CALL( SCIPhashmapCreate(&heurdata->slacktoindivarsmap, SCIPblkmem(scip), nconsindicator) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->indicators, SCIPblkmem(scip), nconsindicator) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->indicopymap, SCIPblkmem(scip), nconsindicator) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->indicopymapback, SCIPblkmem(scip), nconsindicator) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->slackvarlbMap, SCIPblkmem(scip), SCIPgetNOrigVars(scip)) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->slackvarubMap, SCIPblkmem(scip), SCIPgetNOrigVars(scip)) );

   for( i = 0; i < nconsindicator; i++ )
   {
      SCIP_CONS** indicatorconss = SCIPconshdlrGetConss(conshdlrindicator);
      SCIP_CONS* currcons;

      currcons = indicatorconss[i];
      assert(currcons != NULL);

      SCIP_CALL( SCIPhashmapInsert(heurdata->slacktoindivarsmap, SCIPgetSlackVarIndicator(currcons),
            SCIPgetBinaryVarIndicator(currcons)) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->indicators, SCIPgetBinaryVarIndicator(currcons), currcons) );
      SCIP_CALL( SCIPcaptureCons(scip, currcons) );
      SCIP_CALL( SCIPcaptureVar(scip, SCIPgetBinaryVarIndicator(currcons)) );
   }

   /* we introduce slackvariables s+ and s- for each constraint to ensure that the problem is feasible
    * we want to minimize over the sum of these variables, so set the objective to 1 */
   SCIP_CALL( SCIPhashmapCreate(&heurdata->relaxcons, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->relaxconsindi, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&heurdata->slack2var, SCIPblkmem(scip), nvars) );

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(heurdata->integervars), nvars) );
   BMSclearMemoryArray(heurdata->integervars, nvars);
   heurdata->integervarssize = nvars;
   j = 0;

   /* here we relax the variables (or indicator constraints, since indicator variables cannot be relaxed) */
   for( i = 0; i < nvars; ++i )
   {
      var = SCIPvarGetTransVar(vars[i]);
      assert( var != NULL );

      if( ! SCIPvarIsActive(var) )
         continue;

      if( ! SCIPvarIsIntegral(var) )
         continue;

      heurdata->integervars[j++] = vars[i];

      var = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsciptosubscip, var);
      assert( var != NULL );

      /* in this case our variable is an indicator variable */
      if( SCIPhashmapGetImage(heurdata->indicators, SCIPhashmapGetImage(heurdata->varsubsciptoscip, var)) != NULL )
      {
         /* we have to find all the indicator constraints of this variable */
         for( k = 0; k < nconsindicator; k++ )
         {
            SCIP_CONS**  indicatorconss = SCIPconshdlrGetConss(conshdlrindicator);
            SCIP_CONS*   currcons;
            SCIP_VAR*    negatedvar;
            SCIP_VAR*    indicatorbinvar;

            currcons = indicatorconss[k];
            assert(currcons != NULL);

            indicatorbinvar = SCIPgetBinaryVarIndicator(currcons);
            assert(indicatorbinvar != NULL);

            SCIP_CALL( SCIPgetNegatedVar(scip, (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsubsciptoscip, var), &negatedvar) );

            if( indicatorbinvar == SCIPhashmapGetImage(heurdata->varsubsciptoscip, var) || indicatorbinvar == negatedvar )
            {
               /* case that we have a negated variable */
               if( SCIPvarIsNegated(indicatorbinvar) )
               {
                  assert(indicatorbinvar == negatedvar);
                  varobjective = SCIPvarGetObj(negatedvar);
               }
               else
               {
                  assert(indicatorbinvar != negatedvar);
                  varobjective = SCIPvarGetObj(indicatorbinvar);
               }

               varobjective = heurdata->lambdaobj * REALABS(varobjective);

               indicons = currcons;
               assert( indicons != NULL );

               indicons = (SCIP_CONS*)SCIPhashmapGetImage(conssmap, indicons);

               assert( indicons != NULL );
               linindicons = SCIPgetLinearConsIndicator(indicons);

               (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_pos3", SCIPconsGetName(linindicons));
               SCIP_CALL( SCIPcreateVar(heurdata->subscip, &slackvarpos, varname, 0.0, SCIPinfinity(heurdata->subscip),
                     heurdata->lambdaslack *100 + varobjective, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(heurdata->subscip, slackvarpos) );

               (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_neg3", SCIPconsGetName(linindicons));
               SCIP_CALL( SCIPcreateVar(heurdata->subscip, &slackvarneg, varname, 0.0, SCIPinfinity(heurdata->subscip),
                     heurdata->lambdaslack * 100 + varobjective, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(heurdata->subscip, slackvarneg) );

               /* make a copy of the indicator to relax it if this parameter is set true */
               if( heurdata->relaxindicators )
               {
                  SCIP_CONS* imagecons;

                  indicatorbinvar = SCIPgetBinaryVarIndicator(indicons);

                  SCIP_CALL( SCIPgetNegatedVar(heurdata->subscip, indicatorbinvar, &negatedvar) );

                  if( SCIPhashmapGetImage(heurdata->indicopymap, indicatorbinvar) == NULL &&
                     SCIPhashmapGetImage(heurdata->indicopymap, negatedvar) == NULL)
                  {
                     SCIP_Bool negated = FALSE;

                     if (SCIPvarIsNegated(indicatorbinvar))
                     {
                        indicatorbinvar = negatedvar;
                        negated = TRUE;
                     }

                     (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "indicopy_%s", SCIPvarGetName(indicatorbinvar));
                     SCIP_CALL( SCIPcreateVar(heurdata->subscip, &indicatorcopy, varname, SCIPvarGetLbGlobal(indicatorbinvar), SCIPvarGetUbGlobal(indicatorbinvar),
                           SCIPvarGetObj(indicatorbinvar), SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

                     SCIP_CALL( SCIPaddVar(heurdata->subscip, indicatorcopy) );

                     SCIP_CALL( SCIPhashmapInsert(heurdata->indicopymap, indicatorbinvar, indicatorcopy) );
                     SCIP_CALL( SCIPhashmapInsert(heurdata->indicopymapback, indicatorcopy, indicatorbinvar) );
                     SCIP_CALL( SCIPcaptureVar(heurdata->subscip, indicatorbinvar) );

                     (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_pos1", SCIPvarGetName(indicatorbinvar));
                     SCIP_CALL( SCIPcreateVar(heurdata->subscip, &indislackvarpos, varname, 0.0, SCIPinfinity(heurdata->subscip),
                           heurdata->lambdaslack * 100 + varobjective, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
                     SCIP_CALL( SCIPaddVar(heurdata->subscip, indislackvarpos) );

                     (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_neg1", SCIPvarGetName(indicatorbinvar));
                     SCIP_CALL( SCIPcreateVar(heurdata->subscip, &indislackvarneg, varname, 0.0, SCIPinfinity(heurdata->subscip),
                           heurdata->lambdaslack * 100 + varobjective, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
                     SCIP_CALL( SCIPaddVar(heurdata->subscip, indislackvarneg) );

                     /* create linking constraint */
                     (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "linking_%s", SCIPvarGetName(indicatorbinvar));
                     cons = NULL;
                     SCIP_CALL( SCIPcreateConsLinear( heurdata->subscip, &cons, consname, 0, NULL, NULL, 0.0, 0.0,
                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
                     SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, indicatorbinvar, 1.0) );
                     SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, indicatorcopy, -1.0) );
                     SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, indislackvarpos, 1.0) );
                     SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, indislackvarneg, -1.0) );

                     SCIP_CALL( SCIPhashmapInsert(heurdata->relaxconsindi, indicatorbinvar, cons) );
                     SCIP_CALL( SCIPhashmapInsert(heurdata->relaxconsindi, indicatorcopy, cons) );

                     SCIP_CALL( SCIPaddCons(heurdata->subscip, cons) );
                     SCIP_CALL( SCIPcaptureCons(heurdata->subscip, cons) );

                     assert( SCIPhashmapGetImage(heurdata->indicopymap, indicatorbinvar) != NULL );

                     if ( negated )
                     {
                        SCIP_CALL( SCIPgetNegatedVar(heurdata->subscip, indicatorcopy, &indicatorcopy) );
                     }

                     SCIP_CALL( SCIPchgVarType(heurdata->subscip, indicatorbinvar, SCIP_VARTYPE_CONTINUOUS, &feasible) );

                     SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, indislackvarpos, var) );
                     SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, indislackvarneg, var) );
                     SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );
                     SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );
                     SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &indislackvarpos) );
                     SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &indislackvarneg) );
                  }
                  else
                  {
                     if (!SCIPvarIsNegated(indicatorbinvar))
                        indicatorcopy = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->indicopymap, indicatorbinvar);
                     else
                     {
                        negatedvar = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->indicopymap, negatedvar);
                        SCIP_CALL( SCIPgetNegatedVar(heurdata->subscip, negatedvar, &indicatorcopy) );
                     }
                  }

                  cons = NULL;
                  SCIP_CALL( SCIPcreateConsIndicatorLinCons(heurdata->subscip, &cons, SCIPconsGetName(indicons), indicatorcopy,
                        SCIPgetLinearConsIndicator(indicons), SCIPgetSlackVarIndicator(indicons), SCIPconsIsInitial(indicons),
                        SCIPconsIsSeparated(indicons), SCIPconsIsEnforced(indicons), SCIPconsIsChecked(indicons),
                        SCIPconsIsPropagated(indicons), SCIPconsIsLocal(indicons), SCIPconsIsDynamic(indicons),
                        SCIPconsIsRemovable(indicons), SCIPconsIsStickingAtNode(indicons)) );
                  SCIP_CALL( SCIPaddCons(heurdata->subscip, cons) );

                  /* delete old indicator constraints so we can relax the indicator variables */
                  imagecons = (SCIP_CONS*) SCIPhashmapGetImage(heurdata->origsubscipConsMap, (void*)(currcons));
                  assert(imagecons != NULL);
                  SCIP_CALL( SCIPreleaseCons(heurdata->subscip, &imagecons) );
                  SCIP_CALL( SCIPhashmapRemove(heurdata->origsubscipConsMap, currcons) );
                  SCIP_CALL( SCIPhashmapInsert(heurdata->origsubscipConsMap, currcons, cons) );
                  SCIPconsAddUpgradeLocks(SCIPgetLinearConsIndicator(indicons), -1);
                  SCIP_CALL( SCIPdelCons(heurdata->subscip, indicons) );
               }

               SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, slackvarpos, var) );
               SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, slackvarneg, var) );
               SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );
               SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );

               SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, linindicons, slackvarpos, 1.0) );
               SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, linindicons, slackvarneg, -1.0) );
               SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &slackvarpos) );
               SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &slackvarneg) );
            }
         }
         continue;
      }

      if( heurdata->relaxindicators )
      {
         /* relax the old indicator variables*/
         for( k = 0; k < nvars; k++ )
         {
            if( SCIPhashmapGetImage(heurdata->indicators, vars[i]) == NULL )
               continue;

            tmpvar = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsciptosubscip, vars[k]);
            SCIP_CALL( SCIPchgVarType(heurdata->subscip, tmpvar, SCIP_VARTYPE_CONTINUOUS, &feasible) );
            SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, tmpvar, -SCIPinfinity(heurdata->subscip)) );
            SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, tmpvar,  SCIPinfinity(heurdata->subscip)) );
         }

         /* we map all slack variables of indicator constraints to their indicator variables */
         conshdlrindicator  = SCIPfindConshdlr(scip, "indicator");
         nconsindicator = SCIPconshdlrGetNConss(conshdlrindicator);

         /* delete old hashmaps and fill with the new indicators*/
         SCIP_CALL( releaseHashmapEntries(scip, heurdata->slacktoindivarsmap, TRUE) );
         SCIP_CALL( releaseHashmapEntries(scip, heurdata->indicators, FALSE) );
         SCIP_CALL( SCIPhashmapRemoveAll(heurdata->slacktoindivarsmap) );
         SCIP_CALL( SCIPhashmapRemoveAll(heurdata->indicators) );

         /* fill hashmaps with new values */
         for( k = 0; k < nconsindicator; k++ )
         {
            SCIP_CONS** indicatorconss = SCIPconshdlrGetConss(conshdlrindicator);
            SCIP_CONS* currcons;

            currcons = indicatorconss[k];
            assert(currcons != NULL);

            SCIP_CALL( SCIPcaptureVar(scip, SCIPgetBinaryVarIndicator(currcons)) );
            SCIP_CALL( SCIPcaptureCons(scip, currcons) );

            SCIP_CALL( SCIPhashmapInsert(heurdata->slacktoindivarsmap, SCIPgetSlackVarIndicator(currcons),
                  SCIPgetBinaryVarIndicator(currcons)) );
            SCIP_CALL( SCIPhashmapInsert(heurdata->indicators, SCIPgetBinaryVarIndicator(currcons), currcons) );
         }
      }

      /* in this case, we have a normal variable */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "relax_%s", SCIPvarGetName(var));
      cons = NULL;
      SCIP_CALL( SCIPcreateConsLinear( heurdata->subscip, &cons, consname, 0, NULL, NULL, 0.0, 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_pos0", SCIPvarGetName(var));
      SCIP_CALL( SCIPcreateVar( heurdata->subscip, &slackvarpos, varname, 0.0, SCIPinfinity(heurdata->subscip),
            heurdata->lambdaslack * 100, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL,NULL) );
      SCIP_CALL( SCIPaddVar(heurdata->subscip, slackvarpos) );

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_neg0", SCIPvarGetName(var));
      SCIP_CALL( SCIPcreateVar(heurdata->subscip, &slackvarneg, varname, 0.0, SCIPinfinity(heurdata->subscip),
            heurdata->lambdaslack * 100, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL,NULL) );
      SCIP_CALL( SCIPaddVar(heurdata->subscip, slackvarneg) );

      SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, var, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, slackvarpos, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, slackvarneg, -1.0) );

      SCIP_CALL( SCIPaddCons(heurdata->subscip, cons) );

      SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, slackvarpos, var) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, slackvarneg, var) );
      SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );
      SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );
      SCIP_CALL( SCIPhashmapInsert(heurdata->relaxcons, var, cons) );
      SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &slackvarpos) );
      SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &slackvarneg) );

      /* if the var is no indicator, relax it to a continuous variable */
      if( SCIPhashmapGetImage(heurdata->indicators, SCIPhashmapGetImage(heurdata->varsubsciptoscip, var)) == NULL )
      {
         SCIP_CALL( SCIPchgVarType(heurdata->subscip, var, SCIP_VARTYPE_CONTINUOUS, &feasible) );
         SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, var, -SCIPinfinity(heurdata->subscip)) );
         SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, var,  SCIPinfinity(heurdata->subscip)) );
      }
   }

   /* set up relaxation constraints for continous variables */
   if( heurdata->relaxcontvars )
   {
      for( i = 0; i < nvars; ++i )
      {
         var = SCIPvarGetTransVar(vars[i]);
         assert( var != NULL );

         if( ! SCIPvarIsActive(var) )
            continue;

         if( SCIPvarIsIntegral(var) )
            continue;

         if( SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), SCIPvarGetLbGlobal(var)) )
            continue;

         if( (SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), SCIPinfinity(scip))) && (SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), -SCIPinfinity(scip))) )
            continue;

         var = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsciptosubscip, var);
         assert( var != NULL );

         /* in this case, we have a normal variable */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "relax_ub_%s", SCIPvarGetName(var));
         cons = NULL;
         SCIP_CALL( SCIPcreateConsLinear( heurdata->subscip, &cons, consname, 0, NULL, NULL, -SCIPinfinity(heurdata->subscip), SCIPvarGetUbGlobal(var),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_pos2", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVar( heurdata->subscip, &slackvarpos, varname, 0.0, SCIPinfinity(heurdata->subscip),
               heurdata->lambdaslack, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL,NULL) );
         SCIP_CALL( SCIPaddVar(heurdata->subscip, slackvarpos) );

         SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, var, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, slackvarpos, -1.0) );

         SCIP_CALL( SCIPaddCons(heurdata->subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(heurdata->subscip, &cons) );
         SCIP_CALL( SCIPhashmapInsert(heurdata->slackvarubMap, var, slackvarpos) );
         SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, slackvarpos, var) );
         SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "relax_lb_%s", SCIPvarGetName(var));
         cons = NULL;
         SCIP_CALL( SCIPcreateConsLinear( heurdata->subscip, &cons, consname, 0, NULL, NULL, SCIPvarGetLbGlobal(var), SCIPinfinity(heurdata->subscip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );

         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "relax_%s_neg2", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVar( heurdata->subscip, &slackvarneg, varname, 0.0, SCIPinfinity(heurdata->subscip),
               heurdata->lambdaslack, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL,NULL) );
         SCIP_CALL( SCIPaddVar(heurdata->subscip, slackvarneg) );

         SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, var, 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, slackvarneg, 1.0) );

         SCIP_CALL( SCIPaddCons(heurdata->subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(heurdata->subscip, &cons) );
         SCIP_CALL( SCIPhashmapInsert(heurdata->slackvarlbMap, var, slackvarneg) );
         SCIP_CALL( SCIPhashmapInsert(heurdata->slack2var, slackvarneg, var) );
         SCIP_CALL( SCIPcaptureVar(heurdata->subscip, var) );

         SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, var, -SCIPinfinity(heurdata->subscip)) );
         SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, var,  SCIPinfinity(heurdata->subscip)) );
      }
   }

   /* if we have a solution add constraint that the next solution must not be worse than the current one */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "objbound");
   SCIP_CALL( SCIPcreateConsLinear( heurdata->subscip, &cons, consname, 0, NULL, NULL, -SCIPinfinity(scip),
         SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) );
   heurdata->objbound = cons;

   for( i = 0; i < nvars; ++i )
   {
      var = SCIPvarGetTransVar(vars[i]);
      assert( var != NULL );

      if( !SCIPvarIsActive(var) )
         continue;

      subvar = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsciptosubscip, var);
      assert( subvar != NULL );

      SCIP_CALL( SCIPaddCoefLinear(heurdata->subscip, cons, subvar, SCIPvarGetObj(var)) );

      SCIP_CALL( SCIPchgVarObj(heurdata->subscip, subvar, heurdata->lambdaobj * SCIPvarGetObj(subvar) ) );
   }

   SCIP_CALL( SCIPaddCons(heurdata->subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(heurdata->subscip, &cons) );

   /* do not need varsmap and conssmap anymore */
   SCIPhashmapFree(&conssmap);
   SCIPhashmapFree(&varsmap);

   /* enable SCIP output if needed */
   if( heurdata->heurverblevel > 3 )
   {
      SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "display/verblevel", 4) );
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "display/verblevel", 0) );
   }

   heurdata->nintegervars = j;

   return SCIP_OKAY;
}


/** free sub-SCIP data structure */
static
SCIP_RETCODE freeSubSCIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(heurdata->subscip != NULL);

   heurdata->nsubvars = 0;
   heurdata->nvars = 0;

   /* free sub-SCIP */
   SCIP_CALL( SCIPfree(&heurdata->subscip) );

   return SCIP_OKAY;
}


/** create a solution from the values of current nonlinear program */
static
SCIP_RETCODE createSolFromNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL**            sol                 /**< buffer to store solution value; if pointing to NULL a new solution is
                                                created, otherwise values in the given one are overwritten */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR**     subvars;
   SCIP_VAR*      subvar;
   int            i;
   int            nsubvars;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol  != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( *sol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );
   }

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP */
   assert(heurdata->nsubvars <= SCIPgetNOrigVars(heurdata->subscip));

   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, NULL, NULL, NULL, NULL) );

   /* set solution values */
   for( i = 0; i < nsubvars; ++i )
   {
      subvar = subvars[i];
      assert(subvar != NULL);

      subvar = SCIPvarGetTransVar(subvar);

      if( !SCIPvarIsActive(subvar) )
         continue;

      assert(SCIPvarGetNLPSol(subvar) != SCIP_INVALID);/*lint !e777*/
      SCIP_CALL( SCIPsetSolVal(scip, *sol, subvar, SCIPvarGetNLPSol(subvar)) );
   }

   return SCIP_OKAY;
}

#define BIG_VALUE 1E+10

/** method to fix the (relaxed) discrete variables */
static
SCIP_RETCODE fixDiscreteVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from;
                                              *   if NULL, then LP solution is used */
   SCIP_SOL**            transsol            /**< pointer to new created solution with fixed values as solution value */
   )
{
   SCIP_Real  fixval;
   SCIP_VAR*  var;
   SCIP_VAR*  subvar;
   SCIP_CONS* rcons;
   int        i;

   SCIP_CALL( SCIPcreateOrigSol(scip, transsol, NULL) );

   /* fix discrete variables */
   for( i = 0; i < heurdata->nintegervars; i++ )
   {
      var = heurdata->integervars[i];
      assert(var != NULL);

      var = SCIPvarGetTransVar(var);
      assert(var != NULL);
      subvar = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsciptosubscip, var);

      if( subvar == NULL )
         continue;

      if ( SCIPhashmapGetImage(heurdata->indicopymap, subvar) != NULL )
         subvar = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->indicopymap, subvar);

      /* get value of the variables, taking NULL as refpoint gives us the current LP solution,
       * otherwise we get our start point */
      fixval = SCIPgetSolVal(scip, refpoint, heurdata->integervars[i]);

      /* if we do not really have a startpoint, then we should take care that we do not fix variables to very large
       * values - thus, we set to 0.0 here and project on bounds below
       */
      if( REALABS(fixval) > BIG_VALUE && refpoint == NULL && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
         fixval = 0.0;

      /* round fractional variables to the nearest integer,
       * use exact integral value, if the variable is only integral within numerical tolerances
       */
      fixval = SCIPfloor(scip, fixval+0.5);

      /* adjust value to the global bounds of the corresponding SCIP variable */
      fixval = MAX(fixval, SCIPvarGetLbGlobal(heurdata->integervars[i]));  /*lint !e666*/
      fixval = MIN(fixval, SCIPvarGetUbGlobal(heurdata->integervars[i]));  /*lint !e666*/

      SCIP_CALL( SCIPsetSolVal(scip, *transsol, heurdata->integervars[i], fixval) );

      /* adjust the relaxation constraints to the new fixval */
      rcons = (SCIP_CONS*) SCIPhashmapGetImage(heurdata->relaxcons, subvar);

      fixval = MAX(fixval, SCIPvarGetLbGlobal(subvar));/*lint !e666*/
      fixval = MIN(fixval, SCIPvarGetUbGlobal(subvar));/*lint !e666*/
      if( rcons == NULL )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, fixval) );
         SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, fixval) );
         continue;
      }

      SCIP_CALL( SCIPchgLhsLinear(heurdata->subscip, rcons, fixval) );
      SCIP_CALL( SCIPchgRhsLinear(heurdata->subscip, rcons, fixval) );
   }

   return SCIP_OKAY;
}

/** method to free memory before leaving the heuristic or jumping up in the recursion */
static
SCIP_RETCODE freeMemory(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_SOL*             transsol,           /**< sol that has to be freed */
   SCIP_Real*            absranks,           /**< array of absolute rank values */
   SCIP_Real*            ranks,              /**< array of rank values */
   SCIP_VAR**            sortedvars,         /**< array of corresponding variables */
   SCIP_Bool             beforeswitching,    /**< did we call this method before or after switching variables? */
   SCIP_Bool             clearswitchedvars   /**< says if we should clear switchedvars or not */
   )
{
   SCIP_VAR**   subvars;
   SCIP_VAR*    subvar;
   SCIP_VAR*    var;
   SCIP_Real*   val;
   int          nsubvars;
   int          nsubbinvars;
   int          nsubintvars;
   int          i;

   if( clearswitchedvars )
   {
      /* free memory of the solution values in the hashmaps */
      for( i = 0; i < heurdata->nintegervars; i++ )
      {
         var = heurdata->integervars[i];

         if( SCIPhashmapGetImage(heurdata->slacktoindivarsmap, var) != NULL )
            var = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->slacktoindivarsmap, var);

         val = (SCIP_Real*)SCIPhashmapGetImage(heurdata->switchedvars, var);
         if( val != NULL )
         {
            SCIPfreeBlockMemoryArray(heurdata->subscip, &val, 1);
         }

         val = (SCIP_Real*)SCIPhashmapGetImage(heurdata->switchedvars2, var);
         if( val != NULL )
         {
            SCIPfreeBlockMemoryArray(heurdata->subscip, &val, 1);
         }
      }

      SCIP_CALL( SCIPhashmapRemoveAll(heurdata->switchedvars) );
      SCIP_CALL( SCIPhashmapRemoveAll(heurdata->switchedvars2) );
   }

   SCIPfreeBufferArrayNull( scip, &ranks );
   SCIPfreeBufferArrayNull( scip, &absranks );
   SCIPfreeBufferArrayNull( scip, &sortedvars );

   if( transsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &transsol) );
   }

   if( beforeswitching )
   {
      SCIP_CALL( SCIPfreeTransform(heurdata->subscip) );
   }

   /* undo fixing of discrete variables in sub-SCIP */
   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );

   /* set bounds of discrete variables to original values */
   for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
   {
      subvar = subvars[i];
      assert(SCIPvarGetProbindex(subvar) == i);

      var = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsubsciptoscip, subvar);

      if (SCIPhashmapGetImage(heurdata->indicopymapback, subvar) != NULL)
         var = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsubsciptoscip, SCIPhashmapGetImage(heurdata->indicopymapback, subvar));

      assert(var != NULL);

      SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, SCIPvarGetLbGlobal(var)) );
      SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, SCIPvarGetUbGlobal(var)) );
   }

   return SCIP_OKAY;
}

/** computes the ranks, saves them into an array and sorts the variables according to absolute ranks */
static
SCIP_RETCODE computeRanks(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data structure */
   SCIP_Real*            absranks,           /**< array of absolute rank values */
   SCIP_Real*            ranks,              /**< array of rank values */
   SCIP_VAR**            sortedvars          /**< array of corresponding variables */
   )
{
   SCIP_CONSHDLR*   conshdlrindicator;
   SCIP_CONS*       relaxcons;
   SCIP_CONS*       indicons;
   SCIP_CONS*       subcons;
   SCIP_CONS*       transcons;
   SCIP_VAR*        var;
   SCIP_Real*       dualvalue;
   int              nconsindicator;
   int              j;
   int              k;

   conshdlrindicator  = SCIPfindConshdlr(scip, "indicator");
   nconsindicator = SCIPconshdlrGetNConss(conshdlrindicator);

   /* Now we compute the rank of each variable */
   for( j = 0; j < heurdata->nintegervars; j++ )
   {
      sortedvars[j] = heurdata->integervars[j];
      ranks[j] = 0;
      absranks[j] = 0;

      if( sortedvars[j] == NULL )
         break;

      var = SCIPvarGetTransVar(sortedvars[j]);
      assert(var != NULL);

      /* globally fixed variables get rank 0 */
      if (SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)))
      {
         ranks[j] = 0;
         continue;
      }
      else
      {
         var = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsciptosubscip, var);
         assert(var != NULL);
         relaxcons = (SCIP_CONS*)SCIPhashmapGetImage(heurdata->relaxcons, (void*)(var));

         /* get ranks */
         if( relaxcons != NULL )
         {
            SCIP_CALL( SCIPgetTransformedCons(heurdata->subscip, relaxcons, &transcons) );
            dualvalue = (SCIP_Real*)SCIPhashmapGetImage(heurdata->dualvalues, (void*)transcons);

            if( dualvalue == NULL )
               dualvalue = (SCIP_Real*)SCIPhashmapGetImage(heurdata->dualvalues, (void*)(relaxcons));

            if( dualvalue == NULL )
               continue;

            assert(dualvalue != NULL);
            ranks[j] = (*dualvalue);

         }
         else /* if we have an indicator variable */
         {
            assert(ranks[j] == 0.0);

            if (SCIPhashmapGetImage(heurdata->relaxconsindi, (void*)(var)) != NULL)
            {
               subcons = (SCIP_CONS*)SCIPhashmapGetImage(heurdata->relaxconsindi, (void*)(var));

               dualvalue = (SCIP_Real*)SCIPhashmapGetImage(heurdata->dualvalues, (void*)(subcons));

               if( dualvalue == NULL )
                  continue;

               assert(dualvalue != NULL);

               ranks[j] = (*dualvalue);
            }

            /* compute the rank of the indicators, we take the highest dualvalue of an indicator constraint */
            for( k = 0; k < nconsindicator; k++ )
            {
               SCIP_CONS** indicatorconss = SCIPconshdlrGetConss(conshdlrindicator);
               SCIP_CONS* currcons;
               SCIP_VAR* indicatorbinvar;

               currcons = indicatorconss[k];
               assert(currcons != NULL);

               indicatorbinvar = SCIPgetBinaryVarIndicator(currcons);
               assert(indicatorbinvar != NULL);

               if( indicatorbinvar == (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsubsciptoscip, var)
                  || (SCIPvarIsNegated(indicatorbinvar) && indicatorbinvar == SCIPvarGetNegatedVar(var)) )
               {
                  indicons = currcons;
                  assert(indicons != NULL);

                  subcons = (SCIP_CONS*)SCIPhashmapGetImage(heurdata->origsubscipConsMap, (void*)(indicons));
                  assert(subcons != NULL);

                  subcons = SCIPgetLinearConsIndicator(subcons);
                  assert(subcons != NULL);

                  dualvalue = (SCIP_Real*)SCIPhashmapGetImage(heurdata->dualvalues, (void*)(subcons));

                  if( dualvalue == NULL )
                     continue;

                  assert(dualvalue != NULL);

                  if( REALABS(ranks[j]) < REALABS(*dualvalue) )
                     ranks[j] = (*dualvalue);
               }
            }
         }
      }

      /* take the absolute value of each rank */
      absranks[j] = REALABS(ranks[j]);
   }

   SCIPsortDownRealRealPtr(absranks, ranks, (void**)sortedvars, heurdata->nintegervars);

   return SCIP_OKAY;
}

/** compute maximal slack of a variable  */
static
SCIP_Real maximalslack(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_VAR* maxvar;
   SCIP_VAR* subvar;
   SCIP_SOL* bestsol;
   SCIP_Real maxslack;
   int i;
   int nsubvars;
   SCIP_Bool maxslackset;

   /* compute maximal slack */
   nsubvars = SCIPgetNOrigVars(heurdata->subscip);

   /* save information about maximal violation */
   maxvar   = NULL;
   maxslack = -SCIPinfinity(heurdata->subscip);
   maxslackset = FALSE;

   bestsol = SCIPgetBestSol(heurdata->subscip);

   /* search for variable with maximal slack */
   for( i = 0; i < nsubvars; i++ )
   {
      subvar = SCIPgetOrigVars(heurdata->subscip)[i];
      if( subvar == NULL)
         continue;

      /* if variable is slack */
      if( SCIPhashmapGetImage(heurdata->slack2var, subvar) != NULL )
      {
         if( heurdata->isnlp )
         {
            if( maxslack < SCIPvarGetNLPSol(subvar) )
            {
               maxslack = SCIPvarGetNLPSol(subvar);
               maxvar = subvar;
               maxslackset = TRUE;
            }
         }
         else
         {
            assert(bestsol != NULL);
            if( maxslack < SCIPgetSolVal(heurdata->subscip, bestsol, subvar) )
            {
               maxslack = SCIPgetSolVal(heurdata->subscip, bestsol, subvar);
               maxvar = subvar;
               maxslackset = TRUE;
            }
         }
      }
   }

   if( ! maxslackset )
   {
      maxslack = 0;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "could not find a variable with maximal slack!\n");
   }

   assert(maxslack >= 0);

   if( heurdata->heurverblevel > 0 && maxslackset )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "maximum slack: %f %s\n", maxslack, SCIPvarGetName(maxvar));
   }

   return maxslack;
}

/** method called after a solution is found which is feasible in the original problem, stores it and cleans up */
static
SCIP_RETCODE storeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data */
   SCIP_RESULT*          result,             /**< pointer to store result of: did not run, solution found,
                                                no solution found, or fixing is infeasible (cutoff) */
   SCIP_SOL*             transsol,           /**< solution to fix variables */
   SCIP_SOL*             bestsol             /**< solution we create a original scip solution from */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol = NULL;
   SCIP_Bool stored;
   SCIP_Real primalobj;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIP_CALL( createSolFromSubScipSol(scip, heur, &sol, bestsol) );

   /* if this happens, there was an ipopt error - stop the heuristic for there is no good starting point */
   if( heurdata->isnlp && SCIPgetNLPSolstat(heurdata->subscip) > SCIP_NLPSOLSTAT_FEASIBLE )
   {
      *result = SCIP_DIDNOTFIND;
      heurdata->solfound = TRUE;

      /* here we can be sure that we are in the nlp case */
      assert( heurdata->isnlp );
      SCIP_CALL( SCIPfreeSol(heurdata->subscip, &bestsol) );

      SCIP_CALL( freeMemory(scip, heurdata, transsol, NULL, NULL, NULL, TRUE, TRUE) );

      /* don't use the heuristic anymore if IPOPT doesn't give proper solution
       * (normally then this happens in most ipopt runs that may follow)       */
      SCIPheurSetFreq(heur, -1);

      SCIPdebugMsg(scip, "return10 : turn off heuristic, ipopt error\n");

      if( heurdata->heurverblevel > 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "turn off heuristic due to ipopt error");
      }

      return SCIP_OKAY;
   }

   primalobj = SCIPinfinity(scip);

   /* if there is a solution, try to add solution to storage and free it */
   if( sol != NULL )
   {
      primalobj = SCIPsolGetOrigObj(sol);

      if( heurdata->heurverblevel > 0 )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, FALSE, TRUE, &stored) );
      }
      else
      {
         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, FALSE, TRUE, &stored) );
      }
   }
   else
      stored = FALSE;

   if( stored && heurdata->heurverblevel > 1 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "accepted solution\n");

   if( heurdata->isnlp )
      SCIP_CALL( SCIPfreeSol(heurdata->subscip, &bestsol) );

   SCIP_CALL( freeMemory(scip, heurdata, transsol, NULL, NULL, NULL, TRUE, TRUE) );

   if( !stored )
   {
      *result = SCIP_DIDNOTFIND;
      heurdata->solfound = TRUE;

      if( heurdata->heurverblevel >= 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "return9 : found solution that was not stored, objective %f\n", primalobj);/*lint !e644*/
      }

      return SCIP_OKAY;
   }

   heurdata->prevInfeasible = FALSE;
   heurdata->solfound = TRUE;
   *result = SCIP_FOUNDSOL;

   if( heurdata->heurverblevel >= 1 )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "return 9 : found and stored new solution, objective %lf\n", primalobj);

   return SCIP_OKAY;
}

/** main procedure of the dualval heuristic */
SCIP_RETCODE SCIPapplyHeurDualval(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_RESULT*          result,             /**< pointer to store result of: did not run, solution found, no solution
                                              *   found, or fixing is infeasible (cutoff) */
   SCIP_SOL*             refpoint            /**< point to take fixation of discrete variables from; if NULL, then LP
                                              *   solution is used */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_NLROW*  nlrow;
   SCIP_SOL*    transsol;
   SCIP_SOL*    bestsol;
   SCIP_CONS**  subconss;
   SCIP_CONS*   rcons;
   SCIP_VAR**   subvars;
   SCIP_VAR**   sortedvars;
   SCIP_VAR*    var;
   SCIP_VAR*    subvar;
   SCIP_VAR*    v;
   SCIP_RETCODE retcode;
   SCIP_Real*   absranks;
   SCIP_Real*   ranks;
   SCIP_Real*   startpoint;
   SCIP_Real*   dualval;
   SCIP_Real*   lastval;
   SCIP_Real*   seclastval;
   SCIP_Real*   newval;
   SCIP_Real    bound;
   SCIP_Real    maxslack;
   SCIP_Real    objvalue;
   int          i;
   int          k;
   int          nsubvars;
   int          nsubbinvars;
   int          nsubintvars;
   int          nsubconss;
   int          maxequalranks;

   assert(scip != NULL);
   assert(heur != NULL);

   /* dio not run without nlp solver */
   if( SCIPgetNNlpis(scip) <= 0 )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't use the heuristic, if the gap is small so we don't expect to get better solutions than already found */
   if( SCIPgetGap(scip) * 100 < heurdata->mingap )
   {
      SCIPdebugMsg(scip, "return13 : gap is less than mingap\n");
      return SCIP_OKAY;
   }

   /* in the mode 'onlyleaves' don't run the heuristic if we are not in a leaf of the B&B tree */
   if( heurdata->onlyleaves && (SCIPgetNLPBranchCands(scip) != 0 || SCIPgetNPseudoBranchCands(scip) != 0) )
      return SCIP_OKAY;

   /* try to setup subscip if not tried before */
   if( heurdata->subscip == NULL && !heurdata->triedsetupsubscip )
   {
      SCIP_CALL( createSubSCIP(scip, heurdata) );
   }

   /* quit the recursion if we have found a solution */
   if( heurdata->solfound )
   {
      SCIPdebugMsg(scip, "return1 : already found solution \n");
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTRUN;

   /* not initialized */
   if( heurdata->subscip == NULL )
   {
      SCIPdebugMsg(scip, "return2 : subscip is NULL\n");
      return SCIP_OKAY;
   }

   assert(heurdata->nsubvars > 0);
   assert(heurdata->varsubsciptoscip != NULL);

   /* fix discrete variables in sub-SCIP */
   SCIP_CALL( fixDiscreteVars(scip, heurdata, refpoint, &transsol) );
   bound = SCIPgetUpperbound(scip);

   if( heurdata->onlycheaper && !SCIPisInfinity(scip, bound) )
   {
      SCIP_CALL( SCIPchgRhsLinear( heurdata->subscip, heurdata->objbound, bound) );
   }

   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "presolving/maxrounds", 1) );
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "propagating/maxroundsroot", 0) );

   SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 1LL) );
   SCIP_CALL( SCIPpresolve(heurdata->subscip) );

   if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_INFEASIBLE )
   {
      SCIPdebugMsg(scip, "return 4 : subscip is infeasible\n");

      *result = SCIP_DIDNOTFIND;
      heurdata->prevInfeasible = TRUE;
      SCIP_CALL( freeMemory(scip, heurdata, transsol, NULL, NULL, NULL, TRUE, TRUE) );

      return SCIP_OKAY;
   }

   /* If no NLP was constructed, then there were no nonlinearities after presolve.
    * So we increase the nodelimit to 1 and hope that SCIP will find some solution to this probably linear subproblem.
    */
   SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 0LL) );
   retcode = SCIPsolve(heurdata->subscip);
   heurdata->isnlp = TRUE;

   bestsol = NULL;

   /* we have no dualvalues, so give up */
   if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_OPTIMAL)
   {
      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( freeMemory(scip, heurdata, transsol, NULL, NULL, NULL, TRUE, TRUE) );

      return SCIP_OKAY;
   }

   if( ! SCIPisNLPConstructed(heurdata->subscip) && retcode == SCIP_OKAY )
   {
      SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 1LL) );
      SCIP_CALL( SCIPsolve(heurdata->subscip) );
      heurdata->isnlp = FALSE;
      bestsol = SCIPgetBestSol(heurdata->subscip);
   }

   if( heurdata->isnlp )
   {
      /* add non-combinatorial linear constraints from subscip into subNLP */
      SCIP_CALL( addLinearConstraintsToNlp(heurdata->subscip, FALSE, TRUE, heurdata) );

      SCIP_CALL( SCIPallocBufferArray(scip, &startpoint, SCIPgetNNLPVars(heurdata->subscip)) );

      /* set starting values (=refpoint, if not NULL; otherwise LP solution (or pseudo solution)) */
      for( i = 0; i < SCIPgetNNLPVars(heurdata->subscip); ++i )
      {
         SCIP_Real scalar = 1.0;
         SCIP_Real constant = 0.0;

         subvar = SCIPgetNLPVars(heurdata->subscip)[i];

         /* gets corresponding original variable */
         SCIP_CALL( SCIPvarGetOrigvarSum(&subvar, &scalar, &constant) );
         if( subvar == NULL )
         {
            startpoint[i] = constant;
            continue;
         }

         var = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsubsciptoscip, subvar);
         if( var == NULL || REALABS( SCIPgetSolVal(scip, refpoint, var) ) > 1.0e+12 )
         {
            SCIP_Real tmpmax;
            tmpmax = MAX( 0.0, SCIPvarGetLbGlobal(subvar) );/*lint !e666*/
            startpoint[i] = MIN( tmpmax, SCIPvarGetUbGlobal(subvar) );/*lint !e666*/
         }
         else
            /* scalar*subvar+constant corresponds to nlpvar[i], so nlpvar[i] gets value scalar*varval+constant */
            startpoint[i] = scalar * SCIPgetSolVal(scip, refpoint, var) + constant;
      }

      SCIP_CALL( SCIPsetNLPInitialGuess(heurdata->subscip, startpoint) );

      /* don't need startpoint array anymore */
      SCIPfreeBufferArray( scip, &startpoint );

      SCIP_CALL( SCIPsetNLPIntPar(heurdata->subscip, SCIP_NLPPAR_VERBLEVEL, heurdata->nlpverblevel) );

      SCIP_CALL( SCIPsolveNLP(heurdata->subscip) );
      assert(SCIPisNLPConstructed(heurdata->subscip));

      /* in this case there was an error in ipopt, we try to give another startpoint */
      if( SCIPgetNLPSolstat(heurdata->subscip) > SCIP_NLPSOLSTAT_FEASIBLE )
      {
         SCIP_CALL( SCIPsetNLPInitialGuess(heurdata->subscip, NULL) );
         SCIP_CALL( SCIPsolveNLP(heurdata->subscip) );
         assert(SCIPisNLPConstructed(heurdata->subscip));
      }

      nsubconss = SCIPgetNOrigConss(heurdata->subscip);
      subconss = SCIPgetOrigConss(heurdata->subscip);

      /* free memory of all entries and clear the hashmap before filling it */
      for( i = 0; i < nsubconss; i++ )
      {
         dualval = (SCIP_Real*)SCIPhashmapGetImage(heurdata->dualvalues, subconss[i]);
         SCIPfreeBlockMemoryArray(heurdata->subscip, &dualval, 1);
      }
      SCIP_CALL( SCIPhashmapRemoveAll(heurdata->dualvalues) );

      /* save the dualvalues from our nlp solution */
      for( i = 0; i < nsubconss; i++ )
      {
         SCIP_CONS* transcons;

         SCIP_CALL( SCIPgetTransformedCons(heurdata->subscip, subconss[i], &transcons) );

         if( transcons == NULL )
            continue;

         if( SCIPconsGetHdlr(transcons) != SCIPfindConshdlr(heurdata->subscip, "linear") )
            continue;

         nlrow = (SCIP_NLROW*)SCIPhashmapGetImage(heurdata->conss2nlrow, transcons);

         if (nlrow != NULL)
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(heurdata->subscip, &dualval, 1) ); /*lint !e506*/
            *dualval = SCIPnlrowGetDualsol(nlrow);
         }
         else
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(heurdata->subscip, &dualval, 1) ); /*lint !e506*/
            *dualval = 0;
         }

         SCIP_CALL( SCIPhashmapInsert(heurdata->dualvalues, subconss[i], dualval) );
      }

      bestsol = NULL;
      SCIP_CALL( createSolFromNLP(heurdata->subscip, heur, &bestsol) );
   }

   /* if we are infeasible, we can't do anything*/
   if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_INFEASIBLE )
   {
      SCIPdebugMsg(scip, "return4 : the subscip is infeasible\n");

      SCIP_CALL( freeMemory(scip, heurdata, transsol, NULL, NULL, NULL, TRUE, TRUE) );

      return SCIP_OKAY;
   }

   maxslack = maximalslack(scip, heurdata);
   SCIPdebugMsg(scip, "origObj: %f\n", SCIPgetSolOrigObj(heurdata->subscip, bestsol));
   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
   objvalue = 0.0;
   assert(bestsol != NULL);

   /* save information about maximal violation */
   for( i = 0; i < nsubvars; i++ )
   {
      subvar = SCIPgetOrigVars(heurdata->subscip)[i];

      if( SCIPhashmapGetImage(heurdata->slack2var, subvar) == NULL )
         objvalue += SCIPvarGetObj(subvar) * SCIPgetSolVal(heurdata->subscip, bestsol, subvar);
   }

   /* we stop the heuristic if it does not come "closer" to a feasible solution*/
   if( heurdata->forceimprovements )
   {
      if( SCIPisGE(scip, SCIPgetSolOrigObj(heurdata->subscip, bestsol) - objvalue, heurdata->prevobjective) && maxslack > 0 )
      {
         heurdata->nonimprovingRounds++;
         SCIPdebugMsg(scip, "nonimpr rounds %d prevobj %f \n", heurdata->nonimprovingRounds, heurdata->prevobjective);

         /* leave, if we have not improved some iterations*/
         if( heurdata->nonimprovingRounds > heurdata->maxcalls/8 )
         {
            *result = SCIP_DIDNOTFIND;

            if( heurdata->isnlp )
            {
               SCIP_CALL( SCIPfreeSol(heurdata->subscip, &bestsol) );
            }

            SCIP_CALL( freeMemory(scip, heurdata, transsol, NULL, NULL, NULL, TRUE, TRUE) );

            heurdata->solfound = TRUE;
            heurdata->switchdifferent = TRUE;

            SCIPdebugMsg(scip, "return11 : solution did not improve\n");

            return SCIP_OKAY;
         }
      }
   }

   heurdata->prevobjective = SCIPgetSolOrigObj(heurdata->subscip, bestsol) - objvalue;

   /* in this case we found a feasible solution, store it, clean up and stop the heuristic*/
   if( SCIPisFeasLE(heurdata->subscip, maxslack, 0.0) )
      return storeSolution(scip, heur, result, transsol, bestsol);

   SCIP_CALL( SCIPallocBufferArray(scip, &ranks, heurdata->nintegervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedvars, heurdata->nintegervars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &absranks, heurdata->nintegervars) );

   /* compute ranks and sort them in non-increasing order */
   SCIP_CALL( computeRanks(scip, heurdata, absranks, ranks, sortedvars) );

   /* print out the highest ranks */
   if( heurdata->heurverblevel > 1 )
   {
      k = heurdata->rankvalue;

      if( heurdata->nintegervars < heurdata->rankvalue )
         k = heurdata->nintegervars;

      for( i = 0; i < k; i++ )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%i. rank: %f name: %s\n", i, ranks[i], SCIPvarGetName(sortedvars[i]));
      }
   }

   /* free solution */
   if( heurdata->isnlp )
      SCIP_CALL( SCIPfreeSol(heurdata->subscip, &bestsol) );

   /* we don't allow more than a third of the variables to have the same rank */
   maxequalranks = MIN(heurdata->maxequalranks, heurdata->nintegervars/3);

   if( heurdata->maxequalranks >= 0 && SCIPisFeasEQ(heurdata->subscip, REALABS(ranks[0]), REALABS(ranks[maxequalranks])) )
   {
      *result = SCIP_DIDNOTFIND;

      SCIPdebugMsg(scip, "return12 : equal maxranks\n");

      SCIP_CALL( freeMemory(scip, heurdata, transsol, absranks, ranks, sortedvars, TRUE, TRUE ) );
      return SCIP_OKAY;
   }

   /* now we can start switching the variable values */
   SCIP_CALL( SCIPfreeTransform(heurdata->subscip) );

   /* set bounds of fixed discrete variables to original values so we can switch */
   for( k = 0; k < heurdata->nintegervars; ++k )
   {
      var = heurdata->integervars[k];
      if( var == NULL )
         break;

      var = SCIPvarGetTransVar(var);
      subvar = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->varsciptosubscip, var);

      rcons = (SCIP_CONS*) SCIPhashmapGetImage(heurdata->relaxcons, subvar);
      if( rcons != NULL )
         continue;

      assert(var != NULL);
      assert(subvar != NULL);

      if ( SCIPhashmapGetImage(heurdata->indicopymap, subvar) != NULL )
         subvar = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->indicopymap, subvar);

      SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, SCIPvarGetLbGlobal(heurdata->integervars[k])) );
      SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, SCIPvarGetUbGlobal(heurdata->integervars[k])) );
   }

   /* switch variable with maximum ranking if possible */
   for( i = 0; i < heurdata->nintegervars; i++ )
   {
      v = sortedvars[i];
      SCIP_CALL( SCIPallocBlockMemoryArray(heurdata->subscip, &newval, 1) ); /*lint !e506*/

      /* compute the new value of the variable */

      /* if we have an indicator constraint, we turn it off */
      if( SCIPhashmapGetImage(heurdata->slacktoindivarsmap, v) != NULL )
      {
         /* get the indicator var of this constraint */
         v = (SCIP_VAR*)SCIPhashmapGetImage(heurdata->slacktoindivarsmap, v);

         /* set the value to 0 */
         SCIP_CALL( SCIPsetSolVal(scip, transsol, v, 0.0) );
         if( heurdata->heurverblevel > 1 )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Setting value of %s%s to 0\n", SCIPvarIsNegated(v) ? "(negated) " : " ", SCIPvarGetName(v));

         *newval = 0.0;
         SCIP_CALL( SCIPhashmapInsert(heurdata->switchedvars, v, newval) );
      }
      else
      {
         if( ranks[i] > 0 )
         {
            if( SCIPvarIsBinary(v) && SCIPisEQ(scip, 1.0, SCIPgetSolVal(scip, transsol, v)) )
               continue;

            /* ignore fixed vars in input */
            if( SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(v), SCIPvarGetUbGlobal(v)) )
               continue;

            *newval = SCIPgetSolVal(scip, transsol, v) + 1;
         }
         else
         {
            if( SCIPvarIsBinary(v) && SCIPisEQ(scip, 0.0, SCIPgetSolVal(scip, transsol, v)) )
               continue;

            if( SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(v), SCIPvarGetUbGlobal(v)) )
               continue;

            *newval = SCIPgetSolVal(scip, transsol, v) - 1;
         }
      }
      lastval = (SCIP_Real*)SCIPhashmapGetImage(heurdata->switchedvars, v);
      seclastval = (SCIP_Real*)SCIPhashmapGetImage(heurdata->switchedvars2, v);

      /* we don't want to set a variable to a value it already had,or set a binary variable more than once */
      if( (lastval != NULL && (SCIPvarIsBinary(v) || SCIPisFeasEQ(scip, *lastval, *newval))) || (seclastval != NULL && SCIPisFeasEQ(scip, *seclastval, *newval)) )
      {
         SCIPfreeBlockMemoryArray(heurdata->subscip, &newval, 1);
         continue;
      }
      else /* update the switchedvars values, switchedvars2 is the second last and switchedvars the last value */
      {
         if( seclastval != NULL )
            SCIPfreeBlockMemoryArray(heurdata->subscip, &seclastval, 1);

         SCIP_CALL( SCIPhashmapRemove(heurdata->switchedvars2, v) );
         SCIP_CALL( SCIPhashmapInsert(heurdata->switchedvars2, v, lastval) );
         SCIP_CALL( SCIPhashmapRemove(heurdata->switchedvars, v) );
         SCIP_CALL( SCIPhashmapInsert(heurdata->switchedvars, v, newval) );

         if( heurdata->heurverblevel > 1 )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Setting value of %s from %f to %f\n", SCIPvarGetName(v), SCIPgetSolVal(scip, transsol, v), *newval);

         SCIP_CALL( SCIPsetSolVal(scip, transsol, v, *newval) );
      }

      /* if we have exceeded our iterations limit give up without any solution */
      if( heurdata->usedcalls >= heurdata->maxcalls )
      {
         SCIPdebugMsg(scip, "return5 : reached iteration limit\n");

         SCIP_CALL( freeMemory(scip, heurdata, transsol, absranks, ranks, sortedvars, FALSE, TRUE) );
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }

      heurdata->usedcalls++;

      if( heurdata->heurverblevel > 1 )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "----- Total Calls: %d\n", heurdata->usedcalls);

      /* recursive call of the heuristic */
      SCIP_CALL( SCIPapplyHeurDualval(scip, heur, result, transsol) );

      /* just to go up in the recursion */
      if( *result == SCIP_DIDNOTFIND || heurdata->solfound || heurdata->prevInfeasible )
      {
         SCIPdebugMsg(scip, "return6 : go up\n");

         /* here we only go up one step and try another switch (switch the same variables again is forbidden
          * since they are contained in switchedvars) */
         if( heurdata->switchdifferent )
         {
            heurdata->switchdifferent = FALSE;
            heurdata->solfound = FALSE;
            *result = SCIP_DIDNOTRUN;
            heurdata->nonimprovingRounds -= 2;
         }

         if( heurdata->prevInfeasible )
         {
            heurdata->prevInfeasible = FALSE;
            heurdata->solfound = FALSE;
            *result = SCIP_DIDNOTRUN;
            heurdata->nonimprovingRounds++;
         }

         SCIP_CALL( freeMemory(scip, heurdata, transsol, absranks, ranks, sortedvars, FALSE, FALSE) );
         return SCIP_OKAY;
      }
   }

   if( heurdata->subscip == NULL )
   {
      /* something horrible must have happened that we decided to give up completely on this heuristic */
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "return7 : subscip was set NULL\n");

      SCIP_CALL( freeMemory(scip, heurdata, transsol, absranks, ranks, sortedvars, FALSE, TRUE) );
      return SCIP_OKAY;
   }
   assert(!SCIPisTransformed(heurdata->subscip));

   SCIPdebugMsg(scip, "return8 : cannot switch any variable\n");

   SCIP_CALL( freeMemory(scip, heurdata, transsol, absranks, ranks, sortedvars, FALSE, TRUE) );

   *result = SCIP_DIDNOTFIND;
   return SCIP_OKAY;
}


/* Callback methods of primal heuristic */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeDualval)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);

   SCIPfreeBlockMemory(scip, &heurdata);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitDualval)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   /* skip setting up sub-SCIP if heuristic is disabled or we do not want to run the heuristic */
   if( SCIPheurGetFreq(heur) < 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->subscip == NULL);
   assert(!heurdata->triedsetupsubscip);

   /* create sub-SCIP for later use */
   SCIP_CALL( createSubSCIP(scip, heurdata) );

   /* creating sub-SCIP may fail if the solver interfaces did not copy into subscip */
   if( heurdata->subscip == NULL )
      return SCIP_OKAY;

   /* if the heuristic is called at the root node, we want to be called directly after the initial root LP solve */
   if( SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP | HEUR_TIMING);

   SCIP_CALL( SCIPhashmapCreate(&heurdata->dualvalues, SCIPblkmem(scip), 512) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitDualval)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_CONS** subconss;
   SCIP_Real* dualval;
   int i;
   int nsubconss;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &heurdata->integervars, heurdata->integervarssize);

   if( heurdata->subscip != NULL)
   {
      nsubconss = SCIPgetNOrigConss(heurdata->subscip);
      subconss = SCIPgetOrigConss(heurdata->subscip);

      /* free memory of all entries and clear the hashmap before filling it */
      for( i = 0; i < nsubconss; i++ )
      {
         dualval = (SCIP_Real*)SCIPhashmapGetImage(heurdata->dualvalues, subconss[i]);
         SCIPfreeBlockMemoryArrayNull(heurdata->subscip, &dualval, 1);
      }
      SCIP_CALL( SCIPhashmapRemoveAll(heurdata->dualvalues) );
      SCIPhashmapFree(&heurdata->dualvalues);

      if( heurdata->varsciptosubscip != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->varsciptosubscip, TRUE) );

         SCIPhashmapFree(&heurdata->varsciptosubscip);
      }
      if( heurdata->origsubscipConsMap != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->origsubscipConsMap, FALSE) );

         SCIPhashmapFree(&heurdata->origsubscipConsMap);
      }
      if( heurdata->relaxcons != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->relaxcons, FALSE) );

         SCIPhashmapFree(&heurdata->relaxcons);
      }
      if( heurdata->conss2nlrow != NULL )
      {
         SCIP_CALL( releaseHashmapNLPRows(heurdata->subscip, heurdata->conss2nlrow) );

         SCIPhashmapFree(&heurdata->conss2nlrow);
      }
      if( heurdata->slack2var != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->slack2var, TRUE) );

         SCIPhashmapFree(&heurdata->slack2var);
      }
      if( heurdata->indicopymap != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->indicopymap, TRUE) );

         SCIPhashmapFree(&heurdata->indicopymap);
      }
      if( heurdata->indicopymapback != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->indicopymapback, TRUE) );

         SCIPhashmapFree(&heurdata->indicopymapback);
      }
      if( heurdata->relaxconsindi != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->relaxconsindi, FALSE) );

         SCIPhashmapFree(&heurdata->relaxconsindi);
      }
      if( heurdata->slackvarlbMap != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->slackvarlbMap, TRUE) );

         SCIPhashmapFree(&heurdata->slackvarlbMap);
      }
      if( heurdata->slackvarubMap != NULL )
      {
         SCIP_CALL( releaseHashmapEntries(heurdata->subscip, heurdata->slackvarubMap, TRUE) );

         SCIPhashmapFree(&heurdata->slackvarubMap);
      }

      if( heurdata->subscip != NULL )
      {
         SCIP_CALL( freeSubSCIP(scip, heurdata) );
      }
   }

   if( heurdata->varsubsciptoscip != NULL )
   {
      SCIP_CALL( releaseHashmapEntries(scip, heurdata->varsubsciptoscip, TRUE) );

      SCIPhashmapFree(&heurdata->varsubsciptoscip);
   }
   if( heurdata->slacktoindivarsmap != NULL )
   {
      SCIP_CALL( releaseHashmapEntries(scip, heurdata->slacktoindivarsmap, TRUE) );

      SCIPhashmapFree(&heurdata->slacktoindivarsmap);
   }
   if( heurdata->indicators != NULL )
   {
      SCIP_CALL( releaseHashmapEntries(scip, heurdata->indicators, FALSE) );

      SCIPhashmapFree(&heurdata->indicators);
   }
   if( heurdata->switchedvars != NULL )
   {
      SCIPhashmapFree(&heurdata->switchedvars);
   }
   if( heurdata->switchedvars2 != NULL )
   {
      SCIPhashmapFree(&heurdata->switchedvars2);
   }

   /* reset some flags and counters */
   heurdata->triedsetupsubscip = FALSE;
   heurdata->usedcalls = 0;
   heurdata->solfound = FALSE;
   heurdata->prevInfeasible = FALSE;

   assert(heurdata->subscip == NULL);
   assert(heurdata->varsubsciptoscip == NULL);
   assert(heurdata->varsciptosubscip == NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolDualval)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   /* skip setting up sub-SCIP if heuristic is disabled or we do not want to run the heuristic */
   if( SCIPheurGetFreq(heur) < 0 )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* creating sub-SCIP may fail if the solver interfaces did not copy into subscip */
   if( heurdata->subscip == NULL )
      return SCIP_OKAY;

   /* if the heuristic is called at the root node, we want to be called directly after the initial root LP solve */
   if( SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP | HEUR_TIMING);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolDualval)
{
   assert(scip != NULL);
   assert(heur != NULL);

   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecDualval)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* obviously, we did not do anything yet */
   *result = SCIP_DIDNOTRUN;

   /* init data */
   heurdata->usedcalls = 0;
   heurdata->prevInfeasible = FALSE;
   heurdata->solfound = FALSE;
   heurdata->nonimprovingRounds = 0;
   heurdata->prevobjective = INT_MAX;

   SCIP_CALL( SCIPapplyHeurDualval(scip, heur, result, NULL) );

   /* SCIP does not like cutoff as return, so we say didnotfind, since we did not find a solution */
   if( *result == SCIP_CUTOFF )
      *result = SCIP_DIDNOTFIND;

   /* reset timing, if it was changed temporary (at the root node) */
   if( heurtiming != HEUR_TIMING )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/* primal heuristic specific interface methods */

/** creates the dualval primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurDualval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata = NULL;
   SCIP_HEUR* heur = NULL;

   /* create dualval primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* include primal heuristic */

   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecDualval, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeDualval) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitDualval) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitDualval) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolDualval) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolDualval) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/forceimprovements",
         "exit if objective doesn't improve",
         &heurdata->forceimprovements, TRUE, DEFAULT_FORCEIMPROVEMENTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/onlycheaper",
         "add constraint to ensure that discrete vars are improving",
         &heurdata->onlycheaper, TRUE, DEFAULT_ONLYCHEAPER, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/onlyleaves",
         "disable the heuristic if it was not called at a leaf of the B&B tree",
         &heurdata->onlyleaves, FALSE, DEFAULT_ONLYLEAVES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/relaxindicators",
         "relax the indicator variables by introducing continuous copies",
         &heurdata->relaxindicators, FALSE, DEFAULT_RELAXINDICATORS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/relaxcontvars",
         "relax the continous variables",
         &heurdata->relaxcontvars, FALSE, DEFAULT_RELAXCONTVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/heurverblevel",
         "verblevel of the heuristic, default is 0 to display nothing",
         &heurdata->heurverblevel, FALSE, DEFAULT_HEURVERBLEVEL, 0, 4, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nlpverblevel",
         "verblevel of the nlp solver, can be 0 or 1",
         &heurdata->nlpverblevel, FALSE, DEFAULT_NLPVERBLEVEL, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/rankvalue",
         "number of ranks that should be displayed when the heuristic is called",
         &heurdata->rankvalue, FALSE, DEFAULT_RANKVALUE, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxcalls",
         "maximal number of recursive calls of the heuristic (if dynamicdepth is off)",
         &heurdata->maxcalls, FALSE, DEFAULT_MAXCALLS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/dynamicdepth",
         "says if and how the recursion depth is computed at runtime",
         &heurdata->dynamicdepth, FALSE, DEFAULT_DYNAMICDEPTH, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxequalranks",
         "maximal number of variables that may have maximal rank, quit if there are more, turn off by setting -1",
         &heurdata->maxequalranks, FALSE, DEFAULT_MAXEQUALRANKS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/mingap",
         "minimal gap for which we still run the heuristic, if gap is less we return without doing anything",
         &heurdata->mingap, FALSE, DEFAULT_MINGAP, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lambdaslack",
         "value added to objective of slack variables, must not be zero",
         &heurdata->lambdaslack, FALSE, DEFAULT_LAMBDASLACK, 0.1, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/lambdaobj",
         "scaling factor for the objective function",
         &heurdata->lambdaobj, FALSE, DEFAULT_LAMBDAOBJ, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
