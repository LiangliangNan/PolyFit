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

/**@file   concsolver_scip.c
 * @ingroup PARALLEL
 * @brief  implementation of concurrent solver interface for SCIP
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/concsolver_scip.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/syncstore.h"
#include "scip/boundstore.h"
#include <string.h>


/* event handler for synchronization */

#define EVENTHDLR_NAME         "sync"
#define EVENTHDLR_DESC         "event handler for synchronization of concurrent scip sovlers"

/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   int             filterpos;
};

/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeSync)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}



/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitSync)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SYNCSTORE*  syncstore;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);

   if( eventhdlrdata->filterpos < 0 && SCIPsyncstoreIsInitialized(syncstore) )
   {
      /* notify SCIP that your event handler wants to react on synchronization events */
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SYNC, eventhdlr, NULL, &eventhdlrdata->filterpos) );
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitSync)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* notify SCIP that your event handler wants to drop the event type synchronization found */
   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SYNC, eventhdlr, NULL, eventhdlrdata->filterpos) );
      eventhdlrdata->filterpos = -1;
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSync)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);

   SCIP_CALL( SCIPsynchronize(scip) );

   return SCIP_OKAY;
}


/** includes event handler for synchronization found */
static
SCIP_RETCODE includeEventHdlrSync(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLR*     eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata) );
   eventhdlrdata->filterpos = -1;

   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSync, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSync) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSync) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitSync) );

   return SCIP_OKAY;
}

/** data for a concurrent solver type */
struct SCIP_ConcSolverTypeData
{
   SCIP_Bool             loademphasis;       /**< should emphasis settings be loaded whe ncreatig an instance of this concurrent solver */
   SCIP_PARAMEMPHASIS    emphasis;           /**< parameter emphasis that will be loaded if loademphasis is true */
};

/** data for a concurrent solver */
struct SCIP_ConcSolverData
{
   SCIP*                 solverscip;         /**< the concurrent solvers private scip datastructure */
   SCIP_VAR**            vars;               /**< array of variables in the order of the main SCIP's variable array */
   int                   nvars;              /**< number of variables in the above arrays */
};

/** Disable dual reductions that might cut off optimal solutions. Although they keep at least
 *  one optimal solution intact, communicating these bounds may cut off all optimal solutions,
 *  if different optimal solutions were kept in different concurrent solvers. */
static
SCIP_RETCODE disableConflictingDualReductions(
   SCIP*                 scip                /**< SCIP datastructure */
   )
{
   SCIP_Bool commvarbnds;

   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/commvarbnds", &commvarbnds) );

   if( !commvarbnds )
      return SCIP_OKAY;

   SCIP_CALL( SCIPsetBoolParam(scip, "misc/allowdualreds", FALSE) );
   return SCIP_OKAY;
}

/** sets the child selection rule based on the index of the concurrent solver */
static
SCIP_RETCODE setChildSelRule(
   SCIP_CONCSOLVER*      concsolver          /**< the concurrent solver */
   )
{
   SCIP_CONCSOLVERDATA*  data;
   static char childsel[] = { 'h', 'i', 'p', 'r', 'l', 'd', 'u' };

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   SCIP_CALL( SCIPsetCharParam(data->solverscip, "nodeselection/childsel", childsel[SCIPconcsolverGetIdx(concsolver) % 7]) );

   return SCIP_OKAY;
}

/** initialize the concurrent SCIP solver, i.e. setup the copy of the problem and the
 *  mapping of the variables */
static
SCIP_RETCODE initConcsolver(
   SCIP*                 scip,               /**< the main SCIP instance */
   SCIP_CONCSOLVER*      concsolver          /**< the concurrent solver to set up */
   )
{
   int                 i;
   SCIP_VAR**          vars;
   SCIP_Bool           valid;
   SCIP_HASHMAP*       varmapfw;
   SCIP_CONCSOLVERDATA* data;
   int* varperm;

   assert(scip != NULL);
   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   data->nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* create the concurrent solver's SCIP instance and set up the problem */
   SCIP_CALL( SCIPcreate(&data->solverscip) );
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(data->solverscip), data->nvars) );
   SCIP_CALL( SCIPcopy(scip, data->solverscip, varmapfw, NULL, SCIPconcsolverGetName(concsolver), TRUE, FALSE, FALSE, &valid) );
   assert(valid);

   /* allocate memory for the arrays to store the variable mapping */
   SCIP_CALL( SCIPallocBlockMemoryArray(data->solverscip, &data->vars, data->nvars) );
   SCIP_CALL( SCIPallocBufferArray(data->solverscip, &varperm, data->nvars) );

   /* set up the arrays for the variable mapping */
   for( i = 0; i < data->nvars; i++ )
   {
      SCIP_VAR* var;
      var = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);
      varperm[SCIPvarGetIndex(var)] = i;
      data->vars[i] = var;
   }

   /* create the concurrent data structure for the concurrent solver's SCIP */
   /* this assert fails on check/instances/Orbitope/packorb_1-FullIns_3.cip
    * assert(SCIPgetNOrigVars(data->solverscip) == data->nvars);
    * also fails on check/instances/Orbitope/partorb_1-FullIns_3.cip
    * TODO: test if this leads to any problems
    */
   SCIP_CALL( SCIPcreateConcurrent(data->solverscip, concsolver, varperm) );
   SCIPfreeBufferArray(data->solverscip, &varperm);

   /* free the hashmap */
   SCIPhashmapFree(&varmapfw);

   return SCIP_OKAY;
}

/* creates an instance of a concurrent SCIP solver */
static
SCIP_DECL_CONCSOLVERCREATEINST(concsolverScipCreateInstance)
{
   SCIP_CONCSOLVERDATA*     data;
   SCIP_CONCSOLVERTYPEDATA* typedata;
   char*                    prefix;
   char                     filename[SCIP_MAXSTRLEN];
   SCIP_Bool                changechildsel;

   assert(scip != NULL);
   assert(concsolvertype != NULL);
   assert(concsolver != NULL);

   typedata = SCIPconcsolverTypeGetData(concsolvertype);

   SCIP_ALLOC( BMSallocMemory(&data) );
   SCIPconcsolverSetData(concsolver, data);

   SCIP_CALL( initConcsolver(scip, concsolver) );

   /* check if emphasis setting should be loaded */
   if( typedata->loademphasis )
   {
      SCIP_PARAM** params;
      SCIP_PARAM** fixedparams;
      int          nparams;
      int          nfixedparams;
      int          i;

      params = SCIPgetParams(data->solverscip);
      nparams = SCIPgetNParams(data->solverscip);
      SCIP_CALL( SCIPallocBufferArray(data->solverscip, &fixedparams, nparams) );
      nfixedparams = 0;

      /* fix certain parameters before loading emphasis to avoid setting them to default values */
      for( i = 0; i < nparams; ++i )
      {
         const char* paramname;

         paramname = SCIPparamGetName(params[i]);

         if( strncmp(paramname, "limits/", 7) == 0 ||
             strncmp(paramname, "numerics/", 9) == 0 ||
             strncmp(paramname, "memory/", 7) == 0 ||
             strncmp(paramname, "concurrent/sync/", 16) == 0 ||
             strncmp(paramname, "heuristics/sync/", 16) == 0 ||
             strncmp(paramname, "propagating/sync/", 17) == 0 )
         {
            fixedparams[nfixedparams++] = params[i];
            SCIP_CALL( SCIPfixParam(data->solverscip, paramname) );
         }
      }

      SCIP_CALL( SCIPsetEmphasis(data->solverscip, typedata->emphasis, TRUE) );

      for( i = 0; i < nfixedparams; ++i )
         SCIP_CALL( SCIPunfixParam(data->solverscip, SCIPparamGetName(fixedparams[i])) );

      SCIPfreeBufferArray(data->solverscip, &fixedparams);
   }

   /* load settings file if it exists */
   SCIP_CALL( SCIPgetStringParam(scip, "concurrent/paramsetprefix", &prefix) );
   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s%s.set", prefix, SCIPconcsolverGetName(concsolver));

   if( SCIPfileExists(filename) )
   {
      /* load settings file and print info message */
      SCIPinfoMessage(scip, NULL, "reading parameter file <%s> for concurrent solver <%s>\n", filename, SCIPconcsolverGetName(concsolver));
      SCIP_CALL( SCIPreadParams(data->solverscip, filename) );
   }
   else
   {
      /* print message about missing setting files only in verblevel full */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "skipping non existent parameter file <%s> for concurrent solver <%s>\n",
                      filename, SCIPconcsolverGetName(concsolver));
   }

   /* include eventhandler for synchronization */
   SCIP_CALL( includeEventHdlrSync(data->solverscip) );

   /* disable output for subscip */
   SCIP_CALL( SCIPsetIntParam(data->solverscip, "display/verblevel", 0) );

   /* use wall clock time in subscips */
   SCIP_CALL( SCIPsetIntParam(data->solverscip, "timing/clocktype", (int)SCIP_CLOCKTYPE_WALL) );

   /* don't catch ctrlc since already caught in main scip */
   SCIP_CALL( SCIPsetBoolParam(data->solverscip, "misc/catchctrlc", FALSE) );

   /* one solver can do all dual reductions and share them with the other solvers */
   if( SCIPconcsolverGetIdx(concsolver) != 0 )
   {
      SCIP_CALL( disableConflictingDualReductions(data->solverscip) );
   }

   /* set different child selection rules if corresponding parameter is TRUE */
   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/changechildsel", &changechildsel) );
   if( changechildsel )
   {
      SCIP_CALL( setChildSelRule(concsolver) );
   }

   return SCIP_OKAY;
}

/** destroys an instance of a concurrent SCIP solver */
static
SCIP_DECL_CONCSOLVERDESTROYINST(concsolverScipDestroyInstance)
{
   SCIP_CONCSOLVERDATA* data;

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);
   assert(data->solverscip != NULL);

   /* free the array with the variable mapping */
   SCIPfreeBlockMemoryArray(data->solverscip, &data->vars, data->nvars);

   /* free subscip */
   SCIP_CALL( SCIPfree(&data->solverscip) );
   BMSfreeMemory(&data);
   SCIPconcsolverSetData(concsolver, NULL);

   return SCIP_OKAY;
}

/** frees the data of a concurrent solver type */
static
SCIP_DECL_CONCSOLVERTYPEFREEDATA(concsolverTypeScipFreeData)
{
   BMSfreeMemory(data);
}

/** initializes the random and permutation seeds with the given one
 *  and enables permutation of constraints and variables
 */
static
SCIP_DECL_CONCSOLVERINITSEEDS(concsolverScipInitSeeds)
{
   SCIP_CONCSOLVERDATA* data;

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   SCIPinfoMessage(data->solverscip, NULL, "initializing seeds to %d in concurrent solver '%s'\n", (int) seed, SCIPconcsolverGetName(concsolver));

   SCIP_CALL( SCIPsetIntParam(data->solverscip, "randomization/randomseedshift", (int) seed) );
   SCIP_CALL( SCIPsetIntParam(data->solverscip, "randomization/permutationseed", (int) seed) );
   SCIP_CALL( SCIPsetBoolParam(data->solverscip, "randomization/permutevars", TRUE) );
   SCIP_CALL( SCIPsetBoolParam(data->solverscip, "randomization/permuteconss", TRUE) );

   return SCIP_OKAY;
}

/** installs the solving status of this concurrent solver and the solving statistics
 *  into the given SCIP instance
 */
static
SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA(concsolverGetSolvingData)
{
   SCIP_CONCSOLVERDATA* data;
   SCIP_VAR** vars;
   int nvars;
   int nsols;
   SCIP_SOL** sols;
   SCIP_Real* solvals;
   SCIP_HEUR* heur;
   int i;

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);
   assert(data->solverscip != NULL);

   assert(scip != NULL);
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   nsols = SCIPgetNSols(data->solverscip);
   sols = SCIPgetSols(data->solverscip);

   assert(nvars == data->nvars);

   /* allocate buffer array used for translating the solution to the given SCIP */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );

   /* add the solutions to the given SCIP */
   for( i = 0; i < nsols; ++i )
   {
      SCIP_SOL* sol;
      SCIP_Bool stored;
      SCIP_CALL( SCIPgetSolVals(data->solverscip, sols[i], nvars, data->vars, solvals) );

      heur = SCIPsolGetHeur(sols[i]);

      if( heur != NULL )
         heur = SCIPfindHeur(scip, SCIPheurGetName(heur));

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals) );

      SCIP_CALL( SCIPcopySolStats(sols[i], sol) );

      SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );
   }

   /* free the buffer array */
   SCIPfreeBufferArray(scip, &solvals);

   /* copy solving statistics and status from the solver SCIP to the given SCIP */
   SCIP_CALL( SCIPcopyConcurrentSolvingStats(data->solverscip, scip) );

   return SCIP_OKAY;
}

/** start solving the problem until the solving reaches a limit, gets interrupted, or
 *  just finished successfully
 */
static
SCIP_DECL_CONCSOLVEREXEC(concsolverScipExec)
{
   SCIP_CONCSOLVERDATA* data;
   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   /* print info message that solving has started */
   SCIPinfoMessage(data->solverscip, NULL, "starting solve in concurrent solver '%s'\n", SCIPconcsolverGetName(concsolver));

   /* solve */
   SCIP_CALL( SCIPsolve(data->solverscip) );

   /* print info message with status */
   SCIPinfoMessage(data->solverscip, NULL, "concurrent solver '%s' stopped with status ", SCIPconcsolverGetName(concsolver));
   SCIP_CALL( SCIPprintStatus(data->solverscip, NULL) );
   SCIPinfoMessage(data->solverscip, NULL, "\n");

   /* set solving statistics */
   *solvingtime = SCIPgetSolvingTime(data->solverscip);
   *nlpiterations = SCIPgetNLPIterations(data->solverscip);
   *nnodes = SCIPgetNNodes(data->solverscip);

   return SCIP_OKAY;
}

/** stops the concurrent solver as soon as possible */
static
SCIP_DECL_CONCSOLVERSTOP(concsolverScipStop)
{
   SCIP_CONCSOLVERDATA* data;
   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   SCIP_CALL( SCIPinterruptSolve(data->solverscip) );

   return SCIP_OKAY;
}

/** writes new solutions and global boundchanges to the iven synchronization data */
static
SCIP_DECL_CONCSOLVERSYNCWRITE(concsolverScipSyncWrite)
{
   int                    i;
   int                    nsols;
   SCIP_SOL**             sols;
   SCIP_CONCSOLVERDATA*   data;
   SCIP_BOUNDSTORE*       boundstore;
   int                    concsolverid;
   SCIP_STATUS            solverstatus;


   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);
   concsolverid = SCIPconcsolverGetIdx(concsolver);
   solverstatus = SCIPgetStatus(data->solverscip);

   SCIPsyncdataSetStatus(syncdata, solverstatus, concsolverid);
   SCIPsyncdataSetLowerbound(syncdata, SCIPgetDualbound(data->solverscip));
   SCIPsyncdataSetUpperbound(syncdata, SCIPgetPrimalbound(data->solverscip));

   *nsolsshared = 0;

   if( SCIPsyncdataGetStatus(syncdata) != SCIP_STATUS_UNKNOWN )
      return SCIP_OKAY;

   SCIPdebugMessage("syncing in concurrent solver %s\n", SCIPconcsolverGetName(concsolver));

   /* consider at most maxcandsols many solutions, and since
    * the solution array is sorted, we will cosider the best
    * solutions
    */
   nsols = SCIPgetNSols(data->solverscip);
   nsols = MIN(nsols, maxcandsols);
   sols = SCIPgetSols(data->solverscip);

   for( i = 0; i < nsols; ++i )
   {
      if( SCIPIsConcurrentSolNew(data->solverscip, sols[i]) )
      {
         SCIP_Real solobj;
         SCIP_Real* solvals;

         solobj = SCIPgetSolOrigObj(data->solverscip, sols[i]);

         SCIPdebugMessage("adding sol to spi in concurrent solver %s\n", SCIPconcsolverGetName(concsolver));
         SCIPsyncdataGetSolutionBuffer(syncstore, syncdata, solobj, concsolverid, &solvals);

         /* if syncstore has no place for this solution we can stop since the next solution will have
          * a worse objective value and thus won't be accepted either
          */
         if( solvals == NULL )
            break;

         ++(*nsolsshared);
         SCIP_CALL( SCIPgetSolVals(data->solverscip, sols[i], data->nvars, data->vars, solvals) );

         /* if we have added the maximum number of solutions we can also stop */
         if( *nsolsshared == maxsharedsols )
            break;
      }
   }

   boundstore = SCIPgetConcurrentGlobalBoundChanges(data->solverscip);

   if( boundstore != NULL )
      SCIP_CALL( SCIPsyncdataAddBoundChanges(syncstore, syncdata, boundstore) );

   SCIPsyncdataAddMemTotal(syncdata, SCIPgetMemTotal(data->solverscip));

   return SCIP_OKAY;
}

/** reads the solutions and bounds from the given synchronization data */
static
SCIP_DECL_CONCSOLVERSYNCREAD(concsolverScipSyncRead)
{  /*lint --e{715}*/
   int                    i;
   int                    nsols;
   SCIP_Real**            solvals;
   SCIP_CONCSOLVERDATA*   data;
   SCIP_BOUNDSTORE*       boundstore;
   int*                   concsolverids;
   int                    concsolverid;
   int                    nbndchgs;

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   concsolverid = SCIPconcsolverGetIdx(concsolver);

   /* get solutions from synchronization data */
   SCIPsyncdataGetSolutions(syncdata, &solvals, &concsolverids, &nsols);
   *nsolsrecvd = 0;

   for( i = 0; i < nsols; ++i )
   {
      SCIP_SOL* newsol;

      /* do not add own solutions */
      if( concsolverids[i] == concsolverid )
         continue;

      /* solution is from other solver so translate to this solvers variable space
       * and add it to the SCIP
       */
      ++(*nsolsrecvd);
      SCIP_CALL( SCIPcreateOrigSol(data->solverscip, &newsol, NULL) );

      SCIP_CALL( SCIPsetSolVals(data->solverscip, newsol, data->nvars, data->vars, solvals[i]) );
      SCIPdebugMessage("adding solution in concurrent solver %s\n", SCIPconcsolverGetName(concsolver));
      SCIP_CALL( SCIPaddConcurrentSol(data->solverscip, newsol) );
   }

   /* get bound changes from the synchronization data and add it to this
    * concurrent solvers SCIP
    */
   *ntighterbnds = 0;
   *ntighterintbnds = 0;
   boundstore = SCIPsyncdataGetBoundChgs(syncdata);
   nbndchgs = SCIPboundstoreGetNChgs(boundstore);

   for( i = 0; i < nbndchgs; ++i )
   {
      SCIP_VAR*   var;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Real newbound;

      var = data->vars[SCIPboundstoreGetChgVaridx(boundstore, i)];
      boundtype = SCIPboundstoreGetChgType(boundstore, i);
      newbound = SCIPboundstoreGetChgVal(boundstore, i);

      SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );

      /* cannot change bounds of multi-aggregated variables so dont pass this bound-change to the propagator */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
         return SCIP_OKAY;

      /* if bound is not better than also don't pass this bound to the propagator and
       * don't waste memory for storing this boundchange
       */
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisGE(data->solverscip, SCIPvarGetLbGlobal(var), newbound) )
         return SCIP_OKAY;

      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisLE(data->solverscip, SCIPvarGetUbGlobal(var), newbound) )
         return SCIP_OKAY;

      /* bound is better so incremented counters for statistics and pass it to the sync propagator */
      ++(*ntighterbnds);

      if( SCIPvarGetType(var) <= SCIP_VARTYPE_INTEGER )
         ++(*ntighterintbnds);

      SCIP_CALL( SCIPaddConcurrentBndchg(data->solverscip, var, newbound, boundtype) );
   }

   return SCIP_OKAY;
}


/** creates the concurrent SCIP solver plugins and includes them in SCIP */
SCIP_RETCODE SCIPincludeConcurrentScipSolvers(
   SCIP*                 scip                /**< SCIP datastructure */
   )
{
   SCIP_CONCSOLVERTYPEDATA* data;
   assert(scip != NULL);

   /* include concurrent solvers for SCIP for all emphasis settings and without an emphasis setting.
    * For the SCIP without an emphasis setting we set the default preferred priority to 1 and for the other types to 0
    * so that the default concurent solve will use multiple SCIP's using settings as specified by the user in the main SCIP
    */
   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = FALSE;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip", 1.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_DEFAULT;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-default", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_CPSOLVER;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-cpsolver", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_EASYCIP;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-easycip", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_FEASIBILITY;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-feas", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_HARDLP;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-hardlp", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_OPTIMALITY;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-opti", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_COUNTER;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-counter", 0.0,  concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
                                        concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
                                        concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   return SCIP_OKAY;
}
