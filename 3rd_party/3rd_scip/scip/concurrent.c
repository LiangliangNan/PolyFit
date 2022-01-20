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

/**@file   concurrent.c
 * @ingroup PARALLEL
 * @brief  helper functions for concurrent SCIP solvers
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/concurrent.h"
#include "scip/struct_concurrent.h"
#include "scip/concsolver.h"
#include "scip/event.h"
#include "scip/struct_scip.h"
#include "scip/stat.h"
#include "scip/struct_set.h"
#include "scip/struct_primal.h"
#include "scip/struct_stat.h"
#include "scip/struct_sol.h"
#include "scip/struct_prop.h"
#include "scip/struct_heur.h"
#include "scip/struct_sepa.h"
#include "scip/struct_presol.h"
#include "scip/prob.h"
#include "scip/prop_sync.h"
#include "scip/heur_sync.h"
#include "scip/event_globalbnd.h"
#include "scip/scip.h"
#include "scip/syncstore.h"
#include "scip/set.h"
#include "tpi/tpi.h"

/** create concurrent data */
SCIP_RETCODE SCIPcreateConcurrent(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver of given SCIP instance */
   int*                  varperm             /**< permutation of variables for communication */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(concsolver != NULL);
   assert(varperm != NULL);
   assert(scip->concurrent == NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &scip->concurrent) );

   nvars = SCIPgetNOrigVars(scip);
   scip->concurrent->varperm = NULL;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &scip->concurrent->varperm, varperm, nvars) );

   scip->concurrent->concsolver = concsolver;
   scip->concurrent->mainscip = scip;
   scip->concurrent->solidx = 0;
   scip->stat->subscipdepth = 0;

   if( scip->set->parallel_mode == (int) SCIP_PARA_DETERMINISTIC )
   {
      scip->concurrent->dettime = 0.0;
      scip->concurrent->wallclock = NULL;
   }
   else
   {
      SCIP_CALL( SCIPcreateWallClock(scip, &scip->concurrent->wallclock) );
      SCIP_CALL( SCIPstartClock(scip, scip->concurrent->wallclock) );
   }

   assert(SCIPfindHeur(scip, "sync") == NULL);

   SCIP_CALL( SCIPincludeHeurSync(scip) );
   scip->concurrent->heursync = SCIPfindHeur(scip, "sync");

   assert(SCIPfindProp(scip, "sync") == NULL);

   SCIP_CALL( SCIPincludePropSync(scip) );
   scip->concurrent->propsync = SCIPfindProp(scip, "sync");

   scip->concurrent->eventglobalbnd = NULL;
   assert(SCIPfindEventhdlr(scip, "globalbnd") == NULL);

   if( scip->set->concurrent_commvarbnds )
   {
      SCIP_CALL( SCIPincludeEventHdlrGlobalbnd(scip) );
      scip->concurrent->eventglobalbnd = SCIPfindEventhdlr(scip, "globalbnd");
   }

   return SCIP_OKAY;
}

/** get number of initialized concurrent solvers */
int SCIPgetNConcurrentSolvers(
   SCIP*                 scip                /**< SCIP datastructure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nconcsolvers;
}

/** gets the initialized concurrent solvers */
SCIP_CONCSOLVER** SCIPgetConcurrentSolvers(
   SCIP*                 scip                /**< SCIP datastructure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->concsolvers;
}

/** adds an initialized concurrent solver */
SCIP_RETCODE SCIPaddConcurrentSolver(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver of given SCIP instance */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPsetIncludeConcsolver(scip->set, concsolver) );

   return SCIP_OKAY;
}

/** frees concurrent data */
SCIP_RETCODE SCIPfreeConcurrent(
   SCIP*                 scip                /**< SCIP datastructure */
   )
{
   assert(scip != NULL);

   if( scip->concurrent == NULL )
      return SCIP_OKAY;

   assert(scip->concurrent->varperm != NULL);

   /* check if we are the SCIP that is responsible for freeing this concurent struct
    * or just a subscip */
   if( scip->concurrent->mainscip != scip )
   {
      /* we are just a subscip, so don't free the concurrent structure and add the
       * deterministic time that was counted in the subscip but not yet added to the main SCIP */
      scip->concurrent->mainscip->stat->detertimecnt += scip->stat->detertimecnt;
      scip->stat->detertimecnt = 0;
      scip->concurrent = NULL;
   }
   else
   {
      /* we are in the main SCIP so free the concurrent structure */
      if( scip->concurrent->wallclock != NULL )
      {
         SCIP_CALL( SCIPfreeClock(scip, &scip->concurrent->wallclock) );
      }

      SCIPfreeBlockMemoryArray(scip, &scip->concurrent->varperm, SCIPgetNOrigVars(scip));

      SCIPfreeBlockMemory(scip, &scip->concurrent);
   }

   return SCIP_OKAY;
}

/** increments the time counter for synchronization */
SCIP_RETCODE SCIPincrementConcurrentTime(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real             val                 /**< value by which the time counter for synchronization is incremented */
   )
{
   SCIP_Real           syncfreq;
   SCIP*               mainscip;
   SCIP_CLOCK*         wallclock;

   assert(scip != NULL);

   if( scip->concurrent == NULL )
      return SCIP_OKAY;

   syncfreq = SCIPconcsolverGetSyncFreq(scip->concurrent->concsolver);
   wallclock = scip->concurrent->wallclock;
   mainscip = scip->concurrent->mainscip;

   if( wallclock == NULL )
   {
      scip->concurrent->dettime += val;

      if( scip->concurrent->dettime >= syncfreq  )
      {
         SCIP_EVENT* event;
         SCIPconcsolverSetTimeSinceLastSync(scip->concurrent->concsolver, scip->concurrent->dettime);
         scip->concurrent->dettime = 0.0;
         SCIP_CALL( SCIPeventCreateSync(&event, SCIPblkmem(mainscip)) );
         SCIP_CALL( SCIPeventqueueAdd(mainscip->eventqueue, SCIPblkmem(mainscip), mainscip->set,
                                      NULL, NULL, NULL, mainscip->eventfilter, &event) );
      }
   }
   else
   {
      SCIP_Real timesincelastsync;
      timesincelastsync = SCIPgetClockTime(mainscip, wallclock);

      if( timesincelastsync >= syncfreq )
      {
         SCIP_EVENT* event;
         SCIPconcsolverSetTimeSinceLastSync(scip->concurrent->concsolver, timesincelastsync);

         SCIP_CALL( SCIPeventCreateSync(&event, SCIPblkmem(mainscip)) );
         SCIP_CALL( SCIPeventqueueAdd(mainscip->eventqueue, SCIPblkmem(mainscip), mainscip->set,
                                      NULL, NULL, NULL, mainscip->eventfilter, &event) );

         SCIP_CALL( SCIPresetClock(mainscip, wallclock) );
         SCIP_CALL( SCIPstartClock(mainscip, wallclock) );
      }
   }

   return SCIP_OKAY;
}


/** synchronize with other concurrent solvers */
SCIP_RETCODE SCIPsynchronize(
   SCIP*                 scip                /**< SCIP datastructure */
   )
{
   assert(scip != NULL);
   assert(scip->concurrent != NULL);

   SCIP_CALL( SCIPconcsolverSync(scip->concurrent->concsolver, scip->concurrent->mainscip->set) );

   scip->concurrent->mainscip->concurrent->solidx = scip->concurrent->mainscip->stat->solindex;

   if( scip->concurrent->eventglobalbnd != NULL )
      SCIPeventGlobalbndClearBoundChanges(scip->concurrent->eventglobalbnd);

   return SCIP_OKAY;
}

/** disables storing global bound changes */
void SCIPdisableConcurrentBoundStorage(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->concurrent != NULL);

   if( scip->concurrent->eventglobalbnd != NULL )
      SCIPeventGlobalbndDisableBoundStorage(scip->concurrent->eventglobalbnd);
}

/** enables storing global bound changes */
void SCIPenableConcurrentBoundStorage(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->concurrent != NULL);

   if( scip->concurrent->eventglobalbnd != NULL )
      SCIPeventGlobalbndEnableBoundStorage(scip->concurrent->eventglobalbnd);
}

/** gets total memory usage of all concurrent solvers together */
SCIP_Longint SCIPgetConcurrentMemTotal(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Longint memtotal = SCIPgetMemTotal(scip);

   assert(scip != NULL);

   if( scip->concurrent == NULL || scip->concurrent->mainscip != scip || scip->concurrent->concsolver == NULL )
      return memtotal;
   else
   {
      SCIP_Longint concmemtotal = SCIPconcsolverGetMemTotal(scip->concurrent->concsolver);
      return MAX(memtotal, concmemtotal);
   }
}

/** gets the dualbound in the last synchronization */
SCIP_Real SCIPgetConcurrentDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYNCSTORE* syncstore;

   assert(scip != NULL);

   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPsyncstoreGetLastLowerbound(syncstore));
}

/** gets the primalbound in the last synchronization */
SCIP_Real SCIPgetConcurrentPrimalbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYNCSTORE* syncstore;

   assert(scip != NULL);

   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPsyncstoreGetLastUpperbound(syncstore));
}

/** gets the gap in the last synchronization */
SCIP_Real SCIPgetConcurrentGap(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real primalbound;
   SCIP_Real dualbound;

   primalbound = SCIPgetConcurrentPrimalbound(scip);
   dualbound = SCIPgetConcurrentDualbound(scip);

   return SCIPcomputeGap(SCIPepsilon(scip), SCIPinfinity(scip), primalbound, dualbound);
}

/** gives the total number of tightened bounds received from other concurrent solvers */
SCIP_Longint SCIPgetConcurrentNTightenedBnds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip->concurrent != NULL);

   return scip->concurrent->propsync != NULL ? SCIPpropSyncGetNTightenedBnds(scip->concurrent->propsync) : 0;
}

/** gives the total number of tightened bounds for integer variables received from
 *  other concurrent solvers */
SCIP_Longint SCIPgetConcurrentNTightenedIntBnds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip->concurrent != NULL);

   return scip->concurrent->propsync != NULL ? SCIPpropSyncGetNTightenedIntBnds(scip->concurrent->propsync) : 0;
}

/** pass a solution to the given SCIP instance using that was received via synchronization by using
 * the sync heuristic */
SCIP_RETCODE SCIPaddConcurrentSol(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_SOL*             sol                 /**< solution */
   )
{
   assert(scip != NULL);
   assert(scip->concurrent != NULL);
   assert(sol != NULL);

   SCIP_CALL( SCIPheurSyncPassSol(scip, scip->concurrent->heursync, sol) );

   return SCIP_OKAY;
}

/** adds a global boundchange to the given SCIP, by passing it to the sync propagator */
SCIP_RETCODE SCIPaddConcurrentBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for bound */
   SCIP_Real             val,                /**< value of bound */
   SCIP_BOUNDTYPE        bndtype             /**< type of bound */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(scip->concurrent != NULL);
   assert(scip->concurrent->propsync != NULL);

   SCIP_CALL( SCIPpropSyncAddBndchg(scip->concurrent->mainscip, scip->concurrent->propsync, var, val, bndtype) );

   return SCIP_OKAY;
}

/** copy the nodenumber, depth, time, and runnumber of one solution to another one */
SCIP_RETCODE SCIPcopySolStats(
   SCIP_SOL*             source,             /**< source for solution statistics */
   SCIP_SOL*             target              /**< target for solution statistics */
   )
{
   assert(source != NULL);
   assert(target != NULL);

   target->depth = source->depth;
   target->time = source->time;
   target->nodenum = source->nodenum;
   target->runnum = source->runnum;

   return SCIP_OKAY;
}


/** get variable index of original variable that is the same between concurrent solvers */
int SCIPgetConcurrentVaridx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert(scip != NULL);
   assert(scip->concurrent != NULL);
   assert(scip->concurrent->varperm != NULL);
   assert(var != NULL);
   assert(SCIPvarIsOriginal(var));
   assert(SCIPvarGetIndex(var) < SCIPgetNOrigVars(scip));

   return scip->concurrent->varperm[SCIPvarGetIndex(var)];
}

/** is the solution new since the last synchronization point */
SCIP_Bool SCIPIsConcurrentSolNew(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< the solution */
   )
{
   assert(scip != NULL);
   assert(scip->concurrent != NULL);
   assert(sol != NULL);

   return SCIPsolGetIndex(sol) >= scip->concurrent->solidx;
}

/** gets the global lower bound changes since the last synchronization point */
SCIP_BOUNDSTORE* SCIPgetConcurrentGlobalBoundChanges(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->concurrent != NULL);

   if( scip->concurrent->eventglobalbnd != NULL )
      return SCIPeventGlobalbndGetBoundChanges(scip->concurrent->eventglobalbnd);

   return NULL;
}

/** executes the concurrent solver corresponding to the current thread */
static
SCIP_RETCODE execConcsolver(
   void*                 args                /**< SCIP data structure passed in as a void pointer */
   )
{
   SCIP* scip;

   assert(args != NULL);

   scip = (SCIP*) args;

   SCIP_CALL( SCIPconcsolverExec(scip->set->concsolvers[SCIPtpiGetThreadNum()]) );
   SCIP_CALL( SCIPconcsolverSync(scip->set->concsolvers[SCIPtpiGetThreadNum()], scip->set) );

   return SCIP_OKAY;
}

/** start solving in parallel using the given set of concurrent solvers */
SCIP_RETCODE SCIPconcurrentSolve(
   SCIP*                 scip                /**< pointer to scip datastructure */
   )
{
   SCIP_SYNCSTORE*   syncstore;
   int               idx;
   int               jobid;
   int               i;
   SCIP_RETCODE      retcode;
   SCIP_CONCSOLVER** concsolvers;
   int               nconcsolvers;

   assert(scip != NULL);

   syncstore = SCIPgetSyncstore(scip);
   concsolvers = scip->set->concsolvers;
   nconcsolvers = scip->set->nconcsolvers;

   assert(SCIPsyncstoreIsInitialized(syncstore));
   assert(SCIPsyncstoreGetNSolvers(syncstore) == nconcsolvers);

   SCIPsyncstoreSetSolveIsStopped(syncstore, FALSE);
   jobid = SCIPtpiGetNewJobID();

   TPI_PARA
   {
      TPI_SINGLE
      {
         for( i = 0; i < nconcsolvers; ++i )
         {
            /* cppcheck-suppress unassignedVariable */
            SCIP_JOB*         job;
            SCIP_SUBMITSTATUS status;

            SCIP_CALL_ABORT( SCIPtpiCreateJob(&job, jobid, execConcsolver, scip) );
            SCIP_CALL_ABORT( SCIPtpiSumbitJob(job, &status) );

            assert(status == SCIP_SUBMIT_SUCCESS);
         }
      }
   }

   retcode = SCIPtpiCollectJobs(jobid);
   idx = SCIPsyncstoreGetWinner(syncstore);
   assert(idx >= 0 && idx < nconcsolvers);

   SCIP_CALL( SCIPconcsolverGetSolvingData(concsolvers[idx], scip) );

   return retcode;
}

/** copy solving statistics */
SCIP_RETCODE SCIPcopyConcurrentSolvingStats(
   SCIP*                 source,             /**< SCIP data structure */
   SCIP*                 target              /**< target SCIP data structure */
   )
{
   SCIP_Real     tmptime;
   SCIP_HEUR*    heur;
   SCIP_NODE*    root;
   SCIP_PROP*    prop;
   SCIP_SEPA*    sepa;
   SCIP_PRESOL*  presol;
   SCIP_HEUR**   heurs;
   int           nheurs;
   SCIP_PROP**   props;
   int           nprops;
   SCIP_SEPA**   sepas;
   int           nsepas;
   SCIP_PRESOL** presols;
   int           npresols;
   int           i;

   assert(source != NULL);
   assert(target != NULL);

   heurs = SCIPgetHeurs(target);
   nheurs = SCIPgetNHeurs(target);

   for( i = 0; i < nheurs; ++i )
   {
      heur = SCIPfindHeur(source, SCIPheurGetName(heurs[i]));

      if( heur != NULL )
      {
         heurs[i]->nbestsolsfound += heur->nbestsolsfound;
         heurs[i]->ncalls += heur->ncalls;
         heurs[i]->nsolsfound += heur->nsolsfound;
         /* TODO divesets */
         tmptime = SCIPgetClockTime(target, heurs[i]->setuptime);
         tmptime += SCIPgetClockTime(source, heur->setuptime);
         SCIP_CALL( SCIPsetClockTime(target, heurs[i]->setuptime, tmptime) );

         tmptime = SCIPgetClockTime(target, heurs[i]->heurclock);
         tmptime += SCIPgetClockTime(source, heur->heurclock);
         SCIP_CALL( SCIPsetClockTime(target, heurs[i]->heurclock, tmptime) );
      }
   }

   props = SCIPgetProps(target);
   nprops = SCIPgetNProps(target);

   for( i = 0; i < nprops; ++i )
   {
      prop = SCIPfindProp(source, SCIPpropGetName(props[i]));

      if( prop != NULL )
      {
         props[i]->ncalls += prop->ncalls;
         props[i]->nrespropcalls += prop->nrespropcalls;
         props[i]->ncutoffs += prop->ncutoffs;
         props[i]->ndomredsfound += prop->ndomredsfound;

         tmptime = SCIPgetClockTime(target, props[i]->proptime);
         tmptime += SCIPgetClockTime(source, prop->proptime);
         SCIP_CALL( SCIPsetClockTime(target, props[i]->proptime, tmptime) );

         tmptime = SCIPgetClockTime(target, props[i]->sbproptime);
         tmptime += SCIPgetClockTime(source, prop->sbproptime);
         SCIP_CALL( SCIPsetClockTime(target, props[i]->sbproptime, tmptime) );

         tmptime = SCIPgetClockTime(target, props[i]->resproptime);
         tmptime += SCIPgetClockTime(source, prop->resproptime);
         SCIP_CALL( SCIPsetClockTime(target, props[i]->resproptime, tmptime) );

         tmptime = SCIPgetClockTime(target, props[i]->presoltime);
         tmptime += SCIPgetClockTime(source, prop->presoltime);
         SCIP_CALL( SCIPsetClockTime(target, props[i]->presoltime, tmptime) );

         tmptime = SCIPgetClockTime(target, props[i]->setuptime);
         tmptime += SCIPgetClockTime(source, prop->setuptime);
         SCIP_CALL( SCIPsetClockTime(target, props[i]->setuptime, tmptime) );
      }
   }

   presols = SCIPgetPresols(target);
   npresols = SCIPgetNPresols(target);

   for( i = 0; i < npresols; ++i )
   {
      presol = SCIPfindPresol(source, SCIPpresolGetName(presols[i]));

      if( presol != NULL )
      {
         presols[i]->ncalls += presol->ncalls;
         presols[i]->nfixedvars += presol->nfixedvars;
         presols[i]->naggrvars += presol->naggrvars;
         presols[i]->nchgvartypes += presol->nchgvartypes;
         presols[i]->nchgbds += presol->nchgbds;
         presols[i]->naddholes += presol->naddholes;
         presols[i]->ndelconss += presol->ndelconss;
         presols[i]->naddconss += presol->naddconss;
         presols[i]->nupgdconss += presol->nupgdconss;
         presols[i]->nchgcoefs += presol->nchgcoefs;
         presols[i]->nchgsides += presol->nchgsides;
         presols[i]->nfixedvars += presol->nfixedvars;
         presols[i]->nfixedvars += presol->nfixedvars;
         presols[i]->nfixedvars += presol->nfixedvars;

         tmptime = SCIPgetClockTime(target, presols[i]->setuptime);
         tmptime += SCIPgetClockTime(source, presol->setuptime);
         SCIP_CALL( SCIPsetClockTime(target, presols[i]->setuptime, tmptime) );

         tmptime = SCIPgetClockTime(target, presols[i]->presolclock);
         tmptime += SCIPgetClockTime(source, presol->presolclock);
         SCIP_CALL( SCIPsetClockTime(target, presols[i]->presolclock, tmptime) );
      }
   }

   sepas = SCIPgetSepas(target);
   nsepas = SCIPgetNSepas(target);

   for( i = 0; i < nsepas; ++i )
   {
      sepa = SCIPfindSepa(source, SCIPsepaGetName(sepas[i]));

      if( sepa != NULL )
      {
         sepas[i]->lastsepanode = sepa->lastsepanode;
         sepas[i]->ncalls += sepa->ncalls;
         sepas[i]->ncutoffs += sepa->ncutoffs;
         sepas[i]->ncutsfound += sepa->ncutsfound;
         sepas[i]->ncutsapplied += sepa->ncutsapplied;
         sepas[i]->nconssfound += sepa->nconssfound;
         sepas[i]->ndomredsfound += sepa->ndomredsfound;
         sepas[i]->maxbounddist = MAX(sepas[i]->maxbounddist, sepa->maxbounddist);

         tmptime = SCIPgetClockTime(target, sepas[i]->setuptime);
         tmptime += SCIPgetClockTime(source, sepa->setuptime);
         SCIP_CALL( SCIPsetClockTime(target, sepas[i]->setuptime, tmptime) );

         tmptime = SCIPgetClockTime(target, sepas[i]->sepaclock);
         tmptime += SCIPgetClockTime(source, sepa->sepaclock);
         SCIP_CALL( SCIPsetClockTime(target, sepas[i]->sepaclock, tmptime) );
      }
   }

   target->primal->nsolsfound = source->primal->nsolsfound;
   target->primal->nbestsolsfound = source->primal->nbestsolsfound;
   target->primal->nlimsolsfound = source->primal->nlimsolsfound;
   SCIPprobSetDualbound(target->transprob, SCIPprobExternObjval(target->transprob, target->origprob, target->set, SCIPgetDualbound(source)));
   root = SCIPgetRootNode(target);

   if( root != NULL )
   {
      /* in the copied SCIP the dualbound is in the transformed space of the target */
      SCIP_CALL( SCIPupdateNodeLowerbound(target, root, SCIPgetDualbound(source)) );
   }

   target->stat->nlpiterations = source->stat->nlpiterations;
   target->stat->nrootlpiterations = source->stat->nrootlpiterations;
   target->stat->nrootfirstlpiterations = source->stat->nrootfirstlpiterations;
   target->stat->nprimallpiterations = source->stat->nprimallpiterations;
   target->stat->nduallpiterations = source->stat->nduallpiterations;
   target->stat->nlexduallpiterations = source->stat->nlexduallpiterations;
   target->stat->nbarrierlpiterations = source->stat->nbarrierlpiterations;
   target->stat->nprimalresolvelpiterations = source->stat->nprimalresolvelpiterations;
   target->stat->ndualresolvelpiterations = source->stat->ndualresolvelpiterations;
   target->stat->nlexdualresolvelpiterations = source->stat->nlexdualresolvelpiterations;
   target->stat->nnodelpiterations = source->stat->nnodelpiterations;
   target->stat->ninitlpiterations = source->stat->ninitlpiterations;
   target->stat->ndivinglpiterations = source->stat->ndivinglpiterations;
   target->stat->ndivesetlpiterations = source->stat->ndivesetlpiterations;
   target->stat->nsbdivinglpiterations = source->stat->nsbdivinglpiterations;
   target->stat->nsblpiterations = source->stat->nsblpiterations;
   target->stat->nrootsblpiterations = source->stat->nrootsblpiterations;
   target->stat->nconflictlpiterations = source->stat->nconflictlpiterations;
   target->stat->nnodes = source->stat->nnodes;
   target->stat->ninternalnodes = source->stat->ninternalnodes;
   target->stat->nobjleaves = source->stat->nobjleaves;
   target->stat->nfeasleaves = source->stat->nfeasleaves;
   target->stat->ninfeasleaves = source->stat->ninfeasleaves;
   target->stat->ntotalnodes = source->stat->ntotalnodes;
   target->stat->ntotalinternalnodes = source->stat->ntotalinternalnodes;
   target->stat->ncreatednodes = source->stat->ncreatednodes;
   target->stat->ncreatednodesrun = source->stat->ncreatednodesrun;
   target->stat->nactivatednodes = source->stat->nactivatednodes;
   target->stat->ndeactivatednodes = source->stat->ndeactivatednodes;
   target->stat->nearlybacktracks = source->stat->nearlybacktracks;
   target->stat->nnodesaboverefbound = source->stat->nnodesaboverefbound;
   target->stat->nbacktracks = source->stat->nbacktracks;
   target->stat->ndelayedcutoffs = source->stat->ndelayedcutoffs;
   target->stat->nreprops = source->stat->nreprops;
   target->stat->nrepropboundchgs = source->stat->nrepropboundchgs;
   target->stat->nrepropcutoffs = source->stat->nrepropcutoffs;
   target->stat->nlpsolsfound = source->stat->nlpsolsfound;
   target->stat->npssolsfound = source->stat->npssolsfound;
   target->stat->nsbsolsfound = source->stat->nsbsolsfound;
   target->stat->nlpbestsolsfound = source->stat->nlpbestsolsfound;
   target->stat->npsbestsolsfound = source->stat->npsbestsolsfound;
   target->stat->nsbbestsolsfound = source->stat->nsbbestsolsfound;
   target->stat->nexternalsolsfound = source->stat->nexternalsolsfound;
   target->stat->lastdispnode = source->stat->lastdispnode;
   target->stat->lastdivenode = source->stat->lastdivenode;
   target->stat->lastconflictnode = source->stat->lastconflictnode;
   target->stat->bestsolnode = source->stat->bestsolnode;
   target->stat->domchgcount = source->stat->domchgcount;
   target->stat->nboundchgs = source->stat->nboundchgs;
   target->stat->nholechgs = source->stat->nholechgs;
   target->stat->nprobboundchgs = source->stat->nprobboundchgs;
   target->stat->nprobholechgs = source->stat->nprobholechgs;
   target->stat->nsbdowndomchgs = source->stat->nsbdowndomchgs;
   target->stat->nsbupdomchgs = source->stat->nsbupdomchgs;
   target->stat->nsbtimesiterlimhit = source->stat->nsbtimesiterlimhit;
   target->stat->nnodesbeforefirst = source->stat->nnodesbeforefirst;
   target->stat->ninitconssadded = source->stat->ninitconssadded;
   target->stat->firstlpdualbound = SCIPprobExternObjval(target->transprob, target->origprob, target->set, source->stat->firstlpdualbound);
   target->stat->rootlowerbound = SCIPprobExternObjval(source->transprob, source->origprob, source->set, source->stat->rootlowerbound);
   target->stat->vsidsweight = source->stat->vsidsweight;
   target->stat->firstprimalbound = SCIPprobExternObjval(target->transprob, target->origprob, target->set, source->stat->firstprimalbound);
   target->stat->firstprimaltime = source->stat->firstprimaltime;
   target->stat->firstsolgap = source->stat->firstsolgap;
   target->stat->lastsolgap = source->stat->lastsolgap;
   target->stat->primalzeroittime = source->stat->primalzeroittime;
   target->stat->dualzeroittime = source->stat->dualzeroittime;
   target->stat->barrierzeroittime = source->stat->barrierzeroittime;
   target->stat->maxcopytime = MAX(source->stat->maxcopytime, target->stat->maxcopytime);
   target->stat->mincopytime = MIN(source->stat->mincopytime, target->stat->mincopytime);
   target->stat->firstlptime = source->stat->firstlptime;
   target->stat->lastbranchvalue = source->stat->lastbranchvalue;
   target->stat->primaldualintegral = source->stat->primaldualintegral;
   target->stat->previousgap = source->stat->previousgap;
   target->stat->previntegralevaltime = source->stat->previntegralevaltime;
   target->stat->lastprimalbound = SCIPprobExternObjval(source->transprob, source->origprob, source->set, source->stat->lastprimalbound);
   target->stat->lastdualbound = SCIPprobExternObjval(source->transprob, source->origprob, source->set, source->stat->lastdualbound);
   target->stat->lastlowerbound = SCIPprobExternObjval(source->transprob, source->origprob, source->set, source->stat->lastlowerbound);
   target->stat->lastupperbound = SCIPprobExternObjval(source->transprob, source->origprob, source->set, source->stat->lastupperbound);
   target->stat->rootlpbestestimate = source->stat->rootlpbestestimate;
   target->stat->referencebound = source->stat->referencebound;

   /*tmptime = SCIPgetClockTime(target, target->stat->solvingtime);
   tmptime += SCIPgetClockTime(source, source->stat->solvingtime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->solvingtime, tmptime) );*/

   /* TODO */
   tmptime = SCIPgetClockTime(target, target->stat->solvingtimeoverall);
   tmptime += SCIPgetClockTime(source, source->stat->solvingtimeoverall);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->solvingtimeoverall, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->presolvingtime);
   tmptime += SCIPgetClockTime(source, source->stat->presolvingtime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->presolvingtime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->presolvingtimeoverall);
   tmptime += SCIPgetClockTime(source, source->stat->presolvingtimeoverall);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->presolvingtimeoverall, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->primallptime);
   tmptime += SCIPgetClockTime(source, source->stat->primallptime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->primallptime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->duallptime);
   tmptime += SCIPgetClockTime(source, source->stat->duallptime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->duallptime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->lexduallptime);
   tmptime += SCIPgetClockTime(source, source->stat->lexduallptime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->lexduallptime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->barrierlptime);
   tmptime += SCIPgetClockTime(source, source->stat->barrierlptime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->barrierlptime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->divinglptime);
   tmptime += SCIPgetClockTime(source, source->stat->divinglptime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->divinglptime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->strongbranchtime);
   tmptime += SCIPgetClockTime(source, source->stat->strongbranchtime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->strongbranchtime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->conflictlptime);
   tmptime += SCIPgetClockTime(source, source->stat->conflictlptime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->conflictlptime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->lpsoltime);
   tmptime += SCIPgetClockTime(source, source->stat->lpsoltime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->lpsoltime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->pseudosoltime);
   tmptime += SCIPgetClockTime(source, source->stat->pseudosoltime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->pseudosoltime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->sbsoltime);
   tmptime += SCIPgetClockTime(source, source->stat->sbsoltime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->sbsoltime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->nodeactivationtime);
   tmptime += SCIPgetClockTime(source, source->stat->nodeactivationtime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->nodeactivationtime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->nlpsoltime);
   tmptime += SCIPgetClockTime(source, source->stat->nlpsoltime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->nlpsoltime, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->strongpropclock);
   tmptime += SCIPgetClockTime(source, source->stat->strongpropclock);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->strongpropclock, tmptime) );

   tmptime = SCIPgetClockTime(target, target->stat->reoptupdatetime);
   tmptime += SCIPgetClockTime(source, source->stat->reoptupdatetime);
   SCIP_CALL( SCIPsetClockTime(target, target->stat->reoptupdatetime, tmptime) );

   heur = source->stat->firstprimalheur;

   if( heur != NULL )
      target->stat->firstprimalheur = SCIPfindHeur(target, SCIPheurGetName(heur));

   target->stat->status = source->stat->status;
   target->stat->lastbranchdir = source->stat->lastbranchdir;
   target->stat->lastsblpsolstats[0] = source->stat->lastsblpsolstats[0];
   target->stat->lastsblpsolstats[1] = source->stat->lastsblpsolstats[1];
   target->stat->nnz = source->stat->nnz;
   target->stat->lpcount = source->stat->lpcount;
   target->stat->nlps = source->stat->nlps;
   target->stat->nrootlps = source->stat->nrootlps;
   target->stat->nprimallps = source->stat->nprimallps;
   target->stat->nprimalzeroitlps = source->stat->nprimalzeroitlps;
   target->stat->nduallps = source->stat->nduallps;
   target->stat->ndualzeroitlps = source->stat->ndualzeroitlps;
   target->stat->nlexduallps = source->stat->nlexduallps;
   target->stat->nbarrierlps = source->stat->nbarrierlps;
   target->stat->nbarrierzeroitlps = source->stat->nbarrierzeroitlps;
   target->stat->nprimalresolvelps = source->stat->nprimalresolvelps;
   target->stat->ndualresolvelps = source->stat->ndualresolvelps;
   target->stat->nlexdualresolvelps = source->stat->nlexdualresolvelps;
   target->stat->nnodelps = source->stat->nnodelps;
   target->stat->ninitlps = source->stat->ninitlps;
   target->stat->ndivinglps = source->stat->ndivinglps;
   target->stat->ndivesetlps = source->stat->ndivesetlps;
   target->stat->nsbdivinglps = source->stat->nsbdivinglps;
   target->stat->nstrongbranchs = source->stat->nstrongbranchs;
   target->stat->nrootstrongbranchs = source->stat->nrootstrongbranchs;
   target->stat->nconflictlps = source->stat->nconflictlps;
   target->stat->nnlps = source->stat->nnlps;
   target->stat->nisstoppedcalls = source->stat->nisstoppedcalls;
   target->stat->totaldivesetdepth = source->stat->totaldivesetdepth;
   target->stat->ndivesetcalls = source->stat->ndivesetcalls;
   target->stat->nruns = source->stat->nruns;
   target->stat->nconfrestarts = source->stat->nconfrestarts;
   target->stat->nrootboundchgs = source->stat->nrootboundchgs;
   target->stat->nrootboundchgsrun = source->stat->nrootboundchgsrun;
   target->stat->nrootintfixings = source->stat->nrootintfixings;
   target->stat->nrootintfixingsrun = source->stat->nrootintfixingsrun;
   target->stat->prevrunnvars = source->stat->prevrunnvars;
   target->stat->npricerounds = source->stat->npricerounds;
   target->stat->nseparounds = source->stat->nseparounds;
   target->stat->maxdepth = source->stat->maxdepth;
   target->stat->maxtotaldepth = source->stat->maxtotaldepth;
   target->stat->plungedepth = source->stat->plungedepth;
   target->stat->npresolrounds += source->stat->npresolrounds;
   target->stat->npresolroundsfast += source->stat->npresolroundsfast;
   target->stat->npresolroundsmed += source->stat->npresolroundsmed;
   target->stat->npresolroundsext += source->stat->npresolroundsext;
   target->stat->npresolfixedvars += source->stat->npresolfixedvars;
   target->stat->npresolaggrvars += source->stat->npresolaggrvars;
   target->stat->npresolchgvartypes += source->stat->npresolchgvartypes;
   target->stat->npresolchgbds += source->stat->npresolchgbds;
   target->stat->npresoladdholes += source->stat->npresoladdholes;
   target->stat->npresoldelconss += source->stat->npresoldelconss;
   target->stat->npresoladdconss += source->stat->npresoladdconss;
   target->stat->npresolupgdconss += source->stat->npresolupgdconss;
   target->stat->npresolchgcoefs += source->stat->npresolchgcoefs;
   target->stat->npresolchgsides += source->stat->npresolchgsides;
   target->stat->nrunsbeforefirst = source->stat->nrunsbeforefirst;
   target->stat->firstprimaldepth = source->stat->firstprimaldepth;
   target->stat->ncopies += source->stat->ncopies;
   target->stat->nreoptruns = source->stat->nreoptruns;

   /* set the stage but do not set to earlier stage */
   target->set->stage = MAX(source->set->stage, target->set->stage);

   return SCIP_OKAY;
}
