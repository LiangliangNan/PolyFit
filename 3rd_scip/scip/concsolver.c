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

/**@file   concsolver.c
 * @brief  methods for concurrent solvers
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/concsolver.h"
#include "scip/set.h"
#include "scip/scip.h"
#include "scip/concurrent.h"

#include "scip/struct_concsolver.h"
#include "scip/struct_stat.h"
#include "scip/struct_scip.h"
#include "blockmemshell/memory.h"
#include "scip/syncstore.h"
#include "scip/boundstore.h"
#include "scip/clock.h"


/** creates a concurrent solver type */
SCIP_RETCODE SCIPconcsolverTypeCreate(
   SCIP_CONCSOLVERTYPE** concsolvertype,     /**< pointer to concurrent solver data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of concurrent solver */
   SCIP_Real             prefpriodefault,    /**< the default preferred priority of this concurrent solver type */
   SCIP_DECL_CONCSOLVERCREATEINST ((*concsolvercreateinst)),/**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERDESTROYINST ((*concsolverdestroyinst)),/**< data copy method of concurrent solver */
   SCIP_DECL_CONCSOLVERINITSEEDS ((*concsolverinitseeds)),/**< initialize random seeds of concurrent solver */
   SCIP_DECL_CONCSOLVEREXEC ((*concsolverexec)),/**< execution method of concurrent solver */
   SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA ((*concsolvercopysolvdata)),/**< method to copy solving data */
   SCIP_DECL_CONCSOLVERSTOP ((*concsolverstop)),/**< terminate solving in concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCWRITE ((*concsolversyncwrite)),/**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERSYNCREAD ((*concsolversyncread)),/**< synchronization method of concurrent solver */
   SCIP_DECL_CONCSOLVERTYPEFREEDATA ((*concsolvertypefreedata)),/**< method to free data of concurrent solver type */
   SCIP_CONCSOLVERTYPEDATA* data             /**< the concurent solver type's data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(concsolvertype != NULL);
   assert(name != NULL);
   assert(prefpriodefault >= 0.0 && prefpriodefault <= 1.0);

   assert(concsolvercreateinst != NULL);
   assert(concsolverdestroyinst != NULL);
   assert(concsolverexec != NULL);
   assert(concsolvercopysolvdata != NULL);
   assert(concsolverstop != NULL);
   assert(concsolversyncwrite != NULL);
   assert(concsolversyncread != NULL);

   SCIP_ALLOC( BMSallocMemory(concsolvertype) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*concsolvertype)->name, name, strlen(name) + 1) );

   (*concsolvertype)->data = data;
   (*concsolvertype)->ninstances = 0;
   (*concsolvertype)->concsolvercreateinst = concsolvercreateinst;
   (*concsolvertype)->concsolverdestroyinst = concsolverdestroyinst;
   (*concsolvertype)->concsolverinitseeds = concsolverinitseeds;
   (*concsolvertype)->concsolverexec = concsolverexec;
   (*concsolvertype)->concsolvercopysolvdata = concsolvercopysolvdata;
   (*concsolvertype)->concsolverstop = concsolverstop;
   (*concsolvertype)->concsolversyncwrite = concsolversyncwrite;
   (*concsolvertype)->concsolversyncread = concsolversyncread;
   (*concsolvertype)->concsolvertypefreedata = concsolvertypefreedata;

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "concurrent/%s/prefprio", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "the preferred number concurrent solvers of type <%s> with respect to the number of threads", name);
   SCIP_CALL( SCIPsetAddRealParam(set, messagehdlr, blkmem, paramname, paramdesc,
                                  &(*concsolvertype)->prefprio, FALSE, prefpriodefault, 0.0, 1.0,
                                  NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** frees all memory of a concurrent solver type */
void SCIPconcsolverTypeFree(
   SCIP_CONCSOLVERTYPE** concsolvertype      /**< pointer to concurrent solver data structure */
   )
{
   assert(concsolvertype != NULL);
   assert(*concsolvertype != NULL);

   if( (*concsolvertype)->concsolvertypefreedata != NULL )
      (*concsolvertype)->concsolvertypefreedata(&(*concsolvertype)->data);

   BMSfreeMemoryArray(&(*concsolvertype)->name);
   BMSfreeMemory(concsolvertype);
}

/** gets the data of a concurrent solver type */
SCIP_CONCSOLVERTYPEDATA* SCIPconcsolverTypeGetData(
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   )
{
   assert(concsolvertype != NULL);

   return concsolvertype->data;
}

/** sets the data of a concurrent solver type */
void SCIPconcsolverTypeSetData(
   SCIP_CONCSOLVERTYPE*  concsolvertype,     /**< concurrent solver type */
   SCIP_CONCSOLVERTYPEDATA* data             /**< the concurrent solver's data */
   )
{
   assert(concsolvertype != NULL);

   concsolvertype->data = data;
}

/** gets the name of a concurrent solver type */
char* SCIPconcsolverTypeGetName(
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   )
{
   assert(concsolvertype != NULL);

   return concsolvertype->name;
}

/** gets the preferred priority from a concurrent solver type */
SCIP_Real SCIPconcsolverTypeGetPrefPrio(
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   )
{
   assert(concsolvertype != NULL);

   return concsolvertype->prefprio;
}

/** creates an instance of the given concurrent solver type */
SCIP_RETCODE SCIPconcsolverCreateInstance(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVERTYPE*  concsolvertype,     /**< concurrent solver type to create */
   SCIP_CONCSOLVER**     concsolver          /**< pointer to return concurrent solver instance */
   )
{
   char instancename[SCIP_MAXSTRLEN];

   ++concsolvertype->ninstances;
   (void) SCIPsnprintf(instancename, SCIP_MAXSTRLEN, "%s-%i", concsolvertype->name, concsolvertype->ninstances);

   SCIP_ALLOC( BMSallocMemory(concsolver) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*concsolver)->name, instancename, strlen(instancename) + 1) );

   (*concsolver)->type = concsolvertype;

   /* initialize counters for statistics  */
   (*concsolver)->nsolsrecvd = 0;
   (*concsolver)->nsolsshared = 0;
   (*concsolver)->ntighterbnds = 0;
   (*concsolver)->ntighterintbnds = 0;
   SCIP_CALL( SCIPcreateWallClock(set->scip, &(*concsolver)->totalsynctime) );

   /* initialize synchronization fields */
   (*concsolver)->nsyncs = 0;
   (*concsolver)->syncdelay = 0.0;

   /* in deterministic mode use number of nonzeros and variables to get a good initial synchronization frequency
    * in opportunistic mode use the frequency as set by the user
    */
   if( set->parallel_mode == (int) SCIP_PARA_DETERMINISTIC )
      (*concsolver)->syncfreq = 0.01 * set->scip->stat->nnz * SCIPgetNVars(set->scip) * set->concurrent_freqinit;
   else
      (*concsolver)->syncfreq = set->concurrent_freqinit;

   (*concsolver)->syncdata = NULL;

   SCIPdebugMessage("concsolver %s initialized sync freq to %f\n", (*concsolver)->name, (*concsolver)->syncfreq);
   /* register concurrent solver */
   (*concsolver)->idx = SCIPgetNConcurrentSolvers(set->scip);
   SCIP_CALL( concsolvertype->concsolvercreateinst(set->scip, concsolvertype, *concsolver) );
   SCIP_CALL( SCIPaddConcurrentSolver(set->scip, *concsolver) );

   return SCIP_OKAY;
}

/** destroys an instance of the given concurrent solver */
SCIP_RETCODE SCIPconcsolverDestroyInstance(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVER**     concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);
   assert((*concsolver)->type != NULL);
   assert(set != NULL);
   assert((*concsolver)->type->concsolverdestroyinst != NULL);

   SCIP_CALL( (*concsolver)->type->concsolverdestroyinst(set->scip, *concsolver) );
   --(*concsolver)->type->ninstances;

   SCIP_CALL( SCIPfreeClock(set->scip, &(*concsolver)->totalsynctime) );
   BMSfreeMemoryArray(&(*concsolver)->name);

   BMSfreeMemory(concsolver);

   return SCIP_OKAY;
}

/** gets the data of a concurrent solver */
SCIP_CONCSOLVERDATA* SCIPconcsolverGetData(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->data;
}

/** sets the data of a concurrent solver */
void SCIPconcsolverSetData(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP_CONCSOLVERDATA*  data                /**< the concurrent solver's data */
   )
{
   assert(concsolver != NULL);

   concsolver->data = data;
}

/** gets the name of a concurrent solver */
char* SCIPconcsolverGetName(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->name;
}

/** initializes the random seeds of a concurrent solver */
SCIP_RETCODE SCIPconcsolverInitSeeds(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   unsigned int          seed                /**< seed for initializing the solver's internal random seeds */
   )
{
   assert(concsolver != NULL);
   assert(concsolver->type != NULL);

   if( concsolver->type->concsolverinitseeds != NULL )
      SCIP_CALL( concsolver->type->concsolverinitseeds(concsolver, seed) );

   return SCIP_OKAY;
}

/** start the solving process of a concurrent solver */
SCIP_RETCODE SCIPconcsolverExec(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);
   assert(concsolver->type != NULL);
   assert(concsolver->type->concsolverexec != NULL);

   /* set the stopped flag to false */
   concsolver->stopped = FALSE;

   /* then call the execute callback */
   SCIP_CALL( concsolver->type->concsolverexec(concsolver, &concsolver->solvingtime, &concsolver->nlpiterations, &concsolver->nnodes) );

   return SCIP_OKAY;
}

/** gets solving data of concurrent solver and stores it in the given SCIP instance */
SCIP_RETCODE SCIPconcsolverGetSolvingData(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP*                 scip                /**< SCIP datastructure */
   )
{
   assert(concsolver != NULL);
   assert(concsolver->type != NULL);
   assert(concsolver->type->concsolvercopysolvdata != NULL);

   return concsolver->type->concsolvercopysolvdata(concsolver, scip);
}

/** interrupt solving in a concurrent solver */
SCIP_RETCODE SCIPconcsolverStop(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);
   assert(concsolver->type != NULL);
   assert(concsolver->type->concsolverstop != NULL);

   SCIP_CALL( concsolver->type->concsolverstop(concsolver) );

   /* set the stopped flag to true */
   concsolver->stopped = TRUE;

   return SCIP_OKAY;
}

/** let the given concurrent solver synchronize, i.e. pass its own solutions and bounds to
 *  the SPI.
 */
SCIP_RETCODE SCIPconcsolverSync(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_SYNCDATA*   syncdata;
   SCIP_SYNCSTORE*  syncstore;
   int              nsols;
   int              ntighterintbnds;
   int              ntighterbnds;
   SCIP_CONCSOLVERTYPE* concsolvertype;

   assert(concsolver != NULL);
   assert(concsolver->type != NULL);
   assert(concsolver->type->concsolversyncwrite != NULL);
   assert(concsolver->type->concsolversyncread != NULL);

   if( concsolver->stopped )
      return SCIP_OKAY;

   SCIP_CALL( SCIPstartClock(set->scip, concsolver->totalsynctime) );

   concsolvertype = concsolver->type;

   syncstore = SCIPgetSyncstore(set->scip);
   assert(syncstore != NULL);

   SCIP_CALL( SCIPsyncstoreStartSync(syncstore, concsolver->nsyncs, &syncdata) );

   if( syncdata == NULL )
   {
      SCIP_CALL( SCIPstopClock(set->scip, concsolver->totalsynctime) );
      return SCIP_OKAY;
   }

   SCIPdebugMessage("concsolver %s starts sync %lli\n", concsolver->name, concsolver->nsyncs);

   SCIP_CALL( concsolvertype->concsolversyncwrite(concsolver, syncstore, syncdata, set->concurrent_nbestsols, set->concurrent_maxnsols, &nsols) );
   concsolver->nsolsshared += nsols;

   if( SCIPsyncdataGetStatus(syncdata) != SCIP_STATUS_UNKNOWN )
   {
      SCIP_CALL( SCIPconcsolverStop(concsolver) );
   }
   else if( SCIPsyncdataGetNSynced(syncdata) == SCIPsyncstoreGetNSolvers(syncstore) - 1 )
   {
      /* if this is the last concurrent solver that is synchronizing for this synchronization data
       * it will adjust the synchronization frequency using the progress on the gap
       */
      SCIP_Bool lbok;
      SCIP_Bool ubok;
      SCIP_Real progress;
      SCIP_Real prevub;
      SCIP_Real prevlb;
      SCIP_Real newub;
      SCIP_Real newlb;
      SCIP_Real freqfactor;
      SCIP_Real newsyncfreq;
      SCIP_SYNCDATA* prevsync;

      if( concsolver->nsyncs == 0 )
      {
         SCIPsyncdataSetSyncFreq(syncstore, syncdata, concsolver->syncfreq);
      }
      else
      {
         prevsync = SCIPsyncstoreGetSyncdata(syncstore, concsolver->nsyncs - 1);
         assert(SCIPsyncdataGetNSynced(prevsync) == SCIPsyncstoreGetNSolvers(syncstore));

         prevub = SCIPsyncdataGetUpperbound(prevsync);
         prevlb = SCIPsyncdataGetLowerbound(prevsync);
         newub = SCIPsyncdataGetUpperbound(syncdata);
         newlb = SCIPsyncdataGetLowerbound(syncdata);
         lbok = prevlb > -SCIPsetInfinity(set);
         ubok = prevub < SCIPsetInfinity(set);

         if( lbok && ubok )
            progress = SCIPrelDiff(prevub - prevlb, newub - newlb);
         else if( lbok )
            progress = SCIPrelDiff(newlb, prevlb);
         else if( ubok )
            progress = SCIPrelDiff(prevub, newub);
         else if( !SCIPsetIsInfinity(set, -newlb) || !SCIPsetIsInfinity(set, newub) ||
                  SCIPboundstoreGetNChgs(SCIPsyncdataGetBoundChgs(syncdata)) > 0 )
            progress = set->concurrent_targetprogress;
         else
            progress = 0.0;

         /* should not be negative */
         assert(SCIPsetIsGE(set, progress, 0.0));

         if( progress < 0.5 * set->concurrent_targetprogress )
            freqfactor = set->concurrent_freqfactor;
         else if( progress > 2 * set->concurrent_targetprogress )
            freqfactor = 0.5 + 0.5 / set->concurrent_freqfactor;
         else
            freqfactor = 1.0;

         SCIPdebugMessage("syncfreq is %g and freqfactor is %f due to progress %f\n", concsolver->syncfreq, freqfactor, progress);
         newsyncfreq = concsolver->syncfreq * freqfactor;
         SCIPsyncdataSetSyncFreq(syncstore, syncdata, newsyncfreq);
         SCIPdebugMessage("new syncfreq is %g\n", SCIPsyncdataGetSyncFreq(syncdata));
      }
   }

   SCIPdebugMessage("concsolver %s finishing sync %lli\n", concsolver->name, concsolver->nsyncs);

   SCIP_CALL( SCIPsyncstoreFinishSync(syncstore, &syncdata) );
   ++concsolver->nsyncs;

   concsolver->syncdelay += concsolver->timesincelastsync;

   syncdata = SCIPsyncstoreGetNextSyncdata(syncstore, concsolver->syncdata, concsolver->syncfreq, concsolver->nsyncs, &concsolver->syncdelay);

   while( syncdata != NULL )
   {
      SCIP_CALL( SCIPsyncstoreEnsureAllSynced(syncstore, syncdata) );
      concsolver->syncdata = syncdata;
      SCIP_CALL( concsolvertype->concsolversyncread(concsolver, syncstore, syncdata, &nsols, &ntighterbnds, &ntighterintbnds) );
      concsolver->ntighterbnds += ntighterbnds;
      concsolver->ntighterintbnds += ntighterintbnds;
      concsolver->nsolsrecvd += nsols;
      SCIPdebugMessage("syncfreq before reading the next syncdata is %g\n", concsolver->syncfreq);
      concsolver->syncfreq = SCIPsyncdataGetSyncFreq(concsolver->syncdata);
      SCIPdebugMessage("syncfreq after reading the next syncdata is %g\n", concsolver->syncfreq);
      syncdata = SCIPsyncstoreGetNextSyncdata(syncstore, concsolver->syncdata, concsolver->syncfreq, concsolver->nsyncs, &concsolver->syncdelay);
   }

   SCIP_CALL( SCIPstopClock(set->scip, concsolver->totalsynctime) );

   return SCIP_OKAY;
}

/** gets the current synchronization frequency of the concurent solver */
SCIP_Real SCIPconcsolverGetSyncFreq(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->syncfreq;
}

/** gets the total memory used by the concurent solver */
SCIP_Longint SCIPconcsolverGetMemTotal(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->syncdata != NULL ? SCIPsyncdataGetMemTotal(concsolver->syncdata) : 0;
}

/** sets the time elapsed since the last synchronization. Must be set before the synchronization is
 *  started.
 */
void SCIPconcsolverSetTimeSinceLastSync(
   SCIP_CONCSOLVER*      concsolver,         /**< concurrent solver */
   SCIP_Real             time                /**< the time passed since the last synchronization */
   )
{
   assert(concsolver != NULL);

   concsolver->timesincelastsync = time;
}

/** gets the solving time of the concurrent solver */
SCIP_Real SCIPconcsolverGetSolvingTime(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->solvingtime;
}

/** gets the time spent for synchronization for the concurrent solver */
SCIP_Real SCIPconcsolverGetSyncTime(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return SCIPclockGetTime(concsolver->totalsynctime);
}

/** gets the number of lp iterations the concurrent solver used */
SCIP_Longint SCIPconcsolverGetNLPIterations(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->nlpiterations;
}

/** gets the number of branch and bound nodes the concurrent solver used */
SCIP_Longint SCIPconcsolverGetNNodes(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->nnodes;
}

/** gets the number of solutions the concurrent solver received during synchronization */
SCIP_Longint SCIPconcsolverGetNSolsRecvd(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->nsolsrecvd;
}

/** gets the number of solutions the concurrent solver shared during synchronization */
SCIP_Longint SCIPconcsolverGetNSolsShared(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->nsolsshared;
}

/** gets the number of tighter global variable bounds the solver received */
SCIP_Longint SCIPconcsolverGetNTighterBnds(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->ntighterbnds;
}

/** gets the number of tighter global variable bounds of integer variables the solver received */
SCIP_Longint SCIPconcsolverGetNTighterIntBnds(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->ntighterintbnds;
}

/** gets index of concurrent solver */
int SCIPconcsolverGetIdx(
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(concsolver != NULL);

   return concsolver->idx;
}
