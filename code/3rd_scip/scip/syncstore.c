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

/**@file   syncstore.c
 * @ingroup PARALLEL
 * @brief  the function definitions of the synchronization store
 * @author Robert Lion Gottwald
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/pub_message.h"
#include "scip/concsolver.h"
#include "scip/struct_concsolver.h"
#include "scip/prob.h"
#include "scip/scip.h"
#include "blockmemshell/memory.h"
#include "tpi/tpi.h"
#include "scip/struct_syncstore.h"
#include "scip/concurrent.h"
#include "scip/syncstore.h"
#include "scip/boundstore.h"


/** computes the size of the array of synchronization datas, such that
 *  it cannot ever happen that a synchronization data is reused while still
 *  not read by any thread */
static
int getNSyncdata(
   SCIP*                 scip                /**< SCIP main datastructure */
   )
{
   int maxnsyncdelay;
   SCIP_CALL_ABORT( SCIPgetIntParam(scip, "concurrent/sync/maxnsyncdelay", &maxnsyncdelay) );

   return 2 * (maxnsyncdelay + 1);
}

/** creates and captures a new synchronization store */
SCIP_RETCODE SCIPsyncstoreCreate(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to return the created synchronization store */
   )
{
   assert(syncstore != NULL);

   SCIPdebugMessage("SCIPsyncstoreCreate()\n");

   SCIP_ALLOC( BMSallocMemory(syncstore) );

   (*syncstore)->mode = SCIP_PARA_DETERMINISTIC;                      /* initialising the mode */
   (*syncstore)->initialized = FALSE;
   (*syncstore)->syncdata = NULL;
   (*syncstore)->stopped = FALSE;
   (*syncstore)->nuses = 1;
   SCIP_CALL( SCIPtpiInitLock(&(*syncstore)->lock) );

   return SCIP_OKAY;
}

/** releases a synchronization store */
SCIP_RETCODE SCIPsyncstoreRelease(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to the synchronization store */
   )
{
   int references;

   assert(syncstore != NULL);
   if( *syncstore == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPtpiAcquireLock(&(*syncstore)->lock) );
   (*syncstore)->nuses -= 1;
   references = (*syncstore)->nuses;
   SCIP_CALL( SCIPtpiReleaseLock(&(*syncstore)->lock) );

   if( references == 0 )
   {
      if( (*syncstore)->initialized )
      {
         SCIP_CALL( SCIPsyncstoreExit(*syncstore) );
      }

      assert(!(*syncstore)->initialized);
      SCIPtpiDestroyLock(&(*syncstore)->lock);
      BMSfreeMemory(syncstore);
   }
   else
   {
      *syncstore = NULL;
   }

   return SCIP_OKAY;
}

/** captures a synchronization store */
SCIP_RETCODE SCIPsyncstoreCapture(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   SCIP_CALL( SCIPtpiAcquireLock(&syncstore->lock) );

   ++(syncstore->nuses);

   SCIP_CALL( SCIPtpiReleaseLock(&syncstore->lock) );

   return SCIP_OKAY;
}

/** initialize the syncstore for the given SCIP instance */
SCIP_RETCODE SCIPsyncstoreInit(
   SCIP*                 scip                /**< SCIP main datastructure */
   )
{
   SCIP_SYNCSTORE* syncstore;
   int i;
   int j;
   int paramode;

   assert(scip != NULL);
   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);
   syncstore->mainscip = scip;
   syncstore->lastsync = NULL;
   syncstore->nsolvers = SCIPgetNConcurrentSolvers(scip);

   syncstore->ninitvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsols", &syncstore->maxnsols) );
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsyncdelay", &syncstore->maxnsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/minsyncdelay", &syncstore->minsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqinit", &syncstore->syncfreqinit) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqmax", &syncstore->syncfreqmax) );
   syncstore->nsyncdata = getNSyncdata(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &(syncstore->syncdata), syncstore->nsyncdata) );

   for( i = 0; i < syncstore->nsyncdata; ++i )
   {
      syncstore->syncdata[i].syncnum = -1;
      SCIP_CALL( SCIPboundstoreCreate(syncstore->mainscip, &syncstore->syncdata[i].boundstore, syncstore->ninitvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solobj, syncstore->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solsource, syncstore->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols, syncstore->maxnsols) );

      for( j = 0; j < syncstore->maxnsols; ++j )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols[j], syncstore->ninitvars) );
      }

      SCIP_CALL( SCIPtpiInitLock(&(syncstore->syncdata[i].lock)) );
      SCIP_CALL( SCIPtpiInitCondition(&(syncstore->syncdata[i].allsynced)) );
   }

   syncstore->initialized = TRUE;
   syncstore->stopped = FALSE;

   SCIP_CALL( SCIPgetIntParam(scip, "parallel/mode", &paramode) );
   syncstore->mode = (SCIP_PARALLELMODE) paramode;

   SCIP_CALL( SCIPtpiInit(syncstore->nsolvers, INT_MAX, FALSE) );
   SCIP_CALL( SCIPautoselectDisps(scip) );

   if( syncstore->mode == SCIP_PARA_DETERMINISTIC )
   {
      /* in deterministic mode use the number of non-zeros and the number of variables to get a good
       * syncdelay and maximum syncfreq
       */
      syncstore->minsyncdelay *= 0.01 * (SCIPgetNNZs(scip) * SCIPgetNVars(scip));
      syncstore->syncfreqmax *= 0.01 * (SCIPgetNNZs(scip) * SCIPgetNVars(scip));
   }

   return SCIP_OKAY;
}

/** deinitializes the synchronization store */
SCIP_RETCODE SCIPsyncstoreExit(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   int i;
   int j;

   assert(syncstore != NULL);
   assert(syncstore->initialized);

   SCIP_CALL( SCIPtpiExit() );

   for( i = 0; i < syncstore->nsyncdata; ++i )
   {
      SCIPtpiDestroyLock(&(syncstore->syncdata[i].lock));
      SCIPtpiDestroyCondition(&(syncstore->syncdata[i].allsynced));
      SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solobj, syncstore->maxnsols);
      SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solsource, syncstore->maxnsols);
      SCIPboundstoreFree(syncstore->mainscip,  &syncstore->syncdata[i].boundstore);

      for( j = 0; j < syncstore->maxnsols; ++j )
      {
         SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols[j], syncstore->ninitvars);
      }

      SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols, syncstore->maxnsols);
   }

   SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata, syncstore->nsyncdata);

   syncstore->initialized = FALSE;
   syncstore->stopped = FALSE;

   return SCIP_OKAY;
}

/** checks whether the solve-is-stopped flag in the syncstore has been set by any thread */
SCIP_Bool SCIPsyncstoreSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   SCIP_Bool stopped;

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&syncstore->lock) );

   stopped = syncstore->stopped;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&syncstore->lock) );

   return stopped;
}

/** sets the solve-is-stopped flag in the syncstore so that subsequent calls to
 *  SCIPsyncstoreSolveIsStopped will return the given value in any thread
 */
void SCIPsyncstoreSetSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Bool             stopped             /**< flag if the solve is stopped */
   )
{
   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&syncstore->lock) );

   syncstore->stopped = stopped;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&syncstore->lock) );
}

/** gets the upperbound from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastUpperbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? SCIPinfinity(syncstore->mainscip) : syncstore->lastsync->bestupperbound;
}

/** gets the lowerbound from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastLowerbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? -SCIPinfinity(syncstore->mainscip) : syncstore->lastsync->bestlowerbound;
}

/** gets the number of solutions from the last synchronization */
int SCIPsyncstoreGetLastNSols(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0 : syncstore->lastsync->nsols;
}

/** gets the number of boundchanges from the last synchronization */
int SCIPsyncstoreGetLastNBounds(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0 : SCIPboundstoreGetNChgs(syncstore->lastsync->boundstore);
}

/** gets total memory used by all solvers from the last synchronization */
SCIP_Longint SCIPsyncstoreGetLastMemTotal(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0 : syncstore->lastsync->memtotal;
}

/** gets the synchronization frequency from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastSyncfreq(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0.0 : syncstore->lastsync->syncfreq;
}

/** get synchronization data with given number. It is the responsibility of the caller
 *  to only ask for a synchronization number that still exists, which is checked
 *  with an assert in debug mode. */
SCIP_SYNCDATA* SCIPsyncstoreGetSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum             /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   )
{
   int j;

   assert(syncstore != NULL);
   assert(syncstore->initialized);

   j = syncnum % syncstore->nsyncdata;

   /* check if requested syncnumber still exists if in debug mode */
   assert(syncstore->syncdata[j].syncnum == syncnum);

   return &syncstore->syncdata[j];
}

/** get the next synchronization data that should be read and
 *  adjust the delay. Returns NULL if no more data should be read due to minimum delay */
SCIP_SYNCDATA* SCIPsyncstoreGetNextSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq,           /**< the current synchronization frequency */
   SCIP_Longint          writenum,           /**< number of synchronizations the solver has written to */
   SCIP_Real*            delay               /**< pointer holding the current synchronization delay */
   )
{
   SCIP_Real newdelay;
   SCIP_Longint nextsyncnum;

   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(delay != NULL);

   if( syncdata == NULL )
   {
      nextsyncnum = 0;
   }
   else
   {
      if( syncdata->status != SCIP_STATUS_UNKNOWN )
         return NULL;

      nextsyncnum = syncdata->syncnum + 1;
   }

   if( nextsyncnum == writenum )
      return NULL;

   newdelay = *delay - syncfreq;

   /* if the delay would get too small we dont want to read the next syncdata.
    * But due to the limited length of the syncdata array we might need to
    * read this synchronization data anyways which is checked by the second part
    * of the if condition
    */
   if( newdelay < syncstore->minsyncdelay && nextsyncnum >= writenum - syncstore->maxnsyncdelay )
      return NULL;

   *delay = newdelay;
   assert(syncstore->syncdata[nextsyncnum % syncstore->nsyncdata].syncnum == nextsyncnum);

   return &syncstore->syncdata[nextsyncnum % syncstore->nsyncdata];
}

/** ensures that the given synchronization data has been written by
 *  all solvers upon return of this function and blocks the caller if necessary. */
SCIP_RETCODE SCIPsyncstoreEnsureAllSynced(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   /* check if waiting is required, make sure to hold the lock */
   SCIP_CALL( SCIPtpiAcquireLock(&syncdata->lock) );

   while( syncdata->syncedcount < syncstore->nsolvers )
   {
      /* yes, so wait on the condition variable
       * (automatically releases the lock and reacquires it after the waiting)
       */
      SCIP_CALL( SCIPtpiWaitCondition(&syncdata->allsynced, &syncdata->lock) );
   }

   SCIP_CALL( SCIPtpiReleaseLock(&syncdata->lock) );

   return SCIP_OKAY;
}

/** Start synchronization for the given concurrent solver.
 *  Needs to be followed by a call to SCIPsyncstoreFinishSync if
 *  the syncdata that is returned is not NULL
 */
SCIP_RETCODE SCIPsyncstoreStartSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum,            /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   SCIP_SYNCDATA**       syncdata            /**< pointer to return the synchronization data */
   )
{
   int i;

   assert(syncdata != NULL);
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   if( SCIPsyncstoreSolveIsStopped(syncstore) )
   {
      *syncdata = NULL;
      return SCIP_OKAY;
   }

   i = syncnum % syncstore->nsyncdata;
   *syncdata = &syncstore->syncdata[i];
   assert(*syncdata != NULL);

   SCIP_CALL( SCIPtpiAcquireLock(&(*syncdata)->lock) );

   if( (*syncdata)->syncnum != syncnum )
   {
      SCIPboundstoreClear((*syncdata)->boundstore);
      (*syncdata)->nsols = 0;
      (*syncdata)->memtotal = SCIPgetMemTotal(syncstore->mainscip);
      (*syncdata)->syncedcount = 0;
      (*syncdata)->bestupperbound = SCIPinfinity(syncstore->mainscip);
      (*syncdata)->bestlowerbound = -(*syncdata)->bestupperbound;
      (*syncdata)->status = SCIP_STATUS_UNKNOWN;
      (*syncdata)->winner = 0;
      (*syncdata)->syncnum = syncnum;
      (*syncdata)->syncfreq = 0.0;
   }

   return SCIP_OKAY;
}

/** finishes synchronization for the synchronization data */
SCIP_RETCODE SCIPsyncstoreFinishSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA**       syncdata            /**< the synchronization data */
   )
{
   SCIP_Bool printline = FALSE;

   assert(syncdata != NULL);
   assert((*syncdata) != NULL);
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   ++(*syncdata)->syncedcount;

   if( (*syncdata)->syncedcount == syncstore->nsolvers )
   {
      if( (*syncdata)->status != SCIP_STATUS_UNKNOWN )
         SCIPsyncstoreSetSolveIsStopped(syncstore, TRUE);

      syncstore->lastsync = *syncdata;
      printline = TRUE;

      SCIP_CALL( SCIPtpiBroadcastCondition(&(*syncdata)->allsynced) );
   }

   SCIP_CALL( SCIPtpiReleaseLock(&(*syncdata)->lock) );

   if( printline )
   {
      SCIP_CALL( SCIPprintDisplayLine(syncstore->mainscip, NULL, SCIP_VERBLEVEL_HIGH, TRUE) );
   }

   *syncdata = NULL;

   return SCIP_OKAY;
}

/** gets status in synchronization data */
SCIP_STATUS SCIPsyncdataGetStatus(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->status;
}

/** gets the solver that had the best status, or -1 if solve is not stopped yet */
int SCIPsyncstoreGetWinner(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   if( syncstore->lastsync == NULL || syncstore->lastsync->status == SCIP_STATUS_UNKNOWN )
      return -1;

   return syncstore->lastsync->winner;
}

/** how many solvers have already finished synchronizing on this sychronization data */
int SCIPsyncdataGetNSynced(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->syncedcount;
}

/** how many solvers have are running concurrently */
int SCIPsyncstoreGetNSolvers(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->nsolvers;
}


/** read amount total memory used from synchronization data */
SCIP_Longint SCIPsyncdataGetMemTotal(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->memtotal;
}

/** read the synchronization frequency from a synchronization data */
SCIP_Real SCIPsyncdataGetSyncFreq(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->syncfreq;
}

/** read the upperbound stored in a synchronization data */
SCIP_Real SCIPsyncdataGetUpperbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->bestupperbound;
}

/** read the lowerbound stored in a synchronization data */
SCIP_Real SCIPsyncdataGetLowerbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->bestlowerbound;
}

/** read the solutions stored in a synchronization data */
void SCIPsyncdataGetSolutions(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real***          solvalues,          /**< array of buffers containing the solution values */
   int**                 solowner,           /**< array of ownerids of solutions */
   int*                  nsols               /**< pointer to return number of solutions */
   )
{
   assert(syncdata != NULL);
   assert(solvalues != NULL);
   assert(solowner != NULL);
   assert(nsols != NULL);

   *solvalues = syncdata->sols;
   *solowner = syncdata->solsource;
   *nsols = syncdata->nsols;
}

/** read bound changes stored in the synchronization data */
SCIP_BOUNDSTORE* SCIPsyncdataGetBoundChgs(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->boundstore;
}

/** write the synchronization frequency to a synchronization data */
void SCIPsyncdataSetSyncFreq(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq            /**< the synchronization frequency */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);

   syncdata->syncfreq = MIN(syncfreq, syncstore->syncfreqmax);
}

/** set status in the synchronization data */
void SCIPsyncdataSetStatus(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_STATUS           status,             /**< the status */
   int                   solverid            /**< identifier of te solver that has this status */
   )
{
   assert(syncdata != NULL);

   /* check if status is better than current one (closer to SCIP_STATUS_OPTIMAL),
    * break ties by the solverid, and remember the solver wit the best status
    * so that the winner will be selected deterministically
    */
   if( syncdata->status < SCIP_STATUS_OPTIMAL )
   {

      if( status > syncdata->status || (status == syncdata->status && solverid < syncdata->winner) )
      {
         syncdata->status = status;
         syncdata->winner = solverid;
      }
   }
   else if( syncdata->status > SCIP_STATUS_OPTIMAL && status >= SCIP_STATUS_OPTIMAL )
   {
      if( status < syncdata->status || (status == syncdata->status && solverid < syncdata->winner) )
      {
         syncdata->status = status;
         syncdata->winner = solverid;
      }
   }
}

/** adds memory used to the synchronization data */
void SCIPsyncdataAddMemTotal(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Longint          memtotal            /**< the number of bytes used */
   )
{
   assert(syncdata != NULL);

   syncdata->memtotal += memtotal;
}

/** set upperbound to the synchronization data */
void SCIPsyncdataSetUpperbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_Real             upperbound          /**< the upperbound */
   )
{
   assert(syncdata != NULL);

   syncdata->bestupperbound = MIN(syncdata->bestupperbound, upperbound);
}

/** set lowerbound to the synchronization data */
void SCIPsyncdataSetLowerbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the lowerbound should be added to */
   SCIP_Real             lowerbound          /**< the lowerbound */
   )
{
   assert(syncdata != NULL);

   syncdata->bestlowerbound = MAX(syncdata->bestlowerbound, lowerbound);
}

/** gives a buffer to store the solution values, or NULL if solution should not be stored
 *  because there are already better solutions stored.
 */
void SCIPsyncdataGetSolutionBuffer(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Real             solobj,             /**< the objective value of the solution */
   int                   ownerid,            /**< an identifier for the owner of the solution, e.g. the thread number */
   SCIP_Real**           buffer              /**< pointer to return a buffer for the solution values, which must be set
                                              *   if the buffer is not NULL */
   )
{
   int                  pos;
   int                  i;

   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);
   assert(buffer != NULL);

   for( pos = 0; pos < syncdata->nsols; ++pos )
   {
      if( syncdata->solobj[pos] < solobj || (syncdata->solobj[pos] == solobj && ownerid < syncdata->solsource[pos]) ) /*lint !e777*/
         break;
   }

   if( syncdata->nsols < syncstore->maxnsols )
   {
      for( i = syncdata->nsols; i > pos; --i )
      {
         syncdata->solobj[i] = syncdata->solobj[i - 1];
         syncdata->solsource[i] = syncdata->solsource[i - 1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i - 1]);
      }

      ++syncdata->nsols;
   }
   else
   {
      --pos;

      for( i = 0; i < pos; ++i )
      {
         syncdata->solobj[i] = syncdata->solobj[i + 1];
         syncdata->solsource[i] = syncdata->solsource[i + 1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i + 1]);
      }
   }

   if( pos >= 0 )
   {
      syncdata->solobj[pos] = solobj;
      syncdata->solsource[pos] = ownerid;
      *buffer = syncdata->sols[pos];
   }
   else
   {
      *buffer = NULL;
   }
}

/** adds bound changes to the synchronization data */
SCIP_RETCODE SCIPsyncdataAddBoundChanges(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_BOUNDSTORE*      boundstore          /**< bound store containing the bounds to add */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);
   assert(boundstore != NULL);

   SCIP_CALL( SCIPboundstoreMerge(syncstore->mainscip, syncdata->boundstore, boundstore) );

   return SCIP_OKAY;
}

/** is synchronization store initialized */
SCIP_Bool SCIPsyncstoreIsInitialized(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);

   return syncstore->initialized;
}

/** returns the mode of the synchronization store */
SCIP_PARALLELMODE SCIPsyncstoreGetMode(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);

   return syncstore->mode;
}
