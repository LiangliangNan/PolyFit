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

/**@file   syncstore.h
 * @ingroup PARALLEL
 * @brief  the function declarations for the synchronization store
 * @author Robert Lion Gottwald
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SYNCSTORE_H__
#define __SYNCSTORE_H__

#include "scip/def.h"
#include "scip/type_syncstore.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"

/** creates and captures a new synchronization store */
EXTERN
SCIP_RETCODE SCIPsyncstoreCreate(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to return the created synchronization store */
   );

/** releases a synchronization store */
EXTERN
SCIP_RETCODE SCIPsyncstoreRelease(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to the synchronization store */
   );

/** captures a synchronization store */
EXTERN
SCIP_RETCODE SCIPsyncstoreCapture(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** initialize the syncstore for the given SCIP instance */
EXTERN
SCIP_RETCODE SCIPsyncstoreInit(
   SCIP*                 scip                /**< SCIP main datastructure */
   );

/** deinitializes the synchronization store */
EXTERN
SCIP_RETCODE SCIPsyncstoreExit(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** checks whether the solve-is-stopped flag in the syncstore has been set by any thread */
EXTERN
SCIP_Bool SCIPsyncstoreSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** sets the solve-is-stopped flag in the syncstore so that subsequent calls to
 *  SCIPsyncstoreSolveIsStopped will return the given value in any thread
 */
EXTERN
void SCIPsyncstoreSetSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Bool             stopped             /**< flag if the solve is stopped */
   );

/** gets the upperbound from the last synchronization */
EXTERN
SCIP_Real SCIPsyncstoreGetLastUpperbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the lowerbound from the last synchronization */
EXTERN
SCIP_Real SCIPsyncstoreGetLastLowerbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the number of solutions from the last synchronization */
EXTERN
int SCIPsyncstoreGetLastNSols(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the number of boundchanges from the last synchronization */
EXTERN
int SCIPsyncstoreGetLastNBounds(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets total memory used by all solvers from the last synchronization */
EXTERN
SCIP_Longint SCIPsyncstoreGetLastMemTotal(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the synchronization frequency from the last synchronization */
EXTERN
SCIP_Real SCIPsyncstoreGetLastSyncfreq(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** get synchronization data with given number. It is the responsibility of the caller
 *  to only ask for a synchronization number that still exists. */
EXTERN
SCIP_SYNCDATA* SCIPsyncstoreGetSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum             /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   );

/** get the next synchronization data that should be read and
 *  adjust the delay. Returns NULL if no more data should be read due to minimum delay */
EXTERN
SCIP_SYNCDATA* SCIPsyncstoreGetNextSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq,           /**< the current synchronization frequency */
   SCIP_Longint          writenum,           /**< number of synchronizations the solver has written to */
   SCIP_Real*            delay               /**< pointer holding the current synchronization delay */
   );

/** ensures that the given synchronization data has been written by
 *  all solvers upon return of this function and blocks the caller if necessary. */
EXTERN
SCIP_RETCODE SCIPsyncstoreEnsureAllSynced(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** Start synchronization for the given concurrent solver.
 *  Needs to be followed by a call to SCIPsyncstoreFinishSync if
 *  the syncdata that is returned is not NULL
 */
EXTERN
SCIP_RETCODE SCIPsyncstoreStartSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum,            /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   SCIP_SYNCDATA**       syncdata            /**< pointer to return the synchronization data */
   );

/** finishes synchronization for the synchronization data */
EXTERN
SCIP_RETCODE SCIPsyncstoreFinishSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA**       syncdata            /**< the synchronization data */
   );

/** gets status in synchronization data */
EXTERN
SCIP_STATUS SCIPsyncdataGetStatus(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** gets the solver that had the best status, or -1 if solve is not stopped yet */
EXTERN
int SCIPsyncstoreGetWinner(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** how many solvers have already finished synchronizing on this sychronization data */
EXTERN
int SCIPsyncdataGetNSynced(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** how many solvers have are running concurrently */
EXTERN
int SCIPsyncstoreGetNSolvers(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** read amount of memory used from synchronization data */
EXTERN
SCIP_Longint SCIPsyncdataGetMemTotal(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the synchronization frequency from a synchronization data */
EXTERN
SCIP_Real SCIPsyncdataGetSyncFreq(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the upperbound stored in a synchronization data */
EXTERN
SCIP_Real SCIPsyncdataGetUpperbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the lowerbound stored in a synchronization data */
EXTERN
SCIP_Real SCIPsyncdataGetLowerbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the solutions stored in a synchronization data */
EXTERN
void SCIPsyncdataGetSolutions(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real***          solvalues,          /**< pointer to return array of buffers containing the solution values */
   int**                 solowner,           /**< pointer to return array of ownerids of solutions */
   int*                  nsols               /**< pointer to return number of solutions */
   );

/** read bound changes stored in the synchronization data */
EXTERN
SCIP_BOUNDSTORE* SCIPsyncdataGetBoundChgs(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** write the synchronization frequency to a synchronization data */
EXTERN
void SCIPsyncdataSetSyncFreq(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq            /**< the synchronization frequency */
   );

/** set status in the synchronization data */
EXTERN
void SCIPsyncdataSetStatus(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_STATUS           status,             /**< the status */
   int                   solverid            /**< identifier of te solver that has this status */
   );

/** adds memory used to the synchronization data */
EXTERN
void SCIPsyncdataAddMemTotal(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Longint          memtotal            /**< the number of bytes used */
   );

/** set upperbound to the synchronization data */
EXTERN
void SCIPsyncdataSetUpperbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_Real             upperbound          /**< the upperbound */
   );

/** set lowerbound to the synchronization data */
EXTERN
void SCIPsyncdataSetLowerbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the lowerbound should be added to */
   SCIP_Real             lowerbound          /**< the lowerbound */
   );

/** gives a buffer to store the solution values, or NULL if solution should not be stored
 *  because there are already better solutions stored.
 */
EXTERN
void SCIPsyncdataGetSolutionBuffer(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Real             solobj,             /**< the objective value of the solution */
   int                   ownerid,            /**< an identifier for the owner of the solution, e.g. the thread number */
   SCIP_Real**           buffer              /**< pointer to return a buffer for the solution values, which must be set
                                              *   if the buffer is not NULL */
   );

/** adds bound changes to the synchronization data */
EXTERN
SCIP_RETCODE SCIPsyncdataAddBoundChanges(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_BOUNDSTORE*      boundstore          /**< bound store containing the bounds to add */
   );

/** is synchronization store initialized */
EXTERN
SCIP_Bool SCIPsyncstoreIsInitialized(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** returns the mode of the synchronization store */
EXTERN
SCIP_PARALLELMODE SCIPsyncstoreGetMode(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/**@} */
#endif
