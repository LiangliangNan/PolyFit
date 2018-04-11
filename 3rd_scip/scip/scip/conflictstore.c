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

/**@file   conflictstore.c
 * @brief  methods for storing conflicts
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/event.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/prob.h"
#include "scip/reopt.h"
#include "scip/scip.h"
#include "scip/def.h"
#include "scip/cons_linear.h"
#include "scip/struct_conflictstore.h"


#define CONFLICTSTORE_DUALRAYSIZE  100 /* default size of conflict store */
#define CONFLICTSTORE_DUALSOLSIZE   75 /* default size of conflict store */
#define CONFLICTSTORE_MINSIZE     2000 /* default minimal size of a dynamic conflict store */
#define CONFLICTSTORE_MAXSIZE    60000 /* maximal size of a dynamic conflict store (multiplied by 3) */
#define CONFLICTSTORE_SIZE       10000 /* default size of conflict store */
#define CONFLICTSTORE_SORTFREQ      20 /* frequency to resort the conflict array */

/* event handler properties */
#define EVENTHDLR_NAME         "ConflictStore"
#define EVENTHDLR_DESC         "Solution event handler for conflict store."


/* exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecConflictstore)
{/*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_BESTSOLFOUND);

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( SCIPclearConflictStore(scip, event) );
   }

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolConflictstore)
{
   SCIP_Bool cleanboundexceeding;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL( SCIPgetBoolParam(scip, "conflict/cleanboundexceedings", &cleanboundexceeding) );

   if( !cleanboundexceeding )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolConflictstore)
{
   SCIP_Bool cleanboundexceeding;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL( SCIPgetBoolParam(scip, "conflict/cleanboundexceedings", &cleanboundexceeding) );

   if( !cleanboundexceeding )
      return SCIP_OKAY;

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/* comparison method for constraints */
static
SCIP_DECL_SORTPTRCOMP(compareConss)
{
   /*lint --e{715}*/
   SCIP_CONS* cons1 = (SCIP_CONS*)elem1;
   SCIP_CONS* cons2 = (SCIP_CONS*)elem2;

   assert(cons1 != NULL);
   assert(cons2 != NULL);

   if( SCIPconsGetAge(cons1) > SCIPconsGetAge(cons2) + 1e-09 )
      return -1;
   else if ( SCIPconsGetAge(cons1) < SCIPconsGetAge(cons2) - 1e-09 )
      return +1;
   else
#ifdef SCIP_DISABLED_CODE
   /* @todo if both constraints have the same age we prefere the constraint with more non-zeros
    *       this requires a larges change of the callback, passing void-pointer (i.e. a scip
    *       object) would necessary.
    */
   {
      SCIP_Bool success;
      int nvars1;
      int nvars2;

      SCIP_CALL( SCIPgetConsNVars(scip, cons1, &nvars1, &success) );
      assert(success)

      SCIP_CALL( SCIPgetConsNVars(scip, cons2, &nvars2, &success) );
      assert(success)

      if( nvars1 >= nvars2 )
         return -1;
      else
         return +1;
   }
#else
   {
      SCIP_CONSHDLR* conshdlr1 = SCIPconsGetHdlr(cons1);
      SCIP_CONSHDLR* conshdlr2 = SCIPconsGetHdlr(cons2);

      if( strcmp(SCIPconshdlrGetName(conshdlr1), "linear") == strcmp(SCIPconshdlrGetName(conshdlr2), "linear") )
         return 0;
      else if( strcmp(SCIPconshdlrGetName(conshdlr1), "linear") == 0 )
         return -1;
      else
         return +1;
   }
#endif
}

/* initializes the conflict store */
static
SCIP_RETCODE initConflictstore(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob           /**< transformed problem */
   )
{
   assert(conflictstore != NULL);

   /* calculate the maximal size of the conflict store */
   if( conflictstore->maxstoresize == -1 )
   {
      SCIP_CALL( SCIPsetGetIntParam(set, "conflict/maxstoresize", &conflictstore->maxstoresize) );

      /* the size should be dynamic wrt number of variables after presolving */
      if( conflictstore->maxstoresize == -1 )
      {
         int nconss;
         int nvars;

         nconss = SCIPprobGetNConss(transprob);
         nvars = SCIPprobGetNVars(transprob);

         conflictstore->initstoresize = CONFLICTSTORE_MINSIZE;
         conflictstore->initstoresize += 2*nconss;

         if( nvars/2 <= 500 )
            conflictstore->initstoresize += (int) CONFLICTSTORE_MAXSIZE/100;
         else if( nvars/2 <= 5000 )
            conflictstore->initstoresize += (int) CONFLICTSTORE_MAXSIZE/10;
         else
            conflictstore->initstoresize += CONFLICTSTORE_MAXSIZE/2;

         conflictstore->initstoresize = MIN(conflictstore->initstoresize, CONFLICTSTORE_MAXSIZE);
         conflictstore->storesize = conflictstore->initstoresize;
         conflictstore->maxstoresize = (int)(MIN(3.0 * conflictstore->initstoresize, CONFLICTSTORE_MAXSIZE));
      }
      else
      {
         conflictstore->initstoresize = conflictstore->maxstoresize;
         conflictstore->storesize = conflictstore->maxstoresize;
      }

#ifdef NDEBUG
      if( conflictstore->maxstoresize == 0 )
         SCIPsetDebugMsg(set, "usage of conflict pool is disabled.\n");
      else
         SCIPsetDebugMsg(set, "[init,max] size of conflict pool is [%d,%d].\n",
               conflictstore->initstoresize, conflictstore->maxstoresize);
#endif
   }

   return SCIP_OKAY;
}

/** resizes conflict and primal bound arrays to be able to store at least num entries */
static
SCIP_RETCODE conflictstoreEnsureMem(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflictstore != NULL);
   assert(set != NULL);

   /* we do not allocate more memory as allowed */
   if( conflictstore->conflictsize == conflictstore->maxstoresize )
      return SCIP_OKAY;

   if( num > conflictstore->conflictsize )
   {
      int newsize;
#ifndef NDEBUG
      int i;
#endif
      /* initialize the complete data structure */
      if( conflictstore->conflictsize == 0 )
      {
         assert(conflictstore->storesize > 0);

         newsize = MIN(conflictstore->storesize, CONFLICTSTORE_SIZE);
         newsize = MAX(newsize, num);
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->conflicts, newsize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->confprimalbnds, newsize) );
      }
      else
      {
         newsize = SCIPsetCalcMemGrowSize(set, num);
         newsize = MIN(conflictstore->maxstoresize, newsize);
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictstore->conflicts, conflictstore->conflictsize, \
               newsize) );
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictstore->confprimalbnds, conflictstore->conflictsize, \
               newsize) );
      }

#ifndef NDEBUG
      for( i = conflictstore->nconflicts; i < newsize; i++ )
      {
         conflictstore->conflicts[i] = NULL;
         conflictstore->confprimalbnds[i] = -SCIPsetInfinity(set);
      }
#endif
      conflictstore->conflictsize = newsize;
   }
   assert(num <= conflictstore->conflictsize || conflictstore->conflictsize == conflictstore->maxstoresize);

   return SCIP_OKAY;
}

/* increase the dynamic storage if we could not delete enough conflicts
 *
 * we want to have at least set->conf_maxconss free slots in the conflict array, because this is the maximal number
 * of conflicts generated at a node. we increase the size by the minimum of set->conf_maxconss and 1% of the current
 * store size. nevertheless, we don't exceed conflictstore->maxstoresize.
 */
static
void adjustStorageSize(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflictstore != NULL);

   /* increase storage */
   if( conflictstore->storesize - conflictstore->nconflicts <= set->conf_maxconss
      && conflictstore->storesize < conflictstore->maxstoresize )
   {
      SCIP_Real increase = ceil(0.01 * conflictstore->storesize);
      conflictstore->storesize += MIN(set->conf_maxconss, (int)(increase));
      conflictstore->storesize = MIN(conflictstore->storesize, conflictstore->maxstoresize);
   }

   return;
}

/* removes conflict at position pos */
static
SCIP_RETCODE delPosConflict(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL if delete = FALSE */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int                   pos,                /**< position to remove */
   SCIP_Bool             deleteconflict      /**< should the conflict be deleted? */
   )
{
   SCIP_CONS* conflict;
   int lastpos;

   assert(conflictstore != NULL);
   assert(pos >= 0 && pos < conflictstore->nconflicts);

   lastpos = conflictstore->nconflicts-1;
   conflict = conflictstore->conflicts[pos];
   assert(conflict != NULL);

   /* decrease number of conflicts depending an a cutoff bound */
   conflictstore->ncbconflicts -= (SCIPsetIsInfinity(set, REALABS(conflictstore->confprimalbnds[pos])) ? 0 : 1);

#ifdef SCIP_PRINT_DETAILS
   SCIPsetDebugMsg(set, "-> remove conflict at pos=%d with age=%g\n", pos, SCIPconsGetAge(conflict));
#endif

   /* mark the constraint as deleted */
   if( deleteconflict && !SCIPconsIsDeleted(conflict) )
   {
      assert(transprob != NULL);
      SCIP_CALL( SCIPconsDelete(conflictstore->conflicts[pos], blkmem, set, stat, transprob, reopt) );
   }
   SCIP_CALL( SCIPconsRelease(&conflictstore->conflicts[pos], blkmem, set) );

   /* replace with conflict at the last position */
   if( pos < lastpos )
   {
      conflictstore->conflicts[pos] = conflictstore->conflicts[lastpos];
      conflictstore->confprimalbnds[pos] = conflictstore->confprimalbnds[lastpos];
   }

#ifndef NDEBUG
   conflictstore->conflicts[lastpos] = NULL;
   conflictstore->confprimalbnds[lastpos] = -SCIPsetInfinity(set);
#endif

   /* decrease number of conflicts */
   --conflictstore->nconflicts;

   return SCIP_OKAY;
}

/* removes proof based on a dual ray at position pos */
static
SCIP_RETCODE delPosDualray(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL if delete = FALSE */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int                   pos,                /**< position to remove */
   SCIP_Bool             deleteconflict      /**< should the dual ray be deleted? */
   )
{
   SCIP_CONS* dualproof;
   SCIP_Bool success;
   int lastpos;
   int nvars;

   assert(conflictstore != NULL);

   lastpos = conflictstore->ndualrayconfs-1;
   dualproof = conflictstore->dualrayconfs[pos];
   assert(dualproof != NULL);

   /* decrease the number of non-zeros */
   SCIP_CALL( SCIPconsGetNVars(dualproof, set, &nvars, &success) );
   assert(success);
   conflictstore->nnzdualrays -= nvars;

#ifdef SCIP_PRINT_DETAILS
   SCIPsetDebugMsg(set, "-> remove dual proof (ray) at pos=%d age=%g nvars=%d\n", pos, SCIPconsGetAge(dualproof), nvars);
#endif

   /* mark the constraint as deleted */
   if( deleteconflict && !SCIPconsIsDeleted(dualproof) )
   {
      assert(transprob != NULL);
      SCIP_CALL( SCIPconsDelete(dualproof, blkmem, set, stat, transprob, reopt) );
   }
   SCIP_CALL( SCIPconsRelease(&dualproof, blkmem, set) );

   /* replace with dual ray at the last position */
   if( pos < lastpos )
   {
      conflictstore->dualrayconfs[pos] = conflictstore->dualrayconfs[lastpos];

#ifndef NDEBUG
      conflictstore->dualrayconfs[lastpos] = NULL;
#endif
   }

   /* decrease number of dual rays */
   --conflictstore->ndualrayconfs;

   return SCIP_OKAY;
}

/* removes proof based on a dual solution at position pos */
static
SCIP_RETCODE delPosDualsol(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem, or NULL if delete = FALSE */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int                   pos,                /**< position to remove */
   SCIP_Bool             deleteconflict      /**< should the dual ray be deleted? */
   )
{
   SCIP_CONS* dualproof;
   SCIP_Bool success;
   int lastpos;
   int nvars;

   assert(conflictstore != NULL);
   assert(pos >= 0 && pos < conflictstore->ndualsolconfs);

   lastpos = conflictstore->ndualsolconfs-1;
   dualproof = conflictstore->dualsolconfs[pos];
   assert(dualproof != NULL);

   /* decrease the number of non-zeros */
   SCIP_CALL( SCIPconsGetNVars(dualproof, set, &nvars, &success) );
   assert(success);
   conflictstore->nnzdualsols -= nvars;

#ifdef SCIP_PRINT_DETAILS
   SCIPsetDebugMsg(set, "-> remove dual proof (sol) at pos=%d age=%g nvars=%d\n", pos, SCIPconsGetAge(dualproof), nvars);
#endif

   /* mark the constraint as deleted */
   if( deleteconflict && !SCIPconsIsDeleted(dualproof) )
   {
      assert(transprob != NULL);
      SCIP_CALL( SCIPconsDelete(dualproof, blkmem, set, stat, transprob, reopt) );
   }
   SCIP_CALL( SCIPconsRelease(&dualproof, blkmem, set) );

   /* replace with dual ray at the last position */
   if( pos < lastpos )
   {
      conflictstore->dualsolconfs[pos] = conflictstore->dualsolconfs[lastpos];
      conflictstore->dualprimalbnds[pos] = conflictstore->dualprimalbnds[lastpos];
      conflictstore->scalefactors[pos] = conflictstore->scalefactors[lastpos];
      conflictstore->updateside[pos] = conflictstore->updateside[lastpos];

#ifndef NDEBUG
      conflictstore->dualsolconfs[lastpos] = NULL;
      conflictstore->dualprimalbnds[lastpos] = SCIP_UNKNOWN;
      conflictstore->scalefactors[lastpos] = 1.0;
      conflictstore->updateside[lastpos] = FALSE;
#endif
   }

   /* decrease number of dual rays */
   --conflictstore->ndualsolconfs;

   return SCIP_OKAY;
}

/** removes all deleted conflicts from the storage */
static
SCIP_RETCODE cleanDeletedConflicts(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int*                  ndelconfs           /**< pointer to store the number of deleted conflicts */
   )
{
   int i;

   assert(conflictstore != NULL);

   (*ndelconfs) = 0;

   /* we traverse backwards to avoid swapping of pointers */
   for( i = conflictstore->nconflicts-1; i >= 0; i-- )
   {
      assert(conflictstore->conflicts[i] != NULL);

      /* check whether the constraint is already marked as deleted */
      if( SCIPconsIsDeleted(conflictstore->conflicts[i]) )
      {
         /* remove conflict at current position */
         SCIP_CALL( delPosConflict(conflictstore, set, stat, NULL, blkmem, reopt, i, FALSE) );

         ++(*ndelconfs);
      }
   }
   SCIPsetDebugMsg(set, "> removed %d/%d as deleted marked conflicts.\n", *ndelconfs, conflictstore->nconflicts + (*ndelconfs));

   return SCIP_OKAY;
}

/** removes all deleted dual proofs of infeasible LP relaxations from the storage */
static
SCIP_RETCODE cleanDeletedDualrayCons(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int*                  ndelproofs          /**< pointer to store the number of deleted conflicts */
   )
{
   int i;

   assert(conflictstore != NULL);

   (*ndelproofs) = 0;

   /* we traverse backwards to avoid swapping of pointers */
   for( i = conflictstore->ndualrayconfs-1; i >= 0; i-- )
   {
      assert(conflictstore->dualrayconfs[i] != NULL);

      /* check whether the constraint is already marked as deleted */
      if( SCIPconsIsDeleted(conflictstore->dualrayconfs[i]) )
      {
         /* remove proof at current position */
         SCIP_CALL( delPosDualray(conflictstore, set, stat, NULL, blkmem, reopt, i, FALSE) );

         ++(*ndelproofs);
      }
   }

   SCIPsetDebugMsg(set, "> removed %d/%d as deleted marked dual infeasibility proofs.\n",
         *ndelproofs, conflictstore->ndualrayconfs + (*ndelproofs));

   return SCIP_OKAY;
}

/** removes all deleted dual proofs of bound exceeding LP relaxations from the storage */
static
SCIP_RETCODE cleanDeletedDualsolCons(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int*                  ndelproofs          /**< pointer to store the number of deleted conflicts */
   )
{
   int i;

   assert(conflictstore != NULL);

   (*ndelproofs) = 0;

   /* we traverse backwards to avoid swapping of pointers */
   for( i = conflictstore->ndualsolconfs-1; i >= 0; i-- )
   {
      assert(conflictstore->dualsolconfs[i] != NULL);

      /* check whether the constraint is already marked as deleted */
      if( SCIPconsIsDeleted(conflictstore->dualsolconfs[i]) )
      {
         /* remove proof at current position */
         SCIP_CALL( delPosDualsol(conflictstore, set, stat, NULL, blkmem, reopt, i, FALSE) );

         ++(*ndelproofs);
      }
   }

   SCIPsetDebugMsg(set, "> removed %d/%d as deleted marked dual boundexceeding proofs.\n",
         *ndelproofs, conflictstore->ndualrayconfs + (*ndelproofs));

   return SCIP_OKAY;
}

/** cleans up the storage */
static
SCIP_RETCODE conflictstoreCleanUpStorage(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   int ndelconfs;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);

   /* the storage is empty  */
   if( conflictstore->nconflicts == 0 )
      return SCIP_OKAY;
   assert(conflictstore->nconflicts >= 1);

   ndelconfs = 0;

   /* remove all as deleted marked conflicts */
   SCIP_CALL( cleanDeletedConflicts(conflictstore, set, stat, blkmem, reopt, &ndelconfs) );

   /* return if at least one conflict could be deleted */
   if( ndelconfs > 0 )
      goto TERMINATE;

   /* only clean up the storage if it is filled enough */
   if( conflictstore->nconflicts < conflictstore->conflictsize )
      goto TERMINATE;

   /* resort the array regularly */
   if( conflictstore->ncleanups % CONFLICTSTORE_SORTFREQ == 0 )
   {
      /* sort conflict */
      SCIPsortPtrReal((void**)conflictstore->conflicts, conflictstore->confprimalbnds, compareConss, conflictstore->nconflicts);
      assert(SCIPsetIsGE(set, SCIPconsGetAge(conflictstore->conflicts[0]),
            SCIPconsGetAge(conflictstore->conflicts[conflictstore->nconflicts-1])));
   }
   assert(conflictstore->nconflicts > 0);

   if( conflictstore->ncleanups % CONFLICTSTORE_SORTFREQ == 0 )
   {
      /* remove conflict at first position (array is sorted) */
      SCIP_CALL( delPosConflict(conflictstore, set, stat, transprob, blkmem, reopt, 0, TRUE) );
   }
   else
   {
      SCIP_Real maxage;
      int oldest_i;
      int i;

      assert(!SCIPconsIsDeleted(conflictstore->conflicts[0]));

      maxage = SCIPconsGetAge(conflictstore->conflicts[0]);
      oldest_i = 0;

      /* check the first 10% of conflicts and find the oldest */
      for( i = 1; i < 0.1 * conflictstore->nconflicts; i++ )
      {
         assert(!SCIPconsIsDeleted(conflictstore->conflicts[i]));

         if( SCIPconsGetAge(conflictstore->conflicts[i]) > maxage )
         {
            maxage = SCIPconsGetAge(conflictstore->conflicts[i]);
            oldest_i = i;
         }
      }

      /* remove conflict at position oldest_i */
      SCIP_CALL( delPosConflict(conflictstore, set, stat, transprob, blkmem, reopt, oldest_i, TRUE) );
   }
   ++ndelconfs;

   /* adjust size of the storage if we use a dynamic store */
   if( set->conf_maxstoresize == -1 )
      adjustStorageSize(conflictstore, set);
   assert(conflictstore->initstoresize <= conflictstore->storesize);
   assert(conflictstore->storesize <= conflictstore->maxstoresize);

  TERMINATE:

   /* increase the number of clean ups */
   ++conflictstore->ncleanups;

   SCIPsetDebugMsg(set, "clean-up #%lld: removed %d/%d conflicts, %d depending on cutoff bound\n",
         conflictstore->ncleanups, ndelconfs, conflictstore->nconflicts+ndelconfs, conflictstore->ncbconflicts);

   return SCIP_OKAY;
}

/** adds an original conflict constraint to the store
 *
 *  @note the constraint will be only transfered to the storage of the transformed problem after calling
 *        SCIPconflictstoreTransform()
 */
static
SCIP_RETCODE conflictstoreAddOrigConflict(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS*            cons                /**< conflict constraint */
   )
{
   assert(conflictstore != NULL);
   assert(cons != NULL);

   if( conflictstore->origconfs == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->origconfs, CONFLICTSTORE_MINSIZE) );
      conflictstore->origconflictsize = CONFLICTSTORE_MINSIZE;
   }
   else if( conflictstore->norigconfs == conflictstore->origconflictsize )
   {
      int newsize = SCIPsetCalcMemGrowSize(set, conflictstore->origconflictsize+1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictstore->origconfs, conflictstore->origconflictsize, newsize) );
      conflictstore->origconflictsize = newsize;
   }

   SCIPconsCapture(cons);
   conflictstore->origconfs[conflictstore->norigconfs] = cons;
   ++conflictstore->norigconfs;

   return SCIP_OKAY;
}


/** creates conflict store */
SCIP_RETCODE SCIPconflictstoreCreate(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict store */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflictstore != NULL);

   SCIP_ALLOC( BMSallocMemory(conflictstore) );

   (*conflictstore)->conflicts = NULL;
   (*conflictstore)->confprimalbnds = NULL;
   (*conflictstore)->dualprimalbnds = NULL;
   (*conflictstore)->scalefactors = NULL;
   (*conflictstore)->updateside = NULL;
   (*conflictstore)->dualrayconfs = NULL;
   (*conflictstore)->dualsolconfs = NULL;
   (*conflictstore)->origconfs = NULL;
   (*conflictstore)->nnzdualrays = 0;
   (*conflictstore)->nnzdualsols = 0;
   (*conflictstore)->conflictsize = 0;
   (*conflictstore)->origconflictsize = 0;
   (*conflictstore)->nconflicts = 0;
   (*conflictstore)->ndualrayconfs = 0;
   (*conflictstore)->ndualsolconfs = 0;
   (*conflictstore)->norigconfs = 0;
   (*conflictstore)->ncbconflicts = 0;
   (*conflictstore)->nconflictsfound = 0;
   (*conflictstore)->initstoresize = -1;
   (*conflictstore)->storesize = -1;
   (*conflictstore)->maxstoresize = -1;
   (*conflictstore)->ncleanups = 0;
   (*conflictstore)->lastcutoffbound = SCIP_INVALID;
   (*conflictstore)->lastnodenum = -1;
   (*conflictstore)->eventhdlr = SCIPsetFindEventhdlr(set, EVENTHDLR_NAME);

   /* create event handler for LP events */
   if( (*conflictstore)->eventhdlr == NULL )
   {
      SCIP_CALL( SCIPeventhdlrCreate(&(*conflictstore)->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, NULL, NULL,
            NULL, NULL, eventInitsolConflictstore, eventExitsolConflictstore, NULL, eventExecConflictstore, NULL) );
      SCIP_CALL( SCIPsetIncludeEventhdlr(set, (*conflictstore)->eventhdlr) );
   }
   assert((*conflictstore)->eventhdlr != NULL);

   return SCIP_OKAY;
}

/** frees conflict store */
SCIP_RETCODE SCIPconflictstoreFree(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   assert(conflictstore != NULL);
   assert(*conflictstore != NULL);

   /* clear the storage */
   SCIP_CALL( SCIPconflictstoreClear(*conflictstore, blkmem, set, stat, reopt) );

   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->origconfs, (*conflictstore)->origconflictsize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->conflicts, (*conflictstore)->conflictsize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->confprimalbnds, (*conflictstore)->conflictsize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->dualrayconfs, CONFLICTSTORE_DUALRAYSIZE);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->dualsolconfs, CONFLICTSTORE_DUALSOLSIZE);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->dualprimalbnds, CONFLICTSTORE_DUALSOLSIZE);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->scalefactors, CONFLICTSTORE_DUALSOLSIZE);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictstore)->updateside, CONFLICTSTORE_DUALSOLSIZE);
   BMSfreeMemoryNull(conflictstore);

   return SCIP_OKAY;
}

/** clears conflict store */
SCIP_RETCODE SCIPconflictstoreClear(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   int i;

   assert(conflictstore != NULL);

   SCIPsetDebugMsg(set, "clearing conflict store: %d origconfs, %d conflicts, %d dual proofs\n",
         conflictstore->norigconfs, conflictstore->nconflicts, conflictstore->ndualrayconfs + conflictstore->ndualsolconfs);

   /* remove original constraints (if present) */
   if( conflictstore->origconfs != NULL )
   {
      for( i = 0; i < conflictstore->norigconfs; i++ )
      {
         SCIP_CONS* conflict = conflictstore->origconfs[i];
         SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );
      }
      conflictstore->norigconfs = 0;
   }

   /* clean up storage of conflict constraints */
   if( conflictstore->conflicts != NULL )
   {
      /* we traverse in reverse order to avoid swapping of pointers */
      for( i = conflictstore->nconflicts-1; i >= 0; i--)
      {
         SCIP_CALL( delPosConflict(conflictstore, set, stat, NULL, blkmem, reopt, i, FALSE) );
      }
      assert(conflictstore->nconflicts == 0);
   }

   /* clean up storage of proof constraints based on dual rays */
   if( conflictstore->dualrayconfs != NULL )
   {
      /* we traverse in reverse order to avoid swapping of pointers */
      for( i = conflictstore->ndualrayconfs-1; i >= 0; i-- )
      {
         SCIP_CALL( delPosDualray(conflictstore, set, stat, NULL, blkmem, reopt, i, FALSE) );
      }
      assert(conflictstore->ndualrayconfs == 0);
   }

   /* clean up storage of proof constraints based on dual solutions */
   if( conflictstore->dualsolconfs != NULL )
   {
      /* we traverse in reverse order to avoid swapping of pointers */
      for( i = conflictstore->ndualsolconfs-1; i >= 0; i-- )
      {
         SCIP_CALL( delPosDualsol(conflictstore, set, stat, NULL, blkmem, reopt, i, FALSE) );
      }
      assert(conflictstore->ndualsolconfs == 0);
   }

   return SCIP_OKAY;
}

/** cleans up conflict store */
SCIP_RETCODE SCIPconflictstoreClean(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   int ndelconfs;
   int ndeldualray;
   int ndeldualsol;

   assert(conflictstore != NULL);

   SCIPsetDebugMsg(set, "cleaning conflict store: %d conflicts, %d dual proofs\n",
         conflictstore->norigconfs, conflictstore->nconflicts, conflictstore->ndualrayconfs + conflictstore->ndualsolconfs);

   ndelconfs = 0;
   ndeldualray = 0;
   ndeldualsol = 0;

   /* remove all as deleted marked conflicts */
   SCIP_CALL( cleanDeletedConflicts(conflictstore, set, stat, blkmem, reopt, &ndelconfs) );

   /* remove all as deleted marked dual infeasibility proofs */
   SCIP_CALL( cleanDeletedDualrayCons(conflictstore, set, stat, blkmem, reopt, &ndeldualray) );

   /* remove all as deleted marked dual bound exceeding proofs */
   SCIP_CALL( cleanDeletedDualsolCons(conflictstore, set, stat, blkmem, reopt, &ndeldualsol) );

   /* don't update bound exceeding proofs after a restart
    *
    * TODO: check whether we want to delete bound exceeding proofs in general during a restart
    */
   if( SCIPisInRestart(set->scip) )
   {
      int i;
      for( i = 0; i < conflictstore->ndualsolconfs; i++ )
         conflictstore->updateside[i] = FALSE;
   }

   return SCIP_OKAY;
}

/** adds a constraint to the pool of proof constraints based on dual rays
 *
 *  @note this methods captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddDualraycons(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_CONS*            dualproof,          /**< constraint based on a dual ray */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   int nvars;
   SCIP_Bool success;

   assert(conflictstore != NULL);
   assert(conflictstore->ndualrayconfs <= CONFLICTSTORE_DUALRAYSIZE);

   /* mark the constraint to be a conflict */
   SCIPconsMarkConflict(dualproof);

   /* create an array to store constraints based on dual rays */
   if( conflictstore->dualrayconfs == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->dualrayconfs, CONFLICTSTORE_DUALRAYSIZE) );
   }

   /* the store is full, we proceed as follows
    *
    * 1. check whether some constraints are marked as deleted and remove those
    * 2. if no constraint is marked as deleted: remove the oldest
    */
   if( conflictstore->ndualrayconfs == CONFLICTSTORE_DUALRAYSIZE )
   {
      int ndeleted = 0;

      /* remove all as deleted marked dual infeasibility proofs */
      SCIP_CALL( cleanDeletedDualrayCons(conflictstore, set, stat, blkmem, reopt, &ndeleted) );

      /* if we could not remove a dual ray that is already marked as deleted we need to remove the oldest active one */
      if( ndeleted == 0 )
      {
         /* sort dual rays */
         SCIPsortPtr((void**)conflictstore->dualrayconfs, compareConss, conflictstore->ndualrayconfs);
         assert(SCIPsetIsGE(set, SCIPconsGetAge(conflictstore->dualrayconfs[0]),
               SCIPconsGetAge(conflictstore->dualrayconfs[conflictstore->ndualrayconfs-1])));

         SCIP_CALL( delPosDualray(conflictstore, set, stat, transprob, blkmem, reopt, 0, TRUE) );
      }
   }

   /* add the new constraint based on a dual ray at the last position */
   SCIPconsCapture(dualproof);
   conflictstore->dualrayconfs[conflictstore->ndualrayconfs] = dualproof;
   ++conflictstore->ndualrayconfs;

   /* increase the number of non-zeros */
   SCIP_CALL( SCIPconsGetNVars(dualproof, set, &nvars, &success) );
   assert(success);
   conflictstore->nnzdualrays += nvars;

   return SCIP_OKAY;
}

/** adds a constraint to the pool of proof constraints based on dual solutions
 *
 *  @note this methods captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddDualsolcons(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_CONS*            dualproof,          /**< constraint based on a dual solution */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_Real             scale,              /**< scaling factor that needs to be considered when updating the side */
   SCIP_Bool             updateside          /**< should the side be updated if a new incumbent is found */
   )
{
   int nvars;
   SCIP_Bool success;

   assert(conflictstore != NULL);
   assert(conflictstore->ndualsolconfs <= CONFLICTSTORE_DUALSOLSIZE);

   /* mark the constraint to be a conflict */
   SCIPconsMarkConflict(dualproof);

   /* create an array to store constraints based on dual rays */
   if( conflictstore->dualsolconfs == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->dualsolconfs, CONFLICTSTORE_DUALSOLSIZE) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->dualprimalbnds, CONFLICTSTORE_DUALSOLSIZE) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->scalefactors, CONFLICTSTORE_DUALSOLSIZE) );
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->updateside, CONFLICTSTORE_DUALSOLSIZE) );
   }

   /* the store is full, we proceed as follows
    *
    * 1. check whether some constraints are marked as deleted and remove those
    * 2. if no constraint is marked as deleted: remove the oldest
    */
   if( conflictstore->ndualsolconfs == CONFLICTSTORE_DUALSOLSIZE )
   {
      int ndeleted = 0;

      /* remove all as deleted marked dual bound exceeding proofs */
      SCIP_CALL( cleanDeletedDualsolCons(conflictstore, set, stat, blkmem, reopt, &ndeleted) );

      /* if we could not remove a dual proof that is already marked as deleted we need to remove the oldest active one */
      if( ndeleted == 0 )
      {
         /* sort dual rays */
         SCIPsortPtrRealRealInt((void**)conflictstore->dualsolconfs, conflictstore->dualprimalbnds,
               conflictstore->scalefactors, (int*)conflictstore->updateside, compareConss, conflictstore->ndualsolconfs);
         assert(SCIPsetIsGE(set, SCIPconsGetAge(conflictstore->dualsolconfs[0]),
               SCIPconsGetAge(conflictstore->dualsolconfs[conflictstore->ndualsolconfs-1])));

         SCIP_CALL( delPosDualsol(conflictstore, set, stat, transprob, blkmem, reopt, 0, TRUE) );
      }
   }

   /* add the new constraint based on a dual solution at the last position */
   SCIPconsCapture(dualproof);
   conflictstore->dualsolconfs[conflictstore->ndualsolconfs] = dualproof;
   conflictstore->dualprimalbnds[conflictstore->ndualsolconfs] = SCIPgetCutoffbound(set->scip) - SCIPsetSumepsilon(set);
   conflictstore->scalefactors[conflictstore->ndualsolconfs] = scale;
   conflictstore->updateside[conflictstore->ndualsolconfs] = updateside;
   ++conflictstore->ndualsolconfs;

   /* increase the number of non-zeros */
   SCIP_CALL( SCIPconsGetNVars(dualproof, set, &nvars, &success) );
   assert(success);
   conflictstore->nnzdualsols += nvars;

   return SCIP_OKAY;
}

/** adds a conflict to the conflict store
 *
 *  @note this method captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddConflict(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree (or NULL for an original constraint) */
   SCIP_PROB*            transprob,          /**< transformed problem (or NULL for an original constraint) */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_CONS*            cons,               /**< constraint representing the conflict */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             cutoffinvolved,     /**< is a cutoff bound involved in this conflict */
   SCIP_Real             primalbound         /**< primal bound the conflict depend on (or -SCIPinfinity) */
   )
{
   SCIP_Longint curnodenum;
   int nconflicts;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL || SCIPconsIsOriginal(cons));
   assert(transprob != NULL || SCIPconsIsOriginal(cons));
   assert(cons != NULL);
   assert(conftype != SCIP_CONFTYPE_BNDEXCEEDING || cutoffinvolved);
   assert(!cutoffinvolved || !SCIPsetIsInfinity(set, REALABS(primalbound)));

   /* mark the constraint to be a conflict */
   SCIPconsMarkConflict(cons);

   /* add the constraint to a special store */
   if( SCIPconsIsOriginal(cons) )
   {
      assert(SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM);
      SCIP_CALL( conflictstoreAddOrigConflict(conflictstore, set, blkmem, cons) );
      return SCIP_OKAY;
   }

   nconflicts = conflictstore->nconflicts;

   /* initialize the storage */
   if( conflictstore->maxstoresize == -1 )
   {
      SCIP_CALL( initConflictstore(conflictstore, set, transprob) );
   }
   assert(conflictstore->initstoresize >= 0);
   assert(conflictstore->initstoresize <= conflictstore->maxstoresize);

   /* return if conflict pool is disabled */
   if( conflictstore->maxstoresize <= 0 )
      return SCIP_OKAY;

   SCIP_CALL( conflictstoreEnsureMem(conflictstore, set, blkmem, nconflicts+1) );

   /* return if the store has size zero */
   if( conflictstore->conflictsize == 0 )
   {
      assert(conflictstore->maxstoresize == 0);
      return SCIP_OKAY;
   }

   assert(tree != NULL);
   curnodenum = (SCIPtreeGetFocusNode(tree) == NULL ? -1 : SCIPnodeGetNumber(SCIPtreeGetFocusNode(tree)));

   /* clean up the storage if we are at a new node or the storage is full */
   if( conflictstore->lastnodenum != curnodenum || conflictstore->nconflicts == conflictstore->conflictsize )
   {
      SCIP_CALL( conflictstoreCleanUpStorage(conflictstore, set, stat, transprob, blkmem, reopt) );
   }

   /* update the last seen node */
   conflictstore->lastnodenum = curnodenum;

   SCIPconsCapture(cons);
   conflictstore->conflicts[conflictstore->nconflicts] = cons;
   conflictstore->confprimalbnds[conflictstore->nconflicts] = primalbound;
   conflictstore->ncbconflicts += (SCIPsetIsInfinity(set, REALABS(primalbound)) ? 0 : 1);

   ++conflictstore->nconflicts;
   ++conflictstore->nconflictsfound;

#ifdef SCIP_PRINT_DETAILS
   SCIPsetDebugMsg(set, "add conflict <%s> to conflict store at position %d\n", SCIPconsGetName(cons), conflictstore->nconflicts-1);
   SCIPsetDebugMsg(set, " -> conflict type: %d, cutoff involved = %u\n", conftype, cutoffinvolved);
   if( cutoffinvolved )
      SCIPsetDebugMsg(set, " -> current primal bound: %g\n", primalbound);
#endif

   return SCIP_OKAY;
}

/** deletes all conflicts depending on a cutoff bound larger than the given bound */
SCIP_RETCODE SCIPconflictstoreCleanNewIncumbent(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem*/
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_Real             cutoffbound         /**< current cutoff bound */
   )
{
   SCIP_Real improvement;
   int ndelconfs;
   int nchgsides;
   int i;

   assert(conflictstore != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(transprob != NULL);

   /* return if we do not want to use the storage */
   if( set->conf_maxstoresize == 0 )
      return SCIP_OKAY;

   /* return if we do not want to remove conflicts related to an older cutoff bound */
   if( !set->conf_cleanbnddepend )
      return SCIP_OKAY;

   /* there is nothing to clean */
   if( conflictstore->ndualsolconfs == 0 && conflictstore->nconflicts == 0 )
      return SCIP_OKAY;

   /* we can stop whenever we have found a new incumbent but the cutoff bound has not changed */
   if( conflictstore->lastcutoffbound != SCIP_INVALID && SCIPsetIsGE(set, cutoffbound, conflictstore->lastcutoffbound) ) /*lint !e777*/
      return SCIP_OKAY;

   conflictstore->lastcutoffbound = cutoffbound;

   /* calculate scalar to determine whether the old primal bound is worse enough to remove the conflict */
   if( SCIPsetIsPositive(set, cutoffbound) )
      improvement = (1 - set->conf_minimprove);
   else
      improvement = (1 + set->conf_minimprove);

   /* remove all conflicts depending on a primalbound*improvement > cutoffbound
    *
    * note: we cannot remove conflicts that are marked as deleted because at this point in time we would destroy
    *       the internal data structure
    */
   ndelconfs = 0;
   for( i = 0; i < conflictstore->nconflicts; )
   {
      assert(conflictstore->conflicts[i] != NULL);

      /* check if the conflict depends on the cutoff bound */
      if( SCIPsetIsGT(set, improvement * conflictstore->confprimalbnds[i], cutoffbound) )
      {
         /* remove conflict at current position
          *
          * don't increase i because delPosConflict will swap the last pointer to the i-th position
          */
         SCIP_CALL( delPosConflict(conflictstore, set, stat, transprob, blkmem, reopt, i, TRUE) );
         ++ndelconfs;
      }
      else
         /* increase i */
         ++i;
   }
   assert(conflictstore->ncbconflicts >= 0);
   assert(conflictstore->nconflicts >= 0);

   SCIPsetDebugMsg(set, "-> removed %d/%d conflicts, %d depending on cutoff bound\n", ndelconfs,
         conflictstore->nconflicts+ndelconfs, ndelconfs);

   ndelconfs = 0;
   nchgsides = 0;
   /* update all proof constraints based on a dual solution */
   for( i = 0; i < conflictstore->ndualsolconfs; )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONS* dualproof;

      dualproof = conflictstore->dualsolconfs[i];
      assert(dualproof != NULL);

      if( SCIPconsIsDeleted(dualproof) )
      {
         ++i;
         continue;
      }
      if( !conflictstore->updateside[i] || SCIPsetIsLE(set, improvement * conflictstore->dualprimalbnds[i], cutoffbound) )
      {
         ++i;
         continue;
      }
      conshdlr = SCIPconsGetHdlr(dualproof);
      assert(conshdlr != NULL);

      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
      {
         SCIP_Real rhs;
         SCIP_Real newside;

         assert(SCIPsetIsGT(set, conflictstore->dualprimalbnds[i], cutoffbound));

         rhs = SCIPgetRhsLinear(NULL, dualproof);

         if( !SCIPsetIsInfinity(set, rhs) )
         {
            assert(SCIPsetIsInfinity(set, -SCIPgetLhsLinear(NULL, dualproof)));
            assert(SCIPsetIsPositive(set, conflictstore->scalefactors[i]));

            /* get unscaled rhs */
            newside = rhs * conflictstore->scalefactors[i];
            newside -= conflictstore->dualprimalbnds[i];
            newside += cutoffbound - SCIPsetSumepsilon(set);

            /* scale rhs */
            newside /= conflictstore->scalefactors[i];

            SCIP_CALL( SCIPchgRhsLinear(set->scip, dualproof, newside) );
         }
         else
         {
            SCIP_Real lhs = SCIPgetLhsLinear(NULL, dualproof);
            assert(!SCIPsetIsInfinity(set, -lhs));
            assert(SCIPsetIsNegative(set, conflictstore->scalefactors[i]));

            /* get unscaled lhs */
            newside = lhs * conflictstore->scalefactors[i];
            newside += conflictstore->dualprimalbnds[i];
            newside -= (cutoffbound - SCIPsetSumepsilon(set));

            /* scale lhs */
            newside /= conflictstore->scalefactors[i];

            SCIP_CALL( SCIPchgLhsLinear(set->scip, dualproof, newside) );
         }

         ++nchgsides;

         conflictstore->dualprimalbnds[i] = cutoffbound - SCIPsetSumepsilon(set);

         ++i;
      }
      else if( SCIPsetIsGT(set, improvement * conflictstore->dualprimalbnds[i], cutoffbound) )
      {
         /* remove conflict at current position
          *
          * don't increase i because delPosDualsol will swap the last pointer to the i-th position
          */
         SCIP_CALL( delPosDualsol(conflictstore, set, stat, transprob, blkmem, reopt, i, TRUE) );
         ++ndelconfs;
      }
      else
         /* increase i */
         ++i;
   }

   SCIPsetDebugMsg(set, "-> changed %d sides of dual solution constraints\n", nchgsides);
   SCIPsetDebugMsg(set, "-> deleted %d dual solution constraints\n", ndelconfs);

   return SCIP_OKAY;
}

/** returns the maximal size of the conflict pool */
int SCIPconflictstoreGetMaxPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   )
{
   assert(conflictstore != NULL);

   return MIN(conflictstore->storesize, conflictstore->maxstoresize);
}

/** returns the initial size of the conflict pool */
int SCIPconflictstoreGetInitPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   )
{
   assert(conflictstore != NULL);

   return conflictstore->initstoresize;
}

/** returns the number of stored conflicts on the conflict pool
 *
 *  @note the number of active conflicts can be less
 */
int SCIPconflictstoreGetNConflictsInStore(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   )
{
   assert(conflictstore != NULL);

   return conflictstore->nconflicts;
}

/** returns all active conflicts stored in the conflict store */
SCIP_RETCODE SCIPconflictstoreGetConflicts(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_CONS**           conflicts,          /**< array to store conflicts */
   int                   conflictsize,       /**< size of the conflict array */
   int*                  nconflicts          /**< pointer to store the number of conflicts */
   )
{
   int i;

   assert(conflictstore != NULL);

   /* return if the allocated memory is obviously to small */
   if( conflictstore->nconflicts > conflictsize )
   {
      (*nconflicts) = conflictstore->nconflicts;
      return SCIP_OKAY;
   }

   (*nconflicts) = 0;
   for( i = 0; i < conflictstore->nconflicts; i++ )
   {
      SCIP_CONS* conflict;

      conflict = conflictstore->conflicts[i];
      assert(conflict != NULL);

      /* skip deactivated and deleted constraints */
      if( !SCIPconsIsActive(conflict) || SCIPconsIsDeleted(conflict) )
         continue;

      /* count exact number conflicts */
      if( *nconflicts > conflictsize )
         ++(*nconflicts);
      else
      {
         conflicts[*nconflicts] = conflict;
         ++(*nconflicts);
      }
   }

   return SCIP_OKAY;
}

/** transformes all original conflicts into transformed conflicts */
SCIP_RETCODE SCIPconflictstoreTransform(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_REOPT*           reopt               /**< reoptimization data */
   )
{
   int ntransconss;
   int i;

   assert(conflictstore != NULL);
   assert(set != NULL);
   assert(SCIPsetGetStage(set) == SCIP_STAGE_TRANSFORMING);

   /* return if no original constraints are stored */
   if( conflictstore->norigconfs == 0 )
      return SCIP_OKAY;

   ntransconss = 0;

   for( i = 0; i < conflictstore->norigconfs; i++ )
   {
      SCIP_CONS* transcons;

      assert(conflictstore->origconfs[i] != NULL);
      assert(SCIPconsIsOriginal(conflictstore->origconfs[i]));

      transcons = SCIPconsGetTransformed(conflictstore->origconfs[i]);

      if( transcons != NULL )
      {
         SCIP_CALL( SCIPconflictstoreAddConflict(conflictstore, blkmem, set, stat, tree, transprob, reopt, transcons, \
               SCIP_CONFTYPE_UNKNOWN, FALSE, -SCIPsetInfinity(set)) );

         ++ntransconss;
      }

      SCIP_CALL( SCIPconsRelease(&conflictstore->origconfs[i], blkmem, set) );
   }

   SCIPsetDebugMsg(set, "-> transform %d/%d conflicts into transformed space\n", ntransconss, conflictstore->norigconfs);

   conflictstore->norigconfs = 0;

   return SCIP_OKAY;
}

/** returns the average number of non-zeros over all stored dual ray constraints */
SCIP_Real SCIPconflictstoreGetAvgNnzDualInfProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   )
{
   assert(conflictstore != NULL);

   if( conflictstore->ndualrayconfs == 0 )
      return 0.0;
   else
      return (SCIP_Real) conflictstore->nnzdualrays / ((SCIP_Real) conflictstore->ndualrayconfs);
}

/** returns the number of all stored dual ray constraints */
int SCIPconflictstoreGetNDualInfProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   )
{
   assert(conflictstore != NULL);

   return conflictstore->ndualrayconfs;
}

/** returns the average number of non-zeros over all stored boundexceeding proofs */
SCIP_Real SCIPconflictstoreGetAvgNnzDualBndProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   )
{
   assert(conflictstore != NULL);
   assert(conflictstore->ndualsolconfs >= 0);

   if( conflictstore->ndualsolconfs == 0 )
      return 0.0;
   else
      return (SCIP_Real) conflictstore->nnzdualsols / ((SCIP_Real) conflictstore->ndualsolconfs);
}

/** returns the number of all stored boundexceeding proofs */
int SCIPconflictstoreGetNDualBndProofs(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict store */
   )
{
   assert(conflictstore != NULL);

   return conflictstore->ndualsolconfs;
}

