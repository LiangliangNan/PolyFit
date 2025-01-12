/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   benderscut.c
 * @ingroup OTHER_CFILES
 * @brief  methods for Benders' decomposition cut
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/benderscut.h"
#include "scip/reopt.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_benders.h"

#include "scip/struct_benderscut.h"

/* default parameter settings for the Benders' decomposition cuts */
#define SCIP_DEFAULT_ENABLED        TRUE

/** compares two Benders' cuts w. r. to their delay positions and their priority */
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutComp)
{  /*lint --e{715}*/
   SCIP_BENDERSCUT* benderscut1 = (SCIP_BENDERSCUT*)elem1;
   SCIP_BENDERSCUT* benderscut2 = (SCIP_BENDERSCUT*)elem2;

   assert(benderscut1 != NULL);
   assert(benderscut2 != NULL);

   return benderscut2->priority - benderscut1->priority; /* prefer higher priorities */
}

/** comparison method for sorting Benders' cuts w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPbenderscutCompName)
{
   return strcmp(SCIPbenderscutGetName((SCIP_BENDERSCUT*)elem1), SCIPbenderscutGetName((SCIP_BENDERSCUT*)elem2));
}

/** method to call, when the priority of a compression was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdBenderscutPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetBenderscutPriority() to mark the compressions unsorted */
   SCIP_CALL( SCIPsetBenderscutPriority(scip, (SCIP_BENDERSCUT*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given Benders' decomposition cut to a new scip */
SCIP_RETCODE SCIPbenderscutCopyInclude(
   SCIP_BENDERS*         benders,            /**< the Benders' decomposition that the cuts are copied to */
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(benderscut != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( benderscut->benderscutcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including benderscut %s in subscip %p\n", SCIPbenderscutGetName(benderscut), (void*)set->scip);
      SCIP_CALL( benderscut->benderscutcopy(set->scip, benders, benderscut) );
   }

   return SCIP_OKAY;
}

/** internal method for creating a Benders' decomposition structure */
static
SCIP_RETCODE doBenderscutCreate(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT**     benderscut,         /**< pointer to the Benders' decomposition cut data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of the Benders' decomposition cut */
   const char*           desc,               /**< description of the Benders' decomposition cut */
   int                   priority,           /**< priority of the the Benders' decomposition cut */
   SCIP_Bool             islpcut,            /**< indicates whether the cut is generated from the LP solution */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy)),/**< copy method of the Benders' decomposition cut or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree)),/**< destructor of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit)),/**< initialize the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit)),/**< deinitialize the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol)),/**< solving process initialization method of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol)),/**< solving process deinitialization method of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< execution method of the Benders' decomposition cut */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' decomposition cut data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(benderscut != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(benderscutexec != NULL);

   SCIP_ALLOC( BMSallocMemory(benderscut) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benderscut)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*benderscut)->desc, desc, strlen(desc)+1) );
   (*benderscut)->priority = priority;
   (*benderscut)->islpcut = islpcut;
   (*benderscut)->benderscutcopy = benderscutcopy;
   (*benderscut)->benderscutfree = benderscutfree;
   (*benderscut)->benderscutinit = benderscutinit;
   (*benderscut)->benderscutexit = benderscutexit;
   (*benderscut)->benderscutinitsol = benderscutinitsol;
   (*benderscut)->benderscutexitsol = benderscutexitsol;
   (*benderscut)->benderscutexec = benderscutexec;
   (*benderscut)->benderscutdata = benderscutdata;
   SCIP_CALL( SCIPclockCreate(&(*benderscut)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*benderscut)->benderscutclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*benderscut)->ncalls = 0;
   (*benderscut)->nfound = 0;
   (*benderscut)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/benderscut/%s/priority", SCIPbendersGetName(benders), name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of Benders' cut <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*benderscut)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdBenderscutPriority, (SCIP_PARAMDATA*)(*benderscut)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "benders/%s/benderscut/%s/enabled", SCIPbendersGetName(benders), name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
        "is this Benders' decomposition cut method used to generate cuts?", &(*benderscut)->enabled, FALSE,
        SCIP_DEFAULT_ENABLED, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutCreate(
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_BENDERSCUT**     benderscut,         /**< pointer to the Benders' decomposition cut data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of the Benders' decomposition cut */
   const char*           desc,               /**< description of the Benders' decomposition cut */
   int                   priority,           /**< priority of the the Benders' decomposition cut */
   SCIP_Bool             islpcut,            /**< indicates whether the cut is generated from the LP solution */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy)),/**< copy method of the Benders' decomposition cut or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree)),/**< destructor of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit)),/**< initialize the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit)),/**< deinitialize the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol)),/**< solving process initialization method of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol)),/**< solving process deinitialization method of the Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXEC((*benderscutexec)),/**< execution method of the Benders' decomposition cut */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< Benders' decomposition cut data */
   )
{
   assert(benderscut != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(benderscutexec != NULL);

   SCIP_CALL_FINALLY( doBenderscutCreate(benders, benderscut, set, messagehdlr, blkmem, name, desc, priority, islpcut,
         benderscutcopy, benderscutfree, benderscutinit, benderscutexit, benderscutinitsol, benderscutexitsol,
         benderscutexec, benderscutdata), (void)SCIPbenderscutFree(benderscut, set) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutFree(
   SCIP_BENDERSCUT**     benderscut,         /**< pointer to the Benders' decomposition cut data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benderscut != NULL);
   assert(*benderscut != NULL);
   assert(!(*benderscut)->initialized);
   assert(set != NULL);

   /* call destructor of the Benders' decomposition cut */
   if( (*benderscut)->benderscutfree != NULL )
   {
      SCIP_CALL( (*benderscut)->benderscutfree(set->scip, *benderscut) );
   }

   SCIPclockFree(&(*benderscut)->benderscutclock);
   SCIPclockFree(&(*benderscut)->setuptime);
   BMSfreeMemoryArray(&(*benderscut)->name);
   BMSfreeMemoryArray(&(*benderscut)->desc);
   BMSfreeMemory(benderscut);

   return SCIP_OKAY;
}

/** initializes the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutInit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benderscut != NULL);
   assert(set != NULL);

   if( benderscut->initialized )
   {
      SCIPerrorMessage("Benders' decomposition cut <%s> already initialized\n", benderscut->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(benderscut->setuptime);
      SCIPclockReset(benderscut->benderscutclock);

      benderscut->ncalls = 0;
      benderscut->nfound = 0;
   }

   if( benderscut->benderscutinit != NULL )
   {
      /* start timing */
      SCIPclockStart(benderscut->setuptime, set);

      SCIP_CALL( benderscut->benderscutinit(set->scip, benderscut) );

      /* stop timing */
      SCIPclockStop(benderscut->setuptime, set);
   }
   benderscut->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutExit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benderscut != NULL);
   assert(set != NULL);

   if( !benderscut->initialized )
   {
      SCIPerrorMessage("Benders' decomposition cut <%s> not initialized\n", benderscut->name);
      return SCIP_INVALIDCALL;
   }

   if( benderscut->benderscutexit != NULL )
   {
      /* start timing */
      SCIPclockStart(benderscut->setuptime, set);

      SCIP_CALL( benderscut->benderscutexit(set->scip, benderscut) );

      /* stop timing */
      SCIPclockStop(benderscut->setuptime, set);
   }
   benderscut->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs Benders' cut that the branch and bound process is being started */
SCIP_RETCODE SCIPbenderscutInitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benderscut != NULL);
   assert(set != NULL);

   /* call solving process initialization method of the Benders' decomposition cut */
   if( benderscut->benderscutinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benderscut->setuptime, set);

      SCIP_CALL( benderscut->benderscutinitsol(set->scip, benderscut) );

      /* stop timing */
      SCIPclockStop(benderscut->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs Benders' decomposition that the branch and bound process data is being freed */
SCIP_RETCODE SCIPbenderscutExitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(benderscut != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of Benders' decomposition cut */
   if( benderscut->benderscutexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(benderscut->setuptime, set);

      SCIP_CALL( benderscut->benderscutexitsol(set->scip, benderscut) );

      /* stop timing */
      SCIPclockStop(benderscut->setuptime, set);
   }

   return SCIP_OKAY;
}

/** calls execution method of the Benders' decomposition cut */
SCIP_RETCODE SCIPbenderscutExec(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BENDERS*         benders,            /**< Benders' decomposition */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   int                   probnumber,         /**< the number of the subproblem for which the cut is generated */
   SCIP_BENDERSENFOTYPE  type,               /**< the enforcement type calling this function */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_RESULT cutresult;

   assert(benderscut != NULL);
   assert(benderscut->benderscutexec != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   cutresult = SCIP_DIDNOTRUN;

   SCIPsetDebugMsg(set, "executing Benders' decomposition cut <%s>\n", benderscut->name);

   /* start timing */
   SCIPclockStart(benderscut->benderscutclock, set);

   /* call the Benders' decomposition cut if it is enabled */
   if( benderscut->enabled )
   {
      SCIP_CALL( benderscut->benderscutexec(set->scip, benders, benderscut, sol, probnumber, type, &cutresult) );
   }

   /* stop timing */
   SCIPclockStop(benderscut->benderscutclock, set);

   /* evaluate result */
   if( cutresult != SCIP_DIDNOTRUN
      && cutresult != SCIP_DIDNOTFIND
      && cutresult != SCIP_CONSADDED
      && cutresult != SCIP_FEASIBLE
      && cutresult != SCIP_SEPARATED )
   {
      SCIPerrorMessage("execution method of Benders' decomposition cut <%s> returned invalid result <%d>\n",
         benderscut->name, cutresult);
      return SCIP_INVALIDRESULT;
   }

   benderscut->ncalls++;

   if( cutresult == SCIP_CONSADDED || cutresult == SCIP_SEPARATED )
      benderscut->nfound++;

   (*result) = cutresult;

   return SCIP_OKAY;
}

/** gets user data of the Benders' decomposition cut */
SCIP_BENDERSCUTDATA* SCIPbenderscutGetData(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->benderscutdata;
}

/** sets user data of the Benders' decomposition cut; user has to free old data in advance! */
void SCIPbenderscutSetData(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_BENDERSCUTDATA*  benderscutdata      /**< new Benders' decomposition cut user data */
   )
{
   assert(benderscut != NULL);

   benderscut->benderscutdata = benderscutdata;
}

/* new callback setter methods */

/** sets copy callback of the Benders' decomposition cut */
void SCIPbenderscutSetCopy(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTCOPY((*benderscutcopy))/**< copy callback of the Benders' decomposition cut or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(benderscut != NULL);

   benderscut->benderscutcopy = benderscutcopy;
}

/** sets destructor callback of the Benders' decomposition cut */
void SCIPbenderscutSetFree(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTFREE((*benderscutfree))/**< destructor of the Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   benderscut->benderscutfree = benderscutfree;
}

/** sets initialization callback of the Benders' decomposition cut */
void SCIPbenderscutSetInit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINIT((*benderscutinit))/**< initialize the Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   benderscut->benderscutinit = benderscutinit;
}

/** sets deinitialization callback of the Benders' decomposition cut */
void SCIPbenderscutSetExit(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXIT((*benderscutexit))/**< deinitialize the Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   benderscut->benderscutexit = benderscutexit;
}

/** sets solving process initialization callback of the Benders' decomposition cut */
void SCIPbenderscutSetInitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTINITSOL((*benderscutinitsol))/**< solving process initialization callback of the Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   benderscut->benderscutinitsol = benderscutinitsol;
}

/** sets solving process deinitialization callback of Benders' decomposition cut */
void SCIPbenderscutSetExitsol(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_DECL_BENDERSCUTEXITSOL((*benderscutexitsol))/**< solving process deinitialization callback of the Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   benderscut->benderscutexitsol = benderscutexitsol;
}

/** gets name of the Benders' decomposition cut */
const char* SCIPbenderscutGetName(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->name;
}

/** gets description of the Benders' decomposition cut */
const char* SCIPbenderscutGetDesc(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->desc;
}

/** gets priority of the Benders' decomposition cut */
int SCIPbenderscutGetPriority(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->priority;
}

/** sets priority of the Benders' decomposition cut */
void SCIPbenderscutSetPriority(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   int                   priority            /**< new priority of the Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   benderscut->priority = priority;
}

/** gets the number of times, the heuristic was called and tried to find a solution */
SCIP_Longint SCIPbenderscutGetNCalls(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->ncalls;
}

/** gets the number of Benders' cuts found by this Benders' decomposition cut */
SCIP_Longint SCIPbenderscutGetNFound(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->nfound;
}

/** is the Benders' decomposition cut initialized? */
SCIP_Bool SCIPbenderscutIsInitialized(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->initialized;
}

/** gets time in seconds used by this Benders' decomposition cut for setting up */
SCIP_Real SCIPbenderscutGetSetupTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return SCIPclockGetTime(benderscut->setuptime);
}

/** gets time in seconds used in this Benders' decomposition cut */
SCIP_Real SCIPbenderscutGetTime(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return SCIPclockGetTime(benderscut->benderscutclock);
}

/** returns whether the Benders' cut uses the LP information */
SCIP_Bool SCIPbenderscutIsLPCut(
   SCIP_BENDERSCUT*      benderscut          /**< Benders' decomposition cut */
   )
{
   assert(benderscut != NULL);

   return benderscut->islpcut;
}

/** sets the enabled flag of the Benders' decomposition cut method */
void SCIPbenderscutSetEnabled(
   SCIP_BENDERSCUT*      benderscut,         /**< Benders' decomposition cut */
   SCIP_Bool             enabled             /**< flag to indicate whether the Benders' decomposition cut is enabled */
   )
{
   assert(benderscut != NULL);

   benderscut->enabled = enabled;
}
