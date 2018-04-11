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

/**@file   compr.c
 * @brief  methods for tree compressions
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/compr.h"
#include "scip/reopt.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_compr.h"



/** compares two compression methods w. r. to their delay positions and their priority */
SCIP_DECL_SORTPTRCOMP(SCIPcomprComp)
{  /*lint --e{715}*/
   SCIP_COMPR* compr1 = (SCIP_COMPR*)elem1;
   SCIP_COMPR* compr2 = (SCIP_COMPR*)elem2;

   assert(compr1 != NULL);
   assert(compr2 != NULL);

   return compr2->priority - compr1->priority; /* prefer higher priorities */
}

/** comparison method for sorting heuristics w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPcomprCompName)
{
   return strcmp(SCIPcomprGetName((SCIP_COMPR*)elem1), SCIPcomprGetName((SCIP_COMPR*)elem2));
}

/** method to call, when the priority of a compression was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdComprPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetComprPriority() to mark the compressions unsorted */
   SCIP_CALL( SCIPsetComprPriority(scip, (SCIP_COMPR*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given tree compression to a new scip */
SCIP_RETCODE SCIPcomprCopyInclude(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(compr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( compr->comprcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including compr %s in subscip %p\n", SCIPcomprGetName(compr), (void*)set->scip);
      SCIP_CALL( compr->comprcopy(set->scip, compr) );
   }

   return SCIP_OKAY;
}

/** creates a tree compression */
SCIP_RETCODE SCIPcomprCreate(
   SCIP_COMPR**          compr,              /**< pointer to tree compression data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of tree compression */
   const char*           desc,               /**< description of tree compression */
   int                   priority,           /**< priority of the tree compression */
   int                   minnnodes,          /**< minimal number of nodes for calling compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy)),     /**< copy method of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_COMPRFREE   ((*comprfree)),     /**< destructor of tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit)),     /**< initialize tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit)),     /**< deinitialize tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol)), /**< solving process initialization method of tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol)), /**< solving process deinitialization method of tree compression */
   SCIP_DECL_COMPREXEC   ((*comprexec)),     /**< execution method of tree compression */
   SCIP_COMPRDATA*       comprdata           /**< tree compression data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(compr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(comprexec != NULL);

   SCIP_ALLOC( BMSallocMemory(compr) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*compr)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*compr)->desc, desc, strlen(desc)+1) );
   (*compr)->priority = priority;
   (*compr)->minnnodes = minnnodes;
   (*compr)->comprcopy = comprcopy;
   (*compr)->comprfree = comprfree;
   (*compr)->comprinit = comprinit;
   (*compr)->comprexit = comprexit;
   (*compr)->comprinitsol = comprinitsol;
   (*compr)->comprexitsol = comprexitsol;
   (*compr)->comprexec = comprexec;
   (*compr)->comprdata = comprdata;
   SCIP_CALL( SCIPclockCreate(&(*compr)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*compr)->comprclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*compr)->ncalls = 0;
   (*compr)->nfound = 0;
   (*compr)->rate = 0.0;
   (*compr)->initialized = FALSE;
   (*compr)->nnodes = 0;
   (*compr)->loi = 0.0;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "compression/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of compression <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*compr)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdComprPriority, (SCIP_PARAMDATA*)(*compr)) ); /*lint !e740*/
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "compression/%s/minnleaves", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "minimal number of leave nodes for calling tree compression <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*compr)->minnnodes, FALSE, minnnodes, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of tree compression */
SCIP_RETCODE SCIPcomprFree(
   SCIP_COMPR**          compr,              /**< pointer to tree compression data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(compr != NULL);
   assert(*compr != NULL);
   assert(!(*compr)->initialized);
   assert(set != NULL);

   /* call destructor of tree compression */
   if( (*compr)->comprfree != NULL )
   {
      SCIP_CALL( (*compr)->comprfree(set->scip, *compr) );
   }

   SCIPclockFree(&(*compr)->comprclock);
   SCIPclockFree(&(*compr)->setuptime);
   BMSfreeMemoryArray(&(*compr)->name);
   BMSfreeMemoryArray(&(*compr)->desc);
   BMSfreeMemory(compr);

   return SCIP_OKAY;
}

/** initializes tree compression */
SCIP_RETCODE SCIPcomprInit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(compr != NULL);
   assert(set != NULL);

   if( compr->initialized )
   {
      SCIPerrorMessage("tree compression <%s> already initialized\n", compr->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat && !set->reopt_enable )
   {
      SCIPclockReset(compr->setuptime);
      SCIPclockReset(compr->comprclock);

      compr->ncalls = 0;
      compr->nfound = 0;
   }

   if( compr->comprinit != NULL )
   {
      /* start timing */
      SCIPclockStart(compr->setuptime, set);

      SCIP_CALL( compr->comprinit(set->scip, compr) );

      /* stop timing */
      SCIPclockStop(compr->setuptime, set);
   }
   compr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of tree compression */
SCIP_RETCODE SCIPcomprExit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(compr != NULL);
   assert(set != NULL);

   if( !compr->initialized )
   {
      SCIPerrorMessage("tree compression <%s> not initialized\n", compr->name);
      return SCIP_INVALIDCALL;
   }

   if( compr->comprexit != NULL )
   {
      /* start timing */
      SCIPclockStart(compr->setuptime, set);

      SCIP_CALL( compr->comprexit(set->scip, compr) );

      /* stop timing */
      SCIPclockStop(compr->setuptime, set);
   }
   compr->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls execution method of tree compression */
SCIP_RETCODE SCIPcomprExec(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(compr != NULL);
   assert(compr->comprexec != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* do not run if reoptimization data structure is not initialized */
   if( reopt == NULL )
      return SCIP_OKAY;

   /* do not run if the reoptimization tree is not large enough */
   if( SCIPreoptGetNLeaves(reopt, NULL) < compr->minnnodes )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "executing tree compression <%s>\n", compr->name);

   /* start timing */
   SCIPclockStart(compr->comprclock, set);

   /* call external method */
   SCIP_CALL( compr->comprexec(set->scip, compr, result) );

   /* stop timing */
   SCIPclockStop(compr->comprclock, set);

   /* evaluate result */
   if( *result != SCIP_SUCCESS
      && *result != SCIP_DIDNOTFIND
      && *result != SCIP_DIDNOTRUN )
   {
      SCIPerrorMessage("execution method of tree compression <%s> returned invalid result <%d>\n",
         compr->name, *result);
      return SCIP_INVALIDRESULT;
   }

   if( *result != SCIP_DIDNOTRUN )
      compr->ncalls++;

   if( *result == SCIP_SUCCESS )
      compr->nfound++;

   return SCIP_OKAY;
}

/** gets user data of tree compression */
SCIP_COMPRDATA* SCIPcomprGetData(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->comprdata;
}

/** sets user data of tree compression; user has to free old data in advance! */
void SCIPcomprSetData(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_COMPRDATA*       comprdata           /**< new tree compression user data */
   )
{
   assert(compr != NULL);

   compr->comprdata = comprdata;
}

/* new callback setter methods */

/** sets copy callback of tree compression */
void SCIPcomprSetCopy(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy))      /**< copy callback of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(compr != NULL);

   compr->comprcopy = comprcopy;
}

/** sets destructor callback of tree compression */
void SCIPcomprSetFree(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRFREE   ((*comprfree))      /**< destructor of tree compression */
   )
{
   assert(compr != NULL);

   compr->comprfree = comprfree;
}

/** sets initialization callback of tree compression */
void SCIPcomprSetInit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit))      /**< initialize tree compression */
   )
{
   assert(compr != NULL);

   compr->comprinit = comprinit;
}

/** sets deinitialization callback of tree compression */
void SCIPcomprSetExit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit))      /**< deinitialize tree compression */
   )
{
   assert(compr != NULL);

   compr->comprexit = comprexit;
}

/** sets solving process initialization callback of tree compression */
void SCIPcomprSetInitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol))  /**< solving process initialization callback of tree compression */
   )
{
   assert(compr != NULL);

   compr->comprinitsol = comprinitsol;
}

/** sets solving process deinitialization callback of tree compression */
void SCIPcomprSetExitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol))  /**< solving process deinitialization callback of tree compression */
   )
{
   assert(compr != NULL);

   compr->comprexitsol = comprexitsol;
}

/** should the compression be executed at the given depth, number of nodes */
SCIP_Bool SCIPcomprShouldBeExecuted(
   SCIP_COMPR*           compr,              /**< tree compression */
   int                   depth,              /**< depth of current node */
   int                   nnodes              /**< number of open nodes */
   )
{
   assert(compr != NULL);
   assert(depth >= 0);
   assert(nnodes >= 0);

   return nnodes >= compr->minnnodes;
}

/** gets name of tree compression */
const char* SCIPcomprGetName(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->name;
}

/** gets description of tree compression */
const char* SCIPcomprGetDesc(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->desc;
}

/** gets priority of tree compression */
int SCIPcomprGetPriority(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->priority;
}

/** sets priority of tree compression */
void SCIPcomprSetPriority(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the tree compression */
   )
{
   assert(compr != NULL);
   assert(set != NULL);

   compr->priority = priority;
   set->comprssorted = FALSE;
}

/** gets minimal number of nodes for calling tree compression (returns -1, if no node threshold exists) */
int SCIPcomprGetMinNodes(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->minnnodes;
}

/** gets the number of times, the heuristic was called and tried to find a solution */
SCIP_Longint SCIPcomprGetNCalls(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->ncalls;
}

/** gets the number of compressions found by this compression */
SCIP_Longint SCIPcomprGetNFound(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->nfound;
}

/** is tree compression initialized? */
SCIP_Bool SCIPcomprIsInitialized(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return compr->initialized;
}

/** gets time in seconds used in this heuristic for setting up for next stages */
SCIP_Real SCIPcomprGetSetupTime(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return SCIPclockGetTime(compr->setuptime);
}

/** gets time in seconds used in this heuristic */
SCIP_Real SCIPcomprGetTime(
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(compr != NULL);

   return SCIPclockGetTime(compr->comprclock);
}
