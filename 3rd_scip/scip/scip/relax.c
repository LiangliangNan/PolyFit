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

/**@file   relax.c
 * @brief  methods and datastructures for relaxation handlers
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/sol.h"
#include "scip/var.h"
#include "scip/relax.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_relax.h"



/** compares two relaxation handlers w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPrelaxComp)
{  /*lint --e{715}*/
   return ((SCIP_RELAX*)elem2)->priority - ((SCIP_RELAX*)elem1)->priority;
}

/** comparison method for sorting relaxators w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPrelaxCompName)
{
   return strcmp(SCIPrelaxGetName((SCIP_RELAX*)elem1), SCIPrelaxGetName((SCIP_RELAX*)elem2));
}

/** method to call, when the priority of a relaxation handler was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdRelaxPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetRelaxPriority() to mark the relaxs unsorted */
   SCIP_CALL( SCIPsetRelaxPriority(scip, (SCIP_RELAX*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given relaxation handler to a new scip */
SCIP_RETCODE SCIPrelaxCopyInclude(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(relax != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( relax->relaxcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including relaxation handler %s in subscip %p\n", SCIPrelaxGetName(relax), (void*)set->scip);
      SCIP_CALL( relax->relaxcopy(set->scip, relax) );
   }
   return SCIP_OKAY;
}

/** creates a relaxation handler */
SCIP_RETCODE SCIPrelaxCreate(
   SCIP_RELAX**          relax,              /**< pointer to relaxation handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of relaxation handler */
   const char*           desc,               /**< description of relaxation handler */
   int                   priority,           /**< priority of the relaxation handler (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxation handler */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy)),     /**< copy method of relaxation handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxation handler */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxation handler */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxation handler */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxation handler */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxation handler */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< relaxation handler data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(relax != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(relaxexec != NULL);

   SCIP_ALLOC( BMSallocMemory(relax) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*relax)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*relax)->desc, desc, strlen(desc)+1) );
   (*relax)->priority = priority;
   (*relax)->freq = freq;
   (*relax)->relaxcopy = relaxcopy;
   (*relax)->relaxfree = relaxfree;
   (*relax)->relaxinit = relaxinit;
   (*relax)->relaxexit = relaxexit;
   (*relax)->relaxinitsol = relaxinitsol;
   (*relax)->relaxexitsol = relaxexitsol;
   (*relax)->relaxexec = relaxexec;
   (*relax)->relaxdata = relaxdata;
   SCIP_CALL( SCIPclockCreate(&(*relax)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*relax)->relaxclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*relax)->ncalls = 0;
   (*relax)->lastsolvednode = -1;
   (*relax)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "relaxing/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of relaxation handler <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*relax)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
         paramChgdRelaxPriority, (SCIP_PARAMDATA*)(*relax)) ); /*lint !e740*/
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "relaxing/%s/freq", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "frequency for calling relaxation handler <%s> (-1: never, 0: only in root node)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*relax)->freq, FALSE, freq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of relaxation handler */
SCIP_RETCODE SCIPrelaxFree(
   SCIP_RELAX**          relax,              /**< pointer to relaxation handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(*relax != NULL);
   assert(!(*relax)->initialized);
   assert(set != NULL);

   /* call destructor of relaxation handler */
   if( (*relax)->relaxfree != NULL )
   {
      SCIP_CALL( (*relax)->relaxfree(set->scip, *relax) );
   }

   SCIPclockFree(&(*relax)->relaxclock);
   SCIPclockFree(&(*relax)->setuptime);
   BMSfreeMemoryArray(&(*relax)->name);
   BMSfreeMemoryArray(&(*relax)->desc);
   BMSfreeMemory(relax);

   return SCIP_OKAY;
}

/** initializes relaxation handler */
SCIP_RETCODE SCIPrelaxInit(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   if( relax->initialized )
   {
      SCIPerrorMessage("relaxation handler <%s> already initialized\n", relax->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(relax->setuptime);
      SCIPclockReset(relax->relaxclock);
      relax->ncalls = 0;
      relax->lastsolvednode = -1;
   }

   if( relax->relaxinit != NULL )
   {
      /* start timing */
      SCIPclockStart(relax->setuptime, set);

      SCIP_CALL( relax->relaxinit(set->scip, relax) );

      /* stop timing */
      SCIPclockStop(relax->setuptime, set);
   }
   relax->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of relaxation handler */
SCIP_RETCODE SCIPrelaxExit(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   if( !relax->initialized )
   {
      SCIPerrorMessage("relaxation handler <%s> not initialized\n", relax->name);
      return SCIP_INVALIDCALL;
   }

   if( relax->relaxexit != NULL )
   {
      /* start timing */
      SCIPclockStart(relax->setuptime, set);

      SCIP_CALL( relax->relaxexit(set->scip, relax) );

      /* stop timing */
      SCIPclockStop(relax->setuptime, set);
   }
   relax->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs relaxation handler that the branch and bound process is being started */
SCIP_RETCODE SCIPrelaxInitsol(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   /* call solving process initialization method of relaxation handler */
   if( relax->relaxinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(relax->setuptime, set);

      SCIP_CALL( relax->relaxinitsol(set->scip, relax) );

      /* stop timing */
      SCIPclockStop(relax->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs relaxation handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPrelaxExitsol(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of relaxation handler */
   if( relax->relaxexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(relax->setuptime, set);

      SCIP_CALL( relax->relaxexitsol(set->scip, relax) );

      /* stop timing */
      SCIPclockStop(relax->setuptime, set);
   }

   return SCIP_OKAY;
}

/** calls execution method of relaxation handler */
SCIP_RETCODE SCIPrelaxExec(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Real*            lowerbound,         /**< pointer to lower bound computed by the relaxation handler */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(relax != NULL);
   assert(relax->relaxexec != NULL);
   assert(relax->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check, if the relaxation is already solved */
   if( relax->lastsolvednode == stat->ntotalnodes && ! SCIPinProbing(set->scip) )
      return SCIP_OKAY;

   relax->lastsolvednode = stat->ntotalnodes;

   if( (depth == 0 && relax->freq == 0) || (relax->freq > 0 && depth % relax->freq == 0) )
   {
      SCIPsetDebugMsg(set, "executing relaxation handler <%s>\n", relax->name);

      /* start timing */
      SCIPclockStart(relax->relaxclock, set);

      /* call external relaxation method */
      SCIP_CALL( relax->relaxexec(set->scip, relax, lowerbound, result) );

      /* stop timing */
      SCIPclockStop(relax->relaxclock, set);

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_CONSADDED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_SUCCESS
         && *result != SCIP_SUSPENDED
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("execution method of relaxation handler <%s> returned invalid result <%d>\n",
            relax->name, *result);
         return SCIP_INVALIDRESULT;
      }
      if( *result != SCIP_DIDNOTRUN )
      {
         relax->ncalls++;
         stat->relaxcount++;
         if( *result == SCIP_SUSPENDED )
            SCIPrelaxMarkUnsolved(relax);
      }
   }

   return SCIP_OKAY;
}

/** gets user data of relaxation handler */
SCIP_RELAXDATA* SCIPrelaxGetData(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return relax->relaxdata;
}

/** sets user data of relaxation handler; user has to free old data in advance! */
void SCIPrelaxSetData(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< new relaxation handler user data */
   )
{
   assert(relax != NULL);

   relax->relaxdata = relaxdata;
}

/** set copy method of relaxation handler */
void SCIPrelaxSetCopy(
   SCIP_RELAX*           relax,              /**< relaxation handler  */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy))      /**< copy method of relaxation handler */
   )
{
   assert(relax != NULL);

   relax->relaxcopy = relaxcopy;
}

/** set destructor of relaxation handler */
void SCIPrelaxSetFree(
   SCIP_RELAX*           relax,              /**< relaxation handler  */
   SCIP_DECL_RELAXFREE   ((*relaxfree))      /**< destructor of relaxation handler */
   )
{
   assert(relax != NULL);

   relax->relaxfree = relaxfree;
}

/** set initialization method of relaxation handler */
void SCIPrelaxSetInit(
   SCIP_RELAX*           relax,              /**< relaxation handler  */
   SCIP_DECL_RELAXINIT   ((*relaxinit))      /**< initialize relaxation handler */
   )
{
   assert(relax != NULL);

   relax->relaxinit = relaxinit;
}

/** set deinitialization method of relaxation handler */
void SCIPrelaxSetExit(
   SCIP_RELAX*           relax,              /**< relaxation handler  */
   SCIP_DECL_RELAXEXIT   ((*relaxexit))      /**< deinitialize relaxation handler */
   )
{
   assert(relax != NULL);

   relax->relaxexit = relaxexit;
}

/** set solving process initialization method of relaxation handler */
void SCIPrelaxSetInitsol(
   SCIP_RELAX*           relax,              /**< relaxation handler  */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol))   /**< solving process initialization method of relaxation handler */
   )
{
   assert(relax != NULL);

   relax->relaxinitsol = relaxinitsol;
}

/** set solving process deinitialization method of relaxation handler */
void SCIPrelaxSetExitsol(
   SCIP_RELAX*           relax,              /**< relaxation handler  */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol))   /**< solving process deinitialization relaxation handler */
   )
{
   assert(relax != NULL);

   relax->relaxexitsol = relaxexitsol;
}

/** gets name of relaxation handler */
const char* SCIPrelaxGetName(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return relax->name;
}

/** gets description of relaxation handler */
const char* SCIPrelaxGetDesc(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return relax->desc;
}

/** gets priority of relaxation handler */
int SCIPrelaxGetPriority(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return relax->priority;
}

/** sets priority of relaxation handler */
void SCIPrelaxSetPriority(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the relaxation handler */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   relax->priority = priority;
   set->relaxssorted = FALSE;
}

/** gets frequency of relaxation handler */
int SCIPrelaxGetFreq(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return relax->freq;
}

/** gets time in seconds used in this relaxator for setting up for next stages */
SCIP_Real SCIPrelaxGetSetupTime(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return SCIPclockGetTime(relax->setuptime);
}

/** enables or disables all clocks of \p relax, depending on the value of the flag */
void SCIPrelaxEnableOrDisableClocks(
   SCIP_RELAX*           relax,              /**< the relaxation handler for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the relaxation handler be enabled? */
   )
{
   assert(relax != NULL);

   SCIPclockEnableOrDisable(relax->setuptime, enable);
   SCIPclockEnableOrDisable(relax->relaxclock, enable);
}

/** gets time in seconds used in this relaxation handler */
SCIP_Real SCIPrelaxGetTime(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return SCIPclockGetTime(relax->relaxclock);
}

/** gets the total number of times, the relaxation handler was called */
SCIP_Longint SCIPrelaxGetNCalls(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return relax->ncalls;
}

/** is relaxation handler initialized? */
SCIP_Bool SCIPrelaxIsInitialized(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   return relax->initialized;
}

/** returns whether the relaxation was completely solved at the current node */
SCIP_Bool SCIPrelaxIsSolved(
   SCIP_RELAX*           relax,              /**< relaxation handler */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(relax != NULL);
   assert(stat != NULL);

   return (relax->lastsolvednode == stat->ntotalnodes);
}

/** marks the current relaxation unsolved, s.t. the relaxation handler is called again in the next solving round */
void SCIPrelaxMarkUnsolved(
   SCIP_RELAX*           relax               /**< relaxation handler */
   )
{
   assert(relax != NULL);

   relax->lastsolvednode = -1;
}

/* 
 *  methods for the global relaxation data 
 */

/** creates global relaxation data */
SCIP_RETCODE SCIPrelaxationCreate(
   SCIP_RELAXATION**     relaxation,         /**< global relaxation data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(relaxation != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(primal != NULL);
   assert(tree != NULL);

   SCIP_ALLOC( BMSallocMemory(relaxation) );

   (*relaxation)->relaxsolobjval = 0.0;
   (*relaxation)->relaxsolvalid = FALSE;
   (*relaxation)->relaxsolincludeslp = FALSE;
   (*relaxation)->relaxsolzero = TRUE;

   return SCIP_OKAY;
}

/** frees global relaxation data */
SCIP_RETCODE SCIPrelaxationFree(
   SCIP_RELAXATION**     relaxation          /**< global relaxation data */
   )
{
   assert(relaxation != NULL);

   BMSfreeMemory(relaxation);

   return SCIP_OKAY;
}

/** sets the relaxsolzero flag in the relaxation data to the given value */
void SCIPrelaxationSetSolZero(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Bool             iszero              /**< are all values of the relaxation solution set to zero? */
   )
{
   assert(relaxation != NULL);

   relaxation->relaxsolzero = iszero;
}

/** returns whether the global relaxation solution is cleared and all values are set to zero */
SCIP_Bool SCIPrelaxationIsSolZero(
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   )
{
   assert(relaxation != NULL);

   return relaxation->relaxsolzero;
}

/** sets the relaxsolvalid and includeslp flags in the relaxation data to the given values */
void SCIPrelaxationSetSolValid(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Bool             isvalid,            /**< is the stored solution valid? */
   SCIP_Bool             includeslp          /**< does the relaxator contain all cuts in the LP? */
   )
{
   assert(relaxation != NULL);

   relaxation->relaxsolvalid = isvalid;
   relaxation->relaxsolincludeslp = includeslp;
}

/** returns whether the global relaxation solution is valid */
SCIP_Bool SCIPrelaxationIsSolValid(
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   )
{
   assert(relaxation != NULL);

   return relaxation->relaxsolvalid;
}

/** returns whether the global relaxation solution was computed by a relaxator which included all LP cuts */
SCIP_Bool SCIPrelaxationIsLpIncludedForSol(
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   )
{
   assert(relaxation != NULL);

   return relaxation->relaxsolincludeslp;
}

/** sets the objective value of the global relaxation solution */
void SCIPrelaxationSetSolObj(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Real             obj                 /**< objective value */
   )
{
   assert(relaxation != NULL);

   relaxation->relaxsolobjval = obj;
}

/** returns the objective value of the global relaxation solution w.r.t. the transformed problem */
SCIP_Real SCIPrelaxationGetSolObj(
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   )
{
   assert(relaxation != NULL);

   return relaxation->relaxsolobjval;
}

/** adds the given value to the global relaxation solution's objective value */
void SCIPrelaxationSolObjAdd(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Real             val                 /**< value to add to the objective value */
   )
{
   assert(relaxation != NULL);

   relaxation->relaxsolobjval += val;
}

/** updates objective value of current relaxation solution after change of objective coefficient */
void SCIPrelaxationUpdateVarObj(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable with changed objective coefficient */
   SCIP_Real             oldobj,             /**< old objective coefficient */
   SCIP_Real             newobj              /**< new objective coefficient */
   )
{
   SCIP_Real relaxsolval;

   assert(relaxation != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   relaxsolval = SCIPvarGetRelaxSol(var, set);
   relaxation->relaxsolobjval += (newobj - oldobj) * relaxsolval;
}
