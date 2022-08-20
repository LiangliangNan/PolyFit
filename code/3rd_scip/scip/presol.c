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

/**@file   presol.c
 * @brief  methods for presolvers
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 * @todo add maxrounds parameter for single timings, count number of runs of a presolver with given timing
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/pub_misc.h"
#include "scip/presol.h"

#include "scip/struct_presol.h"



/*
 * presolver methods
 */

/** compares two presolvers w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPpresolComp)
{  /*lint --e{715}*/
   return ((SCIP_PRESOL*)elem2)->priority - ((SCIP_PRESOL*)elem1)->priority;
}

/** comparison method for sorting presolvers w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPpresolCompName)
{
   return strcmp(SCIPpresolGetName((SCIP_PRESOL*)elem1), SCIPpresolGetName((SCIP_PRESOL*)elem2));
}

/** method to call, when the priority of a presolver was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdPresolPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPresolPriority() to mark the presols unsorted */
   SCIP_CALL( SCIPsetPresolPriority(scip, (SCIP_PRESOL*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given presolver to a new scip */
SCIP_RETCODE SCIPpresolCopyInclude(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(presol != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( presol->presolcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including presolver %s in subscip %p\n", SCIPpresolGetName(presol), (void*)set->scip);
      SCIP_CALL( presol->presolcopy(set->scip, presol) );
   }
   return SCIP_OKAY;
}

/** creates a presolver */
SCIP_RETCODE SCIPpresolCreate(
   SCIP_PRESOL**         presol,             /**< pointer to store presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of presolver */
   const char*           desc,               /**< description of presolver */
   int                   priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int                   maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   SCIP_PRESOLTIMING     timing,             /**< timing mask of the presolver */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy)),    /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   SCIP_DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   SCIP_DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   SCIP_DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   SCIP_DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   SCIP_DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(presol != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   /* the interface change from delay flags to timings cannot be recognized at compile time: Exit with an appropriate
    * error message
    */
   if( timing < SCIP_PRESOLTIMING_NONE || timing > SCIP_PRESOLTIMING_MAX )
   {
      SCIPmessagePrintError("ERROR: 'PRESOLDELAY'-flag no longer available since SCIP 3.2, use an appropriate "
         "'SCIP_PRESOLTIMING' for <%s> presolver instead.\n", name);

      return SCIP_PARAMETERWRONGVAL;
   }

   SCIP_ALLOC( BMSallocMemory(presol) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*presol)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*presol)->desc, desc, strlen(desc)+1) );
   (*presol)->presolcopy = presolcopy;
   (*presol)->presolfree = presolfree;
   (*presol)->presolinit = presolinit;
   (*presol)->presolexit = presolexit;
   (*presol)->presolinitpre = presolinitpre;
   (*presol)->presolexitpre = presolexitpre;
   (*presol)->presolexec = presolexec;
   (*presol)->presoldata = presoldata;
   SCIP_CALL( SCIPclockCreate(&(*presol)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*presol)->presolclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*presol)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "presolving/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of presolver <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*presol)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4,
         paramChgdPresolPriority, (SCIP_PARAMDATA*)(*presol)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "presolving/%s/maxrounds", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "maximal number of presolving rounds the presolver participates in (-1: no limit)",
         &(*presol)->maxrounds, FALSE, maxrounds, -1, INT_MAX, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "presolving/%s/timing", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing mask of presolver <%s> (%u:FAST, %u:MEDIUM, %u:EXHAUSTIVE, %u:FINAL)",
      name, SCIP_PRESOLTIMING_FAST, SCIP_PRESOLTIMING_MEDIUM, SCIP_PRESOLTIMING_EXHAUSTIVE, SCIP_PRESOLTIMING_FINAL);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)&(*presol)->timing, TRUE, (int)timing, (int) SCIP_PRESOLTIMING_FAST, (int) SCIP_PRESOLTIMING_MAX, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** frees memory of presolver */
SCIP_RETCODE SCIPpresolFree(
   SCIP_PRESOL**         presol,             /**< pointer to presolver data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(*presol != NULL);
   assert(!(*presol)->initialized);
   assert(set != NULL);

   /* call destructor of presolver */
   if( (*presol)->presolfree != NULL )
   {
      SCIP_CALL( (*presol)->presolfree(set->scip, *presol) );
   }

   SCIPclockFree(&(*presol)->presolclock);
   SCIPclockFree(&(*presol)->setuptime);
   BMSfreeMemoryArray(&(*presol)->name);
   BMSfreeMemoryArray(&(*presol)->desc);
   BMSfreeMemory(presol);

   return SCIP_OKAY;
}

/** initializes presolver */
SCIP_RETCODE SCIPpresolInit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(set != NULL);

   if( presol->initialized )
   {
      SCIPerrorMessage("presolver <%s> already initialized\n", presol->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(presol->setuptime);
      SCIPclockReset(presol->presolclock);

      presol->lastnfixedvars = 0;
      presol->lastnaggrvars = 0;
      presol->lastnchgvartypes = 0;
      presol->lastnchgbds = 0;
      presol->lastnaddholes = 0;
      presol->lastndelconss = 0;
      presol->lastnaddconss = 0;
      presol->lastnupgdconss = 0;
      presol->lastnchgcoefs = 0;
      presol->lastnchgsides = 0;
      presol->nfixedvars = 0;
      presol->naggrvars = 0;
      presol->nchgvartypes = 0;
      presol->nchgbds = 0;
      presol->naddholes = 0;
      presol->ndelconss = 0;
      presol->naddconss = 0;
      presol->nupgdconss = 0;
      presol->nchgcoefs = 0;
      presol->nchgsides = 0;
      presol->ncalls = 0;
   }

   /* call initialization method of presolver */
   if( presol->presolinit != NULL )
   {
      /* start timing */
      SCIPclockStart(presol->setuptime, set);

      SCIP_CALL( presol->presolinit(set->scip, presol) );

      /* stop timing */
      SCIPclockStop(presol->setuptime, set);
   }
   presol->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes presolver */
SCIP_RETCODE SCIPpresolExit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(set != NULL);

   if( !presol->initialized )
   {
      SCIPerrorMessage("presolver <%s> not initialized\n", presol->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of presolver */
   if( presol->presolexit != NULL )
   {
      /* start timing */
      SCIPclockStart(presol->setuptime, set);

      SCIP_CALL( presol->presolexit(set->scip, presol) );

      /* stop timing */
      SCIPclockStop(presol->setuptime, set);
   }
   presol->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs presolver that the presolving process is being started */
SCIP_RETCODE SCIPpresolInitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(set != NULL);

   presol->lastnfixedvars = 0;
   presol->lastnaggrvars = 0;
   presol->lastnchgvartypes = 0;
   presol->lastnchgbds = 0;
   presol->lastnaddholes = 0;
   presol->lastndelconss = 0;
   presol->lastnaddconss = 0;
   presol->lastnupgdconss = 0;
   presol->lastnchgcoefs = 0;
   presol->lastnchgsides = 0;

   /* call presolving initialization method of presolver */
   if( presol->presolinitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(presol->setuptime, set);

      SCIP_CALL( presol->presolinitpre(set->scip, presol) );

      /* stop timing */
      SCIPclockStop(presol->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs presolver that the presolving process is finished */
SCIP_RETCODE SCIPpresolExitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(set != NULL);

   /* call presolving deinitialization method of presolver */
   if( presol->presolexitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(presol->setuptime, set);

      SCIP_CALL( presol->presolexitpre(set->scip, presol) );

      /* stop timing */
      SCIPclockStop(presol->setuptime, set);
   }

   return SCIP_OKAY;
}

/** executes presolver */
SCIP_RETCODE SCIPpresolExec(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRESOLTIMING     timing,             /**< current presolving timing */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   int nnewfixedvars;
   int nnewaggrvars;
   int nnewchgvartypes;
   int nnewchgbds;
   int nnewaddholes;
   int nnewdelconss;
   int nnewaddconss;
   int nnewupgdconss;
   int nnewchgcoefs;
   int nnewchgsides;

   assert(presol != NULL);
   assert(presol->presolexec != NULL);
   assert(set != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgvartypes != NULL);
   assert(nchgbds != NULL);
   assert(naddholes != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check number of presolving rounds */
   if( presol->maxrounds >= 0 && nrounds >= presol->maxrounds )
      return SCIP_OKAY;

   /* calculate the number of changes since last call */
   nnewfixedvars = *nfixedvars - presol->lastnfixedvars;
   nnewaggrvars = *naggrvars - presol->lastnaggrvars;
   nnewchgvartypes = *nchgvartypes - presol->lastnchgvartypes;
   nnewchgbds = *nchgbds - presol->lastnchgbds;
   nnewaddholes = *naddholes - presol->lastnaddholes;
   nnewdelconss = *ndelconss - presol->lastndelconss;
   nnewaddconss = *naddconss - presol->lastnaddconss;
   nnewupgdconss = *nupgdconss - presol->lastnupgdconss;
   nnewchgcoefs = *nchgcoefs - presol->lastnchgcoefs;
   nnewchgsides = *nchgsides - presol->lastnchgsides;

   /* remember the number of changes prior to the call of the presolver */
   presol->lastnfixedvars = *nfixedvars;
   presol->lastnaggrvars = *naggrvars;
   presol->lastnchgvartypes = *nchgvartypes;
   presol->lastnchgbds = *nchgbds;
   presol->lastnaddholes = *naddholes;
   presol->lastndelconss = *ndelconss;
   presol->lastnaddconss = *naddconss;
   presol->lastnupgdconss = *nupgdconss;
   presol->lastnchgcoefs = *nchgcoefs;
   presol->lastnchgsides = *nchgsides;

   /* check, if presolver should be called with the current timing */
   if( timing & presol->timing )
   {
      SCIPsetDebugMsg(set, "calling presolver <%s> with timing %u\n", presol->name, timing);

      /* start timing */
      SCIPclockStart(presol->presolclock, set);

      /* call external method */
      SCIP_CALL( presol->presolexec(set->scip, presol, nrounds, timing,
            nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewaddholes,
            nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
            nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
            ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

      /* stop timing */
      SCIPclockStop(presol->presolclock, set);

      /* add/count the new changes */
      presol->nfixedvars += *nfixedvars - presol->lastnfixedvars;
      presol->naggrvars += *naggrvars - presol->lastnaggrvars;
      presol->nchgvartypes += *nchgvartypes - presol->lastnchgvartypes;
      presol->nchgbds += *nchgbds - presol->lastnchgbds;
      presol->naddholes += *naddholes - presol->lastnaddholes;
      presol->ndelconss += *ndelconss - presol->lastndelconss;
      presol->naddconss += *naddconss - presol->lastnaddconss;
      presol->nupgdconss += *nupgdconss - presol->lastnupgdconss;
      presol->nchgcoefs += *nchgcoefs - presol->lastnchgcoefs;
      presol->nchgsides += *nchgsides - presol->lastnchgsides;

      /* check result code of callback method */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_SUCCESS
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("presolver <%s> returned invalid result <%d>\n", presol->name, *result);
         return SCIP_INVALIDRESULT;
      }

      /* increase the number of calls, if the presolver tried to find reductions */
      if( *result != SCIP_DIDNOTRUN )
         ++(presol->ncalls);
   }

   return SCIP_OKAY;
}

/** gets user data of presolver */
SCIP_PRESOLDATA* SCIPpresolGetData(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->presoldata;
}

/** sets user data of presolver; user has to free old data in advance! */
void SCIPpresolSetData(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_PRESOLDATA*      presoldata          /**< new presolver user data */
   )
{
   assert(presol != NULL);

   presol->presoldata = presoldata;
}

/** sets copy method of presolver */
void SCIPpresolSetCopy(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy))     /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(presol != NULL);

   presol->presolcopy = presolcopy;
}

/** sets destructor method of presolver */
void SCIPpresolSetFree(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLFREE  ((*presolfree))     /**< destructor of presolver */
   )
{
   assert(presol != NULL);

   presol->presolfree = presolfree;
}

/** sets initialization method of presolver */
void SCIPpresolSetInit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINIT  ((*presolinit))     /**< initialize presolver */
   )
{
   assert(presol != NULL);

   presol->presolinit = presolinit;
}

/** sets deinitialization method of presolver */
void SCIPpresolSetExit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXIT  ((*presolexit))     /**< deinitialize presolver */
   )
{
   assert(presol != NULL);

   presol->presolexit = presolexit;
}

/** sets solving process initialization method of presolver */
void SCIPpresolSetInitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLINITPRE ((*presolinitpre))/**< solving process initialization method of presolver */
   )
{
   assert(presol != NULL);

   presol->presolinitpre = presolinitpre;
}

/** sets solving process deinitialization method of presolver */
void SCIPpresolSetExitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_DECL_PRESOLEXITPRE ((*presolexitpre))/**< solving process deinitialization method of presolver */
   )
{
   assert(presol != NULL);

   presol->presolexitpre = presolexitpre;
}

/** gets name of presolver */
const char* SCIPpresolGetName(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->name;
}

/** gets description of presolver */
const char* SCIPpresolGetDesc(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->desc;
}

/** gets priority of presolver */
int SCIPpresolGetPriority(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->priority;
}

/** gets round limit of presolver */
int SCIPpresolGetMaxrounds(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->maxrounds;
}

/** sets priority of presolver */
void SCIPpresolSetPriority(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the presolver */
   )
{
   assert(presol != NULL);
   assert(set != NULL);

   presol->priority = priority;
   set->presolssorted = FALSE;
}

/** gets the timing mask of the presolver */
SCIP_PRESOLTIMING SCIPpresolGetTiming(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->timing;
}

/** sets the timing mask of the presolver */
void SCIPpresolSetTiming(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_PRESOLTIMING     timing              /**< timing mask of the presolver */
   )
{
   assert(presol != NULL);

   presol->timing = timing;
}


/** is presolver initialized? */
SCIP_Bool SCIPpresolIsInitialized(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->initialized;
}

/** enables or disables all clocks of \p presol, depending on the value of the flag */
void SCIPpresolEnableOrDisableClocks(
   SCIP_PRESOL*          presol,             /**< the presolver for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the presolver be enabled? */
   )
{
   assert(presol != NULL);

   SCIPclockEnableOrDisable(presol->setuptime, enable);
   SCIPclockEnableOrDisable(presol->presolclock, enable);
}

/** gets time in seconds used in this presolver for setting up for next stages */
SCIP_Real SCIPpresolGetSetupTime(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return SCIPclockGetTime(presol->setuptime);
}

/** gets time in seconds used in this presolver */
SCIP_Real SCIPpresolGetTime(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return SCIPclockGetTime(presol->presolclock);
}

/** gets number of variables fixed in presolver */
int SCIPpresolGetNFixedVars(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nfixedvars;
}

/** gets number of variables aggregated in presolver */
int SCIPpresolGetNAggrVars(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->naggrvars;
}

/** gets number of variable types changed in presolver */
int SCIPpresolGetNChgVarTypes(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgvartypes;
}

/** gets number of bounds changed in presolver */
int SCIPpresolGetNChgBds(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgbds;
}

/** gets number of holes added to domains of variables in presolver */
int SCIPpresolGetNAddHoles(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->naddholes;
}

/** gets number of constraints deleted in presolver */
int SCIPpresolGetNDelConss(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->ndelconss;
}

/** gets number of constraints added in presolver */
int SCIPpresolGetNAddConss(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->naddconss;
}

/** gets number of constraints upgraded in presolver */
int SCIPpresolGetNUpgdConss(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nupgdconss;
}

/** gets number of coefficients changed in presolver */
int SCIPpresolGetNChgCoefs(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgcoefs;
}

/** gets number of constraint sides changed in presolver */
int SCIPpresolGetNChgSides(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgsides;
}

/** gets number of times the presolver was called and tried to find reductions */
int SCIPpresolGetNCalls(
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->ncalls;
}
