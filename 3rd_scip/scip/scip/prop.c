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

/**@file   prop.c
 * @brief  methods and datastructures for propagators
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
#include "scip/var.h"
#include "scip/scip.h"
#include "scip/prop.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_prop.h"


/** compares two propagators w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPpropComp)
{  /*lint --e{715}*/
   return ((SCIP_PROP*)elem2)->priority - ((SCIP_PROP*)elem1)->priority;
}

/** compares two propagators w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPpropCompPresol)
{  /*lint --e{715}*/
   return ((SCIP_PROP*)elem2)->presolpriority - ((SCIP_PROP*)elem1)->presolpriority;
}

/** comparison method for sorting propagators w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPpropCompName)
{
   return strcmp(SCIPpropGetName((SCIP_PROP*)elem1), SCIPpropGetName((SCIP_PROP*)elem2));
}

/** method to call, when the priority of a propagator was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdPropPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPropPriority() to mark the props unsorted */
   SCIP_CALL( SCIPsetPropPriority(scip, (SCIP_PROP*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** method to call, when the presolving priority of a propagator was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdPropPresolPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPropPriority() to mark the props unsorted */
   SCIP_CALL( SCIPsetPropPresolPriority(scip, (SCIP_PROP*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given propagator to a new scip */
SCIP_RETCODE SCIPpropCopyInclude(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(prop != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( prop->propcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including propagator %s in subscip %p\n", SCIPpropGetName(prop), (void*)set->scip);
      SCIP_CALL( prop->propcopy(set->scip, prop) );
   }
   return SCIP_OKAY;
}

/** creates a propagator */
SCIP_RETCODE SCIPpropCreate(
   SCIP_PROP**           prop,               /**< pointer to propagator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagator should be executed */
   int                   presolpriority,     /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the propagator's presolving method */
   SCIP_DECL_PROPCOPY    ((*propcopy)),      /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   SCIP_DECL_PROPINITPRE ((*propinitpre)),   /**< presolving initialization method of propagator */
   SCIP_DECL_PROPEXITPRE ((*propexitpre)),   /**< presolving deinitialization method of propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(prop != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(propexec != NULL);

   /* the interface change from delay flags to timings cannot be recognized at compile time: Exit with an appropriate
    * error message
    */
   if( presoltiming < SCIP_PRESOLTIMING_NONE || presoltiming > SCIP_PRESOLTIMING_MAX )
   {
      SCIPmessagePrintError("ERROR: 'PRESOLDELAY'-flag no longer available since SCIP 3.2, use an appropriate "
         "'SCIP_PRESOLTIMING' for <%s> propagator instead.\n", name);

      return SCIP_PARAMETERWRONGVAL;
   }

   SCIP_ALLOC( BMSallocMemory(prop) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*prop)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*prop)->desc, desc, strlen(desc)+1) );
   (*prop)->priority = priority;
   (*prop)->freq = freq;
   (*prop)->propcopy = propcopy;
   (*prop)->propfree = propfree;
   (*prop)->propinit = propinit;
   (*prop)->propexit = propexit;
   (*prop)->propinitpre = propinitpre;
   (*prop)->propexitpre = propexitpre;
   (*prop)->propinitsol = propinitsol;
   (*prop)->propexitsol = propexitsol;
   (*prop)->proppresol = proppresol;
   (*prop)->propexec = propexec;
   (*prop)->propresprop = propresprop;
   (*prop)->propdata = propdata;
   SCIP_CALL( SCIPclockCreate(&(*prop)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*prop)->proptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*prop)->sbproptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*prop)->resproptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*prop)->presoltime, SCIP_CLOCKTYPE_DEFAULT) );
   (*prop)->ncalls = 0;
   (*prop)->nrespropcalls = 0;
   (*prop)->ncutoffs = 0;
   (*prop)->ndomredsfound = 0;
   (*prop)->wasdelayed = FALSE;
   (*prop)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of propagator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*prop)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4,
         paramChgdPropPriority, (SCIP_PARAMDATA*)(*prop)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/freq", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "frequency for calling propagator <%s> (-1: never, 0: only in root node)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*prop)->freq, FALSE, freq, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/delay", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname,
         "should propagator be delayed, if other propagators found reductions?",
         &(*prop)->delay, TRUE, delay, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/timingmask", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing when propagator should be called (%u:BEFORELP, %u:DURINGLPLOOP, %u:AFTERLPLOOP, %u:ALWAYS))",
      SCIP_PROPTIMING_BEFORELP, SCIP_PROPTIMING_DURINGLPLOOP, SCIP_PROPTIMING_AFTERLPLOOP, SCIP_PROPTIMING_ALWAYS);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)(&(*prop)->timingmask), TRUE, timingmask, (int) SCIP_PROPTIMING_BEFORELP, (int) SCIP_PROPTIMING_ALWAYS, NULL, NULL) ); /*lint !e713*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/presolpriority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "presolving priority of propagator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*prop)->presolpriority, TRUE, presolpriority, INT_MIN/4, INT_MAX/4,
         paramChgdPropPresolPriority, (SCIP_PARAMDATA*)(*prop)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname,
         "maximal number of presolving rounds the propagator participates in (-1: no limit)",
         &(*prop)->maxprerounds, FALSE, presolmaxrounds, -1, INT_MAX, NULL, NULL) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/presoltiming", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "timing mask of the presolving method of propagator <%s> (%u:FAST, %u:MEDIUM, %u:EXHAUSTIVE, %u:FINAL)",
      name, SCIP_PRESOLTIMING_FAST, SCIP_PRESOLTIMING_MEDIUM, SCIP_PRESOLTIMING_EXHAUSTIVE, SCIP_PRESOLTIMING_FINAL);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)&(*prop)->presoltiming, TRUE, (int)presoltiming, (int) SCIP_PRESOLTIMING_NONE, (int) SCIP_PRESOLTIMING_MAX, NULL, NULL) ); /*lint !e740*/


   return SCIP_OKAY;
}

/** calls destructor and frees memory of propagator */
SCIP_RETCODE SCIPpropFree(
   SCIP_PROP**           prop,               /**< pointer to propagator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(*prop != NULL);
   assert(!(*prop)->initialized);
   assert(set != NULL);

   /* call destructor of propagator */
   if( (*prop)->propfree != NULL )
   {
      SCIP_CALL( (*prop)->propfree(set->scip, *prop) );
   }

   SCIPclockFree(&(*prop)->presoltime);
   SCIPclockFree(&(*prop)->resproptime);
   SCIPclockFree(&(*prop)->sbproptime);
   SCIPclockFree(&(*prop)->proptime);
   SCIPclockFree(&(*prop)->setuptime);
   BMSfreeMemoryArray(&(*prop)->desc);
   BMSfreeMemoryArray(&(*prop)->name);
   BMSfreeMemory(prop);

   return SCIP_OKAY;
}

/** initializes propagator */
SCIP_RETCODE SCIPpropInit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   if( prop->initialized )
   {
      SCIPerrorMessage("propagator <%s> already initialized\n", prop->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(prop->proptime);
      SCIPclockReset(prop->sbproptime);
      SCIPclockReset(prop->resproptime);
      SCIPclockReset(prop->presoltime);
      SCIPclockReset(prop->setuptime);

      prop->ncalls = 0;
      prop->nrespropcalls = 0;
      prop->ncutoffs = 0;
      prop->ndomredsfound = 0;
      prop->lastnfixedvars = 0;
      prop->lastnaggrvars = 0;
      prop->lastnchgvartypes = 0;
      prop->lastnchgbds = 0;
      prop->lastnaddholes = 0;
      prop->lastndelconss = 0;
      prop->lastnaddconss = 0;
      prop->lastnupgdconss = 0;
      prop->lastnchgcoefs = 0;
      prop->lastnchgsides = 0;
      prop->nfixedvars = 0;
      prop->naggrvars = 0;
      prop->nchgvartypes = 0;
      prop->nchgbds = 0;
      prop->naddholes = 0;
      prop->ndelconss = 0;
      prop->naddconss = 0;
      prop->nupgdconss = 0;
      prop->nchgcoefs = 0;
      prop->nchgsides = 0;
      prop->npresolcalls = 0;
      prop->wasdelayed = FALSE;
   }

   if( prop->propinit != NULL )
   {
      /* start timing */
      SCIPclockStart(prop->setuptime, set);

      SCIP_CALL( prop->propinit(set->scip, prop) );

      /* stop timing */
      SCIPclockStop(prop->setuptime, set);
   }
   prop->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of propagator */
SCIP_RETCODE SCIPpropExit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   if( !prop->initialized )
   {
      SCIPerrorMessage("propagator <%s> not initialized\n", prop->name);
      return SCIP_INVALIDCALL;
   }

   if( prop->propexit != NULL )
   {
      /* start timing */
      SCIPclockStart(prop->setuptime, set);

      SCIP_CALL( prop->propexit(set->scip, prop) );

      /* stop timing */
      SCIPclockStop(prop->setuptime, set);
   }
   prop->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs propagator that the presolving process is being started */
SCIP_RETCODE SCIPpropInitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   prop->lastnfixedvars = 0;
   prop->lastnaggrvars = 0;
   prop->lastnchgvartypes = 0;
   prop->lastnchgbds = 0;
   prop->lastnaddholes = 0;
   prop->lastndelconss = 0;
   prop->lastnaddconss = 0;
   prop->lastnupgdconss = 0;
   prop->lastnchgcoefs = 0;
   prop->lastnchgsides = 0;
   prop->wasdelayed = FALSE;

   /* call presolving initialization method of propagator */
   if( prop->propinitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(prop->setuptime, set);

      SCIP_CALL( prop->propinitpre(set->scip, prop) );

      /* stop timing */
      SCIPclockStop(prop->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs propagator that the presolving process is finished */
SCIP_RETCODE SCIPpropExitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   /* call presolving deinitialization method of propagator */
   if( prop->propexitpre != NULL )
   {
      /* start timing */
      SCIPclockStart(prop->setuptime, set);

      SCIP_CALL( prop->propexitpre(set->scip, prop) );

      /* stop timing */
      SCIPclockStop(prop->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs propagator that the prop and bound process is being started */
SCIP_RETCODE SCIPpropInitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   /* call solving process initialization method of propagator */
   if( prop->propinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(prop->setuptime, set);

      SCIP_CALL( prop->propinitsol(set->scip, prop) );

      /* stop timing */
      SCIPclockStop(prop->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs propagator that the prop and bound process data is being freed */
SCIP_RETCODE SCIPpropExitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of propagator */
   if( prop->propexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(prop->setuptime, set);

      SCIP_CALL( prop->propexitsol(set->scip, prop, restart) );

      /* stop timing */
      SCIPclockStop(prop->setuptime, set);
   }

   return SCIP_OKAY;
}

/** executes presolving method of propagator */
SCIP_RETCODE SCIPpropPresol(
   SCIP_PROP*            prop,               /**< propagator */
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
   assert(prop != NULL);
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

   if( prop->proppresol == NULL )
      return SCIP_OKAY;

   /* check number of presolving rounds */
   if( prop->maxprerounds >= 0 && nrounds >= prop->maxprerounds )
      return SCIP_OKAY;

   /* check, if presolver should be delayed */
   if( prop->presoltiming & timing )
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

      SCIPsetDebugMsg(set, "calling presolving method of propagator <%s>\n", prop->name);

      /* calculate the number of changes since last call */
      nnewfixedvars = *nfixedvars - prop->lastnfixedvars;
      nnewaggrvars = *naggrvars - prop->lastnaggrvars;
      nnewchgvartypes = *nchgvartypes - prop->lastnchgvartypes;
      nnewchgbds = *nchgbds - prop->lastnchgbds;
      nnewaddholes = *naddholes - prop->lastnaddholes;
      nnewdelconss = *ndelconss - prop->lastndelconss;
      nnewaddconss = *naddconss - prop->lastnaddconss;
      nnewupgdconss = *nupgdconss - prop->lastnupgdconss;
      nnewchgcoefs = *nchgcoefs - prop->lastnchgcoefs;
      nnewchgsides = *nchgsides - prop->lastnchgsides;

      /* remember the number of changes prior to the call of the presolver method of the propagator */
      prop->lastnfixedvars = *nfixedvars;
      prop->lastnaggrvars = *naggrvars;
      prop->lastnchgvartypes = *nchgvartypes;
      prop->lastnchgbds = *nchgbds;
      prop->lastnaddholes = *naddholes;
      prop->lastndelconss = *ndelconss;
      prop->lastnaddconss = *naddconss;
      prop->lastnupgdconss = *nupgdconss;
      prop->lastnchgcoefs = *nchgcoefs;
      prop->lastnchgsides = *nchgsides;

      /* start timing */
      SCIPclockStart(prop->presoltime, set);

      /* call external method */
      SCIP_CALL( prop->proppresol(set->scip, prop, nrounds, timing,
            nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewaddholes,
            nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
            nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
            ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

      /* stop timing */
      SCIPclockStop(prop->presoltime, set);

      /* add/count the new changes */
      prop->nfixedvars += *nfixedvars - prop->lastnfixedvars;
      prop->naggrvars += *naggrvars - prop->lastnaggrvars;
      prop->nchgvartypes += *nchgvartypes - prop->lastnchgvartypes;
      prop->nchgbds += *nchgbds - prop->lastnchgbds;
      prop->naddholes += *naddholes - prop->lastnaddholes;
      prop->ndelconss += *ndelconss - prop->lastndelconss;
      prop->naddconss += *naddconss - prop->lastnaddconss;
      prop->nupgdconss += *nupgdconss - prop->lastnupgdconss;
      prop->nchgcoefs += *nchgcoefs - prop->lastnchgcoefs;
      prop->nchgsides += *nchgsides - prop->lastnchgsides;

      /* check result code of callback method */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_SUCCESS
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("propagator <%s> returned invalid result <%d>\n", prop->name, *result);
         return SCIP_INVALIDRESULT;
      }

      /* increase the number of presolving calls, if the propagator tried to find reductions */
      if( *result != SCIP_DIDNOTRUN )
         ++(prop->npresolcalls);
   }

   return SCIP_OKAY;
}

/** calls execution method of propagator */
SCIP_RETCODE SCIPpropExec(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute propagator even if it is marked to be delayed */
   SCIP_Bool             instrongbranching,  /**< are we currently doing strong branching? */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving process */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(prop != NULL);
   assert(prop->propexec != NULL);
   assert(prop->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( (depth == 0 && prop->freq == 0) || (prop->freq > 0 && depth % prop->freq == 0) )
   {
      if( !prop->delay || execdelayed )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;

         SCIPsetDebugMsg(set, "executing propagator <%s>\n", prop->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;

         /* start timing */
         if( instrongbranching )
            SCIPclockStart(prop->sbproptime, set);
         else
            SCIPclockStart(prop->proptime, set);

         /* call external propagation method */
         SCIP_CALL( prop->propexec(set->scip, prop, proptiming, result) );

         /* stop timing */
         if( instrongbranching )
            SCIPclockStop(prop->sbproptime, set);
         else
            SCIPclockStop(prop->proptime, set);

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
            prop->ncalls++;
         if( *result == SCIP_CUTOFF )
            prop->ncutoffs++;

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         prop->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         prop->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED
            && *result != SCIP_DELAYNODE )
         {
            SCIPerrorMessage("execution method of propagator <%s> returned invalid result <%d>\n", 
               prop->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPsetDebugMsg(set, "propagator <%s> was delayed\n", prop->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether propagator was delayed */
      prop->wasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}

/** resolves the given conflicting bound, that was deduced by the given propagator, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 *
 *  @note it is sufficient to explain the relaxed bound change
 */
SCIP_RETCODE SCIPpropResolvePropagation(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int                   inferinfo,          /**< user inference information attached to the bound change */
   SCIP_BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(prop != NULL);
   assert((inferboundtype == SCIP_BOUNDTYPE_LOWER
         && SCIPgetVarLbAtIndex(set->scip, infervar, bdchgidx, TRUE) > SCIPvarGetLbGlobal(infervar))
      || (inferboundtype == SCIP_BOUNDTYPE_UPPER
         && SCIPgetVarUbAtIndex(set->scip, infervar, bdchgidx, TRUE) < SCIPvarGetUbGlobal(infervar)));
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( prop->propresprop != NULL )
   {

      /* start timing */
      SCIPclockStart(prop->resproptime, set);

      SCIP_CALL( prop->propresprop(set->scip, prop, infervar, inferinfo, inferboundtype, bdchgidx,
            relaxedbd, result) );

      /* stop timing */
      SCIPclockStop(prop->resproptime, set);

      /* update statistic */
      prop->nrespropcalls++;

      /* check result code */
      if( *result != SCIP_SUCCESS && *result != SCIP_DIDNOTFIND )
      {
         SCIPerrorMessage("propagation conflict resolving method of propagator <%s> returned invalid result <%d>\n", 
            prop->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }
   else
   {
      SCIPerrorMessage("propagation conflict resolving method of propagator <%s> is not implemented\n", prop->name);
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** gets user data of propagator */
SCIP_PROPDATA* SCIPpropGetData(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->propdata;
}

/** sets user data of propagator; user has to free old data in advance! */
void SCIPpropSetData(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROPDATA*        propdata            /**< new propagator user data */
   )
{
   assert(prop != NULL);

   prop->propdata = propdata;
}

/** sets copy method of propagator */
void SCIPpropSetCopy(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPCOPY    ((*propcopy))       /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(prop != NULL);

   prop->propcopy = propcopy;
}

/** sets destructor method of propagator */
void SCIPpropSetFree(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPFREE    ((*propfree))       /**< destructor of propagator */
   )
{
   assert(prop != NULL);

   prop->propfree = propfree;
}

/** sets initialization method of propagator */
void SCIPpropSetInit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINIT    ((*propinit))       /**< initialize propagator */
   )
{
   assert(prop != NULL);

   prop->propinit = propinit;
}

/** sets deinitialization method of propagator */
void SCIPpropSetExit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXIT    ((*propexit))       /**< deinitialize propagator */
   )
{
   assert(prop != NULL);

   prop->propexit = propexit;
}

/** sets solving process initialization method of propagator */
void SCIPpropSetInitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITSOL((*propinitsol))     /**< solving process initialization method of propagator */
   )
{
   assert(prop != NULL);

   prop->propinitsol = propinitsol;
}

/** sets solving process deinitialization method of propagator */
void SCIPpropSetExitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol))    /**< solving process deinitialization method of propagator */
   )
{
   assert(prop != NULL);

   prop->propexitsol = propexitsol;
}

/** sets preprocessing initialization method of propagator */
void SCIPpropSetInitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITPRE((*propinitpre))     /**< preprocessing initialization method of propagator */
   )
{
   assert(prop != NULL);

   prop->propinitpre = propinitpre;
}



/** sets preprocessing deinitialization method of propagator */
void SCIPpropSetExitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITPRE((*propexitpre))     /**< preprocessing deinitialization method of propagator */
   )
{
   assert(prop != NULL);

   prop->propexitpre = propexitpre;
}

/** sets presolving method of propagator */
SCIP_RETCODE SCIPpropSetPresol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   int                   presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the propagator's presolving method */
   )
{
   assert(prop != NULL);

   prop->proppresol = proppresol;
   prop->presolpriority = presolpriority;
   /* the interface change from delay flags to timings cannot be recognized at compile time: Exit with an appropriate
    * error message
    */
   if( presoltiming < SCIP_PRESOLTIMING_FAST || presoltiming > SCIP_PRESOLTIMING_MAX )
   {
      SCIPmessagePrintError("ERROR: 'PRESOLDELAY'-flag no longer available since SCIP 3.2, use an appropriate "
         "'SCIP_PRESOLTIMING' for <%s> constraint handler instead.\n", prop->name);

      return SCIP_PARAMETERWRONGVAL;
   }

   prop->presoltiming = presoltiming;
   prop->maxprerounds = presolmaxrounds;

   return SCIP_OKAY;
}

/** sets propagation conflict resolving callback of propagator */
void SCIPpropSetResprop(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop))    /**< propagation conflict resolving callback */
   )
{
   assert(prop != NULL);

   prop->propresprop = propresprop;
}

/** gets name of propagator */
const char* SCIPpropGetName(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->name;
}

/** gets description of propagator */
const char* SCIPpropGetDesc(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->desc;
}

/** gets priority of propagator */
int SCIPpropGetPriority(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->priority;
}

/** gets presolving priority of propagator */
int SCIPpropGetPresolPriority(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->presolpriority;
}

/** sets priority of propagator */
void SCIPpropSetPriority(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the propagator */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   prop->priority = priority;
   set->propssorted = FALSE;
}

/** sets presolving priority of propagator */
void SCIPpropSetPresolPriority(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   presolpriority      /**< new priority of the propagator */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   prop->presolpriority = presolpriority;
   set->propspresolsorted = FALSE;
}

/** gets frequency of propagator */
int SCIPpropGetFreq(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->freq;
}

/** enables or disables all clocks of \p prop, depending on the value of the flag */
void SCIPpropEnableOrDisableClocks(
   SCIP_PROP*            prop,               /**< the propagator for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the propagator be enabled? */
   )
{
   assert(prop != NULL);

   SCIPclockEnableOrDisable(prop->setuptime, enable);
   SCIPclockEnableOrDisable(prop->presoltime, enable);
   SCIPclockEnableOrDisable(prop->proptime, enable);
   SCIPclockEnableOrDisable(prop->resproptime, enable);
   SCIPclockEnableOrDisable(prop->sbproptime, enable);
}

/** gets time in seconds used for setting up this propagator for new stages */
SCIP_Real SCIPpropGetSetupTime(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->setuptime);
}

/** sets frequency of propagator */
void SCIPpropSetFreq(
   SCIP_PROP*            prop,               /**< propagator */
   int                   freq                /**< new frequency of propagator */
   )
{
   assert(prop != NULL);
   assert(freq >= -1);

   prop->freq = freq;
}

/** gets time in seconds used in this propagator for propagation */
SCIP_Real SCIPpropGetTime(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->proptime);
}

/** gets time in seconds used in this propagator for propagation during strong branching */
SCIP_Real SCIPpropGetStrongBranchPropTime(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->sbproptime);
}

/** gets time in seconds used in this propagator for resolve propagation */
SCIP_Real SCIPpropGetRespropTime(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->resproptime);
}

/** gets time in seconds used in this propagator for presolving */
SCIP_Real SCIPpropGetPresolTime(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->presoltime);
}

/** gets the total number of times, the propagator was called */
SCIP_Longint SCIPpropGetNCalls(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ncalls;
}

/** gets the total number of times, the propagator was called for resolving a propagation */
SCIP_Longint SCIPpropGetNRespropCalls(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nrespropcalls;
}

/** gets total number of times, this propagator detected a cutoff */
SCIP_Longint SCIPpropGetNCutoffs(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ncutoffs;
}

/** gets total number of domain reductions found by this propagator */
SCIP_Longint SCIPpropGetNDomredsFound(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ndomredsfound;
}

/** should propagator be delayed, if other propagators found reductions? */
SCIP_Bool SCIPpropIsDelayed(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->delay;
}

/** was propagator delayed at the last call? */
SCIP_Bool SCIPpropWasDelayed(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->wasdelayed;
}

/** is propagator initialized? */
SCIP_Bool SCIPpropIsInitialized(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->initialized;
}

/** gets number of variables fixed during presolving of propagator */
int SCIPpropGetNFixedVars(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nfixedvars;
}

/** gets number of variables aggregated during presolving of propagator  */
int SCIPpropGetNAggrVars(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->naggrvars;
}

/** gets number of variable types changed during presolving of propagator  */
int SCIPpropGetNChgVarTypes(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nchgvartypes;
}

/** gets number of bounds changed during presolving of propagator  */
int SCIPpropGetNChgBds(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nchgbds;
}

/** gets number of holes added to domains of variables during presolving of propagator  */
int SCIPpropGetNAddHoles(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->naddholes;
}

/** gets number of constraints deleted during presolving of propagator */
int SCIPpropGetNDelConss(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ndelconss;
}

/** gets number of constraints added during presolving of propagator */
int SCIPpropGetNAddConss(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->naddconss;
}

/** gets number of constraints upgraded during presolving of propagator  */
int SCIPpropGetNUpgdConss(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nupgdconss;
}

/** gets number of coefficients changed during presolving of propagator */
int SCIPpropGetNChgCoefs(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nchgcoefs;
}

/** gets number of constraint sides changed during presolving of propagator */
int SCIPpropGetNChgSides(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nchgsides;
}

/** gets number of times the propagator was called in presolving and tried to find reductions */
int SCIPpropGetNPresolCalls(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->npresolcalls;
}

/** returns the timing mask of the propagator */
SCIP_PROPTIMING SCIPpropGetTimingmask(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->timingmask;
}

/** does the propagator perform presolving? */
SCIP_Bool SCIPpropDoesPresolve(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return (prop->proppresol != NULL);
}

/** returns the timing mask of the presolving method of the propagator */
SCIP_PRESOLTIMING SCIPpropGetPresolTiming(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->presoltiming;
}

/** sets the timing mask of the presolving method of the propagator */
void SCIPpropSetPresolTiming(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PRESOLTIMING     presoltiming        /** timing mask to be set */
   )
{
   assert(prop != NULL);

   prop->presoltiming = presoltiming;
}
