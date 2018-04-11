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

/**@file   pricer.c
 * @brief  methods for variable pricers
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pricestore.h"
#include "scip/scip.h"
#include "scip/pricer.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"

#include "scip/struct_pricer.h"



/** compares two pricers w. r. to their activity and their priority */
SCIP_DECL_SORTPTRCOMP(SCIPpricerComp)
{  /*lint --e{715}*/
   if( ((SCIP_PRICER*)elem1)->active != ((SCIP_PRICER*)elem2)->active )
      return ((SCIP_PRICER*)elem1)->active ? -1 : +1;
   else
      return ((SCIP_PRICER*)elem2)->priority - ((SCIP_PRICER*)elem1)->priority;
}

/** comparison method for sorting pricers w.r.t. to their name */
SCIP_DECL_SORTPTRCOMP(SCIPpricerCompName)
{
   if( ((SCIP_PRICER*)elem1)->active != ((SCIP_PRICER*)elem2)->active )
      return ((SCIP_PRICER*)elem1)->active ? -1 : +1;
   else
      return strcmp(SCIPpricerGetName((SCIP_PRICER*)elem1), SCIPpricerGetName((SCIP_PRICER*)elem2));
}

/** method to call, when the priority of a pricer was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdPricerPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPricerPriority() to mark the pricers unsorted */
   SCIP_CALL( SCIPsetPricerPriority(scip, (SCIP_PRICER*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given pricer to a new scip */
SCIP_RETCODE SCIPpricerCopyInclude(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_SET*             set,                /**< SCIP_SET of SCIP to copy to */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);
   assert(valid != NULL);
   assert(set->scip != NULL);

   if( pricer->pricercopy != NULL )
   {
      SCIPsetDebugMsg(set, "including pricer %s in subscip %p\n", SCIPpricerGetName(pricer), (void*)set->scip);
      SCIP_CALL( pricer->pricercopy(set->scip, pricer, valid) );
   }
   return SCIP_OKAY;
}

/** creates a variable pricer
 *  To use the variable pricer for solving a problem, it first has to be activated with a call to SCIPactivatePricer().
 */
SCIP_RETCODE SCIPpricerCreate(
   SCIP_PRICER**         pricer,             /**< pointer to variable pricer data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of variable pricer */
   const char*           desc,               /**< description of variable pricer */
   int                   priority,           /**< priority of the variable pricer */
   SCIP_Bool             delay,              /**< should the pricer be delayed until no other pricers or already existing
                                              *   problem variables with negative reduced costs are found */
   SCIP_DECL_PRICERCOPY  ((*pricercopy)),    /**< copy method of pricer or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRICERFREE  ((*pricerfree)),    /**< destructor of variable pricer */
   SCIP_DECL_PRICERINIT  ((*pricerinit)),    /**< initialize variable pricer */
   SCIP_DECL_PRICEREXIT  ((*pricerexit)),    /**< deinitialize variable pricer */
   SCIP_DECL_PRICERINITSOL((*pricerinitsol)),/**< solving process initialization method of variable pricer */
   SCIP_DECL_PRICEREXITSOL((*pricerexitsol)),/**< solving process deinitialization method of variable pricer */
   SCIP_DECL_PRICERREDCOST((*pricerredcost)),/**< reduced cost pricing method of variable pricer for feasible LPs */
   SCIP_DECL_PRICERFARKAS((*pricerfarkas)),  /**< Farkas pricing method of variable pricer for infeasible LPs */
   SCIP_PRICERDATA*      pricerdata          /**< variable pricer data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(pricer != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(pricerredcost != NULL);

   SCIP_ALLOC( BMSallocMemory(pricer) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*pricer)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*pricer)->desc, desc, strlen(desc)+1) );
   (*pricer)->priority = priority;
   (*pricer)->pricercopy = pricercopy;
   (*pricer)->pricerfree = pricerfree;
   (*pricer)->pricerinit = pricerinit;
   (*pricer)->pricerexit = pricerexit;
   (*pricer)->pricerinitsol = pricerinitsol;
   (*pricer)->pricerexitsol = pricerexitsol;
   (*pricer)->pricerredcost = pricerredcost;
   (*pricer)->pricerfarkas = pricerfarkas;
   (*pricer)->pricerdata = pricerdata;
   SCIP_CALL( SCIPclockCreate(&(*pricer)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*pricer)->pricerclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*pricer)->ncalls = 0;
   (*pricer)->nvarsfound = 0;
   (*pricer)->delay = delay;
   (*pricer)->active = FALSE;
   (*pricer)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "pricers/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of pricer <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
                  &(*pricer)->priority, FALSE, priority, INT_MIN/4, INT_MAX/4,
                  paramChgdPricerPriority, (SCIP_PARAMDATA*)(*pricer)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of variable pricer */
SCIP_RETCODE SCIPpricerFree(
   SCIP_PRICER**         pricer,             /**< pointer to variable pricer data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(*pricer != NULL);
   assert(!(*pricer)->initialized);
   assert(set != NULL);

   /* call destructor of variable pricer */
   if( (*pricer)->pricerfree != NULL )
   {
      SCIP_CALL( (*pricer)->pricerfree(set->scip, *pricer) );
   }

   SCIPclockFree(&(*pricer)->pricerclock);
   SCIPclockFree(&(*pricer)->setuptime);
   BMSfreeMemoryArray(&(*pricer)->name);
   BMSfreeMemoryArray(&(*pricer)->desc);
   BMSfreeMemory(pricer);

   return SCIP_OKAY;
}

/** initializes variable pricer */
SCIP_RETCODE SCIPpricerInit(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(pricer->active);
   assert(set != NULL);

   if( pricer->initialized )
   {
      SCIPerrorMessage("variable pricer <%s> already initialized\n", pricer->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(pricer->setuptime);
      SCIPclockReset(pricer->pricerclock);

      pricer->ncalls = 0;
      pricer->nvarsfound = 0;
   }

   if( pricer->pricerinit != NULL )
   {
      /* start timing */
      SCIPclockStart(pricer->setuptime, set);

      SCIP_CALL( pricer->pricerinit(set->scip, pricer) );

      /* stop timing */
      SCIPclockStop(pricer->setuptime, set);
   }
   pricer->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of variable pricer */
SCIP_RETCODE SCIPpricerExit(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(pricer->active);
   assert(set != NULL);

   if( !pricer->initialized )
   {
      SCIPerrorMessage("variable pricer <%s> not initialized\n", pricer->name);
      return SCIP_INVALIDCALL;
   }

   if( pricer->pricerexit != NULL )
   {
      /* start timing */
      SCIPclockStart(pricer->setuptime, set);

      SCIP_CALL( pricer->pricerexit(set->scip, pricer) );

      /* stop timing */
      SCIPclockStop(pricer->setuptime, set);
   }
   pricer->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs variable pricer that the branch and bound process is being started */
SCIP_RETCODE SCIPpricerInitsol(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);

   /* call solving process initialization method of variable pricer */
   if( pricer->pricerinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(pricer->setuptime, set);

      SCIP_CALL( pricer->pricerinitsol(set->scip, pricer) );

      /* stop timing */
      SCIPclockStop(pricer->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs variable pricer that the branch and bound process data is being freed */
SCIP_RETCODE SCIPpricerExitsol(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of variable pricer */
   if( pricer->pricerexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(pricer->setuptime, set);

      SCIP_CALL( pricer->pricerexitsol(set->scip, pricer) );

      /* stop timing */
      SCIPclockStop(pricer->setuptime, set);
   }

   return SCIP_OKAY;
}

/** activates pricer such that it is called in LP solving loop */
SCIP_RETCODE SCIPpricerActivate(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_PROBLEM);

   if( !pricer->active )
   {
      pricer->active = TRUE;
      set->nactivepricers++;
      set->pricerssorted = FALSE;
   }

   return SCIP_OKAY;
}

/** deactivates pricer such that it is no longer called in LP solving loop */
SCIP_RETCODE SCIPpricerDeactivate(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_PROBLEM);

   if( pricer->active )
   {
      pricer->active = FALSE;
      set->nactivepricers--;
      set->pricerssorted = FALSE;
   }

   return SCIP_OKAY;
}

/** calls reduced cost pricing method of variable pricer */
SCIP_RETCODE SCIPpricerRedcost(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem */
   SCIP_Real*            lowerbound,         /**< local lower bound computed by the pricer */
   SCIP_Bool*            stopearly,          /**< should pricing be stopped, although new variables were added? */
   SCIP_RESULT*          result              /**< result of the pricing process */    
   )
{
   int oldnvars;

   assert(pricer != NULL);
   assert(pricer->active);
   assert(pricer->pricerredcost != NULL);
   assert(set != NULL);
   assert(prob != NULL);
   assert(lowerbound != NULL);
   assert(result != NULL);

   SCIPsetDebugMsg(set, "executing reduced cost pricing of variable pricer <%s>\n", pricer->name);

   oldnvars = prob->nvars;

   /* start timing */
   SCIPclockStart(pricer->pricerclock, set);

   /* call external method */
   SCIP_CALL( pricer->pricerredcost(set->scip, pricer, lowerbound, stopearly, result) );

   /* stop timing */
   SCIPclockStop(pricer->pricerclock, set);

   /* evaluate result */
   pricer->ncalls++;
   pricer->nvarsfound += prob->nvars - oldnvars;

   return SCIP_OKAY;
}

/** calls Farkas pricing method of variable pricer */
SCIP_RETCODE SCIPpricerFarkas(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem */
   SCIP_RESULT*          result              /**< result of the pricing process */
   )
{
   int oldnvars;

   assert(pricer != NULL);
   assert(pricer->active);
   assert(set != NULL);
   assert(prob != NULL);

   /* check, if pricer implemented a Farkas pricing algorithm */
   if( pricer->pricerfarkas == NULL )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "executing Farkas pricing of variable pricer <%s>\n", pricer->name);

   oldnvars = prob->nvars;

   /* start timing */
   SCIPclockStart(pricer->pricerclock, set);

   /* call external method */
   SCIP_CALL( pricer->pricerfarkas(set->scip, pricer, result) );

   /* stop timing */
   SCIPclockStop(pricer->pricerclock, set);

   /* evaluate result */
   pricer->ncalls++;
   pricer->nvarsfound += prob->nvars - oldnvars;

   return SCIP_OKAY;
}

/** depending on the LP's solution status, calls reduced cost or Farkas pricing method of variable pricer */
SCIP_RETCODE SCIPpricerExec(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_PRICESTORE*      pricestore,         /**< pricing storage */
   SCIP_Real*            lowerbound,         /**< local lower bound computed by the pricer */
   SCIP_Bool*            stopearly,          /**< should pricing be stopped, although new variables were added? */
   SCIP_RESULT*          result              /**< result of the pricing process */
   )
{
   assert(pricer != NULL);
   assert(lowerbound != NULL);
   assert(stopearly != NULL);
   assert(result != NULL);

   /* set lowerbound, stopearly, and result pointer */
   *lowerbound = - SCIPsetInfinity(set);
   *stopearly = FALSE;
   *result = SCIP_SUCCESS;

   /* check if pricer should be delayed */
   if( pricer->delay && SCIPpricestoreGetNVars(pricestore) > 0 )
      return SCIP_OKAY;

   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      SCIP_CALL( SCIPpricerFarkas(pricer, set, prob, result) );
   }
   else
   {
      *result = SCIP_DIDNOTRUN;
      SCIP_CALL( SCIPpricerRedcost(pricer, set, prob, lowerbound, stopearly, result) );
   }

   return SCIP_OKAY;
}

/** gets user data of variable pricer */
SCIP_PRICERDATA* SCIPpricerGetData(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->pricerdata;
}

/** sets user data of variable pricer; user has to free old data in advance! */
void SCIPpricerSetData(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_PRICERDATA*      pricerdata          /**< new variable pricer user data */
   )
{
   assert(pricer != NULL);

   pricer->pricerdata = pricerdata;
}

/** sets copy callback of pricer */
void SCIPpricerSetCopy(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_DECL_PRICERCOPY  ((*pricercopy))     /**< copy callback of pricer */
   )
{
   assert(pricer != NULL);

   pricer->pricercopy = pricercopy;
}

/** sets destructor callback of pricer */
void SCIPpricerSetFree(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICERFREE  ((*pricerfree))     /**< destructor of pricer */
   )
{
   assert(pricer != NULL);

   pricer->pricerfree = pricerfree;
}

/** sets initialization callback of pricer */
void SCIPpricerSetInit(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICERINIT ((*pricerinit))     /**< initialize pricer */
   )
{
   assert(pricer != NULL);

   pricer->pricerinit = pricerinit;
}

/** sets deinitialization callback of pricer */
void SCIPpricerSetExit(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICEREXIT ((*pricerexit))     /**< deinitialize pricer */
   )
{
   assert(pricer != NULL);

   pricer->pricerexit = pricerexit;
}

/** sets solving process initialization callback of pricer */
void SCIPpricerSetInitsol(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICERINITSOL ((*pricerinitsol))/**< solving process initialization callback of pricer */
   )
{
   assert(pricer != NULL);

   pricer->pricerinitsol = pricerinitsol;
}

/** sets solving process deinitialization callback of pricer */
void SCIPpricerSetExitsol(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICEREXITSOL ((*pricerexitsol))/**< solving process deinitialization callback of pricer */
   )
{
   assert(pricer != NULL);

   pricer->pricerexitsol = pricerexitsol;
}

/** gets name of variable pricer */
const char* SCIPpricerGetName(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->name;
}

/** gets description of variable pricer */
const char* SCIPpricerGetDesc(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->desc;
}

/** gets priority of variable pricer */
int SCIPpricerGetPriority(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->priority;
}

/** sets priority of variable pricer */
void SCIPpricerSetPriority(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the variable pricer */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);

   pricer->priority = priority;
   set->pricerssorted = FALSE;
}

/** gets the number of times, the pricer was called and tried to find a variable with negative reduced costs */
int SCIPpricerGetNCalls(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->ncalls;
}

/** gets the number of variables with negative reduced costs found by this pricer */
int SCIPpricerGetNVarsFound(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->nvarsfound;
}

/** gets time in seconds used in this pricer for setting up for next stages */
SCIP_Real SCIPpricerGetSetupTime(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return SCIPclockGetTime(pricer->setuptime);
}

/** gets time in seconds used in this pricer */
SCIP_Real SCIPpricerGetTime(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return SCIPclockGetTime(pricer->pricerclock);
}

/** enables or disables all clocks of \p pricer, depending on the value of the flag */
void SCIPpricerEnableOrDisableClocks(
   SCIP_PRICER*          pricer,             /**< the pricer for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the pricer be enabled? */
   )
{
   assert(pricer != NULL);

   SCIPclockEnableOrDisable(pricer->setuptime, enable);
   SCIPclockEnableOrDisable(pricer->pricerclock, enable);
}

/** returns whether the given pricer is in use in the current problem */
SCIP_Bool SCIPpricerIsActive(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->active;
}

/** returns whether the pricer should be delayed until no other pricer finds a new variable */
SCIP_Bool SCIPpricerIsDelayed(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->delay;
}

/** is variable pricer initialized? */
SCIP_Bool SCIPpricerIsInitialized(
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->initialized;
}


