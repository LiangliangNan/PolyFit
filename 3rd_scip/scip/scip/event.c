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

/**@file   event.c
 * @brief  methods and datastructures for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/clock.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/primal.h"
#include "scip/branch.h"
#include "scip/pub_message.h"

/* timing the execution methods for event handling takes a lot of time, so it is disabled */
/* #define TIMEEVENTEXEC */


/*
 * Event handler methods
 */

/** copies the given event handler to a new scip */
SCIP_RETCODE SCIPeventhdlrCopyInclude(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( eventhdlr->eventcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including event handler %s in subscip %p\n", SCIPeventhdlrGetName(eventhdlr), (void*)set->scip);
      SCIP_CALL( eventhdlr->eventcopy(set->scip, eventhdlr) );
   }

   return SCIP_OKAY;
}

/** creates an event handler */
SCIP_RETCODE SCIPeventhdlrCreate(
   SCIP_EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   const char*           name,               /**< name of event handler */
   const char*           desc,               /**< description of event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy)),     /**< copy method of event handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit)),     /**< initialize event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialize event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol)),  /**< solving process initialization method of event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol)),  /**< solving process deinitialization method of event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   SCIP_DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   assert(eventhdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(eventexec != NULL);

   SCIP_ALLOC( BMSallocMemory(eventhdlr) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*eventhdlr)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*eventhdlr)->desc, desc, strlen(desc)+1) );
   (*eventhdlr)->eventcopy = eventcopy;
   (*eventhdlr)->eventfree = eventfree;
   (*eventhdlr)->eventinit = eventinit;
   (*eventhdlr)->eventexit = eventexit;
   (*eventhdlr)->eventinitsol = eventinitsol;
   (*eventhdlr)->eventexitsol = eventexitsol;
   (*eventhdlr)->eventdelete = eventdelete;
   (*eventhdlr)->eventexec = eventexec;
   (*eventhdlr)->eventhdlrdata = eventhdlrdata;
   (*eventhdlr)->initialized = FALSE;

   /* create clocks */
   SCIP_CALL( SCIPclockCreate(&(*eventhdlr)->setuptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*eventhdlr)->eventtime, SCIP_CLOCKTYPE_DEFAULT) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of event handler */
SCIP_RETCODE SCIPeventhdlrFree(
   SCIP_EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(*eventhdlr != NULL);
   assert(!(*eventhdlr)->initialized);
   assert(set != NULL);

   /* call destructor of event handler */
   if( (*eventhdlr)->eventfree != NULL )
   {
      SCIP_CALL( (*eventhdlr)->eventfree(set->scip, *eventhdlr) );
   }

   /* free clocks */
   SCIPclockFree(&(*eventhdlr)->eventtime);
   SCIPclockFree(&(*eventhdlr)->setuptime);

   BMSfreeMemoryArray(&(*eventhdlr)->name);
   BMSfreeMemoryArray(&(*eventhdlr)->desc);
   BMSfreeMemory(eventhdlr);

   return SCIP_OKAY;
}

/** initializes event handler */
SCIP_RETCODE SCIPeventhdlrInit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   if( eventhdlr->initialized )
   {
      SCIPerrorMessage("event handler <%s> already initialized\n", eventhdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(eventhdlr->setuptime);
      SCIPclockReset(eventhdlr->eventtime);
   }

   if( eventhdlr->eventinit != NULL )
   {
      /* start timing */
      SCIPclockStart(eventhdlr->setuptime, set);

      SCIP_CALL( eventhdlr->eventinit(set->scip, eventhdlr) );

      /* stop timing */
      SCIPclockStop(eventhdlr->setuptime, set);
   }
   eventhdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of event handler */
SCIP_RETCODE SCIPeventhdlrExit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   if( !eventhdlr->initialized )
   {
      SCIPerrorMessage("event handler <%s> not initialized\n", eventhdlr->name);
      return SCIP_INVALIDCALL;
   }

   if( eventhdlr->eventexit != NULL )
   {
      /* start timing */
      SCIPclockStart(eventhdlr->setuptime, set);

      SCIP_CALL( eventhdlr->eventexit(set->scip, eventhdlr) );

      /* stop timing */
      SCIPclockStop(eventhdlr->setuptime, set);
   }
   eventhdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs event handler that the branch and bound process is being started */
SCIP_RETCODE SCIPeventhdlrInitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   /* call solving process initialization method of event handler */
   if( eventhdlr->eventinitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(eventhdlr->setuptime, set);

      SCIP_CALL( eventhdlr->eventinitsol(set->scip, eventhdlr) );

      /* stop timing */
      SCIPclockStop(eventhdlr->setuptime, set);
   }

   return SCIP_OKAY;
}

/** informs event handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPeventhdlrExitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(eventhdlr != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of event handler */
   if( eventhdlr->eventexitsol != NULL )
   {
      /* start timing */
      SCIPclockStart(eventhdlr->setuptime, set);

      SCIP_CALL( eventhdlr->eventexitsol(set->scip, eventhdlr) );

      /* stop timing */
      SCIPclockStop(eventhdlr->setuptime, set);
   }

   return SCIP_OKAY;
}

/** calls execution method of event handler */
SCIP_RETCODE SCIPeventhdlrExec(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENT*           event,              /**< event to call event handler with */
   SCIP_EVENTDATA*       eventdata           /**< user data for the issued event */
   )
{
   assert(eventhdlr != NULL);
   assert(eventhdlr->eventexec != NULL);
   assert(set != NULL);
   assert(event != NULL);

   SCIPsetDebugMsg(set, "execute event of handler <%s> with event %p of type 0x%" SCIP_EVENTTYPE_FORMAT "\n", eventhdlr->name, (void*)event, event->eventtype);

#ifdef TIMEEVENTEXEC
   /* start timing */
   SCIPclockStart(eventhdlr->eventtime, set);
#endif

   SCIP_CALL( eventhdlr->eventexec(set->scip, eventhdlr, event, eventdata) );

#ifdef TIMEEVENTEXEC
   /* stop timing */
   SCIPclockStop(eventhdlr->eventtime, set);
#endif

   return SCIP_OKAY;
}

/** gets name of event handler */
const char* SCIPeventhdlrGetName(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->name;
}

/** gets user data of event handler */
SCIP_EVENTHDLRDATA* SCIPeventhdlrGetData(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->eventhdlrdata;
}

/** sets user data of event handler; user has to free old data in advance! */
void SCIPeventhdlrSetData(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< new event handler user data */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventhdlrdata = eventhdlrdata;
}

/** sets copy callback for all events of this event handler */
void SCIPeventhdlrSetCopy(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy))      /**< copy callback for events */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventcopy = eventcopy;
}

/** sets destructor callback of this event handler */
void SCIPeventhdlrSetFree(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTFREE   ((*eventfree))      /**< destructor callback of event handler */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventfree = eventfree;
}

/** sets initialization callback of this event handler */
void SCIPeventhdlrSetInit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit))      /**< initialization callback of event handler */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventinit = eventinit;
}

/** sets deinitialization callback of this event handler */
void SCIPeventhdlrSetExit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit))      /**< deinitialization callback of event handler */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventexit = eventexit;
}

/** sets solving process initialization callback of this event handler */
void SCIPeventhdlrSetInitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol))   /**< solving process initialization callback of event handler */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventinitsol = eventinitsol;
}

/** sets solving process deinitialization callback of this event handler */
void SCIPeventhdlrSetExitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol))   /**< solving process deinitialization callback of event handler */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventexitsol = eventexitsol;
}

/** sets callback to free specific event data */
void SCIPeventhdlrSetDelete(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete))    /**< callback to free specific event data */
   )
{
   assert(eventhdlr != NULL);

   eventhdlr->eventdelete = eventdelete;
}

/** is event handler initialized? */
SCIP_Bool SCIPeventhdlrIsInitialized(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return eventhdlr->initialized;
}

/** enables or disables all clocks of \p eventhdlr, depending on the value of the flag */
void SCIPeventhdlrEnableOrDisableClocks(
   SCIP_EVENTHDLR*       eventhdlr,          /**< the event handler for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the event handler be enabled? */
   )
{
   assert(eventhdlr != NULL);

   SCIPclockEnableOrDisable(eventhdlr->setuptime, enable);
   SCIPclockEnableOrDisable(eventhdlr->eventtime, enable);
}

/** gets time in seconds used in this event handler for setting up for next stages */
SCIP_Real SCIPeventhdlrGetSetupTime(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return SCIPclockGetTime(eventhdlr->setuptime);
}

/** gets time in seconds used in this event handler, this measurement is currently disabled so this method will return
 *  0, define TIMEEVENTEXEC in the beginning of this file to enable
 */
SCIP_Real SCIPeventhdlrGetTime(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(eventhdlr != NULL);

   return SCIPclockGetTime(eventhdlr->eventtime);
}



/*
 * Event methods
 */


/** creates a synchronization event */
SCIP_RETCODE SCIPeventCreateSync(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(event != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_SYNC;

   return SCIP_OKAY;
}

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPeventGetType
#undef SCIPeventGetOldobj
#undef SCIPeventGetNewobj
#undef SCIPeventGetOldbound
#undef SCIPeventGetNewbound
#undef SCIPeventGetNode
#undef SCIPeventGetSol
#undef SCIPeventGetRowCol
#undef SCIPeventGetRowOldCoefVal
#undef SCIPeventGetRowNewCoefVal
#undef SCIPeventGetRowOldConstVal
#undef SCIPeventGetRowNewConstVal
#undef SCIPeventGetRowSide
#undef SCIPeventGetRowOldSideVal
#undef SCIPeventGetRowNewSideVal

/** creates an event for an addition of a variable to the problem */
SCIP_RETCODE SCIPeventCreateVarAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that was added to the problem */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_VARADDED;
   (*event)->data.eventvaradded.var = var;

   return SCIP_OKAY;
}

/** creates an event for a deletion of a variable from the problem */
SCIP_RETCODE SCIPeventCreateVarDeleted(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that is to be deleted from the problem */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_VARDELETED;
   (*event)->data.eventvardeleted.var = var;

   return SCIP_OKAY;
}

/** creates an event for a fixing of a variable */
SCIP_RETCODE SCIPeventCreateVarFixed(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that was fixed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_VARFIXED;
   (*event)->data.eventvarfixed.var = var;

   return SCIP_OKAY;
}

/** creates an event for a change in the number of locks of a variable down to zero or one */
SCIP_RETCODE SCIPeventCreateVarUnlocked(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that changed the number of locks */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_VARUNLOCKED;
   (*event)->data.eventvarunlocked.var = var;

   return SCIP_OKAY;
}

/** creates an event for a change in the objective value of a variable */
SCIP_RETCODE SCIPeventCreateObjChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose objective value changed */
   SCIP_Real             oldobj,             /**< old objective value before value changed */
   SCIP_Real             newobj              /**< new objective value after value changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldobj != newobj); /*lint !e777*/

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_OBJCHANGED;
   (*event)->data.eventobjchg.var = var;
   (*event)->data.eventobjchg.oldobj = oldobj;
   (*event)->data.eventobjchg.newobj = newobj;

   return SCIP_OKAY;
}

/** creates an event for a change in the global lower bound of a variable */
SCIP_RETCODE SCIPeventCreateGlbChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldbound != newbound); /*lint !e777*/

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_GLBCHANGED;
   (*event)->data.eventbdchg.var = var;
   (*event)->data.eventbdchg.oldbound = oldbound;
   (*event)->data.eventbdchg.newbound = newbound;

   return SCIP_OKAY;
}

/** creates an event for a change in the global upper bound of a variable */
SCIP_RETCODE SCIPeventCreateGubChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldbound != newbound); /*lint !e777*/

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_GUBCHANGED;
   (*event)->data.eventbdchg.var = var;
   (*event)->data.eventbdchg.oldbound = oldbound;
   (*event)->data.eventbdchg.newbound = newbound;

   return SCIP_OKAY;
}

/** creates an event for a change in the lower bound of a variable */
SCIP_RETCODE SCIPeventCreateLbChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldbound != newbound); /*lint !e777*/

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   if( newbound > oldbound )
      (*event)->eventtype = SCIP_EVENTTYPE_LBTIGHTENED;
   else
      (*event)->eventtype = SCIP_EVENTTYPE_LBRELAXED;
   (*event)->data.eventbdchg.var = var;
   (*event)->data.eventbdchg.oldbound = oldbound;
   (*event)->data.eventbdchg.newbound = newbound;

   return SCIP_OKAY;
}

/** creates an event for a change in the upper bound of a variable */
SCIP_RETCODE SCIPeventCreateUbChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(oldbound != newbound); /*lint !e777*/

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   if( newbound < oldbound )
      (*event)->eventtype = SCIP_EVENTTYPE_UBTIGHTENED;
   else
      (*event)->eventtype = SCIP_EVENTTYPE_UBRELAXED;
   (*event)->data.eventbdchg.var = var;
   (*event)->data.eventbdchg.oldbound = oldbound;
   (*event)->data.eventbdchg.newbound = newbound;

   return SCIP_OKAY;
}

/** creates an event for an addition of a domain hole to a variable */
SCIP_RETCODE SCIPeventCreateGholeAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_GHOLEADDED;
   (*event)->data.eventhole.var = var;
   (*event)->data.eventhole.left = left;
   (*event)->data.eventhole.right = right;

   return SCIP_OKAY;
}

/** creates an event for removing a domain hole of a variable */
SCIP_RETCODE SCIPeventCreateGholeRemoved(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in hole */
   SCIP_Real             right               /**< right bound of open interval in hole */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_GHOLEREMOVED;
   (*event)->data.eventhole.var = var;
   (*event)->data.eventhole.left = left;
   (*event)->data.eventhole.right = right;

   return SCIP_OKAY;
}

/** creates an event for an addition of a domain hole to a variable */
SCIP_RETCODE SCIPeventCreateLholeAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_LHOLEADDED;
   (*event)->data.eventhole.var = var;
   (*event)->data.eventhole.left = left;
   (*event)->data.eventhole.right = right;

   return SCIP_OKAY;
}

/** creates an event for removing a domain hole of a variable */
SCIP_RETCODE SCIPeventCreateLholeRemoved(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in hole */
   SCIP_Real             right               /**< right bound of open interval in hole */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_LHOLEREMOVED;
   (*event)->data.eventhole.var = var;
   (*event)->data.eventhole.left = left;
   (*event)->data.eventhole.right = right;

   return SCIP_OKAY;
}

/** creates an event for an addition to the variable's implications list, clique or variable bounds information */
SCIP_RETCODE SCIPeventCreateImplAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that was fixed */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_IMPLADDED;
   (*event)->data.eventimpladd.var = var;

   return SCIP_OKAY;
}

/** creates an event for the addition of a linear row to the separation storage */
SCIP_RETCODE SCIPeventCreateRowAddedSepa(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was added to the separation storage*/
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_ROWADDEDSEPA;
   (*event)->data.eventrowaddedsepa.row = row;

   return SCIP_OKAY;
}

/** creates an event for the deletion of a linear row from the separation storage */
SCIP_RETCODE SCIPeventCreateRowDeletedSepa(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was deleted from the separation storage */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_ROWDELETEDSEPA;
   (*event)->data.eventrowdeletedsepa.row = row;

   return SCIP_OKAY;
}

/** creates an event for the addition of a linear row to the LP */
SCIP_RETCODE SCIPeventCreateRowAddedLP(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was added to the LP */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_ROWADDEDLP;
   (*event)->data.eventrowaddedlp.row = row;

   return SCIP_OKAY;
}

/** creates an event for the deletion of a linear row from the LP */
SCIP_RETCODE SCIPeventCreateRowDeletedLP(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was deleted from the LP */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_ROWDELETEDLP;
   (*event)->data.eventrowdeletedlp.row = row;

   return SCIP_OKAY;
}

/** creates an event for the change of a coefficient in a linear row */
SCIP_RETCODE SCIPeventCreateRowCoefChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row in which a coefficient changed */
   SCIP_COL*             col,                /**< column which coefficient changed */
   SCIP_Real             oldval,             /**< old value of coefficient */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_ROWCOEFCHANGED;
   (*event)->data.eventrowcoefchanged.row = row;
   (*event)->data.eventrowcoefchanged.col = col;
   (*event)->data.eventrowcoefchanged.oldval = oldval;
   (*event)->data.eventrowcoefchanged.newval = newval;

   return SCIP_OKAY;
}

/** creates an event for the change of a constant in a linear row */
SCIP_RETCODE SCIPeventCreateRowConstChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row in which the constant changed */
   SCIP_Real             oldval,             /**< old value of constant */
   SCIP_Real             newval              /**< new value of constant */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_ROWCONSTCHANGED;
   (*event)->data.eventrowconstchanged.row = row;
   (*event)->data.eventrowconstchanged.oldval = oldval;
   (*event)->data.eventrowconstchanged.newval = newval;

   return SCIP_OKAY;
}

/** creates an event for the change of a side of a linear row */
SCIP_RETCODE SCIPeventCreateRowSideChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row which side has changed */
   SCIP_SIDETYPE         side,               /**< which side has changed */
   SCIP_Real             oldval,             /**< old value of side */
   SCIP_Real             newval              /**< new value of side */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);
   assert(row != NULL);

   /* create event data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, event) );
   (*event)->eventtype = SCIP_EVENTTYPE_ROWSIDECHANGED;
   (*event)->data.eventrowsidechanged.row = row;
   (*event)->data.eventrowsidechanged.side = side;
   (*event)->data.eventrowsidechanged.oldval = oldval;
   (*event)->data.eventrowsidechanged.newval = newval;

   return SCIP_OKAY;
}

/** frees an event */
SCIP_RETCODE SCIPeventFree(
   SCIP_EVENT**          event,              /**< event to free */
   BMS_BLKMEM*           blkmem              /**< block memory buffer */
   )
{
   assert(event != NULL);
   assert(blkmem != NULL);

   BMSfreeBlockMemory(blkmem, event);

   return SCIP_OKAY;
}

/** disables an event */
static
void eventDisable(
   SCIP_EVENT*           event               /**< event to disable */
   )
{
   assert(event != NULL);

   event->eventtype = SCIP_EVENTTYPE_DISABLED;
}

/** gets type of event */
SCIP_EVENTTYPE SCIPeventGetType(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   return event->eventtype;
}

/** sets type of event */
SCIP_RETCODE SCIPeventChgType(
   SCIP_EVENT*           event,              /**< event */
   SCIP_EVENTTYPE        eventtype           /**< new event type */
   )
{
   assert(event != NULL);

   event->eventtype = eventtype;

   return SCIP_OKAY;
}

/** gets variable for a variable event (var added, var deleted, var fixed, objective value or domain change) */
SCIP_VAR* SCIPeventGetVar(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {  
   case SCIP_EVENTTYPE_VARADDED:
      assert(event->data.eventvaradded.var != NULL);
      return event->data.eventvaradded.var;

   case SCIP_EVENTTYPE_VARDELETED:
      assert(event->data.eventvardeleted.var != NULL);
      return event->data.eventvardeleted.var;

   case SCIP_EVENTTYPE_VARFIXED:
      assert(event->data.eventvarfixed.var != NULL);
      return event->data.eventvarfixed.var;

   case SCIP_EVENTTYPE_VARUNLOCKED:
      assert(event->data.eventvarunlocked.var != NULL);
      return event->data.eventvarunlocked.var;

   case SCIP_EVENTTYPE_OBJCHANGED:
      assert(event->data.eventobjchg.var != NULL);
      return event->data.eventobjchg.var;

   case SCIP_EVENTTYPE_GLBCHANGED:
   case SCIP_EVENTTYPE_GUBCHANGED:
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      assert(event->data.eventbdchg.var != NULL);
      return event->data.eventbdchg.var;

   case SCIP_EVENTTYPE_GHOLEADDED:
   case SCIP_EVENTTYPE_GHOLEREMOVED:
   case SCIP_EVENTTYPE_LHOLEADDED:
   case SCIP_EVENTTYPE_LHOLEREMOVED:
      assert(event->data.eventhole.var != NULL);
      return event->data.eventhole.var;

   case SCIP_EVENTTYPE_IMPLADDED:
      assert(event->data.eventimpladd.var != NULL);
      return event->data.eventimpladd.var;

   default:
      SCIPerrorMessage("event does not belong to a variable\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** sets variable for a variable event */
SCIP_RETCODE SCIPeventChgVar(
   SCIP_EVENT*           event,              /**< event */
   SCIP_VAR*             var                 /**< new variable */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {  
   case SCIP_EVENTTYPE_VARADDED:
      assert(event->data.eventvaradded.var != NULL);
      event->data.eventvaradded.var = var;
      break;

   case SCIP_EVENTTYPE_VARDELETED:
      assert(event->data.eventvardeleted.var != NULL);
      event->data.eventvardeleted.var = var;
      break;

   case SCIP_EVENTTYPE_VARFIXED:
      assert(event->data.eventvarfixed.var != NULL);
      event->data.eventvarfixed.var = var;
      break;

   case SCIP_EVENTTYPE_VARUNLOCKED:
      assert(event->data.eventvarunlocked.var != NULL);
      event->data.eventvarunlocked.var = var;
      break;

   case SCIP_EVENTTYPE_OBJCHANGED:
      assert(event->data.eventobjchg.var != NULL);
      event->data.eventobjchg.var = var;
      break;

   case SCIP_EVENTTYPE_GLBCHANGED:
   case SCIP_EVENTTYPE_GUBCHANGED:
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      assert(event->data.eventbdchg.var != NULL);
      event->data.eventbdchg.var = var;
      break;

   case SCIP_EVENTTYPE_GHOLEADDED:
   case SCIP_EVENTTYPE_GHOLEREMOVED:
   case SCIP_EVENTTYPE_LHOLEADDED:
   case SCIP_EVENTTYPE_LHOLEREMOVED:
      assert(event->data.eventhole.var != NULL);
      event->data.eventhole.var = var;
      break;

   case SCIP_EVENTTYPE_IMPLADDED:
      assert(event->data.eventimpladd.var != NULL);
      event->data.eventimpladd.var = var;
      break;

   default:
      SCIPerrorMessage("event does not belong to a variable\n");
      return SCIP_INVALIDDATA;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets old objective value for an objective value change event */
SCIP_Real SCIPeventGetOldobj(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( event->eventtype != SCIP_EVENTTYPE_OBJCHANGED )
   {
      SCIPerrorMessage("event is not an objective value change event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventobjchg.oldobj;
}

/** gets new objective value for an objective value change event */
SCIP_Real SCIPeventGetNewobj(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( event->eventtype != SCIP_EVENTTYPE_OBJCHANGED )
   {
      SCIPerrorMessage("event is not an objective value change event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventobjchg.newobj;
}

/** gets old bound for a bound change event */
SCIP_Real SCIPeventGetOldbound(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {  
   case SCIP_EVENTTYPE_GLBCHANGED:
   case SCIP_EVENTTYPE_GUBCHANGED:
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      return event->data.eventbdchg.oldbound;

   default:
      SCIPerrorMessage("event is not a bound change event\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets new bound for a bound change event */
SCIP_Real SCIPeventGetNewbound(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {  
   case SCIP_EVENTTYPE_GLBCHANGED:
   case SCIP_EVENTTYPE_GUBCHANGED:
   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      return event->data.eventbdchg.newbound;

   default:
      SCIPerrorMessage("event is not a bound change event\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets node for a node or LP event */
SCIP_NODE* SCIPeventGetNode(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & (SCIP_EVENTTYPE_NODEEVENT | SCIP_EVENTTYPE_LPEVENT)) == 0 )
   {
      SCIPerrorMessage("event is neither node nor LP event\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   return event->data.node;
}

/** sets node for a node or LP event */
SCIP_RETCODE SCIPeventChgNode(
   SCIP_EVENT*           event,              /**< event */
   SCIP_NODE*            node                /**< new node */
   )
{
   assert(event != NULL);

   if( (event->eventtype & (SCIP_EVENTTYPE_NODEEVENT | SCIP_EVENTTYPE_LPEVENT)) == 0 )
   {
      SCIPerrorMessage("event is neither node nor LP event\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   event->data.node = node;

   return SCIP_OKAY;
}

/** gets solution for a primal solution event */
SCIP_SOL* SCIPeventGetSol(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_SOLEVENT) == 0 )
   {
      SCIPerrorMessage("event is not a primal solution event\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   return event->data.sol;
}

/** sets solution for a primal solution event */
SCIP_RETCODE SCIPeventChgSol(
   SCIP_EVENT*           event,              /**< event */
   SCIP_SOL*             sol                 /**< new primal solution */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_SOLEVENT) == 0 )
   {
      SCIPerrorMessage("event is not a primal solution event\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;  /*lint !e527*/
   }

   event->data.sol = sol;

   return SCIP_OKAY;
}

/** gets the left bound of open interval in the hole */
SCIP_Real SCIPeventGetHoleLeft(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_HOLECHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a hole added or removed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventhole.left;
}

/** gets the right bound of open interval in the hole */
SCIP_Real SCIPeventGetHoleRight(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_HOLECHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a hole added or removed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventhole.right;
}

/** gets row for a row event */
SCIP_ROW* SCIPeventGetRow(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   switch( event->eventtype )
   {
      case SCIP_EVENTTYPE_ROWADDEDSEPA:
         return event->data.eventrowaddedsepa.row;
      case SCIP_EVENTTYPE_ROWDELETEDSEPA:
         return event->data.eventrowdeletedsepa.row;
      case SCIP_EVENTTYPE_ROWADDEDLP:
         return event->data.eventrowaddedlp.row;
      case SCIP_EVENTTYPE_ROWDELETEDLP:
         return event->data.eventrowdeletedlp.row;
      case SCIP_EVENTTYPE_ROWCOEFCHANGED:
         return event->data.eventrowcoefchanged.row;
      case SCIP_EVENTTYPE_ROWCONSTCHANGED:
         return event->data.eventrowconstchanged.row;
      case SCIP_EVENTTYPE_ROWSIDECHANGED:
         return event->data.eventrowsidechanged.row;
      default:
         SCIPerrorMessage("event does not belong to a row\n");
         SCIPABORT();
         return NULL;  /*lint !e527*/
   }
}

/** gets column for a row change coefficient event */
SCIP_COL* SCIPeventGetRowCol(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWCOEFCHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row coefficient changed event\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   return event->data.eventrowcoefchanged.col;
}

/** gets old coefficient value for a row change coefficient event */
SCIP_Real SCIPeventGetRowOldCoefVal(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWCOEFCHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row coefficient changed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventrowcoefchanged.oldval;
}

/** gets new coefficient value for a row change coefficient event */
SCIP_Real SCIPeventGetRowNewCoefVal(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWCOEFCHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row coefficient changed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventrowcoefchanged.newval;
}

/** gets old constant value for a row change constant event */
SCIP_Real SCIPeventGetRowOldConstVal(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWCONSTCHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row coefficient changed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventrowconstchanged.oldval;
}

/** gets new constant value for a row change constant event */
SCIP_Real SCIPeventGetRowNewConstVal(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWCONSTCHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row coefficient changed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventrowconstchanged.newval;
}

/** gets side for a row change side event */
SCIP_SIDETYPE SCIPeventGetRowSide(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWSIDECHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row side changed event\n");
      SCIPABORT();
      return SCIP_SIDETYPE_LEFT;  /*lint !e527*/
   }

   return event->data.eventrowsidechanged.side;
}

/** gets old side value for a row change side event */
SCIP_Real SCIPeventGetRowOldSideVal(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWSIDECHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row side changed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventrowsidechanged.oldval;
}

/** gets new side value for a row change side event */
SCIP_Real SCIPeventGetRowNewSideVal(
   SCIP_EVENT*           event               /**< event */
   )
{
   assert(event != NULL);

   if( (event->eventtype & SCIP_EVENTTYPE_ROWSIDECHANGED) == 0 )
   {
      SCIPerrorMessage("event is not a row side changed event\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   return event->data.eventrowsidechanged.newval;
}

/** processes event by calling the appropriate event handlers */
SCIP_RETCODE SCIPeventProcess(
   SCIP_EVENT*           event,              /**< event */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data; only needed for objchanged events, or NULL */
   SCIP_LP*              lp,                 /**< current LP data; only needed for obj/boundchanged events, or NULL */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for bound change events, or NULL */
   SCIP_EVENTFILTER*     eventfilter         /**< event filter for global events; not needed for variable specific events */
   )
{
   SCIP_VAR* var;

   assert(event != NULL);
   assert((event->eventtype & SCIP_EVENTTYPE_OBJCHANGED) == 0 || primal != NULL);
   assert((event->eventtype & (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED)) == 0 || lp != NULL);
   assert((event->eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) == 0 || branchcand != NULL);

   SCIPsetDebugMsg(set, "processing event of type 0x%" SCIP_EVENTTYPE_FORMAT "\n", event->eventtype);

   switch( event->eventtype )
   {
   case SCIP_EVENTTYPE_DISABLED:
      break;

   case SCIP_EVENTTYPE_SYNC: /*lint !e30 !e142*/
   case SCIP_EVENTTYPE_VARADDED:
   case SCIP_EVENTTYPE_PRESOLVEROUND:
   case SCIP_EVENTTYPE_NODEFOCUSED:
   case SCIP_EVENTTYPE_NODEFEASIBLE:
   case SCIP_EVENTTYPE_NODEINFEASIBLE:
   case SCIP_EVENTTYPE_NODEBRANCHED:
   case SCIP_EVENTTYPE_FIRSTLPSOLVED:
   case SCIP_EVENTTYPE_LPSOLVED:
   case SCIP_EVENTTYPE_POORSOLFOUND:
   case SCIP_EVENTTYPE_BESTSOLFOUND:
   case SCIP_EVENTTYPE_ROWADDEDSEPA:
   case SCIP_EVENTTYPE_ROWDELETEDSEPA:
   case SCIP_EVENTTYPE_ROWADDEDLP:
   case SCIP_EVENTTYPE_ROWDELETEDLP:
   case SCIP_EVENTTYPE_ROWCOEFCHANGED:
   case SCIP_EVENTTYPE_ROWCONSTCHANGED:
   case SCIP_EVENTTYPE_ROWSIDECHANGED:
      SCIP_CALL( SCIPeventfilterProcess(eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_VARDELETED:
      var = event->data.eventvardeleted.var;
      assert(var != NULL);

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_VARFIXED:
      var = event->data.eventvarfixed.var;
      assert(var != NULL);

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_VARUNLOCKED:
      var = event->data.eventvarunlocked.var;
      assert(var != NULL);

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_OBJCHANGED:
      var = event->data.eventobjchg.var;
      assert(var != NULL);
      assert(var->eventqueueindexobj == -1);
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE);

      /* inform LP about the objective change */
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPcolChgObj(SCIPvarGetCol(var), set, lp, event->data.eventobjchg.newobj) );
         }
         SCIP_CALL( SCIPlpUpdateVarObj(lp, set, var, event->data.eventobjchg.oldobj, event->data.eventobjchg.newobj) );
      }

      /* inform all existing primal solutions about the objective change (only if this is not a temporary change in
       * probing mode)
       */
      if( ! lp->divingobjchg )
      {
         SCIPprimalUpdateVarObj(primal, var, event->data.eventobjchg.oldobj, event->data.eventobjchg.newobj);
      }

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_GLBCHANGED:
      var = event->data.eventbdchg.var;
      assert(var != NULL);

      /* inform LP about global bound change */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         assert(SCIPvarGetProbindex(var) >= 0);
         SCIP_CALL( SCIPlpUpdateVarLbGlobal(lp, set, var, event->data.eventbdchg.oldbound,
               event->data.eventbdchg.newbound) );
      }

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_GUBCHANGED:
      var = event->data.eventbdchg.var;
      assert(var != NULL);

      /* inform LP about global bound change */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         assert(SCIPvarGetProbindex(var) >= 0);
         SCIP_CALL( SCIPlpUpdateVarUbGlobal(lp, set, var, event->data.eventbdchg.oldbound,
               event->data.eventbdchg.newbound) );
      }

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_LBTIGHTENED:
   case SCIP_EVENTTYPE_LBRELAXED:
      var = event->data.eventbdchg.var;
      assert(var != NULL);
      assert(var->eventqueueindexlb == -1);

      /* inform LP about bound change and update branching candidates */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         assert(SCIPvarGetProbindex(var) >= 0);
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPcolChgLb(SCIPvarGetCol(var), set, lp, event->data.eventbdchg.newbound) );
         }
         SCIP_CALL( SCIPlpUpdateVarLb(lp, set, var, event->data.eventbdchg.oldbound,
               event->data.eventbdchg.newbound) );
         SCIP_CALL( SCIPbranchcandUpdateVar(branchcand, set, var) );
      }

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_UBTIGHTENED:
   case SCIP_EVENTTYPE_UBRELAXED:
      var = event->data.eventbdchg.var;
      assert(var != NULL);
      assert(var->eventqueueindexub == -1);

      /* inform LP about bound change and update branching candidates */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         assert(SCIPvarGetProbindex(var) >= 0);
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPcolChgUb(SCIPvarGetCol(var), set, lp, event->data.eventbdchg.newbound) );
         }
         SCIP_CALL( SCIPlpUpdateVarUb(lp, set, var, event->data.eventbdchg.oldbound, 
               event->data.eventbdchg.newbound) );
         SCIP_CALL( SCIPbranchcandUpdateVar(branchcand, set, var) );
      }

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_GHOLEADDED:
   case SCIP_EVENTTYPE_GHOLEREMOVED:
   case SCIP_EVENTTYPE_LHOLEADDED:
   case SCIP_EVENTTYPE_LHOLEREMOVED:
      var = event->data.eventhole.var;
      assert(var != NULL);

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   case SCIP_EVENTTYPE_IMPLADDED:
      var = event->data.eventimpladd.var;
      assert(var != NULL);
      assert(!var->eventqueueimpl);

      /* process variable's event filter */
      SCIP_CALL( SCIPeventfilterProcess(var->eventfilter, set, event) );
      break;

   default:
      SCIPerrorMessage("unknown event type <%d>\n", event->eventtype);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}



/*
 * Event filter methods
 */

/** resizes eventfilter arrays to be able to store at least num entries */
static
SCIP_RETCODE eventfilterEnsureMem(
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of node slots in array */
   )
{
   assert(eventfilter != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);

   if( num > eventfilter->size )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &eventfilter->eventtypes, eventfilter->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &eventfilter->eventhdlrs, eventfilter->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &eventfilter->eventdata, eventfilter->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &eventfilter->nextpos, eventfilter->size, newsize) );
      eventfilter->size = newsize;
   }
   assert(num <= eventfilter->size);

   return SCIP_OKAY;
}

/** creates an event filter */
SCIP_RETCODE SCIPeventfilterCreate(
   SCIP_EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   BMS_BLKMEM*           blkmem              /**< block memory buffer */
   )
{
   assert(eventfilter != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, eventfilter) );
   (*eventfilter)->eventtypes = NULL;
   (*eventfilter)->eventhdlrs = NULL;
   (*eventfilter)->eventdata = NULL;
   (*eventfilter)->nextpos = NULL;
   (*eventfilter)->size = 0;
   (*eventfilter)->len = 0;
   (*eventfilter)->firstfreepos = -1;
   (*eventfilter)->firstdeletedpos = -1;
   (*eventfilter)->eventmask = SCIP_EVENTTYPE_DISABLED;
   (*eventfilter)->delayedeventmask = SCIP_EVENTTYPE_DISABLED;
   (*eventfilter)->delayupdates = FALSE;

   return SCIP_OKAY;
}

/** frees an event filter and the associated event data entries */
SCIP_RETCODE SCIPeventfilterFree(
   SCIP_EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;

   assert(eventfilter != NULL);
   assert(*eventfilter != NULL);
   assert(!(*eventfilter)->delayupdates);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   /* free event data */
   for( i = 0; i < (*eventfilter)->len; ++i )
   {
      if( (*eventfilter)->eventtypes[i] != SCIP_EVENTTYPE_DISABLED )
      {
         assert((*eventfilter)->eventhdlrs[i] != NULL);
         if( (*eventfilter)->eventhdlrs[i]->eventdelete != NULL )
         {
            SCIP_CALL( (*eventfilter)->eventhdlrs[i]->eventdelete(set->scip, (*eventfilter)->eventhdlrs[i],
                  &(*eventfilter)->eventdata[i]) );
         }
      }
   }

   /* free event filter data */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*eventfilter)->eventtypes, (*eventfilter)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*eventfilter)->eventhdlrs, (*eventfilter)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*eventfilter)->eventdata, (*eventfilter)->size);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*eventfilter)->nextpos, (*eventfilter)->size);
   BMSfreeBlockMemory(blkmem, eventfilter);

   return SCIP_OKAY;
}

/** adds element to event filter */
SCIP_RETCODE SCIPeventfilterAdd(
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   int pos;

   assert(eventfilter != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(eventhdlr != NULL);

   if( eventfilter->delayupdates )
   {
      /* insert addition to the end of the arrays;
       * in delayed addition we have to add to the end of the arrays, in order to not destroy the validity of the
       * arrays we are currently iterating over
       */
      SCIP_CALL( eventfilterEnsureMem(eventfilter, blkmem, set, eventfilter->len + 1) );
      pos = eventfilter->len;
      eventfilter->len++;

      /* update delayed event filter mask */
      eventfilter->delayedeventmask |= eventtype;
   }
   else
   {
      if( eventfilter->firstfreepos == -1 )
      {
         /* insert addition to the end of the arrays */
         SCIP_CALL( eventfilterEnsureMem(eventfilter, blkmem, set, eventfilter->len + 1) );
         pos = eventfilter->len;
         eventfilter->len++;
      }
      else
      {
         /* use the first free slot to store the added event filter entry */
         pos = eventfilter->firstfreepos;
         assert(0 <= pos && pos < eventfilter->len);
         assert(eventfilter->eventtypes[pos] == SCIP_EVENTTYPE_DISABLED);
         eventfilter->firstfreepos = eventfilter->nextpos[pos];
         assert(-1 <= eventfilter->firstfreepos && eventfilter->firstfreepos < eventfilter->len);
      }

      /* update event filter mask */
      eventfilter->eventmask |= eventtype;
   }
   assert(0 <= pos && pos < eventfilter->len);

   eventfilter->eventtypes[pos] = eventtype;
   eventfilter->eventhdlrs[pos] = eventhdlr;
   eventfilter->eventdata[pos] = eventdata;
   eventfilter->nextpos[pos] = -2;

   if( filterpos != NULL )
      *filterpos = pos;

   return SCIP_OKAY;
}

/** linear search for the given entry in event filter */
static
int eventfilterSearch(
   SCIP_EVENTFILTER*const eventfilter,       /**< event filter */
   SCIP_EVENTTYPE const  eventtype,          /**< event type */
   SCIP_EVENTHDLR*const  eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*const  eventdata           /**< event data to pass to the event handler for the event processing */
   )
{
   int i;

   assert(eventfilter != NULL);
   assert(eventtype != SCIP_EVENTTYPE_DISABLED);
   assert(eventhdlr != NULL);

   for( i = eventfilter->len - 1; i >= 0; --i )
   {
      if( eventdata == eventfilter->eventdata[i]
         && eventhdlr == eventfilter->eventhdlrs[i]
         && eventtype == eventfilter->eventtypes[i]
         && eventfilter->nextpos[i] == -2 )
         return i;
   }

   return -1;
}

/** deletes element from event filter */
SCIP_RETCODE SCIPeventfilterDel(
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int                   filterpos           /**< position of event filter entry, or -1 if unknown */
   )
{
   assert(eventfilter != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(eventtype != SCIP_EVENTTYPE_DISABLED);
   assert(eventhdlr != NULL);
   assert(-1 <= filterpos && filterpos < eventfilter->len);

   /* search position of event filter entry, if not given by the user */
   if( filterpos == -1 )
      filterpos = eventfilterSearch(eventfilter, eventtype, eventhdlr, eventdata);
   if( filterpos == -1 )
   {
      SCIPerrorMessage("no event for event handler %p with data %p and event mask 0x%x found in event filter %p\n",
         eventhdlr, eventdata, eventtype, eventfilter);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= filterpos && filterpos < eventfilter->len);
   assert(eventfilter->eventtypes[filterpos] == eventtype);
   assert(eventfilter->eventhdlrs[filterpos] == eventhdlr);
   assert(eventfilter->eventdata[filterpos] == eventdata);
   assert(eventfilter->nextpos[filterpos] == -2);

   /* if updates are delayed, insert entry into the list of delayed deletions;
    * otherwise, delete the entry from the filter directly and add the slot to the free list
    */
   if( eventfilter->delayupdates )
   {
      /* append filterpos to the list of deleted entries */
      eventfilter->nextpos[filterpos] = eventfilter->firstdeletedpos;
      eventfilter->firstdeletedpos = filterpos;
   }
   else
   {
      /* disable the entry in the filter and add the slot to the free list */
      assert(eventfilter->nextpos[filterpos] == -2);
      eventfilter->eventtypes[filterpos] = SCIP_EVENTTYPE_DISABLED;
      eventfilter->nextpos[filterpos] = eventfilter->firstfreepos;
      eventfilter->firstfreepos = filterpos;
   }

   return SCIP_OKAY;
}

/** makes the event filter to delay and buffer all updates until eventfilterProcessUpdates() is called */
static
void eventfilterDelayUpdates(
   SCIP_EVENTFILTER*     eventfilter         /**< event filter */
   )
{
   assert(eventfilter != NULL);
   assert(!eventfilter->delayupdates);
   assert(eventfilter->delayedeventmask == SCIP_EVENTTYPE_DISABLED);

   eventfilter->delayupdates = TRUE;
}

/** processes all delayed additions and deletions */
static
void eventfilterProcessUpdates(
   SCIP_EVENTFILTER*     eventfilter         /**< event filter */
   )
{
   int pos;
   int nextpos;

   assert(eventfilter != NULL);
   assert(eventfilter->delayupdates);

   /* move deleted entries into the free list and disable them */
   pos = eventfilter->firstdeletedpos;
   while( pos != -1 )
   {
      assert(0 <= pos && pos < eventfilter->len);
      assert(eventfilter->eventtypes[pos] != SCIP_EVENTTYPE_DISABLED);
      assert(eventfilter->nextpos[pos] >= -1);

      nextpos = eventfilter->nextpos[pos];
      eventfilter->nextpos[pos] = eventfilter->firstfreepos;
      eventfilter->firstfreepos = pos;
      eventfilter->eventtypes[pos] = SCIP_EVENTTYPE_DISABLED;
      pos = nextpos;
   }
   eventfilter->firstdeletedpos = -1;

   /* update event mask */
   eventfilter->eventmask |= eventfilter->delayedeventmask;
   eventfilter->delayedeventmask = SCIP_EVENTTYPE_DISABLED;

   /* mark the event filter updates to be no longer delayed */
   eventfilter->delayupdates = FALSE;
}

/** processes the event with all event handlers with matching filter setting */
SCIP_RETCODE SCIPeventfilterProcess(
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENT*           event               /**< event to process */
   )
{
   SCIP_EVENTTYPE eventtype;
   SCIP_EVENTTYPE* eventtypes;
   SCIP_Bool processed;
   int len;
   int i;

   assert(eventfilter != NULL);
   assert(set != NULL);
   assert(event != NULL);

   SCIPsetDebugMsg(set, "processing event filter %p (len %d, mask 0x%" SCIP_EVENTTYPE_FORMAT ") with event type 0x%" SCIP_EVENTTYPE_FORMAT "\n",
      (void*)eventfilter, eventfilter->len, eventfilter->eventmask, event->eventtype);

   eventtype = event->eventtype;

   /* check, if there may be any event handler for specific event */
   if( (eventtype & eventfilter->eventmask) == 0 )
      return SCIP_OKAY;

   /* delay the updates on this eventfilter, such that changes during event processing to the event filter
    * don't destroy necessary information of the arrays we are currently using
    */
   eventfilterDelayUpdates(eventfilter);

   /* process the event by calling the event handlers */
   processed = FALSE;
   len = eventfilter->len;
   eventtypes = eventfilter->eventtypes;
   for( i = 0; i < len; ++i )
   {
      /* check, if event is applicable for the filter element */
      if( (eventtype & eventtypes[i]) != 0 )
      {
         /* call event handler */
         SCIP_CALL( SCIPeventhdlrExec(eventfilter->eventhdlrs[i], set, event, eventfilter->eventdata[i]) );
         processed = TRUE;
      }
   }

   /* update eventfilter mask, if event was not processed by any event handler */
   if( !processed )
   {
      eventfilter->eventmask &= ~event->eventtype;
      SCIPsetDebugMsg(set, " -> event type 0x%" SCIP_EVENTTYPE_FORMAT " not processed. new mask of event filter %p: 0x%" SCIP_EVENTTYPE_FORMAT "\n",
         event->eventtype, (void*)eventfilter, eventfilter->eventmask);
   }

   /* process delayed events on this eventfilter */
   eventfilterProcessUpdates(eventfilter);

   return SCIP_OKAY;
}



/*
 * Event queue methods
 */

/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPeventqueueIsDelayed

/** resizes events array to be able to store at least num entries */
static
SCIP_RETCODE eventqueueEnsureEventsMem(
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of node slots in array */
   )
{
   assert(eventqueue != NULL);
   assert(set != NULL);

   if( num > eventqueue->eventssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&eventqueue->events, newsize) );
      eventqueue->eventssize = newsize;
   }
   assert(num <= eventqueue->eventssize);

   return SCIP_OKAY;
}

/** creates an event queue */
SCIP_RETCODE SCIPeventqueueCreate(
   SCIP_EVENTQUEUE**     eventqueue          /**< pointer to store the event queue */
   )
{
   assert(eventqueue != NULL);

   SCIP_ALLOC( BMSallocMemory(eventqueue) );
   (*eventqueue)->events = NULL;
   (*eventqueue)->eventssize = 0;
   (*eventqueue)->nevents = 0;
   (*eventqueue)->delayevents = FALSE;

   return SCIP_OKAY;
}

/** frees event queue; there must not be any unprocessed events in the queue! */
SCIP_RETCODE SCIPeventqueueFree(
   SCIP_EVENTQUEUE**     eventqueue          /**< pointer to the event queue */
   )
{
   assert(eventqueue != NULL);
   assert(*eventqueue != NULL);
   assert((*eventqueue)->nevents == 0);

   BMSfreeMemoryArrayNull(&(*eventqueue)->events);
   BMSfreeMemory(eventqueue);

   return SCIP_OKAY;
}

/** appends event to the event queue; sets event to NULL afterwards */
static
SCIP_RETCODE eventqueueAppend(
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENT**          event               /**< pointer to event to append to the queue */
   )
{
   assert(eventqueue != NULL);
   assert(eventqueue->delayevents);
   assert(event != NULL);
   assert(*event != NULL);

   SCIPsetDebugMsg(set, "appending event %p of type 0x%" SCIP_EVENTTYPE_FORMAT " to event queue %p at position %d\n",
      (void*)*event, (*event)->eventtype, (void*)eventqueue, eventqueue->nevents);

   SCIP_CALL( eventqueueEnsureEventsMem(eventqueue, set, eventqueue->nevents+1) );
   eventqueue->events[eventqueue->nevents] = *event;
   eventqueue->nevents++;

   *event = NULL;

   return SCIP_OKAY;
}

/** processes event or adds event to the event queue */
SCIP_RETCODE SCIPeventqueueAdd(
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data; only needed for objchanged events, or NULL */
   SCIP_LP*              lp,                 /**< current LP data; only needed for obj/boundchanged events, or NULL */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for bound change events, or NULL */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events; not needed for variable specific events */
   SCIP_EVENT**          event               /**< pointer to event to add to the queue; will be NULL after queue addition */
   )
{
   SCIP_VAR* var;
   SCIP_EVENT* qevent;
   int pos;

   assert(eventqueue != NULL);
   assert(event != NULL);
   assert(*event != NULL);
   assert(((*event)->eventtype & SCIP_EVENTTYPE_OBJCHANGED) == 0 || primal != NULL);
   assert(((*event)->eventtype & (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED)) == 0 || lp != NULL);
   assert(((*event)->eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) == 0 || branchcand != NULL);

   if( !eventqueue->delayevents )
   {
      /* immediately process event */
      SCIP_CALL( SCIPeventProcess(*event, set, primal, lp, branchcand, eventfilter) );
      SCIP_CALL( SCIPeventFree(event, blkmem) );
   }
   else
   {
      /* delay processing of event by appending it to the event queue */
      SCIPsetDebugMsg(set, "adding event %p of type 0x%" SCIP_EVENTTYPE_FORMAT " to event queue %p\n", (void*)*event, (*event)->eventtype, (void*)eventqueue);

      switch( (*event)->eventtype )
      {
      case SCIP_EVENTTYPE_DISABLED:
         SCIPerrorMessage("cannot add a disabled event to the event queue\n");
         return SCIP_INVALIDDATA;

      case SCIP_EVENTTYPE_SYNC: /*lint !e30 !e142*/
      case SCIP_EVENTTYPE_VARADDED:
      case SCIP_EVENTTYPE_VARDELETED:
      case SCIP_EVENTTYPE_VARFIXED:
      case SCIP_EVENTTYPE_VARUNLOCKED:
      case SCIP_EVENTTYPE_GLBCHANGED:
      case SCIP_EVENTTYPE_GUBCHANGED:
      case SCIP_EVENTTYPE_PRESOLVEROUND:
      case SCIP_EVENTTYPE_NODEFOCUSED:
      case SCIP_EVENTTYPE_NODEFEASIBLE:
      case SCIP_EVENTTYPE_NODEINFEASIBLE:
      case SCIP_EVENTTYPE_NODEBRANCHED:
      case SCIP_EVENTTYPE_FIRSTLPSOLVED:
      case SCIP_EVENTTYPE_LPSOLVED:
      case SCIP_EVENTTYPE_POORSOLFOUND:
      case SCIP_EVENTTYPE_BESTSOLFOUND:
      case SCIP_EVENTTYPE_GHOLEADDED:
      case SCIP_EVENTTYPE_GHOLEREMOVED:
      case SCIP_EVENTTYPE_LHOLEADDED:
      case SCIP_EVENTTYPE_LHOLEREMOVED:
      case SCIP_EVENTTYPE_ROWADDEDSEPA: /* @todo remove previous DELETEDSEPA event */
      case SCIP_EVENTTYPE_ROWDELETEDSEPA: /* @todo remove previous ADDEDSEPA event */
      case SCIP_EVENTTYPE_ROWADDEDLP: /* @todo remove previous DELETEDLP event */
      case SCIP_EVENTTYPE_ROWDELETEDLP: /* @todo remove previous ADDEDLP event */
      case SCIP_EVENTTYPE_ROWCOEFCHANGED: /* @todo merge? */
      case SCIP_EVENTTYPE_ROWCONSTCHANGED: /* @todo merge with previous constchanged event */
      case SCIP_EVENTTYPE_ROWSIDECHANGED: /* @todo merge with previous sidechanged event */
         /* these events cannot (or need not) be merged; just add them to the queue */
         SCIP_CALL( eventqueueAppend(eventqueue, set, event) );
         break;

      case SCIP_EVENTTYPE_OBJCHANGED:
         /* changes in objective value may be merged with older changes in objective value */
         var = (*event)->data.eventobjchg.var;
         assert(var != NULL);
         pos = var->eventqueueindexobj;
         if( pos >= 0 )
         {
            /* the objective value change event already exists -> modify it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_OBJCHANGED);
            assert(qevent->data.eventobjchg.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.eventobjchg.oldobj, qevent->data.eventobjchg.newobj));

            SCIPsetDebugMsg(set, " -> merging OBJ event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               SCIPvarGetName((*event)->data.eventobjchg.var), (*event)->data.eventobjchg.oldobj,
               (*event)->data.eventobjchg.newobj,
               pos, SCIPvarGetName(qevent->data.eventobjchg.var), qevent->data.eventobjchg.oldobj, 
               qevent->data.eventobjchg.newobj);

            qevent->data.eventobjchg.newobj = (*event)->data.eventobjchg.newobj;
            if( qevent->data.eventobjchg.newobj == qevent->data.eventobjchg.oldobj ) /*lint !e777*/
            {
               /* the queued objective value change was reversed -> disable the event in the queue */
               eventDisable(qevent);
               var->eventqueueindexobj = -1;
               SCIPsetDebugMsg(set, " -> event disabled\n");
            }

            /* free the event that is of no use any longer */
            SCIP_CALL( SCIPeventFree(event, blkmem) );
         }
         else
         {
            /* the objective value change event doesn't exist -> add it to the queue, and remember the array index */
            var->eventqueueindexobj = eventqueue->nevents;
            SCIP_CALL( eventqueueAppend(eventqueue, set, event) );
         }
         break;

      case SCIP_EVENTTYPE_LBTIGHTENED:
      case SCIP_EVENTTYPE_LBRELAXED:
         /* changes in lower bound may be merged with older changes in lower bound */
         var = (*event)->data.eventbdchg.var;
         assert(var != NULL);
         pos = var->eventqueueindexlb;
         if( pos >= 0 )
         {
            /* the lower bound change event already exists -> modify it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_LBTIGHTENED || qevent->eventtype == SCIP_EVENTTYPE_LBRELAXED);
            assert(qevent->data.eventbdchg.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.eventbdchg.oldbound, qevent->data.eventbdchg.newbound));

            SCIPsetDebugMsg(set, " -> merging LB event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               SCIPvarGetName((*event)->data.eventbdchg.var), (*event)->data.eventbdchg.oldbound,
               (*event)->data.eventbdchg.newbound,
               pos, SCIPvarGetName(qevent->data.eventbdchg.var), qevent->data.eventbdchg.oldbound, 
               qevent->data.eventbdchg.newbound);

            qevent->data.eventbdchg.newbound = (*event)->data.eventbdchg.newbound;
            /*if( SCIPsetIsLT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            if( qevent->data.eventbdchg.newbound < qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_LBRELAXED;
            /*else if( SCIPsetIsGT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            else if( qevent->data.eventbdchg.newbound > qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_LBTIGHTENED;
            else
            {
               /* the queued bound change was reversed -> disable the event in the queue */
               assert(qevent->data.eventbdchg.newbound == qevent->data.eventbdchg.oldbound); /*lint !e777*/
               eventDisable(qevent);
               var->eventqueueindexlb = -1;
               SCIPsetDebugMsg(set, " -> event disabled\n");
            }

            /* free the event that is of no use any longer */
            SCIP_CALL( SCIPeventFree(event, blkmem) );
         }
         else
         {
            /* the lower bound change event doesn't exist -> add it to the queue, and remember the array index */
            var->eventqueueindexlb = eventqueue->nevents;
            SCIP_CALL( eventqueueAppend(eventqueue, set, event) );
         }
         break;

      case SCIP_EVENTTYPE_UBTIGHTENED:
      case SCIP_EVENTTYPE_UBRELAXED:
         /* changes in upper bound may be merged with older changes in upper bound */
         var = (*event)->data.eventbdchg.var;
         assert(var != NULL);
         pos = var->eventqueueindexub;
         if( pos >= 0 )
         {
            /* the upper bound change event already exists -> modify it accordingly */
            assert(pos < eventqueue->nevents);
            qevent = eventqueue->events[pos];
            assert(qevent != NULL);
            assert(qevent->eventtype == SCIP_EVENTTYPE_UBTIGHTENED || qevent->eventtype == SCIP_EVENTTYPE_UBRELAXED);
            assert(qevent->data.eventbdchg.var == var);
            assert(SCIPsetIsEQ(set, (*event)->data.eventbdchg.oldbound, qevent->data.eventbdchg.newbound));

            SCIPsetDebugMsg(set, " -> merging UB event (<%s>,%g -> %g) with event at position %d (<%s>,%g -> %g)\n",
               SCIPvarGetName((*event)->data.eventbdchg.var), (*event)->data.eventbdchg.oldbound,
               (*event)->data.eventbdchg.newbound,
               pos, SCIPvarGetName(qevent->data.eventbdchg.var), qevent->data.eventbdchg.oldbound,
               qevent->data.eventbdchg.newbound);

            qevent->data.eventbdchg.newbound = (*event)->data.eventbdchg.newbound;
            /*if( SCIPsetIsLT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            if( qevent->data.eventbdchg.newbound < qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_UBTIGHTENED;
            /*else if( SCIPsetIsGT(set, qevent->data.eventbdchg.newbound, qevent->data.eventbdchg.oldbound) )*/
            else if( qevent->data.eventbdchg.newbound > qevent->data.eventbdchg.oldbound )
               qevent->eventtype = SCIP_EVENTTYPE_UBRELAXED;
            else
            {
               /* the queued bound change was reversed -> disable the event in the queue */
               assert(qevent->data.eventbdchg.newbound == qevent->data.eventbdchg.oldbound); /*lint !e777*/
               eventDisable(qevent);
               var->eventqueueindexub = -1;
               SCIPsetDebugMsg(set, " -> event disabled\n");
            }

            /* free the event that is of no use any longer */
            SCIP_CALL( SCIPeventFree(event, blkmem) );
         }
         else
         {
            /* the upper bound change event doesn't exist -> add it to the queue, and remember the array index */
            var->eventqueueindexub = eventqueue->nevents;
            SCIP_CALL( eventqueueAppend(eventqueue, set, event) );
         }
         break;

      case SCIP_EVENTTYPE_IMPLADDED:
         var = (*event)->data.eventimpladd.var;
         assert(var != NULL);
         if( var->eventqueueimpl )
         {
            /* free the event that is of no use any longer */
            SCIP_CALL( SCIPeventFree(event, blkmem) );
         }
         else
         {
            var->eventqueueimpl = TRUE;
            SCIP_CALL( eventqueueAppend(eventqueue, set, event) );
         }
         break;

      default:
         SCIPerrorMessage("unknown event type <%d>\n", (*event)->eventtype);
         return SCIP_INVALIDDATA;
      }
   }

   assert(*event == NULL);

   return SCIP_OKAY;
}

/** marks queue to delay incoming events until a call to SCIPeventqueueProcess() */
SCIP_RETCODE SCIPeventqueueDelay(
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(eventqueue != NULL);
   assert(!eventqueue->delayevents);

   SCIPdebugMessage("event processing is delayed\n");

   eventqueue->delayevents = TRUE;

   return SCIP_OKAY;
}

/** processes all delayed events, marks queue to process events immediately */
SCIP_RETCODE SCIPeventqueueProcess(
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   )
{
   SCIP_EVENT* event;
   int i;

   assert(eventqueue != NULL);
   assert(eventqueue->delayevents);

   SCIPsetDebugMsg(set, "processing %d queued events\n", eventqueue->nevents);

   /* pass events to the responsible event filters
    * During event processing, new events may be raised. We have to loop to the mutable eventqueue->nevents.
    * A loop to something like "nevents = eventqueue->nevents; for(...; i < nevents; ...)" would miss the
    * newly created events. The same holds for eventqueue->events, which can be moved in memory due to
    * memory reallocation in eventqueueAppend().
    */
   for( i = 0; i < eventqueue->nevents; ++i )
   {
      event = eventqueue->events[i];
      assert(event != NULL);

      SCIPsetDebugMsg(set, "processing event %d of %d events in queue: eventtype=0x%" SCIP_EVENTTYPE_FORMAT "\n", i, eventqueue->nevents, event->eventtype);

      /* unmark the event queue index of a variable with changed objective value or bounds, and unmark the event queue
       * member flag of a variable with added implication
       */
      if( (event->eventtype & SCIP_EVENTTYPE_OBJCHANGED) != 0 )
      {
         assert(event->data.eventobjchg.var->eventqueueindexobj == i);
         event->data.eventobjchg.var->eventqueueindexobj = -1;
      }
      else if( (event->eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
      {
         assert(event->data.eventbdchg.var->eventqueueindexlb == i);
         event->data.eventbdchg.var->eventqueueindexlb = -1;
      }
      else if( (event->eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0 )
      {
         assert(event->data.eventbdchg.var->eventqueueindexub == i);
         event->data.eventbdchg.var->eventqueueindexub = -1;
      }
      else if( (event->eventtype & SCIP_EVENTTYPE_IMPLADDED) != 0 )
      {
         assert(event->data.eventimpladd.var->eventqueueimpl);
         event->data.eventimpladd.var->eventqueueimpl = FALSE;
      }

      /* process event */
      SCIP_CALL( SCIPeventProcess(event, set, primal, lp, branchcand, eventfilter) );

      /* free the event immediately, because additionally raised events during event processing
       * can lead to a large event queue
       */
      SCIP_CALL( SCIPeventFree(&eventqueue->events[i], blkmem) );
   }

   assert(i == eventqueue->nevents);
   eventqueue->nevents = 0;
   eventqueue->delayevents = FALSE;

   return SCIP_OKAY;
}

/** returns TRUE iff events of the queue are delayed until the next SCIPeventqueueProcess() call */
SCIP_Bool SCIPeventqueueIsDelayed(
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   assert(eventqueue != NULL);

   return eventqueue->delayevents;
}
