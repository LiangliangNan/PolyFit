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

/**@file   event.h
 * @brief  internal methods for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_H__
#define __SCIP_EVENT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"
#include "scip/type_branch.h"
#include "scip/pub_event.h"

#include "scip/struct_event.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Event handler methods
 */

/** copies the given event handler to a new scip */
extern
SCIP_RETCODE SCIPeventhdlrCopyInclude(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates an event handler */
extern
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
   );

/** calls destructor and frees memory of event handler */
extern
SCIP_RETCODE SCIPeventhdlrFree(
   SCIP_EVENTHDLR**      eventhdlr,          /**< pointer to event handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes event handler */
extern
SCIP_RETCODE SCIPeventhdlrInit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of event handler */
extern
SCIP_RETCODE SCIPeventhdlrExit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for this event */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs event handler that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPeventhdlrInitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs event handler that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPeventhdlrExitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of event handler */
extern
SCIP_RETCODE SCIPeventhdlrExec(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENT*           event,              /**< event to call event handler with */
   SCIP_EVENTDATA*       eventdata           /**< user data for the issued event */
   );

/**
 * callback setter methods of event handlers
 */
/** sets copy callback for all events of this event handler */
extern
void SCIPeventhdlrSetCopy(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy))      /**< copy callback for events */
   );

/** sets destructor callback of this event handler */
extern
void SCIPeventhdlrSetFree(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTFREE   ((*eventfree))      /**< destructor callback of event handler */
   );

/** sets initialization callback of this event handler */
extern
void SCIPeventhdlrSetInit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit))      /**< initialization callback of event handler */
   );

/** sets deinitialization callback of this event handler */
extern
void SCIPeventhdlrSetExit(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit))      /**< deinitialization callback of event handler */
   );

/** sets solving process initialization callback of this event handler */
extern
void SCIPeventhdlrSetInitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol))   /**< solving process initialization callback of event handler */
   );

/** sets solving process deinitialization callback of this event handler */
extern
void SCIPeventhdlrSetExitsol(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol))   /**< solving process deinitialization callback of event handler */
   );

/** sets callback to free specific event data */
extern
void SCIPeventhdlrSetDelete(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete))    /**< callback to free specific event data */
   );

/** enables or disables all clocks of \p eventhdlr, depending on the value of the flag */
extern
void SCIPeventhdlrEnableOrDisableClocks(
   SCIP_EVENTHDLR*       eventhdlr,          /**< the event handler for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the event handler be enabled? */
   );

/*
 * Event methods
 */

/** creates a synchronization event */
extern
SCIP_RETCODE SCIPeventCreateSync(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates an event for an addition of a variable to the problem */
extern
SCIP_RETCODE SCIPeventCreateVarAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that was added to the problem */
   );

/** creates an event for a deletion of a variable from the problem */
extern
SCIP_RETCODE SCIPeventCreateVarDeleted(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that is to be deleted from the problem */
   );

/** creates an event for a fixing of a variable */
extern
SCIP_RETCODE SCIPeventCreateVarFixed(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that was fixed */
   );

/** creates an event for a change in the number of locks of a variable down to zero or one */
extern
SCIP_RETCODE SCIPeventCreateVarUnlocked(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that changed the number of locks */
   );

/** creates an event for a change in the objective value of a variable */
extern
SCIP_RETCODE SCIPeventCreateObjChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose objective value changed */
   SCIP_Real             oldobj,             /**< old objective value before value changed */
   SCIP_Real             newobj              /**< new objective value after value changed */
   );

/** creates an event for a change in the global lower bound of a variable */
extern
SCIP_RETCODE SCIPeventCreateGlbChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   );

/** creates an event for a change in the global upper bound of a variable */
extern
SCIP_RETCODE SCIPeventCreateGubChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   );

/** creates an event for a change in the lower bound of a variable */
extern
SCIP_RETCODE SCIPeventCreateLbChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   );

/** creates an event for a change in the upper bound of a variable */
extern
SCIP_RETCODE SCIPeventCreateUbChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             oldbound,           /**< old bound before bound changed */
   SCIP_Real             newbound            /**< new bound after bound changed */
   );

/** creates an event for an addition of a global domain hole to a variable */
SCIP_RETCODE SCIPeventCreateGholeAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   );

/** creates an event for removing a global domain hole of a variable */
SCIP_RETCODE SCIPeventCreateGholeRemoved(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in hole */
   SCIP_Real             right               /**< right bound of open interval in hole */
   );

/** creates an event for an addition of a local domain hole to a variable */
SCIP_RETCODE SCIPeventCreateLholeAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   );

/** creates an event for removing a local domain hole of a variable */
SCIP_RETCODE SCIPeventCreateLholeRemoved(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable whose bound changed */
   SCIP_Real             left,               /**< left bound of open interval in hole */
   SCIP_Real             right               /**< right bound of open interval in hole */
   );

/** creates an event for an addition to the variable's implications list, clique or variable bounds information */
extern
SCIP_RETCODE SCIPeventCreateImplAdded(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var                 /**< variable that was fixed */
   );

/** creates an event for the addition of a linear row to the separation storage */
extern
SCIP_RETCODE SCIPeventCreateRowAddedSepa(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was added to the separation storage*/
   );

/** creates an event for the deletion of a linear row from the separation storage */
extern
SCIP_RETCODE SCIPeventCreateRowDeletedSepa(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was deleted from the separation storage */
   );

/** creates an event for the addition of a linear row to the LP */
extern
SCIP_RETCODE SCIPeventCreateRowAddedLP(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was added to the LP */
   );

/** creates an event for the deletion of a linear row from the LP */
extern
SCIP_RETCODE SCIPeventCreateRowDeletedLP(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row                 /**< row that was deleted from the LP */
   );

/** creates an event for the change of a coefficient in a linear row */
extern
SCIP_RETCODE SCIPeventCreateRowCoefChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row in which a coefficient changed */
   SCIP_COL*             col,                /**< column which coefficient changed */
   SCIP_Real             oldval,             /**< old value of coefficient */
   SCIP_Real             newval              /**< new value of coefficient */
   );

/** creates an event for the change of a constant in a linear row */
extern
SCIP_RETCODE SCIPeventCreateRowConstChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row in which the constant changed */
   SCIP_Real             oldval,             /**< old value of constant */
   SCIP_Real             newval              /**< new value of constant */
   );

/** creates an event for the change of a side of a linear row */
extern
SCIP_RETCODE SCIPeventCreateRowSideChanged(
   SCIP_EVENT**          event,              /**< pointer to store the event */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_ROW*             row,                /**< row which side has changed */
   SCIP_SIDETYPE         side,               /**< which side has changed */
   SCIP_Real             oldval,             /**< old value of side */
   SCIP_Real             newval              /**< new value of side */
   );

/** frees an event */
extern
SCIP_RETCODE SCIPeventFree(
   SCIP_EVENT**          event,              /**< event to free */
   BMS_BLKMEM*           blkmem              /**< block memory buffer */
   );

/** sets type of event */
extern
SCIP_RETCODE SCIPeventChgType(
   SCIP_EVENT*           event,              /**< event */
   SCIP_EVENTTYPE        eventtype           /**< new event type */
   );

/** sets variable for a variable event */
extern
SCIP_RETCODE SCIPeventChgVar(
   SCIP_EVENT*           event,              /**< event */
   SCIP_VAR*             var                 /**< new variable */
   );

/** sets node for a node or LP event */
extern
SCIP_RETCODE SCIPeventChgNode(
   SCIP_EVENT*           event,              /**< event */
   SCIP_NODE*            node                /**< new node */
   );

/** sets solution for a primal solution event */
extern
SCIP_RETCODE SCIPeventChgSol(
   SCIP_EVENT*           event,              /**< event */
   SCIP_SOL*             sol                 /**< new primal solution */
   );

/** processes event by calling the appropriate event handlers */
extern
SCIP_RETCODE SCIPeventProcess(
   SCIP_EVENT*           event,              /**< event */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data; only needed for objchanged events, or NULL */
   SCIP_LP*              lp,                 /**< current LP data; only needed for obj/boundchanged events, or NULL */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for bound change events, or NULL */
   SCIP_EVENTFILTER*     eventfilter         /**< event filter for global events; not needed for variable specific events */
   );



/*
 * Event filter methods
 */

/** creates an event filter */
extern
SCIP_RETCODE SCIPeventfilterCreate(
   SCIP_EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   BMS_BLKMEM*           blkmem              /**< block memory buffer */
   );

/** frees an event filter and the associated event data entries */
extern
SCIP_RETCODE SCIPeventfilterFree(
   SCIP_EVENTFILTER**    eventfilter,        /**< pointer to store the event filter */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** adds element to event filter */
extern
SCIP_RETCODE SCIPeventfilterAdd(
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   );

/** deletes element from event filter */
extern
SCIP_RETCODE SCIPeventfilterDel(
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int                   filterpos           /**< position of event filter entry, or -1 if unknown */
   );

/** processes the event with all event handlers with matching filter setting */
extern
SCIP_RETCODE SCIPeventfilterProcess(
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENT*           event               /**< event to process */
   );



/*
 * Event queue methods
 */

/** creates an event queue */
extern
SCIP_RETCODE SCIPeventqueueCreate(
   SCIP_EVENTQUEUE**     eventqueue          /**< pointer to store the event queue */
   );

/** frees event queue; there must not be any unprocessed events in the queue! */
extern
SCIP_RETCODE SCIPeventqueueFree(
   SCIP_EVENTQUEUE**     eventqueue          /**< pointer to the event queue */
   );

/** processes event or adds event to the event queue */
extern
SCIP_RETCODE SCIPeventqueueAdd(
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data; only needed for objchanged events, or NULL */
   SCIP_LP*              lp,                 /**< current LP data; only needed for obj/boundchanged events, or NULL */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage; only needed for bound change events, or NULL */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events; not needed for variable specific events */
   SCIP_EVENT**          event               /**< pointer to event to add to the queue; will be NULL after queue addition */
   );

/** marks queue to delay incoming events until a call to SCIPeventqueueProcess() */
extern
SCIP_RETCODE SCIPeventqueueDelay(
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** processes all events in the queue */
extern
SCIP_RETCODE SCIPeventqueueProcess(
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter         /**< event filter for global (not variable dependent) events */
   );

/** returns TRUE iff events of the queue are delayed until the next SCIPeventqueueProcess() call */
extern
SCIP_Bool SCIPeventqueueIsDelayed(
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPeventqueueIsDelayed(eventqueue)       ((eventqueue)->delayevents)

#endif

#ifdef __cplusplus
}
#endif

#endif
