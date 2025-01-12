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

/**@file   struct_event.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_EVENT_H__
#define __SCIP_STRUCT_EVENT_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_event.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data for variable addition events */
struct SCIP_EventVarAdded
{
   SCIP_VAR*             var;                /**< variable that was added to the problem */
};

/** data for variable deletion events */
struct SCIP_EventVarDeleted
{
   SCIP_VAR*             var;                /**< variable that will be deleted from the problem */
};

/** data for variable fixing events */
struct SCIP_EventVarFixed
{
   SCIP_VAR*             var;                /**< variable that was fixed */
};

/** data for locks change events */
struct SCIP_EventVarUnlocked
{
   SCIP_VAR*             var;                /**< variable for which the lock numbers were changed */
};

/** data for objective value change events */
struct SCIP_EventObjChg
{
   SCIP_Real             oldobj;             /**< old objective value before value changed */
   SCIP_Real             newobj;             /**< new objective value after value changed */
   SCIP_VAR*             var;                /**< variable whose objective value changed */
};

/** data for bound change events */
struct SCIP_EventBdChg
{
   SCIP_Real             oldbound;           /**< old bound before bound changed */
   SCIP_Real             newbound;           /**< new bound after bound changed */
   SCIP_VAR*             var;                /**< variable whose bound changed */
};

/** data for domain hole events */
struct SCIP_EventHole
{
   SCIP_Real             left;               /**< left bound of open interval in hole */
   SCIP_Real             right;              /**< right bound of open interval in hole */
   SCIP_VAR*             var;                /**< variable for which a hole was removed */
};

/** data for implication added events */
struct SCIP_EventImplAdd
{
   SCIP_VAR*             var;                /**< variable for which an implication, variable bound, or clique was added */
};

/** data for variable type change events */
struct SCIP_EventTypeChg
{
   SCIP_VARTYPE          oldtype;            /**< old variable type */
   SCIP_VARTYPE          newtype;            /**< new variable type */
   SCIP_VAR*             var;                /**< variable whose type changed */
};

/** data for row addition to separation storage events */
struct SCIP_EventRowAddedSepa
{
   SCIP_ROW*             row;                /**< row that was added to separation storage */
};

/** data for row deletion from separation storage events */
struct SCIP_EventRowDeletedSepa
{
   SCIP_ROW*             row;                /**< row that was deleted from separation storage */
};

/** data for row addition to LP events */
struct SCIP_EventRowAddedLP
{
   SCIP_ROW*             row;                /**< row that was added to the LP */
};

/** data for row deletion from LP events */
struct SCIP_EventRowDeletedLP
{
   SCIP_ROW*             row;                /**< row that was deleted from the LP */
};

/** data for row coefficient change events */
struct SCIP_EventRowCoefChanged
{
   SCIP_ROW*             row;                /**< row which coefficient has changed */
   SCIP_COL*             col;                /**< column which coefficient has changed */
   SCIP_Real             oldval;             /**< old value of coefficient */
   SCIP_Real             newval;             /**< new value of coefficient */
};

/** data for row constant change events */
struct SCIP_EventRowConstChanged
{
   SCIP_ROW*             row;                /**< row which constant has changed */
   SCIP_Real             oldval;             /**< old value of constant */
   SCIP_Real             newval;             /**< new value of constant */
};

/** data for row side change events */
struct SCIP_EventRowSideChanged
{
   SCIP_ROW*             row;                /**< row which side has changed */
   SCIP_SIDETYPE         side;               /**< which side has changed */
   SCIP_Real             oldval;             /**< old value of side */
   SCIP_Real             newval;             /**< new value of side */
};

/** event data structure */
struct SCIP_Event
{
   union
   {
      SCIP_EVENTVARADDED eventvaradded;       /**< data for variable addition events */
      SCIP_EVENTVARDELETED eventvardeleted;   /**< data for variable deletion events */
      SCIP_EVENTVARFIXED eventvarfixed;       /**< data for variable fixing events */
      SCIP_EVENTVARUNLOCKED eventvarunlocked; /**< data for locks change events */
      SCIP_EVENTOBJCHG   eventobjchg;         /**< data for objective value change events */
      SCIP_EVENTBDCHG    eventbdchg;          /**< data for bound change events */
      SCIP_EVENTHOLE     eventhole;           /**< data for domain hole events */
      SCIP_EVENTIMPLADD  eventimpladd;        /**< data for implication added events */
      SCIP_EVENTTYPECHG  eventtypechg;        /**< data for variable type change events */
      SCIP_EVENTROWADDEDSEPA eventrowaddedsepa; /**< data for row addition to separation storage events */
      SCIP_EVENTROWDELETEDSEPA eventrowdeletedsepa; /**< data for row deletion from separation storage events */
      SCIP_EVENTROWADDEDLP eventrowaddedlp;   /**< data for row addition to LP events */
      SCIP_EVENTROWDELETEDLP eventrowdeletedlp; /**< data for row deletion from LP events */
      SCIP_EVENTROWCOEFCHANGED eventrowcoefchanged; /**< data for row coefficient change events */
      SCIP_EVENTROWCONSTCHANGED eventrowconstchanged; /**< data for row constant change events */
      SCIP_EVENTROWSIDECHANGED eventrowsidechanged; /**< data for row side change events */
      SCIP_NODE*         node;                /**< data for node and LP events */
      SCIP_SOL*          sol;                 /**< data for primal solution events */
   } data;
   SCIP_EVENTTYPE        eventtype;          /**< type of event */
};

/** event filter to select events to be processed by an event handler */
struct SCIP_EventFilter
{
   SCIP_EVENTTYPE*       eventtypes;         /**< array with types of event to process; 0 marks a deleted event catch entry */
   SCIP_EVENTHDLR**      eventhdlrs;         /**< array with event handlers to process the event */
   SCIP_EVENTDATA**      eventdata;          /**< array with user data for the issued event */
   int*                  nextpos;            /**< linked lists for free, delayed added and delayed deleted slot positions */
   int                   size;               /**< size of filter arrays (available slots in arrays) */
   int                   len;                /**< number entries in filter arrays (used and deleted) */
   int                   firstfreepos;       /**< first deleted slot; remaining slots are in poslist */
   int                   firstdeletedpos;    /**< first delayed deleted slot; remaining slots are in poslist */
   SCIP_EVENTTYPE        eventmask;          /**< mask for events that are handled by any event handler in the filter */
   SCIP_EVENTTYPE        delayedeventmask;   /**< mask for delayed added events */
   SCIP_Bool             delayupdates;       /**< should additions and deletions to the filter be delayed? */
};

/** event handler */
struct SCIP_Eventhdlr
{
   char*                 name;               /**< name of event handler */
   char*                 desc;               /**< description of event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy));     /**< copy method of event handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_EVENTFREE   ((*eventfree));     /**< destructor of event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit));     /**< initialize event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit));     /**< deinitialize event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol));  /**< solving process initialization method of event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol));  /**< solving process deinitialization method of event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete));   /**< free specific event data */
   SCIP_DECL_EVENTEXEC   ((*eventexec));     /**< execute event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata;      /**< event handler data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this event handler for the next stages */
   SCIP_CLOCK*           eventtime;          /**< time spend in this event handler execution method */
   SCIP_Bool             initialized;        /**< is event handler initialized? */
};

/** event queue to cache events and process them later */
struct SCIP_EventQueue
{
   SCIP_EVENT**          events;             /**< array with queued events */
   int                   eventssize;         /**< number of available slots in events array */
   int                   nevents;            /**< number of events in queue (used slots if events array) */
   SCIP_Bool             delayevents;        /**< should the events be delayed and processed later? */
};

#ifdef __cplusplus
}
#endif

#endif
