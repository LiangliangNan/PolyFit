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

/**@file   type_event.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for managing events
 * @author Tobias Achterberg
 * @author Leona Gottwald
 *
 *  This file defines the interface for event handler implemented in C.
 *
 *  - \ref scip::ObjEventhdlr "C++ wrapper class"
 */

/** @defgroup DEFPLUGINS_EVENT Default event handlers
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default event handlers of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_EVENT_H__
#define __SCIP_TYPE_EVENT_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#if !defined(_MSC_VER) || _MSC_VER > 1600
#ifdef __cplusplus
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#else
#define PRIx64 "llx"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * event types
 */

#define SCIP_EVENTTYPE_DISABLED         UINT64_C(0x000000000)  /**< the event was disabled and has no effect any longer */

/* variable events */
#define SCIP_EVENTTYPE_VARADDED         UINT64_C(0x000000001)  /**< a variable has been added to the transformed problem */
#define SCIP_EVENTTYPE_VARDELETED       UINT64_C(0x000000002)  /**< a variable will be deleted from the transformed problem */
#define SCIP_EVENTTYPE_VARFIXED         UINT64_C(0x000000004)  /**< a variable has been fixed, aggregated, or multi-aggregated */
#define SCIP_EVENTTYPE_VARUNLOCKED      UINT64_C(0x000000008)  /**< the number of rounding locks of a variable was reduced to zero or one */
#define SCIP_EVENTTYPE_OBJCHANGED       UINT64_C(0x000000010)  /**< the objective value of a variable has been changed */
#define SCIP_EVENTTYPE_GLBCHANGED       UINT64_C(0x000000020)  /**< the global lower bound of a variable has been changed */
#define SCIP_EVENTTYPE_GUBCHANGED       UINT64_C(0x000000040)  /**< the global upper bound of a variable has been changed */
#define SCIP_EVENTTYPE_LBTIGHTENED      UINT64_C(0x000000080)  /**< the local lower bound of a variable has been increased */
#define SCIP_EVENTTYPE_LBRELAXED        UINT64_C(0x000000100)  /**< the local lower bound of a variable has been decreased */
#define SCIP_EVENTTYPE_UBTIGHTENED      UINT64_C(0x000000200)  /**< the local upper bound of a variable has been decreased */
#define SCIP_EVENTTYPE_UBRELAXED        UINT64_C(0x000000400)  /**< the local upper bound of a variable has been increased */
#define SCIP_EVENTTYPE_GHOLEADDED       UINT64_C(0x000000800)  /**< a global hole has been added to the hole list of a variable's domain */
#define SCIP_EVENTTYPE_GHOLEREMOVED     UINT64_C(0x000001000)  /**< a global hole has been removed from the hole list of a variable's domain */
#define SCIP_EVENTTYPE_LHOLEADDED       UINT64_C(0x000002000)  /**< a local hole has been added to the hole list of a variable's domain */
#define SCIP_EVENTTYPE_LHOLEREMOVED     UINT64_C(0x000004000)  /**< a local hole has been removed from the hole list of a variable's domain */
#define SCIP_EVENTTYPE_IMPLADDED        UINT64_C(0x000008000)  /**< the variable's implication list, variable bound or clique information was extended */
#define SCIP_EVENTTYPE_TYPECHANGED      UINT64_C(0x000010000)  /**< the type of a variable has changed */

/* presolving events */
#define SCIP_EVENTTYPE_PRESOLVEROUND    UINT64_C(0x000020000)  /**< a presolving round has been finished */

/* node events */
#define SCIP_EVENTTYPE_NODEFOCUSED      UINT64_C(0x000040000)  /**< a node has been focused and is now the focus node */
#define SCIP_EVENTTYPE_NODEFEASIBLE     UINT64_C(0x000080000)  /**< the LP/pseudo solution of the node was feasible */
#define SCIP_EVENTTYPE_NODEINFEASIBLE   UINT64_C(0x000100000)  /**< the focus node has been proven to be infeasible or was bounded */
#define SCIP_EVENTTYPE_NODEBRANCHED     UINT64_C(0x000200000)  /**< the focus node has been solved by branching */
#define SCIP_EVENTTYPE_NODEDELETE       UINT64_C(0x000400000)  /**< a node is about to be deleted from the tree */


/* LP events */
#define SCIP_EVENTTYPE_FIRSTLPSOLVED    UINT64_C(0x000800000)  /**< the node's initial LP was solved */
#define SCIP_EVENTTYPE_LPSOLVED         UINT64_C(0x001000000)  /**< the node's LP was completely solved with cut & price */

/* primal solution events */
#define SCIP_EVENTTYPE_POORSOLFOUND     UINT64_C(0x002000000)  /**< a good enough primal feasible (but not new best) solution was found */
#define SCIP_EVENTTYPE_BESTSOLFOUND     UINT64_C(0x004000000)  /**< a new best primal feasible solution was found */

/* linear row events */
#define SCIP_EVENTTYPE_ROWADDEDSEPA     UINT64_C(0x008000000)  /**< a row has been added to SCIP's separation storage */
#define SCIP_EVENTTYPE_ROWDELETEDSEPA   UINT64_C(0x010000000)  /**< a row has been removed from SCIP's separation storage */
#define SCIP_EVENTTYPE_ROWADDEDLP       UINT64_C(0x020000000)  /**< a row has been added to the LP */
#define SCIP_EVENTTYPE_ROWDELETEDLP     UINT64_C(0x040000000)  /**< a row has been removed from the LP */
#define SCIP_EVENTTYPE_ROWCOEFCHANGED   UINT64_C(0x080000000)  /**< a coefficient of a row has been changed (row specific event) */
#define SCIP_EVENTTYPE_ROWCONSTCHANGED  UINT64_C(0x100000000)  /**< the constant of a row has been changed (row specific event) */
#define SCIP_EVENTTYPE_ROWSIDECHANGED   UINT64_C(0x200000000)  /**< a side of a row has been changed (row specific event) */

/* sync event */
#define SCIP_EVENTTYPE_SYNC             UINT64_C(0x400000000) /**< synchronization event */

/* event masks for variable events */
#define SCIP_EVENTTYPE_GBDCHANGED     (SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED)
#define SCIP_EVENTTYPE_LBCHANGED      (SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_LBRELAXED)
#define SCIP_EVENTTYPE_UBCHANGED      (SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_UBRELAXED)
#define SCIP_EVENTTYPE_BOUNDTIGHTENED (SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_UBTIGHTENED)
#define SCIP_EVENTTYPE_BOUNDRELAXED   (SCIP_EVENTTYPE_LBRELAXED | SCIP_EVENTTYPE_UBRELAXED)
#define SCIP_EVENTTYPE_BOUNDCHANGED   (SCIP_EVENTTYPE_LBCHANGED | SCIP_EVENTTYPE_UBCHANGED)
#define SCIP_EVENTTYPE_GHOLECHANGED   (SCIP_EVENTTYPE_GHOLEADDED | SCIP_EVENTTYPE_GHOLEREMOVED)
#define SCIP_EVENTTYPE_LHOLECHANGED   (SCIP_EVENTTYPE_LHOLEADDED | SCIP_EVENTTYPE_LHOLEREMOVED)
#define SCIP_EVENTTYPE_HOLECHANGED    (SCIP_EVENTTYPE_GHOLECHANGED | SCIP_EVENTTYPE_LHOLECHANGED)
#define SCIP_EVENTTYPE_DOMCHANGED     (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_HOLECHANGED)
#define SCIP_EVENTTYPE_VARCHANGED     (SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED | SCIP_EVENTTYPE_OBJCHANGED \
                                       | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_DOMCHANGED | SCIP_EVENTTYPE_IMPLADDED \
                                       | SCIP_EVENTTYPE_VARDELETED | SCIP_EVENTTYPE_TYPECHANGED)
#define SCIP_EVENTTYPE_VAREVENT       (SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_VARCHANGED | SCIP_EVENTTYPE_TYPECHANGED)

/* event masks for node events */
#define SCIP_EVENTTYPE_NODESOLVED     (SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEINFEASIBLE \
                                       | SCIP_EVENTTYPE_NODEBRANCHED)
#define SCIP_EVENTTYPE_NODEEVENT      (SCIP_EVENTTYPE_NODEFOCUSED | SCIP_EVENTTYPE_NODESOLVED)

/* event masks for LP events */
#define SCIP_EVENTTYPE_LPEVENT        (SCIP_EVENTTYPE_FIRSTLPSOLVED | SCIP_EVENTTYPE_LPSOLVED)

/* event masks for primal solution events */
#define SCIP_EVENTTYPE_SOLFOUND       (SCIP_EVENTTYPE_POORSOLFOUND | SCIP_EVENTTYPE_BESTSOLFOUND)
#define SCIP_EVENTTYPE_SOLEVENT       (SCIP_EVENTTYPE_SOLFOUND)

/* event masks for row events */
#define SCIP_EVENTTYPE_ROWCHANGED     (SCIP_EVENTTYPE_ROWCOEFCHANGED | SCIP_EVENTTYPE_ROWCONSTCHANGED | SCIP_EVENTTYPE_ROWSIDECHANGED)
#define SCIP_EVENTTYPE_ROWEVENT       (SCIP_EVENTTYPE_ROWADDEDSEPA | SCIP_EVENTTYPE_ROWDELETEDSEPA | SCIP_EVENTTYPE_ROWADDEDLP | SCIP_EVENTTYPE_ROWDELETEDLP | SCIP_EVENTTYPE_ROWCHANGED)

typedef uint64_t SCIP_EVENTTYPE;                  /**< type of event (bit field) */
#define SCIP_EVENTTYPE_FORMAT PRIx64

typedef struct SCIP_Eventhdlr SCIP_EVENTHDLR;     /**< event handler for a specific events */
typedef struct SCIP_EventhdlrData SCIP_EVENTHDLRDATA; /**< event handler data */
typedef struct SCIP_Event SCIP_EVENT;             /**< event data structure */
typedef struct SCIP_EventVarAdded SCIP_EVENTVARADDED; /**< data for variable addition events */
typedef struct SCIP_EventVarDeleted SCIP_EVENTVARDELETED; /**< data for variable deletion events */
typedef struct SCIP_EventVarFixed SCIP_EVENTVARFIXED; /**< data for variable fixing events */
typedef struct SCIP_EventVarUnlocked SCIP_EVENTVARUNLOCKED; /**< data for variable unlocked events */
typedef struct SCIP_EventObjChg SCIP_EVENTOBJCHG; /**< data for objective value change events */
typedef struct SCIP_EventBdChg SCIP_EVENTBDCHG;   /**< data for bound change events */
typedef struct SCIP_EventHole SCIP_EVENTHOLE;     /**< data for domain hole events */
typedef struct SCIP_EventImplAdd SCIP_EVENTIMPLADD; /**< data for implication added events */
typedef struct SCIP_EventTypeChg SCIP_EVENTTYPECHG; /**< data for variable type change events */
typedef struct SCIP_EventRowAddedSepa SCIP_EVENTROWADDEDSEPA; /**< data for row addition to sepastorage events */
typedef struct SCIP_EventRowDeletedSepa SCIP_EVENTROWDELETEDSEPA; /**< data for row deletion from sepastorage events */
typedef struct SCIP_EventRowAddedLP SCIP_EVENTROWADDEDLP; /**< data for row addition to LP events */
typedef struct SCIP_EventRowDeletedLP SCIP_EVENTROWDELETEDLP; /**< data for row deletion from LP events */
typedef struct SCIP_EventRowCoefChanged SCIP_EVENTROWCOEFCHANGED; /**< data for row coefficient change events */
typedef struct SCIP_EventRowConstChanged SCIP_EVENTROWCONSTCHANGED; /**< data for row constant change events */
typedef struct SCIP_EventRowSideChanged SCIP_EVENTROWSIDECHANGED; /**< data for row side change events */
typedef struct SCIP_EventData SCIP_EVENTDATA;     /**< locally defined event specific data */
typedef struct SCIP_EventFilter SCIP_EVENTFILTER; /**< event filter to select events to be processed by an event handler */
typedef struct SCIP_EventQueue SCIP_EVENTQUEUE;   /**< event queue to cache events and process them later */

/** copy method for event handler plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 */
#define SCIP_DECL_EVENTCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr)

/** destructor of event handler to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 */
#define SCIP_DECL_EVENTFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr)

/** initialization method of event handler (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 */
#define SCIP_DECL_EVENTINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr)

/** deinitialization method of event handler (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 */
#define SCIP_DECL_EVENTEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr)

/** solving process initialization method of event handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The event handler may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 */
#define SCIP_DECL_EVENTINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr)

/** solving process deinitialization method of event handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The event handler should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 */
#define SCIP_DECL_EVENTEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr)

/** frees specific event data
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 *  - eventdata       : pointer to the event data to free
 */
#define SCIP_DECL_EVENTDELETE(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr, SCIP_EVENTDATA** eventdata)

/** execution method of event handler
 *
 *  Processes the event. The method is called every time an event occurs, for which the event handler
 *  is responsible. Event handlers may declare themselves responsible for events by calling the
 *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
 *  given event handler and event data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - eventhdlr       : the event handler itself
 *  - event           : event to process
 *  - eventdata       : user data for the event
 */
#define SCIP_DECL_EVENTEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_EVENTHDLR* eventhdlr, SCIP_EVENT* event, SCIP_EVENTDATA* eventdata)

#ifdef __cplusplus
}
#endif

#endif
