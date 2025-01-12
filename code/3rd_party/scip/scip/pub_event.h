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

/**@file   pub_event.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for managing events
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_EVENT_H__
#define __SCIP_PUB_EVENT_H__


#include "scip/def.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"

/* In optimized mode, some function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms. For this, we need to include struct_event.h.
 */
#ifdef NDEBUG
#include "scip/struct_event.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Event handler methods
 */

/**@addtogroup PublicEventHandlerMethods
 *
 * @{
 */

/** gets name of event handler */
SCIP_EXPORT
const char* SCIPeventhdlrGetName(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** gets user data of event handler */
SCIP_EXPORT
SCIP_EVENTHDLRDATA* SCIPeventhdlrGetData(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** sets user data of event handler; user has to free old data in advance! */
SCIP_EXPORT
void SCIPeventhdlrSetData(
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< new event handler user data */
   );

/** is event handler initialized? */
SCIP_EXPORT
SCIP_Bool SCIPeventhdlrIsInitialized(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** gets time in seconds used in this event handler for setting up for next stages */
SCIP_EXPORT
SCIP_Real SCIPeventhdlrGetSetupTime(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** gets time in seconds used in this event handler */
SCIP_EXPORT
SCIP_Real SCIPeventhdlrGetTime(
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** @} */

/*
 * Event methods
 */

/**@addtogroup PublicEventMethods
 *
 * @{
 */

/** gets type of event */
SCIP_EXPORT
SCIP_EVENTTYPE SCIPeventGetType(
   SCIP_EVENT*           event               /**< event */
   );

/** gets variable for a variable event (var added, var deleted, var fixed, 
 *  objective value or domain change, domain hole added or removed) */
SCIP_EXPORT
SCIP_VAR* SCIPeventGetVar(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old objective value for an objective value change event */
SCIP_EXPORT
SCIP_Real SCIPeventGetOldobj(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new objective value for an objective value change event */
SCIP_EXPORT
SCIP_Real SCIPeventGetNewobj(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old bound for a bound change event */
SCIP_EXPORT
SCIP_Real SCIPeventGetOldbound(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new bound for a bound change event */
SCIP_EXPORT
SCIP_Real SCIPeventGetNewbound(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old variable type for a variable type change event */
SCIP_EXPORT
SCIP_VARTYPE SCIPeventGetOldtype(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new variable type for a variable type change event */
SCIP_EXPORT
SCIP_VARTYPE SCIPeventGetNewtype(
   SCIP_EVENT*           event               /**< event */
   );

/** gets node for a node or LP event */
SCIP_EXPORT
SCIP_NODE* SCIPeventGetNode(
   SCIP_EVENT*           event               /**< event */
   );

/** gets solution for a primal solution event */
SCIP_EXPORT
SCIP_SOL* SCIPeventGetSol(
   SCIP_EVENT*           event               /**< event */
   );

/** gets the left bound of open interval in the hole */
SCIP_EXPORT
SCIP_Real SCIPeventGetHoleLeft(
   SCIP_EVENT*           event               /**< event */
   );

/** gets the right bound of open interval in the hole */
SCIP_EXPORT
SCIP_Real SCIPeventGetHoleRight(
   SCIP_EVENT*           event               /**< event */
   );

/** gets row for a row event */
SCIP_EXPORT
SCIP_ROW* SCIPeventGetRow(
   SCIP_EVENT*           event               /**< event */
   );

/** gets column for a row change coefficient event */
SCIP_EXPORT
SCIP_COL* SCIPeventGetRowCol(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old coefficient value for a row change coefficient event */
SCIP_EXPORT
SCIP_Real SCIPeventGetRowOldCoefVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new coefficient value for a row change coefficient event */
SCIP_EXPORT
SCIP_Real SCIPeventGetRowNewCoefVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old constant value for a row change constant event */
SCIP_EXPORT
SCIP_Real SCIPeventGetRowOldConstVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new constant value for a row change constant event */
SCIP_EXPORT
SCIP_Real SCIPeventGetRowNewConstVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets side for a row change side event */
SCIP_EXPORT
SCIP_SIDETYPE SCIPeventGetRowSide(
   SCIP_EVENT*           event               /**< event */
   );

/** gets old side value for a row change side event */
SCIP_EXPORT
SCIP_Real SCIPeventGetRowOldSideVal(
   SCIP_EVENT*           event               /**< event */
   );

/** gets new side value for a row change side event */
SCIP_EXPORT
SCIP_Real SCIPeventGetRowNewSideVal(
   SCIP_EVENT*           event               /**< event */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPeventGetType(event)                   ((event)->eventtype)
#define SCIPeventGetOldobj(event)                 ((event)->data.eventobjchg.oldobj)
#define SCIPeventGetNewobj(event)                 ((event)->data.eventobjchg.newobj)
#define SCIPeventGetOldbound(event)               ((event)->data.eventbdchg.oldbound)
#define SCIPeventGetNewbound(event)               ((event)->data.eventbdchg.newbound)
#define SCIPeventGetOldtype(event)                ((event)->data.eventtypechg.oldtype)
#define SCIPeventGetNewtype(event)                ((event)->data.eventtypechg.newtype)
#define SCIPeventGetNode(event)                   ((event)->data.node)
#define SCIPeventGetSol(event)                    ((event)->data.sol)
#define SCIPeventGetRowCol(event)                 ((event)->data.eventrowcoefchanged.col)
#define SCIPeventGetRowOldCoefVal(event)          ((event)->data.eventrowcoefchanged.oldval)
#define SCIPeventGetRowNewCoefVal(event)          ((event)->data.eventrowcoefchanged.newval)
#define SCIPeventGetRowOldConstVal(event)         ((event)->data.eventrowconstchanged.oldval)
#define SCIPeventGetRowNewConstVal(event)         ((event)->data.eventrowconstchanged.newval)
#define SCIPeventGetRowSide(event)                ((event)->data.eventrowsidechanged.side)
#define SCIPeventGetRowOldSideVal(event)          ((event)->data.eventrowsidechanged.oldval)
#define SCIPeventGetRowNewSideVal(event)          ((event)->data.eventrowsidechanged.newval)

#endif

/** @} */

#ifdef __cplusplus
}
#endif

#endif
