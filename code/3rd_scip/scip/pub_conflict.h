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

/**@file   pub_conflict.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for conflict analysis handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CONFLICT_H__
#define __SCIP_PUB_CONFLICT_H__



#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_conflict.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicConflictMethods
 *
 * @{
 */

/** compares two conflict handlers w. r. to their priority */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrComp);

/** comparison method for sorting conflict handler w.r.t. to their name */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrCompName);

/** gets user data of conflict handler */
EXTERN
SCIP_CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** sets user data of conflict handler; user has to free old data in advance! */
EXTERN
void SCIPconflicthdlrSetData(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   );

/** gets name of conflict handler */
EXTERN
const char* SCIPconflicthdlrGetName(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets description of conflict handler */
EXTERN
const char* SCIPconflicthdlrGetDesc(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets priority of conflict handler */
EXTERN
int SCIPconflicthdlrGetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** is conflict handler initialized? */
EXTERN
SCIP_Bool SCIPconflicthdlrIsInitialized(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets time in seconds used in this conflict handler for setting up for next stages */
EXTERN
SCIP_Real SCIPconflicthdlrGetSetupTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** gets time in seconds used in this conflict handler */
EXTERN
SCIP_Real SCIPconflicthdlrGetTime(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
