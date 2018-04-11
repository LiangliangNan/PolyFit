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

/**@file   pub_nodesel.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_NODESEL_H__
#define __SCIP_PUB_NODESEL_H__


#include "scip/def.h"
#include "scip/type_nodesel.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNodeSelectorMethods
 *
 * @{
 */

/** gets name of node selector */
EXTERN
const char* SCIPnodeselGetName(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** gets description of node selector */
EXTERN
const char* SCIPnodeselGetDesc(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** gets priority of node selector in standard mode */
EXTERN
int SCIPnodeselGetStdPriority(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** gets priority of node selector in memory saving mode */
EXTERN
int SCIPnodeselGetMemsavePriority(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** gets user data of node selector */
EXTERN
SCIP_NODESELDATA* SCIPnodeselGetData(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** sets user data of node selector; user has to free old data in advance! */
EXTERN
void SCIPnodeselSetData(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_NODESELDATA*     nodeseldata         /**< new node selector user data */
   );

/** is node selector initialized? */
EXTERN
SCIP_Bool SCIPnodeselIsInitialized(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** gets time in seconds used in this node selector for setting up for next stages */
EXTERN
SCIP_Real SCIPnodeselGetSetupTime(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/** gets time in seconds used in this node selector */
EXTERN
SCIP_Real SCIPnodeselGetTime(
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
