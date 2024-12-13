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

/**@file   compr.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for tree compressions
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_COMPR_H__
#define __SCIP_COMPR_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_reopt.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_compr.h"
#include "scip/pub_compr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given tree compression to a new scip */
SCIP_RETCODE SCIPcomprCopyInclude(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a tree compression */
SCIP_RETCODE SCIPcomprCreate(
   SCIP_COMPR**          compr,              /**< pointer to tree compression data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of tree compression */
   const char*           desc,               /**< description of tree compression */
   int                   priority,           /**< priority of the tree compression */
   int                   minnnodes,          /**< minimal number of nodes for calling compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy)),     /**< copy method of tree compression or NULL if you don't want to copy
                                              *   your plugin into sub-SCIPs */
   SCIP_DECL_COMPRFREE   ((*comprfree)),     /**< destructor of tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit)),     /**< initialize tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit)),     /**< deinitialize tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol)), /**< solving process initialization method of tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol)), /**< solving process deinitialization method of tree compression */
   SCIP_DECL_COMPREXEC   ((*comprexec)),     /**< execution method of tree compression */
   SCIP_COMPRDATA*       comprdata           /**< tree compression data */
   );

/** calls destructor and frees memory of tree compression */
SCIP_RETCODE SCIPcomprFree(
   SCIP_COMPR**          compr,              /**< pointer to tree compression data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes tree compression */
SCIP_RETCODE SCIPcomprInit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of tree compression */
SCIP_RETCODE SCIPcomprExit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs tree compression that the branch and bound process is being started */
SCIP_RETCODE SCIPcomprInitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs tree compression that the branch and bound process data is being freed */
SCIP_RETCODE SCIPcomprExitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of tree compression */
SCIP_RETCODE SCIPcomprExec(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of tree compression */
void SCIPcomprSetPriority(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the tree compression */
   );

/** sets copy callback of tree compression */
void SCIPcomprSetCopy(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy))      /**< copy callback of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor callback of tree compression */
void SCIPcomprSetFree(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRFREE   ((*comprfree))      /**< destructor of tree compression */
   );

/** sets initialization callback of tree compression */
void SCIPcomprSetInit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit))      /**< initialize tree compression */
   );

/** sets deinitialization callback of tree compression */
void SCIPcomprSetExit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit))      /**< deinitialize tree compression */
   );

/** sets solving process initialization callback of tree compression */
void SCIPcomprSetInitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol))  /**< solving process initialization callback of tree compression */
   );

/** sets solving process deinitialization callback of tree compression */
void SCIPcomprSetExitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol))  /**< solving process deinitialization callback of tree compression */
   );

/** should the compression be executed at the given depth, frequency, timing, ... */
SCIP_EXPORT
SCIP_Bool SCIPcomprShouldBeExecuted(
   SCIP_COMPR*           compr,              /**< tree compression */
   int                   depth,              /**< depth of current node */
   int                   nnodes              /**< number of open nodes */
   );

#ifdef __cplusplus
}
#endif

#endif
