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
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_compr.h"
#include "scip/pub_compr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given tree compression to a new scip */
extern
SCIP_RETCODE SCIPcomprCopyInclude(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a tree compression */
extern
SCIP_RETCODE SCIPcomprCreate(
   SCIP_COMPR**          compr,              /**< pointer to tree compression data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of tree compression */
   const char*           desc,               /**< description of tree compression */
   int                   priority,           /**< priority of the tree compression */
   int                   minnnodes,          /**< minimal number of nodes for calling compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy)),     /**< copy method of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_COMPRFREE   ((*comprfree)),     /**< destructor of tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit)),     /**< initialize tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit)),     /**< deinitialize tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol)), /**< solving process initialization method of tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol)), /**< solving process deinitialization method of tree compression */
   SCIP_DECL_COMPREXEC   ((*comprexec)),     /**< execution method of tree compression */
   SCIP_COMPRDATA*       comprdata           /**< tree compression data */
   );

/** calls destructor and frees memory of tree compression */
extern
SCIP_RETCODE SCIPcomprFree(
   SCIP_COMPR**          compr,              /**< pointer to tree compression data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes tree compression */
extern
SCIP_RETCODE SCIPcomprInit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of tree compression */
extern
SCIP_RETCODE SCIPcomprExit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs tree compression that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPcomprInitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs tree compression that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPcomprExitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of tree compression */
extern
SCIP_RETCODE SCIPcomprExec(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of tree compression */
extern
void SCIPcomprSetPriority(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the tree compression */
   );

/** sets copy callback of tree compression */
extern
void SCIPcomprSetCopy(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy))      /**< copy callback of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor callback of tree compression */
extern
void SCIPcomprSetFree(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRFREE   ((*comprfree))      /**< destructor of tree compression */
   );

/** sets initialization callback of tree compression */
extern
void SCIPcomprSetInit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit))      /**< initialize tree compression */
   );

/** sets deinitialization callback of tree compression */
extern
void SCIPcomprSetExit(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit))      /**< deinitialize tree compression */
   );

/** sets solving process initialization callback of tree compression */
extern
void SCIPcomprSetInitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol))  /**< solving process initialization callback of tree compression */
   );

/** sets solving process deinitialization callback of tree compression */
extern
void SCIPcomprSetExitsol(
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol))  /**< solving process deinitialization callback of tree compression */
   );

/** should the compression be executed at the given depth, frequency, timing, ... */
EXTERN
SCIP_Bool SCIPcomprShouldBeExecuted(
   SCIP_COMPR*           compr,              /**< tree compression */
   int                   depth,              /**< depth of current node */
   int                   nnodes              /**< number of open nodes */
   );

#ifdef __cplusplus
}
#endif

#endif
