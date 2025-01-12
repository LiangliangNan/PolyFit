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

/**@file   scip_compr.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for compression plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_COMPR_H__
#define __SCIP_SCIP_COMPR_H__


#include "scip/def.h"
#include "scip/type_compr.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCompressionMethods
 *
 * @{
 */
/** creates a tree compression and includes it in SCIP.
 *
 *  @note method has all compression callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeComprBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCompr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of tree compression */
   const char*           desc,               /**< description of tree compression */
   int                   priority,           /**< priority of the tree compression */
   int                   minnnodes,          /**< minimal number of nodes to call compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy)),     /**< copy method of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_COMPRFREE   ((*comprfree)),     /**< destructor of tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit)),     /**< initialize tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit)),     /**< deinitialize tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol)), /**< solving process initialization method of tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol)), /**< solving process deinitialization method of tree compression */
   SCIP_DECL_COMPREXEC   ((*comprexec)),     /**< execution method of tree compression */
   SCIP_COMPRDATA*       comprdata           /**< tree compression data */
   );

/** creates a tree compression and includes it in SCIP with its most fundamental callbacks.
 *  All non-fundamental (or optional) callbacks
 *  as, e. g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetComprCopy(), SCIPsetComprFree(),
 *  SCIPsetComprInit(), SCIPsetComprExit(), SCIPsetComprInitsol(), and SCIPsetComprExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeCompr() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeComprBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR**          compr,              /**< pointer to tree compression */
   const char*           name,               /**< name of tree compression */
   const char*           desc,               /**< description of tree compression */
   int                   priority,           /**< priority of the tree compression */
   int                   minnnodes,          /**< minimal number of nodes to call the compression */
   SCIP_DECL_COMPREXEC   ((*comprexec)),     /**< execution method of tree compression */
   SCIP_COMPRDATA*       comprdata           /**< tree compression data */
   );

/** sets copy method of tree compression */
SCIP_EXPORT
SCIP_RETCODE SCIPsetComprCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy))      /**< copy method of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of tree compression */
SCIP_EXPORT
SCIP_RETCODE SCIPsetComprFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRFREE   ((*comprfree))      /**< destructor of tree compression */
   );

/** sets initialization method of tree compression */
SCIP_EXPORT
SCIP_RETCODE SCIPsetComprInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit))      /**< initialize tree compression */
   );

/** sets deinitialization method of tree compression */
SCIP_EXPORT
SCIP_RETCODE SCIPsetComprExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit))      /**< deinitialize tree compression */
   );

/** sets solving process initialization method of tree compression */
SCIP_EXPORT
SCIP_RETCODE SCIPsetComprInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol))  /**< solving process initialization method of tree compression */
   );

/** sets solving process deinitialization method of tree compression */
SCIP_EXPORT
SCIP_RETCODE SCIPsetComprExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol))  /**< solving process deinitialization method of tree compression */
   );

/** returns the tree compression of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_COMPR* SCIPfindCompr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of tree compression */
   );

/** returns the array of currently available tree compression */
SCIP_EXPORT
SCIP_COMPR** SCIPgetComprs(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available tree compression */
SCIP_EXPORT
int SCIPgetNCompr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** set the priority of a tree compression method */
SCIP_RETCODE SCIPsetComprPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< compression */
   int                   priority            /**< new priority of the tree compression */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
