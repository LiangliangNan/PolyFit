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

/**@file   scip_nodesel.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for node selector plugins
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

#ifndef __SCIP_SCIP_NODESEL_H__
#define __SCIP_SCIP_NODESEL_H__


#include "scip/def.h"
#include "scip/type_nodesel.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNodeSelectorMethods
 *
 * @{
 */

/** creates a node selector and includes it in SCIP.
 *
 *  @note method has all node selector callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeNodeselBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy)),   /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   SCIP_DECL_NODESELINITSOL((*nodeselinitsol)),/**< solving process initialization method of node selector */
   SCIP_DECL_NODESELEXITSOL((*nodeselexitsol)),/**< solving process deinitialization method of node selector */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   );

/** Creates a node selector and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetNodeselCopy(), SCIPsetNodeselFree(),
 *  SCIPsetNodeselInit(), SCIPsetNodeselExit(), SCIPsetNodeselInitsol(), and SCIPsetNodeselExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeNodesel() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNodeselBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL**        nodesel,            /**< reference to a node selector, or NULL */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   );

/** sets copy method of node selector */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy))    /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of node selector */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELFREE ((*nodeselfree))    /**< destructor of node selector */
   );

/** sets initialization method of node selector */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit))    /**< initialize node selector */
   );

/** sets deinitialization method of node selector */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit))    /**< deinitialize node selector */
   );

/** sets solving process initialization method of node selector */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINITSOL ((*nodeselinitsol))/**< solving process initialization method of node selector */
   );

/** sets solving process deinitialization method of node selector */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXITSOL ((*nodeselexitsol))/**< solving process deinitialization method of node selector */
   );

/** returns the node selector of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_NODESEL* SCIPfindNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of node selector */
   );

/** returns the array of currently available node selectors */
SCIP_EXPORT
SCIP_NODESEL** SCIPgetNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available node selectors */
SCIP_EXPORT
int SCIPgetNNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of a node selector in standard mode */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselStdPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new standard priority of the node selector */
   );

/** sets the priority of a node selector in memory saving mode */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNodeselMemsavePriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new memory saving priority of the node selector */
   );

/** returns the currently used node selector */
SCIP_EXPORT
SCIP_NODESEL* SCIPgetNodesel(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
