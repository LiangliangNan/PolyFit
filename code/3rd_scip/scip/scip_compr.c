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

/**@file   scip_compr.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for compression plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/compr.h"
#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/scip_compr.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates a tree compression and includes it in SCIP.
 *
 *  @note method has all compression callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeComprBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
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
   )
{
   SCIP_COMPR* compr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeCompr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether compression is already present */
   if( SCIPfindCompr(scip, name) != NULL )
   {
      SCIPerrorMessage("compression <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPcomprCreate(&compr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority, minnnodes,
         comprcopy, comprfree, comprinit, comprexit, comprinitsol, comprexitsol, comprexec, comprdata) );

   SCIP_CALL( SCIPsetIncludeCompr(scip->set, compr) );

   return SCIP_OKAY;
}

/** creates a tree compression and includes it in SCIP with its most fundamental callbacks.
 *  All non-fundamental (or optional) callbacks
 *  as, e. g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetComprCopy(), SCIPsetComprFree(),
 *  SCIPsetComprInit(), SCIPsetComprExit(), SCIPsetComprInitsol(), and SCIPsetComprExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeCompr() instead
 */
SCIP_RETCODE SCIPincludeComprBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR**          compr,              /**< pointer to tree compression */
   const char*           name,               /**< name of tree compression */
   const char*           desc,               /**< description of tree compression */
   int                   priority,           /**< priority of the tree compression */
   int                   minnnodes,          /**< minimal number of nodes to call the compression */
   SCIP_DECL_COMPREXEC   ((*comprexec)),     /**< execution method of tree compression */
   SCIP_COMPRDATA*       comprdata           /**< tree compression data */
   )
{
   SCIP_COMPR* comprptr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeComprBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether heuristic is already present */
   if( SCIPfindCompr(scip, name) != NULL )
   {
      SCIPerrorMessage("tree compression <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPcomprCreate(&comprptr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         minnnodes, NULL, NULL, NULL, NULL, NULL, NULL, comprexec, comprdata) );

   assert(comprptr != NULL);

   SCIP_CALL( SCIPsetIncludeCompr(scip->set, comprptr) );

   if( compr != NULL )
      *compr = comprptr;

   return SCIP_OKAY;
}

/* new callback/method setter methods */

/** sets copy method of tree compression */
SCIP_RETCODE SCIPsetComprCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRCOPY   ((*comprcopy))      /**< copy method of tree compression or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetComprCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(compr != NULL);

   SCIPcomprSetCopy(compr, comprcopy);

   return SCIP_OKAY;
}

/** sets destructor method of tree compression */
SCIP_RETCODE SCIPsetComprFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRFREE   ((*comprfree))      /**< destructor of tree compression */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetComprFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(compr != NULL);

   SCIPcomprSetFree(compr, comprfree);

   return SCIP_OKAY;
}

/** sets initialization method of tree compression */
SCIP_RETCODE SCIPsetComprInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINIT   ((*comprinit))      /**< initialize tree compression */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetComprInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(compr != NULL);

   SCIPcomprSetInit(compr, comprinit);

   return SCIP_OKAY;
}

/** sets deinitialization method of tree compression */
SCIP_RETCODE SCIPsetComprExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXIT   ((*comprexit))      /**< deinitialize tree compression */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetComprExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(compr != NULL);

   SCIPcomprSetExit(compr, comprexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of tree compression */
SCIP_RETCODE SCIPsetComprInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPRINITSOL ((*comprinitsol))  /**< solving process initialization method of tree compression */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetComprInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(compr != NULL);

   SCIPcomprSetInitsol(compr, comprinitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of tree compression */
SCIP_RETCODE SCIPsetComprExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< tree compression */
   SCIP_DECL_COMPREXITSOL ((*comprexitsol))  /**< solving process deinitialization method of tree compression */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetComprExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(compr != NULL);

   SCIPcomprSetExitsol(compr, comprexitsol);

   return SCIP_OKAY;
}

/** returns the tree compression of the given name, or NULL if not existing */
SCIP_COMPR* SCIPfindCompr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of tree compression */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindCompr(scip->set, name);
}

/** returns the array of currently available tree compression */
SCIP_COMPR** SCIPgetComprs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortComprs(scip->set);

   return scip->set->comprs;
}

/** returns the number of currently available tree compression */
int SCIPgetNCompr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->ncomprs;
}

/** set the priority of a tree compression method */
SCIP_RETCODE SCIPsetComprPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COMPR*           compr,              /**< compression */
   int                   priority            /**< new priority of the tree compression */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPcomprSetPriority(compr, scip->set, priority);

   return SCIP_OKAY;
}
