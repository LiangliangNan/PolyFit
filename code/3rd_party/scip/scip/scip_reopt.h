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

/**@file   scip_reopt.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for reoptimization
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

#ifndef __SCIP_SCIP_REOPT_H__
#define __SCIP_SCIP_REOPT_H__


#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_reopt.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicReoptimizationMethods
 *
 * @{
 */

/** return the ids of child nodes stored in the reoptimization tree
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetReoptChildIDs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         ids,                /**< array to store the ids of child nodes */
   int                   mem,                /**< allocated memory */
   int*                  nids                /**< number of child nodes */
   );

/** return the ids of all leave nodes store in the reoptimization tree induced by the given node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetReoptLeaveIDs(
   SCIP*                 scip,               /**< SCIP data strcuture */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         ids,                /**< array of ids */
   int                   mem,                /**< allocated memory */
   int*                  nids                /**< number of child nodes */
   );

/** returns the number of nodes in the reoptimization tree induced by @p node; if @p node == NULL, the method
 *  returns the number of nodes of the whole reoptimization tree.
 */
SCIP_EXPORT
int SCIPgetNReoptnodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** returns the number of leave nodes of the subtree induced by @p node; if @p node == NULL, the method
 *  returns the number of leaf nodes of the whole reoptimization tree.
 */
SCIP_EXPORT
int SCIPgetNReoptLeaves(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** gets the node of the reoptimization tree corresponding to the unique @p id */
SCIP_REOPTNODE* SCIPgetReoptnode(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          id                  /**< unique id */
   );

/** add a variable bound change to a given reoptnode
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddReoptnodeBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR*             var,                /**< variable pointer */
   SCIP_Real             bound,              /**< variable bound to add */
   SCIP_BOUNDTYPE        boundtype           /**< bound type of the variable value */
   );

/** set the @p representation as the new search frontier
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsetReoptCompression(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representation,     /**< array of representatives */
   int                   nrepresentatives,   /**< number of representatives */
   SCIP_Bool*            success             /**< pointer to store the result */
   );

/** add stored constraint to a reoptimization node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddReoptnodeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of variable bounds */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of variable boundtypes */
   SCIP_Real             lhs,                /**< lhs of the constraint */
   SCIP_Real             rhs,                /**< rhs of the constraint */
   int                   nvars,              /**< number of variables */
   REOPT_CONSTYPE        constype,           /**< type of the constraint */
   SCIP_Bool             linear              /**< the given constraint has a linear representation */
   );

/** return the branching path stored in the reoptree at ID id */
SCIP_EXPORT
void SCIPgetReoptnodePath(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of variable bounds */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of bound types */
   int                   mem,                /**< allocated memory */
   int*                  nvars,              /**< number of variables */
   int*                  nafterdualvars      /**< number of variables directly after the first based on dual information */
   );

/** initialize a set of empty reoptimization nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPinitRepresentation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives    /**< number of representatives */
   );

/** reset a set of initialized reoptimization nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_RETCODE SCIPresetRepresentation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives    /**< number of representatives */
   );

/** free a set of initialized reoptimization nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeRepresentation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives    /**< number of representatives */
   );

/** reactivate the given @p reoptnode and split them into several nodes if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPapplyReopt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node to reactivate */
   unsigned int          id,                 /**< unique id of the reoptimization node */
   SCIP_Real             estimate,           /**< estimate of the child nodes that should be created */
   SCIP_NODE**           childnodes,         /**< array to store the created child nodes */
   int*                  ncreatedchilds,     /**< pointer to store number of created child nodes */
   int*                  naddedconss,        /**< pointer to store number of generated constraints */
   int                   childnodessize,     /**< available size of childnodes array */
   SCIP_Bool*            success             /**< pointer store the result*/
   );

/** remove the stored information about bound changes based in dual information
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_EXPORT
SCIP_RETCODE SCIPresetReoptnodeDualcons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** splits the root into several nodes and moves the child nodes of the root to one of the created nodes
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsplitReoptRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  ncreatedchilds,     /**< pointer to store the number of created nodes */
   int*                  naddedconss         /**< pointer to store the number added constraints */
   );

/** returns if a node should be reoptimized */
SCIP_EXPORT
SCIP_Bool SCIPreoptimizeNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** deletes the given reoptimization node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdeleteReoptnode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPTNODE**      reoptnode           /**< node of the reoptimization tree */
   );

/** return the similarity between two objective functions */
SCIP_EXPORT
SCIP_Real SCIPgetReoptSimilarity(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   run1,               /**< number of run */
   int                   run2                /**< number of run */
   );

/** check the changes of the variable coefficient in the objective function */
SCIP_EXPORT
void SCIPgetVarCoefChg(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   varidx,             /**< index of variable */
   SCIP_Bool*            negated,            /**< coefficient changed the sign */
   SCIP_Bool*            entering,           /**< coefficient gets non-zero coefficient */
   SCIP_Bool*            leaving             /**< coefficient gets zero coefficient */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
