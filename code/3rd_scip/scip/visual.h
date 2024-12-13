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

/**@file   visual.h
 * @ingroup INTERNALAPI
 * @brief  methods for creating output for visualization tools (VBC, BAK)
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VISUAL_H__
#define __SCIP_VISUAL_H__


#include "scip/def.h"
#include "scip/type_set.h"
#include "scip/type_sol.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_visual.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates visualization data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPvisualCreate(
   SCIP_VISUAL**         visual,             /**< pointer to store the visualization information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** frees visualization data structure */
SCIP_EXPORT
void SCIPvisualFree(
   SCIP_VISUAL**         visual              /**< pointer to store the visualization information */
   );

/** initializes visualization information and creates a file for visualization output */
SCIP_EXPORT
SCIP_RETCODE SCIPvisualInit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** closes the visualization output file */
SCIP_EXPORT
void SCIPvisualExit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** creates a new node entry in the visualization output file */
SCIP_EXPORT
SCIP_RETCODE SCIPvisualNewChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** updates a node entry in the visualization output file */
SCIP_RETCODE SCIPvisualUpdateChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** marks node as solved in visualization output file */
SCIP_EXPORT
void SCIPvisualSolvedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was solved */
   );

/** changes the color of the node to the color of cutoff nodes */
SCIP_EXPORT
void SCIPvisualCutoffNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node, that was cut off */
   SCIP_Bool             infeasible          /**< whether the node is infeasible (otherwise exceeded the cutoff bound) */
   );

/** changes the color of the node to the color of nodes where a conflict constraint was found */
SCIP_EXPORT
void SCIPvisualFoundConflict(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, where the conflict was found */
   );

/** changes the color of the node to the color of nodes that were marked to be repropagated */
SCIP_EXPORT
void SCIPvisualMarkedRepropagateNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was marked to be repropagated */
   );

/** changes the color of the node to the color of repropagated nodes */
SCIP_EXPORT
void SCIPvisualRepropagatedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was repropagated */
   );

/** changes the color of the node to the color of nodes with a primal solution */
SCIP_EXPORT
void SCIPvisualFoundSolution(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node where the solution was found, or NULL */
   SCIP_Bool             bettersol,          /**< the solution was better than the previous ones */
   SCIP_SOL*             sol                 /**< solution that has been found */
   );

/** outputs a new global lower bound to the visualization output file */
SCIP_EXPORT
void SCIPvisualLowerbound(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             lowerbound          /**< new lower bound */
   );

/** outputs a new global upper bound to the visualization output file */
SCIP_EXPORT
void SCIPvisualUpperbound(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             upperbound          /**< new upper bound */
   );

#ifdef __cplusplus
}
#endif

#endif
