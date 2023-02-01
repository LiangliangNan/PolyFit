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
EXTERN
SCIP_RETCODE SCIPvisualCreate(
   SCIP_VISUAL**         visual,             /**< pointer to store the visualization information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** frees visualization data structure */
EXTERN
void SCIPvisualFree(
   SCIP_VISUAL**         visual              /**< pointer to store the visualization information */
   );

/** initializes visualization information and creates a file for visualization output */
EXTERN
SCIP_RETCODE SCIPvisualInit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** closes the visualization output file */
EXTERN
void SCIPvisualExit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** creates a new node entry in the visualization output file */
EXTERN
SCIP_RETCODE SCIPvisualNewChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** updates a node entry in the visualization output file */
extern
SCIP_RETCODE SCIPvisualUpdateChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   );

/** marks node as solved in visualization output file */
EXTERN
void SCIPvisualSolvedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was solved */
   );

/** changes the color of the node to the color of cutoff nodes */
EXTERN
void SCIPvisualCutoffNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node, that was cut off */
   SCIP_Bool             infeasible          /**< whether the node is infeasible (otherwise exceeded the cutoff bound) */
   );

/** changes the color of the node to the color of nodes where a conflict constraint was found */
EXTERN
void SCIPvisualFoundConflict(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, where the conflict was found */
   );

/** changes the color of the node to the color of nodes that were marked to be repropagated */
EXTERN
void SCIPvisualMarkedRepropagateNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was marked to be repropagated */
   );

/** changes the color of the node to the color of repropagated nodes */
EXTERN
void SCIPvisualRepropagatedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was repropagated */
   );

/** changes the color of the node to the color of nodes with a primal solution */
EXTERN
void SCIPvisualFoundSolution(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node where the solution was found, or NULL */
   SCIP_Bool             bettersol,          /**< the solution was better than the previous ones */
   SCIP_SOL*             sol                 /**< solution that has been found */
   );

/** outputs a new global lower bound to the visualization output file */
EXTERN
void SCIPvisualLowerbound(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             lowerbound          /**< new lower bound */
   );

/** outputs a new global upper bound to the visualization output file */
EXTERN
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
