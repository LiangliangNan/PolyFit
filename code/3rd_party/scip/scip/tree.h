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

/**@file   tree.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for branch and bound tree
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TREE_H__
#define __SCIP_TREE_H__


#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/nodesel.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_cons.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_branch.h"
#include "scip/type_prop.h"
#include "scip/type_implics.h"
#include "scip/type_history.h"
#include "scip/type_conflictstore.h"
#include "scip/pub_tree.h"

#ifndef NDEBUG
#include "scip/struct_tree.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Node methods
 */

/** creates a child node of the focus node */
SCIP_RETCODE SCIPnodeCreateChild(
   SCIP_NODE**           node,               /**< pointer to node data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             nodeselprio,        /**< node selection priority of new node */
   SCIP_Real             estimate            /**< estimate for (transformed) objective value of best feasible solution in subtree */
   );

/** frees node */
SCIP_RETCODE SCIPnodeFree(
   SCIP_NODE**           node,               /**< node data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** increases the reference counter of the LP state in the fork or subroot node */
SCIP_RETCODE SCIPnodeCaptureLPIState(
   SCIP_NODE*            node,               /**< fork/subroot node */
   int                   nuses               /**< number to add to the usage counter */
   );

/** decreases the reference counter of the LP state in the fork or subroot node */
SCIP_RETCODE SCIPnodeReleaseLPIState(
   SCIP_NODE*            node,               /**< fork/subroot node */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** installs a child, a sibling, or a leaf node as the new focus node */
SCIP_RETCODE SCIPnodeFocus(
   SCIP_NODE**           node,               /**< pointer to node to focus (or NULL to remove focus); the node
                                              *   is freed, if it was cut off due to a cut off subtree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the given node can be cut off */
   SCIP_Bool             postponed,          /**< was the current focus node postponed? */
   SCIP_Bool             exitsolve           /**< are we in exitsolve stage, so we only need to loose the children */
   );

/** cuts off node and whole sub tree from branch and bound tree */
SCIP_RETCODE SCIPnodeCutoff(
   SCIP_NODE*            node,               /**< node that should be cut off */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** marks node, that propagation should be applied again the next time, a node of its subtree is focused */
void SCIPnodePropagateAgain(
   SCIP_NODE*            node,               /**< node that should be propagated again */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** marks node, that it is completely propagated in the current repropagation subtree level */
void SCIPnodeMarkPropagated(
   SCIP_NODE*            node,               /**< node that should be propagated again */
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** adds constraint locally to the node and captures it; activates constraint, if node is active;
 *  if a local constraint is added to the root node, it is automatically upgraded into a global constraint
 */
SCIP_RETCODE SCIPnodeAddCons(
   SCIP_NODE*            node,               /**< node to add constraint to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONS*            cons                /**< constraint to add */
   );

/** locally deletes constraint at the given node by disabling its separation, enforcing, and propagation capabilities
 *  at the node; captures constraint; disables constraint, if node is active
 */
SCIP_RETCODE SCIPnodeDelCons(
   SCIP_NODE*            node,               /**< node to add constraint to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   );

/** return all bound changes based on constraint propagation; stop saving the bound changes if we reach a branching
 *  decision based on a dual information
 */
void SCIPnodeGetConsProps(
   SCIP_NODE*            node,               /**< node */
   SCIP_VAR**            vars,               /**< array of variables on which constraint propagation triggers a bound change */
   SCIP_Real*            varbounds,          /**< array of bounds set by constraint propagation */
   SCIP_BOUNDTYPE*       varboundtypes,      /**< array of boundtypes set by constraint propagation */
   int*                  nconspropvars,      /**< number of variables on which constraint propagation triggers a bound change
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   conspropvarssize    /**< available slots in arrays */
   );

/** gets all bound changes applied after the first bound change based on dual information.
 *
 *  @note: currently, we can only detect bound changes based in dual information if they arise from strong branching.
 */
void SCIPnodeGetBdChgsAfterDual(
   SCIP_NODE*            node,               /**< node */
   SCIP_VAR**            vars,               /**< array of variables on which the branching has been performed in the parent node */
   SCIP_Real*            varbounds,          /**< array of bounds which the branching in the parent node set */
   SCIP_BOUNDTYPE*       varboundtypes,      /**< array of boundtypes which the branching in the parent node set */
   int                   start,              /**< first free slot in the arrays */
   int*                  nbranchvars,        /**< number of variables on which branching has been performed in the parent node
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   branchvarssize      /**< available slots in arrays */
   );

/** adds bound change with inference information to focus node, child of focus node, or probing node;
 *  if possible, adjusts bound to integral value;
 *  at most one of infercons and inferprop may be non-NULL
 */
SCIP_RETCODE SCIPnodeAddBoundinfer(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   );

/** adds bound change to focus node, or child of focus node, or probing node;
 *  if possible, adjusts bound to integral value
 */
SCIP_RETCODE SCIPnodeAddBoundchg(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_Bool             probingchange       /**< is the bound change a temporary setting due to probing? */
   );

/** adds hole with inference information to focus node, child of focus node, or probing node;
 *  if possible, adjusts bound to integral value;
 *  at most one of infercons and inferprop may be non-NULL
 */
SCIP_RETCODE SCIPnodeAddHoleinfer(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             left,               /**< left bound of open interval defining the hole (left,right) */
   SCIP_Real             right,              /**< right bound of open interval defining the hole (left,right) */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             probingchange,      /**< is the bound change a temporary setting due to probing? */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added, or NULL */
   );

/** adds hole change to focus node, or child of focus node */
SCIP_RETCODE SCIPnodeAddHolechg(
   SCIP_NODE*            node,               /**< node to add bound change to */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             left,               /**< left bound of open interval defining the hole (left,right) */
   SCIP_Real             right,              /**< right bound of open interval defining the hole (left,right) */
   SCIP_Bool             probingchange,      /**< is the bound change a temporary setting due to probing? */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added, or NULL */
   );

/** if given value is larger than the node's lower bound, sets the node's lower bound to the new value */
void SCIPnodeUpdateLowerbound(
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   );

/** updates lower bound of node using lower bound of LP */
SCIP_RETCODE SCIPnodeUpdateLowerboundLP(
   SCIP_NODE*            node,               /**< node to set lower bound for */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp                  /**< LP data */
   );

/** change the node selection priority of the given child */
void SCIPchildChgNodeselPrio(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            child,              /**< child to update the node selection priority */
   SCIP_Real             priority            /**< node selection priority value */
   );


/** sets the node's estimated bound to the new value */
void SCIPnodeSetEstimate(
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newestimate         /**< new estimated bound for the node */
   );

/** propagates implications of binary fixings at the given node triggered by the implication graph and the clique table */
SCIP_RETCODE SCIPnodePropagateImplics(
   SCIP_NODE*            node,               /**< node to propagate implications on */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   );

/** returns all bound changes based on dual information.
 *
 *  currently, this methods works only for bound changes made by strong branching on binary variables. we need this
 *  method to ensure optimality within reoptimization.
 *
 *  since the bound changes made by strong branching are stored as SCIP_BOUNDCHGTYPE_CONSINFER or SCIP_BOUNDCHGTYPE_PROPINFER
 *  with no constraint or propagator, resp., we are are interested in bound changes with these attributes.
 *
 *  all bound changes of type SCIP_BOUNDCHGTYPE_BRANCHING are stored in the beginning of the bound change array, afterwards,
 *  we can find the other two types. thus, we start the search at the end of the list and stop when reaching the first
 *  bound change of type SCIP_BOUNDCHGTYPE_BRANCHING.
 */
void SCIPnodeGetDualBoundchgs(
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR**            vars,               /**< array of variables on which the bound change is based on dual information */
   SCIP_Real*            bounds,             /**< array of bounds which are based on dual information */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array of boundtypes which are based on dual information */
   int*                  nvars,              /**< number of variables on which the bound change is based on dual information
                                              *   if this is larger than the array size, arrays should be reallocated and method
                                              *   should be called again */
   int                   varssize            /**< available slots in arrays */
   );

/** returns the number of bound changes based on dual information.
 *
 *  currently, this methods works only for bound changes made by strong branching on binary variables. we need this
 *  method to ensure optimality within reoptimization.
 *
 *  since the bound changes made by strong branching are stored as SCIP_BOUNDCHGTYPE_CONSINFER or SCIP_BOUNDCHGTYPE_PROPINFER
 *  with no constraint or propagator, resp., we are are interested in bound changes with these attributes.
 *
 *  all bound changes of type SCIP_BOUNDCHGTYPE_BRANCHING are stored in the beginning of the bound change array, afterwards,
 *  we can find the other two types. thus, we start the search at the end of the list and stop when reaching the first
 *  bound change of type SCIP_BOUNDCHGTYPE_BRANCHING.
 */
int SCIPnodeGetNDualBndchgs(
   SCIP_NODE*            node
   );

/*
 * Tree methods
 */

/** creates an initialized tree data structure */
SCIP_RETCODE SCIPtreeCreate(
   SCIP_TREE**           tree,               /**< pointer to tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting leaves in the priority queue */
   );

/** frees tree data structure */
SCIP_RETCODE SCIPtreeFree(
   SCIP_TREE**           tree,               /**< pointer to tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** clears and resets tree data structure and deletes all nodes */
SCIP_RETCODE SCIPtreeClear(
   SCIP_TREE*            tree,               /**< tree data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** creates the root node of the tree and puts it into the leaves queue */
SCIP_RETCODE SCIPtreeCreateRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** creates a temporary presolving root node of the tree and installs it as focus node */
SCIP_RETCODE SCIPtreeCreatePresolvingRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** frees the temporary presolving root and resets tree data structure */
SCIP_RETCODE SCIPtreeFreePresolvingRoot(
   SCIP_TREE*            tree,               /**< tree data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** returns the node selector associated with the given node priority queue */
SCIP_NODESEL* SCIPtreeGetNodesel(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** sets the node selector used for sorting the nodes in the priority queue, and resorts the queue if necessary */
SCIP_RETCODE SCIPtreeSetNodesel(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   );

/** cuts off nodes with lower bound not better than given upper bound */
SCIP_RETCODE SCIPtreeCutoff(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   );

/** constructs the LP relaxation of the focus node */
SCIP_RETCODE SCIPtreeLoadLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool*            initroot            /**< pointer to store whether the root LP relaxation has to be initialized */
   );

/** loads LP state for fork/subroot of the focus node */
SCIP_RETCODE SCIPtreeLoadLPState(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** calculates the node selection priority for moving the given variable's LP value to the given target value;
 *  this node selection priority can be given to the SCIPcreateChild() call
 */
SCIP_Real SCIPtreeCalcNodeselPriority(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_BRANCHDIR        branchdir,          /**< type of branching that was performed: upwards, downwards, or fixed 
                                              * fixed should only be used, when both bounds changed 
                                              */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   );

/** calculates an estimate for the objective of the best feasible solution contained in the subtree after applying the given 
 *  branching; this estimate can be given to the SCIPcreateChild() call
 */
SCIP_Real SCIPtreeCalcChildEstimate(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   );

/** branches on a variable x
 *  if x is a continuous variable, then two child nodes will be created
 *  (x <= x', x >= x')
 *  but if the bounds of x are such that their relative difference is smaller than epsilon,
 *  the variable is fixed to val (if not SCIP_INVALID) or a well chosen alternative in the current node,
 *  i.e., no children are created
 *  if x is not a continuous variable, then:
 *  if solution value x' is fractional, two child nodes will be created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if solution value is integral, the x' is equal to lower or upper bound of the branching
 *  variable and the bounds of x are finite, then two child nodes will be created
 *  (x <= x", x >= x"+1 with x" = floor((lb + ub)/2)),
 *  otherwise (up to) three child nodes will be created
 *  (x <= x'-1, x == x', x >= x'+1)
 *  if solution value is equal to one of the bounds and the other bound is infinite, only two child nodes
 *  will be created (the third one would be infeasible anyway)
 */
SCIP_RETCODE SCIPtreeBranchVar(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on or SCIP_INVALID for branching on current LP/pseudo solution. A branching value is required for branching on continuous variables */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   );

/** branches a variable x using the given domain hole; two child nodes will be created (x <= left, x >= right) */
SCIP_RETCODE SCIPtreeBranchVarHole(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             left,               /**< left side of the domain hole */
   SCIP_Real             right,              /**< right side of the domain hole */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   );

/** n-ary branching on a variable x
 * Branches on variable x such that up to n/2 children are created on each side of the usual branching value.
 * The branching value is selected as in SCIPtreeBranchVar().
 * If n is 2 or the variables local domain is too small for a branching into n pieces, SCIPtreeBranchVar() is called.
 * The parameters minwidth and widthfactor determine the domain width of the branching variable in the child nodes.
 * If n is odd, one child with domain width 'width' and having the branching value in the middle is created.
 * Otherwise, two children with domain width 'width' and being left and right of the branching value are created.
 * Next further nodes to the left and right are created, where width is multiplied by widthfactor with increasing distance from the first nodes.
 * The initial width is calculated such that n/2 nodes are created to the left and to the right of the branching value.
 * If this value is below minwidth, the initial width is set to minwidth, which may result in creating less than n nodes.
 *
 * Giving a large value for widthfactor results in creating children with small domain when close to the branching value
 * and large domain when closer to the current variable bounds. That is, setting widthfactor to a very large value and n to 3
 * results in a ternary branching where the branching variable is mostly fixed in the middle child.
 * Setting widthfactor to 1.0 results in children where the branching variable always has the same domain width
 * (except for one child if the branching value is not in the middle).
 */
SCIP_RETCODE SCIPtreeBranchVarNary(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on or SCIP_INVALID for branching on current LP/pseudo solution.
                                              *   A branching value is required for branching on continuous variables */
   int                   n,                  /**< attempted number of children to be created, must be >= 2 */
   SCIP_Real             minwidth,           /**< minimal domain width in children */
   SCIP_Real             widthfactor,        /**< multiplier for children domain width with increasing distance from val, must be >= 1.0 */
   int*                  nchildren           /**< buffer to store number of created children, or NULL */
   );

/** adds a diving bound change to the tree together with the information if this is a bound change
 *  for the preferred direction or not
 */
SCIP_RETCODE SCIPtreeAddDiveBoundChange(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_VAR*             var,                /**< variable to apply the bound change to */
   SCIP_BRANCHDIR        dir,                /**< direction of the bound change */
   SCIP_Real             value,              /**< value to adjust this variable bound to */
   SCIP_Bool             preferred           /**< is this a bound change for the preferred child? */
   );

/** get the dive bound change data for the preferred or the alternative direction */
void SCIPtreeGetDiveBoundChangeData(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_VAR***           variables,          /**< pointer to store variables for the specified direction */
   SCIP_BRANCHDIR**      directions,         /**< pointer to store the branching directions */
   SCIP_Real**           values,             /**< pointer to store bound change values */
   int*                  ndivebdchgs,        /**< pointer to store the number of dive bound changes */
   SCIP_Bool             preferred           /**< should the dive bound changes for the preferred child be output? */
   );

/** clear the tree dive bound change data structure */
void SCIPtreeClearDiveBoundChanges(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** switches to probing mode and creates a probing root */
SCIP_RETCODE SCIPtreeStartProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_Bool             strongbranching     /**< is the probing mode used for strongbranching? */
   );

/** creates a new probing child node in the probing path */
SCIP_RETCODE SCIPtreeCreateProbingNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** sets the LP state for the current probing node
 *
 *  @note state and norms are stored at the node and later released by SCIP; therefore, the pointers are set
 *        to NULL by the method
 *
 *  @note the pointers to state and norms must not be NULL; however, they may point to a NULL pointer if the
 *        respective information should not be set
 */
SCIP_RETCODE SCIPtreeSetProbingLPState(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_LPISTATE**       lpistate,           /**< pointer to LP state information (like basis information) */
   SCIP_LPINORMS**       lpinorms,           /**< pointer to LP pricing norms information */
   SCIP_Bool             primalfeas,         /**< primal feasibility when LP state information was stored */
   SCIP_Bool             dualfeas            /**< dual feasibility when LP state information was stored */
   );

/** loads the LP state for the current probing node */
SCIP_RETCODE SCIPtreeLoadProbingLPState(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** marks the probing node to have a solved LP relaxation */
SCIP_RETCODE SCIPtreeMarkProbingNodeHasLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** undoes all changes to the problem applied in probing up to the given probing depth;
 *  the changes of the probing node of the given probing depth are the last ones that remain active;
 *  changes that were applied before calling SCIPtreeCreateProbingNode() cannot be undone
 */
SCIP_RETCODE SCIPtreeBacktrackProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_PRIMAL*          primal,             /**< primal data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   );

/** switches back from probing to normal operation mode, frees all nodes on the probing path, restores bounds of all
 *  variables and restores active constraints arrays of focus node
 */
SCIP_RETCODE SCIPtreeEndProbing(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PRIMAL*          primal,             /**< primal LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** stores relaxation solution before diving or probing */
SCIP_RETCODE SCIPtreeStoreRelaxSol(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PROB*            transprob           /**< transformed problem after presolve */
   );

/** restores relaxation solution after diving or probing */
SCIP_RETCODE SCIPtreeRestoreRelaxSol(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_PROB*            transprob           /**< transformed problem after presolve */
   );


/** gets number of children of the focus node  */
int SCIPtreeGetNChildren(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets number of siblings of the focus node  */
int SCIPtreeGetNSiblings(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets number of leaves in the tree (excluding children and siblings of focus nodes) */
int SCIPtreeGetNLeaves(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets number of open nodes in the tree (children + siblings + leaves) */
int SCIPtreeGetNNodes(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns whether the active path goes completely down to the focus node */
SCIP_Bool SCIPtreeIsPathComplete(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns whether the current node is a temporary probing node */
SCIP_Bool SCIPtreeProbing(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns the temporary probing root node, or NULL if the we are not in probing mode */
SCIP_NODE* SCIPtreeGetProbingRoot(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns the current probing depth, i.e. the number of probing sub nodes existing in the probing path */
int SCIPtreeGetProbingDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets focus node of the tree */
SCIP_NODE* SCIPtreeGetFocusNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets depth of focus node in the tree, or -1 if no focus node exists */
int SCIPtreeGetFocusDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns, whether the LP was or is to be solved in the focus node */
SCIP_Bool SCIPtreeHasFocusNodeLP(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** sets mark to solve or to ignore the LP while processing the focus node */
void SCIPtreeSetFocusNodeLP(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             solvelp             /**< should the LP be solved in focus node? */
   );

/** returns whether the LP of the focus node is already constructed */
SCIP_Bool SCIPtreeIsFocusNodeLPConstructed(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns whether the focus node is already solved and only propagated again */
SCIP_Bool SCIPtreeInRepropagation(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets current node of the tree, i.e. the last node in the active path, or NULL if no current node exists */
SCIP_NODE* SCIPtreeGetCurrentNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets depth of current node in the tree, i.e. the length of the active path minus 1, or -1 if no current node exists */
int SCIPtreeGetCurrentDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns, whether the LP was or is to be solved in the current node */
SCIP_Bool SCIPtreeHasCurrentNodeLP(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns the depth of the effective root node (i.e. the first depth level of a node with at least two children) */
int SCIPtreeGetEffectiveRootDepth(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets the root node of the tree */
SCIP_NODE* SCIPtreeGetRootNode(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** returns whether we are in probing and the objective value of at least one column was changed */
SCIP_Bool SCIPtreeProbingObjChanged(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** marks the current probing node to have a changed objective function */
void SCIPtreeMarkProbingObjChanged(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPtreeGetNLeaves(tree)        SCIPnodepqLen((tree)->leaves)
#define SCIPtreeGetNChildren(tree)      ((tree)->nchildren)
#define SCIPtreeGetNSiblings(tree)      ((tree)->nsiblings)
#define SCIPtreeGetNNodes(tree)         \
   (SCIPtreeGetNChildren(tree) + SCIPtreeGetNSiblings(tree) + SCIPtreeGetNLeaves(tree))
#define SCIPtreeIsPathComplete(tree)    ((tree)->focusnode == NULL || (tree)->focusnode->depth < (tree)->pathlen)
#define SCIPtreeProbing(tree)           ((tree)->probingroot != NULL)
#define SCIPtreeGetProbingRoot(tree)    (tree)->probingroot
#define SCIPtreeGetProbingDepth(tree)   (SCIPtreeGetCurrentDepth(tree) - SCIPnodeGetDepth((tree)->probingroot))
#define SCIPtreeGetFocusNode(tree)      (tree)->focusnode
#define SCIPtreeGetFocusDepth(tree)     ((tree)->focusnode != NULL ? (int)(tree)->focusnode->depth : -1)
#define SCIPtreeHasFocusNodeLP(tree)    (tree)->focusnodehaslp
#define SCIPtreeSetFocusNodeLP(tree,solvelp)  ((tree)->focusnodehaslp = solvelp)
#define SCIPtreeIsFocusNodeLPConstructed(tree) (tree)->focuslpconstructed
#define SCIPtreeInRepropagation(tree)   ((tree)->focusnode != NULL \
      && SCIPnodeGetType((tree)->focusnode) == SCIP_NODETYPE_REFOCUSNODE)
#define SCIPtreeGetCurrentNode(tree)    ((tree)->pathlen > 0 ? (tree)->path[(tree)->pathlen-1] : NULL)
#define SCIPtreeGetCurrentDepth(tree)   ((tree)->pathlen-1)
#define SCIPtreeHasCurrentNodeLP(tree)  (SCIPtreeProbing(tree) ? (tree)->probingnodehaslp : SCIPtreeHasFocusNodeLP(tree))
#define SCIPtreeGetEffectiveRootDepth(tree) ((tree)->effectiverootdepth)
#define SCIPtreeGetRootNode(tree)       ((tree)->root)
#define SCIPtreeProbingObjChanged(tree) ((tree)->probingobjchanged)
#define SCIPtreeMarkProbingObjChanged(tree) ((tree)->probingobjchanged = TRUE)

#endif


/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPtreeGetPrioChild(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPtreeGetPrioSibling(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets the best child of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestChild(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the best sibling of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestSibling(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the best leaf from the node queue w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestLeaf(
   SCIP_TREE*            tree                /**< branch and bound tree */
   );

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
SCIP_NODE* SCIPtreeGetBestNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the minimal lower bound of all nodes in the tree */
SCIP_Real SCIPtreeGetLowerbound(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the node with minimal lower bound of all nodes in the tree (child, sibling, or leaf) */
SCIP_NODE* SCIPtreeGetLowerboundNode(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the average lower bound of all nodes in the tree */
SCIP_Real SCIPtreeGetAvgLowerbound(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Real             cutoffbound         /**< global cutoff bound */
   );

/** query if focus node was already branched on */
SCIP_Bool SCIPtreeWasNodeLastBranchParent(
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            node                /**< tree node, or NULL to check focus node */
   );

#ifdef __cplusplus
}
#endif

#endif
