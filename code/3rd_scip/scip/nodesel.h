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

/**@file   nodesel.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for node selectors and node priority queues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NODESEL_H__
#define __SCIP_NODESEL_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_event.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/pub_nodesel.h"

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * node priority queue methods
 */

/** creates node priority queue */
SCIP_RETCODE SCIPnodepqCreate(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   );

/** frees node priority queue, but not the data nodes themselves */
void SCIPnodepqDestroy(
   SCIP_NODEPQ**         nodepq              /**< pointer to a node priority queue */
   );

/** frees node priority queue and all nodes in the queue */
SCIP_RETCODE SCIPnodepqFree(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** deletes all nodes in the node priority queue */
SCIP_RETCODE SCIPnodepqClear(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** returns the node selector associated with the given node priority queue */
SCIP_NODESEL* SCIPnodepqGetNodesel(
   SCIP_NODEPQ*          nodepq              /**< node priority queue */
   );

/** sets the node selector used for sorting the nodes in the queue, and resorts the queue if necessary */
SCIP_RETCODE SCIPnodepqSetNodesel(
   SCIP_NODEPQ**         nodepq,             /**< pointer to a node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector to use for sorting the nodes in the queue */
   );

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
int SCIPnodepqCompare(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node1,              /**< first node to compare */
   SCIP_NODE*            node2               /**< second node to compare */
   );

/** inserts node into node priority queue */
SCIP_RETCODE SCIPnodepqInsert(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to be inserted */
   );

/** removes node from the node priority queue */
SCIP_RETCODE SCIPnodepqRemove(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node to remove */
   );

/** returns the best node of the queue without removing it */
SCIP_NODE* SCIPnodepqFirst(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   );

/** returns the nodes array of the queue */
SCIP_NODE** SCIPnodepqNodes(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   );

/** returns the number of nodes stored in the node priority queue */
int SCIPnodepqLen(
   const SCIP_NODEPQ*    nodepq              /**< node priority queue */
   );

/** gets the minimal lower bound of all nodes in the queue */
SCIP_Real SCIPnodepqGetLowerbound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the node with minimal lower bound of all nodes in the queue */
SCIP_NODE* SCIPnodepqGetLowerboundNode(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets the sum of lower bounds of all nodes in the queue */
SCIP_Real SCIPnodepqGetLowerboundSum(
   SCIP_NODEPQ*          nodepq              /**< node priority queue */
   );

/** free all nodes from the queue that are cut off by the given upper bound */
SCIP_RETCODE SCIPnodepqBound(
   SCIP_NODEPQ*          nodepq,             /**< node priority queue */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             cutoffbound         /**< cutoff bound: all nodes with lowerbound >= cutoffbound are cut off */
   );




/*
 * node selector methods
 */

/** copies the given node selector to a new scip */
SCIP_RETCODE SCIPnodeselCopyInclude(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a node selector */
SCIP_RETCODE SCIPnodeselCreate(
   SCIP_NODESEL**        nodesel,            /**< pointer to store node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
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

/** frees memory of node selector */
SCIP_RETCODE SCIPnodeselFree(
   SCIP_NODESEL**        nodesel,            /**< pointer to node selector data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes node selector */
SCIP_RETCODE SCIPnodeselInit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** deinitializes node selector */
SCIP_RETCODE SCIPnodeselExit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs node selector that the branch and bound process is being started */
SCIP_RETCODE SCIPnodeselInitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs node selector that the branch and bound process data is being freed */
SCIP_RETCODE SCIPnodeselExitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** select next node to be processed */
SCIP_RETCODE SCIPnodeselSelect(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE**           selnode             /**< pointer to store node to be processed next */
   );

/** compares two nodes; returns -1/0/+1 if node1 better/equal/worse than node2 */
int SCIPnodeselCompare(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node1,              /**< first node to compare */
   SCIP_NODE*            node2               /**< second node to compare */
   );

/** sets priority of node selector in standard mode */
void SCIPnodeselSetStdPriority(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the node selector */
   );

/** sets priority of node selector in memory saving mode */
void SCIPnodeselSetMemsavePriority(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the node selector */
   );

/** sets copy method of node selector */
void SCIPnodeselSetCopy(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy))    /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of node selector */
void SCIPnodeselSetFree(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELFREE ((*nodeselfree))    /**< destructor of node selector */
   );

/** sets initialization method of node selector */
void SCIPnodeselSetInit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit))    /**< initialize node selector */
   );

/** sets deinitialization method of node selector */
void SCIPnodeselSetExit(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit))    /**< deinitialize node selector */
   );

/** sets solving process initialization method of node selector */
void SCIPnodeselSetInitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELINITSOL ((*nodeselinitsol))/**< solving process initialization method of node selector */
   );

/** sets solving process deinitialization method of node selector */
void SCIPnodeselSetExitsol(
   SCIP_NODESEL*         nodesel,            /**< node selector */
   SCIP_DECL_NODESELEXITSOL ((*nodeselexitsol))/**< solving process deinitialization method of node selector */
   );

/** enables or disables all clocks of \p nodesel, depending on the value of the flag */
void SCIPnodeselEnableOrDisableClocks(
   SCIP_NODESEL*         nodesel,            /**< the node selector for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the node selector be enabled? */
   );

#ifdef __cplusplus
}
#endif

#endif
