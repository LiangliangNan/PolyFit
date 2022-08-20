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

/**@file   type_tree.h
 * @brief  type definitions for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_TREE_H__
#define __SCIP_TYPE_TREE_H__

#ifdef __cplusplus
extern "C" {
#endif

enum SCIP_NodeType
{
   SCIP_NODETYPE_FOCUSNODE      =  0,   /**< the focus node, whose data is stored in the tree data structure */
   SCIP_NODETYPE_PROBINGNODE    =  1,   /**< temporary child node of the focus or refocused node used for probing */
   SCIP_NODETYPE_SIBLING        =  2,   /**< unsolved sibling of the focus node */
   SCIP_NODETYPE_CHILD          =  3,   /**< unsolved child of the focus node */
   SCIP_NODETYPE_LEAF           =  4,   /**< unsolved leaf of the tree, stored in the tree's queue */
   SCIP_NODETYPE_DEADEND        =  5,   /**< temporary type of focus node, if it was solved completely */
   SCIP_NODETYPE_JUNCTION       =  6,   /**< fork without LP solution */
   SCIP_NODETYPE_PSEUDOFORK     =  7,   /**< fork without LP solution and added rows and columns */
   SCIP_NODETYPE_FORK           =  8,   /**< fork with solved LP and added rows and columns */
   SCIP_NODETYPE_SUBROOT        =  9,   /**< fork with solved LP and arbitrarily changed rows and columns */
   SCIP_NODETYPE_REFOCUSNODE    = 10    /**< junction, fork, or subroot that was refocused for domain propagation */
};
typedef enum SCIP_NodeType SCIP_NODETYPE;         /**< type of node */

typedef struct SCIP_Probingnode SCIP_PROBINGNODE; /**< data for probing nodes */
typedef struct SCIP_Sibling SCIP_SIBLING;         /**< data for sibling nodes */
typedef struct SCIP_Child SCIP_CHILD;             /**< data for child nodes */
typedef struct SCIP_Leaf SCIP_LEAF;               /**< data for leaf nodes */
typedef struct SCIP_Junction SCIP_JUNCTION;       /**< data for junction nodes */
typedef struct SCIP_Pseudofork SCIP_PSEUDOFORK;   /**< data for pseudo fork nodes */
typedef struct SCIP_Fork SCIP_FORK;               /**< data for fork nodes */
typedef struct SCIP_Subroot SCIP_SUBROOT;         /**< data for subroot nodes */
typedef struct SCIP_Node SCIP_NODE;               /**< node data structure */
typedef struct SCIP_PendingBdchg SCIP_PENDINGBDCHG; /**< bound change information for pending bound changes */
typedef struct SCIP_Tree SCIP_TREE;               /**< branch and bound tree */

#ifdef __cplusplus
}
#endif

#endif
