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

/**@file   struct_tree.h
 * @ingroup INTERNALAPI
 * @brief  data structures for branch and bound tree
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_TREE_H__
#define __SCIP_STRUCT_TREE_H__


#include "lpi/type_lpi.h"
#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_history.h"
#include "scip/type_lp.h"
#include "scip/type_nodesel.h"
#include "scip/type_prop.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"


#ifdef __cplusplus
extern "C" {
#endif

/** probing node, possibly with solved LP, where bounds and constraints have been changed,
 *  and rows and columns might have been added
 */
struct SCIP_Probingnode
{
   SCIP_LPISTATE*        lpistate;           /**< LP state information */
   SCIP_LPINORMS*        lpinorms;           /**< LP pricing norms information */
   int                   ninitialcols;       /**< number of LP columns before the node was processed */
   int                   ninitialrows;       /**< number of LP rows before the node was processed */
   int                   ncols;              /**< total number of columns of this node's LP */
   int                   nrows;              /**< total number of rows of this node's LP */
   SCIP_VAR**            origobjvars;        /**< variables whose objective function coefficients have changed */
   SCIP_Real*            origobjvals;        /**< original objective function coefficients */
   int                   nchgdobjs;          /**< number of changed objective coefficients */
   SCIP_Bool             lpwasprimfeas;      /**< primal feasibility of saved LP state information */
   SCIP_Bool             lpwasprimchecked;   /**< primal feasibility check state of saved LP state information  */
   SCIP_Bool             lpwasdualfeas;      /**< dual feasibility of saved LP state information */
   SCIP_Bool             lpwasdualchecked;   /**< dual feasibility check state of saved LP state information  */
};

/** sibling information (should not exceed the size of a pointer) */
struct SCIP_Sibling
{
   int                   arraypos;           /**< position of node in the siblings array */
};

/** child information (should not exceed the size of a pointer) */
struct SCIP_Child
{
   int                   arraypos;           /**< position of node in the children array */
};

/** leaf information (should not exceed the size of a pointer) */
struct SCIP_Leaf
{
   SCIP_NODE*            lpstatefork;        /**< fork/subroot node defining the LP state of the leaf */
};

/** fork without LP solution, where only bounds and constraints have been changed */
struct SCIP_Junction
{
   int                   nchildren;          /**< number of children of this parent node */
};

/** fork without LP solution, where bounds and constraints have been changed, and rows and columns were added */
struct SCIP_Pseudofork
{
   SCIP_COL**            addedcols;          /**< array with pointers to new columns added at this node into the LP */
   SCIP_ROW**            addedrows;          /**< array with pointers to new rows added at this node into the LP */
   int                   naddedcols;         /**< number of columns added at this node */
   int                   naddedrows;         /**< number of rows added at this node */
   int                   nchildren;          /**< number of children of this parent node */
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were added */
struct SCIP_Fork
{
   SCIP_COL**            addedcols;          /**< array with pointers to new columns added at this node into the LP */
   SCIP_ROW**            addedrows;          /**< array with pointers to new rows added at this node into the LP */
   SCIP_LPISTATE*        lpistate;           /**< LP state information */
   SCIP_Real             lpobjval;           /**< the LP objective value for that node, needed to compute the pseudo costs correctly */
   int                   naddedcols;         /**< number of columns added at this node */
   int                   naddedrows;         /**< number of rows added at this node */
   int                   nlpistateref;       /**< number of times, the LP state is needed */
   unsigned int          nchildren:28;       /**< number of children of this parent node */
   unsigned int          lpwasprimfeas:1;    /**< primal feasibility of saved LP state information */
   unsigned int          lpwasprimchecked:1; /**< primal feasibility check state of saved LP state information */
   unsigned int          lpwasdualfeas:1;    /**< dual feasibility of saved LP state information */
   unsigned int          lpwasdualchecked:1; /**< dual feasibility check state of saved LP state information */
};

/** fork with solved LP, where bounds and constraints have been changed, and rows and columns were removed and added */
struct SCIP_Subroot
{
   SCIP_COL**            cols;               /**< array with pointers to the columns in the same order as in the LP */
   SCIP_ROW**            rows;               /**< array with pointers to the rows in the same order as in the LP */
   SCIP_LPISTATE*        lpistate;           /**< LP state information */
   SCIP_Real             lpobjval;           /**< the LP objective value for that node, needed to compute the pseudo costs correctly */
   int                   ncols;              /**< number of columns in the LP */
   int                   nrows;              /**< number of rows in the LP */
   int                   nlpistateref;       /**< number of times, the LP state is needed */
   unsigned int          nchildren:30;       /**< number of children of this parent node */
   unsigned int          lpwasprimfeas:1;    /**< primal feasibility of saved LP state information */
   unsigned int          lpwasprimchecked:1; /**< primal feasibility check state of saved LP state information */
   unsigned int          lpwasdualfeas:1;    /**< dual feasibility of saved LP state information */
   unsigned int          lpwasdualchecked:1; /**< dual feasibility check state of saved LP state information */
};

/** node data structure */
struct SCIP_Node
{
   SCIP_Longint          number;             /**< successively assigned number of the node */
   SCIP_Real             lowerbound;         /**< lower (dual) bound of subtree */
   SCIP_Real             estimate;           /**< estimated value of feasible solution in subtree */
   union
   {
      SCIP_PROBINGNODE*  probingnode;        /**< data for probing nodes */
      SCIP_SIBLING       sibling;            /**< data for sibling nodes */
      SCIP_CHILD         child;              /**< data for child nodes */
      SCIP_LEAF          leaf;               /**< data for leaf nodes */
      SCIP_JUNCTION      junction;           /**< data for junction nodes */
      SCIP_PSEUDOFORK*   pseudofork;         /**< data for pseudo fork nodes */
      SCIP_FORK*         fork;               /**< data for fork nodes */
      SCIP_SUBROOT*      subroot;            /**< data for subroot nodes */
   } data;
   SCIP_NODE*            parent;             /**< parent node in the tree */
   SCIP_CONSSETCHG*      conssetchg;         /**< constraint set changes at this node or NULL */
   SCIP_DOMCHG*          domchg;             /**< domain changes at this node or NULL */
   unsigned int          depth:16;           /**< depth in the tree */
   unsigned int          nodetype:4;         /**< type of node */
   unsigned int          active:1;           /**< is node in the path to the current node? */
   unsigned int          cutoff:1;           /**< should the node and all sub nodes be cut off from the tree? */
   unsigned int          reprop:1;           /**< should propagation be applied again, if the node is on the active path? */
   unsigned int          repropsubtreemark:9;/**< subtree repropagation marker for subtree repropagation */
   unsigned int          reoptid:29;         /**< unique id to identify the node during reoptimization */
   unsigned int          reopttype:3;        /**< node type during reoptimization */
};

/** bound change information for pending bound changes */
struct SCIP_PendingBdchg
{
   SCIP_NODE*            node;               /**< node to add bound change to */
   SCIP_VAR*             var;                /**< variable to change the bounds for */
   SCIP_Real             newbound;           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype;          /**< type of bound: lower or upper bound */
   SCIP_CONS*            infercons;          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop;          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo;          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             probingchange;      /**< is the bound change a temporary setting due to probing? */
};

/** branch and bound tree */
struct SCIP_Tree
{
   SCIP_NODE*            root;               /**< root node of the tree */
   SCIP_NODEPQ*          leaves;             /**< leaves of the tree */
   SCIP_NODE**           path;               /**< array of nodes storing the active path from root to current node, which
                                              *   is usually the focus or a probing node; in case of a cut off, the path
                                              *   may already end earlier */
   SCIP_NODE*            focusnode;          /**< focus node: the node that is stored together with its children and
                                              *   siblings in the tree data structure; the focus node is the currently
                                              *   processed node; it doesn't need to be active all the time, because it
                                              *   may be cut off and the active path stops at the cut off node */
   SCIP_NODE*            focuslpfork;        /**< LP defining pseudofork/fork/subroot of the focus node */
   SCIP_NODE*            focuslpstatefork;   /**< LP state defining fork/subroot of the focus node */
   SCIP_NODE*            focussubroot;       /**< subroot of the focus node's sub tree */
   SCIP_NODE*            probingroot;        /**< root node of the current probing path, or NULL */
   SCIP_NODE**           children;           /**< array with children of the focus node */
   SCIP_NODE**           siblings;           /**< array with siblings of the focus node */
   SCIP_Real*            childrenprio;       /**< array with node selection priorities of children */
   SCIP_Real*            siblingsprio;       /**< array with node selection priorities of siblings */
   SCIP_VAR**            divebdchgvars[2];   /**< two arrays to store variables for branching */
   SCIP_BRANCHDIR*       divebdchgdirs[2];   /**< arrays to hold the directions for diving */
   SCIP_Real*            divebdchgvals[2];   /**< arrays to store bound change values for diving */
   int*                  pathnlpcols;        /**< array with number of LP columns for each problem in active path (except
                                              *   newly added columns of the focus node and the current probing node) */
   int*                  pathnlprows;        /**< array with number of LP rows for each problem in active path (except
                                              *   newly added rows of the focus node and the current probing node) */
   SCIP_LPISTATE*        probinglpistate;    /**< LP state information before probing started */
   SCIP_LPISTATE*        focuslpistate;      /**< LP state information of focus node */
   SCIP_LPINORMS*        probinglpinorms;    /**< LP pricing norms information before probing started */
   SCIP_PENDINGBDCHG*    pendingbdchgs;      /**< array of pending bound changes, or NULL */
   SCIP_Real*            probdiverelaxsol;   /**< array with stored original relaxation solution during diving or probing */
   int                   nprobdiverelaxsol;  /**< size of probdiverelaxsol */
   SCIP_Longint          focuslpstateforklpcount; /**< LP number of last solved LP in current LP state fork, or -1 if unknown */
   SCIP_Longint          lastbranchparentid; /**< last node id/number of branching parent */
   int                   divebdchgsize[2];   /**< holds the two sizes of the dive bound change information */
   int                   ndivebdchanges[2];  /**< current number of stored dive bound changes for the next depth */
   int                   pendingbdchgssize;  /**< size of pendingbdchgs array */
   int                   npendingbdchgs;     /**< number of pending bound changes */
   int                   childrensize;       /**< available slots in children vector */
   int                   nchildren;          /**< number of children of focus node (number of used slots in children vector) */
   int                   siblingssize;       /**< available slots in siblings vector */
   int                   nsiblings;          /**< number of siblings of focus node (number of used slots in siblings vector) */
   int                   pathlen;            /**< length of the current path */
   int                   pathsize;           /**< number of available slots in path arrays */
   int                   effectiverootdepth; /**< first depth with node with at least two children */
   int                   appliedeffectiverootdepth; /**< the effective root depth which was already enforced (that is constraint and bound changes were made global) */
   int                   correctlpdepth;     /**< depth to which current LP data corresponds to LP data of active path */
   int                   cutoffdepth;        /**< depth of first node in active path that is marked being cutoff */
   int                   repropdepth;        /**< depth of first node in active path that has to be propagated again */
   int                   repropsubtreecount; /**< cyclicly increased counter to create markers for subtree repropagation */
   int                   probingsumchgdobjs; /**< number of changed objective coefficients in all probing nodes */
   SCIP_Bool             focusnodehaslp;     /**< is LP being processed in the focus node? */
   SCIP_Bool             probingnodehaslp;   /**< was the LP solved (at least once) in the current probing node? */
   SCIP_Bool             focuslpconstructed; /**< was the LP of the focus node already constructed? */
   SCIP_Bool             cutoffdelayed;      /**< the treeCutoff() call was delayed because of diving and has to be executed */
   SCIP_Bool             probinglpwasflushed;/**< was the LP flushed before we entered the probing mode? */
   SCIP_Bool             probinglpwassolved; /**< was the LP solved before we entered the probing mode? */
   SCIP_Bool             probingloadlpistate;/**< must the LP state be reloaded because of a backtrack in probing? */
   SCIP_Bool             probinglpwasrelax;  /**< was the LP a valid relaxation before we entered the probing mode? */
   SCIP_Bool             probingsolvedlp;    /**< was the LP solved during probing mode, i.e., was SCIPsolveProbingLP() called? */
   SCIP_Bool             forcinglpmessage;   /**< was forcing LP solving message be posted */
   SCIP_Bool             probingobjchanged;  /**< was the objective function changed during probing? */
   SCIP_Bool             sbprobing;          /**< is the probing mode used for strong branching? */
   SCIP_Bool             probinglpwasprimfeas;/**< primal feasibility when probing started */
   SCIP_Bool             probinglpwasprimchecked;/**< primal feasibility has been checked when probing started */
   SCIP_Bool             probinglpwasdualfeas;/**< dual feasibility when probing started */
   SCIP_Bool             probinglpwasdualchecked;/**< dual feasibility has been check when probing started */
   SCIP_Bool             probdiverelaxstored; /**< was a relax solution stored before diving or probing ? */
   SCIP_Bool             probdiverelaxincludeslp; /**< did the stored relaxation solution include all lp cuts ? */
};

#ifdef __cplusplus
}
#endif

#endif
