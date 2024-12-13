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

/**@file   type_nodesel.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for node selectors
 * @author Tobias Achterberg
 */

/** @defgroup DEFPLUGINS_NODESEL Default node selectors
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default node selectors of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_NODESEL_H__
#define __SCIP_TYPE_NODESEL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_tree.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_NodePQ SCIP_NODEPQ;           /**< node priority queue */
typedef struct SCIP_Nodesel SCIP_NODESEL;         /**< node selector data structure */
typedef struct SCIP_NodeselData SCIP_NODESELDATA; /**< node selector specific data */


/** copy method for node selector plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define SCIP_DECL_NODESELCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_NODESEL* nodesel)


/** destructor of node selector to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define SCIP_DECL_NODESELFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_NODESEL* nodesel)

/** initialization method of node selector (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define SCIP_DECL_NODESELINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_NODESEL* nodesel)

/** deinitialization method of node selector (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define SCIP_DECL_NODESELEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_NODESEL* nodesel)

/** solving process initialization method of node selector (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The node selector may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define SCIP_DECL_NODESELINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_NODESEL* nodesel)

/** solving process deinitialization method of node selector (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The node selector should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 */
#define SCIP_DECL_NODESELEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_NODESEL* nodesel)

/** node selection method of node selector
 *
 *  This method is called to select the next leaf of the branch and bound tree to be processed.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 *  - selnode         : pointer to store the selected node
 *
 *  possible return values for *selnode:
 *  - NULL    : problem is solved, because tree is empty
 *  - non-NULL: node to be solved next
 */
#define SCIP_DECL_NODESELSELECT(x) SCIP_RETCODE x (SCIP* scip, SCIP_NODESEL* nodesel, SCIP_NODE** selnode)

/** node comparison method of node selector
 *
 *  This method is called to compare two nodes regarding their order in the node priority queue.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - nodesel         : the node selector itself
 *  - node1           : first node to compare
 *  - node2           : second node to compare
 *
 *  possible return values:
 *  - value < 0: node1 comes before (is better than) node2
 *  - value = 0: both nodes are equally good
 *  - value > 0: node2 comes after (is worse than) node2
 */
#define SCIP_DECL_NODESELCOMP(x) int x (SCIP* scip, SCIP_NODESEL* nodesel, SCIP_NODE* node1, SCIP_NODE* node2)

#ifdef __cplusplus
}
#endif

#endif
