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

/**@file   struct_nodesel.h
 * @ingroup INTERNALAPI
 * @brief  data structures for node selectors and node priority queues
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_NODESEL_H__
#define __SCIP_STRUCT_NODESEL_H__


#include "scip/def.h"
#include "scip/type_tree.h"
#include "scip/type_nodesel.h"

#ifdef __cplusplus
extern "C" {
#endif

/** node priority queue data structure;
 *  the fields lowerboundnode, lowerbound, nlowerbounds and validlowerbound are only used for node selection rules,
 *  that don't store the lowest bound node in the first slot of the queue
 */
struct SCIP_NodePQ
{
   SCIP_Real             lowerboundsum;      /**< sum of lower bounds of all nodes in the queue */
   SCIP_NODESEL*         nodesel;            /**< node selector used for sorting the nodes in the queue */
   SCIP_NODE**           slots;              /**< array of element slots */
   int*                  bfsposs;            /**< position of the slot in the bfs ordered queue */
   int*                  bfsqueue;           /**< queue of slots[] indices sorted by best lower bound */
   int                   len;                /**< number of used element slots */
   int                   size;               /**< total number of available element slots */
};

/** node selector */
struct SCIP_Nodesel
{
   char*                 name;               /**< name of node selector */
   char*                 desc;               /**< description of node selector */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy));   /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_NODESELFREE ((*nodeselfree));   /**< destructor of node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit));   /**< initialize node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit));   /**< deinitialize node selector */
   SCIP_DECL_NODESELINITSOL((*nodeselinitsol));/**< solving process initialization method of node selector */
   SCIP_DECL_NODESELEXITSOL((*nodeselexitsol));/**< solving process deinitialization method of node selector */
   SCIP_DECL_NODESELSELECT((*nodeselselect));/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp));   /**< node comparison method */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this node selector for the next stages */
   SCIP_CLOCK*           nodeseltime;        /**< node selector execution time */
   SCIP_NODESELDATA*     nodeseldata;        /**< node selector data */
   int                   stdpriority;        /**< priority of the node selector in standard mode */
   int                   memsavepriority;    /**< priority of the node selector in memory saving mode */
   SCIP_Bool             initialized;        /**< is node selector initialized? */
};

#ifdef __cplusplus
}
#endif

#endif
