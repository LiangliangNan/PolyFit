/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*                    TCLIQUE --- Algorithm for Maximum Cliques              */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  TCLIQUE is distributed under the terms of the ZIB Academic License.      */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with TCLIQUE; see the file COPYING.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   tclique_coloring.h
 * @brief  coloring part of algorithm for maximum cliques
 * @author Tobias Achterberg
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_COLORING_H__
#define __TCLIQUE_COLORING_H__

#include "blockmemshell/memory.h"
#include "tclique/tclique.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _ITV
{
   int inf;
   int sup;
} ITV;

typedef struct _LIST_ITV
{
   ITV itv;
   struct _LIST_ITV *next;
} LIST_ITV;

typedef struct _NBC
{
   int satdeg;
   LIST_ITV *lcitv;
} NBC;




/** colors the positive weighted nodes of a given set of nodes V with the lowest possible number of colors and
 *  finds a clique in the graph induced by V, an upper bound and an apriori bound for further branching steps */
extern
TCLIQUE_WEIGHT tcliqueColoring( 
   TCLIQUE_GETNNODES((*getnnodes)),          /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_SELECTADJNODES((*selectadjnodes)),/**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   BMS_CHKMEM*           mem,                /**< block memory */
   int*                  buffer,             /**< buffer of size nnodes */
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV,                 /**< number of non-zero weighted nodes for branching */
   NBC*                  gsd,                /**< neighbour color information of all nodes */
   TCLIQUE_Bool*         iscolored,          /**< coloring status of all nodes */
   TCLIQUE_WEIGHT*       apbound,            /**< pointer to store apriori bound of nodes for branching */
   int*                  clique,             /**< buffer for storing the clique */
   int*                  nclique,            /**< pointer to store number of nodes in the clique */
   TCLIQUE_WEIGHT*       weightclique        /**< pointer to store the weight of the clique */
   );

#ifdef __cplusplus
}
#endif

#endif
