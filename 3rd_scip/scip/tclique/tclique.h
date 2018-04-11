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

/**@file   tclique.h
 * @brief  tclique user interface
 * @author Tobias Achterberg
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_H__
#define __TCLIQUE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "tclique/tclique_def.h"

/*
 * Data Types and Structures
 */

typedef int  TCLIQUE_WEIGHT;                 /**< type used for node weights in the graph */
typedef struct TCLIQUE_Graph TCLIQUE_GRAPH;  /**< user defined structure for storing the graph, passed to graph callbacks */
typedef struct TCLIQUE_Data TCLIQUE_DATA;    /**< user defined data to pass to new solution callback method */

#ifndef TCLIQUE_Bool
#define TCLIQUE_Bool unsigned int            /**< type used for boolean values */
#endif
#ifndef TRUE
#define TRUE  1                              /**< boolean value TRUE */
#define FALSE 0                              /**< boolean value FALSE */
#endif

/** return status of the TCLIQUE algorithm */
enum TCLIQUE_Status
{
   TCLIQUE_ERROR           = 0,         /**< an error occurred */
   TCLIQUE_NODELIMIT       = 1,         /**< the node limit was reached */
   TCLIQUE_USERABORT       = 2,         /**< the user call back function aborted the solving process */
   TCLIQUE_OPTIMAL         = 3          /**< the optimal solution was found */
};
typedef enum TCLIQUE_Status TCLIQUE_STATUS;




/*
 * User Callback Methods
 */

/** user callback method which is called whenever a feasible clique was found
 *  input:
 *   - tcliquedata  : user data given to tcliqueMaxClique()
 *   - cliquenodes  : array with nodes of the clique
 *   - ncliquenodes : number of nodes in the clique
 *   - cliqueweight : weight of the clique
 *  output:
 *   - minweight    : new minimal weight for feasible cliques
 *   - acceptsol    : setting TRUE makes clique the new best clique, and updates minweight
 *   - stopsolving  : setting TRUE aborts the search for cliques
 */
#define TCLIQUE_NEWSOL(x) void x (TCLIQUE_DATA* tcliquedata, int* cliquenodes, int ncliquenodes, \
      TCLIQUE_WEIGHT cliqueweight, TCLIQUE_WEIGHT* minweight, TCLIQUE_Bool* acceptsol, TCLIQUE_Bool* stopsolving)

/** user callback method to get number of nodes in the graph
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *  returns:
 *   number of nodes in the graph
 */
#define TCLIQUE_GETNNODES(x) int x (TCLIQUE_GRAPH* tcliquegraph)

/** user callback method to get weights of nodes in the graph
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *  returns:
 *   array of node weights (of length at least equal to the number of nodes in the graph)
 */
#define TCLIQUE_GETWEIGHTS(x) const TCLIQUE_WEIGHT* x (TCLIQUE_GRAPH* tcliquegraph)

/** user callback method to return whether the edge (node1, node2) is in the graph
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *   - node1        : start node of edge (tail)
 *   - node2        : end node of edge (head)
 *  returns:
 *   TRUE if edge is in the graph, FALSE otherwise
 */
#define TCLIQUE_ISEDGE(x) TCLIQUE_Bool x (TCLIQUE_GRAPH* tcliquegraph, int node1, int node2)

/** user callback method to select all nodes from a given set of nodes which are adjacent to a given node
 *  input:
 *   - tcliquegraph : user defined graph data structure given to tcliqueMaxClique()
 *   - node         : node to select adjacent nodes from
 *   - nodes        : array of nodes to select nodes from
 *   - nnodes       : number of nodes in 'nodes' array
 *   - adjnodes     : pointer to memory to store the resulting nodes
 *                    'adjnodes' and 'nodes' may point to the same memory location
 *  output:
 *   - adjnodes     : array of nodes that are contained in 'nodes' and that are adjacent to 'node'
 *  returns:
 *   number of nodes in 'adjnodes'
 */
#define TCLIQUE_SELECTADJNODES(x) int x (TCLIQUE_GRAPH* tcliquegraph, int node, int* nodes, int nnodes, int* adjnodes)




/*
 * Default Graph Implementation: Interface Methods used by the TClique algorithm
 */

/** gets number of nodes in the graph */
EXTERN
TCLIQUE_GETNNODES(tcliqueGetNNodes);

/** gets weight of nodes in the graph */
EXTERN
TCLIQUE_GETWEIGHTS(tcliqueGetWeights);

/** returns, whether the edge (node1, node2) is in the graph */
EXTERN
TCLIQUE_ISEDGE(tcliqueIsEdge);

/** selects all nodes from a given set of nodes which are adjacent to a given node
 *  and returns the number of selected nodes */
EXTERN
TCLIQUE_SELECTADJNODES(tcliqueSelectAdjnodes);




/*
 * Default Graph Implementation: External Interface Methods to access the graph
 */

/** creates graph data structure */
EXTERN
TCLIQUE_Bool tcliqueCreate(
   TCLIQUE_GRAPH**       tcliquegraph        /**< pointer to store graph data structure */
   );

/** frees graph data structure */
EXTERN
void tcliqueFree(
   TCLIQUE_GRAPH**       tcliquegraph        /**< pointer to graph data structure */
   );

/** adds nodes up to the given node number to graph data structure (intermediate nodes have weight 0) */
EXTERN
TCLIQUE_Bool tcliqueAddNode(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   node,               /**< node number to add */
   TCLIQUE_WEIGHT        weight              /**< weight of node to add */
   );

/** changes weight of node in graph data structure */
EXTERN
void tcliqueChangeWeight(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   node,               /**< node to set new weight */
   TCLIQUE_WEIGHT        weight              /**< new weight of node (allready scaled) */
   );

/** adds edge (node1, node2) to graph data structure (node1 and node2 have to be contained in
 *  graph data structure)
 *
 *  New edges are cached, s.t. the graph data structures are not correct until a call to tcliqueFlush();
 *  you have to make sure, that no double edges are inserted.
 */
EXTERN
TCLIQUE_Bool tcliqueAddEdge(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   node1,              /**< start node of edge to add */
   int                   node2               /**< end node of edge to add */
   );

/** inserts all cached edges into the data structures */
EXTERN
TCLIQUE_Bool tcliqueFlush(
   TCLIQUE_GRAPH*        tcliquegraph        /**< graph data structure */
   );

/** loads graph data structure from file */
EXTERN
TCLIQUE_Bool tcliqueLoadFile(
   TCLIQUE_GRAPH**       tcliquegraph,       /**< pointer to store graph data structure */
   const char*           filename,           /**< name of file with graph data */
   double                scaleval,           /**< value to scale weights (only integral part of scaled weights is considered) */
   char*                 probname,           /**< buffer to store the name of the problem */
   int                   sizeofprobname      /**< size of buffer to store the name of the problem */
   );

/** saves graph data structure to file */
EXTERN
TCLIQUE_Bool tcliqueSaveFile(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   const char*           filename,           /**< name of file to create */
   double                scaleval,           /**< value to unscale weights with */
   const char*           probname            /**< name of the problem */
   );

/** gets number of edges in the graph */
EXTERN
int tcliqueGetNEdges(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   );

/** gets degree of nodes in graph */
EXTERN
int* tcliqueGetDegrees(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   );

/** gets adjacent nodes of edges in graph */
EXTERN
int* tcliqueGetAdjnodes(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   );

/** gets pointer to first adjacent edge of given node in graph */
EXTERN
int* tcliqueGetFirstAdjedge(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   int                   node                /**< given node */
   );

/** gets pointer to last adjacent edge of given node in graph */
EXTERN
int* tcliqueGetLastAdjedge(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   int                   node                /**< given node */
   );

/** prints graph data structure */
EXTERN
void tcliquePrintGraph(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   );




/*
 * Interface Methods
 */

/** finds maximum weight clique */
EXTERN
void tcliqueMaxClique(
   TCLIQUE_GETNNODES((*getnnodes)),          /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_ISEDGE((*isedge)),                /**< user function to check for existence of an edge */
   TCLIQUE_SELECTADJNODES((*selectadjnodes)), /**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure that is passed to graph callbacks */
   TCLIQUE_NEWSOL((*newsol)),                /**< user function to call on every new solution */
   TCLIQUE_DATA*         tcliquedata,        /**< user data to pass to new solution callback function */
   int*                  maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*                  nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   TCLIQUE_WEIGHT*       maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   TCLIQUE_WEIGHT        maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used
                                              *   for cliques with at least one fractional node) */
   TCLIQUE_WEIGHT        minweight,          /**< lower bound for weight of generated cliques */
   int                   maxntreenodes,	     /**< maximal number of nodes of b&b tree */
   int                   backtrackfreq,      /**< frequency to backtrack to first level of tree (0: no premature backtracking) */
   int                   maxnzeroextensions, /**< maximal number of zero-valued variables extending the clique */
   int                   fixednode,          /**< node that is forced to be in the clique, or -1; must have positive weight */
   int*                  ntreenodes,         /**< pointer to store the number of used tree nodes (or NULL) */
   TCLIQUE_STATUS*       status              /**< pointer to store the status of the solving call */
   );

#ifdef __cplusplus
}
#endif

#endif
