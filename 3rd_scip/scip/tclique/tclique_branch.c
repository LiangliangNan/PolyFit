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

/**@file   tclique_branch.c
 * @brief  branch and bound part of algorithm for maximum cliques
 * @author Tobias Achterberg
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

#include "tclique/tclique.h"
#include "tclique/tclique_def.h"
#include "tclique/tclique_coloring.h"
#include "blockmemshell/memory.h"



#define CHUNK_SIZE (64)
#define CLIQUEHASH_INITSIZE (1024)




/***************************
 * clique hash table methods
 ***************************/

typedef struct clique CLIQUE;                /**< single element of the clique hash table */

/** single element of the clique hash table */
struct clique
{
   int*                  nodes;              /**< node number of the clique elements */
   int                   nnodes;             /**< length of the clique */
};

typedef struct cliquehash CLIQUEHASH;   /**< hash table for storing cliques */

/** hash table for storing cliques */
struct cliquehash
{
   CLIQUE**              cliques;            /**< elements of the table */
   int                   cliquessize;        /**< size of the table */
   int                   ncliques;           /**< number of cliques stored in the table */
};


/** creates a clique */
static
void createClique(
   CLIQUE**              clique,             /**< pointer to the clique */
   int*                  nodes,              /**< nodes of the clique */
   int                   nnodes              /**< number of nodes in the clique */
   )
{
   int i;

   assert(clique != NULL);

   ALLOC_ABORT( BMSallocMemory(clique) );
   ALLOC_ABORT( BMSallocMemoryArray(&(*clique)->nodes, nnodes) );

   /* sort the nodes into the clique's node array */
   for( i = 0; i < nnodes; ++i )
   {
      int node;
      int j;

      node = nodes[i];
      for( j = i; j > 0 && node < (*clique)->nodes[j-1]; --j )
         (*clique)->nodes[j] = (*clique)->nodes[j-1];
      (*clique)->nodes[j] = node;
   }
   (*clique)->nnodes = nnodes;
}

/** frees a clique */
static
void freeClique(
   CLIQUE**              clique              /**< pointer to the clique */
   )
{
   assert(clique != NULL);
   assert(*clique != NULL);

   BMSfreeMemoryArray(&(*clique)->nodes);
   BMSfreeMemory(clique);
}

/** checks, whether clique1 is a subset of clique2 and returns the following value:
 *   == 0 if clique1 == clique2, or clique1 is contained in clique2,
 *    < 0 if clique1 < clique2, and clique1 is not contained in clique2,
 *    > 0 if clique1 > clique2, and clique1 is not contained in clique2
 */
static
int compSubcliques(
   CLIQUE*               clique1,            /**< first clique to compare */
   CLIQUE*               clique2             /**< second clique to compare */
   )
{
   int pos1;
   int pos2;
   TCLIQUE_Bool clique2smaller;

   assert(clique1 != NULL);
   assert(clique2 != NULL);

   pos1 = 0;
   pos2 = 0;
   clique2smaller = FALSE;
   while( pos1 < clique1->nnodes && pos2 < clique2->nnodes )
   {
      if( clique1->nodes[pos1] < clique2->nodes[pos2] )
      {
         /* clique1 has an element not contained in clique2: clique1 is lex-smaller, if clique2 was not
          * detected earlier to be lex-smaller
          */
         return (clique2smaller ? +1 : -1);
      }
      else if( clique1->nodes[pos1] > clique2->nodes[pos2] )
      {
         /* clique2 has an element not contained in clique1: clique2 is lex-smaller, but probably clique1 is
          * contained in clique2
          */
         pos2++;
         clique2smaller = TRUE;
      }
      else
      {
         pos1++;
         pos2++;
      }
   }

   /* if clique1 has additional elements, it is not contained in clique2 */
   if( pos1 < clique1->nnodes )
      return (clique2smaller ? +1 : -1);

   /* clique1 is contained in clique2 */
   return 0;
}

#ifndef NDEBUG
/** performs an integrity check of the clique hash table */
static
void checkCliquehash(
   CLIQUEHASH*           cliquehash          /**< clique hash table */
   )
{
   int i;

   assert(cliquehash != NULL);

   for( i = 0; i < cliquehash->ncliques-1; ++i )
      assert(compSubcliques(cliquehash->cliques[i], cliquehash->cliques[i+1]) < 0);
}
#else
#define checkCliquehash(cliquehash) /**/
#endif

/** creates a table for storing cliques */
static
void createCliquehash(
   CLIQUEHASH**          cliquehash,         /**< pointer to store the clique hash table */
   int                   tablesize           /**< initial size of the clique hash table */
   )
{
   assert(cliquehash != NULL);
   assert(tablesize > 0);

   ALLOC_ABORT( BMSallocMemory(cliquehash) );
   ALLOC_ABORT( BMSallocMemoryArray(&(*cliquehash)->cliques, tablesize) );
   (*cliquehash)->cliquessize = tablesize;
   (*cliquehash)->ncliques = 0;
}

/** clears the clique hash table and frees all inserted cliques */
static
void clearCliquehash(
   CLIQUEHASH*           cliquehash          /**< clique hash table */
   )
{
   int i;

   assert(cliquehash != NULL);

   /* free the cliques in the table */
   for( i = 0; i < cliquehash->ncliques; ++i )
      freeClique(&cliquehash->cliques[i]);

   cliquehash->ncliques = 0;
}

/** frees the table for storing cliques and all inserted cliques */
static
void freeCliquehash(
   CLIQUEHASH**          cliquehash          /**< pointer to the clique hash table */
   )
{
   assert(cliquehash != NULL);
   assert(*cliquehash != NULL);

   /* free the cliques in the table */
   clearCliquehash(*cliquehash);

   /* free the table data structure */
   BMSfreeMemoryArray(&(*cliquehash)->cliques);
   BMSfreeMemory(cliquehash);
}

/** ensures, that the clique hash table is able to store at least the given number of cliques */
static
void ensureCliquehashSize(
   CLIQUEHASH*           cliquehash,         /**< clique hash table */
   int                   num                 /**< minimal number of cliques to store */
   )
{
   assert(cliquehash != NULL);

   if( num > cliquehash->cliquessize )
   {
      int newsize;

      newsize = 2*cliquehash->cliquessize;
      if( num > newsize )
         newsize = num;

      ALLOC_ABORT( BMSreallocMemoryArray(&cliquehash->cliques, newsize) );
      cliquehash->cliquessize = newsize;
   }
   assert(cliquehash->cliquessize >= num);
}

#ifdef TCLIQUE_DEBUG
/** displayes clique hash table */
static
void printCliquehash(
   CLIQUEHASH*           cliquehash          /**< clique hash table */
   )
{
   int i;

   assert(cliquehash != NULL);

   debugMessage("cliquehash (%d cliques):\n", cliquehash->ncliques);
   for( i = 0; i < cliquehash->ncliques; ++i )
   {
      int j;

      debugPrintf("%d:", i);
      for( j = 0; j < cliquehash->cliques[i]->nnodes; ++j )
         debugPrintf(" %d", cliquehash->cliques[i]->nodes[j]);
      debugPrintf("\n");
   }
}
#endif

/** searches the given clique in the clique hash table and returns whether it (or a stronger clique) exists */
static
TCLIQUE_Bool inCliquehash(
   CLIQUEHASH*           cliquehash,         /**< clique hash table */
   CLIQUE*               clique,             /**< clique to search in the table */
   int*                  insertpos           /**< position where the clique should be inserted in the table */
   )
{
   int left;
   int right;
   int middle;
   int cmp;

   assert(cliquehash != NULL);
   assert(cliquehash->cliquessize > 0);
   assert(clique != NULL);
   assert(insertpos != NULL);

   /* perform a binary search on the clique hash table */
   left = 0;
   right = cliquehash->ncliques-1;
   while( left <= right )
   {
      middle = (left+right)/2;
      cmp = compSubcliques(clique, cliquehash->cliques[middle]);
      if( cmp > 0 )
         left = middle+1;
      else if( cmp < 0 )
         right = middle-1;
      else
      {
         *insertpos = middle;
         return TRUE;
      }
   }

   /* we found the correct insertion position for the clique, but it might be contained in a lex-smaller clique */
   *insertpos = left;
   for( middle = left-1; middle >= 0; --middle )
   {
      cmp = compSubcliques(clique, cliquehash->cliques[middle]);
      assert(cmp >= 0);
      if( cmp == 0 )
         return TRUE;
   }

   return FALSE;
}

/** inserts clique into clique hash table */
static
void insertClique(
   CLIQUEHASH*           cliquehash,         /**< clique hash table */
   CLIQUE*               clique,             /**< clique to search in the table */
   int                   insertpos           /**< position to insert clique into table (returned by inCliquehash()) */
   )
{
   int i;

   assert(cliquehash != NULL);
   assert(clique != NULL);
   assert(0 <= insertpos && insertpos <= cliquehash->ncliques);

   /* allocate memory */
   ensureCliquehashSize(cliquehash, cliquehash->ncliques+1);

   /* insert clique into table */
   for( i = cliquehash->ncliques; i > insertpos; --i )
      cliquehash->cliques[i] = cliquehash->cliques[i-1];
   cliquehash->cliques[insertpos] = clique;
   cliquehash->ncliques++;

   /* check, whether the clique hash table is still sorted */
   checkCliquehash(cliquehash);

   debug(printCliquehash(cliquehash));
}




/****************************
 * clique calculation methods
 ****************************/

/** extends given clique by additional zero-weight nodes of the given node set */
static
void extendCliqueZeroWeight(
   TCLIQUE_SELECTADJNODES((*selectadjnodes)),/**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   int*                  buffer,             /**< buffer of size nnodes */
   int*                  Vzero,              /**< zero weighted nodes */
   int                   nVzero,             /**< number of zero weighted nodes */
   int                   maxnzeroextensions, /**< maximal number of zero-valued variables extending the clique */
   int*                  curcliquenodes,     /**< nodes of the clique */
   int*                  ncurcliquenodes     /**< pointer to store number of nodes in the clique */
   )
{
   int i;
   int* zerocands;
   int nzerocands;
   int nzeroextensions;

   assert(selectadjnodes != NULL);
   assert(buffer != NULL);
   assert(Vzero != NULL);
   assert(curcliquenodes != NULL);
   assert(ncurcliquenodes != NULL);

   debugMessage("extending temporary clique (size %d) with zero-weighted nodes (nVzero=%d)\n", *ncurcliquenodes, nVzero);

   if( maxnzeroextensions == 0 )
      return;

   /* initialize the zero-weighted candidates for clique extension */
   zerocands = buffer;
   BMScopyMemoryArray(zerocands, Vzero, nVzero);
   nzerocands = nVzero;

   /* for each node in the clique, remove non-adjacent nodes from the zero extension candidates */
   for( i = 0; i < *ncurcliquenodes && nzerocands > 0; ++i )
   {
      nzerocands = selectadjnodes(tcliquegraph, curcliquenodes[i], zerocands, nzerocands, zerocands);
   }

   /* put zero-weighted candidates into the clique, and remove non-adjacent nodes from the candidate set */
   nzeroextensions = 0;
   while( nzerocands > 0 )
   {
      /* put first candidate into the clique */
      curcliquenodes[*ncurcliquenodes] = zerocands[0];
      (*ncurcliquenodes)++;
      nzerocands--;
      zerocands++;
      nzeroextensions++;
      if( nzeroextensions >= maxnzeroextensions )
         break;

      /* remove candidates that are not adjacent to the inserted zero-weighted node */
      nzerocands = selectadjnodes(tcliquegraph, curcliquenodes[(*ncurcliquenodes)-1], zerocands, nzerocands, zerocands);
   }
}

/** calls user callback after a new solution was found, that is better than the current incumbent
 *
 *  The callback decides, whether this solution should be accepted as new incumbent, and whether the solution process
 *  should be stopped.
 */
static
void newSolution(
   TCLIQUE_SELECTADJNODES((*selectadjnodes)), /**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   TCLIQUE_NEWSOL((*newsol)),                /**< user function to call on every new solution */
   TCLIQUE_DATA*         tcliquedata,        /**< user data to pass to user callback function */
   CLIQUEHASH*           cliquehash,         /**< clique hash table */
   int*                  buffer,             /**< buffer of size nnodes */
   int*                  Vzero,              /**< zero weighted nodes */
   int                   nVzero,             /**< number of zero weighted nodes */
   int                   maxnzeroextensions, /**< maximal number of zero-valued variables extending the clique */
   int*                  curcliquenodes,     /**< nodes of the new clique */
   int                   ncurcliquenodes,    /**< number of nodes in the new clique */
   TCLIQUE_WEIGHT        curcliqueweight,    /**< weight of the new clique */
   int*                  maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*                  nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   TCLIQUE_WEIGHT*       maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   TCLIQUE_Bool*         stopsolving         /**< pointer to store whether the solving should be stopped */
   )
{
   CLIQUE* clique;
   int insertpos;
   TCLIQUE_Bool acceptsol;

   assert(curcliquenodes != NULL);
   assert(maxcliquenodes != NULL);
   assert(nmaxcliquenodes != NULL);
   assert(maxcliqueweight != NULL);
   assert(curcliqueweight > *maxcliqueweight);
   assert(stopsolving != NULL);
   assert(newsol == NULL || cliquehash != NULL);

   acceptsol = TRUE;
   *stopsolving = FALSE;
   clique = NULL;
   insertpos = 0;

   if( newsol != NULL )
   {
      /* check whether the clique is already stored in the table */
      if( cliquehash->ncliques > 0 )
      {
         createClique(&clique, curcliquenodes, ncurcliquenodes);
         acceptsol = !inCliquehash(cliquehash, clique, &insertpos);
      }
   }

   /* check, if this is a new clique */
   if( acceptsol )
   {
      /* extend the clique with the zero-weighted nodes */
      extendCliqueZeroWeight(selectadjnodes, tcliquegraph, buffer, Vzero, nVzero, maxnzeroextensions,
         curcliquenodes, &ncurcliquenodes);

      if( newsol != NULL )
      {
         /* call user callback method */
         newsol(tcliquedata, curcliquenodes, ncurcliquenodes, curcliqueweight, maxcliqueweight, &acceptsol, stopsolving);

         /* if clique was accepted, clear the clique hash table; otherwise, insert it into the clique hash table, such that
          * the same or a weaker clique is not presented to the user again
          */
         if( acceptsol )
            clearCliquehash(cliquehash);
         else
         {
            /* if the clique was not yet created, do it now */
            if( clique == NULL )
            {
               assert(insertpos == 0);
               assert(cliquehash->ncliques == 0);
               createClique(&clique, curcliquenodes, ncurcliquenodes);
            }

            /* insert clique into clique hash table */
            insertClique(cliquehash, clique, insertpos);
            clique = NULL; /* the clique now belongs to the table */
         }
      }
   }

   /* free the clique, if it was created and not put into the clique hash table */
   if( clique != NULL )
      freeClique(&clique);

   if( acceptsol )
   {
      /* copy the solution to the incumbent */
      BMScopyMemoryArray(maxcliquenodes, curcliquenodes, ncurcliquenodes);
      *nmaxcliquenodes = ncurcliquenodes;
      if( curcliqueweight > *maxcliqueweight )
         *maxcliqueweight = curcliqueweight;
   }

#ifdef TCLIQUE_DEBUG
   debugMessage(" -> clique %s (weight %d):", acceptsol ? "accepted" : "rejected", curcliqueweight);
   {
      int i;
      for( i = 0; i < ncurcliquenodes; ++i )
         debugPrintf(" %d", curcliquenodes[i]);
      debugPrintf("\n");
   }
#endif
}

/** tries to find a clique, if V has only one or two nodes */
static
void reduced(
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_ISEDGE((*isedge)),                /**< user function to check for existence of an edge */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV,                 /**< number of non-zero weighted nodes for branching */
   TCLIQUE_WEIGHT*       apbound,            /**< apriori bound of nodes for branching */
   int*                  tmpcliquenodes,     /**< buffer for storing the temporary clique */
   int*                  ntmpcliquenodes,    /**< pointer to store number of nodes of the temporary clique */
   TCLIQUE_WEIGHT*       tmpcliqueweight     /**< pointer to store weight of the temporary clique */
   )
{
   const TCLIQUE_WEIGHT* weights;

   assert(getweights != NULL);
   assert(isedge != NULL);
   assert(tcliquegraph != NULL);
   assert(V != NULL);
   assert(0 <= nV && nV <= 2);
   assert(apbound != NULL);
   assert(tmpcliquenodes != NULL);
   assert(ntmpcliquenodes != NULL);
   assert(tmpcliqueweight != NULL);

   weights = getweights(tcliquegraph);
   assert(nV == 0 || weights[V[0]] > 0);
   assert(nV <= 1 || weights[V[1]] > 0);

   if( nV >= 1 )
      apbound[0] = weights[V[0]];
   if( nV >= 2 )
      apbound[1] = weights[V[1]];

   /* check if nodes are adjacent */
   if( nV >= 2 && isedge(tcliquegraph, V[0], V[1]) )
   {
      assert(isedge(tcliquegraph, V[1], V[0]));

      /* put nodes into clique */
      tmpcliquenodes[0] = V[0];
      tmpcliquenodes[1] = V[1];
      *ntmpcliquenodes = 2;
      *tmpcliqueweight = weights[V[0]] + weights[V[1]];
      apbound[0] += weights[V[1]];
   }
   else if( nV >= 2 && weights[V[1]] > weights[V[0]] )
   {
      /* put V[1] into clique */
      tmpcliquenodes[0] = V[1];
      *ntmpcliquenodes = 1;
      *tmpcliqueweight = weights[V[1]];
   }
   else if( nV >= 1 )
   {
      /* put V[0] into clique */
      tmpcliquenodes[0] = V[0];
      *ntmpcliquenodes = 1;
      *tmpcliqueweight = weights[V[0]];
   }
   else
   {
      *tmpcliqueweight = 0;
      *ntmpcliquenodes = 0;
   }
}

/** calculates upper bound on remaining subgraph, and heuristically generates a clique */
static
TCLIQUE_WEIGHT boundSubgraph(
   TCLIQUE_GETNNODES((*getnnodes)),          /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_ISEDGE((*isedge)),                /**< user function to check for existence of an edge */
   TCLIQUE_SELECTADJNODES((*selectadjnodes)), /**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   BMS_CHKMEM*           mem,                /**< block memory */
   int*                  buffer,             /**< buffer of size nnodes */
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV,                 /**< number of non-zero weighted nodes for branching */
   NBC*                  gsd,                /**< neighbour color information of all nodes */
   TCLIQUE_Bool*         iscolored,          /**< coloring status of all nodes */
   TCLIQUE_WEIGHT*       apbound,            /**< apriori bound of nodes for branching */
   int*                  tmpcliquenodes,     /**< buffer for storing the temporary clique */
   int*                  ntmpcliquenodes,    /**< pointer to store number of nodes of the temporary clique */
   TCLIQUE_WEIGHT*       tmpcliqueweight     /**< pointer to store weight of the temporary clique */
   )
{
   assert(tmpcliqueweight != NULL);

   /* check if we are in an easy case with at most 2 nodes left */
   if( nV <= 2 )
   {
      /* get 1- or 2-clique and bounds without coloring */
      reduced(getweights, isedge, tcliquegraph, V, nV, apbound, tmpcliquenodes, ntmpcliquenodes, tmpcliqueweight);
      return *tmpcliqueweight;
   }
   else
   {
      /* color the graph induces by nodes of V to get an upper bound for the remaining subgraph */
      return tcliqueColoring(getnnodes, getweights, selectadjnodes, tcliquegraph,
         mem, buffer, V, nV, gsd, iscolored, apbound,
         tmpcliquenodes, ntmpcliquenodes, tmpcliqueweight);
   }
}

/** gets the index of the node of V with the maximum apriori bound; returns -1, if no node exists */
static
int getMaxApBoundIndex(
   int                   nV,                 /**< number of nodes of V */
   TCLIQUE_WEIGHT*       apbound             /**< apriori bound of nodes of V */
   )
{
   TCLIQUE_WEIGHT maxapbound;
   int maxindex;
   int i;

   assert(apbound != NULL);

   maxapbound = 0;
   maxindex = -1;

   for( i = 0 ; i < nV; i++ )
   {
      assert(apbound[i] > 0);
      if( apbound[i] >= maxapbound )
      {
         maxapbound = apbound[i];
         maxindex = i;
      }
   }

   return maxindex;
}

/** gets the index of the node of V with the maximum apriori bound, but ignores nodes with weights
 *  larger than the given maximal weight
 *
 *  Returns -1 if no node with weight smaller or equal than maxweight is found.
 */
static
int getMaxApBoundIndexNotMaxWeight(
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV,                 /**< number of non-zero weighted nodes for branching */
   const TCLIQUE_WEIGHT* apbound,            /**< apriori bound of nodes of V */
   const TCLIQUE_WEIGHT* weights,            /**< weights of nodes */
   TCLIQUE_WEIGHT        maxweight           /**< maximal weight of node to be candidate for selection */
   )
{
   TCLIQUE_WEIGHT maxapbound;
   int maxindex;
   int i;

   assert(apbound != NULL);

   maxapbound = 0;
   maxindex = -1;

   for( i = 0 ; i < nV; i++ )
   {
      assert(apbound[i] > 0);
      assert(weights[V[i]] > 0);

      if( apbound[i] >= maxapbound && weights[V[i]] <= maxweight )
      {
         maxapbound = apbound[i];
         maxindex = i;
      }
   }

   return maxindex;
}

/** branches the searching tree, branching nodes are selected in decreasing order of their apriori bound,
 *  returns the level to which we should backtrack, or INT_MAX for continuing normally
 */
static
int branch(
   TCLIQUE_GETNNODES((*getnnodes)),          /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_ISEDGE((*isedge)),                /**< user function to check for existence of an edge */
   TCLIQUE_SELECTADJNODES((*selectadjnodes)),/**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   TCLIQUE_NEWSOL((*newsol)),                /**< user function to call on every new solution */
   TCLIQUE_DATA*         tcliquedata,        /**< user data to pass to user callback function */
   BMS_CHKMEM*           mem,                /**< block memory */
   CLIQUEHASH*           cliquehash,         /**< clique hash table */
   int*                  buffer,             /**< buffer of size nnodes */
   int                   level,              /**< level of b&b tree */
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV,                 /**< number of non-zero weighted nodes for branching */
   int*                  Vzero,              /**< zero weighted nodes */
   int                   nVzero,             /**< number of zero weighted nodes */
   NBC*                  gsd,                /**< neighbour color information of all nodes */
   TCLIQUE_Bool*         iscolored,          /**< coloring status of all nodes */
   int*                  K,                  /**< nodes from the b&b tree */
   TCLIQUE_WEIGHT        weightK,            /**< weight of the nodes from b&b tree */
   int*                  maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*                  nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   TCLIQUE_WEIGHT*       maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   int*                  curcliquenodes,     /**< pointer to store nodes of currenct clique */
   int*                  ncurcliquenodes,    /**< pointer to store number of nodes in current clique */
   TCLIQUE_WEIGHT*       curcliqueweight,    /**< pointer to store weight of current clique */
   int*                  tmpcliquenodes,     /**< buffer for storing the temporary clique */
   TCLIQUE_WEIGHT        maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used
                                              *   (for cliques with at least one fractional node) */
   int*                  ntreenodes,         /**< pointer to store number of nodes of b&b tree */
   int                   maxntreenodes,      /**< maximal number of nodes of b&b tree */
   int                   backtrackfreq,      /**< frequency to backtrack to first level of tree (0: no premature backtracking) */
   int                   maxnzeroextensions, /**< maximal number of zero-valued variables extending the clique */
   int                   fixednode,          /**< node that is forced to be in the clique, or -1; must have positive weight */
   TCLIQUE_STATUS*       status              /**< pointer to store the status of the solving call */
   )
{
   TCLIQUE_Bool isleaf;
   const TCLIQUE_WEIGHT* weights;
   TCLIQUE_WEIGHT* apbound;
   TCLIQUE_WEIGHT subgraphweight;
   TCLIQUE_WEIGHT weightKold;
   TCLIQUE_WEIGHT tmpcliqueweight;
   int backtracklevel;
   int ntmpcliquenodes;
   int i;

   assert(getnnodes != NULL);
   assert(getweights != NULL);
   assert(selectadjnodes != NULL);
   assert(mem != NULL);
   assert(V != NULL);
   assert(gsd != NULL);
   assert(iscolored != NULL);
   assert(K != NULL);
   assert(maxcliqueweight != NULL);
   assert(curcliquenodes != NULL);
   assert(ncurcliquenodes != NULL);
   assert(curcliqueweight != NULL);
   assert(ntreenodes != NULL);
   assert(maxfirstnodeweight >= 0);
   assert(*ntreenodes >= 0);
   assert(maxntreenodes >= 0);
   assert(status != NULL);

   /* increase the number of nodes, and stop solving, if the node limit is exceeded */
   (*ntreenodes)++;
#ifdef TCLIQUE_DEBUG
   debugMessage("(level %d, treenode %d) maxclique = %d, curclique = %d [mem=%" SCIP_LONGINT_FORMAT " (%" SCIP_LONGINT_FORMAT "), cliques=%d]\n",
      level, *ntreenodes, *maxcliqueweight, *curcliqueweight,
      BMSgetChunkMemoryUsed(mem), BMSgetMemoryUsed(), cliquehash == NULL ? 0 : cliquehash->ncliques);

   debugMessage(" -> current branching (weight %d):", weightK);
   for( i = 0; i < level; ++i )
      debugPrintf(" %d", K[i]);
   debugPrintf("\n");
   debugMessage(" -> branching candidates:");
   for( i = 0; i < nV; ++i )
      debugPrintf(" %d", V[i]);
   debugPrintf("\n");
#endif
   if( *ntreenodes > maxntreenodes )
   {
      *status = TCLIQUE_NODELIMIT;
      return TRUE;
   }

   weights = getweights(tcliquegraph);
   backtracklevel = INT_MAX;
   isleaf = TRUE;

   /* allocate temporary memory for a priori bounds */
   ALLOC_ABORT( BMSallocMemoryArray(&apbound, nV) );
   BMSclearMemoryArray(apbound, nV);

   /* use coloring relaxation to generate an upper bound for the current subtree and a heuristic solution */
   subgraphweight = boundSubgraph(getnnodes, getweights, isedge, selectadjnodes, tcliquegraph,
      mem, buffer, V, nV, gsd, iscolored, apbound,
      tmpcliquenodes, &ntmpcliquenodes, &tmpcliqueweight);

#ifndef NDEBUG
   /* check correctness of V and apbound arrays */
   for( i = 0; i < nV; ++i )
   {
      assert(0 <= V[i] && V[i] < getnnodes(tcliquegraph));
      assert(i == 0 || V[i-1] < V[i]);
      assert(apbound[i] >= 0);
      assert((apbound[i] == 0) == (weights[V[i]] == 0));
   }
#endif

   /* check, whether the heuristic solution is better than the current subtree's solution;
    * if the user wanted to have a fixed variable inside the clique and we are in level 0, we first have to
    * fix this variable in this level (the current clique might not contain the fixed node)
    */
   if( weightK + tmpcliqueweight > *curcliqueweight && (level > 0 || fixednode == -1) )
   {
      /* install the newly generated clique as current clique */
      for( i = 0; i < level; ++i )
         curcliquenodes[i] = K[i];
      for( i = 0; i < ntmpcliquenodes; ++i )
         curcliquenodes[level+i] = tmpcliquenodes[i];
      *ncurcliquenodes = level + ntmpcliquenodes;
      *curcliqueweight = weightK + tmpcliqueweight;

#ifdef TCLIQUE_DEBUG
      debugMessage(" -> new current clique with weight %d at node %d in level %d:",
         *curcliqueweight, *ntreenodes, level);
      for( i = 0; i < *ncurcliquenodes; ++i )
         debugPrintf(" %d", curcliquenodes[i]);
      debugPrintf("\n");
#endif
   }

   /* discard subtree, if the upper bound is not better than the weight of the currently best clique;
    * if only 2 nodes are left, the maximal weighted clique was already calculated in boundSubgraph() and nothing
    * more has to be done;
    * however, if the user wanted to have a fixed node and we are in the first decision level, we have to continue
    */
   if( weightK + subgraphweight > *maxcliqueweight && (nV > 2 || (fixednode >= 0 && level == 0)) )
   {
      int* Vcurrent;
      int nVcurrent;
      int nValive;
      int branchingnode;

      assert(nV > 0);

      /* process current subtree */
      level++;

      /* set up data structures */
      ALLOC_ABORT( BMSallocMemoryArray(&Vcurrent, nV-1) );

      nValive = nV;
      weightKold = weightK;

      debugMessage("============================ branching level %d ===============================\n", level);

      /* branch on the nodes of V by decreasing order of their apriori bound */
      while( backtracklevel >= level && nValive > 0 )
      {
         int branchidx;

         /* check if we meet the backtracking frequency - in this case abort the search until we have reached first level */
         if( level > 1 && backtrackfreq > 0 && (*ntreenodes) % backtrackfreq == 0 )
         {
            backtracklevel = 1;
            break;
         }

         /* get next branching node */
         if( level == 1 && fixednode >= 0 )
         {
            /* select the fixed node as first "branching" candidate */
            for( branchidx = 0; branchidx < nValive && V[branchidx] != fixednode; branchidx++ )
            {}
            assert(branchidx < nValive);
            assert(V[branchidx] == fixednode);
         }
         else if( level == 1 && maxfirstnodeweight > 0 )
            branchidx = getMaxApBoundIndexNotMaxWeight(V, nValive, apbound, weights, maxfirstnodeweight);
         else
            branchidx = getMaxApBoundIndex(nValive, apbound);
         if( branchidx < 0 )
            break;
         assert(0 <= branchidx && branchidx < nValive && nValive <= nV);
         assert(apbound[branchidx] > 0);
         assert(weights[V[branchidx]] > 0);

         /* test a priori bound */
         if( (weightKold + apbound[branchidx]) <= *maxcliqueweight )
            break;

         debugMessage("%d. branching in level %d (treenode %d): bidx=%d, node %d, weight %d, upperbound: %d+%d = %d, maxclique=%d\n",
            nV-nValive+1, level, *ntreenodes, branchidx, V[branchidx], weights[V[branchidx]], weightKold,
            apbound[branchidx], weightKold + apbound[branchidx], *maxcliqueweight);

         /* because we branch on this node, the node is no leaf in the tree */
         isleaf = FALSE;

         /* update the set of nodes from the b&b tree
          *   K = K & {branchingnode}
          */
         branchingnode = V[branchidx];
         K[level-1] = branchingnode;
         weightK = weightKold + weights[branchingnode];

         /* update the set of nodes for branching
          *   V = V \ {branchingnode}
          */
         nValive--;
         for( i = branchidx; i < nValive; ++i )
         {
            V[i] = V[i+1];
            apbound[i] = apbound[i+1];
         }

         /* set the nodes for the next level of b&b tree
          *   Vcurrent = nodes of V, that are adjacent to branchingnode
          */
         nVcurrent = selectadjnodes(tcliquegraph, branchingnode, V, nValive, Vcurrent);

         /* process the selected subtree */
         backtracklevel = branch(getnnodes, getweights, isedge, selectadjnodes, tcliquegraph, newsol, tcliquedata,
            mem, cliquehash, buffer,
            level, Vcurrent, nVcurrent, Vzero, nVzero, gsd, iscolored, K, weightK,
            maxcliquenodes, nmaxcliquenodes, maxcliqueweight,
            curcliquenodes, ncurcliquenodes, curcliqueweight, tmpcliquenodes,
            maxfirstnodeweight, ntreenodes, maxntreenodes, backtrackfreq, maxnzeroextensions, -1, status);

         /* if all other candidates stayed in the candidate list, the current branching was optimal and
          * there is no need to try the remaining ones
          */
         if( nVcurrent == nValive )
         {
            debugMessage("branching on node %d was optimal - ignoring remaining candidates\n", branchingnode);
            nValive = 0;
         }

         /* if we had a fixed node, ignore all other nodes */
         if( fixednode >= 0 )
            nValive = 0;
      }

      debugMessage("========================== branching level %d end =============================\n\n", level);

      /* free data structures */
      BMSfreeMemoryArray(&Vcurrent);
   }

   /* check, whether any branchings have been applied, or if this node is a leaf of the branching tree */
   if( isleaf )
   {
      /* the current clique is the best clique found on the path to this leaf
       * -> check, whether it is an improvement to the currently best clique
       */
      if( *curcliqueweight > *maxcliqueweight )
      {
         TCLIQUE_Bool stopsolving;

         debugMessage("found clique of weight %d at node %d in level %d\n", *curcliqueweight, *ntreenodes, level);
         newSolution(selectadjnodes, tcliquegraph, newsol, tcliquedata, cliquehash, buffer, Vzero, nVzero,
            maxnzeroextensions, curcliquenodes, *ncurcliquenodes, *curcliqueweight,
            maxcliquenodes, nmaxcliquenodes, maxcliqueweight, &stopsolving);

         if( stopsolving )
         {
            debugMessage(" -> solving terminated by callback method\n");
            backtracklevel = 0;
         }
      }

      /* discard the current clique */
      *ncurcliquenodes = 0;
      *curcliqueweight = 0;
   }

#ifdef TCLIQUE_DEBUG
   if( level > backtracklevel )
   {
      debugMessage("premature backtracking after %d nodes - level %d\n", *ntreenodes, level);
   }
#endif

   /* free data structures */
   BMSfreeMemoryArray(&apbound);

   return backtracklevel;
}

/** finds maximum weight clique */
void tcliqueMaxClique(
   TCLIQUE_GETNNODES((*getnnodes)),          /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_ISEDGE((*isedge)),                /**< user function to check for existence of an edge */
   TCLIQUE_SELECTADJNODES((*selectadjnodes)),/**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure that is passed to graph callbacks */
   TCLIQUE_NEWSOL((*newsol)),                /**< user function to call on every new solution */
   TCLIQUE_DATA*         tcliquedata,        /**< user data to pass to new solution callback function */
   int*                  maxcliquenodes,     /**< pointer to store nodes of the maximum weight clique */
   int*                  nmaxcliquenodes,    /**< pointer to store number of nodes in the maximum weight clique */
   TCLIQUE_WEIGHT*       maxcliqueweight,    /**< pointer to store weight of the maximum weight clique */
   TCLIQUE_WEIGHT        maxfirstnodeweight, /**< maximum weight of branching nodes in level 0; 0 if not used
                                              *   for cliques with at least one fractional node) */
   TCLIQUE_WEIGHT        minweight,          /**< lower bound for weight of generated cliques */
   int                   maxntreenodes,      /**< maximal number of nodes of b&b tree */
   int                   backtrackfreq,      /**< frequency to backtrack to first level of tree (0: no premature backtracking) */
   int                   maxnzeroextensions, /**< maximal number of zero-valued variables extending the clique */
   int                   fixednode,          /**< node that is forced to be in the clique, or -1; must have positive weight */
   int*                  ntreenodes,         /**< pointer to store the number of used tree nodes (or NULL) */
   TCLIQUE_STATUS*       status              /**< pointer to store the status of the solving call */
   )
{
   CLIQUEHASH* cliquehash;
   const TCLIQUE_WEIGHT* weights;
   int* buffer;
   int* K;
   int* V;
   int* Vzero;
   int nnodes;
   int nV;
   int nVzero;
   int i;
   BMS_CHKMEM* mem;
   NBC* gsd;
   TCLIQUE_Bool* iscolored;
   int* curcliquenodes;
   int ncurcliquenodes;
   TCLIQUE_WEIGHT curcliqueweight;
   int* tmpcliquenodes;
   int nbbtreenodes;
   int backtracklevel;

   assert(maxcliquenodes != NULL);
   assert(nmaxcliquenodes != NULL);
   assert(maxcliqueweight != NULL);
   assert(maxntreenodes >= 0);
   assert(backtrackfreq >= 0);
   assert(maxnzeroextensions >= 0);
   assert(status != NULL);

   *status = TCLIQUE_OPTIMAL;

   /* use default graph callbacks, if NULL pointers are given */
   if( getnnodes == NULL )
      getnnodes = tcliqueGetNNodes;
   if( getweights == NULL )
      getweights = tcliqueGetWeights;
   if( isedge == NULL )
      isedge = tcliqueIsEdge;
   if( selectadjnodes == NULL )
      selectadjnodes = tcliqueSelectAdjnodes;

   /* get number of nodes */
   nnodes = getnnodes(tcliquegraph);

   debugMessage("calculating maximal weighted clique in graph (%d nodes)\n", nnodes);

   /* set up data structures */
   if( newsol != NULL )
      createCliquehash(&cliquehash, CLIQUEHASH_INITSIZE);
   else
      cliquehash = NULL;
   ALLOC_ABORT( BMSallocMemoryArray(&buffer, nnodes) );
   ALLOC_ABORT( BMSallocMemoryArray(&K, nnodes) );
   ALLOC_ABORT( BMSallocMemoryArray(&V, nnodes) );
   ALLOC_ABORT( BMSallocMemoryArray(&Vzero, nnodes) );
   ALLOC_ABORT( BMSallocMemoryArray(&gsd, nnodes) );
   ALLOC_ABORT( BMSallocMemoryArray(&iscolored, nnodes) );
   ALLOC_ABORT( BMSallocMemoryArray(&curcliquenodes, nnodes) );
   ALLOC_ABORT( BMSallocMemoryArray(&tmpcliquenodes, nnodes) );

   /* set weight and number of nodes of maximum weighted clique */
   *nmaxcliquenodes = 0;
   *maxcliqueweight = minweight-1;
   ncurcliquenodes = 0;
   curcliqueweight = 0;
   nbbtreenodes = 0;

   /* set up V and Vzero */
   weights = getweights(tcliquegraph);
   assert(weights != NULL);
   nV = 0;
   nVzero = 0;
   for( i = 0 ; i <  nnodes; i++ )
   {
      if( weights[i] == 0 )
      {
         Vzero[nVzero] = i;
         nVzero++;
      }
      else
      {
         V[nV] = i;
         nV++;
      }
   }

   /* initialize own memory allocator for coloring */
   mem = BMScreateChunkMemory(sizeof(LIST_ITV), CHUNK_SIZE, -1);

   /* branch to find maximum weight clique */
   backtracklevel = branch(getnnodes, getweights, isedge, selectadjnodes, tcliquegraph, newsol, tcliquedata, mem,
      cliquehash, buffer, 0, V, nV, Vzero, nVzero, gsd, iscolored, K, 0,
      maxcliquenodes, nmaxcliquenodes, maxcliqueweight,
      curcliquenodes, &ncurcliquenodes, &curcliqueweight, tmpcliquenodes,
      maxfirstnodeweight, &nbbtreenodes, maxntreenodes, backtrackfreq, maxnzeroextensions, fixednode, status);

   if ( ntreenodes != NULL )
      *ntreenodes = nbbtreenodes;

   if( backtracklevel != INT_MAX && *status == TCLIQUE_OPTIMAL )
      *status = TCLIQUE_USERABORT;

   /* delete own memory allocator for coloring */
   BMSdestroyChunkMemory(&mem);

   /* free data structures */
   BMSfreeMemoryArray(&tmpcliquenodes);
   BMSfreeMemoryArray(&curcliquenodes);
   BMSfreeMemoryArray(&iscolored);
   BMSfreeMemoryArray(&gsd);
   BMSfreeMemoryArray(&Vzero);
   BMSfreeMemoryArray(&V);
   BMSfreeMemoryArray(&K);
   BMSfreeMemoryArray(&buffer);
   if( newsol != NULL )
      freeCliquehash(&cliquehash);
}
