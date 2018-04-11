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

/**@file   tclique_graph.c
 * @brief  graph data part of algorithm for maximum cliques
 * @author Tobias Achterberg
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "tclique/tclique.h"
#include "tclique/tclique_def.h"
#include "blockmemshell/memory.h"


typedef struct _HEAD_ADJ
{
   int              first;
   int              last;
} HEAD_ADJ;

struct TCLIQUE_Graph
{
   int                   nnodes;             /**< number of nodes in graph */
   int                   nedges;             /**< number of edges in graph */
   TCLIQUE_WEIGHT*       weights;            /**< weight of nodes */
   int*                  degrees;            /**< degree of nodes */
   int*                  adjnodes;           /**< adjacent nodes of edges */
   HEAD_ADJ*             adjedges;           /**< pointer to first and one after last adjacent edge of nodes */
   int                   sizenodes;          /**< size of arrays concerning nodes (weights, degrees and adjedges) */
   int                   sizeedges;          /**< size of arrays concerning edges (adjnodes) */
   int*                  cacheddegrees;      /**< number of adjacent cached edges for each node */
   int*                  cachedorigs;        /**< origin nodes of cached edges */
   int*                  cacheddests;        /**< destination nodes of cached edges */
   int                   ncachededges;       /**< number of cached edges (not yet inserted in all data structures) */
   int                   sizecachededges;    /**< size of arrays concerning cached edges */
}; 




/*
 * Interface Methods used by the TClique algorithm
 */

/** gets number of nodes in the graph */
TCLIQUE_GETNNODES(tcliqueGetNNodes)
{
   assert(tcliquegraph != NULL);

   return tcliquegraph->nnodes;
}

/** gets weight of nodes in the graph */
TCLIQUE_GETWEIGHTS(tcliqueGetWeights)
{
   assert(tcliquegraph != NULL);

   return tcliquegraph->weights;
}

/** returns, whether the edge (node1, node2) is in the graph */
TCLIQUE_ISEDGE(tcliqueIsEdge)
{
   int* currentadjedge;
   int* lastadjedge;
   int tmp;

   assert(tcliquegraph != NULL);
   assert(tcliquegraph->ncachededges == 0);
   assert(0 <= node1 && node1 < tcliquegraph->nnodes);
   assert(0 <= node2 && node2 < tcliquegraph->nnodes);

   if( node1 < node2 )
   {
      tmp = node1;
      node1 = node2;
      node2 = tmp;
   }

   currentadjedge = tcliqueGetFirstAdjedge(tcliquegraph, node1);
   lastadjedge = tcliqueGetLastAdjedge(tcliquegraph, node1);

   if( currentadjedge > lastadjedge || *lastadjedge < node2 )
      return FALSE;

   /* checks if node2 is contained in adjacency list of node1 
    * (list is ordered by adjacent nodes) */
   while( currentadjedge <= lastadjedge ) 
   {
      if( *currentadjedge >= node2 )
      {
         if( *currentadjedge == node2 )
            return TRUE;
         else 
            break;
      }
      currentadjedge++;
   }

   return FALSE;
}

/** selects all nodes from a given set of nodes which are adjacent to a given node
 * and returns the number of selected nodes */
TCLIQUE_SELECTADJNODES(tcliqueSelectAdjnodes)
{
   int nadjnodes;
   int* currentadjedge;
   int* lastadjedge;
   int i;

   assert(tcliquegraph != NULL);
   assert(tcliquegraph->ncachededges == 0);
   assert(0 <= node && node < tcliquegraph->nnodes);
   assert(nnodes == 0 || nodes != NULL);
   assert(adjnodes != NULL);

   nadjnodes = 0;
   currentadjedge = tcliqueGetFirstAdjedge(tcliquegraph, node);
   lastadjedge = tcliqueGetLastAdjedge(tcliquegraph, node);

   /* checks for each node in given set nodes, if it is adjacent to given node 
    * (adjacent nodes are ordered by node index)
    */
   for( i = 0; i < nnodes; i++ )
   {
      assert(0 <= nodes[i] && nodes[i] < tcliquegraph->nnodes);
      assert(i == 0 || nodes[i-1] < nodes[i]);
      for( ; currentadjedge <= lastadjedge; currentadjedge++ )
      {
         if( *currentadjedge >= nodes[i] )
         {
            /* current node is adjacent to given node */
            if( *currentadjedge == nodes[i] )
            {
               adjnodes[nadjnodes] = nodes[i]; 
               nadjnodes++;
            }
            break;
         } 
      }
   }

   return nadjnodes;
}




/*
 * External Interface Methods to access the graph (this can be changed without affecting the TClique algorithm)
 */

/** creates graph data structure */
TCLIQUE_Bool tcliqueCreate(
   TCLIQUE_GRAPH**       tcliquegraph        /**< pointer to store graph data structure */
   )
{
   assert(tcliquegraph != NULL);

   ALLOC_FALSE( BMSallocMemory(tcliquegraph) );

   (*tcliquegraph)->nnodes = 0;
   (*tcliquegraph)->nedges = 0;
   (*tcliquegraph)->weights = NULL;
   (*tcliquegraph)->degrees = NULL;
   (*tcliquegraph)->adjnodes = NULL;
   (*tcliquegraph)->adjedges = NULL;
   (*tcliquegraph)->sizenodes = 0;
   (*tcliquegraph)->sizeedges = 0;
   (*tcliquegraph)->cacheddegrees = NULL;
   (*tcliquegraph)->cachedorigs = NULL;
   (*tcliquegraph)->cacheddests = NULL;
   (*tcliquegraph)->ncachededges = 0;
   (*tcliquegraph)->sizecachededges = 0;

   return TRUE;
}

/** frees graph data structure */
void tcliqueFree(
   TCLIQUE_GRAPH**       tcliquegraph        /**< pointer to graph data structure */
   )
{
   assert(tcliquegraph != NULL);

   if( *tcliquegraph != NULL )
   {
      if ( (*tcliquegraph)->adjedges != NULL )
      {
	 BMSfreeMemoryArray(&(*tcliquegraph)->adjedges);
	 BMSfreeMemoryArray(&(*tcliquegraph)->adjnodes);
	 BMSfreeMemoryArray(&(*tcliquegraph)->degrees);
	 BMSfreeMemoryArray(&(*tcliquegraph)->weights);
      }
      if ( (*tcliquegraph)->cacheddegrees )
      {
	 BMSfreeMemoryArrayNull(&(*tcliquegraph)->cacheddegrees);
	 BMSfreeMemoryArrayNull(&(*tcliquegraph)->cachedorigs);
	 BMSfreeMemoryArrayNull(&(*tcliquegraph)->cacheddests);
      }
      BMSfreeMemory(tcliquegraph);
   }
}

/** ensures, that arrays concerning edges in graph data structure can store at least num entries */
static
TCLIQUE_Bool tcliqueEnsureSizeEdges(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   num                 /**< minimum number of entries concerning edges to store */
   )
{
   assert(tcliquegraph != NULL);

   if( num > tcliquegraph->sizeedges )
   {
      int newsize;

      newsize = 2*tcliquegraph->sizeedges;
      if( newsize < num )
         newsize = num;

      ALLOC_FALSE( BMSreallocMemoryArray(&tcliquegraph->adjnodes, newsize) );
      tcliquegraph->sizeedges = newsize;
   }

   assert(num <= tcliquegraph->sizeedges);

   return TRUE;
}

/** ensures, that arrays concerning cached edges in graph data structure can store at least num entries */
static
TCLIQUE_Bool tcliqueEnsureSizeCachedEdges(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   num                 /**< minimum number of entries concerning cached edges to store */
   )
{
   assert(tcliquegraph != NULL);

   if( num > tcliquegraph->sizecachededges )
   {
      int newsize;

      newsize = 2*tcliquegraph->sizecachededges;
      if( newsize < num )
         newsize = num;

      ALLOC_FALSE( BMSreallocMemoryArray(&tcliquegraph->cachedorigs, newsize) );
      ALLOC_FALSE( BMSreallocMemoryArray(&tcliquegraph->cacheddests, newsize) );
      tcliquegraph->sizecachededges = newsize;
   }

   assert(num <= tcliquegraph->sizecachededges);

   return TRUE;
}

/** ensures, that arrays concerning nodes in graph data structure can store at least num entries */
static
TCLIQUE_Bool tcliqueEnsureSizeNodes(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   num                 /**< minimum number of entries concerning nodes to store */
   )
{
   assert(tcliquegraph != NULL);

   if( !tcliqueEnsureSizeEdges(tcliquegraph, 1) )
      return FALSE;
   assert(tcliquegraph->adjnodes != NULL);

   if( num > tcliquegraph->sizenodes )
   {
      int newsize;
      int i;

      newsize = 2*tcliquegraph->sizenodes;
      if( newsize < num )
         newsize = num;

      ALLOC_FALSE( BMSreallocMemoryArray(&tcliquegraph->weights, newsize) );
      ALLOC_FALSE( BMSreallocMemoryArray(&tcliquegraph->degrees, newsize) );
      ALLOC_FALSE( BMSreallocMemoryArray(&tcliquegraph->adjedges, newsize) );

      for( i = tcliquegraph->sizenodes; i < newsize; i++ )
      {
         tcliquegraph->weights[i] = 0;
         tcliquegraph->degrees[i] = 0;
         tcliquegraph->adjedges[i].first = tcliquegraph->nedges;
         tcliquegraph->adjedges[i].last = tcliquegraph->nedges;
      }

      if( tcliquegraph->ncachededges > 0 )
      {
         assert(tcliquegraph->cacheddegrees != NULL);
         ALLOC_FALSE( BMSreallocMemoryArray(&tcliquegraph->cacheddegrees, newsize) );
         for( i = tcliquegraph->sizenodes; i < newsize; i++ )
            tcliquegraph->cacheddegrees[i] = 0;
      }

      tcliquegraph->sizenodes = newsize;
   }
   assert(num <= tcliquegraph->sizenodes);

   return TRUE;
}


/** adds nodes up to the given node number to graph data structure (intermediate nodes have weight 0) */
TCLIQUE_Bool tcliqueAddNode(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   node,               /**< node number to add */
   TCLIQUE_WEIGHT        weight              /**< weight of node to add */
   )
{
   assert(weight >= 0);

   if( !tcliqueEnsureSizeNodes(tcliquegraph, node + 1) )
      return FALSE;

   tcliquegraph->weights[node] = weight;

   assert(tcliquegraph->degrees[node] == 0);
   assert(tcliquegraph->adjedges[node].first <= tcliquegraph->nedges);
   assert(tcliquegraph->adjedges[node].last == tcliquegraph->adjedges[node].first);
   tcliquegraph->nnodes = MAX(tcliquegraph->nnodes, node+1);

   return TRUE;
}

/** changes weight of node in graph data structure */
void tcliqueChangeWeight(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   node,               /**< node to set new weight */
   TCLIQUE_WEIGHT        weight              /**< new weight of node (allready scaled) */
   )
{
   assert(0 <= node && node < tcliqueGetNNodes(tcliquegraph));
   assert(weight >= 0);

   tcliquegraph->weights[node] = weight;
}

/** adds edge (node1, node2) to graph data structure (node1 and node2 have to be contained in 
 *  graph data structure)
 *
 *  New edges are cached, s.t. the graph data structures are not correct until a call to tcliqueFlush();
 *  you have to make sure, that no double edges are inserted.
 */
TCLIQUE_Bool tcliqueAddEdge(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   int                   node1,              /**< start node of edge to add */
   int                   node2               /**< end node of edge to add */
   )
{
   assert(tcliquegraph != NULL);
   assert(0 <= node1 && node1 < tcliquegraph->nnodes);
   assert(0 <= node2 && node2 < tcliquegraph->nnodes);
   assert(node1 != node2);

   if( !tcliqueEnsureSizeCachedEdges(tcliquegraph, tcliquegraph->ncachededges + 2) )
      return FALSE;

   /* make sure, the array for counting the cached node degrees exists */
   if( tcliquegraph->ncachededges == 0 && tcliquegraph->sizenodes > 0 )
   {
      assert(tcliquegraph->cacheddegrees == NULL);
      ALLOC_FALSE( BMSallocMemoryArray(&tcliquegraph->cacheddegrees, tcliquegraph->sizenodes) );
      BMSclearMemoryArray(tcliquegraph->cacheddegrees, tcliquegraph->sizenodes);
   }
   assert(tcliquegraph->cacheddegrees != NULL);

   /* just remember both new half edges in the cache; the full insertion is done later on demand */
   tcliquegraph->cachedorigs[tcliquegraph->ncachededges] = node1;
   tcliquegraph->cacheddests[tcliquegraph->ncachededges] = node2;
   tcliquegraph->ncachededges++;
   tcliquegraph->cachedorigs[tcliquegraph->ncachededges] = node2;
   tcliquegraph->cacheddests[tcliquegraph->ncachededges] = node1;
   tcliquegraph->ncachededges++;
   tcliquegraph->cacheddegrees[node1]++;
   tcliquegraph->cacheddegrees[node2]++;

   return TRUE;
}

/** inserts all cached edges into the data structures */
TCLIQUE_Bool tcliqueFlush(
   TCLIQUE_GRAPH*        tcliquegraph        /**< graph data structure */
   )
{
   assert(tcliquegraph != NULL);

   /* check, whether there are cached edges */
   if( tcliquegraph->ncachededges > 0 )
   {
      int ninsertedholes;
      int pos;
      int n;
      int i;

      /* reallocate adjnodes array to be able to store all additional edges */
      if( !tcliqueEnsureSizeEdges(tcliquegraph, tcliquegraph->nedges + tcliquegraph->ncachededges) )
         return FALSE;
      assert(tcliquegraph->adjnodes != NULL);
      assert(tcliquegraph->adjedges != NULL);

      /* move the old edges in the adjnodes array, s.t. there is enough free space for the additional edges */
      ninsertedholes = 0;
      pos = tcliquegraph->nedges + tcliquegraph->ncachededges - 1;
      for( n = tcliquegraph->nnodes-1; ; --n ) /* no abort criterion, because at n == 0, the loop is break'ed */
      {
         int olddegree;

         assert(n >= 0);
         assert(tcliquegraph->adjedges[n].last - tcliquegraph->adjedges[n].first == tcliquegraph->degrees[n]);

         /* increase the degree of the node */
         olddegree = tcliquegraph->degrees[n];
         tcliquegraph->degrees[n] += tcliquegraph->cacheddegrees[n];

         /* skip space for new edges */
         pos -= tcliquegraph->cacheddegrees[n];
         ninsertedholes += tcliquegraph->cacheddegrees[n];
         assert(ninsertedholes <= tcliquegraph->ncachededges);
         if( ninsertedholes == tcliquegraph->ncachededges )
            break;
         assert(n > 0);

         /* move old edges */
         for( i = tcliquegraph->adjedges[n].last - 1; i >= tcliquegraph->adjedges[n].first; --i, --pos )
         {
            assert(0 <= i && i < pos && pos < tcliquegraph->nedges + tcliquegraph->ncachededges);
            tcliquegraph->adjnodes[pos] = tcliquegraph->adjnodes[i];
         }

         /* adjust the first and last edge pointers of the node */
         tcliquegraph->adjedges[n].first = pos+1;
         tcliquegraph->adjedges[n].last = pos+1 + olddegree;

         assert(n == tcliquegraph->nnodes-1
            || tcliquegraph->adjedges[n].first + tcliquegraph->degrees[n] == tcliquegraph->adjedges[n+1].first);
      }
      assert(ninsertedholes == tcliquegraph->ncachededges);
      assert(tcliquegraph->adjedges[n].last == pos+1);
#ifndef NDEBUG
      for( --n; n >= 0; --n )
         assert(tcliquegraph->cacheddegrees[n] == 0);
#endif

      /* insert the cached edges into the adjnodes array */
      for( i = 0; i < tcliquegraph->ncachededges; ++i )
      {
         int dest;

         n = tcliquegraph->cachedorigs[i];
         dest = tcliquegraph->cacheddests[i];
         assert(0 <= n && n < tcliquegraph->nnodes);
         assert(0 <= dest && dest < tcliquegraph->nnodes);
         assert(tcliquegraph->adjedges[n].last <= tcliquegraph->nedges + tcliquegraph->ncachededges);
         assert(n == tcliquegraph->nnodes-1 || tcliquegraph->adjedges[n].last <= tcliquegraph->adjedges[n+1].first);
         assert(n == tcliquegraph->nnodes-1
            || tcliquegraph->adjedges[n].first + tcliquegraph->degrees[n] == tcliquegraph->adjedges[n+1].first);

         /* edges of each node must be sorted by increasing destination node number */
         for( pos = tcliquegraph->adjedges[n].last;
              pos > tcliquegraph->adjedges[n].first && dest < tcliquegraph->adjnodes[pos-1]; --pos )
         {
            tcliquegraph->adjnodes[pos] = tcliquegraph->adjnodes[pos-1];
         }
         tcliquegraph->adjnodes[pos] = dest;
         tcliquegraph->adjedges[n].last++;

         assert(n == tcliquegraph->nnodes-1 || tcliquegraph->adjedges[n].last <= tcliquegraph->adjedges[n+1].first);
      }

      /* update the number of edges */
      tcliquegraph->nedges += tcliquegraph->ncachededges;

      /* free the cache */
      BMSfreeMemoryArray(&tcliquegraph->cacheddegrees);
      BMSfreeMemoryArray(&tcliquegraph->cachedorigs);
      BMSfreeMemoryArray(&tcliquegraph->cacheddests);
      tcliquegraph->ncachededges = 0;
      tcliquegraph->sizecachededges = 0;
   }

   /* the cache should now be freed */
   assert(tcliquegraph->ncachededges == 0);
   assert(tcliquegraph->sizecachededges == 0);
   assert(tcliquegraph->cacheddegrees == NULL);
   assert(tcliquegraph->cachedorigs == NULL);
   assert(tcliquegraph->cacheddests == NULL);

#ifndef NDEBUG
   /* check integrity of the data structures */
   {
      int pos;
      int n;

      pos = 0;
      for( n = 0; n < tcliquegraph->nnodes; ++n )
      {
         int i;

         assert(tcliquegraph->adjedges[n].first == pos);
         assert(tcliquegraph->adjedges[n].last == tcliquegraph->adjedges[n].first + tcliquegraph->degrees[n]);

         for( i = tcliquegraph->adjedges[n].first; i < tcliquegraph->adjedges[n].last-1; ++i )
         {
            assert(tcliquegraph->adjnodes[i] < tcliquegraph->adjnodes[i+1]);
         }
         pos = tcliquegraph->adjedges[n].last;
      }
      assert(pos == tcliquegraph->nedges);
   }   
#endif

   return TRUE;
}

/** loads graph data structure from file */
TCLIQUE_Bool tcliqueLoadFile(
   TCLIQUE_GRAPH**       tcliquegraph,       /**< pointer to store graph data structure */
   const char*           filename,           /**< name of file with graph data */
   double                scaleval,           /**< value to scale weights (only integral part of scaled weights is considered) */
   char*                 probname,           /**< buffer to store the name of the problem */
   int                   sizeofprobname      /**< size of buffer to store the name of the problem */
   )
{
   FILE* file;
   double weight;
   int node1;
   int node2;
   int currentnode;
   int i;
   int result;
   char* charresult;
   char* tmp;

   assert(tcliquegraph != NULL);
   assert(scaleval > 0.0);

   /* open file */
   if( (file = fopen(filename, "r")) == NULL )
   {
      if( (file = fopen("default.dat", "r")) == NULL )
      {
         infoMessage("\nCan't open file: %s", filename);
         return FALSE;
      }
   }

   if( !tcliqueCreate(tcliquegraph) )
   {
      fclose(file);
      return FALSE;
   }

   /* set name of problem, copies 'sizeofprobname' characters into probname */
   charresult = fgets(probname, sizeofprobname, file);
   if( charresult == NULL )
   {
      infoMessage("Error while reading probname in file %s", filename);
      fclose(file);
      return FALSE;
   }

   /* allocate temporary memory for skipping rest of problem name */
   BMSallocMemoryArray(&tmp, sizeofprobname +1 );
   if( tmp == NULL )
   {
      infoMessage("[%s:%d] No memory in function call", __FILE__, __LINE__);
      fclose(file);
      return FALSE;
   }

   BMScopyMemoryArray(tmp, probname, sizeofprobname);
   probname[sizeofprobname-1] = '\0';
   tmp[sizeofprobname] = '\0';

   /* continue reading until we reach the end of the problem name */
   while( (int) strlen(tmp) == sizeofprobname && tmp[strlen(tmp)-1] != '\n' )
   {
      charresult = fgets(tmp, sizeofprobname, file);

      if( charresult == NULL )
      {
         infoMessage("Error while reading probname in file %s", filename);
         fclose(file);
         return FALSE;
      }
   }

   /* free temporary memory */
   BMSfreeMemoryArray(&tmp);

   /* set number of nodes and number of edges in graph */
   result = fscanf(file, "%d", &(*tcliquegraph)->nnodes);
   if( result <= 0 )
   {
      infoMessage("Error while reading number of nodes in file %s", filename); 
      fclose(file);
      return FALSE;
   }

   result = fscanf(file, "%d", &(*tcliquegraph)->nedges);
   if( result <= 0 )
   {
      infoMessage("Error while reading number of edges in file %s", filename); 
      fclose(file);
      return FALSE;
   }

   if( (*tcliquegraph)->nnodes < 0 || (*tcliquegraph)->nedges < 0 )
   {
      infoMessage("\nInvalid number of %s (%d) in file: %s", (*tcliquegraph)->nnodes < 0 ? "nodes" : "edges", 
         (*tcliquegraph)->nnodes < 0 ? (*tcliquegraph)->nnodes : (*tcliquegraph)->nedges, filename);
      fclose(file);
      return FALSE;
   }

   /* set data structures for tclique,
    * if an error occured, close the file before returning */
   if( BMSallocMemoryArray(&(*tcliquegraph)->weights, (*tcliquegraph)->nnodes) == NULL )
   {
      infoMessage("Run out of memory while reading file %s", filename); 
      (void) fclose(file);
      return FALSE;
   }

   if( BMSallocMemoryArray(&(*tcliquegraph)->degrees, (*tcliquegraph)->nnodes) == NULL )   
   {
      infoMessage("Run out of memory while reading file %s", filename); 
      (void) fclose(file);
      return FALSE;
   }

   if( BMSallocMemoryArray(&(*tcliquegraph)->adjnodes, (*tcliquegraph)->nedges) == NULL )
   {
      infoMessage("Run out of memory while reading file %s", filename); 
      (void) fclose(file);
      return FALSE;
   }

   if( BMSallocMemoryArray(&(*tcliquegraph)->adjedges, (*tcliquegraph)->nnodes) == NULL )
   {
      infoMessage("Run out of memory while reading file %s", filename); 
      (void) fclose(file);
      return FALSE;
   }

   /* set weights of all nodes (scaled!) */
   for( i = 0; i < (*tcliquegraph)->nnodes; i++ )
   {
      result = fscanf(file, "%lf", &weight);
      if( result <= 0 )
      {
         infoMessage("Error while reading weights of nodes in file %s", filename); 
         fclose(file);
         return FALSE;
      }

      (*tcliquegraph)->weights[i] = (TCLIQUE_WEIGHT)(weight * scaleval);
      assert((*tcliquegraph)->weights[i] >= 0);
   }

   /* set adjacent edges and degree of all nodes */
   currentnode = -1;
   for( i = 0; i < (*tcliquegraph)->nedges; i++ )
   {
      /* read edge (node1, node2) */
      result = fscanf(file, "%d%d", &node1, &node2);
      if( result <= 1 )
      {
         infoMessage("Error while reading edges in file %s", filename); 
         fclose(file);
         return FALSE;
      }

      if( node1 < 0 || node2 < 0 )
      {
         infoMessage("\nInvalid node index (%d) in file: %s", node1 < 0 ? node1 : node2, filename);
         fclose(file);
         return FALSE;
      } 

      /* (node1, node2) is the first adjacent edge of node1 */
      if( node1 != currentnode )
      {
         currentnode = node1;
         (*tcliquegraph)->degrees[currentnode] = 0;
         (*tcliquegraph)->adjedges[currentnode].first = i;
         (*tcliquegraph)->adjedges[currentnode].last = (*tcliquegraph)->adjedges[currentnode].first;
      }
      (*tcliquegraph)->degrees[currentnode]++;
      (*tcliquegraph)->adjnodes[i] = node2;
      (*tcliquegraph)->adjedges[currentnode].last++;
   }

   /* close file */
   fclose(file);

   return TRUE;
}

/** saves graph data structure to file */
TCLIQUE_Bool tcliqueSaveFile(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< graph data structure */
   const char*           filename,           /**< name of file to create */
   double                scaleval,           /**< value to unscale weights with */
   const char*           probname            /**< name of the problem */
   )
{
   FILE* file;
   int i;
   int j;

   assert(tcliquegraph != NULL);
   assert(scaleval > 0.0);

   /* create file */
   if( (file = fopen(filename, "w")) == NULL )
   {
      infoMessage("\nCan't create file: %s", filename);
      return FALSE;
   }

   /* write name of problem, number of nodes and number of edges in graph */
   fprintf(file, "%s\n", probname);
   fprintf(file, "%d\n", tcliquegraph->nnodes);
   fprintf(file, "%d\n", tcliquegraph->nedges);

   /* write weights of all nodes (scaled!) */
   for( i = 0; i < tcliquegraph->nnodes; i++ )
      fprintf(file, "%f\n", (double)tcliquegraph->weights[i]/scaleval);

   /* write edges */
   for( i = 0; i < tcliquegraph->nnodes; i++ )
   {
      for( j = tcliquegraph->adjedges[i].first; j < tcliquegraph->adjedges[i].last; j++ )
         fprintf(file, "%d %d\n", i, tcliquegraph->adjnodes[j]);
   }

   /* close file */
   fclose(file);

   return TRUE;
}

/** gets number of edges in the graph */
int tcliqueGetNEdges(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   )
{
   assert(tcliquegraph != NULL);

   return tcliquegraph->nedges + tcliquegraph->ncachededges;
}

/** gets degree of nodes in graph */
int* tcliqueGetDegrees(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   )
{
   assert(tcliquegraph != NULL);
   assert(tcliquegraph->ncachededges == 0);

   return tcliquegraph->degrees;
}

/** gets adjacent nodes of edges in graph */
int* tcliqueGetAdjnodes(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   )
{
   assert(tcliquegraph != NULL);
   assert(tcliquegraph->ncachededges == 0);

   return tcliquegraph->adjnodes;
}

/** gets pointer to first adjacent edge of given node in graph */
int* tcliqueGetFirstAdjedge(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   int                   node                /**< given node */
   )
{
   HEAD_ADJ* adjedges;
   int* adjnodes;

   assert(tcliquegraph != NULL);
   assert(tcliquegraph->ncachededges == 0);
   assert(0 <= node && node < tcliquegraph->nnodes);

   adjedges = tcliquegraph->adjedges;
   assert(adjedges != NULL);
   assert(adjedges[node].first >= 0);
   assert(adjedges[node].first <= tcliqueGetNEdges(tcliquegraph));

   adjnodes = tcliqueGetAdjnodes(tcliquegraph);
   assert(adjnodes != NULL);

   return &adjnodes[adjedges[node].first];
}

/** gets pointer to last adjacent edge of given node in graph */
int* tcliqueGetLastAdjedge(
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   int                   node                /**< given node */
   )
{
   HEAD_ADJ* adjedges;
   int* adjnodes;
#ifndef NDEBUG
   int* degrees;
#endif

   assert(tcliquegraph != NULL);
   assert(tcliquegraph->ncachededges == 0);
   assert(0 <= node && node < tcliquegraph->nnodes);

   adjedges = tcliquegraph->adjedges;
#ifndef NDEBUG
   degrees = tcliqueGetDegrees(tcliquegraph);
#endif
   assert(adjedges != NULL);
   assert(degrees[node] == 0 || adjedges[node].last-1 >= 0);
   assert(adjedges[node].last-1 <= tcliqueGetNEdges(tcliquegraph));

   assert(adjedges[node].last - adjedges[node].first == degrees[node]);

   adjnodes = tcliqueGetAdjnodes(tcliquegraph);
   assert(adjnodes != NULL);

   return &adjnodes[adjedges[node].last-1];
}

/** prints graph data structure */
void tcliquePrintGraph(
   TCLIQUE_GRAPH*        tcliquegraph        /**< pointer to graph data structure */
   )
{
   const int* weights;
   int* degrees;
   int i;

   assert(tcliquegraph != NULL);
   assert(tcliquegraph->ncachededges == 0);

   degrees = tcliqueGetDegrees(tcliquegraph);
   weights = tcliqueGetWeights(tcliquegraph);

   infoMessage("nnodes=%d, nedges=%d\n", tcliqueGetNNodes(tcliquegraph), tcliqueGetNEdges(tcliquegraph));
   for( i = 0; i < tcliqueGetNNodes(tcliquegraph); i++ )
   {
      int* currentadjedge;
      int* lastadjedge;

      infoMessage("node %d: weight=%d, degree=%d, adjnodes=\n[ ", i, weights[i], degrees[i]);  

      currentadjedge = tcliqueGetFirstAdjedge(tcliquegraph, i);
      lastadjedge = tcliqueGetLastAdjedge(tcliquegraph, i);
      assert(lastadjedge + 1 - currentadjedge == degrees[i]);

      for( ; currentadjedge <= lastadjedge; currentadjedge++ )
      {
	 infoMessage("%d, ", *currentadjedge);
      }
      infoMessage("]\n");
   }
}
