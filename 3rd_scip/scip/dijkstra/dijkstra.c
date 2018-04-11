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

/**@file   dijkstra.c
 * @brief  C implementation of Dijkstra's algorithm
 * @author Thorsten Koch
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "dijkstra.h"


/** Check whether the data structures of the graph are valid. */
DIJKSTRA_Bool dijkstraGraphIsValid(
   const DIJKSTRA_GRAPH* G                   /**< directed graph to be checked */
   )
{
   unsigned int count = 0;
   unsigned int i;
   unsigned int k;

   if ( G == NULL || G->outbeg == NULL || G->outcnt == NULL || G->weight == NULL || G->head == NULL )
      abort();

   for (i = 0; i < G->nodes; ++i)
   {
      for (k = G->outbeg[i]; k < G->outbeg[i] + G->outcnt[i]; ++k)
      {
         if ( G->head[k] >= G->nodes )
            abort();

         if ( G->weight[k] > G->maxweight || G->weight[k] < G->minweight )
            abort();

         ++count;
      }
      if ( G->head[k] != DIJKSTRA_UNUSED )
         abort();

      ++count;
   }
   if ( count > G->arcs )
      abort();

   return TRUE;
}


#ifndef NDEBUG
/** Check whether heap is valid.
 *
 *  @note Sift up/down do not use order, only for the last the changed one is entered.
 */
static
DIJKSTRA_Bool dijkstraHeapIsValid(
   const unsigned int*   entry,              /**< entries of heap */
   const unsigned long long* value,          /**< values in heap */
   const unsigned int*   order,              /**< order of entries */
   const unsigned int    used,               /**< number of used entries */
   const unsigned int    size                /**< size of entry array */
   )
{
   unsigned int i;

   if ( entry == NULL || value == NULL || order == NULL || used  >  size )
      return FALSE;

   /* check heap property */
   for (i = 0; i < used / 2; ++i)
   {
      if ( value[entry[i]] > value[entry[i + i]] )
         return FALSE;
      if ( i + i + 1 < used && value[entry[i]] > value[entry[i + i + 1]] )
         return FALSE;
   }

   return TRUE;
}
#endif


/** Moves an entry down in the vector until the sorting is valid again. */
static
void dijkstraSiftDown(
   unsigned int*         entry,              /**< entries of heap */
   const unsigned long long* value,          /**< values in heap */
   unsigned int*         order,              /**< order of entries */
   unsigned int          used,               /**< number of used entries */
   unsigned int          current             /**< current entry to be sifted */
   )
{
   unsigned long long val;
   unsigned int child;
   unsigned int ent;
   unsigned int e;

   child = current + current;
   ent = entry[current];
   val = value[ent];

   while ( child < used )
   {
      e = entry[child];

      if ( child + 1 < used )
      {
         if ( value[entry[child + 1]] < value[e] )
         {
            ++child;
            e = entry[child];
         }
      }
      if ( value[e] >= val )
         break;

      entry[current] = e;
      order[e] = current;

      current = child;
      child += child;
   }
   entry[current] = ent;
   order[ent] = current;
}


/** Moves an entry up in the vector until the sorting is valid again. */
static
void dijkstraSiftUp(
   unsigned int*         entry,              /**< entries of heap */
   const unsigned long long* value,          /**< values in heap */
   unsigned int*         order,              /**< order of entries */
   unsigned int          current             /**< current entry to be sifted */
   )
{
   unsigned long long val;
   unsigned int parent;
   unsigned int ent;
   unsigned int e;

   ent = entry[current];
   val = value[ent];

   while ( current > 0 )
   {
      parent = current / 2;
      e = entry[parent];

      if ( value[e] <= val )
         break;

      entry[current] = e;
      order[e] = current;
      current = parent;
   }
   entry[current] = ent;
   order[ent] = current;
}


/** Dijkstra's algorithm for shortest paths from a single source using binary heaps */
unsigned int dijkstra(
   const DIJKSTRA_GRAPH* G,                  /**< directed graph */
   unsigned int          source,             /**< source node */
   unsigned long long*   dist,               /**< node distances (allocated by user) */
   unsigned int*         pred,               /**< node predecessors in final shortest path tree (allocated by user) */
   unsigned int*         entry,              /**< temporary storage (for each node - must be allocated by user) */
   unsigned int*         order               /**< temporary storage (for each node - must be allocated by user) */
   )
{
   unsigned long long weight;
   unsigned int iters = 0;
   unsigned int used = 0;
   unsigned int head;
   unsigned int tail;
   unsigned int i;
   unsigned int e;

   assert( dijkstraGraphIsValid(G) );
   assert( source < G->nodes );
   assert( dist != NULL );
   assert( pred != NULL );

   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   /* initialize nodes */
   for (i = 0; i < G->nodes; ++i)
   {
      dist[i] = DIJKSTRA_FARAWAY;
      order[i] = DIJKSTRA_UNUSED;
      pred[i] = DIJKSTRA_UNUSED;
   }

   /* enter source node into heap */
   entry[0] = source;
   order[source] = 0;
   pred[source] = DIJKSTRA_UNUSED;
   dist[source] = 0;

   ++used;

   /* loop while heap is not empty */
   while ( used > 0 )
   {
      /* get next node */
      tail = entry[0];

      /* remove node from heap */
      --used;
      entry[0] = entry[used];
      order[entry[0]] = 0;
      order[tail] = DIJKSTRA_UNUSED;

      dijkstraSiftDown(entry, dist, order, used, 0);

      assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
      assert( entry[used] < G->nodes );

      /* check adjacent nodes */
      for (e = G->outbeg[tail]; G->head[e] != DIJKSTRA_UNUSED; ++e)
      {
         head = G->head[e];
         weight = G->weight[e] + dist[tail];

	 /* Can we improve the current shortest path? */
         if ( dist[head] > weight )
         {
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
            assert( used < G->nodes );
            assert( head <= G->nodes );

            pred[head] = tail;
            dist[head] = weight;

            if ( order[head] == DIJKSTRA_UNUSED )
            {
               assert( head < G->nodes );

               entry[used] = head;
               order[head] = used;

               dijkstraSiftUp(entry, dist, order, used);
               ++used;
            }
            else
            {
               dijkstraSiftUp(entry, dist, order, order[head]);
            }
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

            ++iters;
         }
      }
   }
   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   return iters;
}


/** Dijkstra's algorithm for shortest paths between a pair of nodes using binary heaps */
unsigned int dijkstraPair(
   const DIJKSTRA_GRAPH* G,                  /**< directed graph */
   unsigned int          source,             /**< source node */
   unsigned int          target,             /**< target node */
   unsigned long long*   dist,               /**< node distances (allocated by user) */
   unsigned int*         pred,               /**< node predecessors in final shortest path tree (allocated by user) */
   unsigned int*         entry,              /**< temporary storage (for each node - must be allocated by user) */
   unsigned int*         order               /**< temporary storage (for each node - must be allocated by user) */
   )
{
   unsigned long long weight;
   unsigned int iters = 0;
   unsigned int used = 0;
   unsigned int head;
   unsigned int tail;
   unsigned int i;
   unsigned int e;

   assert( dijkstraGraphIsValid(G) );
   assert( source < G->nodes );
   assert( target < G->nodes );
   assert( dist != NULL );
   assert( pred != NULL );

   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   /* initialize nodes */
   for (i = 0; i < G->nodes; ++i)
   {
      dist[i] = DIJKSTRA_FARAWAY;
      order[i] = DIJKSTRA_UNUSED;
      pred[i] = DIJKSTRA_UNUSED;
   }

   /* enter source node into heap */
   entry[0] = source;
   order[source] = 0;
   pred[source] = DIJKSTRA_UNUSED;
   dist[source] = 0;

   ++used;

   /* loop while heap is not empty */
   while ( used > 0 )
   {
      /* get next node */
      tail = entry[0];

      /* stop if we have found the target node */
      if ( tail == target )
	 break;

      /* remove node from heap */
      --used;
      entry[0] = entry[used];
      order[entry[0]] = 0;
      order[tail] = DIJKSTRA_UNUSED;

      dijkstraSiftDown(entry, dist, order, used, 0);

      assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
      assert( entry[used] < G->nodes );

      /* check adjacent nodes */
      for (e = G->outbeg[tail]; G->head[e] != DIJKSTRA_UNUSED; ++e)
      {
         head = G->head[e];
         weight = G->weight[e] + dist[tail];

	 /* Can we improve the current shortest path? */
         if ( dist[head] > weight )
         {
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
            assert( used < G->nodes );
            assert( head <= G->nodes );

            pred[head] = tail;
            dist[head] = weight;

            if ( order[head] == DIJKSTRA_UNUSED )
            {
               assert( head < G->nodes );

               entry[used] = head;
               order[head] = used;

               dijkstraSiftUp(entry, dist, order, used);
               ++used;
            }
            else
            {
               dijkstraSiftUp(entry, dist, order, order[head]);
            }
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

            ++iters;
         }
      }
   }
   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   return iters;
}


/** Dijkstra's algorithm for shortest paths between a pair of nodes using binary heaps and truncated at cutoff */
unsigned int dijkstraPairCutoff(
   const DIJKSTRA_GRAPH* G,                  /**< directed graph */
   unsigned int          source,             /**< source node */
   unsigned int          target,             /**< target node */
   unsigned long long    cutoff,             /**< if the distance of a node reached this value, we truncate the search */
   unsigned long long*   dist,               /**< node distances (allocated by user) */
   unsigned int*         pred,               /**< node predecessors in final shortest path tree (allocated by user) */
   unsigned int*         entry,              /**< temporary storage (for each node - must be allocated by user) */
   unsigned int*         order               /**< temporary storage (for each node - must be allocated by user) */
   )
{
   unsigned long long weight;
   unsigned int iters = 0;
   unsigned int used = 0;
   unsigned int head;
   unsigned int tail;
   unsigned int i;
   unsigned int e;

   assert( dijkstraGraphIsValid(G) );
   assert( source < G->nodes );
   assert( target < G->nodes );
   assert( dist != NULL );
   assert( pred != NULL );

   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   /* initialize nodes */
   for (i = 0; i < G->nodes; ++i)
   {
      dist[i] = DIJKSTRA_FARAWAY;
      order[i] = DIJKSTRA_UNUSED;
      pred[i] = DIJKSTRA_UNUSED;
   }

   /* enter source node into heap */
   entry[0] = source;
   order[source] = 0;
   pred[source] = DIJKSTRA_UNUSED;
   dist[source] = 0;

   ++used;

   /* loop while heap is not empty */
   while ( used > 0 )
   {
      /* get next node */
      tail = entry[0];

      /* stop if we have found the target node */
      if ( tail == target )
	 break;

      /* remove node from heap */
      --used;
      entry[0] = entry[used];
      order[entry[0]] = 0;
      order[tail] = DIJKSTRA_UNUSED;

      dijkstraSiftDown(entry, dist, order, used, 0);

      assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
      assert( entry[used] < G->nodes );

      /* only work on nodes if their distance is less than the cutoff */
      if ( dist[tail] >= cutoff )
         continue;

      /* check adjacent nodes */
      for (e = G->outbeg[tail]; G->head[e] != DIJKSTRA_UNUSED; ++e)
      {
         head = G->head[e];
         weight = G->weight[e] + dist[tail];

	 /* Can we improve the current shortest path? */
         if ( dist[head] > weight )
         {
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
            assert( used < G->nodes );
            assert( head <= G->nodes );

            pred[head] = tail;
            dist[head] = weight;

            if ( order[head] == DIJKSTRA_UNUSED )
            {
               assert( head < G->nodes );

               entry[used] = head;
               order[head] = used;

               dijkstraSiftUp(entry, dist, order, used);
               ++used;
            }
            else
            {
               dijkstraSiftUp(entry, dist, order, order[head]);
            }
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

            ++iters;
         }
      }
   }
   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   return iters;
}


/** Dijkstra's algorithm for shortest paths between a pair of nodes ignoring nodes, using binary heaps, and truncated at cutoff */
unsigned int dijkstraPairCutoffIgnore(
   const DIJKSTRA_GRAPH* G,                  /**< directed graph */
   unsigned int          source,             /**< source node */
   unsigned int          target,             /**< target node */
   unsigned int*         ignore,             /**< marking nodes to be ignored (if value is nonzero) */
   unsigned long long    cutoff,             /**< if the distance of a node reached this value, we truncate the search */
   unsigned long long*   dist,               /**< node distances (allocated by user) */
   unsigned int*         pred,               /**< node predecessors in final shortest path tree (allocated by user) */
   unsigned int*         entry,              /**< temporary storage (for each node - must be allocated by user) */
   unsigned int*         order               /**< temporary storage (for each node - must be allocated by user) */
   )
{
   unsigned long long weight;
   unsigned int iters = 0;
   unsigned int used = 0;
   unsigned int head;
   unsigned int tail;
   unsigned int i;
   unsigned int e;

   assert( dijkstraGraphIsValid(G) );
   assert( source < G->nodes );
   assert( target < G->nodes );
   assert( dist != NULL );
   assert( pred != NULL );
   assert( ignore[source] == 0 );
   assert( ignore[target] == 0 );

   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   /* initialize nodes */
   for (i = 0; i < G->nodes; ++i)
   {
      dist[i] = DIJKSTRA_FARAWAY;
      order[i] = DIJKSTRA_UNUSED;
      pred[i] = DIJKSTRA_UNUSED;
   }

   /* enter source node into heap */
   entry[0] = source;
   order[source] = 0;
   pred[source] = DIJKSTRA_UNUSED;
   dist[source] = 0;

   ++used;

   /* loop while heap is not empty */
   while ( used > 0 )
   {
      /* get next node */
      tail = entry[0];

      /* stop if we have found the target node */
      if ( tail == target )
	 break;

      /* remove node from heap */
      --used;
      entry[0] = entry[used];
      order[entry[0]] = 0;
      order[tail] = DIJKSTRA_UNUSED;

      dijkstraSiftDown(entry, dist, order, used, 0);

      assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
      assert( entry[used] < G->nodes );

      /* only work on nodes if their distance is less than the cutoff */
      if ( dist[tail] >= cutoff )
         continue;

      /* check adjacent nodes */
      for (e = G->outbeg[tail]; G->head[e] != DIJKSTRA_UNUSED; ++e)
      {
         head = G->head[e];

         /* skip ignored nodes */
         if ( ignore[head] != 0 )
            continue;

         weight = G->weight[e] + dist[tail];

	 /* Can we improve the current shortest path? */
         if ( dist[head] > weight )
         {
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );
            assert( used < G->nodes );
            assert( head <= G->nodes );

            pred[head] = tail;
            dist[head] = weight;

            if ( order[head] == DIJKSTRA_UNUSED )
            {
               assert( head < G->nodes );

               entry[used] = head;
               order[head] = used;

               dijkstraSiftUp(entry, dist, order, used);
               ++used;
            }
            else
            {
               dijkstraSiftUp(entry, dist, order, order[head]);
            }
            assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

            ++iters;
         }
      }
   }
   assert( dijkstraHeapIsValid(entry, dist, order, used, G->nodes) );

   return iters;
}
