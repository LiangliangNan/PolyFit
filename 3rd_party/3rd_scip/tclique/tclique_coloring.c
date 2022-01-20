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

/**@file   tclique_coloring.c
 * @brief  coloring part of algorithm for maximum cliques
 * @author Tobias Achterberg
 * @author Ralf Borndoerfer
 * @author Zoltan Kormos
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "tclique/tclique.h"
#include "tclique/tclique_def.h"
#include "tclique/tclique_coloring.h"
#include "blockmemshell/memory.h"



/** gets index of the uncolored node in a given array of nodes in V with maximum satdeg;
 *  in case of a tie choose node with maximum weight;
 *  if no uncolored node is found, -1 is returned
 */
static
int getMaxSatdegIndex(
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV,                 /**< number of non-zero weighted nodes for branching */
   NBC*                  gsd,                /**< neighbor color information of all nodes */
   TCLIQUE_Bool*         iscolored,          /**< coloring status of all nodes */
   const TCLIQUE_WEIGHT* weights             /**< weight of nodes in grpah */
   )
{   
   TCLIQUE_WEIGHT maxweight;
   int maxsatdeg;
   int maxsatdegindex;
   int i;

   maxweight = -1;
   maxsatdeg = -1;
   maxsatdegindex = -1;

   assert(gsd != NULL);
   assert(iscolored != NULL);

   for( i = 0; i < nV; i++ )
   {
      TCLIQUE_WEIGHT weight;
      int satdeg;

      /* check only uncolored nodes */ 
      if( iscolored[i] ) 
         continue;

      weight = weights[V[i]];
      assert(weight > 0);

      satdeg = gsd[i].satdeg;
      if( satdeg > maxsatdeg || (satdeg == maxsatdeg && weight > maxweight) )
      {
         maxsatdeg = satdeg;
         maxweight = weight;
         maxsatdegindex = i;
      }
   }

   return maxsatdegindex;	
}

/** gets index of the node in a given set of nodes with maximum weight */
static
int getMaxWeightIndex( 
   TCLIQUE_GETNNODES((*getnnodes)),          /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV                  /**< number of non-zero weighted nodes for branching */
   )
{
   const TCLIQUE_WEIGHT* weights;
   TCLIQUE_WEIGHT maxweight;
   int maxweightindex;
   int i;

   assert(getnnodes != NULL);
   assert(getweights != NULL);
   assert(tcliquegraph != NULL);
   assert(nV > 0);

   weights = getweights(tcliquegraph);

   maxweightindex = -1;
   maxweight = 0;

   /* try to improve maxweight */
   for( i = 0 ; i < nV; i++ )
   {
      assert(0 <= V[i] && V[i] < getnnodes(tcliquegraph));
      assert(weights[V[i]] > 0);
      if( weights[V[i]] > maxweight)
      {
         /* node has larger weight */
         maxweight = weights[V[i]];
         maxweightindex = i;
      }
   }
   assert(maxweightindex >= 0);

   return maxweightindex;
}

/** updates the neighbor colors information of a node: updates the list of neighbor color intervals 
 *  by making the union of the existing list and the given list of color intervals, and updates the saturation degree
 */
static
void updateNeighbor(
   BMS_CHKMEM*           mem,                /**< block memory */
   NBC*                  pgsd,               /**< pointer to neighbor color information of node to update */
   LIST_ITV*             pnc                 /**< pointer to given list of color intervals */
   )
{
   LIST_ITV head;
   LIST_ITV* apciv;
   LIST_ITV* pciv;
   LIST_ITV* nciv;

   /* save the pointer to the first element of the list */
   head.next = pgsd->lcitv;
   apciv = &head;
   pciv = apciv->next;

   /* construct the union of the two intervals */
   while( (pnc != NULL) && (pciv != NULL) )
   {
      if( pnc->itv.inf < pciv->itv.inf ) 
      {	
         ALLOC_ABORT( BMSallocChunkMemory(mem, &nciv) );
         nciv->itv = pnc->itv;
         nciv->next = pciv;
         apciv->next = nciv;
         apciv = nciv;

         pnc = pnc->next;	
      }
      else if( pnc->itv.inf <= pciv->itv.sup )
      {
         if( pnc->itv.sup > pciv->itv.sup )
            pciv->itv.sup = pnc->itv.sup;
         pnc = pnc->next;
      }
      else
      {
         apciv = pciv;
         pciv = pciv->next;
      }
   }

   while( pnc != NULL )
   {
      ALLOC_ABORT( BMSallocChunkMemory(mem, &nciv) );
      nciv->itv = pnc->itv;
      nciv->next = NULL;

      apciv->next = nciv;
      apciv = nciv;

      pnc = pnc->next;
   }

   /* try to reduce the number of intervals */
   pgsd->satdeg = 0;
   apciv = head.next;
   while( (pciv = apciv->next) != NULL ) /*lint !e838*/
   {
      if( apciv->itv.sup < (pciv->itv.inf - 1) )
      {
         pgsd->satdeg += apciv->itv.sup - apciv->itv.inf + 1;
         apciv = apciv->next;
      }
      else
      {
         LIST_ITV* tmp;

         if( apciv->itv.sup < pciv->itv.sup )
            apciv->itv.sup = pciv->itv.sup;
         apciv->next = pciv->next;

         /* free data structure for created colorinterval */
         tmp = pciv->next; 
         BMSfreeChunkMemory(mem, &pciv); 
         pciv = tmp; 
      }
   }
   pgsd->satdeg += apciv->itv.sup - apciv->itv.inf + 1;

   /* updates the pointer to the first element of the list */		
   pgsd->lcitv = head.next;
}

/** colors the positive weighted nodes of a given set of nodes V with the lowest possible number of colors and 
 *  finds a clique in the graph induced by V, an upper bound and an apriori bound for further branching steps
 */
TCLIQUE_WEIGHT tcliqueColoring( 
   TCLIQUE_GETNNODES((*getnnodes)),          /**< user function to get the number of nodes */
   TCLIQUE_GETWEIGHTS((*getweights)),        /**< user function to get the node weights */
   TCLIQUE_SELECTADJNODES((*selectadjnodes)),/**< user function to select adjacent edges */
   TCLIQUE_GRAPH*        tcliquegraph,       /**< pointer to graph data structure */
   BMS_CHKMEM*           mem,                /**< block memory */
   int*                  buffer,             /**< buffer of size nnodes */
   int*                  V,                  /**< non-zero weighted nodes for branching */
   int                   nV,                 /**< number of non-zero weighted nodes for branching */
   NBC*                  gsd,                /**< neighbor color information of all nodes */
   TCLIQUE_Bool*         iscolored,          /**< coloring status of all nodes */
   TCLIQUE_WEIGHT*       apbound,            /**< pointer to store apriori bound of nodes for branching */ 
   int*                  clique,             /**< buffer for storing the clique */
   int*                  nclique,            /**< pointer to store number of nodes in the clique */
   TCLIQUE_WEIGHT*       weightclique        /**< pointer to store the weight of the clique */
   )
{
   const TCLIQUE_WEIGHT* weights;
   TCLIQUE_WEIGHT maxsatdegree; 
   TCLIQUE_WEIGHT range;
   TCLIQUE_Bool growclique; 
   int node; 
   int nodeVindex;
   int i;     
   int j;
   LIST_ITV* colorinterval;
   LIST_ITV nwcitv;
   LIST_ITV* pnc;
   LIST_ITV* lcitv;
   LIST_ITV* item;
   LIST_ITV* tmpitem;
   int* workclique;
   int* currentclique;
   int ncurrentclique;
   int weightcurrentclique;
   int* Vadj;
   int nVadj;
   int adjidx;

   assert(getnnodes != NULL);
   assert(getweights != NULL);
   assert(selectadjnodes != NULL);
   assert(buffer != NULL);
   assert(V != NULL);
   assert(nV > 0);
   assert(clique != NULL);
   assert(nclique != NULL);
   assert(weightclique != NULL);
   assert(gsd != NULL);
   assert(iscolored != NULL);

   weights = getweights(tcliquegraph);
   assert(weights != NULL);

   /* initialize maximum weight clique found so far */
   growclique = TRUE;
   *nclique = 0;
   *weightclique = 0;

   /* get node of V with maximum weight */
   nodeVindex = getMaxWeightIndex(getnnodes, getweights, tcliquegraph, V, nV);
   node = V[nodeVindex];
   assert(0 <= node && node < getnnodes(tcliquegraph));
   range = weights[node];
   assert(range > 0);

   /* set up data structures for coloring */
   BMSclearMemoryArray(iscolored, nV); /* new-memory */
   BMSclearMemoryArray(gsd, nV); /* new-memory */
   iscolored[nodeVindex] = TRUE;

   /* color the first node */
   debugMessage("---------------coloring-----------------\n");
   debugMessage("1. node choosen: vindex=%d, vertex=%d, satdeg=%d, range=%d)\n",
      nodeVindex, node, gsd[nodeVindex].satdeg, range);

   /* set apriori bound: apbound(v_i) = satdeg(v_i) + weight(v_i) */
   apbound[nodeVindex] = range;
   assert(apbound[nodeVindex] > 0);

   /* update maximum saturation degree: maxsatdeg = max { satdeg(v_i) + weight(v_i) | v_i in V } */
   maxsatdegree = range;

   debugMessage("-> updated neighbors:\n");

   /* set neighbor color of the adjacent nodes of node */
   Vadj = buffer;
   nVadj = selectadjnodes(tcliquegraph, node, V, nV, Vadj);
   for( i = 0, adjidx = 0; i < nV && adjidx < nVadj; ++i )
   {
      assert(V[i] <= Vadj[adjidx]); /* Vadj is a subset of V */
      if( V[i] == Vadj[adjidx] )
      {
         /* node is adjacent to itself, but we do not need to color it again */
         if( i == nodeVindex )
         {
            /* go to the next node in Vadj */
            adjidx++;
            continue;
         }

         debugMessage("     nodeVindex=%d, node=%d, weight=%d, satdegold=%d  ->  ", 
            i, V[i], weights[V[i]], gsd[i].satdeg); 

         /* sets satdeg for adjacent node */
         gsd[i].satdeg = range;

         /* creates new color interval [1,range] */
         ALLOC_ABORT( BMSallocChunkMemory(mem, &colorinterval) );
         colorinterval->next = NULL;
         colorinterval->itv.inf = 1;
         colorinterval->itv.sup = range;

         /* colorinterval is the first added element of the list of neighborcolors of the adjacent node  */ 
         gsd[i].lcitv = colorinterval;

         /* go to the next node in Vadj */
         adjidx++;

         debugPrintf("satdegnew=%d, nbc=[%d,%d]\n", gsd[i].satdeg, gsd[i].lcitv->itv.inf, gsd[i].lcitv->itv.sup);
      }
   }

   /* set up data structures for the current clique */
   ALLOC_ABORT( BMSallocMemoryArray(&currentclique, nV) );
   workclique = clique;

   /* add node to the current clique */ 
   currentclique[0] = node; 
   ncurrentclique = 1; 
   weightcurrentclique = range; 

   /* color all other nodes of V */
   for( i = 0 ; i < nV-1; i++ )
   {
      assert((workclique == clique) != (currentclique == clique));

      /* selects the next uncolored node to color */
      nodeVindex = getMaxSatdegIndex(V, nV, gsd, iscolored, weights);
      if( nodeVindex == -1 ) /* no uncolored nodes left */
         break;

      node = V[nodeVindex];
      assert(0 <= node && node < getnnodes(tcliquegraph));
      range = weights[node];
      assert(range > 0);
      iscolored[nodeVindex] = TRUE;	

      debugMessage("%d. node choosen: vindex=%d, vertex=%d, satdeg=%d, range=%d, growclique=%u, weight=%d)\n",
         i+2, nodeVindex, node, gsd[nodeVindex].satdeg, range, growclique, weightcurrentclique);

      /* set apriori bound: apbound(v_i) = satdeg(v_i) + weight(v_i) */
      apbound[nodeVindex] = gsd[nodeVindex].satdeg + range;
      assert(apbound[nodeVindex] > 0);

      /* update maximum saturation degree: maxsatdeg = max { satdeg(v_i) + weight(v_i) | v_i in V } */
      if( maxsatdegree < apbound[nodeVindex] )
         maxsatdegree = apbound[nodeVindex];

      /* update clique */
      if( gsd[nodeVindex].satdeg == 0 )
      {
         /* current node is not adjacent to nodes of current clique, 
          * i.e. current clique can not be increased
          */
         debugMessage("current node not adjacend to current clique (weight:%d) -> starting new clique\n", 
            weightcurrentclique);

         /* check, if weight of current clique is larger than weight of maximum weight clique found so far */ 
         if( weightcurrentclique > *weightclique )
         {
            int* tmp;

            /* update maximum weight clique found so far */
            assert((workclique == clique) != (currentclique == clique));
            tmp = workclique;
            *weightclique = weightcurrentclique;
            *nclique = ncurrentclique;
            workclique = currentclique;
            currentclique = tmp;
            assert((workclique == clique) != (currentclique == clique));
         }
         weightcurrentclique = 0;
         ncurrentclique = 0;
         growclique = TRUE;
      }
      if( growclique )
      {
         /* check, if the current node is still adjacent to all nodes in the clique */
         if( gsd[nodeVindex].satdeg == weightcurrentclique )
         {
            assert(ncurrentclique < nV);
            currentclique[ncurrentclique] = node;
            ncurrentclique++; 
            weightcurrentclique += range;
#ifdef TCLIQUE_DEBUG
            {
               int k;
               debugMessage("current clique (size:%d, weight:%d):", ncurrentclique, weightcurrentclique);
               for( k = 0; k < ncurrentclique; ++k )
                  debugPrintf(" %d", currentclique[k]);
               debugPrintf("\n");
            }
#endif
         }
         else
         {
            debugMessage("node satdeg: %d, clique weight: %d -> stop growing clique\n", 
               gsd[nodeVindex].satdeg, weightcurrentclique);
            growclique = FALSE;
         }
      }

      /* search for fitting color intervals for current node */
      pnc = &nwcitv;
      if( gsd[nodeVindex].lcitv == NULL )
      {
         /* current node has no colored neighbors yet: create new color interval [1,range] */
         ALLOC_ABORT( BMSallocChunkMemory(mem, &colorinterval) );
         colorinterval->next = NULL;
         colorinterval->itv.inf = 1;
         colorinterval->itv.sup = range;

         /* add the new colorinterval [1, range] to the list of chosen colorintervals for node */
         pnc->next = colorinterval;
      }
      else
      {
         int tocolor;
         int dif;

         /* current node has colored neighbors */
         tocolor = range;
         lcitv = gsd[nodeVindex].lcitv;

         /* check, if first neighbor color interval [inf, sup] has inf > 1 */
         if( lcitv->itv.inf != 1 )
         {
            /* create new interval [1, min{range, inf}] */ 
            dif =  lcitv->itv.inf - 1 ;
            if( dif > tocolor )
               dif = tocolor;

            ALLOC_ABORT( BMSallocChunkMemory(mem, &colorinterval) );
            colorinterval->next = NULL;
            colorinterval->itv.inf = 1;
            colorinterval->itv.sup = dif;

            tocolor -= dif;
            pnc->next = colorinterval;
            pnc = colorinterval;
         }

         /* as long as node is not colored with all colors, create new color interval by filling 
          * the gaps in the existing neighbor color intervals of the neighbors of node
          */
         while( tocolor > 0 )
         {	
            dif = tocolor;	

            ALLOC_ABORT( BMSallocChunkMemory(mem, &colorinterval) );
            colorinterval->next = NULL;
            colorinterval->itv.inf = lcitv->itv.sup+1;			
            if( lcitv->next != NULL )
            {
               int min;

               min = lcitv->next->itv.inf - lcitv->itv.sup - 1;

               if( dif > min )  
                  dif = min;	
               lcitv = lcitv->next;
            }
            colorinterval->itv.sup = colorinterval->itv.inf + dif - 1;

            tocolor -= dif;
            pnc->next = colorinterval;
            pnc = colorinterval;
         }
      }

      debugMessage("-> updated neighbors:\n"); 

      /* update saturation degree and neighbor colorintervals of all neighbors of node */
      Vadj = buffer;
      nVadj = selectadjnodes(tcliquegraph, node, V, nV, Vadj);
      for( j = 0, adjidx = 0; j < nV && adjidx < nVadj; ++j )
      {
         assert(V[j] <= Vadj[adjidx]); /* Vadj is a subset of V */
         if( V[j] == Vadj[adjidx] )
         {
            if( !iscolored[j] )
            {
               debugMessage("     nodeVindex=%d, node=%d, weight=%d, satdegold=%d  ->  ", 
                  j, V[j], weights[V[j]], gsd[j].satdeg); 
               updateNeighbor(mem, &gsd[j], nwcitv.next);
               debugPrintf("satdegnew=%d, nbc=[%d,%d]\n", gsd[j].satdeg, gsd[j].lcitv->itv.inf, gsd[j].lcitv->itv.sup);
            }

            /* go to the next node in Vadj */
            adjidx++;
         }
      }

      /* free data structure of created colorintervals */
      item = nwcitv.next;
      while( item != NULL )
      {
         tmpitem = item->next;                  
         BMSfreeChunkMemory(mem, &item);       
         item = tmpitem;                        
      }

      /* free data structure of neighbor colorinterval of node just colored */
      item = gsd[nodeVindex].lcitv;
      while( item != NULL )
      {
         tmpitem = item->next;                  
         BMSfreeChunkMemory(mem, &item);       
         item = tmpitem;                        
      }
   }
   assert((workclique == clique) != (currentclique == clique));

   /* update maximum weight clique found so far */
   if( weightcurrentclique > *weightclique )
   {
      int* tmp;

      tmp = workclique;
      *weightclique = weightcurrentclique;
      *nclique = ncurrentclique;
      workclique = currentclique;
      currentclique = tmp;
   }
   assert((workclique == clique) != (currentclique == clique));

   /* move the found clique to the provided clique pointer, if it is not the memory array */
   if( workclique != clique )
   {
      assert(clique == currentclique);
      assert(*nclique <= nV);
      BMScopyMemoryArray(clique, workclique, *nclique);
      currentclique = workclique;
   }

   /* free data structures */
   BMSfreeMemoryArray(&currentclique);

   /* clear chunk memory */
   BMSclearChunkMemory(mem);

   debugMessage("------------coloringend-----------------\n");

   return maxsatdegree;
}
