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

/**@file   heur_octane.c
 * @brief  octane primal heuristic based on Balas, Ceria, Dawande, Margot, and Pataki
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>
#include "scip/heur_octane.h"

#define HEUR_NAME             "octane"
#define HEUR_DESC             "octane primal heuristic for pure {0;1}-problems based on Balas et al."
#define HEUR_DISPCHAR         'O'
#define HEUR_PRIORITY         -1008000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_FMAX          100       /**< {0,1}-points to be checked */
#define DEFAULT_FFIRST        10        /**< {0,1}-points to be generated at first */
#define DEFAULT_USEFRACSPACE  TRUE      /**< use heuristic for the space of fractional variables or for whole space? */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   int                   f_max;              /**< {0,1}-points to be checked */
   int                   f_first;            /**< {0,1}-points to be generated at first in order to check whether restart is necessary */
   int                   lastrule;           /**< last ray selection rule that was performed */
   SCIP_Bool             usefracspace;       /**< use heuristic for the space of fractional variables or for the whole space? */
   SCIP_Bool             useobjray;          /**< should the inner normal of the objective be used as one ray direction? */
   SCIP_Bool             useavgray;          /**< should the average ray of the basic cone be used as one ray direction? */
   SCIP_Bool             usediffray;         /**< should difference between root sol and current LP sol be used as one ray direction? */
   SCIP_Bool             useavgwgtray;       /**< should the weighted average ray of the basic cone be used as one ray direction? */
   SCIP_Bool             useavgnbray;        /**< should the average ray of the nonbasic cone be used as one ray direction? */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
};

/*
 * Local methods
 */


/** tries to insert the facet obtained from facet i flipped in component j into the list of the fmax nearest facets */
static
void tryToInsert(
   SCIP*                 scip,               /**< SCIP data structure                        */
   SCIP_Bool**           facets,             /**< facets got so far                          */
   SCIP_Real*            lambda,             /**< distances of the facets                    */
   int                   i,                  /**< current facet                              */
   int                   j,                  /**< component to flip                          */
   int                   f_max,              /**< maximal number of facets to create         */
   int                   nsubspacevars,      /**< dimension of the fractional space          */
   SCIP_Real             lam,                /**< distance of the current facet              */
   int*                  nfacets             /**< number of facets                           */
   )
{
   SCIP_Bool* lastfacet;
   int k;

   assert(scip != NULL);
   assert(facets != NULL);
   assert(lambda != NULL);
   assert(nfacets != NULL);

   if( SCIPisFeasLE(scip, lam, 0.0) || SCIPisFeasGE(scip, lam, lambda[f_max-1]) )
      return;

   lastfacet = facets[f_max];

   /* shifting lam through lambda, lambda keeps increasingly sorted */
   for( k = f_max; k > 0 && SCIPisFeasGT(scip, lambda[k-1], lam); --k )
   {
      lambda[k] = lambda[k-1];
      facets[k] = facets[k-1];
   }
   assert(i < k && k < f_max );

   /* inserting new facet into list, new facet is facet at position i flipped in coordinate j, new distance lam */
   facets[k] = lastfacet;
   lambda[k] = lam;

   /*lint --e{866}*/
   BMScopyMemoryArray(facets[k], facets[i], nsubspacevars);
   facets[k][j] = !facets[k][j];
   (*nfacets)++;
}

/** constructs a solution from a given facet paying attention to the transformations made at the beginning of OCTANE */
static
SCIP_RETCODE getSolFromFacet(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Bool*            facet,              /**< current facet                         */
   SCIP_SOL*             sol,                /**< solution to create                    */
   SCIP_Bool*            sign,               /**< marker for retransformation            */
   SCIP_VAR**            subspacevars,       /**< pointer to fractional space variables */
   int                   nsubspacevars       /**< dimension of fractional space         */
   )
{
   int v;

   assert(scip != NULL);
   assert(facet != NULL);
   assert(sol != NULL);
   assert(sign != NULL);
   assert(subspacevars != NULL);

   SCIP_CALL( SCIPlinkLPSol(scip, sol) );
   for( v = nsubspacevars - 1; v >= 0; --v )
   {
      /* after permutation, a variable should be set to 1, iff there was no reflection in this coordinate and the hit
       * facet has coordinate + or there was a reflection and the facet has coordinate - */
      if( facet[v] == sign[v] )
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, subspacevars[v], 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, subspacevars[v], 0.0) );
      }
   }

   return SCIP_OKAY;
}

/** generates the direction of the shooting ray as the inner normal of the objective function */
static
SCIP_RETCODE generateObjectiveRay(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Real*            raydirection,       /**< shooting ray                          */
   SCIP_VAR**            subspacevars,       /**< pointer to fractional space variables */
   int                   nsubspacevars       /**< dimension of fractional space         */
   )
{
   int v;

   assert(scip != NULL);
   assert(raydirection != NULL);
   assert(subspacevars != NULL);

   for( v = nsubspacevars - 1; v >= 0; --v )
      raydirection[v] = SCIPvarGetObj(subspacevars[v]);
   return SCIP_OKAY;
}

/** generates the direction of the shooting ray as the difference between the root and the current LP solution */
static
SCIP_RETCODE generateDifferenceRay(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Real*            raydirection,       /**< shooting ray                          */
   SCIP_VAR**            subspacevars,       /**< pointer to fractional space variables */
   int                   nsubspacevars       /**< dimension of fractional space         */
   )
{
   int v;

   assert(scip != NULL);
   assert(raydirection != NULL);
   assert(subspacevars != NULL);

   for( v = nsubspacevars - 1; v >= 0; --v )
      raydirection[v] = SCIPvarGetLPSol(subspacevars[v]) - SCIPvarGetRootSol(subspacevars[v]);

   return SCIP_OKAY;
}


/** generates the direction of the shooting ray as the average of the extreme rays of the basic cone */
static
SCIP_RETCODE generateAverageRay(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Real*            raydirection,       /**< shooting ray                          */
   SCIP_VAR**            subspacevars,       /**< pointer to fractional space variables */
   int                   nsubspacevars,      /**< dimension of fractional space         */
   SCIP_Bool             weighted            /**< should the rays be weighted?          */
   )
{
   SCIP_ROW** rows;
   SCIP_Real** tableaurows;
   SCIP_Real* rownorm;
   SCIP_Real rowweight;
   int** tableaurowinds;                     /* indices of non-zero entries */
   int* ntableaurowinds;                     /* number of non-zero entries */
   SCIP_Bool* usedrowids = NULL;             /* visited row indices */
   int* rowids;                              /* row indices */
   int nrowids = 0;                          /* number of row indices */
   int tableaurowind;
   int nrows;
   int i;
   int j;
   int sparse = -1; /* used to check that either all information is sparse or not sparse */

   assert(scip != NULL);
   assert(raydirection != NULL);
   assert(subspacevars != NULL);
   assert(nsubspacevars > 0);

   /* get data */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert(nrows > 0);
   assert(rows != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &tableaurows, nsubspacevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tableaurowinds, nsubspacevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ntableaurowinds, nsubspacevars) );
   for( j = nsubspacevars - 1; j >= 0; --j )
   {
      /*lint --e{866}*/
      SCIP_CALL( SCIPallocBufferArray(scip, &tableaurows[j], nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tableaurowinds[j], nrows) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &rownorm, nrows) );
   BMSclearMemoryArray(rownorm, nrows);

   /* clear ray */
   BMSclearMemoryArray(raydirection, nsubspacevars);

   /* get the relevant columns of the simplex tableau */
   for( j = nsubspacevars - 1; j >= 0; --j )
   {
      assert(SCIPcolGetLPPos(SCIPvarGetCol(subspacevars[j])) >= 0);
      SCIP_CALL( SCIPgetLPBInvACol(scip, SCIPcolGetLPPos(SCIPvarGetCol(subspacevars[j])), tableaurows[j], tableaurowinds[j], &ntableaurowinds[j]) );

      if( ntableaurowinds[j] == -1 )
      {
         assert(sparse == 0 || sparse == -1);
         sparse = 0;

         for( i = nrows - 1; i >= 0; --i )
            rownorm[i] += tableaurows[j][i] * tableaurows[j][i];
      }
      else
      {
         assert(sparse == 1 || sparse == -1);
         sparse = 1;

         /* allocate temporary memory */
         if( usedrowids == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &rowids, nrows) );
            SCIP_CALL( SCIPallocBufferArray(scip, &usedrowids, nrows) );
            BMSclearMemoryArray(usedrowids, nrows);
         }

         for( i = ntableaurowinds[j] - 1; i >= 0; --i )
         {
            tableaurowind = tableaurowinds[j][i];
            rownorm[tableaurowind] += tableaurows[j][tableaurowind] * tableaurows[j][tableaurowind];
            if( !usedrowids[tableaurowind] )
            {
               usedrowids[tableaurowind] = TRUE;
               rowids[nrowids] = tableaurowind; /*lint !e644*/
               ++nrowids;
               assert(nrowids <= nrows);
            }
         }
      }
   }

   /* compute ray direction (dense) */
   if( sparse == 0 )
   {
      /* take average over all rows of the tableau */
      for( i = nrows - 1; i >= 0; --i )
      {
         if( SCIPisFeasZero(scip, rownorm[i]) )
            continue;
         else
            rownorm[i] = SQRT(rownorm[i]);

         if( weighted )
         {
            rowweight = SCIProwGetDualsol(rows[i]);
            if( SCIPisFeasZero(scip, rowweight) )
               continue;
         }
         else
            rowweight = 1.0;

         for( j = nsubspacevars - 1; j >= 0; --j )
         {
            raydirection[j] += tableaurows[j][i] / (rownorm[i] * rowweight);
            assert( ! SCIPisInfinity(scip, REALABS(raydirection[j])) );
         }
      }
   }
   /* compute ray direction (sparse) */
   else
   {
      SCIP_Real* rowweights;
      int r;
      int k;

      assert(usedrowids != NULL);
      assert(rowids != NULL);
      assert(sparse == 1);
      assert(0 <= nrowids && nrowids <= nrows);

      SCIP_CALL( SCIPallocBufferArray(scip, &rowweights, nrows) );

      /* compute norms of important rows and rowweights */
      for( i = nrowids - 1; i >= 0; --i )
      {
         r = rowids[i];
         assert(0 <= r && r < nrows);
         assert(usedrowids[r]);

         if( SCIPisFeasZero(scip, rownorm[r]) )
         {
            usedrowids[r] = FALSE;
            --nrowids;
            continue;
         }
         else
            rownorm[r] = SQRT(rownorm[r]);

         if( weighted )
         {
            rowweights[r] = SCIProwGetDualsol(rows[r]);
            if( SCIPisFeasZero(scip, rowweights[r]) )
            {
               usedrowids[r] = FALSE;
               --nrowids;
               continue;
            }
         }
         else
            rowweights[r] = 1.0;
      }

      if( nrowids > 0 )
      {
         /* take average over all rows of the tableau */
         for( j = nsubspacevars - 1; j >= 0; --j )
         {
            for( k = ntableaurowinds[j] - 1; k >= 0; --k )
            {
               tableaurowind = tableaurowinds[j][k];

               if( usedrowids[tableaurowind] )
               {
                  raydirection[j] += tableaurows[j][tableaurowind] / (rownorm[tableaurowind] * rowweights[tableaurowind]);
                  assert( ! SCIPisInfinity(scip, REALABS(raydirection[j])) );
               }
            }
         }
      }

      SCIPfreeBufferArray(scip, &rowweights);
      SCIPfreeBufferArray(scip, &usedrowids);
      SCIPfreeBufferArray(scip, &rowids);
   }
   assert(usedrowids == NULL);

   /* free memory */
   SCIPfreeBufferArray(scip, &rownorm);
   for( j = 0; j < nsubspacevars; ++j )
   {
      SCIPfreeBufferArray(scip, &tableaurowinds[j]);
      SCIPfreeBufferArray(scip, &tableaurows[j]);
   }
   SCIPfreeBufferArray(scip, &ntableaurowinds);
   SCIPfreeBufferArray(scip, &tableaurowinds);
   SCIPfreeBufferArray(scip, &tableaurows);

   return SCIP_OKAY;
}


/** generates the direction of the shooting ray as the average of the normalized non-basic vars and rows */
static
SCIP_RETCODE generateAverageNBRay(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Real*            raydirection,       /**< shooting ray                          */
   int*                  fracspace,          /**< index set of fractional variables     */
   SCIP_VAR**            subspacevars,       /**< pointer to fractional space variables */
   int                   nsubspacevars       /**< dimension of fractional space         */
   )
{
   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nrows;
   int ncols;
   int i;

   assert(scip != NULL);
   assert(raydirection != NULL);
   assert(fracspace != NULL);
   assert(subspacevars != NULL);

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* add up non-basic variables */
   for( i = nsubspacevars - 1; i >= 0; --i )
   {
      SCIP_Real solval;

      solval = SCIPvarGetLPSol(subspacevars[i]);

      if( SCIPisFeasEQ(scip, solval, SCIPvarGetLbLocal(subspacevars[i])) )
         raydirection[i] = +1.0;
      else if( SCIPisFeasEQ(scip, solval, SCIPvarGetUbLocal(subspacevars[i])) )
         raydirection[i] = -1.0;
      else
         raydirection[i] = 0.0;
   }

   /* add up non-basic rows */
   for( i = nrows - 1; i >= 0; --i )
   {
      SCIP_Real dualsol;
      SCIP_Real factor;
      SCIP_Real* coeffs;
      SCIP_Real rownorm;
      int j;
      int nnonz;

      dualsol = SCIProwGetDualsol(rows[i]);
      if( SCIPisFeasPositive(scip, dualsol) )
         factor = 1.0;
      else if( SCIPisFeasNegative(scip, dualsol) )
         factor = -1.0;
      else
         continue;

      /* get the row's data */
      coeffs = SCIProwGetVals(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      nnonz = SCIProwGetNNonz(rows[i]);

      rownorm = 0.0;
      for( j = nnonz - 1; j >= 0; --j )
      {
         SCIP_VAR* var;
         var = SCIPcolGetVar(cols[j]);
         if( fracspace[SCIPvarGetProbindex(var)] >= 0 )
            rownorm += coeffs[j] * coeffs[j];
      }

      if( SCIPisFeasZero(scip,rownorm) )
         continue;
      else
      {
         assert(rownorm > 0);
         rownorm = SQRT(rownorm);
      }

      for( j = nnonz - 1; j >= 0; --j )
      {
         SCIP_VAR* var;
         int f;

         var = SCIPcolGetVar(cols[j]);
         f = fracspace[SCIPvarGetProbindex(var)];

         if( f >= 0 )
         {
            raydirection[f] += factor * coeffs[j] / rownorm;
            assert( ! SCIPisInfinity(scip, REALABS(raydirection[f])) );
         }
      }
   }
   return SCIP_OKAY;
}

/** generates the starting point for the shooting ray in original coordinates */
static
SCIP_RETCODE generateStartingPoint(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Real*            rayorigin,          /**< origin of the shooting ray            */
   SCIP_VAR**            subspacevars,       /**< pointer to fractional space variables */
   int                   nsubspacevars       /**< dimension of fractional space         */
   )
{
   int v;

   assert(scip != NULL);
   assert(rayorigin != NULL);
   assert(subspacevars != NULL);

   for( v = nsubspacevars - 1; v >= 0; --v )
      rayorigin[v] = SCIPvarGetLPSol(subspacevars[v]);

   return SCIP_OKAY;
}

/** translates the inner point of the LP to an inner point rayorigin of the unit hyper octahedron and
 *  transforms raydirection and rayorigin by reflections stored in sign
 */
static
void flipCoords(
   SCIP_Real*            rayorigin,          /**< origin of the shooting ray            */
   SCIP_Real*            raydirection,       /**< direction of the shooting ray         */
   SCIP_Bool*            sign,               /**< marker for flipped coordinates        */
   int                   nsubspacevars       /**< dimension of fractional space         */
   )
{
   int v;

   assert(rayorigin != NULL);
   assert(raydirection != NULL);
   assert(sign != NULL);

   for( v = nsubspacevars - 1; v >= 0; --v )
   {
      /* if raydirection[v] is negative, flip its sign */
      if( raydirection[v] < 0 )
      {
         sign[v] = FALSE;
         raydirection[v] *= -1.0;
         rayorigin[v] *= -1.0; /* flip starting point in the same way like raydirection */
      }
      else
         sign[v] = TRUE;
   }
}

/** generates all facets, from which facet i could be obtained by a decreasing + to - flip
 *  or a nonincreasing - to + flip and tests whether they are among the fmax nearest ones
 */
static
void generateNeighborFacets(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Bool**           facets,             /**< facets got so far                     */
   SCIP_Real*            lambda,             /**< distances of the facets               */
   SCIP_Real*            rayorigin,          /**< origin of the shooting ray            */
   SCIP_Real*            raydirection,       /**< direction of the shooting ray         */
   SCIP_Real*            negquotient,        /**< array by which coordinates are sorted */
   int                   nsubspacevars,      /**< dimension of fractional space         */
   int                   f_max,              /**< maximal number of facets to create    */
   int                   i,                  /**< current facet                         */
   int*                  nfacets             /**< number of facets                      */
   )
{
   SCIP_Real p;
   SCIP_Real q;
   SCIP_Real lam;
   int minplus;
   int j;

   assert(scip != NULL);
   assert(facets != NULL);
   assert(facets[i] != NULL);
   assert(lambda != NULL);
   assert(rayorigin != NULL);
   assert(raydirection != NULL);
   assert(negquotient != NULL);
   assert(nfacets != NULL);
   assert(0 <= i && i < f_max);

   /* determine the p and q values of the next facet to fix as a closest one */
   p = 0.5 * nsubspacevars;
   q = 0.0;
   for( j = nsubspacevars - 1; j >= 0; --j )
   {
      if( facets[i][j] )
      {
         p -= rayorigin[j];
         q += raydirection[j];
      }
      else
      {
         p += rayorigin[j];
         q -= raydirection[j];
      }
   }

   /* get the first + entry of the facet */
   minplus = -1;
   for( j = 0; j < nsubspacevars; ++j )
   {
      if( facets[i][j] )
      {
         minplus = j;
         break;
      }
   }

   /* facet (- - ... -) cannot be hit, because raydirection >= 0 */
   assert(minplus >= 0);
   assert(q != 0.0);
   assert(SCIPisFeasEQ(scip, lambda[i], p/q));
   assert(lambda[i] >= 0.0);

   /* reverse search for facets from which the actual facet can be got by a single, decreasing + to - flip */
   /* a facet will be inserted into the queue, iff it is one of the fmax closest ones already found */
   for( j = 0; j < nsubspacevars && !facets[i][j] && SCIPisFeasGT(scip, negquotient[j], lambda[i]); ++j )
   {
      if( SCIPisFeasPositive(scip, q + 2*raydirection[j]) )
      {
         lam = (p - 2*rayorigin[j]) / (q + 2*raydirection[j]);
         tryToInsert(scip, facets, lambda, i, j, f_max, nsubspacevars, lam, nfacets);
      }
   }

   /* reverse search for facets from which the actual facet can be got by a single, nonincreasing - to + flip */
   /* a facet will be inserted into the queue, iff it is one of the fmax closest ones already found */
   for( j = nsubspacevars - 1; j >= 0 && facets[i][j] && SCIPisFeasLE(scip, negquotient[j], lambda[i]); --j )
   {
      if( SCIPisFeasPositive(scip, q - 2*raydirection[j]) )
      {
         lam = (p + 2*rayorigin[j]) / (q - 2*raydirection[j]);
         if( negquotient[minplus] <= lam )
            tryToInsert(scip, facets, lambda, i, j, f_max, nsubspacevars, lam, nfacets);
      }
   }
#ifndef NDEBUG
   for( j = 1; j < f_max; j++)
      assert(SCIPisFeasGE(scip, lambda[j], lambda[j-1]));
#endif
}

/** tests, whether an array is completely zero */
static
SCIP_Bool isZero(
   SCIP*                 scip,               /**< SCIP data structure                   */
   SCIP_Real*            raydirection,       /**< array to be checked                   */
   int                   nsubspacevars       /**< size of array                         */
   )
{
   int v;
   SCIP_Bool iszero;

   assert(scip != NULL);
   assert(raydirection != NULL);
   iszero = TRUE; 
   for( v = nsubspacevars - 1; v >= 0; --v )
   {
      assert(!SCIPisInfinity(scip, raydirection[v]));

      if( !SCIPisFeasZero(scip, raydirection[v]/100) )
         iszero = FALSE;
      else
         raydirection[v] = 0.0;
   }
   return iszero;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyOctane)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurOctane(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeOctane)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitOctane)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->lastrule = 0;
   heurdata->nsuccess = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */

static
SCIP_DECL_HEUREXIT(heurExitOctane)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecOctane)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   SCIP_SOL** first_sols;     /* stores the first ffirst sols in order to check for common violation of a row */

   SCIP_VAR** vars;           /* the variables of the problem */
   SCIP_VAR** fracvars;       /* variables, that are fractional in current LP solution */
   SCIP_VAR** subspacevars;   /* the variables on which the search is performed. Either coinciding with vars or with the
                               * space of all fractional variables of the current LP solution */

   SCIP_Real p;               /* n/2 - <delta,x> ( for some facet delta ) */
   SCIP_Real q;               /* <delta,a> */

   SCIP_Real* rayorigin;      /* origin of the ray, vector x in paper */
   SCIP_Real* raydirection;   /* direction of the ray, vector a in paper */
   SCIP_Real* negquotient;    /* negated quotient of rayorigin and raydirection, vector v in paper */
   SCIP_Real* lambda;         /* stores the distance of the facets (s.b.) to the origin of the ray */

   SCIP_Bool usefracspace;    /* determines whether the search concentrates on fractional variables and fixes integer ones */
   SCIP_Bool cons_viol;       /* used for checking whether a linear constraint is violated by one of the possible solutions */
   SCIP_Bool success;
   SCIP_Bool* sign;           /* signature of the direction of the ray */
   SCIP_Bool** facets;        /* list of extended facets */

   int nvars;            /* number of variables  */
   int nbinvars;         /* number of 0-1-variables */
   int nfracvars;        /* number of fractional variables in current LP solution */
   int nsubspacevars;    /* dimension of the subspace on which the search is performed */
   int nfacets;          /* number of facets hidden by the ray that where already found */
   int i;                /* counter */
   int j;                /* counter */
   int f_max;            /* {0,1}-points to be checked */
   int f_first;          /* {0,1}-points to be generated at first in order to check whether a restart is necessary */
   int r;                /* counter */
   int firstrule;

   int* perm;            /* stores the way in which the coordinates were permuted */
   int* fracspace;       /* maps the variables of the subspace to the original variables */

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );

   /* OCTANE is for use in 0-1 programs only */
   if( nvars != nbinvars )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* don't call heuristic, if it was not successful enough in the past */
   /*lint --e{647}*/
   if( SCIPgetNNodes(scip) % (SCIPheurGetNCalls(heur) / (100 * SCIPheurGetNBestSolsFound(heur) + 10*heurdata->nsuccess + 1) + 1) != 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetLPBranchCands(scip, &fracvars, NULL, NULL, &nfracvars, NULL, NULL) );

   /* don't use integral starting points */
   if( nfracvars == 0 )
      return SCIP_OKAY;

   /* get working pointers from heurdata */
   sol = heurdata->sol;
   assert( sol != NULL );
   f_max = heurdata->f_max;
   f_first = heurdata->f_first;
   usefracspace = heurdata->usefracspace;

   SCIP_CALL( SCIPallocBufferArray(scip, &fracspace, nvars) );

   /* determine the space one which OCTANE should work either as the whole space or as the space of fractional variables */
   if( usefracspace )
   {
      nsubspacevars = nfracvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &subspacevars, nsubspacevars) );
      BMScopyMemoryArray(subspacevars, fracvars, nsubspacevars);
      for( i = nvars - 1; i >= 0; --i )
         fracspace[i] = -1;
      for( i = nsubspacevars - 1; i >= 0; --i )
         fracspace[SCIPvarGetProbindex(subspacevars[i])] = i;
   }
   else
   {
      int currentindex;

      nsubspacevars = nvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &subspacevars, nsubspacevars) );

      /* only copy the variables which are in the current LP */
      currentindex = 0;
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPcolGetLPPos(SCIPvarGetCol(vars[i])) >= 0 )
         {
            subspacevars[currentindex] = vars[i];
            fracspace[i] = currentindex;
            ++currentindex;

         }
         else
         {
            fracspace[i] = -1;
            --nsubspacevars;
         }
      }
   }

   /* nothing to do for empty search space */
   if( nsubspacevars == 0 )
      return SCIP_OKAY;

   assert(0 < nsubspacevars && nsubspacevars <= nvars);

   for( i = 0; i < nsubspacevars; i++)
      assert(fracspace[SCIPvarGetProbindex(subspacevars[i])] == i);

   /* at most 2^(n-1) facets can be hit */
   if( nsubspacevars < 30 )
   {
      /*lint --e{701}*/
      assert(f_max > 0);
      f_max = MIN(f_max, 1 << (nsubspacevars - 1) );
   }

   f_first = MIN(f_first, f_max);

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &rayorigin, nsubspacevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &raydirection, nsubspacevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negquotient, nsubspacevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sign, nsubspacevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nsubspacevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lambda, f_max + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &facets, f_max + 1) );


   for( i = f_max; i >= 0; --i )
   {
      /*lint --e{866}*/
      SCIP_CALL( SCIPallocBufferArray(scip, &facets[i], nsubspacevars) );
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &first_sols, f_first) );

   *result = SCIP_DIDNOTFIND;

   /* starting OCTANE */
   SCIPdebugMsg(scip, "run Octane heuristic on %s variables, which are %d vars, generate at most %d facets, using rule number %d\n",
      usefracspace ? "fractional" : "all", nsubspacevars, f_max, (heurdata->lastrule+1)%5);

   /* generate starting point in original coordinates */
   SCIP_CALL( generateStartingPoint(scip, rayorigin, subspacevars, nsubspacevars) );
   for( i = nsubspacevars - 1; i >= 0; --i )
      rayorigin[i] -= 0.5;

   firstrule = heurdata->lastrule;
   ++firstrule;
   for( r = firstrule; r <= firstrule + 5 && !SCIPisStopped(scip); r++ )
   {
      SCIP_ROW** rows;
      int nrows;

      /* generate shooting ray in original coordinates by certain rules */
      switch(r % 5)
      {
      case 1:
         if( !heurdata->useavgnbray )
            continue;

         SCIP_CALL( generateAverageNBRay(scip, raydirection, fracspace, subspacevars, nsubspacevars) );
         break;
      case 2:
         if( !heurdata->useobjray )
            continue;

         SCIP_CALL( generateObjectiveRay(scip, raydirection, subspacevars, nsubspacevars) );
         break;
      case 3:
         if( !heurdata->usediffray )
            continue;

         SCIP_CALL( generateDifferenceRay(scip, raydirection, subspacevars, nsubspacevars) );
         break;
      case 4:
         if( !heurdata->useavgwgtray || !SCIPisLPSolBasic(scip) )
            continue;

         SCIP_CALL( generateAverageRay(scip, raydirection, subspacevars, nsubspacevars, TRUE) );
         break;
      case 0:
         if( !heurdata->useavgray || !SCIPisLPSolBasic(scip) )
            continue;

         SCIP_CALL( generateAverageRay(scip, raydirection, subspacevars, nsubspacevars, FALSE) );
         break;
      default:
         SCIPerrorMessage("invalid ray rule identifier\n");
         SCIPABORT();
      }

      /* there must be a feasible direction for the shooting ray */
      if( isZero(scip, raydirection, nsubspacevars) )
         continue;

      /* transform coordinates such that raydirection >= 0 */
      flipCoords(rayorigin, raydirection, sign, nsubspacevars);

      for( i = f_max - 1; i >= 0; --i)
         lambda[i] = SCIPinfinity(scip);

      /* calculate negquotient, initialize perm, facets[0], p, and q */
      p = 0.5 * nsubspacevars;
      q = 0.0;
      for( i = nsubspacevars - 1; i >= 0; --i )
      {
         /* calculate negquotient, the ratio of rayorigin and raydirection, paying special attention to the case raydirection[i] == 0 */
         if( SCIPisFeasZero(scip, raydirection[i]) )
         {
            if( rayorigin[i] < 0 )
               negquotient[i] = SCIPinfinity(scip);
            else
               negquotient[i] = -SCIPinfinity(scip);
         }
         else
            negquotient[i] = - (rayorigin[i] / raydirection[i]);

         perm[i] = i;

         /* initialization of facets[0] to the all-one facet with p and q its characteristic values */
         facets[0][i] = TRUE;
         p -= rayorigin[i];
         q += raydirection[i];
      }

      assert(SCIPisPositive(scip, q));

      /* resort the coordinates in nonincreasing order of negquotient */
      SCIPsortDownRealRealRealBoolPtr(negquotient, raydirection, rayorigin, sign, (void**) subspacevars, nsubspacevars);

#ifndef NDEBUG
      for( i = 0; i < nsubspacevars; i++ )
      {
         assert( raydirection[i] >= 0 );
         assert(!SCIPisInfinity(scip, REALABS(raydirection[i])));
      }
      for( i = 1; i < nsubspacevars; i++ )
         assert( negquotient[i - 1] >= negquotient[i] );
#endif
      /* finished initialization */

      /* find the first facet of the octahedron hit by a ray shot from rayorigin into direction raydirection */
      for( i = 0; i < nsubspacevars && negquotient[i] * q > p; ++i )
      {
         facets[0][i] = FALSE;
         p += 2 * rayorigin[i];
         q -= 2 * raydirection[i];
         assert(SCIPisPositive(scip, p));
         assert(SCIPisPositive(scip, q));
      }

      /* avoid dividing by values close to 0.0 */
      if( !SCIPisFeasPositive(scip, q) )
         continue;

      /* assert necessary for flexelint */
      assert(q != 0.0);
      lambda[0] = p / q;

      nfacets = 1;

      /* find the first facets hit by the ray */
      for( i = 0; i < nfacets && i < f_first; ++i)
         generateNeighborFacets(scip, facets, lambda, rayorigin, raydirection, negquotient, nsubspacevars, f_max, i, &nfacets);

      /* construct the first ffirst possible solutions */
      for( i = 0; i < nfacets && i < f_first; ++i )
      {
         SCIP_CALL( SCIPcreateSol(scip, &first_sols[i], heur) );
         SCIP_CALL( getSolFromFacet(scip, facets[i], first_sols[i], sign, subspacevars, nsubspacevars) );
         assert( first_sols[i] != NULL );
      }

      /* try, whether there is a row violated by all of the first ffirst solutions */
      cons_viol = FALSE;
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
      for( i = nrows - 1; i >= 0; --i )
      {
         if( !SCIProwIsLocal(rows[i]) )
         {
            SCIP_COL** cols;
            SCIP_Real constant;
            SCIP_Real lhs;
            SCIP_Real rhs;
            SCIP_Real rowval;
            SCIP_Real* coeffs;
            int nnonzerovars;
            int k;

            /* get the row's data */
            constant = SCIProwGetConstant(rows[i]);
            lhs = SCIProwGetLhs(rows[i]);
            rhs = SCIProwGetRhs(rows[i]);
            coeffs = SCIProwGetVals(rows[i]);
            nnonzerovars = SCIProwGetNNonz(rows[i]);
            cols = SCIProwGetCols(rows[i]);
            rowval = constant;

            for( j = nnonzerovars - 1; j >= 0; --j )
               rowval += coeffs[j] * SCIPgetSolVal(scip, first_sols[0], SCIPcolGetVar(cols[j]));

            /* if the row's lhs is violated by the first sol, test, whether it is violated by the next ones, too */
            if( lhs > rowval )
            {
               cons_viol = TRUE;
               for( k = MIN(f_first, nfacets) - 1; k > 0; --k )
               {
                  rowval = constant;
                  for( j = nnonzerovars - 1; j >= 0; --j )
                     rowval += coeffs[j] * SCIPgetSolVal(scip, first_sols[k], SCIPcolGetVar(cols[j]));
                  if( lhs <= rowval )
                  {
                     cons_viol = FALSE;
                     break;
                  }
               }
            }
            /* dito for the right hand side */
            else if( rhs < rowval )
            {
               cons_viol = TRUE;
               for( k = MIN(f_first, nfacets) - 1; k > 0; --k )
               {
                  rowval = constant;
                  for( j = nnonzerovars - 1; j >= 0; --j )
                     rowval += coeffs[j] * SCIPgetSolVal(scip, first_sols[k], SCIPcolGetVar(cols[j]));
                  if( rhs >= rowval )
                  {
                     cons_viol = FALSE;
                     break;
                  }
               }
            }
            /* break as soon as one row is violated by all of the ffirst solutions */
            if( cons_viol )
               break;
         }
      }


      if( !cons_viol )
      {
         /* if there was no row violated by all solutions, try whether one or more of them are feasible */
         for( i = MIN(f_first, nfacets) - 1; i >= 0; --i )
         {
            assert(first_sols[i] != NULL);
            SCIP_CALL( SCIPtrySol(scip, first_sols[i], FALSE, FALSE, TRUE, FALSE, TRUE, &success) );
            if( success )
               *result = SCIP_FOUNDSOL;
         }
         /* search for further facets and construct and try solutions out of facets fixed as closest ones */
         for( i = f_first; i < f_max; ++i)
         {
            if( i >= nfacets )
               break;
            generateNeighborFacets(scip, facets, lambda, rayorigin, raydirection, negquotient, nsubspacevars, f_max, i, &nfacets);
            SCIP_CALL( getSolFromFacet(scip, facets[i], sol, sign, subspacevars, nsubspacevars) );
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, FALSE, TRUE, &success) );
            if( success )
               *result = SCIP_FOUNDSOL;
         }
      }

      /* finished OCTANE */
      for( i = MIN(f_first, nfacets) - 1; i >= 0; --i )
      {
         SCIP_CALL( SCIPfreeSol(scip, &first_sols[i]) );
      }
   }
   heurdata->lastrule = r;

   if( *result == SCIP_FOUNDSOL )
      ++(heurdata->nsuccess);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &first_sols);
   for( i = 0; i <= f_max; ++i )
      SCIPfreeBufferArray(scip, &facets[i]);
   SCIPfreeBufferArray(scip, &facets);
   SCIPfreeBufferArray(scip, &lambda);
   SCIPfreeBufferArray(scip, &perm);
   SCIPfreeBufferArray(scip, &sign);
   SCIPfreeBufferArray(scip, &negquotient);
   SCIPfreeBufferArray(scip, &raydirection);
   SCIPfreeBufferArray(scip, &rayorigin);
   SCIPfreeBufferArray(scip, &subspacevars);
   SCIPfreeBufferArray(scip, &fracspace);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the octane primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurOctane(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Octane primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecOctane, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyOctane) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeOctane) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitOctane) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitOctane) );

   /* add octane primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/octane/fmax",
         "number of 0-1-points to be tested as possible solutions by OCTANE",
         &heurdata->f_max, TRUE, DEFAULT_FMAX, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/octane/ffirst",
         "number of 0-1-points to be tested at first whether they violate a common row",
         &heurdata->f_first, TRUE, DEFAULT_FFIRST, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/octane/usefracspace",
         "execute OCTANE only in the space of fractional variables (TRUE) or in the full space?",
         &heurdata->usefracspace, TRUE, DEFAULT_USEFRACSPACE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/octane/useobjray",
         "should the inner normal of the objective be used as one ray direction?",
         &heurdata->useobjray, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/octane/useavgray",
         "should the average of the basic cone be used as one ray direction?",
         &heurdata->useavgray, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/octane/usediffray",
         "should the difference between the root solution and the current LP solution be used as one ray direction?",
         &heurdata->usediffray, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/octane/useavgwgtray",
         "should the weighted average of the basic cone be used as one ray direction?",
         &heurdata->useavgwgtray, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/octane/useavgnbray",
         "should the weighted average of the nonbasic cone be used as one ray direction?",
         &heurdata->useavgnbray, TRUE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}
