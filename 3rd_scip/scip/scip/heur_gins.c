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

/**@file   heur_gins.c
 * @brief  LNS heuristic that tries to delimit the search region to a neighborhood in the constraint graph
 * @author Gregor Hendel
 *
 * Graph Induced Neighborhood Search (GINS) is a Large Neighborhood Search Heuristic that attempts to improve
 * an incumbent solution by fixing a suitable percentage of integer variables to the incumbent and
 * solving the resulting, smaller and presumably easier sub-MIP.
 *
 * Its search neighborhoods are based on distances in a bipartite graph \f$G\f$ with the variables and constraints as nodes
 * and an edge between a variable and a constraint, if the variable is part of the constraint.
 * Given an integer \f$k\f$, the \f$k\f$-neighborhood of a variable \f$v\f$ in \f$G\f$ is the set of variables, whose nodes
 * are connected to \f$v\f$ by a path not longer than \f$2 \cdot k\f$. Intuitively, a judiciously chosen neighborhood size
 * allows to consider a local portion of the overall problem.
 *
 * An initial variable selection is made by randomly sampling different neighborhoods across the whole main problem.
 * The neighborhood that offers the largest potential for improvement is selected to become the local search neighborhood,
 * while all variables outside the neighborhood are fixed to their incumbent solution values.
 *
 * GINS also supports a rolling horizon approach, during which several local neighborhoods are considered
 * with increasing distance to the variable selected for the initial sub-problem. The rolling horizon approach ends
 * if no improvement could be found or a sufficient part of the problem component variables has been part of
 * at least one neighborhood.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/heur_gins.h"
#include "scip/pub_misc.h"

#define HEUR_NAME             "gins"
#define HEUR_DESC             "gins works on k-neighborhood in a variable-constraint graph"
#define HEUR_DISPCHAR         'K'
#define HEUR_PRIORITY         -1103000
#define HEUR_FREQ             20
#define HEUR_FREQOFS          8
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODESOFS      500           /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_MAXNODES      5000          /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINIMPROVE    0.01          /**< factor by which Gins should at least improve the incumbent */
#define DEFAULT_MINNODES      50            /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_MINFIXINGRATE 0.66          /**< minimum percentage of integer variables that have to be fixed */
#define DEFAULT_NODESQUOT     0.15          /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_NWAITINGNODES 100           /**< number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_USELPROWS     FALSE         /**< should subproblem be created out of the rows in the LP rows,
                                             **< otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE          /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the
                                             **< cutpool of the original scip be copied to constraints of the subscip */
#define DEFAULT_BESTSOLLIMIT    3            /**< limit on number of improving incumbent solutions in sub-CIP */
#define DEFAULT_FIXCONTVARS FALSE           /**< should continuous variables outside the neighborhoods get fixed? */
#define DEFAULT_POTENTIAL      'r'          /**< the reference point to compute the neighborhood potential: (r)oot or (p)seudo solution */
#define DEFAULT_MAXDISTANCE     3           /**< maximum distance to selected variable to enter the subproblem, or -1 to
                                             *   select the distance that best approximates the minimum fixing rate from below */
#define DEFAULT_RANDSEED       71
#define DEFAULT_RELAXDENSECONSS FALSE       /**< should dense constraints (at least as dense as 1 - minfixingrate) be
                                             *   ignored by connectivity graph? */
#define DEFAULT_USEROLLINGHORIZON TRUE      /**< should the heuristic solve a sequence of sub-MIP's around the first selected variable */
#define DEFAULT_ROLLHORIZONLIMFAC  0.4      /**< limiting percentage for variables already used in sub-SCIPs to terminate rolling
                                             *   horizon approach */
#ifdef SCIP_STATISTIC
#define NHISTOGRAMBINS         10           /* number of bins for histograms */
#endif
/*
 * Data structures
 */

/** rolling horizon data structure to control multiple LNS heuristic runs away from an original source variable
 *
 */
struct RollingHorizon
{
   SCIP_VGRAPH*          variablegraph;      /**< variable graph data structure for breadth-first-search neighborhoods */
   int*                  distances;          /**< distances of the heuristic rolling horizon from the original source
                                              *   variable indexed by probindex */
   SCIP_Bool*            used;               /**< array that represents for every variable whether it has been used
                                              *   in a neighborhood indexed by probindex */
   int                   lastmaxdistance;    /**< the last distance k for a neighborhood, will be decreased
                                              *   during the rolling horizon if the selected neighborhood is too large */
   int                   lastdistance;       /**< last distance from originally selected variable in iteration zero */
   int                   distancessize;      /**< size of the distances and used arrays */
   int                   niterations;        /**< counter for the number of rolling horizon iterations */
   int                   nused;              /**< counts the number variables that have been part of any neighborhood
                                              *   during the rolling horizon approach */
   int                   nnonreachable;      /**< counter for the number of nonreachable variables (distance -1) from
                                              *   the initially selected variable */
};
typedef struct RollingHorizon ROLLINGHORIZON;

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   int                   minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed */
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minimprove;         /**< factor by which Gins should at least improve the incumbent */
   SCIP_Longint          usednodes;          /**< nodes already used by Gins in earlier calls */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             rollhorizonlimfac;  /**< limiting percentage for variables already used in sub-SCIPs to terminate
                                              *   rolling horizon approach */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator                              */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows? */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem? */
   SCIP_Bool             fixcontvars;        /**< should continuous variables outside the neighborhoods get fixed? */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP */
   int                   maxdistance;        /**< maximum distance to selected variable to enter the subproblem, or -1 to
                                              *   select the distance that best approximates the minimum fixing rate from below */
   int                   sumneighborhoodvars;/**< neighborhood variables sum over all seen neighborhoods */
   int                   sumdiscneighborhoodvars; /**< neighborhood discrete variables sum over all seen neighboorhoods */
   int                   nneighborhoods;     /**< number of calculated neighborhoods */
   int                   nsubmips;           /**< counter for the number of sub-MIP's that can be higher than the number of
                                              *   calls of this heuristic */
   SCIP_Bool             relaxdenseconss;    /**< should dense constraints (at least as dense as 1 - minfixingrate) be
                                              *   ignored by connectivity graph? */
   SCIP_Bool             userollinghorizon;  /**< should the heuristic solve a sequence of sub-MIP's around the first
                                              *   selected variable */
   char                  potential;          /**< the reference point to compute the neighborhood potential: (r)oot or
                                              *   (p)seudo solution */
   int                   maxseendistance;    /**< maximum of all distances between two variables */
#ifdef SCIP_STATISTIC
   int                   consvarshist[NHISTOGRAMBINS]; /**< histogram that summarizes the densities of the constraints */
   int                   consdiscvarshist[NHISTOGRAMBINS]; /**< histogram that summarizes the discrete variable densities of the constraints */
   int                   conscontvarshist[NHISTOGRAMBINS]; /**< histogram that summarizes the continuous variable densities of the constraints */
#endif
   int                   nrelaxedconstraints; /**< number of constraints that were relaxed */
   int                   nfailures;           /**< counter for the number of unsuccessful runs of this heuristic */
   SCIP_Longint          nextnodenumber;      /**< the next node number at which the heuristic should be called again */
};

/** represents limits for the sub-SCIP solving process */
struct SolveLimits
{
   SCIP_Longint          nodelimit;          /**< maximum number of solving nodes for the sub-SCIP */
   SCIP_Real             memorylimit;        /**< memory limit for the sub-SCIP */
   SCIP_Real             timelimit;          /**< time limit for the sub-SCIP */
};
typedef struct SolveLimits SOLVELIMITS;

/*
 * Local methods
 */

#ifdef SCIP_STATISTIC
/** resets a histogram */
static
void resetHistogram(
   int*                  histogram           /**< the histogram */
   )
{
   BMSclearMemoryArray(histogram, NHISTOGRAMBINS);
}

/** adds a ratio to the histogram at the right position */
static
void addHistogramEntry(
   int*                  histogram,          /**< the histogram */
   int                   value,              /**< the value */
   int                   basevalue           /**< base value */
   )
{
   SCIP_Real ratio;
   int index;
   assert(value <= basevalue);
   ratio = value/ MAX(1.0, (SCIP_Real)basevalue);

   index = (int)(ratio * NHISTOGRAMBINS);
   ++histogram[index];
}

#endif

/** create a rolling horizon data structure */
static
SCIP_RETCODE rollingHorizonCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   ROLLINGHORIZON**      rollinghorizon      /**< pointer to rolling horizon data structure */
   )
{
   int size;
   assert(scip != NULL);
   assert(rollinghorizon != NULL);

   size = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   SCIP_CALL( SCIPallocBlockMemory(scip, rollinghorizon) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*rollinghorizon)->distances, size) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(*rollinghorizon)->used, size) );
   (*rollinghorizon)->distancessize = size;
   (*rollinghorizon)->variablegraph = NULL;
   (*rollinghorizon)->lastdistance = -1;
   (*rollinghorizon)->niterations = 0;
   (*rollinghorizon)->nused = 0;

   return SCIP_OKAY;
}

/** free a rolling horizon data structure */
static
SCIP_RETCODE rollingHorizonFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ROLLINGHORIZON**      rollinghorizon      /**< pointer to rolling horizon data structure */
   )
{

   assert(scip != NULL);
   assert(rollinghorizon != NULL);
   assert(*rollinghorizon != NULL);

   if( (*rollinghorizon)->variablegraph != NULL )
   {
      SCIPvariableGraphFree(scip, &(*rollinghorizon)->variablegraph);
   }

   SCIPfreeBlockMemoryArray(scip, &(*rollinghorizon)->distances, (*rollinghorizon)->distancessize);
   SCIPfreeBlockMemoryArray(scip, &(*rollinghorizon)->used, (*rollinghorizon)->distancessize);
   SCIPfreeBlockMemory(scip, rollinghorizon);
   return SCIP_OKAY;
}

/** is there potential to run another iteration of the rolling horizon approach? */
static
SCIP_Bool rollingHorizonRunAgain(
   SCIP*                 scip,               /**< SCIP data structure */
   ROLLINGHORIZON*       rollinghorizon,     /**< rolling horizon data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   int maxnused = (int)(heurdata->rollhorizonlimfac * (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip)
         - rollinghorizon->nnonreachable));

   /* run again if a certain percentage of the reachable variables (in the same connected component)
    * was not used in a previous neighborhood
    */
   return (rollinghorizon->nused < maxnused);
}

/** store the distances from the selected variable permanently for the rolling horizon approach */
static
void rollingHorizonStoreDistances(
   SCIP*                 scip,               /**< SCIP data structure */
   ROLLINGHORIZON*       rollinghorizon,     /**< rolling horizon data structure */
   int*                  distances           /**< breadth-first distances indexed by Probindex of variables */
   )
{
   int i;
   BMScopyMemoryArray(rollinghorizon->distances, distances, rollinghorizon->distancessize);
   rollinghorizon->lastdistance = 0;
   rollinghorizon->nnonreachable = 0;

   /* collect number of nonreachable variables */
   for( i = 0; i < rollinghorizon->distancessize; ++i )
   {
      if( distances[i] == -1 )
         ++rollinghorizon->nnonreachable;
   }
}

/** get the potential of a subset of variables (distance to a reference point such as the pseudo-solution or root
 * LP solution)
 */
static
SCIP_Real getPotential(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< length of variable array */
   )
{
   SCIP_Real potential;
   int i;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(sol != NULL);

   if( nvars == 0 )
      return 0.0;

   potential = 0.0;

   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real objdelta;
      SCIP_VAR* var;
      SCIP_Real referencepoint;
      SCIP_Real varobj;

      var = vars[i];
      assert(var != NULL);
      varobj = SCIPvarGetObj(var);

      if( SCIPisZero(scip, varobj) )
         continue;

      /* determine the reference point for potential computation */
      switch( heurdata->potential )
      {
         /* use difference to pseudo solution using the bound in the objective direction */
         case 'p':
            referencepoint = SCIPvarGetObj(var) > 0.0 ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
            break;

         /* use root LP solution difference */
         case 'r':
            referencepoint = SCIPvarGetRootSol(var);
            break;
         default:
            SCIPerrorMessage("Unknown potential computation %c specified\n", heurdata->potential);
            referencepoint = 0.0;
            break;
      }

      if( SCIPisInfinity(scip, REALABS(referencepoint)) )
         continue;

      /* calculate the delta to the variables best bound */
      objdelta = (SCIPgetSolVal(scip, sol, var) - referencepoint) * SCIPvarGetObj(var);
      potential += objdelta;
   }

   return potential;
}

/** is the variable in the current neighborhood which is given by the breadth-first distances from a central variable? */
static
SCIP_Bool isVariableInNeighborhood(
   SCIP_VAR*             var,                /**< problem variable */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   int                   maxdistance         /**< maximum distance (inclusive) to be considered for neighborhoods */
   )
{
   assert(var != NULL);
   assert(distances != NULL);
   assert(maxdistance >= 0);
   assert(SCIPvarGetProbindex(var) >= 0);

   return (distances[SCIPvarGetProbindex(var)] != -1 && distances[SCIPvarGetProbindex(var)] <= maxdistance);
}

/** fixes variables in subproblem based on long breadth-first distances in variable graph */
static
SCIP_RETCODE fixNonNeighborhoodVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   ROLLINGHORIZON*       rollinghorizon,     /**< rolling horizon data structure to save relevant information, or NULL if not needed */
   SCIP_SOL*             sol,                /**< solution in main SCIP for fixing values */
   SCIP_VAR**            vars,               /**< variables in the main SCIP */
   SCIP_VAR**            fixedvars,          /**< buffer to store variables that should be fixed */
   SCIP_Real*            fixedvals,          /**< buffer to store fixing values for fixed variables */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   int                   maxdistance,        /**< maximum distance (inclusive) to be considered for neighborhoods */
   int*                  nfixings            /**< pointer to store number of fixed variables */
   )
{
   int i;
   int nbinvars;
   int nintvars;
   int nvars;
   int nvarstofix;

   SCIP_CALL( SCIPgetVarsData(scip, NULL, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nvarstofix = heurdata->fixcontvars ? nvars : nbinvars + nintvars;
   *nfixings = 0;

   /* change bounds of variables of the subproblem */
   for( i = 0; i < nvarstofix; i++ )
   {
      /* fix all variables that are too far away from this variable according to maxdistance */
      if( distances[i] == -1 || distances[i] > maxdistance )
      {
         SCIP_Real solval;
         SCIP_Real lb;
         SCIP_Real ub;

         solval = SCIPgetSolVal(scip, sol, vars[i]);
         lb = SCIPvarGetLbGlobal(vars[i]);
         ub = SCIPvarGetUbGlobal(vars[i]);
         assert(SCIPisLE(scip, lb, ub));

         /* due to dual reductions, it may happen that the solution value is not in the variable's domain anymore */
         if( SCIPisLT(scip, solval, lb) )
            solval = lb;
         else if( SCIPisGT(scip, solval, ub) )
            solval = ub;

         /* perform the bound change */
         if( !SCIPisInfinity(scip, solval) && !SCIPisInfinity(scip, -solval) )
         {
            fixedvars[*nfixings] = vars[i];
            fixedvals[*nfixings] = solval;
            ++(*nfixings);
         }
      }
      else if( rollinghorizon != NULL && i < nbinvars + nintvars && ! rollinghorizon->used[i] )
      {
         ++rollinghorizon->nused;
         rollinghorizon->used[i] = TRUE;
      }
   }

   if( rollinghorizon != NULL )
   {
      rollinghorizon->lastmaxdistance = maxdistance;
      rollinghorizon->niterations++;
   }

   return SCIP_OKAY;
}

/** determine the maximum allowed distance to stay within the restriction to fix at least minfixingrate many non
 *  neighborhood variables
 */
static
SCIP_RETCODE determineMaxDistance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   int*                  choosevardistance   /**< pointer to store the computed maximum distance */
   )
{
   int* distancescopy;
   int nrelevantdistances;
   int criticalidx;
   int zeropos;
   int nvars;
   int nbinvars;
   int nintvars;

   SCIP_CALL( SCIPgetVarsData(scip, NULL, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   nrelevantdistances = (heurdata->fixcontvars ?  nvars : (nbinvars + nintvars));

   /* copy the relevant distances of either the discrete or all problem variables and sort them */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &distancescopy, distances, nrelevantdistances) );
   SCIPsortInt(distancescopy, nrelevantdistances);

   /* distances can be infinite in the case of multiple connected components; therefore, search for the index of the
    * zero entry, which is the unique representative of the chosen variable in the sorted distances
    */
   zeropos = -1;

   /* TODO: use selection method instead */
   (void)SCIPsortedvecFindInt(distancescopy, 0, nrelevantdistances, &zeropos);
   assert(zeropos >= 0);

   /* determine the critical index to look for an appropriate neighborhood distance, starting from the zero position */
   criticalidx = zeropos + (int)((1.0 - heurdata->minfixingrate) * nrelevantdistances);
   criticalidx = MIN(criticalidx, nrelevantdistances - 1);

   /* determine the maximum breadth-first distance such that the neighborhood size stays small enough (below 1-minfixingrate)*/
   *choosevardistance = distancescopy[criticalidx];

   /* we set the distance to exactly the distance at the critical index. If the distance at critical index is not the
    * last one before the distance increases, we decrease the choosevardistance such that the entire neighborhood
    * fits into the minfixingrate restriction
    */
   if( criticalidx != nrelevantdistances - 1 && distancescopy[criticalidx + 1] == *choosevardistance )
      (*choosevardistance)--;

   /* update the maximum distance */
   heurdata->maxseendistance = MAX(heurdata->maxseendistance, distancescopy[nrelevantdistances - 1]);

   SCIPfreeBufferArray(scip, &distancescopy);

   return SCIP_OKAY;
}

/** gets the average size of a discrete neighborhood over all variables tested */
static
SCIP_Real heurdataAvgDiscreteNeighborhoodSize(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   return heurdata->sumdiscneighborhoodvars / (MAX(1.0, (SCIP_Real)heurdata->nneighborhoods));
}

/** select a good starting variable at the first iteration of a rolling horizon approach */
static
SCIP_RETCODE selectInitialVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_VGRAPH*          vargraph,           /**< variable graph data structure to work on */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   SCIP_VAR**            selvar,             /**< pointer to store the selected variable */
   int*                  selvarmaxdistance   /**< maximal distance k to consider for selected variable neighborhood */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;                          /* original scip variables */
   int nbinvars;
   int nintvars;
   int nvars;
   int nsearched;
   int searchlimit;
   int nintegralvarsleft;
   int nintegralvarsbound;
   int v;
   SCIP_VAR** varscopy;
   int maxdistance;
   SCIP_Real maxpotential;

   assert(vargraph != NULL);
   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(selvar != NULL);
   sol = SCIPgetBestSol(scip);
   assert(sol != NULL);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* copy SCIP variables */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &varscopy, vars, nbinvars + nintvars) );
   nsearched = 0;
   maxpotential = SCIP_REAL_MIN;

   /* determine upper bound on neighborhood size */
   nintegralvarsbound = (int)((1.0 - heurdata->minfixingrate) * (nbinvars + nintvars));

   /* maximum distance from selected variable for breadth-first search (if set to -1, we compute an exhaustive breadth-first
    * search and sort the variables by their distance)
    */
   maxdistance = (heurdata->maxdistance == - 1 ? INT_MAX : heurdata->maxdistance);

   nintegralvarsleft  = nbinvars + nintvars;
   *selvar = NULL;

   /* sort inactive variables to the end of the array */
   for( v = nintegralvarsleft - 1; v >= 0; --v )
   {
      if( ! SCIPvarIsActive(varscopy[v]) )
      {
         varscopy[v] = varscopy[nintegralvarsleft - 1];
         --nintegralvarsleft;
      }
   }

   /* adjust the search limit */
   searchlimit = heurdata->nneighborhoods < 10 ? 5 : (int)(nintegralvarsleft / heurdataAvgDiscreteNeighborhoodSize(heurdata));
   searchlimit = MIN(searchlimit, 10);

   /* multi variable potential: choose different disjoint neighborhoods, compare their potential */
   while( nsearched < searchlimit && nintegralvarsleft > 0 )
   {
      SCIP_VAR** neighborhood;
      SCIP_VAR* choosevar;
      int neighborhoodsize;
      int ndiscvarsneighborhood;
      int choosevardistance;

      ++nsearched;

      /* select a variable to start with randomly, but make sure it is active */
      do
      {
         int idx = SCIPrandomGetInt(heurdata->randnumgen, 0, nintegralvarsleft - 1);
         choosevar = varscopy[idx];
         assert(choosevar != NULL);
         /* sort inactive variables to the end */
         if( SCIPvarGetProbindex(choosevar) < 0 )
         {
            varscopy[idx] = varscopy[nintegralvarsleft - 1];
            --nintegralvarsleft;
         }
      }
      while( SCIPvarGetProbindex(choosevar) < 0 && nintegralvarsleft > 0);

      /* if there was no variable chosen, there are no active variables left */
      if( SCIPvarGetProbindex(choosevar) < 0 )
      {
         SCIPdebugMsg(scip, "No active variable left to perform breadth-first search\n");
         break;
      }

      assert(SCIPvarIsIntegral(choosevar));

      /* get neighborhood storage */
      SCIP_CALL( SCIPallocBufferArray(scip, &neighborhood, nvars) );

      /* determine breadth-first distances from the chosen variable */
      SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, vargraph, &choosevar, 1, distances, maxdistance, INT_MAX, INT_MAX) );

      /* use either automatic or user-defined max-distance for neighborhood in variable constraint graph */
      if( heurdata->maxdistance != -1 )
      {
         choosevardistance = heurdata->maxdistance;
      }
      else
      {
         SCIP_CALL( determineMaxDistance(scip, heurdata, distances, &choosevardistance) );
      }

      ndiscvarsneighborhood = 0;
      neighborhoodsize = 0;

      /* loop over variables and determine neighborhood */
      for( v = nvars - 1; v >= 0; --v )
      {
         SCIP_VAR* currvar;
         currvar = vars[v];

         /* put variable in the neighborhood */
         if( isVariableInNeighborhood(currvar, distances, choosevardistance) )
         {
            neighborhood[neighborhoodsize++] = currvar;

            /* increase discrete variables counter */
            if( SCIPvarIsIntegral(currvar) )
               ++ndiscvarsneighborhood;
         }
      }

      /* check if neighborhood contains too many integer variables in order to satisfy the minimum fixing rate */
      if( ndiscvarsneighborhood >= nintegralvarsbound || ndiscvarsneighborhood <= 1 )
      {
         SCIPdebugMsg(scip, "Too many or too few discrete variables in neighboorhood: %d (%d)\n",
            ndiscvarsneighborhood, nbinvars + nintvars);
      }
      else
      {
         /* compare the neighborhood potential to the best potential found so far */
         SCIP_Real potential = getPotential(scip, heurdata, sol, neighborhood, neighborhoodsize);

         /* big potential, take this variable */
         if( potential > maxpotential )
         {
            maxpotential = potential;
            *selvar = choosevar;
            *selvarmaxdistance = choosevardistance;
         }
      }

      /* sort neighborhood variables to the end so that this neighborhood is not considered in further samples */
      for( v = nintegralvarsleft - 1; v >= 0; --v )
      {
         SCIP_VAR* currvar;
         currvar = vars[v];

         if( isVariableInNeighborhood(currvar, distances, choosevardistance) )
         {
            varscopy[v] = varscopy[nintegralvarsleft - 1];
            --nintegralvarsleft;
         }
      }

      heurdata->sumdiscneighborhoodvars += ndiscvarsneighborhood;
      heurdata->sumneighborhoodvars += neighborhoodsize;
      ++heurdata->nneighborhoods;

      /* free current neighborhood */
      SCIPfreeBufferArray(scip, &neighborhood);
   }

   SCIPfreeBufferArray(scip, &varscopy);

   /* maybe no variable has a positive delta */
   if( !SCIPisPositive(scip, maxpotential) || *selvar == NULL )
   {
      SCIPdebugMsg(scip, "Stopping with maxpotential %15.9f and selected variable %s\n",
            maxpotential, *selvar != NULL ? SCIPvarGetName(*selvar) : "none");
      *selvar = NULL;
   }

   return SCIP_OKAY;
}

/** select the next variable using the rolling horizon */
static
SCIP_RETCODE selectNextVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   ROLLINGHORIZON*       rollinghorizon,     /**< rolling horizon data structure to save relevant information, or NULL if not needed */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   SCIP_VAR**            selvar,             /**< pointer to store the selected variable */
   int*                  selvarmaxdistance   /**< maximal distance k to consider for selected variable neighborhood */
   )
{
   SCIP_VAR** vars;                          /* original scip variables */
   int i;
   int nbinvars;
   int nintvars;
   int minunuseddistance;

   assert(scip != NULL);
   assert(rollinghorizon != NULL);
   assert(distances != NULL);
   assert(selvar != NULL);
   assert(selvarmaxdistance != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   /* loop over the variables that are left and pick the variable with
    * - the smallest, always nondecreasing distance
    * - that was not used before in a neighborhood
    */
   do
   {
      minunuseddistance = INT_MAX;
      *selvarmaxdistance = rollinghorizon->lastmaxdistance;
      *selvar = NULL;
      for( i = 0; i < nbinvars + nintvars && minunuseddistance > rollinghorizon->lastdistance; ++i )
      {
         if( rollinghorizon->distances[i] >= rollinghorizon->lastdistance
               && rollinghorizon->distances[i] < minunuseddistance && ! rollinghorizon->used[i] )
         {
            minunuseddistance = rollinghorizon->distances[i];
            *selvar = vars[i];
         }
      }

      /* if no variable could be selected, we can stop */
      if( *selvar == NULL )
      {
         *selvarmaxdistance = 0;
         return SCIP_OKAY;
      }

      /* determine the distances to other variables from this variable */
      SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, rollinghorizon->variablegraph, selvar, 1, distances, *selvarmaxdistance, INT_MAX, INT_MAX) );

      SCIP_CALL( determineMaxDistance(scip, heurdata, distances, selvarmaxdistance) );

      /* mark this variable as used in order to not find it again */
      if( *selvarmaxdistance == 0 )
      {
         rollinghorizon->used[SCIPvarGetProbindex(*selvar)] = TRUE;
         rollinghorizon->nused++;
         *selvar = NULL;
      }

   } while( rollingHorizonRunAgain(scip, rollinghorizon, heurdata) && (*selvar == NULL || *selvarmaxdistance == 0) );

   /* breadth-first search determines the distances of all variables
    * that are no more than maxdistance away from the start variable
    */
   assert(*selvarmaxdistance <= rollinghorizon->lastmaxdistance);
   *selvarmaxdistance = MIN(*selvarmaxdistance, rollinghorizon->lastmaxdistance);
   rollinghorizon->lastdistance = minunuseddistance;
   rollinghorizon->lastmaxdistance = *selvarmaxdistance;

   return SCIP_OKAY;
}

/** determines the graph-induced variable fixings */
static
SCIP_RETCODE determineVariableFixings(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_VAR**            fixedvars,          /**< buffer to store variables that should be fixed */
   SCIP_Real*            fixedvals,          /**< buffer to store fixing values for fixed variables */
   int*                  nfixings,           /**< pointer to store the number of fixed variables */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   ROLLINGHORIZON*       rollinghorizon,     /**< rolling horizon data structure to save relevant information, or NULL if not needed */
   SCIP_Bool*            success             /**< used to store whether the creation of the subproblem worked */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* sol;                            /* pool of solutions */
   int* distances;
   SCIP_VGRAPH* vargraph;
   SCIP_VAR* selvar;
   int nvars;
   int nbinvars;
   int nintvars;
   int fixthreshold;

   int selvarmaxdistance;

   assert(fixedvars != NULL);
   assert(fixedvals != NULL);
   assert(nfixings != NULL);

   *success = TRUE;
   *nfixings = 0;
   selvarmaxdistance = 0;
   sol = SCIPgetBestSol(scip);
   assert(sol != NULL);

   /* get required data of the original problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   /* create variable graph */
   SCIPdebugMsg(scip, "Creating variable constraint graph\n");

   /* get the saved variable graph, or create a new one */
   if( rollinghorizon != NULL )
   {
      if( rollinghorizon->niterations == 0 )
      {
         SCIP_CALL( SCIPvariableGraphCreate(scip, &rollinghorizon->variablegraph, heurdata->relaxdenseconss, 1.0 - heurdata->minfixingrate, &heurdata->nrelaxedconstraints) );
      }
      else
         assert(rollinghorizon->variablegraph != NULL);

      vargraph = rollinghorizon->variablegraph;
   }
   else
   {
      SCIP_CALL( SCIPvariableGraphCreate(scip, &vargraph, heurdata->relaxdenseconss, 1.0 - heurdata->minfixingrate, &heurdata->nrelaxedconstraints) );
   }

   /* allocate buffer memory to hold distances */
   SCIP_CALL( SCIPallocBufferArray(scip, &distances, nvars) );

   selvar = NULL;

   /* in the first iteration of the rolling horizon approach, we select an initial variable */
   if( rollinghorizon == NULL || rollinghorizon->niterations == 0 )
   {
      SCIP_CALL( selectInitialVariable(scip, heurdata, vargraph, distances, &selvar, &selvarmaxdistance) );

      /* collect and save the distances in the rolling horizon data structure */
      if( selvar != NULL && rollinghorizon != NULL )
      {
         /* collect distances in the variable graph of all variables to the selected variable */
         SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, vargraph, &selvar, 1, distances, INT_MAX, INT_MAX, INT_MAX) );
         rollingHorizonStoreDistances(scip, rollinghorizon, distances);
         rollinghorizon->lastmaxdistance = selvarmaxdistance;
      }
   }
   else
   {
      /* select the next variable, if variables are left */
      SCIP_CALL( selectNextVariable(scip, heurdata, rollinghorizon, distances, &selvar, &selvarmaxdistance) );
   }

   assert(selvar == NULL || selvarmaxdistance > 0);
   if( selvar == NULL )
   {
      *success = FALSE;
   }
   else
   {
      SCIPdebugMsg(scip, "Selected variable <%s> as central variable for a <%d>-neighborhood\n",
         SCIPvarGetName(selvar), selvarmaxdistance);

      /* fix variables that are not in the neighborhood around the selected variable */
      SCIP_CALL( fixNonNeighborhoodVariables(scip, heurdata, rollinghorizon, sol, vars, fixedvars, fixedvals, distances,
            selvarmaxdistance, nfixings) );

      fixthreshold = (int)(heurdata->minfixingrate * (heurdata->fixcontvars ? nvars : (nbinvars + nintvars)));

      /* compare actual number of fixings to limit; if we fixed not enough variables we terminate here;
       * we also terminate if no discrete variables are left
       */
      if( *nfixings < fixthreshold )
      {
         SCIPdebugMsg(scip, "Fixed %d < %d variables in gins heuristic, stopping", *nfixings, fixthreshold);
         *success = FALSE;
      }
   }

   SCIPfreeBufferArray(scip, &distances);
   if( rollinghorizon == NULL )
      SCIPvariableGraphFree(scip, &vargraph);

   return SCIP_OKAY;
}

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_HEUR*            heur,               /**< gins heuristic structure */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** set sub-SCIP solving limits */
static
SCIP_RETCODE setLimits(
   SCIP*                 subscip,            /**< SCIP data structure */
   SOLVELIMITS*          solvelimits         /**< pointer to solving limits data structure */
   )
{
   assert(subscip != NULL);
   assert(solvelimits != NULL);

   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", solvelimits->nodelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", solvelimits->timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", solvelimits->memorylimit) );

   return SCIP_OKAY;
}

/** set up the sub-SCIP */
static
SCIP_RETCODE setupSubScip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP data structure */
   SOLVELIMITS*          solvelimits,        /**< pointer to solving limits data structure */
   SCIP_HEUR*            heur                /**< this heuristic */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real cutoff;
   SCIP_Real upperbound;

   heurdata = SCIPheurGetData(heur);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* disable statistic timing inside sub SCIP */
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", heurdata->bestsollimit) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* enable conflict analysis, disable analysis of boundexceeding LPs, and restrict conflict pool */
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }

   /* speed up sub-SCIP by not checking dual LP feasibility */
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* employ a limit on the number of enforcement rounds in the quadratic constraint handlers; this fixes the issue that
    * sometimes the quadratic constraint handler needs hundreds or thousands of enforcement rounds to determine the
    * feasibility status of a single node without fractional branching candidates by separation (namely for uflquad
    * instances); however, the solution status of the sub-SCIP might get corrupted by this; hence no decutions shall be
    * made for the original SCIP
    */
   if( SCIPfindConshdlr(subscip, "quadratic") != NULL && !SCIPisParamFixed(subscip, "constraints/quadratic/enfolplimit") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "constraints/quadratic/enfolplimit", 10) );
   }

   /* add an objective cutoff */
   assert( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) );

   upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
   if( !SCIPisInfinity(scip, -1.0 * SCIPgetLowerbound(scip)) )
   {
      cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip)
                      + heurdata->minimprove * SCIPgetLowerbound(scip);
   }
   else
   {
      if( SCIPgetUpperbound(scip) >= 0 )
         cutoff = (1 - heurdata->minimprove) * SCIPgetUpperbound(scip);
      else
         cutoff = (1 + heurdata->minimprove) * SCIPgetUpperbound(scip);
   }
   cutoff = MIN(upperbound, cutoff);
   SCIP_CALL(SCIPsetObjlimit(subscip, cutoff));

   /* set solve limits for sub-SCIP */
   SCIP_CALL( setLimits(subscip, solvelimits) );

   return SCIP_OKAY;
}

/** determine limits for a sub-SCIP */
static
SCIP_RETCODE determineLimits(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< this heuristic */
   SOLVELIMITS*          solvelimits,        /**< pointer to solving limits data structure */
   SCIP_Bool*            runagain            /**< can we solve another sub-SCIP with these limits */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real maxnnodesr;
   SCIP_Longint maxnnodes;
   assert(scip != NULL);
   assert(heur != NULL);
   assert(solvelimits != NULL);
   assert(runagain != NULL);

   heurdata = SCIPheurGetData(heur);

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &solvelimits->timelimit) );
   if( !SCIPisInfinity(scip, solvelimits->timelimit) )
      solvelimits->timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &solvelimits->memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, solvelimits->memorylimit) )
   {
      solvelimits->memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      solvelimits->memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( solvelimits->timelimit <= 0.0 || solvelimits->memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
      *runagain = FALSE;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodesr = heurdata->nodesquot * SCIPgetNNodes(scip);

   /* reward gins if it succeeded often, count the setup costs for the sub-MIP as 100 nodes */
   maxnnodesr *= 1.0 + 2.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(heurdata->nsubmips + 1.0);
   maxnnodes = (SCIP_Longint)(maxnnodesr - 100.0 * heurdata->nsubmips);
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   solvelimits->nodelimit = maxnnodes - heurdata->usednodes;
   solvelimits->nodelimit = MIN(solvelimits->nodelimit, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( solvelimits->nodelimit < heurdata->minnodes )
      *runagain = FALSE;

   return SCIP_OKAY;
}

/** updates heurdata after a run of GINS */
static
void updateFailureStatistic(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   /* increase number of failures, calculate next node at which GINS should be called and update actual solutions */
   heurdata->nfailures++;
   heurdata->nextnodenumber = (heurdata->nfailures <= 25
      ? SCIPgetNNodes(scip) + 100*(2LL << heurdata->nfailures) /*lint !e703*/
      : SCIP_LONGINT_MAX);
}

#ifdef SCIP_STATISTIC
/** gets the average neighborhood size of all selected variables */
static
SCIP_Real heurdataAvgNeighborhoodSize(
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   return heurdata->sumneighborhoodvars / (MAX(1.0, (SCIP_Real)heurdata->nneighborhoods));
}

/** prints a histogram */
static
void printHistogram(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  histogram,          /**< histogram values */
   const char*           name                /**< display name of this histogram */
   )
{
   int i;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Gins: %s", name);

   /* write out entries of this histogram */
   for( i = 0; i < NHISTOGRAMBINS; ++i )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " %d", histogram[i]);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");
}
#endif

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyGins)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurGins(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeGins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitGins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->randnumgen == NULL);

   /* initialize data */
   heurdata->usednodes = 0;
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_RANDSEED) );
   heurdata->sumdiscneighborhoodvars = heurdata->sumneighborhoodvars = 0;
   heurdata->nneighborhoods = 0;
   heurdata->maxseendistance = 0;
   heurdata->nsubmips = 0;
   heurdata->nfailures = 0;
   heurdata->nextnodenumber = 0;

#ifdef SCIP_STATISTIC
   resetHistogram(heurdata->conscontvarshist);
   resetHistogram(heurdata->consdiscvarshist);
   resetHistogram(heurdata->conscontvarshist);
#endif

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEUREXIT(heurExitGins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPstatistic(
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Gins: Avg Neighborhood size: %.1f Avg. discrete neighboorhood vars: %.1f\n",
            heurdataAvgNeighborhoodSize(heurdata), heurdataAvgDiscreteNeighborhoodSize(heurdata));
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Gins: Max seen distance %d\n", heurdata->maxseendistance);
      printHistogram(scip, heurdata->consvarshist, "Constraint density histogram");
      printHistogram(scip, heurdata->conscontvarshist, "Constraint continuous density histogram");
      printHistogram(scip, heurdata->consdiscvarshist, "Constraint discrete density histogram");
      )

   SCIPfreeRandom(scip, &heurdata->randnumgen);
   heurdata->randnumgen = NULL;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecGins)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;                  /* heuristic's data */
   SCIP* subscip;                            /* the subproblem created by gins */
   SCIP_VAR** vars;                          /* original problem's variables */
   SCIP_VAR** subvars;                       /* subproblem's variables */
   SCIP_VAR** fixedvars;
   SCIP_Real* fixedvals;
   ROLLINGHORIZON* rollinghorizon;           /* data structure for rolling horizon approach */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */

   int nvars;                                /* number of original problem's variables */
   int i;
   int nfixedvars;
   SOLVELIMITS solvelimits;
   SCIP_Bool runagain;

   SCIP_Bool success;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   /* only call heuristic, if feasible solution is available */
   if( SCIPgetNSols(scip) <= 0 )
      return SCIP_OKAY;

   /* in case of many unsuccessful calls, the nextnodenumber is increased to prevent us from running too often  */
   if( SCIPgetNNodes(scip) < heurdata->nextnodenumber )
      return SCIP_OKAY;

   /* only call heuristic, if the best solution comes from transformed problem */
   assert(SCIPgetBestSol(scip) != NULL);
   if( SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   /* only call heuristic, if enough nodes were processed since last incumbent */
   if( SCIPgetNNodes(scip) - SCIPgetSolNodenum(scip,SCIPgetBestSol(scip)) < heurdata->nwaitingnodes )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if discrete variables are present */
   if( SCIPgetNBinVars(scip) == 0 && SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   runagain = TRUE;

   /* determine solving limits for the sub-SCIP for the first time */
   SCIP_CALL( determineLimits(scip, heur, &solvelimits, &runagain) );

   if( ! runagain )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, nvars) );

   /* only create a rolling horizon data structure if we need to keep it */
   if( heurdata->userollinghorizon )
      SCIP_CALL( rollingHorizonCreate(scip, &rollinghorizon) );
   else
      rollinghorizon = NULL;

   do
   {
      /* create a new problem, by fixing all variables except for a small neighborhood */
      SCIP_CALL( determineVariableFixings(scip, fixedvars, fixedvals, &nfixedvars, heurdata, rollinghorizon, &success) );

      /* terminate if it was not possible to create the subproblem */
      if( !success )
      {
         SCIPdebugMsg(scip, "Could not create the subproblem -> skip call\n");

         /* do not count this as a call to the heuristic */
         *result = SCIP_DIDNOTRUN;

         /* count this as a failure and increase the number of waiting nodes until the next call */
         updateFailureStatistic(scip, heurdata);
         goto TERMINATE;
      }

      /* initializing the subproblem */
      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
      SCIP_CALL( SCIPcreate(&subscip) );
      ++heurdata->nsubmips;

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

      /* create a problem copy as sub SCIP */
      SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "gins", fixedvars, fixedvals, nfixedvars,
            heurdata->uselprows, heurdata->copycuts, &success, NULL) );

      for( i = 0; i < nvars; i++ )
         subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

      /* free hash map */
      SCIPhashmapFree(&varmapfw);

      /* set up limits for the subproblem */
      SCIP_CALL( setupSubScip(scip, subscip, &solvelimits, heur) );

      /* solve the subproblem */
      SCIPdebugMsg(scip, "Solve Gins subMIP\n");

      /* Errors in solving the subproblem should not kill the overall solving process
       * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      SCIP_CALL_ABORT( SCIPsolve(subscip) );

      /* transfer variable statistics from sub-SCIP */
      SCIP_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );

      heurdata->usednodes += SCIPgetNNodes(subscip);

      /* check, whether a solution was found */
      if( SCIPgetNSols(subscip) > 0 )
      {
         SCIP_SOL** subsols;
         int nsubsols;

         /* check, whether a solution was found;
          * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
          */
         nsubsols = SCIPgetNSols(subscip);
         subsols = SCIPgetSols(subscip);
         success = FALSE;
         for( i = 0; i < nsubsols && !success; ++i )
         {
            SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
         }
         if( success )
            *result = SCIP_FOUNDSOL;
      }

      /* free subproblem */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );

      /* check if we want to run another rolling horizon iteration */
      runagain = (*result == SCIP_FOUNDSOL) && heurdata->userollinghorizon;
      if( runagain )
      {
         assert(rollinghorizon != NULL);
         SCIP_CALL( determineLimits(scip, heur, &solvelimits, &runagain ) );
         runagain = runagain && rollingHorizonRunAgain(scip, rollinghorizon, heurdata);
      }

   } while( runagain );

   /* delay the heuristic in case it was not successful */
   if( *result != SCIP_FOUNDSOL )
      updateFailureStatistic(scip, heurdata);

TERMINATE:

   /* only free the rolling horizon data structure if we need to keep it */
   if( heurdata->userollinghorizon )
   {
      SCIP_CALL( rollingHorizonFree(scip, &rollinghorizon) );
   }

   SCIPfreeBufferArray(scip, &fixedvals);
   SCIPfreeBufferArray(scip, &fixedvars);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the gins primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurGins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Gins primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   heurdata->randnumgen = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecGins, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyGins) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeGins) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitGins) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitGins) );

   /* add gins primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nwaitingnodes",
         "number of nodes without incumbent change that heuristic should wait",
         &heurdata->nwaitingnodes, TRUE, DEFAULT_NWAITINGNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "percentage of integer variables that have to be fixed",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, SCIPsumepsilon(scip), 1.0-SCIPsumepsilon(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which " HEUR_NAME " should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/fixcontvars",
         "should continuous variables outside the neighborhoods be fixed?",
         &heurdata->fixcontvars, TRUE, DEFAULT_FIXCONTVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/bestsollimit",
         "limit on number of improving incumbent solutions in sub-CIP",
         &heurdata->bestsollimit, FALSE, DEFAULT_BESTSOLLIMIT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxdistance",
         "maximum distance to selected variable to enter the subproblem, or -1 to select the distance "
         "that best approximates the minimum fixing rate from below",
         &heurdata->maxdistance, FALSE, DEFAULT_MAXDISTANCE, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/potential",
         "the reference point to compute the neighborhood potential: (r)oot or (p)seudo solution",
         &heurdata->potential, TRUE, DEFAULT_POTENTIAL, "pr", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/userollinghorizon",
         "should the heuristic solve a sequence of sub-MIP's around the first selected variable",
         &heurdata->userollinghorizon, TRUE, DEFAULT_USEROLLINGHORIZON, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/relaxdenseconss",
         "should dense constraints (at least as dense as 1 - minfixingrate) be ignored by connectivity graph?",
         &heurdata->relaxdenseconss, TRUE, DEFAULT_RELAXDENSECONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/rollhorizonlimfac",
         "limiting percentage for variables already used in sub-SCIPs to terminate rolling horizon approach",
         &heurdata->rollhorizonlimfac, TRUE, DEFAULT_ROLLHORIZONLIMFAC, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
