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

/**@file   heur_gins.c
 * @ingroup DEFPLUGINS_HEUR
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

#include "blockmemshell/memory.h"
#include "scip/heur_gins.h"
#include "scip/heuristics.h"
#include "scip/pub_dcmp.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include <string.h>

#define HEUR_NAME             "gins"
#define HEUR_DESC             "gins works on k-neighborhood in a variable-constraint graph"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         -1103000
#define HEUR_FREQ             20
#define HEUR_FREQOFS          8
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODESOFS      500            /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_MAXNODES      5000           /**< maximum number of nodes to regard in the subproblem */
#define DEFAULT_MINIMPROVE    0.01           /**< factor by which Gins should at least improve the incumbent */
#define DEFAULT_MINNODES      50             /**< minimum number of nodes to regard in the subproblem */
#define DEFAULT_MINFIXINGRATE 0.66           /**< minimum percentage of integer variables that have to be fixed */
#define DEFAULT_NODESQUOT     0.15           /**< subproblem nodes in relation to nodes of the original problem */
#define DEFAULT_NWAITINGNODES 100            /**< number of nodes without incumbent change that heuristic should wait */
#define DEFAULT_USELPROWS     FALSE          /**< should subproblem be created out of the rows in the LP rows,
                                              *   otherwise, the copy constructors of the constraints handlers are used */
#define DEFAULT_COPYCUTS      TRUE           /**< if DEFAULT_USELPROWS is FALSE, then should all active cuts from the
                                              *   cutpool of the original scip be copied to constraints of the subscip */
#define DEFAULT_BESTSOLLIMIT    3            /**< limit on number of improving incumbent solutions in sub-CIP */
#define DEFAULT_FIXCONTVARS FALSE            /**< should continuous variables outside the neighborhoods get fixed? */
#define DEFAULT_POTENTIAL      'r'           /**< the reference point to compute the neighborhood potential: (r)oot, (l)ocal lp, or (p)seudo solution */
#define DEFAULT_MAXDISTANCE     3            /**< maximum distance to selected variable to enter the subproblem, or -1 to
                                              *   select the distance that best approximates the minimum fixing rate from below */
#define DEFAULT_RANDSEED       71
#define DEFAULT_RELAXDENSECONSS FALSE        /**< should dense constraints (at least as dense as 1 - minfixingrate) be
                                              *   ignored by connectivity graph? */
#define DEFAULT_USEROLLINGHORIZON TRUE       /**< should the heuristic solve a sequence of sub-MIP's around the first selected variable */
#define DEFAULT_ROLLHORIZONLIMFAC  0.4       /**< limiting percentage for variables already used in sub-SCIPs to terminate rolling
                                              *   horizon approach */
#define DEFAULT_USEDECOMP    TRUE            /**< should user decompositions be considered, if available? */
#define DEFAULT_USEDECOMPROLLHORIZON FALSE   /**< should user decompositions be considered for initial selection in rolling horizon, if available? */
#define DEFAULT_USESELFALLBACK  TRUE         /**< should random initial variable selection be used if decomposition was not successful? */
#define DEFAULT_OVERLAP          0.0         /**< overlap of blocks between runs - 0.0: no overlap, 1.0: shift by only 1 block */
#define DEFAULT_CONSECUTIVEBLOCKS TRUE       /**< should blocks be treated consecutively (sorted by ascending label?) */


/*
 * Data structures
 */

/** rolling horizon data structure to control multiple LNS heuristic runs away from an original source variable */
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

/** data structure to enable GINS to solve multiple decompositions in a sequential process */
struct DecompHorizon
{
   SCIP_DECOMP*          decomp;             /**< decomposition data structure used for this horizon */
   SCIP_VAR**            vars;               /**< variables sorted by block indices */
   SCIP_SOL**            lastsolblock;       /**< last solution for which block was part of the sub-SCIP */
   SCIP_Real*            potential;          /**< potential of each block */
   int*                  blocklabels;        /**< sorted block labels of all variable blocks that satisfy the requirements */
   int*                  varblockend;        /**< block end indices in sorted variables array (position of first variable of next block) */
   int*                  ndiscretevars;      /**< number of binary and integer variables in each block */
   int*                  blockindices;       /**< block indices (from 0 to nblocks) with respect to sorting of blocks */
   int*                  nvars;              /**< number of variables (including continuous and implicit integers) in each block */
   SCIP_Bool*            suitable;           /**< TRUE if a block is suitable */
   int                   nsuitableblocks;    /**< the total number of suitable blocks */
   int                   lastblockpos;       /**< last remembered block position (in block indices, i.e., regarding sorting) */
   int                   nblocks;            /**< the number of available variable blocks, only available after initialization */
   int                   memsize;            /**< storage size of the used arrays */
   int                   varsmemsize;        /**< storage size of the vars array */
   int                   overlapinterval[2]; /**< block positions of last interval forbidden by overlap */
   SCIP_Bool             init;               /**< has the decomposition horizon been initialized? */
};
typedef struct DecompHorizon DECOMPHORIZON;

/** data structure to hold elements that are taboo */
struct TabooList
{
   int*                  taboolabels;        /**< labels or indices that are currently taboo */
   int*                  sortedlabels;       /**< array of labels in sorted order for quicker finding */
   int                   memsize;            /**< storage capacity of taboolabels */
   int                   ntaboolabels;       /**< number of elements contained in taboo list */
   SCIP_Bool             needssorting;       /**< has an element been added since the last sort? */
};
typedef struct TabooList TABOOLIST;

/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   int                   maxnodes;           /**< maximum number of nodes to regard in the subproblem */
   int                   minnodes;           /**< minimum number of nodes to regard in the subproblem */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed */
   SCIP_Real             overlap;            /**< overlap of blocks between runs - 0.0: no overlap, 1.0: shift by only 1 block */
   int                   nwaitingnodes;      /**< number of nodes without incumbent change that heuristic should wait */
   SCIP_Real             minimprove;         /**< factor by which Gins should at least improve the incumbent */
   SCIP_Longint          usednodes;          /**< nodes already used by Gins in earlier calls */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem */
   SCIP_Real             rollhorizonlimfac;  /**< limiting percentage for variables already used in sub-SCIPs to terminate
                                              *   rolling horizon approach */
   DECOMPHORIZON*        decomphorizon;      /**< decomposition horizon data structure */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator                              */
   SCIP_SOL*             lastinitsol;        /**< last solution used for selection of initial variable */
   TABOOLIST*            taboolist;          /**< taboo list of block labels that should not be used */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows? */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem? */
   SCIP_Bool             allblocksunsuitable; /**< remember if all blocks are unsuitable w.r.t. the current incumbent solution  */
   SCIP_Bool             fixcontvars;        /**< should continuous variables outside the neighborhoods get fixed? */
   int                   bestsollimit;       /**< limit on number of improving incumbent solutions in sub-CIP */
   int                   maxdistance;        /**< maximum distance to selected variable to enter the subproblem, or -1 to
                                              *   select the distance that best approximates the minimum fixing rate from below */
   int                   sumneighborhoodvars;/**< neighborhood variables sum over all seen neighborhoods */
   int                   sumdiscneighborhoodvars; /**< neighborhood discrete variables sum over all seen neighboorhoods */
   int                   nneighborhoods;     /**< number of calculated neighborhoods */
   int                   nsubmips;           /**< counter for the number of sub-MIP's that can be higher than the number of
                                              *   calls of this heuristic */
   SCIP_Bool             consecutiveblocks;  /**< should blocks be treated consecutively (sorted by ascending label?) */
   SCIP_Bool             relaxdenseconss;    /**< should dense constraints (at least as dense as 1 - minfixingrate) be
                                              *   ignored by connectivity graph? */
   SCIP_Bool             userollinghorizon;  /**< should the heuristic solve a sequence of sub-MIP's around the first
                                              *   selected variable */
   SCIP_Bool             usedecomp;          /**< should user decompositions be considered, if available? */
   SCIP_Bool             usedecomprollhorizon;/**< should user decompositions be considered for initial selection in rolling horizon, if available? */
   SCIP_Bool             useselfallback;     /**< should random initial variable selection be used if decomposition was not successful? */
   char                  potential;          /**< the reference point to compute the neighborhood potential: (r)oot or
                                              *   (p)seudo solution */
   int                   maxseendistance;    /**< maximum of all distances between two variables */
   int                   nrelaxedconstraints; /**< number of constraints that were relaxed */
   int                   nfailures;           /**< counter for the number of unsuccessful runs of this heuristic */
   SCIP_Longint          nextnodenumber;      /**< the next node number at which the heuristic should be called again */
   SCIP_Longint          targetnodes;         /**< number of target nodes, slightly increasing if (stall) node limit is hit */
};

/** represents limits for the sub-SCIP solving process */
struct SolveLimits
{
   SCIP_Longint          nodelimit;          /**< maximum number of solving nodes for the sub-SCIP */
   SCIP_Longint          stallnodelimit;     /**< limit on the number of stall nodes (nodes after last incumbent) */
};
typedef struct SolveLimits SOLVELIMITS;

/*
 * Local methods
 */

/** check if enough fixings have been found */
static
SCIP_Bool checkFixingrate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   nfixings            /**< actual number of fixings */
   )
{
   int fixthreshold;
   int nvars = SCIPgetNVars(scip);
   int nbinvars = SCIPgetNBinVars(scip);
   int nintvars = SCIPgetNIntVars(scip);
   fixthreshold = (int)(heurdata->minfixingrate * (heurdata->fixcontvars ? nvars : (nbinvars + nintvars)));

   /* compare actual number of fixings to limit; if we fixed not enough variables we terminate here;
    * we also terminate if no discrete variables are left
    */
   if( nfixings < fixthreshold )
   {
      SCIPdebugMsg(scip, "Fixed %d < %d variables in gins heuristic, stopping\n", nfixings, fixthreshold);

      return FALSE;
   }
   else
   {
      SCIPdebugMsg(scip, "Fixed enough (%d >= %d) variables in gins heuristic\n", nfixings, fixthreshold);

      return TRUE;
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
            referencepoint = varobj > 0.0 ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
            break;

         /* use root LP solution difference */
         case 'r':
            referencepoint = SCIPvarGetRootSol(var);
            break;

         /* use local LP solution */
         case 'l':
            referencepoint = SCIPgetSolVal(scip, NULL, var);
            break;
         default:
            SCIPerrorMessage("Unknown potential computation %c specified\n", heurdata->potential);
            referencepoint = 0.0;
            break;
      }

      if( SCIPisInfinity(scip, REALABS(referencepoint)) )
         continue;

      /* calculate the delta to the variables best bound */
      objdelta = (SCIPgetSolVal(scip, sol, var) - referencepoint) * varobj;
      potential += objdelta;
   }

   return potential;
}

/** (re)set overlap interval of decomposition horizon */
static
void decompHorizonSetOverlapInterval(
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure */
   int                   leftborder,         /**< left interval border */
   int                   rightborder         /**< right interval border */
   )
{
   assert(decomphorizon != NULL);
   assert(leftborder <= rightborder);

   decomphorizon->overlapinterval[0] = leftborder;
   decomphorizon->overlapinterval[1] = rightborder;
}

/** create a decomp horizon data structure */
static
SCIP_RETCODE decompHorizonCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON**       decomphorizon,      /**< pointer to store decomposition horizon */
   SCIP_DECOMP*          decomp              /**< decomposition in transformed space */
   )
{
   DECOMPHORIZON* decomphorizonptr;
   int nblocks;
   int memsize;

   assert(scip != NULL);
   assert(decomphorizon != NULL);
   assert(decomp != NULL);

   nblocks = SCIPdecompGetNBlocks(decomp);

   assert(nblocks >= 1);

   /* account an additional slot for the border */
   SCIP_CALL( SCIPallocBlockMemory(scip, decomphorizon) );
   decomphorizonptr = *decomphorizon;
   decomphorizonptr->decomp = decomp;
   decomphorizonptr->memsize = memsize = nblocks + 1;

   /* allocate storage for block related information */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decomphorizonptr->blocklabels, memsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decomphorizonptr->varblockend, memsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decomphorizonptr->suitable, memsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decomphorizonptr->ndiscretevars, memsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decomphorizonptr->nvars, memsize) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &decomphorizonptr->lastsolblock, memsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decomphorizonptr->potential, memsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &decomphorizonptr->blockindices, memsize) );

   decomphorizonptr->lastblockpos = INT_MIN; /* cannot use -1 because this is defined for linking variables */

   /* initialize data later */
   decomphorizonptr->init = FALSE;
   decomphorizonptr->vars = NULL;
   decomphorizonptr->varsmemsize = 0;

   return SCIP_OKAY;
}

/** free a decomp horizon data structure */
static
void decompHorizonFree(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON**       decomphorizon       /**< pointer to decomposition horizon that should be freed */
   )
{
   DECOMPHORIZON* decomphorizonptr;

   assert(scip != NULL);
   assert(decomphorizon != NULL);

   /* empty horizon */
   if( *decomphorizon == NULL )
      return;

   decomphorizonptr = *decomphorizon;
   SCIPfreeBlockMemoryArrayNull(scip, &decomphorizonptr->vars, decomphorizonptr->varsmemsize);

   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->blocklabels, decomphorizonptr->memsize);
   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->varblockend, decomphorizonptr->memsize);
   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->suitable, decomphorizonptr->memsize);
   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->ndiscretevars, decomphorizonptr->memsize);
   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->nvars, decomphorizonptr->memsize);
   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->lastsolblock, decomphorizonptr->memsize);
   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->potential, decomphorizonptr->memsize);
   SCIPfreeBlockMemoryArray(scip, &decomphorizonptr->blockindices, decomphorizonptr->memsize);

   SCIPfreeBlockMemory(scip, decomphorizon);

   *decomphorizon = NULL;
}

/** check if another run should be performed within the current decomposition horizon */
static
SCIP_Bool decompHorizonRunAgain(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON*        decomphorizon       /**< decomposition horizon data structure */
   )
{
   assert(scip != NULL);
   assert(decomphorizon != NULL);

   assert(decomphorizon->lastblockpos >= 0);
   assert(decomphorizon->lastblockpos < decomphorizon->nblocks);

   return TRUE;
}

/** return TRUE if the decomposition horizon has already been initialized, FALSE otherwise */
static
SCIP_Bool decompHorizonIsInitialized(
   DECOMPHORIZON*        decomphorizon       /**< decomposition horizon data structure */
   )
{
   assert(decomphorizon != NULL);

   return decomphorizon->init;
}

/** compares two block indices
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(sortIndCompDecompHorizon)
{
   DECOMPHORIZON* decomphorizon = (DECOMPHORIZON*)dataptr;
   SCIP_Real potentialbysize1;
   SCIP_Real potentialbysize2;

   assert(decomphorizon != NULL);
   assert(ind1 >= 0);
   assert(ind2 >= 0);
   assert(ind1 < decomphorizon->nblocks);
   assert(ind2 < decomphorizon->nblocks);

   if( ind1 == ind2 )
      return 0;

   /* linking variables are always sorted up front */
   if( decomphorizon->blocklabels[ind1] == SCIP_DECOMP_LINKVAR )
      return -1;
   else if( decomphorizon->blocklabels[ind2] == SCIP_DECOMP_LINKVAR )
      return 1;

   /* if one of the blocks is not suitable, return the other block */
   if( ! (decomphorizon->suitable[ind1] && decomphorizon->suitable[ind2]) )
   {
      /* prefer the suitable block; break ties based on block position */
      if( decomphorizon->suitable[ind1] )
         return -1;
      else if( decomphorizon->suitable[ind2] )
         return 1;
      else
         return ind1 - ind2;
   }

   assert(decomphorizon->suitable[ind1] && decomphorizon->suitable[ind2]);

   potentialbysize1 = decomphorizon->potential[ind1] / (SCIP_Real)(MAX(decomphorizon->ndiscretevars[ind1], 1.0));
   potentialbysize2 = decomphorizon->potential[ind2] / (SCIP_Real)(MAX(decomphorizon->ndiscretevars[ind2], 1.0));

   /* prefer the block with higher potential */
   if( potentialbysize1 > potentialbysize2 )
      return -1;
   else if( potentialbysize2 > potentialbysize1 )
      return 1;

   /* finally, prefer the block with fewer discrete variables */
   return decomphorizon->ndiscretevars[ind1] - decomphorizon->ndiscretevars[ind2];
}

/** initialize decomposition horizon data structure by setting up data structures and analyzing blocks */
static
SCIP_RETCODE decompHorizonInitialize(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** varscopy;
   int* varlabels;
   int nvars;
   int currblockstart;
   int blockpos;
   int nstblblocks;
   int ncontvarsscip;
   int b;

   SCIP_DECOMP* decomp = decomphorizon->decomp;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(! SCIPdecompIsOriginal(decomp));

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   ncontvarsscip = SCIPgetNContVars(scip) + SCIPgetNImplVars(scip);

   assert(vars != NULL);

   /* get variable labels from decomposition */
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &varscopy, vars, nvars) );

   SCIPdecompGetVarsLabels(decomp, varscopy, varlabels, nvars);

   /*  sort labels and variables */
   SCIPsortIntPtr(varlabels, (void **)varscopy, nvars);
   decomphorizon->vars = varscopy;
   decomphorizon->varsmemsize = nvars;

   blockpos = 0;
   currblockstart = 0;
   nstblblocks = 0;
   /* loop over blocks, and check if they are suitable or not for the improvement heuristic */
   while( currblockstart < nvars )
   {
      int blocklabel;
      int currblockend;
      int ndiscretevars;
      int nfixedvars;
      SCIP_Bool suitable;
      assert(blockpos < decomphorizon->memsize);

      blocklabel = varlabels[currblockstart];
      currblockend = currblockstart;
      ndiscretevars = 0;

      /* determine the block size and the variable types */
      do
      {
         if( SCIPvarGetType(varscopy[currblockend]) < SCIP_VARTYPE_IMPLINT )
            ++ndiscretevars;

         currblockend++;
      }
      while( currblockend < nvars && varlabels[currblockend] == blocklabel );

      if( heurdata->fixcontvars )
         nfixedvars = nvars - (currblockend - currblockstart);
      else
         nfixedvars = nvars - ncontvarsscip - ndiscretevars;

      suitable = nfixedvars > heurdata->minfixingrate * (heurdata->fixcontvars ? nvars : nvars - ncontvarsscip);

      decomphorizon->suitable[blockpos] = suitable;
      decomphorizon->blocklabels[blockpos] = blocklabel;
      decomphorizon->varblockend[blockpos] = currblockend;
      decomphorizon->nvars[blockpos] = currblockend - currblockstart;
      decomphorizon->ndiscretevars[blockpos] = ndiscretevars;

      currblockstart = currblockend;
      nstblblocks += (suitable ? 1 : 0);

      blockpos++;
   }

   /* not necessarily all blocks have variables; store number of available blocks */
   decomphorizon->nblocks = blockpos;
   decomphorizon->nsuitableblocks = nstblblocks;

   /* initialize block indices (identical to blockposition initially) */
   for( b = 0; b < decomphorizon->nblocks; ++b )
      decomphorizon->blockindices[b] = b;

   decompHorizonSetOverlapInterval(decomphorizon, -1, -1);

   SCIPdebugMsg(scip, "Initialized decomposition horizon for %d blocks (%d suitable)\n",
      decomphorizon->nblocks,
      decomphorizon->nsuitableblocks);

   SCIPfreeBufferArray(scip, &varlabels);

   decomphorizon->init = TRUE;

   return SCIP_OKAY;
}

/** get the first block position of the consecutive interval with the highest potential */
static
int decompHorizonGetFirstPosBestPotential(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   maxblocksize        /**< maximum block size in number of variables */
   )
{
   SCIP_SOL* bestsol;
   SCIP_Real intervalpotential;
   int b;
   int nintervalvars;
   int b1;
   int b2;
   int bestpos;
   SCIP_Real maxpotential;
   SCIP_Bool withlinkvars;
   SCIP_Bool linkvarsexist;

   assert(scip != NULL);
   assert(decomphorizon != NULL);
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);

   linkvarsexist = decomphorizon->blocklabels[0] == SCIP_DECOMP_LINKVAR;
   bestpos = 0;

   /* recompute potential of blocks */
   for( b = 0; b < decomphorizon->nblocks; ++b )
   {
      /* unsuitable blocks are left out and should not be contained in an interval */
      if( ! decomphorizon->suitable[b] )
      {
         decomphorizon->potential[b] = SCIP_REAL_MIN;
         continue;
      }

      /* store the potential of this block */
      decomphorizon->potential[b] = getPotential(scip, heurdata, bestsol,
               &decomphorizon->vars[b == 0 ? 0 : decomphorizon->varblockend[b - 1]], decomphorizon->nvars[b]);
   }

   /* sort the blocks such that the suitable blocks with the highest potential come first */
   if( ! heurdata->consecutiveblocks )
   {
      SCIPsortInd(decomphorizon->blockindices, sortIndCompDecompHorizon, (void*)decomphorizon, decomphorizon->nblocks);

      /* best potential block is now at the front (actual block position is retrieved from blockindices */
      SCIPdebugMsg(scip, "New potential based sorting with trailing block: 0 (label %d, potential %.4g)\n",
         decomphorizon->blocklabels[decomphorizon->blockindices[0]], decomphorizon->potential[decomphorizon->blockindices[0]]);

      return 0;
   }

   /* compute the consecutive blocks interval with largest potential */
   b1 = linkvarsexist ? 0 : -1;
   b2 = -1;
   nintervalvars = 0;
   intervalpotential = 0.0;
   maxpotential = 0.0;
   withlinkvars = FALSE;

   while( b1 < decomphorizon->nblocks - 1 )
   {
      int blockindb1;
      int blockindb2;
      ++b1;
      blockindb1 = decomphorizon->blockindices[b1];

      if( ! decomphorizon->suitable[decomphorizon->blockindices[b1]] )
      {
         nintervalvars = 0;
         intervalpotential = 0.0;
         withlinkvars = FALSE;
         b2 = b1;
         continue;
      }

      /* interval starts at b1 */
      if( b2 < b1 )
      {
         nintervalvars = decomphorizon->ndiscretevars[blockindb1];
         assert(nintervalvars < maxblocksize); /* otherwise, it wasn't suitable */
         intervalpotential = decomphorizon->potential[blockindb1];
         withlinkvars = FALSE;
         b2 = b1;
      }
      /* subtract the variables from the previous block */
      else
      {
         int prevblockind;
         assert(b1 > (linkvarsexist ? 1 : 0));
         prevblockind = decomphorizon->blockindices[b1 - 1];
         assert(decomphorizon->suitable[prevblockind]);
         nintervalvars -= decomphorizon->ndiscretevars[prevblockind];
         intervalpotential -= decomphorizon->potential[prevblockind];
      }

      /* check if block allows to include linking variables */
      if( ! withlinkvars && linkvarsexist && decomphorizon->ndiscretevars[0] + decomphorizon->ndiscretevars[blockindb1] < maxblocksize )
      {
         withlinkvars = TRUE;
         nintervalvars = decomphorizon->ndiscretevars[0] + decomphorizon->ndiscretevars[blockindb1];
         b2 = b1;
      }
      else if( withlinkvars && decomphorizon->ndiscretevars[0] + decomphorizon->ndiscretevars[blockindb1] >= maxblocksize )
      {
         withlinkvars = FALSE;
         nintervalvars = decomphorizon->ndiscretevars[blockindb1];
         b2 = b1;
      }

      /* extend the interval by further blocks, if possible */
      while( ++b2 < decomphorizon->nblocks )
      {
         blockindb2 = decomphorizon->blockindices[b2];

         if( ! decomphorizon->suitable[blockindb2] || nintervalvars + decomphorizon->ndiscretevars[blockindb2] >= maxblocksize )
            break;

         nintervalvars += decomphorizon->ndiscretevars[blockindb2];
         intervalpotential += decomphorizon->potential[blockindb2];
      }

      /* store the start of the interval with maximum potential */
      if( intervalpotential > maxpotential )
      {
         bestpos = b1; /* because pos is incremented by 1 again */
         maxpotential = intervalpotential;
      }
   }

   SCIPdebugMsg(scip, "Potential based selection chooses interval starting from block %d with potential %.1g\n",
      bestpos, maxpotential);

   return bestpos;
}

/** has this block been used too recently? */
static
SCIP_Bool decompHorizonBlockUsedRecently(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure */
   int                   blockpos            /**< position of block */
   )
{
   assert(scip != NULL);
   assert(decomphorizon != NULL);
   assert(0 <= blockpos);
   assert(blockpos < decomphorizon->nblocks);

   return (decomphorizon->lastsolblock[decomphorizon->blockindices[blockpos]] == SCIPgetBestSol(scip) ||
      (decomphorizon->overlapinterval[0] <= blockpos && blockpos <= decomphorizon->overlapinterval[1]));
}

/** query the start and end of the next suitable block after the last @p lastblockused
 *
 *  @return TRUE if next suitable block could be found, otherwise FALSE
 */
static
SCIP_Bool decompHorizonNext(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   maxblocksize,       /**< maximum block size in number of variables */
   int*                  blockintervalstart, /**< pointer to store position of first block of interval */
   int*                  blockintervalend,   /**< pointer to store position of last block of interval */
   int*                  nextblocklabel,     /**< pointer to store label of the next suitable block */
   SCIP_Bool*            fixlinkvars         /**< should the linking variables be fixed, as well? */
   )
{
   SCIP_Bool found;
   int pos;
   int firstpos;
   SCIP_SOL* bestsol;
   assert(decomphorizon != NULL);
   assert(blockintervalstart != NULL);
   assert(blockintervalend != NULL);
   assert(nextblocklabel != NULL);
   assert(fixlinkvars != NULL);

   assert(decomphorizon->init);

   if( decomphorizon->nsuitableblocks == 0 )
   {
      return FALSE;
   }

   /* get the last block position that was used by the heuristic. Search for it, and continue with the next block. */
   found = decomphorizon->lastblockpos >= 0;
   firstpos = decomphorizon->lastblockpos;
   assert(! found || (firstpos >= 0 && firstpos < decomphorizon->nblocks));

   bestsol = SCIPgetBestSol(scip);

   /* choose first position based on potential; subtract -1 because we immediately increase it */
   if( ! found || bestsol != decomphorizon->lastsolblock[decomphorizon->blockindices[firstpos]] )
      firstpos = decompHorizonGetFirstPosBestPotential(scip, decomphorizon, heurdata, maxblocksize) - 1;

   /* that's why we subtract 1 from potential based position computation */
   pos = (firstpos + 1) % decomphorizon->nblocks;

   while( pos < decomphorizon->nblocks &&
      (! decomphorizon->suitable[decomphorizon->blockindices[pos]] || decomphorizon->blocklabels[decomphorizon->blockindices[pos]] == SCIP_DECOMP_LINKVAR ||
         decompHorizonBlockUsedRecently(scip, decomphorizon, pos)) )
      pos++;

   if( pos == decomphorizon->nblocks )
   {
      pos = 0;
      while( pos < firstpos &&
         (! decomphorizon->suitable[decomphorizon->blockindices[pos]] || decomphorizon->blocklabels[decomphorizon->blockindices[pos]] == SCIP_DECOMP_LINKVAR ||
            decompHorizonBlockUsedRecently(scip, decomphorizon, pos)) )
         pos++;
   }

   assert(pos == firstpos || (0 <= pos && decomphorizon->nblocks > pos && (decomphorizon->suitable[decomphorizon->blockindices[pos]] || pos == 0)));

   *fixlinkvars = TRUE;
   /* the next suitable block position has been discovered */
   if( pos != firstpos && decomphorizon->suitable[decomphorizon->blockindices[pos]] && !decompHorizonBlockUsedRecently(scip, decomphorizon, pos) )
   {
      int ndiscretevars;
      *nextblocklabel = decomphorizon->blocklabels[decomphorizon->blockindices[pos]];
      *blockintervalstart = pos;
      *blockintervalend = pos;

      ndiscretevars = decomphorizon->ndiscretevars[decomphorizon->blockindices[pos]];
      /* check if linking variable block exceeds maximum block size */
      if( decomphorizon->blocklabels[0] == SCIP_DECOMP_LINKVAR )
      {
         *fixlinkvars = decomphorizon->ndiscretevars[0] + ndiscretevars > maxblocksize;
      }

      /* add linking variables to the block */
      if( !(*fixlinkvars) )
         ndiscretevars += decomphorizon->ndiscretevars[0];

      /* extend the subproblem until maximum target fixing rate is reached */
      while( ++pos < decomphorizon->nblocks && decomphorizon->suitable[decomphorizon->blockindices[pos]] && ndiscretevars + decomphorizon->ndiscretevars[decomphorizon->blockindices[pos]] < maxblocksize )
      {
         *blockintervalend = pos;
         *nextblocklabel = decomphorizon->blocklabels[decomphorizon->blockindices[pos]];
         ndiscretevars += decomphorizon->ndiscretevars[decomphorizon->blockindices[pos]];
      }

      return TRUE;
   }
   else
   {
      /* no next, suitable block exists */
      *blockintervalstart = *blockintervalend = -1;
      *nextblocklabel = -100;

      return FALSE;
   }
}

/** get the variables of this decomposition horizon */
static
SCIP_VAR** decomphorizonGetVars(
   DECOMPHORIZON*        decomphorizon       /**< decomposition horizon data structure */
   )
{
   assert(decomphorizon != NULL);
   assert(decomphorizon->init);

   return decomphorizon->vars;
}

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


/** reset a taboo list */
static
void tabooListReset(
   TABOOLIST*            taboolist           /**< taboo list data structure */
   )
{
   taboolist->ntaboolabels = 0;
   taboolist->needssorting = FALSE;
}

/** create a taboo list data structure */
static
SCIP_RETCODE createTabooList(
   SCIP*                 scip,               /**< SCIP data structure */
   TABOOLIST**           taboolist,          /**< pointer to store taboo list data structure */
   int                   initialsize         /**< initial storage capacity of taboo list */
   )
{
   assert(scip != NULL);
   assert(taboolist != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, taboolist) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*taboolist)->taboolabels, initialsize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*taboolist)->sortedlabels, initialsize) );
   (*taboolist)->memsize = initialsize;
   tabooListReset(*taboolist);

   return SCIP_OKAY;
}

/** free a taboo list data structure */
static
void freeTabooList(
   SCIP*                 scip,               /**< SCIP data structure */
   TABOOLIST**           taboolist           /**< pointer to taboo list data structure that should be freed */
   )
{
   assert(scip != NULL);
   assert(taboolist != NULL);

   if( *taboolist == NULL )
      return;

   SCIPfreeBlockMemoryArray(scip, &(*taboolist)->sortedlabels, (*taboolist)->memsize);
   SCIPfreeBlockMemoryArray(scip, &(*taboolist)->taboolabels, (*taboolist)->memsize);
   SCIPfreeBlockMemory(scip, taboolist);
}

/** add an element to the taboo list */
static
SCIP_RETCODE tabooListAdd(
   SCIP*                 scip,               /**< SCIP data structure */
   TABOOLIST*            taboolist,          /**< taboo list data structure */
   int                   elem                /**< element that should be added to the taboo list */
   )
{
   assert(scip != NULL);
   assert(taboolist != NULL);

   if( taboolist->memsize == taboolist->ntaboolabels )
   {
      int newsize = SCIPcalcMemGrowSize(scip, taboolist->ntaboolabels + 1);
      assert(newsize >= taboolist->ntaboolabels + 1);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &taboolist->taboolabels, taboolist->memsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &taboolist->sortedlabels, taboolist->memsize, newsize) );

      taboolist->memsize = newsize;
   }

   assert(taboolist->ntaboolabels < taboolist->memsize);
   taboolist->taboolabels[taboolist->ntaboolabels++] = elem;

   taboolist->needssorting = TRUE;

   return SCIP_OKAY;
}

/** find an element in the taboo list */
static
SCIP_Bool tabooListFind(
   TABOOLIST*            taboolist,          /**< taboo list data structure */
   int                   elem                /**< element that should be added to the taboo list */
   )
{
   SCIP_Bool found;
   int pos;
   assert(taboolist != NULL);

   if( taboolist->ntaboolabels == 0 )
      return FALSE;

   if( taboolist->needssorting )
   {
      BMScopyMemoryArray(taboolist->sortedlabels, taboolist->taboolabels, taboolist->ntaboolabels);
      SCIPsortInt(taboolist->sortedlabels, taboolist->ntaboolabels);
   }

   found = SCIPsortedvecFindInt(taboolist->sortedlabels, elem, taboolist->ntaboolabels, &pos);

   return found;
}

/** get most recent k elements from taboo list */
static
int* tabooListGetLastK(
   TABOOLIST*            taboolist,          /**< taboo list data structure */
   int                   k                   /**< the number of elements that should be returned */
   )
{
   assert(taboolist != NULL);
   assert(k <= taboolist->ntaboolabels);
   assert(k > 0);

   return &taboolist->taboolabels[taboolist->ntaboolabels - k];
}

/** get number of elements in taboo list */
static
int taboolistgetNElems(
   TABOOLIST*            taboolist           /**< taboo list data structure */
   )
{
   return taboolist->ntaboolabels;
}

/** free a rolling horizon data structure */
static
void rollingHorizonFree(
   SCIP*                 scip,               /**< SCIP data structure */
   ROLLINGHORIZON**      rollinghorizon      /**< pointer to rolling horizon data structure */
   )
{
   assert(scip != NULL);
   assert(rollinghorizon != NULL);

   /* empty rolling horizon */
   if( *rollinghorizon == NULL )
      return;

   if( (*rollinghorizon)->variablegraph != NULL )
   {
      SCIPvariableGraphFree(scip, &(*rollinghorizon)->variablegraph);
   }

   SCIPfreeBlockMemoryArray(scip, &(*rollinghorizon)->distances, (*rollinghorizon)->distancessize);
   SCIPfreeBlockMemoryArray(scip, &(*rollinghorizon)->used, (*rollinghorizon)->distancessize);
   SCIPfreeBlockMemory(scip, rollinghorizon);
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

/** get fixing value of variable */
static
SCIP_Real getFixVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution in main SCIP for fixing values */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real fixval;
   SCIP_Real lb;
   SCIP_Real ub;

   fixval = SCIPgetSolVal(scip, sol, var);
   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   assert(SCIPisLE(scip, lb, ub));

   /* due to dual reductions, it may happen that the solution value is not in the variable's domain anymore */
   if( SCIPisLT(scip, fixval, lb) )
      fixval = lb;
   else if( SCIPisGT(scip, fixval, ub) )
      fixval = ub;

   return fixval;
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
         SCIP_Real fixval;

         fixval = getFixVal(scip, sol, vars[i]);

         /* store variable and value of this fixing */
         if( !SCIPisInfinity(scip, REALABS(fixval)) )
         {
            fixedvars[*nfixings] = vars[i];
            fixedvals[*nfixings] = fixval;
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

/** count occurrences of this label */
static
int countLabel(
   int*                  labels,             /**< sorted array of labels */
   int                   start,              /**< start position */
   int                   nlabels             /**< length of the labels array */
   )
{
   int label = labels[start];
   int end = start;

   assert(labels != NULL);
   assert(start < nlabels);
   assert(start >= 0);

   do
   {
      ++end;
   }
   while( end < nlabels && labels[end] == label );

   return end - start;
}

/** todo select initial variable based on decomposition information, if available */
static
SCIP_RETCODE selectInitialVariableDecomposition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure with variable labels */
   SCIP_VGRAPH*          vargraph,           /**< variable graph data structure to work on */
   int*                  distances,          /**< breadth-first distances indexed by Probindex of variables */
   SCIP_VAR**            selvar,             /**< pointer to store the selected variable */
   int*                  selvarmaxdistance   /**< maximal distance k to consider for selected variable neighborhood */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_VAR** varscopy;
   int* varlabels;
   int* discvaridxs;
   SCIP_Real bestpotential;
   int nbinvars;
   int nintvars;
   int nvars;
   int currblockstart;
   int bestvaridx;
   int cntmessages;
   int nblocks;
   TABOOLIST* taboolist;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(decomp != NULL);
   assert(vargraph != NULL);
   assert(distances != NULL);
   assert(selvar != NULL);
   assert(selvarmaxdistance != NULL);

   sol = SCIPgetBestSol(scip);
   assert(sol != NULL);
   nblocks = SCIPdecompGetNBlocks(decomp);
   /* create a taboo list for recently used block labels at the first initial variable selection */
   if( heurdata->taboolist == NULL )
   {
      SCIPdebugMsg(scip, "Creating taboo list\n");
      SCIP_CALL( createTabooList(scip, &heurdata->taboolist, nblocks) );
   }

   taboolist = heurdata->taboolist;
   assert(taboolist != NULL);

   /* reset taboo list if a new solution has been found since the last initialization call */
   if( sol != heurdata->lastinitsol )
   {
      int nblockstokeep = 3;
      int e;
      int ntaboolistelems;
      ntaboolistelems = taboolistgetNElems(heurdata->taboolist);

      /* keep the last 3 blocks except for border cases of very coarse decomposition, or too few taboo elements */
      if( taboolistgetNElems(heurdata->taboolist) > 0 )
      {
         nblockstokeep = MIN(nblockstokeep, nblocks - 1);
         nblockstokeep = MIN(nblockstokeep, MAX(1, ntaboolistelems - 1));
         nblockstokeep = MAX(nblockstokeep, 0);
      }
      else
         nblockstokeep = 0;

      SCIPdebugMsg(scip, "Resetting taboo list, keeping %d elements\n", nblockstokeep);
      if( nblockstokeep > 0 )
      {
         int* labelstokeep;
         int* taboolistlastk;
         taboolistlastk = tabooListGetLastK(taboolist, nblockstokeep);
         SCIP_CALL( SCIPduplicateBufferArray(scip, &labelstokeep, taboolistlastk, nblockstokeep) );

         tabooListReset(taboolist);

         /* reinstall the last 3 elements into the taboo list */
         for( e = 0; e < nblockstokeep; ++e )
         {
            SCIP_CALL( tabooListAdd(scip, taboolist, labelstokeep[e]) );
         }

         SCIPfreeBufferArray(scip, &labelstokeep);
      }
      else if( ntaboolistelems > 0 )
      {
         tabooListReset(taboolist);
      }

      heurdata->allblocksunsuitable = FALSE;
   }

   *selvar = NULL;
   /* do not continue if, for this solution, all blocks are known to be unsuitable */
   if( heurdata->allblocksunsuitable )
   {
      SCIPdebugMsg(scip, "Skip initial variable selection because all blocks are unsuitable for solution %d\n",
         SCIPsolGetIndex(sol));
      return SCIP_OKAY;
   }

   /* get integer and binary variables from problem and labels for all variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   SCIP_CALL( SCIPduplicateBufferArray(scip, &varscopy, vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &discvaridxs, nvars) );

   SCIPdecompGetVarsLabels(decomp, vars, varlabels, nvars);

   /* sort the variables copy by the labels */
   SCIPsortIntPtr(varlabels, (void **)varscopy, nvars);

   currblockstart = 0;
   bestpotential = 0.0;
   bestvaridx = -1;
   cntmessages = 0;
   /* compute the potential for every block */
   while( currblockstart < nvars )
   {
      int currblockend;
      int v;
      int label = varlabels[currblockstart];
      int ndiscblockvars;
      SCIP_Real potential;

      currblockend = currblockstart + countLabel(varlabels, currblockstart, nvars);

      /* this block was recently used and should be skipped */
      if( tabooListFind(taboolist, label) )
      {
         if( cntmessages++ < 3 )
            SCIPdebugMsg(scip, "Skipping block <%d> from taboo list\n", label);

         currblockstart += currblockend;

         continue;
      }

      /* author bzfhende
       *
       * TODO omit the linking variables from the computation of the potential?
       */
      /* check if block has discrete variables */
      ndiscblockvars = 0;
      for( v = currblockstart; v < currblockend; ++v )
      {
         if( SCIPvarGetType(varscopy[v]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(varscopy[v]) == SCIP_VARTYPE_INTEGER )
            discvaridxs[ndiscblockvars++] = v;
      }

      /* skip potential computation if block has no discrete variables */
      if( ndiscblockvars > 0 )
      {
         potential = getPotential(scip, heurdata, sol, &(varscopy[currblockstart]), currblockend - currblockstart);

         if( potential > bestpotential )
         {
            bestpotential = potential;
            /* randomize the selection of the best variable */
            bestvaridx = discvaridxs[SCIPrandomGetInt(heurdata->randnumgen, 0, ndiscblockvars - 1)];
            assert(bestvaridx >= 0);
         }
      }

      currblockstart += currblockend;
   }

   /* we return the first discrete variable from the block with maximum potential */
   if( bestvaridx >= 0 )
   {
      *selvar = varscopy[bestvaridx];

      /* save block label in taboo list to not use it again too soon */
      SCIP_CALL( tabooListAdd(scip, taboolist, varlabels[bestvaridx]) );

      SCIPdebugMsg(scip, "Select initial variable <%s> from block <%d>\n", SCIPvarGetName(*selvar), varlabels[bestvaridx]);
   }
   else
   {
      SCIPdebugMsg(scip, "Could not find suitable block for variable selection.\n");
      heurdata->allblocksunsuitable = TRUE;
      *selvar = NULL;
   }

   /* use the variable constraint graph to compute distances to all other variables, and store the selvarmaxdistance */
   if( *selvar != NULL )
   {
      SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, vargraph, selvar, 1, distances,
            heurdata->maxdistance == -1 ? INT_MAX : heurdata->maxdistance, INT_MAX, INT_MAX) );

      SCIP_CALL( determineMaxDistance(scip, heurdata, distances, selvarmaxdistance) );

      /* maximum distance is 0, i.e., even the size of the 1-neighborhood of this variable exceeds the fixing rate */
      if( *selvarmaxdistance == 0 )
      {
         SCIPdebugMsg(scip, "1-Neighborhood of variable <%s> too large.\n", SCIPvarGetName(*selvar));
         *selvar = NULL;
      }
   }

   /* remember this solution for the next initial selection */
   heurdata->lastinitsol = sol;

   SCIPfreeBufferArray(scip, &discvaridxs);
   SCIPfreeBufferArray(scip, &varlabels);
   SCIPfreeBufferArray(scip, &varscopy);

   return SCIP_OKAY;
}



/** select a good starting variable at the first iteration of a rolling horizon approach */
static
SCIP_RETCODE selectInitialVariableRandomly(
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
      while( SCIPvarGetProbindex(choosevar) < 0 && nintegralvarsleft > 0 );

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
   }
   while( rollingHorizonRunAgain(scip, rollinghorizon, heurdata) && (*selvar == NULL || *selvarmaxdistance == 0) );

   /* breadth-first search determines the distances of all variables
    * that are no more than maxdistance away from the start variable
    */
   assert(*selvarmaxdistance <= rollinghorizon->lastmaxdistance);
   *selvarmaxdistance = MIN(*selvarmaxdistance, rollinghorizon->lastmaxdistance);
   rollinghorizon->lastdistance = minunuseddistance;
   rollinghorizon->lastmaxdistance = *selvarmaxdistance;

   return SCIP_OKAY;
}

/** mark some of the blocks between currblockstart and currblockend as recently used, depending on overlap */
static
void decompHorizonMarkInterval(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_SOL*             sol,                /**< solution by which some of the blocks should be marked */
   int                   blockstartpos,      /**< current start position of interval */
   int                   blockendpos         /**< current end position (inclusive) of interval */
   )
{
   int nvarsinterval;
   int nvarsstartofinterval;
   int solstamppos;
   int b;
   SCIP_Real overlapcomplement;

   assert(decomphorizon != NULL);
   assert(heurdata != NULL);

   /* is the next block unsuitable or have we reached the end of the blocks? In those cases,
    * we mark all blocks of the current interval; we hence avoid to rerun on a subset of the current subproblem
    * in the next iteration; we achieve this by setting the overlap to 0.0, (its complement to 1.0)
    * such that all blocks are marked
    */
   if( blockendpos == decomphorizon->nblocks - 1 || ! decomphorizon->suitable[decomphorizon->blockindices[blockendpos + 1]] )
      overlapcomplement = 1.0;
   else
      overlapcomplement = 1.0 - heurdata->overlap;

   /* count the total number of variables in the subproblem defining blocks */
   nvarsinterval = 0;
   for( b = blockstartpos; b <= blockendpos; ++b )
      nvarsinterval += decomphorizon->ndiscretevars[decomphorizon->blockindices[b]];

   nvarsstartofinterval = 0;
   /* stamp the first blocks up to the desired overlap by the current incumbent solution */
   for( solstamppos = blockstartpos; solstamppos <= blockendpos; ++solstamppos )
   {
      decomphorizon->lastsolblock[decomphorizon->blockindices[solstamppos]] = sol;
      nvarsstartofinterval += decomphorizon->ndiscretevars[decomphorizon->blockindices[solstamppos]];

      if( nvarsstartofinterval >= overlapcomplement * nvarsinterval )
         break;
   }
   decomphorizon->lastblockpos = solstamppos;
   SCIPdebugMsg(scip, "Blocks %d (label %d)-- %d (label %d) marked with incumbent solution\n",
      blockstartpos, decomphorizon->blocklabels[decomphorizon->blockindices[blockstartpos]],
      solstamppos, decomphorizon->blocklabels[decomphorizon->blockindices[solstamppos]]);

   /* remember the blocks up to the found position as most recent overlap interval */
   decompHorizonSetOverlapInterval(decomphorizon, blockstartpos, solstamppos);
}

/** determine the variable fixings based on a decomposition */
static
SCIP_RETCODE determineVariableFixingsDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure */
   SCIP_VAR**            fixedvars,          /**< buffer to store variables that should be fixed */
   SCIP_Real*            fixedvals,          /**< buffer to store fixing values for fixed variables */
   int*                  nfixings,           /**< pointer to store the number of fixed variables */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool*            success             /**< used to store whether the creation of the subproblem worked */
   )
{
   SCIP_SOL* sol;
   SCIP_Bool hasnext;
   SCIP_Bool fixlinkvars;
   int currblockstart;
   int currblockend;
   int currblocklabel;
   int maxblocksize;

   assert(scip != NULL);
   assert(decomphorizon != NULL);

   sol = SCIPgetBestSol(scip);

   /* initialize the decomposition horizon first for the current variables */
   if( ! decompHorizonIsInitialized(decomphorizon) )
   {
      SCIP_CALL( decompHorizonInitialize(scip, decomphorizon, heurdata) );
   }

   maxblocksize = (int)((1.0 - heurdata->minfixingrate) * (SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip))) - 1;

   /* query the next suitable interval of blocks */
   hasnext = decompHorizonNext(scip, decomphorizon, heurdata, maxblocksize,
            &currblockstart, &currblockend, &currblocklabel, &fixlinkvars);

   if( ! hasnext )
   {
      SCIPdebugMsg(scip, "Could not find a suitable interval that follows %d\n",
               decomphorizon->lastblockpos);

      *success = FALSE;
   }
   else
   {
      /* fix all discrete/continuous variables that are not part of this interval */
      SCIP_VAR** vars;
      int v;
      int startblockposs[] = {fixlinkvars ? 0 : 1, currblockend + 1};
      int endblockposs[] = {currblockstart, decomphorizon->nblocks};
      int p;
      int b;

      SCIPdebugMsg(scip, "Fix %s variables (%scluding linking variables) except blocks %d (label %d) -- %d (label %d)\n",
         heurdata->fixcontvars ? "all" : "discrete",
         fixlinkvars ? "in" : "ex",
         currblockstart, decomphorizon->blocklabels[decomphorizon->blockindices[currblockstart]],
         currblockend, decomphorizon->blocklabels[decomphorizon->blockindices[currblockend]]);

      vars = decomphorizonGetVars(decomphorizon);

      /* iterate over the two blocks left and right of the selected consecutive interval and fix all variables
       *
       * 0, ... b, ... ,[currblockstart, ..., currblockend], currblockend + 1, ..., decomphorizon->nblocks - 1
       */
      for( p = 0; p < 2; ++p )
      {
         /* iterate over all blocks and fix those outside of the blocks interval that defines the current subproblem */
         for( b = startblockposs[p]; b < endblockposs[p]; ++b )
         {
            int blockind = decomphorizon->blockindices[b];
            int varstartpos = blockind == 0 ? 0 : decomphorizon->varblockend[blockind - 1];
            int varendpos = decomphorizon->varblockend[blockind];

            /* fix variables inside of this block */
            for( v = varstartpos; v < varendpos; ++v )
            {
               SCIP_VAR* var = vars[v];

               if( heurdata->fixcontvars || SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
               {
                  SCIP_Real fixval;

                  fixval = getFixVal(scip, sol, var);

                  /* store variable and value of this fixing */
                  if( !SCIPisInfinity(scip, REALABS(fixval)) )
                  {
                     fixedvars[*nfixings] = var;
                     fixedvals[*nfixings] = fixval;
                     ++(*nfixings);
                  }
               }
            }
         }
      }

      *success = checkFixingrate(scip, heurdata, *nfixings);

      decompHorizonMarkInterval(scip, decomphorizon, heurdata, sol, currblockstart, currblockend);
   }

   return SCIP_OKAY; /*lint !e438*/
}

/** choose a decomposition from the store or return NULL if none exists/no decomposition was suitable */
static
SCIP_DECOMP* chooseDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DECOMP** decomps;
   int ndecomps;
   int currdecomp;

   /* TODO coming soon: better selection than last nontrivial decomposition that has been input */
   SCIPgetDecomps(scip, &decomps, &ndecomps, FALSE);
   currdecomp = ndecomps;

   while( --currdecomp >= 0 )
   {
      if( SCIPdecompGetNBlocks(decomps[currdecomp]) > 0 )
         return decomps[currdecomp];
   }

   return NULL;
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
   DECOMPHORIZON*        decomphorizon,      /**< decomposition horizon data structure, or NULL */
   SCIP_Bool*            success             /**< used to store whether the creation of the subproblem worked */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* sol;
   int* distances;
   SCIP_VGRAPH* vargraph;
   SCIP_VAR* selvar;
   int nvars;
   int nbinvars;
   int nintvars;

   int selvarmaxdistance;

   assert(fixedvars != NULL);
   assert(fixedvals != NULL);
   assert(nfixings != NULL);

   *success = TRUE;
   *nfixings = 0;
   selvarmaxdistance = 0;
   sol = SCIPgetBestSol(scip);
   assert(sol != NULL);

   /* determine the variable fixings based on latest user decomposition */
   if( decomphorizon != NULL )
   {
      SCIP_CALL( determineVariableFixingsDecomp(scip, decomphorizon, fixedvars, fixedvals, nfixings, heurdata, success) );

      /* do not use fallback strategy if user parameter does not allow it */
      if( *success || ! heurdata->useselfallback )
         return SCIP_OKAY;
   }

   /* get required data of the source problem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   /* get the saved variable graph, or create a new one */
   if( rollinghorizon != NULL )
   {
      if( rollinghorizon->niterations == 0 )
      {
         /* create variable graph */
         SCIPdebugMsg(scip, "Creating variable constraint graph\n");
         SCIP_CALL( SCIPvariableGraphCreate(scip, &rollinghorizon->variablegraph, heurdata->relaxdenseconss, 1.0 - heurdata->minfixingrate, &heurdata->nrelaxedconstraints) );
      }
      else
         assert(rollinghorizon->variablegraph != NULL);

      vargraph = rollinghorizon->variablegraph;
   }
   else
   {
      /* create variable graph */
      SCIPdebugMsg(scip, "Creating variable constraint graph\n");
      SCIP_CALL( SCIPvariableGraphCreate(scip, &vargraph, heurdata->relaxdenseconss, 1.0 - heurdata->minfixingrate, &heurdata->nrelaxedconstraints) );
   }

   /* allocate buffer memory to hold distances */
   SCIP_CALL( SCIPallocBufferArray(scip, &distances, nvars) );

   selvar = NULL;

   /* in the first iteration of the rolling horizon approach, we select an initial variable */
   if( rollinghorizon == NULL || rollinghorizon->niterations == 0 )
   {
      SCIP_Bool userandomselection = TRUE;

      /* choose the initial variable based on a user decomposition, if available */
      if( heurdata->usedecomprollhorizon )
      {
         SCIP_DECOMP* decomp = chooseDecomp(scip);
         if( decomp != NULL )
         {
            SCIP_CALL( selectInitialVariableDecomposition(scip, heurdata, decomp, vargraph,
                  distances, &selvar, &selvarmaxdistance) );

            userandomselection = (selvar == NULL && heurdata->useselfallback);
         }
      }

      assert(selvar == NULL || !userandomselection);
      /* use random variable selection as fallback strategy, if no user decomposition is available, or the
       * heuristic should not use decomposition
       */
      if( userandomselection )
      {
         SCIP_CALL( selectInitialVariableRandomly(scip, heurdata, vargraph, distances, &selvar, &selvarmaxdistance) );
      }

      /* collect and save the distances in the rolling horizon data structure */
      if( selvar != NULL && rollinghorizon != NULL )
      {
         /* collect distances in the variable graph of all variables to the selected variable */
         SCIP_CALL( SCIPvariablegraphBreadthFirst(scip, vargraph, &selvar, 1, distances, INT_MAX, INT_MAX, INT_MAX) );
         rollingHorizonStoreDistances(rollinghorizon, distances);
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

     *success = checkFixingrate(scip, heurdata, *nfixings);
   }

   SCIPfreeBufferArray(scip, &distances);
   if( rollinghorizon == NULL )
      SCIPvariableGraphFree(scip, &vargraph);

   return SCIP_OKAY;
}

/** set sub-SCIP solving limits */
static
SCIP_RETCODE setLimits(
   SCIP*                 sourcescip,         /**< SCIP data structure */
   SCIP*                 subscip,            /**< SCIP data structure */
   SOLVELIMITS*          solvelimits         /**< pointer to solving limits data structure */
   )
{
   assert(sourcescip != NULL);
   assert(subscip != NULL);
   assert(solvelimits != NULL);

   SCIP_CALL( SCIPcopyLimits(sourcescip, subscip) );

   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", solvelimits->nodelimit) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", solvelimits->stallnodelimit) );

   /* safe long running times caused by lack of dual convergence */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/gap", 0.01) );

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
   SCIP_CALL( setLimits(scip, subscip, solvelimits) );

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
   SCIP_Real confidence;
   SCIP_Longint maxnnodes;
   SCIP_Bool copylimits;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(solvelimits != NULL);
   assert(runagain != NULL);

   heurdata = SCIPheurGetData(heur);

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &copylimits) );

   if( ! copylimits )
      *runagain = FALSE;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   maxnnodesr = heurdata->nodesquot * SCIPgetNNodes(scip);

   /* reward gins if it succeeded often, count the setup costs for the sub-MIP as 100 nodes */
   confidence = (SCIP_Real)SCIPheurGetNCalls(heur);
   confidence = confidence / (confidence + 5.0);
   maxnnodesr *= 1.0 + confidence * 2.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(heurdata->nsubmips + 1.0);
   maxnnodes = (SCIP_Longint)(maxnnodesr - 100.0 * heurdata->nsubmips);
   maxnnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   solvelimits->nodelimit = maxnnodes - heurdata->usednodes;
   solvelimits->nodelimit = MIN(solvelimits->nodelimit, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( solvelimits->nodelimit < heurdata->targetnodes )
      *runagain = FALSE;

   solvelimits->stallnodelimit = heurdata->targetnodes;

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
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, DEFAULT_RANDSEED, TRUE) );
   heurdata->sumdiscneighborhoodvars = heurdata->sumneighborhoodvars = 0;
   heurdata->nneighborhoods = 0;
   heurdata->maxseendistance = 0;
   heurdata->nsubmips = 0;
   heurdata->nfailures = 0;
   heurdata->nextnodenumber = 0;
   heurdata->lastinitsol = NULL;
   heurdata->allblocksunsuitable = FALSE;
   heurdata->taboolist = NULL;
   heurdata->targetnodes = heurdata->minnodes;

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolGins)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* it is likely that the decomposition information changes between runs, we recreate the decomposition horizon */
   decompHorizonFree(scip, &heurdata->decomphorizon);
   assert(heurdata->decomphorizon == NULL);

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

   /* free some data structures that must be reset for a new problem */
   freeTabooList(scip, &heurdata->taboolist);
   SCIPfreeRandom(scip, &heurdata->randnumgen);

   heurdata->taboolist = NULL;
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
   DECOMPHORIZON* decomphorizon;             /* data structure for processing multiple blocks of a decomposition */
   SCIP_DECOMP* decomp;
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

   rollinghorizon = NULL;
   decomp = chooseDecomp(scip);
   if( decomp != NULL && heurdata->usedecomp && heurdata->decomphorizon == NULL )
   {
      SCIP_CALL( decompHorizonCreate(scip, &heurdata->decomphorizon, decomp) );
   }
   decomphorizon = heurdata->decomphorizon;
   /* only create a rolling horizon data structure if we need it */
   if( decomphorizon == NULL && heurdata->userollinghorizon )
   {
      SCIP_CALL( rollingHorizonCreate(scip, &rollinghorizon) );
   }

   do
   {
      SCIP_SOL* oldincumbent;
      SCIP_SOL* newincumbent;

      /* create a new problem, by fixing all variables except for a small neighborhood */
      SCIP_CALL( determineVariableFixings(scip, fixedvars, fixedvals, &nfixedvars, heurdata, rollinghorizon, decomphorizon, &success) );

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

      SCIPdebugMsg(scip, "GINS subscip stats: %.2f sec., %" SCIP_LONGINT_FORMAT " nodes, status:%d\n",
         SCIPgetSolvingTime(subscip), SCIPgetNTotalNodes(subscip), SCIPgetStatus(subscip));

      /* increase target nodes if a (stall) node limit was reached; this will immediately affect the next run */
      if( SCIPgetStatus(subscip) == SCIP_STATUS_NODELIMIT || SCIPgetStatus(subscip) == SCIP_STATUS_STALLNODELIMIT )
      {
         heurdata->targetnodes = (SCIP_Longint)(1.05 * heurdata->targetnodes) + 5L;

         /* but not too far */
         heurdata->targetnodes = MIN(heurdata->targetnodes, heurdata->maxnodes / 2);

         SCIPdebugMsg(scip, "New target nodes after stalling: %" SCIP_LONGINT_FORMAT "\n", heurdata->targetnodes);
      }

      /* transfer variable statistics from sub-SCIP */
      SCIP_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );

      heurdata->usednodes += SCIPgetNNodes(subscip);

      success = FALSE;
      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      oldincumbent = SCIPgetBestSol(scip);

      SCIP_CALL( SCIPtranslateSubSols(scip, subscip, heur, subvars, &success, NULL) );
      if( success )
         *result = SCIP_FOUNDSOL;

      newincumbent = SCIPgetBestSol(scip);

      /* free subproblem */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );

      /* check if we want to run another rolling horizon iteration */
      runagain = success && (newincumbent != oldincumbent) && heurdata->userollinghorizon;
      if( runagain )
      {
         assert(rollinghorizon != NULL || decomphorizon != NULL);
         SCIP_CALL( determineLimits(scip, heur, &solvelimits, &runagain ) );

         if( rollinghorizon != NULL )
            runagain = runagain && rollingHorizonRunAgain(scip, rollinghorizon, heurdata) && (success);
         else if( decomphorizon != NULL )
            runagain = runagain && decompHorizonRunAgain(scip, decomphorizon);
      }
   }
   while( runagain );

   /* delay the heuristic in case it was not successful */
   if( *result != SCIP_FOUNDSOL )
      updateFailureStatistic(scip, heurdata);

TERMINATE:

   /* only free the rolling horizon data structure if we need to keep it */
   if( heurdata->userollinghorizon )
   {
      rollingHorizonFree(scip, &rollinghorizon);
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
   heurdata->decomphorizon = NULL;

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
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolGins) );

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
         "the reference point to compute the neighborhood potential: (r)oot, (l)ocal lp, or (p)seudo solution",
         &heurdata->potential, TRUE, DEFAULT_POTENTIAL, "lpr", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/userollinghorizon",
         "should the heuristic solve a sequence of sub-MIP's around the first selected variable",
         &heurdata->userollinghorizon, TRUE, DEFAULT_USEROLLINGHORIZON, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/relaxdenseconss",
         "should dense constraints (at least as dense as 1 - minfixingrate) be ignored by connectivity graph?",
         &heurdata->relaxdenseconss, TRUE, DEFAULT_RELAXDENSECONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/rollhorizonlimfac",
         "limiting percentage for variables already used in sub-SCIPs to terminate rolling horizon approach",
         &heurdata->rollhorizonlimfac, TRUE, DEFAULT_ROLLHORIZONLIMFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/overlap",
         "overlap of blocks between runs - 0.0: no overlap, 1.0: shift by only 1 block",
         &heurdata->overlap, TRUE, DEFAULT_OVERLAP, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usedecomp",
         "should user decompositions be considered, if available?",
         &heurdata->usedecomp, TRUE, DEFAULT_USEDECOMP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/usedecomprollhorizon",
         "should user decompositions be considered for initial selection in rolling horizon, if available?",
         &heurdata->usedecomprollhorizon, TRUE, DEFAULT_USEDECOMPROLLHORIZON, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/useselfallback",
         "should random initial variable selection be used if decomposition was not successful?",
         &heurdata->useselfallback, TRUE, DEFAULT_USESELFALLBACK, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/consecutiveblocks",
         "should blocks be treated consecutively (sorted by ascending label?)",
         &heurdata->consecutiveblocks, TRUE, DEFAULT_CONSECUTIVEBLOCKS, NULL, NULL) );

   return SCIP_OKAY;
}
