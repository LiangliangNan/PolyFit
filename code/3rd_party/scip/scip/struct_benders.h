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

/**@file   struct_benders.h
 * @ingroup INTERNALAPI
 * @brief  data structures required for Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BENDERS_H__
#define __SCIP_STRUCT_BENDERS_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"

#ifdef __cplusplus
extern "C" {
#endif

struct SCIP_BenderscutCut
{
   SCIP_VAR**            vars;               /**< the variables forming the cut */
   SCIP_Real*            vals;               /**< the coefficients of the variables in the cut */
   SCIP_Real             lhs;                /**< the left hand side of the cut */
   SCIP_Real             rhs;                /**< the right hand side of the cut */
   int                   nvars;              /**< the number of variables in the cut */
};
typedef struct SCIP_BenderscutCut SCIP_BENDERSCUTCUT;

/** Benders' decomposition data */
struct SCIP_Benders
{
   char*                 name;               /**< name of Benders' decomposition */
   char*                 desc;               /**< description of Benders' decomposition */
   SCIP_DECL_BENDERSCOPY ((*benderscopy));   /**< copy method of Benders' decomposition or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BENDERSFREE ((*bendersfree));   /**< destructor of Benders' decomposition */
   SCIP_DECL_BENDERSINIT ((*bendersinit));   /**< initialize Benders' decomposition */
   SCIP_DECL_BENDERSEXIT ((*bendersexit));   /**< deinitialize Benders' decomposition */
   SCIP_DECL_BENDERSINITPRE((*bendersinitpre));/**< presolving initialization method for Benders' decomposition */
   SCIP_DECL_BENDERSEXITPRE((*bendersexitpre));/**< presolving deinitialization method for Benders' decomposition */
   SCIP_DECL_BENDERSINITSOL((*bendersinitsol));/**< solving process initialization method of Benders' decomposition */
   SCIP_DECL_BENDERSEXITSOL((*bendersexitsol));/**< solving process deinitialization method of Benders' decomposition */
   SCIP_DECL_BENDERSGETVAR((*bendersgetvar)); /**< returns the corresponding variable from the master or subproblem */
   SCIP_DECL_BENDERSPRESUBSOLVE((*benderspresubsolve));/**< called prior to the subproblem solving loop */
   SCIP_DECL_BENDERSCREATESUB((*benderscreatesub));/**< creates the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUBCONVEX((*benderssolvesubconvex));/**< the solving method for convex Benders' decomposition subproblems */
   SCIP_DECL_BENDERSSOLVESUB((*benderssolvesub));/**< the solving method for the Benders' decomposition subproblems */
   SCIP_DECL_BENDERSPOSTSOLVE((*benderspostsolve));/**< called after the subproblems are solved. */
   SCIP_DECL_BENDERSFREESUB((*bendersfreesub));/**< the freeing method for the Benders' decomposition subproblems */
   SCIP_DECL_SORTPTRCOMP((*benderssubcomp)); /**< a comparator for defining the solving order of the subproblems */
   SCIP_BENDERSDATA*     bendersdata;        /**< Benders' decomposition local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this Benders' decomposition for the next stages */
   SCIP_CLOCK*           bendersclock;       /**< Benders' decomposition execution time */
   int                   priority;           /**< priority of the Benders' decomposition */
   int                   ncalls;             /**< number of times, this Benders' decomposition was called */
   int                   ncutsfound;         /**< number of cuts found by the Benders' decomposition */
   int                   ntransferred;       /**< number of cuts transferred from sub SCIP to the master SCIP */
   SCIP_Bool             active;             /**< is the Benders' decomposition active? */
   SCIP_Bool             initialized;        /**< is Benders' decomposition initialized? */
   SCIP_Bool             cutlp;              /**< should Benders' cuts be generated for LP solutions? */
   SCIP_Bool             cutpseudo;          /**< should Benders' cuts be generated for pseudo solutions? */
   SCIP_Bool             cutrelax;           /**< should Benders' cuts be generated for relaxation solutions? */
   SCIP_Bool             shareauxvars;       /**< should this Benders' share the highest priority Benders' auxiliary vars */

   /* additional Benders' decomposition parameters */
   SCIP_Bool             transfercuts;       /**< should Benders' cuts generated in LNS heuristics be transferred to the main SCIP instance? */
   SCIP_Bool             lnscheck;           /**< should Benders' decomposition be used in LNS heuristics? */
   int                   lnsmaxdepth;        /**< maximum depth at which the LNS check is performed */
   int                   lnsmaxcalls;        /**< maximum number of Benders' decomposition call in LNS heuristics */
   int                   lnsmaxcallsroot;    /**< maximum number of root node Benders' decomposition call in LNS heuristics */
   SCIP_Bool             cutsasconss;        /**< should the transferred cuts be added as constraints? */
   SCIP_Real             subprobfrac;        /**< fraction of subproblems that are solved in each iteration */
   SCIP_Bool             updateauxvarbound;  /**< should the auxiliary variable lower bound be updated by solving the subproblem? */
   SCIP_Bool             auxvarsimplint;     /**< if subproblem objective is integer, then set the auxiliary variables as implint */
   SCIP_Bool             cutcheck;           /**< should cuts be generated while checking solutions? */
   SCIP_Bool             threadsafe;         /**< has the copy been created requiring thread safety */
   SCIP_Real             solutiontol;        /**< storing the tolerance for optimality in Benders' decomposition */
   int                   numthreads;         /**< the number of threads to use when solving the subproblem */
   SCIP_Bool             execfeasphase;      /**< should a feasibility phase be executed during the root node, i.e.
                                                  adding slack variables to constraints to ensure feasibility */
   SCIP_Real             slackvarcoef;       /**< the initial objective coefficient of the slack variables in the subproblem */
   SCIP_Real             maxslackvarcoef;    /**< the maximal objective coefficient of the slack variables in the subproblem */
   SCIP_Bool             checkconsconvexity; /**< should the constraints of the subproblems be checked for convexity? */

   /* information for heuristics */
   SCIP*                 sourcescip;         /**< the source scip from when the Benders' was copied */
   SCIP_Bool             iscopy;             /**< is the Benders' decomposition struct a copy */
   SCIP_HASHMAP*         mastervarsmap;      /**< hash map for the master variables from the subscip to the master */

   /* the subproblem information */
   SCIP**                subproblems;        /**< the Benders' decomposition subproblems */
   SCIP_VAR**            auxiliaryvars;      /**< the auxiliary variables for the Benders' optimality cuts */
   SCIP_PQUEUE*          subprobqueue;       /**< the priority queue for the subproblems */
   SCIP_SUBPROBLEMSOLVESTAT** solvestat;     /**< storing the solving statistics of all the subproblems */
   SCIP_Real*            subprobobjval;      /**< the objective value of the subproblem in the current iteration */
   SCIP_Real*            bestsubprobobjval;  /**< the best objective value of the subproblem */
   SCIP_Real*            subproblowerbound;  /**< a lower bound on the subproblem - used for the integer cuts */
   int                   naddedsubprobs;     /**< subproblems added to the Benders' decomposition data */
   int                   nsubproblems;       /**< number of subproblems */
   SCIP_BENDERSSUBTYPE*  subprobtype;        /**< the convexity type of the subproblem */
   SCIP_Bool*            subprobisconvex;    /**< is the subproblem convex? This implies that the dual sol can be used for cuts */
   SCIP_Bool*            subprobisnonlinear; /**< does the subproblem contain non-linear constraints */
   int                   nconvexsubprobs;    /**< the number of subproblems that are convex */
   int                   nnonlinearsubprobs; /**< the number of subproblems that are non-linear */
   SCIP_Bool             subprobscreated;    /**< have the subproblems been created for this Benders' decomposition.
                                                  This flag is used when retransforming the problem.*/
   SCIP_Bool*            mastervarscont;     /**< flag to indicate that the master problem variable have been converted
                                               to continuous variables. */
   SCIP_Bool*            subprobsetup;       /**< flag to indicate whether the subproblem has been set up. */
   SCIP_Bool*            indepsubprob;       /**< flag to indicate if a subproblem is independent of the master prob */
   SCIP_Bool*            subprobenabled;     /**< flag to indicate whether the subproblem is enabled */
   int                   nactivesubprobs;    /**< the number of active subproblems */
   SCIP_Bool             freesubprobs;       /**< do the subproblems need to be freed by the Benders' decomposition core? */
   SCIP_Bool             masterisnonlinear;  /**< flag to indicate whether the master problem contains non-linear constraints */

   /* cut strengthening details */
   SCIP_SOL*             corepoint;          /**< the point that is separated for stabilisation */
   SCIP_SOL*             initcorepoint;      /**< the point that was used to initialise the core point */
   SCIP_Real             convexmult;         /**< the multiplier for the convex comb of the LP and sepa point */
   SCIP_Real             perturbeps;         /**< epsilon value to perturb the LP solution */
   int                   noimprovecount;     /**< count of the iterations without improvement */
   int                   noimprovelimit;     /**< limit used to change behaviour of stabilitation */
   SCIP_NODE*            prevnode;           /**< the previous node where the cut strengthening was performed */
   SCIP_Longint          prevnlpiter;        /**< number of LP iters at the previous call of the cut strengthening */
   SCIP_Real             prevlowerbound;     /**< the lowerbound from the previous LP enforcement iteration */
   SCIP_Bool             strengthenenabled;  /**< is the core point cut strengthening enabled */
   char                  strengthenintpoint; /**< where should the strengthening interior point be sourced from ('l'p relaxation, 'f'irst solution, 'i'ncumbent solution, 'r'elative interior point, vector of 'o'nes, vector of 'z'eros)  */
   SCIP_Bool             strengthenround;    /**< flag to indicate whether a cut strengthening round is being performed */
   int                   nstrengthencuts;    /**< the number of strengthened cuts found */
   int                   nstrengthencalls;   /**< the number of calls to the strengthening round */
   int                   nstrengthenfails;   /**< the number of calls to the strengthening round that fail to find cuts */

   /* solving process information */
   int                   npseudosols;        /**< the number of pseudo solutions checked since the last generated cut */
   SCIP_Bool             feasibilityphase;   /**< is the Benders' decomposition in a feasibility phase, i.e. using slack variables */

   /* Bender's cut information */
   SCIP_BENDERSCUT**     benderscuts;        /**< the available Benders' cut algorithms */
   int                   nbenderscuts;       /**< the number of Benders' cut algorithms */
   int                   benderscutssize;    /**< the size of the Benders' cuts algorithms array */
   SCIP_Bool             benderscutssorted;  /**< are the Benders' cuts algorithms sorted by priority */
   SCIP_Bool             benderscutsnamessorted;/**< are the Benders' cuts algorithms sorted by name */

   /* cut storage information */
   SCIP_BENDERSCUTCUT**  storedcuts;         /**< array to store the data required to form a cut/constraint */
   int                   storedcutssize;     /**< the size of the added cuts array */
   int                   nstoredcuts;        /**< the number of the added cuts */

};

/** statistics for solving the subproblems. Used for prioritising the solving of the subproblem */
struct SCIP_SubproblemSolveStat
{
   int                   idx;                /**< the index of the subproblem */
   int                   ncalls;             /**< the number of times this subproblems has been solved */
   SCIP_Real             avgiter;            /**< the average number of LP/NLP iterations performed */
};

/** parameters that are set to solve the subproblem. This will be changed from what the user inputs, so they are stored
 *  and reset after the solving loop. */
struct SCIP_SubproblemParams
{
   SCIP_Real limits_memory;
   SCIP_Real limits_time;
   int cons_linear_propfreq;
   int lp_disablecutoff;
   int lp_scaling;
   int prop_maxrounds;
   int prop_maxroundsroot;
   char lp_initalg;
   char lp_resolvealg;
   SCIP_Bool conflict_enable;
   SCIP_Bool lp_alwaysgetduals;
   SCIP_Bool misc_catchctrlc;
   SCIP_Bool misc_scaleobj;
};
typedef struct SCIP_SubproblemParams SCIP_SUBPROBPARAMS;

#ifdef __cplusplus
}
#endif

#endif
