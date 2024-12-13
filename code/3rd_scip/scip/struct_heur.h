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

/**@file   struct_heur.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_HEUR_H__
#define __SCIP_STRUCT_HEUR_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_heur.h"

#ifdef __cplusplus
extern "C" {
#endif


struct SCIP_DivesetStats
{
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this dive set */
   SCIP_Longint          nlps;               /**< the number of LPs solved by this dive set */
   SCIP_Longint          totaldepth;         /**< the total depth used in this dive set */
   SCIP_Longint          totalsoldepth;      /**< the sum of depths at which this dive set found solutions */
   SCIP_Longint          totalnnodes;        /**< the total number of probing nodes explored by this dive set */
   SCIP_Longint          totalnbacktracks;   /**< the total number of backtracks during the execution of this dive set */
   SCIP_Longint          nsolsfound;         /**< the total number of solutions found */
   SCIP_Longint          nbestsolsfound;     /**< the total number of best solutions found */
   SCIP_Longint          nconflictsfound;    /**< the total number of added conflicts during the execution of this dive set */
   int                   mindepth;           /**< the minimum depth reached by all executions of the dive set */
   int                   maxdepth;           /**< the maximum depth reached by an execution of the dive set */
   int                   minsoldepth;        /**< the minimum depth at which this dive set found a solution */
   int                   maxsoldepth;        /**< the maximum depth at which this dive set found a solution */
   int                   ncalls;             /**< the total number of calls of this dive set */
   int                   nsolcalls;          /**< number of calls with a leaf solution */
};
typedef struct SCIP_DivesetStats SCIP_DIVESETSTATS;

/** common settings for diving heuristics */
struct SCIP_Diveset
{
   SCIP_HEUR*            heur;               /**< the heuristic to which this dive set belongs */
   char*                 name;               /**< name of dive controller, in case that a heuristic has several */
   SCIP_SOL*             sol;                /**< working solution of this dive set */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_DIVESETSTATS*    divesetstats[3];    /**< statistics for individual contexts */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             lpresolvedomchgquot;/**< percentage of immediate domain changes during probing to trigger LP resolve */
   int                   lpsolvefreq;        /**< LP solve frequency for diving heuristics */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   unsigned int          initialseed;        /**< initial seed for the random number generator */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Bool             onlylpbranchcands;  /**< should only LP branching candidates be considered instead of the slower but
                                              *   more general constraint handler diving variable selection? */
   SCIP_Bool             ispublic;           /**< is this dive set publicly available (ie., can be used by other primal heuristics?) */
   SCIP_DIVETYPE         divetypemask;       /**< bit mask that represents the supported dive types by this dive set */
   SCIP_DECL_DIVESETGETSCORE((*divesetgetscore));  /**< method for candidate score and rounding direction */
   SCIP_DECL_DIVESETAVAILABLE((*divesetavailable)); /**< callback to check availability of dive set at the current stage, or NULL if always available */
};

/** primal heuristics data */
struct SCIP_Heur
{
   SCIP_Longint          ncalls;             /**< number of times, this heuristic was called */
   SCIP_Longint          nsolsfound;         /**< number of feasible primal solutions found so far by this heuristic */
   SCIP_Longint          nbestsolsfound;     /**< number of new best primal CIP solutions found so far by this heuristic */
   char*                 name;               /**< name of primal heuristic */
   char*                 desc;               /**< description of primal heuristic */
   SCIP_DECL_HEURCOPY    ((*heurcopy));      /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_HEURFREE    ((*heurfree));      /**< destructor of primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit));      /**< initialize primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit));      /**< deinitialize primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol));   /**< solving process initialization method of primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol));   /**< solving process deinitialization method of primal heuristic */
   SCIP_DECL_HEUREXEC    ((*heurexec));      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata;           /**< primal heuristics local data */
   SCIP_DIVESET**        divesets;           /**< array of diving controllers of this heuristic */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this heuristic for the next stages */
   SCIP_CLOCK*           heurclock;          /**< heuristic execution time */
   int                   priority;           /**< priority of the primal heuristic */
   int                   freq;               /**< frequency for calling primal heuristic */
   int                   freqofs;            /**< frequency offset for calling primal heuristic */
   int                   maxdepth;           /**< maximal depth level to call heuristic at (-1: no limit) */
   int                   delaypos;           /**< position in the delayed heuristics queue, or -1 if not delayed */
   int                   ndivesets;          /**< number of diving controllers of this heuristic */
   SCIP_HEURTIMING       timingmask;         /**< positions in the node solving loop where heuristic should be executed */
   SCIP_Bool             usessubscip;        /**< does the heuristic use a secondary SCIP instance? */
   SCIP_Bool             initialized;        /**< is primal heuristic initialized? */
   char                  dispchar;           /**< display character of primal heuristic */
};

/** variable graph data structure to determine breadth-first distances between variables
 *
 *  the variable graph internally stores a mapping from the variables to the constraints in which they appear.
 *
 *  @see PublicVariableGraphMethods for available methods
 */
struct SCIP_VGraph
{
   SCIP_CONS***          varconss;           /**< constraints of each variable */
   SCIP_HASHTABLE*       visitedconss;       /**< hash table that keeps a record of visited constraints during breadth-first search */
   int*                  nvarconss;          /**< number of constraints for each variable */
   int*                  varconssize;        /**< size array for every varconss entry */
};

#ifdef __cplusplus
}
#endif

#endif
