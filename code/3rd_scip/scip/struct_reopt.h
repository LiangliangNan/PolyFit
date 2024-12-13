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

/**@file   struct_reopt.h
 * @ingroup INTERNALAPI
 * @brief  data structures for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_REOPT_H__
#define __SCIP_STRUCT_REOPT_H__

#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_cons.h"
#include "scip/type_history.h"
#include "scip/type_lp.h"
#include "scip/type_misc.h"
#include "scip/type_reopt.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** nodes of SCIP_SolTree */
struct SCIP_SolNode
{
   SCIP_SOL*             sol;                /**< the stored solution */
   SCIP_SOLNODE*         father;             /**< pointer to the parent node */
   SCIP_SOLNODE*         child;              /**< pointer to left most child node, i.e., node representing the variable
                                               *  with smallest solution value
                                               */
   SCIP_SOLNODE*         sibling;            /**< pointer to next sibling node */
   SCIP_Real             value;              /**< solution value represented by this node */
   SCIP_Bool             updated;            /**< flag if the solution is already updated
                                              *   w.r.t. the new objective function */
#ifndef NDEBUG
   SCIP_VAR*             var;                /**< variable represented by this node */
#endif
};

/** tree for solution */
struct SCIP_SolTree
{
   SCIP_SOLNODE***       sols;               /**< array of arrays of solutions of the reoptimization runs */
   SCIP_SOLNODE*         root;               /**< root node of the solution tree */
   int*                  solssize;           /**< size of sols[x] arrays */
   int*                  nsols;              /**< number of solutions stored in sols[x] array */
};

/** data for constraints to split nodes during reoptimization */
struct SCIP_ReoptConsData
{
   SCIP_VAR**            vars;               /**< array of variables */
   SCIP_Real*            vals;               /**< array of variable coefficients or bounds */
   SCIP_BOUNDTYPE*       boundtypes;         /**< array of variable bounds */
   SCIP_Real             lhs;                /**< left hand side of the constraint */
   SCIP_Real             rhs;                /**< right hand side of the constraint */
   REOPT_CONSTYPE        constype;           /**< type of the constraint */
   SCIP_Bool             linear;             /**< TRUE, iff the constraint is linear, otherwise the constraint is of
                                              *   type bounddisjunction
                                              */
   int                   varssize;           /**< available size in the arrays */
   int                   nvars;              /**< number of entries in the arrays */
};

/** nodes of SCIP_ReoptTree */
struct SCIP_ReoptNode
{
   SCIP_REOPTCONSDATA**  conss;                   /**< array of constraints added to the node, i.e., logic-or constraints */
   SCIP_VAR**            vars;                    /**< variables along the branching path up to the next stored node */
   SCIP_VAR**            afterdualvars;           /**< variables along the branching path after the first decision based on dual information */
   SCIP_REOPTCONSDATA*   dualredscur;             /**< dual reductions that need to be reconstructed the current round */
   SCIP_REOPTCONSDATA*   dualredsnex;             /**< dual reductions that need to be reconstructed the next round */
   SCIP_BOUNDTYPE*       varboundtypes;           /**< boundtypes along the branching path up to the next stored node */
   SCIP_BOUNDTYPE*       afterdualvarboundtypes;  /**< boundtypes along the branching path after the first dual information */
   SCIP_Real*            varbounds;               /**< bounds along the branching path up to the next stored node */
   SCIP_Real*            afterdualvarbounds;      /**< bounds along the branching path after the first decision based on dual information */
   SCIP_Real             lowerbound;              /**< the last lowerbound of this node in the previous iteration */
   SCIP_Bool             dualreds;                /**< flag whether dual reduction were performed */
   int                   nvars;                   /**< number of branching decisions up to the next stored node */
   int                   varssize;                /**< size of allocated memory */
   int                   nafterdualvars;          /**< number of branching decisions after the first dual information */
   int                   afterdualvarssize;       /**< size of allocated memory */
   int                   nchilds;                 /**< number of child nodes */
   int                   allocchildmem;           /**< allocated memory for child nodes */
   int                   nconss;                  /**< number of added constraints */
   int                   consssize;               /**< allocated memory for constraints */
   unsigned int*         childids;                /**< array of child nodes that need to be reoptimized */

   unsigned int          parentID:29;             /**< id of the stored parent node */
   unsigned int          reopttype:3;             /**< reason for storing the node */
};

/* tree to store the current search tree */
struct SCIP_ReoptTree
{
   SCIP_REOPTNODE**      reoptnodes;              /**< array of SCIP_REOPTNODE */
   SCIP_QUEUE*           openids;                 /**< queue of open positions in the reoptnodes array */
   int                   nreoptnodes;             /**< number of saved nodes */
   int                   nfeasnodes;              /**< number of feasible nodes in the current run */
   int                   ntotalfeasnodes;         /**< number of feasible nodes over all runs */
   int                   ninfnodes;               /**< number of (LP-)infeasible nodes in the current run */
   int                   ntotalinfnodes;          /**< number of (LP-)infeasible nodes over all runs */
   int                   nprunednodes;            /**< number of pruned nodes in the current run */
   int                   ntotalprunednodes;       /**< number of pruned nodes over all runs */
   int                   ncutoffreoptnodes;       /**< number of cut off reoptimized nodes in the current run */
   int                   ntotalcutoffreoptnodes;  /**< number of cut off reoptimized nodes over all runs */
   SCIP_Bool             initialized;             /**< is the data structure initialized? */
   unsigned int          reoptnodessize;          /**< size of allocated memory for the reoptnodes array and the openid queue */
};

/** reoptimization data and solution storage */
struct SCIP_Reopt
{
   SCIP_SOL**            prevbestsols;            /**< list of best solutions of all previous rounds */
   SCIP_Real**           objs;                    /**< list of objective coefficients */
   SCIP_HISTORY***       varhistory;              /**< collected variable history */
   SCIP_REOPTCONSDATA**  glbconss;                /**< global constraints that need to be added at the beginning of the next iteration */
   SCIP_REOPTCONSDATA*   dualreds;                /**< dual reductions that probably need to be reconstructed at this node */
   SCIP_REOPTTREE*       reopttree;               /**< data structure to store the current reoptimization search tree */
   SCIP_SOLTREE*         soltree;                 /**< tree to handle all saved solutions */
   SCIP_RANDNUMGEN*      randnumgen;              /**< random number generator */
   SCIP_CLOCK*           savingtime;              /**< time needed to store the nodes */
   SCIP_CONS**           addedconss;              /**< array of added constraints */
   SCIP_Real             simtolastobj;            /**< similarity to the last objective function */
   SCIP_Real             simtofirstobj;           /**< similarity to the first objective function */
   SCIP_Longint          lastbranched;            /**< number of the last branched node */
   SCIP_Longint          lastseennode;            /**< node number of the last caught event */
   int                   nobjvars;                /**< number of variables in the objective function */
   int                   addedconsssize;          /**< size of addedconss array */
   int                   naddedconss;             /**< number of constraints added */
   SCIP_Bool             objhaschanged;           /**< TRUE iff the objective fucntion has changd */
   SCIP_Bool             consadded;               /**< TRUE iff a constraint was added */
   int                   nactiveconss;            /**< number of active constraints stored in activeconss */
   SCIP_CONS**           activeconss;             /**< storage for active constraints */
   int                   nmaxactiveconss;         /**< maximal number of active constraints stored in activeconss */

   /* hashmaps to track global bound reductions and constraints deletion during presolving */
   SCIP_HASHMAP*         glblb;                   /**< global lower bounds after presolving of the first problem */
   SCIP_HASHMAP*         glbub;                   /**< global upper bounds after presolving of the first problem */
   SCIP_HASHSET*         activeconssset;          /**< set of all active constraints after presolving the first problem */

   /* data structure to track decisions based on dual information */
   SCIP_Longint          currentnode;             /**< number of the current node */
   int                   run;                     /**< number of the current reoptimization run */
   int                   runsize;                 /**< allocated memory for runs */
   int                   firstobj;                /**< first non empty objective function */
   int                   noptsolsbyreoptsol;      /**< number of successive optimal solutions found by heur_reoptsols */
   int                   nglbconss;               /**< number of stored global constraints */
   int                   allocmemglbconss;        /**< allocated memory for global constraints */
   int                   ncheckedsols;            /**< number of updated solutions by reoptsols */
   int                   nimprovingsols;          /**< number of improving solutions found by reoptsols */
   int                   nglbrestarts;            /**< number of global restarts */
   int                   ntotallocrestarts;       /**< number of local restarts over all runs */
   int                   nlocrestarts;            /**< number of local restarts in the current iteration */
   int                   firstrestart;            /**< run with the first global restart or -1 of no restart  */
   int                   lastrestart;             /**< run with the last global restart or -1 if no restart */
};

#ifdef __cplusplus
}
#endif

#endif
