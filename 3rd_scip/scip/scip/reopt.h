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

/**@file   reopt.h
 * @ingroup INTERNALAPI
 * @brief  data structures and methods for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_REOPT_H__
#define __SCIP_REOPT_H__


#include "scip/def.h"
#include "scip/pub_reopt.h"
#include "scip/type_primal.h"
#include "scip/type_sepastore.h"
#include "scip/type_cutpool.h"
#include "scip/type_retcode.h"
#include "scip/type_misc.h"
#include "scip/type_reopt.h"
#include "scip/struct_reopt.h"
#include "scip/struct_var.h"
#include "scip/struct_history.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates reopt data */
extern
SCIP_RETCODE SCIPreoptCreate(
   SCIP_REOPT**          reopt,              /**< pointer to reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** frees reopt data */
extern
SCIP_RETCODE SCIPreoptFree(
   SCIP_REOPT**          reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/* release all variables and constraints captured during reoptimization */
extern
SCIP_RETCODE SCIPreoptReleaseData(
   SCIP_REOPT*           reopt,              /**< pointer to reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
    );

/** returns the number of constraints added by the reoptimization plug-in */
extern
int SCIPreoptGetNAddedConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** add a solution to the solution tree */
extern
SCIP_RETCODE SCIPreoptAddSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_Bool             bestsol,            /**< is the current solution an optimal solution? */
   SCIP_Bool*            added,              /**< pointer to store the information if the soltion was added */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of variables */
   int                   run                 /**< number of the current run (1,2,...) */
   );

/** add optimal solution */
extern
SCIP_RETCODE SCIPreoptAddOptSol(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SOL*             sol,                /**< solution to add */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PRIMAL*          origprimal,         /**< original primal */
   SCIP_VAR**            vars,               /**< original problem variables */
   int                   nvars               /**< number of original problem variables */
   );

/** add a run */
extern
SCIP_RETCODE SCIPreoptAddRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data sturcture */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars,          /**< number of original variables */
   int                   size                /**< number of expected solutions */
   );

/** get the number of checked solutions during the reoptimization process */
extern
int SCIPreoptGetNCheckedSols(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** update the number of checked solutions during the reoptimization process */
extern
void SCIPreoptAddNCheckedSols(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int                   ncheckedsols        /**< number of updated solutions */
   );

/** get the number of checked solutions during the reoptimization process */
extern
int SCIPreoptGetNImprovingSols(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** update the number of checked solutions during the reoptimization process */
extern
void SCIPreoptAddNImprovingSols(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int                   nimprovingsols      /**< number of improving solutions */
   );

/** returns number of solutions stored in the solution tree of a given run */
extern
int SCIPreoptGetNSolsRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   int                   run                 /**< number of the run (1,2,..) */
   );

/** returns number of all solutions of all runs */
extern
int SCIPreoptGetNSols(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** return the stored solutions of a given run */
extern
SCIP_RETCODE SCIPreoptGetSolsRun(
   SCIP_REOPT*           reopt,              /**< reopt data */
   int                   run,                /**< number of the run (1,2,...) */
   SCIP_SOL**            sols,               /**< array of solutions to fill */
   int                   solssize,           /**< length of the array */
   int*                  nsols               /**< pointer to store the number of added solutions */
   );

/** returns the number of saved solutions overall runs */
extern
int SCIPreoptGetNSavedSols(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** Check if the reoptimization process should be (locally) restarted.
 *
 *  First, we check whether the current node is the root node, e.g., node == NULL. In this case, we do not need to calculate
 *  the similarity again. We trigger a restart if
 *    1. the objective function has changed too much, or
 *    2. the number of stored nodes is exceeded, or
 *    3. the last n optimal solutions were found by heur_reoptsols (in this case, the stored tree was only needed to
 *       prove the optimality and this can probably be done faster by solving from scratch).
 *
 *  If the current node is different to the root node we calculate the local similarity, i.e., exclude all variables
 *  that are already fixed at the given node.
 */
extern
SCIP_RETCODE SCIPreoptCheckRestart(
   SCIP_REOPT*           reopt,              /**< reopt data */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< current node of the branch and bound tree (or NULL) */
   SCIP_VAR**            transvars,          /**< transformed problem variables */
   int                   ntransvars,         /**< number of transformed problem variables */
   SCIP_Bool*            restart             /**< pointer to store if the reoptimization process should be restarted */
   );

/** returns the similarity to the previous objective function */
extern
SCIP_Real SCIPreoptGetSimToPrevious(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** returns the similarity to the first objective functions */
extern
SCIP_Real SCIPreoptGetSimToFirst(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** return the similarity between two of objective functions of two given runs */
extern
SCIP_Real SCIPreoptGetSimilarity(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   run1,               /**< number of the first run */
   int                   run2,               /**< number of the second run */
   SCIP_VAR**            origvars,           /**< original problem variables */
   int                   norigvars           /**< number of original problem variables */
   );

/** returns the best solution of the last run */
extern
SCIP_SOL* SCIPreoptGetLastBestSol(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** returns the node of the reoptimization tree corresponding to the unique @p id */
extern
SCIP_REOPTNODE* SCIPreoptGetReoptnode(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   unsigned int          id                  /**< unique id */
   );

/** returns the coefficent of variable with index @p idx in run @p run */
extern
SCIP_Real SCIPreoptGetOldObjCoef(
   SCIP_REOPT*           reopt,              /**< reopt data */
   int                   run,                /**< number of the run */
   int                   idx                 /**< problem index of variable */
   );

/** return the best solution of a given run
 *
 *  @note the returned solution is part of the original space.
 */
SCIP_SOL* SCIPreoptGetBestSolRun(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   int                   run                 /**< number of the run (1,2,...) */
   );

/** reset solving specific paramters */
SCIP_RETCODE SCIPreoptReset(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** reset marks of stored solutions to not updated */
extern
void SCIPreoptResetSolMarks(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** returns the number of stored nodes */
extern
int SCIPreoptGetNNodes(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** save information that given node is infeasible */
extern
SCIP_RETCODE SCIPreoptAddInfNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** check the reason for cut off a node and if necessary store the node */
extern
SCIP_RETCODE SCIPreoptCheckCutoff(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memery */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_EVENTTYPE        eventtype,          /**< eventtype */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPSOLSTAT        lpsolstat,          /**< solution status of the LP */
   SCIP_Bool             isrootnode,         /**< the node is the root */
   SCIP_Bool             isfocusnode,        /**< the node is the current focus node */
   SCIP_Real             lowerbound,         /**< lower bound of the node */
   int                   effectiverootdepth  /**< effective root depth */
   );

/** store bound change based on dual information */
extern
SCIP_RETCODE SCIPreoptAddDualBndchg(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   SCIP_VAR*             var,                /**< variables */
   SCIP_Real             newval,             /**< new bound */
   SCIP_Real             oldval              /**< old bound */
   );

/** returns the number of bound changes based on dual information */
extern
int SCIPreoptGetNDualBndchgs(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** returns the number of leaf nodes of the subtree induced by @p node (of the whole tree if node == NULL) */
extern
int SCIPreoptGetNLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node                /**< node of the search tree */
   );

/** returns the child nodes of @p node that need to be reoptimized next or NULL if @p node is a leaf */
extern
SCIP_RETCODE SCIPreoptGetChildIDs(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         childs,             /**< array to store the child ids */
   int                   childssize,         /**< size of the childs array */
   int*                  nchilds             /**< pointer to store the number of child nodes */
   );

/** returns all leaves of the subtree induced by @p node */
extern
SCIP_RETCODE SCIPreoptGetLeaves(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node,               /**< node of the search tree */
   unsigned int*         leaves,             /**< array to the the ids */
   int                   leavessize,         /**< size of leaves array */
   int*                  nleaves             /**< pointer to store the number of leaf node */
   );

/** returns the time needed to store the nodes for reoptimization */
extern
SCIP_Real SCIPreoptGetSavingtime(
   SCIP_REOPT*           reopt               /**< reoptimization data */
   );

/** store a global constraint that should be added at the beginning of the next iteration */
extern
SCIP_RETCODE SCIPreoptAddGlbCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_VAR**            vars,               /**< array to store the variables of the constraint */
   SCIP_Real*            vals,               /**< array to store the coefficients of the variables */
   int                   nvars,              /**< pointer to store the size of the constraints */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** add the stored constraints globally to the problem */
extern
SCIP_RETCODE SCIPreoptApplyGlbConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** add the stored cuts to the separation storage */
SCIP_RETCODE SCIPreoptApplyCuts(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_NODE*            node,               /**< current focus node */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_CUTPOOL*         cutpool,            /**< global cutpool */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_Bool             root                /**< bool whether the current node is the root */
   );

/** check if the LP of the given node should be solved or not */
extern
SCIP_Bool SCIPreoptGetSolveLP(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node of the current search tree */
   );

/** reactivate the given @p reoptnode and split them into several nodes if necessary */
extern
SCIP_RETCODE SCIPreoptApply(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branching tree */
   SCIP_LP*              lp,                 /**< current LP */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree to reactivate */
   unsigned int          id,                 /**< id of the node to reactivate */
   SCIP_Real             estimate,           /**< estimate of the child nodes that should be created */
   SCIP_NODE**           childnodes,         /**< array to store the created child nodes */
   int*                  ncreatedchilds,     /**< pointer to store number of created child nodes */
   int*                  naddedconss,        /**< pointer to store number of generated constraints */
   int                   childnodessize,     /**< available size of childnodes array */
   SCIP_Bool*            success             /**< pointer store the result */
   );

/** delete a node stored in the reoptimization tree */
extern
SCIP_RETCODE SCIPreoptDeleteNode(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_SET*             set,                /**< global SCIP settings */
   unsigned int          id,                 /**< id of the node */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** reset the stored information abound bound changes based on dual information */
extern
SCIP_RETCODE SCIPreoptResetDualBndchgs(
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_NODE*            node,               /**< node of the search tree */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** splits the root into several nodes and moves the child nodes of the root to one of the created nodes */
extern
SCIP_RETCODE SCIPreoptSplitRoot(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int*                  ncreatedchilds,     /**< pointer to store the number of created nodes */
   int*                  naddedconss         /**< pointer to store the number added constraints */
   );

/** reset the complete tree and set the given search frontier */
extern
SCIP_RETCODE SCIPreoptApplyCompression(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE**      representatives,    /**< array of representatives */
   int                   nrepresentatives,   /**< number of representatives */
   SCIP_Bool*            success             /**< pointer to store if the method was successful */
   );

/** add all unprocessed nodes to the reoptimization tree */
extern
SCIP_RETCODE SCIPreoptSaveOpenNodes(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< LP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NODE**           leaves,             /**< array of open leave nodes */
   int                   nleaves,            /**< number of open leave nodes */
   SCIP_NODE**           childs,             /**< array of open children nodes */
   int                   nchilds,            /**< number of open leave nodes */
   SCIP_NODE**           siblings,           /**< array of open sibling nodes */
   int                   nsiblings           /**< number of open leave nodes */
   );

/** merges the variable history of the current run with the stored history */
extern
SCIP_RETCODE SCIPreoptMergeVarHistory(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR**            vars,               /**< original problem variables */
   int                   nvars               /**< number of original problem variables */
   );

/** updates the variable history */
extern
SCIP_RETCODE SCIPreoptUpdateVarHistory(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of variables */
   );

/*
 * methods for reoptnode
 */

/** initialize an empty node */
extern
void SCIPreoptnodeInit(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** reset the given reoptimization node */
extern
SCIP_RETCODE SCIPreoptnodeReset(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_REOPTNODE*       reoptnode           /**< reoptimization node */
   );

/** delete the given reoptimization node */
extern
SCIP_RETCODE SCIPreoptnodeDelete(
   SCIP_REOPTNODE**      reoptnode,          /**< pointer of reoptnode */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** add a variable to a given reoptnode */
extern
SCIP_RETCODE SCIPreoptnodeAddBndchg(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             val,                /**< value of the variable */
   SCIP_BOUNDTYPE        boundtype           /**< boundtype of the variable */
   );

/** add a constraint to a given reoptnode */
extern
SCIP_RETCODE SCIPreoptnodeAddCons(
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reopttree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR**            vars,               /**< variables which are part of the constraint */
   SCIP_Real*            bounds,             /**< bounds of the variables */
   SCIP_BOUNDTYPE*       boundtypes,         /**< boundtypes of the varibales (or NULL is the constraint is a cut) */
   SCIP_Real             lhs,                /**< lhs of the constraint */
   SCIP_Real             rhs,                /**< rhs of the constraint */
   int                   nvars,              /**< number of variables */
   REOPT_CONSTYPE        constype,           /**< type of the constraint */
   SCIP_Bool             linear              /**< the given constraint has a linear representation */
   );

/** return the branching path of the given node in the reoptimization tree */
extern
void SCIPreoptnodeGetPath(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_REOPTNODE*       reoptnode,          /**< node of the reoptimization tree */
   SCIP_VAR**            vars,               /**< array for variables */
   SCIP_Real*            vals,               /**< array for values */
   SCIP_BOUNDTYPE*       boundtypes,         /**< array for bound types */
   int                   varssize,           /**< size of arrays vars, vals, and boundtypes */
   int*                  nbndchgs,           /**< pointer to store the number of bound changes */
   int*                  nbndchgsafterdual   /**< pointer to store the number of bound changes applied after
                                              *  the first dual reduction at the given node */
   );

/** add a constraint to the reoptimization data structure */
extern
SCIP_RETCODE SCIPreoptAddCons(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS*            cons                /**< constraint to add */
   );

/** save global lower and upper bounds
 *
 *  @note this method can only called once, i.e., after fishing presolving of the first problem
 */
extern
SCIP_RETCODE SCIPreoptSaveGlobalBounds(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** save active constraints
 *
 *  @note this method can only called once, i.e., after fishing presolving of the first problem
 */
extern
SCIP_RETCODE SCIPreoptSaveActiveConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** installs global lower and upper bounds */
extern
SCIP_RETCODE SCIPreoptInstallBounds(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** reactivate globally valid constraints that were deactivated and necessary to ensure correctness */
extern
SCIP_RETCODE SCIPreoptResetActiveConss(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic SCIP statistics */
   );


/** returns whether a constraint is necessary to ensure correctness and cannot be deleted */
extern
SCIP_Bool SCIPreoptConsCanBeDeleted(
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_CONS*            cons                /**< problem constraint */
   );

#ifdef __cplusplus
}
#endif

#endif
