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

/**@file   type_cons.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for constraints and constraint handlers
 * @author Tobias Achterberg
 * @author Stefan Heinz
 *
 *  This file defines the interface for constraint handlers implemented in C.
 *
 *  - \ref CONS "Instructions for implementing a constraint handler"
 *  - \ref CONSHDLRS "List of available constraint handlers"
 *  - \ref scip::ObjConshdlr "C++ wrapper class"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CONS_H__
#define __SCIP_TYPE_CONS_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_scip.h"
#include "scip/type_timing.h"
#include "scip/type_heur.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Conshdlr SCIP_CONSHDLR;       /**< constraint handler for a specific constraint type */
typedef struct SCIP_Cons SCIP_CONS;               /**< constraint data structure */
typedef struct SCIP_ConshdlrData SCIP_CONSHDLRDATA; /**< constraint handler data */
typedef struct SCIP_ConsData SCIP_CONSDATA;       /**< locally defined constraint type specific data */
typedef struct SCIP_ConsSetChg SCIP_CONSSETCHG;   /**< tracks additions and removals of the set of active constraints */
typedef struct SCIP_LinConsStats SCIP_LINCONSSTATS; /**< linear constraint classification statistics used for MIPLIB */

/** linear constraint types recognizable */
enum SCIP_LinConstype
{
   SCIP_LINCONSTYPE_EMPTY         =  0,         /**< linear constraints with no variables */
   SCIP_LINCONSTYPE_FREE          =  1,         /**< linear constraints with no finite side */
   SCIP_LINCONSTYPE_SINGLETON     =  2,         /**< linear constraints with a single variable */
   SCIP_LINCONSTYPE_AGGREGATION   =  3,         /**< linear constraints of the type \f$ ax + by = c\f$ */
   SCIP_LINCONSTYPE_PRECEDENCE    =  4,         /**< linear constraints of the type \f$ a x - a y \leq b\f$ where \f$x\f$ and \f$y\f$ must have the same type */
   SCIP_LINCONSTYPE_VARBOUND      =  5,         /**< linear constraints of the form \f$ ax + by \leq c \, x \in \{0,1\} \f$ */
   SCIP_LINCONSTYPE_SETPARTITION  =  6,         /**< linear constraints of the form \f$ \sum x_i = 1\, x_i \in \{0,1\} \forall i \f$ */
   SCIP_LINCONSTYPE_SETPACKING    =  7,         /**< linear constraints of the form \f$ \sum x_i \leq 1\, x_i \in \{0,1\} \forall i \f$ */
   SCIP_LINCONSTYPE_SETCOVERING   =  8,         /**< linear constraints of the form \f$ \sum x_i \geq 1\, x_i \in \{0,1\} \forall i \f$ */
   SCIP_LINCONSTYPE_CARDINALITY   =  9,         /**< linear constraints of the form \f$ \sum x_i = k\, x_i \in \{0,1\} \forall i, \, k\geq 2 \f$ */
   SCIP_LINCONSTYPE_INVKNAPSACK   = 10,         /**< linear constraints of the form \f$ \sum x_i \leq b\, x_i \in \{0,1\} \forall i, \, b\in \mathbb{n} \geq 2 \f$ */
   SCIP_LINCONSTYPE_EQKNAPSACK    = 11,         /**< linear constraints of the form \f$ \sum a_i x_i = b\, x_i \in \{0,1\} \forall i, \, b\in \mathbb{n} \geq 2 \f$ */
   SCIP_LINCONSTYPE_BINPACKING    = 12,         /**< linear constraints of the form \f$ \sum a_i x_i + a x \leq a\, x, x_i \in \{0,1\} \forall i, \, a\in \mathbb{n} \geq 2 \f$ */
   SCIP_LINCONSTYPE_KNAPSACK      = 13,         /**< linear constraints of the form \f$ \sum a_k x_k \leq b\, x_i \in \{0,1\} \forall i, \, b\in \mathbb{n} \geq 2 \f$ */
   SCIP_LINCONSTYPE_INTKNAPSACK   = 14,         /**< linear constraints of the form \f$ \sum a_k x_k \leq b\, x_i \in \mathbb{Z} \forall i, \, b\in \mathbb{n} \f$ */
   SCIP_LINCONSTYPE_MIXEDBINARY   = 15,         /**< linear constraints of the form \f$ \sum a_k x_k + \sum p_j s_j \leq/= b\, x_i \in \{0,1\} \forall i, s_j \in \text{ cont. } \forall j\f$ */
   SCIP_LINCONSTYPE_GENERAL       = 16          /**< general linear constraints with no special structure */
};
typedef enum SCIP_LinConstype SCIP_LINCONSTYPE;

#define SCIP_NLINCONSTYPES ((int)SCIP_LINCONSTYPE_GENERAL+1)

/** copy method for constraint handler plugins (called when SCIP copies plugins)
 *
 *  If the copy process was one to one, the valid pointer can be set to TRUE. Otherwise, this pointer has to be set to
 *  FALSE. If all problem defining objects (constraint handlers and variable pricers) return valid = TRUE for all
 *  their copying calls, SCIP assumes that it is an overall one to one copy of the original instance. In this case any
 *  reductions made in the copied SCIP instance can be transfered to the original SCIP instance. If the valid pointer is
 *  set to TRUE and it was not a one to one copy, it might happen that optimal solutions are cut off.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - valid           : was the copying process valid?
 */
#define SCIP_DECL_CONSHDLRCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_Bool* valid)

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 */
#define SCIP_DECL_CONSFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr)

/** initialization method of constraint handler (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints in transformed problem
 *  - nconss          : number of constraints in transformed problem
 */
#define SCIP_DECL_CONSINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss)

/** deinitialization method of constraint handler (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints in transformed problem
 *  - nconss          : number of constraints in transformed problem
 */
#define SCIP_DECL_CONSEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss)

/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This method is called when the presolving process is about to begin, even if presolving is turned off.
 *  The constraint handler may use this call to initialize its data structures.
 *
 *  Necessary modifications that have to be performed even if presolving is turned off should be done here or in the
 *  presolving deinitialization call (SCIP_DECL_CONSEXITPRE()).
 *
 *  @note Note that the constraint array might contain constraints that were created but not added to the problem.
 *        Constraints that are not added, i.e., for which SCIPconsIsAdded() returns FALSE, cannot be used for problem
 *        reductions.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints in transformed problem
 *  - nconss          : number of constraints in transformed problem
 */
#define SCIP_DECL_CONSINITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss)

/** presolving deinitialization method of constraint handler (called after presolving has been finished)
 *
 *  This method is called after the presolving has been finished, even if presolving is turned off.
 *  The constraint handler may use this call e.g. to clean up or modify its data structures.
 *
 *  Necessary modifications that have to be performed even if presolving is turned off should be done here or in the
 *  presolving initialization call (SCIP_DECL_CONSINITPRE()).
 *
 *  Besides necessary modifications and clean up, no time consuming operations should be performed, especially if the
 *  problem has already been solved.  Use the method SCIPgetStatus(), which in this case returns SCIP_STATUS_OPTIMAL,
 *  SCIP_STATUS_INFEASIBLE, SCIP_STATUS_UNBOUNDED, or SCIP_STATUS_INFORUNBD.
 *
 *  @note Note that the constraint array might contain constraints that were created but not added to the problem.
 *        Constraints that are not added, i.e., for which SCIPconsIsAdded() returns FALSE, cannot be used for problem
 *        reductions.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : final array of constraints in transformed problem
 *  - nconss          : final number of constraints in transformed problem
 */
#define SCIP_DECL_CONSEXITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss)

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The constraint handler may use this call to initialize its branch and bound specific data.
 *
 *  Besides necessary modifications and clean up, no time consuming operations should be performed, especially if the
 *  problem has already been solved.  Use the method SCIPgetStatus(), which in this case returns SCIP_STATUS_OPTIMAL,
 *  SCIP_STATUS_INFEASIBLE, SCIP_STATUS_UNBOUNDED, or SCIP_STATUS_INFORUNBD.
 *
 *  @note Note that the constraint array might contain constraints that were created but not added to the problem.
 *        Constraints that are not added, i.e., for which SCIPconsIsAdded() returns FALSE, cannot be used for problem
 *        reductions.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints of the constraint handler
 *  - nconss          : number of constraints of the constraint handler
 */
#define SCIP_DECL_CONSINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss)

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The constraint handler should use this call to clean up its branch and bound data, in particular to release
 *  all LP rows that he has created or captured.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints of the constraint handler
 *  - nconss          : number of constraints of the constraint handler
 *  - restart         : was this exit solve call triggered by a restart?
 */
#define SCIP_DECL_CONSEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, SCIP_Bool restart)

/** frees specific constraint data
 *
 *  @warning There may exist unprocessed events. For example, a variable's bound may have been already changed, but the
 *           corresponding bound change event was not yet processed.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint belonging to the constraint data
 *  - consdata        : pointer to the constraint data to free
 */
#define SCIP_DECL_CONSDELETE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons, SCIP_CONSDATA** consdata)

/** transforms constraint data into data belonging to the transformed problem
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - sourcecons      : source constraint to transform
 *  - targetcons      : pointer to store created target constraint
 */ 
#define SCIP_DECL_CONSTRANS(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* sourcecons, SCIP_CONS** targetcons)

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved)
 *
 *  Puts the LP relaxations of all "initial" constraints into the LP. The method should put a canonic LP relaxation
 *  of all given constraints to the LP with calls to SCIPaddRow().
 *
 *  @warning It is not guaranteed that the problem is going to be declared infeasible if the infeasible pointer is set
 *           to TRUE. Therefore, it is recommended that users do not end this method prematurely when an infeasiblity
 *           is detected.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *
 *  output:
 *  - infeasible      : pointer to store whether an infeasibility was detected while building the LP
 */
#define SCIP_DECL_CONSINITLP(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, SCIP_Bool* infeasible)

/** separation method of constraint handler for LP solution
 *
 *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
 *  which means that a valid LP solution exists.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - result          : pointer to store the result of the separation call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_NEWROUND   : a cutting plane was generated and a new separation round should immediately start
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
#define SCIP_DECL_CONSSEPALP(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, \
      int nconss, int nusefulconss, SCIP_RESULT* result)

/** separation method of constraint handler for arbitrary primal solution
 *
 *  Separates all constraints of the constraint handler. The method is called outside the LP solution loop (e.g., by
 *  a relaxator or a primal heuristic), which means that there is no valid LP solution.
 *  Instead, the method should produce cuts that separate the given solution.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - sol             : primal solution that should be separated
 *  - result          : pointer to store the result of the separation call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_NEWROUND   : a cutting plane was generated and a new separation round should immediately start
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
#define SCIP_DECL_CONSSEPASOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, \
      int nconss, int nusefulconss, SCIP_SOL* sol, SCIP_RESULT* result)

/** constraint enforcing method of constraint handler for LP solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was solved.
 *  The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or separating the solution with a valid
 *  cutting plane.
 *
 *  The enforcing methods of the active constraint handlers are called in decreasing order of their enforcing
 *  priorities until the first constraint handler returned with the value SCIP_CUTOFF, SCIP_SEPARATED,
 *  SCIP_REDUCEDDOM, SCIP_CONSADDED, or SCIP_BRANCHED.
 *  The integrality constraint handler has an enforcing priority of zero. A constraint handler which can
 *  (or wants) to enforce its constraints only for integral solutions should have a negative enforcing priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to incorporate its own branching strategy even on non-integral
 *  solutions must have an enforcing priority greater than zero (e.g. the SOS-constraint incorporates
 *  SOS-branching on non-integral solutions).
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
 *  be enforced, if no violation was found in the useful constraints.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - solinfeasible   : was the solution already declared infeasible by a constraint handler?
 *  - result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define SCIP_DECL_CONSENFOLP(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nusefulconss, \
      SCIP_Bool solinfeasible, SCIP_RESULT* result)

/** constraint enforcing method of constraint handler for relaxation solutions
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sol             : relaxation solution
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - solinfeasible   : was the solution already declared infeasible by a constraint handler?
 *  - result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the LP
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is true)
 */
#define SCIP_DECL_CONSENFORELAX(x) SCIP_RETCODE x (SCIP* scip, SCIP_SOL* sol, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nusefulconss, \
      SCIP_Bool solinfeasible, SCIP_RESULT* result)

/** constraint enforcing method of constraint handler for pseudo solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was not solved.
 *  The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or adding an additional constraint.
 *  Separation is not possible, since the LP is not processed at the current node. All LP informations like
 *  LP solution, slack values, or reduced costs are invalid and must not be accessed.
 *
 *  Like in the enforcing method for LP solutions, the enforcing methods of the active constraint handlers are
 *  called in decreasing order of their enforcing priorities until the first constraint handler returned with
 *  the value SCIP_CUTOFF, SCIP_REDUCEDDOM, SCIP_CONSADDED, SCIP_BRANCHED, or SCIP_SOLVELP.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
 *  be enforced, if no violation was found in the useful constraints.
 *
 *  If the pseudo solution's objective value is lower than the lower bound of the node, it cannot be feasible
 *  and the enforcing method may skip it's check and set *result to SCIP_DIDNOTRUN. However, it can also process
 *  its constraints and return any other possible result code.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - solinfeasible   : was the solution already declared infeasible by a constraint handler?
 *  - objinfeasible   : is the solution infeasible anyway due to violating lower objective bound?
 *  - result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the LP
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is true)
 */
#define SCIP_DECL_CONSENFOPS(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nusefulconss, \
      SCIP_Bool solinfeasible, SCIP_Bool objinfeasible, SCIP_RESULT* result)

/** feasibility check method of constraint handler for integral solutions
 *
 *  The given solution has to be checked for feasibility.
 *  
 *  The check methods of the active constraint handlers are called in decreasing order of their check
 *  priorities until the first constraint handler returned with the result SCIP_INFEASIBLE.
 *  The integrality constraint handler has a check priority of zero. A constraint handler which can
 *  (or wants) to check its constraints only for integral solutions should have a negative check priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to check feasibility even on non-integral solutions must have a
 *  check priority greater than zero (e.g. if the check is much faster than testing all variables for
 *  integrality).
 *
 *  In some cases, integrality conditions or rows of the current LP don't have to be checked, because their
 *  feasibility is already checked or implicitly given. In these cases, 'checkintegrality' or
 *  'checklprows' is FALSE.
 *
 *  If the solution is not NULL, SCIP should also be informed about the constraint violation with a call to
 *  SCIPupdateSolConsViolation() and additionally SCIPupdateSolLPRowViolation() for every row of the constraint's current
 *  representation in the LP relaxation, if any such rows exist.
 *  As a convenience method, SCIPupdateSolLPConsViolation() can be used if the constraint
 *  is represented completely by a set of LP rows, meaning that the current constraint violation is equal to the maximum
 *  of the contraint violations of the corresponding LP rows.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - sol             : the solution to check feasibility for
 *  - checkintegrality: Has integrality to be checked?
 *  - checklprows     : Do constraints represented by rows in the current LP have to be checked?
 *  - printreason     : Should the reason for the violation be printed?
 *  - completely      : Should all violations be checked?
 *  - result          : pointer to store the result of the feasibility checking call
 *
 *  possible return values for *result:
 *  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define SCIP_DECL_CONSCHECK(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, SCIP_SOL* sol, \
      SCIP_Bool checkintegrality, SCIP_Bool checklprows, SCIP_Bool printreason, SCIP_Bool completely, SCIP_RESULT* result)

/** domain propagation method of constraint handler
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  @note if the constraint handler uses dual information in propagation it is nesassary to check via calling
 *        SCIPallowDualReds and SCIPallowObjProp if dual reductions and propgation with the current cutoff bound, resp.,
 *        are allowed.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - nmarkedconss    : number of constraints which are marked to be definitely propagated
 *  - proptiming      : current point in the node solving loop
 *  - result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched but did not find any domain reductions
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
 *  - SCIP_DELAYNODE  : the current node should be postponed (return value only valid for BEFORELP propagation)
 */
#define SCIP_DECL_CONSPROP(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nusefulconss, \
      int nmarkedconss, SCIP_PROPTIMING proptiming, SCIP_RESULT* result)

/** presolving method of constraint handler
 *
 *  The presolver should go through the variables and constraints and tighten the domains or
 *  constraints. Each tightening should increase the given total number of changes.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nrounds         : number of presolving rounds already done
 *  - presoltiming    : current presolving timing
 *  - nnewfixedvars   : number of variables fixed since the last call to the presolving method
 *  - nnewaggrvars    : number of variables aggregated since the last call to the presolving method
 *  - nnewchgvartypes : number of variable type changes since the last call to the presolving method
 *  - nnewchgbds      : number of variable bounds tightened since the last call to the presolving method
 *  - nnewholes       : number of domain holes added since the last call to the presolving method
 *  - nnewdelconss    : number of deleted constraints since the last call to the presolving method
 *  - nnewaddconss    : number of added constraints since the last call to the presolving method
 *  - nnewupgdconss   : number of upgraded constraints since the last call to the presolving method
 *  - nnewchgcoefs    : number of changed coefficients since the last call to the presolving method
 *  - nnewchgsides    : number of changed left or right hand sides since the last call to the presolving method
 *
 *  @note the counters state the changes since the last call including the changes of this presolving method during its
 *        call
 *
 *  @note if the constraint handler performs dual presolving it is nesassary to check via calling SCIPallowDualReds
 *        if dual reductions are allowed.
 *
 *  input/output:
 *  - nfixedvars      : pointer to count total number of variables fixed of all presolvers
 *  - naggrvars       : pointer to count total number of variables aggregated of all presolvers
 *  - nchgvartypes    : pointer to count total number of variable type changes of all presolvers
 *  - nchgbds         : pointer to count total number of variable bounds tightened of all presolvers
 *  - naddholes       : pointer to count total number of domain holes added of all presolvers
 *  - ndelconss       : pointer to count total number of deleted constraints of all presolvers
 *  - naddconss       : pointer to count total number of added constraints of all presolvers
 *  - nupgdconss      : pointer to count total number of upgraded constraints of all presolvers
 *  - nchgcoefs       : pointer to count total number of changed coefficients of all presolvers
 *  - nchgsides       : pointer to count total number of changed left/right hand sides of all presolvers
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_SUCCESS    : the presolving method found a reduction
 *  - SCIP_DIDNOTFIND : the presolving method searched, but did not find a presolving change
 *  - SCIP_DIDNOTRUN  : the presolving method was skipped
 *  - SCIP_DELAYED    : the presolving method was skipped, but should be called again
 */
#define SCIP_DECL_CONSPRESOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss, int nrounds, \
      SCIP_PRESOLTIMING presoltiming, int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
      int nnewdelconss, int nnewaddconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides, \
      int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes, \
      int* ndelconss, int* naddconss, int* nupgdconss, int* nchgcoefs, int* nchgsides, SCIP_RESULT* result)

/** propagation conflict resolving method of constraint handler
 *
 *  This method is called during conflict analysis. If the constraint handler wants to support conflict analysis,
 *  it should call SCIPinferVarLbCons() or SCIPinferVarUbCons() in domain propagation instead of SCIPchgVarLb() or
 *  SCIPchgVarUb() in order to deduce bound changes on variables.
 *  In the SCIPinferVarLbCons() and SCIPinferVarUbCons() calls, the handler provides the constraint, that deduced the
 *  variable's bound change, and an integer value "inferinfo" that can be arbitrarily chosen.
 *  The propagation conflict resolving method can then be implemented, to provide a "reason" for the bound
 *  changes, i.e., the bounds of variables at the time of the propagation, that forced the constraint to set the
 *  conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
 *  rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
 *  by calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(),
 *  SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), and/or SCIPaddConflictBinvar() in the propagation conflict
 *  resolving method.
 *
 *  For example, the logicor constraint c = "x or y or z" fixes variable z to TRUE (i.e. changes the lower bound of z
 *  to 1.0), if both, x and y, are assigned to FALSE (i.e. if the upper bounds of these variables are 0.0). It uses
 *  SCIPinferVarLbCons(scip, z, 1.0, c, 0) to apply this assignment (an inference information tag is not needed by the
 *  constraint handler and is set to 0).
 *  In the conflict analysis, the constraint handler may be asked to resolve the lower bound change on z with
 *  constraint c, that was applied at a time given by a bound change index "bdchgidx".
 *  With a call to SCIPgetVarLbAtIndex(scip, z, bdchgidx, TRUE), the handler can find out, that the lower bound of
 *  variable z was set to 1.0 at the given point of time, and should call SCIPaddConflictUb(scip, x, bdchgidx) and
 *  SCIPaddConflictUb(scip, y, bdchgidx) to tell SCIP, that the upper bounds of x and y at this point of time were
 *  the reason for the deduction of the lower bound of z.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that deduced the bound change of the conflict variable
 *  - infervar        : the conflict variable whose bound change has to be resolved
 *  - inferinfo       : the user information passed to the corresponding SCIPinferVarLbCons() or SCIPinferVarUbCons() call
 *  - boundtype       : the type of the changed bound (lower or upper bound)
 *  - bdchgidx        : the index of the bound change, representing the point of time where the change took place
 *  - relaxedbd       : the relaxed bound which is sufficient to be explained
 *
 *  output:
 *  - result          : pointer to store the result of the propagation conflict resolving call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the conflicting bound change has been successfully resolved by adding all reason bounds
 *  - SCIP_DIDNOTFIND : the conflicting bound change could not be resolved and has to be put into the conflict set
 *
 *  @note it is sufficient to explain/resolve the relaxed bound
 */
#define SCIP_DECL_CONSRESPROP(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons, \
      SCIP_VAR* infervar, int inferinfo, SCIP_BOUNDTYPE boundtype, SCIP_BDCHGIDX* bdchgidx, SCIP_Real relaxedbd, \
      SCIP_RESULT* result)

/** variable rounding lock method of constraint handler
 *
 *  This method is called, after a constraint is added or removed from the transformed problem.
 *  It should update the rounding locks of all associated variables with calls to SCIPaddVarLocks(),
 *  depending on the way, the variable is involved in the constraint:
 *  - If the constraint may get violated by decreasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the
 *    (positive) constraint infeasible and rounding up is potentially rendering the negation of the constraint
 *    infeasible.
 *  - If the constraint may get violated by increasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlocksneg, nlockspos), saying that rounding up is potentially rendering the
 *    constraint's negation infeasible and rounding up is potentially rendering the constraint itself
 *    infeasible.
 *  - If the constraint may get violated by changing the variable in any direction, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg).
 *
 *  Consider the linear constraint "3x -5y +2z <= 7" as an example. The variable rounding lock method of the
 *  linear constraint handler should call SCIPaddVarLocks(scip, x, nlocksneg, nlockspos), 
 *  SCIPaddVarLocks(scip, y, nlockspos, nlocksneg) and SCIPaddVarLocks(scip, z, nlocksneg, nlockspos) to tell SCIP,
 *  that rounding up of x and z and rounding down of y can destroy the feasibility of the constraint, while rounding
 *  down of x and z and rounding up of y can destroy the feasibility of the constraint's negation "3x -5y +2z > 7".
 *  A linear constraint "2 <= 3x -5y +2z <= 7" should call
 *  SCIPaddVarLocks(scip, ..., nlockspos + nlocksneg, nlockspos + nlocksneg) on all variables, since rounding in both
 *  directions of each variable can destroy both the feasibility of the constraint and it's negation
 *  "3x -5y +2z < 2  or  3x -5y +2z > 7".
 *
 *  If the constraint itself contains other constraints as sub constraints (e.g. the "or" constraint concatenation
 *  "c(x) or d(x)"), the rounding lock methods of these constraints should be called in a proper way.
 *  - If the constraint may get violated by the violation of the sub constraint c, it should call
 *    SCIPaddConsLocks(scip, c, nlockspos, nlocksneg), saying that infeasibility of c may lead to infeasibility of
 *    the (positive) constraint, and infeasibility of c's negation (i.e. feasibility of c) may lead to infeasibility
 *    of the constraint's negation (i.e. feasibility of the constraint).
 *  - If the constraint may get violated by the feasibility of the sub constraint c, it should call
 *    SCIPaddConsLocks(scip, c, nlocksneg, nlockspos), saying that infeasibility of c may lead to infeasibility of
 *    the constraint's negation (i.e. feasibility of the constraint), and infeasibility of c's negation (i.e. feasibility
 *    of c) may lead to infeasibility of the (positive) constraint.
 *  - If the constraint may get violated by any change in the feasibility of the sub constraint c, it should call
 *    SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg).
 *
 *  Consider the or concatenation "c(x) or d(x)". The variable rounding lock method of the or constraint handler
 *  should call SCIPaddConsLocks(scip, c, nlockspos, nlocksneg) and SCIPaddConsLocks(scip, d, nlockspos, nlocksneg)
 *  to tell SCIP, that infeasibility of c and d can lead to infeasibility of "c(x) or d(x)".
 *
 *  As a second example, consider the equivalence constraint "y <-> c(x)" with variable y and constraint c. The
 *  constraint demands, that y == 1 if and only if c(x) is satisfied. The variable lock method of the corresponding
 *  constraint handler should call SCIPaddVarLocks(scip, y, nlockspos + nlocksneg, nlockspos + nlocksneg) and
 *  SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg), because any modification to the
 *  value of y or to the feasibility of c can alter the feasibility of the equivalence constraint.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that should lock rounding of its variables, or NULL if the constraint handler
 *                      does not need constraints
 *  - nlockspos       : number of times, the roundings should be locked for the constraint (may be negative)
 *  - nlocksneg       : number of times, the roundings should be locked for the constraint's negation (may be negative)
 */
#define SCIP_DECL_CONSLOCK(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons, int nlockspos, int nlocksneg)

/** constraint activation notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called after a constraint of the constraint handler was activated. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that has been activated
 */
#define SCIP_DECL_CONSACTIVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons)

/** constraint deactivation notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called before a constraint of the constraint handler is deactivated. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that will be deactivated
 */
#define SCIP_DECL_CONSDEACTIVE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons)

/** constraint enabling notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called after a constraint of the constraint handler was enabled. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that has been enabled
 */
#define SCIP_DECL_CONSENABLE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons)

/** constraint disabling notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called before a constraint of the constraint handler is disabled. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that will be disabled
 */
#define SCIP_DECL_CONSDISABLE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons)

/** variable deletion method of constraint handler
 *
 *  This method is optinal and only of interest if you are using SCIP as a branch-and-price framework. That means, you
 *  are generating new variables during the search. If you are not doing that just define the function pointer to be
 *  NULL.
 *
 *  If this method gets implemented you should iterate over all constraints of the constraint handler and delete all
 *  variables that were marked for deletion by SCIPdelVar().
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints in transformed problem
 *  - nconss          : number of constraints in transformed problem
 */
#define SCIP_DECL_CONSDELVARS(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** conss, int nconss)

/** constraint display method of constraint handler
 *
 *  The constraint handler can store a representation of the constraint into the given text file. Use the method
 *  SCIPinfoMessage() to push a string into the file stream.
 *
 *  @note There are several methods which help to display variables. These are SCIPwriteVarName(), SCIPwriteVarsList(),
 *        SCIPwriteVarsLinearsum(), and SCIPwriteVarsPolynomial().
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that should be displayed
 *  - file            : the text file to store the information into
 */
#define SCIP_DECL_CONSPRINT(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons, FILE* file)

/** constraint copying method of constraint handler
 *
 *  The constraint handler can provide a copy method which copies a constraint from one SCIP data structure into an other
 *  SCIP data structure. If a copy of a constraint is created, the constraint has to be captured. (The capture is usually
 *  already done due to the creation of the constraint).
 *
 *  If the copy process was one to one, the valid pointer can be set to TRUE. Otherwise, you have to set this pointer to
 *  FALSE. In case all problem defining objects (constraint handlers and variable pricers) return a TRUE valid for all
 *  their copying calls, SCIP assumes that it is a overall one to one copy of the original instance. In this case any
 *  reductions made in the copied SCIP instance can be transfered to the original SCIP instance. If the valid pointer is
 *  set to TRUE and it was not a one to one copy, it might happen that optimal solutions are cut off.
 *
 *  To get a copy of a variable in the target SCIP you should use the function SCIPgetVarCopy().
 *
 *  input:
 *  - scip            : target SCIP data structure
 *  - cons            : pointer to store the created target constraint
 *  - name            : name of constraint, or NULL if the name of the source constraint should be used
 *  - sourcescip      : source SCIP data structure
 *  - sourceconshdlr  : source constraint handler of the source SCIP
 *  - sourcecons      : source constraint of the source SCIP
 *  - varmap          : a SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP
 *  - consmap         : a SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP
 *  - initial         : should the LP relaxation of constraint be in the initial LP?
 *  - separate        : should the constraint be separated during LP processing?
 *  - enforce         : should the constraint be enforced during node processing?
 *  - check           : should the constraint be checked for feasibility?
 *  - propagate       : should the constraint be propagated during node processing?
 *  - local           : is constraint only valid locally?
 *  - modifiable      : is constraint modifiable (subject to column generation)?
 *  - dynamic         : is constraint subject to aging?
 *  - removable       : should the relaxation be removed from the LP due to aging or cleanup?
 *  - stickingatnode  : should the constraint always be kept at the node where it was added, even
 *                      if it may be moved to a more global node?
 *  - global          : should a global or a local copy be created?
 *
 *  output:
 *  - valid           : pointer to store whether the copying was valid or not 
 */
#define SCIP_DECL_CONSCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS** cons, const char* name, \
      SCIP* sourcescip, SCIP_CONSHDLR* sourceconshdlr, SCIP_CONS* sourcecons, SCIP_HASHMAP* varmap, SCIP_HASHMAP* consmap, \
      SCIP_Bool initial, SCIP_Bool separate, SCIP_Bool enforce, SCIP_Bool check, SCIP_Bool propagate, \
      SCIP_Bool local, SCIP_Bool modifiable, SCIP_Bool dynamic, SCIP_Bool removable, SCIP_Bool stickingatnode, \
      SCIP_Bool global, SCIP_Bool* valid)

/** constraint parsing method of constraint handler
 *
 *  The constraint handler can provide a callback to parse the output created by the display method
 *  (\ref SCIP_DECL_CONSPRINT) and to create a constraint out of it.
 *
 *  @note For parsing there are several methods which are handy. Have a look at: SCIPparseVarName(),
 *        SCIPparseVarsList(), SCIPparseVarsLinearsum(), SCIPparseVarsPolynomial(), SCIPstrToRealValue(), and
 *        SCIPstrCopySection().
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : pointer to store the created constraint
 *  - name            : name of the constraint
 *  - str             : string to parse
 *  - initial         : should the LP relaxation of constraint be in the initial LP?
 *  - separate        : should the constraint be separated during LP processing?
 *  - enforce         : should the constraint be enforced during node processing?
 *  - check           : should the constraint be checked for feasibility?
 *  - propagate       : should the constraint be propagated during node processing?
 *  - local           : is constraint only valid locally?
 *  - modifiable      : is constraint modifiable (subject to column generation)?
 *  - dynamic         : is constraint subject to aging?
 *  - removable       : should the relaxation be removed from the LP due to aging or cleanup?
 *  - stickingatnode  : should the constraint always be kept at the node where it was added, even
 *                      if it may be moved to a more global node?
 *  output:
 *  - success         : pointer to store whether the parsing was successful or not
 */
#define SCIP_DECL_CONSPARSE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS** cons, \
      const char* name, const char* str, \
      SCIP_Bool initial, SCIP_Bool separate, SCIP_Bool enforce, SCIP_Bool check, SCIP_Bool propagate, SCIP_Bool local, \
      SCIP_Bool modifiable, SCIP_Bool dynamic, SCIP_Bool removable, SCIP_Bool stickingatnode, SCIP_Bool* success)

/** constraint method of constraint handler which returns the variables (if possible)
 *
 *  The constraint handler can (this callback is optional) provide this callback to return the variables which are
 *  involved in that particular constraint. If this is possible, the variables should be copyied into the variables
 *  array and the success pointers has to be set to TRUE. Otherwise the success has to be set FALSE or the callback
 *  should not be implemented.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that should return its variable data
 *  - varssize        : available slots in vars array which is needed to check if the array is large enough
 *
 *  output:
 *  - vars            : array to store/copy the involved variables of the constraint
 *  - success         : pointer to store whether the variables are successfully copied
 */
#define SCIP_DECL_CONSGETVARS(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons, \
      SCIP_VAR** vars, int varssize, SCIP_Bool* success)

/** constraint method of constraint handler which returns the number of variables (if possible)
 *
 *  The constraint handler can (this callback is optional) provide this callback to return the number variable which are
 *  involved in that particular constraint. If this is not possible, the success pointers has to be set to FALSE or the
 *  callback should not be implemented.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : constraint for which the number of variables is wanted
 *
 *  output:
 *  - nvars           : pointer to store the number of variables
 *  - success         : pointer to store whether the constraint successfully returned the number of variables
 */
#define SCIP_DECL_CONSGETNVARS(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_CONS* cons, \
      int* nvars, SCIP_Bool* success)

/** constraint handler method to suggest dive bound changes during the generic diving algorithm
 *
 *  This callback is used inside the various diving heuristics of SCIP and does not affect the normal branching of the
 *  actual search. The constraint handler can provide this callback to render the current solution (even more)
 *  infeasible by suggesting one or several variable bound changes. In fact, since diving heuristics do not necessarily
 *  solve LP relaxations at every probing depth, some of the variable local bounds might already be conflicting with the
 *  solution values.  The solution is rendered infeasible by determining bound changes that should be applied to the
 *  next explored search node via SCIPaddDiveBoundChange().  An alternative in case that the preferred bound change(s)
 *  were detected infeasible must be provided.
 *
 *  The constraint handler must take care to only add bound changes that further shrink the variable domain.
 *
 *  The success pointer must be used to indicate whether the constraint handler succeeded in selecting diving bound
 *  changes. The infeasible pointer should be set to TRUE if the constraint handler found a local infeasibility.  If the
 *  constraint handler needs to select between several candidates, it may use the scoring mechanism of the diveset
 *  argument to control its choice.
 *
 *  This callback is optional.
 *
 *  @note: @p sol is usually the LP relaxation solution unless the caller of the method, usually a diving heuristic,
 *         does not solve LP relaxations at every depth
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - diveset         : diving settings for scoring
 *  - sol             : current diving solution, usually the LP relaxation solution
 *
 *  output:
 *  - success         : pointer to store whether the constraint handler succeeded to determine dive bound changes
 *  - infeasible      : pointer to store whether the constraint handler detected an infeasibility in the local node
 */
#define SCIP_DECL_CONSGETDIVEBDCHGS(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONSHDLR* conshdlr, SCIP_DIVESET* diveset, \
      SCIP_SOL* sol, SCIP_Bool* success, SCIP_Bool* infeasible)

#ifdef __cplusplus
}
#endif

#endif
