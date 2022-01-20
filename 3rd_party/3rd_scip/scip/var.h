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

/**@file   var.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VAR_H__
#define __SCIP_VAR_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_misc.h"
#include "scip/type_history.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_branch.h"
#include "scip/type_cons.h"
#include "scip/pub_var.h"
#include "scip/pub_misc.h"

#ifndef NDEBUG
#include "scip/struct_var.h"
#else
#include "scip/event.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * domain change methods
 */

/** applies single bound change */
extern
SCIP_RETCODE SCIPboundchgApply(
   SCIP_BOUNDCHG*        boundchg,           /**< bound change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   depth,              /**< depth in the tree, where the bound change takes place */
   int                   pos,                /**< position of the bound change in its bound change array */
   SCIP_Bool*            cutoff              /**< pointer to store whether an infeasible bound change was detected */
   );

/** undoes single bound change */
extern
SCIP_RETCODE SCIPboundchgUndo(
   SCIP_BOUNDCHG*        boundchg,           /**< bound change to remove */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** applies domain change to the global problem */
extern
SCIP_RETCODE SCIPdomchgApplyGlobal(
   SCIP_DOMCHG*          domchg,             /**< domain change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether an infeasible domain change was detected */
   );

/** frees domain change data */
extern
SCIP_RETCODE SCIPdomchgFree(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** converts a dynamic domain change data into a static one, using less memory than for a dynamic one */
extern
SCIP_RETCODE SCIPdomchgMakeStatic(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** applies domain change */
extern
SCIP_RETCODE SCIPdomchgApply(
   SCIP_DOMCHG*          domchg,             /**< domain change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   depth,              /**< depth in the tree, where the domain change takes place */
   SCIP_Bool*            cutoff              /**< pointer to store whether an infeasible domain change was detected */
   );

/** undoes domain change */
extern
SCIP_RETCODE SCIPdomchgUndo(
   SCIP_DOMCHG*          domchg,             /**< domain change to remove */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** adds bound change to domain changes */
extern
SCIP_RETCODE SCIPdomchgAddBoundchg(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to change the bounds for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   SCIP_BOUNDCHGTYPE     boundchgtype,       /**< type of bound change: branching decision or inference */
   SCIP_Real             lpsolval,           /**< solval of variable in last LP on path to node, or SCIP_INVALID if unknown */
   SCIP_VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change, or NULL */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change, or NULL */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_BOUNDTYPE        inferboundtype      /**< type of bound for inference var: lower or upper bound */
   );

/** adds hole change to domain changes */
extern
SCIP_RETCODE SCIPdomchgAddHolechg(
   SCIP_DOMCHG**         domchg,             /**< pointer to domain change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HOLELIST**       ptr,                /**< changed list pointer */
   SCIP_HOLELIST*        newlist,            /**< new value of list pointer */
   SCIP_HOLELIST*        oldlist             /**< old value of list pointer */
   );

/** creates an artificial bound change information object with depth = INT_MAX and pos = -1 */
extern
SCIP_RETCODE SCIPbdchginfoCreate(
   SCIP_BDCHGINFO**      bdchginfo,          /**< pointer to store bound change information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_VAR*             var,                /**< active variable that changed the bounds */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound            /**< new value for bound */
   );

/** frees a bound change information object */
extern
void SCIPbdchginfoFree(
   SCIP_BDCHGINFO**      bdchginfo,          /**< pointer to store bound change information */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns the relaxed bound change type */
extern
SCIP_Real SCIPbdchginfoGetRelaxedBound(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict set */
   );

/*
 * methods for variables 
 */

/** creates and captures an original problem variable; an integer variable with bounds
 *  zero and one is automatically converted into a binary variable
 */
extern
SCIP_RETCODE SCIPvarCreateOriginal(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< name of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );

/** creates and captures a loose variable belonging to the transformed problem; an integer variable with bounds
 *  zero and one is automatically converted into a binary variable
 */
extern
SCIP_RETCODE SCIPvarCreateTransformed(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< name of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );

/** copies and captures a variable from source to target SCIP; an integer variable with bounds zero and one is
 *  automatically converted into a binary variable; in case the variable data cannot be copied the variable is not
 *  copied at all
 */
extern
SCIP_RETCODE SCIPvarCopy(
   SCIP_VAR**            var,                /**< pointer to store the target variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_VAR*             sourcevar,          /**< source variable */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool             global              /**< should global or local bounds be used? */
   );

/** parses variable information (in cip format) out of a string; if the parsing process was successful an original
 *  variable is created and captured; if variable is of integral type, fractional bounds are automatically rounded; an
 *  integer variable with bounds zero and one is automatically converted into a binary variable
 */
extern
SCIP_RETCODE SCIPvarParseOriginal(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           str,                /**< string to parse */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata,            /**< user data for this specific variable */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   );

/** parses variable information (in cip format) out of a string; if the parsing process was successful a loose variable
 *  belonging to the transformed problem is created and captured; if variable is of integral type, fractional bounds are
 *  automatically rounded; an integer variable with bounds zero and one is automatically converted into a binary
 *  variable
 */
extern
SCIP_RETCODE SCIPvarParseTransformed(
   SCIP_VAR**            var,                /**< pointer to variable data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           str,                /**< string to parse */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata,            /**< user data for this specific variable */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   );

/** increases usage counter of variable */
extern
void SCIPvarCapture(
   SCIP_VAR*             var                 /**< variable */
   );

/** decreases usage counter of variable, and frees memory if necessary */
extern
SCIP_RETCODE SCIPvarRelease(
   SCIP_VAR**            var,                /**< pointer to variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data (may be NULL, if it's not a column variable) */
   );

/** change variable name */
extern
/** change variable name */
SCIP_RETCODE SCIPvarChgName(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name                /**< name of variable */
   );

/** initializes variable data structure for solving */
extern
void SCIPvarInitSolve(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
extern
SCIP_RETCODE SCIPvarTransform(
   SCIP_VAR*             origvar,            /**< original problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets corresponding transformed variable of an original or negated original variable */
extern
SCIP_RETCODE SCIPvarGetTransformed(
   SCIP_VAR*             origvar,            /**< original problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable, or NULL if not existing yet */
   );

/** converts transformed variable into column variable and creates LP column */
extern
SCIP_RETCODE SCIPvarColumn(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** converts column transformed variable back into loose variable, frees LP column */
extern
SCIP_RETCODE SCIPvarLoose(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp                  /**< current LP data */
   );

/** converts variable into fixed variable */
extern
SCIP_RETCODE SCIPvarFix(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             fixedval,           /**< value to fix variable at */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   );

/** transforms given variables, scalars and constant to the corresponding active variables, scalars and constant
 *
 * If the number of needed active variables is greater than the available slots in the variable array, nothing happens except
 * that the required size is stored in the corresponding variable; hence, if afterwards the required size is greater than the
 * available slots (varssize), nothing happens; otherwise, the active variable representation is stored in the arrays.
 *
 * The reason for this approach is that we cannot reallocate memory, since we do not know how the
 * memory has been allocated (e.g., by a C++ 'new' or SCIP functions).
 */
extern
SCIP_RETCODE SCIPvarGetActiveRepresentatives(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< variable array to get active variables */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   int                   varssize,           /**< available slots in vars and scalars array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   SCIP_Bool             mergemultiples      /**< should multiple occurrences of a var be replaced by a single coeff? */
   );

/** transforms given variable, scalar and constant to the corresponding active, fixed, or
 *  multi-aggregated variable, scalar and constant; if the variable resolves to a fixed variable,
 *  "scalar" will be 0.0 and the value of the sum will be stored in "constant"; a multi-aggregation
 *  with only one active variable (this can happen due to fixings after the multi-aggregation),
 *  is treated like an aggregation; if the multi-aggregation constant is infinite, "scalar" will be 0.0
 */
extern
SCIP_RETCODE SCIPvarGetProbvarSum(
   SCIP_VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   SCIP_Real*            constant            /**< pointer to constant c in sum a*x + c */
   );


/** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later-on */
extern
SCIP_RETCODE SCIPvarFlattenAggregationGraph(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** return for given variables all their active counterparts; all active variables will be pairwise different */
extern
SCIP_RETCODE SCIPvarsGetActiveVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< variable array with given variables and as output all active
                                              *   variables, if enough slots exist
                                              */
   int*                  nvars,              /**< number of given variables, and as output number of active variables,
                                              *   if enough slots exist
                                              */
   int                   varssize,           /**< available slots in vars array */
   int*                  requiredsize        /**< pointer to store the required array size for the active variables */
   );

/** performs second step of SCIPaggregateVars():
 *  the variable to be aggregated is chosen among active problem variables x' and y', preferring a less strict variable
 *  type as aggregation variable (i.e. continuous variables are preferred over implicit integers, implicit integers
 *  or integers over binaries). If none of the variables is continuous, it is tried to find an integer
 *  aggregation (i.e. integral coefficients a'' and b'', such that a''*x' + b''*y' == c''). This can lead to
 *  the detection of infeasibility (e.g. if c'' is fractional), or to a rejection of the aggregation (denoted by
 *  aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 */
extern
SCIP_RETCODE SCIPvarTryAggregateVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             varx,               /**< variable x in equality a*x + b*y == c */
   SCIP_VAR*             vary,               /**< variable y in equality a*x + b*y == c */
   SCIP_Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   SCIP_Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   SCIP_Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   );

/** converts loose variable into aggregated variable */
extern
SCIP_RETCODE SCIPvarAggregate(
   SCIP_VAR*             var,                /**< loose problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             aggvar,             /**< loose variable y in aggregation x = a*y + c */
   SCIP_Real             scalar,             /**< multiplier a in aggregation x = a*y + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a*y + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   );

/** converts variable into multi-aggregated variable */
extern
SCIP_RETCODE SCIPvarMultiaggregate(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   );

/** returns whether variable is not allowed to be multi-aggregated */
extern
SCIP_Bool SCIPvarDoNotMultaggr(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns whether the variable was flagged for deletion from global structures (cliques etc.) */
SCIP_Bool SCIPvarIsMarkedDeleteGlobalStructures(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets negated variable x' = offset - x of problem variable x; the negated variable is created if not yet existing;
 *  the negation offset of binary variables is always 1, the offset of other variables is fixed to lb + ub when the
 *  negated variable is created
 */
extern
SCIP_RETCODE SCIPvarNegate(
   SCIP_VAR*             var,                /**< problem variable to negate */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR**            negvar              /**< pointer to store the negated variable */
   );

/** informs variable that its position in problem's vars array changed */
extern
void SCIPvarSetProbindex(
   SCIP_VAR*             var,                /**< problem variable */
   int                   probindex           /**< new problem index of variable */
   );

/** gives the variable a new name
 *
 *  @note the old pointer is overwritten, which might result in a memory leakage
 */
extern
void SCIPvarSetNamePointer(
   SCIP_VAR*             var,                /**< problem variable */
   const char*           name                /**< new name of variable */
   );

/** informs variable that it will be removed from the problem; adjusts probindex and removes variable from the
 *  implication graph;
 *  If 'final' is TRUE, the thorough implication graph removal is not performed. Instead, only the
 *  variable bounds and implication data structures of the variable are freed. Since in the final removal
 *  of all variables from the transformed problem, this deletes the implication graph completely and is faster
 *  than removing the variables one by one, each time updating all lists of the other variables.
 */
extern
SCIP_RETCODE SCIPvarRemove(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             final               /**< is this the final removal of all problem variables? */
   );

/** marks the variable to be deleted from the problem */
extern
void SCIPvarMarkDeleted(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** marks the variable to not to be multi-aggregated */
extern
SCIP_RETCODE SCIPvarMarkDoNotMultaggr(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** modifies lock numbers for rounding */
extern
SCIP_RETCODE SCIPvarAddLocks(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   int                   addnlocksdown,      /**< increase in number of rounding down locks */
   int                   addnlocksup         /**< increase in number of rounding up locks */
   );

/** changes type of variable; cannot be called, if var belongs to a problem */
extern
SCIP_RETCODE SCIPvarChgType(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_VARTYPE          vartype             /**< new type of variable */
   );

/** changes objective value of variable */
extern
SCIP_RETCODE SCIPvarChgObj(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             newobj              /**< new objective value for variable */
   );

/** adds value to objective value of variable */
extern
SCIP_RETCODE SCIPvarAddObj(
   SCIP_VAR*             var,                /**< variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Real             addobj              /**< additional objective value for variable */
   );

/** changes objective value of variable in current dive */
extern
SCIP_RETCODE SCIPvarChgObjDive(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newobj              /**< new objective value for variable */
   );

/** adjust lower bound to integral value, if variable is integral */
extern
void SCIPvarAdjustLb(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            lb                  /**< pointer to lower bound to adjust */
   );

/** adjust upper bound to integral value, if variable is integral */
extern
void SCIPvarAdjustUb(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real*            ub                  /**< pointer to upper bound to adjust */
   );

/** adjust lower or upper bound to integral value, if variable is integral */
extern
void SCIPvarAdjustBd(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound to adjust */
   SCIP_Real*            bd                  /**< pointer to bound to adjust */
   );

/** changes lower bound of original variable in original problem */
extern
SCIP_RETCODE SCIPvarChgLbOriginal(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** changes upper bound of original variable in original problem */
extern
SCIP_RETCODE SCIPvarChgUbOriginal(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** changes global lower bound of variable; if possible, adjusts bound to integral value;
 *  updates local lower bound if the global bound is tighter
 */
extern
SCIP_RETCODE SCIPvarChgLbGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** changes global upper bound of variable; if possible, adjusts bound to integral value;
 *  updates local upper bound if the global bound is tighter
 */
extern
SCIP_RETCODE SCIPvarChgUbGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** changes global bound of variable; if possible, adjusts bound to integral value;
 *  updates local bound if the global bound is tighter
 */
extern
SCIP_RETCODE SCIPvarChgBdGlobal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Real             newbound,           /**< new bound for variable */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** changes current local lower bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
extern
SCIP_RETCODE SCIPvarChgLbLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** changes current local upper bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
extern
SCIP_RETCODE SCIPvarChgUbLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** changes current local bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
extern
SCIP_RETCODE SCIPvarChgBdLocal(
   SCIP_VAR*             var,                /**< problem variable to change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data, may be NULL for original variables */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage, may be NULL for original variables */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             newbound,           /**< new bound for variable */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** changes lazy lower bound of the variable, this is only possible if the variable is not in the LP yet */
extern
SCIP_RETCODE SCIPvarChgLbLazy(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             lazylb              /**< the lazy lower bound to be set */
   );

/** changes lazy upper bound of the variable, this is only possible if the variable is not in the LP yet */
SCIP_RETCODE SCIPvarChgUbLazy(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             lazylb              /**< the lazy lower bound to be set */
   );

/** changes lower bound of variable in current dive; if possible, adjusts bound to integral value */
extern
SCIP_RETCODE SCIPvarChgLbDive(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** changes upper bound of variable in current dive; if possible, adjusts bound to integral value */
extern
SCIP_RETCODE SCIPvarChgUbDive(
   SCIP_VAR*             var,                /**< problem variable to change */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Real             newbound            /**< new bound for variable */
   );

/** for a multi-aggregated variable, gives the local lower bound computed by adding the local bounds from all aggregation variables
 * this lower bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable is not allowed
 */
extern
SCIP_Real SCIPvarGetMultaggrLbLocal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** for a multi-aggregated variable, gives the local upper bound computed by adding the local bounds from all aggregation variables
 * this upper bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable is not allowed
 */
extern
SCIP_Real SCIPvarGetMultaggrUbLocal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** for a multi-aggregated variable, gives the global lower bound computed by adding the global bounds from all aggregation variables
 * this global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable is not allowed
 */
extern
SCIP_Real SCIPvarGetMultaggrLbGlobal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** for a multi-aggregated variable, gives the global upper bound computed by adding the global bounds from all aggregation variables
 * this upper bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable is not allowed
 */
extern
SCIP_Real SCIPvarGetMultaggrUbGlobal(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** adds a hole to the original domain of the variable*/
extern
SCIP_RETCODE SCIPvarAddHoleOriginal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right               /**< right bound of open interval in new hole */
   );

/** adds a hole to the variable's global domain */
extern
SCIP_RETCODE SCIPvarAddHoleGlobal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right,              /**< right bound of open interval in new hole */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added, or NULL */
   );

/** adds a hole to the variable's current local domain */
extern
SCIP_RETCODE SCIPvarAddHoleLocal(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue, may be NULL for original variables */
   SCIP_Real             left,               /**< left bound of open interval in new hole */
   SCIP_Real             right,              /**< right bound of open interval in new hole */
   SCIP_Bool*            added               /**< pointer to store whether the hole was added, or NULL */
   );

/** resets the global and local bounds of original variable to their original values */
extern
SCIP_RETCODE SCIPvarResetBounds(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** returns at which depth in the tree a bound change was applied to the variable that conflicts with the
 *  given bound; returns -1 if the bound does not conflict with the current local bounds of the variable
 */
extern
int SCIPvarGetConflictingBdchgDepth(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        boundtype,          /**< bound type of the conflicting bound */
   SCIP_Real             bound               /**< conflicting bound */
   );

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  improves the global bounds of the variable and the vlb variable if possible
 */
extern
SCIP_RETCODE SCIPvarAddVlb(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   SCIP_Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   SCIP_Real             vlbconstant,        /**< constant d    in x >= b*z + d */
   SCIP_Bool             transitive,         /**< should transitive closure of implication also be added? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  updates the global bounds of the variable and the vub variable correspondingly
 */
extern
SCIP_RETCODE SCIPvarAddVub(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   SCIP_Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   SCIP_Real             vubconstant,        /**< constant d    in x <= b*z + d */
   SCIP_Bool             transitive,         /**< should transitive closure of implication also be added? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );

/** informs binary variable x about a globally valid implication:  x == 0 or x == 1  ==>  y <= b  or  y >= b;
 *  also adds the corresponding implication or variable bound to the implied variable;
 *  if the implication is conflicting, the variable is fixed to the opposite value;
 *  if the variable is already fixed to the given value, the implication is performed immediately;
 *  if the implication is redundant with respect to the variables' global bounds, it is ignored
 */
extern
SCIP_RETCODE SCIPvarAddImplic(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool             transitive,         /**< should transitive closure of implication also be added? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   );

/** fixes the bounds of a binary variable to the given value, counting bound changes and detecting infeasibility */
extern
SCIP_RETCODE SCIPvarFixBinary(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             value,              /**< value to fix variable to */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   );

/** adds the variable to the given clique and updates the list of cliques the binary variable is member of;
 *  if the variable now appears twice in the clique with the same value, it is fixed to the opposite value;
 *  if the variable now appears twice in the clique with opposite values, all other variables are fixed to
 *  the opposite of the value they take in the clique
 */
extern
SCIP_RETCODE SCIPvarAddClique(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree if in solving stage */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique,             /**< clique the variable should be added to */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to count the number of performed bound changes, or NULL */
   );

/** adds a filled clique to the cliquelists of all corresponding variables */
extern
SCIP_RETCODE SCIPvarsAddClique(
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Bool*            values,             /**< values of the variables in the clique */
   int                   nvars,              /**< number of problem variables */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CLIQUE*          clique              /**< clique that contains all given variables and values */
   );

/** deletes the variable from the given clique and updates the list of cliques the binary variable is member of */
extern
SCIP_RETCODE SCIPvarDelClique(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique              /**< clique the variable should be removed from */
   );

/** adds a clique to the list of cliques of the given binary variable, but does not change the clique
 *  itself
 */
extern
SCIP_RETCODE SCIPvarAddCliqueToList(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique              /**< clique that should be removed from the variable's clique list */
   );

/** deletes a clique from the list of cliques the binary variable is member of, but does not change the clique
 *  itself
 */
extern
SCIP_RETCODE SCIPvarDelCliqueFromList(
   SCIP_VAR*             var,                /**< problem variable  */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             value,              /**< value of the variable in the clique */
   SCIP_CLIQUE*          clique              /**< clique that should be removed from the variable's clique list */
   );

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
extern
SCIP_RETCODE SCIPvarChgBranchFactor(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             branchfactor        /**< factor to weigh variable's branching score with */
   );

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 */
extern
SCIP_RETCODE SCIPvarChgBranchPriority(
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< branching priority of the variable */
   );

/** sets the branch direction of the variable; variables with higher branch direction are always preferred to variables
 *  with lower direction in selection of branching variable
 */
extern
SCIP_RETCODE SCIPvarChgBranchDirection(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   );

/** gets objective value of variable in current SCIP_LP; the value can be different from the bound stored in the variable's own
 *  data due to diving, that operate only on the LP without updating the variables
 */
extern
SCIP_Real SCIPvarGetObjLP(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** gets lower bound of variable in current SCIP_LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
extern
SCIP_Real SCIPvarGetLbLP(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets upper bound of variable in current SCIP_LP; the bound can be different from the bound stored in the variable's own
 *  data due to diving or conflict analysis, that operate only on the LP without updating the variables
 */
extern
SCIP_Real SCIPvarGetUbLP(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns solution value and index of variable lower bound that is closest to the variable's value in the given primal solution
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable lower bound is available
 */
extern
void SCIPvarGetClosestVlb(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   );

/** returns solution value and index of variable upper bound that is closest to the variable's value in the given primal solution;
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is available
 */
extern
void SCIPvarGetClosestVub(
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            closestvub,         /**< pointer to store the value of the closest variable upper bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest variable upper bound */
   );

/** remembers the current solution as root solution in the problem variables */
extern
void SCIPvarStoreRootSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             roothaslp           /**< is the root solution from LP? */
   );

/** updates the current solution as best root solution in the problem variables if it is better */
extern
void SCIPvarUpdateBestRootSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             rootsol,            /**< root solution value */
   SCIP_Real             rootredcost,        /**< root reduced cost */
   SCIP_Real             rootlpobjval        /**< objective value of the root LP */
   );

/** returns the solution value of the problem variables in the relaxation solution */
extern
SCIP_Real SCIPvarGetRelaxSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** returns the solution value of the transformed problem variable in the relaxation solution */
extern
SCIP_Real SCIPvarGetRelaxSolTransVar(
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns for the given binary variable the reduced cost which are given by the variable itself and its implication if
 *  the binary variable is fixed to the given value
 */
extern
SCIP_Real SCIPvarGetImplRedcost(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             varfixing,          /**< FALSE if for x == 0, TRUE for x == 1 */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem, or NULL */
   SCIP_LP*              lp                  /**< current LP data */
   );


/** stores the solution value as relaxation solution in the problem variable */
extern
SCIP_RETCODE SCIPvarSetRelaxSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Real             solval,             /**< solution value in the current relaxation solution */
   SCIP_Bool             updateobj           /**< should the objective value be updated? */
   );

/** stores the solution value as NLP solution in the problem variable */
extern
SCIP_RETCODE SCIPvarSetNLPSol(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             solval              /**< solution value in the current NLP solution */
   );

/** resolves variable to columns and adds them with the coefficient to the row */
extern
SCIP_RETCODE SCIPvarAddToRow(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             val                 /**< value of coefficient */
   );


/** merge two variable histories together; a typical use case is that \p othervar is an image of the target variable
 *  in a SCIP copy. Method should be applied with care, especially because no internal checks are performed whether
 *  the history merge is reasonable
 *
 *  @note Do not use this method if the two variables originate from two SCIP's with different objective functions, since
 *        this corrupts the variable pseudo costs
 *  @note Apply with care; no internal checks are performed if the two variables should be merged
 */
extern
void SCIPvarMergeHistories(
   SCIP_VAR*             targetvar,          /**< the variable that should contain both histories afterwards */
   SCIP_VAR*             othervar,           /**< the variable whose history is to be merged with that of the target variable */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** sets the history of a variable; this method is typacally used within reoptimization to keep and update the variable
 *  history over several iteraions
 */
extern
void SCIPvarSetHistory(
   SCIP_VAR*             var,                /**< variable */
   SCIP_HISTORY*         history,            /**< the history which is to set */
   SCIP_STAT*            stat                /**< problem statistics */
   );

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of
 *  "solvaldelta" in the variable's solution value and resulting change of "objdelta" in the in the LP's objective value
 */
extern
SCIP_RETCODE SCIPvarUpdatePseudocost(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   );

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value */
extern
SCIP_Real SCIPvarGetPseudocost(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPvarGetPseudocostCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
extern
SCIP_Real SCIPvarGetPseudocostCount(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 */
extern
SCIP_Real SCIPvarGetPseudocostCountCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );


/** compares both possible directions for rounding the given solution value and returns the minimum pseudo-costs of the variable */
SCIP_Real SCIPvarGetMinPseudocostScore(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             solval              /**< solution value, e.g., LP solution value */
   );

/** gets the an estimate of the variable's pseudo cost variance in direction \p dir */
extern
SCIP_Real SCIPvarGetPseudocostVariance(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Bool             onlycurrentrun      /**< return pseudo cost variance only for current branch and bound run */
   );

/** calculates a confidence bound for this variable under the assumption of normally distributed pseudo costs
 *
 *  The confidence bound \f$ \theta \geq 0\f$ denotes the interval borders \f$ [X - \theta, \ X + \theta]\f$, which contains
 *  the true pseudo costs of the variable, i.e., the expected value of the normal distribution, with a probability
 *  of 95 %.
 *
 *  @return value of confidence bound for this variable
 */
extern
SCIP_Real SCIPvarCalcPscostConfidenceBound(
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BRANCHDIR        dir,                /**< the branching direction for the confidence bound */
   SCIP_Bool             onlycurrentrun,     /**< should only the current run be taken into account */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for the interval */
   );

/** check if the current pseudo cost relative error in a direction violates the given threshold. The Relative
 *  Error is calculated at a specific confidence level
 */
extern
SCIP_Bool SCIPvarIsPscostRelerrorReliable(
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             threshold,          /**< threshold for relative errors to be considered reliable (enough) */
   SCIP_CONFIDENCELEVEL  clevel              /**< a given confidence level */
   );

/** check if variable pseudo-costs have a significant difference in location. The significance depends on
 *  the choice of \p clevel and on the kind of tested hypothesis. The one-sided hypothesis, which
 *  should be rejected, is that fracy * mu_y >= fracx * mu_x, where mu_y and mu_x denote the
 *  unknown location means of the underlying pseudo-cost distributions of x and y.
 *
 *  This method is applied best if variable x has a better pseudo-cost score than y. The method hypothesizes that y were actually
 *  better than x (despite the current information), meaning that y can be expected to yield branching
 *  decisions as least as good as x in the long run. If the method returns TRUE, the current history information is
 *  sufficient to safely rely on the alternative hypothesis that x yields indeed a better branching score (on average)
 *  than y.
 *
 *  @note The order of x and y matters for the one-sided hypothesis
 *
 *  @note set \p onesided to FALSE if you are not sure which variable is better. The hypothesis tested then reads
 *        fracy * mu_y == fracx * mu_x vs the alternative hypothesis fracy * mu_y != fracx * mu_x.
 *
 *  @return TRUE if the hypothesis can be safely rejected at the given confidence level
 */
extern
SCIP_Bool SCIPvarSignificantPscostDifference(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             varx,               /**< variable x */
   SCIP_Real             fracx,              /**< the fractionality of variable x */
   SCIP_VAR*             vary,               /**< variable y */
   SCIP_Real             fracy,              /**< the fractionality of variable y */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel,             /**< confidence level for rejecting hypothesis */
   SCIP_Bool             onesided            /**< should a one-sided hypothesis y >= x be tested? */
   );

/** tests at a given confidence level whether the variable pseudo-costs only have a small probability to
 *  exceed a \p threshold. This is useful to determine if past observations provide enough evidence
 *  to skip an expensive strong-branching step if there is already a candidate that has been proven to yield an improvement
 *  of at least \p threshold.
 *
 *  @note use \p clevel to adjust the level of confidence. For SCIP_CONFIDENCELEVEL_MIN, the method returns TRUE if
 *        the estimated probability to exceed \p threshold is less than 25 %.
 *
 *  @see  SCIP_Confidencelevel for a list of available levels. The used probability limits refer to the one-sided levels
 *        of confidence.
 *
 *  @return TRUE if the variable pseudo-cost probabilistic model is likely to be smaller than \p threshold
 *          at the given confidence level \p clevel.
 */
extern
SCIP_Bool SCIPvarPscostThresholdProbabilityTest(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable x */
   SCIP_Real             frac,               /**< the fractionality of variable x */
   SCIP_Real             threshold,          /**< the threshold to test against */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for rejecting hypothesis */
   );

/** increases the VSIDS of the variable by the given weight */
extern
SCIP_RETCODE SCIPvarIncVSIDS(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             weight              /**< weight of this update in conflict score */
   );

/** scales the VSIDS of the variable by the given scalar */
extern
SCIP_RETCODE SCIPvarScaleVSIDS(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             scalar              /**< scalar to multiply the conflict scores with */
   );

/** increases the number of active conflicts by one and the overall length of the variable by the given length */
SCIP_RETCODE SCIPvarIncNActiveConflicts(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             length              /**< length of the conflict */
   );

/** gets the number of active conflicts containing this variable in given direction */
SCIP_Longint SCIPvarGetNActiveConflicts(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** gets the number of active conflicts containing this variable in given direction
 *  in the current run
 */
SCIP_Longint SCIPvarGetNActiveConflictsCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/**  gets the average conflict length in given direction due to branching on the variable */
SCIP_Real SCIPvarGetAvgConflictlength(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/**  gets the average conflict length in given direction due to branching on the variable
 *   in the current run
 */
SCIP_Real SCIPvarGetAvgConflictlengthCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** increases the number of branchings counter of the variable */
extern
SCIP_RETCODE SCIPvarIncNBranchings(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   int                   depth               /**< depth at which the bound change took place */
   );

/** increases the inference score of the variable by the given weight */
extern
SCIP_RETCODE SCIPvarIncInferenceSum(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             weight              /**< weight of this update in inference score */
   );

/** increases the cutoff score of the variable by the given weight */
extern
SCIP_RETCODE SCIPvarIncCutoffSum(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_SET*             set,                /**< global SCIP settings, or NULL if the domain value is SCIP_UNKNOWN */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             weight              /**< weight of this update in cutoff score */
   );

/** returns the variable's VSIDS score */
extern
SCIP_Real SCIPvarGetVSIDS_rec(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the variable's VSIDS score only using conflicts of the current run */
extern
SCIP_Real SCIPvarGetVSIDSCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching on the variable in given direction */
extern
SCIP_Real SCIPvarGetAvgInferences(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences found after branching on the variable in given direction
 *  in the current run
 */
extern
SCIP_Real SCIPvarGetAvgInferencesCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of cutoffs found after branching on the variable in given direction */
extern
SCIP_Real SCIPvarGetAvgCutoffs(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run */
extern
SCIP_Real SCIPvarGetAvgCutoffsCurrentRun(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** outputs variable information into file stream */
extern
SCIP_RETCODE SCIPvarPrint(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** includes event handler with given data in variable's event filter */
extern
SCIP_RETCODE SCIPvarCatchEvent(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   );

/** deletes event handler with given data from variable's event filter */
extern
SCIP_RETCODE SCIPvarDropEvent(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler for the event processing */
   int                   filterpos           /**< position of event filter entry returned by SCIPvarCatchEvent(), or -1 */
   );

/** returns the variable's VSIDS score */
extern
SCIP_Real SCIPvarGetVSIDS(
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the position of the bound change index */
extern
int SCIPbdchgidxGetPos(
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index */
   );

/** removes (redundant) cliques, implications and variable bounds of variable from all other variables' implications and variable
 *  bounds arrays, and optionally removes them also from the variable itself
 */
extern
SCIP_RETCODE SCIPvarRemoveCliquesImplicsVbs(
   SCIP_VAR*             var,                /**< problem variable */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             irrelevantvar,      /**< has the variable become irrelevant? */
   SCIP_Bool             onlyredundant,      /**< should only the redundant implications and variable bounds be removed? */
   SCIP_Bool             removefromvar       /**< should the implications and variable bounds be removed from the var itself? */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPvarCatchEvent(var, blkmem, set, eventtype, eventhdlr, eventdata, filterpos) \
   SCIPeventfilterAdd(var->eventfilter, blkmem, set, eventtype, eventhdlr, eventdata, filterpos)
#define SCIPvarDropEvent(var, blkmem, set, eventtype, eventhdlr, eventdata, filterpos) \
   SCIPeventfilterDel(var->eventfilter, blkmem, set, eventtype, eventhdlr, eventdata, filterpos)
#define SCIPvarGetVSIDS(var, stat, dir)    ((var)->varstatus == SCIP_VARSTATUS_LOOSE || (var)->varstatus == SCIP_VARSTATUS_COLUMN ? \
      SCIPhistoryGetVSIDS(var->history, dir)/stat->vsidsweight : SCIPvarGetVSIDS_rec(var, stat, dir))
#define SCIPbdchgidxGetPos(bdchgidx) ((bdchgidx)->pos)

#endif

/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
extern
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyVar);

#ifdef __cplusplus
}
#endif

#endif
