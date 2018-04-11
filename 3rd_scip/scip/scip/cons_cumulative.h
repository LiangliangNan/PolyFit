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

/**@file   cons_cumulative.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for cumulative constraints
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_CUMULATIVE_H__
#define __SCIP_CONS_CUMULATIVE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the constraint handler for cumulative constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrCumulative(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Cumulative Constraints
 *
 * Given:
 * - a set of jobs, represented by their integer start time variables \f$S_j\f$, their array of processing times \f$p_j\f$ and of
 *   their demands \f$d_j\f$.
 * - an integer resource capacity \f$C\f$
 *
 * The cumulative constraint ensures that for each point in time \f$t\f$ \f$\sum_{j: S_j \leq t < S_j + p_j} d_j \leq C\f$ holds.
 *
 * @par
 * Separation:
 * - can be done using binary start time model, see Pritskers, Watters and Wolfe
 * - or by just separating relatively weak cuts on the start time variables
 *
 * @par
 * Propagation:
 * - time tabling, Klein & Scholl (1999)
 * - Edge-finding from Petr Vilim, adjusted and simplified for dynamic repropagation
 *   (2009)
 * - energetic reasoning, see Baptiste, Le Pape, Nuijten (2001)
 *
 * @{
 */

/** creates and captures a cumulative constraint */
EXTERN
SCIP_RETCODE SCIPcreateConsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures an absolute power constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsCumulative(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsCumulative() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity            /**< available cumulative capacity */
   );

/** set the left bound of effective horizon */
EXTERN
SCIP_RETCODE SCIPsetHminCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   int                   hmin                /**< left bound of time axis to be considered */
   );

/** returns the left bound of the effective horizon */
EXTERN
int SCIPgetHminCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );


/** set the right bound of the effective horizon */
EXTERN
SCIP_RETCODE SCIPsetHmaxCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   int                   hmax                /**< right bound of time axis to be considered */
   );

/** returns the right bound of effective horizon */
EXTERN
int SCIPgetHmaxCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the start time variables of the cumulative constraint */
EXTERN
SCIP_VAR** SCIPgetVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the number of start time variables of the cumulative constraint */
EXTERN
int SCIPgetNVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the capacity of the cumulative constraint */
EXTERN
int SCIPgetCapacityCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the durations of the cumulative constraint */
EXTERN
int* SCIPgetDurationsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the demands of the cumulative constraint */
EXTERN
int* SCIPgetDemandsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** check for the given starting time variables with their demands and durations if the cumulative conditions for the
 *  given solution is satisfied
 */
EXTERN
SCIP_RETCODE SCIPcheckCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered */
   int                   hmax,               /**< right bound of time axis to be considered */
   SCIP_Bool*            violated,           /**< pointer to store if the cumulative condition is violated */
   SCIP_CONS*            cons,               /**< constraint which is checked */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   );

/** normalize cumulative condition */
EXTERN
SCIP_RETCODE SCIPnormalizeCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int*                  capacity,           /**< pointer to store the changed cumulative capacity */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   );

/** searches for a time point within the cumulative condition were the cumulative condition can be split */
EXTERN
SCIP_RETCODE SCIPsplitCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   int*                  hmin,               /**< pointer to store the left bound of the effective horizon */
   int*                  hmax,               /**< pointer to store the right bound of the effective horizon */
   int*                  split               /**< point were the cumulative condition can be split */
   );

/** presolve cumulative condition w.r.t. effective horizon by detecting irrelevant variables */
EXTERN
SCIP_RETCODE SCIPpresolveCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int                   hmin,               /**< left bound of time axis to be considered */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_Bool*            downlocks,          /**< array storing if the variable has a down lock, or NULL */
   SCIP_Bool*            uplocks,            /**< array storing if the variable has an up lock, or NULL */
   SCIP_CONS*            cons,               /**< constraint which gets propagated, or NULL */
   SCIP_Bool*            delvars,            /**< array storing the variable which can be deleted from the constraint */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   int*                  nchgsides,          /**< pointer to store the number of changed sides */
   SCIP_Bool*            cutoff              /**< buffer to store whether a cutoff is detected */
   );

/** propagate the given cumulative condition */
EXTERN
SCIP_RETCODE SCIPpropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLTIMING     presoltiming,       /**< current presolving timing */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered */
   int                   hmax,               /**< right bound of time axis to be considered */
   SCIP_CONS*            cons,               /**< constraint which gets propagated */
   int*                  nchgbds,            /**< pointer to store the number of variable bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_Bool*            cutoff              /**< pointer to store if the cumulative condition is violated */
   );

/** resolve propagation w.r.t. the cumulative condition */
EXTERN
SCIP_RETCODE SCIPrespropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound which is sufficient to be explained */
   SCIP_Bool*            explanation,        /**< bool array which marks the variable which are part of the explanation if a cutoff was detected, or NULL */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   );

/** this method visualizes the cumulative structure in GML format */
EXTERN
SCIP_RETCODE SCIPvisualizeConsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< cumulative constraint */
   );

/** solves given cumulative condition as independent sub problem
 *
 *  @note The time and memory limit should be respected.
 *
 *  @note If the problem was solved to the earliest start times (ests) and latest start times (lsts) array contain the
 *        solution values; If the problem was not solved these two arrays contain the global bounds at the time the sub
 *        solver was interrupted.
 *
 *  input:
 *  - njobs           : number of jobs (activities)
 *  - objvals         : array of objective coefficients for each job (linear objective function), or NULL if none
 *  - durations       : array of durations
 *  - demands         : array of demands
 *  - capacity        : cumulative capacity
 *  - hmin            : left bound of time axis to be considered (including hmin)
 *  - hmax            : right bound of time axis to be considered (not including hmax)
 *  - timelimit       : time limit for solving in seconds
 *  - memorylimit     : memory limit for solving in mega bytes (MB)
 *  - maxnodes        : maximum number of branch-and-bound nodes to solve the single cumulative constraint  (-1: no limit)
 *
 *  input/output:
 *  - ests            : array of earliest start times for each job
 *  - lsts            : array of latest start times for each job
 *
 *  output:
 *  - solved          : pointer to store if the problem is solved (to optimality)
 *  - infeasible      : pointer to store if the problem is infeasible
 *  - unbounded       : pointer to store if the problem is unbounded
 *  - error           : pointer to store if an error occurred
 *
 */
#define SCIP_DECL_SOLVECUMULATIVE(x) SCIP_RETCODE x (int njobs, SCIP_Real* ests, SCIP_Real* lsts, SCIP_Real* objvals, \
      int* durations, int* demands, int capacity, int hmin, int hmax, \
      SCIP_Real timelimit, SCIP_Real memorylimit, SCIP_Longint maxnodes, \
      SCIP_Bool* solved, SCIP_Bool* infeasible, SCIP_Bool* unbounded, SCIP_Bool* error)

/** sets method to solve an individual cumulative condition */
EXTERN
SCIP_RETCODE SCIPsetSolveCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_SOLVECUMULATIVE((*solveCumulative)) /**< method to use an individual cumulative condition */
   );

/** solves given cumulative condition as independent sub problem
 *
 *  @note If the problem was solved to the earliest start times (ests) and latest start times (lsts) array contain the
 *        solution values; If the problem was not solved these two arrays contain the global bounds at the time the sub
 *        solver was interrupted.
 */
EXTERN
SCIP_RETCODE SCIPsolveCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   njobs,              /**< number of jobs (activities) */
   SCIP_Real*            ests,               /**< array with the earlier start time for each job */
   SCIP_Real*            lsts,               /**< array with the latest start time for each job */
   SCIP_Real*            objvals,            /**< array of objective coefficients for each job (linear objective function), or NULL if none */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   int                   hmin,               /**< left bound of time axis to be considered (including hmin) */
   int                   hmax,               /**< right bound of time axis to be considered (not including hmax) */
   SCIP_Real             timelimit,          /**< time limit for solving in seconds */
   SCIP_Real             memorylimit,        /**< memory limit for solving in mega bytes (MB) */
   SCIP_Longint          maxnodes,           /**< maximum number of branch-and-bound nodes to solve the single cumulative constraint  (-1: no limit) */
   SCIP_Bool*            solved,             /**< pointer to store if the problem is solved (to optimality) */
   SCIP_Bool*            infeasible,         /**< pointer to store if the problem is infeasible */
   SCIP_Bool*            unbounded,          /**< pointer to store if the problem is unbounded */
   SCIP_Bool*            error               /**< pointer to store if an error occurred */
   );

/** creates the worst case resource profile, that is, all jobs are inserted with the earliest start and latest
 *  completion time
 */
EXTERN
SCIP_RETCODE SCIPcreateWorstCaseProfile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROFILE*         profile,            /**< resource profile */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands             /**< array containing corresponding demands */
   );

/** computes w.r.t. the given worst case resource profile the first time point where the given capacity can be violated */
EXTERN
int SCIPcomputeHmin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROFILE*         profile,            /**< worst case resource profile */
   int                   capacity            /**< capacity to check */
   );

/** computes w.r.t. the given worst case resource profile the first time point where the given capacity is satisfied for sure */
EXTERN
int SCIPcomputeHmax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROFILE*         profile,            /**< worst case profile */
   int                   capacity            /**< capacity to check */
   );


/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
