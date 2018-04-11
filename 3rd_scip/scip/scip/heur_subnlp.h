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

/**@file   heur_subnlp.h
 * @ingroup PRIMALHEURISTICS
 * @brief  NLP local search primal heuristic using sub-SCIPs
 * @author Stefan Vigerske
 *
 * This heuristic applies a NLP local search to a nonlinear CIP after fixing all discrete variables.
 * That is, the CIP is copied, all discrete variables are fixed, presolving is applied,
 * and if the resulting CIP has a nonlinear relaxation, then it is tried to solve this relaxation
 * by an NLP solver.
 * The heuristic only runs if continuous nonlinearities are present (@ref SCIPhasNLPContinuousNonlinearity()).
 *
 * Fixing values for discrete values are either taken from a solution of the LP relaxation which
 * satisfies all integrality constraints, or are provided by SCIPupdateStartpointHeurSubNlp().
 *
 * This heuristic is orthogonal to the undercover heuristic (@ref heur_undercover.h), which fixes
 * variables in a nonlinear CIP in a way that a (possibly mixed-integer) linear subproblem is obtained.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef HEUR_SUBNLP_H_
#define HEUR_SUBNLP_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the NLP local search primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurSubNlp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
  *
  * @{
  */

/** updates the starting point for the NLP heuristic
 * 
 * Is called, for example, by a constraint handler that handles nonlinear constraints when a check on feasibility of a solution fails.
 */
EXTERN
SCIP_RETCODE SCIPupdateStartpointHeurSubNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< subNLP heuristic */
   SCIP_SOL*             solcand,            /**< solution candidate */
   SCIP_Real             violation           /**< constraint violation of solution candidate */
   );

/** main procedure of the subNLP heuristic */
EXTERN
SCIP_RETCODE SCIPapplyHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< pointer to store result of: solution found, no solution found, or fixing is infeasible (cutoff) */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver                                 */
   SCIP_Real             timelimit,          /**< time limit for NLP solver                                      */
   SCIP_Real             minimprove,         /**< desired minimal relative improvement in objective function value */
   SCIP_Longint*         iterused,           /**< buffer to store number of iterations used by NLP solver, or NULL if not of interest */
   SCIP_SOL*             resultsol           /**< a solution where to store found solution values, if any, or NULL if to try adding to SCIP */
   );

/** for a given solution, resolves the corresponding subNLP and updates solution values for continuous variables, if NLP solution is feasible in original problem */
EXTERN
SCIP_RETCODE SCIPresolveSolHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL*             sol,                /**< solution for which to solve NLP, and where to store resolved solution values */
   SCIP_Bool*            success,            /**< buffer where to store whether a feasible solution was found */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver, or -1 for default of NLP heuristic */
   SCIP_Real             timelimit           /**< time limit for NLP solver */
   );

/** adds all known linear constraint to the NLP, if initialized and not done already
 * This function is temporary and will hopefully become obsolete in the near future.
 */ 
EXTERN
SCIP_RETCODE SCIPaddLinearConsToNlpHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints, i.e., linear constraints that involve only discrete variables */
   SCIP_Bool             addcontconss        /**< whether to add continuous    linear constraints, i.e., linear constraints that involve not only discrete variables */
   );

/** gets sub-SCIP used by NLP heuristic, or NULL if none */
EXTERN
SCIP* SCIPgetSubScipHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   );

/** gets mapping of SCIP variables to sub-SCIP variables */
EXTERN
SCIP_VAR** SCIPgetVarMappingScip2SubScipHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   );

/** gets mapping of sub-SCIP variables to SCIP variables */
EXTERN
SCIP_VAR** SCIPgetVarMappingSubScip2ScipHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   );

/** gets startpoint candidate to be used in next call to NLP heuristic, or NULL if none */
EXTERN
SCIP_SOL* SCIPgetStartCandidateHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif /*HEUR_SUBNLP_H_*/
