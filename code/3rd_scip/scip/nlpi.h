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

/**@file   nlpi.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for NLP solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_H__
#define __SCIP_NLPI_H__

#include "scip/type_nlpi.h"
#include "scip/type_misc.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates an NLP solver interface */
SCIP_RETCODE SCIPnlpiCreate(
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying of NLPI, can be NULL */
   SCIP_DECL_NLPIFREE              ((*nlpifree)),               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer)),   /**< get solver pointer, can be NULL */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem)),      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem)),        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer)),  /**< get problem pointer, can be NULL */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars)),            /**< add variables */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints)),     /**< add constraints */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective)),       /**< set objective */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds)),       /**< change variable bounds */
   SCIP_DECL_NLPICHGCONSSIDES      ((*nlpichgconssides)),       /**< change constraint sides */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset)),          /**< delete a set of constraints */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset)),         /**< delete a set of constraints */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs)),     /**< change coefficients in linear part of a constraint or objective */
   SCIP_DECL_NLPICHGEXPR           ((*nlpichgexpr)),            /**< change nonlinear expression a constraint or objective */
   SCIP_DECL_NLPICHGOBJCONSTANT    ((*nlpichgobjconstant)),     /**< change the constant offset in the objective */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess)),    /**< set initial guess, can be NULL */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve)),              /**< solve NLP */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat)),         /**< get solution status */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat)),        /**< get termination status */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution)),        /**< get solution */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics)),      /**< get solve statistics */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
   );

/** sets NLP solver priority */
void SCIPnlpiSetPriority(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   int                   priority            /**< new priority of NLPI */
   );

/** copies an NLPI and includes it into another SCIP instance */
SCIP_RETCODE SCIPnlpiCopyInclude(
   SCIP_NLPI*            sourcenlpi,         /**< the NLP interface to copy */
   SCIP_SET*             targetset           /**< global SCIP settings where to include copy */
   );

/** frees NLPI */
SCIP_RETCODE SCIPnlpiFree(
   SCIP_NLPI**           nlpi,               /**< pointer to NLPI data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes NLPI */
void SCIPnlpiInit(
   SCIP_NLPI*            nlpi                /**< solver interface */
   );

/** gets pointer for NLP solver */
void* SCIPnlpiGetSolverPointer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance, or NULL */
   );

/** creates a problem instance */
SCIP_RETCODE SCIPnlpiCreateProblem(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM**    problem,            /**< problem pointer to store the problem data */
   const char*           name                /**< name of problem, can be NULL */
   );

/** frees a problem instance */
SCIP_RETCODE SCIPnlpiFreeProblem(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM**    problem             /**< pointer where problem instance is stored */
   );

/** gets pointer to solver-internal problem instance */
void* SCIPnlpiGetProblemPointer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   );

/** add variables to nlpi */
SCIP_RETCODE SCIPnlpiAddVars(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      lbs,                /**< lower bounds of variables, can be NULL if -infinity */
   const SCIP_Real*      ubs,                /**< upper bounds of variables, can be NULL if +infinity */
   const char**          varnames            /**< names of variables, can be NULL */
   );

/** add constraints to nlpi */
SCIP_RETCODE SCIPnlpiAddConstraints(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nconss,             /**< number of constraints */
   const SCIP_Real*      lhss,               /**< left hand sides of constraints, can be NULL if -infinity */
   const SCIP_Real*      rhss,               /**< right hand sides of constraints, can be NULL if +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   SCIP_EXPR**           exprs,              /**< expressions for nonlinear part of constraints, entry of array may be NULL in case of no nonlinear part, may be NULL in case of no nonlinear part in any constraint */
   const char**          names               /**< names of constraints, may be NULL or entries may be NULL */
   );

/** sets or overwrites objective, a minimization problem is expected */
SCIP_RETCODE SCIPnlpiSetObjective(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nlins,              /**< number of linear variables */
   const int*            lininds,            /**< variable indices, may be NULL in case of no linear part */
   const SCIP_Real*      linvals,            /**< coefficient values, may be NULL in case of no linear part */
   SCIP_EXPR*            expr,               /**< expression for nonlinear part of objective function, may be NULL in case of no nonlinear part */
   const SCIP_Real       constant            /**< objective value offset */
   );

/** change variable bounds */
SCIP_RETCODE SCIPnlpiChgVarBounds(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   const int             nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   );

/** change constraint sides */
SCIP_RETCODE SCIPnlpiChgConsSides(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   nconss,             /**< number of constraints to change sides */
   const int*            indices,            /**< indices of constraints to change sides */
   const SCIP_Real*      lhss,               /**< new left hand sides */
   const SCIP_Real*      rhss                /**< new right hand sides */
   );

/** delete a set of variables */
SCIP_RETCODE SCIPnlpiDelVarSet(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int*                  dstats,             /**< deletion status of vars; 1 if var should be deleted, 0 if not */
   int                   dstatssize          /**< size of the dstats array */
   );

/** delete a set of constraints */
SCIP_RETCODE SCIPnlpiDelConsSet(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int*                  dstats,             /**< deletion status of constraints; 1 if constraint should be deleted, 0 if not */
   int                   dstatssize          /**< size of the dstats array */
   );

/** changes or adds linear coefficients in a constraint or objective */
SCIP_RETCODE SCIPnlpiChgLinearCoefs(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   idx,                /**< index of constraint or -1 for objective */
   int                   nvals,              /**< number of values in linear constraint to change */
   const int*            varidxs,            /**< indices of variables which coefficient to change */
   const SCIP_Real*      vals                /**< new values for coefficients */
   );

/** change the expression in the nonlinear part */
SCIP_RETCODE SCIPnlpiChgExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   int                   idxcons,            /**< index of constraint or -1 for objective */
   SCIP_EXPR*            expr                /**< new expression for constraint or objective, or NULL to only remove previous tree */
   );

/** change the constant offset in the objective */
SCIP_RETCODE SCIPnlpiChgObjConstant(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real             objconstant         /**< new value for objective constant */
   );

/** sets initial guess for primal variables */
SCIP_RETCODE SCIPnlpiSetInitialGuess(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real*            primalvalues,       /**< initial primal values for variables, or NULL to clear previous values */
   SCIP_Real*            consdualvalues,     /**< initial dual values for constraints, or NULL to clear previous values */
   SCIP_Real*            varlbdualvalues,    /**< initial dual values for variable lower bounds, or NULL to clear previous values */
   SCIP_Real*            varubdualvalues     /**< initial dual values for variable upper bounds, or NULL to clear previous values */
   );

/** tries to solve NLP */
SCIP_RETCODE SCIPnlpiSolve(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPPARAM*        param               /**< solve parameters */
   );

/** gives solution status */
SCIP_NLPSOLSTAT SCIPnlpiGetSolstat(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   );

/** gives termination reason */
SCIP_NLPTERMSTAT SCIPnlpiGetTermstat(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem             /**< problem instance */
   );

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_RETCODE SCIPnlpiGetSolution(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_Real**           primalvalues,       /**< buffer to store pointer to array to primal values, or NULL if not needed */
   SCIP_Real**           consdualvalues,     /**< buffer to store pointer to array to dual values of constraints, or NULL if not needed */
   SCIP_Real**           varlbdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real**           varubdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real*            objval              /**< pointer to store the objective value, or NULL if not needed */
   );

/** gives solve statistics */
SCIP_RETCODE SCIPnlpiGetStatistics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< solver interface */
   SCIP_NLPIPROBLEM*     problem,            /**< problem instance */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_H__ */
