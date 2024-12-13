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

/**@file   nlpioracle.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods to store an NLP and request function, gradient, and Hessian values
 * @author Stefan Vigerske
 *
 * A common part of many NLPIs that takes care of the problem storage and function, gradient, and Hessian evaluation.
 */

#ifndef __SCIP_NLPIORACLE_H__
#define __SCIP_NLPIORACLE_H__

#include "scip/type_message.h"
#include "scip/type_exprinterpret.h"

/**@defgroup NLPIOracle NLPI Oracle
 * @ingroup DataStructures
 * @brief NLP representation used by various NLP solver interface implementations
 *
 *@{
 */


#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_NlpiOracle SCIP_NLPIORACLE; /**< NLPI oracle data structure */

/** creates an NLPIORACLE data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
   );

/** frees an NLPIORACLE data structure */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
   );

/** sets the problem name (used for printing) */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleSetProblemName(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const char*           name                /**< name of problem */
   );

/** gets the problem name, or NULL if none set */
SCIP_EXPORT
const char* SCIPnlpiOracleGetProblemName(
   SCIP_NLPIORACLE*     oracle               /**< pointer to NLPIORACLE data structure */
   );

/** adds variables */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleAddVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to add */
   const SCIP_Real*      lbs,                /**< array with lower bounds of new variables, or NULL if all -infinity */
   const SCIP_Real*      ubs,                /**< array with upper bounds of new variables, or NULL if all +infinity */
   const char**          varnames            /**< array with names of new variables, or NULL if no names should be stored */
   );

/** adds constraints
 *
 *  linear coefficients: row(=constraint) oriented matrix;
 *  quadratic coefficients: row oriented matrix for each constraint
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to add */
   const SCIP_Real*      lhss,               /**< array with left-hand sides of constraints, or NULL if all -infinity */
   const SCIP_Real*      rhss,               /**< array with right-hand sides of constraints, or NULL if all +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   SCIP_EXPR**           exprs,              /**< NULL if no nonlinear parts, otherwise exprs[.] gives nonlinear part,
                                              *   or NULL if no nonlinear part in this constraint */
   const char**          consnames           /**< names of new constraints, or NULL if no names should be stored */
   );

/** sets or overwrites objective, a minimization problem is expected
 *
 *  May change sparsity pattern.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleSetObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real       constant,           /**< constant part of objective */
   int                   nlin,               /**< number of linear variable coefficients */
   const int*            lininds,            /**< indices of linear variables, or NULL if no linear part */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables, or NULL if no linear part */
   SCIP_EXPR*            expr                /**< expression of nonlinear part, or NULL if no nonlinear part */
   );

/** change variable bounds */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< array with indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< array with new lower bounds, or NULL if all should be -infty */
   const SCIP_Real*      ubs                 /**< array with new upper bounds, or NULL if all should be +infty */
   );

/** change constraint sides */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleChgConsSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to change sides */
   const int*            indices,            /**< array with indices of constraints to change sides */
   const SCIP_Real*      lhss,               /**< array with new left-hand sides, or NULL if all should be -infty */
   const SCIP_Real*      rhss                /**< array with new right-hand sides, or NULL if all should be +infty */
   );

/** deletes a set of variables */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleDelVarSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of vars in input (1 if var should be deleted, 0 if not);
                                              *   new position of var in output (-1 if var was deleted) */
   );

/** deletes a set of constraints */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of rows in input (1 if row should be deleted, 0 if not);
                                              *   new position of row in output (-1 if row was deleted) */
   );

/** changes (or adds) linear coefficients in one constraint or objective */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where linear coefficients should be changed, or -1 for objective */
   int                   nentries,           /**< number of coefficients to change */
   const int*            varidxs,            /**< array with indices of variables which coefficients should be changed */
   const SCIP_Real*      newcoefs            /**< array with new coefficients of variables */
   );

/** replaces expression of one constraint or objective */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleChgExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where expression should be changed, or -1 for objective */
   SCIP_EXPR*            expr                /**< new expression, or NULL */
   );

/** changes the constant value in the objective function
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleChgObjConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             objconstant         /**< new value for objective constant */
   );

/** gives the current number of variables */
SCIP_EXPORT
int SCIPnlpiOracleGetNVars(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the current number of constraints */
SCIP_EXPORT
int SCIPnlpiOracleGetNConstraints(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables lower bounds */
SCIP_EXPORT
const SCIP_Real* SCIPnlpiOracleGetVarLbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables upper bounds */
SCIP_EXPORT
const SCIP_Real* SCIPnlpiOracleGetVarUbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables names, or NULL if not set */
SCIP_EXPORT
char** SCIPnlpiOracleGetVarNames(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** indicates whether variable appears nonlinear in any objective or constraint */
SCIP_EXPORT
SCIP_Bool SCIPnlpiOracleIsVarNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   varidx              /**< the variable to check */
   );

/** returns number of linear and nonlinear appearances of variables in objective and constraints */
SCIP_EXPORT
void SCIPnlpiOracleGetVarCounts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           lincounts,          /**< buffer to return pointer to array of counts of linear appearances */
   const int**           nlcounts            /**< buffer to return pointer to array of counts of nonlinear appearances */
   );

/** gives constant term of objective */
SCIP_EXPORT
SCIP_Real SCIPnlpiOracleGetObjectiveConstant(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives left-hand side of a constraint */
SCIP_EXPORT
SCIP_Real SCIPnlpiOracleGetConstraintLhs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   );

/** gives right-hand side of a constraint */
SCIP_EXPORT
SCIP_Real SCIPnlpiOracleGetConstraintRhs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   );

/** gives name of a constraint, may be NULL */
SCIP_EXPORT
char* SCIPnlpiOracleGetConstraintName(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   );

/** indicates whether constraint is nonlinear */
SCIP_EXPORT
SCIP_Bool SCIPnlpiOracleIsConstraintNonlinear(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< index of constraint for which nonlinearity status is returned, or -1 for objective */
   );

/** gives the evaluation capabilities that are shared among all expressions in the problem */
SCIP_EXPORT
SCIP_EXPRINTCAPABILITY SCIPnlpiOracleGetEvalCapability(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** evaluates the objective function in a given point */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            objval              /**< pointer to store objective value */
   );

/** evaluates one constraint function in a given point */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to store constraint value */
   );

/** evaluates all constraint functions in a given point */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            convals             /**< pointer to store constraint values */
   );

/** computes the objective gradient in a given point
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            objval,             /**< pointer to buffer objective value */
   SCIP_Real*            objgrad             /**< pointer to buffer (dense) objective gradient */
   );

/** computes a constraints gradient in a given point
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int             considx,            /**< index of constraint to compute gradient for */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            conval,             /**< pointer to store constraint value */
   SCIP_Real*            congrad             /**< pointer to store (dense) constraint gradient */
   );

/** gets sparsity pattern (rowwise) of Jacobian matrix
 *
 *  Note that internal data is returned in *offset and *col, thus the user does not need to allocate memory there.
 *  Adding or deleting constraints destroys the sparsity structure and make another call to this function necessary.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleGetJacobianSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row,
                                              *   offsets[nconss] gives length of col, can be NULL */
   );

/** evaluates the Jacobian matrix in a given point
 *
 *  The values in the Jacobian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetJacobianSparsity().
 *  The user need to call SCIPnlpiOracleGetJacobianSparsity() at least ones before using this function.
 *
 * @return SCIP_INVALIDDATA, if the Jacobian could not be evaluated (domain error, etc.)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleEvalJacobian(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real*            convals,            /**< pointer to store constraint values, can be NULL */
   SCIP_Real*            jacobi              /**< pointer to store sparse jacobian values */
   );

/** gets sparsity pattern of the Hessian matrix of the Lagrangian
 *
 *  Note that internal data is returned in *offset and *col, thus the user must not to allocate memory there.
 *  Adding or deleting variables, objective, or constraints may destroy the sparsity structure and make another call to this function necessary.
 *  Only elements of the lower left triangle and the diagonal are counted.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row,
                                              *   offsets[nconss] gives length of col, can be NULL */
   );

/** evaluates the Hessian matrix of the Lagrangian in a given point
 *
 *  The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity().
 *  The user must call SCIPnlpiOracleGetHessianLagSparsity() at least ones before using this function.
 *  Only elements of the lower left triangle and the diagonal are computed.
 *
 * @return SCIP_INVALIDDATA, if the Hessian could not be evaluated (domain error, etc.)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx_obj,         /**< has the point x changed since the last call to an objective evaluation function? */
   SCIP_Bool             isnewx_cons,        /**< has the point x changed since the last call to the constraint evaluation function? */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambdas,            /**< array with weights (Lagrangian multipliers) for the constraints */
   SCIP_Real*            hessian             /**< pointer to store sparse hessian values */
   );

/** resets clock that measures evaluation time */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOracleResetEvalTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives time spend in evaluation since last reset of clock
 *
 * Gives 0 if the eval clock is disabled.
 */
SCIP_EXPORT
SCIP_Real SCIPnlpiOracleGetEvalTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** prints the problem to a file. */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOraclePrintProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   );

/** prints the problem to a file in GAMS format
 *
 * If there are variable (equation, resp.) names with more than 9 characters, then variable (equation, resp.) names are prefixed with an unique identifier.
 * This is to make it easier to identify variables solution output in the listing file.
 * Names with more than 64 characters are shorten to 64 letters due to GAMS limits.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
