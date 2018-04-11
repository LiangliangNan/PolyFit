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

/**@file   nlpioracle.h
 * @brief  methods to store an NLP and request function, gradient, and hessian values
 * @author Stefan Vigerske
 *
 * Not a full NLPI, but implements a common part of many NLPIs that takes care
 * of the problem storage and function, gradient, and hessian evaluation.
 */

#ifndef __SCIP_NLPI_ORACLE_H__
#define __SCIP_NLPI_ORACLE_H__

#include "scip/type_message.h"
#include "nlpi/type_nlpi.h"
#include "nlpi/type_exprinterpret.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_NlpiOracle SCIP_NLPIORACLE; /**< NLPI oracle data structure */

/** creates an NLPIORACLE data structure */
extern
SCIP_RETCODE SCIPnlpiOracleCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPIORACLE**     oracle              /**< pointer to store NLPIORACLE data structure */
   );

/** frees an NLPIORACLE data structure */
extern
SCIP_RETCODE SCIPnlpiOracleFree(
   SCIP_NLPIORACLE**     oracle              /**< pointer to NLPIORACLE data structure */
   );

/** sets the value for infinity */
extern
SCIP_RETCODE SCIPnlpiOracleSetInfinity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             infinity            /**< value to use for infinity */
   );

/** gets the value for infinity */
extern
SCIP_Real SCIPnlpiOracleGetInfinity(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** sets the problem name (used for printing) */
extern
SCIP_RETCODE SCIPnlpiOracleSetProblemName(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const char*           name                /**< name of problem */
   );

/** gets the problem name, or NULL if none set */
extern
const char* SCIPnlpiOracleGetProblemName(
   SCIP_NLPIORACLE*     oracle               /**< pointer to NLPIORACLE data structure */
   );

/** adds variables */
extern
SCIP_RETCODE SCIPnlpiOracleAddVars(
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
extern
SCIP_RETCODE SCIPnlpiOracleAddConstraints(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to add */
   const SCIP_Real*      lhss,               /**< array with left-hand sides of constraints, or NULL if all -infinity */
   const SCIP_Real*      rhss,               /**< array with right-hand sides of constraints, or NULL if all +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   const int*            nquadelems,         /**< number of elements in matrix of quadratic part for each constraint,
                                              * may be NULL in case of no quadratic part in any constraint */
   SCIP_QUADELEM* const* quadelems,          /**< quadratic elements specifying quadratic part for each constraint, entry of array may be NULL in case of no quadratic part,
                                              * may be NULL in case of no quadratic part in any constraint */
   int* const*           exprvaridxs,        /**< NULL if no nonquadratic parts, otherwise epxrvaridxs[.] maps variable indices in expression tree to indices in nlp */
   SCIP_EXPRTREE* const* exprtrees,          /**< NULL if no nonquadratic parts, otherwise exprtrees[.] gives nonquadratic part, 
                                              *   or NULL if no nonquadratic part in this constraint */
   const char**          consnames           /**< names of new constraints, or NULL if no names should be stored */
   );

/** sets or overwrites objective, a minimization problem is expected
 * 
 *  May change sparsity pattern.
 */
extern
SCIP_RETCODE SCIPnlpiOracleSetObjective(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real       constant,           /**< constant part of objective */
   int                   nlin,               /**< number of linear variable coefficients */ 
   const int*            lininds,            /**< indices of linear variables, or NULL if no linear part */
   const SCIP_Real*      linvals,            /**< coefficients of linear variables, or NULL if no linear part */
   int                   nquadelems,         /**< number of entries in matrix of quadratic part */
   const SCIP_QUADELEM*  quadelems,          /**< entries in matrix of quadratic part, may be NULL in case of no quadratic part */
   const int*            exprvaridxs,        /**< maps variable indices in expression tree to indices in nlp, or NULL if no nonquadratic part */
   const SCIP_EXPRTREE*  exprtree            /**< expression tree of nonquadratic part, or NULL if no nonquadratic part */
   );

/** change variable bounds */
extern
SCIP_RETCODE SCIPnlpiOracleChgVarBounds(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< array with indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< array with new lower bounds, or NULL if all should be -infty */
   const SCIP_Real*      ubs                 /**< array with new upper bounds, or NULL if all should be +infty */
   );

/** change constraint sides */
extern
SCIP_RETCODE SCIPnlpiOracleChgConsSides(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   nconss,             /**< number of constraints to change sides */
   const int*            indices,            /**< array with indices of constraints to change sides */
   const SCIP_Real*      lhss,               /**< array with new left-hand sides, or NULL if all should be -infty */
   const SCIP_Real*      rhss                /**< array with new right-hand sides, or NULL if all should be +infty */
   );

/** deletes a set of variables */
extern
SCIP_RETCODE SCIPnlpiOracleDelVarSet(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of vars in input (1 if var should be deleted, 0 if not); 
                                              *   new position of var in output (-1 if var was deleted) */
   );

/** deletes a set of constraints */
extern
SCIP_RETCODE SCIPnlpiOracleDelConsSet(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int*                  delstats            /**< array with deletion status of rows in input (1 if row should be deleted, 0 if not); 
                                              *   new position of row in output (-1 if row was deleted) */
   );

/** changes (or adds) linear coefficients in one constraint or objective */
extern
SCIP_RETCODE SCIPnlpiOracleChgLinearCoefs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where linear coefficients should be changed, or -1 for objective */
   int                   nentries,           /**< number of coefficients to change */
   const int*            varidxs,            /**< array with indices of variables which coefficients should be changed */
   const SCIP_Real*      newcoefs            /**< array with new coefficients of variables */
   );

/** changes (or adds) coefficients in the quadratic part of one constraint or objective */
extern
SCIP_RETCODE SCIPnlpiOracleChgQuadCoefs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where quadratic coefficients should be changed, or -1 for objective */
   int                   nquadelems,         /**< number of entries in quadratic constraint to change */
   const SCIP_QUADELEM*  quadelems           /**< new elements in quadratic matrix (replacing already existing ones or adding new ones) */
   );

/** replaces expression tree of one constraint or objective */
extern
SCIP_RETCODE SCIPnlpiOracleChgExprtree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where expression tree should be changed, or -1 for objective */
   const int*            exprvaridxs,        /**< problem indices of variables in expression tree */
   const SCIP_EXPRTREE*  exprtree            /**< new expression tree, or NULL */
   );

/** changes one parameter of expression tree of one constraint or objective
 */
extern
SCIP_RETCODE SCIPnlpiOracleChgExprParam(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint where parameter should be changed in expression tree, or -1 for objective */
   int                   paramidx,           /**< index of parameter */
   SCIP_Real             paramval            /**< new value of parameter */
   );

/** changes the constant value in the objective function
 */
extern
SCIP_RETCODE SCIPnlpiOracleChgObjConstant(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real             objconstant         /**< new value for objective constant */
   );

/** gives the current number of variables */
extern
int SCIPnlpiOracleGetNVars(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the current number of constraints */
extern
int SCIPnlpiOracleGetNConstraints(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables lower bounds */
extern
const SCIP_Real* SCIPnlpiOracleGetVarLbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables upper bounds */
extern
const SCIP_Real* SCIPnlpiOracleGetVarUbs(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives the variables names, or NULL if not set */
extern
char** SCIPnlpiOracleGetVarNames(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** Gives maximum degree of a variable w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
extern
int SCIPnlpiOracleGetVarDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   varidx              /**< the variable for which the degree is returned */
   );

/** Gives maximum degree of all variables w.r.t. objective and all constraints.
 *  The degree of a variable is the degree of the summand where it appears in, and is infinity for nonpolynomial terms.
 */ 
extern
int* SCIPnlpiOracleGetVarDegrees(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** gives left-hand side of a constraint */
extern
SCIP_Real SCIPnlpiOracleGetConstraintLhs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   );

/** gives right-hand side of a constraint */
extern
SCIP_Real SCIPnlpiOracleGetConstraintRhs(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   );

/** gives name of a constraint, may be NULL */
extern
char* SCIPnlpiOracleGetConstraintName(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< constraint index */
   );

/** gives maximum degree of a constraint or objective
 *  The degree is the maximal degree of all summands,, and is infinity for nonpolynomial terms.
 */ 
extern
int SCIPnlpiOracleGetConstraintDegree(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx             /**< index of constraint for which the degree is requested, or -1 for objective */
   );

/** Gives maximum degree over all constraints and the objective (or over all variables, resp.).
 * Thus, if this function returns 0, then the objective and all constraints are constant.
 * If it returns 1, then the problem in linear.
 * If it returns 2, then its a QP, QCP, or QCQP.
 * And if it returns > 2, then it is an NLP.
 */
extern
int SCIPnlpiOracleGetMaxDegree(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** Gives the evaluation capabilities that are shared among all expression trees in the problem. */
extern
SCIP_EXPRINTCAPABILITY SCIPnlpiOracleGetEvalCapability(
   SCIP_NLPIORACLE*      oracle              /**< pointer to NLPIORACLE data structure */
   );

/** evaluates the objective function in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            objval              /**< pointer to store objective value */  
   );

/** evaluates one constraint function in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValue(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   int                   considx,            /**< index of constraint to evaluate */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            conval              /**< pointer to store constraint value */  
   );

/** evaluates all constraint functions in a given point */
extern
SCIP_RETCODE SCIPnlpiOracleEvalConstraintValues(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Real*            convals             /**< pointer to store constraint values */  
   );

/** computes the objective gradient in a given point
 *
 * @return SCIP_INVALIDDATA, if the function or its gradient could not be evaluated (domain error, etc.)
 */
extern
SCIP_RETCODE SCIPnlpiOracleEvalObjectiveGradient(
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
extern
SCIP_RETCODE SCIPnlpiOracleEvalConstraintGradient(
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
extern
SCIP_RETCODE SCIPnlpiOracleGetJacobianSparsity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, 
                                              *   offsets[nconss] gives length of col, can be NULL */
   );

/** evaluates the Jacobi matrix in a given point
 * 
 *  The values in the Jacobi matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetJacobianSparsity.
 *  The user need to call SCIPnlpiOracleGetJacobianSparsity at least ones before using this function.
 *
 *  @return SCIP_INVALIDDATA, if the Jacobian could not be evaluated (domain error, etc.)
 */
extern
SCIP_RETCODE SCIPnlpiOracleEvalJacobian(
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
extern
SCIP_RETCODE SCIPnlpiOracleGetHessianLagSparsity(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const int**           offset,             /**< pointer to store pointer that stores the offsets to each rows sparsity pattern in col, can be NULL */
   const int**           col                 /**< pointer to store pointer that stores the indices of variables that appear in each row, 
                                              *   offsets[nconss] gives length of col, can be NULL */
   );

/** evaluates the Hessian matrix of the Lagrangian in a given point
 * 
 *  The values in the Hessian matrix are returned in the same order as specified by the offset and col arrays obtained by SCIPnlpiOracleGetHessianLagSparsity.
 *  The user must call SCIPnlpiOracleGetHessianLagSparsity at least ones before using this function. 
 *  Only elements of the lower left triangle and the diagonal are computed.
 *
 * @return SCIP_INVALIDDATA, if the Hessian could not be evaluated (domain error, etc.)
 */
extern
SCIP_RETCODE SCIPnlpiOracleEvalHessianLag(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   const SCIP_Real*      x,                  /**< point where to evaluate */
   SCIP_Bool             isnewx,             /**< has the point x changed since the last call to some evaluation function? */
   SCIP_Real             objfactor,          /**< weight for objective function */
   const SCIP_Real*      lambdas,            /**< array with weights (Lagrangian multipliers) for the constraints */ 
   SCIP_Real*            hessian             /**< pointer to store sparse hessian values */  
   );

/** prints the problem to a file. */
extern
SCIP_RETCODE SCIPnlpiOraclePrintProblem(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   );

/** prints the problem to a file in GAMS format
 * If there are variable (equation, resp.) names with more than 9 characters, then variable (equation, resp.) names are prefixed with an unique identifier.
 * This is to make it easier to identify variables solution output in the listing file.
 * Names with more than 64 characters are shorten to 64 letters due to GAMS limits.
 */
extern
SCIP_RETCODE SCIPnlpiOraclePrintProblemGams(
   SCIP_NLPIORACLE*      oracle,             /**< pointer to NLPIORACLE data structure */
   SCIP_Real*            initval,            /**< starting point values for variables or NULL */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to print to, or NULL for standard output */
   );

#ifdef __cplusplus
}
#endif

#endif
