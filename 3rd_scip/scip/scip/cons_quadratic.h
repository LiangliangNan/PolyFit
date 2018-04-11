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

/**@file   cons_quadratic.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for quadratic constraints \f$\textrm{lhs} \leq \sum_{i,j=1}^n a_{i,j} x_ix_j + \sum_{i=1}^n b_i x_i \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_QUADRATIC_H__
#define __SCIP_CONS_QUADRATIC_H__

#include "scip/scip.h"
#include "scip/intervalarith.h"
#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the handler for quadratic constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Quadratic Constraints
 *
 * @{
 *
 * This constraint handler handles constraints of the form
 * \f[
 *   \textrm{lhs} \leq \sum_{i,j=1}^n a_{i,j} x_ix_j + \sum_{i=1}^n b_i x_i \leq \textrm{rhs}
 * \f]
 *
 * Constraints are enforced by separation, domain propagation, and spatial branching.
 *
 * For semidefinite matrices \f$A=(a_{i,j})_{i,j}\f$, cuts based on linearization of \f$\langle x, Ax\rangle\f$ are implemented.
 * For underestimating a non-convex term, McCormick underestimators and secants for univariate concave quadratic terms are implemented.
 * If \f$\langle x, Ax\rangle\f$ is factorable (i.e., can be written as product of two linear functions),
 * specialized separation techniques (e.g., lifted tangent inequalities) that take the constraint sides into account are applied.
 *
 * Branching is performed for variables in nonconvex terms, if the relaxation solution cannot be separated.
 * Further, domain propagation is applied.
 *
 * During presolve, variable products which contain binary variables may be reformulated into linear constraints, thereby introducing new variables.
 *
 * See also
 * @par
 * Timo Berthold and Stefan Heinz and Stefan Vigerske@n
 * <a href="http://dx.doi.org/10.1007/978-1-4614-1927-3">Extending a CIP framework to solve MIQCPs</a>@n
 * In: Jon Lee and Sven Leyffer (eds.),
 *     Mixed-integer nonlinear optimization: Algorithmic advances and applications,
 *     IMA volumes in Mathematics and its Applications, volume 154, 427-444, 2012.
 *
 * @par
 * Stefan Vigerske@n
 * Decomposition of Multistage Stochastic Programs and a Constraint Integer Programming Approach to Mixed-Integer Nonlinear Programming@n
 * PhD Thesis, Humboldt-University Berlin, 2012, submitted.
 *
 * @par
 * Pietro Belotti and Andrew J. Miller and Mahdi Namazifar@n
 * Linear inequalities for bounded products of variables@n
 * SIAG/OPT Views-and-News 22:1, 1-8, 2011.
 */

/** event data for variable bound changes in quadratic constraints */
typedef struct SCIP_QuadVarEventData SCIP_QUADVAREVENTDATA;

/** data structure to store a single term associated to a quadratic variable
 */
struct SCIP_QuadVarTerm
{
   SCIP_VAR*             var;                /**< quadratic variable */
   SCIP_Real             lincoef;            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef;            /**< square coefficient of variable */

   int                   nadjbilin;          /**< number of bilinear terms this variable is involved in */
   int                   adjbilinsize;       /**< size of adjacent bilinear terms array */
   int*                  adjbilin;           /**< indices of associated bilinear terms */

   SCIP_QUADVAREVENTDATA* eventdata;          /**< event data for bound change events */
};
typedef struct SCIP_QuadVarTerm SCIP_QUADVARTERM;

/** data structure to store a single bilinear term (similar to SCIP_QUADELEM)
 * except for temporary reasons, we assume that the index of var1 is smaller than the index of var2
 */
struct SCIP_BilinTerm
{
   SCIP_VAR*             var1;
   SCIP_VAR*             var2;
   SCIP_Real             coef;
};
typedef struct SCIP_BilinTerm SCIP_BILINTERM;

/** storage for a linear row in preparation
 *
 * Uses to assemble data that could eventually make a SCIP_ROW.
 * @note Only one-sided rows are allowed here.
 */
struct SCIP_RowPrep
{
   SCIP_VAR**            vars;               /**< variables */
   SCIP_Real*            coefs;              /**< coefficients of variables */
   int                   nvars;              /**< number of variables (= number of coefficients) */
   int                   varssize;           /**< length of variables array (= lengths of coefficients array) */
   SCIP_Real             side;               /**< side */
   SCIP_SIDETYPE         sidetype;           /**< type of side */
   SCIP_Bool             local;              /**< whether the row is only locally valid (i.e., for the current node) */
   char                  name[SCIP_MAXSTRLEN]; /**< row name */
};
typedef struct SCIP_RowPrep SCIP_ROWPREP;

/** upgrading method for quadratic constraints into more specific constraints
 *
 * the method might upgrade a quadratic constraint into a set of quadratic constraints
 * the caller provided an array upgdconss to store upgrade constraints
 * the length of upgdconss is given by upgdconsssize
 * if an upgrade is not possible, set *nupgdconss to zero
 * if more than upgdconsssize many constraints shall replace cons, the function
 * should return the required number as negated value in *nupgdconss
 * i.e., if cons should be replaced by 3 constraints, the function should set
 * *nupgdconss to -3 and return with SCIP_OKAY
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cons            : the quadratic constraint to upgrade
 *  - nbinlin         : number of binary variables in linear part
 *  - nbinquad        : number of binary variables in quadratic part
 *  - nintlin         : number of integer variables in linear part
 *  - nintquad        : number of integer variables in quadratic part
 *  - nimpllin        : number of implicit integer variables in linear part
 *  - nimplquad       : number of implicit integer variables in quadratic part
 *  - ncontlin        : number of continuous variables in linear part
 *  - ncontquad       : number of continuous variables in quadratic part
 *  - integral        : TRUE iff constraints activity value is always integral
 *  - nupgdconss      : pointer to store number of constraints that replace this constraint
 *  - upgdconss       : array to store constraints that replace this constraint
 *  - upgdconsssize   : length of the provided upgdconss array
 *  - presoltiming    : current presolve timing
 */
#define SCIP_DECL_QUADCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, \
      int nbinlin, int nbinquad, int nintlin, int nintquad, int nimpllin, int nimplquad, int ncontlin, int ncontquad, \
      SCIP_Bool integral, int* nupgdconss, SCIP_CONS** upgdconss, int upgdconsssize, SCIP_PRESOLTIMING presoltiming)

/** includes a quadratic constraint upgrade method into the quadratic constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeQuadconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd)),  /**< method to call for upgrading quadratic constraint */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method be active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   );

/** Creates and captures a quadratic constraint.
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_j z_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< variables in linear part (x_i) or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< coefficients of variables in linear part (b_i) or NULL if nlinvars == 0 */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) or NULL if nquadterms == 0 */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) or NULL if nquadterms == 0 */
   SCIP_Real*            quadcoeffs,         /**< array with coefficients of quadratic terms (a_j) or NULL if nquadterms == 0 */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (l) */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation (u) */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** creates and captures a quadratic constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @see SCIPcreateConsQuadratic() for the default constraint flag configuration
 *

 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs                 /**< right hand side of quadratic equation (u) */
   );

/** creates and captures a quadratic constraint.
 *
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_k v_k w_k \leq u.
 * \f]
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadvarterms,      /**< number of quadratic terms (m) */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms (p) */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             lhs,                /**< constraint left hand side (ell) */
   SCIP_Real             rhs,                /**< constraint right hand side (u) */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

/** creates and captures a quadratic constraint
 *  in its most basic variant, i. e., with all constraint flags set to their default values, which can be set
 *  afterwards using SCIPsetConsFLAGNAME() in scip.h
 *
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_kv_kw_k \leq u.
 * \f]
 *
 *  @see SCIPcreateConsQuadratic2() for the default constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
EXTERN
SCIP_RETCODE SCIPcreateConsBasicQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadvarterms,      /**< number of quadratic terms (m) */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms (p) */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             lhs,                /**< constraint left hand side (ell) */
   SCIP_Real             rhs                 /**< constraint right hand side (u) */
   );

/** Adds a constant to the constraint function, that is, subtracts a constant from both sides */
EXTERN
void SCIPaddConstantQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             constant            /**< constant to subtract from both sides */
   );

/** Adds a linear variable with coefficient to a quadratic constraint.
 */
EXTERN
SCIP_RETCODE SCIPaddLinearVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient of variable */
   );

/** Adds a quadratic variable with linear and square coefficient to a quadratic constraint.
 */
EXTERN
SCIP_RETCODE SCIPaddQuadVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             lincoef,            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef             /**< square coefficient of variable */
   );

/** Adds a linear coefficient for a quadratic variable.
 *
 * Variable will be added with square coefficient 0.0 if not existing yet.
 */
EXTERN
SCIP_RETCODE SCIPaddQuadVarLinearCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to linear coefficient of variable */
   );

/** Adds a square coefficient for a quadratic variable.
 *
 * Variable will be added with linear coefficient 0.0 if not existing yet.
 */
EXTERN
SCIP_RETCODE SCIPaddSquareCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to square coefficient of variable */
   );

/** Adds a bilinear term to a quadratic constraint.
 *
 * Variables will be added with linear and square coefficient 0.0 if not existing yet.
 * If variables are equal, only the square coefficient of the variable is updated.
 */
EXTERN
SCIP_RETCODE SCIPaddBilinTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   );

/** Gets the quadratic constraint as a nonlinear row representation.
 */
EXTERN
SCIP_RETCODE SCIPgetNlRowQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   );

/** Gets the number of variables in the linear part of a quadratic constraint.
 */
EXTERN
int SCIPgetNLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the variables in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
EXTERN
SCIP_VAR** SCIPgetLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the coefficients in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
EXTERN
SCIP_Real* SCIPgetCoefsLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the number of quadratic variable terms of a quadratic constraint.
 */
EXTERN
int SCIPgetNQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the quadratic variable terms of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarTermsQuadratic.
 */
EXTERN
SCIP_QUADVARTERM* SCIPgetQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Ensures that quadratic variable terms are sorted. */
EXTERN
SCIP_RETCODE SCIPsortQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Finds the position of a quadratic variable term for a given variable.
 *
 * @note If the quadratic variable terms have not been sorted before, then a search may reorder the current order of the terms.
 */
EXTERN
SCIP_RETCODE SCIPfindQuadVarTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to search for */
   int*                  pos                 /**< buffer to store position of quadvarterm for var, or -1 if not found */
   );

/** Gets the number of bilinear terms of a quadratic constraint.
 */
EXTERN
int SCIPgetNBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the bilinear terms of a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
EXTERN
SCIP_BILINTERM* SCIPgetBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the left hand side of a quadratic constraint.
 */
EXTERN
SCIP_Real SCIPgetLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the right hand side of a quadratic constraint.
 */
EXTERN
SCIP_Real SCIPgetRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** get index of a variable in linvars that may be decreased without making any other constraint infeasible, or -1 if none */
EXTERN
int SCIPgetLinvarMayDecreaseQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** get index of a variable in linvars that may be increased without making any other constraint infeasible, or -1 if none */
EXTERN
int SCIPgetLinvarMayIncreaseQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Check the quadratic function of a quadratic constraint for its semi-definiteness, if not done yet.
 */
EXTERN
SCIP_RETCODE SCIPcheckCurvatureQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) convex.
 */
EXTERN
SCIP_Bool SCIPisConvexQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) concave.
 */
EXTERN
SCIP_Bool SCIPisConcaveQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the violation of a constraint by a solution. */
EXTERN
SCIP_RETCODE SCIPgetViolationQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< pointer to store violation of constraint */
   );

/** Indicates whether the quadratic constraint is local w.r.t. the current local bounds.
 *
 * That is, checks whether each variable with a square term is fixed and for each bilinear term at least one variable is fixed.
 */
EXTERN
SCIP_Bool SCIPisLinearLocalQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Adds the constraint to an NLPI problem. */
EXTERN
SCIP_RETCODE SCIPaddToNlpiProblemQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< NLPI problem where to add constraint */
   SCIP_HASHMAP*         scipvar2nlpivar,    /**< mapping from SCIP variables to variable indices in NLPI */
   SCIP_Bool             names               /**< whether to pass constraint names to NLPI */
   );

/** sets the left hand side of a quadratic constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint.
 */
EXTERN
SCIP_RETCODE SCIPchgLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   );

/** sets the right hand side of a quadratic constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint.
 */
EXTERN
SCIP_RETCODE SCIPchgRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   );

EXTERN
/** gets the feasibility of the quadratic constraint in the given solution */
SCIP_RETCODE SCIPgetFeasibilityQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_Real*            feasibility         /**< pointer to store the feasibility */
   );

/** gets the activity of the quadratic constraint in the given solution */
EXTERN
SCIP_RETCODE SCIPgetActivityQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_Real*            activity            /**< pointer to store the activity */
   );

/** changes the linear coefficient value for a given quadratic variable in a quadratic constraint data; if not
 *  available, it adds it
 *
 *  @note this is only allowed for original constraints and variables in problem creation stage
 */
EXTERN
SCIP_RETCODE SCIPchgLinearCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< quadratic variable */
   SCIP_Real             coef                /**< new coefficient */
   );

/** changes the square coefficient value for a given quadratic variable in a quadratic constraint data; if not
 *  available, it adds it
 *
 *  @note this is only allowed for original constraints and variables in problem creation stage
 */
EXTERN
SCIP_RETCODE SCIPchgSquareCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< quadratic variable */
   SCIP_Real             coef                /**< new coefficient */
   );

/** changes the bilinear coefficient value for a given quadratic variable in a quadratic constraint data; if not
 *  available, it adds it
 *
 *  @note this is only allowed for original constraints and variables in problem creation stage
 */
EXTERN
SCIP_RETCODE SCIPchgBilinCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var1,               /**< first quadratic variable */
   SCIP_VAR*             var2,               /**< second quadratic variable */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   );

/** returns the total number of bilinear terms that are contained in all quadratic constraints */
int SCIPgetNAllBilinearTermsQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns all bilinear terms that are contained in all quadratic constraints */
EXTERN
SCIP_RETCODE SCIPgetAllBilinearTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR** RESTRICT   x,                  /**< array to store first variable of each bilinear term */
   SCIP_VAR** RESTRICT   y,                  /**< array to second variable of each bilinear term */
   int* RESTRICT         nbilinterms,        /**< buffer to store the total number of bilinear terms */
   int* RESTRICT         nunderests,         /**< array to store the total number of constraints that require to underestimate a bilinear term */
   int* RESTRICT         noverests,          /**< array to store the total number of constraints that require to overestimate a bilinear term */
   SCIP_Real*            maxnonconvexity     /**< estimate of nonconvex eigenvalues of all quadratic constraints containing a bilinear term */
   );

/** adds a globally valid inequality of the form xcoef x <= ycoef y + constant for a bilinear term (x,y)
 *
 *  @note the indices of bilinear terms match with the entries of bilinear terms returned by SCIPgetAllBilinearTermsQuadratic
 */
EXTERN
SCIP_RETCODE SCIPaddBilinearIneqQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   int                   idx,                /**< index of the bilinear term */
   SCIP_Real             xcoef,              /**< x coefficient */
   SCIP_Real             ycoef,              /**< y coefficient */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool*            success             /**< buffer to store whether inequality has been accepted */
   );

/* @} */


#ifdef SCIP_PRIVATE_ROWPREP

/** creates a SCIP_ROWPREP datastructure
 *
 * Initial row represents 0 <= 0.
 */
EXTERN
SCIP_RETCODE SCIPcreateRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep,            /**< buffer to store pointer to rowprep */
   SCIP_SIDETYPE         sidetype,           /**< whether cut will be or lower-equal or larger-equal type */
   SCIP_Bool             local               /**< whether cut will be valid only locally */
);

/** frees a SCIP_ROWPREP datastructure */
EXTERN
void SCIPfreeRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep             /**< pointer that stores pointer to rowprep */
);

/** creates a copy of a SCIP_ROWPREP datastructure */
EXTERN
SCIP_RETCODE SCIPcopyRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        target,             /**< buffer to store pointer of rowprep copy */
   SCIP_ROWPREP*         source              /**< rowprep to copy */
);

/** ensures that rowprep has space for at least given number of additional terms
 *
 * Useful when knowing in advance how many terms will be added.
 */
EXTERN
SCIP_RETCODE SCIPensureRowprepSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   size                /**< number of additional terms for which to alloc space in rowprep */
);

/** prints a rowprep */
EXTERN
void SCIPprintRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be printed */
   FILE*                 file                /**< file to print to, or NULL for stdout */
);

/** adds a term coef*var to a rowprep */
EXTERN
SCIP_RETCODE SCIPaddRowprepTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             coef                /**< coefficient to add */
);

/** adds several terms coef*var to a rowprep */
EXTERN
SCIP_RETCODE SCIPaddRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   nvars,              /**< number of terms to add */
   SCIP_VAR**            vars,               /**< variables to add */
   SCIP_Real*            coefs               /**< coefficients to add */
);

/** adds constant value to side of rowprep */
EXTERN
void SCIPaddRowprepSide(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             side                /**< constant value to be added to side */
);

/** adds constant term to rowprep
 *
 * Substracts constant from side.
 */
EXTERN
void SCIPaddRowprepConstant(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             constant            /**< constant value to be added */
);

#ifdef NDEBUG
#define SCIPaddRowprepSide(rowprep, sideadd)  ((rowprep)->side += (sideadd))
#define SCIPaddRowprepConstant(rowprep, constant)  SCIPaddRowprepSide(rowprep, -(constant))
#endif

/** computes violation of cut in a given solution */
EXTERN
SCIP_Real SCIPgetRowprepViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_SOL*             sol                 /**< solution or NULL for LP solution */
);

/** Merge terms that use same variable and eliminate zero coefficients.
 *
 * Terms are sorted by variable (@see SCIPvarComp) after return.
 */
EXTERN
void SCIPmergeRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep             /**< rowprep to be cleaned up */
);

/* Cleans up and attempts to improve rowprep
 *
 * Drops small or large coefficients if coefrange is too large, if this can be done by relaxing the cut.
 * Scales coefficients up to reach minimal violation, if possible.
 * Scaling is omitted if violation is very small (ROWPREP_SCALEUP_VIOLNONZERO) or
 * maximal coefficient would become huge (ROWPREP_SCALEUP_MAXMAXCOEF).
 * Scales coefficients and side down if they are large and if the minimal violation is still reached.
 * Rounds coefficients close to integral values to integrals, if this can be done by relaxing the cut.
 * Rounds side within epsilon of 0 to 0.0 or +/-1.1*epsilon, whichever relaxes the cut least.
 *
 * After return, the terms in the rowprep will be sorted by absolute value of coefficient, in decreasing order.
 */
EXTERN
SCIP_RETCODE SCIPcleanupRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             maxcoefrange,       /**< maximal allowed coefficients range */
   SCIP_Real             minviol,            /**< minimal absolute violation the row should achieve (w.r.t. sol) */
   SCIP_Real*            coefrange,          /**< buffer to store coefrange of cleaned up cut, or NULL if not of interest */
   SCIP_Real*            viol                /**< buffer to store absolute violation of cleaned up cut in sol, or NULL if not of interest */
);

/** scales a rowprep
 *
 * @return Exponent of actually applied scaling factor, if written as 2^x.
 */
EXTERN
int SCIPscaleRowprep(
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be scaled */
   SCIP_Real             factor              /**< suggested scale factor */
);

/** generates a SCIP_ROW from a rowprep */
EXTERN
SCIP_RETCODE SCIPgetRowprepRowCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
);

/** generates a SCIP_ROW from a rowprep */
EXTERN
SCIP_RETCODE SCIPgetRowprepRowSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_SEPA*            sepa                /**< separator */
);

#endif

/* @} */

#ifdef __cplusplus
}
#endif

#endif
