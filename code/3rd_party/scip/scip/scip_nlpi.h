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

/**@file   scip_nlpi.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for NLPI solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_NLPI_H__
#define __SCIP_SCIP_NLPI_H__

#include "scip/type_nlpi.h"
#include "scip/type_misc.h"
#include "scip/type_lp.h"
#include "blockmemshell/memory.h"
#include "scip/pub_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNLPIInterfaceMethods
 *
 * @{
 */

/** creates an NLPI and includes it into SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlpi(
   SCIP*                           scip,                        /**< SCIP data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying an NLPI, can be NULL */
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

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_NLPI* SCIPfindNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of NLPI */
   );

/** returns the array of currently available NLPIs (sorted by priority) */
SCIP_EXPORT
SCIP_NLPI** SCIPgetNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available NLPIs */
SCIP_EXPORT
int SCIPgetNNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of an NLPI */
SCIP_EXPORT
SCIP_RETCODE SCIPsetNlpiPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of the NLPI */
   );

/** gets internal pointer to NLP solver
 *
 * Depending on the solver interface, a solver pointer may exist for every NLP problem instance.
 * For this case, a NLPI problem can be passed in as well.
 */
SCIP_EXPORT
SCIP_DECL_NLPIGETSOLVERPOINTER(SCIPgetNlpiSolverPointer);

/** creates an empty problem instance */
SCIP_EXPORT
SCIP_DECL_NLPICREATEPROBLEM(SCIPcreateNlpiProblem);

/** frees a problem instance */
SCIP_EXPORT
SCIP_DECL_NLPIFREEPROBLEM(SCIPfreeNlpiProblem);

/** gets internal pointer to solver-internal problem instance */
SCIP_EXPORT
SCIP_DECL_NLPIGETPROBLEMPOINTER(SCIPgetNlpiProblemPointer);

/** add variables to nlpi */
SCIP_EXPORT
SCIP_DECL_NLPIADDVARS(SCIPaddNlpiVars);

/** add constraints to nlpi */
SCIP_EXPORT
SCIP_DECL_NLPIADDCONSTRAINTS(SCIPaddNlpiConstraints);

/** sets or overwrites objective, a minimization problem is expected */
SCIP_EXPORT
SCIP_DECL_NLPISETOBJECTIVE(SCIPsetNlpiObjective);

/** change variable bounds */
SCIP_EXPORT
SCIP_DECL_NLPICHGVARBOUNDS(SCIPchgNlpiVarBounds);

/** change constraint sides */
SCIP_EXPORT
SCIP_DECL_NLPICHGCONSSIDES(SCIPchgNlpiConsSides);

/** delete a set of variables */
SCIP_EXPORT
SCIP_DECL_NLPIDELVARSET(SCIPdelNlpiVarSet);

/** delete a set of constraints */
SCIP_EXPORT
SCIP_DECL_NLPIDELCONSSET(SCIPdelNlpiConsSet);

/** changes or adds linear coefficients in a constraint or objective */
SCIP_EXPORT
SCIP_DECL_NLPICHGLINEARCOEFS(SCIPchgNlpiLinearCoefs);

/** change the expression in the nonlinear part */
SCIP_EXPORT
SCIP_DECL_NLPICHGEXPR(SCIPchgNlpiExpr);

/** change the constant offset in the objective */
SCIP_EXPORT
SCIP_DECL_NLPICHGOBJCONSTANT(SCIPchgNlpiObjConstant);

/** sets initial guess */
SCIP_EXPORT
SCIP_DECL_NLPISETINITIALGUESS(SCIPsetNlpiInitialGuess);

/** try to solve NLP with all parameters given as SCIP_NLPPARAM struct
 *
 * Typical use is
 *
 *     SCIP_NLPPARAM nlparam = { SCIP_NLPPARAM_DEFAULT(scip); }
 *     nlpparam.iterlimit = 42;
 *     SCIP_CALL( SCIPsolveNlpiParam(scip, nlpi, nlpiproblem, nlpparam) );
 *
 * or, in "one" line:
 *
 *     SCIP_CALL( SCIPsolveNlpiParam(scip, nlpi, nlpiproblem,
 *        (SCIP_NLPPARAM){ SCIP_NLPPARAM_DEFAULT(scip), .iterlimit = 42 }) );
 *
 * To get the latter, also \ref SCIPsolveNlpi can be used.
 */
SCIP_EXPORT
SCIP_DECL_NLPISOLVE(SCIPsolveNlpiParam);

/** try to solve NLP with non-default parameters given as optional arguments
 *
 * Typical use is
 *
 *     SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiproblem) );
 *
 * to solve with default parameters.
 * Additionally, one or several values of SCIP_NLPPARAM can be set:
 *
 *     SCIP_CALL( SCIPsolveNlpi(scip, nlpi, nlpiproblem, .iterlimit = 42, .verblevel = 1) );  //lint !e666
 */
/* the problem argument has been made part of the variadic arguments, since ISO C99 requires at least one argument for the "..." part and we want to allow leaving all parameters at default
 * for the same reason, we set the .caller argument, so that macro SCIP_VARARGS_REST will have at least one arg to return
 */
#if !defined(_MSC_VER) || _MSC_VER >= 1800
#define SCIPsolveNlpi(scip, nlpi, ...) \
   SCIPsolveNlpiParam(scip, nlpi, SCIP_VARARGS_FIRST((__VA_ARGS__, ignored)), \
      (SCIP_NLPPARAM){ SCIP_NLPPARAM_DEFAULT_INITS(scip), SCIP_VARARGS_REST(__VA_ARGS__, .caller = __FILE__) })
#else
/* very old MSVC doesn't support C99's designated initializers, so have a version of SCIPsolveNlpi() that just ignores given parameters
 * (compilation of scip_nlpi.c will print a warning)
 */
#define SCIPsolveNlpi(scip, nlpi, ...) \
    SCIPsolveNlpiParam(scip, nlpi, SCIP_VARARGS_FIRST((__VA_ARGS__, ignored)), SCIP_NLPPARAM_DEFAULT_STATIC)
#endif

/** gives solution status */
SCIP_EXPORT
SCIP_DECL_NLPIGETSOLSTAT(SCIPgetNlpiSolstat);

/** gives termination reason */
SCIP_EXPORT
SCIP_DECL_NLPIGETTERMSTAT(SCIPgetNlpiTermstat);

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_EXPORT
SCIP_DECL_NLPIGETSOLUTION(SCIPgetNlpiSolution);

/** gives solve statistics */
SCIP_EXPORT
SCIP_DECL_NLPIGETSTATISTICS(SCIPgetNlpiStatistics);


/**@name Convenience methods to setup and update an NLPI problem using NLROWS
 *
 * These methods can be used, for example, to create a NLPI problem that contains only the convex rows of the SCIP NLP relaxation.
 * @{
 */

/** creates a NLPI problem from given nonlinear rows
 *
 * The function computes for each variable the number of non-linear occurrences and stores it in the nlscore array.
 *
 * @note the first row corresponds always to the cutoff row (even if cutoffbound is SCIPinfinity(scip))
 **/
SCIP_EXPORT
SCIP_RETCODE SCIPcreateNlpiProblemFromNlRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM**    nlpiprob,           /**< buffer to store pointer to created nlpi problem */
   const char*           name,               /**< name to give to problem */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< number of nonlinear rows */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_HASHMAP*         nlrow2idx,          /**< empty hash map to store mapping between variables and indices in nlpiprob, can be NULL */
   SCIP_Real*            nlscore,            /**< array to store the score of each nonlinear variable (NULL if not needed) */
   SCIP_Real             cutoffbound,        /**< cutoff bound */
   SCIP_Bool             setobj,             /**< whether the objective function should be set to one of the SCIP problem */
   SCIP_Bool             onlyconvex          /**< filter only for convex constraints */
   );

/** updates variable bounds and the cutoff row in a NLPI problem
 *
 * The NLPI problem must have been setup by SCIPcreateNlpiProblemFromNlRows().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateNlpiProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_HASHMAP*         var2nlpiidx,        /**< mapping between variables and nlpi indices */
   SCIP_VAR**            nlpivars,           /**< array containing all variables of the nlpi */
   int                   nlpinvars,          /**< total number of nlpi variables */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   );

/** adds SCIP_ROWs to a NLPI problem */
SCIP_EXPORT
SCIP_RETCODE SCIPaddNlpiProblemRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_ROW**            rows,               /**< rows to add */
   int                   nrows               /**< number of rows to add */
   );

/** adds SCIP_NLROWs to a NLPI problem */
SCIP_EXPORT
SCIP_RETCODE SCIPaddNlpiProblemNlRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpiprob */
   SCIP_NLROW**          nlrows,             /**< rows to add */
   int                   nnlrows             /**< number of rows to add */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_SCIP_NLPI_H__ */
