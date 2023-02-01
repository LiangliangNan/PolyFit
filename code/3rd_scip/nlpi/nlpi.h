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

/**@file   nlpi.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for NLPI solver interfaces
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_H__
#define __SCIP_NLPI_H__

#include "nlpi/type_nlpi.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup NLPIS
 *
 * @{
 */

/** compares two NLPIs w.r.t. their priority */
SCIP_DECL_SORTPTRCOMP(SCIPnlpiComp);

/** creates an NLP solver interface */
EXTERN
SCIP_RETCODE SCIPnlpiCreate(
   SCIP_NLPI**                     nlpi,                        /**< pointer to NLP interface data structure */
   const char*                     name,                        /**< name of NLP interface */
   const char*                     description,                 /**< description of NLP interface */
   int                             priority,                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy)),               /**< copying an NLPI */
   SCIP_DECL_NLPIFREE              ((*nlpifree)),               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer)),   /**< get solver pointer */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem)),      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem)),        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer)),  /**< get problem pointer */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars)),            /**< add variables */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints)),     /**< add constraints */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective)),       /**< set objective */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds)),       /**< change variable bounds */
   SCIP_DECL_NLPICHGCONSSIDES      ((*nlpichgconssides)),       /**< change constraint sides */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset)),          /**< delete a set of constraints */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset)),         /**< delete a set of constraints */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs)),     /**< change coefficients in linear part of a constraint or objective */
   SCIP_DECL_NLPICHGQUADCOEFS      ((*nlpichgquadcoefs)),       /**< change coefficients in quadratic part of a constraint or objective */
   SCIP_DECL_NLPICHGEXPRTREE       ((*nlpichgexprtree)),        /**< change nonlinear expression a constraint or objective */
   SCIP_DECL_NLPICHGNONLINCOEF     ((*nlpichgnonlincoef)),      /**< change one parameter in nonlinear expressions of a constraint or objective */
   SCIP_DECL_NLPICHGOBJCONSTANT    ((*nlpichgobjconstant)),     /**< change the constant offset in the objective */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess)),    /**< set initial guess for primal variables */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve)),              /**< solve NLP */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat)),         /**< get solution status */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat)),        /**< get termination status */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution)),        /**< get solution */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics)),      /**< get solve statistics */
   SCIP_DECL_NLPIGETWARMSTARTSIZE  ((*nlpigetwarmstartsize)),   /**< get size for warmstart object buffer */
   SCIP_DECL_NLPIGETWARMSTARTMEMO  ((*nlpigetwarmstartmemo)),   /**< get warmstart object */
   SCIP_DECL_NLPISETWARMSTARTMEMO  ((*nlpisetwarmstartmemo)),   /**< set warmstart object */
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar)),          /**< get value of integer parameter */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar)),          /**< set value of integer parameter */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar)),         /**< get value of floating point parameter */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar)),         /**< set value of floating point parameter */
   SCIP_DECL_NLPIGETSTRINGPAR      ((*nlpigetstringpar)),       /**< get value of string parameter */
   SCIP_DECL_NLPISETSTRINGPAR      ((*nlpisetstringpar)),       /**< set value of string parameter */
   SCIP_DECL_NLPISETMESSAGEHDLR    ((*nlpisetmessagehdlr)),     /**< set message handler */
   SCIP_NLPIDATA*                  nlpidata                     /**< NLP interface local data */
   );

/** copies an NLPI */
EXTERN
SCIP_RETCODE SCIPnlpiCopy(
   BMS_BLKMEM*           blkmem,             /**< block memory in target SCIP */
   SCIP_NLPI*            sourcenlpi,         /**< pointer to NLPI data structure to copy */
   SCIP_NLPI**           targetnlpi          /**< buffer to store pointer to copied NLPI data structure */
   );

/** frees NLPI user data */
EXTERN
SCIP_RETCODE SCIPnlpiFree(
   SCIP_NLPI**           nlpi                /**< pointer to NLPI data structure */
   );

/** gets pointer for NLP solver
 * @return void pointer to solver
 */
EXTERN
void* SCIPnlpiGetSolverPointer(
   SCIP_NLPI*            nlpi                /**< pointer to NLPI datastructure */
   );

/** creates a problem instance */
EXTERN
SCIP_RETCODE SCIPnlpiCreateProblem(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM**    problem,            /**< pointer to store problem data */
   const char*           name                /**< name of problem, can be NULL */
   );

/** frees a problem instance */
EXTERN
SCIP_RETCODE SCIPnlpiFreeProblem(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM**    problem             /**< pointer where problem data is stored */
   );

/** gets pointer to solver-internal problem instance
 * @return void pointer to problem instance
 */
EXTERN
void* SCIPnlpiGetProblemPointer(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer where problem data is stored */
   );

/** add variables to nlpi */
EXTERN
SCIP_RETCODE SCIPnlpiAddVars(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nvars,              /**< number of variables */               
   const SCIP_Real*      lbs,                /**< lower bounds of variables, can be NULL if -infinity */
   const SCIP_Real*      ubs,                /**< upper bounds of variables, can be NULL if +infinity */
   const char**          varnames            /**< varnames names of variables, can be NULL */
   );

/** add constraints to nlpi */
EXTERN
SCIP_RETCODE SCIPnlpiAddConstraints(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI data structure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nconss,             /**< number of added constraints */
   const SCIP_Real*      lhss,               /**< left hand sides of constraints, can be NULL if -infinity */
   const SCIP_Real*      rhss,               /**< right hand sides of constraints, can be NULL if +infinity */
   const int*            nlininds,           /**< number of linear coefficients for each constraint, may be NULL in case of no linear part */
   int* const*           lininds,            /**< indices of variables for linear coefficients for each constraint, may be NULL in case of no linear part */
   SCIP_Real* const*     linvals,            /**< values of linear coefficient for each constraint, may be NULL in case of no linear part */
   const int*            nquadelems,         /**< number of elements in matrix of quadratic part for each constraint,
                                              * may be NULL in case of no quadratic part in any constraint */
   SCIP_QUADELEM* const* quadelems,          /**< quadratic elements specifying quadratic part for each constraint, entry of array may be NULL in case of no quadratic part,
                                              * may be NULL in case of no quadratic part in any constraint */
   int* const*           exprvaridxs,        /**< indices of variables in expression tree, maps variable indices in expression
                                              * tree to indices in nlp, entry of array may be NULL in case of no expression
                                              * tree, may be NULL in case of no expression tree in any constraint */
   SCIP_EXPRTREE* const* exprtrees,          /**< exprtrees expression tree for nonquadratic part of constraints, entry of
                                              * array may be NULL in case of no nonquadratic part, may be NULL in case of no
                                              * nonquadratic part in any constraint */
   const char**          names               /**< names of constraints, may be NULL or entries may be NULL */
   );


/** sets or overwrites objective, a minimization problem is expected */
EXTERN
SCIP_RETCODE SCIPnlpiSetObjective(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nlins,              /**< number of linear variables */
   const int*            lininds,            /**< variable indices, may be NULL in case of no linear part */
   const SCIP_Real*      linvals,            /**< coefficient values, may be NULL in case of no linear part */
   int                   nquadelems,         /**< number of entries in matrix of quadratic part */
   const SCIP_QUADELEM*  quadelems,          /**< entries in matrix of quadratic part, may be NULL in case of no quadratic part */
   const int*            exprvaridxs,        /**< indices of variables in expression tree, maps variable indices in expression
                                              * tree to indices in nlp, may be NULL in case of no expression tree */
   const SCIP_EXPRTREE*  exprtree,           /**< expression tree for nonquadratic part of objective function, may be NULL in
                                              * case of no nonquadratic part */
   const SCIP_Real       constant            /**< objective value offset*/
   );

/** change variable bounds */
EXTERN
SCIP_RETCODE SCIPnlpiChgVarBounds(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nvars,              /**< number of variables to change bounds */
   const int*            indices,            /**< indices of variables to change bounds */
   const SCIP_Real*      lbs,                /**< new lower bounds */
   const SCIP_Real*      ubs                 /**< new upper bounds */
   );

/** change constraint sides */
EXTERN
SCIP_RETCODE SCIPnlpiChgConsSides(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   nconss,             /**< number of constraints to change sides */
   const int*            indices,            /**< indices of constraints to change sides */
   const SCIP_Real*      lhss,               /**< new left hand sides */
   const SCIP_Real*      rhss                /**< new right hand sides */
   );

/** delete a set of variables */
EXTERN
SCIP_RETCODE SCIPnlpiDelVarSet(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int*                  dstats,             /**< deletion status of vars; 1 if var should be deleted, 0 if not; afterwards -1
                                              * if var was deleted */
   int                   dstatssize          /**< size of the dstats array */
   );

/** delete a set of constraints */
EXTERN
SCIP_RETCODE SCIPnlpiDelConsSet(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int*                  dstats,             /**< deletion status of rows; 1 if row should be deleted, 0 if not; afterwards -1
                                              * if row was deleted */
   int                   dstatssize          /**< size of the dstats array */
   );

/** changes or adds linear coefficients in a constraint or objective */
EXTERN
SCIP_RETCODE SCIPnlpiChgLinearCoefs(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   const int             idx,                /**< index of constraint or -1 for objective */
   int                   nvals,              /**< number of values in linear constraint */
   const int*            varidxs,            /**< indices of variable */
   const SCIP_Real*      vals                /**< new values for coefficient */
   );

/** changes or adds coefficients in the quadratic part of a constraint or objective */
EXTERN
SCIP_RETCODE SCIPnlpiChgQuadCoefs(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   idx,                /**< index of constraint or -1 for objective */
   int                   nquadelems,         /**< number of entries in quadratic constraint to change */
   const SCIP_QUADELEM*  quadelems           /**< new elements in quadratic matrix (replacing already existing ones or adding new ones) */
   );

/** change the expression tree in the nonlinear part */
EXTERN
SCIP_RETCODE SCIPnlpiChgExprtree(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   idxcons,            /**< index of constraint or -1 for objective */
   const int*            exprvaridxs,        /**< indices of variables in expression tree, maps variable indices in expression tree to indices in nlp, or NULL */
   const SCIP_EXPRTREE*  exprtree            /**< new expression tree, or NULL for no tree */
   );

/** change the value of one parameter in the nonlinear part */
EXTERN
SCIP_RETCODE SCIPnlpiChgNonlinCoef(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   int                   idxcons,            /**< index of constraint or -1 for objective */
   int                   idxparam,           /**< index of parameter */
   SCIP_Real             value               /**< new value for nonlinear parameter */
   );

/** change the constant offset in the objective */
EXTERN
SCIP_RETCODE SCIPnlpiChgObjConstant(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_Real             objconstant         /**< new value for objective constant */
   );

/** sets initial guess for primal variables */
EXTERN
SCIP_RETCODE SCIPnlpiSetInitialGuess(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_Real*            primalvalues,       /**< initial primal values for variables, or NULL to clear previous values */
   SCIP_Real*            consdualvalues,     /**< initial dual values for constraints, or NULL to clear previous values */
   SCIP_Real*            varlbdualvalues,    /**< initial dual values for variable lower bounds, or NULL to clear previous values */
   SCIP_Real*            varubdualvalues     /**< initial dual values for variable upper bounds, or NULL to clear previous values */
   );

/** tries to solve NLP */
EXTERN
SCIP_RETCODE SCIPnlpiSolve(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   );
   
/** gives solution status
 * @return solution status */
EXTERN
SCIP_NLPSOLSTAT SCIPnlpiGetSolstat(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   );

/** gives termination reason
 * @return termination status */
EXTERN
SCIP_NLPTERMSTAT SCIPnlpiGetTermstat(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   );

/** gives primal and dual solution
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 */
EXTERN
SCIP_RETCODE SCIPnlpiGetSolution(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_Real**           primalvalues,       /**< buffer to store pointer to array to primal values, or NULL if not needed */
   SCIP_Real**           consdualvalues,     /**< buffer to store pointer to array to dual values of constraints, or NULL if not needed */
   SCIP_Real**           varlbdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real**           varubdualvalues,    /**< buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed */
   SCIP_Real*            objval              /**< buffer to store the objective value, or NULL if not needed */
   );

/** gives solve statistics */
EXTERN
SCIP_RETCODE SCIPnlpiGetStatistics(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   );

/** gives required size of a buffer to store a warmstart object */
EXTERN
SCIP_RETCODE SCIPnlpiGetWarmstartSize(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   size_t*               size                /**< pointer to store required size for warmstart buffer */
   );

/** stores warmstart information in buffer */
EXTERN
SCIP_RETCODE SCIPnlpiGetWarmstartMemo(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   void*                 buffer              /**< memory to store warmstart information */
   );

/** sets warmstart information in solver */
EXTERN
SCIP_RETCODE SCIPnlpiSetWarmstartMemo(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   void*                 buffer              /**< warmstart information */
   );

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of NLP */
EXTERN
SCIP_RETCODE SCIPnlpiGetIntPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   );
   
/** sets integer parameter of NLP */
EXTERN
SCIP_RETCODE SCIPnlpiSetIntPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** gets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then gets solver-wide value for infinity */
EXTERN
SCIP_RETCODE SCIPnlpiGetRealPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure, can be NULL only if type == SCIP_NLPPAR_INFINITY */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   );

/** sets floating point parameter of NLP
 * if problem is NULL and type == SCIP_NLPPAR_INFINITY, then sets solver-wide value for infinity */
EXTERN
SCIP_RETCODE SCIPnlpiSetRealPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure, can be NULL only if type == SCIP_NLPPAR_INFINITY */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/** gets string parameter of NLP */
EXTERN
SCIP_RETCODE SCIPnlpiGetStringPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value, the user must not modify the string */
   );

/** sets string parameter of NLP */
EXTERN
SCIP_RETCODE SCIPnlpiSetStringPar(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
   );

/** sets message handler for message output */
EXTERN
SCIP_RETCODE SCIPnlpiSetMessageHdlr(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< pointer to message handler, or NULL to suppress all output */
   );

/** gets data of an NLPI */
SCIP_NLPIDATA* SCIPnlpiGetData(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );
   
/** gets NLP solver name */
const char* SCIPnlpiGetName(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver descriptions */
const char* SCIPnlpiGetDesc(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** gets NLP solver priority */
int SCIPnlpiGetPriority(
   SCIP_NLPI*            nlpi                /**< NLP interface structure */
   );

/** sets NLP solver priority */
void SCIPnlpiSetPriority(
   SCIP_NLPI*            nlpi,               /**< NLP interface structure */
   int                   priority            /**< new priority of NLPI */
   );

/** creates an NLP statistics structure */
SCIP_RETCODE SCIPnlpStatisticsCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
   );

/** frees an NLP statistics structure */
void SCIPnlpStatisticsFree(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_NLPSTATISTICS**  statistics          /**< pointer where to store NLP statistics structure */
   );

/** gets the number of iterations from an NLP statistics structure */
int SCIPnlpStatisticsGetNIterations(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
   );

/** gets the total time from an NLP statistics structure */
SCIP_Real SCIPnlpStatisticsGetTotalTime(
   SCIP_NLPSTATISTICS*   statistics          /**< NLP statistics structure */
   );

/** sets the number of iterations in an NLP statistics structure */
void SCIPnlpStatisticsSetNIterations(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   int                   niterations         /**< number of iterations to store */
   );

/** sets the total time in an NLP statistics structure */
void SCIPnlpStatisticsSetTotalTime(
   SCIP_NLPSTATISTICS*   statistics,         /**< NLP statistics structure */
   SCIP_Real             totaltime           /**< solution time to store */
   );

/** @} */

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_H__ */
