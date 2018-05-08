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

/**@file    nlpi_all.c
 * @ingroup NLPIS
 * @brief   NLP interface that uses all available NLP interfaces
 * @author  Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpi_all.h"
#include "nlpi/nlpi.h"
#include "scip/pub_misc.h"

#include <string.h>

#define NLPI_NAME              "all"                       /* short concise name of solver */
#define NLPI_DESC              "NLP interface that uses all available NLP interfaces" /* description of solver */
#define NLPI_PRIORITY          -3000                       /* priority of NLP solver */

/*
 * Data structures
 */

struct SCIP_NlpiData
{
   SCIP_NLPI**           nlpis;              /**< array containing all nlpis */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   int                   nnlpis;             /**< total number of nlpis */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< message handler */
};

struct SCIP_NlpiProblem
{
   SCIP_NLPIPROBLEM**    nlpiproblems;       /**< array containing all nlpi problems */
   int                   nnlpiproblems;      /**< total number of nlpi problems */
   int                   bestidx;            /**< index of NLP solver with the best solution */
};

#ifdef SCIP_STATISTIC
static int _nnlps = 0;                       /**< number of NLPs that have been solved */
#endif

/*
 * Local methods
 */

/*
 * Callback methods of NLP solver interface
 */

/** copy method of NLP interface (called when SCIP copies plugins)
 *
 * input:
 *  - blkmem block memory in target SCIP
 *  - sourcenlpi the NLP interface to copy
 *  - targetnlpi buffer to store pointer to copy of NLP interface
 */
static
SCIP_DECL_NLPICOPY( nlpiCopyAll )
{
   SCIP_NLPIDATA* sourcedata;
   SCIP_Real infinity;

   assert(sourcenlpi != NULL);
   assert(targetnlpi != NULL);

   sourcedata = SCIPnlpiGetData(sourcenlpi);
   assert(sourcedata != NULL);
   assert(sourcedata->nnlpis > 1);
   assert(sourcedata->nlpis[0] != NULL);

   /* create target nlpis */
   SCIP_CALL( SCIPcreateNlpSolverAll(blkmem, targetnlpi, sourcedata->nlpis, sourcedata->nnlpis) );
   assert(*targetnlpi != NULL);

   SCIP_CALL( SCIPnlpiGetRealPar(sourcedata->nlpis[0], NULL, SCIP_NLPPAR_INFINITY, &infinity) );
   SCIP_CALL( SCIPnlpiSetRealPar(*targetnlpi, NULL, SCIP_NLPPAR_INFINITY, infinity) );
   SCIP_CALL( SCIPnlpiSetMessageHdlr(*targetnlpi, sourcedata->messagehdlr) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data
 *
 * input:
 *  - nlpi datastructure for solver interface
 */
static
SCIP_DECL_NLPIFREE( nlpiFreeAll )
{
   SCIP_NLPIDATA* data;
   int i;

   assert(nlpi != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   for( i = data->nnlpis - 1; i >= 0; --i )
   {
      SCIP_CALL( SCIPnlpiFree(&data->nlpis[i]) );
   }

   BMSfreeBlockMemoryArrayNull(data->blkmem, &data->nlpis, data->nnlpis);
   BMSfreeBlockMemory(data->blkmem, &data);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets pointer for NLP solver
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *
 * return: void pointer to solver
 */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerAll)
{
   assert(nlpi != NULL);

   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** creates a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer to store the problem data
 *  - name name of problem, can be NULL
 */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemAll)
{
   SCIP_NLPIDATA* data;
   int i;

   assert(nlpi    != NULL);
   assert(problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, problem) );
   if( *problem == NULL )
      return SCIP_NOMEMORY;

   /* initialize problem */
   BMSclearMemory((*problem));
   SCIP_ALLOC( BMSallocBlockMemoryArray(data->blkmem, &(*problem)->nlpiproblems, data->nnlpis) );
   (*problem)->nnlpiproblems = data->nnlpis;

   for( i = 0; i < data->nnlpis; ++i )
   {
      assert(data->nlpis[i] != NULL);
      SCIP_CALL( SCIPnlpiCreateProblem(data->nlpis[i], &((*problem)->nlpiproblems[i]), name) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** free a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer where problem data is stored
 */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemAll)
{
   SCIP_NLPIDATA* data;
   int i;

   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   for( i = 0; i < data->nnlpis; ++i )
   {
      assert(data->nlpis[i] != NULL);
      SCIP_CALL( SCIPnlpiFreeProblem(data->nlpis[i], &(*problem)->nlpiproblems[i]) );
   }

   BMSfreeBlockMemoryArrayNull(data->blkmem, &(*problem)->nlpiproblems, data->nnlpis);
   BMSfreeBlockMemory(data->blkmem, problem);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets pointer to solver-internal problem instance
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: void pointer to problem instance
 */
static
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerAll)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);

   return NULL;
}  /*lint !e715*/

/** add variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables
 *  - lbs lower bounds of variables, can be NULL if -infinity
 *  - ubs upper bounds of variables, can be NULL if +infinity
 *  - varnames names of variables, can be NULL
 */
static
SCIP_DECL_NLPIADDVARS( nlpiAddVarsAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiAddVars(nlpidata->nlpis[i], problem->nlpiproblems[i], nvars, lbs, ubs, varnames) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/


/** add constraints
 * quadratic coefficiens: row oriented matrix for each constraint
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - ncons number of added constraints
 *  - lhss left hand sides of constraints
 *  - rhss right hand sides of constraints
 *  - nlininds number of linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - lininds indices of variables for linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - linvals values of linear coefficient for each constraint
 *    may be NULL in case of no linear part
 *  - nquadrows number of columns in matrix of quadratic part for each constraint
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadrowidxs indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds[.] and quadvals[.]
 *    indices are given w.r.t. quadrowidxs., i.e., quadoffsets[.][i] gives the start index of row quadrowidxs[.][i] in quadvals[.]
 *    quadoffsets[.][nquadrows[.]] gives length of quadinds[.] and quadvals[.]
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadinds column indices w.r.t. quadrowidxs, i.e., quadrowidxs[quadinds[.][i]] gives the index of the variable corresponding
 *    to entry i, entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadvals coefficient values
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    entry of array may be NULL in case of no expression tree
 *    may be NULL in case of no expression tree in any constraint
 *  - exprtrees expression tree for nonquadratic part of constraints
 *    entry of array may be NULL in case of no nonquadratic part
 *    may be NULL in case of no nonquadratic part in any constraint
 *  - names of constraints, may be NULL or entries may be NULL
 */
static
SCIP_DECL_NLPIADDCONSTRAINTS( nlpiAddConstraintsAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiAddConstraints(nlpidata->nlpis[i], problem->nlpiproblems[i], ncons, lhss, rhss, nlininds,
            lininds, linvals, nquadelems, quadelems, exprvaridxs, exprtrees, names) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected
 *  May change sparsity pattern.
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nlins number of linear variables
 *  - lininds variable indices
 *    may be NULL in case of no linear part
 *  - linvals coefficient values
 *    may be NULL in case of no linear part
 *  - nquadcols number of columns in matrix of quadratic part
 *  - quadcols indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds and quadvals
 *    quadoffsets[.][nquadcols] gives length of quadinds and quadvals
 *    may be NULL in case of no quadratic part
 *  - quadinds column indices
 *    may be NULL in case of no quadratic part
 *  - quadvals coefficient values
 *    may be NULL in case of no quadratic part
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    may be NULL in case of no expression tree
 *  - exprtree expression tree for nonquadratic part of objective function
 *    may be NULL in case of no nonquadratic part
 *  - constant objective value offset
 */
static
SCIP_DECL_NLPISETOBJECTIVE( nlpiSetObjectiveAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiSetObjective(nlpidata->nlpis[i], problem->nlpiproblems[i], nlins, lininds, linvals, nquadelems,
            quadelems, exprvaridxs, exprtree, constant) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change variable bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables to change bounds
 *  - indices indices of variables to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
static
SCIP_DECL_NLPICHGVARBOUNDS( nlpiChgVarBoundsAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiChgVarBounds(nlpidata->nlpis[i], problem->nlpiproblems[i], nvars, indices, lbs, ubs) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change constraint bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nconss number of constraints to change sides
 *  - indices indices of constraints to change sides
 *  - lhss new left hand sides
 *  - rhss new right hand sides
 */
static
SCIP_DECL_NLPICHGCONSSIDES( nlpiChgConsSidesAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiChgConsSides(nlpidata->nlpis[i], problem->nlpiproblems[i], nconss, indices, lhss, rhss) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** delete a set of variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 *  - size of the dstats array
 *
 * output:
 *  - dstats new position of var, -1 if var was deleted
 */
static
SCIP_DECL_NLPIDELVARSET( nlpiDelVarSetAll )
{
   SCIP_NLPIDATA* nlpidata;
   int* tmpdstats;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_ALLOC( BMSallocBlockMemoryArray(nlpidata->blkmem, &tmpdstats, dstatssize) );

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      if( i < nlpidata->nnlpis -1 )
      {
         /* restore dstats entries */
         BMScopyMemoryArray(tmpdstats, dstats, dstatssize);

         SCIP_CALL( SCIPnlpiDelVarSet(nlpidata->nlpis[i], problem->nlpiproblems[i], tmpdstats, dstatssize) );
      }
      else
      {
         /* NOTE this works only when all dstats array are the same after calling the nlpidelvarset callback
          * As long as all solvers use the SCIP NLPI oracle to store the NLP problem data, this is the case.
          * @TODO Assert that the returned dstats are all the same?
          */
         SCIP_CALL( SCIPnlpiDelVarSet(nlpidata->nlpis[i], problem->nlpiproblems[i], dstats, dstatssize) );
      }
   }

   BMSfreeBlockMemoryArray(nlpidata->blkmem, &tmpdstats, dstatssize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** delete a set of constraints
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of rows; 1 if row should be deleted, 0 if not
 *  - size of the dstats array
 *
 * output:
 *  - dstats new position of row, -1 if row was deleted
 */
static
SCIP_DECL_NLPIDELCONSSET( nlpiDelConstraintSetAll )
{
   SCIP_NLPIDATA* nlpidata;
   int* tmpdstats;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   SCIP_ALLOC( BMSallocBlockMemoryArray(nlpidata->blkmem, &tmpdstats, dstatssize) );

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      if( i < nlpidata->nnlpis - 1 )
      {
         /* restore dstats entries */
         BMScopyMemoryArray(tmpdstats, dstats, dstatssize);

         SCIP_CALL( SCIPnlpiDelConsSet(nlpidata->nlpis[i], problem->nlpiproblems[i], tmpdstats, dstatssize) );
      }
      else
      {
         /* NOTE this works only when all dstats array are the same after calling the nlpidelconsset callback
          * As long as all solvers use the SCIP NLPI oracle to store the NLP problem data, this is the case.
          * @TODO Assert that the returned dstats are all the same?
          */
         SCIP_CALL( SCIPnlpiDelConsSet(nlpidata->nlpis[i], problem->nlpiproblems[i], dstats, dstatssize) );
      }

   }

   BMSfreeBlockMemoryArray(nlpidata->blkmem, &tmpdstats, dstatssize);

   return SCIP_OKAY;
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint to change
 *  - varidxs indices of variables which coefficient to change
 *  - vals new values for coefficients
 */
static
SCIP_DECL_NLPICHGLINEARCOEFS( nlpiChgLinearCoefsAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiChgLinearCoefs(nlpidata->nlpis[i], problem->nlpiproblems[i], idx, nvals, varidxs, vals) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** changes (or adds) coefficients in the quadratic part of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nentries number of entries in quadratic matrix to change
 *  - rows row indices of entries in quadratic matrix where values should be changed
 *  - cols column indices of entries in quadratic matrix where values should be changed
 *  - values new values for entries in quadratic matrix
 */
static
SCIP_DECL_NLPICHGQUADCOEFS( nlpiChgQuadraticCoefsAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiChgQuadCoefs(nlpidata->nlpis[i], problem->nlpiproblems[i], idx, nquadelems, quadelems) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** replaces the expression tree of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp, or NULL
 *  - exprtree new expression tree for constraint or objective, or NULL to only remove previous tree
 */
static
SCIP_DECL_NLPICHGEXPRTREE( nlpiChgExprtreeAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiChgExprtree(nlpidata->nlpis[i], problem->nlpiproblems[i], idxcons, exprvaridxs, exprtree) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change one coefficient in the nonlinear part
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - idxparam index of parameter
 *  - value new value for nonlinear parameter
 *
 * return: Error if parameter does not exist
 */
static
SCIP_DECL_NLPICHGNONLINCOEF( nlpiChgNonlinCoefAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiChgNonlinCoef(nlpidata->nlpis[i], problem->nlpiproblems[i], idxcons, idxparam, value) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** change the constant offset in the objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - objconstant new value for objective constant
 */
static
SCIP_DECL_NLPICHGOBJCONSTANT( nlpiChgObjConstantAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiChgObjConstant(nlpidata->nlpis[i], problem->nlpiproblems[i], objconstant) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets initial guess for primal variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues initial primal values for variables, or NULL to clear previous values
 *  - consdualvalues initial dual values for constraints, or NULL to clear previous values
 *  - varlbdualvalues  initial dual values for variable lower bounds, or NULL to clear previous values
 *  - varubdualvalues  initial dual values for variable upper bounds, or NULL to clear previous values
 */
static
SCIP_DECL_NLPISETINITIALGUESS( nlpiSetInitialGuessAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiSetInitialGuess(nlpidata->nlpis[i], problem->nlpiproblems[i], primalvalues, consdualvalues,
            varlbdualvalues, varubdualvalues) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** tries to solve NLP
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 */
static
SCIP_DECL_NLPISOLVE( nlpiSolveAll )
{
   SCIP_NLPIDATA* nlpidata;
   SCIP_NLPTERMSTAT besttermstat;
   SCIP_NLPSOLSTAT bestsolstat;
   SCIP_Real bestsolval;
   SCIP_Real infinity;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   /* use first solver per default */
   problem->bestidx = 0;

   /* initialize best solution values */
   SCIP_CALL( SCIPnlpiGetRealPar(nlpidata->nlpis[0], problem->nlpiproblems[0], SCIP_NLPPAR_INFINITY, &infinity) );
   besttermstat = SCIP_NLPTERMSTAT_OTHER;
   bestsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   bestsolval = infinity;

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      SCIP_NLPTERMSTAT termstat;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Real solval;
      SCIP_Bool update;

      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      /* solve NLP */
      SCIP_CALL( SCIPnlpiSolve(nlpidata->nlpis[i], problem->nlpiproblems[i]) );

      termstat = SCIPnlpiGetTermstat(nlpidata->nlpis[i], problem->nlpiproblems[i]);
      solstat = SCIPnlpiGetSolstat(nlpidata->nlpis[i], problem->nlpiproblems[i]);
      solval = infinity;
      update = FALSE;

      /* collect solution value */
      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         SCIP_CALL( SCIPnlpiGetSolution(nlpidata->nlpis[i], problem->nlpiproblems[i],
               NULL, NULL, NULL, NULL, &solval) );
         assert(solval != infinity); /*lint !e777*/
      }

      /* better termination status -> update best solver */
      if( termstat < besttermstat )
         update = TRUE;

      /* no feasible solutions have been found so far -> update best solver */
      else if( bestsolstat >= SCIP_NLPSOLSTAT_LOCINFEASIBLE && solstat <= SCIP_NLPSOLSTAT_LOCINFEASIBLE )
         update = TRUE;

      /* use solver with the better solution value */
      else if( solval < bestsolval )
         update = TRUE;

      /* update best solver */
      if( update )
      {
         besttermstat = termstat;
         bestsolstat = solstat;
         bestsolval = solval;
         problem->bestidx = i;
      }

#ifdef SCIP_STATISTIC
      {
         SCIP_NLPSTATISTICS stats;

         SCIP_CALL( SCIPnlpiGetStatistics(nlpidata->nlpis[i], problem->nlpiproblems[i], &stats) );

         SCIPstatisticMessage("%d solver %s termstat %d solstat %d solval %e iters %d time %g\n",
            _nnlps, SCIPnlpiGetName(nlpidata->nlpis[i]), termstat, solstat, solval,
            SCIPnlpStatisticsGetNIterations(&stats), SCIPnlpStatisticsGetTotalTime(&stats));
      }
#endif
   }

#ifdef SCIP_STATISTIC
   ++_nnlps;
#endif

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives solution status
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Solution Status
 */
static
SCIP_DECL_NLPIGETSOLSTAT( nlpiGetSolstatAll )
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* return the solution status of the first nlpi */
   return SCIPnlpiGetSolstat(nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx]);
}

/** gives termination reason
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Termination Status
 */
static
SCIP_DECL_NLPIGETTERMSTAT( nlpiGetTermstatAll )
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* return the solution status of the first nlpi */
   return SCIPnlpiGetTermstat(nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx]);
}

/** gives primal and dual solution values
 *
 * solver can return NULL in dual values if not available
 * but if solver provides dual values for one side of variable bounds, then it must also provide those for the other side
 *
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues buffer to store pointer to array to primal values, or NULL if not needed
 *  - consdualvalues buffer to store pointer to array to dual values of constraints, or NULL if not needed
 *  - varlbdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - varubdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - objval pointer store the objective value, or NULL if not needed
 */
static
SCIP_DECL_NLPIGETSOLUTION( nlpiGetSolutionAll )
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* return the solution status of the first nlpi */
   SCIP_CALL( SCIPnlpiGetSolution(nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx],
         primalvalues, consdualvalues, varlbdualvalues, varubdualvalues, objval) );

   return SCIP_OKAY;
}

/** gives solve statistics
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - statistics pointer to store statistics
 *
 * output:
 *  - statistics solve statistics
 */
static
SCIP_DECL_NLPIGETSTATISTICS( nlpiGetStatisticsAll )
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[problem->bestidx] != NULL);
   assert(problem->nlpiproblems != NULL);
   assert(problem->nlpiproblems[problem->bestidx] != NULL);

   /* collect statistics of the first solver */
   SCIP_CALL( SCIPnlpiGetStatistics(nlpidata->nlpis[problem->bestidx], problem->nlpiproblems[problem->bestidx],
         statistics) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives required size of a buffer to store a warmstart object
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - size pointer to store required size for warmstart buffer
 *
 * output:
 *  - size required size for warmstart buffer
 */
static
SCIP_DECL_NLPIGETWARMSTARTSIZE( nlpiGetWarmstartSizeAll )
{
   SCIPerrorMessage("method of NLP solver is not implemented\n");
   return SCIP_OKAY;
}  /*lint !e715*/

/** stores warmstart information in buffer
 *
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer memory to store warmstart information
 *
 * output:
 *  - buffer warmstart information in solver specific data structure
 */
static
SCIP_DECL_NLPIGETWARMSTARTMEMO( nlpiGetWarmstartMemoAll )
{
   SCIPerrorMessage("method of NLP solver is not implemented\n");
   return SCIP_OKAY;
}  /*lint !e715*/

/** sets warmstart information in solver
 *
 * write warmstart to buffer
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer warmstart information
 */
static
SCIP_DECL_NLPISETWARMSTARTMEMO( nlpiSetWarmstartMemoAll )
{
   SCIPerrorMessage("method of NLP solver is not implemented\n");
   return SCIP_OKAY;
}  /*lint !e715*/

/** gets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival pointer to store the parameter value
 *
 * output:
 *  - ival parameter value
 */
static
SCIP_DECL_NLPIGETINTPAR( nlpiGetIntParAll )
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[0] != NULL);

   /* take the first nlpi */
   SCIP_CALL( SCIPnlpiGetIntPar(nlpidata->nlpis[0], problem->nlpiproblems[0], type, ival) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival parameter value
 */
static
SCIP_DECL_NLPISETINTPAR( nlpiSetIntParAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiSetIntPar(nlpidata->nlpis[i], problem->nlpiproblems[i], type, ival) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval pointer to store the parameter value
 *
 * output:
 *  - dval parameter value
 */
static
SCIP_DECL_NLPIGETREALPAR( nlpiGetRealParAll )
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[0] != NULL);

   /* take the first nlpi */
   SCIP_CALL( SCIPnlpiGetRealPar(nlpidata->nlpis[0], problem->nlpiproblems[0], type, dval) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval parameter value
 */
static
SCIP_DECL_NLPISETREALPAR( nlpiSetRealParAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);

      if( type == SCIP_NLPPAR_INFINITY )
      {
         SCIP_CALL( SCIPnlpiSetRealPar(nlpidata->nlpis[i], NULL, type, dval) );
      }
      else
      {
         assert(problem->nlpiproblems != NULL);
         assert(problem->nlpiproblems[i] != NULL);

         SCIP_CALL( SCIPnlpiSetRealPar(nlpidata->nlpis[i], problem->nlpiproblems[i], type, dval) );
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval pointer to store the string value, the user must not modify the string
 *
 * output:
 *  - sval parameter value
 */
static
SCIP_DECL_NLPIGETSTRINGPAR( nlpiGetStringParAll )
{
   SCIP_NLPIDATA* nlpidata;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);
   assert(nlpidata->nlpis != NULL);
   assert(nlpidata->nlpis[0] != NULL);

   /* take the first nlpi */
   SCIP_CALL( SCIPnlpiGetStringPar(nlpidata->nlpis[0], problem->nlpiproblems[0], type, sval) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval parameter value
 */
static
SCIP_DECL_NLPISETSTRINGPAR( nlpiSetStringParAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);
      assert(problem->nlpiproblems[i] != NULL);

      SCIP_CALL( SCIPnlpiSetStringPar(nlpidata->nlpis[i], problem->nlpiproblems[i], type, sval) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets message handler for message output
 *
 * input:
 *  - nlpi NLP interface structure
 *  - messagehdlr SCIP message handler, or NULL to suppress all output
 */
static
SCIP_DECL_NLPISETMESSAGEHDLR( nlpiSetMessageHdlrAll )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   assert(nlpi != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   nlpidata->messagehdlr = messagehdlr;

   for( i = 0; i < nlpidata->nnlpis; ++i )
   {
      assert(nlpidata->nlpis[i] != NULL);

      SCIP_CALL( SCIPnlpiSetMessageHdlr(nlpidata->nlpis[i], messagehdlr) );
   }

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for All solver */
SCIP_RETCODE SCIPcreateNlpSolverAll(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi,               /**< pointer to buffer for nlpi address */
   SCIP_NLPI**           nlpis,              /**< array containing existing nlpis */
   int                   nnlpis              /**< total number of nlpis */
   )
{
   SCIP_NLPIDATA* nlpidata;
   int i;

   assert(blkmem != NULL);
   assert(nlpi != NULL);
   assert(nlpis != NULL || nnlpis == 0);

   /* the number of nlpis must be >= 2 */
   if( nnlpis < 2 )
   {
      *nlpi = NULL;
      return SCIP_OKAY;
   }
   assert(nlpis != NULL);

   /* create all solver interface data */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &nlpidata) );
   BMSclearMemory(nlpidata);
   nlpidata->blkmem = blkmem;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlpidata->nlpis, nnlpis) );

   /* copy nlpis */
   for( i = 0; i < nnlpis; ++i )
   {
      SCIP_CALL( SCIPnlpiCopy(blkmem, nlpis[i], &nlpidata->nlpis[i]) );
   }
   nlpidata->nnlpis = nnlpis;

   /* create solver interface */
   SCIP_CALL( SCIPnlpiCreate(nlpi,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyAll, nlpiFreeAll, nlpiGetSolverPointerAll,
         nlpiCreateProblemAll, nlpiFreeProblemAll, nlpiGetProblemPointerAll,
         nlpiAddVarsAll, nlpiAddConstraintsAll, nlpiSetObjectiveAll,
         nlpiChgVarBoundsAll, nlpiChgConsSidesAll, nlpiDelVarSetAll, nlpiDelConstraintSetAll,
         nlpiChgLinearCoefsAll, nlpiChgQuadraticCoefsAll, nlpiChgExprtreeAll, nlpiChgNonlinCoefAll,
         nlpiChgObjConstantAll, nlpiSetInitialGuessAll, nlpiSolveAll, nlpiGetSolstatAll, nlpiGetTermstatAll,
         nlpiGetSolutionAll, nlpiGetStatisticsAll,
         nlpiGetWarmstartSizeAll, nlpiGetWarmstartMemoAll, nlpiSetWarmstartMemoAll,
         nlpiGetIntParAll, nlpiSetIntParAll, nlpiGetRealParAll, nlpiSetRealParAll, nlpiGetStringParAll, nlpiSetStringParAll,
         nlpiSetMessageHdlrAll,
         nlpidata) );

   return SCIP_OKAY;
}
