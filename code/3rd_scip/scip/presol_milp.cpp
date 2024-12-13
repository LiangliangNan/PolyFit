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

/**@file   presol_milp.cpp
 * @brief  MILP presolver
 * @author Leona Gottwald
 *
 * Calls the presolve library and communicates (multi-)aggregations, fixings, and bound
 * changes to SCIP by utilizing the postsolve information. Constraint changes can currently
 * only be communicated by deleting all constraints and adding new ones.
 *
 * @todo add infrastructure to SCIP for handling parallel columns
 * @todo better communication of constraint changes by adding more information to the postsolve structure
 * @todo allow to pass additional external locks to the presolve library that are considered when doing reductions
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/presol_milp.h"

#ifndef SCIP_WITH_PAPILO

/** creates the MILP presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolMILP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   return SCIP_OKAY;
}

#else

/* disable some warnings that come up in header files of PAPILOs dependencies */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#pragma GCC diagnostic ignored "-Wredundant-decls"
#endif

#include <assert.h>
#include "scip/cons_linear.h"
#include "scip/pub_matrix.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/scip_general.h"
#include "scip/scip_presol.h"
#include "scip/scip_var.h"
#include "scip/scip_mem.h"
#include "scip/scip_prob.h"
#include "scip/scip_param.h"
#include "scip/scip_cons.h"
#include "scip/scip_numerics.h"
#include "scip/scip_timing.h"
#include "scip/scip_message.h"
#include "scip/scip_randnumgen.h"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/Config.hpp"

#define PRESOL_NAME            "milp"
#define PRESOL_DESC            "MILP specific presolving methods"
#define PRESOL_PRIORITY        9999999 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_MEDIUM /* timing of the presolver (fast, medium, or exhaustive) */

/* default parameter values */
#define DEFAULT_THREADS            1         /**< maximum number of threads presolving may use (0: automatic) */
#define DEFAULT_MAXFILLINPERSUBST  3         /**< maximal possible fillin for substitutions to be considered */
#define DEFAULT_MAXSHIFTPERROW     10        /**< maximal amount of nonzeros allowed to be shifted to make space for substitutions */
#define DEFAULT_DETECTLINDEP       0         /**< should linear dependent equations and free columns be removed? (0: never, 1: for LPs, 2: always) */
#define DEFAULT_MAXBADGESIZE_SEQ   15000     /**< the max badge size in Probing if PaPILO is executed in sequential mode*/
#define DEFAULT_MAXBADGESIZE_PAR   -1        /**< the max badge size in Probing if PaPILO is executed in parallel mode*/
#define DEFAULT_RANDOMSEED         0         /**< the random seed used for randomization of tie breaking */
#define DEFAULT_MODIFYCONSFAC      0.8       /**< modify SCIP constraints when the number of nonzeros or rows is at most this
                                              *   factor times the number of nonzeros or rows before presolving */
#define DEFAULT_MARKOWITZTOLERANCE 0.01      /**< the markowitz tolerance used for substitutions */
#define DEFAULT_HUGEBOUND          1e8       /**< absolute bound value that is considered too huge for activitity based calculations */
#define DEFAULT_ENABLEPARALLELROWS TRUE      /**< should the parallel rows presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEDOMCOL       TRUE      /**< should the dominated column presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEDUALINFER    TRUE      /**< should the dualinfer presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEMULTIAGGR    TRUE      /**< should the multi-aggregation presolver be enabled within the presolve library? */
#define DEFAULT_ENABLEPROBING      TRUE      /**< should the probing presolver be enabled within the presolve library? */
#define DEFAULT_ENABLESPARSIFY     FALSE     /**< should the sparsify presolver be enabled within the presolve library? */
#define DEFAULT_FILENAME_PROBLEM   "-"       /**< default filename to store the instance before presolving */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int lastncols;                            /**< the number of columns from the last call */
   int lastnrows;                            /**< the number of rows from the last call */
   int threads;                              /**< maximum number of threads presolving may use (0: automatic) */
   int maxfillinpersubstitution;             /**< maximal possible fillin for substitutions to be considered */
   int maxbadgesizeseq;                      /**< the max badge size in Probing if PaPILO is called in sequential mode*/
   int maxbadgesizepar;                      /**< the max badge size in Probing if PaPILO is called in parallel mode */
   int maxshiftperrow;                       /**< maximal amount of nonzeros allowed to be shifted to make space for substitutions */
   int detectlineardependency;               /**< should linear dependent equations and free columns be removed? (0: never, 1: for LPs, 2: always) */
   int randomseed;                           /**< the random seed used for randomization of tie breaking */
   SCIP_Bool enablesparsify;                 /**< should the sparsify presolver be enabled within the presolve library? */
   SCIP_Bool enabledomcol;                   /**< should the dominated column presolver be enabled within the presolve library? */
   SCIP_Bool enableprobing;                  /**< should the probing presolver be enabled within the presolve library? */
   SCIP_Bool enabledualinfer;                /**< should the dualinfer presolver be enabled within the presolve library? */
   SCIP_Bool enablemultiaggr;                /**< should the multi-aggregation presolver be enabled within the presolve library? */
   SCIP_Bool enableparallelrows;             /**< should the parallel rows presolver be enabled within the presolve library? */
   SCIP_Real modifyconsfac;                  /**< modify SCIP constraints when the number of nonzeros or rows is at most this
                                              *   factor times the number of nonzeros or rows before presolving */
   SCIP_Real markowitztolerance;             /**< the markowitz tolerance used for substitutions */
   SCIP_Real hugebound;                      /**< absolute bound value that is considered too huge for activitity based calculations */

   char* filename = NULL;                    /**< filename to store the instance before presolving */

};

using namespace papilo;

/*
 * Local methods
 */

/** builds the presolvelib problem datastructure from the matrix */
static
Problem<SCIP_Real> buildProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix              /**< initialized SCIP_MATRIX data structure */
   )
{
   ProblemBuilder<SCIP_Real> builder;

   /* build problem from matrix */
   int nnz = SCIPmatrixGetNNonzs(matrix);
   int ncols = SCIPmatrixGetNColumns(matrix);
   int nrows = SCIPmatrixGetNRows(matrix);
   builder.reserve(nnz, nrows, ncols);

   /* set up columns */
   builder.setNumCols(ncols);
   for(int i = 0; i != ncols; ++i)
   {
      SCIP_VAR* var = SCIPmatrixGetVar(matrix, i);
      SCIP_Real lb = SCIPvarGetLbGlobal(var);
      SCIP_Real ub = SCIPvarGetUbGlobal(var);
      builder.setColLb(i, lb);
      builder.setColUb(i, ub);
      builder.setColLbInf(i, SCIPisInfinity(scip, -lb));
      builder.setColUbInf(i, SCIPisInfinity(scip, ub));
#if PAPILO_VERSION_MAJOR > 2 || (PAPILO_VERSION_MAJOR == 2 && PAPILO_VERSION_MINOR >= 1)
      if ( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
         builder.setColImplInt(i, TRUE);
      else
         builder.setColIntegral(i, SCIPvarIsIntegral(var));
#else
      builder.setColIntegral(i, SCIPvarIsIntegral(var));
#endif
      builder.setObj(i, SCIPvarGetObj(var));
   }

   /* set up rows */
   builder.setNumRows(nrows);
   for(int i = 0; i != nrows; ++i)
   {
      int* rowcols = SCIPmatrixGetRowIdxPtr(matrix, i);
      SCIP_Real* rowvals = SCIPmatrixGetRowValPtr(matrix, i);
      int rowlen = SCIPmatrixGetRowNNonzs(matrix, i);
      builder.addRowEntries(i, rowlen, rowcols, rowvals);

      SCIP_Real lhs = SCIPmatrixGetRowLhs(matrix, i);
      SCIP_Real rhs = SCIPmatrixGetRowRhs(matrix, i);
      builder.setRowLhs(i, lhs);
      builder.setRowRhs(i, rhs);
      builder.setRowLhsInf(i, SCIPisInfinity(scip, -lhs));
      builder.setRowRhsInf(i, SCIPisInfinity(scip, rhs));
   }

   /* init objective offset - the value itself is irrelevant */
   builder.setObjOffset(0);

   return builder.build();
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyMILP)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludePresolMILP(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeMILP)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* data = SCIPpresolGetData(presol);
   assert(data != NULL);

   SCIPpresolSetData(presol, NULL);
   SCIPfreeBlockMemory(scip, &data);
   return SCIP_OKAY;
}

/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitMILP)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* data = SCIPpresolGetData(presol);
   assert(data != NULL);

   data->lastncols = -1;
   data->lastnrows = -1;

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecMILP)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_PRESOLDATA* data;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_Bool infeasible;
   SCIP_Real timelimit;

   *result = SCIP_DIDNOTRUN;

   data = SCIPpresolGetData(presol);

   int nvars = SCIPgetNVars(scip);
   int nconss = SCIPgetNConss(scip);

   /* run only if the problem size reduced by some amount since the last call or if it is the first call */
   if( data->lastncols != -1 && data->lastnrows != -1 &&
       nvars > data->lastncols * 0.85 &&
       nconss > data->lastnrows * 0.85 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, TRUE, &initialized, &complete, &infeasible,
      naddconss, ndelconss, nchgcoefs, nchgbds, nfixedvars) );

   /* if infeasibility was detected during matrix creation, return here */
   if( infeasible )
   {
      if( initialized )
         SCIPmatrixFree(scip, &matrix);

      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* we only work on pure MIPs, also disable to try building the matrix again if it failed once */
   if( !initialized || !complete )
   {
      data->lastncols = 0;
      data->lastnrows = 0;

      if( initialized )
         SCIPmatrixFree(scip, &matrix);

      return SCIP_OKAY;
   }

   /* only allow communication of constraint modifications by deleting all constraints when they have not been upgraded yet */
   SCIP_CONSHDLR* linconshdlr = SCIPfindConshdlr(scip, "linear");
   assert(linconshdlr != NULL);
   bool allowconsmodification = (SCIPconshdlrGetNCheckConss(linconshdlr) == SCIPmatrixGetNRows(matrix));

   Problem<SCIP_Real> problem = buildProblem(scip, matrix);
   Presolve<SCIP_Real> presolve;

   /* store current numbers of aggregations, fixings, and changed bounds for statistics */
   int oldnaggrvars = *naggrvars;
   int oldnfixedvars = *nfixedvars;
   int oldnchgbds = *nchgbds;

   /* important so that SCIP does not throw an error, e.g. when an integer variable is substituted
    * into a knapsack constraint */
   presolve.getPresolveOptions().substitutebinarieswithints = false;

   /* currently these changes cannot be communicated to SCIP correctly since a constraint needs
    * to be modified in the cases where slackvariables are removed from constraints but for the
    * presolve library those look like normal substitution on the postsolve stack */
   presolve.getPresolveOptions().removeslackvars = false;

   /* communicate the SCIP parameters to the presolve libary */
   presolve.getPresolveOptions().maxfillinpersubstitution = data->maxfillinpersubstitution;
   presolve.getPresolveOptions().markowitz_tolerance = data->markowitztolerance;
   presolve.getPresolveOptions().maxshiftperrow = data->maxshiftperrow;
   presolve.getPresolveOptions().hugeval = data->hugebound;

   /* removal of linear dependent equations has only an effect when constraint modifications are communicated */
   presolve.getPresolveOptions().detectlindep = allowconsmodification ? data->detectlineardependency : 0;

   /* communicate the random seed */
   presolve.getPresolveOptions().randomseed = SCIPinitializeRandomSeed(scip, (unsigned int)data->randomseed);

   /* set number of threads to be used for presolve */
   presolve.getPresolveOptions().threads = data->threads;

   /* disable dual reductions that are not permitted */
   if( !complete )
      presolve.getPresolveOptions().dualreds = 0;
   else if( SCIPallowStrongDualReds(scip) )
      presolve.getPresolveOptions().dualreds = 2;
   else if( SCIPallowWeakDualReds(scip) )
      presolve.getPresolveOptions().dualreds = 1;
   else
      presolve.getPresolveOptions().dualreds = 0;

   /* set up the presolvers that shall participate */
   using uptr = std::unique_ptr<PresolveMethod<SCIP_Real>>;

   /* fast presolvers*/
   presolve.addPresolveMethod( uptr( new SingletonCols<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new CoefficientStrengthening<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new ConstraintPropagation<SCIP_Real>() ) );

   /* medium presolver */
   presolve.addPresolveMethod( uptr( new SimpleProbing<SCIP_Real>() ) );
   if( data->enableparallelrows )
      presolve.addPresolveMethod( uptr( new ParallelRowDetection<SCIP_Real>() ) );
   /* todo: parallel cols cannot be handled by SCIP currently
   * addPresolveMethod( uptr( new ParallelColDetection<SCIP_Real>() ) ); */
   presolve.addPresolveMethod( uptr( new SingletonStuffing<SCIP_Real>() ) );
#if PAPILO_VERSION_MAJOR > 2 || (PAPILO_VERSION_MAJOR == 2 && PAPILO_VERSION_MINOR >= 1)
   DualFix<SCIP_Real> *dualfix = new DualFix<SCIP_Real>();
   dualfix->set_fix_to_infinity_allowed(false);
   presolve.addPresolveMethod( uptr( dualfix ) );
#else
   presolve.addPresolveMethod( uptr( new DualFix<SCIP_Real>() ) );
#endif
   presolve.addPresolveMethod( uptr( new FixContinuous<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SimplifyInequalities<SCIP_Real>() ) );
   presolve.addPresolveMethod( uptr( new SimpleSubstitution<SCIP_Real>() ) );

   /* exhaustive presolvers*/
   presolve.addPresolveMethod( uptr( new ImplIntDetection<SCIP_Real>() ) );
   if( data->enabledualinfer )
      presolve.addPresolveMethod( uptr( new DualInfer<SCIP_Real>() ) );
   if( data->enableprobing )
   {
#if PAPILO_VERSION_MAJOR > 2 || (PAPILO_VERSION_MAJOR == 2 && PAPILO_VERSION_MINOR >= 1)
      Probing<SCIP_Real> *probing = new Probing<SCIP_Real>();
      if( presolve.getPresolveOptions().runs_sequential() )
      {
         probing->set_max_badge_size( data->maxbadgesizeseq );
      }
      else
      {
         probing->set_max_badge_size( data->maxbadgesizepar );
      }
      presolve.addPresolveMethod( uptr( probing ) );

#else
      presolve.addPresolveMethod( uptr( new Probing<SCIP_Real>() ) );
      if( data->maxbadgesizeseq != DEFAULT_MAXBADGESIZE_SEQ )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   The parameter 'presolving/milp/maxbadgesizeseq' can only be used with PaPILO 2.1.0 or later versions.\n");

      if( data->maxbadgesizepar != DEFAULT_MAXBADGESIZE_PAR )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   The parameter 'presolving/milp/maxbadgesizepar' can only be used with PaPILO 2.1.0 or later versions.\n");

#endif
   }
   if( data->enabledomcol )
      presolve.addPresolveMethod( uptr( new DominatedCols<SCIP_Real>() ) );
   if( data->enablemultiaggr )
      presolve.addPresolveMethod( uptr( new Substitution<SCIP_Real>() ) );
   if( data->enablesparsify )
      presolve.addPresolveMethod( uptr( new Sparsify<SCIP_Real>() ) );


   /* set tolerances */
   presolve.getPresolveOptions().feastol = SCIPfeastol(scip);
   presolve.getPresolveOptions().epsilon = SCIPepsilon(scip);

   /* adjust output settings of presolve library */
#ifdef SCIP_PRESOLLIB_ENABLE_OUTPUT
   problem.setName(SCIPgetProbName(scip));
#else
   presolve.setVerbosityLevel(VerbosityLevel::kQuiet);
#endif

   /* communicate the time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      presolve.getPresolveOptions().tlim = timelimit - SCIPgetSolvingTime(scip);

   if( 0 != strncmp(data->filename, DEFAULT_FILENAME_PROBLEM, strlen(DEFAULT_FILENAME_PROBLEM)) )
   {
      SCIP_CALL(SCIPwriteTransProblem(scip, data->filename, NULL, FALSE));
   }

   /* call the presolving */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) running MILP presolver\n", SCIPgetSolvingTime(scip));
   int oldnnz = problem.getConstraintMatrix().getNnz();

   /*call presolving without storing information for dual postsolve*/
#if (PAPILO_VERSION_MAJOR >= 2)
   PresolveResult<SCIP_Real> res = presolve.apply(problem, false);
#else
   PresolveResult<SCIP_Real> res = presolve.apply(problem);
#endif
   data->lastncols = problem.getNCols();
   data->lastnrows = problem.getNRows();

   /* evaluate the result */
   switch(res.status)
   {
      case PresolveStatus::kInfeasible:
         *result = SCIP_CUTOFF;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) MILP presolver detected infeasibility\n",
               SCIPgetSolvingTime(scip));
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::kUnbndOrInfeas:
      case PresolveStatus::kUnbounded:
         *result = SCIP_UNBOUNDED;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) MILP presolver detected unboundedness\n",
               SCIPgetSolvingTime(scip));
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::kUnchanged:
         *result = SCIP_DIDNOTFIND;
         data->lastncols = nvars;
         data->lastnrows = nconss;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) MILP presolver found nothing\n",
               SCIPgetSolvingTime(scip));
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      case PresolveStatus::kReduced:
         data->lastncols = problem.getNCols();
         data->lastnrows = problem.getNRows();
         *result = SCIP_SUCCESS;
   }

   /* result indicated success, now populate the changes into the SCIP structures */
   std::vector<SCIP_VAR*> tmpvars;
   std::vector<SCIP_Real> tmpvals;

   /* if the number of nonzeros decreased by a sufficient factor, rather create all constraints from scratch */
   int newnnz = problem.getConstraintMatrix().getNnz();
   bool constraintsReplaced = false;
   if( newnnz == 0 || (allowconsmodification &&
         (problem.getNRows() <= data->modifyconsfac * data->lastnrows ||
          newnnz <= data->modifyconsfac * oldnnz)) )
   {
      int oldnrows = SCIPmatrixGetNRows(matrix);
      int newnrows = problem.getNRows();

      constraintsReplaced = true;

      /* capture constraints that are still present in the problem after presolve */
      for( int i = 0; i < newnrows; ++i )
      {
         SCIP_CONS* c = SCIPmatrixGetCons(matrix, res.postsolve.origrow_mapping[i]);
         SCIP_CALL( SCIPcaptureCons(scip, c) );
      }

      /* delete all constraints */
      *ndelconss += oldnrows;
      *naddconss += newnrows;

      for( int i = 0; i < oldnrows; ++i )
      {
         SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, i)) );
      }

      /* now loop over rows of presolved problem and create them as new linear constraints,
       * then release the old constraint after its name was passed to the new constraint */
      const Vec<RowFlags>& rflags = problem.getRowFlags();
      const auto& consmatrix = problem.getConstraintMatrix();
      for( int i = 0; i < newnrows; ++i )
      {
         auto rowvec = consmatrix.getRowCoefficients(i);
         const int* rowcols = rowvec.getIndices();
         /* SCIPcreateConsBasicLinear() requires a non const pointer */
         SCIP_Real* rowvals = const_cast<SCIP_Real*>(rowvec.getValues());
         int rowlen = rowvec.getLength();

         /* retrieve SCIP compatible left and right hand sides */
         SCIP_Real lhs = rflags[i].test(RowFlag::kLhsInf) ? - SCIPinfinity(scip) : consmatrix.getLeftHandSides()[i];
         SCIP_Real rhs = rflags[i].test(RowFlag::kRhsInf) ? SCIPinfinity(scip) : consmatrix.getRightHandSides()[i];

         /* create variable array matching the value array */
         tmpvars.clear();
         tmpvars.reserve(rowlen);
         for( int j = 0; j < rowlen; ++j )
            tmpvars.push_back(SCIPmatrixGetVar(matrix, res.postsolve.origcol_mapping[rowcols[j]]));

         /* create and add new constraint with name of old constraint */
         SCIP_CONS* oldcons = SCIPmatrixGetCons(matrix, res.postsolve.origrow_mapping[i]);
         SCIP_CONS* cons;
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, SCIPconsGetName(oldcons), rowlen, tmpvars.data(), rowvals, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, cons) );

         /* release old and new constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &oldcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   /* loop over res.postsolve and add all fixed variables and aggregations to scip */
   for( std::size_t i = 0; i != res.postsolve.types.size(); ++i )
   {
      ReductionType type = res.postsolve.types[i];
      int first = res.postsolve.start[i];
      int last = res.postsolve.start[i + 1];

      switch( type )
      {
      case ReductionType::kFixedCol:
      {
         SCIP_Bool infeas;
         SCIP_Bool fixed;
         int col = res.postsolve.indices[first];

         SCIP_VAR* var = SCIPmatrixGetVar(matrix, col);

         SCIP_Real value = res.postsolve.values[first];

         SCIP_CALL( SCIPfixVar(scip, var, value, &infeas, &fixed) );
         *nfixedvars += 1;

         assert(!infeas);
         /* SCIP has different rules for aggregating variables than PaPILO: therefore the variable PaPILO
          * tries to fix now may have been aggregated by SCIP before. Additionally, after aggregation SCIP
          * sometimes performs bound tightening resulting in possible fixings. These cases need to be excluded. */
         assert(fixed || SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED ||
                      SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED);
         break;
      }
/*
 * Dual-postsolving in PaPILO required introducing a postsolve-type for substitution with additional information.
 * Further, the different Substitution-postsolving types store the required postsolving data differently (in different order) in the postsolving stack.
 * Therefore, we need to distinguish how to parse the required data (rowLength, col, side, startRowCoefficients, lastRowCoefficients) from the postsolving stack.
 * If these values are accessed, the procedure is the same for both.
 */
#if (PAPILO_VERSION_MAJOR >= 2)
      case ReductionType::kSubstitutedColWithDual:
#endif
      case ReductionType::kSubstitutedCol:
      {
         int col = 0;
         SCIP_Real side = 0;

         int rowlen = 0;
         int startRowCoefficients = 0;
         int lastRowCoefficients = 0;

         if( type == ReductionType::kSubstitutedCol )
         {
            rowlen = last - first - 1;
            col = res.postsolve.indices[first];
            side = res.postsolve.values[first];

            startRowCoefficients = first + 1;
            lastRowCoefficients = last;
         }
#if (PAPILO_VERSION_MAJOR >= 2)
         if( type == ReductionType::kSubstitutedColWithDual )
         {
            rowlen = (int) res.postsolve.values[first];
            col = res.postsolve.indices[first + 3 + rowlen];
            side = res.postsolve.values[first + 1];

            startRowCoefficients = first + 3;
            lastRowCoefficients = first + 3 + rowlen;

            assert(side == res.postsolve.values[first + 2]);
            assert(res.postsolve.indices[first + 1] == 0);
            assert(res.postsolve.indices[first + 2] == 0);

         }
         assert( type == ReductionType::kSubstitutedCol || type == ReductionType::kSubstitutedColWithDual );
#else
         assert( type == ReductionType::kSubstitutedCol );
#endif
         SCIP_Bool infeas;
         SCIP_Bool aggregated;
         SCIP_Bool redundant = FALSE;
         SCIP_Real constant = 0.0;
         if( rowlen == 2 )
         {
            SCIP_Real updatedSide;
            SCIP_VAR* varx = SCIPmatrixGetVar(matrix, res.postsolve.indices[startRowCoefficients]);
            SCIP_VAR* vary = SCIPmatrixGetVar(matrix, res.postsolve.indices[startRowCoefficients + 1]);
            SCIP_Real scalarx = res.postsolve.values[startRowCoefficients];
            SCIP_Real scalary = res.postsolve.values[startRowCoefficients + 1];

            SCIP_CALL( SCIPgetProbvarSum(scip, &varx, &scalarx, &constant) );
            assert(SCIPvarGetStatus(varx) != SCIP_VARSTATUS_MULTAGGR);

            SCIP_CALL( SCIPgetProbvarSum(scip, &vary, &scalary, &constant) );
            assert(SCIPvarGetStatus(vary) != SCIP_VARSTATUS_MULTAGGR);

            updatedSide = side - constant;

            SCIP_CALL( SCIPaggregateVars(scip, varx, vary, scalarx, scalary, updatedSide, &infeas, &redundant, &aggregated) );
         }
         else
         {
            SCIP_Real colCoef = 0.0;
            SCIP_Real updatedSide;

            for( int j = startRowCoefficients; j < lastRowCoefficients; ++j )
            {
               if( res.postsolve.indices[j] == col )
               {
                  colCoef = res.postsolve.values[j];
                  break;
               }
            }

            tmpvars.clear();
            tmpvals.clear();
            tmpvars.reserve(rowlen);
            tmpvals.reserve(rowlen);

            assert(colCoef != 0.0);
            SCIP_VAR* aggrvar = SCIPmatrixGetVar(matrix, col);

            SCIP_CALL( SCIPgetProbvarSum(scip, &aggrvar, &colCoef, &constant) );
            assert(SCIPvarGetStatus(aggrvar) != SCIP_VARSTATUS_MULTAGGR);

            updatedSide = side - constant;

            for( int j = startRowCoefficients; j < lastRowCoefficients; ++j )
            {
               if( res.postsolve.indices[j] == col )
                  continue;

               tmpvars.push_back(SCIPmatrixGetVar(matrix, res.postsolve.indices[j]));
               tmpvals.push_back(- res.postsolve.values[j] / colCoef);
            }

            SCIP_CALL( SCIPmultiaggregateVar(scip, aggrvar, tmpvars.size(),
               tmpvars.data(), tmpvals.data(), updatedSide / colCoef, &infeas, &aggregated) );
         }

         if( aggregated )
            *naggrvars += 1;
         else if( constraintsReplaced && !redundant )
         {
            /* if the constraints where replaced, we need to add the failed substitution as an equality to SCIP */
            tmpvars.clear();
            tmpvals.clear();
            for( int j = startRowCoefficients; j < lastRowCoefficients; ++j )
            {
               tmpvars.push_back(SCIPmatrixGetVar(matrix, res.postsolve.indices[j]));
               tmpvals.push_back(res.postsolve.values[j]);
            }

            SCIP_CONS* cons;
            String name = fmt::format("{}_failed_aggregation_equality", SCIPvarGetName(SCIPmatrixGetVar(matrix, col)));
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name.c_str(),
               tmpvars.size(), tmpvars.data(), tmpvals.data(), side, side ) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            *naddconss += 1;
         }

         if( infeas )
         {
            *result = SCIP_CUTOFF;
            break;
         }

         break;
      }
      case ReductionType::kParallelCol:
         return SCIP_INVALIDRESULT;
#if (PAPILO_VERSION_MAJOR <= 1 && PAPILO_VERSION_MINOR==0)
#else
      case ReductionType::kFixedInfCol: {
         /* todo: currently SCIP can not handle this kind of reduction (see issue #3391) */
         assert(false);
         if(!constraintsReplaced)
            continue;
         SCIP_Bool infeas;
         SCIP_Bool fixed;
         SCIP_Real value = SCIPinfinity(scip);

         int column = res.postsolve.indices[first];
         bool is_negative_infinity = res.postsolve.values[first] < 0;
         SCIP_VAR* column_variable = SCIPmatrixGetVar(matrix, column);

         if( is_negative_infinity )
         {
            value = -SCIPinfinity(scip);
         }

         SCIP_CALL( SCIPfixVar(scip, column_variable, value, &infeas, &fixed) );
         *nfixedvars += 1;

         assert(!infeas);
         assert(fixed);
         break;
      }
#endif
#if (PAPILO_VERSION_MAJOR >= 2)
      case ReductionType::kVarBoundChange :
      case ReductionType::kRedundantRow :
      case ReductionType::kRowBoundChange :
      case ReductionType::kReasonForRowBoundChangeForcedByRow :
      case ReductionType::kRowBoundChangeForcedByRow :
      case ReductionType::kSaveRow :
      case ReductionType::kReducedBoundsCost :
      case ReductionType::kColumnDualValue :
      case ReductionType::kRowDualValue :
      case ReductionType::kCoefficientChange :
         // dual ReductionTypes should be only calculated for dual reductions and should not appear for MIP
         SCIPerrorMessage("PaPILO: PaPILO should not return dual postsolving reductions in SCIP!!\n");
         SCIPABORT(); /*lint --e{527}*/
         break;
#endif
      default:
         SCIPdebugMsg(scip, "PaPILO returned unknown data type: \n" );
         continue;
      }
   }

   /* tighten bounds of variables that are still present after presolving */
   if( *result != SCIP_CUTOFF )
   {
      VariableDomains<SCIP_Real>& varDomains = problem.getVariableDomains();
      for( int i = 0; i != problem.getNCols(); ++i )
      {
         assert( ! varDomains.flags[i].test(ColFlag::kInactive) );
         SCIP_VAR* var = SCIPmatrixGetVar(matrix, res.postsolve.origcol_mapping[i]);
         if( !varDomains.flags[i].test(ColFlag::kLbInf) )
         {
            SCIP_Bool infeas;
            SCIP_Bool tightened;
            SCIP_CALL( SCIPtightenVarLb(scip, var, varDomains.lower_bounds[i], TRUE, &infeas, &tightened) );

            if( tightened )
               *nchgbds += 1;

            if( infeas )
            {
               *result = SCIP_CUTOFF;
               break;
            }
         }

         if( !varDomains.flags[i].test(ColFlag::kUbInf) )
         {
            SCIP_Bool infeas;
            SCIP_Bool tightened;
            SCIP_CALL( SCIPtightenVarUb(scip, var, varDomains.upper_bounds[i], TRUE, &infeas, &tightened) );

            if( tightened )
               *nchgbds += 1;

            if( infeas )
            {
               *result = SCIP_CUTOFF;
               break;
            }
         }
      }
   }

   /* finish with a final verb message and return */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "   (%.1fs) MILP presolver (%d rounds): %d aggregations, %d fixings, %d bound changes\n",
      SCIPgetSolvingTime(scip), presolve.getStatistics().nrounds, *naggrvars - oldnaggrvars,
      *nfixedvars - oldnfixedvars, *nchgbds - oldnchgbds);

   /* free the matrix */
   assert(initialized);
   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the xyz presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolMILP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

#if defined(PAPILO_VERSION_TWEAK) && PAPILO_VERSION_TWEAK != 0
   String name = fmt::format("PaPILO {}.{}.{}.{}", PAPILO_VERSION_MAJOR, PAPILO_VERSION_MINOR, PAPILO_VERSION_PATCH, PAPILO_VERSION_TWEAK);
#else
   String name = fmt::format("PaPILO {}.{}.{}", PAPILO_VERSION_MAJOR, PAPILO_VERSION_MINOR, PAPILO_VERSION_PATCH);
#endif

#ifdef PAPILO_GITHASH_AVAILABLE
   String desc = fmt::format("parallel presolve for integer and linear optimization (github.com/scipopt/papilo) [GitHash: {}]", PAPILO_GITHASH);
#else
   String desc("parallel presolve for integer and linear optimization (github.com/scipopt/papilo)");
#endif

   /* add external code info for the presolve library */
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, name.c_str(), desc.c_str()) );

   /* create MILP presolver data */
   presoldata = NULL;
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );
   BMSclearMemory(presoldata);

   presol = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecMILP,
         presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyMILP) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeMILP) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitMILP) );

   /* add MILP presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/threads",
         "maximum number of threads presolving may use (0: automatic)",
         &presoldata->threads, FALSE, DEFAULT_THREADS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/maxfillinpersubstitution",
         "maximal possible fillin for substitutions to be considered",
         &presoldata->maxfillinpersubstitution, FALSE, DEFAULT_MAXFILLINPERSUBST, INT_MIN, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/maxshiftperrow",
         "maximal amount of nonzeros allowed to be shifted to make space for substitutions",
         &presoldata->maxshiftperrow, TRUE, DEFAULT_MAXSHIFTPERROW, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/" PRESOL_NAME "/randomseed",
         "the random seed used for randomization of tie breaking",
         &presoldata->randomseed, FALSE, DEFAULT_RANDOMSEED, INT_MIN, INT_MAX, NULL, NULL) );

   if( DependentRows<double>::Enabled )
   {
      SCIP_CALL( SCIPaddIntParam(scip,
            "presolving/" PRESOL_NAME "/detectlineardependency",
            "should linear dependent equations and free columns be removed? (0: never, 1: for LPs, 2: always)",
            &presoldata->detectlineardependency, TRUE, DEFAULT_DETECTLINDEP, 0, 2, NULL, NULL) );
   }
   else
      presoldata->detectlineardependency = 0;

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/" PRESOL_NAME "/modifyconsfac",
         "modify SCIP constraints when the number of nonzeros or rows is at most this factor "
         "times the number of nonzeros or rows before presolving",
         &presoldata->modifyconsfac, FALSE, DEFAULT_MODIFYCONSFAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/" PRESOL_NAME "/markowitztolerance",
         "the markowitz tolerance used for substitutions",
         &presoldata->markowitztolerance, FALSE, DEFAULT_MARKOWITZTOLERANCE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/" PRESOL_NAME "/hugebound",
         "absolute bound value that is considered too huge for activitity based calculations",
         &presoldata->hugebound, FALSE, DEFAULT_HUGEBOUND, 0.0, SCIP_REAL_MAX, NULL, NULL) );

#if PAPILO_VERSION_MAJOR > 2 || (PAPILO_VERSION_MAJOR == 2 && PAPILO_VERSION_MINOR >= 1)
   SCIP_CALL(SCIPaddIntParam(scip, "presolving/" PRESOL_NAME "/maxbadgesizeseq",
         "maximal badge size in Probing in PaPILO if PaPILO is executed in sequential mode",
         &presoldata->maxbadgesizeseq, FALSE, DEFAULT_MAXBADGESIZE_SEQ, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip, "presolving/" PRESOL_NAME "/maxbadgesizepar",
         "maximal badge size in Probing in PaPILO if PaPILO is executed in parallel mode",
         &presoldata->maxbadgesizepar, FALSE, DEFAULT_MAXBADGESIZE_PAR, -1, INT_MAX, NULL, NULL));
#else
   presoldata->maxbadgesizeseq = DEFAULT_MAXBADGESIZE_SEQ;
   presoldata->maxbadgesizepar = DEFAULT_MAXBADGESIZE_PAR;
#endif

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enableparallelrows",
         "should the parallel rows presolver be enabled within the presolve library?",
         &presoldata->enableparallelrows, TRUE, DEFAULT_ENABLEPARALLELROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enabledomcol",
         "should the dominated column presolver be enabled within the presolve library?",
         &presoldata->enabledomcol, TRUE, DEFAULT_ENABLEDOMCOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enabledualinfer",
         "should the dualinfer presolver be enabled within the presolve library?",
         &presoldata->enabledualinfer, TRUE, DEFAULT_ENABLEDUALINFER, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enablemultiaggr",
         "should the multi-aggregation presolver be enabled within the presolve library?",
         &presoldata->enablemultiaggr, TRUE, DEFAULT_ENABLEMULTIAGGR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enableprobing",
         "should the probing presolver be enabled within the presolve library?",
         &presoldata->enableprobing, TRUE, DEFAULT_ENABLEPROBING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/" PRESOL_NAME "/enablesparsify",
         "should the sparsify presolver be enabled within the presolve library?",
         &presoldata->enablesparsify, TRUE, DEFAULT_ENABLESPARSIFY, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "presolving/" PRESOL_NAME "/probfilename",
         "filename to store the problem before MILP presolving starts",
         &presoldata->filename, TRUE, DEFAULT_FILENAME_PROBLEM, NULL, NULL) );

   return SCIP_OKAY;
}

#endif
