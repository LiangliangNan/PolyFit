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

/**@file   lpi_glop.cpp
 * @ingroup LPIS
 * @brief  LP interface for Glop
 * @author Frederic Didier
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* turn off some warnings from includes */
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#pragma GCC diagnostic ignored "-Woverflow"

#include "ortools/base/version.h"
#include "ortools/glop/lp_solver.h"
#include "ortools/glop/revised_simplex.h"
#include "ortools/lp_data/lp_print_utils.h"
#include "ortools/lp_data/lp_data_utils.h"
#include "ortools/lp_data/proto_utils.h"
#include "ortools/util/file_util.h"
#include "ortools/util/stats.h"
#include "ortools/util/time_limit.h"

#include "ortools/base/logging.h"
#include "ortools/base/vlog_is_on.h"

#include "lpi/lpi.h"
#include "scip/pub_message.h"

#include <assert.h>

/* turn warnings on again */
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wpedantic"
#pragma GCC diagnostic warning "-Wignored-qualifiers"
#pragma GCC diagnostic warning "-Wshadow"
#pragma GCC diagnostic warning "-Wnon-virtual-dtor"
#pragma GCC diagnostic warning "-Wctor-dtor-privacy"
#pragma GCC diagnostic warning "-Woverflow"


using operations_research::TimeLimit;
using operations_research::glop::BasisState;
using operations_research::glop::ColIndex;
using operations_research::glop::ColIndexVector;
using operations_research::glop::ConstraintStatus;
using operations_research::glop::ConstraintStatusColumn;
using operations_research::glop::DenseBooleanColumn;
using operations_research::glop::DenseBooleanRow;
using operations_research::glop::DenseColumn;
using operations_research::glop::DenseRow;
using operations_research::glop::SparseColumn;
using operations_research::glop::ScatteredColumn;
using operations_research::glop::ScatteredColumnIterator;
using operations_research::glop::SparseMatrix;
using operations_research::glop::Fractional;
using operations_research::glop::GetProblemStatusString;
using operations_research::glop::ProblemStatus;
using operations_research::glop::RowIndex;
using operations_research::glop::ScatteredRow;
using operations_research::glop::ScatteredRowIterator;
using operations_research::glop::VariableStatus;
using operations_research::glop::VariableStatusRow;
using operations_research::MPModelProto;

/** LP interface */
struct SCIP_LPi
{
   operations_research::glop::LinearProgram*   linear_program;     /**< the linear program */
   operations_research::glop::LinearProgram*   scaled_lp;          /**< scaled linear program */
   operations_research::glop::RevisedSimplex*  solver;             /**< direct reference to the revised simplex, not passing through lp_solver */
   operations_research::glop::GlopParameters*  parameters;         /**< parameters */
   operations_research::glop::LpScalingHelper* scaler;             /**< scaler auxiliary class */

   /* the following is used by SCIPlpiWasSolved() */
   bool                  lp_modified_since_last_solve;
   bool                  lp_time_limit_was_reached;

   /* store the values of some parameters in order to be able to return them */
   bool                  lp_info;            /**< whether additional output is turned on */
   SCIP_PRICING          pricing;            /**< SCIP pricing setting  */
   bool                  from_scratch;       /**< store whether basis is ignored for next solving call */
   int                   numthreads;         /**< number of threads used to solve the LP (0 = automatic) */
   SCIP_Real             conditionlimit;     /**< maximum condition number of LP basis counted as stable (-1.0: no limit) */
   bool                  checkcondition;     /**< Should condition number of LP basis be checked for stability? */
   int                   timing;             /**< type of timer (1 - cpu, 2 - wallclock, 0 - off) */

   /* other data */
   SCIP_Longint          niterations;        /**< number of iterations used */

   /* Temporary vectors allocated here for speed. This gain is non-negligible
    * because in many situations, only a few entries of these vectors are
    * inspected (hypersparsity) and allocating them is in O(num_rows) or
    * O(num_cols) instead of O(num_non_zeros) to read/clear them. */
   ScatteredRow*         tmp_row;            /**< temporary vector */
   ScatteredColumn*      tmp_column;         /**< temporary vector */
};

/** uncomment to turn off scaling */
/* #define NOSCALING */

/** define feasibility check to possibly reoptimize: 0: no check, 1: completely new check, 2: check unscaled variable and activity values */
#ifdef NOSCALING
#define UNSCALEDFEAS_CHECK 0
#else
#define UNSCALEDFEAS_CHECK 2
#endif


/*
 * LP Interface Methods
 */

char* initGlopName( );

static char* glopname = initGlopName( );

char* initGlopName( )
{
   glopname = new char[100];
   (void) snprintf(glopname, 100, "Glop %d.%d", operations_research::OrToolsMajorVersion(), operations_research::OrToolsMinorVersion());
   return glopname;
}

/** gets name and version of LP solver */
const char* SCIPlpiGetSolverName(
   void
   )
{
   return glopname;
}

/** gets description of LP solver (developer, webpage, ...) */
const char* SCIPlpiGetSolverDesc(
   void
   )
{
   return "Glop Linear Solver, developed by Google, part of OR-Tools (developers.google.com/optimization)";
}

/** gets pointer for LP solver - use only with great care */
void* SCIPlpiGetSolverPointer(
   SCIP_LPI*             lpi                 /**< pointer to an LP interface structure */
   )
{
   assert( lpi != NULL );

   SCIPerrorMessage("SCIPlpiGetSolverPointer() has not been implemented yet.\n");

   return NULL;
}

/** pass integrality information to LP solver */
SCIP_RETCODE SCIPlpiSetIntegralityInformation(
   SCIP_LPI*             lpi,                /**< pointer to an LP interface structure */
   int                   ncols,              /**< length of integrality array */
   int*                  intInfo             /**< integrality array (0: continuous, 1: integer). May be NULL iff ncols is 0.  */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( ncols == 0 || ncols == lpi->linear_program->num_variables().value() );

   /* Pass on integrality information (currently not used by Glop). */
   for (ColIndex col(0); col < ColIndex(ncols); ++col)
   {
      assert( intInfo != NULL );
      int info = intInfo[col.value()];
      assert( info == 0 || info == 1 );
      if ( info == 0 )
         lpi->linear_program->SetVariableType(col, operations_research::glop::LinearProgram::VariableType::CONTINUOUS);
      else
         lpi->linear_program->SetVariableType(col, operations_research::glop::LinearProgram::VariableType::INTEGER);
   }

   return SCIP_OKAY;
}

/** informs about availability of a primal simplex solving method */
SCIP_Bool SCIPlpiHasPrimalSolve(
   void
   )
{
   return TRUE;
}

/** informs about availability of a dual simplex solving method */
SCIP_Bool SCIPlpiHasDualSolve(
   void
   )
{
   return TRUE;
}

/** informs about availability of a barrier solving method */
SCIP_Bool SCIPlpiHasBarrierSolve(
   void
   )
{
   return FALSE;
}


/*
 * LPI Creation and Destruction Methods
 */

/**@name LPI Creation and Destruction Methods */
/**@{ */

/** creates an LP problem object */
SCIP_RETCODE SCIPlpiCreate(
   SCIP_LPI**            lpi,                /**< pointer to an LP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert( lpi != NULL );
   assert( name != NULL );

   /* Initilialize memory. */
   SCIP_ALLOC(BMSallocMemory(lpi));
   (*lpi)->linear_program = new operations_research::glop::LinearProgram();
   (*lpi)->scaled_lp = new operations_research::glop::LinearProgram();
   (*lpi)->solver = new operations_research::glop::RevisedSimplex();
   (*lpi)->parameters = new operations_research::glop::GlopParameters();
   (*lpi)->scaler = new operations_research::glop::LpScalingHelper();

   /* Set problem name and objective direction. */
   (*lpi)->linear_program->SetName(std::string(name));
   SCIP_CALL( SCIPlpiChgObjsen(*lpi, objsen) );

   (*lpi)->from_scratch = false;
   (*lpi)->lp_info = false;
   (*lpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*lpi)->lp_modified_since_last_solve = true;
   (*lpi)->lp_time_limit_was_reached = false;
   (*lpi)->conditionlimit = -1.0;
   (*lpi)->checkcondition = false;
   (*lpi)->niterations = 0LL;

   (*lpi)->tmp_row = new ScatteredRow();
   (*lpi)->tmp_column = new ScatteredColumn();

#ifdef NOSCALING
   (*lpi)->parameters->set_use_scaling(false);
#endif

   return SCIP_OKAY;
}

/** deletes an LP problem object */
SCIP_RETCODE SCIPlpiFree(
   SCIP_LPI**            lpi                 /**< pointer to an LP interface structure */
   )
{
   SCIPdebugMessage("SCIPlpiFree\n");

   delete (*lpi)->scaler;
   delete (*lpi)->parameters;
   delete (*lpi)->solver;
   delete (*lpi)->scaled_lp;
   delete (*lpi)->linear_program;

   delete (*lpi)->tmp_row;
   delete (*lpi)->tmp_column;

   BMSfreeMemory(lpi);

   return SCIP_OKAY;
}

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies LP data with column matrix into LP solver */
SCIP_RETCODE SCIPlpiLoadColLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( beg != NULL );
   assert( ind != NULL );
   assert( val != NULL );

   lpi->linear_program->Clear();
   SCIP_CALL( SCIPlpiAddRows(lpi, nrows, lhs, rhs, rownames, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, colnames, nnonz, beg, ind, val) );
   SCIP_CALL( SCIPlpiChgObjsen(lpi, objsen) );

   return SCIP_OKAY;
}

/** adds columns to the LP */
SCIP_RETCODE SCIPlpiAddCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( nnonz >= 0) ;
   assert( ncols >= 0) ;

   SCIPdebugMessage("adding %d columns with %d nonzeros.\n", ncols, nnonz);

   /* @todo add names */
   if ( nnonz > 0 )
   {
      assert( beg != NULL );
      assert( ind != NULL );
      assert( val != NULL );
      assert( ncols > 0 );

#ifndef NDEBUG
      /* perform check that no new rows are added */
      RowIndex num_rows = lpi->linear_program->num_constraints();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( 0 <= ind[j] && ind[j] < num_rows.value() );
         assert( val[j] != 0.0 );
      }
#endif

      int nz = 0;
      for (int i = 0; i < ncols; ++i)
      {
         const ColIndex col = lpi->linear_program->CreateNewVariable();
         lpi->linear_program->SetVariableBounds(col, lb[i], ub[i]);
         lpi->linear_program->SetObjectiveCoefficient(col, obj[i]);
         const int end = (nnonz == 0 || i == ncols - 1) ? nnonz : beg[i + 1];
         while ( nz < end )
         {
            lpi->linear_program->SetCoefficient(RowIndex(ind[nz]), col, val[nz]);
            ++nz;
         }
      }
      assert( nz == nnonz );
   }
   else
   {
      for (int i = 0; i < ncols; ++i)
      {
         const ColIndex col = lpi->linear_program->CreateNewVariable();
         lpi->linear_program->SetVariableBounds(col, lb[i], ub[i]);
         lpi->linear_program->SetObjectiveCoefficient(col, obj[i]);
      }
   }

   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** deletes all columns in the given range from LP */
SCIP_RETCODE SCIPlpiDelCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < lpi->linear_program->num_variables() );

   SCIPdebugMessage("deleting columns %d to %d.\n", firstcol, lastcol);

   const ColIndex num_cols = lpi->linear_program->num_variables();
   DenseBooleanRow columns_to_delete(num_cols, false);
   for (int i = firstcol; i <= lastcol; ++i)
      columns_to_delete[ColIndex(i)] = true;

   lpi->linear_program->DeleteColumns(columns_to_delete);
   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** deletes columns from SCIP_LP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelColset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( dstat != NULL );

   const ColIndex num_cols = lpi->linear_program->num_variables();
   DenseBooleanRow columns_to_delete(num_cols, false);
   int new_index = 0;
   int num_deleted_columns = 0;
   for (ColIndex col(0); col < num_cols; ++col)
   {
      int i = col.value();
      if ( dstat[i] == 1 )
      {
         columns_to_delete[col] = true;
         dstat[i] = -1;
         ++num_deleted_columns;
      }
      else
         dstat[i] = new_index++;
   }
   SCIPdebugMessage("SCIPlpiDelColset: deleting %d columns.\n", num_deleted_columns);
   lpi->linear_program->DeleteColumns(columns_to_delete);
   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** adds rows to the LP */
SCIP_RETCODE SCIPlpiAddRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lhs != NULL );
   assert( rhs != NULL );
   assert( nnonz >= 0) ;
   assert( nrows >= 0) ;

   SCIPdebugMessage("adding %d rows with %d nonzeros.\n", nrows, nnonz);

   /* @todo add names */
   if ( nnonz > 0 )
   {
      assert( beg != NULL );
      assert( ind != NULL );
      assert( val != NULL );
      assert( nrows > 0 );

#ifndef NDEBUG
      /* perform check that no new columns are added - this is likely to be a mistake */
      const ColIndex num_cols = lpi->linear_program->num_variables();
      for (int j = 0; j < nnonz; ++j)
      {
         assert( val[j] != 0.0 );
         assert( 0 <= ind[j] && ind[j] < num_cols.value() );
      }
#endif

      int nz = 0;
      for (int i = 0; i < nrows; ++i)
      {
         const RowIndex row = lpi->linear_program->CreateNewConstraint();
         lpi->linear_program->SetConstraintBounds(row, lhs[i], rhs[i]);
         const int end = (nnonz == 0 || i == nrows - 1) ? nnonz : beg[i + 1];
         while ( nz < end )
         {
            lpi->linear_program->SetCoefficient(row, ColIndex(ind[nz]), val[nz]);
            ++nz;
         }
      }
      assert( nz == nnonz );
   }
   else
   {
      for (int i = 0; i < nrows; ++i)
      {
         const RowIndex row = lpi->linear_program->CreateNewConstraint();
         lpi->linear_program->SetConstraintBounds(row, lhs[i], rhs[i]);
      }
   }

   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** delete rows from LP and update the current basis */
static
void deleteRowsAndUpdateCurrentBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const DenseBooleanColumn& rows_to_delete  /**< array to mark rows that should be deleted */
   )
{
   const RowIndex num_rows = lpi->linear_program->num_constraints();
   const ColIndex num_cols = lpi->linear_program->num_variables();

   /* try to repair basis status if problem size has not changed before */
   BasisState state = lpi->solver->GetState();
   if ( state.statuses.size() == num_cols.value() + num_rows.value() )
   {
      /* Shift the status of the non-deleted rows. Note that if the deleted rows where part of the basis (i.e., constraint
       * not tight), then we should be left with a correct basis afterward. This should be the most common use case in SCIP. */
      ColIndex new_size  = num_cols;
      for (RowIndex row(0); row < num_rows; ++row)
      {
         if ( rows_to_delete[row] )
            continue;
         state.statuses[new_size++] = state.statuses[num_cols + RowToColIndex(row)];
      }
      state.statuses.resize(new_size);
      lpi->solver->LoadStateForNextSolve(state);
   }

   lpi->linear_program->DeleteRows(rows_to_delete);
   lpi->lp_modified_since_last_solve = true;
}

/** deletes all rows in the given range from LP */
SCIP_RETCODE SCIPlpiDelRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( 0 <= firstrow && firstrow <= lastrow && lastrow < lpi->linear_program->num_constraints() );

   const RowIndex num_rows = lpi->linear_program->num_constraints();
   DenseBooleanColumn rows_to_delete(num_rows, false);
   for (int i = firstrow; i <= lastrow; ++i)
      rows_to_delete[RowIndex(i)] = true;

   SCIPdebugMessage("deleting rows %d to %d.\n", firstrow, lastrow);
   deleteRowsAndUpdateCurrentBasis(lpi, rows_to_delete);

   return SCIP_OKAY;
}

/** deletes rows from SCIP_LP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPlpiDelRowset(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   const RowIndex num_rows = lpi->linear_program->num_constraints();
   DenseBooleanColumn rows_to_delete(num_rows, false);
   int new_index = 0;
   int num_deleted_rows = 0;
   for (RowIndex row(0); row < num_rows; ++row)
   {
      int i = row.value();
      if (dstat[i] == 1)
      {
         rows_to_delete[row] = true;
         dstat[i] = -1;
         ++num_deleted_rows;
      }
      else
         dstat[i] = new_index++;
   }

   SCIPdebugMessage("SCIPlpiDelRowset: deleting %d rows.\n", num_deleted_rows);
   deleteRowsAndUpdateCurrentBasis(lpi, rows_to_delete);

   return SCIP_OKAY;
}

/** clears the whole LP */
SCIP_RETCODE SCIPlpiClear(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   SCIPdebugMessage("SCIPlpiClear\n");

   lpi->linear_program->Clear();
   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPlpiChgBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices or NULL if ncols is zero */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds or NULL if ncols is zero */
   const SCIP_Real*      ub                  /**< values for the new upper bounds or NULL if ncols is zero */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( ncols == 0 || (ind != NULL && lb != NULL && ub != NULL) );

   SCIPdebugMessage("changing %d bounds.\n", ncols);
   if ( ncols <= 0 )
      return SCIP_OKAY;

   for (int i = 0; i < ncols; ++i)
   {
      SCIPdebugMessage("  col %d: [%g,%g]\n", ind[i], lb[i], ub[i]);

      if ( SCIPlpiIsInfinity(lpi, lb[i]) )
      {
         SCIPerrorMessage("LP Error: fixing lower bound for variable %d to infinity.\n", ind[i]);
         return SCIP_LPERROR;
      }
      if ( SCIPlpiIsInfinity(lpi, -ub[i]) )
      {
         SCIPerrorMessage("LP Error: fixing upper bound for variable %d to -infinity.\n", ind[i]);
         return SCIP_LPERROR;
      }

      lpi->linear_program->SetVariableBounds(ColIndex(ind[i]), lb[i], ub[i]);
   }
   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPlpiChgSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   if ( nrows <= 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("changing %d sides\n", nrows);

   for (int i = 0; i < nrows; ++i)
      lpi->linear_program->SetConstraintBounds(RowIndex(ind[i]), lhs[i], rhs[i]);

   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPlpiChgCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( 0 <= row && row < lpi->linear_program->num_constraints().value() );
   assert( 0 <= col && col < lpi->linear_program->num_variables().value() );

   SCIPdebugMessage("Set coefficient (%d,%d) to %f.\n", row, col, newval);
   lpi->linear_program->CleanUp();
   assert( lpi->linear_program->IsCleanedUp() );
   lpi->linear_program->SetCoefficient(RowIndex(row), ColIndex(col), newval);

   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPlpiChgObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   assert( lpi != NULL);
   assert( lpi->linear_program != NULL);

   SCIPdebugMessage("changing objective sense to %d\n", objsen);

   switch (objsen)
   {
   case SCIP_OBJSEN_MAXIMIZE:
      lpi->linear_program->SetMaximizationProblem(true);
      break;
   case SCIP_OBJSEN_MINIMIZE:
      lpi->linear_program->SetMaximizationProblem(false);
      break;
   }
   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** changes objective values of columns in the LP */
SCIP_RETCODE SCIPlpiChgObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   const int*            ind,                /**< column indices to change objective value for */
   const SCIP_Real*      obj                 /**< new objective values for columns */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( ind != NULL );
   assert( obj != NULL );

   SCIPdebugMessage("changing %d objective values\n", ncols);

   for (int i = 0; i < ncols; ++i)
      lpi->linear_program->SetObjectiveCoefficient(ColIndex(ind[i]), obj[i]);

   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPlpiScaleRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nnonz;
   int* inds;
   int beg;

   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   SCIPdebugMessage("Scale row %d by %f.\n", row, scaleval);

   /* alloc memory */
   ColIndex num_cols = lpi->linear_program->num_variables();

   SCIP_ALLOC( BMSallocMemoryArray(&inds, num_cols.value()) );
   SCIP_ALLOC( BMSallocMemoryArray(&vals, num_cols.value()) );

   /* get the row */
   SCIP_CALL( SCIPlpiGetRows(lpi, row, row, &lhs, &rhs, &nnonz, &beg, inds, vals) );

   /* scale row coefficients */
   for (int j = 0; j < nnonz; ++j)
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, row, inds[j], vals[j] * scaleval) );
   }
   BMSfreeMemoryArray(&vals);
   BMSfreeMemoryArray(&inds);

   /* scale row sides */
   if ( ! SCIPlpiIsInfinity(lpi, -lhs) )
      lhs *= scaleval;
   else if ( scaleval < 0.0 )
      lhs = SCIPlpiInfinity(lpi);

   if ( ! SCIPlpiIsInfinity(lpi, rhs) )
      rhs *= scaleval;
   else if ( scaleval < 0.0 )
      rhs = -SCIPlpiInfinity(lpi);

   if ( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &lhs, &rhs) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgSides(lpi, 1, &row, &rhs, &lhs) );
   }

   lpi->lp_modified_since_last_solve = true;

   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPlpiScaleCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   SCIP_Real* vals;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   int nnonz;
   int* inds;
   int beg;

   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );

   SCIPdebugMessage("Scale column %d by %f.\n", col, scaleval);

   /* alloc memory */
   RowIndex num_rows = lpi->linear_program->num_constraints();

   SCIP_ALLOC( BMSallocMemoryArray(&inds, num_rows.value()) );
   SCIP_ALLOC( BMSallocMemoryArray(&vals, num_rows.value()) );

   /* get the column */
   SCIP_CALL( SCIPlpiGetCols(lpi, col, col, &lb, &ub, &nnonz, &beg, inds, vals) );

   /* scale column coefficients */
   for (int j = 0; j < nnonz; ++j)
   {
      SCIP_CALL( SCIPlpiChgCoef(lpi, col, inds[j], vals[j] * scaleval) );
   }
   BMSfreeMemoryArray(&vals);
   BMSfreeMemoryArray(&inds);

   /* scale objective value */
   SCIP_CALL( SCIPlpiGetObj(lpi, col, col, &obj) );
   obj *= scaleval;
   SCIP_CALL( SCIPlpiChgObj(lpi, 1, &col, &obj) );

   /* scale bound */
   if ( ! SCIPlpiIsInfinity(lpi, -lb) )
      lb *= scaleval;
   else if ( scaleval < 0.0 )
      lb = SCIPlpiInfinity(lpi);

   if ( ! SCIPlpiIsInfinity(lpi, ub) )
      ub *= scaleval;
   else if ( scaleval < 0.0 )
      ub = -SCIPlpiInfinity(lpi);

   if ( scaleval > 0.0 )
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &lb, &ub) );
   }
   else
   {
      SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &col, &ub, &lb) );
   }

   return SCIP_OKAY;
}


/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the LP */
SCIP_RETCODE SCIPlpiGetNRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( nrows != NULL );

   SCIPdebugMessage("getting number of rows.\n");

   *nrows = lpi->linear_program->num_constraints().value();

   return SCIP_OKAY;
}

/** gets the number of columns in the LP */
SCIP_RETCODE SCIPlpiGetNCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( ncols != NULL );

   SCIPdebugMessage("getting number of columns.\n");

   *ncols = lpi->linear_program->num_variables().value();

   return SCIP_OKAY;
}

/** gets objective sense of the LP */
SCIP_RETCODE SCIPlpiGetObjsen(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( objsen != NULL );

   SCIPdebugMessage("getting objective sense.\n");

   *objsen = lpi->linear_program->IsMaximizationProblem() ? SCIP_OBJSEN_MAXIMIZE : SCIP_OBJSEN_MINIMIZE;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP constraint matrix */
SCIP_RETCODE SCIPlpiGetNNonz(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( nnonz != NULL );

   SCIPdebugMessage("getting number of non-zeros.\n");

   *nnonz = (int) lpi->linear_program->num_entries().value();

   return SCIP_OKAY;
}

/** gets columns from LP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetCols(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get from LP */
   int                   lastcol,            /**< last column to get from LP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < lpi->linear_program->num_variables() );
   assert( (lb != NULL && ub != NULL) || (lb == NULL && ub == NULL) );
   assert( (nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL) );

   const DenseRow& tmplb = lpi->linear_program->variable_lower_bounds();
   const DenseRow& tmpub = lpi->linear_program->variable_upper_bounds();

   if ( nnonz != NULL )
   {
      assert( beg != NULL );
      assert( ind != NULL );
      assert( val != NULL );

      *nnonz = 0;
      int index = 0;
      for (ColIndex col(firstcol); col <= ColIndex(lastcol); ++col, ++index)
      {
         if ( lb != NULL )
            lb[index] = tmplb[col];
         if ( ub != NULL )
            ub[index] = tmpub[col];

         beg[index] = *nnonz;
         const SparseColumn& column = lpi->linear_program->GetSparseColumn(col);
         for (const SparseColumn::Entry& entry : column)
         {
            const RowIndex row = entry.row();
            ind[*nnonz] = row.value();
            val[*nnonz] = entry.coefficient();
            ++(*nnonz);
         }
      }
   }
   else
   {
      int index = 0;
      for (ColIndex col(firstcol); col <= ColIndex(lastcol); ++col, ++index)
      {
         if ( lb != NULL )
            lb[index] = tmplb[col];
         if ( ub != NULL )
            ub[index] = tmpub[col];
      }
   }

   return SCIP_OKAY;
}

/** gets rows from LP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPlpiGetRows(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get from LP */
   int                   lastrow,            /**< last row to get from LP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( 0 <= firstrow && firstrow <= lastrow && lastrow < lpi->linear_program->num_constraints() );
   assert( (lhs == NULL && rhs == NULL) || (rhs != NULL && lhs != NULL) );
   assert( (nnonz != NULL && beg != NULL && ind != NULL && val != NULL) || (nnonz == NULL && beg == NULL && ind == NULL && val == NULL) );

   const DenseColumn& tmplhs = lpi->linear_program->constraint_lower_bounds();
   const DenseColumn& tmprhs = lpi->linear_program->constraint_upper_bounds();

   if ( nnonz != NULL )
   {
      assert( beg != NULL );
      assert( ind != NULL );
      assert( val != NULL );

      const SparseMatrix& matrixtrans = lpi->linear_program->GetTransposeSparseMatrix();

      *nnonz = 0;
      int index = 0;
      for (RowIndex row(firstrow); row <= RowIndex(lastrow); ++row, ++index)
      {
         if ( lhs != NULL )
            lhs[index] = tmplhs[row];
         if ( rhs != NULL )
            rhs[index] = tmprhs[row];

         beg[index] = *nnonz;
         const SparseColumn& column = matrixtrans.column(ColIndex(row.value()));
         for (const SparseColumn::Entry& entry : column)
         {
            const RowIndex rowidx = entry.row();
            ind[*nnonz] = rowidx.value();
            val[*nnonz] = entry.coefficient();
            ++(*nnonz);
         }
      }
   }
   else
   {
      int index = 0;
      for (RowIndex row(firstrow); row <= RowIndex(lastrow); ++row, ++index)
      {
         if ( lhs != NULL )
            lhs[index] = tmplhs[row];
         if ( rhs != NULL )
            rhs[index] = tmprhs[row];
      }
   }

   return SCIP_OKAY;
}

/** gets column names */
SCIP_RETCODE SCIPlpiGetColNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get name from LP */
   int                   lastcol,            /**< last column to get name from LP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for col names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( colnames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstcol && firstcol <= lastcol && lastcol < lpi->linear_program->num_variables() );

   SCIPerrorMessage("SCIPlpiGetColNames() has not been implemented yet.\n");

   return SCIP_NOTIMPLEMENTED;
}

/** gets row names */
SCIP_RETCODE SCIPlpiGetRowNames(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get name from LP */
   int                   lastrow,            /**< last row to get name from LP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) or NULL if namestoragesize is zero */
   char*                 namestorage,        /**< storage for row names or NULL if namestoragesize is zero */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) or NULL if namestoragesize is zero */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( rownames != NULL || namestoragesize == 0 );
   assert( namestorage != NULL || namestoragesize == 0 );
   assert( namestoragesize >= 0 );
   assert( storageleft != NULL );
   assert( 0 <= firstrow && firstrow <= lastrow && lastrow < lpi->linear_program->num_constraints() );

   SCIPerrorMessage("SCIPlpiGetRowNames() has not been implemented yet.\n");

   return SCIP_NOTIMPLEMENTED;
}

/** gets objective coefficients from LP problem object */
SCIP_RETCODE SCIPlpiGetObj(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( firstcol <= lastcol );
   assert( vals != NULL );

   SCIPdebugMessage("getting objective values %d to %d\n", firstcol, lastcol);

   int index = 0;
   for (ColIndex col(firstcol); col <= ColIndex(lastcol); ++col)
   {
      vals[index] = lpi->linear_program->objective_coefficients()[col];
      ++index;
   }

   return SCIP_OKAY;
}

/** gets current bounds from LP problem object */
SCIP_RETCODE SCIPlpiGetBounds(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( firstcol <= lastcol );

   SCIPdebugMessage("getting bounds %d to %d\n", firstcol, lastcol);

   int index = 0;
   for (ColIndex col(firstcol); col <= ColIndex(lastcol); ++col)
   {
      if ( lbs != NULL )
         lbs[index] = lpi->linear_program->variable_lower_bounds()[col];

      if ( ubs != NULL )
         ubs[index] = lpi->linear_program->variable_upper_bounds()[col];

      ++index;
   }

   return SCIP_OKAY;
}

/** gets current row sides from LP problem object */
SCIP_RETCODE SCIPlpiGetSides(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( firstrow <= lastrow );

   SCIPdebugMessage("getting row sides %d to %d\n", firstrow, lastrow);

   int index = 0;
   for (RowIndex row(firstrow); row <= RowIndex(lastrow); ++row)
   {
      if ( lhss != NULL )
         lhss[index] = lpi->linear_program->constraint_lower_bounds()[row];

      if ( rhss != NULL )
         rhss[index] = lpi->linear_program->constraint_upper_bounds()[row];

      ++index;
   }

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPlpiGetCoef(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( val != NULL );

   /* quite slow method: possibly needs linear time if matrix is not sorted */
   const SparseMatrix& matrix = lpi->linear_program->GetSparseMatrix();
   *val = matrix.LookUpValue(RowIndex(row), ColIndex(col));

   return SCIP_OKAY;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** update scaled linear program */
static
void updateScaledLP(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   if ( ! lpi->lp_modified_since_last_solve )
      return;

   lpi->scaled_lp->PopulateFromLinearProgram(*lpi->linear_program);
   lpi->scaled_lp->AddSlackVariablesWhereNecessary(false);

   /* @todo: Avoid doing a copy if there is no scaling. */
   /* @todo: Avoid rescaling if not much changed. */
   if ( lpi->parameters->use_scaling() )
      lpi->scaler->Scale(lpi->scaled_lp);
   else
      lpi->scaler->Clear();
}

/** check primal feasibility */
static
bool checkUnscaledPrimalFeasibility(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->linear_program != NULL );

#if UNSCALEDFEAS_CHECK == 1
   /* get unscaled solution */
   const ColIndex num_cols = lpi->linear_program->num_variables();
   DenseRow unscaledsol(num_cols);
   for (ColIndex col = ColIndex(0); col < num_cols; ++col)
      unscaledsol[col] = lpi->scaler->UnscaleVariableValue(col, lpi->solver->GetVariableValue(col));

   /* if the solution is not feasible w.r.t. absolute tolerances, try to fix it in the unscaled problem */
   const double feastol = lpi->parameters->primal_feasibility_tolerance();
   return lpi->linear_program->SolutionIsLPFeasible(unscaledsol, feastol);

#elif UNSCALEDFEAS_CHECK == 2
   const double feastol = lpi->parameters->primal_feasibility_tolerance();

   /* check bounds of unscaled solution */
   const ColIndex num_cols = lpi->linear_program->num_variables();
   for (ColIndex col = ColIndex(0); col < num_cols; ++col)
   {
      const Fractional val = lpi->scaler->UnscaleVariableValue(col, lpi->solver->GetVariableValue(col));
      const Fractional lb = lpi->linear_program->variable_lower_bounds()[col];
      if ( val < lb - feastol )
         return false;
      const Fractional ub = lpi->linear_program->variable_upper_bounds()[col];
      if ( val > ub + feastol )
         return false;
   }

   /* check activities of unscaled solution */
   const RowIndex num_rows = lpi->linear_program->num_constraints();
   for (RowIndex row(0); row < num_rows; ++row)
   {
      const Fractional val = lpi->scaler->UnscaleConstraintActivity(row, lpi->solver->GetConstraintActivity(row));
      const Fractional lhs = lpi->linear_program->constraint_lower_bounds()[row];
      if ( val < lhs - feastol )
         return false;
      const Fractional rhs = lpi->linear_program->constraint_upper_bounds()[row];
      if ( val > rhs + feastol )
         return false;
   }
#endif

   return true;
}

/** common function between the two LPI Solve() functions */
static
SCIP_RETCODE SolveInternal(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   bool                  recursive,          /**< Is this a recursive call? */
   std::unique_ptr<TimeLimit>& time_limit    /**< time limit */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->parameters != NULL );

   updateScaledLP(lpi);

   lpi->solver->SetParameters(*(lpi->parameters));
   lpi->lp_time_limit_was_reached = false;

   /* possibly ignore warm start information for next solve */
   if ( lpi->from_scratch )
      lpi->solver->ClearStateForNextSolve();

   if ( ! lpi->solver->Solve(*(lpi->scaled_lp), time_limit.get()).ok() )
   {
      return SCIP_LPERROR;
   }
   lpi->lp_time_limit_was_reached = time_limit->LimitReached();
   if ( recursive )
      lpi->niterations += (SCIP_Longint) lpi->solver->GetNumberOfIterations();
   else
      lpi->niterations = (SCIP_Longint) lpi->solver->GetNumberOfIterations();

   SCIPdebugMessage("status=%s  obj=%f  iter=%ld.\n", GetProblemStatusString(lpi->solver->GetProblemStatus()).c_str(),
      lpi->solver->GetObjectiveValue(), lpi->solver->GetNumberOfIterations());

   const ProblemStatus status = lpi->solver->GetProblemStatus();
   if ( (status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::OPTIMAL) && lpi->parameters->use_scaling() )
   {
      if ( ! checkUnscaledPrimalFeasibility(lpi) )
      {
         SCIPdebugMessage("Solution not feasible w.r.t. absolute tolerance %g -> reoptimize.\n", lpi->parameters->primal_feasibility_tolerance());

         /* Re-solve without scaling to try to fix the infeasibility. */
         lpi->parameters->set_use_scaling(false);
         lpi->lp_modified_since_last_solve = true;
         SolveInternal(lpi, true, time_limit);   /* inherit time limit, so used time is not reset; do not change iteration limit for resolve */
         lpi->parameters->set_use_scaling(true);

#ifdef SCIP_DEBUG
         if ( ! checkUnscaledPrimalFeasibility(lpi) )
            SCIPdebugMessage("Solution still not feasible after turning off scaling.\n");
#endif
      }
   }

   lpi->lp_modified_since_last_solve = false;

   return SCIP_OKAY;
}

/** calls primal simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolvePrimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->parameters != NULL );

   SCIPdebugMessage("SCIPlpiSolvePrimal: %d rows, %d cols.\n", lpi->linear_program->num_constraints().value(), lpi->linear_program->num_variables().value());
   std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(*lpi->parameters);
   lpi->niterations = 0;

   lpi->parameters->set_use_dual_simplex(false);
   return SolveInternal(lpi, false, time_limit);
}

/** calls dual simplex to solve the LP */
SCIP_RETCODE SCIPlpiSolveDual(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->parameters != NULL );

   SCIPdebugMessage("SCIPlpiSolveDual: %d rows, %d cols.\n", lpi->linear_program->num_constraints().value(), lpi->linear_program->num_variables().value());
   std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(*lpi->parameters);
   lpi->niterations = 0;

   lpi->parameters->set_use_dual_simplex(true);
   return SolveInternal(lpi, false, time_limit);
}

/** calls barrier or interior point algorithm to solve the LP with crossover to simplex basis */
SCIP_RETCODE SCIPlpiSolveBarrier(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->parameters != NULL );

   SCIPerrorMessage("SCIPlpiSolveBarrier - Not supported.\n");

  return SCIP_NOTIMPLEMENTED;
}

/** start strong branching */
SCIP_RETCODE SCIPlpiStartStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );

   updateScaledLP(lpi);

   /* @todo Save state and do all the branching from there. */
   return SCIP_OKAY;
}

/** end strong branching */
SCIP_RETCODE SCIPlpiEndStrongbranch(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{  /*lint --e{715}*/
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );

   /* @todo Restore the saved state. */
   return SCIP_OKAY;
}

/** determine whether the dual bound is valid */
static
bool IsDualBoundValid(
   ProblemStatus         status              /**< status to be checked */
   )
{
   return status == ProblemStatus::OPTIMAL || status == ProblemStatus::DUAL_FEASIBLE || status == ProblemStatus::DUAL_UNBOUNDED;
}

/** performs strong branching iterations */
static
SCIP_RETCODE strongbranch(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col_index,          /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->scaled_lp != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPdebugMessage("calling strongbranching on variable %d (%d iterations)\n", col_index, itlim);

   /* We work on the scaled problem. */
   const ColIndex col(col_index);
   const Fractional lb = lpi->scaled_lp->variable_lower_bounds()[col];
   const Fractional ub = lpi->scaled_lp->variable_upper_bounds()[col];
   const double value = psol * lpi->scaler->VariableScalingFactor(col);

   /* Configure solver. */

   /* @todo Use the iteration limit once glop supports incrementality. */
   int num_iterations = 0;
   lpi->parameters->set_use_dual_simplex(true);
   lpi->solver->SetParameters(*(lpi->parameters));
   const Fractional eps = lpi->parameters->primal_feasibility_tolerance();

   std::unique_ptr<TimeLimit> time_limit = TimeLimit::FromParameters(*lpi->parameters);

   /* Down branch. */
   const Fractional newub = EPSCEIL(value - 1.0, eps);
   if ( newub >= lb - 0.5 )
   {
      lpi->scaled_lp->SetVariableBounds(col, lb, newub);

      if ( lpi->solver->Solve(*(lpi->scaled_lp), time_limit.get()).ok() )
      {
         num_iterations += (int) lpi->solver->GetNumberOfIterations();
         *down = lpi->solver->GetObjectiveValue();
         *downvalid = IsDualBoundValid(lpi->solver->GetProblemStatus()) ? TRUE : FALSE;

         SCIPdebugMessage("down: itlim=%d col=%d [%f,%f] obj=%f status=%d iter=%ld.\n", itlim, col_index, lb, EPSCEIL(value - 1.0, eps),
            lpi->solver->GetObjectiveValue(), (int) lpi->solver->GetProblemStatus(), lpi->solver->GetNumberOfIterations());
      }
      else
      {
         SCIPerrorMessage("error during solve");
         *down = 0.0;
         *downvalid = FALSE;
      }
   }
   else
   {
      if ( lpi->linear_program->IsMaximizationProblem() )
         *down = lpi->parameters->objective_lower_limit();
      else
         *down = lpi->parameters->objective_upper_limit();
      *downvalid = TRUE;
   }

   /* Up branch. */
   const Fractional newlb = EPSFLOOR(value + 1.0, eps);
   if ( newlb <= ub + 0.5 )
   {
      lpi->scaled_lp->SetVariableBounds(col, newlb, ub);

      if ( lpi->solver->Solve(*(lpi->scaled_lp), time_limit.get()).ok() )
      {
         num_iterations += (int) lpi->solver->GetNumberOfIterations();
         *up = lpi->solver->GetObjectiveValue();
         *upvalid = IsDualBoundValid(lpi->solver->GetProblemStatus()) ? TRUE : FALSE;

         SCIPdebugMessage("up: itlim=%d col=%d [%f,%f] obj=%f status=%d iter=%ld.\n", itlim, col_index, EPSFLOOR(value + 1.0, eps), ub,
            lpi->solver->GetObjectiveValue(), (int) lpi->solver->GetProblemStatus(), lpi->solver->GetNumberOfIterations());
      }
      else
      {
         SCIPerrorMessage("error during solve");
         *up = 0.0;
         *upvalid = FALSE;
      }
   }
   else
   {
      if (lpi->linear_program->IsMaximizationProblem())
         *up = lpi->parameters->objective_lower_limit();
      else
         *up = lpi->parameters->objective_upper_limit();
      *upvalid = TRUE;
   }

   /*  Restore bound. */
   lpi->scaled_lp->SetVariableBounds(col, lb, ub);
   if ( iter != NULL )
      *iter = num_iterations;

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPlpiStrongbranchFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col_index,          /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->scaled_lp != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPdebugMessage("calling strongbranching on fractional variable %d (%d iterations)\n", col_index, itlim);

   SCIP_CALL( strongbranch(lpi, col_index, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPlpiStrongbranchesFrac(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL) ;
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPerrorMessage("SCIPlpiStrongbranchesFrac - not implemented.\n");

   return SCIP_NOTIMPLEMENTED;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPlpiStrongbranchInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIP_CALL( strongbranch(lpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPlpiStrongbranchesInt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL) ;
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );

   SCIPerrorMessage("SCIPlpiStrongbranchesInt - not implemented.\n");

   return SCIP_NOTIMPLEMENTED;
}
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the LP */
SCIP_Bool SCIPlpiWasSolved(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );

   /* @todo Track this to avoid uneeded resolving. */
   return ( ! lpi->lp_modified_since_last_solve );
}

/** gets information about primal and dual feasibility of the current LP solution
 *
 *  The feasibility information is with respect to the last solving call and it is only relevant if SCIPlpiWasSolved()
 *  returns true. If the LP is changed, this information might be invalidated.
 *
 *  Note that @a primalfeasible and @a dualfeasible should only return true if the solver has proved the respective LP to
 *  be feasible. Thus, the return values should be equal to the values of SCIPlpiIsPrimalFeasible() and
 *  SCIPlpiIsDualFeasible(), respectively. Note that if feasibility cannot be proved, they should return false (even if
 *  the problem might actually be feasible).
 */
SCIP_RETCODE SCIPlpiGetSolFeasibility(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store dual feasibility status */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();
   *primalfeasible = (status == ProblemStatus::OPTIMAL || status == ProblemStatus::PRIMAL_FEASIBLE);
   *dualfeasible = (status == ProblemStatus::OPTIMAL || status == ProblemStatus::DUAL_FEASIBLE);

   SCIPdebugMessage("SCIPlpiGetSolFeasibility primal:%d dual:%d\n", *primalfeasible, *dualfeasible);

   return SCIP_OKAY;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiExistsPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return lpi->solver->GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

/** returns TRUE iff LP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPlpiHasPrimalRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return lpi->solver->GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

/** returns TRUE iff LP is proven to be primal unbounded */
SCIP_Bool SCIPlpiIsPrimalUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return lpi->solver->GetProblemStatus() == ProblemStatus::PRIMAL_UNBOUNDED;
}

/** returns TRUE iff LP is proven to be primal infeasible */
SCIP_Bool SCIPlpiIsPrimalInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();

   return status == ProblemStatus::DUAL_UNBOUNDED || status == ProblemStatus::PRIMAL_INFEASIBLE;
}

/** returns TRUE iff LP is proven to be primal feasible */
SCIP_Bool SCIPlpiIsPrimalFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();

   return status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::OPTIMAL;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiExistsDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();

   return status == ProblemStatus::DUAL_UNBOUNDED;
}

/** returns TRUE iff LP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPlpiHasDualRay(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();

   return status == ProblemStatus::DUAL_UNBOUNDED;
}

/** returns TRUE iff LP is proven to be dual unbounded */
SCIP_Bool SCIPlpiIsDualUnbounded(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();
   return status == ProblemStatus::DUAL_UNBOUNDED;
}

/** returns TRUE iff LP is proven to be dual infeasible */
SCIP_Bool SCIPlpiIsDualInfeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();
   return status == ProblemStatus::PRIMAL_UNBOUNDED || status == ProblemStatus::DUAL_INFEASIBLE;
}

/** returns TRUE iff LP is proven to be dual feasible */
SCIP_Bool SCIPlpiIsDualFeasible(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   const ProblemStatus status = lpi->solver->GetProblemStatus();

   return status == ProblemStatus::DUAL_FEASIBLE || status == ProblemStatus::OPTIMAL;
}

/** returns TRUE iff LP was solved to optimality */
SCIP_Bool SCIPlpiIsOptimal(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return lpi->solver->GetProblemStatus() == ProblemStatus::OPTIMAL;
}

/** returns TRUE iff current LP solution is stable
 *
 *  This function should return true if the solution is reliable, i.e., feasible and optimal (or proven
 *  infeasible/unbounded) with respect to the original problem. The optimality status might be with respect to a scaled
 *  version of the problem, but the solution might not be feasible to the unscaled original problem; in this case,
 *  SCIPlpiIsStable() should return false.
 */
SCIP_Bool SCIPlpiIsStable(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   /* For correctness, we need to report "unstable" if Glop was not able to prove optimality because of numerical
    * issues. Currently Glop still reports primal/dual feasible if at the end, one status is within the tolerance but not
    * the other. */
   const ProblemStatus status = lpi->solver->GetProblemStatus();
   if ( (status == ProblemStatus::PRIMAL_FEASIBLE || status == ProblemStatus::DUAL_FEASIBLE) &&
      ! SCIPlpiIsObjlimExc(lpi) && ! SCIPlpiIsIterlimExc(lpi) && ! SCIPlpiIsTimelimExc(lpi))
   {
      SCIPdebugMessage("OPTIMAL not reached and no limit: unstable.\n");
      return FALSE;
   }

   if ( status == ProblemStatus::ABNORMAL || status == ProblemStatus::INVALID_PROBLEM || status == ProblemStatus::IMPRECISE )
      return FALSE;

   /* If we have a regular basis and the condition limit is set, we compute (an upper bound on) the condition number of
    * the basis; everything above the specified threshold is then counted as instable. */
   if ( lpi->checkcondition && (SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi)) )
   {
      SCIP_RETCODE retcode;
      SCIP_Real kappa;

      retcode = SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_ESTIMCONDITION, &kappa);
      if ( retcode != SCIP_OKAY )
         SCIPABORT();
      assert( kappa != SCIP_INVALID ); /*lint !e777*/

      if ( kappa > lpi->conditionlimit )
         return FALSE;
   }

   return TRUE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPlpiIsObjlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return lpi->solver->objective_limit_reached();
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPlpiIsIterlimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->niterations >= (int) lpi->solver->GetNumberOfIterations() );

   int maxiter = (int) lpi->parameters->max_number_of_iterations();
   return maxiter >= 0 && lpi->niterations >= maxiter;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPlpiIsTimelimExc(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return lpi->lp_time_limit_was_reached;
}

/** returns the internal solution status of the solver */
int SCIPlpiGetInternalStatus(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return static_cast<int>(lpi->solver->GetProblemStatus());
}

/** tries to reset the internal status of the LP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPlpiIgnoreInstability(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( success != NULL );

   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPlpiGetObjval(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( objval != NULL );

   *objval = lpi->solver->GetObjectiveValue();

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors for feasible LPs
 *
 *  Before calling this function, the caller must ensure that the LP has been solved to optimality, i.e., that
 *  SCIPlpiIsOptimal() returns true.
 */
SCIP_RETCODE SCIPlpiGetSol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   SCIPdebugMessage("SCIPlpiGetSol\n");
   if ( objval != NULL )
      *objval = lpi->solver->GetObjectiveValue();

   const ColIndex num_cols = lpi->linear_program->num_variables();
   for (ColIndex col(0); col < num_cols; ++col)
   {
      int i = col.value();

      if ( primsol != NULL )
         primsol[i] = lpi->scaler->UnscaleVariableValue(col, lpi->solver->GetVariableValue(col));

      if ( redcost != NULL )
         redcost[i] = lpi->scaler->UnscaleReducedCost(col, lpi->solver->GetReducedCost(col));
   }

   const RowIndex num_rows = lpi->linear_program->num_constraints();
   for (RowIndex row(0); row < num_rows; ++row)
   {
      int j = row.value();

      if ( dualsol != NULL )
         dualsol[j] = lpi->scaler->UnscaleDualValue(row, lpi->solver->GetDualValue(row));

      if ( activity != NULL )
         activity[j] = lpi->scaler->UnscaleConstraintActivity(row, lpi->solver->GetConstraintActivity(row));
   }

   return SCIP_OKAY;
}

/** gets primal ray for unbounded LPs */
SCIP_RETCODE SCIPlpiGetPrimalRay(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( ray != NULL );

   SCIPdebugMessage("SCIPlpiGetPrimalRay\n");

   const ColIndex num_cols = lpi->linear_program->num_variables();
   const DenseRow& primal_ray = lpi->solver->GetPrimalRay();
   for (ColIndex col(0); col < num_cols; ++col)
      ray[col.value()] = lpi->scaler->UnscaleVariableValue(col, primal_ray[col]);

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPlpiGetDualfarkas(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( dualfarkas != NULL );

   SCIPdebugMessage("SCIPlpiGetDualfarkas\n");

   const RowIndex num_rows = lpi->linear_program->num_constraints();
   const DenseColumn& dual_ray = lpi->solver->GetDualRay();
   for (RowIndex row(0); row < num_rows; ++row)
      dualfarkas[row.value()] = -lpi->scaler->UnscaleDualValue(row, dual_ray[row]);  /* reverse sign */

   return SCIP_OKAY;
}

/** gets the number of LP iterations of the last solve call */
SCIP_RETCODE SCIPlpiGetIterations(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( iterations != NULL );

   *iterations = (int) lpi->niterations;

   return SCIP_OKAY;
}

/** gets information about the quality of an LP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The LPI should return SCIP_INVALID for @p quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPlpiGetRealSolQuality(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPSOLQUALITY     qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( quality != NULL );

   SCIPdebugMessage("Requesting solution quality: quality %d\n", qualityindicator);

   switch ( qualityindicator )
   {
   case SCIP_LPSOLQUALITY_ESTIMCONDITION:
      *quality = lpi->solver->GetBasisFactorization().ComputeInfinityNormConditionNumber();
      break;

   case SCIP_LPSOLQUALITY_EXACTCONDITION:
      *quality = lpi->solver->GetBasisFactorization().ComputeInfinityNormConditionNumberUpperBound();
      break;

   default:
      SCIPerrorMessage("Solution quality %d unknown.\n", qualityindicator);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * LP Basis Methods
 */

/**@name LP Basis Methods */
/**@{ */

/** convert Glop variable basis status to SCIP status */
static
SCIP_BASESTAT ConvertGlopVariableStatus(
   VariableStatus        status,             /**< variable status */
   Fractional            rc                  /**< reduced cost of variable */
   )
{
   switch ( status )
   {
   case VariableStatus::BASIC:
      return SCIP_BASESTAT_BASIC;
   case VariableStatus::AT_UPPER_BOUND:
      return SCIP_BASESTAT_UPPER;
   case VariableStatus::AT_LOWER_BOUND:
      return SCIP_BASESTAT_LOWER;
   case VariableStatus::FREE:
      return SCIP_BASESTAT_ZERO;
   case VariableStatus::FIXED_VALUE:
      return rc > 0.0 ? SCIP_BASESTAT_LOWER : SCIP_BASESTAT_UPPER;
   default:
      SCIPerrorMessage("invalid Glop basis status.\n");
      abort();
   }
}

/** convert Glop constraint basis status to SCIP status */
static
SCIP_BASESTAT ConvertGlopConstraintStatus(
   ConstraintStatus      status,             /**< constraint status */
   Fractional            dual                /**< dual variable value */
   )
{
   switch ( status )
   {
   case ConstraintStatus::BASIC:
      return SCIP_BASESTAT_BASIC;
   case ConstraintStatus::AT_UPPER_BOUND:
      return SCIP_BASESTAT_UPPER;
   case ConstraintStatus::AT_LOWER_BOUND:
      return SCIP_BASESTAT_LOWER;
   case ConstraintStatus::FREE:
      return SCIP_BASESTAT_ZERO;
   case ConstraintStatus::FIXED_VALUE:
      return dual > 0.0 ? SCIP_BASESTAT_LOWER : SCIP_BASESTAT_UPPER;
   default:
      SCIPerrorMessage("invalid Glop basis status.\n");
      abort();
   }
}

/** Convert SCIP variable status to Glop status */
static
VariableStatus ConvertSCIPVariableStatus(
   int                   status              /**< SCIP variable status */
   )
{
   switch ( status )
   {
   case SCIP_BASESTAT_BASIC:
      return VariableStatus::BASIC;
   case SCIP_BASESTAT_UPPER:
      return VariableStatus::AT_UPPER_BOUND;
   case SCIP_BASESTAT_LOWER:
      return VariableStatus::AT_LOWER_BOUND;
   case SCIP_BASESTAT_ZERO:
      return VariableStatus::FREE;
   default:
      SCIPerrorMessage("invalid SCIP basis status.\n");
      abort();
   }
}

/** Convert a SCIP constraint status to its corresponding Glop slack VariableStatus.
 *
 *  Note that we swap the upper/lower bounds.
 */
static
VariableStatus ConvertSCIPConstraintStatusToSlackStatus(
   int                   status              /**< SCIP constraint status */
   )
{
   switch ( status )
   {
   case SCIP_BASESTAT_BASIC:
      return VariableStatus::BASIC;
   case SCIP_BASESTAT_UPPER:
      return VariableStatus::AT_LOWER_BOUND;
   case SCIP_BASESTAT_LOWER:
      return VariableStatus::AT_UPPER_BOUND;
   case SCIP_BASESTAT_ZERO:
      return VariableStatus::FREE;
   default:
      SCIPerrorMessage("invalid SCIP basis status.\n");
      abort();
   }
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPlpiGetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   SCIPdebugMessage("SCIPlpiGetBase\n");

   assert( lpi->solver->GetProblemStatus() ==  ProblemStatus::OPTIMAL );

   if ( cstat != NULL )
   {
      const ColIndex num_cols = lpi->linear_program->num_variables();
      for (ColIndex col(0); col < num_cols; ++col)
      {
         int i = col.value();
         cstat[i] = (int) ConvertGlopVariableStatus(lpi->solver->GetVariableStatus(col), lpi->solver->GetReducedCost(col));
      }
   }

   if ( rstat != NULL )
   {
      const RowIndex num_rows = lpi->linear_program->num_constraints();
      for (RowIndex row(0); row < num_rows; ++row)
      {
         int i = row.value();
         rstat[i] = (int) ConvertGlopConstraintStatus(lpi->solver->GetConstraintStatus(row), lpi->solver->GetDualValue(row));
      }
   }

   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPlpiSetBase(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const int*            cstat,              /**< array with column basis status */
   const int*            rstat               /**< array with row basis status */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );

   const ColIndex num_cols = lpi->linear_program->num_variables();
   const RowIndex num_rows = lpi->linear_program->num_constraints();

   assert( cstat != NULL || num_cols == 0 );
   assert( rstat != NULL || num_rows == 0 );

   SCIPdebugMessage("SCIPlpiSetBase\n");

   BasisState state;
   state.statuses.reserve(ColIndex(num_cols.value() + num_rows.value()));

   for (ColIndex col(0); col < num_cols; ++col)
      state.statuses[col] = ConvertSCIPVariableStatus(cstat[col.value()]);

   for (RowIndex row(0); row < num_rows; ++row)
      state.statuses[num_cols + RowToColIndex(row)] = ConvertSCIPConstraintStatusToSlackStatus(cstat[row.value()]);

   lpi->solver->LoadStateForNextSolve(state);

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
SCIP_RETCODE SCIPlpiGetBasisInd(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );
   assert( bind != NULL );

   SCIPdebugMessage("SCIPlpiGetBasisInd\n");

   /* the order is important! */
   const ColIndex num_cols = lpi->linear_program->num_variables();
   const RowIndex num_rows = lpi->linear_program->num_constraints();
   for (RowIndex row(0); row < num_rows; ++row)
   {
      const ColIndex col = lpi->solver->GetBasis(row);
      if (col < num_cols)
         bind[row.value()] = col.value();
      else
      {
         assert( col < num_cols.value() + num_rows.value() );
         bind[row.value()] = -1 - (col - num_cols).value();
      }
   }

   return SCIP_OKAY;
}

/** get row of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvRow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->scaler != NULL );
   assert( coef != NULL );

   lpi->solver->GetBasisFactorization().LeftSolveForUnitRow(ColIndex(r), lpi->tmp_row);
   lpi->scaler->UnscaleUnitRowLeftSolve(lpi->solver->GetBasis(RowIndex(r)), lpi->tmp_row);

   const ColIndex size = lpi->tmp_row->values.size();
   assert( size.value() == lpi->linear_program->num_constraints() );

   /* if we want a sparse vector */
   if ( ninds != NULL && inds != NULL )
   {
      *ninds = 0;
      /* Vectors in Glop might be stored in dense or sparse format depending on the values. If non_zeros are given, we
       * can directly loop over the non_zeros, otherwise we have to collect the nonzeros. */
      if ( ! lpi->tmp_row->non_zeros.empty() )
      {
         ScatteredRowIterator end = lpi->tmp_row->end();
         for (ScatteredRowIterator iter = lpi->tmp_row->begin(); iter != end; ++iter)
         {
            int idx = (*iter).column().value();
            assert( 0 <= idx && idx < lpi->linear_program->num_constraints() );
            coef[idx] = (*iter).coefficient();
            inds[(*ninds)++] = idx;
         }
      }
      else
      {
         /* use dense access to tmp_row */
         const Fractional eps = lpi->parameters->primal_feasibility_tolerance();
         for (ColIndex col(0); col < size; ++col)
         {
            SCIP_Real val = (*lpi->tmp_row)[col];
            if ( fabs(val) >= eps )
            {
               coef[col.value()] = val;
               inds[(*ninds)++] = col.value();
            }
         }
      }
      return SCIP_OKAY;
   }

   /* dense version */
   for (ColIndex col(0); col < size; ++col)
      coef[col.value()] = (*lpi->tmp_row)[col];

   if ( ninds != NULL )
      *ninds = -1;

   return SCIP_OKAY;
}

/** get column of inverse basis matrix B^-1
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvCol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP;
                                              *   you have to call SCIPlpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the LP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef,               /**< pointer to store the coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->scaler != NULL );
   assert( coef != NULL );

   /* we need to loop through the rows to extract the values for column c */
   const ColIndex col(c);
   const RowIndex num_rows = lpi->linear_program->num_constraints();

   /* if we want a sparse vector */
   if ( ninds != NULL && inds != NULL )
   {
      const Fractional eps = lpi->parameters->primal_feasibility_tolerance();

      *ninds = 0;
      for (int row = 0; row < num_rows; ++row)
      {
         lpi->solver->GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row), lpi->tmp_row);
         lpi->scaler->UnscaleUnitRowLeftSolve(lpi->solver->GetBasis(RowIndex(row)), lpi->tmp_row);

         SCIP_Real val = (*lpi->tmp_row)[col];
         if ( fabs(val) >= eps )
         {
            coef[row] = val;
            inds[(*ninds)++] = row;
         }
      }
      return SCIP_OKAY;
   }

   /* dense version */
   for (int row = 0; row < num_rows; ++row)
   {
      lpi->solver->GetBasisFactorization().LeftSolveForUnitRow(ColIndex(row), lpi->tmp_row);
      lpi->scaler->UnscaleUnitRowLeftSolve(lpi->solver->GetBasis(RowIndex(row)), lpi->tmp_row);
      coef[row] = (*lpi->tmp_row)[col];
   }

   if ( ninds != NULL )
      *ninds = -1;

   return SCIP_OKAY;
}

/** get row of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvARow(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPlpiGetBInvRow(), or NULL */
   SCIP_Real*            coef,               /**< vector to return coefficients of the row */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->scaler != NULL );
   assert( coef != NULL );

   /* get row of basis inverse, loop through columns and muliply with matrix */
   lpi->solver->GetBasisFactorization().LeftSolveForUnitRow(ColIndex(r), lpi->tmp_row);
   lpi->scaler->UnscaleUnitRowLeftSolve(lpi->solver->GetBasis(RowIndex(r)), lpi->tmp_row);

   const ColIndex num_cols = lpi->linear_program->num_variables();

   /* if we want a sparse vector */
   if ( ninds != NULL && inds != NULL )
   {
      const Fractional eps = lpi->parameters->primal_feasibility_tolerance();

      *ninds = 0;
      for (ColIndex col(0); col < num_cols; ++col)
      {
         SCIP_Real val = operations_research::glop::ScalarProduct(lpi->tmp_row->values, lpi->linear_program->GetSparseColumn(col));
         if ( fabs(val) >= eps )
         {
            coef[col.value()] = val;
            inds[(*ninds)++] = col.value();
         }
      }
      return SCIP_OKAY;
   }

   /* dense version */
   for (ColIndex col(0); col < num_cols; ++col)
      coef[col.value()] = operations_research::glop::ScalarProduct(lpi->tmp_row->values, lpi->linear_program->GetSparseColumn(col));

   if ( ninds != NULL )
      *ninds = -1;

   return SCIP_OKAY;
}

/** get column of inverse basis matrix times constraint matrix B^-1 * A
 *
 *  @note The LP interface defines slack variables to have coefficient +1. This means that if, internally, the LP solver
 *        uses a -1 coefficient, then rows associated with slacks variables whose coefficient is -1, should be negated;
 *        see also the explanation in lpi.h.
 */
SCIP_RETCODE SCIPlpiGetBInvACol(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef,               /**< vector to return coefficients of the column */
   int*                  inds,               /**< array to store the non-zero indices, or NULL */
   int*                  ninds               /**< pointer to store the number of non-zero indices, or NULL
                                              *   (-1: if we do not store sparsity information) */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( lpi->solver != NULL );
   assert( lpi->scaler != NULL );
   assert( coef != NULL );

   lpi->solver->GetBasisFactorization().RightSolveForProblemColumn(ColIndex(c), lpi->tmp_column);
   lpi->scaler->UnscaleColumnRightSolve(lpi->solver->GetBasisVector(), ColIndex(c), lpi->tmp_column);

   const RowIndex num_rows = lpi->tmp_column->values.size();

   /* if we want a sparse vector */
   if ( ninds != NULL && inds != NULL )
   {
      *ninds = 0;
      /* Vectors in Glop might be stored in dense or sparse format depending on the values. If non_zeros are given, we
       * can directly loop over the non_zeros, otherwise we have to collect the nonzeros. */
      if ( ! lpi->tmp_column->non_zeros.empty() )
      {
         ScatteredColumnIterator end = lpi->tmp_column->end();
         for (ScatteredColumnIterator iter = lpi->tmp_column->begin(); iter != end; ++iter)
         {
            int idx = (*iter).row().value();
            assert( 0 <= idx && idx < num_rows );
            coef[idx] = (*iter).coefficient();
            inds[(*ninds)++] = idx;
         }
      }
      else
      {
         /* use dense access to tmp_column */
         const Fractional eps = lpi->parameters->primal_feasibility_tolerance();
         for (RowIndex row(0); row < num_rows; ++row)
         {
            SCIP_Real val = (*lpi->tmp_column)[row];
            if ( fabs(val) > eps )
            {
               coef[row.value()] = val;
               inds[(*ninds)++] = row.value();
            }
         }
      }
      return SCIP_OKAY;
   }

   /* dense version */
   for (RowIndex row(0); row < num_rows; ++row)
      coef[row.value()] = (*lpi->tmp_column)[row];

   if ( ninds != NULL )
      *ninds = -1;

   return SCIP_OKAY;
}

/**@} */




/*
 * LP State Methods
 */

/**@name LP State Methods */
/**@{ */

/* SCIP_LPiState stores basis information and is implemented by the glop BasisState class.
 * However, because in type_lpi.h there is
 *   typedef struct SCIP_LPiState SCIP_LPISTATE;
 * We cannot just use a typedef here and SCIP_LPiState needs to be a struct.
 */
struct SCIP_LPiState : BasisState {};

/** stores LPi state (like basis information) into lpistate object */
SCIP_RETCODE SCIPlpiGetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpistate != NULL );

   *lpistate = static_cast<SCIP_LPISTATE*>(new BasisState(lpi->solver->GetState()));

   return SCIP_OKAY;
}

/** loads LPi state (like basis information) into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetState()
 */
SCIP_RETCODE SCIPlpiSetState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPISTATE*  lpistate            /**< LPi state information (like basis information), or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpistate != NULL );

   lpi->solver->LoadStateForNextSolve(*lpistate);

   return SCIP_OKAY;
}

/** clears current LPi state (like basis information) of the solver */
SCIP_RETCODE SCIPlpiClearState(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   lpi->solver->ClearStateForNextSolve();

   return SCIP_OKAY;
}

/** frees LPi state information */
SCIP_RETCODE SCIPlpiFreeState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPISTATE**       lpistate            /**< pointer to LPi state information (like basis information) */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );
   assert( lpistate != NULL );

   delete *lpistate;
   *lpistate = NULL;

   return SCIP_OKAY;
}

/** checks, whether the given LP state contains simplex basis information */
SCIP_Bool SCIPlpiHasStateBasis(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPISTATE*        lpistate            /**< LP state information (like basis information), or NULL */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   return lpistate != NULL;
}

/** reads LP state (like basis information from a file */
SCIP_RETCODE SCIPlpiReadState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   SCIPerrorMessage("SCIPlpiReadState - not implemented.\n");

   return SCIP_NOTIMPLEMENTED;
}

/** writes LPi state (i.e. basis information) to a file */
SCIP_RETCODE SCIPlpiWriteState(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert( lpi != NULL );
   assert( lpi->solver != NULL );

   SCIPerrorMessage("SCIPlpiWriteState - not implemented.\n");

   return SCIP_NOTIMPLEMENTED;
}

/**@} */




/*
 * LP Pricing Norms Methods
 */

/**@name LP Pricing Norms Methods */
/**@{ */

/* SCIP_LPiNorms stores norm information so they are not recomputed from one state to the next. */
/* @todo Implement this. */
struct SCIP_LPiNorms {};

/** stores LPi pricing norms information
 *
 *  @todo store primal norms as well?
 */
SCIP_RETCODE SCIPlpiGetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information */
   )
{
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpi->solver != NULL );
   assert( lpinorms != NULL );

   return SCIP_OKAY;
}

/** loads LPi pricing norms into solver; note that the LP might have been extended with additional
 *  columns and rows since the state was stored with SCIPlpiGetNorms()
 */
SCIP_RETCODE SCIPlpiSetNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const SCIP_LPINORMS*  lpinorms            /**< LPi pricing norms information, or NULL */
   )
{
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpi->solver != NULL );
   assert( lpinorms != NULL );

   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPlpiFreeNorms(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_LPINORMS**       lpinorms            /**< pointer to LPi pricing norms information, or NULL */
   )
{
   assert( lpi != NULL );
   assert( blkmem != NULL );
   assert( lpi->solver != NULL );
   assert( lpinorms != NULL );

   return SCIP_OKAY;
}

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of LP */
SCIP_RETCODE SCIPlpiGetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   assert( lpi != NULL );
   assert( lpi->parameters != NULL );

   /* Not (yet) supported by Glop: SCIP_LPPAR_FASTMIP, SCIP_LPPAR_POLISHING, SCIP_LPPAR_REFACTOR */
   switch ( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      *ival = (int) lpi->from_scratch;
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_FROMSCRATCH = %d.\n", *ival);
      break;
   case SCIP_LPPAR_LPINFO:
      *ival = (int) lpi->lp_info;
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_LPINFO = %d.\n", *ival);
      break;
   case SCIP_LPPAR_LPITLIM:
      *ival = (int) lpi->parameters->max_number_of_iterations();
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_LPITLIM = %d.\n", *ival);
      break;
   case SCIP_LPPAR_PRESOLVING:
      *ival = lpi->parameters->use_preprocessing();
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_PRESOLVING = %d.\n", *ival);
      break;
   case SCIP_LPPAR_PRICING:
      *ival = lpi->pricing;
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_PRICING = %d.\n", *ival);
      break;
#ifndef NOSCALING
   case SCIP_LPPAR_SCALING:
      *ival = lpi->parameters->use_scaling();
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_SCALING = %d.\n", *ival);
      break;
#endif
   case SCIP_LPPAR_THREADS:
      *ival = lpi->numthreads;
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_THREADS = %d.\n", *ival);
      break;
   case SCIP_LPPAR_TIMING:
      *ival = lpi->timing;
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_TIMING = %d.\n", *ival);
      break;
   case SCIP_LPPAR_RANDOMSEED:
      *ival = (int) lpi->parameters->random_seed();
      SCIPdebugMessage("SCIPlpiGetIntpar: SCIP_LPPAR_RANDOMSEED = %d.\n", *ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets integer parameter of LP */
SCIP_RETCODE SCIPlpiSetIntpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert( lpi != NULL );
   assert( lpi->parameters != NULL );

   switch ( type )
   {
   case SCIP_LPPAR_FROMSCRATCH:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_FROMSCRATCH -> %d.\n", ival);
      lpi->from_scratch = (bool) ival;
      break;
   case SCIP_LPPAR_LPINFO:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_LPINFO -> %d.\n", ival);
      if ( ival == 0 )
      {
         (void) google::SetVLOGLevel("*", google::GLOG_INFO);
         lpi->lp_info = false;
      }
      else
      {
         (void) google::SetVLOGLevel("*", google::GLOG_ERROR);
         lpi->lp_info = true;
      }
      break;
   case SCIP_LPPAR_LPITLIM:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_LPITLIM -> %d.\n", ival);
      lpi->parameters->set_max_number_of_iterations(ival);
      break;
   case SCIP_LPPAR_PRESOLVING:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_PRESOLVING -> %d.\n", ival);
      lpi->parameters->set_use_preprocessing(ival);
      break;
   case SCIP_LPPAR_PRICING:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_PRICING -> %d.\n", ival);
      lpi->pricing = (SCIP_Pricing)ival;
      switch ( lpi->pricing )
      {
      case SCIP_PRICING_LPIDEFAULT:
      case SCIP_PRICING_AUTO:
      case SCIP_PRICING_PARTIAL:
      case SCIP_PRICING_STEEP:
      case SCIP_PRICING_STEEPQSTART:
         lpi->parameters->set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_STEEPEST_EDGE);
         break;
      case SCIP_PRICING_FULL:
         /* Dantzig does not really fit, but use it anyway */
         lpi->parameters->set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_DANTZIG);
         break;
      case SCIP_PRICING_DEVEX:
         lpi->parameters->set_feasibility_rule(operations_research::glop::GlopParameters_PricingRule_DEVEX);
         break;
      default:
         return SCIP_PARAMETERUNKNOWN;
      }
      break;
#ifndef NOSCALING
   case SCIP_LPPAR_SCALING:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_SCALING -> %d.\n", ival);
      lpi->parameters->set_use_scaling(ival);
      break;
#endif
   case SCIP_LPPAR_THREADS:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_THREADS -> %d.\n", ival);
      assert( ival >= 0 );
      lpi->numthreads = ival;
      if ( ival == 0 )
         lpi->parameters->set_num_omp_threads(1);
      else
         lpi->parameters->set_num_omp_threads(ival);
      break;
   case SCIP_LPPAR_TIMING:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_TIMING -> %d.\n", ival);
      assert( 0 <= ival && ival <= 2 );
      lpi->timing = ival;
      if ( ival == 1 )
         absl::SetFlag(&FLAGS_time_limit_use_usertime, true);
      else
         absl::SetFlag(&FLAGS_time_limit_use_usertime, false);
      break;
   case SCIP_LPPAR_RANDOMSEED:
      SCIPdebugMessage("SCIPlpiSetIntpar: SCIP_LPPAR_RANDOMSEED -> %d.\n", ival);
      assert( ival >= 0 );
      lpi->parameters->set_random_seed(ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** gets floating point parameter of LP */
SCIP_RETCODE SCIPlpiGetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   assert( lpi != NULL );
   assert( lpi->parameters != NULL );

   /* Not (yet) supported by Glop: SCIP_LPPAR_ROWREPSWITCH, SCIP_LPPAR_BARRIERCONVTOL */
   switch ( type )
   {
   case SCIP_LPPAR_FEASTOL:
      *dval = lpi->parameters->primal_feasibility_tolerance();
      SCIPdebugMessage("SCIPlpiGetRealpar: SCIP_LPPAR_FEASTOL = %g.\n", *dval);
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      *dval = lpi->parameters->dual_feasibility_tolerance();
      SCIPdebugMessage("SCIPlpiGetRealpar: SCIP_LPPAR_DUALFEASTOL = %g.\n", *dval);
      break;
   case SCIP_LPPAR_OBJLIM:
      if (lpi->linear_program->IsMaximizationProblem())
         *dval = lpi->parameters->objective_lower_limit();
      else
         *dval = lpi->parameters->objective_upper_limit();
      SCIPdebugMessage("SCIPlpiGetRealpar: SCIP_LPPAR_OBJLIM = %f.\n", *dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      if ( absl::GetFlag(FLAGS_time_limit_use_usertime) )
         *dval = lpi->parameters->max_time_in_seconds();
      else
         *dval = lpi->parameters->max_deterministic_time();
      SCIPdebugMessage("SCIPlpiGetRealpar: SCIP_LPPAR_LPTILIM = %f.\n", *dval);
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      *dval = lpi->conditionlimit;
      break;
#ifdef SCIP_DISABLED_CODE
   /* currently do not apply Markowitz parameter, since the default value does not seem suitable for Glop */
   case SCIP_LPPAR_MARKOWITZ:
      *dval = lpi->parameters->markowitz_singularity_threshold();
      SCIPdebugMessage("SCIPlpiGetRealpar: SCIP_LPPAR_MARKOWITZ = %f.\n", *dval);
      break;
#endif
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of LP */
SCIP_RETCODE SCIPlpiSetRealpar(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_LPPARAM          type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert( lpi != NULL );
   assert( lpi->parameters != NULL );

   switch( type )
   {
   case SCIP_LPPAR_FEASTOL:
      SCIPdebugMessage("SCIPlpiSetRealpar: SCIP_LPPAR_FEASTOL -> %g.\n", dval);
      lpi->parameters->set_primal_feasibility_tolerance(dval);
      break;
   case SCIP_LPPAR_DUALFEASTOL:
      SCIPdebugMessage("SCIPlpiSetRealpar: SCIP_LPPAR_DUALFEASTOL -> %g.\n", dval);
      lpi->parameters->set_dual_feasibility_tolerance(dval);
      break;
   case SCIP_LPPAR_OBJLIM:
      SCIPdebugMessage("SCIPlpiSetRealpar: SCIP_LPPAR_OBJLIM -> %f.\n", dval);
      if (lpi->linear_program->IsMaximizationProblem())
         lpi->parameters->set_objective_lower_limit(dval);
      else
         lpi->parameters->set_objective_upper_limit(dval);
      break;
   case SCIP_LPPAR_LPTILIM:
      SCIPdebugMessage("SCIPlpiSetRealpar: SCIP_LPPAR_LPTILIM -> %f.\n", dval);
      if ( absl::GetFlag(FLAGS_time_limit_use_usertime) )
         lpi->parameters->set_max_time_in_seconds(dval);
      else
         lpi->parameters->set_max_deterministic_time(dval);
      break;
   case SCIP_LPPAR_CONDITIONLIMIT:
      lpi->conditionlimit = dval;
      lpi->checkcondition = (dval >= 0.0);
      break;
#ifdef SCIP_DISABLED_CODE
   /* currently do not apply Markowitz parameter, since the default value does not seem suitable for Glop */
   case SCIP_LPPAR_MARKOWITZ:
      SCIPdebugMessage("SCIPlpiSetRealpar: SCIP_LPPAR_MARKOWITZ -> %f.\n", dval);
      lpi->parameters->set_markowitz_singularity_threshold(dval);
      break;
#endif
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** interrupts the currently ongoing lp solve or disables the interrupt */
SCIP_RETCODE SCIPlpiInterrupt(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Bool             interrupt           /**< TRUE if interrupt should be set, FALSE if it should be disabled */
   )
{
   /*lint --e{715}*/
   assert(lpi != NULL);

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the LP solver */
SCIP_Real SCIPlpiInfinity(
   SCIP_LPI*             lpi                 /**< LP interface structure */
   )
{
   assert( lpi != NULL );
   return std::numeric_limits<SCIP_Real>::infinity();
}

/** checks if given value is treated as infinity in the LP solver */
SCIP_Bool SCIPlpiIsInfinity(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{
   assert( lpi != NULL );

   return val == std::numeric_limits<SCIP_Real>::infinity();
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads LP from a file */
SCIP_RETCODE SCIPlpiReadLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( fname != NULL );

   const std::string filespec(fname);
   MPModelProto proto;
   if ( ! ReadFileToProto(filespec, &proto) )
   {
      SCIPerrorMessage("Could not read <%s>\n", fname);
      return SCIP_READERROR;
   }
   lpi->linear_program->Clear();
   MPModelProtoToLinearProgram(proto, lpi->linear_program);

   return SCIP_OKAY;
}

/** writes LP to a file */
SCIP_RETCODE SCIPlpiWriteLP(
   SCIP_LPI*             lpi,                /**< LP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert( lpi != NULL );
   assert( lpi->linear_program != NULL );
   assert( fname != NULL );

   MPModelProto proto;
   LinearProgramToMPModelProto(*lpi->linear_program, &proto);
   const std::string filespec(fname);
   if ( ! WriteProtoToFile(filespec, proto, operations_research::ProtoWriteFormat::kProtoText, true) )
   {
      SCIPerrorMessage("Could not write <%s>\n", fname);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/**@} */
