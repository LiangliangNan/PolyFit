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

/**@file   pub_matrix.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for matrix
 * @author Dieter Weninger
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MATRIX_H__
#define __SCIP_PUB_MATRIX_H__

#include "scip/def.h"
#include "scip/type_var.h"
#include "scip/type_cons.h"
#include "scip/type_matrix.h"

#ifdef NDEBUG
#include "scip/struct_matrix.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * methods for matrix access
 */

/** get column based start pointer of values */
SCIP_EXPORT
SCIP_Real* SCIPmatrixGetColValPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get column based start pointer of row indices */
SCIP_EXPORT
int* SCIPmatrixGetColIdxPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get the number of non-zero entries of this column */
SCIP_EXPORT
int SCIPmatrixGetColNNonzs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get number of columns of the matrix */
SCIP_EXPORT
int SCIPmatrixGetNColumns(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   );

/** get upper bound of column */
SCIP_EXPORT
SCIP_Real SCIPmatrixGetColUb(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get lower bound of column */
SCIP_EXPORT
SCIP_Real SCIPmatrixGetColLb(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get number of uplocks of column */
SCIP_EXPORT
int SCIPmatrixGetColNUplocks(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get number of downlocks of column */
SCIP_EXPORT
int SCIPmatrixGetColNDownlocks(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get variable pointer of column */
SCIP_EXPORT
SCIP_VAR* SCIPmatrixGetVar(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get name of column/variable */
SCIP_EXPORT
const char* SCIPmatrixGetColName(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get row based start pointer of values */
SCIP_EXPORT
SCIP_Real* SCIPmatrixGetRowValPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get row based start pointer of column indices */
SCIP_EXPORT
int* SCIPmatrixGetRowIdxPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of non-zeros of this row */
SCIP_EXPORT
int SCIPmatrixGetRowNNonzs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get name of row */
SCIP_EXPORT
const char* SCIPmatrixGetRowName(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of rows of the matrix */
SCIP_EXPORT
int SCIPmatrixGetNRows(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   );

/** get left-hand-side of row */
SCIP_EXPORT
SCIP_Real SCIPmatrixGetRowLhs(
   SCIP_MATRIX*          matrix,             /**< matrix instace */
   int                   row                 /**< row index */
   );

/** get right-hand-side of row */
SCIP_EXPORT
SCIP_Real SCIPmatrixGetRowRhs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** flag indicating if right-hand-side of row is infinity */
SCIP_EXPORT
SCIP_Bool SCIPmatrixIsRowRhsInfinity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of non-zeros of matrix */
SCIP_EXPORT
int SCIPmatrixGetNNonzs(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   );

/** get minimal activity of row */
SCIP_EXPORT
SCIP_Real SCIPmatrixGetRowMinActivity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get maximal activity of row */
SCIP_EXPORT
SCIP_Real SCIPmatrixGetRowMaxActivity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of negative infinities present within minimal activity */
SCIP_EXPORT
int SCIPmatrixGetRowNMinActNegInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of positive infinities present within minimal activity */
SCIP_EXPORT
int SCIPmatrixGetRowNMinActPosInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of negative infinities present within maximal activity */
SCIP_EXPORT
int SCIPmatrixGetRowNMaxActNegInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of positive infinities present within maximal activity */
SCIP_EXPORT
int SCIPmatrixGetRowNMaxActPosInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get constraint pointer for constraint representing row */
SCIP_EXPORT
SCIP_CONS* SCIPmatrixGetCons(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get if conflicting uplocks of variable present */
SCIP_EXPORT
SCIP_Bool SCIPmatrixUplockConflict(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get if conflicting downlocks of variable present */
SCIP_EXPORT
SCIP_Bool SCIPmatrixDownlockConflict(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );


#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPmatrixGetColValPtr(matrix,col)             (matrix->colmatval + matrix->colmatbeg[col])
#define SCIPmatrixGetColIdxPtr(matrix,col)             (matrix->colmatind + matrix->colmatbeg[col])
#define SCIPmatrixGetColNNonzs(matrix,col)             (matrix->colmatcnt[col])
#define SCIPmatrixGetNColumns(matrix)                  (matrix->ncols)
#define SCIPmatrixGetColUb(matrix,col)                 (matrix->ub[col])
#define SCIPmatrixGetColLb(matrix,col)                 (matrix->lb[col])
#define SCIPmatrixGetColNUplocks(matrix,col)           (matrix->nuplocks[col])
#define SCIPmatrixGetColNDownlocks(matrix,col)         (matrix->ndownlocks[col])
#define SCIPmatrixGetVar(matrix,col)                   (matrix->vars[col])
#define SCIPmatrixGetColName(matrix,col)               (SCIPvarGetName(matrix->vars[col]))
#define SCIPmatrixGetRowValPtr(matrix,row)             (matrix->rowmatval + matrix->rowmatbeg[row])
#define SCIPmatrixGetRowIdxPtr(matrix,row)             (matrix->rowmatind + matrix->rowmatbeg[row])
#define SCIPmatrixGetRowNNonzs(matrix,row)             (matrix->rowmatcnt[row])
#define SCIPmatrixGetRowName(matrix,row)               (SCIPconsGetName(matrix->cons[row]))
#define SCIPmatrixGetNRows(matrix)                     (matrix->nrows)
#define SCIPmatrixGetRowLhs(matrix,row)                (matrix->lhs[row])
#define SCIPmatrixGetRowRhs(matrix,row)                (matrix->rhs[row])
#define SCIPmatrixIsRowRhsInfinity(matrix,row)         (matrix->isrhsinfinite[row])
#define SCIPmatrixGetNNonzs(matrix)                    (matrix->nnonzs)
#define SCIPmatrixGetRowMinActivity(matrix,row)        (matrix->minactivity[row])
#define SCIPmatrixGetRowMaxActivity(matrix,row)        (matrix->maxactivity[row])
#define SCIPmatrixGetRowNMinActNegInf(matrix,row)      (matrix->minactivityneginf[row])
#define SCIPmatrixGetRowNMinActPosInf(matrix,row)      (matrix->minactivityposinf[row])
#define SCIPmatrixGetRowNMaxActNegInf(matrix,row)      (matrix->maxactivityneginf[row])
#define SCIPmatrixGetRowNMaxActPosInf(matrix,row)      (matrix->maxactivityposinf[row])
#define SCIPmatrixGetCons(matrix,row)                  (matrix->cons[row])

#endif

/** initialize matrix by copying all check constraints
 *
 *  @note Completeness is checked by testing whether all check constraints are from a list of linear constraint handlers
 *        that can be represented.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPmatrixCreate(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX**         matrixptr,          /**< pointer to constraint matrix object to be initialized */
   SCIP_Bool             onlyifcomplete,     /**< should matrix creation be skipped if matrix will not be complete? */
   SCIP_Bool*            initialized,        /**< was the initialization successful? */
   SCIP_Bool*            complete,           /**< are all constraint represented within the matrix? */
   SCIP_Bool*            infeasible,         /**< pointer to return whether problem was detected to be infeasible during matrix creation */
   int*                  naddconss,          /**< pointer to count number of added (linear) constraints during matrix creation */
   int*                  ndelconss,          /**< pointer to count number of deleted specialized linear constraints during matrix creation */
   int*                  nchgcoefs,          /**< pointer to count number of changed coefficients during matrix creation */
   int*                  nchgbds,            /**< pointer to count number of changed bounds during matrix creation */
   int*                  nfixedvars          /**< pointer to count number of fixed variables during matrix creation */
   );

/** frees the constraint matrix */
SCIP_EXPORT
void SCIPmatrixFree(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX**         matrix              /**< constraint matrix object */
   );

/** print one row of the MIP matrix */
SCIP_EXPORT
void SCIPmatrixPrintRow(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row                 /**< row index */
   );

/** detect parallel rows, rhs/lhs are ignored */
SCIP_EXPORT
SCIP_RETCODE SCIPmatrixGetParallelRows(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real*            scale,              /**< scale factors of rows */
   int*                  pclass              /**< parallel row classes */
   );

/** removes the bounds of a column and updates the activities accordingly */
SCIP_EXPORT
void SCIPmatrixRemoveColumnBounds(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX*          matrix,             /**< constraint matrix */
   int                   col                 /**< column variable to remove bounds from */
   );

/** detect parallel columns, obj ignored */
SCIP_EXPORT
SCIP_RETCODE SCIPmatrixGetParallelCols(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real*            scale,              /**< scale factors of cols */
   int*                  pclass,             /**< parallel column classes */
   SCIP_Bool*            varineq             /**< indicating if variable is within an equation */
   );


#ifdef __cplusplus
}
#endif

#endif
