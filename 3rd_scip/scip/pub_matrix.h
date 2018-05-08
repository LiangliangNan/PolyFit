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
EXTERN
SCIP_Real* SCIPmatrixGetColValPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get column based start pointer of row indices */
EXTERN
int* SCIPmatrixGetColIdxPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get the number of non-zero entries of this column */
EXTERN
int SCIPmatrixGetColNNonzs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get number of columns of the matrix */
EXTERN
int SCIPmatrixGetNColumns(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   );

/** get upper bound of column */
EXTERN
SCIP_Real SCIPmatrixGetColUb(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get lower bound of column */
EXTERN
SCIP_Real SCIPmatrixGetColLb(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get number of uplocks of column */
EXTERN
int SCIPmatrixGetColNUplocks(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get number of downlocks of column */
EXTERN
int SCIPmatrixGetColNDownlocks(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get variable pointer of column */
EXTERN
SCIP_VAR* SCIPmatrixGetVar(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get name of column/variable */
EXTERN
const char* SCIPmatrixGetColName(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get row based start pointer of values */
EXTERN
SCIP_Real* SCIPmatrixGetRowValPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get row based start pointer of column indices */
EXTERN
int* SCIPmatrixGetRowIdxPtr(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of non-zeros of this row */
EXTERN
int SCIPmatrixGetRowNNonzs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get name of row */
EXTERN
const char* SCIPmatrixGetRowName(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of rows of the matrix */
EXTERN
int SCIPmatrixGetNRows(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   );

/** get left-hand-side of row */
EXTERN
SCIP_Real SCIPmatrixGetRowLhs(
   SCIP_MATRIX*          matrix,             /**< matrix instace */
   int                   row                 /**< row index */
   );

/** get right-hand-side of row */
EXTERN
SCIP_Real SCIPmatrixGetRowRhs(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** flag indicating if right-hand-side of row is infinity */
EXTERN
SCIP_Bool SCIPmatrixIsRowRhsInfinity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of non-zeros of matrix */
EXTERN
int SCIPmatrixGetNNonzs(
   SCIP_MATRIX*          matrix              /**< matrix instance */
   );

/** get minimal activity of row */
EXTERN
SCIP_Real SCIPmatrixGetRowMinActivity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get maximal activity of row */
EXTERN
SCIP_Real SCIPmatrixGetRowMaxActivity(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of negative infinities present within minimal activity */
EXTERN
int SCIPmatrixGetRowNMinActNegInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of positive infinities present within minimal activity */
EXTERN
int SCIPmatrixGetRowNMinActPosInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of negative infinities present within maximal activity */
EXTERN
int SCIPmatrixGetRowNMaxActNegInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get number of positive infinities present within maximal activity */
EXTERN
int SCIPmatrixGetRowNMaxActPosInf(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get constraint pointer for constraint representing row */
EXTERN
SCIP_CONS* SCIPmatrixGetCons(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   row                 /**< row index */
   );

/** get if conflicting uplocks of variable present */
EXTERN
SCIP_Bool SCIPmatrixUplockConflict(
   SCIP_MATRIX*          matrix,             /**< matrix instance */
   int                   col                 /**< column index */
   );

/** get if conflicting downlocks of variable present */
EXTERN
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
#define SCIPmatrixUplockConflict(matrix,col)           (SCIPvarGetNLocksUp(matrix->vars[col]) == matrix->nuplocks[col] ? FALSE : TRUE)
#define SCIPmatrixDownlockConflict(matrix,col)         (SCIPvarGetNLocksDown(matrix->vars[col]) == matrix->ndownlocks[col] ? FALSE : TRUE)

#endif

/** initialize matrix */
EXTERN
SCIP_RETCODE SCIPmatrixCreate(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_MATRIX**         matrixptr,          /**< pointer to constraint matrix object to be initialized */
   SCIP_Bool*            initialized,        /**< was the initialization successful? */
   SCIP_Bool*            complete            /**< are all constraint represented within the matrix? */
   );

/** frees the constraint matrix */
EXTERN
void SCIPmatrixFree(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX**         matrix              /**< constraint matrix object */
   );

/** print one row of the MIP matrix */
EXTERN
void SCIPmatrixPrintRow(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row                 /**< row index */
   );

/** detect parallel rows, rhs/lhs are ignored */
EXTERN
SCIP_RETCODE SCIPmatrixGetParallelRows(
   SCIP*                 scip,               /**< current SCIP instance */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   SCIP_Real*            scale,              /**< scale factors of rows */
   int*                  pclass              /**< parallel row classes */
   );

/** detect parallel rows, obj ignored */
EXTERN
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
