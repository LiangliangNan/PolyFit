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

/**@file   pub_lp.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for LP management
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_LP_H__
#define __SCIP_PUB_LP_H__


#include "lpi/type_lpi.h"
#include "scip/def.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_sepa.h"
#include "scip/type_var.h"
#include "scip/type_misc.h"

#ifdef NDEBUG
#include "scip/struct_lp.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup PublicColumnMethods
 *
 * @{
 */

/** sorts column entries such that LP rows precede non-LP rows and inside both parts lower row indices precede higher ones
 */
SCIP_EXPORT
void SCIPcolSort(
   SCIP_COL*             col                 /**< column to be sorted */
   );

/** gets objective value of column */
SCIP_EXPORT
SCIP_Real SCIPcolGetObj(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets lower bound of column */
SCIP_EXPORT
SCIP_Real SCIPcolGetLb(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets upper bound of column */
SCIP_EXPORT
SCIP_Real SCIPcolGetUb(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets best bound of column with respect to the objective function */
SCIP_EXPORT
SCIP_Real SCIPcolGetBestBound(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the primal LP solution of a column */
SCIP_EXPORT
SCIP_Real SCIPcolGetPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the minimal LP solution value, this column ever assumed */
SCIP_EXPORT
SCIP_Real SCIPcolGetMinPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the maximal LP solution value, this column ever assumed */
SCIP_EXPORT
SCIP_Real SCIPcolGetMaxPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the basis status of a column in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_ZERO for columns not in the current SCIP_LP
 */
SCIP_EXPORT
SCIP_BASESTAT SCIPcolGetBasisStatus(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets variable this column represents */
SCIP_EXPORT
SCIP_VAR* SCIPcolGetVar(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets unique index of col */
SCIP_EXPORT
int SCIPcolGetIndex(
   SCIP_COL*             col                 /**< LP col */
   );

/** gets probindex of corresponding variable */
SCIP_EXPORT
int SCIPcolGetVarProbindex(
   SCIP_COL*             col                 /**< LP col */
   );

/** returns whether the associated variable is of integral type (binary, integer, implicit integer) */
SCIP_EXPORT
SCIP_Bool SCIPcolIsIntegral(
   SCIP_COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is removable from the LP (due to aging or cleanup) */
SCIP_EXPORT
SCIP_Bool SCIPcolIsRemovable(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets position of column in current LP, or -1 if it is not in LP */
SCIP_EXPORT
int SCIPcolGetLPPos(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets depth in the tree where the column entered the LP, or -1 if it is not in LP */
SCIP_EXPORT
int SCIPcolGetLPDepth(
   SCIP_COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is member of current LP */
SCIP_EXPORT
SCIP_Bool SCIPcolIsInLP(
   SCIP_COL*             col                 /**< LP column */
   );

/** get number of nonzero entries in column vector */
SCIP_EXPORT
int SCIPcolGetNNonz(
   SCIP_COL*             col                 /**< LP column */
   );

/** get number of nonzero entries in column vector, that correspond to rows currently in the SCIP_LP;
 *
 *  @warning This method is only applicable on columns, that are completely linked to their rows (e.g. a column
 *  that is in the current LP and the LP was solved, or a column that was in a solved LP and didn't change afterwards
 */
SCIP_EXPORT
int SCIPcolGetNLPNonz(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets array with rows of nonzero entries */
SCIP_EXPORT
SCIP_ROW** SCIPcolGetRows(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets array with coefficients of nonzero entries */
SCIP_EXPORT
SCIP_Real* SCIPcolGetVals(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
SCIP_EXPORT
SCIP_Longint SCIPcolGetStrongbranchNode(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets number of times, strong branching was applied in current run on the given column */
SCIP_EXPORT
int SCIPcolGetNStrongbranchs(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the age of a column, i.e., the total number of successive times a column was in the LP and was 0.0 in the solution */
SCIP_EXPORT
int SCIPcolGetAge(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets opposite bound type of given bound type */
SCIP_EXPORT
SCIP_BOUNDTYPE SCIPboundtypeOpposite(
   SCIP_BOUNDTYPE        boundtype           /**< type of bound (lower or upper) */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcolGetObj(col)              (col)->obj
#define SCIPcolGetLb(col)               (col)->lb
#define SCIPcolGetUb(col)               (col)->ub
#define SCIPcolGetBestBound(col)        ((col)->obj >= 0.0 ? (col)->lb : (col)->ub)
#define SCIPcolGetPrimsol(col)          ((col)->lppos >= 0 ? (col)->primsol : 0.0)
#define SCIPcolGetMinPrimsol(col)       ((col)->minprimsol)
#define SCIPcolGetMaxPrimsol(col)       ((col)->maxprimsol)
#define SCIPcolGetBasisStatus(col)      ((SCIP_BASESTAT)(col)->basisstatus)
#define SCIPcolGetVar(col)              (col)->var
#define SCIPcolGetIndex(col)            (col)->index
#define SCIPcolIsIntegral(col)          (col)->integral
#define SCIPcolIsRemovable(col)         (col)->removable
#define SCIPcolGetLPPos(col)            (col)->lppos
#define SCIPcolGetLPDepth(col)          (col)->lpdepth
#define SCIPcolIsInLP(col)              ((col)->lppos >= 0)
#define SCIPcolGetNNonz(col)            (col)->len
#define SCIPcolGetNLPNonz(col)          (col)->nlprows
#define SCIPcolGetRows(col)             (col)->rows
#define SCIPcolGetVals(col)             (col)->vals
#define SCIPcolGetStrongbranchNode(col) (col)->sbnode
#define SCIPcolGetNStrongbranchs(col)   (col)->nsbcalls
#define SCIPcolGetAge(col)              (col)->age
#define SCIPboundtypeOpposite(boundtype) \
   ((boundtype) == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER)

#endif

/**@} */



/**@addtogroup PublicRowMethods
 *
 * @{
 */

/** comparison method for sorting rows by non-decreasing index */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIProwComp);

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
SCIP_EXPORT
void SCIProwLock(
   SCIP_ROW*             row                 /**< LP row */
   );

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
SCIP_EXPORT
void SCIProwUnlock(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the scalar product of the coefficient vectors of the two given rows */
SCIP_EXPORT
SCIP_Real SCIProwGetScalarProduct(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   );

/** returns the degree of parallelism between the hyperplanes defined by the two row vectors v, w:
 *  p = |v*w|/(|v|*|w|);
 *  the hyperplanes are parallel, iff p = 1, they are orthogonal, iff p = 0
 */
SCIP_EXPORT
SCIP_Real SCIProwGetParallelism(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   );

/** returns the degree of orthogonality between the hyperplanes defined by the two row vectors v, w:
 *  o = 1 - |v*w|/(|v|*|w|);
 *  the hyperplanes are orthogonal, iff p = 1, they are parallel, iff p = 0
 */
SCIP_EXPORT
SCIP_Real SCIProwGetOrthogonality(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   );

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
SCIP_EXPORT
void SCIProwSort(
   SCIP_ROW*             row                 /**< row to be sorted */
   );

/** get number of nonzero entries in row vector */
SCIP_EXPORT
int SCIProwGetNNonz(
   SCIP_ROW*             row                 /**< LP row */
   );

/** get number of nonzero entries in row vector, that correspond to columns currently in the SCIP_LP;
 *
 *  @warning This method is only applicable on rows, that are completely linked to their columns (e.g. a row
 *  that is in the current LP and the LP was solved, or a row that was in a solved LP and didn't change afterwards
 */
SCIP_EXPORT
int SCIProwGetNLPNonz(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets array with columns of nonzero entries */
SCIP_EXPORT
SCIP_COL** SCIProwGetCols(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets array with coefficients of nonzero entries */
SCIP_EXPORT
SCIP_Real* SCIProwGetVals(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets constant shift of row */
SCIP_EXPORT
SCIP_Real SCIProwGetConstant(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets Euclidean norm of row vector */
SCIP_EXPORT
SCIP_Real SCIProwGetNorm(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets sum norm of row vector (sum of absolute values of coefficients) */
SCIP_EXPORT
SCIP_Real SCIProwGetSumNorm(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the left hand side of the row */
SCIP_EXPORT
SCIP_Real SCIProwGetLhs(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the right hand side of the row */
SCIP_EXPORT
SCIP_Real SCIProwGetRhs(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the dual LP solution of a row */
SCIP_EXPORT
SCIP_Real SCIProwGetDualsol(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the dual Farkas coefficient of a row in an infeasible LP */
SCIP_EXPORT
SCIP_Real SCIProwGetDualfarkas(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the basis status of a row in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_BASIC for rows not in the current SCIP_LP
 */
SCIP_EXPORT
SCIP_BASESTAT SCIProwGetBasisStatus(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the name of the row */
SCIP_EXPORT
const char* SCIProwGetName(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets unique index of row */
SCIP_EXPORT
int SCIProwGetIndex(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets age of row */
SCIP_EXPORT
int SCIProwGetAge(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets rank of row */
SCIP_EXPORT
int SCIProwGetRank(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff the activity of the row (without the row's constant) is always integral in a feasible solution */
SCIP_EXPORT
SCIP_Bool SCIProwIsIntegral(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is only valid locally */
SCIP_EXPORT
SCIP_Bool SCIProwIsLocal(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
SCIP_EXPORT
SCIP_Bool SCIProwIsModifiable(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is removable from the LP (due to aging or cleanup) */
SCIP_EXPORT
SCIP_Bool SCIProwIsRemovable(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns type of origin that created the row */
SCIP_EXPORT
SCIP_ROWORIGINTYPE SCIProwGetOrigintype(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns origin constraint handler that created the row (NULL if not available) */
SCIP_EXPORT
SCIP_CONSHDLR* SCIProwGetOriginConshdlr(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns origin constraint that created the row (NULL if not available) */
SCIP_EXPORT
SCIP_CONS* SCIProwGetOriginCons(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns origin separator that created the row (NULL if not available) */
SCIP_EXPORT
SCIP_SEPA* SCIProwGetOriginSepa(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is member of the global cut pool */
SCIP_EXPORT
SCIP_Bool SCIProwIsInGlobalCutpool(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets position of row in current LP, or -1 if it is not in LP */
SCIP_EXPORT
int SCIProwGetLPPos(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets depth in the tree where the row entered the LP, or -1 if it is not in LP */
SCIP_EXPORT
int SCIProwGetLPDepth(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is member of current LP */
SCIP_EXPORT
SCIP_Bool SCIProwIsInLP(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the number of times that this row has been sharp in an optimal LP solution */
SCIP_EXPORT
SCIP_Longint SCIProwGetActiveLPCount(
   SCIP_ROW*             row                 /**< row */
   );

/** returns the number of LPs since this row has been created */
SCIP_EXPORT
SCIP_Longint SCIProwGetNLPsAfterCreation(
   SCIP_ROW*             row                 /**< row */
   );

/** changes the rank of LP row */
SCIP_EXPORT
void SCIProwChgRank(
   SCIP_ROW*             row,                /**< LP row */
   int                   rank                /**< new value for rank */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIProwGetNNonz(row)            (row)->len
#define SCIProwGetNLPNonz(row)          (row)->nlpcols
#define SCIProwGetCols(row)             (row)->cols
#define SCIProwGetVals(row)             (row)->vals
#define SCIProwGetConstant(row)         (row)->constant
#define SCIProwGetNorm(row)             sqrt((row)->sqrnorm)
#define SCIProwGetSumNorm(row)          (row)->sumnorm
#define SCIProwGetLhs(row)              (row)->lhs
#define SCIProwGetRhs(row)              (row)->rhs
#define SCIProwGetDualsol(row)          ((row)->lppos >= 0 ? (row)->dualsol : 0.0)
#define SCIProwGetDualfarkas(row)       ((row)->lppos >= 0 ? (row)->dualfarkas : 0.0)
#define SCIProwGetBasisStatus(row)      ((SCIP_BASESTAT) (row)->basisstatus)
#define SCIProwGetName(row)             (row)->name
#define SCIProwGetIndex(row)            (row)->index
#define SCIProwGetAge(row)              (row)->age
#define SCIProwGetRank(row)             (row)->rank
#define SCIProwIsIntegral(row)          (row)->integral
#define SCIProwIsLocal(row)             (row)->local
#define SCIProwIsModifiable(row)        (row)->modifiable
#define SCIProwIsRemovable(row)         (row)->removable
#define SCIProwGetOrigintype(row)       (row)->origintype
#define SCIProwGetOriginCons(row)       ((SCIP_CONS*) ((SCIP_ROWORIGINTYPE) row->origintype == SCIP_ROWORIGINTYPE_CONS ? (row)->origin : NULL))
#define SCIProwGetOriginSepa(row)       ((SCIP_SEPA*) ((SCIP_ROWORIGINTYPE) row->origintype == SCIP_ROWORIGINTYPE_SEPA ? (row)->origin : NULL))
#define SCIProwIsInGlobalCutpool(row)   (row)->inglobalcutpool
#define SCIProwGetLPPos(row)            (row)->lppos
#define SCIProwGetLPDepth(row)          (row)->lpdepth
#define SCIProwIsInLP(row)              ((row)->lppos >= 0)
#define SCIProwGetActiveLPCount(row)    ((row)->activeinlpcounter)
#define SCIProwGetNLPsAfterCreation(row) ((row)->nlpsaftercreation)
#define SCIProwChgRank(row, cutrank)    ((row)->rank = (cutrank))

#endif

/**@} */

#ifdef __cplusplus
}
#endif

#endif
