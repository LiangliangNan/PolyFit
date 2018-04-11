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

/**@file   pub_lp.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for LP management
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_LP_H__
#define __SCIP_PUB_LP_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "scip/type_sepa.h"
#include "scip/type_misc.h"
#include "lpi/type_lpi.h"

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
EXTERN
void SCIPcolSort(
   SCIP_COL*             col                 /**< column to be sorted */
   );

/** gets objective value of column */
EXTERN
SCIP_Real SCIPcolGetObj(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets lower bound of column */
EXTERN
SCIP_Real SCIPcolGetLb(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets upper bound of column */
EXTERN
SCIP_Real SCIPcolGetUb(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets best bound of column with respect to the objective function */
EXTERN
SCIP_Real SCIPcolGetBestBound(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the primal LP solution of a column */
EXTERN
SCIP_Real SCIPcolGetPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the minimal LP solution value, this column ever assumed */
EXTERN
SCIP_Real SCIPcolGetMinPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the maximal LP solution value, this column ever assumed */
EXTERN
SCIP_Real SCIPcolGetMaxPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the basis status of a column in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_ZERO for columns not in the current SCIP_LP
 */
EXTERN
SCIP_BASESTAT SCIPcolGetBasisStatus(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets variable this column represents */
EXTERN
SCIP_VAR* SCIPcolGetVar(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets unique index of col */
EXTERN
int SCIPcolGetIndex(
   SCIP_COL*             col                 /**< LP col */
   );

/** returns whether the associated variable is of integral type (binary, integer, implicit integer) */
EXTERN
SCIP_Bool SCIPcolIsIntegral(
   SCIP_COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is removable from the LP (due to aging or cleanup) */
EXTERN
SCIP_Bool SCIPcolIsRemovable(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets position of column in current LP, or -1 if it is not in LP */
EXTERN
int SCIPcolGetLPPos(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets depth in the tree where the column entered the LP, or -1 if it is not in LP */
EXTERN
int SCIPcolGetLPDepth(
   SCIP_COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is member of current LP */
EXTERN
SCIP_Bool SCIPcolIsInLP(
   SCIP_COL*             col                 /**< LP column */
   );

/** get number of nonzero entries in column vector */
EXTERN
int SCIPcolGetNNonz(
   SCIP_COL*             col                 /**< LP column */
   );

/** get number of nonzero entries in column vector, that correspond to rows currently in the SCIP_LP;
 *
 *  @warning This method is only applicable on columns, that are completely linked to their rows (e.g. a column
 *  that is in the current LP and the LP was solved, or a column that was in a solved LP and didn't change afterwards
 */
EXTERN
int SCIPcolGetNLPNonz(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets array with rows of nonzero entries */
EXTERN
SCIP_ROW** SCIPcolGetRows(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets array with coefficients of nonzero entries */
EXTERN
SCIP_Real* SCIPcolGetVals(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
EXTERN
SCIP_Longint SCIPcolGetStrongbranchNode(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets number of times, strong branching was applied in current run on the given column */
EXTERN
int SCIPcolGetNStrongbranchs(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets opposite bound type of given bound type */
EXTERN
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
#define SCIPboundtypeOpposite(boundtype) \
   ((boundtype) == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER)

#endif

/**@} */



/**@addtogroup PublicRowMethods
 *
 * @{
 */

/** comparison method for sorting rows by non-decreasing index */
EXTERN
SCIP_DECL_SORTPTRCOMP(SCIProwComp);

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
EXTERN
void SCIProwLock(
   SCIP_ROW*             row                 /**< LP row */
   );

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
EXTERN
void SCIProwUnlock(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the scalar product of the coefficient vectors of the two given rows */
EXTERN
SCIP_Real SCIProwGetScalarProduct(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   );

/** returns the degree of parallelism between the hyperplanes defined by the two row vectors v, w:
 *  p = |v*w|/(|v|*|w|);
 *  the hyperplanes are parallel, iff p = 1, they are orthogonal, iff p = 0
 */
EXTERN
SCIP_Real SCIProwGetParallelism(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   );

/** returns the degree of orthogonality between the hyperplanes defined by the two row vectors v, w:
 *  o = 1 - |v*w|/(|v|*|w|);
 *  the hyperplanes are orthogonal, iff p = 1, they are parallel, iff p = 0
 */
EXTERN
SCIP_Real SCIProwGetOrthogonality(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   );

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
EXTERN
void SCIProwSort(
   SCIP_ROW*             row                 /**< row to be sorted */
   );

/** get number of nonzero entries in row vector */
EXTERN
int SCIProwGetNNonz(
   SCIP_ROW*             row                 /**< LP row */
   );

/** get number of nonzero entries in row vector, that correspond to columns currently in the SCIP_LP;
 *
 *  @warning This method is only applicable on rows, that are completely linked to their columns (e.g. a row
 *  that is in the current LP and the LP was solved, or a row that was in a solved LP and didn't change afterwards
 */
EXTERN
int SCIProwGetNLPNonz(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets array with columns of nonzero entries */
EXTERN
SCIP_COL** SCIProwGetCols(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets array with coefficients of nonzero entries */
EXTERN
SCIP_Real* SCIProwGetVals(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets constant shift of row */
EXTERN
SCIP_Real SCIProwGetConstant(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets Euclidean norm of row vector */
EXTERN
SCIP_Real SCIProwGetNorm(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets sum norm of row vector (sum of absolute values of coefficients) */
EXTERN
SCIP_Real SCIProwGetSumNorm(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the left hand side of the row */
EXTERN
SCIP_Real SCIProwGetLhs(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the right hand side of the row */
EXTERN
SCIP_Real SCIProwGetRhs(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the dual LP solution of a row */
EXTERN
SCIP_Real SCIProwGetDualsol(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the dual Farkas coefficient of a row in an infeasible LP */
EXTERN
SCIP_Real SCIProwGetDualfarkas(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the basis status of a row in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_BASIC for rows not in the current SCIP_LP
 */
EXTERN
SCIP_BASESTAT SCIProwGetBasisStatus(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the name of the row */
EXTERN
const char* SCIProwGetName(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets unique index of row */
EXTERN
int SCIProwGetIndex(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets age of row */
EXTERN
int SCIProwGetAge(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets rank of row */
EXTERN
int SCIProwGetRank(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff the activity of the row (without the row's constant) is always integral in a feasible solution */
EXTERN
SCIP_Bool SCIProwIsIntegral(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is only valid locally */
EXTERN
SCIP_Bool SCIProwIsLocal(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
EXTERN
SCIP_Bool SCIProwIsModifiable(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is removable from the LP (due to aging or cleanup) */
EXTERN
SCIP_Bool SCIProwIsRemovable(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns type of origin that created the row */
EXTERN
SCIP_ROWORIGINTYPE SCIProwGetOrigintype(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns origin constraint handler that created the row (NULL if not available) */
EXTERN
SCIP_CONSHDLR* SCIProwGetOriginCons(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns origin separator that created the row (NULL if not available) */
EXTERN
SCIP_SEPA* SCIProwGetOriginSepa(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is member of the global cut pool */
EXTERN
SCIP_Bool SCIProwIsInGlobalCutpool(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets position of row in current LP, or -1 if it is not in LP */
EXTERN
int SCIProwGetLPPos(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets depth in the tree where the row entered the LP, or -1 if it is not in LP */
EXTERN
int SCIProwGetLPDepth(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is member of current LP */
EXTERN
SCIP_Bool SCIProwIsInLP(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the number of times that this row has been sharp in an optimal LP solution */
EXTERN
SCIP_Longint SCIProwGetActiveLPCount(
   SCIP_ROW*             row                 /**< row */
   );

/** returns the number of LPs since this row has been created */
EXTERN
SCIP_Longint SCIProwGetNLPsAfterCreation(
   SCIP_ROW*             row                 /**< row */
   );

/** changes the rank of LP row */
EXTERN
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
#define SCIProwGetOriginCons(row)       ((SCIP_CONSHDLR*) ((SCIP_ROWORIGINTYPE) row->origintype == SCIP_ROWORIGINTYPE_CONS ? (row)->origin : NULL))
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
