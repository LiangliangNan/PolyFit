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

/**@file   history.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for branching and inference history
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HISTORY_H__
#define __SCIP_HISTORY_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_history.h"

#ifdef NDEBUG
#include "scip/struct_history.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** creates an empty history entry */
extern
SCIP_RETCODE SCIPhistoryCreate(
   SCIP_HISTORY**        history,            /**< pointer to store branching and inference history */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** frees a history entry */
extern
void SCIPhistoryFree(
   SCIP_HISTORY**        history,            /**< pointer to branching and inference history */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** resets history entry to zero */
extern
void SCIPhistoryReset(
   SCIP_HISTORY*         history             /**< branching and inference history */
   );

/** unites two history entries by adding the values of the second one to the first one */
extern
void SCIPhistoryUnite(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_HISTORY*         addhistory,         /**< history values to add to history */
   SCIP_Bool             switcheddirs        /**< should the history entries be united with switched directories */
   );

/** updates the pseudo costs for a change of "solvaldelta" in the variable's LP solution value and a change of "objdelta"
 *  in the LP's objective value
 */
extern
void SCIPhistoryUpdatePseudocost(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight of this update in pseudo cost sum (added to pscostcount) */
   );


/**@defgroup ValueHistory Value Based History
 * @ingroup INTERNALAPI
 * @brief Value based history methods
 *
 * @{
 */

/** creates an empty value history */
extern
SCIP_RETCODE SCIPvaluehistoryCreate(
   SCIP_VALUEHISTORY**   valuehistory,       /**< pointer to store the value based branching and inference histories */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** frees a value history */
extern
void SCIPvaluehistoryFree(
   SCIP_VALUEHISTORY**   valuehistory,       /**< pointer to value based history */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** finds for the given domain value the history if it does not exist yet it will be created */
extern
SCIP_RETCODE SCIPvaluehistoryFind(
   SCIP_VALUEHISTORY*    valuehistory,       /**< value based history */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             value,              /**< domain value of interest */
   SCIP_HISTORY**        history             /**< pointer to store the history for the given domain value */
   );

/** scales the conflict score values with the given scalar for each value history entry */
extern
void SCIPvaluehistoryScaleVSIDS(
   SCIP_VALUEHISTORY*    valuehistory,       /**< value based history */
   SCIP_Real             scalar              /**< scalar to multiply the conflict scores with */
   );

/**@} */

/** returns the opposite direction of the given branching direction */
extern
SCIP_BRANCHDIR SCIPbranchdirOpposite(
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the expected dual gain for moving the corresponding variable by "solvaldelta" */
extern
SCIP_Real SCIPhistoryGetPseudocost(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** returns the variance of pseudo costs about the mean. */
extern
SCIP_Real SCIPhistoryGetPseudocostVariance(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        direction           /**< direction of variable: 1 for upwards history, 0 for downwards history */
   );

/** returns the (possible fractional) number of (partial) pseudo cost updates performed on this pseudo cost entry in 
 *  the given branching direction
 */
extern
SCIP_Real SCIPhistoryGetPseudocostCount(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns whether the pseudo cost entry is empty in the given branching direction (whether no value was added yet) */
extern
SCIP_Bool SCIPhistoryIsPseudocostEmpty(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** increases the conflict score of the history entry by the given weight */
extern
void SCIPhistoryIncVSIDS(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             weight              /**< weight of this update in conflict score */
   );

 /** scales the conflict score values with the given scalar */
extern
void SCIPhistoryScaleVSIDS(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_Real             scalar              /**< scalar to multiply the conflict scores with */
   );

/** increases the number of active conflicts by one and the overall length of the history entry by the given weight */
extern
void SCIPhistoryIncNActiveConflicts(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_Real             length              /**< length of the conflict */
   );

/** gets the number of active conflicts of the history entry */
extern
SCIP_Longint SCIPhistoryGetNActiveConflicts(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   );

/** gets the average conflict length of the history entry */
extern
SCIP_Real SCIPhistoryGetAvgConflictlength(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction */
   );

/** increases the number of branchings counter */
extern
void SCIPhistoryIncNBranchings(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   int                   depth               /**< depth at which the bound change took place */
   );

/** increases the number of inferences counter */
extern
void SCIPhistoryIncInferenceSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             weight              /**< weight of this update in cutoff score */
   );


/** increases the number of cutoffs counter */
extern
void SCIPhistoryIncCutoffSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Real             weight              /**< weight of this update in cutoff score */
   );

/** get number of branchings counter */
extern
SCIP_Longint SCIPhistoryGetNBranchings(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** get number of inferences counter */
extern
SCIP_Real SCIPhistoryGetInferenceSum(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of inferences per branching */
extern
SCIP_Real SCIPhistoryGetAvgInferences(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average number of cutoffs per branching */
extern
SCIP_Real SCIPhistoryGetAvgCutoffs(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

/** returns the average depth of bound changes due to branching */
extern
SCIP_Real SCIPhistoryGetAvgBranchdepth(
   SCIP_HISTORY*         history,            /**< branching and inference history */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPbranchdirOpposite(dir)                                      \
   ((dir) == SCIP_BRANCHDIR_DOWNWARDS ? SCIP_BRANCHDIR_UPWARDS          \
      : ((dir) == SCIP_BRANCHDIR_UPWARDS ? SCIP_BRANCHDIR_DOWNWARDS : SCIP_BRANCHDIR_AUTO))
#define SCIPhistoryGetPseudocost(history,solvaldelta)                   \
   ( (solvaldelta) >= 0.0 ? (solvaldelta) * ((history)->pscostcount[1] > 0.0 \
      ? (history)->pscostweightedmean[1] : 1.0)      \
      : -(solvaldelta) * ((history)->pscostcount[0] > 0.0               \
         ? (history)->pscostweightedmean[0] : 1.0) )
#define SCIPhistoryGetPseudocostVariance(history, dir)                  \
   ( (history)->pscostcount[dir] >= 1.9 ? 1 / ((history)->pscostcount[dir] - 1)  \
         * ((history)->pscostvariance[dir]) \
         : 0.0)
#define SCIPhistoryGetPseudocostCount(history,dir) ((history)->pscostcount[dir])
#define SCIPhistoryIsPseudocostEmpty(history,dir)  ((history)->pscostcount[dir] == 0.0)
#define SCIPhistoryIncVSIDS(history,dir,weight) (history)->vsids[dir] += (weight)
#define SCIPhistoryScaleVSIDS(history,scalar)  { (history)->vsids[0] *= (scalar); \
      (history)->vsids[1] *= (scalar);  }
#define SCIPhistoryIncNActiveConflicts(history,dir,length) { (history)->nactiveconflicts[dir]++; \
      (history)->conflengthsum[dir] += length; }
#define SCIPhistoryGetNActiveConflicts(history,dir) ((history)->nactiveconflicts[dir])
#define SCIPhistoryGetAvgConflictlength(history,dir) ((history)->conflengthsum[dir] > 0.0 \
      ? (SCIP_Real)(history)->nactiveconflicts[dir]/(SCIP_Real)(history)->conflengthsum[dir] : 0.0)
#define SCIPhistoryIncNBranchings(history,dir,depth) { (history)->nbranchings[dir]++; \
      (history)->branchdepthsum[dir] += depth; }
#define SCIPhistoryIncInferenceSum(history,dir,weight)     (history)->inferencesum[dir] += (weight)
#define SCIPhistoryIncCutoffSum(history,dir,weight)        (history)->cutoffsum[dir] += (weight)
#define SCIPhistoryGetNBranchings(history,dir)     ((history)->nbranchings[dir])
#define SCIPhistoryGetInferenceSum(history,dir)     ((history)->inferencesum[dir])
#define SCIPhistoryGetAvgInferences(history,dir)   ((history)->nbranchings[dir] > 0 \
      ? (SCIP_Real)(history)->inferencesum[dir]/(SCIP_Real)(history)->nbranchings[dir] : 0.0)
#define SCIPhistoryGetCutoffSum(history,dir)        ((history)->cutoffsum[dir])
#define SCIPhistoryGetAvgCutoffs(history,dir)      ((history)->nbranchings[dir] > 0 \
      ? (SCIP_Real)(history)->cutoffsum[dir]/(SCIP_Real)(history)->nbranchings[dir] : 0.0)
#define SCIPhistoryGetAvgBranchdepth(history,dir)  ((history)->nbranchings[dir] > 0 \
      ? (SCIP_Real)(history)->branchdepthsum[dir]/(SCIP_Real)(history)->nbranchings[dir] : 1.0)

#endif

#ifdef __cplusplus
}
#endif

#endif
