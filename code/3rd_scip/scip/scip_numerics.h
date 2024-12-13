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

/**@file   scip_numerics.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for numerical tolerances
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_NUMERICS_H__
#define __SCIP_SCIP_NUMERICS_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/set.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicToleranceMethods
 *
 * @{
 */

/** returns value treated as zero
 *
 *  @return value treated as zero
 */
SCIP_EXPORT
SCIP_Real SCIPepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns value treated as zero for sums of floating point values
 *
 *  @return value treated as zero for sums of floating point values
 */
SCIP_EXPORT
SCIP_Real SCIPsumepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns feasibility tolerance for constraints
 *
 *  @return feasibility tolerance for constraints
 */
SCIP_EXPORT
#ifdef __GNUC__
__attribute__ ((pure))
#endif
SCIP_Real SCIPfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns primal feasibility tolerance of LP solver
 *
 *  @deprecated Please use SCIPgetLPFeastol().
 *
 *  @return primal feasibility tolerance of LP solver
 */
SCIP_DEPRECATED
SCIP_EXPORT
SCIP_Real SCIPlpfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns feasibility tolerance for reduced costs
 *
 *  @return feasibility tolerance for reduced costs
 */
SCIP_EXPORT
#ifdef __GNUC__
__attribute__ ((pure))
#endif
SCIP_Real SCIPdualfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns convergence tolerance used in barrier algorithm
 *
 *  @return convergence tolerance used in barrier algorithm
 */
SCIP_EXPORT
SCIP_Real SCIPbarrierconvtol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return the cutoff bound delta
 *
 *  @return cutoff bound data
 */
SCIP_EXPORT
SCIP_Real SCIPcutoffbounddelta(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return the relaxation primal feasibility tolerance
 *
 *  @see SCIPchgRelaxfeastol
 *  @return relaxfeastol
 */
SCIP_EXPORT
SCIP_Real SCIPrelaxfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the feasibility tolerance for constraints
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgFeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             feastol             /**< new feasibility tolerance for constraints */
   );

/** sets the primal feasibility tolerance of LP solver
 *
 *  @deprecated Please use SCIPsetLPFeastol().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_DEPRECATED
SCIP_RETCODE SCIPchgLpfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpfeastol,          /**< new primal feasibility tolerance of LP solver */
   SCIP_Bool             printnewvalue       /**< should "numerics/lpfeastol = ..." be printed? */
   );

/** sets the feasibility tolerance for reduced costs
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgDualfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dualfeastol         /**< new feasibility tolerance for reduced costs */
   );

/** sets the convergence tolerance used in barrier algorithm
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPchgBarrierconvtol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   );

/** sets the primal feasibility tolerance of relaxations
 *
 * This tolerance value is used by the SCIP core and plugins to tighten then feasibility tolerance on relaxations
 * (especially the LP relaxation) during a solve. It is set to SCIP_INVALID initially, which means that only the
 * feasibility tolerance of the particular relaxation is taken into account. If set to a valid value, however,
 * then this value should be used to reduce the primal feasibility tolerance of a relaxation (thus, use the
 * minimum of relaxfeastol and the relaxations primal feastol).
 *
 * @pre The value of relaxfeastol is reset to SCIP_INVALID when initializing the solve (INITSOL).
 * Therefore, this method can only be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 * @return previous value of relaxfeastol
 */
SCIP_EXPORT
SCIP_Real SCIPchgRelaxfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             relaxfeastol        /**< new primal feasibility tolerance of relaxations */
   );

/** marks that some limit parameter was changed */
SCIP_EXPORT
void SCIPmarkLimitChanged(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns value treated as infinity */
SCIP_EXPORT
SCIP_Real SCIPinfinity(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the minimum value that is regarded as huge and should be handled separately (e.g., in activity
 *  computation)
 */
SCIP_EXPORT
SCIP_Real SCIPgetHugeValue(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** checks, if values are in range of epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) lower than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) greater than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) greater than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) lower than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is (positive) infinite */
SCIP_EXPORT
SCIP_Bool SCIPisInfinity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to be compared against infinity */
   );

/** checks, if value is huge and should be handled separately (e.g., in activity computation) */
SCIP_EXPORT
SCIP_Bool SCIPisHugeValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to be checked whether it is huge */
   );

/** checks, if value is in range epsilon of 0.0 */
SCIP_EXPORT
SCIP_Bool SCIPisZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is greater than epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is lower than -epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is integral within epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
SCIP_EXPORT
SCIP_Bool SCIPisScalingIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val,                /**< unscaled value to check for scaled integrality */
   SCIP_Real             scalar              /**< value to scale val with for checking for integrality */
   );

/** checks, if given fractional part is smaller than epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value + epsilon down to the next integer */
SCIP_EXPORT
SCIP_Real SCIPfloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value - epsilon up to the next integer */
SCIP_EXPORT
SCIP_Real SCIPceil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value to the nearest integer with epsilon tolerance */
SCIP_EXPORT
SCIP_Real SCIPround(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
SCIP_EXPORT
SCIP_Real SCIPfrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to return fractional part for */
   );

/** checks, if values are in range of sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) lower than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisSumLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisSumLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) greater than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisSumGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
SCIP_EXPORT
SCIP_Bool SCIPisSumGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range sumepsilon of 0.0 */
SCIP_EXPORT
SCIP_Bool SCIPisSumZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is greater than sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is lower than -sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if relative difference of values is in range of feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisFeasEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference val1 and val2 is lower than feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisFeasLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisFeasLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than feastol */
SCIP_EXPORT
SCIP_Bool SCIPisFeasGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
SCIP_EXPORT
SCIP_Bool SCIPisFeasGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range feasibility tolerance of 0.0 */
SCIP_EXPORT
SCIP_Bool SCIPisFeasZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is greater than feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisFeasPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is lower than -feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisFeasNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is integral within the LP feasibility bounds */
SCIP_EXPORT
SCIP_Bool SCIPisFeasIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if given fractional part is smaller than feastol */
SCIP_EXPORT
SCIP_Bool SCIPisFeasFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value + feasibility tolerance down to the next integer */
SCIP_EXPORT
SCIP_Real SCIPfeasFloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value - feasibility tolerance up to the next integer */
SCIP_EXPORT
SCIP_Real SCIPfeasCeil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value to the nearest integer in feasibility tolerance */
SCIP_EXPORT
SCIP_Real SCIPfeasRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** returns fractional part of value, i.e. x - floor(x) */
SCIP_EXPORT
SCIP_Real SCIPfeasFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if relative difference of values is in range of dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference val1 and val2 is lower than dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range dual feasibility tolerance of 0.0 */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is greater than dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is lower than -dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if value is integral within the LP dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if given fractional part is smaller than dual feasibility tolerance */
SCIP_EXPORT
SCIP_Bool SCIPisDualfeasFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value + dual feasibility tolerance down to the next integer */
SCIP_EXPORT
SCIP_Real SCIPdualfeasFloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value - dual feasibility tolerance up to the next integer */
SCIP_EXPORT
SCIP_Real SCIPdualfeasCeil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** rounds value to the nearest integer in dual feasibility tolerance */
SCIP_EXPORT
SCIP_Real SCIPdualfeasRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** returns fractional part of value, i.e. x - floor(x) in dual feasibility tolerance */
SCIP_EXPORT
SCIP_Real SCIPdualfeasFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   );

/** checks, if the given new lower bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
SCIP_EXPORT
SCIP_Bool SCIPisLbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   );

/** checks, if the given new upper bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
SCIP_EXPORT
SCIP_Bool SCIPisUbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   );

/** checks, if relative difference of values is in range of epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
SCIP_EXPORT
SCIP_Bool SCIPisRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of values is in range of sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/**! [SnippetCodeStyleNaming] */

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
SCIP_EXPORT
SCIP_Bool SCIPisSumRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

/** converts the given real number representing an integer to an int; in optimized mode the function gets inlined for
 *  performance; in debug mode we check some additional conditions
 */
SCIP_EXPORT
int SCIPconvertRealToInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             real                /**< double bound to convert */
   );

/**! [SnippetCodeStyleNaming] */

/** converts the given real number representing an integer to a long integer; in optimized mode the function gets inlined for
 *  performance; in debug mode we check some additional conditions
 */
SCIP_EXPORT
SCIP_Longint SCIPconvertRealToLongint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             real                /**< double bound to convert */
   );

/** Checks, if an iteratively updated value is reliable or should be recomputed from scratch.
 *  This is useful, if the value, e.g., the activity of a linear constraint or the pseudo objective value, gets a high
 *  absolute value during the optimization process which is later reduced significantly. In this case, the last digits
 *  were canceled out when increasing the value and are random after decreasing it.
 *  We do not consider the cancellations which can occur during increasing the absolute value because they just cannot
 *  be expressed using fixed precision floating point arithmetic, anymore.
 *  In order to get more reliable values, the idea is to always store the last reliable value, where increasing the
 *  absolute of the value is viewed as preserving reliability. Then, after each update, the new absolute value can be
 *  compared against the last reliable one with this method, checking whether it was decreased by a factor of at least
 *  "lp/recompfac" and should be recomputed.
 */
SCIP_EXPORT
SCIP_Bool SCIPisUpdateUnreliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newvalue,           /**< new value after update */
   SCIP_Real             oldvalue            /**< old value, i.e., last reliable value */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPinfinity(scip)                        SCIPsetInfinity((scip)->set)
#define SCIPisInfinity(scip, val)                 SCIPsetIsInfinity((scip)->set, val)
#define SCIPisHugeValue(scip, val)                SCIPsetIsHugeValue((scip)->set, val)
#define SCIPgetHugeValue(scip)                    SCIPsetGetHugeValue((scip)->set)
#define SCIPisEQ(scip, val1, val2)                SCIPsetIsEQ((scip)->set, val1, val2)
#define SCIPisLT(scip, val1, val2)                SCIPsetIsLT((scip)->set, val1, val2)
#define SCIPisLE(scip, val1, val2)                SCIPsetIsLE((scip)->set, val1, val2)
#define SCIPisGT(scip, val1, val2)                SCIPsetIsGT((scip)->set, val1, val2)
#define SCIPisGE(scip, val1, val2)                SCIPsetIsGE((scip)->set, val1, val2)
#define SCIPisZero(scip, val)                     SCIPsetIsZero((scip)->set, val)
#define SCIPisPositive(scip, val)                 SCIPsetIsPositive((scip)->set, val)
#define SCIPisNegative(scip, val)                 SCIPsetIsNegative((scip)->set, val)
#define SCIPisIntegral(scip, val)                 SCIPsetIsIntegral((scip)->set, val)
#define SCIPisScalingIntegral(scip, val, scalar)  SCIPsetIsScalingIntegral((scip)->set, val, scalar)
#define SCIPisFracIntegral(scip, val)             SCIPsetIsFracIntegral((scip)->set, val)
#define SCIPfloor(scip, val)                      SCIPsetFloor((scip)->set, val)
#define SCIPceil(scip, val)                       SCIPsetCeil((scip)->set, val)
#define SCIPround(scip, val)                      SCIPsetRound((scip)->set, val)
#define SCIPfrac(scip, val)                       SCIPsetFrac((scip)->set, val)

#define SCIPisSumEQ(scip, val1, val2)             SCIPsetIsSumEQ((scip)->set, val1, val2)
#define SCIPisSumLT(scip, val1, val2)             SCIPsetIsSumLT((scip)->set, val1, val2)
#define SCIPisSumLE(scip, val1, val2)             SCIPsetIsSumLE((scip)->set, val1, val2)
#define SCIPisSumGT(scip, val1, val2)             SCIPsetIsSumGT((scip)->set, val1, val2)
#define SCIPisSumGE(scip, val1, val2)             SCIPsetIsSumGE((scip)->set, val1, val2)
#define SCIPisSumZero(scip, val)                  SCIPsetIsSumZero((scip)->set, val)
#define SCIPisSumPositive(scip, val)              SCIPsetIsSumPositive((scip)->set, val)
#define SCIPisSumNegative(scip, val)              SCIPsetIsSumNegative((scip)->set, val)

#define SCIPisFeasEQ(scip, val1, val2)            SCIPsetIsFeasEQ((scip)->set, val1, val2)
#define SCIPisFeasLT(scip, val1, val2)            SCIPsetIsFeasLT((scip)->set, val1, val2)
#define SCIPisFeasLE(scip, val1, val2)            SCIPsetIsFeasLE((scip)->set, val1, val2)
#define SCIPisFeasGT(scip, val1, val2)            SCIPsetIsFeasGT((scip)->set, val1, val2)
#define SCIPisFeasGE(scip, val1, val2)            SCIPsetIsFeasGE((scip)->set, val1, val2)
#define SCIPisFeasZero(scip, val)                 SCIPsetIsFeasZero((scip)->set, val)
#define SCIPisFeasPositive(scip, val)             SCIPsetIsFeasPositive((scip)->set, val)
#define SCIPisFeasNegative(scip, val)             SCIPsetIsFeasNegative((scip)->set, val)
#define SCIPisFeasIntegral(scip, val)             SCIPsetIsFeasIntegral((scip)->set, val)
#define SCIPisFeasFracIntegral(scip, val)         SCIPsetIsFeasFracIntegral((scip)->set, val)
#define SCIPfeasFloor(scip, val)                  SCIPsetFeasFloor((scip)->set, val)
#define SCIPfeasCeil(scip, val)                   SCIPsetFeasCeil((scip)->set, val)
#define SCIPfeasRound(scip, val)                  SCIPsetFeasRound((scip)->set, val)
#define SCIPfeasFrac(scip, val)                   SCIPsetFeasFrac((scip)->set, val)

#define SCIPisDualfeasEQ(scip, val1, val2)        SCIPsetIsDualfeasEQ((scip)->set, val1, val2)
#define SCIPisDualfeasLT(scip, val1, val2)        SCIPsetIsDualfeasLT((scip)->set, val1, val2)
#define SCIPisDualfeasLE(scip, val1, val2)        SCIPsetIsDualfeasLE((scip)->set, val1, val2)
#define SCIPisDualfeasGT(scip, val1, val2)        SCIPsetIsDualfeasGT((scip)->set, val1, val2)
#define SCIPisDualfeasGE(scip, val1, val2)        SCIPsetIsDualfeasGE((scip)->set, val1, val2)
#define SCIPisDualfeasZero(scip, val)             SCIPsetIsDualfeasZero((scip)->set, val)
#define SCIPisDualfeasPositive(scip, val)         SCIPsetIsDualfeasPositive((scip)->set, val)
#define SCIPisDualfeasNegative(scip, val)         SCIPsetIsDualfeasNegative((scip)->set, val)
#define SCIPisDualfeasIntegral(scip, val)         SCIPsetIsDualfeasIntegral((scip)->set, val)
#define SCIPisDualfeasFracIntegral(scip, val)     SCIPsetIsDualfeasFracIntegral((scip)->set, val)
#define SCIPdualfeasFloor(scip, val)              SCIPsetDualfeasFloor((scip)->set, val)
#define SCIPdualfeasCeil(scip, val)               SCIPsetDualfeasCeil((scip)->set, val)
#define SCIPdualfeasRound(scip, val)              SCIPsetDualfeasRound((scip)->set, val)
#define SCIPdualfeasFrac(scip, val)               SCIPsetDualfeasFrac((scip)->set, val)

#define SCIPisLbBetter(scip, newlb, oldlb, oldub) SCIPsetIsLbBetter(scip->set, newlb, oldlb, oldub)
#define SCIPisUbBetter(scip, newub, oldlb, oldub) SCIPsetIsUbBetter(scip->set, newub, oldlb, oldub)

#define SCIPisRelEQ(scip, val1, val2)             SCIPsetIsRelEQ((scip)->set, val1, val2)
#define SCIPisRelLT(scip, val1, val2)             SCIPsetIsRelLT((scip)->set, val1, val2)
#define SCIPisRelLE(scip, val1, val2)             SCIPsetIsRelLE((scip)->set, val1, val2)
#define SCIPisRelGT(scip, val1, val2)             SCIPsetIsRelGT((scip)->set, val1, val2)
#define SCIPisRelGE(scip, val1, val2)             SCIPsetIsRelGE((scip)->set, val1, val2)

#define SCIPisSumRelEQ(scip, val1, val2)          SCIPsetIsSumRelEQ((scip)->set, val1, val2)
#define SCIPisSumRelLT(scip, val1, val2)          SCIPsetIsSumRelLT((scip)->set, val1, val2)
#define SCIPisSumRelLE(scip, val1, val2)          SCIPsetIsSumRelLE((scip)->set, val1, val2)
#define SCIPisSumRelGT(scip, val1, val2)          SCIPsetIsSumRelGT((scip)->set, val1, val2)
#define SCIPisSumRelGE(scip, val1, val2)          SCIPsetIsSumRelGE((scip)->set, val1, val2)
#define SCIPconvertRealToInt(scip, real)          ((int)((real) < 0 ? ((real) - 0.5) : ((real) + 0.5)))
#define SCIPconvertRealToLongint(scip, real)      ((SCIP_Longint)((real) < 0 ? ((real) - 0.5) : ((real) + 0.5)))

#define SCIPisUpdateUnreliable(scip, newval, oldval) SCIPsetIsUpdateUnreliable((scip)->set, newval, oldval)

#endif

/** outputs a real number, or "+infinity", or "-infinity" to a file */
SCIP_EXPORT
void SCIPprintReal(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             val,                /**< value to print */
   int                   width,              /**< width of the field */
   int                   precision           /**< number of significant digits printed */
   );

/** parse a real value that was written with SCIPprintReal() */
SCIP_EXPORT
SCIP_Bool SCIPparseReal(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
