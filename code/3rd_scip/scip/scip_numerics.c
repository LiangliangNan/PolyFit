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

/**@file   scip_numerics.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for numerical tolerances
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/scip_numerics.h"
#include "scip/set.h"
#include "scip/struct_lp.h"
#include "scip/struct_scip.h"
#include "scip/scip_lp.h"
#include "scip/scip_message.h"
#include <string.h>


/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPinfinity
#undef SCIPisInfinity
#undef SCIPisEQ
#undef SCIPisLT
#undef SCIPisLE
#undef SCIPisGT
#undef SCIPisGE
#undef SCIPisZero
#undef SCIPisPositive
#undef SCIPisNegative
#undef SCIPisIntegral
#undef SCIPisScalingIntegral
#undef SCIPisFracIntegral
#undef SCIPfloor
#undef SCIPceil
#undef SCIPround
#undef SCIPfrac
#undef SCIPisSumEQ
#undef SCIPisSumLT
#undef SCIPisSumLE
#undef SCIPisSumGT
#undef SCIPisSumGE
#undef SCIPisSumZero
#undef SCIPisSumPositive
#undef SCIPisSumNegative
#undef SCIPisFeasEQ
#undef SCIPisFeasLT
#undef SCIPisFeasLE
#undef SCIPisFeasGT
#undef SCIPisFeasGE
#undef SCIPisFeasZero
#undef SCIPisFeasPositive
#undef SCIPisFeasNegative
#undef SCIPisFeasIntegral
#undef SCIPisFeasFracIntegral
#undef SCIPfeasFloor
#undef SCIPfeasCeil
#undef SCIPfeasRound
#undef SCIPfeasFrac
#undef SCIPisDualfeasEQ
#undef SCIPisDualfeasLT
#undef SCIPisDualfeasLE
#undef SCIPisDualfeasGT
#undef SCIPisDualfeasGE
#undef SCIPisDualfeasZero
#undef SCIPisDualfeasPositive
#undef SCIPisDualfeasNegative
#undef SCIPisDualfeasIntegral
#undef SCIPisDualfeasFracIntegral
#undef SCIPdualfeasFloor
#undef SCIPdualfeasCeil
#undef SCIPdualfeasRound
#undef SCIPdualfeasFrac
#undef SCIPisLbBetter
#undef SCIPisUbBetter
#undef SCIPisRelEQ
#undef SCIPisRelLT
#undef SCIPisRelLE
#undef SCIPisRelGT
#undef SCIPisRelGE
#undef SCIPisSumRelEQ
#undef SCIPisSumRelLT
#undef SCIPisSumRelLE
#undef SCIPisSumRelGT
#undef SCIPisSumRelGE
#undef SCIPconvertRealToInt
#undef SCIPconvertRealToLongint
#undef SCIPisUpdateUnreliable
#undef SCIPisHugeValue
#undef SCIPgetHugeValue

/** returns value treated as zero
 *
 *  @return value treated as zero
 */
SCIP_Real SCIPepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetEpsilon(scip->set);
}

/** returns value treated as zero for sums of floating point values
 *
 *  @return value treated as zero for sums of floating point values
 */
SCIP_Real SCIPsumepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetSumepsilon(scip->set);
}

/** returns feasibility tolerance for constraints
 *
 *  @return feasibility tolerance for constraints
 */
SCIP_Real SCIPfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeastol(scip->set);
}

/** returns primal feasibility tolerance of LP solver
 *
 *  @deprecated Please use SCIPgetLPFeastol().
 *
 *  @return primal feasibility tolerance of LP solver
 */
SCIP_Real SCIPlpfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPgetLPFeastol(scip);
}

/** returns feasibility tolerance for reduced costs
 *
 *  @return feasibility tolerance for reduced costs
 */
SCIP_Real SCIPdualfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetDualfeastol(scip->set);
}

/** returns convergence tolerance used in barrier algorithm
 *
 *  @return convergence tolerance used in barrier algorithm
 */
SCIP_Real SCIPbarrierconvtol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetBarrierconvtol(scip->set);
}

/** return the cutoff bound delta
 *
 *  @return cutoff bound data
 */
SCIP_Real SCIPcutoffbounddelta(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetCutoffbounddelta(scip->set);
}

/** return the relaxation primal feasibility tolerance
 *
 *  @see SCIPchgRelaxfeastol
 *  @return relaxfeastol
 */
SCIP_Real SCIPrelaxfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetRelaxfeastol(scip->set);
}

/** sets the feasibility tolerance for constraints
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgFeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             feastol             /**< new feasibility tolerance for constraints */
   )
{
   assert(scip != NULL);

   /* change the settings */
   SCIP_CALL( SCIPsetSetFeastol(scip->set, scip->lp, feastol) );

   return SCIP_OKAY;
}

/** sets the primal feasibility tolerance of LP solver
 *
 *  @deprecated Please use SCIPsetLPFeastol().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgLpfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lpfeastol,          /**< new primal feasibility tolerance of LP solver */
   SCIP_Bool             printnewvalue       /**< should "numerics/lpfeastol = ..." be printed? */
   )
{
   SCIPsetLPFeastol(scip, lpfeastol);

   if( printnewvalue )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "numerics/lpfeastol = %.15g\n", lpfeastol);
   }

   return SCIP_OKAY;
}

/** sets the feasibility tolerance for reduced costs
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgDualfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dualfeastol         /**< new feasibility tolerance for reduced costs */
   )
{
   assert(scip != NULL);

   /* mark the LP unsolved, if the dual feasibility tolerance was tightened */
   if( scip->lp != NULL && dualfeastol < SCIPsetDualfeastol(scip->set) )
   {
      scip->lp->solved = FALSE;
      scip->lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
   }

   /* change the settings */
   SCIP_CALL( SCIPsetSetDualfeastol(scip->set, dualfeastol) );

   return SCIP_OKAY;
}

/** sets the convergence tolerance used in barrier algorithm
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
SCIP_RETCODE SCIPchgBarrierconvtol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   )
{
   assert(scip != NULL);

   /* mark the LP unsolved, if the convergence tolerance was tightened, and the LP was solved with the barrier algorithm */
   if( scip->lp != NULL && barrierconvtol < SCIPsetBarrierconvtol(scip->set)
      && (scip->lp->lastlpalgo == SCIP_LPALGO_BARRIER || scip->lp->lastlpalgo == SCIP_LPALGO_BARRIERCROSSOVER) )
      scip->lp->solved = FALSE;

   /* change the settings */
   SCIP_CALL( SCIPsetSetBarrierconvtol(scip->set, barrierconvtol) );

   return SCIP_OKAY;
}

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
SCIP_Real SCIPchgRelaxfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             relaxfeastol        /**< new primal feasibility tolerance of relaxations */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPchgRelaxfeastol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPsetSetRelaxfeastol(scip->set, relaxfeastol);
}

/** marks that some limit parameter was changed */
void SCIPmarkLimitChanged(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   /* change the settings */
   SCIPsetSetLimitChanged(scip->set);
}

/** outputs a real number, or "+infinity", or "-infinity" to a file */
void SCIPprintReal(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             val,                /**< value to print */
   int                   width,              /**< width of the field */
   int                   precision           /**< number of significant digits printed */
   )
{
   char s[SCIP_MAXSTRLEN];
   char strformat[SCIP_MAXSTRLEN];

   assert(scip != NULL);

   if( SCIPsetIsInfinity(scip->set, val) )
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "+infinity");
   else if( SCIPsetIsInfinity(scip->set, -val) )
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "-infinity");
   else
   {
      (void) SCIPsnprintf(strformat, SCIP_MAXSTRLEN, "%%.%dg", precision);
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, (const char*)strformat, val);
   }
   (void) SCIPsnprintf(strformat, SCIP_MAXSTRLEN, "%%%ds", width);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, (const char*)strformat, s);
}

/** parse a real value that was written with SCIPprintReal() */
SCIP_Bool SCIPparseReal(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed, otherwise @p str */
   )
{
   char* localstr;

   assert(scip != NULL);
   assert(str != NULL);
   assert(value != NULL);
   assert(endptr != NULL);

   localstr = (char*)str;

   /* ignore white space */
   while(isspace((unsigned char)*localstr))
      ++localstr;

   /* test for a special infinity first */
   if( strncmp(localstr, "+infinity", 9) == 0 )
   {
      *value = SCIPinfinity(scip);
      *endptr = (char*)(localstr + 9);
      return TRUE;
   }
   else if( strncmp(localstr, "-infinity", 9) == 0 )
   {
      *value = -SCIPinfinity(scip);
      *endptr = (char*)(localstr + 9);
      return TRUE;
   }
   else
   {
      /* parse a finite value */
      return SCIPstrToRealValue(str, value, endptr);
   }
}

/** checks, if values are in range of epsilon */
SCIP_Bool SCIPisEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) lower than val2 */
SCIP_Bool SCIPisLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) greater than val2 */
SCIP_Bool SCIPisLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) greater than val2 */
SCIP_Bool SCIPisGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) lower than val2 */
SCIP_Bool SCIPisGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsGE(scip->set, val1, val2);
}

/** returns value treated as infinity */
SCIP_Real SCIPinfinity(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetInfinity(scip->set);
}

/** checks, if value is (positive) infinite */
SCIP_Bool SCIPisInfinity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to be compared against infinity */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsInfinity(scip->set, val);
}

/** checks, if value is huge and should be handled separately (e.g., in activity computation) */
SCIP_Bool SCIPisHugeValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to be checked whether it is huge */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsHugeValue(scip->set, val);
}

/** returns the minimum value that is regarded as huge and should be handled separately (e.g., in activity
 *  computation)
 */
SCIP_Real SCIPgetHugeValue(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetGetHugeValue(scip->set);
}

/** checks, if value is in range epsilon of 0.0 */
SCIP_Bool SCIPisZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsZero(scip->set, val);
}

/** checks, if value is greater than epsilon */
SCIP_Bool SCIPisPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsPositive(scip->set, val);
}

/** checks, if value is lower than -epsilon */
SCIP_Bool SCIPisNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsNegative(scip->set, val);
}

/** checks, if value is integral within epsilon */
SCIP_Bool SCIPisIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsIntegral(scip->set, val);
}

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
SCIP_Bool SCIPisScalingIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val,                /**< unscaled value to check for scaled integrality */
   SCIP_Real             scalar              /**< value to scale val with for checking for integrality */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsScalingIntegral(scip->set, val, scalar);
}

/** checks, if given fractional part is smaller than epsilon */
SCIP_Bool SCIPisFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFracIntegral(scip->set, val);
}

/** rounds value + epsilon down to the next integer */
SCIP_Real SCIPfloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFloor(scip->set, val);
}

/** rounds value - epsilon up to the next integer */
SCIP_Real SCIPceil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetCeil(scip->set, val);
}

/** rounds value to the nearest integer with epsilon tolerance */
SCIP_Real SCIPround(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetRound(scip->set, val);
}

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
SCIP_Real SCIPfrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to return fractional part for */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFrac(scip->set, val);
}

/** checks, if values are in range of sumepsilon */
SCIP_Bool SCIPisSumEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPisSumLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPisSumLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPisSumGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPisSumGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumGE(scip->set, val1, val2);
}

/** checks, if value is in range sumepsilon of 0.0 */
SCIP_Bool SCIPisSumZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumZero(scip->set, val);
}

/** checks, if value is greater than sumepsilon */
SCIP_Bool SCIPisSumPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumPositive(scip->set, val);
}

/** checks, if value is lower than -sumepsilon */
SCIP_Bool SCIPisSumNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumNegative(scip->set, val);
}

/** checks, if relative difference of values is in range of feasibility tolerance */
SCIP_Bool SCIPisFeasEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasEQ(scip->set, val1, val2);
}

/** checks, if relative difference val1 and val2 is lower than feasibility tolerance */
SCIP_Bool SCIPisFeasLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than feasibility tolerance */
SCIP_Bool SCIPisFeasLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than feastol */
SCIP_Bool SCIPisFeasGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
SCIP_Bool SCIPisFeasGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasGE(scip->set, val1, val2);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
SCIP_Bool SCIPisFeasZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasZero(scip->set, val);
}

/** checks, if value is greater than feasibility tolerance */
SCIP_Bool SCIPisFeasPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasPositive(scip->set, val);
}

/** checks, if value is lower than -feasibility tolerance */
SCIP_Bool SCIPisFeasNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasNegative(scip->set, val);
}

/** checks, if value is integral within the LP feasibility bounds */
SCIP_Bool SCIPisFeasIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasIntegral(scip->set, val);
}

/** checks, if given fractional part is smaller than feastol */
SCIP_Bool SCIPisFeasFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasFracIntegral(scip->set, val);
}

/** rounds value + feasibility tolerance down to the next integer */
SCIP_Real SCIPfeasFloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasFloor(scip->set, val);
}

/** rounds value - feasibility tolerance up to the next integer */
SCIP_Real SCIPfeasCeil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasCeil(scip->set, val);
}

/** rounds value to the nearest integer in feasibility tolerance */
SCIP_Real SCIPfeasRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasRound(scip->set, val);
}

/** returns fractional part of value, i.e. x - floor(x) in feasibility tolerance */
SCIP_Real SCIPfeasFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasFrac(scip->set, val);
}

/** checks, if relative difference of values is in range of dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasEQ(scip->set, val1, val2);
}

/** checks, if relative difference val1 and val2 is lower than dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasGE(scip->set, val1, val2);
}

/** checks, if value is in range dual feasibility tolerance of 0.0 */
SCIP_Bool SCIPisDualfeasZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasZero(scip->set, val);
}

/** checks, if value is greater than dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasPositive(scip->set, val);
}

/** checks, if value is lower than -dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasNegative(scip->set, val);
}

/** checks, if value is integral within the LP dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasIntegral(scip->set, val);
}

/** checks, if given fractional part is smaller than dual feasibility tolerance */
SCIP_Bool SCIPisDualfeasFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsDualfeasFracIntegral(scip->set, val);
}

/** rounds value + dual feasibility tolerance down to the next integer */
SCIP_Real SCIPdualfeasFloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetDualfeasFloor(scip->set, val);
}

/** rounds value - dual feasibility tolerance up to the next integer */
SCIP_Real SCIPdualfeasCeil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetDualfeasCeil(scip->set, val);
}

/** rounds value to the nearest integer in dual feasibility tolerance */
SCIP_Real SCIPdualfeasRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetDualfeasRound(scip->set, val);
}

/** returns fractional part of value, i.e. x - floor(x) in dual feasibility tolerance */
SCIP_Real SCIPdualfeasFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetDualfeasFrac(scip->set, val);
}

/** checks, if the given new lower bound is at least min(oldub - oldlb, |oldlb|) times the bound
 *  strengthening epsilon better than the old one
 */
SCIP_Bool SCIPisLbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   assert(scip != NULL);

   return SCIPsetIsLbBetter(scip->set, newlb, oldlb, oldub);
}

/** checks, if the given new upper bound is at least min(oldub - oldlb, |oldub|) times the bound
 *  strengthening epsilon better than the old one
 */
SCIP_Bool SCIPisUbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   assert(scip != NULL);

   return SCIPsetIsUbBetter(scip->set, newub, oldlb, oldub);
}

/** checks, if relative difference of values is in range of epsilon */
SCIP_Bool SCIPisRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelEQ(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is lower than epsilon */
SCIP_Bool SCIPisRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
SCIP_Bool SCIPisRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than epsilon */
SCIP_Bool SCIPisRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
SCIP_Bool SCIPisRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsRelGE(scip->set, val1, val2);
}

/** checks, if relative difference of values is in range of sumepsilon */
SCIP_Bool SCIPisSumRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelEQ(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
SCIP_Bool SCIPisSumRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
SCIP_Bool SCIPisSumRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
SCIP_Bool SCIPisSumRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
SCIP_Bool SCIPisSumRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumRelGE(scip->set, val1, val2);
}

/** converts the given real number representing an integer to an int; in optimized mode the function gets inlined for
 *  performance; in debug mode we check some additional conditions
 */
int SCIPconvertRealToInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             real                /**< double bound to convert */
   )
{
   assert(SCIPisFeasIntegral(scip, real));
   assert(SCIPisFeasEQ(scip, real, (SCIP_Real)(int)(real < 0 ? real - 0.5 : real + 0.5)));
   assert(real < INT_MAX);
   assert(real > INT_MIN);

   return (int)(real < 0 ? (real - 0.5) : (real + 0.5));
}

/** converts the given real number representing an integer to a long integer; in optimized mode the function gets inlined for
 *  performance; in debug mode we check some additional conditions
 */
SCIP_Longint SCIPconvertRealToLongint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             real                /**< double bound to convert */
   )
{
   assert(SCIPisFeasIntegral(scip, real));
   assert(SCIPisFeasEQ(scip, real, (SCIP_Real)(SCIP_Longint)(real < 0 ? real - 0.5 : real + 0.5)));
   assert(real < SCIP_LONGINT_MAX);
   assert(real > SCIP_LONGINT_MIN);

   return (SCIP_Longint)(real < 0 ? (real - 0.5) : (real + 0.5));
}

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
SCIP_Bool SCIPisUpdateUnreliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newvalue,           /**< new value after update */
   SCIP_Real             oldvalue            /**< old value, i.e., last reliable value */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisUpdateUnreliable", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPsetIsUpdateUnreliable(scip->set, newvalue, oldvalue);
}
