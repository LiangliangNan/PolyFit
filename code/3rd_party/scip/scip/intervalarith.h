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

/**@file   intervalarith.h
 * @ingroup PUBLICCOREAPI
 * @brief  interval arithmetics for provable bounds
 * @author Tobias Achterberg
 * @author Stefan Vigerske
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_INTERVALARITH_H__
#define __SCIP_INTERVALARITH_H__


#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@defgroup PublicIntervalArithMethods Interval Arithmetics
 * @ingroup MiscellaneousMethods
 * @brief methods for interval arithmetics
 *
 * @{
 */

/** interval given by infimum and supremum */
struct SCIP_Interval
{
   SCIP_Real             inf;                /**< infimum (lower bound) of interval */
   SCIP_Real             sup;                /**< supremum (upper bound) of interval */
};
typedef struct SCIP_Interval SCIP_INTERVAL;

/** rounding mode of floating point operations (upwards, downwards, nearest, ...)
 *
 * exact values depend on machine and compiler
 */
typedef int SCIP_ROUNDMODE;

/*
 * Interval arithmetic operations
 */

/** returns whether rounding mode control is available */
SCIP_EXPORT
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   );

/** sets rounding mode of floating point operations */
SCIP_EXPORT
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE        roundmode           /**< rounding mode to activate */
   );

/** gets current rounding mode of floating point operations */
SCIP_EXPORT
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   );

/** sets rounding mode of floating point operations to downwards rounding */
SCIP_EXPORT
void SCIPintervalSetRoundingModeDownwards(
   void
   );

/** sets rounding mode of floating point operations to upwards rounding */
SCIP_EXPORT
void SCIPintervalSetRoundingModeUpwards(
   void
   );

/** sets rounding mode of floating point operations to nearest rounding */
SCIP_EXPORT
void SCIPintervalSetRoundingModeToNearest(
   void
   );

/** sets rounding mode of floating point operations to towards zero rounding */
SCIP_EXPORT
void SCIPintervalSetRoundingModeTowardsZero(
   void
   );

/** negates a number in a way that the compiler does not optimize it away */
SCIP_EXPORT
SCIP_Real SCIPintervalNegateReal(
   SCIP_Real             x                   /**< number to negate */
   );

/** returns infimum of interval */
SCIP_EXPORT
SCIP_Real SCIPintervalGetInf(
   SCIP_INTERVAL         interval            /**< interval */
   );

/** returns supremum of interval */
SCIP_EXPORT
SCIP_Real SCIPintervalGetSup(
   SCIP_INTERVAL         interval            /**< interval */
   );

/** stores given value as interval */
SCIP_EXPORT
void SCIPintervalSet(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             value               /**< value to store */
   );

/** stores given infimum and supremum as interval */
SCIP_EXPORT
void SCIPintervalSetBounds(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             inf,                /**< value to store as infimum */
   SCIP_Real             sup                 /**< value to store as supremum */
   );

/** sets interval to empty interval, which will be [1.0, -1.0] */
SCIP_EXPORT
void SCIPintervalSetEmpty(
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   );

/** indicates whether interval is empty, i.e., whether inf > sup */
SCIP_EXPORT
SCIP_Bool SCIPintervalIsEmpty(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** sets interval to entire [-infinity, +infinity] */
SCIP_EXPORT
void SCIPintervalSetEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   );

/** indicates whether interval is entire, i.e., whether inf &le; -infinity and sup &ge; infinity */
SCIP_EXPORT
SCIP_Bool SCIPintervalIsEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** indicates whether interval is positive infinity, i.e., [infinity, infinity] */
SCIP_EXPORT
SCIP_Bool SCIPintervalIsPositiveInfinity(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** indicates whether interval is negative infinity, i.e., [-infinity, -infinity] */
SCIP_EXPORT
SCIP_Bool SCIPintervalIsNegativeInfinity(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

#ifdef NDEBUG

/* In optimized mode, some function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 * With SCIPintervalSetBounds we need to be a bit careful, since i and s could use resultant->inf and resultant->sup,
 * e.g., SCIPintervalSetBounds(&resultant, -resultant->sup, -resultant->inf).
 * So we need to make sure that we first evaluate both terms before setting resultant. 
 */

#define SCIPintervalGetInf(interval)               (interval).inf
#define SCIPintervalGetSup(interval)               (interval).sup
#define SCIPintervalSet(resultant, value)          do { (resultant)->inf = (value); (resultant)->sup = (resultant)->inf; } while( FALSE )
#define SCIPintervalSetBounds(resultant, i, s)     do { SCIP_Real scipintervaltemp; scipintervaltemp = (s); (resultant)->inf = (i); (resultant)->sup = scipintervaltemp; } while( FALSE )
#define SCIPintervalSetEmpty(resultant)            do { (resultant)->inf = 1.0; (resultant)->sup = -1.0; } while( FALSE )
#define SCIPintervalSetEntire(infinity, resultant) do { (resultant)->inf = -(infinity); (resultant)->sup =  (infinity); } while( FALSE )
#define SCIPintervalIsEmpty(infinity, operand)     ( (operand).inf > -(infinity) && (operand).sup < (infinity) && (operand).sup < (operand).inf )
#define SCIPintervalIsEntire(infinity, operand)    ( (operand).inf <= -(infinity) && (operand).sup >= (infinity) )
#define SCIPintervalIsPositiveInfinity(infinity, operand) ( (operand).inf >=  (infinity) && (operand).sup >= (operand).inf )
#define SCIPintervalIsNegativeInfinity(infinity, operand) ( (operand).sup <= -(infinity) && (operand).sup >= (operand).inf )

#endif

/** indicates whether operand1 is contained in operand2 */
SCIP_EXPORT
SCIP_Bool SCIPintervalIsSubsetEQ(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** indicates whether operand1 and operand2 are disjoint */
SCIP_EXPORT
SCIP_Bool SCIPintervalAreDisjoint(
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** indicates whether operand1 and operand2 are disjoint with epsilon tolerance
 *
 * Returns whether minimal (relative) distance of intervals is larger than epsilon.
 * Same as `SCIPintervalIsEmpty(SCIPintervalIntersectEps(operand1, operand2))`.
 */
SCIP_EXPORT
SCIP_Bool SCIPintervalAreDisjointEps(
   SCIP_Real             eps,                /**< epsilon */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** intersection of two intervals */
SCIP_EXPORT
void SCIPintervalIntersect(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** intersection of two intervals with epsilon tolerance
 *
 * If intersection of operand1 and operand2 is empty, but minimal (relative) distance of intervals
 * is at most epsilon, then set resultant to singleton containing the point in operand1
 * that is closest to operand2, i.e.,
 * - `resultant = { operand1.sup }`, if `operand1.sup` < `operand2.inf` and `reldiff(operand2.inf,operand1.sup)` &le; eps
 * - `resultant = { operand1.inf }`, if `operand1.inf` > `operand2.sup` and `reldiff(operand1.inf,operand2.sup)` &le; eps
 * - `resultant` = intersection of `operand1` and `operand2`, otherwise
 */
SCIP_EXPORT
void SCIPintervalIntersectEps(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             eps,                /**< epsilon */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** interval enclosure of the union of two intervals */
SCIP_EXPORT
void SCIPintervalUnify(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and operand2 and stores infimum of result in infimum of resultant */
SCIP_EXPORT
void SCIPintervalAddInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and operand2 and stores supremum of result in supremum of resultant */
SCIP_EXPORT
void SCIPintervalAddSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and operand2 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalAdd(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and scalar operand2 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalAddScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** adds vector operand1 and vector operand2 and stores result in vector resultant */
SCIP_EXPORT
void SCIPintervalAddVectors(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< array of resultant intervals of operation */
   int                   length,             /**< length of arrays */
   SCIP_INTERVAL*        operand1,           /**< array of first operands of operation */
   SCIP_INTERVAL*        operand2            /**< array of second operands of operation */
   );

/** subtracts operand2 from operand1 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalSub(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** subtracts scalar operand2 from operand1 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalSubScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** multiplies operand1 with operand2 and stores infimum of result in infimum of resultant */
SCIP_EXPORT
void SCIPintervalMulInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   );

/** multiplies operand1 with operand2 and stores supremum of result in supremum of resultant */
SCIP_EXPORT
void SCIPintervalMulSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   );

/** multiplies operand1 with operand2 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalMul(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** multiplies operand1 with scalar operand2 and stores infimum of result in infimum of resultant */
SCIP_EXPORT
void SCIPintervalMulScalarInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   );

/** multiplies operand1 with scalar operand2 and stores supremum of result in supremum of resultant */
SCIP_EXPORT
void SCIPintervalMulScalarSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   );

/** multiplies operand1 with scalar operand2 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalMulScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** divides operand1 by operand2 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalDiv(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** divides operand1 by scalar operand2 and stores result in resultant */
SCIP_EXPORT
void SCIPintervalDivScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** computes the scalar product of two vectors of intervals and stores result in resultant */
SCIP_EXPORT
void SCIPintervalScalprod(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first  vector as array of intervals */
   SCIP_INTERVAL*        operand2            /**< second vector as array of intervals */
   );

/** computes the scalar product of a vector of intervals and a vector of scalars and stores infimum of result in infimum of resultant */
SCIP_EXPORT
void SCIPintervalScalprodScalarsInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   );

/** computes the scalar product of a vector of intervals and a vector of scalars and stores supremum of result in supremum of resultant */
SCIP_EXPORT
void SCIPintervalScalprodScalarsSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   );

/** computes the scalar product of a vector of intervals and a vector of scalars and stores result in resultant */
SCIP_EXPORT
void SCIPintervalScalprodScalars(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   );

/** squares operand and stores result in resultant */
SCIP_EXPORT
void SCIPintervalSquare(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores (positive part of) square root of operand in resultant
 * @attention we assume a correctly rounded sqrt(double) function when rounding is to nearest
 */
SCIP_EXPORT
void SCIPintervalSquareRoot(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores operand1 to the power of operand2 in resultant
 * 
 * uses SCIPintervalPowerScalar if operand2 is a scalar, otherwise computes exp(op2*log(op1))
 */
SCIP_EXPORT
void SCIPintervalPower(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** stores operand1 to the power of the scalar operand2 in resultant
 * @attention we assume a correctly rounded pow(double) function when rounding is to nearest
 */
SCIP_EXPORT
void SCIPintervalPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** stores bounds on the power of a scalar operand1 to a scalar operand2 in resultant
 *
 * Both operands need to be finite numbers.
 * Needs to have operand1 &ge; 0 or operand2 integer and needs to have operand2 &ge; 0 if operand1 = 0.
 * @attention we assume a correctly rounded pow(double) function when rounding is to nearest
 */
SCIP_EXPORT
void SCIPintervalPowerScalarScalar(
   SCIP_INTERVAL*        resultant,          /**< resultant of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** computes lower bound on power of a scalar operand1 to an integer operand2
 *
 * Both operands need to be finite numbers.
 * Needs to have operand1 &ge; 0 and need to have operand2 &ge; 0 if operand1 = 0.
 */
SCIP_EXPORT
SCIP_Real SCIPintervalPowerScalarIntegerInf(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   );

/** computes upper bound on power of a scalar operand1 to an integer operand2
 *
 * Both operands need to be finite numbers.
 * Needs to have operand1 &ge; 0 and needs to have operand2 &ge; 0 if operand1 = 0.
 */
SCIP_EXPORT
SCIP_Real SCIPintervalPowerScalarIntegerSup(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   );

/** computes bounds on power of a scalar operand1 to an integer operand2
 *
 * Both operands need to be finite numbers.
 * Needs to have operand1 &ge; 0 and needs to have operand2 &ge; 0 if operand1 = 0.
 */
SCIP_EXPORT
void SCIPintervalPowerScalarInteger(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   );

/** given an interval for the image of a power operation, computes an interval for the origin
 *
 * That is, for \f$y = x^p\f$ with the exponent \f$p\f$ a given scalar and \f$y\f$ = `image` a given interval,
 * computes \f$x \subseteq \text{basedomain}\f$ such that \f$y \in x^p\f$ and such that for all \f$z \in \text{basedomain} \setminus x: z^p \not \in y\f$.
 */
SCIP_EXPORT
void SCIPintervalPowerScalarInverse(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         basedomain,         /**< domain of base */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_INTERVAL         image               /**< interval image of power */
   );

/** stores operand1 to the signed power of the scalar positive operand2 in resultant
 * 
 * The signed power of x w.r.t. an exponent n &ge; 0 is given as \f$\mathrm{sign}(x) |x|^n\f$.
 *
 * @attention we assume correctly rounded sqrt(double) and pow(double) functions when rounding is to nearest
 */
SCIP_EXPORT
void SCIPintervalSignPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** computes the reciprocal of an interval
 */
SCIP_EXPORT
void SCIPintervalReciprocal(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores exponential of operand in resultant
 * @attention we assume a correctly rounded exp(double) function when rounding is to nearest
 */
SCIP_EXPORT
void SCIPintervalExp(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores natural logarithm of operand in resultant
 * @attention we assume a correctly rounded log(double) function when rounding is to nearest
 */
SCIP_EXPORT
void SCIPintervalLog(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores minimum of operands in resultant */
SCIP_EXPORT
void SCIPintervalMin(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** stores maximum of operands in resultant */
SCIP_EXPORT
void SCIPintervalMax(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** stores absolute value of operand in resultant */
SCIP_EXPORT
void SCIPintervalAbs(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores sine value of operand in resultant */
SCIP_EXPORT
void SCIPintervalSin(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores cosine value of operand in resultant */
SCIP_EXPORT
void SCIPintervalCos(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores sign of operand in resultant */
SCIP_EXPORT
void SCIPintervalSign(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores entropy of operand in resultant */
SCIP_EXPORT
void SCIPintervalEntropy(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** computes exact upper bound on \f$ a x^2 + b x \f$ for x in [xlb, xub], b an interval, and a scalar
 * 
 * Uses Algorithm 2.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008).
 */
SCIP_EXPORT
SCIP_Real SCIPintervalQuadUpperBound(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_Real             a,                  /**< coefficient of x^2 */
   SCIP_INTERVAL         b_,                 /**< coefficient of x */
   SCIP_INTERVAL         x                   /**< range of x */
   );

/** stores range of quadratic term in resultant
 * 
 * given scalar a and intervals b and x, computes interval for \f$ a x^2 + b x \f$ */
SCIP_EXPORT
void SCIPintervalQuad(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         xrng                /**< range of x */
   );


/** computes interval with positive solutions of a quadratic equation with interval coefficients
 * 
 * Given intervals a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \in c\f$ within xbnds.
 */
SCIP_EXPORT
void SCIPintervalSolveUnivariateQuadExpressionPositive(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs,                /**< right hand side of equation */
   SCIP_INTERVAL         xbnds               /**< bounds on x */
   );

/** computes interval with negative solutions of a quadratic equation with interval coefficients
 *
 * Given intervals a, b, and c, this function computes an interval that contains all negative solutions of \f$ a x^2 + b x \in c\f$ within xbnds.
 */
SCIP_EXPORT
void SCIPintervalSolveUnivariateQuadExpressionNegative(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs,                /**< right hand side of equation */
   SCIP_INTERVAL         xbnds               /**< bounds on x */
   );

/** computes positive solutions of a quadratic equation with scalar coefficients
 * 
 * Givens scalar a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \geq c\f$ within xbnds.
 * Implements Algorithm 3.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008).
 */
SCIP_EXPORT
void SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_Real             lincoeff,           /**< coefficient of x */
   SCIP_Real             rhs,                /**< right hand side of equation */
   SCIP_INTERVAL         xbnds               /**< bounds on x */
   );

/** solves a quadratic equation with interval coefficients
 *
 * Given intervals a, b and c, this function computes an interval that contains all solutions of \f$ a x^2 + b x \in c\f$ within xbnds.
 */
SCIP_EXPORT
void SCIPintervalSolveUnivariateQuadExpression(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs,                /**< right hand side of equation */
   SCIP_INTERVAL         xbnds               /**< bounds on x */
   );

/** stores range of bivariate quadratic term in resultant
 *
 * Given scalars \f$a_x\f$, \f$a_y\f$, \f$a_{xy}\f$, \f$b_x\f$, and \f$b_y\f$ and intervals for \f$x\f$ and \f$y\f$,
 * computes interval for \f$ a_x x^2 + a_y y^2 + a_{xy} x y + b_x x + b_y y \f$.
 *
 * \attention The operations are not applied rounding-safe here!
 */
SCIP_EXPORT
void SCIPintervalQuadBivar(
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_INTERVAL*        resultant,          /**< buffer where to store result of operation */
   SCIP_Real             ax,                 /**< square coefficient of x */
   SCIP_Real             ay,                 /**< square coefficient of y */
   SCIP_Real             axy,                /**< bilinear coefficients */
   SCIP_Real             bx,                 /**< linear coefficient of x */
   SCIP_Real             by,                 /**< linear coefficient of y */
   SCIP_INTERVAL         xbnds,              /**< bounds on x */
   SCIP_INTERVAL         ybnds               /**< bounds on y */
   );

/** solves a bivariate quadratic equation for the first variable
 *
 * Given scalars \f$a_x\f$, \f$a_y\f$, \f$a_{xy}\f$, \f$b_x\f$ and \f$b_y\f$, and intervals for \f$x\f$, \f$y\f$, and rhs,
 * computes \f$ \{ x \in \mathbf{x} : \exists y \in \mathbf{y} : a_x x^2 + a_y y^2 + a_{xy} x y + b_x x + b_y y \in \mathbf{\mbox{rhs}} \} \f$.
 *
 * \attention the operations are not applied rounding-safe here
 */
SCIP_EXPORT
void SCIPintervalSolveBivariateQuadExpressionAllScalar(
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   SCIP_INTERVAL*        resultant,          /**< buffer where to store result of operation */
   SCIP_Real             ax,                 /**< square coefficient of x */
   SCIP_Real             ay,                 /**< square coefficient of y */
   SCIP_Real             axy,                /**< bilinear coefficients */
   SCIP_Real             bx,                 /**< linear coefficient of x */
   SCIP_Real             by,                 /**< linear coefficient of y */
   SCIP_INTERVAL         rhs,                /**< right-hand-side of equation */
   SCIP_INTERVAL         xbnds,              /**< bounds on x */
   SCIP_INTERVAL         ybnds               /**< bounds on y */
   );

/** propagates a weighted sum of intervals in a given interval
 *
 * Given \f$\text{constant} + \sum_i \text{weights}_i \text{operands}_i \in \text{rhs}\f$,
 * computes possibly tighter interval for each term.
 *
 * @attention Valid values are returned in resultants only if any tightening has been found and no empty interval, that is, function returns with non-zero and `*infeasible` = FALSE.
 *
 * @return Number of terms for which resulting interval is smaller than operand interval.
 */
SCIP_EXPORT
int SCIPintervalPropagateWeightedSum(
   SCIP_Real             infinity,           /**< value for infinity in interval arithmetics */
   int                   noperands,          /**< number of operands (intervals) to propagate */
   SCIP_INTERVAL*        operands,           /**< intervals to propagate */
   SCIP_Real*            weights,            /**< weights of intervals in sum */
   SCIP_Real             constant,           /**< constant in sum */
   SCIP_INTERVAL         rhs,                /**< right-hand-side interval */
   SCIP_INTERVAL*        resultants,         /**< array to store propagated intervals, if any reduction is found at all (check return code and *infeasible) */
   SCIP_Bool*            infeasible          /**< buffer to store if propagation produced empty interval */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
