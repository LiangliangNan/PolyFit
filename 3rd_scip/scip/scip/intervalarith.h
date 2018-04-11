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

/**@file   scip/intervalarith.h
 * @ingroup INTERNALAPI
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

/** interval given by infimum and supremum */
struct SCIP_Interval
{
   SCIP_Real             inf;                /**< infimum (lower bound) of interval */
   SCIP_Real             sup;                /**< supremum (upper bound) of interval */
};
typedef struct SCIP_Interval SCIP_INTERVAL;

/** rounding mode of floating point operations (upwards, downwards, nearest, ...)
 * exact values depend on machine and compiler, so we define a corresponding enum in the header file */
typedef int SCIP_ROUNDMODE;

/*
 * Interval arithmetic operations
 */

/** returns whether rounding mode control is available */
extern
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   );

/** sets rounding mode of floating point operations */
extern
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE        roundmode           /**< rounding mode to activate */
   );

/** gets current rounding mode of floating point operations */
extern
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   );

/** sets rounding mode of floating point operations to downwards rounding */
extern
void SCIPintervalSetRoundingModeDownwards(
   void
   );

/** sets rounding mode of floating point operations to upwards rounding */
extern
void SCIPintervalSetRoundingModeUpwards(
   void
   );

/** sets rounding mode of floating point operations to nearest rounding */
extern
void SCIPintervalSetRoundingModeToNearest(
   void
   );

/** sets rounding mode of floating point operations to towards zero rounding */
extern
void SCIPintervalSetRoundingModeTowardsZero(
   void
   );

/** negates a number in a way that the compiler does not optimize it away */
extern
SCIP_Real SCIPintervalNegateReal(
   SCIP_Real             x                   /**< number to negate */
   );

/** returns infimum of interval */
extern
SCIP_Real SCIPintervalGetInf(
   SCIP_INTERVAL         interval            /**< interval */
   );

/** returns supremum of interval */
extern
SCIP_Real SCIPintervalGetSup(
   SCIP_INTERVAL         interval            /**< interval */
   );

/** stores given value as interval */
extern
void SCIPintervalSet(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             value               /**< value to store */
   );

/** stores given infimum and supremum as interval */
extern
void SCIPintervalSetBounds(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             inf,                /**< value to store as infimum */
   SCIP_Real             sup                 /**< value to store as supremum */
   );

/** sets interval to empty interval, which will be [1.0, -1.0] */
extern
void SCIPintervalSetEmpty(
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   );

/** indicates whether interval is empty, i.e., whether inf > sup */
extern
SCIP_Bool SCIPintervalIsEmpty(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** sets interval to entire [-infinity, +infinity] */
extern
void SCIPintervalSetEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   );

/** indicates whether interval is entire, i.e., whether inf <= -infinity and sup >= infinity */
extern
SCIP_Bool SCIPintervalIsEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** indicates whether interval is positive infinity, i.e., [infinity, infinity] */
extern
SCIP_Bool SCIPintervalIsPositiveInfinity(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** indicates whether interval is negative infinity, i.e., [-infinity, -infinity] */
extern
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
extern
SCIP_Bool SCIPintervalIsSubsetEQ(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** indicates whether operand1 and operand2 are disjoint */
extern
SCIP_Bool SCIPintervalAreDisjoint(
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** intersection of two intervals */
extern
void SCIPintervalIntersect(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** interval enclosure of the union of two intervals */
extern
void SCIPintervalUnify(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and operand2 and stores infimum of result in infimum of resultant */
extern
void SCIPintervalAddInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and operand2 and stores supremum of result in supremum of resultant */
extern
void SCIPintervalAddSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and operand2 and stores result in resultant */
extern
void SCIPintervalAdd(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** adds operand1 and scalar operand2 and stores result in resultant */
extern
void SCIPintervalAddScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** adds vector operand1 and vector operand2 and stores result in vector resultant */
extern
void SCIPintervalAddVectors(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< array of resultant intervals of operation */
   int                   length,             /**< length of arrays */
   SCIP_INTERVAL*        operand1,           /**< array of first operands of operation */
   SCIP_INTERVAL*        operand2            /**< array of second operands of operation */
   );

/** subtracts operand2 from operand1 and stores result in resultant */
extern
void SCIPintervalSub(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** subtracts scalar operand2 from operand1 and stores result in resultant */
extern
void SCIPintervalSubScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** multiplies operand1 with operand2 and stores infimum of result in infimum of resultant */
extern
void SCIPintervalMulInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   );

/** multiplies operand1 with operand2 and stores supremum of result in supremum of resultant */
extern
void SCIPintervalMulSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   );

/** multiplies operand1 with operand2 and stores result in resultant */
extern
void SCIPintervalMul(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** multiplies operand1 with scalar operand2 and stores infimum of result in infimum of resultant */
extern
void SCIPintervalMulScalarInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   );

/** multiplies operand1 with scalar operand2 and stores supremum of result in supremum of resultant */
extern
void SCIPintervalMulScalarSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   );

/** multiplies operand1 with scalar operand2 and stores result in resultant */
extern
void SCIPintervalMulScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** divides operand1 by operand2 and stores result in resultant */
extern
void SCIPintervalDiv(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** divides operand1 by scalar operand2 and stores result in resultant */
extern
void SCIPintervalDivScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** computes the scalar product of two vectors of intervals and stores result in resultant */
extern
void SCIPintervalScalprod(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first  vector as array of intervals */
   SCIP_INTERVAL*        operand2            /**< second vector as array of intervals */
   );

/** computes the scalar product of a vector of intervals and a vector of scalars and stores infimum of result in infimum 
 *  of resultant 
 */
extern
void SCIPintervalScalprodScalarsInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   );

/** computes the scalar product of a vector of intervals and a vector of scalars and stores supremum of result in supremum
 *  of resultant 
 */
extern
void SCIPintervalScalprodScalarsSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   );

/** computes the scalar product of a vector of intervals and a vector of scalars and stores result in resultant */
extern
void SCIPintervalScalprodScalars(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   );

/** squares operand and stores result in resultant */
extern
void SCIPintervalSquare(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores (positive part of) square root of operand in resultant
 * @attention we assume a correctly rounded sqrt(double) function when rounding is to nearest
 */
extern
void SCIPintervalSquareRoot(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores operand1 to the power of operand2 in resultant
 * 
 * uses SCIPintervalPowerScalar if operand2 is a scalar, otherwise computes exp(op2*log(op1))
 */
extern
void SCIPintervalPower(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** stores operand1 to the power of the scalar operand2 in resultant */
extern
void SCIPintervalPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** stores bounds on the power of a scalar operand1 to a scalar operand2 in resultant
 * both operands need to be finite numbers
 * need to have operand1 >= 0 or operand2 integer and need to have operand2 >= 0 if operand1 == 0
 * @attention we assume a correctly rounded pow(double) function when rounding is to nearest
 */
extern
void SCIPintervalPowerScalarScalar(
   SCIP_INTERVAL*        resultant,          /**< resultant of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** computes lower bound on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
extern
SCIP_Real SCIPintervalPowerScalarIntegerInf(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   );

/** computes upper bound on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
extern
SCIP_Real SCIPintervalPowerScalarIntegerSup(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   );

/** computes bounds on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
extern
void SCIPintervalPowerScalarInteger(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   );

/** given an interval for the image of a power operation, computes an interval for the origin
 * that is, for y = x^p with p = exponent a given scalar and y = image a given interval,
 * computes a subinterval x of basedomain such that y in x^p and such that for all z in basedomain less x, z^p not in y
 */
extern
void SCIPintervalPowerScalarInverse(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         basedomain,         /**< domain of base */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_INTERVAL         image               /**< interval image of power */
   );

/** stores operand1 to the signed power of the scalar positive operand2 in resultant 
 * 
 * the signed power of x w.r.t. an exponent n >= 0 is given as sign(x) * abs(x)^n
 *
 * @attention we assume correctly rounded sqrt(double) and pow(double) functions when rounding is to nearest
 */
extern
void SCIPintervalSignPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   );

/** computes the reciprocal of an interval
 */
extern
void SCIPintervalReciprocal(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores exponential of operand in resultant
 * @attention we assume a correctly rounded exp(double) function when rounding is to nearest
 */
extern
void SCIPintervalExp(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores natural logarithm of operand in resultant
 * @attention we assume a correctly rounded log(double) function when rounding is to nearest
 */
extern
void SCIPintervalLog(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores minimum of operands in resultant */
extern
void SCIPintervalMin(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** stores maximum of operands in resultant */
extern
void SCIPintervalMax(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   );

/** stores absolute value of operand in resultant */
extern
void SCIPintervalAbs(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores sine value of operand in resultant
 * NOTE: the operations are not applied rounding-safe here
 */
extern
void SCIPintervalSin(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores cosine value of operand in resultant
 * NOTE: the operations are not applied rounding-safe here
 */
extern
void SCIPintervalCos(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** stores sign of operand in resultant */
extern
void SCIPintervalSign(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   );

/** computes exact upper bound on \f$ a x^2 + b x \f$ for x in [xlb, xub], b an interval, and a scalar
 * 
 * Uses Algorithm 2.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008) */
extern
SCIP_Real SCIPintervalQuadUpperBound(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_Real             a,                  /**< coefficient of x^2 */
   SCIP_INTERVAL         b_,                 /**< coefficient of x */
   SCIP_INTERVAL         x                   /**< range of x */
   );

/** stores range of quadratic term in resultant
 * 
 * given scalar a and intervals b and x, computes interval for \f$ a x^2 + b x \f$ */
extern
void SCIPintervalQuad(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         xrng                /**< range of x */
   );


/** computes interval with positive solutions of a quadratic equation with interval coefficients
 * 
 * Given intervals a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \in c\f$.
 */
extern
void SCIPintervalSolveUnivariateQuadExpressionPositive(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs                 /**< right hand side of equation */
   );

/** computes positive solutions of a quadratic equation with scalar coefficients
 * 
 * Given scalar a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \geq c\f$.
 * Implements Algorithm 3.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008).
 */
extern
void SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_Real             lincoeff,           /**< coefficient of x */
   SCIP_Real             rhs                 /**< right hand side of equation */
   );

/** solves a quadratic equation with interval coefficients
 *
 * Given intervals a, b and c, this function computes an interval that contains all solutions of \f$ a x^2 + b x \in c\f$
 */
extern
void SCIPintervalSolveUnivariateQuadExpression(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs                 /**< right hand side of equation */
   );

/** stores range of bivariate quadratic term in resultant
 * given scalars ax, ay, axy, bx, and by and intervals for x and y, computes interval for \f$ ax x^2 + ay y^2 + axy x y + bx x + by y \f$
 * NOTE: the operations are not applied rounding-safe here
 */
extern
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
 * given scalars ax, ay, axy, bx and by, and intervals for x, y, and rhs,
 * computes \f$ \{ x \in \mathbf{x} : \exists y \in \mathbf{y} : a_x x^2 + a_y y^2 + a_{xy} x y + b_x x + b_y y \in \mathbf{\mbox{rhs}} \} \f$
 * NOTE: the operations are not applied rounding-safe here
 */
extern
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

#ifdef __cplusplus
}
#endif

#endif
