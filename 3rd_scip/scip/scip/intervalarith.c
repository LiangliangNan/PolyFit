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

/**@file   intervalarith.c
 * @brief  interval arithmetics for provable bounds
 * @author Tobias Achterberg
 * @author Stefan Vigerske
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "scip/def.h"
#include "scip/intervalarith.h"
#include "scip/pub_message.h"
#include "scip/misc.h"

#ifdef ROUNDING_FE
#define ROUNDING
/*
 * Linux rounding operations
 */

#include <fenv.h>

/** Linux rounding mode settings */
#define SCIP_ROUND_DOWNWARDS FE_DOWNWARD     /**< round always down */
#define SCIP_ROUND_UPWARDS   FE_UPWARD       /**< round always up */
#define SCIP_ROUND_NEAREST   FE_TONEAREST    /**< round always to nearest */
#define SCIP_ROUND_ZERO      FE_TOWARDZERO   /**< round always towards 0.0 */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return TRUE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE        roundmode           /**< rounding mode to activate */
   )
{
   if( fesetround(roundmode) != 0 )
   {
      SCIPerrorMessage("error setting rounding mode to %d\n", roundmode);
      abort();
   }
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
   return (SCIP_ROUNDMODE)fegetround();
}

#endif



#ifdef ROUNDING_FP
#define ROUNDING
/*
 * OSF rounding operations
 */

#include <float.h>

/** OSF rounding mode settings */
#define SCIP_ROUND_DOWNWARDS FP_RND_RM       /**< round always down */
#define SCIP_ROUND_UPWARDS   FP_RND_RP       /**< round always up */
#define SCIP_ROUND_NEAREST   FP_RND_RN       /**< round always to nearest */
#define SCIP_ROUND_ZERO      FP_RND_RZ       /**< round always towards 0.0 */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return TRUE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE        roundmode           /**< rounding mode to activate */
   )
{
   if( write_rnd(roundmode) != 0 )
   {
      SCIPerrorMessage("error setting rounding mode to %d\n", roundmode);
      abort();
   }
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
   return read_rnd();
}

#endif



#ifdef ROUNDING_MS
#define ROUNDING
/*
 * Microsoft compiler rounding operations
 */

#include <float.h>

/** Microsoft rounding mode settings */
#define SCIP_ROUND_DOWNWARDS RC_DOWN         /**< round always down */
#define SCIP_ROUND_UPWARDS   RC_UP           /**< round always up */
#define SCIP_ROUND_NEAREST   RC_NEAR         /**< round always to nearest */
#define SCIP_ROUND_ZERO      RC_CHOP         /**< round always towards zero */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return TRUE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE        roundmode           /**< rounding mode to activate */
   )
{
   if( (_controlfp(roundmode, _MCW_RC) & _MCW_RC) != roundmode )
   {
      SCIPerrorMessage("error setting rounding mode to %x\n", roundmode);
      abort();
   }
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
   return _controlfp(0, 0) & _MCW_RC;
}
#endif



#ifndef ROUNDING
/*
 * rouding operations not available
 */
#define SCIP_ROUND_DOWNWARDS 0               /**< round always down */
#define SCIP_ROUND_UPWARDS   1               /**< round always up */
#define SCIP_ROUND_NEAREST   2               /**< round always to nearest */
#define SCIP_ROUND_ZERO      3               /**< round always towards zero */

/** returns whether rounding mode control is available */
SCIP_Bool SCIPintervalHasRoundingControl(
   void
   )
{
   return FALSE;
}

/** sets rounding mode of floating point operations */
void SCIPintervalSetRoundingMode(
   SCIP_ROUNDMODE        roundmode           /**< rounding mode to activate */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("setting rounding mode not available - interval arithmetic is invalid!\n");
}

/** gets current rounding mode of floating point operations */
SCIP_ROUNDMODE SCIPintervalGetRoundingMode(
   void
   )
{
   return SCIP_ROUND_NEAREST;
}
#else
#undef ROUNDING
#endif


#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))  /* gcc or icc compiler on x86 32bit or 64bit */

/** gets the negation of a double
 * Do this in a way that the compiler does not "optimize" it away, which usually does not considers rounding modes.
 * However, compiling with -frounding-math would allow to return -x here.
 */
static
double negate(
   /* we explicitely use double here, since I'm not sure the assembler code would work as it for other float's */
   double                x                   /**< number that should be negated */
   )
{
   /* The following line of code is taken from GAOL, http://sourceforge.net/projects/gaol. */
   __asm volatile ("fldl %1; fchs; fstpl %0" : "=m" (x) : "m" (x));
   return x;
}

/* cl or icl compiler on 32bit windows or icl compiler on 64bit windows
 * cl on 64bit windows does not seem to support inline assembler
 */
#elif defined(_MSC_VER) && (defined(__INTEL_COMPILER) || !defined(_M_X64))

/** gets the negation of a double
 * Do this in a way that the compiler does not "optimize" it away, which usually does not considers rounding modes.
 */
static
double negate(
   /* we explicitely use double here, since I'm not sure the assembler code would work as it for other float's */
   double                x                   /**< number that should be negated */
   )
{
   /* The following lines of code are taken from GAOL, http://sourceforge.net/projects/gaol. */
   __asm {
      fld x
         fchs
         fstp x
         }
   return x;
}

#else /* unknown compiler or MSVS 64bit */

/** gets the negation of a double
 *
 * Fallback implementation that calls the negation method from misc.o.
 * Having the implementation in a different object file will hopefully prevent
 * it from being "optimized away".
 */
static
SCIP_Real negate(
   SCIP_Real             x                   /**< number that should be negated */
   )
{
   return SCIPnegateReal(x);
}

#endif


/** sets rounding mode of floating point operations to downwards rounding */
void SCIPintervalSetRoundingModeDownwards(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
}

/** sets rounding mode of floating point operations to upwards rounding */
void SCIPintervalSetRoundingModeUpwards(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
}

/** sets rounding mode of floating point operations to nearest rounding */
void SCIPintervalSetRoundingModeToNearest(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
}

/** sets rounding mode of floating point operations to towards zero rounding */
void SCIPintervalSetRoundingModeTowardsZero(
   void
   )
{
   SCIPintervalSetRoundingMode(SCIP_ROUND_ZERO);
}

/** negates a number in a way that the compiler does not optimize it away */
SCIP_Real SCIPintervalNegateReal(
   SCIP_Real             x                   /**< number to negate */
   )
{
   return negate((double)x);
}

/*
 * Interval arithmetic operations
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPintervalGetInf
#undef SCIPintervalGetSup
#undef SCIPintervalSet
#undef SCIPintervalSetBounds
#undef SCIPintervalSetEmpty
#undef SCIPintervalIsEmpty
#undef SCIPintervalSetEntire
#undef SCIPintervalIsEntire
#undef SCIPintervalIsPositiveInfinity
#undef SCIPintervalIsNegativeInfinity

/** returns infimum of interval */
SCIP_Real SCIPintervalGetInf(
   SCIP_INTERVAL         interval            /**< interval */
   )
{
   return interval.inf;
}

/** returns supremum of interval */
SCIP_Real SCIPintervalGetSup(
   SCIP_INTERVAL         interval            /**< interval */
   )
{
   return interval.sup;
}

/** stores given value as interval */
void SCIPintervalSet(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             value               /**< value to store */
   )
{
   assert(resultant != NULL);

   resultant->inf = value;
   resultant->sup = value;
}

/** stores given infimum and supremum as interval */
void SCIPintervalSetBounds(
   SCIP_INTERVAL*        resultant,          /**< interval to store value into */
   SCIP_Real             inf,                /**< value to store as infimum */
   SCIP_Real             sup                 /**< value to store as supremum */
   )
{
   assert(resultant != NULL);
   assert(inf <= sup);

   resultant->inf = inf;
   resultant->sup = sup;
}

/** sets interval to empty interval, which will be [infinity, -infinity] */
void SCIPintervalSetEmpty(
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   )
{
   assert(resultant != NULL);

   resultant->inf =  1.0;
   resultant->sup = -1.0;
}

/** indicates whether interval is empty, i.e., whether inf > sup */
SCIP_Bool SCIPintervalIsEmpty(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   if( operand.sup >= infinity || operand.inf <= -infinity )
      return FALSE;

   return operand.sup < operand.inf;
}

/** sets interval to entire [-infinity, +infinity] */
void SCIPintervalSetEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant           /**< resultant interval of operation */
   )
{
   assert(resultant != NULL);

   resultant->inf = -infinity;
   resultant->sup =  infinity;
}

/** indicates whether interval is entire, i.e., whether inf <= -infinity and sup >= infinity */
SCIP_Bool SCIPintervalIsEntire(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   return operand.inf <= -infinity && operand.sup >= infinity;
}

/** indicates whether interval is positive infinity, i.e., [infinity, infinity] */
SCIP_Bool SCIPintervalIsPositiveInfinity(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   return operand.inf >=  infinity && operand.sup >= operand.inf;
}

/** indicates whether interval is negative infinity, i.e., [-infinity, -infinity] */
SCIP_Bool SCIPintervalIsNegativeInfinity(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   return operand.sup <= -infinity && operand.inf <= operand.sup;
}

/** indicates whether operand1 is contained in operand2 */
SCIP_Bool SCIPintervalIsSubsetEQ(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   /* the empty interval is contained everywhere */
   if( operand1.inf > operand1.sup )
      return TRUE;

   /* something not-empty is not contained in the empty interval */
   if( operand2.inf > operand2.sup )
      return FALSE;

   return (MAX(-infinity, operand1.inf) >= operand2.inf) &&
      (    MIN( infinity, operand1.sup) <= operand2.sup);
}

/** indicates whether operand1 and operand2 are disjoint */
SCIP_Bool SCIPintervalAreDisjoint(
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   return (operand1.sup < operand2.inf) || (operand2.sup < operand1.inf);
}

/** intersection of two intervals */
void SCIPintervalIntersect(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);

   resultant->inf = MAX(operand1.inf, operand2.inf);
   resultant->sup = MIN(operand1.sup, operand2.sup);
}

/** interval enclosure of the union of two intervals */
void SCIPintervalUnify(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);

   if( operand1.inf > operand1.sup )
   {
      /* operand1 is empty */
      *resultant = operand2;
      return;
   }

   if( operand2.inf > operand2.sup )
   {
      /* operand2 is empty */
      *resultant = operand1;
      return;
   }

   resultant->inf = MIN(operand1.inf, operand2.inf);
   resultant->sup = MAX(operand1.sup, operand2.sup);
}

/** adds operand1 and operand2 and stores infimum of result in infimum of resultant */
void SCIPintervalAddInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS);
   assert(resultant != NULL);

   /* [a,...] + [-inf,...] = [-inf,...] for all a, in particular, [+inf,...] + [-inf,...] = [-inf,...] */
   if( operand1.inf <= -infinity || operand2.inf <= -infinity )
   {
      resultant->inf = -infinity;
   }
   /* [a,...] + [+inf,...] = [+inf,...] for all a > -inf */
   else if( operand1.inf >= infinity || operand2.inf >= infinity )
   {
      resultant->inf = infinity;
   }
   else
   {
      resultant->inf = operand1.inf + operand2.inf;
   }
}

/** adds operand1 and operand2 and stores supremum of result in supremum of resultant */
void SCIPintervalAddSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS);
   assert(resultant != NULL);

   /* [...,b] + [...,+inf] = [...,+inf] for all b, in particular, [...,-inf] + [...,+inf] = [...,+inf] */
   if( operand1.sup >= infinity || operand2.sup >= infinity )
   {
      resultant->sup = infinity;
   }
   /* [...,b] + [...,-inf] = [...,-inf] for all b < +inf */
   else if( operand1.sup <= -infinity || operand2.sup <= -infinity )
   {
      resultant->sup = -infinity;
   }
   else
   {
      resultant->sup = operand1.sup + operand2.sup;
   }
}

/** adds operand1 and operand2 and stores result in resultant */
void SCIPintervalAdd(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalAddInf(infinity, resultant, operand1, operand2);

   /* compute supremum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalAddSup(infinity, resultant, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** adds operand1 and scalar operand2 and stores result in resultant */
void SCIPintervalAddScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));

   roundmode = SCIPintervalGetRoundingMode();

   /* -inf + something >= -inf */
   if( operand1.inf <= -infinity || operand2 <= -infinity )
   {
      resultant->inf = -infinity;
   }
   else if( operand1.inf >= infinity || operand2 >= infinity )
   {
      /* inf + finite = inf, inf + inf = inf */
      resultant->inf = infinity;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      resultant->inf = operand1.inf + operand2;
   }

   /* inf + something <= inf */
   if( operand1.sup >=  infinity || operand2 >= infinity )
   {
      resultant->sup =  infinity;
   }
   else if( operand1.sup <= -infinity || operand2 <= -infinity )
   {
      /* -inf + finite = -inf, -inf + (-inf) = -inf */
      resultant->sup = -infinity;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      resultant->sup = operand1.sup + operand2;
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** adds vector operand1 and vector operand2 and stores result in vector resultant */
void SCIPintervalAddVectors(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< array of resultant intervals of operation */
   int                   length,             /**< length of arrays */
   SCIP_INTERVAL*        operand1,           /**< array of first operands of operation */
   SCIP_INTERVAL*        operand2            /**< array of second operands of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   int i;

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimums of resultant array */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   for( i = 0; i < length; ++i )
   {
      SCIPintervalAddInf(infinity, &resultant[i], operand1[i], operand2[i]);
   }
   /* compute supremums of result array */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   for( i = 0; i < length; ++i )
   {
      SCIPintervalAddSup(infinity, &resultant[i], operand1[i], operand2[i]);
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** subtracts operand2 from operand1 and stores result in resultant */
void SCIPintervalSub(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   roundmode = SCIPintervalGetRoundingMode();

   if( operand1.inf <= -infinity || operand2.sup >=  infinity )
      resultant->inf = -infinity;
   /* [a,b] - [-inf,-inf] = [+inf,+inf] */
   else if( operand1.inf >= infinity || operand2.sup <= -infinity )
   {
      resultant->inf = infinity;
      resultant->sup = infinity;
      return;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      resultant->inf = operand1.inf - operand2.sup;
   }

   if( operand1.sup >=  infinity || operand2.inf <= -infinity )
      resultant->sup =  infinity;
   /* [a,b] - [+inf,+inf] = [-inf,-inf] */
   else if( operand1.sup <= -infinity || operand2.inf >= infinity )
   {
      assert(resultant->inf == -infinity);  /* should be set above, since operand1.inf <= operand1.sup <= -infinity */  /*lint !e777*/
      resultant->sup = -infinity;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      resultant->sup = operand1.sup - operand2.inf;
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** subtracts scalar operand2 from operand1 and stores result in resultant */
void SCIPintervalSubScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIPintervalAddScalar(infinity, resultant, operand1, -operand2);
}

/** multiplies operand1 with operand2 and stores infimum of result in infimum of resultant */
void SCIPintervalMulInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS); 

   if( operand1.inf >= infinity )
   {
      /* operand1 is infinity scalar */
      assert(operand1.sup >= infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand2, infinity);
   }
   else if( operand2.inf >= infinity )
   {
      /* operand2 is infinity scalar */
      assert(operand2.sup >=  infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand1, infinity);
   }
   else if( operand1.sup <= -infinity )
   {
      /* operand1 is -infinity scalar */
      assert(operand1.inf <= -infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand2, -infinity);
   }
   else if( operand2.sup <= -infinity )
   {
      /* operand2 is -infinity scalar */
      assert(operand2.inf <= -infinity);
      SCIPintervalMulScalarInf(infinity, resultant, operand1, -infinity);
   }
   else if( ( operand1.inf <= -infinity && operand2.sup > 0.0 ) 
      || ( operand1.sup > 0.0 && operand2.inf <= -infinity ) 
      || ( operand1.inf < 0.0 && operand2.sup >= infinity ) 
      || ( operand1.sup >= infinity && operand2.inf < 0.0 ) )
   {
      resultant->inf = -infinity;
   }
   else
   {
      SCIP_Real cand1;
      SCIP_Real cand2;
      SCIP_Real cand3;
      SCIP_Real cand4;

      cand1 = operand1.inf * operand2.inf;
      cand2 = operand1.inf * operand2.sup;
      cand3 = operand1.sup * operand2.inf;
      cand4 = operand1.sup * operand2.sup;
      resultant->inf = MIN(MIN(cand1, cand2), MIN(cand3, cand4));
   }
}

/** multiplies operand1 with operand2 and stores supremum of result in supremum of resultant */
void SCIPintervalMulSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS); 

   if( operand1.inf >= infinity )
   {
      /* operand1 is infinity scalar */
      assert(operand1.sup >= infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand2, infinity);
   }
   else if( operand2.inf >= infinity )
   {
      /* operand2 is infinity scalar */
      assert(operand2.sup >=  infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand1, infinity);
   }
   else if( operand1.sup <= -infinity )
   {
      /* operand1 is -infinity scalar */
      assert(operand1.inf <= -infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand2, -infinity);
   }
   else if( operand2.sup <= -infinity )
   {
      /* operand2 is -infinity scalar */
      assert(operand2.inf <= -infinity);
      SCIPintervalMulScalarSup(infinity, resultant, operand1, -infinity);
   }
   else if( ( operand1.inf <= -infinity && operand2.inf < 0.0 ) 
      || ( operand1.inf < 0.0 && operand2.inf <= -infinity ) 
      || ( operand1.sup > 0.0 && operand2.sup >= infinity ) 
      || ( operand1.sup >= infinity && operand2.sup > 0.0 ) )
   {
      resultant->sup =  infinity;
   }
   else
   {
      SCIP_Real cand1;
      SCIP_Real cand2;
      SCIP_Real cand3;
      SCIP_Real cand4;

      cand1 = operand1.inf * operand2.inf;
      cand2 = operand1.inf * operand2.sup;
      cand3 = operand1.sup * operand2.inf;
      cand4 = operand1.sup * operand2.sup;
      resultant->sup = MAX(MAX(cand1, cand2), MAX(cand3, cand4));
   }
}

/** multiplies operand1 with operand2 and stores result in resultant */
void SCIPintervalMul(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation; can be +/-inf */
   SCIP_INTERVAL         operand2            /**< second operand of operation; can be +/-inf */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimum result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalMulInf(infinity, resultant, operand1, operand2);

   /* compute supremum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalMulSup(infinity, resultant, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** multiplies operand1 with scalar operand2 and stores infimum of result in infimum of resultant */
void SCIPintervalMulScalarInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS);
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));

   if( operand2 >= infinity )
   {
      /* result.inf defined by sign of operand1.inf */ 
      if( operand1.inf > 0 )
         resultant->inf = infinity;
      else if( operand1.inf < 0 )
         resultant->inf = -infinity;
      else
         resultant->inf = 0.0;
   }
   else if( operand2 <= -infinity )
   {
      /* result.inf defined by sign of operand1.sup */ 
      if( operand1.sup > 0 )
         resultant->inf = -infinity;
      else if( operand1.sup < 0 )
         resultant->inf = infinity;
      else
         resultant->inf = 0.0;
   }
   else if( operand2 == 0.0 )
   {
      resultant->inf = 0.0;
   }
   else if( operand2 > 0.0 )
   {
      if( operand1.inf <= -infinity )
         resultant->inf = -infinity;
      else if( operand1.inf >= infinity )
         resultant->inf =  infinity;
      else
         resultant->inf = operand1.inf * operand2;
   }
   else
   {
      if( operand1.sup >= infinity )
         resultant->inf = -infinity;
      else if( operand1.sup <= -infinity )
         resultant->inf =  infinity;
      else
         resultant->inf = operand1.sup * operand2;
   }
}

/** multiplies operand1 with scalar operand2 and stores supremum of result in supremum of resultant */
void SCIPintervalMulScalarSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation; can be +/- inf */
   )
{
   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS);
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));

   if( operand2 >= infinity )
   {
      /* result.sup defined by sign of operand1.sup */ 
      if( operand1.sup > 0 )
         resultant->sup = infinity;
      else if( operand1.sup < 0 )
         resultant->sup = -infinity;
      else
         resultant->sup = 0.0;
   }
   else if( operand2 <= -infinity )
   {
      /* result.sup defined by sign of operand1.inf */ 
      if( operand1.inf > 0 )
         resultant->sup = -infinity;
      else if( operand1.inf < 0 )
         resultant->sup = infinity;
      else
         resultant->sup = 0.0;
   }
   else if( operand2 == 0.0 )
   {
      resultant->sup = 0.0;
   }
   else if( operand2 > 0.0 )
   {
      if( operand1.sup >= infinity )
         resultant->sup = infinity;
      else if( operand1.sup <= -infinity )
         resultant->sup = -infinity;
      else
         resultant->sup = operand1.sup * operand2;
   }
   else
   {
      if( operand1.inf <= -infinity )
         resultant->sup = infinity;
      else if( operand1.inf >= infinity )
         resultant->sup = -infinity;
      else
         resultant->sup = operand1.inf * operand2;
   }
}

/** multiplies operand1 with scalar operand2 and stores result in resultant */
void SCIPintervalMulScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));

   roundmode = SCIPintervalGetRoundingMode();

   /* compute infimum result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalMulScalarInf(infinity, resultant, operand1, operand2);

   /* compute supremum of result */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalMulScalarSup(infinity, resultant, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** divides operand1 by operand2 and stores result in resultant */
void SCIPintervalDiv(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL intmed;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   if( operand2.inf <= 0.0 && operand2.sup >= 0.0 )
   {  /* division by [0,0] or interval containing 0 gives [-inf, +inf] */
      resultant->inf = -infinity;
      resultant->sup =  infinity;
      return;
   }

   if( operand1.inf == 0.0 && operand1.sup == 0.0 )
   {  /* division of [0,0] by something nonzero */
      SCIPintervalSet(resultant, 0.0);
      return;
   }

   roundmode = SCIPintervalGetRoundingMode();

   /* division by nonzero: resultant = x * (1/y) */
   if( operand2.sup >=  infinity || operand2.sup <= -infinity )
   {
      intmed.inf = 0.0;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      intmed.inf = 1.0 / operand2.sup;
   }
   if( operand2.inf <= -infinity || operand2.inf >= infinity )
   {
      intmed.sup = 0.0;
   }
   else
   {
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      intmed.sup = 1.0 / operand2.inf;
   }
   SCIPintervalMul(infinity, resultant, operand1, intmed);

   SCIPintervalSetRoundingMode(roundmode);
}

/** divides operand1 by scalar operand2 and stores result in resultant */
void SCIPintervalDivScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));

   roundmode = SCIPintervalGetRoundingMode();

   if( operand2 >= infinity || operand2 <= -infinity )
   {
      /* division by +/-infinity is 0.0 */
      resultant->inf = 0.0;
      resultant->sup = 0.0;
   }
   else if( operand2 > 0.0 )
   {
      if( operand1.inf <= -infinity )
         resultant->inf = -infinity;
      else if( operand1.inf >= infinity )
      {
         /* infinity / + = infinity */
         resultant->inf = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf / operand2;
      }

      if( operand1.sup >= infinity )
         resultant->sup =  infinity;
      else if( operand1.sup <= -infinity )
      {
         /* -infinity / + = -infinity */
         resultant->sup = -infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup / operand2;
      }
   }
   else if( operand2 < 0.0 )
   {
      if( operand1.sup >=  infinity )
         resultant->inf = -infinity;
      else if( operand1.sup <= -infinity )
      {
         /* -infinity / - = infinity */
         resultant->inf = infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.sup / operand2;
      }

      if( operand1.inf <= -infinity )
         resultant->sup = infinity;
      else if( operand1.inf >= infinity )
      {
         /* infinity / - = -infinity */
         resultant->sup = -infinity;
      }
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.inf / operand2;
      }
   }
   else
   { /* division by 0.0 */
      if( operand1.inf >= 0 )
      {
         /* [+,+] / [0,0] = [+inf, +inf] */
         resultant->inf =  infinity;
         resultant->sup =  infinity;
      }
      else if( operand1.sup <= 0 )
      {
         /* [-,-] / [0,0] = [-inf, -inf] */
         resultant->inf = -infinity;
         resultant->sup = -infinity;
      }
      else
      {
         /* [-,+] / [0,0] = [-inf, +inf] */
         resultant->inf = -infinity;
         resultant->sup =  infinity;
      }
      return;
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** computes the scalar product of two vectors of intervals and stores result in resultant */
void SCIPintervalScalprod(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals; can have +/-inf entries */
   SCIP_INTERVAL*        operand2            /**< second vector as array of intervals; can have +/-inf entries */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL prod;
   int i;

   roundmode = SCIPintervalGetRoundingMode();

   resultant->inf = 0.0;
   resultant->sup = 0.0;

   /* compute infimum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   for( i = 0; i < length && resultant->inf > -infinity; ++i )
   {
      SCIPintervalSetEntire(infinity, &prod);
      SCIPintervalMulInf(infinity, &prod, operand1[i], operand2[i]);
      SCIPintervalAddInf(infinity, resultant, *resultant, prod); 
   }
   assert(resultant->sup == 0.0);

   /* compute supremum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   for( i = 0; i < length && resultant->sup < infinity ; ++i )
   {
      SCIPintervalSetEntire(infinity, &prod);
      SCIPintervalMulSup(infinity, &prod, operand1[i], operand2[i]);
      SCIPintervalAddSup(infinity, resultant, *resultant, prod); 
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** computes scalar product of a vector of intervals and a vector of scalars and stores infimum of result in infimum of 
 *  resultant 
 */
void SCIPintervalScalprodScalarsInf(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   )
{
   SCIP_INTERVAL prod;
   int i;

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_DOWNWARDS);

   resultant->inf = 0.0;

   /* compute infimum of resultant */
   SCIPintervalSetEntire(infinity, &prod);
   for( i = 0; i < length && resultant->inf > -infinity; ++i )
   {
      SCIPintervalMulScalarInf(infinity, &prod, operand1[i], operand2[i]);
      assert(prod.sup >= infinity);
      SCIPintervalAddInf(infinity, resultant, *resultant, prod); 
   }
}

/** computes the scalar product of a vector of intervals and a vector of scalars and stores supremum of result in 
 *  supremum of resultant 
 */
void SCIPintervalScalprodScalarsSup(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   )
{
   SCIP_INTERVAL prod;
   int i;

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_UPWARDS);

   resultant->sup = 0.0;

   /* compute supremum of resultant */
   SCIPintervalSetEntire(infinity, &prod);
   for( i = 0; i < length && resultant->sup < infinity; ++i )
   {
      SCIPintervalMulScalarSup(infinity, &prod, operand1[i], operand2[i]);
      assert(prod.inf <= -infinity);
      SCIPintervalAddSup(infinity, resultant, *resultant, prod); 
   }
}

/** computes the scalar product of a vector of intervals and a vector of scalars and stores result in resultant */
void SCIPintervalScalprodScalars(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   int                   length,             /**< length of vectors */
   SCIP_INTERVAL*        operand1,           /**< first vector as array of intervals */
   SCIP_Real*            operand2            /**< second vector as array of scalars; can have +/-inf entries */
   )
{
   SCIP_ROUNDMODE roundmode;

   roundmode = SCIPintervalGetRoundingMode();

   resultant->inf = 0.0;
   resultant->sup = 0.0;

   /* compute infimum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
   SCIPintervalScalprodScalarsInf(infinity, resultant, length, operand1, operand2);
   assert(resultant->sup == 0.0);

   /* compute supremum of resultant */
   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   SCIPintervalScalprodScalarsSup(infinity, resultant, length, operand1, operand2);

   SCIPintervalSetRoundingMode(roundmode);
}

/** squares operand and stores result in resultant */
void SCIPintervalSquare(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   roundmode = SCIPintervalGetRoundingMode();

   if( operand.sup <= 0.0 )
   {  /* operand is left of 0.0 */
      if( operand.sup <= -infinity )
         resultant->inf =  infinity;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand.sup * operand.sup;
      }

      if( operand.inf <= -infinity )
         resultant->sup = infinity;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand.inf * operand.inf;
      }
   }
   else if( operand.inf >= 0.0 )
   {  /* operand is right of 0.0 */
      if( operand.inf >= infinity )
         resultant->inf = infinity;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand.inf * operand.inf;
      }

      if( operand.sup >= infinity )
         resultant->sup = infinity;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand.sup * operand.sup;
      }
   }
   else
   {  /* [-,+]^2 */
      resultant->inf = 0.0;
      if( operand.inf <= -infinity || operand.sup >= infinity )
         resultant->sup = infinity;
      else
      {
         SCIP_Real x;
         SCIP_Real y;

         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         x = operand.inf * operand.inf;
         y = operand.sup * operand.sup;
         resultant->sup = MAX(x, y);
      }
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** stores (positive part of) square root of operand in resultant
 * @attention we assume a correctly rounded sqrt(double) function when rounding is to nearest
 */
void SCIPintervalSquareRoot(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   if( operand.sup < 0.0 )
   {
      SCIPintervalSetEmpty(resultant);
      return;
   }

   if( operand.inf == operand.sup )  /*lint !e777 */
   {
      if( operand.inf >= infinity )
      {
         resultant->inf = infinity;
         resultant->sup = infinity;
      }
      else
      {
         SCIP_Real tmp;

         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
         tmp = sqrt(operand.inf);
         resultant->inf = SCIPnextafter(tmp, SCIP_REAL_MIN);
         resultant->sup = SCIPnextafter(tmp, SCIP_REAL_MAX);
      }

      return;
   }

   if( operand.inf <= 0.0 )
      resultant->inf = 0.0;
   else if( operand.inf >= infinity )
   {
      resultant->inf = infinity;
      resultant->sup = infinity;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
      resultant->inf = SCIPnextafter(sqrt(operand.inf), SCIP_REAL_MIN);
   }

   if( operand.sup >= infinity )
      resultant->sup = infinity;
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
      resultant->sup = SCIPnextafter(sqrt(operand.sup), SCIP_REAL_MAX);
   }
}

/** stores operand1 to the power of operand2 in resultant
 * 
 * uses SCIPintervalPowerScalar if operand2 is a scalar, otherwise computes exp(op2*log(op1))
 */
void SCIPintervalPower(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{  /*lint --e{777}*/
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   if( operand2.inf == operand2.sup )
   {  /* operand is number */
      SCIPintervalPowerScalar(infinity, resultant, operand1, operand2.inf);
      return;
   }

   /* resultant := log(op1) */
   SCIPintervalLog(infinity, resultant, operand1);
   if( SCIPintervalIsEmpty(infinity, *resultant) )
      return;

   /* resultant := op2 * resultant */
   SCIPintervalMul(infinity, resultant, operand2, *resultant);

   /* resultant := exp(resultant) */
   SCIPintervalExp(infinity, resultant, *resultant);
}

/** computes lower bound on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
SCIP_Real SCIPintervalPowerScalarIntegerInf(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_Real result;

   assert(operand1 >= 0.0);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
         return 1.0; /* 0^0 = 1 */
      else
         return 0.0; /* 0^positive = 0 */
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
      return 1.0;

   if( operand2 < 0 )
   {
      /* x^n = 1 / x^(-n) */
      result = SCIPintervalPowerScalarIntegerSup(operand1, -operand2);
      assert(result != 0.0);

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
      result = 1.0 / result;
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      unsigned int n;
      SCIP_Real z;

      roundmode = SCIPintervalGetRoundingMode();

      result = 1.0;
      n = (unsigned int)operand2;
      z = operand1;

      SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);

      /* use a binary exponentiation algorithm:
       * consider the binary representation of n: n = sum_i 2^i d_i with d_i \in {0,1}
       * then x^n = prod_{i:d_i=1} x^(2^i)
       * in the following, we loop over i=1,..., thereby storing x^(2^i) in z
       * whenever d_i is 1, we multiply result with x^(2^i) (the current z)
       * at the last (highest) i with d_i = 1 we stop, thus having x^n stored in result
       *
       * the binary representation of n and bit shifting is used for the loop
       */
      assert(n >= 1);
      do
      {
         if( n & 1 ) /* n is odd (d_i=1), so multiply result with current z (=x^{2^i}) */
         {
            result = result * z;
            n >>= 1;
            if( n == 0 )
               break;
         }
         else
            n >>= 1;
         z = z * z;
      }
      while( TRUE );  /*lint !e506 */

      SCIPintervalSetRoundingMode(roundmode);
   }

   return result;
}

/** computes upper bound on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
SCIP_Real SCIPintervalPowerScalarIntegerSup(
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_Real result;

   assert(operand1 >= 0.0);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
         return 1.0; /* 0^0 = 1 */
      else
         return 0.0; /* 0^positive = 0 */
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
      return 1.0;

   if( operand2 < 0 )
   {
      /* x^n = 1 / x^(-n) */
      result = SCIPintervalPowerScalarIntegerInf(operand1, -operand2);
      assert(result != 0.0);

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
      result = 1.0 / result;
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      unsigned int n;
      SCIP_Real z;

      roundmode = SCIPintervalGetRoundingMode();

      result = 1.0;
      n = (unsigned int)operand2;
      z = operand1;

      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);

      /* use a binary exponentiation algorithm... see comments in SCIPintervalPowerScalarIntegerInf */
      assert(n >= 1);
      do
      {
         if( n&1 )
         {
            result = result * z;
            n >>= 1;
            if( n == 0 )
               break;
         }
         else
            n >>= 1;
         z = z * z;
      }
      while( TRUE );  /*lint !e506 */

      SCIPintervalSetRoundingMode(roundmode);
   }

   return result;
}

/** computes bounds on power of a scalar operand1 to an integer operand2
 * both operands need to be finite numbers
 * need to have operand1 >= 0 and need to have operand2 >= 0 if operand1 == 0
 */
void SCIPintervalPowerScalarInteger(
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   int                   operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(operand1 >= 0.0);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
      {
         SCIPintervalSet(resultant, 1.0); /* 0^0 = 1 */
         return;
      }
      else
      {
         SCIPintervalSet(resultant, 0.0); /* 0^positive = 0 */
         return;
      }
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
   {
      SCIPintervalSet(resultant, 1.0);
      return;
   }

   if( operand2 < 0 )
   {
      /* x^n = 1 / x^(-n) */
      SCIPintervalPowerScalarInteger(resultant, operand1, -operand2);
      assert(resultant->inf > 0.0 || resultant->sup < 0.0);
      SCIPintervalReciprocal(SCIP_REAL_MAX, resultant, *resultant); /* value for infinity does not matter, since there should be no 0.0 in the interval, so just use something large enough */
   }
   else
   {
      unsigned int n;
      SCIP_Real z_inf;
      SCIP_Real z_sup;
      SCIP_Real result_sup;
      SCIP_Real result_inf;

      roundmode = SCIPintervalGetRoundingMode();

      result_inf = 1.0;
      result_sup = 1.0;
      z_inf = operand1;
      z_sup = operand1;
      n = (unsigned int)operand2;

      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);

      /* use a binary exponentiation algorithm... see comments in SCIPintervalPowerScalarIntegerInf
       * we compute lower and upper bounds within the same loop
       * to get correct lower bounds while rounding mode is upwards, we negate arguments */
      assert(n >= 1);
      do
      {
         if( n & 1 )
         {
            result_inf = negate(negate(result_inf) * z_inf);
            result_sup = result_sup * z_sup;
            n >>= 1;
            if( n == 0 )
               break;
         }
         else
            n >>= 1;
         z_inf = negate(negate(z_inf) * z_inf);
         z_sup = z_sup * z_sup;
      }
      while( TRUE );  /*lint !e506 */

      SCIPintervalSetRoundingMode(roundmode);

      resultant->inf = result_inf;
      resultant->sup = result_sup;
      assert(resultant->inf <= resultant->sup);
   }
}

/** stores bounds on the power of a scalar operand1 to a scalar operand2 in resultant
 * both operands need to be finite numbers
 * need to have operand1 >= 0 or operand2 integer and need to have operand2 >= 0 if operand1 == 0
 * @attention we assume a correctly rounded pow(double) function when rounding is to nearest
 */
void SCIPintervalPowerScalarScalar(
   SCIP_INTERVAL*        resultant,          /**< resultant of operation */
   SCIP_Real             operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_Real result;

   assert(resultant != NULL);

   if( operand1 == 0.0 )
   {
      assert(operand2 >= 0);
      if( operand2 == 0 )
      {
         SCIPintervalSet(resultant, 1.0); /* 0^0 = 1 */
         return;
      }
      else
      {
         SCIPintervalSet(resultant, 0.0); /* 0^positive = 0 */
         return;
      }
   }

   /* 1^n = 1, x^0 = 1 */
   if( operand1 == 1.0 || operand2 == 0 )
   {
      SCIPintervalSet(resultant, 1.0);
      return;
   }

   assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
   result = pow(operand1, operand2);

   /* to get safe bounds, get the floating point numbers just below and above result */
   resultant->inf = SCIPnextafter(result, SCIP_REAL_MIN);
   resultant->sup = SCIPnextafter(result, SCIP_REAL_MAX);
}

/** stores operand1 to the power of the scalar operand2 in resultant
 * @attention we assume a correctly rounded pow(double) function when rounding is to nearest
 */
void SCIPintervalPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{  /*lint --e{777}*/
   SCIP_Bool op2isint;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));

   if( operand2 == infinity )
   {
      /* 0^infinity =  0
       * +^infinity =  infinity
       * -^infinity = -infinity
       */
      if( operand1.inf < 0.0 )
         resultant->inf = -infinity;
      else
         resultant->inf = 0.0;
      if( operand1.sup > 0.0 )
         resultant->sup =  infinity;
      else
         resultant->sup = 0.0;
      return;
   }

   if( operand2 == 0.0 )
   { /* special case, since x^0 = 1 for x != 0, but 0^0 = 0 */
      if( operand1.inf == 0.0 && operand1.sup == 0.0 )
      {
         resultant->inf = 0.0;
         resultant->sup = 0.0;
      }
      else if( operand1.inf <= 0.0 || operand1.sup >= 0.0 )
      { /* 0.0 in x gives [0,1] */
         resultant->inf = 0.0;
         resultant->sup = 1.0;
      }
      else
      { /* 0.0 outside x gives [1,1] */
         resultant->inf = 1.0;
         resultant->sup = 1.0;
      }
      return;
   }

   if( operand2 == 1.0 )
   {
      /* x^1 = x */
      *resultant = operand1;
      return;
   }

   op2isint = (ceil(operand2) == operand2);

   if( !op2isint && operand1.inf < 0.0 )
   {  /* x^n with x negative not defined for n not integer*/
      operand1.inf = 0.0;
      if( operand1.sup < operand1.inf )
      {
         SCIPintervalSetEmpty(resultant);
         return;
      }
   }

   if( operand1.inf >= 0.0 )
   {  /* easy case: x^n with x >= 0 */
      if( operand2 >= 0.0 )
      {
         /* inf^+ = inf */
         if( operand1.inf >= infinity )
            resultant->inf = infinity;
         else if( operand1.inf > 0.0 )
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
            resultant->inf = SCIPnextafter(pow(operand1.inf, operand2), SCIP_REAL_MIN);
         }
         else
            resultant->inf = 0.0;

         if( operand1.sup >= infinity )
            resultant->sup = infinity;
         else if( operand1.sup > 0.0 )
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
            resultant->sup = SCIPnextafter(pow(operand1.sup, operand2), SCIP_REAL_MAX);
         }
         else
            resultant->sup = 0.0;
      }
      else
      {
         if( operand1.sup >= infinity )
            resultant->inf = 0.0;
         else if( operand1.sup == 0.0 )
         {
            /* x^(negative even) = infinity for x->0 (from both sides),
             * but x^(negative odd) = -infinity for x->0 from left side */
            if( ceil(operand2/2) == operand2/2 )
               resultant->inf =  infinity;
            else
               resultant->inf = -infinity;
         }
         else
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
            resultant->inf = SCIPnextafter(pow(operand1.sup, operand2), SCIP_REAL_MIN);
         }

         /* 0^(negative) = infinity */
         if( operand1.inf == 0.0 )
            resultant->sup = infinity;
         else
         {
            assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST); /* usually, no-one should have changed rounding mode */
            resultant->sup = SCIPnextafter(pow(operand1.inf, operand2), SCIP_REAL_MAX);
         }
      }
   }
   else if( operand1.sup <= 0.0 )
   {  /* more difficult case: x^n with x < 0; we now know, that n is integer */
      assert(op2isint);
      if( operand2 >= 0.0 && ceil(operand2/2) == operand2/2 )
      {
         /* x^n with n>=2 and even -> x^n is monotonically decreasing for x < 0 */
         if( operand1.sup == -infinity )
            /* (-inf)^n = inf */
            resultant->inf = infinity;
         else
            resultant->inf = SCIPintervalPowerScalarIntegerInf(-operand1.sup, (int)operand2);

         if( operand1.inf <= -infinity )
            resultant->sup = infinity;
         else
            resultant->sup = SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);
      }
      else if( operand2 <= 0.0 && ceil(operand2/2) != operand2/2 )
      {
         /* x^n with n<=-1 and odd -> x^n = 1/x^(-n) is monotonically decreasing for x<0 */
         if( operand1.sup == -infinity )
            /* (-inf)^n = 1/(-inf)^(-n) = 1/(-inf) = 0 */
            resultant->inf = 0.0;
         else if( operand1.sup == 0.0 )
            /* x^n -> -infinity for x->0 from left */
            resultant->inf = -infinity;
         else
            resultant->inf = -SCIPintervalPowerScalarIntegerSup(-operand1.sup, (int)operand2);

         if( operand1.inf <= -infinity )
            /* (-inf)^n = 1/(-inf)^(-n) = 1/(-inf) = 0 */
            resultant->sup = 0.0;
         else if( operand1.inf == 0.0 )
            /* x^n -> infinity for x->0 from right */
            resultant->sup = infinity;
         else
            resultant->sup = -SCIPintervalPowerScalarIntegerInf(-operand1.inf, (int)operand2);
      }
      else if( operand2 >= 0.0 )
      {
         /* x^n with n>0 and odd -> x^n is monotonically increasing for x<0 */
         if( operand1.inf <= -infinity )
            resultant->inf = -infinity;
         else
            resultant->inf = -SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);

         if( operand1.sup <= -infinity )
            resultant->sup = -infinity;
         else
            resultant->sup = -SCIPintervalPowerScalarIntegerInf(-operand1.sup, (int)operand2);
      }
      else
      {
         /* x^n with n<0 and even -> x^n is monotonically increasing for x<0 */
         if( operand1.inf <= -infinity )
            resultant->inf = 0.0;
         else if( operand1.inf == 0.0 )
            /* x^n -> infinity for x->0 from both sides */
            resultant->inf = infinity;
         else
            resultant->inf = SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);

         if( operand1.sup <= -infinity )
            resultant->sup = 0.0;
         else if( operand1.sup == 0.0 )
            /* x^n -> infinity for x->0 from both sides */
            resultant->sup = infinity;
         else
            resultant->sup = SCIPintervalPowerScalarIntegerSup(-operand1.sup, (int)operand2);
      }
      assert(resultant->inf <= resultant->sup || resultant->inf >= infinity || resultant->sup <= -infinity);
   }
   else
   {  /* similar difficult case: x^n with x in [<0, >0], but n is integer */
      assert(op2isint); /* otherwise we had set operand1.inf == 0.0, which was handled in first case */
      if( operand2 >= 0.0 && operand2/2 == ceil(operand2/2) )
      {
         /* n even positive integer */
         resultant->inf = 0.0;
         if( operand1.inf == -infinity || operand1.sup == infinity )
            resultant->sup = infinity;
         else
            resultant->sup = SCIPintervalPowerScalarIntegerSup(MAX(-operand1.inf, operand1.sup), (int)operand2);
      }
      else if( operand2 <= 0.0 && ceil(operand2/2) == operand2/2 )
      {
         /* n even negative integer */
         resultant->sup = infinity;  /* since 0^n = infinity */
         if( operand1.inf == -infinity || operand1.sup == infinity )
            resultant->inf = 0.0;
         else
            resultant->inf = SCIPintervalPowerScalarIntegerInf(MAX(-operand1.inf, operand1.sup), (int)operand2);
      }
      else if( operand2 >= 0.0 )
      {
         /* n odd positive integer, so monotonically increasing function */
         if( operand1.inf == -infinity )
            resultant->inf = -infinity;
         else
            resultant->inf = -SCIPintervalPowerScalarIntegerSup(-operand1.inf, (int)operand2);
         if( operand1.sup == infinity )
            resultant->sup = infinity;
         else
            resultant->sup = SCIPintervalPowerScalarIntegerSup(operand1.sup, (int)operand2);
      }
      else
      {
         /* n odd negative integer:
          * x^n -> -infinity for x->0 from left
          * x^n ->  infinity for x->0 from right */
         resultant->inf = -infinity;
         resultant->sup =  infinity;
      }
   }

   /* if value for infinity is too small, relax intervals so they do not appear empty */
   if( resultant->inf > infinity )
      resultant->inf = infinity;
   if( resultant->sup < -infinity )
      resultant->sup = -infinity;
}

/** given an interval for the image of a power operation, computes an interval for the origin
 * that is, for y = x^p with p = exponent a given scalar and y = image a given interval,
 * computes a subinterval x of basedomain such that y in x^p and such that for all z in basedomain less x, z^p not in y
 */
void SCIPintervalPowerScalarInverse(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         basedomain,         /**< domain of base */
   SCIP_Real             exponent,           /**< exponent */
   SCIP_INTERVAL         image               /**< interval image of power */
   )
{
   SCIP_INTERVAL tmp;
   SCIP_INTERVAL exprecip;

   assert(resultant != NULL);
   assert(image.inf <= image.sup);
   assert(basedomain.inf <= basedomain.sup);

   if( exponent == 0.0 )
   {
      /* exponent is 0.0 */
      if( image.inf <= 1.0 && image.sup >= 1.0 )
      {
         /* 1 in image -> resultant = entire */
         *resultant = basedomain;
      }
      else if( image.inf <= 0.0 && image.sup >= 0.0 )
      {
         /* 0 in image, 1 not in image -> resultant = 0   (provided 0^0 = 0 ???)
          * -> resultant = {0} intersected with basedomain */
         SCIPintervalSetBounds(resultant, MAX(0.0, basedomain.inf), MIN(0.0, basedomain.sup));
      }
      else
      {
         /* 0 and 1 not in image -> resultant = empty */
         SCIPintervalSetEmpty(resultant);
      }
      return;
   }

   /* i = b^e
    *   i >= 0 -> b = i^(1/e) [union -i^(1/e), if e is even]
    *   i < 0, e odd integer -> b = -(-i)^(1/e)
    *   i < 0, e even integer or fractional -> empty
    */

   SCIPintervalSetBounds(&exprecip, exponent, exponent);
   SCIPintervalReciprocal(infinity, &exprecip, exprecip);

   /* invert positive part of image, if any */
   if( image.sup >= 0.0 )
   {
      SCIPintervalSetBounds(&tmp, MAX(image.inf, 0.0), image.sup);
      SCIPintervalPower(infinity, resultant, tmp, exprecip);
      if( basedomain.inf <= -resultant->inf && EPSISINT(exponent, 0.0) && (int)exponent % 2 == 0 )  /*lint !e835 */
      {
         if( basedomain.sup < resultant->inf )
            SCIPintervalSetBounds(resultant, -resultant->sup, -resultant->inf);
         else
            SCIPintervalSetBounds(resultant, -resultant->sup, resultant->sup);
      }

      SCIPintervalIntersect(resultant, *resultant, basedomain);
   }
   else
      SCIPintervalSetEmpty(resultant);

   /* invert negative part of image, if any and if base can take negative value and if exponent is such that negative values are possible */
   if( image.inf < 0.0 && basedomain.inf < 0.0 && EPSISINT(exponent, 0.0) && ((int)exponent % 2 != 0) )  /*lint !e835 */
   {
      SCIPintervalSetBounds(&tmp, MAX(-image.sup, 0.0), -image.inf);
      SCIPintervalPower(infinity, &tmp, tmp, exprecip);
      SCIPintervalSetBounds(&tmp, -tmp.sup, -tmp.inf);
      SCIPintervalIntersect(&tmp, basedomain, tmp);
      SCIPintervalUnify(resultant, *resultant, tmp);
   }
}

/** stores operand1 to the signed power of the scalar positive operand2 in resultant 
 * 
 * the signed power of x w.r.t. an exponent n >= 0 is given as sign(x) * abs(x)^n
 *
 * @attention we assume correctly rounded sqrt(double) and pow(double) functions when rounding is to nearest
 */
void SCIPintervalSignPowerScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_Real             operand2            /**< second operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;
   assert(resultant != NULL);

   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(operand2     >= 0.0);

   if( operand2 == infinity )  /*lint !e777 */
   {
      /* 0^infinity =  0
       * +^infinity =  infinity
       *-+^infinity = -infinity
       */
      if( operand1.inf < 0.0 )
         resultant->inf = -infinity;
      else
         resultant->inf = 0.0;
      if( operand1.sup > 0.0 )
         resultant->sup =  infinity;
      else
         resultant->sup = 0.0;
      return;
   }

   if( operand2 == 0.0 )
   {
      /* special case, since x^0 = 1 for x != 0, but 0^0 = 0 */
      if( operand1.inf < 0.0 )
         resultant->inf = -1.0;
      else if( operand1.inf == 0.0 )
         resultant->inf =  0.0;
      else
         resultant->inf =  1.0;

      if( operand1.sup < 0.0 )
         resultant->sup = -1.0;
      else if( operand1.sup == 0.0 )
         resultant->sup =  0.0;
      else
         resultant->sup =  1.0;

      return;
   }

   if( operand2 == 1.0 )
   { /* easy case that should run fast */
      *resultant = operand1;
      return;
   }

   roundmode = SCIPintervalGetRoundingMode();

   if( operand2 == 2.0 )
   { /* common case where pow can easily be avoided */
      if( operand1.inf <= -infinity )
      {
         resultant->inf = -infinity;
      }
      else if( operand1.inf >= infinity )
      {
         resultant->inf =  infinity;
      }
      else if( operand1.inf > 0.0 )
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = operand1.inf * operand1.inf;
      }
      else
      {
         /* need upwards since we negate result of multiplication */
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->inf = negate(operand1.inf * operand1.inf);
      }

      if( operand1.sup >=  infinity )
      {
         resultant->sup =  infinity;
      }
      else if( operand1.sup <= -infinity )
      {
         resultant->sup = -infinity;
      }
      else if( operand1.sup > 0.0 )
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = operand1.sup * operand1.sup;
      }
      else
      {
         /* need downwards since we negate result of multiplication */
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->sup = negate(operand1.sup * operand1.sup);
      }
      assert(resultant->inf <= resultant->sup);
   }
   else if( operand2 == 0.5 )
   { /* another common case where pow can easily be avoided */
      if( operand1.inf <= -infinity )
         resultant->inf = -infinity;
      else if( operand1.inf >= infinity )
         resultant->inf = infinity;
      else if( operand1.inf >= 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf =  SCIPnextafter(sqrt( operand1.inf), SCIP_REAL_MIN);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf = -SCIPnextafter(sqrt(-operand1.inf), SCIP_REAL_MAX);
      }

      if( operand1.sup >=  infinity )
         resultant->sup =  infinity;
      else if( operand1.sup <= -infinity )
         resultant->sup = -infinity;
      else if( operand1.sup > 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup =  SCIPnextafter(sqrt( operand1.sup), SCIP_REAL_MAX);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup = -SCIPnextafter(sqrt(-operand1.sup), SCIP_REAL_MAX);
      }
      assert(resultant->inf <= resultant->sup);
   }
   else
   {
      if( operand1.inf <= -infinity )
         resultant->inf = -infinity;
      else if( operand1.inf >= infinity )
         resultant->inf =  infinity;
      else if( operand1.inf > 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf =  SCIPnextafter(pow( operand1.inf, operand2), SCIP_REAL_MIN);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->inf = -SCIPnextafter(pow(-operand1.inf, operand2), SCIP_REAL_MAX);
      }

      if( operand1.sup >=  infinity )
         resultant->sup =  infinity;
      else if( operand1.sup <= -infinity )
         resultant->sup = -infinity;
      else if( operand1.sup > 0.0 )
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup =  SCIPnextafter(pow( operand1.sup, operand2), SCIP_REAL_MAX);
      }
      else
      {
         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         resultant->sup = -SCIPnextafter(pow(-operand1.sup, operand2), SCIP_REAL_MIN);
      }
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** computes the reciprocal of an interval
 */
void SCIPintervalReciprocal(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   SCIP_ROUNDMODE roundmode;

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   if( operand.inf == 0.0 && operand.sup == 0.0 )
   { /* 1/0 = [-inf,inf] */
      resultant->inf =  infinity;
      resultant->sup = -infinity;
      return;
   }

   roundmode = SCIPintervalGetRoundingMode();

   if( operand.inf >= 0.0 )
   {  /* 1/x with x >= 0 */
      if( operand.sup >= infinity )
         resultant->inf = 0.0;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = 1.0 / operand.sup;
      }

      if( operand.inf >= infinity )
         resultant->sup = 0.0;
      else if( operand.inf == 0.0 )
         resultant->sup = infinity;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = 1.0 / operand.inf;
      }

      SCIPintervalSetRoundingMode(roundmode);
   }
   else if( operand.sup <= 0.0 )
   {  /* 1/x with x <= 0 */
      if( operand.sup <= -infinity )
         resultant->inf = 0.0;
      else if( operand.sup == 0.0 )
         resultant->inf = -infinity;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
         resultant->inf = 1.0 / operand.sup;
      }

      if( operand.inf <= -infinity )
         resultant->sup = infinity;
      else
      {
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         resultant->sup = 1.0 / operand.inf;
      }
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {  /* 1/x with x in [-,+] is division by zero */
      resultant->inf = -infinity;
      resultant->sup =  infinity;
   }
}

/** stores exponential of operand in resultant
 * @attention we assume a correctly rounded exp(double) function when rounding is to nearest
 */
void SCIPintervalExp(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   if( operand.sup <= -infinity )
   {
      resultant->inf = 0.0;
      resultant->sup = 0.0;
      return;
   }

   if( operand.inf >=  infinity )
   {
      resultant->inf = infinity;
      resultant->sup = infinity;
      return;
   }

   if( operand.inf == operand.sup )  /*lint !e777 */
   {
      if( operand.inf == 0.0 )
      {
         resultant->inf = 1.0;
         resultant->sup = 1.0;
      }
      else
      {
         SCIP_Real tmp;

         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         tmp = exp(operand.inf);
         resultant->inf = tmp > 0.0 ? SCIPnextafter(tmp, SCIP_REAL_MIN) : 0.0;
         assert(resultant->inf >= 0.0);
         resultant->sup = SCIPnextafter(tmp, SCIP_REAL_MAX);

         return;
      }
   }

   if( operand.inf <= -infinity )
   {
      resultant->inf = 0.0;
   }
   else if( operand.inf == 0.0 )
   {
      resultant->inf = 1.0;
   }
   else
   {
      SCIP_Real tmp;

      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      tmp = exp(operand.inf);
      resultant->inf = tmp > 0.0 ? SCIPnextafter(tmp, SCIP_REAL_MIN) : 0.0;
      /* make sure we do not exceed value for infinity, so interval is not declared as empty if inf and sup are both > infinity */
      if( resultant->inf >= infinity )
         resultant->inf = infinity;
   }

   if( operand.sup >=  infinity )
   {
      resultant->sup = infinity;
   }
   else if( operand.sup == 0.0 )
   {
      resultant->sup = 1.0;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      resultant->sup = SCIPnextafter(exp(operand.sup), SCIP_REAL_MAX);
      if( resultant->sup < -infinity )
         resultant->sup = -infinity;
   }
}

/** stores natural logarithm of operand in resultant
 * @attention we assume a correctly rounded log(double) function when rounding is to nearest
 */
void SCIPintervalLog(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   /* if operand.sup == 0.0, we could return -inf in resultant->sup, but that
    * seems of little use and just creates problems somewhere else, e.g., #1230
    */
   if( operand.sup <= 0.0 )
   {
      SCIPintervalSetEmpty(resultant);
      return;
   }

   if( operand.inf == operand.sup )  /*lint !e777 */
   {
      if( operand.sup == 1.0 )
      {
         resultant->inf = 0.0;
         resultant->sup = 0.0;
      }
      else
      {
         SCIP_Real tmp;

         assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
         tmp = log(operand.inf);
         resultant->inf = SCIPnextafter(tmp, SCIP_REAL_MIN);
         resultant->sup = SCIPnextafter(tmp, SCIP_REAL_MAX);
      }

      return;
   }

   if( operand.inf <= 0.0 )
   {
      resultant->inf = -infinity;
   }
   else if( operand.inf == 1.0 )
   {
      resultant->inf = 0.0;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      resultant->inf = SCIPnextafter(log(operand.inf), SCIP_REAL_MIN);
   }

   if( operand.sup >= infinity )
   {
      resultant->sup =  infinity;
   }
   else if( operand.sup == 1.0 )
   {
      resultant->sup = 0.0;
   }
   else
   {
      assert(SCIPintervalGetRoundingMode() == SCIP_ROUND_NEAREST);
      resultant->sup = SCIPnextafter(log(operand.sup), SCIP_REAL_MAX);
   }
}

/** stores minimum of operands in resultant */
void SCIPintervalMin(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   resultant->inf = MIN(operand1.inf, operand2.inf);
   resultant->sup = MIN(operand1.sup, operand2.sup);
}

/** stores maximum of operands in resultant */
void SCIPintervalMax(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand1,           /**< first operand of operation */
   SCIP_INTERVAL         operand2            /**< second operand of operation */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand1));
   assert(!SCIPintervalIsEmpty(infinity, operand2));

   resultant->inf = MAX(operand1.inf, operand2.inf);
   resultant->sup = MAX(operand1.sup, operand2.sup);
}

/** stores absolute value of operand in resultant */
void SCIPintervalAbs(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   if( operand.inf <= 0.0 && operand.sup >= 0.0)
   {
      resultant->inf = 0.0;
      resultant->sup = MAX(-operand.inf, operand.sup);
   }
   else if( operand.inf > 0.0 )
   {
      *resultant = operand;
   }
   else
   {
      resultant->inf = -operand.sup;
      resultant->sup = -operand.inf;
   }
}

/** stores sine value of operand in resultant
 * NOTE: the operations are not applied rounding-safe here
 */
void SCIPintervalSin(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   SCIP_Real intervallen;
   SCIP_Real modinf;
   SCIP_Real modsup;
   SCIP_Real finf;
   SCIP_Real fsup;
   int a;
   int b;
   int nbetween;
   /* first one is 1 so even indices are the maximum points */
   static SCIP_Real extremepoints[] = {0.5*M_PI, 1.5*M_PI, 2.5*M_PI, 3.5*M_PI};

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   intervallen = operand.sup - operand.inf;
   if( intervallen >= 2*M_PI )
   {
      SCIPintervalSetBounds(resultant, -1.0, 1.0);
      return;
   }

   modinf = fmod(operand.inf, 2*M_PI);
   if( modinf < 0.0 )
      modinf += 2*M_PI;
   modsup = modinf + intervallen;

   for( b = 0; ; ++b )
   {
      if( modinf <= extremepoints[b] )
      {
         a = b;
         break;
      }
   }
   for( ; b < 4; ++b )
   {
      if( modsup <= extremepoints[b] )
         break;
   }

   nbetween = b-a;

   if( nbetween > 1 )
   {
      SCIPintervalSetBounds(resultant, -1.0, 1.0);
      return;
   }

   finf = sin(operand.inf);
   fsup = sin(operand.sup);

   if( nbetween == 0 )
   {
      if( a & 1 ) /* next extremepoint is minimum -> decreasing -> finf < fsup */
         SCIPintervalSetBounds(resultant, fsup, finf);
      else
         SCIPintervalSetBounds(resultant, finf, fsup);
   }
   else /* 1 extremepoint in between */
   {
       if( a & 1 ) /* extremepoint is minimum */
          SCIPintervalSetBounds(resultant, -1.0, MAX(finf,fsup));
       else
          SCIPintervalSetBounds(resultant, MIN(finf,fsup), 1.0);
   }

   /* above operations did not take outward rounding into account,
    * so extend the computed interval slightly to increase the chance that it will contain the complete sin(operand)
    */
   if( resultant->inf > -1.0 )
      resultant->inf = MAX(-1.0, resultant->inf - 1e-10 * REALABS(resultant->inf));  /*lint !e666*/
   if( resultant->sup <  1.0 )
      resultant->sup = MIN( 1.0, resultant->sup + 1e-10 * REALABS(resultant->sup));  /*lint !e666*/

   assert(resultant->inf <= resultant->sup);
}

/** stores cosine value of operand in resultant
 * NOTE: the operations are not applied rounding-safe here
 */
void SCIPintervalCos(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   SCIP_Real intervallen;
   SCIP_Real modinf;
   SCIP_Real modsup;
   SCIP_Real finf;
   SCIP_Real fsup;
   int a;
   int b;
   int nbetween;
   /* first one is -1 so even indices are the minimum points */
   static SCIP_Real extremepoints[] = {M_PI, 2*M_PI, 3*M_PI};

   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   intervallen = operand.sup - operand.inf;
   if( intervallen >= 2*M_PI )
   {
      SCIPintervalSetBounds(resultant, -1.0, 1.0);
      return;
   }

   modinf = fmod(operand.inf, 2*M_PI);
   if( modinf < 0.0 )
      modinf += 2*M_PI;
   modsup = modinf + intervallen;

   for( b = 0; ; ++b )
   {
      if( modinf <= extremepoints[b] )
      {
         a = b;
         break;
      }
   }
   for( ; b < 3; ++b )
   {
      if( modsup <= extremepoints[b] )
         break;
   }

   nbetween = b-a;

   if( nbetween > 1 )
   {
      SCIPintervalSetBounds(resultant, -1.0, 1.0);
      return;
   }

   finf = cos(operand.inf);
   fsup = cos(operand.sup);

   if( nbetween == 0 )
   {
      if( a & 1 ) /* next extremepoint is maximum -> increasing -> finf < fsup */
         SCIPintervalSetBounds(resultant, finf, fsup);
      else
         SCIPintervalSetBounds(resultant, fsup, finf);
   }
   else /* 1 extremepoint in between */
   {
       if( a & 1 ) /* extremepoint is maximum */
          SCIPintervalSetBounds(resultant, MIN(finf,fsup), 1.0);
       else
          SCIPintervalSetBounds(resultant, -1.0, MAX(finf,fsup));
   }

   /* above operations did not take outward rounding into account,
    * so extend the computed interval slightly to increase the chance that it will contain the complete cos(operand)
    */
   if( resultant->inf > -1.0 )
      resultant->inf = MAX(-1.0, resultant->inf - 1e-10 * REALABS(resultant->inf));  /*lint !e666*/
   if( resultant->sup <  1.0 )
      resultant->sup = MIN( 1.0, resultant->sup + 1e-10 * REALABS(resultant->sup));  /*lint !e666*/

   assert(resultant->inf <= resultant->sup);
}

/** stores sign of operand in resultant */
void SCIPintervalSign(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         operand             /**< operand of operation */
   )
{
   assert(resultant != NULL);
   assert(!SCIPintervalIsEmpty(infinity, operand));

   if( operand.sup < 0.0 )
   {
      resultant->inf = -1.0;
      resultant->sup = -1.0;
   }
   else if( operand.inf >= 0.0 )
   {
      resultant->inf =  1.0;
      resultant->sup =  1.0;
   }
   else
   {
      resultant->inf = -1.0;
      resultant->sup =  1.0;
   }
}

/** computes exact upper bound on \f$ a x^2 + b x \f$ for x in [xlb, xub], b an interval, and a scalar
 * 
 * Uses Algorithm 2.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008) */
SCIP_Real SCIPintervalQuadUpperBound(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_Real             a,                  /**< coefficient of x^2 */
   SCIP_INTERVAL         b_,                 /**< coefficient of x */
   SCIP_INTERVAL         x                   /**< range of x */
   )
{
   SCIP_Real b;
   SCIP_Real u;

   assert(!SCIPintervalIsEmpty(infinity, x));
   assert(b_.inf <  infinity);
   assert(b_.sup > -infinity);
   assert( x.inf <  infinity);
   assert( x.sup > -infinity);

   /* handle b*x separately */
   if( a == 0.0 )
   {
      if( (b_.inf <= -infinity && x.inf <   0.0     ) ||
         ( b_.inf <   0.0      && x.inf <= -infinity) ||
         ( b_.sup >   0.0      && x.sup >=  infinity) ||
         ( b_.sup >=  infinity && x.sup >   0.0     ) )
      {
         u = infinity;
      }
      else
      {
         SCIP_ROUNDMODE roundmode;
         SCIP_Real cand1, cand2, cand3, cand4;

         roundmode = SCIPintervalGetRoundingMode();
         SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
         cand1 = b_.inf * x.inf;
         cand2 = b_.inf * x.sup;
         cand3 = b_.sup * x.inf;
         cand4 = b_.sup * x.sup;
         u = MAX(MAX(cand1, cand2), MAX(cand3, cand4));

         SCIPintervalSetRoundingMode(roundmode);
      }

      return u;
   }

   if( x.sup <= 0.0 )
   { /* change sign of x: enclose a*x^2 + [-bub, -blb]*(-x) for (-x) in [-xub, -xlb] */
      u = x.sup;
      x.sup = -x.inf;
      x.inf = -u;
      b = -b_.inf;
   }
   else
   {
      b = b_.sup;
   }

   if( x.inf >= 0.0 )
   {  /* upper bound for a*x^2 + b*x */
      SCIP_ROUNDMODE roundmode;
      SCIP_Real s,t;

      if( b >= infinity )
         return infinity;

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);

      u = MAX(x.inf * (a*x.inf + b), x.sup * (a*x.sup + b));
      s = b/2;
      t = s/negate(a);
      if( t > x.inf && negate(2*a)*x.sup > b && s*t > u )
         u = s*t;

      SCIPintervalSetRoundingMode(roundmode);
      return u;
   }
   else
   {
      SCIP_INTERVAL xlow = x;
      SCIP_Real cand1;
      SCIP_Real cand2;
      assert(x.inf < 0.0 && x.sup > 0);

      xlow.sup = 0;  /* so xlow is lower part of interval */ 
      x.inf = 0;     /* so x    is upper part of interval now */
      cand1 = SCIPintervalQuadUpperBound(infinity, a, b_, xlow);
      cand2 = SCIPintervalQuadUpperBound(infinity, a, b_, x);
      return MAX(cand1, cand2);
   }
}

/** stores range of quadratic term in resultant
 * 
 * given scalar a and intervals b and x, computes interval for \f$ a x^2 + b x \f$ */
void SCIPintervalQuad(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         xrng                /**< range of x */
   )
{
   SCIP_Real tmp;

   if( SCIPintervalIsEmpty(infinity, xrng) )
   {
      SCIPintervalSetEmpty(resultant);
      return;
   }
   if( sqrcoeff == 0.0 )
   {
      SCIPintervalMul(infinity, resultant, lincoeff, xrng);
      return;
   }

   resultant->sup =  SCIPintervalQuadUpperBound(infinity,  sqrcoeff, lincoeff, xrng);

   tmp = lincoeff.inf;
   lincoeff.inf = -lincoeff.sup;
   lincoeff.sup = -tmp;
   resultant->inf = -SCIPintervalQuadUpperBound(infinity, -sqrcoeff, lincoeff, xrng);

   assert(resultant->sup >= resultant->inf);
}

/** computes interval with positive solutions of a quadratic equation with interval coefficients
 * 
 * Given intervals a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \in c\f$.
 */
void SCIPintervalSolveUnivariateQuadExpressionPositive(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs                 /**< right hand side of equation */
   )
{
   assert(resultant != NULL);

   /* find x>=0 s.t. a.inf x^2 + b.inf x <= c.sup  -> -a.inf x^2 - b.inf x >= -c.sup */
   if( lincoeff.inf <= -infinity || rhs.sup >= infinity || sqrcoeff.inf <= -infinity )
   {
      resultant->inf = 0.0;
      resultant->sup = infinity;
   }
   else
   {
      SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(infinity, resultant, -sqrcoeff.inf, -lincoeff.inf, -rhs.sup);
      SCIPdebugMessage("solve %g*x^2 + %g*x >= %g gives [%.20f, %.20f]\n", -sqrcoeff.inf, -lincoeff.inf, -rhs.sup, resultant->inf, resultant->sup);
   }

   /* find x>=0 s.t. a.sup x^2 + b.sup x >= c.inf */
   if( lincoeff.sup <  infinity && rhs.inf >  -infinity && sqrcoeff.sup <  infinity )
   {
      SCIP_INTERVAL res2;
      SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(infinity, &res2, sqrcoeff.sup, lincoeff.sup, rhs.inf);
      SCIPdebugMessage("solve %g*x^2 + %g*x >= %g gives [%.20f, %.20f]\n", sqrcoeff.sup, lincoeff.sup, rhs.inf, res2.inf, res2.sup);
      SCIPdebugMessage("intersection of [%.20f, %.20f] and [%.20f, %.20f]", resultant->inf, resultant->sup, res2.inf, res2.sup);
      /* intersect both results */
      SCIPintervalIntersect(resultant, *resultant, res2);
      SCIPdebugPrintf(" gives [%.20f, %.20f]\n", resultant->inf, resultant->sup);
   }
   /* else res2 = [0, infty] */

   if( resultant->inf >= infinity || resultant->sup <= -infinity )
   {
      SCIPintervalSetEmpty(resultant);
   }
}

/** computes positive solutions of a quadratic equation with scalar coefficients
 * 
 * Given scalar a, b, and c, this function computes an interval that contains all positive solutions of \f$ a x^2 + b x \geq c\f$.
 * Implements Algorithm 3.2 from Domes and Neumaier: Constraint propagation on quadratic constraints (2008).
 */
void SCIPintervalSolveUnivariateQuadExpressionPositiveAllScalar(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_Real             sqrcoeff,           /**< coefficient of x^2 */
   SCIP_Real             lincoeff,           /**< coefficient of x */
   SCIP_Real             rhs                 /**< right hand side of equation */
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_Real     b;
   SCIP_Real     delta;
   SCIP_Real     z;

   assert(resultant != NULL);
   assert(sqrcoeff <  infinity);
   assert(sqrcoeff > -infinity);

   resultant->inf = 0.0;
   resultant->sup = infinity;

   roundmode = SCIPintervalGetRoundingMode();

   SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
   b = lincoeff / 2.0;

   if( lincoeff >= 0.0 )
   { /* b >= 0.0 */
      /*SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS); */
      if( rhs > 0.0 )
      { /* b >= 0.0 and c > 0.0 */
         delta = b*b + sqrcoeff*rhs;
         if( delta < 0.0 || (sqrcoeff == 0.0 && lincoeff == 0.0) )
         {
            SCIPintervalSetEmpty(resultant);
         }
         else
         {
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = SCIPnextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
            z += b;
            resultant->inf = negate(negate(rhs)/z);
            if( sqrcoeff < 0.0 )
               resultant->sup = z / negate(sqrcoeff);
         }
      }
      else
      { /* b >= 0.0 and c <= 0.0 */
         if( sqrcoeff < 0.0 )
         {
            delta = b*b + sqrcoeff*rhs;
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = SCIPnextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
            z += b;
            resultant->sup = z / negate(sqrcoeff);
         }
      }
   }
   else
   { /* b < 0.0 */
      if( rhs > 0.0 )
      { /* b < 0.0 and c > 0.0 */
         if( sqrcoeff > 0.0 )
         {
            SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
            delta = b*b + sqrcoeff*rhs;
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = SCIPnextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS);
            z += negate(b);
            resultant->inf = z / sqrcoeff;
         }
         else
         {
            SCIPintervalSetEmpty(resultant);
         }
      }
      else
      { /* b < 0.0 and c <= 0.0 */
         /* SCIPintervalSetRoundingMode(SCIP_ROUND_DOWNWARDS); */
         delta = b*b + sqrcoeff * rhs;
         if( delta >= 0.0 && sqrcoeff <= 0.0 )
         {
            SCIPintervalSetRoundingMode(SCIP_ROUND_NEAREST);
            z = SCIPnextafter(sqrt(delta), SCIP_REAL_MAX);
            SCIPintervalSetRoundingMode(SCIP_ROUND_UPWARDS);
            z += negate(b);
            resultant->sup = negate(rhs/z);
         }
         /* @todo actually we could generate a hole here */
      }
   }

   SCIPintervalSetRoundingMode(roundmode);
}

/** solves a quadratic equation with interval coefficients
 *
 * Given intervals a, b and c, this function computes an interval that contains all solutions of \f$ a x^2 + b x \in c\f$
 */
void SCIPintervalSolveUnivariateQuadExpression(
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        resultant,          /**< resultant interval of operation */
   SCIP_INTERVAL         sqrcoeff,           /**< coefficient of x^2 */
   SCIP_INTERVAL         lincoeff,           /**< coefficient of x */
   SCIP_INTERVAL         rhs                 /**< right hand side of equation */
   )
{
   SCIP_Real tmp;
   SCIP_INTERVAL xpos;
   SCIP_INTERVAL xneg;

   assert(resultant != NULL);

   if( sqrcoeff.inf == 0.0 && sqrcoeff.sup == 0.0 )
   { /* relatively easy case: x \in rhs / lincoeff */
      if( lincoeff.inf == 0.0 && lincoeff.sup == 0.0 )
      { /* equation became 0.0 \in rhs */
         if( rhs.inf <= 0.0 && rhs.sup >= 0.0 )
            SCIPintervalSetEntire(infinity, resultant);
         else
            SCIPintervalSetEmpty(resultant);
      }
      else
         SCIPintervalDiv(infinity, resultant, rhs, lincoeff);
      SCIPdebugMessage("  solving [%g,%g]*x in [%g,%g] gives [%g,%g]\n", SCIPintervalGetInf(lincoeff), SCIPintervalGetSup(lincoeff), SCIPintervalGetInf(rhs), SCIPintervalGetSup(rhs), SCIPintervalGetInf(*resultant), SCIPintervalGetSup(*resultant));
      return;
   }

   if( lincoeff.inf == 0.0 && lincoeff.sup == 0.0 )
   { /* easy case: x \in +/- sqrt(rhs/a) */
      SCIPintervalDiv(infinity, resultant, rhs, sqrcoeff);
      SCIPintervalSquareRoot(infinity, resultant, *resultant);
      resultant->inf = -resultant->sup;
      return;
   }

   /* find all x>=0 such that a*x^2+b*x = c */
   SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &xpos, sqrcoeff, lincoeff, rhs);
   SCIPdebugMessage("  positive solutions of [%g,%g]*x^2 + [%g,%g]*x in [%g,%g] are [%g,%g]\n",
      sqrcoeff.inf, sqrcoeff.sup, lincoeff.inf, lincoeff.sup, rhs.inf, rhs.sup, xpos.inf, xpos.sup);

   tmp = lincoeff.inf;
   lincoeff.inf = -lincoeff.sup;
   lincoeff.sup = -tmp;

   /* find all x>=0 such that a*x^2-b*x = c */
   SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &xneg, sqrcoeff, lincoeff, rhs);
   tmp = xneg.inf;
   xneg.inf = -xneg.sup;
   xneg.sup = -tmp;
   SCIPdebugMessage("  negative solutions of [%g,%g]*x^2 + [%g,%g]*x in [%g,%g] are [%g,%g]\n",
      sqrcoeff.inf, sqrcoeff.sup, lincoeff.inf, lincoeff.sup, rhs.inf, rhs.sup, xneg.inf, xneg.sup);

   SCIPintervalUnify(resultant, xpos, xneg);
   SCIPdebugMessage("  unify gives [%g,%g]\n", SCIPintervalGetInf(*resultant), SCIPintervalGetSup(*resultant));
}

/** stores range of bivariate quadratic term in resultant
 * given scalars ax, ay, axy, bx, and by and intervals for x and y, computes interval for \f$ ax x^2 + ay y^2 + axy x y + bx x + by y \f$
 * NOTE: the operations are not applied rounding-safe here
 */
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
   )
{
   /* we use double double precision and finally widen the computed range by 1e-8% to compensate for not computing rounding-safe here */
   SCIP_Real minval;
   SCIP_Real maxval;
   SCIP_Real val;
   SCIP_Real x;
   SCIP_Real y;
   SCIP_Real denom;

   assert(resultant != NULL);
   assert(xbnds.inf <= xbnds.sup);
   assert(ybnds.inf <= ybnds.sup);

   /* if we are separable, then fall back to use SCIPintervalQuad two times and add */
   if( axy == 0.0 )
   {
      SCIP_INTERVAL tmp;

      SCIPintervalSet(&tmp, bx);
      SCIPintervalQuad(infinity, resultant, ax, tmp, xbnds);

      SCIPintervalSet(&tmp, by);
      SCIPintervalQuad(infinity, &tmp, ay, tmp, ybnds);

      SCIPintervalAdd(infinity, resultant, *resultant, tmp);

      return;
   }

   SCIPintervalSet(resultant, 0.0);

   minval =  infinity;
   maxval = -infinity;

   /* check minima/maxima of expression */
   denom = 4.0 * ax * ay - axy * axy;
   if( REALABS(denom) > 1e-9 )
   {
      x = (axy * by - 2.0 * ay * bx) / denom;
      y = (axy * bx - 2.0 * ax * by) / denom;
      if( xbnds.inf <= x && x <= xbnds.sup && ybnds.inf <= y && y <= ybnds.sup )
      {
         val = (axy * bx * by - ay * bx * bx - ax * by * by) / denom;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);
      }
   }
   else if( REALABS(2.0 * ay * bx - axy * by) <= 1e-9 )
   {
      /* The whole line (x, -bx/axy - (axy/2ay) x) defines an extreme point with value -ay bx^2 / axy^2
       * If x is unbounded, then there is an (x,y) with y in ybnds where the extreme value is assumed.
       * If x is bounded on at least one side, then we can rely that the checks below for x at one of its bounds will check this extreme point.
       */
      if( xbnds.inf <= -infinity && xbnds.sup >= infinity )
      {
         val = -ay * bx * bx / (axy * axy);
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);
      }
   }

   /* check boundary of box xbnds x ybnds */

   if( xbnds.inf <= -infinity )
   {
      /* check value for x -> -infinity */
      if( ax > 0.0 )
         maxval =  infinity;
      else if( ax < 0.0 )
         minval = -infinity;
      else if( ax == 0.0 )
      {
         /* bivar(x,y) tends to -(bx+axy y) * infinity */

         if( ybnds.inf <= -infinity )
            val = (axy < 0.0 ? -infinity : infinity);
         else if( bx + axy * ybnds.inf < 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);

         if( ybnds.sup >= infinity )
            val = (axy < 0.0 ? infinity : -infinity);
         else if( bx + axy * ybnds.sup < 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);
      }
   }
   else
   {
      /* get range of bivar(xbnds.inf, y) for y in ybnds */
      SCIP_INTERVAL tmp;
      SCIP_INTERVAL ycoef;

      SCIPintervalSet(&ycoef, axy * xbnds.inf + by);
      SCIPintervalQuad(infinity, &tmp, ay, ycoef, ybnds);
      SCIPintervalAddScalar(infinity, &tmp, tmp, (SCIP_Real)(ax * xbnds.inf * xbnds.inf + bx * xbnds.inf));
      minval = MIN(tmp.inf, minval);
      maxval = MAX(tmp.sup, maxval);
   }

   if( xbnds.sup >= infinity )
   {
      /* check value for x -> infinity */
      if( ax > 0.0 )
         maxval =  infinity;
      else if( ax < 0.0 )
         minval = -infinity;
      else if( ax == 0.0 )
      {
         /* bivar(x,y) tends to (bx+axy y) * infinity */

         if( ybnds.inf <= -infinity )
            val = (axy > 0.0 ? -infinity : infinity);
         else if( bx + axy * ybnds.inf > 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);

         if( ybnds.sup >= infinity )
            val = (axy > 0.0 ? infinity : -infinity);
         else if( bx + axy * ybnds.sup > 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);
      }
   }
   else
   {
      /* get range of bivar(xbnds.sup, y) for y in ybnds */
      SCIP_INTERVAL tmp;
      SCIP_INTERVAL ycoef;

      SCIPintervalSet(&ycoef, axy * xbnds.sup + by);
      SCIPintervalQuad(infinity, &tmp, ay, ycoef, ybnds);
      SCIPintervalAddScalar(infinity, &tmp, tmp, (SCIP_Real)(ax * xbnds.sup * xbnds.sup + bx * xbnds.sup));
      minval = MIN(tmp.inf, minval);
      maxval = MAX(tmp.sup, maxval);
   }

   if( ybnds.inf <= -infinity )
   {
      /* check value for y -> -infinity */
      if( ay > 0.0 )
         maxval =  infinity;
      else if( ay < 0.0 )
         minval = -infinity;
      else if( ay == 0.0 )
      {
         /* bivar(x,y) tends to -(by+axy x) * infinity */

         if( xbnds.inf <= -infinity )
            val = (axy < 0.0 ? -infinity : infinity);
         else if( by + axy * xbnds.inf < 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);

         if( xbnds.sup >= infinity )
            val = (axy < 0.0 ? infinity : -infinity);
         else if( by + axy * xbnds.sup < 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);
      }
   }
   else
   {
      /* get range of bivar(x, ybnds.inf) for x in xbnds */
      SCIP_INTERVAL tmp;
      SCIP_INTERVAL xcoef;

      SCIPintervalSet(&xcoef, axy * ybnds.inf + bx);
      SCIPintervalQuad(infinity, &tmp, ax, xcoef, xbnds);
      SCIPintervalAddScalar(infinity, &tmp, tmp, ay * ybnds.inf * ybnds.inf + by * ybnds.inf);
      minval = MIN(tmp.inf, minval);
      maxval = MAX(tmp.sup, maxval);
   }

   if( ybnds.sup >= infinity )
   {
      /* check value for y -> infinity */
      if( ay > 0.0 )
         maxval =  infinity;
      else if( ay < 0.0 )
         minval = -infinity;
      else if( ay == 0.0 )
      {
         /* bivar(x,y) tends to (by+axy x) * infinity */

         if( xbnds.inf <= -infinity )
            val = (axy > 0.0 ? -infinity : infinity);
         else if( by + axy * xbnds.inf > 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);

         if( xbnds.sup >= infinity )
            val = (axy > 0.0 ? infinity : -infinity);
         else if( by + axy * xbnds.sup > 0.0 )
            val = infinity;
         else
            val = -infinity;
         minval = MIN(val, minval);
         maxval = MAX(val, maxval);
      }
   }
   else
   {
      /* get range of bivar(x, ybnds.sup) for x in xbnds */
      SCIP_INTERVAL tmp;
      SCIP_INTERVAL xcoef;

      SCIPintervalSet(&xcoef, axy * ybnds.sup + bx);
      SCIPintervalQuad(infinity, &tmp, ax, xcoef, xbnds);
      SCIPintervalAddScalar(infinity, &tmp, tmp, (SCIP_Real)(ay * ybnds.sup * ybnds.sup + by * ybnds.sup));
      minval = MIN(tmp.inf, minval);
      maxval = MAX(tmp.sup, maxval);
   }

   minval -= 1e-10 * REALABS(minval);
   maxval += 1e-10 * REALABS(maxval);
   SCIPintervalSetBounds(resultant, (SCIP_Real)minval, (SCIP_Real)maxval);

   SCIPdebugMessage("range for %gx^2 + %gy^2 + %gxy + %gx + %gy = [%g, %g] for x = [%g, %g], y=[%g, %g]\n",
      ax, ay, axy, bx, by, minval, maxval, xbnds.inf, xbnds.sup, ybnds.inf, ybnds.sup);
}

/** solves a bivariate quadratic equation for the first variable
 * given scalars ax, ay, axy, bx and by, and intervals for x, y, and rhs,
 * computes \f$ \{ x \in \mathbf{x} : \exists y \in \mathbf{y} : a_x x^2 + a_y y^2 + a_{xy} x y + b_x x + b_y y \in \mathbf{\mbox{rhs}} \} \f$
 * NOTE: the operations are not applied rounding-safe here
 */
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
   )
{
   /* we use double double precision and finally widen the computed range by 1e-8% to compensate for not computing rounding-safe here */
   SCIP_Real val;

   assert(resultant != NULL);

   if( axy == 0.0 )
   {
      /* if axy == 0, fall back to SCIPintervalSolveUnivariateQuadExpression */
      SCIP_INTERVAL ytermrng;
      SCIP_INTERVAL sqrcoef;
      SCIP_INTERVAL lincoef;
      SCIP_INTERVAL pos;
      SCIP_INTERVAL neg;

      SCIPintervalSet(&lincoef, by);
      SCIPintervalQuad(infinity, &ytermrng, ay, lincoef, ybnds);
      SCIPintervalSub(infinity, &rhs, rhs, ytermrng);

      SCIPintervalSet(&sqrcoef, ax);

      /* get positive solutions, if of interest */
      if( xbnds.sup >= 0.0 )
      {
         SCIPintervalSet(&lincoef, bx);
         SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &pos, sqrcoef, lincoef, rhs);
         SCIPintervalIntersect(&pos, pos, xbnds);
      }
      else
         SCIPintervalSetEmpty(&pos);

      /* get negative solutions, if of interest */
      if( xbnds.inf < 0.0 )
      {
         SCIPintervalSet(&lincoef, -bx);
         SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &neg, sqrcoef, lincoef, rhs);
         if( !SCIPintervalIsEmpty(infinity, neg) )
         {
            SCIPintervalSetBounds(&neg, -neg.sup, -neg.inf);
            SCIPintervalIntersect(&neg, neg, xbnds);
         }
      }
      else
         SCIPintervalSetEmpty(&neg);

      SCIPintervalUnify(resultant, pos, neg);

      return;
   }

   if( ax < 0.0 )
   {
      SCIP_Real tmp;
      tmp = rhs.inf;
      rhs.inf = -rhs.sup;
      rhs.sup = -tmp;

      SCIPintervalSolveBivariateQuadExpressionAllScalar(infinity, resultant, -ax, -ay, -axy, -bx, -by, rhs, xbnds, ybnds);
      return;
   }
   assert(ax >= 0.0);

   *resultant = xbnds;

   if( ax > 0.0 )
   {
      SCIP_Real sqrtax;
      SCIP_Real minvalleft;
      SCIP_Real maxvalleft;
      SCIP_Real minvalright;
      SCIP_Real maxvalright;
      SCIP_Real ymin;
      SCIP_Real rcoef_y;
      SCIP_Real rcoef_yy;
      SCIP_Real rcoef_const;

      sqrtax = sqrt(ax);

      /* rewrite equation as (sqrt(ax)x + b(y))^2 \in r(rhs,y), where
       * b(y) = (bx + axy y)/(2sqrt(ax)), r(rhs,y) = rhs - ay y^2 - by y + b(y)^2
       *
       * -> r(rhs,y) = bx^2/(4ax) + rhs + (axy bx/(2ax) - by)*y + (axy^2/(4ax) - ay)*y^2
       */
      rcoef_y     = axy * bx  / (2.0*ax) - by;
      rcoef_yy    = axy * axy / (4.0*ax) - ay;
      rcoef_const = bx  * bx  / (4.0*ax);

#define CALCB(y)    ((bx + axy * (y)) / (2.0 * sqrtax))
#define CALCR(c,y)  (rcoef_const + (c) + (rcoef_y + rcoef_yy * (y)) * (y))

      /* check whether r(rhs,y) is always negative */
      if( rhs.sup < infinity )
      {
         SCIP_INTERVAL ycoef;
         SCIP_Real ub;

         SCIPintervalSet(&ycoef, (SCIP_Real)rcoef_y);
         ub = (SCIP_Real)(SCIPintervalQuadUpperBound(infinity, (SCIP_Real)rcoef_yy, ycoef, ybnds) + rhs.sup + rcoef_const);

         if( EPSN(ub, 1e-9) )
         {
            SCIPintervalSetEmpty(resultant);
            return;
         }
         else if( ub < 0.0 )
         {
            /* it looks like there will be no solution (rhs < 0), but we are very close and above operations did not take care of careful rounding
             * thus, we relax rhs a be feasible a bit (-ub would be sufficient, but that would put us exactly onto the boundary)
             * see also #1861
             */
            rhs.sup += -2.0*ub;
         }
      }

      /* we have sqrt(ax)x \in (-sqrt(r(rhs,y))-b(y)) \cup (sqrt(r(rhs,y))-b(y))
       * compute minima and maxima of both functions such that
       *
       * [minvalleft,  maxvalleft ] = -sqrt(r(rhs,y))-b(y)
       * [minvalright, maxvalright] =  sqrt(r(rhs,y))-b(y)
       */

      minvalleft  =  infinity;
      maxvalleft  = -infinity;
      minvalright =  infinity;
      maxvalright = -infinity;

      if( rhs.sup >= infinity )
      {
         /* we can't do much if rhs.sup is infinite
          * but we may do a bit of xbnds isn't too huge and rhs.inf > -infinity
          */
         minvalleft  = -infinity;
         maxvalright =  infinity;
      }

      /* evaluate at lower bound of y, as long as r(rhs,ylb) > 0 */
      if( ybnds.inf <= -infinity )
      {
         /* check limit of +/-sqrt(r(rhs,y))-b(y) for y -> -infty */
         if( !EPSZ(ay, 1e-9) && axy * axy >= 4.0 * ax * ay )
         {
            if( axy < 0.0 )
            {
               minvalleft = -infinity;

               if( ay > 0.0 )
                  minvalright = -infinity;
               else
                  maxvalright = infinity;
            }
            else
            {
               maxvalright = infinity;

               if( ay > 0.0 )
                  maxvalleft = infinity;
               else
                  minvalleft = -infinity;
            }
         }
         else if( !EPSZ(ay, 1e-9) )
         {
            /* here axy * axy < 4 * ax * ay, so need to check for zeros of r(rhs,y), which is done below */
         }
         else
         {
            /* here ay = 0.0, which gives a limit of -by/2 for -sqrt(r(rhs,y))-b(y) */
            minvalleft = -by / 2.0;
            maxvalleft = -by / 2.0;
            /* here ay = 0.0, which gives a limit of +infinity for sqrt(r(rhs,y))-b(y) */
            maxvalright = infinity;
         }
      }
      else
      {
         SCIP_Real b;
         SCIP_Real c;

         b = CALCB(ybnds.inf);

         if( rhs.sup <  infinity )
         {
            c = CALCR(rhs.sup, ybnds.inf);

            if( c > 0.0 )
            {
               SCIP_Real sqrtc;

               sqrtc = sqrt(c);
               minvalleft  = MIN(-sqrtc - b, minvalleft);
               maxvalright = MAX( sqrtc - b, maxvalright);
            }
         }

         if( rhs.inf > -infinity )
         {
            c = CALCR(rhs.inf, ybnds.inf);

            if( c > 0.0 )
            {
               SCIP_Real sqrtc;

               sqrtc = sqrt(c);
               maxvalleft  = MAX(-sqrtc - b, maxvalleft);
               minvalright = MIN( sqrtc - b, minvalright);
            }
         }
      }

      /* evaluate at upper bound of y, as long as r(rhs, yub) > 0 */
      if( ybnds.sup >= infinity )
      {
         /* check limit of +/-sqrt(r(rhs,y))-b(y) for y -> +infty */
         if( !EPSZ(ay, 1e-9) && axy * axy >= 4.0 * ax * ay )
         {
            if( axy > 0.0 )
            {
               minvalleft = -infinity;

               if( ay > 0.0 )
                  minvalright = -infinity;
               else
                  maxvalright = infinity;
            }
            else
            {
               maxvalright = infinity;

               if( ay > 0.0 )
                  maxvalleft = infinity;
               else
                  minvalleft = -infinity;
            }
         }
         else if( !EPSZ(ay, 1e-9) )
         {
            /* here axy * axy < 4 * ax * ay, so need to check for zeros of r(rhs,y), which will happen below */
         }
         else
         {
            /* here ay = 0.0, which gives a limit of -infinity for -sqrt(r(rhs,y))-b(y) */
            minvalleft = -infinity;
            /* here ay = 0.0, which gives a limit of -by/2 for sqrt(r(rhs,y))-b(y) */
            minvalright = MIN(minvalright, -by / 2.0);
            maxvalright = MAX(maxvalright, -by / 2.0);
         }
      }
      else
      {
         SCIP_Real b;
         SCIP_Real c;

         b = CALCB(ybnds.sup);

         if( rhs.sup <  infinity )
         {
            c = CALCR(rhs.sup, ybnds.sup);

            if( c > 0.0 )
            {
               SCIP_Real sqrtc;

               sqrtc = sqrt(c);
               minvalleft  = MIN(-sqrtc - b, minvalleft);
               maxvalright = MAX( sqrtc - b, maxvalright);
            }
         }

         if( rhs.inf > -infinity )
         {
            c = CALCR(rhs.inf, ybnds.sup);

            if( c > 0.0 )
            {
               SCIP_Real sqrtc;

               sqrtc = sqrt(c);
               maxvalleft  = MAX(-sqrtc - b, maxvalleft);
               minvalright = MIN( sqrtc - b, minvalright);
            }
         }
      }

      /* evaluate at ymin = y_{_,+}, if inside ybnds
       * if ay = 0 or 2ay*bx == axy*by, then there is no ymin */
      if( !EPSZ(ay, 1e-9) )
      {
         if( REALABS(axy*axy - 4.0*ax*ay) > 1e-9 )
         {
            SCIP_Real sqrtterm;

            if( rhs.sup < infinity )
            {
               sqrtterm = axy * axy * ay * (ay * bx * bx - axy * bx * by + ax * by * by - axy * axy * rhs.sup + 4.0 * ax * ay * rhs.sup);
               if( !EPSN(sqrtterm, 1e-9) )
               {
                  sqrtterm = sqrt(MAX(sqrtterm, 0.0));
                  /* check first candidate for extreme points of +/-sqrt(rhs(r,y))-b(y) */
                  ymin = axy * ay * bx - 2.0 * ax * ay * by - sqrtterm;
                  ymin /= ay;
                  ymin /= 4.0 * ax * ay - axy * axy;

                  if( ymin > ybnds.inf && ymin < ybnds.sup )
                  {
                     SCIP_Real b;
                     SCIP_Real c;

                     b = CALCB(ymin);
                     c = CALCR(rhs.sup, ymin);

                     if( c > 0.0 )
                     {
                        SCIP_Real sqrtc;

                        sqrtc = sqrt(c);
                        minvalleft  = MIN(-sqrtc - b, minvalleft);
                        maxvalright = MAX( sqrtc - b, maxvalright);
                     }
                  }

                  /* check second candidate for extreme points of +/-sqrt(rhs(r,y))-b(y) */
                  ymin = axy * ay * bx - 2.0 * ax * ay * by + sqrtterm;
                  ymin /= ay;
                  ymin /= 4.0 * ax * ay - axy * axy;

                  if( ymin > ybnds.inf && ymin < ybnds.sup )
                  {
                     SCIP_Real b;
                     SCIP_Real c;

                     b = CALCB(ymin);
                     c = CALCR(rhs.sup, ymin);

                     if( c > 0.0 )
                     {
                        SCIP_Real sqrtc;

                        sqrtc = sqrt(c);
                        minvalleft  = MIN(-sqrtc - b, minvalleft);
                        maxvalright = MAX( sqrtc - b, maxvalright);
                     }
                  }
               }
            }

            if( rhs.inf > -infinity )
            {
               sqrtterm = axy * axy * ay * (ay * bx * bx - axy * bx * by + ax * by * by - axy * axy * rhs.inf + 4.0 * ax * ay * rhs.inf);
               if( !EPSN(sqrtterm, 1e-9) )
               {
                  sqrtterm = sqrt(MAX(sqrtterm, 0.0));
                  /* check first candidate for extreme points of +/-sqrt(r(rhs,y))-b(y) */
                  ymin = axy * ay * bx - 2.0 * ax * ay * by - sqrtterm;
                  ymin /= ay;
                  ymin /= 4.0 * ax * ay - axy * axy;

                  if( ymin > ybnds.inf && ymin < ybnds.sup )
                  {
                     SCIP_Real b;
                     SCIP_Real c;

                     b = CALCB(ymin);
                     c = CALCR(rhs.inf, ymin);

                     if( c > 0.0 )
                     {
                        SCIP_Real sqrtc;

                        sqrtc = sqrt(c);
                        maxvalleft  = MAX(-sqrtc - b, maxvalleft);
                        minvalright = MIN( sqrtc - b, minvalright);
                     }
                  }

                  /* check second candidate for extreme points of +/-sqrt(c(y))-b(y) */
                  ymin = axy * ay * bx - 2.0 * ax * ay * by + sqrtterm;
                  ymin /= ay;
                  ymin /= 4.0 * ax * ay - axy * axy;

                  if( ymin > ybnds.inf && ymin < ybnds.sup )
                  {
                     SCIP_Real b;
                     SCIP_Real c;

                     b = CALCB(ymin);
                     c = CALCR(rhs.inf, ymin);

                     if( c > 0.0 )
                     {
                        SCIP_Real sqrtc;

                        sqrtc = sqrt(c);
                        maxvalleft  = MAX(-sqrtc - b, maxvalleft);
                        minvalright = MIN( sqrtc - b, minvalright);
                     }
                  }
               }
            }

         }
         else if( REALABS(2.0 * ay * bx - axy * by) > 1e-9 )
         {
            if( rhs.sup < infinity )
            {
               ymin = - (4.0 * ay * bx * by - axy * by * by + 4.0 * axy * ay * rhs.sup);
               ymin /= 4.0 * ay;
               ymin /= 2.0 * ay * bx - axy * by;

               if( ymin > ybnds.inf && ymin < ybnds.sup )
               {
                  SCIP_Real b;
                  SCIP_Real c;

                  b = CALCB(ymin);
                  c = CALCR(rhs.sup, ymin);

                  if( c > 0.0 )
                  {
                     SCIP_Real sqrtc;

                     sqrtc = sqrt(c);
                     minvalleft  = MIN(-sqrtc - b, minvalleft);
                     maxvalright = MAX( sqrtc - b, maxvalright);
                  }
               }
            }

            if( rhs.inf > -infinity )
            {
               ymin = - (4.0 * ay * bx * by - axy * by * by + 4.0 * axy * ay * rhs.inf);
               ymin /= 4.0 * ay;
               ymin /= 2.0 * ay * bx - axy * by;

               if( ymin > ybnds.inf && ymin < ybnds.sup )
               {
                  SCIP_Real b;
                  SCIP_Real c;

                  b = CALCB(ymin);
                  c = CALCR(rhs.inf, ymin);

                  if( c > 0.0 )
                  {
                     SCIP_Real sqrtc;

                     sqrtc = sqrt(c);
                     maxvalleft  = MAX(-sqrtc - b, maxvalleft);
                     minvalright = MIN( sqrtc - b, minvalright);
                  }
               }
            }
         }
      }

      /* evaluate the case r(rhs,y) = 0, which is to min/max -b(y) w.r.t. r(rhs,y) = 0, y in ybnds
       * with the above assignments
       *   rcoef_y     = axy * bx  / (2.0*ax) - by;
       *   rcoef_yy    = axy * axy / (4.0*ax) - ay;
       *   rcoef_const = bx  * bx  / (4.0*ax);
       * we have r(rhs,y) = rhs + rcoef_const + rcoef_y * y + rcoef_yy * y^2
       *
       * thus, r(rhs,y) = 0 <-> rcoef_y * y + rcoef_yy * y^2 in -rhs - rcoef_const
       *
       */
      {
         SCIP_INTERVAL rcoef_yy_int;
         SCIP_INTERVAL rcoef_y_int;
         SCIP_INTERVAL rhs2;
         SCIP_Real b;

         /* setup rcoef_yy, rcoef_y and -rhs-rcoef_const as intervals */
         SCIPintervalSet(&rcoef_yy_int, (SCIP_Real)rcoef_yy);
         SCIPintervalSet(&rcoef_y_int, (SCIP_Real)rcoef_y);
         SCIPintervalSetBounds(&rhs2, (SCIP_Real)(-rhs.sup - rcoef_const), (SCIP_Real)(-rhs.inf - rcoef_const));

         /* first find all y >= 0 such that rcoef_y * y + rcoef_yy * y^2 in -rhs2, if ybnds.sup >= 0.0
          * and evaluate -b(y) w.r.t. these values
          */
         if( ybnds.sup >= 0.0 )
         {
            SCIP_INTERVAL ypos;

            SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &ypos, rcoef_yy_int, rcoef_y_int, rhs2);
            SCIPintervalIntersect(&ypos, ypos, ybnds);
            if( !SCIPintervalIsEmpty(infinity, ypos) )
            {
               assert(ypos.inf >= 0.0); /* we computed only positive solutions above */
               b = CALCB(ypos.inf);
               minvalleft  = MIN(minvalleft, -b);
               maxvalleft  = MAX(maxvalleft, -b);
               minvalright = MIN(minvalright, -b);
               maxvalright = MAX(maxvalright, -b);

               if( ypos.sup < infinity )
               {
                  b = CALCB(ypos.sup);
                  minvalleft  = MIN(minvalleft, -b);
                  maxvalleft  = MAX(maxvalleft, -b);
                  minvalright = MIN(minvalright, -b);
                  maxvalright = MAX(maxvalright, -b);
               }
               else
               {
                  /* -b(y) = - (bx + axy * y) / (2.0 * sqrt(ax)) -> -sign(axy)*infinity for y -> infinity */
                  if( axy > 0.0 )
                  {
                     minvalleft  = -infinity;
                     minvalright = -infinity;
                  }
                  else
                  {
                     maxvalleft  =  infinity;
                     maxvalright =  infinity;
                  }
               }
            }
         }

         /* next find all y <= 0 such that rcoef_y * y + rcoef_yy * y^2 in -rhs2, if ybnds.inf < 0.0
          * and evaluate -b(y) w.r.t. these values
          * (the case y fixed to 0 has been handled in the ybnds.sup >= 0 case above)
          */
         if( ybnds.inf < 0.0 )
         {
            SCIP_INTERVAL yneg;

            /* find all y >= 0 such that -rcoef_y * y + rcoef_yy * y^2 in -rhs2 and then negate y */
            SCIPintervalSet(&rcoef_y_int, -(SCIP_Real)rcoef_y);
            SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, &yneg, rcoef_yy_int, rcoef_y_int, rhs2);
            if( !SCIPintervalIsEmpty(infinity, yneg) )
            {
               SCIPintervalSetBounds(&yneg, -yneg.sup, -yneg.inf);
               SCIPintervalIntersect(&yneg, yneg, ybnds);
            }
            if( !SCIPintervalIsEmpty(infinity, yneg) )
            {
               if( yneg.inf > -infinity )
               {
                  b = CALCB(yneg.inf);
                  minvalleft  = MIN(minvalleft,  -b);
                  maxvalleft  = MAX(maxvalleft,  -b);
                  minvalright = MIN(minvalright, -b);
                  maxvalright = MAX(maxvalright, -b);
               }
               else
               {
                  /* -b(y) = - (bx + axy * y) / (2.0 * sqrt(ax)) -> sign(axy)*infinity for y -> -infinity */
                  if( axy > 0.0 )
                  {
                     maxvalleft  =  infinity;
                     maxvalright =  infinity;
                  }
                  else
                  {
                     minvalleft  = -infinity;
                     minvalright = -infinity;
                  }
               }

               assert(yneg.sup <= 0.0); /* we computed only negative solutions above */
               b = CALCB(yneg.sup);
               minvalleft  = MIN(minvalleft,  -b);
               maxvalleft  = MAX(maxvalleft,  -b);
               minvalright = MIN(minvalright, -b);
               maxvalright = MAX(maxvalright, -b);
            }
         }
      }

      if( rhs.inf > -infinity && xbnds.inf > -infinity && EPSGT(xbnds.inf, maxvalleft / sqrtax, 1e-9) )
      {
         /* if sqrt(ax)*x > -sqrt(r(rhs,y))-b(y), then tighten lower bound of sqrt(ax)*x to lower bound of sqrt(r(rhs,y))-b(y)
          * this is only possible if rhs.inf > -infinity, otherwise the value for maxvalleft is not valid (but tightening wouldn't be possible for sure anyway) */
         assert(EPSGE(minvalright, minvalleft, 1e-9)); /* right interval should not be above lower bound of left interval */
         if( minvalright > -infinity )
         {
            assert(minvalright < infinity);
            resultant->inf = (SCIP_Real)(minvalright / sqrtax);
         }
      }
      else
      {
         /* otherwise, tighten lower bound of sqrt(ax)*x to lower bound of -sqrt(r(rhs,y))-b(y) */
         if( minvalleft > -infinity )
         {
            assert(minvalleft < infinity);
            resultant->inf = (SCIP_Real)(minvalleft / sqrtax);
         }
      }

      if( rhs.inf > -infinity && xbnds.sup < infinity && EPSLT(xbnds.sup, minvalright / sqrtax, 1e-9) )
      {
         /* if sqrt(ax)*x < sqrt(r(rhs,y))-b(y), then tighten upper bound of sqrt(ax)*x to upper bound of -sqrt(r(rhs,y))-b(y)
          * this is only possible if rhs.inf > -infinity, otherwise the value for minvalright is not valid (but tightening wouldn't be possible for sure anyway) */
         assert(EPSLE(maxvalleft, maxvalright, 1e-9)); /* left interval should not be above upper bound of right interval */
         if( maxvalleft < infinity )
         {
            assert(maxvalleft > -infinity);
            resultant->sup = (SCIP_Real)(maxvalleft / sqrtax);
         }
      }
      else
      {
         /* otherwise, tighten upper bound of sqrt(ax)*x to upper bound of sqrt(r(rhs,y))-b(y) */
         if( maxvalright < infinity )
         {
            assert(maxvalright > -infinity);
            resultant->sup = (SCIP_Real)(maxvalright / sqrtax);
         }
      }

      resultant->inf -= 1e-10 * REALABS(resultant->inf);
      resultant->sup += 1e-10 * REALABS(resultant->sup);

#undef CALCB
#undef CALCR
   }
   else
   {
      /* case ax == 0 */

      SCIP_Real c;
      SCIP_Real d;
      SCIP_Real ymin;
      SCIP_Real minval;
      SCIP_Real maxval;

      /* consider -bx / axy in ybnds, i.e., bx + axy y can be 0 */
      if( EPSGE(-bx / axy, ybnds.inf, 1e-9) && EPSLE(-bx / axy, ybnds.sup, 1e-9) )
      {
         /* write as (bx + axy y) * x \in (c - ay y^2 - by y)
          * and estimate bx + axy y and c - ay y^2 - by y by intervals independently
          * @todo can we do better, as in the case where bx + axy y is bounded away from 0?
          */
         SCIP_INTERVAL lincoef;
         SCIP_INTERVAL myrhs;
         SCIP_INTERVAL tmp;

         if( xbnds.inf < 0.0 && xbnds.sup > 0.0 )
         {
            /* if (bx + axy y) can be arbitrary small and x be both positive and negative,
             * then nothing we can tighten here
             */
            SCIPintervalSetBounds(resultant, xbnds.inf, xbnds.sup);
            return;
         }

         /* store interval for (bx + axy y) in lincoef */
         SCIPintervalMulScalar(infinity, &lincoef, ybnds, axy);
         SCIPintervalAddScalar(infinity, &lincoef, lincoef, bx);

         /* store interval for (c - ay y^2 - by y) in myrhs */
         SCIPintervalSet(&tmp, by);
         SCIPintervalQuad(infinity, &tmp, ay, tmp, ybnds);
         SCIPintervalSub(infinity, &myrhs, rhs, tmp);

         if( lincoef.inf == 0.0 && lincoef.sup == 0.0 )
         {
            /* equation became 0.0 \in myrhs */
            if( myrhs.inf <= 0.0 && myrhs.sup >= 0.0 )
               SCIPintervalSetEntire(infinity, resultant);
            else
               SCIPintervalSetEmpty(resultant);
         }
         else if( xbnds.inf >= 0.0 )
         {
            SCIP_INTERVAL a_;

            /* need only positive solutions */
            SCIPintervalSet(&a_, 0.0);
            SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, resultant, a_, lincoef, myrhs);
         }
         else
         {
            SCIP_INTERVAL a_;

            assert(xbnds.sup <= 0.0);

            /* need only negative solutions */
            SCIPintervalSet(&a_, 0.0);
            SCIPintervalSetBounds(&lincoef, -lincoef.sup, -lincoef.inf);
            SCIPintervalSolveUnivariateQuadExpressionPositive(infinity, resultant, a_, lincoef, myrhs);
            if( !SCIPintervalIsEmpty(infinity, *resultant) )
               SCIPintervalSetBounds(resultant, -resultant->sup, -resultant->inf);
         }

         SCIPintervalIntersect(resultant, xbnds, *resultant);

         return;
      }

      minval =  infinity;
      maxval = -infinity;

      /* compute a lower bound on x */
      if( bx + axy * (axy > 0.0 ? ybnds.inf : ybnds.sup) > 0.0 )
         c = rhs.inf;
      else
         c = rhs.sup;

      if( c > -infinity && c < infinity )
      {
         if( ybnds.inf <= -infinity )
         {
            /* limit is ay/axy * infinity if ay != 0.0 and -by/axy otherwise */
            if( EPSZ(ay, 1e-9) )
               minval = -by / axy;
            else if( ay * axy < 0.0 )
               minval = -infinity;
         }
         else
         {
            val = (c - ay * ybnds.inf * ybnds.inf - by * ybnds.inf) / (bx + axy * ybnds.inf);
            minval = MIN(val, minval);
         }

         if( ybnds.sup >= infinity )
         {
            /* limit is -ay/axy * infinity if ay != 0.0 and -by/axy otherwise */
            if( EPSZ(ay, 1e-9) )
               minval = MIN(minval, -by / axy);
            else if( ay * axy > 0.0 )
               minval = -infinity;
         }
         else
         {
            val = (c - ay * ybnds.sup * ybnds.sup - by * ybnds.sup) / (bx + axy * ybnds.sup);
            minval = MIN(val, minval);
         }

         if( !EPSZ(ay, 1e-9) )
         {
            d = ay * (ay * bx * bx - axy * (bx * by + axy * c));
            if( !EPSN(d, 1e-9) )
            {
               ymin = -ay * bx + sqrt(MAX(d, 0.0));
               ymin /= axy * ay;

               if( ymin > ybnds.inf && ymin < ybnds.sup )
               {
                  assert(bx + axy * ymin != 0.0);

                  val = (c - ay * ymin * ymin - by * ymin) / (bx + axy * ymin);
                  minval = MIN(val, minval);
               }

               ymin = -ay * bx - sqrt(MAX(d, 0.0));
               ymin /= axy * ay;

               if(ymin > ybnds.inf && ymin < ybnds.sup )
               {
                  assert(bx + axy * ymin != 0.0);

                  val = (c - ay * ymin * ymin - by * ymin) / (bx + axy * ymin);
                  minval = MIN(val, minval);
               }
            }
         }
      }
      else
      {
         minval = -infinity;
      }

      /* compute an upper bound on x */
      if( bx + axy * (axy > 0.0 ? ybnds.inf : ybnds.sup) > 0.0 )
         c = rhs.sup;
      else
         c = rhs.inf;

      if( c > -infinity && c < infinity )
      {
         if( ybnds.inf <= -infinity )
         {
            /* limit is ay/axy * infinity if ay != 0.0 and -by/axy otherwise */
            if( EPSZ(ay, 1e-9) )
               maxval = -by / axy;
            else if( ay * axy > 0.0 )
               maxval = infinity;
         }
         else
         {
            val = (c - ay * ybnds.inf * ybnds.inf - by * ybnds.inf) / (bx + axy * ybnds.inf);
            maxval = MAX(val, maxval);
         }

         if( ybnds.sup >= infinity )
         {
            /* limit is -ay/axy * infinity if ay != 0.0 and -by/axy otherwise */
            if( EPSZ(ay, 1e-9) )
               maxval = MAX(maxval, -by / axy);
            else if( ay * axy < 0.0 )
               maxval = infinity;
         }
         else
         {
            val = (c - ay * ybnds.sup * ybnds.sup - by * ybnds.sup) / (bx + axy * ybnds.sup);
            maxval = MAX(val, maxval);
         }

         if( !EPSZ(ay, 1e-9) )
         {
            d = ay * (ay * bx * bx - axy * (bx * by + axy * c));
            if( !EPSN(d, 1e-9) )
            {
               ymin = ay * bx + sqrt(MAX(d, 0.0));
               ymin /= axy * ay;

               val = (c - ay * ymin * ymin - by * ymin) / (bx + axy * ymin);
               maxval = MAX(val, maxval);

               ymin = ay * bx - sqrt(MAX(d, 0.0));
               ymin /= axy * ay;

               val = (c - ay * ymin * ymin - by * ymin) / (bx + axy * ymin);
               maxval = MAX(val, maxval);
            }
         }
      }
      else
      {
         maxval = infinity;
      }

      if( minval > -infinity )
         resultant->inf = minval - 1e-10 * REALABS(minval);
      else
         resultant->inf = -infinity;
      if( maxval <  infinity )
         resultant->sup = maxval + 1e-10 * REALABS(maxval);
      else
         resultant->sup = infinity;
      SCIPintervalIntersect(resultant, *resultant, xbnds);
   }
}
