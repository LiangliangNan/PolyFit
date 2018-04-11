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

/**@file   intervalarithext.h
 * @brief  C++ extensions to interval arithmetics for provable bounds
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_INTERVALARITHEXT_HPP__
#define __SCIP_INTERVALARITHEXT_HPP__

#include "scip/intervalarith.h"
#include "scip/pub_message.h"

/* In some cases, a user may need the SCIPInterval class to be defined in some namespace.
 * To allow this, the symbol SCIPInterval_NAMESPACE should be defined to the name of that namespace
 * before inclusion of this header file
 * It is in the users responsibility to implement the corresponding symbols from nlpi/intervalarith.c.
 */
#ifdef SCIPInterval_NAMESPACE
namespace SCIPInterval_NAMESPACE {
#endif

/** an interval that extends the SCIP_INTERVAL struct
 *
 *  Used by various methods to allow calculating with intervals as with ordinary numbers.
 */
class SCIPInterval : public SCIP_INTERVAL
{
public:
   /** value to use for infinity
    *
    *  Currently a global variable, thus use with care!
    */
   static SCIP_Real infinity;

   /** default constructor -> gives entire interval */
   SCIPInterval()
   {
      SCIPintervalSetBounds(this, -SCIPInterval::infinity,SCIPInterval::infinity);
   }

   /** constructor for an SCIP_INTERVAL struct */
   SCIPInterval(
      const SCIP_INTERVAL& x                 /**< interval to copy */
      )
   {
      SCIPintervalSetBounds(this, x.inf, x.sup);
   }

   /** constructor for an interval giving pointwise bounds */
   SCIPInterval(
      SCIP_Real          infinum,            /**< lower bound of interval */
      SCIP_Real          supremum            /**< upper bound of interval */
      )
   {
      SCIPintervalSetBounds(this, infinum, supremum);
   }

   /** constructor for a singleton */
   SCIPInterval(
      SCIP_Real          number              /**< number to be represented by interval */
      )
   {
      SCIPintervalSet(this, number);
   }

   /** sets interval bounds */
   void setBounds(
      SCIP_Real          newinf,             /**< lower bound of interval */
      SCIP_Real          newsup              /**< upper bound of interval */
      )
   {
      SCIPintervalSetBounds(this, newinf, newsup);
   }

   /** returns whether this interval is equal to another one */
   bool operator==(
      const SCIP_INTERVAL& y                 /**< interval to compare with */
      ) const
   {
      if( SCIPintervalIsEmpty(SCIPInterval::infinity, *this) && !SCIPintervalIsEmpty(SCIPInterval::infinity, y) )
         return false;
      if( this->inf <= -SCIPInterval::infinity && y.inf > -SCIPInterval::infinity )
         return false;
      if( this->sup >= SCIPInterval::infinity && y.sup < SCIPInterval::infinity )
         return false;
      return (this->inf == y.inf) && (this->sup == y.sup);
   }

   /** returns whether this interval is not equal to another one */
   bool operator!=(
      const SCIP_INTERVAL& y                 /**< interval to compare with */
      ) const
   {
      return !operator==(y);
   }

   /** returns whether this interval is equal to a given number */
   bool operator==(
      const SCIP_Real&   y                   /**< number to compare with */
      ) const
   {
      return ( (inf == y) && (sup == y) ) ||
             ( sup <= -SCIPInterval::infinity && y <= -SCIPInterval::infinity ) ||
             ( inf >= SCIPInterval::infinity && y >= SCIPInterval::infinity );
   }

   /** adds another interval to this one */
   SCIPInterval& operator+=(
      const SCIPInterval& y                  /**< interval to add in */
      )
   {
      SCIPintervalAdd(SCIPInterval::infinity, this, *this, y);
      return *this;
   }

   /** subtracts an interval from this one */
   SCIPInterval& operator-=(
      const SCIPInterval& y                  /**< interval to substract */
      )
   {
      SCIPintervalSub(SCIPInterval::infinity, this, *this, y);
      return *this;
   }

   /** multiplies an interval with this one */
   SCIPInterval& operator*=(
      const SCIPInterval& y                  /**< interval to multiply with */
      )
   {
      SCIPintervalMul(SCIPInterval::infinity, this, *this, y);
      return *this;
   }

   /** divides this interval by another one */
   SCIPInterval& operator/=(
      const SCIPInterval& y                  /**< interval to divide by */
      )
   {
      SCIPintervalDiv(SCIPInterval::infinity, this, *this, y);
      return *this;
   }

   /** assigns a number to this interval */
   SCIPInterval& operator=(
      const SCIP_Real& y                     /**< new value for both interval bounds */
      )
   {
      SCIPintervalSet(this, y);
      return *this;
   }

   /** assign an interval to this interval */
   SCIPInterval& operator=(
      const SCIP_INTERVAL& y                 /**< new value for this interval */
      )
   {
      SCIPintervalSetBounds(this, y.inf, y.sup);
      return *this;
   }

};

/** addition of two intervals */
inline
SCIPInterval operator+(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalAdd(SCIPInterval::infinity, &resultant, x, y);

   return resultant;
}

/** substraction for two intervals */
inline
SCIPInterval operator-(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalSub(SCIPInterval::infinity, &resultant, x, y);

   return resultant;
}

/** negation of an interval */
inline
SCIPInterval operator-(
   const SCIPInterval&   y                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalSetBounds(&resultant, -y.sup, -y.inf);

   return resultant;
}

/** multiplication of two intervals */
inline
SCIPInterval operator*(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalMul(SCIPInterval::infinity, &resultant, x, y);

   return resultant;
}

/** division for two intervals */
inline
SCIPInterval operator/(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalDiv(SCIPInterval::infinity, &resultant, x, y);

   return resultant;
}

/** cosine of an interval */
inline
SCIPInterval cos(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalCos(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** exponential of an interval */
inline
SCIPInterval exp(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalExp(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** natural logarithm of an interval */
inline
SCIPInterval log(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalLog(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** power of an interval to another interval */
inline
SCIPInterval pow(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIPInterval&   y                   /**< second operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalPower(SCIPInterval::infinity, &resultant, x, y);

   return resultant;
}

/** power of an interval to a scalar */
inline
SCIPInterval pow(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIP_Real&      y                   /**< exponent */
   )
{
   SCIPInterval resultant;

   SCIPintervalPowerScalar(SCIPInterval::infinity, &resultant, x, y);

   return resultant;
}

/** signpower of an interval to a scalar */
inline
SCIPInterval signpow(
   const SCIPInterval&   x,                  /**< first operand */
   const SCIP_Real       p                   /**< exponent */
   )
{
   SCIPInterval resultant;

   SCIPintervalSignPowerScalar(SCIPInterval::infinity, &resultant, x, p);

   return resultant;
}

/** sine of an interval */
inline
SCIPInterval sin(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalSin(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** square an interval */
inline
SCIPInterval square(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalSquare(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** square root of an interval */
inline
SCIPInterval sqrt(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalSquareRoot(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** absolute value of an interval */
inline
SCIPInterval abs(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalAbs(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** sign of an interval */
inline
SCIPInterval sign(
   const SCIPInterval&   x                   /**< operand */
   )
{
   SCIPInterval resultant;

   SCIPintervalSign(SCIPInterval::infinity, &resultant, x);

   return resultant;
}

/** macro for easy definition of not implemented interval functions */
#define SCIP_INTERVALARITH_UNDEFFUNC(function)                                  \
inline                                                                          \
SCIPInterval function(                                                          \
   const SCIPInterval&   x                   /**< operand */                    \
   )                                                                            \
{                                                                               \
   SCIPerrorMessage("Error: " #function " not implemented for intervals.\n");   \
   return SCIPInterval();                                                       \
}

SCIP_INTERVALARITH_UNDEFFUNC(tan)
SCIP_INTERVALARITH_UNDEFFUNC(acos)
SCIP_INTERVALARITH_UNDEFFUNC(acosh)
SCIP_INTERVALARITH_UNDEFFUNC(asin)
SCIP_INTERVALARITH_UNDEFFUNC(asinh)
SCIP_INTERVALARITH_UNDEFFUNC(atan)
SCIP_INTERVALARITH_UNDEFFUNC(atanh)
SCIP_INTERVALARITH_UNDEFFUNC(cosh)
SCIP_INTERVALARITH_UNDEFFUNC(sinh)
SCIP_INTERVALARITH_UNDEFFUNC(tanh)
SCIP_INTERVALARITH_UNDEFFUNC(erf)
SCIP_INTERVALARITH_UNDEFFUNC(expm1)
SCIP_INTERVALARITH_UNDEFFUNC(log1p)
#undef SCIP_INTERVALARITH_UNDEFFUNC

#ifdef SCIPInterval_NAMESPACE
} /* namespace SCIPInterval_NAMESPACE */
#endif

#endif
