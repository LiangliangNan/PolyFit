/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxdefines.h
 * @brief Debugging, floating point type and parameter definitions.
 *
 * In optimized code with \c NDEBUG defined, only
 * \ref soplex::SPxOut::INFO1 "INFO1",
 * \ref soplex::SPxOut::INFO2 "INFO2", and
 * \ref soplex::SPxOut::INFO3 "INFO3" are set.
 * If \c NDEBUG is not defined, the code within \#TRACE is used.
 * If \c SOPLEX_DEBUG is defined, the code within
 * \ref soplex::SPxOut::DEBUG "DEBUG" is also used.
 *
 * If \c WITH_LONG_DOUBLE is defined, all Real numbers are of type
 * long double instead of just double.
 */
#ifndef _SPXDEFINES_H_
#define _SPXDEFINES_H_

#include <math.h>
#ifdef _MSC_VER
#include <float.h>
#endif

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <iostream>

namespace soplex
{
#define SOPLEX_VERSION         311
#define SOPLEX_SUBVERSION        0
#define SOPLEX_APIVERSION        1
#define SOPLEX_COPYRIGHT       "Copyright (c) 1996-2018 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)"

/*-----------------------------------------------------------------------------
 * Assertion Macros etc.
 *-----------------------------------------------------------------------------
 */

/**
   \brief Macro to turn some assertions into warnings.

   If both \c NDEBUG and \c WITH_WARNINGS are defined then the failed
   assertion is converted to a warning. In all other cases this macro is
   equivalent to assert().

   @param  prefix  Short string for grepping in source code.
   @param  expr    Expression that must be satisfied.
*/
#if defined (NDEBUG) && defined (WITH_WARNINGS)
#define ASSERT_WARN( prefix, expr )                        \
   if ( !( expr ) )                                        \
      {                                                    \
         std::cerr                                         \
         << prefix                                         \
         << " failed assertion on line " << __LINE__       \
         << " in file " << __FILE__ << ": "                \
         << #expr                                          \
         << std::endl;                                     \
      }
#else // just a normal assert
#define ASSERT_WARN( prefix, expr ) ( assert( expr ) )
#endif



/*-----------------------------------------------------------------------------
 * Debugging Macros etc.
 *-----------------------------------------------------------------------------
 */

/**
   Prints/Executes \p stream with verbosity level \p verbosity, resetting
   the old verbosity level afterwards.
   Usually the parameter \p stream prints something out.
   This is an internal define used by MSG_ERROR, MSG_WARNING, etc.
*/
#ifdef DISABLE_VERBOSITY
#define DO_WITH_TMP_VERBOSITY( verbosity, spxout, do_something ) {}
#define DO_WITH_ERR_VERBOSITY( do_something ) {}
#else
#define DO_WITH_TMP_VERBOSITY( verbosity, spxout, do_something ) \
   {                                                             \
     if( &spxout != NULL )                                       \
     {                                                           \
        if( verbosity <= spxout.getVerbosity() )                 \
        {                                                        \
           const SPxOut::Verbosity  old_verbosity = spxout.getVerbosity(); \
           spxout.setVerbosity( verbosity );                     \
           do_something;                                         \
           spxout.setVerbosity( old_verbosity );                 \
        }                                                        \
     }                                                           \
   }
#define DO_WITH_ERR_VERBOSITY( do_something ) { do_something; }
#endif

/// Prints out message \p x if the verbosity level is at least SPxOut::ERROR.
#define MSG_ERROR(x)            { DO_WITH_ERR_VERBOSITY( x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::WARNING.
#define MSG_WARNING(spxout, x)  { DO_WITH_TMP_VERBOSITY( SPxOut::WARNING, spxout, x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::INFO1.
#define MSG_INFO1(spxout, x)    { DO_WITH_TMP_VERBOSITY( SPxOut::INFO1, spxout, x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::INFO2.
#define MSG_INFO2(spxout, x)    { DO_WITH_TMP_VERBOSITY( SPxOut::INFO2, spxout, x ) }
/// Prints out message \p x if the verbosity level is at least SPxOut::INFO3.
#define MSG_INFO3(spxout, x)    { DO_WITH_TMP_VERBOSITY( SPxOut::INFO3, spxout, x ) }

extern bool msginconsistent(const char* name, const char* file, int line);

#define MSGinconsistent(name) msginconsistent(name, __FILE__, __LINE__)

#if defined(SOPLEX_DEBUG)
// print output in any case, regardless of Param::verbose():
#define MSG_DEBUG(x) { x; }
#else
#define MSG_DEBUG(x) /**/
#endif //!SOPLEX_DEBUG


/*-----------------------------------------------------------------------------
 * multi-thread support
 *-----------------------------------------------------------------------------
 */
// enable the user to compile without thread_local by setting USRCXXFLAGS=-DTHREADLOCAL=""
#if !defined(THREADLOCAL)
#if defined(_MSC_VER) && _MSC_VER < 1900
#define THREADLOCAL
#else
#define THREADLOCAL thread_local
#endif
#endif

/*-----------------------------------------------------------------------------
 * Long double support, Parameters and Epsilons
 *-----------------------------------------------------------------------------
 */

#ifdef WITH_LONG_DOUBLE


typedef long double Real;

#ifndef REAL
#define REAL(x)  x##L
#define REAL_FORMAT "Lf"
#endif
/// default allowed bound violation
#ifndef DEFAULT_BND_VIOL
#define DEFAULT_BND_VIOL   1e-12L
#endif
/// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
#ifndef DEFAULT_EPS_ZERO
#define DEFAULT_EPS_ZERO   1e-28L
#endif
/// epsilon for factorization
#ifndef DEFAULT_EPS_FACTOR
#define DEFAULT_EPS_FACTOR 1e-30L
#endif
/// epsilon for factorization update
#ifndef DEFAULT_EPS_UPDATE
#define DEFAULT_EPS_UPDATE 1e-26L
#endif
#ifndef DEFAULT_EPS_PIVOT
#define DEFAULT_EPS_PIVOT 1e-20L
#endif
///
#define DEFAULT_INFINITY   1e100L


#else

#ifdef WITH_FLOAT

typedef float Real;

#ifndef REAL
#define REAL(x)  x
#define REAL_FORMAT "f"
#endif
/// default allowed bound violation
#ifndef DEFAULT_BND_VIOL
#define DEFAULT_BND_VIOL   1e-1
#endif
/// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
#ifndef DEFAULT_EPS_ZERO
#define DEFAULT_EPS_ZERO   1e-7
#endif
#ifndef DEFAULT_EPS_FACTOR
#define DEFAULT_EPS_FACTOR 1e-7
#endif
#ifndef DEFAULT_EPS_UPDATE
#define DEFAULT_EPS_UPDATE 1e-6
#endif
#ifndef DEFAULT_EPS_PIVOT
#define DEFAULT_EPS_PIVOT 1e-6
#endif
#define DEFAULT_INFINITY   1e100


#else

typedef double Real;

#ifndef REAL
#define REAL(x)  x
#define REAL_FORMAT "lf"
#endif
/// default allowed bound violation
#ifndef DEFAULT_BND_VIOL
#define DEFAULT_BND_VIOL   1e-6
#endif
/// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
#ifndef DEFAULT_EPS_ZERO
#define DEFAULT_EPS_ZERO   1e-16
#endif
#ifndef DEFAULT_EPS_FACTOR
#define DEFAULT_EPS_FACTOR 1e-20
#endif
#ifndef DEFAULT_EPS_UPDATE
#define DEFAULT_EPS_UPDATE 1e-16
#endif
#ifndef DEFAULT_EPS_PIVOT
#define DEFAULT_EPS_PIVOT 1e-10
#endif
#define DEFAULT_INFINITY   1e100

#endif // !WITH_FLOAT
#endif // !WITH_LONG_DOUBLE

#define MAXIMUM(x,y)        ((x)>(y) ? (x) : (y))
#define MINIMUM(x,y)        ((x)<(y) ? (x) : (y))

#define SPX_MAXSTRLEN       1024 /**< maximum string length in SoPlex */

THREADLOCAL extern const Real infinity;

class Param
{
private:

   //------------------------------------
   /**@name Data */
   //@{
   /// default allowed additive zero: 1.0 + EPS_ZERO == 1.0
   THREADLOCAL static Real s_epsilon;
   /// epsilon for factorization
   THREADLOCAL static Real s_epsilon_factorization;
   /// epsilon for factorization update
   THREADLOCAL static Real s_epsilon_update;
   /// epsilon for pivot zero tolerance in factorization
   THREADLOCAL static Real s_epsilon_pivot;
   //@}

public:

   //------------------------------------
   /**@name Access / modification */
   //@{
   ///
   static Real epsilon();
   ///
   static void setEpsilon(Real eps);
   ///
   static Real epsilonFactorization();
   ///
   static void setEpsilonFactorization(Real eps);
   ///
   static Real epsilonUpdate();
   ///
   static void setEpsilonUpdate(Real eps);
   ///
   static Real epsilonPivot();
   ///
   static void setEpsilonPivot(Real eps);
   //@}
};

#ifdef WITH_LONG_DOUBLE
/// returns |a|
inline Real spxAbs(Real a)
{
   return fabsl(a);
}

/// returns square root
inline Real spxSqrt(Real a)
{
   return sqrtl(a);
}

// returns the next representable value after x in the direction of y
#ifndef SOPLEX_LEGACY
inline Real spxNextafter(Real x, Real y)
{
   return nextafterl(x,y);
}
#endif

/// returns x * 2^exp
inline Real spxLdexp(Real x, int exp)
{
   return ldexpl(x,exp);
}

// returns x and exp such that y = x * 2^exp
inline Real spxFrexp(Real y, int* exp)
{
   return frexpl(y, exp);
}
#else
/// returns |a|
inline Real spxAbs(Real a)
{
   return fabs(a);
}

/// returns square root
inline Real spxSqrt(Real a)
{
   return sqrt(a);
}

// returns the next representable value after x in the direction of y
#ifndef SOPLEX_LEGACY
inline Real spxNextafter(Real x, Real y)
{
#ifndef _MSC_VER
   return nextafter(x,y);
#else
   return _nextafter(x,y);
#endif
}
#endif

/// returns x * 2^exp
inline Real spxLdexp(Real x, int exp)
{
   return ldexp(x,exp);
}

// returns x and exp such that y = x * 2^exp
inline Real spxFrexp(Real y, int* exp)
{
   return frexp(y, exp);
}
#endif

/// returns max(|a|,|b|)
inline Real maxAbs(Real a, Real b)
{
   const Real absa = spxAbs(a);
   const Real absb = spxAbs(b);

   return absa > absb ? absa : absb;
}

/// returns (a-b) / max(|a|,|b|,1.0)
inline Real relDiff(Real a, Real b)
{
   return (a - b) / (maxAbs(a, b) > 1.0 ? maxAbs(a, b) : 1.0);
}

/// returns \c true iff |a-b| <= eps
inline bool EQ(Real a, Real b, Real eps = Param::epsilon())
{
   return spxAbs(a - b) <= eps;
}

/// returns \c true iff |a-b| > eps
inline bool NE(Real a, Real b, Real eps = Param::epsilon())
{
   return spxAbs(a - b) > eps;
}

/// returns \c true iff a < b + eps
inline bool LT(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) < -eps;
}

/// returns \c true iff a <= b + eps
inline bool LE(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) < eps;
}

/// returns \c true iff a > b + eps
inline bool GT(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) > eps;
}

/// returns \c true iff a >= b + eps
inline bool GE(Real a, Real b, Real eps = Param::epsilon())
{
   return (a - b) > -eps;
}

/// returns \c true iff |a| <= eps
inline bool isZero(Real a, Real eps = Param::epsilon())
{
   return spxAbs(a) <= eps;
}

/// returns \c true iff |a| > eps
inline bool isNotZero(Real a, Real eps = Param::epsilon())
{
   return spxAbs(a) > eps;
}

/// returns \c true iff |relDiff(a,b)| <= eps
inline bool EQrel(Real a, Real b, Real eps = Param::epsilon())
{
   return spxAbs(relDiff(a, b)) <= eps;
}

/// returns \c true iff |relDiff(a,b)| > eps
inline bool NErel(Real a, Real b, Real eps = Param::epsilon())
{
   return spxAbs(relDiff(a, b)) > eps;
}

/// returns \c true iff relDiff(a,b) <= -eps
inline bool LTrel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) <= -eps;
}

/// returns \c true iff relDiff(a,b) <= eps
inline bool LErel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) <= eps;
}

/// returns \c true iff relDiff(a,b) > eps
inline bool GTrel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) > eps;
}

/// returns \c true iff relDiff(a,b) > -eps
inline bool GErel(Real a, Real b, Real eps = Param::epsilon())
{
   return relDiff(a, b) > -eps;
}

/// safe version of snprintf
inline int spxSnprintf(
   char*                 t,                  /**< target string */
   size_t                len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   )
{
   va_list ap;
   int n;

   assert(t != NULL);
   assert(len > 0);

   va_start(ap, s); /*lint !e826*/

#if defined(_WIN32) || defined(_WIN64)
   n = _vsnprintf(t, len, s, ap);
#else
   n = vsnprintf(t, len, s, ap); /*lint !e571*/
#endif
   va_end(ap);

   if( n < 0 || (size_t) n >= len )
   {
#ifndef NDEBUG
      if( n < 0 )
      {
         MSG_ERROR( std::cerr << "vsnprintf returned " << n << " while reading: " << s << std::endl; )
      }
#endif
      t[len-1] = '\0';
      n = (int) len-1;
   }
   return n;
}

} // namespace soplex
#endif // _SPXDEFINES_H_
