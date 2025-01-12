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

/**@file   def.h
 * @ingroup INTERNALAPI
 * @brief  common defines and data types used in all packages of SCIP
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEF_H__
#define __SCIP_DEF_H__

#ifdef __cplusplus
#define __STDC_LIMIT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

/*
 * include build configuration flags
 */


/*
 * GNU COMPILER VERSION define
 */
#ifdef __GNUC__
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 100                     \
      + __GNUC_MINOR__ * 10                             \
      + __GNUC_PATCHLEVEL__)
#endif
#endif

/*
 * define whether compiler allows variadic macros
 * __STDC_VERSION__ only exists for C code
 * added the extra check using the GCC_VERSION to enable variadic macros also with C++ code with GCC atleast
 *
 */
#if defined(_MSC_VER) || ( __STDC_VERSION__ >= 199901L ) || ( GCC_VERSION >= 480 )
#define SCIP_HAVE_VARIADIC_MACROS 1
#endif

/** get the first parameter and all-but-the-first arguments from variadic arguments
 *
 * normally, SCIP_VARARGS_FIRST_ should be sufficient
 * the SCIP_VARARGS_FIRST_/SCIP_VARARGS_FIRST kludge is to work around a bug in MSVC (https://stackoverflow.com/questions/4750688/how-to-single-out-the-first-parameter-sent-to-a-macro-taking-only-a-variadic-par)
 */
#define SCIP_VARARGS_FIRST_(firstarg, ...) firstarg
#define SCIP_VARARGS_FIRST(args) SCIP_VARARGS_FIRST_ args

/** get all but the first parameter from variadic arguments */
#define SCIP_VARARGS_REST(firstarg, ...) __VA_ARGS__

/*
 * Boolean values
 */

#ifndef SCIP_Bool
#define SCIP_Bool unsigned int               /**< type used for Boolean values */
#ifndef TRUE
#define TRUE  1                              /**< Boolean value TRUE */
#define FALSE 0                              /**< Boolean value FALSE */
#endif
#endif

#ifndef SCIP_Shortbool
#define SCIP_Shortbool uint8_t               /**< type used for Boolean values with less space */
#endif

/*
 * Add some macros for differing functions on Windows
 */
#if defined(_WIN32) || defined(_WIN64)

#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#define getcwd _getcwd
#endif

/*
 * Define the marco SCIP_EXPORT if it is not included from the generated header
 */
#ifndef SCIP_EXPORT
#if defined(_WIN32) || defined(_WIN64)
#define SCIP_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#define SCIP_EXPORT __attribute__((__visibility__("default")))
#else
#define SCIP_EXPORT
#endif
#endif

/* define INLINE */
#ifndef INLINE
#if defined(_WIN32) || defined(__STDC__)
#define INLINE                 __inline
#else
#define INLINE                 inline
#endif
#endif



#include "scip/type_retcode.h"
#include "scip/type_message.h"
#include "scip/pub_message.h"

#ifdef __cplusplus
extern "C" {
#endif


#define SCIP_VERSION                803 /**< SCIP version number (multiplied by 100 to get integer number) */
#define SCIP_SUBVERSION               0 /**< SCIP sub version number */
#define SCIP_APIVERSION             104 /**< SCIP API version number */
#define SCIP_COPYRIGHT   "Copyright (C) 2002-2022 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)"


/*
 * CIP format variable characters
 */

#define SCIP_VARTYPE_BINARY_CHAR 'B'
#define SCIP_VARTYPE_INTEGER_CHAR 'I'
#define SCIP_VARTYPE_IMPLINT_CHAR 'M'
#define SCIP_VARTYPE_CONTINUOUS_CHAR 'C'

/*
 * Long Integer values
 */

#ifndef LLONG_MAX
#define LLONG_MAX        9223372036854775807LL
#define LLONG_MIN        (-LLONG_MAX - 1LL)
#endif

#define SCIP_Longint long long                         /**< type used for long integer values */
#define SCIP_LONGINT_MAX          LLONG_MAX
#define SCIP_LONGINT_MIN          LLONG_MIN
#ifndef SCIP_LONGINT_FORMAT
#if defined(_WIN32) || defined(_WIN64)
#define SCIP_LONGINT_FORMAT           "I64d"
#else
#define SCIP_LONGINT_FORMAT           "lld"
#endif
#endif

/*
 * Floating point values
 */

#define SCIP_Real double                               /**< type used for floating point values */
#define SCIP_REAL_MAX         (SCIP_Real)DBL_MAX
#define SCIP_REAL_MIN        -(SCIP_Real)DBL_MAX
#define SCIP_REAL_FORMAT               "lf"

#define SCIP_DEFAULT_INFINITY         1e+20  /**< default value considered to be infinity */
#define SCIP_DEFAULT_EPSILON          1e-09  /**< default upper bound for floating points to be considered zero */
#define SCIP_DEFAULT_SUMEPSILON       1e-06  /**< default upper bound for sums of floating points to be considered zero */
#define SCIP_DEFAULT_FEASTOL          1e-06  /**< default feasibility tolerance for constraints */
#define SCIP_DEFAULT_CHECKFEASTOLFAC    1.0  /**< default factor to change the feasibility tolerance when testing the best solution for feasibility (after solving process) */
#define SCIP_DEFAULT_LPFEASTOLFACTOR    1.0  /**< default factor w.r.t. primal feasibility tolerance that determines default (and maximal) primal feasibility tolerance of LP solver */
#define SCIP_DEFAULT_DUALFEASTOL      1e-07  /**< default feasibility tolerance for reduced costs */
#define SCIP_DEFAULT_BARRIERCONVTOL   1e-10  /**< default convergence tolerance used in barrier algorithm */
#define SCIP_DEFAULT_BOUNDSTREPS       0.05  /**< default minimal relative improve for strengthening bounds */
#define SCIP_DEFAULT_PSEUDOCOSTEPS    1e-01  /**< default minimal variable distance value to use for pseudo cost updates */
#define SCIP_DEFAULT_PSEUDOCOSTDELTA  1e-04  /**< default minimal objective distance value to use for pseudo cost updates */
#define SCIP_DEFAULT_RECOMPFAC        1e+07  /**< default minimal decrease factor that causes the recomputation of a value (e.g., pseudo objective) instead of an update */
#define SCIP_DEFAULT_HUGEVAL          1e+15  /**< values larger than this are considered huge and should be handled separately (e.g., in activity computation) */
#define SCIP_MAXEPSILON               1e-03  /**< maximum value for any numerical epsilon */
#define SCIP_MINEPSILON               1e-20  /**< minimum value for any numerical epsilon */
#define SCIP_INVALID          (double)1e+99  /**< floating point value is not valid */
#define SCIP_UNKNOWN          (double)1e+98  /**< floating point value is not known (in primal solution) */
#define SCIP_INTERVAL_INFINITY (double)1e+300 /**< infinity value for interval computations */

#define REALABS(x)        (fabs(x))
#define EPSEQ(x,y,eps)    (REALABS((x)-(y)) <= (eps))
#define EPSLT(x,y,eps)    ((x)-(y) < -(eps))
#define EPSLE(x,y,eps)    ((x)-(y) <= (eps))
#define EPSGT(x,y,eps)    ((x)-(y) > (eps))
#define EPSGE(x,y,eps)    ((x)-(y) >= -(eps))
#define EPSZ(x,eps)       (REALABS(x) <= (eps))
#define EPSP(x,eps)       ((x) > (eps))
#define EPSN(x,eps)       ((x) < -(eps))
#define EPSFLOOR(x,eps)   (floor((x)+(eps)))
#define EPSCEIL(x,eps)    (ceil((x)-(eps)))
#define EPSROUND(x,eps)   (ceil((x)-0.5+(eps)))
#define EPSFRAC(x,eps)    ((x)-EPSFLOOR(x,eps))
#define EPSISINT(x,eps)   (EPSFRAC(x,eps) <= (eps))


#ifndef SQR
#define SQR(x)        ((x)*(x))
#define SQRT(x)       (sqrt(x))
#endif

/* platform-dependent specification of the log1p, which is numerically more stable around x = 0.0 */
#ifndef LOG1P
#if defined(_WIN32) || defined(_WIN64)
#define LOG1P(x) (log(1.0+x))
#else
#define LOG1P(x) (log1p(x))
#endif
#endif

#ifndef LOG2
#if defined(_MSC_VER) && (_MSC_VER < 1800)
#define LOG2(x) (log(x) / log(2.0))
#else
#define LOG2(x) log2(x)
#endif
#endif

#ifndef ABS
#define ABS(x)        ((x) >= 0 ? (x) : -(x))
#endif

#ifndef MAX
#define MAX(x, y) ((x) >= (y) ? (x) : (y)) /**< returns maximum of x and y */
#endif

#ifndef MIN
#define MIN(x, y) ((x) <= (y) ? (x) : (y)) /**< returns minimum of x and y */
#endif

#ifndef MAX3
#define MAX3(x, y, z) ((x) >= (y) ? MAX(x, z) : MAX(y, z)) /**< returns maximum of x, y, and z */
#endif

#ifndef MIN3
#define MIN3(x, y, z) ((x) <= (y) ? MIN(x, z) : MIN(y, z)) /**< returns minimum of x, y, and z */
#endif

#ifndef COPYSIGN
#if defined(_MSC_VER) && (_MSC_VER < 1800)
#define COPYSIGN _copysign
#else
#define COPYSIGN copysign
#endif
#endif

/*
 * Pointers
 */

#ifndef NULL
#define NULL ((void*)0)                 /**< zero pointer */
#endif

#ifndef RESTRICT
#if defined(_MSC_VER)
#define RESTRICT __restrict
#else
#ifdef __cplusplus
#define RESTRICT __restrict__
#elif __STDC_VERSION__ >= 199901L
#define RESTRICT restrict
#else
#define RESTRICT
#endif
#endif
#endif

/*
 * Strings
 */

#define SCIP_MAXSTRLEN             1024 /**< maximum string length in SCIP */

/*
 * Memory settings
 */

/* we use SIZE_MAX / 2 to detect negative sizes which got a very large value when casting to size_t */
#define SCIP_MAXMEMSIZE              (SIZE_MAX/2) /**< maximum size of allocated memory (array) */

#define SCIP_HASHSIZE_PARAMS        2048 /**< size of hash table in parameter name tables */
#define SCIP_HASHSIZE_NAMES          500 /**< size of hash table in name tables */
#define SCIP_HASHSIZE_CUTPOOLS       500 /**< size of hash table in cut pools */
#define SCIP_HASHSIZE_CLIQUES        500 /**< size of hash table in clique tables */
#define SCIP_HASHSIZE_NAMES_SMALL    100 /**< size of hash table in name tables for small problems */
#define SCIP_HASHSIZE_CUTPOOLS_SMALL 100 /**< size of hash table in cut pools for small problems */
#define SCIP_HASHSIZE_CLIQUES_SMALL  100 /**< size of hash table in clique tables for small problems */
#define SCIP_HASHSIZE_VBC            500 /**< size of hash map for node -> nodenum mapping used for VBC output */

#define SCIP_DEFAULT_MEM_ARRAYGROWFAC   1.2 /**< memory growing factor for dynamically allocated arrays */
#define SCIP_DEFAULT_MEM_ARRAYGROWINIT    4 /**< initial size of dynamically allocated arrays */

#define SCIP_MEM_NOLIMIT (SCIP_Longint)(SCIP_LONGINT_MAX >> 20)/**< initial size of dynamically allocated arrays */

/*
 * Tree settings
 */

#define SCIP_MAXTREEDEPTH             65534  /**< maximal allowed depth of the branch-and-bound tree */

/*
 * Probing scoring settings
 */

#define SCIP_PROBINGSCORE_PENALTYRATIO    2  /**< ratio for penalizing too small fractionalities in diving heuristics.
                                              *   if the fractional part of a variable is smaller than a given threshold
                                              *   the corresponding score gets penalized. due to numerical troubles
                                              *   we will flip a coin whenever SCIPisEQ(scip, fractionality, threshold)
                                              *   evaluates to true. this parameter defines the chance that this results
                                              *   in penalizing the score, i.e., there is 1:2 chance for penalizing.
                                              */

/*
 * Global debugging settings
 */

/*#define DEBUG*/


/*
 * Defines for handling SCIP return codes
 */

/** this macro is used to stop SCIP in debug mode such that errors can be debugged;
 *
 *  @note In optimized mode this macro has no effect. That means, in case of an error it has to be ensured that code
 *        terminates with an error code or continues safely.
 */
#define SCIPABORT() assert(FALSE) /*lint --e{527} */

#define SCIP_CALL_ABORT_QUIET(x)  do { if( (x) != SCIP_OKAY ) SCIPABORT(); } while( FALSE )
#define SCIP_CALL_QUIET(x)        do { SCIP_RETCODE _restat_; if( (_restat_ = (x)) != SCIP_OKAY ) return _restat_; } while( FALSE )
#define SCIP_ALLOC_ABORT_QUIET(x) do { if( NULL == (x) ) SCIPABORT(); } while( FALSE )
#define SCIP_ALLOC_QUIET(x)       do { if( NULL == (x) ) return SCIP_NOMEMORY; } while( FALSE )

#define SCIP_CALL_ABORT(x) do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _restat_; /*lint -e{506,774}*/                                         \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             SCIPABORT();                                                                     \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_ALLOC_ABORT(x) do                                                                                \
                       {                                                                                      \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             SCIPABORT();                                                                     \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_CALL(x)   do                                                                                     \
                       {                                                                                      \
                          SCIP_RETCODE _restat_; /*lint -e{506,774}*/                                         \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             return _restat_;                                                                 \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_ALLOC(x)  do                                                                                     \
                       {                                                                                      \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_CALL_TERMINATE(retcode, x, TERM)   do                                                            \
                       {                                                                                      \
                          if( ((retcode) = (x)) != SCIP_OKAY )                                                \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", retcode);                      \
                             goto TERM;                                                                       \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_ALLOC_TERMINATE(retcode, x, TERM)   do                                                           \
                       {                                                                                      \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             retcode = SCIP_NOMEMORY;                                                         \
                             goto TERM;                                                                       \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_CALL_FINALLY(x, y)   do                                                                                     \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             (y);                                                                             \
                             return _restat_;                                                                 \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

#define SCIP_UNUSED(x) ((void) (x))

/*
 * Define to mark deprecated API functions
 */

#ifndef SCIP_DEPRECATED
#if defined(_MSC_VER)
#  define SCIP_DEPRECATED __declspec(deprecated)
#elif defined(__GNUC__)
#  define SCIP_DEPRECATED __attribute__ ((deprecated))
#else
#  define SCIP_DEPRECATED
#endif
#endif

#ifdef __cplusplus
}
#endif

#endif
