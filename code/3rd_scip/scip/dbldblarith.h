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

/**@file   dbldblarith.h
 * @brief  defines macros for basic operations in double-double arithmetic giving roughly twice the precision of a double
 * @author Leona Gottwald
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef _SCIP_DBLDBL_ARITH_
#define _SCIP_DBLDBL_ARITH_

#include "math.h"


#ifndef DISABLE_QUADPREC

/* smaller epsilon value for use with quadprecision */
#define QUAD_EPSILON 1e-12

/* convenience macros for nicer usage of double double arithmetic */
#define QUAD_HI(x)  x ## hi
#define QUAD_LO(x)  x ## lo
#define QUAD(x) QUAD_HI(x), QUAD_LO(x)
#define QUAD_MEMBER(x) QUAD_HI(x); QUAD_LO(x)
#define QUAD_TO_DBL(x) ( QUAD_HI(x) + QUAD_LO(x) )
#define QUAD_SCALE(x, a) do { QUAD_HI(x) *= (a); QUAD_LO(x) *= (a); } while(0)
#define QUAD_ASSIGN(a, constant)  do { QUAD_HI(a) = (constant); QUAD_LO(a) = 0.0; } while(0)
#define QUAD_ASSIGN_Q(a, b)  do { QUAD_HI(a) = QUAD_HI(b); QUAD_LO(a) = QUAD_LO(b); } while(0)
#define QUAD_ARRAY_SIZE(size) ((size)*2)
#define QUAD_ARRAY_LOAD(r, a, idx) do { QUAD_HI(r) = (a)[2*(idx)]; QUAD_LO(r) = (a)[2*(idx) + 1]; } while(0)
#define QUAD_ARRAY_STORE(a, idx, x) do { (a)[2*(idx)] = QUAD_HI(x); (a)[2*(idx) + 1] = QUAD_LO(x); } while(0)

/* define all the SCIPquadprec... macros such that they use the SCIPdbldbl... macros that expands the quad precision arguments using the above macros */
#define SCIPquadprecProdDD(r, a, b)  SCIPdbldblProd(QUAD_HI(r), QUAD_LO(r), a, b)
#define SCIPquadprecSquareD(r, a) SCIPdbldblSquare(QUAD_HI(r), QUAD_LO(r), a)
#define SCIPquadprecSumDD(r, a, b) SCIPdbldblSum(QUAD_HI(r), QUAD_LO(r), a, b)
#define SCIPquadprecDivDD(r, a, b) SCIPdbldblDiv(QUAD_HI(r), QUAD_LO(r), a, b)
#define SCIPquadprecSumQD(r, a, b) SCIPdbldblSum21(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), b)
#define SCIPquadprecProdQD(r, a, b) SCIPdbldblProd21(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), b)
#define SCIPquadprecDivDQ(r, a, b) SCIPdbldblDiv12(QUAD_HI(r), QUAD_LO(r), a, QUAD_HI(b), QUAD_LO(b))
#define SCIPquadprecDivQD(r, a, b) SCIPdbldblDiv21(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), b)
#define SCIPquadprecProdQQ(r, a, b) SCIPdbldblProd22(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), QUAD_HI(b), QUAD_LO(b))
#define SCIPquadprecSumQQ(r, a, b) SCIPdbldblSum22(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), QUAD_HI(b), QUAD_LO(b))
#define SCIPquadprecSquareQ(r, a) SCIPdbldblSquare2(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a))
#define SCIPquadprecDivQQ(r, a, b) SCIPdbldblDiv22(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), QUAD_HI(b), QUAD_LO(b))
#define SCIPquadprecSqrtD(r, a) SCIPdbldblSqrt(QUAD_HI(r), QUAD_LO(r), a)
#define SCIPquadprecSqrtQ(r, a) SCIPdbldblSqrt2(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a))
#define SCIPquadprecAbsQ(r, a) SCIPdbldblAbs2(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a))
#define SCIPquadprecFloorQ(r, a) SCIPdbldblFloor2(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a))
#define SCIPquadprecCeilQ(r, a) SCIPdbldblCeil2(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a))
#define SCIPquadprecEpsFloorQ(r, a, eps) SCIPdbldblEpsFloor2(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), eps)
#define SCIPquadprecEpsCeilQ(r, a, eps) SCIPdbldblEpsCeil2(QUAD_HI(r), QUAD_LO(r), QUAD_HI(a), QUAD_LO(a), eps)

#else

/* normal epsilon value if quadprecision is disabled */
#define QUAD_EPSILON 1e-9

/* dummy macros that use normal arithmetic */
#define QUAD_HI(x)  x
#define QUAD_LO(x)  0.0
#define QUAD(x)     x
#define QUAD_MEMBER(x) x
#define QUAD_TO_DBL(x) (x)
#define QUAD_SCALE(x, a) do { (x) *= (a); } while(0)
#define QUAD_ASSIGN(a, constant)  do { (a) = constant; } while(0)
#define QUAD_ASSIGN_Q(a, b)  do { (a) = (b); } while(0)
#define QUAD_ARRAY_SIZE(size) (size)
#define QUAD_ARRAY_LOAD(r, a, idx) do { r = (a)[(idx)]; } while(0)
#define QUAD_ARRAY_STORE(a, idx, x) do { (a)[(idx)] = (x); } while(0)

#define SCIPquadprecProdDD(r, a, b)  do { (r) = (a) * (b); } while(0)
#define SCIPquadprecSquareD(r, a)    do { (r) = (a) * (a); } while(0)
#define SCIPquadprecSumDD(r, a, b)   do { (r) = (a) + (b); } while(0)
#define SCIPquadprecDivDD(r, a, b)   do { (r) = (a) / (b); } while(0)
#define SCIPquadprecSumQD(r, a, b)   do { (r) = (a) + (b); } while(0)
#define SCIPquadprecProdQD(r, a, b)  do { (r) = (a) * (b); } while(0)
#define SCIPquadprecDivDQ(r, a, b)   do { (r) = (a) / (b); } while(0)
#define SCIPquadprecDivQD(r, a, b)   do { (r) = (a) / (b); } while(0)
#define SCIPquadprecProdQQ(r, a, b)  do { (r) = (a) * (b); } while(0)
#define SCIPquadprecSumQQ(r, a, b)   do { (r) = (a) + (b); } while(0)
#define SCIPquadprecSquareQ(r, a)    do { (r) = (a) * (a); } while(0)
#define SCIPquadprecDivQQ(r, a, b)   do { (r) = (a) / (b); } while(0)
#define SCIPquadprecSqrtD(r, a)      do { (r) = sqrt(a); } while(0)
#define SCIPquadprecSqrtQ(r, a)      do { (r) = sqrt(a); } while(0)
#define SCIPquadprecAbsQ(r, a)       do { (r) = fabs(a); } while(0)
#define SCIPquadprecFloorQ(r, a)     do { (r) = floor(a); } while(0)
#define SCIPquadprecCeilQ(r, a)      do { (r) = ceil(a); } while(0)
#define SCIPquadprecEpsFloorQ(r, a, eps) do { (r) = floor((a) + (eps)); } while(0)
#define SCIPquadprecEpsCeilQ(r, a, eps) do { (r) = ceil((a) - (eps)); } while(0)

#endif

#define __SCIPdbldblSplit(rhi, rlo, x) \
    do { \
       const double __tmp_split_dbl = 134217729.0 * (x); \
       (rhi) = __tmp_split_dbl - (__tmp_split_dbl - (x)); \
       (rlo) = (x) - (rhi);\
    } while(0)

/** multiply two floating point numbers, both given by one double, and return the result as two doubles. */
#define SCIPdbldblProd(rhi, rlo, a, b) \
    do { \
        double __tmp_dbldbl_prod_ahi; \
        double __tmp_dbldbl_prod_alo; \
        double __tmp_dbldbl_prod_bhi; \
        double __tmp_dbldbl_prod_blo; \
        __SCIPdbldblSplit(__tmp_dbldbl_prod_ahi, __tmp_dbldbl_prod_alo, a); \
        __SCIPdbldblSplit(__tmp_dbldbl_prod_bhi, __tmp_dbldbl_prod_blo, b); \
        (rhi) = (a) * (b); \
        (rlo) = __tmp_dbldbl_prod_alo * __tmp_dbldbl_prod_blo - \
           ((((rhi) - __tmp_dbldbl_prod_ahi * __tmp_dbldbl_prod_bhi) \
           - __tmp_dbldbl_prod_alo * __tmp_dbldbl_prod_bhi) \
           - __tmp_dbldbl_prod_ahi * __tmp_dbldbl_prod_blo); \
    } while(0)

/** square a floating point number given by one double and return the result as two doubles. */
#define SCIPdbldblSquare(rhi, rlo, a) \
    do { \
        double __tmp_dbldbl_square_ahi; \
        double __tmp_dbldbl_square_alo; \
        __SCIPdbldblSplit(__tmp_dbldbl_square_ahi, __tmp_dbldbl_square_alo, a); \
        (rhi) = (a) * (a); \
        (rlo) = __tmp_dbldbl_square_alo * __tmp_dbldbl_square_alo - \
           ((((rhi) - __tmp_dbldbl_square_ahi * __tmp_dbldbl_square_ahi) \
           - 2.0 * __tmp_dbldbl_square_alo * __tmp_dbldbl_square_ahi)); \
    } while(0)

/** add two floating point numbers, both given by one double, and return the result as two doubles. */
#define SCIPdbldblSum(rhi, rlo, a, b) \
    do { \
        double __tmp1_dbldbl_sum; \
        double __tmp2_dbldbl_sum; \
        __tmp2_dbldbl_sum = (a) + (b); \
        __tmp1_dbldbl_sum = __tmp2_dbldbl_sum - (a); \
        (rlo) = ((a) - (__tmp2_dbldbl_sum - __tmp1_dbldbl_sum)) + ((b) - __tmp1_dbldbl_sum); \
        (rhi) = __tmp2_dbldbl_sum; \
    } while(0)

/** divide two floating point numbers, both given by one double, and return the result as two doubles. */
#define SCIPdbldblDiv(rhi, rlo, a, b) \
    do { \
       double __tmp_dbldbl_div_hi; \
       double __tmp_dbldbl_div_lo; \
       double __estim_dbldbl_div = (a)/(b); \
       SCIPdbldblProd(__tmp_dbldbl_div_hi, __tmp_dbldbl_div_lo, b, __estim_dbldbl_div); \
       SCIPdbldblSum21(__tmp_dbldbl_div_hi, __tmp_dbldbl_div_lo, __tmp_dbldbl_div_hi, __tmp_dbldbl_div_lo, -(a)); \
       __tmp_dbldbl_div_hi /= (b); \
       __tmp_dbldbl_div_lo /= (b); \
       SCIPdbldblSum21(rhi, rlo, -__tmp_dbldbl_div_hi, -__tmp_dbldbl_div_lo, __estim_dbldbl_div); \
    } while(0)

/** add two floating point numbers, the first is given by two doubles, the second is given by one double,
 *  and return the result as two doubles.
 */
#define SCIPdbldblSum21(rhi, rlo, ahi, alo, b) \
   do { \
      double __tmp_dbldbl_sum21_hi; \
      double __tmp_dbldbl_sum21_lo; \
      SCIPdbldblSum(__tmp_dbldbl_sum21_hi, __tmp_dbldbl_sum21_lo, ahi, b); \
      (rlo) = __tmp_dbldbl_sum21_lo + (alo); \
      (rhi) = __tmp_dbldbl_sum21_hi; \
   } while(0)


/** multiply two floating point numbers, the first is given by two doubles, the second is given by one double,
 *  and return the result as two doubles.
 */
#define SCIPdbldblProd21(rhi, rlo, ahi, alo, b) \
    do { \
       double __tmp_dbldbl_prod21_hi; \
       double __tmp_dbldbl_prod21_lo; \
       SCIPdbldblProd(__tmp_dbldbl_prod21_hi, __tmp_dbldbl_prod21_lo, ahi, b); \
       (rlo) = (alo) * (b) + __tmp_dbldbl_prod21_lo; \
       (rhi) = __tmp_dbldbl_prod21_hi; \
    } while(0)

/** divide two floating point numbers, the first is given by one double, the second is given by two doubles,
 *  and return the result as two doubles.
 */
#define SCIPdbldblDiv12(rhi, rlo, a, bhi, blo) \
    do { \
       double __tmp_dbldbl_div12_hi; \
       double __tmp_dbldbl_div12_lo; \
       double __estim_dbldbl_div12 = (a)/(bhi); \
       SCIPdbldblProd21(__tmp_dbldbl_div12_hi, __tmp_dbldbl_div12_lo, bhi, blo, __estim_dbldbl_div12); \
       SCIPdbldblSum21(__tmp_dbldbl_div12_hi, __tmp_dbldbl_div12_lo, __tmp_dbldbl_div12_hi, __tmp_dbldbl_div12_lo, -(a)); \
       __tmp_dbldbl_div12_hi /= (bhi); \
       __tmp_dbldbl_div12_lo /= (bhi); \
       SCIPdbldblSum21(rhi, rlo, -__tmp_dbldbl_div12_hi, -__tmp_dbldbl_div12_lo, __estim_dbldbl_div12); \
    } while(0)


/** divide two floating point numbers, the first is given by two doubles, the second is given by one double,
 *  and return the result as two doubles.
 */
#define SCIPdbldblDiv21(rhi, rlo, ahi, alo, b) \
   do { \
      double __tmp_dbldbl_div21_hi; \
      double __tmp_dbldbl_div21_lo; \
      double __estim_dbldbl_div21_hi; \
      double __estim_dbldbl_div21_lo; \
      __estim_dbldbl_div21_hi = (ahi)/(b); \
      __estim_dbldbl_div21_lo = (alo)/(b); \
      SCIPdbldblProd21(__tmp_dbldbl_div21_hi, __tmp_dbldbl_div21_lo, __estim_dbldbl_div21_hi, __estim_dbldbl_div21_lo, b); \
      SCIPdbldblSum22(__tmp_dbldbl_div21_hi, __tmp_dbldbl_div21_lo, __tmp_dbldbl_div21_hi, __tmp_dbldbl_div21_lo, -(ahi), -(alo)); \
      __tmp_dbldbl_div21_hi /= (b); \
      __tmp_dbldbl_div21_lo /= (b); \
      SCIPdbldblSum22(rhi, rlo, __estim_dbldbl_div21_hi, __estim_dbldbl_div21_lo, -__tmp_dbldbl_div21_hi, -__tmp_dbldbl_div21_lo); \
   } while(0)

/** multiply two floating point numbers, both given by two doubles, and return the result as two doubles. */
#define SCIPdbldblProd22(rhi, rlo, ahi, alo, bhi, blo) \
   do { \
      double __tmp_dbldbl_prod22_hi; \
      double __tmp_dbldbl_prod22_lo; \
      SCIPdbldblProd(__tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, ahi, bhi); \
      SCIPdbldblSum21(__tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, \
                   __tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, (alo) * (bhi)); \
      SCIPdbldblSum21(rhi, rlo, \
                   __tmp_dbldbl_prod22_hi, __tmp_dbldbl_prod22_lo, (ahi) * (blo)); \
   } while(0)

/** add two floating point numbers, both given by two doubles, and return the result as two doubles. */
#define SCIPdbldblSum22(rhi, rlo, ahi, alo, bhi, blo) \
   do { \
      double __tmp_dbldbl_sum22_hi; \
      double __tmp_dbldbl_sum22_lo; \
      SCIPdbldblSum21(__tmp_dbldbl_sum22_hi, __tmp_dbldbl_sum22_lo, ahi, alo, bhi); \
      SCIPdbldblSum21(rhi, rlo, __tmp_dbldbl_sum22_hi, __tmp_dbldbl_sum22_lo, blo); \
   } while(0)

/** square a floating point number given by two doubles and return the result as two doubles. */
#define SCIPdbldblSquare2(rhi, rlo, ahi, alo) \
   do { \
      double __tmp_dbldbl_square2_hi; \
      double __tmp_dbldbl_square2_lo; \
      SCIPdbldblSquare(__tmp_dbldbl_square2_hi, __tmp_dbldbl_square2_lo, (ahi)); \
      SCIPdbldblSum21(rhi, rlo, __tmp_dbldbl_square2_hi, __tmp_dbldbl_square2_lo, 2 * (ahi) * (alo)); \
   } while(0)

/** divide two floating point numbers, both given by two doubles, and return the result as two doubles. */
#define SCIPdbldblDiv22(rhi, rlo, ahi, alo, bhi, blo) \
   do { \
      double __tmp_dbldbl_div22_hi; \
      double __tmp_dbldbl_div22_lo; \
      double __estim_dbldbl_div22_hi = (ahi) / (bhi); \
      double __estim_dbldbl_div22_lo = (alo) / (bhi); \
      SCIPdbldblProd22(__tmp_dbldbl_div22_hi, __tmp_dbldbl_div22_lo, \
                    bhi, blo, __estim_dbldbl_div22_hi, __estim_dbldbl_div22_lo); \
      SCIPdbldblSum22(__tmp_dbldbl_div22_hi, __tmp_dbldbl_div22_lo, \
                   __tmp_dbldbl_div22_hi, __tmp_dbldbl_div22_lo, -(ahi), -(alo)); \
      __tmp_dbldbl_div22_hi /= (bhi); \
      __tmp_dbldbl_div22_lo /= (bhi); \
      SCIPdbldblSum22(rhi, rlo, __estim_dbldbl_div22_hi, __estim_dbldbl_div22_lo, \
                                -__tmp_dbldbl_div22_hi, -__tmp_dbldbl_div22_lo); \
   } while(0)


/** take the square root of a floating point number given by one double and return the result as two doubles. */
#define SCIPdbldblSqrt(rhi, rlo, a) \
   do { \
      double __estim_dbldbl_sqrt = sqrt(a); \
      if( __estim_dbldbl_sqrt != 0.0 ) \
      { \
         SCIPdbldblDiv(rhi, rlo, a, __estim_dbldbl_sqrt); \
         SCIPdbldblSum21(rhi, rlo, rhi, rlo, __estim_dbldbl_sqrt); \
         (rhi) *= 0.5; \
         (rlo) *= 0.5; \
      } \
      else \
      { \
         (rhi) = 0.0; \
         (rlo) = 0.0; \
      } \
   } while(0)


/** take the square root of a floating point number given by two doubles and return the result as two doubles. */
#define SCIPdbldblSqrt2(rhi, rlo, ahi, alo) \
   do { \
      double __estim_dbldbl_sqrt2 = sqrt(ahi + alo); \
      if( __estim_dbldbl_sqrt2 != 0.0 ) \
      { \
         SCIPdbldblDiv21(rhi, rlo, ahi, alo, __estim_dbldbl_sqrt2); \
         SCIPdbldblSum21(rhi, rlo, rhi, rlo, __estim_dbldbl_sqrt2); \
         (rhi) *= 0.5; \
         (rlo) *= 0.5; \
      } \
      else \
      { \
         (rhi) = 0.0; \
         (rlo) = 0.0; \
      } \
   } while(0)

/** compute the absolute value of the floating point number given by two doubles */
#define SCIPdbldblAbs2(rhi, rlo, ahi, alo) \
   do { \
      if( ahi < 0.0 ) \
      { \
         (rhi) = -(ahi); \
         (rlo) = -(alo); \
      } \
      else \
      { \
         (rhi) = (ahi); \
         (rlo) = (alo); \
      } \
   } while(0)

/** compute the floored value of the floating point number given by two doubles */
#define SCIPdbldblFloor2(rhi, rlo, ahi, alo) \
   do { \
      double __tmp_dbldbl_floor; \
      __tmp_dbldbl_floor = floor((ahi) + (alo)); \
      SCIPdbldblSum21(rhi, rlo, ahi, alo, -__tmp_dbldbl_floor); \
      if( ((rhi) - 1.0) + (rlo) < 0.0 && (rhi) + (rlo) >= 0.0 ) \
      { \
         /* floor in double precision was fine */ \
         (rhi) = __tmp_dbldbl_floor; \
         (rlo) = 0.0; \
      } \
      else \
      { \
         /* floor in double precision needs to be corrected */ \
         double __tmp2_dbldbl_floor = floor((rhi) + (rlo)); \
         SCIPdbldblSum(rhi, rlo, __tmp_dbldbl_floor, __tmp2_dbldbl_floor); \
      } \
   } while(0)

/** compute the ceiled value of the floating point number given by two doubles */
#define SCIPdbldblCeil2(rhi, rlo, ahi, alo) \
   do { \
      double __tmp_dbldbl_ceil; \
      __tmp_dbldbl_ceil = ceil((ahi) + (alo)); \
      SCIPdbldblSum21(rhi, rlo, -(ahi), -(alo), __tmp_dbldbl_ceil); \
      if( ((rhi) - 1.0) + (rlo) < 0.0 && (rhi) + (rlo) >= 0.0 ) \
      { \
         /* ceil in double precision was fine */ \
         (rhi) = __tmp_dbldbl_ceil; \
         (rlo) = 0.0; \
      } \
      else \
      { \
         /* ceil in double precision needs to be corrected */ \
         double __tmp2_dbldbl_ceil = floor((rhi) + (rlo)); \
         SCIPdbldblSum(rhi, rlo, __tmp_dbldbl_ceil, -__tmp2_dbldbl_ceil); \
      } \
   } while(0)

/** compute the floored value of the floating point number given by two doubles, add epsilon first for safety */
#define SCIPdbldblEpsFloor2(rhi, rlo, ahi, alo, eps) SCIPdbldblFloor2(rhi, rlo, ahi, (alo) + (eps))

/** compute the ceiled value of the floating point number given by two doubles, subtract epsilon first for safety */
#define SCIPdbldblEpsCeil2(rhi, rlo, ahi, alo, eps) SCIPdbldblCeil2(rhi, rlo, ahi, (alo) - (eps))

#endif
