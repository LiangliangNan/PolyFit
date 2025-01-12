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

/**@file    nlpi_ipopt_dummy.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   dummy Ipopt NLP interface for the case that Ipopt is not available
 * @author  Stefan Vigerske
 * @author  Benjamin Müller
 *
 * This code has been separated from nlpi_ipopt.cpp, so the SCIP build system recognizes it as pure C code,
 * thus the linker does not need to be changed to C++.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/nlpi_ipopt.h"
#include "blockmemshell/memory.h"

/** create solver interface for Ipopt solver and includes it into SCIP, if Ipopt is available */  /*lint -e{715}*/
SCIP_RETCODE SCIPincludeNlpSolverIpopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}  /*lint !e715*/

/** gets string that identifies Ipopt (version number) */
const char* SCIPgetSolverNameIpopt(void)
{
   return "";
}

/** gets string that describes Ipopt */
const char* SCIPgetSolverDescIpopt(void)
{
   return "";
}

/** returns whether Ipopt is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisIpoptAvailableIpopt(void)
{
   return FALSE;
}

/** gives a pointer to the NLPIORACLE object stored in Ipopt-NLPI's NLPI problem data structure */ /*lint -e715*/
void* SCIPgetNlpiOracleIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   )
{
   SCIPerrorMessage("Ipopt not available!\n");
   SCIPABORT();
   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** Calls Lapacks Dsyev routine to compute eigenvalues and eigenvectors of a dense matrix.
 * It's here, because Ipopt is linked against Lapack.
 */ /*lint -e715*/
SCIP_RETCODE SCIPcallLapackDsyevIpopt(
   SCIP_Bool             computeeigenvectors,/**< should also eigenvectors should be computed ? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if computeeigenvectors == TRUE */
   SCIP_Real*            w                   /**< buffer to store eigenvalues (size N) */
   )
{
   SCIPerrorMessage("Ipopt not available, cannot use it's Lapack link!\n");
   return SCIP_ERROR;
}  /*lint !e715*/

/* easier access to the entries of A */
#define ENTRY(i,j) (N * (j) + (i))

/* solves a linear problem of the form Ax = b for a regular 3*3 matrix A */
static
SCIP_RETCODE solveLinearProb3(
   SCIP_Real*            A,                  /**< matrix data on input (size 3*3); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size 3) */
   SCIP_Real*            x,                  /**< buffer to store solution (size 3) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   )
{
   SCIP_Real LU[9];
   SCIP_Real y[3];
   int pivot[3] = {0, 1, 2};
   const int N = 3;
   int k;

   assert(A != NULL);
   assert(b != NULL);
   assert(x != NULL);
   assert(success != NULL);

   *success = TRUE;

   /* copy arrays */
   BMScopyMemoryArray(LU, A, N*N);
   BMScopyMemoryArray(y, b, N);

   /* first step: compute LU factorization */
   for( k = 0; k < N; ++k )
   {
      int p;
      int i;

      p = k;
      for( i = k+1; i < N; ++i )
      {
         if( ABS(LU[ ENTRY(pivot[i],k) ]) > ABS( LU[ ENTRY(pivot[p],k) ]) )
            p = i;
      }

      if( ABS(LU[ ENTRY(pivot[p],k) ]) < 1e-08 )
      {
         /* SCIPerrorMessage("Error in nlpi_ipopt_dummy - matrix is singular!\n"); */
         *success = FALSE;
         return SCIP_OKAY;
      }

      if( p != k )
      {
         int tmp;

         tmp = pivot[k];
         pivot[k] = pivot[p];
         pivot[p] = tmp;
      }

      for( i = k+1; i < N; ++i )
      {
         SCIP_Real m;
         int j;

         m = LU[ ENTRY(pivot[i],k) ] / LU[ ENTRY(pivot[k],k) ];

         for( j = k+1; j < N; ++j )
            LU[ ENTRY(pivot[i],j) ] -= m * LU[ ENTRY(pivot[k],j) ];

         LU[ ENTRY(pivot[i],k) ] = m;
      }
   }

   /* second step: forward substitution */
   y[0] = b[pivot[0]];

   for( k = 1; k < N; ++k )
   {
      SCIP_Real s;
      int j;

      s = b[pivot[k]];
      for( j = 0; j < k; ++j )
      {
         s -= LU[ ENTRY(pivot[k],j) ] * y[j];
      }
      y[k] = s;
   }

   /* third step: backward substitution */
   x[N-1] = y[N-1] / LU[ ENTRY(pivot[N-1],N-1) ];
   for( k = N-2; k >= 0; --k )
   {
      SCIP_Real s;
      int j;

      s = y[k];
      for( j = k+1; j < N; ++j )
      {
         s -= LU[ ENTRY(pivot[k],j) ] * x[j];
      }
      x[k] = s / LU[ ENTRY(pivot[k],k) ];
   }

   return SCIP_OKAY;
}

/* solves a linear problem of the form Ax = b for a regular matrix A */
SCIP_RETCODE SCIPsolveLinearEquationsIpopt(
   int                   N,                  /**< dimension */
   SCIP_Real*            A,                  /**< matrix data on input (size N*N); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size N) */
   SCIP_Real*            x,                  /**< buffer to store solution (size N) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   )
{
   SCIP_Real* LU;
   SCIP_Real* y;
   int* pivot;
   int k;

   SCIP_RETCODE retcode = SCIP_OKAY;

   assert(N > 0);
   assert(A != NULL);
   assert(b != NULL);
   assert(x != NULL);
   assert(success != NULL);

   /* call solveLinearProb3() for performance reasons */
   if( N == 3 )
   {
      SCIP_CALL( solveLinearProb3(A, b, x, success) );
      return SCIP_OKAY;
   }

   *success = TRUE;

   LU = NULL;
   y = NULL;
   pivot = NULL;

   /* copy arrays */
   SCIP_ALLOC_TERMINATE( retcode, BMSduplicateMemoryArray(&LU, A, N*N), TERMINATE ); /*lint !e647*/
   SCIP_ALLOC_TERMINATE( retcode, BMSduplicateMemoryArray(&y, b, N), TERMINATE );
   SCIP_ALLOC_TERMINATE( retcode, BMSallocMemoryArray(&pivot, N), TERMINATE );

   /* initialize values */
   for( k = 0; k < N; ++k )
      pivot[k] = k;

   /* first step: compute LU factorization */
   for( k = 0; k < N; ++k )
   {
      int p;
      int i;

      p = k;
      for( i = k+1; i < N; ++i )
      {
         if( ABS(LU[ ENTRY(pivot[i],k) ]) > ABS( LU[ ENTRY(pivot[p],k) ]) )
            p = i;
      }

      if( ABS(LU[ ENTRY(pivot[p],k) ]) < 1e-08 )
      {
         /* SCIPerrorMessage("Error in nlpi_ipopt_dummy - matrix is singular!\n"); */
         *success = FALSE;
         goto TERMINATE;
      }

      if( p != k )
      {
         int tmp;

         tmp = pivot[k];
         pivot[k] = pivot[p];
         pivot[p] = tmp;
      }

      for( i = k+1; i < N; ++i )
      {
         SCIP_Real m;
         int j;

         m = LU[ ENTRY(pivot[i],k) ] / LU[ ENTRY(pivot[k],k) ];

         for( j = k+1; j < N; ++j )
            LU[ ENTRY(pivot[i],j) ] -= m * LU[ ENTRY(pivot[k],j) ];

         LU[ ENTRY(pivot[i],k) ] = m;
      }
   }

   /* second step: forward substitution */
   y[0] = b[pivot[0]];

   for( k = 1; k < N; ++k )
   {
      SCIP_Real s;
      int j;

      s = b[pivot[k]];
      for( j = 0; j < k; ++j )
      {
         s -= LU[ ENTRY(pivot[k],j) ] * y[j];
      }
      y[k] = s;
   }

   /* third step: backward substitution */
   x[N-1] = y[N-1] / LU[ ENTRY(pivot[N-1],N-1) ];
   for( k = N-2; k >= 0; --k )
   {
      SCIP_Real s;
      int j;

      s = y[k];
      for( j = k+1; j < N; ++j )
      {
         s -= LU[ ENTRY(pivot[k],j) ] * x[j];
      }
      x[k] = s / LU[ ENTRY(pivot[k],k) ];
   }

   TERMINATE:
   /* free arrays */
   BMSfreeMemoryArrayNull(&pivot);
   BMSfreeMemoryArrayNull(&y);
   BMSfreeMemoryArrayNull(&LU);

   return retcode;
}
