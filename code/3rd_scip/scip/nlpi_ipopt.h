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

/**@file    nlpi_ipopt.h
 * @brief   Ipopt NLP interface
 * @ingroup NLPIS
 * @author  Stefan Vigerske
 * @author  Benjamin Müller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_IPOPT_H__
#define __SCIP_NLPI_IPOPT_H__

#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for Ipopt solver and includes it into SCIP, if Ipopt is available
 *
 * @ingroup NLPIIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlpSolverIpopt(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLPIS
 *
 * @{
 */

/** gets string that identifies Ipopt (version number) */
SCIP_EXPORT
const char* SCIPgetSolverNameIpopt(void);

/** gets string that describes Ipopt */
SCIP_EXPORT
const char* SCIPgetSolverDescIpopt(void);

/** returns whether Ipopt is available, i.e., whether it has been linked in */
SCIP_EXPORT
SCIP_Bool SCIPisIpoptAvailableIpopt(void);

/** gives a pointer to the NLPIORACLE object stored in Ipopt-NLPI's NLPI problem data structure */
SCIP_EXPORT
void* SCIPgetNlpiOracleIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   );

/** Calls Lapacks Dsyev routine to compute eigenvalues and eigenvectors of a dense matrix.
 *
 * It's here, because we use Ipopt's C interface to Lapack.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcallLapackDsyevIpopt(
   SCIP_Bool             computeeigenvectors,/**< should also eigenvectors should be computed ? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if computeeigenvectors == TRUE */
   SCIP_Real*            w                   /**< buffer to store eigenvalues (size N) */
   );

/** solves a linear problem of the form Ax = b for a regular matrix A
 *
 *  Calls Lapacks DGETRF routine to calculate a LU factorization and uses this factorization to solve
 *  the linear problem Ax = b.
 *
 * It's here, because we use Ipopt's C interface to Lapack.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsolveLinearEquationsIpopt(
   int                   N,                  /**< dimension */
   SCIP_Real*            A,                  /**< matrix data on input (size N*N); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size N) */
   SCIP_Real*            x,                  /**< buffer to store solution (size N) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_IPOPT_H__ */
