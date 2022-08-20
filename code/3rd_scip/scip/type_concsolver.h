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

/**@file   type_concsolver.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for concurrent solvers
 * @author Robert Lion Gottwald
 *
 *  This file defines the interface for concurrent solvers.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CONCSOLVER_H__
#define __SCIP_TYPE_CONCSOLVER_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_syncstore.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_ConcSolverType SCIP_CONCSOLVERTYPE;         /**< the struct defining a concurrent solver class */
typedef struct SCIP_ConcSolverTypeData SCIP_CONCSOLVERTYPEDATA; /**< concurrent solver class user data */
typedef struct SCIP_ConcSolver SCIP_CONCSOLVER;                   /**< struct for an instance of a concurrent solver */
typedef struct SCIP_ConcSolverData SCIP_CONCSOLVERDATA;           /**< concurrent solver user data */

/** creates a concurrent solver instance
 *
 *  input:
 *  - scip               : SCIP main data structure
 *  - concsolvertype     : type of concurrent solver an instance should be created for
 *  - concsolverinstance : pointer to return concurrent solver instance
 *
 * returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERCREATEINST(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONCSOLVERTYPE* concsolvertype, SCIP_CONCSOLVER* concsolver)

/** destroys a concurrent solver instance
 *
 *  input:
 *  - scip               : SCIP main data structure
 *  - concsolverinstance : concurrent solver instance to destroy
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERDESTROYINST(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONCSOLVER* concsolver)

/** frees data of a concurrent solver type
 *
 *  input:
 *  - scip               : SCIP main data structure
 *  - data               : concurrent solver type data to free
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERTYPEFREEDATA(x) void x (SCIP_CONCSOLVERTYPEDATA** data)

/** initialize random seeds of a concurrent solver
 *
 *  input:
 *  - concsolver      : concurrent solver data structure
 *  - seed            : seed for initializing the solver's internal random seeds
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERINITSEEDS(x) SCIP_RETCODE x (SCIP_CONCSOLVER* concsolver, unsigned int seed)

/** synchronization method of concurrent solver for writing data
 *
 *  Syncronizes with other solvers. The concurrent solver should pass new solutions
 *  and bounds to the syncstore. For the solutions, no more than maxcandsols of the best solution
 *  should be considered for sharing. Additionally a maximum if maxsharedsols should be
 *  passed to the syncstore.
 *
 *  input:
 *  - concsolver      : concurrent solver data structure
 *  - spi             : pointer to the SCIP parallel interface
 *  - syncdata        : concurrent solver data structure
 *  - maxcandsols     : how many of the best solutions should be considered for sharing
 *  - maxsharedsols   : the maximum number of solutions that should be shared
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERSYNCWRITE(x) SCIP_RETCODE x (SCIP_CONCSOLVER* concsolver, SCIP_SYNCSTORE* syncstore, SCIP_SYNCDATA* syncdata, int maxcandsols, int maxsharedsols, int* nsolsshared)

/** synchronization method of concurrent solver for reading data
 *
 *  the concurrent solver should read the solutions and bounds stored in the
 *  given synchronization data
 *
 *  input:
 *  - concsolver      : concurrent solver data structure
 *  - spi             : pointer to the SCIP parallel interface
 *  - syncdata        : concurrent solver data structure
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERSYNCREAD(x) SCIP_RETCODE x (SCIP_CONCSOLVER* concsolver, SCIP_SYNCSTORE* syncstore, SCIP_SYNCDATA* syncdata, int* nsolsrecvd, int* ntighterbnds, int* ntighterintbnds)

/** execution method of concurrent solver
 *
 *  start solving of the problem given during initialization
 *
 *  input:
 *  - concsolver       : concurrent solver data structure
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVEREXEC(x) SCIP_RETCODE x (SCIP_CONCSOLVER* concsolver, SCIP_Real* solvingtime, SCIP_Longint* nlpiterations, SCIP_Longint* nnodes)

/** stop the solving as soon as possible
 *
 *  input:
 *  - concsolver      : concurrent solver data structure
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERSTOP(x) SCIP_RETCODE x (SCIP_CONCSOLVER* concsolver)

/** extract the solving data from the concurrent solver and store it into the SCIP datastructure,
 *  so that this SCIP instance has the optimal solution and reports the correct status and statistics.
 *
 *  input:
 *  - concsolver      : concurrent solver data structure
 *  - scip            : SCIP datastructure
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA(x) SCIP_RETCODE x (SCIP_CONCSOLVER* concsolver, SCIP* scip)


#ifdef __cplusplus
}
#endif

#endif
