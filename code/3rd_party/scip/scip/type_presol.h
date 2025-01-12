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

/**@file   type_presol.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for presolvers
 * @author Tobias Achterberg
 */

/** @defgroup DEFPLUGINS_PRESOL Default Presolvers
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default presolvers of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_PRESOL_H__
#define __SCIP_TYPE_PRESOL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Presol SCIP_PRESOL;           /**< presolver data structure */
typedef struct SCIP_PresolData SCIP_PRESOLDATA;   /**< presolver specific data */


/** copy method for presolver plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** destructor of presolver to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** initialization method of presolver (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** deinitialization method of presolver (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** presolving initialization method of presolver (called when presolving is about to begin)
 *
 *  This method is called when the presolving process is about to begin, even if presolving is turned off.
 *  The presolver may use this call to initialize its data structures.
 *
 *  Necessary modifications that have to be performed even if presolving is turned off should be done here or in the
 *  presolving deinitialization call (SCIP_DECL_PRESOLSEXITPRE()).
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLINITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** presolving deinitialization method of presolver (called after presolving has been finished)
 *
 *  This method is called after the presolving has been finished, even if presolving is turned off.
 *  The presolver may use this call e.g. to clean up or modify its data structures.
 *
 *  Necessary modifications that have to be performed even if presolving is turned off should be done here or in the
 *  presolving initialization call (SCIP_DECL_PRESOLINITPRE()).
 *
 *  Besides necessary modifications and clean up, no time consuming operations should be performed, especially if the
 *  problem has already been solved.  Use the method SCIPgetStatus(), which in this case returns SCIP_STATUS_OPTIMAL,
 *  SCIP_STATUS_INFEASIBLE, SCIP_STATUS_UNBOUNDED, or SCIP_STATUS_INFORUNBD.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 */
#define SCIP_DECL_PRESOLEXITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol)

/** execution method of presolver
 *
 *  The presolver should go through the variables and constraints and tighten the domains or
 *  constraints. Each tightening should increase the given total numbers of changes.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - presol          : the presolver itself
 *  - nrounds         : number of presolving rounds already done
 *  - presoltiming    : current presolving timing
 *  - nnewfixedvars   : number of variables fixed since the last call to the presolver
 *  - nnewaggrvars    : number of variables aggregated since the last call to the presolver
 *  - nnewchgvartypes : number of variable type changes since the last call to the presolver
 *  - nnewchgbds      : number of variable bounds tightened since the last call to the presolver
 *  - nnewholes       : number of domain holes added since the last call to the presolver
 *  - nnewdelconss    : number of deleted constraints since the last call to the presolver
 *  - nnewaddconss    : number of added constraints since the last call to the presolver
 *  - nnewupgdconss   : number of upgraded constraints since the last call to the presolver
 *  - nnewchgcoefs    : number of changed coefficients since the last call to the presolver
 *  - nnewchgsides    : number of changed left or right hand sides since the last call to the presolver
 *
 *  @note the counters state the changes since the last call including the changes of this presolver during its last
 *        last call
 *
 *  @note if the presolver uses dual information it is nesassary to check via calling SCIPallowWeakDualReds and
 *        SCIPallowStrongDualReds if dual reductions are allowed.
 *
 *  input/output:
 *  - nfixedvars      : pointer to total number of variables fixed of all presolvers
 *  - naggrvars       : pointer to total number of variables aggregated of all presolvers
 *  - nchgvartypes    : pointer to total number of variable type changes of all presolvers
 *  - nchgbds         : pointer to total number of variable bounds tightened of all presolvers
 *  - naddholes       : pointer to total number of domain holes added of all presolvers
 *  - ndelconss       : pointer to total number of deleted constraints of all presolvers
 *  - naddconss       : pointer to total number of added constraints of all presolvers
 *  - nupgdconss      : pointer to total number of upgraded constraints of all presolvers
 *  - nchgcoefs       : pointer to total number of changed coefficients of all presolvers
 *  - nchgsides       : pointer to total number of changed left/right hand sides of all presolvers
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_SUCCESS    : the presolver found a reduction
 *  - SCIP_DIDNOTFIND : the presolver searched, but did not find a presolving change
 *  - SCIP_DIDNOTRUN  : the presolver was skipped
 */
#define SCIP_DECL_PRESOLEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_PRESOL* presol, int nrounds, SCIP_PRESOLTIMING presoltiming, \
      int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
      int nnewdelconss, int nnewaddconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides, \
      int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes, \
      int* ndelconss, int* naddconss, int* nupgdconss, int* nchgcoefs, int* nchgsides, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
