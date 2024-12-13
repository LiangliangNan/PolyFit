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

/**@file   type_relax.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_RELAX_H__
#define __SCIP_TYPE_RELAX_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Relax SCIP_RELAX;             /**< relaxator */
typedef struct SCIP_Relaxation SCIP_RELAXATION;   /**< relaxator */
typedef struct SCIP_RelaxData SCIP_RELAXDATA;     /**< locally defined relaxator data */


/** copy method for relaxator plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define SCIP_DECL_RELAXCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax)

/** destructor of relaxator to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define SCIP_DECL_RELAXFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax)

/** initialization method of relaxator (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define SCIP_DECL_RELAXINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax)

/** deinitialization method of relaxator (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define SCIP_DECL_RELAXEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax)

/** solving process initialization method of relaxator (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The relaxator may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define SCIP_DECL_RELAXINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax)

/** solving process deinitialization method of relaxator (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The relaxator should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define SCIP_DECL_RELAXEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax)

/** execution method of relaxator
 *
 *  The method is called in the node processing loop. It solves the current subproblem's relaxation.
 *  Like the LP relaxation, the relaxator should only operate on COLUMN variables.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 *  - lowerbound      : pointer to store a lowerbound for the current node
 *  - result          : pointer to store the result of the relaxation call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated, and the relaxator should not be called again on the
 *                      same relaxation
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced, and the relaxator should not be called again on the same
 *                      relaxation
 *  - SCIP_SEPARATED  : a cutting plane was generated, and the relaxator should not be called again on the same relaxation
 *  - SCIP_SUCCESS    : the relaxator solved the relaxation and should not be called again on the same relaxation
 *  - SCIP_SUSPENDED  : the relaxator interrupted its solving process to wait for additional input (e.g. cutting
 *                      planes); however, it is able to continue the solving in order to improve the dual bound
 *  - SCIP_DIDNOTRUN  : the relaxator was skipped
 */
#define SCIP_DECL_RELAXEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_RELAX* relax, SCIP_Real* lowerbound, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
