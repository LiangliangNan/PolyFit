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

/**@file   type_sepa.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SEPA_H__
#define __SCIP_TYPE_SEPA_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_sol.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Sepa SCIP_SEPA;               /**< separator */
typedef struct SCIP_SepaData SCIP_SEPADATA;       /**< locally defined separator data */


/** copy method for separator plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define SCIP_DECL_SEPACOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa)

/** destructor of separator to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define SCIP_DECL_SEPAFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa)

/** initialization method of separator (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define SCIP_DECL_SEPAINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa)

/** deinitialization method of separator (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define SCIP_DECL_SEPAEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa)

/** solving process initialization method of separator (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The separator may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define SCIP_DECL_SEPAINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa)

/** solving process deinitialization method of separator (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The separator should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define SCIP_DECL_SEPAEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa)

/** LP solution separation method of separator
 *
 *  Searches for cutting planes that separate the current LP solution. The method is called in the LP solving loop,
 *  which means that a valid LP solution exists.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 *  - result          : pointer to store the result of the separation call
 *  - allowlocal      : should the separator allow local cuts?
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_NEWROUND   : a cutting plane was generated and a new separation round should immediately start
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
#define SCIP_DECL_SEPAEXECLP(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa, SCIP_RESULT* result, SCIP_Bool allowlocal)

/** arbitrary primal solution separation method of separator
 *
 *  Searches for cutting planes that separate the given primal solution. The method is called outside the LP solution
 *  loop (e.g., by a relaxator or a primal heuristic), which means that there is no valid LP solution.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 *  - sol             : primal solution that should be separated
 *  - result          : pointer to store the result of the separation call
 *  - allowlocal      : should the separator allow local cuts?
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_NEWROUND   : a cutting plane was generated and a new separation round should immediately start
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
#define SCIP_DECL_SEPAEXECSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_SEPA* sepa, SCIP_SOL* sol, SCIP_RESULT* result, SCIP_Bool allowlocal)

#ifdef __cplusplus
}
#endif

#endif
