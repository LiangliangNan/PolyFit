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

/**@file   type_cutsel.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for cut selectors
 * @author Felipe Serrano
 * @author Mark Turner
 */

/** @defgroup DEFPLUGINS_CUTSEL Default cut selectors
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default cut selectors of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CUTSEL_H__
#define __SCIP_TYPE_CUTSEL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_lp.h"
#include "scip/type_result.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Cutsel SCIP_CUTSEL;         /**< cut selector data structure */
typedef struct SCIP_CutselData SCIP_CUTSELDATA; /**< cut selector specific data */


/** copy method for cut selector plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cutsel          : the cut selector itself
 */
#define SCIP_DECL_CUTSELCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_CUTSEL* cutsel)

/** destructor of cut selector to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cutsel          : the cut selector itself
 */
#define SCIP_DECL_CUTSELFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_CUTSEL* cutsel)

/** initialization method of cut selector (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cutsel          : the cut selector itself
 */
#define SCIP_DECL_CUTSELINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_CUTSEL* cutsel)

/** deinitialization method of cut selector (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cutsel          : the cut selector itself
 */
#define SCIP_DECL_CUTSELEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_CUTSEL* cutsel)

/** solving process initialization method of cut selector (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The cut selector may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cutsel          : the cut selector itself
 */
#define SCIP_DECL_CUTSELINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_CUTSEL* cutsel)

/** solving process deinitialization method of cut selector (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The cut selector should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cutsel          : the cut selector itself
 */
#define SCIP_DECL_CUTSELEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_CUTSEL* cutsel)

/** cut selection method of cut selector
 *
 *  This method is called to select the cuts to be added to the LP.
 *  Forcedcuts must not be changed, and cuts should only be resorted, with the first nselectedcuts of cuts being chosen.
 *  These nselectededcuts are used in addition to the forcedcuts (do not delete nor modify elements, simply resort).
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cutsel          : the cut selector itself
 *  - cuts            : cutting planes to select from
 *  - ncuts           : number of cutting planes to select from (length of cuts)
 *  - forcedcuts      : list of cuts that are forced to be applied (i.e they are going to be selected no matter what)
 *  - nforcedcuts     : number of forced cuts
 *  - root            : are we at the root node?
 *  - maxselectedcuts : maximum number of cuts that can be selected (upper bound for nselectedcuts)
 *  - nselectedcuts   : the first nselectedcuts from cuts are selected
 *  - result          : pointer to store the result of the cut selection call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_SUCCESS    : the cut selection succeeded
 *  - SCIP_DIDNOTFIND : the cut selection did not find good enough cuts to select
 */
#define SCIP_DECL_CUTSELSELECT(x) SCIP_RETCODE x (SCIP* scip, SCIP_CUTSEL* cutsel, SCIP_ROW** cuts, int ncuts, \
      SCIP_ROW** forcedcuts, int nforcedcuts, SCIP_Bool root, int maxnselectedcuts, int* nselectedcuts, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
