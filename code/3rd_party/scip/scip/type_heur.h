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

/**@file   type_heur.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for primal heuristics
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 *  This file defines the interface for primal heuristics implemented in C.
 *
 *  - \ref HEUR "Instructions for implementing a primal heuristic"
 *  - \ref PRIMALHEURISTICS "List of available primal heuristics"
 *  - \ref scip::ObjHeur "C++ wrapper class"
 */

/** @defgroup DEFPLUGINS_HEUR Default Primal Heuristics
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default primal heuristics of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_HEUR_H__
#define __SCIP_TYPE_HEUR_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "scip/type_timing.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** represents different methods for a dive set to explore the next children */
#define SCIP_DIVETYPE_NONE                0x000u  /**< no method specified */
#define SCIP_DIVETYPE_INTEGRALITY         0x001u  /**< use branching on a variable by shrinking the domain in the child nodes */
#define SCIP_DIVETYPE_SOS1VARIABLE        0x002u  /**< branch on a variable solution value by exploiting special-ordered set conflict structure */

typedef unsigned int SCIP_DIVETYPE;

/** context for diving statistics */
enum SCIP_DiveContext
{
   SCIP_DIVECONTEXT_TOTAL  = 0,                   /**< all contexts combined */
   SCIP_DIVECONTEXT_SINGLE = 1,                   /**< single heuristic context */
   SCIP_DIVECONTEXT_ADAPTIVE = 2                  /**< within adaptive diving */
};
typedef enum SCIP_DiveContext SCIP_DIVECONTEXT;


typedef struct SCIP_Heur SCIP_HEUR;               /**< primal heuristic */
typedef struct SCIP_HeurData SCIP_HEURDATA;       /**< locally defined primal heuristic data */
typedef struct SCIP_Diveset SCIP_DIVESET;         /**< common parameters for all diving heuristics */
typedef struct SCIP_VGraph SCIP_VGRAPH;           /**< variable graph data structure to determine breadth-first
                                                    *  distances between variables */

/** commonly used display characters indicating special classes of primal heuristics */
#define SCIP_HEURDISPCHAR_LNS       'L'  /**< a 'L'arge Neighborhood or other local search heuristic */
#define SCIP_HEURDISPCHAR_DIVING    'd'  /**< a 'd'iving heuristic that dives down an auxiliary branch-and-bound path */
#define SCIP_HEURDISPCHAR_ITERATIVE 'i'  /**< an iterative improvement heuristic such as 1-opt or 2-opt */
#define SCIP_HEURDISPCHAR_OBJDIVING 'o'  /**< an 'o'bjective diving or feasibility pump heuristic */
#define SCIP_HEURDISPCHAR_PROP      'p'  /**< a 'p'ropagation heuristic, often applied before branch-and-bound starts */
#define SCIP_HEURDISPCHAR_ROUNDING  'r'  /**< a 'r'ounding heuristic that iteratively tries to round an LP or relaxation solution */
#define SCIP_HEURDISPCHAR_TRIVIAL   't'  /**< a 't'rivial or helper heuristic, usually applied before branch-and-bound starts */

/** copy method for heuristic plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** destructor of primal heuristic to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** initialization method of primal heuristic (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** deinitialization method of primal heuristic (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEUREXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEUREXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** execution method of primal heuristic
 *
 *  Searches for feasible primal solutions. The method is called in the node processing loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 *  - heurtiming      : current point in the node solving loop
 *  - nodeinfeasible  : was the current node already detected to be infeasible?
 *  - result          : pointer to store the result of the heuristic call
 *
 *  possible return values for *result:
 *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
 *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
 *  - SCIP_DIDNOTRUN  : the heuristic was skipped
 *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
 *                      its frequency
 */
#define SCIP_DECL_HEUREXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur, SCIP_HEURTIMING heurtiming, \
      SCIP_Bool nodeinfeasible, SCIP_RESULT* result)


/* callbacks for diving heuristic specific settings */

/** calculate score and preferred rounding direction for the candidate variable; the best candidate maximizes the
 *  score
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - diveset         : diving settings for scoring
 *  - divetype        : represents different methods for a dive set to explore the next children
 *  - cand            : candidate variable for which the score should be determined
 *  - candsol         : solution value of variable in LP relaxation solution
 *  - candsfrac       : fractional part of solution value of variable
 *  - score           : pointer for diving score value - the best candidate maximizes this score
 *  - roundup         : pointer to store whether the preferred rounding direction is upwards
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_DIVESETGETSCORE(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIVESET* diveset, \
   SCIP_DIVETYPE divetype, SCIP_VAR* cand, SCIP_Real candsol, SCIP_Real candsfrac, SCIP_Real* score, SCIP_Bool* roundup)

/**
 * optional callback to check preconditions for diving, e.g., if an incumbent solution is available
 *
 * input:
 *  - scip            : SCIP main data structure
 *  - diveset         : diving settings for scoring
 *
 * output:
 *  - available       : TRUE if diveset can run, otherwise FALSE
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_DIVESETAVAILABLE(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIVESET* diveset, SCIP_Bool* available)

#ifdef __cplusplus
}
#endif

#endif
