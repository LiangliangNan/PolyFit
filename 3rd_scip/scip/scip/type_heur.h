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

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_HEUR_H__
#define __SCIP_TYPE_HEUR_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "scip/type_timing.h"

#ifdef __cplusplus
extern "C" {
#endif

/** represents different methods for a dive set to explore the next children */
#define SCIP_DIVETYPE_NONE                0x000u  /**< no method specified */
#define SCIP_DIVETYPE_INTEGRALITY         0x001u  /**< use branching on a variable by shrinking the domain in the child nodes */
#define SCIP_DIVETYPE_SOS1VARIABLE        0x002u  /**< branch on a variable solution value by exploiting special-ordered set conflict structure */

typedef unsigned int SCIP_DIVETYPE;

typedef struct SCIP_Heur SCIP_HEUR;               /**< primal heuristic */
typedef struct SCIP_HeurData SCIP_HEURDATA;       /**< locally defined primal heuristic data */
typedef struct SCIP_Diveset SCIP_DIVESET;         /**< common parameters for all diving heuristics */
typedef struct SCIP_VGraph SCIP_VGRAPH;           /**< variable graph data structure to determine breadth-first
                                                    *  distances between variables */

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
 *  - cand            : Candidate variable for which the score should be determined
 *  - candsol         : solution value of variable in LP relaxation solution
 *  - candsfrac       : fractional part of solution value of variable
 *  - score           : pointer for diving score value - small scores are preferred
 *  - roundup         : pointer to store whether the preferred rounding direction is upwards
 *
 *  returns SCIP_OKAY if everything worked, otherwise, a suitable error code
 */
#define SCIP_DECL_DIVESETGETSCORE(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIVESET* diveset, \
   SCIP_DIVETYPE divetype, SCIP_VAR* cand, SCIP_Real candsol, SCIP_Real candsfrac, SCIP_Real* score, SCIP_Bool* roundup)

#ifdef __cplusplus
}
#endif

#endif
