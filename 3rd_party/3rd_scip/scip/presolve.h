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

/**@file   presolve.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods commonly used for presolving
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOLVE_H__
#define __SCIP_PRESOLVE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@defgroup PublicSpecialPresolveMethods Special Methods
 * @ingroup PublicPresolverMethods
 * @brief methods commonly used for presolving
 *
 * @{
 */
/** try to reduce the necessary variable in a set of variables with corresponding bounds and boundtypes for which one
 *  must be fulfilled
 *
 *  e.g. a set of logicor or bounddisjunctive constraint variables would be such a set
 *
 *  consider the following set:
 *
 *  x1 >= 1, x2 >= 3, x3 >= 1, x4 <= 0
 *
 *  by (global) implication data (cliques, implications, and variable bounds) we have also the following implications
 *  given:
 *
 *  x1 >= 1 => x3 >= 1
 *  x2 >= 2 => x3 >= 1
 *  x4 <= 0 => x1 >= 1
 *
 *  Because of the last implication x4 is redundant, because x1 >= 1 would also be fulfilled in the variable set, so we
 *  can reduce the set by x4.
 *  Also, the both other implications and x3 >= 1 (in the given variable set) all imply exactly x3 >= 1, so we tighten
 *  the global lower bound of x3 to 1 and the set of variables gets redundant.
 */
EXTERN
SCIP_RETCODE SCIPshrinkDisjunctiveVarSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables array for which at least one must be fulfilled in the
                                              *   following bounds and boundtypes */
   SCIP_Real*            bounds,             /**< bounds array for which at least one must be fulfilled */
   SCIP_Bool*            boundtypes,         /**< boundtypes array (TRUE == SCIP_BOUNDTYPE_UPPER, FALSE == SCIP_BOUNDTYPE_LOWER)
                                              *   for which at least one must be fulfilled */
   SCIP_Bool*            redundants,         /**< array which be filled and then indicate if a variable in the set is redundant */
   int                   nvars,              /**< number of variables */
   int*                  nredvars,           /**< pointer to store how many variables can be removed */
   int*                  nglobalred,         /**< pointer to store number of global reductions on variable bounds found
                                              *   through this set of variables */
   SCIP_Bool*            setredundant,       /**< pointer to store if we found a global reduction on a variable which was part
                                              *   of the given set of variables, this makes this disjunction redundant */
   SCIP_Bool*            glbinfeas,          /**< pointer to store if global infeasibility was detected */
   SCIP_Bool             fullshortening      /**< do we want to try the shortening procedure over the whole set (which might be expensive) */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
