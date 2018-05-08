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

/**@file nodesel_breadthfirst.h
 * @ingroup NODESELECTORS
 * @brief node selector for breadth-first search
 * @author Stefan Heinz
 * @author Gregor Hendel
 *
 * This node selector performs breadth-first search, i.e., it completely evaluates an entire level of the search tree before
 * proceeding to the next level. At one level, nodes are processed in the order they were created by the branching rule.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NODESEL_BREADTHFIRST_H__
#define __SCIP_NODESEL_BREADTHFIRST_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the node selector for breadth first search and includes it in SCIP
 *
 *  @ingroup NodeSelectorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeNodeselBreadthfirst(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
