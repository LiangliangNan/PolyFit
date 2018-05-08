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

/**@file   prop_rootredcost.h
 * @ingroup PROPAGATORS
 * @brief  reduced cost strengthening using root node reduced costs and the cutoff bound
 * @author Tobias Achterberg
 * @author Stefan Heinz
 *
 * This propagator uses the root reduced cost to (globally) propagate against the cutoff bound. The propagator checks if
 * the variables with non-zero root reduced cost can exceed the cutoff bound. If this is the case the corresponding
 * bound can be tightened.
 *
 * The propagate is performed during the search any time a new cutoff bound (primal solution) is found.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_ROOTREDCOST_H__
#define __SCIP_PROP_ROOTREDCOST_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the root node reduced cost strengthening propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePropRootredcost(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
