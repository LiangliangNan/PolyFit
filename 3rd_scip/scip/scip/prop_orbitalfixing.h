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

/**@file   prop_orbitalfixing.h
 * @brief  propagator for orbital fixing
 * @author Marc Pfetsch
 *
 * This propagator implements orbital fixing, a method to exploit symmetries to fix variables in orbits to 0 or 1.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_ORBITALFIXING_H__
#define __SCIP_PROP_ORBITALFIXING_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** creates the orbitalfixing propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePropOrbitalfixing(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
