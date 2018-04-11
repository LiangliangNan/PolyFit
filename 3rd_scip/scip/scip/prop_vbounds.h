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

/**@file   prop_vbounds.h
 * @ingroup PROPAGATORS
 * @brief  variable upper and lower bound propagator
 * @author Stefan Heinz
 * @author Jens Schulz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_VBOUNDS_H__
#define __SCIP_PROP_VBOUNDS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the vbounds propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePropVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PROPAGATORS
  *
  * @{
  */

/** returns TRUE if the propagator has the status that all variable lower and upper bounds are propagated */
EXTERN
SCIP_Bool SCIPisPropagatedVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** performs propagation of variables lower and upper bounds */
EXTERN
SCIP_RETCODE SCIPexecPropVbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_RESULT*          result              /**< pointer to store result */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
