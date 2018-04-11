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

/**@file   concsolver_scip.h
 * @ingroup PARALLEL
 * @brief  implementation of concurrent solver interface for SCIP
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONCSOLVER_SCIP_H__
#define __SCIP_CONCSOLVER_SCIP_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the concurrent SCIP solver plugins and includes them in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeConcurrentScipSolvers(
   SCIP*                 scip                /**< SCIP datastructure */
   );

#ifdef __cplusplus
}
#endif

#endif
