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

/**@file   presol_dualcomp.h
 * @ingroup PRESOLVERS
 * @brief  dual compensation presolver
 * @author Dieter Weninger
 *
 * This presolver looks for variables with
 *         i) objcoef >= 0 and exactly one downlock
 *        ii) objcoef <= 0 and exactly one uplock
 * and fixes the variable in case i) at the lower bound and in case ii) at the
 * upper bound if a combination of singleton continuous variables can compensate
 * the downlock in case i) and the uplock in case ii).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_DUALCOMP_H__
#define __SCIP_PRESOL_DUALCOMP_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dualcomp presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePresolDualcomp(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
