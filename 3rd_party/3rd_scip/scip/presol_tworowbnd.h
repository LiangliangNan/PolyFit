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

/**@file   presol_tworowbnd.h
 * @ingroup PRESOLVERS
 * @brief  do bound tightening by using two rows
 * @author Dieter Weninger
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_TWOROWBND_H__
#define __SCIP_PRESOL_TWOROWBND_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the tworowbnd presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePresolTworowbnd(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
