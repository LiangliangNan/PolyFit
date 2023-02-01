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

/**@file   presol_convertinttobin.h
 * @ingroup PRESOLVERS
 * @brief  presolver that converts integer variables with domain [a,a+1] to binaries
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_CONVERTINTTOBIN_H__
#define __SCIP_PRESOL_CONVERTINTTOBIN_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the convertinttobin presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePresolConvertinttobin(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
