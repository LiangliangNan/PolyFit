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

/**@file   disp_xyz.h
 * @ingroup DISPLAYS
 * @brief  xyz display column
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DISP_XYZ_H__
#define __SCIP_DISP_XYZ_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the xyz display column and includes it in SCIP
 *
 *  @ingroup DisplayIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeDispXyz(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup DISPLAYS
 *
 * @{
 */

/** TODO add further methods to this group for documentation purposes */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
