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

/**@file   table_xyz.h
 * @ingroup TABLES
 * @brief  xyz statistics table
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TABLE_XYZ_H__
#define __SCIP_TABLE_XYZ_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the xyz statistics table and includes it in SCIP
 *
 *  @ingroup TableIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeTableXyz(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup TABLES
 *
 * @{
 */

/** TODO add further methods to this group for documentation purposes */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
