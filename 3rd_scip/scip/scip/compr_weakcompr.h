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

/**@file   compr_weakcompr.h
 * @ingroup COMPRESSION
 * @brief  weakcompr tree compression
 * @author Jakob Witzig
 *
 * template file for tree compression plugins
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_COMPR_WEAKCOMPR_H__
#define __SCIP_COMPR_WEAKCOMPR_H__


#include "scip/scip.h"
#include "scip/type_reopt.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the weakcompr tree compression and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeComprWeakcompr(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
