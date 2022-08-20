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

/**@file   scipcoreplugins.h
 * @ingroup INTERNALAPI
 * @brief  register additional core functionality that is designed as plugins
 * @author Gregor Hendel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIPCOREPLUGINS_H__
#define __SCIP_SCIPCOREPLUGINS_H__

#include "scip/scip.h"

/* include header files here, so that the SCIP core only needs to include the SCIP core plugins
 */
#include "scip/bandit_epsgreedy.h"
#include "scip/bandit_exp3.h"
#include "scip/bandit_ucb.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes core SCIP plugins into SCIP */
extern
SCIP_RETCODE SCIPincludeCorePlugins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
