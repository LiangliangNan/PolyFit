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

/**@file   retcode.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RETCODE_H__
#define __SCIP_RETCODE_H__

#include <stdio.h>

#include "scip/pub_message.h"

#ifdef __cplusplus
extern "C" {
#endif

/** prints error message for return code */
extern
void SCIPretcodePrint(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to write error message */
   SCIP_RETCODE          retcode             /**< SCIP return code causing the error */
   );

/** prints error message for return code via error message */
extern
void SCIPretcodePrintError(
   SCIP_RETCODE          retcode             /**< SCIP return code causing the error */
   );

#ifdef __cplusplus
}
#endif

#endif
