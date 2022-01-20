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

/**@file   message_default.h
 * @ingroup PUBLICMETHODS
 * @brief  default message handler
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MESSAGE_DEFAULT_H__
#define __SCIP_MESSAGE_DEFAULT_H__

#include "scip/def.h"
#include "scip/type_message.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Create default message handler. To free the message handler use SCIPmessagehdlrRelease(). */
EXTERN
SCIP_RETCODE SCIPcreateMessagehdlrDefault(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   const char*           filename,           /**< name of log file, or NULL (stdout) */
   SCIP_Bool             quiet               /**< should screen messages be suppressed? */
   );

#ifdef __cplusplus
}
#endif

#endif
