/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_message.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for message output methods
 * @author Tobias Achterberg
 *
 *  This file defines the interface for message handlers implemented in C.
 *
 *  - \ref scip::ObjMessagehdlr "C++ wrapper class"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_MESSAGE_H__
#define __SCIP_TYPE_MESSAGE_H__


#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** verbosity levels of output */
enum SCIP_VerbLevel
{
   SCIP_VERBLEVEL_NONE    = 0,          /**< only error and warning messages are displayed */
   SCIP_VERBLEVEL_DIALOG  = 1,          /**< only interactive dialogs, errors, and warnings are displayed */
   SCIP_VERBLEVEL_MINIMAL = 2,          /**< only important messages are displayed */
   SCIP_VERBLEVEL_NORMAL  = 3,          /**< standard messages are displayed */
   SCIP_VERBLEVEL_HIGH    = 4,          /**< a lot of information is displayed */
   SCIP_VERBLEVEL_FULL    = 5           /**< all messages are displayed */
};
typedef enum SCIP_VerbLevel SCIP_VERBLEVEL;

typedef struct SCIP_Messagehdlr SCIP_MESSAGEHDLR;           /**< message handler */
typedef struct SCIP_MessagehdlrData SCIP_MESSAGEHDLRDATA;   /**< message handler data */

/** generic messagehandler output function
 *
 *  Should be equal to SCIP_DECL_MESSAGEWARNING, SCIP_DECL_MESSAGEDIALOG, and SCIP_DECL_MESSAGEINFO
 */
#define SCIP_DECL_MESSAGEOUTPUTFUNC(x) void x (SCIP_MESSAGEHDLR* messagehdlr, FILE* file, const char* msg)


/** error message print method
 *
 *  This method is invoked, if SCIP wants to display an error message to the screen or a file.
 *
 *  @note This function is independent of any message handler.
 *
 *  input:
 *  - data            : data pointer
 *  - file            : file stream to print into
 *  - msg             : string to output into the file (or NULL to flush)
 */
#define SCIP_DECL_ERRORPRINTING(x) void x (void* data, FILE* file, const char* msg)

/** warning message print method of message handler
 *
 *  This method is invoked, if SCIP wants to display a warning message to the screen or a file.
 *
 *  input:
 *  - messagehdlr     : the message handler itself
 *  - file            : file stream to print into
 *  - msg             : string to output into the file (or NULL to flush)
 */
#define SCIP_DECL_MESSAGEWARNING(x) void x (SCIP_MESSAGEHDLR* messagehdlr, FILE* file, const char* msg)

/** dialog message print method of message handler
 *
 *  This method is invoked, if SCIP wants to display a dialog message to the screen or a file.
 *
 *  input:
 *  - messagehdlr     : the message handler itself
 *  - file            : file stream to print into
 *  - msg             : string to output into the file (or NULL to flush)
 */
#define SCIP_DECL_MESSAGEDIALOG(x) void x (SCIP_MESSAGEHDLR* messagehdlr, FILE* file, const char* msg)

/** info message print method of message handler
 *
 *  This method is invoked, if SCIP wants to display an information message to the screen or a file.
 *
 *  input:
 *  - messagehdlr     : the message handler itself
 *  - file            : file stream to print into
 *  - msg             : string to output into the file (or NULL to flush)
 */
#define SCIP_DECL_MESSAGEINFO(x) void x (SCIP_MESSAGEHDLR* messagehdlr, FILE* file, const char* msg)

/** destructor of message handler to free message handler data
 *
 *  This method is invoked, if SCIP wants to free a message handler.
 *
 *  input:
 *  - messagehdlr     : the message handler itself
 */
#define SCIP_DECL_MESSAGEHDLRFREE(x) SCIP_RETCODE x (SCIP_MESSAGEHDLR* messagehdlr)

#ifdef __cplusplus
}
#endif

#endif
