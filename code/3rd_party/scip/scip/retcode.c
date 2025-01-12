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

/**@file   retcode.c
 * @ingroup OTHER_CFILES
 * @brief  methods for return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>

#include "scip/retcode.h"

/** prints error message for return code via message handler */
void SCIPretcodePrint(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to write error message */
   SCIP_RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   switch( retcode )
   {
   case SCIP_OKAY:
      SCIPmessageFPrintInfo(messagehdlr, file, "normal termination");
      break;
   case SCIP_ERROR:
      SCIPmessageFPrintInfo(messagehdlr, file, "unspecified error");
      break;
   case SCIP_NOMEMORY:
      SCIPmessageFPrintInfo(messagehdlr, file, "insufficient memory error");
      break;
   case SCIP_READERROR:
      SCIPmessageFPrintInfo(messagehdlr, file, "read error");
      break;
   case SCIP_WRITEERROR:
      SCIPmessageFPrintInfo(messagehdlr, file, "write error");
      break;
   case SCIP_NOFILE:
      SCIPmessageFPrintInfo(messagehdlr, file, "file not found error");
      break;
   case SCIP_FILECREATEERROR:
      SCIPmessageFPrintInfo(messagehdlr, file, "cannot create file");
      break;
   case SCIP_LPERROR:
      SCIPmessageFPrintInfo(messagehdlr, file, "error in LP solver");
      break;
   case SCIP_NOPROBLEM:
      SCIPmessageFPrintInfo(messagehdlr, file, "no problem exists");
      break;
   case SCIP_INVALIDCALL:
      SCIPmessageFPrintInfo(messagehdlr, file, "method cannot be called at this time in solution process");
      break;
   case SCIP_INVALIDDATA:
      SCIPmessageFPrintInfo(messagehdlr, file, "method cannot be called with this type of data");
      break;
   case SCIP_INVALIDRESULT:
      SCIPmessageFPrintInfo(messagehdlr, file, "method returned an invalid result code");
      break;
   case SCIP_PLUGINNOTFOUND:
      SCIPmessageFPrintInfo(messagehdlr, file, "a required plugin was not found");
      break;
   case SCIP_PARAMETERUNKNOWN:
      SCIPmessageFPrintInfo(messagehdlr, file, "the parameter with the given name was not found");
      break;
   case SCIP_PARAMETERWRONGTYPE:
      SCIPmessageFPrintInfo(messagehdlr, file, "the parameter is not of the expected type");
      break;
   case SCIP_PARAMETERWRONGVAL:
      SCIPmessageFPrintInfo(messagehdlr, file, "the value is invalid for the given parameter");
      break;
   case SCIP_KEYALREADYEXISTING:
      SCIPmessageFPrintInfo(messagehdlr, file, "the given key is already existing in table");
      break;
   case SCIP_MAXDEPTHLEVEL:
      SCIPmessageFPrintInfo(messagehdlr, file, "maximal branching depth level exceeded");
      break;
   case SCIP_BRANCHERROR:
      SCIPmessageFPrintInfo(messagehdlr, file, "branching could not be performed (e.g. too large values in variable domain)");
      break;
   case SCIP_NOTIMPLEMENTED:
      SCIPmessageFPrintInfo(messagehdlr, file, "function not implemented");
      break;
   default:
      SCIPmessageFPrintInfo(messagehdlr, file, "unknown error code");
      break;
   }
}

/** prints error message for return code via error message */
void SCIPretcodePrintError(
   SCIP_RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   switch( retcode )
   {
   case SCIP_OKAY:
      SCIPmessagePrintError("normal termination");
      break;
   case SCIP_ERROR:
      SCIPmessagePrintError("unspecified error");
      break;
   case SCIP_NOMEMORY:
      SCIPmessagePrintError("insufficient memory error");
      break;
   case SCIP_READERROR:
      SCIPmessagePrintError("read error");
      break;
   case SCIP_WRITEERROR:
      SCIPmessagePrintError("write error");
      break;
   case SCIP_NOFILE:
      SCIPmessagePrintError("file not found error");
      break;
   case SCIP_FILECREATEERROR:
      SCIPmessagePrintError("cannot create file");
      break;
   case SCIP_LPERROR:
      SCIPmessagePrintError("error in LP solver");
      break;
   case SCIP_NOPROBLEM:
      SCIPmessagePrintError("no problem exists");
      break;
   case SCIP_INVALIDCALL:
      SCIPmessagePrintError("method cannot be called at this time in solution process");
      break;
   case SCIP_INVALIDDATA:
      SCIPmessagePrintError("method cannot be called with this type of data");
      break;
   case SCIP_INVALIDRESULT:
      SCIPmessagePrintError("method returned an invalid result code");
      break;
   case SCIP_PLUGINNOTFOUND:
      SCIPmessagePrintError("a required plugin was not found");
      break;
   case SCIP_PARAMETERUNKNOWN:
      SCIPmessagePrintError("the parameter with the given name was not found");
      break;
   case SCIP_PARAMETERWRONGTYPE:
      SCIPmessagePrintError("the parameter is not of the expected type");
      break;
   case SCIP_PARAMETERWRONGVAL:
      SCIPmessagePrintError("the value is invalid for the given parameter");
      break;
   case SCIP_KEYALREADYEXISTING:
      SCIPmessagePrintError("the given key is already existing in table");
      break;
   case SCIP_MAXDEPTHLEVEL:
      SCIPmessagePrintError("maximal branching depth level exceeded");
      break;
   case SCIP_BRANCHERROR:
      SCIPmessagePrintError("branching could not be performed (e.g. too large values in variable domain)");
      break;
   case SCIP_NOTIMPLEMENTED:
      SCIPmessagePrintError("function not implemented");
      break;
   default:
      SCIPmessagePrintError("unknown error code");
      break;
   }
}
