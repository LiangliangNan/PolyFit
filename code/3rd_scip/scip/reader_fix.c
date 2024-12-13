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

/**@file   reader_fix.c
 * @ingroup DEFPLUGINS_READER
 * @brief  file reader for variable fixings
 * @author Tobias Achterberg
 *
 * This reader allows to read a file containing fixation values for variables of the current problem. Each line of the
 * file should have format
 *
 *    \<variable name\> \<value to fix\>
 *
 * Note that only a subset of the variables may need to appear in the file. Lines with unknown variable names are
 * ignored. The writing functionality is currently not supported.
 *
 * @note The format is equal to the (not xml) solution format of SCIP.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/pub_var.h"
#include "scip/reader_fix.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include "scip/scip_solve.h"
#include "scip/scip_var.h"
#include <string.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <strings.h> /*lint --e{766}*/ /* needed for strncasecmp() */
#endif



#define READER_NAME             "fixreader"
#define READER_DESC             "file reader for variable fixings"
#define READER_EXTENSION        "fix"


/*
 * local methods
 */

/** reads the given solution file */
static
SCIP_RETCODE readSol(
   SCIP*                 scip,               /**< SCIP data structure */   
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_RETCODE retcode;
   SCIP_FILE* file;
   SCIP_Bool error;
   SCIP_Bool unknownvariablemessage;
   int lineno;
   int nfixed;

   assert(scip != NULL);
   assert(filename != NULL);

   /* open input file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }   

   /* read the file */
   error = FALSE;
   unknownvariablemessage = FALSE;
   lineno = 0;
   nfixed = 0;
   while( !SCIPfeof(file) && !error )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char valuestring[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      char format[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real value;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      int nread;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* the lines "solution status: ..." and "objective value: ..." may preceed the solution information */
      if( strncasecmp(buffer, "solution status:", 16) == 0 || strncasecmp(buffer, "objective value:", 16) == 0 )
         continue;

      /* parse the line */
      (void) SCIPsnprintf(format, SCIP_MAXSTRLEN, "%%%ds %%%ds %%%ds\n", SCIP_MAXSTRLEN, SCIP_MAXSTRLEN, SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, varname, valuestring, objstring);
      if( nread < 2 )
      {
         SCIPerrorMessage("invalid input line %d in solution file <%s>: <%s>\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      /* find the variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPwarningMessage(scip, "unknown variable <%s> in line %d of solution file <%s>\n", varname, lineno, filename);
            SCIPwarningMessage(scip, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the value */
      if( strncasecmp(valuestring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(valuestring, "+inf", 4) == 0 || strncasecmp(valuestring, "inf", 3) == 0 )
         value = SCIPinfinity(scip);
      else if( strncasecmp(valuestring, "-inf", 4) == 0 )
         value = -SCIPinfinity(scip);
      else
      {
         /* coverity[secure_coding] */
         nread = sscanf(valuestring, "%lf", &value);
         if( nread != 1 )
         {
            SCIPerrorMessage("invalid solution value <%s> for variable <%s> in line %d of solution file <%s>\n",
               valuestring, varname, lineno, filename);
            error = TRUE;
            break;
         }
      }

      /* fix the variable */
      retcode = SCIPfixVar(scip, var, value, &infeasible, &fixed);
      if( retcode != SCIP_OKAY )
      {
         SCIPerrorMessage("Error fixing variable <%s> to value %.15g in line %d of bounds file <%s>\n",
            varname, value, lineno, filename);
         error = TRUE;
         break;
      }
      if( infeasible )
      {
         SCIPerrorMessage("infeasible solution value of <%s>[%.15g,%.15g] to %.15g in line %d of solution file <%s>\n",
            varname, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), value, lineno, filename);
         error = TRUE;
         break;
      }
      if( fixed )
         nfixed++;
   }

   /* close input file */
   SCIPfclose(file);

   /* display result */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "fixed %d variables from solution file <%s>\n", nfixed, filename);

   if( error )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyFix)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderFix(scip) );

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadFix)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("reading of fixing file is only possible after a problem was created\n");
      return SCIP_READERROR;
   }

   /* free transformed problem, s.t. fixings are applied to the original problem */
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* read (partial) solution from fixing file */
   SCIP_CALL( readSol(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * fix file reader specific interface methods
 */

/** includes the fix file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderFix(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyFix) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadFix) );

   return SCIP_OKAY;
}
