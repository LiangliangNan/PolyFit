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

/**@file   reader_sol.c
 * @ingroup DEFPLUGINS_READER
 * @brief  file reader for primal solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/pub_sol.h"
#include "scip/reader_sol.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_param.h"
#include "scip/scip_reader.h"
#include "scip/scip_sol.h"
#include <string.h>

#define READER_NAME             "solreader"
#define READER_DESC             "file reader for primal solutions"
#define READER_EXTENSION        "sol"


/*
 * Local methods of reader
 */

/** reads a given SCIP solution file, problem has to be transformed in advance */
static
SCIP_RETCODE readSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname,              /**< name of the input file */
   SCIP_Bool             xml                 /**< true, iff the given file is XML */
   )
{
   SCIP_SOL* sol;
   SCIP_Bool error;
   SCIP_Bool partial;
   SCIP_Bool stored;
   SCIP_Bool usevartable;

   assert(scip != NULL);
   assert(fname != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/usevartable", &usevartable) );

   if( !usevartable )
   {
      SCIPerrorMessage("Cannot read solution file if vartable is disabled. Make sure parameter 'misc/usevartable' is set to TRUE.\n");
      return SCIP_READERROR;
   }

   /* create zero solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   SCIP_CALL( SCIPreadSolFile(scip, fname, sol, xml, &partial, &error) );

   if( !error )
   {
      /* add and free the solution */
      if( SCIPisTransformed(scip) )
      {
         SCIP_Bool completely;

         assert(!partial);
         assert(!SCIPsolIsPartial(sol));

         /* use display/allviols to decide whether to print all violations or just the first one */
         SCIP_CALL( SCIPgetBoolParam(scip, "display/allviols", &completely) );

         SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, completely, TRUE, TRUE, TRUE, &stored) );

         /* display result */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "primal solution from solution file <%s> was %s\n",
            fname, stored ? "accepted" : "rejected - solution is infeasible or objective too poor");
      }
      else
      {
         /* add primal solution to solution candidate storage, frees the solution afterwards */
         SCIP_CALL( SCIPaddSolFree(scip, &sol, &stored) );

         /* display result */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "%sprimal solution from solution file <%s> was %s\n",
            partial ? "partial " : "", fname, stored ? "accepted as candidate, will be checked when solving starts" : "rejected - solution objective too poor");
      }

      return SCIP_OKAY;
   }
   else
   {
      /* free solution */
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      return SCIP_READERROR;
   }
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderSol(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader 
 *
 *  In order to determine the type of the file, we have to open it. Thus, it has to be opened
 *  twice. This might be removed, but is likely to not hurt the performance too much.
 */
static
SCIP_DECL_READERREAD(readerReadSol)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   char buffer[SCIP_MAXSTRLEN];

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("reading of solution file is only possible after a problem was created\n");
      return SCIP_READERROR;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
         "primal solution from solution file <%s> was ignored - problem is already solved to optimality\n",
         filename);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* open input file in order to determine type */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* get next line */
   if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
   {
      SCIPerrorMessage("cannot parse file.\n");
      SCIPfclose(file);
      return SCIP_READERROR;
   }
   /* close file */
   SCIPfclose(file);

   /* decide whether it is xml */
   if( SCIPstrAtStart(buffer, "<?xml", (size_t) 5) )
   {
      /* read XML solution and add it to the solution pool */
      SCIP_CALL( readSol(scip, filename, TRUE) );
   }
   else
   {
      /* read the solution and add it to the solution pool */
      SCIP_CALL( readSol(scip, filename, FALSE) );
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * sol file reader specific interface methods
 */

/** includes the sol file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySol) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSol) );

   return SCIP_OKAY;
}

