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

/**@file   reader_rlp.c
 * @brief  RLP file reader (LP format with generic variables and row names)
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/reader_lp.h"
#include "scip/reader_rlp.h"

#define READER_NAME             "rlpreader"
#define READER_DESC             "file reader for MIPs in IBM CPLEX's RLP file format"
#define READER_EXTENSION        "rlp"


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyRlp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderRlp(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadRlp)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadLp(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteRlp)
{  /*lint --e{715}*/
   if( genericnames )
   {
      SCIP_CALL( SCIPwriteLp(scip, file, name, transformed, objsense, objscale, objoffset, vars,
            nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );
   }
   else
   {
      SCIPwarningMessage(scip, "RLP format is LP format with generic variable and constraint names\n");

      if( transformed )
      {
         SCIPwarningMessage(scip, "write transformed problem with generic variable and constraint names\n");
         SCIP_CALL( SCIPprintTransProblem(scip, file, "rlp", TRUE) );
      }
      else
      {
         SCIPwarningMessage(scip, "write original problem with generic variable and constraint names\n");
         SCIP_CALL( SCIPprintOrigProblem(scip, file, "rlp", TRUE) );
      }
      *result = SCIP_SUCCESS;
   }
   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the rlp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderRlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyRlp) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadRlp) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteRlp) );

   return SCIP_OKAY;
}
