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

/**@file   reader_wbo.c
 * @ingroup DEFPLUGINS_READER
 * @brief  WBO file reader (OPB format with weighted constraints)
 * @author Michael Winkler
 */

/* For file format description see opb-reader. */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/pub_reader.h"
#include "scip/reader_opb.h"
#include "scip/reader_wbo.h"
#include "scip/scip_reader.h"
#include <string.h>

#define READER_NAME             "wboreader"
#define READER_DESC             "file reader for pseudoboolean wbo file format"
#define READER_EXTENSION        "wbo"

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyWbo)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderWbo(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadWbo)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadOpb(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteWbo)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPwriteOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
         nvars, nbinvars, nintvars, nimplvars, ncontvars, fixedvars, nfixedvars, conss, nconss, genericnames, result) );

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the wbo file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderWbo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyWbo) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadWbo) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteWbo) );

   return SCIP_OKAY;
}
