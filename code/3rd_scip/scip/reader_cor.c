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

/**@file   reader_cor.c
 * @ingroup DEFPLUGINS_READER
 * @brief  COR file reader (MPS format of the core problem for stochastic programs)
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/pub_reader.h"
#include "scip/reader_cor.h"
#include "scip/reader_mps.h"
#include "scip/scip_mem.h"
#include "scip/scip_reader.h"
#include <string.h>

#define READER_NAME             "correader"
#define READER_DESC             "file reader for CORE problem of stochastic programs in the SMPS file format"
#define READER_EXTENSION        "cor"

#define SCIP_DEFAULT_ARRAYSIZE      100

/** COR reading data */
struct SCIP_ReaderData
{
   const char**          varnames;
   const char**          consnames;
   int                   varnamessize;
   int                   consnamessize;
   int                   nvarnames;
   int                   nconsnames;
   SCIP_Bool             read;
};

/** creates the reader data  */
static
SCIP_RETCODE createReaderdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data structure */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);

   readerdata->read = FALSE;
   readerdata->nvarnames = 0;
   readerdata->nconsnames = 0;
   readerdata->varnamessize = SCIP_DEFAULT_ARRAYSIZE;
   readerdata->consnamessize = SCIP_DEFAULT_ARRAYSIZE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->varnames, readerdata->varnamessize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->consnames, readerdata->consnamessize) );

   return SCIP_OKAY;
}

/** creates the reader data  */
static
SCIP_RETCODE freeReaderdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< the reader data structure */
   )
{
   int i;

   assert(scip != NULL);
   assert(readerdata != NULL);

   for( i = readerdata->nvarnames - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &readerdata->varnames[i], strlen(readerdata->varnames[i]) + 1);

   for( i = readerdata->nconsnames - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &readerdata->consnames[i], strlen(readerdata->consnames[i]) + 1);

   SCIPfreeBlockMemoryArray(scip, &readerdata->consnames, readerdata->consnamessize);
   SCIPfreeBlockMemoryArray(scip, &readerdata->varnames, readerdata->varnamessize);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderCor(scip) );

   return SCIP_OKAY;
}


/** destructor of reader to free user data (called when SCIP is exiting) */
/**! [SnippetReaderFreeCor] */
static
SCIP_DECL_READERFREE(readerFreeCor)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( freeReaderdata(scip, readerdata) );

   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}
/**! [SnippetReaderFreeCor] */


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCor)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadCor(scip, filename, result) );

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the cor file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );
   SCIP_CALL( createReaderdata(scip, readerdata) );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCor) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeCor) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCor) );

   return SCIP_OKAY;
}



/** reads problem from file */
SCIP_RETCODE SCIPreadCor(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPreadMps(scip, reader, filename, result, &readerdata->varnames, &readerdata->consnames,
         &readerdata->varnamessize, &readerdata->consnamessize, &readerdata->nvarnames, &readerdata->nconsnames) );

   if( (*result) == SCIP_SUCCESS )
      readerdata->read = TRUE;

   return SCIP_OKAY;
}

/*
 * Interface method for the tim and sto readers
 */


/** returns whether the COR file has been successfully read. This is used by the TIM and STO readers. */
SCIP_Bool SCIPcorHasRead(
   SCIP_READER*          reader              /**< the file reader itself */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->read;
}

/** returns the number of variable names in the COR problem */
int SCIPcorGetNVarNames(
   SCIP_READER*          reader              /**< the file reader itself */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->nvarnames;
}

/** returns the number of constraint names in the COR problem */
int SCIPcorGetNConsNames(
   SCIP_READER*          reader              /**< the file reader itself */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->nconsnames;
}

/** returns the variable name for the given index */
const char* SCIPcorGetVarName(
   SCIP_READER*          reader,             /**< the file reader itself */
   int                   i                   /**< the index of the variable that is requested */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(i >= 0 && i < readerdata->nvarnames);

   return readerdata->varnames[i];
}

/** returns the constraint name for the given index */
const char* SCIPcorGetConsName(
   SCIP_READER*          reader,             /**< the file reader itself */
   int                   i                   /**< the index of the constraint that is requested */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(i >= 0 && i < readerdata->nconsnames);

   return readerdata->consnames[i];
}
