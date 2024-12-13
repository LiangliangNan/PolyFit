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

/**@file   reader_smps.c
 * @ingroup DEFPLUGINS_READER
 * @brief  SMPS file reader - smps files list the cor, tim and sto files for a single instance
 * @author Stephen J. Maher
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/reader_cor.h"
#include "scip/reader_smps.h"
#include "scip/reader_sto.h"
#include "scip/reader_tim.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include <string.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <strings.h> /*lint --e{766}*/ /* needed for strncasecmp() */
#endif

/*
 * The SMPS reader coordinates the reading of the cor, tim and sto files. The public reading methods from the cor, tim
 * and sto readers are called from the SMPS reader. So, the header files for the cor, tim and sto readers are required.
 */

#define READER_NAME             "smpsreader"
#define READER_DESC             "file reader for core problem of stochastic programs in the SMPS file format"
#define READER_EXTENSION        "smps"

#define SMPS_MAX_LINELEN  1024
#define BLANK              ' '
#define LINEWIDTH           80

#define COR_FILEEXTENSION        "cor"
#define TIM_FILEEXTENSION        "tim"
#define STO_FILEEXTENSION        "sto"

/** enum for the file types that are read by the SMPS reader */
enum SCIP_SmpsFileType
{
   SCIP_SMPSFILETYPE_COR = 0,
   SCIP_SMPSFILETYPE_TIM = 1,
   SCIP_SMPSFILETYPE_STO = 2
};
typedef enum SCIP_SmpsFileType SCIP_SMPSFILETYPE;


/** smps input structure */
struct SmpsInput
{
   SCIP_FILE*            fp;
   int                   lineno;
   SCIP_Bool             haserror;
   char                  buf[SMPS_MAX_LINELEN];
   const char*           f0;
   const char*           f1;
};
typedef struct SmpsInput SMPSINPUT;


/** creates the smps input structure */
static
SCIP_RETCODE smpsinputCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SMPSINPUT**           smpsi,               /**< smps input structure */
   SCIP_FILE*            fp                  /**< file object for the input file */
   )
{
   assert(smpsi != NULL);
   assert(fp != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, smpsi) );

   (*smpsi)->fp          = fp;
   (*smpsi)->lineno      = 0;
   (*smpsi)->haserror    = FALSE;
   (*smpsi)->buf     [0] = '\0';
   (*smpsi)->f0          = NULL;
   (*smpsi)->f1          = NULL;

   return SCIP_OKAY;
}

/** free the smps input structure */
static
void smpsinputFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SMPSINPUT**           smpsi               /**< smps input structure */
   )
{
   SCIPfreeBlockMemory(scip, smpsi);
}

/** return the current value of field 0 */
static
const char* smpsinputField0(
   const SMPSINPUT*      smpsi               /**< smps input structure */
   )
{
   assert(smpsi != NULL);

   return smpsi->f0;
}

/** fill the line from \p pos up to column LINEWIDTH with blanks. */
static
void clearFrom(
   char*                 buf,                /**< buffer to clear */
   unsigned int          pos                 /**< position to start the clearing process */
   )
{
   unsigned int i;

   for(i = pos; i < LINEWIDTH; i++)
      buf[i] = BLANK;
   buf[LINEWIDTH] = '\0';
}

/** read a smps format data line and parse the fields. */
static
SCIP_Bool smpsinputReadLine(
   SMPSINPUT*            smpsi               /**< smps input structure */
   )
{
   unsigned int len;
   unsigned int i;
   SCIP_Bool is_marker;
   SCIP_Bool is_empty;
   char* nexttok;

   do
   {
      smpsi->f0 = smpsi->f1 = 0;
      is_marker = FALSE;

      /* Read until we have not a comment line. */
      do
      {
         smpsi->buf[SMPS_MAX_LINELEN-1] = '\0';
         if( NULL == SCIPfgets(smpsi->buf, (int) sizeof(smpsi->buf), smpsi->fp) )
            return FALSE;
         smpsi->lineno++;
      }
      while( *smpsi->buf == '*' );

      /* Normalize line */
      len = (unsigned int) strlen(smpsi->buf);

      /* replace tabs and new lines by blanks */
      for( i = 0; i < len; i++ )
      {
         if( (smpsi->buf[i] == '\t') || (smpsi->buf[i] == '\n') || (smpsi->buf[i] == '\r') )
            smpsi->buf[i] = BLANK;
      }

      if( len < LINEWIDTH )
         clearFrom(smpsi->buf, len);

      SCIPdebugMessage("line %d: <%s>\n", smpsi->lineno, smpsi->buf);

      assert(strlen(smpsi->buf) >= LINEWIDTH);

      /* Look for new section */
      if( *smpsi->buf != BLANK )
      {
         smpsi->f0 = SCIPstrtok(&smpsi->buf[0], " ", &nexttok);

         assert(smpsi->f0 != 0);

         smpsi->f1 = SCIPstrtok(NULL, " ", &nexttok);

         return TRUE;
      }

      /* check for empty lines */
      is_empty = (smpsi->f0 == NULL && smpsi->f1 == NULL);
   }
   while( is_marker || is_empty );

   return TRUE;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySmps)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderSmps(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSmps)
{  /*lint --e{715}*/
   SCIP_FILE* fp;
   SMPSINPUT* smpsi;
   SCIP_RETCODE retcode = SCIP_OKAY;

   char corfilename[SCIP_MAXSTRLEN];
   char timfilename[SCIP_MAXSTRLEN];
   char stofilename[SCIP_MAXSTRLEN];
   char* tmpfilename;
   char* probname;
   char* fileextension;
   char* fromlastslash;
   char parent[SCIP_MAXSTRLEN];
   size_t parentlen;

   SCIP_Bool hascorfile;
   SCIP_Bool hastimfile;
   SCIP_Bool hasstofile;

   int i;

   assert(scip != NULL);
   assert(filename != NULL);

   /* copy filename */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpfilename, filename, (int)strlen(filename)+1) );

   /* getting the problem name from the SMPS file name */
   SCIPsplitFilename(tmpfilename, NULL, &probname, NULL, NULL);

   fromlastslash = (char*) strrchr(filename, '/');

   if( fromlastslash == NULL )
      parentlen = 0;
   else
      parentlen = strlen(filename) - (strlen(fromlastslash) - 1);

   (void)SCIPstrncpy(parent, filename, (int)parentlen + 1);

   fp = SCIPfopen(filename, "r");
   if( fp == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);

      return SCIP_NOFILE;
   }

   SCIP_CALL( smpsinputCreate(scip, &smpsi, fp) );

   hascorfile = FALSE;
   hastimfile = FALSE;
   hasstofile = FALSE;
   while( smpsinputReadLine(smpsi) )
   {
      char* tmpinput;

      /* copy the input */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpinput, smpsinputField0(smpsi),
            (int)strlen(smpsinputField0(smpsi))+1) ); /*lint !e666*/

      /* get extension from filename */
      SCIPsplitFilename(tmpinput, NULL, NULL, &fileextension, NULL);

      if( strcasecmp(fileextension, COR_FILEEXTENSION) == 0 )
      {
         (void) SCIPsnprintf(corfilename, SCIP_MAXSTRLEN, "%s%s", parent, smpsinputField0(smpsi));
         hascorfile = TRUE;
      }
      else if( strcasecmp(fileextension, TIM_FILEEXTENSION) == 0 )
      {
         (void) SCIPsnprintf(timfilename, SCIP_MAXSTRLEN, "%s%s", parent, smpsinputField0(smpsi));
         hastimfile = TRUE;
      }
      else if( strcasecmp(fileextension, STO_FILEEXTENSION) == 0 )
      {
         (void) SCIPsnprintf(stofilename, SCIP_MAXSTRLEN, "%s%s", parent, smpsinputField0(smpsi));
         hasstofile = TRUE;
      }

      SCIPfreeBufferArray(scip, &tmpinput);
   }

   /* printing errors if the correct files have not been provided */
   if( !hascorfile )
   {
      SCIPerrorMessage("The core file has not been listed in <%s>\n", filename);
   }

   if( !hastimfile )
   {
      SCIPerrorMessage("The tim file has not been listed in <%s>\n", filename);
   }

   if( !hasstofile )
   {
      SCIPerrorMessage("The sto file has not been listed in <%s>\n", filename);
   }

   /* if one of the necessary file has not been provided, then an error will be returned */
   if( !hascorfile || !hastimfile || !hasstofile )
   {
      retcode = SCIP_READERROR;
      goto TERMINATE;
   }

   for( i = 0; i < 3; i++ )
   {
      int nvars;
      int nbinvars;
      int nintvars;
      int nimplintvars;
      int ncontvars;
      SCIP_SMPSFILETYPE type;

      type = (SCIP_SMPSFILETYPE) i;
      switch( type )
      {
         case SCIP_SMPSFILETYPE_COR:
            SCIPinfoMessage(scip, NULL, "reading core file <%s> for problem %s\n", corfilename, probname);
            SCIPinfoMessage(scip, NULL, "============\n");

            /* reading the CORE file */
            SCIP_CALL_TERMINATE( retcode, SCIPreadCor(scip, corfilename, result), TERMINATE );

            /* getting the variable information */
            SCIP_CALL( SCIPgetOrigVarsData(scip, NULL, &nvars, &nbinvars, &nintvars, &nimplintvars, &ncontvars) );
            SCIPinfoMessage(scip, NULL,
               "core problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
               nvars, nbinvars, nintvars, nimplintvars, ncontvars, SCIPgetNOrigConss(scip));
            break;
         case SCIP_SMPSFILETYPE_TIM:
            SCIPinfoMessage(scip, NULL, "reading the time file <%s> for problem %s\n", timfilename, probname);
            SCIPinfoMessage(scip, NULL, "============\n");

            /* reading the TIME file */
            SCIP_CALL_TERMINATE( retcode, SCIPreadTim(scip, timfilename, result), TERMINATE );

            SCIPinfoMessage(scip, NULL, "problem %s has %d stages\n", probname, SCIPtimGetNStages(scip));
            break;
         case SCIP_SMPSFILETYPE_STO:
#ifdef BENDERSBRANCH
            SCIP_Bool usebenders;
#endif

            SCIPinfoMessage(scip, NULL, "read problem <%s>\n", stofilename);
            SCIPinfoMessage(scip, NULL, "============\n");

            /* reading the STO file */
            SCIP_CALL_TERMINATE( retcode, SCIPreadSto(scip, stofilename, result), TERMINATE );

            SCIPinfoMessage(scip, NULL, "problem %s has extended with a total of %d scenarios\n", probname,
               SCIPstoGetNScenarios(scip));

            /* getting the variable information */
            SCIP_CALL( SCIPgetOrigVarsData(scip, NULL, &nvars, &nbinvars, &nintvars, &nimplintvars, &ncontvars) );

            /* if Benders' decomposition is used, the variable will be distributed to a number of subproblems */
#ifdef BENDERSBRANCH
            SCIP_CALL( SCIPgetBoolParam(scip, "reading/sto/usebenders", &usebenders) );
            if( usebenders )
            {
               SCIPinfoMessage(scip, NULL, "Benders' decomposition master problem ");
            }
            else
#endif
            {
               SCIPinfoMessage(scip, NULL, "deterministic equivalent problem ");
            }

            SCIPinfoMessage(scip, NULL,
               "has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
               nvars, nbinvars, nintvars, nimplintvars, ncontvars, SCIPgetNOrigConss(scip));
            break;
         /* coverity[dead_error_begin] */
         default:
            SCIPerrorMessage("This should not happen. Aborting.\n");
            SCIPABORT();
            retcode = SCIP_READERROR;
            goto TERMINATE;
      }

      SCIPinfoMessage(scip, NULL, "\n\n");
   }

   SCIPfclose(fp);

 /* cppcheck-suppress unusedLabel */
TERMINATE:
   smpsinputFree(scip, &smpsi);

   /* freeing buffer array */
   SCIPfreeBufferArray(scip, &tmpfilename);

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   if( retcode == SCIP_NOFILE || retcode == SCIP_READERROR )
      return retcode;

   SCIP_CALL( retcode );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the smps file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSmps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySmps) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSmps) );

   return SCIP_OKAY;
}
