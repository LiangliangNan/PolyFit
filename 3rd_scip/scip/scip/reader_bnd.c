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

/**@file   reader_bnd.c
 * @brief  file reader for variable bounds
 * @author Ambros Gleixner
 * @author Ingmar Vierhaus
 * @author Benjamin Mueller
 *
 * This reader allows to read a file containing new bounds for variables of the current problem.  Each line of the file
 * should have format
 *
 *    \<variable name\> \<lower bound\> \<upper bound\>
 *
 * where infinite bounds can be written as inf, +inf or -inf.  Note that only a subset of the variables may appear in
 * the file.  Lines with unknown variable names are ignored.
 * The writing functionality can be used in problem and transformed stages. Note that in transformed stage,
 * the leading "t_" in the name of a transformed variable will not appear in the output. This way, bounds written in transformed stage
 * can be read again in problem stage.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#include "scip/reader_bnd.h"
#include "scip/pub_misc.h"

#define READER_NAME             "bndreader"
#define READER_DESC             "file reader for variable bounds"
#define READER_EXTENSION        "bnd"

#define DEFAULT_IMPROVEONLY     FALSE        /**< only use improving bounds */


/** BND reader data */
struct SCIP_ReaderData
{
   SCIP_Bool             improveonly;        /**< only use improving bounds */
};


/*
 * Local methods of reader
 */

/** reads a given bound file, problem has to be in problem stage */
static
SCIP_RETCODE readBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname,              /**< name of the input file */
   SCIP_READERDATA*      readerdata          /**< pointer to the data of the reader */
   )
{
   SCIP_FILE* file;
   SCIP_Bool error;
   SCIP_Bool unknownvariablemessage;
   SCIP_Bool usevartable;
   int lineno;

   assert(scip != NULL);
   assert(fname != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/usevartable", &usevartable) );

   if( !usevartable )
   {
      SCIPerrorMessage("Cannot read bounds file if vartable is disabled. Make sure parameter 'misc/usevartable' is set to TRUE.\n");
      return SCIP_READERROR;
   }

   /* open input file */
   file = SCIPfopen(fname, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", fname);
      SCIPprintSysError(fname);
      return SCIP_NOFILE;
   }

   /* read the file */
   error = FALSE;
   unknownvariablemessage = FALSE;
   lineno = 0;
   while( !SCIPfeof(file) && !error )
   {
      char buffer[SCIP_MAXSTRLEN];
      char varname[SCIP_MAXSTRLEN];
      char lbstring[SCIP_MAXSTRLEN];
      char ubstring[SCIP_MAXSTRLEN];
      char format[SCIP_MAXSTRLEN];
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      int nread;
      char* endptr;

      /* get next line */
      if( SCIPfgets(buffer, (int) sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* parse the line */
      (void) SCIPsnprintf(format, SCIP_MAXSTRLEN, "%%%ds %%%ds %%%ds\n", SCIP_MAXSTRLEN, SCIP_MAXSTRLEN, SCIP_MAXSTRLEN);
      (void) sscanf(buffer, format, varname, lbstring, ubstring);

      SCIP_CALL( SCIPparseVarName(scip, buffer, &var, &endptr) );

      /* cppcheck-suppress invalidscanf */
      nread = sscanf(endptr, "%s %s\n", lbstring, ubstring);
      if( nread < 1 )
      {
         SCIPerrorMessage("invalid input line %d in bounds file <%s>: <%s>\n", lineno, fname, buffer);
         error = TRUE;
         break;
      }

      if( var == NULL )
      {
         if( !unknownvariablemessage )
         {
            SCIPwarningMessage(scip, "unable to parse variable name in line %d of bounds file <%s>:\n", lineno, fname);
            SCIPwarningMessage(scip, "line is: %s", buffer);
            SCIPwarningMessage(scip, "  (further unknown variables are ignored)\n");
            unknownvariablemessage = TRUE;
         }
         continue;
      }

      /* cast the lower bound value */
      if( strncasecmp(lbstring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(lbstring, "+inf", 4) == 0 || strncasecmp(lbstring, "inf", 3) == 0 )
         lb = SCIPinfinity(scip);
      else if( strncasecmp(lbstring, "-inf", 4) == 0 )
         lb = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(lbstring, "%lf", &lb);
         if( nread != 1 )
         {
            SCIPerrorMessage("invalid lower bound value <%s> for variable <%s> in line %d of bounds file <%s>\n",
               lbstring, varname, lineno, fname);
            error = TRUE;
            break;
         }
      }

      /* cast the upper bound value */
      if( strncasecmp(ubstring, "inv", 3) == 0 )
         continue;
      else if( strncasecmp(ubstring, "+inf", 4) == 0 || strncasecmp(ubstring, "inf", 3) == 0 )
         ub = SCIPinfinity(scip);
      else if( strncasecmp(ubstring, "-inf", 4) == 0 )
         ub = -SCIPinfinity(scip);
      else
      {
         nread = sscanf(ubstring, "%lf", &ub);
         if( nread != 1 )
         {
            SCIPerrorMessage("invalid lower bound value <%s> for variable <%s> in line %d of bounds file <%s>\n",
               ubstring, varname, lineno, fname);
            error = TRUE;
            break;
         }
      }

      if( readerdata->improveonly )
      {
         if( SCIPisLT(scip, lb, SCIPvarGetLbGlobal(var)) )
         {
            SCIPwarningMessage(scip, "not applying lower bound value %s for variable <%s> in line %d of bounds file %s,"
               " because it does not improve existing bound of %f\n",
               lbstring, SCIPvarGetName(var), lineno, fname, SCIPvarGetLbGlobal(var));
         }
         if( SCIPisGT(scip, ub, SCIPvarGetUbGlobal(var)) )
         {
            SCIPwarningMessage(scip, "not applying upper bound value %s for variable <%s> in line %d of bounds file %s, "
               "because it does not improve existing bound of %f\n",
               ubstring, SCIPvarGetName(var), lineno, fname, SCIPvarGetUbGlobal(var));
         }

         /* collect best variable bounds */
         lb = MAX(lb, SCIPvarGetLbGlobal(var)); /*lint !e666*/
         ub = MIN(ub, SCIPvarGetUbGlobal(var)); /*lint !e666*/
      }

      /* note that we don't need to check if lb > ub in SCIPchgVar{Lb,Ub} */
      SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
   }

   /* close input file */
   SCIPfclose(file);

   /* return error if necessary */
   if ( error )
      return SCIP_READERROR;

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyBnd)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderBnd(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader
 *
 *  In order to determine the type of the file, we have to open it. Thus, it has to be opened
 *  twice. This might be removed, but is likely to not hurt the performance too much.
 */
static
SCIP_DECL_READERREAD(readerReadBnd)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("reading of bounds file is only possible after a problem was created\n");
      return SCIP_READERROR;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("reading of bounds file is only possible during problem creation stage\n");
      return SCIP_READERROR;
   }

   /* read bounds file */
   SCIP_CALL( readBounds(scip, filename, SCIPreaderGetData(reader)) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/** outputs given bounds into a file stream */
static
void printBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   /* print lower bound */
   if( SCIPisInfinity(scip, lb) )
      SCIPmessageFPrintInfo(messagehdlr, file, "+inf ");
   else if( SCIPisInfinity(scip, -lb) )
      SCIPmessageFPrintInfo(messagehdlr, file, "-inf ");
   else
      SCIPmessageFPrintInfo(messagehdlr, file, "%.15" SCIP_REAL_FORMAT " ", lb);

   /* print upper bound */
   if( SCIPisInfinity(scip, ub) )
      SCIPmessageFPrintInfo(messagehdlr, file, "+inf");
   else if( SCIPisInfinity(scip, -ub) )
      SCIPmessageFPrintInfo(messagehdlr, file, "-inf");
   else
      SCIPmessageFPrintInfo(messagehdlr, file, "%.15" SCIP_REAL_FORMAT, ub);
}

/** writes problem to file */
static
SCIP_RETCODE SCIPwriteBnd(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{
   SCIP_MESSAGEHDLR* messagehdlr;
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(result != NULL);

   messagehdlr = SCIPgetMessagehdlr(scip);
   *result = SCIP_SUCCESS;

   if( nvars == 0 )
   {
      SCIPwarningMessage(scip, "Problem has no variables, no bounds written.\n");
      return SCIP_OKAY;
   }

   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      const char* varname;

      var = vars[i];
      assert( var != NULL );
      varname = SCIPvarGetName(var);

      /* strip 't_' from varname */
      if( SCIPvarIsTransformedOrigvar(var) && strncmp(SCIPvarGetName(var), "t_", 2) == 0)
      {
         varname = varname + 2;
      }

      SCIPinfoMessage(scip, file, "<%s> ", varname);

      /* print global bounds for transformed variables, original bounds for original variables */
      if( !SCIPvarIsTransformed(var) )
      {
         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);
      }
      else
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
      }

      /* print bounds into the file */
      printBounds(scip, messagehdlr, file, lb, ub);
      SCIPmessageFPrintInfo(messagehdlr, file, "\n");
   }

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteBnd)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   SCIP_CALL( SCIPwriteBnd(scip, file, vars, nvars, result) );

   return SCIP_OKAY;
}

/** destructor of reader to free reader data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeBnd)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/*
 * bnd file reader specific interface methods
 */

/** includes the bnd file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderBnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyBnd) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadBnd) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteBnd) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeBnd) );

   /* add bnd reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/bndreader/improveonly", "only use improving bounds",
         &readerdata->improveonly, FALSE, DEFAULT_IMPROVEONLY, NULL, NULL) );

   return SCIP_OKAY;
}
