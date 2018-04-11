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

/**@file   reader_ppm.c
 * @brief  file writer for portable pixmap file format (PPM), open with common graphic viewer programs (e.g. xview)
 * @author Michael Winkler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "scip/reader_ppm.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/pub_misc.h"

#define READER_NAME             "ppmreader"
#define READER_DESC             "file writer for portable pixmap file format (PPM), open with common graphic viewer programs (e.g. xview)"
#define READER_EXTENSION        "ppm"

/*
 * Data structures
 */
#define PPM_MAX_LINELEN               71      /**< the maximum length of any line is 70 + '\\0' = 71*/
#define DEFAULT_PPM_RGB_LIMIT        160
#define DEFAULT_PPM_COEF_LIMIT         3
#define DEFAULT_PPM_RGB_RELATIVE    TRUE
#define DEFAULT_PPM_RGB_ASCII       TRUE

/** PPM reading data */
struct SCIP_ReaderData
{
   SCIP_Bool             rgb_relative;
   SCIP_Bool             rgb_ascii;
   int                   rgb_limit;
   int                   coef_limit;
};

/*
 * Local methods (for writing)
 */

/** initializes the reader data */
static
void initReaderdata(
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   assert(readerdata != NULL);

   readerdata->rgb_relative = DEFAULT_PPM_RGB_RELATIVE;
   readerdata->rgb_ascii = DEFAULT_PPM_RGB_ASCII;
   readerdata->rgb_limit = DEFAULT_PPM_RGB_LIMIT;
   readerdata->coef_limit = DEFAULT_PPM_COEF_LIMIT;
}


/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n inrc/scip/reader_ppm.c linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( scalars != NULL );
   assert( nvars != NULL );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert( requiredsize <= *nvars );
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
      }
   }
   return SCIP_OKAY;
}

/** clears the given line buffer */
static
void clearLine(
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of characters in line */
   )
{
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   (*linecnt) = 0;
   linebuffer[0] = '\0';
}

/** ends the given line with '\\0' and prints it to the given file stream */
static
void endLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of characters in line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );

   if( (*linecnt) > 0 )
   {
      linebuffer[(*linecnt)] = '\0';

      if(readerdata->rgb_ascii)
         SCIPinfoMessage(scip, file, "%s", linebuffer);
      else
         SCIPinfoMessage(scip, file, "%s\n", linebuffer);
      clearLine(linebuffer, linecnt);
   }
}

/** appends extension to line and prints it to the give file stream if the line exceeded PPM_PRINTLEN */
static
void appendLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of characters in line */
   const char*           extension           /**< string to extent the line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( extension != NULL );

   if( *linecnt + strlen(extension) > PPM_MAX_LINELEN - 1 )
      endLine(scip, file, readerdata, linebuffer, linecnt);

   /* append extension to linebuffer */
   strncat(linebuffer, extension, PPM_MAX_LINELEN - (unsigned int)(*linecnt) - 1);
   (*linecnt) += (int) strlen(extension);
}


/** calculates the color value for a given coefficient */
static
void calcColorValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_Real             coef,               /**< coefficient to scale */
   int*                  red,                /**< red part */
   int*                  green,              /**< green part */
   int*                  blue,               /**< blue part */
   SCIP_Real             scale               /**< maximal coefficient */
   )
{
   SCIP_Real coeflog;

   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(readerdata->rgb_limit >= 0);
   assert(coef > 0);

   coeflog = SCIPfloor(scip, log10(coef));

   if( !(readerdata->rgb_relative) )
   {
      (*red) = 255;
      (*blue) = readerdata->rgb_limit - (int) (unsigned short) (coef/scale * readerdata->rgb_limit);
      (*green) = *blue;
   }
   else
   {
      if( coeflog >= 0 )
      {
         (*red) = 255;
         if( coeflog >= readerdata->coef_limit )
         {
            (*blue) = 0;
            (*green) = 0;
         }
         else
         {
            (*blue) = readerdata->rgb_limit - (int) (unsigned short) (readerdata->rgb_limit * coeflog/readerdata->coef_limit);
            (*green) = *blue;
         }
      }
      else
      {
         (*blue) = 255;
         coeflog = -1.0*coeflog;
         if( coeflog >= readerdata->coef_limit )
         {
            (*red) = 0;
            (*green) = 0;
         }
         else
         {
            (*red) = (readerdata->rgb_limit) - (int) (unsigned short) ((readerdata->rgb_limit)*coeflog/(readerdata->coef_limit));
            (*green) = *red;
         }
      }
   }
}


/** print row in PPM format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_VAR**            vars,               /**< array of constraint variables */
   SCIP_Real*            vals,               /**< array of constraint values */
   int                   nvars,              /**< number of constraint variables */
   int                   ntotalvars,         /**< number of variables */
   SCIP_Real             maxcoef             /**< maximal coefficient */
   )
{
   int v;
   int i;
   int j;

   int red;
   int green;
   int blue;

   char linebuffer[PPM_MAX_LINELEN];
   int linecnt;
   int varindex;
   int actvarindex;
   int maxvarindex;
   int indexvar = 0;

   char buffer[PPM_MAX_LINELEN];
   const unsigned char max = (unsigned char)255;
   char white[4];

   assert( scip != NULL );
   assert (nvars > 0);
   assert (readerdata != NULL);

   i = 0;
   varindex = -1;
   maxvarindex = 0;

   (void) SCIPsnprintf(white, 4, "%c%c%c", max, max, max);
   clearLine(linebuffer, &linecnt);

   /* calculate maximum index of the variables in this constraint */
   for( v = 0; v < nvars; ++v )
   {
      if( maxvarindex < SCIPvarGetProbindex(vars[v]) )
         maxvarindex = SCIPvarGetProbindex(vars[v]);
   }
   assert(maxvarindex < ntotalvars);

   /* print coefficients */
   for(v = 0; v < nvars; ++v)
   {
      actvarindex = maxvarindex;
      for(j = 0; j < nvars; ++j)
      {
         if( varindex < SCIPvarGetProbindex(vars[j]) && SCIPvarGetProbindex(vars[j]) <= actvarindex )
         {
            actvarindex = SCIPvarGetProbindex(vars[j]);
            indexvar = j;
         }
      }
      varindex = actvarindex;

      /* fill in white points since these variables indices do not exits in this constraint */
      for( ; i < varindex; ++i )
      {
         if(readerdata->rgb_ascii)
            appendLine(scip, file, readerdata, linebuffer, &linecnt, white);
         else
            appendLine(scip, file, readerdata, linebuffer, &linecnt, " 255 255 255 ");
      }

      calcColorValue(scip, readerdata, REALABS(vals[indexvar]), &red, &green, &blue, maxcoef);
      if( readerdata->rgb_ascii )
      {
         if( red == 35 || red == 0 )
            red++;
         if( green==35 || green == 0 )
            green++;
         if( blue==35 || blue == 0 )
            blue++;
         (void) SCIPsnprintf(buffer, PPM_MAX_LINELEN, "%c%c%c", (unsigned char)red, (unsigned char)green, (unsigned char)blue);
      }
      else
         (void) SCIPsnprintf(buffer, PPM_MAX_LINELEN, " %d %d %d ", red, green, blue);

      appendLine(scip, file, readerdata, linebuffer, &linecnt, buffer);
      i++;
   }

   /* fill in white points since these variables indices do not exits in this constraint */
   for( ; i < ntotalvars; ++i )
   {
      if(readerdata->rgb_ascii)
         appendLine(scip, file, readerdata, linebuffer, &linecnt, white);
      else
         appendLine(scip, file, readerdata, linebuffer, &linecnt, " 255 255 255 ");
   }

   endLine(scip, file, readerdata, linebuffer, &linecnt);
}


/** prints given linear constraint information in PPM format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   int                   ncompletevars,      /**< number of variables in whole problem */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Real*            maxcoef,            /**< maximal coefficient */
   SCIP_Bool             printbool           /**< print row or calculate maximum coefficient */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( readerdata != NULL );

   /* duplicate variable and value array */
   nactivevars = nvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   if( vals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

   if( ! readerdata->rgb_relative )
   {
      if( ! printbool )
      {
         for(v = 0; v < nactivevars; ++v)
         {
            if( REALABS(activevals[v]) > *maxcoef)
               *maxcoef = REALABS(activevals[v]);
         }
      }
      else
      {
         assert (*maxcoef > 0);
         /* print constraint */
         printRow(scip, file, readerdata, activevars, activevals, nactivevars, ncompletevars, *maxcoef);
      }
   }
   else
   {
      /* print constraint */
      printRow(scip, file, readerdata, activevars, activevals, nactivevars, ncompletevars, *maxcoef);
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyPpm)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderPpm(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreePpm)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWritePpm)
{  /*lint --e{715}*/

   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPwritePpm(scip, file, name, readerdata, transformed, vars, nvars, conss, nconss, result) );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the ppm file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderPpm(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create ppm reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );
   initReaderdata(readerdata);

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyPpm) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreePpm) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWritePpm) );

   /* add ppm reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/ppmreader/rgbrelativ", "should the coloring values be relativ or absolute",
         &readerdata->rgb_relative, FALSE, DEFAULT_PPM_RGB_RELATIVE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/ppmreader/rgbascii", "should the output format be binary(P6) (otherwise plain(P3) format)",
         &readerdata->rgb_ascii, FALSE, DEFAULT_PPM_RGB_ASCII, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/ppmreader/coefficientlimit",
         "splitting coefficients in this number of intervals",
         &readerdata->coef_limit, FALSE, DEFAULT_PPM_COEF_LIMIT, 3, 16, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/ppmreader/rgblimit",
         "maximal color value",
         &readerdata->rgb_limit, FALSE, DEFAULT_PPM_RGB_LIMIT, 0, 255, NULL, NULL) );

   return SCIP_OKAY;
}


/** writes problem to file */
SCIP_RETCODE SCIPwritePpm(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{  /*lint --e{715}*/
   int c;
   int v;
   int i;

   int linecnt;
   char linebuffer[PPM_MAX_LINELEN];

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;
   int i_max = 1;
   SCIP_Real maxcoef = 0;
   SCIP_Bool printbool = FALSE;

   assert( scip != NULL );
   assert(readerdata != NULL);

   /* print statistics as comment to file */
   if(readerdata->rgb_ascii)
      SCIPinfoMessage(scip, file, "P6\n");
   else
      SCIPinfoMessage(scip, file, "P3\n");
   SCIPinfoMessage(scip, file, "# %s\n", name);
   SCIPinfoMessage(scip, file, "%d %d\n", nvars, nconss);
   SCIPinfoMessage(scip, file, "255\n");

   clearLine(linebuffer, &linecnt);

   if( ! readerdata->rgb_relative )
      i_max = 2;

   for(i = 0; i < i_max; ++i)
   {
      if( i )
      {
         printbool = TRUE;
         SCIPdebugMsgPrint(scip, "Maximal coefficient = %g\n", maxcoef);
      }

      for(c = 0; c < nconss; ++c)
      {
         cons = conss[c];
         assert( cons != NULL);

         /* in case the transformed is written only constraint are posted which are enabled in the current node */
         assert(!transformed || SCIPconsIsEnabled(cons));

         conshdlr = SCIPconsGetHdlr(cons);
         assert( conshdlr != NULL );

         conshdlrname = SCIPconshdlrGetName(conshdlr);
         assert( transformed == SCIPconsIsTransformed(cons) );

         if( strcmp(conshdlrname, "linear") == 0 )
         {
            consvars = SCIPgetVarsLinear(scip, cons);
            nconsvars = SCIPgetNVarsLinear(scip, cons);
            assert( consvars != NULL || nconsvars == 0 );

            if( nconsvars > 0 )
            {
               SCIP_CALL( printLinearCons(scip, file, readerdata, consvars, SCIPgetValsLinear(scip, cons),
                     nconsvars, nvars, transformed, &maxcoef, printbool) );
            }
         }
         else if( strcmp(conshdlrname, "setppc") == 0 )
         {
            consvars = SCIPgetVarsSetppc(scip, cons);
            nconsvars = SCIPgetNVarsSetppc(scip, cons);
            assert( consvars != NULL || nconsvars == 0 );

            if( nconsvars > 0 )
            {
               SCIP_CALL( printLinearCons(scip, file, readerdata, consvars, NULL,
                     nconsvars, nvars, transformed, &maxcoef, printbool) );
            }
         }
         else if( strcmp(conshdlrname, "logicor") == 0 )
         {
            consvars = SCIPgetVarsLogicor(scip, cons);
            nconsvars = SCIPgetNVarsLogicor(scip, cons);
            assert( consvars != NULL || nconsvars == 0 );

            if( nconsvars > 0 )
            {
               SCIP_CALL( printLinearCons(scip, file, readerdata, consvars, NULL,
                     nconsvars, nvars, transformed, &maxcoef, printbool) );
            }
         }
         else if( strcmp(conshdlrname, "knapsack") == 0 )
         {
            SCIP_Longint* weights;

            consvars = SCIPgetVarsKnapsack(scip, cons);
            nconsvars = SCIPgetNVarsKnapsack(scip, cons);
            assert( consvars != NULL || nconsvars == 0 );

            /* copy Longint array to SCIP_Real array */
            weights = SCIPgetWeightsKnapsack(scip, cons);
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
            for( v = 0; v < nconsvars; ++v )
               consvals[v] = (SCIP_Real)weights[v];

            if( nconsvars > 0 )
            {
               SCIP_CALL( printLinearCons(scip, file, readerdata, consvars, consvals, nconsvars, nvars, transformed, &maxcoef, printbool) );
            }

            SCIPfreeBufferArray(scip, &consvals);
         }
         else if( strcmp(conshdlrname, "varbound") == 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

            consvars[0] = SCIPgetVarVarbound(scip, cons);
            consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

            consvals[0] = 1.0;
            consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

            SCIP_CALL( printLinearCons(scip, file, readerdata, consvars, consvals, 2, nvars, transformed, &maxcoef, printbool) );

            SCIPfreeBufferArray(scip, &consvars);
            SCIPfreeBufferArray(scip, &consvals);
         }
         else
         {
            SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
            SCIPinfoMessage(scip, file, "\\ ");
            SCIP_CALL( SCIPprintCons(scip, cons, file) );
            SCIPinfoMessage(scip, file, ";\n");
         }
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
