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

/**@file   reader_pbm.c
 * @brief  file writer for portable bitmap file format (PBM), open with common graphic viewer programs (e.g. xview)
 * @author Alexandra Kraft
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "scip/reader_pbm.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/pub_misc.h"

#define READER_NAME             "pbmreader"
#define READER_DESC             "file writer for portable bitmap file format (PBM), open with common graphic viewer programs (e.g. xview)"
#define READER_EXTENSION        "pbm"

/*
 * Data structures
 */
#define PBM_MAX_LINELEN               71      /**< the maximum length of any line is 70 + '\\0' = 71*/
#define DEFAULT_PBM_BINARY            TRUE    /**< binary is the default format for PBM */
#define DEFAULT_PBM_MAXROWS           1000    /**< allowed maximum of pixel-rows int the picture */
#define DEFAULT_PBM_MAXCOLS           1000    /**< allowed maximum of pixel-columns in the picture */

/** LP reading data */
struct SCIP_ReaderData
{
   SCIP_Bool             binary;             /**< binary is the default format for PBM */
   int                   maxrows;            /**< allowed maximum of pixel-rows int the picture */
   int                   maxcols;            /**< allowed maximum of pixel-columns in the picture */
};

/*
 * Local methods (for writing)
 */

/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(scalars != NULL);
   assert(nvars != NULL);
   assert(constant != NULL);

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &scalars, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
         assert(requiredsize <= *nvars);
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

/** transforms given variables to the corresponding active variables */
static
SCIP_RETCODE getActiveVariables2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars != NULL);

   if( transformed )
   {
      SCIP_CALL( SCIPgetActiveVars(scip, vars, nvars, *nvars, &requiredsize) );

      if( requiredsize > *nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );

         SCIP_CALL( SCIPgetActiveVars(scip, vars, nvars, requiredsize, &requiredsize) );
         assert(requiredsize <= *nvars);
      }
   }
   else
   {
      SCIP_Real scalar;
      SCIP_Real constant;
      for( v = 0; v < *nvars; ++v )
      {
         scalar = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalar, &constant) );
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
   assert(linebuffer != NULL);
   assert(linecnt != NULL);

   (*linecnt) = 0;
   linebuffer[0] = '\0';
}


/** appends a bit to buffer and prints it to the give file stream if we've gather a whole byte */
static
void flushBits(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   unsigned char*        bitcnt,             /**< counts bits until whole byte is gathered */
   unsigned char*        bitbuffer           /**< bit buffer */
   )
{
   assert(scip != NULL);

   if( *bitcnt == 0 )
      return;

   assert(bitbuffer != NULL);
   assert(*bitcnt > 0 && *bitcnt <= 8);

   (*bitbuffer) <<= (8 - *bitcnt); /*lint !e701 !e734*/

   fputc(*bitbuffer, file);

   *bitcnt = 0;
   *bitbuffer = 0;
}


/** appends a bit to buffer and prints it to the given file stream if we've gathered a whole byte */
static
void appendBit(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   unsigned char         bit,                /**< bit to append */
   unsigned char*        bitcnt,             /**< counts bits until whole byte is gathered */
   unsigned char*        bitbuffer           /**< bit buffer */
   )
{
   assert(scip != NULL);
   assert(*bitcnt < 8);
   assert(bitbuffer != NULL);

   (*bitbuffer) = ((*bitbuffer)<<1)|(bit&1); /*lint !e734*/
   *bitcnt += 1;

   if( *bitcnt == 8 )
      flushBits(scip, file, bitcnt, bitbuffer);
}

/** calculates the size of the quadratic matrix, which will correspond to one pixel in the picture */
static
int getSubmatrixSize(
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   int                   nvars,              /**< number of variables */
   int                   nconss              /**< number of constraints */
   )
{
   int sizev;
   int sizeh;
   int res;

   assert(readerdata->maxrows != 0 && readerdata->maxcols != 0);

   if( readerdata->maxrows > nconss )
      readerdata->maxrows = nconss;

   if( readerdata->maxcols > nvars )
      readerdata->maxcols = nvars;

   sizev = (nconss + readerdata->maxrows - 1) / readerdata->maxrows;
   sizeh = (nvars + readerdata->maxcols - 1) / readerdata->maxcols;

   /* both defined with -1 */
   if( readerdata->maxrows == -1 && readerdata->maxcols == -1 )
      res = 1;

   /* only width is defined */
   else if( readerdata->maxrows == -1 && readerdata->maxcols > 0 )
      res = sizeh;

   /* only height is defined */
   else if( readerdata->maxrows > 0 && readerdata->maxcols == -1 )
      res = sizev;

   /* both are defined, use smaller scaling factor */
   else if( sizev > sizeh )
      res = sizev;
   else
      res = sizeh;

   readerdata->maxrows = (nconss + res - 1) / res;
   readerdata->maxcols = (nvars + res - 1) / res;

   return res;
}


/** print row in PBM format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_VAR**            vars,               /**< array of constraint variables */
   int                   conscnt,            /**< current constraint */
   int                   nvars,              /**< number of constraint variables */
   int                   submatrixsize,      /**< size of the submatrices */
   int*                  scaledimage         /**< array of ints that count variables in every submatrix */
   )
{
   int v;
   int i;
   const int y = conscnt / submatrixsize;
   int x;

   assert(scip != NULL);
   assert(nvars > 0);
   assert(readerdata != NULL);

   for( i = 0; i < nvars; ++i )
   {
      v = SCIPvarGetProbindex(vars[i]);
      if( v != -1 )
      {
         x = v / submatrixsize;
         ++(scaledimage[y * readerdata->maxcols + x]);
      }
   }

   return;
}

/** prints given linear constraint information in PBM format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< current constraint */
   int                   conscnt,            /**< counts variables in the constraint */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   int                   submatrixsize,      /**< size of the submatrices */
   int*                  scaledimage         /**< array of ints that count variables in every submatrix */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   SCIP_Real activeconstant = 0.0;
   int nactivevars;
   int v;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(conscnt >= 0);
   assert(readerdata != NULL);

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

   /* print constraint */
   printRow(scip, readerdata, activevars, conscnt, nactivevars, submatrixsize, scaledimage);

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

/* draws the scaled image */
static
void drawScaledImage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   int*                  scaledimage         /**< array of ints that count variables in every submatrix */
   )
{
   int y;
   int x;
   unsigned char bitcnt = 0;
   unsigned char bitbuffer = '\0';

   assert(scip != NULL);
   assert(readerdata != NULL);


   for( y = 0; y < readerdata->maxrows; y++ )
   {
      for( x = 0; x < readerdata->maxcols; x++ )
      {
         unsigned char v = 0;
         if( scaledimage[y*readerdata->maxcols+x] >= 1 )
         {
            v = 1;
         }
         appendBit(scip, file, v, &bitcnt, &bitbuffer);
      }
      flushBits(scip, file, &bitcnt, &bitbuffer);
   }
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyPbm)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderPbm(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreePbm)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** problem reading method of reader */
#define readerReadPbm NULL

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWritePbm)
{  /*lint --e{715}*/

   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( SCIPwritePbm(scip, file, name, readerdata, transformed, nvars, conss, nconss, result) );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the pbm file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderPbm(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create pbm reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );

   /* include pbm reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyPbm, readerFreePbm, readerReadPbm, readerWritePbm, readerdata) );

   /* add pbm reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/pbmreader/binary", "should the output format be binary(P4) (otherwise plain(P1) format)",
         &readerdata->binary, FALSE, DEFAULT_PBM_BINARY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/pbmreader/maxrows", "maximum number of rows in the scaled picture (-1 for no limit)",
         &readerdata->maxrows, FALSE, DEFAULT_PBM_MAXROWS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "reading/pbmreader/maxcols", "maximum number of columns in the scaled picture (-1 for no limit)",
         &readerdata->maxcols, FALSE, DEFAULT_PBM_MAXCOLS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/* writes picture of matrix structure to file */
SCIP_RETCODE SCIPwritePbm(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_READERDATA*      readerdata,         /**< information for reader */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_READERDATA readerdata_copy;
   char linebuffer[PBM_MAX_LINELEN];
   const char* conshdlrname;
   int* scaledimage;
   int submatrixsize;
   int nconsvars;
   int linecnt;
   int c;
   int v;

   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(conss != NULL);

   readerdata_copy = *readerdata;
   submatrixsize = getSubmatrixSize(&readerdata_copy, nvars, nconss);
   readerdata = &readerdata_copy;

   SCIP_CALL( SCIPallocBufferArray(scip, &scaledimage, readerdata_copy.maxrows * readerdata_copy.maxcols) );
   assert(scaledimage != NULL);
   BMSclearMemoryArray(scaledimage, readerdata->maxrows * readerdata->maxcols);

   /* print statistics as comment to file */
   if( readerdata->binary )
      SCIPinfoMessage(scip, file, "P4\n");
   else
      SCIPinfoMessage(scip, file, "P1\n");

   SCIPinfoMessage(scip, file, "# %s\n", name);
   SCIPinfoMessage(scip, file, "%d %d\n", readerdata->maxcols, readerdata->maxrows);

   clearLine(linebuffer, &linecnt);

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);

      /* in case the transformed is written only constraint are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert(transformed == SCIPconsIsTransformed(cons));

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         consvars = SCIPgetVarsLinear(scip, cons);
         nconsvars = SCIPgetNVarsLinear(scip, cons);
         assert(consvars != NULL || nconsvars == 0);

         if( nconsvars > 0 )
         {
            SCIP_CALL( printLinearCons(scip, readerdata, consvars, SCIPgetValsLinear(scip, cons),
                  nconsvars, c, transformed, submatrixsize, scaledimage) );
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);
         assert(consvars != NULL || nconsvars == 0);

         if( nconsvars > 0 )
         {
            SCIP_CALL( printLinearCons(scip, readerdata, consvars, NULL,
                  nconsvars, c, transformed, submatrixsize, scaledimage) );
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         consvars = SCIPgetVarsLogicor(scip, cons);
         nconsvars = SCIPgetNVarsLogicor(scip, cons);
         assert(consvars != NULL || nconsvars == 0);

         if( nconsvars > 0 )
         {
            SCIP_CALL( printLinearCons(scip, readerdata, consvars, NULL,
                  nconsvars, c, transformed, submatrixsize, scaledimage) );
         }
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);
         assert(consvars != NULL || nconsvars == 0);

         /* copy Longint array to SCIP_Real array */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
            consvals[v] = weights[v];

         if( nconsvars > 0 )
         {
            SCIP_CALL( printLinearCons(scip, readerdata, consvars, consvals, nconsvars, c, transformed,
                  submatrixsize, scaledimage) );
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

         SCIP_CALL( printLinearCons(scip, readerdata, consvars, consvals, 2, c, transformed,
               submatrixsize, scaledimage) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else
      {
         SCIP_Bool success;

         consvars = NULL;
         SCIP_CALL( SCIPgetConsNVars(scip, cons, &nconsvars, &success) );

         if( success )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );

            SCIP_CALL( SCIPgetConsVars(scip, cons, consvars, nconsvars, &success) );
         }

         if( success )
         {
            /* retransform given variables to active variables */
            SCIP_CALL( getActiveVariables2(scip, consvars, &nconsvars, transformed) );

            printRow(scip, readerdata, consvars, c, nconsvars, submatrixsize, scaledimage);
         }
         else
         {
            SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
            SCIPinfoMessage(scip, file, "\\ ");
            SCIP_CALL( SCIPprintCons(scip, cons, file) );
         }

         SCIPfreeBufferArrayNull(scip, &consvars);
      }
   }

   drawScaledImage(scip, file, readerdata, scaledimage);

   SCIPfreeBufferArray(scip, &scaledimage);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
