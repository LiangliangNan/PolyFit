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

/**@file   reader_gms.c
 * @brief  GAMS file writer
 * @author Ambros Gleixner
 * @author Stefan Vigerske
 *
 * @todo Check for words reserved for GAMS.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#ifdef WITH_GAMS
#include <sys/stat.h>

#include "gmomcc.h"
#include "gevmcc.h"

#include "reader_gmo.h"
#endif

#include "scip/reader_gms.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_soc.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_indicator.h"
#include "scip/cons_abspower.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_bivariate.h"
#include "scip/pub_misc.h"

#define READER_NAME             "gmsreader"
#ifdef WITH_GAMS
#define READER_DESC             "file writer for MI(NL)(SOC)Ps in GAMS file format"
#else
#define READER_DESC             "file reader and writer for MI(NL)(SOC)Ps in GAMS file format"
#endif
#define READER_EXTENSION        "gms"


#define GMS_MAX_LINELEN      256
#define GMS_MAX_PRINTLEN     256       /**< the maximum length of any line is 255 + '\\0' = 256*/
#define GMS_MAX_NAMELEN      64        /**< the maximum length for any name is 63 + '\\0' = 64 */
#define GMS_PRINTLEN         100
#define GMS_DEFAULT_BIGM     1e+6
#define GMS_DEFAULT_INDICATORREFORM 's'
#define GMS_DEFAULT_SIGNPOWER FALSE

/*
 * Local methods (for writing)
 */

static const char badchars[] = "#*+/-@$[](){}";

/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< vars array to get active variables for */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
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
      SCIPinfoMessage(scip, file, "%s\n", linebuffer);
      clearLine(linebuffer, linecnt);
   }
}

/** appends extension to line and prints it to the give file stream if the
 *  line exceeded the length given in the define GMS_PRINTLEN */
static
void appendLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of characters in line */
   const char*           extension           /**< string to extend the line */
   )
{
   size_t len;
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( extension != NULL );
   assert( strlen(linebuffer) + strlen(extension) < GMS_MAX_PRINTLEN );

   /* NOTE: avoid
    *   sprintf(linebuffer, "%s%s", linebuffer, extension); 
    * because of overlapping memory areas in memcpy used in sprintf.
    */
   len = strlen(linebuffer);
   strncat(linebuffer, extension, GMS_MAX_PRINTLEN - len);

   (*linecnt) += (int) strlen(extension);

   SCIPdebugMsg(scip, "linebuffer <%s>, length = %lu\n", linebuffer, (unsigned long)len);

   if( (*linecnt) > GMS_PRINTLEN )
      endLine(scip, file, linebuffer, linecnt);
}

/** appends extension to line and prints it to the give file stream if the
 *  line exceeded the length given in the define GMS_PRINTLEN
 *  indents the line by some spaces if it is a new line */
static
void appendLineWithIndent(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of characters in line */
   const char*           extension           /**< string to extend the line */
   )
{
   if( *linecnt == 0 )
      /* we start a new line; therefore we indent line */
      appendLine(scip, file, linebuffer, linecnt, "     ");

   appendLine(scip, file, linebuffer, linecnt, extension);
}

/** checks string for occurences of bad symbols and replace those by '_' */
static
void conformName(
   char*                 name                /**< string to adjust */
   )
{
   const char* badchar;

   assert( name != NULL );

   for( badchar = badchars; *badchar; ++badchar )
   {
      char* c = strchr(name, *badchar);

      while( c != NULL )
      {
         assert( *c == *badchar );

         *c = '_';
         c = strchr(c, *badchar);
      }
   }
}

/* print first len-1 characters of name to string s and replace '#', '*', '+', '/', and '-' by '_' if necessary */
static
SCIP_RETCODE printConformName(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 t,                  /**< target string */
   int                   len,                /**< length of t */
   const char*           name                /**< source string or format string */
   )
{
   SCIP_Bool replaceforbiddenchars;

   assert( t != NULL );
   assert( len > 0 );

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/gmsreader/replaceforbiddenchars", &replaceforbiddenchars) );

   (void) SCIPsnprintf(t, len, "%s", name);

   if( replaceforbiddenchars )
      conformName(t);

   return SCIP_OKAY;
}


/* retransform to active variables and print in GAMS format to file stream with surrounding bracket, pre- and suffix */
static
SCIP_RETCODE printActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of characters in line */
   const char*           prefix,             /**< prefix (maybe NULL) */
   const char*           suffix,             /**< suffix (maybe NULL) */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of values (or NULL if all ones) */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int v;
   int closingbracket;

   SCIP_VAR* var;
   char varname[GMS_MAX_NAMELEN];
   char buffer[GMS_MAX_PRINTLEN];
   char ext[GMS_MAX_PRINTLEN];

   SCIP_VAR** activevars = NULL;
   SCIP_Real* activevals = NULL;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( vars != NULL || nvars == 0 );


   if( *linecnt == 0 )
      /* we start a new line; therefore we tab this line */
      appendLine(scip, file, linebuffer, linecnt, "     ");

   if( nvars == 0 )
   {
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s(0)%s", prefix != NULL ? prefix : "", suffix != NULL ? suffix : "");

      appendLine(scip, file, linebuffer, linecnt, buffer);
   }
   else
   {
      nactivevars = nvars;

      /* duplicate variable and value array */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars) );
      if( vals != NULL )
      {
         SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars) );
      }
      else
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

         for( v = 0; v < nactivevars; ++v )
            activevals[v] = 1.0;
      }

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

      assert( nactivevars == 0 || activevals != NULL );

      if( nactivevars == 0 && SCIPisZero(scip, activeconstant) )
      {
         if( *linecnt == 0 )
            /* we start a new line; therefore we tab this line */
            appendLine(scip, file, linebuffer, linecnt, "     ");

         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s(0)%s", prefix != NULL ? prefix : "", suffix != NULL ? suffix : "");

         appendLine(scip, file, linebuffer, linecnt, buffer);
      }
      else
      {
         /* buffer prefix */
         (void) SCIPsnprintf(ext, GMS_MAX_PRINTLEN, "%s(", prefix != NULL ? prefix : "");

         /* find position of closing bracket */
         closingbracket = nactivevars;
         if( SCIPisZero(scip, activeconstant) )
         {
            do
               --closingbracket;
            while( SCIPisZero(scip, activevals[closingbracket]) && closingbracket > 0 );
         }

         /* print active variables */
         for( v = 0; v < nactivevars; ++v )
         {
            var = activevars[v];
            assert( var != NULL );

            if( !SCIPisZero(scip, activevals[v]) )
            {
               if( *linecnt == 0 )
                  /* we start a new line; therefore we tab this line */
                  appendLine(scip, file, linebuffer, linecnt, "     ");

               SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(var)) );

               if( SCIPisEQ(scip, activevals[v], 1.0) )
                  (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s%s%s%s%s", ext, strchr(ext, '(') == NULL ? "+" : "",
                        varname, (v == closingbracket) ? ")" : "", (v == closingbracket && suffix) ? suffix : "");
               else if( SCIPisEQ(scip, activevals[v], -1.0) )
                  (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s-%s%s%s", ext,
                        varname, (v == closingbracket) ? ")" : "", (v == closingbracket && suffix) ? suffix : "");
               else if( strchr(ext, '(') != NULL )
                  (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s%.15g*%s%s%s", ext,
                        activevals[v], varname, (v == closingbracket) ? ")" : "", (v == closingbracket && suffix) ? suffix : "");
               else
                  (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s%+.15g*%s%s%s", ext,
                        activevals[v], varname, (v == closingbracket) ? ")" : "", (v == closingbracket && suffix) ? suffix : "");

               appendLine(scip, file, linebuffer, linecnt, buffer);

               (void) SCIPsnprintf(ext, GMS_MAX_PRINTLEN, (*linecnt == 0) ? "" : " ");
            }
         }

         /* print active constant */
         if( !SCIPisZero(scip, activeconstant) )
         {
            if( *linecnt == 0 )
               /* we start a new line; therefore we tab this line */
               appendLine(scip, file, linebuffer, linecnt, "     ");

            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s%+.15g)%s", ext, activeconstant, suffix ? suffix : "");

            appendLine(scip, file, linebuffer, linecnt, buffer);
         }
         /* nothing has been printed, yet */
         else if( strchr(ext, '(') != NULL )
         {
            if( *linecnt == 0 )
               /* we start a new line; therefore we tab this line */
               appendLine(scip, file, linebuffer, linecnt, "     ");

            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s(0)%s", prefix ? prefix : "", suffix ? suffix : "");

            appendLine(scip, file, linebuffer, linecnt, buffer);
         }
      }

      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &activevars);
      SCIPfreeBufferArray(scip, &activevals);
   }

   return SCIP_OKAY;
}


/* print linear row in GAMS format to file stream (without retransformation to active variables) */
static
SCIP_RETCODE printLinearRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=e=", "=l=", or "=g=") */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of values */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   int v;
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   SCIP_VAR* var;
   char varname[GMS_MAX_NAMELEN];
   char consname[GMS_MAX_NAMELEN + 3]; /* four extra characters for ' ..' */
   char buffer[GMS_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strcmp(type, "=e=") == 0 || strcmp(type, "=l=") == 0 || strcmp(type, "=g=") == 0);
   assert( nvars == 0 || (vars != NULL && vals != NULL) );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   if( strlen(rowname) > 0 || strlen(rownameextension) > 0 )
   {
      (void) SCIPsnprintf(buffer, GMS_MAX_NAMELEN + 3, "%s%s ..", rowname, rownameextension);
      SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN + 3, buffer) );
      appendLine(scip, file, linebuffer, &linecnt, consname);
   }

   /* print coefficients */
   if( nvars == 0 )
   {
      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, "     ");

      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " 0");

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );

      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, "     ");

      SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(var)) );
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %+.15g*%s", vals[v], varname);

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print right hand side */
   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s %.15g;", type, rhs);

   /* we start a new line; therefore we tab this line */
   if( linecnt == 0 )
      appendLine(scip, file, linebuffer, &linecnt, "     ");
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}


/** prints given linear constraint information in GAMS format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int v;
   SCIP_VAR** activevars = NULL;
   SCIP_Real* activevals = NULL;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;

   assert( scip != NULL );
   assert( rowname != NULL );

   /* The GAMS format does not forbid that the variable array is empty */
   assert( nvars == 0 || vars != NULL );

   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   nactivevars = nvars;
   if( nvars > 0 ) 
   {
      /* duplicate variable and value array */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars) );
      if( vals != NULL )
      {
         SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars) );
      }
      else
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

         for( v = 0; v < nactivevars; ++v )
            activevals[v] = 1.0;
      }

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );
   }

   /* print row(s) in GAMS format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* print equality constraint */
      SCIP_CALL( printLinearRow(scip, file, rowname, "", "=e=", 
            nactivevars, activevars, activevals, rhs - activeconstant) );
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printLinearRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", "=g=",
               nactivevars, activevars, activevals, lhs - activeconstant) );
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         SCIP_CALL( printLinearRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "=l=",
               nactivevars, activevars, activevals, rhs - activeconstant) );
      }
   }

   if( nvars > 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &activevars);
      SCIPfreeBufferArray(scip, &activevals);
   }

   return SCIP_OKAY;
}


/* print quadratic row in GAMS format to file stream (performing retransformation to active variables) */
static
SCIP_RETCODE printQuadraticRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=e=", "=l=", or "=g=") */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< variables in linear part */ 
   SCIP_Real*            lincoeffs,          /**< coefficients of variables in linear part */ 
   int                   nquadvarterms,      /**< number of quadratic variable terms */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int t;
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   SCIP_VAR* var;
   char consname[GMS_MAX_NAMELEN + 3]; /* four extra characters for ' ..' */
   char buffer[GMS_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strlen(rowname) > 0 || strlen(rownameextension) > 0 );
   assert( strcmp(type, "=e=") == 0 || strcmp(type, "=l=") == 0 || strcmp(type, "=g=") == 0 );
   assert( nlinvars == 0 || (linvars != NULL && lincoeffs != NULL) );
   assert( nquadvarterms == 0 || quadvarterms != NULL );
   assert( nbilinterms == 0 || bilinterms != NULL );
   assert( nquadvarterms > 0 || nbilinterms == 0 );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   (void) SCIPsnprintf(buffer, GMS_MAX_NAMELEN + 3, "%s%s ..", rowname, rownameextension);
   SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN + 3, buffer) );

   appendLine(scip, file, linebuffer, &linecnt, consname);

   /* print linear terms */
   if( nlinvars > 0 )
   {
      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", " ", nlinvars, linvars, lincoeffs, transformed) );
   }

   /* print linear coefficients of quadratic terms */
   for( t = 0; t < nquadvarterms; ++t )
   {
      var = quadvarterms[t].var;
      assert( var != NULL );

      if( !SCIPisZero(scip, quadvarterms[t].lincoef) )
      {
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%+.15g*", quadvarterms[t].lincoef);

         SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, buffer, NULL, 1, &var, NULL, transformed) );
      }
   }

   /* print square coefficients of quadratic terms */
   for( t = 0; t < nquadvarterms; ++t )
   {
      var = quadvarterms[t].var;
      assert( var != NULL );

      if( !SCIPisZero(scip, quadvarterms[t].sqrcoef) )
      {
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%+.15g*sqr", quadvarterms[t].sqrcoef);

         SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, buffer, NULL, 1, &var, NULL, transformed) );
      }
   }

   /* print bilinear terms */
   for( t = 0; t < nbilinterms; ++t )
   {
      if( !SCIPisZero(scip, bilinterms[t].coef) )
      {
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%+.15g*", bilinterms[t].coef);

         /* print first variable (retransformed to active variables) */
         var = bilinterms[t].var1;
         assert( var != NULL );

         SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, buffer, "", 1, &var, NULL, transformed) );

         /* print second variable (retransformed to active variables) */
         var = bilinterms[t].var2;
         assert( var != NULL );

         SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "*", " ", 1, &var, NULL, transformed) );
      }
   }

   /* print right hand side */
   if( linecnt == 0 )
      /* we start a new line; therefore we tab this line */
      appendLine(scip, file, linebuffer, &linecnt, "     ");

   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s%s %.15g;", (nlinvars == 0 && nquadvarterms == 0) ? "0 " : "", type, rhs);

   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}


/** prints given quadratic constraint information in GAMS format to file stream */
static
SCIP_RETCODE printQuadraticCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< variables in linear part */ 
   SCIP_Real*            lincoeffs,          /**< coefficients of variables in linear part */ 
   int                   nquadvarterms,      /**< number of quadratic variable terms */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   assert( scip != NULL );
   assert( rowname != NULL );
   assert( nlinvars == 0 || (linvars != NULL && lincoeffs != NULL) );
   assert( nquadvarterms == 0 || quadvarterms != NULL );
   assert( nbilinterms == 0 || bilinterms != NULL );
   assert( nquadvarterms > 0 || nbilinterms == 0 );
   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   /* print row(s) in GAMS format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* print equality constraint */
      SCIP_CALL( printQuadraticRow(scip, file, rowname, "", "=e=",
         nlinvars, linvars, lincoeffs,
         nquadvarterms, quadvarterms,
         nbilinterms, bilinterms, rhs, transformed) );
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printQuadraticRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", "=g=",
            nlinvars, linvars, lincoeffs,
            nquadvarterms, quadvarterms,
            nbilinterms, bilinterms, lhs, transformed) );
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         SCIP_CALL( printQuadraticRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "=l=",
            nlinvars, linvars, lincoeffs,
            nquadvarterms, quadvarterms,
            nbilinterms, bilinterms, rhs, transformed) );
      }
   }

   return SCIP_OKAY;
}

/** check GAMS limitations on SOC constraints
 * returns true of constraint can be written as conic equation in GAMS (using equation type =C=)
 */
static
SCIP_Bool isGAMSprintableSOC(
   int                   nlhsvars,           /**< number of variables on left hand side */
   SCIP_VAR**            lhsvars,            /**< variables on left hand side */
   SCIP_Real*            lhscoeffs,          /**< coefficients of variables on left hand side, or NULL if == 1.0 */
   SCIP_Real*            lhsoffsets,         /**< offsets of variables on left hand side, or NULL if == 0.0 */
   SCIP_Real             lhsconstant,        /**< constant on left hand side */
   SCIP_VAR*             rhsvar,             /**< variable on right hand side */
   SCIP_Real             rhscoef,            /**< coefficient of variable on right hand side */
   SCIP_Real             rhsoffset           /**< offset of variable on right hand side */
   )
{
   int i;

   assert(nlhsvars == 0 || lhsvars != NULL);

   if( rhscoef != 1.0 )
      return FALSE;

   if( rhsoffset != 0.0 )
      return FALSE;

   if( rhsvar == NULL )
      return FALSE;

   if( !SCIPvarIsActive(rhsvar) )
      return FALSE;

   if( lhsconstant != 0.0 )
      return FALSE;

   if( nlhsvars < 2 )
      return FALSE;

   for( i = 0; i < nlhsvars; ++i )
   {
      if( lhscoeffs [i] != 1.0 )
         return FALSE;

      if( lhsoffsets[i] != 0.0 )
         return FALSE;

      if( !SCIPvarIsActive(lhsvars[i]) )
         return FALSE;
   }

   return TRUE;
}

/* print second order cone row in GAMS format to file stream (performing retransformation to active variables)
 * The constraints are of the following form:
 * \f[
 *    \left\{ x \;:\; \sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} \leq \alpha_{n+1}\, (x_{n+1}+\beta_{n+1}) \right\}.
 * \f]
 * */
static
SCIP_RETCODE printSOCCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   int                   nlhsvars,           /**< number of variables on left hand side */
   SCIP_VAR**            lhsvars,            /**< variables on left hand side */
   SCIP_Real*            lhscoeffs,          /**< coefficients of variables on left hand side, or NULL if == 1.0 */
   SCIP_Real*            lhsoffsets,         /**< offsets of variables on left hand side, or NULL if == 0.0 */
   SCIP_Real             lhsconstant,        /**< constant on left hand side */
   SCIP_VAR*             rhsvar,             /**< variable on right hand side */
   SCIP_Real             rhscoef,            /**< coefficient of variable on right hand side */
   SCIP_Real             rhsoffset,          /**< offset of variable on right hand side */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   char consname[GMS_MAX_NAMELEN + 3]; /* four extra characters for ' ..' */
   char buffer[GMS_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strlen(rowname) > 0 );
   assert( nlhsvars == 0 || lhsvars != NULL );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   (void) SCIPsnprintf(buffer, GMS_MAX_NAMELEN + 3, "%s ..", rowname);
   SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN + 3, buffer) );

   appendLine(scip, file, linebuffer, &linecnt, consname);

   if( !isGAMSprintableSOC(nlhsvars, lhsvars, lhscoeffs, lhsoffsets, lhsconstant, rhsvar, rhscoef, rhsoffset) )
   {
      int t;

      /* print right-hand side on left */
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "sqr(%.15g +", rhsoffset);

      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, buffer, ")", 1, &rhsvar, &rhscoef, transformed) );

      appendLine(scip, file, linebuffer, &linecnt, " =g= ");

      /* print left-hand side on right */

      if( lhsconstant != 0.0 )
      {
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g", lhsconstant);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      for( t = 0; t < nlhsvars; ++t )
      {
         assert( lhsvars[t] != NULL );

         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "+ sqr(%.15g * (%.15g + ", lhscoeffs ? lhscoeffs[t] : 1.0, lhsoffsets ? lhsoffsets[t] : 0.0);

         SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, buffer, "))", 1, &lhsvars[t], NULL, transformed) );
      }
   }
   else
   {
      /* print right-hand side on left */
      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", " ", 1, &rhsvar, &rhscoef, transformed) );

      appendLine(scip, file, linebuffer, &linecnt, " =c= ");

      /* print left-hand side on right */
      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", " ", nlhsvars, lhsvars, lhscoeffs, transformed) );
   }

   appendLine(scip, file, linebuffer, &linecnt, ";");

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/* print indicator constraint in some GAMS format to file stream (performing retransformation to active variables)
 * The constraints are of the following form:
 * \f[
 *    z = 1 -> s = 0
 * \f]
 * */
static
SCIP_RETCODE printIndicatorCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   SCIP_VAR*             z,                  /**< indicating variable (binary) */
   SCIP_VAR*             s,                  /**< slack variable */
   SCIP_Bool*            sossetdeclr,        /**< buffer to store whether we declared the SOS set for indicator reform */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;
   SCIP_Real coef;
   char indicatorform;

   char consname[GMS_MAX_NAMELEN + 30];
   char buffer[GMS_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strlen(rowname) > 0 );
   assert( z != NULL );
   assert( s != NULL );
   assert( SCIPvarIsBinary(z) );
   assert( sossetdeclr != NULL );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   SCIP_CALL( SCIPgetCharParam(scip, "reading/gmsreader/indicatorreform", &indicatorform) );

   switch( indicatorform )
   {
      case 'b':
      {
         /* print row name */
         (void) SCIPsnprintf(buffer, GMS_MAX_NAMELEN + 3, "%s ..", rowname);
         SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN + 3, buffer) );

         appendLine(scip, file, linebuffer, &linecnt, consname);

         /* write as s <= upperbound(s)*(1-z) or s <= upperbound(s) * negation(z) */
         coef = 1.0;
         SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, NULL, " =l= ", 1, &s, &coef, transformed) );

         coef = SCIPvarGetUbGlobal(s);
         if( SCIPisInfinity(scip, coef) )
         {
            SCIP_CALL( SCIPgetRealParam(scip, "reading/gmsreader/bigmdefault", &coef) );

            SCIPwarningMessage(scip, "do not have upper bound on slack variable <%s> in indicator constraint <%s>, will use M = %g.\n",
               SCIPvarGetName(s), rowname, coef);
         }

         if( SCIPvarIsNegated(z) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, z, &z) );
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "", ";", 1, &z, &coef, transformed) );
         }
         else
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g + ", coef);

            coef = -coef;
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, buffer, ";", 1, &z, &coef, transformed) );
         }

         break;
      }

      case 's':
      {
         /* write as
          * sos1 Variable name_sos(sosset);
          *  name_soseq(sosset).. name_sos(sosset) =e= s$(sameas(sosset,'slack') + z$(sameas(sosset,'bin'));
          */
         coef = 1.0;
         SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN, rowname) );

         /* declare set for SOS1 declarations from reformulation of indicator, if needed */
         if( !*sossetdeclr )
         {
            SCIPinfoMessage(scip, file, " Set sosset / slack, bin /;\n");
            *sossetdeclr = TRUE;
         }

         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "sos1 Variable %s_sos(sosset);", consname);
         appendLine(scip, file, linebuffer, &linecnt, buffer);
         endLine(scip, file, linebuffer, &linecnt);

         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s(sosset).. %s_sos(sosset) =e= ", consname, consname);
         appendLine(scip, file, linebuffer, &linecnt, buffer);
         SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, NULL, "$sameas(sosset,'slack')", 1, &s, &coef, transformed) );
         if( SCIPvarIsNegated(z) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, z, &z) );
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, " + (1-(", "))$sameas(sosset,'bin');", 1, &z, &coef, transformed) );
         }
         else
         {
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, " + ", "$sameas(sosset,'bin');", 1, &z, &coef, transformed) );
         }
         endLine(scip, file, linebuffer, &linecnt);

         break;
      }

      default:
         SCIPerrorMessage("wrong value '%c' for parameter reading/gmsreader/indicatorreform\n", indicatorform);
         return SCIP_ERROR;
   }

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/* print SOS constraint in some GAMS format to file stream (performing retransformation to active variables)
 *
 * write as
 * Set name_sosset /1*nvars/;
 * SOS1/2 Variable name_sosvar(name_sosset);
 * Equation name_sosequ(e1_sosset);
 * name_sosequ(name_sosset).. name_sosvar(e1_sosset) =e=
 * vars[0]$sameas(name_sosset, '1') + vars[1]$sameas(name_sosset, '2') + ... + vars[nvars-1]$sameas(name_sosset, nvars);
 */
static
SCIP_RETCODE printSOSCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   int                   nvars,              /**< number of variables in SOS */
   SCIP_VAR**            vars,               /**< variables in SOS */
   int                   sostype,            /**< type of SOS: 1 or 2 */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;
   SCIP_Real coef;
   int v;

   char consname[GMS_MAX_NAMELEN + 30];
   char buffer[GMS_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strlen(rowname) > 0 );
   assert( vars != NULL || nvars == 0 );
   assert( sostype == 1 || sostype == 2 );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN, rowname) );

   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "Set %s_sosset /1*%d/;", consname, nvars);
   appendLine(scip, file, linebuffer, &linecnt, buffer);
   endLine(scip, file, linebuffer, &linecnt);

   /* explicitly set lower bound of SOS variables to -inf, as GAMS default is 0.0 */
   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " SOS%d Variable %s_sosvar(%s_sosset); %s_sosvar.lo(%s_sosset) = -inf;", sostype, consname, consname, consname, consname);
   appendLine(scip, file, linebuffer, &linecnt, buffer);
   endLine(scip, file, linebuffer, &linecnt);

   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s(%s_sosset).. %s_sosvar(%s_sosset) =e= ", consname, consname, consname, consname);
   appendLine(scip, file, linebuffer, &linecnt, buffer);
   endLine(scip, file, linebuffer, &linecnt);

   coef = 1.0;
   for( v = 0; v < nvars; ++v )
   {
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "$sameas(%s_sosset,'%d')", consname, v+1);
      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, v > 0 ? " + " : NULL, buffer, 1, &vars[v], &coef, transformed) ); /*lint !e613*/
   }
   appendLine(scip, file, linebuffer, &linecnt, ";");
   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/* print signpower row in GAMS format to file stream (performing retransformation to active variables) */
static
SCIP_RETCODE printSignpowerRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=e=", "=l=", or "=g=") */
   SCIP_VAR*             nonlinvar,          /**< nonlinear variable */
   SCIP_VAR*             linvar,             /**< linear variable, may be NULL */
   SCIP_Real             exponent,           /**< exponent of nonlinear variable */
   SCIP_Real             offset,             /**< offset of nonlinear variable */
   SCIP_Real             coeflinear,         /**< coefficient of linear variable */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Bool             signpowerallowed,   /**< allowed to use signpower operator in GAMS? */
   SCIP_Bool*            nsmooth             /**< buffer to store whether we printed a nonsmooth function */
   )
{
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;
   SCIP_Bool nisoddint;
   SCIP_Bool fixedsign;

   char consname[GMS_MAX_NAMELEN + 3]; /* four extra characters for ' ..' */
   char buffer[GMS_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strlen(rowname) > 0 || strlen(rownameextension) > 0 );
   assert( strcmp(type, "=e=") == 0 || strcmp(type, "=l=") == 0 || strcmp(type, "=g=") == 0 );
   assert( nonlinvar != NULL );
   assert( exponent > 1.0 );
   assert( nsmooth != NULL );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   (void) SCIPsnprintf(buffer, GMS_MAX_NAMELEN + 3, "%s%s ..", rowname, rownameextension);
   SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN + 3, buffer) );

   appendLine(scip, file, linebuffer, &linecnt, consname);

   /* print nonlinear term
    * if not signpowerallowed, then signpow(x,n) is printed as x*abs(x) if n == 2, x*(abs(x)**(n-1)) if n is not 2 and not an odd integer, and as power(x,n) if n is an odd integer
    * if signpowerallowed, then signpow(x,n) is printed as power(x,n) if n is an odd integer and as signpower(x,n) otherwiser
    */
   nisoddint = SCIPisIntegral(scip, exponent) && ((int)SCIPfloor(scip, exponent+0.5))%2 == 1;
   fixedsign = !SCIPisNegative(scip, SCIPvarGetLbGlobal(nonlinvar)) || !SCIPisPositive(scip, SCIPvarGetUbGlobal(nonlinvar));
   if( !nisoddint && !fixedsign )
   {
      if( signpowerallowed )
      {
         if( offset != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "signpower(%g ", offset);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", ",", 1, &nonlinvar, NULL, transformed) );
         }
         else
         {
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "signpower(", ",", 1, &nonlinvar, NULL, transformed) );
         }
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%g)", exponent);
         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
      else
      {
         if( offset != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "(%g ", offset);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", ") * ", 1, &nonlinvar, NULL, transformed) );
         }
         else
         {
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, NULL, " * ", 1, &nonlinvar, NULL, transformed) );
         }

         if( exponent == 2.0)
         {
            if( offset != 0.0 )
            {
               (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "abs(%g ", offset);
               appendLine(scip, file, linebuffer, &linecnt, buffer);
               SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", ")", 1, &nonlinvar, NULL, transformed) );
            }
            else
            {
               SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "abs", NULL, 1, &nonlinvar, NULL, transformed) );
            }
         }
         else
         {
            if( offset != 0.0 )
            {
               (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "abs(%g ", offset);
               appendLine(scip, file, linebuffer, &linecnt, buffer);
               SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", ")", 1, &nonlinvar, NULL, transformed) );
            }
            else
            {
               SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "abs", NULL, 1, &nonlinvar, NULL, transformed) );
            }
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "**%g", exponent-1.0);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }
      }
      *nsmooth = TRUE;
   }
   else if( nisoddint || !SCIPisNegative(scip, SCIPvarGetLbGlobal(nonlinvar)) )
   {
      if( exponent == 2.0 )
      {
         if( offset != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "sqr(%g ", offset);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", ")", 1, &nonlinvar, NULL, transformed) );
         }
         else
         {
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "sqr", NULL, 1, &nonlinvar, NULL, transformed) );
         }
      }
      else
      {
         if( offset != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "power(%g ", offset);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", ",", 1, &nonlinvar, NULL, transformed) );
         }
         else
         {
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "power(", ",", 1, &nonlinvar, NULL, transformed) );
         }
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%g)", exponent);
         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
   }
   else
   {
      assert(fixedsign && !SCIPisPositive(scip, SCIPvarGetUbGlobal(nonlinvar)));
      if( exponent == 2.0 )
      {
         if( offset != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "-sqr(%g ", -offset);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "-", ")", 1, &nonlinvar, NULL, transformed) );
         }
         else
         {
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "-sqr(-", ")", 1, &nonlinvar, NULL, transformed) );
         }
      }
      else
      {
         if( offset != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "-power(%g ", -offset);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "-", ",", 1, &nonlinvar, NULL, transformed) );
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%g)", exponent);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }
         else
         {
            SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "-power(-", ",", 1, &nonlinvar, NULL, transformed) );
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%g)", exponent);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }
      }
   }

   /* print linear term */
   if( linvar != NULL )
   {
      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, " +", "", 1, &linvar, &coeflinear, transformed) );
   }

   /* print right hand side */
   if( linecnt == 0 )
   {
      /* we start a new line; therefore we tab this line */
      appendLine(scip, file, linebuffer, &linecnt, "     ");
   }

   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s %.15g;", type, rhs);

   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/* print signpower cons in GAMS format to file stream (performing retransformation to active variables)
 */
static
SCIP_RETCODE printSignpowerCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   SCIP_VAR*             nonlinvar,          /**< nonlinear variable */
   SCIP_VAR*             linvar,             /**< linear variable, may be NULL */
   SCIP_Real             exponent,           /**< exponent of nonlinear variable */
   SCIP_Real             offset,             /**< offset of nonlinear variable */
   SCIP_Real             coeflinear,         /**< coefficient of linear variable */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Bool             signpowerallowed,   /**< allowed to use signpower operator in GAMS? */
   SCIP_Bool*            nsmooth             /**< buffer to store whether we printed a nonsmooth function */
   )
{
   assert( scip != NULL );
   assert( strlen(rowname) > 0 );

   /* print row(s) in GAMS format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* print equality constraint */
      SCIP_CALL( printSignpowerRow(scip, file, rowname, "", "=e=",
         nonlinvar, linvar, exponent, offset, coeflinear, rhs, transformed, signpowerallowed, nsmooth) );
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printSignpowerRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", "=g=",
            nonlinvar, linvar, exponent, offset, coeflinear, lhs, transformed, signpowerallowed, nsmooth) );
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         SCIP_CALL( printSignpowerRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "=l=",
            nonlinvar, linvar, exponent, offset, coeflinear, rhs, transformed, signpowerallowed, nsmooth) );
      }
   }

   return SCIP_OKAY;
}

/* prints expression in GAMS format to file stream */
static
SCIP_RETCODE printExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line buffer of length GMS_MAX_PRINTLEN */
   int*                  linecnt,            /**< number of characters in line so far */
   SCIP_Bool*            nsmooth,            /**< buffer to store whether we printed a nonsmooth function */
   SCIP_Bool             transformed,        /**< expression belongs to transformed constraint? */
   SCIP_EXPR*            expr,               /**< expression to print */
   SCIP_VAR**            exprvars            /**< variables of expression */
   )
{
   char buffer[GMS_MAX_PRINTLEN];

   assert(scip != NULL);
   assert(linebuffer != NULL);
   assert(linecnt != NULL);
   assert(expr != NULL);
   assert(nsmooth != NULL);

   switch( SCIPexprGetOperator(expr) )
   {
      case SCIP_EXPR_VARIDX:
      {
         SCIP_Real one;

         assert(exprvars != NULL);

         one = 1.0;
         SCIP_CALL( printActiveVariables(scip, file, linebuffer, linecnt, "",  "",  1, &exprvars[SCIPexprGetOpIndex(expr)], &one, transformed) );

         break;
      }

      case SCIP_EXPR_PARAM:
      {
         SCIPwarningMessage(scip, "parameterized expression in GAMS writer. GAMS file will not compile.\n");

         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "param%d", SCIPexprGetOpIndex(expr));
         appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);

         break;
      }

      case SCIP_EXPR_CONST:
      {
         if( SCIPexprGetOpReal(expr) < 0.0 )
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "(%.15g)", SCIPexprGetOpReal(expr));
         else
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g", SCIPexprGetOpReal(expr));
         appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);

         break;
      }

      case SCIP_EXPR_PLUS:
      {
         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, " + ");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[1], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      case SCIP_EXPR_MINUS:
      {
         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, " - ");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[1], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      case SCIP_EXPR_MUL:
      {
         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, " * ");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[1], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      case SCIP_EXPR_DIV:
      {
         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, " / ");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[1], exprvars) );
         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      case SCIP_EXPR_REALPOWER:
      {
         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, ")**(%.15g)", SCIPexprGetRealPowerExponent(expr));
         appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
         break;
      }

      case SCIP_EXPR_INTPOWER:
      {
         appendLineWithIndent(scip, file, linebuffer, linecnt, "power(");
         SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, ", %d)", SCIPexprGetIntPowerExponent(expr));
         appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
         break;
      }

      case SCIP_EXPR_SIGNPOWER:
      {
         SCIP_Real exponent;
         SCIP_Bool nisoddint;

         /* signpow(x,y) is printed as x*abs(x) if y == 2, x*(abs(x) ** (y-1)) if y is not 2 and not an odd integer, and as intpower(x,y) if y is an odd integer
          * but if reading/gmsreader/signpower is TRUE, then we print as signpower(x,y), unless y is odd integer
          */
         exponent = SCIPexprGetSignPowerExponent(expr);
         nisoddint = (((SCIP_Real)((int)exponent)) == exponent) && (((int)exponent)%2 == 1);

         if( !nisoddint )
         {
            SCIP_Bool signpowerallowed;

            SCIP_CALL( SCIPgetBoolParam(scip, "reading/gmsreader/signpower", &signpowerallowed) );

            if( signpowerallowed )
            {
               appendLineWithIndent(scip, file, linebuffer, linecnt, " * signpower(");
               SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
               (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, ", %.15g)", exponent);
               appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
            }
            else
            {
               appendLineWithIndent(scip, file, linebuffer, linecnt, "(");
               SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
               appendLineWithIndent(scip, file, linebuffer, linecnt, ")");

               if( exponent == 2.0)
               {
                  appendLineWithIndent(scip, file, linebuffer, linecnt, " * abs(");
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
                  appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
               }
               else
               {
                  appendLineWithIndent(scip, file, linebuffer, linecnt, " * abs(");
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
                  (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, ")**(%g)", SCIPexprGetRealPowerExponent(expr)-1.0);
                  appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
               }
            }
            *nsmooth = TRUE;
         }
         else
         {
            appendLineWithIndent(scip, file, linebuffer, linecnt, " * power(");
            SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, ", %.15g)", exponent);
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
         }

         break;
      }

      case SCIP_EXPR_ABS:
      case SCIP_EXPR_SIGN:
         *nsmooth = TRUE; /*lint -fallthrough*/
      case SCIP_EXPR_SQUARE:
      case SCIP_EXPR_SQRT:
      case SCIP_EXPR_EXP:
      case SCIP_EXPR_LOG:
      case SCIP_EXPR_SIN:
      case SCIP_EXPR_COS:
      case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
      case SCIP_EXPR_MIN:
      case SCIP_EXPR_MAX:
      {
         int i;

         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s(", SCIPexpropGetName(SCIPexprGetOperator(expr)));
         appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);

         for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
         {
            SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[i], exprvars) );
            if( i + 1 < SCIPexprGetNChildren(expr) )
               appendLineWithIndent(scip, file, linebuffer, linecnt, ", ");
         }

         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      case SCIP_EXPR_SUM:
      case SCIP_EXPR_PRODUCT:
      {
         switch( SCIPexprGetNChildren(expr) )
         {
            case 0:
            {
               appendLineWithIndent(scip, file, linebuffer, linecnt, SCIPexprGetOperator(expr) == SCIP_EXPR_SUM ? "0" : "1");
               break;
            }
            case 1:
            {
               SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[0], exprvars) );
               break;
            }
            default:
            {
               int i;
               char opstr[GMS_MAX_PRINTLEN];

               (void) SCIPsnprintf(opstr, GMS_MAX_PRINTLEN, SCIPexprGetOperator(expr) == SCIP_EXPR_SUM ? " + " : " * ");
               appendLineWithIndent(scip, file, linebuffer, linecnt, "(");
               for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
               {
                  if( i > 0 )
                  {
                     appendLineWithIndent(scip, file, linebuffer, linecnt, opstr);
                  }
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[i], exprvars) );
               }
               appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
            }
         }
         break;
      }

      case SCIP_EXPR_LINEAR:
      {
         SCIP_Real constant;
         int i;

         constant = SCIPexprGetLinearConstant(expr);

         if( SCIPexprGetNChildren(expr) == 0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g", constant);
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
            break;
         }

         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");

         if( constant != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g", constant);
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
         }

         for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %+.15g * ", SCIPexprGetLinearCoefs(expr)[i]);
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
            SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[i], exprvars) );
         }

         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      case SCIP_EXPR_QUADRATIC:
      {
         SCIP_Real constant;
         int i;
         SCIP_QUADELEM* quadelems;
         SCIP_Real* lincoefs;

         constant = SCIPexprGetQuadConstant(expr);

         if( SCIPexprGetNChildren(expr) == 0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g", constant);
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
            break;
         }

         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");

         if( constant != 0.0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g", constant);
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
         }

         lincoefs = SCIPexprGetQuadLinearCoefs(expr);
         if( lincoefs != NULL )
            for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
            {
               (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %+.15g * ", lincoefs[i]);
               appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
               SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[i], exprvars) );
            }

         quadelems = SCIPexprGetQuadElements(expr);
         for( i = 0; i < SCIPexprGetNQuadElements(expr); ++i )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %+.15g * ", quadelems[i].coef);
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);

            if( quadelems[i].idx1 == quadelems[i].idx2 )
            {
               appendLineWithIndent(scip, file, linebuffer, linecnt, "sqr(");
               SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[quadelems[i].idx1], exprvars) );
               appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
            }
            else
            {
               SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[quadelems[i].idx1], exprvars) );
               appendLineWithIndent(scip, file, linebuffer, linecnt, " * ");
               SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[quadelems[i].idx2], exprvars) );
            }
         }

         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      case SCIP_EXPR_POLYNOMIAL:
      {
         SCIP_EXPRDATA_MONOMIAL* monomdata;
         SCIP_Real exponent;
         int i;
         int j;

         appendLineWithIndent(scip, file, linebuffer, linecnt, "(");

         if( SCIPexprGetPolynomialConstant(expr) != 0.0 || SCIPexprGetNMonomials(expr) == 0 )
         {
            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%.15g", SCIPexprGetPolynomialConstant(expr));
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
         }

         for( i = 0; i < SCIPexprGetNMonomials(expr); ++i )
         {
            monomdata = SCIPexprGetMonomials(expr)[i];
            assert(monomdata != NULL);

            (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %+.15g", SCIPexprGetMonomialCoef(monomdata));
            appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);

            for( j = 0; j < SCIPexprGetMonomialNFactors(monomdata); ++j )
            {
               appendLineWithIndent(scip, file, linebuffer, linecnt, "*");

               exponent = SCIPexprGetMonomialExponents(monomdata)[j];
               if( exponent == 1.0 )
               {
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[SCIPexprGetMonomialChildIndices(monomdata)[j]], exprvars) );
               }
               else if( exponent == 2.0 )
               {
                  appendLineWithIndent(scip, file, linebuffer, linecnt, "sqr(");
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[SCIPexprGetMonomialChildIndices(monomdata)[j]], exprvars) );
                  appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
               }
               else if( exponent == 0.5 )
               {
                  appendLineWithIndent(scip, file, linebuffer, linecnt, "sqrt(");
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[SCIPexprGetMonomialChildIndices(monomdata)[j]], exprvars) );
                  appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
               }
               else if( ((SCIP_Real)((int)exponent)) == exponent )
               {
                  appendLineWithIndent(scip, file, linebuffer, linecnt, "power(");
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[SCIPexprGetMonomialChildIndices(monomdata)[j]], exprvars) );
                  (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, ", %d)", (int)SCIPround(scip, exponent));
                  appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
               }
               else
               {
                  SCIP_CALL( printExpr(scip, file, linebuffer, linecnt, nsmooth, transformed, SCIPexprGetChildren(expr)[SCIPexprGetMonomialChildIndices(monomdata)[j]], exprvars) );
                  (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " ** %.15g", exponent);
                  appendLineWithIndent(scip, file, linebuffer, linecnt, buffer);
               }
            }
         }

         appendLineWithIndent(scip, file, linebuffer, linecnt, ")");
         break;
      }

      default:
         SCIPerrorMessage("unexpected operand %d in expression\n", SCIPexprGetOperator(expr));
         return SCIP_OKAY;
   } /*lint !e788*/

   return SCIP_OKAY;
}

/* print nonlinear row in GAMS format to file stream */
static
SCIP_RETCODE printNonlinearRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=e=", "=l=", or "=g=") */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< variables in linear part */
   SCIP_Real*            lincoeffs,          /**< coefficients of variables in linear part */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            exprtreecoefs,      /**< expression tree coefficients */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Bool*            nsmooth             /**< buffer to store whether we printed a nonsmooth function */
   )
{
   char linebuffer[GMS_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   char consname[GMS_MAX_NAMELEN + 3]; /* four extra characters for ' ..' */
   char buffer[GMS_MAX_PRINTLEN];

   int i;

   assert( scip != NULL );
   assert( strlen(rowname) > 0 || strlen(rownameextension) > 0 );
   assert( strcmp(type, "=e=") == 0 || strcmp(type, "=l=") == 0 || strcmp(type, "=g=") == 0 );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   (void) SCIPsnprintf(buffer, GMS_MAX_NAMELEN + 3, "%s%s ..", rowname, rownameextension);
   SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN + 3, buffer) );

   appendLine(scip, file, linebuffer, &linecnt, consname);

   /* print nonlinear terms
    */
   for( i = 0; i < nexprtrees; ++i )
   {
      assert(exprtrees[i] != NULL);
      if( exprtreecoefs[i] != 0.0 )
      {
         (void) SCIPsnprintf(buffer, GMS_MAX_NAMELEN + 3, "%+g * (", exprtreecoefs[i]);
         appendLineWithIndent(scip, file, linebuffer, &linecnt, buffer);
         SCIP_CALL( printExpr(scip, file, linebuffer, &linecnt, nsmooth, transformed, SCIPexprtreeGetRoot(exprtrees[i]), SCIPexprtreeGetVars(exprtrees[i])) );
         appendLineWithIndent(scip, file, linebuffer, &linecnt, ")");
      }
   }

   /* print linear terms, do after nonlinear since nonlinear may not print sign in beginning */
   if( nlinvars > 0 )
   {
      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "+", " ", nlinvars, linvars, lincoeffs, transformed) );
   }

   /* print right hand side */
   if( linecnt == 0 )
      /* we start a new line; therefore we tab this line */
      appendLine(scip, file, linebuffer, &linecnt, "     ");

   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s %.15g;", type, rhs);

   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/* print nonlinear row in GAMS format to file stream (performing retransformation to active linear variables)
 * */
static
SCIP_RETCODE printNonlinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< variables in linear part */
   SCIP_Real*            lincoeffs,          /**< coefficients of variables in linear part */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            exprtreecoefs,      /**< expression tree coefficients */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Bool*            nsmooth             /**< buffer to store whether we printed a nonsmooth function */
   )
{
   assert( scip != NULL );
   assert( strlen(rowname) > 0 );

   /* print row(s) in GAMS format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* print equality constraint */
      SCIP_CALL( printNonlinearRow(scip, file, rowname, "", "=e=",
         nlinvars, linvars, lincoeffs, nexprtrees, exprtrees, exprtreecoefs, rhs, transformed, nsmooth) );
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printNonlinearRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", "=g=",
            nlinvars, linvars, lincoeffs, nexprtrees, exprtrees, exprtreecoefs, lhs, transformed, nsmooth) );
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         SCIP_CALL( printNonlinearRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "=l=",
            nlinvars, linvars, lincoeffs, nexprtrees, exprtrees, exprtreecoefs, rhs, transformed, nsmooth) );
      }
   }

   return SCIP_OKAY;
}

/** method check if the variable names are not longer than GMS_MAX_NAMELEN */
static
SCIP_RETCODE checkVarnames(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   nvars               /**< number of variables */
   )
{
   int v;
   SCIP_VAR* var;
   SCIP_Bool replaceforbiddenchars;
   const char* badchar;

   assert( scip != NULL );
   assert( vars != NULL );

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/gmsreader/replaceforbiddenchars", &replaceforbiddenchars) );

   /* check if the variable names contain any of the bad symbols */
   for( badchar = badchars; *badchar; ++badchar )
   {
      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         assert( var != NULL );

         if( strchr(SCIPvarGetName(var), *badchar) != NULL )
         {
            if( replaceforbiddenchars )
            {
               SCIPinfoMessage(scip, NULL, "there is a variable name with symbol '%c', not allowed in GAMS format; all '%c' replaced by '_' (consider using 'write genproblem'/'write gentransproblem').\n", *badchar, *badchar);
            }
            else
            {
               SCIPwarningMessage(scip, "there is a variable name with symbol '%c', not allowed in GAMS format; use 'write genproblem'/'write gentransproblem', or set 'reading/gmsreader/replaceforbiddenchars' to TRUE and risk duplicate variable names.\n", *badchar);
            }

            break;
         }
      }
   }

   /* check if the variable names are too long */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );

      if( strlen(SCIPvarGetName(var)) > GMS_MAX_NAMELEN )
      {
         SCIPwarningMessage(scip, "there is a variable name which has to be cut down to %d characters; GAMS model might be corrupted.\n", 
            GMS_MAX_NAMELEN - 1);
         break;
      }
   }

   return SCIP_OKAY;
}

/** method check if the constraint names are not longer than GMS_MAX_NAMELEN */
static
SCIP_RETCODE checkConsnames(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   int c;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_Bool replaceforbiddenchars;
   const char* badchar;

   assert( scip != NULL );
   assert( conss != NULL );

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/gmsreader/replaceforbiddenchars", &replaceforbiddenchars) );

   /* check if the constraint names contain any of the bad symbols */
   for( badchar = badchars; *badchar; ++badchar )
   {
      for( c = 0; c < nconss; ++c )
      {
         cons = conss[c];
         assert( cons != NULL );

         if( strchr(SCIPconsGetName(cons), *badchar) != NULL )
         {
            if( replaceforbiddenchars )
            {
               SCIPinfoMessage(scip, NULL, "there is a constraint name with symbol '%c', not allowed in GAMS format; all '%c' replaced by '_' (consider using 'write genproblem'/'write gentransproblem').\n", *badchar, *badchar);
            }
            else
            {
               SCIPwarningMessage(scip, "there is a constraint name with symbol '%c', not allowed in GAMS format; use 'write genproblem'/'write gentransproblem', or set 'reading/gmsreader/replaceforbiddenchars' to TRUE and risk duplicate variable names.\n", *badchar);
            }

            break;
         }
      }
   }

   /* check if the constraint names are too long */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );

      /* in case the transformed is written, only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 || strcmp(conshdlrname, "quadratic") == 0 )
      {
         SCIP_Real lhs = strcmp(conshdlrname, "linear") == 0 ? SCIPgetLhsLinear(scip, cons) : SCIPgetLhsQuadratic(scip, cons);
         SCIP_Real rhs = strcmp(conshdlrname, "linear") == 0 ? SCIPgetLhsLinear(scip, cons) : SCIPgetRhsQuadratic(scip, cons);

         if( SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > GMS_MAX_NAMELEN )
         {
            SCIPwarningMessage(scip, "there is a constraint name which has to be cut down to %d characters;\n",
               GMS_MAX_NAMELEN - 1);
            break;
         }
         else if( !SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > GMS_MAX_NAMELEN - 4 )
         {
            SCIPwarningMessage(scip, "there is a constraint name which has to be cut down to %d characters;\n",
               GMS_MAX_NAMELEN - 5);
            break;
         }
      }
      else if( strlen(SCIPconsGetName(conss[c])) > GMS_MAX_NAMELEN )
      {
         SCIPwarningMessage(scip, "there is a constraint name which has to be cut down to %d characters;\n",
            GMS_MAX_NAMELEN - 1);
         break;
      }
   }
   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyGms)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderGms(scip) );

   return SCIP_OKAY;
}

#ifdef WITH_GAMS
/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadGms)
{
   SCIP_RETCODE ret;
   FILE* convertdopt;
   char gamscall[SCIP_MAXSTRLEN];
   char buffer[GMS_SSSIZE];
   int rc;
   gmoHandle_t gmo = NULL;
   gevHandle_t gev = NULL;

   assert(scip != NULL);
   assert(reader != NULL);
   assert(filename != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   ret = SCIP_ERROR;

   /* create temporary directory */
   mkdir("loadgms.tmp", S_IRWXU);

   /* create empty convertd options file */
   convertdopt = fopen("loadgms.tmp/convertd.opt", "w");
   if( convertdopt == NULL )
   {
      SCIPerrorMessage("Could not create convertd options file. Do you have write permissions in execution directory?\n");
      goto TERMINATE;
   }
   fputs(" ", convertdopt);
   fclose(convertdopt);

   /* call GAMS with convertd solver to get compiled model instance in temporary directory */
   SCIPsnprintf(gamscall, SCIP_MAXSTRLEN, WITH_GAMS "/gams %s LP=CONVERTD RMIP=CONVERTD QCP=CONVERTD RMIQCP=CONVERTD NLP=CONVERTD DNLP=CONVERTD RMINLP=CONVERTD CNS=CONVERTD MIP=CONVERTD MIQCP=CONVERTD MINLP=CONVERTD MCP=CONVERTD MPEC=CONVERTD RMPEC=CONVERTD SCRDIR=loadgms.tmp output=loadgms.tmp/listing optdir=loadgms.tmp optfile=1 pf4=0 solprint=0 limcol=0 limrow=0 pc=2 lo=%d",
      filename, SCIPgetVerbLevel(scip) == SCIP_VERBLEVEL_FULL ? 3 : 0);
   SCIPdebugMsg(scip, gamscall);
   rc = system(gamscall);
   if( rc != 0 )
   {
      SCIPerrorMessage("GAMS call returned with code %d, check loadgms.tmp/listing for details.\n", rc);
      /* likely the GAMS model could not be compiled, which we could report as a readerror */
      ret = SCIP_READERROR;
      goto TERMINATE;
   }

   /* initialize GEV library and create GEV */
   if( !gevCreateDD(&gev, WITH_GAMS, buffer, sizeof(buffer)) )
   {
      SCIPerrorMessage(buffer);
      goto TERMINATE;
   }

   /* initialize GMO library and create GMO */
   if( !gmoCreateDD(&gmo, WITH_GAMS, buffer, sizeof(buffer)) )
   {
      SCIPerrorMessage(buffer);
      goto TERMINATE;
   }

   /* load control file */
   if( gevInitEnvironmentLegacy(gev, "loadgms.tmp/gamscntr.dat") )
   {
      SCIPerrorMessage("Could not load control file loadgms.tmp/gamscntr.dat\n");
      goto TERMINATE;
   }

   /* tell GMO about GEV */
   if( gmoRegisterEnvironment(gmo, gev, buffer) )
   {
      SCIPerrorMessage("Error registering GAMS Environment: %s\n", buffer);
      goto TERMINATE;
   }

   /* load GAMS model instance into GMO */
   if( gmoLoadDataLegacy(gmo, buffer) )
   {
      SCIPerrorMessage("Could not load model data.\n");
      goto TERMINATE;
   }

   /* create SCIP problem out of GMO, using the magic from reader_gmo in interfaces/gams */
   SCIP_CALL( SCIPcreateProblemReaderGmo(scip, gmo, NULL, FALSE) );
   *result = SCIP_SUCCESS;

   ret = SCIP_OKAY;

TERMINATE:
   if( gmo != NULL )
      gmoFree(&gmo);
   if( gev != NULL )
      gevFree(&gev);

   /* remove temporary directory content (should have only files and directory itself) */
   if( ret != SCIP_READERROR )
      system("rm loadgms.tmp/* && rmdir loadgms.tmp");

   return ret;
}
#endif

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGms)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPwriteGms(scip, file, name, transformed, objsense, objscale, objoffset, vars,
         nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );

   return SCIP_OKAY;
}

#ifdef WITH_GAMS
/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeGms)
{
   if( gmoLibraryLoaded() )
      gmoLibraryUnload();
   if( gevLibraryLoaded() )
      gevLibraryUnload();

   return SCIP_OKAY;
}
#endif

/*
 * reader specific interface methods
 */

/** includes the gms file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderGms(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyGms) );
#ifdef WITH_GAMS
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadGms) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeGms) );
#endif
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteGms) );

   /* add gms reader parameters for writing routines*/
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/gmsreader/freeints", "have integer variables no upper bound by default (depending on GAMS version)?",
         NULL, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/gmsreader/replaceforbiddenchars", "shall characters '#', '*', '+', '/', and '-' in variable and constraint names be replaced by '_'?",
         NULL, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "reading/gmsreader/bigmdefault", "default M value for big-M reformulation of indicator constraints in case no bound on slack variable is given",
         NULL, FALSE, GMS_DEFAULT_BIGM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
         "reading/gmsreader/indicatorreform", "which reformulation to use for indicator constraints: 'b'ig-M, 's'os1",
         NULL, FALSE, GMS_DEFAULT_INDICATORREFORM, "bs", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/gmsreader/signpower", "is it allowed to use the gams function signpower(x,a)?",
         NULL, FALSE, GMS_DEFAULT_SIGNPOWER, NULL, NULL) );

   return SCIP_OKAY;
}


/** writes problem to gms file */
SCIP_RETCODE SCIPwriteGms(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**            vars,               /**< array with active variables ordered binary, integer, implicit, continuous */
   int                   nvars,              /**< number of active variables in the problem */
   int                   nbinvars,           /**< number of binary variables */
   int                   nintvars,           /**< number of general integer variables */
   int                   nimplvars,          /**< number of implicit integer variables */
   int                   ncontvars,          /**< number of continuous variables */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{
   int c;
   int v;
   int linecnt;
   char linebuffer[GMS_MAX_PRINTLEN];

   char varname[GMS_MAX_NAMELEN];
   char buffer[GMS_MAX_PRINTLEN];

   SCIP_Real* objcoeffs;

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;

   char consname[GMS_MAX_NAMELEN];

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_VAR* var;
   SCIP_VAR* objvar;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool freeints;
   SCIP_Bool nondefbounds;
   SCIP_Bool nlcons;
   SCIP_Bool nqcons;
   SCIP_Bool nsmooth;
   SCIP_Bool discrete;
   SCIP_Bool rangedrow;
   SCIP_Bool indicatorsosdef;
   SCIP_Bool signpowerallowed;
   SCIP_Bool needcomma;

   assert( scip != NULL );
   assert( vars != NULL || nvars == 0 );

   /* check if the variable names are not too long */
   SCIP_CALL( checkVarnames(scip, vars, nvars) );
   /* check if the constraint names are too long */
   SCIP_CALL( checkConsnames(scip, conss, nconss, transformed) );

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/gmsreader/signpower", &signpowerallowed) );

   /* check if the objective is a single continuous variable, so we would not have to introduce an auxiliary variable
    * for GAMS
    */
   objvar = NULL;
   if( objscale == 1.0 && objoffset == 0.0 )
   {
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPvarGetObj(vars[v]) == 0.0 ) /*lint !e613*/
            continue;

         if( objvar == NULL )
         {
            /* first variable with nonzero obj coefficient
             * if not active or having coefficient != 1.0, or being binary/integer, then give up
             */
            if( !SCIPvarIsActive(vars[v]) || SCIPvarGetObj(vars[v]) != 1.0 ||
               SCIPvarGetType(vars[v]) < SCIP_VARTYPE_IMPLINT ) /*lint !e613*/
               break;

            objvar = vars[v]; /*lint !e613*/
         }
         else
         {
            /* second variable with nonzero obj coefficient -> give up */
            objvar = NULL;
            break;
         }
      }
   }

   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "$OFFLISTING\n");
   SCIPinfoMessage(scip, file, "* SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "*   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "*   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "*   Constraints      : %d\n\n", nconss);

   /* print flags */
   SCIPinfoMessage(scip, file, "$MAXCOL %d\n", GMS_MAX_LINELEN - 1);
   SCIPinfoMessage(scip, file, "$OFFDIGIT\n\n");

   /* print variable section */
   SCIPinfoMessage(scip, file, "Variables\n");
   clearLine(linebuffer, &linecnt);

   if( objvar == NULL )
   {
      /* auxiliary objective variable */
      SCIPinfoMessage(scip, file, " objvar%c", nvars > 0 ? ',' : ';');
   }

   /* "model" variables */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v]; /*lint !e613*/
      assert( var != NULL );

      SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(var)) );
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%c", varname, (v < nvars - 1) ? ',' : ';');
      appendLine(scip, file, linebuffer, &linecnt, buffer);

      if( (linecnt > 0 && (v == nbinvars - 1 || v == nbinvars + nintvars - 1 ||
            v == nbinvars + nintvars + nimplvars - 1)) || v == nvars - 1 )
      {
         endLine(scip, file, linebuffer, &linecnt);
         clearLine(linebuffer, &linecnt);
      }
   }

   SCIPinfoMessage(scip, file, "\n");

   /* declare binary variables if present */
   if( nbinvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Binary variables\n");
      clearLine(linebuffer, &linecnt);

      for( v = 0; v < nbinvars; ++v )
      {
         var = vars[v]; /*lint !e613*/

         SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(var)) );
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%s", varname, (v < nbinvars - 1) ? "," : ";");

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      endLine(scip, file, linebuffer, &linecnt);
      SCIPinfoMessage(scip, file, "\n");
   }

   /* declare integer variables if present */
   if( nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Integer variables\n");
      clearLine(linebuffer, &linecnt);

      for( v = 0; v < nintvars; ++v )
      {
         var = vars[nbinvars + v]; /*lint !e613*/

         SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(var)) );
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%s", varname, (v < nintvars - 1) ? "," : ";");

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
      endLine(scip, file, linebuffer, &linecnt);
      SCIPinfoMessage(scip, file, "\n");
   }

   /* print variable bounds */
   SCIPinfoMessage(scip, file, "* Variable bounds\n");
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/gmsreader/freeints", &freeints) );
   nondefbounds = FALSE;

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v]; /*lint !e613*/
      assert( var != NULL );

      SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(var)) );

      if( transformed )
      {
         /* in case the transformed is written only local bounds are posted which are valid in the current node */
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
      }
      else
      {
         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);
      }
      assert( lb <= ub );

      /* fixed */
      if( SCIPisEQ(scip, lb, ub) )
      {
         if( v < nintvars )
            SCIPinfoMessage(scip, file, " %s.fx = %g;\n", varname, SCIPfloor(scip, lb + 0.5));
         else
            SCIPinfoMessage(scip, file, " %s.fx = %.15g;\n", varname, lb);
         nondefbounds = TRUE;

         /* no need to write lower and upper bounds additionally */
         continue;
      }

      /* lower bound */
      if( v < nbinvars + nintvars )
      {
         /* default lower bound of binaries and integers is 0 (also in recent gams versions if pf4=0 is given) */
         if( !SCIPisZero(scip, lb) )
         {
            if( !SCIPisInfinity(scip, -lb) )
               SCIPinfoMessage(scip, file, " %s.lo = %g;\n", varname, SCIPceil(scip, lb));
            else if( freeints )
               SCIPinfoMessage(scip, file, " %s.lo = -inf;\n", varname); /* -inf is allowed when running gams with pf4=0, which we assume if freeints is TRUE */
            else
               SCIPinfoMessage(scip, file, " %s.lo = %g;\n", varname, -SCIPinfinity(scip)); /* sorry, -inf not allowed in gams file here */
            nondefbounds = TRUE;
         }
      }
      else if( v >= nbinvars + nintvars && !SCIPisInfinity(scip, -lb) )
      {
         /* continuous variables are free by default */
         SCIPinfoMessage(scip, file, " %s.lo = %.15g;\n", varname, lb);
         nondefbounds = TRUE;
      }

      /* upper bound */
      if( v < nbinvars )
      {
         if( !SCIPisFeasEQ(scip, ub, 1.0) )
         {
            SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPfeasFloor(scip, ub));
            nondefbounds = TRUE;
         }
      }
      else if( v < nbinvars + nintvars && !freeints )
      {
         /* freeints == FALSE: integer variables have upper bound 100 by default */
         if( !SCIPisFeasEQ(scip, ub, 100.0) )
         {
            if( !SCIPisInfinity(scip, ub) )
               SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPfeasFloor(scip, ub));
            else
               SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPinfinity(scip)); /* sorry, +inf not allowed in gams file here (unless pf4=0) */
            nondefbounds = TRUE;
         }
      }
      else if( v < nbinvars + nintvars && !SCIPisInfinity(scip, ub) )
      {
         /* freeints == TRUE: integer variables have no upper bound by default */
         SCIPinfoMessage(scip, file, " %s.up = %g;\n", varname, SCIPfloor(scip, ub));
         nondefbounds = TRUE;
      }
      else if( v >= nbinvars + nintvars && !SCIPisInfinity(scip, ub) )
      {
         /* continuous variables are free by default */
         SCIPinfoMessage(scip, file, " %s.up = %.15g;\n", varname, ub);
         nondefbounds = TRUE;
      }
   }

   if( !nondefbounds )
      SCIPinfoMessage(scip, file, "* (All other bounds at default value: binary [0,1], integer [%s], continuous [-inf,+inf].)\n", freeints ? "0,+inf" : "0,100");
   SCIPinfoMessage(scip, file, "\n");

   /* print equations section */
   if( nconss > 0 || objvar == NULL )
   {
      SCIPinfoMessage(scip, file, "Equations\n");
      clearLine(linebuffer, &linecnt);
   }
   needcomma = FALSE;

   if( objvar == NULL )
   {
      SCIPinfoMessage(scip, file, " objequ");
      needcomma = TRUE;
   }

   /* declare equations */
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN, SCIPconsGetName(cons)) );
      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      rangedrow = strcmp(conshdlrname, "linear") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsLinear(scip, cons)) && !SCIPisInfinity(scip, SCIPgetRhsLinear(scip, cons))
         && !SCIPisEQ(scip, SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons));
      rangedrow = rangedrow || (strcmp(conshdlrname, "quadratic") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) && !SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons))
         && !SCIPisEQ(scip, SCIPgetLhsQuadratic(scip, cons), SCIPgetRhsQuadratic(scip, cons)));
      rangedrow = rangedrow || (strcmp(conshdlrname, "nonlinear") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsNonlinear(scip, cons)) && !SCIPisInfinity(scip, SCIPgetRhsNonlinear(scip, cons))
         && !SCIPisEQ(scip, SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons)));
      rangedrow = rangedrow || (strcmp(conshdlrname, "abspower") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsAbspower(scip, cons)) && !SCIPisInfinity(scip, SCIPgetRhsAbspower(scip, cons))
         && !SCIPisEQ(scip, SCIPgetLhsAbspower(scip, cons), SCIPgetRhsAbspower(scip, cons)));
      rangedrow = rangedrow || (strcmp(conshdlrname, "bivariate") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsBivariate(scip, cons)) && !SCIPisInfinity(scip, SCIPgetRhsBivariate(scip, cons))
         && !SCIPisEQ(scip, SCIPgetLhsBivariate(scip, cons), SCIPgetRhsBivariate(scip, cons)));
      rangedrow = rangedrow || (strcmp(conshdlrname, "varbound") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsVarbound(scip, cons)) && !SCIPisInfinity(scip, SCIPgetRhsVarbound(scip, cons))
         && !SCIPisEQ(scip, SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons)));

      /* we declare only those constraints which we can print in GAMS format */
      if( strcmp(conshdlrname, "knapsack") != 0 && strcmp(conshdlrname, "logicor") != 0 && strcmp(conshdlrname, "setppc") != 0
          && strcmp(conshdlrname, "linear") != 0 && strcmp(conshdlrname, "quadratic") != 0 && strcmp(conshdlrname, "varbound") != 0
          && strcmp(conshdlrname, "soc") != 0 && strcmp(conshdlrname, "abspower") != 0 && strcmp(conshdlrname, "bivariate") != 0
          && strcmp(conshdlrname, "nonlinear") != 0 && strcmp(conshdlrname, "SOS1") != 0 && strcmp(conshdlrname, "SOS2") != 0
          && strcmp(conshdlrname, "indicator") != 0 )
      {
         SCIPwarningMessage(scip, "Constraint type <%s> not supported. Skip writing constraint <%s>.\n", conshdlrname, SCIPconsGetName(cons));
         continue;
      }

      if( needcomma )
         appendLine(scip, file, linebuffer, &linecnt, ",");

      SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN, SCIPconsGetName(cons)) );
      if( rangedrow )
      {
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s%s%s%s", consname, "_lhs, ", consname, "_rhs");
         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
      else
      {
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s", consname);
         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
      needcomma = TRUE;
   }

   if( nconss > 0 || objvar == NULL )
   {
      (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, ";");
      appendLine(scip, file, linebuffer, &linecnt, buffer);

      endLine(scip, file, linebuffer, &linecnt);
      SCIPinfoMessage(scip, file, "\n");
   }

   if( objvar == NULL )
   {
      /* print objective function equation */
      clearLine(linebuffer, &linecnt);
      if( objoffset != 0.0 )
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " objequ .. objvar =e= %.15g + ", objscale * objoffset);
      else
         (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " objequ .. objvar =e= ");
      appendLine(scip, file, linebuffer, &linecnt, buffer);

      SCIP_CALL( SCIPallocBufferArray(scip, &objcoeffs, nvars) );

      for( v = 0; v < nvars; ++v )
      {
         var = vars[v]; /*lint !e613*/
         assert( var != NULL );

         /* in case the original problem has to be posted the variables have to be either "original" or "negated" */
         assert( transformed || SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED );

         objcoeffs[v] = SCIPisZero(scip, SCIPvarGetObj(var)) ? 0.0 : objscale * SCIPvarGetObj(var);
      }

      SCIP_CALL( printActiveVariables(scip, file, linebuffer, &linecnt, "", ";", nvars, vars, objcoeffs, transformed) );

      SCIPfreeBufferArray(scip, &objcoeffs);
      endLine(scip, file, linebuffer, &linecnt);
      SCIPinfoMessage(scip, file, "\n");
   }

   /* print constraints */
   nlcons = FALSE;
   nqcons = FALSE;
   nsmooth = FALSE;
   discrete = nbinvars > 0 || nintvars > 0;
   indicatorsosdef = FALSE;
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL );

      /* in case the transformed is written, only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      SCIP_CALL( printConformName(scip, consname, GMS_MAX_NAMELEN, SCIPconsGetName(cons)) );
      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
            consvals[v] = (SCIP_Real)weights[v];

         SCIP_CALL( printLinearCons(scip, file, consname, nconsvars, consvars, consvals,
               -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file, consname,
               SCIPgetNVarsLinear(scip, cons), SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
               SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed) );
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, file, consname,
               SCIPgetNVarsLogicor(scip, cons), SCIPgetVarsLogicor(scip, cons), NULL,
               1.0, SCIPinfinity(scip), transformed) );
      }
      else if( strcmp(conshdlrname, "quadratic") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
               SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
               SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
               SCIPgetLhsQuadratic(scip, cons),  SCIPgetRhsQuadratic(scip, cons), transformed) );

         nlcons = TRUE;
      }
      else if( strcmp(conshdlrname, "nonlinear") == 0 )
      {
         /* cons_nonlinear does not have exprtree's at hand during presolve */
         if( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_EXITPRESOLVE
             && SCIPgetExprgraphNonlinear(scip,conshdlr) != NULL )
         {
            SCIP_EXPRTREE* exprtree;
            SCIP_Real coef;

            SCIP_CALL( SCIPexprgraphGetTree(SCIPgetExprgraphNonlinear(scip,conshdlr), SCIPgetExprgraphNodeNonlinear(scip,cons), &exprtree) );
            coef = 1.0;
            SCIP_CALL( printNonlinearCons(scip, file, consname,
                  SCIPgetNLinearVarsNonlinear(scip, cons), SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
                  1, &exprtree, &coef,
                  SCIPgetLhsNonlinear(scip, cons),  SCIPgetRhsNonlinear(scip, cons), transformed, &nsmooth) );

            SCIP_CALL( SCIPexprtreeFree(&exprtree) );
         }
         else
         {
            SCIP_CALL( printNonlinearCons(scip, file, consname,
                  SCIPgetNLinearVarsNonlinear(scip, cons), SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
                  SCIPgetNExprtreesNonlinear(scip, cons), SCIPgetExprtreesNonlinear(scip, cons), SCIPgetExprtreeCoefsNonlinear(scip, cons),
                  SCIPgetLhsNonlinear(scip, cons),  SCIPgetRhsNonlinear(scip, cons), transformed, &nsmooth) );
         }
         nlcons = TRUE;
         nqcons = TRUE;
      }
      else if( strcmp(conshdlrname, "bivariate") == 0 )
      {
         SCIP_EXPRTREE* exprtree;
         SCIP_VAR* linvar;
         SCIP_Real lincoef;
         int exprdegree;
         SCIP_Real one;

         exprtree = SCIPgetExprtreeBivariate(scip, cons);
         assert(exprtree != NULL);

         linvar  = SCIPgetLinearVarBivariate(scip, cons);
         lincoef = SCIPgetLinearCoefBivariate(scip, cons);
         one = 1.0;
         SCIP_CALL( printNonlinearCons(scip, file, consname,
            linvar == NULL ? 0 : 1, &linvar, &lincoef,
            1, &exprtree, &one,
            SCIPgetLhsBivariate(scip, cons),  SCIPgetRhsBivariate(scip, cons), transformed, &nsmooth) );

         SCIP_CALL( SCIPexprtreeGetMaxDegree(exprtree, &exprdegree) );
         if( exprdegree > 1 )
            nlcons = TRUE;
         if( exprdegree > 2)
            nqcons = TRUE;
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printLinearCons(scip, file, consname,
                  nconsvars, consvars, NULL, 1.0, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printLinearCons(scip, file, consname,
                  nconsvars, consvars, NULL, -SCIPinfinity(scip), 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printLinearCons(scip, file, consname,
                  nconsvars, consvars, NULL, 1.0, SCIPinfinity(scip), transformed) );
            break;
         }
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( printLinearCons(scip, file, consname,
               2, consvars, consvals,
               SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "soc") == 0 )
      {
         SCIP_CALL( printSOCCons(scip, file, consname,
            SCIPgetNLhsVarsSOC(scip, cons), SCIPgetLhsVarsSOC(scip, cons), SCIPgetLhsCoefsSOC(scip, cons), SCIPgetLhsOffsetsSOC(scip, cons), SCIPgetLhsConstantSOC(scip, cons),
            SCIPgetRhsVarSOC(scip, cons), SCIPgetRhsCoefSOC(scip, cons), SCIPgetRhsOffsetSOC(scip, cons), transformed) );

         nlcons = nlcons || !isGAMSprintableSOC(SCIPgetNLhsVarsSOC(scip, cons), SCIPgetLhsVarsSOC(scip, cons), SCIPgetLhsCoefsSOC(scip, cons), SCIPgetLhsOffsetsSOC(scip, cons), SCIPgetLhsConstantSOC(scip, cons),
            SCIPgetRhsVarSOC(scip, cons), SCIPgetRhsCoefSOC(scip, cons), SCIPgetRhsOffsetSOC(scip, cons));
      }
      else if( strcmp(conshdlrname, "indicator") == 0 )
      {
         SCIP_CALL( printIndicatorCons(scip, file, consname,
            SCIPgetBinaryVarIndicator(cons), SCIPgetSlackVarIndicator(cons), &indicatorsosdef,
            transformed) );
      }
      else if( strcmp(conshdlrname, "abspower") == 0 )
      {
         SCIP_CALL( printSignpowerCons(scip, file, consname,
            SCIPgetNonlinearVarAbspower(scip, cons), SCIPgetLinearVarAbspower(scip, cons),
            SCIPgetExponentAbspower(scip, cons), SCIPgetOffsetAbspower(scip, cons), SCIPgetCoefLinearAbspower(scip, cons),
            SCIPgetLhsAbspower(scip, cons),  SCIPgetRhsAbspower(scip, cons), transformed, signpowerallowed, &nsmooth) );

         nlcons = TRUE;
         nqcons = TRUE;
      }
      else if( strcmp(conshdlrname, "SOS1") == 0 )
      {
         SCIP_CALL( printSOSCons(scip, file, consname,
            SCIPgetNVarsSOS1(scip, cons), SCIPgetVarsSOS1(scip, cons), 1,
            transformed) );
         discrete = TRUE;
      }
      else if( strcmp(conshdlrname, "SOS2") == 0 )
      {
         SCIP_CALL( printSOSCons(scip, file, consname,
            SCIPgetNVarsSOS2(scip, cons), SCIPgetVarsSOS2(scip, cons), 2,
            transformed) );
         discrete = TRUE;
      }
      else
      {
         SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
         SCIPinfoMessage(scip, file, "* ");
         SCIP_CALL( SCIPprintCons(scip, cons, file) );
         SCIPinfoMessage(scip, file, ";\n");
      }

      SCIPinfoMessage(scip, file, "\n");
   }
   /* if at most quadratic, then cannot have nonsmooth functions */
   assert(nlcons || !nsmooth);

   /* print model creation */
   SCIPinfoMessage(scip, file, "Model m / all /;\n\n");

   /* set some options to reduce listing file size */
   SCIPinfoMessage(scip, file, "option limrow = 0;\n");
   SCIPinfoMessage(scip, file, "option limcol = 0;\n\n");

   /* print solve command */
   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, "%s%s",
         discrete ? "MI" : "", nlcons ? (nqcons ? ((nsmooth && !discrete) ? "DNLP" : "NLP") : "QCP") : (discrete > 0 ? "P" : "LP"));

   if( objvar != NULL )
   {
      SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(objvar)) );
   }

   SCIPinfoMessage(scip, file, "$if not set %s $set %s %s\n", buffer, buffer, buffer);
   SCIPinfoMessage(scip, file, "Solve m using %%%s%% %simizing %s;\n",
         buffer, objsense == SCIP_OBJSENSE_MINIMIZE ? "min" : "max", objvar != NULL ? varname : "objvar");

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
