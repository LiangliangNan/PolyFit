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

/**@file   reader_gms.c
 * @ingroup DEFPLUGINS_READER
 * @brief  GAMS file writer
 * @author Ambros Gleixner
 * @author Stefan Vigerske
 *
 * @todo Check for words reserved for GAMS.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_indicator.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_varbound.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/pub_var.h"
#include "scip/reader_gms.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_reader.h"
#include "scip/scip_var.h"
#include "scip/expr_abs.h"
#include <string.h>


#define READER_NAME             "gmsreader"
#define READER_DESC             "file writer for (MI)(N)LPs in GAMS file format"
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
   SCIP_VAR***           vars,               /**< pointer to vars array to get active variables for */
   SCIP_Real**           scalars,            /**< pointer to scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   int*                  varssize,           /**< pointer to length of vars and scalars array */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( *vars != NULL );
   assert( scalars != NULL );
   assert( *scalars != NULL );
   assert( nvars != NULL );
   assert( varssize != NULL );
   assert( *varssize >= *nvars );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *varssize, constant, &requiredsize, TRUE) );

      if( requiredsize > *varssize )
      {
         *varssize = SCIPcalcMemGrowSize(scip, requiredsize);
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, *varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, scalars, *varssize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *varssize, constant, &requiredsize, TRUE) );
         assert(requiredsize <= *varssize);
      }
   }
   else
   {
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&(*vars)[v], &(*scalars)[v], constant) );
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

/** ends the given line with '\\0' and prints it to the given file stream, with a newline at the end */
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
   assert( 0 <= *linecnt && *linecnt < GMS_MAX_LINELEN );

   if( (*linecnt) > 0 )
   {
      linebuffer[(*linecnt)] = '\0';
      SCIPinfoMessage(scip, file, "%s\n", linebuffer);
      clearLine(linebuffer, linecnt);
   }
}

/** ends the given line with '\\0' and prints it to the given file stream, without a newline at the end */
static
void endLineNoNewline(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt             /**< number of characters in line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( 0 <= *linecnt && *linecnt < GMS_MAX_LINELEN );

   if( (*linecnt) > 0 )
   {
      linebuffer[(*linecnt)] = '\0';
      SCIPinfoMessage(scip, file, "%s", linebuffer);
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
   (void) strncat(linebuffer, extension, GMS_MAX_PRINTLEN - len);

   (*linecnt) += (int) strlen(extension);

   SCIPdebugMsg(scip, "linebuffer <%s>, length = %lu\n", linebuffer, (unsigned long)len);

   if( (*linecnt) > GMS_PRINTLEN )
      endLine(scip, file, linebuffer, linecnt);
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
   int activevarssize;
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
      activevarssize = nactivevars;

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, &activevars, &activevals, &nactivevars, &activevarssize, &activeconstant, transformed) );

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
   char linebuffer[GMS_MAX_PRINTLEN+1] = { '\0' };
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
      assert(vars != NULL);  /* for lint */
      assert(vals != NULL);

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
   int activevarssize;

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
      activevarssize = nactivevars;

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, &activevars, &activevals, &nactivevars, &activevarssize, &activeconstant, transformed) );
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
   char linebuffer[GMS_MAX_PRINTLEN+1] = { '\0' };
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
 * SOS1/2 Variable name_sosvar(name_sosset); name_sosvar.lo(name_sosset) = -inf;
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
   char linebuffer[GMS_MAX_PRINTLEN+1] = { '\0' };
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

/** prints expression in GAMS format to file stream */
static
SCIP_RETCODE printExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line buffer of length GMS_MAX_PRINTLEN */
   int*                  linecnt,            /**< number of characters in line so far */
   SCIP_Bool*            nsmooth,            /**< buffer to store whether we printed a nonsmooth function */
   SCIP_Bool*            nqcons,             /**< buffer to update whether we are still quadratic */
   SCIP_Bool             transformed,        /**< expression belongs to transformed constraint? */
   SCIP_EXPR*            expr                /**< expression to print */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPRITER_STAGE stage;
   int currentchild;
   unsigned int parentprecedence;
   long int fpos;
   SCIP_VAR** activevars = NULL;
   SCIP_Real* activecoefs = NULL;
   int nactivevars;
   int activevarssize;
   SCIP_Real activeconstant = 0.0;
   char varname[GMS_MAX_NAMELEN];
   unsigned int sumprecedence;

   assert(scip != NULL);
   assert(linebuffer != NULL);
   assert(linecnt != NULL);
   assert(nsmooth != NULL);
   assert(nqcons != NULL);
   assert(expr != NULL);

   if( file == NULL )
      file = stdout;

   appendLine(scip, file, linebuffer, linecnt, " ");

   /* store file position at begin of current line */
   fpos = ftell(file) - *linecnt;

   /* print out current buffer, as we print the expression directly to file */
   endLineNoNewline(scip, file, linebuffer, linecnt);

   activevarssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &activevars, activevarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activecoefs, activevarssize) );

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ALLSTAGES);

   sumprecedence = SCIPexprhdlrGetPrecedence(SCIPgetExprhdlrSum(scip));

   while( !SCIPexpriterIsEnd(it) )
   {
      stage = SCIPexpriterGetStageDFS(it);

      if( stage == SCIP_EXPRITER_VISITEDCHILD || stage == SCIP_EXPRITER_VISITINGCHILD )
         currentchild = SCIPexpriterGetChildIdxDFS(it);
      else
         currentchild = -1;

      if( SCIPexpriterGetParentDFS(it) != NULL )
         parentprecedence = SCIPexprhdlrGetPrecedence(SCIPexprGetHdlr(SCIPexpriterGetParentDFS(it)));
      else
         parentprecedence = 0;

      /* print a newline, if we have printed at least GMS_PRINTLEN chars since the last newline */
      if( ftell(file) > fpos + GMS_PRINTLEN )
      {
         SCIPinfoMessage(scip, file, "\n     ");
         /* store file position at begin of current line again; the -5 is for the whitespace we printed already */
         fpos = ftell(file) - 5;
      }

      if( SCIPisExprVar(scip, expr) )
      {
         /* special handler for variables:
          * - map to active variables
          * - pass variable name through printConformName
          */
         if( stage == SCIP_EXPRITER_ENTEREXPR )
         {
            activevars[0] = SCIPgetVarExprVar(expr);
            activecoefs[0] = 1.0;
            nactivevars = 1;

            SCIP_CALL( getActiveVariables(scip, &activevars, &activecoefs, &nactivevars, &activevarssize, &activeconstant, transformed) );

            if( nactivevars == 1 && activecoefs[0] == 1.0 && activeconstant == 0.0 )
            {
               SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(activevars[0])) );
               SCIPinfoMessage(scip, file, "%s", varname);
            }
            else
            {
               SCIP_Bool needsign = FALSE;
               int i;

               /* do as in print of expr_sum: an opening parenthesis may be required */
               if( sumprecedence <= parentprecedence )
                  SCIPinfoMessage(scip, file, "(");

               if( activeconstant != 0.0 )
               {
                  SCIPinfoMessage(scip, file, "%.15g", activeconstant);
                  needsign = TRUE;
               }
               for( i = 0; i < nactivevars; ++i )
               {
                  if( REALABS(activecoefs[i]) != 1.0 )
                  {
                     SCIPinfoMessage(scip, file, needsign ? "%+.15g*" : "%.15g*", activecoefs[i]);
                  }
                  else if( activecoefs[i] == 1.0 && needsign )
                  {
                     SCIPinfoMessage(scip, file, "+");
                  }
                  else if( activecoefs[i] == -1.0 )
                  {
                     SCIPinfoMessage(scip, file, "-");
                  }

                  SCIP_CALL( printConformName(scip, varname, GMS_MAX_NAMELEN, SCIPvarGetName(activevars[0])) );
                  SCIPinfoMessage(scip, file, "%s", varname);

                  needsign = TRUE;

                  /* check whether it is time for a linebreak */
                  if( ftell(file) > fpos + GMS_PRINTLEN )
                  {
                     SCIPinfoMessage(scip, file, "\n     ");
                     fpos = ftell(file) - 5;
                  }
               }

               if( sumprecedence <= parentprecedence )
                  SCIPinfoMessage(scip, file, ")");
            }
         }
      }
      else if( SCIPisExprPower(scip, expr) )
      {
         /* special handler for power */
         SCIP_Real exponent = SCIPgetExponentExprPow(expr);

         if( exponent == 2.0 )
         {
            /* write squares as "sqr(child)" */
            if( stage == SCIP_EXPRITER_ENTEREXPR )
               SCIPinfoMessage(scip, file, "sqr(");
            else if( stage == SCIP_EXPRITER_LEAVEEXPR )
               SCIPinfoMessage(scip, file, ")");
         }
         else if( EPSISINT(exponent, 0.0) )  /*lint !e835*/
         {
            /* write integer powers as "power(child, exponent)" */
            if( stage == SCIP_EXPRITER_ENTEREXPR )
               SCIPinfoMessage(scip, file, "power(");
            else if( stage == SCIP_EXPRITER_LEAVEEXPR )
               SCIPinfoMessage(scip, file, ",%g)", exponent);
            /* if power but not square, then we are no longer quadratic */
            *nqcons = FALSE;
         }
         else if( exponent == 0.5 )
         {
            /* write square roots as "sqrt(child)" */
            if( stage == SCIP_EXPRITER_ENTEREXPR )
               SCIPinfoMessage(scip, file, "sqrt(");
            else if( stage == SCIP_EXPRITER_LEAVEEXPR )
               SCIPinfoMessage(scip, file, ")");
            *nqcons = FALSE;
         }
         else
         {
            /* write any other power as "(child)**exponent" */
            if( stage == SCIP_EXPRITER_ENTEREXPR )
               SCIPinfoMessage(scip, file, "(");
            else if( stage == SCIP_EXPRITER_LEAVEEXPR )
               SCIPinfoMessage(scip, file, exponent >= 0.0 ? ")**%.15g" : ")**(%.15g)", exponent);
            *nqcons = FALSE;
         }
      }
      else
      {
         /* for any other expression, use the print callback of the exprhdlr */
         SCIP_CALL( SCIPcallExprPrint(scip, expr, stage, currentchild, parentprecedence, file) );

         if( !*nsmooth )
         {
            /* check for expression types that require changing modeltype from NLP to DNLP: currently only abs */
            if( SCIPisExprAbs(scip, expr) )
            {
               *nsmooth = TRUE;
               *nqcons = FALSE;
            }
         }

         /* if still quadratic, then check whether expression type is one that cannot occur in quadratics
          * allowed are sum, product, value, var, and power; the latter two were handled above
          */
         if( *nqcons && !SCIPisExprSum(scip, expr) && !SCIPisExprProduct(scip, expr) && !SCIPisExprValue(scip, expr) )
            *nqcons = FALSE;
      }

      expr = SCIPexpriterGetNext(it);
   }

   SCIPfreeExpriter(&it);

   SCIPfreeBufferArray(scip, &activecoefs);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}

/** print nonlinear row in GAMS format to file stream */
static
SCIP_RETCODE printNonlinearRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=e=", "=l=", or "=g=") */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Bool*            nsmooth,            /**< buffer to store whether we printed a nonsmooth function */
   SCIP_Bool*            nqcons              /**< buffer to update whether we are still quadratic */
   )
{
   char linebuffer[GMS_MAX_PRINTLEN+1] = { '\0' };
   int linecnt;
   char consname[GMS_MAX_NAMELEN + 3];
   char buffer[GMS_MAX_PRINTLEN];

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

   SCIP_CALL( printExpr(scip, file, linebuffer, &linecnt, nsmooth, nqcons, transformed, expr) );

   /* print right hand side */
   (void) SCIPsnprintf(buffer, GMS_MAX_PRINTLEN, " %s %.15g;", type, rhs);
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/** print nonlinear row in GAMS format to file stream (performing retransformation to active linear variables) */
static
SCIP_RETCODE printNonlinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Bool*            nsmooth,            /**< buffer to store whether we printed a nonsmooth function */
   SCIP_Bool*            nqcons              /**< buffer to update whether we are still quadratic */
   )
{
   assert( scip != NULL );
   assert( strlen(rowname) > 0 );

   /* print row(s) in GAMS format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* print equality constraint */
      SCIP_CALL( printNonlinearRow(scip, file, rowname, "", "=e=", expr, rhs, transformed, nsmooth, nqcons) );
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printNonlinearRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", "=g=", expr, lhs, transformed, nsmooth, nqcons) );
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         SCIP_CALL( printNonlinearRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "=l=", expr, rhs, transformed, nsmooth, nqcons) );
      }
   }

   if( *nqcons )
   {
      /* if we are still at most quadratic, check whether that is still the case when considering current constraint */
      SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, nqcons) );
      if( *nqcons )
         *nqcons = SCIPexprAreQuadraticExprsVariables(expr);
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

      if( strcmp(conshdlrname, "linear") == 0 || strcmp(conshdlrname, "nonlinear") == 0 )
      {
         SCIP_Real lhs = strcmp(conshdlrname, "linear") == 0 ? SCIPgetLhsLinear(scip, cons) : SCIPgetLhsNonlinear(cons);
         SCIP_Real rhs = strcmp(conshdlrname, "linear") == 0 ? SCIPgetLhsLinear(scip, cons) : SCIPgetRhsNonlinear(cons);

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


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGms)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPwriteGms(scip, file, name, transformed, objsense, objscale, objoffset, vars,
         nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );

   return SCIP_OKAY;
}

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
   char linebuffer[GMS_MAX_PRINTLEN+1];

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
   SCIP_Bool nlcons = FALSE;  /* whether there are nonlinear constraints */
   SCIP_Bool nqcons = TRUE;   /* whether nonlinear constraints are at most quadratic */
   SCIP_Bool nsmooth = FALSE; /* whether there are nonsmooth nonlinear constraints */
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
               SCIPinfoMessage(scip, file, " %s.up = +inf;\n", varname);
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
      rangedrow = rangedrow || (strcmp(conshdlrname, "nonlinear") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsNonlinear(cons)) && !SCIPisInfinity(scip, SCIPgetRhsNonlinear(cons))
         && !SCIPisEQ(scip, SCIPgetLhsNonlinear(cons), SCIPgetRhsNonlinear(cons)));
      rangedrow = rangedrow || (strcmp(conshdlrname, "varbound") == 0
         && !SCIPisInfinity(scip, -SCIPgetLhsVarbound(scip, cons)) && !SCIPisInfinity(scip, SCIPgetRhsVarbound(scip, cons))
         && !SCIPisEQ(scip, SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons)));

      /* we declare only those constraints which we can print in GAMS format */
      if( strcmp(conshdlrname, "knapsack") != 0 && strcmp(conshdlrname, "logicor") != 0 && strcmp(conshdlrname, "setppc") != 0
          && strcmp(conshdlrname, "linear") != 0 && strcmp(conshdlrname, "SOS1") != 0 && strcmp(conshdlrname, "SOS2") != 0
          && strcmp(conshdlrname, "nonlinear") != 0
          && strcmp(conshdlrname, "varbound") != 0
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
      else if( strcmp(conshdlrname, "nonlinear") == 0 )
      {
         SCIP_CALL( printNonlinearCons(scip, file, consname,
            SCIPgetExprNonlinear(cons), SCIPgetLhsNonlinear(cons), SCIPgetRhsNonlinear(cons), transformed, &nsmooth, &nqcons) );
         nlcons = TRUE;
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
      else if( strcmp(conshdlrname, "indicator") == 0 )
      {
         SCIP_CALL( printIndicatorCons(scip, file, consname,
            SCIPgetBinaryVarIndicator(cons), SCIPgetSlackVarIndicator(cons), &indicatorsosdef,
            transformed) );
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
         discrete ? "MI" : "", nlcons ? (nqcons ? "QCP" : ((nsmooth && !discrete) ? "DNLP" : "NLP")) : (discrete > 0 ? "P" : "LP"));

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
