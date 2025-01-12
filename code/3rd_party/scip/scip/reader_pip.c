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

/**@file   reader_pip.c
 * @ingroup DEFPLUGINS_READER
 * @brief  file reader for polynomial mixed-integer programs in PIP format
 * @author Stefan Vigerske
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>

#include "blockmemshell/memory.h"
#include "scip/reader_pip.h"
#include "scip/cons_and.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/expr_sum.h"
#include "scip/expr_var.h"
#include "scip/pub_cons.h"
#include "scip/pub_expr.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nlp.h"
#include "scip/pub_reader.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include "scip/scip_var.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <strings.h> /*lint --e{766}*/ /* needed for strncasecmp() */
#endif

#define READER_NAME             "pipreader"
#define READER_DESC             "file reader for polynomial mixed-integer programs in PIP format"
#define READER_EXTENSION        "pip"


/*
 * Data structures
 */
#define PIP_MAX_LINELEN        65536
#define PIP_MAX_PUSHEDTOKENS   2
#define PIP_INIT_MONOMIALSSIZE 128
#define PIP_INIT_FACTORSSIZE   16
#define PIP_MAX_PRINTLEN       561       /**< the maximum length of any line is 560 + '\\0' = 561*/
#define PIP_MAX_NAMELEN        256       /**< the maximum length for any name is 255 + '\\0' = 256 */
#define PIP_PRINTLEN           100

/** Section in PIP File */
enum PipSection
{
   PIP_START,
   PIP_OBJECTIVE,
   PIP_CONSTRAINTS,
   PIP_BOUNDS,
   PIP_GENERALS,
   PIP_BINARIES,
   PIP_END
};
typedef enum PipSection PIPSECTION;

enum PipExpType
{
   PIP_EXP_NONE,
   PIP_EXP_UNSIGNED,
   PIP_EXP_SIGNED
};
typedef enum PipExpType PIPEXPTYPE;

enum PipSense
{
   PIP_SENSE_NOTHING,
   PIP_SENSE_LE,
   PIP_SENSE_GE,
   PIP_SENSE_EQ
};
typedef enum PipSense PIPSENSE;

/** PIP reading data */
struct PipInput
{
   SCIP_FILE*            file;
   char                  linebuf[PIP_MAX_LINELEN+1];
   char                  probname[PIP_MAX_LINELEN];
   char                  objname[PIP_MAX_LINELEN];
   char*                 token;
   char*                 tokenbuf;
   char*                 pushedtokens[PIP_MAX_PUSHEDTOKENS];
   int                   npushedtokens;
   int                   linenumber;
   int                   linepos;
   PIPSECTION            section;
   SCIP_OBJSENSE         objsense;
   SCIP_Bool             initialconss;       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss;       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamiccols;        /**< should columns be added and removed dynamically to the LP? */
   SCIP_Bool             dynamicrows;        /**< should rows be added and removed dynamically to the LP? */
   SCIP_Bool             haserror;
};
typedef struct PipInput PIPINPUT;

static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = "-+:<>=*^";
static const char commentchars[] = "\\";
static const char namechars[] = "!#$%&;?@_";  /* and characters and numbers */


/*
 * Local methods (for reading)
 */

/** issues an error message and marks the PIP data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   const char*           msg                 /**< error message */
   )
{
   char formatstr[256];

   assert(pipinput != NULL);

   SCIPerrorMessage("Syntax error in line %d: %s ('%s')\n", pipinput->linenumber, msg, pipinput->token);
   if( pipinput->linebuf[strlen(pipinput->linebuf)-1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", pipinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", pipinput->linebuf);
   }
   (void) SCIPsnprintf(formatstr, 256, "         %%%ds\n", pipinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, formatstr, "^");
   pipinput->section  = PIP_END;
   pipinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   return pipinput->haserror;
}

/** returns whether the given character is a token delimiter */
static
SCIP_Bool isDelimChar(
   char                  c                   /**< input character */
   )
{
   return (c == '\0') || (strchr(delimchars, c) != NULL);
}

/** returns whether the given character is a single token */
static
SCIP_Bool isTokenChar(
   char                  c                   /**< input character */
   )
{
   return (strchr(tokenchars, c) != NULL);
}

/** returns whether the current character is member of a value string */
static
SCIP_Bool isValueChar(
   char                  c,                  /**< input character */
   char                  nextc,              /**< next input character */
   SCIP_Bool             firstchar,          /**< is the given character the first char of the token? */
   SCIP_Bool*            hasdot,             /**< pointer to update the dot flag */
   PIPEXPTYPE*           exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit((unsigned char)c) )
      return TRUE;
   else if( (*exptype == PIP_EXP_NONE) && !(*hasdot) && (c == '.') && isdigit((unsigned char)nextc) )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == PIP_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = PIP_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit((unsigned char)nextc) )
      {
         *exptype = PIP_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == PIP_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = PIP_EXP_UNSIGNED;
      return TRUE;
   }

   return FALSE;
}

/** reads the next line from the input file into the line buffer; skips comments;
 *  returns whether a line could be read
 */
static
SCIP_Bool getNextLine(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   int i;

   assert(scip != NULL);  /* for lint */
   assert(pipinput != NULL);

   /* clear the line */
   BMSclearMemoryArray(pipinput->linebuf, PIP_MAX_LINELEN);

   /* read next line */
   pipinput->linepos = 0;
   pipinput->linebuf[PIP_MAX_LINELEN-2] = '\0';
   if( SCIPfgets(pipinput->linebuf, (int) sizeof(pipinput->linebuf), pipinput->file) == NULL )
      return FALSE;
   pipinput->linenumber++;
   if( pipinput->linebuf[PIP_MAX_LINELEN-2] != '\0' )
   {
      SCIPerrorMessage("Error: line %d exceeds %d characters\n", pipinput->linenumber, PIP_MAX_LINELEN-2);
      pipinput->haserror = TRUE;
      return FALSE;
   }
   pipinput->linebuf[PIP_MAX_LINELEN-1] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(pipinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
      }
   }

   return TRUE;
}

/** swaps the addresses of two pointers */
static
void swapPointers(
   char**                pointer1,           /**< first pointer */
   char**                pointer2            /**< second pointer */
   )
{
   char* tmp;

   tmp = *pointer1;
   *pointer1 = *pointer2;
   *pointer2 = tmp;
}

/** reads the next token from the input file into the token buffer; returns whether a token was read */
static
SCIP_Bool getNextToken(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   SCIP_Bool hasdot;
   PIPEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(pipinput != NULL);
   assert(pipinput->linepos < PIP_MAX_LINELEN);

   /* check the token stack */
   if( pipinput->npushedtokens > 0 )
   {
      swapPointers(&pipinput->token, &pipinput->pushedtokens[pipinput->npushedtokens-1]);
      pipinput->npushedtokens--;
      SCIPdebugMsg(scip, "(line %d) read token again: '%s'\n", pipinput->linenumber, pipinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = pipinput->linebuf;
   while( isDelimChar(buf[pipinput->linepos]) )
   {
      if( buf[pipinput->linepos] == '\0' )
      {
         if( !getNextLine(scip, pipinput) )
         {
            pipinput->section = PIP_END;
            SCIPdebugMsg(scip, "(line %d) end of file\n", pipinput->linenumber);
            return FALSE;
         }
         assert(pipinput->linepos == 0);
      }
      else
         pipinput->linepos++;
   }
   assert(pipinput->linepos < PIP_MAX_LINELEN);
   assert(!isDelimChar(buf[pipinput->linepos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = PIP_EXP_NONE;
   if( isValueChar(buf[pipinput->linepos], buf[pipinput->linepos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < PIP_MAX_LINELEN);
         assert(!isDelimChar(buf[pipinput->linepos]));
         pipinput->token[tokenlen] = buf[pipinput->linepos];
         tokenlen++;
         pipinput->linepos++;
      }
      while( isValueChar(buf[pipinput->linepos], buf[pipinput->linepos+1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < PIP_MAX_LINELEN);
         pipinput->token[tokenlen] = buf[pipinput->linepos];
         tokenlen++;
         pipinput->linepos++;
         if( tokenlen == 1 && isTokenChar(pipinput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[pipinput->linepos]) && !isTokenChar(buf[pipinput->linepos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       */
      if( tokenlen >= 1
         && (pipinput->token[tokenlen-1] == '<' || pipinput->token[tokenlen-1] == '>' || pipinput->token[tokenlen-1] == '=')
         && buf[pipinput->linepos] == '=' )
      {
         pipinput->linepos++;
      }
      else if( pipinput->token[tokenlen-1] == '=' && (buf[pipinput->linepos] == '<' || buf[pipinput->linepos] == '>') )
      {
         pipinput->token[tokenlen-1] = buf[pipinput->linepos];
         pipinput->linepos++;
      }
   }
   assert(tokenlen < PIP_MAX_LINELEN);
   pipinput->token[tokenlen] = '\0';

   SCIPdebugMsg(scip, "(line %d) read token: '%s'\n", pipinput->linenumber, pipinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);
   assert(pipinput->npushedtokens < PIP_MAX_PUSHEDTOKENS);

   swapPointers(&pipinput->pushedtokens[pipinput->npushedtokens], &pipinput->token);
   pipinput->npushedtokens++;
}

/** puts the buffered token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushBufferToken(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);
   assert(pipinput->npushedtokens < PIP_MAX_PUSHEDTOKENS);

   swapPointers(&pipinput->pushedtokens[pipinput->npushedtokens], &pipinput->tokenbuf);
   pipinput->npushedtokens++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   swapPointers(&pipinput->token, &pipinput->tokenbuf);
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isNewSection(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   SCIP_Bool iscolon;

   assert(pipinput != NULL);

   /* remember first token by swapping the token buffer */
   swapTokenBuffer(pipinput);

   /* look at next token: if this is a ':', the first token is a name and no section keyword */
   iscolon = FALSE;
   if( getNextToken(scip, pipinput) )
   {
      iscolon = (strcmp(pipinput->token, ":") == 0);
      pushToken(pipinput);
   }

   /* reinstall the previous token by swapping back the token buffer */
   swapTokenBuffer(pipinput);

   /* check for ':' */
   if( iscolon )
      return FALSE;

   if( strcasecmp(pipinput->token, "MINIMIZE") == 0
      || strcasecmp(pipinput->token, "MINIMUM") == 0
      || strcasecmp(pipinput->token, "MIN") == 0 )
   {
      SCIPdebugMsg(scip, "(line %d) new section: OBJECTIVE\n", pipinput->linenumber);
      pipinput->section = PIP_OBJECTIVE;
      pipinput->objsense = SCIP_OBJSENSE_MINIMIZE;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "MAXIMIZE") == 0
      || strcasecmp(pipinput->token, "MAXIMUM") == 0
      || strcasecmp(pipinput->token, "MAX") == 0 )
   {
      SCIPdebugMsg(scip, "(line %d) new section: OBJECTIVE\n", pipinput->linenumber);
      pipinput->section = PIP_OBJECTIVE;
      pipinput->objsense = SCIP_OBJSENSE_MAXIMIZE;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "SUBJECT") == 0 )
   {
      /* check if the next token is 'TO' */
      swapTokenBuffer(pipinput);
      if( getNextToken(scip, pipinput) )
      {
         if( strcasecmp(pipinput->token, "TO") == 0 )
         {
            SCIPdebugMsg(scip, "(line %d) new section: CONSTRAINTS\n", pipinput->linenumber);
            pipinput->section = PIP_CONSTRAINTS;
            return TRUE;
         }
         else
            pushToken(pipinput);
      }
      swapTokenBuffer(pipinput);
   }

   if( strcasecmp(pipinput->token, "SUCH") == 0 )
   {
      /* check if the next token is 'THAT' */
      swapTokenBuffer(pipinput);
      if( getNextToken(scip, pipinput) )
      {
         if( strcasecmp(pipinput->token, "THAT") == 0 )
         {
            SCIPdebugMsg(scip, "(line %d) new section: CONSTRAINTS\n", pipinput->linenumber);
            pipinput->section = PIP_CONSTRAINTS;
            return TRUE;
         }
         else
            pushToken(pipinput);
      }
      swapTokenBuffer(pipinput);
   }

   if( strcasecmp(pipinput->token, "st") == 0
      || strcasecmp(pipinput->token, "S.T.") == 0
      || strcasecmp(pipinput->token, "ST.") == 0 )
   {
      SCIPdebugMsg(scip, "(line %d) new section: CONSTRAINTS\n", pipinput->linenumber);
      pipinput->section = PIP_CONSTRAINTS;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "BOUNDS") == 0
      || strcasecmp(pipinput->token, "BOUND") == 0 )
   {
      SCIPdebugMsg(scip, "(line %d) new section: BOUNDS\n", pipinput->linenumber);
      pipinput->section = PIP_BOUNDS;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "GENERAL") == 0
      || strcasecmp(pipinput->token, "GENERALS") == 0
      || strcasecmp(pipinput->token, "GEN") == 0
      || strcasecmp(pipinput->token, "INTEGER") == 0
      || strcasecmp(pipinput->token, "INTEGERS") == 0
      || strcasecmp(pipinput->token, "INT") == 0 )
   {
      SCIPdebugMsg(scip, "(line %d) new section: GENERALS\n", pipinput->linenumber);
      pipinput->section = PIP_GENERALS;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "BINARY") == 0
      || strcasecmp(pipinput->token, "BINARIES") == 0
      || strcasecmp(pipinput->token, "BIN") == 0 )
   {
      SCIPdebugMsg(scip, "(line %d) new section: BINARIES\n", pipinput->linenumber);
      pipinput->section = PIP_BINARIES;
      return TRUE;
   }

   if( strcasecmp(pipinput->token, "END") == 0 )
   {
      SCIPdebugMsg(scip, "(line %d) new section: END\n", pipinput->linenumber);
      pipinput->section = PIP_END;
      return TRUE;
   }

   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   PIPINPUT*             pipinput,           /**< PIP reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(pipinput != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( pipinput->token[1] == '\0' )
   {
      if( *pipinput->token == '+' )
         return TRUE;
      else if( *pipinput->token == '-' )
      {
         *sign *= -1;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isValue(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(pipinput != NULL);
   assert(value != NULL);

   if( strcasecmp(pipinput->token, "INFINITY") == 0 || strcasecmp(pipinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(pipinput->token, &endptr);
      if( endptr != pipinput->token && *endptr == '\0' )
      {
         *value = val;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is an equation sense */
static
SCIP_Bool isSense(
   PIPINPUT*             pipinput,           /**< PIP reading data */
   PIPSENSE*             sense               /**< pointer to store the equation sense, or NULL */
   )
{
   assert(pipinput != NULL);

   if( strcmp(pipinput->token, "<") == 0 )
   {
      if( sense != NULL )
         *sense = PIP_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(pipinput->token, ">") == 0 )
   {
      if( sense != NULL )
         *sense = PIP_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(pipinput->token, "=") == 0 )
   {
      if( sense != NULL )
         *sense = PIP_SENSE_EQ;
      return TRUE;
   }

   return FALSE;
}

/** returns the variable with the given name, or creates a new variable if it does not exist */
static
SCIP_RETCODE getVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 name,               /**< name of the variable */
   SCIP_Bool             dynamiccols,        /**< should columns be added and removed dynamically to the LP? */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Bool*            created             /**< pointer to store whether a new variable was created, or NULL */
   )
{
   assert(name != NULL);
   assert(var != NULL);

   *var = SCIPfindVar(scip, name);
   if( *var == NULL )
   {
      SCIP_VAR* newvar;

      /* create new variable of the given name */
      SCIPdebugMsg(scip, "creating new variable: <%s>\n", name);
      SCIP_CALL( SCIPcreateVar(scip, &newvar, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
            !dynamiccols, dynamiccols, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, newvar) );
      *var = newvar;

      /* because the variable was added to the problem, it is captured by SCIP and we can safely release it right now
       * without making the returned *var invalid
       */
      SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

      if( created != NULL )
         *created = TRUE;
   }
   else if( created != NULL )
      *created = FALSE;

   return SCIP_OKAY;
}

/** reads the header of the file */
static
SCIP_RETCODE readStart(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   /* everything before first section is treated as comment */
   do
   {
      /* get token */
      if( !getNextToken(scip, pipinput) )
         return SCIP_OKAY;
   }
   while( !isNewSection(scip, pipinput) );

   return SCIP_OKAY;
}

/** ensure that an array of monomials can hold a minimum number of entries */
static
SCIP_RETCODE ensureMonomialsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR***          monomials,          /**< pointer to current array of monomials */
   SCIP_Real**           monomialscoef,      /**< pointer to current array of monomial coefficients */
   int*                  monomialssize,      /**< current size of monomials array at input; new size at exit */
   int                   minnmonomials       /**< required minimal size of monomials array */
   )
{
   int newsize;

   assert(scip != NULL);
   assert(monomials != NULL);
   assert(monomialscoef != NULL);
   assert(monomialssize != NULL);
   assert(*monomials != NULL || *monomialssize == 0);

   if( minnmonomials <= *monomialssize )
      return SCIP_OKAY;

   newsize = SCIPcalcMemGrowSize(scip, minnmonomials);

   if( *monomials != NULL )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, monomials, newsize) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, monomials, newsize) );
   }
   if( *monomialscoef != NULL )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, monomialscoef, newsize) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, monomialscoef, newsize) );
   }
   *monomialssize = newsize;

   return SCIP_OKAY;
}

/** ensure that arrays of exponents and variable indices can hold a minimum number of entries */
static
SCIP_RETCODE ensureFactorsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to current array of variables */
   SCIP_Real**           exponents,          /**< pointer to current array of exponents */
   int*                  factorssize,        /**< current size of arrays at input; new size at exit */
   int                   minnfactors         /**< required minimal size of arrays */
   )
{
   int newsize;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(exponents != NULL);
   assert(factorssize != NULL);
   assert(*exponents != NULL || *factorssize == 0);
   assert(*vars != NULL || *factorssize == 0);

   if( minnfactors <= *factorssize )
      return SCIP_OKAY;

   newsize = SCIPcalcMemGrowSize(scip, minnfactors);

   if( *exponents != NULL )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, exponents, newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, vars, newsize) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, exponents, newsize) );
      SCIP_CALL( SCIPallocBufferArray(scip, vars, newsize) );
   }
   *factorssize = newsize;

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   PIP_MAX_LINELEN */
   SCIP_EXPR**           expr,               /**< pointer to store the constraint function as expression */
   SCIP_Bool*            islinear,           /**< pointer to store polynomial is linear */
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Real coef;
   int coefsign;
   int nextcoefsign;
   int monomialdegree;
   int i;

   SCIP_VAR** vars;
   SCIP_Real constant;

   SCIP_EXPR** monomials;
   SCIP_Real* monomialscoef;
   int monomialssize;
   int nmonomials;

   int nfactors;
   int factorssize;
   SCIP_Real* exponents;

   assert(scip != NULL);
   assert(pipinput != NULL);
   assert(name != NULL);
   assert(expr != NULL);
   assert(islinear != NULL);
   assert(newsection != NULL);

   *name = '\0';
   *expr = NULL;
   *islinear = TRUE;
   *newsection = FALSE;

   /* read the first token, which may be the name of the line */
   if( getNextToken(scip, pipinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(scip, pipinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* remember the token in the token buffer */
      swapTokenBuffer(pipinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(scip, pipinput) )
      {
         if( strcmp(pipinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the line name */
            (void)SCIPstrncpy(name, pipinput->tokenbuf, PIP_MAX_LINELEN);
            SCIPdebugMsg(scip, "(line %d) read constraint name: '%s'\n", pipinput->linenumber, name);
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse them as coefficients */
            pushToken(pipinput);
            pushBufferToken(pipinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it as coefficient */
         pushBufferToken(pipinput);
      }
   }

   /* initialize buffer for storing the monomials */
   monomialssize = PIP_INIT_MONOMIALSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &monomials, monomialssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &monomialscoef, monomialssize) );

   /* initialize buffer for storing the factors in a monomial */
   factorssize = PIP_INIT_FACTORSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &exponents, factorssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, factorssize) );

   /* read the coefficients */
   coefsign = +1;
   nextcoefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   nmonomials = 0;
   nfactors = 0;
   monomialdegree = 0;
   constant = 0.0;
   while( getNextToken(scip, pipinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool issense;
      SCIP_Bool issign;
      SCIP_Bool isnewsection;
      SCIP_Real exponent;

      issign = FALSE;   /* fix compiler warning */
      issense = FALSE;  /* fix lint warning */
      if( (isnewsection = isNewSection(scip, pipinput)) ||  /*lint !e820*/
         (issense = isSense(pipinput, NULL))      ||  /*lint !e820*/
         ((nfactors > 0 || havevalue) && (issign = isSign(pipinput, &nextcoefsign))) )  /*lint !e820*/
      {
         /* finish the current monomial */
         if( nfactors > 0 )
         {
            if( coefsign * coef != 0.0 )
            {
               SCIP_CALL( ensureMonomialsSize(scip, &monomials, &monomialscoef, &monomialssize, nmonomials + 1) );
               SCIP_CALL( SCIPcreateExprMonomial(scip, &monomials[nmonomials], nfactors, vars, exponents, NULL, NULL) );
               monomialscoef[nmonomials] = coefsign * coef;
               ++nmonomials;
            }
         }
         else if( havevalue )
         {
            constant += coefsign * coef;
         }

         if( monomialdegree > 1 )
            *islinear = FALSE;

         /* reset variables */
         nfactors = 0;
         coef = 1.0;
         coefsign = +1;
         havesign = FALSE;
         havevalue = FALSE;
         monomialdegree = 0;

         if( isnewsection )
         {
            *newsection = TRUE;
            break;
         }

         if( issense )
         {
            /* put the sense back onto the token stack */
            pushToken(pipinput);
            break;
         }

         if( issign )
         {
            coefsign = nextcoefsign;
            SCIPdebugMsg(scip, "(line %d) read coefficient sign: %+d\n", pipinput->linenumber, coefsign);
            havesign = TRUE;
            nextcoefsign = +1;
            continue;
         }
      }

      /* check if we read a sign */
      if( isSign(pipinput, &coefsign) )
      {
         SCIPdebugMsg(scip, "(line %d) read coefficient sign: %+d\n", pipinput->linenumber, coefsign);

         if( nfactors > 0 || havevalue )
         {
            syntaxError(scip, pipinput, "sign can only be at beginning of monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }

         havesign = TRUE;
         continue;
      }

      /* check if we are in between factors of a monomial */
      if( strcmp(pipinput->token, "*") == 0 )
      {
         if( nfactors == 0 )
         {
            syntaxError(scip, pipinput, "cannot have '*' before first variable in monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }

         continue;
      }

      /* all but the first monomial need a sign */
      if( nmonomials > 0 && !havesign )
      {
         syntaxError(scip, pipinput, "expected sign ('+' or '-') or sense ('<' or '>')");
         goto TERMINATE_READPOLYNOMIAL;
      }

      /* check if we are at an exponent for the last variable */
      if( strcmp(pipinput->token, "^") == 0 )
      {
         if( !getNextToken(scip, pipinput) || !isValue(scip, pipinput, &exponent) )
         {
            syntaxError(scip, pipinput, "expected exponent value after '^'");
            goto TERMINATE_READPOLYNOMIAL;
         }
         if( nfactors == 0 )
         {
            syntaxError(scip, pipinput, "cannot have '^' before first variable in monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }
         exponents[nfactors-1] = exponent;  /*lint !e530*/
         if( SCIPisIntegral(scip, exponent) && exponent > 0.0 ) /*lint !e530*/
            monomialdegree += (int)exponent - 1; /*lint !e530*//* -1, because we added +1 when we put the variable into varidxs */
         else
            *islinear = FALSE;

         SCIPdebugMsg(scip, "(line %d) read exponent value %g for variable %s\n", pipinput->linenumber, exponent,
            SCIPvarGetName(vars[nfactors-1]));
         continue;
      }

      /* check if we read a value */
      if( isValue(scip, pipinput, &coef) )
      {
         SCIPdebugMsg(scip, "(line %d) read coefficient value: %g with sign %+d\n", pipinput->linenumber, coef, coefsign);

         if( havevalue )
         {
            syntaxError(scip, pipinput, "two consecutive values");
            goto TERMINATE_READPOLYNOMIAL;
         }

         if( nfactors > 0 )
         {
            syntaxError(scip, pipinput, "coefficients can only be at the beginning of a monomial");
            goto TERMINATE_READPOLYNOMIAL;
         }

         havevalue = TRUE;
         continue;
      }

      /* the token is a variable name: get the corresponding variable (or create a new one) */
      SCIP_CALL( getVariable(scip, pipinput->token, pipinput->dynamiccols, &var, NULL) );

      /* ensure that there is enough memory to store all factors */
      SCIP_CALL( ensureFactorsSize(scip, &vars, &exponents, &factorssize, nfactors + 1) );

      /* create and store corresponding variable expression */
      vars[nfactors] = var;
      exponents[nfactors] = 1.0;
      ++nfactors;
      ++monomialdegree;
   }

   if( nfactors > 0 )
   {
      syntaxError(scip, pipinput, "string ended before monomial has finished");
      goto TERMINATE_READPOLYNOMIAL;
   }

   /* create sum expression consisting of all monomial expressions */
   SCIP_CALL( SCIPcreateExprSum(scip, expr, nmonomials, monomials, monomialscoef, constant, NULL, NULL) );

   /* release monomial expressions */
   for( i = 0; i < nmonomials; ++i )
   {
      assert(monomials[i] != NULL);
      SCIP_CALL( SCIPreleaseExpr(scip, &monomials[i]) );
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "read polynomial: ");
   SCIP_CALL( SCIPprintExpr(scip, *expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

 TERMINATE_READPOLYNOMIAL:
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &exponents);
   SCIPfreeBufferArray(scip, &monomialscoef);
   SCIPfreeBufferArray(scip, &monomials);

   return SCIP_OKAY;
}

/** reads the objective section */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   char name[PIP_MAX_LINELEN];
   SCIP_EXPR* expr;
   SCIP_Bool linear;
   SCIP_Bool newsection;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool dynamic;
   SCIP_Bool removable;

   assert(pipinput != NULL);

   /* determine settings; note that reading/{initialconss,dynamicconss,dynamicrows,dynamiccols} apply only to model
    * constraints and variables, not to an auxiliary objective constraint (otherwise it can happen that an auxiliary
    * objective variable is loose with infinite best bound, triggering the problem that an LP that is unbounded because
    * of loose variables with infinite best bound cannot be solved)
    */
   initial = TRUE;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = FALSE;
   removable = FALSE;

   /* read the objective coefficients */
   SCIP_CALL( readPolynomial(scip, pipinput, name, &expr, &linear, &newsection) );
   if( !hasError(pipinput) && expr != NULL )
   {
      SCIP_Real constant = SCIPgetConstantExprSum(expr);

      /* always create a variable that represents the constant; otherwise, this might lead to numerical issues on
       * instances with a relatively large constant, e.g., popdynm* instances
       */
      if( constant != 0.0 )
      {
         SCIP_VAR* objconst;
         SCIP_CALL( SCIPcreateVarBasic(scip, &objconst, "objconst", 1.0, 1.0, constant, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, objconst) );
         SCIP_CALL( SCIPreleaseVar(scip, &objconst) );

         /* remove the constant of the sum expression */
         SCIPsetConstantExprSum(expr, 0.0);
      }

      if( linear )
      {
         int i;

         /* set objective coefficients of variables */
         for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
         {
            SCIP_EXPR* child;
            SCIP_VAR* var;
            SCIP_Real coef;

            child = SCIPexprGetChildren(expr)[i];
            assert(child != NULL);
            assert(SCIPisExprVar(scip, child));

            /* child has to be a variable expression, see SCIPcreateExprMonomial() */
            var = SCIPgetVarExprVar(child);
            assert(var != NULL);
            coef = SCIPgetCoefsExprSum(expr)[i];

            /* adjust the objective coefficient */
            SCIP_CALL( SCIPchgVarObj(scip, var, SCIPvarGetObj(var) + coef) );
         }
      }
      else /* insert dummy variable and constraint to represent the nonlinear objective */
      {
         SCIP_EXPR* nonlinobjvarexpr;
         SCIP_VAR*  nonlinobjvar;
         SCIP_CONS* nonlinobjcons;
         SCIP_Real lhs;
         SCIP_Real rhs;

         SCIP_CALL( SCIPcreateVar(scip, &nonlinobjvar, "nonlinobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, nonlinobjvar) );

         if ( pipinput->objsense == SCIP_OBJSENSE_MINIMIZE )
         {
            lhs = -SCIPinfinity(scip);
            rhs = 0.0;
         }
         else
         {
            lhs = 0.0;
            rhs = SCIPinfinity(scip);
         }

         /* add created objective variable */
         SCIP_CALL( SCIPcreateExprVar(scip, &nonlinobjvarexpr, nonlinobjvar, NULL, NULL) );
         SCIP_CALL( SCIPappendExprSumExpr(scip, expr, nonlinobjvarexpr, -1.0) );
         SCIP_CALL( SCIPreleaseExpr(scip, &nonlinobjvarexpr) );

         /* create nonlinear constraint */
         SCIP_CALL( SCIPcreateConsNonlinear(scip, &nonlinobjcons, "nonlinobj", expr, lhs, rhs, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

         SCIP_CALL( SCIPaddCons(scip, nonlinobjcons) );
         SCIPdebugMsg(scip, "(line %d) added constraint <%s> to represent nonlinear objective: ", pipinput->linenumber, SCIPconsGetName(nonlinobjcons));
         SCIPdebugPrintCons(scip, nonlinobjcons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &nonlinobjcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &nonlinobjvar) );
      }
   }

   /* release expression */
   if( expr != NULL )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }

   return SCIP_OKAY;
}

/** reads the constraints section */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   char name[PIP_MAX_LINELEN];
   SCIP_EXPR* expr;
   SCIP_CONS* cons = NULL;
   SCIP_Bool linear;

   PIPSENSE sense;
   SCIP_Real sidevalue;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool newsection;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   int sidesign;

   assert(pipinput != NULL);

   /* read polynomial */
   SCIP_CALL( readPolynomial(scip, pipinput, name, &expr, &linear, &newsection) );
   if ( hasError(pipinput) )
      return SCIP_READERROR;
   if ( newsection )
   {
      if ( expr != NULL )
         syntaxError(scip, pipinput, "expected constraint sense '<=', '=', or '>='");
      return SCIP_OKAY;
   }

   /* read the constraint sense */
   if ( !getNextToken(scip, pipinput) )
   {
      syntaxError(scip, pipinput, "expected constraint sense.");
      return SCIP_READERROR;
   }
   if ( !isSense(pipinput, &sense) )
   {
      syntaxError(scip, pipinput, "expected constraint sense '<=', '=', or '>='");
      return SCIP_READERROR;
   }

   /* read the right hand side */
   sidesign = +1;
   if ( !getNextToken(scip, pipinput) )
   {
      syntaxError(scip, pipinput, "missing right hand side");
      return SCIP_READERROR;
   }
   if ( isSign(pipinput, &sidesign) )
   {
      if( !getNextToken(scip, pipinput) )
      {
         syntaxError(scip, pipinput, "missing value of right hand side");
         return SCIP_READERROR;
      }
   }
   if ( !isValue(scip, pipinput, &sidevalue) )
   {
      syntaxError(scip, pipinput, "expected value as right hand side");
      return SCIP_READERROR;
   }
   sidevalue *= sidesign;

   /* determine settings */
   initial = pipinput->initialconss;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = pipinput->dynamicconss;
   removable = pipinput->dynamicrows;

   /* assign the left and right hand side, depending on the constraint sense */
   switch ( sense )  /*lint !e530*/
   {
   case PIP_SENSE_GE:
      lhs = sidevalue;
      rhs = SCIPinfinity(scip);
      break;
   case PIP_SENSE_LE:
      lhs = -SCIPinfinity(scip);
      rhs = sidevalue;
      break;
   case PIP_SENSE_EQ:
      lhs = sidevalue;
      rhs = sidevalue;
      break;
   case PIP_SENSE_NOTHING:
   default:
      SCIPerrorMessage("invalid constraint sense <%d>\n", sense);
      return SCIP_INVALIDDATA;
   }

   /* linear constraint function */
   if( linear )
   {
      SCIP_VAR** vars;
      SCIP_Real* coefs;
      SCIP_Real constant;
      int nchildren;
      int i;

      nchildren = SCIPexprGetNChildren(expr);
      constant = SCIPgetConstantExprSum(expr);
      coefs = SCIPgetCoefsExprSum(expr);

      /* allocate memory to store variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nchildren) );

      /* collect variables */
      for( i = 0; i < nchildren; ++i )
      {
         SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
         assert(child != NULL);
         assert(SCIPisExprVar(scip, child));

         vars[i] = SCIPgetVarExprVar(child);
         assert(vars[i] != NULL);
      }

      /* adjust lhs and rhs */
      if( !SCIPisInfinity(scip, -lhs) )
         lhs -= constant;
      if( !SCIPisInfinity(scip, rhs) )
         rhs -= constant;

      /* create linear constraint */
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nchildren, vars, coefs, lhs, rhs, initial, separate, enforce,
         check, propagate, local, modifiable, dynamic, removable, FALSE) );

      /* free memory */
      SCIPfreeBufferArray(scip, &vars);
   }
   else /* nonlinear constraint function */
   {
      SCIP_CALL( SCIPcreateConsNonlinear(scip, &cons, name, expr, lhs, rhs, initial, separate, enforce, check, propagate,
      local, modifiable, dynamic, removable) );
   }

   /* add and release constraint */
   assert(cons != NULL);
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugMsg(scip, "(line %d) created constraint: ", pipinput->linenumber);
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* release expression */
   if( expr != NULL )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }

   return SCIP_OKAY;
}

/** reads the bounds section */
static
SCIP_RETCODE readBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   while( getNextToken(scip, pipinput) )
   {
      SCIP_VAR* var;
      SCIP_Real value;
      SCIP_Real lb;
      SCIP_Real ub;
      int sign;
      SCIP_Bool hassign;
      PIPSENSE leftsense;

      /* check if we reached a new section */
      if( isNewSection(scip, pipinput) )
         return SCIP_OKAY;

      /* default bounds are [0,+inf] */
      lb = 0.0;
      ub = SCIPinfinity(scip);
      leftsense = PIP_SENSE_NOTHING;

      /* check if the first token is a sign */
      sign = +1;
      hassign = isSign(pipinput, &sign);
      if( hassign && !getNextToken(scip, pipinput) )
      {
         syntaxError(scip, pipinput, "expected value");
         return SCIP_OKAY;
      }

      /* the first token must be either a value or a variable name */
      if( isValue(scip, pipinput, &value) )
      {
         /* first token is a value: the second token must be a sense */
         if( !getNextToken(scip, pipinput) || !isSense(pipinput, &leftsense) )
         {
            syntaxError(scip, pipinput, "expected bound sense '<=', '=', or '>='");
            return SCIP_OKAY;
         }

         /* update the bound corresponding to the sense */
         switch( leftsense )
         {
         case PIP_SENSE_GE:
            ub = sign * value;
            break;
         case PIP_SENSE_LE:
            lb = sign * value;
            break;
         case PIP_SENSE_EQ:
            lb = sign * value;
            ub = sign * value;
            break;
         case PIP_SENSE_NOTHING:
         default:
            SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
            return SCIP_INVALIDDATA;
         }
      }
      else if( hassign )
      {
         syntaxError(scip, pipinput, "expected value");
         return SCIP_OKAY;
      }
      else
         pushToken(pipinput);

      /* the next token must be a variable name */
      if( !getNextToken(scip, pipinput) )
      {
         syntaxError(scip, pipinput, "expected variable name");
         return SCIP_OKAY;
      }
      SCIP_CALL( getVariable(scip, pipinput->token, pipinput->dynamiccols, &var, NULL) );

      /* the next token might be another sense, or the word "free" */
      if( getNextToken(scip, pipinput) )
      {
         PIPSENSE rightsense;

         if( isSense(pipinput, &rightsense) )
         {
            /* check, if the senses fit */
            if( leftsense == PIP_SENSE_NOTHING
               || (leftsense == PIP_SENSE_LE && rightsense == PIP_SENSE_LE)
               || (leftsense == PIP_SENSE_GE && rightsense == PIP_SENSE_GE) )
            {
               if( !getNextToken(scip, pipinput) )
               {
                  syntaxError(scip, pipinput, "expected value or sign");
                  return SCIP_OKAY;
               }

               /* check if the next token is a sign */
               sign = +1;
               hassign = isSign(pipinput, &sign);
               if( hassign && !getNextToken(scip, pipinput) )
               {
                  syntaxError(scip, pipinput, "expected value");
                  return SCIP_OKAY;
               }

               /* the next token must be a value */
               if( !isValue(scip, pipinput, &value) )
               {
                  syntaxError(scip, pipinput, "expected value");
                  return SCIP_OKAY;
               }

               /* update the bound corresponding to the sense */
               switch( rightsense )
               {
               case PIP_SENSE_GE:
                  lb = sign * value;
                  break;
               case PIP_SENSE_LE:
                  ub = sign * value;
                  break;
               case PIP_SENSE_EQ:
                  lb = sign * value;
                  ub = sign * value;
                  break;
               case PIP_SENSE_NOTHING:
               default:
                  SCIPerrorMessage("invalid bound sense <%d>\n", leftsense);
                  return SCIP_INVALIDDATA;
               }
            }
            else
            {
               syntaxError(scip, pipinput, "the two bound senses do not fit");
               return SCIP_OKAY;
            }
         }
         else if( strcasecmp(pipinput->token, "FREE") == 0 )
         {
            if( leftsense != PIP_SENSE_NOTHING )
            {
               syntaxError(scip, pipinput, "variable with bound is marked as 'free'");
               return SCIP_OKAY;
            }
            lb = -SCIPinfinity(scip);
            ub = SCIPinfinity(scip);
         }
         else
         {
            /* the token was no sense: push it back to the token stack */
            pushToken(pipinput);
         }
      }

      /* change the bounds of the variable if bounds have been given (do not destroy earlier specification of bounds) */
      if ( lb != 0.0 )
         SCIP_CALL( SCIPchgVarLb(scip, var, lb) );
      /*lint --e{777}*/
      if ( ub != SCIPinfinity(scip) )
         SCIP_CALL( SCIPchgVarUb(scip, var, ub) );
      SCIPdebugMsg(scip, "(line %d) new bounds: <%s>[%g,%g]\n", pipinput->linenumber, SCIPvarGetName(var),
         SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   return SCIP_OKAY;
}

/** reads the generals section */
static
SCIP_RETCODE readGenerals(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   while( getNextToken(scip, pipinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;
      SCIP_Bool infeasible;

      /* check if we reached a new section */
      if( isNewSection(scip, pipinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, pipinput->token, pipinput->dynamiccols, &var, &created) );
      if( created )
      {
         syntaxError(scip, pipinput, "unknown variable in generals section");
         return SCIP_OKAY;
      }

      /* mark the variable to be integral */
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER, &infeasible) );
      /* don't assert feasibility here because the presolver will and should detect a infeasibility */
   }

   return SCIP_OKAY;
}

/** reads the binaries section */
static
SCIP_RETCODE readBinaries(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   assert(pipinput != NULL);

   while( getNextToken(scip, pipinput) )
   {
      SCIP_VAR* var;
      SCIP_Bool created;
      SCIP_Bool infeasible;

      /* check if we reached a new section */
      if( isNewSection(scip, pipinput) )
         return SCIP_OKAY;

      /* the token must be the name of an existing variable */
      SCIP_CALL( getVariable(scip, pipinput->token, pipinput->dynamiccols, &var, &created) );
      if( created )
      {
         syntaxError(scip, pipinput, "unknown variable in binaries section");
         return SCIP_OKAY;
      }

      /* mark the variable to be binary and change its bounds appropriately */
      if( SCIPvarGetLbGlobal(var) < 0.0 )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 0.0) );
      }
      if( SCIPvarGetUbGlobal(var) > 1.0 )
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 1.0) );
      }
      SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
      /* don't assert feasibility here because the presolver will and should detect a infeasibility */
   }

   return SCIP_OKAY;
}

/** reads a PIP file
 */
static
SCIP_RETCODE readPIPFile(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   const char*           filename            /**< name of the input file */
   )
{
   assert(pipinput != NULL);

   /* open file */
   pipinput->file = SCIPfopen(filename, "r");
   if( pipinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* parse the file */
   pipinput->section = PIP_START;
   while( pipinput->section != PIP_END && !hasError(pipinput) )
   {
      switch( pipinput->section )
      {
      case PIP_START:
         SCIP_CALL( readStart(scip, pipinput) );
         break;

      case PIP_OBJECTIVE:
         SCIP_CALL( readObjective(scip, pipinput) );
         break;

      case PIP_CONSTRAINTS:
         SCIP_CALL( readConstraints(scip, pipinput) );
         break;

      case PIP_BOUNDS:
         SCIP_CALL( readBounds(scip, pipinput) );
         break;

      case PIP_GENERALS:
         SCIP_CALL( readGenerals(scip, pipinput) );
         break;

      case PIP_BINARIES:
         SCIP_CALL( readBinaries(scip, pipinput) );
         break;

      case PIP_END: /* this is already handled in the while() loop */
      default:
         SCIPerrorMessage("invalid PIP file section <%d>\n", pipinput->section);
         return SCIP_INVALIDDATA;
      }
   }

   /* close file */
   SCIPfclose(pipinput->file);

   return SCIP_OKAY;
}


/*
 * Local methods (for writing)
 */

/** hash key retrieval function for variables */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff the indices of both variables are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVar)
{  /*lint --e{715}*/
   if ( key1 == key2 )
      return TRUE;
   return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVar)
{  /*lint --e{715}*/
   assert( SCIPvarGetIndex((SCIP_VAR*) key) >= 0 );
   return (unsigned int) SCIPvarGetIndex((SCIP_VAR*) key);
}

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
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
   }
   return SCIP_OKAY;
}

/** checks whether a given expression is a signomial
 *
 * assumes simplified expression
 */
static
SCIP_Bool isExprSignomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);

   if( SCIPisExprVar(scip, expr) || SCIPisExprValue(scip, expr) )
      return TRUE;

   if( SCIPisExprPower(scip, expr) && SCIPisExprVar(scip, SCIPexprGetChildren(expr)[0]) )
      return TRUE;

   if( SCIPisExprProduct(scip, expr) )
   {
      SCIP_EXPR* child;
      int c;

      for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      {
         child = SCIPexprGetChildren(expr)[c];

         if( SCIPisExprVar(scip, child) )
            continue;

         if( SCIPisExprPower(scip, child) && SCIPisExprVar(scip, SCIPexprGetChildren(child)[0]) )
            continue;

         /* the pip format does not allow constants here */

         return FALSE;
      }

      return TRUE;
   }

   return FALSE;
}

/** checks whether a given expression is a sum of signomials (i.e., like a polynomial, but negative and fractional exponents allowed)
 *
 * assumes simplified expression;
 * does not check whether variables in powers with fractional exponent are nonnegative;
 * does not check whether variables in powers with negative exponent are bounded away from zero (the format specification does not require that, too)
 */
static
SCIP_Bool isExprPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int c;

   assert(scip != NULL);
   assert(expr != NULL);

   if( !SCIPisExprSum(scip, expr) )
      return isExprSignomial(scip, expr);

   /* check whether every term of sum is signomial */
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      if( !isExprSignomial(scip, SCIPexprGetChildren(expr)[c]) )
         return FALSE;

   return TRUE;
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
   assert( 0 <= *linecnt && *linecnt < PIP_MAX_PRINTLEN );

   if( (*linecnt) > 0 )
   {
      linebuffer[(*linecnt)] = '\0';
      SCIPinfoMessage(scip, file, "%s\n", linebuffer);
      clearLine(linebuffer, linecnt);
   }
}

/** appends extension to line and prints it to the give file stream if the
 *  line exceeded the length given in the define PIP_PRINTLEN */
static
void appendLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line */
   int*                  linecnt,            /**< number of characters in line */
   const char*           extension           /**< string to extent the line */
   )
{
   assert( scip != NULL );
   assert( linebuffer != NULL );
   assert( linecnt != NULL );
   assert( extension != NULL );
   assert( strlen(linebuffer) + strlen(extension) < PIP_MAX_PRINTLEN );

   /* NOTE: avoid
    *   sprintf(linebuffer, "%s%s", linebuffer, extension);
    * because of overlapping memory areas in memcpy used in sprintf.
    */
   (void) strncat(linebuffer, extension, PIP_MAX_PRINTLEN - strlen(linebuffer));

   (*linecnt) += (int) strlen(extension);

   SCIPdebugMsg(scip, "linebuffer <%s>, length = %lu\n", linebuffer, (unsigned long)strlen(linebuffer));

   if( (*linecnt) > PIP_PRINTLEN )
      endLine(scip, file, linebuffer, linecnt);
}


/** print linear or quadratic row in PIP format to file stream */
static
SCIP_RETCODE printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=", "<=", or ">=") */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficient values */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_EXPR*            quadexpr,           /**< quadratic expression */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed         /**< transformed constraint? */
   )
{
   int v;
   char linebuffer[PIP_MAX_PRINTLEN+1] = { '\0' };
   int linecnt;

   char varname[PIP_MAX_NAMELEN];
   char varname2[PIP_MAX_NAMELEN];
   char consname[PIP_MAX_NAMELEN + 1]; /* an extra character for ':' */
   char buffer[PIP_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strcmp(type, "=") == 0 || strcmp(type, "<=") == 0 || strcmp(type, ">=") == 0 );
   assert( nlinvars == 0 || (linvars != NULL && linvals != NULL) );

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   if ( strlen(rowname) > 0 || strlen(rownameextension) > 0 )
   {
      (void) SCIPsnprintf(consname, PIP_MAX_NAMELEN + 1, "%s%s:", rowname, rownameextension);
      appendLine(scip, file, linebuffer, &linecnt, consname);
   }

   /* print coefficients */
   for( v = 0; v < nlinvars; ++v )
   {
      SCIP_VAR* var;

      assert(linvars != NULL);  /* for lint */
      assert(linvals != NULL);

      var = linvars[v];
      assert( var != NULL );

      /* we start a new line; therefore we tab this line */
      if ( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");

      (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", linvals[v], varname);

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print quadratic part */
   if( quadexpr != NULL )
   {
      SCIP_EXPR** linexprs;
      SCIP_VAR** activevars;
      SCIP_Real* activevals;
      SCIP_Real* lincoefs;
      SCIP_Real constant;
      SCIP_Real activeconstant = 0.0;
      int nbilinexprterms;
      int nactivevars;
      int nquadexprs;
      int nlinexprs;

      /* get data from the quadratic expression */
      SCIPexprGetQuadraticData(quadexpr, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, &nbilinexprterms,
            NULL, NULL);

      /* allocate memory to store active linear variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &activevars, nlinexprs) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, lincoefs, nlinexprs) );
      nactivevars = nlinexprs;

      for( v = 0; v < nlinexprs; ++v )
      {
         assert(linexprs != NULL && linexprs[v] != NULL);
         assert(SCIPisExprVar(scip, linexprs[v]));

         activevars[v] = SCIPgetVarExprVar(linexprs[v]);
         assert(activevars[v] != NULL);
      }

      /* get active variables */
      SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );
      constant += activeconstant;

      /* print linear coefficients of linear variables */
      for( v = 0; v < nactivevars; ++v )
      {
         SCIP_VAR* var;

         assert(activevars != NULL);  /* for lint */
         assert(activevals != NULL);

         var = activevars[v];
         assert( var != NULL );

         /* we start a new line; therefore we tab this line */
         if( linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", activevals[v], varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      /* free memory for active linear variables */
      SCIPfreeBufferArray(scip, &activevals);
      SCIPfreeBufferArray(scip, &activevars);

      /* adjust rhs if there is a constant */
      if( constant != 0.0 && !SCIPisInfinity(scip, rhs) )
         rhs -= constant;

      /* print linear coefficients of quadratic variables */
      for( v = 0; v < nquadexprs; ++v )
      {
         SCIP_EXPR* expr;
         SCIP_VAR* var;
         SCIP_Real lincoef;

         /* get linear coefficient and variable of quadratic term */
         SCIPexprGetQuadraticQuadTerm(quadexpr, v, &expr, &lincoef, NULL, NULL, NULL, NULL);
         assert(expr != NULL);
         assert(SCIPisExprVar(scip, expr));

         var = SCIPgetVarExprVar(expr);
         assert(var != NULL);

         if( lincoef == 0.0 )
            continue;

         /* we start a new line; therefore we tab this line */
         if( linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", lincoef, varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      /* print square terms */
      for( v = 0; v < nquadexprs; ++v )
      {
         SCIP_EXPR* expr;
         SCIP_VAR* var;
         SCIP_Real sqrcoef;

         /* get square coefficient and variable of quadratic term */
         SCIPexprGetQuadraticQuadTerm(quadexpr, v, &expr, NULL, &sqrcoef, NULL, NULL, NULL);
         assert(expr != NULL);
         assert(SCIPisExprVar(scip, expr));

         var = SCIPgetVarExprVar(expr);
         assert(var != NULL);

         if( sqrcoef == 0.0 )
            continue;

         /* we start a new line; therefore we tab this line */
         if( linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s^2", sqrcoef, varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      /* print bilinear terms */
      for( v = 0; v < nbilinexprterms; ++v )
      {
         SCIP_EXPR* expr1;
         SCIP_EXPR* expr2;
         SCIP_VAR* var1;
         SCIP_VAR* var2;
         SCIP_Real bilincoef;

         /* get coefficient and variables of bilinear */
         SCIPexprGetQuadraticBilinTerm(quadexpr, v, &expr1, &expr2, &bilincoef, NULL, NULL);
         assert(expr1 != NULL);
         assert(SCIPisExprVar(scip, expr1));
         assert(expr2 != NULL);
         assert(SCIPisExprVar(scip, expr2));

         var1 = SCIPgetVarExprVar(expr1);
         assert(var1 != NULL);
         var2 = SCIPgetVarExprVar(expr2);
         assert(var2 != NULL);

         /* we start a new line; therefore we tab this line */
         if( linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname,  PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var1));
         (void) SCIPsnprintf(varname2, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var2));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s * %s", bilincoef, varname, varname2);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }
   }

   /* print right hand side */
   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %s %+.15g", type, rhs);

   /* we start a new line; therefore we tab this line */
   if (linecnt == 0 )
      appendLine(scip, file, linebuffer, &linecnt, " ");
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/** print signomial in PIP format to file stream */
static
void printSignomial(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line buffer to append to */
   int*                  linecnt,            /**< count on line buffer use */
   SCIP_EXPR*            expr,               /**< sigomial expression */
   SCIP_Real             coef,               /**< coefficient */
   SCIP_Bool             needsign            /**< whether a sign needs to be ensured */
   )
{
   char buffer[PIP_MAX_PRINTLEN];
   SCIP_EXPR* child;
   int c;

   assert(isExprSignomial(scip, expr));

   if( SCIPisExprProduct(scip, expr) )
      coef *= SCIPgetCoefExprProduct(expr);

   if( SCIPisExprValue(scip, expr) )
      coef *= SCIPgetValueExprValue(expr);

   if( REALABS(coef) != 1.0 )
   {
      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, needsign ? " %+.15g " : " %.15g ", coef);
      appendLine(scip, file, linebuffer, linecnt, buffer);
   }
   else if( coef == 1.0 && needsign )
   {
      appendLine(scip, file, linebuffer, linecnt, " + ");
   }
   else if( coef == -1.0 )
   {
      appendLine(scip, file, linebuffer, linecnt, " - ");
   }
   else
   {
      appendLine(scip, file, linebuffer, linecnt, " ");
   }

   if( SCIPisExprVar(scip, expr) )
   {
      appendLine(scip, file, linebuffer, linecnt, SCIPvarGetName(SCIPgetVarExprVar(expr)));
      return;
   }

   if( SCIPisExprValue(scip, expr) )
   {
      if( REALABS(coef) == 1.0 )
      {
         /* in this case, we will have printed only a sign or space above, so print also a 1.0 */
         appendLine(scip, file, linebuffer, linecnt, "1.0");
      }
      return;
   }

   if( SCIPisExprPower(scip, expr) )
   {
      assert(SCIPisExprVar(scip, SCIPexprGetChildren(expr)[0]));

      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, "%s^%.15g", SCIPvarGetName(SCIPgetVarExprVar(SCIPexprGetChildren(expr)[0])), SCIPgetExponentExprPow(expr));
      appendLine(scip, file, linebuffer, linecnt, buffer);

      return;
   }

   assert(SCIPisExprProduct(scip, expr));
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      child = SCIPexprGetChildren(expr)[c];

      if( c > 0 )
         appendLine(scip, file, linebuffer, linecnt, " ");

      if( SCIPisExprVar(scip, child) )
      {
         appendLine(scip, file, linebuffer, linecnt, SCIPvarGetName(SCIPgetVarExprVar(child)));
         continue;
      }

      assert(SCIPisExprPower(scip, child));
      assert(SCIPisExprVar(scip, SCIPexprGetChildren(child)[0]));

      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, "%s^%.15g", SCIPvarGetName(SCIPgetVarExprVar(SCIPexprGetChildren(child)[0])), SCIPgetExponentExprPow(child));
      appendLine(scip, file, linebuffer, linecnt, buffer);
   }
}

/** print polynomial row in PIP format to file stream */
static
void printRowNl(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=", "<=", or ">=") */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   char consname[PIP_MAX_NAMELEN + 1]; /* an extra character for ':' */
   char buffer[PIP_MAX_PRINTLEN];
   char linebuffer[PIP_MAX_PRINTLEN+1] = { '\0' };
   int linecnt;

   assert(scip != NULL);
   assert(strcmp(type, "=") == 0 || strcmp(type, "<=") == 0 || strcmp(type, ">=") == 0);
   assert(expr != NULL);

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   if( strlen(rowname) > 0 || strlen(rownameextension) > 0 )
   {
      (void) SCIPsnprintf(consname, PIP_MAX_NAMELEN + 1, "%s%s:", rowname, rownameextension);
      appendLine(scip, file, linebuffer, &linecnt, consname);
   }

   if( SCIPisExprSum(scip, expr) )
   {
      int c;
      SCIP_Bool needsign = FALSE;

      if( SCIPgetConstantExprSum(expr) != 0.0 )
      {
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g", SCIPgetConstantExprSum(expr));
         appendLine(scip, file, linebuffer, &linecnt, buffer);

         needsign = TRUE;
      }

      for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      {
         printSignomial(scip, file, linebuffer, &linecnt, SCIPexprGetChildren(expr)[c], SCIPgetCoefsExprSum(expr)[c], needsign);
         needsign = TRUE;
      }
   }
   else
   {
      printSignomial(scip, file, linebuffer, &linecnt, expr, 1.0, FALSE);
   }

   /* print right hand side */
   (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %s %+.15g", type, rhs);

   /* we start a new line; therefore we tab this line */
   if( linecnt == 0 )
      appendLine(scip, file, linebuffer, &linecnt, " ");
   appendLine(scip, file, linebuffer, &linecnt, buffer);

   endLine(scip, file, linebuffer, &linecnt);
}

/** print "and" constraint as row in PIP format to file stream */
static
void printRowAnd(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   SCIP_CONS*            cons                /**< "and" constraint */
   )
{
   char linebuffer[PIP_MAX_PRINTLEN+1] = { '\0' };
   int linecnt;
   int i;

   assert(scip != NULL);
   assert(rowname != NULL);
   assert(cons != NULL);

   clearLine(linebuffer, &linecnt);

   /* start each line with a space */
   appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print row name */
   if( strlen(rowname) > 0 )
   {
      appendLine(scip, file, linebuffer, &linecnt, rowname);
      appendLine(scip, file, linebuffer, &linecnt, ":");
   }

   for( i = 0; i < SCIPgetNVarsAnd(scip, cons); ++i )
   {
      appendLine(scip, file, linebuffer, &linecnt, " ");
      appendLine(scip, file, linebuffer, &linecnt, SCIPvarGetName(SCIPgetVarsAnd(scip, cons)[i]));
   }

   appendLine(scip, file, linebuffer, &linecnt, " - ");
   appendLine(scip, file, linebuffer, &linecnt, SCIPvarGetName(SCIPgetResultantAnd(scip, cons)));

   /* we start a new line; therefore we tab this line */
   if( linecnt == 0 )
      appendLine(scip, file, linebuffer, &linecnt, " ");

   /* print right hand side */
   appendLine(scip, file, linebuffer, &linecnt, " = 0");

   endLine(scip, file, linebuffer, &linecnt);
}

/** prints given (linear or) quadratic constraint information in LP format to file stream */
static
SCIP_RETCODE printQuadraticCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficients values (or NULL if all linear coefficient values are 1) */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_EXPR*            quadexpr,           /**< quadratic expression (or NULL if nlinvars > 0) */
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

   assert( nlinvars == 0 || linvars != NULL );
   assert( quadexpr == NULL || nlinvars == 0);
   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   nactivevars = nlinvars;
   if( nlinvars > 0 )
   {
      /* duplicate variable and value array */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, linvars, nactivevars ) );
      if( linvals != NULL )
      {
         SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, linvals, nactivevars ) );
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

   /* print row(s) in LP format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equal constraint */
      SCIP_CALL( printRow(scip, file, rowname, "", "=", activevars, activevals, nactivevars, quadexpr,
         rhs - activeconstant, transformed) );
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", ">=", activevars,
            activevals, nactivevars, quadexpr, lhs - activeconstant, transformed) );
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         SCIP_CALL( printRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "<=", activevars,
            activevals, nactivevars, quadexpr, rhs - activeconstant, transformed) );
      }
   }

   if( nlinvars > 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &activevars);
      SCIPfreeBufferArray(scip, &activevals);
   }

   return SCIP_OKAY;
}

/** prints given nonlinear constraint information in LP format to file stream */
static
SCIP_RETCODE printNonlinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< name of the row */
   SCIP_EXPR*            expr,               /**< polynomial expression */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   assert(scip != NULL);
   assert(rowname != NULL);
   assert(expr != NULL);
   assert(lhs <= rhs);

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   /* print row(s) in LP format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equal constraint */
      printRowNl(scip, file, rowname, "", "=", expr, rhs);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRowNl(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", ">=", expr, lhs);
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         printRowNl(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "<=", expr, rhs);
      }
   }

   return SCIP_OKAY;
}

/** check whether given variables are aggregated and put them into an array without duplication */
static
SCIP_RETCODE collectAggregatedVars(
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_VAR**            vars,               /**< variable array */
   int*                  nAggregatedVars,    /**< number of aggregated variables on output */
   SCIP_VAR***           aggregatedVars,     /**< array storing the aggregated variables on output */
   SCIP_HASHTABLE**      varAggregated       /**< hashtable for checking duplicates */
   )
{
   int j;

   /* check variables */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VARSTATUS status;
      SCIP_VAR* var;

      var = vars[j];
      status = SCIPvarGetStatus(var);

      /* collect aggregated variables in a list */
      if( status >= SCIP_VARSTATUS_AGGREGATED )
      {
         assert( status == SCIP_VARSTATUS_AGGREGATED ||
            status == SCIP_VARSTATUS_MULTAGGR ||
            status == SCIP_VARSTATUS_NEGATED );

         if ( ! SCIPhashtableExists(*varAggregated, (void*) var) )
         {
            (*aggregatedVars)[(*nAggregatedVars)++] = var;
            SCIP_CALL( SCIPhashtableInsert(*varAggregated, (void*) var) );
         }
      }
   }

   return SCIP_OKAY;
}


/** print aggregated variable-constraints */
static
SCIP_RETCODE printAggregatedCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   int                   nvars,              /**< number of active variables in the problem */
   int                   nAggregatedVars,    /**< number of aggregated variables */
   SCIP_VAR**            aggregatedVars      /**< array storing the aggregated variables */
   )
{
   int j;

   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant = 0.0;
   char consname[PIP_MAX_NAMELEN];

   assert( scip != NULL );

   /* write aggregation constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &activevars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nvars) );

   for (j = 0; j < nAggregatedVars; ++j)
   {
      /* set up list to obtain substitution variables */
      nactivevars = 1;

      activevars[0] = aggregatedVars[j];
      activevals[0] = 1.0;
      activeconstant = 0.0;

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

      activevals[nactivevars] = -1.0;
      activevars[nactivevars] = aggregatedVars[j];
      ++nactivevars;

      /* output constraint */
      (void) SCIPsnprintf(consname, PIP_MAX_NAMELEN, "aggr_%s", SCIPvarGetName(aggregatedVars[j]));
      SCIP_CALL( printRow(scip, file, consname, "", "=", activevars, activevals, nactivevars, NULL, - activeconstant,
         transformed) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

/** returns whether name is valid according to PIP specification
 *
 * Checks these two conditions from http://polip.zib.de/pipformat.php:
 * - Names/labels can contain at most 255 characters.
 * - Name/labels have to consist of the following characters: a-z, A-Z, 0-9, "!", "#", "$", "%", "&", ";", "?", "@", "_". They cannot start with a number.
 *
 * In addition checks that the length is not zero.
 */
static
SCIP_Bool isNameValid(
   const char*           name                /**< name to check */
   )
{
   size_t len;
   size_t i;

   assert(name != NULL);

   len = strlen(name);  /*lint !e613*/
   if( len > (size_t) PIP_MAX_NAMELEN || len == 0 )
      return FALSE;

   /* names cannot start with a number */
   if( isdigit(name[0]) )
      return FALSE;

   for( i = 0; i < len; ++i )
   {
      /* a-z, A-Z, 0-9 are ok */
      if( isalnum(name[i]) )
         continue;

      /* characters in namechars are ok, too */
      if( strchr(namechars, name[i]) != NULL )
         continue;

      return FALSE;
   }

   return TRUE;
}


/** method check if the variable names are valid according to PIP specification */
static
void checkVarnames(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   nvars               /**< number of variables */
   )
{
   int v;

   assert(scip != NULL);
   assert(vars != NULL || nvars == 0);

   /* check if the variable names are not too long and have only characters allowed by PIP */
   for( v = 0; v < nvars; ++v )
   {
      if( !isNameValid(SCIPvarGetName(vars[v])) )
      {
         SCIPwarningMessage(scip, "variable name <%s> is not valid (too long or disallowed characters); PIP might be corrupted\n", SCIPvarGetName(vars[v]));
         return;
      }
   }
}

/** method check if the constraint names are valid according to PIP specification */
static
void checkConsnames(
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

   assert( scip != NULL );
   assert( conss != NULL || nconss == 0 );

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL); /* for lint */
      cons = conss[c];
      assert(cons != NULL );

      /* in case the transformed is written only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( !isNameValid(SCIPconsGetName(cons)) )
      {
         SCIPwarningMessage(scip, "constraint name <%s> is not valid (too long or unallowed characters); PIP might be corrupted\n", SCIPconsGetName(cons));
         return;
      }

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_Real lhs = SCIPgetLhsLinear(scip, cons);
         SCIP_Real rhs = SCIPgetRhsLinear(scip, cons);

         /* for ranged constraints, we need to be able to append _lhs and _rhs to the constraint name, so need additional 4 characters */
         if( !SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > (size_t) PIP_MAX_NAMELEN -  4 )
         {
            SCIPwarningMessage(scip, "name of ranged constraint <%s> has to be cut down to %d characters;\n", SCIPconsGetName(conss[c]),
               PIP_MAX_NAMELEN  - 1);
            return;
         }
      }
   }
}

/** writes problem to file
 * @todo add writing cons_pseudoboolean
 */
SCIP_RETCODE SCIPwritePip(
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
   char linebuffer[PIP_MAX_PRINTLEN+1];

   char varname[PIP_MAX_NAMELEN];
   char buffer[PIP_MAX_PRINTLEN];

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;
   SCIP_CONS** consNonlinear;
   int nConsNonlinear;
   SCIP_CONS** consAnd;
   int nConsAnd;
   char consname[PIP_MAX_NAMELEN];

   SCIP_VAR** aggregatedVars;
   int nAggregatedVars;
   SCIP_HASHTABLE* varAggregated;

   SCIP_VAR** tmpvars;
   int tmpvarssize;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;

   assert( scip != NULL );

   nAggregatedVars = 0;
   nConsNonlinear = 0;
   nConsAnd = 0;

   /* check if the variable names are not to long */
   checkVarnames(scip, vars, nvars);

   /* check if the constraint names are to long */
   checkConsnames(scip, conss, nconss, transformed);

   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "\\ SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "\\   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "\\   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "\\   Constraints      : %d\n", nconss);

   /* print objective sense */
   SCIPinfoMessage(scip, file, "%s\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "Minimize" : "Maximize");

   clearLine(linebuffer, &linecnt);
   appendLine(scip, file, linebuffer, &linecnt, " Obj:");

   for (v = 0; v < nvars; ++v)
   {
      var = vars[v];

#ifndef NDEBUG
      /* in case the original problem has to be posted the variables have to be either "original" or "negated" */
      if ( !transformed )
         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED );
#endif

      if ( SCIPisZero(scip, SCIPvarGetObj(var)) )
         continue;

      /* we start a new line; therefore we tab this line */
      if ( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, "     ");

      (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", objscale * SCIPvarGetObj(var), varname );

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   if( ! SCIPisZero(scip, objoffset) )
   {
      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g", objscale * objoffset);
      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   endLine(scip, file, linebuffer, &linecnt);

   /* print "Subject to" section */
   SCIPinfoMessage(scip, file, "Subject to\n");

   /* collect quadratic, nonlinear, absolute power, and, and bivariate constraints in arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &consNonlinear, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consAnd, nconss) );

   tmpvarssize = SCIPgetNTotalVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, tmpvarssize) );

   for (c = 0; c < nconss; ++c)
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed is written only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      (void) SCIPsnprintf(consname, PIP_MAX_NAMELEN, "%s", SCIPconsGetName(cons));
      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
               NULL, SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed) );
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 1.0, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, -SCIPinfinity(scip), 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 1.0, SCIPinfinity(scip), transformed) );
            break;
         }
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
               NULL, 1.0, SCIPinfinity(scip), transformed) );
      }
      else if ( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
            consvals[v] = (SCIP_Real)weights[v];

         SCIP_CALL( printQuadraticCons(scip, file, consname, consvars, consvals, nconsvars,
               NULL, -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvals);
      }
      else if ( strcmp(conshdlrname, "varbound") == 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( printQuadraticCons(scip, file, consname, consvars, consvals, 2, NULL,
            SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "nonlinear") == 0 )
      {
         SCIP_Bool ispolynomial;
         SCIP_Bool isquadratic;
         SCIP_EXPR* simplifiedexpr = NULL;

         ispolynomial = isExprPolynomial(scip, SCIPgetExprNonlinear(cons));
         if( !ispolynomial )
         {
            /* simplify expression and check again if polynomial
             * simplifying the expr owned by the cons can have undesired sideffects onto the consdata (the varhashmap can get messed up), so we copy first
             */
            SCIP_EXPR* exprcopy;
            SCIP_Bool changed;
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPduplicateExpr(scip, SCIPgetExprNonlinear(cons), &exprcopy, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPsimplifyExpr(scip, exprcopy, &simplifiedexpr, &changed, &infeasible, NULL, NULL) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprcopy) );

            ispolynomial = isExprPolynomial(scip, simplifiedexpr);
         }

         /* nonlinear constraints that are not polynomial cannot be printed as PIP */
         if( !ispolynomial )
         {
            SCIPwarningMessage(scip, "nonlinear constraint <%s> is not polynomial\n", SCIPconsGetName(cons));
            SCIPinfoMessage(scip, file, "\\ ");
            SCIP_CALL( SCIPprintCons(scip, cons, file) );
            SCIPinfoMessage(scip, file, ";\n");

            SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedexpr) );
            return SCIP_OKAY;
         }

         /* check whether constraint is even quadratic
          * (we could also skip this and print as polynomial, but the code exists already)
          */
         SCIP_CALL( SCIPcheckExprQuadratic(scip, simplifiedexpr != NULL ? simplifiedexpr : SCIPgetExprNonlinear(cons), &isquadratic) );
         if( isquadratic )
            isquadratic = SCIPexprAreQuadraticExprsVariables(simplifiedexpr != NULL ? simplifiedexpr : SCIPgetExprNonlinear(cons));

         if( isquadratic )
         {
            SCIP_CALL( printQuadraticCons(scip, file, consname, NULL, NULL, 0, simplifiedexpr != NULL ? simplifiedexpr : SCIPgetExprNonlinear(cons),
               SCIPgetLhsNonlinear(cons), SCIPgetRhsNonlinear(cons), transformed) );
         }
         else
         {
            SCIP_CALL( printNonlinearCons(scip, file, consname, simplifiedexpr != NULL ? simplifiedexpr : SCIPgetExprNonlinear(cons), SCIPgetLhsNonlinear(cons), SCIPgetRhsNonlinear(cons)) );
         }

         consNonlinear[nConsNonlinear++] = cons;

         if( simplifiedexpr != NULL )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedexpr) );
         }
      }
      else if( strcmp(conshdlrname, "and") == 0 )
      {
         printRowAnd(scip, file, consname, cons);

         consAnd[nConsAnd++] = cons;
      }
      else
      {
         SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
         SCIPinfoMessage(scip, file, "\\ ");
         SCIP_CALL( SCIPprintCons(scip, cons, file) );
         SCIPinfoMessage(scip, file, ";\n");
      }
   }

   /* create hashtable for storing aggregated variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &aggregatedVars, nvars) );
   SCIP_CALL( SCIPhashtableCreate(&varAggregated, SCIPblkmem(scip), nvars/10, hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );

   /* check for aggregated variables in nonlinear constraints and output aggregations as linear constraints */
   for( c = 0; c < nConsNonlinear; ++c )
   {
      SCIP_Bool success;
      int ntmpvars;

      /* get variables of the nonlinear constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, consNonlinear[c], &ntmpvars, &success) );
      assert(success);
      if( ntmpvars > tmpvarssize )
      {
         tmpvarssize = SCIPcalcMemGrowSize(scip, ntmpvars);
         SCIP_CALL( SCIPreallocBufferArray(scip, &tmpvars, tmpvarssize) );
      }
      SCIP_CALL( SCIPgetConsVars(scip, consNonlinear[c], tmpvars, tmpvarssize, &success) );
      assert(success);

      SCIP_CALL( collectAggregatedVars(ntmpvars, tmpvars, &nAggregatedVars, &aggregatedVars, &varAggregated) );
   }

   /* check for aggregated variables in and constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsAnd; ++c)
   {
      SCIP_VAR* resultant;

      cons = consAnd[c];

      SCIP_CALL( collectAggregatedVars(SCIPgetNVarsAnd(scip, cons), SCIPgetVarsAnd(scip, cons), &nAggregatedVars, &aggregatedVars, &varAggregated) );

      resultant = SCIPgetResultantAnd(scip, cons);
      SCIP_CALL( collectAggregatedVars(1, &resultant, &nAggregatedVars, &aggregatedVars, &varAggregated) );
   }

   /* print aggregation constraints */
   SCIP_CALL( printAggregatedCons(scip, file, transformed, nvars, nAggregatedVars, aggregatedVars) );

   /* print "Bounds" section */
   SCIPinfoMessage(scip, file, "Bounds\n");
   for (v = 0; v < nvars; ++v)
   {
      var = vars[v];
      assert( var != NULL );
      (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

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

      if ( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
         SCIPinfoMessage(scip, file, " %s free\n", varname);
      else
      {
         /* print lower bound */
         if ( SCIPisInfinity(scip, -lb) )
            SCIPinfoMessage(scip, file, " -inf <= ");
         else
         {
            if ( SCIPisZero(scip, lb) )
            {
               /* variables are nonnegative by default - so we skip these variables */
               if ( SCIPisInfinity(scip, ub) )
                  continue;
               lb = 0.0;
            }

            SCIPinfoMessage(scip, file, " %.15g <= ", lb);
         }
         /* print variable name */
         SCIPinfoMessage(scip, file, "%s", varname);

         /* print upper bound as far this one is not infinity */
         if( !SCIPisInfinity(scip, ub) )
            SCIPinfoMessage(scip, file, " <= %.15g", ub);

         SCIPinfoMessage(scip, file, "\n");
      }
   }

   /* output aggregated variables as 'free' */
   for (v = 0; v < nAggregatedVars; ++v)
   {
      var = aggregatedVars[v];
      assert( var != NULL );
      (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

      SCIPinfoMessage(scip, file, " %s free\n", varname);
   }

   /* free space */
   SCIPfreeBufferArray(scip, &aggregatedVars);
   SCIPhashtableFree(&varAggregated);

   /* print binaries section */
   if ( nbinvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Binaries\n");

      clearLine(linebuffer, &linecnt);

      for (v = 0; v < nvars; ++v)
      {
         var = vars[v];
         assert( var != NULL );

         if ( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
         {
            (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %s", varname);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }
      }

      endLine(scip, file, linebuffer, &linecnt);
   }

   /* print generals section */
   if ( nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "Generals\n");

      for (v = 0; v < nvars; ++v)
      {
         var = vars[v];
         assert( var != NULL );

         if ( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
         {
            (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var) );
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %s", varname);
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }
      }
      endLine(scip, file, linebuffer, &linecnt);
   }

   /* free space */
   SCIPfreeBufferArray(scip, &tmpvars);
   SCIPfreeBufferArray(scip, &consNonlinear);
   SCIPfreeBufferArray(scip, &consAnd);

   /* end of lp format */
   SCIPinfoMessage(scip, file, "%s\n", "End");

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyPip)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderPip(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadPip)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadPip(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWritePip)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPwritePip(scip, file, name, transformed, objsense, objscale, objoffset, vars,
         nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the pip file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderPip(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyPip) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadPip) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWritePip) );

   return SCIP_OKAY;
}


/** reads problem from file */
SCIP_RETCODE SCIPreadPip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{  /*lint --e{715}*/
   PIPINPUT pipinput;
   SCIP_RETCODE retcode;
   int i;

   assert(scip != NULL);  /* for lint */
   assert(reader != NULL);

   /* initialize PIP input data */
   pipinput.file = NULL;
   pipinput.linebuf[0] = '\0';
   pipinput.probname[0] = '\0';
   pipinput.objname[0] = '\0';
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pipinput.token, PIP_MAX_LINELEN) ); /*lint !e506*/
   pipinput.token[0] = '\0';
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pipinput.tokenbuf, PIP_MAX_LINELEN) ); /*lint !e506*/
   pipinput.tokenbuf[0] = '\0';
   for( i = 0; i < PIP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((pipinput.pushedtokens)[i]), PIP_MAX_LINELEN) ); /*lint !e866 !e506*/
   }

   pipinput.npushedtokens = 0;
   pipinput.linenumber = 0;
   pipinput.linepos = 0;
   pipinput.section = PIP_START;
   pipinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   pipinput.haserror = FALSE;

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &(pipinput.initialconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicconss", &(pipinput.dynamicconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &(pipinput.dynamiccols)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &(pipinput.dynamicrows)) );

   /* read the file */
   retcode = readPIPFile(scip, &pipinput, filename);

   /* free dynamically allocated memory */
   for( i = PIP_MAX_PUSHEDTOKENS - 1; i >= 0 ; --i )
   {
      SCIPfreeBlockMemoryArray(scip, &pipinput.pushedtokens[i], PIP_MAX_LINELEN);
   }
   SCIPfreeBlockMemoryArray(scip, &pipinput.tokenbuf, PIP_MAX_LINELEN);
   SCIPfreeBlockMemoryArray(scip, &pipinput.token, PIP_MAX_LINELEN);

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   /* evaluate the result */
   if( pipinput.haserror )
      retcode = SCIP_READERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, pipinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}
