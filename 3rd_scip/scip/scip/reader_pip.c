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

/**@file   reader_pip.c
 * @brief  file reader for polynomial mixed-integer programs in PIP format
 * @author Stefan Vigerske
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/ /* needed for strncasecmp() */
#endif
#include <ctype.h>

#include "scip/reader_pip.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_abspower.h"
#include "scip/cons_and.h"
#include "scip/cons_bivariate.h"
#include "scip/pub_misc.h"

#define READER_NAME             "pipreader"
#define READER_DESC             "file reader for polynomial mixed-integer programs in PIP format"
#define READER_EXTENSION        "pip"


/*
 * Data structures
 */
#define PIP_MAX_LINELEN        65536
#define PIP_MAX_PUSHEDTOKENS   2
#define PIP_INIT_VARSSIZE      256
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
   SCIP_EXPRDATA_MONOMIAL*** monomials,      /**< pointer to current array of monomials */
   int*                  monomialssize,      /**< current size of monomials array at input; new size at exit */
   int                   minnmonomials       /**< required minimal size of monomials array */
   )
{
   int newsize;

   assert(scip != NULL);
   assert(monomials != NULL);
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
   *monomialssize = newsize;

   return SCIP_OKAY;
}

/** ensure that arrays of exponents and variable indices can hold a minimum number of entries */
static
SCIP_RETCODE ensureFactorsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real**           exponents,          /**< pointer to current array of exponents */
   int**                 varidxs,            /**< pointer to current array of variable indices */
   int*                  factorssize,        /**< current size of arrays at input; new size at exit */
   int                   minnfactors         /**< required minimal size of arrays */
   )
{
   int newsize;

   assert(scip != NULL);
   assert(exponents != NULL);
   assert(varidxs != NULL);
   assert(factorssize != NULL);
   assert(*exponents != NULL || *factorssize == 0);
   assert(*varidxs   != NULL || *factorssize == 0);
   assert((*exponents != NULL) == (*varidxs != NULL));

   if( minnfactors <= *factorssize )
      return SCIP_OKAY;

   newsize = SCIPcalcMemGrowSize(scip, minnfactors);

   if( *exponents != NULL )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, exponents, newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, varidxs,   newsize) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, exponents, newsize) );
      SCIP_CALL( SCIPallocBufferArray(scip, varidxs,   newsize) );
   }
   *factorssize = newsize;

   return SCIP_OKAY;
}

/** gives index of variable in vars array, inserts it at the end if not existing yet */
static
SCIP_RETCODE getVariableIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to current array of variables */
   int*                  varssize,           /**< current size of variables array at input; new size at exit */
   int*                  nvars,              /**< number of variables stored in array */
   SCIP_HASHMAP*         varhash,            /**< hashmap variables -> indices */
   SCIP_VAR*             var,                /**< the variable which index we need */
   int*                  varidx              /**< pointer to store index of variable in *vars */
   )
{
   assert(scip != NULL);
   assert(varssize != NULL);
   assert(vars != NULL);
   assert(*vars != NULL || *varssize == 0);
   assert(nvars != NULL);
   assert(*nvars <= *varssize);
   assert(varhash != NULL);
   assert(var != NULL);
   assert(varidx != NULL);

   /* check if we saw this variable before */
   if( SCIPhashmapExists(varhash, (void*)var) )
   {
      *varidx = (int)(size_t)SCIPhashmapGetImage(varhash, (void*)var);
      assert(*varidx >= 0);
      assert(*varidx < *nvars);

      return SCIP_OKAY;
   }

   /* since variable is new, add it to the end of vars array and into hashmap */

   /* ensure enough space in vars array */
   if( *nvars + 1 > *varssize )
   {
      *varssize = SCIPcalcMemGrowSize(scip, *nvars + 1);
      if( *vars == NULL )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, vars, *varssize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, *varssize) );
      }
   }
   assert(*vars != NULL);  /*lint !e613*/

   (*vars)[*nvars] = var;  /*lint !e613*/
   SCIP_CALL( SCIPhashmapInsert(varhash, (void*)var, (void*)(size_t)*nvars) );
   *varidx = *nvars;

   ++*nvars;

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput,           /**< PIP reading data */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   PIP_MAX_LINELEN */
   SCIP_EXPRTREE**       exprtree,           /**< pointer to store constraint function as polynomial expression */
   int*                  degree,             /**< pointer to store degree of polynomial */
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_EXPR* expression;
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Real coef;
   int coefsign;
   int nextcoefsign;
   int monomialdegree;
   SCIP_EXPR** varexprs;
   int i;

   SCIP_VAR** vars;
   int varssize;
   int nvars;
   SCIP_HASHMAP* varhash;

   SCIP_Real constant;

   SCIP_EXPRDATA_MONOMIAL** monomials;
   int monomialssize;
   int nmonomials;

   int nfactors;
   int factorssize;
   SCIP_Real* exponents;
   int* varidxs;

   assert(scip != NULL);
   assert(pipinput != NULL);
   assert(name != NULL);
   assert(exprtree != NULL);
   assert(degree != NULL);
   assert(newsection != NULL);

   *name = '\0';
   *exprtree = NULL;
   *degree = 0;
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
            (void)strncpy(name, pipinput->tokenbuf, PIP_MAX_LINELEN);
            name[PIP_MAX_LINELEN - 1] = '\0';
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

   /* initialize buffer for storing the variables */
   varssize = PIP_INIT_VARSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip), PIP_INIT_VARSSIZE) );

   /* initialize buffer for storing the monomials */
   monomialssize = PIP_INIT_MONOMIALSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &monomials, monomialssize) );

   /* initialize buffer for storing the factors in a monomial */
   factorssize = PIP_INIT_FACTORSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, &exponents, factorssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varidxs,   factorssize) );

   /* read the coefficients */
   coefsign = +1;
   nextcoefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   nmonomials = 0;
   nvars = 0;
   nfactors = 0;
   monomialdegree = 0;
   constant = 0.0;
   while( getNextToken(scip, pipinput) )
   {
      SCIP_VAR* var;
      int varidx;
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
            SCIP_CALL( ensureMonomialsSize(scip, &monomials, &monomialssize, nmonomials + 1) );
            SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip),
                  &monomials[nmonomials], coefsign * coef, nfactors, varidxs, exponents) );
            ++nmonomials;
         }
         else if( havevalue) 
         {
            constant += coefsign * coef;
         }

         if( monomialdegree > *degree )
            *degree = monomialdegree;

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
         exponents[nfactors-1] = exponent;
         if( SCIPisIntegral(scip, exponent) && exponent > 0.0 )
            monomialdegree += (int)exponent - 1; /* -1, because we added +1 when we put the variable into varidxs */
         else
            monomialdegree = SCIP_EXPR_DEGREEINFINITY;

         SCIPdebugMsg(scip, "(line %d) read exponent value %g for variable %s\n", pipinput->linenumber, exponent,
            SCIPvarGetName(vars[varidxs[nfactors-1]]));
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

      /* get the index of the variable in the vars array, or add there if not in it yet */
      SCIP_CALL( getVariableIndex(scip, &vars, &varssize, &nvars, varhash, var, &varidx) );

      SCIP_CALL( ensureFactorsSize(scip, &exponents, &varidxs, &factorssize, nfactors + 1) );

      exponents[nfactors] = 1.0;
      varidxs[nfactors]   = varidx;
      ++nfactors;
      ++monomialdegree;
   }

   if( nfactors > 0 )
   {
      syntaxError(scip, pipinput, "string ended before monomial has finished");
      goto TERMINATE_READPOLYNOMIAL;
   }

   /* create variable expressions */
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexprs[i], SCIP_EXPR_VARIDX, i) );
   }

   /* create polynomial expression, let polynomial take over ownership of monomials */
   SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), &expression, nvars, varexprs,
         nmonomials, monomials, constant, FALSE) );

   SCIPfreeBufferArray(scip, &varexprs);

   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), exprtree, expression, 0, 0, NULL) );
   SCIP_CALL( SCIPexprtreeSetVars(*exprtree, nvars, vars) );

   SCIPdebugMsg(scip, "read polynomial of degree %d: ", *degree);
   SCIPdebug( SCIP_CALL( SCIPexprtreePrintWithNames(*exprtree, SCIPgetMessagehdlr(scip), NULL) ) );
   SCIPdebugMsgPrint(scip, "\n");

 TERMINATE_READPOLYNOMIAL:
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &monomials);
   SCIPfreeBufferArray(scip, &exponents);
   SCIPfreeBufferArray(scip, &varidxs);
   SCIPhashmapFree(&varhash);

   return SCIP_OKAY;
}

/** given an expression tree that holds a polynomial expression of degree at most two,
 * gives the coefficients of the constant, linear, and quadratic part of this expression
 */
static
void getLinearAndQuadraticCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree holding polynomial expression */
   SCIP_Real*            constant,           /**< buffer to store constant monomials */
   int*                  nlinvars,           /**< buffer to store number of linear coefficients */
   SCIP_VAR**            linvars,            /**< array to fill with linear variables */
   SCIP_Real*            lincoefs,           /**< array to fill with coefficients of linear variables */
   int*                  nquadterms,         /**< buffer to store number of quadratic terms */ 
   SCIP_VAR**            quadvars1,          /**< array to fill with quadratic variables */
   SCIP_VAR**            quadvars2,          /**< array to fill with quadratic variables */
   SCIP_Real*            quadcoefs           /**< array to fill with coefficients of quadratic terms */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPRDATA_MONOMIAL** monomials;
   int nmonomials;
   int varidx;
   int i;

   expr = SCIPexprtreeGetRoot(exprtree);
   assert(expr != NULL);
   assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);
   assert(SCIPexprGetNChildren(expr) == SCIPexprtreeGetNVars(exprtree));

   nmonomials = SCIPexprGetNMonomials(expr);
   monomials  = SCIPexprGetMonomials(expr);

   *constant = SCIPexprGetPolynomialConstant(expr);
   *nlinvars = 0;
   *nquadterms = 0;
   for( i = 0; i < nmonomials; ++i )
   {
      assert(SCIPexprGetMonomialNFactors(monomials[i]) >= 0);
      assert(SCIPexprGetMonomialNFactors(monomials[i]) <= 2);
      assert(SCIPexprGetMonomialExponents(monomials[i]) != NULL    || SCIPexprGetMonomialNFactors(monomials[i]) == 0);
      assert(SCIPexprGetMonomialChildIndices(monomials[i]) != NULL || SCIPexprGetMonomialNFactors(monomials[i]) == 0);

      if( SCIPexprGetMonomialNFactors(monomials[i]) == 0 )
      {
         /* constant monomial */
         *constant += SCIPexprGetMonomialCoef(monomials[i]);
      }
      else if( SCIPexprGetMonomialNFactors(monomials[i]) == 1 && SCIPexprGetMonomialExponents(monomials[i])[0] == 1.0 )
      {
         /* linear monomial */
         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */

         lincoefs[*nlinvars] = SCIPexprGetMonomialCoef(monomials[i]);
         linvars[*nlinvars]  = SCIPexprtreeGetVars(exprtree)[varidx];
         ++*nlinvars;
      }
      else if( SCIPexprGetMonomialNFactors(monomials[i]) == 1 )
      {
         /* square monomial */
         assert(SCIPexprGetMonomialExponents(monomials[i])[0] == 2.0);

         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */

         quadcoefs[*nquadterms] = SCIPexprGetMonomialCoef(monomials[i]);
         quadvars1[*nquadterms] = SCIPexprtreeGetVars(exprtree)[varidx];
         quadvars2[*nquadterms] = quadvars1[*nquadterms];
         ++*nquadterms;
      }
      else
      {
         /* bilinear monomial */
         assert(SCIPexprGetMonomialExponents(monomials[i])[0] == 1.0);
         assert(SCIPexprGetMonomialExponents(monomials[i])[1] == 1.0);

         quadcoefs[*nquadterms] = SCIPexprGetMonomialCoef(monomials[i]);

         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */
         quadvars1[*nquadterms] = SCIPexprtreeGetVars(exprtree)[varidx];

         varidx = SCIPexprGetMonomialChildIndices(monomials[i])[1];
         assert(varidx >= 0);
         assert(varidx < SCIPexprtreeGetNVars(exprtree));
         assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */
         quadvars2[*nquadterms] = SCIPexprtreeGetVars(exprtree)[varidx];

         ++*nquadterms;
      }
   }
}

/** reads the objective section */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   char name[PIP_MAX_LINELEN];
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;
   int degree;
   SCIP_Bool newsection;
   int varidx;
   int nmonomials;
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
   SCIP_CALL( readPolynomial(scip, pipinput, name, &exprtree, &degree, &newsection) );
   if( !hasError(pipinput) && exprtree != NULL )
   {
      int i;

      expr = SCIPexprtreeGetRoot(exprtree);
      assert(expr != NULL);
      assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);

      nmonomials = SCIPexprGetNMonomials(expr);

      if( SCIPexprGetPolynomialConstant(expr) != 0.0 )
      {
         SCIP_VAR* objconst;
         SCIP_CALL( SCIPcreateVarBasic(scip, &objconst, "objconst", 1.0, 1.0, SCIPexprGetPolynomialConstant(expr), SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, objconst) );
         SCIP_CALL( SCIPreleaseVar(scip, &objconst) );
      }

      assert(degree >= 0);
      if( degree == 1 )
      {
         SCIP_Real coef;
         SCIP_VAR* var;
         SCIP_EXPRDATA_MONOMIAL** monomials;

         assert(SCIPexprtreeGetVars(exprtree) != NULL);
         assert(SCIPexprGetNChildren(expr) == SCIPexprtreeGetNVars(exprtree));

         monomials  = SCIPexprGetMonomials(expr);

         for( i = 0; i < nmonomials; ++i )
         {
            assert(SCIPexprGetMonomialNFactors(monomials[i]) == 1);
            assert(SCIPexprGetMonomialExponents(monomials[i]) != NULL);
            assert(SCIPexprGetMonomialExponents(monomials[i])[0] == 1.0);
            assert(SCIPexprGetMonomialChildIndices(monomials[i]) != NULL);

            varidx = SCIPexprGetMonomialChildIndices(monomials[i])[0];
            assert(varidx >= 0);
            assert(varidx < SCIPexprGetNChildren(expr));
            assert(SCIPexprGetOperator(SCIPexprGetChildren(expr)[varidx]) == SCIP_EXPR_VARIDX);
            assert(SCIPexprGetOpIndex(SCIPexprGetChildren(expr)[varidx]) == varidx); /* assume that child varidx corresponds to variable varidx */

            coef = SCIPexprGetMonomialCoef(monomials[i]);
            var = SCIPexprtreeGetVars(exprtree)[varidx];

            SCIP_CALL( SCIPchgVarObj(scip, var, SCIPvarGetObj(var) + coef) );
         }
      }
      else if( degree == 2 )
      {
         /* insert dummy variable and constraint to represent quadratic part of objective */

         SCIP_VAR*  quadobjvar;
         SCIP_CONS* quadobjcons;
         SCIP_Real  lhs;
         SCIP_Real  rhs;

         SCIP_Real constant;
         int nlinvars;
         SCIP_VAR** linvars;
         SCIP_Real* lincoefs;
         int nquadterms;
         SCIP_VAR** quadvars1;
         SCIP_VAR** quadvars2;
         SCIP_Real* quadcoefs;

         SCIP_CALL( SCIPallocBufferArray(scip, &linvars,   nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs,  nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, nmonomials) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, nmonomials) );

         getLinearAndQuadraticCoefs(scip, exprtree, &constant, &nlinvars, linvars, lincoefs, &nquadterms, quadvars1, quadvars2, quadcoefs);

         SCIP_CALL( SCIPcreateVar(scip, &quadobjvar, "quadobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, quadobjvar) );

         if ( pipinput->objsense == SCIP_OBJSENSE_MINIMIZE )
         {
            lhs = -SCIPinfinity(scip);
            rhs = -constant;
         }
         else
         {
            lhs = -constant;
            rhs = SCIPinfinity(scip);
         }

         SCIP_CALL( SCIPcreateConsQuadratic(scip, &quadobjcons, "quadobj", nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoefs, lhs, rhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

         SCIP_CALL( SCIPaddLinearVarQuadratic(scip, quadobjcons, quadobjvar, -1.0) );

         SCIP_CALL( SCIPaddCons(scip, quadobjcons) );
         SCIPdebugMsg(scip, "(line %d) added constraint <%s> to represent quadratic objective: ", pipinput->linenumber, SCIPconsGetName(quadobjcons));
         SCIPdebugPrintCons(scip, quadobjcons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &quadobjcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &quadobjvar) );

         /* free memory */
         SCIPfreeBufferArray(scip, &linvars);
         SCIPfreeBufferArray(scip, &lincoefs);
         SCIPfreeBufferArray(scip, &quadvars1);
         SCIPfreeBufferArray(scip, &quadvars2);
         SCIPfreeBufferArray(scip, &quadcoefs);
      }
      else if( degree > 2 )
      {
         /* insert dummy variable and constraint to represent nonlinear part of objective */

         SCIP_VAR*  nonlinobjvar;
         SCIP_CONS* nonlinobjcons;
         SCIP_Real  minusone;
         SCIP_Real  lhs;
         SCIP_Real  rhs;

         SCIP_CALL( SCIPcreateVar(scip, &nonlinobjvar, "nonlinobjvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, nonlinobjvar) );

         minusone = -1.0;

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

         SCIP_CALL( SCIPcreateConsNonlinear(scip, &nonlinobjcons, "nonlinobj", 1, &nonlinobjvar, &minusone, 1, &exprtree, NULL, lhs, rhs,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
         SCIP_CALL( SCIPexprtreeFree(&exprtree) );

         SCIP_CALL( SCIPaddCons(scip, nonlinobjcons) );
         SCIPdebugMsg(scip, "(line %d) added constraint <%s> to represent nonlinear objective: ", pipinput->linenumber, SCIPconsGetName(nonlinobjcons));
         SCIPdebugPrintCons(scip, nonlinobjcons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &nonlinobjcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &nonlinobjvar) );
      }
   }

   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** reads the constraints section
 */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   PIPINPUT*             pipinput            /**< PIP reading data */
   )
{
   char name[PIP_MAX_LINELEN];
   SCIP_CONS* cons;
   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;
   int degree;

   SCIP_Real constant;

   int nlinvars;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;

   int nquadcoefs;
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;

   PIPSENSE sense;
   SCIP_RETCODE retcode = SCIP_OKAY;
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
   int nmonomials;

   assert(pipinput != NULL);

   /* read polynomial */
   SCIP_CALL( readPolynomial(scip, pipinput, name, &exprtree, &degree, &newsection) );
   if ( hasError(pipinput) )
      goto TERMINATE;
   if ( newsection )
   {
      if ( exprtree != NULL )
         syntaxError(scip, pipinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if ( !getNextToken(scip, pipinput) || !isSense(pipinput, &sense) )
   {
      syntaxError(scip, pipinput, "expected constraint sense '<=', '=', or '>='");
      goto TERMINATE;
   }

   /* read the right hand side */
   sidesign = +1;
   if ( !getNextToken(scip, pipinput) )
   {
      syntaxError(scip, pipinput, "missing right hand side");
      goto TERMINATE;
   }
   if ( isSign(pipinput, &sidesign) )
   {
      if( !getNextToken(scip, pipinput) )
      {
         syntaxError(scip, pipinput, "missing value of right hand side");
         goto TERMINATE;
      }
   }
   if ( !isValue(scip, pipinput, &sidevalue) )
   {
      syntaxError(scip, pipinput, "expected value as right hand side");
      goto TERMINATE;
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

   if( degree > 2 )
   {
      /* assign the left and right hand side, depending on the constraint sense */
      switch ( sense )
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

      SCIP_CALL_TERMINATE( retcode, SCIPcreateConsNonlinear(scip, &cons, name, 0, NULL, NULL, 1, &exprtree, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE), TERMINATE );
   }
   else
   {
      expr = SCIPexprtreeGetRoot(exprtree);
      assert(expr != NULL);
      assert(SCIPexprGetOperator(expr) == SCIP_EXPR_POLYNOMIAL);
      nmonomials = SCIPexprGetNMonomials(expr);

      SCIP_CALL( SCIPallocBufferArray(scip, &linvars,   nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs,  nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, nmonomials) );
      SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, nmonomials) );

      getLinearAndQuadraticCoefs(scip, exprtree, &constant, &nlinvars, linvars, lincoefs, &nquadcoefs, quadvars1, quadvars2, quadcoefs);

      /* assign the left and right hand side, depending on the constraint sense */
      switch( sense )
      {
      case PIP_SENSE_GE:
         lhs = sidevalue - constant;
         rhs = SCIPinfinity(scip);
         break;
      case PIP_SENSE_LE:
         lhs = -SCIPinfinity(scip);
         rhs = sidevalue - constant;
         break;
      case PIP_SENSE_EQ:
         lhs = sidevalue - constant;
         rhs = sidevalue - constant;
         break;
      case PIP_SENSE_NOTHING:
      default:
         SCIPerrorMessage("invalid constraint sense <%d>\n", sense);
         return SCIP_INVALIDDATA;
      }

      if( nquadcoefs == 0 )
      {
         retcode = SCIPcreateConsLinear(scip, &cons, name, nlinvars, linvars, lincoefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE);
      }
      else
      {
         retcode = SCIPcreateConsQuadratic(scip, &cons, name, nlinvars, linvars, lincoefs,
            nquadcoefs, quadvars1, quadvars2, quadcoefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable);
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &linvars);
      SCIPfreeBufferArray(scip, &lincoefs);
      SCIPfreeBufferArray(scip, &quadvars1);
      SCIPfreeBufferArray(scip, &quadvars2);
      SCIPfreeBufferArray(scip, &quadcoefs);
   }

   if( retcode == SCIP_OKAY )
   {
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMsg(scip, "(line %d) created constraint: ", pipinput->linenumber);
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

 TERMINATE:
   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   if( hasError(pipinput) )
      retcode = SCIP_READERROR;

   SCIP_CALL( retcode );

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
   strncat(linebuffer, extension, PIP_MAX_PRINTLEN - strlen(linebuffer));

   (*linecnt) += (int) strlen(extension);

   SCIPdebugMsg(scip, "linebuffer <%s>, length = %lu\n", linebuffer, (unsigned long)strlen(linebuffer));

   if( (*linecnt) > PIP_PRINTLEN )
      endLine(scip, file, linebuffer, linecnt);
}


/* print row in PIP format to file stream */
static
void printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=", "<=", or ">=") */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficient values */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nquadvarterms,      /**< number of quadratic variable terms */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   int                   nbilinterms,        /**< number of bilinear terms */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   int v;
   char linebuffer[PIP_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   SCIP_VAR* var;
   char varname[PIP_MAX_NAMELEN];
   char varname2[PIP_MAX_NAMELEN];
   char consname[PIP_MAX_NAMELEN + 1]; /* an extra character for ':' */
   char buffer[PIP_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strcmp(type, "=") == 0 || strcmp(type, "<=") == 0 || strcmp(type, ">=") == 0 );
   assert( nlinvars == 0 || (linvars != NULL && linvals != NULL) );
   assert( nquadvarterms == 0 || quadvarterms != NULL );

   /* if there is a bilinear term, then there need to be at least two quadratic variables */
   assert( nbilinterms == 0 || (bilinterms != NULL && nquadvarterms >= 2) );

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
   if( nquadvarterms > 0 )
   {
      /* print linear coefficients of quadratic variables */
      for( v = 0; v < nquadvarterms; ++v )
      {
         if( quadvarterms[v].lincoef == 0.0 )
            continue;

         /* we start a new line; therefore we tab this line */
         if (linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(quadvarterms[v].var));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", quadvarterms[v].lincoef, varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      /* start quadratic part */

      /* print square terms */
      for( v = 0; v < nquadvarterms; ++v )
      {
         if( quadvarterms[v].sqrcoef == 0.0 )
            continue;

         /* we start a new line; therefore we tab this line */
         if (linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(quadvarterms[v].var));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s^2", quadvarterms[v].sqrcoef, varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
      }

      /* print bilinear terms */
      for( v = 0; v < nbilinterms; ++v )
      {
         /* we start a new line; therefore we tab this line */
         if (linecnt == 0 )
            appendLine(scip, file, linebuffer, &linecnt, " ");

         (void) SCIPsnprintf(varname,  PIP_MAX_NAMELEN, "%s", SCIPvarGetName(bilinterms[v].var1));
         (void) SCIPsnprintf(varname2, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(bilinterms[v].var2));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s * %s", bilinterms[v].coef, varname, varname2);

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
}


/* print row in PIP format to file stream */
static
void printRowNl(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           rowname,            /**< row name */
   const char*           rownameextension,   /**< row name extension */
   const char*           type,               /**< row type ("=", "<=", or ">=") */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficient values */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            exprtreecoefs,      /**< coefficients of expression trees */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   int v;
   int c;
   int e;
   char linebuffer[PIP_MAX_PRINTLEN] = { '\0' };
   int linecnt;

   SCIP_VAR* var;
   char varname[PIP_MAX_NAMELEN];
   char varname2[PIP_MAX_NAMELEN];
   char consname[PIP_MAX_NAMELEN + 1]; /* an extra character for ':' */
   char buffer[PIP_MAX_PRINTLEN];

   assert( scip != NULL );
   assert( strcmp(type, "=") == 0 || strcmp(type, "<=") == 0 || strcmp(type, ">=") == 0 );
   assert( nlinvars == 0 || (linvars != NULL && linvals != NULL) );
   assert( nexprtrees == 0 || exprtrees != NULL );
   assert( nexprtrees == 0 || exprtreecoefs != NULL );

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
      var = linvars[v];
      assert( var != NULL );

      /* we start a new line; therefore we tab this line */
      if ( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");

      (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(var));
      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", linvals[v], varname);

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print nonlinear part */
   for( e = 0; e < nexprtrees; ++e )
   {
      SCIP_VAR** vars;
      SCIP_EXPR* expr;
      SCIP_EXPR** children;
      int nchildren;

      vars = SCIPexprtreeGetVars(exprtrees[e]);
      expr = SCIPexprtreeGetRoot(exprtrees[e]);
      children = SCIPexprGetChildren(expr);
      nchildren = SCIPexprGetNChildren(expr);
      assert(nchildren == 0 || children != NULL);

      /* we start a new line; therefore we tab this line */
      if( linecnt == 0 )
         appendLine(scip, file, linebuffer, &linecnt, " ");

      /* assert that all children of expr correspond to variables */
#ifndef NDEBUG
      for( c = 0; c < nchildren; ++c )
      {
         assert(SCIPexprGetOperator(children[c]) == SCIP_EXPR_VARIDX);
         assert(SCIPexprGetOpIndex(children[c]) >= 0);
         assert(SCIPexprGetOpIndex(children[c]) < SCIPexprtreeGetNVars(exprtrees[e]));
      }
#endif

      switch( SCIPexprGetOperator(expr) )
      {
      case SCIP_EXPR_CONST:
      {
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g", exprtreecoefs[e] * SCIPexprGetOpReal(expr));
         appendLine(scip, file, linebuffer, &linecnt, buffer);

         break;
      }

      case SCIP_EXPR_VARIDX:
      {
         assert(SCIPexprGetOpIndex(expr) >= 0);
         assert(SCIPexprGetOpIndex(expr) < SCIPexprtreeGetNVars(exprtrees[e]));
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(expr)]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", exprtreecoefs[e], varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_PLUS:
      {
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[0])]));
         (void) SCIPsnprintf(varname2, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[1])]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s %+.15g %s", exprtreecoefs[e], varname, exprtreecoefs[e], varname2);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_MINUS:
      {
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[0])]));
         (void) SCIPsnprintf(varname2, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[1])]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s %+.15g %s", exprtreecoefs[e], varname, -exprtreecoefs[e], varname2);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_MUL:
      {
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[0])]));
         (void) SCIPsnprintf(varname2, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[1])]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s %s", exprtreecoefs[e], varname, varname2);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_SQUARE:
      {
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[0])]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s^2", exprtreecoefs[e], varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_SQRT:
      {
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[0])]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s^0.5", exprtreecoefs[e], varname);

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_INTPOWER:
      {
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[0])]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s^%d", exprtreecoefs[e], varname, SCIPexprGetIntPowerExponent(expr));

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_REALPOWER:
      {
         (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[0])]));
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s^%.15g", exprtreecoefs[e], varname, SCIPexprGetRealPowerExponent(expr));

         appendLine(scip, file, linebuffer, &linecnt, buffer);
         break;
      }

      case SCIP_EXPR_SUM:
      {
         for( c = 0; c < nchildren; ++c )
         {
            /* we start a new line; therefore we tab this line */
            if( linecnt == 0 )
               appendLine(scip, file, linebuffer, &linecnt, " ");

            (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[c])]));
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", exprtreecoefs[e], varname);

            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }

         break;
      }

      case SCIP_EXPR_PRODUCT:
      {
         (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g", exprtreecoefs[e]);
         appendLine(scip, file, linebuffer, &linecnt, buffer);

         for( c = 0; c < nchildren; ++c )
         {
            /* we start a new line; therefore we tab this line */
            if( linecnt == 0 )
               appendLine(scip, file, linebuffer, &linecnt, " ");

            (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, " %s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[c])]));
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, varname);

            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }

         break;
      }

      case SCIP_EXPR_LINEAR:
      {
         if( SCIPexprGetLinearConstant(expr) != 0.0 )
         {
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g ", exprtreecoefs[e] * SCIPexprGetLinearConstant(expr));
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }

         for( c = 0; c < nchildren; ++c )
         {
            /* we start a new line; therefore we tab this line */
            if( linecnt == 0 )
               appendLine(scip, file, linebuffer, &linecnt, " ");

            (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[c])]));
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", exprtreecoefs[e] * SCIPexprGetLinearCoefs(expr)[c], varname);

            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }

         break;
      }

      case SCIP_EXPR_QUADRATIC:
      {
         int q;

         if( SCIPexprGetQuadConstant(expr) != 0.0 )
         {
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g ", exprtreecoefs[e] * SCIPexprGetQuadConstant(expr));
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }

         if( SCIPexprGetQuadLinearCoefs(expr) != NULL )
         {
            for( c = 0; c < nchildren; ++c )
            {
               if( SCIPexprGetQuadLinearCoefs(expr)[c] == 0.0 )
                  continue;

               /* we start a new line; therefore we tab this line */
               if( linecnt == 0 )
                  appendLine(scip, file, linebuffer, &linecnt, " ");

               (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[c])]));
               (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", exprtreecoefs[e] * SCIPexprGetQuadLinearCoefs(expr)[c], varname);

               appendLine(scip, file, linebuffer, &linecnt, buffer);
            }
         }

         for( q = 0; q < SCIPexprGetNQuadElements(expr); ++q )
         {
            /* we start a new line; therefore we tab this line */
            if( linecnt == 0 )
               appendLine(scip, file, linebuffer, &linecnt, " ");

            (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[SCIPexprGetQuadElements(expr)[q].idx1])]));

            if( SCIPexprGetQuadElements(expr)[q].idx1 == SCIPexprGetQuadElements(expr)[q].idx2 )
            {
               /* square term */
               (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s^2", exprtreecoefs[e] * SCIPexprGetQuadElements(expr)[q].coef, varname);
            }
            else
            {
               /* bilinear term */
               (void) SCIPsnprintf(varname2, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[SCIPexprGetQuadElements(expr)[q].idx2])]));
               (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s %s", exprtreecoefs[e] * SCIPexprGetQuadElements(expr)[q].coef, varname, varname2);
            }

            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }

         break;
      }

      case SCIP_EXPR_POLYNOMIAL:
      {
         SCIP_EXPRDATA_MONOMIAL* monomial;
         int m;
         int f;

         if( SCIPexprGetPolynomialConstant(expr) != 0.0 )
         {
            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g ", exprtreecoefs[e] * SCIPexprGetPolynomialConstant(expr));
            appendLine(scip, file, linebuffer, &linecnt, buffer);
         }

         for( m = 0; m < SCIPexprGetNMonomials(expr); ++m )
         {
            monomial = SCIPexprGetMonomials(expr)[m];
            assert(monomial != NULL);

            (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g", exprtreecoefs[e] * SCIPexprGetMonomialCoef(monomial));
            appendLine(scip, file, linebuffer, &linecnt, buffer);

            for( f = 0; f < SCIPexprGetMonomialNFactors(monomial); ++f )
            {
               /* we start a new line; therefore we tab this line */
               if( linecnt == 0 )
                  appendLine(scip, file, linebuffer, &linecnt, " ");

               (void) SCIPsnprintf(varname, PIP_MAX_NAMELEN, "%s", SCIPvarGetName(vars[SCIPexprGetOpIndex(children[SCIPexprGetMonomialChildIndices(monomial)[f]])]));
               if( SCIPexprGetMonomialExponents(monomial)[f] != 1.0 )
                  (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %s^%.15g", varname, SCIPexprGetMonomialExponents(monomial)[f]);
               else
                  (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %s", varname);
               appendLine(scip, file, linebuffer, &linecnt, buffer);
            }
         }

         break;
      }

      default:
      {
         /* this should have been caught in SCIPwritePip before */
         SCIPerrorMessage("unsupported operator <%s> in writing of polynomial nonlinear constraint\n", SCIPexpropGetName(SCIPexprGetOperator(expr)));
         return;
      } /*lint !e788*/
      }  /*lint !e788*/
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
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nquadvarterms,      /**< number of quadratic variable terms */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   int                   nbilinterms,        /**< number of bilinear terms */
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
   assert( nquadvarterms == 0 || quadvarterms != NULL );
   assert( nbilinterms == 0 || bilinterms != NULL );

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
      printRow(scip, file, rowname, "", "=", activevars, activevals, nactivevars,
         quadvarterms, nquadvarterms, bilinterms, nbilinterms,
         rhs - activeconstant);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRow(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", ">=",
            activevars, activevals, nactivevars,
            quadvarterms, nquadvarterms, bilinterms, nbilinterms,
            lhs - activeconstant);
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         printRow(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "<=",
            activevars, activevals, nactivevars,
            quadvarterms, nquadvarterms, bilinterms, nbilinterms,
            rhs - activeconstant);
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
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            linvals,            /**< array of linear coefficients values (or NULL if all linear coefficient values are 1) */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            exprtreecoefs,      /**< coefficients of expression trees */
   int                   nexprtrees,         /**< number of expression trees */
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
   assert( nexprtrees == 0 || exprtrees != NULL );
   assert( nexprtrees == 0 || exprtreecoefs != NULL );

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
      printRowNl(scip, file, rowname, "", "=", activevars, activevals, nactivevars,
         exprtrees, exprtreecoefs, nexprtrees,
         rhs - activeconstant);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         printRowNl(scip, file, rowname, SCIPisInfinity(scip, rhs) ? "" : "_lhs", ">=",
            activevars, activevals, nactivevars,
            exprtrees, exprtreecoefs, nexprtrees,
            lhs - activeconstant);
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         printRowNl(scip, file, rowname, SCIPisInfinity(scip, -lhs) ? "" : "_rhs", "<=",
            activevars, activevals, nactivevars,
            exprtrees, exprtreecoefs, nexprtrees,
            rhs - activeconstant);
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

/** check whether given variables are aggregated and put them into an array without duplication */
static
SCIP_RETCODE collectAggregatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
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
      printRow(scip, file, consname, "", "=", activevars, activevals, nactivevars, NULL, 0, NULL, 0, - activeconstant);
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

/** method check if the variable names are not longer than PIP_MAX_NAMELEN */
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

   /* check if the variable names are not to long */
   for( v = 0; v < nvars; ++v )
   {
      if( strlen(SCIPvarGetName(vars[v])) > PIP_MAX_NAMELEN )  /*lint !e613*/
      {
         SCIPwarningMessage(scip, "there is a variable name which has to be cut down to %d characters; LP might be corrupted\n", 
            PIP_MAX_NAMELEN - 1);
         return;
      }
   }
}

/** method check if the constraint names are not longer than PIP_MAX_NAMELEN */
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
   assert( conss != NULL );

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL );

      /* in case the transformed is written only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_Real lhs = SCIPgetLhsLinear(scip, cons);
         SCIP_Real rhs = SCIPgetLhsLinear(scip, cons);

         if( (SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > PIP_MAX_NAMELEN)
            || ( !SCIPisEQ(scip, lhs, rhs) && strlen(SCIPconsGetName(conss[c])) > PIP_MAX_NAMELEN -  4) )
         {
            SCIPwarningMessage(scip, "there is a constraint name which has to be cut down to %d characters;\n",
               PIP_MAX_NAMELEN  - 1);
            return;
         }
      }
      else if( strlen(SCIPconsGetName(conss[c])) > PIP_MAX_NAMELEN )
      {
         SCIPwarningMessage(scip, "there is a constraint name which has to be cut down to %d characters;\n",
            PIP_MAX_NAMELEN  - 1);
         return;
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
   int e;

   int linecnt;
   char linebuffer[PIP_MAX_PRINTLEN];

   char varname[PIP_MAX_NAMELEN];
   char buffer[PIP_MAX_PRINTLEN];

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;
   SCIP_CONS** consQuadratic;
   int nConsQuadratic;
   SCIP_CONS** consNonlinear;
   int nConsNonlinear;
   SCIP_CONS** consAbspower;
   int nConsAbspower;
   SCIP_CONS** consAnd;
   int nConsAnd;
   SCIP_CONS** consBivariate;
   int nConsBivariate;
   char consname[PIP_MAX_NAMELEN];

   SCIP_VAR** aggregatedVars;
   int nAggregatedVars;
   SCIP_HASHTABLE* varAggregated;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int nconsvars;

   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;

   SCIP_EXPRTREE* exprtree;
   SCIP_EXPR* expr;

   assert( scip != NULL );

   nAggregatedVars = 0;
   nConsQuadratic = 0;
   nConsNonlinear = 0;
   nConsAbspower = 0;
   nConsAnd = 0;
   nConsBivariate = 0;

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
   SCIPinfoMessage(scip, file, "\\   Obj. scale       : %.15g\n", objscale);
   SCIPinfoMessage(scip, file, "\\   Obj. offset      : %.15g\n", objoffset);

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
      (void) SCIPsnprintf(buffer, PIP_MAX_PRINTLEN, " %+.15g %s", SCIPvarGetObj(var), varname );

      appendLine(scip, file, linebuffer, &linecnt, buffer);
   }

   endLine(scip, file, linebuffer, &linecnt);

   /* print "Subject to" section */
   SCIPinfoMessage(scip, file, "Subject to\n");

   /* collect quadratic, nonlinear, absolute power, and, and bivariate constraints in arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &consQuadratic, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consNonlinear, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consAbspower, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consAnd, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consBivariate, nconss) );

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
               NULL, 0, NULL, 0, SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed) );
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         switch( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 0, NULL, 0, 1.0, 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 0, NULL, 0, -SCIPinfinity(scip), 1.0, transformed) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printQuadraticCons(scip, file, consname,
                  consvars, NULL, nconsvars, NULL, 0, NULL, 0, 1.0, SCIPinfinity(scip), transformed) );
            break;
         }
      }
      else if ( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
               NULL, 0, NULL, 0, 1.0, SCIPinfinity(scip), transformed) );
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
               NULL, 0, NULL, 0, -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed) );

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

         SCIP_CALL( printQuadraticCons(scip, file, consname, consvars, consvals, 2, NULL, 0, NULL, 0, 
               SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "quadratic") == 0 )
      {
         SCIP_CALL( printQuadraticCons(scip, file, consname,
               SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
               SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
               SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
               SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetLhsQuadratic(scip, cons),
               SCIPgetRhsQuadratic(scip, cons), transformed) );

         consQuadratic[nConsQuadratic++] = cons;
      }
      else if( strcmp(conshdlrname, "nonlinear") == 0 )
      {
         SCIP_Bool ispolynomial;
         int nexprtrees = SCIPgetNExprtreesNonlinear(scip, cons);

         /* check whether expressions are polynomials (assumed simplified exprtrees) */
         ispolynomial = TRUE;
         for( e = 0; e < nexprtrees && ispolynomial; ++e )
         {
            exprtree = SCIPgetExprtreesNonlinear(scip, cons)[e];
            expr = SCIPexprtreeGetRoot(exprtree);
            assert(expr != NULL);

            /* check if operator is something polynomial */
            switch( SCIPexprGetOperator(expr) )
            {
            case SCIP_EXPR_CONST:
            case SCIP_EXPR_VARIDX:
            case SCIP_EXPR_PLUS:
            case SCIP_EXPR_MINUS:
            case SCIP_EXPR_MUL:
            case SCIP_EXPR_SQUARE:
            case SCIP_EXPR_SQRT:
            case SCIP_EXPR_SUM:
            case SCIP_EXPR_PRODUCT:
            case SCIP_EXPR_LINEAR:
            case SCIP_EXPR_QUADRATIC:
               break;

            case SCIP_EXPR_INTPOWER:
            {
               if( SCIPexprGetIntPowerExponent(expr) < 0 )
               {
                  SCIPwarningMessage(scip, "negative exponent %d in intpower in %dth expression tree of constraint <%s> cannot be written in pip format\n", SCIPexprGetIntPowerExponent(expr), e, SCIPconsGetName(cons));
                  ispolynomial = FALSE;
               }

               break;
            }

            case SCIP_EXPR_REALPOWER:
            {
               if( SCIPexprGetRealPowerExponent(expr) < 0.0 )
               {
                  SCIPwarningMessage(scip, "negative exponent %g in realpower in %dth expression tree of constraint <%s> cannot be written in pip format\n", SCIPexprGetRealPowerExponent(expr), e, SCIPconsGetName(cons));
                  ispolynomial = FALSE;
               }

               break;
            }

            case SCIP_EXPR_POLYNOMIAL:
            {
               SCIP_EXPRDATA_MONOMIAL* monomial;
               int m;
               int f;

               for( m = 0; m < SCIPexprGetNMonomials(expr) && ispolynomial; ++m )
               {
                  monomial = SCIPexprGetMonomials(expr)[m];
                  for( f = 0; f < SCIPexprGetMonomialNFactors(monomial); ++f )
                  {
                     if( SCIPexprGetMonomialExponents(monomial)[f] < 0.0 )
                     {
                        SCIPwarningMessage(scip, "negative exponent %g in polynomial in %dth expression tree of constraint <%s> cannot be written in pip format\n", SCIPexprGetMonomialExponents(monomial)[f], e, SCIPconsGetName(cons));
                        ispolynomial = FALSE;
                        break;
                     }
                  }
               }

               break;
            }

            default:
               SCIPwarningMessage(scip, "expression operand <%s> in %dth expression tree of constraint <%s> cannot be written in pip format\n", SCIPexpropGetName(SCIPexprGetOperator(expr)), e, SCIPconsGetName(cons));
               ispolynomial = FALSE;
               break;
            } /*lint !e788*/

            /* check if all children of root expression correspond to variables */
            for( v = 0; v < SCIPexprGetNChildren(expr) && ispolynomial; ++v )
            {
               if( SCIPexprGetOperator(SCIPexprGetChildren(expr)[v]) != SCIP_EXPR_VARIDX )
               {
                  SCIPwarningMessage(scip, "%dth expression tree of constraint <%s> is not simplified, cannot write in pip format\n", e, SCIPconsGetName(cons));
                  ispolynomial = FALSE;
               }
            }
         }

         if( ispolynomial )
         {
            SCIP_CALL( printNonlinearCons(scip, file, consname,
                  SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
                  SCIPgetNLinearVarsNonlinear(scip, cons), SCIPgetExprtreesNonlinear(scip, cons),
                  SCIPgetExprtreeCoefsNonlinear(scip, cons), SCIPgetNExprtreesNonlinear(scip, cons),
                  SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons), transformed) );

            consNonlinear[nConsNonlinear++] = cons;
         }
         else
         {
            SCIPinfoMessage(scip, file, "\\ ");
            SCIP_CALL( SCIPprintCons(scip, cons, file) );
            SCIPinfoMessage(scip, file, ";\n");
         }
      }
      else if( strcmp(conshdlrname, "abspower") == 0 )
      {
         SCIP_VAR* x;
         SCIP_Real xoffset;
         SCIP_Real exponent;
         SCIP_Real treecoef;

         expr = NULL;
         treecoef = 1.0;

         x = SCIPgetNonlinearVarAbspower(scip, cons);
         xoffset = SCIPgetOffsetAbspower(scip, cons);
         exponent = SCIPgetExponentAbspower(scip, cons);

         /* see if we formulate signpower(x+offset,exponent) as usual polynomial */
         if( !SCIPisZero(scip, xoffset) )
         {
            SCIPwarningMessage(scip, "nonzero offset for nonlinear variable in constraint <%s>, cannot write in pip format\n", SCIPconsGetName(cons));
         }
         if( SCIPisIntegral(scip, exponent) && ((int)SCIPround(scip, exponent) % 2 == 1) )
         {
            /* exponent is odd integer, so signpower(x,exponent) = x^exponent */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_INTPOWER, expr, (int)SCIPround(scip, exponent)) );
         }
         else if( SCIPisIntegral(scip, exponent) && ((int)SCIPround(scip, exponent) % 2 == 0) && !SCIPisPositive(scip, SCIPvarGetUbGlobal(x)) )
         {
            /* exponent is even integer and x is negative, so signpower(x,exponent) = -x^exponent */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_INTPOWER, expr, (int)SCIPround(scip, exponent)) );
            treecoef = -1.0;
         }
         else if( !SCIPisNegative(scip, SCIPvarGetLbGlobal(x)) )
         {
            /* x is positive, so signpower(x,exponent) = x^exponent */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_REALPOWER, expr, exponent) );
         }
         else
         {
            SCIPwarningMessage(scip, "cannot formulate signpower(<%s>, %g) in constraint <%s> as polynomial, cannot write in pip format\n", SCIPvarGetName(x), exponent, SCIPconsGetName(cons));
         }

         if( expr != NULL )
         {
            SCIP_VAR* z;
            SCIP_Real zcoef;

            SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 1, 0, NULL) );
            SCIP_CALL( SCIPexprtreeSetVars(exprtree, 1, &x) );

            z = SCIPgetLinearVarAbspower(scip, cons);
            zcoef = SCIPgetCoefLinearAbspower(scip, cons);

            SCIP_CALL( printNonlinearCons(scip, file, consname,
                  &z, &zcoef, 1, &exprtree, &treecoef, 1,
                  SCIPgetLhsAbspower(scip, cons), SCIPgetRhsAbspower(scip, cons), transformed) );

            SCIP_CALL( SCIPexprtreeFree(&exprtree) );

            consAbspower[nConsAbspower++] = cons;
         }
         else
         {
            SCIPinfoMessage(scip, file, "\\ ");
            SCIP_CALL( SCIPprintCons(scip, cons, file) );
            SCIPinfoMessage(scip, file, ";\n");
         }
      }
      else if( strcmp(conshdlrname, "bivariate") == 0 )
      {
         SCIP_Bool ispolynomial;

         /* check whether expression is polynomial (simplified exprtree assumed) */
         ispolynomial = TRUE;
         exprtree = SCIPgetExprtreeBivariate(scip, cons);
         expr = SCIPexprtreeGetRoot(exprtree);
         assert(expr != NULL);

         /* check if operator is something polynomial */
         switch( SCIPexprGetOperator(expr) )
         {
         case SCIP_EXPR_CONST:
         case SCIP_EXPR_VARIDX:
         case SCIP_EXPR_PLUS:
         case SCIP_EXPR_MINUS:
         case SCIP_EXPR_MUL:
         case SCIP_EXPR_SQUARE:
         case SCIP_EXPR_SQRT:
         case SCIP_EXPR_SUM:
         case SCIP_EXPR_PRODUCT:
         case SCIP_EXPR_LINEAR:
         case SCIP_EXPR_QUADRATIC:
            break;

         case SCIP_EXPR_INTPOWER:
         {
            if( SCIPexprGetIntPowerExponent(expr) < 0 )
            {
               SCIPwarningMessage(scip, "negative exponent %d in intpower of constraint <%s> cannot be written in pip format\n", SCIPexprGetIntPowerExponent(expr), SCIPconsGetName(cons));
               ispolynomial = FALSE;
            }

            break;
         }

         case SCIP_EXPR_REALPOWER:
         {
            if( SCIPexprGetRealPowerExponent(expr) < 0.0 )
            {
               SCIPwarningMessage(scip, "negative exponent %g in realpower of constraint <%s> cannot be written in pip format\n", SCIPexprGetRealPowerExponent(expr), SCIPconsGetName(cons));
               ispolynomial = FALSE;
            }

            break;
         }

         case SCIP_EXPR_POLYNOMIAL:
         {
            SCIP_EXPRDATA_MONOMIAL* monomial;
            int m;
            int f;

            for( m = 0; m < SCIPexprGetNMonomials(expr) && ispolynomial; ++m )
            {
               monomial = SCIPexprGetMonomials(expr)[m];
               for( f = 0; f < SCIPexprGetMonomialNFactors(monomial); ++f )
               {
                  if( SCIPexprGetMonomialExponents(monomial)[f] < 0.0 )
                  {
                     SCIPwarningMessage(scip, "negative exponent %g in polynomial of constraint <%s> cannot be written in pip format\n", SCIPexprGetMonomialExponents(monomial)[f], SCIPconsGetName(cons));
                     ispolynomial = FALSE;
                     break;
                  }
               }
            }

            break;
         }

         default:
            SCIPwarningMessage(scip, "expression operand <%s> in constraint <%s> cannot be written in pip format\n", SCIPexpropGetName(SCIPexprGetOperator(expr)), SCIPconsGetName(cons));
            ispolynomial = FALSE;
            break;
         } /*lint !e788*/

         if( ispolynomial )
         {
            /* check if all children of root expression correspond to variables */
            for( v = 0; v < SCIPexprGetNChildren(expr); ++v )
            {
               if( SCIPexprGetOperator(SCIPexprGetChildren(expr)[v]) != SCIP_EXPR_VARIDX )
               {
                  SCIPwarningMessage(scip, "expression tree of constraint <%s> is not simplified, cannot write in pip format\n", SCIPconsGetName(cons));
                  ispolynomial = FALSE;
                  break;
               }
            }
         }

         if( ispolynomial )
         {
            SCIP_VAR* z;
            SCIP_Real zcoef;
            SCIP_Real one;

            z = SCIPgetLinearVarBivariate(scip, cons);
            zcoef = SCIPgetLinearCoefBivariate(scip, cons);

            one = 1.0;
            SCIP_CALL( printNonlinearCons(scip, file, consname,
                  &z, &zcoef, z == NULL ? 0 : 1, &exprtree, &one, 1,
                  SCIPgetLhsBivariate(scip, cons), SCIPgetRhsBivariate(scip, cons), transformed) );

            consBivariate[nConsBivariate++] = cons;
         }
         else
         {
            SCIPinfoMessage(scip, file, "\\ ");
            SCIP_CALL( SCIPprintCons(scip, cons, file) );
            SCIPinfoMessage(scip, file, ";\n");
         }
      }
      else if( strcmp(conshdlrname, "and") == 0 )
      {
         SCIP_EXPR** children;
         SCIP_VAR* resultant;
         SCIP_Real minusone;
         SCIP_Real one;

         /* create expression for product of binaries */
         SCIP_CALL( SCIPallocBufferArray(scip, &children, SCIPgetNVarsAnd(scip, cons)) );
         for( v = 0; v < SCIPgetNVarsAnd(scip, cons); ++v )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[v], SCIP_EXPR_VARIDX, v) );
         }
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PRODUCT, SCIPgetNVarsAnd(scip, cons), children) );
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, SCIPgetNVarsAnd(scip, cons), 0, NULL) );
         SCIP_CALL( SCIPexprtreeSetVars(exprtree, SCIPgetNVarsAnd(scip, cons), SCIPgetVarsAnd(scip, cons)) );

         resultant = SCIPgetResultantAnd(scip, cons);
         minusone = -1.0;

         one = 1.0;
         SCIP_CALL( printNonlinearCons(scip, file, consname, &resultant, &minusone, 1, &exprtree, &one, 1, 0.0, 0.0, transformed) );

         SCIP_CALL( SCIPexprtreeFree(&exprtree) );
         SCIPfreeBufferArray(scip, &children);

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

   /* check for aggregated variables in quadratic parts of quadratic constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsQuadratic; ++c)
   {
      cons = consQuadratic[c];
      for( v = 0; v < SCIPgetNQuadVarTermsQuadratic(scip, cons); ++v )
      {
         SCIP_CALL( collectAggregatedVars(scip, 1, &SCIPgetQuadVarTermsQuadratic(scip, cons)[v].var, 
               &nAggregatedVars, &aggregatedVars, &varAggregated) );
      }         
   }

   /* check for aggregated variables in expression trees of nonlinear constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsNonlinear; ++c)
   {
      cons = consNonlinear[c];
      for( e = 0; e < SCIPgetNExprtreesNonlinear(scip, cons); ++e )
      {
         exprtree = SCIPgetExprtreesNonlinear(scip, cons)[e];
         assert(exprtree != NULL);

         for( v = 0; v < SCIPexprtreeGetNVars(exprtree); ++v )
         {
            SCIP_CALL( collectAggregatedVars(scip, 1, &SCIPexprtreeGetVars(exprtree)[v],
                  &nAggregatedVars, &aggregatedVars, &varAggregated) );
         }
      }
   }

   /* check for aggregated variables in absolute power constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsAbspower; ++c)
   {
      SCIP_VAR* spvars[2];

      cons = consAbspower[c];

      spvars[0] = SCIPgetNonlinearVarAbspower(scip, cons);
      spvars[1] = SCIPgetLinearVarAbspower(scip, cons);
      SCIP_CALL( collectAggregatedVars(scip, 2, spvars, &nAggregatedVars, &aggregatedVars, &varAggregated) );
   }

   /* check for aggregated variables in and constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsAnd; ++c)
   {
      SCIP_VAR* resultant;

      cons = consAnd[c];

      SCIP_CALL( collectAggregatedVars(scip, SCIPgetNVarsAnd(scip, cons), SCIPgetVarsAnd(scip, cons), &nAggregatedVars, &aggregatedVars, &varAggregated) );

      resultant = SCIPgetResultantAnd(scip, cons);
      SCIP_CALL( collectAggregatedVars(scip, 1, &resultant, &nAggregatedVars, &aggregatedVars, &varAggregated) );
   }

   /* check for aggregated variables in bivariate constraints and output aggregations as linear constraints */
   for (c = 0; c < nConsBivariate; ++c)
   {
      SCIP_VAR* z;

      cons = consBivariate[c];

      assert(SCIPexprtreeGetNVars(SCIPgetExprtreeBivariate(scip, cons)) == 2);
      SCIP_CALL( collectAggregatedVars(scip, 2, SCIPexprtreeGetVars(SCIPgetExprtreeBivariate(scip, cons)), &nAggregatedVars, &aggregatedVars, &varAggregated) );

      z = SCIPgetLinearVarBivariate(scip, cons);
      if( z != NULL )
      {
         SCIP_CALL( collectAggregatedVars(scip, 1, &z, &nAggregatedVars, &aggregatedVars, &varAggregated) );
      }
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
   SCIPfreeBufferArray(scip, &consQuadratic);
   SCIPfreeBufferArray(scip, &consNonlinear);
   SCIPfreeBufferArray(scip, &consAbspower);
   SCIPfreeBufferArray(scip, &consAnd);
   SCIPfreeBufferArray(scip, &consBivariate);

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
