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

/**@file   reader_opb.c
 * @brief  pseudo-Boolean file reader (opb format)
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * This file reader parses the @a opb format and is also used by the @a wbo reader for the @a wbo format. For a
 * detailed description of this format see
 *
 * - http://www.cril.univ-artois.fr/PB07/solver_req.html
 * - http://www.cril.univ-artois.fr/PB10/format.pdf
 *
 * The syntax of the input file format can be described by a simple Backus-Naur
 *  form. \<formula\> is the start symbol of this grammar.
 *
 *  \<formula\>::= \<sequence_of_comments\>
 *               [\<objective\>] | [\<softheader\>]
 *               \<sequence_of_comments_or_constraints\>
 *
 *  \<sequence_of_comments\>::= \<comment\> [\<sequence_of_comments\>]
 *  \<comment\>::= "*" \<any_sequence_of_characters_other_than_EOL\> \<EOL\>
 *  \<sequence_of_comments_or_constraints\>::=\<comment_or_constraint\> [\<sequence_of_comments_or_constraints\>]
 *  \<comment_or_constraint\>::=\<comment\>|\<constraint\>
 *
 *  \<objective\>::= "min:" \<zeroOrMoreSpace\> \<sum\>  ";"
 *  \<constraint\>::= \<sum\> \<relational_operator\> \<zeroOrMoreSpace\> \<integer\> \<zeroOrMoreSpace\> ";"
 *
 *  \<sum\>::= \<weightedterm\> | \<weightedterm\> \<sum\>
 *  \<weightedterm\>::= \<integer\> \<oneOrMoreSpace\> \<term\> \<oneOrMoreSpace\>
 *
 *  \<integer\>::= \<unsigned_integer\> | "+" \<unsigned_integer\> | "-" \<unsigned_integer\>
 *  \<unsigned_integer\>::= \<digit\> | \<digit\>\<unsigned_integer\>
 *
 *  \<relational_operator\>::= "\>=" | "="
 *
 *  \<variablename\>::= "x" \<unsigned_integer\>
 *
 *  \<oneOrMoreSpace\>::= " " [\<oneOrMoreSpace\>]
 *  \<zeroOrMoreSpace\>::= [" " \<zeroOrMoreSpace\>]
 *
 *  For linear pseudo-Boolean instances, \<term\> is defined as
 *
 *  \<term\>::=\<variablename\>
 *
 *  For non-linear instances, \<term\> is defined as
 *
 *  \<term\>::= \<oneOrMoreLiterals\>
 *  \<oneOrMoreLiterals\>::= \<literal\> | \<literal\> \<oneOrMoreSpace\> \<oneOrMoreLiterals\>
 *  \<literal\>::= \<variablename\> | "~"\<variablename\>
 *
 * For wbo-files are the following additional/changed things possible.
 *
 *  \<softheader\>::= "soft:" [\<unsigned integer\>] ";"
 *
 *  \<comment_or_constraint\>::=\<comment\>|\<constraint\>|\<softconstraint\>
 *
 *  \<softconstraint\>::= "[" \<zeroOrMoreSpace\> \<unsigned integer\> \<zeroOrMoreSpace\> "]" \<constraint\>
 *
 */

/* Our parser should also be lax by handling variable names and it's possible to read doubles instead of integer and
 * possible some more :). */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/ /* needed for strncasecmp() */
#endif
#include <ctype.h>

#include "scip/cons_and.h"
#include "scip/cons_indicator.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_pseudoboolean.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/pub_misc.h"
#include "scip/reader_opb.h"
#include "scip/debug.h"

#define READER_NAME             "opbreader"
#define READER_DESC             "file reader for pseudo-Boolean problem in opb format"
#define READER_EXTENSION        "opb"

#define GENCONSNAMES            TRUE  /* remove if no constraint names should be generated */
#define LINEAROBJECTIVE         TRUE  /* will all non-linear parts inside the objective function be linearized or will
                                       * an artificial integer variable be created which will represent the objective
                                       * function
                                       */

#define INDICATORVARNAME        "indicatorvar" /* standard part of name for all indicator variables */
#define INDICATORSLACKVARNAME   "indslack"     /* standard part of name for all indicator slack variables; should be the same in cons_indicator */
#define TOPCOSTCONSNAME         "topcostcons"  /* standard name for artificial topcost constraint in wbo problems */

/*
 * Data structures
 */
#define OPB_MAX_LINELEN        65536  /**< size of the line buffer for reading or writing */
#define OPB_MAX_PUSHEDTOKENS   2
#define OPB_INIT_COEFSSIZE     8192

/** Section in OPB File */
enum OpbExpType
{
   OPB_EXP_NONE,
   OPB_EXP_UNSIGNED,
   OPB_EXP_SIGNED
};
typedef enum OpbExpType OPBEXPTYPE;

enum OpbSense
{
   OPB_SENSE_NOTHING,
   OPB_SENSE_LE,
   OPB_SENSE_GE,
   OPB_SENSE_EQ
};
typedef enum OpbSense OPBSENSE;

/** OPB reading data */
struct OpbInput
{
   SCIP_FILE*            file;
   char                  linebuf[OPB_MAX_LINELEN+1];
   char*                 token;
   char*                 tokenbuf;
   char*                 pushedtokens[OPB_MAX_PUSHEDTOKENS];
   int                   npushedtokens;
   int                   linenumber;
   int                   linepos;
   int                   bufpos;
   SCIP_OBJSENSE         objsense;
   SCIP_Bool             comment;
   SCIP_Bool             endline;
   SCIP_Bool             eof;
   SCIP_Bool             haserror;
   int                   nproblemcoeffs;
   SCIP_Bool             wbo;
   SCIP_Real             topcost;
   int                   nindvars;
#if GENCONSNAMES == TRUE
   int                   consnumber;
#endif
};

typedef struct OpbInput OPBINPUT;

static const char commentchars[] = "*";
/*
 * Local methods (for reading)
 */

/** issues an error message and marks the OPB data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   const char*           msg                 /**< error message */
   )
{
   assert(opbinput != NULL);

   SCIPerrorMessage("Syntax error in line %d: %s found <%s>\n", opbinput->linenumber, msg, opbinput->token);
   if( opbinput->linebuf[strlen(opbinput->linebuf)-1] == '\n' )
   {
      SCIPerrorMessage("  input: %s", opbinput->linebuf);
   }
   else
   {
      SCIPerrorMessage("  input: %s\n", opbinput->linebuf);
   }

   opbinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(opbinput != NULL);

   return opbinput->haserror;
}

/** returns whether the given character is a token delimiter */
static
SCIP_Bool isDelimChar(
   char                  c                   /**< input character */
   )
{
   switch (c)
   {
   case ' ':
   case '\f':
   case '\n':
   case '\r':
   case '\t':
   case '\v':
   case '\0':
      return TRUE;
   default:
      return FALSE;
   }
}

/** returns whether the given character is a single token */
static
SCIP_Bool isTokenChar(
   char                  c                   /**< input character */
   )
{
   switch (c)
   {
   case '-':
   case '+':
   case ':':
   case '<':
   case '>':
   case '=':
   case '[':
   case ']':
   case ';':
      return TRUE;
   default:
      return FALSE;
   }
}

/** returns whether the current character is member of a value string */
static
SCIP_Bool isValueChar(
   char                  c,                  /**< input character */
   char                  nextc,              /**< next input character */
   SCIP_Bool             firstchar,          /**< is the given character the first char of the token? */
   SCIP_Bool*            hasdot,             /**< pointer to update the dot flag */
   OPBEXPTYPE*           exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit((unsigned char)c) )
      return TRUE;
   else if( (*exptype == OPB_EXP_NONE) && !(*hasdot) && (c == '.') )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == OPB_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = OPB_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit((unsigned char)nextc) )
      {
         *exptype = OPB_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == OPB_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = OPB_EXP_UNSIGNED;
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
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   int i;

   assert(opbinput != NULL);

   /* if we previously detected a comment we have to parse the remaining line away if there is something left */
   if( !opbinput->endline && opbinput->comment )
   {
      SCIPdebugMsg(scip, "Throwing rest of comment away.\n");

      do
      {
         opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';
         (void)SCIPfgets(opbinput->linebuf, (int) sizeof(opbinput->linebuf), opbinput->file);
      }
      while( opbinput->linebuf[OPB_MAX_LINELEN-2] != '\0' );

      opbinput->comment = FALSE;
      opbinput->endline = TRUE;
   }

   /* clear the line */
   opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';

   /* set line position */
   if( opbinput->endline )
   {
      opbinput->linepos = 0;
      opbinput->linenumber++;
   }
   else
      opbinput->linepos += OPB_MAX_LINELEN - 2;

   if( SCIPfgets(opbinput->linebuf, (int) sizeof(opbinput->linebuf), opbinput->file) == NULL )
      return FALSE;

   opbinput->bufpos = 0;

   if( opbinput->linebuf[OPB_MAX_LINELEN-2] != '\0' )
   {
      char* last;

      /* buffer is full; erase last token since it might be incomplete */
      opbinput->endline = FALSE;
      last = strrchr(opbinput->linebuf, ' ');

      if( last == NULL )
      {
         SCIPwarningMessage(scip, "we read %d character from the file; these might indicates a corrupted input file!",
            OPB_MAX_LINELEN - 2);
         opbinput->linebuf[OPB_MAX_LINELEN-2] = '\0';
         SCIPdebugMsg(scip, "the buffer might be corrupted\n");
      }
      else
      {
         SCIPfseek(opbinput->file, -(long) strlen(last), SEEK_CUR);
         SCIPdebugMsg(scip, "correct buffer, reread the last %ld characters\n", (long) strlen(last));
         *last = '\0';
      }
   }
   else
   {
      /* found end of line */
      opbinput->endline = TRUE;
   }

   opbinput->linebuf[OPB_MAX_LINELEN-1] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

   opbinput->comment = FALSE;

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(opbinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
         opbinput->comment = TRUE;
         break;
      }
   }

   SCIPdebugMsg(scip, "%s\n", opbinput->linebuf);

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
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   SCIP_Bool hasdot;
   OPBEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(opbinput != NULL);
   assert(opbinput->bufpos < OPB_MAX_LINELEN);

   /* check the token stack */
   if( opbinput->npushedtokens > 0 )
   {
      swapPointers(&opbinput->token, &opbinput->pushedtokens[opbinput->npushedtokens-1]);
      opbinput->npushedtokens--;
      SCIPdebugMsg(scip, "(line %d) read token again: '%s'\n", opbinput->linenumber, opbinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = opbinput->linebuf;
   while( isDelimChar(buf[opbinput->bufpos]) )
   {
      if( buf[opbinput->bufpos] == '\0' )
      {
         if( !getNextLine(scip, opbinput) )
         {
            SCIPdebugMsg(scip, "(line %d) end of file\n", opbinput->linenumber);
            return FALSE;
         }
         assert(opbinput->bufpos == 0);
      }
      else
      {
         opbinput->bufpos++;
         opbinput->linepos++;
      }
   }
   assert(opbinput->bufpos < OPB_MAX_LINELEN);
   assert(!isDelimChar(buf[opbinput->bufpos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = OPB_EXP_NONE;
   if( isValueChar(buf[opbinput->bufpos], buf[opbinput->bufpos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < OPB_MAX_LINELEN);
         assert(!isDelimChar(buf[opbinput->bufpos]));
         opbinput->token[tokenlen] = buf[opbinput->bufpos];
         tokenlen++;
         opbinput->bufpos++;
         opbinput->linepos++;
      }
      while( isValueChar(buf[opbinput->bufpos], buf[opbinput->bufpos+1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < OPB_MAX_LINELEN);
         opbinput->token[tokenlen] = buf[opbinput->bufpos];
         tokenlen++;
         opbinput->bufpos++;
         opbinput->linepos++;
         if( tokenlen == 1 && isTokenChar(opbinput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[opbinput->bufpos]) && !isTokenChar(buf[opbinput->bufpos]) );

      /* if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>',
       * replace the token by the inequality sense
       */
      if( tokenlen >= 1
         && (opbinput->token[tokenlen-1] == '<' || opbinput->token[tokenlen-1] == '>' || opbinput->token[tokenlen-1] == '=')
         && buf[opbinput->bufpos] == '=' )
      {
         opbinput->bufpos++;
         opbinput->linepos++;
      }
      else if( opbinput->token[tokenlen-1] == '=' && (buf[opbinput->bufpos] == '<' || buf[opbinput->bufpos] == '>') )
      {
         opbinput->token[tokenlen-1] = buf[opbinput->bufpos];
         opbinput->bufpos++;
         opbinput->linepos++;
      }
   }
   assert(tokenlen < OPB_MAX_LINELEN);
   opbinput->token[tokenlen] = '\0';

   SCIPdebugMsg(scip, "(line %d) read token: '%s'\n", opbinput->linenumber, opbinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(opbinput != NULL);
   assert(opbinput->npushedtokens < OPB_MAX_PUSHEDTOKENS);

   swapPointers(&opbinput->pushedtokens[opbinput->npushedtokens], &opbinput->token);
   opbinput->npushedtokens++;
}

/** puts the buffered token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushBufferToken(
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(opbinput != NULL);
   assert(opbinput->npushedtokens < OPB_MAX_PUSHEDTOKENS);

   swapPointers(&opbinput->pushedtokens[opbinput->npushedtokens], &opbinput->tokenbuf);
   opbinput->npushedtokens++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(opbinput != NULL);

   swapPointers(&opbinput->token, &opbinput->tokenbuf);
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isEndLine(
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(opbinput != NULL);

   if( *(opbinput->token) ==  ';')
      return TRUE;

   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   OPBINPUT*             opbinput,           /**< OPB reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(opbinput != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( strlen(opbinput->token) == 1 )
   {
      assert(opbinput->token[1] == '\0');

      if( *opbinput->token == '+' )
         return TRUE;
      else if( *opbinput->token == '-' )
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
   OPBINPUT*             opbinput,           /**< OPB reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(opbinput != NULL);
   assert(value != NULL);

   if( strcasecmp(opbinput->token, "INFINITY") == 0 || strcasecmp(opbinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(opbinput->token, &endptr);
      if( endptr != opbinput->token && *endptr == '\0' )
      {
         *value = val;
         if( strlen(opbinput->token) > 18 )
            opbinput->nproblemcoeffs++;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is an equation sense */
static
SCIP_Bool isSense(
   OPBINPUT*             opbinput,           /**< OPB reading data */
   OPBSENSE*             sense               /**< pointer to store the equation sense, or NULL */
   )
{
   assert(opbinput != NULL);

   if( strcmp(opbinput->token, "<") == 0 )
   {
      if( sense != NULL )
         *sense = OPB_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(opbinput->token, ">") == 0 )
   {
      if( sense != NULL )
         *sense = OPB_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(opbinput->token, "=") == 0 )
   {
      if( sense != NULL )
         *sense = OPB_SENSE_EQ;
      return TRUE;
   }

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isStartingSoftConstraintWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(scip != NULL);
   assert(opbinput != NULL);

   if( strcmp(opbinput->token, "[") == 0 )
      return TRUE;

   return FALSE;
}

/** returns whether the current token is a value */
static
SCIP_Bool isEndingSoftConstraintWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput            /**< OPB reading data */
   )
{
   assert(scip != NULL);
   assert(opbinput != NULL);

   if( strcmp(opbinput->token, "]") == 0 )
      return TRUE;

   return FALSE;
}

/** create binary variable with given name */
static
SCIP_RETCODE createVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   char*                 name                /**< name for the variable */
   )
{
   SCIP_VAR* newvar;
   SCIP_Bool dynamiccols;
   SCIP_Bool initial;
   SCIP_Bool removable;

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &dynamiccols) );
   initial = !dynamiccols;
   removable = dynamiccols;

   /* create new variable of the given name */
   SCIPdebugMsg(scip, "creating new variable: <%s>\n", name);

   SCIP_CALL( SCIPcreateVar(scip, &newvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
         initial, removable, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, newvar) );
   *var = newvar;

   /* because the variable was added to the problem, it is captured by SCIP and we
    * can safely release it right now without making the returned *var invalid */
   SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

   return SCIP_OKAY;
}

/** returns the variable with the given name, or creates a new variable if it does not exist */
static
SCIP_RETCODE getVariableOrTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   SCIP_VAR***           vars,               /**< pointer to store the variables */
   int*                  nvars,              /**< pointer to store the number of variables */
   int*                  varssize            /**< pointer to store the varsize, if changed (should already be initialized) */
   )
{
   SCIP_Bool negated;
   char* name;

   assert(scip != NULL);
   assert(opbinput != NULL);
   assert(vars != NULL);
   assert(nvars != NULL);
   assert(varssize != NULL);
   assert(*varssize >= 0);

   *nvars = 0;

   name = opbinput->token;
   assert(name != NULL);

   /* parse AND terms */
   while(!isdigit((unsigned char) *name ) && !isTokenChar(*name) && !opbinput->haserror )
   {
      SCIP_VAR* var;

      negated = FALSE;
      if( *name == '~' )
      {
         negated = TRUE;
         ++name;
      }

      var = SCIPfindVar(scip, name);
      if( var == NULL )
      {
         SCIP_CALL( createVariable(scip, &var, name) );
      }

      if( negated )
      {
         SCIP_VAR* negvar;
         SCIP_CALL( SCIPgetNegatedVar(scip, var, &negvar) );

         var = negvar;
      }

      /* reallocated memory */
      if( *nvars == *varssize )
      {
         *varssize = SCIPcalcMemGrowSize(scip, *varssize + 1);
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, *varssize) );
      }

      (*vars)[*nvars] = var;
      ++(*nvars);

      if( !getNextToken(scip, opbinput) )
         opbinput->haserror = TRUE;

      name = opbinput->token;
   }

   /* check if we found at least on variable */
   if( *nvars == 0 )
      syntaxError(scip, opbinput, "expected a variable name");

   pushToken(opbinput);

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readCoefficients(
   SCIP*const            scip,               /**< SCIP data structure */
   OPBINPUT*const        opbinput,           /**< OPB reading data */
   char*const            name,               /**< pointer to store the name of the line; must be at least of size
                                              *   OPB_MAX_LINELEN */
   SCIP_VAR***           linvars,            /**< pointer to store the array with linear variables (must be freed by caller) */
   SCIP_Real**           lincoefs,           /**< pointer to store the array with linear coefficients (must be freed by caller) */
   int*const             nlincoefs,          /**< pointer to store the number of linear coefficients */
   SCIP_VAR****          terms,              /**< pointer to store the array with nonlinear variables (must be freed by caller) */
   SCIP_Real**           termcoefs,          /**< pointer to store the array with nonlinear coefficients (must be freed by caller) */
   int**                 ntermvars,          /**< pointer to store the number of nonlinear variables in the terms (must be freed by caller) */
   int*const             ntermcoefs,         /**< pointer to store the number of nonlinear coefficients */
   SCIP_Bool*const       newsection,         /**< pointer to store whether a new section was encountered */
   SCIP_Bool*const       isNonlinear,        /**< pointer to store if we have a nonlinear constraint */
   SCIP_Bool*const       issoftcons,         /**< pointer to store whether it is a soft constraint (for wbo files) */
   SCIP_Real*const       weight              /**< pointer to store the weight of the soft constraint */
   )
{
   SCIP_VAR** tmpvars;
   SCIP_Real* tmpcoefs;
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Bool haveweightstart;
   SCIP_Bool haveweightend;
   SCIP_Real coef;
   int coefsign;
   int lincoefssize;
   int termcoefssize;
   int tmpvarssize;
   int ntmpcoefs;
   int ntmpvars;

   assert(opbinput != NULL);
   assert(name != NULL);
   assert(linvars != NULL);
   assert(lincoefs != NULL);
   assert(nlincoefs != NULL);
   assert(terms != NULL);
   assert(termcoefs != NULL);
   assert(ntermvars != NULL);
   assert(ntermcoefs != NULL);
   assert(newsection != NULL);

   *linvars = NULL;
   *lincoefs = NULL;
   *terms = NULL;
   *termcoefs = NULL;
   *ntermvars = NULL;
   *name = '\0';
   *nlincoefs = 0;
   *ntermcoefs = 0;
   *newsection = FALSE;
   *isNonlinear = FALSE;
   *issoftcons = FALSE;

   SCIPdebugMsg(scip, "read coefficients\n");

   /* read the first token, which may be the name of the line */
   if( getNextToken(scip, opbinput) )
   {
      /* remember the token in the token buffer */
      swapTokenBuffer(opbinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(scip, opbinput) )
      {
         if( strcmp(opbinput->token, ":") == 0 )
         {
            /* the second token was a colon ':' the first token is a constraint name */
	    (void)SCIPmemccpy(name, opbinput->tokenbuf, '\0', SCIP_MAXSTRLEN);

            name[SCIP_MAXSTRLEN-1] = '\0';
            SCIPdebugMsg(scip, "(line %d) read constraint name: '%s'\n", opbinput->linenumber, name);

            /* all but the first coefficient need a sign */
            if( strcmp(name, "soft") == 0 && (SCIPgetNVars(scip) > 0 || SCIPgetNConss(scip) > 0) )
            {
               syntaxError(scip, opbinput, "Soft top cost line needs to be the first non-comment line, and without any objective function.\n");
               return SCIP_OKAY;
            }
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse them as coefficients */
            SCIPdebugMsg(scip, "(line %d) constraint has no name\n", opbinput->linenumber);
            pushToken(opbinput);
            pushBufferToken(opbinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it as coefficient */
         pushBufferToken(opbinput);
      }
   }
   else
   {
      assert(SCIPfeof( opbinput->file ) );
      opbinput->eof = TRUE;
      return SCIP_OKAY;
   }

   /* initialize buffers for storing the coefficients */
   lincoefssize = OPB_INIT_COEFSSIZE;
   termcoefssize = OPB_INIT_COEFSSIZE;
   tmpvarssize = OPB_INIT_COEFSSIZE;
   SCIP_CALL( SCIPallocBufferArray(scip, linvars, lincoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, lincoefs, lincoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, terms, termcoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, termcoefs, termcoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, ntermvars, termcoefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, tmpvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcoefs, tmpvarssize) );

   /* read the coefficients */
   coefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   haveweightstart = FALSE;
   haveweightend = FALSE;
   ntmpcoefs = 0;
   ntmpvars = 0;
   while( getNextToken(scip, opbinput) && !hasError(opbinput) )
   {
      if( isEndLine(opbinput) )
      {
         *newsection = TRUE;
         goto TERMINATE;
      }

      /* check if we reached an equation sense */
      if( isSense(opbinput, NULL) )
      {
         /* put the sense back onto the token stack */
         pushToken(opbinput);
         goto TERMINATE;
      }

      /* check if we read a sign */
      if( isSign(opbinput, &coefsign) )
      {
         SCIPdebugMsg(scip, "(line %d) read coefficient sign: %+d\n", opbinput->linenumber, coefsign);
         havesign = TRUE;
         continue;
      }

      /* check if we read a value */
      if( isValue(scip, opbinput, &coef) )
      {
         /* all but the first coefficient need a sign */
         if( (*nlincoefs > 0 || *ntermcoefs > 0 || ntmpcoefs > 0) && !havesign )
         {
            syntaxError(scip, opbinput, "expected sign ('+' or '-') or sense ('<' or '>')");
            goto TERMINATE;
         }

         SCIPdebugMsg(scip, "(line %d) read coefficient value: %g with sign %+d\n", opbinput->linenumber, coef, coefsign);
         if( havevalue )
         {
            syntaxError(scip, opbinput, "two consecutive values");
            goto TERMINATE;
         }
         havevalue = TRUE;

         /* if we read a wbo file, the first line should be something like "soft: <weight>;", where weight is a value or nothing */
         if( strcmp(name, "soft") == 0 )
         {
            assert(ntmpcoefs == 0);

            tmpcoefs[ntmpcoefs] = coefsign * coef;
            ++ntmpcoefs;
         }

         continue;
      }

      /* check if we are reading a soft constraint line, it start with "[<weight>]", where weight is a value */
      if( *nlincoefs == 0 && *ntermcoefs == 0 && ntmpcoefs == 0 && !havesign && !havevalue && strcmp(name, "soft") != 0 && isStartingSoftConstraintWeight(scip, opbinput) )
      {
         if( !opbinput->wbo )
         {
            SCIPwarningMessage(scip, "Found in line %d a soft constraint, without having read a starting top-cost line.\n", opbinput->linenumber);
         }
         haveweightstart = TRUE;

         continue;
      }
      if( *nlincoefs == 0 && *ntermcoefs == 0 && ntmpcoefs == 0 && havevalue && haveweightstart && isEndingSoftConstraintWeight(scip, opbinput) )
      {
         *weight = coefsign * coef;
         SCIPdebugMsg(scip, "(line %d) found soft constraint weight: %g\n", opbinput->linenumber, *weight);

         coefsign = +1;
         havesign = FALSE;
         havevalue = FALSE;
         haveweightend = TRUE;
         *issoftcons = TRUE;

         continue;
      }

      /* if we read a '[' we should already read a ']', which indicates that we read a soft constraint,
       * we have a parsing error */
      if( haveweightstart != haveweightend )
      {
         syntaxError(scip, opbinput, "Wrong soft constraint.");
         goto TERMINATE;
      }

      /* if we read the first non-comment line of a wbo file we should never be here */
      if( strcmp(name, "soft") == 0 )
      {
         syntaxError(scip, opbinput, "Wrong soft top cost line.");
         goto TERMINATE;
      }

      /* the token is a variable name: get the corresponding variables (or create a new ones) */
      SCIP_CALL( getVariableOrTerm(scip, opbinput, &tmpvars, &ntmpvars, &tmpvarssize) );

      if( ntmpvars > 1 )
      {
         /* insert non-linear term */
         *isNonlinear = TRUE;

         SCIPdebugMsg(scip, "(line %d) found linear term: %+g", opbinput->linenumber, coefsign * coef);
#ifndef NDEBUG
         {
            int v;
            for( v = 0; v < ntmpvars; ++v )
            {
               SCIPdebugMsgPrint(scip, " %s * ", SCIPvarGetName(tmpvars[v]));
            }
            SCIPdebugMsgPrint(scip, "\n");
         }
#endif
         if( !SCIPisZero(scip, coef) )
         {
            assert(*ntermcoefs <= termcoefssize);
            /* resize the terms, ntermvars, and termcoefs array if needed */
            if( *ntermcoefs == termcoefssize )
            {
               termcoefssize = SCIPcalcMemGrowSize(scip, termcoefssize + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, terms, termcoefssize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, termcoefs, termcoefssize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, ntermvars, termcoefssize) );
            }
            assert(*ntermcoefs < termcoefssize);

            /* get memory for the last term */
            SCIP_CALL( SCIPallocBufferArray(scip, &((*terms)[*ntermcoefs]), ntmpvars) ); /*lint !e866 */

            /* set the number of variable in this term */
            (*ntermvars)[*ntermcoefs] = ntmpvars;

            /* add all variables */
            for( --ntmpvars; ntmpvars >= 0; --ntmpvars )
            {
               (*terms)[*ntermcoefs][ntmpvars] = tmpvars[ntmpvars];
            }
            /* add coefficient */
            (*termcoefs)[*ntermcoefs] = coefsign * coef;

            /***********************/
            if( !SCIPisIntegral(scip, (*termcoefs)[*ntermcoefs]) )
            {
               SCIPwarningMessage(scip, "coefficient %g in line %d not integral.\n", (*termcoefs)[*ntermcoefs], opbinput->linenumber);
            }

            ++(*ntermcoefs);
         }

         /* reset the flags and coefficient value for the next coefficient */
         coefsign = +1;
         coef = 1.0;
         havesign = FALSE;
         havevalue = FALSE;
         ntmpvars = 0;
      }
      else
      {
         assert(ntmpvars == 1);
         /* insert linear term */
         SCIPdebugMsg(scip, "(line %d) found linear term: %+g<%s>\n", opbinput->linenumber, coefsign * coef, SCIPvarGetName(tmpvars[0]));
         if( !SCIPisZero(scip, coef) )
         {
            assert(*nlincoefs <= lincoefssize);
            /* resize the vars and coefs array if needed */
            if( *nlincoefs >= lincoefssize )
            {
               lincoefssize = SCIPcalcMemGrowSize(scip, lincoefssize + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, linvars, lincoefssize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, lincoefs, lincoefssize) );
            }
            assert(*nlincoefs < lincoefssize);

            /* add coefficient */
            (*linvars)[*nlincoefs] = tmpvars[0];
            (*lincoefs)[*nlincoefs] = coefsign * coef;

            /***********************/
            if( !SCIPisIntegral(scip, (*lincoefs)[*nlincoefs]) )
            {
               SCIPwarningMessage(scip, "coefficient %g in line %d not integral.\n", (*lincoefs)[*nlincoefs], opbinput->linenumber);
            }

            ++(*nlincoefs);
         }

         /* reset the flags and coefficient value for the next coefficient */
         coefsign = +1;
         coef = 1.0;
         havesign = FALSE;
         havevalue = FALSE;
         ntmpvars = 0;
      }
   }

 TERMINATE:
   if( !opbinput->haserror )
   {
      /* all variables should be in the right arrays */
      assert(ntmpvars == 0);
      /* the following is only the case if we read topcost's of a wbo file, we need to move this topcost value to the
       * right array */
      if( ntmpcoefs > 0 )
      {
         /* maximal one topcost value is possible */
         assert(ntmpcoefs == 1);
         /* no other coefficient should be found here */
         assert(*nlincoefs == 0 && *ntermcoefs == 0);

         /* copy value */
         (*lincoefs)[*nlincoefs] = tmpcoefs[0];

         /***********************/
         if( !SCIPisIntegral(scip, (*lincoefs)[*nlincoefs]) )
         {
            SCIPwarningMessage(scip, "topcost not integral.\n");
         }

         *nlincoefs = 1;
      }
   }
   /* clear memory */
   SCIPfreeBufferArray(scip, &tmpcoefs);
   SCIPfreeBufferArray(scip, &tmpvars);

   return SCIP_OKAY;
}

/** set the objective section */
static
SCIP_RETCODE setObjective(
   SCIP*const            scip,               /**< SCIP data structure */
   OPBINPUT*const        opbinput,           /**< OPB reading data */
   const char*           sense,              /**< objective sense */
   SCIP_VAR**const       linvars,            /**< array of linear variables */
   SCIP_Real*const       coefs,              /**< array of objective values for linear variables */
   int const             ncoefs,             /**< number of coefficients for linear part */
   SCIP_VAR***const      terms,              /**< array with nonlinear variables */
   SCIP_Real*const       termcoefs,          /**< array of objective values for nonlinear variables */
   int*const             ntermvars,          /**< number of nonlinear variables in the terms */
   int const             ntermcoefs          /**< number of nonlinear coefficients */
   )
{
   assert(scip != NULL);
   assert(opbinput != NULL);
   assert(isEndLine(opbinput));
   assert(ncoefs == 0 || (linvars != NULL && coefs != NULL));
   assert(ntermcoefs == 0 || (terms != NULL && ntermvars != NULL && termcoefs != NULL));

   if( !hasError(opbinput) )
   {
      SCIP_VAR* var;
      int v;
      char name[SCIP_MAXSTRLEN];

      if( strcmp(sense, "max" ) == 0 )
         opbinput->objsense = SCIP_OBJSENSE_MAXIMIZE;

      /* @todo: what todo with non-linear objectives, maybe create the necessary and-constraints and add the arising linear
       * objective (with and-resultants) or add a integer variable to this constraint and put only this variable in the
       * objective, for this we need to expand the pseudo-boolean constraints to handle integer variables
       *
       * integer variant is not implemented
       */
      if( ntermcoefs > 0 )
      {
#if (LINEAROBJECTIVE == TRUE)
         /* all non-linear parts are created as and-constraints, even if the same non-linear part was already part of the objective function */

         SCIP_VAR** vars;
         int nvars;
         int t;
         SCIP_CONS* andcons;

         for( t = 0; t < ntermcoefs; ++t )
         {
            vars = terms[t];
            nvars = ntermvars[t];
            assert(vars != NULL);
            assert(nvars > 1);

            /* create auxiliary variable */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, ARTIFICIALVARNAMEPREFIX"obj_%d", t);
            SCIP_CALL( SCIPcreateVar(scip, &var, name, 0.0, 1.0, termcoefs[t], SCIP_VARTYPE_BINARY,
                  TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );

            /* @todo: check if it is better to change the branching priority for the artificial variables */
#if 1
            /* change branching priority of artificial variable to -1 */
            SCIP_CALL( SCIPchgVarBranchPriority(scip, var, -1) );
#endif

            /* add auxiliary variable to the problem */
            SCIP_CALL( SCIPaddVar(scip, var) );

#ifdef WITH_DEBUG_SOLUTION
            if( SCIPdebugIsMainscip(scip) )
            {
               SCIP_Real val = 0.0;

               for( v = nvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPdebugGetSolVal(scip, vars[v], &val) );
                  assert(SCIPisFeasZero(scip, val) || SCIPisFeasEQ(scip, val, 1.0));

                  if( val < 0.5 )
                     break;
               }
               SCIP_CALL( SCIPdebugAddSolVal(scip, var, (val < 0.5) ? 0.0 : 1.0) );
            }
#endif

            /* @todo: check whether all constraint creation flags are the best option */
            /* create and-constraint */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "obj_andcons_%d", t);
            SCIP_CALL( SCIPcreateConsAnd(scip, &andcons, name, var, nvars, vars,
                  TRUE, TRUE, TRUE, TRUE, TRUE,
                  FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, andcons) );
            SCIPdebugPrintCons(scip, andcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &andcons) );

            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
#else    /* now the integer variant */
         SCIP_CONS* pseudocons;
         SCIP_Real lb;
         SCIP_Real ub;

         lb = 0.0;
         ub = 0.0;

         /* add all non linear coefficients up */
         for( v = 0; v < ntermcoefs; ++v )
         {
            if( termcoefs[v] < 0 )
               lb += termcoefs[v];
            else
               ub += termcoefs[v];
         }
         /* add all linear coefficients up */
         for( v = 0; v < ncoefs; ++v )
         {
            if( coefs[v] < 0 )
               lb += coefs[v];
            else
               ub += coefs[v];
         }
         assert(lb < ub);

         /* create auxiliary variable */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "artificial_int_obj");
         SCIP_CALL( SCIPcreateVar(scip, &var, name, lb, ub, 1.0, SCIP_VARTYPE_INTEGER,
               TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );

         /* @todo: check if it is better to change the branching priority for the artificial variables */
#if 1
         /* change branching priority of artificial variable to -1 */
         SCIP_CALL( SCIPchgVarBranchPriority(scip, var, -1) );
#endif
         /* add auxiliary variable to the problem */
         SCIP_CALL( SCIPaddVar(scip, var) );

#ifdef WITH_DEBUG_SOLUTION
         if( SCIPdebugIsMainscip(scip) )
         {
            SCIP_Real artval = 0.0;
            SCIP_Real val;

            for( t = 0; t < ntermcoefs; ++t )
            {
               vars = terms[t];
               nvars = ntermvars[t];
               assert(vars != NULL);
               assert(nvars > 1);

               for( v = nvars - 1; v >= 0; --v )
               {
                  SCIP_CALL( SCIPdebugGetSolVal(scip, vars[v], &val) );
                  assert(SCIPisFeasZero(scip, val) || SCIPisFeasEQ(scip, val, 1.0));

                  if( val < 0.5 )
                     break;
               }

               artval += (((val < 0.5) ? 0.0 : 1.0) * termcoefs[t]);
            }
            assert(SCIPisFeasLE(scip, lb, artval) && SCIPisFeasGE(scip, ub, artval));

            SCIP_CALL( SCIPdebugAddSolVal(scip, var, artval) );
         }
#endif


         /* create artificial objection function constraint containing the artificial integer variable */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "artificial_obj_cons");
         SCIP_CALL( SCIPcreateConsPseudoboolean(scip, &pseudocons, name, linvars, ncoefs, coefs, terms, ntermcoefs,
               ntermvars, termcoefs, NULL, 0.0, FALSE, var, 0.0, 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE,
               FALSE, FALSE, FALSE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddCons(scip, pseudocons) );
         SCIPdebugPrintCons(scip, pseudocons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &pseudocons) );

         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         return SCIP_OKAY;
#endif
      }
      /* set the objective values */
      for( v = 0; v < ncoefs; ++v )
      {
	 if( SCIPvarIsNegated(linvars[v]) )
	 {
	    SCIP_VAR* negvar = SCIPvarGetNegationVar(linvars[v]);

	    SCIP_CALL( SCIPaddOrigObjoffset(scip, coefs[v]) );
	    SCIP_CALL( SCIPchgVarObj(scip, negvar, SCIPvarGetObj(negvar) - coefs[v]) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPchgVarObj(scip, linvars[v], SCIPvarGetObj(linvars[v]) + coefs[v]) );
	 }
      }
   }

   return SCIP_OKAY;
}

/** reads the constraints section */
static
SCIP_RETCODE readConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   int*                  nNonlinearConss     /**< pointer to store number of nonlinear constraints */
   )
{
   char name[OPB_MAX_LINELEN];
   SCIP_CONS* cons;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   int nlincoefs;
   SCIP_VAR*** terms;
   SCIP_Real* termcoefs;
   int* ntermvars;
   int ntermcoefs;
   OPBSENSE sense;
   SCIP_RETCODE retcode;
   SCIP_Real sidevalue;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool newsection;
   SCIP_Bool initialconss;
   SCIP_Bool dynamicconss;
   SCIP_Bool dynamicrows;
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool modifiable;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   SCIP_Bool isNonlinear;
   int sidesign;
   SCIP_Bool issoftcons;
   SCIP_Real weight;
   SCIP_VAR* indvar;
   char indname[SCIP_MAXSTRLEN];
   int t;

   assert(scip != NULL);
   assert(opbinput != NULL);
   assert(nNonlinearConss != NULL);

   weight = -SCIPinfinity(scip);
   retcode = SCIP_OKAY;

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, opbinput, name, &linvars, &lincoefs, &nlincoefs, &terms, &termcoefs, &ntermvars, &ntermcoefs, &newsection, &isNonlinear, &issoftcons, &weight) );

   if( hasError(opbinput) || opbinput->eof )
      goto TERMINATE;
   if( newsection )
   {
      if( strcmp(name, "min") == 0 || strcmp(name, "max") == 0 )
      {
         if( opbinput->wbo )
         {
            syntaxError(scip, opbinput, "Cannot have an objective function when having soft constraints.\n");
            goto TERMINATE;
         }

         /* set objective function  */
         SCIP_CALL( setObjective(scip, opbinput, name, linvars, lincoefs, nlincoefs, terms, termcoefs, ntermvars, ntermcoefs) );
      }
      else if( strcmp(name, "soft") == 0 )
      {
         /* we have a "weighted boolean optimization"-file(wbo) */
         opbinput->wbo = TRUE;
         if( nlincoefs == 0 )
            opbinput->topcost = SCIPinfinity(scip);
         else
         {
            assert(nlincoefs == 1);
            opbinput->topcost = lincoefs[0];
         }
         SCIPdebugMsg(scip, "Weighted Boolean Optimization problem has topcost of %g\n", opbinput->topcost);
      }
      else if( nlincoefs > 0 )
         syntaxError(scip, opbinput, "expected constraint sense '=' or '>='");
      goto TERMINATE;
   }

   /* read the constraint sense */
   if( !getNextToken(scip, opbinput) || !isSense(opbinput, &sense) )
   {
      syntaxError(scip, opbinput, "expected constraint sense '=' or '>='");
      goto TERMINATE;
   }

   /* read the right hand side */
   sidesign = +1;
   if( !getNextToken(scip, opbinput) )
   {
      syntaxError(scip, opbinput, "missing right hand side");
      goto TERMINATE;
   }
   if( isSign(opbinput, &sidesign) )
   {
      if( !getNextToken(scip, opbinput) )
      {
         syntaxError(scip, opbinput, "missing value of right hand side");
         goto TERMINATE;
      }
   }
   if( !isValue(scip, opbinput, &sidevalue) )
   {
      syntaxError(scip, opbinput, "expected value as right hand side");
      goto TERMINATE;
   }
   sidevalue *= sidesign;

   /* check if we reached the line end */
   if( !getNextToken(scip, opbinput) || !isEndLine(opbinput) )
   {
      syntaxError(scip, opbinput, "expected endline character ';'");
      goto TERMINATE;
   }

   /* assign the left and right hand side, depending on the constraint sense */
   switch( sense )
   {
   case OPB_SENSE_GE:
      lhs = sidevalue;
      rhs = SCIPinfinity(scip);
      break;
   case OPB_SENSE_LE:
      lhs = -SCIPinfinity(scip);
      rhs = sidevalue;
      break;
   case OPB_SENSE_EQ:
      lhs = sidevalue;
      rhs = sidevalue;
      break;
   case OPB_SENSE_NOTHING:
   default:
      SCIPerrorMessage("invalid constraint sense <%d>\n", sense);
      return SCIP_INVALIDDATA;
   }

   /* create and add the linear constraint */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &initialconss) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &dynamicrows) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/" READER_NAME "/dynamicconss", &dynamicconss) );

   initial = initialconss;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   modifiable = FALSE;
   dynamic = FALSE;/*dynamicconss;*/
   removable = dynamicrows;

   /* create corresponding constraint */
   if( issoftcons )
   {
      (void) SCIPsnprintf(indname, SCIP_MAXSTRLEN, INDICATORVARNAME"%d", opbinput->nindvars);
      ++(opbinput->nindvars);
      SCIP_CALL( createVariable(scip, &indvar, indname) );

      assert(!SCIPisInfinity(scip, -weight));
      SCIP_CALL( SCIPchgVarObj(scip, indvar, weight) );
   }
   else
      indvar = NULL;

   if( ntermcoefs > 0 || issoftcons )
   {
#if GENCONSNAMES == TRUE
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pseudoboolean%d", opbinput->consnumber);
      ++(opbinput->consnumber);
#else
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pseudoboolean");
#endif
      retcode = SCIPcreateConsPseudoboolean(scip, &cons, name, linvars, nlincoefs, lincoefs, terms, ntermcoefs,
            ntermvars, termcoefs, indvar, weight, issoftcons, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE);
      if( retcode != SCIP_OKAY )
         goto TERMINATE;
   }
   else
   {
#if GENCONSNAMES == TRUE
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "linear%d", opbinput->consnumber);
      ++(opbinput->consnumber);
#else
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "linear");
#endif
      retcode = SCIPcreateConsLinear(scip, &cons, name, nlincoefs, linvars, lincoefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE);
      if( retcode != SCIP_OKAY )
         goto TERMINATE;
   }

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIPdebugMsg(scip, "(line %d) created constraint: ", opbinput->linenumber);
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   if( isNonlinear )
      ++(*nNonlinearConss);

 TERMINATE:

   /* free memory */
   for( t = ntermcoefs - 1; t >= 0; --t )
      SCIPfreeBufferArrayNull(scip, &(terms[t]));

   SCIPfreeBufferArrayNull(scip, &ntermvars);
   SCIPfreeBufferArrayNull(scip, &termcoefs);
   SCIPfreeBufferArrayNull(scip, &terms);
   SCIPfreeBufferArrayNull(scip, &lincoefs);
   SCIPfreeBufferArrayNull(scip, &linvars);

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}

/** tries to read the first comment line which usually contains information about the max size of "and" products */
static
SCIP_RETCODE getMaxAndConsDim(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_Bool stop;
   char* commentstart;
   char* nproducts;
   int i;

   assert(scip != NULL);
   assert(opbinput != NULL);

   stop = FALSE;
   commentstart = NULL;
   nproducts = NULL;

   do
   {
      if( SCIPfgets(opbinput->linebuf, (int) sizeof(opbinput->linebuf), opbinput->file) == NULL )
      {
         assert(SCIPfeof( opbinput->file ) );
         break;
      }

      /* read characters after comment symbol */
      for( i = 0; commentchars[i] != '\0'; ++i )
      {
         commentstart = strchr(opbinput->linebuf, commentchars[i]);

         /* found a comment line */
         if( commentstart != NULL )
         {
            /* search for "#product= xyz" in comment line, where xyz represents the number of and constraints */
            nproducts = strstr(opbinput->linebuf, "#product= ");
            if( nproducts != NULL )
            {
	       const char delimchars[] = " \t";
               char* pos;

               nproducts += strlen("#product= ");

               pos = strtok(nproducts, delimchars);

               if( pos != NULL )
               {
                  SCIPdebugMsg(scip, "%d products supposed to be in file.\n", atoi(pos));
               }

               pos = strtok (NULL, delimchars);

               if( pos != NULL && strcmp(pos, "sizeproduct=") == 0 )
               {
                  pos = strtok (NULL, delimchars);
                  if( pos != NULL )
                  {
                     SCIPdebugMsg(scip, "sizeproducts = %d\n", atoi(pos));
                  }
               }

               stop = TRUE;
            }
            break;
         }
      }
   }
   while(commentstart != NULL && !stop);

   opbinput->linebuf[0] = '\0';

#if 0 /* following lines should be correct, but it seems that gzseek does not reset the position if standing at the end of a file */
   /* reset filereader pointer to the beginning */
   (void) SCIPfseek(opbinput->file, 0, SEEK_SET);
#else
   SCIPfclose(opbinput->file);
   opbinput->file = SCIPfopen(filename, "r");
#endif

   return SCIP_OKAY;
}

/** reads an OPB file */
static
SCIP_RETCODE readOPBFile(
   SCIP*                 scip,               /**< SCIP data structure */
   OPBINPUT*             opbinput,           /**< OPB reading data */
   const char*           filename            /**< name of the input file */
   )
{
   int nNonlinearConss;
   int i;

   assert(scip != NULL);
   assert(opbinput != NULL);

   /* open file */
   opbinput->file = SCIPfopen(filename, "r");
   if( opbinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* tries to read the first comment line which usually contains information about the max size of "and" products */
   SCIP_CALL( getMaxAndConsDim(scip, opbinput, filename) );

   /* reading additional information about the number of and constraints in comments to avoid reallocating
    * "opbinput.andconss"
    */
   BMSclearMemoryArray(opbinput->linebuf, OPB_MAX_LINELEN);

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   nNonlinearConss = 0;

   while( !SCIPfeof( opbinput->file ) && !hasError(opbinput) )
   {
      SCIP_CALL( readConstraints(scip, opbinput, &nNonlinearConss) );
   }

   /* if we read a wbo file we need to make sure that the top cost won't be exceeded */
   if( opbinput->wbo )
   {
      SCIP_VAR** topcostvars;
      SCIP_Real* topcosts;
      SCIP_VAR** vars;
      int nvars;
      int ntopcostvars;
      SCIP_Longint topcostrhs;
      SCIP_CONS* topcostcons;

      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
      assert(nvars > 0 || vars != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &topcostvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &topcosts, nvars) );

      ntopcostvars = 0;
      for( i = nvars - 1; i >= 0; --i )
         if( !SCIPisZero(scip, SCIPvarGetObj(vars[i])) )
         {
            topcostvars[ntopcostvars] = vars[i];
            topcosts[ntopcostvars] = SCIPvarGetObj(vars[i]);
            ++ntopcostvars;
         }

      if( SCIPisIntegral(scip, opbinput->topcost) )
         topcostrhs = (SCIP_Longint) SCIPfloor(scip, opbinput->topcost - 1);
      else
         topcostrhs = (SCIP_Longint) SCIPfloor(scip, opbinput->topcost);

      SCIP_CALL( SCIPcreateConsLinear(scip, &topcostcons, TOPCOSTCONSNAME, ntopcostvars, topcostvars, topcosts, -SCIPinfinity(scip),
            (SCIP_Real) topcostrhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, topcostcons) );
      SCIPdebugPrintCons(scip, topcostcons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &topcostcons) );

      SCIPfreeBufferArray(scip, &topcosts);
      SCIPfreeBufferArray(scip, &topcostvars);
   }

   /* close file */
   SCIPfclose(opbinput->file);

   return SCIP_OKAY;
}


/*
 * Local methods (for writing)
 */

/** transforms given and constraint variables to the corresponding active or negated variables */
static
SCIP_RETCODE getBinVarsRepresentatives(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< vars array to get active variables for */
   int const             nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Bool const       transformed         /**< transformed constraint? */
   )
{
   SCIP_Bool negated;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );

   if( transformed )
   {
      for( v = nvars - 1; v >= 0; --v )
      {
         /* gets a binary variable that is equal to the given binary variable, and that is either active, fixed, or
          * multi-aggregated, or the negated variable of an active, fixed, or multi-aggregated variable
          */
         SCIP_CALL( SCIPgetBinvarRepresentative( scip, vars[v], &vars[v], &negated) );
      }
   }
   else
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      for( v = nvars - 1; v >= 0; --v )
      {
         scalar = 1.0;
         constant = 0.0;

         /* retransforms given variable, scalar and constant to the corresponding original variable, scalar and constant,
          * if possible; if the retransformation is impossible, NULL is returned as variable
          */
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalar, &constant) );

         if( vars[v] == NULL )
         {
            SCIPdebugMsg(scip, "A variable couldn't retransformed to an original variable.\n");
            return SCIP_INVALIDDATA;
         }
         if( SCIPisEQ(scip, scalar, -1.0) && SCIPisEQ(scip, constant, 1.0) )
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &vars[v]) );
         }
         else
         {
            if( !SCIPisEQ(scip, scalar, 1.0) || !SCIPisZero(scip, constant) )
            {
               SCIPdebugMsg(scip, "A variable couldn't retransformed to an original variable or a negated variable of an original variable (scalar = %g, constant = %g).\n", scalar, constant);
               return SCIP_INVALIDDATA;
            }
         }
      }
   }

   return SCIP_OKAY;
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
         assert( requiredsize <= *nvars );
      }
   }
   else
      for( v = 0; v < *nvars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );

         if( vars[v] == NULL )
            return SCIP_INVALIDDATA;
      }

   return SCIP_OKAY;
}

/* computes all and-resultants and their corresponding constraint variables */
static
SCIP_RETCODE computeAndConstraintInfos(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_Bool const       transformed,        /**< transformed problem? */
   SCIP_VAR***           resvars,            /**< pointer to store all resultant variables */
   int*                  nresvars,           /**< pointer to store the number of all resultant variables */
   SCIP_VAR****          andvars,            /**< pointer to store to all resultant variables their corresponding active( or negated) and-constraint variables */
   int**                 nandvars,           /**< pointer to store the number of all corresponding and-variables to their corresponding resultant variable */
   SCIP_Bool*const       existandconshdlr,   /**< pointer to store whether the and-constrainthandler exists*/
   SCIP_Bool*const       existands           /**< pointer to store if their exists some and-constraints */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);
   assert(resvars != NULL);
   assert(nresvars != NULL);
   assert(andvars != NULL);
   assert(nandvars != NULL);
   assert(existandconshdlr != NULL);
   assert(existands != NULL);

   *resvars = NULL;
   *nandvars = NULL;
   *andvars = NULL;
   *nresvars = 0;

   /* detect all and-resultants */
   conshdlr = SCIPfindConshdlr(scip, "and");
   if( conshdlr != NULL )
   {
      SCIP_CONS** andconss;
      int nandconss;
      int* shouldnotbeinand;
      int a;
      int c;
      int r;
      int v;
      int pos;
      int ncontainedands;

      andconss = NULL;
      nandconss = 0;
      *existandconshdlr = TRUE;

      /* if we write the original problem we need to get the original and constraints */
      if( !transformed )
      {
         SCIP_CONS** origconss;
         int norigconss;

         origconss = SCIPgetOrigConss(scip);
         norigconss = SCIPgetNOrigConss(scip);

         /* allocate memory for all possible and-constraints */
         SCIP_CALL( SCIPallocBufferArray(scip, &andconss, norigconss) );

         /* collect all original and-constraints */
         for( c = norigconss - 1; c >= 0; --c )
         {
            conshdlr = SCIPconsGetHdlr(origconss[c]);
            assert( conshdlr != NULL );

            if( strcmp(SCIPconshdlrGetName(conshdlr), "and") == 0 )
            {
               andconss[nandconss] = origconss[c];
               ++nandconss;
            }
         }
      }
      else
      {
         nandconss = SCIPconshdlrGetNConss(conshdlr);
         andconss = SCIPconshdlrGetConss(conshdlr);
      }

      assert(andconss != NULL || nandconss == 0);

      *nresvars = nandconss;

      if( nandconss > 0 )
      {
         *existands = TRUE;

         assert(andconss != NULL);

         SCIP_CALL( SCIPallocMemoryArray(scip, resvars, *nresvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, andvars, *nresvars) );
         SCIP_CALL( SCIPallocMemoryArray(scip, nandvars, *nresvars) );

         /* collect all and-constraint variables */
         for( c = nandconss - 1; c >= 0; --c )
         {
            SCIP_VAR** scipandvars;

            assert(andconss[c] != NULL);

            scipandvars = SCIPgetVarsAnd(scip, andconss[c]);
            (*nandvars)[c] = SCIPgetNVarsAnd(scip, andconss[c]);
            SCIP_CALL( SCIPduplicateMemoryArray(scip, &((*andvars)[c]), scipandvars, (*nandvars)[c]) );  /*lint !e866 */
            SCIP_CALL( getBinVarsRepresentatives(scip, (*andvars)[c], (*nandvars)[c], transformed) );

            (*resvars)[c] = SCIPgetResultantAnd(scip, andconss[c]);

            assert((*andvars)[c] != NULL && (*nandvars)[c] > 0);
            assert((*resvars)[c] != NULL);
         }

         /* sorted the array */
         SCIPsortPtrPtrInt((void**)(*resvars), (void**)(*andvars), (*nandvars), SCIPvarComp, (*nresvars));
      }
      else
         *existands = FALSE;

      SCIP_CALL( SCIPallocBufferArray(scip, &shouldnotbeinand, *nresvars) );

      /* check that all and-constraints doesn't contain any and-resultants, if they do try to resolve this */
      /* attention: if resolving leads to x = x*y*... , we can't do anything here ( this only means (... >=x and) y >= x, so normally the and-constraint needs to be
         deleted and the inequality from before needs to be added ) */
      assert(*nandvars != NULL || *nresvars == 0);
      for( r = *nresvars - 1; r >= 0; --r )
      {
         ncontainedands = 0;
         shouldnotbeinand[ncontainedands] = r;
         ++ncontainedands;
         v = 0;

         while( v < (*nandvars)[r] )
         {
            if( SCIPsortedvecFindPtr((void**)(*resvars), SCIPvarComp, (*andvars)[r][v], *nresvars, &pos) )
            {
               /* check if the found position "pos" is equal to an already visited and resultant in this constraint,
                * than here could exist a directed cycle
                */
               /* better use tarjan's algorithm
                *        <http://algowiki.net/wiki/index.php?title=Tarjan%27s_algorithm>,
                *        <http://en.wikipedia.org/wiki/Tarjan%E2%80%99s_strongly_connected_components_algorithm>
                * because it could be that the same resultant is part of this and-constraint and than it would fail
                * without no cycle
                * Note1: tarjans standard algorithm doesn't find cycle from one node to the same;
                * Note2: when tarjan's algorithm find a cycle, it's still possible that this cycle is not "real" e.g.
                *        y = y ~y z (z can also be a product) where y = 0 follows and therefor only "0 = z" is necessary
                */
               for( a = ncontainedands - 1; a >= 0; --a )
                  if( shouldnotbeinand[a] == pos )
                  {
                     SCIPwarningMessage(scip, "This should not happen here. The and-constraint with resultant variable: ");
                     SCIP_CALL( SCIPprintVar(scip, (*resvars)[r], NULL) );
                     SCIPwarningMessage(scip, "possible contains a loop with and-resultant:");
                     SCIP_CALL( SCIPprintVar(scip, (*resvars)[pos], NULL) );

                     /* free memory iff necessary */
                     SCIPfreeBufferArray(scip, &shouldnotbeinand);
                     if( !transformed )
                     {
                        SCIPfreeBufferArray(scip, &andconss);
                     }
                     return SCIP_INVALIDDATA;
                  }
               SCIPdebugMsg(scip, "Another and-constraint contains and-resultant:");
               SCIPdebug( SCIP_CALL( SCIPprintVar(scip, (*resvars)[pos], NULL) ) );
               SCIPdebugMsg(scip, "Trying to resolve.\n");

               shouldnotbeinand[ncontainedands] = pos;
               ++ncontainedands;

               /* try to resolve containing ands */

               /* resize array and number of variables */
               (*nandvars)[r] = (*nandvars)[r] + (*nandvars)[pos] - 1;
               SCIP_CALL( SCIPreallocMemoryArray(scip, &((*andvars)[r]), (*nandvars)[r]) ); /*lint !e866 */

               /* copy all variables */
               for( a = (*nandvars)[pos] - 1; a >= 0; --a )
                  (*andvars)[r][(*nandvars)[r] - a - 1] = (*andvars)[pos][a];

               /* check same position with new variable, so we do not increase v */
            }
	    else
	       ++v;
         }
      }
      SCIPfreeBufferArray(scip, &shouldnotbeinand);

      /* free memory iff necessary */
      if( !transformed )
      {
         SCIPfreeBufferArray(scip, &andconss);
      }
   }
   else
   {
      SCIPdebugMsg(scip, "found no and-constraint-handler\n");
      *existands = FALSE;
      *existandconshdlr = FALSE;
   }

   return SCIP_OKAY;
}

/** clears the given line buffer */
static
void clearBuffer(
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
void writeBuffer(
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
      SCIPinfoMessage(scip, file, "%s", linebuffer);
      clearBuffer(linebuffer, linecnt);
   }
}


/** appends extension to line and prints it to the give file stream if the line buffer get full */
static
void appendBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 linebuffer,         /**< line buffer */
   int*                  linecnt,            /**< number of characters in line */
   const char*           extension           /**< string to extent the line */
   )
{
   assert(scip != NULL);
   assert(linebuffer != NULL);
   assert(linecnt != NULL);
   assert(extension != NULL);

   if( (*linecnt) + strlen(extension) >= OPB_MAX_LINELEN - 1 )
      writeBuffer(scip, file, linebuffer, linecnt);

   /* append extension to linebuffer */
   strncat(linebuffer, extension, OPB_MAX_LINELEN - (unsigned int)(*linecnt) - 1);
   (*linecnt) += (int) strlen(extension);
}

/** write objective function */
static
SCIP_RETCODE writeOpbObjective(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_VAR**const       vars,               /**< array with active (binary) variables */
   int const             nvars,              /**< number of active variables in the problem */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_OBJSENSE const   objsense,           /**< objective sense */
   SCIP_Real const       objscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real const       objoffset,          /**< objective offset from bound shifting and fixing */
   char const*const      multisymbol,        /**< the multiplication symbol to use between coefficient and variable */
   SCIP_Bool const       existands,          /**< does some and-constraints exist? */
   SCIP_Bool const       transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   SCIP_VAR* var;
   char linebuffer[OPB_MAX_LINELEN];
   char buffer[OPB_MAX_LINELEN];
   SCIP_Longint mult;
   SCIP_Bool objective;
   int v;
   int linecnt;
   int pos;

   assert(scip != NULL);
   assert(file != NULL);
   assert(vars != NULL || nvars == 0);
   assert(resvars != NULL || nresvars == 0);
   assert(andvars != NULL || nandvars == NULL);
   assert(multisymbol != NULL);

   mult = 1;
   objective = FALSE;

   clearBuffer(linebuffer, &linecnt);

   /* check if a objective function exits and compute the multiplier to
    * shift the coefficients to integers */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v]; /*lint !e613 */

#ifndef NDEBUG
      {
         /* in case the original problem has to be posted the variables have to be either "original" or "negated" */
         if( !transformed )
            assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL ||
               SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED );
      }
#endif

      /* we found a indicator variable so we assume this is a wbo file */
      if( strstr(SCIPvarGetName(var), INDICATORVARNAME) != NULL )
      {
         /* find the topcost linear inequality which gives us the maximal cost which could be violated by our
          * solution, which is an artificial constraint and print this at first
          *
          * @note: only linear constraint handler is enough in problem stage, otherwise it could be any upgraded linear
          *        constraint which handles pure binary variables
          */
         SCIP_CONSHDLR* conshdlr;
         SCIP_CONS* topcostcons;
         SCIP_Bool printed;

         printed = FALSE;
         topcostcons = SCIPfindCons(scip, TOPCOSTCONSNAME);

         if( topcostcons != NULL )
         {
            conshdlr = SCIPconsGetHdlr(topcostcons);
            assert(conshdlr != NULL);

            if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") == 0 )
               (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "soft: %g;\n", SCIPgetRhsLinear(scip, topcostcons));
            else if( strcmp(SCIPconshdlrGetName(conshdlr), "knapsack") == 0 )
               (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "soft: %g;\n", SCIPgetCapacityKnapsack(scip, topcostcons));
            else if( strcmp(SCIPconshdlrGetName(conshdlr), "setppc") == 0 )
               (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "soft: 1;\n");
            else
            {
               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527 */
            }
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);
            writeBuffer(scip, file, linebuffer, &linecnt);
            printed = TRUE;
         }
         /* following works only in transformed stage */
         else
         {
            /* first try linear constraints */
            conshdlr = SCIPfindConshdlr(scip, "linear");

            if( conshdlr != NULL )
            {
               SCIP_CONS** conss;
               int nconss;
               int c;

               conss = SCIPconshdlrGetConss(conshdlr);
               nconss = SCIPconshdlrGetNConss(conshdlr);

               assert(conss != NULL || nconss == 0);

               for( c = 0; c < nconss; ++c )
               {
                  SCIP_VAR** linvars;
                  int nlinvars;
                  int w;
                  SCIP_Bool topcostfound;
                  SCIP_CONS* cons;

                  cons = conss[c]; /*lint !e613 */
                  assert(cons != NULL);

                  linvars = SCIPgetVarsLinear(scip, cons);
                  nlinvars = SCIPgetNVarsLinear(scip, cons);

                  assert(linvars != NULL || nlinvars == 0);
                  topcostfound = FALSE;

                  for( w = 0; w < nlinvars; ++w )
                  {
                     if( strstr(SCIPvarGetName(linvars[w]), INDICATORVARNAME) != NULL ) /*lint !e613 */
                        topcostfound = TRUE;
                     else
                     {
                        assert(!topcostfound);
                        topcostfound = FALSE;
                     }
                  }

                  if( topcostfound )
                  {
                     (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "soft: %g;\n", SCIPgetRhsLinear(scip, cons));
                     appendBuffer(scip, file, linebuffer, &linecnt, buffer);
                     writeBuffer(scip, file, linebuffer, &linecnt);
                     printed = TRUE;
                     break;
                  }
               }
            }

            if( !printed )
            {
               /* second try knapsack constraints */
               conshdlr = SCIPfindConshdlr(scip, "knapsack");

               if( conshdlr != NULL )
               {
                  SCIP_CONS** conss;
                  int nconss;
                  int c;

                  conss = SCIPconshdlrGetConss(conshdlr);
                  nconss = SCIPconshdlrGetNConss(conshdlr);

                  assert(conss != NULL || nconss == 0);

                  for( c = 0; c < nconss; ++c )
                  {
                     SCIP_VAR** topvars;
                     int ntopvars;
                     int w;
                     SCIP_Bool topcostfound;
                     SCIP_CONS* cons;

                     cons = conss[c]; /*lint !e613 */
                     assert(cons != NULL);

                     topvars = SCIPgetVarsKnapsack(scip, cons);
                     ntopvars = SCIPgetNVarsKnapsack(scip, cons);

                     assert(topvars != NULL || ntopvars == 0);
                     topcostfound = FALSE;

                     for( w = 0; w < ntopvars; ++w )
                     {
                        if( strstr(SCIPvarGetName(topvars[w]), INDICATORVARNAME) != NULL ) /*lint !e613 */
                           topcostfound = TRUE;
                        else
                        {
                           assert(!topcostfound);
                           topcostfound = FALSE;
                        }
                     }

                     if( topcostfound )
                     {
                        (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "soft: %g;\n", SCIPgetCapacityKnapsack(scip, cons));
                        appendBuffer(scip, file, linebuffer, &linecnt, buffer);
                        writeBuffer(scip, file, linebuffer, &linecnt);
                        printed = TRUE;
                        break;
                     }
                  }
               }
            }

            if( !printed )
            {
               /* third try setppc constraints */
               conshdlr = SCIPfindConshdlr(scip, "setppc");

               if( conshdlr != NULL )
               {
                  SCIP_CONS** conss;
                  int nconss;
                  int c;

                  conss = SCIPconshdlrGetConss(conshdlr);
                  nconss = SCIPconshdlrGetNConss(conshdlr);

                  assert(conss != NULL || nconss == 0);

                  for( c = 0; c < nconss; ++c )
                  {
                     SCIP_VAR** topvars;
                     int ntopvars;
                     int w;
                     SCIP_Bool topcostfound;
                     SCIP_CONS* cons;

                     cons = conss[c]; /*lint !e613 */
                     assert(cons != NULL);

                     topvars = SCIPgetVarsSetppc(scip, cons);
                     ntopvars = SCIPgetNVarsSetppc(scip, cons);

                     assert(topvars != NULL || ntopvars == 0);
                     topcostfound = FALSE;

                     for( w = 0; w < ntopvars; ++w )
                     {
                        if( strstr(SCIPvarGetName(topvars[w]), INDICATORVARNAME) != NULL ) /*lint !e613 */
                           topcostfound = TRUE;
                        else
                        {
                           assert(!topcostfound);
                           topcostfound = FALSE;
                        }
                     }

                     if( topcostfound )
                     {
                        (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "soft: 1;\n");
                        appendBuffer(scip, file, linebuffer, &linecnt, buffer);
                        writeBuffer(scip, file, linebuffer, &linecnt);
                        printed = TRUE;
                        break;
                     }
                  }
               }
            }
         }

         /* no topcost constraint found, so print empty topcost line, which means there is no upper bound on violated soft constraints */
         if( !printed )
         {
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "soft: ;\n");
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);
            writeBuffer(scip, file, linebuffer, &linecnt);
         }

         return SCIP_OKAY;
      }

      if( !SCIPisZero(scip, SCIPvarGetObj(var)) )
      {
         objective = TRUE;
         while( !SCIPisIntegral(scip, SCIPvarGetObj(var) * mult) )
         {
            assert(mult * 10 > mult);
            mult *= 10;
         }
      }
   }

   if( objective )
   {
      /* there exist a objective function*/
      SCIPinfoMessage(scip, file, "*   Obj. scale       : %.15g\n", objscale * mult);
      SCIPinfoMessage(scip, file, "*   Obj. offset      : %.15g\n", objoffset);

      clearBuffer(linebuffer, &linecnt);

      /* opb format supports only minimization; therefore, a maximization problem has to be converted */
      if( objsense == SCIP_OBJSENSE_MAXIMIZE )
         mult *= -1;

      SCIPdebugMsg(scip, "print objective function multiplied with %" SCIP_LONGINT_FORMAT "\n", mult);

      appendBuffer(scip, file, linebuffer, &linecnt, "min:");

#ifndef NDEBUG
      if( existands )
      {
         int c;
         /* check that these variables are sorted */
         for( c = nresvars - 1; c > 0; --c )
            assert(SCIPvarGetIndex(resvars[c]) >= SCIPvarGetIndex(resvars[c - 1])); /*lint !e613 */
      }
#endif

      for( v = nvars - 1; v >= 0; --v )
      {
         SCIP_Bool negated;
         var = vars[v]; /*lint !e613 */

         assert(var != NULL);

         if( SCIPisZero(scip, SCIPvarGetObj(var)) )
            continue;

         negated = SCIPvarIsNegated(var);

         assert( linecnt != 0 );

         if( SCIPvarGetObj(var) * mult > SCIP_LONGINT_MAX )
         {
            SCIPerrorMessage("Integral objective value to big (mult = %" SCIP_LONGINT_FORMAT ", value = %g, mult*value = %g, printingvalue = %" SCIP_LONGINT_FORMAT ")for printing in opb format.\n", mult, SCIPvarGetObj(var), SCIPvarGetObj(var) * mult, (SCIP_Longint) SCIPround(scip, SCIPvarGetObj(var) * mult));
         }

         /* replace and-resultant with corresponding variables */
         if( existands && SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, var, nresvars, &pos) )
         {
            int a;

            assert(andvars != NULL);
            assert(nandvars != NULL);
            assert(pos >= 0 && nandvars[pos] > 0 && andvars[pos] != NULL);
            assert(andvars[pos][nandvars[pos] - 1] != NULL);

            negated = SCIPvarIsNegated(andvars[pos][nandvars[pos] - 1]);

            /* print and-vars */
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " %+" SCIP_LONGINT_FORMAT "%s%s%s",
               (SCIP_Longint) (SCIPvarGetObj(var) * mult), multisymbol, negated ? "~" : "",
               strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][nandvars[pos] - 1]) : andvars[pos][nandvars[pos] - 1]), "x"));
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);

            for(a = nandvars[pos] - 2; a >= 0; --a )
            {
               negated = SCIPvarIsNegated(andvars[pos][a]);

               (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][a]) : andvars[pos][a]), "x"));
               appendBuffer(scip, file, linebuffer, &linecnt, buffer);
            }
         }
         else
         {
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " %+" SCIP_LONGINT_FORMAT "%s%s%s",
               (SCIP_Longint) (SCIPvarGetObj(var) * mult), multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x"));
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);
         }
      }

      /* and objective function line ends with a ';' */
      appendBuffer(scip, file, linebuffer, &linecnt, " ;\n");
      writeBuffer(scip, file, linebuffer, &linecnt);
   }

   return SCIP_OKAY;
}

/* print maybe non linear row in OPB format to file stream */
static
SCIP_RETCODE printNLRow(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file (or NULL for standard output) */
   char const*const      type,               /**< row type ("=" or ">=") */
   SCIP_VAR**const       vars,               /**< array of variables */
   SCIP_Real const*const vals,               /**< array of values */
   int const             nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_Longint          weight,             /**< if we found a soft constraint this is the weight, otherwise 0 */
   SCIP_Longint*const    mult,               /**< multiplier for the coefficients */
   char const*const      multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR* var;
   char buffer[OPB_MAX_LINELEN];
   char linebuffer[OPB_MAX_LINELEN + 1];
   int v;
   int pos;
   int linecnt;

   assert(scip != NULL);
   assert(strcmp(type, "=") == 0 || strcmp(type, ">=") == 0);
   assert(mult != NULL);
   assert(resvars != NULL);
   assert(nresvars > 0);
   assert(andvars != NULL && nandvars != NULL);

   clearBuffer(linebuffer, &linecnt);

   /* check if all coefficients are internal; if not commentstart multiplier */
   for( v = 0; v < nvars; ++v )
   {
      while( !SCIPisIntegral(scip, vals[v] * (*mult)) )
      {
         if( ABS(*mult) > ABS(*mult * 10) )
            return SCIP_INVALIDDATA;
         (*mult) *= 10;
      }
   }

   while( !SCIPisIntegral(scip, lhs * (*mult)) )
   {
      if( ABS(*mult) > ABS(*mult * 10) )
         return SCIP_INVALIDDATA;
      (*mult) *= 10;
   }

   /* print comment line if we have to multiply the coefficients to get integrals */
   if( ABS(*mult) != 1 )
      SCIPinfoMessage(scip, file, "* the following constraint is multiplied by %" SCIP_LONGINT_FORMAT " to get integral coefficients\n", ABS(*mult) );

#ifndef NDEBUG
   /* check that these variables are sorted */
   for( v = nresvars - 1; v > 0; --v )
      assert(SCIPvarGetIndex(resvars[v]) >= SCIPvarGetIndex(resvars[v - 1]));
#endif

   /* if we have a soft constraint print the weight*/
   if( weight != 0 )
   {
      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "[%+" SCIP_LONGINT_FORMAT "] ", weight);
      appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print coefficients */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Bool negated;

      var = vars[v];
      assert( var != NULL );

      negated = SCIPvarIsNegated(var);

      /* replace and-resultant with corresponding variables */
      if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, var, nresvars, &pos) )
      {
	 int a;

         assert(andvars != NULL);
         assert(nandvars != NULL);
         assert(pos >= 0 && nandvars[pos] > 0 && andvars[pos] != NULL);
         assert(andvars[pos][nandvars[pos] - 1] != NULL);

	 negated = SCIPvarIsNegated(andvars[pos][nandvars[pos] - 1]);

         if( vals[v] * (*mult) > SCIP_LONGINT_MAX )
         {
            SCIPerrorMessage("Integral coefficient to big (mult = %" SCIP_LONGINT_FORMAT ", value = %g, mult*value = %g, printingvalue = %" SCIP_LONGINT_FORMAT ")for printing in opb format.\n", *mult, vals[v], vals[v] * (*mult), (SCIP_Longint) SCIPround(scip, vals[v] * (*mult)));
         }

	 /* print and-vars */
	 (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+" SCIP_LONGINT_FORMAT "%s%s%s",
	    (SCIP_Longint) SCIPround(scip, vals[v] * (*mult)), multisymbol, negated ? "~" : "",
	    strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][nandvars[pos] - 1]) : andvars[pos][nandvars[pos] - 1]), "x") );
	 appendBuffer(scip, file, linebuffer, &linecnt, buffer);

         for(a = nandvars[pos] - 2; a >= 0; --a )
         {
            negated = SCIPvarIsNegated(andvars[pos][a]);

            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][a]) : andvars[pos][a]), "x"));
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);
         }

         appendBuffer(scip, file, linebuffer, &linecnt, " ");
      }
      else
      {
         (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+" SCIP_LONGINT_FORMAT "%s%s%s ",
            (SCIP_Longint) SCIPround(scip, vals[v] * (*mult)), multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x"));
         appendBuffer(scip, file, linebuffer, &linecnt, buffer);
      }
   }

   /* print left hand side */
   if( SCIPisZero(scip, lhs) )
      lhs = 0.0;

   (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s %" SCIP_LONGINT_FORMAT " ;\n", type, (SCIP_Longint) (lhs * (*mult)) );
   appendBuffer(scip, file, linebuffer, &linecnt, buffer);

   writeBuffer(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}


/** prints given maybe non-linear constraint information in OPB format to file stream */
static
SCIP_RETCODE printNonLinearCons(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**const       vars,               /**< array of variables */
   SCIP_Real*const       vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int const             nvars,              /**< number of variables */
   SCIP_Real const       lhs,                /**< left hand side */
   SCIP_Real const       rhs,                /**< right hand side */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int  const*const      nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_Longint          weight,             /**< if we found a soft constraint this is the weight, otherwise 0 */
   SCIP_Bool const       transformed,        /**< transformed constraint? */
   char const*const      multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   SCIP_Real activeconstant;
   SCIP_Longint mult;
   SCIP_RETCODE retcode;
   int v;
   int nactivevars;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nvars > 0);
   assert(lhs <= rhs);
   assert(resvars != NULL);
   assert(nresvars > 0);
   assert(andvars != NULL && nandvars != NULL);

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   activeconstant = 0.0;
   nactivevars = nvars;

   /* duplicate variable and value array */
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

   mult = 1;
   retcode = SCIP_OKAY;

   /* print row(s) in OPB format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equality constraint */
      retcode = printNLRow(scip, file, "=", activevars, activevals, nactivevars, rhs - activeconstant, resvars,
         nresvars, andvars, nandvars, weight, &mult, multisymbol);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         retcode = printNLRow(scip, file, ">=", activevars, activevals, nactivevars, lhs - activeconstant, resvars,
            nresvars, andvars, nandvars, weight, &mult, multisymbol);
      }


      if( !SCIPisInfinity(scip, rhs) )
      {
         mult *= -1;

         /* print inequality ">=" and multiplying all coefficients by -1 */
         retcode = printNLRow(scip, file, ">=", activevars, activevals, nactivevars, rhs - activeconstant, resvars,
            nresvars, andvars, nandvars, weight, &mult, multisymbol);
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return retcode;
}


/* print row in OPB format to file stream */
static
SCIP_RETCODE printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           type,               /**< row type ("=" or ">=") */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Longint          weight,             /**< if we found a soft constraint this is the weight, otherwise 0 */
   SCIP_Longint*         mult,               /**< multiplier for the coefficients */
   const char*           multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR* var;
   char buffer[OPB_MAX_LINELEN];
   char linebuffer[OPB_MAX_LINELEN + 1];
   int v;
   int linecnt;

   assert(scip != NULL);
   assert(strcmp(type, "=") == 0 || strcmp(type, ">=") == 0);
   assert(mult != NULL);

   clearBuffer(linebuffer, &linecnt);

   /* if we found the topcost linear inequality which gives us the maximal cost which could be violated by our solution,
    * we can stop printing because it is an artificial constraint
    */
   if( nvars > 0 && strstr(SCIPvarGetName(vars[0]), INDICATORVARNAME) != NULL )
      return SCIP_OKAY;

   /* check if all coefficients are integral; if not commentstart multiplier */
   for( v = 0; v < nvars; ++v )
   {
      while( !SCIPisIntegral(scip, vals[v] * (*mult)) )
      {
         if( ABS(*mult) > ABS(*mult * 10) )
            return SCIP_INVALIDDATA;
         (*mult) *= 10;
      }
   }

   while( !SCIPisIntegral(scip, lhs * (*mult)) )
   {
      if( ABS(*mult) > ABS(*mult * 10) )
         return SCIP_INVALIDDATA;
      (*mult) *= 10;
   }

   /* print comment line if we have to multiply the coefficients to get integrals */
   if( ABS(*mult) != 1 )
      SCIPinfoMessage(scip, file, "* the following constraint is multiplied by %" SCIP_LONGINT_FORMAT " to get integral coefficients\n", ABS(*mult) );

   /* if we have a soft constraint print the weight*/
   if( weight != 0 )
   {
      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "[%+" SCIP_LONGINT_FORMAT "] ", weight);
      appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print coefficients */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Bool negated;

      var = vars[v];
      assert( var != NULL );

      negated = SCIPvarIsNegated(var);

      if( vals[v] * (*mult) > SCIP_LONGINT_MAX )
      {
         SCIPerrorMessage("Integral coefficient to big (mult = %" SCIP_LONGINT_FORMAT ", value = %g, mult*value = %g, printingvalue = %" SCIP_LONGINT_FORMAT ")for printing in opb format.\n", *mult, vals[v], vals[v] * (*mult), (SCIP_Longint) SCIPround(scip, vals[v] * (*mult)));
      }

      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+" SCIP_LONGINT_FORMAT "%s%s%s ",
         (SCIP_Longint) SCIPround(scip, vals[v] * (*mult)), multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x"));
      appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print left hand side */
   if( SCIPisZero(scip, lhs) )
      lhs = 0.0;

   (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s %" SCIP_LONGINT_FORMAT " ;\n", type, (SCIP_Longint) (lhs * (*mult)) );
   appendBuffer(scip, file, linebuffer, &linecnt, buffer);

   writeBuffer(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}


/** prints given linear constraint information in OPB format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Longint          weight,             /**< if we found a soft constraint this is the weight, otherwise 0 */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   const char*           multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant;
   SCIP_Longint mult;
   SCIP_RETCODE retcode;
   int v;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   activeconstant = 0.0;

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

   mult = 1;
   retcode = SCIP_OKAY;

   /* print row(s) in OPB format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equality constraint */
      retcode = printRow(scip, file, "=", activevars, activevals, nactivevars, rhs - activeconstant, weight, &mult,
         multisymbol);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         retcode = printRow(scip, file, ">=", activevars, activevals, nactivevars, lhs - activeconstant, weight, &mult,
            multisymbol);
      }

      if( !SCIPisInfinity(scip, rhs) )
      {
         mult *= -1;

         /* print inequality ">=" and multiplying all coefficients by -1 */
         retcode = printRow(scip, file, ">=", activevars, activevals, nactivevars, rhs - activeconstant, weight, &mult,
            multisymbol);
      }
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return retcode;
}

/* print row in OPB format to file stream */
static
SCIP_RETCODE printPBRow(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file (or NULL for standard output) */
   const char*           type,               /**< row type ("=" or ">=") */
   SCIP_VAR**const       linvars,            /**< array of variables */
   SCIP_Real*const       linvals,            /**< array of values */
   int const             nlinvars,           /**< number of variables */
   SCIP_VAR***const      termvars,           /**< term array with array of variables to print */
   int*const             ntermvars,          /**< array with number of variables in each term */
   SCIP_Real*const       termvals,           /**< array of coefficient values for non-linear variables */
   int const             ntermvals,          /**< number non-linear variables in the problem */
   SCIP_Bool**const      negatedarrays,      /**< array of arrays to know which variable in a non-linear part is negated */
   SCIP_VAR*const        indvar,             /**< indicator variable, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Longint*         mult,               /**< multiplier for the coefficients */
   const char*           multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR* var;
   char buffer[OPB_MAX_LINELEN];
   char linebuffer[OPB_MAX_LINELEN + 1];
   int v;
   int t;
   int linecnt;

   assert(scip != NULL);
   assert(strcmp(type, "=") == 0 || strcmp(type, ">=") == 0);
   assert(linvars != NULL || nlinvars == 0);
   assert(linvals != NULL || nlinvars == 0);
   assert(termvars != NULL || ntermvals == 0);
   assert(ntermvars != NULL || ntermvals == 0);
   assert(termvals != NULL || ntermvals == 0);
   assert(negatedarrays != NULL || ntermvals == 0);
   assert(mult != NULL);

   clearBuffer(linebuffer, &linecnt);

   /* if we found the topcost linear inequality which gives us the maximal cost which could be violated by our solution,
    * we can stop printing because it is an artificial constraint
    */
   if( ntermvals == 0 && nlinvars > 0 && strstr(SCIPvarGetName(linvars[0]), INDICATORVARNAME) != NULL ) /*lint !e613 */
      return SCIP_OKAY;

   /* check if all linear coefficients are internal; if not commentstart multiplier */
   for( v = 0; v < nlinvars; ++v )
   {
      while( !SCIPisIntegral(scip, linvals[v] * (*mult)) ) /*lint !e613 */
      {
         if( ABS(*mult) > ABS(*mult * 10) )
            return SCIP_INVALIDDATA;
         (*mult) *= 10;
      }
   }

   /* check if all non-linear coefficients are internal; if not commentstart multiplier */
   for( v = 0; v < ntermvals; ++v )
   {
      while( !SCIPisIntegral(scip, termvals[v] * (*mult)) ) /*lint !e613 */
      {
         if( ABS(*mult) > ABS(*mult * 10) )
            return SCIP_INVALIDDATA;
         (*mult) *= 10;
      }
   }

   while( !SCIPisIntegral(scip, lhs * (*mult)) )
   {
      if( ABS(*mult) > ABS(*mult * 10) )
         return SCIP_INVALIDDATA;
      (*mult) *= 10;
   }

   /* print comment line if we have to multiply the coefficients to get integrals */
   if( ABS(*mult) != 1 )
      SCIPinfoMessage(scip, file, "* the following constraint is multiplied by %" SCIP_LONGINT_FORMAT " to get integral coefficients\n", ABS(*mult) );

   /* if indicator variable exist we have a soft constraint */
   if( indvar != NULL )
   {
      SCIP_Real weight;

      weight = SCIPvarGetObj(indvar);
      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "[%+g] ", weight);
      appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print linear part */
   for( v = 0; v < nlinvars; ++v )
   {
      SCIP_Bool negated;

      var = linvars[v]; /*lint !e613 */
      assert(var != NULL);

      negated = SCIPvarIsNegated(var);

      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+" SCIP_LONGINT_FORMAT "%s%s%s ",
         (SCIP_Longint) SCIPround(scip, linvals[v] * (*mult)), multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x")); /*lint !e613 */
      appendBuffer(scip, file, linebuffer, &linecnt, buffer);
   }

   /* print non-linear part */
   for( t = 0; t < ntermvals; ++t )
   {
      (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%+" SCIP_LONGINT_FORMAT, (SCIP_Longint) SCIPround(scip, termvals[t] * (*mult))); /*lint !e613 */
      appendBuffer(scip, file, linebuffer, &linecnt, buffer);

      for( v = 0; v < ntermvars[t]; ++v ) /*lint !e613 */
      {
         SCIP_Bool negated;

         var = termvars[t][v]; /*lint !e613 */
         assert(var != NULL);

         negated = negatedarrays[t][v]; /*lint !e613 */

         (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(var) : var), "x"));
         appendBuffer(scip, file, linebuffer, &linecnt, buffer);
      }
      appendBuffer(scip, file, linebuffer, &linecnt, " ");
   }

   /* print left hand side */
   if( SCIPisZero(scip, lhs) )
      lhs = 0.0;

   (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s %" SCIP_LONGINT_FORMAT " ;\n", type, (SCIP_Longint) (lhs * (*mult)) );
   appendBuffer(scip, file, linebuffer, &linecnt, buffer);

   writeBuffer(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}


/** prints given pseudo boolean constraint information in OPB format to file stream */
static
SCIP_RETCODE printPseudobooleanCons(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_VAR**const       linvars,            /**< array with variables of linear part */
   SCIP_Real*const       linvals,            /**< array of coefficients values of linear part */
   int const             nlinvars,           /**< number variables in linear part of the problem */
   SCIP_VAR***const      termvars,           /**< term array with array of variables to print */
   int*const             ntermvars,          /**< array with number of variables in each term */
   SCIP_Real*const       termvals,           /**< array of coefficient values for non-linear variables */
   int const             ntermvals,          /**< number non-linear variables in the problem */
   SCIP_VAR*const        indvar,             /**< indicator variable, or NULL */
   SCIP_Real const       lhs,                /**< left hand side of constraint */
   SCIP_Real const       rhs,                /**< right hand side of constraint */
   SCIP_Bool             transformed,        /**< should the transformed problem be printed ? */
   const char*           multisymbol         /**< the multiplication symbol to use between coefficient and variable */
   )
{
   SCIP_VAR*** activetermvars;
   SCIP_Bool** negatedarrays;
   SCIP_VAR** activelinvars;
   SCIP_Real* activelinvals;
   int nactivelinvars;
   SCIP_Real activelinconstant;
   SCIP_Longint mult;
   SCIP_RETCODE retcode;
   int v;

   assert(scip != NULL);
   assert(linvars != NULL || nlinvars == 0);
   assert(linvals != NULL || nlinvars == 0);
   assert(termvars != NULL || 0 == ntermvals);
   assert(ntermvars != NULL || 0 == ntermvals);
   assert(termvals != NULL || 0 == ntermvals);
   assert(lhs <= rhs);

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   activelinconstant = 0.0;

   /* duplicate variable and value array for linear part */
   nactivelinvars = nlinvars;
   if( nactivelinvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activelinvars, linvars, nactivelinvars ) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activelinvals, linvals, nactivelinvars ) );

      /* retransform given variables to active variables */
      SCIP_CALL( getActiveVariables(scip, activelinvars, activelinvals, &nactivelinvars, &activelinconstant, transformed) );
   }
   else
   {
      activelinvars = NULL;
      activelinvals = NULL;
   }

   /* create non-linear information for printing */
   if( ntermvals > 0 )
   {
      assert(termvars != NULL);
      assert(ntermvars != NULL);
      assert(termvals != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &activetermvars, ntermvals) );
      SCIP_CALL( SCIPallocBufferArray(scip, &negatedarrays, ntermvals) );
      for( v = ntermvals - 1; v >= 0; --v )
      {
         assert(ntermvars[v] > 0); /*lint !e613 */

         if( transformed )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(activetermvars[v]), ntermvars[v]) ); /*lint !e866 */
            SCIP_CALL( SCIPallocBufferArray(scip, &(negatedarrays[v]), ntermvars[v]) ); /*lint !e866 */

            /* get binary representatives of binary variables in non-linear terms */
            SCIP_CALL( SCIPgetBinvarRepresentatives(scip, ntermvars[v], termvars[v], activetermvars[v], negatedarrays[v]) );
         }
         else
         {
            SCIP_CALL( SCIPduplicateBufferArray(scip, &(activetermvars[v]), termvars[v], ntermvars[v]) ); /*lint !e866 */
            SCIP_CALL( SCIPallocBufferArray(scip, &(negatedarrays[v]), ntermvars[v]) ); /*lint !e866 */
            BMSclearMemoryArray(negatedarrays[v], ntermvars[v]); /*lint !e866 */
         }
      }
   }
   else
   {
      activetermvars = NULL;
      negatedarrays = NULL;
   }

   mult = 1;
   retcode = SCIP_OKAY;

   /* print row(s) in OPB format */
   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equality constraint */
      retcode = printPBRow(scip, file, "=", activelinvars, activelinvals, nactivelinvars, activetermvars,
         ntermvars, termvals, ntermvals, negatedarrays, indvar, rhs - activelinconstant, &mult, multisymbol);
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         retcode = printPBRow(scip, file, ">=", activelinvars, activelinvals, nactivelinvars, activetermvars,
            ntermvars, termvals, ntermvals, negatedarrays, indvar, lhs - activelinconstant, &mult, multisymbol);
      }

      if( !SCIPisInfinity(scip, rhs) )
      {
         mult *= -1;

         /* print inequality ">=" and multiplying all coefficients by -1 */
         retcode = printPBRow(scip, file, ">=", activelinvars, activelinvals, nactivelinvars, activetermvars,
            ntermvars, termvals, ntermvals, negatedarrays, indvar, rhs - activelinconstant, &mult, multisymbol);
      }
   }


   /* free buffers for non-linear arrays */
   if( ntermvals > 0 )
   {
      assert(negatedarrays != NULL);
      assert(activetermvars != NULL);

      for( v = 0; v < ntermvals; ++v )
      {
	 assert(negatedarrays[v] != NULL);
	 assert(activetermvars[v] != NULL);
         SCIPfreeBufferArray(scip, &(negatedarrays[v]));
         SCIPfreeBufferArray(scip, &(activetermvars[v]));
      }
      SCIPfreeBufferArray(scip, &negatedarrays);
      SCIPfreeBufferArray(scip, &activetermvars);
   }

   /* free buffer for linear arrays */
   if( nactivelinvars > 0 )
   {
      SCIPfreeBufferArray(scip, &activelinvars);
      SCIPfreeBufferArray(scip, &activelinvals);
   }

   return retcode;
}


static
SCIP_RETCODE writeOpbConstraints(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_CONS**const      conss,              /**< array with constraints of the problem */
   int const             nconss,             /**< number of constraints in the problem */
   SCIP_VAR**const       vars,               /**< array with active (binary) variables */
   int const             nvars,              /**< number of active variables in the problem */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   char const*const      multisymbol,        /**< the multiplication symbol to use between coefficient and variable */
   SCIP_Bool const       existandconshdlr,   /**< does and-constrainthandler exist? */
   SCIP_Bool const       existands,          /**< does some and-constraints exist? */
   SCIP_Bool const       transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_RETCODE retcode;
   int nconsvars;
   int v, c;
   SCIP_HASHMAP* linconssofindicatorsmap = NULL;
   SCIP_HASHMAP* linconssofpbsmap = NULL;

   assert(scip != NULL);
   assert(file != NULL);
   assert(conss != NULL || nconss == 0);
   assert(vars != NULL || nvars == 0);
   assert(resvars != NULL || nresvars == 0);
   assert(andvars != NULL || nandvars == 0);
   assert(multisymbol != NULL);

   if( transformed )
   {
      conshdlr = SCIPfindConshdlr(scip, "indicator");

      /* find artificial linear constraints which correspond to indicator constraints to avoid double printing */
      if( conshdlr != NULL )
      {
	 SCIP_CONS** indconss;
	 int nindconss;

	 indconss = SCIPconshdlrGetConss(conshdlr);
	 nindconss = SCIPconshdlrGetNConss(conshdlr);
	 assert(indconss != NULL || nindconss == 0);

	 if( nindconss > 0 )
	 {
	    SCIP_CONS* lincons;

	    /* create the linear constraint of indicator constraints hash map */
	    SCIP_CALL( SCIPhashmapCreate(&linconssofindicatorsmap, SCIPblkmem(scip), nindconss) );
	    assert(indconss != NULL);

	    for( c = 0; c < nindconss; ++c )
	    {
	       assert(indconss[c] != NULL);
	       lincons = SCIPgetLinearConsIndicator(indconss[c]);
	       assert(lincons != NULL);

	       /* insert constraint into mapping between */
	       SCIP_CALL( SCIPhashmapInsert(linconssofindicatorsmap, (void*)lincons, (void*)lincons) );
	    }
	 }
      }

      conshdlr = SCIPfindConshdlr(scip, "pseudoboolean");

      /* find artifical linear constraints which correspond to indicator constraints to avoid double printing */
      if( conshdlr != NULL )
      {
	 SCIP_CONS** pbconss;
	 int npbconss;

	 pbconss = SCIPconshdlrGetConss(conshdlr);
	 npbconss = SCIPconshdlrGetNConss(conshdlr);
	 assert(pbconss != NULL || npbconss == 0);

	 if( npbconss > 0 )
	 {
	    SCIP_CONS* lincons;

	    /* create the linear constraint of indicator constraints hash map */
	    SCIP_CALL( SCIPhashmapCreate(&linconssofpbsmap, SCIPblkmem(scip), npbconss) );

	    for( c = 0; c < npbconss; ++c )
	    {
	       assert(pbconss[c] != NULL); /*lint !e613*/
	       lincons = SCIPgetLinearConsPseudoboolean(scip, pbconss[c]); /*lint !e613*/
	       assert(lincons != NULL);

	       /* insert constraint into mapping between */
	       SCIP_CALL( SCIPhashmapInsert(linconssofpbsmap, (void*)lincons, (void*)lincons) );
	    }
	 }
      }
   }
   /* in original space we cannot ask the constraint handler for its constraints, therefore we have to loop over all
    * original to check for artificial linear once
    */
   else
   {
      SCIP_CONS* lincons;
      SCIP_Bool pbhashmapcreated = FALSE;
      SCIP_Bool indhashmapcreated = FALSE;

      /* loop over all constraint for printing */
      for( c = 0; c < nconss; ++c )
      {
	 conshdlr = SCIPconsGetHdlr(conss[c]); /*lint !e613*/
	 assert(conshdlr != NULL);

	 conshdlrname = SCIPconshdlrGetName(conshdlr);

	 if( strcmp(conshdlrname, "pseudoboolean") == 0 )
	 {
	    if( !pbhashmapcreated )
	    {
	       /* create the linear constraint of indicator constraints hash map */
	       SCIP_CALL( SCIPhashmapCreate(&linconssofpbsmap, SCIPblkmem(scip), nconss) );
	       pbhashmapcreated = TRUE;
	    }

	    lincons = SCIPgetLinearConsPseudoboolean(scip, conss[c]); /*lint !e613*/
	    assert(lincons != NULL);

	    /* insert constraint into mapping between */
	    SCIP_CALL( SCIPhashmapInsert(linconssofpbsmap, (void*)lincons, (void*)lincons) );
	 }
	 else if( strcmp(conshdlrname, "indicator") == 0 )
	 {
	    if( !indhashmapcreated )
	    {
	       /* create the linear constraint of indicator constraints hash map */
	       SCIP_CALL( SCIPhashmapCreate(&linconssofindicatorsmap, SCIPblkmem(scip), nconss) );
	       indhashmapcreated = TRUE;
	    }

	    lincons = SCIPgetLinearConsIndicator(conss[c]); /*lint !e613*/
	    assert(lincons != NULL);

	    /* insert constraint into mapping between */
	    SCIP_CALL( SCIPhashmapInsert(linconssofindicatorsmap, (void*)lincons, (void*)lincons) );
	 }
      }
   }

   retcode = SCIP_OKAY;
   cons = NULL;

   /* loop over all constraint for printing */
   for( c = 0; c < nconss && retcode == SCIP_OKAY; ++c )
   {
      SCIP_CONS* artcons;

      artcons = NULL;

      cons = conss[c]; /*lint !e613 */
      assert(cons != NULL);

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert(transformed == SCIPconsIsTransformed(cons));

      /* in case the transformed is written only constraint are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      if( linconssofpbsmap != NULL )
	 artcons = (SCIP_CONS*) SCIPhashmapGetImage(linconssofpbsmap, (void*)cons);
      if( artcons == NULL && linconssofindicatorsmap != NULL )
	 artcons = (SCIP_CONS*) SCIPhashmapGetImage(linconssofindicatorsmap, (void*)cons);

      if( artcons == NULL )
      {
	 if( strcmp(conshdlrname, "linear") == 0 )
	 {
	    if( SCIPgetNVarsLinear(scip, cons) == 0 )
	    {
	       if( SCIPisGT(scip, SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons)) )
	       {
		  SCIPerrorMessage("Cannot print empty violated constraint %s, %g <= %g is not fulfilled\n",
		     SCIPconsGetName(cons), SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons));
	       }
	       continue;
	    }

	    if( existands )
	    {
	       retcode = printNonLinearCons(scip, file,
                  SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), resvars, nresvars, andvars, nandvars,
                  0LL, transformed, multisymbol);
	    }
	    else
	    {
	       retcode = printLinearCons(scip, file,
                  SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), 0LL, transformed, multisymbol);
	    }
	 }
	 else if( strcmp(conshdlrname, "setppc") == 0 )
	 {
	    consvars = SCIPgetVarsSetppc(scip, cons);
	    nconsvars = SCIPgetNVarsSetppc(scip, cons);

	    if( nconsvars == 0 )
	       continue;

	    switch( SCIPgetTypeSetppc(scip, cons) )
	    {
	    case SCIP_SETPPCTYPE_PARTITIONING :
	       if( existands )
	       {
		  retcode = printNonLinearCons(scip, file, consvars, NULL, nconsvars, 1.0, 1.0, resvars, nresvars,
                     andvars, nandvars, 0LL, transformed, multisymbol);
	       }
	       else
	       {
		  retcode = printLinearCons(scip, file,
                     consvars, NULL, nconsvars, 1.0, 1.0, 0LL, transformed, multisymbol);
	       }
	       break;
	    case SCIP_SETPPCTYPE_PACKING :
	       if( existands )
	       {
		  retcode = printNonLinearCons(scip, file,
                     consvars, NULL, nconsvars, -SCIPinfinity(scip), 1.0, resvars, nresvars, andvars, nandvars,
                     0LL, transformed, multisymbol);
	       }
	       else
	       {
		  retcode = printLinearCons(scip, file,
                     consvars, NULL, nconsvars, -SCIPinfinity(scip), 1.0, 0LL, transformed, multisymbol);
	       }
	       break;
	    case SCIP_SETPPCTYPE_COVERING :
	       if( existands )
	       {
		  retcode = printNonLinearCons(scip, file,
                     consvars, NULL, nconsvars, 1.0, SCIPinfinity(scip), resvars, nresvars, andvars, nandvars,
                     0LL, transformed, multisymbol);
	       }
	       else
	       {
		  retcode = printLinearCons(scip, file,
                     consvars, NULL, nconsvars, 1.0, SCIPinfinity(scip), 0LL, transformed, multisymbol);
	       }
	       break;
	    }
	 }
	 else if( strcmp(conshdlrname, "logicor") == 0 )
	 {
	    if( SCIPgetNVarsLogicor(scip, cons) == 0 )
	       continue;

	    if( existands )
	    {
	       retcode = printNonLinearCons(scip, file,
                  SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons), 1.0, SCIPinfinity(scip),
                  resvars, nresvars, andvars, nandvars, 0LL, transformed, multisymbol);
	    }
	    else
	    {
	       retcode = printLinearCons(scip, file,
                  SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
                  1.0, SCIPinfinity(scip), 0LL, transformed, multisymbol);
	    }
	 }
	 else if( strcmp(conshdlrname, "knapsack") == 0 )
	 {
	    SCIP_Longint* weights;

	    consvars = SCIPgetVarsKnapsack(scip, cons);
	    nconsvars = SCIPgetNVarsKnapsack(scip, cons);

	    if( nconsvars == 0 )
	       continue;

	    /* copy Longint array to SCIP_Real array */
	    weights = SCIPgetWeightsKnapsack(scip, cons);
	    SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
	    for( v = 0; v < nconsvars; ++v )
	       consvals[v] = (SCIP_Real)weights[v];

	    if( existands )
	    {
	       retcode = printNonLinearCons(scip, file, consvars, consvals, nconsvars, -SCIPinfinity(scip),
                  (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), resvars, nresvars, andvars, nandvars,
                  0LL, transformed, multisymbol);
	    }
	    else
	    {
	       retcode = printLinearCons(scip, file, consvars, consvals, nconsvars, -SCIPinfinity(scip),
                  (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), 0LL, transformed, multisymbol);
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

	    if( existands )
	    {
	       retcode = printNonLinearCons(scip, file, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
                  SCIPgetRhsVarbound(scip, cons), resvars, nresvars, andvars, nandvars, 0LL, transformed, multisymbol);
	    }
	    else
	    {
	       retcode = printLinearCons(scip, file, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
		     SCIPgetRhsVarbound(scip, cons), 0LL, transformed, multisymbol);
	    }

	    SCIPfreeBufferArray(scip, &consvars);
	    SCIPfreeBufferArray(scip, &consvals);
	 }
	 else if( strcmp(conshdlrname, "pseudoboolean") == 0 )
	 {
	    SCIP_VAR*** termvars;
	    int* ntermvars;
	    int termvarssize;
	    SCIP_CONS** andconss;
	    SCIP_Real* andcoefs ;
	    SCIP_VAR** linvars;
	    SCIP_Real* lincoefs ;
	    int nlinvars;
	    int t;

	    /* get the required array size for the variables array and for the number of variables in each variable array */
	    termvarssize = SCIPgetNAndsPseudoboolean(scip, cons);
	    assert(termvarssize >= 0);

	    /* allocate temporary memory */
	    SCIP_CALL( SCIPallocBufferArray(scip, &andconss, termvarssize) );
	    SCIP_CALL( SCIPallocBufferArray(scip, &termvars, termvarssize) );
	    SCIP_CALL( SCIPallocBufferArray(scip, &andcoefs, termvarssize) );
	    SCIP_CALL( SCIPallocBufferArray(scip, &ntermvars, termvarssize) );

	    /* get all corresponding and-constraints and therefor all variables */
	    SCIP_CALL( SCIPgetAndDatasPseudoboolean(scip, cons, andconss, andcoefs, &termvarssize) );
	    for( t = termvarssize - 1; t >= 0; --t )
	    {
	       termvars[t] = SCIPgetVarsAnd(scip, andconss[t]);
	       ntermvars[t] = SCIPgetNVarsAnd(scip, andconss[t]);
	    }

	    /* gets number of linear variables without artificial terms variables of pseudoboolean constraint */
	    nlinvars = SCIPgetNLinVarsWithoutAndPseudoboolean(scip, cons);

	    /* allocate temporary memory */
	    SCIP_CALL( SCIPallocBufferArray(scip, &linvars, nlinvars) );
	    SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nlinvars) );

	    /* gets linear constraint of pseudoboolean constraint */
	    SCIP_CALL( SCIPgetLinDatasWithoutAndPseudoboolean(scip, cons, linvars, lincoefs, &nlinvars) );

	    retcode = printPseudobooleanCons(scip, file, linvars, lincoefs, nlinvars,
               termvars, ntermvars, andcoefs, termvarssize, SCIPgetIndVarPseudoboolean(scip, cons),
               SCIPgetLhsPseudoboolean(scip, cons), SCIPgetRhsPseudoboolean(scip, cons), transformed, multisymbol);

	    /* free temporary memory */
	    SCIPfreeBufferArray(scip, &lincoefs);
	    SCIPfreeBufferArray(scip, &linvars);
	    SCIPfreeBufferArray(scip, &ntermvars);
	    SCIPfreeBufferArray(scip, &andcoefs);
	    SCIPfreeBufferArray(scip, &termvars);
	    SCIPfreeBufferArray(scip, &andconss);
	 }
	 else if( strcmp(conshdlrname, "indicator") == 0 )
	 {
	    SCIP_CONS* lincons;
	    SCIP_VAR* indvar;
	    SCIP_VAR* slackvar;
	    SCIP_Longint weight;

	    /* get artificial binary indicator variables */
	    indvar = SCIPgetBinaryVarIndicator(cons);
	    assert(indvar != NULL);

	    if( SCIPvarGetStatus(indvar) == SCIP_VARSTATUS_NEGATED )
	    {
	       indvar = SCIPvarGetNegationVar(indvar);
	       assert(indvar != NULL);
	       assert(SCIPvarGetStatus(indvar) != SCIP_VARSTATUS_AGGREGATED && SCIPvarGetStatus(indvar) != SCIP_VARSTATUS_MULTAGGR);

	       /* get the soft cost of this constraint */
	       weight = (SCIP_Longint) SCIPvarGetObj(indvar);
	    }
	    else
	    {
	       assert(SCIPvarGetStatus(indvar) != SCIP_VARSTATUS_AGGREGATED && SCIPvarGetStatus(indvar) != SCIP_VARSTATUS_MULTAGGR);

	       /* get the soft cost of this constraint */
	       weight = -(SCIP_Longint) SCIPvarGetObj(indvar);
	    }

	    /* get artificial slack variable */
	    slackvar = SCIPgetSlackVarIndicator(cons);
	    assert(slackvar != NULL);

	    /* only need to print indicator constraints with weights on their indicator variable */
	    if( weight != 0 )
	    {
	       SCIP_VAR** scipvarslinear;
	       SCIP_Real* scipvalslinear;
	       SCIP_Bool cont;
	       int nonbinarypos;

	       lincons = SCIPgetLinearConsIndicator(cons);
	       assert(lincons != NULL);

	       nconsvars = SCIPgetNVarsLinear(scip, lincons);
	       scipvarslinear = SCIPgetVarsLinear(scip, lincons);
	       scipvalslinear = SCIPgetValsLinear(scip, lincons);

	       /* allocate temporary memory */
	       SCIP_CALL( SCIPduplicateBufferArray(scip, &consvars, scipvarslinear, nconsvars) );
	       SCIP_CALL( SCIPduplicateBufferArray(scip, &consvals, scipvalslinear, nconsvars) );

	       nonbinarypos = -1;
	       cont = FALSE;

	       /* find non-binary variable */
	       for( v = 0; v < nconsvars; ++v )
	       {
		  if( SCIPvarGetType(consvars[v]) != SCIP_VARTYPE_BINARY )
		  {
		     if( consvars[v] == slackvar )
		     {
			assert(nonbinarypos == -1);
			nonbinarypos = v;
		     }
		     else
		     {
			SCIPwarningMessage(scip, "cannot print linear constraint <%s> of indicator constraint <%s> because it has more than one non-binary variable\n", SCIPconsGetName(lincons), SCIPconsGetName(cons) );
			SCIPinfoMessage(scip, file, "* ");
			SCIP_CALL( SCIPprintCons(scip, cons, file) );
			SCIPinfoMessage(scip, file, ";\n");
			cont = TRUE;
			break;
		     }
		  }
	       }

	       /* if we have not found any non-binary variable we do not print the constraint, maybe we should ??? */
	       if( nonbinarypos == -1 )
	       {
		  SCIPwarningMessage(scip, "cannot print linear constraint <%s> of indicator constraint <%s> because it has no slack variable\n", SCIPconsGetName(lincons), SCIPconsGetName(cons) );
		  SCIPinfoMessage(scip, file, "* ");
		  SCIP_CALL( SCIPprintCons(scip, cons, file) );
		  SCIPinfoMessage(scip, file, ";\n");

		  /* free temporary memory */
		  SCIPfreeBufferArray(scip, &consvals);
		  SCIPfreeBufferArray(scip, &consvars);
		  continue;
	       }

	       /* if the constraint has more than two non-binary variables is not printable and we go to the next */
	       if( cont )
	       {
		  /* free temporary memory */
		  SCIPfreeBufferArray(scip, &consvals);
		  SCIPfreeBufferArray(scip, &consvars);
		  continue;
	       }

	       assert(0 <= nonbinarypos && nonbinarypos < nconsvars);

	       /* remove slackvariable in linear constraint for printing */
	       --nconsvars;
	       consvars[nonbinarypos] = consvars[nconsvars];
	       consvals[nonbinarypos] = consvals[nconsvars];

	       if( existands )
	       {
		  retcode = printNonLinearCons(scip, file,
                     consvars, consvals, nconsvars, SCIPgetLhsLinear(scip, lincons),  SCIPgetRhsLinear(scip, lincons),
                     resvars, nresvars, andvars, nandvars,
                     weight, transformed, multisymbol);
	       }
	       else
	       {
		  retcode = printLinearCons(scip, file,
                     consvars, consvals, nconsvars, SCIPgetLhsLinear(scip, lincons), SCIPgetRhsLinear(scip, lincons),
                     weight, transformed, multisymbol);
	       }

	       /* free temporary memory */
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);
	    }
	    else
	    {
	       SCIPwarningMessage(scip, "indicator constraint <%s> will not be printed because the indicator variable has no objective value(= weight of this soft constraint)\n", SCIPconsGetName(cons) );
	       SCIPinfoMessage(scip, file, "* ");
	       SCIP_CALL( SCIPprintCons(scip, cons, file) );
	       SCIPinfoMessage(scip, file, ";\n");
	    }
	 }
	 else if( strcmp(conshdlrname, "and") == 0 )
	 {
	    /* all resultants of the and constraint will be replaced by all corresponding variables of this constraint,
	     * so no and-constraint will be printed directly */
	    assert(existandconshdlr);
	 }
	 else
	 {
	    SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
	    SCIPinfoMessage(scip, file, "* ");
	    SCIP_CALL( SCIPprintCons(scip, cons, file) );
	    SCIPinfoMessage(scip, file, ";\n");
	 }
      }
   }

   if( retcode == SCIP_INVALIDDATA )
   {
      assert(cons != NULL);

      SCIPerrorMessage("Cannot print constraint %s with non-integral coefficient or sides in opb-format\n",
         SCIPconsGetName(cons));
      SCIP_CALL( SCIPprintCons(scip, cons, stderr) );
      SCIPinfoMessage(scip, file, ";\n");
   }

   if( linconssofpbsmap != NULL )
   {
      /* free hash map */
      SCIPhashmapFree(&linconssofpbsmap);
   }
   if( linconssofindicatorsmap != NULL )
   {
      /* free hash map */
      SCIPhashmapFree(&linconssofindicatorsmap);
   }

   return retcode;
}

/* write fixed variables (unless already done because they are an and resultant or and variable) */
static
SCIP_RETCODE writeOpbFixedVars(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_VAR**            vars,               /**< array with active (binary) variables */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_HASHTABLE*const  printedfixing,      /**< hashmap to store if a fixed variable was already printed */
   char const*const      multisymbol,        /**< the multiplication symbol to use between coefficient and variable */
   SCIP_Bool const       transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   char linebuffer[OPB_MAX_LINELEN];
   char buffer[OPB_MAX_LINELEN];
   int linecnt;
   int v;

   assert(scip != NULL);
   assert(file != NULL);
   assert(vars != NULL || nvars == 0);
   assert(printedfixing != NULL);
   assert(multisymbol != NULL);

   clearBuffer(linebuffer, &linecnt);

   /* print variables which are fixed */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool neg = FALSE;

      assert( vars != NULL );
      var = vars[v];

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
      assert(lb > -0.5 && ub < 1.5);
      assert(SCIPisFeasIntegral(scip, lb));
      assert(SCIPisFeasIntegral(scip, ub));

      /* print fixed and-resultants */
      if( lb > 0.5 || ub < 0.5 )
      {
         if( transformed ) {
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &var, &neg) );
         }

         if( SCIPhashtableExists(printedfixing, (void*)var) )
            continue;

         (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "+1%s%s%s = %g ;\n", multisymbol, neg ? "~" : "", strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(var) : var), "x"), lb);
         appendBuffer(scip, file, linebuffer, &linecnt, buffer);

         /* add variable to the hashmap */
         SCIP_CALL( SCIPhashtableInsert(printedfixing, (void*)var) );
      }
   }

   writeBuffer(scip, file, linebuffer, &linecnt);

   return SCIP_OKAY;
}

/* write and constraints of inactive but relevant and-resultants and and variables which are fixed to one */
static
SCIP_RETCODE writeOpbRelevantAnds(
   SCIP*const            scip,               /**< SCIP data structure */
   FILE*const            file,               /**< output file, or NULL if standard output should be used */
   SCIP_VAR**const       resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_HASHTABLE*const  printedfixing,      /**< hashmap to store if a fixed variable was already printed */
   char const*const      multisymbol,        /**< the multiplication symbol to use between coefficient and variable */
   SCIP_Bool const       transformed         /**< TRUE iff problem is the transformed problem */
   )
{
   SCIP_VAR* resvar;
   SCIP_Longint rhslhs;
   char linebuffer[OPB_MAX_LINELEN];
   char buffer[OPB_MAX_LINELEN];
   int linecnt;
   int r, v;

   assert(scip != NULL);
   assert(file != NULL);
   assert(resvars != NULL || nresvars == 0);
   assert(nandvars != NULL || nresvars == 0);
   assert(andvars != NULL || nandvars == NULL);
   assert(multisymbol != NULL);

   clearBuffer(linebuffer, &linecnt);

   /* print and-variables which are fixed */
   /* @todo remove this block here and the hashtable and let writeOpbFixedVars() do the job? */
   for( r = nresvars - 1; r >= 0; --r )
   {
      SCIP_VAR* var;
      SCIP_Bool neg;
      SCIP_Real lb;
      SCIP_Real ub;

      assert( resvars != NULL );
      resvar = resvars[r];

      if( transformed )
      {
         /* in case the transformed is written only local bounds are posted which are valid in the current node */
         lb = SCIPvarGetLbLocal(resvar);
         ub = SCIPvarGetUbLocal(resvar);
      }
      else
      {
         lb = SCIPvarGetLbOriginal(resvar);
         ub = SCIPvarGetUbOriginal(resvar);
      }

      /* print fixed and-resultants */
      if( lb > 0.5 || ub < 0.5 )
      {
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, resvar, &var, &neg) );

         assert(SCIPisFeasIntegral(scip, lb));
         (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "+1%s%s%s = %g ;\n", multisymbol, neg ? "~" : "", strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(var) : var), "x"), lb);
         appendBuffer(scip, file, linebuffer, &linecnt, buffer);

         /* add variable to the hashmap */
         SCIP_CALL( SCIPhashtableInsert(printedfixing, (void*)var) );
      }

      assert( andvars != NULL && nandvars != NULL );
      assert( andvars[r] != NULL || nandvars[r] == 0 );

      /* print fixed and-variables */
      for( v = nandvars[r] - 1; v >= 0; --v ) /*lint !e613 */
      {
         assert( andvars[r] != NULL );
	 assert( andvars[r][v] != NULL );

         if( transformed )
         {
            /* in case the transformed is written only local bounds are posted which are valid in the current node */
            lb = SCIPvarGetLbLocal(andvars[r][v]);
            ub = SCIPvarGetUbLocal(andvars[r][v]);
         }
         else
         {
            lb = SCIPvarGetLbOriginal(andvars[r][v]);
            ub = SCIPvarGetUbOriginal(andvars[r][v]);
         }

         if( lb > 0.5 || ub < 0.5 )
         {
            SCIP_CALL( SCIPgetBinvarRepresentative(scip, andvars[r][v], &var, &neg) ); /*lint !e613 */

            assert(SCIPisFeasIntegral(scip, lb));
            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "+1%s%s%s = %g ;\n", multisymbol, neg ? "~" : "", strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(var) : var), "x"), lb);
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);

            /* add variable to the hashmap */
            SCIP_CALL( SCIPhashtableInsert(printedfixing, (void*)var) );

         }
      }
   }

   /* print and-constraints with fixed and-resultant to zero and all and-constraints with
    * aggregated resultant, otherwise we would loose this information
    */
   for( r = nresvars - 1; r >= 0; --r )
   {
      assert( resvars != NULL );
      resvar = resvars[r];
      rhslhs = (SCIPvarGetUbLocal(resvar) < 0.5) ? 0 : ((SCIPvarGetLbLocal(resvar) > 0.5) ? 1 : -1);

      /* if and resultant is fixed to 0 and at least one and-variable is fixed to zero, we don't print this redundant constraint */
      if( rhslhs == 0 )
      {
         SCIP_Bool cont;

         cont = FALSE;

         assert( andvars != NULL && nandvars != NULL );
         assert( andvars[r] != NULL || nandvars[r] == 0 );

         /* if resultant variable and one other and variable is already zero, so we did not need to print this and
          * constraint because all other variables are free
          */
         for( v = nandvars[r] - 1; v >= 0; --v ) /*lint !e613 */
	 {
            assert( andvars[r] != NULL );
	    assert( andvars[r][v] != NULL );

            if( SCIPvarGetUbLocal(andvars[r][v]) < 0.5 ) /*lint !e613 */
            {
               cont = TRUE;
               break;
            }
	 }

         if( cont )
            continue;
      }
      /* if and resultant is fixed to 1 and all and-variable are fixed to 1 too, we don't print this redundant constraint */
      else if( rhslhs == 1 )
      {
         SCIP_Bool cont;

         cont = TRUE;

         assert( andvars != NULL && nandvars != NULL );
         assert( andvars[r] != NULL || nandvars[r] == 0 );

         /* if all variables are already fixed to one, we do not need to print this and constraint */
         for( v = nandvars[r] - 1; v >= 0; --v )
	 {
            assert( andvars[r] != NULL );
	    assert( andvars[r][v] != NULL );

            if( SCIPvarGetLbLocal(andvars[r][v]) < 0.5 ) /*lint !e613 */
            {
               cont = FALSE;
               break;
	    }
	 }

         if( cont )
            continue;
      }


      /* print and with fixed or aggregated and-resultant */
      /* rhslhs equals to 0 means the and constraint is relevant due to it's not clear on which values the and variables are
       * rhslhs equals to 1 means the and constraint is irrelevant cause all and variables have to be 1 too
       * rhslhs equals to -1 means the and constraint is relevant cause the variable is only aggregated */
      if( !SCIPvarIsActive(resvar) )
      {
         SCIP_VAR* var;
         SCIP_Bool neg;
         SCIP_Bool firstprinted;

         firstprinted = FALSE;

         assert( andvars != NULL && nandvars != NULL );
         assert( andvars[r] != NULL || nandvars[r] == 0 );

         for( v = nandvars[r] - 1; v >= 0; --v )
         {
            assert( andvars[r] != NULL );
            assert( andvars[r][v] != NULL );

            SCIP_CALL( SCIPgetBinvarRepresentative(scip, andvars[r][v], &var, &neg) ); /*lint !e613 */

            (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", (firstprinted) ? multisymbol : "", neg ? "~" : "", strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(var) : var), "x"));
            appendBuffer(scip, file, linebuffer, &linecnt, buffer);

            firstprinted = TRUE;
         }

         /* if the resultant is aggregated we need to print his binary representation */
         if( rhslhs == -1 )
         {
            int pos;

            assert(transformed);

            SCIP_CALL( SCIPgetBinvarRepresentative(scip, resvar, &resvar, &neg) );

#ifndef NDEBUG
            if( neg )
               assert(SCIPvarIsActive(SCIPvarGetNegationVar(resvar)));
            else
               assert(SCIPvarIsActive(resvar));
#endif

            /* replace and-resultant with corresponding variables */
            if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, neg ? SCIPvarGetNegationVar(resvar) : resvar, nresvars, &pos) )
            {
               SCIP_Bool negated;
               int a;

               assert(andvars != NULL);
               assert(nandvars != NULL);
               assert(pos >= 0 && nandvars[pos] > 0 && andvars[pos] != NULL);
               assert(andvars[pos][nandvars[pos] - 1] != NULL);

               negated = SCIPvarIsNegated(andvars[pos][nandvars[pos] - 1]);

               /* print and-vars */
               (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, neg ? " +1%s%s%s" : " -1%s%s%s", multisymbol, negated ? "~" : "",
                  strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][nandvars[pos] - 1]) : andvars[pos][nandvars[pos] - 1]), "x"));
               appendBuffer(scip, file, linebuffer, &linecnt, buffer);

               for(a = nandvars[pos] - 2; a >= 0; --a )
               {
                  negated = SCIPvarIsNegated(andvars[pos][a]);

                  (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, "%s%s%s", multisymbol, negated ? "~" : "", strstr(SCIPvarGetName(negated ? SCIPvarGetNegationVar(andvars[pos][a]) : andvars[pos][a]), "x"));
                  appendBuffer(scip, file, linebuffer, &linecnt, buffer);
               }

               appendBuffer(scip, file, linebuffer, &linecnt, " ");

               if( neg )
                  rhslhs = 1;
               else
                  rhslhs = 0;
            }
            else
            {
               (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " -1%s%s%s", multisymbol, neg ? "~" : "",
                  strstr(SCIPvarGetName(neg ? SCIPvarGetNegationVar(resvar) : resvar), "x"));
               appendBuffer(scip, file, linebuffer, &linecnt, buffer);

               rhslhs = 0;
            }
         }

         /* print rhslhs */
         (void) SCIPsnprintf(buffer, OPB_MAX_LINELEN, " = %" SCIP_LONGINT_FORMAT " ;\n", rhslhs);
         appendBuffer(scip, file, linebuffer, &linecnt, buffer);

         writeBuffer(scip, file, linebuffer, &linecnt);

      }
   }

   return SCIP_OKAY;
}

/* writes problem to file */
static
SCIP_RETCODE writeOpb(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   const char*           name,               /**< problem name */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   SCIP_OBJSENSE         objsense,           /**< objective sense */
   SCIP_Real             objscale,           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objoffset,          /**< objective offset from bound shifting and fixing */
   SCIP_VAR**            vars,               /**< array with active (binary) variables */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_VAR** const      resvars,            /**< array of resultant variables */
   int const             nresvars,           /**< number of resultant variables */
   SCIP_VAR**const*const andvars,            /**< corresponding array of and-variables */
   int const*const       nandvars,           /**< array of numbers of corresponding and-variables */
   SCIP_Bool const       existandconshdlr,   /**< does and-constrainthandler exist? */
   SCIP_Bool const       existands,          /**< does some and-constraints exist? */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{
   char multisymbol[OPB_MAX_LINELEN];
   SCIP_HASHTABLE* printedfixing;
   SCIP_Bool usesymbol;
   SCIP_RETCODE retcode;

   assert( scip != NULL );
   assert( vars != NULL || nvars == 0 );
   assert( conss != NULL || nconss == 0 );
   assert( result != NULL );

   /* check if should use a multipliers symbol star '*' between coefficients and variables */
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/" READER_NAME "/multisymbol", &usesymbol) );
   (void) SCIPsnprintf(multisymbol, OPB_MAX_LINELEN, "%s", usesymbol ? " * " : " ");

   /* print statistics as comment to file */
   SCIPinfoMessage(scip, file, "* SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "*   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "*   Variables        : %d (all binary)\n", nvars);
   SCIPinfoMessage(scip, file, "*   Constraints      : %d\n", nconss);

   /* create a hash table */
   SCIP_CALL( SCIPhashtableCreate(&printedfixing, SCIPblkmem(scip), nvars,
         SCIPvarGetHashkey, SCIPvarIsHashkeyEq, SCIPvarGetHashkeyVal, NULL) );

   /* write objective function */
   SCIP_CALL( writeOpbObjective(scip, file, vars, nvars, resvars, nresvars, andvars, nandvars,
         objsense, objscale, objoffset, multisymbol, existands, transformed) );

   /* write constraints */
   retcode = writeOpbConstraints(scip, file, conss, nconss, vars, nvars, resvars, nresvars, andvars, nandvars,
      multisymbol, existandconshdlr, existands, transformed);

   if( existands && (retcode == SCIP_OKAY) )
   {
      /* write and constraints of inactive but relevant and-resultants and and-variables which are fixed to one
         with no fixed and resultant */
      SCIP_CALL( writeOpbRelevantAnds(scip, file, resvars, nresvars, andvars, nandvars, printedfixing, multisymbol, transformed) );
   }

   /* write fixed variables */
   SCIP_CALL( writeOpbFixedVars(scip, file, vars, nvars, printedfixing, multisymbol, transformed) );

   SCIPhashtableFree(&printedfixing);

   *result = SCIP_SUCCESS;

   return retcode;
}


/*
 * extern methods
 */

/** reads problem from file */
SCIP_RETCODE SCIPreadOpb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{  /*lint --e{715}*/
   OPBINPUT opbinput;
   SCIP_RETCODE retcode;
   int i;

   /* initialize OPB input data */
   opbinput.file = NULL;
   opbinput.linebuf[0] = '\0';
   SCIP_CALL( SCIPallocBufferArray(scip, &opbinput.token, OPB_MAX_LINELEN) );
   opbinput.token[0] = '\0';
   SCIP_CALL( SCIPallocBufferArray(scip, &opbinput.tokenbuf, OPB_MAX_LINELEN) );
   opbinput.tokenbuf[0] = '\0';
   for( i = 0; i < OPB_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(opbinput.pushedtokens[i]), OPB_MAX_LINELEN) ); /*lint !e866 */
   }

   opbinput.npushedtokens = 0;
   opbinput.linenumber = 1;
   opbinput.bufpos = 0;
   opbinput.linepos = 0;
   opbinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   opbinput.comment = FALSE;
   opbinput.endline = FALSE;
   opbinput.eof = FALSE;
   opbinput.haserror = FALSE;
   opbinput.nproblemcoeffs = 0;
   opbinput.wbo = FALSE;
   opbinput.topcost = -SCIPinfinity(scip);
   opbinput.nindvars = 0;
#if GENCONSNAMES == TRUE
   opbinput.consnumber = 0;
#endif

   /* read the file */
   retcode = readOPBFile(scip, &opbinput, filename);

   /* free dynamically allocated memory */
   for( i = OPB_MAX_PUSHEDTOKENS - 1; i >= 0; --i )
   {
      SCIPfreeBufferArrayNull(scip, &(opbinput.pushedtokens[i]));
   }
   SCIPfreeBufferArrayNull(scip, &opbinput.tokenbuf);
   SCIPfreeBufferArrayNull(scip, &opbinput.token);

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   SCIP_CALL( retcode );

   if( opbinput.nproblemcoeffs > 0 )
   {
      SCIPwarningMessage(scip, "there might be <%d> coefficients or weight out of range!\n", opbinput.nproblemcoeffs);
   }

   /* evaluate the result */
   if( opbinput.haserror )
      return SCIP_READERROR;
   else
   {
      /* set objective sense */
      SCIP_CALL( SCIPsetObjsense(scip, opbinput.objsense) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}

/** writes problem to file */
SCIP_RETCODE SCIPwriteOpb(
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
   SCIP_VAR**            fixedvars,          /**< array with fixed variables */
   int                   nfixedvars,         /**< number of fixed and aggregated variables in the problem */
   SCIP_CONS**           conss,              /**< array with constraints of the problem */
   int                   nconss,             /**< number of constraints in the problem */
   SCIP_Bool             genericnames,       /**< should generic variable and constraint names be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file writing call */
   )
{  /*lint --e{715}*/
   SCIP_RETCODE retcode = SCIP_OKAY;

   if( nvars != nbinvars && (nintvars > 0 || SCIPfindConshdlr(scip, "indicator") != NULL
         || ncontvars + nimplvars != SCIPconshdlrGetNConss(SCIPfindConshdlr(scip, "indicator"))) )
   {
      SCIPwarningMessage(scip, "only binary problems can be written in OPB format.\n");
      *result = SCIP_DIDNOTRUN;
   }
   else
   {
      SCIP_VAR*** andvars;
      SCIP_VAR** resvars;
      int* nandvars;
      SCIP_Bool existands;
      SCIP_Bool existandconshdlr;
      int nresvars;
      int v;

      /* computes all and-resultants and their corresponding constraint variables */
      SCIP_CALL( computeAndConstraintInfos(scip, transformed, &resvars, &nresvars, &andvars, &nandvars, &existandconshdlr, &existands) );

      if( genericnames )
      {
#ifndef NDEBUG
         /* check for correct names for opb-format */
         int idx;
         int pos;

         for( v = nvars - 1; v >= 0; --v )
         {
            if( existands )
            {
               /* and variables are artificial */
               if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
                  continue;
            }

            assert(sscanf(SCIPvarGetName(vars[v]), "x%d", &idx) == 1);
         }
#endif
         retcode = writeOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
               nvars, conss, nconss, resvars, nresvars, andvars, nandvars, existandconshdlr, existands, result);
      }
      else
      {
         SCIP_Bool printed;
         int idx;
         int pos;

         printed = FALSE;

         /* check if there are already generic names for all (not fixed variables)*/
         for( v = nvars - 1; v >= 0; --v )
            if( !existands || !SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
            {
               if( sscanf(SCIPvarGetName(vars[v]), transformed ? "t_x%d" : "x%d", &idx) != 1 && strstr(SCIPvarGetName(vars[v]), INDICATORVARNAME) == NULL && strstr(SCIPvarGetName(vars[v]), INDICATORSLACKVARNAME) == NULL )
               {
                  SCIPwarningMessage(scip, "At least following variable name isn't allowed in opb format.\n");
                  SCIP_CALL( SCIPprintVar(scip, vars[v], NULL) );
                  SCIPwarningMessage(scip, "OPB format needs generic variable names!\n");

                  if( transformed )
                  {
                     SCIPwarningMessage(scip, "write transformed problem with generic variable names.\n");
                     SCIP_CALL( SCIPprintTransProblem(scip, file, "opb", TRUE) );
                  }
                  else
                  {
                     SCIPwarningMessage(scip, "write original problem with generic variable names.\n");
                     SCIP_CALL( SCIPprintOrigProblem(scip, file, "opb", TRUE) );
                  }
                  printed = TRUE;
                  break;
               }
            }

         if( !printed )
         {
            /* check if there are already generic names for all (fixed variables)*/
            for( v = nfixedvars - 1; v >= 0; --v )
               if( !existands || !SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
               {
                  if( sscanf(SCIPvarGetName(fixedvars[v]), transformed ? "t_x%d" : "x%d", &idx) != 1 && strstr(SCIPvarGetName(fixedvars[v]), INDICATORVARNAME) == NULL && strstr(SCIPvarGetName(fixedvars[v]), INDICATORSLACKVARNAME) == NULL )
                  {
                     SCIPwarningMessage(scip, "At least following variable name isn't allowed in opb format.\n");
                     SCIP_CALL( SCIPprintVar(scip, fixedvars[v], NULL) );
                     SCIPwarningMessage(scip, "OPB format needs generic variable names!\n");

                     if( transformed )
                     {
                        SCIPwarningMessage(scip, "write transformed problem with generic variable names.\n");
                        SCIP_CALL( SCIPprintTransProblem(scip, file, "opb", TRUE) );
                     }
                     else
                     {
                        SCIPwarningMessage(scip, "write original problem with generic variable names.\n");
                        SCIP_CALL( SCIPprintOrigProblem(scip, file, "opb", TRUE) );
                     }
                     printed = TRUE;
                     break;
                  }
               }
         }

         if( !printed )
         {
#ifndef NDEBUG
            for( v = nvars - 1; v >= 0; --v )
            {
               if( existands )
               {
                  if( SCIPsortedvecFindPtr((void**)resvars, SCIPvarComp, vars[v], nresvars, &pos) )
                     continue;
               }

               assert(sscanf(SCIPvarGetName(vars[v]), transformed ? "t_x%d" : "x%d", &idx) == 1 || strstr(SCIPvarGetName(vars[v]), INDICATORVARNAME) != NULL || strstr(SCIPvarGetName(vars[v]), INDICATORSLACKVARNAME) != NULL );
            }
#endif
            retcode = writeOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
               nvars, conss, nconss, resvars, nresvars, andvars, nandvars, existandconshdlr, existands, result);
         }
      }

      if( existands )
      {
         /* free temporary buffers */
         assert(resvars != NULL);
         assert(andvars != NULL);
         assert(nandvars != NULL);

         for( v = nresvars - 1; v >= 0; --v )
         {
            assert(andvars[v] != NULL);
            SCIPfreeMemoryArray(scip, &andvars[v]);
         }
         SCIPfreeMemoryArray(scip, &nandvars);
         SCIPfreeMemoryArray(scip, &andvars);
         SCIPfreeMemoryArray(scip, &resvars);
      }

      *result = SCIP_SUCCESS;
   }

   if( retcode == SCIP_INVALIDDATA )
      return SCIP_WRITEERROR;

   return retcode;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyOpb)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderOpb(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadOpb)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadOpb(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteOpb)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPwriteOpb(scip, file, name, transformed, objsense, objscale, objoffset, vars,
         nvars, nbinvars, nintvars, nimplvars, ncontvars, fixedvars, nfixedvars, conss, nconss, genericnames, result) );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the opb file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderOpb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyOpb) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadOpb) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteOpb) );

   /* add opb reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/dynamicconss", "should model constraints be subject to aging?",
         NULL, FALSE, FALSE/*TRUE*/, NULL, NULL) ); /* have to be FALSE, otherwise an error might inccur in restart during branch and bound */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/multisymbol", "use '*' between coefficients and variables by writing to problem?",
         NULL, TRUE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
