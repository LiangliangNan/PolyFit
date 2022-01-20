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

/**@file   reader_diff.c
 * @brief  DIFF file reader
 * @author Jakob Witzig
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

#include "scip/reader_diff.h"

#define READER_NAME             "diffreader"
#define READER_DESC             "file reader for changes in the LP file"
#define READER_EXTENSION        "diff"

/*
 * Data structures
 */
#define LP_MAX_LINELEN       65536
#define LP_MAX_PUSHEDTOKENS  2
#define LP_INIT_COEFSSIZE    8192

/** Section in LP File */
enum LpSection
{
   LP_START, LP_OBJECTIVE, LP_END
};
typedef enum LpSection LPSECTION;

enum LpExpType
{
   LP_EXP_NONE, LP_EXP_UNSIGNED, LP_EXP_SIGNED
};
typedef enum LpExpType LPEXPTYPE;

enum LpSense
{
   LP_SENSE_LE, LP_SENSE_GE, LP_SENSE_EQ
};
typedef enum LpSense LPSENSE;

/** LP reading data */
struct LpInput
{
   SCIP_FILE*            file;
   char                  linebuf[LP_MAX_LINELEN+1];
   char                  probname[LP_MAX_LINELEN];
   char                  objname[LP_MAX_LINELEN];
   char*                 token;
   char*                 tokenbuf;
   char*                 pushedtokens[LP_MAX_PUSHEDTOKENS];
   int                   npushedtokens;
   int                   linenumber;
   int                   linepos;
   LPSECTION             section;
   SCIP_OBJSENSE         objsense;
   SCIP_Bool             haserror;
   SCIP_Bool             comment;
   SCIP_Bool             endline;
};
typedef struct LpInput LPINPUT;

static const char commentchars[] = "\\";


/*
 * Local methods (for reading)
 */

/** issues an error message and marks the LP data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   const char*           msg                 /**< error message */
   )
{
   char formatstr[256];

   assert(lpinput != NULL);

   SCIPerrorMessage("Syntax error in line %d ('%s'): %s \n", lpinput->linenumber, lpinput->token, msg);
   if( lpinput->linebuf[strlen(lpinput->linebuf)-1] == '\n' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s", lpinput->linebuf);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "  input: %s\n", lpinput->linebuf);
   }
   (void) SCIPsnprintf(formatstr, 256, "         %%%ds\n", lpinput->linepos);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, (const char*)formatstr, "^");
   lpinput->section  = LP_END;
   lpinput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   return lpinput->haserror;
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
   case '*':
   case '^':
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
   LPEXPTYPE*            exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit((unsigned char)c) )
      return TRUE;
   else if( (*exptype == LP_EXP_NONE) && !(*hasdot) && (c == '.') && isdigit((unsigned char)nextc) )
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == LP_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = LP_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit((unsigned char)nextc) )
      {
         *exptype = LP_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == LP_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = LP_EXP_UNSIGNED;
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
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   int i;

   assert(lpinput != NULL);

   /* if we previously detected a comment we have to parse the remaining line away if there is something left */
   if( !lpinput->endline && lpinput->comment )
   {
      SCIPdebugMsg(scip, "Throwing rest of comment away.\n");

      do
      {
         lpinput->linebuf[LP_MAX_LINELEN-2] = '\0';
         (void)SCIPfgets(lpinput->linebuf, (int) sizeof(lpinput->linebuf), lpinput->file);
      }
      while( lpinput->linebuf[LP_MAX_LINELEN-2] != '\0' );

      lpinput->comment = FALSE;
      lpinput->endline = TRUE;
   }

   /* read next line */
   lpinput->linepos = 0;
   lpinput->linebuf[LP_MAX_LINELEN-2] = '\0';

   if( SCIPfgets(lpinput->linebuf, (int) sizeof(lpinput->linebuf), lpinput->file) == NULL )
   {
      /* clear the line, this is really necessary here! */
      BMSclearMemoryArray(lpinput->linebuf, LP_MAX_LINELEN);

      return FALSE;
   }

   lpinput->linenumber++;

   /* if line is too long for our buffer correct the buffer and correct position in file */
   if( lpinput->linebuf[LP_MAX_LINELEN-2] != '\0' )
   {
      char* last;

      /* buffer is full; erase last token since it might be incomplete */
      lpinput->endline = FALSE;
      last = strrchr(lpinput->linebuf, ' ');

      if( last == NULL )
      {
         SCIPwarningMessage(scip, "we read %d characters from the file; this might indicate a corrupted input file!",
            LP_MAX_LINELEN - 2);
         lpinput->linebuf[LP_MAX_LINELEN-2] = '\0';
         SCIPdebugMsg(scip, "the buffer might be corrupted\n");
      }
      else
      {
         SCIPfseek(lpinput->file, -(long) strlen(last) - 1, SEEK_CUR);
         SCIPdebugMsg(scip, "correct buffer, reread the last %ld characters\n", (long) strlen(last) + 1);
         *last = '\0';
      }
   }
   else
   {
      /* found end of line */
      lpinput->endline = TRUE;
   }
   lpinput->linebuf[LP_MAX_LINELEN-1] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
   lpinput->comment = FALSE;

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(lpinput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */

         lpinput->comment = TRUE;
         break;
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
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   SCIP_Bool hasdot;
   LPEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(lpinput != NULL);
   assert(lpinput->linepos < LP_MAX_LINELEN);

   /* check the token stack */
   if( lpinput->npushedtokens > 0 )
   {
      swapPointers(&lpinput->token, &lpinput->pushedtokens[lpinput->npushedtokens-1]);
      lpinput->npushedtokens--;

      SCIPdebugMsg(scip, "(line %d) read token again: '%s'\n", lpinput->linenumber, lpinput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = lpinput->linebuf;
   while( isDelimChar(buf[lpinput->linepos]) )
   {
      if( buf[lpinput->linepos] == '\0' )
      {
         if( !getNextLine(scip, lpinput) )
         {
            lpinput->section = LP_END;
            SCIPdebugMsg(scip, "(line %d) end of file\n", lpinput->linenumber);
            return FALSE;
         }
         assert(lpinput->linepos == 0);
      }
      else
         lpinput->linepos++;
   }
   assert(lpinput->linepos < LP_MAX_LINELEN);
   assert(!isDelimChar(buf[lpinput->linepos]));

   /* check if the token is a value */
   hasdot = FALSE;
   exptype = LP_EXP_NONE;
   if( isValueChar(buf[lpinput->linepos], buf[lpinput->linepos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < LP_MAX_LINELEN);
         assert(!isDelimChar(buf[lpinput->linepos]));
         lpinput->token[tokenlen] = buf[lpinput->linepos];
         tokenlen++;
         lpinput->linepos++;
      }
      while( isValueChar(buf[lpinput->linepos], buf[lpinput->linepos+1], FALSE, &hasdot, &exptype) );
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < LP_MAX_LINELEN);
         lpinput->token[tokenlen] = buf[lpinput->linepos];
         tokenlen++;
         lpinput->linepos++;
         if( tokenlen == 1 && isTokenChar(lpinput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[lpinput->linepos]) && !isTokenChar(buf[lpinput->linepos]) );

      /* if the token is a power sign '^', skip a following '2'
       * if the token is an equation sense '<', '>', or '=', skip a following '='
       * if the token is an equality token '=' and the next character is a '<' or '>', replace the token by the inequality sense
       */
      if( tokenlen >= 1 && lpinput->token[tokenlen-1] == '^' && buf[lpinput->linepos] == '2' )
      {
         lpinput->linepos++;
      }
      if( tokenlen >= 1
         && (lpinput->token[tokenlen-1] == '<' || lpinput->token[tokenlen-1] == '>' || lpinput->token[tokenlen-1] == '=')
         && buf[lpinput->linepos] == '=' )
      {
         lpinput->linepos++;
      }
      else if( lpinput->token[tokenlen-1] == '=' && (buf[lpinput->linepos] == '<' || buf[lpinput->linepos] == '>') )
      {
         lpinput->token[tokenlen-1] = buf[lpinput->linepos];
         lpinput->linepos++;
      }
   }
   assert(tokenlen < LP_MAX_LINELEN);
   lpinput->token[tokenlen] = '\0';

   SCIPdebugMsg(scip, "(line %d) read token: '%s'\n", lpinput->linenumber, lpinput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);
   assert(lpinput->npushedtokens < LP_MAX_PUSHEDTOKENS);

   swapPointers(&lpinput->pushedtokens[lpinput->npushedtokens], &lpinput->token);
   lpinput->npushedtokens++;
}

/** puts the buffered token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushBufferToken(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);
   assert(lpinput->npushedtokens < LP_MAX_PUSHEDTOKENS);

   swapPointers(&lpinput->pushedtokens[lpinput->npushedtokens], &lpinput->tokenbuf);
   lpinput->npushedtokens++;
}

/** swaps the current token with the token buffer */
static
void swapTokenBuffer(
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   swapPointers(&lpinput->token, &lpinput->tokenbuf);
}

/** checks whether the current token is a section identifier, and if yes, switches to the corresponding section */
static
SCIP_Bool isNewSection(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   SCIP_Bool iscolon;
   size_t len;

   assert(lpinput != NULL);

   /* remember first token by swapping the token buffer */
   swapTokenBuffer(lpinput);

   /* look at next token: if this is a ':', the first token is a name and no section keyword */
   iscolon = FALSE;
   if( getNextToken(scip, lpinput) )
   {
      iscolon = (*lpinput->token == ':');
      pushToken(lpinput);
   }

   /* reinstall the previous token by swapping back the token buffer */
   swapTokenBuffer(lpinput);

   /* check for ':' */
   if( iscolon )
      return FALSE;

   len = strlen(lpinput->token);
   assert(len < LP_MAX_LINELEN);

   /* the section keywords are at least 2 characters up to 8 or exactly 15 characters long */
   if( len > 1 && (len < 9 || len == 15) )
   {
      char token[16];
      int c = 0;

      while( lpinput->token[c] != '\0' )
      {
         token[c] = toupper(lpinput->token[c]); /*lint !e734*/
         ++c;
         assert(c < 16);
      }
      token[c] = '\0';

      if( (len == 3 && strcmp(token, "MIN") == 0)
         || (len == 7 && strcmp(token, "MINIMUM") == 0)
         || (len == 8 && strcmp(token, "MINIMIZE") == 0) )
      {
         SCIPdebugMsg(scip, "(line %d) new section: OBJECTIVE\n", lpinput->linenumber);
         lpinput->section = LP_OBJECTIVE;
         lpinput->objsense = SCIP_OBJSENSE_MINIMIZE;
         return TRUE;
      }

      if( (len == 3 && strcmp(token, "MAX") == 0)
         || (len == 7 && strcmp(token, "MAXIMUM") == 0)
         || (len == 8 && strcmp(token, "MAXIMIZE") == 0) )
      {
         SCIPdebugMsg(scip, "(line %d) new section: OBJECTIVE\n", lpinput->linenumber);
         lpinput->section = LP_OBJECTIVE;
         lpinput->objsense = SCIP_OBJSENSE_MAXIMIZE;
         return TRUE;
      }

      if( len == 3 && strcmp(token, "END") == 0 )
      {
         SCIPdebugMsg(scip, "(line %d) new section: END\n", lpinput->linenumber);
         lpinput->section = LP_END;
         return TRUE;
      }
   }

   return FALSE;
}

/** returns whether the current token is a sign */
static
SCIP_Bool isSign(
   LPINPUT*              lpinput,            /**< LP reading data */
   int*                  sign                /**< pointer to update the sign */
   )
{
   assert(lpinput != NULL);
   assert(sign != NULL);
   assert(*sign == +1 || *sign == -1);

   if( lpinput->token[1] == '\0' )
   {
      if( *lpinput->token == '+' )
         return TRUE;
      else if( *lpinput->token == '-' )
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
   LPINPUT*              lpinput,            /**< LP reading data */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   assert(lpinput != NULL);
   assert(value != NULL);

   if( strcasecmp(lpinput->token, "INFINITY") == 0 || strcasecmp(lpinput->token, "INF") == 0 )
   {
      *value = SCIPinfinity(scip);
      return TRUE;
   }
   else
   {
      double val;
      char* endptr;

      val = strtod(lpinput->token, &endptr);
      if( endptr != lpinput->token && *endptr == '\0' )
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
   LPINPUT*              lpinput,            /**< LP reading data */
   LPSENSE*              sense               /**< pointer to store the equation sense, or NULL */
   )
{
   assert(lpinput != NULL);

   if( strcmp(lpinput->token, "<") == 0 )
   {
      if( sense != NULL )
         *sense = LP_SENSE_LE;
      return TRUE;
   }
   else if( strcmp(lpinput->token, ">") == 0 )
   {
      if( sense != NULL )
         *sense = LP_SENSE_GE;
      return TRUE;
   }
   else if( strcmp(lpinput->token, "=") == 0 )
   {
      if( sense != NULL )
         *sense = LP_SENSE_EQ;
      return TRUE;
   }

   return FALSE;
}

/** returns the variable with the given name, or creates a new variable if it does not exist */
static
SCIP_RETCODE getVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 name,               /**< name of the variable */
   SCIP_VAR**            var                 /**< pointer to store the variable */
   )
{
   assert(name != NULL);
   assert(var != NULL);

   *var = SCIPfindVar(scip, name);

   if( *var == NULL )
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads the header of the file */
static
SCIP_RETCODE readStart(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   assert(lpinput != NULL);

   /* everything before first section is treated as comment */
   do
   {
      /* get token */
      if( !getNextToken(scip, lpinput) )
         return SCIP_OKAY;
   }
   while( !isNewSection(scip, lpinput) );

   return SCIP_OKAY;
}

/** reads an objective or constraint with name and coefficients */
static
SCIP_RETCODE readCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   SCIP_Bool             isobjective,        /**< indicates whether we are currently reading the coefficients of the objective */
   char*                 name,               /**< pointer to store the name of the line; must be at least of size
                                              *   LP_MAX_LINELEN */
   int*                  coefssize,          /**< size of vars and coefs arrays */
   SCIP_VAR***           vars,               /**< pointer to store the array with variables (must be freed by caller) */
   SCIP_Real**           coefs,              /**< pointer to store the array with coefficients (must be freed by caller) */
   int*                  ncoefs,             /**< pointer to store the number of coefficients */
   SCIP_Bool*            newsection          /**< pointer to store whether a new section was encountered */
   )
{
   SCIP_Bool havesign;
   SCIP_Bool havevalue;
   SCIP_Real coef;
   int coefsign;

   assert(lpinput != NULL);
   assert(name != NULL);
   assert(coefssize != NULL);
   assert(vars != NULL);
   assert(coefs != NULL);
   assert(ncoefs != NULL);
   assert(newsection != NULL);

   *coefssize = 0;
   *vars = NULL;
   *coefs = NULL;
   *name = '\0';
   *ncoefs = 0;
   *newsection = FALSE;

   /* read the first token, which may be the name of the line */
   if( getNextToken(scip, lpinput) )
   {
      /* check if we reached a new section */
      if( isNewSection(scip, lpinput) )
      {
         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* remember the token in the token buffer */
      swapTokenBuffer(lpinput);

      /* get the next token and check, whether it is a colon */
      if( getNextToken(scip, lpinput) )
      {
         if( strcmp(lpinput->token, ":") == 0 )
         {
            /* the second token was a colon: the first token is the line name */
            (void)SCIPmemccpy(name, lpinput->tokenbuf, '\0', LP_MAX_LINELEN);

            name[LP_MAX_LINELEN - 1] = '\0';
            SCIPdebugMsg(scip, "(line %d) read constraint name: '%s'\n", lpinput->linenumber, name);
         }
         else
         {
            /* the second token was no colon: push the tokens back onto the token stack and parse them as coefficients */
            pushToken(lpinput);
            pushBufferToken(lpinput);
         }
      }
      else
      {
         /* there was only one token left: push it back onto the token stack and parse it as coefficient */
         pushBufferToken(lpinput);
      }
   }

   /* initialize buffers for storing the coefficients */
   *coefssize = LP_INIT_COEFSSIZE;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, vars, *coefssize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, coefs, *coefssize) );

   /* read the coefficients */
   coefsign = +1;
   coef = 1.0;
   havesign = FALSE;
   havevalue = FALSE;
   *ncoefs = 0;
   while( getNextToken(scip, lpinput) )
   {
      SCIP_VAR* var;

      /* check if we read a sign */
      if( isSign(lpinput, &coefsign) )
      {
         SCIPdebugMsg(scip, "(line %d) read coefficient sign: %+d\n", lpinput->linenumber, coefsign);
         havesign = TRUE;
         continue;
      }

      /* check if we read a value */
      if( isValue(scip, lpinput, &coef) )
      {
         SCIPdebugMsg(scip, "(line %d) read coefficient value: %g with sign %+d\n", lpinput->linenumber, coef, coefsign);
         if( havevalue )
         {
            syntaxError(scip, lpinput, "two consecutive values.");
            return SCIP_OKAY;
         }
         havevalue = TRUE;
         continue;
      }

      /* check if we reached an equation sense */
      if( isSense(lpinput, NULL) )
      {
         if( isobjective )
         {
            syntaxError(scip, lpinput, "no sense allowed in objective");
            return SCIP_OKAY;
         }

         /* put the sense back onto the token stack */
         pushToken(lpinput);
         break;
      }

      /* check if we reached a new section, that will be only allowed when having no current sign and value and if we
       * are not in the quadratic part
       */
      if( (isobjective || (!havevalue && !havesign)) && isNewSection(scip, lpinput) )
      {
         if( havesign && !havevalue )
         {
            SCIPwarningMessage(scip, "skipped single sign %c without value or variable in objective\n", coefsign == 1 ? '+' : '-');
         }
         else if( isobjective && havevalue && !SCIPisZero(scip, coef) )
         {
            SCIPwarningMessage(scip, "constant term %+g in objective is skipped\n", coef * coefsign);
         }

         *newsection = TRUE;
         return SCIP_OKAY;
      }

      /* check if we start a quadratic part */
      if( *lpinput->token ==  '[' )
      {
         syntaxError(scip, lpinput, "diff reader does not support quadratic objective function.");
         return SCIP_READERROR;
      }

      /* all but the first coefficient need a sign */
      if( *ncoefs > 0 && !havesign )
      {
         syntaxError(scip, lpinput, "expected sign ('+' or '-') or sense ('<' or '>').");
         return SCIP_OKAY;
      }

      /* check if the last variable should be squared */
      if( *lpinput->token == '^' )
      {
         syntaxError(scip, lpinput, "diff reader does not support quadratic objective function.");
         return SCIP_READERROR;
      }
      else
      {
         /* the token is a variable name: get the corresponding variable */
         SCIP_CALL( getVariable(scip, lpinput->token, &var) );
      }

      /* insert the linear coefficient */
      SCIPdebugMsg(scip, "(line %d) read linear coefficient: %+g<%s>\n", lpinput->linenumber, coefsign * coef, SCIPvarGetName(var));
      if( !SCIPisZero(scip, coef) )
      {
         /* resize the vars and coefs array if needed */
         if( *ncoefs >= *coefssize )
         {
            int oldcoefssize;

            oldcoefssize = *coefssize;
            *coefssize *= 2;
            *coefssize = MAX(*coefssize, (*ncoefs)+1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, vars, oldcoefssize, *coefssize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, coefs, oldcoefssize, *coefssize) );
         }
         assert(*ncoefs < *coefssize);

         /* add coefficient */
         (*vars)[*ncoefs] = var;
         (*coefs)[*ncoefs] = coefsign * coef;
         (*ncoefs)++;
      }

      /* reset the flags and coefficient value for the next coefficient */
      coefsign = +1;
      coef = 1.0;
      havesign = FALSE;
      havevalue = FALSE;
   }

   return SCIP_OKAY;
}

/** reads the objective section */
static
SCIP_RETCODE readObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput             /**< LP reading data */
   )
{
   char name[LP_MAX_LINELEN];
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Bool newsection;
   int coefssize;
   int ncoefs;

   assert(lpinput != NULL);

   /* read the objective coefficients */
   SCIP_CALL( readCoefficients(scip, lpinput, TRUE, name, &coefssize, &vars, &coefs, &ncoefs, &newsection) );

   /* change the objective function */
   SCIP_CALL( SCIPchgReoptObjective(scip, lpinput->objsense, vars, coefs, ncoefs) );

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &coefs, coefssize);
   SCIPfreeBlockMemoryArrayNull(scip, &vars, coefssize);

   return SCIP_OKAY;
}

/** reads a diff file */
static
SCIP_RETCODE readDiffFile(
   SCIP*                 scip,               /**< SCIP data structure */
   LPINPUT*              lpinput,            /**< LP reading data */
   const char*           filename            /**< name of the input file */
   )
{
   assert(lpinput != NULL);

   /* open file */
   lpinput->file = SCIPfopen(filename, "r");
   if( lpinput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* free transformed problem */
   if( SCIPisReoptEnabled(scip) && SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
   {
      SCIP_CALL( SCIPfreeReoptSolve(scip) );
   }
   else
   {
      SCIP_CALL( SCIPfreeTransform(scip) );
   }

   /* parse the file */
   lpinput->section = LP_START;
   while( lpinput->section != LP_END && !hasError(lpinput) )
   {
      switch( lpinput->section )
      {
      case LP_START:
         SCIP_CALL( readStart(scip, lpinput) );
         break;

      case LP_OBJECTIVE:
         SCIP_CALL( readObjective(scip, lpinput) );
         break;

      case LP_END: /* this is already handled in the while() loop */
      default:
         SCIPerrorMessage("invalid Diff file section <%d>\n", lpinput->section);
         return SCIP_INVALIDDATA;
      }
   }

   /* close file */
   SCIPfclose(lpinput->file);

   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyDiff)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderDiff(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeDiff)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadDiff)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadDiff(scip, reader, filename, result) );

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the lp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderDiff(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, NULL) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyDiff) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeDiff) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadDiff) );

   return SCIP_OKAY;
}


/** reads problem from file */
SCIP_RETCODE SCIPreadDiff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{  /*lint --e{715}*/
   LPINPUT lpinput;
   int i;

   /* initialize LP input data */
   lpinput.file = NULL;
   lpinput.linebuf[0] = '\0';
   lpinput.probname[0] = '\0';
   lpinput.objname[0] = '\0';
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &lpinput.token, LP_MAX_LINELEN) ); /*lint !e506*/
   lpinput.token[0] = '\0';
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &lpinput.tokenbuf, LP_MAX_LINELEN) ); /*lint !e506*/
   lpinput.tokenbuf[0] = '\0';
   for( i = 0; i < LP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(lpinput.pushedtokens[i]), LP_MAX_LINELEN) );  /*lint !e866 !e506*/
   }

   lpinput.npushedtokens = 0;
   lpinput.linenumber = 0;
   lpinput.linepos = 0;
   lpinput.section = LP_START;
   lpinput.objsense = SCIP_OBJSENSE_MINIMIZE;
   lpinput.haserror = FALSE;
   lpinput.comment = FALSE;
   lpinput.endline = FALSE;

   /* read the file */
   SCIP_CALL( readDiffFile(scip, &lpinput, filename) );

   /* free dynamically allocated memory */
   for( i = 0; i < LP_MAX_PUSHEDTOKENS; ++i )
   {
      SCIPfreeBlockMemoryArray(scip, &lpinput.pushedtokens[i], LP_MAX_LINELEN);
   }
   SCIPfreeBlockMemoryArray(scip, &lpinput.tokenbuf, LP_MAX_LINELEN);
   SCIPfreeBlockMemoryArray(scip, &lpinput.token, LP_MAX_LINELEN);

   /* evaluate the result */
   if( lpinput.haserror )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}
