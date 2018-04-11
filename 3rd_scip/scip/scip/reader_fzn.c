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

/**@file   reader_fzn.c
 * @brief  FlatZinc file reader
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * FlatZinc is a low-level solver input language that is the target language for MiniZinc. It is designed to be easy to
 * translate into the form required by a solver. For more details see http://www.g12.cs.mu.oz.au/minizinc/ .
 *
 * @todo Support more general constraint types
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#ifdef ALLDIFFERENT
#include "scip/cons_alldifferent.h"
#endif
#include "scip/cons_and.h"
#include "scip/cons_cumulative.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_or.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"
#include "scip/pub_misc.h"
#include "scip/reader_fzn.h"

#define READER_NAME             "fznreader"
#define READER_DESC             "file reader for FlatZinc format"
#define READER_EXTENSION        "fzn"


#define FZN_BUFFERLEN         8192      /**< size of the line buffer for reading or writing */
#define FZN_MAX_PUSHEDTOKENS  1

/*
 * Data structures
 */

/** number types */
enum FznNumberType
{
   FZN_BOOL,
   FZN_INT,
   FZN_FLOAT
};
typedef enum FznNumberType FZNNUMBERTYPE;

/** Expression type in FlatZinc File */
enum FznExpType
{
   FZN_EXP_NONE,
   FZN_EXP_UNSIGNED,
   FZN_EXP_SIGNED
};
typedef enum FznExpType FZNEXPTYPE;

/* structures to store the dimension information */
struct Dimensions
{
   int*                  lbs;                /**< lower bounds */
   int*                  ubs;                /**< upper bounds */
   int                   ndims;              /**< number of dimensions */
   int                   size;               /**< size of lbs and ubs */
};
typedef struct Dimensions DIMENSIONS;

/** FlatZinc constant */
struct FznConstant
{
   const char*           name;               /**< constant name */
   FZNNUMBERTYPE         type;               /**< constant type */
   SCIP_Real             value;              /**< constant value */
};
typedef struct FznConstant FZNCONSTANT;

/** structure to store information for an array variable */
struct ConstArray
{
   FZNCONSTANT**         constants;          /**< array of constants */
   char*                 name;               /**< name of constant array */
   int                   nconstants;         /**< number of constants */
   FZNNUMBERTYPE         type;               /**< constant type */
};
typedef struct ConstArray CONSTARRAY;

/** structure to store information for an array variable */
struct VarArray
{
   SCIP_VAR**            vars;               /**< variable belonging to the variable array */
   char*                 name;               /**< name of the array variable */
   DIMENSIONS*           info;               /**< dimension information */
   int                   nvars;              /**< number of variables */
   FZNNUMBERTYPE         type;               /**< variable type */
};
typedef struct VarArray VARARRAY;

/** data for FlatZinc reader */
struct SCIP_ReaderData
{
   VARARRAY**            vararrays;          /**< variable arrays to output */
   int                   nvararrays;         /**< number of variables */
   int                   vararrayssize;      /**< size of variable array */
};

/** tries to creates and adds a constraint; sets parameter created to TRUE if method was successful
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - fzninput,       : FZN reading data
 *  - fname,          : functions identifier name
 *  - ftokens,        : function identifier tokens
 *  - nftokens,       : number of function identifier tokes
 *
 *  output
 *  - created         : pointer to store whether a constraint was created or not
 */
#define CREATE_CONSTRAINT(x) SCIP_RETCODE x (SCIP* scip, FZNINPUT* fzninput, const char* fname, char** ftokens, int nftokens, SCIP_Bool* created)


/** FlatZinc reading data */
struct FznInput
{
   SCIP_FILE*            file;
   SCIP_HASHTABLE*       varHashtable;
   SCIP_HASHTABLE*       constantHashtable;
   FZNCONSTANT**         constants;
   char                  linebuf[FZN_BUFFERLEN+1];
   char*                 token;
   char*                 pushedtokens[FZN_MAX_PUSHEDTOKENS];
   int                   npushedtokens;
   int                   linenumber;
   int                   linepos;
   int                   bufpos;
   int                   nconstants;
   int                   sconstants;
   SCIP_OBJSENSE         objsense;
   SCIP_Bool             hasdot;             /**< if the current token is a number, this bool tells if it contains a dot */
   SCIP_Bool             endline;            /**< current buffer contains everything until the line ends */
   SCIP_Bool             comment;            /**< current buffer contains everything until a comment starts */
   SCIP_Bool             haserror;           /**< a error was detected during parsing */
   SCIP_Bool             valid;
   SCIP_Bool             initialconss;       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss;       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamiccols;        /**< should columns be added and removed dynamically to the LP? */
   SCIP_Bool             dynamicrows;        /**< should rows be added and removed dynamically to the LP? */

   VARARRAY**            vararrays;          /**< variable arrays */
   int                   nvararrays;         /**< number of variables */
   int                   vararrayssize;      /**< size of variable array */

   CONSTARRAY**          constarrays;        /**< variable arrays to output */
   int                   nconstarrays;       /**< number of variables */
   int                   constarrayssize;    /**< size of variable array */
};
typedef struct FznInput FZNINPUT;

/** FlatZinc writing data */
struct FznOutput
{
   char*                 varbuffer;          /* buffer for auxiliary variables (float representatives of discrete variables) */
   int                   varbufferlen;       /* current length of the above buffer */
   int                   varbufferpos;       /* current filling position in the above buffer */
   char*                 castbuffer;         /* buffer for int2float conversion constraints */
   int                   castbufferlen;      /* current length of the above buffer */
   int                   castbufferpos;      /* current filling position in the above buffer */
   char*                 consbuffer;         /* buffer for all problem constraints */
   int                   consbufferlen;      /* current length of the above buffer */
   int                   consbufferpos;      /* current filling position in the above buffer */
   int                   ndiscretevars;      /* number of discrete variables in the problem */
   SCIP_Bool*            varhasfloat;        /* array which indicates, whether a discrete variable already has a float representative */
};
typedef struct FznOutput FZNOUTPUT;

static const char delimchars[] = " \f\n\r\t\v";
static const char tokenchars[] = ":<>=;{}[],()";
static const char commentchars[] = "%";

/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVar)
{  /*lint --e{715}*/
   SCIP_VAR* var = (SCIP_VAR*) elem;

   assert(var != NULL);
   return (void*) SCIPvarGetName(var);
}

/** gets the key (i.e. the name) of the flatzinc constant */
static
SCIP_DECL_HASHGETKEY(hashGetKeyConstant)
{  /*lint --e{715}*/
   FZNCONSTANT* constant = (FZNCONSTANT*) elem;

   assert(constant != NULL);
   return (void*) constant->name;
}

/** comparison method for sorting variable arrays  w.r.t. to their name */
static
SCIP_DECL_SORTPTRCOMP(vararraysComp)
{
   return strcmp( ((VARARRAY*)elem1)->name, ((VARARRAY*)elem2)->name );
}


/** frees a given buffer char* array */
static
void freeStringBufferArray(
   SCIP*                 scip,               /**< SCIP data structure */
   char**                array,              /**< buffer array to free */
   int                   nelements           /**< number of elements */
   )
{
   int i;

   for( i = 0; i < nelements; ++i )
   {
      SCIPfreeBufferArray(scip, &array[i]);
   }

   SCIPfreeBufferArray(scip, &array);
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

/** check if the current token is equal to give char */
static
SCIP_Bool isChar(
   const char*           token,              /**< token to be checked */
   char                  c                   /**< char to compare */
   )
{
   if( strlen(token) == 1 && *token == c )
      return TRUE;

   return FALSE;
}

/** check if the current token is Bool expression, this means false or true */
static
SCIP_Bool isBoolExp(
   const char*           name,               /**< name to check */
   SCIP_Bool*            value               /**< pointer to store the Bool value */
   )
{
   /* check the name */
   if( strlen(name) == 4 && strncmp(name, "true", 4) == 0 )
   {
      *value = TRUE;
      return TRUE;
   }
   else if( strlen(name) == 1 && strncmp(name, "1", 1) == 0 )
   {
      /* we also allow 1 as true */
      *value = TRUE;
      return TRUE;
   }
   else if( strlen(name) == 5 && strncmp(name, "false", 5) == 0 )
   {
      *value = FALSE;
      return TRUE;
   }
   else if( strlen(name) == 1 && strncmp(name, "0", 1) == 0 )
   {
      /* we also allow 0 as false */
      *value = FALSE;
      return TRUE;
   }

   return FALSE;
}


/** check if the current token is an identifier, this means [A-Za-z][A-Za-z0-9_]* */
static
SCIP_Bool isIdentifier(
   const char*           name                /**< name to check */
   )
{
   int i;

   /* check if the identifier starts with a letter */
   if( strlen(name) == 0 || !isalpha((unsigned char)name[0]) )
      return FALSE;

   i = 1;
   while( name[i] )
   {
      if( !isalnum((unsigned char)name[i]) && name[i] != '_' )
         return FALSE;
      i++;
   }

   return TRUE;
}

/** returns whether the current character is member of a value string */
static
SCIP_Bool isValueChar(
   char                  c,                  /**< input character */
   char                  nextc,              /**< next input character */
   SCIP_Bool             firstchar,          /**< is the given character the first char of the token? */
   SCIP_Bool*            hasdot,             /**< pointer to update the dot flag */
   FZNEXPTYPE*           exptype             /**< pointer to update the exponent type */
   )
{
   assert(hasdot != NULL);
   assert(exptype != NULL);

   if( isdigit((unsigned char)c) )
      return TRUE;
   else if( firstchar && (c == '+' || c == '-') )
      return TRUE;
   else if( (*exptype == FZN_EXP_NONE) && !(*hasdot) && (c == '.') && (isdigit((unsigned char)nextc)))
   {
      *hasdot = TRUE;
      return TRUE;
   }
   else if( !firstchar && (*exptype == FZN_EXP_NONE) && (c == 'e' || c == 'E') )
   {
      if( nextc == '+' || nextc == '-' )
      {
         *exptype = FZN_EXP_SIGNED;
         return TRUE;
      }
      else if( isdigit((unsigned char)nextc) )
      {
         *exptype = FZN_EXP_UNSIGNED;
         return TRUE;
      }
   }
   else if( (*exptype == FZN_EXP_SIGNED) && (c == '+' || c == '-') )
   {
      *exptype = FZN_EXP_UNSIGNED;
      return TRUE;
   }

   return FALSE;
}

/** compares two token if they are equal */
static
SCIP_Bool equalTokens(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           token1,             /**< first token */
   const char*           token2              /**< second token */
   )
{
   assert(token1 != NULL);
   assert(token2 != NULL);

   if( strlen(token1) != strlen(token2) )
      return FALSE;

   return !strncmp(token1, token2, strlen(token2) );
}

/** reads the next line from the input file into the line buffer; skips comments;
 *  returns whether a line could be read
 */
static
SCIP_Bool getNextLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   int i;

   assert(fzninput != NULL);

   /* if we previously detected a comment we have to parse the remaining line away if there is something left */
   if( !fzninput->endline && fzninput->comment )
   {
      SCIPdebugMsg(scip, "Throwing rest of comment away.\n");

      do
      {
         fzninput->linebuf[FZN_BUFFERLEN-2] = '\0';
         (void)SCIPfgets(fzninput->linebuf, (int) sizeof(fzninput->linebuf), fzninput->file);
      }
      while( fzninput->linebuf[FZN_BUFFERLEN-2] != '\0' );

      fzninput->comment = FALSE;
      fzninput->endline = TRUE;
   }

   /* clear the line */
   BMSclearMemoryArray(fzninput->linebuf, FZN_BUFFERLEN);
   fzninput->linebuf[FZN_BUFFERLEN-2] = '\0';

   /* set line position */
   if( fzninput->endline )
   {
      fzninput->linepos = 0;
      fzninput->linenumber++;
   }
   else
      fzninput->linepos += FZN_BUFFERLEN - 2;

   if( SCIPfgets(fzninput->linebuf, (int) sizeof(fzninput->linebuf), fzninput->file) == NULL )
      return FALSE;

   fzninput->bufpos = 0;

   if( fzninput->linebuf[FZN_BUFFERLEN-2] != '\0' )
   {
      char* last;

      /* buffer is full; erase last token since it might be incomplete */
      fzninput->endline = FALSE;
      last = strrchr(fzninput->linebuf, ' ');

      if( last == NULL )
      {
         SCIPwarningMessage(scip, "we read %d characters from the file; this might indicate a corrupted input file!\n",
            FZN_BUFFERLEN - 2);
         fzninput->linebuf[FZN_BUFFERLEN-2] = '\0';
         SCIPdebugMsg(scip, "the buffer might be corrupted\n");
      }
      else
      {
         SCIPfseek(fzninput->file, -(long) strlen(last), SEEK_CUR);
         SCIPdebugMsg(scip, "correct buffer, reread the last %ld characters\n", (long) strlen(last));
         *last = '\0';
      }
   }
   else
   {
      /* found end of line */
      fzninput->endline = TRUE;
   }

   fzninput->linebuf[FZN_BUFFERLEN-1] = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
   fzninput->comment = FALSE;

   /* skip characters after comment symbol */
   for( i = 0; commentchars[i] != '\0'; ++i )
   {
      char* commentstart;

      commentstart = strchr(fzninput->linebuf, commentchars[i]);
      if( commentstart != NULL )
      {
         *commentstart = '\0';
         *(commentstart+1) = '\0'; /* we want to use lookahead of one char -> we need two \0 at the end */
         fzninput->comment = TRUE;
         break;
      }
   }

   return TRUE;
}


/** reads the next token from the input file into the token buffer; returns whether a token was read */
static
SCIP_Bool getNextToken(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   SCIP_Bool hasdot;
   FZNEXPTYPE exptype;
   char* buf;
   int tokenlen;

   assert(fzninput != NULL);
   assert(fzninput->bufpos < FZN_BUFFERLEN);

   /* if the current line got marked as comment get the next line */
   if( fzninput->comment && !getNextLine(scip, fzninput) )
   {
      SCIPdebugMsg(scip, "(line %d) end of file\n", fzninput->linenumber);
      return FALSE;
   }

   /* check the token stack */
   if( fzninput->npushedtokens > 0 )
   {
      SCIPswapPointers((void**)&fzninput->token, (void**)&fzninput->pushedtokens[fzninput->npushedtokens-1]);
      fzninput->npushedtokens--;
      SCIPdebugMsg(scip, "(line %d) read token again: '%s'\n", fzninput->linenumber, fzninput->token);
      return TRUE;
   }

   /* skip delimiters */
   buf = fzninput->linebuf;
   while( isDelimChar(buf[fzninput->bufpos]) )
   {
      if( buf[fzninput->bufpos] == '\0' )
      {
         if( !getNextLine(scip, fzninput) )
         {
            SCIPdebugMsg(scip, "(line %d) end of file\n", fzninput->linenumber);
            return FALSE;
         }
         assert(fzninput->bufpos == 0);
      }
      else
      {
         fzninput->bufpos++;
         fzninput->linepos++;
      }
   }
   assert(fzninput->bufpos < FZN_BUFFERLEN);
   assert(!isDelimChar(buf[fzninput->bufpos]));

   hasdot = FALSE;
   exptype = FZN_EXP_NONE;

   if( buf[fzninput->bufpos] == '.' && buf[fzninput->bufpos+1] == '.')
   {
      /* found <..> which only occurs in Ranges and is a "keyword" */
      tokenlen = 2;
      fzninput->bufpos += 2;
      fzninput->linepos += 2;
      fzninput->token[0] = '.';
      fzninput->token[1] = '.';
   }
   else if( isValueChar(buf[fzninput->bufpos], buf[fzninput->bufpos+1], TRUE, &hasdot, &exptype) )
   {
      /* read value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < FZN_BUFFERLEN);
         assert(!isDelimChar(buf[fzninput->bufpos]));
         fzninput->token[tokenlen] = buf[fzninput->bufpos];
         tokenlen++;
         fzninput->bufpos++;
         fzninput->linepos++;
         assert(fzninput->bufpos < FZN_BUFFERLEN);
      }
      while( isValueChar(buf[fzninput->bufpos], buf[fzninput->bufpos+1], FALSE, &hasdot, &exptype) );

      fzninput->hasdot = hasdot;
   }
   else
   {
      /* read non-value token */
      tokenlen = 0;
      do
      {
         assert(tokenlen < FZN_BUFFERLEN);
         fzninput->token[tokenlen] = buf[fzninput->bufpos];
         tokenlen++;
         fzninput->bufpos++;
         fzninput->linepos++;

         /* check for annotations */
         if(tokenlen == 1 && fzninput->token[0] == ':' && buf[fzninput->bufpos] == ':')
         {
            fzninput->token[tokenlen] = buf[fzninput->bufpos];
            tokenlen++;
            fzninput->bufpos++;
            fzninput->linepos++;
            break;
         }

         if( tokenlen == 1 && isTokenChar(fzninput->token[0]) )
            break;
      }
      while( !isDelimChar(buf[fzninput->bufpos]) && !isTokenChar(buf[fzninput->bufpos]) );
   }

   assert(tokenlen < FZN_BUFFERLEN);
   fzninput->token[tokenlen] = '\0';

   SCIPdebugMsg(scip, "(line %d) read token: '%s'\n", fzninput->linenumber, fzninput->token);

   return TRUE;
}

/** puts the current token on the token stack, such that it is read at the next call to getNextToken() */
static
void pushToken(
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   assert(fzninput != NULL);
   assert(fzninput->npushedtokens < FZN_MAX_PUSHEDTOKENS);

   SCIPswapPointers((void**)&fzninput->pushedtokens[fzninput->npushedtokens], (void**)&fzninput->token);
   fzninput->npushedtokens++;
}

/** checks whether the current token is a semicolon which closes a statement */
static
SCIP_Bool isEndStatement(
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   assert(fzninput != NULL);

   return isChar(fzninput->token, ';');
}

/** returns whether the current token is a value */
static
SCIP_Bool isValue(
   const char*           token,              /**< token to check */
   SCIP_Real*            value               /**< pointer to store the value (unchanged, if token is no value) */
   )
{
   double val;
   char* endptr;

   assert(value != NULL);

   val = strtod(token, &endptr);
   if( endptr != token && *endptr == '\0' )
   {
      *value = val;
      return TRUE;
   }

   return FALSE;
}

/*
 * Local methods (for reading)
 */

/** issues an error message and marks the FlatZinc data to have errors */
static
void syntaxError(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           msg                 /**< error message */
   )
{
   assert(fzninput != NULL);

   SCIPerrorMessage("Syntax error in line %d: %s found <%s>\n", fzninput->linenumber, msg, fzninput->token);
   SCIPerrorMessage("  input: %s\n", fzninput->linebuf);

   fzninput->haserror = TRUE;
}

/** returns whether a syntax error was detected */
static
SCIP_Bool hasError(
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   assert(fzninput != NULL);

   return (fzninput->haserror || !fzninput->valid);
}

/** create reader data */
static
SCIP_RETCODE readerdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA**     readerdata          /**< pointer to reader data */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, readerdata) );

   (*readerdata)->vararrays = NULL;
   (*readerdata)->nvararrays = 0;
   (*readerdata)->vararrayssize = 0;

   return SCIP_OKAY;
}

/** ensure the size if the variable array */
static
SCIP_RETCODE ensureVararrySize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   int nvararrays;
   int vararrayssize;

   nvararrays = readerdata->nvararrays;
   vararrayssize = readerdata->vararrayssize;

   if( vararrayssize == nvararrays )
   {
      if( vararrayssize == 0 )
      {
         vararrayssize = 100;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->vararrays, vararrayssize) );
      }
      else
      {
         vararrayssize *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &readerdata->vararrays, readerdata->vararrayssize, vararrayssize) );
      }
   }

   readerdata->vararrayssize = vararrayssize;

   return SCIP_OKAY;
}

/** ensure the size if the variable array */
static
SCIP_RETCODE ensureVararrySizeFznInput(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   int nvararrays;
   int vararrayssize;

   nvararrays = fzninput->nvararrays;
   vararrayssize = fzninput->vararrayssize;

   if( vararrayssize == nvararrays )
   {
      if( vararrayssize == 0 )
      {
         vararrayssize = 100;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &fzninput->vararrays, vararrayssize) );
      }
      else
      {
         vararrayssize *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &fzninput->vararrays, fzninput->vararrayssize, vararrayssize) );
      }
   }

   fzninput->vararrayssize = vararrayssize;

   return SCIP_OKAY;
}

/** ensure the size if the variable array */
static
SCIP_RETCODE ensureConstarrySizeFznInput(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   int nconstarrays;
   int constarrayssize;

   nconstarrays = fzninput->nconstarrays;
   constarrayssize = fzninput->constarrayssize;

   if( constarrayssize == nconstarrays )
   {
      if( constarrayssize == 0 )
      {
         constarrayssize = 100;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &fzninput->constarrays, constarrayssize) );
      }
      else
      {
         constarrayssize *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &fzninput->constarrays, fzninput->constarrayssize, constarrayssize) );
      }
   }

   fzninput->constarrayssize = constarrayssize;

   return SCIP_OKAY;
}

/** print given value in FlatZinc format to given stream */
static
void printValue(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             value,              /**< value to print */
   FZNNUMBERTYPE         type                /**< FlatZinc number type */
   )
{
   switch( type )
   {
   case FZN_BOOL:
      if( value < 0.5 )
         SCIPinfoMessage(scip, file, "false");
      else
         SCIPinfoMessage(scip, file, "true");
      break;
   case FZN_INT:
   {
      SCIP_Longint longvalue;
      longvalue = SCIPconvertRealToLongint(scip, value);
      SCIPinfoMessage(scip, file, "%" SCIP_LONGINT_FORMAT "", longvalue);
      break;
   }
   case FZN_FLOAT:
      if( SCIPisIntegral(scip, value) )
      {
         printValue(scip, file, value, FZN_INT);

         /* add a ".0" to be type save */
         SCIPinfoMessage(scip, file, ".0");
      }
      else
      {
         SCIPinfoMessage(scip, file, "%.1f", value);
      }
      break;
   }
}

/*
 * Local methods (for VARARRAY)
 */

/** free dimension structure */
static
SCIP_RETCODE copyDimensions(
   SCIP*                 scip,               /**< SCIP data structure */
   DIMENSIONS**          target,             /**< pointer to dimension target structure */
   DIMENSIONS*           source              /**< dimension source */
   )
{
   if( source != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, target) );

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*target)->lbs, source->lbs, source->ndims) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*target)->ubs, source->ubs, source->ndims) );
      (*target)->ndims = source->ndims;
      (*target)->size = source->ndims;
   }
   else
      *target = NULL;

   return SCIP_OKAY;
}

/** create variable array data structure */
static
SCIP_RETCODE createVararray(
   SCIP*                 scip,               /**< SCIP data structure */
   VARARRAY**            vararray,           /**< pointer to variable array */
   const char*           name,               /**< name of the variable array */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   nvars,              /**< number of variables */
   FZNNUMBERTYPE         type,               /**< variable type */
   DIMENSIONS*           info                /**< dimension information for output */
   )
{
   /* allocate memory for the new vararray struct */
   SCIP_CALL( SCIPallocBlockMemory(scip, vararray) );

   /* copy variable pointers */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*vararray)->vars, vars, nvars) );

   /* copy variable array name */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*vararray)->name, name, strlen(name)+1) );

   SCIP_CALL( copyDimensions(scip, &(*vararray)->info, info) );

   (*vararray)->nvars = nvars;
   (*vararray)->type = type;

   return SCIP_OKAY;
}

/** free dimension structure */
static
void freeDimensions(
   SCIP*                 scip,               /**< SCIP data structure */
   DIMENSIONS**          dim                 /**< pointer to dimension structure */
   )
{
   if( *dim != NULL )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*dim)->lbs, (*dim)->size);
      SCIPfreeBlockMemoryArrayNull(scip, &(*dim)->ubs, (*dim)->size);
      SCIPfreeBlockMemory(scip, dim);
   }
}

/** free variable array data structure */
static
void freeVararray(
   SCIP*                 scip,               /**< SCIP data structure */
   VARARRAY**            vararray            /**< pointer to variable array */
   )
{
   freeDimensions(scip, &(*vararray)->info);

   SCIPfreeBlockMemoryArray(scip, &(*vararray)->name, strlen((*vararray)->name) + 1);
   SCIPfreeBlockMemoryArray(scip, &(*vararray)->vars, (*vararray)->nvars);

   SCIPfreeBlockMemory(scip, vararray);
}

/** searches the variable array data base if a constant array exists with the given name; if it exists it is returned */
static
VARARRAY* findVararray(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name                /**< variable array name */
   )
{
   VARARRAY* vararray;
   int c;

   /* search in constants array list for a constants array with the given name */
   for( c = 0; c < fzninput->nvararrays; ++c )
   {
      vararray = fzninput->vararrays[c];

      if( equalTokens(scip, name, vararray->name) )
         return vararray;
   }

   return NULL;
}

/*
 * Local methods (for CONSTARRAY)
 */

/** create constant array data structure */
static
SCIP_RETCODE createConstarray(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTARRAY**          constarray,         /**< pointer to constant array */
   const char*           name,               /**< name of the variable array */
   FZNCONSTANT**         constants,          /**< array of constants */
   int                   nconstants,         /**< number of constants */
   FZNNUMBERTYPE         type                /**< constant type */
   )
{
   SCIPdebugMsg(scip, "create constant array <%s>\n", name);

   /* allocate memory for the new constarray struct */
   SCIP_CALL( SCIPallocBlockMemory(scip, constarray) );

   /* copy constant values */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*constarray)->constants, constants, nconstants) );

   /* copy constant array name */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*constarray)->name, name, strlen(name)+1) );

   (*constarray)->nconstants = nconstants;
   (*constarray)->type = type;

   return SCIP_OKAY;
}

/** free constant array data structure */
static
void freeConstarray(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTARRAY**          constarray          /**< pointer to constant array */
   )
{
   SCIPdebugMsg(scip, "free constant array <%s>\n", (*constarray)->name);

   /* free variable pointers */
   SCIPfreeBlockMemoryArray(scip, &(*constarray)->constants, (*constarray)->nconstants);

   /* free variable array name */
   SCIPfreeBlockMemoryArray(scip, &(*constarray)->name, strlen((*constarray)->name) + 1);

   /* allocate memory for the new vararray struct */
   SCIPfreeBlockMemory(scip, constarray);
}

/** searches the constant array data base if a constant array exists with the given name; if it exists it is returned */
static
CONSTARRAY* findConstarray(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name                /**< constant array name */
   )
{
   CONSTARRAY* constarray;
   int c;

   /* search in constants array list for a constants array with the given name */
   for( c = 0; c < fzninput->nconstarrays; ++c )
   {
      constarray = fzninput->constarrays[c];

      if( equalTokens(scip, name, constarray->name) )
         return constarray;
   }

   return NULL;
}

/** add variable to the reader data */
static
SCIP_RETCODE readerdataAddOutputvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   SCIP_VAR*             var,                /**< variable to add to the reader data */
   FZNNUMBERTYPE         type                /**< variable type */
   )
{
   DIMENSIONS* info;
   const char* name;
   VARARRAY* vararray;
   int nvararrays;

   nvararrays = readerdata->nvararrays;

   SCIP_CALL( ensureVararrySize(scip, readerdata) );
   assert(nvararrays < readerdata->vararrayssize);

   /* get variable name */
   name = SCIPvarGetName(var);

   /* allocate memory for the new vararray struct */
   SCIP_CALL( SCIPallocBlockMemory(scip, &vararray) );

   /* copy variable pointers */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vararray->vars, &var, 1) );

   /* copy variable array name */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vararray->name, name, strlen(name)+1) );

   SCIP_CALL( SCIPallocBlockMemory(scip, &info) );
   info->lbs = NULL;
   info->ubs = NULL;
   info->ndims = 0;
   info->size = 0;

   vararray->info = info;
   vararray->nvars = 1;
   vararray->type = type;

   readerdata->vararrays[nvararrays] = vararray;
   readerdata->nvararrays++;

   return SCIP_OKAY;
}

/** add variable to the reader data */
static
SCIP_RETCODE readerdataAddOutputvararray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   const char*           name,               /**< name of the variable array */
   SCIP_VAR**            vars,               /**< array of variable to add to the reader data */
   int                   nvars,              /**< number of variables */
   FZNNUMBERTYPE         type,               /**< variable type */
   DIMENSIONS*           info                /**< dimension information for output */
   )
{
   VARARRAY* vararray;
   int nvararrays;

   nvararrays = readerdata->nvararrays;

   SCIP_CALL( ensureVararrySize(scip, readerdata) );
   assert(nvararrays < readerdata->vararrayssize);

   /* create variable array data structure */
   SCIP_CALL( createVararray(scip, &vararray, name, vars, nvars, type, info) );

   readerdata->vararrays[nvararrays] = vararray;
   readerdata->nvararrays++;

   return SCIP_OKAY;
}

/** add variable to the input data */
static
SCIP_RETCODE fzninputAddVararray(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< name of the variable array */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   nvars,              /**< number of variables */
   FZNNUMBERTYPE         type,               /**< variable type */
   DIMENSIONS*           info                /**< dimension information for output */
   )
{
   VARARRAY* vararray;
   int nvararrays;

   nvararrays = fzninput->nvararrays;

   SCIP_CALL( ensureVararrySizeFznInput(scip, fzninput) );
   assert(nvararrays < fzninput->vararrayssize);

   /* create variable array data structure */
   SCIP_CALL( createVararray(scip, &vararray, name, vars, nvars, type, info) );

   fzninput->vararrays[nvararrays] = vararray;
   fzninput->nvararrays++;

   return SCIP_OKAY;
}

/** add variable to the reader data */
static
SCIP_RETCODE fzninputAddConstarray(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< name of the variable array */
   FZNCONSTANT**         constants,          /**< array of constants */
   int                   nconstants,         /**< number of constants */
   FZNNUMBERTYPE         type                /**< variable type */
   )
{
   CONSTARRAY* constarray;
   int nconstarrays;

   nconstarrays = fzninput->nconstarrays;

   SCIP_CALL( ensureConstarrySizeFznInput(scip, fzninput) );
   assert(nconstarrays < fzninput->constarrayssize);

   /* create constant array structure */
   SCIP_CALL( createConstarray(scip, &constarray, name, constants, nconstants, type) );

   fzninput->constarrays[nconstarrays] = constarray;
   fzninput->nconstarrays++;

   return SCIP_OKAY;
}

/** creates, adds, and releases a linear constraint */
static
SCIP_RETCODE createQuadraticCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation (u) */
   SCIP_Bool             initialconss,       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss,       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamicrows         /**< should rows be added and removed dynamically to the LP? */
   )
{
   SCIP_CONS* cons;

   SCIP_CALL( SCIPcreateConsQuadratic(scip, &cons, name, nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoefs, lhs, rhs,
         initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows) );

   SCIPdebugPrintCons(scip, cons, NULL);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** creates, adds, and releases a linear constraint */
static
SCIP_RETCODE createLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initialconss,       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss,       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamicrows         /**< should rows be added and removed dynamically to the LP? */
   )
{
   SCIP_CONS* cons;

   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nvars, vars, vals, lhs, rhs,
         initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, dynamicconss, dynamicrows, FALSE) );

   SCIPdebugPrintCons(scip, cons, NULL);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/** create a linking between the two given identifiers */
static
SCIP_RETCODE createLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           consname,           /**< name of constraint */
   const char*           name1,              /**< name of first identifier */
   const char*           name2,              /**< name of second identifier */
   SCIP_Real             lhs,                /**< left hand side of the linking */
   SCIP_Real             rhs                 /**< right hand side of the linking */
   )
{
   SCIP_VAR** vars;
   SCIP_Real vals[] = {0.0,0.0};
   SCIP_Real value1;
   SCIP_Real value2;
   int nvars;

   nvars = 0;
   value1 = 0.0;
   value2 = 0.0;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );

   vars[nvars] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) name1);
   if( vars[nvars] != NULL )
   {
      vals[nvars] = 1.0;
      nvars++;
   }
   else if( !isValue(name1, &value1) )
   {
      FZNCONSTANT* constant;

      constant = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, (char*) name1);
      assert(constant != NULL);

      value1 = constant->value;
   }

   vars[nvars] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) name2);
   if( vars[nvars] != NULL )
   {
      vals[nvars] = -1.0;
      nvars++;
   }
   else if( !isValue(name2, &value2) )
   {
      FZNCONSTANT* constant;

      constant = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, (char*) name2);
      assert(constant != NULL);

      value2 = constant->value;
   }

   if( !SCIPisInfinity(scip, -lhs) )
      lhs += (value2 - value1);

   if( !SCIPisInfinity(scip, rhs) )
      rhs += (value2 - value1);

   SCIP_CALL( createLinearCons(scip, consname, nvars, vars, vals, lhs, rhs, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** parse array index expression */
static
void parseArrayIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   int*                  idx                 /**< pointer to store the array index */
   )
{
   SCIP_Real value;

   assert( isChar(fzninput->token, '[') );

   /* parse array index expression */
   if( !getNextToken(scip, fzninput) || isEndStatement(fzninput) )
   {
      syntaxError(scip, fzninput, "expecting array index expression");
      return;
   }

   if( isIdentifier(fzninput->token) )
   {
      FZNCONSTANT* constant;

      /* identifier has to be one of a constant */
      constant = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, fzninput->token);

      if( constant == NULL )
         syntaxError(scip, fzninput, "unknown index name");
      else
      {
         assert(constant->type == FZN_INT);
         *idx = (int) constant->value;
      }
   }
   else if( isValue(fzninput->token, &value) )
   {
      assert( fzninput->hasdot == FALSE );
      *idx = (int) value;
   }
   else
      syntaxError(scip, fzninput, "expecting array index expression");
}

/** unroll assignment if it is an array access one */
static
void flattenAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   char*                 assignment          /**< assignment to unroll */
   )
{
   assert(scip != NULL);
   assert(fzninput != NULL);

   SCIPdebugMsg(scip, "parse assignment expression\n");

   if( !getNextToken(scip, fzninput) || isEndStatement(fzninput) )
   {
      syntaxError(scip, fzninput, "expecting more tokens");
      return;
   }

   if( isIdentifier(fzninput->token) )
   {
      char name[FZN_BUFFERLEN];
      int idx;

      (void) SCIPsnprintf(name, FZN_BUFFERLEN, "%s", fzninput->token);

      if( !getNextToken(scip, fzninput) )
      {
         syntaxError(scip, fzninput, "expecting at least a semicolon to close the statement");
         return;
      }

      /* check if it is an array access expression */
      if( isChar(fzninput->token, '[') )
      {
         idx = -1;
         parseArrayIndex(scip, fzninput, &idx);

         assert(idx >= 0);

         if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, ']') )
         {
            syntaxError(scip, fzninput, "expecting token <]>");
            return;
         }

         /* put constant name or variable name together */
         (void) SCIPsnprintf(assignment, FZN_BUFFERLEN, "%s[%d]", name, idx);
      }
      else
      {
         (void) SCIPsnprintf(assignment, FZN_BUFFERLEN, "%s", name);

         /* push the current token back for latter evaluations */
         pushToken(fzninput);
      }
   }
   else
      (void) SCIPsnprintf(assignment, FZN_BUFFERLEN, "%s", fzninput->token);
}

/** computes w.r.t. to the given side value and relation the left and right side for a SCIP linear constraint */
static
void computeLinearConsSides(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< name of the relation */
   SCIP_Real             sidevalue,          /**< parsed side value */
   SCIP_Real*            lhs,                /**< pointer to left hand side */
   SCIP_Real*            rhs                 /**< pointer to right hand side */
   )
{
   SCIPdebugMsg(scip, "check relation <%s>\n", name);

   /* compute left and right hand side of the linear constraint */
   if( equalTokens(scip, name, "eq") )
   {
      *lhs = sidevalue;
      *rhs = sidevalue;
   }
   else if( equalTokens(scip, name, "ge") )
   {
      *lhs = sidevalue;
   }
   else if( equalTokens(scip, name, "le") )
   {
      *rhs = sidevalue;
   }
   else if( equalTokens(scip, name, "gt") )
   {
      /* greater than only works if there are not continuous variables are involved */
      *lhs = sidevalue + 1.0;
   }
   else if( equalTokens(scip, name, "lt") )
   {
      /* less than only works if there are not continuous variables are involved */
      *rhs = sidevalue - 1.0;
   }
   else
      syntaxError(scip, fzninput, "unknown relation in constraint identifier name");

   SCIPdebugMsg(scip, "lhs = %g, rhs = %g\n", *lhs, *rhs);
}

/** parse a list of elements which is separates by a comma */
static
SCIP_RETCODE parseList(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   char***               elements,           /**< pointer to char* array for storing the elements of the list */
   int*                  nelements,          /**< pointer to store the number of elements */
   int                   selements           /**< size of the elements char* array */
   )
{
   char assignment[FZN_BUFFERLEN];
   assert(selements > 0);

   /* check if the list is not empty */
   if( getNextToken(scip, fzninput) && !isChar(fzninput->token, ']') )
   {
      /* push back token */
      pushToken(fzninput);

      /* loop through the array */
      do
      {
         if(selements == *nelements)
         {
            selements *= 2;
            SCIP_CALL( SCIPreallocBufferArray(scip, elements, selements) );
         }

         /* parse and flatten assignment */
         flattenAssignment(scip, fzninput, assignment);

         if( hasError(fzninput) )
            break;

         /* store assignment */
         SCIP_CALL( SCIPduplicateBufferArray(scip, &(*elements)[(*nelements)], assignment, (int) strlen(assignment) + 1) ); /*lint !e866*/

         (*nelements)++;
      }
      while( getNextToken(scip, fzninput) && isChar(fzninput->token, ',') );
   }
   else
   {
      SCIPdebugMsg(scip, "list is empty\n");
   }


   /* push back ']' which closes the list */
   pushToken(fzninput);

   return SCIP_OKAY;
}

/** parse range expression */
static
void parseRange(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   FZNNUMBERTYPE*        type,               /**< pointer to store the number type */
   SCIP_Real*            lb,                 /**< pointer to store the lower bound */
   SCIP_Real*            ub                  /**< pointer to store the upper bound */
   )
{
   if( !getNextToken(scip, fzninput) )
   {
      syntaxError(scip, fzninput, "expected left side of range");
      return;
   }

   /* current token should be the lower bound */
   if( !isValue(fzninput->token, lb) )
   {
      syntaxError(scip, fzninput, "expected lower bound value");
      return;
   }

   /* check if we have a float notation or an integer notation which defines the type of the variable */
   if( fzninput->hasdot || !SCIPisIntegral(scip, *lb) )
      *type = FZN_FLOAT;
   else
      *type = FZN_INT;

   /* parse next token which should be <..> */
   if( !getNextToken(scip, fzninput) || !equalTokens(scip, fzninput->token, "..") )
   {
      syntaxError(scip, fzninput, "expected <..>");
      return;
   }

   /* parse upper bound */
   if( !getNextToken(scip, fzninput) || !isValue(fzninput->token, ub) )
   {
      syntaxError(scip, fzninput, "expected upper bound value");
      return;
   }

   /* check if upper bound notation fits which lower bound notation */
   if( fzninput->hasdot != (*type == FZN_FLOAT) )
   {
      SCIPwarningMessage(scip, "lower bound and upper bound mismatch in value type, assume %s variable type\n",
         fzninput->hasdot ? "an integer" : "a continuous");
   }
}

/** parse dimension information */
static
SCIP_RETCODE parseOutputDimensioninfo(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   DIMENSIONS**          info                /**< pointer to store the output dimension information if one */
   )
{
   FZNNUMBERTYPE type;
   SCIP_Real lb;
   SCIP_Real ub;
   int nelements;
   int size;

   nelements = 0;
   size = 100;

   SCIP_CALL( SCIPallocBlockMemory(scip, info) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*info)->lbs, size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*info)->ubs, size) );
   (*info)->size = size;

   /* check for bracket */
   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, '(') )
   {
      syntaxError(scip, fzninput, "expecting  <(> after <output_array>");
      return SCIP_OKAY;
   }

   while( getNextToken(scip, fzninput) && !isChar(fzninput->token, ']') )
   {
      parseRange(scip, fzninput, &type, &lb, &ub);

      if( fzninput->haserror )
         return SCIP_OKAY;

      assert(type == FZN_INT);

      if( nelements == size )
      {
         size *= 2;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(*info)->lbs, (*info)->size, size) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(*info)->ubs, (*info)->size, size) );
         (*info)->size = size;
      }

      /* we assume integer bounds */
      (*info)->lbs[nelements] = (int) lb;
      (*info)->ubs[nelements] = (int) ub;
      nelements++;
   }

   (*info)->ndims = nelements;

   /* check for colon */
   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, ')') )
      syntaxError(scip, fzninput, "expecting  <)>");

   return SCIP_OKAY;
}

/** parse identifier name without annotations */
static
SCIP_RETCODE parseName(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   char*                 name,               /**< pointer to store the name */
   SCIP_Bool*            output,             /**< pointer to store if the name has the annotations to output */
   DIMENSIONS**          info                /**< pointer to store the output dimension information if one */
   )
{
   if( output != NULL )
      (*output) = FALSE;

   /* check for colon */
   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, ':') )
   {
      syntaxError(scip, fzninput, "expecting colon <:>");
      return SCIP_OKAY;
   }

   /* parse identifier name */
   if( !getNextToken(scip, fzninput) || !isIdentifier(fzninput->token) )
   {
      syntaxError(scip, fzninput, "expecting identifier name");
      return SCIP_OKAY;
   }

   /* copy identifier name */
   (void)SCIPsnprintf(name, FZN_BUFFERLEN-1, (const char*)fzninput->token);

   /* search for an assignment; therefore, skip annotations */
   do
   {
      if( !getNextToken(scip, fzninput) )
      {
         syntaxError(scip, fzninput, "expected at least a semicolon to close statement");
         return SCIP_OKAY;
      }

      /* check if the name has the annotation to be part of the output */
      if( equalTokens(scip, fzninput->token, "output_var") && output != NULL )
         (*output) = TRUE;
      else if( equalTokens(scip, fzninput->token, "output_array") && output != NULL)
      {
         (*output) = TRUE;
         assert(info != NULL);
         SCIP_CALL( parseOutputDimensioninfo(scip, fzninput, info) );
      }

      if( isEndStatement(fzninput) )
         break;
   }
   while( !isChar(fzninput->token, '=') );

   /* push back '=' or ';' */
   pushToken(fzninput);

   return SCIP_OKAY;
}

/** parse variable/constant (array) type (integer, float, bool, or set) */
static
void parseType(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   FZNNUMBERTYPE*        type,               /**< pointer to store the number type */
   SCIP_Real*            lb,                 /**< pointer to store the lower bound */
   SCIP_Real*            ub                  /**< pointer to store the lower bound */
   )
{
   if( !getNextToken(scip, fzninput) || isEndStatement(fzninput) )
   {
      syntaxError(scip, fzninput, "missing token");
      return;
   }

   *lb = -SCIPinfinity(scip);
   *ub = SCIPinfinity(scip);

   /* parse variable type or bounds */
   if( equalTokens(scip, fzninput->token, "bool") )
   {
      *type = FZN_BOOL;
      *lb = 0.0;
      *ub = 1.0;
   }
   else if( equalTokens(scip, fzninput->token, "float") )
      *type = FZN_FLOAT;
   else if( equalTokens(scip, fzninput->token, "int") )
      *type = FZN_INT;
   else if( equalTokens(scip, fzninput->token, "set") || isChar(fzninput->token, '{') )
   {
      SCIPwarningMessage(scip, "sets are not supported yet\n");
      fzninput->valid = FALSE;
      return;
   }
   else
   {
      /* the type is not explicitly given; it is given through the a range
       * expression; therefore, push back the current token since it
       * belongs to the range expression */
      pushToken(fzninput);
      parseRange(scip, fzninput, type, lb, ub);

      if( fzninput->haserror )
         return;
   }

   SCIPdebugMsg(scip, "range =  [%g,%g]\n", *lb, *ub);

   assert(*lb <= *ub);
}

/** applies assignment */
static
SCIP_RETCODE applyVariableAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   SCIP_VAR*             var,                /**< variable to assign something */
   FZNNUMBERTYPE         type,               /**< number type */
   const char*           assignment          /**< assignment */
   )
{
   FZNCONSTANT* constant;
   SCIP_VAR* linkVar;
   SCIP_Bool boolvalue;
   SCIP_Real realvalue;
   SCIP_Real fixvalue;
   SCIP_Real vals[] = {1.0,-1.0};

   linkVar = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) assignment);
   constant = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, (char*) assignment);

   realvalue = SCIP_INVALID;
   boolvalue = FALSE;

   if( linkVar == NULL )
   {
      if( isBoolExp(assignment, &boolvalue) && type == FZN_BOOL )
         fixvalue = (SCIP_Real) boolvalue;
      else if( isValue(assignment, &realvalue) && type != FZN_BOOL )
         fixvalue = realvalue;
      else if( constant != NULL )
         fixvalue = constant->value;
      else
      {
         syntaxError(scip, fzninput, "assignment is not recognizable");
         return SCIP_OKAY;
      }

      /* create fixing constraint */
      SCIP_CALL( createLinearCons(scip, "fixing", 1, &var, vals, fixvalue, fixvalue, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );
   }
   else
   {
      SCIP_VAR** vars;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
      vars[0] = var;
      vars[1] = linkVar;

      SCIP_CALL( createLinearCons(scip, "link", 2, vars, vals, 0.0, 0.0, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );

      SCIPfreeBufferArray(scip, &vars);
   }

   return SCIP_OKAY;
}

/** applies constant assignment expression */
static
SCIP_RETCODE createConstantAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNCONSTANT**         constant,           /**< pointer to constant */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< constant name */
   FZNNUMBERTYPE         type,               /**< number type */
   const char*           assignment          /**< assignment to apply */
   )
{
   SCIP_Bool boolvalue;
   SCIP_Real realvalue;
   SCIP_Real value;

   (*constant) = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, (char*) assignment);
   realvalue = SCIP_INVALID;
   boolvalue = FALSE;

   if( *constant != NULL )
   {
      /* check if the constant type fits */
      if( type != (*constant)->type )
      {
         syntaxError(scip, fzninput, "type error");
         return SCIP_OKAY;
      }

      value = (*constant)->value;
   }
   else if( isBoolExp(assignment, &boolvalue) && type == FZN_BOOL )
   {
      value = (SCIP_Real) boolvalue;
   }
   else if( isValue(assignment, &realvalue) && type != FZN_BOOL )
   {
      value = realvalue;
   }
   else
   {
      syntaxError(scip, fzninput, "assignment is not recognizable");
      return SCIP_OKAY;
   }

   /* get buffer memory for FZNCONSTANT struct */
   SCIP_CALL( SCIPallocBuffer(scip, constant) );

   (*constant)->type = type;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &(*constant)->name, name, (int) strlen(name) + 1) );
   (*constant)->value = value;

   /* store constant */
   if( fzninput->sconstants == fzninput->nconstants )
   {
      assert(fzninput->sconstants > 0);
      fzninput->sconstants *= 2;
      SCIP_CALL( SCIPreallocBufferArray(scip, &fzninput->constants, fzninput->sconstants) );
   }

   assert(fzninput->sconstants > fzninput->nconstants);
   fzninput->constants[fzninput->nconstants] = *constant;
   fzninput->nconstants++;

   SCIP_CALL( SCIPhashtableInsert(fzninput->constantHashtable, (void*) (*constant)) );

   return SCIP_OKAY;
}

/** parse array type ( (i) variable or constant; (ii) integer, float, bool, or set) */
static
void parseArrayType(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   SCIP_Bool*            isvararray,         /**< pointer to store if it is a variable or constant array */
   FZNNUMBERTYPE*        type,               /**< pointer to store number type */
   SCIP_Real*            lb,                 /**< pointer to store the lower bound */
   SCIP_Real*            ub                  /**< pointer to store the lower bound */
   )
{
   if( !getNextToken(scip, fzninput) || !equalTokens(scip, fzninput->token, "of") )
   {
      syntaxError(scip, fzninput, "expected keyword  <of>");
      return;
   }

   if( !getNextToken(scip, fzninput) )
   {
      syntaxError(scip, fzninput, "expected more tokens");
      return;
   }

   /* check if it is a variable or constant array */
   if( equalTokens(scip, fzninput->token, "var") )
      *isvararray = TRUE;
   else
   {
      /* push token back since it belongs to the type declaration */
      pushToken(fzninput);
      *isvararray = FALSE;
   }

   /* pares array type and range */
   parseType(scip, fzninput, type, lb, ub);
}

/** parse an array assignment */
static
SCIP_RETCODE parseArrayAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   char***               elements,           /**< pointer to string array to store the parsed elements */
   int*                  nelements,          /**< pointer to store the number of parsed elements */
   int                   selements           /**< size of the string array elements */
   )
{
   assert(scip != NULL);
   assert(fzninput != NULL);
   assert(*nelements >= 0);
   assert(selements >= *nelements);

   /* check for opening brackets */
   if( !getNextToken(scip, fzninput) ||  !isChar(fzninput->token, '[') )
   {
      syntaxError(scip, fzninput, "expected token <[>");
      return SCIP_OKAY;
   }

   SCIP_CALL( parseList(scip, fzninput, elements, nelements, selements) );

   if( hasError(fzninput) )
      return SCIP_OKAY;

   /* check for closing brackets */
   if( !getNextToken(scip, fzninput) ||  !isChar(fzninput->token, ']') )
      syntaxError(scip, fzninput, "expected token <]>");

   return SCIP_OKAY;
}

/** parse array dimension */
static
void parseArrayDimension(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   int*                  nelements           /**< pointer to store the size of the array */
   )
{
   FZNNUMBERTYPE type;
   SCIP_Real left;
   SCIP_Real right;

   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, '[') )
   {
      syntaxError(scip, fzninput, "expected token <[> for array dimension");
      return;
   }

   /* get array dimension */
   parseRange(scip, fzninput, &type, &left, &right);

   if( fzninput->haserror )
      return;

   if( type != FZN_INT || left != 1.0  || right <= 0.0 )
   {
      syntaxError(scip, fzninput, "invalid array dimension format");
      return;
   }

   *nelements = (int) right;

   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, ']') )
   {
      syntaxError(scip, fzninput, "expected token <]> for array dimension");
      return;
   }
}

/** creates and adds a variable to SCIP and stores it for latter use  in fzninput structure */
static
SCIP_RETCODE createVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   SCIP_VAR**            var,                /**< pointer to hold the created variable, or NULL */
   const char*           name,               /**< name of the variable */
   SCIP_Real             lb,                 /**< lower bound of the variable */
   SCIP_Real             ub,                 /**< upper bound of the variable */
   FZNNUMBERTYPE         type                /**< number type */
   )
{
   SCIP_VAR* varcopy;
   SCIP_VARTYPE vartype;

   assert(scip != NULL);
   assert(fzninput != NULL);
   assert(lb <= ub);

   switch(type)
   {
   case FZN_BOOL:
      vartype = SCIP_VARTYPE_BINARY;
      break;
   case FZN_INT:
      vartype = SCIP_VARTYPE_INTEGER;
      break;
   case FZN_FLOAT:
      vartype = SCIP_VARTYPE_CONTINUOUS;
      break;
   default:
      syntaxError(scip, fzninput, "unknown variable type");
      return SCIP_OKAY;
   }

   /* create variable */
   SCIP_CALL( SCIPcreateVar(scip, &varcopy, name, lb, ub, 0.0, vartype, !(fzninput->dynamiccols), fzninput->dynamiccols, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, varcopy) );

   SCIPdebugMsg(scip, "created variable\n");
   SCIPdebug( SCIP_CALL( SCIPprintVar(scip, varcopy, NULL) ) );

   /* variable name should not exist before */
   assert(SCIPhashtableRetrieve(fzninput->varHashtable, varcopy) == NULL);

   /* insert variable into the hashmap for later use in the constraint section */
   SCIP_CALL( SCIPhashtableInsert(fzninput->varHashtable, varcopy) );

   /* copy variable pointer before releasing the variable to keep the pointer to the variable */
   if( var != NULL )
      *var = varcopy;

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &varcopy) );

   return SCIP_OKAY;
}


/** parse variable array assignment and create the variables */
static
SCIP_RETCODE parseVariableArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< array name */
   int                   nvars,              /**< number of variables */
   FZNNUMBERTYPE         type,               /**< number type */
   SCIP_Real             lb,                 /**< lower bound of the variables */
   SCIP_Real             ub,                 /**< lower bound of the variables */
   DIMENSIONS*           info                /**< dimension information */
   )
{
   SCIP_VAR** vars;
   char varname[FZN_BUFFERLEN];
   int v;

   /* create variables and add them to the problem */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   for( v = 0; v < nvars; ++v )
   {
      (void) SCIPsnprintf(varname, FZN_BUFFERLEN, "%s[%d]", name, v + 1);

      /* create variable */
      SCIP_CALL( createVariable(scip, fzninput, &vars[v], varname, lb, ub, type) );
   }

   if( !getNextToken(scip, fzninput) )
   {
      syntaxError(scip, fzninput, "expected semicolon");
   }
   else
   {
      if( isChar(fzninput->token, '=') )
      {
         char** assigns;
         int nassigns;

         SCIP_CALL( SCIPallocBufferArray(scip, &assigns, nvars) );
         nassigns = 0;

         SCIP_CALL( parseArrayAssignment(scip, fzninput, &assigns, &nassigns, nvars) );

         if(!hasError(fzninput) )
         {
            for( v = 0; v < nvars && !hasError(fzninput); ++v )
            {
               /* parse and apply assignment */
               SCIP_CALL( applyVariableAssignment(scip, fzninput, vars[v], type, assigns[v]) );
            }
         }

         freeStringBufferArray(scip, assigns, nassigns);
      }
      else
      {
         /* push back the ';' */
         assert( isEndStatement(fzninput) );
         pushToken(fzninput);
      }

      if( info != NULL )
      {
         SCIP_CALL( readerdataAddOutputvararray(scip, readerdata, name, vars, nvars, type, info) );
      }

      /* add variable information to fzninput since this array name might be used later in the fzn file */
      SCIP_CALL( fzninputAddVararray(scip, fzninput, name, vars, nvars, type, info) );
   }

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** parse constant array assignment and create the constants */
static
SCIP_RETCODE parseConstantArray(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< array name */
   int                   nconstants,         /**< number of constants */
   FZNNUMBERTYPE         type                /**< number type */
   )
{
   FZNCONSTANT** constants;
   char** assigns;
   char constantname[FZN_BUFFERLEN];
   int nassigns;
   int c;

   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, '=') )
   {
      syntaxError(scip, fzninput, "expected token <=>");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &assigns, nconstants) );
   SCIP_CALL( SCIPallocBufferArray(scip, &constants, nconstants) );
   nassigns = 0;

   SCIP_CALL( parseArrayAssignment(scip, fzninput, &assigns, &nassigns, nconstants) );

   if( !hasError(fzninput) )
   {
      for( c = 0; c < nconstants; ++c )
      {
         (void) SCIPsnprintf(constantname, FZN_BUFFERLEN, "%s[%d]", name, c + 1);
         SCIP_CALL( createConstantAssignment(scip, &constants[c], fzninput, constantname, type, assigns[c]) );
      }

      /* add variable information to fzninput since this array name might be used later in the fzn file */
      SCIP_CALL( fzninputAddConstarray(scip, fzninput, name, constants, nconstants, type) );
   }

   SCIPfreeBufferArray(scip, &constants);
   freeStringBufferArray(scip, assigns, nassigns);

   return SCIP_OKAY;
}

/** parse predicate expression */
static
SCIP_RETCODE parsePredicate(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   /* mark predicate expression as comment such that it gets skipped */
   fzninput->comment = TRUE;

   return SCIP_OKAY;
}

/** parse array expression */
static
SCIP_RETCODE parseArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   FZNNUMBERTYPE type;
   DIMENSIONS* info;
   int nelements;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool isvararray;
   SCIP_Bool output;
   char name[FZN_BUFFERLEN];

   assert(scip != NULL);
   assert(fzninput != NULL);

   info = NULL;
   isvararray = FALSE;
   nelements = -1;

   SCIPdebugMsg(scip, "parse array expression\n");

   /* parse array dimension */
   parseArrayDimension(scip, fzninput, &nelements);
   assert(hasError(fzninput) || nelements > 0);

   if( hasError(fzninput) )
      return SCIP_OKAY;

   /* parse array type ( (i) variable or constant; (ii) integer, float, bool, or set) */
   parseArrayType(scip, fzninput, &isvararray, &type, &lb, &ub);

   if( hasError(fzninput) )
      return SCIP_OKAY;

   /* parse array name */
   SCIP_CALL( parseName(scip, fzninput, name, &output, &info) );
   assert(!output || info != NULL);

   if( hasError(fzninput) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "found <%s> array named <%s> of type <%s> and size <%d> with bounds [%g,%g] (output %u)\n",
      isvararray ? "variable" : "constant", name,
      type == FZN_BOOL ? "bool" : type == FZN_INT ? "integer" : "float", nelements, lb, ub, output);

   if( isvararray )
      SCIP_CALL( parseVariableArray(scip, readerdata, fzninput, name, nelements, type, lb, ub, info) );
   else
      SCIP_CALL( parseConstantArray(scip, fzninput, name, nelements, type) );

   freeDimensions(scip, &info);

   return SCIP_OKAY;
}

/** parse variable expression */
static
SCIP_RETCODE parseVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   SCIP_VAR* var;
   FZNNUMBERTYPE type;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool output;
   char assignment[FZN_BUFFERLEN];
   char name[FZN_BUFFERLEN];

   assert(scip != NULL);
   assert(fzninput != NULL);

   SCIPdebugMsg(scip, "parse variable expression\n");

   /* pares variable type and range */
   parseType(scip, fzninput, &type, &lb, &ub);

   if( hasError(fzninput) )
      return SCIP_OKAY;

   /* parse variable name without annotations */
   SCIP_CALL( parseName(scip, fzninput, name, &output, NULL) );

   if( hasError(fzninput) )
      return SCIP_OKAY;

   assert(type == FZN_BOOL || type == FZN_INT || type == FZN_FLOAT);

   /* create variable */
   SCIP_CALL( createVariable(scip, fzninput, &var, name, lb, ub, type) );

   /* check if the variable should be part of the output */
   if( output )
   {
      SCIP_CALL( readerdataAddOutputvar(scip, readerdata, var, type) );
   }

   if( !getNextToken(scip, fzninput) )
   {
      syntaxError(scip, fzninput, "expected semicolon");
      return SCIP_OKAY;
   }

   if( isChar(fzninput->token, '=') )
   {
      /* parse and flatten assignment */
      flattenAssignment(scip, fzninput, assignment);

      /* apply assignment */
      SCIP_CALL( applyVariableAssignment(scip, fzninput, var, type, assignment) );
   }
   else
      pushToken(fzninput);

   return SCIP_OKAY;
}

/** parse constant expression */
static
SCIP_RETCODE parseConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   FZNNUMBERTYPE         type                /**< constant type */
   )
{
   FZNCONSTANT* constant;
   char name[FZN_BUFFERLEN];
   char assignment[FZN_BUFFERLEN];

   assert(scip != NULL);
   assert(fzninput != NULL);
   assert(type == FZN_INT || type == FZN_FLOAT || type == FZN_BOOL);

   SCIPdebugMsg(scip, "parse constant expression\n");

   /* parse name of the constant */
   SCIP_CALL( parseName(scip, fzninput, name, NULL, NULL) );

   if( hasError(fzninput) )
      return SCIP_OKAY;

   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, '=') )
   {
      syntaxError(scip, fzninput, "expected token <=>");
      return SCIP_OKAY;
   }

   /* the assignment has to be an other constant or a suitable value */
   flattenAssignment(scip, fzninput, assignment);

   /* applies constant assignment and creates constant */
   SCIP_CALL( createConstantAssignment(scip, &constant, fzninput, name, type, assignment) );

   return SCIP_OKAY;
}

/** evaluates current token as constant */
static
void parseValue(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   SCIP_Real*            value,              /**< pointer to store value */
   const char*           assignment          /**< assignment to parse a value */
   )
{
   if( isValue(assignment, value) )
      return;

   /* if it is an identifier name, it has to belong to a constant or fixed variable */
   if( isIdentifier(assignment) )
   {
      FZNCONSTANT* constant;

      /* identifier has to be one of a constant */
      constant = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, (char*) assignment);

      if( constant == NULL )
      {
         SCIP_VAR* var;

         var = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) assignment);

         if( var == NULL )
            syntaxError(scip, fzninput, "unknown constant name");
         else
         {
            if( SCIPisEQ(scip, SCIPvarGetLbOriginal(var), SCIPvarGetUbOriginal(var)) )
               (*value) = SCIPvarGetLbOriginal(var);
            else
               syntaxError(scip, fzninput, "expected fixed variable");
         }
      }
      else
         (*value) = constant->value;
   }
   else
      syntaxError(scip, fzninput, "expected constant expression");
}

/** parse array expression containing constants */
static
SCIP_RETCODE parseConstantArrayAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   SCIP_Real**           vals,               /**< pointer to value array */
   int*                  nvals,              /**< pointer to store the number if values */
   int                   sizevals            /**< size of the vals array */
   )
{
   int c;

   assert(*nvals <= sizevals);

   /* check for next token */
   if( !getNextToken(scip, fzninput) )
   {
      syntaxError(scip, fzninput, "expected constant array");
      return SCIP_OKAY;
   }

   /* check if an array is given explicitly */
   if( isChar(fzninput->token, '[') )
   {
      char** elements;
      SCIP_Real value;
      int nelements;

      SCIP_CALL( SCIPallocBufferArray(scip, &elements, sizevals) );
      nelements = 0;

      value = 0.0;

      /* push back '[' which closes the list */
      pushToken(fzninput);

      /* pares array assignment */
      SCIP_CALL( parseArrayAssignment(scip, fzninput, &elements, &nelements, sizevals) );

      if( sizevals <= *nvals + nelements )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vals, *nvals + nelements) );
      }

      for( c = 0; c < nelements && !hasError(fzninput); ++c )
      {
         parseValue(scip, fzninput, &value, elements[c]);
         assert(!hasError(fzninput));

         (*vals)[(*nvals)] = value;
         (*nvals)++;
      }

      freeStringBufferArray(scip, elements, nelements);
   }
   else
   {
      /* array is not given explicitly; therefore, check constant array data base if the given constant array name was
       * parsed before
       */

      CONSTARRAY* constarray;

      constarray = findConstarray(scip, fzninput, fzninput->token);

      if( constarray != NULL )
      {
         /* ensure variable array size */
         if( sizevals <= *nvals + constarray->nconstants )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, vals, *nvals + constarray->nconstants) );
         }

         for( c = 0; c < constarray->nconstants; ++c )
         {
            (*vals)[(*nvals)] = constarray->constants[c]->value;
            (*nvals)++;
         }
      }
      else
      {
         /* there is no constant array with the given name; therefore check the variable array data base if such an
          * array exist with fixed variables
          */

         VARARRAY* vararray;

         vararray = findVararray(scip, fzninput, fzninput->token);

         if( vararray == NULL )
         {
            syntaxError(scip, fzninput, "unknown constants array name");
         }
         else
         {
            /* ensure variable array size */
            if( sizevals <= *nvals + vararray->nvars )
            {
               SCIP_CALL( SCIPreallocBufferArray(scip, vals, *nvals + vararray->nvars) );
            }

            for( c = 0; c < vararray->nvars; ++c )
            {
               SCIP_VAR* var;

               var = vararray->vars[c];
               assert(var != NULL);

               if( SCIPisEQ(scip, SCIPvarGetLbOriginal(var), SCIPvarGetUbOriginal(var)) )
               {
                  (*vals)[(*nvals)] = SCIPvarGetLbOriginal(var);
                  (*nvals)++;
               }
               else
               {
                  syntaxError(scip, fzninput, "variable array contains unfixed variable");
                  break;
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** parse array expression containing variables */
static
SCIP_RETCODE parseVariableArrayAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   SCIP_VAR***           vars,               /**< pointer to variable array */
   int*                  nvars,              /**< pointer to store the number if variable */
   int                   sizevars            /**< size of the variable array */
   )
{
   int v;

   assert(*nvars <= sizevars);

   /* check for next token */
   if( !getNextToken(scip, fzninput) )
   {
      syntaxError(scip, fzninput, "expected constant array");
      return SCIP_OKAY;
   }

   if( isChar(fzninput->token, '[') )
   {
      char** elements;
      int nelements;

      SCIP_CALL( SCIPallocBufferArray(scip, &elements, sizevars) );
      nelements = 0;

      /* push back '[' which closes the list */
      pushToken(fzninput);

      SCIP_CALL( parseArrayAssignment(scip, fzninput, &elements, &nelements, sizevars) );

      if( sizevars <= *nvars + nelements )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, vars, *nvars + nelements) );
      }

      for( v = 0; v < nelements; ++v )
      {
         (*vars)[(*nvars)] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, elements[v]);

         if( (*vars)[(*nvars)] == NULL )
         {
            /* since the given element does not correspond to a variable name
             * it might be the case that it is a constant which can be seen as
             * as a fixed variable
             */

            FZNCONSTANT* constant;
            SCIP_Real value;

            constant = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, (char*) elements[v]);

            if( constant != NULL )
            {
               assert(constant->type == FZN_FLOAT);
               value = constant->value;
            }
            else if(!isValue(elements[v], &value) )
            {
               char* tmptoken;

               tmptoken = fzninput->token;
               fzninput->token = elements[v];
               syntaxError(scip, fzninput, "expected variable name or constant");

               fzninput->token = tmptoken;
               break;
            }

            /* create a fixed variable */
            SCIP_CALL( createVariable(scip, fzninput, &(*vars)[*nvars], elements[v], value, value, FZN_FLOAT) );
         }

         (*nvars)++;
      }

      freeStringBufferArray(scip, elements, nelements);
   }
   else
   {
      VARARRAY* vararray;

      vararray = findVararray(scip, fzninput, fzninput->token);

      if( vararray != NULL )
      {
         assert(vararray != NULL);

         /* ensure variable array size */
         if( sizevars <= *nvars + vararray->nvars )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, vars, *nvars + vararray->nvars) );
         }

         for( v = 0; v < vararray->nvars; ++v )
         {
            (*vars)[(*nvars)] = vararray->vars[v];
            (*nvars)++;
         }
      }
      else
         syntaxError(scip, fzninput, "unknown variable array name");
   }

   return SCIP_OKAY;
}

/** parse linking statement */
static
SCIP_RETCODE parseQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name                /**< name of constraint */
   )
{
   char** elements;
   int nelements;

   SCIP_CALL( SCIPallocBufferArray(scip, &elements, 3) );
   nelements = 0;

   /* parse the list of three elements */
   SCIP_CALL( parseList(scip, fzninput, &elements, &nelements, 3) );
   assert(nelements == 3);

   if( !hasError(fzninput) )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real rhs;
      int v;

      rhs = 0.0;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, 3) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, 3) );

      for( v = 0; v < 3; ++v )
      {
         /* collect variable if constraint identifier is a variable */
         vars[v] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) elements[v]);

         /* parse the numeric value otherwise */
         if( vars[v] == NULL )
         {
            parseValue(scip, fzninput, &vals[v], elements[v]);
            assert(!hasError(fzninput));
         }
         else
            vals[v] = SCIP_INVALID;
      }

      /* the first two identifiers are proper variables => the constraints is indeed quadratic */
      if( vars[0] != NULL && vars[1] != NULL )
      {
         SCIP_Real quadval;
         quadval = 1.0;

         /* we might have an additional linear term or just a constant */
         if( vars[2] != NULL )
         {
            SCIP_Real linval;
            linval = -1.0;

            SCIP_CALL( createQuadraticCons(scip, name, 1, &vars[2], &linval, 1, &vars[0], &vars[1], &quadval, rhs, rhs,
                  fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );
         }
         else
         {
            rhs += vals[2];
            SCIP_CALL( createQuadraticCons(scip, name, 0, NULL, NULL, 1, &vars[0], &vars[1], &quadval, rhs, rhs,
                  fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows));
         }
      }
      else if( vars[0] != NULL || vars[1] != NULL )
      {
         int nvars;
         nvars = 1;

         /* the left hand side of x*y = z is linear (but not constant) */
         if( vars[0] == NULL )
            SCIPswapPointers((void**)&vars[0], (void**)&vars[1]);
         else
            SCIPswapPointers((void**)&vals[0], (void**)&vals[1]);

         /* after swapping, the variable and the coefficient should stand in front */
         assert(vars[0] != NULL && vals[0] != SCIP_INVALID ); /*lint !e777*/

         /* the right hand side might be a variable or a constant */
         if( vars[2] != NULL )
         {
            SCIPswapPointers((void**)&vars[1], (void**)&vars[2]);
            vals[1] = -1.0;
            nvars++;
         }
         else
         {
            assert(vals[2] != SCIP_INVALID); /*lint !e777*/
            rhs += vals[2];
         }

         SCIP_CALL( createLinearCons(scip, name, nvars, vars, vals, rhs, rhs, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );
      }
      else
      {
         /* the left hand side of x*y = z is constant */
         assert(vals[0] != SCIP_INVALID && vals[1] != SCIP_INVALID); /*lint !e777*/

         rhs = rhs - vals[0]*vals[1];

         /* the right hand side might be a variable or a constant */
         if( vars[2] != NULL )
         {
            SCIP_Real val;
            val = -1.0;
            SCIP_CALL( createLinearCons(scip, name, 1, &vars[2], &val, rhs, rhs, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );
         }
         else
         {
            assert(vals[2] != SCIP_INVALID); /*lint !e777*/
            rhs += vals[2];
            SCIP_CALL( createLinearCons(scip, name, 0, NULL, NULL, rhs, rhs, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );
         }
      }

      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* free elements array */
   freeStringBufferArray(scip, elements, nelements);

   return SCIP_OKAY;
}

/** parse aggregation statement (plus, minus, negate) */
static
SCIP_RETCODE parseAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< name of constraint */
   const char*           type                /**< linear constraint type */
   )
{
   /* here we take care of the three expression
    *
    * - int_plus(x1,x2,x3)   -> x1 + x2 == x3
    * - int_minus(x1,x2,x3)  -> x1 - x2 == x3
    * - int_negate(x1,x2)    -> x1 + x2 == 0
    */
   char** elements;
   int nelements;

   SCIP_CALL( SCIPallocBufferArray(scip, &elements, 3) );
   nelements = 0;

   /* parse the list of three elements */
   SCIP_CALL( parseList(scip, fzninput, &elements, &nelements, 3) );
   assert(nelements == 3 || nelements == 2);

   if( !hasError(fzninput) )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real value;
      SCIP_Real rhs;
      int nvars;

      nvars = 0;
      rhs = 0.0;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, 3) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, 3) );

      /* parse first element */
      vars[nvars] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) elements[0]);
      if( vars[nvars] == NULL )
      {
         parseValue(scip, fzninput, &value, elements[0]);
         assert(!hasError(fzninput));

         rhs -= value;
      }
      else
      {
         vals[nvars] = 1.0;
         nvars++;
      }

      /* parse second element */
      vars[nvars] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) elements[1]);
      if( vars[nvars] == NULL )
      {
         parseValue(scip, fzninput, &value, elements[1]);
         assert(!hasError(fzninput));

         if( equalTokens(scip, type, "minus") )
            rhs += value;
         else
            rhs -= value;
      }
      else
      {
         if( equalTokens(scip, type, "minus") )
         {
            /* in case of minus the second element get a -1.0 as coefficient */
            vals[nvars] = -1.0;
         }
         else
            vals[nvars] = 1.0;

         nvars++;
      }

      if( !equalTokens(scip, type, "negate") )
      {
         /* parse third element in case of "minus" or "plus" */
         vars[nvars] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) elements[2]);
         if( vars[nvars] == NULL )
         {
            parseValue(scip, fzninput, &value, elements[2]);
            assert(!hasError(fzninput));

            rhs += value;
         }
         else
         {
            vals[nvars] = -1.0;
            nvars++;
         }
      }

      SCIP_CALL( createLinearCons(scip, name, nvars, vars, vals, rhs, rhs, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );

      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* free elements array */
   freeStringBufferArray(scip, elements, nelements);
   return SCIP_OKAY;
}

/** parse linking statement */
static
SCIP_RETCODE parseLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           name,               /**< name of constraint */
   const char*           type,               /**< linear constraint type */
   SCIP_Real             sidevalue           /**< side value of constraint */
   )
{
   char** names;
   SCIP_Real lhs = SCIP_INVALID;
   SCIP_Real rhs = SCIP_INVALID;
   int nnames;

   nnames = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &names, 2) );

   SCIP_CALL( parseList(scip, fzninput, &names, &nnames, 2) );
   assert(nnames == 2);

   if( hasError(fzninput) )
      goto TERMINATE;

   /* compute left and right side */
   computeLinearConsSides(scip, fzninput, type, sidevalue, &lhs, &rhs);

   if( hasError(fzninput) )
      goto TERMINATE;

   SCIP_CALL( createLinking(scip, fzninput, name, names[0], names[1], lhs, rhs) );

 TERMINATE:
   freeStringBufferArray(scip, names, nnames);

   return SCIP_OKAY;
}

/** creates a linear constraint for an array operation */
static
CREATE_CONSTRAINT(createCoercionOpCons)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(fzninput != NULL);

   /* check if the function identifier name is array operation */
   if( !equalTokens(scip, fname, "int2float") && !equalTokens(scip, fname, "bool2int") )
      return SCIP_OKAY;

   SCIP_CALL( parseLinking(scip, fzninput, fname, "eq", 0.0) );

   *created = TRUE;

   return SCIP_OKAY;
}

/** creates a linear constraint for an array operation */
static
CREATE_CONSTRAINT(createSetOpCons)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(fzninput != NULL);

   /* check if the function identifier name is array operation */
   if( !equalTokens(scip, ftokens[0], "set") )
      return SCIP_OKAY;

   fzninput->valid = FALSE;
   SCIPwarningMessage(scip, "set operation are not supported yet\n");

   return SCIP_OKAY;
}

/** creates linear constraint for an array operation */
static
CREATE_CONSTRAINT(createArrayOpCons)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(fzninput != NULL);

   /* check if the function identifier name is array operation */
   if( !equalTokens(scip, ftokens[0], "array") )
      return SCIP_OKAY;

   fzninput->valid = FALSE;
   SCIPwarningMessage(scip, "array operation are not supported yet\n");

   return SCIP_OKAY;
}

/** creates a linear constraint for a logical operation */
static
CREATE_CONSTRAINT(createLogicalOpCons)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(fzninput != NULL);

   /* check if the function identifier name is array operation */
   if(nftokens < 2)
      return SCIP_OKAY;

   if(equalTokens(scip, ftokens[0], "bool") && nftokens == 2 )
   {
      char** elements;
      int nelements;

      /* the bool_eq constraint is processed in createComparisonOpCons() */
      if( equalTokens(scip, ftokens[1], "eq") || equalTokens(scip, ftokens[1], "ge") || equalTokens(scip, ftokens[1], "le")
         || equalTokens(scip, ftokens[1], "lt") || equalTokens(scip, ftokens[1], "gt") )
         return SCIP_OKAY;

      SCIP_CALL( SCIPallocBufferArray(scip, &elements, 3) );
      nelements = 0;

      SCIP_CALL( parseList(scip, fzninput, &elements, &nelements, 3) );

      if( !hasError(fzninput) )
      {
         SCIP_CONS* cons;
         SCIP_VAR** vars;
         int v;
         int nvars;

         if( equalTokens(scip, ftokens[1], "ne") || equalTokens(scip, ftokens[1], "not") )
            nvars = 2;
         else
            nvars = 3;

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

         /* collect variable if constraint identifier is a variable */
         for( v = 0; v < nvars; ++v )
         {
            vars[v] = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) elements[v]);

            if( vars[v] == NULL )
            {
               syntaxError(scip, fzninput, "unknown variable identifier name");
               goto TERMINATE;
            }
         }


         if( equalTokens(scip, ftokens[1], "ne" ) || equalTokens(scip, ftokens[1], "not") )
         {
            SCIP_Real vals[] = {1.0, 1.0};

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, fname, 2, vars, vals, 1.0, 1.0,
                  fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );

            *created = TRUE;
         }
         else if( equalTokens(scip, ftokens[1], "or" ) )
         {
            SCIP_CALL( SCIPcreateConsOr(scip, &cons, fname, vars[2], 2, vars,
                  fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );

            *created = TRUE;
         }
         else if( equalTokens(scip, ftokens[1], "and") )
         {
            SCIP_CALL( SCIPcreateConsAnd(scip, &cons, fname, vars[2], 2, vars,
                  fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );

            *created = TRUE;
         }
         else if( equalTokens(scip, ftokens[1], "xor") )
         {
            /* swap resultant to front */
            SCIPswapPointers((void**)&vars[0], (void**)&vars[2]);

            SCIP_CALL( SCIPcreateConsXor(scip, &cons, fname, FALSE, 3, vars,
                  fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );

            *created = TRUE;
         }
         else
         {
            fzninput->valid = FALSE;
            SCIPwarningMessage(scip, "logical operation <%s> is not supported yet\n", fname);
            goto TERMINATE;
         }

         SCIPdebugPrintCons(scip, cons, NULL);

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      TERMINATE:
         SCIPfreeBufferArray(scip, &vars);
      }

      /* free elements array */
      freeStringBufferArray(scip, elements, nelements);

   }
   else if(equalTokens(scip, ftokens[1], "bool") && nftokens == 3 )
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_VAR* resvar;
      int nvars;
      char** elements;
      int nelements;
      int size;

      if( !equalTokens(scip, ftokens[2], "or" ) && !equalTokens(scip, ftokens[2], "and" ) )
      {
         fzninput->valid = FALSE;
         SCIPwarningMessage(scip, "logical operation <%s> is not supported yet\n", fname);
         return SCIP_OKAY;
      }

      size = 10;
      nvars = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &elements, 1) );
      nelements = 0;

      SCIPdebugMsg(scip, "found and constraint <%s>\n", fname);

      /* parse operand variable array */
      SCIP_CALL( parseVariableArrayAssignment(scip, fzninput, &vars, &nvars, size) );

      /* check error and for the comma between the variable array and side value */
      if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
      {
         if( hasError(fzninput) )
            syntaxError(scip, fzninput, "unexpected error in fzn input");
         else
            syntaxError(scip, fzninput, "expected token <,>");

         goto TERMINATE2;
      }

      /* parse resultant variable array */
      SCIP_CALL( parseList(scip, fzninput, &elements, &nelements, 1) );
      resvar = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) elements[0]);

      /* check error and for the comma between the variable array and side value */
      if( hasError(fzninput) || resvar == NULL )
      {
         if( hasError(fzninput) )
            syntaxError(scip, fzninput, "unexpected error in fzn input");
         else
            syntaxError(scip, fzninput, "unknown variable identifier name");
         goto TERMINATE2;
      }

      /* create the constraint */
      if( equalTokens(scip, ftokens[2], "or" ) )
      {
         SCIP_CALL( SCIPcreateConsOr(scip, &cons, fname, resvar, nvars, vars,
               fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );
      }
      else
      {
         assert( equalTokens(scip, ftokens[2], "and") );

         SCIP_CALL( SCIPcreateConsAnd(scip, &cons, fname, resvar, nvars, vars,
               fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );
      }

      SCIPdebugPrintCons(scip, cons, NULL);
      *created = TRUE;

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   TERMINATE2:
      /* free elements array */
      freeStringBufferArray(scip, elements, nelements);
      SCIPfreeBufferArray(scip, &vars);
   }
   else if( equalTokens(scip, ftokens[1], "bool") )
   {
      fzninput->valid = FALSE;
      SCIPwarningMessage(scip, "logical operation <%s> is not supported yet\n", fname);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** creates a linear constraint for a comparison operation */
static
CREATE_CONSTRAINT(createComparisonOpCons)
{  /*lint --e{715}*/
   char assignment[FZN_BUFFERLEN];

   assert(scip != NULL);
   assert(fzninput != NULL);

   /* check if the function name ends of "reif" (reified constraint) which SCIP does not support yet */
   if( equalTokens(scip, ftokens[nftokens - 1], "reif") )
   {
      SCIPwarningMessage(scip, "reified constraints are not supported\n");
      fzninput->valid = FALSE;
      return SCIP_OKAY;
   }

   /* the last token can be
    * 'eq' -- equal
    * 'ne' -- not equal
    * 'lt' -- less than
    * 'gt' -- greater than
    * 'le' -- less or equal than
    * 'ge' -- greater or equal than
    *         => these are comparison constraints
    * 'plus'   -- addition
    * 'minus'  -- subtraction
    * 'negate' -- negation
    *             => these are aggregation constraints
    * 'times' -- multiplication
    *            => this is a nonlinear constraint
    */
   if( strlen(ftokens[nftokens - 1]) != 2 && nftokens != 2 )
      return SCIP_OKAY;

   /* check if any sets are involved in the constraint */
   if( equalTokens(scip, ftokens[0], "set") )
   {
      SCIPwarningMessage(scip, "constraints using sets are not supported\n");
      fzninput->valid = FALSE;
      return SCIP_OKAY;
   }

   /* check if the constraint is a 'not equal' one */
   if( equalTokens(scip, ftokens[nftokens - 1], "ne") )
   {
      SCIPwarningMessage(scip, "constraints with 'not equal' relation are not supported\n");
      fzninput->valid = FALSE;
      return SCIP_OKAY;
   }

   /* check if the constraint contains float variable and coefficients and '<' or '>' relation */
   if( equalTokens(scip, ftokens[0], "float") &&
      (equalTokens(scip, ftokens[nftokens - 1], "lt") || equalTokens(scip, ftokens[nftokens - 1], "gt") ) )
   {
      SCIPwarningMessage(scip, "constraints with '<' or '>' relation and continuous variables are not supported\n");
      fzninput->valid = FALSE;
      return SCIP_OKAY;
   }

   if( equalTokens(scip, ftokens[1], "lin") )
   {
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real sidevalue;
      int nvars;
      int nvals;
      int size;

      assert(nftokens == 3);

      size = 10;
      nvars = 0;
      nvals = 0;
      sidevalue = SCIP_INVALID;

      SCIP_CALL( SCIPallocBufferArray(scip, &vars, size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, size) );

      SCIPdebugMsg(scip, "found linear constraint <%s>\n", fname);

      /* pares coefficients array */
      SCIP_CALL( parseConstantArrayAssignment(scip, fzninput, &vals, &nvals, size) );

      /* check error and for the comma between the coefficient and variable array */
      if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
      {
         if( !hasError(fzninput) )
            syntaxError(scip, fzninput, "expected token <,>");

         goto TERMINATE;
      }

      /* pares variable array */
      SCIP_CALL( parseVariableArrayAssignment(scip, fzninput, &vars, &nvars, size) );

      /* check error and for the comma between the variable array and side value */
      if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
      {
         if( !hasError(fzninput) )
            syntaxError(scip, fzninput, "expected token <,>");

         goto TERMINATE;
      }

      /* pares sidevalue */
      flattenAssignment(scip, fzninput, assignment);
      parseValue(scip, fzninput, &sidevalue, assignment);

      if( !hasError(fzninput) )
      {
         SCIP_Real lhs = -SCIPinfinity(scip);
         SCIP_Real rhs = SCIPinfinity(scip);

         assert(sidevalue != SCIP_INVALID); /*lint !e777*/

         /* compute left and right side */
         computeLinearConsSides(scip, fzninput, ftokens[2], sidevalue, &lhs, &rhs);

         if( hasError(fzninput) )
            goto TERMINATE;

         SCIP_CALL( createLinearCons(scip, fname, nvars, vars, vals, lhs, rhs, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );
      }

   TERMINATE:
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals);
   }
   else if( equalTokens(scip, ftokens[1], "minus") || equalTokens(scip, ftokens[1], "plus") || equalTokens(scip, ftokens[1], "negate") )
   {
      assert(nftokens == 2);
      SCIP_CALL( parseAggregation(scip, fzninput, fname, ftokens[1]) );
   }
   else if( equalTokens(scip, ftokens[1], "eq") || equalTokens(scip, ftokens[1], "le") || equalTokens(scip, ftokens[1], "ge")
      || equalTokens(scip, ftokens[1], "lt") || equalTokens(scip, ftokens[1], "gt") )
   {
      assert(nftokens == 2);
      SCIP_CALL( parseLinking(scip, fzninput, fname, ftokens[1], 0.0) );
   }
   else if( equalTokens(scip, ftokens[1], "times") )
   {
      assert(nftokens == 2);
      SCIP_CALL( parseQuadratic(scip, fzninput, fname) );
   }
   else
   {
      syntaxError(scip, fzninput, "unknown constraint type");
   }

   *created = TRUE;

   return SCIP_OKAY;
}

/** creates an alldifferent constraint */
static
CREATE_CONSTRAINT(createAlldifferentOpCons)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
#ifdef ALLDIFFERENT
   SCIP_CONS* cons;
#endif
   int nvars;
   int size;

   assert(scip != NULL);
   assert(fzninput != NULL);

   /* check if the function identifier name is array operation */
   if( !equalTokens(scip, ftokens[0], "all") || !equalTokens(scip, ftokens[1], "different") )
      return SCIP_OKAY;

   size = 10;
   nvars = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, size) );

   SCIPdebugMsg(scip, "parse alldifferent expression\n");

   /* pares variable array */
   SCIP_CALL( parseVariableArrayAssignment(scip, fzninput, &vars, &nvars, size) );

#ifdef ALLDIFFERENT
   /* create alldifferent constraint */
   SCIP_CALL( SCIPcreateConsAlldifferent(scip, &cons, fname, nvars, vars,
         fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );

   SCIPdebugPrintCons(scip, cons, NULL);

   /* add and release the constraint to the problem */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   *created = TRUE;
#endif

   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** creates an alldifferent constraint */
static
CREATE_CONSTRAINT(createCumulativeOpCons)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* vals = NULL;
   int* durations = NULL;
   int* demands = NULL;
   SCIP_Real val;
   int capacity;
   char assignment[FZN_BUFFERLEN];

   int nvars;
   int ndurations;
   int ndemads;
   int size;
   int i;

   assert(scip != NULL);
   assert(fzninput != NULL);

   /* check if the function identifier name is array operation */
   if( !equalTokens(scip, ftokens[0], "cumulative") )
      return SCIP_OKAY;

   size = 10;
   nvars = 0;
   ndurations = 0;
   ndemads = 0;

   SCIPdebugMsg(scip, "parse cumulative expression\n");

   /* pares start time variable array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, size) );
   SCIP_CALL( parseVariableArrayAssignment(scip, fzninput, &vars, &nvars, size) );

   /* check error and for the comma between the variable array and side value */
   if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
   {
      if( !hasError(fzninput) )
         syntaxError(scip, fzninput, "expected token <,>");

      goto TERMINATE;
   }

   /* pares job duration array */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, size) );
   SCIP_CALL( parseConstantArrayAssignment(scip, fzninput, &vals, &ndurations, size) );

   SCIP_CALL( SCIPallocBufferArray(scip, &durations, ndurations) );
   for( i = 0; i < ndurations; ++i )
      durations[i] = (int)vals[i];

   /* check error and for the comma between the variable array and side value */
   if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
   {
      if( !hasError(fzninput) )
         syntaxError(scip, fzninput, "expected token <,>");

      goto TERMINATE;
   }

   /* pares job demand array */
   SCIP_CALL( parseConstantArrayAssignment(scip, fzninput, &vals, &ndemads, size) );

   SCIP_CALL( SCIPallocBufferArray(scip, &demands, ndemads) );
   for( i = 0; i < ndemads; ++i )
      demands[i] = (int)vals[i];

   /* check error and for the comma between the variable array and side value */
   if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
   {
      if( !hasError(fzninput) )
         syntaxError(scip, fzninput, "expected token <,>");

      goto TERMINATE;
   }

   /* parse cumulative capacity */
   flattenAssignment(scip, fzninput, assignment);
   parseValue(scip, fzninput, &val, assignment);
   assert(!hasError(fzninput));

   capacity = (int)val;

   assert(nvars == ndurations);
   assert(nvars == ndemads);

   /* create cumulative constraint */
   SCIP_CALL( SCIPcreateConsCumulative(scip, &cons, fname, nvars, vars, durations, demands, capacity,
         fzninput->initialconss, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, fzninput->dynamicconss, fzninput->dynamicrows, FALSE) );

   SCIPdebugPrintCons(scip, cons, NULL);

   /* add and release the constraint to the problem */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   assert(!hasError(fzninput));
   *created = TRUE;

 TERMINATE:
   /* free buffers */
   SCIPfreeBufferArrayNull(scip, &demands);
   SCIPfreeBufferArrayNull(scip, &durations);
   SCIPfreeBufferArrayNull(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/* function pointer array containing all function which can create a constraint */
static CREATE_CONSTRAINT((*constypes[])) =  {
   createCoercionOpCons,
   createSetOpCons,
   createLogicalOpCons,
   createArrayOpCons,
   createComparisonOpCons,
   createAlldifferentOpCons,
   createCumulativeOpCons
};

/** size of the function pointer array */
static int nconstypes = 7;


/** parse constraint expression */
static
SCIP_RETCODE parseConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   SCIP_VAR* var;
   char* tokens[4];
   char* token;
   char* nexttoken;
   char name[FZN_BUFFERLEN];
   char fname[FZN_BUFFERLEN];
   SCIP_Bool created;
   int ntokens;
   int i;
   int c;

   assert(scip != NULL);
   assert(fzninput != NULL);

   SCIPdebugMsg(scip, "parse constraint expression\n");

   /* get next token already flatten */
   flattenAssignment(scip, fzninput, name);

   /* check if constraint identifier is a variable */
   var = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) name);

   if( var != NULL )
   {
      SCIP_Real vals[] = {1.0};

      /* create fixing constraint */
      SCIP_CALL( createLinearCons(scip, "fixing", 1, &var, vals, 1.0, 1.0, fzninput->initialconss, fzninput->dynamicconss, fzninput->dynamicrows) );
      return SCIP_OKAY;
   }

   /* check constraint identifier name */
   if( !isIdentifier(name) )
   {
      syntaxError(scip, fzninput, "expected constraint identifier name");
      return SCIP_OKAY;
   }

   /* check if we have a opening parenthesis */
   if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, '(') )
   {
      syntaxError(scip, fzninput, "expected token <(>");
      return SCIP_OKAY;
   }

   /* copy function name */
   (void) SCIPsnprintf(fname, FZN_BUFFERLEN, "%s", name);

   /* truncate the function identifier name in separate tokens */
   token = SCIPstrtok(name, "_", &nexttoken);
   ntokens = 0;
   while( token != NULL )
   {
      if( ntokens == 4 )
         break;

      SCIP_CALL( SCIPduplicateBufferArray(scip, &(tokens[ntokens]), token, (int) strlen(token) + 1) ); /*lint !e866*/
      ntokens++;

      token = SCIPstrtok(NULL, "_", &nexttoken);
   }

   assert(token == NULL || tokens[0] != NULL); /*lint !e771*/
   for( i = 0; i < ntokens; ++i )
   {
      SCIPdebugMsgPrint(scip, "%s ", tokens[i]);
   }
   SCIPdebugMsgPrint(scip, "\n");

   created = FALSE;

   /* loop over all methods which can create a constraint */
   for( c = 0; c < nconstypes && !created && !hasError(fzninput); ++c )
   {
      SCIP_CALL( constypes[c](scip, fzninput, fname, tokens, ntokens, &created) );
   }

   /* check if a constraint was created */
   if( !hasError(fzninput) && !created )
   {
      fzninput->valid = FALSE;
      SCIPwarningMessage(scip, "constraint <%s> is not supported yet\n", fname);
   }

   /* free memory */
   for( i = ntokens - 1; i >= 0 ; --i )
   {
      SCIPfreeBufferArray(scip, &tokens[i]);
   }

   /* check for the closing parenthesis */
   if( !hasError(fzninput) && ( !getNextToken(scip, fzninput) || !isChar(fzninput->token, ')')) )
      syntaxError(scip, fzninput, "expected token <)>");

   return SCIP_OKAY;
}

/** parse solve item expression */
static
SCIP_RETCODE parseSolveItem(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNINPUT*             fzninput            /**< FZN reading data */
   )
{
   assert(scip != NULL);
   assert(fzninput != NULL);

   SCIPdebugMsg(scip, "parse solve item expression\n");

   if( !getNextToken(scip, fzninput) )
   {
      syntaxError(scip, fzninput, "expected solving specification");
      return SCIP_OKAY;
   }

   /* check for annotations */
   if( equalTokens(scip, fzninput->token, "::") )
   {
      /* skip the annotation */
      do
      {
         if( !getNextToken(scip, fzninput) )
            syntaxError(scip, fzninput, "expected more tokens");
      }
      while( !equalTokens(scip, fzninput->token, "satisfy")
         && !equalTokens(scip, fzninput->token, "minimize")
         && !equalTokens(scip, fzninput->token, "maximize") );
   }

   if( equalTokens(scip, fzninput->token, "satisfy") )
   {
      SCIPdebugMsg(scip, "detected a satisfiability problem\n");
   }
   else
   {
      SCIP_VAR* var;
      FZNCONSTANT* constant;
      char name[FZN_BUFFERLEN];

      if( equalTokens(scip, fzninput->token, "minimize") )
      {
         fzninput->objsense = SCIP_OBJSENSE_MINIMIZE;
         SCIPdebugMsg(scip, "detected a minimization problem\n");
      }
      else
      {
         assert(equalTokens(scip, fzninput->token, "maximize"));
         fzninput->objsense = SCIP_OBJSENSE_MAXIMIZE;
         SCIPdebugMsg(scip, "detected a maximization problem\n");
      }

      /* parse objective coefficients */

      /* parse and flatten assignment */
      flattenAssignment(scip, fzninput, name);

      var = (SCIP_VAR*) SCIPhashtableRetrieve(fzninput->varHashtable, (char*) name);
      constant = (FZNCONSTANT*) SCIPhashtableRetrieve(fzninput->constantHashtable, (char*) name);

      if( var != NULL )
      {
         SCIP_CALL(SCIPchgVarObj(scip, var, 1.0) );
      }
      else if( constant != NULL )
      {
         SCIPdebugMsg(scip, "optimizing a constant is equal to a satisfiability problem!\n");
      }
      else if( equalTokens(scip, name, "int_float_lin") )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int nvars;
         int nvals;
         int size;
         int v;

         nvars = 0;
         nvals = 0;
         size = 10;

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, size) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, size) );

         SCIPdebugMsg(scip, "found linear objective\n");

         if( !getNextToken(scip, fzninput) || !isChar(fzninput->token, '(') )
         {
            syntaxError(scip, fzninput, "expected token <(>");
            goto TERMINATE;
         }

         /* pares coefficients array for integer variables */
         SCIP_CALL( parseConstantArrayAssignment(scip, fzninput, &vals, &nvals, size) );

         /* check error and for the comma between the coefficient and variable array */
         if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
         {
            if( !hasError(fzninput) )
               syntaxError(scip, fzninput, "expected token <,>");

            goto TERMINATE;
         }

         /* pares coefficients array for continuous variables */
         SCIP_CALL( parseConstantArrayAssignment(scip, fzninput, &vals, &nvals, MAX(size, nvals)) );

         /* check error and for the comma between the coefficient and variable array */
         if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
         {
            if( !hasError(fzninput) )
               syntaxError(scip, fzninput, "expected token <,>");

            goto TERMINATE;
         }

         /* pares integer variable array */
         SCIP_CALL( parseVariableArrayAssignment(scip, fzninput, &vars, &nvars, size) );

         /* check error and for the comma between the variable array and side value */
         if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ',') )
         {
            if( !hasError(fzninput) )
               syntaxError(scip, fzninput, "expected token <,>");

            goto TERMINATE;
         }

         assert(nvars <= nvals);

         /* pares continuous variable array */
         SCIP_CALL( parseVariableArrayAssignment(scip, fzninput, &vars, &nvars, MAX(size, nvars)) );

         /* check error and for the ')' */
         if( hasError(fzninput) || !getNextToken(scip, fzninput) || !isChar(fzninput->token, ')') )
         {
            if( !hasError(fzninput) )
               syntaxError(scip, fzninput, "expected token <)>");

            goto TERMINATE;
         }

         assert( nvars == nvals );

         for( v = 0; v < nvars; ++v )
         {
            SCIP_CALL(SCIPchgVarObj(scip, vars[v], vals[v]) );
         }

      TERMINATE:
         SCIPfreeBufferArray(scip, &vars);
         SCIPfreeBufferArray(scip, &vals);
      }
      else
      {
         syntaxError(scip, fzninput, "unknown identifier expression for a objective function");
      }
   }

   return SCIP_OKAY;
}

/** reads a FlatZinc model */
static
SCIP_RETCODE readFZNFile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READERDATA*      readerdata,         /**< reader data */
   FZNINPUT*             fzninput,           /**< FZN reading data */
   const char*           filename            /**< name of the input file */
   )
{
   assert(scip != NULL);
   assert(readerdata != NULL);
   assert(fzninput != NULL);

   /* open file */
   fzninput->file = SCIPfopen(filename, "r");
   if( fzninput->file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* create two auxiliary variable for true and false values */
   SCIP_CALL( createVariable(scip, fzninput, NULL, "true", 1.0, 1.0, FZN_BOOL) );
   SCIP_CALL( createVariable(scip, fzninput, NULL, "false", 0.0, 0.0, FZN_BOOL) );

   /* parse through statements one-by-one */
   while( !SCIPfeof( fzninput->file ) && !hasError(fzninput) )
   {
      /* read the first token (keyword) of a new statement */
      if( getNextToken(scip, fzninput) )
      {
         if( equalTokens(scip, fzninput->token, "predicate") )
         {
            /* parse array expression containing constants or variables */
            SCIP_CALL( parsePredicate(scip, fzninput) );
         }
         else if( equalTokens(scip, fzninput->token, "array") )
         {
            /* parse array expression containing constants or variables */
            SCIP_CALL( parseArray(scip, readerdata, fzninput) );
         }
         else if( equalTokens(scip, fzninput->token, "constraint") )
         {
            /* parse a constraint */
            SCIP_CALL( parseConstraint(scip, fzninput) );
         }
         else if( equalTokens(scip, fzninput->token, "int") )
         {
            /* parse an integer constant */
            SCIP_CALL( parseConstant(scip, fzninput, FZN_INT) );
         }
         else if( equalTokens(scip, fzninput->token, "float") )
         {
            /* parse a float constant */
            SCIP_CALL( parseConstant(scip, fzninput, FZN_FLOAT) );
         }
         else if( equalTokens(scip, fzninput->token, "bool") )
         {
            /* parse a bool constant */
            SCIP_CALL( parseConstant(scip, fzninput, FZN_BOOL) );
         }
         else if( equalTokens(scip, fzninput->token, "set") )
         {
            /* deal with sets */
            SCIPwarningMessage(scip, "sets are not supported yet\n");
            fzninput->valid = FALSE;
            break;
         }
         else if( equalTokens(scip, fzninput->token, "solve") )
         {
            /* parse solve item (objective sense and objective function) */
            SCIP_CALL( parseSolveItem(scip, fzninput) );
         }
         else if( equalTokens(scip, fzninput->token, "var") )
         {
            /* parse variables */
            SCIP_CALL( parseVariable(scip, readerdata, fzninput) );
         }
         else if( equalTokens(scip, fzninput->token, "output") )
         {
            /* the output section is the last section in the flatzinc model and can be skipped */
            SCIPdebugMsg(scip, "skip output section\n");
            break;
         }
         else
         {
            FZNNUMBERTYPE type;
            SCIP_Real lb;
            SCIP_Real ub;

            /* check if the new statement starts with a range expression
             * which indicates a constant; therefore, push back the current token
             * since it belongs to the range expression */
            pushToken(fzninput);

            /* parse range to detect constant type */
            parseRange(scip, fzninput, &type, &lb, &ub);

            if( hasError(fzninput) )
               break;

            /* parse the remaining constant statement */
            SCIP_CALL( parseConstant(scip, fzninput, type) );

            if( hasError(fzninput) )
            {
               SCIPwarningMessage(scip, "unknown keyword <%s> skip statement\n", fzninput->token);
               SCIPABORT();
               return SCIP_OKAY;  /*lint !e527*/
            }
         }

         if( hasError(fzninput) )
            break;

         /* if the current statement got marked as comment continue with the next line */
         if( fzninput->comment )
            continue;

         /* each statement should be closed with a semicolon */
         if( !getNextToken(scip, fzninput) )
            syntaxError(scip, fzninput, "expected semicolon");

         /* check for annotations */
         if( equalTokens(scip, fzninput->token, "::") )
         {
            /* skip the annotation */
            do
            {
               if( !getNextToken(scip, fzninput) )
                  syntaxError(scip, fzninput, "expected more tokens");
            }
            while( !isEndStatement(fzninput) );
         }

         if( !isEndStatement(fzninput) )
            syntaxError(scip, fzninput, "expected semicolon");
      }
   }

   /* close file */
   SCIPfclose(fzninput->file);

   if( hasError(fzninput) )
   {
      SCIP_CALL( SCIPfreeProb(scip) );

      /* create empty problem */
      SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   }
   else
   {
      SCIP_CALL( SCIPsetObjsense(scip, fzninput->objsense) );
   }

   return SCIP_OKAY;
}


/*
 * Local methods (for writing)
 */


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
   int requiredsize;                         /* number of active variables */
   int v;

   assert( scip != NULL );
   assert( scalars != NULL );
   assert( nvars != NULL );
   assert( vars != NULL || *nvars == 0 );
   assert( constant != NULL );

   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

      /* avoid overflow by reallocation */
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
         assert(vars != NULL);
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &scalars[v], constant) );
      }

   return SCIP_OKAY;
}

/** ends the given line with '\\0' and prints it to the given file stream */
static
void writeBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   char*                 buffer,             /**< line */
   int                   bufferpos           /**< number of characters in buffer */
   )
{
   assert( scip != NULL );
   assert( buffer != NULL );

   if( bufferpos > 0 )
   {
      buffer[bufferpos] = '\0';

      SCIPinfoMessage(scip, file, "%s", buffer);
   }
}

/** appends extension to line and prints it to the give file stream if the line buffer get full */
static
SCIP_RETCODE appendBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   char**                buffer,             /**< buffer which should be extended */
   int*                  bufferlen,          /**< length of the buffer */
   int*                  bufferpos,          /**< current position in  the buffer */
   const char*           extension           /**< string to extend the line */
   )
{
   int newpos;

   assert( scip != NULL );
   assert( buffer != NULL );
   assert( bufferlen != NULL );
   assert( bufferpos != NULL );
   assert( extension != NULL );

   /* avoid overflow by reallocation */
   newpos = (*bufferpos) + (int)strlen(extension);
   if( newpos >= (*bufferlen) )
   {
      *bufferlen = MAX( newpos, 2*(*bufferlen) );

      SCIP_CALL( SCIPreallocBufferArray(scip, buffer, (*bufferlen)));
   }

   /* append extension to linebuffer */
   (void)strncpy((*buffer)+(*bufferpos), extension, (size_t)strlen(extension));
   *bufferpos = newpos;

   return SCIP_OKAY;
}

/* Writes a real value to a string with full precision, if fractional and adds a ".0" if integral */
static
void flattenFloat(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val,                /**< value to flatten */
   char*                 buffer              /**< string buffer to print in */
   )
{
   if( SCIPisIntegral(scip, val) )
      (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%.1f", SCIPround(scip, val));
   else
      (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%+.15g", val);

}

/* print row in FZN format to file stream */
static
SCIP_RETCODE printRow(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNOUTPUT*            fznoutput,          /**< output data structure containing the buffers to write to */
   const char*           type,               /**< row type ("eq", "le" or "ge") */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvars,              /**< number of variables */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             hasfloats           /**< are there continuous variables or coefficients in the constraint? */
   )
{
   SCIP_VAR* var;                            /* some variable */
   int v;                                    /* variable counter */
   char buffer[FZN_BUFFERLEN];
   char buffy[FZN_BUFFERLEN];

   assert( scip != NULL );
   assert( strcmp(type, "eq") == 0 || strcmp(type, "le") == 0 || strcmp(type, "ge") == 0 );

   /* Add a constraint of type float_lin or int_lin, depending on whether there are continuous variables or coefficients */
   SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos),"constraint ") );
   if( hasfloats )
      (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "float_lin_%s([",type);
   else
      (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "int_lin_%s([",type);
   SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos),buffer) );

   /* print all coefficients but the last one */
   for( v = 0; v < nvars-1; ++v )
   {
      if( hasfloats )
      {
         flattenFloat(scip, vals[v], buffy);
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%s, ", buffy);
      }
      else
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%.f, ", vals[v]);
      SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos), buffer) );
   }

   /* print last coefficient */
   if( nvars > 0 )
   {
      if( hasfloats )
      {
         flattenFloat(scip, vals[nvars-1], buffy);
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%s", buffy);
      }
      else
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%.f", vals[nvars-1]);

      SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos), buffer) );
   }

   SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos), "], [") );

   /* print all variables but the last one */
   for( v = 0; v < nvars-1; ++v )
   {
      var = vars[v];
      assert( var != NULL );

      if( hasfloats )
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%s%s, ", SCIPvarGetName(var), SCIPvarGetProbindex(var) < fznoutput->ndiscretevars ? "_float" : "");
      else
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%s, ", SCIPvarGetName(var) );
      SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos), buffer) );
   }

   /* print last variable */
   if( nvars > 0 )
   {
      if( hasfloats )
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%s%s",SCIPvarGetName(vars[nvars-1]),
            SCIPvarGetProbindex(vars[nvars-1]) < fznoutput->ndiscretevars ? "_float" : "");
      else
         (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%s",SCIPvarGetName(vars[nvars-1]));

      SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos),buffer) );
   }

   SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos), "], ") );

   /* print right hand side */
   if( SCIPisZero(scip, rhs) )
      rhs = 0.0;

   if( hasfloats )
   {
      flattenFloat(scip, rhs, buffy);
      (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%s);\n",buffy);
   }
   else
      (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "%.f);\n",rhs);
   SCIP_CALL( appendBuffer(scip, &(fznoutput->consbuffer), &(fznoutput->consbufferlen), &(fznoutput->consbufferpos),buffer) );

   return SCIP_OKAY;
}

/** prints given linear constraint information in FZN format to file stream */
static
SCIP_RETCODE printLinearCons(
   SCIP*                 scip,               /**< SCIP data structure */
   FZNOUTPUT*            fznoutput,          /**< output data structure containing the buffers to write to */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SCIP_Bool             mayhavefloats       /**< may there be continuous variables in the constraint? */
   )
{
   SCIP_VAR** activevars;                    /* active problem variables of a constraint */
   SCIP_Real* activevals;                    /* coefficients in the active representation */

   SCIP_Real activeconstant;                 /* offset (e.g., due to fixings) in the active representation */
   int nactivevars;                          /* number of active problem variables */
   int v;                                    /* variable counter */

   char buffer[FZN_BUFFERLEN];
   SCIP_Bool hasfloats;

   assert( scip != NULL );
   assert( vars != NULL || nvars == 0 );
   assert( fznoutput != NULL );
   assert( lhs <= rhs );

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   /* duplicate variable and value array */
   nactivevars = nvars;
   hasfloats = FALSE;
   activevars = NULL;
   activeconstant = 0.0;

   if( vars != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   }

   if( vals != NULL )
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, activevars, activevals, &nactivevars, &activeconstant, transformed) );

   /* If there may be continuous variables or coefficients in the constraint, scan for them */
   if( mayhavefloats )
   {
      /* fractional sides trigger a constraint to be of float type */
      if( !SCIPisInfinity(scip, -lhs) )
         hasfloats = hasfloats || !SCIPisIntegral(scip, lhs-activeconstant);
      if( !SCIPisInfinity(scip, rhs) )
         hasfloats = hasfloats || !SCIPisIntegral(scip, rhs-activeconstant);

      /* any continuous variable or fractional variable coefficient triggers a constraint to be of float type */
      for( v = 0; v < nactivevars && !hasfloats; v++ )
      {
         SCIP_VAR* var;

         assert(activevars != 0);
         var = activevars[v];

         hasfloats = hasfloats || (SCIPvarGetType(var) != SCIP_VARTYPE_BINARY &&  SCIPvarGetType(var) != SCIP_VARTYPE_INTEGER);
         hasfloats = hasfloats || !SCIPisIntegral(scip, activevals[v]);
      }

      /* If the constraint has to be written as float type, all discrete variables need to have a float counterpart */
      if( hasfloats )
      {
         for( v = 0; v < nactivevars; v++ )
         {
            SCIP_VAR* var;
            int idx;

            assert(activevars != 0);
            var = activevars[v];
            idx = SCIPvarGetProbindex(var);
            assert( idx >= 0);

            /* If there was no float representation of the variable before, add an auxiliary variable and a conversion constraint */
            if( idx < fznoutput->ndiscretevars && !fznoutput->varhasfloat[idx] )
            {
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

               (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "var float: %s_float;\n", SCIPvarGetName(var));
               SCIP_CALL( appendBuffer(scip, &(fznoutput->varbuffer), &(fznoutput->varbufferlen), &(fznoutput->varbufferpos),buffer) );

               (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "constraint int2float(%s, %s_float);\n", SCIPvarGetName(var), SCIPvarGetName(var));
               SCIP_CALL( appendBuffer(scip, &(fznoutput->castbuffer), &(fznoutput->castbufferlen), &(fznoutput->castbufferpos),buffer) );

               fznoutput->varhasfloat[idx] = TRUE;
            }
         }
      }
   }


   if( SCIPisEQ(scip, lhs, rhs) )
   {
      assert( !SCIPisInfinity(scip, rhs) );

      /* equality constraint */
      SCIP_CALL( printRow(scip, fznoutput, "eq", activevars, activevals, nactivevars, rhs - activeconstant, hasfloats) );
   }
   else
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* print inequality ">=" */
         SCIP_CALL( printRow(scip, fznoutput, "ge", activevars, activevals, nactivevars, lhs - activeconstant, hasfloats) );
      }

      if( !SCIPisInfinity(scip, rhs) )
      {
         /* print inequality "<=" */
         SCIP_CALL( printRow(scip, fznoutput, "le", activevars, activevals, nactivevars, rhs - activeconstant, hasfloats) );
      }
   }

   /* free buffer arrays */
   if( activevars != NULL )
      SCIPfreeBufferArray(scip, &activevars);
   SCIPfreeBufferArray(scip, &activevals);

   return SCIP_OKAY;
}

/* writes problem to a flatzinc conform file, including introduction of several auxiliary variables and constraints */
static
SCIP_RETCODE writeFzn(
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
   FZNOUTPUT fznoutput;                   /* data structure for writing in fzn format */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS* cons;
   const char* conshdlrname;
   SCIP_VAR** consvars;                   /* variables of a specific constraint */
   SCIP_VAR* var;
   SCIP_BOUNDTYPE* boundtypes;            /* indicates whether to which side the variables are bounded */
   SCIP_Real* consvals;                   /* coefficients of a specific constraint */

   int* boundedvars;                      /* variables which are bounded to exactly one side */
   int* floatobjvars;                     /* discrete variables which have a fractional objective coefficient */
   int* intobjvars;                       /* discrete variables which have an integral objective coefficient */

   SCIP_Real lb;                          /* lower bound of some variable */
   SCIP_Real ub;                          /* upper bound of some variable */

   int nboundedvars;                      /* number of variables which are bounded to exactly one side */
   int nconsvars;                         /* number of variables appearing in a specific constraint */
   int nfloatobjvars;                     /* number of discrete variables which have a fractional objective coefficient */
   int nintobjvars;                       /* number of discrete variables which have an integral objective coefficient */
   int c;                                 /* counter for the constraints */
   int v;                                 /* counter for the variables */
   const int ndiscretevars = nbinvars+nintvars;  /* number of discrete variables */

   char varname[SCIP_MAXSTRLEN];          /* buffer for storing variable names */
   char buffer[FZN_BUFFERLEN];            /* buffer for storing auxiliary variables and constraints */
   char buffy[FZN_BUFFERLEN];

   assert( scip != NULL );

   /* print problem statistics as comment to file */
   SCIPinfoMessage(scip, file, "%% SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "%% Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "%% Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "%% Constraints      : %d\n", nconss);

   SCIP_CALL( SCIPallocBufferArray(scip, &boundedvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, nvars) );
   nboundedvars = 0;

   if( nvars > 0 )
      SCIPinfoMessage(scip, file, "\n%%%%%%%%%%%% Problem variables %%%%%%%%%%%%\n");

   /* write all (active) problem variables */
   for( v = 0; v < nvars; v++ )
   {
      var = vars[v];
      assert( var != NULL );
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", SCIPvarGetName(var) );

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

      /* If a variable is bounded to both sides, the bounds are added to the declaration,
       * for variables bounded to exactly one side, an auxiliary constraint will be added later-on.
       */
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
      {
         SCIP_Bool fixed;
         fixed = FALSE;

         if( SCIPisEQ(scip, lb, ub) )
            fixed = TRUE;

         if( v < ndiscretevars )
         {
            assert( SCIPisFeasIntegral(scip, lb) && SCIPisFeasIntegral(scip, ub) );

            if( fixed )
               SCIPinfoMessage(scip, file, "var int: %s = %.f;\n", varname, lb);
            else
               SCIPinfoMessage(scip, file, "var %.f..%.f: %s;\n", lb, ub, varname);
         }
         else
         {
            /* Real valued bounds have to be made type conform */
            if( fixed )
            {
               flattenFloat(scip, lb, buffy);
               SCIPinfoMessage(scip, file, "var float: %s = %s;\n", varname, buffy);
            }
            else
            {
               char buffy2[FZN_BUFFERLEN];

               flattenFloat(scip, lb, buffy);
               flattenFloat(scip, ub, buffy2);
               SCIPinfoMessage(scip, file,  "var %s..%s: %s;\n", buffy, buffy2, varname);
            }
         }
      }
      else
      {
         assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
         assert( v >= nbinvars );

         /* declare the variable without any bound */
         if( v < ndiscretevars )
            SCIPinfoMessage(scip, file, "var int: %s;\n", varname);
         else
            SCIPinfoMessage(scip, file, "var float: %s;\n", varname);

         /* if there is a bound, store the variable and its boundtype for adding a corresponding constraint later-on */
         if( SCIPisInfinity(scip, ub) )
         {
            boundedvars[nboundedvars] = v;
            boundtypes[nboundedvars] = SCIP_BOUNDTYPE_LOWER;
            nboundedvars++;
         }
         if( SCIPisInfinity(scip, -lb) )
         {
            boundedvars[nboundedvars] = v;
            boundtypes[nboundedvars] = SCIP_BOUNDTYPE_UPPER;
            nboundedvars++;
         }
      }
   }

   /* set up the datastructures for the auxiliary int2float variables, the casting constraints and the problem constraints */
   fznoutput.ndiscretevars = ndiscretevars;
   fznoutput.varbufferpos = 0;
   fznoutput.consbufferpos = 0;
   fznoutput.castbufferpos = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &fznoutput.varhasfloat, ndiscretevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fznoutput.varbuffer, FZN_BUFFERLEN) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fznoutput.castbuffer, FZN_BUFFERLEN) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fznoutput.consbuffer, FZN_BUFFERLEN) );
   fznoutput.consbufferlen = FZN_BUFFERLEN;
   fznoutput.varbufferlen = FZN_BUFFERLEN;
   fznoutput.castbufferlen = FZN_BUFFERLEN;

   for( v = 0; v < ndiscretevars; v++ )
      fznoutput.varhasfloat[v] = FALSE;
   fznoutput.varbuffer[0] = '\0';
   fznoutput.consbuffer[0] = '\0';
   fznoutput.castbuffer[0] = '\0';

   /* output all problem constraints */
   for( c = 0; c < nconss; c++ )
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed is written only constraint are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( transformed == SCIPconsIsTransformed(cons) );

      /* By now, only linear, setppc, logicor, knapsack, and varbound constraints can be written.
       * Since they are all linearizable, a linear representation of them is written.
       */
      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, &fznoutput,
               SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
               SCIPgetLhsLinear(scip, cons),  SCIPgetRhsLinear(scip, cons), transformed, TRUE) );
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         consvars = SCIPgetVarsSetppc(scip, cons);
         nconsvars = SCIPgetNVarsSetppc(scip, cons);

         /* Setppc constraints only differ in their lhs/rhs (+- INF or 1) */
         switch( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            SCIP_CALL( printLinearCons(scip, &fznoutput,
                  consvars, NULL, nconsvars, 1.0, 1.0, transformed, FALSE) );
            break;
         case SCIP_SETPPCTYPE_PACKING :
            SCIP_CALL( printLinearCons(scip, &fznoutput,
                  consvars, NULL, nconsvars, -SCIPinfinity(scip), 1.0, transformed, FALSE) );
            break;
         case SCIP_SETPPCTYPE_COVERING :
            SCIP_CALL( printLinearCons(scip, &fznoutput,
                  consvars, NULL, nconsvars, 1.0, SCIPinfinity(scip), transformed, FALSE) );
            break;
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( printLinearCons(scip, &fznoutput,
               SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons),
               1.0, SCIPinfinity(scip), transformed, FALSE) );
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;

         consvars = SCIPgetVarsKnapsack(scip, cons);
         nconsvars = SCIPgetNVarsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
         for( v = 0; v < nconsvars; ++v )
            consvals[v] = (SCIP_Real)weights[v];

         SCIP_CALL( printLinearCons(scip, &fznoutput, consvars, consvals, nconsvars, -SCIPinfinity(scip),
               (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), transformed, FALSE) );

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

         /* Varbound constraints always consist of exactly two variables */
         SCIP_CALL( printLinearCons(scip, &fznoutput,
               consvars, consvals, 2,
               SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons), transformed, TRUE) );

         SCIPfreeBufferArray(scip, &consvars);
         SCIPfreeBufferArray(scip, &consvals);
      }
      else if(  strcmp(conshdlrname, "cumulative") == 0 )
      {
         int* intvals;

         consvars = SCIPgetVarsCumulative(scip, cons);
         nconsvars = SCIPgetNVarsCumulative(scip, cons);

         SCIP_CALL( appendBuffer(scip, &(fznoutput.consbuffer), &(fznoutput.consbufferlen), &(fznoutput.consbufferpos), "cumulative([") );

         for( v = 0; v < nconsvars; ++v )
         {
            if( v < nconsvars - 1)
               (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s, ", SCIPvarGetName(consvars[v]) );
            else
               (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s", SCIPvarGetName(consvars[v]) );

            SCIP_CALL( appendBuffer(scip, &(fznoutput.consbuffer), &(fznoutput.consbufferlen), &(fznoutput.consbufferpos), varname) );
         }

         SCIP_CALL( appendBuffer(scip, &(fznoutput.consbuffer), &(fznoutput.consbufferlen), &(fznoutput.consbufferpos), "], [") );

         intvals = SCIPgetDurationsCumulative(scip, cons);

         for( v = 0; v < nconsvars; ++v )
         {
            if( v < nconsvars - 1)
               (void) SCIPsnprintf(buffy, SCIP_MAXSTRLEN, "%d, ", intvals[v] );
            else
               (void) SCIPsnprintf(buffy, SCIP_MAXSTRLEN, "%d", intvals[v] );

            SCIP_CALL( appendBuffer(scip, &(fznoutput.consbuffer), &(fznoutput.consbufferlen), &(fznoutput.consbufferpos), buffy) );
         }

         SCIP_CALL( appendBuffer(scip, &(fznoutput.consbuffer), &(fznoutput.consbufferlen), &(fznoutput.consbufferpos), "], [") );

         intvals = SCIPgetDemandsCumulative(scip, cons);

         for( v = 0; v < nconsvars; ++v )
         {
            if( v < nconsvars - 1)
               (void) SCIPsnprintf(buffy, SCIP_MAXSTRLEN, "%d, ", intvals[v] );
            else
               (void) SCIPsnprintf(buffy, SCIP_MAXSTRLEN, "%d", intvals[v] );

            SCIP_CALL( appendBuffer(scip, &(fznoutput.consbuffer), &(fznoutput.consbufferlen), &(fznoutput.consbufferpos), buffy) );
         }
         (void) SCIPsnprintf(buffy, SCIP_MAXSTRLEN, "], %d);\n", SCIPgetCapacityCumulative(scip, cons) );

         SCIP_CALL( appendBuffer(scip, &(fznoutput.consbuffer), &(fznoutput.consbufferlen), &(fznoutput.consbufferpos), buffy) );
      }
      else
      {
         SCIPwarningMessage(scip, "constraint handler <%s> cannot print flatzinc format\n", conshdlrname );
      }
   }

   SCIP_CALL( SCIPallocBufferArray(scip,&intobjvars,ndiscretevars) );
   SCIP_CALL( SCIPallocBufferArray(scip,&floatobjvars,nvars) );
   nintobjvars = 0;
   nfloatobjvars = 0;

   /* scan objective function: Which variables have to be put to the float part, which to the int part? */
   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real obj;

      var = vars[v];
      obj = SCIPvarGetObj(var);

      if( !SCIPisZero(scip,obj) )
      {
         /* only discrete variables with integral objective coefficient will be put to the int part of the objective */
         if( v < ndiscretevars && SCIPisIntegral(scip, objscale*obj) )
         {
            intobjvars[nintobjvars] = v;
            SCIPdebugMsg(scip, "variable <%s> at pos <%d,%d> has an integral obj: %f=%f*%f\n",
               SCIPvarGetName(var), nintobjvars, v, obj, objscale, SCIPvarGetObj(var));
            nintobjvars++;
         }
         else
         {
            /* if not happened yet, introduce an auxiliary variable for discrete variables with fractional coefficients */
            if( v < ndiscretevars && !fznoutput.varhasfloat[v] )
            {
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

               (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "var float: %s_float;\n", SCIPvarGetName(var));
               SCIP_CALL( appendBuffer(scip, &(fznoutput.varbuffer), &(fznoutput.varbufferlen), &(fznoutput.varbufferpos),buffer) );

               (void) SCIPsnprintf(buffer, FZN_BUFFERLEN, "constraint int2float(%s, %s_float);\n", SCIPvarGetName(var), SCIPvarGetName(var));
               SCIP_CALL( appendBuffer(scip, &(fznoutput.castbuffer), &(fznoutput.castbufferlen), &(fznoutput.castbufferpos),buffer) );

               fznoutput.varhasfloat[v] = TRUE;
            }

            floatobjvars[nfloatobjvars] = v;
            nfloatobjvars++;
         }
      }
   }

   /* output all created auxiliary variables (float equivalents of discrete variables) */
   if( fznoutput.varbufferpos > 0 )
   {
      SCIPinfoMessage(scip, file, "\n%%%%%%%%%%%% Auxiliary variables %%%%%%%%%%%%\n");
      writeBuffer(scip, file, fznoutput.varbuffer, fznoutput.varbufferpos );
   }

   /* output all int2float casting/conversion constraints */
   if( fznoutput.castbufferpos > 0 )
   {
      SCIPinfoMessage(scip, file, "\n%%%%%%%%%%%% Variable conversions %%%%%%%%%%%%\n");
      writeBuffer(scip, file, fznoutput.castbuffer, fznoutput.castbufferpos );
   }

   if( nboundedvars > 0 )
      SCIPinfoMessage(scip, file, "\n%%%%%%%%%%%% Variable bounds %%%%%%%%%%%%\n");

   /* output all bounds of variables with exactly one bound*/
   for( v = 0; v < nboundedvars; v++ )
   {
      var = vars[boundedvars[v]];

      if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
      {
         if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
            SCIPinfoMessage(scip, file,"constraint int_ge(%s, %.f);\n",SCIPvarGetName(var),
               transformed ? SCIPvarGetLbLocal(var) : SCIPvarGetLbOriginal(var));
         else
         {
            assert( boundtypes[v] == SCIP_BOUNDTYPE_UPPER );
            SCIPinfoMessage(scip, file,"constraint int_le(%s, %.f);\n",SCIPvarGetName(var),
               transformed ? SCIPvarGetUbLocal(var) : SCIPvarGetUbOriginal(var));
         }
      }
      else
      {
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT || SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

         if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
         {
            flattenFloat(scip, transformed ? SCIPvarGetLbLocal(var) : SCIPvarGetLbOriginal(var), buffy);
            SCIPinfoMessage(scip, file,"constraint float_ge(%s, %s);\n", SCIPvarGetName(var), buffy);
         }
         else
         {
            assert( boundtypes[v] == SCIP_BOUNDTYPE_UPPER );
            flattenFloat(scip, transformed ? SCIPvarGetUbLocal(var) : SCIPvarGetUbOriginal(var), buffy);
            SCIPinfoMessage(scip, file,"constraint float_le(%s, %s);\n",SCIPvarGetName(var), buffy);
         }
      }
   }

   /* output all problem constraints */
   if( fznoutput.consbufferpos > 0 )
   {
      SCIPinfoMessage(scip, file, "\n%%%%%%%%%%%% Problem constraints %%%%%%%%%%%%\n");
      writeBuffer(scip, file, fznoutput.consbuffer, fznoutput.consbufferpos );
   }

   SCIPinfoMessage(scip, file, "\n%%%%%%%%%%%% Objective function %%%%%%%%%%%%\n");

   /* If there is at least one variable in the objective function write down the optimization problem, else declare it to be a satisfiability problem */
   if( nintobjvars > 0 || nfloatobjvars > 0 )
   {
      SCIPinfoMessage(scip, file, "solve %s int_float_lin([", objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize" );

      /* first array: coefficients (in float representation) of discrete variables with integral objective coefficient */
      for( v = 0; v < nintobjvars; v++ )
      {
         SCIP_Real obj;
         var = vars[intobjvars[v]];
         obj = objscale*SCIPvarGetObj(var);
         SCIPdebugMsg(scip, "variable <%s> at pos <%d,%d> has an integral obj: %f=%f*%f\n", SCIPvarGetName(var), v, intobjvars[v], obj, objscale, SCIPvarGetObj(var));

         assert( SCIPisIntegral(scip, obj) );
         flattenFloat(scip, obj, buffy);
         SCIPinfoMessage(scip, file, "%s%s", buffy, v < nintobjvars-1 ? ", " : "" );
      }

      /* second array: all other objective coefficients */
      SCIPinfoMessage(scip, file, "], [");
      for( v = 0; v < nfloatobjvars; v++ )
      {
         SCIP_Real obj;
         obj = objscale*SCIPvarGetObj(vars[floatobjvars[v]]);
         flattenFloat(scip, obj, buffy);
         assert( !SCIPisIntegral(scip, obj) || SCIPvarGetType(vars[floatobjvars[v]]) == SCIP_VARTYPE_CONTINUOUS
            || SCIPvarGetType(vars[floatobjvars[v]]) == SCIP_VARTYPE_IMPLINT);
         SCIPinfoMessage(scip, file, "%s%s", buffy, v < nfloatobjvars-1 ? ", " : "" );
      }

      /* potentially add an objective offset */
      if( !SCIPisZero(scip, objoffset) )
      {
         flattenFloat(scip, objoffset, buffy);
         SCIPinfoMessage(scip, file, "%s%s", nfloatobjvars == 0 ? "" : ", ", buffy );
      }

      /* third array: all discrete variables with integral objective coefficient */
      SCIPinfoMessage(scip, file, "], [");
      for( v = 0; v < nintobjvars; v++ )
         SCIPinfoMessage(scip, file, "%s%s", SCIPvarGetName(vars[intobjvars[v]]), v < nintobjvars-1 ? ", " : "" );

      /* fourth array: all other variables with nonzero objective coefficient */
      SCIPinfoMessage(scip, file, "], [");
      for( v = 0; v < nfloatobjvars; v++ )
         SCIPinfoMessage(scip, file, "%s%s%s", SCIPvarGetName(vars[floatobjvars[v]]), floatobjvars[v] < ndiscretevars ? "_float" : "", v < nfloatobjvars-1 ? ", " : "" );

      /* potentially add a 1.0 for the objective offset */
      if( !SCIPisZero(scip, objoffset) )
         SCIPinfoMessage(scip, file, "%s%.1f", nfloatobjvars == 0 ? "" : ", ", 1.0 );
      SCIPinfoMessage(scip, file, "]);\n");
   }
   else
      SCIPinfoMessage(scip, file, "solve satisfy;\n");

   /* free all memory */
   SCIPfreeBufferArray(scip, &fznoutput.castbuffer);
   SCIPfreeBufferArray(scip, &fznoutput.consbuffer);
   SCIPfreeBufferArray(scip, &fznoutput.varbuffer);

   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &boundedvars);
   SCIPfreeBufferArray(scip, &floatobjvars);
   SCIPfreeBufferArray(scip, &intobjvars);
   SCIPfreeBufferArray(scip, &fznoutput.varhasfloat);

   *result = SCIP_SUCCESS;
   return  SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyFzn)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderFzn(scip) );

   return SCIP_OKAY;
}


/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeFzn)
{
   SCIP_READERDATA* readerdata;
   int v;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   /* free all variable array elements */
   for( v = 0; v < readerdata->nvararrays; ++v )
   {
      freeVararray(scip, &readerdata->vararrays[v]);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &readerdata->vararrays, readerdata->vararrayssize);

   /* free reader data */
   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadFzn)
{  /*lint --e{715}*/

   FZNINPUT fzninput;
   int i;

   /* initialize FZN input data */
   fzninput.file = NULL;
   fzninput.linebuf[0] = '\0';
   SCIP_CALL( SCIPallocBufferArray(scip, &fzninput.token, FZN_BUFFERLEN) );
   fzninput.token[0] = '\0';

   for( i = 0; i < FZN_MAX_PUSHEDTOKENS; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(fzninput.pushedtokens[i]), FZN_BUFFERLEN) ); /*lint !e866*/
   }

   fzninput.npushedtokens = 0;
   fzninput.linenumber = 1;
   fzninput.bufpos = 0;
   fzninput.linepos = 0;
   fzninput.objsense = SCIP_OBJSENSE_MINIMIZE;
   fzninput.endline = FALSE;
   fzninput.comment = FALSE;
   fzninput.haserror = FALSE;
   fzninput.valid = TRUE;
   fzninput.vararrays = NULL;
   fzninput.nvararrays = 0;
   fzninput.vararrayssize = 0;
   fzninput.constarrays = NULL;
   fzninput.nconstarrays = 0;
   fzninput.constarrayssize = 0;

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &(fzninput.initialconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicconss", &(fzninput.dynamicconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &(fzninput.dynamiccols)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &(fzninput.dynamicrows)) );

   SCIP_CALL( SCIPhashtableCreate(&fzninput.varHashtable, SCIPblkmem(scip), SCIP_HASHSIZE_NAMES,
         hashGetKeyVar, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );

   SCIP_CALL( SCIPhashtableCreate(&fzninput.constantHashtable, SCIPblkmem(scip), SCIP_HASHSIZE_NAMES,
         hashGetKeyConstant, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fzninput.constants, 10) );

   fzninput.nconstants = 0;
   fzninput.sconstants = 10;

   /* read the file */
   SCIP_CALL( readFZNFile(scip, SCIPreaderGetData(reader), &fzninput, filename) );

   /* free dynamically allocated memory */
   SCIPfreeBufferArrayNull(scip, &fzninput.token);
   for( i = 0; i < FZN_MAX_PUSHEDTOKENS; ++i )
   {
      SCIPfreeBufferArrayNull(scip, &fzninput.pushedtokens[i]);
   }

   /* free buffer memory */
   for( i = 0; i < fzninput.nconstants; ++i )
   {
      SCIPfreeBufferArray(scip, &fzninput.constants[i]->name);
      SCIPfreeBuffer(scip, &fzninput.constants[i]);
   }
   SCIPhashtableFree(&fzninput.varHashtable);
   SCIPhashtableFree(&fzninput.constantHashtable);
   SCIPfreeBufferArray(scip, &fzninput.constants);

   /* free variable arrays */
   for( i = 0; i < fzninput.nvararrays; ++i )
   {
      freeVararray(scip, &fzninput.vararrays[i]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(fzninput.vararrays), fzninput.vararrayssize);

   /* free constant arrays */
   for( i = 0; i < fzninput.nconstarrays; ++i )
   {
      freeConstarray(scip, &(fzninput.constarrays[i]));
   }
   SCIPfreeBlockMemoryArrayNull(scip, &fzninput.constarrays, fzninput.constarrayssize);

   /* evaluate the result */
   if( fzninput.haserror )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteFzn)
{  /*lint --e{715}*/
   if( genericnames )
   {
      SCIP_CALL( writeFzn(scip, file, name, transformed, objsense, objscale, objoffset, vars,
            nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );
   }
   else
   {
      int i;
      SCIP_Bool legal;

      legal = TRUE;

      /* Scan whether all variable names are flatzinc conform */
      for( i = 0; i < nvars; i++ )
      {
         const char* varname;
         size_t length;

         varname = SCIPvarGetName(vars[i]);
         length = strlen(varname);
         legal = isIdentifier(varname);
         if( !legal )
         {
            SCIPwarningMessage(scip, "The name of variable <%d>: \"%s\" is not conform to the fzn standard.\n", i, varname);
            break;
         }

         if( length >= 7 )
            legal = (strncmp(&varname[length-6],"_float",6) != 0);
         if( !legal )
         {
            SCIPwarningMessage(scip, "The name of variable <%d>: \"%s\" ends with \"_float\" which is not supported.\n", i, varname);
            break;
         }
      }

      /* If there is at least one name, which is not conform, use generic names */
      if( legal )
      {
         SCIP_CALL( writeFzn(scip, file, name, transformed, objsense, objscale, objoffset, vars,
               nvars, nbinvars, nintvars, nimplvars, ncontvars, conss, nconss, result) );
      }
      else if( transformed )
      {
         SCIPwarningMessage(scip, "Write transformed problem with generic variable names.\n");
         SCIP_CALL( SCIPprintTransProblem(scip, file, "fzn", TRUE) );
      }
      else
      {
         SCIPwarningMessage(scip, "Write original problem with generic variable names.\n");
         SCIP_CALL( SCIPprintOrigProblem(scip, file, "fzn", TRUE) );
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * reader specific interface methods
 */

/** includes the fzn file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderFzn(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create fzn reader data */
   SCIP_CALL( readerdataCreate(scip, &readerdata) );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyFzn) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeFzn) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadFzn) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteFzn) );

   return SCIP_OKAY;
}

/** print given solution in Flatzinc format w.r.t. the output annotation */
SCIP_RETCODE SCIPprintSolReaderFzn(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   SCIP_VAR** vars;
   VARARRAY** vararrays;
   DIMENSIONS* info;
   VARARRAY* vararray;
   FZNNUMBERTYPE type;
   SCIP_Real solvalue;
   int nvararrays;
   int nvars;
   int i;
   int v;

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   vararrays = readerdata->vararrays;
   nvararrays = readerdata->nvararrays;

   /* sort variable arrays */
   SCIPsortPtr((void**)vararrays, vararraysComp, nvararrays);

   for( i = 0; i < nvararrays; ++i )
   {
      vararray = vararrays[i];
      info = vararray->info;
      vars = vararray->vars;
      nvars = vararray->nvars;
      type =  vararray->type;

      if( info->ndims == 0 )
      {
         solvalue = SCIPgetSolVal(scip, sol, vars[0]);

         SCIPinfoMessage(scip, file, "%s = ", vararray->name);

         printValue(scip, file, solvalue, type);

         SCIPinfoMessage(scip, file, ";\n");
      }
      else
      {

         SCIPinfoMessage(scip, file, "%s = array%dd(", vararray->name, info->ndims);

         for( v = 0; v < info->ndims; ++v )
         {
            SCIPinfoMessage(scip, file, "%d..%d, ", info->lbs[v], info->ubs[v]);
         }

         SCIPinfoMessage(scip, file, "[");

         for( v = 0; v < nvars; ++v )
         {
            if( v > 0)
               SCIPinfoMessage(scip, file, ", ");

            solvalue = SCIPgetSolVal(scip, sol, vars[v]);
            printValue(scip, file, solvalue, type);
         }

         SCIPinfoMessage(scip, file, "]);\n");
      }
   }

   SCIPinfoMessage(scip, file, "----------\n");

   return SCIP_OKAY;
}
