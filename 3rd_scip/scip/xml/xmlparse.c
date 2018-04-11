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

/**@file   xmldef.h
 * @brief  declarations for XML parsing
 * @author Thorsten Koch
 * @author Marc Pfetsch
 *
 * If SPEC_LIKE_SPACE_HANDLING is not defined, all LF,CR will be changed into spaces and from a
 * sequence of spaces only one will be used.
 *
 * @todo Implement possibility to avoid the construction of parsing information for certain tags
 * (and their children). For solution files this would avoid parsing the constraints section.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <blockmemshell/memory.h>

#include "xml.h"
#include "xmldef.h"


#include <sys/types.h>
#ifdef WITH_ZLIB
#if defined(_WIN32) || defined(_WIN64)
#define R_OK _A_RDONLY
#define access _access
#include <io.h>
#else
#include <unistd.h>
#endif
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>


#define NAME_EXT_SIZE 128
#define ATTR_EXT_SIZE 4096
#define DATA_EXT_SIZE 4096
#define LINE_BUF_SIZE 8192

#define xmlError(a, b) xmlErrmsg(a, b, FALSE, __FILE__, __LINE__)


/* forward declarations */
typedef struct parse_stack_struct PSTACK;
typedef struct parse_pos_struct   PPOS;

/** state of the parser */
enum parse_state_enum
{
   XML_STATE_ERROR,
   XML_STATE_BEFORE,
   XML_STATE_IN_TAG,
   XML_STATE_PCDATA,
   XML_STATE_EOF
};
typedef enum   parse_state_enum   PSTATE;

/** Stack as a (singly) linked list. The top element is the current node. */
struct parse_stack_struct
{
   XML_NODE*             node;
   PSTACK*               next;
};

/** Store the current position in the file and the state of the parser. */
struct parse_pos_struct
{
   const char*           filename;
   FPTYPE                fp;
   char                  buf[LINE_BUF_SIZE];
   int                   pos;
   int                   lineno;
   int                   nextsym;
   int                   lastsym;
   PSTATE                state;
   PSTACK*               top;
};


/** output error message with corresponding line and position */
static void xmlErrmsg(
   PPOS*                 ppos,
   const char*           msg,
   XML_Bool              msg_only,
   const char*           file,
   int                   line
   )
{
#ifndef NDEBUG
   int ret;
   assert( ppos != NULL );

   if ( ! msg_only )
   {
      ret = fprintf(stderr, "%s(%d) Error in file %s line %d\n", file, line, ppos->filename, ppos->lineno);
      assert(ret >= 0);

      ret = fprintf(stderr, "%s", ppos->buf);
      assert(ret >= 0);

      if ( strchr(ppos->buf, '\n') == NULL )
      {
         int retc;

         retc = fputc('\n', stderr);
         assert(retc != EOF);
      }

      ret = fprintf(stderr, "%*s\n", ppos->pos, "^");
      assert(ret >= 0);
   }
   ret = fprintf(stderr, "%s\n\n", msg);
   assert(ret >= 0);

#else

   if ( ! msg_only )
   {
      (void) fprintf(stderr, "%s(%d) Error in file %s line %d\n", file, line, ppos->filename, ppos->lineno);

      (void) fprintf(stderr, "%s", ppos->buf);

      if ( strchr(ppos->buf, '\n') == NULL )
      {
         (void) fputc('\n', stderr);
      }

      (void) fprintf(stderr, "%*s\n", ppos->pos, "^");
   }
   (void) fprintf(stderr, "%s\n\n", msg);
#endif
}


/** Push new element on the parse stack.
 *
 *  TRUE if it worked, FAILURE otherwise.
 */
static
XML_Bool pushPstack(
   PPOS*                 ppos,
   XML_NODE*             node
   )
{
   PSTACK* p;

   assert(ppos != NULL);
   assert(node != NULL);

   debugMessage("Pushing %s\n", node->name);

   ALLOC_FALSE( BMSallocMemory(&p) );
   assert(p != NULL);

   p->node   = node;
   p->next   = ppos->top;
   ppos->top = p;

   return TRUE;
}

/** returns top element on stack (which has to be present) */
static XML_NODE* topPstack(
   const PPOS*           ppos
   )
{
   assert(ppos      != NULL);
   assert(ppos->top != NULL);

   return ppos->top->node;
}

/** remove top element from stack and deletes it
 *
 *  TRUE if ok, FALSE otherwise
 */
static
XML_Bool popPstack(
   PPOS*                 ppos                /**< input stream position */
   )
{
   PSTACK* p;
   XML_Bool result;

   assert(ppos != NULL);

   if ( ppos->top == NULL )
   {
      xmlError(ppos, "Stack underflow");
      result = FALSE;
   }
   else
   {
      result = TRUE;
      p = ppos->top;
      ppos->top = p->next;

      debugMessage("Poping %s\n", p->node->name);
      BMSfreeMemory(&p);
   }
   return result;
}

/** remove complete stack */
static
void clearPstack(
   PPOS*                 ppos
   )
{
   assert(ppos != NULL);

   while ( ppos->top != NULL )
      (void) popPstack(ppos);
}

/** Returns the next character from the input buffer and fills the buffer if it is empty (similar to fgetc()). */
static
int mygetc(
   PPOS*                 ppos
   )
{
   assert(ppos      != NULL);
   assert(ppos->fp  != NULL);
   assert(ppos->pos <  LINE_BUF_SIZE);

   if ( ppos->buf[ppos->pos] == '\0' )
   {
#ifdef SCIP_DISABLED_CODE
      /* the low level function gzread/fread used below seem to be faster */
      if ( NULL == FGETS(ppos->buf, sizeof(ppos->buf), ppos->fp) )
         return EOF;
#else
      size_t len = (size_t) FREAD(ppos->buf, sizeof(ppos->buf) - 1, ppos->fp); /*lint !e571 !e747*/

      if( len == 0 || len > sizeof(ppos->buf) - 1 )
         return EOF;

      ppos->buf[len] = '\0';
#endif
      ppos->pos = 0;
   }
   return (unsigned char)ppos->buf[ppos->pos++];
}


#ifdef SPEC_LIKE_SPACE_HANDLING
/** Read input from fp_in.
 *
 * If there is a LF, CR, CR/LF, or LF/CR it returns exactly on LF.  Also counts the number of
 * characters.
 */
static
int getsymbol(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos != NULL);

   if ( ppos->nextsym == 0 )
      c = mygetc(ppos);
   else
   {
      c = ppos->nextsym;
      ppos->nextsym = 0;
   }
   assert(ppos->nextsym == 0);

   if (((c == '\n') && (ppos->lastsym == '\r')) || ((c == '\r') && (ppos->lastsym == '\n')))
      c = mygetc(ppos);

   ppos->lastsym = c;

   if ( c == '\r' )
      c = '\n';

   if ( c == '\n' )
      ++ppos->lineno;

   return c;
}
#else
/** Read input from fp_in (variant).
 *
 *  Here we convert all LF or CR into SPACE and return maximally one SPACE after the other.
 *
 *  @note This function counts lines differently. On systems that have only one '\\r' as line feed
 *  (MAC) it does not count correctly.
 */
static
int getsymbol(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos != NULL);

   do
   {
      if ( ppos->nextsym == 0 )
         c = mygetc(ppos);
      else
      {
         c = ppos->nextsym;
         ppos->nextsym = 0;
      }
      assert(ppos->nextsym == 0);

      if ( c == '\n' )
         ++ppos->lineno;

      if ((c == '\n') || (c == '\r'))
         c = ' ';
   } while((c == ' ') && (ppos->lastsym == c));

   ppos->lastsym = c;

   debugMessage("[%c]\n", c);

   return c;
}
#endif

/** Reinserts a character into the input stream */
static
void ungetsymbol(
   PPOS*                 ppos,
   int                   c
   )
{
   assert(ppos          != NULL);
   assert(ppos->nextsym == 0);

   ppos->nextsym = c;
}

/** Skip all spaces and return the next non-space character or EOF */
static
int skipSpace(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos != NULL);

   do
   {
      c = getsymbol(ppos);
   }
   while(isspace(c));

   return c;
}

/** Get name of a TAG or attribute from the input stream.
 *
 *  Either it returns a pointer to allocated memory which contains the name or it returns NULL if
 *  there is some error.
 */
static
char* getName(
   PPOS*                 ppos
   )
{
   char* name = NULL;
   size_t size = 0;
   size_t len  = 0;
   int c;

   assert(ppos != NULL);

   c = getsymbol(ppos);

   if ( ! isalpha(c) && (c != '_') && (c != ':') )
   {
      xmlError(ppos, "Name starting with illegal charater");
      return NULL;
   }

   /* The following is wrong: Here almost all characters that we casted to unicode are feasible */
   while ( isalnum(c) || (c == '_') || (c == ':') || (c == '.') || (c == '-') )
   {
      if ( len + 1 >= size )
      {
         size += NAME_EXT_SIZE;

         if ( name == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&name, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&name, size) );
         }
      }
      assert(name != NULL);
      assert(size > len);

      name[len++] = (char)c;

      c = getsymbol(ppos);
   }
   if ( c != EOF )
      ungetsymbol(ppos, c);

   assert(name != NULL);

   if ( len == 0 )
   {
      BMSfreeMemoryArray(&name);
      name = NULL;
   }
   else
      name[len] = '\0';

   return name;
}

/** Read the value of an attribute from the input stream.
 *
 *  The value has to be between two " or ' (the other character is then valid as well). The function
 *  returns a pointer to allocated memory containing the value or it returns NULL in case of an
 *  error.
 */
static
char* getAttrval(
   PPOS*                 ppos
   )
{
   char*  attr = NULL;
   int    c;
   int    stop;
   size_t len = 0;
   size_t size = 0;

   assert(ppos != NULL);

   /* The following is not allowed according to the specification (the value has to be directly
    * after the equation sign). */
   c = skipSpace(ppos);

   if ( (c != '"') && (c != '\'') )
   {
      xmlError(ppos, "Atribute value does not start with \" or \'");
      return NULL;
   }
   stop = c;

   for(;;)
   {
      if ( len == size )
      {
         size += ATTR_EXT_SIZE;

         if ( attr == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&attr, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&attr, size) );
         }
      }
      assert(attr != NULL);
      assert(size >  len);

      c = getsymbol(ppos);

      if ( (c == stop) || (c == EOF) )
         break;

      attr[len++] = (char)c;
   }

   if ( c != EOF )
      attr[len] = '\0';
   else
   {
      BMSfreeMemoryArray(&attr);
      attr = NULL;
   }
   return attr;
}

/** Skip comment
 *
 *  Return FALSE if an error occurs.
 */
static
XML_Bool doComment(
   PPOS*                 ppos
   )
{
   XML_Bool result = TRUE;
   int c;
   int state = 0;

   assert(ppos != NULL);

   for(;;)
   {
      c = getsymbol(ppos);

      if ( c == EOF )
         break;

      if ( (c == '>') && (state >= 2) )
         break;

      state = (c == '-') ? state + 1 : 0;
   }
   if ( c == EOF )
   {
      xmlError(ppos, "Unexpected EOF in comment");
      result = FALSE;
   }
   return result;
}

/** Handles a CDATA section.
 *
 *  Returns a pointer to allocated memory containing the data of this section or NULL in case of an
 *  error.
 */
static
char* doCdata(
   PPOS*                 ppos
   )
{
   char* data  = NULL;
   size_t size  = 0;
   size_t len   = 0;
   int state = 0;
   int c;

   assert(ppos != NULL);

   for(;;)
   {
      c = getsymbol(ppos);

      if ( c == EOF )
         break;

      if ( c == ']' )
         state++;
      else
         if ( (c == '>') && (state >= 2) )
            break;
         else
            state = 0;

      if ( len == size )
      {
         size += DATA_EXT_SIZE;

         if ( data == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&data, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&data, size) );
         }
      }
      assert(data != NULL);
      assert(size >  len);

      data[len++] = (char)c;
   }
   assert(data != NULL);

   /*lint --e{527}*/
   if ( c != EOF )
   {
      assert(len  >= 2);
      assert(data != NULL);

      data[len - 2] = '\0'; /*lint !e413*/
   }
   else
   {
      BMSfreeMemoryArray(&data);
      data = NULL;
      xmlError(ppos, "Unexpected EOF in CDATA");
   }
   return data;
}

/** Handle processing instructions (skipping) */
static
void handlePi(
   PPOS*                 ppos
   )
{
   int c;

   assert(ppos        != NULL);
   assert(ppos->state == XML_STATE_BEFORE);

   do
   {
      c = getsymbol(ppos);
   }
   while ( (c != EOF) && (c != '>') );

   if ( c != EOF )
      ppos->state = XML_STATE_PCDATA;
   else
   {
      xmlError(ppos, "Unexpected EOF in PI");
      ppos->state = XML_STATE_ERROR;
   }
}

/** Handles declarations that start with a <!.
 *
 *  This includes comments. Does currenlty not work very well, because of DTDs.
 */
static
void handleDecl(
   PPOS*                 ppos
   )
{
   enum XmlSection
   {
      IS_COMMENT,
      IS_ATTLIST,
      IS_DOCTYPE,
      IS_ELEMENT,
      IS_ENTITY,
      IS_NOTATION,
      IS_CDATA
   };
   typedef enum XmlSection XMLSECTION;

   static struct
   {
      const char* name;
      XMLSECTION  what;
   } key[] =
   {
      { "--",       IS_COMMENT  },
      { "ATTLIST",  IS_ATTLIST  },
      { "DOCTYPE",  IS_DOCTYPE  },
      { "ELEMENT",  IS_ELEMENT  },
      { "ENTITY",   IS_ENTITY   },
      { "NOTATION", IS_NOTATION },
      { "[CDATA[",  IS_CDATA    }
   };
   XML_NODE* node;
   char* data;
   int c;
   int k = 0;
   int beg = 0;
   int end;

   assert(ppos        != NULL);
   assert(ppos->state == XML_STATE_BEFORE);

   end = (int) (sizeof(key) / sizeof(key[0])) - 1;
   do
   {
      c = getsymbol(ppos);

      for(; (beg <= end) && (c != key[beg].name[k]); beg++)
         ;
      for(; (end >= beg) && (c != key[end].name[k]); end--)
         ;
      k++;
   } while(beg < end);

   if ( beg != end )
   {
      xmlError(ppos, "Unknown declaration");

      while ( (c != EOF) && (c != '>') )
         c = getsymbol(ppos);
   }
   else
   {
      assert(beg == end);
      assert(beg <  (int)(sizeof(key) / sizeof(*key)));

      switch(key[beg].what)
      {
      case IS_COMMENT :
         if ( ! doComment(ppos) )
            ppos->state = XML_STATE_ERROR;
         break;
      case IS_CDATA :
         if ( (data = doCdata(ppos)) == NULL )
            ppos->state = XML_STATE_ERROR;
         else
         {
            if ( NULL == (node = xmlNewNode("#CDATA", ppos->lineno)) )
            {
               xmlError(ppos, "Can't create new node");
               ppos->state = XML_STATE_ERROR;
            }
            else
            {
               BMSduplicateMemoryArray(&node->data, data, strlen(data)+1);
               BMSfreeMemoryArray(&data);
               xmlAppendChild(topPstack(ppos), node);
            }
         }
         break;
      case IS_ATTLIST :
      case IS_ELEMENT :
      case IS_NOTATION :
      case IS_ENTITY :
      case IS_DOCTYPE :
         break;
      default :
         abort();
      }
   }
}

/** Handle end tag */
static
void handleEndtag(
   PPOS*                 ppos
   )
{
   char* name;
   int   c;

   assert(ppos != NULL);

   if ( (name = getName(ppos)) == NULL )
      xmlError(ppos, "Missing name in endtag");
   else
   {
      c = skipSpace(ppos);

      if ( c != '>' )
      {
         xmlError(ppos, "Missing '>' in endtag");
         ppos->state = XML_STATE_ERROR;
      }
      else
      {
         if ( strcmp(name, topPstack(ppos)->name) )
         {
            xmlError(ppos, "Name of endtag does not match starttag");
            ppos->state = XML_STATE_ERROR;
         }
         else
         {
            if ( popPstack(ppos) )
               ppos->state = XML_STATE_PCDATA;
            else
               ppos->state = XML_STATE_ERROR;
         }
      }

      BMSfreeMemoryArray(&name);
   }
}

/** Handle start tag */
static
void handleStarttag(
   PPOS*                 ppos
   )
{
   XML_NODE* node;
   char* name;

   assert(ppos != NULL);

   name = getName(ppos);
   if ( name == NULL )
   {
      xmlError(ppos, "Missing name in tagstart");
      ppos->state = XML_STATE_ERROR;
   }
   else
   {
      node = xmlNewNode(name, ppos->lineno);
      if ( node == NULL )
      {
         xmlError(ppos, "Can't create new node");
         ppos->state = XML_STATE_ERROR;
      }
      else
      {
         xmlAppendChild(topPstack(ppos), node);

         if ( pushPstack(ppos, node) )
            ppos->state = XML_STATE_IN_TAG;
         else
            ppos->state = XML_STATE_ERROR;
      }
      BMSfreeMemoryArray(&name);
   }
}

/** Checks for next tag */
static
void procBefore(
   PPOS*                 ppos                /**< input stream position */
   )
{
   int c;

   assert(ppos        != NULL);
   assert(ppos->state == XML_STATE_BEFORE);

   c = skipSpace(ppos);

   if ( c != '<' )
   {
      xmlError(ppos, "Expecting '<'");
      ppos->state = XML_STATE_ERROR;
   }
   else
   {
      c = getsymbol(ppos);

      switch(c)
      {
      case EOF :
         xmlError(ppos, "Unexpected EOF");
         ppos->state = XML_STATE_ERROR;
         break;
      case '!' :
         handleDecl(ppos);
         break;
      case '?' :
         handlePi(ppos);
         break;
      case '/' :
         handleEndtag(ppos);
         break;
      default :
         ungetsymbol(ppos, c);
         handleStarttag(ppos);
         break;
      }
   }
}

/** Process tag */
static
void procInTag(
   PPOS*                 ppos                /**< input stream position */
   )
{
   XML_ATTR* attr;
   int     c;
   XML_Bool empty = FALSE;
   char*   name;
   char*   value;

   assert(ppos        != NULL);
   assert(ppos->state == XML_STATE_IN_TAG);

   c = skipSpace(ppos);

   if ( (c == '/') || (c == '>') || (c == EOF) )
   {
      if ( c == '/' )
      {
         empty = TRUE;
         c = getsymbol(ppos);
      }

      if ( c == EOF )
      {
         xmlError(ppos, "Unexpected EOF while in a tag");
         ppos->state = XML_STATE_ERROR;
      }

      if ( c == '>' )
      {
         ppos->state = XML_STATE_PCDATA;

         if (empty && ! popPstack(ppos))
            ppos->state = XML_STATE_ERROR;
      }
      else
      {
         xmlError(ppos, "Expected tag end marker '>'");
         ppos->state = XML_STATE_ERROR;
      }
   }
   else
   {
      ungetsymbol(ppos, c);

      name = getName(ppos);
      if ( name == NULL )
      {
         xmlError(ppos, "No name for attribute");
         ppos->state = XML_STATE_ERROR;
      }
      else
      {
         c = skipSpace(ppos);

         if ( (c != '=') || ((value = getAttrval(ppos)) == NULL) )
         {
            xmlError(ppos, "Missing attribute value");
            ppos->state = XML_STATE_ERROR;
            BMSfreeMemoryArray(&name);
         }
         else
         {
            attr = xmlNewAttr(name, value);
            if ( attr == NULL )
            {
               xmlError(ppos, "Can't create new attribute");
               ppos->state = XML_STATE_ERROR;
            }
            else
            {
               xmlAddAttr(topPstack(ppos), attr);
            }
            BMSfreeMemoryArray(&name);
            BMSfreeMemoryArray(&value);
         }
      }
   }
}

/* Handles PCDATA */
static
void procPcdata(
   PPOS*                 ppos                /**< input stream position */
   )
{
   XML_NODE* node;
   char*   data   = NULL;
   size_t  size   = 0;
   size_t  len    = 0;
   int     c;

   assert(ppos        != NULL);
   assert(ppos->state == XML_STATE_PCDATA);

#ifndef SPEC_LIKE_SPACE_HANDLING
   c = skipSpace(ppos);
   if ( c != EOF )
      ungetsymbol(ppos, c);
#endif
   c = getsymbol(ppos);

   while ( (c != EOF) && (c != '<') )
   {
      if ( len + 1 >= size ) /* leave space for terminating '\0' */
      {
         size += DATA_EXT_SIZE;

         if ( data == NULL )
         {
            ALLOC_ABORT( BMSallocMemoryArray(&data, size) );
         }
         else
         {
            ALLOC_ABORT( BMSreallocMemoryArray(&data, size) );
         }
      }
      assert(data != NULL);
      assert(size > len + 1);

      data[len++] = (char)c;

      c = getsymbol(ppos);
   }
   if ( data == NULL )
   {
      if ( c == EOF )
         ppos->state = XML_STATE_EOF;
      else
      {
         assert(c == '<');
         ppos->state = XML_STATE_BEFORE;
         ungetsymbol(ppos, c);
      }
   }
   else
   {
      assert(len < size);
      data[len] = '\0';

      if ( c == EOF )
         ppos->state = XML_STATE_ERROR;
      else
      {
         ungetsymbol(ppos, c);

         node = xmlNewNode("#PCDATA", ppos->lineno);
         if ( node == NULL )
         {
            xmlError(ppos, "Can't create new node");
            ppos->state = XML_STATE_ERROR;
         }
         else
         {
            BMSduplicateMemoryArray(&node->data, data, strlen(data)+1);
            xmlAppendChild(topPstack(ppos), node);
            ppos->state = XML_STATE_BEFORE;
         }
      }

      BMSfreeMemoryArray(&data);
   }
}

/** Parse input stream */
static
XML_Bool xmlParse(
   PPOS*                 ppos                /**< input stream position */
   )
{
   XML_Bool ok = TRUE;

   while (ok)
   {
      debugMessage("state=%d\n", ppos->state);

      switch (ppos->state)
      {
      case XML_STATE_BEFORE :
         procBefore(ppos);
         break;
      case XML_STATE_IN_TAG :
         procInTag(ppos);
         break;
      case XML_STATE_PCDATA :
         procPcdata(ppos);
         break;
      case XML_STATE_EOF :
         ok = FALSE;
         break;
      case XML_STATE_ERROR :
         ok = FALSE;
         break;
      default :
         xmlError(ppos, "Internal Error, illegal state");
         ok = FALSE;
      }
   }
   return (ppos->state == XML_STATE_EOF);
}

/** Parse file */
XML_NODE* xmlProcess(
   const char*           filename            /**< XML file name */
   )
{
   PPOS      ppos;
   XML_NODE* node = NULL;
   XML_ATTR* attr;
   XML_Bool  result = FALSE;
   char*     myfilename;
   size_t    filenamelen;

   /* allocate space and copy filename (possibly modified below) in two steps in order to satisfy valgrind */
   assert( filename != NULL );
   filenamelen = strlen(filename);
   if ( BMSallocMemoryArray(&myfilename, filenamelen + 5) == NULL )
      return NULL;
   BMScopyMemoryArray(myfilename, filename, filenamelen + 1);

#ifdef WITH_ZLIB
   if ( access(filename, R_OK) != 0 )
   {
      strcat(myfilename, ".gz");

      /* If .gz also does not work, revert to the old name
       * to get a better error message.
       */
      if ( access(myfilename, R_OK) != 0 )
         strcpy(myfilename, filename);
   }
#endif
   ppos.fp = FOPEN(myfilename, "r");
   if ( ppos.fp == NULL )
      perror(myfilename);
   else
   {
      ppos.filename = myfilename;
      ppos.buf[0]   = '\0';
      ppos.pos      = 0;
      ppos.lineno   = 1;
      ppos.nextsym  = 0;
      ppos.lastsym  = 0;
      ppos.state    = XML_STATE_BEFORE;
      ppos.top      = NULL;

      node = xmlNewNode("#ROOT", ppos.lineno);
      if ( node == NULL )
      {
         xmlError(&ppos, "Can't create new node");
      }
      else
      {
         attr = xmlNewAttr("filename", myfilename);
         if ( attr == NULL )
            xmlError(&ppos, "Can't create new attribute");
         else
         {
            xmlAddAttr(node, attr);

            /* push root node on stack and start to process */
            if ( pushPstack(&ppos, node) )
            {
               result = xmlParse(&ppos);

               clearPstack(&ppos);
            }
         }
      }

      if ( ! result && (node != NULL) )
      {
         xmlErrmsg(&ppos, "Parsing error, processing stopped", TRUE, __FILE__, __LINE__);
         xmlFreeNode(node);
         node = NULL;
      }
      if ( FCLOSE(ppos.fp) )
         perror(myfilename);
   }
   BMSfreeMemoryArray(&myfilename);

   return node;
}






/*----------------------------------------------------------------------------------------------*/


/** create new node */
XML_NODE* xmlNewNode(
   const char*           name,
   int                   lineno
   )
{
   XML_NODE* n = NULL;

   assert(name != NULL);

   if ( BMSallocMemory(&n) != NULL )
   {
      BMSclearMemory(n);
      BMSduplicateMemoryArray(&n->name, name, strlen(name)+1);
      n->lineno = lineno;
   }
   return n;
}

/** create new attribute */
XML_ATTR* xmlNewAttr(
   const char*           name,
   const char*           value
   )
{
   XML_ATTR* a = NULL;

   assert(name  != NULL);
   assert(value != NULL);

   if ( BMSallocMemory(&a) != NULL )
   {
      BMSclearMemory(a);
      BMSduplicateMemoryArray(&a->name, name, strlen(name)+1);
      BMSduplicateMemoryArray(&a->value, value, strlen(value)+1);
   }
   return a;
}

/** add attribute */
void xmlAddAttr(
   XML_NODE*             n,
   XML_ATTR*             a
   )
{
   assert(n != NULL);
   assert(a != NULL);

   a->next = n->attrlist;
   n->attrlist = a;
}

/** append child node */
void xmlAppendChild(
   XML_NODE*             parent,
   XML_NODE*             child
   )
{
   assert(parent != NULL);
   assert(child  != NULL);

   child->parent = parent;
   child->prevsibl = parent->lastchild;
   child->nextsibl = NULL;
   parent->lastchild = child;

   if ( child->prevsibl != NULL )
      child->prevsibl->nextsibl = child;

   if ( parent->firstchild == NULL )
      parent->firstchild = child;
}

/** free attribute */
static
void xmlFreeAttr(
   XML_ATTR*             attr
   )
{
   XML_ATTR* a;

   /* Note: use an iterative implementation instead of a recursive one; the latter is much slower for large instances
    * and might overflow the heap. */
   a = attr;
   while (a != NULL)
   {
      XML_ATTR* b;
      b = a->next;

      assert(a->name  != NULL);
      assert(a->value != NULL);

      BMSfreeMemoryArray(&a->name);
      BMSfreeMemoryArray(&a->value);
      BMSfreeMemory(&a);
      a = b;
   }
}

/** free node */
void xmlFreeNode(
   XML_NODE*             node
   )
{
   XML_NODE* n;

   if ( node == NULL )
      return;

   /* Free data from back to front (because free is faster this way). */
   /* Note: use an iterative implementation instead of a recursive one; the latter is much slower for large instances
    * and might overflow the heap. */
   n = node->lastchild;
   while ( n != NULL )
   {
      XML_NODE* m;
      m = n->prevsibl;
      xmlFreeNode(n);
      n = m;
   }

   xmlFreeAttr(node->attrlist);

   if ( node->data != NULL )
   {
      BMSfreeMemoryArray(&node->data);
   }
   assert(node->name != NULL);

   BMSfreeMemoryArray(&node->name);
   BMSfreeMemory(&node);
}

/** output node */
void xmlShowNode(
   const XML_NODE*       root
   )
{
   const XML_NODE* n;
   const XML_ATTR* a;

   assert(root != NULL);

   for (n = root; n != NULL; n = n->nextsibl)
   {
      infoMessage("Name: %s\n", n->name);
      infoMessage("Line: %d\n", n->lineno);
      infoMessage("Data: %s\n", (n->data != NULL) ? n->data : "***");

      for (a = n->attrlist; a != NULL; a = a->next)
         infoMessage("Attr: %s = [%s]\n", a->name, a->value);

      if ( n->firstchild != NULL )
      {
         infoMessage("->\n");
         xmlShowNode(n->firstchild);
         infoMessage("<-\n");
      }
   }
}

/** get attribute value */
const char* xmlGetAttrval(
   const XML_NODE*       node,
   const char*           name
   )
{
   XML_ATTR* a;

   assert(node != NULL);
   assert(name != NULL);

   for (a = node->attrlist; a != NULL; a = a->next)
   {
      if ( ! strcmp(name, a->name) )
         break;
   }

#ifdef SCIP_DEBUG
   if (a == NULL)
      infoMessage("Error: Attribute %s in TAG <%s> not found\n", name, node->name);
#endif

   return (a == NULL) ? NULL : a->value;
}

/** return first node */
const XML_NODE* xmlFirstNode(
   const XML_NODE*       node,
   const char*           name
   )
{
   const XML_NODE* n;

   assert(node != NULL);
   assert(name != NULL);

   for (n = node; n != NULL; n = n->nextsibl)
   {
      if ( ! strcmp(name, n->name) )
         break;
   }

   return n;
}

/** return next node */
const XML_NODE* xmlNextNode(
   const XML_NODE*       node,
   const char*           name
   )
{
   assert(node != NULL);
   assert(name != NULL);

   return (node->nextsibl == NULL) ? NULL : xmlFirstNode(node->nextsibl, name);
}

/** find node */
const XML_NODE* xmlFindNode(
   const XML_NODE*       node,
   const char*           name
   )
{
   const XML_NODE* n;
   const XML_NODE* r;

   assert(node != NULL);
   assert(name != NULL);

   if ( ! strcmp(name, node->name) )
      return node;

   for (n = node->firstchild; n != NULL; n = n->nextsibl)
   {
      r = xmlFindNode(n, name);
      if ( r != NULL )
         return r;
   }

   return NULL;
}

/** find node with bound on the depth */
const XML_NODE* xmlFindNodeMaxdepth(
   const XML_NODE*       node,               /**< current node - use start node to begin */
   const char*           name,               /**< name of tag to search for */
   int                   depth,              /**< current depth - start with 0 for root */
   int                   maxdepth            /**< maximal depth */
   )
{
   const XML_NODE* n;
   const XML_NODE* r;

   assert(node != NULL);
   assert(name != NULL);

   if ( ! strcmp(name, node->name) )
      return node;

   if ( depth < maxdepth )
   {
      for (n = node->firstchild; n != NULL; n = n->nextsibl)
      {
         r = xmlFindNodeMaxdepth(n, name, depth+1, maxdepth);
         if ( r != NULL )
            return r;
      }
   }

   return NULL;
}

/** return next sibling */
const XML_NODE* xmlNextSibl(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->nextsibl;
}

/** return previous sibling */
const XML_NODE* xmlPrevSibl(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->prevsibl;
}

/** return first child */
const XML_NODE* xmlFirstChild(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->firstchild;
}

/** return last child */
const XML_NODE* xmlLastChild(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->lastchild;
}

/** return name of node */
const char* xmlGetName(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->name;
}

/** get line number */
int xmlGetLine(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->lineno;
}

/** get data */
const char* xmlGetData(
   const XML_NODE*       node
   )
{
   assert(node != NULL);

   return node->data;
}

/** find PCDATA */
const char* xmlFindPcdata(
   const XML_NODE*       node,
   const char*           name
   )
{
   const XML_NODE* n;

   assert(node != NULL);
   assert(name != NULL);

   n = xmlFindNode(node, name);
   if ( n == NULL )
      return NULL;

   if ( ! strcmp(n->firstchild->name, "#PCDATA") )
      return n->firstchild->data;

   return NULL;
}
