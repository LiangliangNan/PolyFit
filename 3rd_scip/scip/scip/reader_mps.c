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

/**@file   reader_mps.c
 * @brief  (extended) MPS file reader
 * @author Thorsten Koch
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Stefan Heinz
 * @author Stefan Vigerske
 * @author Michael Winkler
 *
 * This reader/writer handles MPS files in extended MPS format, as it
 * is used by CPLEX. In the extended format the limits on variable
 * name lengths and coefficients are considerably relaxed. The columns
 * in the format are then separated by whitespaces.
 *
 * @todo Check whether constructing the names for aggregated constraint yields name clashes (aggrXXX).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/reader_mps.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_and.h"
#include "scip/cons_sos1.h"
#include "scip/cons_sos2.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_soc.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/pub_misc.h"

#define READER_NAME             "mpsreader"
#define READER_DESC             "file reader for MIQPs in IBM's Mathematical Programming System format"
#define READER_EXTENSION        "mps"

#define DEFAULT_LINEARIZE_ANDS         TRUE  /**< should possible \"and\" constraint be linearized when writing the mps file? */
#define DEFAULT_AGGRLINEARIZATION_ANDS TRUE  /**< should an aggregated linearization for and constraints be used? */

/*
 * mps reader internal methods
 */

#define MPS_MAX_LINELEN  1024
#define MPS_MAX_NAMELEN   256
#define MPS_MAX_VALUELEN   26
#define MPS_MAX_FIELDLEN   20

#define PATCH_CHAR    '_'
#define BLANK         ' '

/** MPS reading data */
struct SCIP_ReaderData
{
   SCIP_Bool             linearizeands;
   SCIP_Bool             aggrlinearizationands;
};

/** enum containing all mps sections */
enum MpsSection
{
   MPS_NAME,
   MPS_OBJSEN,
   MPS_OBJNAME,
   MPS_ROWS,
   MPS_USERCUTS,
   MPS_LAZYCONS,
   MPS_COLUMNS,
   MPS_RHS,
   MPS_RANGES,
   MPS_BOUNDS,
   MPS_SOS,
   MPS_QUADOBJ,
   MPS_QMATRIX,
   MPS_QCMATRIX,
   MPS_INDICATORS,
   MPS_ENDATA
};
typedef enum MpsSection MPSSECTION;

/** mps input structure */
struct MpsInput
{
   MPSSECTION            section;
   SCIP_FILE*            fp;
   int                   lineno;
   SCIP_OBJSENSE         objsense;
   SCIP_Bool             haserror;
   char                  buf[MPS_MAX_LINELEN];
   const char*           f0;
   const char*           f1;
   const char*           f2;
   const char*           f3;
   const char*           f4;
   const char*           f5;
   char                  probname[MPS_MAX_NAMELEN];
   char                  objname [MPS_MAX_NAMELEN];
   SCIP_Bool             initialconss;       /**< should model constraints be marked as initial? */
   SCIP_Bool             dynamicconss;       /**< should model constraints be subject to aging? */
   SCIP_Bool             dynamiccols;        /**< should columns be added and removed dynamically to the LP? */
   SCIP_Bool             dynamicrows;        /**< should rows be added and removed dynamically to the LP? */
   SCIP_Bool             isinteger;
   SCIP_Bool             isnewformat;
};
typedef struct MpsInput MPSINPUT;

/** sparse matrix representation */
struct SparseMatrix
{
   SCIP_Real*            values;             /**< matrix element */
   SCIP_VAR**            columns;            /**< corresponding variables */
   const char**          rows;               /**< corresponding constraint names */ 
   int                   nentries;           /**< number of elements in the arrays */
   int                   sentries;           /**< number of slots in the arrays */
};
typedef struct SparseMatrix SPARSEMATRIX;

/** struct for mapping cons names to numbers */
struct ConsNameFreq
{
   const char*           consname;           /**< name of the constraint */
   int                   freq;               /**< how often we have seen the name */
};
typedef struct ConsNameFreq CONSNAMEFREQ;

/** creates the mps input structure */
static
SCIP_RETCODE mpsinputCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   MPSINPUT**            mpsi,               /**< mps input structure */
   SCIP_FILE*            fp                  /**< file object for the input file */
   )
{
   assert(mpsi != NULL);
   assert(fp != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, mpsi) );

   (*mpsi)->section     = MPS_NAME;
   (*mpsi)->fp          = fp;
   (*mpsi)->lineno      = 0;
   (*mpsi)->objsense    = SCIP_OBJSENSE_MINIMIZE;
   (*mpsi)->haserror    = FALSE;
   (*mpsi)->isinteger   = FALSE;
   (*mpsi)->isnewformat = FALSE;
   (*mpsi)->buf     [0] = '\0';
   (*mpsi)->probname[0] = '\0';
   (*mpsi)->objname [0] = '\0';
   (*mpsi)->f0          = NULL;
   (*mpsi)->f1          = NULL;
   (*mpsi)->f2          = NULL;
   (*mpsi)->f3          = NULL;
   (*mpsi)->f4          = NULL;
   (*mpsi)->f5          = NULL;

   SCIP_CALL( SCIPgetBoolParam(scip, "reading/initialconss", &((*mpsi)->initialconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicconss", &((*mpsi)->dynamicconss)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamiccols", &((*mpsi)->dynamiccols)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "reading/dynamicrows", &((*mpsi)->dynamicrows)) );

   return SCIP_OKAY;
}

/** free the mps input structure */
static
void mpsinputFree(
   SCIP*                 scip,               /**< SCIP data structure */
   MPSINPUT**            mpsi                /**< mps input structure */
   )
{
   SCIPfreeBlockMemory(scip, mpsi);
}

/** returns the current section */
static
MPSSECTION mpsinputSection(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->section;
}

/** return the current value of field 0 */
static
const char* mpsinputField0(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->f0;
}

/** return the current value of field 1 */
static
const char* mpsinputField1(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->f1;
}

/** return the current value of field 2 */
static
const char* mpsinputField2(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->f2;
}

/** return the current value of field 3 */
static
const char* mpsinputField3(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->f3;
}

/** return the current value of field 4 */
static
const char* mpsinputField4(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->f4;
}

/** return the current value of field 5 */
static
const char* mpsinputField5(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->f5;
}

#if 0
/** returns the problem name */
static
const char* mpsinputProbname(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->probname;
}
#endif

/** returns the objective name */
static
const char* mpsinputObjname(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->objname;
}

/** returns the objective sense */
static
SCIP_OBJSENSE mpsinputObjsense(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->objsense;
}

/** returns if an error was detected */
static
SCIP_Bool mpsinputHasError(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->haserror;
}

/** returns the value of the Bool "is integer" in the mps input */
static
SCIP_Bool mpsinputIsInteger(
   const MPSINPUT*       mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   return mpsi->isinteger;
}

/** set the section in the mps input structure to given section */
static
void mpsinputSetSection(
   MPSINPUT*             mpsi,               /**< mps input structure */
   MPSSECTION            section             /**< section that is set */
   )
{
   assert(mpsi != NULL);

   mpsi->section = section;
}

/** set the problem name in the mps input structure to given problem name */
static
void mpsinputSetProbname(
   MPSINPUT*             mpsi,               /**< mps input structure */
   const char*           probname            /**< name of the problem to set */
   )
{
   assert(mpsi     != NULL);
   assert(probname != NULL);
   assert(strlen(probname) < sizeof(mpsi->probname));

   (void)SCIPmemccpy(mpsi->probname, probname, '\0', MPS_MAX_NAMELEN - 1);
}

/** set the objective name in the mps input structure to given objective name */
static
void mpsinputSetObjname(
   MPSINPUT*             mpsi,               /**< mps input structure */
   const char*           objname             /**< name of the objective function to set */
   )
{
   assert(mpsi != NULL);
   assert(objname != NULL);
   assert(strlen(objname) < sizeof(mpsi->objname));

   (void)SCIPmemccpy(mpsi->objname, objname, '\0', MPS_MAX_NAMELEN - 1);
}

/** set the objective sense in the mps input structure to given objective sense */
static
void mpsinputSetObjsense(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP_OBJSENSE         sense               /**< sense of the objective function */
   )
{
   assert(mpsi != NULL);

   mpsi->objsense = sense;
}

static
void mpsinputSyntaxerror(
   MPSINPUT*             mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   SCIPerrorMessage("Syntax error in line %d\n", mpsi->lineno);
   mpsi->section  = MPS_ENDATA;
   mpsi->haserror = TRUE;
}

/** method post a ignore message  */
static
void mpsinputEntryIgnored(
   SCIP*                 scip,               /**< SCIP data structure */
   MPSINPUT*             mpsi,               /**< mps input structure */
   const char*           what,               /**< what get ignored */
   const char*           what_name,          /**< name of that object */
   const char*           entity,             /**< entity */
   const char*           entity_name,        /**< entity name */
   SCIP_VERBLEVEL        verblevel           /**< SCIP verblevel for this message */
   )
{
   assert(mpsi        != NULL);
   assert(what        != NULL);
   assert(what_name   != NULL);
   assert(entity      != NULL);
   assert(entity_name != NULL);

   SCIPverbMessage(scip, verblevel, NULL,
      "Warning line %d: %s \"%s\" for %s \"%s\" ignored\n", mpsi->lineno, what, what_name, entity, entity_name);
}

/** fill the line from \p pos up to column 80 with blanks. */
static
void clearFrom(
   char*                 buf,                /**< buffer to clear */
   unsigned int          pos                 /**< position to start the clearing process */
   )
{
   unsigned int i;

   for(i = pos; i < 80; i++)
      buf[i] = BLANK;
   buf[80] = '\0';
}

/** change all blanks inside a field to #PATCH_CHAR. */
static
void patchField(
   char*                 buf,                /**< buffer to patch */
   int                   beg,                /**< position to begin */
   int                   end                 /**< position to end */
   )
{
   int i;

   while( (beg <= end) && (buf[end] == BLANK) )
      end--;

   while( (beg <= end) && (buf[beg] == BLANK) )
      beg++;

   for( i = beg; i <= end; i++ )
      if( buf[i] == BLANK )
         buf[i] = PATCH_CHAR;
}

/** read a mps format data line and parse the fields. */
static
SCIP_Bool mpsinputReadLine(
   MPSINPUT*             mpsi                /**< mps input structure */
   )
{
   unsigned int len;
   unsigned int i;
   int space;
   char* s;
   SCIP_Bool is_marker;
   SCIP_Bool is_empty;
   char* nexttok;

   do
   {
      mpsi->f0 = mpsi->f1 = mpsi->f2 = mpsi->f3 = mpsi->f4 = mpsi->f5 = 0;
      is_marker = FALSE;

      /* Read until we have not a comment line. */
      do
      {
         mpsi->buf[MPS_MAX_LINELEN-1] = '\0';
         if( NULL == SCIPfgets(mpsi->buf, (int) sizeof(mpsi->buf), mpsi->fp) )
            return FALSE;
         mpsi->lineno++;
      }
      while( *mpsi->buf == '*' );

      /* Normalize line */
      len = (unsigned int) strlen(mpsi->buf);

      for( i = 0; i < len; i++ )
         if( (mpsi->buf[i] == '\t') || (mpsi->buf[i] == '\n') || (mpsi->buf[i] == '\r') )
            mpsi->buf[i] = BLANK;

      if( len < 80 )
         clearFrom(mpsi->buf, len);

      SCIPdebugMessage("line %d: <%s>\n", mpsi->lineno, mpsi->buf);

      assert(strlen(mpsi->buf) >= 80);

      /* Look for new section */
      if( *mpsi->buf != BLANK )
      {
         mpsi->f0 = SCIPstrtok(&mpsi->buf[0], " ", &nexttok);

         assert(mpsi->f0 != 0);

         mpsi->f1 = SCIPstrtok(NULL, " ", &nexttok);

         return TRUE;
      }

      /* If we decide to use the new format we never revert this decision */
      if( !mpsi->isnewformat )
      {
         /* Test for fixed format comments */
         if( (mpsi->buf[14] == '$') && (mpsi->buf[13] == ' ') )
            clearFrom(mpsi->buf, 14);
         else if( (mpsi->buf[39] == '$') && (mpsi->buf[38] == ' ') )
            clearFrom(mpsi->buf, 39);

         /* Test for fixed format */
         space = mpsi->buf[12] | mpsi->buf[13]
            | mpsi->buf[22] | mpsi->buf[23]
            | mpsi->buf[36] | mpsi->buf[37] | mpsi->buf[38]
            | mpsi->buf[47] | mpsi->buf[48]
            | mpsi->buf[61] | mpsi->buf[62] | mpsi->buf[63];

         if( space == BLANK )
         {
            /* Now we have space at the right positions.
             * But are there also the non space where they
             * should be ?
             */
            SCIP_Bool number;

            number = isdigit((unsigned char)mpsi->buf[24]) || isdigit((unsigned char)mpsi->buf[25])
               || isdigit((unsigned char)mpsi->buf[26]) || isdigit((unsigned char)mpsi->buf[27])
               || isdigit((unsigned char)mpsi->buf[28]) || isdigit((unsigned char)mpsi->buf[29])
               || isdigit((unsigned char)mpsi->buf[30]) || isdigit((unsigned char)mpsi->buf[31])
               || isdigit((unsigned char)mpsi->buf[32]) || isdigit((unsigned char)mpsi->buf[33])
               || isdigit((unsigned char)mpsi->buf[34]) || isdigit((unsigned char)mpsi->buf[35]);

            /* len < 14 is handle ROW lines with embedded spaces
             * in the names correctly
             */
            if( number || len < 14 )
            {
               /* We assume fixed format, so we patch possible embedded spaces. */
               patchField(mpsi->buf,  4, 12);
               patchField(mpsi->buf, 14, 22);
               patchField(mpsi->buf, 39, 47);
            }
            else
            {
               if( mpsi->section == MPS_COLUMNS || mpsi->section == MPS_RHS
                  || mpsi->section == MPS_RANGES  || mpsi->section == MPS_BOUNDS )
                  mpsi->isnewformat = TRUE;
            }
         }
         else
         {
            mpsi->isnewformat = TRUE;
         }
      }
      s = &mpsi->buf[1];

      /* At this point it is not clear if we have a indicator field.
       * If there is none (e.g. empty) f1 will be the first name field.
       * If there is one, f2 will be the first name field.
       *
       * Initially comment marks '$' are only allowed in the beginning
       * of the 2nd and 3rd name field. We test all fields but the first.
       * This makes no difference, since if the $ is at the start of a value
       * field, the line will be erroneous anyway.
       */
      do
      {
         if( NULL == (mpsi->f1 = SCIPstrtok(s, " ", &nexttok)) )
            break;

         if( (NULL == (mpsi->f2 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f2 == '$') )
         {
            mpsi->f2 = 0;
            break;
         }
         if( !strcmp(mpsi->f2, "'MARKER'") )
            is_marker = TRUE;

         if( (NULL == (mpsi->f3 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f3 == '$') )
         {
            mpsi->f3 = 0;
            break;
         }
         if( is_marker )
         {
            if( !strcmp(mpsi->f3, "'INTORG'") )
               mpsi->isinteger = TRUE;
            else if( !strcmp(mpsi->f3, "'INTEND'") )
               mpsi->isinteger = FALSE;
            else
               break; /* unknown marker */
         }
         if( !strcmp(mpsi->f3, "'MARKER'") )
            is_marker = TRUE;

         if( (NULL == (mpsi->f4 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f4 == '$') )
         {
            mpsi->f4 = 0;
            break;
         }
         if( is_marker )
         {
            if( !strcmp(mpsi->f4, "'INTORG'") )
               mpsi->isinteger = TRUE;
            else if( !strcmp(mpsi->f4, "'INTEND'") )
               mpsi->isinteger = FALSE;
            else
               break; /* unknown marker */
         }
         if( (NULL == (mpsi->f5 = SCIPstrtok(NULL, " ", &nexttok))) || (*mpsi->f5 == '$') )
            mpsi->f5 = 0;
      }
      while( FALSE );

      /* check for empty lines */
      is_empty = (mpsi->f0 == NULL && mpsi->f1 == NULL);
   }
   while( is_marker || is_empty );

   return TRUE;
}

/** Insert \p str as field 4 and shift all other fields up. */
static
void mpsinputInsertField4(
   MPSINPUT*             mpsi,               /**< mps input structure */
   const char*           str                 /**< str to insert */
   )
{
   assert(mpsi != NULL);
   assert(str != NULL);

   mpsi->f5 = mpsi->f4;
   mpsi->f4 = str;
}

/** Insert \p name as field 1 or 2 and shift all other fields up. */
static
void mpsinputInsertName(
   MPSINPUT*             mpsi,               /**< mps input structure */
   const char*           name,               /**< name to insert */
   SCIP_Bool             second              /**< insert as second field? */
   )
{
   assert(mpsi != NULL);
   assert(name != NULL);

   mpsi->f5 = mpsi->f4;
   mpsi->f4 = mpsi->f3;
   mpsi->f3 = mpsi->f2;

   if( second )
      mpsi->f2 = name;
   else
   {
      mpsi->f2 = mpsi->f1;
      mpsi->f1 = name;
   }
}

/** Process NAME section. */
static
SCIP_RETCODE readName(
   SCIP*                 scip,               /**< SCIP data structure */
   MPSINPUT*             mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   SCIPdebugMsg(scip, "read problem name\n");

   /* This has to be the Line with the NAME section. */
   if( !mpsinputReadLine(mpsi) || mpsinputField0(mpsi) == NULL || strcmp(mpsinputField0(mpsi), "NAME") )
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   /* Sometimes the name is omitted. */
   mpsinputSetProbname(mpsi, (mpsinputField1(mpsi) == 0) ? "_MPS_" : mpsinputField1(mpsi));

   /* This hat to be a new section */
   if( !mpsinputReadLine(mpsi) || (mpsinputField0(mpsi) == NULL) )
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   if( !strncmp(mpsinputField0(mpsi), "ROWS", 4) )
      mpsinputSetSection(mpsi, MPS_ROWS);
   else if( !strncmp(mpsinputField0(mpsi), "USERCUTS", 8) )
      mpsinputSetSection(mpsi, MPS_USERCUTS);
   else if( !strncmp(mpsinputField0(mpsi), "LAZYCONS", 8) )
      mpsinputSetSection(mpsi, MPS_LAZYCONS);
   else if( !strncmp(mpsinputField0(mpsi), "OBJSEN", 6) )
      mpsinputSetSection(mpsi, MPS_OBJSEN);
   else if( !strncmp(mpsinputField0(mpsi), "OBJNAME", 7) )
      mpsinputSetSection(mpsi, MPS_OBJNAME);
   else
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** Process OBJSEN section. This Section is a CPLEX extension. */
static
SCIP_RETCODE readObjsen(
   SCIP*                 scip,               /**< SCIP data structure */
   MPSINPUT*             mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   SCIPdebugMsg(scip, "read objective sense\n");

   /* This has to be the Line with MIN or MAX. */
   if( !mpsinputReadLine(mpsi) || (mpsinputField1(mpsi) == NULL) )
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   if( !strncmp(mpsinputField1(mpsi), "MIN", 3) )
      mpsinputSetObjsense(mpsi, SCIP_OBJSENSE_MINIMIZE);
   else if( !strncmp(mpsinputField1(mpsi), "MAX", 3) )
      mpsinputSetObjsense(mpsi, SCIP_OBJSENSE_MAXIMIZE);
   else
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   /* Look for ROWS, USERCUTS, LAZYCONS, or OBJNAME Section */
   if( !mpsinputReadLine(mpsi) || mpsinputField0(mpsi) == NULL )
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   if( !strcmp(mpsinputField0(mpsi), "ROWS") )
      mpsinputSetSection(mpsi, MPS_ROWS);
   else if( !strcmp(mpsinputField0(mpsi), "USERCUTS") )
      mpsinputSetSection(mpsi, MPS_USERCUTS);
   else if( !strcmp(mpsinputField0(mpsi), "LAZYCONS") )
      mpsinputSetSection(mpsi, MPS_LAZYCONS);
   else if( !strcmp(mpsinputField0(mpsi), "OBJNAME") )
      mpsinputSetSection(mpsi, MPS_OBJNAME);
   else
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** Process OBJNAME section. This Section is a CPLEX extension. */
static
SCIP_RETCODE readObjname(
   SCIP*                 scip,               /**< SCIP data structure */
   MPSINPUT*             mpsi                /**< mps input structure */
   )
{
   assert(mpsi != NULL);

   SCIPdebugMsg(scip, "read objective name\n");

   /* This has to be the Line with the name. */
   if( !mpsinputReadLine(mpsi) || mpsinputField1(mpsi) == NULL )
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   mpsinputSetObjname(mpsi, mpsinputField1(mpsi));

   /* Look for ROWS, USERCUTS, or LAZYCONS Section */
   if( !mpsinputReadLine(mpsi) || mpsinputField0(mpsi) == NULL )
   {
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }
   if( !strcmp(mpsinputField0(mpsi), "ROWS") )
      mpsinputSetSection(mpsi, MPS_ROWS);
   else if( !strcmp(mpsinputField0(mpsi), "USERCUTS") )
      mpsinputSetSection(mpsi, MPS_USERCUTS);
   else if( !strcmp(mpsinputField0(mpsi), "LAZYCONS") )
      mpsinputSetSection(mpsi, MPS_LAZYCONS);
   else
      mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/** Process ROWS, USERCUTS, or LAZYCONS section. */
static
SCIP_RETCODE readRows(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIPdebugMsg(scip, "read rows\n");

   while( mpsinputReadLine(mpsi) )
   {
      if( mpsinputField0(mpsi) != NULL )
      {
         if( !strcmp(mpsinputField0(mpsi), "ROWS") )
            mpsinputSetSection(mpsi, MPS_ROWS);
         else if( !strcmp(mpsinputField0(mpsi), "USERCUTS") )
            mpsinputSetSection(mpsi, MPS_USERCUTS);
         else if( !strcmp(mpsinputField0(mpsi), "LAZYCONS") )
            mpsinputSetSection(mpsi, MPS_LAZYCONS);
         else if( !strcmp(mpsinputField0(mpsi), "COLUMNS") )
            mpsinputSetSection(mpsi, MPS_COLUMNS);
         else
            mpsinputSyntaxerror(mpsi);

         return SCIP_OKAY;
      }

      if( *mpsinputField1(mpsi) == 'N' )
      {
         if( *mpsinputObjname(mpsi) == '\0' )
            mpsinputSetObjname(mpsi, mpsinputField2(mpsi));
         else
            mpsinputEntryIgnored(scip, mpsi, "row", mpsinputField2(mpsi), "objective function", "N", SCIP_VERBLEVEL_NORMAL);
      }
      else
      {
         SCIP_CONS* cons;
         SCIP_Bool initial;
         SCIP_Bool separate;
         SCIP_Bool enforce;
         SCIP_Bool check;
         SCIP_Bool propagate;
         SCIP_Bool local;
         SCIP_Bool modifiable;
         SCIP_Bool dynamic;
         SCIP_Bool removable;

         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons != NULL )
            break;

         initial = mpsi->initialconss && (mpsinputSection(mpsi) == MPS_ROWS);
         separate = TRUE;
         enforce = (mpsinputSection(mpsi) != MPS_USERCUTS);
         check = (mpsinputSection(mpsi) != MPS_USERCUTS);
         propagate = TRUE;
         local = FALSE;
         modifiable = FALSE;
         dynamic = mpsi->dynamicconss;
         removable = mpsi->dynamicrows || (mpsinputSection(mpsi) == MPS_USERCUTS);

         switch(*mpsinputField1(mpsi))
         {
         case 'G' :
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, mpsinputField2(mpsi), 0, NULL, NULL, 0.0, SCIPinfinity(scip),
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
            break;
         case 'E' :
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, mpsinputField2(mpsi), 0, NULL, NULL, 0.0, 0.0,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
            break;
         case 'L' :
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, mpsinputField2(mpsi), 0, NULL, NULL, -SCIPinfinity(scip), 0.0,
                  initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, FALSE) );
            break;
         default :
            mpsinputSyntaxerror(mpsi);
            return SCIP_OKAY;
         }
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/** Process COLUMNS section. */
static
SCIP_RETCODE readCols(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char          colname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS*    cons;
   SCIP_VAR*     var;
   SCIP_Real     val;

   SCIPdebugMsg(scip, "read columns\n");

   var = NULL;
   while( mpsinputReadLine(mpsi) )
   {
      if( mpsinputField0(mpsi) != 0 )
      {
         if( strcmp(mpsinputField0(mpsi), "RHS") )
            break;

         /* add the last variable to the problem */
         if( var != NULL )
         {
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
         assert(var == NULL);

         mpsinputSetSection(mpsi, MPS_RHS);
         return SCIP_OKAY;
      }
      if( mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL )
         break;

      /* new column? */
      if( strcmp(colname, mpsinputField1(mpsi)) )
      {
         /* add the last variable to the problem */
         if( var != NULL )
         {
            SCIP_CALL( SCIPaddVar(scip, var) );
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
         assert(var == NULL);

	 (void)SCIPmemccpy(colname, mpsinputField1(mpsi), '\0', MPS_MAX_NAMELEN - 1);

         if( mpsinputIsInteger(mpsi) )
         {
            /* for integer variables, default bounds are 0 <= x < 1(not +infinity, like it is for continuous variables), and default cost is 0 */
            SCIP_CALL( SCIPcreateVar(scip, &var, colname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
                  !mpsi->dynamiccols, mpsi->dynamiccols, NULL, NULL, NULL, NULL, NULL) );
         }
         else
         {
            /* for continuous variables, default bounds are 0 <= x, and default cost is 0 */
            SCIP_CALL( SCIPcreateVar(scip, &var, colname, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
                  !mpsi->dynamiccols, mpsi->dynamiccols, NULL, NULL, NULL, NULL, NULL) );
         }
      }
      assert(var != NULL);

      val = atof(mpsinputField3(mpsi));

      if( !strcmp(mpsinputField2(mpsi), mpsinputObjname(mpsi)) )
      {
         SCIP_CALL( SCIPchgVarObj(scip, var, val) );
      }
      else
      {
         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons == NULL )
            mpsinputEntryIgnored(scip, mpsi, "Column", mpsinputField1(mpsi), "row", mpsinputField2(mpsi), SCIP_VERBLEVEL_FULL);
         else if( !SCIPisZero(scip, val) )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, val) );
         }
      }
      if( mpsinputField5(mpsi) != NULL )
      {
         assert(mpsinputField4(mpsi) != NULL);

         val = atof(mpsinputField5(mpsi));

         if( !strcmp(mpsinputField4(mpsi), mpsinputObjname(mpsi)) )
         {
            SCIP_CALL( SCIPchgVarObj(scip, var, val) );
         }
         else
         {
            cons = SCIPfindCons(scip, mpsinputField4(mpsi));
            if( cons == NULL )
               mpsinputEntryIgnored(scip, mpsi, "Column", mpsinputField1(mpsi), "row", mpsinputField4(mpsi), SCIP_VERBLEVEL_FULL);
            else if( !SCIPisZero(scip, val) )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, val) );
            }
         }
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/** Process RHS section. */
static
SCIP_RETCODE readRhs(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char        rhsname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS*  cons;
   SCIP_Real   lhs;
   SCIP_Real   rhs;
   SCIP_Real   val;

   SCIPdebugMsg(scip, "read right hand sides\n");

   while( mpsinputReadLine(mpsi) )
   {
      if( mpsinputField0(mpsi) != NULL )
      {
         if( !strcmp(mpsinputField0(mpsi), "RANGES") )
            mpsinputSetSection(mpsi, MPS_RANGES);
         else if( !strcmp(mpsinputField0(mpsi), "BOUNDS") )
            mpsinputSetSection(mpsi, MPS_BOUNDS);
         else if( !strcmp(mpsinputField0(mpsi), "SOS") )
            mpsinputSetSection(mpsi, MPS_SOS);
         else if( !strcmp(mpsinputField0(mpsi), "QMATRIX") )
            mpsinputSetSection(mpsi, MPS_QMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "QUADOBJ") )
            mpsinputSetSection(mpsi, MPS_QUADOBJ);
         else if( !strcmp(mpsinputField0(mpsi), "QCMATRIX") )
            mpsinputSetSection(mpsi, MPS_QCMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "INDICATORS") )
            mpsinputSetSection(mpsi, MPS_INDICATORS);
         else if( !strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
         else
            break;
         return SCIP_OKAY;
      }
      if( (mpsinputField2(mpsi) != NULL && mpsinputField3(mpsi) == NULL)
         || (mpsinputField4(mpsi) != NULL && mpsinputField5(mpsi) == NULL) )
      {
         SCIPwarningMessage(scip, "reading rhs section, a field is missing, assuming that the vector name is the missing one(, row identfier <%s>)\n", mpsinputField2(mpsi));

         mpsinputInsertName(mpsi, "_RHS_", FALSE);
      }

      if( mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL )
         break;

      if( *rhsname == '\0' )
	 (void)SCIPmemccpy(rhsname, mpsinputField1(mpsi), '\0', MPS_MAX_NAMELEN - 1);

      if( !strcmp(rhsname, mpsinputField1(mpsi)) )
      {
         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons == NULL )
         {
            /* the rhs of the objective row is treated as objective constant */
            if( !strcmp(mpsinputField2(mpsi), mpsinputObjname(mpsi)) )
            {
               val = atof(mpsinputField3(mpsi));
               SCIP_CALL( SCIPaddOrigObjoffset(scip, -val) );
            }
            else
               mpsinputEntryIgnored(scip, mpsi, "RHS", mpsinputField1(mpsi), "row", mpsinputField2(mpsi), SCIP_VERBLEVEL_NORMAL);
         }
         else
         {
            val = atof(mpsinputField3(mpsi));

            /* find out the row sense */
            lhs = SCIPgetLhsLinear(scip, cons);
            rhs = SCIPgetRhsLinear(scip, cons);
            if( SCIPisInfinity(scip, -lhs) )
            {
               /* lhs = -infinity -> lower or equal */
               assert(SCIPisZero(scip, rhs));
               SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
            }
            else if( SCIPisInfinity(scip, rhs) )
            {
               /* rhs = +infinity -> greater or equal */
               assert(SCIPisZero(scip, lhs));
               SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
            }
            else
            {
               /* lhs > -infinity, rhs < infinity -> equality */
               assert(SCIPisZero(scip, lhs));
               assert(SCIPisZero(scip, rhs));
               SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
               SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
            }
            SCIPdebugMsg(scip, "RHS <%s> lhs: %g  rhs: %g  val: <%22.12g>\n", mpsinputField2(mpsi), lhs, rhs, val);
         }
         if( mpsinputField5(mpsi) != NULL )
         {
            cons = SCIPfindCons(scip, mpsinputField4(mpsi));
            if( cons == NULL )
            {
               /* the rhs of the objective row is treated as objective constant */
               if( !strcmp(mpsinputField2(mpsi), mpsinputObjname(mpsi)) )
               {
                  val = atof(mpsinputField3(mpsi));
                  SCIP_CALL( SCIPaddOrigObjoffset(scip, -val) );
               }
               else
                  mpsinputEntryIgnored(scip, mpsi, "RHS", mpsinputField1(mpsi), "row", mpsinputField4(mpsi), SCIP_VERBLEVEL_NORMAL);
            }
            else
            {
               val = atof(mpsinputField5(mpsi));

               /* find out the row sense */
               lhs = SCIPgetLhsLinear(scip, cons);
               rhs = SCIPgetRhsLinear(scip, cons);
               if( SCIPisInfinity(scip, -lhs) )
               {
                  /* lhs = -infinity -> lower or equal */
                  assert(SCIPisZero(scip, rhs));
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
               }
               else if( SCIPisInfinity(scip, rhs) )
               {
                  /* rhs = +infinity -> greater or equal */
                  assert(SCIPisZero(scip, lhs));
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
               }
               else
               {
                  /* lhs > -infinity, rhs < infinity -> equality */
                  assert(SCIPisZero(scip, lhs));
                  assert(SCIPisZero(scip, rhs));
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, val) );
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, val) );
               }
               SCIPdebugMsg(scip, "RHS <%s> lhs: %g  rhs: %g  val: <%22.12g>\n", mpsinputField4(mpsi), lhs, rhs, val);
            }
         }
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/** Process RANGES section */
static
SCIP_RETCODE readRanges(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char        rngname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS*  cons;
   SCIP_Real   lhs;
   SCIP_Real   rhs;
   SCIP_Real   val;

   SCIPdebugMsg(scip, "read ranges\n");

   while( mpsinputReadLine(mpsi) )
   {
      if( mpsinputField0(mpsi) != NULL )
      {
         if( !strcmp(mpsinputField0(mpsi), "BOUNDS") )
            mpsinputSetSection(mpsi, MPS_BOUNDS);
         else if( !strcmp(mpsinputField0(mpsi), "SOS") )
            mpsinputSetSection(mpsi, MPS_SOS);
         else if( !strcmp(mpsinputField0(mpsi), "QMATRIX") )
            mpsinputSetSection(mpsi, MPS_QMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "QUADOBJ") )
            mpsinputSetSection(mpsi, MPS_QUADOBJ);
         else if( !strcmp(mpsinputField0(mpsi), "QCMATRIX") )
            mpsinputSetSection(mpsi, MPS_QCMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "INDICATORS") )
            mpsinputSetSection(mpsi, MPS_INDICATORS);
         else if( !strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
         else
            break;
         return SCIP_OKAY;
      }
      if( (mpsinputField2(mpsi) != NULL && mpsinputField3(mpsi) == NULL)
         || (mpsinputField4(mpsi) != NULL && mpsinputField5(mpsi) == NULL) )
      {
         SCIPwarningMessage(scip, "reading ranged section, a field is missing, assuming that the vector name is the missing one(, row identfier <%s>)\n", mpsinputField2(mpsi));

         mpsinputInsertName(mpsi, "_RNG_", FALSE);
      }

      if( mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL )
         break;

      if( *rngname == '\0' )
	 (void)SCIPmemccpy(rngname, mpsinputField1(mpsi), '\0', MPS_MAX_NAMELEN - 1);

      /* The rules are:
       * Row Sign   LHS             RHS
       * ----------------------------------------
       *  G   +/-   rhs             rhs + |range|
       *  L   +/-   rhs - |range|   rhs
       *  E   +     rhs             rhs + range
       *  E   -     rhs + range     rhs
       * ----------------------------------------
       */
      if( !strcmp(rngname, mpsinputField1(mpsi)) )
      {
         cons = SCIPfindCons(scip, mpsinputField2(mpsi));
         if( cons == NULL )
            mpsinputEntryIgnored(scip, mpsi, "Range", mpsinputField1(mpsi), "row", mpsinputField2(mpsi), SCIP_VERBLEVEL_NORMAL);
         else
         {
            val = atof(mpsinputField3(mpsi));

            /* find out the row sense */
            lhs = SCIPgetLhsLinear(scip, cons);
            rhs = SCIPgetRhsLinear(scip, cons);
            if( SCIPisInfinity(scip, -lhs) )
            {
               /* lhs = -infinity -> lower or equal */
               SCIP_CALL( SCIPchgLhsLinear(scip, cons, rhs - REALABS(val)) );
            }
            else if( SCIPisInfinity(scip, rhs) )
            {
               /* rhs = +infinity -> greater or equal */
               SCIP_CALL( SCIPchgRhsLinear(scip, cons, lhs + REALABS(val)) );
            }
            else
            {
               /* lhs > -infinity, rhs < infinity -> equality */
               assert(SCIPisEQ(scip, lhs, rhs));
               if( val >= 0.0 )
               {
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, rhs + val) );
               }
               else
               {
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, lhs + val) );
               }
            }
         }
         if( mpsinputField5(mpsi) != NULL )
         {
            cons = SCIPfindCons(scip, mpsinputField4(mpsi));
            if( cons == NULL )
               mpsinputEntryIgnored(scip, mpsi, "Range", mpsinputField1(mpsi), "row", mpsinputField4(mpsi), SCIP_VERBLEVEL_NORMAL);
            else
            {
               val = atof(mpsinputField5(mpsi));

               /* find out the row sense */
               lhs = SCIPgetLhsLinear(scip, cons);
               rhs = SCIPgetRhsLinear(scip, cons);
               if( SCIPisInfinity(scip, -lhs) )
               {
                  /* lhs = -infinity -> lower or equal */
                  SCIP_CALL( SCIPchgLhsLinear(scip, cons, rhs - REALABS(val)) );
               }
               else if( SCIPisInfinity(scip, rhs) )
               {
                  /* rhs = +infinity -> greater or equal */
                  SCIP_CALL( SCIPchgRhsLinear(scip, cons, lhs + REALABS(val)) );
               }
               else
               {
                  /* lhs > -infinity, rhs < infinity -> equality */
                  assert(SCIPisEQ(scip, lhs, rhs));
                  if( val >= 0.0 )
                  {
                     SCIP_CALL( SCIPchgRhsLinear(scip, cons, rhs + val) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPchgLhsLinear(scip, cons, lhs + val) );
                  }
               }
            }
         }
      }
   }
   mpsinputSyntaxerror(mpsi);

   return SCIP_OKAY;
}

/** Process BOUNDS section. */
static
SCIP_RETCODE readBounds(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char        bndname[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_VAR*   var;
   SCIP_RETCODE retcode;
   SCIP_Real   val;
   SCIP_Bool   shifted;

   SCIP_VAR** semicont;
   int nsemicont;
   int semicontsize;

   retcode = SCIP_OKAY;

   semicont = NULL;
   nsemicont = 0;
   semicontsize = 0;

   SCIPdebugMsg(scip, "read bounds\n");

   while( mpsinputReadLine(mpsi) )
   {
      if( mpsinputField0(mpsi) != 0 )
      {
         if( !strcmp(mpsinputField0(mpsi), "SOS") )
            mpsinputSetSection(mpsi, MPS_SOS);
         else if( !strcmp(mpsinputField0(mpsi), "QMATRIX") )
            mpsinputSetSection(mpsi, MPS_QMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "QUADOBJ") )
            mpsinputSetSection(mpsi, MPS_QUADOBJ);
         else if( !strcmp(mpsinputField0(mpsi), "QCMATRIX") )
            mpsinputSetSection(mpsi, MPS_QCMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "INDICATORS") )
            mpsinputSetSection(mpsi, MPS_INDICATORS);
         else if( !strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
         else
            break;
         goto READBOUNDS_FINISH;
      }

      shifted = FALSE;

      /* Is the value field used ? */
      if( !strcmp(mpsinputField1(mpsi), "LO")  /* lower bound given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "UP")  /* upper bound given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "FX")  /* fixed value given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "LI")  /* CPLEX extension: lower bound of integer variable given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "UI")  /* CPLEX extension: upper bound of integer variable given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "SC")  /* CPLEX extension: semi-continuous variable, upper bound given in field 4 */
         || !strcmp(mpsinputField1(mpsi), "SI") )/* CPLEX extension: semi-integer variable, upper bound given in field 4 */
      {
         if( mpsinputField3(mpsi) != NULL && mpsinputField4(mpsi) == NULL )
         {
            int l;

            /* check what might be missing, if field 3 is a number the bound name might be missing */
            for( l = (int) strlen(mpsinputField3(mpsi)) - 1; l >= 0; --l )
            {
               if( mpsinputField3(mpsi)[l] != '.' && !isdigit(mpsinputField3(mpsi)[l]) )
                  break;
            }

            /* the bound name?! is missing */
            if( l < 0 )
            {
               SCIPwarningMessage(scip, "in bound section a name for value <%s> might be missing\n", mpsinputField3(mpsi));

               mpsinputInsertName(mpsi, "_BND_", TRUE);
               shifted = TRUE;
            }
            /* the bound is be missing */
            else
            {
               SCIPwarningMessage(scip, "in bound section a value for column <%s> is missing, assuming 0.0\n", mpsinputField3(mpsi));

               mpsinputInsertField4(mpsi, "0.0");
               shifted = TRUE;
            }
         }
      }
      else if( !strcmp(mpsinputField1(mpsi), "FR") /* free variable */
         || !strcmp(mpsinputField1(mpsi), "MI")    /* lower bound is minus infinity */
         || !strcmp(mpsinputField1(mpsi), "PL")    /* upper bound is plus infinity */
         || !strcmp(mpsinputField1(mpsi), "BV") )  /* CPLEX extension: binary variable */
      {
         if( mpsinputField2(mpsi) != NULL && mpsinputField3(mpsi) == NULL )
         {
            SCIPwarningMessage(scip, "in bound section a name for a column is missing\n");

            mpsinputInsertName(mpsi, "_BND_", TRUE);
            shifted = TRUE;
         }
      }
      else
      {
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      if( mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL || mpsinputField3(mpsi) == NULL )
         break;

      if( *bndname == '\0' )
	 (void)SCIPmemccpy(bndname, mpsinputField2(mpsi), '\0', MPS_MAX_NAMELEN - 1);

      /* Only read the first Bound in section */
      if( !strcmp(bndname, mpsinputField2(mpsi)) )
      {
         SCIP_VARTYPE oldvartype;
         SCIP_Bool infeasible;

         var = SCIPfindVar(scip, mpsinputField3(mpsi));
         /* if variable did not appear in columns section before, then it may still come in later sections (QCMATRIX, QMATRIX, SOS, ...)
          * thus add it as continuous variables, which has default bounds 0.0 <= x, and default cost 0.0 */
         if( var == NULL )
         {
            SCIP_VAR* varcpy;

            SCIP_CALL( SCIPcreateVar(scip, &var, mpsinputField3(mpsi), 0.0, SCIPinfinity(scip), 0.0, 
                  SCIP_VARTYPE_CONTINUOUS, !mpsi->dynamiccols, mpsi->dynamiccols, NULL, NULL, NULL, NULL, NULL) );

            SCIP_CALL( SCIPaddVar(scip, var) );
            varcpy = var;
            SCIP_CALL( SCIPreleaseVar(scip, &varcpy) );
            /* mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField3(mpsi), "bound", bndname, SCIP_VERBLEVEL_NORMAL); */
         }
         assert(var != NULL);

         if( mpsinputField4(mpsi) == NULL )
            val = 0.0;
         else
            val = atof(mpsinputField4(mpsi));

         /* remember variable type */
         oldvartype = SCIPvarGetType(var);

         /* If a bound of a binary variable is given, the variable is converted into an integer variable
          * with default bounds 0 <= x <= infinity before applying the bound. Note that integer variables
          * are by default assumed to be binary, but an explicit lower bound of 0 turns them into integer variables.
          * Only if the upper bound is explicitly set to 1, we leave the variable as a binary one.
          */
         if( oldvartype == SCIP_VARTYPE_BINARY && !((mpsinputField1(mpsi)[0] == 'U' ||
                  (mpsinputField1(mpsi)[0] == 'F' && mpsinputField1(mpsi)[1] == 'X')) && SCIPisFeasEQ(scip, val, 1.0))
            && !(mpsinputField1(mpsi)[0] == 'F' && mpsinputField1(mpsi)[1] == 'X'&& SCIPisFeasEQ(scip, val, 0.0)) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER, &infeasible) );
            assert(!infeasible);

            oldvartype =  SCIP_VARTYPE_INTEGER;
            SCIP_CALL( SCIPchgVarUb(scip, var, SCIPinfinity(scip)) );
         }

         /* switch variable type to continuous before applying the bound, this is necessary for stupid non-integral
          * bounds on general variables, which even might lead to infeasibility
          */
         if( oldvartype != SCIP_VARTYPE_CONTINUOUS )
         {
            assert(SCIP_VARTYPE_CONTINUOUS >= SCIP_VARTYPE_IMPLINT && SCIP_VARTYPE_IMPLINT >= SCIP_VARTYPE_INTEGER
               && SCIP_VARTYPE_INTEGER >= SCIP_VARTYPE_BINARY);
            /* relaxing variable type */
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_CONTINUOUS, &infeasible) );
         }
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

         switch( mpsinputField1(mpsi)[0] )
         {
         case 'L':
            if( !SCIPisZero(scip, SCIPvarGetLbGlobal(var)) && SCIPisLT(scip, val, SCIPvarGetLbGlobal(var)) )
            {
               SCIPwarningMessage(scip, "Relaxing already defined lower bound %g of variable <%s> to %g not allowed.\n", SCIPvarGetLbGlobal(var), SCIPvarGetName(var), val);
            }

            SCIP_CALL( SCIPchgVarLb(scip, var, val) );

            if( mpsinputField1(mpsi)[1] == 'I' ) /* CPLEX extension (Integer Bound) */
            {
               if( !SCIPisFeasIntegral(scip, val) )
               {
                  SCIPwarningMessage(scip, "variable <%s> declared as integral has a non-integral lower bound (%.14g) -> if feasible, bounds will be adjusted\n", SCIPvarGetName(var), val);
               }
               SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER, &infeasible) );
               /* don't assert feasibility here because the presolver will and should detect a infeasibility */
            }
            else if( oldvartype < SCIP_VARTYPE_CONTINUOUS )
            {
               if( !SCIPisFeasIntegral(scip, val) )
               {
                  SCIPwarningMessage(scip, "variable <%s> declared as integral has a non-integral lower bound (%.14g) -> if feasible, bounds will be adjusted\n", SCIPvarGetName(var), val);
               }
            }

            break;
         case 'U':
            if( SCIPisGT(scip, val, SCIPvarGetUbGlobal(var)) )
            {
               SCIPwarningMessage(scip, "Relaxing already defined upper bound %g of variable <%s> to %g not allowed.\n", SCIPvarGetUbGlobal(var), SCIPvarGetName(var), val);
            }

            SCIP_CALL( SCIPchgVarUb(scip, var, val) );
            if( mpsinputField1(mpsi)[1] == 'I' ) /* CPLEX extension (Integer Bound) */
            {
               if( !SCIPisFeasIntegral(scip, val) )
               {
                  SCIPwarningMessage(scip, "variable <%s> declared as integral has a non-integral upper bound (%.14g) -> if feasible, bounds will be adjusted\n", SCIPvarGetName(var), val);
               }

               SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER, &infeasible) );
               /* don't assert feasibility here because the presolver will and should detect an infeasibility */
            }
            else if( oldvartype < SCIP_VARTYPE_CONTINUOUS )
            {
               if( !SCIPisFeasIntegral(scip, val) )
               {
                  SCIPwarningMessage(scip, "variable <%s> declared as integral has a non-integral upper bound (%.14g) -> if feasible, bounds will be adjusted\n", SCIPvarGetName(var), val);
               }
            }
            break;
         case 'S':
            assert(mpsinputField1(mpsi)[1] == 'C' || mpsinputField1(mpsi)[1] == 'I'); /* semi-continuous or semi-integer (CPLEX extension) */
            /* remember that variable is semi-continuous/-integer */
            if( semicontsize <= nsemicont )
            {
               semicontsize = SCIPcalcMemGrowSize(scip, nsemicont+1);
               if( semicont == NULL )
               {
                  SCIP_CALL( SCIPallocBufferArray(scip, &semicont, semicontsize) );
               }
               else
               {
                  SCIP_CALL( SCIPreallocBufferArray(scip, &semicont, semicontsize) );
               }
            }
            assert(semicont != NULL);
            semicont[nsemicont] = var;
            ++nsemicont;

            if( mpsinputField1(mpsi)[1] == 'I' ) /* variable is semi-integer, hence change its type to integer (the "semi" part will be handled below) */
            {
               SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_INTEGER, &infeasible) );
               /* don't assert feasibility here because the presolver will and should detect an infeasibility */
            }

            /* if both bounds are infinite anyway, we do not need to print a warning or change the bound */
            if( !SCIPisInfinity(scip, val) || !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
            {
               if( SCIPisGT(scip, val, SCIPvarGetUbGlobal(var)) )
               {
                  SCIPwarningMessage(scip, "Relaxing already defined upper bound %g of variable <%s> to %g not allowed.\n", SCIPvarGetUbGlobal(var), SCIPvarGetName(var), val);
               }

               SCIP_CALL( SCIPchgVarUb(scip, var, val) );
            }
            break;
         case 'F':
            if( mpsinputField1(mpsi)[1] == 'X' )
            {
               SCIP_CALL( SCIPchgVarLb(scip, var, val) );
               SCIP_CALL( SCIPchgVarUb(scip, var, val) );
            }
            else
            {
               SCIP_CALL( SCIPchgVarLb(scip, var, -SCIPinfinity(scip)) );
               SCIP_CALL( SCIPchgVarUb(scip, var, +SCIPinfinity(scip)) );
            }
            break;
         case 'M':
            SCIP_CALL( SCIPchgVarLb(scip, var, -SCIPinfinity(scip)) );
            break;
         case 'P':
            SCIP_CALL( SCIPchgVarUb(scip, var, +SCIPinfinity(scip)) );
            break;
         case 'B' : /* CPLEX extension (Binary) */
            SCIP_CALL( SCIPchgVarLb(scip, var, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, var, 1.0) );
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            /* don't assert feasibility here because the presolver will and should detect a infeasibility */
            break;
         default:
            mpsinputSyntaxerror(mpsi);
            return SCIP_OKAY;
         }

         /* switch variable type back to old type if necessary */
         if( oldvartype < SCIPvarGetType(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, oldvartype, &infeasible) );
         }
      }
      else
      {
         /* check for syntax error */
         assert(*bndname != '\0');
         if( strcmp(bndname, mpsinputField3(mpsi)) == 0 && shifted )
         {
            mpsinputSyntaxerror(mpsi);
            return SCIP_OKAY;
         }

         mpsinputEntryIgnored(scip, mpsi, "bound", mpsinputField2(mpsi), "variable", mpsinputField3(mpsi), SCIP_VERBLEVEL_NORMAL);
      }
   }
   mpsinputSyntaxerror(mpsi);


 READBOUNDS_FINISH:
   if( nsemicont > 0 )
   {
      SCIP_CONS* cons;
      SCIP_VAR* vars[2];
      SCIP_BOUNDTYPE boundtypes[2];
      SCIP_Real bounds[2];
      char name[SCIP_MAXSTRLEN];
      SCIP_Real oldlb;
      int i;

      assert(semicont != NULL);

      /* add bound disjunction constraints for semi-continuous and semi-integer variables */
      for( i = 0; i < nsemicont; ++i )
      {
         var = semicont[i];
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

         oldlb = SCIPvarGetLbGlobal(var);
         assert(oldlb >= 0.0);

         /* if no bound was specified (which we assume if we see lower bound 0.0),
          * then the default lower bound for a semi-continuous variable is 1.0 */
         if( oldlb == 0.0 )
            oldlb = 1.0;

         /* change the lower bound to 0.0 */
         SCIP_CALL( SCIPchgVarLb(scip, var, 0.0) );

         /* add a bound disjunction constraint to say var <= 0.0 or var >= oldlb */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "semicont_%s", SCIPvarGetName(var));

         vars[0] = var;
         vars[1] = var;
         boundtypes[0] = SCIP_BOUNDTYPE_UPPER;
         boundtypes[1] = SCIP_BOUNDTYPE_LOWER;
         bounds[0] = 0.0;
         bounds[1] = oldlb;

         retcode = SCIPcreateConsBounddisjunction(scip, &cons, name, 2, vars, boundtypes, bounds,
            !mpsi->dynamiccols, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, mpsi->dynamicconss, mpsi->dynamiccols, FALSE);

         if( retcode != SCIP_OKAY )
            break;

         SCIP_CALL( SCIPaddCons(scip, cons) );

         SCIPdebugMsg(scip, "add bound disjunction constraint for semi-continuity/-integrality of <%s>:\n\t", SCIPvarGetName(var));
         SCIPdebugPrintCons(scip, cons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   SCIPfreeBufferArrayNull(scip, &semicont);

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}


/** Process SOS section.
 *
 *  We read the SOS section, which is a nonstandard section introduced by CPLEX.
 *
 *  @note Currently we do not support the standard way of specifying SOS constraints via markers.
 */
static
SCIP_RETCODE readSOS(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   char name[MPS_MAX_NAMELEN] = { '\0' };
   SCIP_CONS* cons = NULL;
   int consType = -1;
   int cnt = 0;

   SCIPdebugMsg(scip, "read SOS constraints\n");

   /* standard settings for SOS constraints: */
   initial = mpsi->initialconss;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   dynamic = mpsi->dynamicconss;
   removable = mpsi->dynamicrows;

   /* loop through section */
   while( mpsinputReadLine(mpsi) )
   {
      int type = -1;

      /* check if next section is found */
      if( mpsinputField0(mpsi) != NULL )
      {
         if( !strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
         else if( !strcmp(mpsinputField0(mpsi), "QMATRIX") )
            mpsinputSetSection(mpsi, MPS_QMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "QUADOBJ") )
            mpsinputSetSection(mpsi, MPS_QUADOBJ);
         else if( !strcmp(mpsinputField0(mpsi), "QCMATRIX") )
            mpsinputSetSection(mpsi, MPS_QCMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "INDICATORS") )
            mpsinputSetSection(mpsi, MPS_INDICATORS);
         break;
      }
      if( mpsinputField1(mpsi) == NULL )
      {
         SCIPerrorMessage("empty data in a non-comment line.\n");
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      /* check for new SOS set */
      if( strcmp(mpsinputField1(mpsi), "S1") == 0 )
         type = 1;
      if( strcmp(mpsinputField1(mpsi), "S2") == 0 )
         type = 2;

      /* add last constraint and create a new one */
      if( type > 0 )
      {
         assert( type == 1 || type == 2 );
         if( cons != NULL )
         {
            /* add last constraint */
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIPdebugMsg(scip, "(line %d) added constraint <%s>: ", mpsi->lineno, SCIPconsGetName(cons));
            SCIPdebugPrintCons(scip, cons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }

         /* check name */
         if( mpsinputField2(mpsi) != NULL )
            (void)SCIPmemccpy(name, mpsinputField2(mpsi), '\0', MPS_MAX_NAMELEN - 1);
         else
         {
            /* create new name */
            (void) SCIPsnprintf(name, MPS_MAX_NAMELEN, "SOS%d", ++cnt);
         }

         /* create new SOS constraint */
         if( type == 1 )
         {
            /* we do not know the name of the constraint */
            SCIP_CALL( SCIPcreateConsSOS1(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
                  local, dynamic, removable, FALSE) );
         }
         else
         {
            assert( type == 2 );
            SCIP_CALL( SCIPcreateConsSOS2(scip, &cons, name, 0, NULL, NULL, initial, separate, enforce, check, propagate,
                  local, dynamic, removable, FALSE) );
         }
         consType = type;
         SCIPdebugMsg(scip, "created constraint <%s> of type %d.\n", name, type);
         /* note: we ignore the priorities! */
      }
      else
      {
         /* otherwise we are in the section given variables */
         SCIP_VAR* var;
         SCIP_Real weight;
         char* endptr;

         if( consType != 1 && consType != 2 )
         {
            SCIPerrorMessage("missing SOS type specification.\n");
            mpsinputSyntaxerror(mpsi);
            return SCIP_OKAY;
         }

         /* get variable */
         var = SCIPfindVar(scip, mpsinputField1(mpsi));
         if( var == NULL )
         {
            /* ignore unknown variables - we would not know the type anyway */
            mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField1(mpsi), "SOS", name, SCIP_VERBLEVEL_NORMAL);
         }
         else
         {
            /* get weight */
            weight = strtod(mpsinputField2(mpsi), &endptr);
            if( endptr == mpsinputField2(mpsi) || *endptr != '\0' )
            {
               SCIPerrorMessage("weight for variable <%s> not specified.\n", mpsinputField1(mpsi));
               mpsinputSyntaxerror(mpsi);
               return SCIP_OKAY;
            }

            /* add variable and weight */
            assert( consType == 1 || consType == 2 );
            switch( consType )
            {
            case 1: 
               SCIP_CALL( SCIPaddVarSOS1(scip, cons, var, weight) );
               break;
            case 2: 
               SCIP_CALL( SCIPaddVarSOS2(scip, cons, var, weight) );
               break;
            default: 
               SCIPerrorMessage("unknown SOS type: <%d>\n", type); /* should not happen */
               SCIPABORT();
               return SCIP_INVALIDDATA;  /*lint !e527*/
            }
            SCIPdebugMsg(scip, "added variable <%s> with weight %g.\n", SCIPvarGetName(var), weight);
         }
         /* check other fields */
         if( (mpsinputField3(mpsi) != NULL && *mpsinputField3(mpsi) != '\0' ) ||
            (mpsinputField4(mpsi) != NULL && *mpsinputField4(mpsi) != '\0' ) ||
            (mpsinputField5(mpsi) != NULL && *mpsinputField5(mpsi) != '\0' ) )
         {
            SCIPwarningMessage(scip, "ignoring data in fields 3-5 <%s> <%s> <%s>.\n",
               mpsinputField3(mpsi), mpsinputField4(mpsi), mpsinputField5(mpsi));
         }
      }
   }

   if( cons != NULL )
   {
      /* add last constraint */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMsg(scip, "(line %d) added constraint <%s>: ", mpsi->lineno, SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}


/** Process QMATRIX or QUADOBJ section.
 *
 *  - We read the QMATRIX or QUADOBJ section, which is a nonstandard section introduced by CPLEX.
 *  - We create a quadratic constraint for this matrix and add a variable to the objective to
 *    represent the value of the QMATRIX.
 *  - For a QMATRIX, we expect that both lower and upper diagonal elements are given and every
 *    coefficient has to be divided by 2.0.
 *  - For a QUADOBJ, we expect that only the upper diagonal elements are given and thus only
 *    coefficients on the diagonal have to be divided by 2.0.
 */
static
SCIP_RETCODE readQMatrix(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP_Bool             isQuadObj,          /**< whether we actually read a QUADOBJ section */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;
   SCIP_RETCODE retcode;
   int cnt  = 0; /* number of qmatrix elements processed so far */
   int size;     /* size of quad* arrays */

   SCIPdebugMsg(scip, "read %s objective\n", isQuadObj ? "QUADOBJ" : "QMATRIX");

   retcode = SCIP_OKAY;

   size = 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, size) );

   /* loop through section */
   while( mpsinputReadLine(mpsi) )
   {
      /* otherwise we are in the section given variables */
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real coef;

      /* check if next section is found */
      if( mpsinputField0(mpsi) != NULL )
      {
         if( !strcmp(mpsinputField0(mpsi), "QCMATRIX") )
            mpsinputSetSection(mpsi, MPS_QCMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "INDICATORS") )
            mpsinputSetSection(mpsi, MPS_INDICATORS);
         else if( !strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
         break;
      }
      if( mpsinputField1(mpsi) == NULL && mpsinputField2(mpsi) == NULL )
      {
         SCIPerrorMessage("empty data in a non-comment line.\n");
         mpsinputSyntaxerror(mpsi);
         SCIPfreeBufferArray(scip, &quadvars1);
         SCIPfreeBufferArray(scip, &quadvars2);
         SCIPfreeBufferArray(scip, &quadcoefs);
         return SCIP_OKAY;
      }

      /* get first variable */
      var1 = SCIPfindVar(scip, mpsinputField1(mpsi));
      if( var1 == NULL )
      {
         /* ignore unknown variables - we would not know the type anyway */
         mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField1(mpsi), "QMatrix", "QMATRIX", SCIP_VERBLEVEL_NORMAL);
      }
      else
      {
         int k;
         for( k = 1; k <= 2; ++k )
         {
            /* get second variable */
            var2 = SCIPfindVar(scip, k == 1 ? mpsinputField2(mpsi) : mpsinputField4(mpsi));
            if( var2 == NULL )
            {
               /* ignore unknown variables - we would not know the type anyway */
               mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField2(mpsi), "QMatrix", "QMATRIX", SCIP_VERBLEVEL_NORMAL);
            }
            else
            {
               const char* field;
               char* endptr;

               /* get coefficient */
               field = (k == 1 ? mpsinputField3(mpsi) :  mpsinputField5(mpsi));
               coef = strtod(field, &endptr);
               if( endptr == field || *endptr != '\0' )
               {
                  SCIPerrorMessage("coefficient of term <%s>*<%s> not specified.\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
                  mpsinputSyntaxerror(mpsi);
                  SCIPfreeBufferArray(scip, &quadvars1);
                  SCIPfreeBufferArray(scip, &quadvars2);
                  SCIPfreeBufferArray(scip, &quadcoefs);
                  return SCIP_OKAY;
               }

               /* store variables and coefficient */
               if( cnt >= size )
               {
                  int newsize = SCIPcalcMemGrowSize(scip, size+1);
                  assert(newsize > size);
                  SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars1, newsize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars2, newsize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &quadcoefs, newsize) );
                  size = newsize;
               }
               assert(cnt < size);
               quadvars1[cnt] = var1;
               quadvars2[cnt] = var2;
               quadcoefs[cnt] = coef;

               /* diagonal elements have to be divided by 2.0
                * in a QMATRIX section also off-diagonal have to be divided by 2.0, since both lower and upper diagonal elements are given
                */
               if( var1 == var2 || !isQuadObj )
                  quadcoefs[cnt] /= 2.0;
               ++cnt;

               SCIPdebugMsg(scip, "stored term %g*<%s>*<%s>.\n", coef, SCIPvarGetName(var1), SCIPvarGetName(var2));
            }

            if( mpsinputField4(mpsi) == NULL || *mpsinputField4(mpsi) == '\0' )
               break;

            if( mpsinputField5(mpsi) == NULL || *mpsinputField5(mpsi) == '\0' )
            {
               /* ignore unknown variables - we would not know the type anyway */
               mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField4(mpsi), "QMatrix", "QMATRIX", SCIP_VERBLEVEL_NORMAL);
               break;
            }
         }
      }
   }

   /* add constraint */
   if( cnt )
   {
      SCIP_Bool  initial, separate, enforce, check, propagate;
      SCIP_Bool  local, modifiable, dynamic, removable;
      SCIP_CONS* cons = NULL;
      SCIP_VAR*  qmatrixvar = NULL;
      SCIP_Real  lhs, rhs;
      SCIP_Real  minusone = -1.0;

      /* determine settings; note that reading/{initialconss,dynamicconss,dynamicrows,dynamiccols} apply only to model
       * constraints and variables, not to an auxiliary objective constraint (otherwise it can happen that an auxiliary
       * objective variable is loose with infinite best bound, triggering the problem that an LP that is unbounded
       * because of loose variables with infinite best bound cannot be solved)
       */
      initial    = TRUE;
      separate   = TRUE;
      enforce    = TRUE;
      check      = TRUE;
      propagate  = TRUE;
      local      = FALSE;
      modifiable = FALSE;
      dynamic    = FALSE;
      removable  = FALSE;

      SCIP_CALL( SCIPcreateVar(scip, &qmatrixvar, "qmatrixvar", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, qmatrixvar) );

      if( mpsinputObjsense(mpsi) == SCIP_OBJSENSE_MINIMIZE )
      {
         lhs = -SCIPinfinity(scip);
         rhs = 0.0;
      }
      else
      {
         lhs = 0.0;
         rhs = SCIPinfinity(scip);
      }

      retcode = SCIPcreateConsQuadratic(scip, &cons, "qmatrix", 1, &qmatrixvar, &minusone, cnt, quadvars1, quadvars2, quadcoefs, lhs, rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable);

      if( retcode == SCIP_OKAY )
      {
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIPdebugMsg(scip, "(line %d) added constraint <%s>: ", mpsi->lineno, SCIPconsGetName(cons));
         SCIPdebugPrintCons(scip, cons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         SCIP_CALL( SCIPreleaseVar(scip, &qmatrixvar) );
      }
   }
   else
   {
      SCIPwarningMessage(scip, "%s section has no entries.\n", isQuadObj ? "QUADOBJ" : "QMATRIX");
   }

   SCIPfreeBufferArray(scip, &quadvars1);
   SCIPfreeBufferArray(scip, &quadvars2);
   SCIPfreeBufferArray(scip, &quadcoefs);

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}


/** Process QCMATRIX section.
 *
 *  We read the QCMATRIX section, which is a nonstandard section introduced by CPLEX.
 *
 *  We replace the corresponding linear constraint by a quadratic constraint which contains the
 *  original linear constraint plus the quadratic part specified in the QCMATRIX.
 */
static
SCIP_RETCODE readQCMatrix(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONS* lincons; /* the linear constraint that was added for the corresponding row */
   SCIP_VAR** quadvars1;
   SCIP_VAR** quadvars2;
   SCIP_Real* quadcoefs;
   SCIP_RETCODE retcode;
   int cnt  = 0; /* number of qcmatrix elements processed so far */
   int size;     /* size of quad* arrays */

   if( mpsinputField1(mpsi) == NULL )
   {
      SCIPerrorMessage("no row name in QCMATRIX line.\n");
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   retcode = SCIP_OKAY;

   SCIPdebugMsg(scip, "read QCMATRIX section for row <%s>\n", mpsinputField1(mpsi));

   lincons = SCIPfindCons(scip, mpsinputField1(mpsi));
   if( lincons == NULL )
   {
      SCIPerrorMessage("no row under name <%s> processed so far.\n");
      mpsinputSyntaxerror(mpsi);
      return SCIP_OKAY;
   }

   size = 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, size) );

   /* loop through section */
   while( mpsinputReadLine(mpsi) )
   {
      /* otherwise we are in the section given variables */
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real coef;

      /* check if next section is found */
      if( mpsinputField0(mpsi) != NULL )
      {
         if( !strcmp(mpsinputField0(mpsi), "QMATRIX") )
            mpsinputSetSection(mpsi, MPS_QMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "QUADOBJ") )
            mpsinputSetSection(mpsi, MPS_QUADOBJ);
         else if( !strcmp(mpsinputField0(mpsi), "QCMATRIX") )
            mpsinputSetSection(mpsi, MPS_QCMATRIX);
         else if( !strcmp(mpsinputField0(mpsi), "INDICATORS") )
            mpsinputSetSection(mpsi, MPS_INDICATORS);
         else if( !strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
         break;
      }
      if( mpsinputField1(mpsi) == NULL && mpsinputField2(mpsi) == NULL )
      {
         SCIPerrorMessage("empty data in a non-comment line.\n");
         mpsinputSyntaxerror(mpsi);

         goto TERMINATE;
      }

      /* get first variable */
      var1 = SCIPfindVar(scip, mpsinputField1(mpsi));
      if( var1 == NULL )
      {
         /* ignore unknown variables - we would not know the type anyway */
         mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField1(mpsi), "QCMatrix", SCIPconsGetName(lincons), SCIP_VERBLEVEL_NORMAL);
      }
      else
      {
         /* get second variable */
         var2 = SCIPfindVar(scip, mpsinputField2(mpsi));
         if( var2 == NULL )
         {
            /* ignore unknown variables - we would not know the type anyway */
            mpsinputEntryIgnored(scip, mpsi, "column", mpsinputField2(mpsi), "QCMatrix", SCIPconsGetName(lincons), SCIP_VERBLEVEL_NORMAL);
         }
         else
         {
            char* endptr;
            /* get coefficient */
            coef = strtod(mpsinputField3(mpsi), &endptr);
            if( endptr == mpsinputField3(mpsi) || *endptr != '\0' )
            {
               SCIPerrorMessage("coefficient of term <%s>*<%s> not specified.\n", mpsinputField1(mpsi), mpsinputField2(mpsi));
               mpsinputSyntaxerror(mpsi);

               goto TERMINATE;
            }

            /* store variables and coefficient */
            if( cnt >= size )
            {
               int newsize = SCIPcalcMemGrowSize(scip, size+1);
               assert(newsize > size);
               SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars1, newsize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &quadvars2, newsize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, &quadcoefs, newsize) );
               size = newsize;
            }
            assert(cnt < size);
            quadvars1[cnt] = var1;
            quadvars2[cnt] = var2;
            quadcoefs[cnt] = coef;
            ++cnt;

            SCIPdebugMsg(scip, "stored term %g*<%s>*<%s>.\n", coef, SCIPvarGetName(var1), SCIPvarGetName(var2));

            /* check other fields */
            if( (mpsinputField4(mpsi) != NULL && *mpsinputField4(mpsi) != '\0' ) ||
               (mpsinputField5(mpsi) != NULL && *mpsinputField5(mpsi) != '\0' ) )
            {
               SCIPwarningMessage(scip, "ignoring data in fields 4 and 5 <%s> <%s>.\n", mpsinputField4(mpsi), mpsinputField5(mpsi));
            }
         }
      }
   }

   /* replace linear constraint by quadratic constraint */
   if( cnt )
   {
      SCIP_CONS* cons = NULL;

      retcode = SCIPcreateConsQuadratic(scip, &cons, SCIPconsGetName(lincons),
            SCIPgetNVarsLinear(scip, lincons), SCIPgetVarsLinear(scip, lincons), SCIPgetValsLinear(scip, lincons),
            cnt, quadvars1, quadvars2, quadcoefs, SCIPgetLhsLinear(scip, lincons), SCIPgetRhsLinear(scip, lincons),
            SCIPconsIsInitial(lincons), SCIPconsIsSeparated(lincons), SCIPconsIsEnforced(lincons), SCIPconsIsChecked(lincons),
            SCIPconsIsPropagated(lincons), SCIPconsIsLocal(lincons), SCIPconsIsModifiable(lincons), SCIPconsIsDynamic(lincons),
            SCIPconsIsRemovable(lincons));

      if( retcode != SCIP_OKAY )
         goto TERMINATE;

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMsg(scip, "(line %d) added constraint <%s>: ", mpsi->lineno, SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      SCIP_CALL( SCIPdelCons(scip, lincons) );
   }
   else
   {
      SCIPwarningMessage(scip, "QCMATRIX section has no entries.\n");
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &quadcoefs);
   SCIPfreeBufferArray(scip, &quadvars2);
   SCIPfreeBufferArray(scip, &quadvars1);

   SCIP_CALL( retcode );

   return SCIP_OKAY;
}


/** Process INDICATORS section.
 *
 *  We read the INDICATORS section, which is a nonstandard section introduced by CPLEX.
 *  Note that CPLEX does not allow ranged rows.
 *
 *  If the linear constraints are equations or ranged rows, we generate two indicator
 *  constraints.
 *
 *  The section has to come after the QMATRIX* sections.
 */
static
SCIP_RETCODE readIndicators(
   MPSINPUT*             mpsi,               /**< mps input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool initial;
   SCIP_Bool separate;
   SCIP_Bool enforce;
   SCIP_Bool check;
   SCIP_Bool propagate;
   SCIP_Bool local;
   SCIP_Bool dynamic;
   SCIP_Bool removable;
   SCIP_Bool stickingatnode;
   char name[MPS_MAX_NAMELEN] = { '\0' };

   SCIPdebugMsg(scip, "read INDICATORS constraints\n");

   /* standard settings for indicator constraints: */
   initial = mpsi->initialconss;
   separate = TRUE;
   enforce = TRUE;
   check = TRUE;
   propagate = TRUE;
   local = FALSE;
   dynamic = mpsi->dynamicconss;
   removable = mpsi->dynamicrows;
   stickingatnode = FALSE;

   /* loop through section */
   while( mpsinputReadLine(mpsi) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_VARTYPE slackvartype;
      SCIP_CONS* cons;
      SCIP_CONS* lincons;
      SCIP_VAR* binvar;
      SCIP_VAR* slackvar;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real sign;
      SCIP_VAR** linvars;
      SCIP_Real* linvals;
      int nlinvars;
      int i;

      /* check if next section is found */
      if( mpsinputField0(mpsi) != NULL )
      {
         if( !strcmp(mpsinputField0(mpsi), "ENDATA") )
            mpsinputSetSection(mpsi, MPS_ENDATA);
         break;
      }
      if( mpsinputField1(mpsi) == NULL || mpsinputField2(mpsi) == NULL )
      {
         SCIPerrorMessage("empty data in a non-comment line.\n");
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      /* check for new indicator constraint */
      if( strcmp(mpsinputField1(mpsi), "IF") != 0 )
      {
         SCIPerrorMessage("Indicator constraints need to be introduced by 'IF' in column 1.\n");
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      /* get linear constraint (row) */
      lincons = SCIPfindCons(scip, mpsinputField2(mpsi));
      if( lincons == NULL )
      {
         SCIPerrorMessage("row <%s> does not exist.\n", mpsinputField2(mpsi));
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      /* check whether constraint is really linear */
      conshdlr = SCIPconsGetHdlr(lincons);
      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear") != 0 )
      {
         SCIPerrorMessage("constraint <%s> is not linear.\n", mpsinputField2(mpsi));
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      /* get binary variable */
      binvar = SCIPfindVar(scip, mpsinputField3(mpsi));
      if( binvar == NULL )
      {
         SCIPerrorMessage("binary variable <%s> does not exist.\n", mpsinputField3(mpsi));
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      /* check type */
      if( SCIPvarGetType(binvar) != SCIP_VARTYPE_BINARY )
      {
         SCIPerrorMessage("variable <%s> is not binary.\n", mpsinputField3(mpsi));
         mpsinputSyntaxerror(mpsi);
         return SCIP_OKAY;
      }

      /* check whether we need the negated variable */
      if( mpsinputField4(mpsi) != NULL )
      {
         if( *mpsinputField4(mpsi) == '0' )
         {
            SCIP_VAR* var;
            SCIP_CALL( SCIPgetNegatedVar(scip, binvar, &var) );
            binvar = var;
            assert( binvar != NULL );
         }
         else
         {
            if( *mpsinputField4(mpsi) != '1' )
            {
               SCIPerrorMessage("binary variable <%s> can only take values 0/1 (%s).\n", mpsinputField3(mpsi), mpsinputField4(mpsi));
               mpsinputSyntaxerror(mpsi);
               return SCIP_OKAY;
            }
         }
      }

      /* check lhs/rhs */
      lhs = SCIPgetLhsLinear(scip, lincons);
      rhs = SCIPgetRhsLinear(scip, lincons);
      nlinvars = SCIPgetNVarsLinear(scip, lincons);
      linvars = SCIPgetVarsLinear(scip, lincons);
      linvals = SCIPgetValsLinear(scip, lincons);

      sign = -1.0;
      if( !SCIPisInfinity(scip, -lhs) )
      {
         if( SCIPisInfinity(scip, rhs) )
            sign = 1.0;
         else
         {
            /* create second indicator constraint */
            SCIP_VAR** vars;
            SCIP_Real* vals;
            SCIP_RETCODE retcode;

            SCIP_CALL( SCIPallocBufferArray(scip, &vars, nlinvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vals, nlinvars) );
            for( i = 0; i < nlinvars; ++i )
            {
               vars[i] = linvars[i];
               vals[i] = -linvals[i];
            }

            /* create new name */
            (void) SCIPsnprintf(name, MPS_MAX_NAMELEN, "indlhs_%s", SCIPconsGetName(lincons));

            /* create indicator constraint */
            retcode = SCIPcreateConsIndicator(scip, &cons, name, binvar, nlinvars, vars, vals, -lhs,
               initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode);

            if( retcode == SCIP_OKAY )
            {
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIPdebugMsg(scip, "created indicator constraint <%s>\n", mpsinputField2(mpsi));
               SCIPdebugPrintCons(scip, cons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }

            SCIPfreeBufferArray(scip, &vals);
            SCIPfreeBufferArray(scip, &vars);

            SCIP_CALL( retcode );
         }
      }

      /* check if slack variable can be made implicitly integer */
      slackvartype = SCIP_VARTYPE_IMPLINT;
      for (i = 0; i < nlinvars; ++i)
      {
         if( ! SCIPvarIsIntegral(linvars[i]) || ! SCIPisIntegral(scip, linvals[i]) )
         {
            slackvartype = SCIP_VARTYPE_CONTINUOUS;
            break;
         }
      }

      /* create slack variable */
      if ( ! SCIPisInfinity(scip, -lhs) )
         (void) SCIPsnprintf(name, MPS_MAX_NAMELEN, "indslack_indrhs_%s", SCIPconsGetName(lincons));
      else
         (void) SCIPsnprintf(name, MPS_MAX_NAMELEN, "indslack_%s", SCIPconsGetName(lincons));
      SCIP_CALL( SCIPcreateVar(scip, &slackvar, name, 0.0, SCIPinfinity(scip), 0.0, slackvartype, TRUE, FALSE,
            NULL, NULL, NULL, NULL, NULL) );

      /* add slack variable */
      SCIP_CALL( SCIPaddVar(scip, slackvar) );
      SCIP_CALL( SCIPaddCoefLinear(scip, lincons, slackvar, sign) );

      /* correct linear constraint and create new name */
      if ( ! SCIPisInfinity(scip, -lhs) )
      {
         /* we have added lhs above and only need the rhs */
         SCIP_CALL( SCIPchgLhsLinear(scip, lincons, -SCIPinfinity(scip) ) );
         (void) SCIPsnprintf(name, MPS_MAX_NAMELEN, "indrhs_%s", SCIPconsGetName(lincons));
      }
      else
         (void) SCIPsnprintf(name, MPS_MAX_NAMELEN, "ind_%s", SCIPconsGetName(lincons));

      /* create indicator constraint */
      SCIP_CALL( SCIPcreateConsIndicatorLinCons(scip, &cons, name, binvar, lincons, slackvar,
            initial, separate, enforce, check, propagate, local, dynamic, removable, stickingatnode) );

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebugMsg(scip, "created indicator constraint <%s>", mpsinputField2(mpsi));
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}


/** Read LP in "MPS File Format".
 *
 *  A specification of the MPS format can be found at
 *
 *  http://plato.asu.edu/ftp/mps_format.txt,
 *  ftp://ftp.caam.rice.edu/pub/people/bixby/miplib/mps_format,
 *
 *  and in the
 *
 *  CPLEX Reference Manual
 *
 *  This routine should read all valid MPS format files.
 *  What it will not do, is to find all cases where a file is ill formed.
 *  If this happens it may complain and read nothing or read "something".
 */
static
SCIP_RETCODE readMps(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* fp;
   MPSINPUT* mpsi;
   SCIP_RETCODE retcode;
   SCIP_Bool error = TRUE;

   assert(scip != NULL);
   assert(filename != NULL);

   fp = SCIPfopen(filename, "r");
   if( fp == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   SCIP_CALL( mpsinputCreate(scip, &mpsi, fp) );

   SCIP_CALL_TERMINATE( retcode, readName(scip, mpsi), TERMINATE );

   SCIP_CALL_TERMINATE( retcode, SCIPcreateProb(scip, mpsi->probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL), TERMINATE );

   if( mpsinputSection(mpsi) == MPS_OBJSEN )
   {
      SCIP_CALL_TERMINATE( retcode, readObjsen(scip, mpsi), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_OBJNAME )
   {
      SCIP_CALL_TERMINATE( retcode, readObjname(scip, mpsi), TERMINATE );
   }
   while( mpsinputSection(mpsi) == MPS_ROWS
      || mpsinputSection(mpsi) == MPS_USERCUTS
      || mpsinputSection(mpsi) == MPS_LAZYCONS )
   {
      SCIP_CALL_TERMINATE( retcode, readRows(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_COLUMNS )
   {
      SCIP_CALL_TERMINATE( retcode, readCols(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_RHS )
   {
      SCIP_CALL_TERMINATE( retcode, readRhs(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_RANGES )
   {
      SCIP_CALL_TERMINATE( retcode, readRanges(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_BOUNDS )
   {
      SCIP_CALL_TERMINATE( retcode, readBounds(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_SOS )
   {
      SCIP_CALL_TERMINATE( retcode, readSOS(mpsi, scip), TERMINATE );
   }
   while( mpsinputSection(mpsi) == MPS_QCMATRIX )
   {
      SCIP_CALL_TERMINATE( retcode, readQCMatrix(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_QMATRIX )
   {
      SCIP_CALL_TERMINATE( retcode, readQMatrix(mpsi, FALSE, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_QUADOBJ )
   {
      SCIP_CALL_TERMINATE( retcode, readQMatrix(mpsi, TRUE, scip), TERMINATE );
   }
   while( mpsinputSection(mpsi) == MPS_QCMATRIX )
   {
      SCIP_CALL_TERMINATE( retcode, readQCMatrix(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) == MPS_INDICATORS )
   {
      SCIP_CALL_TERMINATE( retcode, readIndicators(mpsi, scip), TERMINATE );
   }
   if( mpsinputSection(mpsi) != MPS_ENDATA )
      mpsinputSyntaxerror(mpsi);

   SCIPfclose(fp);

   error = mpsinputHasError(mpsi);

   if( !error )
   {
      SCIP_CALL_TERMINATE( retcode, SCIPsetObjsense(scip, mpsinputObjsense(mpsi)), TERMINATE );
   }

 TERMINATE:
   mpsinputFree(scip, &mpsi);

   if( error )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}

/*
 * local methods for writing problem
 */

/** gets the key (i.e. the name) of the given namefreq */
static
SCIP_DECL_HASHGETKEY(hashGetKeyNamefreq)
{  /*lint --e{715}*/
   CONSNAMEFREQ* consnamefreq = (CONSNAMEFREQ*)elem;

   assert(consnamefreq != NULL);
   assert(consnamefreq->consname != NULL);

   return (void*)consnamefreq->consname;
}

/** returns TRUE iff both keys (i.e. strings) are equal up to max length*/
static
SCIP_DECL_HASHKEYEQ(hashKeyEqString)
{  /*lint --e{715}*/
   const char* string1 = (const char*)key1;
   const char* string2 = (const char*)key2;

   return (strncmp(string1, string2, MPS_MAX_NAMELEN - 1) == 0);
}

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
   if( key1 == key2 )
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


/** computes the field width such that the output file is nicely arranged */
static
unsigned int computeFieldWidth(
   unsigned int          width               /**< required width */
   )
{
   width = MAX(8u, width);
   return MIN(MPS_MAX_FIELDLEN, width);
}


/** output two strings in columns 1 and 2 with computed widths */
static
void printRecord(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           col1,               /**< column 1 */
   const char*           col2,               /**< column 2 */
   unsigned int          maxnamelen          /**< maximum name length */
   )
{
   unsigned int fieldwidth;
   char format[32];

   assert( scip != NULL );
   assert( col1 != NULL );
   assert( col2 != NULL );
   assert( strlen(col1) < MPS_MAX_NAMELEN );
   assert( strlen(col2) < MPS_MAX_VALUELEN );
   assert( maxnamelen > 0 );

   fieldwidth = computeFieldWidth(maxnamelen);
   (void) SCIPsnprintf(format, 32," %%-%ds %%%ds ", fieldwidth, MPS_MAX_VALUELEN - 1);

   SCIPinfoMessage(scip, file, (const char *)format, col1, col2);
}

/** output two strings in columns 1 (width 2) and 2 (width 8) */
static
void printStart(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           col1,               /**< column 1 */
   const char*           col2,               /**< column 2 */
   int                   maxnamelen          /**< maximum name length (-1 if irrelevant) */
   )
{
   unsigned int fieldwidth;
   char format[32];

   assert( scip != NULL );
   assert( col1 != NULL );
   assert( col2 != NULL );
   assert( strlen(col1) <= 2 );
   assert( strlen(col2) < MPS_MAX_NAMELEN );
   assert( maxnamelen == -1 || maxnamelen > 0 );

   if( maxnamelen < 0 )
   {
      /* format does not matter */
      (void) SCIPsnprintf(format, 32, " %%-2.2s %%-s ");
   }
   else
   {
      fieldwidth = computeFieldWidth((unsigned int) maxnamelen);
      (void) SCIPsnprintf(format, 32, " %%-2.2s %%-%ds ", fieldwidth);
   }

   SCIPinfoMessage(scip, file, (const char*)format, col1, col2);
}

/** prints the given data as column entry */
static
void printEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           varname,            /**< variable name */
   const char*           consname,           /**< constraint name */
   SCIP_Real             value,              /**< value to display */
   int*                  recordcnt,          /**< pointer to store the number of records per line */
   unsigned int          maxnamelen          /**< maximum name length */
   )
{
   char valuestr[MPS_MAX_VALUELEN] = { '\0' };

   assert( scip != NULL );
   assert( recordcnt != NULL );
   assert( *recordcnt >= 0 && *recordcnt < 2 );

   (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", value);

   if( *recordcnt == 0 )
   {
      /* start new line with an empty first column and the variable name in the second column */
      printStart(scip, file, "", varname, (int) maxnamelen);
      *recordcnt = 0;
   }

   printRecord(scip, file, consname, valuestr, maxnamelen);
   (*recordcnt)++;

   if( *recordcnt == 2 )
   {
      /* each line can have at most two records */
      SCIPinfoMessage(scip, file, "\n");
      *recordcnt = 0;
   }
}

/** prints the constraint type to file stream */
static
void printRowType(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   const char*           name                /**< constraint name */
   )
{
   char rowtype[2];

   assert( scip != NULL );
   assert( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) );
   assert( SCIPisGT(scip, rhs, lhs) || SCIPisEQ(scip, lhs, rhs) );
   assert( name != NULL );

   if( SCIPisEQ(scip, lhs, rhs) )
      (void) SCIPsnprintf(rowtype, 2, "%s", "E");
   else
   {
      /* in case the right hand side and the left hand side are not infinity we print a
       * less or equal constraint and put the right hand side in the RHS section and the
       * left hand side (hidden) in the RANGE section */
      if( !SCIPisInfinity(scip, rhs) )
         (void) SCIPsnprintf(rowtype, 2, "%s", "L");
      else
      {
         assert( !SCIPisInfinity(scip, -lhs) );
         (void) SCIPsnprintf(rowtype, 2, "%s", "G");
      }
   }

   printStart(scip, file, rowtype, name, -1);
   SCIPinfoMessage(scip, file, "\n");
}


/** initializes the sparse matrix */
static
SCIP_RETCODE initializeMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SPARSEMATRIX**        matrix,             /**< pointer to sparse matrix containing the entries */
   int                   slots               /**< number of slots */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, matrix) );
   (*matrix)->nentries = 0;
   (*matrix)->sentries = slots;
   SCIP_CALL( SCIPallocBufferArray(scip, &(*matrix)->values, (*matrix)->sentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*matrix)->columns, (*matrix)->sentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*matrix)->rows, (*matrix)->sentries) );

   return SCIP_OKAY;
}

/** this method takes care that the required capacity is available in the sparse matrix */
static
SCIP_RETCODE checkSparseMatrixCapacity(
   SCIP*                 scip,               /**< SCIP data structure */
   SPARSEMATRIX*         matrix,             /**< sparse matrix for storing the coefficient */
   int                   capacity            /**< needed capacity */
   )
{
   if( matrix->nentries + capacity >= matrix->sentries )
   {
      matrix->sentries = matrix->sentries * 2 + capacity;
      SCIP_CALL( SCIPreallocBufferArray(scip, &matrix->values, matrix->sentries) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &matrix->columns, matrix->sentries) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &matrix->rows, matrix->sentries) );
   }
   return SCIP_OKAY;
}

/** frees the sparse matrix */
static
void freeMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SPARSEMATRIX*         matrix              /**< sparse matrix to free */
   )
{
   SCIPfreeBufferArray(scip, &matrix->rows);
   SCIPfreeBufferArray(scip, &matrix->columns);
   SCIPfreeBufferArray(scip, &matrix->values);

   SCIPfreeBuffer(scip, &matrix);
}


/** computes the coefficient for the given variables and linear constraint information */
static
SCIP_RETCODE getLinearCoeffs(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           consname,           /**< name of the constraint */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_Real*            vals,               /**< array of coefficients values (or NULL if all coefficient values are 1) */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             transformed,        /**< transformed constraint? */
   SPARSEMATRIX*         matrix,             /**< sparse matrix for storing the coefficient */
   SCIP_Real*            rhs                 /**< pointer to right hand side */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   SCIP_Real activeconstant = 0.0;

   int nactivevars;
   int requiredsize;
   int v;

   assert( scip != NULL );
   assert( nvars == 0 || vars != NULL );
   assert( !SCIPisInfinity(scip, *rhs) );
   assert( matrix != NULL );

   /* if the variables array contains no variables, then return without
    * doing any thing; The MPS format and LP format do not forbid this
    * situation */
   if( nvars == 0 ) 
      return SCIP_OKAY;

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
   if( transformed )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, activevars, activevals, &nactivevars, nactivevars, &activeconstant, &requiredsize, TRUE) );

      if( requiredsize > nactivevars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &activevars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &activevals, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, activevars, activevals, &nactivevars, requiredsize, &activeconstant, &requiredsize, TRUE) );
         assert( requiredsize <= nactivevars );
      }
   }
   else
   {
      for( v = 0; v < nactivevars; ++v )
      {
         SCIP_CALL( SCIPvarGetOrigvarSum(&activevars[v], &activevals[v], &activeconstant) );

         /* negated variables with an original counterpart may also be returned by SCIPvarGetOrigvarSum();
          * make sure we get the original variable in that case
          */
         if( SCIPvarGetStatus(activevars[v]) == SCIP_VARSTATUS_NEGATED )
         {
            activevars[v] = SCIPvarGetNegatedVar(activevars[v]);
            activevals[v] *= -1.0;
            activeconstant += 1.0;
         }
      }
   }

   /* copy the (matrix) row into the sparse matrix */
   SCIP_CALL( checkSparseMatrixCapacity(scip, matrix, nactivevars) );
   assert( matrix->nentries + nactivevars < matrix->sentries );

   for( v = 0; v < nactivevars; ++v )
   {
      matrix->values[matrix->nentries] = activevals[v];
      matrix->columns[matrix->nentries] = activevars[v];
      matrix->rows[matrix->nentries] = consname;
      matrix->nentries++;
   }

   /* adjust right hand side */
   (*rhs) -= activeconstant;

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevals);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}


/** check whether given variables are aggregated and put them into an array without duplication */
static
SCIP_RETCODE collectAggregatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars,              /**< number of active variables in the problem */
   SCIP_VAR***           aggvars,            /**< pointer to array storing the aggregated variables on output */
   int*                  naggvars,           /**< pointer to number of aggregated variables on output */
   int*                  saggvars,           /**< pointer to number of slots in aggvars array */
   SCIP_HASHTABLE*       varAggregated       /**< hashtable for checking duplicates */
   )
{
   int v;

   assert( scip != NULL );
   assert( aggvars != NULL );
   assert( naggvars != NULL );
   assert( saggvars != NULL );

   /* check variables */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VARSTATUS status;
      SCIP_VAR* var;

      var = vars[v];
      status = SCIPvarGetStatus(var);

      /* collect aggregated variables in a list */
      if( status >= SCIP_VARSTATUS_AGGREGATED )
      {
         assert( status == SCIP_VARSTATUS_AGGREGATED || status == SCIP_VARSTATUS_MULTAGGR || status == SCIP_VARSTATUS_NEGATED );
         assert( varAggregated != NULL );

         if( ! SCIPhashtableExists(varAggregated, (void*) var) )
         {
            /* possibly enlarge array */
            if ( *saggvars <= *naggvars )
            {
               int newsize;
               newsize = SCIPcalcMemGrowSize(scip, *naggvars + 1);
               assert( newsize > *saggvars );
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggvars, *saggvars, newsize) );
               *saggvars = newsize;
            }

            (*aggvars)[*naggvars] = var;
            (*naggvars)++;
            SCIP_CALL( SCIPhashtableInsert(varAggregated, (void*) var) );
            assert( *naggvars <= *saggvars );
         }
      }
   }
   return SCIP_OKAY;
}


/** method check if the variable names are not longer than MPS_MAX_NAMELEN - 1*/
static
SCIP_RETCODE checkVarnames(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int                   nvars,              /**< number of variables */
   unsigned int*         maxnamelen,         /**< pointer to store the maximum name length */
   const char***         varnames,           /**< pointer to array of variable names */
   SCIP_HASHMAP**        varnameHashmap      /**< pointer to hash map storing variable, variable name mapping */      
   )
{
   int v;
   int faulty;
   char* varname;
   SCIP_VAR* var;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( maxnamelen != NULL );

   faulty = 0;

   /* allocate memory */
   SCIP_CALL( SCIPhashmapCreate(varnameHashmap, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, varnames, nvars) );

   /* check if the variable names are not to long */
   for( v = 0; v < nvars; ++v )
   {
      size_t l;

      var = vars[v];
      assert( var != NULL );

      l = strlen(SCIPvarGetName(var));

      if( l >= MPS_MAX_NAMELEN )
      {
         faulty++;
         (*maxnamelen) = MPS_MAX_NAMELEN - 1;
      }
      else
      {
         (*maxnamelen) = MAX(*maxnamelen, (unsigned int) l);
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &varname, (int) *maxnamelen + 1) );
      (void) SCIPsnprintf(varname, (int)(*maxnamelen) + 1, "%s", SCIPvarGetName(var) );

      /* insert variable with variable name into hash map */
      assert( !SCIPhashmapExists(*varnameHashmap, var) );
      SCIP_CALL( SCIPhashmapInsert(*varnameHashmap, var, (void*) varname) );

      (*varnames)[v] = varname;
   }

   if( faulty > 0 )
   {
      SCIPwarningMessage(scip, "there are %d variable names which have to be cut down to %d characters; LP might be corrupted\n",
         faulty, MPS_MAX_NAMELEN - 1);
   }
   return SCIP_OKAY;
}

/** method check if the constraint names are not longer than MPS_MAX_NAMELEN - 1 */
static
SCIP_RETCODE checkConsnames(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of all constraints */
   int                   nconss,             /**< number of all constraints */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   unsigned int*         maxnamelen,         /**< pointer to store the maximum name length */
   const char***         consnames,          /**< pointer to array of constraint names */
   SCIP_Bool*            error               /**< pointer to store whether all constraint names exist */
   )
{
   SCIP_HASHTABLE* consfreq;
   CONSNAMEFREQ* consnamefreqs;
   SCIP_CONS* cons;
   char* consname;
   int i;

   assert(scip != NULL);
   assert(maxnamelen != NULL);

   *error = FALSE;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &consnamefreqs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, consnames, nconss) );
   SCIP_CALL( SCIPhashtableCreate(&consfreq, SCIPblkmem(scip), SCIP_HASHSIZE_NAMES,
         hashGetKeyNamefreq, hashKeyEqString, SCIPhashKeyValString, NULL) );

   for( i = 0; i < nconss; ++i )
   {
      CONSNAMEFREQ* consnamefreq;
      size_t l;
      int freq;

      cons = conss[i];
      assert( cons != NULL );

      /* in case the transformed problem is written, only constraints are posted which are enabled in the current node */
      assert(!transformed || SCIPconsIsEnabled(cons));

      l = strlen(SCIPconsGetName(cons));

      if( l == 0 )
      {
         SCIPwarningMessage(scip, "At least one name of a constraint is empty, so file will be written with generic names.\n");
         *error = TRUE;

         goto TERMINATE;
      }

      consnamefreqs[i].consname = SCIPconsGetName(cons);
      consnamefreqs[i].freq = 0;
      freq = 0;

      /* check for duplicate names */
      if( NULL != (consnamefreq = (CONSNAMEFREQ *)SCIPhashtableRetrieve(consfreq, (void*)SCIPconsGetName(cons))) )
      {
         consnamefreq->freq += 1;
         consnamefreqs[i] = *consnamefreq;
         freq = consnamefreq->freq;
      }
      SCIP_CALL( SCIPhashtableInsert(consfreq, (void*)(&consnamefreqs[i])) );

      /* the new length is the length of the old name + a '_' and the freq number which has floor(log10(freq)) + 1 characters */
      if( freq > 0 )
         l = l + 1 + (size_t)log10((SCIP_Real) freq) + 1;

      if( l >= MPS_MAX_NAMELEN )
      {
         SCIPwarningMessage(scip, "Constraints have duplicate name and are too long to fix, so file will be written with generic names.\n");
         *error = TRUE;

         goto TERMINATE;
      }

      (*maxnamelen) = MAX(*maxnamelen, (unsigned int) l);

      SCIP_CALL( SCIPallocBufferArray(scip, &consname, (int) l + 1) );
      if( freq > 0 )
         (void) SCIPsnprintf(consname, (int)l + 1, "%s_%d", SCIPconsGetName(cons), freq);
      else
         (void) SCIPsnprintf(consname, (int)l + 1, "%s", SCIPconsGetName(cons));

      (*consnames)[i] = consname;
   }

TERMINATE:
   SCIPfreeBufferArray(scip, &consnamefreqs);

   if( *error )
   {
      --i;  /*lint !e445*/
      for( ; i >= 0; --i)  /*lint !e445*/
      {
         SCIPfreeBufferArray(scip, &((*consnames)[i]));
      }
      SCIPfreeBufferArray(scip, consnames);
   }

   SCIPhashtableFree(&consfreq);

   return SCIP_OKAY;
}


/** outputs the COLUMNS section of the MPS format */
static
void printColumnSection(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   SPARSEMATRIX*         matrix,             /**< sparse matrix containing the entries */
   SCIP_HASHMAP*         varnameHashmap,     /**< map from SCIP_VAR* to variable name */
   SCIP_HASHTABLE*       indicatorSlackHash, /**< hashtable containing slack variables from indicators (or NULL) */
   unsigned int          maxnamelen          /**< maximum name length */
   )
{
   SCIP_Bool intSection;
   SCIP_VAR* var;
   const char* varname;
   SCIP_Real value;
   int v;
   int recordcnt;

   /* sort sparse matrix w.r.t. the variable indices */
   SCIPsortPtrPtrReal((void**) matrix->columns, (void**) matrix->rows, matrix->values, SCIPvarComp, matrix->nentries);

   /* print COLUMNS section */
   SCIPinfoMessage(scip, file, "COLUMNS\n");

   intSection = FALSE;

   for( v = 0; v < matrix->nentries; )
   {
      var = matrix->columns[v];
      assert( var != NULL );

      /* skip slack variables in output */
      if( indicatorSlackHash != NULL && SCIPhashtableExists(indicatorSlackHash, var) )
      {
         ++v;
         continue;
      }

      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && intSection )
      {
         /* end integer section in MPS format */
         printStart(scip, file, "", "INTEND", (int) maxnamelen);
         printRecord(scip, file, "'MARKER'", "", maxnamelen);
         printRecord(scip, file, "'INTEND'", "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n", maxnamelen);
         intSection = FALSE;
      }
      else if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !intSection )
      {
         /* start integer section in MPS format */
         printStart(scip, file, "", "INTSTART", (int) maxnamelen);
         printRecord(scip, file, "'MARKER'", "", maxnamelen);
         printRecord(scip, file, "'INTORG'", "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         intSection = TRUE;
      }

      SCIPdebugMsg(scip, "create entries for variable <%s>\n", SCIPvarGetName(var));

      /* record count; there are at most two records per line */
      recordcnt = 0;

      /* get variable name */
      assert ( SCIPhashmapExists(varnameHashmap, var) );
      varname = (const char*) SCIPhashmapGetImage(varnameHashmap, var);

      /* output all entries of the same variable */
      do
      {
         value = matrix->values[v];

         /* print record to file */
         printEntry(scip, file, varname, matrix->rows[v], value, &recordcnt, maxnamelen);
         v++;
      }
      while( v < matrix->nentries && var == matrix->columns[v] );

      if( recordcnt == 1 )
         SCIPinfoMessage(scip, file, "\n");
   }
   /* end integer section, if the columns sections ends with integer variables */
   if( intSection )
   {
      /* end integer section in MPS format */
      printStart(scip, file, "", "INTEND", (int) maxnamelen);
      printRecord(scip, file, "'MARKER'", "", maxnamelen);
      printRecord(scip, file, "'INTEND'", "", maxnamelen);
      SCIPinfoMessage(scip, file, "\n", maxnamelen);
   }
}


/** outputs the right hand side section */
static
void printRhsSection(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   int                   nconss,             /**< number of constraints */
   const char**          consnames,          /**< constraint names */
   SCIP_Real*            rhss,               /**< right hand side array */
   unsigned int          maxnamelen          /**< maximum name length */
   )
{
   int recordcnt = 0;
   int c;

   assert( rhss != NULL );

   SCIPinfoMessage(scip, file, "RHS\n");
   SCIPdebugMsg(scip, "start printing RHS section\n");

   /* take care of the linear constraints */
   for( c = 0; c < nconss; ++c )
   {
      /* skip all constraints which have a right hand side of infinity */
      if( SCIPisInfinity(scip, rhss[c]) )
         continue;

      assert(consnames[c] != NULL);

      printEntry(scip, file, "RHS", consnames[c], rhss[c], &recordcnt, maxnamelen);
   }

   if( recordcnt == 1 )
      SCIPinfoMessage(scip, file, "\n");
}


/** outputs the range section */
static
void printRangeSection(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   SCIP_CONS**           conss,              /**< constraint array */
   int                   nconss,             /**< number of constraints */
   const char**          consnames,          /**< constraint names */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   unsigned int          maxnamelen          /**< maximum name length */
   )
{
   int c;
   int recordcnt = 0;

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   SCIP_CONS* cons;
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIPinfoMessage(scip, file, "RANGES\n");
   SCIPdebugMsg(scip, "start printing RANGES section\n");

   for( c = 0; c < nconss; ++c  )
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed problems is written only constraint are posted which are enabled in the current node;
       * the conss array should only contain relevant constraints
       */
      assert( !transformed || SCIPconsIsEnabled(cons) );

      assert( consnames[c] != NULL );

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         lhs = SCIPgetLhsLinear(scip, cons);
         rhs = SCIPgetRhsLinear(scip, cons);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         lhs = SCIPgetLhsVarbound(scip, cons);
         rhs = SCIPgetRhsVarbound(scip, cons);
      }
      else
         continue;

      if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, rhs, lhs) )
      {
         assert( SCIPisGT(scip, rhs, lhs) );
         printEntry(scip, file, "RANGE", consnames[c], rhs - lhs, &recordcnt, maxnamelen);
      }
   }
   if(recordcnt == 1 )
      SCIPinfoMessage(scip, file, "\n");
}

/** print bound section name */
static
void printBoundSectionName(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file, or NULL if standard output should be used */
   )
{
   SCIPinfoMessage(scip, file, "BOUNDS\n");
   SCIPdebugMsg(scip, "start printing BOUNDS section\n");
}

/** output bound section */
static
void printBoundSection(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL if standard output should be used */
   SCIP_VAR**            vars,               /**< active variables */
   int                   nvars,              /**< number of active variables */
   SCIP_VAR**            aggvars,            /**< needed aggregated variables */
   int                   naggvars,           /**< number of aggregated variables */
   SCIP_VAR**            fixvars,            /**< all fixed variables (or NULL if nfixvars is 0) */
   int                   nfixvars,           /**< number of fixed variables */
   SCIP_Bool             transformed,        /**< TRUE iff problem is the transformed problem */
   const char**          varnames,           /**< array with variable names */
   SCIP_HASHTABLE*       indicatorSlackHash, /**< hashtable containing slack variables from indicators (or NULL) */
   unsigned int          maxnamelen          /**< maximum name length */
   )
{
   int v;
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool sectionName;
   const char* varname;
   char valuestr[MPS_MAX_VALUELEN] = { '\0' };

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nfixvars == 0 || fixvars != NULL);

   sectionName = FALSE;

   /* output the active variables */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert( var != NULL );

      /* skip slack variables in output */
      if( indicatorSlackHash != NULL && SCIPhashtableExists(indicatorSlackHash, var) )
         continue;

      /* get variable name */
      varname = varnames[v];
      assert(strncmp(varname, SCIPvarGetName(var), maxnamelen) == 0);

      if( transformed )
      {
         /* in case the transformed is written only local bounds are posted
          * which are valid in the current node */
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
      }
      else
      {
         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);
      }

      /* take care of binary variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         if( !sectionName )
         {
            printBoundSectionName(scip, file);
            sectionName = TRUE;
         }

         if( !SCIPisFeasZero(scip, lb) || !SCIPisFeasEQ(scip, ub, 1.0) )
         {
            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", lb);
            printStart(scip, file, "LO", "Bound", (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n");

            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", ub);
            printStart(scip, file, "UP", "Bound", (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
         }
         else
         {
            printStart(scip, file, "BV", "Bound", (int) maxnamelen);
            printRecord(scip, file, varname, "", maxnamelen);
         }
         SCIPinfoMessage(scip, file, "\n");

         continue;
      }

      /* take care of free variables */
      if( SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub) )
      {
         if( !sectionName )
         {
            printBoundSectionName(scip, file);
            sectionName = TRUE;
         }

         /* variable is free */
         printStart(scip, file, "FR", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* take care of fixed variables */
      if( SCIPisEQ(scip, lb, ub) )
      {
         if( !sectionName )
         {
            printBoundSectionName(scip, file);
            sectionName = TRUE;
         }

         /* variable is fixed */
         (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", lb);
         printStart(scip, file, "FX", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, valuestr, maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
         continue;
      }

      /* print lower bound */
      if( SCIPisInfinity(scip, -lb) )
      {
         if( !sectionName )
         {
            printBoundSectionName(scip, file);
            sectionName = TRUE;
         }

         /* the free variables are processed above */
         assert( !SCIPisInfinity(scip, ub) );
         printStart(scip, file, "MI", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
      else
      {
         if( SCIPisZero(scip, lb) )
         {
            lb = 0.0;
         }
         else
         {
            if( !sectionName )
            {
               printBoundSectionName(scip, file);
               sectionName = TRUE;
            }

            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", lb);
            printStart(scip, file, "LO", "Bound", (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n");
         }
      }

      /* print upper bound, infinity has to be printed for integer (!) variables, because during
       * reading an mps file no upper bound of an integer variable means that the upper bound will
       * be set to 1 instead of +infinity (like it is for continuous variables) */
      if( SCIPisInfinity(scip, ub) )
      {
         if( !sectionName )
         {
            printBoundSectionName(scip, file);
            sectionName = TRUE;
         }

         /* the free variables are processed above */
         assert( !SCIPisInfinity(scip, -lb) );
         printStart(scip, file, "PL", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
      else
      {
         if( !sectionName )
         {
            printBoundSectionName(scip, file);
            sectionName = TRUE;
         }

         (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", ub);
         printStart(scip, file, "UP", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, valuestr, maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
   }

   /* output aggregated variables as 'free', except if they are binary */
   for( v = 0; v < naggvars; ++v )
   {
      if( !sectionName )
      {
         printBoundSectionName(scip, file);
         sectionName = TRUE;
      }

      var = aggvars[v];
      assert( var != NULL );

      /* get variable name */
      varname = varnames[nvars + v];
      assert(strncmp(varname, SCIPvarGetName(var), maxnamelen) == 0);

      /* take care of binary variables */
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         printStart(scip, file, "BV", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
      else
      {
         /* variable is free */
         printStart(scip, file, "FR", "Bound", (int) maxnamelen);
         printRecord(scip, file, varname, "", maxnamelen);
         SCIPinfoMessage(scip, file, "\n");
      }
   }

   /* output all fixed variables */
   for( v = 0; v < nfixvars; ++v )
   {
      /* we should print the transformed problem, otherwise no fixed variable should exists */
      assert(transformed);
      assert(fixvars != NULL && fixvars[v] != NULL);

      /* cppcheck-suppress nullPointer */
      var = fixvars[v];

      assert(var != NULL);
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);

      /* get variable name */
      varname = varnames[nvars + naggvars + v];
      assert(strncmp(varname, SCIPvarGetName(var), maxnamelen) == 0);

      /* only local bounds are posted which are valid in the current node */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      assert(SCIPisEQ(scip, lb, ub));

      if( !sectionName )
      {
         printBoundSectionName(scip, file);
         sectionName = TRUE;
      }

      /* print fixed variable */
      (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", lb);
      printStart(scip, file, "FX", "Bound", (int) maxnamelen);
      printRecord(scip, file, varname, valuestr, maxnamelen);
      SCIPinfoMessage(scip, file, "\n");
   }
}


/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
/**! [SnippetReaderCopyMps] */
static
SCIP_DECL_READERCOPY(readerCopyMps)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderMps(scip) );

   return SCIP_OKAY;
}
/**! [SnippetReaderCopyMps] */

/** destructor of reader to free user data (called when SCIP is exiting) */
/**! [SnippetReaderFreeMps] */
static
SCIP_DECL_READERFREE(readerFreeMps)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   SCIPfreeBlockMemory(scip, &readerdata);

   return SCIP_OKAY;
}
/**! [SnippetReaderFreeMps] */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadMps)
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   retcode = readMps(scip, filename);

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   if( retcode == SCIP_NOFILE || retcode == SCIP_READERROR )
      return retcode;

   SCIP_CALL( retcode );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteMps)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;
   int naddrows;
   int faulty = 0;
   int c;
   int v;
   int k;
   char* namestr;

   SCIP_CONS* cons = NULL;
   const char* consname;
   const char** consnames;

   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;

   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real* rhss;
   SCIP_Real value;

   SCIP_VAR* var = NULL;
   const char* varname;
   const char** varnames;

   char valuestr[MPS_MAX_VALUELEN] = { '\0' };

   SCIP_CONS** consIndicator;
   SCIP_CONS** consSOS1;
   SCIP_CONS** consSOS2;
   SCIP_CONS** consQuadratic;
   SCIP_CONS** consSOC;
   int nConsIndicator;
   int nConsSOS1;
   int nConsSOS2;
   int nConsQuadratic;
   int nConsSOC;

   SCIP_HASHMAP* varnameHashmap;           /* hash map from SCIP_VAR* to variable name */
   SPARSEMATRIX* matrix;

   SCIP_VAR** aggvars;
   int naggvars = 0;
   int saggvars;
   SCIP_HASHTABLE* varFixedHash;
   SCIP_HASHTABLE* indicatorSlackHash;

   SCIP_VAR** fixvars = NULL;
   int nfixvars = 0;

   SCIP_VAR** consvars;
   int nconsvars;
   SCIP_Real* vals;
   SCIP_Longint* weights;

   SCIP_Bool needRANGES;
   unsigned int maxnamelen;

   SCIP_Bool error;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   needRANGES = FALSE;
   maxnamelen = 0;
   nConsSOS1 = 0;
   nConsSOS2 = 0;
   nConsQuadratic = 0;
   nConsSOC = 0;
   nConsIndicator = 0;

   /* check if the constraint names are too long and build the constraint names */
   SCIP_CALL( checkConsnames(scip, conss, nconss, transformed, &maxnamelen, &consnames, &error) );
   if( error )
   {
      /* call writing with generic names */
      if( transformed )
      {
         SCIPwarningMessage(scip, "write transformed problem with generic variable and constraint names\n");
         SCIP_CALL( SCIPprintTransProblem(scip, file, "mps", TRUE) );
      }
      else
      {
         SCIPwarningMessage(scip, "write original problem with generic variable and constraint names\n");
         SCIP_CALL( SCIPprintOrigProblem(scip, file, "mps", TRUE) );
      }
      *result = SCIP_SUCCESS;

      return SCIP_OKAY;
   }

   /* check if the variable names are not too long and build the "variable" -> "variable name" hash map */
   SCIP_CALL( checkVarnames(scip, vars, nvars, &maxnamelen, &varnames, &varnameHashmap) );

   /* collect SOS, quadratic, and indicator constraints in array for later output */
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS1, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOS2, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consQuadratic, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consSOC, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consIndicator, nconss) );

   /* nfixedvars counts all variables with status SCIP_VARSTATUS_FIXED, SCIP_VARSTATUS_AGGREGATED, SCIP_VARSTATUS_MULTAGGR, but not SCIP_VARSTATUS_NEGATED */
   saggvars = nfixedvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &aggvars, saggvars) );

   /* create hashtable for storing aggregated variables */
   if( nfixedvars > 0 )
   {
      SCIP_CALL( SCIPhashtableCreate(&varFixedHash, SCIPblkmem(scip), nfixedvars, hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   }
   else
      varFixedHash = NULL;

   if( nvars > 0 )
   {
      SCIP_CALL( SCIPhashtableCreate(&indicatorSlackHash, SCIPblkmem(scip), nvars, hashGetKeyVar, hashKeyEqVar, hashKeyValVar, NULL) );
   }
   else
      indicatorSlackHash = NULL;

   /* initialize sparse matrix */
   SCIP_CALL( initializeMatrix(scip, &matrix, (nvars * 2 + nfixedvars)) );
   assert( matrix->sentries >= nvars );

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   naddrows = 0;

   /* determine and-constraints and printing format to resize necessary arrays */
   if( readerdata->linearizeands )
   {
      SCIP_CONSHDLR* andconshdlr = SCIPfindConshdlr(scip, "and");

      if( andconshdlr != NULL )
      {
         /* need to check for and-constraints, note that in the original problem you cannot get the number of
          * and-constraints by one call */
         for( c = nconss - 1; c >= 0; --c )
         {
            conshdlr = SCIPconsGetHdlr(conss[c]);
            assert(conshdlr != NULL);

            conshdlrname = SCIPconshdlrGetName(conshdlr);

            if( strcmp(conshdlrname, "and") == 0 )
            {
               if( readerdata->aggrlinearizationands )
                  ++naddrows;
               else
                  naddrows += SCIPgetNVarsAnd(scip, conss[c]);
            }
         }
         assert(naddrows >= 0);

         if( naddrows > 0 )
         {
            /* resize consnames vector */
            SCIP_CALL( SCIPreallocBufferArray(scip, &consnames, nconss + naddrows) );
         }
      }
   }

   /* initialize rhs vector */
   SCIP_CALL( SCIPallocBufferArray(scip, &rhss, nconss + naddrows) );

   /* print statistics as comment to file stream */
   SCIPinfoMessage(scip, file, "* SCIP STATISTICS\n");
   SCIPinfoMessage(scip, file, "*   Problem name     : %s\n", name);
   SCIPinfoMessage(scip, file, "*   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   SCIPinfoMessage(scip, file, "*   Constraints      : %d\n", nconss);
   SCIPinfoMessage(scip, file, "*   Obj. scale       : %.15g\n", objscale);
   SCIPinfoMessage(scip, file, "*   Obj. offset      : %.15g\n", objoffset);

   /* print NAME of the problem */
   SCIPinfoMessage(scip, file, "%-14s%s\n", "NAME", name);

   /* print OBJSENSE of the problem */
   SCIPinfoMessage(scip, file, "OBJSENSE\n");
   SCIPinfoMessage(scip, file, "%s\n", objsense == SCIP_OBJSENSE_MAXIMIZE ? "  MAX" : "  MIN");

   /* start ROWS section */
   SCIPinfoMessage(scip, file, "ROWS\n");

   /* print row type for the objective function */
   printStart(scip, file, "N", "Obj", -1);
   SCIPinfoMessage(scip, file, "\n");

   /* first fill the matrix with the objective coefficients */
   for( v = 0; v < nvars; ++v )
   {
      /* take care of the objective entry */
      var = vars[v];
      value = SCIPvarGetObj(var);

      /* we also want to add integer variables to the columns section, even if the objective value is 0, because it
       * might happen that they only exist in non-linear constraints, which leads to no other line in the column section
       * and therefore do not mark the variable as an integer
       */
      if( !SCIPisZero(scip, value) || SCIPvarGetType(var) < SCIP_VARTYPE_IMPLINT || ((SCIPvarGetNLocksDown(var) == 0) && (SCIPvarGetNLocksUp(var) == 0)) )
      {
         assert( matrix->nentries < matrix->sentries );

         matrix->values[matrix->nentries] = value;
         matrix->columns[matrix->nentries] = var;
         matrix->rows[matrix->nentries] = "Obj";
         matrix->nentries++;
      }
   }

   /* loop over all constraints */
   k = nconss;
   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert( cons != NULL);

      /* in case the transformed problems is written only constraint are posted which are enabled in the current node;
       * the conss array should only contain relevant constraints
       */
      assert( !transformed || SCIPconsIsEnabled(cons) );

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);

      /* construct constraint name */
      consname = consnames[c];

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         lhs = SCIPgetLhsLinear(scip, cons);
         rhs = SCIPgetRhsLinear(scip, cons);

         /* there is nothing to do if the left hand side is minus infinity and the right side is infinity */
         if( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) )
         {
            if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
               needRANGES = TRUE;

            /* print row entry */
            printRowType(scip, file, lhs, rhs, consname);

            if( SCIPisInfinity(scip, rhs) )
               rhss[c] = lhs;
            else
               rhss[c] = rhs;

            assert( !SCIPisInfinity(scip, rhss[c]) );

            /* compute column entries */
            SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
                  SCIPgetNVarsLinear(scip, cons), transformed, matrix, &rhss[c]) );
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         /* print row entry */
         switch( SCIPgetTypeSetppc(scip, cons) )
         {
         case SCIP_SETPPCTYPE_PARTITIONING :
            printRowType(scip, file, 1.0, 1.0, consname);
            break;
         case SCIP_SETPPCTYPE_PACKING :
            printRowType(scip, file, -SCIPinfinity(scip), 1.0, consname);
            break;
         case SCIP_SETPPCTYPE_COVERING :
            printRowType(scip, file, 1.0, SCIPinfinity(scip), consname);
            break;
         }

         rhss[c] = 1.0;

         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetVarsSetppc(scip, cons), NULL, SCIPgetNVarsSetppc(scip, cons), transformed, matrix, &rhss[c]) );
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         /* print row entry */
         printRowType(scip, file, 1.0, SCIPinfinity(scip), consname);

         rhss[c] = 1.0;

         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetVarsLogicor(scip, cons), NULL, SCIPgetNVarsLogicor(scip, cons), transformed, matrix, &rhss[c]) );
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         int i;

         /* print row entry */
         printRowType(scip, file, -SCIPinfinity(scip), (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons), consname);

         nconsvars = SCIPgetNVarsKnapsack(scip, cons);
         weights = SCIPgetWeightsKnapsack(scip, cons);

         /* copy Longint array to SCIP_Real array */
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, nconsvars ) );
         for( i = 0; i < nconsvars; ++i )
            vals[i] = (SCIP_Real)weights[i];

         rhss[c] = (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons);

         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetVarsKnapsack(scip, cons), vals, nconsvars, transformed, matrix, &rhss[c]) );

         SCIPfreeBufferArray(scip, &vals);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         lhs = SCIPgetLhsVarbound(scip, cons);
         rhs = SCIPgetRhsVarbound(scip, cons);

         /* there is nothing to do if the left hand side is minus infinity and the right side is infinity */
         if( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) )
         {
            if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
               needRANGES = TRUE;

            /* print row entry */
            printRowType(scip, file, lhs, rhs, consname);

            /* allocate memory */
            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

            consvars[0] = SCIPgetVarVarbound(scip, cons);
            consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

            vals[0] = 1.0;
            vals[1] = SCIPgetVbdcoefVarbound(scip, cons);

            if( SCIPisInfinity(scip, rhs) )
               rhss[c] = lhs;
            else
               rhss[c] = rhs;

            assert( !SCIPisInfinity(scip, rhss[c]) );

            /* compute column entries */
            SCIP_CALL( getLinearCoeffs(scip, consname, consvars, vals, 2, transformed, matrix, &rhss[c]) );

            SCIPfreeBufferArray(scip, &vals);
            SCIPfreeBufferArray(scip, &consvars);
         }
      }
      else if( strcmp(conshdlrname, "indicator") == 0 )
      {
         SCIP_VAR* slackvar;
         SCIP_VAR* binvar;

         /* store slack variable in hash */
         slackvar = SCIPgetSlackVarIndicator(cons);
         assert( slackvar != NULL );
         assert( indicatorSlackHash != NULL );
         assert( !SCIPhashtableExists(indicatorSlackHash, (void*) slackvar) );
         SCIP_CALL( SCIPhashtableInsert(indicatorSlackHash, (void*) slackvar) );

         /* if slackvariable is aggregated, we store it in the list of aggregated variables */
         if ( SCIPvarGetStatus(slackvar) == SCIP_VARSTATUS_AGGREGATED )
         {
            SCIP_CALL( collectAggregatedVars(scip, &slackvar, 1, &aggvars, &naggvars, &saggvars, varFixedHash) );
         }

         /* store aggregated variables */
         binvar = SCIPgetBinaryVarIndicator(cons);
         if( SCIPvarIsNegated(binvar) )
            binvar = SCIPvarGetNegatedVar(binvar);
         assert( binvar != NULL );
         SCIP_CALL( collectAggregatedVars(scip, &binvar, 1, &aggvars, &naggvars, &saggvars, varFixedHash) );

         /* indicator constraint do not have a right hand side; mark this with SCIPinfinity(scip) */
         rhss[c] = SCIPinfinity(scip);

         /* store constraint */
         consIndicator[nConsIndicator++] = cons;
         continue;
      }
      else if( strcmp(conshdlrname, "SOS1") == 0 )
      {
         /* store constraint */
         consSOS1[nConsSOS1++] = cons;

         /* check for aggregated variables in SOS1 constraints for later output
          * of aggregations as linear constraints */
         consvars = SCIPgetVarsSOS1(scip, cons);
         nconsvars = SCIPgetNVarsSOS1(scip, cons);

         /* SOS constraint do not have a right hand side; mark this with SCIPinfinity(scip) */
         rhss[c] = SCIPinfinity(scip);

         SCIP_CALL( collectAggregatedVars(scip, consvars, nconsvars, &aggvars, &naggvars, &saggvars, varFixedHash) );
      }
      else if( strcmp(conshdlrname, "SOS2") == 0 )
      {
         /* store constraint */
         consSOS2[nConsSOS2++] = cons;

         /* check for aggregated variables in SOS2 constraints for later output aggregations as linear constraints */
         consvars = SCIPgetVarsSOS2(scip, cons);
         nconsvars = SCIPgetNVarsSOS2(scip, cons);

         /* SOS constraint do not have a right hand side; mark this with SCIPinfinity(scip) */
         rhss[c] = SCIPinfinity(scip);

         SCIP_CALL( collectAggregatedVars(scip, consvars, nconsvars, &aggvars, &naggvars, &saggvars, varFixedHash) );
      }
      else if( strcmp(conshdlrname, "quadratic") == 0 )
      {
         SCIP_VAR** quadvars;
         SCIP_Real* quadvarlincoefs;
         int j;

         /* store constraint */
         consQuadratic[nConsQuadratic++] = cons;

         /* collect linear coefficients of quadratic part */
         SCIP_CALL( SCIPallocBufferArray(scip, &quadvars, SCIPgetNQuadVarTermsQuadratic(scip, cons)) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadvarlincoefs, SCIPgetNQuadVarTermsQuadratic(scip, cons)) );
         for( j = 0; j < SCIPgetNQuadVarTermsQuadratic(scip, cons); ++j )
         {
            quadvars[j]        = SCIPgetQuadVarTermsQuadratic(scip, cons)[j].var;
            quadvarlincoefs[j] = SCIPgetQuadVarTermsQuadratic(scip, cons)[j].lincoef;
         }

         lhs = SCIPgetLhsQuadratic(scip, cons);
         rhs = SCIPgetRhsQuadratic(scip, cons);

         /* there is nothing to do if the left hand side is minus infinity and the right side is infinity */
         if( !SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs) )
         {
            if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && !SCIPisEQ(scip, lhs, rhs) )
               needRANGES = TRUE;

            /* print row entry */
            printRowType(scip, file, lhs, rhs, consname);

            if( SCIPisInfinity(scip, rhs) )
               rhss[c] = lhs;
            else
               rhss[c] = rhs;

            assert( !SCIPisInfinity(scip, rhss[c]) );

            /* compute column entries for linear part */
            SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
                  SCIPgetNLinearVarsQuadratic(scip, cons), transformed, matrix, &rhss[c]) );

            /* compute column entries for linear part in quadratic part */
            SCIP_CALL( getLinearCoeffs(scip, consname, quadvars, quadvarlincoefs, SCIPgetNQuadVarTermsQuadratic(scip, cons),
                  transformed, matrix, &rhss[c]) );
         }

         /* check for aggregated variables in quadratic part of quadratic constraints for later output of
          * aggregations as linear constraints */
         consvars = quadvars;
         nconsvars = SCIPgetNQuadVarTermsQuadratic(scip, cons);

         SCIP_CALL( collectAggregatedVars(scip, consvars, nconsvars, &aggvars, &naggvars, &saggvars, varFixedHash) );

         SCIPfreeBufferArray(scip, &quadvars);
         SCIPfreeBufferArray(scip, &quadvarlincoefs);
      }
      else if( strcmp(conshdlrname, "soc") == 0 )
      {
         /* SOC constraints are of the form lhsconstant + sum_i (lhscoef_i*(lhsvar_i+lhsoffset_i))^2 <= (rhscoef*(rhsvar+rhsoffset))^2 */
         SCIP_Real* lincoefs;
         SCIP_Real  coef;
         SCIP_Real  offset;

         /* store constraint */
         consSOC[nConsSOC++] = cons;

         consvars  = SCIPgetLhsVarsSOC(scip, cons);
         nconsvars = SCIPgetNLhsVarsSOC(scip, cons);

         rhs = -SCIPgetLhsConstantSOC(scip, cons);

         /* offsets on lhs give linear coefficients that need to be processed here */
         SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nconsvars) );

         for( v = 0; v < nconsvars; ++v )
         {
            offset = SCIPgetLhsOffsetsSOC(scip, cons)[v];
            coef = SCIPgetLhsCoefsSOC(scip, cons)[v];

            lincoefs[v] = 2 * offset * coef * coef;
            rhs -= offset * offset * coef * coef;
         }

         SCIP_CALL( getLinearCoeffs(scip, consname, SCIPgetLhsVarsSOC(scip, cons), lincoefs, nconsvars, transformed, matrix, &rhs) );

         SCIPfreeBufferArray(scip, &lincoefs);

         /* if there is an offsets on rhs, then we have linear a coefficient that need to be processed here */
         if( SCIPgetRhsOffsetSOC(scip, cons) != 0.0 )
         {
            SCIP_VAR* rhsvar;
            SCIP_Real lincoef;

            coef   = SCIPgetRhsCoefSOC(scip, cons);
            offset = SCIPgetRhsOffsetSOC(scip, cons);
            rhsvar = SCIPgetRhsVarSOC(scip, cons);
            lincoef = -2 * offset * coef * coef;
            rhs += offset * offset * coef * coef;

            SCIP_CALL( getLinearCoeffs(scip, consname, &rhsvar, &lincoef, 1, transformed, matrix, &rhs) );
         }

         assert(!SCIPisInfinity(scip, ABS(rhs)));

         /* print row entry */
         printRowType(scip, file, -SCIPinfinity(scip), rhs, consname);

         rhss[c] = rhs;

         /* check for aggregated variables in for later output of aggregations as linear constraints */
         SCIP_CALL( collectAggregatedVars(scip, consvars, nconsvars, &aggvars, &naggvars, &saggvars, varFixedHash) );
         var = SCIPgetRhsVarSOC(scip, cons);
         SCIP_CALL( collectAggregatedVars(scip, &var, 1, &aggvars, &naggvars, &saggvars, varFixedHash) );
      }
      else if( strcmp(conshdlrname, "and") == 0 )
      {
         if( readerdata->linearizeands )
         {
            SCIP_VAR** rowvars;
            SCIP_VAR** operands;
            SCIP_VAR* resultant;
            SCIP_Real* rowvals;
            char* rowname;
            int nrowvars;
            int l;
            int n;

            nrowvars = SCIPgetNVarsAnd(scip, cons);
            operands = SCIPgetVarsAnd(scip, cons);
            resultant = SCIPgetResultantAnd(scip, cons);

            /* allocate buffer array */
            SCIP_CALL( SCIPallocBufferArray(scip, &rowvars, nrowvars + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &rowvals, nrowvars + 1) );

            /* get length of constraint name */
            l = (int) strlen(consname);

            /* the tight relaxtion, number of and-constraint operands rows */
            if( !readerdata->aggrlinearizationands )
            {
               rowvars[0] = resultant;
               rowvals[0] = 1.0;
               rowvals[1] = -1.0;

               /* compute maximal length for rowname */
               n = (int) log10((double)nrowvars) + 1 + l;

               /* assure maximal allowed value */
               if( n >= MPS_MAX_NAMELEN )
                  n = MPS_MAX_NAMELEN - 1;

               /* update maxnamelen */
               maxnamelen = MAX(maxnamelen, (unsigned int) n);

               /* print operator rows */
               for( v = 0; v < nrowvars; ++v )
               {
                  /* compute maximal length for rowname */
                  if( v == 0 )
                     n = 2;
                  else
                     n = (int) log10((double)v) + 2;
                  n += l;

                  /* assure maximal allowed value */
                  if( n >= MPS_MAX_NAMELEN )
                  {
                     n = MPS_MAX_NAMELEN - 1;
                     ++faulty;
                  }

                  /* need memory for additional row */
                  SCIP_CALL( SCIPallocBufferArray(scip, &rowname, n + 1) );

                  assert(k < nconss + naddrows);
                  consnames[k] = rowname;

                  (void) SCIPsnprintf(rowname, n + 1, "%s_%d", consname, v);
                  rowvars[1] = operands[v];

                  /* print row entry */
                  printRowType(scip, file, -SCIPinfinity(scip), 0.0, rowname);

                  rhss[k] = 0.0;

                  /* compute column entries */
                  SCIP_CALL( getLinearCoeffs(scip, rowname, rowvars, rowvals, 2, transformed, matrix, &rhss[k]) );
                  ++k;
               }
            }

            /* prepare for next row */
            for( v = nrowvars - 1; v >= 0; --v )
            {
               rowvars[v] = operands[v];
               rowvals[v] = -1.0;
            }

            rowvars[nrowvars] = resultant;

            /* the weak relaxtion, only one constraint */
            if( readerdata->aggrlinearizationands )
            {
               /* compute maximal length for rowname */
               n = l + 3;

               /* assure maximal allowed value */
               if( n >= MPS_MAX_NAMELEN )
               {
                  n = MPS_MAX_NAMELEN - 1;
                  ++faulty;
               }

               /* update maxnamelen */
               maxnamelen = MAX(maxnamelen, (unsigned int) n);

               /* need memory for additional row */
               SCIP_CALL( SCIPallocBufferArray(scip, &rowname, n + 1) );

               assert(k < nconss + naddrows);
               consnames[k] = rowname;

               /* adjust rowname of constraint */
               (void) SCIPsnprintf(rowname, n + 1, "%s_op", consname);

               rowvals[nrowvars] = (SCIP_Real) nrowvars;

               /* print row entry */
               printRowType(scip, file, -SCIPinfinity(scip), 0.0, rowname);

               rhss[k] = 0.0;

               /* compute column entries */
               SCIP_CALL( getLinearCoeffs(scip, rowname, rowvars, rowvals, nrowvars + 1, transformed, matrix, &rhss[k]) );

               SCIPdebugMsg(scip, "%g, %g\n", rowvals[1], rhss[k]);
               ++k;
            }

            rowvals[nrowvars] = 1.0;

            /* print row entry */
            printRowType(scip, file, -nrowvars + 1.0, SCIPinfinity(scip), consname);

            rhss[c] = -nrowvars + 1.0;

            /* compute column entries */
            SCIP_CALL( getLinearCoeffs(scip, consname, rowvars, rowvals, nrowvars + 1, transformed, matrix, &rhss[c]) );

            /* free buffer array */
            SCIPfreeBufferArray(scip, &rowvals);
            SCIPfreeBufferArray(scip, &rowvars);
         }
         else
         {
            /* and constraint printing not enabled; mark this with SCIPinfinity(scip) */
            rhss[c] = SCIPinfinity(scip);

            SCIPwarningMessage(scip, "change parameter \"reading/" READER_NAME "/linearize-and-constraints\" to TRUE to print and-constraints\n");
         }
      }
      else
      {
         /* unknown constraint type; mark this with SCIPinfinity(scip) */
         rhss[c] = SCIPinfinity(scip);

         SCIPwarningMessage(scip, "constraint handler <%s> cannot print requested format\n", conshdlrname );
      }
   }

   if( faulty > 0 )
   {
      SCIPwarningMessage(scip, "there are %d and-constraint-rownames which have to be cut down to %d characters; MPS file might be corrupted\n",
         faulty, MPS_MAX_NAMELEN - 1);
   }

   /* free hash table */
   if( varFixedHash != NULL )
      SCIPhashtableFree(&varFixedHash);

   if( indicatorSlackHash != NULL && nConsIndicator == 0 )
   {
      SCIPhashtableFree(&indicatorSlackHash);
      assert( indicatorSlackHash == NULL );
   }

   if( naggvars > 0 )
   {
      /* construct variables name of the needed aggregated variables and the constraint names for the aggregation constraints */

      /* realloc memory */
      SCIP_CALL( SCIPreallocBufferArray(scip, &consnames, nconss + naddrows + naggvars) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &rhss, nconss + naddrows + naggvars) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &varnames, nvars + naggvars) );

      for( c = 0; c < naggvars; ++c )
      {
         size_t l;

         /* create variable name */
         var = aggvars[c];

         l = strlen(SCIPvarGetName(var));
         if( l >= MPS_MAX_NAMELEN )
            maxnamelen = MPS_MAX_NAMELEN - 1;
         else
            maxnamelen = MAX(maxnamelen, (unsigned int) l);

         SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );
         (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

         /* insert variable with variable name into hash map */
         varnames[nvars + c] = namestr;
         assert( !SCIPhashmapExists(varnameHashmap, var) );
         SCIP_CALL( SCIPhashmapInsert(varnameHashmap, var, (void*) namestr) );

         /* output row type (it is an equation) */
         SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) ); /* note that namestr above is freed via varnames */
         (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "aggr_%s", SCIPvarGetName(var));
         printRowType(scip, file, 1.0, 1.0, namestr);

         l = strlen(namestr);
         maxnamelen = MAX(maxnamelen, (unsigned int) l);
         consnames[nconss + naddrows + c] = namestr;
         rhss[nconss + naddrows + c] = 0.0;

         /* compute column entries */
         SCIP_CALL( getLinearCoeffs(scip, namestr, &(aggvars[c]), NULL, 1, transformed, matrix, &rhss[nconss + naddrows + c]) );

         /* add the aggregated variables to the sparse matrix */
         SCIP_CALL( checkSparseMatrixCapacity(scip, matrix, 1) );
         matrix->values[matrix->nentries] = -1.0;
         matrix->columns[matrix->nentries] = aggvars[c];
         matrix->rows[matrix->nentries] = namestr;
         matrix->nentries++;
      }
   }

   /* collect also fixed variables, because they might not be removed from all constraints */
   /* @todo only collect fixed variables in the non-linear constraint types, where they (could not be)/(were not) removed */
   if( nfixedvars > 0 )
   {
      int startpos = nvars + naggvars;
      /* construct variables name of fixed variables */

      /* realloc memory */
      SCIP_CALL( SCIPreallocBufferArray(scip, &varnames, startpos + nfixedvars) );

      /* allocate memory for fixed variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &fixvars, nfixedvars) );

      for( v = nfixedvars - 1; v >= 0; --v )
      {
         /* create variable name */
         var = fixedvars[v];

	 if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
	 {
            size_t l;
	    l = strlen(SCIPvarGetName(var));
	    if( l >= MPS_MAX_NAMELEN )
	       maxnamelen = MPS_MAX_NAMELEN - 1;
	    else
	       maxnamelen = MAX(maxnamelen, (unsigned int) l);

	    SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );
	    (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPvarGetName(var) );

	    varnames[startpos + nfixvars] = namestr;
	    fixvars[nfixvars] = var;
	    ++nfixvars;

	    /* insert variable with variable name into hash map */
	    assert(!SCIPhashmapExists(varnameHashmap, var));
	    SCIP_CALL( SCIPhashmapInsert(varnameHashmap, var, (void*) namestr) );

	    /* add the fixed variables to the sparse matrix, needed for columns section */
	    SCIP_CALL( checkSparseMatrixCapacity(scip, matrix, 1) );
	    matrix->values[matrix->nentries] = 0.0;
	    matrix->columns[matrix->nentries] = var;
	    matrix->rows[matrix->nentries] = "Obj";
	    matrix->nentries++;
	 }
      }
   }

   /* output COLUMNS section */
   printColumnSection(scip, file, matrix, varnameHashmap, indicatorSlackHash, maxnamelen);

   /* output RHS section */
   printRhsSection(scip, file, nconss + naddrows +naggvars, consnames, rhss, maxnamelen);

   /* output RANGES section */
   if( needRANGES )
      printRangeSection(scip, file, conss, nconss, consnames, transformed, maxnamelen);

   /* output BOUNDS section */
   printBoundSection(scip, file, vars, nvars, aggvars, naggvars, fixvars, nfixvars, transformed, varnames, indicatorSlackHash, maxnamelen);

   if( nfixedvars > 0 )
   {
      SCIPfreeBufferArray(scip, &fixvars);
   }

   /* print SOS section */
   if( nConsSOS1 > 0 || nConsSOS2 > 0 )
   {
      SCIP_Real* sosweights;

      SCIPinfoMessage(scip, file, "SOS\n");
      SCIPdebugMsg(scip, "start printing SOS section\n");

      SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );

      /* first output SOS1 constraints */
      for( c = 0; c < nConsSOS1; ++c )
      {
         cons = consSOS1[c];
         consvars = SCIPgetVarsSOS1(scip, cons);
         nconsvars = SCIPgetNVarsSOS1(scip, cons);
         sosweights = SCIPgetWeightsSOS1(scip, cons);
         (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         printStart(scip, file, "S1", namestr, -1);
         SCIPinfoMessage(scip, file, "\n");

         for( v = 0; v < nconsvars; ++v )
         {
            /* get variable name */
            assert ( SCIPhashmapExists(varnameHashmap, consvars[v]) );
            varname = (const char*) SCIPhashmapGetImage(varnameHashmap, consvars[v]);

            printStart(scip, file, "", varname, (int) maxnamelen);

            if( sosweights != NULL )
               (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", sosweights[v]);
            else
               (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25d ", v);

            SCIPinfoMessage(scip, file, "%25s\n", valuestr);
         }
      }

      /* next output SOS2 constraints */
      for( c = 0; c < nConsSOS2; ++c )
      {
         cons = consSOS2[c];
         consvars = SCIPgetVarsSOS2(scip, cons);
         nconsvars = SCIPgetNVarsSOS2(scip, cons);
         sosweights = SCIPgetWeightsSOS2(scip, cons);
         (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         printStart(scip, file, "S2", namestr, -1);
         SCIPinfoMessage(scip, file, "\n");

         for( v = 0; v < nconsvars; ++v )
         {
            /* get variable name */
            assert ( SCIPhashmapExists(varnameHashmap, consvars[v]) );
            varname = (const char*) SCIPhashmapGetImage(varnameHashmap, consvars[v]);

            printStart(scip, file, "", varname, (int) maxnamelen);

            if( sosweights != NULL )
               (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", sosweights[v]);
            else
               (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25d ", v);

            SCIPinfoMessage(scip, file, "%25s\n", valuestr);
         }
      }
      SCIPfreeBufferArray(scip, &namestr);
   }

   /* print QCMATRIX sections for quadratic constraints
    * in difference to a quadratic term in the objective function, the quadratic part is not divided by 2 here
    */
   if( nConsQuadratic > 0 )
   {
      SCIP_QUADVARTERM* quadvarterms;
      SCIP_BILINTERM*   bilinterms;
      const char* varname2;
      int nbilin;

      SCIPdebugMsg(scip, "start printing QCMATRIX sections for quadratic constraints\n");
      SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );

      for( c = 0; c < nConsQuadratic; ++c )
      {
         cons = consQuadratic[c];
         nconsvars = SCIPgetNQuadVarTermsQuadratic(scip, cons);
         quadvarterms = SCIPgetQuadVarTermsQuadratic(scip, cons);
         bilinterms = SCIPgetBilinTermsQuadratic(scip, cons);
         nbilin = SCIPgetNBilinTermsQuadratic(scip, cons);

         (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );

         SCIPinfoMessage(scip, file, "QCMATRIX %s\n", namestr);

         /* print x^2 terms */
         for( v = 0; v < nconsvars; ++v )
         {
            if( quadvarterms[v].sqrcoef == 0.0 )
               continue;

            /* get variable name */
            assert ( SCIPhashmapExists(varnameHashmap, quadvarterms[v].var) );
            varname = (const char*) SCIPhashmapGetImage(varnameHashmap, quadvarterms[v].var);

            /* get coefficient as string */
            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", quadvarterms[v].sqrcoef);

            /* print "x x coeff" line */
            printStart(scip, file, "", varname, (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n", valuestr);
         }

         /* print bilinear terms; CPLEX format expects a symmetric matrix with all coefficients specified,
          * i.e., we have to split bilinear coefficients into two off diagonal elements */
         for( v = 0; v < nbilin; ++v )
         {
            if( bilinterms[v].coef == 0.0 )
               continue;

            /* get name of first variable */
            assert ( SCIPhashmapExists(varnameHashmap, bilinterms[v].var1) );
            varname = (const char*) SCIPhashmapGetImage(varnameHashmap, bilinterms[v].var1);

            /* get name of second variable */
            assert ( SCIPhashmapExists(varnameHashmap, bilinterms[v].var2) );
            varname2 = (const char*) SCIPhashmapGetImage(varnameHashmap, bilinterms[v].var2);

            /* get coefficient as string */
            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", 0.5*bilinterms[v].coef);

            /* print "x y coeff/2" line */
            printStart(scip, file, "", varname, (int) maxnamelen);
            printRecord(scip, file, varname2, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n", valuestr);

            /* print "y x coeff/2" line */
            printStart(scip, file, "", varname2, (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n", valuestr);
         }
      }

      SCIPfreeBufferArray(scip, &namestr);
   }

   /* print QCMATRIX sections for second order cone constraints */
   if( nConsSOC > 0 )
   {
      SCIP_Real* coefs;

      SCIPdebugMsg(scip, "start printing QCMATRIX sections for soc constraints\n");
      SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );

      for( c = 0; c < nConsSOC; ++c )
      {
         cons = consSOC[c];
         consvars = SCIPgetLhsVarsSOC(scip, cons);
         nconsvars = SCIPgetNLhsVarsSOC(scip, cons);
         coefs = SCIPgetLhsCoefsSOC(scip, cons);

         (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "%s", SCIPconsGetName(cons) );
         SCIPinfoMessage(scip, file, "QCMATRIX %s\n", namestr);

         /* print alpha_i^2 x_i^2 terms */
         for( v = 0; v < nconsvars; ++v )
         {
            if( coefs[v] == 0.0 )
               continue;

            /* get variable name */
            assert ( SCIPhashmapExists(varnameHashmap, consvars[v]) );
            varname = (const char*) SCIPhashmapGetImage(varnameHashmap, consvars[v]);

            /* get coefficient^2 as string */
            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", coefs[v]*coefs[v]);

            /* print "x x coeff" line */
            printStart(scip, file, "", varname, (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n", valuestr);
         }

         /* print -(alpha_{n+1} x_{n+1})^2 term */

         /* get variable name */
         var = SCIPgetRhsVarSOC(scip, cons);
         assert ( SCIPhashmapExists(varnameHashmap, var) );
         varname = (const char*) SCIPhashmapGetImage(varnameHashmap, var);

         /* get -coefficient^2 as string */
         (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25.15g", -SCIPgetRhsCoefSOC(scip, cons)*SCIPgetRhsCoefSOC(scip, cons));

         /* print "x x coeff" line */
         printStart(scip, file, "", varname, (int) maxnamelen);
         printRecord(scip, file, varname, valuestr, maxnamelen);
         SCIPinfoMessage(scip, file, "\n", valuestr);
      }

      SCIPfreeBufferArray(scip, &namestr);
   }

   /* print indicator section */
   if( nConsIndicator > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &namestr, MPS_MAX_NAMELEN) );

      SCIPinfoMessage(scip, file, "INDICATORS\n");
      SCIPdebugMsg(scip, "start printing INDICATOR section\n");

      /* output each indicator constraint */
      for( c = 0; c < nConsIndicator; ++c )
      {
         SCIP_CONS* lincons;
         SCIP_VAR* slackvar;
         SCIP_VAR* binvar;

         cons = consIndicator[c];
         binvar = SCIPgetBinaryVarIndicator(cons);
         lincons = SCIPgetLinearConsIndicator(cons);
         slackvar = SCIPgetSlackVarIndicator(cons);

         /* create variable and value strings */
         if( SCIPvarIsNegated(binvar) )
         {
            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25d", 0);
            assert( SCIPvarGetNegatedVar(binvar) != NULL );
            assert( SCIPhashmapExists(varnameHashmap, SCIPvarGetNegatedVar(binvar)) );
            varname = (const char*) SCIPhashmapGetImage(varnameHashmap, SCIPvarGetNegatedVar(binvar));
         }
         else
         {
            (void) SCIPsnprintf(valuestr, MPS_MAX_VALUELEN, "%25d", 1);
            assert ( SCIPhashmapExists(varnameHashmap, binvar) );
            varname = (const char*) SCIPhashmapGetImage(varnameHashmap, binvar);
         }

         /* write records */
         if ( SCIPvarGetStatus(slackvar) == SCIP_VARSTATUS_AGGREGATED )
         {
            /* for aggregated variables output name of aggregating constraint */
            (void) SCIPsnprintf(namestr, MPS_MAX_NAMELEN, "aggr_%s", SCIPvarGetName(slackvar));
            printStart(scip, file, "IF", namestr, (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n");
         }
         else
         {
            printStart(scip, file, "IF", SCIPconsGetName(lincons), (int) maxnamelen);
            printRecord(scip, file, varname, valuestr, maxnamelen);
            SCIPinfoMessage(scip, file, "\n");
         }
      }
      SCIPfreeBufferArray(scip, &namestr);
   }

   /* free matrix data structure */
   freeMatrix(scip, matrix);

   /* free slackvar hashtable */
   if( indicatorSlackHash != NULL )
      SCIPhashtableFree(&indicatorSlackHash);

   /* free variable hashmap */
   SCIPhashmapFree(&varnameHashmap);

   SCIPfreeBlockMemoryArray(scip, &aggvars, saggvars);
   SCIPfreeBufferArray(scip, &rhss);

   /* free buffer arrays for SOS1, SOS2, and quadratic */
   SCIPfreeBufferArray(scip, &consIndicator);
   SCIPfreeBufferArray(scip, &consSOC);
   SCIPfreeBufferArray(scip, &consQuadratic);
   SCIPfreeBufferArray(scip, &consSOS2);
   SCIPfreeBufferArray(scip, &consSOS1);

   /* free variable and constraint name array */
   for( v = nvars + naggvars + nfixvars - 1; v >= 0; --v )
      SCIPfreeBufferArray(scip, &varnames[v]);
   SCIPfreeBufferArray(scip, &varnames);

   for( c = nconss + naddrows + naggvars - 1; c >= 0; --c )
      SCIPfreeBufferArray(scip, &consnames[c]);
   SCIPfreeBufferArray(scip, &consnames);

   /* print end of data line */
   SCIPinfoMessage(scip, file, "ENDATA");

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * mps file reader specific interface methods
 */

/** includes the mps file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderMps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyMps) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeMps) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadMps) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteMps) );

   /* add lp-reader parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/linearize-and-constraints",
         "should possible \"and\" constraint be linearized when writing the mps file?",
         &readerdata->linearizeands, TRUE, DEFAULT_LINEARIZE_ANDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/" READER_NAME "/aggrlinearization-ands",
         "should an aggregated linearization for and constraints be used?",
         &readerdata->aggrlinearizationands, TRUE, DEFAULT_AGGRLINEARIZATION_ANDS, NULL, NULL) );

   return SCIP_OKAY;
}
