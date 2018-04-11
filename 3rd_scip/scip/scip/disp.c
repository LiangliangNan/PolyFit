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

/**@file   disp.c
 * @brief  methods and datastructures for displaying runtime statistics
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/scip.h"
#include "scip/disp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/syncstore.h"
#include "scip/struct_disp.h"



/*
 * display column methods
 */

/** parameter change information method to autoselect display columns again */
SCIP_DECL_PARAMCHGD(SCIPparamChgdDispActive)
{  /*lint --e{715}*/
   /* automatically select the now active display columns */
   SCIP_CALL( SCIPautoselectDisps(scip) );

   return SCIP_OKAY;
}

/** copies the given display to a new scip */
SCIP_RETCODE SCIPdispCopyInclude(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(disp != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( disp->dispcopy != NULL )
   {
      SCIPsetDebugMsg(set, "including display column %s in subscip %p\n", SCIPdispGetName(disp), (void*)set->scip);
      SCIP_CALL( disp->dispcopy(set->scip, disp) );
   }
   return SCIP_OKAY;
}

/** creates a display column */
SCIP_RETCODE SCIPdispCreate(
   SCIP_DISP**           disp,               /**< pointer to store display column */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of display column */
   const char*           desc,               /**< description of display column */
   const char*           header,             /**< head line of display column */
   SCIP_DISPSTATUS       dispstatus,         /**< display activation status of display column */
   SCIP_DECL_DISPCOPY    ((*dispcopy)),      /**< copy method of display column or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   SCIP_DECL_DISPINIT    ((*dispinit)),      /**< initialize display column */
   SCIP_DECL_DISPEXIT    ((*dispexit)),      /**< deinitialize display column */
   SCIP_DECL_DISPINITSOL ((*dispinitsol)),   /**< solving process initialization method of display column */
   SCIP_DECL_DISPEXITSOL ((*dispexitsol)),   /**< solving process deinitialization method of display column */
   SCIP_DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   SCIP_DISPDATA*        dispdata,           /**< display column data */
   int                   width,              /**< width of display column (no. of chars used) */
   int                   priority,           /**< priority of display column */
   int                   position,           /**< relative position of display column */
   SCIP_Bool             stripline           /**< should the column be separated with a line from its right neighbor? */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(disp != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(header != NULL);
   assert(dispoutput != NULL);
   assert(width >= 0);

   SCIP_ALLOC( BMSallocMemory(disp) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*disp)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*disp)->desc, desc, strlen(desc)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*disp)->header, header, strlen(header)+1) );
   (*disp)->dispstatus = dispstatus;
   (*disp)->dispcopy = dispcopy;
   (*disp)->dispfree = dispfree;
   (*disp)->dispinit = dispinit;
   (*disp)->dispexit = dispexit;
   (*disp)->dispinitsol = dispinitsol;
   (*disp)->dispexitsol = dispexitsol;
   (*disp)->dispoutput = dispoutput;
   (*disp)->dispdata = dispdata;
   (*disp)->width = width;
   (*disp)->priority = priority;
   (*disp)->position = position;
   (*disp)->stripline = stripline;
   (*disp)->initialized = FALSE;
   (*disp)->active = (dispstatus == SCIP_DISPSTATUS_ON);
   (*disp)->mode = SCIP_DISPMODE_DEFAULT;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "display/%s/active", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "display activation status of display column <%s> (0: off, 1: auto, 2:on)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, messagehdlr, blkmem, paramname, paramdesc,
         (int*)(&(*disp)->dispstatus), FALSE, (int)dispstatus, 0, 2, SCIPparamChgdDispActive, NULL) );

   return SCIP_OKAY;
}

/** frees memory of display column */
SCIP_RETCODE SCIPdispFree(
   SCIP_DISP**           disp,               /**< pointer to display column data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(disp != NULL);
   assert(*disp != NULL);
   assert(!(*disp)->initialized);
   assert(set != NULL);

   /* call destructor of display column */
   if( (*disp)->dispfree != NULL )
   {
      SCIP_CALL( (*disp)->dispfree(set->scip, *disp) );
   }

   BMSfreeMemoryArray(&(*disp)->name);
   BMSfreeMemoryArray(&(*disp)->desc);
   BMSfreeMemoryArray(&(*disp)->header);
   BMSfreeMemory(disp);

   return SCIP_OKAY;
}

/** initializes display column */
SCIP_RETCODE SCIPdispInit(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(disp != NULL);
   assert(set != NULL);

   if( disp->initialized )
   {
      SCIPerrorMessage("display column <%s> already initialized\n", disp->name);
      return SCIP_INVALIDCALL;
   }

   if( disp->dispinit != NULL )
   {
      SCIP_CALL( disp->dispinit(set->scip, disp) );
   }
   disp->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes display column */
SCIP_RETCODE SCIPdispExit(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(disp != NULL);
   assert(set != NULL);

   if( !disp->initialized )
   {
      SCIPerrorMessage("display column <%s> not initialized\n", disp->name);
      return SCIP_INVALIDCALL;
   }

   if( disp->dispexit != NULL )
   {
      SCIP_CALL( disp->dispexit(set->scip, disp) );
   }
   disp->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs display column that the branch and bound process is being started */
SCIP_RETCODE SCIPdispInitsol(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(disp != NULL);
   assert(set != NULL);

   /* call solving process initialization method of display column */
   if( disp->dispinitsol != NULL )
   {
      SCIP_CALL( disp->dispinitsol(set->scip, disp) );
   }

   return SCIP_OKAY;
}

/** informs display column that the branch and bound process data is being freed */
SCIP_RETCODE SCIPdispExitsol(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(disp != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of display column */
   if( disp->dispexitsol != NULL )
   {
      SCIP_CALL( disp->dispexitsol(set->scip, disp) );
   }

   return SCIP_OKAY;
}

/** output display column to screen */
SCIP_RETCODE SCIPdispOutput(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_SET*             set,                /**< global SCIP settings */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(disp != NULL);
   assert(disp->dispoutput != NULL);
   assert(set != NULL);

   SCIP_CALL( disp->dispoutput(set->scip, disp, file) );

   return SCIP_OKAY;
}

/** gets user data of display column */
SCIP_DISPDATA* SCIPdispGetData(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->dispdata;
}

/** sets user data of display column; user has to free old data in advance! */
void SCIPdispSetData(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_DISPDATA*        dispdata            /**< new display column user data */
   )
{
   assert(disp != NULL);

   disp->dispdata = dispdata;
}

/** gets name of display column */
const char* SCIPdispGetName(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->name;
}

/** gets description of display column */
const char* SCIPdispGetDesc(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->desc;
}

/** gets head line of display column */
const char* SCIPdispGetHeader(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->header;
}

/** gets width of display column */
int SCIPdispGetWidth(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->width;
}

/** gets priority of display column */
int SCIPdispGetPriority(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->priority;
}

/** gets position of display column */
int SCIPdispGetPosition(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->position;
}

/** gets status of display column */
SCIP_DISPSTATUS SCIPdispGetStatus(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->dispstatus;
}

/** is display column initialized? */
SCIP_Bool SCIPdispIsInitialized(
   SCIP_DISP*            disp                /**< display column */
   )
{
   assert(disp != NULL);

   return disp->initialized;
}

/** prints one line of output with the active display columns */
SCIP_RETCODE SCIPdispPrintLine(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             forcedisplay,       /**< should the line be printed without regarding frequency? */
   SCIP_Bool             endline             /**< should the line be terminated with a newline symbol? */
   )
{
   assert(set != NULL);
   assert(set->disp_freq >= -1);
   assert(set->disp_headerfreq >= -1);
   assert(stat != NULL);

   if( (SCIP_VERBLEVEL)set->disp_verblevel < SCIP_VERBLEVEL_NORMAL || set->disp_freq == -1 )
      return SCIP_OKAY;

   if( forcedisplay
      || (stat->nnodes != stat->lastdispnode
         && set->disp_freq > 0
         && (stat->nnodes % set->disp_freq == 0 || stat->nnodes == 1)) )
   {
      int i;
      int j;
      SCIP_Bool stripline;

      /* display header line */
      if( (set->disp_headerfreq == 0 && stat->ndisplines == 0)
         || (set->disp_headerfreq > 0 && stat->ndisplines % set->disp_headerfreq == 0) )
      {
         int fillspace;

         stripline = FALSE;
         for( i = 0; i < set->ndisps; ++i )
         {
            assert(set->disps[i] != NULL);
            if( set->disps[i]->active )
            {
               if( stripline )
                  SCIPmessageFPrintInfo(messagehdlr, file, "|");
               fillspace = set->disps[i]->width - (int)strlen(set->disps[i]->header);
               for( j = 0; j < (fillspace)/2; ++j )
                  SCIPmessageFPrintInfo(messagehdlr, file, " ");
               SCIPmessageFPrintInfo(messagehdlr, file, "%s", (const char*)set->disps[i]->header);
               for( j = 0; j < (fillspace+1)/2; ++j )
                  SCIPmessageFPrintInfo(messagehdlr, file, " ");
               stripline = set->disps[i]->stripline;
            }
         }
         SCIPmessageFPrintInfo(messagehdlr, file, "\n");
      }

      /* display node information line */
      stripline = FALSE;
      for( i = 0; i < set->ndisps; ++i )
      {
         assert(set->disps[i] != NULL);
         if( set->disps[i]->active )
         {
            if( stripline )
               SCIPmessageFPrintInfo(messagehdlr, file, "|");
            SCIP_CALL( SCIPdispOutput(set->disps[i], set, file) );
            stripline = set->disps[i]->stripline;
         }
      }
      if( endline )
      {
         SCIPmessageFPrintInfo(messagehdlr, file, "\n");
      }
      fflush(stdout);

      stat->lastdispnode = stat->nnodes;
      stat->ndisplines++;
   }

   return SCIP_OKAY;
}

/** comparison method for display columns */
static
SCIP_DECL_SORTPTRCOMP(dispComp)
{  /*lint --e{715}*/
   return ((SCIP_DISP*)elem2)->priority - ((SCIP_DISP*)elem1)->priority;
}

/** activates all display lines fitting in the display w.r. to priority */
SCIP_RETCODE SCIPdispAutoActivate(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_DISP** disps;
   SCIP_SYNCSTORE* syncstore;
   SCIP_DISPMODE mode;
   int totalwidth;
   int width;
   int i;

   assert(set != NULL);

   syncstore = SCIPgetSyncstore(set->scip);
   assert(syncstore != NULL);

   /* sort display columns w.r. to their priority */
   SCIP_ALLOC( BMSduplicateMemoryArray(&disps, set->disps, set->ndisps) );
   SCIPsortPtr((void**)disps, dispComp, set->ndisps);

   totalwidth = 0;

   if( SCIPsyncstoreIsInitialized(syncstore) )
      mode = SCIP_DISPMODE_CONCURRENT;
   else
      mode = SCIP_DISPMODE_DEFAULT;

   /* first activate all columns with display status ON */
   for( i = 0; i < set->ndisps; ++i )
   {
      width = disps[i]->width;
      if( disps[i]->stripline )
         width++;
      if( disps[i]->dispstatus == SCIP_DISPSTATUS_ON && (disps[i]->mode & mode) )
      {
         disps[i]->active = TRUE;
         totalwidth += width;
      }
      else
         disps[i]->active = FALSE;
   }

   /* beginning with highest priority display column, activate AUTO columns as long as it fits into display width */
   for( i = 0; i < set->ndisps; ++i )
   {
      if( disps[i]->dispstatus == SCIP_DISPSTATUS_AUTO )
      {
         assert(!disps[i]->active);

         width = disps[i]->width;
         if( disps[i]->stripline )
            width++;
         if( totalwidth + width <= set->disp_width && (disps[i]->mode & mode) )
         {
            disps[i]->active = TRUE;
            totalwidth += width;
         }
      }
   }

   /* free temporary memory */
   BMSfreeMemoryArray(&disps);

   return SCIP_OKAY;
}

/** changes the display column mode */
void SCIPdispChgMode(
   SCIP_DISP*            disp,               /**< display column */
   SCIP_DISPMODE         mode                /**< the display column mode */
   )
{
   disp->mode = mode;
}

static
const char decpowerchar[] = {' ', 'k', 'M', 'G', 'T', 'P', 'E'};
#define MAXDECPOWER 6

/** displays a long integer in decimal form fitting in a given width */
void SCIPdispLongint(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   SCIP_Longint          val,                /**< value to display */
   int                   width               /**< width to fit into */
   )
{
   assert(width >= 1);

   if( width == 1 )
   {
      if( val < 0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "-");
      else if( val < 10 )
         SCIPmessageFPrintInfo(messagehdlr, file, "%" SCIP_LONGINT_FORMAT, val);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "+");
   }
   else
   {
      char format[SCIP_MAXSTRLEN];
      SCIP_Longint maxval;
      int decpower;
      int i;

      maxval = 1;
      for( i = 0; i < width-1; ++i )
         maxval *= 10;
      if( val < 0 )
         maxval /= 10;
      decpower = 0;
      while( ABS(val) >= maxval && decpower < MAXDECPOWER )
      {
         decpower++;
         val /= 1000;
      }
      (void) SCIPsnprintf(format, SCIP_MAXSTRLEN, "%%%d" SCIP_LONGINT_FORMAT "%c", width-1, decpowerchar[decpower]);

      if( width == 2 && val < 0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "-%c", decpowerchar[decpower]);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, (const char*)format, val);
   }
}

/** displays an integer in decimal form fitting in a given width */
void SCIPdispInt(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   int                   val,                /**< value to display */
   int                   width               /**< width to fit into */
   )
{
   SCIPdispLongint(messagehdlr, file, (SCIP_Longint)val, width);
}


static
const char timepowerchar[] = {'s', 'm', 'h', 'd', 'y'};
const SCIP_Real timepowerval[] = {1.0, 60.0, 60.0, 24.0, 365.0};
#define MAXTIMEPOWER 4

/** displays a time value fitting in a given width */
void SCIPdispTime(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< output stream */
   SCIP_Real             val,                /**< value in seconds to display */
   int                   width               /**< width to fit into */
   )
{
   assert(width >= 1);

   if( width == 1 )
   {
      if( val < 0.0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "-");
      else if( val < 10.0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "%.0f", val);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "+");
   }
   else
   {
      char format[SCIP_MAXSTRLEN];
      SCIP_Longint maxval;
      int timepower;
      int i;

      maxval = 1;
      for( i = 0; i < width-1; ++i )
         maxval *= 10;
      if( val < 0.0 )
         maxval /= 10;
      timepower = 0;
      while( REALABS(val) + 0.5 >= maxval && timepower < MAXTIMEPOWER )
      {
         timepower++;
         val /= timepowerval[timepower];
      }
      if( REALABS(val) + 0.05 < maxval/100.0 )
         (void) SCIPsnprintf(format, SCIP_MAXSTRLEN, "%%%d.1f%c", width-1, timepowerchar[timepower]);
      else
         (void) SCIPsnprintf(format, SCIP_MAXSTRLEN, "%%%d.0f%c", width-1, timepowerchar[timepower]);

      if( width == 2 && val < 0.0 )
         SCIPmessageFPrintInfo(messagehdlr, file, "-%c", timepowerchar[timepower]);
      else
         SCIPmessageFPrintInfo(messagehdlr, file, (const char*)format, val);
   }
}
