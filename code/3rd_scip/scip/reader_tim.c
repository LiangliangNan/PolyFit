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

/**@file   reader_tim.c
 * @ingroup DEFPLUGINS_READER
 * @brief  TIM file reader - the stage information for a stochastic programming instance in SMPS format
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_cons.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_reader.h"
#include "scip/reader_cor.h"
#include "scip/reader_tim.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_reader.h"
#include <string.h>

#define READER_NAME             "timreader"
#define READER_DESC             "file reader for the TIME file of a stochastic program in SMPS format"
#define READER_EXTENSION        "tim"

/*
 * tim reader internal methods
 */

#define TIM_MAX_LINELEN       1025
#define TIM_MAX_NAMELEN        256
#define TIM_DEFAULT_STAGESIZE   10
#define TIM_DEFAULT_ARRAYSIZE  100

#define BLANK         ' '

struct TimStage
{
   SCIP_VAR**            vars;
   SCIP_CONS**           conss;
   SCIP_HASHMAP*         varnametovar;
   SCIP_HASHMAP*         consnametocons;
   int                   nvars;
   int                   nconss;
   int                   varssize;
   int                   conssize;
};
typedef struct TimStage TIMSTAGE;

/** TIM reading data */
struct SCIP_ReaderData
{
   SCIP_Bool             read;               /**< flag to indicate whether the time file has been read */
   int                   nstages;            /**< the number of stages in the stochastic program */
   const char**          stagestartvars;     /**< the variables that start each stage */
   const char**          stagestartcons;     /**< the constraints that start each stage */
   const char**          stagenames;         /**< the name of the stage */
   TIMSTAGE**            stages;             /**< the stages for the stochastic program */
};

/** enum containing all tim sections */
enum TimSection
{
   TIM_TIME,
   TIM_PERIODS,
   TIM_ENDATA
};
typedef enum TimSection TIMSECTION;

/** tim input structure */
struct TimInput
{
   TIMSECTION            section;
   SCIP_FILE*            fp;
   int                   lineno;
   SCIP_Bool             haserror;
   char                  buf[TIM_MAX_LINELEN];
   const char*           f0;
   const char*           f1;
   const char*           f2;
   const char*           f3;
   char                  probname[TIM_MAX_NAMELEN];
   const char**          stagestartvars;
   const char**          stagestartcons;
   const char**          stagenames;
   int                   nstages;
   int                   stagesize;
};
typedef struct TimInput TIMINPUT;

/** adds the variable to the given stage */
static
SCIP_RETCODE addVariableToStage(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMSTAGE*             stage,              /**< the stage structure */
   const char*           varname             /**< the name of the variable to add to the stage */
   )
{
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(stage != NULL);

   var = SCIPfindVar(scip, varname);

   if( var == NULL )
   {
      SCIPwarningMessage(scip, "This is an error. All variables should in the problem.\n");
      return SCIP_OKAY;
   }

   /* adding the variable to the hashmap */
   SCIP_CALL( SCIPhashmapInsert(stage->varnametovar, (void*) varname, var) );

   /* adding the variable to the variable storage */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &stage->vars, &stage->varssize, stage->nvars + 1) );
   stage->vars[stage->nvars] = var;
   stage->nvars++;

   return SCIP_OKAY;
}

/** adds the constraint to the given stage */
static
SCIP_RETCODE addConstraintToStage(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMSTAGE*             stage,              /**< the stage structure */
   const char*           consname            /**< the name of the constraint to add to the stage */
   )
{
   SCIP_CONS* cons;

   assert(scip != NULL);
   assert(stage != NULL);

   cons = SCIPfindCons(scip, consname);

   if( cons == NULL )
   {
      SCIPwarningMessage(scip, "This is an error. All constraints should in the problem.\n");
      return SCIP_OKAY;
   }

   /* adding the constraint to the hashmap */
   SCIP_CALL( SCIPhashmapInsert(stage->consnametocons, (void*) consname, cons) );

   /* adding the constraint to the constraint storage */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &stage->conss, &stage->conssize, stage->nconss + 1) );
   stage->conss[stage->nconss] = cons;
   stage->nconss++;

   return SCIP_OKAY;
}

/** creates the stage data */
static
SCIP_RETCODE createStages(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the reader structure */
   SCIP_READER*          correader           /**< the reader structure for the core file */
   )
{
   SCIP_READERDATA* readerdata;
   int stage;
   int i;

   assert(scip != NULL);
   assert(reader != NULL);
   assert(correader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   stage = 0;

   /* assigning the variables to the stages */
   for( i = 0; i < SCIPcorGetNVarNames(correader); i++ )
   {
      /* the first variable in the var names list should be the start of the first stage */
      assert((stage == 0 && i == 0 && strcmp(SCIPcorGetVarName(correader, i), readerdata->stagestartvars[stage]) == 0)
         || i > 0);
      /* checking whether the next stage has been found */
      if( i > 0 && stage < readerdata->nstages - 1
         && strcmp(SCIPcorGetVarName(correader, i), readerdata->stagestartvars[stage + 1]) == 0 )
         stage++;

      /* adding the variable to the stage */
      SCIP_CALL( addVariableToStage(scip, readerdata->stages[stage], SCIPcorGetVarName(correader, i)) );
   }

   stage = 0;

   /* assigning the constraint to the stages */
   for( i = 0; i < SCIPcorGetNConsNames(correader); i++ )
   {
      /* checking whether the next stage has been found */
      if( i > 0 && stage < readerdata->nstages - 1
         && strcmp(SCIPcorGetConsName(correader, i), readerdata->stagestartcons[stage + 1]) == 0 )
         stage++;

      /* adding the consiable to the stage */
      SCIP_CALL( addConstraintToStage(scip, readerdata->stages[stage], SCIPcorGetConsName(correader, i)) );
   }

   return SCIP_OKAY;
}

/** creates the reader data for the time input data */
static
SCIP_RETCODE createReaderdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the reader structure */
   TIMINPUT*             timi                /**< tim input structure */
   )
{
   SCIP_READERDATA* readerdata;
   int hashmapsize;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(reader != NULL);
   assert(timi != NULL);

   readerdata = SCIPreaderGetData(reader);

   assert(readerdata != NULL);

   /* getting the total number of variables in the problem. The hash maps will be of size nvars/nstages. */
   nvars = SCIPgetNVars(scip);

   readerdata->read = TRUE;
   readerdata->nstages = timi->nstages;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->stagestartvars, readerdata->nstages) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->stagestartcons, readerdata->nstages) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->stagenames, readerdata->nstages) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->stages, readerdata->nstages) );

   for( i = 0; i < readerdata->nstages; i++ )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &readerdata->stagestartvars[i],
            timi->stagestartvars[i], strlen(timi->stagestartvars[i]) + 1) );  /*lint !e866*/
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &readerdata->stagestartcons[i],
            timi->stagestartcons[i], strlen(timi->stagestartcons[i]) + 1) );  /*lint !e866*/
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &readerdata->stagenames[i],
            timi->stagenames[i], strlen(timi->stagenames[i]) + 1) );          /*lint !e866*/

      /* creating the data for the stages */
      SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata->stages[i]) );        /*lint !e866*/
      readerdata->stages[i]->nvars = 0;
      readerdata->stages[i]->nconss = 0;
      readerdata->stages[i]->varssize = TIM_DEFAULT_ARRAYSIZE;
      readerdata->stages[i]->conssize = TIM_DEFAULT_ARRAYSIZE;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->stages[i]->vars, readerdata->stages[i]->varssize) );      /*lint !e866*/
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &readerdata->stages[i]->conss, readerdata->stages[i]->conssize) );     /*lint !e866*/

      /* creating the hashmaps */
      hashmapsize = (int) SCIPceil(scip, (SCIP_Real) nvars/(SCIP_Real) readerdata->nstages);
      SCIP_CALL( SCIPhashmapCreate(&readerdata->stages[i]->varnametovar, SCIPblkmem(scip), hashmapsize) );
      SCIP_CALL( SCIPhashmapCreate(&readerdata->stages[i]->consnametocons, SCIPblkmem(scip), hashmapsize) );
   }

   return SCIP_OKAY;
}

/** free the reader data */
static
void freeReaderdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader              /**< the reader structure */
   )
{
   SCIP_READERDATA* readerdata;
   int i;

   assert(scip != NULL);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);

   assert(readerdata != NULL);

   /* only free the reader data is a file has been read */
   if( readerdata->read )
   {
      for( i = 0; i < readerdata->nstages; i++ )
      {
         /* freeing the hashmaps */
         SCIPhashmapFree(&readerdata->stages[i]->consnametocons);
         SCIPhashmapFree(&readerdata->stages[i]->varnametovar);

         SCIPfreeBlockMemoryArray(scip, &readerdata->stagestartvars[i], strlen(readerdata->stagestartvars[i]) + 1);
         SCIPfreeBlockMemoryArray(scip, &readerdata->stagestartcons[i], strlen(readerdata->stagestartcons[i]) + 1);
         SCIPfreeBlockMemoryArray(scip, &readerdata->stagenames[i], strlen(readerdata->stagenames[i]) + 1);

         /* freeing the memory for the stage data */
         SCIPfreeBlockMemoryArray(scip, &readerdata->stages[i]->vars, readerdata->stages[i]->varssize);
         SCIPfreeBlockMemoryArray(scip, &readerdata->stages[i]->conss, readerdata->stages[i]->conssize);
         SCIPfreeBlockMemory(scip, &readerdata->stages[i]);    /*lint !e866*/
      }

      SCIPfreeBlockMemoryArray(scip, &readerdata->stages, readerdata->nstages);
      SCIPfreeBlockMemoryArray(scip, &readerdata->stagenames, readerdata->nstages);
      SCIPfreeBlockMemoryArray(scip, &readerdata->stagestartcons, readerdata->nstages);
      SCIPfreeBlockMemoryArray(scip, &readerdata->stagestartvars, readerdata->nstages);
   }

   SCIPfreeBlockMemory(scip, &readerdata);
}


/** creates the tim input structure */
static
SCIP_RETCODE timinputCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMINPUT**            timi,               /**< tim input structure */
   SCIP_FILE*            fp                  /**< file object for the input file */
   )
{
   assert(timi != NULL);
   assert(fp != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, timi) );

   (*timi)->section     = TIM_TIME;
   (*timi)->fp          = fp;
   (*timi)->lineno      = 0;
   (*timi)->haserror    = FALSE;
   (*timi)->buf     [0] = '\0';
   (*timi)->probname[0] = '\0';
   (*timi)->f0          = NULL;
   (*timi)->f1          = NULL;
   (*timi)->f2          = NULL;
   (*timi)->f3          = NULL;
   (*timi)->nstages     = 0;
   (*timi)->stagesize   = TIM_DEFAULT_STAGESIZE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*timi)->stagestartvars, (*timi)->stagesize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*timi)->stagestartcons, (*timi)->stagesize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*timi)->stagenames, (*timi)->stagesize) );

   return SCIP_OKAY;
}

/** free the tim input structure */
static
void timinputFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMINPUT**            timi                /**< tim input structure */
   )
{
   int i;

   for( i = 0; i < (*timi)->nstages; i++ )
   {
      SCIPfreeBlockMemoryArray(scip, &(*timi)->stagestartvars[i], strlen((*timi)->stagestartvars[i]) + 1);
      SCIPfreeBlockMemoryArray(scip, &(*timi)->stagestartcons[i], strlen((*timi)->stagestartcons[i]) + 1);
      SCIPfreeBlockMemoryArray(scip, &(*timi)->stagenames[i], strlen((*timi)->stagenames[i]) + 1);
   }

   SCIPfreeBlockMemoryArray(scip, &(*timi)->stagestartvars, (*timi)->stagesize);
   SCIPfreeBlockMemoryArray(scip, &(*timi)->stagestartcons, (*timi)->stagesize);
   SCIPfreeBlockMemoryArray(scip, &(*timi)->stagenames, (*timi)->stagesize);

   SCIPfreeBlockMemory(scip, timi);
}

/** returns the current section */
static
TIMSECTION timinputSection(
   const TIMINPUT*       timi                /**< tim input structure */
   )
{
   assert(timi != NULL);

   return timi->section;
}

/** return the current value of field 0 */
static
const char* timinputField0(
   const TIMINPUT*       timi                /**< tim input structure */
   )
{
   assert(timi != NULL);

   return timi->f0;
}

/** return the current value of field 1 */
static
const char* timinputField1(
   const TIMINPUT*       timi                /**< tim input structure */
   )
{
   assert(timi != NULL);

   return timi->f1;
}

/** return the current value of field 2 */
static
const char* timinputField2(
   const TIMINPUT*       timi                /**< tim input structure */
   )
{
   assert(timi != NULL);

   return timi->f2;
}

/** return the current value of field 3 */
static
const char* timinputField3(
   const TIMINPUT*       timi                /**< tim input structure */
   )
{
   assert(timi != NULL);

   return timi->f3;
}

/** returns if an error was detected */
static
SCIP_Bool timinputHasError(
   const TIMINPUT*       timi                /**< tim input structure */
   )
{
   assert(timi != NULL);

   return timi->haserror;
}

/** set the section in the tim input structure to given section */
static
void timinputSetSection(
   TIMINPUT*             timi,               /**< tim input structure */
   TIMSECTION            section             /**< section that is set */
   )
{
   assert(timi != NULL);

   timi->section = section;
}

/** set the problem name in the tim input structure to given problem name */
static
void timinputSetProbname(
   TIMINPUT*             timi,               /**< tim input structure */
   const char*           probname            /**< name of the problem to set */
   )
{
   assert(timi     != NULL);
   assert(probname != NULL);
   assert(strlen(probname) < sizeof(timi->probname));

   (void)SCIPmemccpy(timi->probname, probname, '\0', TIM_MAX_NAMELEN - 1);
}

/** set the problem var name that starts a stage in the tim input structure to given objective name */
static
SCIP_RETCODE timinputSetStageStartVar(
   TIMINPUT*             timi,               /**< tim input structure */
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           varname,            /**< name of the variable that starts the stage */
   int                   stagenum            /**< the stage number the variable starts */
   )
{
   assert(timi != NULL);
   assert(varname != NULL);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &timi->stagestartvars[stagenum], varname, strlen(varname) + 1) );  /*lint !e866*/

   return SCIP_OKAY;
}

/** set the problem constraint name that starts a stage in the tim input structure to given objective name */
static
SCIP_RETCODE timinputSetStageStartCons(
   TIMINPUT*             timi,               /**< tim input structure */
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           consname,           /**< name of the constraint that starts the stage */
   int                   stagenum            /**< the stage number the constraint starts */
   )
{
   assert(timi != NULL);
   assert(consname != NULL);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &timi->stagestartcons[stagenum], consname, strlen(consname) + 1) );   /*lint !e866*/

   return SCIP_OKAY;
}

/** set the stage name in the tim input structure to given objective name */
static
SCIP_RETCODE timinputSetStageName(
   TIMINPUT*             timi,               /**< tim input structure */
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           stagename,          /**< name of the stage */
   int                   stagenum            /**< the stage number the constraint starts */
   )
{
   assert(timi != NULL);
   assert(stagename != NULL);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &timi->stagenames[stagenum], stagename, strlen(stagename) + 1) );  /*lint !e866*/

   return SCIP_OKAY;
}

static
void timinputSyntaxerror(
   TIMINPUT*             timi                /**< tim input structure */
   )
{
   assert(timi != NULL);

   SCIPerrorMessage("Syntax error in line %d\n", timi->lineno);
   timi->section  = TIM_ENDATA;
   timi->haserror = TRUE;
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

/** read a tim format data line and parse the fields. */
static
SCIP_Bool timinputReadLine(
   TIMINPUT*             timi                /**< tim input structure */
   )
{
   unsigned int len;
   unsigned int i;
   char* s;
   SCIP_Bool is_empty;
   char* nexttok;

   do
   {
      timi->f0 = timi->f1 = timi->f2 = timi->f3 = 0;

      /* Read until we have not a comment line. */
      do
      {
         timi->buf[TIM_MAX_LINELEN-1] = '\0';
         if( NULL == SCIPfgets(timi->buf, (int) sizeof(timi->buf), timi->fp) )
            return FALSE;
         timi->lineno++;
      }
      while( *timi->buf == '*' );   /* coverity[a_loop_bound] */

      /* Normalize line */
      len = (unsigned int) strlen(timi->buf);

      for( i = 0; i < len; i++ )
         if( (timi->buf[i] == '\t') || (timi->buf[i] == '\n') || (timi->buf[i] == '\r') )
            timi->buf[i] = BLANK;

      if( len < 80 )
         clearFrom(timi->buf, len);

      SCIPdebugMessage("line %d: <%s>\n", timi->lineno, timi->buf);

      assert(strlen(timi->buf) >= 80);

      /* Look for new section */
      if( *timi->buf != BLANK )
      {
         timi->f0 = SCIPstrtok(&timi->buf[0], " ", &nexttok);

         assert(timi->f0 != 0);

         timi->f1 = SCIPstrtok(NULL, " ", &nexttok);

         return TRUE;
      }

      s = &timi->buf[1];

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
         if( NULL == (timi->f1 = SCIPstrtok(s, " ", &nexttok)) )
            break;

         if( (NULL == (timi->f2 = SCIPstrtok(NULL, " ", &nexttok))) || (*timi->f2 == '$') )
         {
            timi->f2 = 0;
            break;
         }

         if( (NULL == (timi->f3 = SCIPstrtok(NULL, " ", &nexttok))) || (*timi->f3 == '$') )
         {
            timi->f3 = 0;
            break;
         }
      }
      while( FALSE );

      /* check for empty lines */
      is_empty = (timi->f0 == NULL && timi->f1 == NULL);
   }
   while( is_empty );

   return TRUE;
}

/** Process TIME section. */
static
SCIP_RETCODE readTime(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMINPUT*             timi                /**< tim input structure */
   )
{
   SCIPdebugMsg(scip, "read problem name from TIME section\n");

   /* This has to be the Line with the TIME section. */
   if( !timinputReadLine(timi) || timinputField0(timi) == NULL || strcmp(timinputField0(timi), "TIME") )
   {
      timinputSyntaxerror(timi);
      return SCIP_OKAY;
   }

   /* Sometimes the name is omitted. */
   timinputSetProbname(timi, (timinputField1(timi) == 0) ? "_TIM_" : timinputField1(timi));

   /* This has to be a new section */
   /* coverity[tainted_data] */
   if( !timinputReadLine(timi) || (timinputField0(timi) == NULL) )
   {
      timinputSyntaxerror(timi);
      return SCIP_OKAY;
   }

   if( strncmp(timinputField0(timi), "PERIODS", 7) == 0 )
      timinputSetSection(timi, TIM_PERIODS);
   else
   {
      timinputSyntaxerror(timi);
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** Process PERIODS section. */
static
SCIP_RETCODE readPeriods(
   TIMINPUT*             timi,               /**< tim input structure */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIPdebugMsg(scip, "read Periods\n");

   /* coverity[tainted_data_sink_lv_call] */
   /* coverity[tainted_data] */
   while( timinputReadLine(timi) )
   {
      if( timinputField0(timi) != NULL )
      {
         if( strcmp(timinputField0(timi), "PERIODS") == 0 )
            timinputSetSection(timi, TIM_PERIODS);
         else if( strcmp(timinputField0(timi), "ENDATA") == 0 )
            timinputSetSection(timi, TIM_ENDATA);
         else
            timinputSyntaxerror(timi);

         return SCIP_OKAY;
      }

      if( timi->nstages + 1 >= timi->stagesize )
      {
         SCIP_CALL( SCIPensureBlockMemoryArray(scip, &timi->stagestartvars, &timi->stagesize, timi->nstages + 1) );
         SCIP_CALL( SCIPensureBlockMemoryArray(scip, &timi->stagestartcons, &timi->stagesize, timi->nstages + 1) );
         SCIP_CALL( SCIPensureBlockMemoryArray(scip, &timi->stagenames, &timi->stagesize, timi->nstages + 1) );
      }

      SCIP_CALL( timinputSetStageStartVar(timi, scip, timinputField1(timi), timi->nstages) );
      SCIP_CALL( timinputSetStageStartCons(timi, scip, timinputField2(timi), timi->nstages) );
      SCIP_CALL( timinputSetStageName(timi, scip, timinputField3(timi), timi->nstages) );

      timi->nstages++;
   }
   timinputSyntaxerror(timi);

   return SCIP_OKAY;
}


/** Read time data for the SMPS file format. */
static
SCIP_RETCODE readTim(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< the file reader itself */
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* fp;
   TIMINPUT* timi;
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

   SCIP_CALL_FINALLY( timinputCreate(scip, &timi, fp), SCIPfclose(fp) );

   SCIP_CALL_TERMINATE( retcode, readTime(scip, timi), TERMINATE );

   while( timinputSection(timi) == TIM_PERIODS )
   {
      /* coverity[tainted_data] */
      SCIP_CALL_TERMINATE( retcode, readPeriods(timi, scip), TERMINATE );
   }
   if( timinputSection(timi) != TIM_ENDATA )
      timinputSyntaxerror(timi);

   error = timinputHasError(timi);

   if( !error )
   {
      SCIP_CALL_TERMINATE( retcode, createReaderdata(scip, reader, timi), TERMINATE );
   }

 /* cppcheck-suppress unusedLabel */
 TERMINATE:
   timinputFree(scip, &timi);
   SCIPfclose(fp);

   if( error )
      return SCIP_READERROR;
   else
      return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyTim)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderTim(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeTim)
{
   freeReaderdata(scip, reader);

   return SCIP_OKAY;
}

/** reads the stage information for a stochastic programming instance in SMPS format */
static
SCIP_DECL_READERREAD(readerReadTim)
{  /*lint --e{715}*/
   SCIP_READER* correader;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   correader = SCIPfindReader(scip, "correader");

   if( correader == NULL )
   {
      SCIPwarningMessage(scip, "It is necessary to include the \"cor\" reader\n");
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* checking whether the cor file has been read */
   if( !SCIPcorHasRead(correader) )
   {
      SCIPwarningMessage(scip, "The core file must be read before the time and stochastic files.\n");
      (*result) = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPreadTim(scip, filename, result) );

   return SCIP_OKAY;
}


/*
 * tim file reader specific interface methods
 */

/** includes the tim file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderTim(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &readerdata) );
   readerdata->read = FALSE;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyTim) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeTim) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadTim) );

   return SCIP_OKAY;
}


/** reads problem from file */
SCIP_RETCODE SCIPreadTim(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{
   SCIP_READER* reader;
   SCIP_RETCODE retcode;
   SCIP_READERDATA* readerdata;

   assert(scip != NULL);
   assert(result != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   retcode = readTim(scip, reader, filename);

   if( retcode == SCIP_PLUGINNOTFOUND )
      retcode = SCIP_READERROR;

   if( retcode == SCIP_NOFILE || retcode == SCIP_READERROR )
      return retcode;

   SCIP_CALL( retcode );

   /* creating the stages */
   SCIP_CALL( createStages(scip, reader, SCIPfindReader(scip, "correader")) );

   /* setting the read flag to TRUE */
   readerdata = SCIPreaderGetData(reader);
   readerdata->read = TRUE;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/*
 * Interface methods for the cor and sto files
 */

/* return whether the tim file has been read */
SCIP_Bool SCIPtimHasRead(
   SCIP_READER*          reader              /**< the file reader itself */
   )
{
   SCIP_READERDATA* readerdata;

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->read;
}


/* returns the number of stages */
int SCIPtimGetNStages(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->nstages;
}

/* returns the name for a given stage */
const char* SCIPtimGetStageName(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(stagenum >= 0 && stagenum < readerdata->nstages);

   return readerdata->stagenames[stagenum];
}

/* returns the stage name for a given constraint name */
const char* SCIPtimConsGetStageName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           consname            /**< the constraint to search for */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   int stagenum;
   int i;
   int j;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   /* looping over all stages to find the provided constraint */
   stagenum = -1;
   for( i = 0; i < readerdata->nstages; i++ )
   {
      for( j = 0; j < readerdata->stages[i]->nconss; j++ )
      {
         if( strcmp(SCIPconsGetName(readerdata->stages[i]->conss[j]), consname) == 0 )
         {
            stagenum = i;
            break;
         }
      }

      if( stagenum >= 0 )
         break;
   }
   assert(stagenum >= 0 && stagenum < readerdata->nstages);

   return readerdata->stagenames[stagenum];
}

/* returns the number for a given stage */
int SCIPtimFindStage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           stage               /**< the name of the requested stage */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   int i;
   int stagenum;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   stagenum = -1;
   for( i = 0; i < readerdata->nstages; i++ )
   {
      if( strcmp(readerdata->stagenames[i], stage) == 0 )
      {
         stagenum = i;
         break;
      }
   }

   if( stagenum < 0 )
   {
      SCIPerrorMessage("Stage <%s> was not found in the TIM file. Check the SMPS files (COR, TIM and STO)\n", stage);
      SCIPABORT();
   }

   return stagenum;
}

/* returns the array of variables for a given stage */
SCIP_VAR** SCIPtimGetStageVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(stagenum >= 0 && stagenum < readerdata->nstages);

   return readerdata->stages[stagenum]->vars;
}

/* returns an array of constraints for a given stage */
SCIP_CONS** SCIPtimGetStageConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(stagenum >= 0 && stagenum < readerdata->nstages);

   return readerdata->stages[stagenum]->conss;
}

/* returns the number of variables for a given stage */
int SCIPtimGetStageNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(stagenum >= 0 && stagenum < readerdata->nstages);

   return readerdata->stages[stagenum]->nvars;
}

/* returns the number of constraints for a given stage */
int SCIPtimGetStageNConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   stagenum            /**< the number of the requested stage */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, READER_NAME);

   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(stagenum >= 0 && stagenum < readerdata->nstages);

   return readerdata->stages[stagenum]->nconss;
}
