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

/**@file   table.c
 * @brief  methods and datastructures for displaying statistics tables
 * @author Tristan Gally
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
#include "scip/table.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/syncstore.h"
#include "scip/struct_table.h"



/*
 * statistics table methods
 */

/** copies the given statistics table to a new scip */
SCIP_RETCODE SCIPtableCopyInclude(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(table != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( table->tablecopy != NULL )
   {
      SCIPsetDebugMsg(set, "including statistics table %s in subscip %p\n", SCIPtableGetName(table), (void*)set->scip);
      SCIP_CALL( table->tablecopy(set->scip, table) );
   }
   return SCIP_OKAY;
}

/** creates a statistics table */
SCIP_RETCODE SCIPtableCreate(
   SCIP_TABLE**          table,              /**< pointer to store statistics table */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of statistics table */
   const char*           desc,               /**< description of statistics table */
   SCIP_Bool             active,             /**< should the table be activated by default? */
   SCIP_DECL_TABLECOPY   ((*tablecopy)),     /**< copy method of statistics table or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_TABLEFREE   ((*tablefree)),     /**< destructor of statistics table */
   SCIP_DECL_TABLEINIT   ((*tableinit)),     /**< initialize statistics table */
   SCIP_DECL_TABLEEXIT   ((*tableexit)),     /**< deinitialize statistics table */
   SCIP_DECL_TABLEINITSOL ((*tableinitsol)), /**< solving process initialization method of statistics table */
   SCIP_DECL_TABLEEXITSOL ((*tableexitsol)), /**< solving process deinitialization method of statistics table */
   SCIP_DECL_TABLEOUTPUT ((*tableoutput)),   /**< output method */
   SCIP_TABLEDATA*       tabledata,          /**< display statistics table */
   int                   position,           /**< position of statistics table */
   SCIP_STAGE            earlieststage       /**< output of the statistics table is only printed from this stage onwards */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(table != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(tableoutput != NULL);

   SCIP_ALLOC( BMSallocMemory(table) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*table)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*table)->desc, desc, strlen(desc)+1) );
   (*table)->tablecopy = tablecopy;
   (*table)->tablefree = tablefree;
   (*table)->tableinit = tableinit;
   (*table)->tableexit = tableexit;
   (*table)->tableinitsol = tableinitsol;
   (*table)->tableexitsol = tableexitsol;
   (*table)->tableoutput = tableoutput;
   (*table)->tabledata = tabledata;
   (*table)->position = position;
   (*table)->earlieststage = earlieststage;
   (*table)->initialized = FALSE;
   (*table)->active = active;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "table/%s/active", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "is statistics table <%s> active", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, messagehdlr, blkmem, paramname, paramdesc,
         &(*table)->active, FALSE, active, NULL, NULL) );

   return SCIP_OKAY;
}

/** frees memory of statistics table */
SCIP_RETCODE SCIPtableFree(
   SCIP_TABLE**          table,              /**< pointer to statistics table data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(table != NULL);
   assert(*table != NULL);
   assert(!(*table)->initialized);
   assert(set != NULL);

   /* call destructor of statistics table */
   if( (*table)->tablefree != NULL )
   {
      SCIP_CALL( (*table)->tablefree(set->scip, *table) );
   }

   BMSfreeMemoryArray(&(*table)->name);
   BMSfreeMemoryArray(&(*table)->desc);
   BMSfreeMemory(table);

   return SCIP_OKAY;
}

/** initializes statistics table */
SCIP_RETCODE SCIPtableInit(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(table != NULL);
   assert(set != NULL);

   if( table->initialized )
   {
      SCIPerrorMessage("statistics table <%s> already initialized\n", table->name);
      return SCIP_INVALIDCALL;
   }

   if( table->tableinit != NULL )
   {
      SCIP_CALL( table->tableinit(set->scip, table) );
   }
   table->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes statistics table */
SCIP_RETCODE SCIPtableExit(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(table != NULL);
   assert(set != NULL);

   if( !table->initialized )
   {
      SCIPerrorMessage("statistics table <%s> not initialized\n", table->name);
      return SCIP_INVALIDCALL;
   }

   if( table->tableexit != NULL )
   {
      SCIP_CALL( table->tableexit(set->scip, table) );
   }
   table->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs statistics table that the branch and bound process is being started */
SCIP_RETCODE SCIPtableInitsol(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(table != NULL);
   assert(set != NULL);

   /* call solving process initialization method of statistics table */
   if( table->tableinitsol != NULL )
   {
      SCIP_CALL( table->tableinitsol(set->scip, table) );
   }

   return SCIP_OKAY;
}

/** informs statistics table that the branch and bound process data is being freed */
SCIP_RETCODE SCIPtableExitsol(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(table != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of statistics table */
   if( table->tableexitsol != NULL )
   {
      SCIP_CALL( table->tableexitsol(set->scip, table) );
   }

   return SCIP_OKAY;
}

/** output statistics table to screen */
SCIP_RETCODE SCIPtableOutput(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_SET*             set,                /**< global SCIP settings */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(table != NULL);
   assert(table->tableoutput != NULL);
   assert(set != NULL);

   SCIP_CALL( table->tableoutput(set->scip, table, file) );

   return SCIP_OKAY;
}

/** gets user data of statistics table */
SCIP_TABLEDATA* SCIPtableGetData(
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(table != NULL);

   return table->tabledata;
}

/** sets user data of statistics table; user has to free old data in advance! */
void SCIPtableSetData(
   SCIP_TABLE*           table,              /**< statistics table */
   SCIP_TABLEDATA*       tabledata           /**< new statistics table user data */
   )
{
   assert(table != NULL);

   table->tabledata = tabledata;
}

/** gets name of statistics table */
const char* SCIPtableGetName(
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(table != NULL);

   return table->name;
}

/** gets description of statistics table */
const char* SCIPtableGetDesc(
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(table != NULL);

   return table->desc;
}

/** gets position of statistics table */
int SCIPtableGetPosition(
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(table != NULL);

   return table->position;
}

/** gets earliest stage of statistics table */
SCIP_STAGE SCIPtableGetEarliestStage(
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(table != NULL);

   return table->earlieststage;
}

/** is statistics table currently active? */
SCIP_Bool SCIPtableIsActive(
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(table != NULL);

   return table->active;
}

/** is statistics table initialized? */
SCIP_Bool SCIPtableIsInitialized(
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(table != NULL);

   return table->initialized;
}
