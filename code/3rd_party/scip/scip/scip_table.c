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

/**@file   scip_table.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for statistics table plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/scip_table.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/table.h"


/** creates a statistics table and includes it in SCIP */
SCIP_RETCODE SCIPincludeTable(
   SCIP*                 scip,               /**< SCIP data structure */
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
   SCIP_TABLEDATA*       tabledata,          /**< statistics table data */
   int                   position,           /**< position of statistics table */
   SCIP_STAGE            earlieststage       /**< output of the statistics table is only printed from this stage onwards */
   )
{
   SCIP_TABLE* table;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeTable", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether statistics table is already present */
   if( SCIPfindTable(scip, name) != NULL )
   {
      SCIPerrorMessage("statistics table <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtableCreate(&table, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, active, tablecopy,
         tablefree, tableinit, tableexit, tableinitsol, tableexitsol, tableoutput, tabledata,
         position, earlieststage) );
   SCIP_CALL( SCIPsetIncludeTable(scip->set, table) );

   return SCIP_OKAY;
}

/** returns the statistics table of the given name, or NULL if not existing */
SCIP_TABLE* SCIPfindTable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of statistics table */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindTable(scip->set, name);
}

/** returns the array of currently available statistics tables */
SCIP_TABLE** SCIPgetTables(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->tables;
}

/** returns the number of currently available statistics tables */
int SCIPgetNTables(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->ntables;
}
