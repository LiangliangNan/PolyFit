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

/**@file   table_xyz.c
 * @ingroup OTHER_CFILES
 * @brief  xyz statistics table
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/table_xyz.h"


#define TABLE_NAME              "xyz"
#define TABLE_DESC              "xyz statistics table"
#define TABLE_POSITION          -1000           /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE    SCIP_STAGE_INIT /**< output of the statistics table is only printed from this stage onwards */




/*
 * Data structures
 */

/* TODO: fill in the necessary statistics table data */

/** statistics table data */
struct SCIP_TableData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of statistics table
 */

/* TODO: Implement all necessary statistics table methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for statistics table plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_TABLECOPY(tableCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableCopyXyz NULL
#endif


/** destructor of statistics table to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_TABLEFREE(tableFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableFreeXyz NULL
#endif


/** initialization method of statistics table (called after problem was transformed) */
#if 0
static
SCIP_DECL_TABLEINIT(tableInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableInitXyz NULL
#endif


/** deinitialization method of statistics table (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_TABLEEXIT(tableExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableExitXyz NULL
#endif


/** solving process initialization method of statistics table (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_TABLEINITSOL(tableInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableInitsolXyz NULL
#endif


/** solving process deinitialization method of statistics table (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_TABLEEXITSOL(tableExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableExitsolXyz NULL
#endif


/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * statistics table specific interface methods
 */

/** creates the xyz statistics table and includes it in SCIP */
SCIP_RETCODE SCIPincludeTableXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;

   /* create xyz statistics table data */
   tabledata = NULL;
   /* TODO: (optional) create statistics table specific data here */

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME, TABLE_DESC, TRUE,
         tableCopyXyz, tableFreeXyz, tableInitXyz, tableExitXyz,
         tableInitsolXyz, tableExitsolXyz, tableOutputXyz,
         tabledata, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   /* add xyz statistics table parameters */
   /* TODO: (optional) add statistics table specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
