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

/**@file   disp_xyz.c
 * @ingroup DEFPLUGINS_DISP
 * @brief  xyz display column
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/disp_xyz.h"


#define DISP_NAME               "xyz"                
#define DISP_DESC               "xyz display column"
#define DISP_HEADER             "xyz" 
#define DISP_WIDTH              14      /**< the width of the display column */
#define DISP_PRIORITY           110000  /**< the priority of the display column */
#define DISP_POSITION           30100   /**< the relative position of the display column */
#define DISP_STRIPLINE          TRUE    /**< the default for whether the display column should be separated 
                                         *   with a line from its right neighbor */




/*
 * Data structures
 */

/* TODO: fill in the necessary display column data */

/** display column data */
struct SCIP_DispData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of display column
 */

/* TODO: Implement all necessary display column methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for dialog plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_DISPCOPY(dispCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispCopyXyz NULL
#endif

/** destructor of display column to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_DISPFREE(dispFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispFreeXyz NULL
#endif


/** initialization method of display column (called after problem was transformed) */
#if 0
static
SCIP_DECL_DISPINIT(dispInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispInitXyz NULL
#endif


/** deinitialization method of display column (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_DISPEXIT(dispExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitXyz NULL
#endif


/** solving process initialization method of display column (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_DISPINITSOL(dispInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispInitsolXyz NULL
#endif


/** solving process deinitialization method of display column (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_DISPEXITSOL(dispExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitsolXyz NULL
#endif




/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * display column specific interface methods
 */

/** creates the xyz display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDispXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DISPDATA* dispdata;

   /* create xyz display column data */
   dispdata = NULL;
   /* TODO: (optional) create display column specific data here */

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, SCIP_DISPSTATUS_AUTO, 
         dispCopyXyz,
         dispFreeXyz, dispInitXyz, dispExitXyz, 
         dispInitsolXyz, dispExitsolXyz, dispOutputXyz, 
         dispdata, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

   /* add xyz display column parameters */
   /* TODO: (optional) add display column specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
