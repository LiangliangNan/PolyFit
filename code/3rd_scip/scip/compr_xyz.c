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

/**@file   compr_xyz.c
 * @ingroup OTHER_CFILES
 * @brief  xyz tree compression
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/compr_xyz.h"


#define COMPR_NAME             "xyz"
#define COMPR_DESC             "tree compression template"
#define COMPR_PRIORITY         0
#define COMPR_MINNNODES        0


/*
 * Data structures
 */

/* TODO: fill in the necessary compression data */

/** tree compression data */
struct SCIP_ComprData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of tree compression
 */

/* TODO: Implement all necessary tree compression methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for tree compression plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_COMPRCOPY(comprCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz tree compression not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define comprCopyXyz NULL
#endif

/** destructor of tree compression to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_COMPRFREE(comprFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz tree compression not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define comprFreeXyz NULL
#endif


/** initialization method of tree compression (called after problem was transformed) */
#if 0
static
SCIP_DECL_COMPRINIT(comprInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz tree compression not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define comprInitXyz NULL
#endif


/** deinitialization method of tree compression (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_COMPREXIT(comprExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz tree compression not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define comprExitXyz NULL
#endif


/** solving process initialization method of tree compression (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_COMPRINITSOL(comprInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz tree compression not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define comprInitsolXyz NULL
#endif


/** solving process deinitialization method of tree compression (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_COMPREXITSOL(comprExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz tree compression not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define comprExitsolXyz NULL
#endif


/** execution method of tree compression */
static
SCIP_DECL_COMPREXEC(comprExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz tree compression not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/*
 * tree compression specific interface methods
 */

/** creates the xyz tree compression and includes it in SCIP */
SCIP_RETCODE SCIPincludeComprXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_COMPRDATA* comprdata;
   SCIP_COMPR* compr;

   /* create xyz tree compression data */
   comprdata = NULL;

   compr = NULL;

   /* include tree compression */
#if 0
   /* use SCIPincludeCompr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeCompr(scip, COMPR_NAME, COMPR_DESC, COMPR_PRIORITY, COMPR_MINNNODES,
         comprCopyXyz, comprFreeXyz, comprInitXyz, comprExitXyz, comprInitsolXyz, comprExitsolXyz, comprExecXyz,
         comprdata) );
#else
   /* use SCIPincludeComprBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeComprBasic(scip, &compr, COMPR_NAME, COMPR_DESC,
         COMPR_PRIORITY, COMPR_MINNNODES, comprExecXyz, comprdata) );

   assert(compr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetComprCopy(scip, compr, comprCopyXyz) );
   SCIP_CALL( SCIPsetComprFree(scip, compr, comprFreeXyz) );
   SCIP_CALL( SCIPsetComprInit(scip, compr, comprInitXyz) );
   SCIP_CALL( SCIPsetComprExit(scip, compr, comprExitXyz) );
   SCIP_CALL( SCIPsetComprInitsol(scip, compr, comprInitsolXyz) );
   SCIP_CALL( SCIPsetComprExitsol(scip, compr, comprExitsolXyz) );
#endif

   /* add xyz tree compression parameters */
   /* TODO: (optional) add tree compression specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
