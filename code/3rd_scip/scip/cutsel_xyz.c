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

/**@file   cutsel_xyz.c
 * @ingroup DEFPLUGINS_CUTSEL
 * @brief  xyz cut selector
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cutsel_xyz.h"


#define CUTSEL_NAME              "xyz"
#define CUTSEL_DESC              "cut selector template"
#define CUTSEL_PRIORITY                 0


/*
 * Data structures
 */

/* TODO: fill in the necessary cut selector data */

/** cut selector data */
struct SCIP_CutselData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of cut selector
 */

/* TODO: Implement all necessary cut selector methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for cut selector plugin (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CUTSELCOPY(cutselCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define cutselCopyXyz NULL
#endif


/** destructor of cut selector to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CUTSELFREE(cutselFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define cutselFreeXyz NULL
#endif


/** initialization method of cut selector (called after problem was transformed) */
#if 0
static
SCIP_DECL_CUTSELINIT(cutselInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define cutselInitXyz NULL
#endif


/** deinitialization method of cut selector (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CUTSELEXIT(cutselExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define cutselExitXyz NULL
#endif


/** solving process initialization method of cut selector (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CUTSELINITSOL(cutselInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define cutselInitsolXyz NULL
#endif


/** solving process deinitialization method of cut selector (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CUTSELEXITSOL(cutselExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define cutselExitsolXyz NULL
#endif


/** cut selection method of cut selector */
static
SCIP_DECL_CUTSELSELECT(cutselSelectXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz cut selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/*
 * cut selector specific interface methods
 */

/** creates the xyz cut selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeCutselXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CUTSELDATA* cutseldata = NULL;
   SCIP_CUTSEL* cutsel = NULL;

   /* create xyz cut selector data */

   /* TODO: (optional) create cut selector specific data here */

   /* include cut selector */
#if 0
   /* use SCIPincludeCutsel() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeCutsel(scip, CUTSEL_NAME, CUTSEL_DESC, CUTSEL_PRIORITY,
         cutselCopyXyz, cutselFreeXyz, cutselInitXyz, cutselExitXyz, cutselInitsolXyz, cutselExitsolXyz, cutselSelectXyz,
         cutseldata) );
#else
   /* use SCIPincludeCutselBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independently of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeCutselBasic(scip, &cutsel, CUTSEL_NAME, CUTSEL_DESC, CUTSEL_PRIORITY, cutselSelectXyz,
         cutseldata) );

   assert(cutsel != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetCutselCopy(scip, cutsel, cutselCopyXyz) );
   SCIP_CALL( SCIPsetCutselFree(scip, cutsel, cutselFreeXyz) );
   SCIP_CALL( SCIPsetCutselInit(scip, cutsel, cutselInitXyz) );
   SCIP_CALL( SCIPsetCutselExit(scip, cutsel, cutselExitXyz) );
   SCIP_CALL( SCIPsetCutselInitsol(scip, cutsel, cutselInitsolXyz) );
   SCIP_CALL( SCIPsetCutselExitsol(scip, cutsel, cutselExitsolXyz) );
#endif

   /* add xyz cut selector parameters */
   /* TODO: (optional) add cut selector specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
