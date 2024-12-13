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

/**@file   benderscut_xyz.c
 * @ingroup OTHER_CFILES
 * @brief  xyz Benders' decomposition cut
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/benderscut_xyz.h"


#define BENDERSCUT_NAME             "xyz"
#define BENDERSCUT_DESC             "Benders' decomposition cut template"
#define BENDERSCUT_PRIORITY         0
#define BENDERSCUT_LPCUT            TRUE


/*
 * Data structures
 */

/* TODO: fill in the necessary compression data */

/** Benders' decomposition cut data */
struct SCIP_BenderscutData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of the Benders' decomposition cut
 */

/* TODO: Implement all necessary Benders' decomposition cut methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for the Benders' decomposition cut plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCUTCOPY(benderscutCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition cut not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutCopyXyz NULL
#endif

/** destructor of the Benders' decomposition cut to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_BENDERSCUTFREE(benderscutFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition cut not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutFreeXyz NULL
#endif


/** initialization method of the Benders' decomposition cut (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSCUTINIT(benderscutInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition cut not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutInitXyz NULL
#endif


/** deinitialization method of the Benders' decomposition cut (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSCUTEXIT(benderscutExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition cut not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutExitXyz NULL
#endif


/** solving process initialization method of the Benders' decomposition cut (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSCUTINITSOL(benderscutInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition cut not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutInitsolXyz NULL
#endif


/** solving process deinitialization method of the Benders' decomposition cut (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSCUTEXITSOL(benderscutExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition cut not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define benderscutExitsolXyz NULL
#endif


/** execution method of the Benders' decomposition cut */
static
SCIP_DECL_BENDERSCUTEXEC(benderscutExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition cut not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/*
 * Benders' decomposition cut specific interface methods
 */

/** creates the xyz Benders' decomposition cut and includes it in SCIP */
SCIP_RETCODE SCIPincludeBenderscutXyz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;
   SCIP_BENDERSCUT* benderscut;

   assert(benders != NULL);

   /* create xyz Benders' decomposition cut data */
   benderscutdata = NULL;

   benderscut = NULL;

   /* include the Benders' decomposition cut */
#if 0
   /* use SCIPincludeBenderscut() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenderscut(scip, benders, BENDERSCUT_NAME, BENDERSCUT_DESC, BENDERSCUT_PRIORITY,
         BENDERSCUT_LPCUT, benderscutCopyXyz, benderscutFreeXyz, benderscutInitXyz, benderscutExitXyz,
         benderscutInitsolXyz, benderscutExitsolXyz, benderscutExecXyz, benderscutdata) );
#else
   /* use SCIPincludeBenderscutBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenderscutBasic(scip, benders, &benderscut, BENDERSCUT_NAME, BENDERSCUT_DESC,
         BENDERSCUT_PRIORITY, BENDERSCUT_LPCUT, benderscutExecXyz, benderscutdata) );

   assert(benderscut != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBenderscutCopy(scip, benderscut, benderscutCopyXyz) );
   SCIP_CALL( SCIPsetBenderscutFree(scip, benderscut, benderscutFreeXyz) );
   SCIP_CALL( SCIPsetBenderscutInit(scip, benderscut, benderscutInitXyz) );
   SCIP_CALL( SCIPsetBenderscutExit(scip, benderscut, benderscutExitXyz) );
   SCIP_CALL( SCIPsetBenderscutInitsol(scip, benderscut, benderscutInitsolXyz) );
   SCIP_CALL( SCIPsetBenderscutExitsol(scip, benderscut, benderscutExitsolXyz) );
#endif

   /* add xyz Benders' decomposition cut parameters */
   /* TODO: (optional) add the Benders' decomposition cut specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
