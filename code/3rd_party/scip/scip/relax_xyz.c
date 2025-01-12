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

/**@file   relax_xyz.c
 * @ingroup OTHER_CFILES
 * @brief  xyz relaxator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/relax_xyz.h"


#define RELAX_NAME             "xyz"
#define RELAX_DESC             "relaxator template"
#define RELAX_PRIORITY         0
#define RELAX_FREQ             1




/*
 * Data structures
 */

/* TODO: fill in the necessary relaxator data */

/** relaxator data */
struct SCIP_RelaxData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of relaxator
 */

/* TODO: Implement all necessary relaxator methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for relaxator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_RELAXCOPY(relaxCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxCopyXyz NULL
#endif

/** destructor of relaxator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_RELAXFREE(relaxFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxFreeXyz NULL
#endif


/** initialization method of relaxator (called after problem was transformed) */
#if 0
static
SCIP_DECL_RELAXINIT(relaxInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitXyz NULL
#endif


/** deinitialization method of relaxator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_RELAXEXIT(relaxExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitXyz NULL
#endif


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_RELAXINITSOL(relaxInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitsolXyz NULL
#endif


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitsolXyz NULL
#endif


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * relaxator specific interface methods
 */

/** creates the xyz relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create xyz relaxator data */
   relaxdata = NULL;
   /* TODO: (optional) create relaxator specific data here */

   relax = NULL;

   /* include relaxator */
#if 0
   /* use SCIPincludeRelax() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, RELAX_INCLUDESLP,
         relaxCopyXyz, relaxFreeXyz, relaxInitXyz, relaxExitXyz, relaxInitsolXyz, relaxExitsolXyz, relaxExecXyz,
         relaxdata) );
#else
   /* use SCIPincludeRelaxBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecXyz, relaxdata) );

   assert(relax != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopyXyz) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeXyz) );
   SCIP_CALL( SCIPsetRelaxInit(scip, relax, relaxInitXyz) );
   SCIP_CALL( SCIPsetRelaxExit(scip, relax, relaxExitXyz) );
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitsolXyz) );
   SCIP_CALL( SCIPsetRelaxExitsol(scip, relax, relaxExitsolXyz) );
#endif

   /* add xyz relaxator parameters */
   /* TODO: (optional) add relaxator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
