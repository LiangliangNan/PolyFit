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

/**@file   prop_xyz.c
 * @ingroup DEFPLUGINS_PROP
 * @brief  xyz propagator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/prop_xyz.h"

/* fundamental propagator properties */
#define PROP_NAME              "xyz"
#define PROP_DESC              "propagator template"
#define PROP_PRIORITY                 0 /**< propagator priority */ 
#define PROP_FREQ                    10 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */

/* optional propagator properties
 * TODO: remove those properties if propagator does not support presolving */
#define PROP_PRESOL_PRIORITY          0 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_MEDIUM /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */


/*
 * Data structures
 */

/* TODO: fill in the necessary propagator data */

/** propagator data */
struct SCIP_PropData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of propagator
 */

/* TODO: Implement all necessary propagator methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for propagator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PROPCOPY(propCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propCopyXyz NULL
#endif

/** destructor of propagator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PROPFREE(propFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propFreeXyz NULL
#endif


/** initialization method of propagator (called after problem was transformed) */
#if 0
static
SCIP_DECL_PROPINIT(propInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propInitXyz NULL
#endif


/** deinitialization method of propagator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PROPEXIT(propExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propExitXyz NULL
#endif


/** presolving initialization method of propagator (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PROPINITPRE(propInitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propInitpreXyz NULL
#endif


/** presolving deinitialization method of propagator (called after presolving has been finished) */
#if 0
static
SCIP_DECL_PROPEXITPRE(propExitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propExitpreXyz NULL
#endif


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_PROPINITSOL(propInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propInitsolXyz NULL
#endif


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_PROPEXITSOL(propExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propExitsolXyz NULL
#endif


/** presolving method of propagator */
#if 0
static
SCIP_DECL_PROPPRESOL(propPresolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propPresolXyz NULL
#endif


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the xyz propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create xyz propagator data */
   propdata = NULL;
   /* TODO: (optional) create propagator specific data here */

   prop = NULL;

   /* include propagator */
#if 0
   /* use SCIPincludeProp() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY,
         PROP_TIMING, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING,
         propCopyXyz, propFreeXyz, propInitXyz, propExitXyz, propInitpreXyz, propExitpreXyz,
         propInitsolXyz, propExitsolXyz, propPresolXyz, propExecXyz, propRespropXyz,
         propdata) );
#else
   /* use SCIPincludePropBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecXyz, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyXyz) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeXyz) );
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitXyz) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitXyz) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolXyz) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolXyz) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreXyz) );
   SCIP_CALL( SCIPsetPropExitpre(scip, prop, propExitpreXyz) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolXyz, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropXyz) );
#endif

   /* add xyz propagator parameters */
   /* TODO: (optional) add propagator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
