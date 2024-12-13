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

/**@file   benders_xyz.c
 * @ingroup OTHER_CFILES
 * @brief  xyz Benders' decomposition algorithm
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/benders_xyz.h"
#include "scip/pub_benders.h"


#define BENDERS_NAME                "xyz"
#define BENDERS_DESC                "Benders' decomposition template"
#define BENDERS_PRIORITY            0
#define BENDERS_CUTLP            TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO        TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX         TRUE   /**< should Benders' cut be generated for relaxation solutions */
#define BENDERS_SHAREAUXVARS    FALSE   /**< should this Benders' share the highest priority Benders' aux vars */




/*
 * Data structures
 */

/* TODO: fill in the necessary Benders' decomposition data */

/** Benders' decomposition data */
struct SCIP_BendersData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods for Benders' decomposition
 */

/* TODO: Implement all necessary Benders' decomposition methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for benders plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCOPY(bendersCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decompostion not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   (*valid) = TRUE;

   return SCIP_OKAY;
}
#else
#define bendersCopyXyz NULL
#endif

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_BENDERSFREE(bendersFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersFreeXyz NULL
#endif


/** initialization method of Benders' decomposition (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSINIT(bendersInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitXyz NULL
#endif


/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSEXIT(bendersExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitXyz NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 */
#if 0
static
SCIP_DECL_BENDERSINITPRE(bendersInitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitpreXyz NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_BENDERSEXITPRE(bendersExitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitpreXyz NULL
#endif


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSINITSOL(bendersInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitsolXyz NULL
#endif


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitsolXyz NULL
#endif


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must create the SCIP instance for the subproblem and add the required variables and constraints. In
 *  addition, the settings required for the solving the problem must be set here. However, some settings will be
 *  overridden by the standard solving method included in the Benders' decomposition framework. If a special solving
 *  method is desired, the user can implement the bendersSolvesubXyz callback.
 *
 *  If the user defines a subproblem solving method in the bendersSolvesubXyz that doesn't involve a SCIP instance,
 *  then the user must explicitly specify the subproblem type. This is necessary because the dual solutions from convex
 *  problems can be used to generate cuts. The classical Benders' optimality and feasibility cuts require that the
 *  subproblems are convex. The subproblem type is specified by calling SCIPbendersSetSubproblemType. The available
 *  subproblem types are defined in SCIP_BENDERSSUBTYPE.
 */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** called before the subproblem solve for Benders' decomposition */
#if 0
static
SCIP_DECL_BENDERSPRESUBSOLVE(bendersPresubsolveXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersPresubsolveXyz NULL
#endif

/** the solving method for a convex subproblem for Benders' decomposition. In this method the subproblem is setup with
 *  the given solution and then solved.
 *  NOTE: if either bendersSolvesubconvexXyz or bendersSolvesubXyz callbacks are implemented then the bendersFreesubXyz
 *  callback must be implemented
 */
#if 0
static
SCIP_DECL_BENDERSSOLVESUBCONVEX(bendersSolvesubconvexXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersSolvesubconvexXyz NULL
#endif

/** the subproblem solving method for Benders' decomposition. In this method the subproblem is setup with the given
 *  solution and then solved.
 *  NOTE: if either bendersSolvesubconvexXyz or bendersSolvesubXyz callbacks are implemented then the bendersFreesubXyz
 *  callback must be implemented
 */
#if 0
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersSolvesubXyz NULL
#endif

#if 0
/** the post-solve method for Benders' decomposition */
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersPostsolveXyz NULL
#endif


/** the subproblem freeing method for Benders' decomposition. This is called between subproblem solves to clear the
 *  solving data. Generally this will only require a call to SCIPfreeTransform. However, depending on the problem it
 *  could additional freeing methods.
 *  NOTE: the bendersFreesubXyz callback must be implemented if the bendersSolvesubXyz is implemented */
#if 0
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersFreesubXyz NULL
#endif




/*
 * Benders' decomposition specific interface methods
 */

/** creates the xyz Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   /* create xyz Benders' decomposition data */
   bendersdata = NULL;
   /* TODO: (optional) create Benders' decomposition specific data here */

   benders = NULL;

   /* include Benders' decomposition */
#if 0
   /* use SCIPincludeBenders() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenders(scip, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, BENDERS_CUTLP, BENDERS_CUTPSEUDO,
         BENDERS_CUTRELAX, bendersCopyXyz, bendersFreeXyz, bendersInitXyz, bendersExitXyz, bendersInitpreXyz,
         bendersExitpreXyz, bendersInitsolXyz, bendersExitsolXyz, bendersGetvarXyz, bendersPresubsolveXyz,
         bendersCreatesubXyz, bendersSolvesubconvexXyz, bendersSolvesubXyz, bendersPostsolveXyz, bendersFreesubXyz,
         bendersdata) );
#else
   /* use SCIPincludeBendersBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, BENDERS_CUTLP,
         BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, BENDERS_SHAREAUXVARS, bendersGetvarXyz, bendersCreatesubXyz,
         bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyXyz) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeXyz) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitXyz) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitXyz) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreXyz) );
   SCIP_CALL( SCIPsetBendersExitpre(scip, benders, bendersExitpreXyz) );
   SCIP_CALL( SCIPsetBendersInitsol(scip, benders, bendersInitsolXyz) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolXyz) );
   SCIP_CALL( SCIPsetBendersPresubsolve(scip, benders, bendersPresubsolveXyz) );
   SCIP_CALL( SCIPsetBendersSolveAndFreesub(scip, benders, bendersSolvesubconvexXyz, bendersSolvesubXyz,
         bendersFreesubXyz) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveXyz) );
#endif

   /* OPTIONAL: including the default cuts for Benders' decomposition */
#if 0
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );
#endif

   /* add xyz Benders' decomposition parameters */
   /* TODO: (optional) add Benders' decomposition specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
