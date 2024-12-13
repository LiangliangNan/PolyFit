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

/**@file   objbenders.cpp
 * @brief  C++ wrapper for the Benders' decomposition plugins
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objbenders.h"




/*
 * Data structures
 */

/** Benders' decomposition data */
struct SCIP_BendersData
{
   scip::ObjBenders*     objbenders;         /**< the Benders' decomposition object */
   SCIP_Bool             deleteobject;       /**< should the Benders' decomposition object be deleted when benders is freed? */
};




/*
 * Callback methods of the Benders' decomposition framework
 */

extern "C"
{

/** copy method for Benders' decomposition plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BENDERSCOPY(bendersCopyObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);
   assert(bendersdata->objbenders->scip_ != scip);

   if( bendersdata->objbenders->iscloneable() )
   {
      scip::ObjBenders* newobjbenders;
      newobjbenders = dynamic_cast<scip::ObjBenders*> (bendersdata->objbenders->clone(scip));

      /* call include method of Benders' decomposition object */
      SCIP_CALL( SCIPincludeObjBenders(scip, newobjbenders, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BENDERSFREE(bendersFreeObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);
   assert(bendersdata->objbenders->scip_ == scip);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_free(scip, benders) );

   /* free benders object */
   if( bendersdata->deleteobject )
      delete bendersdata->objbenders;

   /* free benders data */
   delete bendersdata;
   SCIPbendersSetData(benders, NULL); /*lint !e64*/

   return SCIP_OKAY;
}


/** initialization method of Benders' decomposition (called after problem was transformed) */
static
SCIP_DECL_BENDERSINIT(bendersInitObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);
   assert(bendersdata->objbenders->scip_ == scip);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_init(scip, benders) );

   return SCIP_OKAY;
}


/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
static
SCIP_DECL_BENDERSEXIT(bendersExitObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_exit(scip, benders) );

   return SCIP_OKAY;
}


/** presolving initialization method of Benders' decomposition (called when presolving is about to begin) */
static
SCIP_DECL_BENDERSINITPRE(bendersInitpreObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_initpre(scip, benders) );

   return SCIP_OKAY;
}


/** presolving deinitialization method of Benders' decomposition (called after presolving has been finished) */
static
SCIP_DECL_BENDERSEXITPRE(bendersExitpreObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_exitpre(scip, benders) );

   return SCIP_OKAY;
}


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
static
SCIP_DECL_BENDERSINITSOL(bendersInitsolObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_initsol(scip, benders) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_exitsol(scip, benders) );

   return SCIP_OKAY;
}


/** method that is called to create the subproblem and register it with the Benders' decomposition structure. */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_createsub(scip, benders, probnumber) );

   return SCIP_OKAY;
}


/** methods called prior to solving the subproblems */
static
SCIP_DECL_BENDERSPRESUBSOLVE(bendersPresubsolveObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_presubsolve(scip, benders, sol, type, checkint, infeasible, auxviol,
         skipsolve, result) );

   return SCIP_OKAY;
}


/** method called to solve the convex relaxation of an individual subproblem of the Benders' decomposition */
static
SCIP_DECL_BENDERSSOLVESUBCONVEX(bendersSolvesubconvexObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_solvesubconvex(scip, benders, sol, probnumber, onlyconvexcheck, objective,
      result) );

   return SCIP_OKAY;
}


/** method called to solve an individual subproblem of the Benders' decomposition */
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_solvesub(scip, benders, sol, probnumber, objective, result) );

   return SCIP_OKAY;
}


/** method called after the subproblems are solved in the Benders' decomposition algorithm */
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_postsolve(scip, benders, sol, type, mergecands, npriomergecands,
      nmergecands, checkint, infeasible, merged) );

   return SCIP_OKAY;
}


/** frees an individual subproblem. Called in each iteration of the Benders' decomposition algorithm */
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_freesub(scip, benders, probnumber) );

   return SCIP_OKAY;
}


/** callback method to retrieve the master (subproblem) variable corresponding to the input subproblem (master) variable */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarObj)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);
   assert(bendersdata->objbenders != NULL);

   /* call virtual method of benders object */
   SCIP_CALL( bendersdata->objbenders->scip_getvar(scip, benders, var, mappedvar, probnumber) );

   return SCIP_OKAY;
}


}


/*
 * Benders' decomposition specific interface methods
 */

/** creates the Benders' decomposition for the given Benders' decomposition object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjBenders*     objbenders,         /**< Benders' decomposition object */
   SCIP_Bool             deleteobject        /**< should the Benders' decomposition object be deleted when benders is freed? */
   )
{
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   assert(objbenders != NULL);
   assert(objbenders->scip_ == scip);

   /* create obj Benders' decomposition data */
   bendersdata = new SCIP_BENDERSDATA;
   bendersdata->objbenders = objbenders;
   bendersdata->deleteobject = deleteobject;

   /* include Benders' decomposition */
   SCIP_CALL( SCIPincludeBenders(scip, objbenders->scip_name_, objbenders->scip_desc_,
      objbenders->scip_priority_, objbenders->scip_cutlp_, objbenders->scip_cutpseudo_,
      objbenders->scip_cutrelax_, objbenders->scip_shareauxvars_, bendersCopyObj, bendersFreeObj, bendersInitObj,
      bendersExitObj, bendersInitpreObj, bendersExitpreObj, bendersInitsolObj, bendersExitsolObj, bendersGetvarObj,
      bendersCreatesubObj, bendersPresubsolveObj, bendersSolvesubconvexObj, bendersSolvesubObj, bendersPostsolveObj,
      bendersFreesubObj, bendersdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the benders object of the given name, or 0 if not existing */
scip::ObjBenders* SCIPfindObjBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of Benders' decomposition */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSDATA* bendersdata;

   benders = SCIPfindBenders(scip, name);
   if( benders == NULL )
      return 0;

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   return bendersdata->objbenders;
}

/** returns the benders object for the given Benders' decomposition */
scip::ObjBenders* SCIPgetObjBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   )
{
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   return bendersdata->objbenders;
}
