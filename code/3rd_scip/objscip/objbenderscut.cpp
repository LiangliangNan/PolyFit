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

/**@file   objbenderscut.cpp
 * @brief  C++ wrapper for the Benders' decomposition cut plugins
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objbenderscut.h"




/*
 * Data structures
 */

/** Benders' decomposition data */
struct SCIP_BenderscutData
{
   scip::ObjBenderscut*  objbenderscut;      /**< the Benders' decomposition object */
   SCIP_Bool             deleteobject;       /**< should the Benders' decomposition object be deleted when benders is freed? */
};




/*
 * Callback methods of the Benders' decomposition framework
 */

extern "C"
{

/** copy method for Benders' decomposition cuts (called when SCIP copies plugins) */
static
SCIP_DECL_BENDERSCUTCOPY(benderscutCopyObj)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert(scip != NULL);

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->objbenderscut != NULL);
   assert(benderscutdata->objbenderscut->scip_ != scip);

   if( benderscutdata->objbenderscut->iscloneable() )
   {
      scip::ObjBenderscut* newobjbenderscut;
      newobjbenderscut = dynamic_cast<scip::ObjBenderscut*> (benderscutdata->objbenderscut->clone(scip));

      /* call include method of Benders' decomposition object */
      SCIP_CALL( SCIPincludeObjBenderscut(scip, SCIPgetObjBenders(scip, benders), newobjbenderscut, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of Benders' decomposition cut to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BENDERSCUTFREE(benderscutFreeObj)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->objbenderscut != NULL);
   assert(benderscutdata->objbenderscut->scip_ == scip);

   /* call virtual method of benderscut object */
   SCIP_CALL( benderscutdata->objbenderscut->scip_free(scip, benderscut) );

   /* free benderscut object */
   if( benderscutdata->deleteobject )
      delete benderscutdata->objbenderscut;

   /* free benderscut data */
   delete benderscutdata;
   SCIPbenderscutSetData(benderscut, NULL); /*lint !e64*/

   return SCIP_OKAY;
}


/** initialization method of Benders' decomposition cut (called after problem was transformed) */
static
SCIP_DECL_BENDERSCUTINIT(benderscutInitObj)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->objbenderscut != NULL);
   assert(benderscutdata->objbenderscut->scip_ == scip);

   /* call virtual method of benderscut object */
   SCIP_CALL( benderscutdata->objbenderscut->scip_init(scip, benderscut) );

   return SCIP_OKAY;
}


/** deinitialization method of Benders' decomposition cut (called before transformed problem is freed) */
static
SCIP_DECL_BENDERSCUTEXIT(benderscutExitObj)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->objbenderscut != NULL);

   /* call virtual method of benderscut object */
   SCIP_CALL( benderscutdata->objbenderscut->scip_exit(scip, benderscut) );

   return SCIP_OKAY;
}


/** solving process initialization method of Benders' decomposition cut (called when branch and bound process is about to begin) */
static
SCIP_DECL_BENDERSCUTINITSOL(benderscutInitsolObj)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->objbenderscut != NULL);

   /* call virtual method of benderscut object */
   SCIP_CALL( benderscutdata->objbenderscut->scip_initsol(scip, benderscut) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of Benders' decomposition cut (called before branch and bound process data is freed) */
static
SCIP_DECL_BENDERSCUTEXITSOL(benderscutExitsolObj)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->objbenderscut != NULL);

   /* call virtual method of benderscut object */
   SCIP_CALL( benderscutdata->objbenderscut->scip_exitsol(scip, benderscut) );

   return SCIP_OKAY;
}


/** execution method of the Benders' decomposition cuts */
static
SCIP_DECL_BENDERSCUTEXEC(benderscutExecObj)
{  /*lint --e{715}*/
   SCIP_BENDERSCUTDATA* benderscutdata;

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);
   assert(benderscutdata->objbenderscut != NULL);

   /* call virtual method of benderscut object */
   SCIP_CALL( benderscutdata->objbenderscut->scip_exec(scip, benders, benderscut, sol, probnumber, type, result) );

   return SCIP_OKAY;
}


}


/*
 * Benders' decomposition cut specific interface methods
 */

/** creates the Benders' decomposition cut for the given Benders' decomposition cut object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjBenderscut(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjBenders*     objbenders,         /**< Benders' decomposition object */
   scip::ObjBenderscut*  objbenderscut,      /**< Benders' decomposition cut object */
   SCIP_Bool             deleteobject        /**< should the Benders' decomposition cut object be deleted when benderscut is freed? */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert(scip != NULL);
   assert(objbenderscut != NULL);
   assert(objbenderscut->scip_ == scip);

   /* create obj Benderscut' decomposition data */
   benderscutdata = new SCIP_BENDERSCUTDATA;
   benderscutdata->objbenderscut = objbenderscut;
   benderscutdata->deleteobject = deleteobject;

   benders = SCIPfindBenders(scip, objbenders->scip_name_);
   assert(benders != NULL);

   /* include Benderscut' decomposition */
   SCIP_CALL( SCIPincludeBenderscut(scip, benders, objbenderscut->scip_name_, objbenderscut->scip_desc_,
         objbenderscut->scip_priority_, objbenderscut->scip_islpcut_, benderscutCopyObj, benderscutFreeObj,
         benderscutInitObj, benderscutExitObj, benderscutInitsolObj, benderscutExitsolObj, benderscutExecObj,
         benderscutdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the benderscut object of the given name, or 0 if not existing */
scip::ObjBenderscut* SCIPfindObjBenderscut(
   scip::ObjBenders*     objbenders,         /**< Benders' decomposition object */
   const char*           name                /**< name of Benderscut' decomposition */
   )
{
   SCIP_BENDERS* benders;
   SCIP_BENDERSCUT* benderscut;
   SCIP_BENDERSCUTDATA* benderscutdata;

   benders = SCIPfindBenders(objbenders->scip_, objbenders->scip_name_);
   assert(benders != NULL);

   benderscut = SCIPfindBenderscut(benders, name);
   if( benderscut == NULL )
      return 0;

   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);

   return benderscutdata->objbenderscut;
}

/** returns the benderscut object for the given Benderscut' decomposition */
scip::ObjBenderscut* SCIPgetObjBenderscut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERSCUT*      benderscut          /**< Benderscut' decomposition */
   )
{
   SCIP_BENDERSCUTDATA* benderscutdata;

   assert(scip != NULL);
   benderscutdata = SCIPbenderscutGetData(benderscut);
   assert(benderscutdata != NULL);

   return benderscutdata->objbenderscut;
}
