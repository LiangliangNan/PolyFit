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

/**@file   objcutsel.cpp
 * @brief  C++ wrapper for cut selectors
 * @author Felipe Serrano
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objcutsel.h"




/*
 * Data structures
 */

/** cut selector data */
struct SCIP_CutselData
{
   scip::ObjCutsel*      objcutsel;          /**< cut selector object */
   SCIP_Bool             deleteobject;       /**< should the cut selector object be deleted when cut selector is freed? */
};




/*
 * Callback methods of cut selector
 */

extern "C"
{

/** copy method for cut selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CUTSELCOPY(cutselCopyObj)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   assert(scip != NULL);

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->objcutsel != NULL);
   assert(cutseldata->objcutsel->scip_ != scip);

   if( cutseldata->objcutsel->iscloneable() )
   {
      scip::ObjCutsel* newobjcutsel;
      newobjcutsel = dynamic_cast<scip::ObjCutsel*> (cutseldata->objcutsel->clone(scip));

      /* call include method of cut selector object */
      SCIP_CALL( SCIPincludeObjCutsel(scip, newobjcutsel, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of cut selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_CUTSELFREE(cutselFreeObj)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->objcutsel != NULL);
   assert(cutseldata->objcutsel->scip_ == scip);

   /* call virtual method of cutsel object */
   SCIP_CALL( cutseldata->objcutsel->scip_free(scip, cutsel) );

   /* free cutsel object */
   if( cutseldata->deleteobject )
      delete cutseldata->objcutsel;

   /* free cutsel data */
   delete cutseldata;
   SCIPcutselSetData(cutsel, NULL); /*lint !e64*/

   return SCIP_OKAY;
}


/** initialization method of cut selector (called after problem was transformed) */
static
SCIP_DECL_CUTSELINIT(cutselInitObj)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->objcutsel != NULL);
   assert(cutseldata->objcutsel->scip_ == scip);

   /* call virtual method of cutsel object */
   SCIP_CALL( cutseldata->objcutsel->scip_init(scip, cutsel) );

   return SCIP_OKAY;
}


/** deinitialization method of cut selector (called before transformed problem is freed) */
static
SCIP_DECL_CUTSELEXIT(cutselExitObj)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->objcutsel != NULL);

   /* call virtual method of cutsel object */
   SCIP_CALL( cutseldata->objcutsel->scip_exit(scip, cutsel) );

   return SCIP_OKAY;
}


/** solving process initialization method of cut selector (called when branch and bound process is about to begin) */
static
SCIP_DECL_CUTSELINITSOL(cutselInitsolObj)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->objcutsel != NULL);

   /* call virtual method of cutsel object */
   SCIP_CALL( cutseldata->objcutsel->scip_initsol(scip, cutsel) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of cut selector (called before branch and bound process data is freed) */
static
SCIP_DECL_CUTSELEXITSOL(cutselExitsolObj)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->objcutsel != NULL);

   /* call virtual method of cutsel object */
   SCIP_CALL( cutseldata->objcutsel->scip_exitsol(scip, cutsel) );

   return SCIP_OKAY;
}


/** cut selection method of cut selector */
static
SCIP_DECL_CUTSELSELECT(cutselSelectObj)
{  /*lint --e{715}*/
   SCIP_CUTSELDATA* cutseldata;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);
   assert(cutseldata->objcutsel != NULL);

   /* call virtual method of cutsel object */
   SCIP_CALL( cutseldata->objcutsel->scip_select(scip, cutsel, cuts, ncuts, forcedcuts, nforcedcuts,
      root, maxnselectedcuts, nselectedcuts, result) );

   return SCIP_OKAY;
}
}



/*
 * cut selector specific interface methods
 */

/** creates the cut selector for the given cut selector object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjCutsel(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjCutsel*      objcutsel,          /**< cut selector object */
   SCIP_Bool             deleteobject        /**< should the cut selector object be deleted when cut selector is freed? */
   )
{
   SCIP_CUTSELDATA* cutseldata;

   assert(scip != NULL);
   assert(objcutsel != NULL);

   /* create cut selector data */
   cutseldata = new SCIP_CUTSELDATA;
   cutseldata->objcutsel = objcutsel;
   cutseldata->deleteobject = deleteobject;

   /* include cut selector */
   SCIP_CALL( SCIPincludeCutsel(scip, objcutsel->scip_name_, objcutsel->scip_desc_,
         objcutsel->scip_priority_,
         cutselCopyObj,
         cutselFreeObj, cutselInitObj, cutselExitObj,
         cutselInitsolObj, cutselExitsolObj, cutselSelectObj,
         cutseldata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the cutsel object of the given name, or 0 if not existing */
scip::ObjCutsel* SCIPfindObjCutsel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of cut selector */
   )
{
   SCIP_CUTSEL* cutsel;
   SCIP_CUTSELDATA* cutseldata;

   cutsel = SCIPfindCutsel(scip, name);
   if( cutsel == NULL )
      return 0;

   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);

   return cutseldata->objcutsel;
}

/** returns the cutsel object for the given cut selector */
scip::ObjCutsel* SCIPgetObjCutsel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL*          cutsel              /**< cut selector */
   )
{
   SCIP_CUTSELDATA* cutseldata;

   assert(scip != NULL);
   cutseldata = SCIPcutselGetData(cutsel);
   assert(cutseldata != NULL);

   return cutseldata->objcutsel;
}
