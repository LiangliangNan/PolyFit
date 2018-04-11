/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objnodesel.cpp
 * @brief  C++ wrapper for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objnodesel.h"




/*
 * Data structures
 */

/** node selector data */
struct SCIP_NodeselData
{
   scip::ObjNodesel*     objnodesel;         /**< node selector object */
   SCIP_Bool             deleteobject;       /**< should the node selector object be deleted when node selector is freed? */
};




/*
 * Callback methods of node selector
 */

extern "C"
{

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;
   
   assert(scip != NULL);
   
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);
   assert(nodeseldata->objnodesel->scip_ != scip);

   if( nodeseldata->objnodesel->iscloneable() )
   {
      scip::ObjNodesel* newobjnodesel;
      newobjnodesel = dynamic_cast<scip::ObjNodesel*> (nodeseldata->objnodesel->clone(scip));

      /* call include method of node selector object */
      SCIP_CALL( SCIPincludeObjNodesel(scip, newobjnodesel, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_NODESELFREE(nodeselFreeObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);
   assert(nodeseldata->objnodesel->scip_ == scip);

   /* call virtual method of nodesel object */
   SCIP_CALL( nodeseldata->objnodesel->scip_free(scip, nodesel) );

   /* free nodesel object */
   if( nodeseldata->deleteobject )
      delete nodeseldata->objnodesel;

   /* free nodesel data */
   delete nodeseldata;
   SCIPnodeselSetData(nodesel, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of node selector (called after problem was transformed) */
static
SCIP_DECL_NODESELINIT(nodeselInitObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);
   assert(nodeseldata->objnodesel->scip_ == scip);

   /* call virtual method of nodesel object */
   SCIP_CALL( nodeseldata->objnodesel->scip_init(scip, nodesel) );

   return SCIP_OKAY;
}


/** deinitialization method of node selector (called before transformed problem is freed) */
static
SCIP_DECL_NODESELEXIT(nodeselExitObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   SCIP_CALL( nodeseldata->objnodesel->scip_exit(scip, nodesel) );

   return SCIP_OKAY;
}


/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
static
SCIP_DECL_NODESELINITSOL(nodeselInitsolObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   SCIP_CALL( nodeseldata->objnodesel->scip_initsol(scip, nodesel) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of node selector (called before branch and bound process data is freed) */
static
SCIP_DECL_NODESELEXITSOL(nodeselExitsolObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   SCIP_CALL( nodeseldata->objnodesel->scip_exitsol(scip, nodesel) );

   return SCIP_OKAY;
}


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   SCIP_CALL( nodeseldata->objnodesel->scip_select(scip, nodesel, selnode) );

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompObj)
{  /*lint --e{715}*/
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   assert(nodeseldata->objnodesel != NULL);

   /* call virtual method of nodesel object */
   return nodeseldata->objnodesel->scip_comp(scip, nodesel, node1, node2);
}
}



/*
 * node selector specific interface methods
 */

/** creates the node selector for the given node selector object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjNodesel*     objnodesel,         /**< node selector object */
   SCIP_Bool             deleteobject        /**< should the node selector object be deleted when node selector is freed? */
   )
{
   SCIP_NODESELDATA* nodeseldata;

   assert(scip != NULL);
   assert(objnodesel != NULL);

   /* create node selector data */
   nodeseldata = new SCIP_NODESELDATA;
   nodeseldata->objnodesel = objnodesel;
   nodeseldata->deleteobject = deleteobject;

   /* include node selector */
   SCIP_CALL( SCIPincludeNodesel(scip, objnodesel->scip_name_, objnodesel->scip_desc_,
         objnodesel->scip_stdpriority_, objnodesel->scip_memsavepriority_,
         nodeselCopyObj,
         nodeselFreeObj, nodeselInitObj, nodeselExitObj,
         nodeselInitsolObj, nodeselExitsolObj, nodeselSelectObj, nodeselCompObj,
         nodeseldata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the nodesel object of the given name, or 0 if not existing */
scip::ObjNodesel* SCIPfindObjNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of node selector */
   )
{
   SCIP_NODESEL* nodesel;
   SCIP_NODESELDATA* nodeseldata;

   nodesel = SCIPfindNodesel(scip, name);
   if( nodesel == NULL )
      return 0;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   return nodeseldata->objnodesel;
}
   
/** returns the nodesel object for the given node selector */
scip::ObjNodesel* SCIPgetObjNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   SCIP_NODESELDATA* nodeseldata;

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   return nodeseldata->objnodesel;
}
