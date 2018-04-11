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

/**@file   objrelax.cpp
 * @brief  C++ wrapper for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objrelax.h"




/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   scip::ObjRelax*       objrelax;           /**< relaxator object */
   SCIP_Bool             deleteobject;       /**< should the relaxator object be deleted when relaxator is freed? */
};




/*
 * Callback methods of relaxator
 */

extern "C"
{

/** copy method for relaxator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopyObj)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;
   
   assert(scip != NULL);
   
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);
   assert(relaxdata->objrelax->scip_ != scip);

   if( relaxdata->objrelax->iscloneable() )
   {
      scip::ObjRelax* newobjrelax;
      newobjrelax = dynamic_cast<scip::ObjRelax*> (relaxdata->objrelax->clone(scip));

      /* call include method of relaxator object */
      SCIP_CALL( SCIPincludeObjRelax(scip, newobjrelax, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeObj)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);
   assert(relaxdata->objrelax->scip_ == scip);

   /* call virtual method of relax object */
   SCIP_CALL( relaxdata->objrelax->scip_free(scip, relax) );

   /* free relax object */
   if( relaxdata->deleteobject )
      delete relaxdata->objrelax;

   /* free relax data */
   delete relaxdata;
   SCIPrelaxSetData(relax, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of relaxator (called after problem was transformed) */
static
SCIP_DECL_RELAXINIT(relaxInitObj)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);
   assert(relaxdata->objrelax->scip_ == scip);

   /* call virtual method of relax object */
   SCIP_CALL( relaxdata->objrelax->scip_init(scip, relax) );

   return SCIP_OKAY;
}


/** deinitialization method of relaxator (called before transformed problem is freed) */
static
SCIP_DECL_RELAXEXIT(relaxExitObj)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   SCIP_CALL( relaxdata->objrelax->scip_exit(scip, relax) );

   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolObj)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   SCIP_CALL( relaxdata->objrelax->scip_initsol(scip, relax) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolObj)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   SCIP_CALL( relaxdata->objrelax->scip_exitsol(scip, relax) );

   return SCIP_OKAY;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecObj)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   assert(relaxdata->objrelax != NULL);

   /* call virtual method of relax object */
   SCIP_CALL( relaxdata->objrelax->scip_exec(scip, relax, lowerbound, result) );

   return SCIP_OKAY;
}
}



/*
 * relaxator specific interface methods
 */

/** creates the relaxator for the given relaxator object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjRelax*       objrelax,           /**< relaxator object */
   SCIP_Bool             deleteobject        /**< should the relaxator object be deleted when relaxator is freed? */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(objrelax != NULL);

   /* create relaxator data */
   relaxdata = new SCIP_RELAXDATA;
   relaxdata->objrelax = objrelax;
   relaxdata->deleteobject = deleteobject;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelax(scip, objrelax->scip_name_, objrelax->scip_desc_,
         objrelax->scip_priority_, objrelax->scip_freq_, relaxCopyObj,
         relaxFreeObj, relaxInitObj, relaxExitObj,
         relaxInitsolObj, relaxExitsolObj, relaxExecObj,
         relaxdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the relax object of the given name, or 0 if not existing */
scip::ObjRelax* SCIPfindObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxator */
   )
{
   SCIP_RELAX* relax;
   SCIP_RELAXDATA* relaxdata;

   relax = SCIPfindRelax(scip, name);
   if( relax == NULL )
      return 0;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->objrelax;
}
   
/** returns the relax object for the given relaxator */
scip::ObjRelax* SCIPgetObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->objrelax;
}
