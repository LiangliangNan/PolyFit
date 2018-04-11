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

/**@file   objprop.cpp
 * @brief  C++ wrapper for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objprop.h"




/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   scip::ObjProp*        objprop;            /**< propagator object */
   SCIP_Bool             deleteobject;       /**< should the propagator object be deleted when propagator is freed? */
};




/*
 * Callback methods of propagator
 */

extern "C"
{

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   
   assert(scip != NULL);
   
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);
   assert(propdata->objprop->scip_ != scip);

   if( propdata->objprop->iscloneable() )
   {
      scip::ObjProp* newobjprop;
      newobjprop = dynamic_cast<scip::ObjProp*> (propdata->objprop->clone(scip));

      /* call include method of propagator object */
      SCIP_CALL( SCIPincludeObjProp(scip, newobjprop, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);
   assert(propdata->objprop->scip_ == scip);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_free(scip, prop) );

   /* free prop object */
   if( propdata->deleteobject )
      delete propdata->objprop;

   /* free prop data */
   delete propdata;
   SCIPpropSetData(prop, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);
   assert(propdata->objprop->scip_ == scip);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_init(scip, prop) );

   return SCIP_OKAY;
}


/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_exit(scip, prop) );

   return SCIP_OKAY;
}


/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_initpre(scip, prop) );

   return SCIP_OKAY;
}


/** presolving deinitialization method of propagator (called after presolving has been finished) */
static
SCIP_DECL_PROPEXITPRE(propExitpreObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_exitpre(scip, prop) );

   return SCIP_OKAY;
}


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_initsol(scip, prop) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_exitsol(scip, prop, restart) );

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_presol(scip, prop, nrounds, presoltiming,
         nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
         nnewdelconss, nnewaddconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
         nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
         ndelconss, naddconss, nupgdconss, nchgcoefs, nchgsides, result) );

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_exec(scip, prop, proptiming, result) );

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropObj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->objprop != NULL);

   /* call virtual method of prop object */
   SCIP_CALL( propdata->objprop->scip_resprop(scip, prop, infervar, inferinfo, boundtype, bdchgidx, relaxedbd, result) );

   return SCIP_OKAY;
}
}



/*
 * propagator specific interface methods
 */

/** creates the propagator for the given propagator object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjProp*        objprop,            /**< propagator object */
   SCIP_Bool             deleteobject        /**< should the propagator object be deleted when propagator is freed? */
   )
{
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(objprop != NULL);

   /* create propagator data */
   propdata = new SCIP_PROPDATA;
   propdata->objprop = objprop;
   propdata->deleteobject = deleteobject;

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, objprop->scip_name_, objprop->scip_desc_,
         objprop->scip_priority_, objprop->scip_freq_, objprop->scip_delay_,
         objprop->scip_timingmask_, objprop->scip_presol_priority_, objprop->scip_presol_maxrounds_, objprop->scip_presol_timing_,
         propCopyObj, propFreeObj, propInitObj, propExitObj, propInitpreObj, propExitpreObj, propInitsolObj, propExitsolObj,
         propPresolObj, propExecObj, propRespropObj,
         propdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the prop object of the given name, or 0 if not existing */
scip::ObjProp* SCIPfindObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of propagator */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   prop = SCIPfindProp(scip, name);
   if( prop == NULL )
      return 0;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   return propdata->objprop;
}
   
/** returns the prop object for the given propagator */
scip::ObjProp* SCIPgetObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop                /**< propagator */
   )
{
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   return propdata->objprop;
}
