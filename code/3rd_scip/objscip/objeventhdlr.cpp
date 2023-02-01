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

/**@file   objeventhdlr.cpp
 * @brief  C++ wrapper for event handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objeventhdlr.h"




/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   scip::ObjEventhdlr*   objeventhdlr;       /**< event handler object */
   SCIP_Bool             deleteobject;       /**< should the event handler object be deleted when eventhdlristic is freed? */
};




/*
 * Callback methods of event handler
 */

extern "C"
{

/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventhdlrCopyObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   
   assert(scip != NULL);
   
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);
   assert(eventhdlrdata->objeventhdlr->scip_ != scip);

   if( eventhdlrdata->objeventhdlr->iscloneable() )
   {
      scip::ObjEventhdlr*  newobjeventhdlr;
      newobjeventhdlr = dynamic_cast<scip::ObjEventhdlr*> (eventhdlrdata->objeventhdlr->clone(scip));

      /* call include method of event handler object */
      SCIP_CALL( SCIPincludeObjEventhdlr(scip, newobjeventhdlr, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventhdlrFreeObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);
   assert(eventhdlrdata->objeventhdlr->scip_ == scip);

   /* call virtual method of eventhdlr object */
   SCIP_CALL( eventhdlrdata->objeventhdlr->scip_free(scip, eventhdlr) );

   /* free eventhdlr object */
   if( eventhdlrdata->deleteobject )
      delete eventhdlrdata->objeventhdlr;

   /* free eventhdlr data */
   delete eventhdlrdata;
   SCIPeventhdlrSetData(eventhdlr, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventhdlrInitObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);
   assert(eventhdlrdata->objeventhdlr->scip_ == scip);

   /* call virtual method of eventhdlr object */
   SCIP_CALL( eventhdlrdata->objeventhdlr->scip_init(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventhdlrExitObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   SCIP_CALL( eventhdlrdata->objeventhdlr->scip_exit(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventhdlrInitsolObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   SCIP_CALL( eventhdlrdata->objeventhdlr->scip_initsol(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventhdlrExitsolObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   SCIP_CALL( eventhdlrdata->objeventhdlr->scip_exitsol(scip, eventhdlr) );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_EVENTDELETE(eventhdlrDeleteObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   SCIP_CALL( eventhdlrdata->objeventhdlr->scip_delete(scip, eventhdlr, eventdata) );

   return SCIP_OKAY;
}


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventhdlrExecObj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->objeventhdlr != NULL);

   /* call virtual method of eventhdlr object */
   SCIP_CALL( eventhdlrdata->objeventhdlr->scip_exec(scip, eventhdlr, event, eventdata) );

   return SCIP_OKAY;
}
}



/*
 * event handler specific interface methods
 */

/** creates the event handler for the given event handler object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjEventhdlr*   objeventhdlr,       /**< event handler object */
   SCIP_Bool             deleteobject        /**< should the event handler object be deleted when eventhdlristic is freed? */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(objeventhdlr != NULL);

   /* create event handler data */
   eventhdlrdata = new SCIP_EVENTHDLRDATA;
   eventhdlrdata->objeventhdlr = objeventhdlr;
   eventhdlrdata->deleteobject = deleteobject;

   /* include event handler */
   SCIP_CALL( SCIPincludeEventhdlr(scip, objeventhdlr->scip_name_, objeventhdlr->scip_desc_,
         eventhdlrCopyObj,
         eventhdlrFreeObj, eventhdlrInitObj, eventhdlrExitObj,
         eventhdlrInitsolObj, eventhdlrExitsolObj, eventhdlrDeleteObj, eventhdlrExecObj,
         eventhdlrdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the eventhdlr object of the given name, or 0 if not existing */
scip::ObjEventhdlr* SCIPfindObjEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of event handler */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlr = SCIPfindEventhdlr(scip, name);
   if( eventhdlr == NULL )
      return 0;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return eventhdlrdata->objeventhdlr;
}
   
/** returns the eventhdlr object for the given event handler */
scip::ObjEventhdlr* SCIPgetObjEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return eventhdlrdata->objeventhdlr;
}
