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
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_globalbnd.c
 * @brief  eventhandler for storing all global bound changes
 * @author Robert Lion Gottwald
 *
 * the bound changes are stored so that they can be shared with other threads
 * in a concurrent solve.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_globalbnd.h"
#include "scip/pub_misc.h"
#include "scip/concurrent.h"
#include "scip/boundstore.h"
#include "scip/syncstore.h"
#include <string.h>

#define EVENTHDLR_NAME         "globalbnd"
#define EVENTHDLR_DESC         "event handler for globalbnd event"


/*
 * Data structures
 */



/** event handler data */
struct SCIP_EventhdlrData
{
   int                    filterpos;
   SCIP_Bool              storebounds;
   SCIP_BOUNDSTORE*       boundstore;
};

/*
 * Local methods
 */

/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeGlobalbnd)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitGlobalbnd)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->filterpos < 0 && SCIPgetSubscipDepth(scip) == 0 && SCIPsyncstoreIsInitialized(SCIPgetSyncstore(scip)) )
   {
      int        i;
      int        nvars;
      SCIP_VAR** vars;
      SCIPdebugMsg(scip, "catching events in " EVENTHDLR_NAME " eventhdlr\n");
      /* notify SCIP that this event handler wants to react on global bound change events */
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
      eventhdlrdata->storebounds = TRUE;
      SCIP_CALL( SCIPboundstoreCreate(scip, &eventhdlrdata->boundstore, SCIPgetNOrigVars(scip)) );

      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, &eventhdlrdata->filterpos) );
      for( i = 0; i < nvars ; ++i )
      {
         SCIP_CALL( SCIPcatchVarEvent(scip, vars[i], SCIP_EVENTTYPE_GBDCHANGED, eventhdlr, NULL, NULL) );
      }
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitGlobalbnd)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* notify SCIP that your event handler wants to drop the event type var added */
   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, eventhdlrdata->filterpos) );
      eventhdlrdata->filterpos = -1;
      SCIPboundstoreFree(scip, &eventhdlrdata->boundstore);
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecGlobalbnd)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_VAR*           var;
   SCIP_Real           newbound;
   SCIP_BOUNDTYPE      boundtype;
   SCIP_Real           constant;
   SCIP_Real           scalar;
   SCIPdebugMsg(scip, "exec method of eventhdlr " EVENTHDLR_NAME "\n");
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   var = SCIPeventGetVar(event);
   switch( SCIPeventGetType(event) )
   {
      case SCIP_EVENTTYPE_VARADDED:
         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GBDCHANGED, eventhdlr, NULL, NULL) );
         return SCIP_OKAY;
      case SCIP_EVENTTYPE_GLBCHANGED:
         boundtype = SCIP_BOUNDTYPE_LOWER;
         break;
      case SCIP_EVENTTYPE_GUBCHANGED:
         boundtype = SCIP_BOUNDTYPE_UPPER;
         break;
      default:
         SCIPABORT();
         return SCIP_ERROR; /*lint !e527*/
   }

   if( !eventhdlrdata->storebounds )
      return SCIP_OKAY;

   newbound = SCIPeventGetNewbound(event);
   constant = 0.0;
   scalar = 1.0;
   SCIP_CALL( SCIPvarGetOrigvarSum(&var, &scalar, &constant) );
   if( var != NULL )
   {
      int varidx;

      varidx = SCIPgetConcurrentVaridx(scip, var);

      boundtype = scalar < 0.0 ? SCIPboundtypeOpposite(boundtype) : boundtype;
      newbound = (newbound - constant) / scalar;

      SCIP_CALL( SCIPboundstoreAdd(scip, eventhdlrdata->boundstore, varidx, newbound, boundtype) );
   }
   return SCIP_OKAY;
}

/** creates event handler for globalbnd event */
SCIP_RETCODE SCIPincludeEventHdlrGlobalbnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create globalbnd event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   eventhdlrdata->filterpos = -1;
   eventhdlr = NULL;

   /* include event handler into SCIP */

   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecGlobalbnd, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeGlobalbnd) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitGlobalbnd) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitGlobalbnd) );

   return SCIP_OKAY;
}


/** gets the global bound changes stored in the eventhandler */
SCIP_BOUNDSTORE* SCIPeventGlobalbndGetBoundChanges(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(eventhdlr != NULL);
   assert(strcmp(EVENTHDLR_NAME, SCIPeventhdlrGetName(eventhdlr)) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return eventhdlrdata->boundstore;
}

/** enables storing of bound changes */
void SCIPeventGlobalbndEnableBoundStorage(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(eventhdlr != NULL);
   assert(strcmp(EVENTHDLR_NAME, SCIPeventhdlrGetName(eventhdlr)) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->storebounds = TRUE;
}

/** disables storing of bound changes */
void SCIPeventGlobalbndDisableBoundStorage(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(eventhdlr != NULL);
   assert(strcmp(EVENTHDLR_NAME, SCIPeventhdlrGetName(eventhdlr)) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->storebounds = FALSE;
}

/** clears all bound changes stored in the eventhandler */
void SCIPeventGlobalbndClearBoundChanges(
   SCIP_EVENTHDLR*       eventhdlr           /**< the globalbound eventhandler */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(eventhdlr != NULL);
   assert(strcmp(EVENTHDLR_NAME, SCIPeventhdlrGetName(eventhdlr)) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPboundstoreClear(eventhdlrdata->boundstore);
}
