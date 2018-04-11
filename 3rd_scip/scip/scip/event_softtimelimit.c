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

/**@file   event_softtimelimit.c
 * @brief  eventhdlr for soft time limit
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "event_softtimelimit.h"

#include <string.h>

#define EVENTHDLR_NAME         "softtimelimit"
#define EVENTHDLR_DESC         "event handler for soft time limit"

/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Real softtimelimit;
   int filterpos;
};

/*
 * Callback methods of event handler
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
/**! [SnippetEventCopySofttimelimit] */
static
SCIP_DECL_EVENTCOPY(eventCopySofttimelimit)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* call inclusion method of event handler */
   SCIP_CALL( SCIPincludeEventHdlrSofttimelimit(scip) );

   return SCIP_OKAY;
}
/**! [SnippetEventCopySofttimelimit] */

/** destructor of event handler to free user data (called when SCIP is exiting) */
/**! [SnippetEventFreeSofttimelimit] */
static
SCIP_DECL_EVENTFREE(eventFreeSofttimelimit)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}
/**! [SnippetEventFreeSofttimelimit] */



/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitSofttimelimit)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->filterpos < 0 && !SCIPisNegative(scip, eventhdlrdata->softtimelimit) )
   {
      /* notify SCIP that your event handler wants to react on the event type best solution found */
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, &eventhdlrdata->filterpos) );
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitSofttimelimit)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* notify SCIP that your event handler wants to drop the event type best solution found */
   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, eventhdlrdata->filterpos) );
      eventhdlrdata->filterpos = -1;
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSofttimelimit)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_Real timelimit;

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPdebugMsg(scip, "exec method of event handler for soft time limit\n");

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   if( eventhdlrdata->softtimelimit < timelimit )
   {
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", eventhdlrdata->softtimelimit) );
   }

   /* notify SCIP that your event handler wants to drop the event type best solution found */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, eventhdlrdata->filterpos) );
   eventhdlrdata->filterpos = -1;

   /* print best solution value */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "changed time limit to %.1f after first solution was found\n",
      eventhdlrdata->softtimelimit);

   return SCIP_OKAY;
}

/** includes event handler for best solution found */
SCIP_RETCODE SCIPincludeEventHdlrSofttimelimit(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata) );
   eventhdlrdata->filterpos = -1;

   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSofttimelimit, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopySofttimelimit) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSofttimelimit) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSofttimelimit) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitSofttimelimit) );

   SCIP_CALL( SCIPaddRealParam(scip, "limits/softtime",
         "soft time limit which should be applied after first solution was found (-1.0: disabled)",
         &eventhdlrdata->softtimelimit, FALSE, -1.0, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
