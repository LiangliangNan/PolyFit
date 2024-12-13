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

/**@file   event_xyz.c
 * @ingroup DEFPLUGINS_EVENT
 * @brief  eventhdlr for xyz event
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_xyz.h"

#define EVENTHDLR_NAME         "xyz"
#define EVENTHDLR_DESC         "event handler for xyz event"


/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */

/** event handler data */
struct SCIP_EventhdlrData
{
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */



/*
 * Callback methods of event handler
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_EVENTCOPY(eventCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventCopyXyz NULL
#endif

/** destructor of event handler to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_EVENTFREE(eventFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventFreeXyz NULL
#endif

/** initialization method of event handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_EVENTINIT(eventInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventInitXyz NULL
#endif

/** deinitialization method of event handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_EVENTEXIT(eventExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitXyz NULL
#endif

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_EVENTINITSOL(eventInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventInitsolXyz NULL
#endif

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_EVENTEXITSOL(eventExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitsolXyz NULL
#endif

/** frees specific event data */
#if 0
static
SCIP_DECL_EVENTDELETE(eventDeleteXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventDeleteXyz NULL
#endif

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** creates event handler for xyz event */
SCIP_RETCODE SCIPincludeEventHdlrXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create xyz event handler data */
   eventhdlrdata = NULL;
   /* TODO: (optional) create event handler specific data here */

   eventhdlr = NULL;

   /* include event handler into SCIP */
#if 0
   /* use SCIPincludeEventhdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventCopyXyz,
         eventFreeXyz, eventInitXyz, eventExitXyz, 
         eventInitsolXyz, eventExitsolXyz, eventDeleteXyz, eventExecXyz,
         eventhdlrdata) );
#else
   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecXyz, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyXyz) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeXyz) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitXyz) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitXyz) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolXyz) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolXyz) );
   SCIP_CALL( SCIPsetEventhdlrDelete(scip, eventhdlr, eventDeleteXyz) );
#endif

   /* add xyz event handler parameters */
   /* TODO: (optional) add event handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
