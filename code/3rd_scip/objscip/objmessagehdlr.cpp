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

/**@file   objmessagehdlr.cpp
 * @brief  C++ wrapper for message handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objmessagehdlr.h"




/*
 * Data structures
 */

/** message handler data */
struct SCIP_MessagehdlrData
{
   scip::ObjMessagehdlr* objmessagehdlr;     /**< message handler object */
   SCIP_Bool             deleteobject;       /**< should the message handler object be deleted when message handler is freed? */
};




/*
 * Callback methods of file reader
 */

extern "C"
{
/** error message print method of message handler */
static
SCIP_DECL_ERRORPRINTING(messagehdlrErrorObj)
{
   assert( data != 0 );

   scip::ObjMessagehdlr* objmessagehdlr = static_cast<scip::ObjMessagehdlr*>(data);

   objmessagehdlr->scip_error(0, file, msg);
}

/** warning message print method of message handler */
static
SCIP_DECL_MESSAGEWARNING(messagehdlrWarningObj)
{  /*lint --e{715}*/
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL && messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   messagehdlrdata->objmessagehdlr->scip_warning(messagehdlr, file, msg);
}


/** dialog message print method of message handler */
static
SCIP_DECL_MESSAGEDIALOG(messagehdlrDialogObj)
{  /*lint --e{715}*/
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL && messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   messagehdlrdata->objmessagehdlr->scip_dialog(messagehdlr, file, msg);
}


/** info message print method of message handler */
static
SCIP_DECL_MESSAGEINFO(messagehdlrInfoObj)
{  /*lint --e{715}*/
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL && messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   messagehdlrdata->objmessagehdlr->scip_info(messagehdlr, file, msg);
}

/** destructor of message handler to free message handler data */
static
SCIP_DECL_MESSAGEHDLRFREE(messagehdlrFree)
{  /*lint --e{715}*/
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL && messagehdlrdata->objmessagehdlr != NULL);

   /* call virtual method of messagehdlr object */
   SCIP_CALL( messagehdlrdata->objmessagehdlr->scip_free(messagehdlr) );

   /* free message handler object */
   if( messagehdlrdata->deleteobject )
      delete messagehdlrdata->objmessagehdlr;

   /* free message handler data */
   delete messagehdlrdata;
   SCIP_CALL( SCIPmessagehdlrSetData(messagehdlr, NULL) ); /*lint !e64*/

   return SCIP_OKAY;
}
}



/*
 * message handler specific interface methods
 */

/** creates the message handler for the given message handler object */
SCIP_RETCODE SCIPcreateObjMessagehdlr(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   scip::ObjMessagehdlr* objmessagehdlr,     /**< message handler object */
   SCIP_Bool             deleteobject        /**< should the message handler object be deleted when message handler is freed? */
   )
{
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;
   SCIP_RETCODE retcode;

   /* create file messagehdlr data */
   messagehdlrdata = new SCIP_MESSAGEHDLRDATA;
   messagehdlrdata->objmessagehdlr = objmessagehdlr;
   messagehdlrdata->deleteobject = deleteobject;

   /* create message handler */
   retcode = SCIPmessagehdlrCreate(messagehdlr, objmessagehdlr->scip_bufferedoutput_, (const char*)NULL, FALSE,
      messagehdlrWarningObj, messagehdlrDialogObj, messagehdlrInfoObj,
      messagehdlrFree, messagehdlrdata); /*lint !e429*/

   if( retcode != SCIP_OKAY )
   {
      /* free message handler object */
      if( messagehdlrdata->deleteobject )
         delete messagehdlrdata->objmessagehdlr;

      delete messagehdlrdata;
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY; /*lint !e429 !e593*/
}

/** returns the message handler object for the given message handler */
scip::ObjMessagehdlr* SCIPgetObjMessagehdlr(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;

   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL && messagehdlrdata->objmessagehdlr != NULL);

   return messagehdlrdata->objmessagehdlr;
}

/** set static error output function to the corresponding function of message handler */
void SCIPsetStaticErrorPrintingMessagehdlr(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert( messagehdlr != NULL );

   SCIPmessageSetErrorPrinting(messagehdlrErrorObj, (void*) SCIPgetObjMessagehdlr(messagehdlr));
}
