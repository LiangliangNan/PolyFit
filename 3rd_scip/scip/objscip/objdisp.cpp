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

/**@file   objdisp.cpp
 * @brief  C++ wrapper for display column
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objdisp.h"




/*
 * Data structures
 */

/** display column data */
struct SCIP_DispData
{
   scip::ObjDisp*        objdisp;            /**< display column object */
   SCIP_Bool             deleteobject;       /**< should the display column object be deleted when display column is freed? */
};




/*
 * Callback methods of display column
 */

extern "C"
{

/** copy method for display column plugins (called when SCIP copies plugins) */
static
SCIP_DECL_DISPCOPY(dispCopyObj)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;
   
   assert(scip != NULL);
   
   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);
   assert(dispdata->objdisp != NULL);
   assert(dispdata->objdisp->scip_ != scip);

   if( dispdata->objdisp->iscloneable() )
   {
      scip::ObjDisp*  newobjdisp;
      newobjdisp = dynamic_cast<scip::ObjDisp*> (dispdata->objdisp->clone(scip));

      /* call include method of display column object */
      SCIP_CALL( SCIPincludeObjDisp(scip, newobjdisp, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of display column to free user data (called when SCIP is exiting) */
static
SCIP_DECL_DISPFREE(dispFreeObj)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);
   assert(dispdata->objdisp != NULL);
   assert(dispdata->objdisp->scip_ == scip);

   /* call virtual method of display column object */
   SCIP_CALL( dispdata->objdisp->scip_free(scip, disp) );

   /* free display column object */
   if( dispdata->deleteobject )
      delete dispdata->objdisp;

   /* free display column data */
   delete dispdata;
   SCIPdispSetData(disp, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** initialization method of display column (called after problem was transformed) */
static
SCIP_DECL_DISPINIT(dispInitObj)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);
   assert(dispdata->objdisp != NULL);
   assert(dispdata->objdisp->scip_ == scip);

   /* call virtual method of display column object */
   SCIP_CALL( dispdata->objdisp->scip_init(scip, disp) );

   return SCIP_OKAY;
}


/** deinitialization method of display column (called before transformed problem is freed) */
static
SCIP_DECL_DISPEXIT(dispExitObj)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);
   assert(dispdata->objdisp != NULL);

   /* call virtual method of display column object */
   SCIP_CALL( dispdata->objdisp->scip_exit(scip, disp) );

   return SCIP_OKAY;
}


/** solving process initialization method of display column (called when branch and bound process is about to begin) */
static
SCIP_DECL_DISPINITSOL(dispInitsolObj)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);
   assert(dispdata->objdisp != NULL);

   /* call virtual method of display column object */
   SCIP_CALL( dispdata->objdisp->scip_initsol(scip, disp) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of display column (called before branch and bound process data is freed) */
static
SCIP_DECL_DISPEXITSOL(dispExitsolObj)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);
   assert(dispdata->objdisp != NULL);

   /* call virtual method of display column object */
   SCIP_CALL( dispdata->objdisp->scip_exitsol(scip, disp) );

   return SCIP_OKAY;
}


/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputObj)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);
   assert(dispdata->objdisp != NULL);

   /* call virtual method of display column object */
   SCIP_CALL( dispdata->objdisp->scip_output(scip, disp, file) );

   return SCIP_OKAY;
}
}



/*
 * display column specific interface methods
 */

/** creates the display column for the given display column object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjDisp*        objdisp,            /**< display column object */
   SCIP_Bool             deleteobject        /**< should the display column object be deleted when display column is freed? */
   )
{
   SCIP_DISPDATA* dispdata;

   assert(scip != NULL);
   assert(objdisp != NULL);
   
   /* create display column data */
   dispdata = new SCIP_DISPDATA;
   dispdata->objdisp = objdisp;
   dispdata->deleteobject = deleteobject;

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, objdisp->scip_name_, objdisp->scip_desc_, 
         objdisp->scip_header_, SCIP_DISPSTATUS_AUTO,
         dispCopyObj,
         dispFreeObj, dispInitObj, dispExitObj, dispInitsolObj, 
         dispExitsolObj, dispOutputObj, dispdata, objdisp->scip_width_, objdisp->scip_priority_, objdisp->scip_position_, 
         objdisp->scip_stripline_) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the display column object of the given name, or 0 if not existing */
scip::ObjDisp* SCIPfindObjDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of display column */
   )
{
   SCIP_DISP* disp;
   SCIP_DISPDATA* dispdata;

   disp = SCIPfindDisp(scip, name);
   if( disp == NULL )
      return 0;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);

   return dispdata->objdisp;
}
   
/** returns the display column object for the given display column */
scip::ObjDisp* SCIPgetObjDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DISP*            disp                /**< display column */
   )
{
   SCIP_DISPDATA* dispdata;

   dispdata = SCIPdispGetData(disp);
   assert(dispdata != NULL);

   return dispdata->objdisp;
}
