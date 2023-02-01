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

/**@file   objdialog.cpp
 * @brief  C++ wrapper for dialogs
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objdialog.h"




/*
 * Data structures
 */

/** dialog data */
struct SCIP_DialogData
{
   scip::ObjDialog*      objdialog;          /**< dialog object */
   SCIP_Bool             deleteobject;       /**< should the dialog object be deleted when dialog is freed? */
};




/*
 * Callback methods of dialog
 */

extern "C"
{

/** copy method for dialog plugins (called when SCIP copies plugins) */
static
SCIP_DECL_DIALOGCOPY(dialogCopyObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;
   
   assert(scip != 0);
   
   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != 0);
   assert(dialogdata->objdialog != 0);
   assert(dialogdata->objdialog->scip_ != scip);

   if( dialogdata->objdialog->iscloneable() )
   {
      scip::ObjDialog*  newobjdialog;
      newobjdialog = dynamic_cast<scip::ObjDialog*> (dialogdata->objdialog->clone(scip));

      /* call include method of dialog object */
      SCIP_CALL( SCIPincludeObjDialog(scip, newobjdialog, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of dialog to free user data (called when SCIP is exiting) */
static
SCIP_DECL_DIALOGFREE(dialogFreeObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != 0);
   assert(dialogdata->objdialog != 0);
   assert(dialogdata->objdialog->scip_ == scip);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_free(scip, dialog) );

   /* free dialog object */
   if( dialogdata->deleteobject )
      delete dialogdata->objdialog;

   /* free dialog data */
   delete dialogdata;
   SCIPdialogSetData(dialog, 0); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** description output method of dialog */
static
SCIP_DECL_DIALOGDESC(dialogDescObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != 0);
   assert(dialogdata->objdialog != 0);
   assert(dialogdata->objdialog->scip_ == scip);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_desc(scip, dialog) );

   return SCIP_OKAY;
}

/** execution method of dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecObj)
{  /*lint --e{715}*/
   SCIP_DIALOGDATA* dialogdata;

   dialogdata = SCIPdialogGetData(dialog);
   assert(dialogdata != 0);
   assert(dialogdata->objdialog != 0);

   /* call virtual method of dialog object */
   SCIP_CALL( dialogdata->objdialog->scip_exec(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}
}



/*
 * dialog specific interface methods
 */

/** creates the dialog for the given dialog object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjDialog*      objdialog,          /**< dialog object */
   SCIP_Bool             deleteobject        /**< should the dialog object be deleted when dialog is freed? */
   )
{/*lint --e{429} */
   SCIP_DIALOG* parentdialog;

   assert(scip != 0);
   assert(objdialog != 0);

   /* get parent dialog */
   parentdialog = SCIPgetRootDialog(scip);
   assert(parentdialog != 0);
   /* TODO: (optional) change parent dialog from root dialog to another existing dialog (needs to be a menu) */

   /* create, include, and release dialog */
   if( !SCIPdialogHasEntry(parentdialog, objdialog->scip_name_) )
   {
      SCIP_DIALOGDATA* dialogdata; /*lint !e593*/
      SCIP_DIALOG* dialog;
      SCIP_RETCODE retcode;

      dialog = 0;

      /* create dialog data */
      dialogdata = new SCIP_DIALOGDATA;
      dialogdata->objdialog = objdialog;
      dialogdata->deleteobject = deleteobject;

      retcode = SCIPincludeDialog(scip, &dialog, dialogCopyObj, dialogExecObj, dialogDescObj, dialogFreeObj,
         objdialog->scip_name_, objdialog->scip_desc_, objdialog->scip_issubmenu_, dialogdata);
      if( retcode != SCIP_OKAY )
      {
          delete dialogdata;
          SCIP_CALL( retcode );
      }
      SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );  /*lint !e593*/
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );                /*lint !e593*/
   }  /*lint !e593*/

   return SCIP_OKAY;  /*lint !e593*/
}
