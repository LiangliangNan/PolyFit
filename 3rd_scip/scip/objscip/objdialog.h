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

/**@file   objdialog.h
 * @brief  C++ wrapper for dialogs
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJDIALOG_H__
#define __SCIP_OBJDIALOG_H__

#include <cstring>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/**
 *  @brief C++ wrapper for dialogs
 *
 *  This class defines the interface for dialogs implemented in C++. Note that there is a pure virtual function (this
 *  function has to be implemented). This function is: scip_exec().
 *
 * - \ref DIALOG "Instructions for implementing a dialog"
 * - \ref DIALOGS "List of available dialogs"
 *  - \ref type_dialog.h "Corresponding C interface"
 */
class ObjDialog : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the dialog */
   char* scip_name_;

   /** description of the dialog */
   char* scip_desc_;

   /** default for whether the dialog is a menu */
   const SCIP_Bool scip_issubmenu_;

   /** default constructor */
   ObjDialog(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of the dialog */
      const char*        desc,               /**< description of the dialog */
      SCIP_Bool          issubmenu           /**< default for whether the dialog is a menu */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_issubmenu_(issubmenu)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjDialog()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of dialog to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_DIALOGFREE(x) in @ref type_dialog.h
    */
   virtual SCIP_DECL_DIALOGFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** description output method of dialog
    *
    *  @see SCIP_DECL_DIALOGDESC(x) in @ref type_dialog.h
    */
   virtual SCIP_DECL_DIALOGDESC(scip_desc)
   {  /*lint --e{715}*/
      SCIPdialogMessage(scip, 0, "%s", scip_desc_);
      return SCIP_OKAY;
   }

   /** execution method of dialog
    *
    *  @see SCIP_DECL_DIALOGEXEC(x) in @ref type_dialog.h
    */
   virtual SCIP_DECL_DIALOGEXEC(scip_exec) = 0;
};

} /* namespace scip */



/** creates the dialog for the given dialog object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyDialog* mydialog = new MyDialog(...);
 *       SCIP_CALL( SCIPincludeObjDialog(scip, &mydialog, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mydialog;    // delete dialog AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjDialog(scip, new MyDialog(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyDialog is called here
 */
EXTERN
SCIP_RETCODE SCIPincludeObjDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjDialog*      objdialog,          /**< dialog object */
   SCIP_Bool             deleteobject        /**< should the dialog object be deleted when dialog is freed? */
   );

#endif
