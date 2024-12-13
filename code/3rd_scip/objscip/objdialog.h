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

/**@file   objdialog.h
 * @brief  C++ wrapper for dialogs
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJDIALOG_H__
#define __SCIP_OBJDIALOG_H__

#include <cstring>
#include <utility>

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

   /** copy constructor */
   ObjDialog(const ObjDialog& o) : ObjDialog(o.scip_, o.scip_name_, o.scip_desc_, o.scip_issubmenu_) {}

   /** move constructor */
   ObjDialog(ObjDialog&& o) : scip_(o.scip_), scip_name_(0), scip_desc_(0), scip_issubmenu_(o.scip_issubmenu_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjDialog()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjDialog& operator=(const ObjDialog& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjDialog& operator=(ObjDialog&& o) = delete;

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
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjDialog*      objdialog,          /**< dialog object */
   SCIP_Bool             deleteobject        /**< should the dialog object be deleted when dialog is freed? */
   );

#endif
