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

/**@file   objmessagehdlr.h
 * @brief  C++ wrapper for message handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJMESSAGEHDLR_H__
#define __SCIP_OBJMESSAGEHDLR_H__

#include "scip/scip.h"

namespace scip
{

/**
 *  @brief C++ wrapper for message handlers
 *
 *  This class defines the interface for message handlers implemented in C++. Note that all functions are pure virtual
 *  (these functions have to be implemented).
 *
 *  - \ref type_message.h "Corresponding C interface"
 */
class ObjMessagehdlr
{
public:
   /** should the output be buffered up to the next newline? */
   const SCIP_Bool scip_bufferedoutput_;

   /** default constructor */
   explicit ObjMessagehdlr(
      SCIP_Bool          bufferedoutput      /**< should the output be buffered up to the next newline? */
      )
      : scip_bufferedoutput_(bufferedoutput)
   {
   }

   /** destructor */
   virtual ~ObjMessagehdlr()
   {
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjMessagehdlr& operator=(const ObjMessagehdlr& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjMessagehdlr& operator=(ObjMessagehdlr&& o) = delete;

   /** error message print method of message handler
    *
    *  @note This function can be activated by calling SCIPsetStaticErrorPrintingMessagehdlr().
    *
    *  @see SCIP_DECL_ERRORPRINTING(x) in @ref type_message.h
    */ /*lint -e715*/
   virtual void scip_error(
      SCIP_MESSAGEHDLR*  messagehdlr,        /**< the message handler itself */
      FILE*              file,               /**< file stream to print into (NULL for stderr) */
      const char*        msg                 /**< string to output into the file (or NULL to flush) */
      )
   { /*lint --e{715}*/

   }

   /** warning message print method of message handler
    *
    *  @see SCIP_DECL_MESSAGEWARNING(x) in @ref type_message.h
    */ /*lint -e715*/
   virtual SCIP_DECL_MESSAGEWARNING(scip_warning)
   { /*lint --e{715}*/

   }

   /** dialog message print method of message handler
    *
    *  @see SCIP_DECL_MESSAGEDIALOG(x) in @ref type_message.h
    */ /*lint -e715*/
   virtual SCIP_DECL_MESSAGEDIALOG(scip_dialog)
   { /*lint --e{715}*/

   }

   /** info message print method of message handler
    *
    *  @see SCIP_DECL_MESSAGEINFO(x) in @ref type_message.h
    */ /*lint -e715*/
   virtual SCIP_DECL_MESSAGEINFO(scip_info)
   { /*lint --e{715}*/

   }

   /** destructor of message handler to free message handler data
    *
    *  @see SCIP_DECL_MESSAGEHDLRFREE(x) in @ref type_message.h
    */ /*lint -e715*/
   virtual SCIP_DECL_MESSAGEHDLRFREE(scip_free)
   { /*lint --e{715}*/
      return SCIP_OKAY;
   }
};

} /* namespace scip */



/** creates the message handler for the given message handler object */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateObjMessagehdlr(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   scip::ObjMessagehdlr* objmessagehdlr,     /**< message handler object */
   SCIP_Bool             deleteobject        /**< should the message handler object be deleted when message handler is freed? */
   );

/** returns the message handler object for the given message handler */
SCIP_EXPORT
scip::ObjMessagehdlr* SCIPgetObjMessagehdlr(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** set static error output function to the corresponding function of message handler */
SCIP_EXPORT
void SCIPsetStaticErrorPrintingMessagehdlr(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

#endif
