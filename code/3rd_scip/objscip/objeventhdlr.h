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

/**@file   objeventhdlr.h
 * @brief  C++ wrapper for event handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJEVENTHDLR_H__
#define __SCIP_OBJEVENTHDLR_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for event handlers
 *
 *  This class defines the interface for eventt handlers implemented in C++. Note that there is a pure virtual function
 *  (this function has to be implemented). This function is: scip_exec().
 *
 *  - \ref EVENT "Instructions for implementing an event handler"
 *  - \ref type_event.h "Corresponding C interface"
 */
class ObjEventhdlr : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the event handler */
   char* scip_name_;

   /** description of the event handler */
   char* scip_desc_;

   /** default constructor */
   ObjEventhdlr(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of event handler */
      const char*        desc                /**< description of event handler */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjEventhdlr(const ObjEventhdlr& o) : ObjEventhdlr(o.scip_, o.scip_name_, o.scip_desc_) {}

   /** move constructor */
   ObjEventhdlr(ObjEventhdlr&& o) : scip_(o.scip_), scip_name_(0), scip_desc_(0)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjEventhdlr()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjEventhdlr& operator=(const ObjEventhdlr& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjEventhdlr& operator=(ObjEventhdlr&& o) = delete;

   /** destructor of event handler to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_EVENTFREE(x) in @ref type_event.h
    */
   virtual SCIP_DECL_EVENTFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of event handler (called after problem was transformed)
    *
    *  @see SCIP_DECL_EVENTINIT(x) in @ref type_event.h
    */
   virtual SCIP_DECL_EVENTINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of event handler (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_EVENTEXIT(x) in @ref type_event.h
    */
   virtual SCIP_DECL_EVENTEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of event handler (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_EVENTINITSOL(x) in @ref type_event.h
    */
   virtual SCIP_DECL_EVENTINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of event handler (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_EVENTEXITSOL(x) in @ref type_event.h
    */
   virtual SCIP_DECL_EVENTEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** frees specific constraint data
    *
    *  @see SCIP_DECL_EVENTDELETE(x) in @ref type_event.h
    */
   virtual SCIP_DECL_EVENTDELETE(scip_delete)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** execution method of event handler
    *
    *  @see SCIP_DECL_EVENTEXEC(x) in @ref type_event.h
    */
   virtual SCIP_DECL_EVENTEXEC(scip_exec) = 0;
};

} /* namespace scip */



/** creates the event handler for the given event handler object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyEventhdlr* myeventhdlr = new MyEventhdlr(...);
 *       SCIP_CALL( SCIPincludeObjEventhdlr(scip, &myeventhdlr, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myeventhdlr;    // delete eventhdlr AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjEventhdlr(scip, new MyEventhdlr(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyEventhdlr is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjEventhdlr*   objeventhdlr,       /**< event handler object */
   SCIP_Bool             deleteobject        /**< should the event handler object be deleted when eventhdlristic is freed? */
   );

/** returns the eventhdlr object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjEventhdlr* SCIPfindObjEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of event handler */
   );

/** returns the eventhdlr object for the given event handler */
SCIP_EXPORT
scip::ObjEventhdlr* SCIPgetObjEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   );

#endif
