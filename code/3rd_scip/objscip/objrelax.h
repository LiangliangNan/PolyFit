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

/**@file   objrelax.h
 * @brief  C++ wrapper for relaxation handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJRELAX_H__
#define __SCIP_OBJRELAX_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for relaxation handlers
 *
 *  This class defines the interface for relaxation handlers implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_exec().
 *
 *  - \ref RELAX "Instructions for implementing a relaxation handler"
 *  - \ref type_relax.h "Corresponding C interface"
 */
class ObjRelax : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the relaxator */
   char* scip_name_;

   /** description of the relaxator */
   char* scip_desc_;

   /** default priority of the relaxator (negative: call after LP, non-negative: call before LP) */
   const int scip_priority_;

   /** frequency for calling relaxator */
   const int scip_freq_;

   /** does the relaxator contain all cuts in the LP? */
   const SCIP_Bool scip_includeslp_;

   /** default constructor */
   ObjRelax(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of relaxator */
      const char*        desc,               /**< description of relaxator */
      int                priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
      int                freq,               /**< frequency for calling relaxator */
      SCIP_Bool          includeslp          /**< Does the relaxator contain all cuts in the LP? */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_freq_(freq),
        scip_includeslp_(includeslp)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjRelax(const ObjRelax& o)
       : ObjRelax(o.scip_, o.scip_name_, o.scip_desc_, o.scip_priority_, o.scip_priority_, o.scip_includeslp_)
   {
   }

   /** move constructor */
   ObjRelax(ObjRelax&& o)
       : scip_(o.scip_),
         scip_name_(0),
         scip_desc_(0),
         scip_priority_(o.scip_priority_),
         scip_freq_(o.scip_freq_),
         scip_includeslp_(o.scip_includeslp_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjRelax()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjRelax& operator=(const ObjRelax& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjRelax& operator=(ObjRelax&& o) = delete;

   /** destructor of relaxator to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_RELAXFREE(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of relaxator (called after problem was transformed)
    *
    *  @see SCIP_DECL_RELAXINIT(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of relaxator (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_RELAXEXIT(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of relaxator (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_RELAXINITSOL(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of relaxator (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_RELAXEXITSOL(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** execution method of relaxator
    *
    *  @see SCIP_DECL_RELAXEXEC(x) in @ref type_relax.h
    */
   virtual SCIP_DECL_RELAXEXEC(scip_exec) = 0;
};

} /* namespace scip */



/** creates the relaxator for the given relaxator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyRelax* myrelax = new MyRelax(...);
 *       SCIP_CALL( SCIPincludeObjRelax(scip, &myrelax, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myrelax;    // delete relax AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjRelax(scip, new MyRelax(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyRelax is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjRelax*  objrelax,           /**< relaxator object */
   SCIP_Bool             deleteobject        /**< should the relaxator object be deleted when relaxator is freed? */
   );

/** returns the relax object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjRelax* SCIPfindObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxator */
   );

/** returns the relax object for the given relaxator */
SCIP_EXPORT
scip::ObjRelax* SCIPgetObjRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< relaxator */
   );

#endif
