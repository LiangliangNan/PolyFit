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

/**@file   objprop.h
 * @brief  C++ wrapper for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPROP_H__
#define __SCIP_OBJPROP_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for propagators
 *
 *  This class defines the interface for propagators implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_exec().
 *
 *  - \ref PROP "Instructions for implementing a propagator"
 *  - \ref PROPAGATORS "List of available propagators"
 *  - \ref type_prop.h "Corresponding C interface"
 */
class ObjProp : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the propagator */
   char* scip_name_;

   /** description of the propagator */
   char* scip_desc_;

   /** default priority of the propagator */
   const int scip_priority_;

   /** frequency for calling propagator */
   const int scip_freq_;

   /** should propagator be delayed, if other propagators found reductions? */
   const SCIP_Bool scip_delay_;

   /** positions in the node solving loop where propagator should be executed */
   const SCIP_PROPTIMING scip_timingmask_;

   /** default presolving priority of the propagator */
   const int scip_presol_priority_;

   /** frequency for calling propagator */
   const int scip_presol_maxrounds_;

   /**< timing mask of the propagator's presolving method */
   const SCIP_PRESOLTIMING scip_presol_timing_;


   /** default constructor */
   ObjProp(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of propagator */
      const char*        desc,               /**< description of propagator */
      int                priority,           /**< priority of the propagator */
      int                freq,               /**< frequency for calling propagator */
      SCIP_Bool          delay,              /**< should propagator be delayed, if other propagators found reductions? */
      SCIP_PROPTIMING    timingmask,         /**< positions in the node solving loop where propagator should be executed */
      int                presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
      int                presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
      SCIP_PRESOLTIMING  presoltiming        /**< timing mask of the propagator's presolving method */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_freq_(freq),
        scip_delay_(delay),
        scip_timingmask_(timingmask),
        scip_presol_priority_(presolpriority),
        scip_presol_maxrounds_(presolmaxrounds),
        scip_presol_timing_(presoltiming)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjProp(const ObjProp& o)
       : ObjProp(o.scip_, o.scip_name_, o.scip_desc_, o.scip_priority_, o.scip_freq_, o.scip_delay_, o.scip_timingmask_,
                 o.scip_presol_priority_, o.scip_presol_maxrounds_, o.scip_presol_timing_)
   {
   }

   /** move constructor */
   ObjProp(ObjProp&& o)
       : scip_(o.scip_),
         scip_name_(0),
         scip_desc_(0),
         scip_priority_(o.scip_priority_),
         scip_freq_(o.scip_freq_),
         scip_delay_(o.scip_delay_),
         scip_timingmask_(o.scip_timingmask_),
         scip_presol_priority_(o.scip_presol_priority_),
         scip_presol_maxrounds_(o.scip_presol_maxrounds_),
         scip_presol_timing_(o.scip_presol_timing_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjProp()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjProp& operator=(const ObjProp& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjProp& operator=(ObjProp&& o) = delete;

   /** destructor of propagator to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_PROPFREE(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of propagator (called after problem was transformed)
    *
    *  @see SCIP_DECL_PROPINIT(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of propagator (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_PROPEXIT(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving initialization method of propagator (called when presolving is about to begin)
    *
    *  @see SCIP_DECL_PROPINITPRE(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPINITPRE(scip_initpre)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving deinitialization method of propagator (called after presolving has been finished)
    *
    *  @see SCIP_DECL_PROPEXITPRE(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPEXITPRE(scip_exitpre)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of propagator (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_PROPINITSOL(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of propagator (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_PROPEXITSOL(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving method of propagator
    *
    *  @see SCIP_DECL_PROPPRESOL(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPPRESOL(scip_presol)
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** execution method of propagator
    *
    *  @see SCIP_DECL_PROPEXEC(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPEXEC(scip_exec) = 0;

   /** propagation conflict resolving method of propagator
    *
    *  @see SCIP_DECL_PROPRESPROP(x) in @ref type_prop.h
    */
   virtual SCIP_DECL_PROPRESPROP(scip_resprop)
   {  /*lint --e{715}*/

      /* set result pointer to indicate the propagation was not resolved */
      assert(result != NULL);
      (*result) = SCIP_DIDNOTFIND;

      return SCIP_OKAY;
   }
};

} /* namespace scip */



/** creates the propagator for the given propagator object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyProp* myprop = new MyProp(...);
 *       SCIP_CALL( SCIPincludeObjProp(scip, &myprop, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myprop;    // delete prop AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjProp(scip, new MyProp(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyProp is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjProp*        objprop,            /**< propagator object */
   SCIP_Bool             deleteobject        /**< should the propagator object be deleted when propagator is freed? */
   );

/** returns the prop object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjProp* SCIPfindObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of propagator */
   );

/** returns the prop object for the given propagator */
SCIP_EXPORT
scip::ObjProp* SCIPgetObjProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop                /**< propagator */
   );

#endif
