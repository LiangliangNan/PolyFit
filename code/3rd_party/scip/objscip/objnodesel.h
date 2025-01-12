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

/**@file   objnodesel.h
 * @brief  C++ wrapper for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJNODESEL_H__
#define __SCIP_OBJNODESEL_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for primal heuristics
 *
 *  This class defines the interface for node selectors implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_comp().
 *
 *  - \ref NODESEL "Instructions for implementing a  node selector"
 *  - \ref NODESELECTORS "List of available node selectors"
 *  - \ref type_nodesel.h "Corresponding C interface"
 */
class ObjNodesel : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the node selector */
   char* scip_name_;

   /** description of the node selector */
   char* scip_desc_;

   /** priority of the node selector in standard mode */
   const int scip_stdpriority_;

   /** priority of the node selector in memory saving mode */
   const int scip_memsavepriority_;

   /** default constructor */
   ObjNodesel(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of node selector */
      const char*        desc,               /**< description of node selector */
      int                stdpriority,        /**< priority of the node selector in standard mode */
      int                memsavepriority     /**< priority of the node selector in memory saving mode */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_stdpriority_(stdpriority),
        scip_memsavepriority_(memsavepriority)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjNodesel(const ObjNodesel& o)
       : ObjNodesel(o.scip_, o.scip_name_, o.scip_desc_, o.scip_stdpriority_, o.scip_memsavepriority_)
   {
   }

   /** move constructor */
   ObjNodesel(ObjNodesel&& o)
       : scip_(o.scip_),
         scip_name_(0),
         scip_desc_(0),
         scip_stdpriority_(o.scip_stdpriority_),
         scip_memsavepriority_(o.scip_memsavepriority_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjNodesel()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjNodesel& operator=(const ObjNodesel& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjNodesel& operator=(ObjNodesel&& o) = delete;

   /** destructor of node selector to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_NODESELFREE(x) in @ref type_nodesel.h
    */
   virtual SCIP_DECL_NODESELFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of node selector (called after problem was transformed)
    *
    *  @see SCIP_DECL_NODESELINIT(x) in @ref type_nodesel.h
    */
   virtual SCIP_DECL_NODESELINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of node selector (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_NODESELEXIT(x) in @ref type_nodesel.h
    */
   virtual SCIP_DECL_NODESELEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of node selector (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_NODESELINITSOL(x) in @ref type_nodesel.h
    */
   virtual SCIP_DECL_NODESELINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of node selector (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_NODESELEXITSOL(x) in @ref type_nodesel.h
    */
   virtual SCIP_DECL_NODESELEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** node selection method of node selector
    *
    *  @see SCIP_DECL_NODESELSELECT(x) in @ref type_nodesel.h
    */
   virtual SCIP_DECL_NODESELSELECT(scip_select) = 0;

   /** node comparison method of node selector
    *
    *  @see SCIP_DECL_NODESELCOMP(x) in @ref type_nodesel.h
    */
   virtual SCIP_DECL_NODESELCOMP(scip_comp) = 0;
};

} /* namespace scip */



/** creates the node selector for the given node selector object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyNodesel* mynodesel = new MyNodesel(...);
 *       SCIP_CALL( SCIPincludeObjNodesel(scip, &mynodesel, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mynodesel;    // delete nodesel AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjNodesel(scip, new MyNodesel(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyNodesel is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjNodesel*     objnodesel,         /**< node selector object */
   SCIP_Bool             deleteobject        /**< should the node selector object be deleted when node selector is freed? */
   );

/** returns the nodesel object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjNodesel* SCIPfindObjNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of node selector */
   );

/** returns the nodesel object for the given node selector */
SCIP_EXPORT
scip::ObjNodesel* SCIPgetObjNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel             /**< node selector */
   );

#endif
