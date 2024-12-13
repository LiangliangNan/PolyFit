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

/**@file   objpricer.h
 * @brief  C++ wrapper for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJPRICER_H__
#define __SCIP_OBJPRICER_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objprobcloneable.h"

namespace scip
{

/** @brief C++ wrapper for variable pricer
 *
 *  This class defines the interface for variable pricer implemented in C++. Note that there is a pure virtual
 *  function (this function has to be implemented). This function is: scip_redcost().
 *
 *  - \ref PRICER "Instructions for implementing a variable pricer"
 *  - \ref type_pricer.h "Corresponding C interface"
 */
class ObjPricer : public ObjProbCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the variable pricer */
   char* scip_name_;

   /** description of the variable pricer */
   char* scip_desc_;

   /** default priority of the variable pricer */
   const int scip_priority_;

   /** should the pricer be delayed until no other pricers or already existing problem variables with negative reduced
    *  costs are found?
    */
   const SCIP_Bool scip_delay_;

   /** default constructor */
   ObjPricer(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of variable pricer */
      const char*        desc,               /**< description of variable pricer */
      int                priority,           /**< priority of the variable pricer */
      SCIP_Bool          delay               /**< should the pricer be delayed until no other pricers or already existing
                                              *   problem variables with negative reduced costs are found?
                                              *   if this is set to FALSE it may happen that the pricer produces columns
                                              *   that already exist in the problem (which are also priced in by the
                                              *   default problem variable pricing in the same round)
                                              */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_priority_(priority),
        scip_delay_(delay)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjPricer(const ObjPricer& o) : ObjPricer(o.scip_, o.scip_name_, o.scip_desc_, o.scip_priority_, o.scip_delay_) {}

   /** move constructor */
   ObjPricer(ObjPricer&& o)
       : scip_(o.scip_), scip_name_(0), scip_desc_(0), scip_priority_(o.scip_priority_), scip_delay_(o.scip_delay_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjPricer()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjPricer& operator=(const ObjPricer& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjPricer& operator=(ObjPricer&& o) = delete;

   /** destructor of variable pricer to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_PRICERFREE(x) in @ref type_pricer.h
    */
   virtual SCIP_DECL_PRICERFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of variable pricer (called after problem was transformed)
    *
    *  @see SCIP_DECL_PRICERINIT(x) in @ref type_pricer.h
    */
   virtual SCIP_DECL_PRICERINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of variable pricer (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_PRICEREXIT(x) in @ref type_pricer.h
    */
   virtual SCIP_DECL_PRICEREXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of variable pricer (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_PRICERINITSOL(x) in @ref type_pricer.h
    */
   virtual SCIP_DECL_PRICERINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of variable pricer (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_PRICEREXITSOL(x) in @ref type_pricer.h
    */
   virtual SCIP_DECL_PRICEREXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** reduced cost pricing method of variable pricer for feasible LPs
    *
    *  @see SCIP_DECL_PRICERREDCOST(x) in @ref type_pricer.h
    */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost) = 0;

   /** farkas pricing method of variable pricer for infeasible LPs
    *
    *  @see SCIP_DECL_PRICERFARKAS(x) in @ref type_pricer.h
    */
   virtual SCIP_DECL_PRICERFARKAS(scip_farkas)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
};

} /* namespace scip */



/** creates the variable pricer for the given variable pricer object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyPricer* mypricer = new MyPricer(...);
 *       SCIP_CALL( SCIPincludeObjPricer(scip, &mypricer, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mypricer;    // delete pricer AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjPricer(scip, new MyPricer(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyPricer is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjPricer*      objpricer,          /**< variable pricer object */
   SCIP_Bool             deleteobject        /**< should the pricer object be deleted when pricer is freed? */
   );

/** returns the variable pricer object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjPricer* SCIPfindObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable pricer */
   );

/** returns the variable pricer object for the given pricer */
SCIP_EXPORT
scip::ObjPricer* SCIPgetObjPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer              /**< pricer */
   );

#endif
