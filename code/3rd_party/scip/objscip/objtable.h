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

/**@file   objtable.h
 * @brief  C++ wrapper for statistics tables
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJTABLE_H__
#define __SCIP_OBJTABLE_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"
namespace scip
{

/**
 *  @brief C++ wrapper for statistics tables
 *
 *  This class defines the interface for statistics tables implemented in C++. Note that there is a pure virtual function
 *  (this function has to be implemented). This function is: scip_output().
 *
 * - \ref TABLE "Instructions for implementing a statistics table"
 * - \ref TABLES "List of available statistics tables"
 * - \ref type_table.h "Corresponding C interface"
 */
class ObjTable : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the statistics tables */
   char* scip_name_;

   /** description of the statistics table */
   char* scip_desc_;

   /** position of the statistics table */
   const int scip_position_;

   /** output of the statistics table is only printed from this stage onwards */
   SCIP_STAGE scip_earlieststage_;

   /** default constructor */
   ObjTable(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of statistics table */
      const char*        desc,               /**< description of statistics table */
      int                position,           /**< position of statistics table */
      SCIP_STAGE         earlieststage       /**< output of the statistics table is only printed from this stage onwards */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_position_(position),
        scip_earlieststage_(earlieststage)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** copy constructor */
   ObjTable(const ObjTable& o) : ObjTable(o.scip_, o.scip_name_, o.scip_desc_, o.scip_position_, o.scip_earlieststage_)
   {
   }

   /** move constructor */
   ObjTable(ObjTable&& o)
       : scip_(o.scip_),
         scip_name_(0),
         scip_desc_(0),
         scip_position_(o.scip_position_),
         scip_earlieststage_(o.scip_earlieststage_)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
   }

   /** destructor */
   virtual ~ObjTable()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjTable& operator=(const ObjTable& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjTable& operator=(ObjTable&& o) = delete;

   /** destructor of statistics table to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_TABLEFREE(x) in @ref type_disp.h
    */
   virtual SCIP_DECL_TABLEFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** initialization method of statistics table (called after problem was transformed)
    *
    *  @see SCIP_DECL_TABLEINIT(x) in @ref type_table.h
    */
   virtual SCIP_DECL_TABLEINIT(scip_init)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** deinitialization method of statistics table (called before transformed problem is freed)
    *
    *  @see SCIP_DECL_TABLEEXIT(x) in @ref type_table.h
    */
   virtual SCIP_DECL_TABLEEXIT(scip_exit)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process initialization method of statistics table (called when branch and bound process is about to begin)
    *
    *  @see SCIP_DECL_TABLEINITSOL(x) in @ref type_table.h
    */
   virtual SCIP_DECL_TABLEINITSOL(scip_initsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** solving process deinitialization method of statistics table (called before branch and bound process data is freed)
    *
    *  @see SCIP_DECL_TABLEEXITSOL(x) in @ref type_table.h
    */
   virtual SCIP_DECL_TABLEEXITSOL(scip_exitsol)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** output method of statistics table to output file stream 'file'
    *
    *  @see SCIP_DECL_TABLEOUTPUT(x) in @ref type_table.h
    */
   virtual SCIP_DECL_TABLEOUTPUT(scip_output) = 0;
};

} /* namespace scip */



/** creates the statistics table for the given statistics table object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyTable* mytable = new MyTable(...);
 *       SCIP_CALL( SCIPincludeObjTable(scip, &mytable, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete mytable;    // delete table AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjTable(scip, new MyTable(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyTable is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjTable(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjTable*       objtable,           /**< statistics table object */
   SCIP_Bool             deleteobject        /**< should the statistics table object be deleted when statistics table is freed? */
   );

/** returns the statistics table object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjTable* SCIPfindObjTable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of statistics table */
   );

/** returns the statistics table object for the given statistics table */
SCIP_EXPORT
scip::ObjTable* SCIPgetObjTable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_TABLE*           table               /**< statistics table */
   );

#endif
