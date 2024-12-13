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

/**@file   objreader.h
 * @brief  C++ wrapper for file readers and writers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJREADER_H__
#define __SCIP_OBJREADER_H__

#include <cstring>
#include <utility>

#include "scip/scip.h"
#include "objscip/objcloneable.h"

namespace scip
{

/** @brief C++ wrapper for file readers and writers
 *
 *  This class defines the interface for file readers and writers implemented in C++.
 *
 *  - \ref READER "Instructions for implementing a file reader and writer"
 *  - \ref FILEREADERS "List of available file readers and writers"
 *  - \ref type_reader.h "Corresponding C interface"
 */
class ObjReader : public ObjCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the file reader */
   char* scip_name_;

   /** description of the file reader */
   char* scip_desc_;

   /** file extension that reader processes */
   char* scip_extension_;

   /** default constructor */
   ObjReader(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of file reader */
      const char*        desc,               /**< description of file reader */
      const char*        extension           /**< file extension that reader processes */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_extension_(0)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_extension_, extension, std::strlen(extension)+1) );
   }

   /** copy constructor */
   ObjReader(const ObjReader& o) : ObjReader(o.scip_, o.scip_name_, o.scip_desc_, o.scip_extension_) {}

   /** move constructor */
   ObjReader(ObjReader&& o) : scip_(o.scip_), scip_name_(0), scip_desc_(0), scip_extension_(0)
   {
      std::swap(scip_name_, o.scip_name_);
      std::swap(scip_desc_, o.scip_desc_);
      std::swap(scip_extension_, o.scip_extension_);
   }

   /** destructor */
   virtual ~ObjReader()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
      SCIPfreeMemoryArray(scip_, &scip_extension_);
   }

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjReader& operator=(const ObjReader& o) = delete;

   /** assignment of polymorphic classes causes slicing and is therefore disabled. */
   ObjReader& operator=(ObjReader&& o) = delete;

   /** destructor of file reader to free user data (called when SCIP is exiting)
    *
    *  @see SCIP_DECL_READERFREE(x) in @ref type_reader.h
    */
   virtual SCIP_DECL_READERFREE(scip_free)
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** problem reading method of reader
    *
    *  @see SCIP_DECL_READERREAD(x) in @ref type_reader.h
    */
   virtual SCIP_DECL_READERREAD(scip_read)
   {  /*lint --e{715}*/

      /* set result pointer to indicate that the reading was not performed */
      assert(result != NULL);
      (*result) = SCIP_DIDNOTRUN;

      return SCIP_OKAY;
   }

   /** problem writing method of reader; NOTE: if the parameter "genericnames" is TRUE, then
    *  SCIP already set all variable and constraint names to generic names; therefore, this
    *  method should always use SCIPvarGetName() and SCIPconsGetName(); 
    *
    *  @see SCIP_DECL_READERWRITE(x) in @ref type_reader.h
    */
   virtual SCIP_DECL_READERWRITE(scip_write)
   {  /*lint --e{715}*/

      /* set result pointer to indicate that the writing was not performed */
      assert(result != NULL);
      (*result) = SCIP_DIDNOTRUN;

      return SCIP_OKAY;
   }
};

} /* namespace scip */



/** creates the file reader for the given file reader object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyReader* myreader = new MyReader(...);
 *       SCIP_CALL( SCIPincludeObjReader(scip, &myreader, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myreader;    // delete reader AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjReader(scip, new MyReader(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyReader is called here
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjReader*      objreader,          /**< file reader object */
   SCIP_Bool             deleteobject        /**< should the reader object be deleted when reader is freed? */
   );

/** returns the reader object of the given name, or 0 if not existing */
SCIP_EXPORT
scip::ObjReader* SCIPfindObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of file reader */
   );

/** returns the reader object for the given file reader */
SCIP_EXPORT
scip::ObjReader* SCIPgetObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader              /**< file reader */
   );

#endif
