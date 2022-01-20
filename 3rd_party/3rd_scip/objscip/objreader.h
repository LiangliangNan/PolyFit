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

/**@file   objreader.h
 * @brief  C++ wrapper for file readers and writers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJREADER_H__
#define __SCIP_OBJREADER_H__

#include <cstring>

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

   /** destructor */
   virtual ~ObjReader()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
      SCIPfreeMemoryArray(scip_, &scip_extension_);
   }

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
EXTERN
SCIP_RETCODE SCIPincludeObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjReader*      objreader,          /**< file reader object */
   SCIP_Bool             deleteobject        /**< should the reader object be deleted when reader is freed? */
   );

/** returns the reader object of the given name, or 0 if not existing */
EXTERN
scip::ObjReader* SCIPfindObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of file reader */
   );

/** returns the reader object for the given file reader */
EXTERN
scip::ObjReader* SCIPgetObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader              /**< file reader */
   );

#endif
