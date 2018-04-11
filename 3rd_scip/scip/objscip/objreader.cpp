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

/**@file   objreader.cpp
 * @brief  C++ wrapper for file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objreader.h"




/*
 * Data structures
 */

/** file reader data */
struct SCIP_ReaderData
{
   scip::ObjReader*      objreader;          /**< file reader object */
   SCIP_Bool             deleteobject;       /**< should the reader object be deleted when reader is freed? */
};




/*
 * Callback methods of file reader
 */

extern "C"
{

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyObj)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;
   
   assert(scip != NULL);
   
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(readerdata->objreader != NULL);
   assert(readerdata->objreader->scip_ != scip);

   if( readerdata->objreader->iscloneable() )
   {
      scip::ObjReader* newobjreader;
      newobjreader = dynamic_cast<scip::ObjReader*> (readerdata->objreader->clone(scip));

      /* call include method of reader object */
      SCIP_CALL( SCIPincludeObjReader(scip, newobjreader, TRUE) );
   }

   return SCIP_OKAY;
}

/** destructor of file reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeObj)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(readerdata->objreader != NULL);
   assert(readerdata->objreader->scip_ == scip);

   /* call virtual method of reader object */
   SCIP_CALL( readerdata->objreader->scip_free(scip, reader) );

   /* free reader object */
   if( readerdata->deleteobject )
      delete readerdata->objreader;

   /* free reader data */
   delete readerdata;
   SCIPreaderSetData(reader, NULL); /*lint !e64*/
   
   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadObj)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(readerdata->objreader != NULL);
   assert(readerdata->objreader->scip_ == scip);

   /* call virtual method of reader object */
   SCIP_CALL( readerdata->objreader->scip_read(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteObj)
{  /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   assert(readerdata->objreader != NULL);

   /* call virtual method of reader object */
   SCIP_CALL( readerdata->objreader->scip_write(scip, reader, file, name, probdata, transformed, 
         objsense, objscale, objoffset, 
         vars, nvars, nbinvars, nintvars, nimplvars, ncontvars, fixedvars, nfixedvars, startnvars,
         conss, nconss, maxnconss, startnconss, genericnames, result) );
   
   return SCIP_OKAY;
}
}


/*
 * file reader specific interface methods
 */

/** creates the file reader for the given file reader object and includes it in SCIP */
SCIP_RETCODE SCIPincludeObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjReader*      objreader,          /**< file reader object */
   SCIP_Bool             deleteobject        /**< should the reader object be deleted when reader is freed? */
   )
{
   SCIP_READERDATA* readerdata;

   assert(scip != NULL);
   assert(objreader != NULL);

   /* create file reader data */
   readerdata = new SCIP_READERDATA;
   readerdata->objreader = objreader;
   readerdata->deleteobject = deleteobject;

   /* include file reader */
   SCIP_CALL( SCIPincludeReader(scip, objreader->scip_name_, objreader->scip_desc_, objreader->scip_extension_,
         readerCopyObj,
         readerFreeObj, readerReadObj, readerWriteObj, readerdata) ); /*lint !e429*/

   return SCIP_OKAY; /*lint !e429*/
}

/** returns the reader object of the given name, or 0 if not existing */
scip::ObjReader* SCIPfindObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of file reader */
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;

   reader = SCIPfindReader(scip, name);
   if( reader == NULL )
      return 0;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->objreader;
}
   
/** returns the reader object for the given file reader */
scip::ObjReader* SCIPgetObjReader(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader              /**< file reader */
   )
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   return readerdata->objreader;
}
