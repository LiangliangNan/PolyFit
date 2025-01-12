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

/**@file   scip_reader.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for reader plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/reader.h"
#include "scip/scip_reader.h"
#include "scip/set.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"

/** creates a reader and includes it in SCIP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note method has all reader callbacks as arguments and is thus changed every time a new callback is added
 *        in future releases; consider using SCIPincludeReaderBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of reader */
   const char*           desc,               /**< description of reader */
   const char*           extension,          /**< file extension that reader processes */
   SCIP_DECL_READERCOPY  ((*readercopy)),    /**< copy method of reader or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   SCIP_DECL_READERREAD  ((*readerread)),    /**< read method */
   SCIP_DECL_READERWRITE ((*readerwrite)),   /**< write method */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   SCIP_READER* reader;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeReader", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether reader is already present */
   if( SCIPfindReader(scip, name) != NULL )
   {
      SCIPerrorMessage("reader <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPreaderCreate(&reader, scip->set, name, desc, extension, readercopy, readerfree, readerread,
      readerwrite, readerdata) );
   SCIP_CALL( SCIPsetIncludeReader(scip->set, reader) );

   return SCIP_OKAY;
}

/** creates a reader and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see
 *  SCIPsetReaderCopy(), SCIPsetReaderFree(), SCIPsetReaderRead(), SCIPsetReaderWrite().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeReader() instead
 */
SCIP_RETCODE SCIPincludeReaderBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER**         readerptr,          /**< reference to reader pointer, or NULL */
   const char*           name,               /**< name of reader */
   const char*           desc,               /**< description of reader */
   const char*           extension,          /**< file extension that reader processes */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   SCIP_READER* reader;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeReaderBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether reader is already present */
   if( SCIPfindReader(scip, name) != NULL )
   {
      SCIPerrorMessage("reader <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPreaderCreate(&reader, scip->set, name, desc, extension, NULL, NULL, NULL, NULL, readerdata) );
   SCIP_CALL( SCIPsetIncludeReader(scip->set, reader) );

   if( readerptr != NULL )
      *readerptr = reader;

   return SCIP_OKAY;
}

/** set copy method of reader
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetReaderCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERCOPY  ((*readercopy))     /**< copy method of reader or NULL if you don't want to copy your plugin into sub-SCIPs */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetReaderCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPreaderSetCopy(reader, readercopy);

   return SCIP_OKAY;
}

/** set deinitialization method of reader
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetReaderFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERFREE  ((*readerfree))     /**< destructor of reader */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetReaderFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPreaderSetFree(reader, readerfree);

   return SCIP_OKAY;
}

/** set read method of reader
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetReaderRead(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERREAD  ((*readerread))     /**< read method of reader */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetReaderRead", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPreaderSetRead(reader, readerread);

   return SCIP_OKAY;
}

/** set write method of reader
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetReaderWrite(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERWRITE ((*readerwrite))    /**< write method of reader */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetReaderWrite", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPreaderSetWrite(reader, readerwrite);

   return SCIP_OKAY;
}

/** returns the reader of the given name, or NULL if not existing */
SCIP_READER* SCIPfindReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of reader */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindReader(scip->set, name);
}

/** returns the array of currently available readers */
SCIP_READER** SCIPgetReaders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->readers;
}

/** returns the number of currently available readers */
int SCIPgetNReaders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nreaders;
}
