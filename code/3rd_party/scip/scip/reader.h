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

/**@file   reader.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_H__
#define __SCIP_READER_H__


#include "scip/def.h"
#include "scip/type_prob.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_reader.h"
#include "scip/pub_reader.h"

#ifdef __cplusplus
extern "C" {
#endif


/** copies the given reader to a new scip */
SCIP_RETCODE SCIPreaderCopyInclude(
   SCIP_READER*          reader,             /**< reader */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a reader */
SCIP_RETCODE SCIPreaderCreate(
   SCIP_READER**         reader,             /**< pointer to store reader */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of reader */
   const char*           desc,               /**< description of reader */
   const char*           extension,          /**< file extension that reader processes */
   SCIP_DECL_READERCOPY  ((*readercopy)),    /**< copy method of reader or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   SCIP_DECL_READERREAD  ((*readerread)),    /**< read method */
   SCIP_DECL_READERWRITE ((*readerwrite)),   /**< write method */
   SCIP_READERDATA*      readerdata          /**< reader data */
   );

/** frees memory of reader */
SCIP_RETCODE SCIPreaderFree(
   SCIP_READER**         reader,             /**< pointer to reader data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** reads problem data from file with given reader or returns SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPreaderRead(
   SCIP_READER*          reader,             /**< reader */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           filename,           /**< name of the input file */
   const char*           extension,          /**< extension of the input file name */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** writes problem data to file with given reader or returns SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPreaderWrite(
   SCIP_READER*          reader,             /**< reader */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           format,             /**< file format (or NULL) */
   SCIP_Bool             genericnames,       /**< using generic variable and constraint names? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** gets time in seconds used in this reader for reading */
SCIP_Real SCIPreaderGetReadingTime(
   SCIP_READER*          reader              /**< reader */
   );

/** enables or disables all clocks of \p reader, depending on the value of the flag */
void SCIPreaderEnableOrDisableClocks(
   SCIP_READER*          reader,             /**< the reader for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks be enabled? */
   );

/** resets reading time of reader */
SCIP_RETCODE SCIPreaderResetReadingTime(
   SCIP_READER*          reader              /**< reader */
   );

/** sets copy method of reader */
void SCIPreaderSetCopy(
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERCOPY  ((*readercopy))     /**< copy method of reader or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor of reader */
void SCIPreaderSetFree(
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERFREE  ((*readerfree))     /**< destructor of reader */
   );

/** sets read method of reader */
void SCIPreaderSetRead(
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERREAD  ((*readerread))     /**< read method */
   );

/** sets write method of reader */
void SCIPreaderSetWrite(
   SCIP_READER*          reader,             /**< reader */
   SCIP_DECL_READERWRITE ((*readerwrite))    /**< write method */
   );

#ifdef __cplusplus
}
#endif

#endif
