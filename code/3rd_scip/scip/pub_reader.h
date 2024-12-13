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

/**@file   pub_reader.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for input file readers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_READER_H__
#define __SCIP_PUB_READER_H__


#include "scip/def.h"
#include "scip/type_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicReaderMethods
 *
 * @{
 */

/** gets user data of reader */
SCIP_EXPORT
SCIP_READERDATA* SCIPreaderGetData(
   SCIP_READER*          reader              /**< reader */
   );

/** sets user data of reader; user has to free old data in advance! */
SCIP_EXPORT
void SCIPreaderSetData(
   SCIP_READER*          reader,             /**< reader */
   SCIP_READERDATA*      readerdata          /**< new reader user data */
   );

/** gets name of reader */
SCIP_EXPORT
const char* SCIPreaderGetName(
   SCIP_READER*          reader              /**< reader */
   );

/** gets description of reader */
SCIP_EXPORT
const char* SCIPreaderGetDesc(
   SCIP_READER*          reader              /**< reader */
   );

/** gets file extension of reader */
SCIP_EXPORT
const char* SCIPreaderGetExtension(
   SCIP_READER*          reader              /**< reader */
   );

/** return whether the reader can read files */
SCIP_EXPORT
SCIP_Bool SCIPreaderCanRead(
   SCIP_READER*          reader              /**< reader */
   );

/** return whether the reader can write files */
SCIP_EXPORT
SCIP_Bool SCIPreaderCanWrite(
   SCIP_READER*          reader              /**< reader */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
