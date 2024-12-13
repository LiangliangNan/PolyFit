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

/**@file   pub_fileio.h
 * @ingroup PUBLICCOREAPI
 * @brief  wrapper functions to map file i/o to standard or zlib file i/o
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_FILEIO_H__
#define __SCIP_PUB_FILEIO_H__

#include <stddef.h>
#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_File SCIP_FILE;          /**< file data structure */

SCIP_EXPORT SCIP_FILE* SCIPfopen(const char *path, const char *mode);
SCIP_EXPORT SCIP_FILE* SCIPfdopen(int fildes, const char *mode);
SCIP_EXPORT size_t SCIPfread(void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream);
SCIP_EXPORT size_t SCIPfwrite(const void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream);
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
SCIP_EXPORT int SCIPfprintf(SCIP_FILE *stream, const char *format, ...);
SCIP_EXPORT int SCIPfputc(int c, SCIP_FILE *stream);
SCIP_EXPORT int SCIPfputs(const char *s, SCIP_FILE *stream);
SCIP_EXPORT int SCIPfgetc(SCIP_FILE *stream);
SCIP_EXPORT char* SCIPfgets(char *s, int size, SCIP_FILE *stream);
SCIP_EXPORT int SCIPfflush(SCIP_FILE *stream);
SCIP_EXPORT int SCIPfseek(SCIP_FILE *stream, long offset, int whence);
SCIP_EXPORT void SCIPrewind(SCIP_FILE *stream);
SCIP_EXPORT long SCIPftell(SCIP_FILE *stream);
SCIP_EXPORT int SCIPfeof(SCIP_FILE *stream);
SCIP_EXPORT int SCIPfclose(SCIP_FILE *fp);

#ifdef __cplusplus
}
#endif

#endif
