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
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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

EXTERN SCIP_FILE* SCIPfopen(const char *path, const char *mode);
EXTERN SCIP_FILE* SCIPfdopen(int fildes, const char *mode);
EXTERN size_t SCIPfread(void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream);
EXTERN size_t SCIPfwrite(const void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream);
EXTERN int SCIPfprintf(SCIP_FILE *stream, const char *format, ...);
EXTERN int SCIPfputc(int c, SCIP_FILE *stream);
EXTERN int SCIPfputs(const char *s, SCIP_FILE *stream);
EXTERN int SCIPfgetc(SCIP_FILE *stream);
EXTERN char* SCIPfgets(char *s, int size, SCIP_FILE *stream);
EXTERN int SCIPfflush(SCIP_FILE *stream);
EXTERN int SCIPfseek(SCIP_FILE *stream, long offset, int whence);
EXTERN void SCIPrewind(SCIP_FILE *stream);
EXTERN long SCIPftell(SCIP_FILE *stream);
EXTERN int SCIPfeof(SCIP_FILE *stream);
EXTERN int SCIPfclose(SCIP_FILE *fp);

#ifdef __cplusplus
}
#endif

#endif
