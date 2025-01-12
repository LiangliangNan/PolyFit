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

/**@file   xmldef.h
 * @brief  definitions for XML parsing
 * @author Thorsten Koch
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_XMLDEF_H__
#define __SCIP_XMLDEF_H__

#ifdef __cplusplus
extern "C" {
#endif


#ifndef XML_Bool
#define XML_Bool unsigned int           /**< type used for boolean values */
#endif

#ifndef TRUE
#define TRUE  1                         /**< boolean value TRUE */
#define FALSE 0                         /**< boolean value FALSE */
#endif


#ifdef SCIP_WITH_ZLIB
#include <zlib.h>

#define FOPEN(file, mode)    gzopen(file, mode)
#define FCLOSE(fp)           gzclose(fp)
#define FGETS(buf, len, fp)  gzgets(fp, buf, len) /*lint !e755 */
#define FREAD(buf, len, fp)  gzread(fp, buf, len)
#define FPTYPE               gzFile
#else
#define FOPEN(file, mode)    fopen(file, mode)
#define FCLOSE(fp)           fclose(fp)
#define FGETS(buf, len, fp)  fgets(buf, len, fp) /*lint !e755 */
#define FREAD(buf, len, fp)  fread(buf, 1, len, fp)
#define FPTYPE               FILE*
#endif /* SCIP_WITH_ZLIB */


#ifndef ALLOC_ABORT
#define ALLOC_ABORT(x) do                                               \
   {                                                                    \
      if( NULL == (x) )                                                 \
      {                                                                 \
         printf("[%s:%d] No memory in function call\n", __FILE__, __LINE__); \
         abort();                                                       \
      }                                                                 \
   }                                                                    \
   while( FALSE )
#endif

#ifndef ALLOC_FALSE
#define ALLOC_FALSE(x)  do                                              \
   {                                                                    \
      if( NULL == (x) )                                                 \
      {                                                                 \
         printf("[%s:%d] No memory in function call\n", __FILE__, __LINE__); \
         return FALSE;                                                  \
      }                                                                 \
   }                                                                    \
   while( FALSE )
#endif

#ifdef XML_DEBUG
#define debug(x)                        x
#define debugMessage                    printf("[%s:%d] debug: ", __FILE__, __LINE__); printf
#define debugPrintf                     printf
#else
#define debug(x)                        /**/
#define debugMessage                    while( FALSE ) printf
#define debugPrintf                     while( FALSE ) printf
#endif

#ifndef infoMessage
#define infoMessage printf
#endif

#ifdef __cplusplus
}
#endif

#endif
