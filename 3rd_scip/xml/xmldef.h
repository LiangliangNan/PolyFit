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


#ifdef WITH_ZLIB
#include <zlib.h>

#define FOPEN(file, mode)    gzopen(file, mode)
#define FCLOSE(fp)           gzclose(fp)
#define FGETS(buf, len, fp)  gzgets(fp, buf, len)
#define FREAD(buf, len, fp)  gzread(fp, buf, len)
#define FPTYPE               gzFile
#else
#define FOPEN(file, mode)    fopen(file, mode)
#define FCLOSE(fp)           fclose(fp)
#define FGETS(buf, len, fp)  fgets(buf, len, fp)
#define FREAD(buf, len, fp)  fread(buf, 1, len, fp)
#define FPTYPE               FILE*
#endif /* WITH_ZLIB */


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
