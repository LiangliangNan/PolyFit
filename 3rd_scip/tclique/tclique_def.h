/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*                    TCLIQUE --- Algorithm for Maximum Cliques              */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  TCLIQUE is distributed under the terms of the ZIB Academic License.      */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with TCLIQUE; see the file COPYING.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   tclique_def.h
 * @brief  tclique defines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_DEF_H__
#define __TCLIQUE_DEF_H__

#ifdef WITH_SCIPDEF
#include "scip/def.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

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

#ifndef debug
#ifdef TCLIQUE_DEBUG
#define debug(x)                        x
#define debugMessage                    printf("[%s:%d] debug: ", __FILE__, __LINE__); printf
#define debugPrintf                     printf
#else
#define debug(x)                        /**/
#define debugMessage                    while( FALSE ) printf
#define debugPrintf                     while( FALSE ) printf
#endif
#endif

#ifndef infoMessage
#define infoMessage printf
#endif

#ifndef MAX
#define MAX(x,y)  ((x) >= (y) ? (x) : (y))
#endif

#ifdef __cplusplus
}
#endif

#endif
