/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*              TCLIQUE --- Algorithm for Maximum Cliques                    */
/*                                                                           */
/*  Copyright 1996-2022 Zuse Institute Berlin                                */
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
/*  along with TCLIQUE; see the file LICENSE.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   tclique_def.h
 * @brief  tclique defines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TCLIQUE_DEF_H__
#define __TCLIQUE_DEF_H__

/*
 * include build configuration flags
 */
#ifndef NO_CONFIG_HEADER
#include "scip/config.h"
#include "scip/scip_export.h"
#endif

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
