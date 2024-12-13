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

/**@file   type_sepastore.h
 * @brief  type definitions for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SEPASTORE_H__
#define __SCIP_TYPE_SEPASTORE_H__

#ifdef __cplusplus
extern "C" {
#endif

/** possible settings for specifying the solution for which cuts are selected */
enum SCIP_Efficiacychoice
{
   SCIP_EFFICIACYCHOICE_LP    = 0,           /**< use LP solution to base efficacy on */
   SCIP_EFFICIACYCHOICE_RELAX = 1,           /**< use relaxation solution to base efficacy on */
   SCIP_EFFICIACYCHOICE_NLP   = 2            /**< use NLP solution to base efficacy on */
};
typedef enum SCIP_Efficiacychoice SCIP_EFFICIACYCHOICE;

typedef struct SCIP_SepaStore SCIP_SEPASTORE;     /**< storage for separated variables */

#ifdef __cplusplus
}
#endif

#endif
