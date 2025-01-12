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

/**@file   scip_randnumgen.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for random numbers
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_RANDNUMGEN_H__
#define __SCIP_SCIP_RANDNUMGEN_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup RandomNumbers
 *
 * @{
 */

/** creates and initializes a random number generator
 *
 *  @note The initial seed is changed using SCIPinitializeRandomSeed()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   unsigned int          initialseed,        /**< initial random seed */
   SCIP_Bool             useglobalseed       /**< should SCIP's global seed be used to initialise the supplied seed? */
   );

/** frees a random number generator */
SCIP_EXPORT
void SCIPfreeRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen          /**< random number generator */
   );

/** initializes a random number generator with a given seed
 *
 *  @note The seed is changed using SCIPinitializeRandomSeed()
 */
SCIP_EXPORT
void SCIPsetRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   unsigned int          seed                /**< new random seed */
   );


/** modifies an initial seed value with the global shift of random seeds */
SCIP_EXPORT
unsigned int SCIPinitializeRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          initialseedvalue    /**< initial seed value to be modified */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
