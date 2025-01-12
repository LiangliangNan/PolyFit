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

/**@file   scip_randnumgen.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for random numbers
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/misc.h"
#include "scip/pub_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_randnumgen.h"
#include "scip/set.h"
#include "scip/struct_scip.h"

/** creates and initializes a random number generator
 *
 *  @note The initial seed is changed using SCIPinitializeRandomSeed()
 */
SCIP_RETCODE SCIPcreateRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen,         /**< random number generator */
   unsigned int          initialseed,        /**< initial random seed */
   SCIP_Bool             useglobalseed       /**< should the supplied seed be initialized by SCIP's global seed shift? */
   )
{
   unsigned int modifiedseed;

   assert(scip != NULL);
   assert(randnumgen != NULL);

   if( useglobalseed )
      modifiedseed = SCIPinitializeRandomSeed(scip, initialseed);
   else
      modifiedseed = initialseed;

   SCIP_CALL( SCIPrandomCreate(randnumgen, SCIPblkmem(scip), modifiedseed) );

   return SCIP_OKAY;
}

/** frees a random number generator */
void SCIPfreeRandom(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN**     randnumgen          /**< random number generator */
   )
{
   assert(scip != NULL);
   assert(randnumgen != NULL);

   SCIPrandomFree(randnumgen, SCIPblkmem(scip));
}

/** initializes a random number generator with a given start seed
 *
 *  @note The seed is changed using SCIPinitializeRandomSeed()
 */
void SCIPsetRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   unsigned int          seed                /**< new random seed */
   )
{
   unsigned int modifiedseed;

   assert(scip != NULL);
   assert(randnumgen != NULL);

   modifiedseed = SCIPinitializeRandomSeed(scip, seed);

   SCIPrandomSetSeed(randnumgen, modifiedseed);
}

/** modifies an initial seed value with the global shift of random seeds */
unsigned int SCIPinitializeRandomSeed(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          initialseedvalue    /**< initial seed value to be modified */
   )
{
   assert(scip != NULL);

   return SCIPsetInitializeRandomSeed(scip->set, initialseedvalue);
}
