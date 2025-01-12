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

/**@file   scip_bandit.c
 * @ingroup OTHER_CFILES
 * @brief  public functions for bandit algorithms
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/bandit.h"
#include "scip/pub_message.h"
#include "scip/scip_bandit.h"
#include "scip/scip_mem.h"
#include "scip/scip_randnumgen.h"
#include "scip/set.h"
#include "scip/struct_scip.h"
#include "scip/pub_message.h"
#include "scip/scip_bandit.h"
#include "scip/scip_mem.h"
#include "scip/scip_randnumgen.h"
#include "scip/set.h"
#include "scip/struct_scip.h"

/** includes a bandit algorithm virtual function table  */
SCIP_RETCODE SCIPincludeBanditvtable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDITVTABLE**   banditvtable,       /**< bandit algorithm virtual function table */
   const char*           name,               /**< a name for the algorithm represented by this vtable */
   SCIP_DECL_BANDITFREE  ((*banditfree)),    /**< callback to free bandit specific data structures */
   SCIP_DECL_BANDITSELECT((*banditselect)),  /**< selection callback for bandit selector */
   SCIP_DECL_BANDITUPDATE((*banditupdate)),  /**< update callback for bandit algorithms */
   SCIP_DECL_BANDITRESET ((*banditreset))    /**< update callback for bandit algorithms */
   )
{
   SCIP_BANDITVTABLE* vtableptr;

   assert(scip != NULL);
   assert(banditvtable != NULL);

   if( SCIPfindBanditvtable(scip, name) != NULL )
   {
      SCIPerrorMessage("bandit VTable <%s> already included.\n", name);
       return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbanditvtableCreate(&vtableptr, name,
         banditfree, banditselect, banditupdate, banditreset) );

   SCIP_CALL( SCIPsetIncludeBanditvtable(scip->set, vtableptr) );

   *banditvtable = vtableptr;

   return SCIP_OKAY;
}

/** returns the bandit virtual function table of the given name, or NULL if not existing */
SCIP_BANDITVTABLE* SCIPfindBanditvtable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of bandit algorithm virtual function table */
   )
{
   assert(scip != NULL);

   return SCIPsetFindBanditvtable(scip->set, name);
}

/** reset the bandit algorithm */
SCIP_RETCODE SCIPresetBandit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT*          bandit,             /**< pointer to bandit algorithm data structure */
   SCIP_Real*            priorities,         /**< priorities for every action, or NULL if not needed */
   unsigned int          seed                /**< initial random seed for bandit selection */
   )
{
   assert(scip != NULL);
   assert(bandit != NULL);

   SCIP_CALL( SCIPbanditReset(SCIPbuffer(scip), bandit, priorities, SCIPinitializeRandomSeed(scip, seed)) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of bandit algorithm */
SCIP_RETCODE SCIPfreeBandit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         bandit              /**< pointer to bandit algorithm data structure */
   )
{
   assert(scip != NULL);
   assert(bandit != NULL);
   assert(*bandit != NULL);

   SCIP_CALL( SCIPbanditFree(SCIPblkmem(scip), bandit) );

   return SCIP_OKAY;
}
