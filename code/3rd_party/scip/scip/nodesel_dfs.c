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

/**@file   nodesel_dfs.c
 * @ingroup DEFPLUGINS_NODESEL
 * @brief  node selector for depth first search
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nodesel_dfs.h"
#include "scip/pub_message.h"
#include "scip/pub_nodesel.h"
#include "scip/pub_tree.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_tree.h"
#include <string.h>

#define NODESEL_NAME             "dfs"
#define NODESEL_DESC             "depth first search"
#define NODESEL_STDPRIORITY           0
#define NODESEL_MEMSAVEPRIORITY  100000


/*
 * Callback methods
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyDfs)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );

   return SCIP_OKAY;
}


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectDfs)
{  /*lint --e{715}*/
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   *selnode = SCIPgetPrioChild(scip);
   if( *selnode == NULL )
   {
      *selnode = SCIPgetPrioSibling(scip);
      if( *selnode == NULL )
      {
         SCIPdebugMsg(scip, "select best leaf\n");
         *selnode = SCIPgetBestLeaf(scip);
      }

      SCIPdebugMsg(scip, "select best sibling leaf\n");
   }

   return SCIP_OKAY;
}


/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompDfs)
{  /*lint --e{715}*/
   int depth1;
   int depth2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   depth1 = SCIPnodeGetDepth(node1);
   depth2 = SCIPnodeGetDepth(node2);
   if( depth1 > depth2 )
      return -1;
   else if( depth1 < depth2 )
      return +1;
   else
   {
      SCIP_Real lowerbound1;
      SCIP_Real lowerbound2;

      lowerbound1 = SCIPnodeGetLowerbound(node1);
      lowerbound2 = SCIPnodeGetLowerbound(node2);
      if( lowerbound1 < lowerbound2 )
         return -1;
      else if( lowerbound1 > lowerbound2 )
         return +1;
      else
         return 0;
   }
}


/*
 * dfs specific interface methods
 */

/** creates the node selector for depth first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselDfs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESEL* nodesel;

   /* include node selector */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
          nodeselSelectDfs, nodeselCompDfs, NULL) );

   assert(nodesel != NULL);

   SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyDfs) );

   return SCIP_OKAY;
}
