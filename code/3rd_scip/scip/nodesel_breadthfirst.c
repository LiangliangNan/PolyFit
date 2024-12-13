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

/**@file nodesel_breadthfirst.h
 * @ingroup DEFPLUGINS_NODESEL
 * @ingroup NODESELECTORS
 * @brief node selector for breadth-first search
 * @author Stefan Heinz
 * @author Gregor Hendel
 *
 * This node selector performs breadth-first search, i.e., it completely evaluates an entire level of the search tree before
 * proceeding to the next level. At one level, nodes are processed in the order they were created by the branching rule.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/nodesel_breadthfirst.h"
#include "scip/pub_message.h"
#include "scip/pub_nodesel.h"
#include "scip/pub_tree.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_tree.h"
#include <string.h>

#define NODESEL_NAME             "breadthfirst"
#define NODESEL_DESC             "breadth first search"
#define NODESEL_STDPRIORITY        -10000
#define NODESEL_MEMSAVEPRIORITY  -1000000

/*
 * Callback methods
 */

/** copy method for node selector plugins (called when SCIP copies plugins) */
static
SCIP_DECL_NODESELCOPY(nodeselCopyBreadthfirst)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);

   /* call inclusion method of node selector */
   SCIP_CALL( SCIPincludeNodeselBreadthfirst(scip) );

   return SCIP_OKAY;
}

/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectBreadthfirst)
{  /*lint --e{715}*/
   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   /* siblings come before leaves at the same level. Sometimes it can occur that no leaves are left except for children */
   *selnode = SCIPgetBestSibling(scip);
   if( *selnode == NULL )
   {
      *selnode = SCIPgetBestLeaf(scip);
      if( *selnode == NULL )
         *selnode=SCIPgetBestChild(scip);
   }
   if( *selnode != NULL )
   {
      SCIPdebugMsg(scip, "Selecting next node number %" SCIP_LONGINT_FORMAT " at depth %d\n", SCIPnodeGetNumber(*selnode), SCIPnodeGetDepth(*selnode));
   }

   return SCIP_OKAY;
}


/** node comparison method of breadth first search: nodes with lower depth are preferred; in case of a tie, the node
 *  which was created earlier (and therefore has a smaller node number) is preferred */
static
SCIP_DECL_NODESELCOMP(nodeselCompBreadthfirst)
{  /*lint --e{715}*/
   int depth1;
   int depth2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   depth1 = SCIPnodeGetDepth(node1);
   depth2 = SCIPnodeGetDepth(node2);

   /* if depths differ, prefer node with smaller depth */
   if( depth1 < depth2 )
      return -1;
   else if( depth1 > depth2 )
      return +1;
   else
   {
      /* depths are equal; prefer node with smaller number */
      SCIP_Longint number1;
      SCIP_Longint number2;

      number1 = SCIPnodeGetNumber(node1);
      number2 = SCIPnodeGetNumber(node2);
      assert(number1 != number2);

      if( number1 < number2 )
         return -1;
      else
         return +1;
   }
}

/*
 * breadth first specific interface methods
 */

/** creates the node selector for breadth first search and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselBreadthfirst(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESEL* nodesel;

   /* include node selector */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
         nodeselSelectBreadthfirst, nodeselCompBreadthfirst, NULL) );

   assert(nodesel != NULL);

   /* set non-fundamental callback functions via setter functions */
   SCIP_CALL ( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyBreadthfirst) );

   return SCIP_OKAY;
}
