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

/**@file   type_reopt.h
 * @brief  type definitions for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_REOPT_H__
#define __SCIP_TYPE_REOPT_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Reopt SCIP_REOPT;             /**< reopt data */

typedef struct SCIP_SolTree SCIP_SOLTREE;         /**< tree to check solutions */

typedef struct SCIP_SolNode SCIP_SOLNODE;         /**< nodes of SCIP_SOLTREE */

typedef struct SCIP_ReoptTree SCIP_REOPTTREE;     /**< data structure to store the search tree */

typedef struct SCIP_ReoptNode SCIP_REOPTNODE;     /**< nodes of SCIP_REOPTTREE */

typedef struct SCIP_ReoptNode SCIP_REPRESENTATIVE;/**< representatives of the search frontier */

typedef struct SCIP_ReoptConsData SCIP_REOPTCONSDATA; /**< data for constraints to handle dual information \
                                                        *  within (mixed) binary programs
                                                        */

/* type of nodes during reoptimization */
enum SCIP_ReoptType
{
   SCIP_REOPTTYPE_NONE        = 0,                /**< node is not part of the reoptimizationtree */
   SCIP_REOPTTYPE_TRANSIT     = 1,                /**< node is only needed for reconstructing the tree */
   SCIP_REOPTTYPE_INFSUBTREE  = 2,                /**< node contains dual reductions which leed to LP infeasibility */
   SCIP_REOPTTYPE_STRBRANCHED = 3,                /**< node contains dual reductions */
   SCIP_REOPTTYPE_LOGICORNODE = 4,                /**< node contains additional constraints */
   SCIP_REOPTTYPE_LEAF        = 5,                /**< node is a leaf node */
   SCIP_REOPTTYPE_PRUNED      = 6,                /**< node is a leaf node and pruned by boudning */
   SCIP_REOPTTYPE_FEASIBLE    = 7                 /**< node is a leaf node and has an integral optimal LP solution */
};
typedef enum SCIP_ReoptType SCIP_REOPTTYPE;       /**< type nodes during reoptimization */

enum Reopt_ConsType
{
   REOPT_CONSTYPE_INFSUBTREE   = 0,               /**< constraint cutoffs an LP infeasible subtree */
   REOPT_CONSTYPE_DUALREDS     = 1,               /**< constraint reconstructs dual reductions */
   REOPT_CONSTYPE_CUT          = 2,               /**< constraint representing a cut, e.g., to separate a solution */
   REOPT_CONSTYPE_UNKNOWN      = 3                /**< constraint was added by SCIP, e.g., a (local) conflict */
};
typedef enum Reopt_ConsType REOPT_CONSTYPE;       /**< tye of constraunts added during reoptimization */

#ifdef __cplusplus
}
#endif

#endif
