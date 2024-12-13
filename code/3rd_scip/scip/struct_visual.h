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

/**@file   struct_visual.h
 * @ingroup INTERNALAPI
 * @brief  data structures for output for visualization tools (VBC, BAK)
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_VISUAL_H__
#define __SCIP_STRUCT_VISUAL_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_visual.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** visual data structure */
struct SCIP_Visual
{
   FILE*                 vbcfile;            /**< file to store VBC information */
   FILE*                 bakfile;            /**< file to store BAK information */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< message handler to use */
   SCIP_HASHMAP*         nodenum;            /**< hash map for mapping nodes to node numbers */
   SCIP_Longint          timestep;           /**< time step counter for non real time output */
   SCIP_NODE*            lastnode;           /**< last node that was colored */
   SCIP_VBCCOLOR         lastcolor;          /**< last color that was used */
   SCIP_Bool             userealtime;        /**< should the real solving time be used instead of a time step counter? */
   SCIP_Real             lastlowerbound;     /**< last lower bound that was output */
};

#ifdef __cplusplus
}
#endif

#endif
