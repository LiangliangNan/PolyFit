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

/**@file   type_visual.h
 * @brief  type definitions for output for visualization tools (VBC, BAK)
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_VISUAL_H__
#define __SCIP_TYPE_VISUAL_H__

#ifdef __cplusplus
extern "C" {
#endif

/** node colors in VBC output:
 *   1: indian red
 *   2: green
 *   3: light gray
 *   4: red
 *   5: blue
 *   6: black
 *   7: light pink
 *   8: cyan
 *   9: dark green
 *  10: brown
 *  11: orange
 *  12: yellow
 *  13: pink
 *  14: purple
 *  15: light blue
 *  16: muddy green
 *  17: white
 *  18: light grey
 *  19: light grey
 *  20: light grey
 */
enum SCIP_VBCColor
{
   SCIP_VBCCOLOR_UNSOLVED   =  3,       /**< color for newly created, unsolved nodes */
   SCIP_VBCCOLOR_SOLVED     =  2,       /**< color for solved nodes */
   SCIP_VBCCOLOR_CUTOFF     =  4,       /**< color for nodes that were cut off */
   SCIP_VBCCOLOR_CONFLICT   = 15,       /**< color for nodes where a conflict constraint was found */
   SCIP_VBCCOLOR_MARKREPROP = 11,       /**< color for nodes that were marked to be repropagated */
   SCIP_VBCCOLOR_REPROP     = 12,       /**< color for repropagated nodes */
   SCIP_VBCCOLOR_SOLUTION   = 14,       /**< color for solved nodes, where a solution has been found */
   SCIP_VBCCOLOR_NONE       = -1        /**< color should not be changed */
};
typedef enum SCIP_VBCColor SCIP_VBCCOLOR;


typedef struct SCIP_Visual SCIP_VISUAL;      /**< VBC Tool data structure */

#ifdef __cplusplus
}
#endif

#endif
