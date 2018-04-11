/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
};

#ifdef __cplusplus
}
#endif

#endif
