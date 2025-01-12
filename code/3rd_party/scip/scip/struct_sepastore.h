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

/**@file   struct_sepastore.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for storing separated cuts
 * @author Tobias Achterberg
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SEPASTORE_H__
#define __SCIP_STRUCT_SEPASTORE_H__


#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_sepastore.h"

#ifdef __cplusplus
extern "C" {
#endif

/** storage for separated cuts */
struct SCIP_SepaStore
{
   SCIP_ROW**            cuts;               /**< array with separated cuts sorted by score */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator used for tie breaking */
   int                   cutssize;           /**< size of cuts and score arrays */
   int                   ncuts;              /**< number of separated cuts (max. is set->sepa_maxcuts) */
   int                   nforcedcuts;        /**< number of forced separated cuts (first positions in cuts array) */
   int                   ncutsadded;         /**< total number of cuts added so far */
   int                   ncutsaddedviapool;  /**< total number of cuts added from cutpool */
   int                   ncutsaddeddirect;   /**< total number of cuts added directly */
   int                   ncutsfoundround;    /**< number of cuts found so far in this separation round */
   int                   ncutsapplied;       /**< total number of cuts applied to the LP */
   SCIP_Bool             initiallp;          /**< is the separation storage currently being filled with the initial LP rows? */
   SCIP_Bool             forcecuts;          /**< should the cuts be used despite the number of cuts parameter limit? */
};

#ifdef __cplusplus
}
#endif

#endif
