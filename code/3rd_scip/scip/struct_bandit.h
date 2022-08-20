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

/**@file   struct_bandit.h
 * @ingroup INTERNALAPI
 * @brief  data structures for bandit selection algorithms
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BANDIT_H__
#define __SCIP_STRUCT_BANDIT_H__


#include "scip/def.h"
#include "misc.h"
#include "scip/type_clock.h"
#include "scip/type_bandit.h"

#ifdef __cplusplus
extern "C" {
#endif

/** virtual function table for bandit selection algorithms */
struct SCIP_BanditVTable
{
   const char*           name;               /**< name of the represented bandit algorithm */
   SCIP_DECL_BANDITFREE  ((*banditfree));    /**< callback to free bandit specific data structures */
   SCIP_DECL_BANDITSELECT((*banditselect));  /**< selection callback for bandit selector */
   SCIP_DECL_BANDITUPDATE((*banditupdate));  /**< update callback for bandit algorithms */
   SCIP_DECL_BANDITRESET ((*banditreset));   /**< update callback for bandit algorithms */
};

/** data structure for bandit algorithms */
struct SCIP_Bandit
{
   SCIP_BANDITVTABLE*    vtable;             /**< virtual function table for callbacks */
   SCIP_RANDNUMGEN*      rng;                /**< random number generator for randomized selection */
   int                   nactions;           /**< the number of actions to select from */
   SCIP_BANDITDATA*      data;               /**< specific data for bandit algorithm implementations */
};
#ifdef __cplusplus
}
#endif

#endif
