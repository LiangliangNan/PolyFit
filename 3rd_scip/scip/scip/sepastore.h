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

/**@file   sepastore.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing separated cuts
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPASTORE_H__
#define __SCIP_SEPASTORE_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_implics.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_sepastore.h"
#include "scip/type_branch.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates separation storage */
extern
SCIP_RETCODE SCIPsepastoreCreate(
   SCIP_SEPASTORE**      sepastore           /**< pointer to store separation storage */
   );

/** frees separation storage */
extern
SCIP_RETCODE SCIPsepastoreFree(
   SCIP_SEPASTORE**      sepastore           /**< pointer to store separation storage */
   );

/** informs separation storage, that the setup of the initial LP starts now */
extern
void SCIPsepastoreStartInitialLP(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** informs separation storage, that the setup of the initial LP is now finished */
extern
void SCIPsepastoreEndInitialLP(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** informs separation storage, that the following cuts should be used in any case */
extern
void SCIPsepastoreStartForceCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** informs separation storage, that the following cuts should no longer be used in any case */
extern
void SCIPsepastoreEndForceCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** adds cut to separation storage and captures it;
 *  if the cut should be forced to enter the LP, an infinite score has to be used
 */
extern
SCIP_RETCODE SCIPsepastoreAddCut(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut,           /**< should the cut be forced to enter the LP? */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_Bool*            infeasible          /**< pointer to store whether the cut is infeasible */
   );

/** adds cuts to the LP and clears separation storage */
extern
SCIP_RETCODE SCIPsepastoreApplyCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_EFFICIACYCHOICE  efficiacychoice,    /**< type of solution to base efficiacy computation on */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   );

/** clears the separation storage without adding the cuts to the LP */
extern
SCIP_RETCODE SCIPsepastoreClearCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp                  /**< LP data */
   );

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP */
extern
SCIP_RETCODE SCIPsepastoreRemoveInefficaciousCuts(
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_EFFICIACYCHOICE  efficiacychoice     /**< type of solution to base efficiacy computation on */
   );

/** indicates whether a cut is applicable
 *
 * A cut is applicable if it is modifiable, not a bound change, or a bound change that changes bounds by at least epsilon.
 */
extern
SCIP_Bool SCIPsepastoreIsCutApplicable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             cut                 /**< cut to check */
   );

/** get cuts in the separation storage */
extern
SCIP_ROW** SCIPsepastoreGetCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** get number of cuts in the separation storage */
extern
int SCIPsepastoreGetNCuts(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** get total number of cuts found so far */
extern
int SCIPsepastoreGetNCutsFound(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** get number of cuts found so far in current separation round */
extern
int SCIPsepastoreGetNCutsFoundRound(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

/** get total number of cuts applied to the LPs */
extern
int SCIPsepastoreGetNCutsApplied(
   SCIP_SEPASTORE*       sepastore           /**< separation storage */
   );

#ifdef __cplusplus
}
#endif

#endif
