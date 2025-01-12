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

/**@file   heur.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_H__
#define __SCIP_HEUR_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_primal.h"
#include "scip/type_heur.h"
#include "scip/pub_heur.h"
#include "scip/stat.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create a set of diving heuristic settings */
SCIP_RETCODE SCIPdivesetCreate(
   SCIP_DIVESET**        divesetptr,         /**< pointer to the freshly created diveset */
   SCIP_HEUR*            heur,               /**< the heuristic to which this dive setting belongs */
   const char*           name,               /**< name for the diveset, or NULL if the name of the heuristic should be used */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   SCIP_Real             minreldepth,        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth,        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot,      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   SCIP_Real             maxdiveubquot,      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot,     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol, /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol,/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             lpresolvedomchgquot,/**< percentage of immediate domain changes during probing to trigger LP resolve */
   int                   lpsolvefreq,        /**< LP solve frequency for (0: only if enough domain reductions are found by propagation)*/
   int                   maxlpiterofs,       /**< additional number of allowed LP iterations */
   unsigned int          initialseed,        /**< initial seed for random number generation */
   SCIP_Bool             backtrack,          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Bool             onlylpbranchcands,  /**< should only LP branching candidates be considered instead of the slower but
                                              *   more general constraint handler diving variable selection? */
   SCIP_Bool             ispublic,           /**< is this dive set publicly available (ie., can be used by other primal heuristics?) */
   SCIP_DIVETYPE         divetypemask,       /**< bit mask that represents the supported dive types by this dive set */
   SCIP_DECL_DIVESETGETSCORE((*divesetgetscore)), /**< method for candidate score and rounding direction */
   SCIP_DECL_DIVESETAVAILABLE((*divesetavailable)) /**< callback to check availability of dive set at the current stage, or NULL if always available */
   );

/** resets diving settings counters */
SCIP_RETCODE SCIPdivesetReset(
   SCIP_DIVESET*         diveset,            /**< diveset to be reset */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** update diveset statistics and global diveset statistics */
void SCIPdivesetUpdateStats(
   SCIP_DIVESET*         diveset,            /**< diveset to be reset */
   SCIP_STAT*            stat,               /**< global SCIP statistics */
   int                   depth,              /**< the depth reached this time */
   int                   nprobingnodes,      /**< the number of probing nodes explored this time */
   int                   nbacktracks,        /**< the number of backtracks during probing this time */
   SCIP_Longint          nsolsfound,         /**< number of new solutions found this time */
   SCIP_Longint          nbestsolsfound,     /**< number of new best solutions found this time */
   SCIP_Longint          nconflictsfound,    /**< number of new conflicts found this time */
   SCIP_Bool             leavesol,           /**< has the diving heuristic reached a feasible leaf */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** get the candidate score and preferred rounding direction for a candidate variable */
SCIP_RETCODE SCIPdivesetGetScore(
   SCIP_DIVESET*         diveset,            /**< general diving settings */
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_DIVETYPE         divetype,           /**< the type of diving that should be applied */
   SCIP_VAR*             divecand,           /**< the candidate for which the branching direction is requested */
   SCIP_Real             divecandsol,        /**< LP solution value of the candidate */
   SCIP_Real             divecandfrac,       /**< fractionality of the candidate */
   SCIP_Real*            candscore,          /**< pointer to store the candidate score */
   SCIP_Bool*            roundup             /**< pointer to store whether preferred direction for diving is upwards */
   );

/** check specific preconditions for diving, e.g., if an incumbent solution is available */
SCIP_RETCODE SCIPdivesetIsAvailable(
   SCIP_DIVESET*         diveset,            /**< diving heuristic settings */
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_Bool*            available           /**< pointer to store if the diving can run at the current solving stage */
   );

/** update diveset LP statistics, should be called after every LP solved by this diving heuristic */
void SCIPdivesetUpdateLPStats(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_STAT*            stat,               /**< global SCIP statistics */
   SCIP_Longint          niterstoadd,        /**< additional number of LP iterations to be added */
   SCIP_DIVECONTEXT      divecontext         /**< context for diving statistics */
   );

/** copies the given primal heuristic to a new scip */
SCIP_RETCODE SCIPheurCopyInclude(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a primal heuristic */
SCIP_RETCODE SCIPheurCreate(
   SCIP_HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   SCIP_HEURTIMING       timingmask,         /**< positions in the node solving loop where heuristic should be executed */
   SCIP_Bool             usessubscip,        /**< does the heuristic use a secondary SCIP instance? */
   SCIP_DECL_HEURCOPY    ((*heurcopy)),      /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol)),   /**< solving process initialization method of primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol)),   /**< solving process deinitialization method of primal heuristic */
   SCIP_DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   );

/** calls destructor and frees memory of primal heuristic */
SCIP_RETCODE SCIPheurFree(
   SCIP_HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** initializes primal heuristic */
SCIP_RETCODE SCIPheurInit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of primal heuristic */
SCIP_RETCODE SCIPheurExit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs primal heuristic that the branch and bound process is being started */
SCIP_RETCODE SCIPheurInitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs primal heuristic that the branch and bound process data is being freed */
SCIP_RETCODE SCIPheurExitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** should the heuristic be executed at the given depth, frequency, timing, ... */
SCIP_Bool SCIPheurShouldBeExecuted(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   depth,              /**< depth of current node */
   int                   lpstateforkdepth,   /**< depth of the last node with solved LP */
   SCIP_HEURTIMING       heurtiming,         /**< current point in the node solving process */
   SCIP_Bool*            delayed             /**< pointer to store whether the heuristic should be delayed */
   );

/** calls execution method of primal heuristic */
SCIP_RETCODE SCIPheurExec(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   int                   depth,              /**< depth of current node */
   int                   lpstateforkdepth,   /**< depth of the last node with solved LP */
   SCIP_HEURTIMING       heurtiming,         /**< current point in the node solving process */
   SCIP_Bool             nodeinfeasible,     /**< was the current node already detected to be infeasible? */
   int*                  ndelayedheurs,      /**< pointer to count the number of delayed heuristics */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of primal heuristic */
void SCIPheurSetPriority(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the primal heuristic */
   );

/** sets copy callback of primal heuristic */
void SCIPheurSetCopy(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURCOPY    ((*heurcopy))       /**< copy callback of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor callback of primal heuristic */
void SCIPheurSetFree(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURFREE    ((*heurfree))       /**< destructor of primal heuristic */
   );

/** sets initialization callback of primal heuristic */
void SCIPheurSetInit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit))       /**< initialize primal heuristic */
   );

/** sets deinitialization callback of primal heuristic */
void SCIPheurSetExit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit))       /**< deinitialize primal heuristic */
   );

/** sets solving process initialization callback of primal heuristic */
void SCIPheurSetInitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol))    /**< solving process initialization callback of primal heuristic */
   );

/** sets solving process deinitialization callback of primal heuristic */
void SCIPheurSetExitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol))    /**< solving process deinitialization callback of primal heuristic */
   );

/** enables or disables all clocks of \p heur, depending on the value of the flag */
void SCIPheurEnableOrDisableClocks(
   SCIP_HEUR*            heur,               /**< the heuristic for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the heuristic be enabled? */
   );

#ifdef __cplusplus
}
#endif

#endif
