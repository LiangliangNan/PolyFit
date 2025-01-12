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

/**@file   scip_conflict.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for conflict handler plugins and conflict analysis
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/conflict.h"
#include "scip/debug.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/scip_conflict.h"
#include "scip/scip_tree.h"
#include "scip/set.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_var.h"

/** creates a conflict handler and includes it in SCIP
 *
 *  @note method has all conflict handler callbacks as arguments and is thus changed every time a new
 *        callback is added
 *        in future releases; consider using SCIPincludeConflicthdlrBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_RETCODE SCIPincludeConflicthdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy)),  /**< copy method of conflict handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   )
{
   SCIP_CONFLICTHDLR* conflicthdlr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeConflicthdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether conflict handler is already present */
   if( SCIPfindConflicthdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("conflict handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconflicthdlrCreate(&conflicthdlr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         conflictcopy,
         conflictfree, conflictinit, conflictexit, conflictinitsol, conflictexitsol, conflictexec,
         conflicthdlrdata) );
   SCIP_CALL( SCIPsetIncludeConflicthdlr(scip->set, conflicthdlr) );

   return SCIP_OKAY;
}

/** creates a conflict handler and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions SCIPsetConflicthdlrCopy(), SCIPsetConflicthdlrFree(),
 *  SCIPsetConflicthdlrInit(), SCIPsetConflicthdlrExit(), SCIPsetConflicthdlrInitsol(),
 *  and SCIPsetConflicthdlrExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeConflicthdlr() instead
 */
SCIP_RETCODE SCIPincludeConflicthdlrBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR**   conflicthdlrptr,    /**< reference to a conflict handler pointer, or NULL */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   )
{
   SCIP_CONFLICTHDLR* conflicthdlr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeConflicthdlrBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether conflict handler is already present */
   if( SCIPfindConflicthdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("conflict handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconflicthdlrCreate(&conflicthdlr, scip->set, scip->messagehdlr, scip->mem->setmem, name, desc, priority,
         NULL, NULL, NULL, NULL, NULL, NULL, conflictexec, conflicthdlrdata) );
   SCIP_CALL( SCIPsetIncludeConflicthdlr(scip->set, conflicthdlr) );

   if( conflicthdlrptr != NULL )
      *conflicthdlrptr = conflicthdlr;

   return SCIP_OKAY;
}

/** set copy method of conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy))   /**< copy method of conflict handler */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConflicthdlrCopy", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconflicthdlrSetCopy(conflicthdlr, conflictcopy);

   return SCIP_OKAY;
}

/** set destructor of conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTFREE((*conflictfree))   /**< destructor of conflict handler */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConflicthdlrFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconflicthdlrSetFree(conflicthdlr, conflictfree);

   return SCIP_OKAY;
}

/** set initialization method of conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit))   /**< initialize conflict handler */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConflicthdlrInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconflicthdlrSetInit(conflicthdlr, conflictinit);

   return SCIP_OKAY;
}

/** set deinitialization method of conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit))   /**< deinitialize conflict handler */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConflicthdlrExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconflicthdlrSetExit(conflicthdlr, conflictexit);

   return SCIP_OKAY;
}

/** set solving process initialization method of conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol))/**< solving process initialization method of conflict handler */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConflicthdlrInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconflicthdlrSetInitsol(conflicthdlr, conflictinitsol);

   return SCIP_OKAY;
}

/** set solving process deinitialization method of conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol))/**< solving process deinitialization method of conflict handler */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetConflicthdlrExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPconflicthdlrSetExitsol(conflicthdlr, conflictexitsol);

   return SCIP_OKAY;
}

/** returns the conflict handler of the given name, or NULL if not existing */
SCIP_CONFLICTHDLR* SCIPfindConflicthdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of conflict handler */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindConflicthdlr(scip->set, name);
}

/** returns the array of currently available conflict handlers */
SCIP_CONFLICTHDLR** SCIPgetConflicthdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortConflicthdlrs(scip->set);

   return scip->set->conflicthdlrs;
}

/** returns the number of currently available conflict handlers */
int SCIPgetNConflicthdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nconflicthdlrs;
}

/** sets the priority of a conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   int                   priority            /**< new priority of the conflict handler */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPconflicthdlrSetPriority(conflicthdlr, scip->set, priority);

   return SCIP_OKAY;
}

/** return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *  conflict analysis since it will not be applied
 *
 *  @return return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *          conflict analysis since it will not be applied
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_Bool SCIPisConflictAnalysisApplicable(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisConflictAnalysisApplicable", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return (SCIPgetDepth(scip) > 0 && SCIPconflictApplicable(scip->set));
}

/** initializes the conflict analysis by clearing the conflict candidate queue; this method must be called before you
 *  enter the conflict variables by calling SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPinitConflictAnalysis(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFTYPE         conftype,           /**< type of conflict */
   SCIP_Bool             iscutoffinvolved    /**< is the current cutoff bound involved? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPinitConflictAnalysis", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictInit(scip->conflict, scip->set, scip->stat, scip->transprob, conftype, iscutoffinvolved) );

   return SCIP_OKAY;
}

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictLb() should be called for each lower bound
 *      that led to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictLb() should be called
 *      for each lower bound, whose current assignment led to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPaddConflictLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose lower bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflictLb", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, SCIP_BOUNDTYPE_LOWER, bdchgidx) );

   return SCIP_OKAY;
}

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage
 *  with the additional information of a relaxed lower bound; this relaxed lower bound is the one which would be enough
 *  to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedLb() should be called for each (relaxed) lower bound
 *      that led to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelexedLb() should be called
 *      for each (relaxed) lower bound, whose current assignment led to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPaddConflictRelaxedLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose lower bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Real             relaxedlb           /**< the relaxed lower bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflictRelaxedLb", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPconflictAddRelaxedBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, SCIP_BOUNDTYPE_LOWER, bdchgidx, relaxedlb) );

   return SCIP_OKAY;
}

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictUb() should be called for each upper bound that
 *      led to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictUb() should be called for
 *      each upper bound, whose current assignment led to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPaddConflictUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflictUb", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, SCIP_BOUNDTYPE_UPPER, bdchgidx) );

   return SCIP_OKAY;
}

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage
 *  with the additional information of a relaxed upper bound; this relaxed upper bound is the one which would be enough
 *  to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedUb() should be called for each (relaxed) upper
 *      bound that led to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelaxedUb() should be
 *      called for each (relaxed) upper bound, whose current assignment led to the deduction of the given conflict
 *      bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPaddConflictRelaxedUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Real             relaxedub           /**< the relaxed upper bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflictRelaxedUb", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPconflictAddRelaxedBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, SCIP_BOUNDTYPE_UPPER, bdchgidx, relaxedub) );

   return SCIP_OKAY;
}

/** adds lower or upper bound of variable at the time of the given bound change index to the conflict analysis' candidate
 *  storage; this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBd() should be called for each bound
 *      that led to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBd() should be called
 *      for each bound, whose current assignment led to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPaddConflictBd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the conflicting bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflictBd", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, boundtype, bdchgidx) );

   return SCIP_OKAY;
}

/** adds lower or upper bound of variable at the time of the given bound change index to the conflict analysis'
 *  candidate storage; with the additional information of a relaxed upper bound; this relaxed upper bound is the one
 *  which would be enough to explain a certain bound change;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictRelaxedBd() should be called for each (relaxed)
 *      bound that led to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictRelaxedBd() should be
 *      called for each (relaxed) bound, whose current assignment led to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPaddConflictRelaxedBd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the conflicting bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Real             relaxedbd           /**< the relaxed bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflictRelaxedBd", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPconflictAddRelaxedBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, boundtype, bdchgidx, relaxedbd) );

   return SCIP_OKAY;
}

/** adds changed bound of fixed binary variable to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBinvar() should be called for each fixed binary
 *      variable that led to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBinvar() should be called
 *      for each binary variable, whose current fixing led to the deduction of the given conflict bound.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPaddConflictBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< binary variable whose changed bound should be added to conflict queue */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddConflictBinvar", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(var->scip == scip);
   assert(SCIPvarIsBinary(var));

   if( SCIPvarGetLbLocal(var) > 0.5 )
   {
      SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, SCIP_BOUNDTYPE_LOWER, NULL) );
   }
   else if( SCIPvarGetUbLocal(var) < 0.5 )
   {
      SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->mem->probmem, scip->set, scip->stat, var, SCIP_BOUNDTYPE_UPPER, NULL) );
   }

   return SCIP_OKAY;
}

/** checks if the given variable is already part of the current conflict set or queued for resolving with the same or
 *  even stronger bound
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPisConflictVarUsed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the conflicting bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   SCIP_Bool*            used                /**< pointer to store if the variable is already used */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisConflictVarUsed", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPconflictIsVarUsed(scip->conflict, var, scip->set, boundtype, bdchgidx, used);
}

/** returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *  bound
 *
 *  @return returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *          bound
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_Real SCIPgetConflictVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetConflictVarLb", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPconflictGetVarLb(scip->conflict, var);
}

/** returns the conflict upper bound if the variable is present in the current conflict set; otherwise minus global
 *  upper bound
 *
 *  @return returns the conflict upper bound if the variable is present in the current conflict set; otherwise minus global
 *          upper bound
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_Real SCIPgetConflictVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetConflictVarUb", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPconflictGetVarUb(scip->conflict, var);
}

/** analyzes conflict bounds that were added after a call to SCIPinitConflictAnalysis() with calls to
 *  SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(),
 *  SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar(); on success, calls the conflict
 *  handlers to create a conflict constraint out of the resulting conflict set; the given valid depth must be a depth
 *  level, at which the conflict set defined by calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), and SCIPaddConflictBinvar() is
 *  valid for the whole subtree; if the conflict was found by a violated constraint, use SCIPanalyzeConflictCons()
 *  instead of SCIPanalyzeConflict() to make sure, that the correct valid depth is used
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPanalyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPanalyzeConflict", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictAnalyze(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
         scip->transprob, scip->tree, validdepth, success) );

   return SCIP_OKAY;
}

/** analyzes conflict bounds that were added with calls to SCIPaddConflictLb(), SCIPaddConflictUb(),
 *  SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or
 *  SCIPaddConflictBinvar(); on success, calls the conflict handlers to create a conflict constraint out of the
 *  resulting conflict set; the given constraint must be the constraint that detected the conflict, i.e. the constraint
 *  that is infeasible in the local bounds of the initial conflict set (defined by calls to SCIPaddConflictLb(),
 *  SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(),
 *  SCIPaddConflictRelaxedBd(), and SCIPaddConflictBinvar())
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note SCIP stage does not get changed
 */
SCIP_RETCODE SCIPanalyzeConflictCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that detected the conflict */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPanalyzeConflictCons", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPconsIsGlobal(cons) )
   {
      SCIP_CALL( SCIPconflictAnalyze(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->tree, 0, success) );
   }
   else if( SCIPconsIsActive(cons) )
   {
      SCIP_CALL( SCIPconflictAnalyze(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->tree, SCIPconsGetValidDepth(cons), success) );
   }

   return SCIP_OKAY;
}
