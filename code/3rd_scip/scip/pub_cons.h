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

/**@file   pub_cons.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for managing constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CONS_H__
#define __SCIP_PUB_CONS_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_cons.h"

#ifdef NDEBUG
#include "scip/struct_cons.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Constraint handler methods
 */

/**@addtogroup PublicConshdlrMethods
 *
 * @{
 */

/** compares two constraint handlers w. r. to their separation priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompSepa);

/** compares two constraint handlers w. r. to their enforcing priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo);

/** compares two constraint handlers w. r. to their feasibility check priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPconshdlrCompCheck);

/** gets name of constraint handler */
SCIP_EXPORT
const char* SCIPconshdlrGetName(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets description of constraint handler */
SCIP_EXPORT
const char* SCIPconshdlrGetDesc(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets user data of constraint handler */
SCIP_EXPORT
SCIP_CONSHDLRDATA* SCIPconshdlrGetData(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** sets user data of constraint handler; user has to free old data in advance! */
SCIP_EXPORT
void SCIPconshdlrSetData(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< new constraint handler user data */
   );

/** sets all separation related callbacks of the constraint handler */
SCIP_EXPORT
void SCIPconshdlrSetSepa(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   SCIP_Bool             delaysepa           /**< should separation method be delayed, if other separators found cuts? */
   );

/** sets both the propagation callback and the propagation frequency of the constraint handler */
SCIP_EXPORT
void SCIPconshdlrSetProp(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask          /**< positions in the node solving loop where propagators should be executed */
   );

/** sets the relaxation enforcement method of the constraint handler */
SCIP_EXPORT
void SCIPconshdlrSetEnforelax(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_DECL_CONSENFORELAX ((*consenforelax)) /**< constraint copying method */
   );

/** gets array with constraints of constraint handler; the first SCIPconshdlrGetNActiveConss() entries are the active
 *  constraints, the last SCIPconshdlrGetNConss() - SCIPconshdlrGetNActiveConss() constraints are deactivated
 *
 *  @note A constraint is active if it is global and was not removed or it was added locally (in that case the local
 *        flag is TRUE) and the current node belongs to the corresponding sub tree.
 */ 
SCIP_EXPORT
SCIP_CONS** SCIPconshdlrGetConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets array with enforced constraints of constraint handler; this is local information */
SCIP_EXPORT
SCIP_CONS** SCIPconshdlrGetEnfoConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets array with checked constraints of constraint handler; this is local information */
SCIP_EXPORT
SCIP_CONS** SCIPconshdlrGetCheckConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets array with delayed update constraints
 *
 * @attention Usually, there should be no need to access this array. Use this only if you are absolutely sure what you are doing.
 */
SCIP_EXPORT
SCIP_CONS** SCIPconshdlrGetUpdateConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of existing transformed constraints of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of enforced constraints of constraint handler; this is local information */
SCIP_EXPORT
int SCIPconshdlrGetNEnfoConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of checked constraints of constraint handler; this is local information */
SCIP_EXPORT
int SCIPconshdlrGetNCheckConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of active constraints of constraint handler
 *
 *  @note A constraint is active if it is global and was not removed or it was added locally (in that case the local
 *        flag is TRUE) and the current node belongs to the corresponding sub tree.
 */ 
SCIP_EXPORT
int SCIPconshdlrGetNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of enabled constraints of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNEnabledConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints that have delayed updates */
SCIP_EXPORT
int SCIPconshdlrGetNUpdateConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for setting up this constraint handler for new stages */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetSetupTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for presolving in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetPresolTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for separation in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetSepaTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for LP enforcement in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetEnfoLPTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for pseudo enforcement in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetEnfoPSTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for relaxation enforcement in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetEnfoRelaxTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for propagation in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetPropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for propagation in this constraint handler during strong branching */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetStrongBranchPropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for feasibility checking in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetCheckTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for resolving propagation in this constraint handler */
SCIP_EXPORT
SCIP_Real SCIPconshdlrGetRespropTime(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's separation method */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNSepaCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's LP enforcing method */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNEnfoLPCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's pseudo enforcing method */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNEnfoPSCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's relaxation enforcing method */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNEnfoRelaxCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's propagation method */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNPropCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's checking method */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNCheckCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's resolve propagation method */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNRespropCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of times, this constraint handler detected a cutoff */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNCutoffs(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of cuts found by this constraint handler */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNCutsFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of cuts found by this constraint handler applied to lp */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNCutsApplied(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of additional constraints added by this constraint handler */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNConssFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of domain reductions found by this constraint handler */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNDomredsFound(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of children created by this constraint handler */
SCIP_EXPORT
SCIP_Longint SCIPconshdlrGetNChildren(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets maximum number of active constraints of constraint handler existing at the same time */
SCIP_EXPORT
int SCIPconshdlrGetMaxNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets initial number of active constraints of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetStartNActiveConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variables fixed in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNFixedVars(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variables aggregated in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNAggrVars(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variable types changed in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNChgVarTypes(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of bounds changed in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNChgBds(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of holes added to domains of variables in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNAddHoles(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints deleted in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNDelConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints added in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNAddConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints upgraded in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNUpgdConss(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of coefficients changed in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNChgCoefs(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraint sides changed in presolving method of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetNChgSides(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of times the presolving method of the constraint handler was called and tried to find reductions */
SCIP_EXPORT
int SCIPconshdlrGetNPresolCalls(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets separation priority of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetSepaPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets enforcing priority of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetEnfoPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets checking priority of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetCheckPriority(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets separation frequency of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetSepaFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets propagation frequency of constraint handler */
SCIP_EXPORT
int SCIPconshdlrGetPropFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets frequency of constraint handler for eager evaluations in separation, propagation and enforcement */
SCIP_EXPORT
int SCIPconshdlrGetEagerFreq(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** needs constraint handler a constraint to be called? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrNeedsCons(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** does the constraint handler perform presolving? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrDoesPresolve(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** should separation method be delayed, if other separators found cuts? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrIsSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** should propagation method be delayed, if other propagators found reductions? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrIsPropagationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** was LP separation method delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrWasLPSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** was primal solution separation method delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrWasSolSeparationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** was propagation method delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrWasPropagationDelayed(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** is constraint handler initialized? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrIsInitialized(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** does the constraint handler have a copy function? */
SCIP_EXPORT
SCIP_Bool SCIPconshdlrIsClonable(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** returns the timing mask of the propagation method of the constraint handler */
SCIP_EXPORT
SCIP_PROPTIMING SCIPconshdlrGetPropTiming(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/*
 * Methods for constraint change sets
 */
/** gets added constraints data for a constraint set change */
SCIP_EXPORT
void SCIPconssetchgGetAddedConsData(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to get data from */
   SCIP_CONS***          conss,              /**< reference to constraints array added in the conssetchg, or NULL */
   int*                  nconss              /**< reference to store the size of the constraints array, or NULL */
   );

/** sets the timing mask of the propagation method of the constraint handler */
SCIP_EXPORT
void SCIPconshdlrSetPropTiming(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_PROPTIMING       proptiming          /**< timing mask to be set */
   );


/** returns the timing mask of the presolving method of the constraint handler */
SCIP_EXPORT
SCIP_PRESOLTIMING SCIPconshdlrGetPresolTiming(
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** sets the timing mask of the presolving method of the constraint handler */
SCIP_EXPORT
void SCIPconshdlrSetPresolTiming(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_PRESOLTIMING     presoltiming        /** timing mask to be set */
   );

/** @} */

/*
 * Constraint methods
 */

/**@addtogroup PublicConstraintMethods
 *
 * @{
 */


/** returns the name of the constraint 
 *
 *  @note to change the name of a constraint, use SCIPchgConsName() from scip.h
 */
SCIP_EXPORT
const char* SCIPconsGetName(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the position of constraint in the corresponding handler's conss array */
SCIP_EXPORT
int SCIPconsGetPos(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the constraint handler of the constraint */
SCIP_EXPORT
SCIP_CONSHDLR* SCIPconsGetHdlr(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the constraint data field of the constraint */
SCIP_EXPORT
SCIP_CONSDATA* SCIPconsGetData(
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets number of times, the constraint is currently captured */
SCIP_EXPORT
int SCIPconsGetNUses(
   SCIP_CONS*            cons                /**< constraint */
   );

/** for an active constraint, returns the depth in the tree at which the constraint was activated */
SCIP_EXPORT
int SCIPconsGetActiveDepth(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns the depth in the tree at which the constraint is valid; returns INT_MAX, if the constraint is local
 *  and currently not active
 */
SCIP_EXPORT
int SCIPconsGetValidDepth(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is active in the current node */
SCIP_EXPORT
SCIP_Bool SCIPconsIsActive(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint has to be deactivated in update phase */
SCIP_EXPORT
SCIP_Bool SCIPconsIsUpdatedeactivate(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is enabled in the current node */
SCIP_EXPORT
SCIP_Bool SCIPconsIsEnabled(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's separation is enabled in the current node */
SCIP_EXPORT
SCIP_Bool SCIPconsIsSeparationEnabled(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's propagation is enabled in the current node */
SCIP_EXPORT
SCIP_Bool SCIPconsIsPropagationEnabled(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is deleted or marked to be deleted */
SCIP_EXPORT
SCIP_Bool SCIPconsIsDeleted(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is marked obsolete */
SCIP_EXPORT
SCIP_Bool SCIPconsIsObsolete(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is marked as a conflict */
SCIP_EXPORT
SCIP_Bool SCIPconsIsConflict(
   SCIP_CONS*            cons                /**< constraint */
   );

/** gets age of constraint */
SCIP_EXPORT
SCIP_Real SCIPconsGetAge(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff the LP relaxation of constraint should be in the initial LP */
SCIP_EXPORT
SCIP_Bool SCIPconsIsInitial(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be separated during LP processing */
SCIP_EXPORT
SCIP_Bool SCIPconsIsSeparated(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be enforced during node processing */
SCIP_EXPORT
SCIP_Bool SCIPconsIsEnforced(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be checked for feasibility */
SCIP_EXPORT
SCIP_Bool SCIPconsIsChecked(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns whether the constraint is marked for propagation */
SCIP_EXPORT
SCIP_Bool SCIPconsIsMarkedPropagate(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be propagated during node processing */
SCIP_EXPORT
SCIP_Bool SCIPconsIsPropagated(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is globally valid */
SCIP_EXPORT
SCIP_Bool SCIPconsIsGlobal(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is only locally valid or not added to any (sub)problem */
SCIP_EXPORT
SCIP_Bool SCIPconsIsLocal(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is modifiable (subject to column generation) */
SCIP_EXPORT
SCIP_Bool SCIPconsIsModifiable(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is subject to aging */
SCIP_EXPORT
SCIP_Bool SCIPconsIsDynamic(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
SCIP_EXPORT
SCIP_Bool SCIPconsIsRemovable(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint's relaxation should be removed from the LP due to aging or cleanup */
SCIP_EXPORT
SCIP_Bool SCIPconsIsStickingAtNode(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint belongs to the global problem */
SCIP_EXPORT
SCIP_Bool SCIPconsIsInProb(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is belonging to original space */
SCIP_EXPORT
SCIP_Bool SCIPconsIsOriginal(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is belonging to transformed space */
SCIP_EXPORT
SCIP_Bool SCIPconsIsTransformed(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff roundings for variables in constraint are locked */
SCIP_EXPORT
SCIP_Bool SCIPconsIsLockedPos(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff roundings for variables in constraint's negation are locked */
SCIP_EXPORT
SCIP_Bool SCIPconsIsLockedNeg(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff roundings for variables in constraint or in constraint's negation are locked */
SCIP_EXPORT
SCIP_Bool SCIPconsIsLocked(
   SCIP_CONS*            cons                /**< constraint */
   );

/** get number of times the roundings for variables in constraint are locked */
SCIP_EXPORT
int SCIPconsGetNLocksPos(
   SCIP_CONS*            cons                /**< constraint */
   );

/** get number of times the roundings for variables in constraint's negation are locked */
SCIP_EXPORT
int SCIPconsGetNLocksNeg(
   SCIP_CONS*            cons                /**< constraint */
   );

/** returns TRUE iff roundings of the given locktype for variables in constraint are locked */
SCIP_EXPORT
SCIP_Bool SCIPconsIsLockedTypePos(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_LOCKTYPE         locktype            /**< variable lock type */
   );

/** returns TRUE iff roundings of the given locktype for variables in constraint are locked */
SCIP_EXPORT
SCIP_Bool SCIPconsIsLockedTypeNeg(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_LOCKTYPE         locktype            /**< variable lock type */
   );

/** returns TRUE iff roundings of the given locktype for variables in constraint or in constraint's negation are locked */
SCIP_EXPORT
SCIP_Bool SCIPconsIsLockedType(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_LOCKTYPE         locktype            /**< variable lock type */
   );

/** get number of times the roundings of given locktype for variables in constraint are locked */
SCIP_EXPORT
int SCIPconsGetNLocksTypePos(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_LOCKTYPE         locktype            /**< variable lock type */
   );

/** get number of times the roundings of given locktype for variables in constraint's negation are locked */
SCIP_EXPORT
int SCIPconsGetNLocksTypeNeg(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_LOCKTYPE         locktype            /**< variable lock type */
   );

/** returns if the constraint was already added to a SCIP instance */
SCIP_EXPORT
SCIP_Bool SCIPconsIsAdded(
   SCIP_CONS*            cons                /**< constraint */
   );

/** adds locks to (dis-)allow upgrading of constraint */
SCIP_EXPORT
void SCIPconsAddUpgradeLocks(
   SCIP_CONS*            cons,               /**< constraint to add locks */
   int                   nlocks              /**< number of locks to add */
   );

/** gets number of locks against upgrading the constraint, 0 means this constraint can be upgraded */
SCIP_EXPORT
int SCIPconsGetNUpgradeLocks(
   SCIP_CONS*            cons                /**< constraint */
   );

#ifdef NDEBUG

/* In optimized mode, the function calls are overwritten by defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPconsGetName(cons)           (cons)->name
#define SCIPconsGetPos(cons)            (cons)->consspos
#define SCIPconsGetHdlr(cons)           (cons)->conshdlr
#define SCIPconsGetData(cons)           (cons)->consdata
#define SCIPconsGetNUses(cons)          (cons)->nuses
#define SCIPconsGetActiveDepth(cons)    (cons)->activedepth
#define SCIPconsGetValidDepth(cons)     (!(cons)->local ? 0     \
      : !SCIPconsIsActive(cons) ? INT_MAX                       \
      : (cons)->validdepth == -1 ? SCIPconsGetActiveDepth(cons) \
      : (cons)->validdepth)
#define SCIPconsIsActive(cons)          ((cons)->updateactivate || ((cons)->active && !(cons)->updatedeactivate))
#define SCIPconsIsEnabled(cons)         ((cons)->updateenable || ((cons)->enabled && !(cons)->updatedisable))
#define SCIPconsIsSeparationEnabled(cons)                               \
   (SCIPconsIsEnabled(cons) && ((cons)->updatesepaenable || ((cons)->sepaenabled && !(cons)->updatesepadisable)))
#define SCIPconsIsPropagationEnabled(cons)                              \
   (SCIPconsIsEnabled(cons) && ((cons)->updatepropenable || ((cons)->propenabled && !(cons)->updatepropdisable)))
#define SCIPconsIsDeleted(cons)         ((cons)->deleted)
#define SCIPconsIsObsolete(cons)        ((cons)->updateobsolete || (cons)->obsolete)
#define SCIPconsIsConflict(cons)        ((cons)->conflict)
#define SCIPconsGetAge(cons)            (cons)->age
#define SCIPconsIsInitial(cons)         (cons)->initial
#define SCIPconsIsSeparated(cons)       (cons)->separate
#define SCIPconsIsEnforced(cons)        (cons)->enforce
#define SCIPconsIsChecked(cons)         (cons)->check
#define SCIPconsIsMarkedPropagate(cons) ((cons)->updatemarkpropagate || ((cons)->markpropagate && !(cons)->updateunmarkpropagate))
#define SCIPconsIsPropagated(cons)      (cons)->propagate
#define SCIPconsIsGlobal(cons)          !(cons)->local
#define SCIPconsIsLocal(cons)           (cons)->local
#define SCIPconsIsModifiable(cons)      (cons)->modifiable
#define SCIPconsIsDynamic(cons)         (cons)->dynamic
#define SCIPconsIsRemovable(cons)       (cons)->removable
#define SCIPconsIsStickingAtNode(cons)  (cons)->stickingatnode
#define SCIPconsIsInProb(cons)          ((cons)->addconssetchg == NULL && (cons)->addarraypos >= 0)
#define SCIPconsIsOriginal(cons)        (cons)->original
#define SCIPconsIsTransformed(cons)     !(cons)->original
#define SCIPconsIsLockedPos(cons)       ((cons)->nlockspos[SCIP_LOCKTYPE_MODEL] > 0)
#define SCIPconsIsLockedNeg(cons)       ((cons)->nlocksneg[SCIP_LOCKTYPE_MODEL] > 0)
#define SCIPconsIsLocked(cons)          ((cons)->nlockspos[SCIP_LOCKTYPE_MODEL] > 0 || (cons)->nlocksneg[SCIP_LOCKTYPE_MODEL] > 0)
#define SCIPconsGetNLocksPos(cons)      ((cons)->nlockspos[SCIP_LOCKTYPE_MODEL])
#define SCIPconsGetNLocksNeg(cons)      ((cons)->nlocksneg[SCIP_LOCKTYPE_MODEL])
#define SCIPconsIsLockedTypePos(cons, locktype)  ((cons)->nlockspos[locktype] > 0)
#define SCIPconsIsLockedTypeNeg(cons, locktype)  ((cons)->nlocksneg[locktype] > 0)
#define SCIPconsIsLockedType(cons, locktype)     ((cons)->nlockspos[locktype] > 0 || (cons)->nlocksneg[locktype] > 0)
#define SCIPconsGetNLocksTypePos(cons, locktype) ((cons)->nlockspos[locktype])
#define SCIPconsGetNLocksTypeNeg(cons, locktype) ((cons)->nlocksneg[locktype])
#define SCIPconsIsAdded(cons)           ((cons)->addarraypos >= 0)
#define SCIPconsGetNUpgradeLocks(cons)  ((cons)->nupgradelocks)

#endif

/** @} */

/**@addtogroup PublicProblemMethods
 *
 * public methods to query linear constraint classification statistics
 *
 * @{
 */

/** create linear constraint statistics */
SCIP_EXPORT
SCIP_RETCODE SCIPlinConsStatsCreate(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_LINCONSSTATS**   linconsstats        /**< pointer to linear constraint classification statistics */
   );

/** free linear constraint statistics */
SCIP_EXPORT
void SCIPlinConsStatsFree(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_LINCONSSTATS**   linconsstats        /**< pointer to linear constraint classification statistics */
   );

/** resets linear constraint statistics */
SCIP_EXPORT
void SCIPlinConsStatsReset(
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint classification statistics */
   );

/** returns the number of occurrences of a specific type of linear constraint */
SCIP_EXPORT
int SCIPlinConsStatsGetTypeCount(
   SCIP_LINCONSSTATS*    linconsstats,       /**< linear constraint classification statistics */
   SCIP_LINCONSTYPE      linconstype         /**< linear constraint type */
   );

/** returns the total number of classified constraints */
SCIP_EXPORT
int SCIPlinConsStatsGetSum(
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint classification statistics */
   );

/** increases the number of occurrences of a specific type of linear constraint */
SCIP_EXPORT
void SCIPlinConsStatsIncTypeCount(
   SCIP_LINCONSSTATS*    linconsstats,       /**< linear constraint classification statistics */
   SCIP_LINCONSTYPE      linconstype,        /**< linear constraint type */
   int                   increment           /**< positive increment */
   );

/** print linear constraint classification statistics */
SCIP_EXPORT
void SCIPprintLinConsStats(
   SCIP*                 scip,               /**< scip data structure */
   FILE*                 file,               /**< file handle or NULL to print to standard out */
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint classification statistics */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
