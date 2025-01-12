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

/**@file   type_prop.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for propagators
 * @author Tobias Achterberg
 */

/** @defgroup DEFPLUGINS_PROP Default Propagators
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default propagators of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_PROP_H__
#define __SCIP_TYPE_PROP_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"
#include "scip/type_timing.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Prop SCIP_PROP;               /**< propagator */
typedef struct SCIP_PropData SCIP_PROPDATA;       /**< locally defined propagator data */


/** copy method for propagator plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 */
#define SCIP_DECL_PROPCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop)

/** destructor of propagator to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 */
#define SCIP_DECL_PROPFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop)

/** initialization method of propagator (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 */
#define SCIP_DECL_PROPINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop)

/** deinitialization method of propagator (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 */
#define SCIP_DECL_PROPEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop)

/** presolving initialization method of propagator (called when presolving is about to begin)
 *
 *  This method is called when the presolving process is about to begin, even if presolving is turned off.  The
 *  propagator may use this call to initialize its presolving data, before the presolving process begins.
 *
 *  Necessary modifications that have to be performed even if presolving is turned off should be done here or in the
 *  presolving deinitialization call (SCIP_DECL_PROPEXITPRE()).
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 */
#define SCIP_DECL_PROPINITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop)

/** presolving deinitialization method of propagator (called after presolving has been finished)
 *
 *  This method is called after the presolving has been finished, even if presolving is turned off.
 *  The propagator may use this call e.g. to clean up its presolving data.
 *
 *  Besides necessary modifications and clean up, no time consuming operations should be performed, especially if the
 *  problem has already been solved.  Use the method SCIPgetStatus(), which in this case returns SCIP_STATUS_OPTIMAL,
 *  SCIP_STATUS_INFEASIBLE, SCIP_STATUS_UNBOUNDED, or SCIP_STATUS_INFORUNBD.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 */
#define SCIP_DECL_PROPEXITPRE(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop)

/** solving process initialization method of propagator (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The propagator may use this call to initialize its branch and bound specific data.
 *
 *  Besides necessary modifications and clean up, no time consuming operations should be performed, especially if the
 *  problem has already been solved.  Use the method SCIPgetStatus(), which in this case returns SCIP_STATUS_OPTIMAL,
 *  SCIP_STATUS_INFEASIBLE, SCIP_STATUS_UNBOUNDED, or SCIP_STATUS_INFORUNBD.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 */
#define SCIP_DECL_PROPINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop)

/** solving process deinitialization method of propagator (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The propagator should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 *  - restart         : was this exit solve call triggered by a restart?
 */
#define SCIP_DECL_PROPEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop, SCIP_Bool restart)

/** presolving method of propagator
 *
 *  The presolver should go through the variables and constraints and tighten the domains or
 *  constraints. Each tightening should increase the given total numbers of changes.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 *  - nrounds         : number of presolving rounds already done
 *  - presoltiming    : current presolving timing
 *  - nnewfixedvars   : number of variables fixed since the last call to the presolving method
 *  - nnewaggrvars    : number of variables aggregated since the last call to the presolving method
 *  - nnewchgvartypes : number of variable type changes since the last call to the presolving method
 *  - nnewchgbds      : number of variable bounds tightened since the last call to the presolving method
 *  - nnewholes       : number of domain holes added since the last call to the presolving method
 *  - nnewdelconss    : number of deleted constraints since the last call to the presolving method
 *  - nnewaddconss    : number of added constraints since the last call to the presolving method
 *  - nnewupgdconss   : number of upgraded constraints since the last call to the presolving method
 *  - nnewchgcoefs    : number of changed coefficients since the last call to the presolving method
 *  - nnewchgsides    : number of changed left or right hand sides since the last call to the presolving method
 *
 *  @note the counters state the changes since the last call including the changes of this presolving method during its
 *        last call
 *
 *  @note if the propagator uses dual information for presolving it is nesassary to check via calling SCIPallowWeakDualReds
 *        or SCIPallowStrongDualReds if dual reductions are allowed.
 *
 *  input/output:
 *  - nfixedvars      : pointer to total number of variables fixed of all presolvers
 *  - naggrvars       : pointer to total number of variables aggregated of all presolvers
 *  - nchgvartypes    : pointer to total number of variable type changes of all presolvers
 *  - nchgbds         : pointer to total number of variable bounds tightened of all presolvers
 *  - naddholes       : pointer to total number of domain holes added of all presolvers
 *  - ndelconss       : pointer to total number of deleted constraints of all presolvers
 *  - naddconss       : pointer to total number of added constraints of all presolvers
 *  - nupgdconss      : pointer to total number of upgraded constraints of all presolvers
 *  - nchgcoefs       : pointer to total number of changed coefficients of all presolvers
 *  - nchgsides       : pointer to total number of changed left/right hand sides of all presolvers
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_SUCCESS    : the presolving method found a reduction
 *  - SCIP_DIDNOTFIND : the presolving method searched, but did not find a presolving change
 *  - SCIP_DIDNOTRUN  : the presolving method was skipped
 *  - SCIP_DELAYED    : the presolving method was skipped, but should be called again
 */
#define SCIP_DECL_PROPPRESOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop, int nrounds, SCIP_PRESOLTIMING presoltiming,  \
      int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
      int nnewdelconss, int nnewaddconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides, \
      int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes, \
      int* ndelconss, int* naddconss, int* nupgdconss, int* nchgcoefs, int* nchgsides, SCIP_RESULT* result)

/** execution method of propagator
 *
 *  Searches for domain propagations. The method is called in the node processing loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 *  - proptiming      : current point in the node solving loop
 *  - result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the current node is infeasible for the current domains
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched, but did not find a domain reduction
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
 *  - SCIP_DELAYNODE  : the current node should be postponed (return value only valid for BEFORELP propagation)
 */
#define SCIP_DECL_PROPEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop,  SCIP_PROPTIMING proptiming, SCIP_RESULT* result)


/** propagation conflict resolving method of propagator
 *
 *  This method is called during conflict analysis. If the propagator wants to support conflict analysis,
 *  it should call SCIPinferVarLbProp() or SCIPinferVarUbProp() in domain propagation instead of SCIPchgVarLb() or
 *  SCIPchgVarUb() in order to deduce bound changes on variables.
 *  In the SCIPinferVarLbProp() and SCIPinferVarUbProp() calls, the propagator provides a pointer to itself
 *  and an integer value "inferinfo" that can be arbitrarily chosen.
 *  The propagation conflict resolving method can then be implemented, to provide a "reasons" for the bound
 *  changes, i.e. the bounds of variables at the time of the propagation, that forced the propagator to set the
 *  conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
 *  rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
 *  by calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(),
 *  SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), and/or SCIPaddConflictBinvar() in the propagation conflict
 *  resolving method.
 *
 *  See the description of the propagation conflict resolving method of constraint handlers for further details.
 *
 *  @note if the propagtor uses dual information it is nesassary to check via calling SCIPallowWeakDualReds and
 *        SCIPallowStrongDualReds if dual reductions and propgation with the current cutoff bound, resp., are allowed.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - prop            : the propagator itself
 *  - infervar        : the conflict variable whose bound change has to be resolved
 *  - inferinfo       : the user information passed to the corresponding SCIPinferVarLbProp() or SCIPinferVarUbProp() call
 *  - boundtype       : the type of the changed bound (lower or upper bound)
 *  - bdchgidx        : the index of the bound change, representing the point of time where the change took place
 *  - relaxedbd       : the relaxed bound which is sufficient to be explained
 *
 *  output:
 *  - result          : pointer to store the result of the propagation conflict resolving call
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the conflicting bound change has been successfully resolved by adding all reason bounds
 *  - SCIP_DIDNOTFIND : the conflicting bound change could not be resolved and has to be put into the conflict set
 *
 *  @note it is sufficient to explain/resolve the relaxed bound
 */
#define SCIP_DECL_PROPRESPROP(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROP* prop, SCIP_VAR* infervar, int inferinfo, \
      SCIP_BOUNDTYPE boundtype, SCIP_BDCHGIDX* bdchgidx, SCIP_Real relaxedbd, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
