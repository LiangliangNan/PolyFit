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

/**@file   prop.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_H__
#define __SCIP_PROP_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prop.h"
#include "scip/pub_prop.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given propagator to a new scip */
extern
SCIP_RETCODE SCIPpropCopyInclude(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a propagator */
extern
SCIP_RETCODE SCIPpropCreate(
   SCIP_PROP**           prop,               /**< pointer to propagator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagator should be executed */
   int                   presolpriority,     /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the propagator's presolving method */
   SCIP_DECL_PROPCOPY    ((*propcopy)),      /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   SCIP_DECL_PROPINITPRE ((*propinitpre)),   /**< presolving initialization method of propagator */
   SCIP_DECL_PROPEXITPRE ((*propexitpre)),   /**< presolving deinitialization method of propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   );

/** calls destructor and frees memory of propagator */
extern
SCIP_RETCODE SCIPpropFree(
   SCIP_PROP**           prop,               /**< pointer to propagator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes propagator */
extern
SCIP_RETCODE SCIPpropInit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of propagator */
extern
SCIP_RETCODE SCIPpropExit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs propagator that the presolving process is being started */
extern
SCIP_RETCODE SCIPpropInitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs propagator that the presolving is finished */
extern
SCIP_RETCODE SCIPpropExitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs propagator that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPpropInitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs propagator that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPpropExitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   );

/** executes presolving method of propagator */
extern
SCIP_RETCODE SCIPpropPresol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRESOLTIMING     timing,             /**< current presolving timing */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls execution method of propagator */
extern
SCIP_RETCODE SCIPpropExec(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute propagator even if it is marked to be delayed */
   SCIP_Bool             instrongbranching,  /**< are we currently doing strong branching? */
   SCIP_PROPTIMING       proptiming,         /**< current point in the node solving process */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** resolves the given conflicting bound, that was deduced by the given propagator, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  SCIPaddConflictRelaxedLb(), SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), or SCIPaddConflictBinvar();
 *
 *  @note it is sufficient to explain the relaxed bound change
 */
extern
SCIP_RETCODE SCIPpropResolvePropagation(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int                   inferinfo,          /**< user inference information attached to the bound change */
   SCIP_BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   SCIP_Real             relaxedbd,          /**< the relaxed bound */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of propagator */
extern
void SCIPpropSetPriority(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the propagator */
   );

/** sets presolving priority of propagator */
extern
void SCIPpropSetPresolPriority(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   presolpriority      /**< new priority of the propagator */
   );

/** sets copy method of propagator */
extern
void SCIPpropSetCopy(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPCOPY    ((*propcopy))       /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of propagator */
extern
void SCIPpropSetFree(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPFREE    ((*propfree))       /**< destructor of propagator */
   );

/** sets initialization method of propagator */
extern
void SCIPpropSetInit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINIT    ((*propinit))       /**< initialize propagator */
   );

/** sets deinitialization method of propagator */
extern
void SCIPpropSetExit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXIT    ((*propexit))       /**< deinitialize propagator */
   );

/** sets solving process initialization method of propagator */
extern
void SCIPpropSetInitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITSOL((*propinitsol))     /**< solving process initialization method of propagator */
   );

/** sets solving process deinitialization method of propagator */
extern
void SCIPpropSetExitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol))    /**< solving process deinitialization method of propagator */
   );

/** sets preprocessing initialization method of propagator */
extern
void SCIPpropSetInitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPINITPRE((*propinitpre))     /**< preprocessing initialization method of propagator */
   );

/** sets preprocessing deinitialization method of propagator */
extern
void SCIPpropSetExitpre(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPEXITPRE((*propexitpre))     /**< preprocessing deinitialization method of propagator */
   );

/** sets presolving method of propagator */
extern
SCIP_RETCODE SCIPpropSetPresol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   int                   presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the propagator's presolving method */
   );

/** sets propagation conflict resolving callback of propagator */
extern
void SCIPpropSetResprop(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop))    /**< propagation conflict resolving callback */
   );

/** enables or disables all clocks of \p prop, depending on the value of the flag */
extern
void SCIPpropEnableOrDisableClocks(
   SCIP_PROP*            prop,               /**< the propagator for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the propagator be enabled? */
   );

#ifdef __cplusplus
}
#endif

#endif
