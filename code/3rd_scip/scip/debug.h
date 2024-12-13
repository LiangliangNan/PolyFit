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

/**@file   debug.h
 * @ingroup INTERNALAPI
 * @brief  methods for debugging
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DEBUG_H__
#define __SCIP_DEBUG_H__

/** uncomment this define to activate debugging the LP interface  */
/* #define SCIP_DEBUG_LP_INTERFACE */

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef WITH_DEBUG_SOLUTION
#include "blockmemshell/memory.h"
#include "scip/type_cons.h"
#include "scip/type_lp.h"
#include "scip/type_misc.h"
#include "scip/type_set.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** solution data for debugging purposes */
typedef struct SCIP_DebugSolData SCIP_DEBUGSOLDATA;

#ifdef WITH_DEBUG_SOLUTION

/** creates debug solution data */
SCIP_RETCODE SCIPdebugSolDataCreate(
   SCIP_DEBUGSOLDATA**   debugsoldata        /**< pointer to debug solution data */
   );

/** frees the debug solution */
SCIP_RETCODE SCIPdebugFreeSol(
   SCIP_SET*             set
   );

/** resets the data structure after restart */
SCIP_RETCODE SCIPdebugReset(
   SCIP_SET*             set
   );

/** frees debugging data for the particular instance */
SCIP_RETCODE SCIPdebugFreeDebugData(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees all debugging data */
SCIP_RETCODE SCIPdebugFree(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** checks for validity of the debugging solution in given constraints */
SCIP_RETCODE SCIPdebugCheckConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to check for validity */
   int                   nconss              /**< number of given constraints */
   );

/** checks whether given row is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< row to check for validity */
   );

/** checks whether given global lower bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lb                  /**< lower bound */
   );

/** checks whether given global upper bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             ub                  /**< upper bound */
   );

/** checks whether given local bound implication is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckInference(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< local node where this bound change was applied */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** informs solution debugger, that the given node will be freed */
SCIP_RETCODE SCIPdebugRemoveNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node that will be freed */
   );

/** checks whether global lower bound does not exceed debuging solution value */
SCIP_RETCODE SCIPdebugCheckGlobalLowerbound(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** checks whether local lower bound does not exceed debuging solution value */
SCIP_RETCODE SCIPdebugCheckLocalLowerbound(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node that will be freed */
   );

/** checks whether given variable bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckVbound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable x in x <= b*z + d  or  x >= b*z + d */
   SCIP_BOUNDTYPE        vbtype,             /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             vbvar,              /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbcoef,             /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbconstant          /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   );

/** checks whether given implication is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckImplic(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound           /**< bound b    in implication y <= b or y >= b */
   );

/** checks whether given (multi)-aggregation is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckAggregation(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VAR**            aggrvars,           /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   int                   naggrvars           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   );

/** check whether given clique is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckClique(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< binary variables in the clique: at most one can be set to the given value */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars               /**< number of variables in the clique */
   );

/** checks whether given conflict is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckConflict(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node where the conflict clause is added */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change informations of the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos         /**< number of bound changes in the conflict set */
   );

/** checks whether given conflict graph frontier is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckConflictFrontier(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node where the conflict clause is added */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change info which got resolved, or NULL */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change informations of the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos,        /**< number of bound changes in the conflict set */
   SCIP_PQUEUE*          bdchgqueue,         /**< unprocessed conflict bound changes */
   SCIP_PQUEUE*          forcedbdchgqueue    /**< unprocessed conflict bound changes that must be resolved */
   );

/** creates the debugging propagator and includes it in SCIP */
SCIP_RETCODE SCIPdebugIncludeProp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds a solution value for a new variable in the transformed problem that has no original counterpart
 * a value can only be set if no value has been set for this variable before
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdebugAddSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which to add a value */
   SCIP_Real             val                 /**< solution value for variable */
   );

/** gets pointer to the debug solution */
SCIP_EXPORT
SCIP_RETCODE SCIPdebugGetSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol                 /**< buffer to store pointer to the debug solution */
   );

/** gets value for a variable in the debug solution
 *
 * if no value is stored for the variable, gives 0.0
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdebugGetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which to get the value */
   SCIP_Real*            val                 /**< buffer to store solution value */
   );

/** check whether the debugging solution is valid in the current node */
SCIP_EXPORT
SCIP_RETCODE SCIPdebugSolIsValidInSubtree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            isvalidinsubtree    /**< pointer to store whether the solution is valid in the current
                                              *   subtree
                                              */
   );

/** checks whether SCIP data structure is the main SCIP (the one for which debugging is enabled) */
SCIP_EXPORT
SCIP_Bool SCIPdebugIsMainscip(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** enabling solution debugging mechanism */
SCIP_EXPORT
void SCIPdebugSolEnable(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** disabling solution debugging mechanism */
SCIP_EXPORT
void SCIPdebugSolDisable(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** check if solution debugging mechanism is enabled */
SCIP_EXPORT
SCIP_Bool SCIPdebugSolIsEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** check if SCIP is compiled with WITH_DEBUG_SOLUTION */
SCIP_EXPORT
SCIP_Bool SCIPwithDebugSol(void);

#else

#define SCIPdebugSolDataCreate(debugsoldata) SCIP_OKAY
#define SCIPdebugFreeSol(set) SCIP_OKAY
#define SCIPdebugReset(set) SCIP_OKAY
#define SCIPdebugFreeDebugData(set) SCIP_OKAY
#define SCIPdebugFree(set) SCIP_OKAY
#define SCIPdebugCheckConss(scip,conss,nconss) SCIP_OKAY
#define SCIPdebugCheckRow(set,row) SCIP_OKAY
#define SCIPdebugCheckLbGlobal(scip,var,lb) SCIP_OKAY
#define SCIPdebugCheckUbGlobal(scip,var,ub) SCIP_OKAY
#define SCIPdebugCheckInference(blkmem,set,node,var,newbound,boundtype) SCIP_OKAY
#define SCIPdebugRemoveNode(blkmem,set,node) SCIP_OKAY
#define SCIPdebugCheckGlobalLowerbound(blkmem,set) SCIP_OKAY
#define SCIPdebugCheckLocalLowerbound(blkmem,set,node) SCIP_OKAY
#define SCIPdebugCheckVbound(set,var,vbtype,vbvar,vbcoef,vbconstant) SCIP_OKAY
#define SCIPdebugCheckImplic(set,var,varfixing,implvar,impltype,implbound) SCIP_OKAY
#define SCIPdebugCheckAggregation(set,var,aggrvars,scalars,constant,naggrvars) SCIP_OKAY
#define SCIPdebugCheckClique(set,vars,values,nvars) SCIP_OKAY
#define SCIPdebugCheckConflict(blkmem,set,node,bdchginfos,relaxedbds,nliterals) SCIP_OKAY
#define SCIPdebugCheckConflictFrontier(blkmem,set,node,bdchginfo,bdchginfos,relaxedbds,nliterals,bdchgqueue,forcedbdchgqueue) SCIP_OKAY
#define SCIPdebugIncludeProp(scip) SCIP_OKAY
#define SCIPdebugAddSolVal(scip,var,val) SCIP_OKAY
#define SCIPdebugGetSolVal(scip,var,val) SCIP_OKAY
#define SCIPdebugSolIsValidInSubtree(scip,isvalidinsubtree) SCIP_OKAY
#define SCIPdebugSolEnable(scip) /**/
#define SCIPdebugSolDisable(scip) /**/
#define SCIPdebugSolIsEnabled(scip) FALSE
#define SCIPwithDebugSol(void) FALSE

#endif


/* 
 * debug method for LP interface, to check if the LP interface works correct 
 */
#ifdef SCIP_DEBUG_LP_INTERFACE 

/* check if the coef is the r-th line of the inverse matrix B^-1; this is
 * the case if (coef * B) is the r-th unit vector */
SCIP_RETCODE SCIPdebugCheckBInvRow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   );

#else

#define SCIPdebugCheckBInvRow(scip,r,coef) SCIP_OKAY

#endif

/** checks, if SCIP is in one of the feasible stages */
#ifndef NDEBUG

SCIP_RETCODE SCIPcheckStage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           method,             /**< method that was called */
   SCIP_Bool             init,               /**< may method be called in the INIT stage? */
   SCIP_Bool             problem,            /**< may method be called in the PROBLEM stage? */
   SCIP_Bool             transforming,       /**< may method be called in the TRANSFORMING stage? */
   SCIP_Bool             transformed,        /**< may method be called in the TRANSFORMED stage? */
   SCIP_Bool             initpresolve,       /**< may method be called in the INITPRESOLVE stage? */
   SCIP_Bool             presolving,         /**< may method be called in the PRESOLVING stage? */
   SCIP_Bool             exitpresolve,       /**< may method be called in the EXITPRESOLE stage? */
   SCIP_Bool             presolved,          /**< may method be called in the PRESOLVED stage? */
   SCIP_Bool             initsolve,          /**< may method be called in the INITSOLVE stage? */
   SCIP_Bool             solving,            /**< may method be called in the SOLVING stage? */
   SCIP_Bool             solved,             /**< may method be called in the SOLVED stage? */
   SCIP_Bool             exitsolve,          /**< may method be called in the EXITSOLVE stage? */
   SCIP_Bool             freetrans,          /**< may method be called in the FREETRANS stage? */
   SCIP_Bool             freescip            /**< may method be called in the FREE stage? */
   );
#else

#define SCIPcheckStage(scip,method,init,problem,transforming,transformed,initpresolve,presolving,exitpresolve,presolved, \
   initsolve,solving,solved,exitsolve,freetrans,freescip) SCIP_OKAY

#endif

#ifdef __cplusplus
}
#endif

#endif
