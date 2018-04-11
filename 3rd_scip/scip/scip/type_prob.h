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

/**@file   type_prob.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_PROB_H__
#define __SCIP_TYPE_PROB_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** objective sense: minimization or maximization */
enum SCIP_Objsense
{
   SCIP_OBJSENSE_MAXIMIZE = -1,         /**< maximization of objective function */
   SCIP_OBJSENSE_MINIMIZE = +1          /**< minimization of objective function (the default) */
};
typedef enum SCIP_Objsense SCIP_OBJSENSE;

typedef struct SCIP_Prob SCIP_PROB;               /**< main problem to solve */
typedef struct SCIP_ProbData SCIP_PROBDATA;       /**< user problem data set by the reader */


/** frees user data of original problem (called when the original problem is freed)
 *
 *  This method should free the user data of the original problem.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - probdata        : pointer to the user problem data to free
 */
#define SCIP_DECL_PROBDELORIG(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROBDATA** probdata)

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 *
 *  Because the original problem and the user data of the original problem should not be
 *  modified during the solving process, a transformed problem is created as a copy of
 *  the original problem. If the user problem data is never modified during the solving
 *  process anyways, it is enough to simple copy the user data's pointer. This is the
 *  default implementation, which is used when a NULL is given as PROBTRANS method.
 *  If the user data may be modified during the solving process (e.g. during preprocessing),
 *  the PROBTRANS method must be given and has to copy the user problem data to a different
 *  memory location.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sourcedata      : source problem data to transform
 *  - targetdata      : pointer to store created transformed problem data
 */
#define SCIP_DECL_PROBTRANS(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROBDATA* sourcedata, SCIP_PROBDATA** targetdata)

/** frees user data of transformed problem (called when the transformed problem is freed)
 *
 *  This method has to be implemented, if the PROBTRANS method is not a simple pointer
 *  copy operation like in the default PROBTRANS implementation. It should free the
 *  user data of the transformed problem, that was created in the PROBTRANS method.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - probdata        : pointer to the user problem data to free
 */
#define SCIP_DECL_PROBDELTRANS(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROBDATA** probdata)

/** solving process initialization method of transformed data (called before the branch and bound process begins)
 *
 *  This method is called before the branch and bound process begins and can be used to initialize user problem
 *  data that depends for example on the number of active problem variables, because these are now fixed.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - probdata        : user problem data
 */
#define SCIP_DECL_PROBINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROBDATA* probdata)

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed)
 *
 *  This method is called before the branch and bound data is freed and should be used to free all data that
 *  was allocated in the solving process initialization method. The user has to make sure, that all LP rows associated
 *  to the transformed user problem data are released.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - probdata        : user problem data
 *  - restart         : was this exit solve call triggered by a restart?
 */
#define SCIP_DECL_PROBEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_PROBDATA* probdata, SCIP_Bool restart)

/** copies user data of source SCIP for the target SCIP
 *
 *  This method should copy the problem data of the source SCIP and create a target problem data for (target)
 *  SCIP. Implementing this callback is optional. If the copying process was successful the target SCIP gets this
 *  problem data assigned. In case the result pointer is set to SCIP_DIDNOTRUN the target SCIP will have no problem data
 *  at all.
 *
 *  The variable map and the constraint map can be used via the function SCIPgetVarCopy() and SCIPgetConsCopy(),
 *  respectively, to get for certain variables and constraints of the source SCIP the counter parts in the target
 *  SCIP. You should be very carefully in using these two methods since they could lead to an infinite loop due to
 *  recursion.
 *
 *  input:
 *  - scip            : target SCIP data structure
 *  - sourcescip      : source SCIP main data structure
 *  - sourcedata      : source user problem data
 *  - varmap,         : a hashmap which stores the mapping of source variables to corresponding target variables
 *  - consmap,        : a hashmap which stores the mapping of source constraints to corresponding target constraints
 *  - targetdata      : pointer to the target user problem data to create
 *  - global          : create a global or a local copy? 
 *
 *  output:
 *  - result          : pointer to store the result of the call
 *
 *  possible return values for *result:
 *  - SCIP_DIDNOTRUN  : the copying process was not performed 
 *  - SCIP_SUCCESS    : the copying process was successfully performed
 */
#define SCIP_DECL_PROBCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP* sourcescip, SCIP_PROBDATA* sourcedata, \
      SCIP_HASHMAP* varmap, SCIP_HASHMAP* consmap, SCIP_PROBDATA** targetdata, SCIP_Bool global, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
