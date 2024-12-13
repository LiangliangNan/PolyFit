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

/**@file   struct_nlpi.h
 * @ingroup INTERNALAPI
 * @brief  data definitions for an NLP solver interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_NLPI_H__
#define __SCIP_STRUCT_NLPI_H__

#include "scip/def.h"
#include "scip/type_nlpi.h"
#include "scip/type_clock.h"

#ifdef __cplusplus
extern "C" {
#endif

/** NLP interface data */
struct SCIP_Nlpi
{
   char*                           name;                        /**< name of NLP solver */
   char*                           description;                 /**< description of NLP solver */
   int                             priority;                    /**< priority of NLP interface */
   SCIP_DECL_NLPICOPY              ((*nlpicopy));               /**< copy an NLPI */
   SCIP_DECL_NLPIFREE              ((*nlpifree));               /**< free NLPI user data */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer));   /**< get solver pointer */
   SCIP_DECL_NLPICREATEPROBLEM     ((*nlpicreateproblem));      /**< create a new problem instance */
   SCIP_DECL_NLPIFREEPROBLEM       ((*nlpifreeproblem));        /**< free a problem instance */
   SCIP_DECL_NLPIGETPROBLEMPOINTER ((*nlpigetproblempointer));  /**< get problem pointer */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars));            /**< add variables to a problem */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints));     /**< add constraints to a problem  */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective));       /**< set objective of a problem  */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds));       /**< change variable bounds in a problem  */
   SCIP_DECL_NLPICHGCONSSIDES      ((*nlpichgconssides));       /**< change constraint sides in a problem  */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset));          /**< delete a set of variables from a problem  */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset));         /**< delete a set of constraints from a problem */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs));     /**< change coefficients in linear part of a constraint or objective */
   SCIP_DECL_NLPICHGEXPR           ((*nlpichgexpr));            /**< change nonlinear expression a constraint or objective */
   SCIP_DECL_NLPICHGOBJCONSTANT    ((*nlpichgobjconstant));     /**< change the constant offset in the objective */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess));    /**< set initial guess */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve));              /**< solve a problem */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat));         /**< get solution status for a problem  */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat));        /**< get termination status for a problem  */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution));        /**< get solution of a problem  */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics));      /**< get solve statistics for a problem  */
   SCIP_NLPIDATA*                  nlpidata;                    /**< NLP interface local data */

   /* statistics */
   int                             nproblems;                   /**< number of problems created */
   int                             nsolves;                     /**< number of solves */
   SCIP_CLOCK*                     problemtime;                 /**< time spend in problem setup and modification */
   SCIP_Real                       solvetime;                   /**< time spend in solve as reported by solver */
   SCIP_Real                       evaltime;                    /**< time spend in function evaluation during solve */
   SCIP_Longint                    niter;                       /**< total number of iterations */
   int                             ntermstat[SCIP_NLPTERMSTAT_OTHER+1]; /**< number of times a specific termination status occurred */
   int                             nsolstat[SCIP_NLPSOLSTAT_UNKNOWN+1]; /**< number of times a specific solution status occurred */
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_NLPI_H__ */
