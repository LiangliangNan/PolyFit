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

/**@file    exprinterpret_none.c
 * @brief   function definitions for nonexisting expression interpreter to resolve linking references
 * @ingroup DEFPLUGINS_EXPRINT
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/exprinterpret.h"

/** gets name and version of expression interpreter */
const char* SCIPexprintGetName(
   void
   )
{
   return "NONE";
}  /*lint !e715*/

/** gets descriptive text of expression interpreter */
const char* SCIPexprintGetDesc(
   void
   )
{
   return "dummy expression interpreter which sole purpose it is to resolve linking symbols";
}  /*lint !e715*/

/** gets capabilities of expression interpreter (using bitflags) */
SCIP_EXPRINTCAPABILITY SCIPexprintGetCapability(
   void
   )
{
   return SCIP_EXPRINTCAPABILITY_NONE;
}  /*lint !e715*/

/** creates an expression interpreter object */
SCIP_RETCODE SCIPexprintCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
   )
{
   *exprint = (SCIP_EXPRINT*)1u;  /* some code checks that a non-NULL pointer is returned here, even though it may not point anywhere */

   return SCIP_OKAY;
}  /*lint !e715*/

/** frees an expression interpreter object */
SCIP_RETCODE SCIPexprintFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
   )
{
   *exprint = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

/** compiles an expression and returns interpreter-specific data for expression
 *
 * can be called again with existing exprintdata if expression has been changed
 *
 * @attention *exprintdata needs to be initialized to NULL at first call
 * @attention the expression is assumed to use varidx expressions instead of var expressions
 */
SCIP_RETCODE SCIPexprintCompile(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA**    exprintdata         /**< buffer to store pointer to compiled data */
   )
{
   return SCIP_OKAY;
}  /*lint !e715*/

/** frees interpreter data for expression */
SCIP_RETCODE SCIPexprintFreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA**    exprintdata         /**< pointer to pointer to compiled data to be freed */
   )
{
   assert(exprintdata  != NULL);
   assert(*exprintdata == NULL);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives the capability to evaluate an expression by the expression interpreter
 *
 * In cases of user-given expressions, higher order derivatives may not be available for the user-expression,
 * even if the expression interpreter could handle these. This method allows to recognize that, e.g., the
 * Hessian for an expression is not available because it contains a user expression that does not provide
 * Hessians.
 */
SCIP_EXPRINTCAPABILITY SCIPexprintGetExprCapability(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata         /**< interpreter-specific data for expression */
   )
{
   return SCIP_EXPRINTCAPABILITY_NONE;
} /*lint !e715*/

/** evaluates an expression */
SCIP_RETCODE SCIPexprintEval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value of expression */
   )
{
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}  /*lint !e715*/

/** computes value and gradient of an expression */
SCIP_RETCODE SCIPexprintGrad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient */
   )
{
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}  /*lint !e715*/

/** gives sparsity pattern of lower-triangular part of Hessian
 *
 * Since the AD code might need to do a forward sweep, variable values need to be passed in here.
 *
 * Result will have `(*colidxs)[i] <= (*rowidixs)[i]` for `i=0..*nnz`.
 */
SCIP_RETCODE SCIPexprintHessianSparsity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables */
   int**                 rowidxs,            /**< buffer to return array with row indices of Hessian elements */
   int**                 colidxs,            /**< buffer to return array with column indices of Hessian elements */
   int*                  nnz                 /**< buffer to return length of arrays */
   )
{
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}  /*lint !e715*/

/** computes value and Hessian of an expression
 *
 * Returned arrays `rowidxs` and `colidxs` and number of elements `nnz` are the same as given by SCIPexprintHessianSparsity().
 * Returned array `hessianvals` will contain the corresponding Hessian elements.
 */
SCIP_RETCODE SCIPexprintHessian(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRINTDATA*     exprintdata,        /**< interpreter-specific data for expression */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   int**                 rowidxs,            /**< buffer to return array with row indices of Hessian elements */
   int**                 colidxs,            /**< buffer to return array with column indices of Hessian elements */
   SCIP_Real**           hessianvals,        /**< buffer to return array with Hessian elements */
   int*                  nnz                 /**< buffer to return length of arrays */
   )
{
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}  /*lint !e715*/
