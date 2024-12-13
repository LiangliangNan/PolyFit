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

/**@file   scip_var.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for SCIP variables
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

#include "blockmemshell/memory.h"
#include "lpi/lpi.h"
#include "scip/branch.h"
#include "scip/clock.h"
#include "scip/conflict.h"
#include "scip/debug.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/relax.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/struct_lp.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"
#include <ctype.h>


/** creates and captures problem variable; if variable is of integral type, fractional bounds are automatically rounded;
 *  an integer variable with bounds zero and one is automatically converted into a binary variable;
 *
 *  @warning When doing column generation and the original problem is a maximization problem, notice that SCIP will
 *           transform the problem into a minimization problem by multiplying the objective function by -1.  Thus, the
 *           original objective function value of variables created during the solving process has to be multiplied by
 *           -1, too.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the variable gets captured, hence at one point you have to release it using the method SCIPreleaseVar()
 */
SCIP_RETCODE SCIPcreateVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(lb <= ub);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* forbid infinite objective function values */
   if( SCIPisInfinity(scip, REALABS(obj)) )
   {
      SCIPerrorMessage("invalid objective function value: value is infinite\n");
      return SCIP_INVALIDDATA;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPvarCreateOriginal(var, scip->mem->probmem, scip->set, scip->stat,
            name, lb, ub, obj, vartype, initial, removable, vardelorig, vartrans, vardeltrans, varcopy, vardata) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPvarCreateTransformed(var, scip->mem->probmem, scip->set, scip->stat,
            name, lb, ub, obj, vartype, initial, removable, vardelorig, vartrans, vardeltrans, varcopy, vardata) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** creates and captures problem variable with optional callbacks and variable data set to NULL, which can be set
 *  afterwards using SCIPvarSetDelorigData(), SCIPvarSetTransData(),
 *  SCIPvarSetDeltransData(), SCIPvarSetCopy(), and SCIPvarSetData(); sets variable flags initial=TRUE
 *  and removable = FALSE, which can be adjusted by using SCIPvarSetInitial() and SCIPvarSetRemovable(), resp.;
 *  if variable is of integral type, fractional bounds are automatically rounded;
 *  an integer variable with bounds zero and one is automatically converted into a binary variable;
 *
 *  @warning When doing column generation and the original problem is a maximization problem, notice that SCIP will
 *           transform the problem into a minimization problem by multiplying the objective function by -1.  Thus, the
 *           original objective function value of variables created during the solving process has to be multiplied by
 *           -1, too.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the variable gets captured, hence at one point you have to release it using the method SCIPreleaseVar()
 */
SCIP_RETCODE SCIPcreateVarBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype             /**< type of variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateVarBasic", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcreateVar(scip, var, name, lb, ub, obj, vartype, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** outputs the variable name to the file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPwriteVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR*             var,                /**< variable to output */
   SCIP_Bool             type                /**< should the variable type be also posted */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteVarName", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* print variable name */
   if( SCIPvarIsNegated(var) )
   {
      SCIP_VAR* negatedvar;

      SCIP_CALL( SCIPgetNegatedVar(scip, var, &negatedvar) );
      SCIPinfoMessage(scip, file, "<~%s>", SCIPvarGetName(negatedvar));
   }
   else
   {
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(var));
   }

   if( type )
   {
      /* print variable type */
      SCIPinfoMessage(scip, file, "[%c]",
         SCIPvarGetType(var) == SCIP_VARTYPE_BINARY ? SCIP_VARTYPE_BINARY_CHAR :
         SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER ? SCIP_VARTYPE_INTEGER_CHAR :
         SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT ? SCIP_VARTYPE_IMPLINT_CHAR : SCIP_VARTYPE_CONTINUOUS_CHAR);
   }

   return SCIP_OKAY;
}

/** print the given list of variables to output stream separated by the given delimiter character;
 *
 *  i. e. the variables x1, x2, ..., xn with given delimiter ',' are written as: \<x1\>, \<x2\>, ..., \<xn\>;
 *
 *  the method SCIPparseVarsList() can parse such a string
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The printing process is done via the message handler system.
 */
SCIP_RETCODE SCIPwriteVarsList(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR**            vars,               /**< variable array to output */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             type,               /**< should the variable type be also posted */
   char                  delimiter           /**< character which is used for delimitation */
   )
{
   int v;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteVarsList", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( v > 0 )
      {
         SCIPinfoMessage(scip, file, "%c", delimiter);
      }

      /* print variable name */
      SCIP_CALL( SCIPwriteVarName(scip, file, vars[v], type) );
   }

   return SCIP_OKAY;
}

/** print the given variables and coefficients as linear sum in the following form
 *  c1 \<x1\> + c2 \<x2\>   ... + cn \<xn\>
 *
 *  This string can be parsed by the method SCIPparseVarsLinearsum().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The printing process is done via the message handler system.
 */
SCIP_RETCODE SCIPwriteVarsLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR**            vars,               /**< variable array to output */
   SCIP_Real*            vals,               /**< array of coefficients or NULL if all coefficients are 1.0 */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             type                /**< should the variable type be also posted */
   )
{
   int v;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteVarsLinearsum", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( vals != NULL )
      {
         if( vals[v] == 1.0 )
         {
            if( v > 0 )
               SCIPinfoMessage(scip, file, " +");
         }
         else if( vals[v] == -1.0 )
            SCIPinfoMessage(scip, file, " -");
         else
            SCIPinfoMessage(scip, file, " %+.15g", vals[v]);
      }
      else if( nvars > 0 )
         SCIPinfoMessage(scip, file, " +");

      /* print variable name */
      SCIP_CALL( SCIPwriteVarName(scip, file, vars[v], type) );
   }

   return SCIP_OKAY;
}

/** print the given terms as signomial in the following form
 *  c1 \<x11\>^e11 \<x12\>^e12 ... \<x1n\>^e1n + c2 \<x21\>^e21 \<x22\>^e22 ... + ... + cn \<xn1\>^en1 ...
 *
 *  This string can be parsed by the method SCIPparseVarsPolynomial().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The printing process is done via the message handler system.
 */
SCIP_RETCODE SCIPwriteVarsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR***           monomialvars,       /**< arrays with variables for each monomial */
   SCIP_Real**           monomialexps,       /**< arrays with variable exponents, or NULL if always 1.0 */
   SCIP_Real*            monomialcoefs,      /**< array with monomial coefficients */
   int*                  monomialnvars,      /**< array with number of variables for each monomial */
   int                   nmonomials,         /**< number of monomials */
   SCIP_Bool             type                /**< should the variable type be also posted */
   )
{
   int i;
   int v;

   assert(scip != NULL);
   assert(monomialvars  != NULL || nmonomials == 0);
   assert(monomialcoefs != NULL || nmonomials == 0);
   assert(monomialnvars != NULL || nmonomials == 0);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteVarsPolynomial", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( nmonomials == 0 )
   {
      SCIPinfoMessage(scip, file, " 0 ");
      return SCIP_OKAY;
   }

   for( i = 0; i < nmonomials; ++i )
   {
      if( monomialcoefs[i] == 1.0 ) /*lint !e613*/
      {
         if( i > 0 )
            SCIPinfoMessage(scip, file, " +");
      }
      else if( monomialcoefs[i] == -1.0 ) /*lint !e613*/
         SCIPinfoMessage(scip, file, " -");
      else
         SCIPinfoMessage(scip, file, " %+.15g", monomialcoefs[i]); /*lint !e613*/

      assert(monomialvars[i] != NULL || monomialnvars[i] == 0); /*lint !e613*/

      for( v = 0; v < monomialnvars[i]; ++v ) /*lint !e613*/
      {
         SCIP_CALL( SCIPwriteVarName(scip, file, monomialvars[i][v], type) ); /*lint !e613*/
         if( monomialexps != NULL && monomialexps[i] != NULL && monomialexps[i][v] != 1.0 )
         {
            SCIPinfoMessage(scip, file, "^%.15g", monomialexps[i][v]);
         }
      }
   }

   return SCIP_OKAY;
}

/** parses variable information (in cip format) out of a string; if the parsing process was successful a variable is
 *  created and captured; if variable is of integral type, fractional bounds are automatically rounded; an integer
 *  variable with bounds zero and one is automatically converted into a binary variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPparseVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to store the problem variable */
   const char*           str,                /**< string to parse */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata,            /**< user data for this specific variable */
   char**                endptr,             /**< pointer to store the final string position if successful */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   )
{
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPparseVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPvarParseOriginal(var, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
            str, initial, removable, varcopy, vardelorig, vartrans, vardeltrans, vardata, endptr, success) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPvarParseTransformed(var, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
            str, initial, removable, varcopy, vardelorig, vartrans, vardeltrans, vardata, endptr, success) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** parses the given string for a variable name and stores the variable in the corresponding pointer if such a variable
 *  exits and returns the position where the parsing stopped
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPparseVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            var,                /**< pointer to store the problem variable, or NULL if it does not exit */
   char**                endptr              /**< pointer to store the final string position if successful */
   )
{
   char varname[SCIP_MAXSTRLEN];

   assert(str != NULL);
   assert(var != NULL);
   assert(endptr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPparseVarName", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPstrCopySection(str, '<', '>', varname, SCIP_MAXSTRLEN, endptr);
   assert(*endptr != NULL);

   if( *varname == '\0' )
   {
      SCIPerrorMessage("invalid variable name string given: could not find '<'\n");
      return SCIP_INVALIDDATA;
   }

   /* check if we have a negated variable */
   if( *varname == '~' )
   {
      SCIPdebugMsg(scip, "parsed negated variable name <%s>\n", &varname[1]);

      /* search for the variable and ignore '~' */
      (*var) = SCIPfindVar(scip, &varname[1]);

      if( *var  != NULL )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, *var, var) );
      }
   }
   else
   {
      SCIPdebugMsg(scip, "parsed variable name <%s>\n", varname);

      /* search for the variable */
      (*var) = SCIPfindVar(scip, varname);
   }

   str = *endptr;

   /* skip additional variable type marker */
   if( *str == '[' && (str[1] == SCIP_VARTYPE_BINARY_CHAR || str[1] == SCIP_VARTYPE_INTEGER_CHAR ||
       str[1] == SCIP_VARTYPE_IMPLINT_CHAR || str[1] == SCIP_VARTYPE_CONTINUOUS_CHAR )  && str[2] == ']' )
      (*endptr) += 3;

   return SCIP_OKAY;
}

/** parse the given string as variable list (here ',' is the delimiter)) (\<x1\>, \<x2\>, ..., \<xn\>) (see
 *  SCIPwriteVarsList() ); if it was successful, the pointer success is set to TRUE
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist.
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
SCIP_RETCODE SCIPparseVarsList(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variable */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successful */
   char                  delimiter,          /**< character which is used for delimitation */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successful or not */
   )
{
   SCIP_VAR** tmpvars;
   SCIP_VAR* var;
   int ntmpvars = 0;
   int v;

   assert( nvars != NULL );
   assert( requiredsize != NULL );
   assert( endptr != NULL );
   assert( success != NULL );

   SCIP_CALL( SCIPcheckStage(scip, "SCIPparseVarsList", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* allocate buffer memory for temporary storing the parsed variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, varssize) );

   (*success) = TRUE;

   do
   {
      *endptr = (char*)str;

      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, str, &var, endptr) );

      if( var == NULL )
      {
         SCIPdebugMsg(scip, "variable with name <%s> does not exist\n", SCIPvarGetName(var));
         (*success) = FALSE;
         break;
      }

      /* store the variable in the tmp array */
      if( ntmpvars < varssize )
         tmpvars[ntmpvars] = var;

      ntmpvars++;

      str = *endptr;

      while( isspace((unsigned char)*str) )
         str++;
   }
   while( *str == delimiter );

   *endptr = (char*)str;

   /* if all variable name searches were successful and the variable array has enough slots, copy the collected variables */
   if( (*success) && ntmpvars <= varssize )
   {
      for( v = 0; v < ntmpvars; ++v )
         vars[v] = tmpvars[v];

      (*nvars) = ntmpvars;
   }
   else
      (*nvars) = 0;

   (*requiredsize) = ntmpvars;

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &tmpvars);

   return SCIP_OKAY;
}

/** parse the given string as linear sum of variables and coefficients (c1 \<x1\> + c2 \<x2\> + ... + cn \<xn\>)
 *  (see SCIPwriteVarsLinearsum() ); if it was successful, the pointer success is set to TRUE
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note The pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist.
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
SCIP_RETCODE SCIPparseVarsLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variables */
   SCIP_Real*            vals,               /**< array to store the parsed coefficients */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successful */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successful or not */
   )
{
   SCIP_VAR*** monomialvars;
   SCIP_Real** monomialexps;
   SCIP_Real*  monomialcoefs;
   int*        monomialnvars;
   int         nmonomials;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPparseVarsLinearsum", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(scip != NULL);
   assert(str != NULL);
   assert(vars != NULL || varssize == 0);
   assert(vals != NULL || varssize == 0);
   assert(nvars != NULL);
   assert(requiredsize != NULL);
   assert(endptr != NULL);
   assert(success != NULL);

   *requiredsize = 0;

   SCIP_CALL( SCIPparseVarsPolynomial(scip, str, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, &nmonomials, endptr, success) );

   if( !*success )
   {
      assert(nmonomials == 0); /* SCIPparseVarsPolynomial should have freed all buffers, so no need to call free here */
      return SCIP_OKAY;
   }

   /* check if linear sum is just "0" */
   if( nmonomials == 1 && monomialnvars[0] == 0 && monomialcoefs[0] == 0.0 )
   {
      *nvars = 0;
      *requiredsize = 0;

      SCIPfreeParseVarsPolynomialData(scip, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, nmonomials);

      return SCIP_OKAY;
   }

   *nvars = nmonomials;
   *requiredsize = nmonomials;

   /* if we have enough slots in the variables array, copy variables over */
   if( varssize >= nmonomials )
   {
      int v;

      for( v = 0; v < nmonomials; ++v )
      {
         if( monomialnvars[v] == 0 )
         {
            SCIPerrorMessage("constant in linear sum\n");
            *success = FALSE;
            break;
         }
         if( monomialnvars[v] > 1 || monomialexps[v][0] != 1.0 )
         {
            SCIPerrorMessage("nonlinear monomial in linear sum\n");
            *success = FALSE;
            break;
         }
         assert(monomialnvars[v]   == 1);
         assert(monomialvars[v][0] != NULL);
         assert(monomialexps[v][0] == 1.0);

         vars[v] = monomialvars[v][0]; /*lint !e613*/
         vals[v] = monomialcoefs[v]; /*lint !e613*/
      }
   }

   SCIPfreeParseVarsPolynomialData(scip, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, nmonomials);

   return SCIP_OKAY;
}

/** parse the given string as signomial of variables and coefficients
 *  (c1 \<x11\>^e11 \<x12\>^e12 ... \<x1n\>^e1n + c2 \<x21\>^e21 \<x22\>^e22 ... + ... + cn \<xn1\>^en1 ...)
 *  (see SCIPwriteVarsPolynomial()); if it was successful, the pointer success is set to TRUE
 *
 *  The user has to call SCIPfreeParseVarsPolynomialData(scip, monomialvars, monomialexps,
 *  monomialcoefs, monomialnvars, *nmonomials) short after SCIPparseVarsPolynomial to free all the
 *  allocated memory again.
 *
 *  Parsing is stopped at the end of string (indicated by the \\0-character) or when no more monomials
 *  are recognized.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPparseVarsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int*                  nmonomials,         /**< pointer to store number of parsed monomials */
   char**                endptr,             /**< pointer to store the final string position if successful */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successful or not */
   )
{
   typedef enum
   {
      SCIPPARSEPOLYNOMIAL_STATE_BEGIN,       /* we are at the beginning of a monomial */
      SCIPPARSEPOLYNOMIAL_STATE_INTERMED,    /* we are in between the factors of a monomial */
      SCIPPARSEPOLYNOMIAL_STATE_COEF,        /* we parse the coefficient of a monomial */
      SCIPPARSEPOLYNOMIAL_STATE_VARS,        /* we parse monomial variables */
      SCIPPARSEPOLYNOMIAL_STATE_EXPONENT,    /* we parse the exponent of a variable */
      SCIPPARSEPOLYNOMIAL_STATE_END,         /* we are at the end the polynomial */
      SCIPPARSEPOLYNOMIAL_STATE_ERROR        /* a parsing error occured */
   } SCIPPARSEPOLYNOMIAL_STATES;

   SCIPPARSEPOLYNOMIAL_STATES state;
   int monomialssize;

   /* data of currently parsed monomial */
   int varssize;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Real* exponents;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(str != NULL);
   assert(monomialvars != NULL);
   assert(monomialexps != NULL);
   assert(monomialnvars != NULL);
   assert(monomialcoefs != NULL);
   assert(nmonomials != NULL);
   assert(endptr != NULL);
   assert(success != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPparseVarsPolynomial", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *success = FALSE;
   *nmonomials = 0;
   monomialssize = 0;
   *monomialvars = NULL;
   *monomialexps = NULL;
   *monomialcoefs = NULL;
   *monomialnvars = NULL;

   /* initialize state machine */
   state = SCIPPARSEPOLYNOMIAL_STATE_BEGIN;
   varssize = 0;
   nvars = 0;
   vars = NULL;
   exponents = NULL;
   coef = SCIP_INVALID;

   SCIPdebugMsg(scip, "parsing polynomial from '%s'\n", str);

   while( *str && state != SCIPPARSEPOLYNOMIAL_STATE_END && state != SCIPPARSEPOLYNOMIAL_STATE_ERROR )
   {
      /* skip white space */
      while( isspace((unsigned char)*str) )
         str++;

      assert(state != SCIPPARSEPOLYNOMIAL_STATE_END);

      switch( state )
      {
      case SCIPPARSEPOLYNOMIAL_STATE_BEGIN:
      {
         if( coef != SCIP_INVALID  ) /*lint !e777*/
         {
            SCIPdebugMsg(scip, "push monomial with coefficient <%g> and <%d> vars\n", coef, nvars);

            /* push previous monomial */
            if( monomialssize <= *nmonomials )
            {
               monomialssize = SCIPcalcMemGrowSize(scip, *nmonomials+1);

               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialvars,  *nmonomials, monomialssize) );
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialexps,  *nmonomials, monomialssize) );
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialnvars, *nmonomials, monomialssize) );
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialcoefs, *nmonomials, monomialssize) );
            }

            if( nvars > 0 )
            {
               SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*monomialvars)[*nmonomials], vars, nvars) ); /*lint !e866*/
               SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*monomialexps)[*nmonomials], exponents, nvars) ); /*lint !e866*/
            }
            else
            {
               (*monomialvars)[*nmonomials] = NULL;
               (*monomialexps)[*nmonomials] = NULL;
            }
            (*monomialcoefs)[*nmonomials] = coef;
            (*monomialnvars)[*nmonomials] = nvars;
            ++*nmonomials;

            nvars = 0;
            coef = SCIP_INVALID;
         }

         if( *str == '<' )
         {
            /* there seem to come a variable at the beginning of a monomial
             * so assume the coefficient is 1.0
             */
            state = SCIPPARSEPOLYNOMIAL_STATE_VARS;
            coef = 1.0;
         }
         else if( *str == '-' || *str == '+' || isdigit(*str) )
            state = SCIPPARSEPOLYNOMIAL_STATE_COEF;
         else
            state = SCIPPARSEPOLYNOMIAL_STATE_END;

         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_INTERMED:
      {
         if( *str == '<' )
         {
            /* there seem to come another variable */
            state = SCIPPARSEPOLYNOMIAL_STATE_VARS;
         }
         else if( *str == '-' || *str == '+' || isdigit(*str) )
         {
            /* there seem to come a coefficient, which means the next monomial */
            state = SCIPPARSEPOLYNOMIAL_STATE_BEGIN;
         }
         else /* since we cannot detect the symbols we stop parsing the polynomial */
            state = SCIPPARSEPOLYNOMIAL_STATE_END;

         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_COEF:
      {
         if( *str == '+' && !isdigit(str[1]) )
         {
            /* only a plus sign, without number */
            coef =  1.0;
            ++str;
         }
         else if( *str == '-' && !isdigit(str[1]) )
         {
            /* only a minus sign, without number */
            coef = -1.0;
            ++str;
         }
         else if( SCIPstrToRealValue(str, &coef, endptr) )
         {
            str = *endptr;
         }
         else
         {
            SCIPerrorMessage("could not parse number in the beginning of '%s'\n", str);
            state = SCIPPARSEPOLYNOMIAL_STATE_ERROR;
            break;
         }

         /* after the coefficient we go into the intermediate state, i.e., expecting next variables */
         state = SCIPPARSEPOLYNOMIAL_STATE_INTERMED;  /*lint !e838*/

         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_VARS:
      {
         SCIP_VAR* var;

         assert(*str == '<');

         /* parse variable name */
         SCIP_CALL( SCIPparseVarName(scip, str, &var, endptr) );

         /* check if variable name was parsed */
         if( *endptr == str )
         {
            state = SCIPPARSEPOLYNOMIAL_STATE_END;
            break;
         }

         if( var == NULL )
         {
            SCIPerrorMessage("did not find variable in the beginning of %s\n", str);
            state = SCIPPARSEPOLYNOMIAL_STATE_ERROR;
            break;
         }

         /* add variable to vars array */
         if( nvars + 1 > varssize )
         {
            varssize = SCIPcalcMemGrowSize(scip, nvars+1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &vars,      nvars, varssize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &exponents, nvars, varssize) );
         }
         assert(vars != NULL);
         assert(exponents != NULL);

         vars[nvars] = var;
         exponents[nvars] = 1.0;
         ++nvars;

         str = *endptr;

         if( *str == '^' )
            state = SCIPPARSEPOLYNOMIAL_STATE_EXPONENT;
         else
            state = SCIPPARSEPOLYNOMIAL_STATE_INTERMED;

         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_EXPONENT:
      {
         assert(*str == '^');
         assert(nvars > 0); /* we should be in a monomial that has already a variable */
         assert(exponents != NULL);
         ++str;

         if( !SCIPstrToRealValue(str, &exponents[nvars-1], endptr) )
         {
            SCIPerrorMessage("could not parse number in the beginning of '%s'\n", str);
            state = SCIPPARSEPOLYNOMIAL_STATE_ERROR;
            break;
         }
         str = *endptr;

         /* after the exponent we go into the intermediate state, i.e., expecting next variables */
         state = SCIPPARSEPOLYNOMIAL_STATE_INTERMED;  /*lint !e838*/
         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_END:
      case SCIPPARSEPOLYNOMIAL_STATE_ERROR:
      default:
         SCIPerrorMessage("unexpected state\n");
         return SCIP_READERROR;
      }
   }

   /* set end pointer */
   *endptr = (char*)str;

   /* check state at end of string */
   switch( state )
   {
   case SCIPPARSEPOLYNOMIAL_STATE_BEGIN:
   case SCIPPARSEPOLYNOMIAL_STATE_END:
   case SCIPPARSEPOLYNOMIAL_STATE_INTERMED:
   {
      if( coef != SCIP_INVALID ) /*lint !e777*/
      {
         /* push last monomial */
         SCIPdebugMsg(scip, "push monomial with coefficient <%g> and <%d> vars\n", coef, nvars);
         if( monomialssize <= *nmonomials )
         {
            monomialssize = *nmonomials+1;
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialvars,  *nmonomials, monomialssize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialexps,  *nmonomials, monomialssize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialnvars, *nmonomials, monomialssize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialcoefs, *nmonomials, monomialssize) );
         }

         if( nvars > 0 )
         {
            /* shrink vars and exponents array to needed size and take over ownership */
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &vars, varssize, nvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &exponents, varssize, nvars) );
            (*monomialvars)[*nmonomials] = vars;
            (*monomialexps)[*nmonomials] = exponents;
            vars = NULL;
            exponents = NULL;
         }
         else
         {
            (*monomialvars)[*nmonomials] = NULL;
            (*monomialexps)[*nmonomials] = NULL;
         }
         (*monomialcoefs)[*nmonomials] = coef;
         (*monomialnvars)[*nmonomials] = nvars;
         ++*nmonomials;
      }

      *success = TRUE;
      break;
   }

   case SCIPPARSEPOLYNOMIAL_STATE_COEF:
   case SCIPPARSEPOLYNOMIAL_STATE_VARS:
   case SCIPPARSEPOLYNOMIAL_STATE_EXPONENT:
   {
      SCIPerrorMessage("unexpected parsing state at end of polynomial string\n");
   }
   /*lint -fallthrough*/
   case SCIPPARSEPOLYNOMIAL_STATE_ERROR:
      assert(!*success);
      break;
   }

   /* free memory to store current monomial, if still existing */
   SCIPfreeBlockMemoryArrayNull(scip, &vars, varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &exponents, varssize);

   if( *success && *nmonomials > 0 )
   {
      /* shrink arrays to required size, so we do not need to keep monomialssize around */
      assert(*nmonomials <= monomialssize);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialvars,  monomialssize, *nmonomials) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialexps,  monomialssize, *nmonomials) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialnvars, monomialssize, *nmonomials) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, monomialcoefs, monomialssize, *nmonomials) );

      /* SCIPwriteVarsPolynomial(scip, NULL, *monomialvars, *monomialexps, *monomialcoefs, *monomialnvars, *nmonomials, FALSE); */
   }
   else
   {
      /* in case of error, cleanup all data here */
      SCIPfreeParseVarsPolynomialData(scip, monomialvars, monomialexps, monomialcoefs, monomialnvars, *nmonomials);
      *nmonomials = 0;
   }

   return SCIP_OKAY;
}

/** frees memory allocated when parsing a signomial from a string
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
void SCIPfreeParseVarsPolynomialData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int                   nmonomials          /**< pointer to store number of parsed monomials */
   )
{
   int i;

   assert(scip != NULL);
   assert(monomialvars  != NULL);
   assert(monomialexps  != NULL);
   assert(monomialcoefs != NULL);
   assert(monomialnvars != NULL);
   assert((*monomialvars  != NULL) == (nmonomials > 0));
   assert((*monomialexps  != NULL) == (nmonomials > 0));
   assert((*monomialcoefs != NULL) == (nmonomials > 0));
   assert((*monomialnvars != NULL) == (nmonomials > 0));

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPfreeParseVarsPolynomialData", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( nmonomials == 0 )
      return;

   for( i = nmonomials - 1; i >= 0; --i )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*monomialexps)[i], (*monomialnvars)[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*monomialvars)[i], (*monomialnvars)[i]);
   }

   SCIPfreeBlockMemoryArray(scip, monomialcoefs, nmonomials);
   SCIPfreeBlockMemoryArray(scip, monomialnvars, nmonomials);
   SCIPfreeBlockMemoryArray(scip, monomialexps, nmonomials);
   SCIPfreeBlockMemoryArray(scip, monomialvars, nmonomials);
}

/** increases usage counter of variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPcaptureVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to capture */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcaptureVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   assert(var->scip == scip);

   SCIPvarCapture(var);

   return SCIP_OKAY;
}

/** decreases usage counter of variable, if the usage pointer reaches zero the variable gets freed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note the pointer of the variable will be NULLed
 */
SCIP_RETCODE SCIPreleaseVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var                 /**< pointer to variable */
   )
{
   assert(var != NULL);
   assert(*var != NULL);
   assert((*var)->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPreleaseVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      if( !SCIPvarIsTransformed(*var) && (*var)->nuses == 1 && (*var)->data.original.transvar != NULL )
      {
         SCIPerrorMessage("cannot release last use of original variable while associated transformed variable exists\n");
         return SCIP_INVALIDCALL;
      }
      SCIP_CALL( SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** changes the name of a variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_PROBLEM
 *
 *  @note to get the current name of a variable, use SCIPvarGetName() from pub_var.h
 */
SCIP_RETCODE SCIPchgVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable */
   const char*           name                /**< new name of constraint */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarName", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   assert( var->scip == scip );

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("variable names can only be changed in problem creation stage\n");
      SCIPABORT();
      return SCIP_INVALIDCALL; /*lint !e527*/
   }

   /* remove variable's name from the namespace if the variable was already added */
   if( SCIPvarGetProbindex(var) != -1 )
   {
      SCIP_CALL( SCIPprobRemoveVarName(scip->origprob, var) );
   }

   /* change variable name */
   SCIP_CALL( SCIPvarChgName(var, SCIPblkmem(scip), name) );

   /* add variable's name to the namespace if the variable was already added */
   if( SCIPvarGetProbindex(var) != -1 )
   {
      SCIP_CALL( SCIPprobAddVarName(scip->origprob, var) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPtransformVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get/create transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   assert(transvar != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtransformVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPvarIsTransformed(var) )
   {
      *transvar = var;
      SCIPvarCapture(*transvar);
   }
   else
   {
      SCIP_CALL( SCIPvarTransform(var, scip->mem->probmem, scip->set, scip->stat, scip->origprob->objsense, transvar) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed variables for an array of variables;
 *  if a variable of the array is not yet transformed, a new transformed variable for this variable is created;
 *  it is possible to call this method with vars == transvars
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPtransformVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get/create transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get/create transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || transvars != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtransformVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarIsTransformed(vars[v]) )
      {
         transvars[v] = vars[v];
         SCIPvarCapture(transvars[v]);
      }
      else
      {
         SCIP_CALL( SCIPvarTransform(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->origprob->objsense,
               &transvars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed variable of a given variable;
 *  returns NULL as transvar, if transformed variable is not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetTransformedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   assert(transvar != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetTransformedVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPvarIsTransformed(var) )
      *transvar = var;
   else
   {
      SCIP_CALL( SCIPvarGetTransformed(var, scip->mem->probmem, scip->set, scip->stat, transvar) );
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed variables for an array of variables;
 *  stores NULL in a transvars slot, if the transformed variable is not yet existing;
 *  it is possible to call this method with vars == transvars, but remember that variables that are not
 *  yet transformed will be replaced with NULL
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetTransformedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || transvars != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetTransformedVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarIsTransformed(vars[v]) )
         transvars[v] = vars[v];
      else
      {
         SCIP_CALL( SCIPvarGetTransformed(vars[v], scip->mem->probmem, scip->set, scip->stat, &transvars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** gets negated variable x' = lb + ub - x of variable x; negated variable is created, if not yet existing;
 *  in difference to \ref SCIPcreateVar, the negated variable must not be released (unless captured explicitly)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetNegatedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get negated variable for */
   SCIP_VAR**            negvar              /**< pointer to store the negated variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNegatedVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   assert( var->scip == scip );

   SCIP_CALL( SCIPvarNegate(var, scip->mem->probmem, scip->set, scip->stat, negvar) );

   return SCIP_OKAY;
}

/** gets negated variables x' = lb + ub - x of variables x; negated variables are created, if not yet existing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetNegatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get negated variables for */
   SCIP_VAR**            vars,               /**< array of variables to get negated variables for */
   SCIP_VAR**            negvars             /**< array to store the negated variables */
   )
{
   int v;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNegatedVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarNegate(vars[v], scip->mem->probmem, scip->set, scip->stat, &(negvars[v])) );
   }

   return SCIP_OKAY;
}

/** gets a binary variable that is equal to the given binary variable, and that is either active, fixed, or
 *  multi-aggregated, or the negated variable of an active, fixed, or multi-aggregated variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPgetBinvarRepresentative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to get binary representative for */
   SCIP_VAR**            repvar,             /**< pointer to store the binary representative */
   SCIP_Bool*            negated             /**< pointer to store whether the negation of an active variable was returned */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(repvar != NULL);
   assert(negated != NULL);
   assert(var->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetBinvarRepresentative", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* get the active representative of the given variable */
   *repvar = var;
   *negated = FALSE;
   SCIP_CALL( SCIPvarGetProbvarBinary(repvar, negated) );

   /* negate the representative, if it corresponds to the negation of the given variable */
   if( *negated )
   {
      SCIP_CALL( SCIPgetNegatedVar(scip, *repvar, repvar) );
   }

   return SCIP_OKAY;
}

/** gets binary variables that are equal to the given binary variables, and which are either active, fixed, or
 *  multi-aggregated, or the negated variables of active, fixed, or multi-aggregated variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPgetBinvarRepresentatives(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of binary variables to get representatives for */
   SCIP_VAR**            vars,               /**< binary variables to get binary representatives for */
   SCIP_VAR**            repvars,            /**< array to store the binary representatives */
   SCIP_Bool*            negated             /**< array to store whether the negation of an active variable was returned */
   )
{
   int v;

   assert(scip != NULL);
   assert(vars != NULL || nvars == 0);
   assert(repvars != NULL || nvars == 0);
   assert(negated != NULL || nvars == 0);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetBinvarRepresentatives", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( nvars == 0 )
      return SCIP_OKAY;

   /* get the active representative of the given variable */
   BMScopyMemoryArray(repvars, vars, nvars);
   BMSclearMemoryArray(negated, nvars);
   SCIP_CALL( SCIPvarsGetProbvarBinary(&repvars, &negated, nvars) );

   /* negate the representatives, if they correspond to the negation of the given variables */
   for( v = nvars - 1; v >= 0; --v )
      if( negated[v] )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, repvars[v], &(repvars[v])) );
      }

   return SCIP_OKAY;
}

/** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later on
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPflattenVarAggregationGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   SCIP_CALL( SCIPcheckStage(scip, "SCIPflattenVarAggregationGraph", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarFlattenAggregationGraph(var, scip->mem->probmem, scip->set, scip->eventqueue) );

   return SCIP_OKAY;
}

/** Transforms a given linear sum of variables, that is a_1*x_1 + ... + a_n*x_n + c into a corresponding linear sum of
 *  active variables, that is b_1*y_1 + ... + b_m*y_m + d.
 *
 *  If the number of needed active variables is greater than the available slots in the variable array, nothing happens
 *  except that the required size is stored in the corresponding variable (requiredsize). Otherwise, the active variable
 *  representation is stored in the variable array, scalar array and constant.
 *
 *  The reason for this approach is that we cannot reallocate memory, since we do not know how the memory has been
 *  allocated (e.g., by a C++ 'new' or SCIP functions).
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The resulting linear sum is stored into the given variable array, scalar array, and constant. That means the
 *        given entries are overwritten.
 *
 *  @note That method can be used to convert a single variables into variable space of active variables. Therefore call
 *        the method with the linear sum 1.0*x + 0.0.
 */
SCIP_RETCODE SCIPgetProbvarLinearSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array x_1, ..., x_n in the linear sum which will be
                                              *   overwritten by the variable array y_1, ..., y_m in the linear sum
                                              *   w.r.t. active variables */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum which will be overwritten to the
                                              *   scalars b_1, ..., b_m in the linear sum of the active variables  */
   int*                  nvars,              /**< pointer to number of variables in the linear sum which will be
                                              *   overwritten by the number of variables in the linear sum corresponding
                                              *   to the active variables */
   int                   varssize,           /**< available slots in vars and scalars array which is needed to check if
                                              *   the array are large enough for the linear sum w.r.t. active
                                              *   variables */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c which
                                              *   will chnage to constant d in the linear sum b_1*y_1 + ... + b_m*y_m +
                                              *   d w.r.t. the active variables */
   int*                  requiredsize,       /**< pointer to store the required array size for the linear sum w.r.t. the
                                              *   active variables */
   SCIP_Bool             mergemultiples      /**< should multiple occurrences of a var be replaced by a single coeff? */
   )
{
   assert( scip != NULL );
   assert( nvars != NULL );
   assert( vars != NULL || *nvars == 0 );
   assert( scalars != NULL || *nvars == 0 );
   assert( constant != NULL );
   assert( requiredsize != NULL );
   assert( *nvars <= varssize );

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetProbvarLinearSum", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPvarGetActiveRepresentatives(scip->set, vars, scalars, nvars, varssize, constant, requiredsize, mergemultiples) );

   return SCIP_OKAY;
}

/** transforms given variable, scalar and constant to the corresponding active, fixed, or
 *  multi-aggregated variable, scalar and constant; if the variable resolves to a fixed variable,
 *  "scalar" will be 0.0 and the value of the sum will be stored in "constant"; a multi-aggregation
 *  with only one active variable (this can happen due to fixings after the multi-aggregation),
 *  is treated like an aggregation; if the multi-aggregation constant is infinite, "scalar" will be 0.0
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetProbvarSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to problem variable x in sum a*x + c */
   SCIP_Real*            scalar,             /**< pointer to scalar a in sum a*x + c */
   SCIP_Real*            constant            /**< pointer to constant c in sum a*x + c */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(scalar != NULL);
   assert(constant != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetProbvarSum", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPvarGetProbvarSum(var, scip->set, scalar, constant) );

   return SCIP_OKAY;
}

/** return for given variables all their active counterparts; all active variables will be pairwise different
 *  @note It does not hold that the first output variable is the active variable for the first input variable.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPgetActiveVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array with given variables and as output all active
					      *   variables, if enough slots exist
					      */
   int*                  nvars,              /**< number of given variables, and as output number of active variables,
					      *   if enough slots exist
					      */
   int                   varssize,           /**< available slots in vars array */
   int*                  requiredsize        /**< pointer to store the required array size for the active variables */
   )
{
   assert(scip != NULL);
   assert(nvars != NULL);
   assert(vars != NULL || *nvars == 0);
   assert(varssize >= *nvars);
   assert(requiredsize != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetActiveVars", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   SCIP_CALL( SCIPvarsGetActiveVars(scip->set, vars, nvars, varssize, requiredsize) );

   return SCIP_OKAY;
}

/** returns the reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 *
 *  @note The return value of this method should be used carefully if the dual feasibility check was explictely disabled.
 */
SCIP_Real SCIPgetVarRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( var->scip == scip );

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPgetVarRedcost(scip, var->data.original.transvar);

   case SCIP_VARSTATUS_COLUMN:
      return SCIPgetColRedcost(scip, SCIPvarGetCol(var));

   case SCIP_VARSTATUS_LOOSE:
      return SCIP_INVALID;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      return 0.0;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns the implied reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 *
 *  @note The return value of this method should be used carefully if the dual feasibility check was explictely disabled.
 */
SCIP_Real SCIPgetVarImplRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get reduced costs, should be a column in current node LP */
   SCIP_Bool             varfixing           /**< FALSE if for x == 0, TRUE for x == 1 */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   assert( var->scip == scip );

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPgetVarImplRedcost(scip, var->data.original.transvar, varfixing);

   case SCIP_VARSTATUS_COLUMN:
      return SCIPvarGetImplRedcost(var, scip->set, varfixing, scip->stat, scip->transprob, scip->lp);

   case SCIP_VARSTATUS_LOOSE:
      return SCIP_INVALID;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      return 0.0;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}


/** returns the Farkas coefficient of the variable in the current node's LP relaxation;
 *  the current node has to have an infeasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetVarFarkasCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(var->scip == scip);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPgetVarFarkasCoef(scip,var->data.original.transvar);

   case SCIP_VARSTATUS_COLUMN:
      return SCIPgetColFarkasCoef(scip,SCIPvarGetCol(var));

   case SCIP_VARSTATUS_LOOSE:
      return SCIP_INVALID;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      return 0.0;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return 0.0; /*lint !e527*/
   }
}

/** returns lower bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
SCIP_Real SCIPgetVarLbAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   SCIP_VARSTATUS varstatus;
   SCIP_BDCHGINFO* bdchginfo;
   assert(var != NULL);

   varstatus = SCIPvarGetStatus(var);

   /* get bounds of attached variables */
   switch( varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      return SCIPgetVarLbAtIndex(scip, var->data.original.transvar, bdchgidx, after);

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      if( bdchgidx == NULL )
         return SCIPvarGetLbLocal(var);
      else
      {
         bdchginfo = SCIPvarGetLbchgInfo(var, bdchgidx, after);
         if( bdchginfo != NULL )
            return SCIPbdchginfoGetNewbound(bdchginfo);
         else
            return var->glbdom.lb;
      }

   case SCIP_VARSTATUS_FIXED:
      return var->glbdom.lb;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_Real lb;

         lb = SCIPgetVarLbAtIndex(scip, var->data.aggregate.var, bdchgidx, after);

         /* a > 0 -> get lower bound of y */
         if( SCIPisInfinity(scip, -lb) )
            return -SCIPinfinity(scip);
         else if( SCIPisInfinity(scip, lb) )
            return SCIPinfinity(scip);
         else
            return var->data.aggregate.scalar * lb + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         SCIP_Real ub;

         ub = SCIPgetVarUbAtIndex(scip, var->data.aggregate.var, bdchgidx, after);

         /* a < 0 -> get upper bound of y */
         if( SCIPisInfinity(scip, -ub) )
            return SCIPinfinity(scip);
         else if( SCIPisInfinity(scip, ub) )
            return -SCIPinfinity(scip);
         else
            return var->data.aggregate.scalar * ub + var->data.aggregate.constant;
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         SCIPABORT();
         return SCIP_INVALID; /*lint !e527*/
      }

   case SCIP_VARSTATUS_MULTAGGR:
      /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
      if ( var->data.multaggr.nvars == 1 )
      {
         assert(var->data.multaggr.vars != NULL);
         assert(var->data.multaggr.scalars != NULL);
         assert(var->data.multaggr.vars[0] != NULL);

         if( var->data.multaggr.scalars[0] > 0.0 )
         {
            SCIP_Real lb;

            lb = SCIPgetVarLbAtIndex(scip, var->data.multaggr.vars[0], bdchgidx, after);

            /* a > 0 -> get lower bound of y */
            if( SCIPisInfinity(scip, -lb) )
               return -SCIPinfinity(scip);
            else if( SCIPisInfinity(scip, lb) )
               return SCIPinfinity(scip);
            else
               return var->data.multaggr.scalars[0] * lb + var->data.multaggr.constant;
         }
         else if( var->data.multaggr.scalars[0] < 0.0 )
         {
            SCIP_Real ub;

            ub = SCIPgetVarUbAtIndex(scip, var->data.multaggr.vars[0], bdchgidx, after);

            /* a < 0 -> get upper bound of y */
            if( SCIPisInfinity(scip, -ub) )
               return SCIPinfinity(scip);
            else if( SCIPisInfinity(scip, ub) )
               return -SCIPinfinity(scip);
            else
               return var->data.multaggr.scalars[0] * ub + var->data.multaggr.constant;
         }
         else
         {
            SCIPerrorMessage("scalar is zero in multi-aggregation\n");
            SCIPABORT();
            return SCIP_INVALID; /*lint !e527*/
         }
      }
      SCIPerrorMessage("cannot get the bounds of a multi-aggregated variable.\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPgetVarUbAtIndex(scip, var->negatedvar, bdchgidx, after);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** returns upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
SCIP_Real SCIPgetVarUbAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   SCIP_VARSTATUS varstatus;
   SCIP_BDCHGINFO* bdchginfo;
   assert(var != NULL);

   varstatus = SCIPvarGetStatus(var);

   /* get bounds of attached variables */
   switch( varstatus )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      assert(var->data.original.transvar != NULL);
      return SCIPgetVarUbAtIndex(scip, var->data.original.transvar, bdchgidx, after);

   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
      if( bdchgidx == NULL )
         return SCIPvarGetUbLocal(var);
      else
      {
         bdchginfo = SCIPvarGetUbchgInfo(var, bdchgidx, after);
         if( bdchginfo != NULL )
            return SCIPbdchginfoGetNewbound(bdchginfo);
         else
            return var->glbdom.ub;
      }

   case SCIP_VARSTATUS_FIXED:
      return var->glbdom.ub;

   case SCIP_VARSTATUS_AGGREGATED: /* x = a*y + c  ->  y = (x-c)/a */
      assert(var->data.aggregate.var != NULL);
      if( var->data.aggregate.scalar > 0.0 )
      {
         SCIP_Real ub;

         ub = SCIPgetVarUbAtIndex(scip, var->data.aggregate.var, bdchgidx, after);

         /* a > 0 -> get lower bound of y */
         if( SCIPisInfinity(scip, -ub) )
            return -SCIPinfinity(scip);
         else if( SCIPisInfinity(scip, ub) )
            return SCIPinfinity(scip);
         else
            return var->data.aggregate.scalar * ub + var->data.aggregate.constant;
      }
      else if( var->data.aggregate.scalar < 0.0 )
      {
         SCIP_Real lb;

         lb = SCIPgetVarLbAtIndex(scip, var->data.aggregate.var, bdchgidx, after);

         /* a < 0 -> get upper bound of y */
         if ( SCIPisInfinity(scip, -lb) )
            return SCIPinfinity(scip);
         else if ( SCIPisInfinity(scip, lb) )
            return -SCIPinfinity(scip);
         else
            return var->data.aggregate.scalar * lb + var->data.aggregate.constant;
      }
      else
      {
         SCIPerrorMessage("scalar is zero in aggregation\n");
         SCIPABORT();
         return SCIP_INVALID; /*lint !e527*/
      }

   case SCIP_VARSTATUS_MULTAGGR:
      /* handle multi-aggregated variables depending on one variable only (possibly caused by SCIPvarFlattenAggregationGraph()) */
      if ( var->data.multaggr.nvars == 1 )
      {
         assert(var->data.multaggr.vars != NULL);
         assert(var->data.multaggr.scalars != NULL);
         assert(var->data.multaggr.vars[0] != NULL);

         if( var->data.multaggr.scalars[0] > 0.0 )
         {
            SCIP_Real ub;

            ub = SCIPgetVarUbAtIndex(scip, var->data.multaggr.vars[0], bdchgidx, after);

            /* a > 0 -> get lower bound of y */
            if ( SCIPisInfinity(scip, -ub) )
               return -SCIPinfinity(scip);
            else if ( SCIPisInfinity(scip, ub) )
               return SCIPinfinity(scip);
            else
               return var->data.multaggr.scalars[0] * ub + var->data.multaggr.constant;
         }
         else if( var->data.multaggr.scalars[0] < 0.0 )
         {
            SCIP_Real lb;

            lb = SCIPgetVarLbAtIndex(scip, var->data.multaggr.vars[0], bdchgidx, after);

            /* a < 0 -> get upper bound of y */
            if ( SCIPisInfinity(scip, -lb) )
               return SCIPinfinity(scip);
            else if ( SCIPisInfinity(scip, lb) )
               return -SCIPinfinity(scip);
            else
               return var->data.multaggr.scalars[0] * lb + var->data.multaggr.constant;
         }
         else
         {
            SCIPerrorMessage("scalar is zero in multi-aggregation\n");
            SCIPABORT();
            return SCIP_INVALID; /*lint !e527*/
         }
      }
      SCIPerrorMessage("cannot get the bounds of a multiple aggregated variable.\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/

   case SCIP_VARSTATUS_NEGATED: /* x' = offset - x  ->  x = offset - x' */
      assert(var->negatedvar != NULL);
      assert(SCIPvarGetStatus(var->negatedvar) != SCIP_VARSTATUS_NEGATED);
      assert(var->negatedvar->negatedvar == var);
      return var->data.negate.constant - SCIPgetVarLbAtIndex(scip, var->negatedvar, bdchgidx, after);

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** returns lower or upper bound of variable directly before or after the bound change given by the bound change index
 *  was applied
 */
SCIP_Real SCIPgetVarBdAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
      return SCIPgetVarLbAtIndex(scip, var, bdchgidx, after);
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      return SCIPgetVarUbAtIndex(scip, var, bdchgidx, after);
   }
}

/** returns whether the binary variable was fixed at the time given by the bound change index */
SCIP_Bool SCIPgetVarWasFixedAtIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index representing time on path to current node */
   SCIP_Bool             after               /**< should the bound change with given index be included? */
   )
{
   assert(var != NULL);
   assert(SCIPvarIsBinary(var));

   /* check the current bounds first in order to decide at which bound change information we have to look
    * (which is expensive because we have to follow the aggregation tree to the active variable)
    */
   return ((SCIPvarGetLbLocal(var) > 0.5 && SCIPgetVarLbAtIndex(scip, var, bdchgidx, after) > 0.5)
      || (SCIPvarGetUbLocal(var) < 0.5 && SCIPgetVarUbAtIndex(scip, var, bdchgidx, after) < 0.5));
}

/** gets solution value for variable in current node
 *
 *  @return solution value for variable in current node
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetVarSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get solution value for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   assert( var->scip == scip );

   return SCIPvarGetSol(var, SCIPtreeHasCurrentNodeLP(scip->tree));
}

/** gets solution values of multiple variables in current node
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetVarSols(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarSols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPvarGetLPSol(vars[v]);
   }
   else
   {
      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPvarGetPseudoSol(vars[v]);
   }

   return SCIP_OKAY;
}

/** sets the solution value of all variables in the global relaxation solution to zero
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPclearRelaxSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< relaxator data structure */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPclearRelaxSolVals", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* update the responsible relax pointer */
   SCIPrelaxationSetSolRelax(scip->relaxation, relax);

   /* the relaxation solution is already cleared */
   if( SCIPrelaxationIsSolZero(scip->relaxation) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_CALL( SCIPvarSetRelaxSol(vars[v], scip->set, scip->relaxation, 0.0, FALSE) );
   }

   SCIPrelaxationSetSolObj(scip->relaxation, 0.0);
   SCIPrelaxationSetSolZero(scip->relaxation, TRUE);

   return SCIP_OKAY;
}

/** sets the value of the given variable in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  You can use SCIPclearRelaxSolVals() to set all values to zero, initially;
 *  after setting all solution values, you have to call SCIPmarkRelaxSolValid()
 *  to inform SCIP that the stored solution is valid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This method incrementally updates the objective value of the relaxation solution. If the whole solution
 *        should be updated, using SCIPsetRelaxSolVals() instead or calling SCIPclearRelaxSolVals() before setting
 *        the first value to reset the solution and the objective value to 0 may help the numerics.
 */
SCIP_RETCODE SCIPsetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxator data structure */
   SCIP_VAR*             var,                /**< variable to set value for */
   SCIP_Real             val                 /**< solution value of variable */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetRelaxSolVal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarSetRelaxSol(var, scip->set, scip->relaxation, val, TRUE) );

   if( val != 0.0 )
      SCIPrelaxationSetSolZero(scip->relaxation, FALSE);
   SCIPrelaxationSetSolValid(scip->relaxation, FALSE, FALSE);
   SCIPrelaxationSetSolRelax(scip->relaxation, relax);

   return SCIP_OKAY;
}

/** sets the values of the given variables in the global relaxation solution and informs SCIP about the validity
 *  and whether the solution can be enforced via linear cuts;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  the solution is automatically cleared, s.t. all other variables get value 0.0
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetRelaxSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxator data structure */
   int                   nvars,              /**< number of variables to set relaxation solution value for */
   SCIP_VAR**            vars,               /**< array with variables to set value for */
   SCIP_Real*            vals,               /**< array with solution values of variables */
   SCIP_Bool             includeslp          /**< does the relaxator contain all cuts in the LP? */
   )
{
   int v;

   assert(scip != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetRelaxSolVals", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPclearRelaxSolVals(scip, relax) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_CALL( SCIPvarSetRelaxSol(vars[v], scip->set, scip->relaxation, vals[v], TRUE) );
   }

   SCIPrelaxationSetSolZero(scip->relaxation, FALSE);
   SCIPrelaxationSetSolValid(scip->relaxation, TRUE, includeslp);
   SCIPrelaxationSetSolRelax(scip->relaxation, relax);

   return SCIP_OKAY;
}

/** sets the values of the variables in the global relaxation solution to the values in the given primal solution
 *  and informs SCIP about the validity and whether the solution can be enforced via linear cuts;
 *  the relaxation solution can be filled by the relaxation handlers and might be used by heuristics and for separation
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetRelaxSolValsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxator data structure */
   SCIP_SOL*             sol,                /**< primal relaxation solution */
   SCIP_Bool             includeslp          /**< does the relaxator contain all cuts in the LP? */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   int v;

   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetRelaxSolValsSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* alloc buffer array for solution values of the variables and get the values */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, vals) );

   SCIP_CALL( SCIPclearRelaxSolVals(scip, relax) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_CALL( SCIPvarSetRelaxSol(vars[v], scip->set, scip->relaxation, vals[v], FALSE) );
   }

   SCIPrelaxationSetSolObj(scip->relaxation, SCIPsolGetObj(sol, scip->set, scip->transprob, scip->origprob));

   SCIPrelaxationSetSolZero(scip->relaxation, FALSE);
   SCIPrelaxationSetSolValid(scip->relaxation, TRUE, includeslp);
   SCIPrelaxationSetSolRelax(scip->relaxation, relax);

   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

/** returns whether the relaxation solution is valid
 *
 *  @return TRUE, if the relaxation solution is valid; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Bool SCIPisRelaxSolValid(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisRelaxSolValid", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIPrelaxationIsSolValid(scip->relaxation);
}

/** informs SCIP that the relaxation solution is valid and whether the relaxation can be enforced through linear cuts
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPmarkRelaxSolValid(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxator data structure that set the current relaxation solution */
   SCIP_Bool             includeslp          /**< does the relaxator contain all cuts in the LP? */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPmarkRelaxSolValid", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPrelaxationSetSolValid(scip->relaxation, TRUE, includeslp);
   SCIPrelaxationSetSolRelax(scip->relaxation, relax);

   return SCIP_OKAY;
}

/** informs SCIP, that the relaxation solution is invalid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPmarkRelaxSolInvalid(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPmarkRelaxSolInvalid", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPrelaxationSetSolValid(scip->relaxation, FALSE, FALSE);

   return SCIP_OKAY;
}

/** gets the relaxation solution value of the given variable
 *
 *  @return the relaxation solution value of the given variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(var->scip == scip);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRelaxSolVal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("Relaxation Solution is not valid!\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }

   return SCIPvarGetRelaxSol(var, scip->set);
}

/** gets the relaxation solution objective value
 *
 *  @return the objective value of the relaxation solution
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetRelaxSolObj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRelaxSolObj", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("Relaxation Solution is not valid!\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }

   return SCIPrelaxationGetSolObj(scip->relaxation);
}

/** determine which branching direction should be evaluated first by strong branching
 *
 *  @return TRUE iff strong branching should first evaluate the down child
 *
 */
SCIP_Bool SCIPisStrongbranchDownFirst(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to determine the branching direction on */
   )
{
   switch( scip->set->branch_firstsbchild )
   {
      case 'u':
         return FALSE;
      case 'd':
         return TRUE;
      case 'a':
         return (SCIPvarGetNLocksDown(var) > SCIPvarGetNLocksUp(var));
      default:
         assert(scip->set->branch_firstsbchild == 'h');
         return (SCIPgetVarAvgCutoffs(scip, var, SCIP_BRANCHDIR_DOWNWARDS) > SCIPgetVarAvgCutoffs(scip, var, SCIP_BRANCHDIR_UPWARDS));
   }
}

/** start strong branching - call before any strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note if propagation is enabled, strong branching is not done directly on the LP, but probing nodes are created
 *        which allow to perform propagation but also creates some overhead
 */
SCIP_RETCODE SCIPstartStrongbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             enablepropagation   /**< should propagation be done before solving the strong branching LP? */
   )
{
   assert( scip != NULL );
   SCIP_CALL( SCIPcheckStage(scip, "SCIPstartStrongbranch", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(!SCIPinProbing(scip));

   SCIPdebugMsg(scip, "starting strong branching mode%s: lpcount=%" SCIP_LONGINT_FORMAT "\n", enablepropagation ? " with propagation" : "", scip->stat->lpcount - scip->stat->nsbdivinglps);

   /* start probing mode to allow propagation before solving the strong branching LPs; if no propagation should be done,
    * start the strong branching mode in the LP interface
    */
   if( enablepropagation )
   {
      if( SCIPtreeProbing(scip->tree) )
      {
         SCIPerrorMessage("cannot start strong branching with propagation while in probing mode\n");
         return SCIP_INVALIDCALL;
      }

      if( scip->lp != NULL && SCIPlpDiving(scip->lp) )
      {
         SCIPerrorMessage("cannot start strong branching with propagation while in diving mode\n");
         return SCIP_INVALIDCALL;
      }

      /* other then in SCIPstartProbing(), we do not disable collecting variable statistics during strong branching;
       * we cannot disable it, because the pseudo costs would not be updated, otherwise,
       * and reliability branching would end up doing strong branching all the time
       */
      SCIP_CALL( SCIPtreeStartProbing(scip->tree, scip->mem->probmem, scip->set, scip->lp, scip->relaxation, scip->transprob, TRUE) );

      /* inform the LP that the current probing mode is used for strong branching */
      SCIPlpStartStrongbranchProbing(scip->lp);
   }
   else
   {
      SCIP_CALL( SCIPlpStartStrongbranch(scip->lp) );
   }

   /* reset local strong branching info */
   scip->stat->lastsblpsolstats[0] = scip->stat->lastsblpsolstats[1] = SCIP_LPSOLSTAT_NOTSOLVED;

   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPendStrongbranch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   SCIP_CALL( SCIPcheckStage(scip, "SCIPendStrongbranch", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* depending on whether the strong branching mode was started with propagation enabled or not, we end the strong
    * branching probing mode or the LP strong branching mode
    */
   if( SCIPtreeProbing(scip->tree) )
   {
      SCIP_NODE* node;
      SCIP_DOMCHG* domchg;
      SCIP_VAR** boundchgvars;
      SCIP_Real* bounds;
      SCIP_BOUNDTYPE* boundtypes;
      int nboundchgs;
      int nbnds;
      int i;

      /* collect all bound changes deducted during probing, which were applied at the probing root and apply them to the
       * focusnode
       */
      node = SCIPgetCurrentNode(scip);
      assert(SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);
      assert(SCIPgetProbingDepth(scip) == 0);

      domchg = SCIPnodeGetDomchg(node);
      nboundchgs = SCIPdomchgGetNBoundchgs(domchg);

      SCIP_CALL( SCIPallocBufferArray(scip, &boundchgvars, nboundchgs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nboundchgs) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, nboundchgs) );

      for( i = 0, nbnds = 0; i < nboundchgs; ++i )
      {
         SCIP_BOUNDCHG* boundchg;

         boundchg = SCIPdomchgGetBoundchg(domchg, i);

         /* ignore redundant bound changes */
         if( SCIPboundchgIsRedundant(boundchg) )
            continue;

         boundchgvars[nbnds] = SCIPboundchgGetVar(boundchg);
         bounds[nbnds] = SCIPboundchgGetNewbound(boundchg);
         boundtypes[nbnds] = SCIPboundchgGetBoundtype(boundchg);
         ++nbnds;
      }

      SCIPdebugMsg(scip, "ending strong branching with probing: %d bound changes collected\n", nbnds);

      /* inform the LP that the probing mode is not used for strong branching anymore */
      SCIPlpEndStrongbranchProbing(scip->lp);

      /* switch back from probing to normal operation mode and restore variables and constraints to focus node */
      SCIP_CALL( SCIPtreeEndProbing(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat,
         scip->transprob, scip->origprob, scip->lp, scip->relaxation, scip->primal,
         scip->branchcand, scip->eventqueue, scip->eventfilter, scip->cliquetable) );

      /* apply the collected bound changes */
      for( i = 0; i < nbnds; ++i )
      {
         if( boundtypes[i] == SCIP_BOUNDTYPE_LOWER )
         {
            SCIPdebugMsg(scip, "apply probing lower bound change <%s> >= %.9g\n", SCIPvarGetName(boundchgvars[i]), bounds[i]);
            SCIP_CALL( SCIPchgVarLb(scip, boundchgvars[i], bounds[i]) );
         }
         else
         {
            SCIPdebugMsg(scip, "apply probing upper bound change <%s> <= %.9g\n", SCIPvarGetName(boundchgvars[i]), bounds[i]);
            SCIP_CALL( SCIPchgVarUb(scip, boundchgvars[i], bounds[i]) );
         }
      }

      SCIPfreeBufferArray(scip, &boundtypes);
      SCIPfreeBufferArray(scip, &bounds);
      SCIPfreeBufferArray(scip, &boundchgvars);
   }
   else
   {
      SCIPdebugMsg(scip, "ending strong branching\n");

      SCIP_CALL( SCIPlpEndStrongbranch(scip->lp) );
   }

   return SCIP_OKAY;
}

/** analyze the strong branching for the given variable; that includes conflict analysis for infeasible branches and
 *  storing of root reduced cost information
 */
static
SCIP_RETCODE analyzeStrongbranch(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to analyze */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict          /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   )
{
   SCIP_COL* col;
   SCIP_Bool downcutoff;
   SCIP_Bool upcutoff;

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   downcutoff = col->sbdownvalid && SCIPsetIsGE(scip->set, col->sbdown, scip->lp->cutoffbound);
   upcutoff = col->sbupvalid && SCIPsetIsGE(scip->set, col->sbup, scip->lp->cutoffbound);

   if( downinf != NULL )
      *downinf = downcutoff;
   if( upinf != NULL )
      *upinf = upcutoff;

   /* analyze infeasible strong branching sub problems:
    * because the strong branching's bound change is necessary for infeasibility, it cannot be undone;
    * therefore, infeasible strong branchings on non-binary variables will not produce a valid conflict constraint
    */
   if( scip->set->conf_enable && scip->set->conf_usesb && scip->set->nconflicthdlrs > 0
      && SCIPvarIsBinary(var) && SCIPtreeGetCurrentDepth(scip->tree) > 0 )
   {
      if( (downcutoff && SCIPsetFeasCeil(scip->set, col->primsol-1.0) >= col->lb - 0.5)
         || (upcutoff && SCIPsetFeasFloor(scip->set, col->primsol+1.0) <= col->ub + 0.5) )
      {
         assert(downconflict != NULL);
         assert(upconflict   != NULL);
         SCIP_CALL( SCIPconflictAnalyzeStrongbranch(scip->conflict, scip->conflictstore, scip->mem->probmem, scip->set, scip->stat,
               scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, col, downconflict, upconflict) );
      }
   }

   /* the strong branching results can be used to strengthen the root reduced cost information which is used for example
    * to propagate against the cutoff bound
    *
    * @note Ignore the results if the LP solution of the down (up) branch LP is smaller which should not happened by
    *       theory but can arise due to numerical issues.
    */
   if( SCIPtreeGetCurrentDepth(scip->tree) == 0 && SCIPvarIsBinary(var) && SCIPlpIsDualReliable(scip->lp) )
   {
      SCIP_Real lpobjval;

      assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

      lpobjval =  SCIPlpGetObjval(scip->lp, scip->set, scip->transprob);

      if( col->sbdownvalid && SCIPsetFeasCeil(scip->set, col->primsol-1.0) >= col->lb - 0.5 && lpobjval < col->sbdown )
         SCIPvarUpdateBestRootSol(var, scip->set, SCIPvarGetUbGlobal(var), -(col->sbdown - lpobjval), lpobjval);
      if( col->sbupvalid && SCIPsetFeasFloor(scip->set, col->primsol+1.0) <= col->ub + 0.5 && lpobjval < col->sbup )
         SCIPvarUpdateBestRootSol(var, scip->set, SCIPvarGetLbGlobal(var), col->sbup - lpobjval,  lpobjval);
   }

   return SCIP_OKAY;
}

/** gets strong branching information on column variable with fractional value
 *
 *  Before calling this method, the strong branching mode must have been activated by calling SCIPstartStrongbranch();
 *  after strong branching was done for all candidate variables, the strong branching mode must be ended by
 *  SCIPendStrongbranch(). Since this method does not apply domain propagation before strongbranching,
 *  propagation should not be enabled in the SCIPstartStrongbranch() call.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetVarStrongbranchFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Bool             idempotent,         /**< should scip's state remain the same after the call (statistics, column states...), or should it be updated ? */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL* col;
   SCIP_Real localdown;
   SCIP_Real localup;
   SCIP_Bool localdownvalid;
   SCIP_Bool localupvalid;

   assert(scip != NULL);
   assert(var != NULL);
   assert(lperror != NULL);
   assert(!SCIPtreeProbing(scip->tree)); /* we should not be in strong branching with propagation mode */
   assert(var->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarStrongbranchFrac", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( downvalid != NULL )
      *downvalid = FALSE;
   if( upvalid != NULL )
      *upvalid = FALSE;
   if( downinf != NULL )
      *downinf = FALSE;
   if( upinf != NULL )
      *upinf = FALSE;
   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
   {
      SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
      return SCIP_OKAY;
   }

   /* call strong branching for column with fractional value */
   SCIP_CALL( SCIPcolGetStrongbranch(col, FALSE, scip->set, scip->stat, scip->transprob, scip->lp, itlim, !idempotent, !idempotent,
         &localdown, &localup, &localdownvalid, &localupvalid, lperror) );

   /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
    * declare the sub nodes infeasible
    */
   if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
   {
      if( !idempotent ) /*lint !e774*/
      {
         SCIP_CALL( analyzeStrongbranch(scip, var, downinf, upinf, downconflict, upconflict) );
      }
      else
      {
         if( downinf != NULL )
            *downinf = localdownvalid && SCIPsetIsGE(scip->set, localdown, scip->lp->cutoffbound);
         if( upinf != NULL )
            *upinf = localupvalid && SCIPsetIsGE(scip->set, localup, scip->lp->cutoffbound);
      }
   }

   if( down != NULL )
      *down = localdown;
   if( up != NULL )
      *up = localup;
   if( downvalid != NULL )
      *downvalid = localdownvalid;
   if( upvalid != NULL )
      *upvalid = localupvalid;

   return SCIP_OKAY;
}

/** create, solve, and evaluate a single strong branching child (for strong branching with propagation) */
static
SCIP_RETCODE performStrongbranchWithPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   SCIP_Bool             down,               /**< do we regard the down child? */
   SCIP_Bool             firstchild,         /**< is this the first of the two strong branching children? */
   SCIP_Bool             propagate,          /**< should domain propagation be performed? */
   SCIP_Real             newbound,           /**< new bound to apply at the strong branching child */
   int                   itlim,              /**< iteration limit for strong branchings */
   int                   maxproprounds,      /**< maximum number of propagation rounds (-1: no limit, -2: parameter
                                              *   settings) */
   SCIP_Real*            value,              /**< stores dual bound for strong branching child */
   SCIP_Bool*            valid,              /**< stores whether the returned value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Longint*         ndomreductions,     /**< pointer to store the number of domain reductions found, or NULL */
   SCIP_Bool*            conflict,           /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible strong branching child, or NULL */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_VAR**            vars,               /**< active problem variables */
   int                   nvars,              /**< number of active problem variables */
   SCIP_Real*            newlbs,             /**< array to store valid lower bounds for all active variables, or NULL */
   SCIP_Real*            newubs,             /**< array to store valid upper bounds for all active variables, or NULL */
   SCIP_Bool*            foundsol,           /**< pointer to store whether a primal solution was found during strong branching */
   SCIP_Bool*            cutoff              /**< pointer to store whether the strong branching child is infeasible */
   )
{
   SCIP_Longint ndomreds;

   assert(value != NULL);
   assert(foundsol != NULL);
   assert(cutoff != NULL);
   assert(lperror != NULL);
   assert(valid != NULL ? !(*valid) : TRUE);

   *foundsol = FALSE;
   *cutoff = FALSE;
   *lperror = FALSE;

   /* check whether the strong branching child is already infeasible due to the bound change */
   if( down )
   {
      /* the down branch is infeasible due to the branching bound change; since this means that solval is not within the
       * bounds, this should only happen if previous strong branching calls on other variables detected bound changes which
       * are valid for and were already applied at the probing root
       */
      if( newbound < SCIPvarGetLbLocal(var) - 0.5 )
      {
         *value = SCIPinfinity(scip);

         if( valid != NULL )
            *valid = TRUE;

         /* bound changes are applied in SCIPendStrongbranch(), which can be seen as a conflict constraint */
         if( conflict != NULL )
            *conflict = TRUE;

         *cutoff = TRUE;

         return SCIP_OKAY;
      }
   }
   else
   {
      /* the up branch is infeasible due to the branching bound change; since this means that solval is not within the
       * bounds, this should only happen if previous strong branching calls on other variables detected bound changes which
       * are valid for and were already applied at the probing root
       */
      if( newbound > SCIPvarGetUbLocal(var) + 0.5 )
      {
         *value = SCIPinfinity(scip);

         if( valid != NULL )
            *valid = TRUE;

         /* bound changes are applied in SCIPendStrongbranch(), which can be seen as a conflict constraint */
         if( conflict != NULL )
            *conflict = TRUE;

         *cutoff = TRUE;

         return SCIP_OKAY;
      }
   }

   /* we need to ensure that we can create at least one new probing node without exceeding the maximal tree depth */
   if( SCIP_MAXTREEDEPTH > SCIPtreeGetProbingDepth(scip->tree) )
   {
      /* create a new probing node for the strong branching child and apply the new bound for the variable */
      SCIP_CALL( SCIPnewProbingNode(scip) );

      if( down )
      {
         assert(SCIPisGE(scip, newbound, SCIPvarGetLbLocal(var)));
         if( SCIPisLT(scip, newbound, SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPchgVarUbProbing(scip, var, newbound) );
         }
      }
      else
      {
         assert(SCIPisLE(scip, newbound, SCIPvarGetUbLocal(var)));
         if( SCIPisGT(scip, newbound, SCIPvarGetLbLocal(var)) )
         {
            SCIP_CALL( SCIPchgVarLbProbing(scip, var, newbound) );
         }
      }
   }
   else
   {
      if( valid != NULL )
         *valid = FALSE;

      *cutoff = FALSE;

      if( conflict != NULL )
         *conflict = FALSE;

      return SCIP_OKAY;
   }

   /* propagate domains at the probing node */
   if( propagate )
   {
      /* start time measuring */
      SCIPclockStart(scip->stat->strongpropclock, scip->set);

      ndomreds = 0;
      SCIP_CALL( SCIPpropagateProbing(scip, maxproprounds, cutoff, &ndomreds) );

      /* store number of domain reductions in strong branching */
      if( down )
         SCIPstatAdd(scip->stat, scip->set, nsbdowndomchgs, ndomreds);
      else
         SCIPstatAdd(scip->stat, scip->set, nsbupdomchgs, ndomreds);

      if( ndomreductions != NULL )
         *ndomreductions = ndomreds;

      /* stop time measuring */
      SCIPclockStop(scip->stat->strongpropclock, scip->set);

      if( *cutoff )
      {
         *value = SCIPinfinity(scip);

         if( valid != NULL )
            *valid = TRUE;

         SCIPdebugMsg(scip, "%s branch of var <%s> detected infeasible during propagation\n",
            down ? "down" : "up", SCIPvarGetName(var));
      }
   }

   /* if propagation did not already detect infeasibility, solve the probing LP */
   if( !(*cutoff) )
   {
      SCIP_CALL( SCIPsolveProbingLP(scip, itlim, lperror, cutoff) );
      assert(SCIPisLPRelax(scip));

      if( *cutoff )
      {
         assert(!(*lperror));

         *value = SCIPinfinity(scip);

         if( valid != NULL )
            *valid = TRUE;

         SCIPdebugMsg(scip, "%s branch of var <%s> detected infeasible in LP solving: status=%d\n",
            down ? "down" : "up", SCIPvarGetName(var), SCIPgetLPSolstat(scip));
      }
      else if( !(*lperror) )
      {
         /* save the lp solution status */
         scip->stat->lastsblpsolstats[down ? 0 : 1] = SCIPgetLPSolstat(scip);

         switch( SCIPgetLPSolstat(scip) )
         {
         case SCIP_LPSOLSTAT_OPTIMAL:
         {
            *value = SCIPgetLPObjval(scip);
            assert(SCIPisLT(scip, *value, SCIPgetCutoffbound(scip)));

            SCIPdebugMsg(scip, "probing LP solved to optimality, objective value: %16.9g\n", *value);

            if( valid != NULL )
               *valid = TRUE;

            /* check the strong branching LP solution for feasibility */
            SCIP_CALL( SCIPtryStrongbranchLPSol(scip, foundsol, cutoff) );
            break;
         }
         case SCIP_LPSOLSTAT_ITERLIMIT:
            ++scip->stat->nsbtimesiterlimhit;
            /*lint -fallthrough*/
         case SCIP_LPSOLSTAT_TIMELIMIT:
         {
            /* use LP value as estimate */
            SCIP_LPI* lpi;
            SCIP_Real objval;
            SCIP_Real looseobjval;

            SCIPdebugMsg(scip, "probing LP hit %s limit\n", SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_ITERLIMIT ? "iteration" : "time");

            /* we access the LPI directly, because when a time limit was hit, we cannot access objective value and dual
             * feasibility using the SCIPlp... methods; we should try to avoid direct calls to the LPI, but this is rather
             * uncritical here, because we are immediately after the SCIPsolveProbingLP() call, because we access the LPI
             * read-only, and we check SCIPlpiWasSolved() first
             */
            SCIP_CALL( SCIPgetLPI(scip, &lpi) );

            if( SCIPlpiWasSolved(lpi) )
            {
               SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
               looseobjval = SCIPlpGetLooseObjval(scip->lp, scip->set, scip->transprob);

               /* the infinity value in the LPI should not be smaller than SCIP's infinity value */
               assert(!SCIPlpiIsInfinity(lpi, objval) || SCIPisInfinity(scip, objval));

               /* we use SCIP's infinity value here because a value larger than this is counted as infeasible by SCIP */
               if( SCIPisInfinity(scip, objval) )
                  *value = SCIPinfinity(scip);
               else if( SCIPisInfinity(scip, -looseobjval) )
                  *value = -SCIPinfinity(scip);
               else
                  *value = objval + looseobjval;

               if( SCIPlpiIsDualFeasible(lpi) )
               {
                  if( valid != NULL )
                     *valid = TRUE;

                  if( SCIPisGE(scip, *value, SCIPgetCutoffbound(scip)) )
                     *cutoff = TRUE;
               }
            }
            break;
         }
         case SCIP_LPSOLSTAT_ERROR:
         case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
            *lperror = TRUE;
            break;
         case SCIP_LPSOLSTAT_NOTSOLVED: /* should only be the case for *cutoff = TRUE or *lperror = TRUE */
         case SCIP_LPSOLSTAT_OBJLIMIT: /* in this case, *cutoff should be TRUE and we should not get here */
         case SCIP_LPSOLSTAT_INFEASIBLE: /* in this case, *cutoff should be TRUE and we should not get here */
         default:
            SCIPerrorMessage("invalid LP solution status <%d>\n", SCIPgetLPSolstat(scip));
            return SCIP_INVALIDDATA;
         }  /*lint !e788*/
      }

      /* If columns are missing in the LP, the cutoff flag may be wrong. Therefore, we need to set it and the valid pointer
       * to false here.
       */
      if( (*cutoff) && !SCIPallColsInLP(scip) )
      {
         *cutoff = FALSE;
      }

#ifndef NDEBUG
      if( *lperror )
      {
         SCIPdebugMsg(scip, "error during strong branching probing LP solving: status=%d\n", SCIPgetLPSolstat(scip));
      }
#endif
   }

   /* if the subproblem was feasible, we store the local bounds of the variables after propagation and (possibly)
    * conflict analysis
    * @todo do this after propagation? should be able to get valid bounds more often, but they might be weaker
    */
   if( !(*cutoff) && newlbs != NULL)
   {
      int v;

      assert(newubs != NULL);

      /* initialize the newlbs and newubs to the current local bounds */
      if( firstchild )
      {
         for( v = 0; v < nvars; ++v )
         {
            newlbs[v] = SCIPvarGetLbLocal(vars[v]);
            newubs[v] = SCIPvarGetUbLocal(vars[v]);
         }
      }
      /* update newlbs and newubs: take the weaker of the already stored bounds and the current local bounds */
      else
      {
         for( v = 0; v < nvars; ++v )
         {
            SCIP_Real lb = SCIPvarGetLbLocal(vars[v]);
            SCIP_Real ub = SCIPvarGetUbLocal(vars[v]);

            newlbs[v] = MIN(newlbs[v], lb);
            newubs[v] = MAX(newubs[v], ub);
         }
      }
   }

   /* revert all changes at the probing node */
   SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

   return SCIP_OKAY;
}

/** gets strong branching information with previous domain propagation on column variable
 *
 *  Before calling this method, the strong branching mode must have been activated by calling SCIPstartStrongbranch();
 *  after strong branching was done for all candidate variables, the strong branching mode must be ended by
 *  SCIPendStrongbranch(). Since this method applies domain propagation before strongbranching, propagation has to be be
 *  enabled in the SCIPstartStrongbranch() call.
 *
 *  Before solving the strong branching LP, domain propagation can be performed. The number of propagation rounds
 *  can be specified by the parameter @p maxproprounds.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @warning When using this method, LP banching candidates and solution values must be copied beforehand, because
 *           they are updated w.r.t. the strong branching LP solution.
 */
SCIP_RETCODE SCIPgetVarStrongbranchWithPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   SCIP_Real             solval,             /**< value of the variable in the current LP solution */
   SCIP_Real             lpobjval,           /**< LP objective value of the current LP solution */
   int                   itlim,              /**< iteration limit for strong branchings */
   int                   maxproprounds,      /**< maximum number of propagation rounds (-1: no limit, -2: parameter
                                              *   settings) */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Longint*         ndomredsdown,       /**< pointer to store the number of domain reductions down, or NULL */
   SCIP_Longint*         ndomredsup,         /**< pointer to store the number of domain reductions up, or NULL */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   SCIP_Real*            newlbs,             /**< array to store valid lower bounds for all active variables, or NULL */
   SCIP_Real*            newubs              /**< array to store valid upper bounds for all active variables, or NULL */
   )
{
   SCIP_COL* col;
   SCIP_VAR** vars;
   SCIP_Longint oldniters;
   SCIP_Real newub;
   SCIP_Real newlb;
   SCIP_Bool propagate;
   SCIP_Bool cutoff;
   SCIP_Bool downchild;
   SCIP_Bool firstchild;
   SCIP_Bool foundsol;
   SCIP_Bool downvalidlocal;
   SCIP_Bool upvalidlocal;
   SCIP_Bool allcolsinlp;
   SCIP_Bool enabledconflict;
   int oldnconflicts;
   int nvars;

   assert(scip != NULL);
   assert(var != NULL);
   assert(SCIPvarIsIntegral(var));
   assert(down != NULL);
   assert(up != NULL);
   assert(lperror != NULL);
   assert((newlbs != NULL) == (newubs != NULL));
   assert(SCIPinProbing(scip));
   assert(var->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarStrongbranchWithPropagation", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether propagation should be performed */
   propagate = (maxproprounds != 0 && maxproprounds != -3);

   /* Check, if all existing columns are in LP.
    * If this is not the case, we may still return that the up and down dual bounds are valid, because the branching
    * rule should not apply them otherwise.
    * However, we must not set the downinf or upinf pointers to TRUE based on the dual bound, because we cannot
    * guarantee that this node can be cut off.
    */
   allcolsinlp = SCIPallColsInLP(scip);

   /* if maxproprounds is -2, change it to 0, which for the following calls means using the parameter settings */
   if( maxproprounds == -2 )
      maxproprounds = 0;

   *down = lpobjval;
   *up = lpobjval;
   if( downvalid != NULL )
      *downvalid = FALSE;
   if( upvalid != NULL )
      *upvalid = FALSE;
   if( downinf != NULL )
      *downinf = FALSE;
   if( upinf != NULL )
      *upinf = FALSE;
   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;
   if( ndomredsdown != NULL )
      *ndomredsdown = 0;
   if( ndomredsup != NULL )
      *ndomredsup = 0;

   *lperror = FALSE;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   scip->stat->lastsblpsolstats[0] = scip->stat->lastsblpsolstats[1] = SCIP_LPSOLSTAT_NOTSOLVED;

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
      return SCIP_OKAY;
   }

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
   {
      SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   newlb = SCIPfeasFloor(scip, solval + 1.0);
   newub = SCIPfeasCeil(scip, solval - 1.0);

   SCIPdebugMsg(scip, "strong branching on var <%s>: solval=%g, lb=%g, ub=%g\n", SCIPvarGetName(var), solval,
      SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

   /* the up branch is infeasible due to the branching bound change; since this means that solval is not within the
    * bounds, this should only happen if previous strong branching calls on other variables detected bound changes which
    * are valid for and were already applied at the probing root
    */
   if( newlb > SCIPvarGetUbLocal(var) + 0.5 )
   {
      *up = SCIPinfinity(scip);

      if( upinf != NULL )
         *upinf = TRUE;

      if( upvalid != NULL )
         *upvalid = TRUE;

      /* bound changes are applied in SCIPendStrongbranch(), which can be seen as a conflict constraint */
      if( upconflict != NULL )
         *upconflict = TRUE;

      SCIPcolSetStrongbranchData(col, scip->set, scip->stat, scip->lp, lpobjval, solval,
         *down, *up, FALSE, TRUE, 0LL, INT_MAX);

      /* we do not regard the down branch; its valid pointer stays set to FALSE */
      return SCIP_OKAY;
   }

   /* the down branch is infeasible due to the branching bound change; since this means that solval is not within the
    * bounds, this should only happen if previous strong branching calls on other variables detected bound changes which
    * are valid for and were already applied at the probing root
    */
   if( newub < SCIPvarGetLbLocal(var) - 0.5 )
   {
      *down = SCIPinfinity(scip);

      if( downinf != NULL )
         *downinf = TRUE;

      if( downvalid != NULL )
         *downvalid = TRUE;

      /* bound changes are applied in SCIPendStrongbranch(), which can be seen as a conflict constraint */
      if( downconflict != NULL )
         *downconflict = TRUE;

      SCIPcolSetStrongbranchData(col, scip->set, scip->stat, scip->lp, lpobjval, solval,
         *down, *up, TRUE, FALSE, 0LL, INT_MAX);

      /* we do not regard the up branch; its valid pointer stays set to FALSE */
      return SCIP_OKAY;
   }

   /* We now do strong branching by creating the two potential child nodes as probing nodes and solving them one after
    * the other. We will stop when the first child is detected infeasible, saving the effort we would need for the
    * second child. Since empirically, the up child tends to be infeasible more often, we do strongbranching first on
    * the up branch.
    */
   oldniters = scip->stat->nsbdivinglpiterations;
   firstchild = TRUE;
   cutoff = FALSE;

   /* switch conflict analysis according to usesb parameter */
   enabledconflict = scip->set->conf_enable;
   scip->set->conf_enable = (scip->set->conf_enable && scip->set->conf_usesb);

   /* @todo: decide the branch to look at first based on the cutoffs in previous calls? */
   downchild = SCIPisStrongbranchDownFirst(scip, var);

   downvalidlocal = FALSE;
   upvalidlocal = FALSE;

   do
   {
      oldnconflicts = SCIPconflictGetNConflicts(scip->conflict);

      if( downchild )
      {
         SCIP_CALL( performStrongbranchWithPropagation(scip, var, downchild, firstchild, propagate, newub, itlim, maxproprounds,
               down, &downvalidlocal, ndomredsdown, downconflict, lperror, vars, nvars, newlbs, newubs, &foundsol, &cutoff) );

         /* check whether a new solutions rendered the previous child infeasible */
         if( foundsol && !firstchild && allcolsinlp )
         {
            if( SCIPisGE(scip, *up, SCIPgetCutoffbound(scip)) )
            {
               if( upinf != NULL )
                  *upinf = TRUE;
            }
         }

         /* check for infeasibility */
         if( cutoff )
         {
            if( downinf != NULL )
               *downinf = TRUE;

            if( downconflict != NULL &&
               (SCIPvarGetLbLocal(var) > newub + 0.5 || SCIPconflictGetNConflicts(scip->conflict) > oldnconflicts) )
            {
               *downconflict = TRUE;
            }

            if( !scip->set->branch_forceall )
            {
               /* if this is the first call, we do not regard the up branch, its valid pointer is initially set to FALSE */
               break;
            }
         }
      }
      else
      {
         SCIP_CALL( performStrongbranchWithPropagation(scip, var, downchild, firstchild, propagate, newlb, itlim, maxproprounds,
               up, &upvalidlocal, ndomredsup, upconflict, lperror, vars, nvars, newlbs, newubs, &foundsol, &cutoff) );

         /* check whether a new solutions rendered the previous child infeasible */
         if( foundsol && !firstchild && allcolsinlp )
         {
            if( SCIPisGE(scip, *down, SCIPgetCutoffbound(scip)) )
            {
               if( downinf != NULL )
                  *downinf = TRUE;
            }
         }

         /* check for infeasibility */
         if( cutoff )
         {
            if( upinf != NULL )
               *upinf = TRUE;

            assert(upinf == NULL || (*upinf) == TRUE);

            if( upconflict != NULL &&
               (SCIPvarGetUbLocal(var) < newlb - 0.5 || SCIPconflictGetNConflicts(scip->conflict) > oldnconflicts) )
            {
               *upconflict = TRUE;
            }

            if( !scip->set->branch_forceall )
            {
               /* if this is the first call, we do not regard the down branch, its valid pointer is initially set to FALSE */
               break;
            }
         }
      }

      downchild = !downchild;
      firstchild = !firstchild;
   }
   while( !firstchild );

   /* set strong branching information in column */
   if( *lperror )
   {
      SCIPcolInvalidateStrongbranchData(col, scip->set, scip->stat, scip->lp);
   }
   else
   {
      SCIPcolSetStrongbranchData(col, scip->set, scip->stat, scip->lp, lpobjval, solval,
         *down, *up, downvalidlocal, upvalidlocal, scip->stat->nsbdivinglpiterations - oldniters, itlim);
   }

   if( downvalid != NULL )
      *downvalid = downvalidlocal;
   if( upvalid != NULL )
      *upvalid = upvalidlocal;

   scip->set->conf_enable = enabledconflict;

   return SCIP_OKAY;   /*lint !e438*/
}

/** gets strong branching information on column variable x with integral LP solution value (val); that is, the down branch
 *  is (val -1.0) and the up brach ins (val +1.0)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note If the integral LP solution value is the lower or upper bound of the variable, the corresponding branch will be
 *        marked as infeasible. That is, the valid pointer and the infeasible pointer are set to TRUE.
 */
SCIP_RETCODE SCIPgetVarStrongbranchInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Bool             idempotent,         /**< should scip's state remain the same after the call (statistics, column states...), or should it be updated ? */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL* col;
   SCIP_Real localdown;
   SCIP_Real localup;
   SCIP_Bool localdownvalid;
   SCIP_Bool localupvalid;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarStrongbranchInt", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(lperror != NULL);
   assert(var->scip == scip);

   if( downvalid != NULL )
      *downvalid = FALSE;
   if( upvalid != NULL )
      *upvalid = FALSE;
   if( downinf != NULL )
      *downinf = FALSE;
   if( upinf != NULL )
      *upinf = FALSE;
   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
   {
      SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
      return SCIP_OKAY;
   }

   /* call strong branching for column */
   SCIP_CALL( SCIPcolGetStrongbranch(col, TRUE, scip->set, scip->stat, scip->transprob, scip->lp, itlim, !idempotent, !idempotent,
         &localdown, &localup, &localdownvalid, &localupvalid, lperror) );

   /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
    * declare the sub nodes infeasible
    */
   if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
   {
      if( !idempotent ) /*lint !e774*/
      {
         SCIP_CALL( analyzeStrongbranch(scip, var, downinf, upinf, downconflict, upconflict) );
      }
      else
      {
         if( downinf != NULL )
            *downinf = localdownvalid && SCIPsetIsGE(scip->set, localdown, scip->lp->cutoffbound);
         if( upinf != NULL )
            *upinf = localupvalid && SCIPsetIsGE(scip->set, localup, scip->lp->cutoffbound);
      }
   }

   if( down != NULL )
      *down = localdown;
   if( up != NULL )
      *up = localup;
   if( downvalid != NULL )
      *downvalid = localdownvalid;
   if( upvalid != NULL )
      *upvalid = localupvalid;

   return SCIP_OKAY;
}

/** gets strong branching information on column variables with fractional values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetVarsStrongbranchesFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables to get strong branching values for */
   int                   nvars,              /**< number of variables */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching variables down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching variables up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< array to store whether the downward branches are infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< array to store whether the upward branches are infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< array to store whether conflict constraints were created for
                                              *   infeasible downward branches, or NULL */
   SCIP_Bool*            upconflict,         /**< array to store whether conflict constraints were created for
                                              *   infeasible upward branches, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL** cols;
   int j;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarsStrongbranchesFrac", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( lperror != NULL );
   assert( vars != NULL );

   /* set up data */
   cols = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &cols, nvars) );
   assert(cols != NULL);
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      SCIP_COL* col;

      if( downvalid != NULL )
         downvalid[j] = FALSE;
      if( upvalid != NULL )
         upvalid[j] = FALSE;
      if( downinf != NULL )
         downinf[j] = FALSE;
      if( upinf != NULL )
         upinf[j] = FALSE;
      if( downconflict != NULL )
         downconflict[j] = FALSE;
      if( upconflict != NULL )
         upconflict[j] = FALSE;

      var = vars[j];
      assert( var != NULL );
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      cols[j] = col;

      if( !SCIPcolIsInLP(col) )
      {
         SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
   }
   else
   {
      /* call strong branching for columns with fractional value */
      SCIP_CALL( SCIPcolGetStrongbranches(cols, nvars, FALSE, scip->set, scip->stat, scip->transprob, scip->lp, itlim,
            down, up, downvalid, upvalid, lperror) );

      /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
       * declare the sub nodes infeasible
       */
      if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
      {
         for( j = 0; j < nvars; ++j )
         {
            SCIP_CALL( analyzeStrongbranch(scip, vars[j], (downinf != NULL) ? (&(downinf[j])) : NULL,
                  (upinf != NULL) ? (&(upinf[j])) : NULL, (downconflict != NULL) ? (&(downconflict[j])) : NULL,
                  (upconflict != NULL) ? (&(upconflict[j])) : NULL) );
         }
      }
   }
   SCIPfreeBufferArray(scip, &cols);

   return SCIP_OKAY;
}

/** gets strong branching information on column variables with integral values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetVarsStrongbranchesInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables to get strong branching values for */
   int                   nvars,              /**< number of variables */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching variables down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching variables up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< array to store whether the downward branches are infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< array to store whether the upward branches are infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< array to store whether conflict constraints were created for
                                              *   infeasible downward branches, or NULL */
   SCIP_Bool*            upconflict,         /**< array to store whether conflict constraints were created for
                                              *   infeasible upward branches, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL** cols;
   int j;

   assert(lperror != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarsStrongbranchesInt", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( vars != NULL );

   /* set up data */
   cols = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &cols, nvars) );
   assert(cols != NULL);
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      SCIP_COL* col;

      if( downvalid != NULL )
         downvalid[j] = FALSE;
      if( upvalid != NULL )
         upvalid[j] = FALSE;
      if( downinf != NULL )
         downinf[j] = FALSE;
      if( upinf != NULL )
         upinf[j] = FALSE;
      if( downconflict != NULL )
         downconflict[j] = FALSE;
      if( upconflict != NULL )
         upconflict[j] = FALSE;

      var = vars[j];
      assert( var != NULL );
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      cols[j] = col;

      if( !SCIPcolIsInLP(col) )
      {
         SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
   }
   else
   {
      /* call strong branching for columns */
      SCIP_CALL( SCIPcolGetStrongbranches(cols, nvars, TRUE, scip->set, scip->stat, scip->transprob, scip->lp, itlim,
            down, up, downvalid, upvalid, lperror) );

      /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
       * declare the sub nodes infeasible
       */
      if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
      {
         for( j = 0; j < nvars; ++j )
         {
            SCIP_CALL( analyzeStrongbranch(scip, vars[j], (downinf != NULL) ? (&(downinf[j])) : NULL,
                  (upinf != NULL) ? (&(upinf[j])) : NULL, (downconflict != NULL) ? (&(downconflict[j])) : NULL,
                  (upconflict != NULL) ? (&(upconflict[j])) : NULL) );
         }
      }
   }
   SCIPfreeBufferArray(scip, &cols);

   return SCIP_OKAY;
}

/** get LP solution status of last strong branching call (currently only works for strong branching with propagation) */
SCIP_LPSOLSTAT SCIPgetLastStrongbranchLPSolStat(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        branchdir           /**< branching direction for which LP solution status is requested */
   )
{
   assert(NULL != scip);
   assert(branchdir == SCIP_BRANCHDIR_DOWNWARDS || branchdir == SCIP_BRANCHDIR_UPWARDS);

   return scip->stat->lastsblpsolstats[branchdir == SCIP_BRANCHDIR_DOWNWARDS ? 0 : 1];
}

/** gets strong branching information on COLUMN variable of the last SCIPgetVarStrongbranch() call;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given variable;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPgetVarStrongbranchLast(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get last strong branching values for */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Real*            solval,             /**< stores LP solution value of variable at the last strong branching call, or NULL */
   SCIP_Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarStrongbranchLast", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable\n");
      return SCIP_INVALIDDATA;
   }

   SCIPcolGetStrongbranchLast(SCIPvarGetCol(var), down, up, downvalid, upvalid, solval, lpobjval);

   return SCIP_OKAY;
}

/** sets strong branching information for a column variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetVarStrongbranchData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to set last strong branching values for */
   SCIP_Real             lpobjval,           /**< objective value of the current LP */
   SCIP_Real             primsol,            /**< primal solution value of the column in the current LP */
   SCIP_Real             down,               /**< dual bound after branching column down */
   SCIP_Real             up,                 /**< dual bound after branching column up */
   SCIP_Bool             downvalid,          /**< is the returned down value a valid dual bound? */
   SCIP_Bool             upvalid,            /**< is the returned up value a valid dual bound? */
   SCIP_Longint          iter,               /**< total number of strong branching iterations */
   int                   itlim               /**< iteration limit applied to the strong branching call */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetVarStrongbranchData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot set strong branching information on non-COLUMN variable\n");
      return SCIP_INVALIDDATA;
   }

   SCIPcolSetStrongbranchData(SCIPvarGetCol(var), scip->set, scip->stat, scip->lp, lpobjval, primsol,
      down, up, downvalid, upvalid, iter, itlim);

   return SCIP_OKAY;
}

/** rounds the current solution and tries it afterwards; if feasible, adds it to storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPtryStrongbranchLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            foundsol,           /**< stores whether solution was feasible and good enough to keep */
   SCIP_Bool*            cutoff              /**< stores whether solution was cutoff due to exceeding the cutoffbound */
   )
{
   assert(scip != NULL);
   assert(foundsol != NULL);
   assert(cutoff != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtryStrongbranchLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->set->branch_checksbsol )
   {
      SCIP_SOL* sol;
      SCIP_Bool rounded = TRUE;
      SCIP_Real value = SCIPgetLPObjval(scip);
      SCIP_Longint oldnbestsolsfound = scip->primal->nbestsolsfound;

      /* start clock for strong branching solutions */
      SCIPclockStart(scip->stat->sbsoltime, scip->set);

      SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL) );
      SCIPsolSetStrongbranching(sol);

      /* try to round the strong branching solution */
      if( scip->set->branch_roundsbsol )
      {
         SCIP_CALL( SCIProundSol(scip, sol, &rounded) );
      }

      /* check the solution for feasibility if rounding worked well (or was not tried) */
      if( rounded )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, FALSE, TRUE, FALSE, foundsol) );
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }

      if( *foundsol )
      {
         SCIPdebugMsg(scip, "found new solution in strong branching\n");

         scip->stat->nsbsolsfound++;

         if( scip->primal->nbestsolsfound != oldnbestsolsfound )
         {
            scip->stat->nsbbestsolsfound++;
         }

         if( SCIPisGE(scip, value, SCIPgetCutoffbound(scip)) )
            *cutoff = TRUE;
      }

      /* stop clock for strong branching solutions */
      SCIPclockStop(scip->stat->sbsoltime, scip->set);
   }
   return SCIP_OKAY;
}


/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given variable, or -1 if strong branching was never applied to the variable in current run
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetVarStrongbranchNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarStrongbranchNode", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return -1;

   return SCIPcolGetStrongbranchNode(SCIPvarGetCol(var));
}

/** if strong branching was already applied on the variable at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this variable was applied;
 *  if strong branching was not yet applied on the variable at the current node, returns INT_MAX
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetVarStrongbranchLPAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get strong branching LP age for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarStrongbranchLPAge", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return SCIP_LONGINT_MAX;

   return SCIPcolGetStrongbranchLPAge(SCIPvarGetCol(var), scip->stat);
}

/** gets number of times, strong branching was applied in current run on the given variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetVarNStrongbranchs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarNStrongbranchs", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return 0;

   return SCIPcolGetNStrongbranchs(SCIPvarGetCol(var));
}

/** adds given values to lock numbers of type @p locktype of variable for rounding
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPaddVarLocksType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_LOCKTYPE         locktype,           /**< type of the variable locks */
   int                   nlocksdown,         /**< modification in number of rounding down locks */
   int                   nlocksup            /**< modification in number of rounding up locks */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarLocksType", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );

   assert( var->scip == scip );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      SCIP_CALL( SCIPvarAddLocks(var, scip->mem->probmem, scip->set, scip->eventqueue, locktype, nlocksdown, nlocksup) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** adds given values to lock numbers of variable for rounding
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note This method will always add variable locks of type model
 *
 *  @note It is recommented to use SCIPaddVarLocksType()
 */
SCIP_RETCODE SCIPaddVarLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   nlocksdown,         /**< modification in number of rounding down locks */
   int                   nlocksup            /**< modification in number of rounding up locks */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarLocks", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPaddVarLocksType(scip, var, SCIP_LOCKTYPE_MODEL, nlocksdown, nlocksup) );

   return SCIP_OKAY;
}

/** add locks of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be locked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be locked in upwards direction? */
   )
{
   int nlocksdown[NLOCKTYPES];
   int nlocksup[NLOCKTYPES];
   int i;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPlockVarCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );

   assert( var->scip == scip );

   for( i = 0; i < NLOCKTYPES; i++ )
   {
      nlocksdown[i] = 0;
      nlocksup[i] = 0;

      if( SCIPconsIsLockedTypePos(cons, (SCIP_LOCKTYPE) i) )
      {
         if( lockdown )
            ++nlocksdown[i];
         if( lockup )
            ++nlocksup[i];
      }
      if( SCIPconsIsLockedTypeNeg(cons, (SCIP_LOCKTYPE) i) )
      {
         if( lockdown )
            ++nlocksup[i];
         if( lockup )
            ++nlocksdown[i];
      }
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      for( i = 0; i < NLOCKTYPES; i++ )
      {
         if( nlocksdown[i] == 0 && nlocksup[i] == 0 )
            continue;

         SCIP_CALL( SCIPvarAddLocks(var, scip->mem->probmem, scip->set, scip->eventqueue, (SCIP_LOCKTYPE) i, nlocksdown[i], nlocksup[i]) );
      }
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** remove locks of type @p locktype of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_RETCODE SCIPunlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be unlocked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be unlocked in upwards direction? */
   )
{
   int nlocksdown[NLOCKTYPES];
   int nlocksup[NLOCKTYPES];
   int i;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPunlockVarCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE) );

   assert( var->scip == scip );

   for( i = 0; i < NLOCKTYPES; i++ )
   {
      nlocksdown[i] = 0;
      nlocksup[i] = 0;

      if( SCIPconsIsLockedTypePos(cons, (SCIP_LOCKTYPE) i) )
      {
         if( lockdown )
            ++nlocksdown[i];
         if( lockup )
            ++nlocksup[i];
      }
      if( SCIPconsIsLockedTypeNeg(cons, (SCIP_LOCKTYPE) i) )
      {
         if( lockdown )
            ++nlocksup[i];
         if( lockup )
            ++nlocksdown[i];
      }
   }
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
      for( i = 0; i < NLOCKTYPES; i++ )
      {
         if( nlocksdown[i] == 0 && nlocksup[i] == 0 )
            continue;

         SCIP_CALL( SCIPvarAddLocks(var, scip->mem->probmem, scip->set, scip->eventqueue, (SCIP_LOCKTYPE)  i, -nlocksdown[i], -nlocksup[i]) );
      }
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** changes variable's objective value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_RETCODE SCIPchgVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarObj", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   /* forbid infinite objective values */
   if( SCIPisInfinity(scip, REALABS(newobj)) )
   {
      SCIPerrorMessage("invalid objective value: objective value is infinite\n");
      return SCIP_INVALIDDATA;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgObj(var, scip->mem->probmem, scip->set, scip->origprob, scip->primal, scip->lp, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
      SCIP_CALL( SCIPvarChgObj(var, scip->mem->probmem, scip->set,  scip->transprob, scip->primal, scip->lp, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** adds value to variable's objective value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_RETCODE SCIPaddVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             addobj              /**< additional objective value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarObj", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarAddObj(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob, scip->primal,
            scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue, addobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
      SCIP_CALL( SCIPvarAddObj(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob, scip->primal,
            scip->tree, scip->reopt, scip->lp, scip->eventfilter, scip->eventqueue, addobj) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 *
 *  @return adjusted lower bound for the given variable; the bound of the variable is not changed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Real SCIPadjustedVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             lb                  /**< lower bound value to adjust */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPadjustedVarLb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &lb);

   return lb;
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) upper bound value;
 *  does not change the bounds of the variable
 *
 *  @return adjusted upper bound for the given variable; the bound of the variable is not changed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
SCIP_Real SCIPadjustedVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             ub                  /**< upper bound value to adjust */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPadjustedVarUb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIPvarAdjustUb(var, scip->set, &ub);

   return ub;
}

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLb", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* ignore tightenings of lower bounds to +infinity during solving process */
   if( SCIPisInfinity(scip, newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var));
#endif
      return SCIP_OKAY;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVED:
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable,
               var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, newbound,
            SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUb", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* ignore tightenings of upper bounds to -infinity during solving process */
   if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
         SCIPvarGetUbLocal(var));
#endif
      return SCIP_OKAY;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVED:
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
               scip->cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarLbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLbNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( node == NULL )
   {
      SCIP_CALL( SCIPchgVarLb(scip, var, newbound) );
   }
   else
   {
      SCIPvarAdjustLb(var, scip->set, &newbound);

      /* ignore tightenings of lower bounds to +infinity during solving process */
      if( SCIPisInfinity(scip, newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
#ifndef NDEBUG
         SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
            SCIPvarGetLbLocal(var));
#endif
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPnodeAddBoundchg(node, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
            SCIP_BOUNDTYPE_LOWER, FALSE) );
   }

   return SCIP_OKAY;
}

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarUbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUbNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( node == NULL )
   {
      SCIP_CALL( SCIPchgVarUb(scip, var, newbound) );
   }
   else
   {
      SCIPvarAdjustUb(var, scip->set, &newbound);

      /* ignore tightenings of upper bounds to -infinity during solving process */
      if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
#ifndef NDEBUG
         SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
            SCIPvarGetUbLocal(var));
#endif
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPnodeAddBoundchg(node, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
            SCIP_BOUNDTYPE_UPPER, FALSE) );
   }

   return SCIP_OKAY;
}

/** changes global lower bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLbGlobal", FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* ignore tightenings of lower bounds to +infinity during solving process */
   if( SCIPisInfinity(scip, newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var));
#endif
      return SCIP_OKAY;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
            SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** changes global upper bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUbGlobal", FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* ignore tightenings of upper bounds to -infinity during solving process */
   if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
         SCIPvarGetUbLocal(var));
#endif
      return SCIP_OKAY;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
            SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** changes lazy lower bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  Lazy bounds are bounds that are already enforced by constraints and the objective function.
 *  Setting a lazy lower bound has the consequence that for variables which lower bound equals the lazy lower bound,
 *  the lower bound does not need to be passed on to the LP solver.
 *  This is especially useful in a column generation (branch-and-price) setting.
 *
 *  @attention If the variable has a global lower bound below lazylb, then the global lower bound is tightened to
 *     lazylb by a call to SCIPchgVarLbGlobal().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarLbLazy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lazylb              /**< the lazy lower bound to be set */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLbLazy", FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPisGT(scip, lazylb, SCIPvarGetLbGlobal(var)) )
   {
      SCIP_CALL( SCIPchgVarLbGlobal(scip, var, lazylb) );
   }

   SCIP_CALL( SCIPvarChgLbLazy(var, scip->set, lazylb) );

   return SCIP_OKAY;
}

/** changes lazy upper bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  Lazy bounds are bounds that are already enforced by constraints and the objective function.
 *  Setting a lazy upper bound has the consequence that for variables which upper bound equals the lazy upper bound,
 *  the upper bound does not need to be passed on to the LP solver.
 *  This is especially useful in a column generation (branch-and-price) setting.
 *
 *  @attention If the variable has a global upper bound above lazyub, then the global upper bound is tightened to
 *     lazyub by a call to SCIPchgVarUbGlobal().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarUbLazy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lazyub              /**< the lazy lower bound to be set */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUbLazy", FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPisLT(scip, lazyub, SCIPvarGetUbGlobal(var)) )
   {
      SCIP_CALL( SCIPchgVarUbGlobal(scip, var, lazyub) );
   }

   SCIP_CALL( SCIPvarChgUbLazy(var, scip->set, lazyub) );

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtightenVarLb", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   /** @todo if needed provide pending local/global bound changes that will be flushed after leaving diving mode (as in struct_tree.h) */
   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM || !SCIPinDive(scip));

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* ignore tightenings of lower bounds to +infinity during solving process */
   if( SCIPisInfinity(scip, newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPcomputeVarLbLocal(scip, var);
   ub = SCIPcomputeVarUbLocal(scip, var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( (force && SCIPsetIsLE(scip->set, newbound, lb)) || (!force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;
   case SCIP_STAGE_TRANSFORMED:
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;
   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable,
            var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* check whether the lower bound improved */
   if( tightened != NULL && lb < SCIPcomputeVarLbLocal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPtightenVarUb", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /** @todo if needed provide pending local/global bound changes that will be flushed after leaving diving mode (as in struct_tree.h) */
   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM || !SCIPinDive(scip));

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* ignore tightenings of upper bounds to -infinity during solving process */
   if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
         SCIPvarGetUbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPcomputeVarLbLocal(scip, var);
   ub = SCIPcomputeVarUbLocal(scip, var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( (force && SCIPsetIsGE(scip->set, newbound, ub)) || (!force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;
   case SCIP_STAGE_TRANSFORMED:
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;
   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* check whether the upper bound improved */
   if( tightened != NULL && ub > SCIPcomputeVarUbLocal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** fixes variable in preprocessing or in the current node, if the new bound is tighter (w.r.t. bound strengthening
 *  epsilon) than the current bound; if possible, adjusts bound to integral value; the given inference constraint is
 *  stored, such that the conflict analysis is able to find out the reason for the deduction of the bound change
 *
 *  @note In presolving stage when not in probing mode the variable will be fixed directly, otherwise this method
 *        changes first the lowerbound by calling SCIPinferVarLbCons and second the upperbound by calling
 *        SCIPinferVarUbCons
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarFixCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             fixedval,           /**< new value for fixation */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferVarFixCons", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( tightened != NULL )
      *tightened = FALSE;

   /* in presolving case we take the shortcut to directly fix the variables */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && SCIPtreeGetCurrentDepth(scip->tree) == 0 )
   {
      SCIP_Bool fixed;

      SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter,
            scip->eventqueue, scip->cliquetable, fixedval, infeasible, &fixed) );

      if( tightened != NULL )
	 *tightened = fixed;
   }
   /* otherwise we use the lb and ub methods */
   else
   {
      SCIP_Bool lbtightened;

      SCIP_CALL( SCIPinferVarLbCons(scip, var, fixedval, infercons, inferinfo, force, infeasible, &lbtightened) );

      if( ! (*infeasible) )
      {
	 SCIP_CALL( SCIPinferVarUbCons(scip, var, fixedval, infercons, inferinfo, force, infeasible, tightened) );

	 if( tightened != NULL )
	    *tightened |= lbtightened;
      }
   }

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarLbCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferVarLbCons", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* ignore tightenings of lower bounds to +infinity during solving process */
   if( SCIPisInfinity(scip, newbound)  && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( (force && SCIPsetIsLE(scip->set, newbound, lb)) || (!force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, infercons, NULL, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* check whether the lower bound improved */
   if( tightened != NULL && lb < SCIPcomputeVarLbLocal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarUbCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferVarUbCons", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* ignore tightenings of upper bounds to -infinity during solving process */
   if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
         SCIPvarGetUbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( (force && SCIPsetIsGE(scip->set, newbound, ub)) || (!force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, infercons, NULL, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* check whether the upper bound improved */
   if( tightened != NULL && ub > SCIPcomputeVarUbLocal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPinferBinvarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to fix */
   SCIP_Bool             fixedval,           /**< value to fix binary variable to */
   SCIP_CONS*            infercons,          /**< constraint that deduced the fixing */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(SCIPvarIsBinary(var));
   assert(fixedval == TRUE || fixedval == FALSE);
   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferBinvarCons", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsEQ(scip->set, lb, 0.0) || SCIPsetIsEQ(scip->set, lb, 1.0));
   assert(SCIPsetIsEQ(scip->set, ub, 0.0) || SCIPsetIsEQ(scip->set, ub, 1.0));
   assert(SCIPsetIsLE(scip->set, lb, ub));

   /* check, if variable is already fixed */
   if( (lb > 0.5) || (ub < 0.5) )
   {
      *infeasible = (fixedval == (lb < 0.5));

      return SCIP_OKAY;
   }

   /* apply the fixing */
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      {
         SCIP_Bool fixed;

         SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue,
               scip->cliquetable, (SCIP_Real)fixedval, infeasible, &fixed) );
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
               scip->cliquetable, var, 1.0, SCIP_BOUNDTYPE_LOWER, infercons, NULL, inferinfo, FALSE) );
      }
      else
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
               scip->cliquetable, var, 0.0, SCIP_BOUNDTYPE_UPPER, infercons, NULL, inferinfo, FALSE) );
      }
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** fixes variable in preprocessing or in the current node, if the new bound is tighter (w.r.t. bound strengthening
 *  epsilon) than the current bound; if possible, adjusts bound to integral value; the given inference constraint is
 *  stored, such that the conflict analysis is able to find out the reason for the deduction of the bound change
 *
 *  @note In presolving stage when not in probing mode the variable will be fixed directly, otherwise this method
 *        changes first the lowerbound by calling SCIPinferVarLbProp and second the upperbound by calling
 *        SCIPinferVarUbProp
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarFixProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             fixedval,           /**< new value for fixation */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferVarFixProp", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( tightened != NULL )
      *tightened = FALSE;

   /* in presolving case we take the shortcut to directly fix the variables */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && SCIPtreeGetCurrentDepth(scip->tree) == 0 )
   {
      SCIP_Bool fixed;

      SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue,
            scip->cliquetable, fixedval, infeasible, &fixed) );

      if( tightened != NULL )
	 *tightened = fixed;
   }
   /* otherwise we use the lb and ub methods */
   else
   {
      SCIP_Bool lbtightened;

      SCIP_CALL( SCIPinferVarLbProp(scip, var, fixedval, inferprop, inferinfo, force, infeasible, &lbtightened) );

      if( ! (*infeasible) )
      {
	 SCIP_CALL( SCIPinferVarUbProp(scip, var, fixedval, inferprop, inferinfo, force, infeasible, tightened) );

	 if( tightened != NULL )
	    *tightened |= lbtightened;
      }
   }

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarLbProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferVarLbProp", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* ignore tightenings of lower bounds to +infinity during solving process */
   if( SCIPisInfinity(scip, newbound)  && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( (!force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub))
      || SCIPsetIsLE(scip->set, newbound, lb) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, newbound, SCIP_BOUNDTYPE_LOWER, NULL, inferprop, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* check whether the lower bound improved */
   if( tightened != NULL && lb < SCIPcomputeVarLbLocal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarUbProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferVarUbProp", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* ignore tightenings of upper bounds to -infinity during solving process */
   if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
         SCIPvarGetUbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( (!force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub))
      || SCIPsetIsGE(scip->set, newbound, ub) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue,
            scip->cliquetable, var, newbound, SCIP_BOUNDTYPE_UPPER, NULL, inferprop, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* check whether the upper bound improved */
   if( tightened != NULL && ub > SCIPcomputeVarUbLocal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPinferBinvarProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to fix */
   SCIP_Bool             fixedval,           /**< value to fix binary variable to */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the fixing */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(SCIPvarIsBinary(var));
   assert(fixedval == TRUE || fixedval == FALSE);
   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPinferBinvarProp", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsEQ(scip->set, lb, 0.0) || SCIPsetIsEQ(scip->set, lb, 1.0));
   assert(SCIPsetIsEQ(scip->set, ub, 0.0) || SCIPsetIsEQ(scip->set, ub, 1.0));
   assert(SCIPsetIsLE(scip->set, lb, ub));

   /* check, if variable is already fixed */
   if( (lb > 0.5) || (ub < 0.5) )
   {
      *infeasible = (fixedval == (lb < 0.5));

      return SCIP_OKAY;
   }

   /* apply the fixing */
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      {
         SCIP_Bool fixed;

         SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue,
               scip->cliquetable, (SCIP_Real)fixedval, infeasible, &fixed) );
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, 1.0,
               SCIP_BOUNDTYPE_LOWER, NULL, inferprop, inferinfo, FALSE) );
      }
      else
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, 0.0,
               SCIP_BOUNDTYPE_UPPER, NULL, inferprop, inferinfo, FALSE) );
      }
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes global lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current global bound; if possible, adjusts bound to integral value;
 *  also tightens the local bound, if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtightenVarLbGlobal", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* ignore tightenings of lower bounds to +infinity during solving process */
   if( SCIPisInfinity(scip, newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore lower bound tightening for %s from %e to +infinity\n", SCIPvarGetName(var),
         SCIPvarGetLbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   assert(scip->set->stage == SCIP_STAGE_PROBLEM || SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   /* bound changes of less than epsilon are ignored by SCIPvarChgLb or raise an assert in SCIPnodeAddBoundinfer,
    * so don't apply them even if force is set
    */
   if( SCIPsetIsEQ(scip->set, lb, newbound) || (!force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
            SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* coverity: unreachable code */
   if( tightened != NULL && lb < SCIPcomputeVarLbGlobal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes global upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current global bound; if possible, adjusts bound to integral value;
 *  also tightens the local bound, if the global bound is better than the local bound
 *
 *  @warning If SCIP is in presolving stage, it can happen that the internal variable array (which can be accessed via
 *           SCIPgetVars()) gets resorted.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note During presolving, an integer variable whose bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPtightenVarUbGlobal", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* ignore tightenings of upper bounds to -infinity during solving process */
   if( SCIPisInfinity(scip, -newbound) && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
#ifndef NDEBUG
      SCIPwarningMessage(scip, "ignore upper bound tightening for %s from %e to -infinity\n", SCIPvarGetName(var),
         SCIPvarGetUbLocal(var));
#endif
      return SCIP_OKAY;
   }

   /* get current bounds */
   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   assert(scip->set->stage == SCIP_STAGE_PROBLEM || SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   /* bound changes of less than epsilon are ignored by SCIPvarChgUb or raise an assert in SCIPnodeAddBoundinfer,
    * so don't apply them even if force is set
    */
   if( SCIPsetIsEQ(scip->set, ub, newbound) || (!force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
               scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
               SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( (SCIP_VARTYPE)var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, var, newbound,
            SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   /* coverity: unreachable code */
   if( tightened != NULL && ub > SCIPcomputeVarUbGlobal(scip, var) )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/* some simple variable functions implemented as defines */
#undef SCIPcomputeVarLbGlobal
#undef SCIPcomputeVarUbGlobal
#undef SCIPcomputeVarLbLocal
#undef SCIPcomputeVarUbLocal

/** for a multi-aggregated variable, returns the global lower bound computed by adding the global bounds from all aggregation variables
 *
 *  This global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbGlobal.
 *
 *  @return the global lower bound computed by adding the global bounds from all aggregation variables
 */
SCIP_Real SCIPcomputeVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrLbGlobal(var, scip->set);
   else
      return SCIPvarGetLbGlobal(var);
}

/** for a multi-aggregated variable, returns the global upper bound computed by adding the global bounds from all aggregation variables
 *
 *  This global bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbGlobal
 *
 *  @return the global upper bound computed by adding the global bounds from all aggregation variables
 */
SCIP_Real SCIPcomputeVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrUbGlobal(var, scip->set);
   else
      return SCIPvarGetUbGlobal(var);
}

/** for a multi-aggregated variable, returns the local lower bound computed by adding the local bounds from all aggregation variables
 *
 *  This local bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbLocal.
 *
 *  @return the local lower bound computed by adding the global bounds from all aggregation variables
 */
SCIP_Real SCIPcomputeVarLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrLbLocal(var, scip->set);
   else
      return SCIPvarGetLbLocal(var);
}

/** for a multi-aggregated variable, returns the local upper bound computed by adding the local bounds from all aggregation variables
 *
 *  This local bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is not updated if bounds of aggregation variables are changing
 *  calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbLocal.
 *
 *  @return the local upper bound computed by adding the global bounds from all aggregation variables
 */
SCIP_Real SCIPcomputeVarUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrUbLocal(var, scip->set);
   else
      return SCIPvarGetUbLocal(var);
}

/** for a multi-aggregated variable, gives the global lower bound computed by adding the global bounds from all
 *  aggregation variables, this global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPgetVarMultaggrLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   return SCIPvarGetMultaggrLbGlobal(var, scip->set);
}

/** for a multi-aggregated variable, gives the global upper bound computed by adding the global bounds from all
 *  aggregation variables, this upper bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPgetVarMultaggrUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   return SCIPvarGetMultaggrUbGlobal(var, scip->set);
}

/** for a multi-aggregated variable, gives the local lower bound computed by adding the local bounds from all
 *  aggregation variables, this lower bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPgetVarMultaggrLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   return SCIPvarGetMultaggrLbLocal(var, scip->set);
}

/** for a multi-aggregated variable, gives the local upper bound computed by adding the local bounds from all
 *  aggregation variables, this upper bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is
 *  not updated if bounds of aggregation variables are changing
 *
 *  calling this function for a non-multi-aggregated variable is not allowed
 */
SCIP_Real SCIPgetVarMultaggrUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
   return SCIPvarGetMultaggrUbLocal(var, scip->set);
}

/** returns solution value and index of variable lower bound that is closest to the variable's value in the given primal
 *  solution or current LP solution if no primal solution is given; returns an index of -1 if no variable lower bound is
 *  available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetVarClosestVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarClosestVlb", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPvarGetClosestVlb(var, sol, scip->set, scip->stat, closestvlb, closestvlbidx);

   return SCIP_OKAY;
}

/** returns solution value and index of variable upper bound that is closest to the variable's value in the given primal solution;
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetVarClosestVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvub,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest variable lower bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetVarClosestVub", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPvarGetClosestVub(var, sol, scip->set, scip->stat, closestvub, closestvubidx);

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  if z is non-continuous and 1/b not too small, the corresponding valid upper/lower bound
 *  z <= (x-d)/b or z >= (x-d)/b (depending on the sign of of b) is added, too;
 *  improves the global bounds of the variable and the vlb variable if possible
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   SCIP_Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   SCIP_Real             vlbconstant,        /**< constant d    in x >= b*z + d */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   int nlocalbdchgs;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarVlb", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarAddVlb(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob, scip->tree,
         scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, vlbvar, vlbcoef, vlbconstant,
         TRUE, infeasible, &nlocalbdchgs) );

   *nbdchgs = nlocalbdchgs;

   /* if x is not continuous we add a variable bound for z; do not add it if cofficient would be too small or we already
    * detected infeasibility
    */
   if( !(*infeasible) && SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPisZero(scip, 1.0/vlbcoef) )
   {
      if( vlbcoef > 0.0 )
      {
         /* if b > 0, we have a variable upper bound: x >= b*z + d  =>  z <= (x-d)/b */
         SCIP_CALL( SCIPvarAddVub(vlbvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var, 1.0/vlbcoef,
               -vlbconstant/vlbcoef, TRUE, infeasible, &nlocalbdchgs) );
      }
      else
      {
         /* if b < 0, we have a variable lower bound: x >= b*z + d  =>  z >= (x-d)/b */
         SCIP_CALL( SCIPvarAddVlb(vlbvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var, 1.0/vlbcoef,
               -vlbconstant/vlbcoef, TRUE, infeasible, &nlocalbdchgs) );
      }
      *nbdchgs += nlocalbdchgs;
   }

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  if z is non-continuous and 1/b not too small, the corresponding valid lower/upper bound
 *  z >= (x-d)/b or z <= (x-d)/b (depending on the sign of of b) is added, too;
 *  improves the global bounds of the variable and the vlb variable if possible
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   SCIP_Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   SCIP_Real             vubconstant,        /**< constant d    in x <= b*z + d */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   int nlocalbdchgs;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarVub", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarAddVub(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob, scip->tree,
         scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, vubvar, vubcoef, vubconstant, TRUE,
         infeasible, &nlocalbdchgs) );

   *nbdchgs = nlocalbdchgs;

   /* if x is not continuous we add a variable bound for z; do not add it if cofficient would be too small or we already
    * detected infeasibility
    */
   if( !(*infeasible) && SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPisZero(scip, 1.0/vubcoef) )
   {
      if( vubcoef > 0.0 )
      {
         /* if b < 0, we have a variable lower bound: x >= b*z + d  =>  z >= (x-d)/b */
         SCIP_CALL( SCIPvarAddVlb(vubvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var, 1.0/vubcoef,
               -vubconstant/vubcoef, TRUE, infeasible, &nlocalbdchgs) );
      }
      else
      {
         /* if b > 0, we have a variable upper bound: x >= b*z + d  =>  z <= (x-d)/b */
         SCIP_CALL( SCIPvarAddVub(vubvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var, 1.0/vubcoef,
               -vubconstant/vubcoef, TRUE, infeasible, &nlocalbdchgs) );
      }
      *nbdchgs += nlocalbdchgs;
   }

   return SCIP_OKAY;
}

/** informs binary variable x about a globally valid implication:  x == 0 or x == 1  ==>  y <= b  or  y >= b;
 *  also adds the corresponding implication or variable bound to the implied variable;
 *  if the implication is conflicting, the variable is fixed to the opposite value;
 *  if the variable is already fixed to the given value, the implication is performed immediately;
 *  if the implication is redundant with respect to the variables' global bounds, it is ignored
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarImplication(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER)
                                              *                          or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   SCIP_VAR* implprobvar;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarImplication", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(infeasible != NULL);
   *infeasible = FALSE;

   if ( nbdchgs != NULL )
      *nbdchgs = 0;

   if( !SCIPvarIsBinary(var) )
   {
      SCIPerrorMessage("can't add implication for nonbinary variable\n");
      return SCIP_INVALIDDATA;
   }

   implprobvar = SCIPvarGetProbvar(implvar);
   /* transform implication containing two binary variables to a clique; the condition ensures that the active representative
    * of implvar is actually binary
    */
   if( SCIPvarIsBinary(implvar) && (SCIPvarIsActive(implvar) || (implprobvar != NULL && SCIPvarIsBinary(implprobvar))) )
   {
      assert(SCIPisFeasEQ(scip, implbound, 1.0) || SCIPisFeasZero(scip, implbound));
      assert((impltype == SCIP_BOUNDTYPE_UPPER) == SCIPisFeasZero(scip, implbound));

      /* only add clique if implication is not redundant with respect to global bounds of the implication variable */
      if( (impltype == SCIP_BOUNDTYPE_LOWER && SCIPvarGetLbGlobal(implvar) < 0.5) ||
          (impltype == SCIP_BOUNDTYPE_UPPER && SCIPvarGetUbGlobal(implvar) > 0.5) )
      {
         SCIP_VAR* vars[2];
         SCIP_Bool vals[2];

         vars[0] = var;
         vars[1] = implvar;
         vals[0] = varfixing;
         vals[1] = (impltype == SCIP_BOUNDTYPE_UPPER);

         SCIP_CALL( SCIPaddClique(scip, vars, vals, 2, FALSE, infeasible, nbdchgs) );
      }

      return SCIP_OKAY;
   }

   /* the implication graph can only handle 'real' binary (SCIP_VARTYPE_BINARY) variables, therefore we transform the
    * implication in variable bounds, (lowerbound of y will be abbreviated by lby, upperbound equivlaent) the follwing
    * four cases are:
    *
    * 1. (x >= 1 => y >= b) => y >= (b - lby) * x + lby
    * 2. (x >= 1 => y <= b) => y <= (b - uby) * x + uby
    * 3. (x <= 0 => y >= b) => y >= (lby - b) * x + b
    * 4. (x <= 0 => y <= b) => y <= (uby - b) * x + b
    */
   if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
   {
      SCIP_Real lby;
      SCIP_Real uby;

      lby = SCIPvarGetLbGlobal(implvar);
      uby = SCIPvarGetUbGlobal(implvar);

      if( varfixing == TRUE )
      {
         if( impltype == SCIP_BOUNDTYPE_LOWER )
         {
            /* we return if the lower bound is infinity */
            if( SCIPisInfinity(scip, -lby) )
               return SCIP_OKAY;

            SCIP_CALL( SCIPvarAddVlb(implvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
                  scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var,
                  implbound - lby, lby, TRUE, infeasible, nbdchgs) );
         }
         else
         {
            /* we return if the upper bound is infinity */
            if( SCIPisInfinity(scip, uby) )
               return SCIP_OKAY;

            SCIP_CALL( SCIPvarAddVub(implvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
                  scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var,
                  implbound - uby, uby, TRUE, infeasible, nbdchgs) );
         }
      }
      else
      {
         if( impltype == SCIP_BOUNDTYPE_LOWER )
         {
            /* we return if the lower bound is infinity */
            if( SCIPisInfinity(scip, -lby) )
               return SCIP_OKAY;

            SCIP_CALL( SCIPvarAddVlb(implvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
                  scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var,
                  lby - implbound, implbound, TRUE, infeasible, nbdchgs) );
         }
         else
         {
            /* we return if the upper bound is infinity */
            if( SCIPisInfinity(scip, uby) )
               return SCIP_OKAY;

            SCIP_CALL( SCIPvarAddVub(implvar, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
                  scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, var,
                  uby - implbound, implbound, TRUE, infeasible, nbdchgs) );
         }
      }
   }
   else
   {
      SCIP_CALL( SCIPvarAddImplic(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventqueue, varfixing, implvar, impltype,
            implbound, TRUE, infeasible, nbdchgs) );
   }

   return SCIP_OKAY;
}

/** adds a clique information to SCIP, stating that at most one of the given binary variables can be set to 1;
 *  if a variable appears twice in the same clique, the corresponding implications are performed
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddClique(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars,              /**< number of variables in the clique */
   SCIP_Bool             isequation,         /**< is the clique an equation or an inequality? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddClique", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( nbdchgs != NULL )
      *nbdchgs = 0;

   if( nvars > 1 )
   {
      /* add the clique to the clique table */
      SCIP_CALL( SCIPcliquetableAdd(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, vars, values, nvars, isequation,
            infeasible, nbdchgs) );
   }

   return SCIP_OKAY;
}

/** relabels the given labels in-place in an increasing fashion: the first seen label is 0, the next label 1, etc...
 *
 *  @note every label equal to -1 is treated as a previously unseen, unique label and gets a new ordered label.
 */
static
SCIP_RETCODE relabelOrderConsistent(
   SCIP*const            scip,               /**< SCIP data structure */
   int*                  labels,             /**< current labels that will be overwritten */
   int const             nlabels,            /**< number of variables in the clique */
   int*                  nclasses            /**< pointer to store the total number of distinct labels */
   )
{
   SCIP_HASHMAP* classidx2newlabel;

   int classidx;
   int i;

   SCIP_CALL( SCIPhashmapCreate(&classidx2newlabel, SCIPblkmem(scip), nlabels) );

   classidx = 0;

   /* loop over labels to create local class indices that obey the variable order */
   for( i = 0; i < nlabels; ++i )
   {
      int currentlabel = labels[i];
      int localclassidx;

      /* labels equal to -1 are stored as singleton classes */
      if( currentlabel == -1 )
      {
         ++classidx;
         localclassidx = classidx;
      }
      else
      {
         assert(currentlabel >= 0);
         /* look up the class index image in the hash map; if it is not stored yet, new class index is created and stored */
         if( !SCIPhashmapExists(classidx2newlabel, (void*)(size_t)currentlabel) )
         {
            ++classidx;
            localclassidx = classidx;
            SCIP_CALL( SCIPhashmapInsertInt(classidx2newlabel, (void*)(size_t)currentlabel, classidx) ); /*lint !e571*/
         }
         else
         {
            localclassidx = SCIPhashmapGetImageInt(classidx2newlabel, (void*)(size_t)currentlabel); /*lint !e571*/
         }
      }
      assert(localclassidx - 1 >= 0);
      assert(localclassidx - 1 <= i);

      /* indices start with zero, but we have an offset of 1 because we cannot store 0 in a hashmap */
      labels[i] = localclassidx - 1;
   }

   assert(classidx > 0);
   assert(classidx <= nlabels);
   *nclasses = classidx;

   SCIPhashmapFree(&classidx2newlabel);

   return SCIP_OKAY;
}

/** sort the variables w.r.t. the given labels; thereby ensure the current order of the variables with the same label. */
static
SCIP_RETCODE labelSortStable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array */
   int*                  classlabels,        /**< array that contains a class label for every variable */
   SCIP_VAR**            sortedvars,         /**< array to store variables after stable sorting */
   int*                  sortedindices,      /**< array to store indices of sorted variables in the original vars array */
   int*                  classesstartposs,   /**< starting position array for each label class (must have size nclasses + 1) */
   int                   nvars,              /**< size of the vars arrays */
   int                   nclasses            /**< number of label classes */
   )
{
   SCIP_VAR*** varpointers;
   int** indexpointers;
   int* classcount;

   int nextpos;
   int c;
   int v;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(sortedindices != NULL);
   assert(classesstartposs != NULL);

   assert(nvars == 0 || vars != NULL);

   if( nvars == 0 )
      return SCIP_OKAY;

   assert(classlabels != NULL);
   assert(nclasses > 0);

   /* we first count all class cardinalities and allocate temporary memory for a bucket sort */
   SCIP_CALL( SCIPallocBufferArray(scip, &classcount, nclasses) );
   BMSclearMemoryArray(classcount, nclasses);

   /* first we count for each class the number of elements */
   for( v = nvars - 1; v >= 0; --v )
   {
      assert(0 <= classlabels[v] && classlabels[v] < nclasses);
      ++(classcount[classlabels[v]]);
   }

#ifndef NDEBUG
   BMSclearMemoryArray(sortedvars, nvars);
   BMSclearMemoryArray(sortedindices, nvars);
#endif
   SCIP_CALL( SCIPallocBufferArray(scip, &varpointers, nclasses) );
   SCIP_CALL( SCIPallocBufferArray(scip, &indexpointers, nclasses) );

   nextpos = 0;
   /* now we initialize all start pointers for each class, so they will be ordered */
   for( c = 0; c < nclasses; ++c )
   {
      /* to reach the goal that all variables of each class will be standing next to each other we will initialize the
       * starting pointers for each class by adding the cardinality of each class to the last class starting pointer
       * e.g. class1 has 4 elements and class2 has 3 elements then the starting pointer for class1 will be the pointer
       *      to sortedvars[0], the starting pointer to class2 will be the pointer to sortedvars[4] and to class3 it will be
       *      the pointer to sortedvars[7]
       */
      varpointers[c] = (SCIP_VAR**) (sortedvars + nextpos);
      indexpointers[c] = (int*) (sortedindices + nextpos);
      classesstartposs[c] = nextpos;
      assert(classcount[c] > 0);
      nextpos += classcount[c];
      assert(nextpos > 0);
   }
   assert(nextpos == nvars);
   classesstartposs[c] = nextpos;

   /* now we copy all variables to the right order */
   for( v = 0; v < nvars; ++v )
   {
      /* copy variable itself to the right position */
      *(varpointers[classlabels[v]]) = vars[v];  /*lint !e613*/
      ++(varpointers[classlabels[v]]);

      /* copy index */
      *(indexpointers[classlabels[v]]) = v;
      ++(indexpointers[classlabels[v]]);
   }

/* in debug mode, we ensure the correctness of the mapping */
#ifndef NDEBUG
   for( v = 0; v < nvars; ++v )
   {
      assert(sortedvars[v] != NULL);
      assert(sortedindices[v] >= 0);

      /* assert that the sorted indices map back to the correct variable in the original order */
      assert(vars[sortedindices[v]] == sortedvars[v]);
   }
#endif

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &indexpointers);
   SCIPfreeBufferArray(scip, &varpointers);
   SCIPfreeBufferArray(scip, &classcount);

   return SCIP_OKAY;
}


/* calculate clique partition for a maximal amount of comparisons on variables due to expensive algorithm
 * @todo: check for a good value, maybe it's better to check parts of variables
 */
#define MAXNCLIQUEVARSCOMP 1000000

/** calculates a partition of the given set of binary variables into cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same clique;
 *  the first variable is always assigned to clique 0, and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique at most 1 variables can be set to TRUE in a feasible solution;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
static
SCIP_RETCODE calcCliquePartitionGreedy(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   SCIP_Bool*const       values,             /**< clique value (TRUE or FALSE) for each variable in the clique */
   int const             nvars,              /**< number of variables in the array */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   )
{
   SCIP_VAR** cliquevars;
   SCIP_Bool* cliquevalues;
   int i;
   int maxncliquevarscomp;
   int ncliquevars;

   /* allocate temporary memory for storing the variables of the current clique */
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &cliquevars, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &cliquevalues, nvars) );

   /* initialize the cliquepartition array with -1 */
   for( i = nvars - 1; i >= 0; --i )
      cliquepartition[i] = -1;

   maxncliquevarscomp = (int) MIN(nvars * (SCIP_Longint)nvars, MAXNCLIQUEVARSCOMP);
   /* calculate the clique partition */
   *ncliques = 0;
   for( i = 0; i < nvars; ++i )
   {
      if( cliquepartition[i] == -1 )
      {
         int j;

         /* variable starts a new clique */
         cliquepartition[i] = *ncliques;
         cliquevars[0] = vars[i];
         cliquevalues[0] = values[i];
         ncliquevars = 1;

         /* if variable is not active (multi-aggregated or fixed), it cannot be in any clique */
         if( SCIPvarIsActive(vars[i]) && SCIPvarGetNCliques(vars[i], values[i]) > 0 )
         {
            /* greedily fill up the clique */
            for( j = i+1; j < nvars; ++j )
            {
               /* if variable is not active (multi-aggregated or fixed), it cannot be in any clique */
               if( cliquepartition[j] == -1 && SCIPvarIsActive(vars[j]) )
               {
                  int k;

                  /* check if every variable in the current clique can be extended by tmpvars[j] */
                  for( k = ncliquevars - 1; k >= 0; --k )
                  {
                     if( !SCIPvarsHaveCommonClique(vars[j], values[j], cliquevars[k], cliquevalues[k], FALSE) )
                        break;
                  }

                  if( k == -1 )
                  {
                     /* put the variable into the same clique */
                     cliquepartition[j] = cliquepartition[i];
                     cliquevars[ncliquevars] = vars[j];
                     cliquevalues[ncliquevars] = values[j];
                     ++ncliquevars;
                  }
               }
            }
         }

         /* this clique is finished */
         ++(*ncliques);
      }
      assert(cliquepartition[i] >= 0 && cliquepartition[i] < i+1);

      /* break if we reached the maximal number of comparisons */
      if( i * nvars > maxncliquevarscomp )
         break;
   }
   /* if we had to many variables fill up the cliquepartition and put each variable in a separate clique */
   for( ; i < nvars; ++i )
   {
      if( cliquepartition[i] == -1 )
      {
         cliquepartition[i] = *ncliques;
         ++(*ncliques);
      }
   }

   SCIPsetFreeBufferArray(scip->set, &cliquevalues);
   SCIPsetFreeBufferArray(scip->set, &cliquevars);

   return SCIP_OKAY;
}

/** calculates a partition of the given set of binary variables into cliques; takes into account independent clique components
 *
 *  The algorithm performs the following steps:
 *  - recomputes connected components of the clique table, if necessary
 *  - computes a clique partition for every connected component greedily.
 *  - relabels the resulting clique partition such that it satisfies the description below
 *
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same clique;
 *  the first variable is always assigned to clique 0, and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique at most 1 variables can be set to TRUE in a feasible solution;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcalcCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   )
{
   SCIP_VAR** tmpvars;

   SCIP_VAR** sortedtmpvars;
   SCIP_Bool* tmpvalues;
   SCIP_Bool* sortedtmpvalues;
   int* componentlabels;
   int* sortedindices;
   int* componentstartposs;
   int i;
   int c;

   int ncomponents;

   assert(scip != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || cliquepartition != NULL);
   assert(ncliques != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcalcCliquePartition", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( nvars == 0 )
   {
      *ncliques = 0;
      return SCIP_OKAY;
   }

   /* early abort if no cliques are present */
   if( SCIPgetNCliques(scip) == 0 )
   {
      for( i = 0; i < nvars; ++i )
         cliquepartition[i] = i;

      *ncliques = nvars;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &tmpvalues, nvars) );
   SCIP_CALL( SCIPsetDuplicateBufferArray(scip->set, &tmpvars, vars, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &componentlabels, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &sortedindices, nvars) );

   /* initialize the tmpvalues array */
   for( i = nvars - 1; i >= 0; --i )
   {
      tmpvalues[i] = TRUE;
      cliquepartition[i] = -1;
   }

   /* get corresponding active problem variables */
   SCIP_CALL( SCIPvarsGetProbvarBinary(&tmpvars, &tmpvalues, nvars) );

   ncomponents = -1;

   /* update clique components if necessary */
   if( SCIPcliquetableNeedsComponentUpdate(scip->cliquetable) )
   {
      SCIP_VAR** allvars;
      int nallbinvars;
      int nallintvars;
      int nallimplvars;

      SCIP_CALL( SCIPgetVarsData(scip, &allvars, NULL, &nallbinvars, &nallintvars, &nallimplvars, NULL) );

      SCIP_CALL( SCIPcliquetableComputeCliqueComponents(scip->cliquetable, scip->set, SCIPblkmem(scip), allvars, nallbinvars, nallintvars, nallimplvars) );
   }

   assert(!SCIPcliquetableNeedsComponentUpdate(scip->cliquetable));

   /* store the global clique component labels */
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPvarIsActive(tmpvars[i]) )
         componentlabels[i] = SCIPcliquetableGetVarComponentIdx(scip->cliquetable, tmpvars[i]);
      else
         componentlabels[i] = -1;
   }

   /* relabel component labels order consistent as prerequisite for a stable sort */
   SCIP_CALL( relabelOrderConsistent(scip, componentlabels, nvars, &ncomponents) );
   assert(ncomponents >= 1);
   assert(ncomponents <= nvars);

   /* allocate storage array for the starting positions of the components */
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &componentstartposs, ncomponents + 1) );

   /* stable sort the variables w.r.t. the component labels so that we can restrict the quadratic algorithm to the components */
   if( ncomponents > 1 )
   {
      SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &sortedtmpvars, nvars) );
      SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &sortedtmpvalues, nvars) );
      SCIP_CALL( labelSortStable(scip, tmpvars, componentlabels, sortedtmpvars, sortedindices, componentstartposs, nvars, ncomponents) );

      /* reassign the tmpvalues with respect to the sorting */
      for( i = 0; i < nvars; ++i )
      {
         assert(tmpvars[sortedindices[i]] == sortedtmpvars[i]);
         sortedtmpvalues[i] = tmpvalues[sortedindices[i]];
      }
   }
   else
   {
      /* if we have only one large connected component, skip the stable sorting and prepare the data differently */
      sortedtmpvars = tmpvars;
      sortedtmpvalues = tmpvalues;
      componentstartposs[0] = 0;
      componentstartposs[1] = nvars;

      /* sorted indices are the identity */
      for( i = 0; i < nvars; ++i )
         sortedindices[i] = i;
   }

   *ncliques = 0;
   /* calculate a greedy clique partition for each connected component */
   for( c = 0; c < ncomponents; ++c )
   {
      int* localcliquepartition;
      int nlocalcliques;
      int ncomponentvars;
      int l;

      /* extract the number of variables in this connected component */
      ncomponentvars = componentstartposs[c + 1] - componentstartposs[c];
      nlocalcliques = 0;

      /* allocate necessary memory to hold the intermediate component clique partition */
      SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &localcliquepartition, ncomponentvars) );

      /* call greedy clique algorithm for all component variables */
      SCIP_CALL( calcCliquePartitionGreedy(scip, &(sortedtmpvars[componentstartposs[c]]), &(sortedtmpvalues[componentstartposs[c]]),
            ncomponentvars, localcliquepartition, &nlocalcliques) );

      assert(nlocalcliques >= 1);
      assert(nlocalcliques <= ncomponentvars);

      /* store the obtained clique partition with an offset of ncliques for the original variables */
      for( l = componentstartposs[c]; l < componentstartposs[c + 1]; ++l )
      {
         int origvaridx = sortedindices[l];
         assert(cliquepartition[origvaridx] == -1);
         assert(localcliquepartition[l - componentstartposs[c]] <= l - componentstartposs[c]);
         cliquepartition[origvaridx] = localcliquepartition[l - componentstartposs[c]] + (*ncliques);
      }
      *ncliques += nlocalcliques;

      /* free the local clique partition */
      SCIPsetFreeBufferArray(scip->set, &localcliquepartition);
   }

   /* except in the two trivial cases, we have to ensure the order consistency of the partition indices */
   if( ncomponents > 1 && ncomponents < nvars )
   {
      int partitionsize;
      SCIP_CALL( relabelOrderConsistent(scip, cliquepartition, nvars, &partitionsize) );

      assert(partitionsize == *ncliques);
   }

   if( ncomponents > 1 )
   {
      SCIPsetFreeBufferArray(scip->set, &sortedtmpvalues);
      SCIPsetFreeBufferArray(scip->set, &sortedtmpvars);
   }

   /* use the greedy algorithm as a whole to verify the result on small number of variables */
#ifdef SCIP_DISABLED_CODE
   {
      int* debugcliquepartition;
      int ndebugcliques;

      SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &debugcliquepartition, nvars) );

      /* call greedy clique algorithm for all component variables */
      SCIP_CALL( calcCliquePartitionGreedy(scip, tmpvars, tmpvalues, nvars, debugcliquepartition, &ndebugcliques) );

      /* loop and compare the traditional greedy clique with  */
      for( i = 0; i < nvars; ++i )
         assert(i * nvars > MAXNCLIQUEVARSCOMP || cliquepartition[i] == debugcliquepartition[i]);

      SCIPsetFreeBufferArray(scip->set, &debugcliquepartition);
   }
#endif

   /* free temporary memory */
   SCIPsetFreeBufferArray(scip->set, &componentstartposs);
   SCIPsetFreeBufferArray(scip->set, &sortedindices);
   SCIPsetFreeBufferArray(scip->set, &componentlabels);
   SCIPsetFreeBufferArray(scip->set, &tmpvars);
   SCIPsetFreeBufferArray(scip->set, &tmpvalues);

   return SCIP_OKAY;
}

/** calculates a partition of the given set of binary variables into negated cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same negated clique;
 *  the first variable is always assigned to clique 0 and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique with n_c variables at least n_c-1 variables can be set to TRUE in a feasible solution;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcalcNegatedCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   )
{
   SCIP_VAR** negvars;
   int v;

   assert(scip != NULL);
   assert(cliquepartition != NULL || nvars == 0);
   assert(ncliques != NULL);

   if( nvars == 0 )
   {
      *ncliques = 0;
      return SCIP_OKAY;
   }
   assert(vars != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &negvars, nvars) );

   /* get all negated variables */
   for( v = nvars - 1; v >= 0; --v )
   {
      SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &(negvars[v])) );
   }

   /* calculate cliques on negated variables, which are "negated" cliques on normal variables array */
   SCIP_CALL( SCIPcalcCliquePartition( scip, negvars, nvars, cliquepartition, ncliques) );

   /* free temporary memory */
   SCIPsetFreeBufferArray(scip->set, &negvars);

   return SCIP_OKAY;
}


/** force SCIP to clean up all cliques; cliques do not get automatically cleaned up after presolving. Use
 *  this method to prevent inactive variables in cliques when retrieved via SCIPgetCliques()
 *
 *  @return SCIP_OKAY if everything worked, otherwise a suitable error code is passed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPcleanupCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            infeasible          /**< pointer to store if cleanup detected infeasibility */
   )
{
   int nlocalbdchgs;
   SCIP_Bool globalinfeasibility;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPcleanupCliques", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   globalinfeasibility = FALSE;
   nlocalbdchgs = 0;
   SCIP_CALL( SCIPcliquetableCleanup(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
         scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, &nlocalbdchgs,
         &globalinfeasibility) );

   if( infeasible != NULL )
      *infeasible = globalinfeasibility;

   if( globalinfeasibility )
      scip->stat->status = SCIP_STATUS_INFEASIBLE;

   return SCIP_OKAY;
}

/** gets the number of cliques in the clique table
 *
 *  @return number of cliques in the clique table
 *
 *  @note cliques do not get automatically cleaned up after presolving. Use SCIPcleanupCliques()
 *  to prevent inactive variables in cliques when retrieved via SCIPgetCliques(). This might reduce the number of cliques
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNCliques(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNCliques", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcliquetableGetNCliques(scip->cliquetable);
}

/** gets the number of cliques created so far by the cliquetable
 *
 *  @return number of cliques created so far by the cliquetable
 *
 *  @note cliques do not get automatically cleaned up after presolving. Use SCIPcleanupCliques()
 *  to prevent inactive variables in cliques when retrieved via SCIPgetCliques(). This might reduce the number of cliques
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNCliquesCreated(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNCliquesCreated", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcliquetableGetNCliquesCreated(scip->cliquetable);
}

/** gets the array of cliques in the clique table
 *
 *  @return array of cliques in the clique table
 *
 *  @note cliques do not get automatically cleaned up after presolving. Use SCIPcleanupCliques()
 *  to prevent inactive variables in cliques when retrieved via SCIPgetCliques(). This might reduce the number of cliques
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_CLIQUE** SCIPgetCliques(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetCliques", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPcliquetableGetCliques(scip->cliquetable);
}

/** returns whether there is a clique that contains both given variable/value pairs;
 *  the variables must be active binary variables;
 *  if regardimplics is FALSE, only the cliques in the clique table are looked at;
 *  if regardimplics is TRUE, both the cliques and the implications of the implication graph are regarded
 *
 *  @return TRUE, if there is a clique that contains both variable/clique pairs; FALSE, otherwise
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note a variable with it's negated variable are NOT! in a clique
 *  @note a variable with itself are in a clique
 */
SCIP_Bool SCIPhaveVarsCommonClique(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_Bool             value1,             /**< value of first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Bool             value2,             /**< value of second variable */
   SCIP_Bool             regardimplics       /**< should the implication graph also be searched for a clique? */
   )
{
   assert(scip != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(SCIPvarIsActive(var1));
   assert(SCIPvarIsActive(var2));
   assert(SCIPvarIsBinary(var1));
   assert(SCIPvarIsBinary(var2));

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPhaveVarsCommonClique", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* if both variables together have more cliques then actual cliques exist, then they have a common clique (in debug
    * mode we check this for correctness), otherwise we need to call the pairwise comparison method for these variables
    */
#ifndef NDEBUG
   assert((SCIPvarGetNCliques(var1, value1) + SCIPvarGetNCliques(var2, value2) > SCIPcliquetableGetNCliques(scip->cliquetable)) ? SCIPvarsHaveCommonClique(var1, value1, var2, value2, FALSE) : TRUE);
#endif

   return (SCIPvarGetNCliques(var1, value1) + SCIPvarGetNCliques(var2, value2) > SCIPcliquetableGetNCliques(scip->cliquetable)
      || SCIPvarsHaveCommonClique(var1, value1, var2, value2, regardimplics));
}

/** writes the clique graph to a gml file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @note there can be duplicated arcs in the output file
 *
 *  If @p writenodeweights is true, only nodes corresponding to variables that have a fractional value and only edges
 *  between such nodes are written.
 */
SCIP_RETCODE SCIPwriteCliqueGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname,              /**< name of file */
   SCIP_Bool             writenodeweights    /**< should we write weights of nodes? */
   )
{
   FILE* gmlfile;
   SCIP_HASHMAP* nodehashmap;
   SCIP_CLIQUE** cliques;
   SCIP_VAR** clqvars;
   SCIP_VAR** allvars;
   SCIP_Bool* clqvalues;
   char nodename[SCIP_MAXSTRLEN];
   int nallvars;
   int nbinvars;
   int nintvars;
   int nimplvars;
   int ncliques;
   int c;
   int v1;
   int v2;
   int id1;
   int id2;

   assert(scip != NULL);
   assert(fname != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPwriteCliqueGraph", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* get all active variables */
   SCIP_CALL( SCIPgetVarsData(scip, &allvars, &nallvars, &nbinvars, &nintvars, &nimplvars, NULL) );

   /* no possible variables for cliques exist */
   if( nbinvars + nimplvars == 0 )
      return SCIP_OKAY;

   ncliques = SCIPgetNCliques(scip);

   /* no cliques and do not wont to check for binary implications */
   if( ncliques == 0 )
      return SCIP_OKAY;

   /* open gml file */
   gmlfile = fopen(fname, "w");

   if( gmlfile == NULL )
   {
      SCIPerrorMessage("cannot open graph file <%s>\n", fname);
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   /* create the hash map */
   SCIP_CALL_FINALLY( SCIPhashmapCreate(&nodehashmap, SCIPblkmem(scip), nbinvars+nimplvars), fclose(gmlfile) );

   /* write starting of gml file */
   SCIPgmlWriteOpening(gmlfile, TRUE);

   cliques = SCIPgetCliques(scip);

   /* write nodes and arcs for all cliques */
   for( c = ncliques - 1; c >= 0; --c )
   {
      clqvalues = SCIPcliqueGetValues(cliques[c]);
      clqvars = SCIPcliqueGetVars(cliques[c]);

      for( v1 = SCIPcliqueGetNVars(cliques[c]) - 1; v1 >= 0; --v1 )
      {
	 id1 = clqvalues[v1] ? SCIPvarGetProbindex(clqvars[v1]) : (nallvars + SCIPvarGetProbindex(clqvars[v1]));

	 /* if corresponding node was not added yet, add it */
	 if( !SCIPhashmapExists(nodehashmap, (void*)(size_t)id1) )
	 {
            assert(id1 >= 0);
	    SCIP_CALL_FINALLY( SCIPhashmapInsertInt(nodehashmap, (void*)(size_t)id1, 1), fclose(gmlfile) ); /*lint !e571*/

	    (void) SCIPsnprintf(nodename, SCIP_MAXSTRLEN, "%s%s", (id1 >= nallvars ? "~" : ""), SCIPvarGetName(clqvars[v1]));

            /* write new gml node for new variable */
            if ( writenodeweights )
            {
               if ( ! SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, clqvars[v1])) )
                  SCIPgmlWriteNodeWeight(gmlfile, (unsigned int)id1, nodename, NULL, NULL, NULL, SCIPgetSolVal(scip, NULL, clqvars[v1]));
            }
            else
            {
               SCIPgmlWriteNode(gmlfile, (unsigned int)id1, nodename, NULL, NULL, NULL);
            }
	 }

	 for( v2 = SCIPcliqueGetNVars(cliques[c]) - 1; v2 >= 0; --v2 )
	 {
	    if( v1 == v2 )
	       continue;

	    id2 = clqvalues[v2] ? SCIPvarGetProbindex(clqvars[v2]) : (nallvars + SCIPvarGetProbindex(clqvars[v2]));

	    /* if corresponding node was not added yet, add it */
	    if( !SCIPhashmapExists(nodehashmap, (void*)(size_t)id2) )
	    {
               assert(id2 >= 0);
	       SCIP_CALL_FINALLY( SCIPhashmapInsertInt(nodehashmap, (void*)(size_t)id2, 1), fclose(gmlfile) ); /*lint !e571*/

	       (void) SCIPsnprintf(nodename, SCIP_MAXSTRLEN, "%s%s", (id2 >= nallvars ? "~" : ""), SCIPvarGetName(clqvars[v2]));

	       /* write new gml node for new variable */
               if ( writenodeweights )
               {
                  if ( ! SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, clqvars[v2])) )
                     SCIPgmlWriteNodeWeight(gmlfile, (unsigned int)id2, nodename, NULL, NULL, NULL, SCIPgetSolVal(scip, NULL, clqvars[v2]));
               }
               else
               {
                  SCIPgmlWriteNode(gmlfile, (unsigned int)id2, nodename, NULL, NULL, NULL);
               }
            }

	    /* write gml arc between resultant and operand */
            if ( ! writenodeweights || ! SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, clqvars[v2])) )
               SCIPgmlWriteArc(gmlfile, (unsigned int)id1, (unsigned int)id2, NULL, NULL);
	 }
      }
   }

   /* free the hash map */
   SCIPhashmapFree(&nodehashmap);

   SCIPgmlWriteClosing(gmlfile);
   fclose(gmlfile);

   return SCIP_OKAY;
}

/** Removes (irrelevant) variable from all its global structures, i.e. cliques, implications and variable bounds.
 *  This is an advanced method which should be used with care.
 *
 *  @return SCIP_OKAY if everything worked, otherwise a suitable error code is passed
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPremoveVarFromGlobalStructures(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to remove from global structures */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPremoveVarFromGlobalStructures", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* mark the variable as deletable from global structures - This is necessary for the delayed clean up of cliques */
   SCIPvarMarkDeleteGlobalStructures(var);

   /* remove variable from all its cliques, implications, and variable bounds */
   SCIP_CALL( SCIPvarRemoveCliquesImplicsVbs(var, SCIPblkmem(scip), scip->cliquetable, scip->set, TRUE, FALSE, TRUE) );

   return SCIP_OKAY;
}

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarChgBranchFactor(var, scip->set, branchfactor) );

   return SCIP_OKAY;
}

/** scales the branch factor of the variable with the given value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPscaleVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             scale               /**< factor to scale variable's branching factor with */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPscaleVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarChgBranchFactor(var, scip->set, scale * SCIPvarGetBranchFactor(var)) );

   return SCIP_OKAY;
}

/** adds the given value to the branch factor of the variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             addfactor           /**< value to add to the branch factor of the variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarChgBranchFactor(var, scip->set, addfactor + SCIPvarGetBranchFactor(var)) );

   return SCIP_OKAY;
}

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *
 * @note the default branching priority is 0
 */
SCIP_RETCODE SCIPchgVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< branch priority of the variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( SCIPisTransformed(scip)  )
   {
      assert(scip->branchcand != NULL);

      /* inform the pseudo branch candidates that the branch priority changes and change the branch priority */
      SCIP_CALL( SCIPbranchcandUpdateVarBranchPriority(scip->branchcand, scip->set, var, branchpriority) );
   }
   else
   {
      /* change the branching priority of the variable */
      SCIP_CALL( SCIPvarChgBranchPriority(var, branchpriority) );
   }

   return SCIP_OKAY;
}

/** changes the branch priority of the variable to the given value, if it is larger than the current priority
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPupdateVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< new branch priority of the variable, if it is larger than current priority */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( branchpriority > SCIPvarGetBranchPriority(var) )
   {
      SCIP_CALL( SCIPvarChgBranchPriority(var, branchpriority) );
   }

   return SCIP_OKAY;
}

/** adds the given value to the branch priority of the variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   addpriority         /**< value to add to the branch priority of the variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPvarChgBranchPriority(var, addpriority + SCIPvarGetBranchPriority(var)) );

   return SCIP_OKAY;
}

/** sets the branch direction of the variable (-1: prefer downwards branch, 0: automatic selection, +1: prefer upwards
 *  branch)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarBranchDirection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarBranchDirection", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   SCIP_CALL( SCIPvarChgBranchDirection(var, branchdirection) );

   return SCIP_OKAY;
}

/** tightens the variable bounds due to a new variable type */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_VARTYPE          vartype,            /**< new type of variable */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected (, due to
                                              *   integrality condition of the new variable type) */
   )
{
   assert(scip != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM || SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
   assert(scip->set->stage == SCIP_STAGE_PROBLEM || SCIPvarIsTransformed(var));
   assert(var->scip == scip);

   *infeasible = FALSE;

   /* adjusts bounds if the variable type changed form continuous to non-continuous (integral) */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && vartype != SCIP_VARTYPE_CONTINUOUS )
   {
      SCIP_Bool tightened;

      /* we adjust variable bounds to integers first, since otherwise a later bound tightening with a fractional old
       * bound may give an assert because SCIP expects non-continuous variables to have non-fractional bounds
       *
       * we adjust bounds with a fractionality within [eps,feastol] only if the resulting bound change is a bound
       * tightening, because relaxing bounds may not be allowed
       */
      if( !SCIPisFeasIntegral(scip, SCIPvarGetLbGlobal(var)) ||
         (!SCIPisIntegral(scip, SCIPvarGetLbGlobal(var)) && SCIPvarGetLbGlobal(var) < SCIPfeasCeil(scip, SCIPvarGetLbGlobal(var))) ||
         (!SCIPsetIsEQ(scip->set, SCIPvarGetLbGlobal(var), SCIPfeasCeil(scip, SCIPvarGetLbGlobal(var))) &&
          SCIPvarGetLbGlobal(var) < SCIPfeasCeil(scip, SCIPvarGetLbGlobal(var)))
        )
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, SCIPfeasCeil(scip, SCIPvarGetLbGlobal(var)), TRUE, infeasible, &tightened) );
         if( *infeasible )
            return SCIP_OKAY;

         /* the only reason for not applying a forced boundchange is when the new bound is reduced because the variables upper bound is below the new bound
          * in a concrete case, lb == ub == 100.99999001; even though within feastol of 101, the lower bound cannot be tighented to 101 due to the upper bound
          */
         assert(tightened || SCIPisFeasLE(scip, SCIPvarGetUbGlobal(var), SCIPfeasCeil(scip, SCIPvarGetLbGlobal(var))));
      }
      if( !SCIPisFeasIntegral(scip, SCIPvarGetUbGlobal(var)) ||
         (!SCIPisIntegral(scip, SCIPvarGetUbGlobal(var)) && SCIPvarGetUbGlobal(var) > SCIPfeasFloor(scip, SCIPvarGetUbGlobal(var)))
        )
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, SCIPfeasFloor(scip, SCIPvarGetUbGlobal(var)), TRUE, infeasible, &tightened) );
         if( *infeasible )
            return SCIP_OKAY;

         assert(tightened || SCIPisFeasGE(scip, SCIPvarGetLbGlobal(var), SCIPfeasFloor(scip, SCIPvarGetUbGlobal(var))));
      }
   }

   return SCIP_OKAY;
}

/** changes type of variable in the problem;
 *
 *  @warning This type change might change the variable array returned from SCIPgetVars() and SCIPgetVarsData();
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *
 *  @note If SCIP is already beyond the SCIP_STAGE_PROBLEM and a original variable is passed, the variable type of the
 *        corresponding transformed variable is changed; the type of the original variable does not change
 *
 *  @note If the type changes from a continuous variable to a non-continuous variable the bounds of the variable get
 *        adjusted w.r.t. to integrality information
 */
SCIP_RETCODE SCIPchgVarType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_VARTYPE          vartype,            /**< new type of variable */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected (, due to
                                              *   integrality condition of the new variable type) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarType", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(var != NULL);
   assert(var->scip == scip);

   if( SCIPvarIsNegated(var) )
   {
      SCIPdebugMsg(scip, "upgrading type of negated variable <%s> from %d to %d\n", SCIPvarGetName(var), SCIPvarGetType(var), vartype);
      var = SCIPvarGetNegationVar(var);
   }
#ifndef NDEBUG
   else
   {
      if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM )
      {
         SCIPdebugMsg(scip, "upgrading type of variable <%s> from %d to %d\n", SCIPvarGetName(var), SCIPvarGetType(var), vartype);
      }
   }
#endif

   /* change variable type */
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));

      /* first adjust the variable due to new integrality information */
      SCIP_CALL( tightenBounds(scip, var, vartype, infeasible) );

      /* second change variable type */
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         SCIP_CALL( SCIPprobChgVarType(scip->origprob, scip->mem->probmem, scip->set, scip->primal, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, var, vartype) );
      }
      else
      {
         SCIP_CALL( SCIPvarChgType(var, scip->mem->probmem, scip->set, scip->primal, scip->lp,
            scip->eventqueue, vartype) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPvarIsTransformed(var) )
      {
         SCIP_VAR* transvar;

         SCIP_CALL( SCIPgetTransformedVar(scip, var, &transvar) );
         assert(transvar != NULL);

         /* recall method with transformed variable */
         SCIP_CALL( SCIPchgVarType(scip, transvar, vartype, infeasible) );
         return SCIP_OKAY;
      }

      /* first adjust the variable due to new integrality information */
      SCIP_CALL( tightenBounds(scip, var, vartype, infeasible) );

      /* second change variable type */
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         SCIP_CALL( SCIPprobChgVarType(scip->transprob, scip->mem->probmem, scip->set, scip->primal, scip->lp,
            scip->branchcand, scip->eventqueue, scip->cliquetable, var, vartype) );
      }
      else
      {
         SCIP_CALL( SCIPvarChgType(var, scip->mem->probmem, scip->set, scip->primal, scip->lp,
            scip->eventqueue, vartype) );
      }
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** in problem creation and solving stage, both bounds of the variable are set to the given value;
 *  in presolving stage, the variable is converted into a fixed variable, and bounds are changed respectively;
 *  conversion into a fixed variable changes the vars array returned from SCIPgetVars() and SCIPgetVarsData(),
 *  and also renders arrays returned from the SCIPvarGetImpl...() methods invalid
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPfixVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to fix */
   SCIP_Real             fixedval,           /**< value to fix variable to */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   )
{
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(fixed != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPfixVar", FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   *fixed = FALSE;

   /* in the problem creation stage, modify the bounds as requested, independently from the current bounds */
   if( scip->set->stage != SCIP_STAGE_PROBLEM )
   {
      if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPsetIsFeasIntegral(scip->set, fixedval))
         || SCIPsetIsFeasLT(scip->set, fixedval, SCIPvarGetLbLocal(var))
         || SCIPsetIsFeasGT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
      {
         *infeasible = !SCIPsetIsFeasEQ(scip->set, fixedval, SCIPvarGetLbLocal(var));
         return SCIP_OKAY;
      }
   }
   else
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* in the problem creation stage, modify the bounds as requested, independently from the current bounds;
       * we have to make sure, that the order of the bound changes does not intermediately produce an invalid
       * interval lb > ub
       */
      if( fixedval <= SCIPvarGetLbLocal(var) )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, fixedval) );
         SCIP_CALL( SCIPchgVarUb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, fixedval) );
         SCIP_CALL( SCIPchgVarLb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      {
         SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue,
               scip->cliquetable, fixedval, infeasible, fixed) );
         return SCIP_OKAY;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      if( SCIPsetIsFeasGT(scip->set, fixedval, SCIPvarGetLbLocal(var)) )
      {
         if( SCIPsetIsFeasGT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
         {
            *infeasible = TRUE;
            return SCIP_OKAY;
         }
         else
         {
            SCIP_CALL( SCIPchgVarLb(scip, var, fixedval) );
            *fixed = TRUE;
         }
      }
      if( SCIPsetIsFeasLT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
      {
         if( SCIPsetIsFeasLT(scip->set, fixedval, SCIPvarGetLbLocal(var)) )
         {
            *infeasible = TRUE;
            return SCIP_OKAY;
         }
         else
         {
            SCIP_CALL( SCIPchgVarUb(scip, var, fixedval) );
            *fixed = TRUE;
         }
      }
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** From a given equality a*x + b*y == c, aggregates one of the variables and removes it from the set of
 *  active problem variables. This changes the vars array returned from SCIPgetVars() and SCIPgetVarsData(),
 *  and also renders the arrays returned from the SCIPvarGetImpl...() methods for the two variables invalid.
 *  In the first step, the equality is transformed into an equality with active problem variables
 *  a'*x' + b'*y' == c'. If x' == y', this leads to the detection of redundancy if a' == -b' and c' == 0,
 *  of infeasibility, if a' == -b' and c' != 0, or to a variable fixing x' == c'/(a'+b') (and possible
 *  infeasibility) otherwise.
 *  In the second step, the variable to be aggregated is chosen among x' and y', prefering a less strict variable
 *  type as aggregation variable (i.e. continuous variables are preferred over implicit integers, implicit integers
 *  over integers, and integers over binaries). If none of the variables is continuous, it is tried to find an integer
 *  aggregation (i.e. integral coefficients a'' and b'', such that a''*x' + b''*y' == c''). This can lead to
 *  the detection of infeasibility (e.g. if c'' is fractional), or to a rejection of the aggregation (denoted by
 *  aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 *
 *  The output flags have the following meaning:
 *  - infeasible: the problem is infeasible
 *  - redundant:  the equality can be deleted from the constraint set
 *  - aggregated: the aggregation was successfully performed (the variables were not aggregated before)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_PRESOLVING
 */
SCIP_RETCODE SCIPaggregateVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             varx,               /**< variable x in equality a*x + b*y == c */
   SCIP_VAR*             vary,               /**< variable y in equality a*x + b*y == c */
   SCIP_Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   SCIP_Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   SCIP_Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            redundant,          /**< pointer to store whether the equality is (now) redundant */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_Real constantx;
   SCIP_Real constanty;

   assert(infeasible != NULL);
   assert(redundant != NULL);
   assert(aggregated != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaggregateVars", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   *redundant = FALSE;
   *aggregated = FALSE;

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot aggregate variables during probing\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);

   /* do not perform aggregation if it is globally deactivated */
   if( scip->set->presol_donotaggr )
      return SCIP_OKAY;

   /* get the corresponding equality in active problem variable space:
    * transform both expressions "a*x + 0" and "b*y + 0" into problem variable space
    */
   constantx = 0.0;
   constanty = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&varx, scip->set, &scalarx, &constantx) );
   SCIP_CALL( SCIPvarGetProbvarSum(&vary, scip->set, &scalary, &constanty) );

   /* we cannot aggregate multi-aggregated variables */
   if( SCIPvarGetStatus(varx) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(vary) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   /* move the constant to the right hand side to acquire the form "a'*x' + b'*y' == c'" */
   rhs -= (constantx + constanty);

   /* if a scalar is zero, treat the variable as fixed-to-zero variable */
   if( SCIPsetIsZero(scip->set, scalarx) )
      varx = NULL;
   if( SCIPsetIsZero(scip->set, scalary) )
      vary = NULL;

   /* capture the special cases that less than two variables are left, due to resolutions to a fixed variable or
    * to the same active variable
    */
   if( varx == NULL && vary == NULL )
   {
      /* both variables were resolved to fixed variables */
      *infeasible = !SCIPsetIsZero(scip->set, rhs);
      *redundant = TRUE;
   }
   else if( varx == NULL )
   {
      assert(SCIPsetIsZero(scip->set, scalarx));
      assert(!SCIPsetIsZero(scip->set, scalary));

      /* variable x was resolved to fixed variable: variable y can be fixed to c'/b' */
      SCIP_CALL( SCIPvarFix(vary, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue,
            scip->cliquetable, rhs/scalary, infeasible, aggregated) );
      *redundant = TRUE;
   }
   else if( vary == NULL )
   {
      assert(SCIPsetIsZero(scip->set, scalary));
      assert(!SCIPsetIsZero(scip->set, scalarx));

      /* variable y was resolved to fixed variable: variable x can be fixed to c'/a' */
      SCIP_CALL( SCIPvarFix(varx, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue,
            scip->cliquetable, rhs/scalarx, infeasible, aggregated) );
      *redundant = TRUE;
   }
   else if( varx == vary )
   {
      /* both variables were resolved to the same active problem variable: this variable can be fixed */
      scalarx += scalary;
      if( SCIPsetIsZero(scip->set, scalarx) )
      {
         /* left hand side of equality is zero: equality is potentially infeasible */
         *infeasible = !SCIPsetIsZero(scip->set, rhs);
      }
      else
      {
         /* sum of scalars is not zero: fix variable x' == y' to c'/(a'+b') */
         SCIP_CALL( SCIPvarFix(varx, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
               scip->primal, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventfilter, scip->eventqueue,
               scip->cliquetable, rhs/scalarx, infeasible, aggregated) );
      }
      *redundant = TRUE;
   }
   else
   {
      /* both variables are different active problem variables, and both scalars are non-zero: try to aggregate them */
      SCIP_CALL( SCIPvarTryAggregateVars(scip->set, scip->mem->probmem, scip->stat, scip->transprob, scip->origprob,
            scip->primal, scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventfilter,
            scip->eventqueue, varx, vary, scalarx, scalary, rhs, infeasible, aggregated) );
      *redundant = *aggregated;
   }

   return SCIP_OKAY;
}

/** converts variable into multi-aggregated variable; this changes the variable array returned from
 *  SCIPgetVars() and SCIPgetVarsData();
 *
 *  @warning The integrality condition is not checked anymore on the multi-aggregated variable. You must not
 *           multi-aggregate an integer variable without being sure, that integrality on the aggregation variables
 *           implies integrality on the aggregated variable.
 *
 *  The output flags have the following meaning:
 *  - infeasible: the problem is infeasible
 *  - aggregated: the aggregation was successfully performed (the variables were not aggregated before)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can only be called if @p scip is in stage \ref SCIP_STAGE_PRESOLVING
 */
SCIP_RETCODE SCIPmultiaggregateVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable x to aggregate */
   int                   naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPmultiaggregateVar", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(var->scip == scip);

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot multi-aggregate variables during probing\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);

   SCIP_CALL( SCIPvarMultiaggregate(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->primal, scip->tree, scip->reopt, scip->lp, scip->cliquetable, scip->branchcand, scip->eventfilter,
         scip->eventqueue, naggvars, aggvars, scalars, constant, infeasible, aggregated) );

   return SCIP_OKAY;
}

/** returns whether aggregation of variables is not allowed */
SCIP_Bool SCIPdoNotAggr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->set->presol_donotaggr;
}

/** returns whether multi-aggregation is disabled */
SCIP_Bool SCIPdoNotMultaggr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->set->presol_donotmultaggr;
}

/** returns whether variable is not allowed to be aggregated */
SCIP_Bool SCIPdoNotAggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable x to aggregate */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(var->scip == scip);

   return scip->set->presol_donotaggr || SCIPvarDoNotAggr(var);
}

/** returns whether variable is not allowed to be multi-aggregated */
SCIP_Bool SCIPdoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable x to aggregate */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(var->scip == scip);

   return scip->set->presol_donotmultaggr || SCIPvarDoNotMultaggr(var);
}

/** returns whether dual reductions are allowed during propagation and presolving
 *
 *  @deprecated Please use SCIPallowStrongDualReds()
 */
SCIP_Bool SCIPallowDualReds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return !scip->set->reopt_enable && scip->set->misc_allowstrongdualreds;
}

/** returns whether strong dual reductions are allowed during propagation and presolving
 *
 *  @note A reduction is called strong dual, if it may discard feasible/optimal solutions, but leaves at least one
 *        optimal solution intact. Often such reductions are based on analyzing the objective function and variable
 *        locks.
 */
SCIP_Bool SCIPallowStrongDualReds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return !scip->set->reopt_enable && scip->set->misc_allowstrongdualreds;
}

/** returns whether propagation w.r.t. current objective is allowed
 *
 *  @deprecated Please use SCIPallowWeakDualReds()
 */
SCIP_Bool SCIPallowObjProp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return !scip->set->reopt_enable && scip->set->misc_allowweakdualreds;
}

/** returns whether weak dual reductions are allowed during propagation and presolving
 *
 *  @note A reduction is called weak dual, if it may discard feasible solutions, but leaves at all optimal solutions
 *        intact. Often such reductions are based on analyzing the objective function, reduced costs, and/or dual LPs.
 */
SCIP_Bool SCIPallowWeakDualReds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return !scip->set->reopt_enable && scip->set->misc_allowweakdualreds;
}

/** marks the variable that it must not be aggregated
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *
 *  @note There exists no "unmark" method since it has to be ensured that if a plugin requires that a variable is not
 *        aggregated that this is will be the case.
 */
SCIP_RETCODE SCIPmarkDoNotAggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to delete */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(var->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPmarkDoNotAggrVar", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPvarMarkDoNotAggr(var) );

   return SCIP_OKAY;
}

/** marks the variable that it must not be multi-aggregated
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *
 *  @note There exists no "unmark" method since it has to be ensured that if a plugin requires that a variable is not
 *        multi-aggregated that this is will be the case.
 */
SCIP_RETCODE SCIPmarkDoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to delete */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(var->scip == scip);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPmarkDoNotMultaggrVar", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE) );

   SCIP_CALL( SCIPvarMarkDoNotMultaggr(var) );

   return SCIP_OKAY;
}

/** enables the collection of statistics for a variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPenableVarHistory(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPenableVarHistory", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPstatEnableVarHistory(scip->stat);
}

/** disables the collection of any statistic for a variable
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPdisableVarHistory(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPdisableVarHistory", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPstatDisableVarHistory(scip->stat);
}

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of "solvaldelta" in the
 *  variable's solution value and resulting change of "objdelta" in the in the LP's objective value;
 *  the update is ignored, if the objective value difference is infinite
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPupdateVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateVarPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPsetIsInfinity(scip->set, 2*objdelta) ) /* differences  infinity - eps  should also be treated as infinity */
   {
      if( scip->set->branch_divingpscost || (!scip->lp->diving && !SCIPtreeProbing(scip->tree)) )
      {
         SCIP_CALL( SCIPvarUpdatePseudocost(var, scip->set, scip->stat, solvaldelta, objdelta, weight) );
      }
   }

   return SCIP_OKAY;
}

/** gets the variable's pseudo cost value for the given change of the variable's LP value
 *
 *  @return the variable's pseudo cost value for the given change of the variable's LP value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   assert( var->scip == scip );

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostVal", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarGetPseudocost(var, scip->stat, solvaldelta);
}

/** gets the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostValCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   assert( var->scip == scip );

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostValCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarGetPseudocostCurrentRun(var, scip->stat, solvaldelta);
}

/** gets the variable's pseudo cost value for the given direction
 *
 *  @return the variable's pseudo cost value for the given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocost", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert(var->scip == scip);

   return SCIPvarGetPseudocost(var, scip->stat, dir == SCIP_BRANCHDIR_DOWNWARDS ? -1.0 : 1.0);
}

/** gets the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert(var->scip == scip);

   return SCIPvarGetPseudocostCurrentRun(var, scip->stat, dir == SCIP_BRANCHDIR_DOWNWARDS ? -1.0 : 1.0);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction
 *
 *  @return the variable's (possible fractional) number of pseudo cost updates for the given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostCount", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert(var->scip == scip);

   return SCIPvarGetPseudocostCount(var, dir);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostCountCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert(var->scip == scip);

   return SCIPvarGetPseudocostCountCurrentRun(var, dir);
}

/** get pseudo cost variance of the variable, either for entire solve or only for current branch and bound run
 *
 *  @return returns the (corrected) variance of pseudo code information collected so far.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostVariance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Bool             onlycurrentrun      /**< only for pseudo costs of current branch and bound run */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostVariance", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   assert(var->scip == scip);

   return SCIPvarGetPseudocostVariance(var, dir, onlycurrentrun);
}

/** calculates a confidence bound for this variable under the assumption of normally distributed pseudo costs
 *
 *  The confidence bound \f$ \theta \geq 0\f$ denotes the interval borders \f$ [X - \theta, \ X + \theta]\f$, which contains
 *  the true pseudo costs of the variable, i.e., the expected value of the normal distribution, with a probability
 *  of 2 * clevel - 1.
 *
 *  @return value of confidence bound for this variable
 */
SCIP_Real SCIPcalculatePscostConfidenceBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_BRANCHDIR        dir,                /**< the branching direction for the confidence bound */
   SCIP_Bool             onlycurrentrun,     /**< should only the current run be taken into account */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for the interval */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPcalculatePscostConfidenceBound", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarCalcPscostConfidenceBound(var, scip->set, dir, onlycurrentrun, clevel);
}

/** check if variable pseudo-costs have a significant difference in location. The significance depends on
 *  the choice of \p clevel and on the kind of tested hypothesis. The one-sided hypothesis, which
 *  should be rejected, is that fracy * mu_y >= fracx * mu_x, where mu_y and mu_x denote the
 *  unknown location means of the underlying pseudo-cost distributions of x and y.
 *
 *  This method is applied best if variable x has a better pseudo-cost score than y. The method hypothesizes that y were actually
 *  better than x (despite the current information), meaning that y can be expected to yield branching
 *  decisions as least as good as x in the long run. If the method returns TRUE, the current history information is
 *  sufficient to safely rely on the alternative hypothesis that x yields indeed a better branching score (on average)
 *  than y.
 *
 *  @note The order of x and y matters for the one-sided hypothesis
 *
 *  @note set \p onesided to FALSE if you are not sure which variable is better. The hypothesis tested then reads
 *        fracy * mu_y == fracx * mu_x vs the alternative hypothesis fracy * mu_y != fracx * mu_x.
 *
 *  @return TRUE if the hypothesis can be safely rejected at the given confidence level
 */
SCIP_Bool SCIPsignificantVarPscostDifference(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             varx,               /**< variable x */
   SCIP_Real             fracx,              /**< the fractionality of variable x */
   SCIP_VAR*             vary,               /**< variable y */
   SCIP_Real             fracy,              /**< the fractionality of variable y */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel,             /**< confidence level for rejecting hypothesis */
   SCIP_Bool             onesided            /**< should a one-sided hypothesis y >= x be tested? */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPsignificantVarPscostDifference", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarSignificantPscostDifference(scip->set, scip->stat, varx, fracx, vary, fracy, dir, clevel, onesided);
}

/** tests at a given confidence level whether the variable pseudo-costs only have a small probability to
 *  exceed a \p threshold. This is useful to determine if past observations provide enough evidence
 *  to skip an expensive strong-branching step if there is already a candidate that has been proven to yield an improvement
 *  of at least \p threshold.
 *
 *  @note use \p clevel to adjust the level of confidence. For SCIP_CONFIDENCELEVEL_MIN, the method returns TRUE if
 *        the estimated probability to exceed \p threshold is less than 25 %.
 *
 *  @see  SCIP_Confidencelevel for a list of available levels. The used probability limits refer to the one-sided levels
 *        of confidence.
 *
 *  @return TRUE if the variable pseudo-cost probabilistic model is likely to be smaller than \p threshold
 *          at the given confidence level \p clevel.
 */
SCIP_Bool SCIPpscostThresholdProbabilityTest(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable x */
   SCIP_Real             frac,               /**< the fractionality of variable x */
   SCIP_Real             threshold,          /**< the threshold to test against */
   SCIP_BRANCHDIR        dir,                /**< branching direction */
   SCIP_CONFIDENCELEVEL  clevel              /**< confidence level for rejecting hypothesis */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPpscostThresholdProbabilityTest", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarPscostThresholdProbabilityTest(scip->set, scip->stat, var, frac, threshold, dir, clevel);
}

/** check if the current pseudo cost relative error in a direction violates the given threshold. The Relative
 *  Error is calculated at a specific confidence level
 *
 *  @return TRUE if relative error in variable pseudo costs is smaller than \p threshold
 */
SCIP_Bool SCIPisVarPscostRelerrorReliable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable in question */
   SCIP_Real             threshold,          /**< threshold for relative errors to be considered reliable (enough) */
   SCIP_CONFIDENCELEVEL  clevel              /**< a given confidence level */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisVarPscostRelerrorReliable", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarIsPscostRelerrorReliable(var, scip->set, scip->stat, threshold, clevel);
}

/** gets the variable's pseudo cost score value for the given LP solution value
 *
 *  @return the variable's pseudo cost score value for the given LP solution value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   )
{
   SCIP_Real downsol;
   SCIP_Real upsol;
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   downsol = SCIPsetFeasCeil(scip->set, solval-1.0);
   upsol = SCIPsetFeasFloor(scip->set, solval+1.0);
   pscostdown = SCIPvarGetPseudocost(var, scip->stat, downsol-solval);
   pscostup = SCIPvarGetPseudocost(var, scip->stat, upsol-solval);

   return SCIPbranchGetScore(scip->set, var, pscostdown, pscostup);
}

/** gets the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 *
 *  @return the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarPseudocostScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   )
{
   SCIP_Real downsol;
   SCIP_Real upsol;
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarPseudocostScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   downsol = SCIPsetFeasCeil(scip->set, solval-1.0);
   upsol = SCIPsetFeasFloor(scip->set, solval+1.0);
   pscostdown = SCIPvarGetPseudocostCurrentRun(var, scip->stat, downsol-solval);
   pscostup = SCIPvarGetPseudocostCurrentRun(var, scip->stat, upsol-solval);

   return SCIPbranchGetScore(scip->set, var, pscostdown, pscostup);
}

/** returns the variable's VSIDS value
 *
 *  @return the variable's VSIDS value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarVSIDS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarVSIDS", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( dir != SCIP_BRANCHDIR_DOWNWARDS && dir != SCIP_BRANCHDIR_UPWARDS )
   {
      SCIPerrorMessage("invalid branching direction %d when asking for VSIDS value\n", dir);
      return SCIP_INVALID;
   }

   return SCIPvarGetVSIDS(var, scip->stat, dir);
}

/** returns the variable's VSIDS value only using conflicts of the current run
 *
 *  @return the variable's VSIDS value only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarVSIDSCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarVSIDSCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( dir != SCIP_BRANCHDIR_DOWNWARDS && dir != SCIP_BRANCHDIR_UPWARDS )
   {
      SCIPerrorMessage("invalid branching direction %d when asking for VSIDS value\n", dir);
      return SCIP_INVALID;
   }

   return SCIPvarGetVSIDSCurrentRun(var, scip->stat, dir);
}

/** returns the variable's conflict score value
 *
 *  @return the variable's conflict score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarConflictScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarConflictScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   downscore = SCIPvarGetVSIDS(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetVSIDS(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's conflict score value only using conflicts of the current run
 *
 *  @return the variable's conflict score value only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarConflictScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarConflictScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   downscore = SCIPvarGetVSIDSCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetVSIDSCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's conflict length score
 *
 *  @return the variable's conflict length score
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarConflictlengthScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarConflictlengthScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   downscore = SCIPvarGetAvgConflictlength(var, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetAvgConflictlength(var, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's conflict length score only using conflicts of the current run
 *
 *  @return the variable's conflict length score only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarConflictlengthScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarConflictlengthScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   downscore = SCIPvarGetAvgConflictlengthCurrentRun(var, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetAvgConflictlengthCurrentRun(var, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's average conflict length
 *
 *  @return the variable's average conflict length
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgConflictlength(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgConflictlength", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPvarGetAvgConflictlength(var, dir);
}

/** returns the variable's average  conflict length only using conflicts of the current run
 *
 *  @return the variable's average conflict length only using conflicts of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgConflictlengthCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgConflictlengthCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPvarGetAvgConflictlengthCurrentRun(var, dir);
}

/** returns the average number of inferences found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of inferences found after branching on the variable in given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgInferences", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPvarGetAvgInferences(var, scip->stat, dir);
}

/** returns the average number of inferences found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of inferences found after branching on the variable in given direction in the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgInferencesCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, dir);
}

/** returns the variable's average inference score value
 *
 *  @return the variable's average inference score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgInferenceScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real inferdown;
   SCIP_Real inferup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgInferenceScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   inferdown = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, inferdown, inferup);
}

/** returns the variable's average inference score value only using inferences of the current run
 *
 *  @return the variable's average inference score value only using inferences of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgInferenceScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real inferdown;
   SCIP_Real inferup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgInferenceScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   inferdown = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, inferdown, inferup);
}

/** initializes the upwards and downwards pseudocosts, conflict scores, conflict lengths, inference scores, cutoff scores
 *  of a variable to the given values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPinitVarBranchStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which should be initialized */
   SCIP_Real             downpscost,         /**< value to which pseudocosts for downwards branching should be initialized */
   SCIP_Real             uppscost,           /**< value to which pseudocosts for upwards branching should be initialized */
   SCIP_Real             downvsids,          /**< value to which VSIDS score for downwards branching should be initialized */
   SCIP_Real             upvsids,            /**< value to which VSIDS score for upwards branching should be initialized */
   SCIP_Real             downconflen,        /**< value to which conflict length score for downwards branching should be initialized */
   SCIP_Real             upconflen,          /**< value to which conflict length score for upwards branching should be initialized */
   SCIP_Real             downinfer,          /**< value to which inference counter for downwards branching should be initialized */
   SCIP_Real             upinfer,            /**< value to which inference counter for upwards branching should be initialized */
   SCIP_Real             downcutoff,         /**< value to which cutoff counter for downwards branching should be initialized */
   SCIP_Real             upcutoff            /**< value to which cutoff counter for upwards branching should be initialized */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPinitVarBranchStats", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(downpscost >= 0.0 && uppscost >= 0.0);
   assert(downvsids >= 0.0 && upvsids >= 0.0);
   assert(downconflen >= 0.0 && upconflen >= 0.0);
   assert(downinfer >= 0.0 && upinfer >= 0.0);
   assert(downcutoff >= 0.0 && upcutoff >= 0.0);

   if( !SCIPisFeasZero(scip, downpscost) || !SCIPisFeasZero(scip, downvsids)
      || !SCIPisFeasZero(scip, downinfer) || !SCIPisFeasZero(scip, downcutoff) )
   {
      SCIP_CALL( SCIPvarIncNBranchings(var, NULL, NULL, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, SCIP_UNKNOWN, 1) );
      SCIP_CALL( SCIPvarUpdatePseudocost(var, scip->set, scip->stat, -1.0, downpscost, 1.0) );
      SCIP_CALL( SCIPvarIncInferenceSum(var,  NULL, NULL, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, SCIP_UNKNOWN, downinfer) );
      SCIP_CALL( SCIPvarIncVSIDS(var, NULL, scip->set, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, SCIP_UNKNOWN, downvsids) );
      SCIP_CALL( SCIPvarIncCutoffSum(var, NULL, NULL, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, SCIP_UNKNOWN, downcutoff) );
   }

   if( !SCIPisFeasZero(scip, downconflen) )
   {
      SCIP_CALL( SCIPvarIncNActiveConflicts(var, NULL, NULL, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, SCIP_UNKNOWN, downconflen) );
   }

   if( !SCIPisFeasZero(scip, uppscost) || !SCIPisFeasZero(scip, upvsids)
      || !SCIPisFeasZero(scip, upinfer) || !SCIPisFeasZero(scip, upcutoff) )
   {
      SCIP_CALL( SCIPvarIncNBranchings(var, NULL, NULL, scip->stat, SCIP_BRANCHDIR_UPWARDS, SCIP_UNKNOWN, 1) );
      SCIP_CALL( SCIPvarUpdatePseudocost(var, scip->set, scip->stat, 1.0, uppscost, 1.0) );
      SCIP_CALL( SCIPvarIncInferenceSum(var, NULL, NULL, scip->stat, SCIP_BRANCHDIR_UPWARDS, SCIP_UNKNOWN, upinfer) );
      SCIP_CALL( SCIPvarIncVSIDS(var, NULL, scip->set, scip->stat, SCIP_BRANCHDIR_UPWARDS, SCIP_UNKNOWN, upvsids) );
      SCIP_CALL( SCIPvarIncCutoffSum(var, NULL, NULL, scip->stat, SCIP_BRANCHDIR_UPWARDS, SCIP_UNKNOWN, upcutoff) );
   }

   if( !SCIPisFeasZero(scip, upconflen) )
   {
      SCIP_CALL( SCIPvarIncNActiveConflicts(var, NULL, NULL, scip->stat, SCIP_BRANCHDIR_UPWARDS, SCIP_UNKNOWN, upconflen) );
   }

   return SCIP_OKAY;
}

/** initializes the upwards and downwards conflict scores, conflict lengths, inference scores, cutoff scores of a
 *  variable w.r.t. a value by the given values (SCIP_VALUEHISTORY)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPinitVarValueBranchStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which should be initialized */
   SCIP_Real             value,              /**< domain value, or SCIP_UNKNOWN */
   SCIP_Real             downvsids,          /**< value to which VSIDS score for downwards branching should be initialized */
   SCIP_Real             upvsids,            /**< value to which VSIDS score for upwards branching should be initialized */
   SCIP_Real             downconflen,        /**< value to which conflict length score for downwards branching should be initialized */
   SCIP_Real             upconflen,          /**< value to which conflict length score for upwards branching should be initialized */
   SCIP_Real             downinfer,          /**< value to which inference counter for downwards branching should be initialized */
   SCIP_Real             upinfer,            /**< value to which inference counter for upwards branching should be initialized */
   SCIP_Real             downcutoff,         /**< value to which cutoff counter for downwards branching should be initialized */
   SCIP_Real             upcutoff            /**< value to which cutoff counter for upwards branching should be initialized */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPinitVarValueBranchStats", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(downvsids >= 0.0 && upvsids >= 0.0);
   assert(downconflen >= 0.0 && upconflen >= 0.0);
   assert(downinfer >= 0.0 && upinfer >= 0.0);
   assert(downcutoff >= 0.0 && upcutoff >= 0.0);

   if( !SCIPisFeasZero(scip, downvsids) || !SCIPisFeasZero(scip, downinfer) || !SCIPisFeasZero(scip, downcutoff) )
   {
      SCIP_CALL( SCIPvarIncNBranchings(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, value, 1) );
      SCIP_CALL( SCIPvarIncInferenceSum(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, value, downinfer) );
      SCIP_CALL( SCIPvarIncVSIDS(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, value, downvsids) );
      SCIP_CALL( SCIPvarIncCutoffSum(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, value, downcutoff) );
   }

   if( !SCIPisFeasZero(scip, downconflen) )
   {
      SCIP_CALL( SCIPvarIncNActiveConflicts(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, value, downconflen) );
   }

   if( !SCIPisFeasZero(scip, upvsids) || !SCIPisFeasZero(scip, upinfer) || !SCIPisFeasZero(scip, upcutoff) )
   {
      SCIP_CALL( SCIPvarIncNBranchings(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_UPWARDS, value, 1) );
      SCIP_CALL( SCIPvarIncInferenceSum(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_UPWARDS, value, upinfer) );
      SCIP_CALL( SCIPvarIncVSIDS(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_UPWARDS, value, upvsids) );
      SCIP_CALL( SCIPvarIncCutoffSum(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_UPWARDS, value, upcutoff) );
   }

   if( !SCIPisFeasZero(scip, upconflen) )
   {
      SCIP_CALL( SCIPvarIncNActiveConflicts(var, SCIPblkmem(scip), scip->set, scip->stat, SCIP_BRANCHDIR_UPWARDS, value, upconflen) );
   }

   return SCIP_OKAY;
}

/** returns the average number of cutoffs found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of cutoffs found after branching on the variable in given direction
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgCutoffs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPvarGetAvgCutoffs(var, scip->stat, dir);
}

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 *
 *  @return the average number of cutoffs found after branching on the variable in given direction in the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgCutoffsCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   return SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, dir);
}

/** returns the variable's average cutoff score value
 *
 *  @return the variable's average cutoff score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgCutoffScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   cutoffdown = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, cutoffdown, cutoffup);
}

/** returns the variable's average cutoff score value, only using cutoffs of the current run
 *
 *  @return the variable's average cutoff score value, only using cutoffs of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   cutoffdown = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, cutoffdown, cutoffup);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor
 *
 *  @return the variable's average inference/cutoff score value
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgInferenceCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   SCIP_Real avginferdown;
   SCIP_Real avginferup;
   SCIP_Real avginfer;
   SCIP_Real inferdown;
   SCIP_Real inferup;
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgInferenceCutoffScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   avginferdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   avginferup = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);
   avginfer = (avginferdown + avginferup)/2.0;
   inferdown = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);
   cutoffdown = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var,
      inferdown + cutoffweight * avginfer * cutoffdown, inferup + cutoffweight * avginfer * cutoffup);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor, only using inferences and cutoffs of the current run
 *
 *  @return the variable's average inference/cutoff score value, only using inferences and cutoffs of the current run
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetVarAvgInferenceCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   SCIP_Real avginferdown;
   SCIP_Real avginferup;
   SCIP_Real avginfer;
   SCIP_Real inferdown;
   SCIP_Real inferup;
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetVarAvgInferenceCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   avginferdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   avginferup = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);
   avginfer = (avginferdown + avginferup)/2.0;
   inferdown = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);
   cutoffdown = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var,
      inferdown + cutoffweight * avginfer * cutoffdown, inferup + cutoffweight * avginfer * cutoffup);
}

/** outputs variable information to file stream via the message system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed
 */
SCIP_RETCODE SCIPprintVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPvarPrint(var, scip->set, scip->messagehdlr, file) );

   return SCIP_OKAY;
}
