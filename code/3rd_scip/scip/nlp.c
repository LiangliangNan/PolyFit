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

/**@file   nlp.c
 * @brief  NLP management methods and datastructures
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 *
 *  In NLP management, we have to differ between the current NLP and the NLPI problem
 *  stored in the NLP solver. All NLP methods affect the current NLP only.
 *  Before solving the current NLP with the NLP solver, the NLP solvers data
 *  has to be updated to the current NLP with a call to nlpFlush().
 *
 *  @todo handle linear rows from LP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/intervalarith.h"
#include "scip/clock.h"
#include "scip/nlp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/event.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "nlpi/nlpi.h"
#include "nlpi/pub_expr.h"
#include "nlpi/struct_expr.h"
#include "scip/struct_nlp.h"
/* to get value of parameter "nlp/solver" and nlpis array and to get access to set->lp for releasing a variable */
#include "scip/struct_set.h"
/* to get nlp, set, ... in event handling */
#include "scip/struct_scip.h"

/* defines */

#define EVENTHDLR_NAME   "nlpEventHdlr"      /**< name of NLP event handler that catches variable events */
#define EVENTHDLR_DESC   "handles all events necessary for maintaining NLP data"  /**< description of NLP event handler */
#define ADDNAMESTONLPI   0                   /**< whether to give variable and row names to NLPI */

#ifdef __cplusplus
extern "C" {
#endif

/* avoid inclusion of scip.h */
BMS_BLKMEM* SCIPblkmem(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

/*
 * forward declarations
 */

/** NLP event handler execution method */
static
SCIP_DECL_EVENTEXEC( eventExecNlp );

/** announces, that a row of the NLP was modified
 * adjusts status of current solution
 * calling method has to ensure that change is passed to the NLPI!
 */
static
SCIP_RETCODE nlpRowChanged(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row which was changed */
   );

/*
 * public expression tree methods
 */

/** returns variables of expression tree */
SCIP_VAR** SCIPexprtreeGetVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(tree != NULL);

   return (SCIP_VAR**)tree->vars;
}

/** stores array of variables in expression tree */
SCIP_RETCODE SCIPexprtreeSetVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
   )
{
   assert(tree != NULL);
   assert(vars != NULL || nvars == 0);

   if( nvars == 0 )
   {
      BMSfreeBlockMemoryArrayNull(tree->blkmem, &tree->vars, tree->nvars);
      tree->nvars = 0;
   }
   else if( tree->vars != NULL )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars, nvars) );
      BMScopyMemoryArray(tree->vars, (void**)vars, nvars);
   }
   else
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(tree->blkmem, &tree->vars, (void**)vars, nvars) );
   }

   tree->nvars = nvars;

   assert(tree->vars != NULL || tree->nvars == 0);

   return SCIP_OKAY;
}

/** adds variables to the expression tree variables array */
SCIP_RETCODE SCIPexprtreeAddVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
   )
{
   assert(tree != NULL);
   assert(vars != NULL || nvars == 0);
   assert(tree->vars != NULL || tree->nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   if( tree->nvars == 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(tree->blkmem, &tree->vars, (void**)vars, nvars) );
      tree->nvars = nvars;
      return SCIP_OKAY;
   }

   /* append vars to tree->vars array */
   SCIP_ALLOC( BMSreallocBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars, tree->nvars + nvars) );
   BMScopyMemoryArray(&tree->vars[tree->nvars], (void**)vars, nvars);  /*lint !e866*/
   tree->nvars += nvars;

   return SCIP_OKAY;
}

/** prints an expression tree using variable names from variables array */
SCIP_RETCODE SCIPexprtreePrintWithNames(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file for printing, or NULL for stdout */
   )
{
   const char** varnames;
   int i;

   assert(tree != NULL);

   if( tree->nvars == 0 )
   {
      SCIPexprtreePrint(tree, messagehdlr, file, NULL, NULL);
      return SCIP_OKAY;
   }

   assert(tree->vars != NULL);

   SCIP_ALLOC( BMSallocMemoryArray(&varnames, tree->nvars) );
   for( i = 0; i < tree->nvars; ++i )
      varnames[i] = SCIPvarGetName((SCIP_VAR*)tree->vars[i]);

   SCIPexprtreePrint(tree, messagehdlr, file, varnames, NULL);

   BMSfreeMemoryArray(&varnames);

   return SCIP_OKAY;
}

/** searches the variables array of an expression tree for a variable and returns its position, or -1 if not found
 * Note that this is an O(n) operation!
 */
int SCIPexprtreeFindVar(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_VAR*             var                 /**< variable to search for */
   )
{
   int i;

   assert(tree != NULL);
   assert(var  != NULL);

   for( i = 0; i < tree->nvars; ++i )
      if( (SCIP_VAR*)tree->vars[i] == var )
         return i;

   return -1;
}

/** removes fixed variables from an expression tree, so that at exit all variables are active */
SCIP_RETCODE SCIPexprtreeRemoveFixedVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool*            changed,            /**< buffer to store whether the tree was changed, i.e., whether there was a fixed variable */
   int*                  varpos,             /**< array of length at least tree->nvars to store new indices of previously existing variables in expression tree, or -1 if variable was removed; set to NULL if not of interest */
   int*                  newvarsstart        /**< buffer to store index in tree->vars array where new variables begin, or NULL if not of interest */
   )
{
   SCIP_HASHMAP* varhash;
   int i;
   int j;
   int nvarsold;
   SCIP_VAR* var;
   SCIP_Real scalar;
   SCIP_Real constant;
   SCIP_EXPR** replaceexprs;
   SCIP_Bool havefixedvar;
   int idx;
   int* newpos;
   int offset;

   assert(tree != NULL);
   assert(tree->vars != NULL || tree->nvars == 0);
   assert(changed != NULL);

   *changed = FALSE;
   if( newvarsstart != NULL )
      *newvarsstart = tree->nvars;

   if( tree->nvars == 0 )
      return SCIP_OKAY;

   /* create hash map from variable to indices in tree->vars and check if there is a non-fixed variable */
   havefixedvar = FALSE;
   SCIP_CALL( SCIPhashmapCreate(&varhash, tree->blkmem, tree->nvars) );
   for( i = 0; i < tree->nvars; ++i )
   {
      /* it's not possible to add a variable twice to the varhash map */
      if( SCIPhashmapExists(varhash, tree->vars[i]) )
         continue;

      SCIP_CALL( SCIPhashmapInsert(varhash, tree->vars[i], (void*)(size_t)i) );

      if( !SCIPvarIsActive((SCIP_VAR*)tree->vars[i]) )
         havefixedvar = TRUE;
   }

   if( !havefixedvar )
   {
      /* nothing to do */
      if( varpos != NULL )
         for( i = 0; i < tree->nvars; ++i )
            varpos[i] = i;
      SCIPhashmapFree(&varhash);
      return SCIP_OKAY;
   }

   /* we will do something */
   *changed = TRUE;

   nvarsold = tree->nvars;

   /* array to store expressions that replace a variable expression in the tree */
   SCIP_ALLOC( BMSallocBlockMemoryArray(tree->blkmem, &replaceexprs, nvarsold) );
   BMSclearMemoryArray(replaceexprs, nvarsold);

   /* construct for each nonactive variable an expression that replaces this variable in the tree */
   for( i = 0; i < nvarsold; ++i )
   {
      var = (SCIP_VAR*)tree->vars[i];

      if( SCIPvarIsActive(var) )
         continue;

      scalar   = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &scalar, &constant) );

      if( scalar == 0.0 )
      {
         /* variable is fixed, thus replace by constant expression in tree */
         SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_CONST, constant) );
         continue;
      }

      if( SCIPvarIsActive(var) )
      {
         /* variable was aggregated or negated, thus replace by scalar * var + constant */
         if( !SCIPhashmapExists(varhash, var) )
         {
            /* var not in tree yet, so add it */
            SCIP_CALL( SCIPexprtreeAddVars(tree, 1, &var) );
            idx = tree->nvars - 1;
            SCIP_CALL( SCIPhashmapInsert(varhash, (void*)var, (void*)(size_t)idx) );
         }
         else
         {
            idx = (int)(size_t) SCIPhashmapGetImage(varhash, (void*)var);
         }
         assert(idx >= 0 && idx < tree->nvars);
         assert((SCIP_VAR*)tree->vars[idx] == var);

         SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_VARIDX, idx) );
         if( scalar != 1.0 || constant != 0.0 )
         {
            /* multiply by scalar and add constant -> linear expression */
            SCIP_CALL( SCIPexprCreateLinear(tree->blkmem, &replaceexprs[i], 1, &replaceexprs[i], &scalar, constant) );
         }
         continue;
      }

      {
         SCIP_EXPR** children;
         SCIP_Real*  coefs;
         int         nchildren;
         SCIP_VAR*   mvar;
         SCIP_Real   mscalar;

         /* var is now multi-aggregated, thus replace by scalar * (multaggrconst + sum_j multaggrscalar_j*multaggrvar_j) + constant */
         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR );

         /* allocate array for children and coefficients */
         SCIP_ALLOC( BMSallocBlockMemoryArray(tree->blkmem, &children, SCIPvarGetMultaggrNVars(var)) );  /*lint !e666 */
         SCIP_ALLOC( BMSallocBlockMemoryArray(tree->blkmem, &coefs,    SCIPvarGetMultaggrNVars(var)) );  /*lint !e666 */
         nchildren = 0;

         /* linear part
          * turn each variable in SCIPvarGetMultaggrVars(var) into an active or multi-aggregated one and add corresponding term to summands */
         for( j = 0; j < SCIPvarGetMultaggrNVars(var); ++j )
         {
            mvar      = SCIPvarGetMultaggrVars(var)[j];
            mscalar   = scalar * SCIPvarGetMultaggrScalars(var)[j];
            SCIP_CALL( SCIPvarGetProbvarSum(&mvar, set, &mscalar, &constant) );

            /* if variable mvar is fixed, constant has been added to constant and we can continue */
            if( mscalar == 0.0 )
               continue;

            assert(SCIPvarIsActive(mvar) || SCIPvarGetStatus(mvar) == SCIP_VARSTATUS_MULTAGGR);

            /* add mvar to tree, if not in tree yet */
            if( !SCIPhashmapExists(varhash, mvar) )
            {
               /* var not in tree yet, so add it */
               SCIP_CALL( SCIPexprtreeAddVars(tree, 1, &mvar) );
               idx = tree->nvars - 1;
               SCIP_CALL( SCIPhashmapInsert(varhash, (void*)mvar, (void*)(size_t)idx) );
            }
            else
            {
               idx = (int)(size_t) SCIPhashmapGetImage(varhash, (void*)mvar);
            }
            assert(idx >= 0 && idx < tree->nvars);
            assert((SCIP_VAR*)tree->vars[idx] == mvar);

            SCIP_CALL( SCIPexprCreate(tree->blkmem, &children[nchildren], SCIP_EXPR_VARIDX, idx) );
            coefs[nchildren] = mscalar;
            ++nchildren;
         }

         /* constant part */
         constant += scalar * SCIPvarGetMultaggrConstant(var);

         if( nchildren == 0 )
         {
            /* somehow all aggregation variables were fixed */
            SCIP_CALL( SCIPexprCreate(tree->blkmem, &replaceexprs[i], SCIP_EXPR_CONST, constant) );
         }
         else if( nchildren == 1 && constant == 0.0 )
         {
            /* somehow everything collapsed to one summand -> use that one for replaceexprs[i]*/
            replaceexprs[i] = children[0];
         }
         else
         {
            /* set replaceexprs[i] to linear expression in children */
            SCIP_CALL( SCIPexprCreateLinear(tree->blkmem, &replaceexprs[i], nchildren, children, coefs, constant) );
         }

         BMSfreeBlockMemoryArray(tree->blkmem, &children, SCIPvarGetMultaggrNVars(var));
         BMSfreeBlockMemoryArray(tree->blkmem, &coefs,    SCIPvarGetMultaggrNVars(var));
      }
   }

   /* replace variables in tree by assembled expressions */
   SCIP_CALL( SCIPexprtreeSubstituteVars(tree, replaceexprs) );
   /* free replaceexprs */
   for( i = 0; i < nvarsold; ++i )
      if( replaceexprs[i] != NULL )
         SCIPexprFreeDeep(tree->blkmem, &replaceexprs[i]);
   BMSfreeBlockMemoryArray(tree->blkmem, &replaceexprs, nvarsold);

   /* the varhash is not needed anymore */
   SCIPhashmapFree(&varhash);

   /* remove inactive variables from vars array and recompute variable indices */
   SCIP_ALLOC( BMSallocBlockMemoryArray(tree->blkmem, &newpos, tree->nvars) );
   offset = 0;
   for( i = 0; i < tree->nvars; ++i )
   {
      if( SCIPvarIsActive((SCIP_VAR*)tree->vars[i]) || i >= nvarsold )
      {
         /* a new variable need to be either active or multi-aggregated */
         assert(i < nvarsold || SCIPvarIsActive((SCIP_VAR*)tree->vars[i]) || SCIPvarGetStatus((SCIP_VAR*)tree->vars[i]) == SCIP_VARSTATUS_MULTAGGR);
         newpos[i] = i - offset;
      }
      else
      {
         /* non-active variable are removed */
         newpos[i] = -1;
         ++offset;
      }
      if( varpos != NULL && i < nvarsold )
         varpos[i] = newpos[i];
   }
   if( newvarsstart != NULL )
      *newvarsstart -= offset;

   /* update indices in tree */
   SCIPexprReindexVars(tree->root, newpos);

   /* move variable in expression tree vars array
    * check if there is a fixed variable left */
   havefixedvar = FALSE;
   for( i = 0; i < tree->nvars; ++i )
   {
      if( newpos[i] == -1 )
      {
         /* variable was removed */
         assert(!SCIPvarIsActive((SCIP_VAR*)tree->vars[i]));
         continue;
      }
      /* variable is moved */
      tree->vars[newpos[i]] = tree->vars[i];
      if( !SCIPvarIsActive((SCIP_VAR*)tree->vars[i]) )
         havefixedvar = TRUE;
   }

   /* free newpos array; resize vars array */
   BMSfreeBlockMemoryArray(tree->blkmem, &newpos, tree->nvars);
   if( offset < tree->nvars )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars, tree->nvars - offset) );
      tree->nvars -= offset;
   }
   else
   {
      /* all variables were removed */
      BMSfreeBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars);
      tree->nvars = 0;
   }

   if( havefixedvar )
   {
      /* if there are still fixed variables left, then this are newly added multi-aggregated variables
       * it is then save to call this function recursively, since the original active variables should not be moved,
       * i.e., varpos and *newvarsstart will remain valid
       */
      SCIP_Bool gotchange;

      SCIP_CALL( SCIPexprtreeRemoveFixedVars(tree, set, &gotchange, NULL, NULL) );
      assert(gotchange);
   }

   return SCIP_OKAY;
}

/*
 * private NLP nonlinear row methods
 */

/** announces, that the given linear coefficient in the constraint matrix changed */
static
SCIP_RETCODE nlrowLinearCoefChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var,                /**< variable which coefficient changed */
   SCIP_Real             coef,               /**< new coefficient of variable, 0.0 if deleted */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);
   assert(var   != NULL);

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      /* update NLPI problem, if row is in NLPI already */
      if( nlrow->nlpiindex >= 0 )
      {
         int idx;

         /* get index of variable in NLPI */
         assert(SCIPhashmapExists(nlp->varhash, var));
         idx = (int)(size_t)SCIPhashmapGetImage(nlp->varhash, var);
         assert(idx >= 0 && idx < nlp->nvars);

         idx = nlp->varmap_nlp2nlpi[idx];
         assert(idx >= 0 && idx < nlp->nvars_solver);

         /* change coefficient in NLPI problem */
         SCIP_CALL( SCIPnlpiChgLinearCoefs(nlp->solver, nlp->problem, nlrow->nlpiindex, 1, &idx, &coef) );
      }
   }

   return SCIP_OKAY;
}

/** announces, that an element in the quadratic part of a nonlinear row changed */
static
SCIP_RETCODE nlrowQuadElemChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_QUADELEM         quadelem,           /**< new element (variable indices and new values), quadelem.coef == 0 if it was deleted */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);
   assert(quadelem.idx1 >= 0);
   assert(quadelem.idx1 < nlrow->nquadvars);
   assert(quadelem.idx2 >= 0);
   assert(quadelem.idx2 < nlrow->nquadvars);

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      /* update NLPI problem, if row is in NLPI already */
      if( nlrow->nlpiindex >= 0 )
      {
         SCIP_QUADELEM elem;

         /* get NLPI index of first variable */
         assert(nlrow->quadvars[quadelem.idx1] != NULL);
         assert(SCIPhashmapExists(nlp->varhash, nlrow->quadvars[quadelem.idx1]));
         elem.idx1 = (int)(size_t)SCIPhashmapGetImage(nlp->varhash, nlrow->quadvars[quadelem.idx1]);
         assert(elem.idx1 >= 0 && elem.idx1 < nlp->nvars);

         elem.idx1 = nlp->varmap_nlp2nlpi[elem.idx1];
         assert(elem.idx1 >= 0 && elem.idx1 < nlp->nvars_solver);

         /* get NLPI index of second variable */
         assert(nlrow->quadvars[quadelem.idx2] != NULL);
         assert(SCIPhashmapExists(nlp->varhash, nlrow->quadvars[quadelem.idx2]));
         elem.idx2 = (int)(size_t)SCIPhashmapGetImage(nlp->varhash, nlrow->quadvars[quadelem.idx2]);
         assert(elem.idx2 >= 0 && elem.idx2 < nlp->nvars);

         elem.idx2 = nlp->varmap_nlp2nlpi[elem.idx2];
         assert(elem.idx2 >= 0 && elem.idx2 < nlp->nvars_solver);

         /* make sure idx1 <= idx2 */
         if( elem.idx1 > elem.idx2 )
         {
            int tmp;
            tmp = elem.idx2;
            elem.idx2 = elem.idx1;
            elem.idx1 = tmp;
         }

         elem.coef = quadelem.coef;

         /* change coefficient in NLPI problem */
         SCIP_CALL( SCIPnlpiChgQuadCoefs(nlp->solver, nlp->problem, nlrow->nlpiindex, 1, &elem) );
      }
   }

   return SCIP_OKAY;
}

/** announces, that an expression tree changed */
static
SCIP_RETCODE nlrowExprtreeChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      if( nlrow->nlpiindex >= 0 )
      {
         /* change expression tree in NLPI problem */
         int* nlinidxs;

         /* get indices of variables in expression tree part of row */
         if( nlrow->exprtree != NULL )
         {
            int i;
            int n;
            SCIP_VAR* var;

            n = SCIPexprtreeGetNVars(nlrow->exprtree);
            assert(n == 0 || SCIPexprtreeGetVars(nlrow->exprtree) != NULL);

            SCIP_CALL( SCIPsetAllocBufferArray(set, &nlinidxs, n) );

            for( i = 0; i < n; ++i )
            {
               var = SCIPexprtreeGetVars(nlrow->exprtree)[i];
               assert(var != NULL);
               assert(SCIPvarIsActive(var)); /* at this point, there should be only active variables in the row */

               assert(SCIPhashmapExists(nlp->varhash, var));
               nlinidxs[i] = nlp->varmap_nlp2nlpi[(size_t) (void*) SCIPhashmapGetImage(nlp->varhash, var)];
            }

            SCIP_CALL( SCIPnlpiChgExprtree(nlp->solver, nlp->problem, nlrow->nlpiindex, nlinidxs, nlrow->exprtree) );

            SCIPsetFreeBufferArray(set, &nlinidxs);
         }
         else
         {
            SCIP_CALL( SCIPnlpiChgExprtree(nlp->solver, nlp->problem, nlrow->nlpiindex, NULL, NULL) );
         }
      }
   }

   return SCIP_OKAY;
}

/** announces, that a parameter in an expression tree has changed */
static
SCIP_RETCODE nlrowExprtreeParamChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   paramidx,           /**< index of parameter which has changed, or -1 if all changed */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->exprtree != NULL);
   assert(paramidx >= -1);
   assert(paramidx <  SCIPexprtreeGetNParams(nlrow->exprtree));

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      if( nlrow->nlpiindex >= 0 )
      {
         if( paramidx >= 0 )
         {
            /* change coefficient in NLPI problem */
            SCIP_CALL( SCIPnlpiChgNonlinCoef(nlp->solver, nlp->problem, nlrow->nlpiindex, paramidx, SCIPexprtreeGetParamVals(nlrow->exprtree)[paramidx]) );
         }
         else
         {
            SCIP_Real* paramvals;
            int i;
            int n;

            /* change all coefficients in NLPI problem */
            n = SCIPexprtreeGetNParams(nlrow->exprtree);
            paramvals = SCIPexprtreeGetParamVals(nlrow->exprtree);
            for( i = 0; i < n; ++i )
            {
               SCIP_CALL( SCIPnlpiChgNonlinCoef(nlp->solver, nlp->problem, nlrow->nlpiindex, i, paramvals[i]) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** notifies nonlinear row, that its sides were changed */
static
SCIP_RETCODE nlrowSideChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      if( nlrow->nlpiindex >= 0 )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

         /* change sides in NLPI problem */
         lhs = nlrow->lhs;
         rhs = nlrow->rhs;
         if( !SCIPsetIsInfinity(set, -lhs) )
            lhs -= nlrow->constant;
         if( !SCIPsetIsInfinity(set,  rhs) )
            rhs -= nlrow->constant;

         SCIP_CALL( SCIPnlpiChgConsSides(nlp->solver, nlp->problem, 1, &nlrow->nlpiindex, &lhs, &rhs) );
      }
   }

   return SCIP_OKAY;
}

/** notifies nonlinear row, that its constant was changed */
static
SCIP_RETCODE nlrowConstantChanged(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlrow != NULL);

   nlrow->activity = SCIP_INVALID;
   nlrow->validactivitynlp = -1;
   nlrow->pseudoactivity = SCIP_INVALID;
   nlrow->validpsactivitydomchg = -1;
   nlrow->minactivity = SCIP_INVALID;
   nlrow->maxactivity = SCIP_INVALID;
   nlrow->validactivitybdsdomchg = -1;

   if( nlrow->nlpindex >= 0 )
   {
      assert(nlp != NULL);

      /* notify NLP that row has changed */
      SCIP_CALL( nlpRowChanged(nlp, set, stat, nlrow) );

      if( nlrow->nlpiindex >= 0 )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

         lhs = nlrow->lhs;
         rhs = nlrow->rhs;
         if( !SCIPsetIsInfinity(set, -lhs) )
            lhs -= nlrow->constant;
         if( !SCIPsetIsInfinity(set,  rhs) )
            rhs -= nlrow->constant;

         /* change sides in NLPI problem */
         SCIP_CALL( SCIPnlpiChgConsSides(nlp->solver, nlp->problem, 1, &nlrow->nlpiindex, &lhs, &rhs) );
      }
   }

   return SCIP_OKAY;
}

/** sorts linear part of row entries such that lower variable indices precede higher ones */
static
void nlrowSortLinear(
   SCIP_NLROW*           nlrow               /**< nonlinear row to be sorted */
   )
{
   assert(nlrow != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( nlrow->linvarssorted )
      return;

   /* sort linear coefficients */
   SCIPsortPtrReal((void**)nlrow->linvars, nlrow->lincoefs, SCIPvarComp, nlrow->nlinvars);

   nlrow->linvarssorted = TRUE;
}

/** searches linear variable in nonlinear row, returns position in linvars vector or -1 if not found */
static
int nlrowSearchLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be searched in */
   SCIP_VAR*             var                 /**< variable to be searched for */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);

   if( nlrow->nlinvars == 0 )
      return -1;

   nlrowSortLinear(nlrow);
   if( !SCIPsortedvecFindPtr((void**)nlrow->linvars, SCIPvarComp, (void*)var, nlrow->nlinvars, &pos) )
      return -1;

   return pos;
}

/** moves a coefficient in a nonlinear row to a different place, and updates all corresponding data structures */
static
void nlrowMoveLinearCoef(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= oldpos && oldpos < nlrow->nlinvars);
   assert(0 <= newpos && newpos < nlrow->nlinvars);
   assert(nlrow->linvars[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   nlrow->linvars[newpos]  = nlrow->linvars[oldpos];
   nlrow->lincoefs[newpos] = nlrow->lincoefs[oldpos];

   /* update sorted flags */
   nlrow->linvarssorted = FALSE;
}

/** adds a previously non existing linear coefficient to a nonlinear row */
static
SCIP_RETCODE nlrowAddLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value of coefficient */
   )
{
   int pos;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(var    != NULL);
   assert(!SCIPsetIsZero(set, coef));

   /* assert that only active variables are added once the row is in the NLP */
   assert(nlrow->nlpindex == -1 || SCIPvarIsActive(var) );

   SCIP_CALL( SCIPnlrowEnsureLinearSize(nlrow, blkmem, set, nlrow->nlinvars+1) );
   assert(nlrow->linvars  != NULL);
   assert(nlrow->lincoefs != NULL);

   pos = nlrow->nlinvars;
   nlrow->nlinvars++;

   /* insert the variable */
   nlrow->linvars [pos] = var;
   nlrow->lincoefs[pos] = coef;

   SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, var, coef, nlp) );

   /* update sorted flag */
   if( pos > 0 && SCIPvarCompare(nlrow->linvars[pos-1], nlrow->linvars[pos]) > 0 )
      nlrow->linvarssorted = FALSE;

   SCIPsetDebugMsg(set, "added linear coefficient %g * <%s> at position %d to nonlinear row <%s>\n",
      coef, SCIPvarGetName(var), pos, nlrow->name);

   return SCIP_OKAY;
}

/** adds a linear coefficient to a nonlinear row
 * if the variable exists in the linear part of the row already, the coefficients are added
 * otherwise the variable is added to the row */
static
SCIP_RETCODE nlrowAddToLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef,               /**< value of coefficient */
   SCIP_Bool             removefixed         /**< whether to disaggregate var before adding */
   )
{
   int pos;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(var    != NULL);

   if( removefixed && !SCIPvarIsActive(var) )
   {
      SCIP_Real constant;

      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &coef, &constant) );
      if( constant != 0.0 )
      {
         nlrow->constant += constant;
         SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
      }

      if( !SCIPvarIsActive(var) )
      {
         int j;

         /* if var is still not active, then it is multi-aggregated */
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

         if( SCIPvarGetMultaggrConstant(var) != 0.0 )
         {
            nlrow->constant += coef * SCIPvarGetMultaggrConstant(var);
            SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
         }

         for( j = 0; j < SCIPvarGetMultaggrNVars(var); ++j )
         {
            SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var)[j], SCIPvarGetMultaggrScalars(var)[j] * coef, TRUE) );
         }

         return SCIP_OKAY;
      }
   }
   assert(!removefixed || SCIPvarIsActive(var));

   if( SCIPsetIsZero(set, coef) )
      return SCIP_OKAY;

   pos = nlrowSearchLinearCoef(nlrow, var);

   if( pos == -1 )
   {
      /* add as new coefficient */
      SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, var, coef) );
   }
   else
   {
      assert(pos >= 0);
      assert(pos <  nlrow->nlinvars);
      assert(nlrow->linvars[pos] == var);

      /* add to previously existing coefficient */
      nlrow->lincoefs[pos] += coef;
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
SCIP_RETCODE nlrowDelLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_VAR* var;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] != NULL);

   var = nlrow->linvars[pos];

   /* move last coefficient to position of empty slot (should set sorted flag to FALSE, if not last variable was deleted) */
   nlrowMoveLinearCoef(nlrow, nlrow->nlinvars-1, pos);
   nlrow->nlinvars--;
   assert(pos == nlrow->nlinvars || nlrow->linvarssorted == FALSE);

   SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, var, 0.0, nlp) );

   return SCIP_OKAY;
}

/** changes a coefficient at given position of a nonlinear row */
static
SCIP_RETCODE nlrowChgLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos,                /**< position in row vector to change */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] != NULL);

   if( SCIPsetIsZero(set, coef) )
   {
      /* delete existing coefficient */
      SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, stat, nlp, pos) );
   }
   else if( !SCIPsetIsEQ(set, nlrow->lincoefs[pos], coef) )
   {
      /* change existing coefficient */
      nlrow->lincoefs[pos] = coef;
      SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, nlrow->linvars[pos], coef, nlp) );
   }

   return SCIP_OKAY;
}

/** sets up the variable hash for quadratic variables, if the number of variables exceeds some given threshold */
static
SCIP_RETCODE nlrowSetupQuadVarsHash(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;
   assert(blkmem != NULL);
   assert(nlrow  != NULL);
   assert(nlrow->quadvarshash == NULL);

   if( nlrow->nquadvars < 3 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapCreate(&nlrow->quadvarshash, blkmem, nlrow->nquadvars) );
   assert(nlrow->quadvarshash != NULL);

   for( i = 0; i < nlrow->nquadvars; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(nlrow->quadvarshash, (void*)nlrow->quadvars[i], (void*)(size_t)i) );
   }

   return SCIP_OKAY;
}

/** sorts quadratic part of row entries */
static
void nlrowSortQuadElem(
   SCIP_NLROW*           nlrow               /**< nonlinear row to be sorted */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->quadelems != NULL);

   /* check, if row is already sorted in the LP part, or if the sorting should be delayed */
   if( nlrow->quadelemssorted )
      return;

   /* sort quadratic elements */
   SCIPquadelemSort(nlrow->quadelems, nlrow->nquadelems);

   nlrow->quadelemssorted = TRUE;
}

/** searches quadratic elements in nonlinear row, returns position of given index pair in quadelems array or -1 if not found */
static
int nlrowSearchQuadElem(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be searched in */
   int                   idx1,               /**< index of first  variable to be searched for */
   int                   idx2                /**< index of second variable to be searched for */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(idx1 >= 0);
   assert(idx1 <  nlrow->nquadvars);
   assert(idx2 >= 0);
   assert(idx2 <  nlrow->nquadvars);

   nlrowSortQuadElem(nlrow);
   if( !SCIPquadelemSortedFind(nlrow->quadelems, idx1, idx2, nlrow->nquadelems, &pos) )
      pos = -1;

   return pos;
}

/** moves a quadratic element in a nonlinear row to a different place, and updates all corresponding data structures */
static
void nlrowMoveQuadElement(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   oldpos,             /**< old position of coefficient */
   int                   newpos              /**< new position of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= oldpos && oldpos < nlrow->nquadelems);
   assert(0 <= newpos && newpos < nlrow->nquadelems);

   if( oldpos == newpos )
      return;

   nlrow->quadelems[newpos] = nlrow->quadelems[oldpos];

   /* update sorted flags */
   nlrow->quadelemssorted = FALSE;
}

/** adds a previously non existing quadratic element to a nonlinear row */
static
SCIP_RETCODE nlrowAddQuadElement(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_QUADELEM         elem                /**< quadratic element to add */
   )
{
   int pos;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(elem.idx1 >= 0);
   assert(elem.idx1 <  nlrow->nquadvars);
   assert(elem.idx2 >= 0);
   assert(elem.idx2 <  nlrow->nquadvars);

   if( SCIPsetIsZero(set, elem.coef) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPnlrowEnsureQuadElementsSize(nlrow, blkmem, set, nlrow->nquadelems+1) );
   assert(nlrow->quadelems != NULL);

   pos = nlrow->nquadelems;
   nlrow->nquadelems++;

   /* insert the element */
   nlrow->quadelems[pos] = elem;

   /* notify row and NLP */
   SCIP_CALL( nlrowQuadElemChanged(nlrow, set, stat, elem, nlp) );

   /* update sorted flag */
   if( pos > 0 )
      nlrow->quadelemssorted = FALSE;

   SCIPsetDebugMsg(set, "added quadratic element %g * <%s> * <%s> at position %d to nonlinear row <%s>\n",
      elem.coef, SCIPvarGetName(nlrow->quadvars[elem.idx1]), SCIPvarGetName(nlrow->quadvars[elem.idx2]), pos, nlrow->name);

   return SCIP_OKAY;
}

/** deletes coefficient at given position from row */
static
SCIP_RETCODE nlrowDelQuadElemPos(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos                 /**< position in row vector to delete */
   )
{
   SCIP_QUADELEM elem;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(0 <= pos && pos < nlrow->nquadelems);

   SCIPsetDebugMsg(set, "delete quad element (%d,%d) at pos %d\n", nlrow->quadelems[pos].idx1, nlrow->quadelems[pos].idx2, pos);

   elem = nlrow->quadelems[pos];

   /* move last coefficient to position of empty slot (should set sorted flag to FALSE, if not last element was deleted) */
   nlrowMoveQuadElement(nlrow, nlrow->nquadelems-1, pos);
   nlrow->nquadelems--;
   assert(pos == nlrow->nquadelems || nlrow->quadelemssorted == FALSE);

   /* notify row and NLP */
   elem.coef = 0.0;
   SCIP_CALL( nlrowQuadElemChanged(nlrow, set, stat, elem, nlp) );

   return SCIP_OKAY;
}

/** changes a coefficient at given position of quadratic element in nonlinear row */
static
SCIP_RETCODE nlrowChgQuadElemPos(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos,                /**< position in quadratic elements array to change */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   assert(nlrow != NULL);
   assert(0 <= pos && pos < nlrow->nquadelems);

   SCIPsetDebugMsg(set, "change quad element (%d,%d) at pos %d to %g\n", nlrow->quadelems[pos].idx1, nlrow->quadelems[pos].idx2, pos, coef);

   if( SCIPsetIsZero(set, coef) )
   {
      /* delete existing coefficient */
      SCIP_CALL( nlrowDelQuadElemPos(nlrow, set, stat, nlp, pos) );
   }
   else if( !SCIPsetIsEQ(set, nlrow->quadelems[pos].coef, coef) )
   {
      /* change existing coefficient */
      nlrow->quadelems[pos].coef = coef;
      SCIP_CALL( nlrowQuadElemChanged(nlrow, set, stat, nlrow->quadelems[pos], nlp) );
   }

   return SCIP_OKAY;
}

/** calculates minimal and maximal activity of row w.r.t. the variable's bounds */
static
SCIP_RETCODE nlrowCalcActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   SCIP_Real inf;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL bounds;
   int i;

   assert(nlrow != NULL);
   assert(set   != NULL);
   assert(stat  != NULL);

   inf = SCIPsetInfinity(set);

   /* calculate activity bounds */
   SCIPintervalSet(&activity, nlrow->constant);
   for( i = 0; i < nlrow->nlinvars && !SCIPintervalIsEntire(inf, activity); ++i )
   {
      SCIPintervalSetBounds(&bounds, SCIPvarGetLbLocal(nlrow->linvars[i]), SCIPvarGetUbLocal(nlrow->linvars[i]));
      SCIPintervalMulScalar(inf, &bounds, bounds, nlrow->lincoefs[i]);
      SCIPintervalAdd(inf, &activity, activity, bounds);
   }

   /* @todo make sure quadelems is sorted */
   for( i = 0; i < nlrow->nquadelems && !SCIPintervalIsEntire(inf, activity); )
   {
      SCIP_Real a;
      SCIP_INTERVAL b, tmp;
      int idx1;

      idx1 = nlrow->quadelems[i].idx1;
      SCIPintervalSetBounds(&bounds, SCIPvarGetLbLocal(nlrow->quadvars[idx1]), SCIPvarGetUbLocal(nlrow->quadvars[idx1]));

      /* for x_i*(a*x_i + sum_j b_jx_j) we assemble a and sum_j b_jx_j */
      a = 0.0;
      SCIPintervalSet(&b, 0.0);
      do
      {
         if( nlrow->quadelems[i].idx1 == nlrow->quadelems[i].idx2 )
         {
            a = nlrow->quadelems[i].coef;
         }
         else
         {
            SCIPintervalSetBounds(&tmp, SCIPvarGetLbLocal(nlrow->quadvars[nlrow->quadelems[i].idx2]), SCIPvarGetUbLocal(nlrow->quadvars[nlrow->quadelems[i].idx2]));
            SCIPintervalMulScalar(inf, &tmp, tmp, nlrow->quadelems[i].coef);
            SCIPintervalAdd(inf, &b, b, tmp);
         }
         ++i;
      }
      while( i < nlrow->nquadvars && idx1 == nlrow->quadelems[i].idx1 );

      /* compute bounds for a*x_i^2 + b*x_i and add to activity bounds */
      SCIPintervalQuad(inf, &bounds, a, b, bounds);
      SCIPintervalAdd(inf, &activity, activity, bounds);
   }

   if( nlrow->exprtree != NULL && !SCIPintervalIsEntire(inf, activity))
   {
      SCIP_INTERVAL* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPsetAllocBufferArray(set, &varvals, n) );

      for( i = 0; i < n; ++i )
      {
         SCIPintervalSetBounds(&varvals[i], SCIPvarGetLbLocal(SCIPexprtreeGetVars(nlrow->exprtree)[i]), SCIPvarGetUbLocal(SCIPexprtreeGetVars(nlrow->exprtree)[i]));
      }

      SCIP_CALL( SCIPexprtreeEvalInt(nlrow->exprtree, inf, varvals, &bounds) );
      SCIPintervalAdd(inf, &activity, activity, bounds);

      SCIPsetFreeBufferArray(set, &varvals);
   }

   nlrow->minactivity = SCIPintervalGetInf(activity);
   nlrow->maxactivity = SCIPintervalGetSup(activity);

   nlrow->validactivitybdsdomchg = stat->domchgcount;

   return SCIP_OKAY;
}

/** makes sure that there is no fixed variable at position pos of the linear part of a nonlinear row
 * a fixed variable is replaced with the corresponding constant or disaggregated term
 */
static
SCIP_RETCODE nlrowRemoveFixedLinearCoefPos(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   pos                 /**< position of variable in linear variables array */
   )
{
   SCIP_Real oldconstant;
   SCIP_VAR* var;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(pos >= 0);
   assert(pos <  nlrow->nlinvars);

   var = nlrow->linvars[pos];

   if( SCIPvarIsActive(var) )
      return SCIP_OKAY;

   oldconstant = nlrow->constant;

   /* replace fixed, aggregated, or negated variable */
   SCIP_CALL( SCIPvarGetProbvarSum( &nlrow->linvars[pos], set, &nlrow->lincoefs[pos], &nlrow->constant) );

   /* if var had been fixed, entry should be removed from row */
   if( nlrow->lincoefs[pos] == 0.0 )
   {
      nlrowMoveLinearCoef(nlrow, nlrow->nlinvars-1, pos);
      nlrow->nlinvars--;

      if( pos < nlrow->nlinvars )
      {
         SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, pos) );
      }

      return SCIP_OKAY;
   }
   nlrow->linvarssorted = FALSE;

   /* notify nlrow that coefficient of var is now 0.0 in row */
   SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, var, 0.0, nlp) );

   /* notify nlrow that constant of row has changed */
   if( oldconstant != nlrow->constant )  /*lint !e777*/
      SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );

   if( SCIPvarIsActive(nlrow->linvars[pos]) )
   {
      /* if var was aggregated or negated, notify nlrow about new coefficient */
      SCIP_CALL( nlrowLinearCoefChanged(nlrow, set, stat, nlrow->linvars[pos], nlrow->lincoefs[pos], nlp) );
   }
   else
   {
      SCIP_Real coef;
      int i;

      /* if not removed or active, the new variable should be multi-aggregated */
      assert(SCIPvarGetStatus(nlrow->linvars[pos]) == SCIP_VARSTATUS_MULTAGGR);

      var  = nlrow->linvars[pos];
      coef = nlrow->lincoefs[pos];

      /* remove the variable from the row */
      SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, stat, nlp, pos) );

      /* add multi-aggregated term to row */
      if( SCIPvarGetMultaggrConstant(var) != 0.0 )
      {
         nlrow->constant += coef * SCIPvarGetMultaggrConstant(var);
         SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
      }
      SCIP_CALL( SCIPnlrowEnsureLinearSize(nlrow, blkmem, set, nlrow->nlinvars + SCIPvarGetMultaggrNVars(var)) );
      for( i = 0; i < SCIPvarGetMultaggrNVars(var); ++i )
      {
         SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var)[i], coef * SCIPvarGetMultaggrScalars(var)[i]) );
         assert(SCIPvarGetMultaggrVars(var)[i] == nlrow->linvars[nlrow->nlinvars-1]);
         if( !SCIPvarIsActive(SCIPvarGetMultaggrVars(var)[i]) )
         {
            /* if newly added variable is fixed, replace it now */
            SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, nlrow->nlinvars-1) );
         }
      }

      /* due to nlrowDelLinearCoefPos, an inactive variable may have moved to position pos
       * if that is the case, call ourself recursively
       */
      if( pos < nlrow->nlinvars && !SCIPvarIsActive(nlrow->linvars[pos]) )
      {
         SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, pos) );
      }
   }

   return SCIP_OKAY;
}

/** removes fixed variables from the linear part of a nonlinear row */
static
SCIP_RETCODE nlrowRemoveFixedLinearCoefs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   int i;
   int oldlen;

   assert(nlrow != NULL);
   assert(nlrow->linvars != NULL || nlrow->nlinvars == 0);

   oldlen = nlrow->nlinvars;
   for( i = 0; i < MIN(oldlen, nlrow->nlinvars); ++i )
   {
      assert(nlrow->linvars[i] != NULL);
      SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, i) );
   }

   return SCIP_OKAY;
}

/** removes fixed quadratic variables of a nonlinear row by replacing them with the corresponding constant or disaggregated terms */
static
SCIP_RETCODE nlrowRemoveFixedQuadVars(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   int i;
   int nvarsold;
   SCIP_Bool* used;
   SCIP_QUADELEM elem;
   SCIP_QUADELEM newelem;
   int idx2;
   SCIP_Bool havechange;

   SCIP_VAR* var1;
   SCIP_Real coef1;
   SCIP_Real constant1;
   SCIP_VAR* var2;
   SCIP_Real coef2;
   SCIP_Real constant2;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);

   if( nlrow->nquadvars == 0 )
      return SCIP_OKAY;

   SCIPsetDebugMsg(set, "removing fixed quadratic variables from nlrow\n");

   nvarsold = nlrow->nquadvars;
   havechange = FALSE;

   /* allocate array to count number of uses for each variable */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &used, nlrow->nquadvars) );
   BMSclearMemoryArray(used, nlrow->nquadvars);

   i = 0;
   while( i < nlrow->nquadelems )
   {
      elem = nlrow->quadelems[i];

      assert(elem.idx1 < nlrow->nquadvars);
      assert(elem.idx2 < nlrow->nquadvars);
      if( SCIPvarIsActive(nlrow->quadvars[elem.idx1]) && SCIPvarIsActive(nlrow->quadvars[elem.idx2]) )
      {
         /* both variables of quadratic element are active
          * thus, we just remember that we saw them and can continue with the next element
          */
         if( elem.idx1 < nvarsold )
            used[elem.idx1] = TRUE;
         if( elem.idx2 < nvarsold )
            used[elem.idx2] = TRUE;
         ++i;
         continue;
      }

      SCIPsetDebugMsg(set, "removing fixed quadratic variables from %dth element %g <%s> <%s>\n",
         i, elem.coef, SCIPvarGetName(nlrow->quadvars[elem.idx1]), SCIPvarGetName(nlrow->quadvars[elem.idx2]));

      /* if one of the variable is not active, we remove the element and insert new disaggregated ones */
      SCIP_CALL( nlrowDelQuadElemPos(nlrow, set, stat, nlp, i) );
      havechange = TRUE;

      var1 = nlrow->quadvars[elem.idx1];
      var2 = nlrow->quadvars[elem.idx2];
      coef1 = 1.0;
      coef2 = 1.0;
      constant1 = 0.0;
      constant2 = 0.0;

      SCIP_CALL( SCIPvarGetProbvarSum(&var1, set, &coef1, &constant1) );
      SCIP_CALL( SCIPvarGetProbvarSum(&var2, set, &coef2, &constant2) );

      if( coef1 == 0.0 && coef2 == 0.0 )
      {
         /* both variables were fixed, so we may add a constant term and continue */
         if( constant1 != 0.0 && constant2 != 0.0 )
         {
            nlrow->constant += elem.coef * constant1 * constant2;
            SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
         }
         continue;
      }

      if( coef1 == 0.0 )
      {
         /* only the first variable was fixed, so we may add a linear term
          * elem.coef * x * y -> elem.coef * constant1 * (coef2 * var2 + constant2) */
         if( constant1 != 0.0 )
         {
            SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, var2, elem.coef * constant1 * coef2, TRUE) );
            if( constant2 != 0.0 )
            {
               nlrow->constant += elem.coef * constant1 * constant2;
               SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
            }
         }
         /* continue with next element that is at position i now */
         continue;
      }

      if( coef2 == 0.0 )
      {
         /* only the second variable was fixed, so we may add a linear term
          * elem.coef * x * y -> elem.coef * (coef1 * var1 + constant1) * constant2 */
         if( constant2 != 0.0 )
         {
            SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, var1, elem.coef * coef1 * constant2, TRUE) );
            if( constant1 != 0.0 )
            {
               nlrow->constant += elem.coef * constant1 * constant2;
               SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
            }
         }
         /* continue with next element that is at position i now */
         continue;
      }

      if( var1 == var2 && !SCIPvarIsActive(var1) )
      {
         SCIP_Real tmp;
         int* multaggrvaridxs;
         int j, k;

         assert(SCIPvarGetStatus(var1) == SCIP_VARSTATUS_MULTAGGR);
         assert(coef1 == coef2);  /*lint !e777*/
         assert(constant1 == constant2);  /*lint !e777*/
         /* square term which variable is multi-aggregated
          * elem.coef * x^2 -> elem.coef * (coef1 * (multaggrconstant + sum_i multaggrscalar_i*multaggrvar_i) + constant1)^2
          *    = elem.coef * ( (coef1 * multaggrconstant + constant1)^2 +
          *                    2 * (coef1 * multaggrconstant + constant1) * coef1 * (sum_j multaggrscalar_j*multaggrvar_j) +
          *                    coef1^2 * (sum_{j,k} multaggrscalar_j*multaggrscalar_k*multaggrvar_j*multaggrvar_k)
          *                  )
          */

         /* add constant part */
         tmp = coef1 * SCIPvarGetMultaggrConstant(var1) + constant1;
         if( tmp != 0.0 )
         {
            nlrow->constant += elem.coef * tmp * tmp;
            SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
         }

         /* add linear part */
         if( constant1 != 0.0 || SCIPvarGetMultaggrConstant(var1) != 0.0 )
         {
            for( j = 0; j < SCIPvarGetMultaggrNVars(var1); ++j )
            {
               SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var1)[j],
                     2.0 * elem.coef * (coef1 * SCIPvarGetMultaggrConstant(var1) + constant1) * coef1 * SCIPvarGetMultaggrScalars(var1)[j], TRUE) );
            }
         }

         /* setup array with indices of multi-aggregated variables in quadvars */
         SCIP_CALL( SCIPsetAllocBufferArray(set, &multaggrvaridxs, SCIPvarGetMultaggrNVars(var1)) );
         for( j = 0; j < SCIPvarGetMultaggrNVars(var1); ++j )
         {
            multaggrvaridxs[j] = SCIPnlrowSearchQuadVar(nlrow, SCIPvarGetMultaggrVars(var1)[j]);
            if( multaggrvaridxs[j] == -1 )
            {
               /* variable multaggrvar_j not existing in quadvars array yet, so add it */
               SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, blkmem, set, SCIPvarGetMultaggrVars(var1)[j]) );
               multaggrvaridxs[j] = nlrow->nquadvars-1;
            }
            assert(nlrow->quadvars[multaggrvaridxs[j]] == SCIPvarGetMultaggrVars(var1)[j]);
         }

         /* add quadratic elements elem.coef * coef1^2 * (sum_{j,k} multaggrscalar_j*multaggrscalar_k*multaggrvar_j*multaggrvar_k) */
         for( j = 0; j < SCIPvarGetMultaggrNVars(var1); ++j )
         {
            /* bilinear terms */
            for( k = 0; k < j; ++k )
            {
               newelem.idx1 = MIN(multaggrvaridxs[j], multaggrvaridxs[k]);
               newelem.idx2 = MAX(multaggrvaridxs[j], multaggrvaridxs[k]);
               newelem.coef = 2 * elem.coef * coef1 * coef1 * SCIPvarGetMultaggrScalars(var1)[j] * SCIPvarGetMultaggrScalars(var1)[k];
               SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, blkmem, set, stat, nlp, newelem) );
            }

            /* square term */
            newelem.idx1 = multaggrvaridxs[j];
            newelem.idx2 = multaggrvaridxs[j];
            newelem.coef = elem.coef * coef1 * coef1 * SCIPvarGetMultaggrScalars(var1)[j] * SCIPvarGetMultaggrScalars(var1)[j];
            SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, blkmem, set, stat, nlp, newelem) );
         }

         SCIPsetFreeBufferArray(set, &multaggrvaridxs);

         /* continue with next element that is at position i now */
         continue;
      }

      assert(var1 != NULL);
      assert(var2 != NULL);
      if( SCIPvarIsActive(var1) && !SCIPvarIsActive(var2) )
      {
         /* if the second variable is multi-aggregated, but the first one is not, swap both terms */
         SCIP_VAR* tmpvar;
         SCIP_Real tmpcoef;
         SCIP_Real tmpconstant;

         tmpvar      = var1;
         tmpcoef     = coef1;
         tmpconstant = constant1;
         var2      = var1;
         coef2     = coef1;
         constant2 = constant1;
         var1      = tmpvar;
         coef1     = tmpcoef;
         constant1 = tmpconstant;
      }

      if( !SCIPvarIsActive(var1) )
      {
         SCIP_Real tmp;
         int j;

         assert(SCIPvarGetStatus(var1) == SCIP_VARSTATUS_MULTAGGR);

         /* the first variable is multi-aggregated, add a constant and sequences of linear and quadratic terms:
          * elem.coef * x * y -> elem.coef * (coef1 * (multaggrconstant + sum_i multaggrscalar_i*multaggrvar_i) + constant1) * (coef2 * var2 + constant2)
          *    = elem.coef * ( (coef1 * multaggrconstant + constant1) * constant2 +
          *                    (coef1 * multaggrconstant + constant1) * coef2 * var2 +
          *                    (coef1 * (sum_j multaggrscalar_j*multaggrvar_j)) * constant2 +
          *                    (coef1 * (sum_j multaggrscalar_j*multaggrvar_j)) * coef2 * var2
          *                  )
          */

         /* add constant part */
         tmp = elem.coef * (coef1 * SCIPvarGetMultaggrConstant(var1) + constant1) * constant2;
         if( tmp != 0.0 )
         {
            nlrow->constant += tmp;
            SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
         }

         /* add linear part */
         SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, var2, elem.coef * (coef1 * SCIPvarGetMultaggrConstant(var1) + constant1) * coef2, TRUE) );
         if( constant2 != 0.0 )
         {
            for( j = 0; j < SCIPvarGetMultaggrNVars(var1); ++j )
            {
               SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var1)[j], elem.coef * coef1 * SCIPvarGetMultaggrScalars(var1)[j] * constant2, TRUE) );
            }
         }

         /* get index of var2 in quadvars array */
         idx2 = SCIPnlrowSearchQuadVar(nlrow, var2);
         if( idx2 == -1 )
         {
            /* variable var2 not existing in quadvars array yet, so add it */
            SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, blkmem, set, var2) );
            idx2 = nlrow->nquadvars-1;
            assert(nlrow->quadvars[idx2] == var2);
         }

         /* add quadratic elements elem.coef * coef1 * (sum_j multaggrscalar_j*multaggrvar_j) * coef2 * var2 */
         for( j = 0; j < SCIPvarGetMultaggrNVars(var1); ++j )
         {
            newelem.idx1 = SCIPnlrowSearchQuadVar(nlrow, SCIPvarGetMultaggrVars(var1)[j]);
            if( newelem.idx1 == -1 )
            {
               /* variable not existing in quadvars array yet, so add it */
               SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, blkmem, set, SCIPvarGetMultaggrVars(var1)[j]) );
               newelem.idx1 = nlrow->nquadvars-1;
               assert(nlrow->quadvars[newelem.idx1] == SCIPvarGetMultaggrVars(var1)[j]);
            }

            newelem.idx2 = idx2;

            /* swap indices if newelem.idx1 <= newelem.idx2 */
            if( newelem.idx1 > idx2 )
            {
               newelem.idx2 = newelem.idx1;
               newelem.idx1 = idx2;
            }

            newelem.coef = elem.coef * coef1 * coef2 * SCIPvarGetMultaggrScalars(var1)[j];

            SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, blkmem, set, stat, nlp, newelem) );

            /* continue with next element that is at position i now */
            continue;
         }
      }

      assert(SCIPvarIsActive(var1));
      assert(SCIPvarIsActive(var2));
      /* add elem.coef * (coef1 * var1 + constant1) * (coef2 * var2 + constant2) */
      /* add constant part */
      if( constant1 != 0.0 && constant2 != 0.0 )
      {
         nlrow->constant += elem.coef * constant1 * constant2;
         SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
      }
      /* add linear coefficients */
      SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, var1, elem.coef * coef1 * constant2, TRUE) );
      SCIP_CALL( nlrowAddToLinearCoef(nlrow, blkmem, set, stat, nlp, var2, elem.coef * coef2 * constant1, TRUE) );
      /* get index of var1 in quadvars array */
      newelem.idx1 = SCIPnlrowSearchQuadVar(nlrow, var1);
      if( newelem.idx1 == -1 )
      {
         /* variable var2 not existing in quadvars array yet, so add it */
         SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, blkmem, set, var1) );
         newelem.idx1 = nlrow->nquadvars-1;
         assert(nlrow->quadvars[newelem.idx1] == var1);
      }
      /* get index of var2 in quadvars array */
      newelem.idx2 = SCIPnlrowSearchQuadVar(nlrow, var2);
      if( newelem.idx2 == -1 )
      {
         /* variable var2 not existing in quadvars array yet, so add it */
         SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, blkmem, set, var2) );
         newelem.idx2 = nlrow->nquadvars-1;
         assert(nlrow->quadvars[newelem.idx2] == var2);
      }
      /* make sure idx1 <= idx2 */
      if( newelem.idx1 > newelem.idx2 )
      {
         idx2 = newelem.idx2;
         newelem.idx2 = newelem.idx1;
         newelem.idx1 = idx2;
      }
      newelem.coef = elem.coef * coef1 * coef2;
      /* add new quadratic element */
      SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, blkmem, set, stat, nlp, newelem) );

      /* continue with next element that is at position i now */
   }

   /* clean up unused variables */
   if( nlrow->nquadelems == 0 )
   {
      /* the complete quadratic part was fixed or linearized, so we just free up all memory */
      BMSfreeBlockMemoryArray(blkmem, &nlrow->quadvars, nlrow->quadvarssize);
      if( nlrow->quadvarshash != NULL )
         SCIPhashmapFree(&nlrow->quadvarshash);
      BMSfreeBlockMemoryArray(blkmem, &nlrow->quadelems, nlrow->quadelemssize);
      nlrow->nquadvars = 0;
      nlrow->quadvarssize = 0;
      nlrow->nquadelems = 0;
      nlrow->quadelemssize = 0;
      nlrow->quadelemssorted = TRUE;
   }
   else if( havechange )
   {
      /* something had changed, so we likely have quadratic variables to remove */
      int* newpos;
      int offset;

      /* compute new positions of variables in quadvars array */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &newpos, nlrow->nquadvars) );

      offset = 0;
      for( i = 0; i < nvarsold; ++i )
      {
         /* previously existing variables should either be active or not used anymore */
         assert(!used[i] || SCIPvarIsActive(nlrow->quadvars[i]));

         if( !used[i] )
         {
            /* variable has been removed */
            newpos[i] = -1;
            ++offset;
         }
         else
         {
            /* variable will move to position i-offset */
            newpos[i] = i-offset;
         }
      }
      for( ; i < nlrow->nquadvars; ++i )
      {
         if( !SCIPvarIsActive(nlrow->quadvars[i]) )
         {
            /* it can have happened that a new quadratic variable was added that is not active (when multiplying two multi-aggregations)
             * in this case, the variable was only temporarily used and should not be used anymore (this is asserted in the next for-loop below),
             * thus we can remove it
             */
            newpos[i] = -1;
            ++offset;
         }
         else
         {
            /* variable will move to position i-offset */
            newpos[i] = i-offset;
         }
      }

      /* adjust variable indices in quadratic elements */
      for( i = 0; i < nlrow->nquadelems; ++i )
      {
         assert(newpos[nlrow->quadelems[i].idx1] >= 0);
         assert(newpos[nlrow->quadelems[i].idx2] >= 0);
         nlrow->quadelems[i].idx1 = newpos[nlrow->quadelems[i].idx1];
         nlrow->quadelems[i].idx2 = newpos[nlrow->quadelems[i].idx2];
         assert(nlrow->quadelems[i].idx1 <= nlrow->quadelems[i].idx2); /* the way we shrink the quadvars array, variables should stay in the same relative position to each other */
      }

      /* move variables in quadvars array and update quadvarshash */
      for( i = 0; i < nlrow->nquadvars; ++i )
      {
         if( newpos[i] == -1 )
         {
            if( nlrow->quadvarshash != NULL )
            {
               SCIP_CALL( SCIPhashmapRemove(nlrow->quadvarshash, (void*)nlrow->quadvars[i]) );
            }
         }
         else
         {
            nlrow->quadvars[newpos[i]] = nlrow->quadvars[i];
            if( nlrow->quadvarshash != NULL )
            {
               SCIP_CALL( SCIPhashmapSetImage(nlrow->quadvarshash, (void*)nlrow->quadvars[i], (void*)(size_t)newpos[i]) );
            }
         }
      }
      nlrow->nquadvars -= offset;

      SCIPsetFreeBufferArray(set, &newpos);
   }

   SCIPsetFreeBufferArray(set, &used);

   SCIPsetDebugMsg(set, "finished removing fixed quadratic variables\n");

   return SCIP_OKAY;
}

/** removes fixed variables from expression tree of a nonlinear row */
static
SCIP_RETCODE nlrowRemoveFixedExprtreeVars(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   SCIP_Bool changed;

   if( nlrow->exprtree == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPexprtreeRemoveFixedVars(nlrow->exprtree, set, &changed, NULL, NULL) );
   if( changed )
   {
      SCIP_CALL( nlrowExprtreeChanged(nlrow, set, stat, nlp) );
   }

   if( SCIPexprtreeGetNVars(nlrow->exprtree) == 0 && SCIPexprtreeGetNParams(nlrow->exprtree) == 0 )
   {
      /* if expression tree is constant and not parameterized now, remove it */
      SCIP_Real exprval;
      SCIP_CALL( SCIPexprtreeEval(nlrow->exprtree, NULL, &exprval) );
      SCIP_CALL( SCIPnlrowChgConstant(nlrow, set, stat, nlp, nlrow->constant + exprval) );

      SCIP_CALL( SCIPexprtreeFree(&nlrow->exprtree) );
   }

   return SCIP_OKAY;
}

/** removes fixed variable from nonlinear row */
static
SCIP_RETCODE nlrowRemoveFixedVar(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< variable that had been fixed */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);
   assert(!SCIPvarIsActive(var));

   /* search for variable in linear part and remove if existing */
   pos = nlrowSearchLinearCoef(nlrow, var);
   if( pos >= 0 )
   {
      SCIP_CALL( nlrowRemoveFixedLinearCoefPos(nlrow, blkmem, set, stat, nlp, pos) );
   }

   /* search for variable in quadratic part and remove all fixed quadratic variables if existing */
   pos = SCIPnlrowSearchQuadVar(nlrow, var);
   if( pos >= 0 )
   {
      SCIP_CALL( nlrowRemoveFixedQuadVars(nlrow, blkmem, set, stat, nlp) );
   }

   /* search for variable in non-quadratic part and remove all fixed variables in expression tree if existing */
   if( nlrow->exprtree != NULL && SCIPexprtreeFindVar(nlrow->exprtree, var) >= 0 )
   {
      SCIP_CALL( nlrowRemoveFixedExprtreeVars(nlrow, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/*
 * public NLP nonlinear row methods
 */

/** create a new nonlinear row
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreate(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   int                   nquadvars,          /**< number of variables in quadratic terms */
   SCIP_VAR**            quadvars,           /**< variables in quadratic terms, or NULL if nquadvars == 0 */
   int                   nquadelems,         /**< number of entries in quadratic term matrix */
   SCIP_QUADELEM*        quadelems,          /**< elements of quadratic term matrix, or NULL if nquadelems == 0 */
   SCIP_EXPRTREE*        exprtree,           /**< expression tree, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_EXPRCURV         curvature           /**< curvature of the nonlinear row */
   )
{
#ifndef NDEBUG
   int i;
#endif

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(name   != NULL);
   assert(!SCIPsetIsInfinity(set, ABS(constant)));
   assert(nlinvars   == 0 || linvars   != NULL);
   assert(nlinvars   == 0 || lincoefs  != NULL);
   assert(nquadvars  == 0 || quadvars  != NULL);
   assert(nquadelems == 0 || quadelems != NULL);
   assert(nquadelems == 0 || nquadvars > 0);
   assert(SCIPsetIsRelLE(set, lhs, rhs));

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, nlrow) );

   /* constant part */
   assert(!SCIPsetIsInfinity(set, REALABS(constant)));
   (*nlrow)->constant = constant;

#ifndef NDEBUG
   for( i = 0; i < nlinvars; ++i )
   {
      assert(linvars[i] != NULL);
      assert(!SCIPsetIsInfinity(set, REALABS(lincoefs[i])));
   }
#endif

   /* linear part */
   (*nlrow)->nlinvars = nlinvars;
   (*nlrow)->linvarssize = nlinvars;
   if( nlinvars > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->linvars,  linvars,  nlinvars) );
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->lincoefs, lincoefs, nlinvars) );
      (*nlrow)->linvarssorted = FALSE;
   }
   else
   {
      (*nlrow)->linvars  = NULL;
      (*nlrow)->lincoefs = NULL;
      (*nlrow)->linvarssorted = TRUE;
   }

   /* quadratic variables */
#ifndef NDEBUG
   for( i = 0; i < nquadvars; ++i )
      assert(quadvars[i] != NULL);
#endif

   (*nlrow)->nquadvars    = nquadvars;
   (*nlrow)->quadvarssize = nquadvars;
   (*nlrow)->quadvarshash = NULL;
   if( nquadvars > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->quadvars, quadvars, nquadvars) );
      SCIP_CALL( nlrowSetupQuadVarsHash(*nlrow, blkmem) );
   }
   else
   {
      (*nlrow)->quadvars = NULL;
   }

   /* quadratic elements */
#ifndef NDEBUG
   for( i = 0; i < nquadelems; ++i )
   {
      assert(quadelems[i].idx1 >= 0);
      assert(quadelems[i].idx1 <  nquadvars);
      assert(quadelems[i].idx2 >= 0);
      assert(quadelems[i].idx2 <  nquadvars);
      assert(quadelems[i].idx1 <= quadelems[i].idx2);
      assert(!SCIPsetIsInfinity(set, REALABS(quadelems[i].coef)));
   }
#endif

   (*nlrow)->nquadelems = nquadelems;
   (*nlrow)->quadelemssize = nquadelems;
   if( nquadelems > 0 )
   {
      assert(nquadvars > 0);
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->quadelems, quadelems, nquadelems) );
      (*nlrow)->quadelemssorted = FALSE;
   }
   else
   {
      (*nlrow)->quadelems       = NULL;
      (*nlrow)->quadelemssorted = TRUE;
   }

   /* non-quadratic part */
   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeCopy( blkmem, &(*nlrow)->exprtree, exprtree) );
   }
   else
   {
      (*nlrow)->exprtree = NULL;
   }

   /* left and right hand sides, asserted above that lhs is relatively less equal than rhs */
   (*nlrow)->lhs = MIN(lhs, rhs);
   (*nlrow)->rhs = MAX(lhs, rhs);

   /* miscellaneous */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlrow)->name, name, strlen(name)+1) );
   (*nlrow)->activity = SCIP_INVALID;
   (*nlrow)->validactivitynlp = FALSE;
   (*nlrow)->pseudoactivity = SCIP_INVALID;
   (*nlrow)->validpsactivitydomchg = FALSE;
   (*nlrow)->minactivity = SCIP_INVALID;
   (*nlrow)->maxactivity = SCIP_INVALID;
   (*nlrow)->validactivitybdsdomchg = FALSE;
   (*nlrow)->nlpindex = -1;
   (*nlrow)->nlpiindex = -1;
   (*nlrow)->nuses = 0;
   (*nlrow)->dualsol = 0.0;
   (*nlrow)->curvature = curvature;

   /* capture the nonlinear row */
   SCIPnlrowCapture(*nlrow);

   return SCIP_OKAY;
}

/** create a nonlinear row that is a copy of a given row
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreateCopy(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           sourcenlrow         /**< nonlinear row to copy */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(sourcenlrow != NULL);

   SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, sourcenlrow->name,
         sourcenlrow->constant,
         sourcenlrow->nlinvars, sourcenlrow->linvars, sourcenlrow->lincoefs,
         sourcenlrow->nquadvars, sourcenlrow->quadvars, sourcenlrow->nquadelems, sourcenlrow->quadelems,
         sourcenlrow->exprtree,
         sourcenlrow->lhs, sourcenlrow->rhs, sourcenlrow->curvature) );

   (*nlrow)->linvarssorted          = sourcenlrow->linvarssorted;
   (*nlrow)->quadelemssorted        = sourcenlrow->quadelemssorted;
   (*nlrow)->activity               = sourcenlrow->activity;
   (*nlrow)->validactivitynlp       = sourcenlrow->validactivitynlp;
   (*nlrow)->pseudoactivity         = sourcenlrow->pseudoactivity;
   (*nlrow)->validpsactivitydomchg  = sourcenlrow->validpsactivitydomchg;
   (*nlrow)->minactivity            = sourcenlrow->minactivity;
   (*nlrow)->maxactivity            = sourcenlrow->maxactivity;
   (*nlrow)->validactivitybdsdomchg = sourcenlrow->validactivitybdsdomchg;

   return SCIP_OKAY;
}

/** create a new nonlinear row from a linear row
 * the new row is already captured
 */
SCIP_RETCODE SCIPnlrowCreateFromRow(
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< the linear row to copy */
   )
{
   int rownz;

   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(row    != NULL);

   rownz = SCIProwGetNNonz(row);

   if( rownz > 1 )
   {
      SCIP_VAR** rowvars;
      int i;

      SCIP_CALL( SCIPsetAllocBufferArray(set, &rowvars, rownz) );

      for( i = 0; i < rownz; ++i )
      {
         rowvars[i] = SCIPcolGetVar(SCIProwGetCols(row)[i]);
         assert(rowvars[i] != NULL);
      }

      SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, SCIProwGetName(row),
            SCIProwGetConstant(row),
            rownz, rowvars, SCIProwGetVals(row),
            0, NULL, 0, NULL,
            NULL,
            SCIProwGetLhs(row), SCIProwGetRhs(row),
            SCIP_EXPRCURV_LINEAR) );

      SCIPsetFreeBufferArray(set, &rowvars);
   }
   else if( rownz == 1 )
   {
      SCIP_VAR* rowvar;

      rowvar = SCIPcolGetVar(SCIProwGetCols(row)[0]);

      SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, SCIProwGetName(row),
            SCIProwGetConstant(row),
            1, &rowvar, SCIProwGetVals(row),
            0, NULL, 0, NULL,
            NULL,
            SCIProwGetLhs(row), SCIProwGetRhs(row),
            SCIP_EXPRCURV_LINEAR) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowCreate(nlrow, blkmem, set, SCIProwGetName(row),
            SCIProwGetConstant(row),
            0, NULL, NULL,
            0, NULL, 0, NULL,
            NULL,
            SCIProwGetLhs(row), SCIProwGetRhs(row),
            SCIP_EXPRCURV_LINEAR) );
   }

   return SCIP_OKAY;   
}

/** frees a nonlinear row */
SCIP_RETCODE SCIPnlrowFree(
   SCIP_NLROW**          nlrow,              /**< pointer to NLP row */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);
   assert(nlrow  != NULL);
   assert(*nlrow != NULL);
   assert((*nlrow)->nuses == 0);
   assert((*nlrow)->nlpindex == -1);
   assert((*nlrow)->nlpiindex == -1);

   /* linear part */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->linvars,   (*nlrow)->linvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->lincoefs,  (*nlrow)->linvarssize);

   /* quadratic part */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->quadvars,  (*nlrow)->quadvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlrow)->quadelems, (*nlrow)->quadelemssize);
   if( (*nlrow)->quadvarshash != NULL )
      SCIPhashmapFree(&(*nlrow)->quadvarshash);

   /* non-quadratic part */
   if( (*nlrow)->exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&(*nlrow)->exprtree) );
   }

   /* miscellaneous */
   BMSfreeBlockMemoryArray(blkmem, &(*nlrow)->name, strlen((*nlrow)->name)+1);

   BMSfreeBlockMemory(blkmem, nlrow);

   return SCIP_OKAY;
}

/** output nonlinear row to file stream */
SCIP_RETCODE SCIPnlrowPrint(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   assert(nlrow != NULL);

   /* print row name */
   if( nlrow->name != NULL && nlrow->name[0] != '\0' )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "%s: ", nlrow->name);
   }

   /* print left hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, "%.15g <= ", nlrow->lhs);

   /* print constant */
   SCIPmessageFPrintInfo(messagehdlr, file, "%.15g ", nlrow->constant);

   /* print linear coefficients */
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);
      assert(SCIPvarGetName(nlrow->linvars[i]) != NULL);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g<%s> ", nlrow->lincoefs[i], SCIPvarGetName(nlrow->linvars[i]));
   }

   /* print quadratic elements */
   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      assert(SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx1]) != NULL);
      assert(SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx2]) != NULL);
      if( nlrow->quadelems[i].idx1 == nlrow->quadelems[i].idx2 )
         SCIPmessageFPrintInfo(messagehdlr, file, "%+.15gsqr(<%s>) ", nlrow->quadelems[i].coef, SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx1]));
      else
         SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g<%s><%s> ", nlrow->quadelems[i].coef, SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx1]), SCIPvarGetName(nlrow->quadvars[nlrow->quadelems[i].idx2]));
   }

   /* print non-quadratic part */
   if( nlrow->exprtree != NULL )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, " + ");
      SCIP_CALL( SCIPexprtreePrintWithNames(nlrow->exprtree, messagehdlr, file) );
   }

   /* print right hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, "<= %.15g\n", nlrow->rhs);

   return SCIP_OKAY;
}

/** increases usage counter of NLP nonlinear row */
void SCIPnlrowCapture(
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nuses >= 0);

   SCIPdebugMessage("capture nonlinear row <%s> with nuses=%d\n", nlrow->name, nlrow->nuses);
   nlrow->nuses++;
}

/** decreases usage counter of NLP nonlinear row */
SCIP_RETCODE SCIPnlrowRelease(
   SCIP_NLROW**          nlrow,              /**< nonlinear row to free */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(blkmem != NULL);
   assert(nlrow  != NULL);
   assert(*nlrow != NULL);
   assert((*nlrow)->nuses >= 1);

   SCIPsetDebugMsg(set, "release nonlinear row <%s> with nuses=%d\n", (*nlrow)->name, (*nlrow)->nuses);
   (*nlrow)->nuses--;
   if( (*nlrow)->nuses == 0 )
   {
      SCIP_CALL( SCIPnlrowFree(nlrow, blkmem) );
   }

   *nlrow = NULL;

   return SCIP_OKAY;
} /*lint !e715*/

/** ensures, that linear coefficient array of nonlinear row can store at least num entries */
SCIP_RETCODE SCIPnlrowEnsureLinearSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nlinvars <= nlrow->linvarssize);

   if( num > nlrow->linvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->linvars,    nlrow->linvarssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->lincoefs,   nlrow->linvarssize, newsize) );
      nlrow->linvarssize = newsize;
   }
   assert(num <= nlrow->linvarssize);

   return SCIP_OKAY;
}

/** adds a previously non existing linear coefficient to an NLP nonlinear row */
SCIP_RETCODE SCIPnlrowAddLinearCoef(
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   /* if row is in NLP already, make sure that only active variables are added */
   if( nlrow->nlpindex >= 0 )
   {
      SCIP_Real constant;

      /* get corresponding active or multi-aggregated variable */
      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, set, &val, &constant) );

      /* add constant */
      SCIP_CALL( SCIPnlrowChgConstant(nlrow, set, stat, nlp, nlrow->constant + constant) );

      if( val == 0.0 )
         /* var has been fixed */
         return SCIP_OKAY;

      if( !SCIPvarIsActive(var) )
      {
         /* var should be multi-aggregated, so call this function recursively */
         int i;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);
         for( i = 0; i < SCIPvarGetMultaggrNVars(var); ++i )
         {
            SCIP_CALL( SCIPnlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, SCIPvarGetMultaggrVars(var)[i], SCIPvarGetMultaggrScalars(var)[i] * val) );
         }
         return SCIP_OKAY;
      }

      /* var is active, so can go on like normal */
   }

   SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, var, val) );

   return SCIP_OKAY;
}

/** deletes linear coefficient from nonlinear row */
SCIP_RETCODE SCIPnlrowDelLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var                 /**< coefficient to be deleted */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);

   /* if the row is in the NLP already, we can only have active variables, so var should also be active; in non-debug mode, one gets an error below */
   assert(nlrow->nlpindex == -1 || SCIPvarIsActive(var) );

   /* search the position of the variable in the row's variable vector */
   pos = nlrowSearchLinearCoef(nlrow, var);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for variable <%s> doesn't exist in nonlinear row <%s>\n", SCIPvarGetName(var), nlrow->name);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < nlrow->nlinvars);
   assert(nlrow->linvars[pos] == var);

   /* delete the variable from the row's variable vector */
   SCIP_CALL( nlrowDelLinearCoefPos(nlrow, set, stat, nlp, pos) );

   return SCIP_OKAY;
}

/** changes or adds a linear coefficient to a nonlinear row */
SCIP_RETCODE SCIPnlrowChgLinearCoef(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(nlp != NULL);
   assert(var != NULL);

   /* search the position of the variable in the row's linvars vector */
   pos = nlrowSearchLinearCoef(nlrow, var);

   /* check, if column already exists in the row's linear variables vector */
   if( pos == -1 )
   {
      if( !SCIPsetIsZero(set, coef) )
      {
         /* add previously not existing coefficient */
         SCIP_CALL( nlrowAddLinearCoef(nlrow, blkmem, set, stat, nlp, var, coef) );
      }
   }
   else
   {
      /* change the coefficient in the row */
      SCIP_CALL( nlrowChgLinearCoefPos(nlrow, set, stat, nlp, pos, coef) );
   }

   return SCIP_OKAY;
}

/** ensures, that quadratic variables array of nonlinear row can store at least num entries */
SCIP_RETCODE SCIPnlrowEnsureQuadVarsSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nquadvars <= nlrow->quadvarssize);

   if( num > nlrow->quadvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->quadvars, nlrow->quadvarssize, newsize) );
      nlrow->quadvarssize = newsize;
   }
   assert(num <= nlrow->quadvarssize);

   return SCIP_OKAY;
}

/** adds variable to quadvars array of row */
SCIP_RETCODE SCIPnlrowAddQuadVar(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable to search for */
   )
{
   assert(blkmem != NULL);
   assert(nlrow  != NULL);
   assert(var    != NULL);

   /* assert that only active variables are added once the row is in the NLP */
   assert(nlrow->nlpindex == -1 || SCIPvarIsActive(var) );

   /* assert that variable has not been added already */
   assert(SCIPnlrowSearchQuadVar(nlrow, var) == -1);

   SCIP_CALL( SCIPnlrowEnsureQuadVarsSize(nlrow, blkmem, set, nlrow->nquadvars+1) );
   nlrow->quadvars[nlrow->nquadvars] = var;
   nlrow->nquadvars++;

   if( nlrow->quadvarshash == NULL )
   {
      SCIP_CALL( nlrowSetupQuadVarsHash(nlrow, blkmem) );
   }
   else
   {
      SCIP_CALL( SCIPhashmapInsert(nlrow->quadvarshash, (void*)var, (void*)(size_t)(nlrow->nquadvars-1)) );
   }
   assert(SCIPnlrowSearchQuadVar(nlrow, var) == nlrow->nquadvars-1);

   return SCIP_OKAY;
}

/** ensures, that quadratic elements array of nonlinear row can store at least num entries */
SCIP_RETCODE SCIPnlrowEnsureQuadElementsSize(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlrow != NULL);
   assert(nlrow->nquadelems <= nlrow->quadelemssize);

   if( num > nlrow->quadelemssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlrow->quadelems, nlrow->quadelemssize, newsize) );
      nlrow->quadelemssize = newsize;
   }
   assert(num <= nlrow->quadelemssize);

   return SCIP_OKAY;
}

/** adds a previously non existing quadratic element to an NLP nonlinear row */
SCIP_RETCODE SCIPnlrowAddQuadElement(
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_QUADELEM         elem                /**< quadratic element to add */
   )
{
   SCIP_CALL( nlrowAddQuadElement(nlrow, blkmem, set, stat, nlp, elem) );

   return SCIP_OKAY;
}

/** deletes quadratic element from nonlinear row */
SCIP_RETCODE SCIPnlrowDelQuadElement(
   SCIP_NLROW*           nlrow,              /**< nonlinear row to be changed */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   idx1,               /**< index of first variable in element */
   int                   idx2                /**< index of second variable in element */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(idx1 >= 0);
   assert(idx1 <  nlrow->nquadvars);
   assert(idx2 >= 0);
   assert(idx2 <  nlrow->nquadvars);
   assert(idx1 <= idx2);

   /* search the position of the variable in the row's variable vector */
   pos = nlrowSearchQuadElem(nlrow, idx1, idx2);
   if( pos == -1 )
   {
      SCIPerrorMessage("coefficient for index pair (idx1, idx2) doesn't exist in nonlinear row <%s>\n", idx1, idx2, nlrow->name);
      return SCIP_INVALIDDATA;
   }
   assert(0 <= pos && pos < nlrow->nquadelems);

   /* delete the element from the row's quadratic elements array */
   SCIP_CALL( nlrowDelQuadElemPos(nlrow, set, stat, nlp, pos) );

   return SCIP_OKAY;
}

/** changes or adds a quadratic element to a nonlinear row */
SCIP_RETCODE SCIPnlrowChgQuadElem(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_QUADELEM         elem                /**< new quadratic element */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(nlp != NULL);

   /* search the position of the variable in the row's linvars vector */
   pos = nlrowSearchQuadElem(nlrow, elem.idx1, elem.idx2);

   if( pos == -1 )
   {
      /* add previously not existing element */
      SCIP_CALL( nlrowAddQuadElement(nlrow, blkmem, set, stat, nlp, elem) );
   }
   else
   {
      /* change the coefficient in the row */
      SCIP_CALL( nlrowChgQuadElemPos(nlrow, set, stat, nlp, pos, elem.coef) );
   }

   return SCIP_OKAY;
}

/** replaces an expression tree in nonlinear row */
SCIP_RETCODE SCIPnlrowChgExprtree(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_EXPRTREE*        exprtree            /**< new expression tree */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);

   /* free previous expression tree */
   if( nlrow->exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&nlrow->exprtree) );
      assert(nlrow->exprtree == NULL);
   }

   /* adds new expression tree */
   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeCopy(blkmem, &nlrow->exprtree, exprtree) );

      /* if row is already in NLP, ensure that exprtree has only active variables */
      if( nlrow->nlpindex >= 0 )
      {
         SCIP_Bool dummy;
         SCIP_CALL( SCIPexprtreeRemoveFixedVars(nlrow->exprtree, set, &dummy, NULL, NULL) );
      }
   }

   /* notify row about the change */
   SCIP_CALL( nlrowExprtreeChanged(nlrow, set, stat, nlp) );

   return SCIP_OKAY;
}

/** changes a parameter in an expression of a nonlinear row */
SCIP_RETCODE SCIPnlrowChgExprtreeParam(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   int                   paramidx,           /**< index of parameter in expression tree's parameter array */
   SCIP_Real             paramval            /**< new value of parameter */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(nlrow->exprtree != NULL);

   SCIPexprtreeSetParamVal(nlrow->exprtree, paramidx, paramval);

   /* notify row about the change */
   SCIP_CALL( nlrowExprtreeParamChanged(nlrow, set, stat, paramidx, nlp) );

   return SCIP_OKAY;
}

/** changes all parameters in an expression of a nonlinear row */
SCIP_RETCODE SCIPnlrowChgExprtreeParams(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            paramvals           /**< new values of parameters */
   )
{
   assert(nlrow  != NULL);
   assert(blkmem != NULL);
   assert(nlrow->exprtree != NULL);

   SCIP_CALL( SCIPexprtreeSetParams(nlrow->exprtree, SCIPexprtreeGetNParams(nlrow->exprtree), paramvals) );

   /* notify row about the change */
   SCIP_CALL( nlrowExprtreeParamChanged(nlrow, set, stat, -1, nlp) );

   return SCIP_OKAY;
}

/** changes constant of nonlinear row */
SCIP_RETCODE SCIPnlrowChgConstant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             constant            /**< new constant */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->constant, constant) )
   {
      nlrow->constant = constant;
      SCIP_CALL( nlrowConstantChanged(nlrow, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/** changes left hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgLhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->lhs, lhs) )
   {
      nlrow->lhs = lhs;
      SCIP_CALL( nlrowSideChanged(nlrow, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/** changes right hand side of nonlinear row */
SCIP_RETCODE SCIPnlrowChgRhs(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   assert(nlrow != NULL);

   if( !SCIPsetIsEQ(set, nlrow->rhs, rhs) )
   {
      nlrow->rhs = rhs;
      SCIP_CALL( nlrowSideChanged(nlrow, set, stat, nlp) );
   }

   return SCIP_OKAY;
}

/** removes (or substitutes) all fixed, negated, aggregated, multi-aggregated variables from the linear, quadratic, and non-quadratic terms of a nonlinear row */
SCIP_RETCODE SCIPnlrowRemoveFixedVars(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   SCIP_CALL( nlrowRemoveFixedLinearCoefs(nlrow, blkmem, set, stat, nlp) );
   SCIP_CALL( nlrowRemoveFixedQuadVars(nlrow, blkmem, set, stat, nlp) );
   SCIP_CALL( nlrowRemoveFixedExprtreeVars(nlrow, set, stat, nlp) );

   return SCIP_OKAY;
}

/** recalculates the current activity of a nonlinear row */
SCIP_RETCODE SCIPnlrowRecalcNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   SCIP_Real val1, val2;
   int i;
   int previdx1;

   assert(nlrow  != NULL);
   assert(stat   != NULL);
   assert(nlp    != NULL);

   if( nlp->solstat > SCIP_NLPSOLSTAT_LOCINFEASIBLE )
   {
      SCIPerrorMessage("do not have NLP solution for computing NLP activity\n");
      return SCIP_ERROR;
   }

   nlrow->activity = nlrow->constant;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);
      assert(SCIPvarGetNLPSol(nlrow->linvars[i]) < SCIP_INVALID);

      nlrow->activity += nlrow->lincoefs[i] * SCIPvarGetNLPSol(nlrow->linvars[i]);
   }

   val1 = 0.0; /* for lint */
   previdx1 = -1;
   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      /* if first index of quadelems is the same as in last round, val1 is still up to date */
      if( previdx1 != nlrow->quadelems[i].idx1 )
      {
         previdx1 = nlrow->quadelems[i].idx1;
         val1 = SCIPvarGetNLPSol(nlrow->quadvars[previdx1]);
         assert(val1 < SCIP_INVALID);

         if( val1 == 0.0 )
            continue;
      }

      val2 = SCIPvarGetNLPSol(nlrow->quadvars[nlrow->quadelems[i].idx2]);
      assert(val2 < SCIP_INVALID);

      nlrow->activity += nlrow->quadelems[i].coef * val1 * val2;
   }

   if( nlrow->exprtree != NULL )
   {
      SCIP_Real* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPsetAllocBufferArray(set, &varvals, n) );

      for( i = 0; i < n; ++i )
      {
         varvals[i] = SCIPvarGetNLPSol(SCIPexprtreeGetVars(nlrow->exprtree)[i]);
      }

      SCIP_CALL( SCIPexprtreeEval(nlrow->exprtree, varvals, &val1) );
      nlrow->activity += val1;

      SCIPsetFreeBufferArray(set, &varvals);
   }

   nlrow->validactivitynlp = stat->nnlps;

   return SCIP_OKAY;
}

/** returns the activity of a nonlinear row in the current NLP solution */
SCIP_RETCODE SCIPnlrowGetNLPActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   assert(nlrow  != NULL);
   assert(stat   != NULL);
   assert(activity != NULL);

   assert(nlrow->validactivitynlp <= stat->nnlps);

   if( nlrow->validactivitynlp != stat->nnlps )
   {
      SCIP_CALL( SCIPnlrowRecalcNLPActivity(nlrow, set, stat, nlp) );
   }
   assert(nlrow->validactivitynlp == stat->nnlps);
   assert(nlrow->activity < SCIP_INVALID);

   *activity = nlrow->activity;

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row in the current NLP solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetNLPFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_Real activity;

   assert(nlrow != NULL);
   assert(feasibility != NULL);

   SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, set, stat, nlp, &activity) );
   *feasibility = MIN(nlrow->rhs - activity, activity - nlrow->lhs);

   return SCIP_OKAY;
}

/** calculates the current pseudo activity of a nonlinear row */
SCIP_RETCODE SCIPnlrowRecalcPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Real val1, val2;
   int i;

   assert(nlrow  != NULL);
   assert(stat   != NULL);

   nlrow->pseudoactivity = nlrow->constant;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);

      val1 = SCIPvarGetBestBoundLocal(nlrow->linvars[i]);
      nlrow->pseudoactivity += nlrow->lincoefs[i] * val1;
   }

   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      val1 = SCIPvarGetBestBoundLocal(nlrow->quadvars[nlrow->quadelems[i].idx1]);
      if( val1 == 0.0 )
         continue;

      val2 = SCIPvarGetBestBoundLocal(nlrow->quadvars[nlrow->quadelems[i].idx2]);
      nlrow->pseudoactivity += nlrow->quadelems[i].coef * val1 * val2;
   }

   if( nlrow->exprtree != NULL )
   {
      SCIP_Real* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPsetAllocBufferArray(set, &varvals, n) );

      for( i = 0; i < n; ++i )
         varvals[i] = SCIPvarGetBestBoundLocal(SCIPexprtreeGetVars(nlrow->exprtree)[i]);

      SCIP_CALL( SCIPexprtreeEval(nlrow->exprtree, varvals, &val1) );
      nlrow->pseudoactivity += val1;

      SCIPsetFreeBufferArray(set, &varvals);
   }

   nlrow->validpsactivitydomchg = stat->domchgcount;

   return SCIP_OKAY;
}

/** returns the pseudo activity of a nonlinear row in the current pseudo solution */
SCIP_RETCODE SCIPnlrowGetPseudoActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudoactivity      /**< buffer to store pseudo activity value */
   )
{
   assert(nlrow != NULL);
   assert(stat  != NULL);
   assert(pseudoactivity != NULL);
   assert(nlrow->validpsactivitydomchg <= stat->domchgcount);

   /* check, if pseudo activity has to be calculated */
   if( nlrow->validpsactivitydomchg != stat->domchgcount )
   {
      SCIP_CALL( SCIPnlrowRecalcPseudoActivity(nlrow, set, stat) );
   }
   assert(nlrow->validpsactivitydomchg == stat->domchgcount);
   assert(nlrow->pseudoactivity < SCIP_INVALID);

   *pseudoactivity = nlrow->pseudoactivity;

   return SCIP_OKAY;
}

/** returns the pseudo feasibility of a nonlinear row in the current pseudo solution: negative value means infeasibility */
SCIP_RETCODE SCIPnlrowGetPseudoFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudofeasibility   /**< buffer to store pseudo feasibility value */
   )
{
   SCIP_Real pseudoactivity;

   assert(nlrow != NULL);
   assert(stat  != NULL);
   assert(pseudofeasibility != NULL);

   SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, set, stat, &pseudoactivity) );
   *pseudofeasibility = MIN(nlrow->rhs - pseudoactivity, pseudoactivity - nlrow->lhs);

   return SCIP_OKAY;
}

/** returns the activity of a nonlinear row for a given solution */
SCIP_RETCODE SCIPnlrowGetSolActivity(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   SCIP_Real inf;
   SCIP_Real val1, val2;
   int i;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(activity != NULL);

   *activity = nlrow->constant;
   for( i = 0; i < nlrow->nlinvars; ++i )
   {
      assert(nlrow->linvars[i] != NULL);

      val1 = SCIPsolGetVal(sol, set, stat, nlrow->linvars[i]);
      if( val1 == SCIP_UNKNOWN ) /*lint !e777*/
      {
         *activity = SCIP_INVALID;
         return SCIP_OKAY;
      }
      *activity += nlrow->lincoefs[i] * val1;
   }

   for( i = 0; i < nlrow->nquadelems; ++i )
   {
      val1 = SCIPsolGetVal(sol, set, stat, nlrow->quadvars[nlrow->quadelems[i].idx1]);
      if( val1 == SCIP_UNKNOWN ) /*lint !e777*/
      {
         *activity = SCIP_INVALID;
         return SCIP_OKAY;
      }
      if( val1 == 0.0 )
         continue;

      val2 = SCIPsolGetVal(sol, set, stat, nlrow->quadvars[nlrow->quadelems[i].idx2]);
      if( val2 == SCIP_UNKNOWN ) /*lint !e777*/
      {
         *activity = SCIP_INVALID;
         return SCIP_OKAY;
      }
      *activity += nlrow->quadelems[i].coef * val1 * val2;
   }

   if( nlrow->exprtree != NULL )
   {
      SCIP_Real* varvals;
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);

      SCIP_CALL( SCIPsetAllocBufferArray(set, &varvals, n) );

      for( i = 0; i < n; ++i )
      {
         varvals[i] = SCIPsolGetVal(sol, set, stat, SCIPexprtreeGetVars(nlrow->exprtree)[i]);
         if( varvals[i] == SCIP_UNKNOWN ) /*lint !e777*/
         {
            *activity = SCIP_INVALID;
            SCIPsetFreeBufferArray(set, &varvals);
            return SCIP_OKAY;
         }
      }

      SCIP_CALL( SCIPexprtreeEval(nlrow->exprtree, varvals, &val1) );
      *activity += val1;

      SCIPsetFreeBufferArray(set, &varvals);
   }

   inf = SCIPsetInfinity(set);
   *activity = MAX(*activity, -inf);
   *activity = MIN(*activity, +inf);

   return SCIP_OKAY;
}

/** returns the feasibility of a nonlinear row for the given solution */
SCIP_RETCODE SCIPnlrowGetSolFeasibility(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_Real activity;

   assert(nlrow != NULL);
   assert(feasibility != NULL);

   SCIP_CALL( SCIPnlrowGetSolActivity(nlrow, set, stat, sol, &activity) );

   *feasibility = MIN(nlrow->rhs - activity, activity - nlrow->lhs);

   return SCIP_OKAY;
}

/** returns the minimal activity of a nonlinear row w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowGetActivityBounds(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   )
{
   assert(nlrow != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(nlrow->validactivitybdsdomchg <= stat->domchgcount);

   /* check, if activity bounds has to be calculated */
   if( nlrow->validactivitybdsdomchg != stat->domchgcount )
   {
      SCIP_CALL( nlrowCalcActivityBounds(nlrow, set, stat) );
   }
   assert(nlrow->validactivitybdsdomchg == stat->domchgcount);
   assert(nlrow->minactivity < SCIP_INVALID);
   assert(nlrow->maxactivity < SCIP_INVALID);

   if( minactivity != NULL )
      *minactivity = nlrow->minactivity;
   if( maxactivity != NULL )
      *maxactivity = nlrow->maxactivity;

   return SCIP_OKAY;
}

/** returns whether the nonlinear row is redundant w.r.t. the variables' bounds */
SCIP_RETCODE SCIPnlrowIsRedundant(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_Bool*            isredundant         /**< buffer to store whether row is redundant */
   )
{
   SCIP_Real minactivity;
   SCIP_Real maxactivity;

   assert(nlrow != NULL);
   assert(set != NULL);
   assert(isredundant != NULL);

   SCIP_CALL( SCIPnlrowGetActivityBounds(nlrow, set, stat, &minactivity, &maxactivity) );

   *isredundant = TRUE;
   if( (!SCIPsetIsInfinity(set, -nlrow->lhs) && SCIPsetIsFeasLT(set, minactivity, nlrow->lhs)) ||
      ( !SCIPsetIsInfinity(set,  nlrow->rhs) && SCIPsetIsFeasGT(set, maxactivity, nlrow->rhs)) )
      *isredundant = FALSE;

   return SCIP_OKAY;
}

/** gets constant */
SCIP_Real SCIPnlrowGetConstant(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->constant;
}

/** gets number of variables of linear part */
int SCIPnlrowGetNLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlinvars;
}

/** gets array with variables of linear part */
SCIP_VAR** SCIPnlrowGetLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->linvars;
}

/** gets array with coefficients in linear part */
SCIP_Real* SCIPnlrowGetLinearCoefs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->lincoefs;
}

/** gets number of quadratic variables in quadratic part */
int SCIPnlrowGetNQuadVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nquadvars;
}

/** gets quadratic variables in quadratic part */
SCIP_VAR** SCIPnlrowGetQuadVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->quadvars;
}

/** gives position of variable in quadvars array of row, or -1 if not found */
int SCIPnlrowSearchQuadVar(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_VAR*             var                 /**< variable to search for */
   )
{
   int pos;

   assert(nlrow != NULL);
   assert(var   != NULL);

   if( nlrow->quadvarshash != NULL )
   {
      pos = SCIPhashmapExists(nlrow->quadvarshash, var) ? (int)(size_t)SCIPhashmapGetImage(nlrow->quadvarshash, var) : -1;
   }
   else
   {
      for( pos = nlrow->nquadvars-1; pos >= 0; --pos )
         if( nlrow->quadvars[pos] == var )
            break;
   }

   assert(pos == -1 || (pos < nlrow->nquadvars && nlrow->quadvars[pos] == var));

   return pos;
}

/** gets number of quadratic elements in quadratic part */
int SCIPnlrowGetNQuadElems(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nquadelems;
}

/** gets quadratic elements in quadratic part */
SCIP_QUADELEM* SCIPnlrowGetQuadElems(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->quadelems;
}

/** gets array with coefficients in linear part */
void SCIPnlrowGetQuadData(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int*                  nquadvars,          /**< buffer to store number of variables in quadratic term, or NULL if not of interest */
   SCIP_VAR***           quadvars,           /**< buffer to store pointer to array of variables in quadratic term, or NULL if not of interest */
   int*                  nquadelems,         /**< buffer to store number of entries in quadratic term, or NULL if not of interest */
   SCIP_QUADELEM**       quadelems           /**< buffer to store pointer to array of entries in quadratic term, or NULL if not of interest */
   )
{
   assert(nlrow != NULL);

   if( nquadvars  != NULL )
      *nquadvars  = nlrow->nquadvars;
   if( quadvars   != NULL )
      *quadvars   = nlrow->quadvars;
   if( nquadelems != NULL )
      *nquadelems = nlrow->nquadelems;
   if( quadelems  != NULL )
      *quadelems  = nlrow->quadelems;
}

/** gets expression tree */
SCIP_EXPRTREE* SCIPnlrowGetExprtree(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->exprtree;
}

/** returns the left hand side of a nonlinear row */
SCIP_Real SCIPnlrowGetLhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->lhs;
}

/** returns the right hand side of a nonlinear row */
SCIP_Real SCIPnlrowGetRhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->rhs;
}

/** returns the curvature of a nonlinear row */
SCIP_EXPRCURV SCIPnlrowGetCurvature(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);
   return nlrow->curvature;
}

/** sets the curvature of a nonlinear row */
void SCIPnlrowSetCurvature(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_EXPRCURV         curvature           /**< curvature of NLP row */
   )
{
   assert(nlrow != NULL);
   nlrow->curvature = curvature;
}

/** returns the name of a nonlinear row */
const char* SCIPnlrowGetName(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->name;
}

/** gets position of a nonlinear row in current NLP, or -1 if not in NLP */
int SCIPnlrowGetNLPPos(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpindex;
}

/** returns TRUE iff row is member of current NLP */
SCIP_Bool SCIPnlrowIsInNLP(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpindex != -1;
}

/** gets the dual NLP solution of a nlrow
 * for a ranged constraint, the dual value is positive if the right hand side is active and negative if the left hand side is active
 */
SCIP_Real SCIPnlrowGetDualsol(
   SCIP_NLROW*           nlrow               /**< NLP row */
   )
{
   assert(nlrow != NULL);

   return nlrow->nlpiindex >= 0 ? nlrow->dualsol : 0.0;
}

/*
 * private NLP methods
 */

/** announces, that a row of the NLP was modified
 * adjusts status of current solution
 * calling method has to ensure that change is passed to the NLPI!
 */
static
SCIP_RETCODE nlpRowChanged(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row which was changed */
   )
{
   assert(nlp != NULL);
   assert(nlrow != NULL);
   assert(!nlp->indiving);
   assert(nlrow->nlpindex >= 0);

   /* nlrow is a row in the NLP, so changes effect feasibility */
   /* if we have a feasible NLP solution and it satisfies the modified row, then it is still feasible
    * if the NLP was globally or locally infeasible or unbounded, then this may not be the case anymore
    */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      SCIP_Real feasibility;
      SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, set, stat, nlp, &feasibility) );
      if( !SCIPsetIsFeasNegative(set, feasibility) )
         nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
      else
         nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
   }
   else
   {
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
   }

   return SCIP_OKAY;
}

/** adds a set of nonlinear rows to the NLP and captures them */
static
SCIP_RETCODE nlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of nonlinear rows to add */
   SCIP_NLROW**          nlrows              /**< nonlinear rows to add */
   )
{
#ifndef NDEBUG
   int i;
#endif
   int j;
   SCIP_NLROW* nlrow;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlrows != NULL || nnlrows == 0);
   assert(!nlp->indiving);

   SCIP_CALL( SCIPnlpEnsureNlRowsSize(nlp, blkmem, set, nlp->nnlrows + nnlrows) );

   for( j = 0; j < nnlrows; ++j )
   {
      nlrow = nlrows[j];  /*lint !e613*/

      /* assert that row is not in NLP (or even NLPI) yet */
      assert(nlrow->nlpindex == -1);
      assert(nlrow->nlpiindex == -1);

      /* make sure there are only active variables in row */
      SCIP_CALL( SCIPnlrowRemoveFixedVars(nlrow, blkmem, set, stat, nlp) );

#ifndef NDEBUG
      /* assert that variables of row are in NLP */
      for( i = 0; i < nlrow->nlinvars; ++i )
         assert(SCIPhashmapExists(nlp->varhash, nlrow->linvars[i]));

      for( i = 0; i < nlrow->nquadvars; ++i )
         assert(SCIPhashmapExists(nlp->varhash, nlrow->quadvars[i]));

      if( nlrow->exprtree )
      {
         int n;

         n = SCIPexprtreeGetNVars(nlrow->exprtree);
         assert(SCIPexprtreeGetVars(nlrow->exprtree) != NULL || n == 0);

         for( i = 0; i < n; ++i )
            assert(SCIPhashmapExists(nlp->varhash, SCIPexprtreeGetVars(nlrow->exprtree)[i]));
      }
#endif

      /* add row to NLP and capture it */
      nlp->nlrows[nlp->nnlrows + j] = nlrow;
      nlrow->nlpindex = nlp->nnlrows + j;

      SCIPnlrowCapture(nlrow);

      /* if we have a feasible NLP solution and it satisfies the new solution, then it is still feasible
       * if the NLP was globally or locally infeasible, then it stays that way
       * if the NLP was unbounded, then this may not be the case anymore
       */
      if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         SCIP_Real feasibility;
         SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, set, stat, nlp, &feasibility) );
         if( !SCIPsetIsFeasNegative(set, feasibility) )
            nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
         else
            nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
      }
      else if( nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      {
         nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
      }
   }

   nlp->nnlrows += nnlrows;
   nlp->nunflushednlrowadd += nnlrows;

   return SCIP_OKAY;
}

/** moves a nonlinear row to a different place, and updates all corresponding data structures */
static
void nlpMoveNlrow(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   int                   oldpos,             /**< old position of nonlinear row */
   int                   newpos              /**< new position of nonlinear row */
   )
{
   assert(nlp != NULL);
   assert(0 <= oldpos && oldpos < nlp->nnlrows);
   assert(0 <= newpos && newpos < nlp->nnlrows);
   assert(nlp->nlrows[oldpos] != NULL);

   if( oldpos == newpos )
      return;

   nlp->nlrows[newpos] = nlp->nlrows[oldpos];
   nlp->nlrows[newpos]->nlpindex = newpos;
}

/** deletes nonlinear row with given position from NLP */
static
SCIP_RETCODE nlpDelNlRowPos(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   pos                 /**< position of nonlinear row that is to be removed */
   )
{
   SCIP_NLROW* nlrow;

   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(pos >= 0);
   assert(pos < nlp->nnlrows);
   assert(!nlp->indiving);

   nlrow = nlp->nlrows[pos];
   assert(nlrow != NULL);
   assert(nlrow->nlpindex == pos);

   /* if row is in NLPI, then mark that it has to be removed in the next flush
    * if row was not in NLPI yet, then we have one unflushed nlrow addition less */
   if( nlrow->nlpiindex >= 0 )
   {
      assert(nlrow->nlpiindex < nlp->nnlrows_solver);
      nlp->nlrowmap_nlpi2nlp[nlrow->nlpiindex] = -1;
      nlrow->nlpiindex = -1;
      ++nlp->nunflushednlrowdel;
   }
   else
   {
      assert(nlrow->nlpiindex == -1);
      --nlp->nunflushednlrowadd;
   }

   /* move NLP row from the end to pos and mark nlrow to be not in NLP anymore */
   nlpMoveNlrow(nlp, nlp->nnlrows-1, pos);
   nlrow->nlpindex = -1;

   /* forget about restriction */
   SCIP_CALL( SCIPnlrowRelease(&nlrow, blkmem, set) );
   --nlp->nnlrows;

   if( nlp->solstat < SCIP_NLPSOLSTAT_LOCOPT )
      nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
   else if( nlp->solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;

   return SCIP_OKAY;
}

/** updates bounds on a variable in the NLPI problem */
static
SCIP_RETCODE nlpUpdateVarBounds(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable which bounds have changed */
   SCIP_Bool             tightened           /**< whether the bound change was a bound tightening */
   )
{
   int pos;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* original variable bounds are ignored during diving
    * (all variable bounds are reset to their current value in exitDiving) */
   if( nlp->indiving )
      return SCIP_OKAY;

   /* get position of variable in NLP */
   pos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);

   /* if variable not in NLPI yet, nothing to do */
   if( nlp->varmap_nlp2nlpi[pos] == -1 )
      return SCIP_OKAY;

   /* update bounds in NLPI problem */
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   pos = nlp->varmap_nlp2nlpi[pos];
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, 1, &pos, &lb, &ub) );

   /* if we have a feasible NLP solution and it satisfies the new bounds, then it is still feasible
    * if the NLP was globally or locally infeasible and we tightened a bound, then it stays that way
    * if the NLP was unbounded and we tightened a bound, then this may not be the case anymore
    */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
   {
      if( !tightened ||
         ((SCIPsetIsInfinity(set, -lb) || SCIPsetIsFeasLE(set, lb, SCIPvarGetNLPSol(var))) &&
          (SCIPsetIsInfinity(set,  ub) || SCIPsetIsFeasGE(set, ub, SCIPvarGetNLPSol(var)))) )
         nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
      else
         nlp->solstat = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
   }
   else if( !tightened || nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
   {
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;
   }

   return SCIP_OKAY;
}

/** updates coefficient of a variable in the objective */
static
SCIP_RETCODE nlpUpdateObjCoef(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_VAR*             var                 /**< variable which bounds have changed */
   )
{
   int pos;
   int objidx;
   SCIP_Real coef;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* if the objective in the NLPI is not up to date, then we do not need to do something here */
   if( !nlp->objflushed )
      return SCIP_OKAY;

   /* original objective is ignored during diving
    * we just need to remember that at end of diving we have to flush the objective */
   if( nlp->indiving )
   {
      nlp->objflushed = FALSE;
      return SCIP_OKAY;
   }

   /* get position of variable in NLP and objective coefficient */
   pos  = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);
   assert(nlp->varmap_nlp2nlpi[pos] == -1 || nlp->solver != NULL);

   /* actually we only need to remember flushing the objective if we also have an NLPI */
   if( nlp->solver == NULL )
      return SCIP_OKAY;

   coef = SCIPvarGetObj(var);

   /* if variable not in NLPI yet, then we only need to remember to update the objective after variable additions were flushed */
   if( nlp->varmap_nlp2nlpi[pos] == -1 && coef != 0.0 )
   {
      nlp->objflushed = FALSE;

      return SCIP_OKAY;
   }

   /* if we are here, then the objective in the NLPI is up to date,
    * we keep it this way by changing the coefficient of var in the NLPI problem objective */
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   pos = nlp->varmap_nlp2nlpi[pos];
   objidx = -1;
   SCIP_CALL( SCIPnlpiChgLinearCoefs(nlp->solver, nlp->problem, objidx, 1, &pos, &coef) );

   /* if we had a solution and it was locally (or globally) optimal, then now we can only be sure that it is still feasible */
   if( nlp->solstat < SCIP_NLPSOLSTAT_FEASIBLE )
      nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;

   return SCIP_OKAY;
}

/** adds new variables to the NLP */
static
SCIP_RETCODE nlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variable to add to NLP */
   )
{
   int i;
   SCIP_VAR* var;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(vars   != NULL || nvars == 0);
   assert(!nlp->indiving || nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPnlpEnsureVarsSize(nlp, blkmem, set, nlp->nvars + nvars) );
   assert(nlp->sizevars >= nlp->nvars + nvars);

   for( i = 0; i < nvars; ++i )
   {
      var = vars[i];  /*lint !e613*/

      assert(SCIPvarIsTransformed(var));
      assert(SCIPvarIsActive(var));
      assert(!SCIPhashmapExists(nlp->varhash, var));

      SCIPvarCapture(var);

      nlp->vars[nlp->nvars+i]            = var;
      nlp->varmap_nlp2nlpi[nlp->nvars+i] = -1;
      SCIP_CALL( SCIPhashmapInsert(nlp->varhash, var, (void*) (size_t) (nlp->nvars+i)) );

      nlp->varlbdualvals[nlp->nvars+i]   = 0.0;
      nlp->varubdualvals[nlp->nvars+i]   = 0.0;

      /* update objective, if necessary (new variables have coefficient 0.0 anyway) */
      if( SCIPvarGetObj(var) != 0.0 )
      {
         SCIP_CALL( nlpUpdateObjCoef(nlp, var) );
      }

      /* let's keep the previous initial guess and set it for the new variable to the best bound
       * (since there can be no row that uses this variable yet, this seems a good guess) */
      if( nlp->haveinitguess )
      {
         assert(nlp->initialguess != NULL);

         nlp->initialguess[nlp->nvars+i] = SCIPvarGetBestBoundLocal(var);
      }

      /* if we have a feasible NLP solution, then it remains feasible
       * but we have to update the objective function
       */
      if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         SCIP_CALL( SCIPvarSetNLPSol(var, set, SCIPvarGetBestBoundLocal(var)) );
         nlp->primalsolobjval += SCIPvarGetObj(var) * SCIPvarGetBestBoundLocal(var);
         nlp->solstat = SCIP_NLPSOLSTAT_FEASIBLE;
      }

      /* catch events on variable */
      SCIP_CALL( SCIPvarCatchEvent(var, blkmem, set, \
            SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED, \
            nlp->eventhdlr, (SCIP_EVENTDATA*)nlp, NULL) ); /* @todo should store event filter position in nlp? */
   }

   nlp->nvars += nvars;
   nlp->nunflushedvaradd += nvars;

   return SCIP_OKAY;
}

/** moves a variable to a different place, and updates all corresponding data structures */
static
SCIP_RETCODE nlpMoveVar(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   int                   oldpos,             /**< old position of variable */
   int                   newpos              /**< new position of variable */
   )
{
   int nlpipos;

   assert(nlp != NULL);
   assert(0 <= oldpos && oldpos < nlp->nvars);
   assert(0 <= newpos && newpos < nlp->nvars);
   assert(nlp->vars[oldpos] != NULL);

   if( oldpos == newpos )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapSetImage(nlp->varhash, nlp->vars[oldpos], (void*) (size_t) newpos) );
   nlp->vars[newpos]            = nlp->vars[oldpos];
   nlp->varmap_nlp2nlpi[newpos] = nlp->varmap_nlp2nlpi[oldpos];
   nlp->varlbdualvals[newpos]   = nlp->varlbdualvals[oldpos];
   nlp->varubdualvals[newpos]   = nlp->varubdualvals[oldpos];
   if( nlp->initialguess != NULL )
      nlp->initialguess[newpos] = nlp->initialguess[oldpos];

   nlpipos = nlp->varmap_nlp2nlpi[newpos];
   if( nlpipos > 0 )
      nlp->varmap_nlpi2nlp[nlpipos] = newpos;

   return SCIP_OKAY;
}

/** deletes variable with given position from NLP */
static
SCIP_RETCODE nlpDelVarPos(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed if a column-variable is freed */
   int                   pos                 /**< position of nonlinear row that is to be removed */
   )
{
   SCIP_VAR* var;
#ifndef NDEBUG
   int i;
#endif
   int nlpipos;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(pos >= 0);
   assert(pos < nlp->nvars);
   assert(!nlp->indiving);

   var = nlp->vars[pos];
   assert(var != NULL);

#ifndef NDEBUG
   /* assert that variable is not used by any nonlinear row */
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      int j;
      SCIP_NLROW* nlrow;

      nlrow = nlp->nlrows[i];
      assert(nlrow != NULL);

      /* use nlrowSearchLinearCoef only if already sorted, since otherwise we may change the solving process slightly */
      if( nlrow->linvarssorted )
         assert( nlrowSearchLinearCoef(nlrow, var) == -1 );
      else
         for( j = 0; j < nlrow->nlinvars; ++j )
            assert( nlrow->linvars[j] != var );

      assert( SCIPnlrowSearchQuadVar(nlrow, var) == -1);

      assert(nlrow->exprtree == NULL || SCIPexprtreeFindVar(nlrow->exprtree, var) == -1);
   }
#endif

   /* if we had a feasible solution, then adjust objective function value
    * if NLP was unbounded before, then maybe it is not anymore */
   if( nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      nlp->primalsolobjval -= SCIPvarGetObj(var) * SCIPvarGetNLPSol(var);
   else if( nlp->solstat == SCIP_NLPSOLSTAT_UNBOUNDED )
      nlp->solstat = SCIP_NLPSOLSTAT_UNKNOWN;

   /* if variable is in NLPI problem, mark that we have to remember to delete it there
    * if it was not in the NLPI yet, then we have one unflushed var addition less now */
   nlpipos = nlp->varmap_nlp2nlpi[pos];
   if( nlpipos >= 0 )
   {
      assert(nlpipos < nlp->nvars_solver);

      nlp->varmap_nlpi2nlp[nlpipos] = -1;
      ++nlp->nunflushedvardel;
   }
   else
      --nlp->nunflushedvaradd;

   /* drop events on variable */
   SCIP_CALL( SCIPvarDropEvent(var, blkmem, set, \
         SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_OBJCHANGED, \
         nlp->eventhdlr, (SCIP_EVENTDATA*)nlp, -1) );

   /* move variable from end to pos */
   SCIP_CALL( nlpMoveVar(nlp, nlp->nvars-1, pos) );

   /* forget about variable */
   SCIP_CALL( SCIPhashmapRemove(nlp->varhash, var) );
   SCIP_CALL( SCIPvarRelease(&var, blkmem, set, eventqueue, lp) );
   --nlp->nvars;

   return SCIP_OKAY;
}

/** notifies NLP that a variable was fixed, so it is removed from objective, all rows, and the NLP variables */
static
SCIP_RETCODE nlpRemoveFixedVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed to release variable */
   SCIP_VAR*             var                 /**< variable that has been fixed */
   )
{
   int i;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(!SCIPvarIsActive(var));
   assert(!nlp->indiving);
   assert(SCIPhashmapExists(nlp->varhash, var));

   /* remove var from all rows */
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIP_CALL( nlrowRemoveFixedVar(nlp->nlrows[i], blkmem, set, stat, nlp, var) );
   }

   /* remove variable from NLP */
   SCIP_CALL( SCIPnlpDelVar(nlp, blkmem, set, eventqueue, lp, var) );

   return SCIP_OKAY;
}

/** creates arrays with NLPI variable indices of variables in a nonlinear row */
static
SCIP_RETCODE nlpSetupNlpiIndices(
   SCIP_NLP*             nlp,                /**< NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   int**                 linidxs,            /**< buffer to store pointer to NLPI indices of linear variables */
   SCIP_QUADELEM**       quadelems,          /**< buffer to store pointer to quadratic elements w.r.t. NLPI indices */
   int**                 nlinidxs            /**< buffer to store pointer to NLPI indices of nonlinear variables */
   )
{
   int i;
   SCIP_VAR* var;

   assert(nlp    != NULL);
   assert(set    != NULL);
   assert(nlrow  != NULL);
   assert(linidxs   != NULL);
   assert(quadelems != NULL);
   assert(nlinidxs  != NULL);

   /* get indices of variables in linear part of row */
   if( nlrow->nlinvars > 0 )
   {
      assert(nlrow->linvars  != NULL);
      assert(nlrow->lincoefs != NULL);

      SCIP_CALL( SCIPsetAllocBufferArray(set, linidxs, nlrow->nlinvars) );

      for( i = 0; i < nlrow->nlinvars; ++i )
      {
         var = nlrow->linvars[i];
         assert(var != NULL);
         assert(SCIPvarIsActive(var)); /* at this point, there should be only active variables in the row */

         assert(SCIPhashmapExists(nlp->varhash, var));
         (*linidxs)[i] = nlp->varmap_nlp2nlpi[(size_t) (void*) SCIPhashmapGetImage(nlp->varhash, var)];
         assert((*linidxs)[i] >= 0);
      }
   }
   else
      *linidxs = NULL;

   /* get indices of variables in quadratic part of row */
   if( nlrow->nquadvars > 0 )
   {
      int* quadvarsidx;

      assert(nlrow->quadvars    != NULL);
      assert(nlrow->nquadelems  > 0);
      assert(nlrow->quadelems   != NULL);

      /* compute mapping of variable indices quadratic term -> NLPI */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &quadvarsidx, nlrow->nquadvars) );
      for( i = 0; i < nlrow->nquadvars; ++i )
      {
         var = nlrow->quadvars[i];
         assert(var != NULL);
         assert(SCIPvarIsActive(var)); /* at this point, there should be only active variables in the row */

         assert(SCIPhashmapExists(nlp->varhash, var));
         quadvarsidx[i] = nlp->varmap_nlp2nlpi[(size_t) (void*) SCIPhashmapGetImage(nlp->varhash, var)];
      }

      /* compute quad elements using NLPI indices */
      SCIP_CALL( SCIPsetAllocBufferArray(set, quadelems, nlrow->nquadelems) );
      for( i = 0; i < nlrow->nquadelems; ++i )
      {
         assert(nlrow->quadelems[i].idx1 >= 0);
         assert(nlrow->quadelems[i].idx1 < nlrow->nquadvars);
         assert(nlrow->quadelems[i].idx2 >= 0);
         assert(nlrow->quadelems[i].idx2 < nlrow->nquadvars);

         (*quadelems)[i].idx1 = quadvarsidx[nlrow->quadelems[i].idx1];
         (*quadelems)[i].idx2 = quadvarsidx[nlrow->quadelems[i].idx2];
         if( (*quadelems)[i].idx1 > (*quadelems)[i].idx2 )
         {
            int tmp = (*quadelems)[i].idx1;
            (*quadelems)[i].idx1 = (*quadelems)[i].idx2;
            (*quadelems)[i].idx2 = tmp;
         }
         (*quadelems)[i].coef = nlrow->quadelems[i].coef;
      }

      SCIPsetFreeBufferArray(set, &quadvarsidx);
   }
   else
      *quadelems = NULL;

   /* get indices of variables in expression tree part of row */
   if( nlrow->exprtree != NULL )
   {
      int n;

      n = SCIPexprtreeGetNVars(nlrow->exprtree);
      assert(n == 0 || SCIPexprtreeGetVars(nlrow->exprtree) != NULL);

      SCIP_CALL( SCIPsetAllocBufferArray(set, nlinidxs, n) );

      for( i = 0; i < n; ++i )
      {
         var = SCIPexprtreeGetVars(nlrow->exprtree)[i];
         assert(var != NULL);
         assert(SCIPvarIsActive(var)); /* at this point, there should be only active variables in the row */

         assert(SCIPhashmapExists(nlp->varhash, var));
         (*nlinidxs)[i] = nlp->varmap_nlp2nlpi[(size_t) (void*) SCIPhashmapGetImage(nlp->varhash, var)];
      }
   }
   else
      *nlinidxs = NULL;

   return SCIP_OKAY;
}

/** ensures, that NLPI variables array of NLP can store at least num entries */
static
SCIP_RETCODE nlpEnsureVarsSolverSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nvars_solver <= nlp->sizevars_solver);

   if( num > nlp->sizevars_solver )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varmap_nlpi2nlp, nlp->sizevars_solver, newsize) );

      nlp->sizevars_solver = newsize;
   }
   assert(num <= nlp->sizevars_solver);

   return SCIP_OKAY;
}

/** ensures, that NLPI nonlinear rows array of NLP can store at least num entries */
static
SCIP_RETCODE nlpEnsureNlRowsSolverSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nnlrows_solver <= nlp->sizenlrows_solver);

   if( num > nlp->sizenlrows_solver )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->nlrowmap_nlpi2nlp, nlp->sizenlrows_solver, newsize) );

      nlp->sizenlrows_solver = newsize;
   }
   assert(num <= nlp->sizenlrows_solver);

   return SCIP_OKAY;
}

/** deletes rows from the NLPI problem that have been marked as to remove */
static
SCIP_RETCODE nlpFlushNlRowDeletions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int         j;
   int         c;      /* counts the number of rows to delete */
   int*        rowset; /* marks which rows to delete and stores new indices */
   SCIP_NLROW* nlrow;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushednlrowdel >= 0);
   assert(!nlp->indiving);

   if( nlp->nunflushednlrowdel == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending removals of nonlinear rows */
      for( j = 0; j < nlp->nnlrows_solver; ++j )
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* create marker which rows have to be deleted */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rowset, nlp->nnlrows_solver) );
   c = 0;
   for( j = 0; j < nlp->nnlrows_solver; ++j )
   {
      if( nlp->nlrowmap_nlpi2nlp[j] == -1 )
      {
         rowset[j] = 1;
         ++c;
      }
      else
         rowset[j] = 0;
   }
   assert(c == nlp->nunflushednlrowdel);

   /* remove rows from NLPI problem */
   SCIP_CALL( SCIPnlpiDelConsSet(nlp->solver, nlp->problem, rowset, nlp->nnlrows_solver) );

   /* update NLPI row indices */
   for( j = 0; j < nlp->nnlrows_solver; ++j )
   {
      assert(rowset[j] <= j); /* we assume that the NLP solver did not move a row behind its previous position!! */
      if( rowset[j] < 0 )
      {
         /* assert that row was marked as deleted */
         assert(nlp->nlrowmap_nlpi2nlp[j] == -1);
      }
      else if( rowset[j] < j )
      {
         /* nlrow at position j moved (forward) to position rowset[j] */
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
         assert(nlp->nlrowmap_nlpi2nlp[j] < nlp->nnlrows);

         nlrow = nlp->nlrows[nlp->nlrowmap_nlpi2nlp[j]];
         assert(nlrow->nlpiindex == j);

         /* there should be no row at the new position already */
         assert(nlp->nlrowmap_nlpi2nlp[rowset[j]] == -1);

         nlrow->nlpiindex = rowset[j];
         nlp->nlrowmap_nlpi2nlp[rowset[j]] = nlrow->nlpindex;
      }
      else
      {
         /* row j stays at position j */
         assert(nlp->nlrowmap_nlpi2nlp[j] >= 0);
         assert(nlp->nlrowmap_nlpi2nlp[j] < nlp->nnlrows);
         assert(nlp->nlrows[nlp->nlrowmap_nlpi2nlp[j]]->nlpiindex == j);
      }
   }
   nlp->nnlrows_solver -= c;
   nlp->nunflushednlrowdel = 0;

   /* cleanup */
   SCIPsetFreeBufferArray(set, &rowset);

   return SCIP_OKAY;
}

/** deletes variables from the NLPI problem that have been marked as to remove
 * assumes that there are no pending row deletions (nlpFlushNlRowDeletions should be called first)
 */
static
SCIP_RETCODE nlpFlushVarDeletions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int  i;
   int  c;      /* counter on number of variables to remove in solver */
   int* colset; /* marks which variables to delete and stores new indices */

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushedvardel >= 0);
   assert(nlp->nunflushednlrowdel == 0);
   assert(!nlp->indiving);

   if( nlp->nunflushedvardel == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending removals of variables */
      for( i = 0; i < nlp->nvars_solver; ++i )
         assert(nlp->varmap_nlpi2nlp[i] >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* create marker which variables have to be deleted */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &colset, nlp->nvars_solver) );
   c = 0;
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      if( nlp->varmap_nlpi2nlp[i] == -1 )
      {
         colset[i] = 1;
         ++c;
      }
      else
         colset[i] = 0;
   }
   assert(c == nlp->nunflushedvardel);

   /* delete variables from NLPI problem */
   SCIP_CALL( SCIPnlpiDelVarSet(nlp->solver, nlp->problem, colset, nlp->nvars_solver) );

   /* update NLPI variable indices */
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      assert(colset[i] <= i); /* we assume that the NLP solver did not move a variable behind its previous position!! */
      if( colset[i] < 0 )
      {
         /* assert that variable was marked as deleted */
         assert(nlp->varmap_nlpi2nlp[i] == -1);
      }
      else if( colset[i] < i)
      {
         /* variable at position i moved (forward) to position colset[i] */
         int varpos;

         varpos = nlp->varmap_nlpi2nlp[i]; /* position of variable i in NLP */
         assert(varpos >= 0);
         assert(varpos < nlp->nvars);
         assert(nlp->varmap_nlp2nlpi[varpos] == i);

         /* there should be no variable at the new position already */
         assert(nlp->varmap_nlpi2nlp[colset[i]] == -1);

         nlp->varmap_nlp2nlpi[varpos] = colset[i];
         nlp->varmap_nlpi2nlp[colset[i]] = varpos;
      }
      else
      {
         /* variable i stays at position i */
         assert(nlp->varmap_nlpi2nlp[i] >= 0);
         assert(nlp->varmap_nlpi2nlp[i] < nlp->nvars);
         assert(nlp->varmap_nlp2nlpi[nlp->varmap_nlpi2nlp[i]] == i);
      }
   }

   nlp->nvars_solver -= c;
   nlp->nunflushedvardel = 0;

   /* cleanup */
   SCIPsetFreeBufferArray(set, &colset);

   return SCIP_OKAY;
}

/** adds nonlinear rows to NLPI problem that have been added to NLP before
 * assumes that there are no pending variable additions or deletions (nlpFlushVarDeletions and nlpFlushVarAdditions should be called first) */
static
SCIP_RETCODE nlpFlushNlRowAdditions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int c, i;
   SCIP_NLROW* nlrow;
   SCIP_Real*  lhss;
   SCIP_Real*  rhss;
   int*        nlinvars;
   int**       linidxs;
   SCIP_Real** lincoefs;
   int*        nquadelems;
   SCIP_QUADELEM** quadelems;
   int**       nlidxs;
   SCIP_EXPRTREE** exprtrees;
   const char** names;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushednlrowadd >= 0);
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->nunflushedvardel == 0);
   assert(!nlp->indiving);

   if( nlp->nunflushednlrowadd == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending additions of variables */
      for( i = 0; i < nlp->nnlrows; ++i )
         assert(nlp->nlrows[i]->nlpiindex >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( nlpEnsureNlRowsSolverSize(nlp, blkmem, set, nlp->nnlrows_solver + nlp->nunflushednlrowadd) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, &lhss,        nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rhss,        nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &nlinvars,    nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &linidxs,     nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lincoefs,    nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &nquadelems,  nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &quadelems,   nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &nlidxs,      nlp->nunflushednlrowadd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &exprtrees,   nlp->nunflushednlrowadd) );
#if ADDNAMESTONLPI
   SCIP_CALL( SCIPsetAllocBufferArray(set, &names,       nlp->nunflushednlrowadd) );
#else
   names = NULL;
#endif

   c = 0;
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      nlrow = nlp->nlrows[i];
      assert(nlrow != NULL);

      /* skip nonlinear rows already in NLPI problem */
      if( nlrow->nlpiindex >= 0 )
         continue;
      assert(c < nlp->nunflushednlrowadd);

      /* get indices in NLPI */
      SCIP_CALL( nlpSetupNlpiIndices(nlp, set, nlrow, &linidxs[c], &quadelems[c], &nlidxs[c]) );
      assert(linidxs[c]   != NULL || nlrow->nlinvars  == 0);
      assert(quadelems[c] != NULL || nlrow->nquadvars == 0);
      assert(nlidxs[c]    != NULL || nlrow->exprtree  == NULL);

      nlp->nlrowmap_nlpi2nlp[nlp->nnlrows_solver+c] = i;
      nlrow->nlpiindex = nlp->nnlrows_solver+c;

      lhss[c] = nlrow->lhs;
      rhss[c] = nlrow->rhs;
      if( nlrow->constant != 0.0 )
      {
         if( !SCIPsetIsInfinity(set, -nlrow->lhs) )
            lhss[c] -= nlrow->constant;
         if( !SCIPsetIsInfinity(set,  nlrow->rhs) )
            rhss[c] -= nlrow->constant;
      }
      if( rhss[c] < lhss[c] )
      {
         assert(SCIPsetIsEQ(set, lhss[c], rhss[c]));
         rhss[c] = lhss[c];
      }

      nlinvars[c] = nlrow->nlinvars;
      lincoefs[c] = nlrow->lincoefs;

      nquadelems[c] = nlrow->nquadelems;

      exprtrees[c]  = nlrow->exprtree;

#if ADDNAMESTONLPI
      names[c]      = nlrow->name;
#endif

      ++c;

#ifdef NDEBUG
      /* have c vars to add already, there can be no more */
      if( c == nlp->nunflushednlrowadd )
         break;
#endif
   }
   assert(c == nlp->nunflushednlrowadd);

   nlp->nnlrows_solver += c;

   SCIP_CALL( SCIPnlpiAddConstraints(nlp->solver, nlp->problem, c, lhss, rhss,
         nlinvars, linidxs, lincoefs,
         nquadelems, quadelems,
         nlidxs, exprtrees,
         names) );

   for( c = 0; c < nlp->nunflushednlrowadd; ++c )
   {
      if( linidxs[c] != NULL )
         SCIPsetFreeBufferArray(set, &linidxs[c]);
      if( quadelems[c] != NULL )
         SCIPsetFreeBufferArray(set, &quadelems[c]);
      if( nlidxs[c] != NULL )
         SCIPsetFreeBufferArray(set, &nlidxs[c]);
   }

#if ADDNAMESTONLPI
   SCIPsetFreeBufferArray(set, &names);
#endif
   SCIPsetFreeBufferArray(set, &lhss);
   SCIPsetFreeBufferArray(set, &rhss);
   SCIPsetFreeBufferArray(set, &nlinvars);
   SCIPsetFreeBufferArray(set, &linidxs);
   SCIPsetFreeBufferArray(set, &lincoefs);
   SCIPsetFreeBufferArray(set, &nquadelems);
   SCIPsetFreeBufferArray(set, &quadelems);
   SCIPsetFreeBufferArray(set, &nlidxs);
   SCIPsetFreeBufferArray(set, &exprtrees);

   nlp->nunflushednlrowadd = 0;

   return SCIP_OKAY;
}


/** adds variables to NLPI problem that have been added to NLP before
 * may set nlp->objflushed to FALSE if a variable with nonzero obj.coefficient is added to the NLPI problem */
static
SCIP_RETCODE nlpFlushVarAdditions(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i, c;
   SCIP_Real*  lbs;
   SCIP_Real*  ubs;
   const char** names;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushedvaradd >= 0);
   assert(!nlp->indiving);

   if( nlp->nunflushedvaradd == 0 )
   {
#ifndef NDEBUG
      /* check that there are really no pending additions of variables */
      for( i = 0; i < nlp->nvars; ++i )
         assert(nlp->varmap_nlp2nlpi[i] >= 0);
#endif
      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( nlpEnsureVarsSolverSize(nlp, blkmem, set, nlp->nvars_solver + nlp->nunflushedvaradd) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbs,   nlp->nunflushedvaradd) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubs,   nlp->nunflushedvaradd) );
#if ADDNAMESTONLPI
   SCIP_CALL( SCIPsetAllocBufferArray(set, &names, nlp->nunflushedvaradd) );
#else
   names = NULL;
#endif

   c = 0;
   for( i = 0; i < nlp->nvars; ++i )
   {
      /* skip variables already in NLPI problem */
      if( nlp->varmap_nlp2nlpi[i] >= 0 )
         continue;
      assert(c < nlp->nunflushedvaradd);

      nlp->varmap_nlpi2nlp[nlp->nvars_solver+c] = i;
      nlp->varmap_nlp2nlpi[i] = nlp->nvars_solver+c;
      lbs[c]   = SCIPvarGetLbLocal(nlp->vars[i]);
      ubs[c]   = SCIPvarGetUbLocal(nlp->vars[i]);
#if ADDNAMESTONLPI
      names[c] = SCIPvarGetName(nlp->vars[i]);
#endif
      ++c;

      /* if the new variable has a nonzero objective coefficient, then the objective need to be updated */
      if( !SCIPsetIsZero(set, SCIPvarGetObj(nlp->vars[i])) )
         nlp->objflushed = FALSE;

#ifdef NDEBUG
      /* have c vars to add already, there can be no more */
      if( c == nlp->nunflushedvaradd )
         break;
#endif
   }
   assert(c == nlp->nunflushedvaradd);

   nlp->nvars_solver += c;

   SCIP_CALL( SCIPnlpiAddVars(nlp->solver, nlp->problem, c, lbs, ubs, names) );

#if ADDNAMESTONLPI
   SCIPsetFreeBufferArray(set, &names);
#endif
   SCIPsetFreeBufferArray(set, &lbs);
   SCIPsetFreeBufferArray(set, &ubs);

   nlp->nunflushedvaradd = 0;

   return SCIP_OKAY;
}

/** updates the objective in the NLPI problem, if necessary
 * assumes that there are no unflushed variable additions or deletions (nlpFlushVarDeletions and nlpFlushVarAdditions should be called first)
 */
static
SCIP_RETCODE nlpFlushObjective(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int*       linindices;
   SCIP_Real* lincoefs;
   SCIP_Real  coef;
   int        i;
   int        nz;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->nunflushedvardel == 0);
   assert(!nlp->indiving);

   if( nlp->objflushed )
      return SCIP_OKAY;

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* assemble coefficients */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &linindices, nlp->nvars_solver) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lincoefs,   nlp->nvars_solver) );

   nz = 0;
   for( i = 0; i < nlp->nvars_solver; ++i )
   {
      assert(nlp->varmap_nlpi2nlp[i] >= 0); /* there should be no variable deletions pending */

      coef = SCIPvarGetObj(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
      if( SCIPsetIsZero(set, coef) )
         continue;

      linindices[nz] = i;
      lincoefs[nz]   = coef;
      ++nz;
   }

   SCIP_CALL( SCIPnlpiSetObjective(nlp->solver, nlp->problem,
         nz, linindices, lincoefs,
         0, NULL,
         NULL, NULL,
         0.0) );

   SCIPsetFreeBufferArray(set, &linindices);
   SCIPsetFreeBufferArray(set, &lincoefs);

   nlp->objflushed = TRUE;

   return SCIP_OKAY;
}

/** solves the NLP, assuming it has been flushed already
 *
 *  is used also to solve diving NLP
 */
static
SCIP_RETCODE nlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( nlp->solver == NULL )
   {
      SCIPmessagePrintWarning(messagehdlr, "Attempted to solve NLP, but no solver available.\n");

      nlp->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      nlp->termstat = SCIP_NLPTERMSTAT_OTHER;

      return SCIP_OKAY;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* set initial guess, if available */
   if( nlp->haveinitguess )
   {
      /* @todo should we not set it if we had set it already? (initguessflushed...) */
      SCIP_Real* initialguess_solver;
      int nlpidx;

      assert(nlp->initialguess != NULL);

      SCIP_CALL( SCIPsetAllocBufferArray(set, &initialguess_solver, nlp->nvars_solver) );

      for( i = 0; i < nlp->nvars_solver; ++i )
      {
         nlpidx = nlp->varmap_nlpi2nlp[i];
         assert(nlpidx >= 0);
         assert(nlpidx < nlp->nvars);

         initialguess_solver[i] = nlp->initialguess[nlpidx];
      }
      SCIP_CALL( SCIPnlpiSetInitialGuess(nlp->solver, nlp->problem, initialguess_solver, NULL, NULL, NULL) );

      SCIPsetFreeBufferArray(set, &initialguess_solver);
   }

   /* set NLP tolerances to current SCIP primal and dual feasibility tolerance */
   SCIP_CALL( SCIPnlpiSetRealPar(nlp->solver, nlp->problem, SCIP_NLPPAR_FEASTOL, SCIPsetFeastol(set)) );
   SCIP_CALL( SCIPnlpiSetRealPar(nlp->solver, nlp->problem, SCIP_NLPPAR_RELOBJTOL, SCIPsetDualfeastol(set)) );

   /* let NLP solver do his work */
   SCIPclockStart(stat->nlpsoltime, set);

   SCIP_CALL( SCIPnlpiSolve(nlp->solver, nlp->problem) );

   SCIPclockStop(stat->nlpsoltime, set);
   ++stat->nnlps;

   nlp->termstat = SCIPnlpiGetTermstat(nlp->solver, nlp->problem);
   nlp->solstat  = SCIPnlpiGetSolstat(nlp->solver, nlp->problem);
   switch( nlp->solstat )
   {
   case SCIP_NLPSOLSTAT_GLOBOPT:
   case SCIP_NLPSOLSTAT_LOCOPT:
   case SCIP_NLPSOLSTAT_FEASIBLE:
   case SCIP_NLPSOLSTAT_LOCINFEASIBLE:
   {
      SCIP_Real* primalvals;
      SCIP_Real* nlrowdualvals;
      SCIP_Real* varlbdualvals;
      SCIP_Real* varubdualvals;

      primalvals    = NULL;
      nlrowdualvals = NULL;
      varlbdualvals = NULL;
      varubdualvals = NULL;

      /* get NLP solution */
      SCIP_CALL( SCIPnlpiGetSolution(nlp->solver, nlp->problem, &primalvals, &nlrowdualvals, &varlbdualvals, &varubdualvals, NULL) );
      assert(primalvals != NULL || nlp->nvars == 0);
      assert((varlbdualvals != NULL) == (varubdualvals != NULL)); /* if there are duals for one bound, then there should also be duals for the other bound */

      /* store solution primal values in variable and evaluate objective function */
      if( nlp->indiving && nlp->divingobj != NULL )
      {
         for( i = 0; i < nlp->nvars; ++i )
         {
            SCIP_CALL( SCIPvarSetNLPSol(nlp->vars[i], set, primalvals[nlp->varmap_nlp2nlpi[i]]) );  /*lint !e613 */
         }

         /* evaluate modified diving objective */
         SCIP_CALL( SCIPnlrowGetNLPActivity(nlp->divingobj, set, stat, nlp, &nlp->primalsolobjval) );
      }
      else
      {
         /* evaluate SCIP objective function */
         nlp->primalsolobjval = 0.0;
         for( i = 0; i < nlp->nvars; ++i )
         {
            SCIP_Real solval = primalvals[nlp->varmap_nlp2nlpi[i]];  /*lint !e613 */

            /* do a quick assert that variable bounds are satisfied, if feasibility is claimed */
            assert(SCIPsetIsFeasGE(set, solval, SCIPvarGetLbLocal(nlp->vars[i])) || nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE);
            assert(SCIPsetIsFeasLE(set, solval, SCIPvarGetUbLocal(nlp->vars[i])) || nlp->solstat > SCIP_NLPSOLSTAT_FEASIBLE);

            SCIP_CALL( SCIPvarSetNLPSol(nlp->vars[i], set, solval) );  /*lint !e613 */
            nlp->primalsolobjval += SCIPvarGetObj(nlp->vars[i]) * solval;  /*lint !e613 */
         }
      }

      /* store solution dual values in nlrows and variables */
      for( i = 0; i < nlp->nnlrows; ++i )
      {
         assert(nlp->nlrows[i]->nlpiindex >= 0); /* NLP was flushed before solve, so all nlrows should be in there */

         nlp->nlrows[i]->dualsol = nlrowdualvals != NULL ? nlrowdualvals[nlp->nlrows[i]->nlpiindex] : 0.0;

         /* SCIPsetDebugMsg(set, "dual of nlrow <%s> = %g\n", nlp->nlrows[i]->name, nlp->nlrows[i]->dualsol); */
      }
      assert(nlp->varlbdualvals != NULL || nlp->nvars == 0);
      assert(nlp->varubdualvals != NULL || nlp->nvars == 0);
      if( varlbdualvals != NULL )
      {
         for( i = 0; i < nlp->nvars; ++i )
         {
            assert(nlp->varmap_nlp2nlpi[i] >= 0); /* NLP was flushed before solve, so all vars should be in there */

            nlp->varlbdualvals[i] = varlbdualvals[nlp->varmap_nlp2nlpi[i]];
            nlp->varubdualvals[i] = varubdualvals[nlp->varmap_nlp2nlpi[i]];

            /* SCIPsetDebugMsg(set, "duals of var <%s> = %g %g\n", SCIPvarGetName(nlp->vars[i]), nlp->varlbdualvals[i], nlp->varubdualvals[i]); */
         }
      }
      else if( nlp->nvars > 0 )
      {
         BMSclearMemoryArray(nlp->varlbdualvals, nlp->nvars);
         BMSclearMemoryArray(nlp->varubdualvals, nlp->nvars);
      }

      break;
   }
   default:
      nlp->primalsolobjval = SCIP_INVALID;
      break;
   } /*lint !e788*/

   return SCIP_OKAY;
}

/** assembles list of fractional variables in last NLP solution */
static
SCIP_RETCODE nlpCalcFracVars(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(nlp->validfracvars <= stat->nnlps);
   assert(SCIPnlpHasSolution(nlp));

   SCIPsetDebugMsg(set, "calculating NLP fractional variables: validfracvars=%" SCIP_LONGINT_FORMAT ", nnlps=%" SCIP_LONGINT_FORMAT "\n", nlp->validfracvars, stat->nnlps);

   if( nlp->solstat > SCIP_NLPSOLSTAT_LOCINFEASIBLE )
   {
      nlp->nfracvars     = 0;
      nlp->npriofracvars = 0;
      nlp->validfracvars = stat->nnlps;

      SCIPsetDebugMsg(set, "NLP globally infeasible, unbounded, or worse -> no solution values -> no fractional variables\n");
      return SCIP_OKAY;
   }

   /* check, if the current NLP fractional variables array is invalid */
   if( nlp->validfracvars < stat->nnlps )
   {
      SCIP_VAR* var;
      SCIP_Real primsol;
      SCIP_Real frac;
      int branchpriority;
      int insertpos;
      int maxpriority;
      int i;

      SCIPsetDebugMsg(set, " -> recalculating NLP fractional variables\n");

      if( nlp->fracvarssize == 0 )
      {
         assert(nlp->fracvars     == NULL);
         assert(nlp->fracvarssol  == NULL);
         assert(nlp->fracvarsfrac == NULL);
         nlp->fracvarssize = 5;
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->fracvars,     nlp->fracvarssize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->fracvarssol,  nlp->fracvarssize) );
         SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &nlp->fracvarsfrac, nlp->fracvarssize) );
      }

      maxpriority = INT_MIN;
      nlp->nfracvars = 0;
      nlp->npriofracvars = 0;
      for( i = 0; i < nlp->nvars; ++i )
      {
         var = nlp->vars[i];
         assert(var != NULL);

         primsol = SCIPvarGetNLPSol(var);
         assert(primsol < SCIP_INVALID);

         /* consider only binary and integer variables */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY && SCIPvarGetType(var) != SCIP_VARTYPE_INTEGER )
            continue;

         /* ignore fixed variables (due to numerics, it is possible, that the NLP solution of a fixed integer variable
          * (with large fixed value) is fractional in terms of absolute feasibility measure)
          */
         if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
            continue;

         /* check, if the LP solution value is fractional */
         frac = SCIPsetFeasFrac(set, primsol);

         /* The fractionality should not be smaller than -feastol, however, if the primsol is large enough
          * and close to an integer, fixed precision floating point arithmetic might give us values slightly
          * smaller than -feastol. Originally, the "frac >= -feastol"-check was within SCIPsetIsFeasFracIntegral(),
          * however, we relaxed it to "frac >= -2*feastol" and have the stricter check here for small-enough primsols.
          */
         assert(SCIPsetIsGE(set, frac, -SCIPsetFeastol(set)) || (primsol > 1e14 * SCIPsetFeastol(set)));

         if( SCIPsetIsFeasFracIntegral(set, frac) )
            continue;

         /* ensure enough space in fracvars arrays */
         if( nlp->fracvarssize <= nlp->nfracvars )
         {
            int newsize;

            newsize = SCIPsetCalcMemGrowSize(set, nlp->nfracvars + 1);
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->fracvars,     nlp->fracvarssize, newsize) );
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->fracvarssol,  nlp->fracvarssize, newsize) );
            SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->fracvarsfrac, nlp->fracvarssize, newsize) );
            nlp->fracvarssize = newsize;
         }
         assert(nlp->nfracvars < nlp->fracvarssize);
         assert(nlp->fracvars     != NULL);
         assert(nlp->fracvarssol  != NULL);
         assert(nlp->fracvarsfrac != NULL);

         /* insert candidate in candidate list */
         branchpriority = SCIPvarGetBranchPriority(var);
         insertpos = nlp->nfracvars;
         nlp->nfracvars++;
         if( branchpriority > maxpriority )
         {
            /* candidate has higher priority than the current maximum:
             * move it to the front and declare it to be the single best candidate
             */
            if( insertpos != 0 )
            {
               nlp->fracvars[insertpos]     = nlp->fracvars[0];
               nlp->fracvarssol[insertpos]  = nlp->fracvarssol[0];
               nlp->fracvarsfrac[insertpos] = nlp->fracvarsfrac[0];
               insertpos = 0;
            }
            nlp->npriofracvars = 1;
            maxpriority = branchpriority;
         }
         else if( branchpriority == maxpriority )
         {
            /* candidate has equal priority as the current maximum:
             * move away the first non-maximal priority candidate, move the current candidate to the correct
             * slot (binaries first) and increase the number of maximal priority candidates
             */
            if( insertpos != nlp->npriofracvars )
            {
               nlp->fracvars[insertpos]     = nlp->fracvars[nlp->npriofracvars];
               nlp->fracvarssol[insertpos]  = nlp->fracvarssol[nlp->npriofracvars];
               nlp->fracvarsfrac[insertpos] = nlp->fracvarsfrac[nlp->npriofracvars];
               insertpos = nlp->npriofracvars;
            }
            ++nlp->npriofracvars;
         }
         nlp->fracvars[insertpos]     = var;
         nlp->fracvarssol[insertpos]  = primsol;
         nlp->fracvarsfrac[insertpos] = frac;

         SCIPsetDebugMsg(set, " -> candidate %d: var=<%s>, sol=%g, frac=%g, prio=%d (max: %d) -> pos %d\n",
            nlp->nfracvars, SCIPvarGetName(var), primsol, frac, branchpriority, maxpriority, insertpos);
      }

      nlp->validfracvars = stat->nnlps;
   }
   assert(0 <= nlp->npriofracvars);
   assert(nlp->npriofracvars <= nlp->nfracvars);

   SCIPsetDebugMsg(set, " -> %d fractional variables (%d of maximal priority)\n", nlp->nfracvars, nlp->npriofracvars);

   return SCIP_OKAY;
}

/** event handling for variable events */
static
SCIP_DECL_EVENTEXEC(eventExecNlp)
{
   SCIP_EVENTTYPE etype;
   SCIP_VAR*      var;

   assert(scip      != NULL);
   assert(eventhdlr != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);

   assert((SCIP_NLP*)eventdata == scip->nlp);

   etype = SCIPeventGetType(event);
   var   = SCIPeventGetVar(event);

   if( SCIP_EVENTTYPE_VARADDED & etype )
   {
      SCIPdebugMessage("-> handling varadd event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( SCIPnlpAddVar(scip->nlp, SCIPblkmem(scip), scip->set, var) );
   }
   else if( SCIP_EVENTTYPE_VARDELETED & etype )
   {
      SCIPdebugMessage("-> handling vardel event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( SCIPnlpDelVar(scip->nlp, SCIPblkmem(scip), scip->set, scip->eventqueue, scip->lp, var) );
   }
   else if( SCIP_EVENTTYPE_VARFIXED & etype )
   {
      /* variable was fixed, aggregated, or multi-aggregated */
      SCIPdebugMessage("-> handling variable fixation event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( nlpRemoveFixedVar(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, scip->eventqueue, scip->lp, var) );
   }
   else if( SCIP_EVENTTYPE_BOUNDCHANGED & etype )
   {
      SCIPdebugMessage("-> handling bound changed event %" SCIP_EVENTTYPE_FORMAT ", variable <%s>\n", etype, SCIPvarGetName(var) );
      SCIP_CALL( nlpUpdateVarBounds(scip->nlp, scip->set, var, (SCIP_Bool)(SCIP_EVENTTYPE_BOUNDTIGHTENED & etype)) );
   }
   else if( SCIP_EVENTTYPE_OBJCHANGED & etype )
   {
      SCIPdebugMessage("-> handling objchg event, variable <%s>\n", SCIPvarGetName(var) );
      SCIP_CALL( nlpUpdateObjCoef(scip->nlp, var) );
   }
   else
   {
      SCIPerrorMessage("unexpected event %d on variable <%s>\n", etype, SCIPvarGetName(var) );
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/*
 * public NLP methods
 */

/** includes NLP specific plugins (e.g., event handler) and parameters */
SCIP_RETCODE SCIPnlpInclude(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(set != NULL);
   assert(set->stage == SCIP_STAGE_INIT);

   /* check whether event handler is already present */
   if( SCIPsetFindEventhdlr(set, EVENTHDLR_NAME) != NULL )
   {
      SCIPerrorMessage("event handler <" EVENTHDLR_NAME "> already included.\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, eventExecNlp, NULL) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(set, eventhdlr) );

   return SCIP_OKAY;
} /*lint !e715*/

/** construct a new empty NLP */
SCIP_RETCODE SCIPnlpCreate(
   SCIP_NLP**            nlp,                /**< NLP handler, call by reference */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   const char*           name,               /**< problem name */
   int                   nvars_estimate      /**< an estimate on the number of variables that may be added to the NLP later */
   )
{
   assert(nlp  != NULL);
   assert(blkmem != NULL);
   assert(set  != NULL);
   assert(stat != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocMemory(nlp) );

   /* select NLP solver (if any available) and setup problem */
   if( set->nnlpis > 0 )
   {
      assert(set->nlp_solver != NULL);
      if( set->nlp_solver[0] == '\0' )
      { /* take solver with highest priority */
         assert(set->nlpis != NULL);

         /* sort the NLPIs if necessary */
         if( !set->nlpissorted )
            SCIPsetSortNlpis(set);

         (*nlp)->solver = set->nlpis[0];
      }
      else
      { /* find user specified NLP solver */
         (*nlp)->solver = SCIPsetFindNlpi(set, set->nlp_solver);
         if( (*nlp)->solver == NULL )
         {
            SCIPerrorMessage("Selected NLP solver <%s> not available.\n", set->nlp_solver);
            return SCIP_PLUGINNOTFOUND;
         }
      }
      assert((*nlp)->solver != NULL);
      SCIP_CALL( SCIPnlpiCreateProblem((*nlp)->solver, &(*nlp)->problem, "scip_nlp") );
   }
   else
   {
      /* maybe someone wanna use the NLP just to collect nonlinearities, but is not necessarily interesting on solving
       * so we allow this and just continue */
      (*nlp)->solver = NULL;
      (*nlp)->problem = NULL;
   }

   /* status */
   (*nlp)->nunflushedvaradd   = 0;
   (*nlp)->nunflushedvardel   = 0;
   (*nlp)->nunflushednlrowadd = 0;
   (*nlp)->nunflushednlrowdel = 0;
   (*nlp)->isrelax    = TRUE;
   (*nlp)->indiving   = FALSE;

   /* variables in problem and NLPI problem */
   (*nlp)->nvars = 0;
   (*nlp)->sizevars = 0;
   (*nlp)->vars = NULL;
   SCIP_CALL( SCIPhashmapCreate(&(*nlp)->varhash, blkmem, nvars_estimate) );

   (*nlp)->nvars_solver = 0;
   (*nlp)->sizevars_solver = 0;
   (*nlp)->varmap_nlp2nlpi = NULL;
   (*nlp)->varmap_nlpi2nlp = NULL;

   /* nonlinear rows in problem and NLPI problem */
   (*nlp)->nnlrows = 0;
   (*nlp)->sizenlrows = 0;
   (*nlp)->nlrows = NULL;

   (*nlp)->nnlrows_solver = 0;
   (*nlp)->sizenlrows_solver = 0;
   (*nlp)->nlrowmap_nlpi2nlp = NULL;

   /* objective function */
   (*nlp)->objflushed = TRUE;
   (*nlp)->divingobj = NULL;

   /* initial guess */
   (*nlp)->haveinitguess = FALSE;
   (*nlp)->initialguess = NULL;

   /* solution of NLP */
   (*nlp)->primalsolobjval = SCIP_INVALID;
   (*nlp)->solstat         = SCIP_NLPSOLSTAT_UNKNOWN;
   (*nlp)->termstat        = SCIP_NLPTERMSTAT_OTHER;
   (*nlp)->varlbdualvals   = NULL;
   (*nlp)->varubdualvals   = NULL;

   /* event handling: catch variable addition and deletion events */
   (*nlp)->eventhdlr = SCIPsetFindEventhdlr(set, EVENTHDLR_NAME);
   if( (*nlp)->eventhdlr == NULL )
   {
      SCIPerrorMessage("NLP eventhandler <" EVENTHDLR_NAME "> not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   SCIP_CALL( SCIPeventfilterAdd(set->scip->eventfilter, blkmem, set,
         SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_VARDELETED,
         (*nlp)->eventhdlr, (SCIP_EVENTDATA*)(*nlp), &(*nlp)->globalfilterpos) );

   /* fractional variables in last NLP solution */
   (*nlp)->fracvars     = NULL;
   (*nlp)->fracvarssol  = NULL;
   (*nlp)->fracvarsfrac = NULL;
   (*nlp)->nfracvars     = 0;
   (*nlp)->npriofracvars = 0;
   (*nlp)->fracvarssize  = 0;
   (*nlp)->validfracvars = -1;

   /* miscellaneous */
   SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*nlp)->name, name, strlen(name)+1) );

   return SCIP_OKAY;
}

/** frees NLP data object */
SCIP_RETCODE SCIPnlpFree(
   SCIP_NLP**            nlp,                /**< pointer to NLP data object */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   )
{
   assert(nlp    != NULL);
   assert(*nlp   != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   /* drop fractional variables */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->fracvars,     (*nlp)->fracvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->fracvarssol,  (*nlp)->fracvarssize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->fracvarsfrac, (*nlp)->fracvarssize);

   /* drop global events (variable addition and deletion) */
   SCIP_CALL( SCIPeventfilterDel(set->scip->eventfilter, blkmem, set,
         SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_VARDELETED,
         (*nlp)->eventhdlr, (SCIP_EVENTDATA*)(*nlp), (*nlp)->globalfilterpos) );

   SCIP_CALL( SCIPnlpReset(*nlp, blkmem, set, eventqueue, lp) );
   assert((*nlp)->nnlrows == 0);
   assert((*nlp)->nnlrows_solver == 0);
   assert((*nlp)->nvars == 0);
   assert((*nlp)->nvars_solver == 0);
   assert((*nlp)->initialguess == NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*nlp)->name, strlen((*nlp)->name)+1);

   /* free nonlinear rows arrays */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->nlrowmap_nlpi2nlp, (*nlp)->sizenlrows_solver);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->nlrows, (*nlp)->sizenlrows);

   /* free variables arrays */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varmap_nlp2nlpi, (*nlp)->sizevars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varmap_nlpi2nlp, (*nlp)->sizevars_solver);
   SCIPhashmapFree(&(*nlp)->varhash);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->vars, (*nlp)->sizevars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varlbdualvals, (*nlp)->sizevars);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*nlp)->varubdualvals, (*nlp)->sizevars);

   /* free NLPI problem */
   if( (*nlp)->problem != NULL )
   {
      SCIP_CALL( SCIPnlpiFreeProblem((*nlp)->solver, &(*nlp)->problem) );
   }

   /* free NLP data structure */
   BMSfreeMemory(nlp);

   return SCIP_OKAY;
}

/** resets the NLP to the empty NLP by removing all variables and rows from NLP,
 *  releasing all rows, and flushing the changes to the NLP solver
 */
SCIP_RETCODE SCIPnlpReset(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< SCIP LP, needed for releasing variables */
   )
{
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   if( nlp->indiving )
   {
      SCIP_CALL( SCIPnlpEndDive(nlp, blkmem, set) );
   }

   nlp->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   nlp->termstat = SCIP_NLPTERMSTAT_OTHER;

   BMSfreeBlockMemoryArrayNull(blkmem, &nlp->initialguess, nlp->nvars);
   nlp->haveinitguess = FALSE;

   for(i = nlp->nnlrows - 1; i >= 0; --i)
   {
      SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, i) );
   }

   for(i = nlp->nvars - 1; i >= 0; --i)
   {
      SCIP_CALL( nlpDelVarPos(nlp, blkmem, set, eventqueue, lp, i) );
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set) );

   return SCIP_OKAY;
}

/** currently a dummy function that always returns TRUE */
SCIP_Bool SCIPnlpHasCurrentNodeNLP(
   SCIP_NLP*             nlp                 /**< NLP data */
   )
{
   return TRUE;
} /*lint !e715*/

/** ensures, that variables array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureVarsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nvars <= nlp->sizevars);

   if( num > nlp->sizevars )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->vars,            nlp->sizevars, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varmap_nlp2nlpi, nlp->sizevars, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varlbdualvals,   nlp->sizevars, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->varubdualvals,   nlp->sizevars, newsize) );
      if( nlp->initialguess != NULL )
      {
         SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->initialguess, nlp->sizevars, newsize) );
      }

      nlp->sizevars = newsize;
   }
   assert(num <= nlp->sizevars);

   return SCIP_OKAY;
}

/** adds a variable to the NLP and captures the variable */
SCIP_RETCODE SCIPnlpAddVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPvarIsTransformed(var));
   assert(!SCIPhashmapExists(nlp->varhash, var));

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add variable during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddVars(nlp, blkmem, set, 1, &var) );

   return SCIP_OKAY;
}

/** adds a set of variables to the NLP and captures the variables */
SCIP_RETCODE SCIPnlpAddVars(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables to add */
   SCIP_VAR**            vars                /**< variables to add */
   )
{
   assert(nlp != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(vars != NULL || nvars == 0);

   if( nlp->indiving && nvars > 0)
   {
      SCIPerrorMessage("cannot add variables during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddVars(nlp, blkmem, set, nvars, vars) );

   return SCIP_OKAY;
}

/** deletes a variable from the NLP and releases the variable */
SCIP_RETCODE SCIPnlpDelVar(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< SCIP LP, needed to release variable */
   SCIP_VAR*             var                 /**< variable */
   )
{
   int varpos;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(var    != NULL);

   if( !SCIPhashmapExists(nlp->varhash, var) )
   {
      SCIPerrorMessage("variable <%s> not found in NLP, cannot delete\n", SCIPvarGetName(var));
      return SCIP_ERROR;
   }

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot delete variable during NLP diving\n");
      return SCIP_ERROR;
   }

   varpos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);

   SCIP_CALL( nlpDelVarPos(nlp, blkmem, set, eventqueue, lp, varpos) );

   return SCIP_OKAY;
}

/** ensures, that nonlinear rows array of NLP can store at least num entries */
SCIP_RETCODE SCIPnlpEnsureNlRowsSize(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlp->nnlrows <= nlp->sizenlrows);

   if( num > nlp->sizenlrows )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &nlp->nlrows, nlp->sizenlrows, newsize) );

      nlp->sizenlrows = newsize;
   }
   assert(num <= nlp->sizenlrows);

   return SCIP_OKAY;
}

/** adds a nonlinear row to the NLP and captures it
 * all variables of the row need to be present in the NLP */
SCIP_RETCODE SCIPnlpAddNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   assert(nlp   != NULL);
   assert(nlrow != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add row during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddNlRows(nlp, blkmem, set, stat, 1, &nlrow) );

   return SCIP_OKAY;
}

/** adds nonlinear rows to the NLP and captures them
 * all variables of the row need to be present in the NLP */
SCIP_RETCODE SCIPnlpAddNlRows(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   int                   nnlrows,            /**< number of rows to add */
   SCIP_NLROW**          nlrows              /**< rows to add */
   )
{
   assert(nlp    != NULL);
   assert(nlrows != NULL || nnlrows == 0);

   if( nnlrows == 0 )
      return SCIP_OKAY;

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot add rows during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpAddNlRows(nlp, blkmem, set, stat, nnlrows, nlrows) );

   return SCIP_OKAY;
}

/** deletes a nonlinear row from the NLP
 * does nothing if nonlinear row is not in NLP */
SCIP_RETCODE SCIPnlpDelNlRow(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(nlrow  != NULL);

   /* if row not in NLP, nothing to do */
   if( nlrow->nlpindex == -1 )
      return SCIP_OKAY;

   assert(nlrow->nlpindex >= 0);
   assert(nlrow->nlpindex < nlp->nnlrows);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot delete row during NLP diving\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, nlrow->nlpindex) );

   return SCIP_OKAY;
}

/** applies all cached changes to the NLP solver */
SCIP_RETCODE SCIPnlpFlush(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot flush NLP during NLP diving\n");
      return SCIP_ERROR;
   }

   /* flush removals of nonlinear rows and variables */
   SCIP_CALL( nlpFlushNlRowDeletions(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushVarDeletions(nlp, blkmem, set) );
   assert(nlp->nunflushednlrowdel == 0);
   assert(nlp->nunflushedvardel   == 0);

   /* flush addition of variables, objective, and addition of rows */
   SCIP_CALL( nlpFlushVarAdditions(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushObjective(nlp, blkmem, set) );
   SCIP_CALL( nlpFlushNlRowAdditions(nlp, blkmem, set) );
   assert(nlp->nunflushedvaradd == 0);
   assert(nlp->objflushed == TRUE);
   assert(nlp->nunflushednlrowadd == 0);

   assert(nlp->nvars   == nlp->nvars_solver);
   assert(nlp->nnlrows == nlp->nnlrows_solver);

   return SCIP_OKAY;
}

/** solves the NLP */
SCIP_RETCODE SCIPnlpSolve(
   SCIP_NLP*             nlp,                /**< NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot solve NLP during NLP diving (use SCIPsolveDiveNLP)\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set) );

   SCIP_CALL( nlpSolve(nlp, blkmem, set, messagehdlr, stat) );

   return SCIP_OKAY;
}

/** gets objective value of current NLP */
SCIP_Real SCIPnlpGetObjval(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->primalsolobjval;
}

/** gives current pseudo objective value */
SCIP_RETCODE SCIPnlpGetPseudoObjval(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real*            pseudoobjval        /**< buffer to store pseudo objective value */
   )
{
   assert(nlp != NULL);
   assert(pseudoobjval != NULL);

   if( nlp->divingobj != NULL )
   {
      assert(nlp->indiving);
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlp->divingobj, set, stat, pseudoobjval) );
   }
   else
   {
      int i;

      *pseudoobjval = 0.0;
      for( i = 0; i < nlp->nvars; ++i )
         *pseudoobjval += SCIPvarGetObj(nlp->vars[i]) * SCIPvarGetBestBoundLocal(nlp->vars[i]);
   }

   return SCIP_OKAY;
}

/** gets fractional variables of last NLP solution along with solution values and fractionalities
 */
SCIP_RETCODE SCIPnlpGetFracVars(
   SCIP_NLP*             nlp,                /**< NLP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR***           fracvars,           /**< pointer to store the array of NLP fractional variables, or NULL */
   SCIP_Real**           fracvarssol,        /**< pointer to store the array of NLP fractional variables solution values, or NULL */
   SCIP_Real**           fracvarsfrac,       /**< pointer to store the array of NLP fractional variables fractionalities, or NULL */
   int*                  nfracvars,          /**< pointer to store the number of NLP fractional variables , or NULL */
   int*                  npriofracvars       /**< pointer to store the number of NLP fractional variables with maximal branching priority, or NULL */
   )
{
   assert(nlp != NULL);

   SCIP_CALL( nlpCalcFracVars(nlp, blkmem, set, stat) );
   assert(nlp->fracvars     != NULL);
   assert(nlp->fracvarssol  != NULL);
   assert(nlp->fracvarsfrac != NULL);

   if( fracvars != NULL )
      *fracvars = nlp->fracvars;
   if( fracvarssol != NULL )
      *fracvarssol = nlp->fracvarssol;
   if( fracvarsfrac != NULL )
      *fracvarsfrac = nlp->fracvarsfrac;
   if( nfracvars != NULL )
      *nfracvars = nlp->nfracvars;
   if( npriofracvars != NULL )
      *npriofracvars = nlp->npriofracvars;

   return SCIP_OKAY;
}

/** removes all redundant nonlinear rows */
SCIP_RETCODE SCIPnlpRemoveRedundantNlRows(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_NLPSOLSTAT solstatus;
   SCIP_Bool isredundant;
   int i;

   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(set    != NULL);
   assert(stat   != NULL);

   if( nlp->nnlrows == 0 )
      return SCIP_OKAY;

   if( nlp->indiving )
   {
      SCIPerrorMessage("cannot remove redundant rows during NLP diving\n");
      return SCIP_ERROR;
   }

   /* removing redundant rows should not change the solution status, so we reset it at the end */
   solstatus = nlp->solstat;

   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIP_CALL( SCIPnlrowIsRedundant(nlp->nlrows[i], set, stat, &isredundant) );
      if( isredundant )
      {
         SCIP_CALL( nlpDelNlRowPos(nlp, blkmem, set, i) );
      }
   }

   nlp->solstat = solstatus;

   return SCIP_OKAY;
}

/** set initial guess (approximate primal solution) for next solve
 *
 *  array initguess must be NULL or have length at least SCIPnlpGetNVars()
 */
SCIP_RETCODE SCIPnlpSetInitialGuess(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_Real*            initguess           /**< new initial guess, or NULL to clear previous one */
   )
{
   assert(nlp    != NULL);
   assert(blkmem != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   /* if user wants to let NLP solver choose start point, then invalidate current initial guess both in NLP and in NLPI */
   if( initguess == NULL )
   {
      nlp->haveinitguess = FALSE;
      SCIP_CALL( SCIPnlpiSetInitialGuess(nlp->solver, nlp->problem, NULL, NULL, NULL, NULL) );
      return SCIP_OKAY;
   }

   if( nlp->initialguess != NULL )
   {
      BMScopyMemoryArray(nlp->initialguess, initguess, nlp->nvars);
   }
   else
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &nlp->initialguess, initguess, nlp->nvars) );
   }
   nlp->haveinitguess = TRUE;

   return SCIP_OKAY;
}

/** writes NLP to a file */
SCIP_RETCODE SCIPnlpWrite(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           fname               /**< file name */
   )
{
   FILE* file;
   int i;

   assert(nlp != NULL);

   if( fname != NULL )
   {
      file = fopen(fname, "w");
      if( file == NULL )
      {
         SCIPerrorMessage("could not open file <%s> for writing\n", fname);
         return SCIP_FILECREATEERROR;
      }
   }
   else
      file = stdout;

   SCIPmessageFPrintInfo(messagehdlr, file, "STATISTICS\n");
   SCIPmessageFPrintInfo(messagehdlr, file, "  NLP name: %s\n", nlp->name);
   SCIPmessageFPrintInfo(messagehdlr, file, "  Variables: %d\n", nlp->nvars);
   SCIPmessageFPrintInfo(messagehdlr, file, "  Rows: %d\n", nlp->nnlrows);

   SCIPmessageFPrintInfo(messagehdlr, file, "VARIABLES\n");
   for( i = 0; i < nlp->nvars; ++i )
   {
      SCIP_CALL( SCIPvarPrint(nlp->vars[i], set, messagehdlr, file) );
   }

   SCIPmessageFPrintInfo(messagehdlr, file, "NONLINEAR ROWS\n");
   for( i = 0; i < nlp->nnlrows; ++i )
   {
      SCIPmessageFPrintInfo(messagehdlr, file, "  ");
      SCIP_CALL( SCIPnlrowPrint(nlp->nlrows[i], messagehdlr, file) );
   }

   if( fname != NULL )
   {
      fclose(file);
   }

   return SCIP_OKAY;
}

/** gets array with variables of the NLP */
SCIP_VAR** SCIPnlpGetVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->vars;
}

/** gets current number of variables in NLP */
int SCIPnlpGetNVars(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nvars;
}

/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var */
SCIP_RETCODE SCIPnlpGetVarsNonlinearity(
   SCIP_NLP*             nlp,                /**< current NLP data */
   int*                  nlcount             /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   )
{
   SCIP_NLROW* nlrow;
   int varidx;
   int i;
   int c;

   assert(nlp != NULL);
   assert(nlcount != NULL || nlp->nvars == 0);

   BMSclearMemoryArray(nlcount, nlp->nvars);

   for( c = 0; c < nlp->nnlrows; ++c )
   {
      nlrow = nlp->nlrows[c];
      assert(nlrow != NULL);

      for( i = 0; i < nlrow->nquadvars; ++i )
      {
         assert(SCIPhashmapExists(nlp->varhash, (void*)nlrow->quadvars[i]));
         varidx = (int)(size_t) SCIPhashmapGetImage(nlp->varhash, (void*)nlrow->quadvars[i]);
         assert(varidx < nlp->nvars);
         ++nlcount[varidx];  /*lint !e613 */
      }

      if( nlrow->exprtree != NULL )
      {
         SCIP_VAR** exprtreevars;
         int nexprtreevars;

         exprtreevars = SCIPexprtreeGetVars(nlrow->exprtree);
         nexprtreevars = SCIPexprtreeGetNVars(nlrow->exprtree);
         assert(exprtreevars != NULL || nexprtreevars == 0);
         for( i = 0; i < nexprtreevars; ++i )
         {
            assert(SCIPhashmapExists(nlp->varhash, (void*)exprtreevars[i]));  /*lint !e613 */

            /* skip variables that also appear in quadratic part, so they are not counted twice */
            if( nlrow->quadvarshash != NULL && SCIPhashmapExists(nlrow->quadvarshash, (void*)exprtreevars[i]) )  /*lint !e613 */
               continue;

            varidx = (int)(size_t) SCIPhashmapGetImage(nlp->varhash, (void*)exprtreevars[i]);  /*lint !e613 */
            assert(varidx < nlp->nvars);
            ++nlcount[varidx];  /*lint !e613 */
         }
      }
   }

   return SCIP_OKAY;
}


/** indicates whether there exists a row that contains a continuous variable in a nonlinear term
 *
 * @note The method may have to touch every row and nonlinear term to compute its result.
 */
SCIP_Bool SCIPnlpHasContinuousNonlinearity(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   SCIP_NLROW* nlrow;
   int c;
   int i;

   assert(nlp != NULL);

   for( c = 0; c < nlp->nnlrows; ++c )
   {
      nlrow = nlp->nlrows[c];
      assert(nlrow != NULL);

      for( i = 0; i < nlrow->nquadvars; ++i )
         if( SCIPvarGetType(nlrow->quadvars[i]) == SCIP_VARTYPE_CONTINUOUS )
            return TRUE;

      if( nlrow->exprtree != NULL )
      {
         SCIP_VAR** exprtreevars;
         int nexprtreevars;

         exprtreevars = SCIPexprtreeGetVars(nlrow->exprtree);
         nexprtreevars = SCIPexprtreeGetNVars(nlrow->exprtree);
         assert(exprtreevars != NULL || nexprtreevars == 0);

         for( i = 0; i < nexprtreevars; ++i )
            if( SCIPvarGetType(exprtreevars[i]) == SCIP_VARTYPE_CONTINUOUS ) /*lint !e613*/
               return TRUE;
      }
   }

   return FALSE;
}

/** gives dual solution values associated with lower bounds of NLP variables */
SCIP_Real* SCIPnlpGetVarsLbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->varlbdualvals;
}

/** gives dual solution values associated with upper bounds of NLP variables */
SCIP_Real* SCIPnlpGetVarsUbDualsol(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->varubdualvals;
}

/** gets array with nonlinear rows of the NLP */
SCIP_NLROW** SCIPnlpGetNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nlrows;
}

/** gets current number of nonlinear rows in NLP */
int SCIPnlpGetNNlRows(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->nnlrows;
}

/** gets the NLP solver interface */
SCIP_NLPI* SCIPnlpGetNLPI(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solver;
}

/** gets the NLP problem in the solver interface */
SCIP_NLPIPROBLEM* SCIPnlpGetNLPIProblem(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->problem;
}

/** indicates whether NLP is currently in diving mode */
SCIP_Bool SCIPnlpIsDiving(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->indiving;
}

/** gets solution status of current NLP */
SCIP_NLPSOLSTAT SCIPnlpGetSolstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solstat;
}

/** gets termination status of last NLP solve */
SCIP_NLPTERMSTAT SCIPnlpGetTermstat(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->termstat;
}

/** gives statistics (number of iterations, solving time, ...) of last NLP solve */
SCIP_RETCODE SCIPnlpGetStatistics(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   assert(nlp != NULL);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);
   assert(statistics != NULL);

   SCIP_CALL( SCIPnlpiGetStatistics(nlp->solver, nlp->problem, statistics) );

   return SCIP_OKAY;
}

/** indicates whether a feasible solution for the current NLP is available
 * thus, returns whether the solution status <= feasible  */
SCIP_Bool SCIPnlpHasSolution(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   assert(nlp != NULL);

   return nlp->solstat <= SCIP_NLPSOLSTAT_FEASIBLE;
}

/** gets integer parameter of NLP */
SCIP_RETCODE SCIPnlpGetIntPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);
   assert(ival != NULL);

   SCIP_CALL( SCIPnlpiGetIntPar(nlp->solver, nlp->problem, type, ival) );

   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_RETCODE SCIPnlpSetIntPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( SCIPnlpiSetIntPar(nlp->solver, nlp->problem, type, ival) );

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpGetRealPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);
   assert(dval != NULL);

   SCIP_CALL( SCIPnlpiGetRealPar(nlp->solver, nlp->problem, type, dval) );

   return SCIP_OKAY;
}

/** sets floating point parameter of NLP */
SCIP_RETCODE SCIPnlpSetRealPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( SCIPnlpiSetRealPar(nlp->solver, nlp->problem, type, dval) );

   return SCIP_OKAY;
}

/** gets string parameter of NLP */
SCIP_RETCODE SCIPnlpGetStringPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value */
   )
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);
   assert(sval != NULL);

   SCIP_CALL( SCIPnlpiGetStringPar(nlp->solver, nlp->problem, type, sval) );

   return SCIP_OKAY;
}

/** sets string parameter of NLP */
SCIP_RETCODE SCIPnlpSetStringPar(
   SCIP_NLP*             nlp,                /**< pointer to NLP datastructure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
   )
{
   assert(nlp  != NULL);
   assert(nlp->solver  != NULL);
   assert(nlp->problem != NULL);

   SCIP_CALL( SCIPnlpiSetStringPar(nlp->solver, nlp->problem, type, sval) );

   return SCIP_OKAY;
}

/*
 * NLP diving methods
 */

/** signals start of diving */
SCIP_RETCODE SCIPnlpStartDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(nlp != NULL);

   if( nlp->indiving )
   {
      SCIPerrorMessage("NLP is already in diving mode\n");
      return SCIP_ERROR;
   }

   if( nlp->solver == NULL )
   {
      /* In diving mode we do not cache changes but put them directly in the NLPI problem, which does not exist if there is no solver.
       * So we forbid diving of no solver is available. */
      SCIPerrorMessage("Cannot start diving if no NLP solver is available\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlpFlush(nlp, blkmem, set) );

   nlp->indiving = TRUE;

   return SCIP_OKAY;
}

/** resets the bound and objective changes made during diving and disables diving mode */
SCIP_RETCODE SCIPnlpEndDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   int* varidx;
   SCIP_Real* varlb;
   SCIP_Real* varub;

   assert(nlp != NULL);
   assert(set != NULL);
   assert(nlp->nvars == nlp->nvars_solver);

   if( !nlp->indiving )
   {
      SCIPerrorMessage("NLP not in diving mode, cannot end dive\n");
      return SCIP_ERROR;
   }

   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* reset variable bounds in NLPI problem to their current values */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varidx, nlp->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varlb,  nlp->nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varub,  nlp->nvars) );
   for( i = 0; i < nlp->nvars; ++i )
   {
      varidx[i] = i;
      varlb[i] = SCIPvarGetLbLocal(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
      varub[i] = SCIPvarGetUbLocal(nlp->vars[nlp->varmap_nlpi2nlp[i]]);
   }

   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, nlp->nvars, varidx, varlb, varub) );

   SCIPsetFreeBufferArray(set, &varidx);
   SCIPsetFreeBufferArray(set, &varlb);
   SCIPsetFreeBufferArray(set, &varub);

   /* clear diving objective, if one was used (i.e., if SCIPnlpChgVarObjDive had been called)
    * the objective in the NLPI will be reset in the next flush */
   if( nlp->divingobj != NULL )
   {
      SCIP_CALL( SCIPnlrowRelease(&nlp->divingobj, blkmem, set) );
      assert(nlp->divingobj == NULL);
      assert(nlp->objflushed == FALSE);
   }

   /* we do not have a valid solution anymore */
   nlp->solstat  = SCIP_NLPSOLSTAT_UNKNOWN;
   nlp->termstat = SCIP_NLPTERMSTAT_OTHER;
   nlp->primalsolobjval = SCIP_INVALID;

   nlp->indiving = FALSE;

   return SCIP_OKAY;
}

/** changes coefficient of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarObjDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new linear coefficient of variable in objective */
   )
{
   int pos;
   int objidx;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));
   assert(nlp->indiving);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* get position of variable in NLPI problem */
   pos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);
   pos = nlp->varmap_nlp2nlpi[pos];
   assert(pos >= 0);

   /* set coefficient in NLPI problem objective */
   objidx = -1;
   SCIP_CALL( SCIPnlpiChgLinearCoefs(nlp->solver, nlp->problem, objidx, 1, &pos, &coef) );

   /* create an nlrow that holds the diving objective, if not done yet */
   if( nlp->divingobj == NULL )
   {
      SCIP_Real* coefs;
      int        i;

      SCIP_CALL( SCIPsetAllocBufferArray(set, &coefs, nlp->nvars) );
      for( i = 0; i < nlp->nvars; ++i )
         coefs[i] = SCIPvarGetObj(nlp->vars[i]);

      SCIP_CALL( SCIPnlrowCreate(&nlp->divingobj, blkmem, set, "divingobj",
            0.0,
            nlp->nvars, nlp->vars, coefs,
            0, NULL, 0, NULL,
            NULL,
            -SCIPsetInfinity(set), SCIPsetInfinity(set),
            SCIP_EXPRCURV_LINEAR) );

      SCIPsetFreeBufferArray(set, &coefs);
   }
   assert(nlp->divingobj != NULL);

   /* modify coefficient in diving objective */
   SCIP_CALL( SCIPnlrowChgLinearCoef(nlp->divingobj, blkmem, set, stat, nlp, var, coef) );

   /* remember that we have to store objective after diving ended */
   nlp->objflushed = FALSE;

   return SCIP_OKAY;
}

/** changes bounds of variable in diving NLP */
SCIP_RETCODE SCIPnlpChgVarBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             lb,                 /**< new lower bound of variable */
   SCIP_Real             ub                  /**< new upper bound of variable */
   )
{
   int pos;

   assert(nlp != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(nlp->varhash, var));
   assert(nlp->indiving);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   /* get position of variable in NLPI problem */
   pos = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, var);
   pos = nlp->varmap_nlp2nlpi[pos];
   assert(pos >= 0);

   /* set new bounds in NLPI */
   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, 1, &pos, &lb, &ub) );

   return SCIP_OKAY;
}

/** changes bounds of a set of variables in diving NLP */
SCIP_RETCODE SCIPnlpChgVarsBoundsDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nvars,              /**< number of variables which bounds to change */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds of variables */
   SCIP_Real*            ubs                 /**< new upper bounds of variables */
   )
{
   int i;
   int* poss;

   assert(nlp  != NULL);
   assert(vars != NULL || nvars == 0);
   assert(nlp->indiving);
   assert(lbs  != NULL || nvars == 0);
   assert(ubs  != NULL || nvars == 0);
   assert(nlp->solver != NULL);
   assert(nlp->problem != NULL);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPsetAllocBufferArray(set, &poss, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      assert(SCIPhashmapExists(nlp->varhash, vars[i]));  /*lint !e613*/

      /* get position of variable in NLPI problem */
      poss[i] = (int) (size_t) SCIPhashmapGetImage(nlp->varhash, vars[i]);   /*lint !e613*/
      poss[i] = nlp->varmap_nlp2nlpi[poss[i]];
      assert(poss[i] >= 0);
   }

   /* set new bounds in NLPI */
   SCIP_CALL( SCIPnlpiChgVarBounds(nlp->solver, nlp->problem, nvars, poss, lbs, ubs) );

   SCIPsetFreeBufferArray(set, &poss);

   return SCIP_OKAY;
}

/** returns whether the objective function has been changed during diving */
SCIP_Bool SCIPnlpIsDivingObjChanged(
   SCIP_NLP*             nlp                 /**< current NLP data */
   )
{
   return nlp->divingobj != NULL;
}

/** solves diving NLP */
SCIP_RETCODE SCIPnlpSolveDive(
   SCIP_NLP*             nlp,                /**< current NLP data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_CALL( nlpSolve(nlp, blkmem, set, messagehdlr, stat) );

   return SCIP_OKAY;
}

