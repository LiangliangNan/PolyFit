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

/**@file   nlhdlr_perspective.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  perspective nonlinear handler
 * @author Ksenia Bestuzheva
 */

#include <string.h>

#include "scip/nlhdlr_perspective.h"
#include "scip/cons_nonlinear.h"
#include "scip/scip_sol.h"
#include "scip/pub_misc_rowprep.h"
#include "scip/nlhdlr.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "perspective"
#define NLHDLR_DESC               "perspective handler for expressions"
#define NLHDLR_DETECTPRIORITY     -20   /**< detect last so that to make use of what other handlers detected */
#define NLHDLR_ENFOPRIORITY       125   /**< enforce first because perspective cuts are always stronger */

#define DEFAULT_MAXPROPROUNDS     1     /**< maximal number of propagation rounds in probing */
#define DEFAULT_MINDOMREDUCTION   0.1   /**< minimal relative reduction in a variable's domain for applying probing */
#define DEFAULT_MINVIOLPROBING    1e-05 /**< minimal violation w.r.t. auxiliary variables for applying probing */
#define DEFAULT_PROBINGONLYINSEPA TRUE  /**< whether to do probing only in separation loop */
#define DEFAULT_PROBINGFREQ       1     /**< probing frequency (-1 - no probing, 0 - root node only) */
#define DEFAULT_CONVEXONLY        FALSE /**< whether perspective cuts are added only for convex expressions */
#define DEFAULT_TIGHTENBOUNDS     TRUE  /**< whether variable semicontinuity is used to tighten variable bounds */
#define DEFAULT_ADJREFPOINT       TRUE  /**< whether to adjust the reference point if indicator is not 1 */

/*
 * Data structures
 */

/** data structure to store information of a semicontinuous variable
 *
 * For a variable x (not stored in the struct), this stores the data of nbnds implications
 *   bvars[i] = 0 -> x = vals[i]
 *   bvars[i] = 1 -> lbs[i] <= x <= ubs[i]
 * where bvars[i] are binary variables.
 */
struct SCVarData
{
   SCIP_Real*            vals0;              /**< values of the variable when the corresponding bvars[i] = 0 */
   SCIP_Real*            lbs1;               /**< global lower bounds of the variable when the corresponding bvars[i] = 1 */
   SCIP_Real*            ubs1;               /**< global upper bounds of the variable when the corresponding bvars[i] = 1 */
   SCIP_VAR**            bvars;              /**< the binary variables on which the variable domain depends */
   int                   nbnds;              /**< number of suitable on/off bounds the var has */
   int                   bndssize;           /**< size of the arrays */
};
typedef struct SCVarData SCVARDATA;

/** nonlinear handler expression data
 *
 * For an expression expr (not stored in the struct), this stores the data of nindicators implications
 *   indicators[i] = 0 -> expr = exprvals[0]
 * where indicators[i] is an indicator (binary) variable, corresponding to some bvars entry in SCVarData.
 *
 * Also stores the variables the expression depends on.
 */
struct SCIP_NlhdlrExprData
{
   SCIP_Real*            exprvals0;          /**< 'off' values of the expression for each indicator variable */
   SCIP_VAR**            vars;               /**< expression variables (both original and auxiliary) */
   int                   nvars;              /**< total number of variables in the expression */
   int                   varssize;           /**< size of the vars array */
   SCIP_VAR**            indicators;         /**< all indicator variables for the expression */
   int                   nindicators;        /**< number of indicator variables */
};

/** nonlinear handler data */
struct SCIP_NlhdlrData
{
   SCIP_HASHMAP*         scvars;             /**< maps semicontinuous variables to their on/off bounds (SCVarData) */

   /* parameters */
   int                   maxproprounds;      /**< maximal number of propagation rounds in probing */
   SCIP_Real             mindomreduction;    /**< minimal relative reduction in a variable's domain for applying probing */
   SCIP_Real             minviolprobing;     /**< minimal violation w.r.t. auxiliary variables for applying probing */
   SCIP_Bool             probingonlyinsepa;  /**< whether to do probing only in separation loop */
   int                   probingfreq;        /**< if and when to do probing */
   SCIP_Bool             convexonly;         /**< whether perspective cuts are added only for convex expressions */
   SCIP_Bool             tightenbounds;      /**< whether variable semicontinuity is used to tighten variable bounds */
   SCIP_Bool             adjrefpoint;        /**< whether to adjust the reference point if indicator is not 1 */
};

/*
 * Local methods
 */

/*
 * Helper methods for working with nlhdlrExprData
 */

/** frees nlhdlrexprdata structure */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nlhdlr expression data */
   )
{
   int v;

   if( nlhdlrexprdata->nindicators != 0 )
   {
      assert(nlhdlrexprdata->indicators != NULL);
      for( v = nlhdlrexprdata->nindicators - 1; v >= 0; --v )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(nlhdlrexprdata->indicators[v])) );
      }
      SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->indicators), nlhdlrexprdata->nindicators);
      SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->exprvals0), nlhdlrexprdata->nindicators);
   }

   for( v = nlhdlrexprdata->nvars - 1; v >= 0; --v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(nlhdlrexprdata->vars[v])) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &nlhdlrexprdata->vars, nlhdlrexprdata->varssize);

   return SCIP_OKAY;
}

/* remove an indicator from nlhdlr expression data */
static
SCIP_RETCODE removeIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlexprdata,         /**< nlhdlr expression data */
   int                   pos                 /**< position of the indicator */
   )
{
   int i;

   assert(pos >= 0 && pos < nlexprdata->nindicators);

   SCIP_CALL( SCIPreleaseVar(scip, &nlexprdata->indicators[pos]) );
   for( i = pos; i < nlexprdata->nindicators - 1; ++i )
   {
      nlexprdata->indicators[i] = nlexprdata->indicators[i+1];
   }

   --nlexprdata->nindicators;

   return SCIP_OKAY;
}

/** adds an auxiliary variable to the vars array in nlhdlrexprdata */
static
SCIP_RETCODE addAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_HASHMAP*         auxvarmap,          /**< hashmap linking auxvars to positions in nlhdlrexprdata->vars */
   SCIP_VAR*             auxvar              /**< variable to be added */
   )
{
   int pos;
   int newsize;

   assert(nlhdlrexprdata != NULL);
   assert(auxvar != NULL);

   pos = SCIPhashmapGetImageInt(auxvarmap, (void*) auxvar);

   if( pos != INT_MAX )
      return SCIP_OKAY;

   /* ensure size */
   if( nlhdlrexprdata->nvars + 1 > nlhdlrexprdata->varssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, nlhdlrexprdata->nvars + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->vars, nlhdlrexprdata->varssize, newsize) );
      nlhdlrexprdata->varssize = newsize;
   }
   assert(nlhdlrexprdata->nvars + 1 <= nlhdlrexprdata->varssize);

   nlhdlrexprdata->vars[nlhdlrexprdata->nvars] = auxvar;
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPhashmapSetImageInt(auxvarmap, (void*) auxvar, nlhdlrexprdata->nvars) );
   ++(nlhdlrexprdata->nvars);

   return SCIP_OKAY;
}

/*
 * Semicontinuous variable methods
 */

/** adds an indicator to the data of a semicontinuous variable */
static
SCIP_RETCODE addSCVarIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCVARDATA*            scvdata,            /**< semicontinuous variable data */
   SCIP_VAR*             indicator,          /**< indicator to be added */
   SCIP_Real             val0,               /**< value of the variable when indicator == 0 */
   SCIP_Real             lb1,                /**< lower bound of the variable when indicator == 1 */
   SCIP_Real             ub1                 /**< upper bound of the variable when indicator == 1 */
   )
{
   int newsize;
   int i;
   SCIP_Bool found;
   int pos;

   assert(scvdata != NULL);
   assert(indicator != NULL);

   /* find the position where to insert */
   if( scvdata->bvars == NULL )
   {
      assert(scvdata->nbnds == 0 && scvdata->bndssize == 0);
      found = FALSE;
      pos = 0;
   }
   else
   {
      found = SCIPsortedvecFindPtr((void**)scvdata->bvars, SCIPvarComp, (void*)indicator, scvdata->nbnds, &pos);
   }

   if( found )
      return SCIP_OKAY;

   /* ensure sizes */
   if( scvdata->nbnds + 1 > scvdata->bndssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, scvdata->nbnds + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->bvars, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->lbs1, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->ubs1, scvdata->bndssize, newsize) );
      scvdata->bndssize = newsize;
   }
   assert(scvdata->nbnds + 1 <= scvdata->bndssize);
   assert(scvdata->bvars != NULL);

   /* move entries if needed */
   for( i = scvdata->nbnds; i > pos; --i )
   {
      /* coverity[forward_null] */
      scvdata->bvars[i] = scvdata->bvars[i-1];
      scvdata->vals0[i] = scvdata->vals0[i-1];
      scvdata->lbs1[i] = scvdata->lbs1[i-1];
      scvdata->ubs1[i] = scvdata->ubs1[i-1];
   }

   scvdata->bvars[pos] = indicator;
   scvdata->vals0[pos] = val0;
   scvdata->lbs1[pos] = lb1;
   scvdata->ubs1[pos] = ub1;
   ++scvdata->nbnds;

   return SCIP_OKAY;
}

/** find scvardata of var and position of indicator in it
 *
 *  If indicator is not there, returns NULL.
 */
static
SCVARDATA* getSCVarDataInd(
   SCIP_HASHMAP*         scvars,             /**< hashmap linking variables to scvardata */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             indicator,          /**< indicator variable */
   int*                  pos                 /**< pointer to store the position of indicator */
   )
{
   SCIP_Bool exists;
   SCVARDATA* scvdata;

   assert(var != NULL);
   assert(scvars != NULL);
   assert(indicator != NULL);

   scvdata = (SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)var);
   if( scvdata != NULL )
   {
      /* look for the indicator variable */
      exists = SCIPsortedvecFindPtr((void**)scvdata->bvars, SCIPvarComp, (void*)indicator, scvdata->nbnds, pos);
      if( !exists )
         return NULL;

      return scvdata;
   }

   return NULL;
}

/** checks if a variable is semicontinuous and, if needed, updates the scvars hashmap
 *
 * A variable \f$x\f$ is semicontinuous if its bounds depend on at least one binary variable called the indicator,
 * and indicator = 0 &rArr; \f$x = x^0\f$ for some real constant \f$x^0\f$.
 */
static
SCIP_RETCODE varIsSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable to check */
   SCIP_HASHMAP*         scvars,             /**< semicontinuous variable information */
   SCIP_Bool*            result              /**< buffer to store whether var is semicontinuous */
   )
{
   SCIP_Real lb0;
   SCIP_Real ub0;
   SCIP_Real lb1;
   SCIP_Real ub1;
   SCIP_Real glb;
   SCIP_Real gub;
   SCIP_Bool exists;
   int c;
   int pos;
   SCIP_VAR** vlbvars;
   SCIP_VAR** vubvars;
   SCIP_Real* vlbcoefs;
   SCIP_Real* vubcoefs;
   SCIP_Real* vlbconstants;
   SCIP_Real* vubconstants;
   int nvlbs;
   int nvubs;
   SCVARDATA* scvdata;
   SCIP_VAR* bvar;

   assert(scip != NULL);
   assert(var != NULL);
   assert(scvars != NULL);
   assert(result != NULL);

   scvdata = (SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)var);
   if( scvdata != NULL )
   {
      *result = TRUE;
      return SCIP_OKAY;
   }

   vlbvars = SCIPvarGetVlbVars(var);
   vubvars = SCIPvarGetVubVars(var);
   vlbcoefs = SCIPvarGetVlbCoefs(var);
   vubcoefs = SCIPvarGetVubCoefs(var);
   vlbconstants = SCIPvarGetVlbConstants(var);
   vubconstants = SCIPvarGetVubConstants(var);
   nvlbs = SCIPvarGetNVlbs(var);
   nvubs = SCIPvarGetNVubs(var);
   glb = SCIPvarGetLbGlobal(var);
   gub = SCIPvarGetUbGlobal(var);

   *result = FALSE;

   /* Scan through lower bounds; for each binary vlbvar save the corresponding lb0 and lb1.
    * Then check if there is an upper bound with this vlbvar and save ub0 and ub1.
    * If the found bounds imply that the var value is fixed to some val0 when vlbvar = 0,
    * save vlbvar and val0 to scvdata.
    */
   for( c = 0; c < nvlbs; ++c )
   {
      if( SCIPvarGetType(vlbvars[c]) != SCIP_VARTYPE_BINARY )
         continue;

      SCIPdebugMsg(scip, "var <%s>[%f, %f] lower bound: %f <%s> %+f", SCIPvarGetName(var), glb, gub, vlbcoefs[c], SCIPvarGetName(vlbvars[c]), vlbconstants[c]);

      bvar = vlbvars[c];

      lb0 = MAX(vlbconstants[c], glb);
      lb1 = MAX(vlbconstants[c] + vlbcoefs[c], glb);

      /* look for bvar in vubvars */
      if( vubvars != NULL )
         exists = SCIPsortedvecFindPtr((void**)vubvars, SCIPvarComp, bvar, nvubs, &pos);
      else
         exists = FALSE;
      if( exists )
      { /*lint --e{644}*/
         SCIPdebugMsgPrint(scip, ", upper bound: %f <%s> %+f", vubcoefs[pos], SCIPvarGetName(vubvars[pos]), vubconstants[pos]); /*lint !e613*/

         /* save the upper bounds */
         ub0 = MIN(vubconstants[pos], gub);
         ub1 = MIN(vubconstants[pos] + vubcoefs[pos], gub);
      }
      else
      {
         /* if there is no upper bound with vubvar = bvar, use global var bounds */
         ub0 = gub;
         ub1 = gub;
      }

      /* the 'off' domain of a semicontinuous var should reduce to a single point and be different from the 'on' domain */
      SCIPdebugMsgPrint(scip, " -> <%s> in [%f, %f] (off), [%f, %f] (on)\n", SCIPvarGetName(var), lb0, ub0, lb1, ub1);
      if( SCIPisEQ(scip, lb0, ub0) && (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
      {
         if( scvdata == NULL )
         {
            SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
         }

         SCIP_CALL( addSCVarIndicator(scip, scvdata, bvar, lb0, lb1, ub1) );
      }
   }

   /* look for vubvars that have not been processed yet */
   assert(vubvars != NULL || nvubs == 0);
   for( c = 0; c < nvubs; ++c )
   {
      /* coverity[forward_null] */
      if( SCIPvarGetType(vubvars[c]) != SCIP_VARTYPE_BINARY )  /*lint !e613*/
         continue;

      bvar = vubvars[c];  /*lint !e613*/

      /* skip vars that are in vlbvars */
      if( vlbvars != NULL && SCIPsortedvecFindPtr((void**)vlbvars, SCIPvarComp, bvar, nvlbs, &pos) )
         continue;

      SCIPdebugMsg(scip, "var <%s>[%f, %f] upper bound: %f <%s> %+f",
         SCIPvarGetName(var), glb, gub, vubcoefs[c], SCIPvarGetName(vubvars[c]), vubconstants[c]);  /*lint !e613*/

      lb0 = glb;
      lb1 = glb;
      ub0 = MIN(vubconstants[c], gub);
      ub1 = MIN(vubconstants[c] + vubcoefs[c], gub);

      /* the 'off' domain of a semicontinuous var should reduce to a single point and be different from the 'on' domain */
      SCIPdebugMsgPrint(scip, " -> <%s> in [%f, %f] (off), [%f, %f] (on)\n", SCIPvarGetName(var), lb0, ub0, lb1, ub1);
      if( SCIPisEQ(scip, lb0, ub0) && (!SCIPisEQ(scip, lb0, lb1) || !SCIPisEQ(scip, ub0, ub1)) )
      {
         if( scvdata == NULL )
         {
            SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
         }

         SCIP_CALL( addSCVarIndicator(scip, scvdata, bvar, lb0, lb1, ub1) );
      }
   }

   if( scvdata != NULL )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "var <%s> has global bounds [%f, %f] and the following on/off bounds:\n", SCIPvarGetName(var), glb, gub);
      for( c = 0; c < scvdata->nbnds; ++c )
      {
         SCIPdebugMsg(scip, " c = %d, bvar <%s>: val0 = %f\n", c, SCIPvarGetName(scvdata->bvars[c]), scvdata->vals0[c]);
      }
#endif
      SCIP_CALL( SCIPhashmapInsert(scvars, var, scvdata) );
      *result = TRUE;
   }

   return SCIP_OKAY;
}

/*
 * Semicontinuous expression methods
 */

/* checks if an expression is semicontinuous
 *
 * An expression is semicontinuous if all of its nonlinear variables are semicontinuous
 * and share at least one common indicator variable
 */
static
SCIP_RETCODE exprIsSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool*            res                 /**< buffer to store whether the expression is semicontinuous */
   )
{
   int v;
   SCIP_Bool var_is_sc;
   SCVARDATA* scvdata;
   SCIP_VAR* var;
   int nindicators;
   int nbnds0;
   int c;
   SCIP_VAR** indicators;
   SCIP_Bool* nonlinear;

   *res = FALSE;

   /* constant expression is not semicontinuous; variable expressions are of no interest here */
   if( nlhdlrexprdata->nvars == 0 )
      return SCIP_OKAY;

   indicators = NULL;
   nindicators = 0;
   nbnds0 = 0;

   if( SCIPisExprSum(scip, expr) )
   {
      SCIP_EXPRITER* it;
      SCIP_EXPR* child;
      SCIP_EXPR* curexpr;
      int pos;
      SCIP_Bool issc;

      /* sums are treated separately because if there are variables that are non-semicontinuous but
       * appear only linearly, we still want to apply perspective to expr
       */

      SCIP_CALL( SCIPallocClearBufferArray(scip, &nonlinear, nlhdlrexprdata->nvars) );
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );

      for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      {
         child = SCIPexprGetChildren(expr)[c];

         if( SCIPisExprVar(scip, child) )
         {
            var = SCIPgetVarExprVar(child);

            /* save information on semicontinuity of child */
            SCIP_CALL( varIsSemicontinuous(scip, var, nlhdlrdata->scvars, &var_is_sc) );

            /* since child is a variable, go on regardless of the value of var_is_sc */
            continue;
         }

         issc = TRUE;

         SCIP_CALL( SCIPexpriterInit(it, child, SCIP_EXPRITER_DFS, FALSE) );
         curexpr = SCIPexpriterGetCurrent(it);

         /* all nonlinear terms of a sum should be semicontinuous in original variables */
         while( !SCIPexpriterIsEnd(it) )
         {
            assert(curexpr != NULL);

            if( SCIPisExprVar(scip, curexpr) )
            {
               var = SCIPgetVarExprVar(curexpr);

               if( !SCIPvarIsRelaxationOnly(var) )
               {
                  SCIP_CALL( varIsSemicontinuous(scip, var, nlhdlrdata->scvars, &var_is_sc) );

                  /* mark the variable as nonlinear */
                  (void) SCIPsortedvecFindPtr((void**) nlhdlrexprdata->vars, SCIPvarComp, (void*) var, nlhdlrexprdata->nvars,
                        &pos);
                  assert(0 <= pos && pos < nlhdlrexprdata->nvars);
                  nonlinear[pos] = TRUE;

                  if( !var_is_sc )
                  {
                     /* non-semicontinuous child which is (due to a previous check) not a var ->
                      * expr is non-semicontinuous
                      */
                     issc = FALSE;
                     break;
                  }
               }
            }
            curexpr = SCIPexpriterGetNext(it);
         }

         if( !issc )
         {
            SCIPfreeExpriter(&it);
            goto TERMINATE;
         }
      }
      SCIPfreeExpriter(&it);
   }
   else
   {
      /* non-sum expression */
      nonlinear = NULL;

      /* all variables of a non-sum on/off expression should be semicontinuous */
      for( v = 0; v < nlhdlrexprdata->nvars; ++v )
      {
         SCIP_CALL( varIsSemicontinuous(scip, nlhdlrexprdata->vars[v], nlhdlrdata->scvars, &var_is_sc) );
         if( !var_is_sc )
            return SCIP_OKAY;
      }
   }

   /* look for common binary variables for all variables of the expression */

   SCIPdebugMsg(scip, "Array intersection for var <%s>\n", SCIPvarGetName(nlhdlrexprdata->vars[0]));
   for( v = 0; v < nlhdlrexprdata->nvars; ++v )
   {
      SCIPdebugMsg(scip, "%s; \n", SCIPvarGetName(nlhdlrexprdata->vars[v]));

      if( nonlinear != NULL && !nonlinear[v] )
         continue;

      scvdata = (SCVARDATA*)SCIPhashmapGetImage(nlhdlrdata->scvars, (void*) nlhdlrexprdata->vars[v]);

      /* we should have exited earlier if there is a nonlinear non-semicontinuous variable */
      assert(scvdata != NULL);

      if( indicators == NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &indicators, scvdata->bvars, scvdata->nbnds) );
         nbnds0 = scvdata->nbnds;
         nindicators = nbnds0;
      }
      else
      {
         SCIPcomputeArraysIntersectionPtr((void**)indicators, nindicators, (void**)scvdata->bvars, scvdata->nbnds,
               SCIPvarComp, (void**)indicators, &nindicators);
      }

      /* if we have found out that the intersection is empty, expr is not semicontinuous */
      if( indicators != NULL && nindicators == 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &indicators, nbnds0);
         goto TERMINATE;
      }
   }

   /* this can happen if all children are linear vars and none are semicontinuous */
   if( indicators == NULL )
   {
      goto TERMINATE;
   }
   assert(nindicators > 0 && nindicators <= nbnds0);

   if( nindicators < nbnds0 )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &indicators, nbnds0, nindicators) );
   }

   for( v = 0; v < nindicators; ++v )
   {
      SCIP_CALL( SCIPcaptureVar(scip, indicators[v]) );
   }
   nlhdlrexprdata->indicators = indicators;
   nlhdlrexprdata->nindicators = nindicators;
   *res = TRUE;

 TERMINATE:
   SCIPfreeBufferArrayNull(scip, &nonlinear);

   return SCIP_OKAY;
}

/** computes the 'off' value of the expression and the 'off' values of
  * semicontinuous auxiliary variables for each indicator variable
  */
static
SCIP_RETCODE computeOffValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPRITER* it;
   SCIP_SOL* sol;
   int i;
   int v;
   int norigvars;
   SCIP_Real* origvals0;
   SCIP_VAR** origvars;
   SCVARDATA* scvdata;
   SCIP_VAR* auxvar;
   SCIP_EXPR* curexpr;
   SCIP_HASHMAP* auxvarmap;
   SCIP_Bool hasnonsc;
   int pos;

   assert(expr != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nlhdlrexprdata->exprvals0), nlhdlrexprdata->nindicators) );
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &origvals0, nlhdlrexprdata->nvars) );
   SCIP_CALL( SCIPhashmapCreate(&auxvarmap, SCIPblkmem(scip), 10) );
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &origvars, nlhdlrexprdata->vars, nlhdlrexprdata->nvars) );
   norigvars = nlhdlrexprdata->nvars;

   for( i = nlhdlrexprdata->nindicators - 1; i >= 0; --i )
   {
      hasnonsc = FALSE;

      /* set sol to the off value of all expr vars for this indicator */
      for( v = 0; v < norigvars; ++v )
      {
         /* set vals0[v] = 0 if var is non-sc with respect to indicators[i] - then it will not
          * contribute to exprvals0[i] since any non-sc var must be linear
          */
         scvdata = getSCVarDataInd(nlhdlrdata->scvars, origvars[v], nlhdlrexprdata->indicators[i], &pos);
         if( scvdata == NULL )
         {
            origvals0[v] = 0.0;
            hasnonsc = TRUE;
         }
         else
         {
            origvals0[v] = scvdata->vals0[pos];
         }
      }
      SCIP_CALL( SCIPsetSolVals(scip, sol, norigvars, origvars, origvals0) );
      SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0L) );

      if( SCIPexprGetEvalValue(expr) == SCIP_INVALID ) /*lint !e777*/
      {
         SCIPdebugMsg(scip, "expression evaluation failed for %p, removing indicator %s\n",
                             (void*)expr, SCIPvarGetName(nlhdlrexprdata->indicators[i]));
         /* TODO should we fix the indicator variable to 1? */
         /* since the loop is backwards, this only modifies the already processed part of nlhdlrexprdata->indicators */
         SCIP_CALL( removeIndicator(scip, nlhdlrexprdata, i) );
         continue;
      }

      nlhdlrexprdata->exprvals0[i] = SCIPexprGetEvalValue(expr);

      /* iterate through the expression and create scvdata for aux vars */
      SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
      curexpr = SCIPexpriterGetCurrent(it);

      while( !SCIPexpriterIsEnd(it) )
      {
         auxvar = SCIPgetExprAuxVarNonlinear(curexpr);

         if( auxvar != NULL )
         {
            SCIP_Bool issc = TRUE;
#ifndef NDEBUG
            SCIP_EXPR** childvarexprs;
            int nchildvarexprs;
            SCIP_VAR* var;
#endif

            if( hasnonsc )
            {
               /* expr is a sum with non-semicontinuous linear terms. Therefore, curexpr might be
                * non-semicontinuous. In that case the auxvar is also non-semicontinuous, so
                * we will skip on/off bounds computation.
                */
               if( SCIPisExprVar(scip, curexpr) )
               {
                  /* easy case: curexpr is a variable, can check semicontinuity immediately */
                  scvdata = getSCVarDataInd(nlhdlrdata->scvars, SCIPgetVarExprVar(curexpr),
                        nlhdlrexprdata->indicators[i], &pos);
                  issc = scvdata != NULL;
               }
               else if( SCIPisExprSum(scip, curexpr) && curexpr == expr )
               {
                  /* if expr itself is a sum, this is an exception since a sum with nonlinear terms is
                   * allowed to have both semicontinuous and non-semicontinuous variables; we skip it here
                   * and then analyse it term by term
                   */
                  issc = FALSE;
               }

#ifndef NDEBUG
               if( !SCIPisExprVar(scip, curexpr) && (!SCIPisExprSum(scip, curexpr) || curexpr != expr) )
               {
                  /* curexpr is a non-variable expression and does not fit the sum special case,
                   * so it belongs to the non-linear part of expr.
                   * Since the non-linear part of expr must be semicontinuous with respect to
                   * nlhdlrexprdata->indicators[i], curexpr must be semicontinuous
                   */

                  SCIP_CALL( SCIPallocBufferArray(scip, &childvarexprs, norigvars) );
                  SCIP_CALL( SCIPgetExprVarExprs(scip, curexpr, childvarexprs, &nchildvarexprs) );

                  /* all nonlinear variables of a sum on/off term should be semicontinuous */
                  for( v = 0; v < nchildvarexprs; ++v )
                  {
                     var = SCIPgetVarExprVar(childvarexprs[v]);
                     scvdata = getSCVarDataInd(nlhdlrdata->scvars, var, nlhdlrexprdata->indicators[i], &pos);
                     assert(scvdata != NULL);

                     SCIP_CALL( SCIPreleaseExpr(scip, &childvarexprs[v]) );
                  }

                  SCIPfreeBufferArray(scip, &childvarexprs);
               }
#endif
            }

            if( issc )
            {
               /* we know that all vars are semicontinuous with respect to exprdata->indicators; it remains to:
                * - get or create the scvardata structure for auxvar
                * - if had to create scvardata, add it to scvars hashmap
                * - add the indicator and the off value (= curexpr's off value) to scvardata
                */
               scvdata = (SCVARDATA*) SCIPhashmapGetImage(nlhdlrdata->scvars, (void*)auxvar);
               if( scvdata == NULL )
               {
                  SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->bvars,  nlhdlrexprdata->nindicators) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->vals0, nlhdlrexprdata->nindicators) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->lbs1, nlhdlrexprdata->nindicators) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->ubs1, nlhdlrexprdata->nindicators) );
                  scvdata->bndssize = nlhdlrexprdata->nindicators;
                  SCIP_CALL( SCIPhashmapInsert(nlhdlrdata->scvars, auxvar, scvdata) );
               }

               SCIP_CALL( addSCVarIndicator(scip, scvdata, nlhdlrexprdata->indicators[i],
                     SCIPexprGetEvalValue(curexpr), SCIPvarGetLbGlobal(auxvar), SCIPvarGetUbGlobal(auxvar)) );
            }

            SCIP_CALL( addAuxVar(scip, nlhdlrexprdata, auxvarmap, auxvar) );
         }

         curexpr = SCIPexpriterGetNext(it);
      }
   }

   SCIPfreeExpriter(&it);
   SCIPhashmapFree(&auxvarmap);
   SCIPfreeBufferArray(scip, &origvars);
   SCIPfreeBufferArray(scip, &origvals0);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   return SCIP_OKAY;
}

/*
 * Probing and bound tightening methods
 */

/** go into probing and set some variable bounds */
static
SCIP_RETCODE startProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_VAR*             indicator,          /**< indicator variable */
   SCIP_VAR**            probingvars,        /**< array of vars whose bounds we will change in probing */
   SCIP_INTERVAL*        probingdoms,        /**< array of intervals to which bounds of probingvars will be changed in probing */
   int                   nprobingvars,       /**< number of probing vars */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_SOL**            solcopy,            /**< buffer for a copy of sol before going into probing; if *solcopy == sol, then copy is created */
   SCIP_Bool*            cutoff_probing      /**< pointer to store whether indicator == 1 is infeasible */
   )
{
   int v;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Bool propagate;

   propagate = SCIPgetDepth(scip) == 0;

   /* if a copy of sol has not been created yet, then create one now and copy the relevant var values from sol,
    * because sol can change after SCIPstartProbing, e.g., when linked to the LP solution
    */
   if( *solcopy == sol )
   {
      SCIP_CALL( SCIPcreateSol(scip, solcopy, NULL) );
      for( v = 0; v < nlhdlrexprdata->nvars; ++v )
      {
         SCIP_CALL( SCIPsetSolVal(scip, *solcopy, nlhdlrexprdata->vars[v], SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[v])) );
      }
      for( v = 0; v < nlhdlrexprdata->nindicators; ++v )
      {
         SCIP_CALL( SCIPsetSolVal(scip, *solcopy, nlhdlrexprdata->indicators[v], SCIPgetSolVal(scip, sol, nlhdlrexprdata->indicators[v])) );
      }
   }

   /* go into probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* create a probing node */
   SCIP_CALL( SCIPnewProbingNode(scip) );

   /* set indicator to 1 */
   SCIP_CALL( SCIPchgVarLbProbing(scip, indicator, 1.0) );

   /* apply stored bounds */
   for( v = 0; v < nprobingvars; ++v )
   {
      newlb = SCIPintervalGetInf(probingdoms[v]);
      newub = SCIPintervalGetSup(probingdoms[v]);

      if( SCIPisGT(scip, newlb, SCIPvarGetLbLocal(probingvars[v])) || (newlb >= 0.0 && SCIPvarGetLbLocal(probingvars[v]) < 0.0) )
      {
         SCIP_CALL( SCIPchgVarLbProbing(scip, probingvars[v], newlb) );
      }
      if( SCIPisLT(scip, newub, SCIPvarGetUbLocal(probingvars[v])) || (newub <= 0.0 && SCIPvarGetUbLocal(probingvars[v]) > 0.0) )
      {
         SCIP_CALL( SCIPchgVarUbProbing(scip, probingvars[v], newub) );
      }
   }

   if( propagate )
   {
      SCIP_Longint ndomreds;

      SCIP_CALL( SCIPpropagateProbing(scip, nlhdlrdata->maxproprounds, cutoff_probing, &ndomreds) );
   }

   return SCIP_OKAY;
}

/** analyse on/off bounds on a variable
 *
 * analyses for
 * 1. tightening bounds in probing for indicator = 1,
 * 2. fixing indicator / detecting cutoff if one or both states are infeasible,
 * 3. tightening local bounds if indicator is fixed.
 *
 * `probinglb` and `probingub` are only set if `doprobing` is TRUE.
 * They are either set to bounds that should be used in probing or to `SCIP_INVALID` if bounds on
 * `var` shouldn't be changed in probing.
 */
static
SCIP_RETCODE analyseVarOnoffBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             indicator,          /**< indicator variable */
   SCIP_Bool             indvalue,           /**< indicator value for which the bounds are applied */
   SCIP_Bool*            infeas,             /**< pointer to store whether infeasibility has been detected */
   SCIP_Real*            probinglb,          /**< pointer to store the lower bound to be applied in probing */
   SCIP_Real*            probingub,          /**< pointer to store the upper bound to be applied in probing */
   SCIP_Bool             doprobing,          /**< whether we currently consider to go into probing */
   SCIP_Bool*            reduceddom          /**< pointer to store whether any variables were fixed */
   )
{
   SCVARDATA* scvdata;
   int pos;
   SCIP_Real sclb;
   SCIP_Real scub;
   SCIP_Real loclb;
   SCIP_Real locub;
   SCIP_Bool bndchgsuccess;

   assert(var != NULL);
   assert(indicator != NULL);
   assert(infeas != NULL);
   assert(reduceddom != NULL);

   /* shouldn't be called if indicator is fixed to !indvalue */
   assert((indvalue && SCIPvarGetUbLocal(indicator) > 0.5) || (!indvalue && SCIPvarGetLbLocal(indicator) < 0.5));

   *infeas = FALSE;
   *reduceddom = FALSE;
   scvdata = getSCVarDataInd(nlhdlrdata->scvars, var, indicator, &pos);
   if( doprobing )
   {
      assert(probinglb != NULL);
      assert(probingub != NULL);

      *probinglb = SCIP_INVALID;
      *probingub = SCIP_INVALID;
   }

   /* nothing to do for non-semicontinuous variables */
   if( scvdata == NULL )
   {
      return SCIP_OKAY;
   }

   sclb = indvalue ? scvdata->lbs1[pos] : scvdata->vals0[pos];
   scub = indvalue ? scvdata->ubs1[pos] : scvdata->vals0[pos];
   loclb = SCIPvarGetLbLocal(var);
   locub = SCIPvarGetUbLocal(var);

   /* nothing to do for fixed variables */
   if( SCIPisEQ(scip, loclb, locub) )
      return SCIP_OKAY;

   /* use a non-redundant lower bound */
   if( SCIPisGT(scip, sclb, SCIPvarGetLbLocal(var)) || (sclb >= 0.0 && loclb < 0.0) )
   {
      /* first check for infeasibility */
      if( SCIPisFeasGT(scip, sclb, SCIPvarGetUbLocal(var)) )
      {
         SCIP_CALL( SCIPfixVar(scip, indicator, indvalue ? 0.0 : 1.0, infeas, &bndchgsuccess) );
         *reduceddom += bndchgsuccess;
         if( *infeas )
         {
            return SCIP_OKAY;
         }
      }
      else if( nlhdlrdata->tightenbounds &&
              (SCIPvarGetUbLocal(indicator) <= 0.5 || SCIPvarGetLbLocal(indicator) >= 0.5) )
      {
         /* indicator is fixed; due to a previous check, here it can only be fixed to indvalue;
          * therefore, sclb is valid for the current node
          */

         if( indvalue == 0 )
         {
            assert(sclb == scub); /*lint !e777*/
            SCIP_CALL( SCIPfixVar(scip, var, sclb, infeas, &bndchgsuccess) );
         }
         else
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, sclb, FALSE, infeas, &bndchgsuccess) );
         }
         *reduceddom += bndchgsuccess;
         if( *infeas )
         {
            return SCIP_OKAY;
         }
      }
   }

   /* use a non-redundant upper bound */
   if( SCIPisLT(scip, scub, SCIPvarGetUbLocal(var)) || (scub <= 0.0 && locub > 0.0) )
   {
      /* first check for infeasibility */
      if( SCIPisFeasLT(scip, scub, SCIPvarGetLbLocal(var)) )
      {
         SCIP_CALL( SCIPfixVar(scip, indicator, indvalue ? 0.0 : 1.0, infeas, &bndchgsuccess) );
         *reduceddom += bndchgsuccess;
         if( *infeas )
         {
            return SCIP_OKAY;
         }
      }
      else if( nlhdlrdata->tightenbounds &&
              (SCIPvarGetUbLocal(indicator) <= 0.5 || SCIPvarGetLbLocal(indicator) >= 0.5) )
      {
         /* indicator is fixed; due to a previous check, here it can only be fixed to indvalue;
          * therefore, scub is valid for the current node
          */

         if( indvalue == 0 )
         {
            assert(sclb == scub); /*lint !e777*/
            SCIP_CALL( SCIPfixVar(scip, var, scub, infeas, &bndchgsuccess) );
         }
         else
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, scub, FALSE, infeas, &bndchgsuccess) );
         }
         *reduceddom += bndchgsuccess;
         if( *infeas )
         {
            return SCIP_OKAY;
         }
      }
   }

   /* If a bound change has been found and indvalue == TRUE, try to use the new bounds.
    * This is only done for indvalue == TRUE since this is where enfo asks other nlhdlrs to estimate,
    * and at indicator == FALSE we already only have a single point
    */
   if( doprobing && indvalue && (((scub - sclb) / (locub - loclb)) <= 1.0 - nlhdlrdata->mindomreduction ||
       (sclb >= 0.0 && loclb < 0.0) || (scub <= 0.0 && locub > 0.0)) )
   {
      *probinglb = sclb;
      *probingub = scub;
   }

   SCIPdebugMsg(scip, "%s in [%g, %g] instead of [%g, %g] (vals0 = %g)\n", SCIPvarGetName(var), sclb, scub,
                SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), scvdata->vals0[pos]);

   return SCIP_OKAY;
}

/** looks for bound tightenings to be applied either in the current node or in probing
 *
 * Loops through both possible values of indicator and calls analyseVarOnoffBounds().
 * Might update the `*doprobing` flag by setting it to `FALSE` if:
 * - indicator is fixed or
 * - analyseVarOnoffBounds() hasn't found a sufficient improvement at indicator==1.
 *
 * If `*doprobing==TRUE`, stores bounds suggested by analyseVarOnoffBounds() in order to apply them in probing together
 * with the fixing `indicator=1`.
 */
static
SCIP_RETCODE analyseOnoffBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_VAR*             indicator,          /**< indicator variable */
   SCIP_VAR***           probingvars,        /**< array to store variables whose bounds will be changed in probing */
   SCIP_INTERVAL**       probingdoms,        /**< array to store bounds to be applied in probing */
   int*                  nprobingvars,       /**< pointer to store number of vars whose bounds will be changed in probing */
   SCIP_Bool*            doprobing,          /**< pointer to the flag telling whether we want to do probing */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   int v;
   SCIP_VAR* var;
   SCIP_Bool infeas;
   int b;
   SCIP_Real probinglb = SCIP_INVALID;
   SCIP_Real probingub = SCIP_INVALID;
   SCIP_Bool changed;
   SCIP_Bool reduceddom;

   assert(indicator != NULL);
   assert(nprobingvars != NULL);
   assert(doprobing != NULL);
   assert(result != NULL);

   changed = FALSE;

   /* no probing if indicator already fixed */
   if( SCIPvarGetUbLocal(indicator) <= 0.5 || SCIPvarGetLbLocal(indicator) >= 0.5 )
   {
      *doprobing = FALSE;
   }

   /* consider each possible value of indicator */
   for( b = 0; b < 2; ++b )
   {
      for( v = 0; v < nlhdlrexprdata->nvars; ++v )
      {
         /* nothing left to do if indicator is already fixed to !indvalue
          * (checked in the inner loop since analyseVarOnoff bounds might fix the indicator)
          */
         if( (b == 1 && SCIPvarGetUbLocal(indicator) <= 0.5) || (b == 0 && SCIPvarGetLbLocal(indicator) >= 0.5) )
         {
            *doprobing = FALSE;
            break;
         }

         var = nlhdlrexprdata->vars[v];

         SCIP_CALL( analyseVarOnoffBounds(scip, nlhdlrdata, var, indicator, b == 1, &infeas, &probinglb,
               &probingub, *doprobing, &reduceddom) );

         if( infeas )
         {
            *result = SCIP_CUTOFF;
            *doprobing = FALSE;
            return SCIP_OKAY;
         }
         else if( reduceddom )
         {
            *result = SCIP_REDUCEDDOM;
         }

         if( !(*doprobing) )
            continue;

         /* if bounds to be applied in probing have been found, store them */
         if( probinglb != SCIP_INVALID ) /*lint !e777*/
         {
            assert(probingub != SCIP_INVALID); /*lint !e777*/

            SCIP_CALL( SCIPreallocBufferArray(scip, probingvars, *nprobingvars + 1) );
            SCIP_CALL( SCIPreallocBufferArray(scip, probingdoms, *nprobingvars + 1) );
            (*probingvars)[*nprobingvars] = var;
            (*probingdoms)[*nprobingvars].inf = probinglb;
            (*probingdoms)[*nprobingvars].sup = probingub;
            ++*nprobingvars;

            changed = TRUE;
         }
      }
   }

   if( !changed )
   {
      *doprobing = FALSE;
   }

   return SCIP_OKAY;
}

/** saves local bounds on all expression variables, including auxiliary variables, obtained from propagating
 * indicator == 1 to the corresponding SCVARDATA (should only be used in the root node)
 */
static
SCIP_RETCODE tightenOnBounds(
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nlhdlr expression data */
   SCIP_HASHMAP*         scvars,             /**< hashmap with semicontinuous variables */
   SCIP_VAR*             indicator           /**< indicator variable */
   )
{
   int v;
   SCIP_VAR* var;
   SCVARDATA* scvdata;
   int pos;
   SCIP_Real lb;
   SCIP_Real ub;

   for( v = 0; v < nlhdlrexprdata->nvars; ++v )
   {
      var = nlhdlrexprdata->vars[v];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      scvdata = getSCVarDataInd(scvars, var, indicator, &pos);

      if( scvdata != NULL )
      {
         scvdata->lbs1[pos] = MAX(scvdata->lbs1[pos], lb);
         scvdata->ubs1[pos] = MIN(scvdata->ubs1[pos], ub);
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrPerspective)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrPerspective(targetscip) );

   return SCIP_OKAY;
}


/** callback to free data of handler */
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataPerspective)
{ /*lint --e{715}*/
   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}


/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataPerspective)
{  /*lint --e{715}*/
   SCIP_CALL( freeNlhdlrExprData(scip, *nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** callback to be called in initialization */
#if 0
static
SCIP_DECL_NLHDLRINIT(nlhdlrInitPerspective)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}
#endif

/** callback to be called in deinitialization */
static
SCIP_DECL_NLHDLREXIT(nlhdlrExitPerspective)
{  /*lint --e{715}*/
   SCIP_HASHMAPENTRY* entry;
   SCVARDATA* data;
   int c;
   SCIP_NLHDLRDATA* nlhdlrdata;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->scvars != NULL )
   {
      for( c = 0; c < SCIPhashmapGetNEntries(nlhdlrdata->scvars); ++c )
      {
         entry = SCIPhashmapGetEntry(nlhdlrdata->scvars, c);
         if( entry != NULL )
         {
            data = (SCVARDATA*) SCIPhashmapEntryGetImage(entry);
            SCIPfreeBlockMemoryArray(scip, &data->ubs1, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->lbs1, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->vals0, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->bvars, data->bndssize);
            SCIPfreeBlockMemory(scip, &data);
         }
      }
      SCIPhashmapFree(&nlhdlrdata->scvars);
      assert(nlhdlrdata->scvars == NULL);
   }

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree
 *
 *  We are looking for expressions g(x), where x is a vector of semicontinuous variables that all share at least one
 *  indicator variable.
 */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectPerspective)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_EXPR** varexprs;
   SCIP_Bool success = FALSE;
   int i;
   SCIP_Bool hassepabelow = FALSE;
   SCIP_Bool hassepaabove = FALSE;
   SCIP_Bool hasnondefault = FALSE;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(participating != NULL);
   assert(enforcing != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrdata != NULL);

   /* do not run if we will have no auxvar to add a cut for */
   if( SCIPgetExprNAuxvarUsesNonlinear(expr) == 0 )
      return SCIP_OKAY;

   if( SCIPgetNBinVars(scip) == 0 )
   {
      SCIPdebugMsg(scip, "problem has no binary variables, not running perspective detection\n");
      return SCIP_OKAY;
   }

   for( i = 0; i < SCIPgetExprNEnfosNonlinear(expr); ++i )
   {
      SCIP_NLHDLR* nlhdlr2;
      SCIP_NLHDLR_METHOD nlhdlr2participates;
      SCIP_Bool sepabelowusesactivity;
      SCIP_Bool sepaaboveusesactivity;
      SCIPgetExprEnfoDataNonlinear(expr, i, &nlhdlr2, NULL, &nlhdlr2participates, &sepabelowusesactivity, &sepaaboveusesactivity, NULL);

      if( (nlhdlr2participates & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 )
         continue;

      if( !SCIPnlhdlrHasEstimate(nlhdlr2) )
         continue;

      if( strcmp(SCIPnlhdlrGetName(nlhdlr2), "default") != 0 )
         hasnondefault = TRUE;

      /* If we are supposed to run only on convex expressions, than check whether there is a nlhdlr
       * that participates in separation without using activity for it. Otherwise, check for
       * participation regardless of activity usage.
       */
      if( (nlhdlr2participates & SCIP_NLHDLR_METHOD_SEPABELOW) && (!nlhdlrdata->convexonly || !sepabelowusesactivity) )
         hassepabelow = TRUE;

      if( (nlhdlr2participates & SCIP_NLHDLR_METHOD_SEPAABOVE) && (!nlhdlrdata->convexonly || !sepaaboveusesactivity) )
         hassepaabove = TRUE;
   }

   /* If a sum expression is handled only by default nlhdlr, then all the children will have auxiliary vars.
    * Since the sum will then be linear in auxiliary variables, perspective can't improve anything for it
    */
   if( SCIPisExprSum(scip, expr) && !hasnondefault )
   {
      SCIPdebugMsg(scip, "sum expr only has default exprhdlr, not running perspective detection\n");
      return SCIP_OKAY;
   }

   /* If no other nlhdlr separates, neither does perspective (if convexonly, only separation
    * without using activity counts)
    */
   if( !hassepabelow && !hassepaabove )
   {
      SCIPdebugMsg(scip, "no nlhdlr separates without using activity, not running perspective detection\n");
      return SCIP_OKAY;
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Called perspective detect, expr = %p: ", (void*)expr);
   SCIPprintExpr(scip, expr, NULL);
   SCIPdebugMsgPrint(scip, "\n");
#endif

   /* allocate memory */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   if( nlhdlrdata->scvars == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&(nlhdlrdata->scvars), SCIPblkmem(scip), SCIPgetNVars(scip)) );
   }

   /* save varexprs to nlhdlrexprdata */
   SCIP_CALL( SCIPgetExprNVars(scip, expr, &(*nlhdlrexprdata)->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, (*nlhdlrexprdata)->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, (*nlhdlrexprdata)->nvars) );
   (*nlhdlrexprdata)->varssize = (*nlhdlrexprdata)->nvars;
   SCIP_CALL( SCIPgetExprVarExprs(scip, expr, varexprs, &(*nlhdlrexprdata)->nvars) );
   for( i = 0; i < (*nlhdlrexprdata)->nvars; ++i )
   {
      (*nlhdlrexprdata)->vars[i] = SCIPgetVarExprVar(varexprs[i]);
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[i]) );
      SCIP_CALL( SCIPcaptureVar(scip, (*nlhdlrexprdata)->vars[i]) );
   }
   SCIPsortPtr((void**) (*nlhdlrexprdata)->vars, SCIPvarComp, (*nlhdlrexprdata)->nvars);
   SCIPfreeBufferArray(scip, &varexprs);

   /* check if expr is semicontinuous and save indicator variables */
   SCIP_CALL( exprIsSemicontinuous(scip, nlhdlrdata, *nlhdlrexprdata, expr, &success) );

   if( success )
   {
      assert(*nlhdlrexprdata != NULL);
      assert((*nlhdlrexprdata)->nindicators > 0);

      if( hassepaabove )
         *participating |= SCIP_NLHDLR_METHOD_SEPAABOVE;
      if( hassepabelow )
         *participating |= SCIP_NLHDLR_METHOD_SEPABELOW;

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "detected an on/off expr: ");
      SCIPprintExpr(scip, expr, NULL);
      SCIPinfoMessage(scip, NULL, "\n");
#endif
   }
   else
   {
      assert(*nlhdlrexprdata != NULL);
      SCIP_CALL( nlhdlrFreeExprDataPerspective(scip, nlhdlr, expr, nlhdlrexprdata) );
   }

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxPerspective)
{ /*lint --e{715}*/
   int e;
   SCIP_Real maxdiff;
   SCIP_Real auxvarvalue;
   SCIP_Real enfoauxval;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(auxvalue != NULL);

   auxvarvalue = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr));
   maxdiff = 0.0;
   *auxvalue = auxvarvalue;

   /* use the auxvalue from one of the other nlhdlrs that estimates for this expr: take the one that is farthest
    * from the current value of auxvar
    */
   for( e = 0; e < SCIPgetExprNEnfosNonlinear(expr); ++e )
   {
      SCIP_NLHDLR* nlhdlr2;
      SCIP_NLHDLREXPRDATA* nlhdlr2exprdata;
      SCIP_NLHDLR_METHOD nlhdlr2participation;

      SCIPgetExprEnfoDataNonlinear(expr, e, &nlhdlr2, &nlhdlr2exprdata, &nlhdlr2participation, NULL, NULL, NULL);

      /* skip nlhdlr that do not participate or do not provide estimate */
      if( (nlhdlr2participation & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 || !SCIPnlhdlrHasEstimate(nlhdlr2) )
         continue;

      SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr2, expr, nlhdlr2exprdata, &enfoauxval, sol) );

      SCIPsetExprEnfoAuxValueNonlinear(expr, e, enfoauxval);

      if( REALABS(enfoauxval - auxvarvalue) > maxdiff && enfoauxval != SCIP_INVALID ) /*lint !e777*/
      {
         maxdiff = REALABS(enfoauxval - auxvarvalue);
         *auxvalue = enfoauxval;
      }
   }

   return SCIP_OKAY;
}

/** separation initialization method of a nonlinear handler */
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaPerspective)
{ /*lint --e{715}*/
   int sindicators;

   sindicators = nlhdlrexprdata->nindicators;

   /* compute 'off' values of expr and subexprs (and thus auxvars too) */
   SCIP_CALL( computeOffValues(scip, SCIPnlhdlrGetData(nlhdlr), nlhdlrexprdata, expr) );

   /* some indicator variables might have been removed if evaluation failed, check how many remain */
   if( nlhdlrexprdata->nindicators == 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->indicators, sindicators);
      SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->exprvals0, sindicators);
   }
   else if( nlhdlrexprdata->nindicators < sindicators )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->indicators, sindicators, nlhdlrexprdata->nindicators) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->exprvals0, sindicators, nlhdlrexprdata->nindicators) );
   }

   return SCIP_OKAY;
}


#if 0
/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
static
SCIP_DECL_NLHDLREXITSEPA(nlhdlrExitSepaPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#endif

/** nonlinear handler enforcement callback
 *
 * "Perspectivies" cuts produced by other nonlinear handlers.
 *
 * Suppose that we want to separate \f$x\f$ from the set \f$\{ x : g(x) \leq 0\}\f$.
 * If \f$g(x) = g^0\f$ if indicator \f$z = 0\f$, and a cut is given by \f$\sum_i a_ix_i + c \leq \text{aux}\f$, where \f$x_i = x_i^0\f$ if \f$z = 0\f$ for all \f$i\f$,
 * then the "perspectivied" cut is \f[\sum_i a_ix_i + c + (1 - z)\,(g^0 - c - \sum_i a_ix_i^0) \leq \text{aux}.\f]
 * This ensures that at \f$z = 1\f$, the new cut is equivalent to the given cut, and at \f$z = 0\f$ it reduces to \f$g^0 \leq \text{aux}\f$.
 */
static
SCIP_DECL_NLHDLRENFO(nlhdlrEnfoPerspective)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* auxvar;
   int i;
   int j;
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_Real cst0;
   SCIP_VAR* indicator;
   SCIP_PTRARRAY* rowpreps2;
   SCIP_PTRARRAY* rowpreps;
   int nrowpreps;
   SCIP_SOL* solcopy;
   SCIP_Bool doprobing;
   SCIP_BOOLARRAY* addedbranchscores2;
   SCIP_Bool stop;
   int nenfos;
   int* enfoposs;
   SCIP_SOL* soladj;
   int pos;
   SCVARDATA* scvdata;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "enforcement method of perspective nonlinear handler called for expr %p: ", (void*)expr);
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, " at\n");
   for( i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      SCIPinfoMessage(scip, NULL, "%s = %g\n", SCIPvarGetName(nlhdlrexprdata->vars[i]),
              SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[i]));
   }
   SCIPinfoMessage(scip, NULL, "%s = %g", SCIPvarGetName(SCIPgetExprAuxVarNonlinear(expr)),
           SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr)));
#endif

   assert(scip != NULL);
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrdata != NULL);

   if( nlhdlrexprdata->nindicators == 0 )
   {
      /* we might have removed all indicators in initsepa */
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   auxvar = SCIPgetExprAuxVarNonlinear(expr);
   assert(auxvar != NULL);

   /* detect should have picked only those expressions for which at least one other nlhdlr can enforce */
   assert(SCIPgetExprNEnfosNonlinear(expr) > 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &enfoposs, SCIPgetExprNEnfosNonlinear(expr) - 1) );

   doprobing = FALSE;
   nenfos = 0;
   soladj = NULL;

   /* find suitable nlhdlrs and check if there is enough violation to do probing */
   for( j = 0; j < SCIPgetExprNEnfosNonlinear(expr); ++j )
   {
      SCIP_NLHDLR* nlhdlr2;
      SCIP_NLHDLREXPRDATA* nlhdlr2exprdata;
      SCIP_NLHDLR_METHOD nlhdlr2participate;
      SCIP_Real nlhdlr2auxvalue;
      SCIP_Real violation;
      SCIP_Bool violbelow;
      SCIP_Bool violabove;
      SCIP_Bool sepausesactivity = FALSE;

      SCIPgetExprEnfoDataNonlinear(expr, j, &nlhdlr2, &nlhdlr2exprdata, &nlhdlr2participate, !overestimate ? &sepausesactivity : NULL, overestimate ? &sepausesactivity: NULL, &nlhdlr2auxvalue);  /*lint !e826*/

      if( nlhdlr2 == nlhdlr )
         continue;

      /* if nlhdlr2 cannot estimate, then cannot use it */
      if( !SCIPnlhdlrHasEstimate(nlhdlr2) )
         continue;

      /* if nlhdlr2 does not participate in the separation on the desired side (overestimate), then skip it */
      if( (nlhdlr2participate & (overestimate ? SCIP_NLHDLR_METHOD_SEPAABOVE : SCIP_NLHDLR_METHOD_SEPABELOW)) == 0 )
         continue;

      /* if only working on convex-looking expressions, then skip nlhdlr if it uses activity for estimates */
      if( nlhdlrdata->convexonly && sepausesactivity )
         continue;

      /* evalaux should have called evalaux of nlhdlr2 by now
       * check whether handling the violation for nlhdlr2 requires under- or overestimation and this fits to
       * overestimate flag
       */
      SCIP_CALL( SCIPgetExprAbsAuxViolationNonlinear(scip, expr, nlhdlr2auxvalue, sol, &violation, &violbelow,
            &violabove) );
      assert(violation >= 0.0);

      if( (overestimate && !violabove) || (!overestimate && !violbelow) )
         continue;

      /* if violation is small, cuts would likely be weak - skip perspectification */
      if( !allowweakcuts && violation < SCIPfeastol(scip) )
         continue;

      enfoposs[nenfos] = j;
      ++nenfos;

      /* enable probing if tightening the domain could be useful for nlhdlr and violation is above threshold */
      if( sepausesactivity && violation >= nlhdlrdata->minviolprobing )
         doprobing = TRUE;
   }

   if( nenfos == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      SCIPfreeBufferArray(scip, &enfoposs);
      return SCIP_OKAY;
   }

   /* check probing frequency against depth in b&b tree */
   if( nlhdlrdata->probingfreq == -1 || (nlhdlrdata->probingfreq == 0 && SCIPgetDepth(scip) != 0) ||
      (nlhdlrdata->probingfreq > 0 && SCIPgetDepth(scip) % nlhdlrdata->probingfreq != 0)  )
      doprobing = FALSE;

   /* if addbranchscores is TRUE, then we can assume to be in enforcement and not in separation */
   if( nlhdlrdata->probingonlyinsepa && addbranchscores )
      doprobing = FALSE;

   /* disable probing if already being in probing or if in a subscip */
   if( SCIPinProbing(scip) || SCIPgetSubscipDepth(scip) != 0 )
      doprobing = FALSE;

   nrowpreps = 0;
   *result = SCIP_DIDNOTFIND;
   solcopy = sol;
   stop = FALSE;

   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps2) );
   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );
   SCIP_CALL( SCIPcreateBoolarray(scip, &addedbranchscores2) );

   /* build cuts for every indicator variable */
   for( i = 0; i < nlhdlrexprdata->nindicators && !stop; ++i )
   {
      int v;
      int minidx;
      int maxidx;
      int r;
      SCIP_VAR** probingvars;
      SCIP_INTERVAL* probingdoms;
      int nprobingvars;
      SCIP_Bool doprobingind;
      SCIP_Real indval;
      SCIP_Real solval;
      SCIP_Bool adjrefpoint;

      indicator = nlhdlrexprdata->indicators[i];
      probingvars = NULL;
      probingdoms = NULL;
      nprobingvars = 0;
      doprobingind = doprobing;
      solval = SCIPgetSolVal(scip, solcopy, indicator);
      adjrefpoint = nlhdlrdata->adjrefpoint && !SCIPisFeasEQ(scip, solval, 1.0);

      SCIP_CALL( analyseOnoffBounds(scip, nlhdlrdata, nlhdlrexprdata, indicator, &probingvars, &probingdoms,
            &nprobingvars, &doprobingind, result) );

      /* don't add perspective cuts for fixed indicators since there is no use for perspectivy */
      if( SCIPvarGetLbLocal(indicator) >= 0.5 )
      {
         assert(!doprobingind);
         continue;
      }
      if( SCIPvarGetUbLocal(indicator) <= 0.5 )
      { /* this case is stronger as it implies that everything is fixed;
         * therefore we are now happy
         */
         assert(!doprobingind);
         goto TERMINATE;
      }

      if( doprobingind )
      {
         SCIP_Bool propagate;
         SCIP_Bool cutoff_probing;
         SCIP_Bool cutoff;
         SCIP_Bool fixed;

#ifndef NDEBUG
         SCIP_Real* solvals;
         SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nlhdlrexprdata->nvars) );
         for( v = 0; v < nlhdlrexprdata->nvars; ++v )
         {
            solvals[v] = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[v]);
         }
#endif

         propagate = SCIPgetDepth(scip) == 0;

         SCIP_CALL( startProbing(scip, nlhdlrdata, nlhdlrexprdata, indicator, probingvars, probingdoms, nprobingvars,
               sol, &solcopy, &cutoff_probing) );

#ifndef NDEBUG
         for( v = 0; v < nlhdlrexprdata->nvars; ++v )
         {
            assert(solvals[v] == SCIPgetSolVal(scip, solcopy, nlhdlrexprdata->vars[v])); /*lint !e777*/
         }
         SCIPfreeBufferArray(scip, &solvals);
#endif

         if( propagate )
         { /* we are in the root node and startProbing did propagation */
            /* probing propagation might have detected infeasibility */
            if( cutoff_probing )
            {
               /* indicator == 1 is infeasible -> set indicator to 0 */
               SCIPfreeBufferArrayNull(scip, &probingvars);
               SCIPfreeBufferArrayNull(scip, &probingdoms);

               SCIP_CALL( SCIPendProbing(scip) );

               SCIP_CALL( SCIPfixVar(scip, indicator, 0.0, &cutoff, &fixed) );

               if( cutoff )
               {
                  *result = SCIP_CUTOFF;
                  goto TERMINATE;
               }

               continue;
            }

            /* probing propagation in the root node can provide better on/off bounds */
            SCIP_CALL( tightenOnBounds(nlhdlrexprdata, nlhdlrdata->scvars, indicator) );
         }
      }

      if( adjrefpoint )
      {
         /* make sure that when we adjust the point, we don't divide by something too close to 0.0 */
         indval = MAX(solval, 0.1);

         /* create an adjusted point x^adj = (x* - x0) / z* + x0 */
         SCIP_CALL( SCIPcreateSol(scip, &soladj, NULL) );
         for( v = 0; v < nlhdlrexprdata->nvars; ++v )
         {
            if( SCIPvarGetStatus(nlhdlrexprdata->vars[v]) == SCIP_VARSTATUS_FIXED )
               continue;

            scvdata = getSCVarDataInd(nlhdlrdata->scvars, nlhdlrexprdata->vars[v], indicator, &pos);

            /* a non-semicontinuous variable must be linear in expr; skip it */
            if( scvdata == NULL )
               continue;

            SCIP_CALL( SCIPsetSolVal(scip, soladj, nlhdlrexprdata->vars[v],
                  (SCIPgetSolVal(scip, solcopy, nlhdlrexprdata->vars[v]) - scvdata->vals0[pos]) / indval
                  + scvdata->vals0[pos]) );
         }
         for( v = 0; v < nlhdlrexprdata->nindicators; ++v )
         {
            if( SCIPvarGetStatus(nlhdlrexprdata->indicators[v]) == SCIP_VARSTATUS_FIXED )
               continue;

            SCIP_CALL( SCIPsetSolVal(scip, soladj, nlhdlrexprdata->indicators[v],
                  SCIPgetSolVal(scip, solcopy, nlhdlrexprdata->indicators[v])) );
         }
         if( SCIPvarGetStatus(auxvar) != SCIP_VARSTATUS_FIXED )
         {
            SCIP_CALL( SCIPsetSolVal(scip, soladj, auxvar, SCIPgetSolVal(scip, solcopy, auxvar)) );
         }
      }

      /* use cuts from every suitable nlhdlr */
      for( j = 0; j < nenfos; ++j )
      {
         SCIP_Bool addedbranchscores2j;
         SCIP_NLHDLR* nlhdlr2;
         SCIP_NLHDLREXPRDATA* nlhdlr2exprdata;
         SCIP_Real nlhdlr2auxvalue;
         SCIP_Bool success2;

         SCIPgetExprEnfoDataNonlinear(expr, enfoposs[j], &nlhdlr2, &nlhdlr2exprdata, NULL, NULL, NULL, &nlhdlr2auxvalue);
         assert(SCIPnlhdlrHasEstimate(nlhdlr2) && nlhdlr2 != nlhdlr);

         SCIPdebugMsg(scip, "asking nonlinear handler %s to %sestimate\n", SCIPnlhdlrGetName(nlhdlr2), overestimate ? "over" : "under");

         /* ask the nonlinear handler for an estimator */
         if( adjrefpoint )
         {
            SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr2, expr, nlhdlr2exprdata, &nlhdlr2auxvalue, soladj) );

            /* coverity[copy_paste_error] */
            SCIP_CALL( SCIPnlhdlrEstimate(scip, conshdlr, nlhdlr2, expr,
                  nlhdlr2exprdata, soladj,
                  nlhdlr2auxvalue, overestimate, SCIPgetSolVal(scip, solcopy, auxvar),
                  addbranchscores, rowpreps2, &success2, &addedbranchscores2j) );
         }
         else
         {
            SCIP_CALL( SCIPnlhdlrEstimate(scip, conshdlr, nlhdlr2, expr,
                  nlhdlr2exprdata, solcopy,
                  nlhdlr2auxvalue, overestimate, SCIPgetSolVal(scip, solcopy, auxvar),
                  addbranchscores, rowpreps2, &success2, &addedbranchscores2j) );
         }

         minidx = SCIPgetPtrarrayMinIdx(scip, rowpreps2);
         maxidx = SCIPgetPtrarrayMaxIdx(scip, rowpreps2);

         assert((success2 && minidx <= maxidx) || (!success2 && minidx > maxidx));

         /* perspectivy all cuts from nlhdlr2 and add them to rowpreps */
         for( r = minidx; r <= maxidx; ++r )
         {
            SCIP_Real maxcoef;
            SCIP_Real* rowprepcoefs;
            SCIP_VAR** rowprepvars;

            rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps2, r);
            assert(rowprep != NULL);

#ifdef SCIP_DEBUG
            SCIPinfoMessage(scip, NULL, "rowprep for expr ");
            SCIPprintExpr(scip, expr, NULL);
            SCIPinfoMessage(scip, NULL, "rowprep before perspectivy is: \n");
            SCIPprintRowprep(scip, rowprep, NULL);
#endif

            /* given a rowprep: sum aixi + sum biyi + c, where xi are semicontinuous variables and yi are
             * non-semicontinuous variables (which appear in expr linearly, which detect must have ensured),
             * perspectivy the semicontinuous part by adding (1-z)(g0 - c - sum aix0i) (the constant is
             * treated as belonging to the semicontinuous part)
             */

            /* we want cst0 = g0 - c - sum aix0i; first add g0 - c */
            cst0 = nlhdlrexprdata->exprvals0[i] + SCIProwprepGetSide(rowprep);

            maxcoef = 0.0;
            rowprepcoefs = SCIProwprepGetCoefs(rowprep);
            rowprepvars = SCIProwprepGetVars(rowprep);

            for( v = 0; v < SCIProwprepGetNVars(rowprep); ++v )
            {
               if( REALABS( rowprepcoefs[v]) > maxcoef )
               {
                  maxcoef = REALABS(rowprepcoefs[v]);
               }

               scvdata = getSCVarDataInd(nlhdlrdata->scvars, rowprepvars[v], indicator, &pos);

               /* a non-semicontinuous variable must be linear in expr; skip it */
               if( scvdata == NULL )
                  continue;

               cst0 -= rowprepcoefs[v] * scvdata->vals0[pos];
            }

            /* only perspectivy when the absolute value of cst0 is not too small
             * TODO on ex1252a there was cst0=0 - ok to still use the cut?
            */
            if( cst0 == 0.0 || maxcoef / REALABS(cst0) <= 10.0 / SCIPfeastol(scip) )
            {
               /* update the rowprep by adding cst0 - cst0*z */
               SCIProwprepAddConstant(rowprep, cst0);
               SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, indicator, -cst0));
            }
            else
            {
               SCIPfreeRowprep(scip, &rowprep);
               continue;
            }

            SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0));

            SCIPdebugMsg(scip, "rowprep after perspectivy is: \n");
#ifdef SCIP_DEBUG
            SCIPprintRowprep(scip, rowprep, NULL);
#endif

            SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, nrowpreps, rowprep) );
            SCIP_CALL( SCIPsetBoolarrayVal(scip, addedbranchscores2, nrowpreps, addedbranchscores2j) );
            ++nrowpreps;
         }

         SCIP_CALL( SCIPclearPtrarray(scip, rowpreps2) );
      }

      if( adjrefpoint )
      {
         SCIP_CALL( SCIPfreeSol(scip, &soladj) );
      }

      if( doprobingind )
      {
         SCIP_CALL( SCIPendProbing(scip) );
      }

      /* add all cuts found for indicator i */
      for( r = SCIPgetPtrarrayMinIdx(scip, rowpreps); r <= SCIPgetPtrarrayMaxIdx(scip, rowpreps) && !stop; ++r )
      {
         SCIP_RESULT resultr;

#ifdef SCIP_DEBUG
         SCIPprintRowprep(scip, rowprep, NULL);
#endif
         rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, r);
         resultr = SCIP_DIDNOTFIND;

         (void) strcat(SCIProwprepGetName(rowprep), "_persp_indicator_");
         (void) strcat(SCIProwprepGetName(rowprep), SCIPvarGetName(indicator));

         SCIP_CALL( SCIPprocessRowprepNonlinear(scip, nlhdlr, cons, expr, rowprep, overestimate, auxvar, auxvalue,
               allowweakcuts, SCIPgetBoolarrayVal(scip, addedbranchscores2, r), addbranchscores, solcopy, &resultr) );

         if( resultr == SCIP_SEPARATED )
            *result = SCIP_SEPARATED;
         else if( resultr == SCIP_CUTOFF )
         {
            *result = SCIP_CUTOFF;
            stop = TRUE;
         }
         else if( resultr == SCIP_BRANCHED )
         {
            if( *result != SCIP_SEPARATED && *result != SCIP_REDUCEDDOM )
               *result = SCIP_BRANCHED;
         }
         else if( resultr != SCIP_DIDNOTFIND )
         {
            SCIPerrorMessage("estimate called by perspective nonlinear handler returned invalid result <%d>\n", resultr);
            return SCIP_INVALIDRESULT;
         }
      }

      /* free all rowpreps for indicator i */
      for( r = SCIPgetPtrarrayMinIdx(scip, rowpreps); r <= SCIPgetPtrarrayMaxIdx(scip, rowpreps); ++r )
      {
         rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, r);
         SCIPfreeRowprep(scip, &rowprep);
      }

      SCIPfreeBufferArrayNull(scip, &probingvars);
      SCIPfreeBufferArrayNull(scip, &probingdoms);
      SCIP_CALL( SCIPclearPtrarray(scip, rowpreps) );
   }

TERMINATE:
   SCIP_CALL( SCIPfreeBoolarray(scip, &addedbranchscores2) );
   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps2) );
   if( solcopy != sol )
   {
      SCIP_CALL( SCIPfreeSol(scip, &solcopy) );
   }
   SCIPfreeBufferArray(scip, &enfoposs);

   return SCIP_OKAY;
}


/*
 * nonlinear handler specific interface methods
 */

/** includes perspective nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrPerspective(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   /* create nonlinear handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   BMSclearMemory(nlhdlrdata);

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectPerspective, nlhdlrEvalauxPerspective, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/maxproprounds",
           "maximal number of propagation rounds in probing",
           &nlhdlrdata->maxproprounds, FALSE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/mindomreduction",
           "minimal relative reduction in a variable's domain for applying probing",
           &nlhdlrdata->mindomreduction, FALSE, DEFAULT_MINDOMREDUCTION, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/minviolprobing",
           "minimal violation w.r.t. auxiliary variables for applying probing",
           &nlhdlrdata->minviolprobing, FALSE, DEFAULT_MINVIOLPROBING, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/probingonlyinsepa",
           "whether to do probing only in separation",
           &nlhdlrdata->probingonlyinsepa, FALSE, DEFAULT_PROBINGONLYINSEPA, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "nlhdlr/" NLHDLR_NAME "/probingfreq",
           "probing frequency (-1 - no probing, 0 - root node only)",
           &nlhdlrdata->probingfreq, FALSE, DEFAULT_PROBINGFREQ, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/convexonly",
           "whether perspective cuts are added only for convex expressions",
           &nlhdlrdata->convexonly, FALSE, DEFAULT_CONVEXONLY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/tightenbounds",
           "whether variable semicontinuity is used to tighten variable bounds",
           &nlhdlrdata->tightenbounds, FALSE, DEFAULT_TIGHTENBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/adjrefpoint",
           "whether to adjust the reference point",
           &nlhdlrdata->adjrefpoint, FALSE, DEFAULT_ADJREFPOINT, NULL, NULL) );

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrPerspective);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrdataPerspective);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataPerspective);
   SCIPnlhdlrSetInitExit(nlhdlr, NULL, nlhdlrExitPerspective);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaPerspective, nlhdlrEnfoPerspective, NULL, NULL);

   return SCIP_OKAY;
}
