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

/**@file   expr.c
 * @ingroup OTHER_CFILES
 * @brief  functions for algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#include <assert.h>
#include <ctype.h>

#include "scip/expr.h"
#include "scip/struct_expr.h"
#include "scip/pub_misc.h"
#include "scip/clock.h"
#include "scip/set.h"
#include "scip/pub_var.h"
#include "scip/sol.h"
#include "scip/tree.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/nlpi_ipopt.h" /* for LAPACK */

/*lint -e440*/
/*lint -e441*/
/*lint -e777*/

/*
 * Data structures
 */

/** printing to file data */
struct SCIP_ExprPrintData
{
   FILE*                 file;               /**< file to print to */
   SCIP_EXPRITER*        iterator;           /**< iterator to use */
   SCIP_Bool             closefile;          /**< whether file need to be closed when finished printing */
   SCIP_HASHMAP*         leaveexprs;         /**< hashmap storing leave (no children) expressions */
   SCIP_EXPRPRINT_WHAT   whattoprint;        /**< flags that indicate what to print for each expression */
};

/*
 * Local methods
 */

/** frees an expression */
static
SCIP_RETCODE freeExpr(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr                /**< pointer to free the expression */
   )
{
   assert(expr != NULL);
   assert(*expr != NULL);
   assert((*expr)->nuses == 1);
   assert((*expr)->quaddata == NULL);
   assert((*expr)->ownerdata == NULL);

   /* free children array, if any */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*expr)->children, (*expr)->childrensize);

   BMSfreeBlockMemory(blkmem, expr);
   assert(*expr == NULL);

   return SCIP_OKAY;
}

/*
 * quadratic representation of expression
 */

/** first time seen quadratically and
 * seen before linearly --> --nlinterms; assign 2; ++nquadterms
 * not seen before linearly --> assing 1; ++nquadterms
 *
 * seen before --> assign += 1
 */
static
SCIP_RETCODE quadDetectProcessExpr(
   SCIP_EXPR*            expr,               /**< the expression */
   SCIP_HASHMAP*         seenexpr,           /**< hash map */
   int*                  nquadterms,         /**< number of quadratic terms */
   int*                  nlinterms           /**< number of linear terms */
   )
{
   if( SCIPhashmapExists(seenexpr, (void*)expr) )
   {
      int nseen = SCIPhashmapGetImageInt(seenexpr, (void*)expr);

      if( nseen < 0 )
      {
         /* only seen linearly before */
         assert(nseen == -1);

         --*nlinterms;
         ++*nquadterms;
         SCIP_CALL( SCIPhashmapSetImageInt(seenexpr, (void*)expr, 2) );
      }
      else
      {
         assert(nseen > 0);
         SCIP_CALL( SCIPhashmapSetImageInt(seenexpr, (void*)expr, nseen + 1) );
      }
   }
   else
   {
      ++*nquadterms;
      SCIP_CALL( SCIPhashmapInsertInt(seenexpr, (void*)expr, 1) );
   }

   return SCIP_OKAY;
}

/** returns a quadexprterm that contains the expr
 *
 * it either finds one that already exists or creates a new one
 */
static
SCIP_RETCODE quadDetectGetQuadexprterm(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< the expression */
   SCIP_HASHMAP*         expr2idx,           /**< map: expr to index in quadexpr->quadexprterms */
   SCIP_HASHMAP*         seenexpr,           /**< map: expr to number of times it was seen */
   SCIP_QUADEXPR*        quadexpr,           /**< data of quadratic representation of expression */
   SCIP_QUADEXPR_QUADTERM** quadexprterm     /**< buffer to store quadexprterm */
   )
{
   assert(expr != NULL);
   assert(expr2idx != NULL);
   assert(quadexpr != NULL);
   assert(quadexprterm != NULL);

   if( SCIPhashmapExists(expr2idx, (void*)expr) )
   {
      *quadexprterm = &quadexpr->quadexprterms[SCIPhashmapGetImageInt(expr2idx, (void*)expr)];
      assert((*quadexprterm)->expr == expr);
   }
   else
   {
      SCIP_CALL( SCIPhashmapInsertInt(expr2idx, expr, quadexpr->nquadexprs) );
      *quadexprterm = &quadexpr->quadexprterms[quadexpr->nquadexprs];
      ++quadexpr->nquadexprs;

      (*quadexprterm)->expr = expr;
      (*quadexprterm)->sqrcoef = 0.0;
      (*quadexprterm)->sqrexpr = NULL;
      (*quadexprterm)->lincoef = 0.0;
      (*quadexprterm)->nadjbilin = 0;
      (*quadexprterm)->adjbilinsize = SCIPhashmapGetImageInt(seenexpr, (void*)expr);
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*quadexprterm)->adjbilin, (*quadexprterm)->adjbilinsize) );
   }

   return SCIP_OKAY;
}


/** evaluate and forward-differentiate expression
 *
 * also initializes derivative and bardot to 0.0
 */
static
SCIP_RETCODE evalAndDiff(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   SCIP_Longint          soltag,             /**< tag that uniquely identifies the solution (with its values), or 0. */
   SCIP_SOL*             direction           /**< direction for directional derivative */
   )
{
   SCIP_EXPRITER* it;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);

   /* assume we'll get a domain error, so we don't have to get this expr back if we abort the iteration
    * if there is no domain error, then we will overwrite the evalvalue in the last leaveexpr stage
    */
   expr->evalvalue = SCIP_INVALID;
   expr->evaltag = soltag;
   expr->dot = SCIP_INVALID;

   /* start a new difftag */
   ++stat->exprlastdifftag;

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      /* evaluate expression only if necessary */
      if( soltag == 0 || expr->evaltag != soltag )
      {
         SCIP_CALL( SCIPexprhdlrEvalExpr(expr->exprhdlr, set, NULL, expr, &expr->evalvalue, NULL, sol) );

         expr->evaltag = soltag;
      }

      if( expr->evalvalue == SCIP_INVALID )
         break;

      if( expr->difftag != stat->exprlastdifftag )
      {
         /* compute forward diff */
         SCIP_CALL( SCIPexprhdlrFwDiffExpr(expr->exprhdlr, set, expr, &expr->dot, direction) );

         if( expr->dot == SCIP_INVALID )
            break;

         expr->derivative = 0.0;
         expr->bardot = 0.0;
         expr->difftag = stat->exprlastdifftag;
      }
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}


/*
 * Public methods
 */

/* Undo the defines from pub_expr.h, which exist if NDEBUG is defined. */
#ifdef NDEBUG
#undef SCIPexprhdlrSetCopyFreeHdlr
#undef SCIPexprhdlrSetCopyFreeData
#undef SCIPexprhdlrSetPrint
#undef SCIPexprhdlrSetParse
#undef SCIPexprhdlrSetCurvature
#undef SCIPexprhdlrSetMonotonicity
#undef SCIPexprhdlrSetIntegrality
#undef SCIPexprhdlrSetHash
#undef SCIPexprhdlrSetCompare
#undef SCIPexprhdlrSetDiff
#undef SCIPexprhdlrSetIntEval
#undef SCIPexprhdlrSetSimplify
#undef SCIPexprhdlrSetReverseProp
#undef SCIPexprhdlrSetEstimate
#undef SCIPexprhdlrGetName
#undef SCIPexprhdlrGetDescription
#undef SCIPexprhdlrGetPrecedence
#undef SCIPexprhdlrGetData
#undef SCIPexprhdlrHasPrint
#undef SCIPexprhdlrHasBwdiff
#undef SCIPexprhdlrHasFwdiff
#undef SCIPexprhdlrHasIntEval
#undef SCIPexprhdlrHasEstimate
#undef SCIPexprhdlrHasInitEstimates
#undef SCIPexprhdlrHasSimplify
#undef SCIPexprhdlrHasCurvature
#undef SCIPexprhdlrHasMonotonicity
#undef SCIPexprhdlrHasReverseProp
#undef SCIPexprhdlrGetNCreated
#undef SCIPexprhdlrGetNIntevalCalls
#undef SCIPexprhdlrGetIntevalTime
#undef SCIPexprhdlrGetNReversepropCalls
#undef SCIPexprhdlrGetReversepropTime
#undef SCIPexprhdlrGetNCutoffs
#undef SCIPexprhdlrGetNDomainReductions
#undef SCIPexprhdlrIncrementNDomainReductions
#undef SCIPexprhdlrGetNEstimateCalls
#undef SCIPexprhdlrGetEstimateTime
#undef SCIPexprhdlrGetNBranchings
#undef SCIPexprhdlrIncrementNBranchings
#undef SCIPexprhdlrGetNSimplifyCalls
#undef SCIPexprhdlrGetSimplifyTime
#undef SCIPexprhdlrGetNSimplifications
#endif

/** create expression handler */
SCIP_RETCODE SCIPexprhdlrCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRHDLR**       exprhdlr,           /**< buffer where to store created expression handler */
   const char*           name,               /**< name of expression handler (must not be NULL) */
   const char*           desc,               /**< description of expression handler (can be NULL) */
   unsigned int          precedence,         /**< precedence of expression operation (used for printing) */
   SCIP_DECL_EXPREVAL((*eval)),              /**< point evaluation callback (must not be NULL) */
   SCIP_EXPRHDLRDATA*    data                /**< data of expression handler (can be NULL) */
   )
{
   assert(exprhdlr != NULL);
   assert(name != NULL);

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, exprhdlr) );

   SCIP_ALLOC( BMSduplicateMemoryArray(&(*exprhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_ALLOC( BMSduplicateMemoryArray(&(*exprhdlr)->desc, desc, strlen(desc)+1) );
   }

   (*exprhdlr)->precedence = precedence;
   (*exprhdlr)->eval = eval;
   (*exprhdlr)->data = data;

   /* create clocks */
   SCIP_CALL( SCIPclockCreate(&(*exprhdlr)->estimatetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*exprhdlr)->intevaltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*exprhdlr)->proptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*exprhdlr)->simplifytime, SCIP_CLOCKTYPE_DEFAULT) );

   return SCIP_OKAY;
}

/** frees expression handler */
SCIP_RETCODE SCIPexprhdlrFree(
   SCIP_EXPRHDLR**       exprhdlr,           /**< pointer to expression handler to be freed */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   if( (*exprhdlr)->freehdlr != NULL )
   {
      SCIP_CALL( (*exprhdlr)->freehdlr(set->scip, *exprhdlr, &(*exprhdlr)->data) );
   }

   /* free clocks */
   SCIPclockFree(&(*exprhdlr)->simplifytime);
   SCIPclockFree(&(*exprhdlr)->intevaltime);
   SCIPclockFree(&(*exprhdlr)->proptime);
   SCIPclockFree(&(*exprhdlr)->estimatetime);

   BMSfreeMemoryArrayNull(&(*exprhdlr)->desc);
   BMSfreeMemoryArray(&(*exprhdlr)->name);

   BMSfreeBlockMemory(blkmem, exprhdlr);

   return SCIP_OKAY;
}

/**@addtogroup PublicExprHandlerMethods
 * @{
 */

/** set the expression handler callbacks to copy and free an expression handler */
void SCIPexprhdlrSetCopyFreeHdlr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYHDLR((*copyhdlr)),      /**< handler copy callback (can be NULL) */
   SCIP_DECL_EXPRFREEHDLR((*freehdlr))       /**< handler free callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->copyhdlr = copyhdlr;
   exprhdlr->freehdlr = freehdlr;
}

/** set the expression handler callbacks to copy and free expression data */
void SCIPexprhdlrSetCopyFreeData(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOPYDATA((*copydata)),      /**< expression data copy callback (can be NULL for expressions
                                                  without data) */
   SCIP_DECL_EXPRFREEDATA((*freedata))       /**< expression data free callback (can be NULL if data does not
                                                  need to be freed) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->copydata = copydata;
   exprhdlr->freedata = freedata;
}

/** set the print callback of an expression handler */
void SCIPexprhdlrSetPrint(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPRINT((*print))             /**< print callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->print = print;
}

/** set the parse callback of an expression handler */
void SCIPexprhdlrSetParse(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRPARSE((*parse))             /**< parse callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->parse = parse;
}

/** set the curvature detection callback of an expression handler */
void SCIPexprhdlrSetCurvature(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCURVATURE((*curvature))     /**< curvature detection callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->curvature = curvature;
}

/** set the monotonicity detection callback of an expression handler */
void SCIPexprhdlrSetMonotonicity(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRMONOTONICITY((*monotonicity)) /**< monotonicity detection callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->monotonicity = monotonicity;
}

/** set the integrality detection callback of an expression handler */
void SCIPexprhdlrSetIntegrality(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEGRALITY((*integrality)) /**< integrality detection callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->integrality = integrality;
}

/** set the hash callback of an expression handler */
void SCIPexprhdlrSetHash(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRHASH((*hash))               /**< hash callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->hash = hash;
}

/** set the compare callback of an expression handler */
void SCIPexprhdlrSetCompare(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRCOMPARE((*compare))         /**< compare callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->compare = compare;
}

/** set differentiation callbacks of an expression handler */
void SCIPexprhdlrSetDiff(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRBWDIFF((*bwdiff)),          /**< backward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRFWDIFF((*fwdiff)),          /**< forward derivative evaluation callback (can be NULL) */
   SCIP_DECL_EXPRBWFWDIFF((*bwfwdiff))       /**< backward-forward derivative evaluation callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->bwdiff = bwdiff;
   exprhdlr->fwdiff = fwdiff;
   exprhdlr->bwfwdiff = bwfwdiff;
}

/** set the interval evaluation callback of an expression handler */
void SCIPexprhdlrSetIntEval(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINTEVAL((*inteval))         /**< interval evaluation callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->inteval = inteval;
}

/** set the simplify callback of an expression handler */
void SCIPexprhdlrSetSimplify(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRSIMPLIFY((*simplify))       /**< simplify callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->simplify = simplify;
}

/** set the reverse propagation callback of an expression handler */
void SCIPexprhdlrSetReverseProp(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->reverseprop = reverseprop;
}

/** set the estimation callbacks of an expression handler */
void SCIPexprhdlrSetEstimate(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_DECL_EXPRINITESTIMATES((*initestimates)), /**< initial estimators callback (can be NULL) */
   SCIP_DECL_EXPRESTIMATE((*estimate))       /**< estimator callback (can be NULL) */
   )
{
   assert(exprhdlr != NULL);

   exprhdlr->initestimates = initestimates;
   exprhdlr->estimate = estimate;
}

/** gives the name of an expression handler */
const char* SCIPexprhdlrGetName(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->name;
}

/** gives the description of an expression handler (can be NULL) */
const char* SCIPexprhdlrGetDescription(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->desc;
}

/** gives the precedence of an expression handler */
unsigned int SCIPexprhdlrGetPrecedence(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->precedence;
}

/** gives the data of an expression handler */
SCIP_EXPRHDLRDATA* SCIPexprhdlrGetData(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->data;
}

/** returns whether expression handler implements the print callback */
SCIP_Bool SCIPexprhdlrHasPrint(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->print != NULL;
}

/** returns whether expression handler implements the backward differentiation callback */
SCIP_Bool SCIPexprhdlrHasBwdiff(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->bwdiff != NULL;
}

/** returns whether expression handler implements the forward differentiation callback */
SCIP_Bool SCIPexprhdlrHasFwdiff(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->fwdiff != NULL;
}

/** returns whether expression handler implements the interval evaluation callback */
SCIP_Bool SCIPexprhdlrHasIntEval(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->inteval != NULL;
}

/** returns whether expression handler implements the estimator callback */
SCIP_Bool SCIPexprhdlrHasEstimate(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->estimate != NULL;
}

/** returns whether expression handler implements the initial estimators callback */
SCIP_Bool SCIPexprhdlrHasInitEstimates(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->initestimates != NULL;
}

/** returns whether expression handler implements the simplification callback */
SCIP_Bool SCIPexprhdlrHasSimplify(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->simplify != NULL;
}

/** returns whether expression handler implements the curvature callback */
SCIP_Bool SCIPexprhdlrHasCurvature(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->curvature != NULL;
}

/** returns whether expression handler implements the monotonicity callback */
SCIP_Bool SCIPexprhdlrHasMonotonicity(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->monotonicity != NULL;
}

/** returns whether expression handler implements the reverse propagation callback */
SCIP_Bool SCIPexprhdlrHasReverseProp(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->reverseprop != NULL;
}

/** compares two expression handler w.r.t. their name */
SCIP_DECL_SORTPTRCOMP(SCIPexprhdlrComp)
{
   return strcmp(((SCIP_EXPRHDLR*)elem1)->name, ((SCIP_EXPRHDLR*)elem2)->name);
}

/** gets number of times an expression has been created with given expression handler */
unsigned int SCIPexprhdlrGetNCreated(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->ncreated;
}

/** gets number of times the interval evaluation callback was called */
SCIP_Longint SCIPexprhdlrGetNIntevalCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->nintevalcalls;
}

/** gets time spend in interval evaluation callback */
SCIP_Real SCIPexprhdlrGetIntevalTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return SCIPclockGetTime(exprhdlr->intevaltime);
}

/** gets number of times the reverse propagation callback was called */
SCIP_Longint SCIPexprhdlrGetNReversepropCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->npropcalls;
}

/** gets time spend in reverse propagation callback */
SCIP_Real SCIPexprhdlrGetReversepropTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return SCIPclockGetTime(exprhdlr->proptime);
}

/** gets number of times an empty interval was found in reverse propagation */
SCIP_Longint SCIPexprhdlrGetNCutoffs(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->ncutoffs;
}

/** gets number of times a bound reduction was found in reverse propagation (and accepted by caller) */
SCIP_Longint SCIPexprhdlrGetNDomainReductions(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->ndomreds;
}

/** increments the domain reductions count of an expression handler */
void SCIPexprhdlrIncrementNDomainReductions(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   int                   nreductions         /**< number of reductions to add to counter */
   )
{
   assert(exprhdlr != NULL);
   assert(nreductions >= 0);

   exprhdlr->ndomreds += nreductions;
}

/** gets number of times the estimation callback was called */
SCIP_Longint SCIPexprhdlrGetNEstimateCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->nestimatecalls;
}

/** gets time spend in estimation callback */
SCIP_Real SCIPexprhdlrGetEstimateTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return SCIPclockGetTime(exprhdlr->estimatetime);
}

/** gets number of times branching candidates reported by of this expression handler were used to
 * assemble branching candidates
 *
 * that is, how often did we consider branching on a child of this expression
 */
SCIP_Longint SCIPexprhdlrGetNBranchings(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->nbranchscores;
}

/** increments the branching candidates count of an expression handler */
void SCIPexprhdlrIncrementNBranchings(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   ++exprhdlr->nbranchscores;
}

/** gets number of times the simplify callback was called */
SCIP_Longint SCIPexprhdlrGetNSimplifyCalls(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->nsimplifycalls;
}

/** gets time spend in simplify callback */
SCIP_Real SCIPexprhdlrGetSimplifyTime(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return SCIPclockGetTime(exprhdlr->simplifytime);
}

/** gets number of times the simplify callback found a simplification */
SCIP_Longint SCIPexprhdlrGetNSimplifications(
   SCIP_EXPRHDLR*        exprhdlr            /**< expression handler */
   )
{
   assert(exprhdlr != NULL);

   return exprhdlr->nsimplified;
}

/** @} */

/** copies the given expression handler to a new scip */
SCIP_RETCODE SCIPexprhdlrCopyInclude(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             targetset           /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(exprhdlr != NULL);
   assert(targetset != NULL);
   assert(targetset->scip != NULL);

   if( exprhdlr->copyhdlr != NULL )
   {
      SCIPsetDebugMsg(targetset, "including expression handler <%s> in subscip %p\n",
            SCIPexprhdlrGetName(exprhdlr), (void*)targetset->scip);
      SCIP_CALL( exprhdlr->copyhdlr(targetset->scip, exprhdlr) );
   }
   else
   {
      SCIPsetDebugMsg(targetset, "expression handler <%s> cannot be copied to subscip %p due "
            "to missing copyhdlr callback\n", SCIPexprhdlrGetName(exprhdlr), (void*)targetset->scip);
   }

   return SCIP_OKAY;
}

/** initialization of expression handler (resets statistics) */
void SCIPexprhdlrInit(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(exprhdlr != NULL);

   if( set->misc_resetstat )
   {
      exprhdlr->ncreated = 0;
      exprhdlr->nestimatecalls = 0;
      exprhdlr->nintevalcalls = 0;
      exprhdlr->npropcalls = 0;
      exprhdlr->ncutoffs = 0;
      exprhdlr->ndomreds = 0;
      exprhdlr->nbranchscores = 0;
      exprhdlr->nsimplifycalls = 0;
      exprhdlr->nsimplified = 0;

      SCIPclockReset(exprhdlr->estimatetime);
      SCIPclockReset(exprhdlr->intevaltime);
      SCIPclockReset(exprhdlr->proptime);
      SCIPclockReset(exprhdlr->simplifytime);
   }
}

/** calls the print callback of an expression handler
 *
 * The method prints an expression.
 * It is called while iterating over the expression graph at different stages.
 *
 * @see SCIP_DECL_EXPRPRINT
 */
SCIP_RETCODE SCIPexprhdlrPrintExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRITER_STAGE   stage,              /**< stage of expression iteration */
   int                   currentchild,       /**< index of current child if in stage visitingchild or visitedchild */
   unsigned int          parentprecedence,   /**< precedence of parent */
   FILE*                 file                /**< the file to print to */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(messagehdlr != NULL);

   if( SCIPexprhdlrHasPrint(exprhdlr) )
   {
      SCIP_CALL( exprhdlr->print(set->scip, expr, stage, currentchild, parentprecedence, file) );
   }
   else
   {
      /* default: <hdlrname>(<child1>, <child2>, ...) */
      switch( stage )
      {
         case SCIP_EXPRITER_ENTEREXPR :
         {
            SCIPmessageFPrintInfo(messagehdlr, file, "%s", SCIPexprhdlrGetName(expr->exprhdlr));
            if( expr->nchildren > 0 )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, "(");
            }
            break;
         }

         case SCIP_EXPRITER_VISITEDCHILD :
         {
            assert(currentchild >= 0);
            assert(currentchild < expr->nchildren);
            if( currentchild < expr->nchildren-1 )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, ", ");
            }
            else
            {
               SCIPmessageFPrintInfo(messagehdlr, file, ")");
            }

            break;
         }

         case SCIP_EXPRITER_VISITINGCHILD :
         case SCIP_EXPRITER_LEAVEEXPR :
         default:
            break;
      }
   }

   return SCIP_OKAY;
}

/** calls the parse callback of an expression handler
 *
 * The method parses an expression.
 * It should be called when parsing an expression and an operator with the expr handler name is found.
 *
 * @see SCIP_DECL_EXPRPARSE
 */
SCIP_RETCODE SCIPexprhdlrParseExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           string,             /**< string containing expression to be parse */
   const char**          endstring,          /**< buffer to store the position of string after parsing */
   SCIP_EXPR**           expr,               /**< buffer to store the parsed expression */
   SCIP_Bool*            success,            /**< buffer to store whether the parsing was successful or not */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);

   *expr = NULL;

   if( exprhdlr->parse == NULL )
   {
      /* TODO we could just look for a comma separated list of operands and try to initialize the expr with this one?
       * That would be sufficient for sin, cos, exp, log, abs, for example.
       */
      SCIPdebugMessage("Expression handler <%s> has no parsing method.\n", SCIPexprhdlrGetName(exprhdlr));
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* give control to exprhdlr's parser */
   SCIP_CALL( exprhdlr->parse(set->scip, exprhdlr, string, endstring, expr, success, ownercreate, ownercreatedata) );

   assert(*success || (*expr == NULL));

   return SCIP_OKAY;
}

/** calls the curvature check callback of an expression handler
 *
 * @see SCIP_DECL_EXPRCURVATURE
 */
SCIP_RETCODE SCIPexprhdlrCurvatureExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to check the curvature for */
   SCIP_EXPRCURV         exprcurvature,      /**< desired curvature of this expression */
   SCIP_Bool*            success,            /**< buffer to store whether the desired curvature be obtained */
   SCIP_EXPRCURV*        childcurv           /**< array to store required curvature for each child */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(success != NULL);

   *success = FALSE;

   if( exprhdlr->curvature != NULL )
   {
      SCIP_CALL( exprhdlr->curvature(set->scip, expr, exprcurvature, success, childcurv) );
   }

   return SCIP_OKAY;
}

/** calls the monotonicity check callback of an expression handler
 *
 * @see SCIP_DECL_EXPRMONOTONICITY
 */
SCIP_RETCODE SCIPexprhdlrMonotonicityExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to check the monotonicity for */
   int                   childidx,           /**< index of the considered child expression */
   SCIP_MONOTONE*        result              /**< buffer to store the monotonicity */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(result != NULL);

   *result = SCIP_MONOTONE_UNKNOWN;

   /* check whether the expression handler implements the monotonicity callback */
   if( exprhdlr->monotonicity != NULL )
   {
      SCIP_CALL( exprhdlr->monotonicity(set->scip, expr, childidx, result) );
   }

   return SCIP_OKAY;
}

/** calls the integrality check callback of an expression handler
 *
 * @see SCIP_DECL_EXPRINTEGRALITY
 */
SCIP_RETCODE SCIPexprhdlrIntegralityExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to check integrality for */
   SCIP_Bool*            isintegral          /**< buffer to store whether expression is integral */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(isintegral != NULL);

   *isintegral = FALSE;

   /* check whether the expression handler implements the monotonicity callback */
   if( exprhdlr->integrality != NULL )
   {
      SCIP_CALL( exprhdlr->integrality(set->scip, expr, isintegral) );
   }

   return SCIP_OKAY;
}

/** calls the hash callback of an expression handler
 *
 * The method hashes an expression by taking the hashes of its children into account.
 *
 * @see SCIP_DECL_EXPRHASH
 */
SCIP_RETCODE SCIPexprhdlrHashExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be hashed */
   unsigned int*         hashkey,            /**< buffer to store the hash value */
   unsigned int*         childrenhashes      /**< array with hash values of children */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL || expr->nchildren == 0);

   if( expr->exprhdlr->hash != NULL )
   {
      SCIP_CALL( expr->exprhdlr->hash(set->scip, expr, hashkey, childrenhashes) );
   }
   else
   {
      int i;

      /* compute initial hash from expression handler name if callback is not implemented
       * this can lead to more collisions and thus a larger number of expensive expression compare calls
       */
      *hashkey = 0;
      for( i = 0; expr->exprhdlr->name[i] != '\0'; i++ )
         *hashkey += (unsigned int) expr->exprhdlr->name[i]; /*lint !e571*/

      *hashkey = SCIPcalcFibHash((SCIP_Real)*hashkey);

      /* now make use of the hashkeys of the children */
      for( i = 0; i < expr->nchildren; ++i )
         *hashkey ^= childrenhashes[i];
   }

   return SCIP_OKAY;
}

/** calls the compare callback of an expression handler
 *
 * The method receives two expressions, expr1 and expr2, and returns
 * - -1 if expr1 < expr2,
 * - 0  if expr1 = expr2,
 * - 1  if expr1 > expr2.
 *
 * @see SCIP_DECL_EXPRCOMPARE
 */
int SCIPexprhdlrCompareExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr1,              /**< first expression in comparison */
   SCIP_EXPR*            expr2               /**< second expression in comparison */
   )
{
   int i;

   assert(expr1 != NULL);
   assert(expr2 != NULL);
   assert(expr1->exprhdlr == expr2->exprhdlr);

   if( expr1->exprhdlr->compare != NULL )
   {
      /* enforces OR1-OR4 */
      return expr1->exprhdlr->compare(set->scip, expr1, expr2);
   }

   /* enforces OR5: default comparison method of expressions of the same type:
    * expr1 < expr2 if and only if expr1_i = expr2_i for all i < k and expr1_k < expr2_k.
    * if there is no such k, use number of children to decide
    * if number of children is equal, both expressions are equal
    * @note: Warning, this method doesn't know about expression data. So if your expressions have special data,
    * you must implement the compare callback: SCIP_DECL_EXPRCOMPARE
    */
   for( i = 0; i < expr1->nchildren && i < expr2->nchildren; ++i )
   {
      int compareresult = SCIPexprCompare(set, expr1->children[i], expr2->children[i]);
      if( compareresult != 0 )
         return compareresult;
   }

   return expr1->nchildren == expr2->nchildren ? 0 : expr1->nchildren < expr2->nchildren ? -1 : 1;
}

/** calls the evaluation callback of an expression handler
 *
 * The method evaluates an expression by taking the values of its children into account.
 *
 * Further, allows to evaluate w.r.t. given expression and children values instead of those stored in children expressions.
 *
 * @see SCIP_DECL_EXPREVAL
 */
SCIP_RETCODE SCIPexprhdlrEvalExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            val,                /**< buffer to store value of expression */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_SOL*             sol                 /**< solution that is evaluated (can be NULL) */
   )
{
   SCIP_Real* origvals = NULL;

   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(exprhdlr->eval != NULL);
   assert(val != NULL);

   /* temporarily overwrite the evalvalue in all children with values from childrenvals */
   if( childrenvals != NULL && expr->nchildren > 0 )
   {
      int c;

      assert(bufmem != NULL);

      SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &origvals, expr->nchildren) );

      for( c = 0; c < expr->nchildren; ++c )
      {
         origvals[c] = expr->children[c]->evalvalue;
         expr->children[c]->evalvalue = childrenvals[c];
      }
   }

   /* call expression eval callback */
   SCIP_CALL( exprhdlr->eval(set->scip, expr, val, sol) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*val) )
      *val = SCIP_INVALID;

   /* restore original evalvalues in children */
   if( origvals != NULL )
   {
      int c;
      for( c = 0; c < expr->nchildren; ++c )
         expr->children[c]->evalvalue = origvals[c];

      BMSfreeBufferMemoryArray(bufmem, &origvals);
   }

   return SCIP_OKAY;
}

/** calls the backward derivative evaluation callback of an expression handler
 *
 * The method should compute the partial derivative of expr w.r.t its child at childidx.
 * That is, it returns
 * \f[
 *   \frac{\partial \text{expr}}{\partial \text{child}_{\text{childidx}}}
 * \f]
 *
 * Further, allows to differentiate w.r.t. given expression and children values instead of those stored in children expressions.
 *
 * @see SCIP_DECL_EXPRBWDIFF
 */
SCIP_RETCODE SCIPexprhdlrBwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   int                   childidx,           /**< index of the child */
   SCIP_Real*            derivative,         /**< buffer to store the partial derivative w.r.t. the i-th children */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should be used */
   SCIP_Real             exprval             /**< value for expression, used only if childrenvals is not NULL */
   )
{
   SCIP_Real* origchildrenvals;
   SCIP_Real origexprval;
   int c;

   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(derivative != NULL);

   if( exprhdlr->bwdiff == NULL )
   {
      *derivative = SCIP_INVALID;
      return SCIP_OKAY;
   }

   if( childrenvals != NULL )
   {
      /* temporarily overwrite the evalvalue in all children and expr with values from childrenvals and exprval, resp. */
      if( expr->nchildren > 0 )
      {
         assert(bufmem != NULL);
         SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &origchildrenvals, expr->nchildren) );

         for( c = 0; c < expr->nchildren; ++c )
         {
            origchildrenvals[c] = expr->children[c]->evalvalue;
            expr->children[c]->evalvalue = childrenvals[c];
         }
      }

      origexprval = expr->evalvalue;
      expr->evalvalue = exprval;
   }

   SCIP_CALL( expr->exprhdlr->bwdiff(set->scip, expr, childidx, derivative) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*derivative) )
      *derivative = SCIP_INVALID;

   /* restore original evalvalues in children */
   if( childrenvals != NULL )
   {
      if( expr->nchildren > 0 )
      {
         for( c = 0; c < expr->nchildren; ++c )
            expr->children[c]->evalvalue = origchildrenvals[c];  /*lint !e644*/

         BMSfreeBufferMemoryArray(bufmem, &origchildrenvals);
      }

      expr->evalvalue = origexprval;   /*lint !e644*/
   }

   return SCIP_OKAY;
}

/** calls the forward differentiation callback of an expression handler
 *
 * @see SCIP_DECL_EXPRFWDIFF
 */
SCIP_RETCODE SCIPexprhdlrFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_Real*            dot,                /**< buffer to store derivative value */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions) */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(dot != NULL);

   if( exprhdlr->fwdiff == NULL )
   {
      *dot = SCIP_INVALID;
      return SCIP_OKAY;
   }

   SCIP_CALL( exprhdlr->fwdiff(set->scip, expr, dot, direction) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*dot) )
      *dot = SCIP_INVALID;

   return SCIP_OKAY;
}

/** calls the evaluation and forward-differentiation callback of an expression handler
 *
 * The method evaluates an expression by taking the values of its children into account.
 * The method differentiates an expression by taking the values and directional derivatives of its children into account.
 *
 * Further, allows to evaluate and differentiate w.r.t. given values for children instead of those stored in children expressions.
 *
 * It probably doesn't make sense to call this function for a variable-expression if sol and/or direction are not given.
 *
 * @see SCIP_DECL_EXPREVAL
 * @see SCIP_DECL_EXPRFWDIFF
 */
SCIP_RETCODE SCIPexprhdlrEvalFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory, can be NULL if childrenvals is NULL */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_Real*            val,                /**< buffer to store value of expression */
   SCIP_Real*            dot,                /**< buffer to store derivative value */
   SCIP_Real*            childrenvals,       /**< values for children, or NULL if values stored in children should
                                                  be used */
   SCIP_SOL*             sol,                /**< solution that is evaluated (can be NULL) */
   SCIP_Real*            childrendirs,       /**< directional derivatives for children, or NULL if dot-values stored
                                                  in children should be used */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions, can
                                                  be NULL if childrendirs is given) */
   )
{
   SCIP_Real origval;
   SCIP_Real* origvals = NULL;
   SCIP_Real* origdots = NULL;

   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(exprhdlr->eval != NULL);
   assert(val != NULL);
   assert(dot != NULL);

   /* temporarily overwrite the evalvalue in all children with values from childrenvals */
   if( childrenvals != NULL && expr->nchildren > 0 )
   {
      int c;

      assert(bufmem != NULL);

      SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &origvals, expr->nchildren) );

      for( c = 0; c < expr->nchildren; ++c )
      {
         origvals[c] = expr->children[c]->evalvalue;
         expr->children[c]->evalvalue = childrenvals[c];
      }
   }

   /* temporarily overwrite the dot in all children with values from childrendirs */
   if( childrendirs != NULL && expr->nchildren > 0 )
   {
      int c;

      assert(bufmem != NULL);

      SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &origdots, expr->nchildren) );

      for( c = 0; c < expr->nchildren; ++c )
      {
         origdots[c] = expr->children[c]->dot;
         expr->children[c]->dot = childrendirs[c];
      }
   }

   /* remember original value */
   origval = expr->evalvalue;

   /* call expression eval callback */
   SCIP_CALL( exprhdlr->eval(set->scip, expr, val, sol) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*val) )
      *val = SCIP_INVALID;

   /* temporarily overwrite evalvalue of expr, since some exprhdlr (e.g., product) access this value in fwdiff */
   expr->evalvalue = *val;

   /* call forward-differentiation callback (if available) */
   SCIP_CALL( SCIPexprhdlrFwDiffExpr(exprhdlr, set, expr, dot, direction) );

   /* restore original value */
   expr->evalvalue = origval;

   /* restore original dots in children */
   if( origdots != NULL )
   {
      int c;
      for( c = 0; c < expr->nchildren; ++c )
         expr->children[c]->dot = origdots[c];

      BMSfreeBufferMemoryArray(bufmem, &origdots);
   }

   /* restore original evalvalues in children */
   if( origvals != NULL )
   {
      int c;
      for( c = 0; c < expr->nchildren; ++c )
         expr->children[c]->evalvalue = origvals[c];

      BMSfreeBufferMemoryArray(bufmem, &origvals);
   }

   return SCIP_OKAY;
}

/** calls the evaluation callback for Hessian directions (backward over forward) of an expression handler
 *
 * @see SCIP_DECL_EXPRBWFWDIFF
 */
SCIP_RETCODE SCIPexprhdlrBwFwDiffExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   int                   childidx,           /**< index of the child */
   SCIP_Real*            bardot,             /**< buffer to store derivative value */
   SCIP_SOL*             direction           /**< direction of the derivative (useful only for var expressions) */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(childidx >= 0);
   assert(childidx < expr->nchildren);
   assert(bardot != NULL);

   if( exprhdlr->bwfwdiff == NULL )
   {
      *bardot = SCIP_INVALID;
      return SCIP_OKAY;
   }

   SCIP_CALL( expr->exprhdlr->bwfwdiff(set->scip, expr, childidx, bardot, direction) );

   /* if there was some evaluation error (e.g., overflow) that hasn't been caught yet, then do so now */
   if( !SCIPisFinite(*bardot) )
      *bardot = SCIP_INVALID;

   return SCIP_OKAY;
}

/** calls the interval evaluation callback of an expression handler
 *
 * @see SCIP_DECL_EXPRINTEVAL
 */
SCIP_RETCODE SCIPexprhdlrIntEvalExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_INTERVAL*        interval,           /**< buffer where to store interval */
   SCIP_DECL_EXPR_INTEVALVAR((*intevalvar)), /**< callback to be called when interval-evaluating a variable */
   void*                 intevalvardata      /**< data to be passed to intevalvar callback */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(interval != NULL);

   if( exprhdlr->inteval != NULL )
   {
      SCIPclockStart(exprhdlr->intevaltime, set);
      SCIP_CALL( exprhdlr->inteval(set->scip, expr, interval, intevalvar, intevalvardata) );
      SCIPclockStop(exprhdlr->intevaltime, set);

      ++exprhdlr->nintevalcalls;
   }

   return SCIP_OKAY;
}

/** calls the estimator callback of an expression handler
 *
 * @see SCIP_DECL_EXPRESTIMATE
 */
SCIP_RETCODE SCIPexprhdlrEstimateExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be estimated */
   SCIP_INTERVAL*        localbounds,        /**< current bounds for children */
   SCIP_INTERVAL*        globalbounds,       /**< global bounds for children */
   SCIP_Real*            refpoint,           /**< children values for the reference point where to estimate */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_Real             targetvalue,        /**< a value that the estimator shall exceed, can be +/-infinity */
   SCIP_Real*            coefs,              /**< array to store coefficients of estimator */
   SCIP_Real*            constant,           /**< buffer to store constant part of estimator */
   SCIP_Bool*            islocal,            /**< buffer to store whether estimator is valid locally only */
   SCIP_Bool*            success,            /**< buffer to indicate whether an estimator could be computed */
   SCIP_Bool*            branchcand          /**< array to indicate which children (not) to consider for branching */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(coefs != NULL);
   assert(islocal != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( exprhdlr->estimate != NULL )
   {
      SCIPclockStart(exprhdlr->estimatetime, set);
      SCIP_CALL( exprhdlr->estimate(set->scip, expr, localbounds, globalbounds, refpoint, overestimate, targetvalue,
            coefs, constant, islocal, success, branchcand) );
      SCIPclockStop(exprhdlr->estimatetime, set);

      /* update statistics */
      ++exprhdlr->nestimatecalls;
   }

   return SCIP_OKAY;
}

/** calls the intitial estimators callback of an expression handler
 *
 * @see SCIP_DECL_EXPRINITESTIMATES
 */
SCIP_RETCODE SCIPexprhdlrInitEstimatesExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to be estimated */
   SCIP_INTERVAL*        bounds,             /**< bounds for children */
   SCIP_Bool             overestimate,       /**< whether the expression shall be overestimated or underestimated */
   SCIP_Real*            coefs[SCIP_EXPR_MAXINITESTIMATES], /**< buffer to store coefficients of computed estimators */
   SCIP_Real             constant[SCIP_EXPR_MAXINITESTIMATES], /**< buffer to store constant of computed estimators */
   int*                  nreturned           /**< buffer to store number of estimators that have been computed */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(nreturned != NULL);

   *nreturned = 0;

   if( exprhdlr->initestimates )
   {
      SCIPclockStart(expr->exprhdlr->estimatetime, set);
      SCIP_CALL( exprhdlr->initestimates(set->scip, expr, bounds, overestimate, coefs, constant, nreturned) );
      SCIPclockStop(expr->exprhdlr->estimatetime, set);

      ++exprhdlr->nestimatecalls;
   }

   return SCIP_OKAY;
}

/** calls the simplification callback of an expression handler
 *
 * @see SCIP_DECL_EXPRSIMPLIFY
 */
SCIP_RETCODE SCIPexprhdlrSimplifyExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to simplify */
   SCIP_EXPR**           simplifiedexpr,     /**< buffer to store the simplified expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(simplifiedexpr != NULL);

   if( exprhdlr->simplify != NULL )
   {
      SCIPclockStart(expr->exprhdlr->simplifytime, set);
      SCIP_CALL( exprhdlr->simplify(set->scip, expr, simplifiedexpr, ownercreate, ownercreatedata) );
      SCIPclockStop(expr->exprhdlr->simplifytime, set);

      /* update statistics */
      ++exprhdlr->nsimplifycalls;
      if( expr != *simplifiedexpr )
         ++exprhdlr->nsimplified;
   }
   else
   {
      *simplifiedexpr = expr;

      /* if an expression handler doesn't implement simplify, we assume that it is already simplified
       * we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created
       */
      SCIPexprCapture(expr);
   }

   return SCIP_OKAY;
}

/** calls the reverse propagation callback of an expression handler
 *
 * The method propagates given bounds over the children of an expression.
 *
 * @see SCIP_DECL_EXPRREVERSEPROP
 */
SCIP_RETCODE SCIPexprhdlrReversePropExpr(
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to propagate */
   SCIP_INTERVAL         bounds,             /**< the bounds on the expression that should be propagated */
   SCIP_INTERVAL*        childrenbounds,     /**< array to store computed bounds for children, initialized with
                                                  current activity */
   SCIP_Bool*            infeasible          /**< buffer to store whether a children bounds were propagated to
                                                  an empty interval */
   )
{
   assert(exprhdlr != NULL);
   assert(set != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr == exprhdlr);
   assert(childrenbounds != NULL || expr->nchildren == 0);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   if( exprhdlr->reverseprop != NULL )
   {
      SCIPclockStart(exprhdlr->proptime, set);
      SCIP_CALL( exprhdlr->reverseprop(set->scip, expr, bounds, childrenbounds, infeasible) );
      SCIPclockStop(exprhdlr->proptime, set);

      /* update statistics */
      if( *infeasible )
         ++expr->exprhdlr->ncutoffs;
      ++expr->exprhdlr->npropcalls;
   }

   return SCIP_OKAY;
}

/**@name Expression Methods */
/**@{ */

/* from expr.h */

#ifdef NDEBUG
#undef SCIPexprCapture
#undef SCIPexprIsVar
#undef SCIPexprIsValue
#undef SCIPexprIsSum
#undef SCIPexprIsProduct
#undef SCIPexprIsPower
#endif

/** creates and captures an expression with given expression data and children */
SCIP_RETCODE SCIPexprCreate(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data (expression assumes ownership) */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children (can be NULL if nchildren is 0) */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   int c;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(children != NULL || nchildren == 0);
   assert(exprdata == NULL || exprhdlr->copydata != NULL); /* copydata must be available if there is expression data */
   assert(exprdata == NULL || exprhdlr->freedata != NULL); /* freedata must be available if there is expression data */

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, expr) );

   (*expr)->exprhdlr = exprhdlr;
   (*expr)->exprdata = exprdata;
   (*expr)->activitytag = -1;  /* to be less than initial domchgcount */
   (*expr)->curvature = SCIP_EXPRCURV_UNKNOWN;

   /* initialize activity to entire interval */
   SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &(*expr)->activity);

   if( nchildren > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*expr)->children, children, nchildren) );
      (*expr)->nchildren = nchildren;
      (*expr)->childrensize = nchildren;

      for( c = 0; c < nchildren; ++c )
         SCIPexprCapture((*expr)->children[c]);
   }

   SCIPexprCapture(*expr);

   ++exprhdlr->ncreated;

   /* initializes the ownerdata */
   if( ownercreate != NULL )
   {
      SCIP_CALL( ownercreate(set->scip, *expr, &(*expr)->ownerdata, &(*expr)->ownerfree, &(*expr)->ownerprint,
            &(*expr)->ownerevalactivity, ownercreatedata) );
   }

   return SCIP_OKAY;
}

/** appends child to the children list of expr */
SCIP_RETCODE SCIPexprAppendChild(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR*            child               /**< expression to be appended */
   )
{
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(child != NULL);
   assert(expr->nchildren <= expr->childrensize);

   if( expr->nchildren == expr->childrensize )
   {
      expr->childrensize = SCIPsetCalcMemGrowSize(set, expr->nchildren+1);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &expr->children, expr->nchildren, expr->childrensize) );
   }

   expr->children[expr->nchildren] = child;
   ++expr->nchildren;

   /* capture child */
   SCIPexprCapture(child);

   return SCIP_OKAY;
}

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured, unless they are the same (=same pointer)
 */
SCIP_RETCODE SCIPexprReplaceChild(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression where a child is going to be replaced */
   int                   childidx,           /**< index of child being replaced */
   SCIP_EXPR*            newchild            /**< the new child */
   )
{
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(newchild != NULL);
   assert(childidx >= 0);
   assert(childidx < expr->nchildren);

   /* do nothing if child is not changing */
   if( newchild == expr->children[childidx] )
      return SCIP_OKAY;

   /* capture new child (do this before releasing the old child in case there are equal */
   SCIPexprCapture(newchild);

   SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &(expr->children[childidx])) );
   expr->children[childidx] = newchild;

   return SCIP_OKAY;
}

/** remove all children of expr */
SCIP_RETCODE SCIPexprRemoveChildren(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int c;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);

   for( c = 0; c < expr->nchildren; ++c )
   {
      assert(expr->children[c] != NULL);
      SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &(expr->children[c])) );
   }

   expr->nchildren = 0;

   return SCIP_OKAY;
}

/** copies an expression including subexpressions
 *
 * @note If copying fails due to an expression handler not being available in the targetscip, then *targetexpr will be set to NULL.
 *
 * For all or some expressions, a mapping to an existing expression can be specified via the mapexpr callback.
 * The mapped expression (including its children) will not be copied in this case and its ownerdata will not be touched.
 * If, however, the mapexpr callback returns NULL for the targetexpr, then the expr will be copied in the usual way.
 */
SCIP_RETCODE SCIPexprCopy(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             targetset,          /**< global SCIP settings data structure where target expression will live */
   SCIP_STAT*            targetstat,         /**< dynamic problem statistics in target SCIP */
   BMS_BLKMEM*           targetblkmem,       /**< block memory in target SCIP */
   SCIP_EXPR*            sourceexpr,         /**< expression to be copied */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to copy of source expression */
   SCIP_DECL_EXPR_MAPEXPR((*mapexpr)),       /**< expression mapping function, or NULL for creating new expressions */
   void*                 mapexprdata,        /**< data of expression mapping function */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPRITER_USERDATA expriteruserdata;
   SCIP_EXPR* expr;
   SCIP* sourcescip = set->scip;        /* SCIP data structure corresponding to source expression */
   SCIP* targetscip = targetset->scip;  /* SCIP data structure where target expression will live */

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(targetset != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);
   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, sourceexpr, SCIP_EXPRITER_DFS, TRUE) );  /*TODO use FALSE, i.e., don't duplicate common subexpr? */
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_VISITEDCHILD);

   expr = sourceexpr;
   while( !SCIPexpriterIsEnd(it) )
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_ENTEREXPR :
         {
            /* create expr that will hold the copy */
            SCIP_EXPRHDLR* targetexprhdlr;
            SCIP_EXPRDATA* targetexprdata;
            SCIP_EXPR* exprcopy = NULL;

            if( mapexpr != NULL )
            {
               SCIP_CALL( mapexpr(targetscip, &exprcopy, sourcescip, expr, ownercreate, ownercreatedata, mapexprdata) );
               if( exprcopy != NULL )
               {
                  /* map callback gave us an expression to use for the copy */
                  /* store targetexpr */
                  expriteruserdata.ptrval = exprcopy;
                  SCIPexpriterSetCurrentUserData(it, expriteruserdata);

                  /* skip subexpression (assume that exprcopy is a complete copy) and continue */
                  expr = SCIPexpriterSkipDFS(it);
                  continue;
               }
            }

            /* get the exprhdlr of the target scip */
            if( targetscip != sourcescip )
            {
               targetexprhdlr = SCIPsetFindExprhdlr(targetset, expr->exprhdlr->name);

               if( targetexprhdlr == NULL )
               {
                  /* expression handler not in target scip (probably did not have a copy callback) -> abort */
                  expriteruserdata.ptrval = NULL;
                  SCIPexpriterSetCurrentUserData(it, expriteruserdata);

                  expr = SCIPexpriterSkipDFS(it);
                  continue;
               }
            }
            else
            {
               targetexprhdlr = expr->exprhdlr;
            }
            assert(targetexprhdlr != NULL);

            /* copy expression data */
            if( expr->exprdata != NULL )
            {
               assert(expr->exprhdlr->copydata != NULL);
               SCIP_CALL( expr->exprhdlr->copydata(targetscip, targetexprhdlr, &targetexprdata, sourcescip, expr) );
            }
            else
            {
               targetexprdata = NULL;
            }

            /* create in targetexpr an expression of the same type as expr, but without children for now */
            SCIP_CALL( SCIPexprCreate(targetset, targetblkmem, &exprcopy, targetexprhdlr, targetexprdata, 0, NULL,
                  ownercreate, ownercreatedata) );

            /* store targetexpr */
            expriteruserdata.ptrval = exprcopy;
            SCIPexpriterSetCurrentUserData(it, expriteruserdata);

            break;
         }

         case SCIP_EXPRITER_VISITEDCHILD :
         {
            /* just visited child so a copy of himself should be available; append it */
            SCIP_EXPR* exprcopy;
            SCIP_EXPR* childcopy;

            exprcopy = (SCIP_EXPR*)SCIPexpriterGetCurrentUserData(it).ptrval;

            /* get copy of child */
            childcopy = (SCIP_EXPR*)SCIPexpriterGetChildUserDataDFS(it).ptrval;
            if( childcopy == NULL )
            {
               /* abort */
               /* release exprcopy (should free also the already copied children) */
               SCIP_CALL( SCIPexprRelease(targetset, targetstat, targetblkmem, (SCIP_EXPR**)&exprcopy) );

               expriteruserdata.ptrval = NULL;
               SCIPexpriterSetCurrentUserData(it, expriteruserdata);

               expr = SCIPexpriterSkipDFS(it);
               continue;
            }

            /* append child to exprcopy */
            SCIP_CALL( SCIPexprAppendChild(targetset, targetblkmem, exprcopy, childcopy) );

            /* release childcopy (still captured by exprcopy) */
            SCIP_CALL( SCIPexprRelease(targetset, targetstat, targetblkmem, &childcopy) );

            break;
         }

         default:
            /* we should never be called in this stage */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

   /* the target expression should be stored in the userdata of the sourceexpr (can be NULL if aborted) */
   *targetexpr = (SCIP_EXPR*)SCIPexpriterGetExprUserData(it, sourceexpr).ptrval;

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** duplicates the given expression without its children */
SCIP_RETCODE SCIPexprDuplicateShallow(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store (shallow) duplicate of expr */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call on expression copy to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRDATA* exprdatacopy = NULL;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(copyexpr != NULL);

   /* copy expression data */
   if( expr->exprdata != NULL )
   {
      assert(expr->exprhdlr->copydata != NULL);
      SCIP_CALL( expr->exprhdlr->copydata(set->scip, expr->exprhdlr, &exprdatacopy, set->scip, expr) );
   }

   /* create expression with same handler and copied data, but without children */
   SCIP_CALL( SCIPexprCreate(set, blkmem, copyexpr, expr->exprhdlr, exprdatacopy, 0, NULL, ownercreate,
         ownercreatedata) );

   return SCIP_OKAY;
}

/** captures an expression (increments usage count) */
void SCIPexprCapture(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   ++expr->nuses;
}

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_RETCODE SCIPexprRelease(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           rootexpr            /**< pointer to expression */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;

   assert(rootexpr != NULL);
   assert(*rootexpr != NULL);
   assert((*rootexpr)->nuses > 0);

   if( (*rootexpr)->nuses > 1 )
   {
      --(*rootexpr)->nuses;
      *rootexpr = NULL;

      return SCIP_OKAY;
   }

   /* handle the root expr separately: free ownerdata, quaddata, and exprdata first */

   /* call ownerfree callback, if given
    * we intentially call this also if ownerdata is NULL, so owner can be notified without storing data
    */
   if( (*rootexpr)->ownerfree != NULL )
   {
      SCIP_CALL( (*rootexpr)->ownerfree(set->scip, *rootexpr, &(*rootexpr)->ownerdata) );
      assert((*rootexpr)->ownerdata == NULL);
   }

   /* free quadratic info */
   SCIPexprFreeQuadratic(blkmem, *rootexpr);

   /* free expression data */
   if( (*rootexpr)->exprdata != NULL )
   {
      assert((*rootexpr)->exprhdlr->freedata != NULL);
      SCIP_CALL( (*rootexpr)->exprhdlr->freedata(set->scip, *rootexpr) );
   }

   /* now release and free children, where no longer in use */
   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, *rootexpr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_VISITEDCHILD);
   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it) ; )
   {
      /* expression should be used by its parent and maybe by the iterator (only the root!)
       * in VISITEDCHILD we assert that expression is only used by its parent
       */
      assert(expr != NULL);
      assert(0 <= expr->nuses && expr->nuses <= 2);

      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITINGCHILD :
         {
            /* check whether a child needs to be visited (nuses == 1)
             * if not, then we still have to release it
             */
            SCIP_EXPR* child;

            child = SCIPexpriterGetChildExprDFS(it);
            if( child->nuses > 1 )
            {
               /* child is not going to be freed: just release it */
               SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &child) );
               expr = SCIPexpriterSkipDFS(it);
               continue;
            }

            assert(child->nuses == 1);

            /* free child's quaddata, ownerdata, and exprdata when entering child */
            if( child->ownerfree != NULL )
            {
               SCIP_CALL( child->ownerfree(set->scip, child, &child->ownerdata) );
               assert(child->ownerdata == NULL);
            }

            /* free quadratic info */
            SCIPexprFreeQuadratic(blkmem, child);

            /* free expression data */
            if( child->exprdata != NULL )
            {
               assert(child->exprhdlr->freedata != NULL);
               SCIP_CALL( child->exprhdlr->freedata(set->scip, child) );
               assert(child->exprdata == NULL);
            }

            break;
         }

         case SCIP_EXPRITER_VISITEDCHILD :
         {
            /* free child after visiting it */
            SCIP_EXPR* child;

            child = SCIPexpriterGetChildExprDFS(it);
            /* child should only be used by its parent */
            assert(child->nuses == 1);

            /* child should have no data associated */
            assert(child->exprdata == NULL);

            /* free child expression */
            SCIP_CALL( freeExpr(blkmem, &child) );
            expr->children[SCIPexpriterGetChildIdxDFS(it)] = NULL;

            break;
         }

         default:
            SCIPABORT(); /* we should never be called in this stage */
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

   SCIPexpriterFree(&it);

   /* handle the root expr separately: free its children and itself here */
   SCIP_CALL( freeExpr(blkmem, rootexpr) );

   return SCIP_OKAY;
}

/** returns whether an expression is a variable expression */
SCIP_Bool SCIPexprIsVar(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(set != NULL);
   assert(expr != NULL);

   return expr->exprhdlr == set->exprhdlrvar;
}

/** returns whether an expression is a value expression */
SCIP_Bool SCIPexprIsValue(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(set != NULL);
   assert(expr != NULL);

   return expr->exprhdlr == set->exprhdlrval;
}

/** returns whether an expression is a sum expression */
SCIP_Bool SCIPexprIsSum(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(set != NULL);
   assert(expr != NULL);

   return expr->exprhdlr == set->exprhdlrsum;
}

/** returns whether an expression is a product expression */
SCIP_Bool SCIPexprIsProduct(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(set != NULL);
   assert(expr != NULL);

   return expr->exprhdlr == set->exprhdlrproduct;
}

/** returns whether an expression is a power expression */
SCIP_Bool SCIPexprIsPower(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(set != NULL);
   assert(expr != NULL);

   return expr->exprhdlr == set->exprhdlrpow;
}

/** print an expression as info-message */
SCIP_RETCODE SCIPexprPrint(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPR*            expr                /**< expression to be printed */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPRITER_STAGE stage;
   int currentchild;
   unsigned int parentprecedence;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ALLSTAGES);

   while( !SCIPexpriterIsEnd(it) )
   {
      assert(expr->exprhdlr != NULL);
      stage = SCIPexpriterGetStageDFS(it);

      if( stage == SCIP_EXPRITER_VISITEDCHILD || stage == SCIP_EXPRITER_VISITINGCHILD )
         currentchild = SCIPexpriterGetChildIdxDFS(it);
      else
         currentchild = -1;

      if( SCIPexpriterGetParentDFS(it) != NULL )
         parentprecedence = SCIPexprhdlrGetPrecedence(SCIPexprGetHdlr(SCIPexpriterGetParentDFS(it)));
      else
         parentprecedence = 0;

      SCIP_CALL( SCIPexprhdlrPrintExpr(expr->exprhdlr, set, messagehdlr, expr, stage, currentchild,
            parentprecedence, file) );

      expr = SCIPexpriterGetNext(it);
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format to a give FILE* pointer */
SCIP_RETCODE SCIPexprPrintDotInit(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRPRINTDATA**  printdata,          /**< buffer to store dot printing data */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPRPRINT_WHAT   whattoprint         /**< info on what to print for each expression */
   )
{
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(printdata != NULL);

   if( file == NULL )
      file = stdout;

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, printdata) );

   (*printdata)->file = file;
   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &(*printdata)->iterator) );
   (*printdata)->closefile = FALSE;
   (*printdata)->whattoprint = whattoprint;
   SCIP_CALL( SCIPhashmapCreate(&(*printdata)->leaveexprs, blkmem, 100) );

   fputs("strict digraph exprgraph {\n", file);
   fputs("node [fontcolor=white, style=filled, rankdir=LR]\n", file);

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_RETCODE SCIPexprPrintDotInit2(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRPRINTDATA**  printdata,          /**< buffer to store dot printing data */
   const char*           filename,           /**< name of file to print to */
   SCIP_EXPRPRINT_WHAT   whattoprint         /**< info on what to print for each expression */
   )
{
   FILE* f;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(printdata != NULL);
   assert(filename != NULL);

   f = fopen(filename, "w");
   if( f == NULL )
   {
      SCIPerrorMessage("could not open file <%s> for writing\n", filename);  /* error code would be in errno */
      return SCIP_FILECREATEERROR;
   }

   SCIP_CALL_FINALLY( SCIPexprPrintDotInit(set, stat, blkmem, printdata, f, whattoprint),
      fclose(f) );
   (*printdata)->closefile = TRUE;

   return SCIP_OKAY;
} /*lint !e429*/

/** main part of printing an expression in dot format */
SCIP_RETCODE SCIPexprPrintDot(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPRPRINTDATA*   printdata,          /**< data as initialized by \ref SCIPprintExprDotInit() */
   SCIP_EXPR*            expr                /**< expression to be printed */
   )
{
   SCIP_Real color;
   int c;

   assert(set != NULL);
   assert(printdata != NULL);
   assert(expr != NULL);
   assert(expr->exprhdlr != NULL);

   SCIP_CALL( SCIPexpriterInit(printdata->iterator, expr, SCIP_EXPRITER_DFS, FALSE) );

   while( !SCIPexpriterIsEnd(printdata->iterator) )
   {
      /* print expression as dot node */

      if( expr->nchildren == 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(printdata->leaveexprs, (void*)expr, NULL) );
      }

      /* make up some color from the expression type (it's name) */
      color = 0.0;
      for( c = 0; expr->exprhdlr->name[c] != '\0'; ++c )
         color += (tolower(expr->exprhdlr->name[c]) - 'a') / 26.0;
      color = SCIPsetFrac(set, color);
      fprintf(printdata->file, "n%p [fillcolor=\"%g,%g,%g\", label=\"", (void*)expr, color, color, color);

      if( printdata->whattoprint & SCIP_EXPRPRINT_EXPRHDLR )
      {
         fprintf(printdata->file, "%s\\n", expr->exprhdlr->name);
      }

      if( printdata->whattoprint & SCIP_EXPRPRINT_EXPRSTRING )
      {
         SCIP_CALL( SCIPexprhdlrPrintExpr(expr->exprhdlr, set, messagehdlr, expr, SCIP_EXPRITER_ENTEREXPR, -1, 0,
               printdata->file) );
         for( c = 0; c < expr->nchildren; ++c )
         {
            SCIP_CALL( SCIPexprhdlrPrintExpr(expr->exprhdlr, set, messagehdlr, expr, SCIP_EXPRITER_VISITINGCHILD,
                  c, 0, printdata->file) );
            fprintf(printdata->file, "c%d", c);
            SCIP_CALL( SCIPexprhdlrPrintExpr(expr->exprhdlr, set, messagehdlr, expr, SCIP_EXPRITER_VISITEDCHILD,
                  c, 0, printdata->file) );
         }
         SCIP_CALL( SCIPexprhdlrPrintExpr(expr->exprhdlr, set, messagehdlr, expr, SCIP_EXPRITER_LEAVEEXPR, -1, 0,
               printdata->file) );

         fputs("\\n", printdata->file);
      }

      if( printdata->whattoprint & SCIP_EXPRPRINT_NUSES )
      {
         /* print number of uses */
         fprintf(printdata->file, "%d uses\\n", expr->nuses);
      }

      if( printdata->whattoprint & SCIP_EXPRPRINT_OWNER )
      {
         /* print ownerdata */
         if( expr->ownerprint != NULL )
         {
            SCIP_CALL( expr->ownerprint(set->scip, printdata->file, expr, expr->ownerdata) );
         }
         else if( expr->ownerdata != NULL )
         {
            fprintf(printdata->file, "owner=%p\\n", (void*)expr->ownerdata);
         }
      }

      if( printdata->whattoprint & SCIP_EXPRPRINT_EVALVALUE )
      {
         /* print eval value */
         fprintf(printdata->file, "val=%g", expr->evalvalue);

         if( (printdata->whattoprint & SCIP_EXPRPRINT_EVALTAG) == SCIP_EXPRPRINT_EVALTAG )
         {
            /* print also eval tag */
            fprintf(printdata->file, " (%" SCIP_LONGINT_FORMAT ")", expr->evaltag);
         }
         fputs("\\n", printdata->file);
      }

      if( printdata->whattoprint & SCIP_EXPRPRINT_ACTIVITY )
      {
         /* print activity */
         fprintf(printdata->file, "[%g,%g]", expr->activity.inf, expr->activity.sup);

         if( (printdata->whattoprint & SCIP_EXPRPRINT_ACTIVITYTAG) == SCIP_EXPRPRINT_ACTIVITYTAG )
         {
            /* print also activity eval tag */
            fprintf(printdata->file, " (%" SCIP_LONGINT_FORMAT ")", expr->activitytag);
         }
         fputs("\\n", printdata->file);
      }

      fputs("\"]\n", printdata->file);  /* end of label and end of node */

      /* add edges from expr to its children */
      for( c = 0; c < expr->nchildren; ++c )
         fprintf(printdata->file, "n%p -> n%p [label=\"c%d\"]\n", (void*)expr, (void*)expr->children[c], c);

      expr = SCIPexpriterGetNext(printdata->iterator);
   }

   return SCIP_OKAY;
}

/** finishes printing of expressions in dot format */
SCIP_RETCODE SCIPexprPrintDotFinal(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRPRINTDATA**  printdata           /**< buffer where dot printing data has been stored */
   )
{
   SCIP_EXPR* expr;
   SCIP_HASHMAPENTRY* entry;
   FILE* file;
   int i;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(printdata != NULL);
   assert(*printdata != NULL);

   file = (*printdata)->file;
   assert(file != NULL);

   /* iterate through all entries of the map */
   fputs("{rank=same;", file);
   for( i = 0; i < SCIPhashmapGetNEntries((*printdata)->leaveexprs); ++i )
   {
      entry = SCIPhashmapGetEntry((*printdata)->leaveexprs, i);

      if( entry != NULL )
      {
         expr = (SCIP_EXPR*) SCIPhashmapEntryGetOrigin(entry);
         assert(expr != NULL);
         assert(expr->nchildren == 0);

         fprintf(file, " n%p", (void*)expr);
      }
   }
   fprintf(file, "}\n");

   fprintf(file, "}\n");

   SCIPhashmapFree(&(*printdata)->leaveexprs);

   SCIPexpriterFree(&(*printdata)->iterator);

   if( (*printdata)->closefile )
      fclose((*printdata)->file);

   BMSfreeBlockMemory(blkmem, printdata);

   return SCIP_OKAY;
}

/** prints structure of an expression a la Maple's dismantle */
SCIP_RETCODE SCIPexprDismantle(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPR*            expr                /**< expression to dismantle */
   )
{
   SCIP_EXPRITER* it;
   int depth = -1;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(messagehdlr != NULL);
   assert(expr != NULL);

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);

   for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_ENTEREXPR:
         {
            int nspaces;

            ++depth;
            nspaces = 3 * depth;

            /* use depth of expression to align output */
            SCIPmessageFPrintInfo(messagehdlr, file, "%*s[%s]: ", nspaces, "", expr->exprhdlr->name);

            if( SCIPexprIsVar(set, expr) )
            {
               SCIP_VAR* var;

               var = SCIPgetVarExprVar(expr);
               SCIPmessageFPrintInfo(messagehdlr, file, "%s in [%g, %g]", SCIPvarGetName(var), SCIPvarGetLbLocal(var),
                  SCIPvarGetUbLocal(var));
            }
            else if( SCIPexprIsSum(set, expr) )
               SCIPmessageFPrintInfo(messagehdlr, file, "%g", SCIPgetConstantExprSum(expr));
            else if( SCIPexprIsProduct(set, expr) )
               SCIPmessageFPrintInfo(messagehdlr, file, "%g", SCIPgetCoefExprProduct(expr));
            else if( SCIPexprIsValue(set, expr) )
               SCIPmessageFPrintInfo(messagehdlr, file, "%g", SCIPgetValueExprValue(expr));
            else if( SCIPexprIsPower(set, expr) || strcmp(expr->exprhdlr->name, "signpower") == 0)
               SCIPmessageFPrintInfo(messagehdlr, file, "%g", SCIPgetExponentExprPow(expr));

            SCIPmessageFPrintInfo(messagehdlr, file, "\n");

            if( expr->ownerprint != NULL )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, "%*s   ", nspaces, "");
               SCIP_CALL( expr->ownerprint(set->scip, file, expr, expr->ownerdata) );
            }

            break;
         }

         case SCIP_EXPRITER_VISITINGCHILD:
         {
            int nspaces = 3 * depth;

            if( SCIPexprIsSum(set, expr) )
            {
               SCIPmessageFPrintInfo(messagehdlr, file, "%*s   ", nspaces, "");
               SCIPmessageFPrintInfo(messagehdlr, file, "[coef]: %g\n", SCIPgetCoefsExprSum(expr)[SCIPexpriterGetChildIdxDFS(it)]);
            }

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR:
         {
            --depth;
            break;
         }

         default:
            /* shouldn't be here */
            SCIPABORT();
            break;
      }
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** evaluate an expression in a point
 *
 * Iterates over expressions to also evaluate children, if necessary.
 * Value can be received via SCIPexprGetEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPexprGetEvalTag().
 */
SCIP_RETCODE SCIPexprEval(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_EXPRITER* it;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);

   /* if value is up-to-date, then nothing to do */
   if( soltag != 0 && expr->evaltag == soltag )
      return SCIP_OKAY;

   /* assume we'll get a domain error, so we don't have to get this expr back if we abort the iteration
    * if there is no domain error, then we will overwrite the evalvalue in the last leaveexpr stage
    */
   expr->evalvalue = SCIP_INVALID;
   expr->evaltag = soltag;

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);

   while( !SCIPexpriterIsEnd(it) )
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITINGCHILD :
         {
            SCIP_EXPR* child;

            if( soltag == 0 )
               break;

            /* check whether child has been evaluated for that solution already */
            child = SCIPexpriterGetChildExprDFS(it);
            if( soltag == child->evaltag )
            {
               if( child->evalvalue == SCIP_INVALID )
                  goto TERMINATE;

               /* skip this child
                * this already returns the next one, so continue with loop
                */
               expr = SCIPexpriterSkipDFS(it);
               continue;
            }

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR :
         {
            SCIP_CALL( SCIPexprhdlrEvalExpr(expr->exprhdlr, set, NULL , expr, &expr->evalvalue, NULL, sol) );
            expr->evaltag = soltag;

            if( expr->evalvalue == SCIP_INVALID )
               goto TERMINATE;

            break;
         }

         default :
            /* we should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

TERMINATE:
   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** evaluates gradient of an expression for a given point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_RETCODE SCIPexprEvalGradient(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr,           /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_EXPR* child;
   SCIP_Real derivative;
   SCIP_Longint difftag;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(rootexpr != NULL);

   /* ensure expression is evaluated */
   SCIP_CALL( SCIPexprEval(set, stat, blkmem, rootexpr, sol, soltag) );

   /* check if expression could not be evaluated */
   if( SCIPexprGetEvalValue(rootexpr) == SCIP_INVALID )
   {
      rootexpr->derivative = SCIP_INVALID;
      return SCIP_OKAY;
   }

   if( SCIPexprIsValue(set, rootexpr) )
   {
      rootexpr->derivative = 0.0;
      return SCIP_OKAY;
   }

   difftag = ++(stat->exprlastdifftag);

   rootexpr->derivative = 1.0;
   rootexpr->difftag = difftag;

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      assert(expr->evalvalue != SCIP_INVALID);

      child = SCIPexpriterGetChildExprDFS(it);
      assert(child != NULL);

      /* reset the value of the partial derivative w.r.t. a variable expression if we see it for the first time */
      if( child->difftag != difftag && SCIPexprIsVar(set, child) )
         child->derivative = 0.0;

      /* update differentiation tag of the child */
      child->difftag = difftag;

      /* call backward differentiation callback */
      if( SCIPexprIsValue(set, child) )
      {
         derivative = 0.0;
      }
      else
      {
         derivative = SCIP_INVALID;
         SCIP_CALL( SCIPexprhdlrBwDiffExpr(expr->exprhdlr, set, NULL, expr, SCIPexpriterGetChildIdxDFS(it),
               &derivative, NULL, 0.0) );

         if( derivative == SCIP_INVALID )
         {
            rootexpr->derivative = SCIP_INVALID;
            break;
         }
      }

      /* update partial derivative stored in the child expression
       * for a variable, we have to sum up the partial derivatives of the root w.r.t. this variable over all parents
       * for other intermediate expressions, we only store the partial derivative of the root w.r.t. this expression
       */
      if( !SCIPexprIsVar(set, child) )
         child->derivative = expr->derivative * derivative;
      else
         child->derivative += expr->derivative * derivative;
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** evaluates Hessian-vector product of an expression for a given point and direction
 *
 * Evaluates children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffGradientDirNonlinear()
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_RETCODE SCIPexprEvalHessianDir(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr,           /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   SCIP_Longint          soltag,             /**< tag that uniquely identifies the solution (with its values), or 0. */
   SCIP_SOL*             direction           /**< direction */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_EXPR* child;
   SCIP_Real derivative;
   SCIP_Real hessiandir;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(rootexpr != NULL);

   if( SCIPexprIsValue(set, rootexpr) )
   {
      rootexpr->dot = 0.0;
      rootexpr->bardot = 0.0;
      return SCIP_OKAY;
   }

   /* evaluate expression and directional derivative */
   SCIP_CALL( evalAndDiff(set, stat, blkmem, rootexpr, sol, soltag, direction) );

   if( rootexpr->evalvalue == SCIP_INVALID )
   {
      rootexpr->derivative = SCIP_INVALID;
      rootexpr->bardot = SCIP_INVALID;
      return SCIP_OKAY;
   }

   rootexpr->derivative = 1.0;

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);

   /* compute reverse diff and bardots: i.e. hessian times direction */
   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      assert(expr->evalvalue != SCIP_INVALID);

      child = SCIPexpriterGetChildExprDFS(it);
      assert(child != NULL);

      /* call backward and forward-backward differentiation callback */
      if( SCIPexprIsValue(set, child) )
      {
         derivative = 0.0;
         hessiandir = 0.0;
      }
      else
      {
         derivative = SCIP_INVALID;
         hessiandir = SCIP_INVALID;
         SCIP_CALL( SCIPexprhdlrBwDiffExpr(expr->exprhdlr, set, NULL, expr, SCIPexpriterGetChildIdxDFS(it),
               &derivative, NULL, SCIP_INVALID) );
         SCIP_CALL( SCIPexprhdlrBwFwDiffExpr(expr->exprhdlr, set, expr, SCIPexpriterGetChildIdxDFS(it),
               &hessiandir, NULL) );

         if( derivative == SCIP_INVALID || hessiandir == SCIP_INVALID )
         {
            rootexpr->derivative = SCIP_INVALID;
            rootexpr->bardot = SCIP_INVALID;
            break;
         }
      }

      /* update partial derivative and hessian stored in the child expression
       * for a variable, we have to sum up the partial derivatives of the root w.r.t. this variable over all parents
       * for other intermediate expressions, we only store the partial derivative of the root w.r.t. this expression
       */
      if( !SCIPexprIsVar(set, child) )
      {
         child->derivative = expr->derivative * derivative;
         child->bardot = expr->bardot * derivative + expr->derivative * hessiandir;
      }
      else
      {
         child->derivative += expr->derivative * derivative;
         child->bardot += expr->bardot * derivative + expr->derivative * hessiandir;
      }
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is no longer uptodate.
 * If the expr owner provided a evalactivity-callback, then call this.
 * Otherwise, loop over descendants and compare activitytag with stat's domchgcount, i.e.,
 * whether some bound was changed since last evaluation, to check whether exprhdlrs INTEVAL should be called.
 *
 * @note If expression is set to be integral, then activities are tightened to integral values.
 *   Thus, ensure that the integrality information is valid (if set to TRUE; the default (FALSE) is always ok).
 */
SCIP_RETCODE SCIPexprEvalActivity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr            /**< expression */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;

   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(rootexpr != NULL);

   if( rootexpr->ownerevalactivity != NULL )
   {
      /* call owner callback for activity-eval */
      SCIP_CALL( rootexpr->ownerevalactivity(set->scip, rootexpr, rootexpr->ownerdata) );

      return SCIP_OKAY;
   }

   /* fallback if no callback is given */

   assert(rootexpr->activitytag <= stat->domchgcount);

   /* if value is up-to-date, then nothing to do */
   if( rootexpr->activitytag == stat->domchgcount )
   {
#ifdef DEBUG_PROP
      SCIPsetDebugMsg(set, "activitytag of root expr equals domchgcount (%u), skip evalactivity\n", stat->domchgcount);
#endif

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it);  )
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITINGCHILD :
         {
            /* skip child if it has been evaluated already */
            SCIP_EXPR* child;

            child = SCIPexpriterGetChildExprDFS(it);
            if( child->activitytag == stat->domchgcount )
            {
               expr = SCIPexpriterSkipDFS(it);
               continue;
            }

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR :
         {
            /* we should not have entered this expression if its activity was already uptodate */
            assert(expr->activitytag < stat->domchgcount);

            /* reset activity to entire if invalid, so we can use it as starting point below */
            SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &expr->activity);

#ifdef DEBUG_PROP
            SCIPsetDebugMsg(set, "interval evaluation of expr %p ", (void*)expr);
            SCIP_CALL( SCIPprintExpr(set->scip, expr, NULL) );
            SCIPsetDebugMsgPrint(set, "\n");
#endif

            /* call the inteval callback of the exprhdlr */
            SCIP_CALL( SCIPexprhdlrIntEvalExpr(expr->exprhdlr, set, expr, &expr->activity, NULL, NULL) );
#ifdef DEBUG_PROP
            SCIPsetDebugMsg(set, " exprhdlr <%s>::inteval = [%.20g, %.20g]", expr->exprhdlr->name, expr->activity.inf,
                  expr->activity.sup);
#endif

            /* if expression is integral, then we try to tighten the interval bounds a bit
             * this should undo the addition of some unnecessary safety added by use of nextafter() in interval
             * arithmetics, e.g., when doing pow() it would be ok to use ceil() and floor(), but for safety we
             * use SCIPceil and SCIPfloor for now the default intevalVar does not relax variables, so can omit
             * expressions without children (constants should be ok, too)
             */
            if( expr->isintegral && expr->nchildren > 0 )
            {
               if( expr->activity.inf > -SCIP_INTERVAL_INFINITY )
                  expr->activity.inf = SCIPsetCeil(set, expr->activity.inf);
               if( expr->activity.sup <  SCIP_INTERVAL_INFINITY )
                  expr->activity.sup = SCIPsetFloor(set, expr->activity.sup);
#ifdef DEBUG_PROP
               SCIPsetDebugMsg(set, " applying integrality: [%.20g, %.20g]\n", expr->activity.inf, expr->activity.sup);
#endif
            }

            /* mark activity as empty if either the lower/upper bound is above/below +/- SCIPinfinity()
             * TODO this is a problem if dual-presolve fixed a variable to +/- infinity
             */
            if( SCIPsetIsInfinity(set, expr->activity.inf) || SCIPsetIsInfinity(set, -expr->activity.sup) )
            {
               SCIPsetDebugMsg(set, "treat activity [%g,%g] as empty as beyond infinity\n", expr->activity.inf, expr->activity.sup);
               SCIPintervalSetEmpty(&expr->activity);
            }

            /* remember that activity is uptodate now */
            expr->activitytag = stat->domchgcount;

            break;
         }

         default:
            /* you should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** compare expressions
 *
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note The given expressions are assumed to be simplified.
 */
int SCIPexprCompare(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr1,              /**< first expression */
   SCIP_EXPR*            expr2               /**< second expression */
   )
{
   SCIP_EXPRHDLR* exprhdlr1;
   SCIP_EXPRHDLR* exprhdlr2;
   int retval;

   exprhdlr1 = expr1->exprhdlr;
   exprhdlr2 = expr2->exprhdlr;

   /* expressions are of the same kind/type; use compare callback or default method */
   if( exprhdlr1 == exprhdlr2 )
      return SCIPexprhdlrCompareExpr(set, expr1, expr2);

   /* expressions are of different kind/type */
   /* enforces OR6 */
   if( SCIPexprIsValue(set, expr1) )
   {
      return -1;
   }
   /* enforces OR12 */
   if( SCIPexprIsValue(set, expr2) )
      return -SCIPexprCompare(set, expr2, expr1);

   /* enforces OR7 */
   if( SCIPexprIsSum(set, expr1) )
   {
      int compareresult;
      int nchildren;

      nchildren = expr1->nchildren;
      compareresult = SCIPexprCompare(set, expr1->children[nchildren-1], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* "base" of the largest expression of the sum is equal to expr2, coefficient might tell us that
       * expr2 is larger */
      if( SCIPgetCoefsExprSum(expr1)[nchildren-1] < 1.0 )
         return -1;

      /* largest expression of sum is larger or equal than expr2 => expr1 > expr2 */
      return 1;
   }
   /* enforces OR12 */
   if( SCIPexprIsSum(set, expr2) )
      return -SCIPexprCompare(set, expr2, expr1);

   /* enforces OR8 */
   if( SCIPexprIsProduct(set, expr1) )
   {
      int compareresult;
      int nchildren;

      nchildren = expr1->nchildren;
      compareresult = SCIPexprCompare(set, expr1->children[nchildren-1], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* largest expression of product is larger or equal than expr2 => expr1 > expr2 */
      return 1;
   }
   /* enforces OR12 */
   if( SCIPexprIsProduct(set, expr2) )
      return -SCIPexprCompare(set, expr2, expr1);

   /* enforces OR9 */
   if( SCIPexprIsPower(set, expr1) )
   {
      int compareresult;

      compareresult = SCIPexprCompare(set, expr1->children[0], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* base equal to expr2, exponent might tell us that expr2 is larger */
      if( SCIPgetExponentExprPow(expr1) < 1.0 )
         return -1;

      /* power expression is larger => expr1 > expr2 */
      return 1;
   }
   /* enforces OR12 */
   if( SCIPexprIsPower(set, expr2) )
      return -SCIPexprCompare(set, expr2, expr1);

   /* enforces OR10 */
   if( SCIPexprIsVar(set, expr1) )
      return -1;
   /* enforces OR12 */
   if( SCIPexprIsVar(set, expr2) )
      return -SCIPexprCompare(set, expr2, expr1);

   /* enforces OR11 */
   retval = strcmp(SCIPexprhdlrGetName(exprhdlr1), SCIPexprhdlrGetName(exprhdlr2));
   return retval == 0 ? 0 : retval < 0 ? -1 : 1;
}

/** simplifies an expression
 *
 * @see SCIPsimplifyExpr
 */
SCIP_RETCODE SCIPexprSimplify(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            rootexpr,           /**< expression to be simplified */
   SCIP_EXPR**           simplified,         /**< buffer to store simplified expression */
   SCIP_Bool*            changed,            /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*            infeasible,         /**< buffer to store whether infeasibility has been detected */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPRITER* it;

   assert(rootexpr != NULL);
   assert(simplified != NULL);
   assert(changed != NULL);
   assert(infeasible != NULL);

   /* simplify bottom up
    * when leaving an expression it simplifies it and stores the simplified expr in its iterators expression data
    * after the child was visited, it is replaced with the simplified expr
    */
   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );  /* TODO can we set allowrevisited to FALSE?*/
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITEDCHILD | SCIP_EXPRITER_LEAVEEXPR);

   *changed = FALSE;
   *infeasible = FALSE;
   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITEDCHILD:
         {
            SCIP_EXPR* newchild;
            SCIP_EXPR* child;

            newchild = (SCIP_EXPR*)SCIPexpriterGetChildUserDataDFS(it).ptrval;
            child = SCIPexpriterGetChildExprDFS(it);
            assert(newchild != NULL);

            /* if child got simplified, replace it with the new child */
            if( newchild != child )
            {
               SCIP_CALL( SCIPexprReplaceChild(set, stat, blkmem, expr, SCIPexpriterGetChildIdxDFS(it), newchild) );
            }

            /* we do not need to hold newchild anymore */
            SCIP_CALL( SCIPexprRelease(set, stat, blkmem, &newchild) );

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR:
         {
            SCIP_EXPR* refexpr = NULL;
            SCIP_EXPRITER_USERDATA iterdata;

            /* TODO we should do constant folding (handle that all children are value-expressions) here in a generic way
             * instead of reimplementing it in every handler
             */

            /* use simplification of expression handlers */
            SCIP_CALL( SCIPexprhdlrSimplifyExpr(expr->exprhdlr, set, expr, &refexpr, ownercreate, ownercreatedata) );
            assert(refexpr != NULL);
            if( expr != refexpr )
               *changed = TRUE;

            iterdata.ptrval = (void*) refexpr;
            SCIPexpriterSetCurrentUserData(it, iterdata);

            break;
         }

         default:
            SCIPABORT(); /* we should never be called in this stage */
            break;
      }
   }

   *simplified = (SCIP_EXPR*)SCIPexpriterGetExprUserData(it, rootexpr).ptrval;
   assert(*simplified != NULL);

   SCIPexpriterFree(&it);

   return SCIP_OKAY;
}

/** checks whether an expression is quadratic
 *
 * An expression is quadratic if it is either a power expression with exponent 2.0, a product of two expressions,
 * or a sum of terms where at least one is a square or a product of two.
 *
 * Use \ref SCIPexprGetQuadraticData to get data about the representation as quadratic.
 */
SCIP_RETCODE SCIPexprCheckQuadratic(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool*            isquadratic         /**< buffer to store result */
   )
{
   SCIP_HASHMAP* expr2idx;
   SCIP_HASHMAP* seenexpr = NULL;
   int nquadterms = 0;
   int nlinterms = 0;
   int nbilinterms = 0;
   int c;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(isquadratic != NULL);

   if( expr->quadchecked )
   {
      *isquadratic = expr->quaddata != NULL;
      return SCIP_OKAY;
   }
   assert(expr->quaddata == NULL);

   expr->quadchecked = TRUE;
   *isquadratic = FALSE;

   /* check if expression is a quadratic expression */
   SCIPsetDebugMsg(set, "checking if expr %p is quadratic\n", (void*)expr);

   /* handle single square term */
   if( SCIPexprIsPower(set, expr) && SCIPgetExponentExprPow(expr) == 2.0 )
   {
      SCIPsetDebugMsg(set, "expr %p looks like square: fill data structures\n", (void*)expr);
      SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, &expr->quaddata) );

      expr->quaddata->nquadexprs = 1;
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &expr->quaddata->quadexprterms, 1) );
      expr->quaddata->quadexprterms[0].expr = expr->children[0];
      expr->quaddata->quadexprterms[0].sqrcoef = 1.0;

      expr->quaddata->allexprsarevars = SCIPexprIsVar(set, expr->quaddata->quadexprterms[0].expr);

      *isquadratic = TRUE;
      return SCIP_OKAY;
   }

   /* handle single bilinear term */
   if( SCIPexprIsProduct(set, expr) && SCIPexprGetNChildren(expr) == 2 )
   {
      SCIPsetDebugMsg(set, "expr %p looks like bilinear product: fill data structures\n", (void*)expr);
      SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, &expr->quaddata) );
      expr->quaddata->nquadexprs = 2;

      SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &expr->quaddata->quadexprterms, 2) );
      expr->quaddata->quadexprterms[0].expr = SCIPexprGetChildren(expr)[0];
      expr->quaddata->quadexprterms[0].nadjbilin = 1;
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &expr->quaddata->quadexprterms[0].adjbilin, 1) );
      expr->quaddata->quadexprterms[0].adjbilin[0] = 0;

      expr->quaddata->quadexprterms[1].expr = SCIPexprGetChildren(expr)[1];
      expr->quaddata->quadexprterms[1].nadjbilin = 1;
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &expr->quaddata->quadexprterms[1].adjbilin, 1) );
      expr->quaddata->quadexprterms[1].adjbilin[0] = 0;

      expr->quaddata->nbilinexprterms = 1;
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &expr->quaddata->bilinexprterms, 1) );
      expr->quaddata->bilinexprterms[0].expr1 = SCIPexprGetChildren(expr)[1];
      expr->quaddata->bilinexprterms[0].expr2 = SCIPexprGetChildren(expr)[0];
      expr->quaddata->bilinexprterms[0].coef = 1.0;

      expr->quaddata->allexprsarevars = SCIPexprIsVar(set, expr->quaddata->quadexprterms[0].expr)
            && SCIPexprIsVar(set, expr->quaddata->quadexprterms[1].expr);

      *isquadratic = TRUE;
      return SCIP_OKAY;
   }

   /* neither a sum, nor a square, nor a bilinear term */
   if( !SCIPexprIsSum(set, expr) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapCreate(&seenexpr, blkmem, 2*SCIPexprGetNChildren(expr)) );
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      SCIP_EXPR* child;

      child = SCIPexprGetChildren(expr)[c];
      assert(child != NULL);

      if( SCIPexprIsPower(set, child) && SCIPgetExponentExprPow(child) == 2.0 ) /* quadratic term */
      {
         SCIP_CALL( quadDetectProcessExpr(SCIPexprGetChildren(child)[0], seenexpr, &nquadterms, &nlinterms) );
      }
      else if( SCIPexprIsProduct(set, child) && SCIPexprGetNChildren(child) == 2 ) /* bilinear term */
      {
         ++nbilinterms;
         SCIP_CALL( quadDetectProcessExpr(SCIPexprGetChildren(child)[0], seenexpr, &nquadterms, &nlinterms) );
         SCIP_CALL( quadDetectProcessExpr(SCIPexprGetChildren(child)[1], seenexpr, &nquadterms, &nlinterms) );
      }
      else
      {
         /* first time seen linearly --> assign -1; ++nlinterms
          * not first time --> assign +=1;
          */
         if( SCIPhashmapExists(seenexpr, (void*)child) )
         {
            assert(SCIPhashmapGetImageInt(seenexpr, (void*)child) > 0);

            SCIP_CALL( SCIPhashmapSetImageInt(seenexpr, (void*)child, SCIPhashmapGetImageInt(seenexpr, (void*)child) + 1) );
         }
         else
         {
            ++nlinterms;
            SCIP_CALL( SCIPhashmapInsertInt(seenexpr, (void*)child, -1) );
         }
      }
   }

   if( nquadterms == 0 )
   {
      /* only linear sum */
      SCIPhashmapFree(&seenexpr);
      return SCIP_OKAY;
   }

   SCIPsetDebugMsg(set, "expr %p looks quadratic: fill data structures\n", (void*)expr);

   /* expr2idx maps expressions to indices; if index > 0, it is its index in the linexprs array, otherwise -index-1 is
    * its index in the quadexprterms array
    */
   SCIP_CALL( SCIPhashmapCreate(&expr2idx, blkmem, nquadterms + nlinterms) );

   /* allocate memory, etc */
   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, &expr->quaddata) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &expr->quaddata->quadexprterms, nquadterms) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &expr->quaddata->linexprs, nlinterms) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &expr->quaddata->lincoefs, nlinterms) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &expr->quaddata->bilinexprterms, nbilinterms) );

   expr->quaddata->constant = SCIPgetConstantExprSum(expr);

   expr->quaddata->allexprsarevars = TRUE;
   /* for every term of the sum-expr */
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      SCIP_EXPR* child;
      SCIP_Real coef;

      child = SCIPexprGetChildren(expr)[c];
      coef = SCIPgetCoefsExprSum(expr)[c];

      assert(child != NULL);
      assert(coef != 0.0);

      if( SCIPexprIsPower(set, child) && SCIPgetExponentExprPow(child) == 2.0 ) /* quadratic term */
      {
         SCIP_QUADEXPR_QUADTERM* quadexprterm;
         assert(SCIPexprGetNChildren(child) == 1);

         child = SCIPexprGetChildren(child)[0];
         assert(SCIPhashmapGetImageInt(seenexpr, (void *)child) > 0);

         SCIP_CALL( quadDetectGetQuadexprterm(blkmem, child, expr2idx, seenexpr, expr->quaddata, &quadexprterm) );
         assert(quadexprterm->expr == child);
         quadexprterm->sqrcoef = coef;
         quadexprterm->sqrexpr = SCIPexprGetChildren(expr)[c];

         if( expr->quaddata->allexprsarevars )
            expr->quaddata->allexprsarevars = SCIPexprIsVar(set, quadexprterm->expr);
      }
      else if( SCIPexprIsProduct(set, child) && SCIPexprGetNChildren(child) == 2 ) /* bilinear term */
      {
         SCIP_QUADEXPR_BILINTERM* bilinexprterm;
         SCIP_QUADEXPR_QUADTERM* quadexprterm;
         SCIP_EXPR* expr1;
         SCIP_EXPR* expr2;

         assert(SCIPgetCoefExprProduct(child) == 1.0);

         expr1 = SCIPexprGetChildren(child)[0];
         expr2 = SCIPexprGetChildren(child)[1];
         assert(expr1 != NULL && expr2 != NULL);

         bilinexprterm = &expr->quaddata->bilinexprterms[expr->quaddata->nbilinexprterms];

         bilinexprterm->coef = coef;
         if( SCIPhashmapGetImageInt(seenexpr, (void*)expr1) >= SCIPhashmapGetImageInt(seenexpr, (void*)expr2) )
         {
            bilinexprterm->expr1 = expr1;
            bilinexprterm->expr2 = expr2;
         }
         else
         {
            bilinexprterm->expr1 = expr2;
            bilinexprterm->expr2 = expr1;
         }
         bilinexprterm->prodexpr = child;

         SCIP_CALL( quadDetectGetQuadexprterm(blkmem, expr1, expr2idx, seenexpr, expr->quaddata, &quadexprterm) );
         assert(quadexprterm->expr == expr1);
         quadexprterm->adjbilin[quadexprterm->nadjbilin] = expr->quaddata->nbilinexprterms;
         ++quadexprterm->nadjbilin;

         if( expr->quaddata->allexprsarevars )
            expr->quaddata->allexprsarevars = SCIPexprIsVar(set, quadexprterm->expr);

         SCIP_CALL( quadDetectGetQuadexprterm(blkmem, expr2, expr2idx, seenexpr, expr->quaddata, &quadexprterm) );
         assert(quadexprterm->expr == expr2);
         quadexprterm->adjbilin[quadexprterm->nadjbilin] = expr->quaddata->nbilinexprterms;
         ++quadexprterm->nadjbilin;

         if( expr->quaddata->allexprsarevars )
            expr->quaddata->allexprsarevars = SCIPexprIsVar(set, quadexprterm->expr);

         ++expr->quaddata->nbilinexprterms;

         /* store position of second factor in quadexprterms */
         bilinexprterm->pos2 = SCIPhashmapGetImageInt(expr2idx, (void*)bilinexprterm->expr2);
      }
      else /* linear term */
      {
         if( SCIPhashmapGetImageInt(seenexpr, (void*)child) < 0 )
         {
            assert(SCIPhashmapGetImageInt(seenexpr, (void*)child) == -1);

            /* expression only appears linearly */
            expr->quaddata->linexprs[expr->quaddata->nlinexprs] = child;
            expr->quaddata->lincoefs[expr->quaddata->nlinexprs] = coef;
            expr->quaddata->nlinexprs++;

            if( expr->quaddata->allexprsarevars )
               expr->quaddata->allexprsarevars = SCIPexprIsVar(set, child);
         }
         else
         {
            /* expression appears non-linearly: set lin coef */
            SCIP_QUADEXPR_QUADTERM* quadexprterm;
            assert(SCIPhashmapGetImageInt(seenexpr, (void*)child) > 0);

            SCIP_CALL( quadDetectGetQuadexprterm(blkmem, child, expr2idx, seenexpr, expr->quaddata, &quadexprterm) );
            assert(quadexprterm->expr == child);
            quadexprterm->lincoef = coef;

            if( expr->quaddata->allexprsarevars )
               expr->quaddata->allexprsarevars = SCIPexprIsVar(set, quadexprterm->expr);
         }
      }
   }
   assert(expr->quaddata->nquadexprs == nquadterms);
   assert(expr->quaddata->nlinexprs == nlinterms);
   assert(expr->quaddata->nbilinexprterms == nbilinterms);

   SCIPhashmapFree(&seenexpr);
   SCIPhashmapFree(&expr2idx);

   *isquadratic = TRUE;

   return SCIP_OKAY;
}

/** frees information on quadratic representation of an expression
 *
 * Reverts SCIPexprCheckQuadratic().
 * Before doing changes to an expression, it can be useful to call this function.
 */
void SCIPexprFreeQuadratic(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int i;
   int n;

   assert(blkmem != NULL);
   assert(expr != NULL);

   expr->quadchecked = FALSE;

   if( expr->quaddata == NULL )
      return;

   n = expr->quaddata->nquadexprs;

   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->linexprs, expr->quaddata->nlinexprs);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->lincoefs, expr->quaddata->nlinexprs);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->bilinexprterms, expr->quaddata->nbilinexprterms);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->eigenvalues, n);
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->eigenvectors, n * n);  /*lint !e647*/

   for( i = 0; i < n; ++i )
   {
      BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->quadexprterms[i].adjbilin,
         expr->quaddata->quadexprterms[i].adjbilinsize);
   }
   BMSfreeBlockMemoryArrayNull(blkmem, &expr->quaddata->quadexprterms, n);

   BMSfreeBlockMemory(blkmem, &expr->quaddata);
}

/** Checks the curvature of the quadratic function stored in quaddata
 *
 * For this, it builds the matrix Q of quadratic coefficients and computes its eigenvalues using LAPACK.
 * If Q is
 * - semidefinite positive -> curv is set to convex,
 * - semidefinite negative -> curv is set to concave,
 * - otherwise -> curv is set to unknown.
 *
 * If `assumevarfixed` is given and some expressions in quadratic terms correspond to variables present in
 * this hashmap, then the corresponding rows and columns are ignored in the matrix Q.
 */
SCIP_RETCODE SCIPexprComputeQuadraticCurvature(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_EXPRCURV*        curv,               /**< pointer to store the curvature of quadratics */
   SCIP_HASHMAP*         assumevarfixed,     /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   SCIP_Bool             storeeigeninfo      /**< whether the eigenvalues and eigenvectors should be stored */
   )
{
   SCIP_QUADEXPR* quaddata;
   SCIP_HASHMAP* expr2matrix;
   double* matrix;
   double* alleigval;
   int nvars;
   int nn;
   int n;
   int i;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(bufmem != NULL);
   assert(messagehdlr != NULL);
   assert(expr != NULL);
   assert(curv != NULL);

   quaddata = expr->quaddata;
   assert(quaddata != NULL);

   /* do not store eigen information if we are not considering full matrix */
   if( assumevarfixed != NULL )
      storeeigeninfo = FALSE;

   if( quaddata->eigeninfostored || (quaddata->curvaturechecked && !storeeigeninfo) )
   {
      *curv = quaddata->curvature;
      /* if we are convex or concave on the full set of variables, then we will also be so on a subset */
      if( assumevarfixed == NULL || quaddata->curvature != SCIP_EXPRCURV_UNKNOWN )
         return SCIP_OKAY;
   }
   assert(quaddata->curvature == SCIP_EXPRCURV_UNKNOWN || assumevarfixed != NULL
         || (storeeigeninfo && !quaddata->eigeninfostored));

   *curv = SCIP_EXPRCURV_UNKNOWN;

   n  = quaddata->nquadexprs;

   /* do not check curvature if nn will be too large
    * we want nn * sizeof(real) to fit into an unsigned int, so n must be <= sqrt(unit_max/sizeof(real))
    * sqrt(2*214748364/8) = 7327.1475350234
    */
   if( n > 7000 )
   {
      SCIPmessageFPrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL, NULL,
            "number of quadratic variables is too large (%d) to check the curvature\n", n);
      return SCIP_OKAY;
   }

   /* TODO do some simple tests first; like diagonal entries don't change sign, etc */

   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   nn = n * n;
   assert(nn > 0);
   assert((unsigned)nn < UINT_MAX / sizeof(SCIP_Real));

   if( storeeigeninfo )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &quaddata->eigenvalues, n));
      SCIP_ALLOC( BMSallocClearBlockMemoryArray(blkmem, &quaddata->eigenvectors, nn));

      alleigval = quaddata->eigenvalues;
      matrix = quaddata->eigenvectors;
   }
   else
   {
      SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &alleigval, n) );
      SCIP_ALLOC( BMSallocClearBufferMemoryArray(bufmem, &matrix, nn) );
   }

   SCIP_CALL( SCIPhashmapCreate(&expr2matrix, blkmem, n) );

   /* fill matrix's diagonal */
   nvars = 0;
   for( i = 0; i < n; ++i )
   {
      SCIP_QUADEXPR_QUADTERM quadexprterm;

      quadexprterm = quaddata->quadexprterms[i];

      assert(!SCIPhashmapExists(expr2matrix, (void*)quadexprterm.expr));

      /* skip expr if it is a variable mentioned in assumevarfixed */
      if( assumevarfixed != NULL && SCIPexprIsVar(set, quadexprterm.expr)
          && SCIPhashmapExists(assumevarfixed, (void*)SCIPgetVarExprVar(quadexprterm.expr)) )
         continue;

      if( quadexprterm.sqrcoef == 0.0 && ! storeeigeninfo )
      {
         assert(quadexprterm.nadjbilin > 0);
         /* SCIPdebugMsg(scip, "var <%s> appears in bilinear term but is not squared
          * --> indefinite quadratic\n", SCIPvarGetName(quadexprterm.var)); */
         goto CLEANUP;
      }

      matrix[nvars * n + nvars] = quadexprterm.sqrcoef;

      /* remember row of variable in matrix */
      SCIP_CALL( SCIPhashmapInsert(expr2matrix, (void *)quadexprterm.expr, (void *)(size_t)nvars) );
      nvars++;
   }

   /* fill matrix's upper-diagonal */
   for( i = 0; i < quaddata->nbilinexprterms; ++i )
   {
      SCIP_QUADEXPR_BILINTERM bilinexprterm;
      int col;
      int row;

      bilinexprterm = quaddata->bilinexprterms[i];

      /* each factor should have been added to expr2matrix unless it corresponds to a variable mentioned in assumevarfixed */
      assert(SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr1)
            || (assumevarfixed != NULL && SCIPexprIsVar(set, bilinexprterm.expr1)
            && SCIPhashmapExists(assumevarfixed, (void*)SCIPgetVarExprVar(bilinexprterm.expr1))));
      assert(SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr2)
            || (assumevarfixed != NULL && SCIPexprIsVar(set, bilinexprterm.expr2)
            && SCIPhashmapExists(assumevarfixed, (void*)SCIPgetVarExprVar(bilinexprterm.expr2))));

      /* skip bilinear terms where at least one of the factors should be assumed to be fixed
       * (i.e., not present in expr2matrix map) */
      if( !SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr1)
          || !SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr2) )
         continue;

      row = (int)(size_t)SCIPhashmapGetImage(expr2matrix, bilinexprterm.expr1);
      col = (int)(size_t)SCIPhashmapGetImage(expr2matrix, bilinexprterm.expr2);

      assert(row != col);

      if( row < col )
         matrix[row * n + col] = bilinexprterm.coef / 2.0;
      else
         matrix[col * n + row] = bilinexprterm.coef / 2.0;
   }

   /* compute eigenvalues */
   if( SCIPcallLapackDsyevIpopt(storeeigeninfo, n, matrix, alleigval) != SCIP_OKAY )
   {
      SCIPmessagePrintWarning(messagehdlr, "Failed to compute eigenvalues of quadratic coefficient "
            "matrix --> don't know curvature\n");
      goto CLEANUP;
   }

   /* check convexity */
   if( !SCIPsetIsNegative(set, alleigval[0]) )
      *curv = SCIP_EXPRCURV_CONVEX;
   else if( !SCIPsetIsPositive(set, alleigval[n-1]) )
      *curv = SCIP_EXPRCURV_CONCAVE;

CLEANUP:
   SCIPhashmapFree(&expr2matrix);

   if( !storeeigeninfo )
   {
      BMSfreeBufferMemoryArray(bufmem, &matrix);
      BMSfreeBufferMemoryArray(bufmem, &alleigval);
   }
   else
   {
      assert(!quaddata->eigeninfostored);
      quaddata->eigeninfostored = TRUE;
   }

   /* if checked convexity on full Q matrix, then remember it
    * if indefinite on submatrix, then it will also be indefinite on full matrix, so can remember that, too */
   if( assumevarfixed == NULL || *curv == SCIP_EXPRCURV_UNKNOWN )
   {
      quaddata->curvature = *curv;
      quaddata->curvaturechecked = TRUE;
   }

   return SCIP_OKAY;
}


/* from pub_expr.h */

#ifdef NDEBUG
#undef SCIPexprGetNUses
#undef SCIPexprGetNChildren
#undef SCIPexprGetChildren
#undef SCIPexprGetHdlr
#undef SCIPexprGetData
#undef SCIPexprSetData
#undef SCIPexprGetOwnerData
#undef SCIPexprGetEvalValue
#undef SCIPexprGetEvalTag
#undef SCIPexprGetDerivative
#undef SCIPexprGetDot
#undef SCIPexprGetBardot
#undef SCIPexprGetDiffTag
#undef SCIPexprGetActivity
#undef SCIPexprGetActivityTag
#undef SCIPexprSetActivity
#undef SCIPexprGetCurvature
#undef SCIPexprSetCurvature
#undef SCIPexprIsIntegral
#undef SCIPexprSetIntegrality
#undef SCIPexprAreQuadraticExprsVariables
#endif

/** gets the number of times the expression is currently captured */
int SCIPexprGetNUses(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nuses;
}

/** gives the number of children of an expression */
int SCIPexprGetNChildren(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nchildren;
}

/** gives the children of an expression (can be NULL if no children) */
SCIP_EXPR** SCIPexprGetChildren(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->children;
}

/** gets the expression handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_EXPRHDLR* SCIPexprGetHdlr(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprhdlr;
}

/** gets the expression data of an expression */
SCIP_EXPRDATA* SCIPexprGetData(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprdata;
}

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler only.
 */
void SCIPexprSetData(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRDATA*        exprdata            /**< expression data to be set (can be NULL) */
   )
{
   assert(expr != NULL);
   assert(exprdata == NULL || expr->exprhdlr->copydata != NULL);  /* copydata must be available if there is expression data */
   assert(exprdata == NULL || expr->exprhdlr->freedata != NULL);  /* freedata must be available if there is expression data */

   expr->exprdata = exprdata;
}

/** gets the data that the owner of an expression has stored in an expression */
SCIP_EXPR_OWNERDATA* SCIPexprGetOwnerData(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->ownerdata;
}

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error)
 *
 * @see SCIPevalExpr to evaluate the expression at a given solution.
 */
SCIP_Real SCIPexprGetEvalValue(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evalvalue;
}

/** gives the evaluation tag from the last evaluation, or 0
 *
 * @see SCIPevalExpr
 */
SCIP_Longint SCIPexprGetEvalTag(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evaltag;
}

/** returns the derivative stored in an expression (or SCIP_INVALID if there was an evaluation error)
 *
 * @see SCIPevalExprGradient
 */
SCIP_Real SCIPexprGetDerivative(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->derivative;
}

/** gives the value of directional derivative from the last evaluation of a directional derivative of
 * expression (or SCIP_INVALID if there was an error)
 *
 * @see SCIPevalExprHessianDir
 */
SCIP_Real SCIPexprGetDot(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->dot;
}

/** gives the value of directional derivative from the last evaluation of a directional derivative of
 * derivative of root (or SCIP_INVALID if there was an error)
 *
 * @see SCIPevalExprHessianDir
 */
SCIP_Real SCIPexprGetBardot(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->bardot;
}

/** returns the difftag stored in an expression
 *
 * can be used to check whether partial derivative value is valid
 *
 * @see SCIPevalExprGradient
 */
SCIP_Longint SCIPexprGetDiffTag(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->difftag;
}

/** returns the activity that is currently stored for an expression
 *
 * @see SCIPevalExprActivity
 */
SCIP_INTERVAL SCIPexprGetActivity(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->activity;
}

/** returns the tag associated with the activity of the expression
 *
 * It can depend on the owner of the expression how to interpret this tag.
 * SCIPevalExprActivity() compares with `stat->domchgcount`.
 *
 * @see SCIPevalExprActivity
 */
SCIP_Longint SCIPexprGetActivityTag(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->activitytag;
}

/** set the activity with tag for an expression */
void SCIPexprSetActivity(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_INTERVAL         activity,           /**< new activity */
   SCIP_Longint          activitytag         /**< tag associated with activity */
   )
{
   assert(expr != NULL);

   expr->activity = activity;
   expr->activitytag = activitytag;
}

/** returns the curvature of an expression
 *
 *  @note Call SCIPcomputeExprCurvature() before calling this function.
 */
SCIP_EXPRCURV SCIPexprGetCurvature(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->curvature;
}

/** sets the curvature of an expression */
void SCIPexprSetCurvature(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRCURV         curvature           /**< curvature of the expression */
   )
{
   assert(expr != NULL);

   expr->curvature = curvature;
}

/** returns whether an expression is integral */
SCIP_Bool SCIPexprIsIntegral(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->isintegral;
}

/** sets the integrality flag of an expression */
void SCIPexprSetIntegrality(
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool             isintegral          /**< integrality of the expression */
   )
{
   assert(expr != NULL);

   expr->isintegral = isintegral;
}

/** gives the coefficients and expressions that define a quadratic expression
 *
 * It can return the constant part, the number, arguments, and coefficients of the purely linear part
 * and the number of quadratic terms and bilinear terms.
 * Note that for arguments that appear in the quadratic part, a linear coefficient is
 * stored with the quadratic term.
 * Use SCIPexprGetQuadraticQuadTerm() and SCIPexprGetQuadraticBilinTerm()
 * to access the data for a quadratic or bilinear term.
 *
 * It can also return the eigenvalues and the eigenvectors of the matrix \f$Q\f$ when the quadratic is written
 * as \f$x^T Q x + b^T x + c^T y + d\f$, where \f$c^T y\f$ defines the purely linear part.
 * Note, however, that to have access to them one needs to call SCIPcomputeExprQuadraticCurvature()
 * with `storeeigeninfo=TRUE`. If the eigen information was not stored or it failed to be computed,
 * `eigenvalues` and `eigenvectors` will be set to NULL.
 *
 * This function returns pointers to internal data in linexprs and lincoefs.
 * The user must not change this data.
 *
 * @attention SCIPcheckExprQuadratic() needs to be called first to check whether expression is quadratic and initialize the data of the quadratic representation.
 */
void SCIPexprGetQuadraticData(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_Real*            constant,           /**< buffer to store constant term, or NULL */
   int*                  nlinexprs,          /**< buffer to store number of expressions that appear linearly, or NULL */
   SCIP_EXPR***          linexprs,           /**< buffer to store pointer to array of expressions that appear linearly,
                                                  or NULL */
   SCIP_Real**           lincoefs,           /**< buffer to store pointer to array of coefficients of expressions that
                                                  appear linearly, or NULL */
   int*                  nquadexprs,         /**< buffer to store number of expressions in quadratic terms, or NULL */
   int*                  nbilinexprs,        /**< buffer to store number of bilinear expressions terms, or NULL */
   SCIP_Real**           eigenvalues,        /**< buffer to store pointer to array of eigenvalues of Q, or NULL */
   SCIP_Real**           eigenvectors        /**< buffer to store pointer to array of eigenvectors of Q, or NULL */
   )
{
   SCIP_QUADEXPR* quaddata;

   assert(expr != NULL);

   quaddata = expr->quaddata;
   assert(quaddata != NULL);

   if( constant != NULL )
      *constant = quaddata->constant;
   if( nlinexprs != NULL )
      *nlinexprs = quaddata->nlinexprs;
   if( linexprs != NULL )
      *linexprs = quaddata->linexprs;
   if( lincoefs != NULL )
      *lincoefs = quaddata->lincoefs;
   if( nquadexprs != NULL )
      *nquadexprs = quaddata->nquadexprs;
   if( nbilinexprs != NULL )
      *nbilinexprs = quaddata->nbilinexprterms;
   if( eigenvalues != NULL )
      *eigenvalues = quaddata->eigenvalues;
   if( eigenvectors != NULL )
      *eigenvectors = quaddata->eigenvectors;
}

/** gives the data of a quadratic expression term
 *
 * For a term \f$a \cdot \text{expr}^2 + b \cdot \text{expr} + \sum_i (c_i \cdot \text{expr} \cdot \text{otherexpr}_i)\f$, returns
 * `expr`, \f$a\f$, \f$b\f$, the number of summands, and indices of bilinear terms in the quadratic expressions `bilinexprterms`.
 *
 * This function returns pointers to internal data in adjbilin.
 * The user must not change this data.
 */
void SCIPexprGetQuadraticQuadTerm(
   SCIP_EXPR*            quadexpr,           /**< quadratic expression */
   int                   termidx,            /**< index of quadratic term */
   SCIP_EXPR**           expr,               /**< buffer to store pointer to argument expression (the 'x') of this term,
                                                  or NULL */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of variable, or NULL */
   SCIP_Real*            sqrcoef,            /**< buffer to store square coefficient of variable, or NULL */
   int*                  nadjbilin,          /**< buffer to store number of bilinear terms this variable is involved in,
                                                  or NULL */
   int**                 adjbilin,           /**< buffer to store pointer to indices of associated bilinear terms, or NULL */
   SCIP_EXPR**           sqrexpr             /**< buffer to store pointer to square expression (the 'x^2') of this term
                                                  or NULL if no square expression, or NULL */
   )
{
   SCIP_QUADEXPR_QUADTERM* quadexprterm;

   assert(quadexpr != NULL);
   assert(quadexpr->quaddata != NULL);
   assert(quadexpr->quaddata->quadexprterms != NULL);
   assert(termidx >= 0);
   assert(termidx < quadexpr->quaddata->nquadexprs);

   quadexprterm = &quadexpr->quaddata->quadexprterms[termidx];

   if( expr != NULL )
      *expr = quadexprterm->expr;
   if( lincoef != NULL )
      *lincoef = quadexprterm->lincoef;
   if( sqrcoef != NULL )
      *sqrcoef = quadexprterm->sqrcoef;
   if( nadjbilin != NULL )
      *nadjbilin = quadexprterm->nadjbilin;
   if( adjbilin != NULL )
      *adjbilin = quadexprterm->adjbilin;
   if( sqrexpr != NULL )
      *sqrexpr = quadexprterm->sqrexpr;
}

/** gives the data of a bilinear expression term
 *
 * For a term a*expr1*expr2, returns expr1, expr2, a, and
 * the position of the quadratic expression term that uses expr2 in the quadratic expressions `quadexprterms`.
 */
void SCIPexprGetQuadraticBilinTerm(
   SCIP_EXPR*            expr,               /**< quadratic expression */
   int                   termidx,            /**< index of bilinear term */
   SCIP_EXPR**           expr1,              /**< buffer to store first factor, or NULL */
   SCIP_EXPR**           expr2,              /**< buffer to store second factor, or NULL */
   SCIP_Real*            coef,               /**< buffer to coefficient, or NULL */
   int*                  pos2,               /**< buffer to position of expr2 in quadexprterms array of quadratic
                                                  expression, or NULL */
   SCIP_EXPR**           prodexpr            /**< buffer to store pointer to expression that is product if first
                                                  and second factor, or NULL */
   )
{
   SCIP_QUADEXPR_BILINTERM* bilinexprterm;

   assert(expr != NULL);
   assert(expr->quaddata != NULL);
   assert(expr->quaddata->bilinexprterms != NULL);
   assert(termidx >= 0);
   assert(termidx < expr->quaddata->nbilinexprterms);

   bilinexprterm = &expr->quaddata->bilinexprterms[termidx];

   if( expr1 != NULL )
      *expr1 = bilinexprterm->expr1;
   if( expr2 != NULL )
      *expr2 = bilinexprterm->expr2;
   if( coef != NULL )
      *coef = bilinexprterm->coef;
   if( pos2 != NULL )
      *pos2 = bilinexprterm->pos2;
   if( prodexpr != NULL )
      *prodexpr = bilinexprterm->prodexpr;
}

/** returns whether all expressions that are used in a quadratic expression are variable expressions
 *
 * @return TRUE iff all `linexprs` and `quadexprterms[.].expr` are variable expressions
 */
SCIP_Bool SCIPexprAreQuadraticExprsVariables(
   SCIP_EXPR*            expr                /**< quadratic expression */
   )
{
   assert(expr != NULL);
   assert(expr->quaddata != NULL);

   return expr->quaddata->allexprsarevars;
}

/**@} */
