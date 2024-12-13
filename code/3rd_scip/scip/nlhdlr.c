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

/**@file   nlhdlr.c
 * @ingroup OTHER_CFILES
 * @brief  functions for nonlinearity handlers of nonlinear constraint handler
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#include <assert.h>

#include "scip/pub_nlhdlr.h"
#include "scip/nlhdlr.h"
#include "scip/struct_nlhdlr.h"
#include "scip/scip_timing.h"
#include "scip/scip_mem.h"
#include "scip/scip_param.h"
#include "scip/scip_message.h"
#include "scip/pub_misc.h"

/**@addtogroup PublicNlhdlrInterfaceMethods
 * @{
 */

#ifdef NDEBUG
/* Undo the defines from pub_nlhdlr.h, which exist if NDEBUG is defined. */
#undef SCIPnlhdlrSetCopyHdlr
#undef SCIPnlhdlrSetFreeHdlrData
#undef SCIPnlhdlrSetFreeExprData
#undef SCIPnlhdlrSetInitExit
#undef SCIPnlhdlrSetProp
#undef SCIPnlhdlrSetSepa
#undef SCIPnlhdlrGetName
#undef SCIPnlhdlrGetDesc
#undef SCIPnlhdlrGetDetectPriority
#undef SCIPnlhdlrGetEnfoPriority
#undef SCIPnlhdlrIsEnabled
#undef SCIPnlhdlrGetData
#undef SCIPnlhdlrHasIntEval
#undef SCIPnlhdlrHasReverseProp
#undef SCIPnlhdlrHasInitSepa
#undef SCIPnlhdlrHasExitSepa
#undef SCIPnlhdlrHasEnfo
#undef SCIPnlhdlrHasEstimate
#endif

/** sets the copy handler callback of a nonlinear handler */
void SCIPnlhdlrSetCopyHdlr(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRCOPYHDLR((*copy))         /**< copy callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->copyhdlr = copy;
}

/** sets the nonlinear handler callback to free the nonlinear handler data */
void SCIPnlhdlrSetFreeHdlrData(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEHDLRDATA((*freehdlrdata)) /**< handler free callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->freehdlrdata = freehdlrdata;
}

/** sets the nonlinear handler callback to free expression specific data of nonlinear handler */
void SCIPnlhdlrSetFreeExprData(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEEXPRDATA((*freeexprdata)) /**< nonlinear handler expression data free callback
                                                      (can be NULL if data does not need to be freed) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->freeexprdata = freeexprdata;
}

/** sets the initialization and deinitialization callback of a nonlinear handler */
void SCIPnlhdlrSetInitExit(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINIT((*init)),            /**< initialization callback (can be NULL) */
   SCIP_DECL_NLHDLREXIT((*exit_))            /**< deinitialization callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->init = init;
   nlhdlr->exit = exit_;
}

/** sets the propagation callbacks of a nonlinear handler */
void SCIPnlhdlrSetProp(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINTEVAL((*inteval)),      /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);

   nlhdlr->inteval = inteval;
   nlhdlr->reverseprop = reverseprop;
}

/** sets the enforcement callbacks of a nonlinear handler */
void SCIPnlhdlrSetSepa(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINITSEPA((*initsepa)),    /**< separation initialization callback (can be NULL) */
   SCIP_DECL_NLHDLRENFO((*enfo)),            /**< enforcement callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_NLHDLRESTIMATE((*estimate)),    /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_NLHDLREXITSEPA((*exitsepa))     /**< separation deinitialization callback (can be NULL) */
   )
{
   assert(nlhdlr != NULL);
   assert(enfo != NULL || estimate != NULL);

   nlhdlr->initsepa = initsepa;
   nlhdlr->enfo = enfo;
   nlhdlr->estimate = estimate;
   nlhdlr->exitsepa = exitsepa;
}

/** gives name of nonlinear handler */
const char* SCIPnlhdlrGetName(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->name;
}

/** gives description of nonlinear handler, can be NULL */
const char* SCIPnlhdlrGetDesc(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->desc;
}

/** gives detection priority of nonlinear handler */
int SCIPnlhdlrGetDetectPriority(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->detectpriority;
}

/** gives enforcement priority of nonlinear handler */
int SCIPnlhdlrGetEnfoPriority(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->enfopriority;
}

/** returns whether nonlinear handler is enabled */
SCIP_Bool SCIPnlhdlrIsEnabled(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->enabled;
}

/** gives handler data of nonlinear handler */
SCIP_NLHDLRDATA* SCIPnlhdlrGetData(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);

   return nlhdlr->data;
}

/** returns whether nonlinear handler implements the interval evaluation callback */
SCIP_Bool SCIPnlhdlrHasIntEval(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->inteval != NULL;
}

/** returns whether nonlinear handler implements the reverse propagation callback */
SCIP_Bool SCIPnlhdlrHasReverseProp(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->reverseprop != NULL;
}

/** returns whether nonlinear handler implements the separation initialization callback */
SCIP_Bool SCIPnlhdlrHasInitSepa(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->initsepa != NULL;
}

/** returns whether nonlinear handler implements the separation deinitialization callback */
SCIP_Bool SCIPnlhdlrHasExitSepa(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->exitsepa != NULL;
}

/** returns whether nonlinear handler implements the enforcement callback */
SCIP_Bool SCIPnlhdlrHasEnfo(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->enfo != NULL;
}

/** returns whether nonlinear handler implements the estimator callback */
SCIP_Bool SCIPnlhdlrHasEstimate(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   return nlhdlr->estimate != NULL;
}

/** compares two nonlinear handlers by detection priority
 *
 * if handlers have same detection priority, then compare by name
 */
SCIP_DECL_SORTPTRCOMP(SCIPnlhdlrComp)
{
   SCIP_NLHDLR* h1;
   SCIP_NLHDLR* h2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   h1 = (SCIP_NLHDLR*)elem1;
   h2 = (SCIP_NLHDLR*)elem2;

   if( h1->detectpriority != h2->detectpriority )
      return h1->detectpriority - h2->detectpriority;

   return strcmp(h1->name, h2->name);
}

#ifdef SCIP_DISABLED_CODE
/** compares nonlinear handler by enforcement priority
 *
 * if handlers have same enforcement priority, then compare by detection priority, then by name
 */
SCIP_DECL_SORTPTRCOMP(SCIPnlhdlrCompEnfo)
{
   SCIP_NLHDLR* h1;
   SCIP_NLHDLR* h2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   h1 = (SCIP_NLHDLR*)elem1;
   h2 = (SCIP_NLHDLR*)elem2;

   if( h1->enfopriority != h2->enfopriority )
      return h1->enfopriority - h2->enfopriority;

   if( h1->detectpriority != h2->detectpriority )
      return h1->detectpriority - h2->detectpriority;

   return strcmp(h1->name, h2->name);
}
#endif

/** @} */

/* nlhdlr private API functions from nlhdlr.h */

#ifndef NDEBUG
#undef SCIPnlhdlrResetNDetectionslast
#undef SCIPnlhdlrIncrementNCutoffs
#undef SCIPnlhdlrIncrementNSeparated
#endif

/** creates a nonlinear handler */
SCIP_RETCODE SCIPnlhdlrCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlr,             /**< buffer to store pointer to created nonlinear handler */
   const char*           name,               /**< name of nonlinear handler (must not be NULL) */
   const char*           desc,               /**< description of nonlinear handler (can be NULL) */
   int                   detectpriority,     /**< detection priority of nonlinear handler */
   int                   enfopriority,       /**< enforcement priority of nonlinear handler */
   SCIP_DECL_NLHDLRDETECT((*detect)),        /**< structure detection callback of nonlinear handler */
   SCIP_DECL_NLHDLREVALAUX((*evalaux)),      /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_NLHDLRDATA*      nlhdlrdata          /**< data of nonlinear handler (can be NULL) */
   )
{
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(name != NULL);
   assert(detect != NULL);
   assert(evalaux != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*nlhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL_FINALLY( SCIPduplicateMemoryArray(scip, &(*nlhdlr)->desc, desc, strlen(desc)+1),
         SCIPfreeMemoryArray(scip, &(*nlhdlr)->name) );
   }

   (*nlhdlr)->detectpriority = detectpriority;
   (*nlhdlr)->enfopriority = enfopriority;
   (*nlhdlr)->data = nlhdlrdata;
   (*nlhdlr)->detect = detect;
   (*nlhdlr)->evalaux = evalaux;

   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->detecttime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->enfotime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->proptime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*nlhdlr)->intevaltime) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "nlhdlr/%s/enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, paramname, "should this nonlinear handler be used",
      &(*nlhdlr)->enabled, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}

/** frees a nonlinear handler */
SCIP_RETCODE SCIPnlhdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlr              /**< pointer to nonlinear handler to be freed */
   )
{
   assert(nlhdlr != NULL);
   assert(*nlhdlr != NULL);

   if( (*nlhdlr)->freehdlrdata != NULL )
   {
      SCIP_CALL( (*nlhdlr)->freehdlrdata(scip, *nlhdlr, &(*nlhdlr)->data) );
   }

   /* free clocks */
   SCIP_CALL( SCIPfreeClock(scip, &(*nlhdlr)->detecttime) );
   SCIP_CALL( SCIPfreeClock(scip, &(*nlhdlr)->enfotime) );
   SCIP_CALL( SCIPfreeClock(scip, &(*nlhdlr)->proptime) );
   SCIP_CALL( SCIPfreeClock(scip, &(*nlhdlr)->intevaltime) );

   SCIPfreeMemory(scip, &(*nlhdlr)->name);
   SCIPfreeMemoryNull(scip, &(*nlhdlr)->desc);

   SCIPfreeBlockMemory(scip, nlhdlr);

   return SCIP_OKAY;
}

/** call the handler copy callback of a nonlinear handler */
SCIP_DECL_NLHDLRCOPYHDLR(SCIPnlhdlrCopyhdlr)
{
   /* TODO for now just don't copy disabled nlhdlr, a clean way would probably be to first copy and disable then */
   if( sourcenlhdlr->copyhdlr != NULL && sourcenlhdlr->enabled )
   {
      SCIP_CALL( sourcenlhdlr->copyhdlr(targetscip, targetconshdlr, sourceconshdlr, sourcenlhdlr) );
   }

   return SCIP_OKAY;
}

/** call the free expression specific data callback of a nonlinear handler */
SCIP_DECL_NLHDLRFREEEXPRDATA(SCIPnlhdlrFreeexprdata)
{
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   if( nlhdlr->freeexprdata != NULL )
   {
      SCIP_CALL( nlhdlr->freeexprdata(scip, nlhdlr, expr, nlhdlrexprdata) );
      assert(*nlhdlrexprdata == NULL);
   }

   return SCIP_OKAY;
}

/** call the initialization callback of a nonlinear handler */
SCIP_DECL_NLHDLRINIT(SCIPnlhdlrInit)
{
   assert(nlhdlr != NULL);

   nlhdlr->nenfocalls = 0;
   nlhdlr->nintevalcalls = 0;
   nlhdlr->npropcalls = 0;
   nlhdlr->nseparated = 0;
   nlhdlr->ncutoffs = 0;
   nlhdlr->ndomreds = 0;
   nlhdlr->nbranchscores = 0;
   nlhdlr->ndetections = 0;
   nlhdlr->ndetectionslast = 0;

   SCIP_CALL( SCIPresetClock(scip, nlhdlr->detecttime) );
   SCIP_CALL( SCIPresetClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( SCIPresetClock(scip, nlhdlr->proptime) );
   SCIP_CALL( SCIPresetClock(scip, nlhdlr->intevaltime) );

   if( nlhdlr->init != NULL )
   {
      SCIP_CALL( nlhdlr->init(scip, nlhdlr) );
   }

   return SCIP_OKAY;
}

/** call the deinitialization callback of a nonlinear handler */
SCIP_DECL_NLHDLREXIT(SCIPnlhdlrExit)
{
   assert(nlhdlr != NULL);

   if( nlhdlr->exit != NULL )
   {
      SCIP_CALL( nlhdlr->exit(scip, nlhdlr) );
   }

   return SCIP_OKAY;
}

/** call the detect callback of a nonlinear handler */
SCIP_DECL_NLHDLRDETECT(SCIPnlhdlrDetect)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->detect != NULL);
   assert(nlhdlr->detecttime != NULL);
   assert(participating != NULL);

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->detecttime) );
   SCIP_CALL( nlhdlr->detect(scip, conshdlr, nlhdlr, expr, cons, enforcing, participating, nlhdlrexprdata) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->detecttime) );

   if( *participating != SCIP_NLHDLR_METHOD_NONE )
   {
      ++nlhdlr->ndetections;
      ++nlhdlr->ndetectionslast;
   }

   return SCIP_OKAY;
}

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLREVALAUX(SCIPnlhdlrEvalaux)
{
   assert(nlhdlr != NULL);
   assert(nlhdlr->evalaux != NULL);

   SCIP_CALL( nlhdlr->evalaux(scip, nlhdlr, expr, nlhdlrexprdata, auxvalue, sol) );

   return SCIP_OKAY;
}

/** call the interval evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLRINTEVAL(SCIPnlhdlrInteval)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->intevaltime != NULL);

   if( nlhdlr->inteval != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->intevaltime) );
      SCIP_CALL( nlhdlr->inteval(scip, nlhdlr, expr, nlhdlrexprdata, interval, intevalvar, intevalvardata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->intevaltime) );

      ++nlhdlr->nintevalcalls;
   }

   return SCIP_OKAY;
}

/** call the reverse propagation callback of a nonlinear handler */
SCIP_DECL_NLHDLRREVERSEPROP(SCIPnlhdlrReverseprop)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->proptime != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   if( nlhdlr->reverseprop == NULL )
   {
      *infeasible = FALSE;
      *nreductions = 0;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->proptime) );
   SCIP_CALL( nlhdlr->reverseprop(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, bounds, infeasible, nreductions) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->proptime) );

   /* update statistics */
   nlhdlr->ndomreds += *nreductions;
   if( *infeasible )
      ++nlhdlr->ncutoffs;
   ++nlhdlr->npropcalls;

   return SCIP_OKAY;
}

/** call the separation initialization callback of a nonlinear handler */
SCIP_DECL_NLHDLRINITSEPA(SCIPnlhdlrInitsepa)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(infeasible != NULL);

   if( nlhdlr->initsepa == NULL )
   {
      *infeasible = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->initsepa(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, overestimate, underestimate, infeasible) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   ++nlhdlr->nenfocalls;
   if( *infeasible )
      ++nlhdlr->ncutoffs;

   return SCIP_OKAY;
}

/** call the separation deinitialization callback of a nonlinear handler */
SCIP_DECL_NLHDLREXITSEPA(SCIPnlhdlrExitsepa)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);

   if( nlhdlr->exitsepa != NULL )
   {
      SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
      SCIP_CALL( nlhdlr->exitsepa(scip, nlhdlr, expr, nlhdlrexprdata) );
      SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );
   }

   return SCIP_OKAY;
}

/** call the enforcement callback of a nonlinear handler */
SCIP_DECL_NLHDLRENFO(SCIPnlhdlrEnfo)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(result != NULL);

   if( nlhdlr->enfo == NULL )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      /* we should get EXACTLY the same value from calling evalaux with the same solution as before */
      assert(auxvalue == auxvaluetest);  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->enfo(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue,
         overestimate, allowweakcuts, separated, addbranchscores, result) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;
   switch( *result )
   {
      case SCIP_SEPARATED :
         ++nlhdlr->nseparated;
         break;
      case SCIP_BRANCHED:
         ++nlhdlr->nbranchscores;
         break;
      case SCIP_CUTOFF:
         ++nlhdlr->ncutoffs;
         break;
      case SCIP_REDUCEDDOM:
         ++nlhdlr->ndomreds;
         break;
      default: ;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** call the estimator callback of a nonlinear handler */
SCIP_DECL_NLHDLRESTIMATE(SCIPnlhdlrEstimate)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlr->enfotime != NULL);
   assert(success != NULL);
   assert(addedbranchscores != NULL);

   if( nlhdlr->estimate == NULL )
   {
      *success = FALSE;
      *addedbranchscores = FALSE;
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   /* check that auxvalue is correct by reevaluating */
   {
      SCIP_Real auxvaluetest;
      SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, nlhdlrexprdata, &auxvaluetest, sol) );
      /* we should get EXACTLY the same value from calling evalaux with the same solution as before */
      assert(auxvalue == auxvaluetest);  /*lint !e777*/
   }
#endif

   SCIP_CALL( SCIPstartClock(scip, nlhdlr->enfotime) );
   SCIP_CALL( nlhdlr->estimate(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate, targetvalue, addbranchscores, rowpreps, success, addedbranchscores) );
   SCIP_CALL( SCIPstopClock(scip, nlhdlr->enfotime) );

   /* update statistics */
   ++nlhdlr->nenfocalls;

   return SCIP_OKAY;
}

/** reset number of detections counter for last round */
void SCIPnlhdlrResetNDetectionslast(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);
   nlhdlr->ndetectionslast = 0;
}

/** increments number of cutoffs in statistics */
void SCIPnlhdlrIncrementNCutoffs(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);
   ++nlhdlr->ncutoffs;
}

/** increments number of separations in statistics */
void SCIPnlhdlrIncrementNSeparated(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   )
{
   assert(nlhdlr != NULL);
   ++nlhdlr->nseparated;
}

/** print statistics for nonlinear handlers */
void SCIPnlhdlrPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlrs,            /**< nonlinear handlers */
   int                   nnlhdlrs,           /**< number of nonlinear handlers */
   FILE*                 file                /**< file handle, or NULL for standard out */
   )
{
   int i;

   SCIPinfoMessage(scip, file, "Nlhdlrs            : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
      "Detects", "DetectAll", "DetectTime",
      "#IntEval", "IntEvalTi",
      "#RevProp", "RevPropTi", "DomReds", "Cutoffs",
      "#Enforce", "EnfoTime", "Cuts", "Branching");

   for( i = 0; i < nnlhdlrs; ++i )
   {
      /* skip disabled nlhdlr */
      if( !nlhdlrs[i]->enabled )
         continue;

      SCIPinfoMessage(scip, file, "  %-17s:", nlhdlrs[i]->name);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->ndetectionslast);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->ndetections);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlrs[i]->detecttime));

      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->nintevalcalls);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlrs[i]->intevaltime));

      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->npropcalls);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlrs[i]->proptime));
      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->ndomreds);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->ncutoffs);

      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->nenfocalls);
      SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, nlhdlrs[i]->enfotime));
      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->nseparated);
      SCIPinfoMessage(scip, file, " %10lld", nlhdlrs[i]->nbranchscores);

      SCIPinfoMessage(scip, file, "\n");
   }
}
