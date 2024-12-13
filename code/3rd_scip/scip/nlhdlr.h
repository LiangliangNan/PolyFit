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

/**@file   nlhdlr.h
 * @ingroup INTERNALAPI
 * @brief  private functions of nonlinear handlers of nonlinear constraints
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_NLHDLR_H_
#define SCIP_NLHDLR_H_

#include "scip/pub_nlhdlr.h"

#ifndef NDEBUG
#include "scip/struct_nlhdlr.h"
#endif

#ifdef __cplusplus
extern "C" {
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
   );

/** frees a nonlinear handler */
SCIP_RETCODE SCIPnlhdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlr              /**< pointer to nonlinear handler to be freed */
   );

/** call the handler copy callback of a nonlinear handler */
SCIP_DECL_NLHDLRCOPYHDLR(SCIPnlhdlrCopyhdlr);

/** call the free expression specific data callback of a nonlinear handler */
SCIP_DECL_NLHDLRFREEEXPRDATA(SCIPnlhdlrFreeexprdata);

/** call the initialization callback of a nonlinear handler */
SCIP_DECL_NLHDLRINIT(SCIPnlhdlrInit);

/** call the deinitialization callback of a nonlinear handler */
SCIP_DECL_NLHDLREXIT(SCIPnlhdlrExit);

/** call the detect callback of a nonlinear handler */
SCIP_DECL_NLHDLRDETECT(SCIPnlhdlrDetect);

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLREVALAUX(SCIPnlhdlrEvalaux);

/** call the interval evaluation callback of a nonlinear handler */
SCIP_DECL_NLHDLRINTEVAL(SCIPnlhdlrInteval);

/** call the reverse propagation callback of a nonlinear handler */
SCIP_DECL_NLHDLRREVERSEPROP(SCIPnlhdlrReverseprop);

/** call the separation initialization callback of a nonlinear handler */
SCIP_DECL_NLHDLRINITSEPA(SCIPnlhdlrInitsepa);

/** call the separation deinitialization callback of a nonlinear handler */
SCIP_DECL_NLHDLREXITSEPA(SCIPnlhdlrExitsepa);

/** call the enforcement callback of a nonlinear handler */
SCIP_DECL_NLHDLRENFO(SCIPnlhdlrEnfo);

/** call the estimator callback of a nonlinear handler */
SCIP_DECL_NLHDLRESTIMATE(SCIPnlhdlrEstimate);

/** reset number of detections counter for last round */
void SCIPnlhdlrResetNDetectionslast(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** increments number of cutoffs in statistics */
void SCIPnlhdlrIncrementNCutoffs(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** increments number of separations in statistics */
void SCIPnlhdlrIncrementNSeparated(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** print statistics for nonlinear handlers */
void SCIPnlhdlrPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlrs,            /**< nonlinear handlers */
   int                   nnlhdlrs,           /**< number of nonlinear handlers */
   FILE*                 file                /**< file handle, or NULL for standard out */
   );

#ifndef NDEBUG
#define SCIPnlhdlrResetNDetectionslast(nlhdlr)  (nlhdlr)->ndetectionslast = 0
#define SCIPnlhdlrIncrementNCutoffs(nlhdlr)     ++(nlhdlr)->ncutoffs
#define SCIPnlhdlrIncrementNSeparated(nlhdlr)   ++(nlhdlr)->nseparated
#endif

#ifdef __cplusplus
}
#endif

#endif /* SCIP_NLHDLR_H_ */
