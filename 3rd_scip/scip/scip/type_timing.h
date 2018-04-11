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

/**@file   type_timing.h
 * @brief  timing definitions for SCIP
 * @author Timo Berthold
 * @author Matthias Miltenberger
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_TIMING_H__
#define __SCIP_TYPE_TIMING_H__

#include "scip/def.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** presolving execution timing flags
 *
 *  @note: in contrast to other timings, these timings need to start from 0x0002u in order to avoid confusion with
 *         the now obsolete '(presol)delay' boolean flag used until SCIP version 3.1.1
 */
#define SCIP_PRESOLTIMING_NONE            0x002u  /**< presolving disabled */
#define SCIP_PRESOLTIMING_FAST            0x004u  /**< timing for fast presolving methods */
#define SCIP_PRESOLTIMING_MEDIUM          0x008u  /**< timing for more expensive presolving methods */
#define SCIP_PRESOLTIMING_EXHAUSTIVE      0x010u  /**< timing for most expensive presolving methods */
#define SCIP_PRESOLTIMING_FINAL           0x020u  /**< timing for final presolving methods */

/** call presolver in every timing */
#define SCIP_PRESOLTIMING_ALWAYS (SCIP_PRESOLTIMING_FAST | SCIP_PRESOLTIMING_MEDIUM | SCIP_PRESOLTIMING_EXHAUSTIVE )
#define SCIP_PRESOLTIMING_MAX    (SCIP_PRESOLTIMING_FAST | SCIP_PRESOLTIMING_MEDIUM | SCIP_PRESOLTIMING_EXHAUSTIVE | SCIP_PRESOLTIMING_FINAL)

typedef unsigned int SCIP_PRESOLTIMING;


/** propagation execution timing flags */
#define SCIP_PROPTIMING_BEFORELP          0x001u  /**< call propagator before LP is solved */
#define SCIP_PROPTIMING_DURINGLPLOOP      0x002u  /**< call propagator after each LP solving during cut-and-price loop */
#define SCIP_PROPTIMING_AFTERLPLOOP       0x004u  /**< call propagator after the cut-and-price loop was finished */
#define SCIP_PROPTIMING_AFTERLPNODE       0x008u  /**< call propagator after the processing of a node with solved LP was
                                                   *   finished */

/** call propagator regardless of current status */
#define SCIP_PROPTIMING_ALWAYS (SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_DURINGLPLOOP | SCIP_PROPTIMING_AFTERLPLOOP | SCIP_PROPTIMING_AFTERLPNODE )

typedef unsigned int SCIP_PROPTIMING;


/** heuristics execution timing flags */
#define SCIP_HEURTIMING_BEFORENODE        0x001u  /**< call heuristic before the processing of the node starts */
#define SCIP_HEURTIMING_DURINGLPLOOP      0x002u  /**< call heuristic after each LP solving during cut-and-price loop */
#define SCIP_HEURTIMING_AFTERLPLOOP       0x004u  /**< call heuristic after the cut-and-price loop was finished */
#define SCIP_HEURTIMING_AFTERLPNODE       0x008u  /**< call heuristic after the processing of a node with solved LP was
                                                   *   finished */
#define SCIP_HEURTIMING_AFTERPSEUDONODE   0x010u  /**< call heuristic after the processing of a node without solved LP was
                                                   *   finished */
#define SCIP_HEURTIMING_AFTERLPPLUNGE     0x020u  /**< call heuristic after the processing of the last node in the current
                                                   *   plunge was finished, and only if the LP was solved for this node */
#define SCIP_HEURTIMING_AFTERPSEUDOPLUNGE 0x040u  /**< call heuristic after the processing of the last node in the current
                                                   *   plunge was finished, and only if the LP was not solved for this node */
#define SCIP_HEURTIMING_DURINGPRICINGLOOP 0x080u  /**< call heuristic during pricing loop */
#define SCIP_HEURTIMING_BEFOREPRESOL      0x100u  /**< call heuristic before presolving */
#define SCIP_HEURTIMING_DURINGPRESOLLOOP  0x200u  /**< call heuristic during presolving loop */
#define SCIP_HEURTIMING_AFTERPROPLOOP     0x400u  /**< call heuristic after propagation which is performed before solving the LP */
/* it turned out that a heuristic timing DURINGPROPLOOP causes severe troubles with the resolving of propagations */

/** call heuristic after the processing of a node was finished */
#define SCIP_HEURTIMING_AFTERNODE (SCIP_HEURTIMING_AFTERLPNODE | SCIP_HEURTIMING_AFTERPSEUDONODE)

/** call heuristic after the processing of the last node in the current plunge was finished */
#define SCIP_HEURTIMING_AFTERPLUNGE (SCIP_HEURTIMING_AFTERLPPLUNGE | SCIP_HEURTIMING_AFTERPSEUDOPLUNGE)

typedef unsigned int SCIP_HEURTIMING;

#ifdef __cplusplus
}
#endif

#endif
