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

/**@file    prop_genvbounds.h
 * @ingroup PROPAGATORS
 * @brief   generalized variable bounds propagator
 * @author  Stefan Weltge
 * @author  Ambros Gleixner
 *
 *  A generalized variable bound is a linear inequality of the form
 *  \f[
 *     c \, x_i \geq \sum (a_j \, x_j) + d \cdot \mbox{primal\_bound} + \mbox{const},
 *  \f]
 *  where \f$c\f$ is either 1 or -1 and \f$primal\_bound\f$ is an upper bound on the optimal objective
 *  value, which may improve during the solving process. In SCIP, generalized variable bounds are
 *  used for providing bounds on the LHS's variable \f$x_i\f$. If the above inequality is valid, the
 *  following bounds, depending on \f$x_i\f$'s coefficient, are also valid:
 *  \f[
 *     c = 1   \qquad\Rightarrow\qquad   x_i \geq  \mbox{minactivity}(\sum a_j \, x_j)
 *                                       + d \cdot \mbox{primal\_bound} + \mbox{const}
 *  \f]
 *  \f[
 *     c = -1  \qquad\Rightarrow\qquad   x_i \leq - \mbox{minactivity}(\sum a_j \, x_j)
 *                                       - d \cdot \mbox{primal\_bound} - \mbox{const}.
 *  \f]
 *
 *  Note that for feasible problems, \f$d \leq 0\f$ must hold. If \f$d < 0\f$ a decrease of the
 *  primal bound causes an improvement of the provided bound. Similarly, if \f$a_j > 0\f$ (\f$< 0\f$), a
 *  tightened lower (upper) bound of a variable \f$x_j\f$ also yields a better bound for \f$x_i\f$.
 *
 *  The genvbounds propagator sorts its stored generalized variable bounds topologically in the
 *  following order: A generalized variable bound A (\f$c\, x_i \geq \ldots\f$) preceeds a
 *  generalized variable bound B if the left-hand side variable of A appears in the right-hand side
 *  of B with sign of its coefficient equal to c; i.e., if A is propagated and tightens the
 *  corresponding bound of x_i, then the minactivity on the right-hand side of B increases. We
 *  assume that this order is acyclic for the generalized variable bounds added. Under this
 *  condition, propagating the generalized variable bounds in a topological order ensures that all
 *  propagations are found in one round.
 *
 *  Both global and local propagation is applied: If the primal bound improves, generalized variable bounds with a
 *  nonzero coefficient d are enforced in order to tighten global bounds using the global variable bounds for computing
 *  the minactivity. Independently, the genvbounds propagator catches events SCIP_EVENTTYPE_LBTIGHTENED and
 *  SCIP_EVENTTYPE_UBTIGHTENED, i.e., locally tightened bounds of variables that occur in the right-hand sides of
 *  generalized variable bounds, in order to perform an efficient local propagation when called.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_GENVBOUNDS_H__
#define __SCIP_PROP_GENVBOUNDS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PROPAGATORS
  *
  * @{
  */

/** adds a generalized variable bound to the genvbounds propagator; if there is already a genvbound for the bound
 *  "boundtype" of variable "var", it will be replaced
 */
EXTERN
SCIP_RETCODE SCIPgenVBoundAdd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            genvboundprop,      /**< genvbound propagator */
   SCIP_VAR**            vars,               /**< array of RHSs variables */
   SCIP_VAR*             var,                /**< LHSs variable */
   SCIP_Real*            coefs,              /**< array of coefficients for the RHSs variables */
   int                   ncoefs,             /**< size of coefs array */
   SCIP_Real             coefprimalbound,    /**< nonpositive value of the primal bounds multiplier */
   SCIP_Real             constant,           /**< constant term */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound provided by the genvbound */
   );

/* @} */


/** creates the genvbounds propagator and includes it in SCIP
 *
 * @ingroup PropagatorIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludePropGenvbounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
