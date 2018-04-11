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

/**@file   type_objprobcloneable.h
 * @ingroup TYPEDEFINITIONS
 * @brief  function type definitions for clonable classes which define problem data
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_OBJPROBCLONEABLE_H__
#define __SCIP_TYPE_OBJPROBCLONEABLE_H__

/** clone method which will be used to copy constraint handler and variable pricer objects
 *
 *  input:
 *  - scip            : SCIP main data structure
 *
 *  output:
 *  - valid           : pointer to store whether to copy is valid w.r.t. copying dual reductions
 */
#define SCIP_DECL_OBJPROBCLONE(x) x (SCIP* scip, SCIP_Bool* valid) const
#define SCIP_DECL_CONSHDLRCLONE(x) x (SCIP* scip, SCIP_Bool* valid) const
#define SCIP_DECL_PRICERCLONE(x) x (SCIP* scip, SCIP_Bool* valid) const

/** returns whether the plugin object is copyable
 *
 *  return value      : whether object is copyable
 */
#define SCIP_DECL_OBJPROBISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_CONSHDLRISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_PRICERISCLONEABLE(x) SCIP_Bool x (void) const

#endif
