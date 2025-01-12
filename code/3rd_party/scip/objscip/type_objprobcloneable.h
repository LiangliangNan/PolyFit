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
#define SCIP_DECL_BENDERSCLONE(x) x (SCIP* scip, SCIP_Bool* valid) const

/** returns whether the plugin object is copyable
 *
 *  return value      : whether object is copyable
 */
#define SCIP_DECL_OBJPROBISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_CONSHDLRISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_PRICERISCLONEABLE(x) SCIP_Bool x (void) const
#define SCIP_DECL_BENDERSISCLONEABLE(x) SCIP_Bool x (void) const

#endif
