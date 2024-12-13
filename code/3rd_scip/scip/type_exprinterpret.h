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

/**@file   type_exprinterpret.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for expression interpreter
 * @author Stefan Vigerske
 */

/** @defgroup DEFPLUGINS_EXPRINT Default expression interpreter
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c/.cpp files) of the default expression handlers of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_EXPRINTERPRET_H__
#define __SCIP_TYPE_EXPRINTERPRET_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_ExprInt     SCIP_EXPRINT;      /**< an expression interpreter */
typedef struct SCIP_ExprIntData SCIP_EXPRINTDATA;  /**< data of an expression interpreter */
typedef unsigned int            SCIP_EXPRINTCAPABILITY; /**< type of expression interpreter capability */

#define SCIP_EXPRINTCAPABILITY_NONE         0x00000000  /**< the expression interpreter is capable of nothing */
#define SCIP_EXPRINTCAPABILITY_FUNCVALUE    0x00000001  /**< the expression interpreter is able to compute a function value in a point */
#define SCIP_EXPRINTCAPABILITY_GRADIENT     0x00000010  /**< the expression interpreter is able to compute a gradient in a point */
#define SCIP_EXPRINTCAPABILITY_HESSIAN      0x00000100  /**< the expression interpreter is able to compute a full hessian in a point */
#define SCIP_EXPRINTCAPABILITY_ALL          (SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_HESSIAN)

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_EXPRINTERPRET_H__ */
