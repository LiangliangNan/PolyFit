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

/**@file   type_exprinterpret.h
 * @brief  type definitions for expression interpreter
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
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
#define SCIP_EXPRINTCAPABILITY_INTFUNCVALUE 0x00000002  /**< the expression interpreter is able to compute an interval function value */
#define SCIP_EXPRINTCAPABILITY_GRADIENT     0x00000010  /**< the expression interpreter is able to compute a gradient in a point */
#define SCIP_EXPRINTCAPABILITY_INTGRADIENT  0x00000020  /**< the expression interpreter is able to compute an interval gradient */
#define SCIP_EXPRINTCAPABILITY_HESSIAN      0x00000100  /**< the expression interpreter is able to compute a full hessian in a point */
#define SCIP_EXPRINTCAPABILITY_INTHESSIAN   0x00000200  /**< the expression interpreter is able to compute an interval hessian */
#define SCIP_EXPRINTCAPABILITY_ALL          (SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_INTFUNCVALUE | SCIP_EXPRINTCAPABILITY_GRADIENT | SCIP_EXPRINTCAPABILITY_INTGRADIENT | SCIP_EXPRINTCAPABILITY_HESSIAN | SCIP_EXPRINTCAPABILITY_INTHESSIAN)

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_EXPRINTERPRET_H__ */
