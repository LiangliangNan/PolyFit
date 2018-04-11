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

/**@file   type_symmetry.h
 * @brief  type definitions for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SYMMETRY_H_
#define __SCIP_TYPE_SYMMETRY_H_

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** symmetry type specification */
#define SYM_SPEC_INTEGER                UINT32_C(0x00000001)  /**< need symmetries for integer variables only */
#define SYM_SPEC_BINARY                 UINT32_C(0x00000002)  /**< need symmetries for binary variables only */
#define SYM_SPEC_REAL                   UINT32_C(0x00000004)  /**< need symmetries also for continuous variables */

typedef uint32_t SYM_SPEC;              /**< types of variables handled by symmetry */

/** define sense of rhs */
enum SYM_Rhssense
{
   SYM_SENSE_UNKOWN     = 0,                 /**< unknown sense */
   SYM_SENSE_INEQUALITY = 1,                 /**< linear inequality */
   SYM_SENSE_EQUATION   = 2,                 /**< linear equation */
   SYM_SENSE_XOR        = 3,                 /**< XOR constraint */
   SYM_SENSE_AND        = 4,                 /**< AND constraint */
   SYM_SENSE_OR         = 5                  /**< OR constrant */
};
typedef enum SYM_Rhssense SYM_RHSSENSE;

/* type of symmetry handling codes */
#define SYM_HANDLETYPE_NONE             UINT32_C(0x00000000)  /**< no symmetry handling */
#define SYM_HANDLETYPE_SYMBREAK         UINT32_C(0x00000001)  /**< symmetry breaking inequalities */
#define SYM_HANDLETYPE_ORBITALFIXING    UINT32_C(0x00000002)  /**< orbital fixing */

typedef uint32_t SYM_HANDLETYPE;        /**< type of symmetry handling */

typedef struct SYM_Vartype SYM_VARTYPE;      /**< data of variables that are considered to be equivalent */
typedef struct SYM_Matrixdata SYM_MATRIXDATA;/**< data for symmetry group computation */

#ifdef __cplusplus
}
#endif

#endif
