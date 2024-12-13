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

/**@file   type_symmetry.h
 * @brief  type definitions for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SYMMETRY_H_
#define __SCIP_TYPE_SYMMETRY_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** symmetry type specification */
#define SYM_SPEC_INTEGER                UINT32_C(0x00000001)  /**< need symmetries for integer variables only */
#define SYM_SPEC_BINARY                 UINT32_C(0x00000002)  /**< need symmetries for binary variables only */
#define SYM_SPEC_REAL                   UINT32_C(0x00000004)  /**< need symmetries also for continuous variables */

typedef uint32_t SYM_SPEC;              /**< types of variables handled by symmetry */

/** symmetry timings */
#define SYM_COMPUTETIMING_BEFOREPRESOL    0  /**< compute symmetries before presolving */
#define SYM_COMPUTETIMING_DURINGPRESOL    1  /**< compute symmetries during presolving */
#define SYM_COMPUTETIMING_AFTERPRESOL     2  /**< compute symmetries after presolving */

/** define sense of rhs */
enum SYM_Rhssense
{
   SYM_SENSE_UNKOWN     = 0,                 /**< unknown sense */
   SYM_SENSE_INEQUALITY = 1,                 /**< linear inequality */
   SYM_SENSE_EQUATION   = 2,                 /**< linear equation */
   SYM_SENSE_XOR        = 3,                 /**< XOR constraint */
   SYM_SENSE_AND        = 4,                 /**< AND constraint */
   SYM_SENSE_OR         = 5,                 /**< OR constrant */
   SYM_SENSE_BOUNDIS_TYPE_1 = 6,             /**< bounddisjunction type 1 */
   SYM_SENSE_BOUNDIS_TYPE_2 = 7              /**< bounddisjunction type 2 */
};
typedef enum SYM_Rhssense SYM_RHSSENSE;

/* type of symmetry handling codes */
#define SYM_HANDLETYPE_NONE             UINT32_C(0x00000000)  /**< no symmetry handling */
#define SYM_HANDLETYPE_SYMBREAK         UINT32_C(0x00000001)  /**< symmetry breaking inequalities */
#define SYM_HANDLETYPE_ORBITALFIXING    UINT32_C(0x00000002)  /**< orbital fixing */
#define SYM_HANDLETYPE_SST              UINT32_C(0x00000004)  /**< Schreier Sims cuts */
#define SYM_HANDLETYPE_SYMCONS (SYM_HANDLETYPE_SYMBREAK | SYM_HANDLETYPE_SST)

typedef uint32_t SYM_HANDLETYPE;        /**< type of symmetry handling */

typedef struct SYM_Vartype SYM_VARTYPE;      /**< data of variables that are considered to be equivalent */
typedef struct SYM_Optype SYM_OPTYPE;        /**< data of operators that are considered to be equivalent */
typedef struct SYM_Consttype SYM_CONSTTYPE;  /**< data of constants that are considered to be equivalent */
typedef struct SYM_Rhstype SYM_RHSTYPE;      /**< data of constraint sides that are considered to be equivalent */
typedef struct SYM_Matrixdata SYM_MATRIXDATA;/**< data for symmetry group computation on linear constraints */
typedef struct SYM_Exprdata SYM_EXPRDATA;    /**< data for symmetry group computation on nonlinear constraints */

/** selection rules for leaders in SST cuts */
enum SCIP_LeaderRule
{
   SCIP_LEADERRULE_FIRSTINORBIT        = 0,       /**< first var in orbit */
   SCIP_LEADERRULE_LASTINORBIT         = 1,       /**< last var in orbit */
   SCIP_LEADERRULE_MAXCONFLICTSINORBIT = 2,       /**< var with most conflicting vars in its orbit */
   SCIP_LEADERRULE_MAXCONFLICTS        = 3        /**< var with most conflicting vars in problem */
};
typedef enum SCIP_LeaderRule SCIP_LEADERRULE;

/** tie breaks for leader rule based on the leader's orbit */
enum SCIP_LeaderTiebreakRule
{
   SCIP_LEADERTIEBREAKRULE_MINORBIT            = 0,    /**< orbit of minimum size */
   SCIP_LEADERTIEBREAKRULE_MAXORBIT            = 1,    /**< orbit of maximum size */
   SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT = 2     /**< orbit with maximum number of vars in conflict with leader */
};

/** variable types for leader in Schreier Sims cuts */
enum SCIP_SSTType
{
   SCIP_SSTTYPE_BINARY                 = 1,    /**< binary variables */
   SCIP_SSTTYPE_INTEGER                = 2,    /**< integer variables */
   SCIP_SSTTYPE_IMPLINT                = 4,    /**< implicitly integer variables */
   SCIP_SSTTYPE_CONTINUOUS             = 8     /**< continuous variables */
};

typedef enum SCIP_SSTType SCIP_SSTTYPE;

/** type of orbitope constraint: full, packing, or partitioning orbitope */
enum SCIP_OrbitopeType
{
   SCIP_ORBITOPETYPE_FULL         = 0,       /**< constraint is a full orbitope constraint:         rowsum(x) unrestricted */
   SCIP_ORBITOPETYPE_PARTITIONING = 1,       /**< constraint is a partitioning orbitope constraint: rowsum(x) == 1 */
   SCIP_ORBITOPETYPE_PACKING      = 2        /**< constraint is a packing orbitope constraint:      rowsum(x) <= 1 */
};
typedef enum SCIP_OrbitopeType SCIP_ORBITOPETYPE;

/** conditions to recompute symmetries after a restart */
enum SCIP_RecomputesymType
{
   SCIP_RECOMPUTESYM_NEVER         = 0,       /**< never recompute symmetries */
   SCIP_RECOMPUTESYM_ALWAYS        = 1,       /**< always recompute symmetries */
   SCIP_RECOMPUTESYM_OFFOUNDRED    = 2        /**< only if orbital fixing found a reduction in previous run */
};
typedef enum SCIP_RecomputesymType SCIP_RECOMPUTESYMTYPE;


#ifdef __cplusplus
}
#endif

#endif
