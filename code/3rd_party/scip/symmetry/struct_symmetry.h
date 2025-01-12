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

/**@file   struct_symmetry.h
 * @brief  structs for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SYMMETRY_H_
#define __SCIP_STRUCT_SYMMETRY_H_

#include "scip/scip.h"
#include "symmetry/type_symmetry.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data of variables that are considered to be equivalent */
struct SYM_Vartype
{
   SCIP_Real             obj;                /**< objective of variable */
   SCIP_Real             lb;                 /**< lower bound of variable */
   SCIP_Real             ub;                 /**< upper bound of variable */
   SCIP_VARTYPE          type;               /**< type of variable */
   int                   nconss;             /**< number of conss a variable is contained in */
   int                   color;              /**< store color */
};

/** data of operators that are considered to be equivalent */
struct SYM_Optype
{
   SCIP_EXPR*            expr;               /**< the underlying expression */
   int                   level;              /**< level of operator in its expression tree */
   int                   color;              /**< store color */
};

/** data of constants that are considered to be equivalent */
struct SYM_Consttype
{
   SCIP_Real             value;              /**< value of constant */
   int                   color;              /**< store color */
};

/** data of coefficients that are considered to be equivalent */
struct SYM_Rhstype
{
   SCIP_Real             lhs;                /**< value of left-hand-side */
   SCIP_Real             rhs;                /**< value of right-hand-side */
   int                   color;              /**< store color */
};

/** data for symmetry group computation on linear constraints */
struct SYM_Matrixdata
{
   SCIP_Real*            matcoef;            /**< nonzero coefficients appearing in the matrix */
   SCIP_Real*            rhscoef;            /**< rhs coefficients */
   SYM_RHSSENSE*         rhssense;           /**< sense of rhs */
   int*                  matrhsidx;          /**< indices of rhs corresponding to matrix entries */
   int*                  matvaridx;          /**< indices of variables for matrix entries */
   int*                  matidx;             /**< indices in mat(rhs/var)idx array corresponding to matrix coefficients */
   int*                  rhsidx;             /**< indices in rhstype array corresponding to rhs coefficients */
   int*                  permvarcolors;      /**< array for storing the colors of the individual variables */
   int*                  matcoefcolors;      /**< array for storing the colors of all matrix coefficients */
   int*                  rhscoefcolors;      /**< array for storing the colors of all rhs coefficients */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nmatcoef;           /**< number of coefficients in matrix */
   int                   nrhscoef;           /**< number of coefficients in rhs */
   int                   nmaxmatcoef;        /**< maximal number of matrix coefficients (will be increase on demand) */
   int                   nuniquevars;        /**< number of unique variable types */
   int                   nuniquerhs;         /**< number of unique rhs types */
   int                   nuniquemat;         /**< number of unique matrix coefficients */
};

/** data for symmetry group computation on nonlinear constraints */
struct SYM_Exprdata
{
   int                   nuniqueconstants;   /**< number of unique constants */
   int                   nuniqueoperators;   /**< number of unique operators */
   int                   nuniquecoefs;       /**< number of unique coefficients */
};

#ifdef __cplusplus
}
#endif

#endif
