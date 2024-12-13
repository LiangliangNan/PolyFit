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

/**@file   struct_prob.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PROB_H__
#define __SCIP_STRUCT_PROB_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_cons.h"

#ifdef __cplusplus
extern "C" {
#endif

/** main problem to solve */
struct SCIP_Prob
{
   SCIP_Real             objoffset;          /**< objective offset from bound shifting and fixing (fixed vars result) */
   SCIP_Real             objscale;           /**< scalar applied to objective function; external objective value is
                                              *   extobj = objsense * objscale * (intobj + objoffset) */
   SCIP_Real             objlim;             /**< objective limit as external value (original problem space) */
   SCIP_Real             dualbound;          /**< dual bound as external value (original problem space) which is given or update during presolving */
   char*                 name;               /**< problem name */
   SCIP_DECL_PROBCOPY    ((*probcopy));      /**< copies user data if you want to copy it to a subscip, or NULL */
   SCIP_DECL_PROBDELORIG ((*probdelorig));   /**< frees user data of original problem */
   SCIP_DECL_PROBTRANS   ((*probtrans));     /**< creates user data of transformed problem by transforming original user data */
   SCIP_DECL_PROBDELTRANS((*probdeltrans));  /**< frees user data of transformed problem */
   SCIP_DECL_PROBINITSOL ((*probinitsol));   /**< solving process initialization method of transformed data */
   SCIP_DECL_PROBEXITSOL ((*probexitsol));   /**< solving process deinitialization method of transformed data */
   SCIP_PROBDATA*        probdata;           /**< user problem data set by the reader */
   SCIP_HASHTABLE*       varnames;           /**< hash table storing variable's names */
   SCIP_VAR**            vars;               /**< array with active variables ordered binary, integer, implicit, continuous */
   SCIP_VAR**            fixedvars;          /**< array with fixed and aggregated variables */
   SCIP_VAR**            deletedvars;        /**< array to temporarily store deleted variables */
   SCIP_HASHTABLE*       consnames;          /**< hash table storing constraints' names */
   SCIP_CONS**           conss;              /**< array with constraints of the problem */
   int                   varssize;           /**< available slots in vars array */
   int                   nvars;              /**< number of active variables in the problem (used slots in vars array) */
   int                   nbinvars;           /**< number of binary variables */
   int                   nintvars;           /**< number of general integer variables */
   int                   nimplvars;          /**< number of implicit integer variables */
   int                   ncontvars;          /**< number of continuous variables */
   int                   ncolvars;           /**< number of variables with attached column information */
   int                   fixedvarssize;      /**< available slots in fixedvars array */
   int                   nfixedvars;         /**< number of fixed and aggregated variables in the problem */
   int                   deletedvarssize;    /**< available slots in deletedvars array */
   int                   ndeletedvars;       /**< number of deleted variables in the problem */
   int                   nobjvars;           /**< number of variables with a non-zero objective coefficient */
   int                   consssize;          /**< available slots in conss array */
   int                   nconss;             /**< number of constraints in the problem (number of used slots in conss array) */
   int                   maxnconss;          /**< maximum number of constraints existing at the same time */
   int                   startnvars;         /**< number of variables existing when problem solving started */
   int                   startnconss;        /**< number of constraints existing when problem solving started */
   SCIP_OBJSENSE         objsense;           /**< objective sense of the original problem */
   SCIP_Bool             objisintegral;      /**< is objective value always integral for feasible solutions? */
   SCIP_Bool             transformed;        /**< TRUE iff problem is the transformed problem */
   SCIP_Bool             nlpenabled;         /**< marks whether an NLP relaxation should be constructed */
   SCIP_Bool             permuted;           /**< TRUE iff the problem is already permuted */
   SCIP_Bool             conscompression;    /**< TRUE for problems for which constraint compression on a set of fixed variables is desired */
};

#ifdef __cplusplus
}
#endif

#endif
