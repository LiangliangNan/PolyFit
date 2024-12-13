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

/**@file   struct_sepa.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SEPA_H__
#define __SCIP_STRUCT_SEPA_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_sepa.h"

#ifdef __cplusplus
extern "C" {
#endif

/** separators data */
struct SCIP_Sepa
{
   SCIP_Longint          lastsepanode;       /**< last (total) node where this separator was called */
   SCIP_Longint          ncalls;             /**< number of times, this separator was called */
   SCIP_Longint          nrootcalls;         /**< number of times, this separator was called at the root */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this separator */
   SCIP_Longint          ncutsfound;         /**< number of cutting planes found so far by this separator */
   SCIP_Longint          ncutsadded;         /**< number of cutting planes added to sepastore equal to
                                              *   the sum of added cuts via pool and direct.*/
   SCIP_Longint          ncutsaddedviapool;  /**< number of cutting planes added from cutpool */
   SCIP_Longint          ncutsaddeddirect;   /**< number of cutting planes added directly */
   SCIP_Longint          ncutsappliedviapool;/**< number of cutting planes applied to LP via cutpool */
   SCIP_Longint          ncutsapplieddirect; /**< number of cutting planes applied to LP directly from sepastore */
   SCIP_Longint          nconssfound;        /**< number of additional constraints added by this separator */
   SCIP_Longint          ndomredsfound;      /**< number of domain reductions found so far by this separator */
   SCIP_Real             maxbounddist;       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   char*                 name;               /**< name of separator */
   char*                 desc;               /**< description of separator */
   SCIP_DECL_SEPACOPY    ((*sepacopy));      /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_SEPAFREE    ((*sepafree));      /**< destructor of separator */
   SCIP_DECL_SEPAINIT    ((*sepainit));      /**< initialize separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit));      /**< deinitialize separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol));   /**< solving process initialization method of separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol));   /**< solving process deinitialization method of separator */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp));    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol));   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata;           /**< separators local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this separator for the next stages */
   SCIP_CLOCK*           sepaclock;          /**< separation time */
   int                   priority;           /**< priority of the separator */
   int                   freq;               /**< frequency for calling separator */
   int                   ncallsatnode;       /**< number of times, this separator was called at the current node */
   int                   ncutsfoundatnode;   /**< number of cutting planes found at the current node */
   int                   expbackoff;         /**< base for exponential increase of frequency at which the separator is called */
   SCIP_Bool             usessubscip;        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay;              /**< should separator be delayed, if other separators found cuts? */
   SCIP_Bool             lpwasdelayed;       /**< was the LP separation delayed at the last call? */
   SCIP_Bool             solwasdelayed;      /**< was the solution separation delayed at the last call? */
   SCIP_Bool             initialized;        /**< is separator initialized? */
   SCIP_Bool             isparentsepa;       /**< is separator a parent separator that create cuts of child separators? */
   struct SCIP_Sepa*     parentsepa;         /**< pointer to parent separator or NULL */
};

#ifdef __cplusplus
}
#endif

#endif
