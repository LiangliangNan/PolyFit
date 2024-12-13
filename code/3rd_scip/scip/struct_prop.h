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

/**@file   struct_prop.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for propagators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_PROP_H__
#define __SCIP_STRUCT_PROP_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_prop.h"

#ifdef __cplusplus
extern "C" {
#endif

/** propagators data */
struct SCIP_Prop
{
   SCIP_Longint          ncalls;             /**< number of times, this propagator was called */
   SCIP_Longint          nrespropcalls;      /**< number of times, the resolve propagation was called */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this propagator */
   SCIP_Longint          ndomredsfound;      /**< number of domain reductions found so far by this propagator */
   char*                 name;               /**< name of propagator */
   char*                 desc;               /**< description of propagator */
   SCIP_DECL_PROPCOPY    ((*propcopy));      /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PROPFREE    ((*propfree));      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit));      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit));      /**< deinitialize propagator */
   SCIP_DECL_PROPINITPRE ((*propinitpre));   /**< presolving initialization method of propagator */
   SCIP_DECL_PROPEXITPRE ((*propexitpre));   /**< presolving deinitialization method of propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol));   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol));   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol));    /**< presolving method of propagator */
   SCIP_DECL_PROPEXEC    ((*propexec));      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop));   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata;           /**< propagators local data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this propagator for the next stages */
   SCIP_CLOCK*           proptime;           /**< time used for propagation of this propagator */
   SCIP_CLOCK*           sbproptime;         /**< time used for propagation of this propagator during strong branching */
   SCIP_CLOCK*           resproptime;        /**< time used for resolve propagation of this propagator */
   SCIP_CLOCK*           presoltime;         /**< time used for presolving of this propagator */
   int                   priority;           /**< priority of the propagator for propagation */
   int                   freq;               /**< frequency for calling propagator */
   SCIP_PROPTIMING       timingmask;         /**< positions in the node solving loop where propagator should be executed */
   SCIP_PRESOLTIMING     presoltiming;       /**< timing mask of the presolving method of the propagator */
   int                   presolpriority;     /**< priority of the presolving of the propagator */
   int                   maxprerounds;       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   int                   lastnfixedvars;     /**< number of variables fixed before the last call to the propagator */
   int                   lastnaggrvars;      /**< number of variables aggregated in presolving before the last call to the propagator */
   int                   lastnchgvartypes;   /**< number of variable type changes in presolving before the last call to the propagator */
   int                   lastnchgbds;        /**< number of variable bounds tightened in presolving before the last call to the propagator */
   int                   lastnaddholes;      /**< number of domain holes added in presolving before the last call to the propagator */
   int                   lastndelconss;      /**< number of deleted constraints in presolving before the last call to the propagator */
   int                   lastnaddconss;      /**< number of added constraints in presolving before the last call to the propagator */
   int                   lastnupgdconss;     /**< number of upgraded constraints in presolving before the last call to the propagator */
   int                   lastnchgcoefs;      /**< number of changed coefficients in presolving before the last call to the propagator */
   int                   lastnchgsides;      /**< number of changed left or right hand sides in presolving before the last call to the propagator */
   int                   nfixedvars;         /**< total number of variables fixed by this propagator in presolving */
   int                   naggrvars;          /**< total number of variables aggregated by this propagator in presolving */
   int                   nchgvartypes;       /**< total number of variable type changes by this propagator in presolving */
   int                   nchgbds;            /**< total number of variable bounds tightened by this propagator in presolving */
   int                   naddholes;          /**< total number of domain holes added by this propagator in presolving */
   int                   ndelconss;          /**< total number of deleted constraints by this propagator in presolving */
   int                   naddconss;          /**< total number of added constraints by this propagator in presolving */
   int                   nupgdconss;         /**< total number of upgraded constraints by this propagator in presolving */
   int                   nchgcoefs;          /**< total number of changed coefficients by this propagator in presolving */
   int                   nchgsides;          /**< total number of changed left or right hand sides by this propagator in presolving */
   int                   npresolcalls;       /**< number of times the propagator was called in presolving and tried to find reductions */
   SCIP_Bool             delay;              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_Bool             wasdelayed;         /**< was the propagator delayed at the last call? */
   SCIP_Bool             initialized;        /**< is propagator initialized? */
};

#ifdef __cplusplus
}
#endif

#endif
