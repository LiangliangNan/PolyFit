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

/**@file   struct_conflict.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONFLICT_H__
#define __SCIP_STRUCT_CONFLICT_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_misc.h"
#include "scip/type_var.h"
#include "scip/type_conflict.h"
#include "lpi/type_lpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** conflict handler */
struct SCIP_Conflicthdlr
{
   char*                 name;               /**< name of conflict handler */
   char*                 desc;               /**< description of conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy));  /**< copy method of conflict handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONFLICTFREE((*conflictfree));  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit));  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit));  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol));/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol));/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec));  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata;  /**< conflict handler data */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this conflict handler for the next stages */
   SCIP_CLOCK*           conflicttime;       /**< conflict handler execution time */
   int                   priority;           /**< priority of the conflict handler */
   SCIP_Bool             initialized;        /**< is conflict handler initialized? */
};

/** set of conflicting bound changes */
struct SCIP_ConflictSet
{
   SCIP_BDCHGINFO**      bdchginfos;         /**< bound change informations of the conflict set */
   SCIP_BDCHGINFO*       confbdchginfo;      /**< a bound change at the conflict depth */
   SCIP_Real*            relaxedbds;         /**< array of relaxed bounds which are efficient for a valid conflict */
   SCIP_Real             confrelaxedbd;      /**< relaxed bound belonging the the bound change at the conflict depth */
   int*                  sortvals;           /**< aggregated var index/bound type values for sorting */
   int                   bdchginfossize;     /**< size of bdchginfos array */
   int                   nbdchginfos;        /**< number of bound change informations in the conflict set */
   int                   validdepth;         /**< depth in the tree where the conflict set is valid */
   int                   insertdepth;        /**< depth level where constraint should be added */
   int                   conflictdepth;      /**< depth in the tree where the conflict set yields a conflict */
   int                   repropdepth;        /**< depth at which the conflict set triggers a deduction */
   unsigned int          repropagate:1;      /**< should the conflict constraint trigger a repropagation? */
   unsigned int          depthcalced:1;      /**< are the conflict and repropagation depth calculated? */
   unsigned int          sorted:1;           /**< is the conflict set sorted */
   unsigned int          usescutoffbound:1;  /**< is the conflict based on the cutoff bound? */
   unsigned int          hasrelaxonlyvar:1;  /**< is one of the bound change informations using a relaxation-only variable */
   SCIP_CONFTYPE         conflicttype;       /**< conflict type: unknown, infeasible LP, bound exceeding LP, propagation */
};

/** set of conflicting bound changes */
struct SCIP_ProofSet
{
   SCIP_Real*            vals;
   int*                  inds;
   SCIP_Real             rhs;
   int                   nnz;
   int                   size;
   int                   validdepth;
   SCIP_CONFTYPE         conflicttype;       /**< conflict type: unknown, infeasible LP, bound exceeding LP */
};

/** set of LP bound change */
struct SCIP_LPBdChgs
{
   int*                  bdchginds;          /**< array of column indices */
   SCIP_Real*            bdchglbs;           /**< array of lower bounds */
   SCIP_Real*            bdchgubs;           /**< array of upper bounds */
   int*                  bdchgcolinds;       /**< array of ???????????? */
   SCIP_Bool*            usedcols;           /**< array to mark if a column is used */
   int                   nbdchgs;            /**< number of stored LP bound changes */
};

/** conflict analysis data structure */
struct SCIP_Conflict
{
   SCIP_Longint          nglbchgbds;         /**< total number of applied global bound changes */
   SCIP_Longint          nappliedglbconss;   /**< total number of conflict constraints added globally to the problem */
   SCIP_Longint          nappliedglbliterals;/**< total number of literals in globally applied conflict constraints */
   SCIP_Longint          nlocchgbds;         /**< total number of applied local bound changes */
   SCIP_Longint          nappliedlocconss;   /**< total number of conflict constraints added locally to the problem */
   SCIP_Longint          nappliedlocliterals;/**< total number of literals in locally applied conflict constraints */
   SCIP_Longint          npropcalls;         /**< number of calls to propagation conflict analysis */
   SCIP_Longint          npropsuccess;       /**< number of calls yielding at least one conflict constraint */
   SCIP_Longint          npropconfconss;     /**< number of valid conflict constraints detected in propagation conflict analysis */
   SCIP_Longint          npropconfliterals;  /**< total number of literals in valid propagation conflict constraints */
   SCIP_Longint          npropreconvconss;   /**< number of reconvergence constraints detected in propagation conflict analysis */
   SCIP_Longint          npropreconvliterals;/**< total number of literals in valid propagation reconvergence constraints */
   SCIP_Longint          ninflpcalls;        /**< number of calls to infeasible LP conflict analysis */
   SCIP_Longint          ninflpsuccess;      /**< number of calls yielding at least one conflict constraint */
   SCIP_Longint          ninflpconfconss;    /**< number of valid conflict constraints detected in infeasible LP conflict
                                              *   analysis */
   SCIP_Longint          ninflpconfliterals; /**< total number of literals in valid infeasible LP conflict constraints */
   SCIP_Longint          ninflpreconvconss;  /**< number of reconvergence constraints detected in infeasible LP conflict
                                              *   analysis */
   SCIP_Longint          ninflpreconvliterals; /**< total number of literals in valid infeasible LP reconvergence
                                                *   constraints */
   SCIP_Longint          ninflpiterations;   /**< total number of LP iterations used in infeasible LP conflict analysis */
   SCIP_Longint          nboundlpcalls;      /**< number of calls to bound exceeding LP conflict analysis */
   SCIP_Longint          nboundlpsuccess;    /**< number of calls yielding at least one conflict constraint */
   SCIP_Longint          nboundlpconfconss;  /**< number of valid conflict constraints detected in bound exceeding LP
                                              *   conflict analysis */
   SCIP_Longint          nboundlpconfliterals; /**< total number of literals in valid bound exceeding LP conflict
                                                *   constraints */
   SCIP_Longint          nboundlpreconvconss;/**< number of reconvergence constraints detected in bound exceeding LP
                                              *   conflict analysis */
   SCIP_Longint          nboundlpreconvliterals; /**< total number of literals in valid bound exceeding LP reconvergence
                                                  *   constraints */
   SCIP_Longint          nboundlpiterations; /**< total number of LP iterations used in bound exceeding LP conflict
                                              *   analysis */
   SCIP_Longint          nsbcalls;           /**< number of calls to infeasible strong branching conflict analysis */
   SCIP_Longint          nsbsuccess;         /**< number of calls yielding at least one conflict constraint */
   SCIP_Longint          nsbconfconss;       /**< number of conflict constraints detected in strong branching conflict analysis */
   SCIP_Longint          nsbconfliterals;    /**< total number of literals in valid strong branching conflict constraints */
   SCIP_Longint          nsbreconvconss;     /**< number of reconvergence constraints detected in strong branch conflict analysis */
   SCIP_Longint          nsbreconvliterals;  /**< total number of literals in valid strong branching reconvergence constraints */
   SCIP_Longint          nsbiterations;      /**< total number of LP iterations used in strong branching conflict analysis */
   SCIP_Longint          npseudocalls;       /**< number of calls to pseudo solution conflict analysis */
   SCIP_Longint          npseudosuccess;     /**< number of calls yielding at least one conflict constraint */
   SCIP_Longint          npseudoconfconss;   /**< number of valid conflict constraints detected in pseudo sol conflict analysis */
   SCIP_Longint          npseudoconfliterals;/**< total number of literals in valid pseudo solution conflict constraints */
   SCIP_Longint          npseudoreconvconss; /**< number of reconvergence constraints detected in pseudo sol conflict analysis */
   SCIP_Longint          npseudoreconvliterals;/**< total number of literals in valid pseudo solution reconvergence constraints */
   SCIP_Longint          ndualproofsinfglobal;/**< number of globally added dual proof constraints derived from infeasible LP */
   SCIP_Longint          ndualproofsinflocal;/**< number of locally added dual proof constraints derived from infeasible LP */
   SCIP_Longint          ndualproofsinfsuccess;/**< number of successfully dual proof analysis calls for infeasible LPs */
   SCIP_Longint          dualproofsinfnnonzeros;/**< number of non-zeros over all accepted dual proof constraints derived from infeasible LP */
   SCIP_Longint          ndualproofsbndglobal;/**< number of globally added dual proof constraints derived from bound exceeding LP */
   SCIP_Longint          ndualproofsbndlocal;/**< number of locally added dual proof constraints derived from bound exceeding LP */
   SCIP_Longint          ndualproofsbndsuccess;/**< number of successfully dual proof analysis calls for bound exceeding LPs */
   SCIP_Longint          dualproofsbndnnonzeros;/**< number of non-zeros over all accepted dual proof constraints derived from bound exceeding LPs */

   SCIP_CLOCK*           dIBclock;           /**< time used for detect implied bounds */

   SCIP_CLOCK*           propanalyzetime;    /**< time used for propagation conflict analysis */
   SCIP_CLOCK*           inflpanalyzetime;   /**< time used for infeasible LP conflict analysis */
   SCIP_CLOCK*           boundlpanalyzetime; /**< time used for bound exceeding LP conflict analysis */
   SCIP_CLOCK*           sbanalyzetime;      /**< time used for strong branching LP conflict analysis */
   SCIP_CLOCK*           pseudoanalyzetime;  /**< time used for pseudo solution conflict analysis */
   SCIP_PQUEUE*          bdchgqueue;         /**< unprocessed conflict bound changes */
   SCIP_PQUEUE*          forcedbdchgqueue;   /**< unprocessed conflict bound changes that must be resolved */
   SCIP_PROOFSET*        proofset;           /**< proof sets found at the current node */
   SCIP_PROOFSET**       proofsets;          /**< proof sets found at the current node */
   SCIP_CONFLICTSET*     conflictset;        /**< bound changes resembling the current conflict set */
   SCIP_CONFLICTSET**    conflictsets;       /**< conflict sets found at the current node */
   SCIP_Real*            conflictsetscores;  /**< score values of the conflict sets found at the current node */
   SCIP_BDCHGINFO**      tmpbdchginfos;      /**< temporarily created bound change information data */
   int                   conflictsetssize;   /**< size of conflictsets array */
   int                   nconflictsets;      /**< number of available conflict sets (used slots in conflictsets array) */
   int                   proofsetssize;      /**< size of proofsets array */
   int                   nproofsets;         /**< number of available proof sets (used slots in proofsets array) */
   int                   tmpbdchginfossize;  /**< size of tmpbdchginfos array */
   int                   ntmpbdchginfos;     /**< number of temporary created bound change information data */
   int                   count;              /**< conflict set counter to label binary conflict variables with */
};

#ifdef __cplusplus
}
#endif

#endif
