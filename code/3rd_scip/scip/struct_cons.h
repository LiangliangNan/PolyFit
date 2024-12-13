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

/**@file   struct_cons.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for constraints and constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONS_H__
#define __SCIP_STRUCT_CONS_H__

#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_cons.h"


#ifdef __cplusplus
extern "C" {
#endif

/** constraint data structure */
struct SCIP_Cons
{
   SCIP_Real             age;                /**< age of constraint: number of successive times, the constraint was irrelevant */
   char*                 name;               /**< name of the constraint */
   SCIP_CONSHDLR*        conshdlr;           /**< constraint handler for this constraint */
   SCIP_CONSDATA*        consdata;           /**< data for this specific constraint */
   SCIP_CONS*            transorigcons;      /**< for original constraints: associated transformed constraint or NULL,
                                              *   for transformed constraints: associated original constraint or NULL */
   SCIP_CONSSETCHG*      addconssetchg;      /**< constraint change that added constraint to current subproblem, or NULL if
                                              *   constraint is from global problem */
   int                   addarraypos;        /**< position of constraint in the conssetchg's/prob's addedconss/conss array */
   int                   consspos;           /**< position of constraint in the handler's conss array */
   int                   initconsspos;       /**< position of constraint in the handler's initconss array */
   int                   sepaconsspos;       /**< position of constraint in the handler's sepaconss array */
   int                   enfoconsspos;       /**< position of constraint in the handler's enfoconss array */
   int                   checkconsspos;      /**< position of constraint in the handler's checkconss array */
   int                   propconsspos;       /**< position of constraint in the handler's propconss array */
   int                   nlockspos[NLOCKTYPES]; /**< array of times, the constraint locked rounding of its variables */
   int                   nlocksneg[NLOCKTYPES]; /**< array of times, the constraint locked vars for the constraint's negation */
   int                   activedepth;        /**< depth level of constraint activation (-2: inactive, -1: problem constraint) */
   int                   validdepth;         /**< depth level where constraint is valid (-1: equals activedepth) */
   int                   nuses;              /**< number of times, this constraint is referenced */
   unsigned int          initial:1;          /**< TRUE iff LP relaxation of constraint should be in initial LP, if possible */
   unsigned int          separate:1;         /**< TRUE iff constraint should be separated during LP processing */
   unsigned int          enforce:1;          /**< TRUE iff constraint should be enforced during node processing */
   unsigned int          check:1;            /**< TRUE iff constraint should be checked for feasibility */
   unsigned int          propagate:1;        /**< TRUE iff constraint should be propagated during node processing */
   unsigned int          sepaenabled:1;      /**< TRUE iff constraint should be separated in the next separation call */
   unsigned int          propenabled:1;      /**< TRUE iff constraint should be propagated in the next propagation call */
   unsigned int          local:1;            /**< TRUE iff constraint is only valid locally */
   unsigned int          modifiable:1;       /**< TRUE iff constraint is modifiable (subject to column generation) */
   unsigned int          dynamic:1;          /**< TRUE iff constraint is subject to aging */
   unsigned int          removable:1;        /**< TRUE iff relaxation should be removed from the LP due to aging or cleanup */
   unsigned int          stickingatnode:1;   /**< TRUE iff the node should always be kept at the node where it was added */
   unsigned int          original:1;         /**< TRUE iff constraint belongs to original problem */
   unsigned int          deleteconsdata:1;   /**< TRUE iff constraint data has to be deleted if constraint is freed */
   unsigned int          active:1;           /**< TRUE iff constraint is active in the current node; a constraint is
                                              *   active if it is global and was not removed during presolving or it was
                                              *   added locally (in that case the local flag is TRUE) and the current
                                              *   node belongs to the corresponding sub tree
                                              */
   unsigned int          conflict:1;         /**< TRUE iff constraint is a conflict */
   unsigned int          enabled:1;          /**< TRUE iff constraint is enforced, separated, and propagated in current node */
   unsigned int          obsolete:1;         /**< TRUE iff constraint is too seldomly used and therefore obsolete */
   unsigned int          markpropagate:1;    /**< TRUE iff constraint is marked to be propagated in the next round */
   unsigned int          deleted:1;          /**< TRUE iff constraint was globally deleted */
   unsigned int          update:1;           /**< TRUE iff constraint has to be updated in update phase */
   unsigned int          updateinsert:1;     /**< TRUE iff constraint has to be inserted in the conss array */
   unsigned int          updateactivate:1;   /**< TRUE iff constraint has to be activated in update phase */
   unsigned int          updatedeactivate:1; /**< TRUE iff constraint has to be deactivated in update phase */
   unsigned int          updateenable:1;     /**< TRUE iff constraint has to be enabled in update phase */
   unsigned int          updatedisable:1;    /**< TRUE iff constraint has to be disabled in update phase */
   unsigned int          updatesepaenable:1; /**< TRUE iff constraint's separation has to be enabled in update phase */
   unsigned int          updatesepadisable:1;/**< TRUE iff constraint's separation has to be disabled in update phase */
   unsigned int          updatepropenable:1; /**< TRUE iff constraint's propagation has to be enabled in update phase */
   unsigned int          updatepropdisable:1;/**< TRUE iff constraint's propagation has to be disabled in update phase */
   unsigned int          updateobsolete:1;   /**< TRUE iff obsolete status of constraint has to be updated in update phase */
   unsigned int          updatefree:1;       /**< TRUE iff constraint has to be freed in update phase */
   unsigned int          updateactfocus:1;   /**< TRUE iff delayed constraint activation happened at focus node */
   unsigned int          updatemarkpropagate:1;/**< TRUE iff constraint has to be marked to be propagated in update phase */
   unsigned int          updateunmarkpropagate:1;/**< TRUE iff constraint has to be unmarked to be propagated in update phase */
   unsigned int          nupgradelocks:28;   /**< number of times, a constraint is locked against an upgrade
                                              *   (e.g. linear -> logicor), 0 means a constraint can be upgraded */
#ifndef NDEBUG
   SCIP*                 scip;               /**< SCIP data structure */
#endif
};

/** tracks additions and removals of the set of active constraints */
struct SCIP_ConsSetChg
{
   SCIP_CONS**           addedconss;         /**< constraints added to the set of active constraints */
   SCIP_CONS**           disabledconss;      /**< constraints disabled in the set of active constraints */
   int                   addedconsssize;     /**< size of added constraints array */
   int                   naddedconss;        /**< number of added constraints */
   int                   disabledconsssize;  /**< size of disabled constraints array */
   int                   ndisabledconss;     /**< number of disabled constraints */
};

/** constraint handler */
struct SCIP_Conshdlr
{
   SCIP_Longint          nsepacalls;         /**< number of times, the separator was called */
   SCIP_Longint          nenfolpcalls;       /**< number of times, the LP enforcer was called */
   SCIP_Longint          nenfopscalls;       /**< number of times, the pseudo enforcer was called */
   SCIP_Longint          nenforelaxcalls;    /**< number of times, the relaxation enforcer was called */
   SCIP_Longint          npropcalls;         /**< number of times, the propagator was called */
   SCIP_Longint          ncheckcalls;        /**< number of times, the feasibility check was called */
   SCIP_Longint          nrespropcalls;      /**< number of times, the resolve propagation was called */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this constraint handler */
   SCIP_Longint          ncutsfound;         /**< number of cuts found by this constraint handler */
   SCIP_Longint          ncutsapplied;       /**< number of cuts found by this constraint handler applied to lp */
   SCIP_Longint          nconssfound;        /**< number of additional constraints added by this constraint handler */
   SCIP_Longint          ndomredsfound;      /**< number of domain reductions found so far by this constraint handler */
   SCIP_Longint          nchildren;          /**< number of children the constraint handler created during branching */
   SCIP_Longint          lastpropdomchgcount;/**< last bound change number, where the domain propagation was called */
   SCIP_Longint          storedpropdomchgcount;/**< bound change number, where the domain propagation was called last before starting probing */
   SCIP_Longint          lastenfolpdomchgcount;/**< last bound change number, where the LP enforcement was called */
   SCIP_Longint          lastenfopsdomchgcount;/**< last bound change number, where the pseudo enforcement was called */
   SCIP_Longint          lastenforelaxdomchgcount;/**< last bound change number, where the relaxation enforcement was called */
   SCIP_Longint          lastenfolpnode;     /**< last node at which the LP enforcement was called */
   SCIP_Longint          lastenfopsnode;     /**< last node at which the pseudo enforcement was called */
   SCIP_Longint          lastenforelaxnode;  /**< last node at which the relaxation enforcement was called */
   SCIP_RESULT           lastenfolpresult;   /**< result of last LP enforcement call */
   SCIP_RESULT           lastenfopsresult;   /**< result of last pseudo enforcement call */
   SCIP_RESULT           lastenforelaxresult;/**< result of last relaxation enforcement call */
   SCIP_Real             ageresetavg;        /**< exp. decaying weighted average of constraint ages at moment of age reset */
   char*                 name;               /**< name of constraint handler */
   char*                 desc;               /**< description of constraint handler */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy));  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSFREE    ((*consfree));      /**< destructor of constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit));      /**< initialize constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit));      /**< deinitialize constraint handler */
   SCIP_DECL_CONSINITPRE ((*consinitpre));   /**< presolving initialization method of constraint handler */
   SCIP_DECL_CONSEXITPRE ((*consexitpre));   /**< presolving deinitialization method of constraint handler */
   SCIP_DECL_CONSINITSOL ((*consinitsol));   /**< solving process initialization method of constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol));   /**< solving process deinitialization method of constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete));    /**< free specific constraint data */
   SCIP_DECL_CONSTRANS   ((*constrans));     /**< transform constraint data into data belonging to the transformed problem */
   SCIP_DECL_CONSINITLP  ((*consinitlp));    /**< initialize LP with relaxations of "initial" constraints */
   SCIP_DECL_CONSSEPALP  ((*conssepalp));    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol));   /**< separate cutting planes for arbitrary primal solution */
   SCIP_DECL_CONSENFOLP  ((*consenfolp));    /**< enforcing constraints for LP solutions */
   SCIP_DECL_CONSENFORELAX ((*consenforelax)); /**< enforcing constraints for relaxation solutions */
   SCIP_DECL_CONSENFOPS  ((*consenfops));    /**< enforcing constraints for pseudo solutions */
   SCIP_DECL_CONSCHECK   ((*conscheck));     /**< check feasibility of primal solution */
   SCIP_DECL_CONSPROP    ((*consprop));      /**< propagate variable domains */
   SCIP_DECL_CONSPRESOL  ((*conspresol));    /**< presolving method */
   SCIP_DECL_CONSRESPROP ((*consresprop));   /**< propagation conflict resolving method */
   SCIP_DECL_CONSLOCK    ((*conslock));      /**< variable rounding lock method */
   SCIP_DECL_CONSACTIVE  ((*consactive));    /**< activation notification method */
   SCIP_DECL_CONSDEACTIVE((*consdeactive));  /**< deactivation notification method */
   SCIP_DECL_CONSENABLE  ((*consenable));    /**< enabling notification method */
   SCIP_DECL_CONSDISABLE ((*consdisable));   /**< disabling notification method */
   SCIP_DECL_CONSDELVARS ((*consdelvars));   /**< variable deletion method */
   SCIP_DECL_CONSPRINT   ((*consprint));     /**< constraint display method */
   SCIP_DECL_CONSCOPY    ((*conscopy));      /**< constraint copying method */
   SCIP_DECL_CONSPARSE   ((*consparse));     /**< constraint parsing method */
   SCIP_DECL_CONSGETVARS ((*consgetvars));   /**< constraint get variables method */
   SCIP_DECL_CONSGETNVARS((*consgetnvars));  /**< constraint get number of variable method */
   SCIP_DECL_CONSGETDIVEBDCHGS((*consgetdivebdchgs)); /**< constraint handler diving solution enforcement method */
   SCIP_CONSHDLRDATA*    conshdlrdata;       /**< constraint handler data */
   SCIP_CONS**           conss;              /**< array with all transformed constraints, active ones preceed inactive
                                              *   ones; a constraint is active if it is global and was not removed
                                              *   during presolving or it was added locally (in that case the local flag
                                              *   is TRUE) and the current node belongs to the corresponding sub tree */
   SCIP_CONS**           initconss;          /**< array with active constraints that must enter the LP with their initial representation */
   SCIP_CONS**           sepaconss;          /**< array with active constraints that must be separated during LP processing */
   SCIP_CONS**           enfoconss;          /**< array with active constraints that must be enforced during node processing */
   SCIP_CONS**           checkconss;         /**< array with active constraints that must be checked for feasibility */
   SCIP_CONS**           propconss;          /**< array with active constraints that must be propagated during node processing */
   SCIP_CONS**           storedpropconss;    /**< array to store constraints that were marked for propagation before
                                              *   starting probing mode
                                              */
   SCIP_CONS**           updateconss;        /**< array with constraints that changed and have to be update in the handler */
   SCIP_CLOCK*           setuptime;          /**< time spend for setting up this constraint handler for the next stages */
   SCIP_CLOCK*           presoltime;         /**< time used for presolving of this constraint handler */
   SCIP_CLOCK*           sepatime;           /**< time used for separation of this constraint handler */
   SCIP_CLOCK*           enfolptime;         /**< time used for LP enforcement of this constraint handler */
   SCIP_CLOCK*           enfopstime;         /**< time used for pseudo enforcement of this constraint handler */
   SCIP_CLOCK*           enforelaxtime;      /**< time used for relaxation enforcement of this constraint handler */
   SCIP_CLOCK*           proptime;           /**< time used for propagation of this constraint handler */
   SCIP_CLOCK*           sbproptime;         /**< time used for propagation of this constraint handler during strong branching */
   SCIP_CLOCK*           checktime;          /**< time used for feasibility check of this constraint handler */
   SCIP_CLOCK*           resproptime;        /**< time used for resolve propagation of this constraint handler */
   SCIP_Longint          lastsepalpcount;    /**< last LP number, where the separations was called */
   SCIP_Longint          lastenfolplpcount;  /**< last LP number, where the LP enforcement was called */
   SCIP_Longint          lastenforelaxrelaxcount; /**< last relax number, where the relax enforcement was called */
   int                   sepapriority;       /**< priority of the constraint handler for separation */
   int                   enfopriority;       /**< priority of the constraint handler for constraint enforcing */
   int                   checkpriority;      /**< priority of the constraint handler for checking infeasibility */
   int                   sepafreq;           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   propfreq;           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int                   eagerfreq;          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int                   maxprerounds;       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   int                   consssize;          /**< size of conss array */
   int                   nconss;             /**< total number of constraints */
   int                   nactiveconss;       /**< total number of active constraints */
   int                   maxnactiveconss;    /**< maximal number of active constraints existing at the same time */
   int                   startnactiveconss;  /**< number of active constraints existing when problem solving started */
   int                   initconsssize;      /**< size of initconss array */
   int                   ninitconss;         /**< number of active constraints that must enter the LP */
   int                   ninitconsskept;     /**< number of active constraints that must enter the LP, but were not initialized at
                                              *   their valid node, so that they have to be initialized at every node at which they
                                              *   are active; these constraints come first in the initconss array */
   int                   sepaconsssize;      /**< size of sepaconss array */
   int                   nsepaconss;         /**< number of active constraints that may be separated during LP processing */
   int                   nusefulsepaconss;   /**< number of non-obsolete active constraints that should be separated */
   int                   enfoconsssize;      /**< size of enfoconss array */
   int                   nenfoconss;         /**< number of active constraints that must be enforced during node processing */
   int                   nusefulenfoconss;   /**< number of non-obsolete active constraints that must be enforced */
   int                   checkconsssize;     /**< size of checkconss array */
   int                   ncheckconss;        /**< number of active constraints that must be checked for feasibility */
   int                   nusefulcheckconss;  /**< number of non-obsolete active constraints that must be checked */
   int                   propconsssize;      /**< size of propconss array */
   int                   npropconss;         /**< number of active constraints that may be propagated during node processing */
   int                   nmarkedpropconss;   /**< number of active constraints which are marked to be propagated in the next round */
   int                   nusefulpropconss;   /**< number of non-obsolete active constraints that should be propagated */
   int                   storedpropconsssize;/**< size of array for storing away marked propagation constraints */
   int                   storednmarkedpropconss;/**< number of marked propagation constraints that are stored away */
   int                   updateconsssize;    /**< size of updateconss array */
   int                   nupdateconss;       /**< number of update constraints */
   int                   nenabledconss;      /**< total number of enabled constraints of the handler */
   int                   lastnusefulpropconss;/**< number of already propagated useful constraints on current domains */
   int                   lastnusefulsepaconss;/**< number of already separated useful constraints on current solution */
   int                   lastnusefulenfoconss;/**< number of already enforced useful constraints on current solution */
   int                   lastnfixedvars;     /**< number of variables fixed before the last call to the presolver */
   int                   lastnaggrvars;      /**< number of variables aggregated before the last call to the presolver */
   int                   lastnchgvartypes;   /**< number of variable type changes before the last call to the presolver */
   int                   lastnchgbds;        /**< number of variable bounds tightened before the last call to the presolver */
   int                   lastnaddholes;      /**< number of domain holes added before the last call to the presolver */
   int                   lastndelconss;      /**< number of deleted constraints before the last call to the presolver */
   int                   lastnaddconss;      /**< number of added constraints before the last call to the presolver */
   int                   lastnupgdconss;     /**< number of upgraded constraints before the last call to the presolver */
   int                   lastnchgcoefs;      /**< number of changed coefficients before the last call to the presolver */
   int                   lastnchgsides;      /**< number of changed left or right hand sides before the last call to the presolver */
   int                   nfixedvars;         /**< total number of variables fixed by this presolver */
   int                   naggrvars;          /**< total number of variables aggregated by this presolver */
   int                   nchgvartypes;       /**< total number of variable type changes by this presolver */
   int                   nchgbds;            /**< total number of variable bounds tightened by this presolver */
   int                   naddholes;          /**< total number of domain holes added by this presolver */
   int                   ndelconss;          /**< total number of deleted constraints by this presolver */
   int                   naddconss;          /**< total number of added constraints by this presolver */
   int                   nupgdconss;         /**< total number of upgraded constraints by this presolver */
   int                   nchgcoefs;          /**< total number of changed coefficients by this presolver */
   int                   nchgsides;          /**< total number of changed left or right hand sides by this presolver */
   int                   npresolcalls;       /**< number of times the constraint handler was called in presolving and tried to find reductions */
   int                   delayupdatecount;   /**< must the updates of the constraint arrays be delayed until processUpdates()? */
   SCIP_Bool             delaysepa;          /**< should separation method be delayed, if other separators found cuts? */
   SCIP_Bool             delayprop;          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_Bool             needscons;          /**< should the constraint handler be skipped, if no constraints are available? */
   SCIP_Bool             sepalpwasdelayed;   /**< was the LP separation method delayed at the last call? */
   SCIP_Bool             sepasolwasdelayed;  /**< was the SOL separation method delayed at the last call? */
   SCIP_Bool             propwasdelayed;     /**< was the propagation method delayed at the last call? */
   SCIP_Bool             initialized;        /**< is constraint handler initialized? */
   SCIP_Bool             duringsepa;         /**< is the constraint handler currently performing separation? */
   SCIP_Bool             duringprop;         /**< is the constraint handler currently performing propagation? */
   SCIP_PROPTIMING       proptiming;         /**< positions in the node solving loop where propagation method of constraint handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming;       /**< timing mask of the constraint handler's presolving method */
};

/** linear constraint classification statistics used for MIPLIB */
struct SCIP_LinConsStats
{
   int                   counter[SCIP_NLINCONSTYPES]; /**< count statistics per type of linear constraint */
   int                   sum;                         /**< sum of all counters */
};

#ifdef __cplusplus
}
#endif

#endif
