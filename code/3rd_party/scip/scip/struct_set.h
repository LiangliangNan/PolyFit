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

/**@file   struct_set.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SET_H__
#define __SCIP_STRUCT_SET_H__


#include "scip/def.h"
#include "scip/message.h"
#include "scip/type_bandit.h"
#include "scip/type_set.h"
#include "scip/type_clock.h"
#include "scip/type_paramset.h"
#include "scip/type_event.h"
#include "scip/type_scip.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_cons.h"
#include "scip/type_disp.h"
#include "scip/type_dialog.h"
#include "scip/type_heur.h"
#include "scip/type_compr.h"
#include "scip/type_nodesel.h"
#include "scip/type_presol.h"
#include "scip/type_pricer.h"
#include "scip/type_reader.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_cutsel.h"
#include "scip/type_table.h"
#include "scip/type_prop.h"
#include "scip/type_nlpi.h"
#include "scip/type_concsolver.h"
#include "scip/type_benders.h"
#include "scip/type_expr.h"
#include "scip/debug.h"

#ifdef __cplusplus
extern "C" {
#endif

/** global SCIP settings */
struct SCIP_Set
{
   SCIP_STAGE            stage;              /**< SCIP operation stage */
   SCIP*                 scip;               /**< very ugly: pointer to scip main data structure for callback methods */
   SCIP_PARAMSET*        paramset;           /**< set of parameters */
   BMS_BUFMEM*           buffer;             /**< memory buffers for short living temporary objects */
   BMS_BUFMEM*           cleanbuffer;        /**< memory buffers for short living temporary objects init. to all zero */
   SCIP_READER**         readers;            /**< file readers */
   SCIP_PRICER**         pricers;            /**< variable pricers */
   SCIP_CONSHDLR**       conshdlrs;          /**< constraint handlers (sorted by check priority) */
   SCIP_CONSHDLR**       conshdlrs_sepa;     /**< constraint handlers (sorted by separation priority) */
   SCIP_CONSHDLR**       conshdlrs_enfo;     /**< constraint handlers (sorted by enforcement priority) */
   SCIP_CONSHDLR**       conshdlrs_include;  /**< constraint handlers (sorted by inclusion order) */
   SCIP_CONFLICTHDLR**   conflicthdlrs;      /**< conflict handlers */
   SCIP_PRESOL**         presols;            /**< presolvers */
   SCIP_RELAX**          relaxs;             /**< relaxators */
   SCIP_SEPA**           sepas;              /**< separators */
   SCIP_CUTSEL**         cutsels;            /**< cut selectors */
   SCIP_PROP**           props;              /**< propagators */
   SCIP_PROP**           props_presol;       /**< propagators (sorted by presol priority) */
   SCIP_HEUR**           heurs;              /**< primal heuristics */
   SCIP_COMPR**          comprs;             /**< tree compressions */
   SCIP_EVENTHDLR**      eventhdlrs;         /**< event handlers */
   SCIP_NODESEL**        nodesels;           /**< node selectors */
   SCIP_NODESEL*         nodesel;            /**< currently used node selector, or NULL if invalid */
   SCIP_BRANCHRULE**     branchrules;        /**< branching rules */
   SCIP_DISP**           disps;              /**< display columns */
   SCIP_TABLE**          tables;             /**< statistics tables */
   SCIP_DIALOG**         dialogs;            /**< dialogs */
   SCIP_EXPRHDLR**       exprhdlrs;          /**< expression handlers */
   SCIP_EXPRHDLR*        exprhdlrvar;        /**< expression handler for variables (for quick access) */
   SCIP_EXPRHDLR*        exprhdlrval;        /**< expression handler for constant values (for quick access) */
   SCIP_EXPRHDLR*        exprhdlrsum;        /**< expression handler for sums (for quick access) */
   SCIP_EXPRHDLR*        exprhdlrproduct;    /**< expression handler for products (for quick access) */
   SCIP_EXPRHDLR*        exprhdlrpow;        /**< expression handler for power (for quick access) */
   SCIP_NLPI**           nlpis;              /**< interfaces to NLP solvers */
   SCIP_CONCSOLVERTYPE** concsolvertypes;    /**< concurrent solver types */
   SCIP_CONCSOLVER**     concsolvers;        /**< the concurrent solvers used for solving */
   SCIP_BENDERS**        benders;            /**< the data structures managing the Benders' decomposition algorithm */
   SCIP_DEBUGSOLDATA*    debugsoldata;       /**< data for debug solutions */
   SCIP_BANDITVTABLE**   banditvtables;      /**< virtual function tables for bandit algorithms */
   char**                extcodenames;       /**< names of externals codes */
   char**                extcodedescs;       /**< descriptions of external codes */
   int                   nreaders;           /**< number of file readers */
   int                   readerssize;        /**< size of readers array */
   int                   npricers;           /**< number of variable pricers */
   int                   nactivepricers;     /**< number of variable pricers used in the current problem */
   int                   pricerssize;        /**< size of pricers array */
   int                   nconshdlrs;         /**< number of constraint handlers */
   int                   conshdlrssize;      /**< size of conshdlrs array */
   int                   nconflicthdlrs;     /**< number of conflict handlers */
   int                   conflicthdlrssize;  /**< size of conflicthdlrs array */
   int                   npresols;           /**< number of presolvers */
   int                   presolssize;        /**< size of presols array */
   int                   nrelaxs;            /**< number of relaxators */
   int                   relaxssize;         /**< size of relaxs array */
   int                   nsepas;             /**< number of separators */
   int                   sepassize;          /**< size of sepas array */
   int                   ncutsels;           /**< number of cut selectors */
   int                   cutselssize;        /**< size of cutsels array */
   int                   nprops;             /**< number of propagators */
   int                   propssize;          /**< size of props array */
   int                   nheurs;             /**< number of primal heuristics */
   int                   heurssize;          /**< size of heurs array */
   int                   ncomprs;            /**< number of tree compressions */
   int                   comprssize;         /**< size of comprs array */
   int                   neventhdlrs;        /**< number of event handlers */
   int                   eventhdlrssize;     /**< size of eventhdlrs array */
   int                   nnodesels;          /**< number of node selectors */
   int                   nodeselssize;       /**< size of nodesels array */
   int                   nbranchrules;       /**< number of branching rules */
   int                   branchrulessize;    /**< size of branchrules array */
   int                   ndisps;             /**< number of display columns */
   int                   dispssize;          /**< size of disps array */
   int                   ntables;            /**< number of statistics tables */
   int                   tablessize;         /**< size of tables array */
   int                   ndialogs;           /**< number of dialogs */
   int                   dialogssize;        /**< size of dialogs array */
   int                   nexprhdlrs;         /**< number of expression handlers */
   int                   exprhdlrssize;      /**< size of expression handlers array */
   int                   nnlpis;             /**< number of NLPIs */
   int                   nlpissize;          /**< size of NLPIs array */
   int                   nconcsolvertypes;   /**< number of concurrent solver types */
   int                   concsolvertypessize;/**< size of concurrent solver types array */
   int                   nconcsolvers;       /**< number of concurrent solvers used for solving */
   int                   concsolverssize;    /**< size of concurrent solvers array */
   int                   nbenders;           /**< number of Benders' decomposition algorithms */
   int                   nactivebenders;     /**< number of Benders' decomposition algorithms that are used */
   int                   benderssize;        /**< size of Benders' decomposition algorithms array */
   int                   nextcodes;          /**< number of external codes */
   int                   extcodessize;       /**< size of external code arrays */
   int                   nbanditvtables;     /**< number of bandit algorithm virtual function tables */
   int                   banditvtablessize;  /**< size of banditvtables array */
   SCIP_Bool             pricerssorted;      /**< are the pricers sorted by activity and priority? */
   SCIP_Bool             pricersnamesorted;  /**< are the pricers sorted by name? */
   SCIP_Bool             conflicthdlrssorted;/**< are the conflict handlers sorted by priority? */
   SCIP_Bool             conflicthdlrsnamesorted;/**< are the conflict handlers sorted by name? */
   SCIP_Bool             presolssorted;      /**< are the presolvers sorted by priority? */
   SCIP_Bool             presolsnamesorted;  /**< are the presolvers sorted by name? */
   SCIP_Bool             relaxssorted;       /**< are the relaxators sorted by priority? */
   SCIP_Bool             relaxsnamesorted;   /**< are the relaxators sorted by name? */
   SCIP_Bool             sepassorted;        /**< are the separators sorted by priority? */
   SCIP_Bool             sepasnamesorted;    /**< are the separators sorted by name? */
   SCIP_Bool             cutselssorted;      /**< are the cutsels sorted by priority? */
   SCIP_Bool             propssorted;        /**< are the propagators sorted by priority? */
   SCIP_Bool             propspresolsorted;  /**< are the propagators in prop_presol sorted? */
   SCIP_Bool             propsnamesorted;    /**< are the propagators sorted by name? */
   SCIP_Bool             heurssorted;        /**< are the heuristics sorted by priority? */
   SCIP_Bool             heursnamesorted;    /**< are the heuristics sorted by name? */
   SCIP_Bool             comprssorted;       /**< are the compressions sorted by priority? */
   SCIP_Bool             comprsnamesorted;   /**< are the compressions sorted by name? */
   SCIP_Bool             branchrulessorted;  /**< are the branching rules sorted by priority? */
   SCIP_Bool             branchrulesnamesorted;/**< are the branching rules sorted by name? */
   SCIP_Bool             tablessorted;       /**< are the tables sorted by position? */
   SCIP_Bool             exprhdlrssorted;    /**< are the expression handlers sorted by name? */
   SCIP_Bool             nlpissorted;        /**< are the NLPIs sorted by priority? */
   SCIP_Bool             benderssorted;      /**< are the Benders' algorithms sorted by activity and priority? */
   SCIP_Bool             bendersnamesorted;  /**< are the Benders' algorithms sorted by name? */
   SCIP_Bool             limitchanged;       /**< marks whether any of the limit parameters was changed */
   SCIP_Bool             subscipsoff;        /**< marks whether the sub-SCIPs have been deactivated */

   /* branching settings */
   char                  branch_scorefunc;   /**< branching score function ('s'um, 'p'roduct, 'q'uotient) */
   char                  branch_firstsbchild;/**< child node to be regarded first during strong branching (only with propagation): 'u'p child, 'd'own child, 'h'istory-based, or 'a'utomatic */
   SCIP_Real             branch_scorefac;    /**< branching score factor to weigh downward and upward gain prediction
                                              *   in sum score function */
   SCIP_Bool             branch_preferbinary;/**< should branching on binary variables be preferred? */
   SCIP_Real             branch_clamp;       /**< minimal fractional distance of branching point to a continuous variable' bounds; a value of 0.5 leads to branching always in the middle of a bounded domain */
   SCIP_Real             branch_midpull;     /**< fraction by which to move branching point of a continuous variable towards the middle of the domain; a value of 1.0 leads to branching always in the middle of the domain */
   SCIP_Real             branch_midpullreldomtrig; /**< multiply midpull by relative domain width if the latter is below this value */
   char                  branch_lpgainnorm;  /**< strategy for normalizing LP gain when updating pseudo costs of continuous variables */
   SCIP_Bool             branch_delaypscost; /**< whether to delay pseudo costs updates for continuous variables to after separation */
   SCIP_Bool             branch_divingpscost;/**< should pseudo costs be updated also in diving and probing mode? */
   SCIP_Bool             branch_forceall;    /**< should all strong branching children be regarded even if
                                              *   one is detected to be infeasible? (only with propagation) */
   SCIP_Bool             branch_checksbsol;  /**< should LP solutions during strong branching with propagation be checked for feasibility? */
   SCIP_Bool             branch_roundsbsol;  /**< should LP solutions during strong branching with propagation be rounded? (only when checksbsol=TRUE) */
   SCIP_Bool             branch_sumadjustscore; /**< score adjustment near zero by \b adding epsilon (TRUE) or using maximum (FALSE) */

   /* conflict analysis settings */
   SCIP_Real             conf_maxvarsfac;    /**< maximal fraction of variables involved in a conflict constraint */
   int                   conf_minmaxvars;    /**< minimal absolute maximum of variables involved in a conflict constraint */
   int                   conf_maxlploops;    /**< maximal number of LP resolving loops during conflict analysis
                                              *   (-1: no limit) */
   int                   conf_lpiterations;  /**< maximal number of LP iterations in each LP resolving loop
                                              *   (-1: no limit) */
   int                   conf_fuiplevels;    /**< number of depth levels up to which first UIP's are used in conflict
                                              *   analysis (-1: use All-FirstUIP rule) */
   int                   conf_interconss;    /**< maximal number of intermediate conflict constraints generated in conflict
                                              *   graph (-1: use every intermediate constraint) */
   int                   conf_maxconss;      /**< maximal number of conflict constraints accepted at an infeasible node
                                              *   (-1: use all generated conflict constraints) */
   int                   conf_maxstoresize;  /**< maximal size of conflict store */
   int                   conf_reconvlevels;  /**< number of depth levels up to which UIP reconvergence constraints are
                                              *   generated (-1: generate reconvergence constraints in all depth levels) */
   SCIP_Bool             conf_enable;        /**< should conflict analysis be enabled? */
   SCIP_Bool             conf_cleanbnddepend;/**< should conflicts related to an old cutoff bound be removed? */
   SCIP_Bool             conf_useprop;       /**< should propagation conflict analysis be used? (uses conflict graph only) */
   char                  conf_useinflp;      /**< should infeasible LP conflict analysis be used?
                                              *   ('o'ff, 'c'onflict graph, 'd'ual ray, 'b'oth conflict graph and dual ray)
                                              */
   char                  conf_useboundlp;    /**< should bound exceeding LP conflict analysis be used?
                                              *   ('o'ff, 'c'onflict graph, 'd'ual ray, 'b'oth conflict graph and dual ray)
                                              */
   SCIP_Bool             conf_usesb;         /**< should infeasible/bound exceeding strong branching conflict analysis be
                                              *   used? */
   SCIP_Bool             conf_usepseudo;     /**< should pseudo solution conflict analysis be used? */
   SCIP_Bool             conf_prefinfproof;  /**< prefer infeasibility proof to boundexceeding proof */
   SCIP_Bool             conf_preferbinary;  /**< should binary conflicts be preferred? */
   SCIP_Bool             conf_allowlocal;    /**< should conflict constraints be generated that are only valid locally? */
   SCIP_Bool             conf_settlelocal;   /**< should conflict constraints be attached only to the local subtree where
                                              *   they can be useful? */
   SCIP_Bool             conf_repropagate;   /**< should earlier nodes be repropagated in order to replace branching
                                              *   decisions by deductions? */
   SCIP_Bool             conf_keepreprop;    /**< should constraints be kept for repropagation even if they are too long? */
   SCIP_Bool             conf_separate;      /**< should the conflict constraints be separated? */
   SCIP_Bool             conf_dynamic;       /**< should the conflict constraints be subject to aging? */
   SCIP_Bool             conf_removable;     /**< should the conflict's relaxations be subject to LP aging and cleanup? */
   SCIP_Real             conf_depthscorefac; /**< score factor for depth level in bound relaxation heuristic */
   SCIP_Real             conf_proofscorefac; /**< score factor for contribution to infeasibility proof in bound relaxation heuristic */
   SCIP_Real             conf_uplockscorefac;/**< score factor for number of up locks in bound relaxation heuristic */
   SCIP_Real             conf_downlockscorefac;/**< score factor for number of down locks in bound relaxation heuristic */
   SCIP_Real             conf_scorefac;      /**< factor to decrease importance of variables' earlier conflict scores */
   int                   conf_restartnum;    /**< number of successful conflict analysis calls that trigger a restart
                                              *   (0: disable conflict restarts) */
   SCIP_Real             conf_restartfac;    /**< factor to increase restartnum with after each restart */
   SCIP_Bool             conf_ignorerelaxedbd;/**< should relaxed bounds be ignored? */
   int                   conf_maxvarsdetectimpliedbounds;/**< maximal number of variables to try to detect global bound
                                                          *   implications and shorten the whole conflict set (0:
                                                          *   disabled )
                                                          */
   SCIP_Bool             conf_fullshortenconflict;/**< try to shorten the whole conflict set or terminate early
                                                   *   (depending on the 'maxvarsdetectimpliedbounds' parameter)
                                                   */
   SCIP_Real             conf_conflictweight;/**< the weight the VSIDS score is weight by updating the VSIDS for a
                                              *   variable if it is part of a conflict
                                              */
   SCIP_Real             conf_conflictgraphweight; /**< the weight the VSIDS score is weight by updating the VSIDS for a
                                                    *   variable if it is part of a conflict graph
                                                    */
   SCIP_Real             conf_weightsize;    /**< weight of the size of a conflict used in score calculation */
   SCIP_Real             conf_weightrepropdepth;/**< weight of the prepropagtion depth of a conflict used in score calculation */
   SCIP_Real             conf_weightvaliddepth;/**< weight of the valid depth of a conflict used in score calculation */
   SCIP_Bool             conf_sepaaltproofs;      /**< separate valid inequalities from dualray proofs */
   SCIP_Real             conf_minimprove;    /**< minimal improvement of primal bound to remove conflicts depending on
                                              *   a previous incumbent.
                                              */
   SCIP_Bool             conf_uselocalrows;  /**< use local rows to construct infeasibility proofs */

   /* constraint settings */
   int                   cons_agelimit;      /**< maximum age an unnecessary constraint can reach before it is deleted
                                              *   (0: dynamic, -1: disable aging) */
   int                   cons_obsoleteage;   /**< age of a constraint after which it is marked obsolete
                                              *   (0: dynamic, -1: disable obsoletion) */
   SCIP_Bool             cons_disableenfops; /**< should enforcement of pseudo solution be disabled? */

   /* display settings */
   SCIP_VERBLEVEL        disp_verblevel;     /**< verbosity level of output */
   int                   disp_width;         /**< maximal number of characters in a node information line */
   int                   disp_freq;          /**< frequency for displaying node information lines */
   int                   disp_headerfreq;    /**< frequency for displaying header lines (every n'th node information line) */
   SCIP_Bool             disp_lpinfo;        /**< should the LP solver display status messages? */
   SCIP_Bool             disp_allviols;      /**< display all violations of the best solution after the solving process finished? */
   SCIP_Bool             disp_relevantstats; /**< should the relevant statistics be displayed at the end of solving? */

   /* heuristics settings */
   SCIP_Bool             heur_useuctsubscip; /**< should setting of common subscip parameters include the activation of the UCT node selector? */

   /* history settings */
   SCIP_Bool             history_valuebased; /**< should statistics be collected for variable domain value pairs? */
   SCIP_Bool             history_allowmerge; /**< should variable histories be merged from sub-SCIPs whenever possible? */
   SCIP_Bool             history_allowtransfer; /**< should variable histories be transferred to initialize SCIP copies? */

   /* limit settings */
   SCIP_Real             limit_time;         /**< maximal time in seconds to run */
   SCIP_Real             limit_memory;       /**< maximal memory usage in MB */
   SCIP_Real             limit_gap;          /**< solving stops, if the given gap is reached */
   SCIP_Real             limit_absgap;       /**< solving stops, if the absolute difference between primal and dual bound
                                              *   reaches this value */
   SCIP_Longint          limit_nodes;        /**< maximal number of nodes to process (-1: no limit) */
   SCIP_Longint          limit_totalnodes;   /**< maximal number of total nodes (incl. restarts) to process (-1: no limit) */
   SCIP_Longint          limit_stallnodes;   /**< solving stops, if the given number of nodes was processed since the
                                              *   last improvement of the primal solution value (-1: no limit) */
   int                   limit_solutions;    /**< solving stops, if the given number of solutions were found (-1: no limit) */
   int                   limit_bestsol;      /**< solving stops, if the given number of solution improvements were found
                                              *   (-1: no limit) */
   int                   limit_maxsol;       /**< maximal number of solutions to store in the solution storage */
   int                   limit_maxorigsol;   /**< maximal number of solutions candidates to store in the solution storage of the original problem */
   int                   limit_restarts;     /**< solving stops, if the given number of restarts was triggered (-1: no limit) */
   int                   limit_autorestartnodes;/**< nodes to trigger automatic restart */

   SCIP_Bool             istimelimitfinite;  /**< is the time limit finite */

   /* LP settings */
   int                   lp_solvefreq;       /**< frequency for solving LP at the nodes (-1: never; 0: only root LP) */
   SCIP_Longint          lp_iterlim;         /**< iteration limit for each single LP solve; -1: no limit */
   SCIP_Longint          lp_rootiterlim;     /**< iteration limit for initial root LP solve; -1: no limit */
   int                   lp_solvedepth;      /**< maximal depth for solving LP at the nodes (-1: no depth limit) */
   char                  lp_initalgorithm;   /**< LP algorithm for solving initial LP relaxations ('s'implex, 'b'arrier,
                                              *   barrier with 'c'rossover) */
   char                  lp_resolvealgorithm;/**< LP algorithm for resolving LP relaxations if a starting basis exists
                                              *   ('s'implex, 'b'arrier, barrier with 'c'rossover) */
   char                  lp_pricing;         /**< LP pricing strategy ('a'uto, 'f'ull pricing, 's'teepest edge pricing,
                                              *   'q'uickstart steepest edge pricing, 'd'evex pricing) */
   SCIP_Bool             lp_clearinitialprobinglp;/**< should lp state be cleared at the end of probing mode when LP
                                              *   was initially unsolved, e.g., when called right after presolving? */
   SCIP_Bool             lp_resolverestore;  /**< should the LP be resolved to restore the state at start of diving (if
                                              *   FALSE we buffer the solution values)? */
   SCIP_Bool             lp_freesolvalbuffers; /**< should the buffers for storing LP solution values during diving be
                                              *   freed at end of diving? */
   int                   lp_colagelimit;     /**< maximum age a column can reach before it is deleted from the SCIP_LP
                                              *   (-1: don't delete columns due to aging) */
   int                   lp_rowagelimit;     /**< maximum age a row can reach before it is deleted from the LP 
                                              *   (-1: don't delete rows due to aging) */
   SCIP_Bool             lp_cleanupcols;     /**< should new non-basic columns be removed after LP solving? */
   SCIP_Bool             lp_cleanupcolsroot; /**< should new non-basic columns be removed after root LP solving? */
   SCIP_Bool             lp_cleanuprows;     /**< should new basic rows be removed after LP solving? */
   SCIP_Bool             lp_cleanuprowsroot; /**< should new basic rows be removed after root LP solving? */
   SCIP_Bool             lp_checkstability;  /**< should LP solver's return status be checked for stability? */
   SCIP_Real             lp_conditionlimit;  /**< maximum condition number of LP basis counted as stable (-1.0: no check) */
   SCIP_Real             lp_markowitz;       /**< minimal Markowitz threshold to control sparsity/stability in LU factorization */
   SCIP_Bool             lp_checkprimfeas;   /**< should LP solutions be checked for primal feasibility, resolving LP when numerical troubles occur? */
   SCIP_Bool             lp_checkdualfeas;   /**< should LP solutions be checked for dual feasibility, resolving LP when numerical troubles occur? */
   SCIP_Bool             lp_checkfarkas;     /**< should infeasibility proofs from the LP be checked? */
   int                   lp_fastmip;         /**< which FASTMIP setting of LP solver should be used? 0: off, 1: medium, 2: full */
   int                   lp_scaling;         /**< LP scaling (0: none, 1: normal, 2: aggressive) */
   SCIP_Bool             lp_presolving;      /**< should presolving of LP solver be used? */
   SCIP_Bool             lp_lexdualalgo;     /**< should the lexicographic dual algorithm be used? */
   SCIP_Bool             lp_lexdualrootonly; /**< should the lexicographic dual algorithm be applied only at the root node */
   int                   lp_lexdualmaxrounds;/**< maximum number of rounds in the lexicographic dual algorithm */
   SCIP_Bool             lp_lexdualbasic;    /**< choose fractional basic variables in lexicographic dual algorithm */
   SCIP_Bool             lp_lexdualstalling; /**< turn on the lex dual algorithm only when stalling? */
   int                   lp_disablecutoff;   /**< disable the cutoff bound in the LP solver? (0: enabled, 1: disabled, 2: auto) */
   SCIP_Real             lp_rowrepswitch;    /**< simplex algorithm shall use row representation of the basis
                                              *   if number of rows divided by number of columns exceeds this value */
   int                   lp_threads;         /**< number of threads used for solving the LP (0: automatic) */
   SCIP_Real             lp_resolveiterfac;  /**< factor of average LP iterations that is used as LP iteration limit
                                              *   for LP resolve (-1: unlimited) */
   int                   lp_resolveitermin;  /**< minimum number of iterations that are allowed for LP resolve */
   int                   lp_solutionpolishing;/**< LP solution polishing method (0: disabled, 1: only root, 2: always, 3: auto) */
   int                   lp_refactorinterval;/**< LP refactorization interval (0: automatic) */
   SCIP_Bool             lp_alwaysgetduals;  /**< should the dual solution always be collected for LP solutions. */

   /* NLP settings */
   SCIP_Bool             nlp_disable;        /**< should the NLP be disabled even if a constraint handler enabled it? */
   char*                 nlp_solver;         /**< name of NLP solver to use */

   /* memory settings */
   SCIP_Real             mem_savefac;        /**< fraction of maximal memory usage resulting in switch to memory saving mode */
   SCIP_Real             mem_arraygrowfac;   /**< memory growing factor for dynamically allocated arrays */
   SCIP_Real             mem_treegrowfac;    /**< memory growing factor for tree array */
   SCIP_Real             mem_pathgrowfac;    /**< memory growing factor for path array */
   int                   mem_arraygrowinit;  /**< initial size of dynamically allocated arrays */
   int                   mem_treegrowinit;   /**< initial size of tree array */
   int                   mem_pathgrowinit;   /**< initial size of path array */

   /* miscellaneous settings */
   SCIP_Bool             misc_catchctrlc;    /**< should the CTRL-C interrupt be caught by SCIP? */
   SCIP_Bool             misc_usevartable;   /**< should a hashtable be used to map from variable names to variables? */
   SCIP_Bool             misc_useconstable;  /**< should a hashtable be used to map from constraint names to constraints? */
   SCIP_Bool             misc_usesmalltables;/**< should smaller hashtables be used? yields better performance for small problems with about 100 variables */
   SCIP_Bool             misc_exactsolve;    /**< should the problem be solved exactly (with proven dual bounds)? */
   SCIP_Bool             misc_resetstat;     /**< should the statistics be reset if the transformed problem is freed
                                              *   otherwise the statistics get reset after original problem is freed (in
                                              *   case of bender decomposition this parameter should be set to FALSE and
                                              *   therefore can be used to collect statistics over all runs) */
   SCIP_Bool             misc_improvingsols; /**< should only solutions be checked which improve the primal bound */
   SCIP_Bool             misc_printreason;   /**< should the reason be printed if a given start solution is infeasible? */
   SCIP_Bool             misc_estimexternmem;/**< should the usage of external memory be estimated? */
   SCIP_Bool             misc_avoidmemout;   /**< try to avoid running into memory limit by restricting plugins like heuristics? */
   SCIP_Bool             misc_transorigsols; /**< should SCIP try to transfer original solutions to the transformed space (after presolving)? */
   SCIP_Bool             misc_transsolsorig; /**< should SCIP try to transfer transformed solutions to the original space (after solving)? */
   SCIP_Bool             misc_calcintegral;  /**< should SCIP calculate the primal dual integral value which may require
                                              *   a large number of additional clock calls (and decrease the performance)? */
   SCIP_Bool             misc_finitesolstore;/**< should SCIP try to remove infinite fixings from solutions copied to the solution store? */
   SCIP_Bool             misc_outputorigsol; /**< should the best solution be transformed to the orignal space and be output in command line run? */
   SCIP_Bool             misc_allowstrongdualreds; /**< should strong dual reductions be allowed in propagation and presolving? */
   SCIP_Bool             misc_allowweakdualreds;  /**< should weak dual reductions be allowed in propagation and presolving? */
   SCIP_Real             misc_referencevalue;/**< objective value for reference purposes */
   int                   misc_usesymmetry;   /**< bitset describing used symmetry handling technique (0: off; 1: polyhedral (orbitopes and/or symresacks);
                                              *   2: orbital fixing; 3: orbitopes and orbital fixing; 4: Schreier Sims cuts; 5: Schreier Sims cuts and
                                              *   symresacks) */
   char*                 misc_debugsol;      /**< path to a debug solution */
   SCIP_Bool             misc_scaleobj;      /**< should the objective function be scaled? */
   SCIP_Bool             misc_showdivingstats;/**< should detailed statistics for diving heuristics be shown? */

   /* randomization parameters */
   int                   random_randomseedshift;/**< global shift of all random seeds in the plugins, this will have no impact on the permutation and LP seeds */
   int                   random_permutationseed;/**< seed value for permuting the problem after reading/transformation
                                                 *   (0: no permutation) */
   int                   random_randomseed;     /**< random seed for LP solver, e.g. for perturbations in the simplex (0: LP default) */
   SCIP_Bool             random_permuteconss;   /**< should order of constraints be permuted (depends on permutationseed)? */
   SCIP_Bool             random_permutevars;    /**< should order of variables be permuted (depends on permutationseed)? */

   /* node selection settings */
   char                  nodesel_childsel;   /**< child selection rule ('d'own, 'u'p, 'p'seudo costs, 'i'nference, 'l'p value,
                                              *   'r'oot LP value difference, 'h'brid inference/root LP value difference) */

   /* numerical settings */
   SCIP_Real             num_infinity;       /**< values larger than this are considered infinity */
   SCIP_Real             num_epsilon;        /**< absolute values smaller than this are considered zero */
   SCIP_Real             num_sumepsilon;     /**< absolute values of sums smaller than this are considered zero */
   SCIP_Real             num_feastol;        /**< feasibility tolerance for constraints */
   SCIP_Real             num_checkfeastolfac;/**< factor to change the feasibility tolerance when testing the best
                                              *   solution for feasibility (after solving process) */
   SCIP_Real             num_lpfeastolfactor;/**< factor w.r.t. primal feasibility tolerance that determines default (and maximal) primal feasibility tolerance of LP solver (user parameter, see also num_relaxfeastol) */
   SCIP_Real             num_dualfeastol;    /**< feasibility tolerance for reduced costs */
   SCIP_Real             num_barrierconvtol; /**< convergence tolerance used in barrier algorithm */
   SCIP_Real             num_boundstreps;    /**< minimal improve for strengthening bounds */
   SCIP_Real             num_pseudocosteps;  /**< minimal variable distance value to use for pseudo cost updates */
   SCIP_Real             num_pseudocostdelta;/**< minimal objective distance value to use for pseudo cost updates */
   SCIP_Real             num_recompfac;      /**< minimal decrease factor that causes the recomputation of a value
                                              *   (e.g., pseudo objective) instead of an update */
   SCIP_Real             num_hugeval;        /**< values larger than this are considered huge and should be handled
                                              *   separately (e.g., in activity computation) */
   SCIP_Real             num_relaxfeastol;   /**< primal feasibility tolerance for relaxations (set by core or plugins, not a parameter) */

   /* presolving settings */
   SCIP_Real             presol_abortfac;    /**< abort presolve, if l.t. this frac of the problem was changed in last round */
   int                   presol_maxrounds;   /**< maximal number of presolving rounds (-1: unlimited) */
   int                   presol_maxrestarts; /**< maximal number of restarts (-1: unlimited) */
   SCIP_Real             presol_clqtablefac; /**< limit on number of entries in clique table relative to number of problem nonzeros */
   SCIP_Real             presol_restartfac;  /**< fraction of integer variables that were fixed in the root node
                                              *   triggering a restart with preprocessing after root node evaluation */
   SCIP_Real             presol_immrestartfac;/**< fraction of integer variables that were fixed in the root node triggering an
                                               *   immediate restart with preprocessing */
   SCIP_Real             presol_subrestartfac;/**< fraction of integer variables that were globally fixed during the
                                               *   solving process triggering a restart with preprocessing */
   SCIP_Real             presol_restartminred;/**< minimal fraction of integer variables removed after restart to allow for
                                               *   an additional restart */
   SCIP_Bool             presol_donotmultaggr;/**< should multi-aggregation of variables be forbidden? */
   SCIP_Bool             presol_donotaggr;    /**< should aggregation of variables be forbidden? */

   /* pricing settings */
   SCIP_Real             price_abortfac;     /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
   int                   price_maxvars;      /**< maximal number of variables priced in per pricing round */
   int                   price_maxvarsroot;  /**< maximal number of priced variables at the root node */
   SCIP_Bool             price_delvars;      /**< should variables created at the current node be deleted when the node is solved
                                              *   in case they are not present in the LP anymore? */
   SCIP_Bool             price_delvarsroot;  /**< should variables created at the root node be deleted when the root is solved
                                              *   in case they are not present in the LP anymore? */

   /* Decomposition settings */
   SCIP_Bool             decomp_benderslabels; /**< should the variables be labeled for the application of Benders'
                                                *   decomposition */
   SCIP_Bool             decomp_applybenders;  /**< if a decomposition exists, should Benders' decomposition be applied*/
   int                   decomp_maxgraphedge;  /**< maximum number of edges in block graph computation (-1: no limit, 0: disable block graph computation) */
   SCIP_Bool             decomp_disablemeasures; /**< disable expensive measures */

   /* Benders' decomposition settings */
   SCIP_Real             benders_soltol;     /**< the tolerance for checking optimality in Benders' decomposition */
   SCIP_Bool             benders_cutlpsol;   /**< should cuts be generated from the solution to the LP relaxation? */
   SCIP_Bool             benders_copybenders;/**< should Benders' decomposition be copied for sub-SCIPs? */

   /* propagation settings */
   int                   prop_maxrounds;     /**< maximal number of propagation rounds per node (-1: unlimited) */
   int                   prop_maxroundsroot; /**< maximal number of propagation rounds in the root node (-1: unlimited) */
   SCIP_Bool             prop_abortoncutoff; /**< should propagation be aborted immediately? setting this to FALSE could
                                              *   help conflict analysis to produce more conflict constraints */

   /* reoptimization settings */
   SCIP_Real             reopt_objsimsol;    /**< similarity of two objective functions to reuse stored solutions. */
   SCIP_Real             reopt_objsimrootlp; /**< similarity of two sequential objective function to disable solving the
                                              *   root LP.
                                              */
   SCIP_Real             reopt_objsimdelay;  /**< minimum similarity for using reoptimization of the search tree. */
   char                  reopt_varorderinterdiction; /** use the 'd'efault or a 'r'andom variable order for interdiction
                                                      *  branching when applying the reoptimization
                                                      */
   int                   reopt_forceheurrestart; /**< force a restart if the last n optimal solutions were found by
                                                  *   heuristic reoptsols
                                                  */
   int                   reopt_maxcutage;    /**< maximal age of cuts to use them in reoptimization */
   int                   reopt_maxdiffofnodes;/**< maximal number of bound changes between two stored nodes on one path */
   int                   reopt_maxsavednodes;/**< maximal number of saved nodes */
   int                   reopt_solvelp;      /**< strategy for solving the LP at nodes from reoptimization */
   int                   reopt_solvelpdiff;  /**< maximal number of bound changes at node to skip solving the LP */
   int                   reopt_savesols;     /**< number of best solutions which should be saved for the following runs.
                                              *   (-1: save all)
                                              */
   SCIP_Bool             reopt_commontimelimit;/**< time limit over all reoptimization rounds? */
   SCIP_Bool             reopt_enable;       /**< enable reoptimization */
   SCIP_Bool             reopt_reducetofrontier; /**< delete stored nodes which were not reoptimized */
   SCIP_Bool             reopt_saveconsprop; /**< save constraint propagations */
   SCIP_Bool             reopt_sbinit;       /**< try to fix variables before reoptimizing by probing like strong
                                              *   branching
                                              */
   SCIP_Bool             reopt_shrinkinner;  /**< replace branched inner nodes by their child nodes, if the number of
                                              *   bound changes is not to large
                                              */
   SCIP_Bool             reopt_sepaglbinfsubtrees;/**< save global constraints to separate infeasible subtrees */
   SCIP_Bool             reopt_sepabestsol;  /**< separate only the best solution, i.e., for constrained shortest path */
   SCIP_Bool             reopt_storevarhistory;/**< use variable history of the previous solve if the objective function
                                                *   has changed only slightly
                                                */
   SCIP_Bool             reopt_usepscost;    /**< reuse pseudo costs if the objective function changed only slightly */
   SCIP_Bool             reopt_usecuts;      /**< reoptimize cuts found at the root node */
   SCIP_Bool             reopt_usesplitcons; /**< use constraints to reconstruct the subtree pruned be dual reduction
                                              *   when reactivating the node
                                              */

   /* separation settings */
   SCIP_Real             sepa_maxbounddist;  /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying separation
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_Real             sepa_maxlocalbounddist;/**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying local separation
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_Real             sepa_maxcoefratio;  /**< maximal ratio between coefficients in strongcg, cmir, and flowcover cuts */
   SCIP_Real             sepa_maxcoefratiofacrowprep; /**< maximal ratio between coefficients (as factor of 1/feastol) to ensure in rowprep cleanup */
   SCIP_Real             sepa_minefficacy;   /**< minimal efficacy for a cut to enter the LP */
   SCIP_Real             sepa_minefficacyroot; /**< minimal efficacy for a cut to enter the LP in the root node */
   SCIP_Real             sepa_minortho;      /**< minimal orthogonality for a cut to enter the LP */
   SCIP_Real             sepa_minorthoroot;  /**< minimal orthogonality for a cut to enter the LP in the root node */
   SCIP_Real             sepa_minactivityquot; /**< minimum cut activity quotient to convert cuts into constraints
                                                *   during a restart (0.0: all cuts are converted) */
   char                  sepa_orthofunc;     /**< function used for calc. scalar prod. in orthogonality test ('e'uclidean, 'd'iscrete) */
   char                  sepa_efficacynorm;  /**< row norm to use for efficacy calculation ('e'uclidean, 'm'aximum, 's'um,
                                              *   'd'iscrete) */
   char                  sepa_cutselrestart; /**< cut selection during restart ('a'ge, activity 'q'uotient) */
   char                  sepa_cutselsubscip; /**< cut selection for sub SCIPs  ('a'ge, activity 'q'uotient) */
   SCIP_Bool             sepa_filtercutpoolrel; /**< should cutpool separate only cuts with high relative efficacy? */
   int                   sepa_maxruns;       /**< maximal number of runs for which separation is enabled (-1: unlimited) */
   int                   sepa_maxrounds;     /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   sepa_maxroundsroot; /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   sepa_maxroundsrootsubrun; /**< maximal number of separation rounds in the root node of a subsequent run (-1: unlimited) */
   int                   sepa_maxaddrounds;  /**< maximal additional number of separation rounds in subsequent price-and-cut
                                              *   loops (-1: no additional restriction) */
   int                   sepa_maxstallrounds;/**< maximal number of consecutive separation rounds without objective
                                              *   or integrality improvement (-1: no additional restriction) */
   int                   sepa_maxstallroundsroot;/**< maximal number of consecutive separation rounds without objective
                                              *   or integrality improvement (-1: no additional restriction) */
   int                   sepa_maxcuts;       /**< maximal number of cuts separated per separation round */
   int                   sepa_maxcutsroot;   /**< maximal number of separated cuts at the root node */
   int                   sepa_cutagelimit;   /**< maximum age a cut can reach before it is deleted from the global cut pool */
   int                   sepa_poolfreq;      /**< separation frequency for the global cut pool */

   /* parallel settings */
   int                   parallel_mode;      /**< the mode for the parallel implementation. 0: opportunistic or
                                              *   1: deterministic */
   int                   parallel_minnthreads;/**< the minimum number of threads used for parallel code */
   int                   parallel_maxnthreads;/**< the maximum number of threads used for parallel code */

   /* concurrent solver settings */
   SCIP_Bool             concurrent_changeseeds;    /**< change the seeds in the different solvers? */
   SCIP_Bool             concurrent_changechildsel; /**< change the child selection rule in different solvers? */
   SCIP_Bool             concurrent_commvarbnds;    /**< should the concurrent solvers communicate global variable bound changes? */
   SCIP_Bool             concurrent_presolvebefore; /**< should the problem be presolved before it is copied to the concurrent solvers? */
   int                   concurrent_initseed;       /**< the seed for computing the concurrent solver seeds */
   SCIP_Real             concurrent_freqinit;       /**< initial frequency of synchronization */
   SCIP_Real             concurrent_freqmax;        /**< maximal frequency of synchronization */
   SCIP_Real             concurrent_freqfactor;     /**< factor by which the frequency of synchronization changes */
   SCIP_Real             concurrent_targetprogress; /**< when adapting the synchronization frequency this value is the targeted
                                                     *   relative difference by which the absolute gap decreases per synchronization */
   int                   concurrent_maxnsols;       /**< maximum number of solutions that will get stored in one synchronization */
   int                   concurrent_nbestsols;      /**< number of best solutions that should be considered for synchronization */
   int                   concurrent_maxnsyncdelay;  /**< max number of synchronizations before data is used */
   SCIP_Real             concurrent_minsyncdelay;   /**< min offset before synchronization data is used */
   char*                 concurrent_paramsetprefix; /**< path prefix for parameter setting files of concurrent solver scip-custom */

   /* timing settings */
   SCIP_CLOCKTYPE        time_clocktype;     /**< default clock type to use */
   SCIP_Bool             time_enabled;       /**< is timing enabled? */
   SCIP_Bool             time_reading;       /**< belongs reading time to solving time? */
   SCIP_Bool             time_rareclockcheck;/**< should clock checks of solving time be performed less frequently (might exceed time limit slightly) */
   SCIP_Bool             time_statistictiming;  /**< should timing for statistic output be enabled? */
   SCIP_Bool             time_nlpieval;      /**< should time for evaluation in NLP solves be measured? */

   /* tree compression parameters (for reoptimization) */
   SCIP_Bool             compr_enable;       /**< should automatic tree compression after presolving be enabled? (only for reoptimization) */
   SCIP_Real             compr_time;         /**< maximum time to run tree compression heuristics */

   /* visualization settings */
   char*                 visual_vbcfilename; /**< name of the VBC tool output file, or - if no VBC output should be created */
   char*                 visual_bakfilename; /**< name of the BAK tool output file, or - if no BAK output should be created */
   SCIP_Bool             visual_realtime;    /**< should the real solving time be used instead of time step counter in visualization? */
   SCIP_Bool             visual_dispsols;    /**< should the node where solutions are found be visualized? */
   SCIP_Bool             visual_displb;      /**< should lower bound information be visualized? */
   SCIP_Bool             visual_objextern;   /**< should be output the external value of the objective? */

   /* Reading */
   SCIP_Bool             read_initialconss;  /**< should model constraints be marked as initial? */
   SCIP_Bool             read_dynamicconss;  /**< should model constraints be subject to aging? */
   SCIP_Bool             read_dynamiccols;   /**< should columns be added and removed dynamically to the LP? */
   SCIP_Bool             read_dynamicrows;   /**< should rows be added and removed dynamically to the LP? */

   /* Writing */
   SCIP_Bool             write_allconss;     /**< should all constraints be written (including the redundant constraints)? */
   SCIP_Bool             write_printzeros;   /**< should variables set to zero be printed? */
   int                   write_genoffset;    /**< when writing the problem with generic names, we start with index
                                              *   0; using this parameter we can change the starting index to be
                                              *   different */
};

#ifdef __cplusplus
}
#endif

#endif
