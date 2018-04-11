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

/**@file   struct_stat.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for problem statistics
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Gerald Gamrath
 * @author Marc Pfetsch
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_STAT_H__
#define __SCIP_STRUCT_STAT_H__


#include "scip/def.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "scip/type_visual.h"
#include "scip/type_history.h"
#include "scip/type_var.h"
#include "scip/type_lp.h"
#include "scip/type_heur.h"
#include "scip/type_relax.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** problem and runtime specific statistics */
struct SCIP_Stat
{
   SCIP_REGRESSION*      regressioncandsobjval;/**< linear regression of pairs (nbranchcands, lpobjval) for every node */
   SCIP_Longint          nlpiterations;      /**< total number of LP iterations */
   SCIP_Longint          nrootlpiterations;  /**< total number of LP iterations in root node */
   SCIP_Longint          nrootfirstlpiterations;/**< number of LP iterations for first LP solved at the root node */
   SCIP_Longint          nprimallpiterations;/**< number of iterations in primal simplex */
   SCIP_Longint          nduallpiterations;  /**< number of iterations in dual simplex */
   SCIP_Longint          nlexduallpiterations;/**< number of iterations in lexicographic dual simplex */
   SCIP_Longint          nbarrierlpiterations;/**< number of iterations in barrier algorithm */
   SCIP_Longint          nprimalresolvelpiterations;  /**< number of primal LP iterations with advanced start basis */
   SCIP_Longint          ndualresolvelpiterations;    /**< number of dual LP iterations with advanced start basis */
   SCIP_Longint          nlexdualresolvelpiterations; /**< number of lexicographic dual LP iterations with advanced start basis */
   SCIP_Longint          nnodelpiterations;  /**< number of iterations for totally solving node relaxations */
   SCIP_Longint          ninitlpiterations;  /**< number of iterations for solving nodes' initial relaxations */
   SCIP_Longint          ndivinglpiterations;/**< number of iterations in diving and probing */
   SCIP_Longint          ndivesetlpiterations; /**< total number of LP iterations performed by divesets */
   SCIP_Longint          nsbdivinglpiterations;/**< number of iterations in probing mode for strong branching */
   SCIP_Longint          nsblpiterations;    /**< number of simplex iterations used in strong branching */
   SCIP_Longint          nrootsblpiterations;/**< number of simplex iterations used in strong branching at the root node */
   SCIP_Longint          nconflictlpiterations;/**< number of simplex iterations used in conflict analysis */
   SCIP_Longint          nnodes;             /**< number of nodes processed in current run (including focus node) */
   SCIP_Longint          ninternalnodes;     /**< number of nodes processed in current run where a branching was performed */
   SCIP_Longint          nobjleaves;         /**< number of leaf nodes processed that reached the cutoff bound */
   SCIP_Longint          nfeasleaves;        /**< number of leaf nodes processed with feasible relaxation solution */
   SCIP_Longint          ninfeasleaves;      /**< number of infeasible leaf nodes processed */
   SCIP_Longint          ntotalnodes;        /**< total number of nodes processed in all runs (including focus node) */
   SCIP_Longint          ntotalinternalnodes;/**< total number of nodes processed in all runs where a branching was performed */
   SCIP_Longint          ntotalnodesmerged;  /**< total number of nodes added ot the statistics of the main SCIP so far (see SCIPmergeStatistics) */
   SCIP_Longint          ncreatednodes;      /**< total number of nodes created */
   SCIP_Longint          ncreatednodesrun;   /**< number of nodes created in current run */
   SCIP_Longint          nactivatednodes;    /**< number of times, a node got activated in current run */
   SCIP_Longint          ndeactivatednodes;  /**< number of times, a node got deactivated in current run */
   SCIP_Longint          nearlybacktracks;   /**< counter for early switches (if children dual bound is below reference value) */
   SCIP_Longint          nnodesaboverefbound;/**< counter for the number of focus nodes exceeding the reference bound */
   SCIP_Longint          nbacktracks;        /**< number of times, the new node was chosen from the leaves queue */
   SCIP_Longint          ndelayedcutoffs;    /**< number of times, the selected node was from a cut off subtree */
   SCIP_Longint          nreprops;           /**< number of times, a solved node is repropagated again */
   SCIP_Longint          nrepropboundchgs;   /**< number of bound changes generated in repropagating nodes */
   SCIP_Longint          nrepropcutoffs;     /**< number of times, a repropagated node was cut off */
   SCIP_Longint          nlpsolsfound;       /**< number of CIP-feasible LP solutions found so far */
   SCIP_Longint          nrelaxsolsfound;    /**< number of CIP-feasible relaxation solutions found so far */
   SCIP_Longint          npssolsfound;       /**< number of CIP-feasible pseudo solutions found so far */
   SCIP_Longint          nsbsolsfound;       /**< number of CIP-feasible solutions found during strong branching so far */
   SCIP_Longint          nlpbestsolsfound;   /**< number of new best CIP-feasible LP solutions found so far */
   SCIP_Longint          nrelaxbestsolsfound;/**< number of new best CIP-feasible relaxation solutions found so far */
   SCIP_Longint          npsbestsolsfound;   /**< number of new best CIP-feasible pseudo solutions found so far */
   SCIP_Longint          nsbbestsolsfound;   /**< number of new best CIP-feasible solutions found during strong branching so far */
   SCIP_Longint          nexternalsolsfound; /**< number of externally given CIP-feasible solutions (or new solutions found when transforming old ones) */
   SCIP_Longint          lastdispnode;       /**< last node for which an information line was displayed */
   SCIP_Longint          lastdivenode;       /**< last node where LP diving was applied */
   SCIP_Longint          lastconflictnode;   /**< last node where conflict analysis was applied */
   SCIP_Longint          bestsolnode;        /**< node number where the last incumbent solution was found */
   SCIP_Longint          domchgcount;        /**< internal counter, where all domain changes are counted */
   SCIP_Longint          nboundchgs;         /**< total number of bound changes generated in the tree */
   SCIP_Longint          nholechgs;          /**< total number of hole changes generated in the tree */
   SCIP_Longint          nprobboundchgs;     /**< total number of bound changes generated in the tree during probing */
   SCIP_Longint          nprobholechgs;      /**< total number of hole changes generated in the tree  during probing */
   SCIP_Longint          nsbdowndomchgs;     /**< total number of domain changes generated at down children during strong branching */
   SCIP_Longint          nsbupdomchgs;       /**< total number of domain changes generated at up children during strong branching */
   SCIP_Longint          nsbtimesiterlimhit; /**< total number of times that the strong branching iteration limit was hit */
   SCIP_Longint          nnodesbeforefirst;  /**< number of nodes before first primal solution */
   SCIP_Longint          ninitconssadded;    /**< total number of initial constraints added during the solve */
   SCIP_Longint          nactiveconssadded;  /**< total number of active constraints added */
   SCIP_Longint          externmemestim;     /**< estimation of external memory usage, e.g., by LP solver */
   SCIP_Real             avgnnz;             /**< average number of nonzeros per constraint in presolved problem */
   SCIP_Real             firstlpdualbound;   /**< dual bound of root node computed by first LP solve (without cuts) */
   SCIP_Real             rootlowerbound;     /**< lower bound of root node */
   SCIP_Real             vsidsweight;        /**< current weight to use for updating VSIDS in history */
   SCIP_Real             firstprimalbound;   /**< objective value of first primal solution */
   SCIP_Real             firstprimaltime;    /**< time (in seconds) needed for first primal solution */
   SCIP_Real             firstsolgap;        /**< solution gap when first solution is found */
   SCIP_Real             lastsolgap;         /**< solution gap when last solution is found */
   SCIP_Real             primalzeroittime;   /**< time used in primal simplex calls without iterations */
   SCIP_Real             dualzeroittime;     /**< time used in dual simplex calls without iterations */
   SCIP_Real             barrierzeroittime;  /**< time used in barrier calls without iterations */
   SCIP_Real             maxcopytime;        /**< maxmimal time needed for copying a problem */
   SCIP_Real             mincopytime;        /**< minimal time needed for copying a problem */
   SCIP_Real             firstlptime;        /**< time needed to solve the very first LP in the root node */
   SCIP_Real             lastbranchvalue;    /**< domain value of the last branching */
   SCIP_Real             primaldualintegral; /**< current primal-dual integral value */
   SCIP_Real             previousgap;        /**< primal dual gap preceding the current gap */
   SCIP_Real             previntegralevaltime;/**< last time of primal-dual integral evaluation */
   SCIP_Real             lastprimalbound;    /**< last (non-infinite) primal bound (in transformed space) for integral evaluation */
   SCIP_Real             lastdualbound;      /**< last (non-infinite) dual bound (in transformed space) for integral evaluation */
   SCIP_Real             lastlowerbound;     /**< last lower bound (in transformed space) for integral evaluation */
   SCIP_Real             lastupperbound;     /**< last upper bound (in transformed space) for integral evaluation */
   SCIP_Real             rootlpbestestimate; /**< best-estimate for final root LP solution that changes with every pseudo-cost update */
   SCIP_Real             referencebound;     /**< objective bound for reference purposes */
   SCIP_Real             bestefficacy;       /**< best efficacy of global pool cut seen so far */
   SCIP_Real             minefficacyfac;     /**< factor of best efficacy to use as min efficacy */
   SCIP_Real             detertimecnt;       /**< internal counter for deterministic time */
   SCIP_CLOCK*           solvingtime;        /**< total time used for solving (including presolving) the current problem */
   SCIP_CLOCK*           solvingtimeoverall; /**< total time used for solving (including presolving) during reoptimization */
   SCIP_CLOCK*           presolvingtime;     /**< total time used for presolving the current problem */
   SCIP_CLOCK*           presolvingtimeoverall;/**< total time used for presolving during reoptimization */
   SCIP_CLOCK*           primallptime;       /**< primal LP solution time */
   SCIP_CLOCK*           duallptime;         /**< dual LP solution time */
   SCIP_CLOCK*           lexduallptime;      /**< lexicographic dual LP solution time */
   SCIP_CLOCK*           barrierlptime;      /**< barrier LP solution time */
   SCIP_CLOCK*           divinglptime;       /**< diving and probing LP solution time */
   SCIP_CLOCK*           strongbranchtime;   /**< strong branching time */
   SCIP_CLOCK*           conflictlptime;     /**< conflict analysis LP solution time */
   SCIP_CLOCK*           lpsoltime;          /**< time needed for storing feasible LP solutions */
   SCIP_CLOCK*           relaxsoltime;       /**< time needed for storing feasible relaxation solutions */
   SCIP_CLOCK*           pseudosoltime;      /**< time needed for storing feasible pseudo solutions */
   SCIP_CLOCK*           sbsoltime;          /**< time needed for searching and storing feasible strong branching solutions */
   SCIP_CLOCK*           nodeactivationtime; /**< time needed for path switching and activating nodes */
   SCIP_CLOCK*           nlpsoltime;         /**< time needed for solving NLPs */
   SCIP_CLOCK*           copyclock;          /**< time needed for copying problems */
   SCIP_CLOCK*           strongpropclock;    /**< time needed for propagation during strong branching */
   SCIP_CLOCK*           reoptupdatetime;    /**< time needed for storing and recreating nodes and solutions for reoptimization */
   SCIP_HISTORY*         glbhistory;         /**< global history information over all variables */
   SCIP_HISTORY*         glbhistorycrun;     /**< global history information over all variables for current run */
   SCIP_VAR*             lastbranchvar;      /**< last variable, that was branched on */
   SCIP_VISUAL*          visual;             /**< visualization information */
   SCIP_HEUR*            firstprimalheur;    /**< heuristic which found the first primal solution */
   SCIP_STATUS           status;             /**< SCIP solving status */
   SCIP_BRANCHDIR        lastbranchdir;      /**< direction of the last branching */
   SCIP_LPSOLSTAT        lastsblpsolstats[2];/**< last LP solving statuses for variable strong branching */
   SCIP_Longint          nnz;                /**< number of nonzeros in presolved problem */
   SCIP_Longint          lpcount;            /**< internal counter, where all lp calls are counted; this includes the restored lps after diving and probing */
   SCIP_Longint          relaxcount;         /**< internal counter, where all relax calls are counted */
   SCIP_Longint          nlps;               /**< total number of LPs solved with at least 1 iteration */
   SCIP_Longint          nrootlps;           /**< number of LPs solved at the root node with at least 1 iteration */
   SCIP_Longint          nprimallps;         /**< number of primal LPs solved with at least 1 iteration */
   SCIP_Longint          nprimalzeroitlps;   /**< number of primal LPs with 0 iterations */
   SCIP_Longint          nduallps;           /**< number of dual LPs solved with at least 1 iteration */
   SCIP_Longint          ndualzeroitlps;     /**< number of dual LPs with 0 iterations */
   SCIP_Longint          nlexduallps;        /**< number of lexicographic dual LPs solved */
   SCIP_Longint          nbarrierlps;        /**< number of barrier LPs solved with at least 1 iteration */
   SCIP_Longint          nbarrierzeroitlps;  /**< number of barrier LPs with 1 iteration */
   SCIP_Longint          nprimalresolvelps;  /**< number of primal LPs solved with advanced start basis and at least 1 iteration */
   SCIP_Longint          ndualresolvelps;    /**< number of dual LPs solved with advanced start basis and at least 1 iteration */
   SCIP_Longint          nlexdualresolvelps; /**< number of lexicographic dual LPs solved with advanced start basis and at least 1 iteration */
   SCIP_Longint          nnodelps;           /**< number of LPs solved for node relaxations */
   SCIP_Longint          ninitlps;           /**< number of LPs solved for nodes' initial relaxations */
   SCIP_Longint          ndivinglps;         /**< number of LPs solved during diving and probing */
   SCIP_Longint          ndivesetlps;        /**< total number of diveset LPs */
   SCIP_Longint          nsbdivinglps;       /**< number of LPs solved during strong branching probing mode */
   SCIP_Longint          nnumtroublelpmsgs;  /**< number of messages about numerical trouble in LP on verblevel HIGH or lower */
   SCIP_Longint          nstrongbranchs;     /**< number of strong branching calls */
   SCIP_Longint          nrootstrongbranchs; /**< number of strong branching calls at the root node */
   SCIP_Longint          nconflictlps;       /**< number of LPs solved during conflict analysis */
   SCIP_Longint          nnlps;              /**< number of NLPs solved */
   SCIP_Longint          nisstoppedcalls;    /**< number of calls to SCIPsolveIsStopped() */
   SCIP_Longint          totaldivesetdepth;  /**< the total probing depth over all diveset calls */
   int                   subscipdepth;       /**< depth of current scip instance (increased by each copy call) */
   int                   ndivesetcalls;      /**< total number of diveset diving calls */
   int                   nruns;              /**< number of branch and bound runs on current problem, including current run */
   int                   ncutpoolfails;      /**< number of fails in a cutpool to separate efficacious cuts */
   int                   nconfrestarts;      /**< number of restarts performed due to conflict analysis */
   int                   nrootboundchgs;     /**< total number of bound changes generated in the root node */
   int                   nrootboundchgsrun;  /**< total number of bound changes generated in the root node of current run */
   int                   nrootintfixings;    /**< total number of global fixings of integer variables */
   int                   nrootintfixingsrun; /**< total number of global fixings of integer variables of current run */
   int                   prevrunnvars;       /**< number of variables in the previous run */
   int                   nvaridx;            /**< number of used variable indices */
   int                   ncolidx;            /**< number of used column indices */
   int                   nrowidx;            /**< number of used row indices */
   int                   marked_nvaridx;     /**< number of used variable indices before solving started */
   int                   marked_ncolidx;     /**< number of used column indices before solving started */
   int                   marked_nrowidx;     /**< number of used row indices before solving started */
   int                   npricerounds;       /**< number of pricing rounds performed in current node */
   int                   nseparounds;        /**< number of separation rounds performed in current node */
   int                   nincseparounds;     /**< number of separation rounds performed in current node that increased the maximum number of rows in the LP */
   int                   ndisplines;         /**< number of displayed information lines */
   int                   maxdepth;           /**< maximal depth of all processed nodes in current run */
   int                   maxtotaldepth;      /**< maximal depth of all processed nodes over all runs */
   int                   plungedepth;        /**< current plunging depth (successive times, a child was selected as next node) */
   int                   nactiveconss;       /**< total number of currently active constraints */
   int                   nenabledconss;      /**< total number of currently enabled constraints */
   int                   nimplications;      /**< total number of implications stored in the implication graph */
   int                   npresolrounds;      /**< number of presolving rounds in current run */
   int                   npresolroundsfast;  /**< number of fast presolving rounds in current run */
   int                   npresolroundsmed;   /**< number of medium presolving rounds in current run */
   int                   npresolroundsext;   /**< number of exhaustive presolving rounds in current run */
   int                   npresolfixedvars;   /**< number of presolving fixings in current run */
   int                   npresolaggrvars;    /**< number of presolving aggregations in current run */
   int                   npresolchgvartypes; /**< number of presolving variable type changes in current run */
   int                   npresolchgbds;      /**< number of presolving bound changes in current run */
   int                   npresoladdholes;    /**< number of presolving hole additions in current run */
   int                   npresoldelconss;    /**< number of presolving constraint deletions in current run */
   int                   npresoladdconss;    /**< number of presolving constraint additions in current run */
   int                   npresolupgdconss;   /**< number of presolving constraint upgrades in current run */
   int                   npresolchgcoefs;    /**< number of presolving coefficient changes in current run */
   int                   npresolchgsides;    /**< number of presolving side changes in current run */
   int                   lastnpresolfixedvars;/**< number of presolving fixings before presolving round */
   int                   lastnpresolaggrvars;/**< number of presolving aggregations before presolving round */
   int                   lastnpresolchgvartypes;/**< number of presolving variable type changes before presolving round */
   int                   lastnpresolchgbds;  /**< number of presolving bound changes before presolving round */
   int                   lastnpresoladdholes;/**< number of presolving hole additions before presolving round */
   int                   lastnpresoldelconss;/**< number of presolving constraint deletions before presolving round */
   int                   lastnpresoladdconss;/**< number of presolving constraint additions before presolving round */
   int                   lastnpresolupgdconss;/**< number of presolving constraint upgrades before presolving round */
   int                   lastnpresolchgcoefs;/**< number of presolving coefficient changes before presolving round */
   int                   lastnpresolchgsides;/**< number of presolving side changes before presolving round */
#ifdef SCIP_DISABLED_CODE
   int                   lastnpresolimplications;/**< number of implications before presolving round */
   int                   lastnpresolcliques; /**< number of cliques before presolving round */
#endif
   int                   solindex;           /**< consecutively numbered solution index */
   int                   nrunsbeforefirst;   /**< number of runs until first primal solution */
   int                   firstprimaldepth;   /**< depth in which first primal solution was found */
   int                   ncopies;            /**< counter how often SCIPcopy() was performed */
   int                   nreoptruns;         /**< number of reoptimization runs */
   int                   nclockskipsleft;    /**< how many times the timing should be skipped in SCIPsolveIsStopped() */
   SCIP_Bool             memsavemode;        /**< should algorithms be switched to memory saving mode? */
   SCIP_Bool             userinterrupt;      /**< has the user asked to interrupt the solving process? */
   SCIP_Bool             userrestart;        /**< has the user asked to restart the solving process? */
   SCIP_Bool             inrestart;          /**< are we currently restarting the system? */
   SCIP_Bool             collectvarhistory;  /**< should variable history statistics be collected */
   SCIP_Bool             performpresol;      /**< indicates whether presolving is enabled */
   SCIP_Bool             branchedunbdvar;    /**< indicates whether branching on an unbounded variable has been performed */
   SCIP_Bool             disableenforelaxmsg;/**< was disable enforelax message printed? */
};

#ifdef __cplusplus
}
#endif

#endif
