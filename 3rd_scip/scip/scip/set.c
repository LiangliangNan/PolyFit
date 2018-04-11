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

/**@file   set.c
 * @brief  methods for global SCIP settings
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 * @todo Functions like SCIPsetFeastol() are misleading (it seems that the feasibility tolerance can be set).
 *       Rename all functions starting with SCIPsetXXX, e.g., SCIPsetGetFeastol() and SCIPsetSetFeastol().
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/bandit.h"
#include "scip/branch.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/disp.h"
#include "scip/dialog.h"
#include "scip/heur.h"
#include "scip/concsolver.h"
#include "scip/compr.h"
#include "scip/nodesel.h"
#include "scip/presol.h"
#include "scip/pricer.h"
#include "scip/reader.h"
#include "scip/relax.h"
#include "scip/sepa.h"
#include "scip/table.h"
#include "scip/prop.h"
#include "nlpi/nlpi.h"
#include "scip/struct_scip.h" /* for SCIPsetPrintDebugMessage() */

/*
 * Default settings
 */


/* Branching */

#define SCIP_DEFAULT_BRANCH_SCOREFUNC       'p' /**< branching score function ('s'um, 'p'roduct) */
#define SCIP_DEFAULT_BRANCH_SCOREFAC      0.167 /**< branching score factor to weigh downward and upward gain prediction
                                                 *   in sum score function */
#define SCIP_DEFAULT_BRANCH_PREFERBINARY  FALSE /**< should branching on binary variables be preferred? */
#define SCIP_DEFAULT_BRANCH_CLAMP           0.2 /**< minimal fractional distance of branching point to a continuous variable'
                                                 *   bounds; a value of 0.5 leads to branching always in the middle of a bounded domain */
#define SCIP_DEFAULT_BRANCH_LPGAINNORMALIZE 's' /**< strategy for normalizing LP gain when updating pseudo costs of continuous variables */
#define SCIP_DEFAULT_BRANCH_DELAYPSCOST    TRUE /**< should updating pseudo costs of continuous variables be delayed to after separation */
#define SCIP_DEFAULT_BRANCH_DIVINGPSCOST   TRUE /**< should pseudo costs be updated also in diving and probing mode? */
#define SCIP_DEFAULT_BRANCH_FORCEALL      FALSE /**< should all strong branching children be regarded even if
                                                 *   one is detected to be infeasible? (only with propagation) */
#define SCIP_DEFAULT_BRANCH_FIRSTSBCHILD    'a' /**< child node to be regarded first during strong branching (only with propagation): 'u'p child, 'd'own child, 'h'istory-based, or 'a'utomatic */
#define SCIP_DEFAULT_BRANCH_CHECKSBSOL     TRUE /**< should LP solutions during strong branching with propagation be checked for feasibility? */
#define SCIP_DEFAULT_BRANCH_ROUNDSBSOL     TRUE /**< should LP solutions during strong branching with propagation be rounded? (only when checksbsol=TRUE) */
#define SCIP_DEFAULT_BRANCH_SUMADJUSTSCORE FALSE /**< score adjustment near zero by adding epsilon (TRUE) or using maximum (FALSE) */

/* Tree Compression */

#define SCIP_DEFAULT_COMPR_ENABLE         FALSE /**< should automatic tree compression in reoptimization after presolving be enabled? */


/* Conflict Analysis (general) */

#define SCIP_DEFAULT_CONF_ENABLE           TRUE /**< conflict analysis be enabled? */
#define SCIP_DEFAULT_CONF_MAXVARSFAC       0.15 /**< maximal fraction of variables involved in a conflict constraint */
#define SCIP_DEFAULT_CONF_MINMAXVARS          0 /**< minimal absolute maximum of variables involved in a conflict constraint */
#define SCIP_DEFAULT_CONF_MAXLPLOOPS          2 /**< maximal number of LP resolving loops during conflict analysis
                                                 *   (-1: no limit) */
#define SCIP_DEFAULT_CONF_LPITERATIONS       10 /**< maximal number of LP iterations in each LP resolving loop
                                                 *   (-1: no limit) */
#define SCIP_DEFAULT_CONF_USEPROP          TRUE /**< should propagation conflict analysis be used? */
#define SCIP_DEFAULT_CONF_USEINFLP          'b' /**< should infeasible LP conflict analysis be used?
                                                 *   ('o'ff, 'c'onflict graph, 'd'ual ray, 'b'oth conflict graph and dual ray)
                                                 */
#define SCIP_DEFAULT_CONF_USEBOUNDLP        'b' /**< should bound exceeding LP conflict analysis be used?
                                                 *   ('o'ff, 'c'onflict graph, 'd'ual ray, 'b'oth conflict graph and dual solution)
                                                 */
#define SCIP_DEFAULT_CONF_USESB            TRUE /**< should infeasible/bound exceeding strong branching conflict analysis
                                                 *   be used? */
#define SCIP_DEFAULT_CONF_USEPSEUDO        TRUE /**< should pseudo solution conflict analysis be used? */
#define SCIP_DEFAULT_CONF_PREFINFPROOF     TRUE /**< prefer infeasibility proof to boundexceeding proof */
#define SCIP_DEFAULT_CONF_SEPARATE         TRUE /**< should the conflict constraints be separated? */
#define SCIP_DEFAULT_CONF_DYNAMIC          TRUE /**< should the conflict constraints be subject to aging? */


/* Conflict Analysis (conflict graph) */

#define SCIP_DEFAULT_CONF_MAXSTORESIZE    10000 /**< maximal size of the conflict pool */
#define SCIP_DEFAULT_CONF_RECONVLEVELS       -1 /**< number of depth levels up to which UIP reconvergence constraints are
                                                 *   generated (-1: generate reconvergence constraints in all depth levels) */
#define SCIP_DEFAULT_CONF_CLEANBNDDEPEND   TRUE /**< should conflicts based on an old cutoff bound removed? */
#define SCIP_DEFAULT_CONF_FUIPLEVELS         -1 /**< number of depth levels up to which first UIP's are used in conflict
                                                 *   analysis (-1: use All-FirstUIP rule) */
#define SCIP_DEFAULT_CONF_INTERCONSS         -1 /**< maximal number of intermediate conflict constraints generated in
                                                 *   conflict graph (-1: use every intermediate constraint) */
#define SCIP_DEFAULT_CONF_MAXCONSS           10 /**< maximal number of conflict constraints accepted at an infeasible node
                                                 *   (-1: use all generated conflict constraints) */
#define SCIP_DEFAULT_CONF_PREFERBINARY    FALSE /**< should binary conflicts be preferred? */
#define SCIP_DEFAULT_CONF_ALLOWLOCAL       TRUE /**< should conflict constraints be generated that are only valid locally? */
#define SCIP_DEFAULT_CONF_SETTLELOCAL     FALSE /**< should conflict constraints be attached only to the local subtree
                                                 *   where they can be useful? */
#define SCIP_DEFAULT_CONF_REPROPAGATE      TRUE /**< should earlier nodes be repropagated in order to replace branching
                                                 *   decisions by deductions? */
#define SCIP_DEFAULT_CONF_KEEPREPROP       TRUE /**< should constraints be kept for repropagation even if they are too long? */
#define SCIP_DEFAULT_CONF_REMOVEABLE       TRUE /**< should the conflict's relaxations be subject to LP aging and cleanup? */
#define SCIP_DEFAULT_CONF_DEPTHSCOREFAC     1.0 /**< score factor for depth level in bound relaxation heuristic of LP analysis */
#define SCIP_DEFAULT_CONF_PROOFSCOREFAC     1.0 /**< score factor for impact on acticity in bound relaxation heuristic of LP analysis */
#define SCIP_DEFAULT_CONF_UPLOCKSCOREFAC    0.0 /**< score factor for up locks in bound relaxation heuristic of LP analysis */
#define SCIP_DEFAULT_CONF_DOWNLOCKSCOREFAC  0.0 /**< score factor for down locks in bound relaxation heuristic of LP analysis */
#define SCIP_DEFAULT_CONF_SCOREFAC         0.98 /**< factor to decrease importance of variables' earlier conflict scores */
#define SCIP_DEFAULT_CONF_RESTARTNUM          0 /**< number of successful conflict analysis calls that trigger a restart
                                                 *   (0: disable conflict restarts) */
#define SCIP_DEFAULT_CONF_RESTARTFAC        1.5 /**< factor to increase restartnum with after each restart */
#define SCIP_DEFAULT_CONF_IGNORERELAXEDBD FALSE /**< should relaxed bounds be ignored? */
#define SCIP_DEFAULT_CONF_MAXVARSDETECTIMPLIEDBOUNDS 250 /**< maximal number of variables to try to detect global bound implications and shorten the whole conflict set (0: disabled) */
#define SCIP_DEFAULT_CONF_FULLSHORTENCONFLICT TRUE /**< try to shorten the whole conflict set or terminate early (depending on the 'maxvarsdetectimpliedbounds' parameter) */
#define SCIP_DEFAULT_CONF_CONFLITWEIGHT     0.0 /**< the weight the VSIDS score is weight by updating the VSIDS for a variable if it is part of a conflict */
#define SCIP_DEFAULT_CONF_CONFLITGRAPHWEIGHT 1.0/**< the weight the VSIDS score is weight by updating the VSIDS for a variable if it is part of a conflict graph */
#define SCIP_DEFAULT_CONF_WEIGHTSIZE      0.001 /**< weight of the size of a conflict used in score calculation */
#define SCIP_DEFAULT_CONF_WEIGHTREPROPDEPTH 0.1 /**< weight of the repropagation depth of a conflict used in score calculation */
#define SCIP_DEFAULT_CONF_WEIGHTVALIDDEPTH  1.0 /**< weight of the valid depth of a conflict used in score calculation */
#define SCIP_DEFAULT_CONF_MINIMPROVE       0.05 /**< minimal improvement of primal bound to remove conflicts based on a previous incumbent */

/* Conflict Analysis (dual ray) */

#define SCIP_DEFAULT_CONF_SEPAALTPROOFS   FALSE /**< apply cut generating functions to construct alternative proofs */

/* Constraints */

#define SCIP_DEFAULT_CONS_AGELIMIT            0 /**< maximum age an unnecessary constraint can reach before it is deleted
                                                 *   (0: dynamic adjustment, -1: constraints are never deleted) */
#define SCIP_DEFAULT_CONS_OBSOLETEAGE        -1 /**< age of a constraint after which it is marked obsolete
                                                 *   (0: dynamic adjustment, -1: constraints are never marked obsolete) */
#define SCIP_DEFAULT_CONS_DISABLEENFOPS   FALSE /**< should enforcement of pseudo solution be disabled? */


/* Display */

#define SCIP_DEFAULT_DISP_VERBLEVEL SCIP_VERBLEVEL_HIGH /**< verbosity level of output */
#define SCIP_DEFAULT_DISP_WIDTH             139 /**< maximal number of characters in a node information line */
#define SCIP_DEFAULT_DISP_FREQ              100 /**< frequency for displaying node information lines */
#define SCIP_DEFAULT_DISP_HEADERFREQ         15 /**< frequency for displaying header lines (every n'th node info line) */
#define SCIP_DEFAULT_DISP_LPINFO          FALSE /**< should the LP solver display status messages? */
#define SCIP_DEFAULT_DISP_ALLVIOLS        FALSE /**< display all violations of the best solution after the solving process finished? */

/* History */

#define SCIP_DEFAULT_HISTORY_VALUEBASED   FALSE /**< should statistics be collected for variable domain value pairs */
#define SCIP_DEFAULT_HISTORY_ALLOWMERGE   FALSE /**< should variable histories be merged from sub-SCIPs whenever possible? */
#define SCIP_DEFAULT_HISTORY_ALLOWTRANSFER FALSE /**< should variable histories be transferred to initialize SCIP copies? */

/* Limits */

#define SCIP_DEFAULT_LIMIT_TIME           1e+20 /**< maximal time in seconds to run */
#define SCIP_DEFAULT_LIMIT_MEMORY SCIP_MEM_NOLIMIT/**< maximal memory usage in MB */
#define SCIP_DEFAULT_LIMIT_GAP              0.0 /**< solving stops, if the gap is below the given value */
#define SCIP_DEFAULT_LIMIT_ABSGAP           0.0 /**< solving stops, if the absolute difference between primal and dual
                                                 *   bound reaches this value */
#define SCIP_DEFAULT_LIMIT_NODES           -1LL /**< maximal number of nodes to process (-1: no limit) */
#define SCIP_DEFAULT_LIMIT_STALLNODES      -1LL /**< solving stops, if the given number of nodes was processed since the
                                                 *   last improvement of the primal solution value (-1: no limit) */
#define SCIP_DEFAULT_LIMIT_SOLUTIONS         -1 /**< solving stops, if given number of sols were found (-1: no limit) */
#define SCIP_DEFAULT_LIMIT_BESTSOL           -1 /**< solving stops, if given number of solution improvements were found
                                                 *   (-1: no limit) */
#define SCIP_DEFAULT_LIMIT_MAXSOL           100 /**< maximal number of solutions to store in the solution storage */
#define SCIP_DEFAULT_LIMIT_MAXORIGSOL        10 /**< maximal number of solutions candidates to store in the solution storage of the original problem */
#define SCIP_DEFAULT_LIMIT_RESTARTS          -1 /**< solving stops, if the given number of restarts was triggered (-1: no limit) */
#define SCIP_DEFAULT_LIMIT_AUTORESTARTNODES  -1 /**< if solve exceeds this number of nodes, an automatic restart is triggered (-1: no automatic restart)*/


/* LP */

#define SCIP_DEFAULT_LP_SOLVEFREQ             1 /**< frequency for solving LP at the nodes; -1: never; 0: only root LP */
#define SCIP_DEFAULT_LP_ITERLIM            -1LL /**< iteration limit for each single LP solve; -1: no limit */
#define SCIP_DEFAULT_LP_ROOTITERLIM        -1LL /**< iteration limit for initial root LP solve; -1: no limit */
#define SCIP_DEFAULT_LP_SOLVEDEPTH           -1 /**< maximal depth for solving LPs (-1: no depth limit) */
#define SCIP_DEFAULT_LP_INITALGORITHM       's' /**< LP algorithm for solving initial LP relaxations ('s'implex, 'b'arrier,
                                                 *   barrier with 'c'rossover) */
#define SCIP_DEFAULT_LP_RESOLVEALGORITHM    's' /**< LP algorithm for resolving LP relaxations if a starting basis exists
                                                 *   ('s'implex, 'b'arrier, barrier with 'c'rossover) */
#define SCIP_DEFAULT_LP_PRICING             'l' /**< LP pricing strategy ('l'pi default, 'a'uto, 'f'ull pricing, 'p'artial,
                                                 *   's'teepest edge pricing, 'q'uickstart steepest edge pricing,
                                                 *   'd'evex pricing) */
#define SCIP_DEFAULT_LP_CLEARINITIALPROBINGLP TRUE/**< should lp state be cleared at the end of probing mode when lp
                                                   *   was initially unsolved, e.g., when called right after presolving? */
#define SCIP_DEFAULT_LP_RESOLVERESTORE    FALSE /**< should the LP be resolved to restore the state at start of diving (if FALSE we buffer the solution values)? */
#define SCIP_DEFAULT_LP_FREESOLVALBUFFERS FALSE /**< should the buffers for storing LP solution values during diving be freed at end of diving? */
#define SCIP_DEFAULT_LP_COLAGELIMIT          10 /**< maximum age a dynamic column can reach before it is deleted from SCIP_LP
                                                 *   (-1: don't delete columns due to aging) */
#define SCIP_DEFAULT_LP_ROWAGELIMIT          10 /**< maximum age a dynamic row can reach before it is deleted from SCIP_LP
                                                 *   (-1: don't delete rows due to aging) */
#define SCIP_DEFAULT_LP_CLEANUPCOLS       FALSE /**< should new non-basic columns be removed after LP solving? */
#define SCIP_DEFAULT_LP_CLEANUPCOLSROOT   FALSE /**< should new non-basic columns be removed after root LP solving? */
#define SCIP_DEFAULT_LP_CLEANUPROWS        TRUE /**< should new basic rows be removed after LP solving? */
#define SCIP_DEFAULT_LP_CLEANUPROWSROOT    TRUE /**< should new basic rows be removed after root LP solving? */
#define SCIP_DEFAULT_LP_CHECKSTABILITY     TRUE /**< should LP solver's return status be checked for stability? */
#define SCIP_DEFAULT_LP_CONDITIONLIMIT     -1.0 /**< maximum condition number of LP basis counted as stable (-1.0: no limit) */
#define SCIP_DEFAULT_LP_CHECKPRIMFEAS      TRUE /**< should LP solutions be checked for primal feasibility to resolve LP at numerical troubles? */
#define SCIP_DEFAULT_LP_CHECKDUALFEAS      TRUE /**< should LP solutions be checked for dual feasibility to resolve LP at numerical troubles? */
#define SCIP_DEFAULT_LP_FASTMIP               1 /**< should FASTMIP setting of LP solver be used? */
#define SCIP_DEFAULT_LP_SCALING               1 /**< LP scaling (0: none, 1: normal, 2: aggressive) */
#define SCIP_DEFAULT_LP_PRESOLVING         TRUE /**< should presolving of LP solver be used? */
#define SCIP_DEFAULT_LP_LEXDUALALGO       FALSE /**< should the dual lexicographic algorithm be used? */
#define SCIP_DEFAULT_LP_LEXDUALROOTONLY    TRUE /**< should the lexicographic dual algorithm be applied only at the root node */
#define SCIP_DEFAULT_LP_LEXDUALMAXROUNDS      2 /**< maximum number of rounds in the dual lexicographic algorithm */
#define SCIP_DEFAULT_LP_LEXDUALBASIC      FALSE /**< choose fractional basic variables in lexicographic dual algorithm */
#define SCIP_DEFAULT_LP_LEXDUALSTALLING    TRUE /**< turn on the lex dual algorithm only when stalling? */
#define SCIP_DEFAULT_LP_DISABLECUTOFF         2 /**< disable the cutoff bound in the LP solver? (0: enabled, 1: disabled, 2: auto) */
#define SCIP_DEFAULT_LP_ROWREPSWITCH        1.2 /**< simplex algorithm shall use row representation of the basis
                                                 *   if number of rows divided by number of columns exceeds this value */
#define SCIP_DEFAULT_LP_THREADS               0 /**< number of threads used for solving the LP (0: automatic) */
#define SCIP_DEFAULT_LP_RESOLVEITERFAC     -1.0 /**< factor of average LP iterations that is used as LP iteration limit
                                                 *   for LP resolve (-1.0: unlimited) */
#define SCIP_DEFAULT_LP_RESOLVEITERMIN     1000 /**< minimum number of iterations that are allowed for LP resolve */
#define SCIP_DEFAULT_LP_SOLUTIONPOLISHING     3 /**< LP solution polishing method (0: disabled, 1: only root, 2: always, 3: auto) */
#define SCIP_DEFAULT_LP_REFACTORINTERVAL      0 /**< LP refactorization interval (0: automatic) */

/* NLP */

#define SCIP_DEFAULT_NLP_SOLVER              "" /**< name of NLP solver to use, or "" if solver should be chosen by priority */
#define SCIP_DEFAULT_NLP_DISABLE          FALSE /**< should the NLP be always disabled? */


/* Memory */

#define SCIP_DEFAULT_MEM_SAVEFAC            0.8 /**< fraction of maximal mem usage when switching to memory saving mode */
#define SCIP_DEFAULT_MEM_TREEGROWFAC        2.0 /**< memory growing factor for tree array */
#define SCIP_DEFAULT_MEM_PATHGROWFAC        2.0 /**< memory growing factor for path array */
#define SCIP_DEFAULT_MEM_TREEGROWINIT     65536 /**< initial size of tree array */
#define SCIP_DEFAULT_MEM_PATHGROWINIT       256 /**< initial size of path array */


/* Miscellaneous */

#define SCIP_DEFAULT_MISC_CATCHCTRLC       TRUE /**< should the CTRL-C interrupt be caught by SCIP? */
#define SCIP_DEFAULT_MISC_USEVARTABLE      TRUE /**< should a hashtable be used to map from variable names to variables? */
#define SCIP_DEFAULT_MISC_USECONSTABLE     TRUE /**< should a hashtable be used to map from constraint names to constraints? */
#define SCIP_DEFAULT_MISC_USESMALLTABLES  FALSE /**< should smaller hashtables be used? yields better performance for small problems with about 100 variables */
#define SCIP_DEFAULT_MISC_EXACTSOLVE      FALSE /**< should the problem be solved exactly (with proven dual bounds)? */
#define SCIP_DEFAULT_MISC_RESETSTAT        TRUE /**< should the statistics be reset if the transformed problem is
                                                 *   freed otherwise the statistics get reset after original problem is
                                                 *   freed (in case of Benders decomposition this parameter should be set
                                                 *   to FALSE and therefore can be used to collect statistics over all
                                                 *   runs) */
#define SCIP_DEFAULT_MISC_IMPROVINGSOLS   FALSE /**< should only solutions be checked which improve the primal bound */
#define SCIP_DEFAULT_MISC_PRINTREASON      TRUE /**< should the reason be printed if a given start solution is infeasible? */
#define SCIP_DEFAULT_MISC_ESTIMEXTERNMEM   TRUE /**< should the usage of external memory be estimated? */
#define SCIP_DEFAULT_MISC_TRANSORIGSOLS    TRUE /**< should SCIP try to transfer original solutions to the transformed space (after presolving)? */
#define SCIP_DEFAULT_MISC_TRANSSOLSORIG    TRUE /**< should SCIP try to transfer transformed solutions to the original space (after solving)? */
#define SCIP_DEFAULT_MISC_CALCINTEGRAL     TRUE /**< should SCIP calculate the primal dual integral? */
#define SCIP_DEFAULT_MISC_FINITESOLSTORE  FALSE /**< should SCIP try to remove infinite fixings from solutions copied to the solution store? */
#define SCIP_DEFAULT_MISC_OUTPUTORIGSOL    TRUE /**< should the best solution be transformed to the orignal space and be output in command line run? */
#define SCIP_DEFAULT_MISC_ALLOWDUALREDS    TRUE /**< should dual reductions in propagation methods and presolver be allowed? */
#define SCIP_DEFAULT_MISC_ALLOWOBJPROP     TRUE /**< should propagation to the current objective be allowed in propagation methods? */
#define SCIP_DEFAULT_MISC_REFERENCEVALUE   1e99 /**< objective value for reference purposes */
#define SCIP_DEFAULT_MISC_USESYMMETRY         2 /**< used symmetry handling technique (0: off; 1: polyhedral; 2: orbital fixing) */


#ifdef WITH_DEBUG_SOLUTION
#define SCIP_DEFAULT_MISC_DEBUGSOLUTION     "-" /**< path to a debug solution */
#endif

/* Randomization */
#define SCIP_DEFAULT_RANDOM_RANDSEEDSHIFT     0 /**< global shift of all random seeds in the plugins, this will have no impact on the permutation and LP seeds */
#define SCIP_DEFAULT_RANDOM_PERMUTATIONSEED   0 /**< seed value for permuting the problem after reading/transformation (0: no permutation) */
#define SCIP_DEFAULT_RANDOM_LPSEED            0 /**< random seed for LP solver, e.g. for perturbations in the simplex (0: LP default) */
#define SCIP_DEFAULT_RANDOM_PERMUTECONSS   TRUE /**< should order of constraints be permuted (depends on permutationseed)? */
#define SCIP_DEFAULT_RANDOM_PERMUTEVARS   FALSE /**< should order of variables be permuted (depends on permutationseed)? */


/* Node Selection */

#define SCIP_DEFAULT_NODESEL_CHILDSEL       'h' /**< child selection rule ('d'own, 'u'p, 'p'seudo costs, 'i'nference, 'l'p value,
                                                 *   'r'oot LP value difference, 'h'brid inference/root LP value difference) */


/* Presolving */

#define SCIP_DEFAULT_PRESOL_ABORTFAC      8e-04 /**< abort presolve, if at most this fraction of the problem was changed
                                                 *   in last presolve round */
#define SCIP_DEFAULT_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds (-1: unlimited, 0: off) */
#define SCIP_DEFAULT_PRESOL_MAXRESTARTS      -1 /**< maximal number of restarts (-1: unlimited) */
#define SCIP_DEFAULT_PRESOL_RESTARTFAC    0.025 /**< fraction of integer variables that were fixed in the root node
                                                 *   triggering a restart with preprocessing after root node evaluation */
#define SCIP_DEFAULT_PRESOL_IMMRESTARTFAC  0.10 /**< fraction of integer variables that were fixed in the root node triggering an
                                                 *   immediate restart with preprocessing */
#define SCIP_DEFAULT_PRESOL_SUBRESTARTFAC  1.00 /**< fraction of integer variables that were globally fixed during the
                                                 *   solving process triggering a restart with preprocessing */
#define SCIP_DEFAULT_PRESOL_RESTARTMINRED  0.10 /**< minimal fraction of integer variables removed after restart to allow
                                                 *   for an additional restart */
#define SCIP_DEFAULT_PRESOL_DONOTMULTAGGR FALSE /**< should multi-aggregation of variables be forbidden? */
#define SCIP_DEFAULT_PRESOL_DONOTAGGR     FALSE /**< should aggregation of variables be forbidden? */


/* Pricing */

#define SCIP_DEFAULT_PRICE_ABORTFAC         2.0 /**< pricing is aborted, if fac * price_maxvars pricing candidates were
                                                 *   found */
#define SCIP_DEFAULT_PRICE_MAXVARS          100 /**< maximal number of variables priced in per pricing round */
#define SCIP_DEFAULT_PRICE_MAXVARSROOT     2000 /**< maximal number of priced variables at the root node */
#define SCIP_DEFAULT_PRICE_DELVARS        FALSE /**< should variables created at the current node be deleted when the node is solved
                                                 *   in case they are not present in the LP anymore? */
#define SCIP_DEFAULT_PRICE_DELVARSROOT    FALSE /**< should variables created at the root node be deleted when the root is solved
                                                 *   in case they are not present in the LP anymore? */

/* Reoptimization */

#define SCIP_DEFAULT_REOPT_OBJSIMSOL       -1.0 /**< re-use stored solutions only if the similarity of the new and the old objective
                                                     function is greater or equal than this value */
#define SCIP_DEFAULT_REOPT_OBJSIMROOTLP     0.8 /**< similarity of two sequential objective function to disable solving the root LP. */
#define SCIP_DEFAULT_REOPT_OBJSIMDELAY     -1.0 /**< start reoptimizing the search if the new objective function has similarity of
                                                 *   at least SCIP_DEFAULT_REOPT_DELAY w.r.t. the previous objective function. */
#define SCIP_DEFAULT_REOPT_VARORDERINTERDICTION 'd' /** use 'd'efault, 'r'andom or 'i'nference score for variable
                                                     *  ordering when performing interdiction branching during
                                                     *  reoptimization of nodes
                                                     */
#define SCIP_DEFAULT_REOPT_MAXSAVEDNODES  INT_MAX/**< maximum number of saved nodes */
#define SCIP_DEFAULT_REOPT_MAXDIFFOFNODES INT_MAX/**< maximum number of bound changes of two ancestor nodes
                                                  *  such that the path get not shrunk */
#define SCIP_DEFAULT_REOPT_FORCEHEURRESTART   3 /**< force a restart if the last n optimal solutions are found by
                                                 *   reoptsols heuristic
                                                 */
#define SCIP_DEFAULT_REOPT_SAVESOLS      INT_MAX/**< save n best solutions found so far. */
#define SCIP_DEFAULT_REOPT_SOLVELP            1 /**< strategy for solving the LP at nodes from reoptimization */
#define SCIP_DEFAULT_REOPT_SOLVELPDIFF        1 /**< difference of path length between two ancestor nodes to solve the LP */
#define SCIP_DEFAULT_REOPT_ENABLE         FALSE /**< enable reoptimization */
#define SCIP_DEFAULT_REOPT_SEPAGLBINFSUBTREES TRUE/**< save global constraints to separate infeasible subtrees */
#define SCIP_DEFAULT_REOPT_SEPABESTSOL    FALSE /**< separate the optimal solution, e.g., for solving constraint shortest
                                                 *   path problems
                                                 */
#define SCIP_DEFAULT_REOPT_STOREVARHISTOTY FALSE/**< use the variable history of the previous solve if the objective
                                                 *   function has changed only slightly
                                                 */
#define SCIP_DEFAULT_REOPT_USEPSCOST      FALSE /**< re-use pseudo costs of the objective function changed only slightly */
#define SCIP_DEFAULT_REOPT_COMMONTIMELIMIT FALSE/**< is the given time limit for all reoptimization round? */
#define SCIP_DEFAULT_REOPT_SHRINKINNER     TRUE /**< replace branched transit nodes by their child nodes, if the number
                                                 *   of bound changes is not to large
                                                 */
#define SCIP_DEFAULT_REOPT_STRONGBRANCHINIT TRUE/**< try to fix variables before reoptimizing by probing like strong
                                                 *   branching
                                                 */
#define SCIP_DEFAULT_REOPT_REDUCETOFRONTIER TRUE/**< delete stored nodes which were not reoptimized */
#define SCIP_DEFAULT_REOPT_SAVECONSPROP     FALSE/**< save constraint propagation */
#define SCIP_DEFAULT_REOPT_USESPLITCONS    TRUE /**< use constraints to reconstruct the subtree pruned be dual reduction
                                                 *   when reactivating the node
                                                 */
#define SCIP_DEFAULT_REOPT_USECUTS        FALSE /**< reoptimize cuts found at the root node */
#define SCIP_DEFAULT_REOPT_MAXCUTAGE          0 /**< maximal age of a cut to be use for reoptimization */

/* Propagating */

#define SCIP_DEFAULT_PROP_MAXROUNDS         100 /**< maximal number of propagation rounds per node (-1: unlimited) */
#define SCIP_DEFAULT_PROP_MAXROUNDSROOT    1000 /**< maximal number of propagation rounds in root node (-1: unlimited) */
#define SCIP_DEFAULT_PROP_ABORTONCUTOFF    TRUE /**< should propagation be aborted immediately? setting this to FALSE could
                                                 *   help conflict analysis to produce more conflict constraints */


/* Separation */

#define SCIP_DEFAULT_SEPA_MAXBOUNDDIST      1.0 /**< maximal relative distance from current node's dual bound to primal
                                                 *   bound compared to best node's dual bound for applying separation
                                                 *   (0.0: only on current best node, 1.0: on all nodes) */
#define SCIP_DEFAULT_SEPA_MAXLOCALBOUNDDIST 0.0 /**< maximal relative distance from current node's dual bound to primal
                                                 *   bound compared to best node's dual bound for applying local separation
                                                 *   (0.0: only on current best node, 1.0: on all nodes) */
#define SCIP_DEFAULT_SEPA_MAXCOEFRATIO     1e+4 /**< maximal ratio between coefficients in strongcg, cmir, and flowcover cuts */
#define SCIP_DEFAULT_SEPA_MINEFFICACY      1e-4 /**< minimal efficacy for a cut to enter the LP */
#define SCIP_DEFAULT_SEPA_MINEFFICACYROOT  1e-4 /**< minimal efficacy for a cut to enter the LP in the root node */
#define SCIP_DEFAULT_SEPA_MINORTHO         0.90 /**< minimal orthogonality for a cut to enter the LP */
#define SCIP_DEFAULT_SEPA_MINORTHOROOT     0.90 /**< minimal orthogonality for a cut to enter the LP in the root node */
#define SCIP_DEFAULT_SEPA_OBJPARALFAC       0.1 /**< factor to scale objective parallelism of cut in score calculation */
#define SCIP_DEFAULT_SEPA_INTSUPPORTFAC     0.1 /**< factor to scale integral support of cut in score calculation */
#define SCIP_DEFAULT_SEPA_ORTHOFUNC         'e' /**< function used for calc. scalar prod. in orthogonality test ('e'uclidean, 'd'iscrete) */
#define SCIP_DEFAULT_SEPA_EFFICACYNORM      'e' /**< row norm to use for efficacy calculation ('e'uclidean, 'm'aximum,
                                                 *   's'um, 'd'iscrete) */
#define SCIP_DEFAULT_SEPA_CUTSELRESTART     'a' /**< cut selection during restart ('a'ge, activity 'q'uotient) */
#define SCIP_DEFAULT_SEPA_CUTSELSUBSCIP     'a' /**< cut selection for sub SCIPs  ('a'ge, activity 'q'uotient) */
#define SCIP_DEFAULT_SEPA_MAXRUNS            -1 /**< maximal number of runs for which separation is enabled (-1: unlimited) */
#define SCIP_DEFAULT_SEPA_MAXROUNDS          -1 /**< maximal number of separation rounds per node (-1: unlimited) */
#define SCIP_DEFAULT_SEPA_MAXROUNDSROOT      -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define SCIP_DEFAULT_SEPA_MAXROUNDSROOTSUBRUN -1 /**< maximal number of separation rounds in the root node of a subsequent run (-1: unlimited) */
#define SCIP_DEFAULT_SEPA_MAXADDROUNDS        1 /**< maximal additional number of separation rounds in subsequent
                                                 *   price-and-cut loops (-1: no additional restriction) */
#define SCIP_DEFAULT_SEPA_MAXSTALLROUNDSROOT 10 /**< maximal number of consecutive separation rounds without objective
                                                 *   or integrality improvement in the root node (-1: no additional restriction) */
#define SCIP_DEFAULT_SEPA_MAXSTALLROUNDS      1 /**< maximal number of consecutive separation rounds without objective
                                                 *   or integrality improvement in local nodes (-1: no additional restriction) */
#define SCIP_DEFAULT_SEPA_MAXINCROUNDS       20 /**< maximal number of consecutive separation rounds that increase the size of the LP relaxation per node (-1: unlimited) */
#define SCIP_DEFAULT_SEPA_MAXCUTS           100 /**< maximal number of cuts separated per separation round */
#define SCIP_DEFAULT_SEPA_MAXCUTSROOT      2000 /**< maximal separated cuts at the root node */
#define SCIP_DEFAULT_SEPA_CUTAGELIMIT        80 /**< maximum age a cut can reach before it is deleted from global cut pool
                                                 *   (-1: cuts are never deleted from the global cut pool) */
#define SCIP_DEFAULT_SEPA_POOLFREQ           10 /**< separation frequency for the global cut pool */
#define SCIP_DEFAULT_SEPA_MINACTIVITYQUOT   0.8 /**< minimum cut activity quotient to convert cuts into constraints
                                                 *   during a restart (0.0: all cuts are converted) */

/* Parallel */
#define SCIP_DEFAULT_PARALLEL_MODE               1     /**< the mode for the parallel implementation. Either 0: opportunistic or
                                                        *   1: deterministic */
#define SCIP_DEFAULT_PARALLEL_MINNTHREADS        1     /**< the minimum number of threads used in parallel code */
#define SCIP_DEFAULT_PARALLEL_MAXNTHREADS        8     /**< the maximum number of threads used in parallel code */

/* Concurrent solvers */
#define SCIP_DEFAULT_CONCURRENT_CHANGESEEDS     TRUE /**< should the concurrent solvers use different random seeds? */
#define SCIP_DEFAULT_CONCURRENT_CHANGECHILDSEL  TRUE /**< should the concurrent solvers use different child selection rules? */
#define SCIP_DEFAULT_CONCURRENT_COMMVARBNDS     TRUE /**< should the concurrent solvers communicate variable bounds? */
#define SCIP_DEFAULT_CONCURRENT_PRESOLVEBEFORE  TRUE /**< should the problem be presolved before it is copied to the concurrent solvers? */
#define SCIP_DEFAULT_CONCURRENT_INITSEED     5131912 /**< the seed used to initialize the random seeds for the concurrent solvers */
#define SCIP_DEFAULT_CONCURRENT_FREQINIT        10.0 /**< initial frequency of synchronization with other threads
                                                      *   (fraction of time required for solving the root LP) */
#define SCIP_DEFAULT_CONCURRENT_FREQMAX         10.0 /**< maximal frequency of synchronization with other threads
                                                      *   (fraction of time required for solving the root LP) */
#define SCIP_DEFAULT_CONCURRENT_FREQFACTOR       1.5 /**< factor by which the frequency of synchronization is changed */
#define SCIP_DEFAULT_CONCURRENT_TARGETPROGRESS 0.001 /**< when adapting the synchronization frequency this value is the targeted
                                                       *   relative difference by which the absolute gap decreases per synchronization */
#define SCIP_DEFAULT_CONCURRENT_MAXNSOLS           3 /**< maximum number of solutions that will be shared in a single synchronization */
#define SCIP_DEFAULT_CONCURRENT_MAXNSYNCDELAY      7 /**< maximum number of synchronizations before reading is enforced regardless of delay */
#define SCIP_DEFAULT_CONCURRENT_MINSYNCDELAY    10.0 /**< minimum delay before synchronization data is read */
#define SCIP_DEFAULT_CONCURRENT_NBESTSOLS         10 /**< how many of the N best solutions should be considered for synchronization */
#define SCIP_DEFAULT_CONCURRENT_PARAMSETPREFIX    "" /**< path prefix for parameter setting files of concurrent solvers */


/* Timing */

#define SCIP_DEFAULT_TIME_CLOCKTYPE  SCIP_CLOCKTYPE_CPU  /**< default clock type for timing */
#define SCIP_DEFAULT_TIME_ENABLED          TRUE /**< is timing enabled? */
#define SCIP_DEFAULT_TIME_READING         FALSE /**< belongs reading time to solving time? */
#define SCIP_DEFAULT_TIME_RARECLOCKCHECK  FALSE /**< should clock checks of solving time be performed less frequently (might exceed time limit slightly) */
#define SCIP_DEFAULT_TIME_STATISTICTIMING  TRUE /**< should timing for statistic output be enabled? */


/* visualization output */

#define SCIP_DEFAULT_VISUAL_VBCFILENAME     "-" /**< name of the VBC tool output file, or "-" if no VBC tool output should be created */
#define SCIP_DEFAULT_VISUAL_BAKFILENAME     "-" /**< name of the BAK tool output file, or "-" if no BAK tool output should be created */
#define SCIP_DEFAULT_VISUAL_REALTIME       TRUE /**< should the real solving time be used instead of a time step counter in visualization? */
#define SCIP_DEFAULT_VISUAL_DISPSOLS      FALSE /**< should the node where solutions are found be visualized? */
#define SCIP_DEFAULT_VISUAL_OBJEXTERN      TRUE /**< should be output the external value of the objective? */


/* Reading */

#define SCIP_DEFAULT_READ_INITIALCONSS     TRUE /**< should model constraints be marked as initial? */
#define SCIP_DEFAULT_READ_DYNAMICCONSS     TRUE /**< should model constraints be subject to aging? */
#define SCIP_DEFAULT_READ_DYNAMICCOLS     FALSE /**< should columns be added and removed dynamically to the LP? */
#define SCIP_DEFAULT_READ_DYNAMICROWS     FALSE /**< should rows be added and removed dynamically to the LP? */
#define SCIP_DEFAULT_WRITE_GENNAMES_OFFSET    0 /**< when writing the problem with generic names, we start with index
                                                 *   0; using this parameter we can change the starting index to be
                                                 *   different */


/* Writing */

#define SCIP_DEFAULT_WRITE_ALLCONSS       FALSE /**< should all constraints be written (including the redundant constraints)? */
#define SCIP_DEFAULT_PRINTZEROS           FALSE /**< should variables set to zero be printed? */



/** calculate memory size for dynamically allocated arrays */
static
int calcGrowSize(
   int                   initsize,           /**< initial size of array */
   SCIP_Real             growfac,            /**< growing factor of array */
   int                   num                 /**< minimum number of entries to store */
   )
{
   int size;

   assert(initsize >= 0);
   assert(growfac >= 1.0);
   assert(num >= 0);

   if( growfac == 1.0 )
      size = MAX(initsize, num);
   else
   {
      int oldsize;

      /* calculate the size with this loop, such that the resulting numbers are always the same (-> block memory) */
      initsize = MAX(initsize, 4);
      size = initsize;
      oldsize = size - 1;

      /* second condition checks against overflow */
      while( size < num && size > oldsize )
      {
         oldsize = size;
         size = (int)(growfac * size + initsize);
      }

      /* if an overflow happened, set the correct value */
      if( size <= oldsize )
         size = num;
   }

   assert(size >= initsize);
   assert(size >= num);

   return size;
}


/** information method for a parameter change of feastol */
static
SCIP_DECL_PARAMCHGD(paramChgdFeastol)
{  /*lint --e{715}*/
   SCIP_Real newfeastol;

   newfeastol = SCIPparamGetReal(param);

   /* change the feastol through the SCIP call in order to adjust lpfeastol if necessary */
   SCIP_CALL( SCIPchgFeastol(scip, newfeastol) );

   return SCIP_OKAY;
}

/** information method for a parameter change of lpfeastol */
static
SCIP_DECL_PARAMCHGD(paramChgdLpfeastol)
{  /*lint --e{715}*/
   SCIP_Real newlpfeastol;

   newlpfeastol = SCIPparamGetReal(param);

   /* change the lpfeastol through the SCIP call in order to mark the LP unsolved and control that it does not exceed
    * SCIP's feastol
    */
   SCIP_CALL( SCIPchgLpfeastol(scip, newlpfeastol, FALSE) );

   return SCIP_OKAY;
}

/** information method for a parameter change of dualfeastol */
static
SCIP_DECL_PARAMCHGD(paramChgdDualfeastol)
{  /*lint --e{715}*/
   SCIP_Real newdualfeastol;

   newdualfeastol = SCIPparamGetReal(param);

   /* change the dualfeastol through the SCIP call in order to mark the LP unsolved */
   SCIP_CALL( SCIPchgDualfeastol(scip, newdualfeastol) );

   return SCIP_OKAY;
}

/** information method for a parameter change of barrierconvtol */
static
SCIP_DECL_PARAMCHGD(paramChgdBarrierconvtol)
{  /*lint --e{715}*/
   SCIP_Real newbarrierconvtol;

   newbarrierconvtol = SCIPparamGetReal(param);

   /* change the barrierconvtol through the SCIP call in order to mark the LP unsolved */
   SCIP_CALL( SCIPchgBarrierconvtol(scip, newbarrierconvtol) );

   return SCIP_OKAY;
}

/** information method for a parameter change of infinity value */
static
SCIP_DECL_PARAMCHGD(paramChgInfinity)
{  /*lint --e{715}*/
   SCIP_Real infinity;

   infinity = SCIPparamGetReal(param);

   /* Check that infinity value of LP-solver is at least as large as the one used in SCIP. This is necessary, because we
    * transfer SCIP infinity values to the ones by the LPI, but not the converse. */
   if ( scip->lp != NULL && scip->lp->lpi != NULL && infinity > SCIPlpiInfinity(scip->lp->lpi) )
   {
      SCIPerrorMessage("The infinity value of the LP solver has to be at least as large as the one of SCIP.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   return SCIP_OKAY;
}

/** parameter change information method to autoselect display columns again */
static
SCIP_DECL_PARAMCHGD(SCIPparamChgdDispWidth)
{  /*lint --e{715}*/
   /* automatically select the new active display columns */
   SCIP_CALL( SCIPautoselectDisps(scip) );

   return SCIP_OKAY;
}

/** parameter change information method that some limit was changed */
static
SCIP_DECL_PARAMCHGD(SCIPparamChgdLimit)
{  /*lint --e{715}*/

   SCIPmarkLimitChanged(scip);
   return SCIP_OKAY;
}

/** information method for a parameter change of mem_arraygrowfac */
static
SCIP_DECL_PARAMCHGD(paramChgdArraygrowfac)
{  /*lint --e{715}*/
   SCIP_Real newarraygrowfac;

   newarraygrowfac = SCIPparamGetReal(param);

   /* change arraygrowfac */
   BMSsetBufferMemoryArraygrowfac(SCIPbuffer(scip), newarraygrowfac);
   BMSsetBufferMemoryArraygrowfac(SCIPcleanbuffer(scip), newarraygrowfac);

   return SCIP_OKAY;
}

/** information method for a parameter change of mem_arraygrowinit */
static
SCIP_DECL_PARAMCHGD(paramChgdArraygrowinit)
{  /*lint --e{715}*/
   int newarraygrowinit;

   newarraygrowinit = SCIPparamGetInt(param);

   /* change arraygrowinit */
   BMSsetBufferMemoryArraygrowinit(SCIPbuffer(scip), newarraygrowinit);
   BMSsetBufferMemoryArraygrowinit(SCIPcleanbuffer(scip), newarraygrowinit);

   return SCIP_OKAY;
}

/** information method for a parameter change of reopt_enable */
static
SCIP_DECL_PARAMCHGD(paramChgdEnableReopt)
{  /*lint --e{715}*/

   /* create or deconstruct the reoptimization data structures */

   SCIP_CALL( SCIPenableReoptimization(scip, SCIPparamGetBool(param)) );

   return SCIP_OKAY;
}

/** information method for a parameter change of usesymmetry */
static
SCIP_DECL_PARAMCHGD(paramChgdUsesymmetry)
{  /*lint --e{715}*/

   if ( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_SOLVED )
   {
      if ( SCIPparamGetInt(param) > 0 )
      {
         SCIPerrorMessage("Cannot turn on symmetry handling during (pre)solving.\n");
      }
   }

   return SCIP_OKAY;
}

/** set parameters for reoptimization */
SCIP_RETCODE SCIPsetSetReoptimizationParams(
   SCIP_SET*             set,                /**< SCIP data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert(set != NULL);
   assert(messagehdlr != NULL);

   if( set->reopt_enable )
   {
      /* disable some parts of conflict analysis */
      SCIP_CALL( SCIPsetSetCharParam(set, messagehdlr, "conflict/useboundlp", 'o') );
      SCIP_CALL( SCIPsetSetBoolParam(set, messagehdlr, "conflict/usepseudo", FALSE) );

      /* TODO check wheather multi aggregation can be enabled in reoptimization */
      if( SCIPsetIsParamFixed(set, "presolving/donotmultaggr") )
      {
         SCIP_CALL( SCIPsetChgParamFixed(set, "presolving/donotmultaggr", FALSE) );
      }
      SCIP_CALL( SCIPsetSetBoolParam(set, messagehdlr, "presolving/donotmultaggr", TRUE) );

      if( SCIPsetIsParamFixed(set, "branching/nodereopt/priority") )
      {
         SCIP_CALL( SCIPsetChgParamFixed(set, "branching/nodereopt/priority", FALSE) );
      }
      SCIP_CALL( SCIPsetSetIntParam(set, messagehdlr, "branching/nodereopt/priority", INT_MAX/4) );
   }
   else
   {
      /* disable conflict analysis */
      if( SCIPsetIsParamFixed(set, "conflict/enable") )
      {
         SCIP_CALL( SCIPsetChgParamFixed(set, "conflict/enable", FALSE) );
      }
      SCIP_CALL( SCIPsetResetParam(set, messagehdlr, "conflict/enable") );

      /* TODO check wheather multi aggregation can be enabled in reoptimization */
      if( SCIPsetIsParamFixed(set, "presolving/donotmultaggr") )
      {
         SCIP_CALL( SCIPsetChgParamFixed(set, "presolving/donotmultaggr", FALSE) );
      }
      SCIP_CALL( SCIPsetResetParam(set, messagehdlr, "presolving/donotmultaggr") );

      /* set priority to defeault */
      if( SCIPsetFindBranchrule(set, "nodereopt") != NULL )
      {
         if( SCIPsetIsParamFixed(set, "branching/nodereopt/priority") )
         {
            SCIP_CALL( SCIPsetChgParamFixed(set, "branching/nodereopt/priority", FALSE) );
         }
         SCIP_CALL( SCIPsetResetParam(set, messagehdlr, "branching/nodereopt/priority") );
      }
   }

   return SCIP_OKAY;
}

/** enable or disable all plugin timers depending on the value of the flag \p enabled */
void SCIPsetEnableOrDisablePluginClocks(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_Bool             enabled             /**< should plugin clocks be enabled? */
   )
{
   int i;

   assert(set != NULL);

   /* go through all plugin types and enable or disable their respective clocks */
   for( i = set->nreaders - 1; i >= 0; --i )
      SCIPreaderEnableOrDisableClocks(set->readers[i], enabled);

   for( i = set->npricers - 1; i >= 0; --i )
      SCIPpricerEnableOrDisableClocks(set->pricers[i], enabled);

   for( i = set->nconshdlrs - 1; i >= 0; --i )
      SCIPconshdlrEnableOrDisableClocks(set->conshdlrs[i], enabled);

   for( i = set->nconflicthdlrs - 1; i >= 0; --i )
      SCIPconflicthdlrEnableOrDisableClocks(set->conflicthdlrs[i], enabled);

   for( i = set->npresols - 1; i >= 0; --i )
      SCIPpresolEnableOrDisableClocks(set->presols[i], enabled);

   for( i = set->nrelaxs - 1; i >= 0; --i )
      SCIPrelaxEnableOrDisableClocks(set->relaxs[i], enabled);

   for( i = set->nsepas - 1; i >= 0; --i )
      SCIPsepaEnableOrDisableClocks(set->sepas[i], enabled);

   for( i = set->nprops - 1; i >= 0; --i )
      SCIPpropEnableOrDisableClocks(set->props[i], enabled);

   for( i = set->nheurs - 1; i >= 0; --i )
      SCIPheurEnableOrDisableClocks(set->heurs[i], enabled);

   for( i = set->neventhdlrs - 1; i >= 0; --i )
      SCIPeventhdlrEnableOrDisableClocks(set->eventhdlrs[i], enabled);

   for( i = set->nnodesels - 1; i >= 0; --i )
      SCIPnodeselEnableOrDisableClocks(set->nodesels[i], enabled);

   for( i = set->nbranchrules - 1; i >= 0; --i )
      SCIPbranchruleEnableOrDisableClocks(set->branchrules[i], enabled);
}

/* method to be invoked when the parameter timing/statistictiming is changed */
static
SCIP_DECL_PARAMCHGD(paramChgdStatistictiming)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPenableOrDisableStatisticTiming(scip) );

   return SCIP_OKAY;
}

/** copies plugins from sourcescip to targetscip; in case that a constraint handler which does not need constraints
 *  cannot be copied, valid will return FALSE. All plugins can declare that, if their copy process failed, the
 *  copied SCIP instance might not represent the same problem semantics as the original.
 *  Note that in this case dual reductions might be invalid. */
SCIP_RETCODE SCIPsetCopyPlugins(
   SCIP_SET*             sourceset,          /**< source SCIP_SET data structure */
   SCIP_SET*             targetset,          /**< target SCIP_SET data structure */
   SCIP_Bool             copyreaders,        /**< should the file readers be copied */
   SCIP_Bool             copypricers,        /**< should the variable pricers be copied */
   SCIP_Bool             copyconshdlrs,      /**< should the constraint handlers be copied */
   SCIP_Bool             copyconflicthdlrs,  /**< should the conflict handlers be copied */
   SCIP_Bool             copypresolvers,     /**< should the presolvers be copied */
   SCIP_Bool             copyrelaxators,     /**< should the relaxators be copied */
   SCIP_Bool             copyseparators,     /**< should the separators be copied */
   SCIP_Bool             copypropagators,    /**< should the propagators be copied */
   SCIP_Bool             copyheuristics,     /**< should the heuristics be copied */
   SCIP_Bool             copyeventhdlrs,     /**< should the event handlers be copied */
   SCIP_Bool             copynodeselectors,  /**< should the node selectors be copied */
   SCIP_Bool             copybranchrules,    /**< should the branchrules be copied */
   SCIP_Bool             copydisplays,       /**< should the display columns be copied */
   SCIP_Bool             copydialogs,        /**< should the dialogs be copied */
   SCIP_Bool             copytables,         /**< should the statistics tables be copied */
   SCIP_Bool             copynlpis,          /**< should the NLP interfaces be copied */
   SCIP_Bool*            allvalid            /**< pointer to store whether all plugins were validly copied */
   )
{
   int p;
   SCIP_Bool valid;

   assert(sourceset != NULL);
   assert(targetset != NULL);
   assert(sourceset != targetset);
   assert(allvalid != NULL);

   *allvalid = TRUE;

   /* copy all reader plugins */
   if( copyreaders && sourceset->readers != NULL )
   {
      for( p = sourceset->nreaders - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPreaderCopyInclude(sourceset->readers[p], targetset) );
      }
   }

   /* copy all variable pricer plugins */
   if( copypricers && sourceset->pricers != NULL )
   {
      for( p = sourceset->npricers - 1; p >= 0; --p )
      {
         valid = FALSE;
         SCIP_CALL( SCIPpricerCopyInclude(sourceset->pricers[p], targetset, &valid) );
         *allvalid = *allvalid && valid;
         if( SCIPpricerIsActive(sourceset->pricers[p]) )
         {
            SCIP_CALL( SCIPactivatePricer(targetset->scip, targetset->pricers[p]) );
         }
      }
   }

   /* copy all constraint handler plugins */
   if( copyconshdlrs && sourceset->conshdlrs_include != NULL )
   {
      /* copy them in order they were added to the sourcescip
       *
       * @note we only have to set the valid pointer to FALSE in case that a constraint handler, which does not need
       *       constraints, does not copy; in the case that a constraint handler does not copy and it needs constraint
       *       we will detect later that the problem is not valid if a constraint of that type exits
       */
      for( p = 0; p < sourceset->nconshdlrs; ++p )
      {
         if( SCIPconshdlrIsClonable(sourceset->conshdlrs_include[p]) )
         {
            valid = FALSE;
            SCIP_CALL( SCIPconshdlrCopyInclude(sourceset->conshdlrs_include[p], targetset, &valid) );
            *allvalid = *allvalid && valid;
            SCIPsetDebugMsg(sourceset, "Copying conshdlr <%s> was%s valid.\n", SCIPconshdlrGetName(sourceset->conshdlrs_include[p]), valid ? "" : " not");
         }
         else if( !SCIPconshdlrNeedsCons(sourceset->conshdlrs_include[p]) )
         {
            SCIPsetDebugMsg(sourceset, "Copying Conshdlr <%s> without constraints not valid.\n", SCIPconshdlrGetName(sourceset->conshdlrs_include[p]));
            *allvalid = FALSE;
         }
      }
   }

   /* copy all conflict handler plugins */
   if( copyconflicthdlrs && sourceset->conflicthdlrs != NULL )
   {
      for( p = sourceset->nconflicthdlrs - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPconflicthdlrCopyInclude(sourceset->conflicthdlrs[p], targetset) );
      }
   }

   /* copy all presolver plugins */
   if( copypresolvers && sourceset->presols != NULL )
   {
      for( p = sourceset->npresols - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPpresolCopyInclude(sourceset->presols[p], targetset) );
      }
   }


   /* copy all relaxator plugins */
   if( copyrelaxators && sourceset->relaxs != NULL )
   {
      for( p = sourceset->nrelaxs - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPrelaxCopyInclude(sourceset->relaxs[p], targetset) );
      }
   }


   /* copy all separator plugins */
   if( copyseparators && sourceset->sepas != NULL )
   {
      for( p = sourceset->nsepas - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPsepaCopyInclude(sourceset->sepas[p], targetset) );
      }
   }

   /* copy all propagators plugins */
   if( copypropagators && sourceset->props != NULL )
   {
      for( p = sourceset->nprops - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPpropCopyInclude(sourceset->props[p], targetset) );
      }
   }

   /* copy all primal heuristics plugins */
   if( copyheuristics && sourceset->heurs != NULL )
   {
      for( p = sourceset->nheurs - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPheurCopyInclude(sourceset->heurs[p], targetset) );
      }
   }

   /* copy all event handler plugins */
   if( copyeventhdlrs && sourceset->eventhdlrs != NULL )
   {
      for( p = sourceset->neventhdlrs - 1; p >= 0; --p )
      {
         /* @todo: the copying process of event handlers is currently not checked for consistency */
         SCIP_CALL( SCIPeventhdlrCopyInclude(sourceset->eventhdlrs[p], targetset) );
      }
   }


   /* copy all node selector plugins */
   if( copynodeselectors && sourceset->nodesels != NULL )
   {
      for( p = sourceset->nnodesels - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPnodeselCopyInclude(sourceset->nodesels[p], targetset) );
      }
   }

   /* copy all branchrule plugins */
   if( copybranchrules && sourceset->branchrules != NULL )
   {
      for( p = sourceset->nbranchrules - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPbranchruleCopyInclude(sourceset->branchrules[p], targetset) );
      }
   }


   /* copy all display plugins */
   if( copydisplays && sourceset->disps != NULL )
   {
      for( p = sourceset->ndisps - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPdispCopyInclude(sourceset->disps[p], targetset) );
      }
   }


   /* copy all dialog plugins */
   if( copydialogs && sourceset->dialogs != NULL )
   {
      for( p = sourceset->ndialogs - 1; p >= 0; --p )
      {
         /* @todo: the copying process of dialog handlers is currently not checked for consistency */
         SCIP_CALL( SCIPdialogCopyInclude(sourceset->dialogs[p], targetset) );
      }
   }

   /* copy all table plugins */
   if( copytables && sourceset->tables != NULL )
   {
      for( p = sourceset->ntables - 1; p >= 0; --p )
      {
         SCIP_CALL( SCIPtableCopyInclude(sourceset->tables[p], targetset) );
      }
   }

   /* copy all NLP interfaces */
   if( copynlpis && sourceset->nlpis != NULL )
   {
      for( p = sourceset->nnlpis - 1; p >= 0; --p )
      {
         SCIP_NLPI* nlpicopy;

         SCIP_CALL( SCIPnlpiCopy(SCIPblkmem(targetset->scip), sourceset->nlpis[p], &nlpicopy) );
         SCIP_CALL( SCIPincludeNlpi(targetset->scip, nlpicopy) );
      }
   }

   return SCIP_OKAY;
}

/** copies parameters from sourcescip to targetscip */
SCIP_RETCODE SCIPsetCopyParams(
   SCIP_SET*             sourceset,          /**< source SCIP_SET data structure */
   SCIP_SET*             targetset,          /**< target SCIP_SET data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler of target SCIP */
   )
{
   assert(sourceset != NULL);
   assert(targetset != NULL);
   assert(sourceset != targetset);
   assert(targetset->scip != NULL);

   SCIP_CALL( SCIPparamsetCopyParams(sourceset->paramset, targetset->paramset, targetset, messagehdlr) );

   return SCIP_OKAY;
}

/** creates global SCIP settings */
SCIP_RETCODE SCIPsetCreate(
   SCIP_SET**            set,                /**< pointer to SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(set != NULL);
   assert(scip != NULL);

   SCIP_ALLOC( BMSallocMemory(set) );

   (*set)->stage = SCIP_STAGE_INIT;
   (*set)->scip = scip;
   (*set)->buffer = SCIPbuffer(scip);
   (*set)->cleanbuffer = SCIPcleanbuffer(scip);

   SCIP_CALL( SCIPparamsetCreate(&(*set)->paramset, blkmem) );

   (*set)->readers = NULL;
   (*set)->nreaders = 0;
   (*set)->readerssize = 0;
   (*set)->pricers = NULL;
   (*set)->npricers = 0;
   (*set)->nactivepricers = 0;
   (*set)->pricerssize = 0;
   (*set)->pricerssorted = FALSE;
   (*set)->pricersnamesorted = FALSE;
   (*set)->conshdlrs = NULL;
   (*set)->conshdlrs_sepa = NULL;
   (*set)->conshdlrs_enfo = NULL;
   (*set)->conshdlrs_include = NULL;
   (*set)->nconshdlrs = 0;
   (*set)->conshdlrssize = 0;
   (*set)->conflicthdlrs = NULL;
   (*set)->nconflicthdlrs = 0;
   (*set)->conflicthdlrssize = 0;
   (*set)->conflicthdlrssorted = FALSE;
   (*set)->conflicthdlrsnamesorted = FALSE;

   (*set)->debugsoldata = NULL;
   SCIP_CALL( SCIPdebugSolDataCreate(&(*set)->debugsoldata) ); /*lint !e506 !e774*/

   (*set)->presols = NULL;
   (*set)->npresols = 0;
   (*set)->presolssize = 0;
   (*set)->presolssorted = FALSE;
   (*set)->presolsnamesorted = FALSE;
   (*set)->relaxs = NULL;
   (*set)->nrelaxs = 0;
   (*set)->relaxssize = 0;
   (*set)->relaxssorted = FALSE;
   (*set)->relaxsnamesorted = FALSE;
   (*set)->sepas = NULL;
   (*set)->nsepas = 0;
   (*set)->sepassize = 0;
   (*set)->sepassorted = FALSE;
   (*set)->sepasnamesorted = FALSE;
   (*set)->props = NULL;
   (*set)->props_presol = NULL;
   (*set)->nprops = 0;
   (*set)->propssize = 0;
   (*set)->propssorted = FALSE;
   (*set)->propspresolsorted = FALSE;
   (*set)->propsnamesorted = FALSE;
   (*set)->concsolvertypes = NULL;
   (*set)->nconcsolvertypes = 0;
   (*set)->concsolvertypessize = 0;
   (*set)->concsolvers = NULL;
   (*set)->nconcsolvers = 0;
   (*set)->concsolverssize = 0;
   (*set)->concurrent_paramsetprefix = NULL;
   (*set)->heurs = NULL;
   (*set)->nheurs = 0;
   (*set)->heurssize = 0;
   (*set)->heurssorted = FALSE;
   (*set)->heursnamesorted = FALSE;
   (*set)->comprs = NULL;
   (*set)->ncomprs = 0;
   (*set)->comprssize = 0;
   (*set)->comprssorted = FALSE;
   (*set)->comprsnamesorted = FALSE;
   (*set)->eventhdlrs = NULL;
   (*set)->neventhdlrs = 0;
   (*set)->eventhdlrssize = 0;
   (*set)->nodesels = NULL;
   (*set)->nnodesels = 0;
   (*set)->nodeselssize = 0;
   (*set)->nodesel = NULL;
   (*set)->branchrules = NULL;
   (*set)->nbranchrules = 0;
   (*set)->branchrulessize = 0;
   (*set)->branchrulessorted = FALSE;
   (*set)->branchrulesnamesorted = FALSE;
   (*set)->banditvtables = NULL;
   (*set)->banditvtablessize = 0;
   (*set)->nbanditvtables = 0;
   (*set)->disps = NULL;
   (*set)->ndisps = 0;
   (*set)->dispssize = 0;
   (*set)->tables = NULL;
   (*set)->ntables = 0;
   (*set)->tablessize = 0;
   (*set)->tablessorted = FALSE;
   (*set)->dialogs = NULL;
   (*set)->ndialogs = 0;
   (*set)->dialogssize = 0;
   (*set)->nlpis = NULL;
   (*set)->nnlpis = 0;
   (*set)->nlpissize = 0;
   (*set)->nlpissorted = FALSE;
   (*set)->limitchanged = FALSE;
   (*set)->extcodenames = NULL;
   (*set)->extcodedescs = NULL;
   (*set)->nextcodes = 0;
   (*set)->extcodessize = 0;
   (*set)->visual_vbcfilename = NULL;
   (*set)->visual_bakfilename = NULL;
   (*set)->nlp_solver = NULL;
   (*set)->nlp_disable = FALSE;
   (*set)->num_relaxfeastol = SCIP_INVALID;
   (*set)->misc_debugsol = NULL;

   /* the default time limit is infinite */
   (*set)->istimelimitfinite = FALSE;

   /* branching parameters */
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "branching/scorefunc",
         "branching score function ('s'um, 'p'roduct, 'q'uotient)",
         &(*set)->branch_scorefunc, TRUE, SCIP_DEFAULT_BRANCH_SCOREFUNC, "spq",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "branching/scorefac",
         "branching score factor to weigh downward and upward gain prediction in sum score function",
         &(*set)->branch_scorefac, TRUE, SCIP_DEFAULT_BRANCH_SCOREFAC, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "branching/preferbinary",
         "should branching on binary variables be preferred?",
         &(*set)->branch_preferbinary, FALSE, SCIP_DEFAULT_BRANCH_PREFERBINARY,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "branching/clamp",
         "minimal relative distance of branching point to bounds when branching on a continuous variable",
         &(*set)->branch_clamp, FALSE, SCIP_DEFAULT_BRANCH_CLAMP, 0.0, 0.5,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "branching/lpgainnormalize",
         "strategy for normalization of LP gain when updating pseudocosts of continuous variables (divide by movement of 'l'p value, reduction in 'd'omain width, or reduction in domain width of 's'ibling)",
         &(*set)->branch_lpgainnorm, FALSE, SCIP_DEFAULT_BRANCH_LPGAINNORMALIZE, "dls",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "branching/delaypscostupdate",
         "should updating pseudo costs for continuous variables be delayed to the time after separation?",
         &(*set)->branch_delaypscost, FALSE, SCIP_DEFAULT_BRANCH_DELAYPSCOST,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "branching/divingpscost",
         "should pseudo costs be updated also in diving and probing mode?",
         &(*set)->branch_divingpscost, FALSE, SCIP_DEFAULT_BRANCH_DIVINGPSCOST,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "branching/forceallchildren",
         "should all strong branching children be regarded even if one is detected to be infeasible? (only with propagation)",
         &(*set)->branch_forceall, TRUE, SCIP_DEFAULT_BRANCH_FORCEALL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "branching/firstsbchild",
         "child node to be regarded first during strong branching (only with propagation): 'u'p child, 'd'own child, 'h'istory-based, or 'a'utomatic",
         &(*set)->branch_firstsbchild, TRUE, SCIP_DEFAULT_BRANCH_FIRSTSBCHILD, "aduh",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "branching/checksol",
         "should LP solutions during strong branching with propagation be checked for feasibility?",
         &(*set)->branch_checksbsol, TRUE, SCIP_DEFAULT_BRANCH_CHECKSBSOL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "branching/roundsbsol",
         "should LP solutions during strong branching with propagation be rounded? (only when checksbsol=TRUE)",
         &(*set)->branch_roundsbsol, TRUE, SCIP_DEFAULT_BRANCH_ROUNDSBSOL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "branching/sumadjustscore",
         "score adjustment near zero by adding epsilon (TRUE) or using maximum (FALSE)",
         &(*set)->branch_sumadjustscore, TRUE, SCIP_DEFAULT_BRANCH_SUMADJUSTSCORE,
         NULL, NULL) );

   /* tree compression parameters */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "compression/enable",
         "should automatic tree compression after the presolving be enabled?",
         &(*set)->compr_enable, TRUE, SCIP_DEFAULT_COMPR_ENABLE,
         NULL, NULL) );

   /* conflict analysis parameters */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/enable",
         "should conflict analysis be enabled?",
         &(*set)->conf_enable, FALSE, SCIP_DEFAULT_CONF_ENABLE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/cleanboundexceedings",
         "should conflicts based on an old cutoff bound be removed from the conflict pool after improving the primal bound?",
         &(*set)->conf_cleanbnddepend, TRUE, SCIP_DEFAULT_CONF_CLEANBNDDEPEND,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/useprop",
         "should propagation conflict analysis be used?",
         &(*set)->conf_useprop, FALSE, SCIP_DEFAULT_CONF_USEPROP,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "conflict/useinflp",
         "should infeasible LP conflict analysis be used? ('o'ff, 'c'onflict graph, 'd'ual ray, 'b'oth conflict graph and dual ray)",
         &(*set)->conf_useinflp, FALSE, SCIP_DEFAULT_CONF_USEINFLP, "ocdb",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "conflict/useboundlp",
         "should bound exceeding LP conflict analysis be used? ('o'ff, 'c'onflict graph, 'd'ual ray, 'b'oth conflict graph and dual ray)",
         &(*set)->conf_useboundlp, FALSE, SCIP_DEFAULT_CONF_USEBOUNDLP, "ocdb",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/usesb",
         "should infeasible/bound exceeding strong branching conflict analysis be used?",
         &(*set)->conf_usesb, FALSE, SCIP_DEFAULT_CONF_USESB,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/usepseudo",
         "should pseudo solution conflict analysis be used?",
         &(*set)->conf_usepseudo, FALSE, SCIP_DEFAULT_CONF_USEPSEUDO,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/maxvarsfac",
         "maximal fraction of variables involved in a conflict constraint",
         &(*set)->conf_maxvarsfac, TRUE, SCIP_DEFAULT_CONF_MAXVARSFAC, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/minmaxvars",
         "minimal absolute maximum of variables involved in a conflict constraint",
         &(*set)->conf_minmaxvars, TRUE, SCIP_DEFAULT_CONF_MINMAXVARS, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/maxlploops",
         "maximal number of LP resolving loops during conflict analysis (-1: no limit)",
         &(*set)->conf_maxlploops, TRUE, SCIP_DEFAULT_CONF_MAXLPLOOPS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/lpiterations",
         "maximal number of LP iterations in each LP resolving loop (-1: no limit)",
         &(*set)->conf_lpiterations, TRUE, SCIP_DEFAULT_CONF_LPITERATIONS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/fuiplevels",
         "number of depth levels up to which first UIP's are used in conflict analysis (-1: use All-FirstUIP rule)",
         &(*set)->conf_fuiplevels, TRUE, SCIP_DEFAULT_CONF_FUIPLEVELS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/interconss",
         "maximal number of intermediate conflict constraints generated in conflict graph (-1: use every intermediate constraint)",
         &(*set)->conf_interconss, TRUE, SCIP_DEFAULT_CONF_INTERCONSS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/reconvlevels",
         "number of depth levels up to which UIP reconvergence constraints are generated (-1: generate reconvergence constraints in all depth levels)",
         &(*set)->conf_reconvlevels, TRUE, SCIP_DEFAULT_CONF_RECONVLEVELS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/maxconss",
         "maximal number of conflict constraints accepted at an infeasible node (-1: use all generated conflict constraints)",
         &(*set)->conf_maxconss, TRUE, SCIP_DEFAULT_CONF_MAXCONSS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/maxstoresize",
         "maximal size of conflict store (-1: auto, 0: disable storage)",
         &(*set)->conf_maxstoresize, TRUE, SCIP_DEFAULT_CONF_MAXSTORESIZE, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/preferbinary",
         "should binary conflicts be preferred?",
         &(*set)->conf_preferbinary, FALSE, SCIP_DEFAULT_CONF_PREFERBINARY,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/prefinfproof",
         "prefer infeasibility proof to boundexceeding proof",
         &(*set)->conf_prefinfproof, TRUE, SCIP_DEFAULT_CONF_PREFINFPROOF,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/allowlocal",
         "should conflict constraints be generated that are only valid locally?",
         &(*set)->conf_allowlocal, TRUE, SCIP_DEFAULT_CONF_ALLOWLOCAL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/settlelocal",
         "should conflict constraints be attached only to the local subtree where they can be useful?",
         &(*set)->conf_settlelocal, TRUE, SCIP_DEFAULT_CONF_SETTLELOCAL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/repropagate",
         "should earlier nodes be repropagated in order to replace branching decisions by deductions?",
         &(*set)->conf_repropagate, TRUE, SCIP_DEFAULT_CONF_REPROPAGATE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/keepreprop",
         "should constraints be kept for repropagation even if they are too long?",
         &(*set)->conf_keepreprop, TRUE, SCIP_DEFAULT_CONF_KEEPREPROP,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/separate",
         "should the conflict constraints be separated?",
         &(*set)->conf_separate, TRUE, SCIP_DEFAULT_CONF_SEPARATE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/dynamic",
         "should the conflict constraints be subject to aging?",
         &(*set)->conf_dynamic, TRUE, SCIP_DEFAULT_CONF_DYNAMIC,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/removable",
         "should the conflict's relaxations be subject to LP aging and cleanup?",
         &(*set)->conf_removable, TRUE, SCIP_DEFAULT_CONF_REMOVEABLE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/graph/depthscorefac",
         "score factor for depth level in bound relaxation heuristic",
         &(*set)->conf_depthscorefac, TRUE, SCIP_DEFAULT_CONF_DEPTHSCOREFAC, SCIP_REAL_MIN, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/proofscorefac",
         "score factor for impact on acticity in bound relaxation heuristic",
         &(*set)->conf_proofscorefac, TRUE, SCIP_DEFAULT_CONF_PROOFSCOREFAC, SCIP_REAL_MIN, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/uplockscorefac",
         "score factor for up locks in bound relaxation heuristic",
         &(*set)->conf_uplockscorefac, TRUE, SCIP_DEFAULT_CONF_UPLOCKSCOREFAC, SCIP_REAL_MIN, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/downlockscorefac",
         "score factor for down locks in bound relaxation heuristic",
         &(*set)->conf_downlockscorefac, TRUE, SCIP_DEFAULT_CONF_DOWNLOCKSCOREFAC, SCIP_REAL_MIN, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/scorefac",
         "factor to decrease importance of variables' earlier conflict scores",
         &(*set)->conf_scorefac, TRUE, SCIP_DEFAULT_CONF_SCOREFAC, 1e-6, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/restartnum",
         "number of successful conflict analysis calls that trigger a restart (0: disable conflict restarts)",
         &(*set)->conf_restartnum, FALSE, SCIP_DEFAULT_CONF_RESTARTNUM, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/restartfac",
         "factor to increase restartnum with after each restart",
         &(*set)->conf_restartfac, FALSE, SCIP_DEFAULT_CONF_RESTARTFAC, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/ignorerelaxedbd",
         "should relaxed bounds be ignored?",
         &(*set)->conf_ignorerelaxedbd, TRUE, SCIP_DEFAULT_CONF_IGNORERELAXEDBD,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "conflict/maxvarsdetectimpliedbounds",
         "maximal number of variables to try to detect global bound implications and shorten the whole conflict set (0: disabled)",
         &(*set)->conf_maxvarsdetectimpliedbounds, TRUE, SCIP_DEFAULT_CONF_MAXVARSDETECTIMPLIEDBOUNDS, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/fullshortenconflict",
         "try to shorten the whole conflict set or terminate early (depending on the 'maxvarsdetectimpliedbounds' parameter)",
         &(*set)->conf_fullshortenconflict, TRUE, SCIP_DEFAULT_CONF_FULLSHORTENCONFLICT,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/conflictweight",
         "the weight the VSIDS score is weight by updating the VSIDS for a variable if it is part of a conflict",
         &(*set)->conf_conflictweight, FALSE, SCIP_DEFAULT_CONF_CONFLITWEIGHT, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/conflictgraphweight",
         "the weight the VSIDS score is weight by updating the VSIDS for a variable if it is part of a conflict graph",
         &(*set)->conf_conflictgraphweight, FALSE, SCIP_DEFAULT_CONF_CONFLITGRAPHWEIGHT, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/minimprove",
         "minimal improvement of primal bound to remove conflicts based on a previous incumbent",
         &(*set)->conf_minimprove, TRUE, SCIP_DEFAULT_CONF_MINIMPROVE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/weightsize",
         "weight of the size of a conflict used in score calculation",
         &(*set)->conf_weightsize, TRUE, SCIP_DEFAULT_CONF_WEIGHTSIZE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/weightrepropdepth",
         "weight of the repropagation depth of a conflict used in score calculation",
         &(*set)->conf_weightrepropdepth, TRUE, SCIP_DEFAULT_CONF_WEIGHTREPROPDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "conflict/weightvaliddepth",
         "weight of the valid depth of a conflict used in score calculation",
         &(*set)->conf_weightvaliddepth, TRUE, SCIP_DEFAULT_CONF_WEIGHTVALIDDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "conflict/sepaaltproofs",
         "apply cut generating functions to construct alternative proofs",
         &(*set)->conf_sepaaltproofs, FALSE, SCIP_DEFAULT_CONF_SEPAALTPROOFS,
         NULL, NULL) );

   /* constraint parameters */
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "constraints/agelimit",
         "maximum age an unnecessary constraint can reach before it is deleted (0: dynamic, -1: keep all constraints)",
         &(*set)->cons_agelimit, TRUE, SCIP_DEFAULT_CONS_AGELIMIT, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "constraints/obsoleteage",
         "age of a constraint after which it is marked obsolete (0: dynamic, -1 do not mark constraints obsolete)",
         &(*set)->cons_obsoleteage, TRUE, SCIP_DEFAULT_CONS_OBSOLETEAGE, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "constraints/disableenfops",
         "should enforcement of pseudo solution be disabled?",
         &(*set)->cons_disableenfops, TRUE, SCIP_DEFAULT_CONS_DISABLEENFOPS,
         NULL, NULL) );

   /* display parameters */
   assert(sizeof(int) == sizeof(SCIP_VERBLEVEL));
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "display/verblevel",
         "verbosity level of output",
         (int*)&(*set)->disp_verblevel, FALSE, (int)SCIP_DEFAULT_DISP_VERBLEVEL,
         (int)SCIP_VERBLEVEL_NONE, (int)SCIP_VERBLEVEL_FULL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "display/width",
         "maximal number of characters in a node information line",
         &(*set)->disp_width, FALSE, SCIP_DEFAULT_DISP_WIDTH, 0, INT_MAX,
         SCIPparamChgdDispWidth, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "display/freq",
         "frequency for displaying node information lines",
         &(*set)->disp_freq, FALSE, SCIP_DEFAULT_DISP_FREQ, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "display/headerfreq",
         "frequency for displaying header lines (every n'th node information line)",
         &(*set)->disp_headerfreq, FALSE, SCIP_DEFAULT_DISP_HEADERFREQ, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "display/lpinfo",
         "should the LP solver display status messages?",
         &(*set)->disp_lpinfo, FALSE, SCIP_DEFAULT_DISP_LPINFO,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "display/allviols",
         "display all violations for a given start solution / the best solution after the solving process?",
         &(*set)->disp_allviols, FALSE, SCIP_DEFAULT_DISP_ALLVIOLS,
         NULL, NULL) );

   /* history parameters */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "history/valuebased",
         "should statistics be collected for variable domain value pairs?",
         &(*set)->history_valuebased, FALSE, SCIP_DEFAULT_HISTORY_VALUEBASED,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "history/allowmerge",
         "should variable histories be merged from sub-SCIPs whenever possible?",
         &(*set)->history_allowmerge, FALSE, SCIP_DEFAULT_HISTORY_ALLOWMERGE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "history/allowtransfer",
         "should variable histories be transferred to initialize SCIP copies?",
         &(*set)->history_allowtransfer, FALSE, SCIP_DEFAULT_HISTORY_ALLOWTRANSFER,
         NULL, NULL) );

   /* limit parameters */
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "limits/time",
         "maximal time in seconds to run",
         &(*set)->limit_time, FALSE, SCIP_DEFAULT_LIMIT_TIME, 0.0, SCIP_DEFAULT_LIMIT_TIME,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddLongintParam(*set, messagehdlr, blkmem,
         "limits/nodes",
         "maximal number of nodes to process (-1: no limit)",
         &(*set)->limit_nodes, FALSE, SCIP_DEFAULT_LIMIT_NODES, -1LL, SCIP_LONGINT_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddLongintParam(*set, messagehdlr, blkmem,
         "limits/totalnodes",
         "maximal number of total nodes (incl. restarts) to process (-1: no limit)",
         &(*set)->limit_totalnodes, FALSE, SCIP_DEFAULT_LIMIT_NODES, -1LL, SCIP_LONGINT_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddLongintParam(*set, messagehdlr, blkmem,
         "limits/stallnodes",
         "solving stops, if the given number of nodes was processed since the last improvement of the primal solution value (-1: no limit)",
         &(*set)->limit_stallnodes, FALSE, SCIP_DEFAULT_LIMIT_STALLNODES, -1LL, SCIP_LONGINT_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "limits/memory",
         "maximal memory usage in MB; reported memory usage is lower than real memory usage!",
         &(*set)->limit_memory, FALSE, SCIP_DEFAULT_LIMIT_MEMORY, 0.0, SCIP_MEM_NOLIMIT,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "limits/gap",
         "solving stops, if the relative gap = |primal - dual|/MIN(|dual|,|primal|) is below the given value",
         &(*set)->limit_gap, FALSE, SCIP_DEFAULT_LIMIT_GAP, 0.0, SCIP_REAL_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "limits/absgap",
         "solving stops, if the absolute gap = |primalbound - dualbound| is below the given value",
         &(*set)->limit_absgap, FALSE, SCIP_DEFAULT_LIMIT_ABSGAP, 0.0, SCIP_REAL_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "limits/solutions",
         "solving stops, if the given number of solutions were found (-1: no limit)",
         &(*set)->limit_solutions, FALSE, SCIP_DEFAULT_LIMIT_SOLUTIONS, -1, INT_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "limits/bestsol",
         "solving stops, if the given number of solution improvements were found (-1: no limit)",
         &(*set)->limit_bestsol, FALSE, SCIP_DEFAULT_LIMIT_BESTSOL, -1, INT_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "limits/maxsol",
         "maximal number of solutions to store in the solution storage",
         &(*set)->limit_maxsol, FALSE, SCIP_DEFAULT_LIMIT_MAXSOL, 1, INT_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "limits/maxorigsol",
         "maximal number of solutions candidates to store in the solution storage of the original problem",
         &(*set)->limit_maxorigsol, FALSE, SCIP_DEFAULT_LIMIT_MAXORIGSOL, 0, INT_MAX,
         SCIPparamChgdLimit, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "limits/restarts",
         "solving stops, if the given number of restarts was triggered (-1: no limit)",
         &(*set)->limit_restarts, FALSE, SCIP_DEFAULT_LIMIT_RESTARTS, -1, INT_MAX,
         SCIPparamChgdLimit, NULL) );

   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "limits/autorestartnodes",
         "if solve exceeds this number of nodes for the first time, an automatic restart is triggered (-1: no automatic restart)",
         &(*set)->limit_autorestartnodes, FALSE, SCIP_DEFAULT_LIMIT_AUTORESTARTNODES, -1, INT_MAX,
         SCIPparamChgdLimit, NULL) );

   /* LP parameters */
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/solvefreq",
         "frequency for solving LP at the nodes (-1: never; 0: only root LP)",
         &(*set)->lp_solvefreq, FALSE, SCIP_DEFAULT_LP_SOLVEFREQ, -1, SCIP_MAXTREEDEPTH,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddLongintParam(*set, messagehdlr, blkmem,
         "lp/iterlim",
         "iteration limit for each single LP solve (-1: no limit)",
         &(*set)->lp_iterlim, TRUE, SCIP_DEFAULT_LP_ITERLIM, -1LL, SCIP_LONGINT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddLongintParam(*set, messagehdlr, blkmem,
         "lp/rootiterlim",
         "iteration limit for initial root LP solve (-1: no limit)",
         &(*set)->lp_rootiterlim, TRUE, SCIP_DEFAULT_LP_ROOTITERLIM, -1LL, SCIP_LONGINT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/solvedepth",
         "maximal depth for solving LP at the nodes (-1: no depth limit)",
         &(*set)->lp_solvedepth, FALSE, SCIP_DEFAULT_LP_SOLVEDEPTH, -1, SCIP_MAXTREEDEPTH,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "lp/initalgorithm",
         "LP algorithm for solving initial LP relaxations (automatic 's'implex, 'p'rimal simplex, 'd'ual simplex, 'b'arrier, barrier with 'c'rossover)",
         &(*set)->lp_initalgorithm, FALSE, SCIP_DEFAULT_LP_INITALGORITHM, "spdbc",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "lp/resolvealgorithm",
         "LP algorithm for resolving LP relaxations if a starting basis exists (automatic 's'implex, 'p'rimal simplex, 'd'ual simplex, 'b'arrier, barrier with 'c'rossover)",
         &(*set)->lp_resolvealgorithm, FALSE, SCIP_DEFAULT_LP_RESOLVEALGORITHM, "spdbc",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "lp/pricing",
         "LP pricing strategy ('l'pi default, 'a'uto, 'f'ull pricing, 'p'artial, 's'teepest edge pricing, 'q'uickstart steepest edge pricing, 'd'evex pricing)",
         &(*set)->lp_pricing, FALSE, SCIP_DEFAULT_LP_PRICING, "lafpsqd",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/clearinitialprobinglp",
         "should lp state be cleared at the end of probing mode when lp was initially unsolved, e.g., when called right after presolving?",
         &(*set)->lp_clearinitialprobinglp, TRUE, SCIP_DEFAULT_LP_CLEARINITIALPROBINGLP,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/resolverestore",
         "should the LP be resolved to restore the state at start of diving (if FALSE we buffer the solution values)?",
         &(*set)->lp_resolverestore, TRUE, SCIP_DEFAULT_LP_RESOLVERESTORE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/freesolvalbuffers",
         "should the buffers for storing LP solution values during diving be freed at end of diving?",
         &(*set)->lp_freesolvalbuffers, TRUE, SCIP_DEFAULT_LP_FREESOLVALBUFFERS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/colagelimit",
         "maximum age a dynamic column can reach before it is deleted from the LP (-1: don't delete columns due to aging)",
         &(*set)->lp_colagelimit, TRUE, SCIP_DEFAULT_LP_COLAGELIMIT, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/rowagelimit",
         "maximum age a dynamic row can reach before it is deleted from the LP (-1: don't delete rows due to aging)",
         &(*set)->lp_rowagelimit, TRUE, SCIP_DEFAULT_LP_ROWAGELIMIT, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/cleanupcols",
         "should new non-basic columns be removed after LP solving?",
         &(*set)->lp_cleanupcols, TRUE, SCIP_DEFAULT_LP_CLEANUPCOLS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/cleanupcolsroot",
         "should new non-basic columns be removed after root LP solving?",
         &(*set)->lp_cleanupcolsroot, TRUE, SCIP_DEFAULT_LP_CLEANUPCOLSROOT,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/cleanuprows",
         "should new basic rows be removed after LP solving?",
         &(*set)->lp_cleanuprows, TRUE, SCIP_DEFAULT_LP_CLEANUPROWS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/cleanuprowsroot",
         "should new basic rows be removed after root LP solving?",
         &(*set)->lp_cleanuprowsroot, TRUE, SCIP_DEFAULT_LP_CLEANUPROWSROOT,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/checkstability",
         "should LP solver's return status be checked for stability?",
         &(*set)->lp_checkstability, TRUE, SCIP_DEFAULT_LP_CHECKSTABILITY,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "lp/conditionlimit",
         "maximum condition number of LP basis counted as stable (-1.0: no limit)",
         &(*set)->lp_conditionlimit, TRUE, SCIP_DEFAULT_LP_CONDITIONLIMIT, -1.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/checkprimfeas",
         "should LP solutions be checked for primal feasibility, resolving LP when numerical troubles occur?",
         &(*set)->lp_checkprimfeas, TRUE, SCIP_DEFAULT_LP_CHECKPRIMFEAS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/checkdualfeas",
         "should LP solutions be checked for dual feasibility, resolving LP when numerical troubles occur?",
         &(*set)->lp_checkdualfeas, TRUE, SCIP_DEFAULT_LP_CHECKDUALFEAS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/fastmip",
         "which FASTMIP setting of LP solver should be used? 0: off, 1: low",
         &(*set)->lp_fastmip, TRUE, SCIP_DEFAULT_LP_FASTMIP, 0, 1,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/scaling",
         "LP scaling (0: none, 1: normal, 2: aggressive)",
         &(*set)->lp_scaling, TRUE, SCIP_DEFAULT_LP_SCALING, 0, 2,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/presolving",
         "should presolving of LP solver be used?",
         &(*set)->lp_presolving, TRUE, SCIP_DEFAULT_LP_PRESOLVING,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/lexdualalgo",
         "should the lexicographic dual algorithm be used?",
         &(*set)->lp_lexdualalgo, TRUE, SCIP_DEFAULT_LP_LEXDUALALGO,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/lexdualrootonly",
         "should the lexicographic dual algorithm be applied only at the root node",
         &(*set)->lp_lexdualrootonly, TRUE, SCIP_DEFAULT_LP_LEXDUALROOTONLY,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/lexdualmaxrounds",
         "maximum number of rounds in the lexicographic dual algorithm (-1: unbounded)",
         &(*set)->lp_lexdualmaxrounds, TRUE, SCIP_DEFAULT_LP_LEXDUALMAXROUNDS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/lexdualbasic",
         "choose fractional basic variables in lexicographic dual algorithm?",
         &(*set)->lp_lexdualbasic, TRUE, SCIP_DEFAULT_LP_LEXDUALBASIC,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "lp/lexdualstalling",
         "turn on the lex dual algorithm only when stalling?",
         &(*set)->lp_lexdualstalling, TRUE, SCIP_DEFAULT_LP_LEXDUALSTALLING,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/disablecutoff",
         "disable the cutoff bound in the LP solver? (0: enabled, 1: disabled, 2: auto)",
         &(*set)->lp_disablecutoff, TRUE, SCIP_DEFAULT_LP_DISABLECUTOFF,
         0, 2, NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "lp/rowrepswitch",
         "simplex algorithm shall use row representation of the basis if number of rows divided by number of columns exceeds this value (-1.0 to disable row representation)",
         &(*set)->lp_rowrepswitch, TRUE, SCIP_DEFAULT_LP_ROWREPSWITCH, -1.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/threads",
         "number of threads used for solving the LP (0: automatic)",
         &(*set)->lp_threads, TRUE, SCIP_DEFAULT_LP_THREADS, 0, 64,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "lp/resolveiterfac",
         "factor of average LP iterations that is used as LP iteration limit for LP resolve (-1: unlimited)",
         &(*set)->lp_resolveiterfac, TRUE, SCIP_DEFAULT_LP_RESOLVEITERFAC, -1.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/resolveitermin",
         "minimum number of iterations that are allowed for LP resolve",
         &(*set)->lp_resolveitermin, TRUE, SCIP_DEFAULT_LP_RESOLVEITERMIN, 1, INT_MAX,
         NULL, NULL) );

   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/solutionpolishing",
         "LP solution polishing method (0: disabled, 1: only root, 2: always, 3: auto)",
         &(*set)->lp_solutionpolishing, TRUE, SCIP_DEFAULT_LP_SOLUTIONPOLISHING, 0, 3,
         NULL, NULL) );

   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "lp/refactorinterval",
         "LP refactorization interval (0: auto)",
         &(*set)->lp_refactorinterval, TRUE, SCIP_DEFAULT_LP_REFACTORINTERVAL, 0, INT_MAX,
         NULL, NULL) );

   /* NLP parameters */
   SCIP_CALL( SCIPsetAddStringParam(*set, messagehdlr, blkmem,
         "nlp/solver",
         "solver to use for solving NLPs; leave empty to select NLPI with highest priority",
         &(*set)->nlp_solver, FALSE, SCIP_DEFAULT_NLP_SOLVER,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "nlp/disable",
         "should the NLP relaxation be always disabled (also for NLPs/MINLPs)?",
         &(*set)->nlp_disable, FALSE, SCIP_DEFAULT_NLP_DISABLE,
         NULL, NULL) );

   /* memory parameters */
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "memory/savefac",
         "fraction of maximal memory usage resulting in switch to memory saving mode",
         &(*set)->mem_savefac, FALSE, SCIP_DEFAULT_MEM_SAVEFAC, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "memory/arraygrowfac",
         "memory growing factor for dynamically allocated arrays",
         &(*set)->mem_arraygrowfac, TRUE, SCIP_DEFAULT_MEM_ARRAYGROWFAC, 1.0, 10.0,
         paramChgdArraygrowfac, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "memory/arraygrowinit",
         "initial size of dynamically allocated arrays",
         &(*set)->mem_arraygrowinit, TRUE, SCIP_DEFAULT_MEM_ARRAYGROWINIT, 0, INT_MAX,
         paramChgdArraygrowinit, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "memory/treegrowfac",
         "memory growing factor for tree array",
         &(*set)->mem_treegrowfac, TRUE, SCIP_DEFAULT_MEM_TREEGROWFAC, 1.0, 10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "memory/treegrowinit",
         "initial size of tree array",
         &(*set)->mem_treegrowinit, TRUE, SCIP_DEFAULT_MEM_TREEGROWINIT, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "memory/pathgrowfac",
         "memory growing factor for path array",
         &(*set)->mem_pathgrowfac, TRUE, SCIP_DEFAULT_MEM_PATHGROWFAC, 1.0, 10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "memory/pathgrowinit",
         "initial size of path array",
         &(*set)->mem_pathgrowinit, TRUE, SCIP_DEFAULT_MEM_PATHGROWINIT, 0, INT_MAX,
         NULL, NULL) );

   /* miscellaneous parameters */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/catchctrlc",
         "should the CTRL-C interrupt be caught by SCIP?",
         &(*set)->misc_catchctrlc, FALSE, SCIP_DEFAULT_MISC_CATCHCTRLC,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/usevartable",
         "should a hashtable be used to map from variable names to variables?",
         &(*set)->misc_usevartable, FALSE, SCIP_DEFAULT_MISC_USEVARTABLE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/useconstable",
         "should a hashtable be used to map from constraint names to constraints?",
         &(*set)->misc_useconstable, FALSE, SCIP_DEFAULT_MISC_USECONSTABLE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/usesmalltables",
         "should smaller hashtables be used? yields better performance for small problems with about 100 variables",
         &(*set)->misc_usesmalltables, FALSE, SCIP_DEFAULT_MISC_USESMALLTABLES,
         NULL, NULL) );
#if 0 /**@todo activate exactsolve parameter and finish implementation of solving MIPs exactly */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/exactsolve",
         "should the problem be solved exactly (with proven dual bounds)?",
         &(*set)->misc_exactsolve, FALSE, SCIP_DEFAULT_MISC_EXACTSOLVE,
         NULL, NULL) );
#else
   (*set)->misc_exactsolve = SCIP_DEFAULT_MISC_EXACTSOLVE;
#endif

   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/resetstat",
         "should the statistics be reset if the transformed problem is freed (in case of a Benders decomposition this parameter should be set to FALSE)",
         &(*set)->misc_resetstat, FALSE, SCIP_DEFAULT_MISC_RESETSTAT,
         NULL, NULL) );

   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/improvingsols",
         "should only solutions be checked which improve the primal bound",
         &(*set)->misc_improvingsols, FALSE, SCIP_DEFAULT_MISC_IMPROVINGSOLS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/printreason",
         "should the reason be printed if a given start solution is infeasible",
         &(*set)->misc_printreason, FALSE, SCIP_DEFAULT_MISC_PRINTREASON,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/estimexternmem",
         "should the usage of external memory be estimated?",
         &(*set)->misc_estimexternmem, FALSE, SCIP_DEFAULT_MISC_ESTIMEXTERNMEM,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/transorigsols",
         "should SCIP try to transfer original solutions to the transformed space (after presolving)?",
         &(*set)->misc_transorigsols, FALSE, SCIP_DEFAULT_MISC_TRANSORIGSOLS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "misc/transsolsorig",
         "should SCIP try to transfer transformed solutions to the original space (after solving)?",
         &(*set)->misc_transsolsorig, FALSE, SCIP_DEFAULT_MISC_TRANSSOLSORIG,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
            "misc/calcintegral",
            "should SCIP calculate the primal dual integral value?",
            &(*set)->misc_calcintegral, FALSE, SCIP_DEFAULT_MISC_CALCINTEGRAL,
            NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
            "misc/finitesolutionstore",
            "should SCIP try to remove infinite fixings from solutions copied to the solution store?",
            &(*set)->misc_finitesolstore, FALSE, SCIP_DEFAULT_MISC_FINITESOLSTORE,
            NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
            "misc/outputorigsol",
            "should the best solution be transformed to the orignal space and be output in command line run?",
            &(*set)->misc_outputorigsol, FALSE, SCIP_DEFAULT_MISC_OUTPUTORIGSOL,
            NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
            "misc/allowdualreds",
            "should dual reductions in propagation methods and presolver be allowed?",
            &(*set)->misc_allowdualreds, FALSE, SCIP_DEFAULT_MISC_ALLOWDUALREDS,
            NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
            "misc/allowobjprop",
            "should propagation to the current objective be allowed in propagation methods?",
            &(*set)->misc_allowobjprop, FALSE, SCIP_DEFAULT_MISC_ALLOWOBJPROP,
            NULL, NULL) );

   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "misc/referencevalue",
         "objective value for reference purposes",
         &(*set)->misc_referencevalue, FALSE, SCIP_DEFAULT_MISC_REFERENCEVALUE, SCIP_REAL_MIN, SCIP_REAL_MAX,
         NULL, NULL) );

#ifdef WITH_DEBUG_SOLUTION
   SCIP_CALL( SCIPsetAddStringParam(*set, messagehdlr, blkmem,
         "misc/debugsol",
         "path to a debug solution",
         &(*set)->misc_debugsol, FALSE, SCIP_DEFAULT_MISC_DEBUGSOLUTION,
         NULL, NULL) );
#endif

   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "misc/usesymmetry",
         "used symmetry handling technique (0: off; 1: polyhedral; 2: orbital fixing)",
         &(*set)->misc_usesymmetry, FALSE, SCIP_DEFAULT_MISC_USESYMMETRY, 0, 2,
         paramChgdUsesymmetry, NULL) );

   /* randomization parameters */
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "randomization/randomseedshift",
         "global shift of all random seeds in the plugins and the LP random seed",
         &(*set)->random_randomseedshift, FALSE, SCIP_DEFAULT_RANDOM_RANDSEEDSHIFT, 0, INT_MAX,
         NULL, NULL) );

   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "randomization/permutationseed",
         "seed value for permuting the problem after reading/transformation (0: no permutation)",
         &(*set)->random_permutationseed, FALSE, SCIP_DEFAULT_RANDOM_PERMUTATIONSEED, 0, INT_MAX,
         NULL, NULL) );

   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "randomization/permuteconss",
         "should order of constraints be permuted (depends on permutationseed)?",
         &(*set)->random_permuteconss, TRUE, SCIP_DEFAULT_RANDOM_PERMUTECONSS,
         NULL, NULL) );

   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "randomization/permutevars",
         "should order of variables be permuted (depends on permutationseed)?",
         &(*set)->random_permutevars, TRUE, SCIP_DEFAULT_RANDOM_PERMUTEVARS,
         NULL, NULL) );


   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "randomization/lpseed",
         "random seed for LP solver, e.g. for perturbations in the simplex (0: LP default)",
         &(*set)->random_randomseed, FALSE, SCIP_DEFAULT_RANDOM_LPSEED, 0, INT_MAX,
         NULL, NULL) );

   /* node selection */
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "nodeselection/childsel",
         "child selection rule ('d'own, 'u'p, 'p'seudo costs, 'i'nference, 'l'p value, 'r'oot LP value difference, 'h'ybrid inference/root LP value difference)",
         &(*set)->nodesel_childsel, FALSE, SCIP_DEFAULT_NODESEL_CHILDSEL, "dupilrh",
         NULL, NULL) );

   /* numerical parameters */
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/infinity",
         "values larger than this are considered infinity",
         &(*set)->num_infinity, FALSE, SCIP_DEFAULT_INFINITY, 1e+10, SCIP_INVALID/10.0,
         paramChgInfinity, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/epsilon",
         "absolute values smaller than this are considered zero",
         &(*set)->num_epsilon, FALSE, SCIP_DEFAULT_EPSILON, SCIP_MINEPSILON, SCIP_MAXEPSILON,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/sumepsilon",
         "absolute values of sums smaller than this are considered zero",
         &(*set)->num_sumepsilon, FALSE, SCIP_DEFAULT_SUMEPSILON, SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/feastol",
         "feasibility tolerance for constraints",
         &(*set)->num_feastol, FALSE, SCIP_DEFAULT_FEASTOL, SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON,
         paramChgdFeastol, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/checkfeastolfac",
         "feasibility tolerance factor; for checking the feasibility of the best solution",
         &(*set)->num_checkfeastolfac, FALSE, SCIP_DEFAULT_CHECKFEASTOLFAC, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/lpfeastol",
         "primal feasibility tolerance of LP solver",
         &(*set)->num_lpfeastol, FALSE, SCIP_DEFAULT_LPFEASTOL, SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON,
         paramChgdLpfeastol, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/dualfeastol",
         "feasibility tolerance for reduced costs in LP solution",
         &(*set)->num_dualfeastol, FALSE, SCIP_DEFAULT_DUALFEASTOL, SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON,
         paramChgdDualfeastol, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/barrierconvtol",
         "LP convergence tolerance used in barrier algorithm",
         &(*set)->num_barrierconvtol, TRUE, SCIP_DEFAULT_BARRIERCONVTOL, SCIP_MINEPSILON*1e+03, SCIP_MAXEPSILON,
         paramChgdBarrierconvtol, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/boundstreps",
         "minimal relative improve for strengthening bounds",
         &(*set)->num_boundstreps, TRUE, SCIP_DEFAULT_BOUNDSTREPS, SCIP_MINEPSILON*1e+03, SCIP_INVALID/10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/pseudocosteps",
         "minimal variable distance value to use for branching pseudo cost updates",
         &(*set)->num_pseudocosteps, TRUE, SCIP_DEFAULT_PSEUDOCOSTEPS, SCIP_MINEPSILON*1e+03, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/pseudocostdelta",
         "minimal objective distance value to use for branching pseudo cost updates",
         &(*set)->num_pseudocostdelta, TRUE, SCIP_DEFAULT_PSEUDOCOSTDELTA, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/recomputefac",
         "minimal decrease factor that causes the recomputation of a value (e.g., pseudo objective) instead of an update",
         &(*set)->num_recompfac, TRUE, SCIP_DEFAULT_RECOMPFAC, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "numerics/hugeval",
         "values larger than this are considered huge and should be handled separately (e.g., in activity computation)",
         &(*set)->num_hugeval, TRUE, SCIP_DEFAULT_HUGEVAL, 0.0, SCIP_INVALID/10.0,
         NULL, NULL) );

   /* presolving parameters */
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "presolving/maxrounds",
         "maximal number of presolving rounds (-1: unlimited, 0: off)",
         &(*set)->presol_maxrounds, FALSE, SCIP_DEFAULT_PRESOL_MAXROUNDS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "presolving/abortfac",
         "abort presolve, if at most this fraction of the problem was changed in last presolve round",
         &(*set)->presol_abortfac, TRUE, SCIP_DEFAULT_PRESOL_ABORTFAC, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "presolving/maxrestarts",
         "maximal number of restarts (-1: unlimited)",
         &(*set)->presol_maxrestarts, FALSE, SCIP_DEFAULT_PRESOL_MAXRESTARTS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "presolving/restartfac",
         "fraction of integer variables that were fixed in the root node triggering a restart with preprocessing after root node evaluation",
         &(*set)->presol_restartfac, TRUE, SCIP_DEFAULT_PRESOL_RESTARTFAC, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "presolving/immrestartfac",
         "fraction of integer variables that were fixed in the root node triggering an immediate restart with preprocessing",
         &(*set)->presol_immrestartfac, TRUE, SCIP_DEFAULT_PRESOL_IMMRESTARTFAC, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "presolving/subrestartfac",
         "fraction of integer variables that were globally fixed during the solving process triggering a restart with preprocessing",
         &(*set)->presol_subrestartfac, TRUE, SCIP_DEFAULT_PRESOL_SUBRESTARTFAC, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "presolving/restartminred",
         "minimal fraction of integer variables removed after restart to allow for an additional restart",
         &(*set)->presol_restartminred, TRUE, SCIP_DEFAULT_PRESOL_RESTARTMINRED, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "presolving/donotmultaggr",
         "should multi-aggregation of variables be forbidden?",
         &(*set)->presol_donotmultaggr, TRUE, SCIP_DEFAULT_PRESOL_DONOTMULTAGGR,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "presolving/donotaggr",
         "should aggregation of variables be forbidden?",
         &(*set)->presol_donotaggr, TRUE, SCIP_DEFAULT_PRESOL_DONOTAGGR,
         NULL, NULL) );


   /* pricing parameters */
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "pricing/maxvars",
         "maximal number of variables priced in per pricing round",
         &(*set)->price_maxvars, FALSE, SCIP_DEFAULT_PRICE_MAXVARS, 1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "pricing/maxvarsroot",
         "maximal number of priced variables at the root node",
         &(*set)->price_maxvarsroot, FALSE, SCIP_DEFAULT_PRICE_MAXVARSROOT, 1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "pricing/abortfac",
         "pricing is aborted, if fac * pricing/maxvars pricing candidates were found",
         &(*set)->price_abortfac, FALSE, SCIP_DEFAULT_PRICE_ABORTFAC, 1.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "pricing/delvars",
         "should variables created at the current node be deleted when the node is solved in case they are not present in the LP anymore?",
         &(*set)->price_delvars, FALSE, SCIP_DEFAULT_PRICE_DELVARS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "pricing/delvarsroot",
         "should variables created at the root node be deleted when the root is solved in case they are not present in the LP anymore?",
         &(*set)->price_delvarsroot, FALSE, SCIP_DEFAULT_PRICE_DELVARSROOT,
         NULL, NULL) );

   /* propagation parameters */
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "propagating/maxrounds",
         "maximal number of propagation rounds per node (-1: unlimited)",
         &(*set)->prop_maxrounds, FALSE, SCIP_DEFAULT_PROP_MAXROUNDS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "propagating/maxroundsroot",
         "maximal number of propagation rounds in the root node (-1: unlimited)",
         &(*set)->prop_maxroundsroot, FALSE, SCIP_DEFAULT_PROP_MAXROUNDSROOT, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "propagating/abortoncutoff",
         "should propagation be aborted immediately? setting this to FALSE could help conflict analysis to produce more conflict constraints",
         &(*set)->prop_abortoncutoff, FALSE, SCIP_DEFAULT_PROP_ABORTONCUTOFF,
         NULL, NULL) );

   /* reoptimization */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/enable",
         "should reoptimization used?",
         &(*set)->reopt_enable, FALSE, SCIP_DEFAULT_REOPT_ENABLE,
         paramChgdEnableReopt, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "reoptimization/maxsavednodes",
         "maximal number of saved nodes",
         &(*set)->reopt_maxsavednodes, TRUE, SCIP_DEFAULT_REOPT_MAXSAVEDNODES, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "reoptimization/maxdiffofnodes",
         "maximal number of bound changes between two stored nodes on one path",
         &(*set)->reopt_maxdiffofnodes, TRUE, SCIP_DEFAULT_REOPT_MAXDIFFOFNODES, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/globalcons/sepainfsubtrees",
         "save global constraints to separate infeasible subtrees.",
         &(*set)->reopt_sepaglbinfsubtrees, FALSE, SCIP_DEFAULT_REOPT_SEPAGLBINFSUBTREES,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/sepabestsol",
         "separate the optimal solution, i.e., for constrained shortest path",
         &(*set)->reopt_sepabestsol, TRUE, SCIP_DEFAULT_REOPT_SEPABESTSOL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/storevarhistory",
         "use variable history of the previous solve if the objctive function has changed only slightly",
         &(*set)->reopt_storevarhistory, TRUE, SCIP_DEFAULT_REOPT_STOREVARHISTOTY,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/usepscost",
         "re-use pseudo costs if the objective function changed only slightly ",
         &(*set)->reopt_usepscost, TRUE, SCIP_DEFAULT_REOPT_USEPSCOST,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "reoptimization/solvelp",
         "at which reopttype should the LP be solved? (1: transit, 3: strong branched, 4: w/ added logicor, 5: only leafs).",
         &(*set)->reopt_solvelp, TRUE, SCIP_DEFAULT_REOPT_SOLVELP, 1, 5,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "reoptimization/solvelpdiff",
         "maximal number of bound changes at node to skip solving the LP",
         &(*set)->reopt_solvelpdiff, TRUE, SCIP_DEFAULT_REOPT_SOLVELPDIFF, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "reoptimization/savesols",
         "number of best solutions which should be saved for the following runs. (-1: save all)",
         &(*set)->reopt_savesols, TRUE, SCIP_DEFAULT_REOPT_SAVESOLS, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "reoptimization/objsimrootLP",
         "similarity of two sequential objective function to disable solving the root LP.",
         &(*set)->reopt_objsimrootlp, TRUE, SCIP_DEFAULT_REOPT_OBJSIMROOTLP, -1.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "reoptimization/objsimsol",
         "similarity of two objective functions to re-use stored solutions",
         &(*set)->reopt_objsimsol, TRUE, SCIP_DEFAULT_REOPT_OBJSIMSOL, -1.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "reoptimization/delay",
         "minimum similarity for using reoptimization of the search tree.",
         &(*set)->reopt_objsimdelay, TRUE, SCIP_DEFAULT_REOPT_OBJSIMDELAY, -1.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/commontimelimit",
         "time limit over all reoptimization rounds?.",
         &(*set)->reopt_commontimelimit, TRUE, SCIP_DEFAULT_REOPT_COMMONTIMELIMIT,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/shrinkinner",
         "replace branched inner nodes by their child nodes, if the number of bound changes is not to large",
         &(*set)->reopt_shrinkinner, TRUE, SCIP_DEFAULT_REOPT_SHRINKINNER,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/strongbranchinginit",
         "try to fix variables at the root node before reoptimizing by probing like strong branching",
         &(*set)->reopt_sbinit, TRUE, SCIP_DEFAULT_REOPT_STRONGBRANCHINIT,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/reducetofrontier",
         "delete stored nodes which were not reoptimized",
         &(*set)->reopt_reducetofrontier, TRUE, SCIP_DEFAULT_REOPT_REDUCETOFRONTIER,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "reoptimization/forceheurrestart",
         "force a restart if the last n optimal solutions were found by heuristic reoptsols",
         &(*set)->reopt_forceheurrestart, TRUE, SCIP_DEFAULT_REOPT_FORCEHEURRESTART, 1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/saveconsprop",
         "save constraint propagations",
         &(*set)->reopt_saveconsprop, TRUE, SCIP_DEFAULT_REOPT_SAVECONSPROP,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/usesplitcons", "use constraints to reconstruct the subtree pruned be dual reduction when reactivating the node",
         &(*set)->reopt_usesplitcons, TRUE, SCIP_DEFAULT_REOPT_USESPLITCONS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "reoptimization/varorderinterdiction", "use 'd'efault, 'r'andom or a variable ordering based on 'i'nference score for interdiction branching used during reoptimization",
         &(*set)->reopt_varorderinterdiction, TRUE, SCIP_DEFAULT_REOPT_VARORDERINTERDICTION, "dir",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reoptimization/usecuts",
         "reoptimize cuts found at the root node",
         &(*set)->reopt_usecuts, TRUE, SCIP_DEFAULT_REOPT_USECUTS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "reoptimization/maxcutage",
         "maximal age of a cut to be use for reoptimization",
         &(*set)->reopt_maxcutage, TRUE, SCIP_DEFAULT_REOPT_MAXCUTAGE, 0, INT_MAX,
         NULL, NULL) );

   /* separation parameters */
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/maxbounddist",
         "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying separation (0.0: only on current best node, 1.0: on all nodes)",
         &(*set)->sepa_maxbounddist, FALSE, SCIP_DEFAULT_SEPA_MAXBOUNDDIST, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/maxlocalbounddist",
         "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying local separation (0.0: only on current best node, 1.0: on all nodes)",
         &(*set)->sepa_maxlocalbounddist, FALSE, SCIP_DEFAULT_SEPA_MAXLOCALBOUNDDIST, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/maxcoefratio",
         "maximal ratio between coefficients in strongcg, cmir, and flowcover cuts",
         &(*set)->sepa_maxcoefratio, FALSE, SCIP_DEFAULT_SEPA_MAXCOEFRATIO, 1.0, SCIP_INVALID/10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/minefficacy",
         "minimal efficacy for a cut to enter the LP",
         &(*set)->sepa_minefficacy, FALSE, SCIP_DEFAULT_SEPA_MINEFFICACY, 0.0, SCIP_INVALID/10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/minefficacyroot",
         "minimal efficacy for a cut to enter the LP in the root node",
         &(*set)->sepa_minefficacyroot, FALSE, SCIP_DEFAULT_SEPA_MINEFFICACYROOT, 0.0, SCIP_INVALID/10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/minortho",
         "minimal orthogonality for a cut to enter the LP",
         &(*set)->sepa_minortho, FALSE, SCIP_DEFAULT_SEPA_MINORTHO, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/minorthoroot",
         "minimal orthogonality for a cut to enter the LP in the root node",
         &(*set)->sepa_minorthoroot, FALSE, SCIP_DEFAULT_SEPA_MINORTHOROOT, 0.0, 1.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/objparalfac",
         "factor to scale objective parallelism of cut in separation score calculation",
         &(*set)->sepa_objparalfac, TRUE, SCIP_DEFAULT_SEPA_OBJPARALFAC, 0.0, SCIP_INVALID/10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "separating/intsupportfac",
         "factor to scale integral support of cut in separation score calculation",
         &(*set)->sepa_intsupportfac, TRUE, SCIP_DEFAULT_SEPA_INTSUPPORTFAC, 0.0, SCIP_INVALID/10.0,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
           "separating/minactivityquot",
           "minimum cut activity quotient to convert cuts into constraints during a restart (0.0: all cuts are converted)",
           &(*set)->sepa_minactivityquot, FALSE, SCIP_DEFAULT_SEPA_MINACTIVITYQUOT, 0.0, 1.0,
           NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "separating/orthofunc",
         "function used for calc. scalar prod. in orthogonality test ('e'uclidean, 'd'iscrete)",
         &(*set)->sepa_orthofunc, TRUE, SCIP_DEFAULT_SEPA_ORTHOFUNC, "ed",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "separating/efficacynorm",
         "row norm to use for efficacy calculation ('e'uclidean, 'm'aximum, 's'um, 'd'iscrete)",
         &(*set)->sepa_efficacynorm, TRUE, SCIP_DEFAULT_SEPA_EFFICACYNORM, "emsd",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "separating/cutselrestart",
         "cut selection during restart ('a'ge, activity 'q'uotient)",
         &(*set)->sepa_cutselrestart, TRUE, SCIP_DEFAULT_SEPA_CUTSELRESTART, "aq",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddCharParam(*set, messagehdlr, blkmem,
         "separating/cutselsubscip",
         "cut selection for sub SCIPs  ('a'ge, activity 'q'uotient)",
         &(*set)->sepa_cutselsubscip, TRUE, SCIP_DEFAULT_SEPA_CUTSELSUBSCIP, "aq",
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxruns",
         "maximal number of runs for which separation is enabled (-1: unlimited)",
         &(*set)->sepa_maxruns, TRUE, SCIP_DEFAULT_SEPA_MAXRUNS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &(*set)->sepa_maxrounds, FALSE, SCIP_DEFAULT_SEPA_MAXROUNDS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxroundsroot",
         "maximal number of separation rounds in the root node (-1: unlimited)",
         &(*set)->sepa_maxroundsroot, FALSE, SCIP_DEFAULT_SEPA_MAXROUNDSROOT, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxroundsrootsubrun",
         "maximal number of separation rounds in the root node of a subsequent run (-1: unlimited)",
         &(*set)->sepa_maxroundsrootsubrun, TRUE, SCIP_DEFAULT_SEPA_MAXROUNDSROOTSUBRUN, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxaddrounds",
         "maximal additional number of separation rounds in subsequent price-and-cut loops (-1: no additional restriction)",
         &(*set)->sepa_maxaddrounds, TRUE, SCIP_DEFAULT_SEPA_MAXADDROUNDS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxstallrounds",
         "maximal number of consecutive separation rounds without objective or integrality improvement in local nodes (-1: no additional restriction)",
         &(*set)->sepa_maxstallrounds, FALSE, SCIP_DEFAULT_SEPA_MAXSTALLROUNDS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxstallroundsroot",
         "maximal number of consecutive separation rounds without objective or integrality improvement in the root node (-1: no additional restriction)",
         &(*set)->sepa_maxstallroundsroot, FALSE, SCIP_DEFAULT_SEPA_MAXSTALLROUNDSROOT, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxincrounds",
         "maximal number of consecutive separation rounds that increase the size of the LP relaxation per node (-1: unlimited)",
         &(*set)->sepa_maxincrounds, FALSE, SCIP_DEFAULT_SEPA_MAXINCROUNDS, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxcuts",
         "maximal number of cuts separated per separation round (0: disable local separation)",
         &(*set)->sepa_maxcuts, FALSE, SCIP_DEFAULT_SEPA_MAXCUTS, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/maxcutsroot",
         "maximal number of separated cuts at the root node (0: disable root node separation)",
         &(*set)->sepa_maxcutsroot, FALSE, SCIP_DEFAULT_SEPA_MAXCUTSROOT, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/cutagelimit",
         "maximum age a cut can reach before it is deleted from the global cut pool, or -1 to keep all cuts",
         &(*set)->sepa_cutagelimit, TRUE, SCIP_DEFAULT_SEPA_CUTAGELIMIT, -1, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "separating/poolfreq",
         "separation frequency for the global cut pool (-1: disable global cut pool, 0: only separate pool at the root)",
         &(*set)->sepa_poolfreq, FALSE, SCIP_DEFAULT_SEPA_POOLFREQ, -1, SCIP_MAXTREEDEPTH,
         NULL, NULL) );

   /* parallel parameters */
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "parallel/mode",
         "parallel optimisation mode, 0: opportunistic or 1: deterministic.",
         &(*set)->parallel_mode, FALSE, SCIP_DEFAULT_PARALLEL_MODE, 0, 1,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "parallel/minnthreads",
         "the minimum number of threads used during parallel solve",
         &(*set)->parallel_minnthreads, FALSE, SCIP_DEFAULT_PARALLEL_MINNTHREADS, 0, 64,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "parallel/maxnthreads",
         "the maximum number of threads used during parallel solve",
         &(*set)->parallel_maxnthreads, FALSE, SCIP_DEFAULT_PARALLEL_MAXNTHREADS, 0, 64,
         NULL, NULL) );

   /* concurrent solver parameters */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "concurrent/changeseeds",
         "set different random seeds in each concurrent solver?",
         &(*set)->concurrent_changeseeds, FALSE, SCIP_DEFAULT_CONCURRENT_CHANGESEEDS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "concurrent/changechildsel",
         "use different child selection rules in each concurrent solver?",
         &(*set)->concurrent_changechildsel, FALSE, SCIP_DEFAULT_CONCURRENT_CHANGECHILDSEL,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "concurrent/commvarbnds",
         "should the concurrent solvers communicate global variable bound changes?",
         &(*set)->concurrent_commvarbnds, FALSE, SCIP_DEFAULT_CONCURRENT_COMMVARBNDS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "concurrent/presolvebefore",
         "should the problem be presolved before it is copied to the concurrent solvers?",
         &(*set)->concurrent_presolvebefore, FALSE, SCIP_DEFAULT_CONCURRENT_PRESOLVEBEFORE,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "concurrent/initseed",
         "maximum number of solutions that will be shared in a one synchronization",
         &(*set)->concurrent_initseed, FALSE, SCIP_DEFAULT_CONCURRENT_INITSEED, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "concurrent/sync/freqinit",
         "initial frequency of synchronization with other threads",
         &(*set)->concurrent_freqinit, FALSE, SCIP_DEFAULT_CONCURRENT_FREQINIT, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
      SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "concurrent/sync/freqmax",
         "maximal frequency of synchronization with other threads",
         &(*set)->concurrent_freqmax, FALSE, SCIP_DEFAULT_CONCURRENT_FREQMAX, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "concurrent/sync/freqfactor",
         "factor by which the frequency of synchronization is changed",
         &(*set)->concurrent_freqfactor, FALSE, SCIP_DEFAULT_CONCURRENT_FREQFACTOR, 1.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "concurrent/sync/targetprogress",
         "when adapting the synchronization frequency this value is the targeted relative difference by which the absolute gap decreases per synchronization",
         &(*set)->concurrent_targetprogress, FALSE, SCIP_DEFAULT_CONCURRENT_TARGETPROGRESS, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "concurrent/sync/maxnsols",
         "maximum number of solutions that will be shared in a single synchronization",
         &(*set)->concurrent_maxnsols, FALSE, SCIP_DEFAULT_CONCURRENT_MAXNSOLS, 0, 1000,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "concurrent/sync/maxnsyncdelay",
         "maximum number of synchronizations before reading is enforced regardless of delay",
         &(*set)->concurrent_maxnsyncdelay, TRUE, SCIP_DEFAULT_CONCURRENT_MAXNSYNCDELAY, 0, 100,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddRealParam(*set, messagehdlr, blkmem,
         "concurrent/sync/minsyncdelay",
         "minimum delay before synchronization data is read",
         &(*set)->concurrent_minsyncdelay, FALSE, SCIP_DEFAULT_CONCURRENT_MINSYNCDELAY, 0.0, SCIP_REAL_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "concurrent/sync/nbestsols",
         "how many of the N best solutions should be considered for synchronization?",
         &(*set)->concurrent_nbestsols, FALSE, SCIP_DEFAULT_CONCURRENT_NBESTSOLS, 0, INT_MAX,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddStringParam(*set, messagehdlr, blkmem,
         "concurrent/paramsetprefix",
         "path prefix for parameter setting files of concurrent solvers",
         &(*set)->concurrent_paramsetprefix, FALSE, SCIP_DEFAULT_CONCURRENT_PARAMSETPREFIX,
         NULL, NULL) );

   /* timing parameters */
   assert(sizeof(int) == sizeof(SCIP_CLOCKTYPE));
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "timing/clocktype",
         "default clock type (1: CPU user seconds, 2: wall clock time)",
         (int*)&(*set)->time_clocktype, FALSE, (int)SCIP_DEFAULT_TIME_CLOCKTYPE, 1, 2,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "timing/enabled",
         "is timing enabled?",
         &(*set)->time_enabled, FALSE, SCIP_DEFAULT_TIME_ENABLED,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "timing/reading",
         "belongs reading time to solving time?",
         &(*set)->time_reading, FALSE, SCIP_DEFAULT_TIME_READING,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "timing/rareclockcheck",
         "should clock checks of solving time be performed less frequently (note: time limit could be exceeded slightly)",
         &(*set)->time_rareclockcheck, FALSE, SCIP_DEFAULT_TIME_RARECLOCKCHECK,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "timing/statistictiming",
         "should timing for statistic output be performed?",
         &(*set)->time_statistictiming, FALSE, SCIP_DEFAULT_TIME_STATISTICTIMING,
         paramChgdStatistictiming, NULL) );

   /* visualization parameters */
   SCIP_CALL( SCIPsetAddStringParam(*set, messagehdlr, blkmem,
         "visual/vbcfilename",
         "name of the VBC tool output file, or - if no VBC tool output should be created",
         &(*set)->visual_vbcfilename, FALSE, SCIP_DEFAULT_VISUAL_VBCFILENAME,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddStringParam(*set, messagehdlr, blkmem,
         "visual/bakfilename",
         "name of the BAK tool output file, or - if no BAK tool output should be created",
         &(*set)->visual_bakfilename, FALSE, SCIP_DEFAULT_VISUAL_BAKFILENAME,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "visual/realtime",
         "should the real solving time be used instead of a time step counter in visualization?",
         &(*set)->visual_realtime, FALSE, SCIP_DEFAULT_VISUAL_REALTIME,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "visual/dispsols",
         "should the node where solutions are found be visualized?",
         &(*set)->visual_dispsols, FALSE, SCIP_DEFAULT_VISUAL_DISPSOLS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "visual/objextern",
         "should be output the external value of the objective?",
         &(*set)->visual_objextern, FALSE, SCIP_DEFAULT_VISUAL_OBJEXTERN,
         NULL, NULL) );

   /* Reading parameters */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reading/initialconss",
         "should model constraints be marked as initial?",
         &(*set)->read_initialconss, FALSE, SCIP_DEFAULT_READ_INITIALCONSS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reading/dynamicconss",
         "should model constraints be subject to aging?",
         &(*set)->read_dynamicconss, FALSE, SCIP_DEFAULT_READ_DYNAMICCONSS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reading/dynamiccols",
         "should columns be added and removed dynamically to the LP?",
         &(*set)->read_dynamiccols, FALSE, SCIP_DEFAULT_READ_DYNAMICCOLS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "reading/dynamicrows",
         "should rows be added and removed dynamically to the LP?",
         &(*set)->read_dynamicrows, FALSE, SCIP_DEFAULT_READ_DYNAMICROWS,
         NULL, NULL) );

   /* Writing parameters */
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "write/allconss",
         "should all constraints be written (including the redundant constraints)?",
         &(*set)->write_allconss, FALSE, SCIP_DEFAULT_WRITE_ALLCONSS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddBoolParam(*set, messagehdlr, blkmem,
         "write/printzeros",
         "should variables set to zero be printed?",
         &(*set)->write_printzeros, FALSE, SCIP_DEFAULT_PRINTZEROS,
         NULL, NULL) );
   SCIP_CALL( SCIPsetAddIntParam(*set, messagehdlr, blkmem,
         "write/genericnamesoffset",
         "when writing a generic problem the index for the first variable should start with?",
         &(*set)->write_genoffset, FALSE, SCIP_DEFAULT_WRITE_GENNAMES_OFFSET, 0, INT_MAX/2,
         NULL, NULL) );

   return SCIP_OKAY;
}

/** frees global SCIP settings */
SCIP_RETCODE SCIPsetFree(
   SCIP_SET**            set,                /**< pointer to SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(set != NULL);

   if( *set == NULL )
      return SCIP_OKAY;

   /* free parameter set */
   SCIPparamsetFree(&(*set)->paramset, blkmem);

   /* free file readers */
   for( i = 0; i < (*set)->nreaders; ++i )
   {
      SCIP_CALL( SCIPreaderFree(&(*set)->readers[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->readers);

   /* free variable pricers */
   for( i = 0; i < (*set)->npricers; ++i )
   {
      SCIP_CALL( SCIPpricerFree(&(*set)->pricers[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->pricers);

   /* free constraint handlers */
   for( i = 0; i < (*set)->nconshdlrs; ++i )
   {
      SCIP_CALL( SCIPconshdlrFree(&(*set)->conshdlrs[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->conshdlrs);
   BMSfreeMemoryArrayNull(&(*set)->conshdlrs_sepa);
   BMSfreeMemoryArrayNull(&(*set)->conshdlrs_enfo);
   BMSfreeMemoryArrayNull(&(*set)->conshdlrs_include);

   /* free conflict handlers */
   for( i = 0; i < (*set)->nconflicthdlrs; ++i )
   {
      SCIP_CALL( SCIPconflicthdlrFree(&(*set)->conflicthdlrs[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->conflicthdlrs);

   /* free presolvers */
   for( i = 0; i < (*set)->npresols; ++i )
   {
      SCIP_CALL( SCIPpresolFree(&(*set)->presols[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->presols);

   /* free relaxators */
   for( i = 0; i < (*set)->nrelaxs; ++i )
   {
      SCIP_CALL( SCIPrelaxFree(&(*set)->relaxs[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->relaxs);

   /* free separators */
   for( i = 0; i < (*set)->nsepas; ++i )
   {
      SCIP_CALL( SCIPsepaFree(&(*set)->sepas[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->sepas);

   /* free propagators */
   for( i = 0; i < (*set)->nprops; ++i )
   {
      SCIP_CALL( SCIPpropFree(&(*set)->props[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->props);
   BMSfreeMemoryArrayNull(&(*set)->props_presol);

   /* free primal heuristics */
   for( i = 0; i < (*set)->nheurs; ++i )
   {
      SCIP_CALL( SCIPheurFree(&(*set)->heurs[i], *set, blkmem) );
   }
   BMSfreeMemoryArrayNull(&(*set)->heurs);

   /* free tree compressions */
   for( i = 0; i < (*set)->ncomprs; ++i )
   {
      SCIP_CALL( SCIPcomprFree(&(*set)->comprs[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->comprs);

   /* free event handlers */
   for( i = 0; i < (*set)->neventhdlrs; ++i )
   {
      SCIP_CALL( SCIPeventhdlrFree(&(*set)->eventhdlrs[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->eventhdlrs);

   /* free node selectors */
   for( i = 0; i < (*set)->nnodesels; ++i )
   {
      SCIP_CALL( SCIPnodeselFree(&(*set)->nodesels[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->nodesels);

   /* free branching methods */
   for( i = 0; i < (*set)->nbranchrules; ++i )
   {
      SCIP_CALL( SCIPbranchruleFree(&(*set)->branchrules[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->branchrules);

   /* free statistics tables */
   for( i = 0; i < (*set)->ntables; ++i )
   {
      SCIP_CALL( SCIPtableFree(&(*set)->tables[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->tables);

   /* free display columns */
   for( i = 0; i < (*set)->ndisps; ++i )
   {
      SCIP_CALL( SCIPdispFree(&(*set)->disps[i], *set) );
   }
   BMSfreeMemoryArrayNull(&(*set)->disps);

   /* free dialogs */
   BMSfreeMemoryArrayNull(&(*set)->dialogs);

   /* free NLPIs */
   for( i = 0; i < (*set)->nnlpis; ++i )
   {
      SCIP_CALL( SCIPnlpiFree(&(*set)->nlpis[i]) );
   }
   BMSfreeMemoryArrayNull(&(*set)->nlpis);

   /* free concsolvers */
   SCIP_CALL( SCIPsetFreeConcsolvers(*set) );

   /* free concsolvers types */
   for( i = 0; i < (*set)->nconcsolvertypes; ++i )
   {
      SCIPconcsolverTypeFree(&(*set)->concsolvertypes[i]);
   }
   BMSfreeMemoryArrayNull(&(*set)->concsolvertypes);


   /* free information on external codes */
   for( i = 0; i < (*set)->nextcodes; ++i )
   {
      BMSfreeMemoryArrayNull(&(*set)->extcodenames[i]);
      BMSfreeMemoryArrayNull(&(*set)->extcodedescs[i]);
   }
   BMSfreeMemoryArrayNull(&(*set)->extcodenames);
   BMSfreeMemoryArrayNull(&(*set)->extcodedescs);

   /* free virtual tables of bandit algorithms */
   for( i = 0; i < (*set)->nbanditvtables; ++i )
   {
      SCIPbanditvtableFree(&(*set)->banditvtables[i]);
   }
   BMSfreeMemoryArrayNull(&(*set)->banditvtables);

   BMSfreeMemory(set);

   return SCIP_OKAY;
}

/** returns current stage of SCIP */
SCIP_STAGE SCIPsetGetStage(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->stage;
}

/** creates a SCIP_Bool parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetAddBool(set->paramset, messagehdlr, blkmem, name, desc, valueptr, isadvanced,
         defaultvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates an int parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetAddInt(set->paramset, messagehdlr, blkmem, name, desc, valueptr, isadvanced,
         defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a SCIP_Longint parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetAddLongint(set->paramset, messagehdlr, blkmem, name, desc, valueptr, isadvanced,
         defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a SCIP_Real parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetAddReal(set->paramset, messagehdlr, blkmem, name, desc, valueptr, isadvanced,
         defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetAddChar(set->paramset, messagehdlr, blkmem, name, desc, valueptr, isadvanced,
         defaultvalue, allowedvalues, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPsetAddStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetAddString(set->paramset, messagehdlr, blkmem, name, desc, valueptr, isadvanced,
         defaultvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** gets the fixing status value of an existing parameter */
SCIP_Bool SCIPsetIsParamFixed(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of the parameter */
   )
{
   assert(set != NULL);

   return SCIPparamsetIsFixed(set->paramset, name);
}

/** returns the pointer to the SCIP parameter with the given name */
SCIP_PARAM* SCIPsetGetParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of the parameter */
   )
{
   assert(set != NULL);

   return SCIPparamsetGetParam(set->paramset, name);
}

/** gets the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetGetBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool*            value               /**< pointer to store the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetGetBool(set->paramset, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Int parameter */
SCIP_RETCODE SCIPsetGetIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   int*                  value               /**< pointer to store the value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetGetInt(set->paramset, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPsetGetLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint*         value               /**< pointer to store the value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetGetLongint(set->paramset, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPsetGetRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Real*            value               /**< pointer to store the value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetGetReal(set->paramset, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Char parameter */
SCIP_RETCODE SCIPsetGetCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   char*                 value               /**< pointer to store the value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetGetChar(set->paramset, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing String parameter */
SCIP_RETCODE SCIPsetGetStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   char**                value               /**< pointer to store the value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetGetString(set->paramset, name, value) );

   return SCIP_OKAY;
}

/** changes the fixing status of an existing parameter */
SCIP_RETCODE SCIPsetChgParamFixed(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             fixed               /**< new fixing status of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetFix(set->paramset, name, fixed) );

   return SCIP_OKAY;
}

/** changes the value of an existing parameter */
SCIP_RETCODE SCIPsetSetParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   void*                 value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSet(set->paramset, set, messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetChgBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(set != NULL);

   retcode = SCIPparamSetBool(param, set, messagehdlr, value, FALSE, TRUE);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetSetBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetBool(set->paramset, set, messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** sets the default value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetSetDefaultBoolParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             defaultvalue        /**< new default value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetDefaultBool(set->paramset, name, defaultvalue) );

   return SCIP_OKAY;
}


/** changes the value of an existing Int parameter */
SCIP_RETCODE SCIPsetChgIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   int                   value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(set != NULL);
   assert(param != NULL);

   retcode = SCIPparamSetInt(param, set, messagehdlr, value, FALSE, TRUE);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing Int parameter */
SCIP_RETCODE SCIPsetSetIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   int                   value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetInt(set->paramset, set, messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** changes the default value of an existing Int parameter */
SCIP_RETCODE SCIPsetSetDefaultIntParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   int                   defaultvalue        /**< new default value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetDefaultInt(set->paramset, name, defaultvalue) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPsetChgLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(set != NULL);
   assert(param != NULL);

   retcode = SCIPparamSetLongint(param, set, messagehdlr, value, FALSE, TRUE);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPsetSetLongintParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetLongint(set->paramset, set, messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPsetChgRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   SCIP_Real             value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(set != NULL);
   assert(param != NULL);

   retcode = SCIPparamSetReal(param, set, messagehdlr, value, FALSE, TRUE);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPsetSetRealParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   SCIP_Real             value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetReal(set->paramset, set, messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Char parameter */
SCIP_RETCODE SCIPsetChgCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   char                  value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(set != NULL);
   assert(param != NULL);

   retcode = SCIPparamSetChar(param, set, messagehdlr, value, FALSE, TRUE);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing Char parameter */
SCIP_RETCODE SCIPsetSetCharParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   char                  value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetChar(set->paramset, set, messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing String parameter */
SCIP_RETCODE SCIPsetChgStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAM*           param,              /**< parameter */
   const char*           value               /**< new value of the parameter */
   )
{
   SCIP_RETCODE retcode;

   assert(set != NULL);
   assert(param != NULL);

   retcode = SCIPparamSetString(param, set, messagehdlr, value, TRUE);

   if( retcode != SCIP_PARAMETERWRONGVAL )
   {
      SCIP_CALL( retcode );
   }

   return retcode;
}

/** changes the value of an existing String parameter */
SCIP_RETCODE SCIPsetSetStringParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name,               /**< name of the parameter */
   const char*           value               /**< new value of the parameter */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetSetString(set->paramset, set, messagehdlr, name, value) );

   return SCIP_OKAY;
}

/** reads parameters from a file */
SCIP_RETCODE SCIPsetReadParams(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename            /**< file name */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetRead(set->paramset, set, messagehdlr, filename) );

   return SCIP_OKAY;
}

/** writes all parameters in the parameter set to a file */
SCIP_RETCODE SCIPsetWriteParams(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   )
{
   assert(set != NULL);

   SCIP_CALL( SCIPparamsetWrite(set->paramset, messagehdlr, filename, comments, onlychanged) );

   return SCIP_OKAY;
}

/** resets a single parameters to its default value */
SCIP_RETCODE SCIPsetResetParam(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   const char*           name                /**< name of the parameter */
   )
{
   SCIP_CALL( SCIPparamsetSetToDefault(set->paramset, set, messagehdlr, name) );

   return SCIP_OKAY;
}

/** resets all parameters to their default values */
SCIP_RETCODE SCIPsetResetParams(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_CALL( SCIPparamsetSetToDefaults(set->paramset, set, messagehdlr) );

   return SCIP_OKAY;
}

/** sets parameters to
 *
 *  - \ref SCIP_PARAMEMPHASIS_DEFAULT to use default values (see also SCIPsetResetParams())
 *  - \ref SCIP_PARAMEMPHASIS_COUNTER to get feasible and "fast" counting process
 *  - \ref SCIP_PARAMEMPHASIS_CPSOLVER to get CP like search (e.g. no LP relaxation)
 *  - \ref SCIP_PARAMEMPHASIS_EASYCIP to solve easy problems fast
 *  - \ref SCIP_PARAMEMPHASIS_FEASIBILITY to detect feasibility fast
 *  - \ref SCIP_PARAMEMPHASIS_HARDLP to be capable to handle hard LPs
 *  - \ref SCIP_PARAMEMPHASIS_OPTIMALITY to prove optimality fast
 *  - \ref SCIP_PARAMEMPHASIS_PHASEFEAS to find feasible solutions during a 3 phase solution process
 *  - \ref SCIP_PARAMEMPHASIS_PHASEIMPROVE to find improved solutions during a 3 phase solution process
 *  - \ref SCIP_PARAMEMPHASIS_PHASEPROOF to proof optimality during a 3 phase solution process
 */
SCIP_RETCODE SCIPsetSetEmphasis(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMEMPHASIS    paramemphasis,      /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( SCIPparamsetSetEmphasis(set->paramset, set, messagehdlr, paramemphasis, quiet) );

   return SCIP_OKAY;
}

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 */
SCIP_RETCODE SCIPsetSetSubscipsOff(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( SCIPparamsetSetToSubscipsOff(set->paramset, set, messagehdlr, quiet) );

   return SCIP_OKAY;
}

/** sets heuristic parameters values to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spent on heuristics is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristics are called more aggressively
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 */
SCIP_RETCODE SCIPsetSetHeuristics(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( SCIPparamsetSetHeuristics(set->paramset, set, messagehdlr, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** sets presolving parameters to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spent on presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving is more aggressive
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 */
SCIP_RETCODE SCIPsetSetPresolving(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( SCIPparamsetSetPresolving(set->paramset, set, messagehdlr, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** sets separating parameters to
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spent on separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that separating is more aggressive
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 */
SCIP_RETCODE SCIPsetSetSeparating(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( SCIPparamsetSetSeparating(set->paramset, set, messagehdlr, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** returns the array of all available SCIP parameters */
SCIP_PARAM** SCIPsetGetParams(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return SCIPparamsetGetParams(set->paramset);
}

/** returns the total number of all available SCIP parameters */
int SCIPsetGetNParams(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return SCIPparamsetGetNParams(set->paramset);
}

/** inserts file reader in file reader list */
SCIP_RETCODE SCIPsetIncludeReader(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_READER*          reader              /**< file reader */
   )
{
   assert(set != NULL);
   assert(reader != NULL);

   if( set->nreaders >= set->readerssize )
   {
      set->readerssize = SCIPsetCalcMemGrowSize(set, set->nreaders+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->readers, set->readerssize) );
   }
   assert(set->nreaders < set->readerssize);

   set->readers[set->nreaders] = reader;
   set->nreaders++;

   return SCIP_OKAY;
}

/** returns the file reader of the given name, or NULL if not existing */
SCIP_READER* SCIPsetFindReader(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of file reader */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nreaders; ++i )
   {
      if( strcmp(SCIPreaderGetName(set->readers[i]), name) == 0 )
         return set->readers[i];
   }

   return NULL;
}

/** inserts variable pricer in variable pricer list */
SCIP_RETCODE SCIPsetIncludePricer(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   assert(set != NULL);
   assert(pricer != NULL);

   if( set->npricers >= set->pricerssize )
   {
      set->pricerssize = SCIPsetCalcMemGrowSize(set, set->npricers+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->pricers, set->pricerssize) );
   }
   assert(set->npricers < set->pricerssize);

   set->pricers[set->npricers] = pricer;
   set->npricers++;
   set->pricerssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the variable pricer of the given name, or NULL if not existing */
SCIP_PRICER* SCIPsetFindPricer(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of variable pricer */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->npricers; ++i )
   {
      if( strcmp(SCIPpricerGetName(set->pricers[i]), name) == 0 )
         return set->pricers[i];
   }

   return NULL;
}

/** sorts pricers by priorities */
void SCIPsetSortPricers(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->pricerssorted )
   {
      SCIPsortPtr((void**)set->pricers, SCIPpricerComp, set->npricers);
      set->pricerssorted = TRUE;
      set->pricersnamesorted = FALSE;
   }
}

/** sorts pricers by name */
void SCIPsetSortPricersName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->pricersnamesorted )
   {
      SCIPsortPtr((void**)set->pricers, SCIPpricerCompName, set->npricers);
      set->pricerssorted = FALSE;
      set->pricersnamesorted = TRUE;
   }
}

/** inserts constraint handler in constraint handler list */
SCIP_RETCODE SCIPsetIncludeConshdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   int priority;
   int i;

   assert(set != NULL);
   assert(conshdlr != NULL);
   assert(!SCIPconshdlrIsInitialized(conshdlr));

   /* allocate memory */
   if( set->nconshdlrs >= set->conshdlrssize )
   {
      set->conshdlrssize = SCIPsetCalcMemGrowSize(set, set->nconshdlrs+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->conshdlrs, set->conshdlrssize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&set->conshdlrs_sepa, set->conshdlrssize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&set->conshdlrs_enfo, set->conshdlrssize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&set->conshdlrs_include, set->conshdlrssize) );
   }
   assert(set->nconshdlrs < set->conshdlrssize);

   /* sort constraint handler into conshdlrs array sorted by check priority */
   priority = SCIPconshdlrGetCheckPriority(conshdlr);
   for( i = set->nconshdlrs; i > 0 && SCIPconshdlrGetCheckPriority(set->conshdlrs[i-1]) < priority; --i )
   {
      set->conshdlrs[i] = set->conshdlrs[i-1];
   }
   set->conshdlrs[i] = conshdlr;

   /* sort constraint handler into conshdlrs_sepa array sorted by sepa priority */
   priority = SCIPconshdlrGetSepaPriority(conshdlr);
   for( i = set->nconshdlrs; i > 0 && SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i-1]) < priority; --i )
   {
      set->conshdlrs_sepa[i] = set->conshdlrs_sepa[i-1];
   }
   set->conshdlrs_sepa[i] = conshdlr;

   /* sort constraint handler into conshdlrs_enfo array sorted by enfo priority */
   priority = SCIPconshdlrGetEnfoPriority(conshdlr);
   for( i = set->nconshdlrs; i > 0 && SCIPconshdlrGetEnfoPriority(set->conshdlrs_enfo[i-1]) < priority; --i )
   {
      set->conshdlrs_enfo[i] = set->conshdlrs_enfo[i-1];
   }
   set->conshdlrs_enfo[i] = conshdlr;

   /* add constraint handler into conshdlrs_include array sorted by inclusion order */
   set->conshdlrs_include[set->nconshdlrs] = conshdlr;

   set->nconshdlrs++;

   return SCIP_OKAY;
}

/** reinserts a constraint handler with modified sepa priority into the sepa priority sorted array */
void SCIPsetReinsertConshdlrSepaPrio(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler to be reinserted */
   int                   oldpriority         /**< the old separation priority of constraint handler */
   )
{
   int newpriority;
   int newpos;
   int i;
   assert(set != NULL);
   assert(conshdlr != NULL);

   newpriority = SCIPconshdlrGetSepaPriority(conshdlr);
   newpos = -1;

   /* search for the old position of constraint handler; determine its new position at the same time */
   if( newpriority > oldpriority )
   {
      i = 0;
      while( i < set->nconshdlrs &&
            strcmp(SCIPconshdlrGetName(set->conshdlrs_sepa[i]), SCIPconshdlrGetName(conshdlr)) != 0 )
      {
         int priorityatpos;

         priorityatpos = SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i]);
         assert(priorityatpos >= oldpriority);

         /* current index is the position to insert the constraint handler */
         if( newpriority > priorityatpos && newpos == -1 )
            newpos = i;

         ++i;
      }
      assert(i < set->nconshdlrs);

      /* constraint must change its position in array */
      if( newpos != -1 )
      {
         /* shift all constraint handlers between old and new position by one, and insert constraint handler */
         for( ; i > newpos; --i )
         {
            set->conshdlrs_sepa[i] = set->conshdlrs_sepa[i-1];
         }
         set->conshdlrs_sepa[newpos] = conshdlr;
      }

   }
   else if( newpriority < oldpriority )
   {
      i = set->nconshdlrs - 1;
      while( i >= 0 &&
                  strcmp(SCIPconshdlrGetName(set->conshdlrs_sepa[i]), SCIPconshdlrGetName(conshdlr)) != 0 )
      {
         int priorityatpos;

         priorityatpos = SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i]);
         assert(priorityatpos <= oldpriority);

         /* current index is the position to insert the constraint handler */
         if( newpriority < priorityatpos && newpos == -1 )
            newpos = i;

         --i;
      }
      assert(i >= 0);

      /* constraint must change its position in array */
      if( newpos != -1 )
      {
         /* shift all constraint handlers between old and new position by one, and insert constraint handler */
         for(; i < newpos; ++i )
         {
            set->conshdlrs_sepa[i] = set->conshdlrs_sepa[i + 1];
         }
         set->conshdlrs_sepa[newpos] = conshdlr;
      }
#ifndef NDEBUG
      for( i = 0; i < set->nconshdlrs - 1; ++i )
         assert(SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i])
               >= SCIPconshdlrGetSepaPriority(set->conshdlrs_sepa[i + 1]));
#endif
   }
}

/** returns the constraint handler of the given name, or NULL if not existing */
SCIP_CONSHDLR* SCIPsetFindConshdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of constraint handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nconshdlrs; ++i )
   {
      if( strcmp(SCIPconshdlrGetName(set->conshdlrs[i]), name) == 0 )
         return set->conshdlrs[i];
   }

   return NULL;
}

/** inserts conflict handler in conflict handler list */
SCIP_RETCODE SCIPsetIncludeConflicthdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(set != NULL);
   assert(conflicthdlr != NULL);
   assert(!SCIPconflicthdlrIsInitialized(conflicthdlr));

   if( set->nconflicthdlrs >= set->conflicthdlrssize )
   {
      set->conflicthdlrssize = SCIPsetCalcMemGrowSize(set, set->nconflicthdlrs+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->conflicthdlrs, set->conflicthdlrssize) );
   }
   assert(set->nconflicthdlrs < set->conflicthdlrssize);

   set->conflicthdlrs[set->nconflicthdlrs] = conflicthdlr;
   set->nconflicthdlrs++;
   set->conflicthdlrssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the conflict handler of the given name, or NULL if not existing */
SCIP_CONFLICTHDLR* SCIPsetFindConflicthdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of conflict handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      if( strcmp(SCIPconflicthdlrGetName(set->conflicthdlrs[i]), name) == 0 )
         return set->conflicthdlrs[i];
   }

   return NULL;
}

/** sorts conflict handlers by priorities */
void SCIPsetSortConflicthdlrs(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->conflicthdlrssorted )
   {
      SCIPsortPtr((void**)set->conflicthdlrs, SCIPconflicthdlrComp, set->nconflicthdlrs);
      set->conflicthdlrssorted = TRUE;
      set->conflicthdlrsnamesorted = FALSE;
   }
}

/** sorts conflict handlers by name */
void SCIPsetSortConflicthdlrsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->conflicthdlrsnamesorted )
   {
      SCIPsortPtr((void**)set->conflicthdlrs, SCIPconflicthdlrCompName, set->nconflicthdlrs);
      set->conflicthdlrssorted = FALSE;
      set->conflicthdlrsnamesorted = TRUE;
   }
}

/** inserts presolver in presolver list */
SCIP_RETCODE SCIPsetIncludePresol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRESOL*          presol              /**< presolver */
   )
{
   assert(set != NULL);
   assert(presol != NULL);

   if( set->npresols >= set->presolssize )
   {
      set->presolssize = SCIPsetCalcMemGrowSize(set, set->npresols+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->presols, set->presolssize) );
   }
   assert(set->npresols < set->presolssize);

   set->presols[set->npresols] = presol;
   set->npresols++;
   set->presolssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the presolver of the given name, or NULL if not existing */
SCIP_PRESOL* SCIPsetFindPresol(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of presolver */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->npresols; ++i )
   {
      if( strcmp(SCIPpresolGetName(set->presols[i]), name) == 0 )
         return set->presols[i];
   }

   return NULL;
}

/** sorts presolvers by priorities */
void SCIPsetSortPresols(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->presolssorted )
   {
      SCIPsortPtr((void**)set->presols, SCIPpresolComp, set->npresols);
      set->presolssorted = TRUE;
      set->presolsnamesorted = FALSE;
   }
}

/** sorts presolvers by name */
void SCIPsetSortPresolsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->presolsnamesorted )
   {
      SCIPsortPtr((void**)set->presols, SCIPpresolCompName, set->npresols);
      set->presolssorted = FALSE;
      set->presolsnamesorted = TRUE;
   }
}

/** inserts relaxator in relaxator list */
SCIP_RETCODE SCIPsetIncludeRelax(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(set != NULL);
   assert(relax != NULL);
   assert(!SCIPrelaxIsInitialized(relax));

   if( set->nrelaxs >= set->relaxssize )
   {
      set->relaxssize = SCIPsetCalcMemGrowSize(set, set->nrelaxs+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->relaxs, set->relaxssize) );
   }
   assert(set->nrelaxs < set->relaxssize);

   set->relaxs[set->nrelaxs] = relax;
   set->nrelaxs++;
   set->relaxssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the relaxator of the given name, or NULL if not existing */
SCIP_RELAX* SCIPsetFindRelax(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of relaxator */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nrelaxs; ++i )
   {
      if( strcmp(SCIPrelaxGetName(set->relaxs[i]), name) == 0 )
         return set->relaxs[i];
   }

   return NULL;
}

/** sorts relaxators by priorities */
void SCIPsetSortRelaxs(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->relaxssorted )
   {
      SCIPsortPtr((void**)set->relaxs, SCIPrelaxComp, set->nrelaxs);
      set->relaxssorted = TRUE;
      set->relaxsnamesorted = FALSE;
   }
}

/** sorts relaxators by priorities */
void SCIPsetSortRelaxsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->relaxsnamesorted )
   {
      SCIPsortPtr((void**)set->relaxs, SCIPrelaxCompName, set->nrelaxs);
      set->relaxssorted = FALSE;
      set->relaxsnamesorted = TRUE;
   }
}

/** inserts separator in separator list */
SCIP_RETCODE SCIPsetIncludeSepa(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(set != NULL);
   assert(sepa != NULL);
   assert(!SCIPsepaIsInitialized(sepa));

   if( set->nsepas >= set->sepassize )
   {
      set->sepassize = SCIPsetCalcMemGrowSize(set, set->nsepas+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->sepas, set->sepassize) );
   }
   assert(set->nsepas < set->sepassize);

   set->sepas[set->nsepas] = sepa;
   set->nsepas++;
   set->sepassorted = FALSE;

   return SCIP_OKAY;
}

/** returns the separator of the given name, or NULL if not existing */
SCIP_SEPA* SCIPsetFindSepa(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of separator */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nsepas; ++i )
   {
      if( strcmp(SCIPsepaGetName(set->sepas[i]), name) == 0 )
         return set->sepas[i];
   }

   return NULL;
}

/** sorts separators by priorities */
void SCIPsetSortSepas(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->sepassorted )
   {
      SCIPsortPtr((void**)set->sepas, SCIPsepaComp, set->nsepas);
      set->sepassorted = TRUE;
      set->sepasnamesorted = FALSE;
   }
}

/** sorts separators by name */
void SCIPsetSortSepasName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->sepasnamesorted )
   {
      SCIPsortPtr((void**)set->sepas, SCIPsepaCompName, set->nsepas);
      set->sepassorted = FALSE;
      set->sepasnamesorted = TRUE;
   }
}

/** inserts propagator in propagator list */
SCIP_RETCODE SCIPsetIncludeProp(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(set != NULL);
   assert(prop != NULL);
   assert(!SCIPpropIsInitialized(prop));

   if( set->nprops >= set->propssize )
   {
      set->propssize = SCIPsetCalcMemGrowSize(set, set->nprops+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->props, set->propssize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&set->props_presol, set->propssize) );
   }
   assert(set->nprops < set->propssize);

   set->props[set->nprops] = prop;
   set->props_presol[set->nprops] = prop;
   set->nprops++;
   set->propssorted = FALSE;
   set->propspresolsorted = FALSE;

   return SCIP_OKAY;
}

/** returns the propagator of the given name, or NULL if not existing */
SCIP_PROP* SCIPsetFindProp(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of propagator */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nprops; ++i )
   {
      if( strcmp(SCIPpropGetName(set->props[i]), name) == 0 )
         return set->props[i];
   }

   return NULL;
}

/** sorts propagators by priorities */
void SCIPsetSortProps(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->propssorted )
   {
      SCIPsortPtr((void**)set->props, SCIPpropComp, set->nprops);
      set->propssorted = TRUE;
      set->propsnamesorted = FALSE;
   }
}

/** sorts propagators by priorities for presolving */
void SCIPsetSortPropsPresol(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->propspresolsorted )
   {
      SCIPsortPtr((void**)set->props_presol, SCIPpropCompPresol, set->nprops);
      set->propspresolsorted = TRUE;
      set->propsnamesorted = FALSE;
   }
}

/** sorts propagators w.r.t. names */
void SCIPsetSortPropsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->propsnamesorted )
   {
      SCIPsortPtr((void**)set->props, SCIPpropCompName, set->nprops);
      set->propssorted = FALSE;
      set->propsnamesorted = TRUE;
   }
}

/** inserts bandit virtual function table into set */
SCIP_RETCODE SCIPsetIncludeBanditvtable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BANDITVTABLE*    banditvtable        /**< bandit algorithm virtual function table */
   )
{
   assert(set != NULL);
   assert(banditvtable != NULL);

   if( set->nbanditvtables >= set->banditvtablessize )
   {
      int newsize = SCIPsetCalcMemGrowSize(set, set->nbanditvtables + 1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->banditvtables, newsize) );
      set->banditvtablessize = newsize;
   }

   assert(set->nbanditvtables < set->banditvtablessize);
   set->banditvtables[set->nbanditvtables++] = banditvtable;

   return SCIP_OKAY;
}

/** returns the bandit virtual function table of the given name, or NULL if not existing */
SCIP_BANDITVTABLE* SCIPsetFindBanditvtable(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of bandit algorithm virtual function table */
   )
{
   int b;

   assert(set != NULL);
   assert(name != NULL);

   /* search for a bandit v table of the given name */
   for( b = 0; b < set->nbanditvtables; ++b )
   {
      if( strcmp(name, SCIPbanditvtableGetName(set->banditvtables[b])) == 0 )
         return set->banditvtables[b];
   }

   return NULL;
}

/** inserts concurrent solver type into the concurrent solver type list */
SCIP_RETCODE SCIPsetIncludeConcsolverType(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVERTYPE*  concsolvertype      /**< concurrent solver type */
   )
{
   assert(set != NULL);
   assert(concsolvertype != NULL);

   if( set->nconcsolvertypes >= set->concsolvertypessize )
   {
      set->concsolvertypessize = SCIPsetCalcMemGrowSize(set, set->nconcsolvertypes + 1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->concsolvertypes, set->concsolvertypessize) );
   }
   assert(set->nconcsolvertypes < set->concsolvertypessize);

   set->concsolvertypes[set->nconcsolvertypes] = concsolvertype;
   set->nconcsolvertypes++;

   return SCIP_OKAY;
}

/** returns the concurrent solver type with the given name, or NULL if not existing */
SCIP_CONCSOLVERTYPE* SCIPsetFindConcsolverType(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of concurrent solver type */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nconcsolvertypes; ++i )
   {
      if( strcmp(SCIPconcsolverTypeGetName(set->concsolvertypes[i]), name) == 0 )
         return set->concsolvertypes[i];
   }

   return NULL;
}

/** inserts concurrent solver into the concurrent solver list */
SCIP_RETCODE SCIPsetIncludeConcsolver(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONCSOLVER*      concsolver          /**< concurrent solver */
   )
{
   assert(set != NULL);
   assert(concsolver != NULL);

   if( set->nconcsolvers >= set->concsolverssize )
   {
      set->concsolverssize = SCIPsetCalcMemGrowSize(set, set->nconcsolvers + 1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->concsolvers, set->concsolverssize) );
   }
   assert(set->nconcsolvers < set->concsolverssize);

   set->concsolvers[set->nconcsolvers] = concsolver;
   assert(set->nconcsolvers == SCIPconcsolverGetIdx(concsolver));

   set->nconcsolvers++;

   return SCIP_OKAY;
}

/** frees all concurrent solvers in the concurrent solver list */
SCIP_RETCODE SCIPsetFreeConcsolvers(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int i;
   assert(set != NULL);

   /* call user callback for each concurrent solver */
   for( i = 0; i < set->nconcsolvers; ++i )
   {
      SCIP_CALL( SCIPconcsolverDestroyInstance(set, &set->concsolvers[i]) );
   }

   /* set size and number to zero and free the concurent solver array */
   set->nconcsolvers = 0;
   set->concsolverssize = 0;
   BMSfreeMemoryArrayNull(&set->concsolvers);

   return SCIP_OKAY;
}

/** inserts primal heuristic in primal heuristic list */
SCIP_RETCODE SCIPsetIncludeHeur(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_HEUR*            heur                /**< primal heuristic */
   )
{
   assert(set != NULL);
   assert(heur != NULL);
   assert(!SCIPheurIsInitialized(heur));

   if( set->nheurs >= set->heurssize )
   {
      set->heurssize = SCIPsetCalcMemGrowSize(set, set->nheurs+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->heurs, set->heurssize) );
   }
   assert(set->nheurs < set->heurssize);

   set->heurs[set->nheurs] = heur;
   set->nheurs++;
   set->heurssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the primal heuristic of the given name, or NULL if not existing */
SCIP_HEUR* SCIPsetFindHeur(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of primal heuristic */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nheurs; ++i )
   {
      if( strcmp(SCIPheurGetName(set->heurs[i]), name) == 0 )
         return set->heurs[i];
   }

   return NULL;
}

/** sorts heuristics by priorities */
void SCIPsetSortHeurs(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->heurssorted )
   {
      SCIPsortPtr((void**)set->heurs, SCIPheurComp, set->nheurs);
      set->heurssorted = TRUE;
      set->heursnamesorted = FALSE;
   }
}

/** sorts heuristics by names */
void SCIPsetSortHeursName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->heursnamesorted )
   {
      SCIPsortPtr((void**)set->heurs, SCIPheurCompName, set->nheurs);
      set->heurssorted = FALSE;
      set->heursnamesorted = TRUE;
   }
}

/** inserts tree compression in tree compression list */
SCIP_RETCODE SCIPsetIncludeCompr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_COMPR*           compr               /**< tree compression */
   )
{
   assert(set != NULL);
   assert(compr != NULL);
   assert(!SCIPcomprIsInitialized(compr));

   if( set->ncomprs >= set->comprssize )
   {
      set->comprssize = SCIPsetCalcMemGrowSize(set, set->ncomprs+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->comprs, set->comprssize) );
   }
   assert(set->ncomprs < set->comprssize);

   set->comprs[set->ncomprs] = compr;
   set->ncomprs++;
   set->comprssorted = FALSE;

   return SCIP_OKAY;
}

/** returns the tree compression of the given name, or NULL if not existing */
SCIP_COMPR* SCIPsetFindCompr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of tree compression */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->ncomprs; ++i )
   {
      if( strcmp(SCIPcomprGetName(set->comprs[i]), name) == 0 )
         return set->comprs[i];
   }

   return NULL;
}

/** sorts compressions by priorities */
void SCIPsetSortComprs(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->comprssorted )
   {
      SCIPsortPtr((void**)set->comprs, SCIPcomprComp, set->ncomprs);
      set->comprssorted = TRUE;
      set->comprsnamesorted = FALSE;
   }
}

/** sorts heuristics by names */
void SCIPsetSortComprsName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->comprsnamesorted )
   {
      SCIPsortPtr((void**)set->comprs, SCIPcomprCompName, set->ncomprs);
      set->comprssorted = FALSE;
      set->comprsnamesorted = TRUE;
   }
}

/** inserts event handler in event handler list */
SCIP_RETCODE SCIPsetIncludeEventhdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(set != NULL);
   assert(eventhdlr != NULL);
   assert(!SCIPeventhdlrIsInitialized(eventhdlr));

   if( set->neventhdlrs >= set->eventhdlrssize )
   {
      set->eventhdlrssize = SCIPsetCalcMemGrowSize(set, set->neventhdlrs+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->eventhdlrs, set->eventhdlrssize) );
   }
   assert(set->neventhdlrs < set->eventhdlrssize);

   set->eventhdlrs[set->neventhdlrs] = eventhdlr;
   set->neventhdlrs++;

   return SCIP_OKAY;
}

/** returns the event handler of the given name, or NULL if not existing */
SCIP_EVENTHDLR* SCIPsetFindEventhdlr(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of event handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->neventhdlrs; ++i )
   {
      if( strcmp(SCIPeventhdlrGetName(set->eventhdlrs[i]), name) == 0 )
         return set->eventhdlrs[i];
   }

   return NULL;
}

/** inserts node selector in node selector list */
SCIP_RETCODE SCIPsetIncludeNodesel(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODESEL*         nodesel             /**< node selector */
   )
{
   int i;
   int nodeselstdprio;

   assert(set != NULL);
   assert(nodesel != NULL);
   assert(!SCIPnodeselIsInitialized(nodesel));

   if( set->nnodesels >= set->nodeselssize )
   {
      set->nodeselssize = SCIPsetCalcMemGrowSize(set, set->nnodesels+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->nodesels, set->nodeselssize) );
   }
   assert(set->nnodesels < set->nodeselssize);

   nodeselstdprio = SCIPnodeselGetStdPriority(nodesel);

   for( i = set->nnodesels; i > 0 && nodeselstdprio > SCIPnodeselGetStdPriority(set->nodesels[i-1]); --i )
      set->nodesels[i] = set->nodesels[i-1];

   set->nodesels[i] = nodesel;
   set->nnodesels++;

   return SCIP_OKAY;
}

/** returns the node selector of the given name, or NULL if not existing */
SCIP_NODESEL* SCIPsetFindNodesel(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of event handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nnodesels; ++i )
   {
      if( strcmp(SCIPnodeselGetName(set->nodesels[i]), name) == 0 )
         return set->nodesels[i];
   }

   return NULL;
}

/** returns node selector with highest priority in the current mode */
SCIP_NODESEL* SCIPsetGetNodesel(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(set != NULL);
   assert(stat != NULL);

   /* check, if old node selector is still valid */
   if( set->nodesel == NULL && set->nnodesels > 0 )
   {
      int i;

      set->nodesel = set->nodesels[0];

      /* search highest priority node selector */
      if( stat->memsavemode )
      {
         for( i = 1; i < set->nnodesels; ++i )
         {
            if( SCIPnodeselGetMemsavePriority(set->nodesels[i]) > SCIPnodeselGetMemsavePriority(set->nodesel) )
               set->nodesel = set->nodesels[i];
         }
      }
      else
      {
         for( i = 1; i < set->nnodesels; ++i )
         {
            if( SCIPnodeselGetStdPriority(set->nodesels[i]) > SCIPnodeselGetStdPriority(set->nodesel) )
               set->nodesel = set->nodesels[i];
         }
      }
   }

   return set->nodesel;
}

/** inserts branching rule in branching rule list */
SCIP_RETCODE SCIPsetIncludeBranchrule(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(set != NULL);
   assert(branchrule != NULL);
   assert(!SCIPbranchruleIsInitialized(branchrule));

   if( set->nbranchrules >= set->branchrulessize )
   {
      set->branchrulessize = SCIPsetCalcMemGrowSize(set, set->nbranchrules+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->branchrules, set->branchrulessize) );
   }
   assert(set->nbranchrules < set->branchrulessize);

   set->branchrules[set->nbranchrules] = branchrule;
   set->nbranchrules++;
   set->branchrulessorted = FALSE;

   return SCIP_OKAY;
}

/** returns the branching rule of the given name, or NULL if not existing */
SCIP_BRANCHRULE* SCIPsetFindBranchrule(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of event handler */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nbranchrules; ++i )
   {
      if( strcmp(SCIPbranchruleGetName(set->branchrules[i]), name) == 0 )
         return set->branchrules[i];
   }

   return NULL;
}

/** sorts branching rules by priorities */
void SCIPsetSortBranchrules(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->branchrulessorted )
   {
      SCIPsortPtr((void**)set->branchrules, SCIPbranchruleComp, set->nbranchrules);
      set->branchrulessorted = TRUE;
      set->branchrulesnamesorted = FALSE;
   }
}

/** sorts branching rules by priorities */
void SCIPsetSortBranchrulesName(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->branchrulesnamesorted )
   {
      SCIPsortPtr((void**)set->branchrules, SCIPbranchruleCompName, set->nbranchrules);
      set->branchrulessorted = FALSE;
      set->branchrulesnamesorted = TRUE;
   }
}

/** inserts display column in display column list */
SCIP_RETCODE SCIPsetIncludeDisp(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DISP*            disp                /**< display column */
   )
{
   int i;
   int disppos;

   assert(set != NULL);
   assert(disp != NULL);
   assert(!SCIPdispIsInitialized(disp));

   if( set->ndisps >= set->dispssize )
   {
      set->dispssize = SCIPsetCalcMemGrowSize(set, set->ndisps+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->disps, set->dispssize) );
   }
   assert(set->ndisps < set->dispssize);

   disppos = SCIPdispGetPosition(disp);

   for( i = set->ndisps; i > 0 && disppos < SCIPdispGetPosition(set->disps[i-1]); --i )
   {
      set->disps[i] = set->disps[i-1];
   }
   set->disps[i] = disp;
   set->ndisps++;

   return SCIP_OKAY;
}

/** returns the display column of the given name, or NULL if not existing */
SCIP_DISP* SCIPsetFindDisp(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of display */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->ndisps; ++i )
   {
      if( strcmp(SCIPdispGetName(set->disps[i]), name) == 0 )
         return set->disps[i];
   }

   return NULL;
}

/** inserts statistics table in statistics table list */
SCIP_RETCODE SCIPsetIncludeTable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TABLE*           table               /**< statistics table */
   )
{
   assert(set != NULL);
   assert(table != NULL);
   assert(!SCIPtableIsInitialized(table));

   if( set->ntables >= set->tablessize )
   {
      set->tablessize = SCIPsetCalcMemGrowSize(set, set->ntables+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->tables, set->tablessize) );
   }
   assert(set->ntables < set->tablessize);

   /* we insert in arbitrary order and sort once before printing statistics */
   set->tables[set->ntables] = table;
   set->ntables++;
   set->tablessorted = FALSE;

   return SCIP_OKAY;
}

/** returns the statistics table of the given name, or NULL if not existing */
SCIP_TABLE* SCIPsetFindTable(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of statistics table */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->ntables; ++i )
   {
      if( strcmp(SCIPtableGetName(set->tables[i]), name) == 0 )
         return set->tables[i];
   }

   return NULL;
}

/** inserts dialog in dialog list */
SCIP_RETCODE SCIPsetIncludeDialog(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(set != NULL);
   assert(dialog != NULL);

   if( set->ndialogs >= set->dialogssize )
   {
      set->dialogssize = SCIPsetCalcMemGrowSize(set, set->ndialogs+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->dialogs, set->dialogssize) );
   }
   assert(set->ndialogs < set->dialogssize);

   set->dialogs[set->ndialogs] = dialog;
   set->ndialogs++;

   return SCIP_OKAY;
}

/** returns if the dialog already exists */
SCIP_Bool SCIPsetExistsDialog(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   int i;

   assert(set != NULL);

   if( dialog == NULL )
      return FALSE;

   for( i = 0; i < set->ndialogs; ++i )
   {
      if( set->dialogs[i] == dialog )
         return TRUE;
   }

   return FALSE;
}

/** inserts NLPI in NLPI list */
SCIP_RETCODE SCIPsetIncludeNlpi(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi                /**< NLPI */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);

   if( set->nnlpis >= set->nlpissize )
   {
      set->nlpissize = SCIPsetCalcMemGrowSize(set, set->nnlpis+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->nlpis, set->nlpissize) );
   }
   assert(set->nnlpis < set->nlpissize);

   set->nlpis[set->nnlpis] = nlpi;
   set->nnlpis++;
   set->nlpissorted = FALSE;

   return SCIP_OKAY;
}

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_NLPI* SCIPsetFindNlpi(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name                /**< name of NLPI */
   )
{
   int i;

   assert(set != NULL);
   assert(name != NULL);

   for( i = 0; i < set->nnlpis; ++i )
   {
      if( strcmp(SCIPnlpiGetName(set->nlpis[i]), name) == 0 )
         return set->nlpis[i];
   }

   return NULL;
}

/** sorts NLPIs by priorities */
void SCIPsetSortNlpis(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( !set->nlpissorted )
   {
      SCIPsortPtr((void**)set->nlpis, SCIPnlpiComp, set->nnlpis);
      set->nlpissorted = TRUE;
   }
}

/** set priority of an NLPI */
void SCIPsetSetPriorityNlpi(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of NLPI */
   )
{
   assert(set != NULL);
   assert(nlpi != NULL);

   SCIPnlpiSetPriority(nlpi, priority);
   set->nlpissorted = FALSE;
}

/** inserts information about an external code in external codes list */
SCIP_RETCODE SCIPsetIncludeExternalCode(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of external code */
   const char*           description         /**< description of external code, can be NULL */
   )
{
   assert(set  != NULL);
   assert(name != NULL);

   if( set->nextcodes >= set->extcodessize )
   {
      set->extcodessize = SCIPsetCalcMemGrowSize(set, set->nextcodes+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&set->extcodenames, set->extcodessize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&set->extcodedescs, set->extcodessize) );
   }
   assert(set->nextcodes < set->extcodessize);

   BMSduplicateMemoryArray(&(set->extcodenames[set->nextcodes]), name, (int) (strlen(name)+1));  /*lint !e866*/
   if( description != NULL )
   {
      BMSduplicateMemoryArray(&(set->extcodedescs[set->nextcodes]), description, (int) (strlen(description)+1));  /*lint !e866*/
   }
   else
   {
      set->extcodedescs[set->nextcodes] = NULL;
   }
   set->nextcodes++;

   return SCIP_OKAY;
}

/** calls init methods of all plugins */
SCIP_RETCODE SCIPsetInitPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   int i;

   assert(set != NULL);

   /* active variable pricers */
   SCIPsetSortPricers(set);
   for( i = 0; i < set->nactivepricers; ++i )
   {
      SCIP_CALL( SCIPpricerInit(set->pricers[i], set) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      SCIP_CALL( SCIPconshdlrInit(set->conshdlrs[i], blkmem, set, stat) );
   }

   /* conflict handlers */
   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      SCIP_CALL( SCIPconflicthdlrInit(set->conflicthdlrs[i], set) );
   }

   /* presolvers */
   for( i = 0; i < set->npresols; ++i )
   {
      SCIP_CALL( SCIPpresolInit(set->presols[i], set) );
   }

   /* relaxators */
   for( i = 0; i < set->nrelaxs; ++i )
   {
      SCIP_CALL( SCIPrelaxInit(set->relaxs[i], set) );
   }

   /* separators */
   for( i = 0; i < set->nsepas; ++i )
   {
      SCIP_CALL( SCIPsepaInit(set->sepas[i], set) );
   }

   /* propagators */
   for( i = 0; i < set->nprops; ++i )
   {
      SCIP_CALL( SCIPpropInit(set->props[i], set) );
   }

   /* primal heuristics */
   for( i = 0; i < set->nheurs; ++i )
   {
      SCIP_CALL( SCIPheurInit(set->heurs[i], set) );
   }

   /* tree compression */
   for( i = 0; i < set->ncomprs; ++i )
   {
      SCIP_CALL( SCIPcomprInit(set->comprs[i], set) );
   }

   /* event handlers */
   for( i = 0; i < set->neventhdlrs; ++i )
   {
      SCIP_CALL( SCIPeventhdlrInit(set->eventhdlrs[i], set) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      SCIP_CALL( SCIPnodeselInit(set->nodesels[i], set) );
   }

   /* branching rules */
   for( i = 0; i < set->nbranchrules; ++i )
   {
      SCIP_CALL( SCIPbranchruleInit(set->branchrules[i], set) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      SCIP_CALL( SCIPdispInit(set->disps[i], set) );
   }
   SCIP_CALL( SCIPdispAutoActivate(set) );

   /* statistics tables */
   for( i = 0; i < set->ntables; ++i )
   {
      SCIP_CALL( SCIPtableInit(set->tables[i], set) );
   }

   return SCIP_OKAY;
}

/** calls exit methods of all plugins */
SCIP_RETCODE SCIPsetExitPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   int i;

   assert(set != NULL);

   /* active variable pricers */
   SCIPsetSortPricers(set);
   for( i = 0; i < set->nactivepricers; ++i )
   {
      SCIP_CALL( SCIPpricerExit(set->pricers[i], set) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      SCIP_CALL( SCIPconshdlrExit(set->conshdlrs[i], blkmem, set, stat) );
   }

   /* conflict handlers */
   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      SCIP_CALL( SCIPconflicthdlrExit(set->conflicthdlrs[i], set) );
   }

   /* presolvers */
   for( i = 0; i < set->npresols; ++i )
   {
      SCIP_CALL( SCIPpresolExit(set->presols[i], set) );
   }

   /* relaxators */
   for( i = 0; i < set->nrelaxs; ++i )
   {
      SCIP_CALL( SCIPrelaxExit(set->relaxs[i], set) );
   }

   /* separators */
   for( i = 0; i < set->nsepas; ++i )
   {
      SCIP_CALL( SCIPsepaExit(set->sepas[i], set) );
   }

   /* propagators */
   for( i = 0; i < set->nprops; ++i )
   {
      SCIP_CALL( SCIPpropExit(set->props[i], set) );
   }

   /* primal heuristics */
   for( i = 0; i < set->nheurs; ++i )
   {
      SCIP_CALL( SCIPheurExit(set->heurs[i], set) );
   }

   /* tree compression */
   for( i = 0; i < set->ncomprs; ++i )
   {
      SCIP_CALL( SCIPcomprExit(set->comprs[i], set) );
   }

   /* event handlers */
   for( i = 0; i < set->neventhdlrs; ++i )
   {
      SCIP_CALL( SCIPeventhdlrExit(set->eventhdlrs[i], set) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      SCIP_CALL( SCIPnodeselExit(set->nodesels[i], set) );
   }

   /* branching rules */
   for( i = 0; i < set->nbranchrules; ++i )
   {
      SCIP_CALL( SCIPbranchruleExit(set->branchrules[i], set) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      SCIP_CALL( SCIPdispExit(set->disps[i], set) );
   }

   /* statistics tables */
   for( i = 0; i < set->ntables; ++i )
   {
      SCIP_CALL( SCIPtableExit(set->tables[i], set) );
   }

   return SCIP_OKAY;
}

/** calls initpre methods of all plugins */
SCIP_RETCODE SCIPsetInitprePlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   int i;

   assert(set != NULL);

   /* inform presolvers that the presolving is abound to begin */
   for( i = 0; i < set->npresols; ++i )
   {
      SCIP_CALL( SCIPpresolInitpre(set->presols[i], set) );
   }

   /* inform propagators that the presolving is abound to begin */
   for( i = 0; i < set->nprops; ++i )
   {
      SCIP_CALL( SCIPpropInitpre(set->props[i], set) );
   }

   /* inform constraint handlers that the presolving is abound to begin */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      SCIP_CALL( SCIPconshdlrInitpre(set->conshdlrs[i], blkmem, set, stat) );
   }

   return SCIP_OKAY;
}

/** calls exitpre methods of all plugins */
SCIP_RETCODE SCIPsetExitprePlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   int i;

   assert(set != NULL);

   /* inform presolvers that the presolving is abound to begin */
   for( i = 0; i < set->npresols; ++i )
   {
      SCIP_CALL( SCIPpresolExitpre(set->presols[i], set) );
   }

   /* inform propagators that the presolving is abound to begin */
   for( i = 0; i < set->nprops; ++i )
   {
      SCIP_CALL( SCIPpropExitpre(set->props[i], set) );
   }

   /* inform constraint handlers that the presolving is abound to begin */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      SCIP_CALL( SCIPconshdlrExitpre(set->conshdlrs[i], blkmem, set, stat) );
   }

   return SCIP_OKAY;
}

/** calls initsol methods of all plugins */
SCIP_RETCODE SCIPsetInitsolPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   int i;

   assert(set != NULL);

   /* reset SCIP-defined feasibility tolerance for relaxations
    * if this is invalid, then only the relaxation specific feasibility tolerance,
    * e.g., numerics/lpfeastol is applied
    * SCIP plugins or core may set num_relaxfeastol to request a
    * tighter feasibility tolerance, though
    * see also documentation of SCIPchgRelaxfeastol
    */
   set->num_relaxfeastol = SCIP_INVALID;

   /* active variable pricers */
   SCIPsetSortPricers(set);
   for( i = 0; i < set->nactivepricers; ++i )
   {
      SCIP_CALL( SCIPpricerInitsol(set->pricers[i], set) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      SCIP_CALL( SCIPconshdlrInitsol(set->conshdlrs[i], blkmem, set, stat) );
   }

   /* conflict handlers */
   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      SCIP_CALL( SCIPconflicthdlrInitsol(set->conflicthdlrs[i], set) );
   }

   /* relaxators */
   for( i = 0; i < set->nrelaxs; ++i )
   {
      SCIP_CALL( SCIPrelaxInitsol(set->relaxs[i], set) );
   }

   /* separators */
   for( i = 0; i < set->nsepas; ++i )
   {
      SCIP_CALL( SCIPsepaInitsol(set->sepas[i], set) );
   }

   /* propagators */
   for( i = 0; i < set->nprops; ++i )
   {
      SCIP_CALL( SCIPpropInitsol(set->props[i], set) );
   }

   /* primal heuristics */
   for( i = 0; i < set->nheurs; ++i )
   {
      SCIP_CALL( SCIPheurInitsol(set->heurs[i], set) );
   }

   /* event handlers */
   for( i = 0; i < set->neventhdlrs; ++i )
   {
      SCIP_CALL( SCIPeventhdlrInitsol(set->eventhdlrs[i], set) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      SCIP_CALL( SCIPnodeselInitsol(set->nodesels[i], set) );
   }

   /* branching rules */
   for( i = 0; i < set->nbranchrules; ++i )
   {
      SCIP_CALL( SCIPbranchruleInitsol(set->branchrules[i], set) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      SCIP_CALL( SCIPdispInitsol(set->disps[i], set) );
   }

   /* statistics tables */
   for( i = 0; i < set->ntables; ++i )
   {
      SCIP_CALL( SCIPtableInitsol(set->tables[i], set) );
   }

   return SCIP_OKAY;
}

/** calls exitsol methods of all plugins */
SCIP_RETCODE SCIPsetExitsolPlugins(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   )
{
   int i;

   assert(set != NULL);

   /* active variable pricers */
   SCIPsetSortPricers(set);
   for( i = 0; i < set->nactivepricers; ++i )
   {
      SCIP_CALL( SCIPpricerExitsol(set->pricers[i], set) );
   }

   /* constraint handlers */
   for( i = 0; i < set->nconshdlrs; ++i )
   {
      SCIP_CALL( SCIPconshdlrExitsol(set->conshdlrs[i], blkmem, set, stat, restart) );
   }

   /* conflict handlers */
   for( i = 0; i < set->nconflicthdlrs; ++i )
   {
      SCIP_CALL( SCIPconflicthdlrExitsol(set->conflicthdlrs[i], set) );
   }

   /* relaxators */
   for( i = 0; i < set->nrelaxs; ++i )
   {
      SCIP_CALL( SCIPrelaxExitsol(set->relaxs[i], set) );
   }

   /* separators */
   for( i = 0; i < set->nsepas; ++i )
   {
      SCIP_CALL( SCIPsepaExitsol(set->sepas[i], set) );
   }

   /* propagators */
   for( i = 0; i < set->nprops; ++i )
   {
      SCIP_CALL( SCIPpropExitsol(set->props[i], set, restart) );
   }

   /* primal heuristics */
   for( i = 0; i < set->nheurs; ++i )
   {
      SCIP_CALL( SCIPheurExitsol(set->heurs[i], set) );
   }

   /* event handlers */
   for( i = 0; i < set->neventhdlrs; ++i )
   {
      SCIP_CALL( SCIPeventhdlrExitsol(set->eventhdlrs[i], set) );
   }

   /* node selectors */
   for( i = 0; i < set->nnodesels; ++i )
   {
      SCIP_CALL( SCIPnodeselExitsol(set->nodesels[i], set) );
   }

   /* branching rules */
   for( i = 0; i < set->nbranchrules; ++i )
   {
      SCIP_CALL( SCIPbranchruleExitsol(set->branchrules[i], set) );
   }

   /* display columns */
   for( i = 0; i < set->ndisps; ++i )
   {
      SCIP_CALL( SCIPdispExitsol(set->disps[i], set) );
   }

   /* statistics tables */
   for( i = 0; i < set->ntables; ++i )
   {
      SCIP_CALL( SCIPtableExitsol(set->tables[i], set) );
   }

   return SCIP_OKAY;
}

/** calculate memory size for dynamically allocated arrays */
int SCIPsetCalcMemGrowSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->mem_arraygrowinit, set->mem_arraygrowfac, num);
}

/** calculate memory size for tree array */
int SCIPsetCalcTreeGrowSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->mem_treegrowinit, set->mem_treegrowfac, num);
}

/** calculate memory size for path array */
int SCIPsetCalcPathGrowSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   return calcGrowSize(set->mem_pathgrowinit, set->mem_pathgrowfac, num);
}

/** sets verbosity level for message output */
SCIP_RETCODE SCIPsetSetVerbLevel(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VERBLEVEL        verblevel           /**< verbosity level for message output */
   )
{
   assert(set != NULL);

   if( verblevel > SCIP_VERBLEVEL_FULL )
   {
      SCIPerrorMessage("invalid verbosity level <%d>, maximum is <%d>\n", verblevel, SCIP_VERBLEVEL_FULL);
      return SCIP_INVALIDCALL;
   }

   set->disp_verblevel = verblevel;

   return SCIP_OKAY;
}

/** sets feasibility tolerance */
SCIP_RETCODE SCIPsetSetFeastol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             feastol             /**< new feasibility tolerance */
   )
{
   assert(set != NULL);

   set->num_feastol = feastol;

   /* the feasibility tolerance of the LP solver should never be larger than SCIP's feasibility tolerance; if necessary,
    * decrease it; use the SCIP change method in order to mark the LP unsolved
    */
   if( SCIPsetFeastol(set) < SCIPsetLpfeastol(set) )
   {
      SCIPsetDebugMsg(set, "decreasing lpfeastol along with feastol to %g\n", SCIPsetFeastol(set));
      SCIP_CALL( SCIPchgLpfeastol(set->scip, SCIPsetFeastol(set), TRUE) );
   }

   return SCIP_OKAY;
}

/** sets primal feasibility tolerance of LP solver */
SCIP_RETCODE SCIPsetSetLpfeastol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             lpfeastol,          /**< new primal feasibility tolerance of LP solver */
   SCIP_Bool             printnewvalue       /**< should "numerics/lpfeastol = ..." be printed? */
   )
{
   SCIP_RETCODE retcode;

   assert(set != NULL);

   retcode = SCIP_OKAY;

   /* the feasibility tolerance of the LP solver should never be larger than SCIP's feasibility tolerance; if this is
    * tried, we correct it to feastol; note that when we are called, e.g., by paramChgdLpfeastol, lpfeastol has already
    * been modified and so we cannot leave the lpfeastol value unchanged; if we would not return SCIP_PARAMETERWRONGVAL
    * in this case, the interactive shell would print the incorrect value to be set
    */
   if( lpfeastol > SCIPsetFeastol(set) )
   {
      SCIPerrorMessage("LP feasibility tolerance must be at least as tight as SCIP's feasibility tolerance\n");

      retcode = SCIP_PARAMETERWRONGVAL;
      printnewvalue = TRUE;

      set->num_lpfeastol = SCIPsetFeastol(set);
   }
   else
      set->num_lpfeastol = lpfeastol;

   if( printnewvalue )
   {
      SCIPverbMessage(set->scip, SCIP_VERBLEVEL_HIGH, NULL, "numerics/lpfeastol = %.15g\n", SCIPsetLpfeastol(set));
   }

   return retcode;
}

/** sets feasibility tolerance for reduced costs in LP solution */
SCIP_RETCODE SCIPsetSetDualfeastol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             dualfeastol         /**< new reduced costs feasibility tolerance */
   )
{
   assert(set != NULL);

   set->num_dualfeastol = dualfeastol;

   return SCIP_OKAY;
}

/** sets LP convergence tolerance used in barrier algorithm */
SCIP_RETCODE SCIPsetSetBarrierconvtol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   )
{
   assert(set != NULL);

   set->num_barrierconvtol = barrierconvtol;

   return SCIP_OKAY;
}

/** sets primal feasibility tolerance for relaxations (relaxfeastol)
 *
 * @note Set to SCIP_INVALID to apply relaxation-specific feasibility tolerance only.
 *
 * @return Previous value of relaxfeastol.
 */
SCIP_Real SCIPsetSetRelaxfeastol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             relaxfeastol        /**< new primal feasibility tolerance for relaxations, or SCIP_INVALID */
   )
{
   SCIP_Real oldval;

   assert(set != NULL);
   assert(relaxfeastol >= 0.0);

   oldval = set->num_relaxfeastol;
   set->num_relaxfeastol = relaxfeastol;

   return oldval;
}

/** marks that some limit parameter was changed */
void SCIPsetSetLimitChanged(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   set->limitchanged = TRUE;

   set->istimelimitfinite = (set->limit_time < SCIP_DEFAULT_LIMIT_TIME);
}

/** returns the maximal number of variables priced into the LP per round */
int SCIPsetGetPriceMaxvars(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   assert(set != NULL);

   if( root )
      return set->price_maxvarsroot;
   else
      return set->price_maxvars;
}

/** returns the maximal number of cuts separated per round */
int SCIPsetGetSepaMaxcuts(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             root                /**< are we at the root node? */
   )
{
   assert(set != NULL);

   if( root )
      return set->sepa_maxcutsroot;
   else
      return set->sepa_maxcuts;
}

/** returns user defined objective value (in original space) for reference purposes */
SCIP_Real SCIPsetGetReferencevalue(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(NULL != set);

   return set->misc_referencevalue;
}


/** returns debug solution data */
SCIP_DEBUGSOLDATA* SCIPsetGetDebugSolData(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->debugsoldata;
}


/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPsetInfinity
#undef SCIPsetEpsilon
#undef SCIPsetSumepsilon
#undef SCIPsetFeastol
#undef SCIPsetLpfeastol
#undef SCIPsetDualfeastol
#undef SCIPsetBarrierconvtol
#undef SCIPsetPseudocosteps
#undef SCIPsetPseudocostdelta
#undef SCIPsetCutoffbounddelta
#undef SCIPsetRelaxfeastol
#undef SCIPsetRecompfac
#undef SCIPsetIsEQ
#undef SCIPsetIsLT
#undef SCIPsetIsLE
#undef SCIPsetIsGT
#undef SCIPsetIsGE
#undef SCIPsetIsInfinity
#undef SCIPsetIsZero
#undef SCIPsetIsPositive
#undef SCIPsetIsNegative
#undef SCIPsetIsIntegral
#undef SCIPsetIsScalingIntegral
#undef SCIPsetIsFracIntegral
#undef SCIPsetFloor
#undef SCIPsetCeil
#undef SCIPsetRound
#undef SCIPsetFrac
#undef SCIPsetIsSumEQ
#undef SCIPsetIsSumLT
#undef SCIPsetIsSumLE
#undef SCIPsetIsSumGT
#undef SCIPsetIsSumGE
#undef SCIPsetIsSumZero
#undef SCIPsetIsSumPositive
#undef SCIPsetIsSumNegative
#undef SCIPsetSumFloor
#undef SCIPsetSumCeil
#undef SCIPsetSumRound
#undef SCIPsetSumFrac
#undef SCIPsetIsFeasEQ
#undef SCIPsetIsFeasLT
#undef SCIPsetIsFeasLE
#undef SCIPsetIsFeasGT
#undef SCIPsetIsFeasGE
#undef SCIPsetIsFeasZero
#undef SCIPsetIsFeasPositive
#undef SCIPsetIsFeasNegative
#undef SCIPsetIsFeasIntegral
#undef SCIPsetIsFeasFracIntegral
#undef SCIPsetFeasFloor
#undef SCIPsetFeasCeil
#undef SCIPsetFeasRound
#undef SCIPsetFeasFrac
#undef SCIPsetIsDualfeasEQ
#undef SCIPsetIsDualfeasLT
#undef SCIPsetIsDualfeasLE
#undef SCIPsetIsDualfeasGT
#undef SCIPsetIsDualfeasGE
#undef SCIPsetIsDualfeasZero
#undef SCIPsetIsDualfeasPositive
#undef SCIPsetIsDualfeasNegative
#undef SCIPsetIsDualfeasIntegral
#undef SCIPsetIsDualfeasFracIntegral
#undef SCIPsetDualfeasFloor
#undef SCIPsetDualfeasCeil
#undef SCIPsetDualfeasRound
#undef SCIPsetDualfeasFrac
#undef SCIPsetIsLbBetter
#undef SCIPsetIsUbBetter
#undef SCIPsetIsEfficacious
#undef SCIPsetIsRelEQ
#undef SCIPsetIsRelLT
#undef SCIPsetIsRelLE
#undef SCIPsetIsRelGT
#undef SCIPsetIsRelGE
#undef SCIPsetIsSumRelEQ
#undef SCIPsetIsSumRelLT
#undef SCIPsetIsSumRelLE
#undef SCIPsetIsSumRelGT
#undef SCIPsetIsSumRelGE
#undef SCIPsetIsUpdateUnreliable
#undef SCIPsetInitializeRandomSeed
#undef SCIPsetIsHugeValue
#undef SCIPsetGetHugeValue

/** returns value treated as infinity */
SCIP_Real SCIPsetInfinity(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_infinity;
}

/** returns the minimum value that is regarded as huge and should be handled separately (e.g., in activity
 *  computation)
 */
SCIP_Real SCIPsetGetHugeValue(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_hugeval;
}

/** returns value treated as zero */
SCIP_Real SCIPsetEpsilon(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_epsilon;
}

/** returns value treated as zero for sums of floating point values */
SCIP_Real SCIPsetSumepsilon(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_sumepsilon;
}

/** returns feasibility tolerance for constraints */
SCIP_Real SCIPsetFeastol(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_feastol;
}

/** returns feasibility tolerance for reduced costs */
SCIP_Real SCIPsetDualfeastol(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_dualfeastol;
}

/** returns primal feasibility tolerance of LP solver given as minimum of lpfeastol option and relaxfeastol */
SCIP_Real SCIPsetLpfeastol(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( set->num_relaxfeastol != SCIP_INVALID ) /*lint !e777*/
      return MIN(set->num_relaxfeastol, set->num_lpfeastol);

   return set->num_lpfeastol;
}

/** returns convergence tolerance used in barrier algorithm */
SCIP_Real SCIPsetBarrierconvtol(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_barrierconvtol;
}

/** returns minimal variable distance value to use for pseudo cost updates */
SCIP_Real SCIPsetPseudocosteps(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_pseudocosteps;
}

/** returns minimal minimal objective distance value to use for pseudo cost updates */
SCIP_Real SCIPsetPseudocostdelta(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_pseudocostdelta;
}

/** return the delta to use for computing the cutoff bound for integral objectives */
SCIP_Real SCIPsetCutoffbounddelta(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_Real feastol;

   assert(set != NULL);

   feastol = SCIPsetFeastol(set);

   return MIN(100.0 * feastol, 0.0001);
}

/** return the primal feasibility tolerance for relaxations */
SCIP_Real SCIPsetRelaxfeastol(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_relaxfeastol;
}

/** returns minimal decrease factor that causes the recomputation of a value
 *  (e.g., pseudo objective) instead of an update */
SCIP_Real SCIPsetRecompfac(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   return set->num_recompfac;
}

/** checks, if value is (positive) infinite */
SCIP_Bool SCIPsetIsInfinity(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be compared against infinity */
   )
{
   assert(set != NULL);

   return (val >= set->num_infinity);
}

/** checks, if value is huge and should be handled separately (e.g., in activity computation) */
SCIP_Bool SCIPsetIsHugeValue(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to be checked whether it is huge */
   )
{
   assert(set != NULL);

   return (val >= set->num_hugeval);
}

/** checks, if values are in range of epsilon */
SCIP_Bool SCIPsetIsEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSEQ(val1, val2, set->num_epsilon);
}

/** checks, if val1 is (more than epsilon) lower than val2 */
SCIP_Bool SCIPsetIsLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSLT(val1, val2, set->num_epsilon);
}

/** checks, if val1 is not (more than epsilon) greater than val2 */
SCIP_Bool SCIPsetIsLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSLE(val1, val2, set->num_epsilon);
}

/** checks, if val1 is (more than epsilon) greater than val2 */
SCIP_Bool SCIPsetIsGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSGT(val1, val2, set->num_epsilon);
}

/** checks, if val1 is not (more than epsilon) lower than val2 */
SCIP_Bool SCIPsetIsGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSGE(val1, val2, set->num_epsilon);
}

/** checks, if value is in range epsilon of 0.0 */
SCIP_Bool SCIPsetIsZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSZ(val, set->num_epsilon);
}

/** checks, if value is greater than epsilon */
SCIP_Bool SCIPsetIsPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSP(val, set->num_epsilon);
}

/** checks, if value is lower than -epsilon */
SCIP_Bool SCIPsetIsNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSN(val, set->num_epsilon);
}

/** checks, if value is integral within epsilon */
SCIP_Bool SCIPsetIsIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSISINT(val, set->num_epsilon);
}

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
SCIP_Bool SCIPsetIsScalingIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val,                /**< unscaled value to check for scaled integrality */
   SCIP_Real             scalar              /**< value to scale val with for checking for integrality */
   )
{
   SCIP_Real scaledeps;

   assert(set != NULL);

   scaledeps = REALABS(scalar);
   scaledeps = MAX(scaledeps, 1.0);
   scaledeps *= set->num_epsilon;

   return EPSISINT(scalar*val, scaledeps);
}

/** checks, if given fractional part is smaller than epsilon */
SCIP_Bool SCIPsetIsFracIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);
   assert(SCIPsetIsGE(set, val, -set->num_epsilon));
   assert(SCIPsetIsLE(set, val, 1.0+set->num_epsilon));

   return (val <= set->num_epsilon);
}

/** rounds value + feasibility tolerance down to the next integer in epsilon tolerance */
SCIP_Real SCIPsetFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSFLOOR(val, set->num_epsilon);
}

/** rounds value - feasibility tolerance up to the next integer in epsilon tolerance */
SCIP_Real SCIPsetCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSCEIL(val, set->num_epsilon);
}

/** rounds value to the nearest integer in epsilon tolerance */
SCIP_Real SCIPsetRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSROUND(val, set->num_epsilon);
}

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
SCIP_Real SCIPsetFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to return fractional part for */
   )
{
   assert(set != NULL);

   return EPSFRAC(val, set->num_epsilon);
}

/** checks, if values are in range of sumepsilon */
SCIP_Bool SCIPsetIsSumEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSEQ(val1, val2, set->num_sumepsilon);
}

/** checks, if val1 is (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPsetIsSumLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSLT(val1, val2, set->num_sumepsilon);
}

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPsetIsSumLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSLE(val1, val2, set->num_sumepsilon);
}

/** checks, if val1 is (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPsetIsSumGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSGT(val1, val2, set->num_sumepsilon);
}

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPsetIsSumGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return EPSGE(val1, val2, set->num_sumepsilon);
}

/** checks, if value is in range sumepsilon of 0.0 */
SCIP_Bool SCIPsetIsSumZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSZ(val, set->num_sumepsilon);
}

/** checks, if value is greater than sumepsilon */
SCIP_Bool SCIPsetIsSumPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSP(val, set->num_sumepsilon);
}

/** checks, if value is lower than -sumepsilon */
SCIP_Bool SCIPsetIsSumNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSN(val, set->num_sumepsilon);
}

/** rounds value + sumepsilon tolerance down to the next integer */
SCIP_Real SCIPsetSumFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSFLOOR(val, set->num_sumepsilon);
}

/** rounds value - sumepsilon tolerance up to the next integer */
SCIP_Real SCIPsetSumCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSCEIL(val, set->num_sumepsilon);
}

/** rounds value to the nearest integer in sumepsilon tolerance */
SCIP_Real SCIPsetSumRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSROUND(val, set->num_sumepsilon);
}

/** returns fractional part of value, i.e. x - floor(x) in sumepsilon tolerance */
SCIP_Real SCIPsetSumFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSFRAC(val, set->num_sumepsilon);
}

/** checks, if relative difference of values is in range of feastol */
SCIP_Bool SCIPsetIsFeasEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSZ(diff, set->num_feastol);
}

/** checks, if relative difference of val1 and val2 is lower than feastol */
SCIP_Bool SCIPsetIsFeasLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSN(diff, set->num_feastol);
}

/** checks, if relative difference of val1 and val2 is not greater than feastol */
SCIP_Bool SCIPsetIsFeasLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSP(diff, set->num_feastol);
}

/** checks, if relative difference of val1 and val2 is greater than feastol */
SCIP_Bool SCIPsetIsFeasGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSP(diff, set->num_feastol);
}

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
SCIP_Bool SCIPsetIsFeasGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSN(diff, set->num_feastol);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
SCIP_Bool SCIPsetIsFeasZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSZ(val, set->num_feastol);
}

/** checks, if value is greater than feasibility tolerance */
SCIP_Bool SCIPsetIsFeasPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSP(val, set->num_feastol);
}

/** checks, if value is lower than -feasibility tolerance */
SCIP_Bool SCIPsetIsFeasNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSN(val, set->num_feastol);
}

/** checks, if value is integral within the feasibility bounds */
SCIP_Bool SCIPsetIsFeasIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSISINT(val, set->num_feastol);
}

/** checks, if given fractional part is smaller than feastol */
SCIP_Bool SCIPsetIsFeasFracIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);
   assert(SCIPsetIsGE(set, val, -2*set->num_feastol));
   assert(SCIPsetIsLE(set, val, 1.0+set->num_feastol));

   return (val <= set->num_feastol);
}

/** rounds value + feasibility tolerance down to the next integer in feasibility tolerance */
SCIP_Real SCIPsetFeasFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSFLOOR(val, set->num_feastol);
}

/** rounds value - feasibility tolerance up to the next integer in feasibility tolerance */
SCIP_Real SCIPsetFeasCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSCEIL(val, set->num_feastol);
}

/** rounds value to the nearest integer in feasibility tolerance */
SCIP_Real SCIPsetFeasRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSROUND(val, set->num_feastol);
}

/** returns fractional part of value, i.e. x - floor(x) in feasibility tolerance */
SCIP_Real SCIPsetFeasFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSFRAC(val, set->num_feastol);
}

/** checks, if relative difference of values is in range of dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSZ(diff, set->num_dualfeastol);
}

/** checks, if relative difference of val1 and val2 is lower than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSN(diff, set->num_dualfeastol);
}

/** checks, if relative difference of val1 and val2 is not greater than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSP(diff, set->num_dualfeastol);
}

/** checks, if relative difference of val1 and val2 is greater than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSP(diff, set->num_dualfeastol);
}

/** checks, if relative difference of val1 and val2 is not lower than -dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSN(diff, set->num_dualfeastol);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
SCIP_Bool SCIPsetIsDualfeasZero(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSZ(val, set->num_dualfeastol);
}

/** checks, if value is greater than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasPositive(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSP(val, set->num_dualfeastol);
}

/** checks, if value is lower than -dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasNegative(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSN(val, set->num_dualfeastol);
}

/** checks, if value is integral within the dual feasibility bounds */
SCIP_Bool SCIPsetIsDualfeasIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSISINT(val, set->num_dualfeastol);
}

/** checks, if given fractional part is smaller than dual feasibility tolerance */
SCIP_Bool SCIPsetIsDualfeasFracIntegral(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);
   assert(SCIPsetIsGE(set, val, -set->num_dualfeastol));
   assert(SCIPsetIsLE(set, val, 1.0+set->num_dualfeastol));

   return (val <= set->num_dualfeastol);
}

/** rounds value + dual feasibility tolerance down to the next integer */
SCIP_Real SCIPsetDualfeasFloor(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSFLOOR(val, set->num_dualfeastol);
}

/** rounds value - dual feasibility tolerance up to the next integer */
SCIP_Real SCIPsetDualfeasCeil(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSCEIL(val, set->num_dualfeastol);
}

/** rounds value to the nearest integer in dual feasibility tolerance */
SCIP_Real SCIPsetDualfeasRound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSROUND(val, set->num_dualfeastol);
}

/** returns fractional part of value, i.e. x - floor(x) in dual feasibility tolerance */
SCIP_Real SCIPsetDualfeasFrac(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(set != NULL);

   return EPSFRAC(val, set->num_dualfeastol);
}

/** checks, if the given new lower bound is at least min(oldub - oldlb, |oldlb|) times the bound
 *  strengthening epsilon better than the old one or the change in the lower bound would fix the
 *  sign of the variable
 */
SCIP_Bool SCIPsetIsLbBetter(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   SCIP_Real eps;

   assert(set != NULL);
   assert(SCIPsetIsLE(set, oldlb, oldub));

   /* if lower bound is moved to 0 or higher, always accept bound change */
   if( oldlb < 0.0 && newlb >= 0.0 )
      return TRUE;

   eps = REALABS(oldlb);
   eps = MIN(oldub - oldlb, eps);
   return EPSGT(newlb, oldlb, set->num_boundstreps * MAX(eps, 1e-3));
}

/** checks, if the given new upper bound is at least min(oldub - oldlb, |oldub|) times the bound
 *  strengthening epsilon better than the old one or the change in the upper bound would fix the
 *  sign of the variable
 */
SCIP_Bool SCIPsetIsUbBetter(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   SCIP_Real eps;

   assert(set != NULL);
   assert(SCIPsetIsLE(set, oldlb, oldub));

   /* if upper bound is moved to 0 or lower, always accept bound change */
   if( oldub > 0.0 && newub <= 0.0 )
      return TRUE;

   eps = REALABS(oldub);
   eps = MIN(oldub - oldlb, eps);
   return EPSLT(newub, oldub, set->num_boundstreps * MAX(eps, 1e-3));
}

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy */
SCIP_Bool SCIPsetIsEfficacious(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             root,               /**< should the root's minimal cut efficacy be used? */
   SCIP_Real             efficacy            /**< efficacy of the cut */
   )
{
   assert(set != NULL);

   if( root )
      return EPSP(efficacy, set->sepa_minefficacyroot);
   else
      return EPSP(efficacy, set->sepa_minefficacy);
}

/** checks, if relative difference of values is in range of epsilon */
SCIP_Bool SCIPsetIsRelEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSZ(diff, set->num_epsilon);
}

/** checks, if relative difference of val1 and val2 is lower than epsilon */
SCIP_Bool SCIPsetIsRelLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSN(diff, set->num_epsilon);
}

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
SCIP_Bool SCIPsetIsRelLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSP(diff, set->num_epsilon);
}

/** checks, if relative difference of val1 and val2 is greater than epsilon */
SCIP_Bool SCIPsetIsRelGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSP(diff, set->num_epsilon);
}

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
SCIP_Bool SCIPsetIsRelGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSN(diff, set->num_epsilon);
}

/** checks, if relative difference of values is in range of sumepsilon */
SCIP_Bool SCIPsetIsSumRelEQ(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSZ(diff, set->num_sumepsilon);
}

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
SCIP_Bool SCIPsetIsSumRelLT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSN(diff, set->num_sumepsilon);
}

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
SCIP_Bool SCIPsetIsSumRelLE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSP(diff, set->num_sumepsilon);
}

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
SCIP_Bool SCIPsetIsSumRelGT(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return EPSP(diff, set->num_sumepsilon);
}

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
SCIP_Bool SCIPsetIsSumRelGE(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   SCIP_Real diff;

   assert(set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPsetIsInfinity(set, val1) || !SCIPsetIsInfinity(set, val2))
         && (!SCIPsetIsInfinity(set, -val1) || !SCIPsetIsInfinity(set, -val2)))
      || val1 == val2 );    /*lint !e777*/

   diff = SCIPrelDiff(val1, val2);

   return !EPSN(diff, set->num_sumepsilon);
}

/** Checks, if an iteratively updated value is reliable or should be recomputed from scratch.
 *  This is useful, if the value, e.g., the activity of a linear constraint or the pseudo objective value, gets a high
 *  absolute value during the optimization process which is later reduced significantly. In this case, the last digits
 *  were canceled out when increasing the value and are random after decreasing it.
 *  We dot not consider the cancellations which can occur during increasing the absolute value because they just cannot
 *  be expressed using fixed precision floating point arithmetic, anymore.
 *  The idea to get more reliable values is to always store the last reliable value, where increasing the absolute of
 *  the value is viewed as preserving reliability. Then, after each update, the new absolute value can be compared
 *  against the last reliable one with this method, checking whether it was decreased by a factor of at least
 *  "lp/recompfac" and should be recomputed.
 */
SCIP_Bool SCIPsetIsUpdateUnreliable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             newvalue,           /**< new value after update */
   SCIP_Real             oldvalue            /**< old value, i.e., last reliable value */
   )
{
   SCIP_Real quotient;

   assert(set != NULL);

   quotient = ABS(oldvalue) / MAX(ABS(newvalue), set->num_epsilon);

   return quotient >= set->num_recompfac;
}

/** prints a debug message */
void SCIPsetPrintDebugMessage(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline,         /**< line in the source file where the function was called */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   int subscipdepth = 0;
   SCIP* scip;
   va_list ap;

   assert( sourcefile != NULL );
   assert( set != NULL );

   scip = set->scip;
   assert( scip != NULL );

   if ( scip->stat != NULL )
      subscipdepth = scip->stat->subscipdepth;

   if ( subscipdepth > 0 )
      SCIPmessageFPrintInfo(scip->messagehdlr, NULL, "%d: [%s:%d] debug: ", subscipdepth, sourcefile, sourceline);
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, NULL, "[%s:%d] debug: ", sourcefile, sourceline);

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(scip->messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a debug message without precode */
void SCIPsetDebugMessagePrint(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   assert( set != NULL );
   assert( set->scip != NULL );

   va_start(ap, formatstr); /*lint !e838*/
   SCIPmessageVFPrintInfo(set->scip->messagehdlr, NULL, formatstr, ap);
   va_end(ap);
}

/** modifies an initial seed value with the global shift of random seeds */
int SCIPsetInitializeRandomSeed(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   initialseedvalue    /**< initial seed value to be modified */
   )
{
   assert(set != NULL);

   return (initialseedvalue + set->random_randomseedshift);
}
