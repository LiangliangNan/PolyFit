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

/**@file   cons_nonlinear.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for nonlinear constraints specified by algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef SCIP_DEBUG
#define ENFO_LOGGING
#endif

/* enable to get log output for enforcement */
/* #define ENFO_LOGGING */
/* define to get enforcement logging into file */
/* #define ENFOLOGFILE "consexpr_enfo.log" */

/* define to get more debug output from domain propagation */
/* #define DEBUG_PROP */

/*lint -e440*/
/*lint -e441*/
/*lint -e528*/
/*lint -e666*/
/*lint -e777*/
/*lint -e866*/

#include <assert.h>
#include <ctype.h>

#include "scip/cons_nonlinear.h"
#include "scip/nlhdlr.h"
#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_value.h"
#include "scip/expr_pow.h"
#include "scip/nlhdlr_convex.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/nlpi_ipopt.h"  /* for SCIPsolveLinearEquationsIpopt */
#include "scip/debug.h"
#include "scip/dialog_default.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "nonlinear"
#define CONSHDLR_DESC          "handler for nonlinear constraints specified by algebraic expressions"
#define CONSHDLR_ENFOPRIORITY       -60 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000010 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_ALWAYS /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */

/* properties of the nonlinear constraint handler statistics table */
#define TABLE_NAME_NONLINEAR           "cons_nonlinear"
#define TABLE_DESC_NONLINEAR           "nonlinear constraint handler statistics"
#define TABLE_POSITION_NONLINEAR       14600                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_NONLINEAR SCIP_STAGE_TRANSFORMED /**< output of the statistics table is only printed from this stage onwards */

/* properties of the nonlinear handler statistics table */
#define TABLE_NAME_NLHDLR              "nlhdlr"
#define TABLE_DESC_NLHDLR              "nonlinear handler statistics"
#define TABLE_POSITION_NLHDLR          14601                  /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE_NLHDLR    SCIP_STAGE_PRESOLVING  /**< output of the statistics table is only printed from this stage onwards */

#define DIALOG_NAME            "nlhdlrs"
#define DIALOG_DESC            "display nonlinear handlers"
#define DIALOG_ISSUBMENU          FALSE

#define VERTEXPOLY_MAXPERTURBATION      1e-3 /**< maximum perturbation */
#define VERTEXPOLY_USEDUALSIMPLEX       TRUE /**< use dual or primal simplex algorithm? */
#define VERTEXPOLY_RANDNUMINITSEED  20181029 /**< seed for random number generator, which is used to move points away from the boundary */
#define VERTEXPOLY_ADJUSTFACETFACTOR     1e1 /**< adjust resulting facets in checkRikun() up to a violation of this value times lpfeastol */

#define BRANCH_RANDNUMINITSEED      20191229 /**< seed for random number generator, which is used to select from several similar good branching candidates */

#define BILIN_MAXNAUXEXPRS                10 /**< maximal number of auxiliary expressions per bilinear term */

/** translate from one value of infinity to another
 *
 *  if val is &ge; infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/** translates x to 2^x for non-negative integer x */
#define POWEROFTWO(x) (0x1u << (x))

#ifdef ENFO_LOGGING
#define ENFOLOG(x) if( SCIPgetSubscipDepth(scip) == 0 && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL ) { x }
FILE* enfologfile = NULL;
#else
#define ENFOLOG(x)
#endif

/*
 * Data structures
 */

/** enforcement data of an expression */
typedef struct
{
   SCIP_NLHDLR*          nlhdlr;             /**< nonlinear handler */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata;     /**< data of nonlinear handler */
   SCIP_NLHDLR_METHOD nlhdlrparticipation; /**< methods where nonlinear handler participates */
   SCIP_Bool             issepainit;         /**< was the initsepa callback of nlhdlr called */
   SCIP_Real             auxvalue;           /**< auxiliary value of expression w.r.t. currently enforced solution */
   SCIP_Bool             sepabelowusesactivity;/**< whether sepabelow uses activity of some expression */
   SCIP_Bool             sepaaboveusesactivity;/**< whether sepaabove uses activity of some expression */
} EXPRENFO;

/** data stored by constraint handler in an expression that belongs to a nonlinear constraint */
struct SCIP_Expr_OwnerData
{
   SCIP_CONSHDLR*        conshdlr;           /** nonlinear constraint handler */

   /* locks and monotonicity */
   int                   nlockspos;          /**< positive locks counter */
   int                   nlocksneg;          /**< negative locks counter */
   SCIP_MONOTONE*        monotonicity;       /**< array containing monotonicity of expression w.r.t. each child */
   int                   monotonicitysize;   /**< length of monotonicity array */

   /* propagation (in addition to activity that is stored in expr) */
   SCIP_INTERVAL         propbounds;         /**< bounds to propagate in reverse propagation */
   unsigned int          propboundstag;      /**< tag to indicate whether propbounds are valid for the current propagation rounds */
   SCIP_Bool             inpropqueue;        /**< whether expression is queued for propagation */

   /* enforcement of expr == auxvar (or expr <= auxvar, or expr >= auxvar) */
   EXPRENFO**            enfos;              /**< enforcements */
   int                   nenfos;             /**< number of enforcements, or -1 if not initialized */
   unsigned int          lastenforced;       /**< last enforcement round where expression was enforced successfully */
   unsigned int          nactivityusesprop;  /**< number of nonlinear handlers whose activity computation (or domain propagation) depends on the activity of the expression */
   unsigned int          nactivityusessepa;  /**< number of nonlinear handlers whose separation (estimate or enfo) depends on the activity of the expression */
   unsigned int          nauxvaruses;        /**< number of nonlinear handlers whose separation uses an auxvar in the expression */
   SCIP_VAR*             auxvar;             /**< auxiliary variable used for outer approximation cuts */

   /* branching */
   SCIP_Real             violscoresum;       /**< sum of violation scores for branching stored for this expression */
   SCIP_Real             violscoremax;       /**< max of violation scores for branching stored for this expression */
   int                   nviolscores;        /**< number of violation scores stored for this expression */
   unsigned int          violscoretag;       /**< tag to decide whether a violation score of an expression needs to be initialized */

   /* additional data for variable expressions (TODO move into sub-struct?) */
   SCIP_CONS**           conss;              /**< constraints in which this variable appears */
   int                   nconss;             /**< current number of constraints in conss */
   int                   consssize;          /**< length of conss array */
   SCIP_Bool             consssorted;        /**< is the array of constraints sorted */

   int                   filterpos;          /**< position of eventdata in SCIP's event filter, -1 if not catching events */
};

/** constraint data for nonlinear constraints */
struct SCIP_ConsData
{
   /* data that defines the constraint: expression and sides */
   SCIP_EXPR*            expr;               /**< expression that represents this constraint */
   SCIP_Real             lhs;                /**< left-hand side */
   SCIP_Real             rhs;                /**< right-hand side */

   /* variables */
   SCIP_EXPR**           varexprs;           /**< array containing all variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */
   SCIP_Bool             catchedevents;      /**< do we catch events on variables? */

   /* constraint violation */
   SCIP_Real             lhsviol;            /**< violation of left-hand side by current solution */
   SCIP_Real             rhsviol;            /**< violation of right-hand side by current solution */
   SCIP_Real             gradnorm;           /**< norm of gradient of constraint function in current solution (if evaluated) */
   SCIP_Longint          gradnormsoltag;     /**< tag of solution used that gradnorm corresponds to */

   /* status flags */
   unsigned int          ispropagated:1;     /**< did we propagate the current bounds already? */
   unsigned int          issimplified:1;     /**< did we simplify the expression tree already? */

   /* locks */
   int                   nlockspos;          /**< number of positive locks */
   int                   nlocksneg;          /**< number of negative locks */

   /* repair infeasible solutions */
   SCIP_VAR*             linvardecr;         /**< variable that may be decreased without making any other constraint infeasible, or NULL if none */
   SCIP_VAR*             linvarincr;         /**< variable that may be increased without making any other constraint infeasible, or NULL if none */
   SCIP_Real             linvardecrcoef;     /**< linear coefficient of linvardecr */
   SCIP_Real             linvarincrcoef;     /**< linear coefficient of linvarincr */

   /* miscellaneous */
   SCIP_EXPRCURV         curv;               /**< curvature of the root expression w.r.t. the original variables */
   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */
   int                   consindex;          /**< an index of the constraint that is unique among all expr-constraints in this SCIP instance and is constant */
};

/** constraint upgrade method */
typedef struct
{
   SCIP_DECL_NONLINCONSUPGD((*consupgd));    /**< method to call for upgrading nonlinear constraint */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
} CONSUPGRADE;

/** constraint handler data */
struct SCIP_ConshdlrData
{
   /* nonlinear handler */
   SCIP_NLHDLR**         nlhdlrs;            /**< nonlinear handlers */
   int                   nnlhdlrs;           /**< number of nonlinear handlers */
   int                   nlhdlrssize;        /**< size of nlhdlrs array */
   SCIP_Bool             indetect;           /**< whether we are currently in detectNlhdlr */
   SCIP_Bool             registerusesactivitysepabelow; /**< a flag that is used only during \ref @detectNlhdlr() */
   SCIP_Bool             registerusesactivitysepaabove; /**< a flag that is used only during \ref @detectNlhdlr() */

   /* constraint upgrades */
   CONSUPGRADE**         consupgrades;       /**< constraint upgrade methods for specializing nonlinear constraints */
   int                   consupgradessize;   /**< size of consupgrades array */
   int                   nconsupgrades;      /**< number of constraint upgrade methods */

   /* other plugins */
   SCIP_EVENTHDLR*       eventhdlr;          /**< handler for variable bound change events */
   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subnlp heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the trysol heuristic, if available */

   /* tags and counters */
   int                   auxvarid;           /**< unique id for the next auxiliary variable */
   SCIP_Longint          curboundstag;       /**< tag indicating current variable bounds */
   SCIP_Longint          lastboundrelax;     /**< tag when bounds where most recently relaxed */
   SCIP_Longint          lastvaractivitymethodchange; /**< tag when method used to evaluate activity of variables changed last */
   unsigned int          enforound;          /**< total number of enforcement calls, including current one */
   int                   lastconsindex;      /**< last used consindex, plus one */

   /* activity intervals and domain propagation */
   SCIP_DECL_EXPR_INTEVALVAR((*intevalvar)); /**< method currently used for activity calculation of variable expressions */
   SCIP_Bool             globalbounds;       /**< whether global variable bounds should be used for activity calculation */
   SCIP_QUEUE*           reversepropqueue;   /**< expression queue to be used in reverse propagation, filled by SCIPtightenExprIntervalNonlinear */
   SCIP_Bool             forceboundtightening; /**< whether bound change passed to SCIPtightenExprIntervalNonlinear should be forced */
   unsigned int          curpropboundstag;   /**< tag indicating current propagation rounds, to match with expr->propboundstag */

   /* parameters */
   int                   maxproprounds;      /**< limit on number of propagation rounds for a set of constraints within one round of SCIP propagation */
   SCIP_Bool             propauxvars;        /**< whether to check bounds of all auxiliary variable to seed reverse propagation */
   char                  varboundrelax;      /**< strategy on how to relax variable bounds during bound tightening */
   SCIP_Real             varboundrelaxamount; /**< by how much to relax variable bounds during bound tightening */
   SCIP_Real             conssiderelaxamount; /**< by how much to relax constraint sides during bound tightening */
   SCIP_Real             vp_maxperturb;      /**< maximal relative perturbation of reference point */
   SCIP_Real             vp_adjfacetthreshold; /**< adjust computed facet up to a violation of this value times lpfeastol */
   SCIP_Bool             vp_dualsimplex;     /**< whether to use dual simplex instead of primal simplex for facet computing LP */
   SCIP_Bool             reformbinprods;     /**< whether to reformulate products of binary variables during presolving */
   SCIP_Bool             reformbinprodsand;  /**< whether to use the AND constraint handler for reformulating binary products */
   int                   reformbinprodsfac;  /**< minimum number of terms to reformulate bilinear binary products by factorizing variables (<= 1: disabled) */
   SCIP_Bool             forbidmultaggrnlvar; /**< whether to forbid multiaggregation of variables that appear in a nonlinear term of a constraint */
   SCIP_Bool             tightenlpfeastol;   /**< whether to tighten LP feasibility tolerance during enforcement, if it seems useful */
   SCIP_Bool             propinenforce;      /**< whether to (re)run propagation in enforcement */
   SCIP_Real             weakcutthreshold;   /**< threshold for when to regard a cut from an estimator as weak */
   SCIP_Real             strongcutmaxcoef;   /**< "strong" cuts will be scaled to have their maximal coef in [1/strongcutmaxcoef,strongcutmaxcoef] */
   SCIP_Bool             strongcutefficacy;  /**< consider efficacy requirement when deciding whether a cut is "strong" */
   SCIP_Bool             forcestrongcut;     /**< whether to force "strong" cuts in enforcement */
   SCIP_Real             enfoauxviolfactor;  /**< an expression will be enforced if the "auxiliary" violation is at least enfoauxviolfactor times the "original" violation */
   SCIP_Real             weakcutminviolfactor; /**< retry with weak cuts for constraints with violation at least this factor of maximal violated constraints */
   char                  rownotremovable;    /**< whether to make rows to be non-removable in the node where they are added (can prevent some cycling): 'o'ff, in 'e'nforcement only, 'a'lways */
   char                  violscale;          /**< method how to scale violations to make them comparable (not used for feasibility check) */
   char                  checkvarlocks;      /**< whether variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction) */
   int                   branchauxmindepth;  /**< from which depth on to allow branching on auxiliary variables */
   SCIP_Bool             branchexternal;     /**< whether to use external branching candidates for branching */
   SCIP_Real             branchhighviolfactor; /**< consider a constraint highly violated if its violation is >= this factor * maximal violation among all constraints */
   SCIP_Real             branchhighscorefactor; /**< consider a variable branching score high if its branching score >= this factor * maximal branching score among all variables */
   SCIP_Real             branchviolweight;   /**< weight by how much to consider the violation assigned to a variable for its branching score */
   SCIP_Real             branchdualweight;   /**< weight by how much to consider the dual values of rows that contain a variable for its branching score */
   SCIP_Real             branchpscostweight; /**< weight by how much to consider the pseudo cost of a variable for its branching score */
   SCIP_Real             branchdomainweight; /**< weight by how much to consider the domain width in branching score */
   SCIP_Real             branchvartypeweight;/**< weight by how much to consider variable type in branching score */
   char                  branchscoreagg;     /**< how to aggregate several branching scores given for the same expression ('a'verage, 'm'aximum, or 's'um) */
   char                  branchviolsplit;    /**< method used to split violation in expression onto variables ('u'niform, 'm'idness of solution, 'd'omain width, 'l'ogarithmic domain width) */
   SCIP_Real             branchpscostreliable; /**< minimum pseudo-cost update count required to consider pseudo-costs reliable */
   char                  linearizeheursol;   /**< whether tight linearizations of nonlinear constraints should be added to cutpool when some heuristics finds a new solution ('o'ff, on new 'i'ncumbents, on 'e'very solution) */
   SCIP_Bool             assumeconvex;       /**< whether to assume that any constraint is convex */

   /* statistics */
   SCIP_Longint          nweaksepa;          /**< number of times we used "weak" cuts for enforcement */
   SCIP_Longint          ntightenlp;         /**< number of times we requested solving the LP with a smaller feasibility tolerance when enforcing */
   SCIP_Longint          ndesperatetightenlp; /**< number of times we requested solving the LP with a smaller feasibility tolerance when enforcing because we didn't know anything better */
   SCIP_Longint          ndesperatebranch;   /**< number of times we branched on some variable because normal enforcement was not successful */
   SCIP_Longint          ndesperatecutoff;   /**< number of times we cut off a node in enforcement because no branching candidate could be found */
   SCIP_Longint          nforcelp;           /**< number of times we forced solving the LP when enforcing a pseudo solution */
   SCIP_CLOCK*           canonicalizetime;   /**< time spend for canonicalization */
   SCIP_Longint          ncanonicalizecalls; /**< number of times we called canonicalization */

   /* facets of envelops of vertex-polyhedral functions */
   SCIP_RANDNUMGEN*      vp_randnumgen;      /**< random number generator used to perturb reference point */
   SCIP_LPI*             vp_lp[SCIP_MAXVERTEXPOLYDIM+1];  /**< LPs used to compute facets for functions of different dimension */

   /* hashing of bilinear terms */
   SCIP_HASHTABLE*       bilinhashtable;     /**< hash table for bilinear terms */
   SCIP_CONSNONLINEAR_BILINTERM* bilinterms; /**< bilinear terms */
   int                   nbilinterms;        /**< total number of bilinear terms */
   int                   bilintermssize;     /**< size of bilinterms array */
   int                   bilinmaxnauxexprs;  /**< maximal number of auxiliary expressions per bilinear term */

   /* branching */
   SCIP_RANDNUMGEN*      branchrandnumgen;   /**< random number generated used in branching variable selection */
   char                  branchpscostupdatestrategy; /**< value of parameter branching/lpgainnormalize */

   /* misc */
   SCIP_Bool             checkedvarlocks;    /**< whether variables contained in a single constraint have been already considered */
   SCIP_HASHMAP*         var2expr;           /**< hashmap to map SCIP variables to variable-expressions */
   int                   newsoleventfilterpos; /**< filter position of new solution event handler, if caught */
};

/** branching candidate with various scores */
typedef struct
{
   SCIP_EXPR*            expr;               /**< expression that holds branching candidate */
   SCIP_Real             auxviol;            /**< aux-violation score of candidate */
   SCIP_Real             domain;             /**< domain score of candidate */
   SCIP_Real             dual;               /**< dual score of candidate */
   SCIP_Real             pscost;             /**< pseudo-cost score of candidate */
   SCIP_Real             vartype;            /**< variable type score of candidate */
   SCIP_Real             weighted;           /**< weighted sum of other scores, see scoreBranchingCandidates() */
} BRANCHCAND;

/*
 * Local methods
 */

/* forward declaration */
static
SCIP_RETCODE forwardPropExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_EXPR*            rootexpr,           /**< expression */
   SCIP_Bool             tightenauxvars,     /**< should the bounds of auxiliary variables be tightened? */
   SCIP_Bool*            infeasible,         /**< buffer to store whether the problem is infeasible (NULL if not needed) */
   int*                  ntightenings        /**< buffer to store the number of auxiliary variable tightenings (NULL if not needed) */
   );

/** frees auxiliary variables of expression, if any */
static
SCIP_RETCODE freeAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression which auxvar to free, if any */
   )
{
   SCIP_EXPR_OWNERDATA* mydata;

   assert(scip != NULL);
   assert(expr != NULL);

   mydata = SCIPexprGetOwnerData(expr);
   assert(mydata != NULL);

   if( mydata->auxvar == NULL )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "remove auxiliary variable <%s> for expression %p\n", SCIPvarGetName(mydata->auxvar), (void*)expr);

   /* remove variable locks
    * as this is a relaxation-only variable, no other plugin should use it for deducing any type of reductions or cutting planes
    */
   SCIP_CALL( SCIPaddVarLocks(scip, mydata->auxvar, -1, -1) );

   /* release auxiliary variable */
   SCIP_CALL( SCIPreleaseVar(scip, &mydata->auxvar) );
   assert(mydata->auxvar == NULL);

   return SCIP_OKAY;
}

/** frees data used for enforcement of expression, that is, nonlinear handlers
 *
 * can also clear indicators whether expr needs enforcement methods, that is,
 * free an associated auxiliary variable and reset the nactivityuses counts
 */
static
SCIP_RETCODE freeEnfoData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression whose enforcement data will be released */
   SCIP_Bool             freeauxvar          /**< whether aux var should be released and activity usage counts be reset */
   )
{
   SCIP_EXPR_OWNERDATA* mydata;
   int e;

   mydata = SCIPexprGetOwnerData(expr);
   assert(mydata != NULL);

   if( freeauxvar )
   {
      /* free auxiliary variable */
      SCIP_CALL( freeAuxVar(scip, expr) );
      assert(mydata->auxvar == NULL);

      /* reset count on activity and auxvar usage */
      mydata->nactivityusesprop = 0;
      mydata->nactivityusessepa = 0;
      mydata->nauxvaruses = 0;
   }

   /* free data stored by nonlinear handlers */
   for( e = 0; e < mydata->nenfos; ++e )
   {
      SCIP_NLHDLR* nlhdlr;

      assert(mydata->enfos[e] != NULL);

      nlhdlr = mydata->enfos[e]->nlhdlr;
      assert(nlhdlr != NULL);

      if( mydata->enfos[e]->issepainit )
      {
         /* call the separation deinitialization callback of the nonlinear handler */
         SCIP_CALL( SCIPnlhdlrExitsepa(scip, nlhdlr, expr, mydata->enfos[e]->nlhdlrexprdata) );
         mydata->enfos[e]->issepainit = FALSE;
      }

      /* free nlhdlr exprdata, if there is any and there is a method to free this data */
      if( mydata->enfos[e]->nlhdlrexprdata != NULL )
      {
         SCIP_CALL( SCIPnlhdlrFreeexprdata(scip, nlhdlr, expr, &mydata->enfos[e]->nlhdlrexprdata) );
         assert(mydata->enfos[e]->nlhdlrexprdata == NULL);
      }

      /* free enfo data */
      SCIPfreeBlockMemory(scip, &mydata->enfos[e]);
   }

   /* free array with enfo data */
   SCIPfreeBlockMemoryArrayNull(scip, &mydata->enfos, mydata->nenfos);

   /* we need to look at this expression in detect again */
   mydata->nenfos = -1;

   return SCIP_OKAY;
}

/** callback that frees data that this conshdlr stored in an expression */
static
SCIP_DECL_EXPR_OWNERFREE(exprownerFree)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(ownerdata != NULL);
   assert(*ownerdata != NULL);

   /* expression should not be locked anymore */
   assert((*ownerdata)->nlockspos == 0);
   assert((*ownerdata)->nlocksneg == 0);

   SCIP_CALL( freeEnfoData(scip, expr, TRUE) );

   /* expression should not be enforced anymore */
   assert((*ownerdata)->nenfos <= 0);
   assert((*ownerdata)->auxvar == NULL);

   if( SCIPisExprVar(scip, expr) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_VAR* var;

      /* there should be no constraints left that still use this variable */
      assert((*ownerdata)->nconss == 0);
      /* thus, there should also be no variable event catched (via this exprhdlr) */
      assert((*ownerdata)->filterpos == -1);

      SCIPfreeBlockMemoryArrayNull(scip, &(*ownerdata)->conss, (*ownerdata)->consssize);

      /* update var2expr hashmap in conshdlrdata */
      conshdlrdata = SCIPconshdlrGetData((*ownerdata)->conshdlr);
      assert(conshdlrdata != NULL);

      var = SCIPgetVarExprVar(expr);
      assert(var != NULL);

      /* remove var -> expr map from hashmap if present
       *  (if no variable-expression stored for var hashmap, then the var hasn't been used in any constraint, so do nothing
       *   if variable-expression stored for var is different, then also do nothing)
       */
      if( SCIPhashmapGetImage(conshdlrdata->var2expr, var) == (void*)expr )
      {
         SCIP_CALL( SCIPhashmapRemove(conshdlrdata->var2expr, var) );
      }
   }

   SCIPfreeBlockMemory(scip, ownerdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPR_OWNERPRINT(exprownerPrint)
{  /*lint --e{715}*/
   assert(ownerdata != NULL);

   /* print nl handlers associated to expr */
   if( ownerdata->nenfos > 0 )
   {
      int i;
      SCIPinfoMessage(scip, file, "   {");

      for( i = 0; i < ownerdata->nenfos; ++i )
      {
         SCIPinfoMessage(scip, file, "%s:", SCIPnlhdlrGetName(ownerdata->enfos[i]->nlhdlr));
         if( ownerdata->enfos[i]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_ACTIVITY )
            SCIPinfoMessage(scip, file, "a");
         if( ownerdata->enfos[i]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABELOW )
            SCIPinfoMessage(scip, file, "u");
         if( ownerdata->enfos[i]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPAABOVE )
            SCIPinfoMessage(scip, file, "o");
         if( i < ownerdata->nenfos-1 )
            SCIPinfoMessage(scip, file, ", ");
      }

      SCIPinfoMessage(scip, file, "}");
   }

   /* print aux var associated to expr */
   if( ownerdata->auxvar != NULL )
   {
      SCIPinfoMessage(scip, file, "  (<%s> in [%g, %g])", SCIPvarGetName(ownerdata->auxvar), SCIPvarGetLbLocal(ownerdata->auxvar), SCIPvarGetUbLocal(ownerdata->auxvar));
   }
   SCIPinfoMessage(scip, file, "\n");

   return SCIP_OKAY;
}

/** possibly reevaluates and then returns the activity of the expression
 *
 * Reevaluate activity if currently stored is not up to date (some bound was changed since last evaluation).
 */
static
SCIP_DECL_EXPR_OWNEREVALACTIVITY(exprownerEvalactivity)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(ownerdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPexprGetActivityTag(expr) < conshdlrdata->curboundstag )
   {
      /* update activity of expression */
      SCIP_CALL( forwardPropExpr(scip, ownerdata->conshdlr, expr, FALSE, NULL, NULL) );

      assert(SCIPexprGetActivityTag(expr) == conshdlrdata->curboundstag);
   }

   return SCIP_OKAY;
}

/** callback that creates data that this conshdlr wants to store in an expression */
static
SCIP_DECL_EXPR_OWNERCREATE(exprownerCreate)
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(ownerdata != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, ownerdata) );
   (*ownerdata)->nenfos = -1;
   (*ownerdata)->conshdlr = (SCIP_CONSHDLR*)ownercreatedata;

   if( SCIPisExprVar(scip, expr) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_VAR* var;

      (*ownerdata)->filterpos = -1;

      /* add to var2expr hashmap if not having expr for var yet */

      conshdlrdata = SCIPconshdlrGetData((*ownerdata)->conshdlr);
      assert(conshdlrdata != NULL);

      var = SCIPgetVarExprVar(expr);

      if( !SCIPhashmapExists(conshdlrdata->var2expr, (void*)var) )
      {
         /* store the variable expression in the hashmap */
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->var2expr, (void*)var, (void*)expr) );
      }
      else
      {
         /* if expr was just created, then it shouldn't already be stored as image of var */
         assert(SCIPhashmapGetImage(conshdlrdata->var2expr, (void*)var) != (void*)expr);
      }
   }
   else
   {
      /* just so that we can use filterpos to recognize whether an expr is a varexpr if not having a SCIP pointer around */
      (*ownerdata)->filterpos = -2;
   }

   *ownerfree = exprownerFree;
   *ownerprint = exprownerPrint;
   *ownerevalactivity = exprownerEvalactivity;

   return SCIP_OKAY;
}

/** creates a variable expression or retrieves from hashmap in conshdlr data */
static
SCIP_RETCODE createExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_VAR*             var                 /**< variable to be stored */
   )
{
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(var != NULL);

   /* get variable expression representing the given variable if there is one already */
   *expr = (SCIP_EXPR*) SCIPhashmapGetImage(SCIPconshdlrGetData(conshdlr)->var2expr, (void*) var);

   if( *expr == NULL )
   {
      /* create a new variable expression; this also captures the expression */
      SCIP_CALL( SCIPcreateExprVar(scip, expr, var, exprownerCreate, (void*)conshdlr) );
      assert(*expr != NULL);
      /* exprownerCreate should have added var->expr to var2expr */
      assert(SCIPhashmapGetImage(SCIPconshdlrGetData(conshdlr)->var2expr, (void*)var) == (void*)*expr);
   }
   else
   {
      /* only capture already existing expr to get a consistent uses-count */
      SCIPcaptureExpr(*expr);
   }

   return SCIP_OKAY;
}

/* map var exprs to var-expr from var2expr hashmap */
static
SCIP_DECL_EXPR_MAPEXPR(mapexprvar)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr = (SCIP_CONSHDLR*)mapexprdata;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);
   assert(*targetexpr == NULL);
   assert(mapexprdata != NULL);

   /* do not provide map if not variable */
   if( !SCIPisExprVar(sourcescip, sourceexpr) )
      return SCIP_OKAY;

   SCIP_CALL( createExprVar(targetscip, conshdlr, targetexpr, SCIPgetVarExprVar(sourceexpr)) );

   return SCIP_OKAY;
}

/* map var exprs to var-expr from var2expr hashmap corresponding to transformed var */
static
SCIP_DECL_EXPR_MAPEXPR(mapexprtransvar)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr = (SCIP_CONSHDLR*)mapexprdata;
   SCIP_VAR* var;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);
   assert(*targetexpr == NULL);
   assert(mapexprdata != NULL);

   /* do not provide map if not variable */
   if( !SCIPisExprVar(sourcescip, sourceexpr) )
      return SCIP_OKAY;

   var = SCIPgetVarExprVar(sourceexpr);
   assert(var != NULL);

   /* transform variable */
   SCIP_CALL( SCIPgetTransformedVar(sourcescip, var, &var) );
   assert(var != NULL);

   SCIP_CALL( createExprVar(targetscip, conshdlr, targetexpr, var) );

   return SCIP_OKAY;
}

/** stores all variable expressions into a given constraint */
static
SCIP_RETCODE storeVarExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int varexprssize;
   int i;

   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already */
   if( consdata->varexprs != NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs == NULL);
   assert(consdata->nvarexprs == 0);

   /* get an upper bound on number of variable expressions */
   if( consdata->issimplified )
   {
      /* if simplified, then we should have removed inactive variables and replaced common subexpressions,
       * so we cannot have more variable expression than the number of active variables
       */
      varexprssize = SCIPgetNVars(scip);
   }
   else
   {
      SCIP_CALL( SCIPgetExprNVars(scip, consdata->expr, &varexprssize) );
   }

   /* create array to store all variable expressions */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->varexprs, varexprssize) );

   SCIP_CALL( SCIPgetExprVarExprs(scip, consdata->expr, consdata->varexprs, &(consdata->nvarexprs)) );
   assert(varexprssize >= consdata->nvarexprs);

   /* shrink array if there are less variables in the expression than in the problem */
   if( varexprssize > consdata->nvarexprs )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->varexprs, varexprssize, consdata->nvarexprs) );
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->var2expr != NULL);

   /* ensure that for every variable an entry exists in the var2expr hashmap
    * when removing duplicate subexpressions it can happen that a var->varexpr map was removed from the hashmap
    */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      if( !SCIPhashmapExists(conshdlrdata->var2expr, SCIPgetVarExprVar(consdata->varexprs[i])) )
      {
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->var2expr, SCIPgetVarExprVar(consdata->varexprs[i]), consdata->varexprs[i]) );
      }
   }

   return SCIP_OKAY;
}

/** frees all variable expression stored in storeVarExprs() */
static
SCIP_RETCODE freeVarExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int i;

   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already*/
   if( consdata->varexprs == NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* release variable expressions */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);
      SCIP_CALL( SCIPreleaseExpr(scip, &consdata->varexprs[i]) );
      assert(consdata->varexprs[i] == NULL);
   }

   /* free variable expressions */
   SCIPfreeBlockMemoryArrayNull(scip, &consdata->varexprs, consdata->nvarexprs);
   consdata->varexprs = NULL;
   consdata->nvarexprs = 0;

   return SCIP_OKAY;
}

/** interval evaluation of variables as used in bound tightening
 *
 * Returns slightly relaxed local variable bounds of a variable as interval.
 * Does not relax beyond integer values, thus does not relax bounds on integer variables at all.
 */
static
SCIP_DECL_EXPR_INTEVALVAR(intEvalVarBoundTightening)
{
   SCIP_INTERVAL interval;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)intevalvardata;
   assert(conshdlrdata != NULL);

   if( conshdlrdata->globalbounds )
   {
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
   }
   else
   {
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
   }
   assert(lb <= ub);  /* SCIP should ensure that variable bounds are not contradicting */

   /* implicit integer variables may have non-integer bounds, apparently (run space25a) */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
   {
      lb = EPSROUND(lb, 0.0); /*lint !e835*/
      ub = EPSROUND(ub, 0.0); /*lint !e835*/
   }

   /* integer variables should always have integral bounds in SCIP */
   assert(EPSFRAC(lb, 0.0) == 0.0 || !SCIPvarIsIntegral(var));  /*lint !e835*/
   assert(EPSFRAC(ub, 0.0) == 0.0 || !SCIPvarIsIntegral(var));  /*lint !e835*/

   switch( conshdlrdata->varboundrelax )
   {
      case 'n' : /* no relaxation */
         break;

      case 'a' : /* relax by absolute value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         if( !SCIPisInfinity(scip, -lb) )
         {
            /* reduce lb by epsilon, or to the next integer value, which ever is larger */
            SCIP_Real bnd = floor(lb);
            lb = MAX(bnd, lb - conshdlrdata->varboundrelaxamount);
         }

         if( !SCIPisInfinity(scip, ub) )
         {
            /* increase ub by epsilon, or to the next integer value, which ever is smaller */
            SCIP_Real bnd = ceil(ub);
            ub = MIN(bnd, ub + conshdlrdata->varboundrelaxamount);
         }

         break;
      }

      case 'b' : /* relax always by absolute value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         if( !SCIPisInfinity(scip, -lb) )
            lb -= conshdlrdata->varboundrelaxamount;

         if( !SCIPisInfinity(scip, ub) )
            ub += conshdlrdata->varboundrelaxamount;

         break;
      }

      case 'r' : /* relax by relative value */
      {
         /* do not look at integer variables, they already have integral bounds, so wouldn't be relaxed */
         if( SCIPvarIsIntegral(var) )
            break;

         /* relax bounds by epsilon*max(1,|bnd|), instead of just epsilon as in case 'a', thus we trust the first log(epsilon) digits
          * however, when domains get small, relaxing can excessively weaken bound tightening, thus do only fraction of |ub-lb| if that is smaller
          * further, do not relax beyond next integer value
          */
         if( !SCIPisInfinity(scip, -lb) )
         {
            SCIP_Real bnd = floor(lb);
            lb = MAX(bnd, lb - MIN(conshdlrdata->varboundrelaxamount * MAX(1.0, REALABS(lb)), 0.001 * REALABS(ub-lb)));
         }

         if( !SCIPisInfinity(scip, ub) )
         {
            SCIP_Real bnd = ceil(ub);
            ub = MIN(bnd, ub + MIN(conshdlrdata->varboundrelaxamount * MAX(1.0, REALABS(ub)), 0.001 * REALABS(ub-lb)));
         }

         break;
      }

      default :
      {
         SCIPerrorMessage("Unsupported value '%c' for varboundrelax option.\n", conshdlrdata->varboundrelax);
         SCIPABORT();
         break;
      }
   }

   /* convert SCIPinfinity() to SCIP_INTERVAL_INFINITY */
   lb = -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb);
   ub =  infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, ub);
   assert(lb <= ub);

   SCIPintervalSetBounds(&interval, lb, ub);

   return interval;
}

/** compares two nonlinear constraints by its index
 *
 * Usable as compare operator in array sort functions.
 */
static
SCIP_DECL_SORTPTRCOMP(compIndexConsNonlinear)
{
   SCIP_CONSDATA* consdata1 = SCIPconsGetData((SCIP_CONS*)elem1);
   SCIP_CONSDATA* consdata2 = SCIPconsGetData((SCIP_CONS*)elem2);

   assert(consdata1 != NULL);
   assert(consdata2 != NULL);

   return consdata1->consindex - consdata2->consindex;
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{  /*lint --e{715}*/
   SCIP_EVENTTYPE eventtype;
   SCIP_EXPR* expr;
   SCIP_EXPR_OWNERDATA* ownerdata;

   eventtype = SCIPeventGetType(event);
   assert(eventtype & (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED));

   assert(eventdata != NULL);
   expr = (SCIP_EXPR*) eventdata;
   assert(SCIPisExprVar(scip, expr));

   SCIPdebugMsg(scip, "  exec event %" SCIP_EVENTTYPE_FORMAT " for variable <%s> (local [%g,%g], global [%g,%g])\n", eventtype,
         SCIPvarGetName(SCIPeventGetVar(event)),
         SCIPvarGetLbLocal(SCIPeventGetVar(event)), SCIPvarGetUbLocal(SCIPeventGetVar(event)),
         SCIPvarGetLbGlobal(SCIPeventGetVar(event)), SCIPvarGetUbGlobal(SCIPeventGetVar(event)));

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   /* we only catch varevents for variables in constraints, so there should be constraints */
   assert(ownerdata->nconss > 0);
   assert(ownerdata->conss != NULL);

   /* notify constraints that use this variable expression (expr) to repropagate and possibly resimplify
    * - propagation can only find something new if a bound was tightened
    * - simplify can only find something new if a var is fixed (or maybe a bound is tightened)
    *   and we look at global changes (that is, we are not looking at boundchanges in probing)
    */
   if( eventtype & (SCIP_EVENTTYPE_BOUNDTIGHTENED | SCIP_EVENTTYPE_VARFIXED) )
   {
      SCIP_CONSDATA* consdata;
      int c;

      for( c = 0; c < ownerdata->nconss; ++c )
      {
         assert(ownerdata->conss[c] != NULL);
         consdata = SCIPconsGetData(ownerdata->conss[c]);

         /* if bound tightening, then mark constraints to be propagated again
          * TODO we could try be more selective here and only trigger a propagation if a relevant bound has changed,
          *   that is, we don't need to repropagate x + ... <= rhs if only the upper bound of x has been tightened
          *   the locks don't help since they are not available separately for each constraint
          */
         if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
         {
            consdata->ispropagated = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for propagate\n", SCIPconsGetName(ownerdata->conss[c]));
         }

         /* if still in presolve (but not probing), then mark constraints to be unsimplified */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && !SCIPinProbing(scip) )
         {
            consdata->issimplified = FALSE;
            SCIPdebugMsg(scip, "  marked <%s> for simplify\n", SCIPconsGetName(ownerdata->conss[c]));
         }
      }
   }

   /* update curboundstag, lastboundrelax, and expr activity */
   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_INTERVAL activity;

      conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
      assert(conshdlrdata != NULL);

      /* increase tag on bounds */
      ++conshdlrdata->curboundstag;
      assert(conshdlrdata->curboundstag > 0);

      /* remember also if we relaxed bounds now */
      if( eventtype & SCIP_EVENTTYPE_BOUNDRELAXED )
         conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;

      /* update the activity of the var-expr here immediately
       * (we could call expr->activity = intevalvar(var, consdhlr) directly, but then the exprhdlr statistics are not updated)
       */
      SCIP_CALL( SCIPcallExprInteval(scip, expr, &activity, conshdlrdata->intevalvar, conshdlrdata) );
      /* activity = conshdlrdata->intevalvar(scip, SCIPgetVarExprVar(expr), conshdlrdata); */
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, "  var-exprhdlr::inteval = [%.20g, %.20g]\n", activity.inf, activity.sup);
#endif
      SCIPexprSetActivity(expr, activity, conshdlrdata->curboundstag);
   }

   return SCIP_OKAY;
}

/** registers event handler to catch variable events on variable
 *
 * Additionally, the given constraint is stored in the ownerdata of the variable-expression.
 * When an event occurs, all stored constraints are notified.
 */
static
SCIP_RETCODE catchVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_EXPR*            expr,               /**< variable expression */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;

   assert(eventhdlr != NULL);
   assert(expr != NULL);
   assert(SCIPisExprVar(scip, expr));
   assert(cons != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

#ifndef NDEBUG
   /* assert that constraint does not double-catch variable */
   {
      int i;
      for( i = 0; i < ownerdata->nconss; ++i )
         assert(ownerdata->conss[i] != cons);
   }
#endif

   /* append cons to ownerdata->conss */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ownerdata->conss, &ownerdata->consssize, ownerdata->nconss + 1) );
   ownerdata->conss[ownerdata->nconss++] = cons;
   /* we're not capturing the constraint here to avoid circular references */

   /* updated sorted flag */
   if( ownerdata->nconss <= 1 )
      ownerdata->consssorted = TRUE;
   else if( ownerdata->consssorted )
      ownerdata->consssorted = compIndexConsNonlinear(ownerdata->conss[ownerdata->nconss-2], ownerdata->conss[ownerdata->nconss-1]) > 0;

   /* catch variable events, if not done so yet (first constraint) */
   if( ownerdata->filterpos < 0 )
   {
      SCIP_EVENTTYPE eventtype;

      assert(ownerdata->nconss == 1);

      eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;

      SCIP_CALL( SCIPcatchVarEvent(scip, SCIPgetVarExprVar(expr), eventtype, eventhdlr, (SCIP_EVENTDATA*)expr, &ownerdata->filterpos) );
      assert(ownerdata->filterpos >= 0);
   }

   return SCIP_OKAY;
}

/** catch variable events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   int i;

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* check if we have catched variable events already */
   if( consdata->catchedevents )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->intevalvar == intEvalVarBoundTightening);

   SCIPdebugMsg(scip, "catchVarEvents for %s\n", SCIPconsGetName(cons));

   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      expr = consdata->varexprs[i];

      assert(expr != NULL);
      assert(SCIPisExprVar(scip, expr));

      SCIP_CALL( catchVarEvent(scip, eventhdlr, expr, cons) );

      /* from now on, activity of var-expr will usually be updated in processVarEvent if variable bound is changing
       * since we just registered this eventhdlr, we should make sure that the activity is also up to date now
       */
      if( SCIPexprGetActivityTag(expr) < conshdlrdata->curboundstag )
      {
         SCIP_INTERVAL activity;
         SCIP_CALL( SCIPcallExprInteval(scip, expr, &activity, intEvalVarBoundTightening, conshdlrdata) );
         /* activity = intEvalVarBoundTightening(scip, SCIPgetVarExprVar(expr), conshdlrdata); */
         SCIPexprSetActivity(expr, activity, conshdlrdata->curboundstag);
#ifdef DEBUG_PROP
         SCIPdebugMsg(scip, "var-exprhdlr::inteval for var <%s> = [%.20g, %.20g]\n", SCIPvarGetName(SCIPgetVarExprVar(expr)), activity.inf, activity.sup);
#endif
      }
   }

   consdata->catchedevents = TRUE;

   return SCIP_OKAY;
}

/** unregisters event handler to catch variable events on variable
 *
 * The given constraint is removed from the constraints array in the ownerdata of the variable-expression.
 * If this was the last constraint, then the event handler is unregistered for this variable.
 */
static
SCIP_RETCODE dropVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_EXPR*            expr,               /**< variable expression */
   SCIP_CONS*            cons                /**< expr constraint */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   int pos;

   assert(eventhdlr != NULL);
   assert(expr != NULL);
   assert(SCIPisExprVar(scip, expr));
   assert(cons != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   assert(ownerdata->nconss > 0);

   if( ownerdata->conss[ownerdata->nconss-1] == cons )
   {
      pos = ownerdata->nconss-1;
   }
   else
   {
      if( !ownerdata->consssorted )
      {
         SCIPsortPtr((void**)ownerdata->conss, compIndexConsNonlinear, ownerdata->nconss);
         ownerdata->consssorted = TRUE;
      }

      if( !SCIPsortedvecFindPtr((void**)ownerdata->conss, compIndexConsNonlinear, cons, ownerdata->nconss, &pos) )
      {
         SCIPerrorMessage("Constraint <%s> not in constraint array of expression for variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(SCIPgetVarExprVar(expr)));
         return SCIP_ERROR;
      }
      assert(pos >= 0 && pos < ownerdata->nconss);
   }
   assert(ownerdata->conss[pos] == cons);

   /* move last constraint into position of removed constraint */
   if( pos < ownerdata->nconss-1 )
   {
      ownerdata->conss[pos] = ownerdata->conss[ownerdata->nconss-1];
      ownerdata->consssorted = FALSE;
   }
   --ownerdata->nconss;

   /* drop variable events if that was the last constraint */
   if( ownerdata->nconss == 0 )
   {
      SCIP_EVENTTYPE eventtype;

      assert(ownerdata->filterpos >= 0);

      eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;

      SCIP_CALL( SCIPdropVarEvent(scip, SCIPgetVarExprVar(expr), eventtype, eventhdlr, (SCIP_EVENTDATA*)expr, ownerdata->filterpos) );
      ownerdata->filterpos = -1;
   }

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to drop bound change events */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if we have catched variable events already */
   if( !consdata->catchedevents )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   SCIPdebugMsg(scip, "dropVarEvents for %s\n", SCIPconsGetName(cons));

   for( i = consdata->nvarexprs - 1; i >= 0; --i )
   {
      assert(consdata->varexprs[i] != NULL);

      SCIP_CALL( dropVarEvent(scip, eventhdlr, consdata->varexprs[i], cons) );
   }

   consdata->catchedevents = FALSE;

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *
 * @attention Use copyexpr=FALSE only if expr is already "owned" by conshdlr, that is, if expressions were created with exprownerCreate() and ownerdata passed in the last two arguments
 */
static
SCIP_RETCODE createCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             copyexpr,           /**< whether to copy the expression or reuse the given expr (capture it) */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(conshdlr != NULL);
   assert(expr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( local && SCIPgetDepth(scip) != 0 )
   {
      SCIPerrorMessage("Locally valid nonlinear constraints are not supported, yet.\n");
      return SCIP_INVALIDCALL;
   }

   /* TODO we should allow for non-initial nonlinear constraints */
   if( !initial )
   {
      SCIPerrorMessage("Non-initial nonlinear constraints are not supported, yet.\n");
      return SCIP_INVALIDCALL;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &consdata) );

   if( copyexpr )
   {
      /* copy expression, thereby map variables expressions to already existing variables expressions in var2expr map, or augment var2expr map */
      SCIP_CALL( SCIPduplicateExpr(scip, expr, &consdata->expr, mapexprvar, conshdlr, exprownerCreate, (void*)conshdlr) );
   }
   else
   {
      consdata->expr = expr;
      SCIPcaptureExpr(consdata->expr);
   }
   consdata->lhs = lhs;
   consdata->rhs = rhs;
   consdata->consindex = conshdlrdata->lastconsindex++;
   consdata->curv = SCIP_EXPRCURV_UNKNOWN;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   return SCIP_OKAY;
}

/** returns absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z &le; f(x) and sets `violover` to TRUE.
 * If there are positive locks, then return the violation of z &ge; f(x) and sets `violunder` to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z = f(x).
 * If f could not be evaluated, then return SCIPinfinity() and set both `violover` and `violunder` to TRUE.
 *
 * @note This does not reevaluate the violation, but assumes that the expression has been evaluated
 */
static
SCIP_Real getExprAbsOrigViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_Real auxvarvalue;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   assert(ownerdata->auxvar != NULL);

   if( SCIPexprGetEvalValue(expr) == SCIP_INVALID )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = TRUE;
      return SCIPinfinity(scip);
   }

   auxvarvalue = SCIPgetSolVal(scip, sol, ownerdata->auxvar);

   if( ownerdata->nlocksneg > 0 && auxvarvalue > SCIPexprGetEvalValue(expr) )
   {
      if( violunder != NULL )
         *violunder = FALSE;
      if( violover != NULL )
         *violover = TRUE;
      return auxvarvalue - SCIPexprGetEvalValue(expr);
   }

   if( ownerdata->nlockspos > 0 && SCIPexprGetEvalValue(expr) > auxvarvalue )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = FALSE;
      return SCIPexprGetEvalValue(expr) - auxvarvalue;
   }

   if( violunder != NULL )
      *violunder = FALSE;
   if( violover != NULL )
      *violover = FALSE;
   return 0.0;
}

/** returns absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then return the violation of z &le; f(w) and sets `violover` to TRUE.
 * If there are positive locks, then return the violation of z &ge; f(w) and sets `violunder` to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z = f(w).
 * If f could not be evaluated, then return SCIPinfinity() and set both `violover` and `violunder` to TRUE.
 *
 * @note This does not reevaluate the violation, but assumes that f(w) is passed in with auxvalue.
 */
static
SCIP_Real getExprAbsAuxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_Real auxvarvalue;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   assert(ownerdata->auxvar != NULL);

   if( auxvalue == SCIP_INVALID )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = TRUE;
      return SCIPinfinity(scip);
   }

   auxvarvalue = SCIPgetSolVal(scip, sol, ownerdata->auxvar);

   if( ownerdata->nlocksneg > 0 && auxvarvalue > auxvalue )
   {
      if( violunder != NULL )
         *violunder = FALSE;
      if( violover != NULL )
         *violover = TRUE;
      return auxvarvalue - auxvalue;
   }

   if( ownerdata->nlockspos > 0 && auxvalue > auxvarvalue )
   {
      if( violunder != NULL )
         *violunder = TRUE;
      if( violover != NULL )
         *violover = FALSE;
      return auxvalue - auxvarvalue;
   }

   if( violunder != NULL )
      *violunder = FALSE;
   if( violover != NULL )
      *violover = FALSE;

   return 0.0;
}

/** computes violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPevalExpr(scip, consdata->expr, sol, soltag) );
   activity = SCIPexprGetEvalValue(consdata->expr);

   /* consider constraint as violated if it is undefined in the current point */
   if( activity == SCIP_INVALID )
   {
      consdata->lhsviol = SCIPinfinity(scip);
      consdata->rhsviol = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   /* compute violations */
   consdata->lhsviol = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : consdata->lhs  - activity;
   consdata->rhsviol = SCIPisInfinity(scip,  consdata->rhs) ? -SCIPinfinity(scip) : activity - consdata->rhs;

   return SCIP_OKAY;
}

/** returns absolute violation of a constraint
 *
 * @note This does not reevaluate the violation, but assumes that computeViolation() has been called before.
 */
static
SCIP_Real getConsAbsViolation(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return MAX3(0.0, consdata->lhsviol, consdata->rhsviol);
}

/** computes relative violation of a constraint
 *
 * @note This does not reevaluate the violation, but assumes that computeViolation() has been called before.
 */
static
SCIP_RETCODE getConsRelViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            viol,               /**< buffer to store violation */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0 */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real scale;

   assert(cons != NULL);
   assert(viol != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *viol = getConsAbsViolation(cons);

   if( conshdlrdata->violscale == 'n' )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, *viol) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( conshdlrdata->violscale == 'a' )
   {
      scale = MAX(1.0, REALABS(SCIPexprGetEvalValue(consdata->expr)));

      /* consider value of side that is violated for scaling, too */
      if( consdata->lhsviol > 0.0 && REALABS(consdata->lhs) > scale )
      {
         assert(!SCIPisInfinity(scip, -consdata->lhs));
         scale = REALABS(consdata->lhs);
      }
      else if( consdata->rhsviol > 0.0 && REALABS(consdata->rhs) > scale )
      {
         assert(!SCIPisInfinity(scip,  consdata->rhs));
         scale = REALABS(consdata->rhs);
      }

      *viol /= scale;
      return SCIP_OKAY;
   }

   /* if not 'n' or 'a', then it has to be 'g' at the moment */
   assert(conshdlrdata->violscale == 'g');
   if( soltag == 0L || consdata->gradnormsoltag != soltag )
   {
      /* we need the varexprs to conveniently access the gradient */
      SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

      /* update cached value of norm of gradient */
      consdata->gradnorm = 0.0;

      /* compute gradient */
      SCIP_CALL( SCIPevalExprGradient(scip, consdata->expr, sol, soltag) );

      /* gradient evaluation error -> no scaling */
      if( SCIPexprGetDerivative(consdata->expr) != SCIP_INVALID )
      {
         int i;
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_Real deriv;

            assert(SCIPexprGetDiffTag(consdata->expr) == SCIPexprGetDiffTag(consdata->varexprs[i]));
            deriv = SCIPexprGetDerivative(consdata->varexprs[i]);
            if( deriv == SCIP_INVALID )
            {
               /* SCIPdebugMsg(scip, "gradient evaluation error for component %d\n", i); */
               consdata->gradnorm = 0.0;
               break;
            }

            consdata->gradnorm += deriv*deriv;
         }
      }
      consdata->gradnorm = sqrt(consdata->gradnorm);
      consdata->gradnormsoltag = soltag;
   }

   *viol /= MAX(1.0, consdata->gradnorm);

   return SCIP_OKAY;
}

/** returns whether constraint is currently violated
 *
 * @note This does not reevaluate the violation, but assumes that computeViolation() has been called before.
 */
static
SCIP_Bool isConsViolated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   return getConsAbsViolation(cons) > SCIPfeastol(scip);
}

/** checks for a linear variable that can be increased or decreased without harming feasibility */
static
void findUnlockedLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int poslock;
   int neglock;
   int i;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->linvarincr = NULL;
   consdata->linvardecr = NULL;
   consdata->linvarincrcoef = 0.0;
   consdata->linvardecrcoef = 0.0;

   /* root expression is not a sum -> no unlocked linear variable available */
   if( !SCIPisExprSum(scip, consdata->expr) )
      return;

   for( i = 0; i < SCIPexprGetNChildren(consdata->expr); ++i )
   {
      SCIP_EXPR* child;

      child = SCIPexprGetChildren(consdata->expr)[i];
      assert(child != NULL);

      /* check whether the child is a variable expression */
      if( SCIPisExprVar(scip, child) )
      {
         SCIP_VAR* var = SCIPgetVarExprVar(child);
         SCIP_Real coef = SCIPgetCoefsExprSum(consdata->expr)[i];

         if( coef > 0.0 )
         {
            poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         }
         else
         {
            poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
            neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         }
         SCIPdebugMsg(scip, "child <%s> locks: %d %d\n", SCIPvarGetName(var), SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL), SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL));

         if( SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL) - neglock == 0 )
         {
            /* for a*x + f(y) \in [lhs, rhs], we can decrease x without harming other constraints
             * if we have already one candidate, then take the one where the loss in the objective function is less
             */
            if( (consdata->linvardecr == NULL) ||
               (SCIPvarGetObj(consdata->linvardecr) / consdata->linvardecrcoef > SCIPvarGetObj(var) / coef) )
            {
               consdata->linvardecr = var;
               consdata->linvardecrcoef = coef;
            }
         }

         if( SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL) - poslock == 0 )
         {
            /* for a*x + f(y) \in [lhs, rhs], we can increase x without harm
             * if we have already one candidate, then take the one where the loss in the objective function is less
             */
            if( (consdata->linvarincr == NULL) ||
               (SCIPvarGetObj(consdata->linvarincr) / consdata->linvarincrcoef > SCIPvarGetObj(var) / coef) )
            {
               consdata->linvarincr = var;
               consdata->linvarincrcoef = coef;
            }
         }
      }
   }

   assert(consdata->linvarincr == NULL || consdata->linvarincrcoef != 0.0);
   assert(consdata->linvardecr == NULL || consdata->linvardecrcoef != 0.0);

   if( consdata->linvarincr != NULL )
   {
      SCIPdebugMsg(scip, "may increase <%s> to become feasible\n", SCIPvarGetName(consdata->linvarincr));
   }
   if( consdata->linvardecr != NULL )
   {
      SCIPdebugMsg(scip, "may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->linvardecr));
   }
}

/** Given a solution where every nonlinear constraint is either feasible or can be made feasible by
 *  moving a linear variable, construct the corresponding feasible solution and pass it to the trysol heuristic.
 *
 *  The method assumes that this is always possible and that not all constraints are feasible already.
 */
static
SCIP_RETCODE proposeFeasibleSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to process */
   SCIP_Bool*            success             /**< buffer to store whether we succeeded to construct a solution that satisfies all provided constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* newsol;
   int c;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(success != NULL);

   *success = FALSE;

   /* don't propose new solutions if not in presolve or solving */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( sol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );
   SCIPdebugMsg(scip, "attempt to make solution from <%s> feasible by shifting linear variable\n",
      sol != NULL ? (SCIPsolGetHeur(sol) != NULL ? SCIPheurGetName(SCIPsolGetHeur(sol)) : "tree") : "LP");

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      SCIP_Real viol = 0.0;
      SCIP_Real delta;
      SCIP_Real gap;

      assert(consdata != NULL);

      /* get absolute violation and sign */
      if( consdata->lhsviol > SCIPfeastol(scip) )
         viol = consdata->lhsviol; /* lhs - activity */
      else if( consdata->rhsviol > SCIPfeastol(scip) )
         viol = -consdata->rhsviol; /* rhs - activity */
      else
         continue; /* constraint is satisfied */

      if( consdata->linvarincr != NULL &&
         ((viol > 0.0 && consdata->linvarincrcoef > 0.0) || (viol < 0.0 && consdata->linvarincrcoef < 0.0)) )
      {
         SCIP_VAR* var = consdata->linvarincr;

         /* compute how much we would like to increase var */
         delta = viol / consdata->linvarincrcoef;
         assert(delta > 0.0);

         /* if var has an upper bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            gap = SCIPvarGetUbGlobal(var) - SCIPgetSolVal(scip, newsol, var);
            delta = MIN(MAX(0.0, gap), delta);
         }
         if( SCIPisPositive(scip, delta) )
         {
            /* if variable is integral, round delta up so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPceil(scip, delta);

            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy lhs-violation %g of cons <%s>\n",
               SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));  /*lint !e613*/

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->linvarincrcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->linvardecr != NULL &&
         ((viol > 0.0 && consdata->linvardecrcoef < 0.0) || (viol < 0.0 && consdata->linvardecrcoef > 0.0)) )
      {
         SCIP_VAR* var = consdata->linvardecr;

         /* compute how much we would like to decrease var */
         delta = viol / consdata->linvardecrcoef;
         assert(delta < 0.0);

         /* if var has a lower bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
         {
            gap = SCIPgetSolVal(scip, newsol, var) - SCIPvarGetLbGlobal(var);
            delta = MAX(MIN(0.0, gap), delta);
         }
         if( SCIPisNegative(scip, delta) )
         {
            /* if variable is integral, round delta down so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPfloor(scip, delta);
            SCIP_CALL( SCIPincSolVal(scip, newsol, consdata->linvardecr, delta) );
            /*lint --e{613} */
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy rhs-violation %g of cons <%s>\n",
               SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->linvardecrcoef * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      /* still here... so probably we could not make constraint feasible due to variable bounds, thus give up */
      break;
   }

   /* if we have a solution that should satisfy all quadratic constraints and has a better objective than the current upper bound,
    * then pass it to the trysol heuristic
    */
   if( c == nconss && (SCIPisInfinity(scip, SCIPgetUpperbound(scip)) || SCIPisSumLT(scip, SCIPgetSolTransObj(scip, newsol), SCIPgetUpperbound(scip))) )
   {
      SCIPdebugMsg(scip, "pass solution with objective val %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      assert(conshdlrdata->trysolheur != NULL);
      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );

      *success = TRUE;
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** adds globally valid tight estimators in a given solution as cut to cutpool
 *
 * Called by addTightEstimatorCuts() for a specific expression, nlhdlr, and estimate-direction (over or under).
 */
static
SCIP_RETCODE addTightEstimatorCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_EXPR*            expr,               /**< expression */
   EXPRENFO*             exprenfo,           /**< expression enfo data, e.g., nlhdlr to use */
   SCIP_SOL*             sol,                /**< reference point where to estimate */
   SCIP_Bool             overestimate,       /**< whether to overestimate */
   SCIP_PTRARRAY*        rowpreps            /**< array for rowpreps */
   )
{
   SCIP_Bool estimatesuccess = FALSE;
   SCIP_Bool branchscoresuccess = FALSE;
   int minidx;
   int maxidx;
   int r;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprenfo != NULL);
   assert(rowpreps != NULL);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   %sestimate using nlhdlr <%s> for expr %p (%s)\n",
      overestimate ? "over" : "under", SCIPnlhdlrGetName(exprenfo->nlhdlr), (void*)expr, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr))); )

   SCIP_CALL( SCIPnlhdlrEstimate(scip, conshdlr, exprenfo->nlhdlr, expr, exprenfo->nlhdlrexprdata, sol,
      exprenfo->auxvalue, overestimate, overestimate ? SCIPinfinity(scip) : -SCIPinfinity(scip), FALSE, rowpreps, &estimatesuccess, &branchscoresuccess) );

   minidx = SCIPgetPtrarrayMinIdx(scip, rowpreps);
   maxidx = SCIPgetPtrarrayMaxIdx(scip, rowpreps);
   assert(estimatesuccess == (minidx <= maxidx));

   if( !estimatesuccess )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s failed\n", SCIPnlhdlrGetName(exprenfo->nlhdlr)); )
      return SCIP_OKAY;
   }

   for( r = minidx; r <= maxidx; ++r )
   {
      SCIP_ROWPREP* rowprep;
      SCIP_ROW* row;
      SCIP_Real estimateval;
      int i;

      rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, r);
      assert(rowprep != NULL);
      assert(SCIProwprepGetSidetype(rowprep) == (overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT));

      /* if estimators is only local valid, then skip */
      if( SCIProwprepIsLocal(rowprep) )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    skip local estimator\n"); )
         SCIPfreeRowprep(scip, &rowprep);
         continue;
      }

      /* compute value of estimator */
      estimateval = -SCIProwprepGetSide(rowprep);
      for( i = 0; i < SCIProwprepGetNVars(rowprep); ++i )
         estimateval += SCIProwprepGetCoefs(rowprep)[i] * SCIPgetSolVal(scip, sol, SCIProwprepGetVars(rowprep)[i]);

      /* if estimator value is not tight (or even "more than tight", e.g., when estimating in integer vars), then skip */
      if( (overestimate && !SCIPisFeasLE(scip, estimateval, SCIPexprGetEvalValue(expr))) ||
         (!overestimate && !SCIPisFeasGE(scip, estimateval, SCIPexprGetEvalValue(expr))) )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    skip non-tight estimator with value %g, expr value %g\n", estimateval, SCIPexprGetEvalValue(expr)); )
         SCIPfreeRowprep(scip, &rowprep);
         continue;
      }

      /* complete estimator to cut and clean it up */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetExprAuxVarNonlinear(expr), -1.0) );
      SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIPinfinity(scip), &estimatesuccess) );

      /* if cleanup failed or rowprep is local now, then skip */
      if( !estimatesuccess || SCIProwprepIsLocal(rowprep) )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    skip after cleanup failed or made estimator locally valid\n"); )
         SCIPfreeRowprep(scip, &rowprep);
         continue;
      }

      /* generate row and add to cutpool */
      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    adding cut ");
      SCIP_CALL( SCIPprintRow(scip, row, enfologfile) ); )

      SCIP_CALL( SCIPaddPoolCut(scip, row) );
      /* SCIPnlhdlrIncrementNSeparated(nlhdlr); */

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      SCIPfreeRowprep(scip, &rowprep);
   }

   SCIP_CALL( SCIPclearPtrarray(scip, rowpreps) );

   return SCIP_OKAY;
}

/** adds globally valid tight estimators in a given solution as cuts to cutpool
 *
 * Essentially we want to ensure that the LP relaxation is tight in the new solution, if possible.
 * For convex constraints, we would achieve this by linearizing.
 * To avoid checking explicitly for convexity, we compute estimators via any nlhdlr that didn't say it would
 * use bound information and check whether the estimator is tight.
 *
 * Since linearization may happen in auxiliary variables, we ensure that auxiliary variables are set
 * to the eval-value of its expression, i.e., we change sol so it is also feasible in the extended formulation.
 */
static
SCIP_RETCODE addTightEstimatorCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol                 /**< reference point where to estimate */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Longint soltag;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_PTRARRAY* rowpreps;
   int c, e;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "add tight estimators in new solution from <%s> to cutpool\n", SCIPheurGetName(SCIPsolGetHeur(sol))); )

   /* TODO probably we just evaluated all expressions when checking the sol before it was added
    * would be nice to recognize this and skip reevaluating
    */
   soltag = SCIPgetExprNewSoltag(scip);

   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   for( c = 0; c < nconss; ++c )
   {
      /* skip constraints that are not enabled or deleted or have separation disabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) || !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* TODO we could remember for which constraints there is a chance that we would add anything,
       * i.e., there is some convex-like expression, and skip other constraints
       */

      ENFOLOG(
      {
         int i;
         SCIPinfoMessage(scip, enfologfile, " constraint ");
         SCIP_CALL( SCIPprintCons(scip, conss[c], enfologfile) );
         SCIPinfoMessage(scip, enfologfile, "\n and point\n");
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_VAR* var;
            var = SCIPgetVarExprVar(consdata->varexprs[i]);
            SCIPinfoMessage(scip, enfologfile, "  %-10s = %15g bounds: [%15g,%15g]\n", SCIPvarGetName(var),
                  SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         }
      })

      SCIP_CALL( SCIPevalExpr(scip, consdata->expr, sol, soltag) );
      assert(SCIPexprGetEvalValue(consdata->expr) != SCIP_INVALID);

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         SCIP_EXPR_OWNERDATA* ownerdata;

         ownerdata = SCIPexprGetOwnerData(expr);
         assert(ownerdata != NULL);

         /* we can only generate a cut from an estimator if there is an auxvar */
         if( ownerdata->auxvar == NULL )
            continue;

         /* set value for auxvar in sol to value of expr, in case it is used to compute estimators higher up of this expression */
         assert(SCIPexprGetEvalTag(expr) == soltag);
         assert(SCIPexprGetEvalValue(expr) != SCIP_INVALID);
         SCIP_CALL( SCIPsetSolVal(scip, sol, ownerdata->auxvar, SCIPexprGetEvalValue(expr)) );

         /* generate cuts from estimators of each nonlinear handler that provides estimates */
         for( e = 0; e < ownerdata->nenfos; ++e )
         {
            SCIP_NLHDLR* nlhdlr;

            nlhdlr = ownerdata->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* skip nlhdlr that does not implement estimate (so it does enfo) */
            if( !SCIPnlhdlrHasEstimate(nlhdlr) )
               continue;

            /* skip nlhdlr that does not participate in separation or looks like it would give only locally-valid estimators
             * (because it uses activities on vars/auxvars)
             */
            if( ((ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPAABOVE) == 0 || ownerdata->enfos[e]->sepaaboveusesactivity) &&
                ((ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABELOW) == 0 || ownerdata->enfos[e]->sepabelowusesactivity) )
               continue;

            /* skip nlhdlr_default on sum, as the estimator doesn't depend on the reference point (expr is linear in auxvars) */
            if( SCIPisExprSum(scip, expr) && strcmp(SCIPnlhdlrGetName(nlhdlr), "default") == 0 )
               continue;

            /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables, since some nlhdlr expect this before their estimate is called */
            SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata, &ownerdata->enfos[e]->auxvalue, sol) );
            ENFOLOG(
               SCIPinfoMessage(scip, enfologfile, "  expr ");
               SCIPprintExpr(scip, expr, enfologfile);
               SCIPinfoMessage(scip, enfologfile, " (%p): evalvalue %.15g auxvarvalue %.15g, nlhdlr <%s> auxvalue: %.15g\n",
                  (void*)expr, SCIPexprGetEvalValue(expr), SCIPgetSolVal(scip, sol, ownerdata->auxvar), SCIPnlhdlrGetName(nlhdlr), ownerdata->enfos[e]->auxvalue);
            )
            /* due to setting values of auxvars to expr values in sol, the auxvalue should equal to expr evalvalue */
            assert(SCIPisEQ(scip, ownerdata->enfos[e]->auxvalue, SCIPexprGetEvalValue(expr)));

            /* if nlhdlr wants to be called for overestimate and does not use local bounds, then call estimate of nlhdlr */
            if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPAABOVE) && !ownerdata->enfos[e]->sepaaboveusesactivity )
            {
               SCIP_CALL( addTightEstimatorCut(scip, conshdlr, conss[c], expr, ownerdata->enfos[e], sol, TRUE, rowpreps) );
            }

            /* if nlhdlr wants to be called for underestimate and does not use local bounds, then call estimate of nlhdlr */
            if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABELOW) && !ownerdata->enfos[e]->sepabelowusesactivity )
            {
               SCIP_CALL( addTightEstimatorCut(scip, conshdlr, conss[c], expr, ownerdata->enfos[e], sol, FALSE, rowpreps) );
            }
         }
      }
   }

   SCIPfreeExpriter(&it);
   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* sol;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   if( SCIPconshdlrGetNConss(conshdlr) == 0 )
      return SCIP_OKAY;

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* we are only interested in solution coming from some heuristic other than trysol, but not from the tree
    * the reason for ignoring trysol solutions is that they may come ~~from an NLP solve in sepalp, where we already added linearizations, or are~~
    * from the tree, but postprocessed via proposeFeasibleSolution
    */
   if( SCIPsolGetHeur(sol) == NULL || SCIPsolGetHeur(sol) == conshdlrdata->trysolheur )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "caught new sol event %" SCIP_EVENTTYPE_FORMAT " from heur <%s>\n", SCIPeventGetType(event), SCIPheurGetName(SCIPsolGetHeur(sol)));

   SCIP_CALL( addTightEstimatorCuts(scip, conshdlr, SCIPconshdlrGetConss(conshdlr), SCIPconshdlrGetNConss(conshdlr), sol) );

   return SCIP_OKAY;
}

/** tightens the bounds of the auxiliary variable associated with an expression (or original variable if being a variable-expression) according to given bounds
 *
 *  The given bounds may very well be the exprs activity (when called from forwardPropExpr()), but can also be some
 *  tighter bounds (when called from SCIPtightenExprIntervalNonlinear()).
 *
 *  Nothing will happen if SCIP is not in presolve or solve.
 */
static
SCIP_RETCODE tightenAuxVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_EXPR*            expr,               /**< expression whose auxvar is to be tightened */
   SCIP_INTERVAL         bounds,             /**< bounds to be used for tightening (must not be empty) */
   SCIP_Bool*            cutoff,             /**< buffer to store whether a cutoff was detected */
   int*                  ntightenings        /**< buffer to add the total number of tightenings, or NULL */
   )
{
   SCIP_VAR* var;
   SCIP_Bool tightenedlb;
   SCIP_Bool tightenedub;
   SCIP_Bool force;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);

   /* the given bounds must not be empty (we could cope, but we shouldn't be called in this situation) */
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bounds));

   *cutoff = FALSE;

   /* do not tighten variable in problem stage (important for unittests)
    * TODO put some kind of #ifdef UNITTEST around this
    */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) > SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   var = SCIPgetExprAuxVarNonlinear(expr);
   if( var == NULL )
      return SCIP_OKAY;

   /* force tightening if conshdlrdata says so or it would mean fixing the variable */
   force = SCIPconshdlrGetData(conshdlr)->forceboundtightening || SCIPisEQ(scip, bounds.inf, bounds.sup);

   /* try to tighten lower bound of (auxiliary) variable */
   SCIP_CALL( SCIPtightenVarLb(scip, var, bounds.inf, force, cutoff, &tightenedlb) );
   if( tightenedlb )
   {
      if( ntightenings != NULL )
         ++*ntightenings;
      SCIPdebugMsg(scip, "tightened lb on auxvar <%s> to %.15g (forced:%u)\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), force);
   }
   if( *cutoff )
   {
      SCIPdebugMsg(scip, "cutoff when tightening lb on auxvar <%s> to %.15g\n", SCIPvarGetName(var), bounds.inf);
      return SCIP_OKAY;
   }

   /* try to tighten upper bound of (auxiliary) variable */
   SCIP_CALL( SCIPtightenVarUb(scip, var, bounds.sup, force, cutoff, &tightenedub) );
   if( tightenedub )
   {
      if( ntightenings != NULL )
         ++*ntightenings;
      SCIPdebugMsg(scip, "tightened ub on auxvar <%s> to %.15g (forced:%u)\n", SCIPvarGetName(var), SCIPvarGetUbLocal(var), force);
   }
   if( *cutoff )
   {
      SCIPdebugMsg(scip, "cutoff when tightening ub on auxvar <%s> to %.15g\n", SCIPvarGetName(var), bounds.sup);
      return SCIP_OKAY;
   }

   /* TODO expr->activity should have been reevaluated now due to boundchange-events, but it used to relax bounds
    * that seems unnecessary and we could easily undo this here, e.g.,
    * if( tightenedlb ) expr->activity.inf = bounds.inf
    */

   return SCIP_OKAY;
}

/** propagate bounds of the expressions in a given expression tree (that is, updates activity intervals)
 *  and tries to tighten the bounds of the auxiliary variables accordingly
 */
static
SCIP_RETCODE forwardPropExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_EXPR*            rootexpr,           /**< expression */
   SCIP_Bool             tightenauxvars,     /**< should the bounds of auxiliary variables be tightened? */
   SCIP_Bool*            infeasible,         /**< buffer to store whether the problem is infeasible (NULL if not needed) */
   int*                  ntightenings        /**< buffer to store the number of auxiliary variable tightenings (NULL if not needed) */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(rootexpr != NULL);

   if( infeasible != NULL )
      *infeasible = FALSE;
   if( ntightenings != NULL )
      *ntightenings = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if value is valid and empty, then we cannot improve, so do nothing */
   if( SCIPexprGetActivityTag(rootexpr) >= conshdlrdata->lastboundrelax && SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(rootexpr)) )
   {
      SCIPdebugMsg(scip, "stored activity of root expr is empty and valid (activitytag >= lastboundrelax (%" SCIP_LONGINT_FORMAT ")), skip forwardPropExpr -> cutoff\n", conshdlrdata->lastboundrelax);

      if( infeasible != NULL )
         *infeasible = TRUE;

      /* just update tag to curboundstag */
      SCIPexprSetActivity(rootexpr, SCIPexprGetActivity(rootexpr), conshdlrdata->curboundstag);

      return SCIP_OKAY;
   }

   /* if value is up-to-date, then nothing to do */
   if( SCIPexprGetActivityTag(rootexpr) == conshdlrdata->curboundstag )
   {
      SCIPdebugMsg(scip, "activitytag of root expr equals curboundstag (%" SCIP_LONGINT_FORMAT "), skip forwardPropExpr\n", conshdlrdata->curboundstag);

      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(rootexpr))); /* handled in previous if() */

      return SCIP_OKAY;
   }

   ownerdata = SCIPexprGetOwnerData(rootexpr);
   assert(ownerdata != NULL);

   /* if activity of rootexpr is not used, but expr participated in detect (nenfos >= 0), then we do nothing
    * it seems wrong to be called for such an expression (unless we are in detect at the moment), so I add a SCIPABORT()
    * during detect, we are in some in-between state where we may want to eval activity
    * on exprs that we did not notify about their activity usage
    */
   if( ownerdata->nenfos >= 0 && ownerdata->nactivityusesprop == 0 && ownerdata->nactivityusessepa == 0 && !conshdlrdata->indetect)
   {
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, "root expr activity is not used but enfo initialized, skip inteval\n");
#endif
      SCIPABORT();
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it);  )
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITINGCHILD :
         {
            /* skip child if it has been evaluated already */
            SCIP_EXPR* child;

            child = SCIPexpriterGetChildExprDFS(it);
            if( conshdlrdata->curboundstag == SCIPexprGetActivityTag(child) )
            {
               if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(child)) && infeasible != NULL )
                  *infeasible = TRUE;

               expr = SCIPexpriterSkipDFS(it);
               continue;
            }

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR :
         {
            SCIP_INTERVAL activity;

            /* we should not have entered this expression if its activity was already up to date */
            assert(SCIPexprGetActivityTag(expr) < conshdlrdata->curboundstag);

            ownerdata = SCIPexprGetOwnerData(expr);
            assert(ownerdata != NULL);

            /* for var exprs where varevents are catched, activity is updated immediately when the varbound has been changed
             * so we can assume that the activity is up to date for all these variables
             * UNLESS we changed the method used to evaluate activity of variable expressions
             *   or we currently use global bounds (varevents are catched for local bound changes only)
             */
            if( SCIPisExprVar(scip, expr) && ownerdata->filterpos >= 0 &&
                SCIPexprGetActivityTag(expr) >= conshdlrdata->lastvaractivitymethodchange && !conshdlrdata->globalbounds )
            {
#ifndef NDEBUG
               SCIP_INTERVAL exprhdlrinterval;

               SCIP_CALL( SCIPcallExprInteval(scip, expr, &exprhdlrinterval, conshdlrdata->intevalvar, conshdlrdata) );
               assert(SCIPisRelEQ(scip, exprhdlrinterval.inf, SCIPexprGetActivity(expr).inf));
               assert(SCIPisRelEQ(scip, exprhdlrinterval.sup, SCIPexprGetActivity(expr).sup));
#endif
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, "skip interval evaluation of expr for var <%s> [%g,%g]\n", SCIPvarGetName(SCIPgetVarExprVar(expr)), SCIPexprGetActivity(expr).inf, SCIPexprGetActivity(expr).sup);
#endif
               SCIPexprSetActivity(expr, SCIPexprGetActivity(expr), conshdlrdata->curboundstag);

               break;
            }

            if( SCIPexprGetActivityTag(expr) < conshdlrdata->lastboundrelax )
            {
               /* start with entire activity if current one is invalid */
               SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &activity);
            }
            else if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(expr)) )
            {
               /* If already empty, then don't try to compute even better activity.
                * If cons_nonlinear were alone, then we should have noted that we are infeasible
                * so an assert(infeasible == NULL || *infeasible) should work here.
                * However, after reporting a cutoff due to expr->activity being empty,
                * SCIP may wander to a different node and call propagation again.
                * If no bounds in a nonlinear constraint have been relaxed when switching nodes
                * (so expr->activitytag >= conshdlrdata->lastboundrelax), then
                * we will still have expr->activity being empty, but will have forgotten
                * that we found infeasibility here before (!2221#note_134120).
                * Therefore we just set *infeasibility=TRUE here and stop.
                */
               if( infeasible != NULL )
                  *infeasible = TRUE;
               SCIPdebugMsg(scip, "expr %p already has empty activity -> cutoff\n", (void*)expr);
               break;
            }
            else
            {
               /* start with current activity, since it is valid */
               activity = SCIPexprGetActivity(expr);
            }

            /* if activity of expr is not used, but expr participated in detect (nenfos >= 0), then do nothing */
            if( ownerdata->nenfos >= 0 && ownerdata->nactivityusesprop == 0 && ownerdata->nactivityusessepa == 0 && !conshdlrdata->indetect )
            {
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, "expr %p activity is not used but enfo initialized, skip inteval\n", (void*)expr);
#endif
               break;
            }

#ifdef DEBUG_PROP
            SCIPdebugMsg(scip, "interval evaluation of expr %p ", (void*)expr);
            SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
            SCIPdebugMsgPrint(scip, ", current activity = [%.20g, %.20g]\n", SCIPexprGetActivity(expr).inf, SCIPexprGetActivity(expr).sup);
#endif

            /* run interval eval of nonlinear handlers or expression handler */
            if( ownerdata->nenfos > 0 )
            {
               SCIP_NLHDLR* nlhdlr;
               SCIP_INTERVAL nlhdlrinterval;
               int e;

               /* for expressions with enforcement, nlhdlrs take care of interval evaluation */
               for( e = 0; e < ownerdata->nenfos && !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, activity); ++e )
               {
                  /* skip nlhdlr if it does not want to participate in activity computation */
                  if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_ACTIVITY) == 0 )
                     continue;

                  nlhdlr = ownerdata->enfos[e]->nlhdlr;
                  assert(nlhdlr != NULL);

                  /* skip nlhdlr if it does not provide interval evaluation (so it may only provide reverse propagation) */
                  if( !SCIPnlhdlrHasIntEval(nlhdlr) )
                     continue;

                  /* let nlhdlr evaluate current expression */
                  nlhdlrinterval = activity;
                  SCIP_CALL( SCIPnlhdlrInteval(scip, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata,
                     &nlhdlrinterval, conshdlrdata->intevalvar, conshdlrdata) );
#ifdef DEBUG_PROP
                  SCIPdebugMsg(scip, " nlhdlr <%s>::inteval = [%.20g, %.20g]", SCIPnlhdlrGetName(nlhdlr), nlhdlrinterval.inf, nlhdlrinterval.sup);
#endif

                  /* update activity by intersecting with computed activity */
                  SCIPintervalIntersectEps(&activity, SCIPepsilon(scip), activity, nlhdlrinterval);
#ifdef DEBUG_PROP
                  SCIPdebugMsgPrint(scip, " -> new activity: [%.20g, %.20g]\n", activity.inf, activity.sup);
#endif
               }
            }
            else
            {
               /* for node without enforcement (before or during detect), call the callback of the exprhdlr directly */
               SCIP_INTERVAL exprhdlrinterval = activity;
               SCIP_CALL( SCIPcallExprInteval(scip, expr, &exprhdlrinterval, conshdlrdata->intevalvar, conshdlrdata) );
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " exprhdlr <%s>::inteval = [%.20g, %.20g]", SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), exprhdlrinterval.inf, exprhdlrinterval.sup);
#endif

               /* update expr->activity by intersecting with computed activity */
               SCIPintervalIntersectEps(&activity, SCIPepsilon(scip), activity, exprhdlrinterval);
#ifdef DEBUG_PROP
               SCIPdebugMsgPrint(scip, " -> new activity: [%.20g, %.20g]\n", activity.inf, activity.sup);
#endif
            }

            /* if expression is integral, then we try to tighten the interval bounds a bit
             * this should undo the addition of some unnecessary safety added by use of nextafter() in interval arithmetics, e.g., when doing pow()
             * it would be ok to use ceil() and floor(), but for safety we use SCIPceil and SCIPfloor for now
             * do this only if using boundtightening-inteval and not in redundancy check (there we really want to relax all variables)
             * boundtightening-inteval does not relax integer variables, so can omit expressions without children
             * (constants should be ok, too)
             */
            if( SCIPexprIsIntegral(expr) && conshdlrdata->intevalvar == intEvalVarBoundTightening && SCIPexprGetNChildren(expr) > 0 )
            {
               if( activity.inf > -SCIP_INTERVAL_INFINITY )
                  activity.inf = SCIPceil(scip, activity.inf);
               if( activity.sup <  SCIP_INTERVAL_INFINITY )
                  activity.sup = SCIPfloor(scip, activity.sup);
#ifdef DEBUG_PROP
               SCIPdebugMsg(scip, " applying integrality: [%.20g, %.20g]\n", activity.inf, activity.sup);
#endif
            }

            /* mark the current node to be infeasible if either the lower/upper bound is above/below +/- SCIPinfinity()
             * TODO this is a problem if dual-presolve fixed a variable to +/- infinity
             */
            if( SCIPisInfinity(scip, activity.inf) || SCIPisInfinity(scip, -activity.sup) )
            {
               SCIPdebugMsg(scip, "cut off due to activity [%g,%g] beyond infinity\n", activity.inf, activity.sup);
               SCIPintervalSetEmpty(&activity);
            }

            /* now finally store activity in expr */
            SCIPexprSetActivity(expr, activity, conshdlrdata->curboundstag);

            if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, activity) )
            {
               if( infeasible != NULL )
                  *infeasible = TRUE;
            }
            else if( tightenauxvars && ownerdata->auxvar != NULL )
            {
               SCIP_Bool tighteninfeasible;

               SCIP_CALL( tightenAuxVarBounds(scip, conshdlr, expr, activity, &tighteninfeasible, ntightenings) );
               if( tighteninfeasible )
               {
                  if( infeasible != NULL )
                     *infeasible = TRUE;
                  SCIPintervalSetEmpty(&activity);
                  SCIPexprSetActivity(expr, activity, conshdlrdata->curboundstag);
               }
            }

            break;
         }

         default:
            /* you should never be here */
            SCIPerrorMessage("unexpected iterator stage\n");
            SCIPABORT();
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** returns whether intersecting `oldinterval` with `newinterval` would provide a properly smaller interval
 *
 * If `subsetsufficient` is TRUE, then the intersection being smaller than oldinterval is sufficient.
 *
 * If `subsetsufficient` is FALSE, then we require
 *  - a change from an unbounded interval to a bounded one, or
 *  - or a change from an unfixed (width > epsilon) to a fixed interval, or
 *  - a minimal tightening of one of the interval bounds as defined by SCIPis{Lb,Ub}Better().
 */
static
SCIP_Bool isIntervalBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             subsetsufficient,   /**< whether the intersection being a proper subset of oldinterval is sufficient */
   SCIP_INTERVAL         newinterval,        /**< new interval */
   SCIP_INTERVAL         oldinterval         /**< old interval */
   )
{
   assert(scip != NULL);
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newinterval));
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, oldinterval));

   if( subsetsufficient )
      /* oldinterval \cap newinterval < oldinterval iff not oldinterval is subset of newinterval */
      return !SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, oldinterval, newinterval);

   /* check whether lower bound of interval becomes finite */
   if( oldinterval.inf <= -SCIP_INTERVAL_INFINITY && newinterval.inf > -SCIP_INTERVAL_INFINITY )
      return TRUE;

   /* check whether upper bound of interval becomes finite */
   if( oldinterval.sup >=  SCIP_INTERVAL_INFINITY && newinterval.sup >  SCIP_INTERVAL_INFINITY )
      return TRUE;

   /* check whether intersection will have width <= epsilon, if oldinterval doesn't have yet */
   if( !SCIPisEQ(scip, oldinterval.inf, oldinterval.sup) && SCIPisEQ(scip, MAX(oldinterval.inf, newinterval.inf), MIN(oldinterval.sup, newinterval.sup)) )
      return TRUE;

   /* check whether lower bound on interval will be better by SCIP's quality measures for boundchanges */
   if( SCIPisLbBetter(scip, newinterval.inf, oldinterval.inf, oldinterval.sup) )
      return TRUE;

   /* check whether upper bound on interval will be better by SCIP's quality measures for boundchanges */
   if( SCIPisUbBetter(scip, newinterval.sup, oldinterval.inf, oldinterval.sup) )
      return TRUE;

   return FALSE;
}

/** propagates bounds for each sub-expression in the `reversepropqueue` by starting from the root expressions
 *
 *  The expression will be traversed in breadth first search by using this queue.
 *
 *  @note Calling this function requires feasible intervals for each sub-expression; this is guaranteed by calling
 *  forwardPropExpr() before calling this function.
 *
 *  @note Calling this function with `*infeasible` = TRUE will only empty the queue.
 */
static
SCIP_RETCODE reversePropQueue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_Bool*            infeasible,         /**< buffer to update whether an expression's bounds were propagated to an empty interval */
   int*                  ntightenings        /**< buffer to store the number of (variable) tightenings */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* expr;
   SCIP_EXPR_OWNERDATA* ownerdata;

   assert(infeasible != NULL);
   assert(ntightenings != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *ntightenings = 0;

   /* main loop that calls reverse propagation for expressions on the queue
    * when reverseprop finds a tightening for an expression, then that expression is added to the queue (within the reverseprop call)
    */
   while( !SCIPqueueIsEmpty(conshdlrdata->reversepropqueue) && !(*infeasible) )
   {
      SCIP_INTERVAL propbounds;
      int e;

      expr = (SCIP_EXPR*) SCIPqueueRemove(conshdlrdata->reversepropqueue);
      assert(expr != NULL);

      ownerdata = SCIPexprGetOwnerData(expr);
      assert(ownerdata != NULL);

      assert(ownerdata->inpropqueue);
      /* mark that the expression is not in the queue anymore */
      ownerdata->inpropqueue = FALSE;

      /* since the expr was in the propagation queue, the propbounds should belong to current propagation and should not be empty
       * (propbounds being entire doesn't make much sense, so assert this for now, too, but that could be removed)
       */
      assert(ownerdata->propboundstag == conshdlrdata->curpropboundstag);
      assert(!SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, ownerdata->propbounds));
      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, ownerdata->propbounds));

      /* this intersects propbounds with activity and auxvar bounds
       * I doubt this would be much helpful, since propbounds are already subset of activity and we also propagate
       * auxvar bounds separately, so disabling this for now
       */
#ifdef SCIP_DISABLED_CODE
      propbounds = SCIPgetExprBoundsNonlinear(scip, expr);
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, propbounds) )
      {
         *infeasible = TRUE;
         break;
      }
#else
      propbounds = ownerdata->propbounds;
#endif

      if( ownerdata->nenfos > 0 )
      {
         /* for nodes with enforcement, call reverse propagation callbacks of nlhdlrs */
         for( e = 0; e < ownerdata->nenfos && !*infeasible; ++e )
         {
            SCIP_NLHDLR* nlhdlr;
            int nreds;

            /* skip nlhdlr if it does not want to participate in activity computation */
            if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_ACTIVITY) == 0 )
               continue;

            nlhdlr = ownerdata->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* call the reverseprop of the nlhdlr */
#ifdef SCIP_DEBUG
            SCIPdebugMsg(scip, "call reverse propagation for ");
            SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
            SCIPdebugMsgPrint(scip, " in [%g,%g] using nlhdlr <%s>\n", propbounds.inf, propbounds.sup, SCIPnlhdlrGetName(nlhdlr));
#endif

            nreds = 0;
            SCIP_CALL( SCIPnlhdlrReverseprop(scip, conshdlr, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata, propbounds, infeasible, &nreds) );
            assert(nreds >= 0);
            *ntightenings += nreds;
         }
      }
      else if( SCIPexprhdlrHasReverseProp(SCIPexprGetHdlr(expr)) )
      {
         /* if expr without enforcement (before detect), call reverse propagation callback of exprhdlr directly */
         SCIP_INTERVAL* childrenbounds;
         int c;

#ifdef SCIP_DEBUG
         SCIPdebugMsg(scip, "call reverse propagation for ");
         SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
         SCIPdebugMsgPrint(scip, " in [%g,%g] using exprhdlr <%s>\n", SCIPexprGetActivity(expr).inf, SCIPexprGetActivity(expr).sup, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));
#endif

         /* if someone added an expr without nlhdlr into the reversepropqueue, then this must be because its enfo hasn't
          * been initialized in detectNlhdlr yet (nenfos < 0)
          */
         assert(ownerdata->nenfos < 0);

         SCIP_CALL( SCIPallocBufferArray(scip, &childrenbounds, SCIPexprGetNChildren(expr)) );
         for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
            childrenbounds[c] = SCIPgetExprBoundsNonlinear(scip, SCIPexprGetChildren(expr)[c]);

         /* call the reverseprop of the exprhdlr */
         SCIP_CALL( SCIPcallExprReverseprop(scip, expr, propbounds, childrenbounds, infeasible) );

         if( !*infeasible )
            for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
            {
               SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, SCIPexprGetChildren(expr)[c], childrenbounds[c], infeasible, ntightenings) );
            }

         SCIPfreeBufferArray(scip, &childrenbounds);
      }
   }

   /* reset inpropqueue for all remaining expr's in queue (can happen in case of early stop due to infeasibility) */
   while( !SCIPqueueIsEmpty(conshdlrdata->reversepropqueue) )
   {
      expr = (SCIP_EXPR*) SCIPqueueRemove(conshdlrdata->reversepropqueue);
      assert(expr != NULL);

      ownerdata = SCIPexprGetOwnerData(expr);
      assert(ownerdata != NULL);

      /* mark that the expression is not in the queue anymore */
      ownerdata->inpropqueue = FALSE;
   }

   return SCIP_OKAY;
}

/** calls domain propagation for a given set of constraints
 *
 *  The algorithm alternates calls of forward and reverse propagation.
 *  Forward propagation ensures that activity of expressions is up to date.
 *  Reverse propagation tries to derive tighter variable bounds by reversing the activity computation, using the constraints
 *  [lhs,rhs] interval as starting point.
 *
 *  The propagation algorithm works as follows:
 *   1. apply forward propagation (update activities) for all constraints not marked as propagated
 *   2. if presolve or propauxvars is disabled: collect expressions for which the constraint sides provide tighter bounds
 *      if solve and propauxvars is enabled: collect expressions for which auxvars (including those in root exprs)
 *      provide tighter bounds
 *   3. apply reverse propagation to all collected expressions; don't explore
 *      sub-expressions which have not changed since the beginning of the propagation loop
 *   4. if we have found enough tightenings go to 1, otherwise leave propagation loop
 *
 *  @note After calling forward propagation for a constraint, we mark this constraint as propagated. This flag might be
 *  reset during the reverse propagation when we find a bound tightening of a variable expression contained in the
 *  constraint. Resetting this flag is done in the EVENTEXEC callback of the event handler
 *
 *  TODO should we distinguish between expressions where activity information is used for separation and those where not,
 *    e.g., try less to propagate on convex constraints?
 */
static
SCIP_RETCODE propConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds             /**< buffer to add the number of changed bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_Bool cutoff = FALSE;
   SCIP_INTERVAL conssides;
   int ntightenings;
   int roundnr;
   SCIP_EXPRITER* revpropcollectit = NULL;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(*nchgbds >= 0);

   /* no constraints to propagate */
   if( nconss == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->intevalvar == intEvalVarBoundTightening);
   assert(!conshdlrdata->globalbounds);

   *result = SCIP_DIDNOTFIND;
   roundnr = 0;

   /* tightenAuxVarBounds() needs to know whether boundtightenings are to be forced */
   conshdlrdata->forceboundtightening = force;

   /* invalidate all propbounds (probably not needed) */
   ++conshdlrdata->curpropboundstag;

   /* create iterator that we will use if we need to look at all auxvars */
   if( conshdlrdata->propauxvars )
   {
      SCIP_CALL( SCIPcreateExpriter(scip, &revpropcollectit) );
   }

   /* main propagation loop */
   do
   {
      SCIPdebugMsg(scip, "start propagation round %d\n", roundnr);

      assert(SCIPqueueIsEmpty(conshdlrdata->reversepropqueue));

      /* apply forward propagation (update expression activities)
       * and add promising root expressions into queue for reversepropagation
       */
      for( i = 0; i < nconss; ++i )
      {
         consdata = SCIPconsGetData(conss[i]);
         assert(consdata != NULL);

         /* skip deleted, non-active, or propagation-disabled constraints */
         if( SCIPconsIsDeleted(conss[i]) || !SCIPconsIsActive(conss[i]) || !SCIPconsIsPropagationEnabled(conss[i]) )
            continue;

         /* skip already propagated constraints, i.e., constraints where no (original) variable has changed and thus
          * activity didn't change
          */
         if( consdata->ispropagated )
            continue;

         /* update activities in expression */
         SCIPdebugMsg(scip, "call forwardPropExpr() for constraint <%s> (round %d): ", SCIPconsGetName(conss[i]), roundnr);
         SCIPdebugPrintCons(scip, conss[i], NULL);

         ntightenings = 0;
         SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, TRUE, &cutoff, &ntightenings) );
         assert(cutoff || !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(consdata->expr)));

         if( cutoff )
         {
            SCIPdebugMsg(scip, " -> cutoff in forwardPropExpr (due to domain error or auxvar tightening) of constraint <%s>\n", SCIPconsGetName(conss[i]));
            *result = SCIP_CUTOFF;
            break;
         }

         ownerdata = SCIPexprGetOwnerData(consdata->expr);

         /* TODO for a constraint that only has an auxvar for consdata->expr (e.g., convex quadratic), we could also just do the if(TRUE)-branch */
         if( !conshdlrdata->propauxvars || ownerdata->auxvar == NULL )
         {
            /* check whether constraint sides (relaxed by epsilon) or auxvar bounds provide a tightening
             *   (if we have auxvar (not in presolve), then bounds of the auxvar are initially set to constraint sides,
             *   so taking auxvar bounds is enough)
             */
            if( ownerdata->auxvar == NULL )
            {
               /* relax sides by SCIPepsilon() and handle infinite sides */
               SCIP_Real lhs = SCIPisInfinity(scip, -consdata->lhs) ? -SCIP_INTERVAL_INFINITY : consdata->lhs - conshdlrdata->conssiderelaxamount;
               SCIP_Real rhs = SCIPisInfinity(scip,  consdata->rhs) ?  SCIP_INTERVAL_INFINITY : consdata->rhs + conshdlrdata->conssiderelaxamount;
               SCIPintervalSetBounds(&conssides, lhs, rhs);
            }
            else
            {
               conssides = intEvalVarBoundTightening(scip, ownerdata->auxvar, (void*)conshdlrdata);
            }
            SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, consdata->expr, conssides, &cutoff, &ntightenings) );
         }
         else
         {
            /* check whether bounds of any auxvar used in constraint provides a tightening
             *   (for the root expression, bounds of auxvar are initially set to constraint sides)
             * but skip exprs that have an auxvar, but do not participate in propagation
             */
            SCIP_EXPR* expr;

            assert(revpropcollectit != NULL);
            SCIP_CALL( SCIPexpriterInit(revpropcollectit, consdata->expr, SCIP_EXPRITER_BFS, FALSE) );
            for( expr = SCIPexpriterGetCurrent(revpropcollectit); !SCIPexpriterIsEnd(revpropcollectit) && !cutoff; expr = SCIPexpriterGetNext(revpropcollectit) )
            {
               ownerdata = SCIPexprGetOwnerData(expr);
               assert(ownerdata != NULL);

               if( ownerdata->auxvar == NULL )
                  continue;

               if( ownerdata->nactivityusesprop == 0 && ownerdata->nactivityusessepa == 0 )
                  continue;

               conssides = intEvalVarBoundTightening(scip, ownerdata->auxvar, (void*)conshdlrdata);
               SCIP_CALL( SCIPtightenExprIntervalNonlinear(scip, expr, conssides, &cutoff, &ntightenings) );
            }
         }

         if( cutoff )
         {
            SCIPdebugMsg(scip, " -> cutoff after intersect with conssides of constraint <%s>\n", SCIPconsGetName(conss[i]));
            *result = SCIP_CUTOFF;
            break;
         }

         assert(ntightenings >= 0);
         if( ntightenings > 0 )
         {
            *nchgbds += ntightenings;
            *result = SCIP_REDUCEDDOM;
         }

         /* mark constraint as propagated; this will be reset via the event system when we find a variable tightening */
         consdata->ispropagated = TRUE;
      }

      /* apply backward propagation (if cutoff is TRUE, then this call empties the queue) */
      SCIP_CALL( reversePropQueue(scip, conshdlr, &cutoff, &ntightenings) );
      assert(ntightenings >= 0);
      assert(SCIPqueueIsEmpty(conshdlrdata->reversepropqueue));

      if( cutoff )
      {
         SCIPdebugMsg(scip, " -> cutoff\n");
         *result = SCIP_CUTOFF;
         break;
      }

      if( ntightenings > 0 )
      {
         *nchgbds += ntightenings;
         *result = SCIP_REDUCEDDOM;
      }
   }
   while( ntightenings > 0 && ++roundnr < conshdlrdata->maxproprounds );

   if( conshdlrdata->propauxvars )
   {
      SCIPfreeExpriter(&revpropcollectit);
   }

   conshdlrdata->forceboundtightening = FALSE;

   /* invalidate propbounds in all exprs, so noone accidentally uses them outside propagation */
   ++conshdlrdata->curpropboundstag;

   return SCIP_OKAY;
}

/** calls the reverseprop callbacks of all nlhdlrs in all expressions in all constraints using activity as bounds
 *
 * This is meant to propagate any domain restrictions on functions onto variable bounds, if possible.
 *
 * Assumes that activities are still valid and curpropboundstag does not need to be increased.
 * Therefore, a good place to call this function is immediately after propConss() or after forwardPropExpr() if outside propagation.
 */
static
SCIP_RETCODE propExprDomains(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds             /**< buffer to add the number of changed bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_Bool cutoff = FALSE;
   int ntightenings;
   int c;
   int e;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(*nchgbds >= 0);

   assert(SCIPconshdlrGetData(conshdlr)->intevalvar == intEvalVarBoundTightening);
   assert(!SCIPconshdlrGetData(conshdlr)->globalbounds);
   assert(SCIPqueueIsEmpty(SCIPconshdlrGetData(conshdlr)->reversepropqueue));

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( c = 0; c < nconss && !cutoff; ++c )
   {
      /* skip deleted, non-active, or propagation-disabled constraints */
      if( SCIPconsIsDeleted(conss[c]) || !SCIPconsIsActive(conss[c]) || !SCIPconsIsPropagationEnabled(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it) && !cutoff; expr = SCIPexpriterGetNext(it) )
      {
         ownerdata = SCIPexprGetOwnerData(expr);
         assert(ownerdata != NULL);

         /* call reverseprop for those nlhdlr that participate in this expr's activity computation
          * this will propagate the current activity
          */
         for( e = 0; e < ownerdata->nenfos; ++e )
         {
            SCIP_NLHDLR* nlhdlr;
            assert(ownerdata->enfos[e] != NULL);

            nlhdlr = ownerdata->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);
            if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_ACTIVITY) == 0 )
               continue;

            SCIPdebugMsg(scip, "propExprDomains calling reverseprop for expression %p [%g,%g]\n", (void*)expr,
                  SCIPexprGetActivity(expr).inf, SCIPexprGetActivity(expr).sup);
            ntightenings = 0;
            SCIP_CALL( SCIPnlhdlrReverseprop(scip, conshdlr, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata,
                    SCIPexprGetActivity(expr), &cutoff, &ntightenings) );

            if( cutoff )
            {
               /* stop everything if we detected infeasibility */
               SCIPdebugMsg(scip, "detect infeasibility for constraint <%s> during reverseprop()\n", SCIPconsGetName(conss[c]));
               *result = SCIP_CUTOFF;
               break;
            }

            assert(ntightenings >= 0);
            if( ntightenings > 0 )
            {
               *nchgbds += ntightenings;
               *result = SCIP_REDUCEDDOM;
            }
         }
      }
   }

   /* apply backward propagation (if cutoff is TRUE, then this call empties the queue) */
   SCIP_CALL( reversePropQueue(scip, conshdlr, &cutoff, &ntightenings) );
   assert(ntightenings >= 0);

   if( cutoff )
   {
      SCIPdebugMsg(scip, " -> cutoff\n");
      *result = SCIP_CUTOFF;
   }
   else if( ntightenings > 0 )
   {
      *nchgbds += ntightenings;
      *result = SCIP_REDUCEDDOM;
   }

   SCIPfreeExpriter(&it);

   /* invalidate propbounds in all exprs, so noone accidentally uses them outside propagation */
   ++SCIPconshdlrGetData(conshdlr)->curpropboundstag;

   return SCIP_OKAY;
}

/** propagates variable locks through expression and adds locks to variables */
static
SCIP_RETCODE propagateLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int                   nlockspos,          /**< number of positive locks */
   int                   nlocksneg           /**< number of negative locks */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_EXPRITER* it;
   SCIP_EXPRITER_USERDATA ituserdata;

   assert(expr != NULL);

   /* if no locks, then nothing to propagate */
   if( nlockspos == 0 && nlocksneg == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_VISITINGCHILD | SCIP_EXPRITER_LEAVEEXPR);
   assert(SCIPexpriterGetCurrent(it) == expr); /* iterator should not have moved */

   /* store locks in root node */
   ituserdata.intvals[0] = nlockspos;
   ituserdata.intvals[1] = nlocksneg;
   SCIPexpriterSetCurrentUserData(it, ituserdata);

   while( !SCIPexpriterIsEnd(it) )
   {
      /* collect locks */
      ituserdata = SCIPexpriterGetCurrentUserData(it);
      nlockspos = ituserdata.intvals[0];
      nlocksneg = ituserdata.intvals[1];

      ownerdata = SCIPexprGetOwnerData(expr);

      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_ENTEREXPR:
         {
            if( SCIPisExprVar(scip, expr) )
            {
               /* if a variable, then also add nlocksneg/nlockspos via SCIPaddVarLocks() */
               SCIP_CALL( SCIPaddVarLocks(scip, SCIPgetVarExprVar(expr), nlocksneg, nlockspos) );
            }

            /* add locks to expression */
            ownerdata->nlockspos += nlockspos;
            ownerdata->nlocksneg += nlocksneg;

            /* add monotonicity information if expression has been locked for the first time */
            if( ownerdata->nlockspos == nlockspos && ownerdata->nlocksneg == nlocksneg && SCIPexprGetNChildren(expr) > 0
               && SCIPexprhdlrHasMonotonicity(SCIPexprGetHdlr(expr)) )
            {
               int i;

               assert(ownerdata->monotonicity == NULL);
               assert(ownerdata->monotonicitysize == 0);

               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ownerdata->monotonicity, SCIPexprGetNChildren(expr)) );
               ownerdata->monotonicitysize = SCIPexprGetNChildren(expr);

               /* store the monotonicity for each child */
               for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
               {
                  SCIP_CALL( SCIPcallExprMonotonicity(scip, expr, i, &ownerdata->monotonicity[i]) );
               }
            }
            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR :
         {
            /* remove monotonicity information if expression has been unlocked */
            if( ownerdata->nlockspos == 0 && ownerdata->nlocksneg == 0 && ownerdata->monotonicity != NULL )
            {
               assert(ownerdata->monotonicitysize > 0);
               /* keep this assert for checking whether someone changed an expression without updating locks properly */
               assert(ownerdata->monotonicitysize == SCIPexprGetNChildren(expr));

               SCIPfreeBlockMemoryArray(scip, &ownerdata->monotonicity, ownerdata->monotonicitysize);
               ownerdata->monotonicitysize = 0;
            }
            break;
         }

         case SCIP_EXPRITER_VISITINGCHILD :
         {
            SCIP_MONOTONE monotonicity;

            /* get monotonicity of child */
            /* NOTE: the monotonicity stored in an expression might be different from the result obtained by
             * SCIPcallExprMonotonicity
             */
            monotonicity = ownerdata->monotonicity != NULL ? ownerdata->monotonicity[SCIPexpriterGetChildIdxDFS(it)] : SCIP_MONOTONE_UNKNOWN;

            /* compute resulting locks of the child expression */
            switch( monotonicity )
            {
               case SCIP_MONOTONE_INC:
                  ituserdata.intvals[0] = nlockspos;
                  ituserdata.intvals[1] = nlocksneg;
                  break;
               case SCIP_MONOTONE_DEC:
                  ituserdata.intvals[0] = nlocksneg;
                  ituserdata.intvals[1] = nlockspos;
                  break;
               case SCIP_MONOTONE_UNKNOWN:
                  ituserdata.intvals[0] = nlockspos + nlocksneg;
                  ituserdata.intvals[1] = nlockspos + nlocksneg;
                  break;
               case SCIP_MONOTONE_CONST:
                  ituserdata.intvals[0] = 0;
                  ituserdata.intvals[1] = 0;
                  break;
            }
            /* set locks in child expression */
            SCIPexpriterSetChildUserData(it, ituserdata);

            break;
         }

         default :
            /* you should never be here */
            SCIPABORT();
            break;
      }

      expr = SCIPexpriterGetNext(it);
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** main function for adding locks to expressions and variables
 *
 * Locks for a nonlinear constraint are used to update locks for all sub-expressions and variables.
 * Locks of expressions depend on the monotonicity of expressions w.r.t. their children, e.g.,
 * consider the constraint \f$x^2 \leq 1\f$ with \f$x \in [-2,-1]\f$ implies an up-lock for the root
 * expression (pow) and a down-lock for its child \f$x\f$ because \f$x^2\f$ is decreasing on [-2,-1].
 * Since the monotonicity (and thus the locks) might also depend on variable bounds, the function remembers
 * the computed monotonicity information of each expression until all locks of an expression have been removed,
 * which implies that updating the monotonicity information during the next locking of this expression does not
 * break existing locks.
 *
 * @note When modifying the structure of an expression, e.g., during simplification, it is necessary to remove all
 *       locks from an expression and repropagating them after the structural changes have been applied.
 *       Because of existing common sub-expressions, it might be necessary to remove the locks of all constraints
 *       to ensure that an expression is unlocked (see canonicalizeConstraints() for an example)
 */
static
SCIP_RETCODE addLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   int                   nlockspos,          /**< number of positive rounding locks */
   int                   nlocksneg           /**< number of negative rounding locks */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( nlockspos == 0 && nlocksneg == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* no constraint sides -> nothing to lock */
   if( SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, -consdata->lhs) )
      return SCIP_OKAY;

   /* remember locks */
   consdata->nlockspos += nlockspos;
   consdata->nlocksneg += nlocksneg;

   assert(consdata->nlockspos >= 0);
   assert(consdata->nlocksneg >= 0);

   /* compute locks for lock propagation */
   if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlockspos + nlocksneg, nlockspos + nlocksneg));
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlockspos, nlocksneg));
   }
   else
   {
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      SCIP_CALL( propagateLocks(scip, consdata->expr, nlocksneg, nlockspos));
   }

   return SCIP_OKAY;
}

/** create a nonlinear row representation of a nonlinear constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* better curvature info will be set in initSolve() just before nlrow is added to NLP */
   SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
         0, NULL, NULL, NULL, consdata->lhs, consdata->rhs, SCIP_EXPRCURV_UNKNOWN) );

   if( SCIPisExprSum(scip, consdata->expr) )
   {
      /* if root is a sum, then split into linear and nonlinear terms */
      SCIP_EXPR* nonlinpart;
      SCIP_EXPR* child;
      SCIP_Real* coefs;
      int i;

      coefs = SCIPgetCoefsExprSum(consdata->expr);

      /* constant term of sum */
      SCIP_CALL( SCIPchgNlRowConstant(scip, consdata->nlrow, SCIPgetConstantExprSum(consdata->expr)) );

      /* a sum-expression that will hold the nonlinear terms and be passed to the nlrow eventually */
      SCIP_CALL( SCIPcreateExprSum(scip, &nonlinpart, 0, NULL, NULL, 0.0, exprownerCreate, (void*)SCIPconsGetHdlr(cons)) );

      for( i = 0; i < SCIPexprGetNChildren(consdata->expr); ++i )
      {
         child = SCIPexprGetChildren(consdata->expr)[i];
         if( SCIPisExprVar(scip, child) )
         {
            /* linear term */
            SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, SCIPgetVarExprVar(child), coefs[i]) );
         }
         else
         {
            /* nonlinear term */
            SCIP_CALL( SCIPappendExprSumExpr(scip, nonlinpart, child, coefs[i]) );
         }
      }

      if( SCIPexprGetNChildren(nonlinpart) > 0 )
      {
         /* add expression to nlrow (this will make a copy) */
         SCIP_CALL( SCIPsetNlRowExpr(scip, consdata->nlrow, nonlinpart) );
      }
      SCIP_CALL( SCIPreleaseExpr(scip, &nonlinpart) );
   }
   else
   {
      SCIP_CALL( SCIPsetNlRowExpr(scip, consdata->nlrow, consdata->expr) );
   }

   return SCIP_OKAY;
}

/** compares enfodata by enforcement priority of nonlinear handler
 *
 * If handlers have same enforcement priority, then compare by detection priority, then by name.
 */
static
SCIP_DECL_SORTPTRCOMP(enfodataCmp)
{
   SCIP_NLHDLR* h1;
   SCIP_NLHDLR* h2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   h1 = ((EXPRENFO*)elem1)->nlhdlr;
   h2 = ((EXPRENFO*)elem2)->nlhdlr;

   assert(h1 != NULL);
   assert(h2 != NULL);

   if( SCIPnlhdlrGetEnfoPriority(h1) != SCIPnlhdlrGetEnfoPriority(h2) )
      return SCIPnlhdlrGetEnfoPriority(h1) - SCIPnlhdlrGetEnfoPriority(h2);

   if( SCIPnlhdlrGetDetectPriority(h1) != SCIPnlhdlrGetDetectPriority(h2) )
      return SCIPnlhdlrGetDetectPriority(h1) - SCIPnlhdlrGetDetectPriority(h2);

   return strcmp(SCIPnlhdlrGetName(h1), SCIPnlhdlrGetName(h2));
}

/** install nlhdlrs in one expression */
static
SCIP_RETCODE detectNlhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression for which to run detection routines */
   SCIP_CONS*            cons                /**< constraint for which expr == consdata->expr, otherwise NULL */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_NLHDLR_METHOD enforcemethodsallowed;
   SCIP_NLHDLR_METHOD enforcemethods;
   SCIP_NLHDLR_METHOD enforcemethodsnew;
   SCIP_NLHDLR_METHOD nlhdlrenforcemethods;
   SCIP_NLHDLR_METHOD nlhdlrparticipating;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata;
   int enfossize;  /* allocated length of expr->enfos array */
   int h;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->auxvarid >= 0);
   assert(!conshdlrdata->indetect);

   /* there should be no enforcer yet and detection should not even have considered expr yet */
   assert(ownerdata->nenfos < 0);
   assert(ownerdata->enfos == NULL);

   /* check which enforcement methods are required by setting flags in enforcemethods for those that are NOT required
    * - if no auxiliary variable is used, then do not need sepabelow or sepaabove
    * - if auxiliary variable is used, but nobody positively (up) locks expr -> only need to enforce expr >= auxvar -> no need for underestimation
    * - if auxiliary variable is used, but nobody negatively (down) locks expr -> only need to enforce expr <= auxvar -> no need for overestimation
    * - if no one uses activity, then do not need activity methods
    */
   enforcemethods = SCIP_NLHDLR_METHOD_NONE;
   if( ownerdata->nauxvaruses == 0 )
      enforcemethods |= SCIP_NLHDLR_METHOD_SEPABOTH;
   else
   {
      if( ownerdata->nlockspos == 0 )  /* no need for underestimation */
         enforcemethods |= SCIP_NLHDLR_METHOD_SEPABELOW;
      if( ownerdata->nlocksneg == 0 )  /* no need for overestimation */
         enforcemethods |= SCIP_NLHDLR_METHOD_SEPAABOVE;
   }
   if( ownerdata->nactivityusesprop == 0 && ownerdata->nactivityusessepa == 0 )
      enforcemethods |= SCIP_NLHDLR_METHOD_ACTIVITY;

   /* it doesn't make sense to have been called on detectNlhdlr, if the expr isn't used for anything */
   assert(enforcemethods != SCIP_NLHDLR_METHOD_ALL);

   /* all methods that have not been flagged above are the ones that we want to be handled by nlhdlrs */
   enforcemethodsallowed = ~enforcemethods & SCIP_NLHDLR_METHOD_ALL;

   ownerdata->nenfos = 0;
   enfossize = 2;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ownerdata->enfos, enfossize) );
   conshdlrdata->indetect = TRUE;

   SCIPdebugMsg(scip, "detecting nlhdlrs for %s expression %p (%s); requiring%s%s%s\n",
      cons != NULL ? "root" : "non-root", (void*)expr, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)),
      (enforcemethods & SCIP_NLHDLR_METHOD_SEPABELOW) != 0 ? "" : " sepabelow",
      (enforcemethods & SCIP_NLHDLR_METHOD_SEPAABOVE) != 0 ? "" : " sepaabove",
      (enforcemethods & SCIP_NLHDLR_METHOD_ACTIVITY) != 0 ? "" : " activity");

   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
   {
      SCIP_NLHDLR* nlhdlr;

      nlhdlr = conshdlrdata->nlhdlrs[h];
      assert(nlhdlr != NULL);

      /* skip disabled nlhdlrs */
      if( !SCIPnlhdlrIsEnabled(nlhdlr) )
         continue;

      /* call detect routine of nlhdlr */
      nlhdlrexprdata = NULL;
      enforcemethodsnew = enforcemethods;
      nlhdlrparticipating = SCIP_NLHDLR_METHOD_NONE;
      conshdlrdata->registerusesactivitysepabelow = FALSE;  /* SCIPregisterExprUsageNonlinear() as called by detect may set this to TRUE */
      conshdlrdata->registerusesactivitysepaabove = FALSE;  /* SCIPregisterExprUsageNonlinear() as called by detect may set this to TRUE */
      /* coverity[forward_null] */
      SCIP_CALL( SCIPnlhdlrDetect(scip, ownerdata->conshdlr, nlhdlr, expr, cons, &enforcemethodsnew, &nlhdlrparticipating, &nlhdlrexprdata) );

      /* nlhdlr might have claimed more than needed: clean up sepa flags */
      nlhdlrparticipating &= enforcemethodsallowed;

      /* detection is only allowed to augment to nlhdlrenforcemethods, so previous enforcemethods must still be set */
      assert((enforcemethodsnew & enforcemethods) == enforcemethods);

      /* Because of the previous assert, nlhdlrenforcenew ^ enforcemethods are the methods enforced by this nlhdlr.
       * They are also cleaned up here to ensure that only the needed methods are claimed.
       */
      nlhdlrenforcemethods = (enforcemethodsnew ^ enforcemethods) & enforcemethodsallowed;

      /* nlhdlr needs to participate for the methods it is enforcing */
      assert((nlhdlrparticipating & nlhdlrenforcemethods) == nlhdlrenforcemethods);

      if( nlhdlrparticipating == SCIP_NLHDLR_METHOD_NONE )
      {
         /* nlhdlr might not have detected anything, or all set flags might have been removed by
          * clean up; in the latter case, we may need to free nlhdlrexprdata */

         /* free nlhdlr exprdata, if there is any and there is a method to free this data */
         if( nlhdlrexprdata != NULL )
         {
            SCIP_CALL( SCIPnlhdlrFreeexprdata(scip, nlhdlr, expr, &nlhdlrexprdata) );
         }
         /* nlhdlr cannot have added an enforcement method if it doesn't participate (actually redundant due to previous asserts) */
         assert(nlhdlrenforcemethods == SCIP_NLHDLR_METHOD_NONE);

         SCIPdebugMsg(scip, "nlhdlr <%s> detect unsuccessful\n", SCIPnlhdlrGetName(nlhdlr));

         continue;
      }

      SCIPdebugMsg(scip, "nlhdlr <%s> detect successful; sepabelow: %s, sepaabove: %s, activity: %s\n",
         SCIPnlhdlrGetName(nlhdlr),
         ((nlhdlrenforcemethods & SCIP_NLHDLR_METHOD_SEPABELOW) != 0) ? "enforcing" : ((nlhdlrparticipating & SCIP_NLHDLR_METHOD_SEPABELOW) != 0) ? "participating" : "no",
         ((nlhdlrenforcemethods & SCIP_NLHDLR_METHOD_SEPAABOVE) != 0) ? "enforcing" : ((nlhdlrparticipating & SCIP_NLHDLR_METHOD_SEPAABOVE) != 0) ? "participating" : "no",
         ((nlhdlrenforcemethods & SCIP_NLHDLR_METHOD_ACTIVITY) != 0) ? "enforcing" : ((nlhdlrparticipating & SCIP_NLHDLR_METHOD_ACTIVITY) != 0) ? "participating" : "no");

      /* store nlhdlr and its data */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ownerdata->enfos, &enfossize, ownerdata->nenfos+1) );
      SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[ownerdata->nenfos]) );
      ownerdata->enfos[ownerdata->nenfos]->nlhdlr = nlhdlr;
      ownerdata->enfos[ownerdata->nenfos]->nlhdlrexprdata = nlhdlrexprdata;
      ownerdata->enfos[ownerdata->nenfos]->nlhdlrparticipation = nlhdlrparticipating;
      ownerdata->enfos[ownerdata->nenfos]->issepainit = FALSE;
      ownerdata->enfos[ownerdata->nenfos]->sepabelowusesactivity = conshdlrdata->registerusesactivitysepabelow;
      ownerdata->enfos[ownerdata->nenfos]->sepaaboveusesactivity = conshdlrdata->registerusesactivitysepaabove;
      ownerdata->nenfos++;

      /* update enforcement flags */
      enforcemethods = enforcemethodsnew;
   }

   conshdlrdata->indetect = FALSE;

   /* stop if an enforcement method is missing but we are already in solving stage
    * (as long as the expression provides its callbacks, the default nlhdlr should have provided all enforcement methods)
    */
   if( enforcemethods != SCIP_NLHDLR_METHOD_ALL && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIPerrorMessage("no nonlinear handler provided some of the required enforcement methods\n");
      return SCIP_ERROR;
   }

   assert(ownerdata->nenfos > 0);

   /* sort nonlinear handlers by enforcement priority, in decreasing order */
   if( ownerdata->nenfos > 1 )
      SCIPsortDownPtr((void**)ownerdata->enfos, enfodataCmp, ownerdata->nenfos);

   /* resize enfos array to be nenfos long */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &ownerdata->enfos, enfossize, ownerdata->nenfos) );

   return SCIP_OKAY;
}

/** detect nlhdlrs that can handle the expressions */
static
SCIP_RETCODE detectNlhdlrs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints for which to run nlhdlr detect */
   int                   nconss              /**< total number of constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_EXPRITER* it;
   int i;

   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE || SCIPgetStage(scip) == SCIP_STAGE_SOLVING);  /* should only be called in presolve or initsolve or consactive */

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, TRUE) );

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPgetDepth(scip) != 0 )
   {
      /* ensure that activities are recomputed w.r.t. the global variable bounds if CONSACTIVE is called in a local node;
       * for example, this happens if globally valid nonlinear constraints are added during the tree search
       */
      SCIPincrementCurBoundsTagNonlinear(conshdlr, TRUE);
      conshdlrdata->globalbounds = TRUE;
      conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   }

   for( i = 0; i < nconss; ++i )
   {
      assert(conss != NULL && conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      /* if a constraint is separated, we currently need it to be initial, too
       * this is because INITLP will create the auxiliary variables that are used for any separation
       * TODO we may relax this with a little more programming effort when required, see also TODO in INITLP
       */
      assert((!SCIPconsIsSeparated(conss[i]) && !SCIPconsIsEnforced(conss[i])) || SCIPconsIsInitial(conss[i]));

      ownerdata = SCIPexprGetOwnerData(consdata->expr);
      assert(ownerdata != NULL);

      /* because of common sub-expressions it might happen that we already detected a nonlinear handler and added it to the expr
       * then we would normally skip to run DETECT again
       * HOWEVER: most likely we have been running DETECT with cons == NULL, which may interest less nlhdlrs
       * thus, if expr is the root expression, we rerun DETECT
       */
      if( ownerdata->nenfos > 0 )
      {
         SCIP_CALL( freeEnfoData(scip, consdata->expr, FALSE) );
         assert(ownerdata->nenfos < 0);
      }

      /* if constraint will be enforced, and we are in solve, then ensure auxiliary variable for root expression
       *   this way we can treat the root expression like any other expression when enforcing via separation
       * if constraint will be propagated, then register activity usage of root expression
       * this can trigger a call to forwardPropExpr, for which we better have the indetect flag set
       */
      conshdlrdata->indetect = TRUE;
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, consdata->expr,
         SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && (SCIPconsIsSeparated(conss[i]) || SCIPconsIsEnforced(conss[i])),
         SCIPconsIsPropagated(conss[i]),
         FALSE, FALSE) );
      conshdlrdata->indetect = FALSE;

      /* compute integrality information for all subexpressions */
      SCIP_CALL( SCIPcomputeExprIntegrality(scip, consdata->expr) );

      /* run detectNlhdlr on all expr where required */
      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         ownerdata = SCIPexprGetOwnerData(expr);
         assert(ownerdata != NULL);

         /* skip exprs that we already looked at */
         if( ownerdata->nenfos >= 0 )
            continue;

         /* if there is use of the auxvar, then someone requires that
          *   auxvar == expr (or auxvar >= expr or auxvar <= expr) or we are at the root expression (expr==consdata->expr)
          *   thus, we need to find nlhdlrs that separate or estimate
          * if there is use of the activity, then there is someone requiring that
          *   activity of this expression is updated; this someone would also benefit from better bounds on the activity of this expression
          *   thus, we need to find nlhdlrs that do interval-evaluation
          */
         if( ownerdata->nauxvaruses > 0 || ownerdata->nactivityusesprop > 0 || ownerdata->nactivityusessepa > 0 )
         {
            SCIP_CALL( detectNlhdlr(scip, expr, expr == consdata->expr ? conss[i] : NULL) );

            assert(ownerdata->nenfos >= 0);
         }
         else
         {
            /* remember that we looked at this expression during detectNlhdlrs
             * even though we have not actually run detectNlhdlr, because no nlhdlr showed interest in this expr,
             * in some situations (forwardPropExpr, to be specific) we will have to distinguish between exprs for which
             * we have not initialized enforcement yet (nenfos < 0) and expressions which are just not used in enforcement (nenfos == 0)
             */
            ownerdata->nenfos = 0;
         }
      }

      /* include this constraint into the next propagation round because the added nlhdlr may do find tighter bounds now */
      if( SCIPconsIsPropagated(conss[i]) )
         consdata->ispropagated = FALSE;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPgetDepth(scip) != 0 )
   {
      /* ensure that the local bounds are used again when reevaluating the expressions later;
       * this is only needed if CONSACTIVE is called in a local node (see begin of this function)
       */
      SCIPincrementCurBoundsTagNonlinear(conshdlr, FALSE);
      conshdlrdata->globalbounds = FALSE;
      conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   }
   else
   {
      /* ensure that all activities (except for var-exprs) are reevaluated since better methods may be available now */
      SCIPincrementCurBoundsTagNonlinear(conshdlr, FALSE);
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** initializes (pre)solving data of constraints
 *
 * This initializes data in a constraint that is used for separation, propagation, etc, and assumes that expressions will
 * not be modified.
 * In particular, this function
 * - runs the detection method of nlhldrs
 * - looks for unlocked linear variables
 * - checks curvature (if not in presolve)
 * - creates and add row to NLP (if not in presolve)
 *
 * This function can be called in presolve and solve and can be called several times with different sets of constraints,
 * e.g., it should be called in INITSOL and for constraints that are added during solve.
 */
static
SCIP_RETCODE initSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   int c;

   for( c = 0; c < nconss; ++c )
   {
      /* check for a linear variable that can be increase or decreased without harming feasibility */
      findUnlockedLinearVar(scip, conss[c]);

      if( SCIPgetStage(scip) == SCIP_STAGE_INITSOLVE || SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      {
         SCIP_CONSDATA* consdata;
         SCIP_Bool success = FALSE;

         consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
         assert(consdata != NULL);
         assert(consdata->expr != NULL);

         if( !SCIPconshdlrGetData(conshdlr)->assumeconvex )
         {
            /* call the curvature detection algorithm of the convex nonlinear handler
             * Check only for those curvature that may result in a convex inequality, i.e.,
             * whether f(x) is concave when f(x) >= lhs and/or f(x) is convex when f(x) <= rhs.
             * Also we can assume that we are nonlinear, so do not check for convex if already concave.
             */
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               SCIP_CALL( SCIPhasExprCurvature(scip, consdata->expr, SCIP_EXPRCURV_CONCAVE, &success, NULL) );
               if( success )
                  consdata->curv = SCIP_EXPRCURV_CONCAVE;
            }
            if( !success && !SCIPisInfinity(scip, consdata->rhs) )
            {
               SCIP_CALL( SCIPhasExprCurvature(scip, consdata->expr, SCIP_EXPRCURV_CONVEX, &success, NULL) );
               if( success )
                  consdata->curv = SCIP_EXPRCURV_CONVEX;
            }
         }
         else
         {
            if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
            {
               SCIPwarningMessage(scip, "Nonlinear constraint <%s> has finite left- and right-hand side, but constraints/nonlinear/assumeconvex is enabled.\n", SCIPconsGetName(conss[c]));
               consdata->curv = SCIP_EXPRCURV_LINEAR;
            }
            else
            {
               consdata->curv = !SCIPisInfinity(scip, consdata->rhs) ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_CONCAVE;
            }
         }
         SCIPdebugMsg(scip, "root curvature of constraint %s = %d\n", SCIPconsGetName(conss[c]), consdata->curv);

         /* add nlrow representation to NLP, if NLP had been constructed */
         if( SCIPisNLPConstructed(scip) && SCIPconsIsActive(conss[c]) )
         {
            if( consdata->nlrow == NULL )
            {
               SCIP_CALL( createNlRow(scip, conss[c]) );
               assert(consdata->nlrow != NULL);
            }
            SCIPnlrowSetCurvature(consdata->nlrow, consdata->curv);
            SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
         }
      }
   }

   /* register non linear handlers */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}

/** deinitializes (pre)solving data of constraints
 *
 * This removes the initialization data created in initSolve().
 *
 * This function can be called in presolve and solve.
 *
 * TODO At the moment, it should not be called for a constraint if there are other constraints
 * that use the same expressions but still require their nlhdlr.
 * We should probably only decrement the auxvar and activity usage for the root expr and then
 * proceed as in detectNlhdlrs(), i.e., free enfo data only where none is used.
 */
static
SCIP_RETCODE deinitSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_CONSDATA* consdata;
   SCIP_Bool rootactivityvalid;
   int c;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   /* call deinitialization callbacks of expression and nonlinear handlers
    * free nonlinear handlers information from expressions
    * remove auxiliary variables and nactivityuses counts from expressions
    */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      /* check and remember whether activity in root is valid */
      rootactivityvalid = SCIPexprGetActivityTag(consdata->expr) >= SCIPconshdlrGetData(conshdlr)->lastboundrelax;

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         SCIPdebugMsg(scip, "exitsepa and free nonlinear handler data for expression %p\n", (void*)expr);

         /* remove nonlinear handlers in expression and their data and auxiliary variables; reset activityusage count */
         SCIP_CALL( freeEnfoData(scip, expr, TRUE) );

         /* remove quadratic info */
         SCIPfreeExprQuadratic(scip, expr);

         if( rootactivityvalid )
         {
            /* ensure activity is valid if consdata->expr activity is valid
             * this is mainly to ensure that we do not leave invalid activities in parts of the expression tree where activity was not used,
             * e.g., an expr's activity was kept up to date by a nlhdlr, but without using some childs activity
             * so this childs activity would be invalid, which can generate confusion
             */
            SCIP_CALL( SCIPevalExprActivity(scip, expr) );
         }
      }

      if( consdata->nlrow != NULL )
      {
         /* remove row from NLP, if still in solving
          * if we are in exitsolve, the whole NLP will be freed anyway
          */
         if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
         {
            SCIP_CALL( SCIPdelNlRow(scip, consdata->nlrow) );
         }

         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }

      /* forget about linear variables that can be increased or decreased without harming feasibility */
      consdata->linvardecr = NULL;
      consdata->linvarincr = NULL;

      /* forget about curvature */
      consdata->curv = SCIP_EXPRCURV_UNKNOWN;
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** helper method to decide whether a given expression is product of at least two binary variables */
static
SCIP_Bool isBinaryProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   int i;

   assert(expr != NULL);

   /* check whether the expression is a product */
   if( !SCIPisExprProduct(scip, expr) )
      return FALSE;

   /* don't consider products with a coefficient != 1 and products with a single child
    * simplification will take care of this expression later
    */
   if( SCIPexprGetNChildren(expr) <= 1 || SCIPgetCoefExprProduct(expr) != 1.0 )
      return FALSE;

   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child;
      SCIP_VAR* var;
      SCIP_Real ub;
      SCIP_Real lb;

      child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      if( !SCIPisExprVar(scip, child) )
         return FALSE;

      var = SCIPgetVarExprVar(child);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      /* check whether variable is integer and has [0,1] as variable bounds */
      if( !SCIPvarIsIntegral(var) || !SCIPisEQ(scip, lb, 0.0) || !SCIPisEQ(scip, ub, 1.0) )
         return FALSE;
   }

   return TRUE;
}

/** helper method to collect all bilinear binary product terms */
static
SCIP_RETCODE getBilinearBinaryTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            sumexpr,            /**< sum expression */
   SCIP_VAR**            xs,                 /**< array to collect first variable of each bilinear binary product */
   SCIP_VAR**            ys,                 /**< array to collect second variable of each bilinear binary product */
   int*                  childidxs,          /**< array to store the index of the child of each stored bilinear binary product */
   int*                  nterms              /**< pointer to store the total number of bilinear binary terms */
   )
{
   int i;

   assert(sumexpr != NULL);
   assert(SCIPisExprSum(scip, sumexpr));
   assert(xs != NULL);
   assert(ys != NULL);
   assert(childidxs != NULL);
   assert(nterms != NULL);

   *nterms = 0;

   for( i = 0; i < SCIPexprGetNChildren(sumexpr); ++i )
   {
      SCIP_EXPR* child;

      child = SCIPexprGetChildren(sumexpr)[i];
      assert(child != NULL);

      if( SCIPexprGetNChildren(child) == 2 && isBinaryProduct(scip, child) )
      {
         SCIP_VAR* x = SCIPgetVarExprVar(SCIPexprGetChildren(child)[0]);
         SCIP_VAR* y = SCIPgetVarExprVar(SCIPexprGetChildren(child)[1]);

         assert(x != NULL);
         assert(y != NULL);

         if( x != y )
         {
            xs[*nterms] = x;
            ys[*nterms] = y;
            childidxs[*nterms] = i;
            ++(*nterms);
         }
      }
   }

   return SCIP_OKAY;
}

/** helper method to reformulate \f$x_i \sum_j c_{ij} x_j\f$ */
static
SCIP_RETCODE reformulateFactorizedBinaryQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             facvar,             /**< variable that has been factorized */
   SCIP_VAR**            vars,               /**< variables of sum_j c_ij x_j */
   SCIP_Real*            coefs,              /**< coefficients of sum_j c_ij x_j */
   int                   nvars,              /**< total number of variables in sum_j c_ij x_j */
   SCIP_EXPR**           newexpr,            /**< pointer to store the new expression */
   int*                  naddconss           /**< pointer to update the total number of added constraints (might be NULL) */
   )
{
   SCIP_VAR* auxvar;
   SCIP_CONS* newcons;
   SCIP_Real minact = 0.0;
   SCIP_Real maxact = 0.0;
   SCIP_Bool integral = TRUE;
   char name [SCIP_MAXSTRLEN];
   int i;

   assert(facvar != NULL);
   assert(vars != NULL);
   assert(nvars > 1);
   assert(newexpr != NULL);

   /* compute minimum and maximum activity of sum_j c_ij x_j */
   /* TODO could compute minact and maxact for facvar=0 and facvar=1 separately, taking implied bounds into account, allowing for possibly tighter big-M's below */
   for( i = 0; i < nvars; ++i )
   {
      minact += MIN(coefs[i], 0.0);
      maxact += MAX(coefs[i], 0.0);
      integral = integral && SCIPisIntegral(scip, coefs[i]);
   }
   assert(minact <= maxact);

   /* create and add auxiliary variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateVarBasic(scip, &auxvar, name, minact, maxact, 0.0, integral ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, auxvar) );

   /* create and add z - maxact x <= 0 */
   if( !SCIPisZero(scip, maxact) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_1", SCIPconsGetName(cons), SCIPvarGetName(facvar));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &newcons, name, auxvar, facvar, -maxact, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      if( naddconss != NULL )
         ++(*naddconss);
   }

   /* create and add  0 <= z - minact x */
   if( !SCIPisZero(scip, minact) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_2", SCIPconsGetName(cons), SCIPvarGetName(facvar));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &newcons, name, auxvar, facvar, -minact, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
      if( naddconss != NULL )
         ++(*naddconss);
   }

   /* create and add minact <= sum_j c_j x_j - z + minact x_i */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_3", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &newcons, name, nvars, vars, coefs, minact, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCoefLinear(scip, newcons, auxvar, -1.0) );
   if( !SCIPisZero(scip, minact) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, facvar, minact) );
   }
   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   if( naddconss != NULL )
      ++(*naddconss);

   /* create and add sum_j c_j x_j - z + maxact x_i <= maxact */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_4", SCIPconsGetName(cons), SCIPvarGetName(facvar));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &newcons, name, nvars, vars, coefs, -SCIPinfinity(scip), maxact) );
   SCIP_CALL( SCIPaddCoefLinear(scip, newcons, auxvar, -1.0) );
   if( !SCIPisZero(scip, maxact) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, facvar, maxact) );
   }
   SCIP_CALL( SCIPaddCons(scip, newcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   if( naddconss != NULL )
      ++(*naddconss);

   /* create variable expression */
   SCIP_CALL( createExprVar(scip, conshdlr, newexpr, auxvar) );

   /* release auxvar */
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

   return SCIP_OKAY;
}

/** helper method to generate an expression for a sum of products of binary variables; note that the method captures the generated expression */
static
SCIP_RETCODE getFactorizedBinaryQuadraticExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_EXPR*            sumexpr,            /**< expression */
   int                   minterms,           /**< minimum number of terms in a the sum of x_i sum_j c_j x_j */
   SCIP_EXPR**           newexpr,            /**< pointer to store the expression that represents the binary quadratic */
   int*                  naddconss           /**< pointer to update the total number of added constraints (might be NULL) */
   )
{
   SCIP_EXPR** exprs = NULL;
   SCIP_VAR** tmpvars = NULL;
   SCIP_VAR** vars = NULL;
   SCIP_VAR** xs = NULL;
   SCIP_VAR** ys = NULL;
   SCIP_Real* exprcoefs = NULL;
   SCIP_Real* tmpcoefs = NULL;
   SCIP_Real* sumcoefs;
   SCIP_Bool* isused  = NULL;
   int* childidxs = NULL;
   int* count = NULL;
   int nchildren;
   int nexprs = 0;
   int nterms;
   int nvars;
   int ntotalvars;
   int i;

   assert(sumexpr != NULL);
   assert(minterms > 1);
   assert(newexpr != NULL);

   *newexpr = NULL;

   /* check whether sumexpr is indeed a sum */
   if( !SCIPisExprSum(scip, sumexpr) )
      return SCIP_OKAY;

   nchildren = SCIPexprGetNChildren(sumexpr);
   sumcoefs = SCIPgetCoefsExprSum(sumexpr);
   nvars = SCIPgetNVars(scip);
   ntotalvars = SCIPgetNTotalVars(scip);

   /* check whether there are enough terms available */
   if( nchildren < minterms )
      return SCIP_OKAY;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &xs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ys, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &childidxs, nchildren) );

   /* collect all bilinear binary product terms */
   SCIP_CALL( getBilinearBinaryTerms(scip, sumexpr, xs, ys, childidxs, &nterms) );

   /* check whether there are enough terms available */
   if( nterms < minterms )
      goto TERMINATE;

   /* store how often each variable appears in a bilinear binary product */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, SCIPgetVars(scip), nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &count, ntotalvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &isused, nchildren) );

   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exprcoefs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, MIN(nterms, nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcoefs, MIN(nterms, nvars)) );

   for( i = 0; i < nterms; ++i )
   {
      int xidx;
      int yidx;

      assert(xs[i] != NULL);
      assert(ys[i] != NULL);

      xidx = SCIPvarGetIndex(xs[i]);
      assert(xidx < ntotalvars);
      yidx = SCIPvarGetIndex(ys[i]);
      assert(yidx < ntotalvars);

      ++count[xidx];
      ++count[yidx];

      SCIPdebugMsg(scip, "increase counter for %s to %d\n", SCIPvarGetName(xs[i]), count[xidx]);
      SCIPdebugMsg(scip, "increase counter for %s to %d\n", SCIPvarGetName(ys[i]), count[yidx]);
   }

   /* sort variables; don't change order of count array because it depends on problem indices */
   {
      int* tmpcount;

      SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpcount, count, nvars) );
      SCIPsortDownIntPtr(tmpcount, (void**)vars, nvars);
      SCIPfreeBufferArray(scip, &tmpcount);
   }

   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* facvar = vars[i];
      int ntmpvars = 0;
      int j;

      /* skip candidate if there are not enough terms left */
      if( count[SCIPvarGetIndex(vars[i])] < minterms )
         continue;

      SCIPdebugMsg(scip, "consider facvar = %s with count = %d\n", SCIPvarGetName(facvar), count[SCIPvarGetIndex(vars[i])]);

      /* collect variables for x_i * sum_j c_ij x_j */
      for( j = 0; j < nterms; ++j )
      {
         int childidx = childidxs[j];
         assert(childidx >= 0 && childidx < nchildren);

         if( !isused[childidx] && (xs[j] == facvar || ys[j] == facvar) )
         {
            SCIP_Real coef;
            int xidx;
            int yidx;

            coef = sumcoefs[childidx];
            assert(coef != 0.0);

            /* collect corresponding variable */
            tmpvars[ntmpvars] = (xs[j] == facvar) ? ys[j] : xs[j];
            tmpcoefs[ntmpvars] = coef;
            ++ntmpvars;

            /* update counters */
            xidx = SCIPvarGetIndex(xs[j]);
            assert(xidx < ntotalvars);
            yidx = SCIPvarGetIndex(ys[j]);
            assert(yidx < ntotalvars);
            --count[xidx];
            --count[yidx];
            assert(count[xidx] >= 0);
            assert(count[yidx] >= 0);

            /* mark term to be used */
            isused[childidx] = TRUE;
         }
      }
      assert(ntmpvars >= minterms);
      assert(SCIPvarGetIndex(facvar) < ntotalvars);
      assert(count[SCIPvarGetIndex(facvar)] == 0); /* facvar should not appear in any other bilinear term */

      /* create required constraints and store the generated expression */
      SCIP_CALL( reformulateFactorizedBinaryQuadratic(scip, conshdlr, cons, facvar, tmpvars, tmpcoefs, ntmpvars, &exprs[nexprs], naddconss) );
      exprcoefs[nexprs] = 1.0;
      ++nexprs;
   }

   /* factorization was only successful if at least one expression has been generated */
   if( nexprs > 0 )
   {
      int nexprsold = nexprs;

      /* add all children of the sum that have not been used */
      for( i = 0; i < nchildren; ++i )
      {
         if( !isused[i] )
         {
            exprs[nexprs] = SCIPexprGetChildren(sumexpr)[i];
            exprcoefs[nexprs] = sumcoefs[i];
            ++nexprs;
         }
      }

      /* create a new sum expression */
      SCIP_CALL( SCIPcreateExprSum(scip, newexpr, nexprs, exprs, exprcoefs, SCIPgetConstantExprSum(sumexpr), exprownerCreate, (void*)conshdlr) );

      /* release all expressions that have been generated by reformulateFactorizedBinaryQuadratic() */
      for( i = 0; i < nexprsold; ++i )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &exprs[i]) );
      }
   }

TERMINATE:
   /* free memory */
   SCIPfreeBufferArrayNull(scip, &tmpcoefs);
   SCIPfreeBufferArrayNull(scip, &tmpvars);
   SCIPfreeBufferArrayNull(scip, &exprcoefs);
   SCIPfreeBufferArrayNull(scip, &exprs);
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &isused);
   SCIPfreeBufferArrayNull(scip, &count);
   SCIPfreeBufferArray(scip, &childidxs);
   SCIPfreeBufferArray(scip, &ys);
   SCIPfreeBufferArray(scip, &xs);

   return SCIP_OKAY;
}

/** helper method to create an AND constraint or varbound constraints for a given binary product expression */
static
SCIP_RETCODE getBinaryProductExprDo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_EXPR*            prodexpr,           /**< product expression */
   SCIP_EXPR**           newexpr,            /**< pointer to store the expression that represents the product */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   SCIP_Bool             empathy4and         /**< whether to use an AND constraint, if possible */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS* cons;
   SCIP_Real* coefs;
   SCIP_VAR* w;
   char* name;
   int nchildren;
   int i;

   assert(conshdlr != NULL);
   assert(prodexpr != NULL);
   assert(SCIPisExprProduct(scip, prodexpr));
   assert(newexpr != NULL);

   nchildren = SCIPexprGetNChildren(prodexpr);
   assert(nchildren >= 2);

   /* memory to store the variables of the variable expressions (+1 for w) and their name */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nchildren + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nchildren + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &name, nchildren * (SCIP_MAXSTRLEN + 1) + 20) );

   /* prepare the names of the variable and the constraints */
   strcpy(name, "binreform");
   for( i = 0; i < nchildren; ++i )
   {
      vars[i] = SCIPgetVarExprVar(SCIPexprGetChildren(prodexpr)[i]);
      coefs[i] = 1.0;
      assert(vars[i] != NULL);
      (void) strcat(name, "_");
      (void) strcat(name, SCIPvarGetName(vars[i]));
   }

   /* create and add variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIPdebugMsg(scip, "  created auxiliary variable %s\n", name);

   /* use variable bound constraints if it is a bilinear product and there is no empathy for an AND constraint */
   if( nchildren == 2 && !empathy4and )
   {
      SCIP_VAR* x = vars[0];
      SCIP_VAR* y = vars[1];

      assert(x != NULL);
      assert(y != NULL);
      assert(x != y);

      /* create and add x - w >= 0 */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_1", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, name, x, w, -1.0, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* create and add y - w >= 0 */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_2", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, name, y, w, -1.0, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* create and add x + y - w <= 1 */
      vars[2] = w;
      coefs[2] = -1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "binreform_%s_%s_3", SCIPvarGetName(x), SCIPvarGetName(y));
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 3, vars, coefs, -SCIPinfinity(scip), 1.0) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* update number of added constraints */
      if( naddconss != NULL )
         *naddconss += 3;
   }
   else
   {
      /* create, add, and release AND constraint */
      SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, name, w, nchildren, vars) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIPdebugMsg(scip, "  create AND constraint\n");

      /* update number of added constraints */
      if( naddconss != NULL )
         *naddconss += 1;
   }

   /* create variable expression */
   SCIP_CALL( createExprVar(scip, conshdlr, newexpr, w) );

   /* release created variable */
   SCIP_CALL( SCIPreleaseVar(scip, &w) );

   /* free memory */
   SCIPfreeBufferArray(scip, &name);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** helper method to generate an expression for the product of binary variables; note that the method captures the generated expression */
static
SCIP_RETCODE getBinaryProductExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_HASHMAP*         exprmap,            /**< map to remember generated variables for visited product expressions */
   SCIP_EXPR*            prodexpr,           /**< product expression */
   SCIP_EXPR**           newexpr,            /**< pointer to store the expression that represents the product */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to update the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nchildren;

   assert(prodexpr != NULL);
   assert(newexpr != NULL);

   *newexpr = NULL;

   /* only consider products of binary variables */
   if( !isBinaryProduct(scip, prodexpr) )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   nchildren = SCIPexprGetNChildren(prodexpr);
   assert(nchildren >= 2);

   /* check whether there is already an expression that represents the product */
   if( SCIPhashmapExists(exprmap, (void*)prodexpr) )
   {
      *newexpr = (SCIP_EXPR*) SCIPhashmapGetImage(exprmap, (void*)prodexpr);
      assert(*newexpr != NULL);

      /* capture expression */
      SCIPcaptureExpr(*newexpr);
   }
   else
   {
      SCIPdebugMsg(scip, "  product expression %p has been considered for the first time\n", (void*)prodexpr);

      if( nchildren == 2 )
      {
         SCIP_CLIQUE** xcliques;
         SCIP_VAR* x;
         SCIP_VAR* y;
         SCIP_Bool found_clique = FALSE;
         int c;

         /* get variables from the product expression */
         x = SCIPgetVarExprVar(SCIPexprGetChildren(prodexpr)[0]);
         assert(x != NULL);
         y = SCIPgetVarExprVar(SCIPexprGetChildren(prodexpr)[1]);
         assert(y != NULL);
         assert(x != y);

         /* first try to find a clique containing both variables */
         xcliques = SCIPvarGetCliques(x, TRUE);

         /* look in cliques containing x */
         for( c = 0; c < SCIPvarGetNCliques(x, TRUE); ++c )
         {
            if( SCIPcliqueHasVar(xcliques[c], y, TRUE) ) /* x + y <= 1 => x*y = 0 */
            {
               /* create zero value expression */
               SCIP_CALL( SCIPcreateExprValue(scip, newexpr, 0.0, exprownerCreate, (void*)conshdlr) );

               if( nchgcoefs != NULL )
                  *nchgcoefs += 1;

               found_clique = TRUE;
               break;
            }

            if( SCIPcliqueHasVar(xcliques[c], y, FALSE) ) /* x + (1-y) <= 1 => x*y = x */
            {
               /* create variable expression for x */
               SCIP_CALL( createExprVar(scip, conshdlr, newexpr, x) );

               if( nchgcoefs != NULL )
                  *nchgcoefs += 2;

               found_clique = TRUE;
               break;
            }
         }

         if( !found_clique )
         {
            xcliques = SCIPvarGetCliques(x, FALSE);

            /* look in cliques containing complement of x */
            for( c = 0; c < SCIPvarGetNCliques(x, FALSE); ++c )
            {
               if( SCIPcliqueHasVar(xcliques[c], y, TRUE) ) /* (1-x) + y <= 1 => x*y = y */
               {
                  /* create variable expression for y */
                  SCIP_CALL( createExprVar(scip, conshdlr, newexpr, y) );

                  if( nchgcoefs != NULL )
                     *nchgcoefs += 1;

                  found_clique = TRUE;
                  break;
               }

               if( SCIPcliqueHasVar(xcliques[c], y, FALSE) ) /* (1-x) + (1-y) <= 1 => x*y = x + y - 1 */
               {
                  /* create sum expression */
                  SCIP_EXPR* sum_children[2];
                  SCIP_Real sum_coefs[2];
                  SCIP_CALL( createExprVar(scip, conshdlr, &sum_children[0], x) );
                  SCIP_CALL( createExprVar(scip, conshdlr, &sum_children[1], y) );
                  sum_coefs[0] = 1.0;
                  sum_coefs[1] = 1.0;
                  SCIP_CALL( SCIPcreateExprSum(scip, newexpr, 2, sum_children, sum_coefs, -1.0, exprownerCreate, (void*)conshdlr) );

                  SCIP_CALL( SCIPreleaseExpr(scip, &sum_children[0]) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &sum_children[1]) );

                  if( nchgcoefs != NULL )
                     *nchgcoefs += 3;

                  found_clique = TRUE;
                  break;
               }
            }
         }

         /* if the variables are not in a clique, do standard linearization */
         if( !found_clique )
         {
            SCIP_CALL( getBinaryProductExprDo(scip, conshdlr, prodexpr, newexpr, naddconss, conshdlrdata->reformbinprodsand) );
         }
      }
      else
      {
         /* linearize binary product using an AND constraint because nchildren > 2 */
         SCIP_CALL( getBinaryProductExprDo(scip, conshdlr, prodexpr, newexpr, naddconss, conshdlrdata->reformbinprodsand) );
      }

      /* hash variable expression */
      SCIP_CALL( SCIPhashmapInsert(exprmap, (void*)prodexpr, *newexpr) );
   }

   return SCIP_OKAY;
}

/** helper function to replace binary products in a given constraint */
static
SCIP_RETCODE replaceBinaryProducts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_HASHMAP*         exprmap,            /**< map to remember generated variables for visited product expressions */
   SCIP_EXPRITER*        it,                 /**< expression iterator */
   int*                  naddconss,          /**< pointer to update the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to update the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(exprmap != NULL);
   assert(it != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   SCIPdebugMsg(scip, "  check constraint %s\n", SCIPconsGetName(cons));

   for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      SCIP_EXPR* newexpr = NULL;
      SCIP_EXPR* childexpr;
      int childexpridx;

      childexpridx = SCIPexpriterGetChildIdxDFS(it);
      assert(childexpridx >= 0 && childexpridx < SCIPexprGetNChildren(expr));
      childexpr = SCIPexpriterGetChildExprDFS(it);
      assert(childexpr != NULL);

      /* try to factorize variables in a sum expression that contains several products of binary variables */
      if( conshdlrdata->reformbinprodsfac > 1 )
      {
         SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, cons, childexpr, conshdlrdata->reformbinprodsfac, &newexpr, naddconss) );
      }

      /* try to create an expression that represents a product of binary variables */
      if( newexpr == NULL )
      {
         SCIP_CALL( getBinaryProductExpr(scip, conshdlr, exprmap, childexpr, &newexpr, naddconss, nchgcoefs) );
      }

      if( newexpr != NULL )
      {
         assert(naddconss == NULL || *naddconss > 0 || nchgcoefs == NULL || *nchgcoefs > 0);

         /* replace product expression */
         SCIP_CALL( SCIPreplaceExprChild(scip, expr, childexpridx, newexpr) );

         /* note that the expression has been captured by getBinaryProductExpr and SCIPreplaceExprChild */
         SCIP_CALL( SCIPreleaseExpr(scip, &newexpr) );

         /* mark the constraint to not be simplified anymore */
         consdata->issimplified = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** reformulates products of binary variables during presolving in the following way:
 *
 * Let \f$\sum_{i,j} Q_{ij} x_i x_j\f$ be a subexpression that only contains binary variables.
 * Each term \f$x_i x_j\f$ is reformulated with the help of an extra (implicit integer) variable \f$z_{ij}\f$ in {0,1}:
 * \f[
 *    z_{ij} \leq x_i, \qquad z_{ij} \leq x_j, \qquad x_i + x_j - z_{ij} \leq 1.
 * \f]
 *
 * Before reformulating \f$x_i x_j\f$ in this way, it is checked whether there is a clique that contains \f$x_i\f$ and \f$x_j\f$.
 * These cliques allow for a better reformulation. There are four cases:
 *
 *    1. \f$x_i + x_j \leq 1\f$ implies that \f$x_i x_j = 0\f$
 *    2. \f$x_i + (1 - x_j) \leq 1\f$ implies \f$x_i x_j = x_i\f$
 *    3. \f$(1 - x_i) + x_j \leq 1\f$ implies \f$x_i x_j = x_j\f$
 *    4. \f$(1 - x_i) + (1 - x_j) \leq 1\f$ implies \f$x_i x_j = x_i + x_j - 1\f$
 *
 * The reformulation using \f$z_{ij}\f$ or the cliques is implemented in getBinaryProductExpr().
 *
 * Introducing too many extra variables and constraints can have a negative impact on the performance (e.g., due to
 * slow probing). For this reason, it is checked in getFactorizedBinaryQuadraticExpr() whether \f$\sum_{i,j} Q_{ij} x_i x_j\f$
 * contains large (&ge; `reformbinprodsfac` parameter) lower sums of the form \f$x_i \sum_j Q_{ij} x_j\f$.
 * Such a lower sum is reformulated with only one extra variable w_i:
 * \f{align}{
 *    \text{maxact} & := \sum_j \max(0, Q_{ij}), \\
 *    \text{minact} & := \sum_j \min(0, Q_{ij}), \\
 *    \text{minact}\, x_i & \leq w_i, \\
 *    w_i &\leq \text{maxact}\, x_i, \\
 *    \text{minact} &\leq \sum_j Q_{ij} x_j - w_i + \text{minact}\, x_i \\
 *    \text{maxact} &\geq \sum_j Q_{ij} x_j - w_i + \text{maxact}\, x_i
 * \f}
 * We mark \f$w_i\f$ to be implicit integer if all \f$Q_{ij}\f$ are integer. After each replacement of a lower sum, it
 * is checked whether there are enough terms left to factorize other binary variables. Lower sums with a larger number
 * of terms are prioritized.
 */
static
SCIP_RETCODE presolveBinaryProducts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< total number of constraints */
   int*                  naddconss,          /**< pointer to store the total number of added constraints (might be NULL) */
   int*                  nchgcoefs           /**< pointer to store the total number of changed coefficients (might be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_HASHMAP* exprmap;
   SCIP_EXPRITER* it;
   int c;

   assert(conshdlr != NULL);

   /* no nonlinear constraints or binary variables -> skip */
   if( nconss == 0 || SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;
   assert(conss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create expression hash map */
   SCIP_CALL( SCIPhashmapCreate(&exprmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   /* create expression iterator */
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);

   SCIPdebugMsg(scip, "call presolveBinaryProducts()\n");

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_EXPR* newexpr = NULL;

      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* try to reformulate the root expression */
      if( conshdlrdata->reformbinprodsfac > 1 )
      {
         SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, conss[c], consdata->expr, conshdlrdata->reformbinprodsfac, &newexpr, naddconss) );
      }

      /* release the root node if another expression has been found */
      if( newexpr != NULL )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );
         consdata->expr = newexpr;

         /* mark constraint to be not simplified anymore */
         consdata->issimplified = FALSE;
      }

      /* replace each product of binary variables separately */
      SCIP_CALL( replaceBinaryProducts(scip, conshdlr, conss[c], exprmap, it, naddconss, nchgcoefs) );
   }

   /* free memory */
   SCIPhashmapFree(&exprmap);
   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** scales the sides of the constraint \f$\ell \leq \sum_i c_i f_i(x) \leq r\f$.
 *
 *  Let \f$n_+\f$ the number of positive coefficients \f$c_i\f$ and \f$n_-\f$ be the number of negative coefficients.
 *  Then scale by -1 if
 *  - \f$n_+ < n_-\f$, or
 *  - \f$n_+ = n_-\f$ and \f$r = \infty\f$.
 */
static
SCIP_RETCODE scaleConsSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool*            changed             /**< buffer to store if the expression of cons changed */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPisExprSum(scip, consdata->expr) )
   {
      SCIP_Real* coefs;
      SCIP_Real constant;
      int nchildren;
      int counter = 0;

      coefs = SCIPgetCoefsExprSum(consdata->expr);
      constant = SCIPgetConstantExprSum(consdata->expr);
      nchildren = SCIPexprGetNChildren(consdata->expr);

      /* handle special case when constraint is l <= -f(x) <= r and f(x) not a sum: simplfy ensures f is not a sum */
      if( nchildren == 1 && constant == 0.0 && coefs[0] == -1.0 )
      {
         SCIP_EXPR* expr;
         expr = consdata->expr;

         consdata->expr = SCIPexprGetChildren(expr)[0];
         assert(!SCIPisExprSum(scip, consdata->expr));

         SCIPcaptureExpr(consdata->expr);

         SCIPswapReals(&consdata->lhs, &consdata->rhs);
         consdata->lhs = -consdata->lhs;
         consdata->rhs = -consdata->rhs;

         SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
         *changed = TRUE;
         return SCIP_OKAY;
      }

      /* compute n_+ - n_i */
      for( i = 0; i < nchildren; ++i )
         counter += coefs[i] > 0 ? 1 : -1;

      if( counter < 0 || (counter == 0 && SCIPisInfinity(scip, consdata->rhs)) )
      {
         SCIP_EXPR* expr;
         SCIP_Real* newcoefs;

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &newcoefs, nchildren) );

         for( i = 0; i < nchildren; ++i )
            newcoefs[i] = -coefs[i];

         /* create a new sum expression */
         SCIP_CALL( SCIPcreateExprSum(scip, &expr, nchildren, SCIPexprGetChildren(consdata->expr), newcoefs, -constant, exprownerCreate, (void*)conshdlr) );

         /* replace expression in constraint data and scale sides */
         SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );
         consdata->expr = expr;
         SCIPswapReals(&consdata->lhs, &consdata->rhs);
         consdata->lhs = -consdata->lhs;
         consdata->rhs = -consdata->rhs;

         /* free memory */
         SCIPfreeBufferArray(scip, &newcoefs);

         *changed = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** forbid multiaggrations of variables that appear nonlinear in constraints */
static
SCIP_RETCODE forbidNonlinearVariablesMultiaggration(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_EXPRITER* it;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* expr;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   if( !SCIPconshdlrGetData(conshdlr)->forbidmultaggrnlvar )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* if root expression is sum, then forbid multiaggregation only for variables that are not in linear terms of sum,
       *   i.e., skip children of sum that are variables
       */
      if( SCIPisExprSum(scip, consdata->expr) )
      {
         int i;
         SCIP_EXPR* child;
         for( i = 0; i < SCIPexprGetNChildren(consdata->expr); ++i )
         {
            child = SCIPexprGetChildren(consdata->expr)[i];

            /* skip variable expression, as they correspond to a linear term */
            if( SCIPisExprVar(scip, child) )
               continue;

            for( expr = SCIPexpriterRestartDFS(it, child); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
               if( SCIPisExprVar(scip, expr) )
               {
                  SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, SCIPgetVarExprVar(expr)) );
               }
         }
      }
      else
      {
         for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
            if( SCIPisExprVar(scip, expr) )
            {
               SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, SCIPgetVarExprVar(expr)) );
            }
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** simplifies expressions and replaces common subexpressions for a set of constraints
 * @todo put the constant to the constraint sides
 */
static
SCIP_RETCODE canonicalizeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< total number of constraints */
   SCIP_PRESOLTIMING     presoltiming,       /**< presolve timing (SCIP_PRESOLTIMING_ALWAYS if not in presolving) */
   SCIP_Bool*            infeasible,         /**< buffer to store whether infeasibility has been detected */
   int*                  ndelconss,          /**< counter to add number of deleted constraints, or NULL */
   int*                  naddconss,          /**< counter to add number of added constraints, or NULL */
   int*                  nchgcoefs           /**< counter to add number of changed coefficients, or NULL */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int* nlockspos;
   int* nlocksneg;
   SCIP_Bool havechange;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(infeasible != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* update number of canonicalize calls */
   ++(conshdlrdata->ncanonicalizecalls);

   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->canonicalizetime) );

   *infeasible = FALSE;

   /* set havechange to TRUE in the first call of canonicalize; otherwise we might not replace common subexpressions */
   havechange = conshdlrdata->ncanonicalizecalls == 1;

   /* free nonlinear handlers information from expressions */  /* TODO can skip this in first presolve round */
   SCIP_CALL( deinitSolve(scip, conshdlr, conss, nconss) );

   /* allocate memory for storing locks of each constraint */
   SCIP_CALL( SCIPallocBufferArray(scip, &nlockspos, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksneg, nconss) );

   /* unlock all constraints */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* remember locks */
      nlockspos[i] = consdata->nlockspos;
      nlocksneg[i] = consdata->nlocksneg;

      /* remove locks */
      SCIP_CALL( addLocks(scip, conss[i], -consdata->nlockspos, -consdata->nlocksneg) );
      assert(consdata->nlockspos == 0);
      assert(consdata->nlocksneg == 0);
   }

#ifndef NDEBUG
   /* check whether all locks of each expression have been removed */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_EXPR* expr;
      SCIP_EXPRITER* it;

      SCIP_CALL( SCIPcreateExpriter(scip, &it) );

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      SCIP_CALL( SCIPexpriterInit(it, consdata->expr, SCIP_EXPRITER_RTOPOLOGIC, TRUE) );
      for( expr = consdata->expr; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         assert(expr != NULL);
         assert(SCIPexprGetOwnerData(expr)->nlocksneg == 0);
         assert(SCIPexprGetOwnerData(expr)->nlockspos == 0);
      }
      SCIPfreeExpriter(&it);
   }
#endif

   /* reformulate products of binary variables */
   if( conshdlrdata->reformbinprods && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING
      && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) )
   {
      int tmpnaddconss = 0;
      int tmpnchgcoefs = 0;

      /* call this function before simplification because expressions might not be simplified after reformulating
       * binary products; the detection of some nonlinear handlers might assume that expressions are simplified
       */
      SCIP_CALL( presolveBinaryProducts(scip, conshdlr, conss, nconss, &tmpnaddconss, &tmpnchgcoefs) );

      /* update counters */
      if( naddconss != NULL )
         *naddconss += tmpnaddconss;
      if( nchgcoefs != NULL )
         *nchgcoefs += tmpnchgcoefs;

      /* check whether at least one expression has changed */
      if( tmpnaddconss + tmpnchgcoefs > 0 )
         havechange = TRUE;
   }

   for( i = 0; i < nconss; ++i )
   {
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* call simplify for each expression */
      if( !consdata->issimplified && consdata->expr != NULL )
      {
         SCIP_EXPR* simplified;
         SCIP_Bool changed;

         changed = FALSE;
         SCIP_CALL( SCIPsimplifyExpr(scip, consdata->expr, &simplified, &changed, infeasible, exprownerCreate, (void*)conshdlr) );
         consdata->issimplified = TRUE;

         if( changed )
            havechange = TRUE;

         /* If root expression changed, then we need to take care updating the locks as well (the consdata is the one holding consdata->expr "as a child").
          * If root expression did not change, some subexpression may still have changed, but the locks were taking care of in the corresponding SCIPreplaceExprChild() call.
          */
         if( simplified != consdata->expr )
         {
            assert(changed);

            /* release old expression */
            SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );

            /* store simplified expression */
            consdata->expr = simplified;
         }
         else
         {
            /* The simplify captures simplified in any case, also if nothing has changed.
             * Therefore, we have to release it here.
             */
            SCIP_CALL( SCIPreleaseExpr(scip, &simplified) );
         }

         if( *infeasible )
            break;

         /* scale constraint sides */
         SCIP_CALL( scaleConsSides(scip, conshdlr, conss[i], &changed) );

         if( changed )
            havechange = TRUE;

         /* handle constant root expression; either the problem is infeasible or the constraint is redundant */
         if( SCIPisExprValue(scip, consdata->expr) )
         {
            SCIP_Real value = SCIPgetValueExprValue(consdata->expr);
            if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasNegative(scip, value - consdata->lhs)) ||
                (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasPositive(scip, value - consdata->rhs)) )
            {
               SCIPdebugMsg(scip, "<%s> with constant expression found infeasible\n", SCIPconsGetName(conss[i]));
               SCIPdebugPrintCons(scip, conss[i], NULL);
               *infeasible = TRUE;
               break;
            }
            else
            {
               SCIP_CALL( addLocks(scip, conss[i], nlockspos[i], nlocksneg[i]) );
               SCIP_CALL( SCIPdelCons(scip, conss[i]) );
               if( ndelconss != NULL )
                  ++*ndelconss;
               havechange = TRUE;
            }
         }
      }
   }

   /* replace common subexpressions */
   if( havechange && !*infeasible )
   {
      SCIP_CONS** consssorted;
      SCIP_EXPR** rootexprs;
      SCIP_Bool replacedroot;

      SCIP_CALL( SCIPallocBufferArray(scip, &rootexprs, nconss) );
      for( i = 0; i < nconss; ++i )
         rootexprs[i] = SCIPconsGetData(conss[i])->expr;

      SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, rootexprs, nconss, &replacedroot) );

      /* update pointer to root expr in constraints, if any has changed
       * SCIPreplaceCommonSubexpressions will have released the old expr and captures the new one
       */
      if( replacedroot )
         for( i = 0; i < nconss; ++i )
            SCIPconsGetData(conss[i])->expr = rootexprs[i];

      SCIPfreeBufferArray(scip, &rootexprs);

      /* TODO this is a possibly expensive way to update the variable expressions stored inside an expression which might have
       * been changed after simplification; now we completely recollect all variable expression and variable events
       */

      /* Each variable stores the constraints for which it catched varbound events sorted by the constraint index.
       * Thus, for performance reasons, it is better to call dropVarEvents in descending order of constraint index.
       */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consssorted, conss, nconss) );
      SCIPsortPtr((void**)consssorted, compIndexConsNonlinear, nconss);

      for( i = nconss-1; i >= 0; --i )
      {
         assert(i == 0 || compIndexConsNonlinear((void*)consssorted[i-1], (void*)consssorted[i]) < 0);
         if( SCIPconsIsDeleted(consssorted[i]) )
            continue;

         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
         SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(consssorted[i])) );
      }
      for( i = 0; i < nconss; ++i )
      {
         if( SCIPconsIsDeleted(consssorted[i]) )
            continue;

         SCIP_CALL( storeVarExprs(scip, conshdlr, SCIPconsGetData(consssorted[i])) );
         SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
      }

      SCIPfreeBufferArray(scip, &consssorted);

      /* forbid multiaggregation for nonlinear variables again (in case new variables appeared now)
       * a multiaggregation of a nonlinear variable can yield to a large increase in expressions due to
       * expanding terms in simplify, e.g. ,(sum_i x_i)^2, so we just forbid these
       */
      SCIP_CALL( forbidNonlinearVariablesMultiaggration(scip, conshdlr, conss, nconss) );
   }

   /* restore locks */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsIsDeleted(conss[i]) )
         continue;

      SCIP_CALL( addLocks(scip, conss[i], nlockspos[i], nlocksneg[i]) );
   }

   /* run nlhdlr detect if in presolving stage (that is, not in exitpre)
    * TODO can we skip this in presoltiming fast?
    */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && !*infeasible )
   {
      /* reset one of the number of detections counter to count only current presolving round */
      for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
         SCIPnlhdlrResetNDetectionslast(conshdlrdata->nlhdlrs[i]);

      SCIP_CALL( initSolve(scip, conshdlr, conss, nconss) );
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &nlocksneg);
   SCIPfreeBufferArray(scip, &nlockspos);

   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->canonicalizetime) );

   return SCIP_OKAY;
}

/** merges constraints that have the same root expression */
static
SCIP_RETCODE presolveMergeConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            success             /**< pointer to store whether at least one constraint could be deleted */
   )
{
   SCIP_HASHMAP* expr2cons;
   SCIP_Bool* updatelocks;
   int* nlockspos;
   int* nlocksneg;
   int c;

   assert(success != NULL);

   *success = FALSE;

   /* not enough constraints available */
   if( nconss <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapCreate(&expr2cons, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &updatelocks, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlockspos, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksneg, nconss) );

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      /* ignore deleted constraints */
      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* add expression to the hash map if not seen so far */
      if( !SCIPhashmapExists(expr2cons, (void*)consdata->expr) )
      {
         SCIP_CALL( SCIPhashmapInsertInt(expr2cons, (void*)consdata->expr, c) );
      }
      else
      {
         SCIP_CONSDATA* imgconsdata;
         int idx;

         idx = SCIPhashmapGetImageInt(expr2cons, (void*)consdata->expr);
         assert(idx >= 0 && idx < nconss);

         imgconsdata = SCIPconsGetData(conss[idx]);
         assert(imgconsdata != NULL);
         assert(imgconsdata->expr == consdata->expr);

         SCIPdebugMsg(scip, "merge constraint %g <= %s <= %g with %g <= %s <= %g\n", consdata->lhs,
            SCIPconsGetName(conss[c]), consdata->rhs, imgconsdata->lhs, SCIPconsGetName(conss[idx]), imgconsdata->rhs);

         /* check whether locks need to be updated */
         if( !updatelocks[idx] && ((SCIPisInfinity(scip, -imgconsdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs))
            || (SCIPisInfinity(scip, imgconsdata->rhs) && !SCIPisInfinity(scip, consdata->rhs))) )
         {
            nlockspos[idx] = imgconsdata->nlockspos;
            nlocksneg[idx] = imgconsdata->nlocksneg;
            SCIP_CALL( addLocks(scip, conss[idx], -imgconsdata->nlockspos, -imgconsdata->nlocksneg) );
            updatelocks[idx] = TRUE;
         }

         /* update constraint sides */
         imgconsdata->lhs = MAX(imgconsdata->lhs, consdata->lhs);
         imgconsdata->rhs = MIN(imgconsdata->rhs, consdata->rhs);

         /* delete constraint */
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         *success = TRUE;
      }
   }

   /* restore locks of updated constraints */
   if( *success )
   {
      for( c = 0; c < nconss; ++c )
      {
         if( updatelocks[c] )
         {
            SCIP_CALL( addLocks(scip, conss[c], nlockspos[c], nlocksneg[c]) );
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &nlocksneg);
   SCIPfreeBufferArray(scip, &nlockspos);
   SCIPfreeBufferArray(scip, &updatelocks);
   SCIPhashmapFree(&expr2cons);

   return SCIP_OKAY;
}

/** interval evaluation of variables as used in redundancy check
 *
 * Returns local variable bounds of a variable, relaxed by feastol, as interval.
 */
static
SCIP_DECL_EXPR_INTEVALVAR(intEvalVarRedundancyCheck)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL interval;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(var != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)intevalvardata;
   assert(conshdlrdata != NULL);

   if( conshdlrdata->globalbounds )
   {
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
   }
   else
   {
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
   }
   assert(lb <= ub);  /* can SCIP ensure by now that variable bounds are not contradicting? */

   /* relax variable bounds, if there are bounds and variable is not fixed
    * (actually some assert complains if trying SCIPisRelEQ if both bounds are at different infinity)
    */
   if( !(SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub)) && !SCIPisRelEQ(scip, lb, ub) )
   {
      if( !SCIPisInfinity(scip, -lb) )
         lb -= SCIPfeastol(scip);

      if( !SCIPisInfinity(scip, ub) )
         ub += SCIPfeastol(scip);
   }

   /* convert SCIPinfinity() to SCIP_INTERVAL_INFINITY */
   lb = -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb);
   ub =  infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  ub);
   assert(lb <= ub);

   SCIPintervalSetBounds(&interval, lb, ub);

   return interval;
}

/** removes constraints that are always feasible or very simple
 *
 * Checks whether the activity of constraint functions is a subset of the constraint sides (relaxed by feastol).
 * To compute the activity, we use forwardPropExpr(), but relax variable bounds by feastol, because solutions to be checked
 * might violate variable bounds by up to feastol, too.
 * This is the main reason why the redundancy check is not done in propConss(), which relaxes variable bounds by epsilon only.
 *
 * Also removes constraints of the form lhs &le; variable &le; rhs.
 *
 * @todo it would be sufficient to check constraints for which we know that they are not currently violated by a valid solution
 *
 * @note This could should not run during solving, because the forwardProp takes the bounds of auxiliary variables into account.
 * For the root expression, these bounds are already set to the constraint sides, so that the activity of every expression
 * would appear as if the constraint is redundant.
 */
static
SCIP_RETCODE presolveRedundantConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_Bool*            cutoff,             /**< pointer to store whether infeasibility has been identified */
   int*                  ndelconss,          /**< buffer to add the number of deleted constraints */
   int*                  nchgbds             /**< buffer to add the number of variable bound tightenings */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL activity;
   SCIP_INTERVAL sides;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgbds != NULL);

   /* no constraints to check */
   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* increase curboundstag and set lastvaractivitymethodchange
    * we do this here to trigger a reevaluation of all variable bounds, since we will relax variable bounds
    * for the redundancy check differently than for domain propagation
    * we also update lastboundrelax to ensure activites of all expressions are indeed reevaluated
    */
   ++conshdlrdata->curboundstag;
   assert(conshdlrdata->curboundstag > 0);
   conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
   conshdlrdata->intevalvar = intEvalVarRedundancyCheck;

   SCIPdebugMsg(scip, "checking %d constraints for redundancy\n", nconss);

   *cutoff = FALSE;
   for( i = 0; i < nconss; ++i )
   {
      if( !SCIPconsIsActive(conss[i]) || SCIPconsIsDeleted(conss[i]) )
         continue;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* handle constant expressions separately: either the problem is infeasible or the constraint is redundant */
      if( SCIPisExprValue(scip, consdata->expr) )
      {
         SCIP_Real value = SCIPgetValueExprValue(consdata->expr);

         if( (!SCIPisInfinity(scip, -consdata->lhs) && value < consdata->lhs - SCIPfeastol(scip)) ||
             (!SCIPisInfinity(scip,  consdata->rhs) && value > consdata->rhs + SCIPfeastol(scip)) )
         {
            SCIPdebugMsg(scip, "constant constraint <%s> is infeasible: %g in [%g,%g] ", SCIPconsGetName(conss[i]), value, consdata->lhs, consdata->rhs);
            *cutoff = TRUE;

            goto TERMINATE;
         }

         SCIPdebugMsg(scip, "constant constraint <%s> is redundant: %g in [%g,%g] ", SCIPconsGetName(conss[i]), value, consdata->lhs, consdata->rhs);

         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      /* handle variable expressions separately: tighten variable bounds to constraint sides, then remove constraint (now redundant) */
      if( SCIPisExprVar(scip, consdata->expr) )
      {
         SCIP_VAR* var;
         SCIP_Bool tightened;

         var = SCIPgetVarExprVar(consdata->expr);
         assert(var != NULL);

         SCIPdebugMsg(scip, "variable constraint <%s> can be made redundant: <%s>[%g,%g] in [%g,%g]\n", SCIPconsGetName(conss[i]), SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), consdata->lhs, consdata->rhs);

         /* ensure that variable bounds are within constraint sides */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, consdata->lhs, TRUE, cutoff, &tightened) );

            if( tightened )
               ++*nchgbds;

            if( *cutoff )
               goto TERMINATE;
         }

         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, consdata->rhs, TRUE, cutoff, &tightened) );

            if( tightened )
               ++*nchgbds;

            if( *cutoff )
               goto TERMINATE;
         }

         /* delete the (now) redundant constraint locally */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      /* reevaluate expression activity, now using intEvalVarRedundancyCheck
       * we relax variable bounds by feastol here, as solutions that are checked later can also violate
       * variable bounds by up to feastol
       * (relaxing fixed variables seems to be too much, but they would be removed by presolve soon anyway)
       */
      SCIPdebugMsg(scip, "call forwardPropExpr() for constraint <%s>: ", SCIPconsGetName(conss[i]));
      SCIPdebugPrintCons(scip, conss[i], NULL);

      SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, FALSE, cutoff, NULL) );
      assert(*cutoff || !SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(consdata->expr)));

      /* it is unlikely that we detect infeasibility by doing forward propagation */
      if( *cutoff )
      {
         SCIPdebugMsg(scip, " -> cutoff\n");
         goto TERMINATE;
      }

      assert(SCIPexprGetActivityTag(consdata->expr) == conshdlrdata->curboundstag);
      activity = SCIPexprGetActivity(consdata->expr);

      /* relax sides by feastol
       * we could accept every solution that violates constraints up to feastol as redundant, so this is the most permissive we can be
       */
      SCIPintervalSetBounds(&sides,
         SCIPisInfinity(scip, -consdata->lhs) ? -SCIP_INTERVAL_INFINITY : consdata->lhs - SCIPfeastol(scip),
         SCIPisInfinity(scip,  consdata->rhs) ?  SCIP_INTERVAL_INFINITY : consdata->rhs + SCIPfeastol(scip));

      if( SCIPintervalIsSubsetEQ(SCIP_INTERVAL_INFINITY, activity, sides) )
      {
         SCIPdebugMsg(scip, " -> redundant: activity [%g,%g] within sides [%g,%g]\n", activity.inf, activity.sup, consdata->lhs, consdata->rhs);

         SCIP_CALL( SCIPdelConsLocal(scip, conss[i]) );
         ++*ndelconss;

         continue;
      }

      SCIPdebugMsg(scip, " -> not redundant: activity [%g,%g] not within sides [%g,%g]\n", activity.inf, activity.sup, consdata->lhs, consdata->rhs);
   }

TERMINATE:
   /* make sure all activities are reevaluated again, since we relaxed bounds in a different way */
   ++conshdlrdata->curboundstag;
   conshdlrdata->lastvaractivitymethodchange = conshdlrdata->curboundstag;
   conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
   conshdlrdata->intevalvar = intEvalVarBoundTightening;

   return SCIP_OKAY;
}

/** tries to automatically convert a nonlinear constraint into a more specific and more specialized constraint */
static
SCIP_RETCODE presolveUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_Bool*            upgraded,           /**< buffer to store whether constraint was upgraded */
   int*                  nupgdconss,         /**< buffer to increase if constraint was upgraded */
   int*                  naddconss           /**< buffer to increase with number of additional constraints created during upgrade */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS** upgdconss;
   int upgdconsssize;
   int nupgdconss_;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsModifiable(cons));
   assert(upgraded   != NULL);
   assert(nupgdconss != NULL);
   assert(naddconss  != NULL);

   *upgraded = FALSE;

   nupgdconss_ = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if there are no upgrade methods, we can stop */
   if( conshdlrdata->nconsupgrades == 0 )
      return SCIP_OKAY;

   upgdconsssize = 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &upgdconss, upgdconsssize) );

   /* call the upgrading methods */
   SCIPdebugMsg(scip, "upgrading nonlinear constraint <%s> (up to %d upgrade methods): ", SCIPconsGetName(cons), conshdlrdata->nconsupgrades);
   SCIPdebugPrintCons(scip, cons, NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nconsupgrades; ++i )
   {
      if( !conshdlrdata->consupgrades[i]->active )
         continue;

      assert(conshdlrdata->consupgrades[i]->consupgd != NULL);

      SCIP_CALL( conshdlrdata->consupgrades[i]->consupgd(scip, cons, consdata->nvarexprs, &nupgdconss_, upgdconss, upgdconsssize) );

      while( nupgdconss_ < 0 )
      {
         /* upgrade function requires more memory: resize upgdconss and call again */
         assert(-nupgdconss_ > upgdconsssize);
         upgdconsssize = -nupgdconss_;
         SCIP_CALL( SCIPreallocBufferArray(scip, &upgdconss, -nupgdconss_) );

         SCIP_CALL( conshdlrdata->consupgrades[i]->consupgd(scip, cons, consdata->nvarexprs, &nupgdconss_, upgdconss, upgdconsssize) );

         assert(nupgdconss_ != 0);
      }

      if( nupgdconss_ > 0 )
      {
         /* got upgrade */
         int j;

         SCIPdebugMsg(scip, " -> upgraded to %d constraints:\n", nupgdconss_);

         /* add the upgraded constraints to the problem and forget them */
         for( j = 0; j < nupgdconss_; ++j )
         {
            SCIPdebugMsgPrint(scip, "\t");
            SCIPdebugPrintCons(scip, upgdconss[j], NULL);

            SCIP_CALL( SCIPaddCons(scip, upgdconss[j]) );      /*lint !e613*/
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconss[j]) ); /*lint !e613*/
         }

         /* count the first upgrade constraint as constraint upgrade and the remaining ones as added constraints */
         *nupgdconss += 1;
         *naddconss += nupgdconss_ - 1;
         *upgraded = TRUE;

         /* delete upgraded constraint */
         SCIPdebugMsg(scip, "delete constraint <%s> after upgrade\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );

         break;
      }
   }

   SCIPfreeBufferArray(scip, &upgdconss);

   return SCIP_OKAY;
}

/** returns whether the variable of a given variable expression is a candidate for presolveSingleLockedVars(), i.e.,
 *  the variable is only contained in a single nonlinear constraint, has no objective coefficient, has finite
 *  variable bounds, and is not binary
 */
static
SCIP_Bool isSingleLockedCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< variable expression */
   )
{
   SCIP_VAR* var;
   SCIP_EXPR_OWNERDATA* ownerdata;

   assert(SCIPisExprVar(scip, expr));

   var = SCIPgetVarExprVar(expr);
   assert(var != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   return SCIPvarGetNLocksDownType(var, SCIP_LOCKTYPE_MODEL) == ownerdata->nlocksneg
      && SCIPvarGetNLocksUpType(var, SCIP_LOCKTYPE_MODEL) == ownerdata->nlockspos
      && ownerdata->nconss == 1 && SCIPisZero(scip, SCIPvarGetObj(var))
      && !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var))
      && SCIPvarGetType(var) != SCIP_VARTYPE_BINARY
      && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
}

/** removes all variable expressions that are contained in a given expression from a hash map */
static
SCIP_RETCODE removeSingleLockedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRITER*        it,                 /**< expression iterator */
   SCIP_HASHMAP*         exprcands           /**< map to hash variable expressions */
   )
{
   SCIP_EXPR* e;

   for( e = SCIPexpriterRestartDFS(it, expr); !SCIPexpriterIsEnd(it); e = SCIPexpriterGetNext(it) )
   {
      if( SCIPisExprVar(scip, e) && SCIPhashmapExists(exprcands, (void*)e) )
      {
         SCIP_CALL( SCIPhashmapRemove(exprcands, (void*)e) );
      }
   }

   return SCIP_OKAY;
}

/** presolving method to fix a variable \f$x_i\f$ to one of its bounds if the variable is only contained in a single
 *  nonlinear constraint g(x) &le; rhs (&ge; lhs) if g() is concave (convex) in \f$x_i\f$
 *
 *  If a continuous variable has bounds [0,1], then the variable type is changed to be binary.
 *  Otherwise, a bound disjunction constraint is added.
 *
 *  @todo the same reduction can be applied if g(x) is not concave, but monotone in \f$x_i\f$ for g(x) &le; rhs
 *  @todo extend this to cases where a variable can appear in a monomial with an exponent, essentially relax
 *    g(x) to \f$\sum_i [a_i,b_i] x^{p_i}\f$ for a single variable \f$x\f$ and try to conclude montonicity or convexity/concavity
 *    on this (probably have one or two flags per variable and update this whenever another \f$x^{p_i}\f$ is found)
 *  @todo reduction should also be applicable if variable appears in the objective with the right sign (sign such that opt is at boundary)
 */
static
SCIP_RETCODE presolveSingleLockedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   int*                  nchgvartypes,       /**< pointer to store the total number of changed variable types */
   int*                  naddconss,          /**< pointer to store the total number of added constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether problem is infeasible */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR** singlelocked;
   SCIP_HASHMAP* exprcands;
   SCIP_Bool hasbounddisj;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int nsinglelocked = 0;
   int i;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(nchgvartypes != NULL);
   assert(naddconss != NULL);
   assert(infeasible != NULL);

   *nchgvartypes = 0;
   *naddconss = 0;
   *infeasible = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* only consider constraints with one finite side */
   if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
      return SCIP_OKAY;

   /* only consider sum expressions */
   if( !SCIPisExprSum(scip, consdata->expr) )
      return SCIP_OKAY;

   /* remember which side is finite */
   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   /* allocate memory */
   SCIP_CALL( SCIPhashmapCreate(&exprcands, SCIPblkmem(scip), consdata->nvarexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &singlelocked, consdata->nvarexprs) );

   /* check all variable expressions for single locked variables */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);

      if( isSingleLockedCand(scip, consdata->varexprs[i]) )
      {
         SCIP_CALL( SCIPhashmapInsert(exprcands, (void*)consdata->varexprs[i], NULL) );
         singlelocked[nsinglelocked++] = consdata->varexprs[i];
      }
   }
   SCIPdebugMsg(scip, "found %d single locked variables for constraint %s\n", nsinglelocked, SCIPconsGetName(cons));

   if( nsinglelocked > 0 )
   {
      SCIP_EXPR** children;
      SCIP_EXPRITER* it;
      int nchildren;

      children = SCIPexprGetChildren(consdata->expr);
      nchildren = SCIPexprGetNChildren(consdata->expr);

      /* create iterator */
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
      SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
      SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR);

      for( i = 0; i < nchildren; ++i )
      {
         SCIP_EXPR* child;
         SCIP_Real coef;

         child = children[i];
         assert(child != NULL);
         coef = SCIPgetCoefsExprSum(consdata->expr)[i];

         /* ignore linear terms */
         if( SCIPisExprVar(scip, child) )
            continue;

         /* consider products prod_j f_j(x); ignore f_j(x) if it is a single variable, otherwise iterate through the
          * expression that represents f_j and remove each variable expression from exprcands
          */
         else if( SCIPisExprProduct(scip, child) )
         {
            int j;

            for( j = 0; j < SCIPexprGetNChildren(child); ++j )
            {
               SCIP_EXPR* grandchild = SCIPexprGetChildren(child)[j];

               if( !SCIPisExprVar(scip, grandchild) )
               {
                  /* mark all variable expressions that are contained in the expression */
                  SCIP_CALL( removeSingleLockedVars(scip, grandchild, it, exprcands) );
               }
            }
         }
         /* fixing a variable x to one of its bounds is only valid for ... +x^p >= lhs or ... -x^p <= rhs if p = 2k
          * for an integer k >= 1
          */
         else if( SCIPisExprPower(scip, child) )
         {
            SCIP_EXPR* grandchild = SCIPexprGetChildren(child)[0];
            SCIP_Real exponent = SCIPgetExponentExprPow(child);
            SCIP_Bool valid;

            /* check for even integral exponent */
            valid = exponent > 1.0 && fmod(exponent, 2.0) == 0.0;

            if( !valid || !SCIPisExprVar(scip, grandchild) || (hasrhs && coef > 0.0) || (haslhs && coef < 0.0) )
            {
               /* mark all variable expressions that are contained in the expression */
               SCIP_CALL( removeSingleLockedVars(scip, grandchild, it, exprcands) );
            }
         }
         /* all other cases cannot be handled */
         else
         {
            /* mark all variable expressions that are contained in the expression */
            SCIP_CALL( removeSingleLockedVars(scip, child, it, exprcands) );
         }
      }

      /* free expression iterator */
      SCIPfreeExpriter(&it);
   }

   /* check whether the bound disjunction constraint handler is available */
   hasbounddisj = SCIPfindConshdlr(scip, "bounddisjunction") != NULL;

   /* fix variable to one of its bounds by either changing its variable type or adding a disjunction constraint */
   for( i = 0; i < nsinglelocked; ++i )
   {
      /* only consider expressions that are still contained in the exprcands map */
      if( SCIPhashmapExists(exprcands, (void*)singlelocked[i]) )
      {
         SCIP_CONS* newcons;
         SCIP_VAR* vars[2];
         SCIP_BOUNDTYPE boundtypes[2];
         SCIP_Real bounds[2];
         char name[SCIP_MAXSTRLEN];
         SCIP_VAR* var;

         var = SCIPgetVarExprVar(singlelocked[i]);
         assert(var != NULL);
         SCIPdebugMsg(scip, "found single locked variable %s in [%g,%g] that can be fixed to one of its bounds\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));

         /* try to change the variable type to binary */
         if( conshdlrdata->checkvarlocks == 't' && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0) && SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0) )
         {
            assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            ++(*nchgvartypes);

            if( *infeasible )
            {
               SCIPdebugMsg(scip, "detect infeasibility after changing variable type of <%s>\n", SCIPvarGetName(var));
               break;
            }
         }
         /* add bound disjunction constraint if bounds of the variable are finite */
         else if( hasbounddisj && !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            vars[0] = var;
            vars[1] = var;
            boundtypes[0] = SCIP_BOUNDTYPE_LOWER;
            boundtypes[1] = SCIP_BOUNDTYPE_UPPER;
            bounds[0] = SCIPvarGetUbGlobal(var);
            bounds[1] = SCIPvarGetLbGlobal(var);

            SCIPdebugMsg(scip, "add bound disjunction constraint for %s\n", SCIPvarGetName(var));

            /* create, add, and release bound disjunction constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quadvarbnddisj_%s", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &newcons, name, 2, vars, boundtypes, bounds, TRUE, TRUE,
               TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, newcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
            ++(*naddconss);
         }
      }
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &singlelocked);
   SCIPhashmapFree(&exprcands);

   return SCIP_OKAY;
}

/** presolving method to check if there is a single linear continuous variable that can be made implicit integer */
static
SCIP_RETCODE presolveImplint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< nonlinear constraints */
   int                   nconss,             /**< total number of nonlinear constraints */
   int*                  nchgvartypes,       /**< pointer to update the total number of changed variable types */
   SCIP_Bool*            infeasible          /**< pointer to store whether problem is infeasible */
   )
{
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nchgvartypes != NULL);
   assert(infeasible != NULL);

   *infeasible = FALSE;

   /* nothing can be done if there are no binary and integer variables available */
   if( SCIPgetNBinVars(scip) == 0 && SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   /* no continuous var can be made implicit-integer if there are no continuous variables */
   if( SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_EXPR** children;
      int nchildren;
      SCIP_Real* coefs;
      SCIP_EXPR* cand = NULL;
      SCIP_Real candcoef = 0.0;
      int i;

      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* the constraint must be an equality constraint */
      if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
         continue;

      /* the root expression needs to be a sum expression */
      if( !SCIPisExprSum(scip, consdata->expr) )
         continue;

      children = SCIPexprGetChildren(consdata->expr);
      nchildren = SCIPexprGetNChildren(consdata->expr);

      /* the sum expression must have at least two children
       * (with one child, we would look for a coef*x = constant, which is presolved away anyway)
       */
      if( nchildren <= 1 )
         continue;

      coefs = SCIPgetCoefsExprSum(consdata->expr);

      /* find first continuous variable and get value of its coefficient */
      for( i = 0; i < nchildren; ++i )
      {
         if( !SCIPisExprVar(scip, children[i]) || SCIPvarIsIntegral(SCIPgetVarExprVar(children[i])) )
            continue;

         candcoef = coefs[i];
         assert(candcoef != 0.0);

         /* lhs/rhs - constant divided by candcoef must be integral
          * if not, break with cand == NULL, so give up
          */
         if( SCIPisIntegral(scip, (consdata->lhs - SCIPgetConstantExprSum(consdata->expr)) / candcoef) )
            cand = children[i];

         break;
      }

      /* no suitable continuous variable found */
      if( cand == NULL )
         continue;

      /* check whether all other coefficients are integral when diving by candcoef and all other children are integral */
      for( i = 0; i < nchildren; ++i )
      {
         if( children[i] == cand )
            continue;

         /* child i must be integral */
         if( !SCIPexprIsIntegral(children[i]) )
         {
            cand = NULL;
            break;
         }

         /* coefficient of child i must be integral if diving by candcoef */
         if( !SCIPisIntegral(scip, coefs[i] / candcoef) )  /*lint !e414*/
         {
            cand = NULL;
            break;
         }
      }

      if( cand == NULL )
         continue;

      SCIPdebugMsg(scip, "make variable <%s> implicit integer due to constraint <%s>\n",
         SCIPvarGetName(SCIPgetVarExprVar(cand)), SCIPconsGetName(conss[c]));

      /* change variable type */
      SCIP_CALL( SCIPchgVarType(scip, SCIPgetVarExprVar(cand), SCIP_VARTYPE_IMPLINT, infeasible) );

      if( *infeasible )
         return SCIP_OKAY;

      /* mark expression as being integral (as would be done by expr_var.c in the next round of updating integrality info) */
      SCIPexprSetIntegrality(cand, TRUE);
   }

   return SCIP_OKAY;
}

/** creates auxiliary variable for a given expression
 *
 * @note for a variable expression it does nothing
 * @note this function can only be called in stage SCIP_STAGE_SOLVING
 */
static
SCIP_RETCODE createAuxVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VARTYPE vartype;
   SCIP_INTERVAL activity;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   assert(ownerdata->nauxvaruses > 0);

   /* if we already have auxvar, then do nothing */
   if( ownerdata->auxvar != NULL )
      return SCIP_OKAY;

   /* if expression is a variable-expression, then do nothing */
   if( SCIPisExprVar(scip, expr) )
      return SCIP_OKAY;

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
   {
      SCIPerrorMessage("it is not possible to create auxiliary variables during stage=%d\n", SCIPgetStage(scip));
      return SCIP_INVALIDCALL;
   }

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->auxvarid >= 0);

   /* it doesn't harm much to have an auxvar for a constant, as this can be handled well by the default hdlr,
    * but it usually indicates a missing simplify
    * if we find situations where we need to have an auxvar for a constant, then remove this assert
    */
   assert(!SCIPisExprValue(scip, expr));

   /* create and capture auxiliary variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "auxvar_%s_%d", SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), conshdlrdata->auxvarid);
   ++conshdlrdata->auxvarid;

   /* type of auxiliary variable depends on integrality information of the expression */
   vartype = SCIPexprIsIntegral(expr) ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS;

   /* get activity of expression to initialize variable bounds, if something valid is available (evalActivity was called in initSepa) */
   if( SCIPexprGetActivityTag(expr) >= conshdlrdata->lastboundrelax )
   {
      activity = SCIPexprGetActivity(expr);
      /* we cannot handle a domain error here at the moment, but it seems unlikely that it could occur
       * if it appear, then we could change code to handle this properly, but for now we just ensure that we continue correctly
       * and abort in debug mode only
       */
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, activity) )
      {
         SCIPABORT();
         SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &activity);
      }
   }
   else
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &activity);

   /* if root node, then activity is globally valid, so use it to initialize the global bounds of the auxvar
    * otherwise, we create var without bounds here and use activity to set local bounds below (needs to be after adding var)
    */
   if( SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &ownerdata->auxvar, name, MAX(-SCIPinfinity(scip), activity.inf), MIN(SCIPinfinity(scip), activity.sup), 0.0, vartype) );
   }
   else
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &ownerdata->auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, vartype) );
   }

   /* mark the auxiliary variable to be added for the relaxation only
    * this prevents SCIP to create linear constraints from cuts or conflicts that contain auxiliary variables,
    * or to copy the variable to a subscip
    */
   SCIPvarMarkRelaxationOnly(ownerdata->auxvar);

   SCIP_CALL( SCIPaddVar(scip, ownerdata->auxvar) );

   SCIPdebugMsg(scip, "added auxiliary variable <%s> [%g,%g] for expression %p\n", SCIPvarGetName(ownerdata->auxvar), SCIPvarGetLbGlobal(ownerdata->auxvar), SCIPvarGetUbGlobal(ownerdata->auxvar), (void*)expr);

   /* add variable locks in both directions
    * TODO should be sufficient to lock only according to expr->nlockspos/neg,
    *   but then we need to also update the auxvars locks when the expr locks change
    */
   SCIP_CALL( SCIPaddVarLocks(scip, ownerdata->auxvar, 1, 1) );

#ifdef WITH_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      /* store debug solution value of auxiliary variable
       * assumes that expression has been evaluated in debug solution before
       */
      SCIP_CALL( SCIPdebugAddSolVal(scip, ownerdata->auxvar, SCIPexprGetEvalValue(expr)) );
   }
#endif

   if( SCIPgetDepth(scip) > 0 )
   {
      /* initialize local bounds to (locally valid) activity */
      SCIP_Bool cutoff;
      SCIP_CALL( tightenAuxVarBounds(scip, ownerdata->conshdlr, expr, activity, &cutoff, NULL) );
      assert(!cutoff);  /* should not happen as activity wasn't empty and variable is new */
   }

   return SCIP_OKAY;
}

/** initializes separation for constraint
 *
 * - ensures that activities are up to date in all expressions
 * - creates auxiliary variables where required
 * - calls propExprDomains() to possibly tighten auxvar bounds
 * - calls separation initialization callback of nlhdlrs
 */
static
SCIP_RETCODE initSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether the problem is infeasible or not */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_RESULT result;
   SCIP_VAR* auxvar;
   int nreductions = 0;
   int c, e;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(infeasible != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* start with new propbounds (just to be sure, should not be needed) */
   ++conshdlrdata->curpropboundstag;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   /* first ensure activities are up to date and create auxvars */
   *infeasible = FALSE;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_SOL* debugsol;

         SCIP_CALL( SCIPdebugGetSol(scip, &debugsol) );

         if( debugsol != NULL ) /* it can be compiled WITH_DEBUG_SOLUTION, but still no solution given */
         {
            /* evaluate expression in debug solution, so we can set the solution value of created auxiliary variables
             * in createAuxVar()
             */
            SCIP_CALL( SCIPevalExpr(scip, consdata->expr, debugsol, 0) );
         }
      }
#endif

      /* ensure we have a valid activity for auxvars and propExprDomains() call below */
      SCIP_CALL( SCIPevalExprActivity(scip, consdata->expr) );

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         if( SCIPexprGetOwnerData(expr)->nauxvaruses > 0 )
         {
            SCIP_CALL( createAuxVar(scip, expr) );
         }
      }

      auxvar = SCIPexprGetOwnerData(consdata->expr)->auxvar;
      if( auxvar != NULL )
      {
         SCIPdebugMsg(scip, "tighten auxvar <%s> bounds using constraint sides [%g,%g]\n",
               SCIPvarGetName(auxvar), consdata->lhs, consdata->rhs);
         /* change the bounds of the auxiliary variable of the root node to [lhs,rhs] */
         SCIP_CALL( SCIPtightenVarLb(scip, auxvar, consdata->lhs, TRUE, infeasible, NULL) );
         if( *infeasible )
         {
            SCIPdebugMsg(scip, "infeasibility detected while tightening auxvar lb (%g) using lhs of constraint (%g)\n", SCIPvarGetLbLocal(auxvar), consdata->lhs);
            break;
         }

         SCIP_CALL( SCIPtightenVarUb(scip, auxvar, consdata->rhs, TRUE, infeasible, NULL) );
         if( *infeasible )
         {
            SCIPdebugMsg(scip, "infeasibility detected while tightening auxvar ub (%g) using rhs of constraint (%g)\n", SCIPvarGetUbLocal(auxvar), consdata->rhs);
            break;
         }
      }
   }

   /* now run a special version of reverseprop to ensure that important bound information (like function domains) is stored in bounds of auxvars,
    * since sometimes they cannot be recovered from activity evaluation even after some rounds of domain propagation
    * (e.g., log(x*y), which becomes log(w), w=x*y
    *  log(w) implies w >= 0, but we may not be able to derive bounds on x and y such that w >= 0 is ensured)
    */
   SCIP_CALL( propExprDomains(scip, conshdlr, conss, nconss, &result, &nreductions) );
   if( result == SCIP_CUTOFF )
      *infeasible = TRUE;

   /* now call initsepa of nlhdlrs
    * TODO skip if !SCIPconsIsInitial(conss[c]) ?
    *   but at the moment, initSepa() is called from INITLP anyway, so we have SCIPconsIsInitial(conss[c]) anyway
    */
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   for( c = 0; c < nconss && !*infeasible; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->expr != NULL);

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it) && !*infeasible; expr = SCIPexpriterGetNext(it) )
      {
         SCIP_EXPR_OWNERDATA* ownerdata;

         ownerdata = SCIPexprGetOwnerData(expr);
         assert(ownerdata != NULL);

         if( ownerdata->nauxvaruses == 0 )
            continue;

         for( e = 0; e < ownerdata->nenfos; ++e )
         {
            SCIP_NLHDLR* nlhdlr;
            SCIP_Bool underestimate;
            SCIP_Bool overestimate;
            assert(ownerdata->enfos[e] != NULL);

            /* skip if initsepa was already called, e.g., because this expression is also part of a constraint
             * which participated in a previous initSepa() call
             */
            if( ownerdata->enfos[e]->issepainit )
               continue;

            /* only call initsepa if it will actually separate */
            if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 )
               continue;

            nlhdlr = ownerdata->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* only init sepa if there is an initsepa callback */
            if( !SCIPnlhdlrHasInitSepa(nlhdlr) )
               continue;

            /* check whether expression needs to be under- or overestimated */
            overestimate = ownerdata->nlocksneg > 0;
            underestimate = ownerdata->nlockspos > 0;
            assert(underestimate || overestimate);

            SCIPdebugMsg(scip, "initsepa under=%u over=%u for expression %p\n", underestimate, overestimate, (void*)expr);

            /* call the separation initialization callback of the nonlinear handler */
            SCIP_CALL( SCIPnlhdlrInitsepa(scip, conshdlr, conss[c], nlhdlr, expr,
               ownerdata->enfos[e]->nlhdlrexprdata, overestimate, underestimate, infeasible) );
            ownerdata->enfos[e]->issepainit = TRUE;

            if( *infeasible )
            {
               /* stop everything if we detected infeasibility */
               SCIPdebugMsg(scip, "detect infeasibility for constraint %s during initsepa()\n", SCIPconsGetName(conss[c]));
               break;
            }
         }
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** returns whether we are ok to branch on auxiliary variables
 *
 * Currently returns whether depth of node in B&B tree is at least value of constraints/nonlinear/branching/aux parameter.
 */
static
SCIP_Bool branchAuxNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->branchauxmindepth <= SCIPgetDepth(scip);
}

/** gets weight of variable when splitting violation score onto several variables in an expression */
static
SCIP_Real getViolSplitWeight(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expr constraint handler */
   SCIP_VAR*             var,                /**< variable */
   SCIP_SOL*             sol                 /**< current solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   switch( conshdlrdata->branchviolsplit )
   {
      case 'u' :  /* uniform: everyone gets the same score */
         return 1.0;

      case 'm' :  /* midness of solution: 0.5 if in middle of domain, 0.05 if close to lower or upper bound */
      {
         SCIP_Real weight;
         weight = MIN(SCIPgetSolVal(scip, sol, var) - SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var) - SCIPgetSolVal(scip, sol, var)) / (SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var));
         return MAX(0.05, weight);
      }

      case 'd' :  /* domain width */
         return SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);

      case 'l' :  /* logarithmic domain width: log-scale if width is below 0.1 or above 10, otherwise actual width */
      {
         SCIP_Real width = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);
         assert(width > 0.0);
         if( width > 10.0 )
            return 10.0*log10(width);
         if( width < 0.1 )
            return 0.1/(-log10(width));
         return width;
      }

      default :
         SCIPerrorMessage("invalid value for parameter constraints/expr/branching/violsplit");
         SCIPABORT();
         return SCIP_INVALID;
   }
}

/** adds violation-branching score to a set of expressions, thereby distributing the score
 *
 * Each expression must either be a variable expression or have an aux-variable.
 *
 * If unbounded variables are present, each unbounded var gets an even score.
 * If no unbounded variables, then parameter constraints/nonlinear/branching/violsplit decides weight for each var.
 */
static
void addExprsViolScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           exprs,              /**< expressions where to add branching score */
   int                   nexprs,             /**< number of expressions */
   SCIP_Real             violscore,          /**< violation-branching score to add to expression */
   SCIP_SOL*             sol,                /**< current solution */
   SCIP_Bool*            success             /**< buffer to store whether at least one violscore was added */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR* var;
   SCIP_Real weight;
   SCIP_Real weightsum = 0.0; /* sum of weights over all candidates with bounded domain */
   int nunbounded = 0;  /* number of candidates with unbounded domain */
   int i;

   assert(exprs != NULL);
   assert(nexprs > 0);
   assert(success != NULL);

   if( nexprs == 1 )
   {
      SCIPaddExprViolScoreNonlinear(scip, exprs[0], violscore);
      SCIPdebugMsg(scip, "add score %g to <%s>[%g,%g]\n", violscore,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(exprs[0])), SCIPvarGetLbLocal(SCIPgetExprAuxVarNonlinear(exprs[0])), SCIPvarGetUbLocal(SCIPgetExprAuxVarNonlinear(exprs[0])));
      *success = TRUE;
      return;
   }

   conshdlr = SCIPexprGetOwnerData(exprs[0])->conshdlr;

   for( i = 0; i < nexprs; ++i )
   {
      var = SCIPgetExprAuxVarNonlinear(exprs[i]);
      assert(var != NULL);

      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
         ++nunbounded;
      else if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         weightsum += getViolSplitWeight(scip, conshdlr, var, sol);
   }

   *success = FALSE;
   for( i = 0; i < nexprs; ++i )
   {
      var = SCIPgetExprAuxVarNonlinear(exprs[i]);
      assert(var != NULL);

      if( nunbounded > 0 )
      {
         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
         {
            SCIPaddExprViolScoreNonlinear(scip, exprs[i], violscore / nunbounded);
            SCIPdebugMsg(scip, "add score %g (%g%% of %g) to <%s>[%g,%g]\n", violscore / nunbounded,
               100.0/nunbounded, violscore,
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
            *success = TRUE;
         }
      }
      else if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
      {
         assert(weightsum > 0.0);

         weight = getViolSplitWeight(scip, conshdlr, var, sol);
         SCIPaddExprViolScoreNonlinear(scip, exprs[i], violscore * weight / weightsum);
         SCIPdebugMsg(scip, "add score %g (%g%% of %g) to <%s>[%g,%g]\n", violscore * weight / weightsum,
            100*weight / weightsum, violscore,
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         *success = TRUE;
      }
      else
      {
         SCIPdebugMsg(scip, "skip score for fixed variable <%s>[%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      }
   }
}

/** adds violation-branching score to children of expression for given auxiliary variables
 *
 * Iterates over the successors of `expr` to find expressions that are associated with one of the given auxiliary variables.
 * Adds violation-branching scores to all found exprs by means of SCIPaddExprsViolScoreNonlinear().
 *
 * @note This method may modify the given auxvars array by means of sorting.
 */
static
SCIP_RETCODE addExprViolScoresAuxVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression where to start searching */
   SCIP_Real             violscore,          /**< violation score to add to expression */
   SCIP_VAR**            auxvars,            /**< auxiliary variables for which to find expression */
   int                   nauxvars,           /**< number of auxiliary variables */
   SCIP_SOL*             sol,                /**< current solution (NULL for the LP solution) */
   SCIP_Bool*            success             /**< buffer to store whether at least one violscore was added */
   )
{
   SCIP_EXPRITER* it;
   SCIP_VAR* auxvar;
   SCIP_EXPR** exprs;
   int nexprs;
   int pos;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(auxvars != NULL);
   assert(success != NULL);

   /* sort variables to make lookup below faster */
   SCIPsortPtr((void**)auxvars, SCIPvarComp, nauxvars);

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_BFS, FALSE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nauxvars) );
   nexprs = 0;

   for( expr = SCIPexpriterGetNext(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      auxvar = SCIPgetExprAuxVarNonlinear(expr);
      if( auxvar == NULL )
         continue;

      /* if auxvar of expr is contained in auxvars array, add branching score to expr */
      if( SCIPsortedvecFindPtr((void**)auxvars, SCIPvarComp, auxvar, nauxvars, &pos) )
      {
         assert(auxvars[pos] == auxvar);

         SCIPdebugMsg(scip, "adding branchingscore for expr %p with auxvar <%s>\n", (void*)expr, SCIPvarGetName(auxvar));
         exprs[nexprs++] = expr;

         if( nexprs == nauxvars )
            break;
      }
   }

   SCIPfreeExpriter(&it);

   if( nexprs > 0 )
   {
      SCIP_CALL( SCIPaddExprsViolScoreNonlinear(scip, exprs, nexprs, violscore, sol, success) );
   }
   else
      *success = FALSE;

   SCIPfreeBufferArray(scip, &exprs);

   return SCIP_OKAY;
}

/** registers all unfixed variables in violated constraints as branching candidates */
static
SCIP_RETCODE registerBranchingCandidatesAllUnfixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int c;
   int i;

   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nnotify != NULL);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* consider only violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      /* register all variables that have not been fixed yet */
      assert(consdata->varexprs != NULL);
      for( i = 0; i < consdata->nvarexprs; ++i )
      {
         var = SCIPgetVarExprVar(consdata->varexprs[i]);
         assert(var != NULL);

         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, getConsAbsViolation(conss[c]), SCIP_INVALID) );
            ++(*nnotify);
         }
      }
   }

   return SCIP_OKAY;
}

/** registers all variables in violated constraints with branching scores as external branching candidates */
static
SCIP_RETCODE registerBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            success             /**< buffer to store whether at least one branching candidate was added */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it = NULL;
   int c;

   assert(conshdlr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( branchAuxNonlinear(scip, conshdlr) )
   {
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
      SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   }

   /* register external branching candidates */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->varexprs != NULL);

      /* consider only violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      if( !branchAuxNonlinear(scip, conshdlr) )
      {
         int i;

         /* if not branching on auxvars, then violation-branching scores will have been added to original variables
          * only, so we can loop over variable expressions
          */
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_Real violscore;
            SCIP_Real lb;
            SCIP_Real ub;
            SCIP_VAR* var;

            violscore = SCIPgetExprViolScoreNonlinear(consdata->varexprs[i]);

            /* skip variable expressions that do not have a violation score */
            if( violscore == 0.0 )
               continue;

            var = SCIPgetVarExprVar(consdata->varexprs[i]);
            assert(var != NULL);

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);

            /* consider variable for branching if it has not been fixed yet */
            if( !SCIPisEQ(scip, lb, ub) )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " add variable <%s>[%g,%g] as extern branching candidate with score %g\n", SCIPvarGetName(var), lb, ub, violscore); )
               SCIP_CALL( SCIPaddExternBranchCand(scip, var, violscore, SCIP_INVALID) );
               *success = TRUE;
            }
            else
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
            }

            /* invalidate violscore-tag, so that we do not register variables that appear in multiple constraints
             * several times as external branching candidate, see SCIPgetExprViolScoreNonlinear()
             */
            SCIPexprGetOwnerData(consdata->varexprs[i])->violscoretag = 0;
         }
      }
      else
      {
         SCIP_EXPR* expr;
         SCIP_VAR* var;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real violscore;

         for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
         {
            violscore = SCIPgetExprViolScoreNonlinear(expr);
            if( violscore == 0.0 )
               continue;

            /* if some nlhdlr added a branching score for this expression, then it considered this expression as a
             * variable, so this expression should either be an original variable or have an auxiliary variable
             */
            var = SCIPgetExprAuxVarNonlinear(expr);
            assert(var != NULL);

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);

            /* consider variable for branching if it has not been fixed yet */
            if( !SCIPisEQ(scip, lb, ub) )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " add variable <%s>[%g,%g] as extern branching candidate with score %g\n", SCIPvarGetName(var), lb, ub, violscore); )

               SCIP_CALL( SCIPaddExternBranchCand(scip, var, violscore, SCIP_INVALID) );
               *success = TRUE;
            }
            else
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
            }
         }
      }
   }

   if( it != NULL )
      SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** collect branching candidates from violated constraints
 *
 * Fills array with expressions that serve as branching candidates.
 * Collects those expressions that have a branching score assigned and stores the score in the auxviol field of the
 * branching candidate.
 *
 * If branching on aux-variables is allowed, then iterate through expressions of violated constraints, otherwise iterate
 * through variable-expressions only.
 */
static
SCIP_RETCODE collectBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Real             maxrelconsviol,     /**< maximal scaled constraint violation */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_Longint          soltag,             /**< tag of solution */
   BRANCHCAND*           cands,              /**< array where to store candidates, must be at least SCIPgetNVars() long */
   int*                  ncands              /**< number of candidates found */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it = NULL;
   int c;
   int attempt;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cands != NULL);
   assert(ncands != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( branchAuxNonlinear(scip, conshdlr) )
   {
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
      SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   }

   *ncands = 0;
   for( attempt = 0; attempt < 2; ++attempt )
   {
      /* collect branching candidates from violated constraints
       * in the first attempt, consider only constraints with large violation
       * in the second attempt, consider all remaining violated constraints
       */
      for( c = 0; c < nconss; ++c )
      {
         SCIP_Real consviol;

         assert(conss != NULL && conss[c] != NULL);

         /* consider only violated constraints */
         if( !isConsViolated(scip, conss[c]) )
            continue;

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         assert(consdata->varexprs != NULL);

         SCIP_CALL( getConsRelViolation(scip, conss[c], &consviol, sol, soltag) );

         if( attempt == 0 && consviol < conshdlrdata->branchhighviolfactor * maxrelconsviol )
            continue;
         else if( attempt == 1 && consviol >= conshdlrdata->branchhighviolfactor * maxrelconsviol )
            continue;

         if( !branchAuxNonlinear(scip, conshdlr) )
         {
            int i;

            /* if not branching on auxvars, then violation-branching scores will be available for original variables
             * only, so we can loop over variable expressions
             * unfortunately, we don't know anymore which constraint contributed the violation-branching score to the
             * variable, therefore we invalidate the score of a variable after processing it.
             */
            for( i = 0; i < consdata->nvarexprs; ++i )
            {
               SCIP_Real lb;
               SCIP_Real ub;

               /* skip variable expressions that do not have a valid violation score */
               if( conshdlrdata->enforound != SCIPexprGetOwnerData(consdata->varexprs[i])->violscoretag )
                  continue;

               var = SCIPgetVarExprVar(consdata->varexprs[i]);
               assert(var != NULL);

               lb = SCIPvarGetLbLocal(var);
               ub = SCIPvarGetUbLocal(var);

               /* skip already fixed variable */
               if( SCIPisEQ(scip, lb, ub) )
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
                  continue;
               }

               assert(*ncands + 1 < SCIPgetNVars(scip));
               cands[*ncands].expr = consdata->varexprs[i];
               cands[*ncands].auxviol = SCIPgetExprViolScoreNonlinear(consdata->varexprs[i]);
               ++(*ncands);

               /* invalidate violscore-tag, so that we do not register variables that appear in multiple constraints
                * several times as external branching candidate */
               SCIPexprGetOwnerData(consdata->varexprs[i])->violscoretag = 0;
            }
         }
         else
         {
            SCIP_EXPR* expr;
            SCIP_Real lb;
            SCIP_Real ub;

            for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
            {
               if( SCIPexprGetOwnerData(expr)->violscoretag != conshdlrdata->enforound )
                  continue;

               /* if some nlhdlr added a branching score for this expression, then it considered this expression as
                * variables, so this expression should either be an original variable or have an auxiliary variable
                */
               var = SCIPgetExprAuxVarNonlinear(expr);
               assert(var != NULL);

               lb = SCIPvarGetLbLocal(var);
               ub = SCIPvarGetUbLocal(var);

               /* skip already fixed variable */
               if( SCIPisEQ(scip, lb, ub) )
               {
                  ENFOLOG( SCIPinfoMessage(scip, enfologfile, " skip fixed variable <%s>[%.15g,%.15g]\n", SCIPvarGetName(var), lb, ub); )
                  continue;
               }

               assert(*ncands + 1 < SCIPgetNVars(scip));
               cands[*ncands].expr = expr;
               cands[*ncands].auxviol = SCIPgetExprViolScoreNonlinear(expr);
               ++(*ncands);
            }
         }
      }

      /* if we have branching candidates, then we don't need another attempt */
      if( *ncands > 0 )
         break;
   }

   if( it != NULL )
      SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** computes a branching score for a variable that reflects how important branching on this variable would be for
 * improving the dual bound from the LP relaxation
 *
 * Assume the Lagrangian for the current LP is something of the form
 *   L(x,z,lambda) = c'x + sum_i lambda_i (a_i'x - z_i + b_i) + ...
 * where x are the original variables, z the auxiliary variables,
 * and a_i'x - z_i + b_i <= 0 are the rows of the LP.
 *
 * Assume that a_i'x + b_i <= z_i was derived from some nonlinear constraint f(x) <= z and drop index i.
 * If we could have used not only an estimator, but the actual function f(x), then this would
 * have contributed lambda*(f(x) - z) to the Lagrangian function (though the value of z would be different).
 * Using a lot of handwaving, we claim that
 *   lambda_i * (f(x) - a_i'x + b_i)
 * is a value that can be used to quantity how much improving the estimator a'x + b <= z could change the dual bound.
 * If an estimator depended on local bounds, then it could be improved by branching.
 * We use row-is-local as proxy for estimator-depending-on-lower-bounds.
 *
 * To score a variable, we then sum the values lambda_i * (f(x) - a_i'x + b_i) for all rows in which the variable appears.
 * To scale, we divide by the LP objective value (if >1).
 *
 * TODO if we branch only on original variables, we neglect here estimators that are build on auxiliary variables;
 *     these are affected by the bounds on original variables indirectly (through forward-propagation)
 *
 * TODO if we branch also on auxiliary variables, then separating z from the x-variables in the row a'x+b <= z should happen;
 *     in effect, we should go from the row to the expression for which it was generated and consider only variables that
 *     would also be branching candidates
 */
static
SCIP_Real getDualBranchscore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_COL* col;
   SCIP_ROW** rows;
   int nrows;
   int r;
   SCIP_Real dualscore;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(var != NULL);

   /* if LP not solved, then the dual branching score is not available */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return 0.0;

   /* if var is not in the LP, then the dual branching score is not available */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return 0.0;

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
      return 0.0;

   nrows = SCIPcolGetNLPNonz(col);  /* TODO there is a big warning on when not to use this method; is the check for SCIPcolIsInLP sufficient? */
   rows = SCIPcolGetRows(col);

   /* SCIPinfoMessage(scip, enfologfile, " dualscoring <%s>\n", SCIPvarGetName(var)); */

   /* aggregate duals from all rows from consexpr with non-zero dual
    * TODO: this is a quick-and-dirty implementation, and not used by default
    *   in the long run, this should be either removed or replaced by a proper implementation
    */
   dualscore = 0.0;
   for( r = 0; r < nrows; ++r )
   {
      SCIP_Real estimategap;
      const char* estimategapstr;

      /* rows from cuts that may be replaced by tighter ones after branching are the interesting ones
       * these would typically be local, unless they are created at the root node
       * so not check for local now, but trust that estimators that do not improve after branching will have an estimategap of 0
      if( !SCIProwIsLocal(rows[r]) )
         continue;
       */
      if( SCIProwGetOriginConshdlr(rows[r]) != conshdlr )
         continue;
      if( SCIPisZero(scip, SCIProwGetDualsol(rows[r])) )
         continue;

      estimategapstr = strstr(SCIProwGetName(rows[r]), "_estimategap=");
      if( estimategapstr == NULL ) /* gap not stored, maybe because it was 0 */
         continue;
      estimategap = atof(estimategapstr + 13);
      assert(estimategap >= 0.0);
      if( !SCIPisFinite(estimategap) || SCIPisHugeValue(scip, estimategap) )
         estimategap = SCIPgetHugeValue(scip);

      /* SCIPinfoMessage(scip, enfologfile, "  row <%s> contributes %g*|%g|: ", SCIProwGetName(rows[r]), estimategap, SCIProwGetDualsol(rows[r]));
      SCIP_CALL( SCIPprintRow(scip, rows[r], enfologfile) ); */

      dualscore += estimategap * REALABS(SCIProwGetDualsol(rows[r]));
   }

   /* divide by optimal value of LP for scaling */
   dualscore /= MAX(1.0, REALABS(SCIPgetLPObjval(scip)));

   return dualscore;
}

/** computes branching scores (including weighted score) for a set of candidates
 *
 * For each candidate in the array, compute and store the various branching scores (violation, pseudo-costs, vartype, domainwidth).
 * For pseudo-costs, it's possible that the score is not available, in which case cands[c].pscost will be set to SCIP_INVALID.
 *
 * For each score, compute the maximum over all candidates.
 *
 * Then compute for each candidate a "weighted" score using the weights as specified by parameters
 * and the scores as previously computed, but scale each score to be in [0,1], i.e., divide each score by the maximum
 * score of all candidates.
 * Further divide by the sum of all weights where a score was available (even if the score was 0).
 *
 * For example:
 * - Let variable x have violation-score 10.0 and pseudo-cost-score 5.0.
 * - Let variable y have violation-score 12.0 but no pseudo-cost-score (because it hasn't yet been branched on sufficiently often).
 * - Assuming violation is weighted by 2.0 and pseudo-costs are weighted by 3.0.
 * - Then the weighted scores for x will be (2.0 * 10.0/12.0 + 3.0 * 5.0/5.0) / (2.0 + 3.0) = 0.9333.
 *   The weighted score for y will be (2.0 * 12.0/12.0) / 2.0 = 1.0.
 */
static
void scoreBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BRANCHCAND*           cands,              /**< branching candidates */
   int                   ncands,             /**< number of candidates */
   SCIP_SOL*             sol                 /**< solution to enforce (NULL for the LP solution) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   BRANCHCAND maxscore;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cands != NULL);
   assert(ncands > 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* initialize counts to 0 */
   memset(&maxscore, 0, sizeof(BRANCHCAND));

   for( c = 0; c < ncands; ++c )
   {
      if( conshdlrdata->branchviolweight > 0.0 )
      {
         /* cands[c].auxviol was set in collectBranchingCandidates, so only update maxscore here */
         maxscore.auxviol = MAX(maxscore.auxviol, cands[c].auxviol);
      }

      if( conshdlrdata->branchdomainweight > 0.0 )
      {
         SCIP_Real domainwidth;
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         /* get domain width, taking infinity at 1e20 on purpose */
         domainwidth = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);

         /* domain-score is going to be log(2*infinity / domainwidth) if domain width >= 1
          * and log(2 * infinity *  MAX(epsilon, domainwidth)) for domain width < 1
          * the idea is to penalize very large and very small domains
          */
         if( domainwidth >= 1.0 )
            cands[c].domain = log10(2 * SCIPinfinity(scip) / domainwidth);
         else
            cands[c].domain = log10(2 * SCIPinfinity(scip) * MAX(SCIPepsilon(scip), domainwidth));

         maxscore.domain = MAX(cands[c].domain, maxscore.domain);
      }
      else
         cands[c].domain = 0.0;

      if( conshdlrdata->branchdualweight > 0.0 )
      {
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         cands[c].dual = getDualBranchscore(scip, conshdlr, var);
         maxscore.dual = MAX(cands[c].dual, maxscore.dual);
      }

      if( conshdlrdata->branchpscostweight > 0.0 && SCIPgetNObjVars(scip) > 0 )
      {
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
            cands[c].pscost = SCIP_INVALID;
         else
         {
            SCIP_Real brpoint;
            SCIP_Real pscostdown;
            SCIP_Real pscostup;
            char strategy;

            /* decide how to compute pseudo-cost scores
             * this should be consistent with the way how pseudo-costs are updated in the core, which is decided by
             * branching/lpgainnormalize for continuous variables and move in LP-value for non-continuous variables
             */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               strategy = conshdlrdata->branchpscostupdatestrategy;
            else
               strategy = 'l';

            brpoint = SCIPgetBranchingPoint(scip, var, SCIP_INVALID);

            /* branch_relpscost deems pscosts as reliable, if the pseudo-count is at least something between 1 and 4
             * or it uses some statistical tests involving SCIPisVarPscostRelerrorReliable
             * For here, I use a simple #counts >= branchpscostreliable.
             * TODO use SCIPgetVarPseudocostCount() instead?
             */
            if( SCIPgetVarPseudocostCountCurrentRun(scip, var, SCIP_BRANCHDIR_DOWNWARDS) >= conshdlrdata->branchpscostreliable )
            {
               switch( strategy )
               {
                  case 's' :
                     pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPvarGetUbLocal(var) - SCIPadjustedVarLb(scip, var, brpoint)));
                     break;
                  case 'd' :
                     pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPadjustedVarUb(scip, var, brpoint) - SCIPvarGetLbLocal(var)));
                     break;
                  case 'l' :
                     if( SCIPisInfinity(scip, SCIPgetSolVal(scip, sol, var)) )
                        pscostdown = SCIP_INVALID;
                     else if( SCIPgetSolVal(scip, sol, var) <= SCIPadjustedVarUb(scip, var, brpoint) )
                        pscostdown = SCIPgetVarPseudocostVal(scip, var, 0.0);
                     else
                        pscostdown = SCIPgetVarPseudocostVal(scip, var, -(SCIPgetSolVal(scip, NULL, var) - SCIPadjustedVarUb(scip, var, brpoint)));
                     break;
                  default :
                     SCIPerrorMessage("pscost update strategy %c unknown\n", strategy);
                     pscostdown = SCIP_INVALID;
               }
            }
            else
               pscostdown = SCIP_INVALID;

            if( SCIPgetVarPseudocostCountCurrentRun(scip, var, SCIP_BRANCHDIR_UPWARDS) >= conshdlrdata->branchpscostreliable )
            {
               switch( strategy )
               {
                  case 's' :
                     pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPadjustedVarUb(scip, var, brpoint) - SCIPvarGetLbLocal(var));
                     break;
                  case 'd' :
                     pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPvarGetUbLocal(var) - SCIPadjustedVarLb(scip, var, brpoint));
                     break;
                  case 'l' :
                     if( SCIPisInfinity(scip, -SCIPgetSolVal(scip, sol, var)) )
                        pscostup = SCIP_INVALID;
                     else if( SCIPgetSolVal(scip, NULL, var) >= SCIPadjustedVarLb(scip, var, brpoint) )
                        pscostup = SCIPgetVarPseudocostVal(scip, var, 0.0);
                     else
                        pscostup = SCIPgetVarPseudocostVal(scip, var, SCIPadjustedVarLb(scip, var, brpoint) - SCIPgetSolVal(scip, NULL, var) );
                     break;
                  default :
                     SCIPerrorMessage("pscost update strategy %c unknown\n", strategy);
                     pscostup = SCIP_INVALID;
               }
            }
            else
               pscostup = SCIP_INVALID;

            /* TODO if both are valid, we get pscostdown*pscostup, but does this compare well with vars were only pscostdown or pscostup is used?
             * maybe we should use (pscostdown+pscostup)/2 or sqrt(pscostdown*pscostup) ?
             */
            if( pscostdown == SCIP_INVALID && pscostup == SCIP_INVALID )
               cands[c].pscost = SCIP_INVALID;
            else if( pscostdown == SCIP_INVALID )
               cands[c].pscost = pscostup;
            else if( pscostup == SCIP_INVALID )
               cands[c].pscost = pscostdown;
            else
               cands[c].pscost = SCIPgetBranchScore(scip, NULL, pscostdown, pscostup);  /* pass NULL for var to avoid multiplication with branch-factor */
         }

         if( cands[c].pscost != SCIP_INVALID )
            maxscore.pscost = MAX(cands[c].pscost, maxscore.pscost);
      }

      if( conshdlrdata->branchvartypeweight > 0.0 )
      {
         SCIP_VAR* var;

         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         assert(var != NULL);

         switch( SCIPvarGetType(var) )
         {
            case SCIP_VARTYPE_BINARY :
               cands[c].vartype = 1.0;
               break;
            case SCIP_VARTYPE_INTEGER :
               cands[c].vartype = 0.1;
               break;
            case SCIP_VARTYPE_IMPLINT :
               cands[c].vartype = 0.01;
               break;
            case SCIP_VARTYPE_CONTINUOUS :
            default:
               cands[c].vartype = 0.0;
         }
         maxscore.vartype = MAX(cands[c].vartype, maxscore.vartype);
      }
   }

   /* now compute a weighted score for each candidate from the single scores
    * the single scores are scaled to be in [0,1] for this
    */
   for( c = 0; c < ncands; ++c )
   {
      SCIP_Real weightsum;

      ENFOLOG(
         SCIP_VAR* var;
         var = SCIPgetExprAuxVarNonlinear(cands[c].expr);
         SCIPinfoMessage(scip, enfologfile, " scoring <%8s>[%7.1g,%7.1g]:(", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         )

      cands[c].weighted = 0.0;
      weightsum = 0.0;

      if( maxscore.auxviol > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchviolweight * cands[c].auxviol / maxscore.auxviol;
         weightsum += conshdlrdata->branchviolweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(viol)", conshdlrdata->branchviolweight, cands[c].auxviol / maxscore.auxviol); )
      }

      if( maxscore.domain > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchdomainweight * cands[c].domain / maxscore.domain;
         weightsum += conshdlrdata->branchdomainweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(domain)", conshdlrdata->branchdomainweight, cands[c].domain / maxscore.domain); )
      }

      if( maxscore.dual > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchdualweight * cands[c].dual / maxscore.dual;
         weightsum += conshdlrdata->branchdualweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(dual)", conshdlrdata->branchdualweight, cands[c].dual / maxscore.dual); )
      }

      if( maxscore.pscost > 0.0 )
      {
         /* use pseudo-costs only if available */
         if( cands[c].pscost != SCIP_INVALID )
         {
            cands[c].weighted += conshdlrdata->branchpscostweight * cands[c].pscost / maxscore.pscost;
            weightsum += conshdlrdata->branchpscostweight;

            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%7.2g(pscost)", conshdlrdata->branchpscostweight, cands[c].pscost / maxscore.pscost); )
         }
         else
         {
            /* do not add pscostscore, if not available, also do not add into weightsum */
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " +0.0*    n/a(pscost)"); )
         }
      }

      if( maxscore.vartype > 0.0 )
      {
         cands[c].weighted += conshdlrdata->branchvartypeweight * cands[c].vartype / maxscore.vartype;
         weightsum += conshdlrdata->branchvartypeweight;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %+g*%6.2g(vartype)", conshdlrdata->branchvartypeweight, cands[c].vartype / maxscore.vartype); )
      }
      assert(weightsum > 0.0);  /* we should have got at least one valid score */
      cands[c].weighted /= weightsum;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " ) / %g = %g\n", weightsum, cands[c].weighted); )
   }
}

/** compare two branching candidates by their weighted score
 *
 * if weighted score is equal, use variable index of (aux)var
 */
static
SCIP_DECL_SORTINDCOMP(branchcandCompare)
{
   BRANCHCAND* cands = (BRANCHCAND*)dataptr;

   if( cands[ind1].weighted != cands[ind2].weighted )
      return cands[ind1].weighted < cands[ind2].weighted ? -1 : 1;
   else
      return SCIPvarGetIndex(SCIPgetExprAuxVarNonlinear(cands[ind1].expr)) - SCIPvarGetIndex(SCIPgetExprAuxVarNonlinear(cands[ind2].expr));
}

/** do branching or register branching candidates */
static
SCIP_RETCODE branching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Real             maxrelconsviol,     /**< maximal scaled constraint violation */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_Longint          soltag,             /**< tag of solution */
   SCIP_RESULT*          result              /**< pointer to store the result of branching */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   BRANCHCAND* cands;
   int ncands;
   SCIP_VAR* var;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;

   assert(conshdlr != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->branchexternal )
   {
      /* just register branching candidates as external */
      SCIP_Bool success;

      SCIP_CALL( registerBranchingCandidates(scip, conshdlr, conss, nconss, &success) );
      if( success )
         *result = SCIP_INFEASIBLE;

      return SCIP_OKAY;
   }

   /* collect branching candidates and their auxviol-score */
   SCIP_CALL( SCIPallocBufferArray(scip, &cands, SCIPgetNVars(scip)) );
   SCIP_CALL( collectBranchingCandidates(scip, conshdlr, conss, nconss, maxrelconsviol, sol, soltag, cands, &ncands) );

   /* if no unfixed branching candidate in all violated constraint, then it's probably numerics that prevented us to separate or decide a cutoff
    * we will return here and let the fallbacks in consEnfo() decide how to proceed
    */
   if( ncands == 0 )
      goto TERMINATE;

   if( ncands > 1 )
   {
      /* if there are more than one candidate, then compute scores and select */
      int* perm;
      int c;
      int left;
      int right;
      SCIP_Real threshold;

      /* compute additional scores on branching candidates and weighted score */
      scoreBranchingCandidates(scip, conshdlr, cands, ncands, sol);

      /* sort candidates by weighted score */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, ncands) );
      SCIPsortDown(perm, branchcandCompare, (void*)cands, ncands);

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %d branching candidates <%s>(%g)...<%s>(%g)\n", ncands,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[0]].expr)), cands[perm[0]].weighted,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[ncands - 1]].expr)), cands[perm[ncands - 1]].weighted); )

      /* binary search to find first low-scored (score below branchhighscorefactor * maximal-score)  candidate */
      left = 0;
      right = ncands - 1;
      threshold = conshdlrdata->branchhighscorefactor * cands[perm[0]].weighted;
      while( left < right )
      {
         int mid = (left + right) / 2;
         if( cands[perm[mid]].weighted >= threshold )
            left = mid + 1;
         else
            right = mid;
      }
      assert(left <= ncands);

      if( left < ncands )
      {
         if( cands[perm[left]].weighted >= threshold )
         {
            assert(left + 1 == ncands || cands[perm[left + 1]].weighted < threshold);
            ncands = left + 1;
         }
         else
         {
            assert(cands[perm[left]].weighted < threshold);
            ncands = left;
         }
      }
      assert(ncands > 0);

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " %d branching candidates <%s>(%g)...<%s>(%g) after removing low scores\n", ncands,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[0]].expr)), cands[perm[0]].weighted,
         SCIPvarGetName(SCIPgetExprAuxVarNonlinear(cands[perm[ncands - 1]].expr)), cands[perm[ncands - 1]].weighted); )

      if( ncands > 1 )
      {
         /* choose at random from candidates 0..ncands-1 */
         if( conshdlrdata->branchrandnumgen == NULL )
         {
            SCIP_CALL( SCIPcreateRandom(scip, &conshdlrdata->branchrandnumgen, BRANCH_RANDNUMINITSEED, TRUE) );
         }
         c = SCIPrandomGetInt(conshdlrdata->branchrandnumgen, 0, ncands - 1);
         var = SCIPgetExprAuxVarNonlinear(cands[perm[c]].expr);
      }
      else
         var = SCIPgetExprAuxVarNonlinear(cands[perm[0]].expr);

      SCIPfreeBufferArray(scip, &perm);
   }
   else
   {
      var = SCIPgetExprAuxVarNonlinear(cands[0].expr);
   }
   assert(var != NULL);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " branching on variable <%s>[%g,%g]\n", SCIPvarGetName(var),
            SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)); )

   SCIP_CALL( SCIPbranchVarVal(scip, var, SCIPgetBranchingPoint(scip, var, SCIP_INVALID), &downchild, &eqchild,
            &upchild) );
   if( downchild != NULL || eqchild != NULL || upchild != NULL )
      *result = SCIP_BRANCHED;
   else
      /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
      *result = SCIP_REDUCEDDOM;

 TERMINATE:
   SCIPfreeBufferArray(scip, &cands);

   return SCIP_OKAY;
}

/** call enforcement or estimate callback of nonlinear handler
 *
 * Calls the enforcement callback, if available.
 * Otherwise, calls the estimate callback, if available, and constructs a cut from the estimator.
 *
 * If cut is weak, but estimator is not tight, tries to add branching candidates.
 */
static
SCIP_RETCODE enforceExprNlhdlr(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler data of expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_Bool             separated,          /**< whether another nonlinear handler already added a cut for this expression */
   SCIP_Bool             allowweakcuts,      /**< whether we allow for weak cuts */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement (and not just separation) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   assert(result != NULL);

   /* call enforcement callback of the nlhdlr */
   SCIP_CALL( SCIPnlhdlrEnfo(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
            allowweakcuts, separated, inenforcement, result) );

   /* if it was not running (e.g., because it was not available) or did not find anything, then try with estimator callback */
   if( *result != SCIP_DIDNOTRUN && *result != SCIP_DIDNOTFIND )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    sepa of nlhdlr %s succeeded with result %d\n",
               SCIPnlhdlrGetName(nlhdlr), *result); )
      return SCIP_OKAY;
   }
   else
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    sepa of nlhdlr <%s> did not succeed with result %d\n", SCIPnlhdlrGetName(nlhdlr), *result); )
   }

   *result = SCIP_DIDNOTFIND;

   /* now call the estimator callback of the nlhdlr */
   if( SCIPnlhdlrHasEstimate(nlhdlr) )
   {
      SCIP_VAR* auxvar;
      SCIP_Bool sepasuccess = FALSE;
      SCIP_Bool branchscoresuccess = FALSE;
      SCIP_PTRARRAY* rowpreps;
      int minidx;
      int maxidx;
      int r;
      SCIP_ROWPREP* rowprep;

      SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );

      auxvar = SCIPgetExprAuxVarNonlinear(expr);
      assert(auxvar != NULL);

      SCIP_CALL( SCIPnlhdlrEstimate(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, overestimate,
               SCIPgetSolVal(scip, sol, auxvar), inenforcement, rowpreps, &sepasuccess, &branchscoresuccess) );

      minidx = SCIPgetPtrarrayMinIdx(scip, rowpreps);
      maxidx = SCIPgetPtrarrayMaxIdx(scip, rowpreps);

      assert((sepasuccess && minidx <= maxidx) || (!sepasuccess && minidx > maxidx));

      if( !sepasuccess )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s failed\n",
                  SCIPnlhdlrGetName(nlhdlr)); )
      }

      for( r = minidx; r <= maxidx; ++r )
      {
         rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, r);

         assert(rowprep != NULL);
         assert(SCIProwprepGetSidetype(rowprep) == (overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT));

         /* complete estimator to cut */
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0) );

         /* add the cut and/or branching scores */
         SCIP_CALL( SCIPprocessRowprepNonlinear(scip, nlhdlr, cons, expr, rowprep, overestimate, auxvar,
               auxvalue, allowweakcuts, branchscoresuccess, inenforcement, sol, result) );

         SCIPfreeRowprep(scip, &rowprep);
      }

      SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   }

   return SCIP_OKAY;
}

/** tries to enforce violation in an expression by separation, bound tightening, or finding a branching candidate
 *
 * if not inenforcement, then we should be called by consSepa(), and thus only try separation
 */
static
SCIP_RETCODE enforceExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Longint          soltag,             /**< tag of solution */
   SCIP_Bool             allowweakcuts,      /**< whether we allow weak cuts */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement (and not just separation) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_Real origviol;
   SCIP_Bool underestimate;
   SCIP_Bool overestimate;
   SCIP_Real auxviol;
   SCIP_Bool auxunderestimate;
   SCIP_Bool auxoverestimate;
   SCIP_RESULT hdlrresult;
   int e;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   assert(ownerdata->auxvar != NULL);  /* there must be a variable attached to the expression in order to construct a cut here */

   *result = SCIP_DIDNOTFIND;

   /* make sure that this expression has been evaluated */
   SCIP_CALL( SCIPevalExpr(scip, expr, sol, soltag) );

   /* decide whether under- or overestimate is required and get amount of violation */
   origviol = getExprAbsOrigViolation(scip, expr, sol, &underestimate, &overestimate);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* no sufficient violation w.r.t. the original variables -> skip expression */
   if( !overestimate && !underestimate )
   {
      return SCIP_OKAY;
   }

   /* check aux-violation w.r.t. each nonlinear handlers and try to enforce when there is a decent violation */
   for( e = 0; e < ownerdata->nenfos; ++e )
   {
      SCIP_NLHDLR* nlhdlr;

      /* skip nlhdlr that do not want to participate in any separation */
      if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 )
         continue;

      nlhdlr = ownerdata->enfos[e]->nlhdlr;
      assert(nlhdlr != NULL);

      /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
      SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata, &ownerdata->enfos[e]->auxvalue, sol) );
      ENFOLOG(
         SCIPinfoMessage(scip, enfologfile, "  expr ");
         SCIPprintExpr(scip, expr, enfologfile);
         SCIPinfoMessage(scip, enfologfile, " (%p): evalvalue %.15g auxvarvalue %.15g [%.15g,%.15g], nlhdlr <%s> " \
            "auxvalue: %.15g\n", (void*)expr, SCIPexprGetEvalValue(expr), SCIPgetSolVal(scip, sol, ownerdata->auxvar),
            SCIPexprGetActivity(expr).inf, SCIPexprGetActivity(expr).sup, SCIPnlhdlrGetName(nlhdlr), ownerdata->enfos[e]->auxvalue);
      )

      /* TODO if expr is root of constraint (consdata->expr == expr),
       * then compare auxvalue with constraint sides instead of auxvarvalue, as the former is what actually matters
       * that is, if auxvalue is good enough for the constraint to be satisfied, but when looking at evalvalue we see
       * the the constraint is violated, then some of the auxvars that nlhdlr uses is not having a good enough value,
       * so we should enforce in these auxiliaries first
       * if changing this here, we must also adapt analyzeViolation()
       */

      auxviol = getExprAbsAuxViolation(scip, expr, ownerdata->enfos[e]->auxvalue, sol, &auxunderestimate, &auxoverestimate);
      assert(auxviol >= 0.0);

      /* if aux-violation is much smaller than orig-violation, then better enforce further down in the expression first */
      if( !SCIPisInfinity(scip, auxviol) && auxviol < conshdlrdata->enfoauxviolfactor * origviol )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   skip enforce using nlhdlr <%s> for expr %p (%s) with " \
                  "auxviolation %g << origviolation %g under:%d over:%d\n", SCIPnlhdlrGetName(nlhdlr), (void*)expr,
                  SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), auxviol, origviol, underestimate, overestimate); )

         /* TODO should we do expr->lastenforced = conshdlrdata->enforound even though we haven't enforced, but only decided not to enforce? */
         continue;
      }

      /* if aux-violation is small (below feastol) and we look only for strong cuts, then it's unlikely to give a strong cut, so skip it */
      if( !allowweakcuts && auxviol < SCIPfeastol(scip) )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   skip enforce using nlhdlr <%s> for expr %p (%s) with tiny " \
                  "auxviolation %g under:%d over:%d\n", SCIPnlhdlrGetName(nlhdlr), (void*)expr, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), auxviol,
                  underestimate, overestimate); )

            /* TODO should we do expr->lastenforced = conshdlrdata->enforound even though we haven't enforced, but only decided not to enforce? */
         continue;
      }

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "   enforce using nlhdlr <%s> for expr %p (%s) with auxviolation " \
               "%g origviolation %g under:%d over:%d weak:%d\n", SCIPnlhdlrGetName(nlhdlr), (void*)expr, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)),
               auxviol, origviol, underestimate, overestimate, allowweakcuts); )

      /* if we want to overestimate and violation w.r.t. auxiliary variables is also present on this side and nlhdlr
       * wants to be called for separation on this side, then call separation of nlhdlr
       */
      if( overestimate && auxoverestimate && (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPAABOVE) != 0 )
      {
         /* call the separation or estimation callback of the nonlinear handler for overestimation */
         hdlrresult = SCIP_DIDNOTFIND;
         SCIP_CALL( enforceExprNlhdlr(scip, conshdlr, cons, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata, sol,
            ownerdata->enfos[e]->auxvalue, TRUE, *result == SCIP_SEPARATED, allowweakcuts, inenforcement, &hdlrresult) );

         if( hdlrresult == SCIP_CUTOFF )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    found a cutoff -> stop separation\n"); )
            *result = SCIP_CUTOFF;
            ownerdata->lastenforced = conshdlrdata->enforound;
            break;
         }

         if( hdlrresult == SCIP_SEPARATED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by cut\n", SCIPnlhdlrGetName(nlhdlr)); )
            *result = SCIP_SEPARATED;
            ownerdata->lastenforced = conshdlrdata->enforound;
            /* TODO or should we give other nlhdlr another chance? (also #3070) */
            break;
         }

         if( hdlrresult == SCIP_REDUCEDDOM )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by boundchange\n", SCIPnlhdlrGetName(nlhdlr)); )
            *result = SCIP_REDUCEDDOM;
            ownerdata->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_BRANCHED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> added branching candidate\n", SCIPnlhdlrGetName(nlhdlr)); )
            assert(inenforcement);

            /* separation and domain reduction takes precedence over branching */
            assert(*result == SCIP_DIDNOTFIND || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED);
            if( *result == SCIP_DIDNOTFIND )
               *result = SCIP_BRANCHED;
            ownerdata->lastenforced = conshdlrdata->enforound;
         }
      }

      /* if we want to underestimate and violation w.r.t. auxiliary variables is also present on this side and nlhdlr
       * wants to be called for separation on this side, then call separation of nlhdlr
       */
      if( underestimate && auxunderestimate && (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABELOW) != 0 )
      {
         /* call the separation or estimation callback of the nonlinear handler for underestimation */
         hdlrresult = SCIP_DIDNOTFIND;
         SCIP_CALL( enforceExprNlhdlr(scip, conshdlr, cons, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata, sol,
            ownerdata->enfos[e]->auxvalue, FALSE, *result == SCIP_SEPARATED, allowweakcuts, inenforcement, &hdlrresult) );

         if( hdlrresult == SCIP_CUTOFF )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    found a cutoff -> stop separation\n"); )
            *result = SCIP_CUTOFF;
            ownerdata->lastenforced = conshdlrdata->enforound;
            break;
         }

         if( hdlrresult == SCIP_SEPARATED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by cut\n", SCIPnlhdlrGetName(nlhdlr)); )
            *result = SCIP_SEPARATED;
            ownerdata->lastenforced = conshdlrdata->enforound;
            /* TODO or should we give other nlhdlr another chance? (also #3070) */
            break;
         }

         if( hdlrresult == SCIP_REDUCEDDOM )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> separating the current solution by boundchange\n", SCIPnlhdlrGetName(nlhdlr)); )
            *result = SCIP_REDUCEDDOM;
            ownerdata->lastenforced = conshdlrdata->enforound;
            /* TODO or should we always just stop here? */
         }

         if( hdlrresult == SCIP_BRANCHED )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    nlhdlr <%s> added branching candidate\n", SCIPnlhdlrGetName(nlhdlr)); )
            assert(inenforcement);

            /* separation takes precedence over branching */
            assert(*result == SCIP_DIDNOTFIND || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED);
            if( *result == SCIP_DIDNOTFIND )
               *result = SCIP_BRANCHED;
            ownerdata->lastenforced = conshdlrdata->enforound;
         }
      }
   }

   return SCIP_OKAY;
}

/** helper function to enforce a single constraint */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_Longint          soltag,             /**< tag of solution */
   SCIP_EXPRITER*        it,                 /**< expression iterator that we can just use here */
   SCIP_Bool             allowweakcuts,      /**< whether to allow weak cuts in this round */
   SCIP_Bool             inenforcement,      /**< whether to we are in enforcement, and not just separation */
   SCIP_RESULT*          result,             /**< pointer to update with result of the enforcing call */
   SCIP_Bool*            success             /**< buffer to store whether some enforcement took place */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* expr;

   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(it != NULL);
   assert(result != NULL);
   assert(success != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPexprGetOwnerData(consdata->expr)->nenfos >= 0);

   *success = FALSE;

   if( inenforcement && !consdata->ispropagated )
   {
      /* If there are boundchanges that haven't been propagated to activities yet, then do this now and update bounds of
       * auxiliary variables, since some nlhdlr/exprhdlr may look at auxvar bounds or activities
       * (TODO: nlhdlr tells us now whether they do and so we could skip).
       * For now, update bounds of auxiliary variables only if called from enforcement, since updating auxvar bounds in
       * separation doesn't seem to be right (it would be ok if the boundchange cuts off the current LP solution by a
       * nice amount, but if not, we may just add a boundchange that doesn't change the dual bound much and could
       * confuse the stalling check for how long to do separation).
       */
      SCIP_Bool infeasible;
      int ntightenings;

      SCIP_CALL( forwardPropExpr(scip, conshdlr, consdata->expr, inenforcement, &infeasible, &ntightenings) );
      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      /* if we tightened an auxvar bound, we better communicate that */
      if( ntightenings > 0 )
         *result = SCIP_REDUCEDDOM;
   }

   for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      SCIP_EXPR_OWNERDATA* ownerdata;
      SCIP_RESULT resultexpr;

      ownerdata = SCIPexprGetOwnerData(expr);
      assert(ownerdata != NULL);

      /* we can only enforce if there is an auxvar to compare with */
      if( ownerdata->auxvar == NULL )
         continue;

      assert(ownerdata->lastenforced <= conshdlrdata->enforound);
      if( ownerdata->lastenforced == conshdlrdata->enforound )
      {
         ENFOLOG(
            SCIPinfoMessage(scip, enfologfile, "  skip expr ");
            SCIPprintExpr(scip, expr, enfologfile);
            SCIPinfoMessage(scip, enfologfile, " as already enforced in this enforound\n");
         )
         *success = TRUE;
         continue;
      }

      SCIP_CALL( enforceExpr(scip, conshdlr, cons, expr, sol, soltag, allowweakcuts, inenforcement, &resultexpr) );

      /* if not enforced, then we must not have found a cutoff, cut, domain reduction, or branchscore */
      assert((ownerdata->lastenforced == conshdlrdata->enforound) == (resultexpr != SCIP_DIDNOTFIND));
      if( ownerdata->lastenforced == conshdlrdata->enforound )
         *success = TRUE;

      if( resultexpr == SCIP_CUTOFF )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if( resultexpr == SCIP_SEPARATED )
         *result = SCIP_SEPARATED;

      if( resultexpr == SCIP_REDUCEDDOM && *result != SCIP_SEPARATED )
         *result = SCIP_REDUCEDDOM;

      if( resultexpr == SCIP_BRANCHED && *result != SCIP_SEPARATED && *result != SCIP_REDUCEDDOM )
         *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** try to separate violated constraints and, if in enforcement, register branching scores
 *
 * Sets result to
 * - SCIP_DIDNOTFIND, if nothing of the below has been done
 * - SCIP_CUTOFF, if node can be cutoff,
 * - SCIP_SEPARATED, if a cut has been added,
 * - SCIP_REDUCEDDOM, if a domain reduction has been found,
 * - SCIP_BRANCHED, if branching has been done,
 * - SCIP_REDUCEDDOM, if a variable got fixed (in an attempt to branch on it),
 * - SCIP_INFEASIBLE, if external branching candidates were registered
 */
static
SCIP_RETCODE enforceConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_Longint          soltag,             /**< tag of solution */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, and not just separation */
   SCIP_Real             maxrelconsviol,     /**< largest scaled violation among all violated expr-constraints, only used if in enforcement */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRITER* it;
   SCIP_Bool consenforced;  /* whether any expression in constraint could be enforced */
   int c;

   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* increase tag to tell whether branching scores in expression belong to this sweep
    * and which expressions have already been enforced in this sweep
    * (we also want to distinguish sepa rounds, so this need to be here and not in consEnfo)
    */
   ++(conshdlrdata->enforound);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, TRUE) );

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      /* skip constraints that are not enabled or deleted */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      /* skip constraints that have separation disabled if we are only in separation */
      if( !inenforcement && !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;

      /* skip non-violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      ENFOLOG(
      {
         SCIP_CONSDATA* consdata;
         int i;
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         SCIPinfoMessage(scip, enfologfile, " constraint ");
         SCIP_CALL( SCIPprintCons(scip, conss[c], enfologfile) );
         SCIPinfoMessage(scip, enfologfile, "\n with viol %g and point\n", getConsAbsViolation(conss[c]));
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_VAR* var;
            var = SCIPgetVarExprVar(consdata->varexprs[i]);
            SCIPinfoMessage(scip, enfologfile, "  %-10s = %15g bounds: [%15g,%15g]\n", SCIPvarGetName(var),
                  SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         }
      })

      SCIP_CALL( enforceConstraint(scip, conshdlr, conss[c], sol, soltag, it, FALSE, inenforcement, result, &consenforced) );

      if( *result == SCIP_CUTOFF )
         break;

      if( !consenforced && inenforcement )
      {
         SCIP_Real viol;

         SCIP_CALL( getConsRelViolation(scip, conss[c], &viol, sol, soltag) );
         if( viol > conshdlrdata->weakcutminviolfactor * maxrelconsviol )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, " constraint <%s> could not be enforced, try again with weak "\
                     "cuts allowed\n", SCIPconsGetName(conss[c])); )

            SCIP_CALL( enforceConstraint(scip, conshdlr, conss[c], sol, soltag, it, TRUE, inenforcement, result, &consenforced) );

            if( consenforced )
               ++conshdlrdata->nweaksepa;  /* TODO maybe this should not be counted per constraint, but per enforcement round? */

            if( *result == SCIP_CUTOFF )
               break;
         }
      }
   }

   SCIPfreeExpriter(&it);

   ENFOLOG( if( enfologfile != NULL ) fflush( enfologfile); )

   /* if having branching scores, then propagate them from expressions with children to variable expressions */
   if( *result == SCIP_BRANCHED )
   {
      /* having result set to branched here means only that we have branching candidates, we still need to do the actual
       * branching
       */
      SCIP_CALL( branching(scip, conshdlr, conss, nconss, maxrelconsviol, sol, soltag, result) );

      /* branching should either have branched: result == SCIP_BRANCHED,
       * or fixed a variable: result == SCIP_REDUCEDDOM,
       * or have registered external branching candidates: result == SCIP_INFEASIBLE,
       * or have not done anything: result == SCIP_DIDNOTFIND
       */
      assert(*result == SCIP_BRANCHED || *result == SCIP_REDUCEDDOM || *result == SCIP_INFEASIBLE || *result == SCIP_DIDNOTFIND);
   }

   ENFOLOG( if( enfologfile != NULL ) fflush( enfologfile); )

   return SCIP_OKAY;
}

/** collect (and print (if debugging enfo)) information on violation in expressions
 *
 * assumes that constraint violations have been computed
 */
static
SCIP_RETCODE analyzeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Longint          soltag,             /**< tag of solution */
   SCIP_Real*            maxabsconsviol,     /**< buffer to store maximal absolute violation of constraints */
   SCIP_Real*            maxrelconsviol,     /**< buffer to store maximal relative violation of constraints */
   SCIP_Real*            minauxviol,         /**< buffer to store minimal (nonzero) violation of auxiliaries */
   SCIP_Real*            maxauxviol,         /**< buffer to store maximal violation of auxiliaries (violation in "extended formulation") */
   SCIP_Real*            maxvarboundviol     /**< buffer to store maximal violation of variable bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;
   SCIP_Real v;
   int c;

   assert(conss != NULL || nconss == 0);
   assert(maxabsconsviol != NULL);
   assert(maxrelconsviol != NULL);
   assert(maxauxviol != NULL);
   assert(maxvarboundviol != NULL);

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   *maxabsconsviol = 0.0;
   *maxrelconsviol = 0.0;
   *minauxviol = SCIPinfinity(scip);
   *maxauxviol = 0.0;
   *maxvarboundviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* skip constraints that are not enabled, deleted, or have separation disabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) || !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      v = getConsAbsViolation(conss[c]);
      *maxabsconsviol = MAX(*maxabsconsviol, v);

      /* skip non-violated constraints */
      if( !isConsViolated(scip, conss[c]) )
         continue;

      SCIP_CALL( getConsRelViolation(scip, conss[c], &v, sol, soltag) );
      *maxrelconsviol = MAX(*maxrelconsviol, v);

      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         SCIP_EXPR_OWNERDATA* ownerdata;
         SCIP_Real auxvarvalue;
         SCIP_Real auxvarlb;
         SCIP_Real auxvarub;
         SCIP_Bool violunder;
         SCIP_Bool violover;
         SCIP_Real origviol;
         SCIP_Real auxviol;
         int e;

         ownerdata = SCIPexprGetOwnerData(expr);
         assert(ownerdata != NULL);

         if( ownerdata->auxvar == NULL )
         {
            /* check violation of variable bounds of original variable */
            if( SCIPisExprVar(scip, expr) )
            {
               SCIP_VAR* var;
               var = SCIPgetVarExprVar(expr);
               auxvarvalue = SCIPgetSolVal(scip, sol, var);
               auxvarlb = SCIPvarGetLbLocal(var);
               auxvarub = SCIPvarGetUbLocal(var);

               origviol = 0.0;
               if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
                  origviol = auxvarlb - auxvarvalue;
               else if( auxvarub < auxvarvalue && !SCIPisInfinity(scip, auxvarub) )
                  origviol = auxvarvalue - auxvarub;
               if( origviol <= 0.0 )
                  continue;

               *maxvarboundviol = MAX(*maxvarboundviol, origviol);

               ENFOLOG(
               SCIPinfoMessage(scip, enfologfile, "var <%s>[%.15g,%.15g] = %.15g", SCIPvarGetName(var), auxvarlb, auxvarub, auxvarvalue);
               if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
                  SCIPinfoMessage(scip, enfologfile, " var >= lb violated by %g", auxvarlb - auxvarvalue);
               if( auxvarub < auxvarvalue && !SCIPisInfinity(scip,  auxvarub) )
                  SCIPinfoMessage(scip, enfologfile, " var <= ub violated by %g", auxvarvalue - auxvarub);
               SCIPinfoMessage(scip, enfologfile, "\n");
               )
            }

            continue;
         }

         auxvarvalue = SCIPgetSolVal(scip, sol, ownerdata->auxvar);
         auxvarlb = SCIPvarGetLbLocal(ownerdata->auxvar);
         auxvarub = SCIPvarGetUbLocal(ownerdata->auxvar);

         /* check violation of variable bounds of auxiliary variable */
         if( auxvarlb - auxvarvalue > *maxvarboundviol && !SCIPisInfinity(scip, -auxvarlb) )
            *maxvarboundviol = auxvarlb - auxvarvalue;
         else if( auxvarvalue - auxvarub > *maxvarboundviol && !SCIPisInfinity(scip,  auxvarub) )
            *maxvarboundviol = auxvarvalue - auxvarub;

         origviol = getExprAbsOrigViolation(scip, expr, sol, &violunder, &violover);

         ENFOLOG(
         if( origviol > 0.0 || auxvarlb > auxvarvalue || auxvarub < auxvarvalue )
         {
            SCIPinfoMessage(scip, enfologfile, "expr ");
            SCIP_CALL( SCIPprintExpr(scip, expr, enfologfile) );
            SCIPinfoMessage(scip, enfologfile, " (%p)[%.15g,%.15g] = %.15g\n", (void*)expr, SCIPexprGetActivity(expr).inf, SCIPexprGetActivity(expr).sup, SCIPexprGetEvalValue(expr));

            SCIPinfoMessage(scip, enfologfile, "  auxvar <%s>[%.15g,%.15g] = %.15g", SCIPvarGetName(ownerdata->auxvar), auxvarlb, auxvarub, auxvarvalue);
            if( origviol > 0.0 )
               SCIPinfoMessage(scip, enfologfile, " auxvar %s expr violated by %g", violunder ? ">=" : "<=", origviol);
            if( auxvarlb > auxvarvalue && !SCIPisInfinity(scip, -auxvarlb) )
               SCIPinfoMessage(scip, enfologfile, " auxvar >= auxvar's lb violated by %g", auxvarlb - auxvarvalue);
            if( auxvarub < auxvarvalue && !SCIPisInfinity(scip,  auxvarub) )
               SCIPinfoMessage(scip, enfologfile, " auxvar <= auxvar's ub violated by %g", auxvarvalue - auxvarub);
            SCIPinfoMessage(scip, enfologfile, "\n");
         }
         )

         /* no violation w.r.t. the original variables -> skip expression */
         if( origviol == 0.0 )
            continue;

         /* compute aux-violation for each nonlinear handlers */
         for( e = 0; e < ownerdata->nenfos; ++e )
         {
            SCIP_NLHDLR* nlhdlr;

            /* eval in auxvars is only defined for nlhdrs that separate; there might not even be auxvars otherwise */
            if( (ownerdata->enfos[e]->nlhdlrparticipation & SCIP_NLHDLR_METHOD_SEPABOTH) == 0 )
               continue;

            nlhdlr = ownerdata->enfos[e]->nlhdlr;
            assert(nlhdlr != NULL);

            /* evaluate the expression w.r.t. the nlhdlrs auxiliary variables */
            SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, ownerdata->enfos[e]->nlhdlrexprdata, &ownerdata->enfos[e]->auxvalue, sol) );

            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "  nlhdlr <%s> = %.15g", SCIPnlhdlrGetName(nlhdlr), ownerdata->enfos[e]->auxvalue); )

            auxviol = getExprAbsAuxViolation(scip, expr, ownerdata->enfos[e]->auxvalue, sol, &violunder, &violover);

            if( auxviol > 0.0 )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, " auxvar %s nlhdlr-expr violated by %g", violover ? "<=" : ">=", auxviol); )
               *maxauxviol = MAX(*maxauxviol, auxviol);
               *minauxviol = MIN(*minauxviol, auxviol);
            }
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "\n"); )
         }
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
} /*lint !e715*/

/** enforcement of constraints called by enfolp and enforelax */
static
SCIP_RETCODE consEnfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real maxabsconsviol;
   SCIP_Real maxrelconsviol;
   SCIP_Real minauxviol;
   SCIP_Real maxauxviol;
   SCIP_Real maxvarboundviol;
   SCIP_Longint soltag;
   int nnotify;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlr != NULL);

   soltag = SCIPgetExprNewSoltag(scip);

   *result = SCIP_FEASIBLE;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
         *result = SCIP_INFEASIBLE;
   }

   if( *result == SCIP_FEASIBLE )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: all expr-constraints feasible, skip enforcing\n",
               SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )
      return SCIP_OKAY;
   }

   SCIP_CALL( analyzeViolation(scip, conss, nconss, sol, soltag, &maxabsconsviol, &maxrelconsviol,
            &minauxviol, &maxauxviol, &maxvarboundviol) );

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: enforcing constraints with max conssviol=%e (rel=%e), "\
            "auxviolations in %g..%g, variable bounds violated by at most %g, LP feastol=%e\n",
            SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), maxabsconsviol, maxrelconsviol, minauxviol, maxauxviol,
            maxvarboundviol, SCIPgetLPFeastol(scip)); )

   assert(maxvarboundviol <= SCIPgetLPFeastol(scip));

   /* try to propagate */
   if( conshdlrdata->propinenforce )
   {
      SCIP_RESULT propresult;
      int nchgbds = 0;

      SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds) );

      if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
      {
         *result = propresult;
         return SCIP_OKAY;
      }
   }

   /* tighten the LP tolerance if violation in variables bounds is larger than aux-violation (max |expr - auxvar| over
    * all violated expr/auxvar in violated constraints)
    */
   if( conshdlrdata->tightenlpfeastol && maxvarboundviol > maxauxviol && SCIPisPositive(scip, SCIPgetLPFeastol(scip)) &&
         sol == NULL )
   {
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxvarboundviol / 2.0, SCIPgetLPFeastol(scip) / 2.0)));
      ++conshdlrdata->ntightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " variable bound violation %g larger than auxiliary violation %g, "\
               "reducing LP feastol to %g\n", maxvarboundviol, maxauxviol, SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   /* tighten the LP tolerance if violation in auxiliaries is below LP feastol, as we could have problems to find a cut
    * with violation above LP tolerance (especially when auxviolation is below 10*eps = ROWPREP_SCALEUP_VIOLNONZERO in misc_rowprep.c)
    */
   if( conshdlrdata->tightenlpfeastol && maxauxviol < SCIPgetLPFeastol(scip) && SCIPisPositive(scip, SCIPgetLPFeastol(scip)) && sol == NULL )
   {
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), maxauxviol/2.0));
      ++conshdlrdata->ntightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " auxiliary violation %g below LP feastol, reducing LP feastol to %g\n", maxauxviol, SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   SCIP_CALL( enforceConstraints(scip, conshdlr, conss, nconss, sol, soltag, TRUE, maxrelconsviol, result) );

   if( *result == SCIP_CUTOFF || *result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM || *result == SCIP_BRANCHED ||
         *result == SCIP_INFEASIBLE )
      return SCIP_OKAY;

   assert(*result == SCIP_DIDNOTFIND);

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " could not enforce violation %g in regular ways, LP feastol=%g, "\
            "becoming desperate now...\n", maxabsconsviol, SCIPgetLPFeastol(scip)); )

   if( conshdlrdata->tightenlpfeastol && SCIPisPositive(scip, maxvarboundviol) && SCIPisPositive(scip, SCIPgetLPFeastol(scip)) && sol == NULL )
   {
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxvarboundviol / 2.0, SCIPgetLPFeastol(scip) / 2.0)));
      ++conshdlrdata->ntightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " variable bounds are violated by more than eps, reduced LP "\
               "feasibility tolerance to %g\n", SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   if( conshdlrdata->tightenlpfeastol && SCIPisPositive(scip, maxauxviol) && SCIPisPositive(scip,
            SCIPgetLPFeastol(scip)) && sol == NULL )
   {
      /* try whether tighten the LP feasibility tolerance could help
       * maybe it is just some cut that hasn't been taken into account sufficiently
       * in the next enforcement round, we would then also allow even weaker cuts, as we want a minimal cut violation of LP's feastol
       * unfortunately, we do not know the current LP solution primal infeasibility, so sometimes this just repeats without effect
       * until the LP feastol reaches epsilon
       * (this is similar to the "tighten the LP tolerance if violation in auxiliaries is below LP feastol..." case above, but applies
       * when maxauxviol is above LP feastol)
       */
      SCIPsetLPFeastol(scip, MAX(SCIPepsilon(scip), MIN(maxauxviol / 2.0, SCIPgetLPFeastol(scip) / 10.0)));
      ++conshdlrdata->ndesperatetightenlp;

      *result = SCIP_SOLVELP;

      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " reduced LP feasibility tolerance to %g and hope\n", SCIPgetLPFeastol(scip)); )

      return SCIP_OKAY;
   }

   /* try to propagate, if not tried above TODO(?) allow to disable this as well */
   if( !conshdlrdata->propinenforce )
   {
      SCIP_RESULT propresult;
      int nchgbds = 0;

      SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds) );

      if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
      {
         *result = propresult;
         return SCIP_OKAY;
      }
   }

   /* could not find branching candidates even when looking at minimal violated (>eps) expressions
    * now look if we find any unfixed variable that we could still branch on
    */
   SCIP_CALL( registerBranchingCandidatesAllUnfixed(scip, conshdlr, conss, nconss, &nnotify) );

   if( nnotify > 0 )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, " registered %d unfixed variables as branching candidates\n", nnotify); )
      ++conshdlrdata->ndesperatebranch;

      *result = SCIP_INFEASIBLE;  /* enforceConstraints may have changed it to SCIP_DIDNOTFIND */

      return SCIP_OKAY;
   }

   /* if everything is fixed in violated constraints, then let's cut off the node
    * - bound tightening with all vars fixed should prove cutoff, but interval arithmetic overestimates and so the
    *   result may not be conclusive (when constraint violations are small)
    * - if tightenlpfeastol=FALSE, then the LP solution that we try to enforce here may just not be within bounds
    *   sufficiently (see st_e40)
    * - but if the LP solution is really within bounds and since variables are fixed, cutting off the node is actually
    *   not "desperate", but a pretty obvious thing to do
    */
   ENFOLOG( SCIPinfoMessage(scip, enfologfile, " enforcement with max. violation %g failed; cutting off node\n", maxabsconsviol); )
   *result = SCIP_CUTOFF;

   /* it's only "desperate" if the LP solution does not coincide with variable fixings (should we use something tighter than epsilon here?) */
   if( !SCIPisZero(scip, maxvarboundviol) )
      ++conshdlrdata->ndesperatecutoff;

   return SCIP_OKAY;
}

/** separation for all violated constraints to be used by SEPA callbacks */
static
SCIP_RETCODE consSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_Longint soltag;
   SCIP_Bool haveviol = FALSE;
   int c;

   *result = SCIP_DIDNOTFIND;

   soltag = SCIPgetExprNewSoltag(scip);

   /* compute violations */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);

      /* skip constraints that are not enabled, deleted, or have separation disabled */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) || !SCIPconsIsSeparationEnabled(conss[c]) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
         haveviol = TRUE;
   }

   /* if none of our constraints are violated, don't attempt separation */
   if( !haveviol )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: skip separation of non-violated constraints\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )
      return SCIP_OKAY;
   }

   ENFOLOG( SCIPinfoMessage(scip, enfologfile, "node %lld: separation\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip))); )

   /* call separation */
   SCIP_CALL( enforceConstraints(scip, conshdlr, conss, nconss, sol, soltag, FALSE, SCIP_INVALID, result) );

   return SCIP_OKAY;
}

/** hash key retrieval function for bilinear term entries */
static
SCIP_DECL_HASHGETKEY(bilinearTermsGetHashkey)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int idx;

   conshdlrdata = (SCIP_CONSHDLRDATA*)userptr;
   assert(conshdlrdata != NULL);

   idx = ((int)(size_t)elem) - 1;
   assert(idx >= 0 && idx < conshdlrdata->nbilinterms);

   return (void*)&conshdlrdata->bilinterms[idx];
}

/** returns TRUE iff the bilinear term entries are equal */
static
SCIP_DECL_HASHKEYEQ(bilinearTermsIsHashkeyEq)
{  /*lint --e{715}*/
   SCIP_CONSNONLINEAR_BILINTERM* entry1;
   SCIP_CONSNONLINEAR_BILINTERM* entry2;

   /* get corresponding entries */
   entry1 = (SCIP_CONSNONLINEAR_BILINTERM*)key1;
   entry2 = (SCIP_CONSNONLINEAR_BILINTERM*)key2;
   assert(entry1->x != NULL && entry1->y != NULL);
   assert(entry2->x != NULL && entry2->y != NULL);
   assert(SCIPvarCompare(entry1->x, entry1->y) < 1);
   assert(SCIPvarCompare(entry2->x, entry2->y) < 1);

   return entry1->x == entry2->x && entry1->y == entry2->y;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(bilinearTermsGetHashkeyVal)
{  /*lint --e{715}*/
   SCIP_CONSNONLINEAR_BILINTERM* entry;

   entry = (SCIP_CONSNONLINEAR_BILINTERM*)key;
   assert(entry->x != NULL && entry->y != NULL);
   assert(SCIPvarCompare(entry->x, entry->y) < 1);

   return SCIPhashTwo(SCIPvarGetIndex(entry->x), SCIPvarGetIndex(entry->y));
}

/** compare two auxiliary expressions
 *
 *  Compares auxiliary variables, followed by coefficients, and then constants.
 */
static
SCIP_DECL_SORTPTRCOMP(auxexprComp)
{
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr1 = (SCIP_CONSNONLINEAR_AUXEXPR*)elem1;
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr2 = (SCIP_CONSNONLINEAR_AUXEXPR*)elem2;
   int compvars;
   int i;

   /* compare the auxiliary variables */
   compvars = SCIPvarCompare(auxexpr1->auxvar, auxexpr2->auxvar); /* TODO can one of these be NULL? */

   if( compvars != 0 )
      return compvars;

   /* compare the coefficients and constants */
   for( i = 0; i < 3; ++i )
   {
      if( auxexpr1->coefs[i] != auxexpr2->coefs[i] )
         return auxexpr1->coefs[i] < auxexpr2->coefs[i] ? -1 : 1;
   }

   return auxexpr1->cst < auxexpr2->cst ? -1 : auxexpr1->cst == auxexpr2->cst ? 0 : 1;
}

/* add an auxiliary expression to a bilinear term */
static
SCIP_RETCODE bilinTermAddAuxExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< nonlinear constraint handler data */
   SCIP_CONSNONLINEAR_BILINTERM* term,       /**< bilinear term */
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr,      /**< auxiliary expression to add */
   SCIP_Bool*            added               /**< pointer to store whether auxexpr has been added */
   )
{
   SCIP_Bool found;
   int pos;
   int i;

   *added = FALSE;

   /* check if auxexpr has already been added to term */
   if( term->nauxexprs == 0 )
   {
      found = FALSE;
      pos = 0;
   }
   else
   {
      found = SCIPsortedvecFindPtr((void**)term->aux.exprs, auxexprComp, auxexpr, term->nauxexprs, &pos);
   }

   if( !found )
   {
      if( term->nauxexprs >= conshdlrdata->bilinmaxnauxexprs )
         return SCIP_OKAY;

      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &term->aux.exprs, &term->auxexprssize, term->nauxexprs + 1) );
      assert(term->auxexprssize >= term->nauxexprs + 1);

      /* insert expression at the correct position */
      for( i = term->nauxexprs; i > pos; --i )
      {
         term->aux.exprs[i] = term->aux.exprs[i-1];
      }
      term->aux.exprs[pos] = auxexpr;
      ++(term->nauxexprs);
      *added = TRUE;
   }
   else
   {
      term->aux.exprs[pos]->underestimate |= auxexpr->underestimate;
      term->aux.exprs[pos]->overestimate  |= auxexpr->overestimate;
   }

   return SCIP_OKAY;
}

/** iterates through all expressions of all nonlinear constraints and adds the corresponding bilinear terms to the hash table */
static
SCIP_RETCODE bilinearTermsInsertAll(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< nonlinear constraints */
   int                   nconss              /**< total number of nonlinear constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPRITER* it;
   int c;

   assert(conss != NULL || nconss == 0);

   if( nconss == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether the bilinear terms have been stored already */
   if( conshdlrdata->bilinterms != NULL )
      return SCIP_OKAY;

   /* create and initialize iterator */
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR);

   /* iterate through all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      SCIP_EXPR* expr;

      assert(conss != NULL && conss[c] != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* iterate through all expressions */
      for( expr = SCIPexpriterRestartDFS(it, consdata->expr); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      {
         SCIP_EXPR** children = SCIPexprGetChildren(expr);
         SCIP_VAR* x = NULL;
         SCIP_VAR* y = NULL;

         /* check whether the expression is of the form f(..)^2 */
         if( SCIPisExprPower(scip, expr) && SCIPgetExponentExprPow(expr) == 2.0 )
         {
            x = SCIPgetExprAuxVarNonlinear(children[0]);
            y = x;
         }
         /* check whether the expression is of the form f(..) * g(..) */
         else if( SCIPisExprProduct(scip, expr) && SCIPexprGetNChildren(expr) == 2 )
         {
            x = SCIPgetExprAuxVarNonlinear(children[0]);
            y = SCIPgetExprAuxVarNonlinear(children[1]);
         }

         /* add variables to the hash table */
         if( x != NULL && y != NULL )
         {
            SCIP_CALL( SCIPinsertBilinearTermExistingNonlinear(scip, conshdlr, x, y, SCIPgetExprAuxVarNonlinear(expr),
               SCIPgetExprNLocksPosNonlinear(expr), SCIPgetExprNLocksNegNonlinear(expr)) );
         }
      }
   }

   /* release iterator */
   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** store x, y and the locks in a new bilinear term */
static
SCIP_RETCODE bilinearTermsInsertEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_VAR*             x,                  /**< the first variable */
   SCIP_VAR*             y,                  /**< the second variable */
   int                   nlockspos,          /**< number of positive locks of the bilinear term */
   int                   nlocksneg,          /**< number of negative locks of the bilinear term */
   int*                  idx,                /**< pointer to store the position of the term in bilinterms array */
   SCIP_Bool             existing            /**< whether the term exists explicitly in the problem */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM* term;

   assert(conshdlr != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(nlockspos >= 0);
   assert(nlocksneg >= 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* ensure that x.index <= y.index */
   if( SCIPvarCompare(x, y) == 1 )
   {
      SCIPswapPointers((void**)&x, (void**)&y);
   }
   assert(SCIPvarCompare(x, y) < 1);

   *idx = SCIPgetBilinTermIdxNonlinear(conshdlr, x, y);

   /* update or create the term */
   if( *idx >= 0 )
   { /* the term has already been added */
      assert(conshdlrdata->bilinterms[*idx].x == x);
      assert(conshdlrdata->bilinterms[*idx].y == y);

      /* get term and add locks */
      term = &conshdlrdata->bilinterms[*idx];
      assert(existing <= term->existing); /* implicit terms are added after existing ones */
      term->nlockspos += nlockspos;
      term->nlocksneg += nlocksneg;
   }
   else
   { /* this is the first time we encounter this product */
      /* ensure size of bilinterms array */
      SCIP_CALL( SCIPensureBlockMemoryArray(scip, &conshdlrdata->bilinterms, &conshdlrdata->bilintermssize, conshdlrdata->nbilinterms + 1) );

      *idx = conshdlrdata->nbilinterms;

      /* get term and set values in the created bilinear term */
      term = &conshdlrdata->bilinterms[*idx];
      assert(term != NULL);
      term->x = x;
      term->y = y;
      term->nauxexprs = 0;
      term->auxexprssize = 0;
      term->nlockspos = nlockspos;
      term->nlocksneg = nlocksneg;
      term->existing = existing;
      if( existing )
         term->aux.var = NULL;
      else
         term->aux.exprs = NULL;

      /* increase the total number of bilinear terms */
      ++(conshdlrdata->nbilinterms);

      /* save to the hashtable */
      if( conshdlrdata->bilinhashtable == NULL )
      {
         SCIP_CALL( SCIPhashtableCreate(&conshdlrdata->bilinhashtable, SCIPblkmem(scip), conshdlrdata->nbilinterms,
               bilinearTermsGetHashkey, bilinearTermsIsHashkeyEq, bilinearTermsGetHashkeyVal,
               (void*)conshdlrdata) );
      }
      assert(conshdlrdata->bilinhashtable != NULL);

      /* insert the index of the bilinear term into the hash table; note that the index of the i-th element is (i+1)
       * because zero can not be inserted into hash table
       */
      SCIP_CALL( SCIPhashtableInsert(conshdlrdata->bilinhashtable, (void*)(size_t)(*idx + 1)) ); /*lint !e571 !e776*/

      /* capture product variables */
      SCIP_CALL( SCIPcaptureVar(scip, x) );
      SCIP_CALL( SCIPcaptureVar(scip, y) );
   }

   return SCIP_OKAY;
}

/** frees array of bilinear terms and hash table */
static
SCIP_RETCODE bilinearTermsFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   int i;
   int j;

   assert(conshdlrdata != NULL);

   /* check whether bilinear terms have been stored */
   if( conshdlrdata->bilinterms == NULL )
   {
      assert(conshdlrdata->bilinterms == NULL);
      assert(conshdlrdata->nbilinterms == 0);
      assert(conshdlrdata->bilintermssize == 0);

      return SCIP_OKAY;
   }

   /* release variables */
   for( i = 0; i < conshdlrdata->nbilinterms; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].y) );
      SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].x) );

      for( j = 0; j < conshdlrdata->bilinterms[i].nauxexprs; ++j )
      {
         if( conshdlrdata->bilinterms[i].aux.exprs[j]->auxvar != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].aux.exprs[j]->auxvar) );
         }
         SCIPfreeBlockMemory(scip, &(conshdlrdata->bilinterms[i].aux.exprs[j]));
      }

      if( conshdlrdata->bilinterms[i].nauxexprs > 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &(conshdlrdata->bilinterms[i].aux.exprs), conshdlrdata->bilinterms[i].auxexprssize);
         continue;
      }

      /* the rest is for simple terms with a single auxvar */

      /* it might be that there is a bilinear term without a corresponding auxiliary variable */
      if( conshdlrdata->bilinterms[i].aux.var != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->bilinterms[i].aux.var) );
      }
   }

   /* free hash table */
   if( conshdlrdata->bilinhashtable != NULL )
   {
      SCIPhashtableFree(&conshdlrdata->bilinhashtable);
   }

   /* free bilinterms array; reset counters */
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bilinterms, conshdlrdata->bilintermssize);
   conshdlrdata->nbilinterms = 0;
   conshdlrdata->bilintermssize = 0;

   return SCIP_OKAY;
}

/*
 * vertex polyhedral separation
 */

/** builds LP used to compute facets of the convex envelope of vertex-polyhedral functions */
static
SCIP_RETCODE buildVertexPolyhedralSeparationLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of (unfixed) variables in vertex-polyhedral functions */
   SCIP_LPI**            lp                  /**< pointer to store created LP */
   )
{
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* val;
   int* beg;
   int* ind;
   unsigned int nnonz;
   unsigned int ncols;
   unsigned int nrows;
   unsigned int i;
   unsigned int k;

   assert(scip != NULL);
   assert(lp != NULL);
   assert(nvars > 0);
   assert(nvars <= SCIP_MAXVERTEXPOLYDIM);

   SCIPdebugMsg(scip, "Building LP for computing facets of convex envelope of vertex-polyhedral function\n");

   /* create lpi to store the LP */
   SCIP_CALL( SCIPlpiCreate(lp, SCIPgetMessagehdlr(scip), "facet finding LP", SCIP_OBJSEN_MINIMIZE) );

   nrows = (unsigned int)nvars + 1;
   ncols = POWEROFTWO((unsigned int)nvars);
   nnonz = (ncols * (nrows + 1)) / 2;

   /* allocate necessary memory; set obj, lb, and ub to zero */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &obj, ncols) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, nnonz) );

   /* calculate nonzero entries in the LP */
   for( i = 0, k = 0; i < ncols; ++i )
   {
      int row;
      unsigned int a;

      /* an upper bound of 1.0 is implied by the last row, but I presume that LP solvers prefer unbounded variables */
      ub[i] = SCIPlpiInfinity(*lp);

      SCIPdebugMsg(scip, "col %u starts at position %u\n", i, k);
      beg[i] = (int)k;
      row = 0;

      /* iterate through the bit representation of i */
      a = 1;
      while( a <= i )
      {
         if( (a & i) != 0 )
         {
            val[k] = 1.0;
            ind[k] = row;

            SCIPdebugMsg(scip, " val[%d][%u] = 1 (position  %u)\n", row, i, k);

            ++k;
         }

         a <<= 1;
         ++row;
         assert(0 <= row && row <= SCIP_MAXVERTEXPOLYDIM);
         assert(POWEROFTWO(row) == a);
      }

      /* put 1 as a coefficient for sum_{i} \lambda_i = 1 row (last row) */
      val[k] = 1.0;
      ind[k] = (int)nrows - 1;
      ++k;
      SCIPdebugMsg(scip, " val[%u][%u] = 1 (position  %u)\n", nrows - 1, i, k);
   }
   assert(k == nnonz);

   /* load all data into LP interface
    * we can assume nrows (=nvars+1) <= ncols (=2^nvars), so we can pass lb as dummy lhs and rhs
    */
   assert(nrows <= ncols);
   SCIP_CALL( SCIPlpiLoadColLP(*lp, SCIP_OBJSEN_MINIMIZE,
      (int)ncols, obj, lb, ub, NULL,
      (int)nrows, lb, lb, NULL,
      (int)nnonz, beg, ind, val) );

   /* for the last row, we can set the rhs to 1.0 already */
   ind[0] = (int)nrows - 1;
   val[0] = 1.0;
   SCIP_CALL( SCIPlpiChgSides(*lp, 1, ind, val, val) );

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** the given facet might not be a valid under(over)estimator, because of numerics and bad fixings; we compute \f$
 * \max_{v \in V} f(v) - (\alpha v + \beta) \f$ (\f$\max_{v \in V} \alpha v + \beta - f(v) \f$) where \f$ V \f$ is the
 * set of vertices of the domain
 */
static
SCIP_Real computeVertexPolyhedralMaxFacetError(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             overestimate,       /**< whether we check for an over or underestimator */
   SCIP_Real*            funvals,            /**< array containing the evaluation of the function at all corners, length: 2^nvars */
   SCIP_Real*            box,                /**< box for which facet was computed, length: 2*nallvars */
   int                   nallvars,           /**< number of all variables */
   int                   nvars,              /**< number of unfixed variables */
   int*                  nonfixedpos,        /**< indices of unfixed variables, length: nvars */
   SCIP_Real*            facetcoefs,         /**< current facet candidate's coefficients, length: nallvars */
   SCIP_Real             facetconstant       /**< current facet candidate's constant, length: nallvars */
   )
{
   SCIP_Real maxerror;
   SCIP_Real facetval;
   SCIP_Real funval;
   SCIP_Real error;
   unsigned int i;
   unsigned int ncorners;
   unsigned int prev;

   assert(scip != NULL);
   assert(funvals != NULL);
   assert(box != NULL);
   assert(nonfixedpos != NULL);
   assert(facetcoefs != NULL);

   ncorners = POWEROFTWO(nvars);
   maxerror = 0.0;

   /* check the origin (all variables at lower bound) */
   facetval = facetconstant;
   for( i = 0; i < (unsigned int) nallvars; ++i )
      facetval += facetcoefs[i] * box[2*i];

   /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
   funval = funvals[0];
   if( overestimate )
      error = funval - facetval;
   else
      error = facetval - funval;

   /* update maximum error */
   maxerror = MAX(error, maxerror);

   prev = 0;
   for( i = 1; i < ncorners; ++i )
   {
      unsigned int gray;
      unsigned int diff;
      unsigned int pos;
      int origpos;

      gray = i ^ (i >> 1);
      diff = gray ^ prev;

      /* compute position of unique 1 of diff */
      pos = 0;
      while( (diff >>= 1) != 0 )
         ++pos;
      assert(pos < (unsigned int)nvars);

      origpos = nonfixedpos[pos];

      if( gray > prev )
         facetval += facetcoefs[origpos] * (box[2*origpos+1] - box[2*origpos]);
      else
         facetval -= facetcoefs[origpos] * (box[2*origpos+1] - box[2*origpos]);

      /* compute largest/smallest possible value of function, depending on whether we are over/under-estimating */
      funval = funvals[gray];
      if( overestimate )
         error = funval - facetval;
      else
         error = facetval - funval;

      /* update  maximum error */
      maxerror = MAX(error, maxerror);

      prev = gray;
   }

   SCIPdebugMsg(scip, "maximum error of facet: %2.8e\n", maxerror);

   return maxerror;
}

/** computes a facet of the convex or concave envelope of a vertex polyhedral function by solving an LP */  /*lint -e{715}*/
static
SCIP_RETCODE computeVertexPolyhedralFacetLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   int*                  nonfixedpos,        /**< indices of nonfixed variables */
   SCIP_Real*            funvals,            /**< values of function in all corner points (w.r.t. nonfixed variables) */
   int                   nvars,              /**< number of nonfixed variables */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an zero'ed array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
   )
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_LPI* lp;
   SCIP_Real* aux; /* used to transform x^* and then to store LP solution */
   int* inds;
   int ncols;
   int nrows;
   int i;
   SCIP_Real facetvalue;
   SCIP_Real mindomwidth;
   SCIP_RETCODE lpsolveretcode;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(xstar != NULL);
   assert(box != NULL);
   assert(nonfixedpos != NULL);
   assert(funvals != NULL);
   assert(nvars >= 0);
   assert(nvars <= SCIP_MAXVERTEXPOLYDIM);
   assert(success != NULL);
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->vp_randnumgen == NULL && conshdlrdata->vp_maxperturb > 0.0 )
   {
      SCIP_CALL( SCIPcreateRandom(scip, &conshdlrdata->vp_randnumgen, VERTEXPOLY_RANDNUMINITSEED, TRUE) );
   }

   /* construct an LP for this size, if not having one already */
   if( conshdlrdata->vp_lp[nvars] == NULL )
   {
      SCIP_CALL( buildVertexPolyhedralSeparationLP(scip, nvars, &conshdlrdata->vp_lp[nvars]) );
   }
   lp = conshdlrdata->vp_lp[nvars];
   assert(lp != NULL);

   /* get number of cols and rows of separation lp */
   SCIP_CALL( SCIPlpiGetNCols(lp, &ncols) );
   SCIP_CALL( SCIPlpiGetNRows(lp, &nrows) );

   /* number of columns should equal the number of corners = 2^nvars */
   assert(ncols == (int)POWEROFTWO(nvars));

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &aux, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, ncols) );

   /*
    * set up the described LP on the transformed space
    */

   for( i = 0; i < ncols; ++i )
      inds[i] = i;

   /* compute T^-1(x^*), i.e. T^-1(x^*)_i = (x^*_i - lb_i)/(ub_i - lb_i) */
   mindomwidth = 2*SCIPinfinity(scip);
   for( i = 0; i < nrows-1; ++i )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;
      int varpos;

      assert(i < nvars);

      varpos = nonfixedpos[i];
      lb = box[2 * varpos];
      ub = box[2 * varpos + 1];
      solval = xstar[varpos];

      if( ub - lb < mindomwidth )
         mindomwidth = ub - lb;

      /* explicitly handle solution which violate bounds of variables (this can happen because of tolerances) */
      if( solval <= lb )
         aux[i] = 0.0;
      else if( solval >= ub )
         aux[i] = 1.0;
      else
         aux[i] = (solval - lb) / (ub - lb);

      /* perturb point to hopefully obtain a facet of the convex envelope */
      if( conshdlrdata->vp_maxperturb > 0.0 )
      {
         assert(conshdlrdata->vp_randnumgen != NULL);

         if( aux[i] == 1.0 )
            aux[i] -= SCIPrandomGetReal(conshdlrdata->vp_randnumgen, 0.0, conshdlrdata->vp_maxperturb);
         else if( aux[i] == 0.0 )
            aux[i] += SCIPrandomGetReal(conshdlrdata->vp_randnumgen, 0.0, conshdlrdata->vp_maxperturb);
         else
         {
            SCIP_Real perturbation;

            perturbation = MIN( aux[i], 1.0 - aux[i] ) / 2.0;
            perturbation = MIN( perturbation, conshdlrdata->vp_maxperturb );
            aux[i] += SCIPrandomGetReal(conshdlrdata->vp_randnumgen, -perturbation, perturbation);
         }
         assert(0.0 < aux[i] && aux[i] < 1.0);
      }

      SCIPdebugMsg(scip, "LP row %d in [%e, %e]\n", i, aux[i], aux[i]);
   }

   /* update LP */
   SCIP_CALL( SCIPlpiChgObj(lp, ncols, inds, funvals) );
   SCIP_CALL( SCIPlpiChgSides(lp, nrows-1, inds, aux, aux) );
   SCIP_CALL( SCIPlpiChgObjsen(lp, overestimate ? SCIP_OBJSEN_MAXIMIZE : SCIP_OBJSEN_MINIMIZE) );

   /* we can stop the LP solve if will not meet the target value anyway, but only if xstar hasn't been perturbed */
   if( conshdlrdata->vp_maxperturb == 0.0 && !SCIPisInfinity(scip, REALABS(targetvalue)) )
   {
      SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_OBJLIM, targetvalue) );
   }
   /* set an iteration limit so we do not run forever */
   SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_LPITLIM, 100*ncols) );
   /* since we work with the dual of the LP, primal feastol determines how much we want the computed facet to be the best possible one */
   SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_FEASTOL, SCIPfeastol(scip)) );
   /* since we work with the dual of the LP, dual feastol determines validity of the facet
    * if some ub-lb is small, we need higher accuracy, since below we divide coefs by ub-lb (we moved and scaled the box)
    * thus, we set the dual feastol to be between SCIPepsilon and SCIPfeastol
    */
   SCIP_CALL( SCIPlpiSetRealpar(lp, SCIP_LPPAR_DUALFEASTOL, MIN(SCIPfeastol(scip), MAX(SCIPepsilon(scip), mindomwidth * SCIPfeastol(scip)))) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPlpiSetIntpar(lp, SCIP_LPPAR_LPINFO, 1) );
#endif

   /*
    * solve the LP and store the resulting facet for the transformed space
    */
   if( conshdlrdata->vp_dualsimplex )
   {
      lpsolveretcode = SCIPlpiSolveDual(lp);
   }
   else
   {
      lpsolveretcode = SCIPlpiSolvePrimal(lp);
   }
   if( lpsolveretcode == SCIP_LPERROR )
   {
      SCIPdebugMsg(scip, "LP error, aborting.\n");
      goto CLEANUP;
   }
   SCIP_CALL( lpsolveretcode );

   /* any dual feasible solution should provide a valid estimator (and a dual optimal one a facet) */
   if( !SCIPlpiIsDualFeasible(lp) )
   {
      SCIPdebugMsg(scip, "LP not solved to dual feasibility, aborting.\n");
      goto CLEANUP;
   }

   /* get dual solution (facet of convex envelope); again, we have to be careful since the LP can have more rows and
    * columns than needed, in particular, \bar \beta is the last dual multiplier
    */
   SCIP_CALL( SCIPlpiGetSol(lp, NULL, NULL, aux, NULL, NULL) );

   for( i = 0; i < nvars; ++i )
      facetcoefs[nonfixedpos[i]] = aux[i];
   /* last dual multiplier is the constant */
   *facetconstant = aux[nrows - 1];

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "facet for the transformed problem: ");
   for( i = 0; i < nallvars; ++i )
   {
      SCIPdebugMsgPrint(scip, "%3.4e * x%d + ", facetcoefs[i], i);
   }
   SCIPdebugMsgPrint(scip, "%3.4e\n", *facetconstant);
#endif

   /*
    * transform the facet to original space and compute value at x^*, i.e., alpha x + beta
    */

   SCIPdebugMsg(scip, "facet in orig. space: ");

   facetvalue = 0.0;
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      int varpos;

      varpos = nonfixedpos[i];
      lb = box[2 * varpos];
      ub = box[2 * varpos + 1];
      assert(!SCIPisEQ(scip, lb, ub));

      /* alpha_i := alpha_bar_i / (ub_i - lb_i) */
      facetcoefs[varpos] = facetcoefs[varpos] / (ub - lb);

      /* beta = beta_bar - sum_i alpha_i * lb_i */
      *facetconstant -= facetcoefs[varpos] * lb;

      /* evaluate */
      facetvalue += facetcoefs[varpos] * xstar[varpos];

      SCIPdebugMsgPrint(scip, "%3.4e * x%d + ", facetcoefs[varpos], varpos);
   }
   SCIPdebugMsgPrint(scip, "%3.4e ", *facetconstant);

   /* add beta to the facetvalue: at this point in the code, facetvalue = g(x^*) */
   facetvalue += *facetconstant;

   SCIPdebugMsgPrint(scip, "has value %g, target = %g\n", facetvalue, targetvalue);

    /* if overestimate, then we want facetvalue < targetvalue
    * if underestimate, then we want facetvalue > targetvalue
    * if none holds, give up
    * so maybe here we should check against the minimal violation
    */
   if( overestimate == (facetvalue > targetvalue) )
   {
      SCIPdebugMsg(scip, "missed the target, facetvalue %g targetvalue %g, overestimate=%u\n", facetvalue, targetvalue, overestimate);
      goto CLEANUP;
   }

   /* if we made it until here, then we have a nice facet */
   *success = TRUE;

CLEANUP:
   /* free allocated memory */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &aux);

   return SCIP_OKAY;
}

/** computes a facet of the convex or concave envelope of a univariant vertex polyhedral function
 *
 * In other words, compute the line that passes through two given points.
 */
static
SCIP_RETCODE computeVertexPolyhedralFacetUnivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             left,               /**< left coordinate */
   SCIP_Real             right,              /**< right coordinate */
   SCIP_Real             funleft,            /**< value of function in left coordinate */
   SCIP_Real             funright,           /**< value of function in right coordinate */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoef,          /**< buffer to store coefficient of facet defining inequality */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
   )
{
   assert(scip != NULL);
   assert(SCIPisLE(scip, left, right));
   assert(!SCIPisInfinity(scip, -left));
   assert(!SCIPisInfinity(scip, right));
   assert(SCIPisFinite(funleft) && funleft != SCIP_INVALID);
   assert(SCIPisFinite(funright) && funright != SCIP_INVALID);
   assert(success != NULL);
   assert(facetcoef != NULL);
   assert(facetconstant != NULL);

   *facetcoef = (funright - funleft) / (right - left);
   *facetconstant = funleft - *facetcoef * left;

   *success = TRUE;

   return SCIP_OKAY;
}

/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 *
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
static
SCIP_RETCODE computeHyperplaneThreePoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a1,                 /**< first coordinate of a */
   SCIP_Real             a2,                 /**< second coordinate of a */
   SCIP_Real             a3,                 /**< third coordinate of a */
   SCIP_Real             b1,                 /**< first coordinate of b */
   SCIP_Real             b2,                 /**< second coordinate of b */
   SCIP_Real             b3,                 /**< third coordinate of b */
   SCIP_Real             c1,                 /**< first coordinate of c */
   SCIP_Real             c2,                 /**< second coordinate of c */
   SCIP_Real             c3,                 /**< third coordinate of c */
   SCIP_Real*            alpha,              /**< coefficient of first coordinate */
   SCIP_Real*            beta,               /**< coefficient of second coordinate */
   SCIP_Real*            gamma_,             /**< coefficient of third coordinate */
   SCIP_Real*            delta               /**< constant right-hand side */
   )
{
   assert(scip != NULL);
   assert(alpha != NULL);
   assert(beta  != NULL);
   assert(gamma_ != NULL);
   assert(delta != NULL);

   *alpha  = -b3*c2 + a3*(-b2+c2) + a2*(b3-c3) + b2*c3;
   *beta   = -(-b3*c1 + a3*(-b1+c1) + a1*(b3-c3) + b1*c3);
   *gamma_ = -a2*b1 + a1*b2 + a2*c1 - b2*c1 - a1*c2 + b1*c2;
   *delta  = -a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3;

   /* SCIPdebugMsg(scip, "alpha: %g beta: %g gamma: %g delta: %g\n", *alpha, *beta, *gamma_, *delta); */

   if( SCIPisInfinity(scip, REALABS(*gamma_ * a3)) ||
      SCIPisInfinity(scip, REALABS(*gamma_ * b3)) ||
      SCIPisInfinity(scip, REALABS(*gamma_ * c3)) )
   {
      SCIPdebugMsg(scip, "activity above SCIP infinity\n");
      *delta  = 0.0;
      *alpha  = 0.0;
      *beta   = 0.0;
      *gamma_ = 0.0;
      return SCIP_OKAY;
   }

   /* check if hyperplane contains all three points (necessary because of numerical troubles) */
   if( !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3) ||
      !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3) ||
      !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3) )
   {
      SCIP_Real m[9];
      SCIP_Real rhs[3];
      SCIP_Real x[3];
      SCIP_Bool success;

      /*
      SCIPdebugMsg(scip, "a = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", a1, a2, a3, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3, SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3));
      SCIPdebugMsg(scip, "b = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", b1, b2, b3, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3, SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3));
      SCIPdebugMsg(scip, "c = (%g,%g,%g) hyperplane: %g rhs %g EQdelta: %d\n", c1, c2, c3, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3, SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3));
      */

      /* initialize matrix column-wise */
      m[0] = a1;
      m[1] = b1;
      m[2] = c1;
      m[3] = a2;
      m[4] = b2;
      m[5] = c2;
      m[6] = a3;
      m[7] = b3;
      m[8] = c3;

      rhs[0] = 1.0;
      rhs[1] = 1.0;
      rhs[2] = 1.0;

      SCIPdebugMsg(scip, "numerical troubles - try to solve the linear system via an LU factorization\n");

      /* solve the linear problem */
      SCIP_CALL( SCIPsolveLinearEquationsIpopt(3, m, rhs, x, &success) );

      *delta  = rhs[0];
      *alpha  = x[0];
      *beta   = x[1];
      *gamma_ = x[2];

      /* set all coefficients to zero if one of the points is not contained in the hyperplane; this ensures that we do
       * not add a cut to SCIP and that all assertions are trivially fulfilled
       */
      if( !success || !SCIPisRelEQ(scip, *alpha * a1 + *beta * a2 - *delta, -*gamma_ * a3) ||
         !SCIPisRelEQ(scip, *alpha * b1 + *beta * b2 - *delta, -*gamma_ * b3) ||
         !SCIPisRelEQ(scip, *alpha * c1 + *beta * c2 - *delta, -*gamma_ * c3) ) /*lint !e774*/
      {
         SCIPdebugMsg(scip, "could not resolve numerical difficulties\n");
         *delta  = 0.0;
         *alpha  = 0.0;
         *beta   = 0.0;
         *gamma_ = 0.0;
      }
   }

   if( *gamma_ < 0.0 )
   {
      *alpha  = -*alpha;
      *beta   = -*beta;
      *gamma_ = -*gamma_;
      *delta  = -*delta;
   }

   return SCIP_OKAY;
}

/** computes a facet of the convex or concave envelope of a bivariate vertex polyhedral function */
static
SCIP_RETCODE computeVertexPolyhedralFacetBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_Real             p1[2],              /**< first vertex of box */
   SCIP_Real             p2[2],              /**< second vertex of box */
   SCIP_Real             p3[2],              /**< third vertex of box */
   SCIP_Real             p4[2],              /**< forth vertex of box */
   SCIP_Real             p1val,              /**< value in p1 */
   SCIP_Real             p2val,              /**< value in p2 */
   SCIP_Real             p3val,              /**< value in p3 */
   SCIP_Real             p4val,              /**< value in p4 */
   SCIP_Real             xstar[2],           /**< point to be separated */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least 2 */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
   )
{
   SCIP_Real alpha, beta, gamma_, delta;
   SCIP_Real xstarval, candxstarval = 0.0;
   int leaveout;

   assert(scip != NULL);
   assert(success != NULL);
   assert(SCIPisFinite(p1val) && p1val != SCIP_INVALID);
   assert(SCIPisFinite(p2val) && p2val != SCIP_INVALID);
   assert(SCIPisFinite(p3val) && p3val != SCIP_INVALID);
   assert(SCIPisFinite(p4val) && p4val != SCIP_INVALID);
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   /* if we want an underestimator, flip f(x,y), i.e., do as if we compute an overestimator for -f(x,y) */
   if( !overestimate )
   {
      p1val = -p1val;
      p2val = -p2val;
      p3val = -p3val;
      p4val = -p4val;
      targetvalue = -targetvalue;
   }

   SCIPdebugMsg(scip, "p1 = (%g, %g), f(p1) = %g\n", p1[0], p1[1], p1val);
   SCIPdebugMsg(scip, "p2 = (%g, %g), f(p2) = %g\n", p2[0], p2[1], p2val);
   SCIPdebugMsg(scip, "p3 = (%g, %g), f(p3) = %g\n", p3[0], p3[1], p3val);
   SCIPdebugMsg(scip, "p4 = (%g, %g), f(p4) = %g\n", p4[0], p4[1], p4val);

   /* Compute coefficients alpha, beta, gamma (>0), delta such that
    *   alpha*x + beta*y + gamma*z = delta
    * is satisfied by at least three of the corner points (p1,f(p1)), ..., (p4,f(p4)) and
    * the fourth corner point lies below this hyperplane.
    * Since we assume that f is vertex-polyhedral, we then know that all points (x,y,f(x,y)) are below this hyperplane, i.e.,
    *    alpha*x + beta*y - delta <= -gamma * f(x,y),
    * or, equivalently,
    *   -alpha/gamma*x - beta/gamma*y + delta/gamma >= f(x,y).
    */
   for( leaveout = 1; leaveout <= 4; ++leaveout )
   {
      switch( leaveout)
      {
         case 1 :
            /* get hyperplane through p2, p3, p4 */
            SCIP_CALL( computeHyperplaneThreePoints(scip, p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p1, then go to next candidate */
            if( alpha * p1[0] + beta * p1[1] + gamma_ * p1val - delta > 0.0 )
               continue;
            break;

         case 2 :
            /* get hyperplane through p1, p3, p4 */
            SCIP_CALL( computeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p2, then go to next candidate */
            if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val - delta > 0.0 )
               continue;
            break;

         case 3 :
            /* get hyperplane through p1, p2, p4 */
            SCIP_CALL( computeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1], p4val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p3, then go to next candidate */
            if( alpha * p3[0] + beta * p3[1] + gamma_ * p3val - delta > 0.0 )
               continue;
            break;

         case 4 :
            /* get hyperplane through p1, p2, p3 */
            SCIP_CALL( computeHyperplaneThreePoints(scip, p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val,
               &alpha, &beta, &gamma_, &delta) );
            /* if not underestimating in p4, then stop */
            if( alpha * p4[0] + beta * p4[1] + gamma_ * p4val - delta > 0.0 )
               continue;
            break;

         default: /* only for lint */
            alpha = SCIP_INVALID;
            beta = SCIP_INVALID;
            gamma_ =  SCIP_INVALID;
            delta = SCIP_INVALID;
            break;
      }

      /* check if bad luck: should not happen if numerics are fine */
      if( SCIPisZero(scip, gamma_) )
         continue;
      assert(!SCIPisNegative(scip, gamma_));

      /* if coefficients become tiny because division by gamma makes them < SCIPepsilon(scip), then skip, too */
      if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
         ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta/gamma_)) )
         continue;

      SCIPdebugMsg(scip, "alpha = %g, beta = %g, gamma = %g, delta = %g\n", alpha, beta, gamma_, delta);

      /* value of hyperplane candidate in xstar */
      xstarval = -alpha/gamma_ * xstar[0] -beta/gamma_ * xstar[1] + delta/gamma_;

      /* if reaching target and first or better than previous candidate, then update */
      if( xstarval <= targetvalue && (!*success || xstarval < candxstarval) )
      {
         /* flip hyperplane */
         if( !overestimate )
            gamma_ = -gamma_;

         facetcoefs[0] = -alpha/gamma_;
         facetcoefs[1] = -beta/gamma_;
         *facetconstant = delta/gamma_;

         *success = TRUE;
         candxstarval = xstarval;
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLR*     targetconshdlr;
   SCIP_CONSHDLRDATA* sourceconshdlrdata;
   int                i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(valid != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* create basic data of constraint handler and include it to scip */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );

   targetconshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(targetconshdlr != NULL);
   assert(targetconshdlr != conshdlr);

   sourceconshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(sourceconshdlrdata != NULL);

   /* copy nonlinear handlers */
   for( i = 0; i < sourceconshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CALL( SCIPnlhdlrCopyhdlr(scip, targetconshdlr, conshdlr, sourceconshdlrdata->nlhdlrs[i]) );
   }

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free nonlinear handlers */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CALL( SCIPnlhdlrFree(scip, &conshdlrdata->nlhdlrs[i]) );
      assert(conshdlrdata->nlhdlrs[i] == NULL);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->nlhdlrs, conshdlrdata->nlhdlrssize);
   conshdlrdata->nlhdlrssize = 0;

   /* free upgrade functions */
   for( i = 0; i < conshdlrdata->nconsupgrades; ++i )
   {
      assert(conshdlrdata->consupgrades[i] != NULL);
      SCIPfreeBlockMemory(scip, &conshdlrdata->consupgrades[i]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->consupgrades, conshdlrdata->consupgradessize);

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->canonicalizetime) );

   SCIPqueueFree(&conshdlrdata->reversepropqueue);

   if( conshdlrdata->vp_randnumgen != NULL )
      SCIPfreeRandom(scip, &conshdlrdata->vp_randnumgen);

   /* free LPs used to construct facets of envelops of vertex-polyhedral functions */
   for( i = 0; i <= SCIP_MAXVERTEXPOLYDIM; ++i )
   {
      if( conshdlrdata->vp_lp[i] != NULL )
      {
         SCIP_CALL( SCIPlpiFree(&conshdlrdata->vp_lp[i]) );
      }
   }

   assert(conshdlrdata->branchrandnumgen == NULL);

   assert(SCIPhashmapGetNElements(conshdlrdata->var2expr) == 0);
   SCIPhashmapFree(&conshdlrdata->var2expr);

   SCIPfreeBlockMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* make sure current activity tags in expressions are invalid, because we start catching variable events only now */
   conshdlrdata->lastboundrelax = ++conshdlrdata->curboundstag;
   /* set to 1 so it is larger than initial value of lastenforound in exprs */
   conshdlrdata->enforound = 1;
   /* reset numbering for auxiliary variables */
   conshdlrdata->auxvarid = 0;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( storeVarExprs(scip, conshdlr, SCIPconsGetData(conss[i])) );
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[i]) );
   }

   /* sort nonlinear handlers by detection priority, in decreasing order */
   if( conshdlrdata->nnlhdlrs > 1 )
      SCIPsortDownPtr((void**)conshdlrdata->nlhdlrs, SCIPnlhdlrComp, conshdlrdata->nnlhdlrs);

   /* get heuristics for later use */
   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   /* reset statistics in nonlinear handlers (TODO only if misc/resetstat == TRUE) and call nlhdlrInit */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CALL( SCIPnlhdlrInit(scip, conshdlrdata->nlhdlrs[i]) );
   }

   /* reset statistics in constraint handler */
   conshdlrdata->nweaksepa = 0;
   conshdlrdata->ntightenlp = 0;
   conshdlrdata->ndesperatebranch = 0;
   conshdlrdata->ndesperatecutoff = 0;
   conshdlrdata->ndesperatetightenlp = 0;
   conshdlrdata->nforcelp = 0;
   SCIP_CALL( SCIPresetClock(scip, conshdlrdata->canonicalizetime) );
   conshdlrdata->ncanonicalizecalls = 0;

#ifdef ENFOLOGFILE
   ENFOLOG( enfologfile = fopen(ENFOLOGFILE, "w"); )
#endif

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS** consssorted;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( nconss > 0 )
   {
      /* for better performance of dropVarEvents, we sort by index, descending */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consssorted, conss, nconss) );
      SCIPsortDownPtr((void**)consssorted, compIndexConsNonlinear, nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, consssorted[i]) );
         SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(consssorted[i])) );
      }

      SCIPfreeBufferArray(scip, &consssorted);
   }

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   if( conshdlrdata->vp_randnumgen != NULL )
      SCIPfreeRandom(scip, &conshdlrdata->vp_randnumgen);

   /* free LPs used to construct facets of envelops of vertex-polyhedral functions */
   for( i = 0; i <= SCIP_MAXVERTEXPOLYDIM; ++i )
   {
      if( conshdlrdata->vp_lp[i] != NULL )
      {
         SCIP_CALL( SCIPlpiFree(&conshdlrdata->vp_lp[i]) );
      }
   }

   if( conshdlrdata->branchrandnumgen != NULL )
      SCIPfreeRandom(scip, &conshdlrdata->branchrandnumgen);

   /* deinitialize nonlinear handlers */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_CALL( SCIPnlhdlrExit(scip, conshdlrdata->nlhdlrs[i]) );
   }

   ENFOLOG(
   if( enfologfile != NULL )
   {
      fclose(enfologfile);
      enfologfile = NULL;
   })

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_CONSINITPRE(consInitpreNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreNonlinear NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreNonlinear)
{  /*lint --e{715}*/
   SCIP_Bool infeasible;

   if( nconss == 0 )
      return SCIP_OKAY;

   /* skip some extra work if already known to be infeasible */
   if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE )
      return SCIP_OKAY;

   /* simplify constraints and replace common subexpressions */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, nconss, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );

   /* currently SCIP does not offer to communicate this,
    * but at the moment this can only become true if canonicalizeConstraints called detectNlhdlrs (which it doesn't do in EXITPRESOLVE stage)
    * or if a constraint expression became constant
    */
   assert(!infeasible);

   /* tell SCIP that we have something nonlinear */
   SCIPenableNLP(scip);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   /* skip remaining initializations if we have solved already
    * if infeasibility was found by our boundtightening, then curvature check may also fail as some exprhdlr (e.g., pow)
    * assumes nonempty activities in expressions
    */
   switch( SCIPgetStatus(scip) )
   {
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_INFEASIBLE:
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_INFORUNBD:
         return SCIP_OKAY;
      default: ;
   } /*lint !e788 */

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* reset one of the number of detections counter to count only current round */
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
      SCIPnlhdlrResetNDetectionslast(conshdlrdata->nlhdlrs[i]);

   SCIP_CALL( initSolve(scip, conshdlr, conss, nconss) );

   /* catch new solution event */
   if( nconss != 0 && conshdlrdata->linearizeheursol != 'o' )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME "_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPcatchEvent(scip, conshdlrdata->linearizeheursol == 'i' ? SCIP_EVENTTYPE_BESTSOLFOUND : SCIP_EVENTTYPE_SOLFOUND,
         eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   /* check that branching/lpgainnormalize is set to a known value if pseudo-costs are used in branching */
   if( conshdlrdata->branchpscostweight > 0.0 )
   {
      SCIP_CALL( SCIPgetCharParam(scip, "branching/lpgainnormalize", &(conshdlrdata->branchpscostupdatestrategy)) );
      if( strchr("lds", conshdlrdata->branchpscostupdatestrategy) == NULL )
      {
         SCIPerrorMessage("branching/lpgainnormalize strategy %c unknown\n", conshdlrdata->branchpscostupdatestrategy);
         SCIPABORT();
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( deinitSolve(scip, conshdlr, conss, nconss) );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free hash table for bilinear terms */
   SCIP_CALL( bilinearTermsFree(scip, conshdlrdata) );

   /* reset flag to allow another call of presolSingleLockedVars() after a restart */
   conshdlrdata->checkedvarlocks = FALSE;

   /* drop catching new solution event, if catched before */
   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME "_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, conshdlrdata->linearizeheursol == 'i' ? SCIP_EVENTTYPE_BESTSOLFOUND : SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteNonlinear)
{  /*lint --e{715}*/
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->expr != NULL);

   /* constraint locks should have been removed */
   assert((*consdata)->nlockspos == 0);
   assert((*consdata)->nlocksneg == 0);

   /* free variable expressions */
   SCIP_CALL( freeVarExprs(scip, *consdata) );

   SCIP_CALL( SCIPreleaseExpr(scip, &(*consdata)->expr) );

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransNonlinear)
{  /*lint --e{715}*/
   SCIP_EXPR* targetexpr;
   SCIP_CONSDATA* sourcedata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( SCIPduplicateExpr(scip, sourcedata->expr, &targetexpr, mapexprtransvar, conshdlr, exprownerCreate, (void*)conshdlr) );
   assert(targetexpr != NULL);  /* SCIPduplicateExpr cannot fail */

   /* create transformed cons (only captures targetexpr, no need to copy again) */
   SCIP_CALL( createCons(scip, conshdlr, targetcons, SCIPconsGetName(sourcecons),
      targetexpr, sourcedata->lhs, sourcedata->rhs, FALSE,
      SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
      SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
      SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
      SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons)) );

   /* release target expr */
   SCIP_CALL( SCIPreleaseExpr(scip, &targetexpr) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpNonlinear)
{  /*lint --e{715}*/
   /* create auxiliary variables and call separation initialization callbacks of the expression handlers
    * TODO if we ever want to allow constraints that are separated but not initial, then we need to call initSepa also
    *   during SEPALP, ENFOLP, etc, whenever a constraint may be separated the first time
    *   for now, there is an assert in detectNlhdlrs to require initial if separated
    */
   SCIP_CALL( initSepa(scip, conshdlr, conss, nconss, infeasible) );

   /* collect all bilinear terms for which an auxvar is present
    * TODO this will only do something for the first call of initlp after initsol, because it cannot handle
    * addition (and removal?) of constraints during solve
    * this is typically the majority of constraints, but the method should be made more flexible
    */
   SCIP_CALL( bilinearTermsInsertAll(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpNonlinear)
{  /*lint --e{715}*/
   SCIP_CALL( consSepa(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolNonlinear)
{  /*lint --e{715}*/
   SCIP_CALL( consSepa(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpNonlinear)
{  /*lint --e{715}*/
   SCIP_CALL( consEnfo(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxNonlinear)
{  /*lint --e{715}*/
   SCIP_CALL( consEnfo(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsNonlinear)
{  /*lint --e{715}*/
   SCIP_RESULT propresult;
   SCIP_Longint soltag;
   int nchgbds;
   int nnotify;
   int c;

   soltag = SCIPgetExprNewSoltag(scip);

   *result = SCIP_FEASIBLE;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], NULL, soltag) );

      if( isConsViolated(scip, conss[c]) )
         *result = SCIP_INFEASIBLE;
   }

   if( *result == SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* try to propagate
    * TODO obey propinenfo parameter, but we need something to recognize cutoff
    */
   nchgbds = 0;
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, TRUE, &propresult, &nchgbds) );

   if( (propresult == SCIP_CUTOFF) || (propresult == SCIP_REDUCEDDOM) )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* register all unfixed variables in all violated constraints as branching candidates */
   SCIP_CALL( registerBranchingCandidatesAllUnfixed(scip, conshdlr, conss, nconss, &nnotify) );
   if( nnotify > 0 )
   {
      SCIPdebugMsg(scip, "registered %d external branching candidates\n", nnotify);

      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "could not find branching candidates, forcing to solve LP\n");
   *result = SCIP_SOLVELP;
   ++SCIPconshdlrGetData(conshdlr)->nforcelp;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   SCIP_Bool          maypropfeasible;
   SCIP_Longint       soltag;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;
   soltag = SCIPgetExprNewSoltag(scip);
   maxviol = 0.0;
   maypropfeasible = conshdlrdata->trysolheur != NULL && SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED
      && SCIPgetStage(scip) <= SCIP_STAGE_SOLVING;

   /* check nonlinear constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      if( isConsViolated(scip, conss[c]) )
      {
         *result = SCIP_INFEASIBLE;
         maxviol = MAX(maxviol, getConsAbsViolation(conss[c]));

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /* print reason for infeasibility */
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( consdata->lhsviol > SCIPfeastol(scip) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhsviol);
            }
            if( consdata->rhsviol > SCIPfeastol(scip) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", consdata->rhsviol);
            }
         }
         else if( (conshdlrdata->subnlpheur == NULL || sol == NULL) && !maypropfeasible && !completely )
         {
            /* if we don't want to pass to subnlp heuristic and don't need to print reasons, then can stop checking here */
            return SCIP_OKAY;
         }

         /* do not try to shift linear variables if violation is at infinity (leads to setting variable to infinity in solution, which is not allowed) */
         if( maypropfeasible && SCIPisInfinity(scip, getConsAbsViolation(conss[c])) )
            maypropfeasible = FALSE;

         if( maypropfeasible )
         {
            if( consdata->lhsviol > SCIPfeastol(scip) )
            {
               /* check if there is a variable which may help to get the left hand side satisfied
                * if there is no such variable, then we cannot get feasible
                */
               if( !(consdata->linvarincr != NULL && consdata->linvarincrcoef > 0.0) &&
                   !(consdata->linvardecr != NULL && consdata->linvardecrcoef < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(consdata->rhsviol > SCIPfeastol(scip));
               /* check if there is a variable which may help to get the right hand side satisfied
                * if there is no such variable, then we cannot get feasible
                */
               if( !(consdata->linvarincr != NULL && consdata->linvarincrcoef < 0.0) &&
                   !(consdata->linvardecr != NULL && consdata->linvardecrcoef > 0.0) )
                  maypropfeasible = FALSE;
            }
         }
      }
   }

   if( *result == SCIP_INFEASIBLE && maypropfeasible )
   {
      SCIP_Bool success;

      SCIP_CALL( proposeFeasibleSolution(scip, conshdlr, conss, nconss, sol, &success) );

      /* do not pass solution to NLP heuristic if we made it feasible this way */
      if( success )
         return SCIP_OKAY;
   }

   if( *result == SCIP_INFEASIBLE && conshdlrdata->subnlpheur != NULL && sol != NULL && !SCIPisInfinity(scip, maxviol) )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropNonlinear)
{  /*lint --e{715}*/
   int nchgbds = 0;

   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, FALSE, result, &nchgbds) );
   assert(nchgbds >= 0);

   /* TODO would it make sense to check for redundant constraints? */

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible;
   int c;

   *result = SCIP_DIDNOTFIND;

   if( nconss == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* simplify constraints and replace common subexpressions, reinit nlhdlrs */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, nconss, presoltiming, &infeasible, ndelconss, naddconss, nchgcoefs) );
   if( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* merge constraints with the same root expression */
   if( presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE )
   {
      SCIP_Bool success;

      SCIP_CALL( presolveMergeConss(scip, conss, nconss, &success) );
      if( success )
         *result = SCIP_SUCCESS;
   }

   /* propagate constraints */
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, FALSE, result, nchgbds) );
   if( *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   /* propagate function domains (TODO integrate with simplify?) */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) || nrounds == 0 )
   {
      SCIP_RESULT localresult;
      SCIP_CALL( propExprDomains(scip, conshdlr, conss, nconss, &localresult, nchgbds) );
      if( localresult == SCIP_CUTOFF )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( localresult == SCIP_REDUCEDDOM )
         *result = SCIP_REDUCEDDOM;
   }

   /* check for redundant constraints, remove constraints that are a value expression */
   SCIP_CALL( presolveRedundantConss(scip, conshdlr, conss, nconss, &infeasible, ndelconss, nchgbds) );
   if( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* try to upgrade constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_Bool upgraded;

      /* skip inactive and deleted constraints */
      if( SCIPconsIsDeleted(conss[c]) || !SCIPconsIsActive(conss[c]) )
         continue;

      SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss) );
   }

   /* try to change continuous variables that appear linearly to be implicit integer */
   if( presoltiming & SCIP_PRESOLTIMING_MEDIUM )
   {
      SCIP_CALL( presolveImplint(scip, conshdlr, conss, nconss, nchgvartypes, &infeasible) );

      if( infeasible )
      {
         SCIPdebugMsg(scip, "presolveImplint() detected infeasibility\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* fix variables that are contained in only one nonlinear constraint to their upper or lower bounds, if possible */
   if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) && SCIPisPresolveFinished(scip)
      && !conshdlrdata->checkedvarlocks && conshdlrdata->checkvarlocks != 'd' )
   {
      /* run this presolving technique only once because we don't want to generate identical bound disjunction
       * constraints multiple times
       */
      conshdlrdata->checkedvarlocks = TRUE;

      for( c = 0; c < nconss; ++c )
      {
         int tmpnchgvartypes = 0;
         int tmpnaddconss = 0;

         SCIP_CALL( presolveSingleLockedVars(scip, conshdlr, conss[c], &tmpnchgvartypes, &tmpnaddconss, &infeasible) );
         SCIPdebugMsg(scip, "presolSingleLockedVars() for %s: nchgvartypes=%d naddconss=%d infeas=%u\n",
            SCIPconsGetName(conss[c]), tmpnchgvartypes, tmpnaddconss, infeasible);

         if( infeasible )
         {
            SCIPdebugMsg(scip, "presolSingleLockedVars() detected infeasibility\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         (*nchgvartypes) += tmpnchgvartypes;
         (*naddconss) += tmpnaddconss;
      }
   }

   if( *ndelconss > 0 || *nchgbds > 0 || *nupgdconss > 0 || *naddconss > 0 || *nchgvartypes > 0 )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_CONSRESPROP(consRespropNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropNonlinear NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_Bool reinitsolve = FALSE;

   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   ownerdata = SCIPexprGetOwnerData(consdata->expr);

   /* check whether we need to initSolve again because
    * - we have enfo initialized (nenfos >= 0)
    * - and locks appeared (going from zero to nonzero) or disappeared (going from nonzero to zero) now
    */
   if( ownerdata->nenfos >= 0 )
   {
      if( (consdata->nlockspos == 0) != (nlockspos == 0) )
         reinitsolve = TRUE;
      if( (consdata->nlocksneg == 0) != (nlocksneg == 0) )
         reinitsolve = TRUE;
   }

   if( reinitsolve )
   {
      SCIP_CALL( deinitSolve(scip, conshdlr, &cons, 1) );
   }

   /* add locks */
   SCIP_CALL( addLocks(scip, cons, nlockspos, nlocksneg) );

   if( reinitsolve )
   {
      SCIP_CALL( initSolve(scip, conshdlr, &cons, 1) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* simplify root expression if the constraint has been added after presolving */
   if( SCIPgetStage(scip) > SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_Bool replacedroot;

      if( !consdata->issimplified )
      {
         SCIP_EXPR* simplified;
         SCIP_Bool changed;

         /* simplify constraint */
         SCIP_CALL( SCIPsimplifyExpr(scip, consdata->expr, &simplified, &changed, &infeasible, exprownerCreate, (void*)conshdlr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );
         assert(simplified != NULL);
         consdata->expr = simplified;
         consdata->issimplified = TRUE;
      }

      /* ensure each variable is represented by one variable expression only (need this for storeVarExprs() with simplified=TRUE below) */
      SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &consdata->expr, 1, &replacedroot) );
      assert(!replacedroot);  /* root expression cannot have been equal to one of its subexpressions */

      /* ensure that varexprs in consdata->expr are the one from var2expr hashmap */
      {
         SCIP_CONSHDLRDATA* conshdlrdata;
         SCIP_EXPRITER* it;
         SCIP_EXPR* expr;

         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);

         SCIP_CALL( SCIPcreateExpriter(scip, &it) );
         SCIP_CALL( SCIPexpriterInit(it, consdata->expr, SCIP_EXPRITER_DFS, FALSE) );
         SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITINGCHILD);
         for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
         {
            SCIP_EXPR* child;
            SCIP_EXPR* hashmapexpr;

            child = SCIPexpriterGetChildExprDFS(it);
            if( !SCIPisExprVar(scip, child) )
               continue;

            /* check which expression is stored in the hashmap for the var of child */
            hashmapexpr = (SCIP_EXPR*)SCIPhashmapGetImage(conshdlrdata->var2expr, SCIPgetVarExprVar(child));
            /* if a varexpr exists already in the hashmap, but it is child, then replace child by the one in the hashmap */
            if( hashmapexpr != NULL && hashmapexpr != child )
            {
               SCIP_CALL( SCIPreplaceExprChild(scip, expr, SCIPexpriterGetChildIdxDFS(it), hashmapexpr) );
            }
         }
         SCIPfreeExpriter(&it);
      }
   }

   /* store variable expressions */
   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );
   }

   /* add manually locks to constraints that are not checked for feasibility */
   if( !SCIPconsIsChecked(cons) )
   {
      assert(consdata->nlockspos == 0);
      assert(consdata->nlocksneg == 0);

      SCIP_CALL( addLocks(scip, cons, 1, 0) );
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_INITPRESOLVE && !infeasible )
   {
      SCIP_CALL( initSolve(scip, conshdlr, &cons, 1) );
   }

   /* TODO deal with infeasibility */
   assert(!infeasible);

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) < SCIP_STAGE_EXITSOLVE )
   {
      SCIP_CALL( deinitSolve(scip, conshdlr, &cons, 1) );
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(cons)) );
   }

   /* remove locks that have been added in consActiveExpr() */
   if( !SCIPconsIsChecked(cons) )
   {
      SCIP_CALL( addLocks(scip, cons, -1, 0) );

      assert(SCIPconsGetData(cons)->nlockspos == 0);
      assert(SCIPconsGetData(cons)->nlocksneg == 0);
   }

   return SCIP_OKAY;
}


/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}


/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** variable deletion of constraint handler */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_CONSDELVARS(consDelvarsNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsNonlinear NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   /* print left hand side for ranged constraints */
   if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);
   }

   /* print expression */
   SCIP_CALL( SCIPprintExpr(scip, consdata->expr, file) );

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* targetconshdlr;
   SCIP_EXPR* targetexpr = NULL;
   SCIP_CONSDATA* sourcedata;

   assert(cons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   targetconshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(targetconshdlr != NULL);

   SCIP_CALL( SCIPcopyExpr(sourcescip, scip, sourcedata->expr, &targetexpr, exprownerCreate, (void*)targetconshdlr, varmap, consmap, global, valid) );

   if( targetexpr == NULL )
      *valid = FALSE;

   *cons = NULL;
   if( *valid )
   {
      /* create copy (only capture targetexpr, no need to copy again) */
      SCIP_CALL( createCons(scip, targetconshdlr, cons, name != NULL ? name : SCIPconsGetName(sourcecons),
         targetexpr, sourcedata->lhs, sourcedata->rhs, FALSE,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
   }

   if( targetexpr != NULL )
   {
      /* release target expr */
      SCIP_CALL( SCIPreleaseExpr(scip, &targetexpr) );
   }

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseNonlinear)
{  /*lint --e{715}*/
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   const char* endptr;
   SCIP_EXPR* consexprtree;

   SCIPdebugMsg(scip, "cons_nonlinear::consparse parsing %s\n", str);

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = FALSE;

   /* return if string empty */
   if( !*str )
      return SCIP_OKAY;

   endptr = str;

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   /* parse constraint to get lhs, rhs, and expression in between (from cons_linear.c::consparse, but parsing whole string first, then getting expression) */

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, (char**)&endptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the beginning of the expression and not a left-hand-side */
         lhs = -SCIPinfinity(scip);
      }
      else
      {
         /* it was indeed a left-hand-side, so continue parsing after it */
         str = endptr + 2;

         /* ignore whitespace */
         while( isspace((unsigned char)*str) )
            ++str;
      }
   }

   SCIPdebugMsg(scip, "str should start at beginning of expr: %s\n", str);

   /* parse expression: so far we did not allocate memory, so can just return in case of readerror */
   SCIP_CALL( SCIPparseExpr(scip, &consexprtree, str, &str, exprownerCreate, (void*)conshdlr) );

   /* check for left or right hand side */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for free constraint */
   if( strncmp(str, "[free]", 6) == 0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIPerrorMessage("cannot have left hand side and [free] status \n");
         SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
         return SCIP_OKAY;
      }
      *success = TRUE;
   }
   else
   {
      switch( *str )
      {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, (char**)&endptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );
            return SCIP_OKAY;
      }
   }

   /* create constraint */
   SCIP_CALL( createCons(scip, conshdlr, cons, name,
      consexprtree, lhs, rhs, FALSE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );
   assert(*cons != NULL);

   SCIP_CALL( SCIPreleaseExpr(scip, &consexprtree) );

   SCIPdebugMsg(scip, "created nonlinear constraint: <%s>\n", SCIPconsGetName(*cons));

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* store variable expressions if not done so far */
   SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

   /* check whether array is too small in order to store all variables */
   if( varssize < consdata->nvarexprs )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      vars[i] = SCIPgetVarExprVar(consdata->varexprs[i]);
      assert(vars[i] != NULL);
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* store variable expressions if not done so far */
   SCIP_CALL( storeVarExprs(scip, conshdlr, consdata) );

   *nvars = consdata->nvarexprs;
   *success = TRUE;

   return SCIP_OKAY;
}

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsNonlinear NULL
#endif

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputNonlinear)
{ /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* print statistics for constraint handler */
   SCIPinfoMessage(scip, file, "Nonlinear Conshdlr : %10s %10s %10s %10s %10s %10s %10s\n", "WeakSepa", "TightenLP", "DespTghtLP", "DespBranch", "DespCutoff", "ForceLP", "CanonTime");
   SCIPinfoMessage(scip, file, "  enforce%-10s:", "");
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nweaksepa);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ntightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatetightenlp);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatebranch);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->ndesperatecutoff);
   SCIPinfoMessage(scip, file, " %10lld", conshdlrdata->nforcelp);
   SCIPinfoMessage(scip, file, "\n");
   SCIPinfoMessage(scip, file, "  presolve%-9s: %-65s", "", "");
   SCIPinfoMessage(scip, file, " %10.2f", SCIPgetClockTime(scip, conshdlrdata->canonicalizetime));
   SCIPinfoMessage(scip, file, "\n");

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputNlhdlr)
{ /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   /* skip nlhdlr table if there never were active nonlinear constraints */
   if( SCIPconshdlrGetMaxNActiveConss(conshdlr) == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* print statistics for nonlinear handlers */
   SCIPnlhdlrPrintStatistics(scip, conshdlrdata->nlhdlrs, conshdlrdata->nnlhdlrs, file);

   return SCIP_OKAY;
}

/** execution method of display nlhdlrs dialog */
static
SCIP_DECL_DIALOGEXEC(dialogExecDisplayNlhdlrs)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   /* add dialog to history of dialogs that have been executed */
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* display list of nonlinear handler */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " nonlinear handler  enabled  detectprio  enforceprio  description\n");
   SCIPdialogMessage(scip, NULL, " -----------------  -------  ----------  -----------  -----------\n");
   for( i = 0; i < conshdlrdata->nnlhdlrs; ++i )
   {
      SCIP_NLHDLR* nlhdlr = conshdlrdata->nlhdlrs[i];
      assert(nlhdlr != NULL);

      SCIPdialogMessage(scip, NULL, " %-17s ", SCIPnlhdlrGetName(nlhdlr));
      SCIPdialogMessage(scip, NULL, " %7s ", SCIPnlhdlrIsEnabled(nlhdlr) ? "yes" : "no");
      SCIPdialogMessage(scip, NULL, " %10d ", SCIPnlhdlrGetDetectPriority(nlhdlr));
      SCIPdialogMessage(scip, NULL, " %11d ", SCIPnlhdlrGetEnfoPriority(nlhdlr));
      SCIPdialogMessage(scip, NULL, " %s", SCIPnlhdlrGetDesc(nlhdlr));
      SCIPdialogMessage(scip, NULL, "\n");
   }
   SCIPdialogMessage(scip, NULL, "\n");

   /* next dialog will be root dialog again */
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/*
 * constraint handler specific interface methods
 */

/** creates the handler for nonlinear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_DIALOG* parentdialog;

   /* create nonlinear constraint handler data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->intevalvar = intEvalVarBoundTightening;
   conshdlrdata->curboundstag = 1;
   conshdlrdata->lastboundrelax = 1;
   conshdlrdata->curpropboundstag = 1;
   conshdlrdata->newsoleventfilterpos = -1;
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->canonicalizetime) );
   SCIP_CALL( SCIPqueueCreate(&conshdlrdata->reversepropqueue, 100, 2.0) );
   SCIP_CALL( SCIPhashmapCreate(&conshdlrdata->var2expr, SCIPblkmem(scip), 100) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyNonlinear,
         consFreeNonlinear, consInitNonlinear, consExitNonlinear,
         consInitpreNonlinear, consExitpreNonlinear, consInitsolNonlinear, consExitsolNonlinear,
         consDeleteNonlinear, consTransNonlinear, consInitlpNonlinear,
         consSepalpNonlinear, consSepasolNonlinear, consEnfolpNonlinear, consEnforelaxNonlinear, consEnfopsNonlinear, consCheckNonlinear,
         consPropNonlinear, consPresolNonlinear, consRespropNonlinear, consLockNonlinear,
         consActiveNonlinear, consDeactiveNonlinear,
         consEnableNonlinear, consDisableNonlinear, consDelvarsNonlinear,
         consPrintNonlinear, consCopyNonlinear, consParseNonlinear,
         consGetVarsNonlinear, consGetNVarsNonlinear, consGetDiveBdChgsNonlinear, conshdlrdata) );

   /* add nonlinear constraint handler parameters */
   /* TODO organize into more subcategories */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a set of constraints within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/propauxvars",
         "whether to check bounds of all auxiliary variable to seed reverse propagation",
         &conshdlrdata->propauxvars, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/varboundrelax",
         "strategy on how to relax variable bounds during bound tightening: relax (n)ot, relax by (a)bsolute value, relax always by a(b)solute value, relax by (r)relative value",
         &conshdlrdata->varboundrelax, TRUE, 'r', "nabr", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/varboundrelaxamount",
         "by how much to relax variable bounds during bound tightening if strategy 'a', 'b', or 'r'",
         &conshdlrdata->varboundrelaxamount, TRUE, SCIPepsilon(scip), 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/conssiderelaxamount",
         "by how much to relax constraint sides during bound tightening",
         &conshdlrdata->conssiderelaxamount, TRUE, SCIPepsilon(scip), 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/vpmaxperturb",
         "maximal relative perturbation of reference point when computing facet of envelope of vertex-polyhedral function (dim>2)",
         &conshdlrdata->vp_maxperturb, TRUE, VERTEXPOLY_MAXPERTURBATION, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/vpadjfacetthresh",
         "adjust computed facet of envelope of vertex-polyhedral function up to a violation of this value times LP feasibility tolerance",
         &conshdlrdata->vp_adjfacetthreshold, TRUE, VERTEXPOLY_ADJUSTFACETFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/vpdualsimplex",
         "whether to use dual simplex instead of primal simplex for LP that computes facet of vertex-polyhedral function",
         &conshdlrdata->vp_dualsimplex, TRUE, VERTEXPOLY_USEDUALSIMPLEX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/bilinmaxnauxexprs",
           "maximal number of auxiliary expressions per bilinear term",
           &conshdlrdata->bilinmaxnauxexprs, FALSE, BILIN_MAXNAUXEXPRS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprods",
         "whether to reformulate products of binary variables during presolving",
         &conshdlrdata->reformbinprods, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsand",
         "whether to use the AND constraint handler for reformulating binary products",
         &conshdlrdata->reformbinprodsand, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/reformbinprodsfac",
         "minimum number of terms to reformulate bilinear binary products by factorizing variables (<= 1: disabled)",
         &conshdlrdata->reformbinprodsfac, FALSE, 50, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forbidmultaggrnlvar",
         "whether to forbid multiaggregation of nonlinear variables",
         &conshdlrdata->forbidmultaggrnlvar, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/tightenlpfeastol",
         "whether to tighten LP feasibility tolerance during enforcement, if it seems useful",
         &conshdlrdata->tightenlpfeastol, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/propinenforce",
         "whether to (re)run propagation in enforcement",
         &conshdlrdata->propinenforce, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/weakcutthreshold",
         "threshold for when to regard a cut from an estimator as weak (lower values allow more weak cuts)",
         &conshdlrdata->weakcutthreshold, TRUE, 0.2, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/strongcutmaxcoef",
         "\"strong\" cuts will be scaled to have their maximal coef in [1/strongcutmaxcoef,strongcutmaxcoef]",
         &conshdlrdata->strongcutmaxcoef, TRUE, 1000.0, 1.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/strongcutefficacy",
         "consider efficacy requirement when deciding whether a cut is \"strong\"",
         &conshdlrdata->strongcutefficacy, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/forcestrongcut",
         "whether to force \"strong\" cuts in enforcement",
         &conshdlrdata->forcestrongcut, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/enfoauxviolfactor",
         "an expression will be enforced if the \"auxiliary\" violation is at least this factor times the \"original\" violation",
         &conshdlrdata->enfoauxviolfactor, TRUE, 0.01, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/weakcutminviolfactor",
         "retry enfo of constraint with weak cuts if violation is least this factor of maximal violated constraints",
         &conshdlrdata->weakcutminviolfactor, TRUE, 0.5, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/rownotremovable",
         "whether to make rows to be non-removable in the node where they are added (can prevent some cycling): 'o'ff, in 'e'nforcement only, 'a'lways",
         &conshdlrdata->rownotremovable, TRUE, 'o', "oea", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/violscale",
         "method how to scale violations to make them comparable (not used for feasibility check): (n)one, (a)ctivity and side, norm of (g)radient",
         &conshdlrdata->violscale, TRUE, 'n', "nag", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/checkvarlocks",
         "whether variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction)",
         &conshdlrdata->checkvarlocks, TRUE, 't', "bdt", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/branching/aux",
         "from which depth on in the tree to allow branching on auxiliary variables (variables added for extended formulation)",
         &conshdlrdata->branchauxmindepth, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/branching/external",
         "whether to use external branching candidates and branching rules for branching",
         &conshdlrdata->branchexternal, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/highviolfactor",
         "consider a constraint highly violated if its violation is >= this factor * maximal violation among all constraints",
         &conshdlrdata->branchhighviolfactor, FALSE, 0.0, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/highscorefactor",
         "consider a variable branching score high if its branching score >= this factor * maximal branching score among all variables",
         &conshdlrdata->branchhighscorefactor, FALSE, 0.9, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/violweight",
         "weight by how much to consider the violation assigned to a variable for its branching score",
         &conshdlrdata->branchviolweight, FALSE, 1.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/dualweight",
         "weight by how much to consider the dual values of rows that contain a variable for its branching score",
         &conshdlrdata->branchdualweight, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/pscostweight",
         "weight by how much to consider the pseudo cost of a variable for its branching score",
         &conshdlrdata->branchpscostweight, FALSE, 1.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/domainweight",
         "weight by how much to consider the domain width in branching score",
         &conshdlrdata->branchdomainweight, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/vartypeweight",
         "weight by how much to consider variable type (continuous: 0, binary: 1, integer: 0.1, impl-integer: 0.01) in branching score",
         &conshdlrdata->branchvartypeweight, FALSE, 0.5, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branching/scoreagg",
         "how to aggregate several branching scores given for the same expression: 'a'verage, 'm'aximum, 's'um",
         &conshdlrdata->branchscoreagg, FALSE, 's', "ams", NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branching/violsplit",
         "method used to split violation in expression onto variables: 'u'niform, 'm'idness of solution, 'd'omain width, 'l'ogarithmic domain width",
         &conshdlrdata->branchviolsplit, FALSE, 'm', "umdl", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/branching/pscostreliable",
         "minimum pseudo-cost update count required to consider pseudo-costs reliable",
         &conshdlrdata->branchpscostreliable, FALSE, 2.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/linearizeheursol",
         "whether tight linearizations of nonlinear constraints should be added to cutpool when some heuristics finds a new solution ('o'ff, on new 'i'ncumbents, on 'e'very solution)",
         &conshdlrdata->linearizeheursol, FALSE, 'o', "oie", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/assumeconvex",
         "whether to assume that any constraint is convex",
         &conshdlrdata->assumeconvex, FALSE, FALSE, NULL, NULL) );

   /* include handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, CONSHDLR_NAME "_boundchange",
         "signals a bound change to a nonlinear constraint", processVarEvent, NULL) );
   assert(conshdlrdata->eventhdlr != NULL);

   /* include tables for statistics */
   assert(SCIPfindTable(scip, TABLE_NAME_NONLINEAR) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_NONLINEAR, TABLE_DESC_NONLINEAR, FALSE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputNonlinear,
         NULL, TABLE_POSITION_NONLINEAR, TABLE_EARLIEST_STAGE_NONLINEAR) );

   assert(SCIPfindTable(scip, TABLE_NAME_NLHDLR) == NULL);
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_NLHDLR, TABLE_DESC_NLHDLR, TRUE,
         NULL, NULL, NULL, NULL, NULL, NULL, tableOutputNlhdlr,
         NULL, TABLE_POSITION_NLHDLR, TABLE_EARLIEST_STAGE_NLHDLR) );

   /* create, include, and release display nlhdlrs dialog */
   if( SCIPgetRootDialog(scip) != NULL && SCIPdialogFindEntry(SCIPgetRootDialog(scip), "display", &parentdialog) == 1 )
   {
      SCIP_DIALOG* dialog;

      assert(parentdialog != NULL);
      assert(!SCIPdialogHasEntry(parentdialog, DIALOG_NAME));

      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, dialogExecDisplayNlhdlrs, NULL, NULL,
            DIALOG_NAME, DIALOG_DESC, DIALOG_ISSUBMENU, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, CONSHDLR_NAME "_newsolution", "handles the event that a new primal solution has been found",
         processNewSolutionEvent, NULL) );

   return SCIP_OKAY;
}

/** includes a nonlinear constraint upgrade method into the nonlinear constraint handler */
SCIP_RETCODE SCIPincludeConsUpgradeNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_NONLINCONSUPGD((*nlconsupgd)),  /**< method to call for upgrading nonlinear constraint */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   CONSUPGRADE*       consupgrade;
   char               paramname[SCIP_MAXSTRLEN];
   char               paramdesc[SCIP_MAXSTRLEN];
   int                i;

   assert(conshdlrname != NULL );
   assert(nlconsupgd != NULL);

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether upgrade method exists already */
   for( i = conshdlrdata->nconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->consupgrades[i]->consupgd == nlconsupgd )
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage(scip, "Try to add already known upgrade method for constraint handler <%s>.\n", conshdlrname);
#endif
         return SCIP_OKAY;
      }
   }

   /* create a nonlinear constraint upgrade data object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consupgrade) );
   consupgrade->consupgd = nlconsupgd;
   consupgrade->priority = priority;
   consupgrade->active   = active;

   /* insert nonlinear constraint upgrade method into constraint handler data */
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &conshdlrdata->consupgrades, &conshdlrdata->consupgradessize, conshdlrdata->nconsupgrades+1) );
   assert(conshdlrdata->nconsupgrades+1 <= conshdlrdata->consupgradessize);

   for( i = conshdlrdata->nconsupgrades; i > 0 && conshdlrdata->consupgrades[i-1]->priority < consupgrade->priority; --i )
      conshdlrdata->consupgrades[i] = conshdlrdata->consupgrades[i-1];
   assert(0 <= i && i <= conshdlrdata->nconsupgrades);
   conshdlrdata->consupgrades[i] = consupgrade;
   conshdlrdata->nconsupgrades++;

   /* adds parameter to turn on and off the upgrading step */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/" CONSHDLR_NAME "/upgrade/%s", conshdlrname);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable nonlinear upgrading for constraint handler <%s>", conshdlrname);
   SCIP_CALL( SCIPaddBoolParam(scip,
         paramname, paramdesc,
         &consupgrade->active, FALSE, active, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsNonlinear() call, if you don't need all the information */
   SCIP_CONSHDLR* conshdlr;

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( createCons(scip, conshdlr, cons, name, expr, lhs, rhs, TRUE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint with all its constraint flags set to their default values
 *
 *  All flags can be set via SCIPconsSetFLAGNAME-methods.
 *
 *  @see SCIPcreateConsNonlinear() for information about the basic constraint flag configuration.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*            expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name, expr, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a quadratic nonlinear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsQuadraticNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_EXPR* expr;

   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadterms == 0 || (quadvars1 != NULL && quadvars2 != NULL && quadcoefs != NULL));

   /* get nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create quadratic expression */
   SCIP_CALL( SCIPcreateExprQuadratic(scip, &expr, nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoefs, exprownerCreate, (void*)conshdlr) );
   assert(expr != NULL);

   /* create nonlinear constraint */
   SCIP_CALL( createCons(scip, conshdlr, cons, name, expr, lhs, rhs, FALSE,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   /* release quadratic expression (captured by constraint now) */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** creates and captures a quadratic nonlinear constraint with all its constraint flags set to their default values
 *
 *  All flags can be set via SCIPconsSetFLAGNAME-methods.
 *
 *  @see SCIPcreateConsQuadraticNonlinear() for information about the basic constraint flag configuration.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicQuadraticNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs                 /**< right hand side of quadratic equation */
   )
{
   SCIP_CALL( SCIPcreateConsQuadraticNonlinear(scip, cons, name, nlinvars, linvars, lincoefs, nquadterms, quadvars1, quadvars2, quadcoefs, lhs, rhs,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint that is a second-order cone constraint with all its constraint flags set to their default values
 *
 * \f$\sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} \leq \alpha_{n+1}\, (x_{n+1}+\beta_{n+1})\f$
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSOCNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables on left hand side of constraint (n) */
   SCIP_VAR**            vars,               /**< array with variables on left hand side (x_i) */
   SCIP_Real*            coefs,              /**< array with coefficients of left hand side variables (alpha_i), or NULL if all 1.0 */
   SCIP_Real*            offsets,            /**< array with offsets of variables (beta_i), or NULL if all 0.0 */
   SCIP_Real             constant,           /**< constant on left hand side (gamma) */
   SCIP_VAR*             rhsvar,             /**< variable on right hand side of constraint (x_{n+1}) */
   SCIP_Real             rhscoeff,           /**< coefficient of variable on right hand side (alpha_{n+1}) */
   SCIP_Real             rhsoffset           /**< offset of variable on right hand side (beta_{n+1}) */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPR* lhssum;
   SCIP_EXPR* terms[2];
   SCIP_Real termcoefs[2];
   int i;

   assert(vars != NULL || nvars == 0);

   SCIP_CALL( SCIPcreateExprSum(scip, &lhssum, 0, NULL, NULL, constant, NULL, NULL) );  /* gamma */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_EXPR* varexpr;
      SCIP_EXPR* powexpr;

      SCIP_CALL( SCIPcreateExprVar(scip, &varexpr, vars[i], NULL, NULL) );   /* x_i */
      if( offsets != NULL && offsets[i] != 0.0 )
      {
         SCIP_EXPR* sum;
         SCIP_CALL( SCIPcreateExprSum(scip, &sum, 1, &varexpr, NULL, offsets[i], NULL, NULL) );  /* x_i + beta_i */
         SCIP_CALL( SCIPcreateExprPow(scip, &powexpr, sum, 2.0, NULL, NULL) );   /* (x_i + beta_i)^2 */
         SCIP_CALL( SCIPreleaseExpr(scip, &sum) );
      }
      else
      {
         SCIP_CALL( SCIPcreateExprPow(scip, &powexpr, varexpr, 2.0, NULL, NULL) );  /* x_i^2 */
      }

      SCIP_CALL( SCIPappendExprSumExpr(scip, lhssum, powexpr, coefs != NULL ? coefs[i]*coefs[i] : 1.0) );  /* + alpha_i^2 (x_i + beta_i)^2 */
      SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );
      SCIP_CALL( SCIPreleaseExpr(scip, &powexpr) );
   }

   SCIP_CALL( SCIPcreateExprPow(scip, &terms[0], lhssum, 0.5, NULL, NULL) );  /* sqrt(...) */
   SCIP_CALL( SCIPreleaseExpr(scip, &lhssum) );
   termcoefs[0] = 1.0;

   SCIP_CALL( SCIPcreateExprVar(scip, &terms[1], rhsvar, NULL, NULL) );  /* x_{n+1} */
   termcoefs[1] = -rhscoeff;

   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 2, terms, termcoefs, 0.0, NULL, NULL) );  /* sqrt(...) - alpha_{n+1}x_{n_1} */

   SCIP_CALL( SCIPreleaseExpr(scip, &terms[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &terms[0]) );

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, cons, name, expr, -SCIPinfinity(scip), rhscoeff * rhsoffset) );

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** creates and captures a signpower nonlinear constraint with all its constraint flags set to their default values
 *
 * \f$\textrm{lhs} \leq \textrm{sign}(x+a) |x+a|^n + c z \leq \textrm{rhs}\f$
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSignpowerNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             x,                  /**< nonlinear variable x in constraint */
   SCIP_VAR*             z,                  /**< linear variable z in constraint */
   SCIP_Real             exponent,           /**< exponent n of |x+offset|^n term in constraint */
   SCIP_Real             xoffset,            /**< offset in |x+offset|^n term in constraint */
   SCIP_Real             zcoef,              /**< coefficient of z in constraint */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_EXPR* xexpr;
   SCIP_EXPR* terms[2];
   SCIP_Real coefs[2];
   SCIP_EXPR* sumexpr;

   assert(x != NULL);
   assert(z != NULL);

   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   if( xoffset != 0.0 )
   {
      SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 1, &xexpr, NULL, xoffset, NULL, NULL) ); /* x + xoffset */
      SCIP_CALL( SCIPcreateExprSignpower(scip, &terms[0], sumexpr, exponent, NULL, NULL) ); /* signpow(x + xoffset, exponent) */

      SCIP_CALL( SCIPreleaseExpr(scip,  &sumexpr) );
   }
   else
   {
      SCIP_CALL( SCIPcreateExprSignpower(scip, &terms[0], xexpr, exponent, NULL, NULL) );  /* signpow(x, exponent) */
   }
   coefs[0] = 1.0;

   SCIP_CALL( SCIPcreateExprVar(scip, &terms[1], z, NULL, NULL) );
   coefs[1] = zcoef;

   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 2, terms, coefs, 0.0, NULL, NULL) );  /* signpowexpr + zcoef * z */

   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, cons, name, sumexpr, lhs, rhs) );

   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &terms[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &terms[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );

   return SCIP_OKAY;
}

/** gets tag indicating current local variable bounds */
SCIP_Longint SCIPgetCurBoundsTagNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->curboundstag;
}

/** gets the `curboundstag` from the last time where variable bounds were relaxed */
SCIP_Longint SCIPgetLastBoundRelaxTagNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->lastboundrelax;
}

/** increments `curboundstag` and resets `lastboundrelax` in constraint handler data
 *
 * @attention This method is not intended for normal use.
 *   These tags are maintained by the event handler for variable bound change events.
 *   This method is used by some unittests.
 */
void SCIPincrementCurBoundsTagNonlinear(
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_Bool             boundrelax          /**< indicates whether a bound was relaxed, i.e., lastboundrelax should be set too */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ++conshdlrdata->curboundstag;
   assert(conshdlrdata->curboundstag > 0);

   if( boundrelax )
      conshdlrdata->lastboundrelax = conshdlrdata->curboundstag;
}

/** returns the hashmap that is internally used to map variables to their corresponding variable expressions */
SCIP_HASHMAP* SCIPgetVarExprHashmapNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->var2expr;
}

/** processes a rowprep for cut addition and maybe report branchscores */
SCIP_RETCODE SCIPprocessRowprepNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler which provided the estimator */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_ROWPREP*         rowprep,            /**< cut to be added */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             allowweakcuts,      /**< whether we should only look for "strong" cuts, or anything that separates is fine */
   SCIP_Bool             branchscoresuccess, /**< whether the estimator generation generated branching scores */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, or only in separation */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_Real cutviol;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real auxvarvalue = SCIP_INVALID;
   SCIP_Bool sepasuccess;
   SCIP_Real estimateval = SCIP_INVALID;
   SCIP_Real mincutviolation;

   assert(nlhdlr != NULL);
   assert(cons != NULL);
   assert(expr != NULL);
   assert(rowprep != NULL);
   assert(auxvar != NULL);
   assert(result != NULL);

   /* decide on minimal violation of cut */
   if( sol == NULL )
      mincutviolation = SCIPgetLPFeastol(scip);  /* we enforce an LP solution */
   else
      mincutviolation = SCIPfeastol(scip);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   sepasuccess = TRUE;

   cutviol = SCIPgetRowprepViolation(scip, rowprep, sol, NULL);
   if( cutviol > 0.0 )
   {
      auxvarvalue = SCIPgetSolVal(scip, sol, auxvar);

      /* check whether cut is weak (if f(x) not defined, then it's never weak) */
      if( !allowweakcuts && auxvalue != SCIP_INVALID )
      {
         /* let the estimator be c'x-b, the auxvar is z (=auxvarvalue), and the expression is f(x) (=auxvalue)
          * then if we are underestimating and since the cut is violated, we should have z <= c'x-b <= f(x)
          * cutviol is c'x-b - z, so estimator value is c'x-b = z + cutviol
          * if the estimator value (c'x-b) is too close to z (auxvarvalue), when compared to f(x) (auxvalue),
          * then let's call this a weak cut that is, it's a weak cut if c'x-b <= z + weakcutthreshold * (f(x)-z)
          *   <->   c'x-b - z <= weakcutthreshold * (f(x)-z)
          *
          * if we are overestimating, we have z >= c'x-b >= f(x)
          * cutviol is z - (c'x-b), so estimator value is c'x-b = z - cutviol
          * it's weak if c'x-b >= f(x) + (1-weakcutthreshold) * (z - f(x))
          *   <->   c'x-b - z >= weakcutthreshold * (f(x)-z)
          *
          * when linearizing convex expressions, then we should have c'x-b = f(x), so they would never be weak
          */
         if( (!overestimate && ( cutviol <= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) ||
             ( overestimate && (-cutviol >= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded, but cut is too "\
                           "weak: auxvarvalue %g estimateval %g auxvalue %g (over %d)\n",
                                     SCIPnlhdlrGetName(nlhdlr), auxvarvalue,
                                     auxvarvalue + (overestimate ? -cutviol : cutviol), auxvalue, overestimate); )
            sepasuccess = FALSE;
         }
      }

      /* save estimator value for later, see long comment above why this gives the value for c'x-b */
      estimateval = auxvarvalue + (!overestimate ? cutviol : -cutviol);
   }
   else
   {
      sepasuccess = FALSE;
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded, but cut does not "\
                     "separate\n", SCIPnlhdlrGetName(nlhdlr)); )
   }

   /* clean up estimator */
   if( sepasuccess )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    estimate of nlhdlr %s succeeded: auxvarvalue %g "\
                     "estimateval %g auxvalue %g (over %d)\n    ", SCIPnlhdlrGetName(nlhdlr), auxvarvalue,
                               auxvarvalue + (overestimate ? -cutviol : cutviol), auxvalue, overestimate);
                       SCIPprintRowprep(scip, rowprep, enfologfile); )

      /* if not allowweakcuts, then do not attempt to get cuts more violated by scaling them up,
       * instead, may even scale them down, that is, scale so that max coef is close to 1
       */
      if( !allowweakcuts )
      {
         SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, conshdlrdata->strongcutmaxcoef, &sepasuccess) );

         if( !sepasuccess )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup cut failed due to bad numerics\n"); )
         }
         else
         {
            cutviol = SCIPgetRowprepViolation(scip, rowprep, sol, &sepasuccess);
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup succeeded, violation = %g and %sreliable, "\
                           "min requ viol = %g\n", cutviol, sepasuccess ? "" : "not ", mincutviolation); )
            if( sepasuccess )
               sepasuccess = cutviol > mincutviolation;
         }

         if( sepasuccess && auxvalue != SCIP_INVALID )
         {
            /* check whether cut is weak now
             * auxvar z may now have a coefficient due to scaling (down) in cleanup - take this into account when
             * reconstructing estimateval from cutviol (TODO improve or remove?)
             */
            SCIP_Real auxvarcoef = 0.0;
            int i;

            /* get absolute value of coef of auxvar in row - this makes the whole check here more expensive than
             * it should be...
             */
            for( i = 0; i < SCIProwprepGetNVars(rowprep); ++i )
            {
               if( SCIProwprepGetVars(rowprep)[i] == auxvar )
               {
                  auxvarcoef = REALABS(SCIProwprepGetCoefs(rowprep)[i]);
                  break;
               }
            }

            if( auxvarcoef == 0.0 ||
                (!overestimate && ( cutviol / auxvarcoef <= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) ||
                ( overestimate && (-cutviol / auxvarcoef >= conshdlrdata->weakcutthreshold * (auxvalue - auxvarvalue))) )
            {
               ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cut is too weak after cleanup: auxvarvalue %g estimateval %g auxvalue %g (over %d)\n",
                  auxvarvalue, auxvarvalue + (overestimate ? -cutviol : cutviol) / auxvarcoef, auxvalue, overestimate); )
               sepasuccess = FALSE;
            }
         }
      }
      else
      {
         /* TODO if violations are really tiny, then maybe handle special (decrease LP feastol, for example) */

         /* if estimate didn't report branchscores explicitly, then consider branching on those children for
          * which the following cleanup changes coefficients (we had/have this in expr_sum this way)
          */
         if( !branchscoresuccess )
            SCIProwprepRecordModifications(rowprep);

         SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, mincutviolation, &cutviol, &sepasuccess) );

         if( !sepasuccess )
         {
            ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cleanup failed, %d coefs modified, cutviol %g\n",
                                     SCIProwprepGetNModifiedVars(rowprep), cutviol); )
         }

         /* if cleanup left us with a useless cut, then consider branching on variables for which coef were
          * changed
          */
         if( !sepasuccess && !branchscoresuccess && SCIProwprepGetNModifiedVars(rowprep) > 0 )
         {
            SCIP_Real violscore;

#ifdef BRSCORE_ABSVIOL
            violscore = getExprAbsAuxViolation(scip, expr, auxvalue, sol, NULL, NULL);
#else
            SCIP_CALL( SCIPgetExprRelAuxViolationNonlinear(scip, expr, auxvalue, sol, &violscore, NULL, NULL) );
#endif
            SCIP_CALL( addExprViolScoresAuxVars(scip, expr, violscore, SCIProwprepGetModifiedVars(rowprep), SCIProwprepGetNModifiedVars(rowprep), sol, &branchscoresuccess) );

            /* addConsExprExprBranchScoresAuxVars can fail if the only vars for which the coef was changed
             * - were fixed,
             * - are this expr's auxvar (I don't think it makes sense to branch on that one (would it?)), or
             * - if a variable in the rowprep is not in expr (can happen with indicator added by perspective)
             * the first case came up again in #3085 and I don't see how to exclude this in the assert,
             * so I'm disabling the assert for now
             */
            /* assert(branchscoresuccess || (rowprep->nmodifiedvars == 1 && rowprep->modifiedvars[0] == auxvar) ||
                  strcmp(SCIPnlhdlrGetName(nlhdlr), "perspective")==0); */
         }
      }
   }

   /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
   if( sepasuccess )
   {
      SCIP_ROW* row;

      if( conshdlrdata->branchdualweight > 0.0 )
      {
         /* store remaining gap |f(x)-estimateval| in row name, which could be used in getDualBranchscore
          * skip if gap is zero
          */
         if( auxvalue == SCIP_INVALID )
            strcat(SCIProwprepGetName(rowprep), "_estimategap=inf");
         else if( !SCIPisEQ(scip, auxvalue, estimateval) )
         {
            char gap[40];
            (void) SCIPsnprintf(gap, 40, "_estimategap=%g", REALABS(auxvalue - estimateval));
            strcat(SCIProwprepGetName(rowprep), gap);
         }
      }

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

      if( !allowweakcuts && conshdlrdata->strongcutefficacy && !SCIPisCutEfficacious(scip, sol, row) )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cut efficacy %g is too low (minefficacy=%g)\n",
                                  SCIPgetCutEfficacy(scip, sol, row), SCIPgetSepaMinEfficacy(scip)); )
      }
      else if( !SCIPisCutApplicable(scip, row) )
      {
         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    cut not applicable (e.g., cut is boundchange below eps)\n"); )
      }
      else
      {
         SCIP_Bool infeasible;

         ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    adding cut ");
           SCIP_CALL( SCIPprintRow(scip, row, enfologfile) ); )

         /* I take !allowweakcuts as equivalent for having a strong cut (we usually have allowweakcuts=TRUE only
          * if we haven't found strong cuts before)
          */
         SCIP_CALL( SCIPaddRow(scip, row, conshdlrdata->forcestrongcut && !allowweakcuts && inenforcement, &infeasible) );

         /* mark row as not removable from LP for current node (this can prevent some cycling) */
         if( conshdlrdata->rownotremovable == 'a' || (conshdlrdata->rownotremovable == 'e' && inenforcement) )
            SCIPmarkRowNotRemovableLocal(scip, row);

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIPnlhdlrIncrementNCutoffs(nlhdlr);
         }
         else
         {
            *result = SCIP_SEPARATED;
            SCIPnlhdlrIncrementNSeparated(nlhdlr);
         }
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }
   else if( branchscoresuccess )
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    separation with estimate of nlhdlr %s failed, but "\
                     "branching candidates added\n", SCIPnlhdlrGetName(nlhdlr)); )

      /* well, not branched, but addConsExprExprViolScoresAuxVars() added scores to (aux)variables and that makes the
       * expressions eligible for branching candidate, see enforceConstraints() and branching()
       */
      *result = SCIP_BRANCHED;
   }
   else
   {
      ENFOLOG( SCIPinfoMessage(scip, enfologfile, "    separation with estimate of nlhdlr %s failed and no "\
                     "branching candidates%s\n", SCIPnlhdlrGetName(nlhdlr), (allowweakcuts && inenforcement) ?
                                                                                    " (!)" : ""); )
   }

   return SCIP_OKAY;
}

/** returns whether all nonlinear constraints are assumed to be convex */
SCIP_Bool SCIPassumeConvexNonlinear(
   SCIP_CONSHDLR*        conshdlr
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->assumeconvex;
}

/** collects all bilinear terms for a given set of constraints
 *
 * @attention This method should only be used for unit tests that depend on SCIPgetBilinTermsNonlinear(),
 *       SCIPgetBilinTermNonlinear() or SCIPgetBilinTermIdxNonlinear().
 */
SCIP_RETCODE SCIPcollectBilinTermsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_CONS**           conss,              /**< nonlinear constraints */
   int                   nconss              /**< total number of nonlinear constraints */
   )
{
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( bilinearTermsInsertAll(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}

/** returns the total number of bilinear terms that are contained in all nonlinear constraints
 *
 *  @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 */
int SCIPgetNBilinTermsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nbilinterms;
}

/** returns all bilinear terms that are contained in all nonlinear constraints
 *
 * @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @note The value of the auxiliary variable of a bilinear term might be NULL, which indicates that the term does not have an auxiliary variable.
 */
SCIP_CONSNONLINEAR_BILINTERM* SCIPgetBilinTermsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->bilinterms;
}

/** returns the index of the bilinear term representing the product of the two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns -1 if the variables do not appear bilinearly.
 */
int SCIPgetBilinTermIdxNonlinear(
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y                   /**< second variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM entry;
   int idx;

   assert(conshdlr != NULL);
   assert(x != NULL);
   assert(y != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->bilinhashtable == NULL )
   {
      return -1;
   }

   /* ensure that x.index <= y.index */
   if( SCIPvarCompare(x, y) == 1 )
   {
      SCIPswapPointers((void**)&x, (void**)&y);
   }
   assert(SCIPvarCompare(x, y) < 1);

   /* use a new entry to find the image in the bilinear hash table */
   entry.x = x;
   entry.y = y;
   idx = (int)(size_t)SCIPhashtableRetrieve(conshdlrdata->bilinhashtable, (void*)&entry) - 1;
   assert(idx >= -1 && idx < conshdlrdata->nbilinterms);
   assert(idx < 0 || conshdlrdata->bilinterms[idx].x == x);
   assert(idx < 0 || conshdlrdata->bilinterms[idx].y == y);

   return idx;
}

/** returns the bilinear term that represents the product of two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns NULL if the variables do not appear bilinearly.
 */
SCIP_CONSNONLINEAR_BILINTERM* SCIPgetBilinTermNonlinear(
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y                   /**< second variable */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int idx;

   assert(conshdlr != NULL);
   assert(x != NULL);
   assert(y != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   idx = SCIPgetBilinTermIdxNonlinear(conshdlr, x, y);
   assert(idx >= -1 && idx < conshdlrdata->nbilinterms);

   if( idx >= 0 )
   {
      return &conshdlrdata->bilinterms[idx];
   }

   return NULL;
}

/** evaluates an auxiliary expression for a bilinear term */
SCIP_Real SCIPevalBilinAuxExprNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable of the bilinear term */
   SCIP_VAR*             y,                  /**< second variable of the bilinear term */
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr,      /**< auxiliary expression */
   SCIP_SOL*             sol                 /**< solution at which to evaluate (can be NULL) */
   )
{
   assert(scip != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(auxexpr != NULL);
   assert(auxexpr->auxvar != NULL);

   return auxexpr->cst + auxexpr->coefs[0] * SCIPgetSolVal(scip, sol, auxexpr->auxvar) +
          auxexpr->coefs[1] * SCIPgetSolVal(scip, sol, x) + auxexpr->coefs[2] * SCIPgetSolVal(scip, sol, y);
}

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_RETCODE SCIPinsertBilinearTermExistingNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   int                   nlockspos,          /**< number of positive expression locks */
   int                   nlocksneg           /**< number of negative expression locks */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM* term;
   int idx;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( bilinearTermsInsertEntry(scip, conshdlr, x, y, nlockspos, nlocksneg, &idx, TRUE) );

   term = &conshdlrdata->bilinterms[idx];
   assert(term != NULL);
   assert(term->nauxexprs == 0);  /* existing terms should be added before implicit terms */
   assert(term->aux.var == NULL); /* there should not already be an auxvar, that is, existing terms should exist only once (common subexprs should have been eliminated) */

   /* store and capture auxiliary variable */
   if( auxvar != NULL )
   {
      term->aux.var = auxvar;
      SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   }

   return SCIP_OKAY;
}

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_RETCODE SCIPinsertBilinearTermImplicitNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   SCIP_Real             coefx,              /**< coefficient of x in the auxiliary expression */
   SCIP_Real             coefy,              /**< coefficient of y in the auxiliary expression */
   SCIP_Real             coefaux,            /**< coefficient of auxvar in the auxiliary expression */
   SCIP_Real             cst,                /**< constant of the auxiliary expression */
   SCIP_Bool             overestimate        /**< whether the auxiliary expression overestimates the bilinear product */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSNONLINEAR_BILINTERM* term;
   SCIP_CONSNONLINEAR_AUXEXPR* auxexpr;
   int idx;
   int nlockspos;
   int nlocksneg;
   SCIP_Bool added;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nlockspos = overestimate ? 1 : 0;
   nlocksneg = overestimate ? 0 : 1;

   SCIP_CALL( bilinearTermsInsertEntry(scip, conshdlr, x, y, nlockspos, nlocksneg, &idx, FALSE) );

   term = &conshdlrdata->bilinterms[idx];
   assert(term != NULL);
   assert(SCIPvarCompare(term->x, term->y) < 1);

   if( term->existing && term->nauxexprs == 0 && term->aux.var != NULL )
   {
      SCIP_CONSNONLINEAR_AUXEXPR* auxvarexpr;
      /* this is the case where we are adding an implicitly defined relation for a product that has already
       * been explicitly defined; convert auxvar into an auxexpr */

      /* nothing to do if we aren't allowed to add more than one auxexpr per term */
      if( conshdlrdata->bilinmaxnauxexprs <= 1 )
         return SCIP_OKAY;

      SCIP_CALL( SCIPallocBlockMemory(scip, &auxvarexpr) );
      auxvarexpr->cst = 0.0;
      auxvarexpr->coefs[0] = 1.0;
      auxvarexpr->coefs[1] = 0.0;
      auxvarexpr->coefs[2] = 0.0;
      auxvarexpr->auxvar = term->aux.var;
      auxvarexpr->underestimate = term->nlocksneg > 0;
      auxvarexpr->overestimate = term->nlockspos > 0;

      /* before we were working with term->aux.var; now aux.var has been saved and aux.exprs can be initialised to NULL */
      term->aux.exprs = NULL;

      SCIP_CALL( bilinTermAddAuxExpr(scip, conshdlrdata, term, auxvarexpr, &added) );

      /* since there were no auxexprs before and we've already checked for bilinmaxnauxexprs, auxvarexpr should always be added */
      assert(added);
   }

   /* create and add auxexpr */
   SCIP_CALL( SCIPallocBlockMemory(scip, &auxexpr) );
   auxexpr->underestimate = !overestimate;
   auxexpr->overestimate = overestimate;
   auxexpr->auxvar = auxvar;
   auxexpr->coefs[0] = coefaux;
   if( term->x == x )
   {
      assert(term->y == y);
      auxexpr->coefs[1] = coefx;
      auxexpr->coefs[2] = coefy;
   }
   else
   {
      assert(term->x == y);
      assert(term->y == x);
      auxexpr->coefs[1] = coefy;
      auxexpr->coefs[2] = coefx;
   }
   auxexpr->cst = cst;
   SCIP_CALL( bilinTermAddAuxExpr(scip, conshdlrdata, term, auxexpr, &added) );

   if( !added )
   {
      SCIPfreeBlockMemory(scip, &auxexpr);
   }
   else if( auxvar != NULL )
   { /* capture auxiliary variable */
      SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   }

   return SCIP_OKAY;
}

/* replication of long comment on SCIPcomputeFacetVertexPolyhedralNonlinear() in cons_nonlinear.h omitted here */
SCIP_RETCODE SCIPcomputeFacetVertexPolyhedralNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_DECL_VERTEXPOLYFUN((*function)),     /**< pointer to vertex polyhedral function */
   void*                 fundata,            /**< data for function evaluation (can be NULL) */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
   )
{
   SCIP_Real* corner;
   SCIP_Real* funvals;
   int* nonfixedpos;
   SCIP_Real maxfaceterror;
   int nvars; /* number of nonfixed variables */
   unsigned int ncorners;
   unsigned int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(function != NULL);
   assert(xstar != NULL);
   assert(box != NULL);
   assert(success != NULL);
   assert(facetcoefs != NULL);
   assert(facetconstant != NULL);

   *success = FALSE;

   /* identify fixed variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &nonfixedpos, nallvars) );
   nvars = 0;
   for( j = 0; j < nallvars; ++j )
   {
      if( SCIPisRelEQ(scip, box[2 * j], box[2 * j + 1]) )
         continue;
      nonfixedpos[nvars] = j;
      nvars++;
   }

   /* if all variables are fixed, then we could provide something trivial, but that wouldn't be the job of separation
    * if too many variables are not fixed, then we do nothing currently
    */
   if( nvars == 0 || nvars > SCIP_MAXVERTEXPOLYDIM )
   {
      SCIPwarningMessage(scip, "SCIPcomputeFacetVertexPolyhedralNonlinear() called with %d nonfixed variables. Must be between [1,%d].\n", nvars, SCIP_MAXVERTEXPOLYDIM);
      SCIPfreeBufferArray(scip, &nonfixedpos);
      return SCIP_OKAY;
   }

   /* compute f(v^i) for each corner v^i of [l,u] */
   ncorners = POWEROFTWO(nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &funvals, ncorners) );
   SCIP_CALL( SCIPallocBufferArray(scip, &corner, nallvars) );
   for( j = 0; j < nallvars; ++j )
   {
      if( SCIPisRelEQ(scip, box[2 * j], box[2 * j + 1]) )
         corner[j] = (box[2 * j] + box[2 * j + 1]) / 2.0;
   }
   for( i = 0; i < ncorners; ++i )
   {
      SCIPdebugMsg(scip, "corner %u: ", i);
      for( j = 0; j < nvars; ++j )
      {
         int varpos = nonfixedpos[j];
         /* if j'th bit of row index i is set, then take upper bound on var j, otherwise lower bound var j
          * we check this by shifting i for j positions to the right and checking whether the last bit is set
          */
         if( (i >> j) & 0x1 )
            corner[varpos] = box[2 * varpos + 1]; /* ub of var */
         else
            corner[varpos] = box[2 * varpos ]; /* lb of var */
         SCIPdebugMsgPrint(scip, "%g, ", corner[varpos]);
         assert(!SCIPisInfinity(scip, REALABS(corner[varpos])));
      }

      funvals[i] = function(corner, nallvars, fundata);

      SCIPdebugMsgPrint(scip, "obj = %e\n", funvals[i]);

      if( funvals[i] == SCIP_INVALID || SCIPisInfinity(scip, REALABS(funvals[i])) )
      {
         SCIPdebugMsg(scip, "cannot compute underestimator; function value at corner is too large %g\n", funvals[i]);
         goto CLEANUP;
      }
   }

   /* clear coefs array; below we only fill in coefs for nonfixed variables */
   BMSclearMemoryArray(facetcoefs, nallvars);

   if( nvars == 1 )
   {
      SCIP_CALL( computeVertexPolyhedralFacetUnivariate(scip, box[2 * nonfixedpos[0]], box[2 * nonfixedpos[0] + 1], funvals[0], funvals[1], success, &facetcoefs[nonfixedpos[0]], facetconstant) );

      /* check whether target has been missed */
      if( *success && overestimate == (*facetconstant + facetcoefs[nonfixedpos[0]] * xstar[nonfixedpos[0]] > targetvalue) )
      {
         SCIPdebugMsg(scip, "computed secant, but missed target %g (facetvalue=%g, overestimate=%u)\n", targetvalue, *facetconstant + facetcoefs[nonfixedpos[0]] * xstar[nonfixedpos[0]], overestimate);
         *success = FALSE;
      }
   }
   else if( nvars == 2 )
   {
      int idx1 = nonfixedpos[0];
      int idx2 = nonfixedpos[1];
      SCIP_Real p1[2] = { box[2*idx1],   box[2*idx2]   }; /* corner 0: 0>>0 & 0x1 = 0, 0>>1 & 0x1 = 0 */
      SCIP_Real p2[2] = { box[2*idx1+1], box[2*idx2]   }; /* corner 1: 1>>0 & 0x1 = 1, 1>>1 & 0x1 = 0 */
      SCIP_Real p3[2] = { box[2*idx1],   box[2*idx2+1] }; /* corner 2: 2>>0 & 0x1 = 0, 2>>1 & 0x1 = 1 */
      SCIP_Real p4[2] = { box[2*idx1+1], box[2*idx2+1] }; /* corner 3: 3>>0 & 0x1 = 1, 3>>1 & 0x1 = 1 */
      SCIP_Real xstar2[2] = { xstar[idx1], xstar[idx2] };
      SCIP_Real coefs[2] = { 0.0, 0.0 };

      SCIP_CALL( computeVertexPolyhedralFacetBivariate(scip, overestimate, p1, p2, p3, p4, funvals[0], funvals[1], funvals[2], funvals[3], xstar2, targetvalue, success, coefs, facetconstant) );

      facetcoefs[idx1] = coefs[0];
      facetcoefs[idx2] = coefs[1];
   }
   else
   {
      SCIP_CALL( computeVertexPolyhedralFacetLP(scip, conshdlr, overestimate, xstar, box, nallvars, nonfixedpos, funvals, nvars, targetvalue, success, facetcoefs, facetconstant) );
   }
   if( !*success )
   {
      SCIPdebugMsg(scip, "no success computing facet, %d vars\n", nvars);
      goto CLEANUP;
   }

   /*
    *  check and adjust facet with the algorithm of Rikun et al.
    */

   maxfaceterror = computeVertexPolyhedralMaxFacetError(scip, overestimate, funvals, box, nallvars, nvars, nonfixedpos, facetcoefs, *facetconstant);

   /* adjust constant part of the facet by maxerror to make it a valid over/underestimator (not facet though) */
   if( maxfaceterror > 0.0 )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Real midval;
      SCIP_Real feastol;

      feastol = SCIPgetStage(scip) == SCIP_STAGE_SOLVING ? SCIPgetLPFeastol(scip) : SCIPfeastol(scip);

      /* evaluate function in middle point to get some idea for a scaling */
      for( j = 0; j < nvars; ++j )
         corner[nonfixedpos[j]] = (box[2 * nonfixedpos[j]] + box[2 * nonfixedpos[j] + 1]) / 2.0;
      midval = function(corner, nallvars, fundata);
      if( midval == SCIP_INVALID )
         midval = 1.0;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* there seem to be numerical problems if the error is too large; in this case we reject the facet */
      if( maxfaceterror > conshdlrdata->vp_adjfacetthreshold * feastol * fabs(midval) )
      {
         SCIPdebugMsg(scip, "ignoring facet due to instability, it cuts off a vertex by %g (midval=%g).\n", maxfaceterror, midval);
         *success = FALSE;
         goto CLEANUP;
      }

      SCIPdebugMsg(scip, "maximum facet error %g (midval=%g), adjust constant to make cut valid!\n", maxfaceterror, midval);

      if( overestimate )
         *facetconstant += maxfaceterror;
      else
         *facetconstant -= maxfaceterror;
   }

   /* if we made it until here, then we have a nice facet */
   assert(*success);

CLEANUP:
   /* free allocated memory */
   SCIPfreeBufferArray(scip, &corner);
   SCIPfreeBufferArray(scip, &funvals);
   SCIPfreeBufferArray(scip, &nonfixedpos);

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** returns the expression of the given nonlinear constraint */
SCIP_EXPR* SCIPgetExprNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->expr;
}

/** gets the left hand side of a nonlinear constraint */
SCIP_Real SCIPgetLhsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets the right hand side of a nonlinear constraint */
SCIP_Real SCIPgetRhsNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets the nonlinear constraint as a nonlinear row representation. */
SCIP_RETCODE SCIPgetNlRowNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow == NULL )
   {
      SCIP_CALL( createNlRow(scip, cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** returns the curvature of the expression of a given nonlinear constraint
 *
 * @note The curvature information is computed during CONSINITSOL.
 */
SCIP_EXPRCURV SCIPgetCurvatureNonlinear(
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->curv;
}

/** checks whether expression of constraint can be represented as quadratic form
 *
 * Only sets `*isquadratic` to TRUE if the whole expression is quadratic (in the non-extended formulation) and non-linear.
 * That is, the expression in each \ref SCIP_QUADEXPR_QUADTERM will be a variable expressions and
 * \ref SCIPgetVarExprVar() can be used to retrieve the variable.
 */
SCIP_RETCODE SCIPcheckQuadraticNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Bool*            isquadratic         /**< buffer to store whether constraint is quadratic */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(isquadratic != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   /* check whether constraint expression is quadratic in extended formulation */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, consdata->expr, isquadratic) );

   /* if not quadratic in non-extended formulation, then do indicate quadratic */
   if( *isquadratic )
      *isquadratic = SCIPexprAreQuadraticExprsVariables(consdata->expr);

   return SCIP_OKAY;
}

/** changes left-hand-side of a nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_RETCODE SCIPchgLhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left-hand-side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("SCIPchgLhsNonlinear can only be called in problem stage.\n");
      return SCIP_INVALIDCALL;
   }

   /* we should have an original constraint */
   assert(SCIPconsIsOriginal(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->lhs == lhs )
      return SCIP_OKAY;

   consdata->lhs = lhs;

   /* not sure we care about any of these flags for original constraints */
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** changes right-hand-side of a nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_RETCODE SCIPchgRhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right-hand-side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("SCIPchgLhsNonlinear can only be called in problem stage.\n");
      return SCIP_INVALIDCALL;
   }

   /* we should have an original constraint */
   assert(SCIPconsIsOriginal(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rhs == rhs )
      return SCIP_OKAY;

   consdata->rhs = rhs;

   /* not sure we care about any of these flags for original constraints */
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** changes expression of a nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_RETCODE SCIPchgExprNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_EXPR*            expr                /**< new expression */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(expr != NULL);

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("SCIPchgExprNonlinear can only be called in problem stage.\n");
      return SCIP_INVALIDCALL;
   }

   /* we should have an original constraint */
   assert(SCIPconsIsOriginal(cons));

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   /* we should not have collected additional data for the expr
    * if some of these asserts fail, we may have to remove it and add some code to keep information up to date
    */
   assert(consdata->nvarexprs == 0);
   assert(consdata->varexprs == NULL);
   assert(!consdata->catchedevents);

   SCIP_CALL( SCIPreleaseExpr(scip, &consdata->expr) );

   /* copy expression, thereby map variables expressions to already existing variables expressions in var2expr map, or augment var2expr map */
   SCIP_CALL( SCIPduplicateExpr(scip, expr, &consdata->expr, mapexprvar, conshdlr, exprownerCreate, (void*)conshdlr) );

   /* not sure we care about any of these flags for original constraints */
   consdata->curv = SCIP_EXPRCURV_UNKNOWN;
   consdata->issimplified = FALSE;
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** adds coef * var to nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_RETCODE SCIPaddLinearVarNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* varexpr;

   assert(scip != NULL);
   assert(cons != NULL);

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("SCIPaddLinearVarNonlinear can only be called in problem stage.\n");
      return SCIP_INVALIDCALL;
   }

   /* we should have an original constraint */
   assert(SCIPconsIsOriginal(cons));

   if( coef == 0.0 )
      return SCIP_OKAY;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   /* we should not have collected additional data for it
    * if some of these asserts fail, we may have to remove it and add some code to keep information up to date
    */
   assert(consdata->nvarexprs == 0);
   assert(consdata->varexprs == NULL);
   assert(!consdata->catchedevents);

   SCIP_CALL( createExprVar(scip, conshdlr, &varexpr, var) );

   /* append to sum, if consdata->expr is sum and not used anywhere else */
   if( SCIPexprGetNUses(consdata->expr) == 1 && SCIPisExprSum(scip, consdata->expr) )
   {
      SCIP_CALL( SCIPappendExprSumExpr(scip, consdata->expr, varexpr, coef) );
   }
   else
   {
      /* create new expression = 1 * consdata->expr + coef * var */
      SCIP_EXPR* children[2] = { consdata->expr, varexpr };
      SCIP_Real coefs[2] = { 1.0, coef };

      SCIP_CALL( SCIPcreateExprSum(scip, &consdata->expr, 2, children, coefs, 0.0, exprownerCreate, (void*)conshdlr) );

      /* release old root expr */
      SCIP_CALL( SCIPreleaseExpr(scip, &children[0]) );
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );

   /* not sure we care about any of these flags for original constraints */
   consdata->issimplified = FALSE;
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** adds coef * expr to nonlinear constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_RETCODE SCIPaddExprNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             coef                /**< coefficient */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR* exprowned;

   assert(scip != NULL);
   assert(cons != NULL);

   if( SCIPgetStage(scip) != SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("SCIPaddLinearVarNonlinear can only be called in problem stage.\n");
      return SCIP_INVALIDCALL;
   }

   /* we should have an original constraint */
   assert(SCIPconsIsOriginal(cons));

   if( coef == 0.0 )
      return SCIP_OKAY;

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->expr != NULL);

   /* we should not have collected additional data for it
    * if some of these asserts fail, we may have to remove it and add some code to keep information up to date
    */
   assert(consdata->nvarexprs == 0);
   assert(consdata->varexprs == NULL);
   assert(!consdata->catchedevents);

   /* copy expression, thereby map variables expressions to already existing variables expressions in var2expr map, or augment var2expr map */
   SCIP_CALL( SCIPduplicateExpr(scip, expr, &exprowned, mapexprvar, conshdlr, exprownerCreate, (void*)conshdlr) );

   /* append to sum, if consdata->expr is sum and not used anywhere else */
   if( SCIPexprGetNUses(consdata->expr) == 1 && SCIPisExprSum(scip, consdata->expr) )
   {
      SCIP_CALL( SCIPappendExprSumExpr(scip, consdata->expr, exprowned, coef) );
   }
   else
   {
      /* create new expression = 1 * consdata->expr + coef * var */
      SCIP_EXPR* children[2] = { consdata->expr, exprowned };
      SCIP_Real coefs[2] = { 1.0, coef };

      SCIP_CALL( SCIPcreateExprSum(scip, &consdata->expr, 2, children, coefs, 0.0, exprownerCreate, (void*)conshdlr) );

      /* release old root expr */
      SCIP_CALL( SCIPreleaseExpr(scip, &children[0]) );
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &exprowned) );

   /* not sure we care about any of these flags for original constraints */
   consdata->issimplified = FALSE;
   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** gets absolute violation of nonlinear constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * If this value is at most SCIPfeastol(), the constraint would be considered feasible.
 */
SCIP_RETCODE SCIPgetAbsViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   )
{
   assert(cons != NULL);
   assert(viol != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, 0L) );
   *viol = getConsAbsViolation(cons);

   return SCIP_OKAY;
}

/** gets scaled violation of nonlinear constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * The scaling that is applied to the absolute violation of the constraint
 * depends on the setting of parameter constraints/nonlinear/violscale.
 */
SCIP_RETCODE SCIPgetRelViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   )
{
   assert(cons != NULL);
   assert(viol != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, 0L) );
   SCIP_CALL( getConsRelViolation(scip, cons, viol, sol, 0L) );

   return SCIP_OKAY;
}

/** returns a variable that appears linearly that may be decreased without making any other constraint infeasible */
void SCIPgetLinvarMayDecreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != NULL);

   /* check for a linear variable that can be increased or decreased without harming feasibility */
   findUnlockedLinearVar(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *var = consdata->linvardecr;
   *coef = consdata->linvardecrcoef;
}

/** returns a variable that appears linearly that may be increased without making any other constraint infeasible */
void SCIPgetLinvarMayIncreaseNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != NULL);

   /* check for a linear variable that can be increased or decreased without harming feasibility */
   findUnlockedLinearVar(scip, cons);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *var = consdata->linvarincr;
   *coef = consdata->linvarincrcoef;
}


/*
 * Methods for Expressions in Nonlinear Constraints
 */

/** returns the number of positive rounding locks of an expression */
int SCIPgetExprNLocksPosNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetOwnerData(expr) != NULL);

   return SCIPexprGetOwnerData(expr)->nlockspos;
}

/** returns the number of negative rounding locks of an expression */
int SCIPgetExprNLocksNegNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetOwnerData(expr) != NULL);

   return SCIPexprGetOwnerData(expr)->nlocksneg;
}

/** returns the variable used for linearizing a given expression (return value might be NULL)
 *
 * @note for variable expression it returns the corresponding variable
 */
SCIP_VAR* SCIPgetExprAuxVarNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   return ownerdata->filterpos >= -1 ? SCIPgetVarExprVar(expr) : ownerdata->auxvar;
}

/** returns the number of enforcements for an expression */
int SCIPgetExprNEnfosNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetOwnerData(expr) != NULL);

   return SCIPexprGetOwnerData(expr)->nenfos;
}

/** returns the data for one of the enforcements of an expression */
void SCIPgetExprEnfoDataNonlinear(
   SCIP_EXPR*            expr,               /**< expression */
   int                   idx,                /**< position of enforcement in enfos array */
   SCIP_NLHDLR**         nlhdlr,             /**< buffer to store nlhldr */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< buffer to store nlhdlr data for expression, or NULL */
   SCIP_NLHDLR_METHOD* nlhdlrparticipation, /**< buffer to store methods where nonlinear handler participates, or NULL */
   SCIP_Bool*            sepabelowusesactivity, /**< buffer to store whether sepabelow uses activity of some expression, or NULL */
   SCIP_Bool*            sepaaboveusesactivity, /**< buffer to store whether sepaabove uses activity of some expression, or NULL */
   SCIP_Real*            auxvalue            /**< buffer to store current auxvalue, or NULL */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   assert(idx >= 0);
   assert(idx < ownerdata->nenfos);
   assert(ownerdata->enfos[idx] != NULL);
   assert(nlhdlr != NULL);

   *nlhdlr = ownerdata->enfos[idx]->nlhdlr;

   if( nlhdlrexprdata != NULL )
      *nlhdlrexprdata = ownerdata->enfos[idx]->nlhdlrexprdata;

   if( nlhdlrparticipation != NULL )
      *nlhdlrparticipation = ownerdata->enfos[idx]->nlhdlrparticipation;

   if( sepabelowusesactivity != NULL )
      *sepabelowusesactivity = ownerdata->enfos[idx]->sepabelowusesactivity;

   if( sepaaboveusesactivity != NULL )
      *sepaaboveusesactivity = ownerdata->enfos[idx]->sepaaboveusesactivity;

   if( auxvalue != NULL )
      *auxvalue = ownerdata->enfos[idx]->auxvalue;
}

/** sets the auxiliary value of expression for one of the enforcements of an expression */
void SCIPsetExprEnfoAuxValueNonlinear(
   SCIP_EXPR*            expr,               /**< expression */
   int                   idx,                /**< position of enforcement in enfos array */
   SCIP_Real             auxvalue            /**< the new value of auxval */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   assert(idx >= 0);
   assert(idx < ownerdata->nenfos);
   assert(ownerdata->enfos[idx] != NULL);

   ownerdata->enfos[idx]->auxvalue = auxvalue;
}

/** number of nonlinear handlers whose activity computation and propagation methods depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
unsigned int SCIPgetExprNPropUsesActivityNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetOwnerData(expr) != NULL);

   return SCIPexprGetOwnerData(expr)->nactivityusesprop;
}

/** number of nonlinear handlers whose separation methods (estimate or enforcement) depend on the activity of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
unsigned int SCIPgetExprNSepaUsesActivityNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetOwnerData(expr) != NULL);

   return SCIPexprGetOwnerData(expr)->nactivityusessepa;
}

/** number of nonlinear handlers whose separation methods (estimate or enforcement) use auxiliary variable of the expression
 *
 * @note This method can only be used after the detection methods of the nonlinear handlers have been called.
 */
unsigned int SCIPgetExprNAuxvarUsesNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(expr != NULL);
   assert(SCIPexprGetOwnerData(expr) != NULL);

   return SCIPexprGetOwnerData(expr)->nauxvaruses;
}

/** method to be called by a nlhdlr during NLHDLRDETECT to notify an expression that it will be used
 *
 * - if `useauxvar` is enabled, then ensures that an auxiliary variable will be created in INITLP
 * - if `useactivityforprop` or `useactivityforsepa{below,above}` is enabled, then ensured that activity will be updated for `expr`
 * - if `useactivityforprop` is enabled, then increments the count returned by SCIPgetExprNPropUsesActivityNonlinear()
 * - if `useactivityforsepa{below,above}` is enabled, then increments the count returned by SCIPgetExprNSepaUsesActivityNonlinear()
 *   and also increments this count for all variables in the expression.
 *
 * The distinction into `useactivityforprop` and `useactivityforsepa{below,above}` is to recognize variables which domain influences
 * under/overestimators. Domain propagation routines (like OBBT) may invest more work for these variables.
 * The distinction into `useactivityforsepabelow` and `useactivityforsepaabove` is to recognize whether a nlhdlr that called this method
 * will use activity of `expr` in enfomethod \ref SCIP_NLHDLR_METHOD_SEPABELOW or \ref SCIP_NLHDLR_METHOD_SEPAABOVE.
 */
SCIP_RETCODE SCIPregisterExprUsageNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool             useauxvar,          /**< whether an auxiliary variable will be used for estimate or cut generation */
   SCIP_Bool             useactivityforprop, /**< whether activity of expr will be used by domain propagation or activity calculation (inteval) */
   SCIP_Bool             useactivityforsepabelow, /**< whether activity of expr will be used by underestimation */
   SCIP_Bool             useactivityforsepaabove  /**< whether activity of expr will be used by overestimation */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   /* do not store auxvar request for variable expressions */
   if( useauxvar && SCIPisExprVar(scip, expr) )
      useauxvar = FALSE;

   if( ownerdata->nenfos >= 0 &&
      ( (ownerdata->nactivityusesprop == 0 && ownerdata->nactivityusessepa == 0 && (useactivityforprop || useactivityforsepabelow || useactivityforsepaabove)) ||
        (ownerdata->nauxvaruses == 0 && useauxvar)
      ) )
   {
      /* if we already have ran detect of nlhdlrs on expr (nenfos >= 0), then we need to rerun detection if
       * we require additional enforcement methods, that is,
       * - activity of expr was not used before but will be used now, or
       * - auxiliary variable of expr was not required before but will be used now
       */
      SCIP_CALL( freeEnfoData(scip, expr, FALSE) );
   }

   if( useauxvar )
      ++ownerdata->nauxvaruses;

   if( useactivityforprop )
      ++ownerdata->nactivityusesprop;

   if( useactivityforsepabelow || useactivityforsepaabove )
      ++ownerdata->nactivityusessepa;

   /* remember that SCIPregisterExprUsageNonlinear() has been called with useactivityforsepa{below,above}=TRUE; this
    * information is used in detectNlhdlr()
    */
   if( useactivityforsepabelow )
      SCIPconshdlrGetData(ownerdata->conshdlr)->registerusesactivitysepabelow = TRUE;
   if( useactivityforsepaabove )
      SCIPconshdlrGetData(ownerdata->conshdlr)->registerusesactivitysepaabove = TRUE;

   if( useactivityforprop )
   {
      /* if activity will be used for propagation, then make sure there is a valid activity
       * this way, we can do a reversepropcall after detectNlhdlr
       */
      SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   }

   /* increase the nactivityusedsepa counter for all variables used in the given expression */
   if(( useactivityforsepabelow || useactivityforsepaabove) && SCIPexprGetNChildren(expr) > 0 )
   {
      SCIP_EXPRITER* it;

      /* create and initialize iterator */
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );
      SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

      for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
         if( SCIPisExprVar(scip, expr) )
            ++SCIPexprGetOwnerData(expr)->nactivityusessepa;

      /* free iterator */
      SCIPfreeExpriter(&it);
   }

   return SCIP_OKAY;
}

/** computes absolute violation for auxvar relation in an expression w.r.t. original variables
 *
 * Assume the expression is f(x), where x are original (i.e., not auxiliary) variables.
 * Assume that f(x) is associated with auxiliary variable z.
 *
 * If there are negative locks, then returns the violation of z &le; f(x) and sets `violover` to TRUE.
 * If there are positive locks, then returns the violation of z &ge; f(x) and sets `violunder` to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z = f(x).
 *
 * If necessary, f is evaluated in the given solution. If that fails (domain error),
 * then `viol` is set to SCIPinfinity() and both `violover` and `violunder` are set to TRUE.
 */
SCIP_RETCODE SCIPgetExprAbsOrigViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Longint          soltag,             /**< tag of solution */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(x) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(x) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* make sure expression has been evaluated */
   SCIP_CALL( SCIPevalExpr(scip, expr, sol, soltag) );

   /* get violation from internal method */
   *viol = getExprAbsOrigViolation(scip, expr, sol, violunder, violover);

   return SCIP_OKAY;
}

/** computes absolute violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * If there are negative locks, then returns the violation of z &le; f(w) and sets `violover` to TRUE.
 * If there are positive locks, then returns the violation of z &ge; f(w) and sets `violunder` to TRUE.
 * Of course, if there both negative and positive locks, then return the violation of z = f(w).
 *
 * If the given value of f(w) is SCIP_INVALID, then `viol` is set to SCIPinfinity() and
 * both `violover` and `violunder` are set to TRUE.
 */
SCIP_RETCODE SCIPgetExprAbsAuxViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* get violation from internal method */
   *viol = getExprAbsAuxViolation(scip, expr, auxvalue, sol, violunder, violover);

   return SCIP_OKAY;
}


/** computes relative violation for auxvar relation in an expression w.r.t. auxiliary variables
 *
 * Assume the expression is f(w), where w are auxiliary variables that were introduced by some nlhdlr.
 * Assume that f(w) is associated with auxiliary variable z.
 *
 * Taking the absolute violation from SCIPgetExprAbsAuxViolationNonlinear(), this function returns
 * the absolute violation divided by max(1,|f(w)|).
 *
 * If the given value of f(w) is SCIP_INVALID, then `viol` is set to SCIPinfinity() and
 * both `violover` and `violunder` are set to TRUE.
 */
SCIP_RETCODE SCIPgetExprRelAuxViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             auxvalue,           /**< the value of f(w) */
   SCIP_SOL*             sol,                /**< solution that has been evaluated */
   SCIP_Real*            viol,               /**< buffer to store computed violation */
   SCIP_Bool*            violunder,          /**< buffer to store whether z >= f(w) is violated, or NULL */
   SCIP_Bool*            violover            /**< buffer to store whether z <= f(w) is violated, or NULL */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(viol != NULL);

   /* get violation from internal method */
   *viol = getExprAbsAuxViolation(scip, expr, auxvalue, sol, violunder, violover);

   if( !SCIPisInfinity(scip, *viol) )
   {
      assert(auxvalue != SCIP_INVALID);
      /* TODO maybe we should rather use max(eps,|auxvalue|)? */
      *viol /= MAX(1.0, REALABS(auxvalue));
   }

   return SCIP_OKAY;
}

/** returns bounds on the expression
 *
 * This gives an intersection of bounds from
 * - activity calculation (SCIPexprGetActivity()), if valid,
 * - auxiliary variable, if present,
 * - stored by SCIPtightenExprIntervalNonlinear() during domain propagation
 *
 * @note The returned interval can be empty!
 */
SCIP_INTERVAL SCIPgetExprBoundsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL bounds;

   assert(scip != NULL);
   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);

   /* SCIPdebugMsg(scip, "get bounds expr %p:", expr); */

   /* start with propbounds if they belong to current propagation */
   if( ownerdata->propboundstag == conshdlrdata->curpropboundstag )
   {
      bounds = ownerdata->propbounds;
      /* SCIPdebugMsgPrint(scip, " propbounds [%.15g,%.15g]", ownerdata->propbounds.inf, ownerdata->propbounds.sup); */
   }
   else
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &bounds);

   if( SCIPexprGetActivityTag(expr) >= conshdlrdata->lastboundrelax )
   {
      /* apply propbounds to expr activity, but ensure it's not-empty if very close disjoint intervals */
      /* SCIPdebugMsgPrint(scip, " activity [%.15g,%.15g]", expr->activity.inf, expr->activity.sup); */
      SCIPintervalIntersectEps(&bounds, SCIPepsilon(scip), SCIPexprGetActivity(expr), bounds);
   }

   if( ownerdata->auxvar != NULL )
   {
      /* apply auxiliary variable bounds to bounds */
      SCIP_INTERVAL auxvarbounds;

      auxvarbounds = conshdlrdata->intevalvar(scip, ownerdata->auxvar, conshdlrdata);
      /* SCIPdebugMsgPrint(scip, " auxvar [%.15g,%.15g]", auxvarbounds.inf, auxvarbounds.sup); */
      SCIPintervalIntersectEps(&bounds, SCIPepsilon(scip), bounds, auxvarbounds);
   }

   /* SCIPdebugMsgPrint(scip, " -> [%.15g,%.15g]\n", bounds.inf, bounds.sup); */

   return bounds;
}

/** informs the expression about new bounds that can be used for reverse-propagation and to tighten bounds of
 * corresponding (auxiliary) variable (if any)
 *
 * @attention this function should only be called during domain propagation in cons_nonlinear
 */
SCIP_RETCODE SCIPtightenExprIntervalNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be tightened */
   SCIP_INTERVAL         newbounds,          /**< new bounds for the expression */
   SCIP_Bool*            cutoff,             /**< buffer to store whether a cutoff was detected */
   int*                  ntightenings        /**< buffer to add the total number of tightenings, or NULL */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);
   assert(ownerdata->conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);

   /* the code below assumes that current activity is valid
    * if it turns out that we cannot ensure that, then we should change code
    */
   assert(SCIPexprGetActivityTag(expr) >= conshdlrdata->lastboundrelax || SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(expr)));
   assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPexprGetActivity(expr)));

   *cutoff = FALSE;

#ifdef DEBUG_PROP
   SCIPdebugMsg(scip, "Trying to tighten bounds of expr ");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPdebugMsgPrint(scip, " with activity [%.15g,%.15g] to [%.15g,%.15g] (force=%d)\n", SCIPexprGetActivity(expr).inf, SCIPexprGetActivity(expr).sup, newbounds.inf, newbounds.sup, conshdlrdata->forceboundtightening);
#endif

   if( SCIPexprIsIntegral(expr) )
   {
      /* apply integrality to new bounds
       * it should be ok to use normal ceil() and floor(), but for safety, we use SCIPceil and SCIPfloor for now
       */
      if( newbounds.inf > -SCIP_INTERVAL_INFINITY )
         newbounds.inf = SCIPceil(scip, newbounds.inf);
      if( newbounds.sup <  SCIP_INTERVAL_INFINITY )
         newbounds.sup = SCIPfloor(scip, newbounds.sup);
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, " applied integrality: [%.15g,%.15g]\n", newbounds.inf, newbounds.sup);
#endif
   }

   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      SCIPdebugMsg(scip, " cut off due to new bounds being empty\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* treat the new bounds as empty if either the lower/upper bound is above/below +/- SCIPinfinity() */
   if( SCIPisInfinity(scip, newbounds.inf) || SCIPisInfinity(scip, -newbounds.sup) )
   {
      SCIPdebugMsg(scip, " cut off due to new bounds being beyond infinity\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* tighten newbounds w.r.t. existing expr->propbounds or activity */
   if( ownerdata->propboundstag == conshdlrdata->curpropboundstag )
   {
      /* if already having propbounds in expr, then tighten newbounds by propbounds */
      SCIPintervalIntersectEps(&newbounds, SCIPepsilon(scip), ownerdata->propbounds, newbounds);
   }
   else
   {
      /* first time we have propbounds for expr in this propagation rounds:
       * intersect with activity (though don't let it become empty if very close intervals)
       */
      SCIPintervalIntersectEps(&newbounds, SCIPepsilon(scip), SCIPexprGetActivity(expr), newbounds);
   }
#ifdef DEBUG_PROP
   SCIPdebugMsg(scip, " applied %s: [%.20g,%.20g]\n", ownerdata->propboundstag == conshdlrdata->curpropboundstag ? "previous propbounds" : "activity", newbounds.inf, newbounds.sup);
#endif

   /* check if the new bounds lead to an empty interval */
   if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, newbounds) )
   {
      SCIPdebugMsg(scip, " cut off due to empty intersection with previous propbounds or activity\n");

      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* if expr is not constant or variable, then store newbounds in expr->propbounds
    * - for constant, the intersection with activity should have been sufficient to determine infeasibilty
    * - for variable, the tightenAuxVarBounds call below should be suffient to have to new bounds acknowledged
    */
   if( SCIPexprGetNChildren(expr) > 0 )
   {
      ownerdata->propbounds = newbounds;
      ownerdata->propboundstag = conshdlrdata->curpropboundstag;
   }

   /* if updated propbounds do not allow a sufficient tightening, then do not consider adding to queue for reverse
    * propagation or update of auxvar bounds
    * TODO? if we first had a considerable tightening and then only get small tightenings under the same
    *   curpropboundstag, then these will still be considered as isIntervalBetter, since we compare with activity here and
    *   not with the propbounds as set in the beginning; I'm not sure, though, that comparing always with previous
    *   propbounds would be better, since a number of small updates to propbounds could eventually lead to a considerable
    *   one or should we not even update propbounds to newbounds if the update is small?
    */
   if( !isIntervalBetter(scip, conshdlrdata->forceboundtightening, newbounds, SCIPexprGetActivity(expr)) )
   {
#ifdef DEBUG_PROP
      SCIPdebugMsg(scip, " new bounds [%g,%g] for expr %p not sufficiently tighter than activity -- not adding to propqueue or tightening auxvar\n", newbounds.inf, newbounds.sup, (void*)expr);
#endif
      return SCIP_OKAY;
   }

   if( SCIPexprGetNChildren(expr) > 0 && !ownerdata->inpropqueue && (ownerdata->nactivityusesprop > 0 || ownerdata->nactivityusessepa > 0 || ownerdata->nenfos < 0) )
   {
      /* add expression to propagation queue if not there yet and not var or constant and
       * if it should have a nlhdlr with a reverseprop callback or nlhdlrs are not initialized yet (nenfos < 0)
       */
#ifdef DEBUG_PROP
         SCIPdebugMsg(scip, " insert expr <%p> (%s) into reversepropqueue\n", (void*)expr, SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)));
#endif
         SCIP_CALL( SCIPqueueInsert(conshdlrdata->reversepropqueue, expr) );
         ownerdata->inpropqueue = TRUE;
   }

   /* update bounds on variable or auxiliary variable */
   SCIP_CALL( tightenAuxVarBounds(scip, ownerdata->conshdlr, expr, newbounds, cutoff, ntightenings) );

   return SCIP_OKAY;
}

/** mark constraints that include this expression to be propagated again
 *
 * This can be used by, e.g., nlhdlrs, to trigger a new propagation of constraints without
 * a change of variable bounds, e.g., because new information on the expression is available
 * that could potentially lead to tighter expression activity values.
 *
 * Note, that this call marks also constraints for propagation which only share some variable
 * with this expression.
 */
SCIP_RETCODE SCIPmarkExprPropagateNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression to propagate again */
   )
{
   SCIP_EXPRITER* it;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR_OWNERDATA* ownerdata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   SCIPincrementCurBoundsTagNonlinear(ownerdata->conshdlr, FALSE);

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

   for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      if( !SCIPisExprVar(scip, expr) )
         continue;

      ownerdata = SCIPexprGetOwnerData(expr);
      assert(ownerdata != NULL);

      for( c = 0; c < ownerdata->nconss; ++c )
      {
         consdata = SCIPconsGetData(ownerdata->conss[c]);
         assert(consdata != NULL);
         consdata->ispropagated = FALSE;
      }
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}

/** adds violation-branching score to an expression
 *
 * Adds a score to the expression-specific violation-branching score, thereby marking it as branching candidate.
 * The expression must either be a variable expression or have an aux-variable.
 * In the latter case, branching on auxiliary variables must have been enabled.
 * In case of doubt, use SCIPaddExprsViolScoreNonlinear(). Roughly, the difference between these functions is that the current
 * function adds `violscore` to the expression directly, while SCIPaddExprsViolScoreNonlinear() will split the
 * violation score among all the given expressions according to parameter constraints/nonlinear/branching/violsplit.
 *
 * @see SCIPaddExprsViolScoreNonlinear()
 */
void SCIPaddExprViolScoreNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression where to add branching score */
   SCIP_Real             violscore           /**< violation score to add to expression */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(violscore >= 0.0);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);

   /* if not allowing to branch on auxvars, then expr must be a var-expr */
   assert(branchAuxNonlinear(scip, ownerdata->conshdlr) || SCIPisExprVar(scip, expr));
   /* if allowing to branch on auxvars, then expr must be a var-expr or have an auxvar */
   assert(!branchAuxNonlinear(scip, ownerdata->conshdlr) || SCIPisExprVar(scip, expr) || ownerdata->auxvar != NULL);

   /* reset branching score if we are in a different enfo round */
   if( ownerdata->violscoretag != conshdlrdata->enforound )
   {
      ownerdata->violscoresum = violscore;
      ownerdata->violscoremax = violscore;
      ownerdata->nviolscores = 1;
      ownerdata->violscoretag = conshdlrdata->enforound;
      return;
   }

   ownerdata->violscoresum += violscore;
   if( violscore > ownerdata->violscoremax )
      ownerdata->violscoremax = violscore;
   ++ownerdata->nviolscores;
}

/** adds violation-branching score to a set of expressions, distributing the score among all the expressions
 *
 * Each expression must either be a variable expression or have an aux-variable.
 * If branching on aux-variables is disabled, then the violation branching score will be distributed among all
 * variables present in `exprs`.
 */
SCIP_RETCODE SCIPaddExprsViolScoreNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           exprs,              /**< expressions where to add branching score */
   int                   nexprs,             /**< number of expressions */
   SCIP_Real             violscore,          /**< violation score to add to expression */
   SCIP_SOL*             sol,                /**< current solution */
   SCIP_Bool*            success             /**< buffer to store whether at least one violscore was added */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR** varexprs;
   SCIP_EXPR* e;
   int nvars;
   int varssize;
   int i;

   assert(exprs != NULL || nexprs == 0);
   assert(success != NULL);

   if( nexprs == 0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* if allowing to branch on auxiliary variables, then call internal addConsExprExprsViolScore immediately */
   if( branchAuxNonlinear(scip, SCIPexprGetOwnerData(exprs[0])->conshdlr) )
   {
      addExprsViolScore(scip, exprs, nexprs, violscore, sol, success);
      return SCIP_OKAY;
   }

   /* if not allowing to branch on aux vars, then create new array containing var expressions that exprs depend on */
   nvars = 0;
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, varssize) );

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );
   SCIP_CALL( SCIPexpriterInit(it, NULL, SCIP_EXPRITER_DFS, FALSE) );

   for( i = 0; i < nexprs; ++i )
   {
      for( e = SCIPexpriterRestartDFS(it, exprs[i]); !SCIPexpriterIsEnd(it); e = SCIPexpriterGetNext(it) )
      {
         assert(e != NULL);

         if( SCIPisExprVar(scip, e) )
         {
            /* add variable expression to vars array */
            if( varssize == nvars )
            {
               varssize = SCIPcalcMemGrowSize(scip, nvars + 1);
               SCIP_CALL( SCIPreallocBufferArray(scip, &varexprs, varssize) );
            }
            assert(varssize > nvars);

            varexprs[nvars++] = e;
         }
      }
   }

   SCIPfreeExpriter(&it);

   addExprsViolScore(scip, varexprs, nvars, violscore, sol, success);

   SCIPfreeBufferArray(scip, &varexprs);

   return SCIP_OKAY;
}

/** gives violation-branching score stored in expression, or 0.0 if no valid score has been stored */
SCIP_Real SCIPgetExprViolScoreNonlinear(
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->enforound != ownerdata->violscoretag )
      return 0.0;

   if( ownerdata->nviolscores == 0 )
      return 0.0;

   switch( conshdlrdata->branchscoreagg )
   {
      case 'a' :
         /* average */
         return ownerdata->violscoresum / ownerdata->nviolscores;

      case 'm' :
         /* maximum */
         return ownerdata->violscoremax;

      case 's' :
         /* sum */
         return ownerdata->violscoresum;

      default:
         SCIPerrorMessage("Invalid value %c for branchscoreagg parameter\n", conshdlrdata->branchscoreagg);
         SCIPABORT();
         return SCIP_INVALID;
   }
}

/** returns the partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @see SCIPexprGetDerivative()
 */
SCIP_Real SCIPgetExprPartialDiffNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< root expression of constraint used in the last SCIPevalExprGradient() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* varexpr;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(var != NULL);

   /* return 0.0 for value expression */
   if( SCIPisExprValue(scip, expr) )
   {
      assert(SCIPexprGetDerivative(expr) == 0.0);
      return 0.0;
   }

   /* check if an error occurred during the last SCIPevalExprGradient() call */
   if( SCIPexprGetDerivative(expr) == SCIP_INVALID )
      return SCIP_INVALID;

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);

   /* use variable to expressions mapping which is stored in the constraint handler data */
   assert(SCIPhashmapExists(conshdlrdata->var2expr, var));

   varexpr = (SCIP_EXPR*)SCIPhashmapGetImage(conshdlrdata->var2expr, var);
   assert(varexpr != NULL);
   assert(SCIPisExprVar(scip, varexpr));

   /* use difftag to decide whether the variable belongs to the expression */
   return (SCIPexprGetDiffTag(expr) != SCIPexprGetDiffTag(varexpr)) ? 0.0 : SCIPexprGetDerivative(varexpr);
}

/** returns the var's coordinate of Hu partial derivative of an expression w.r.t. a variable (or SCIP_INVALID if there was an evaluation error)
 *
 * @see SCIPexprGetBardot()
 */
SCIP_Real SCIPgetExprPartialDiffGradientDirNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< root expression of constraint used in the last SCIPevalExprHessianDir() call */
   SCIP_VAR*             var                 /**< variable (needs to be in the expression) */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EXPR* varexpr;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(var != NULL);

   /* return 0.0 for value expression */
   if( SCIPisExprValue(scip, expr) )
      return 0.0;

   /* check if an error occurred during the last SCIPevalExprHessianDir() call */
   if( SCIPexprGetBardot(expr) == SCIP_INVALID )
      return SCIP_INVALID;

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(ownerdata->conshdlr);
   assert(conshdlrdata != NULL);

   /* use variable to expressions mapping which is stored in the constraint handler data;
    * if this fails it means that we are asking for the var's component of H*u for a var
    * that doesn't appear in any nonlinear constraint, so maybe we can also just return 0.0
    */
   assert(SCIPhashmapExists(conshdlrdata->var2expr, var));

   varexpr = (SCIP_EXPR*)SCIPhashmapGetImage(conshdlrdata->var2expr, var);
   assert(varexpr != NULL);
   assert(SCIPisExprVar(scip, varexpr));

   /* use difftag to decide whether the variable belongs to the expression */
   return (SCIPexprGetDiffTag(expr) != SCIPexprGetDiffTag(varexpr)) ? 0.0 : SCIPexprGetBardot(varexpr);
}

/** evaluates quadratic term in a solution w.r.t. auxiliary variables
 *
 * \note This requires that for every expr used in the quadratic data, a variable or auxiliary variable is available.
 */
SCIP_Real SCIPevalExprQuadraticAuxNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   )
{
   SCIP_Real auxvalue;
   int nlinexprs;
   SCIP_Real* lincoefs;
   SCIP_EXPR** linexprs;
   int nquadexprs;
   int nbilinexprs;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);

   SCIPexprGetQuadraticData(expr, &auxvalue, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, &nbilinexprs, NULL, NULL);

   /* linear terms */
   for( i = 0; i < nlinexprs; ++i )
   {
      assert(SCIPgetExprAuxVarNonlinear(linexprs[i]) != NULL);
      auxvalue += lincoefs[i] * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(linexprs[i]));
   }

   /* quadratic terms */
   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_EXPR* quadexprterm;
      SCIP_Real lincoef;
      SCIP_Real sqrcoef;
      SCIP_Real solval;

      SCIPexprGetQuadraticQuadTerm(expr, i, &quadexprterm, &lincoef, &sqrcoef, NULL, NULL, NULL);

      assert(SCIPgetExprAuxVarNonlinear(quadexprterm) != NULL);

      solval = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(quadexprterm));
      auxvalue += (lincoef + sqrcoef * solval) * solval;
   }

   /* bilinear terms */
   for( i = 0; i < nbilinexprs; ++i )
   {
      SCIP_EXPR* expr1;
      SCIP_EXPR* expr2;
      SCIP_Real coef;

      SCIPexprGetQuadraticBilinTerm(expr, i, &expr1, &expr2, &coef, NULL, NULL);

      assert(SCIPgetExprAuxVarNonlinear(expr1) != NULL);
      assert(SCIPgetExprAuxVarNonlinear(expr2) != NULL);
      auxvalue += coef * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr1)) * SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(expr2));
   }

   return auxvalue;
}

/**@addtogroup PublicNlhdlrInterfaceMethods
 * @{
 */

/** creates a nonlinear handler and includes it into the nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlr,             /**< buffer where to store nonlinear handler */
   const char*           name,               /**< name of nonlinear handler (must not be NULL) */
   const char*           desc,               /**< description of nonlinear handler (can be NULL) */
   int                   detectpriority,     /**< detection priority of nonlinear handler */
   int                   enfopriority,       /**< enforcement priority of nonlinear handler */
   SCIP_DECL_NLHDLRDETECT((*detect)),        /**< structure detection callback of nonlinear handler */
   SCIP_DECL_NLHDLREVALAUX((*evalaux)),      /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_NLHDLRDATA*      nlhdlrdata          /**< data of nonlinear handler (can be NULL) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(detect != NULL);

   /* find myself */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create nlhdlr */
   SCIP_CALL( SCIPnlhdlrCreate(scip, nlhdlr, name, desc, detectpriority, enfopriority, detect, evalaux, nlhdlrdata) );

   /* include into constraint handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &conshdlrdata->nlhdlrs, &conshdlrdata->nlhdlrssize, conshdlrdata->nnlhdlrs+1) );

   conshdlrdata->nlhdlrs[conshdlrdata->nnlhdlrs] = *nlhdlr;
   ++conshdlrdata->nnlhdlrs;

   /* sort nonlinear handlers by detection priority, in decreasing order
    * will happen in INIT, so only do when called late
    */
   if( SCIPgetStage(scip) > SCIP_STAGE_INIT && conshdlrdata->nnlhdlrs > 1 )
      SCIPsortDownPtr((void**)conshdlrdata->nlhdlrs, SCIPnlhdlrComp, conshdlrdata->nnlhdlrs);

   return SCIP_OKAY;
}

/** get number of nonlinear handler */
int SCIPgetNNlhdlrsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nnlhdlrs;
}

/** get nonlinear handlers */
SCIP_NLHDLR** SCIPgetNlhdlrsNonlinear(
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nlhdlrs;
}

/** returns a nonlinear handler of a given name (or NULL if not found) */
SCIP_NLHDLR* SCIPfindNlhdlrNonlinear(
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   const char*           name                /**< name of nonlinear handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int h;

   assert(conshdlr != NULL);
   assert(name != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
      if( strcmp(SCIPnlhdlrGetName(conshdlrdata->nlhdlrs[h]), name) == 0 )
         return conshdlrdata->nlhdlrs[h];

   return NULL;
}

/** gives expression data that a given nonlinear handler stored in an expression
 *
 * Returns NULL if expr has not been detected by nlhdlr or nlhdlr did not store data.
 */
SCIP_NLHDLREXPRDATA* SCIPgetNlhdlrExprDataNonlinear(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPR_OWNERDATA* ownerdata;
   int e;

   assert(nlhdlr != NULL);
   assert(expr != NULL);

   ownerdata = SCIPexprGetOwnerData(expr);
   assert(ownerdata != NULL);

   for( e = 0; e < ownerdata->nenfos; ++e )
      if( ownerdata->enfos[e]->nlhdlr == nlhdlr )
         return ownerdata->enfos[e]->nlhdlrexprdata;

   return NULL;
}

/** @} */
