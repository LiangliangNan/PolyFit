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

/**@file   cons_quadratic.c
 * @brief  constraint handler for quadratic constraints \f$\textrm{lhs} \leq \sum_{i,j=1}^n a_{i,j} x_i x_j + \sum_{i=1}^n b_i x_i \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * 
 * @todo SCIP might fix linear variables on +/- infinity; remove them in presolve and take care later
 * @todo round constraint sides to integers if all coefficients and variables are (impl.) integer
 * @todo constraints in one variable should be replaced by linear or bounddisjunction constraint
 * @todo check if some quadratic terms appear in several constraints and try to simplify (e.g., nous1)
 * @todo skip separation in enfolp if for current LP (check LP id) was already separated
 * @todo watch unbounded variables to enable/disable propagation
 * @todo sort order in bilinvar1/bilinvar2 such that the var which is involved in more terms is in bilinvar1, and use this info propagate and AddLinearReform
 * @todo underestimate for multivariate concave quadratic terms as in cons_nonlinear
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h> /* for strcmp */ 
#include <ctype.h>  /* for isspace */
#include <math.h>

#define SCIP_PRIVATE_ROWPREP

#include "scip/cons_nonlinear.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_and.h"
#include "scip/cons_varbound.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/intervalarith.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"
#include "nlpi/nlpi.h"
#include "nlpi/nlpi_ipopt.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "quadratic"
#define CONSHDLR_DESC          "quadratic constraints of the form lhs <= b' x + x' A x <= rhs"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -50 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING            SCIP_PRESOLTIMING_ALWAYS /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */
#define NONLINCONSUPGD_PRIORITY   40000 /**< priority of upgrading nonlinear constraints */
#define INITLPMAXVARVAL          1000.0 /**< maximal absolute value of variable for still generating a linearization cut at that point in initlp */

/* Activating this define enables reformulation of bilinear terms x*y with implications from x to y into linear terms.
 * However, implications are not enforced by SCIP. Thus, if, e.g., the used implication was derived from this constraint and we then reformulate the constraint,
 * then the implication may not be enforced in a solution.
 * This issue need to be fixed before this feature can be enabled.
 */
/* #define CHECKIMPLINBILINEAR */

/* enable new propagation for bivariate quadratic terms */
#define PROPBILINNEW

/* epsilon for differentiating between a boundary and interior point */
#define INTERIOR_EPS 1e-1

/* scaling factor for gauge function */
#define GAUGESCALE 0.99999

#define ROWPREP_SCALEUP_VIOLNONZERO    (minviol / 10.0)            /**< minimal violation for considering up-scaling of rowprep (we want to avoid upscaling very small violations) */
#define ROWPREP_SCALEUP_MINVIOLFACTOR  2.0                         /**< scale up will target a violation of ~MINVIOLFACTOR*minviol, where minviol is given by caller */
#define ROWPREP_SCALEUP_MAXMINCOEF     (1.0 / SCIPfeastol(scip))   /**< scale up only if min. coef is below this number (before scaling) */
#define ROWPREP_SCALEUP_MAXMAXCOEF     SCIPgetHugeValue(scip)      /**< scale up only if max. coef will not exceed this number by scaling */
#define ROWPREP_SCALEUP_MAXSIDE        SCIPgetHugeValue(scip)      /**< scale up only if side will not exceed this number by scaling */
#define ROWPREP_SCALEDOWN_MINMAXCOEF   (1.0 / SCIPfeastol(scip))   /**< scale down if max. coef is at least this number (before scaling) */
#define ROWPREP_SCALEDOWN_MINCOEF      SCIPfeastol(scip)           /**< scale down only if min. coef does not drop below this number by scaling */

/*
 * Data structures
 */

/** eventdata for variable bound change events in quadratic constraints */
struct SCIP_QuadVarEventData
{
   SCIP_CONS*            cons;               /**< the constraint */
   int                   varidx;             /**< the index of the variable which bound change is caught, positive for linear variables, negative for quadratic variables */
   int                   filterpos;          /**< position of eventdata in SCIP's event filter */
};

/** Data of a quadratic constraint. */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */

   int                   nlinvars;           /**< number of linear variables */
   int                   linvarssize;        /**< length of linear variable arrays */
   SCIP_VAR**            linvars;            /**< linear variables */
   SCIP_Real*            lincoefs;           /**< coefficients of linear variables */
   SCIP_QUADVAREVENTDATA** lineventdata;       /**< eventdata for bound change of linear variable */

   int                   nquadvars;          /**< number of variables in quadratic terms */
   int                   quadvarssize;       /**< length of quadratic variable terms arrays */
   SCIP_QUADVARTERM*     quadvarterms;       /**< array with quadratic variable terms */

   int                   nbilinterms;        /**< number of bilinear terms */
   int                   bilintermssize;     /**< length of bilinear term arrays */
   SCIP_BILINTERM*       bilinterms;         /**< bilinear terms array */
   int*                  bilintermsidx;      /**< unique index of each bilinear term xy in the bilinestimators array of the constraint handler data */

   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */

   unsigned int          linvarssorted:1;    /**< are the linear variables already sorted? */
   unsigned int          linvarsmerged:1;    /**< are equal linear variables already merged? */
   unsigned int          quadvarssorted:1;   /**< are the quadratic variables already sorted? */
   unsigned int          quadvarsmerged:1;   /**< are equal quadratic variables already merged? */
   unsigned int          bilinsorted:1;      /**< are the bilinear terms already sorted? */
   unsigned int          bilinmerged:1;      /**< are equal bilinear terms (and bilinear terms with zero coefficient) already merged? */

   unsigned int          isconvex:1;         /**< is quadratic function is convex ? */
   unsigned int          isconcave:1;        /**< is quadratic function is concave ? */
   unsigned int          iscurvchecked:1;    /**< is quadratic function checked on convexity or concavity ? */
   unsigned int          isremovedfixings:1; /**< did we removed fixed/aggr/multiaggr variables ? */
   unsigned int          ispropagated:1;     /**< was the constraint propagated with respect to the current bounds ? */
   unsigned int          ispresolved:1;      /**< did we checked for possibilities of upgrading or implicit integer variables ? */
   unsigned int          initialmerge:1;     /**< did we perform an initial merge and clean in presolving yet ? */
#ifdef CHECKIMPLINBILINEAR
   unsigned int          isimpladded:1;      /**< has there been an implication added for a binary variable in a bilinear term? */
#endif
   unsigned int          isgaugeavailable:1; /**< is the gauge function computed? */
   unsigned int          isedavailable:1;    /**< is the eigen decomposition of A available? */

   SCIP_Real             minlinactivity;     /**< sum of minimal activities of all linear terms with finite minimal activity */
   SCIP_Real             maxlinactivity;     /**< sum of maximal activities of all linear terms with finite maximal activity */
   int                   minlinactivityinf;  /**< number of linear terms with infinite minimal activity */
   int                   maxlinactivityinf;  /**< number of linear terms with infinity maximal activity */
   SCIP_INTERVAL         quadactivitybounds; /**< bounds on the activity of the quadratic term, if up to date, otherwise empty interval */
   SCIP_Real             activity;           /**< activity of quadratic function w.r.t. current solution */
   SCIP_Real             lhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */

   int                   linvar_maydecrease; /**< index of a variable in linvars that may be decreased without making any other constraint infeasible, or -1 if none */
   int                   linvar_mayincrease; /**< index of a variable in linvars that may be increased without making any other constraint infeasible, or -1 if none */

   SCIP_VAR**            sepaquadvars;       /**< variables corresponding to quadvarterms to use in separation, only available in solving stage */
   int*                  sepabilinvar2pos;   /**< position of second variable in bilinear terms to use in separation, only available in solving stage */
   SCIP_Real             lincoefsmin;        /**< minimal absolute value of coefficients in linear part, only available in solving stage */
   SCIP_Real             lincoefsmax;        /**< maximal absolute value of coefficients in linear part, only available in solving stage */

   SCIP_Real*            factorleft;         /**< coefficients of left factor if constraint function is factorable */
   SCIP_Real*            factorright;        /**< coefficients of right factor if constraint function is factorable */

   SCIP_Real*            gaugecoefs;         /**< coefficients of the gauge function */
   SCIP_Real             gaugeconst;         /**< constant of the gauge function */
   SCIP_Real*            interiorpoint;      /**< interior point of the region defined by the convex function */
   SCIP_Real             interiorpointval;   /**< function value at interior point */

   SCIP_Real*            eigenvalues;        /**< eigenvalues of A */
   SCIP_Real*            eigenvectors;       /**< orthonormal eigenvectors of A; if A = P D P^T, then eigenvectors is P^T */
   SCIP_Real*            bp;                 /**< stores b * P where b are the linear coefficients of the quadratic vars */
   SCIP_Real             maxnonconvexity;    /**< nonconvexity measure: estimate on largest absolute value of nonconvex eigenvalues */

   SCIP_Bool             isdisaggregated;    /**< has the constraint already been disaggregated? if might happen that more disaggreation would be potentially
                                                  possible, but we reached the maximum number of sparsity components during presolveDisaggregate() */
};

/** quadratic constraint update method */
struct SCIP_QuadConsUpgrade
{
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd));  /**< method to call for upgrading quadratic constraint */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
};
typedef struct SCIP_QuadConsUpgrade SCIP_QUADCONSUPGRADE; /**< quadratic constraint update method */

/** structure to store everything needed for using linear inequalities to improve upon the McCormick relaxation */
struct BilinearEstimator
{
   SCIP_VAR*             x;                 /**< first variable */
   SCIP_VAR*             y;                 /**< second variable */
   SCIP_Real             inequnderest[6];   /**< at most two inequalities that can be used to underestimate xy; stored as (xcoef,ycoef,constant) with xcoef x <= ycoef y + constant */
   SCIP_Real             ineqoverest[6];    /**< at most two inequalities that can be used to overestimate xy; stored as (xcoef,ycoef,constant) with xcoef x <= ycoef y + constant */
   SCIP_Real             maxnonconvexity;   /**< estimate on largest absolute value of nonconvex eigenvalues of all quadratic constraint containing xy */
   int                   ninequnderest;     /**< total number of inequalities for underestimating xy */
   int                   nineqoverest;      /**< total number of inequalities for overestimating xy */
   int                   nunderest;         /**< number of constraints that require to underestimate xy */
   int                   noverest;          /**< number of constraints that require to overestimate xy */

   SCIP_Real             lastimprfac;       /**< last achieved improvement factor */
};
typedef struct BilinearEstimator BILINESTIMATOR;

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   replacebinaryprodlength; /**< length of linear term which when multiplied with a binary variable is replaced by an auxiliary variable and an equivalent linear formulation */
   int                   empathy4and;        /**< how much empathy we have for using the AND constraint handler: 0 avoid always; 1 use sometimes; 2 use as often as possible */
   SCIP_Bool             binreforminitial;   /**< whether to make constraints added due to replacing products with binary variables initial */
   SCIP_Bool             binreformbinaryonly;/**< whether to consider only binary variables when reformulating products with binary variables */
   SCIP_Real             binreformmaxcoef;   /**< factor on 1/feastol to limit coefficients and coef range in linear constraints created by binary reformulation */
   SCIP_Real             cutmaxrange;        /**< maximal range (maximal coef / minimal coef) of a cut in order to be added to LP */
   SCIP_Bool             linearizeheursol;   /**< whether linearizations of convex quadratic constraints should be added to cutpool when some heuristics finds a new solution */
   SCIP_Bool             checkcurvature;     /**< whether functions should be checked for convexity/concavity */
   SCIP_Bool             checkfactorable;    /**< whether functions should be checked to be factorable */
   char                  checkquadvarlocks;  /**< whether quadratic variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction) */
   SCIP_Bool             linfeasshift;       /**< whether to make solutions in check feasible if possible */
   int                   maxdisaggrsize;     /**< maximum number of components when disaggregating a quadratic constraint (<= 1: off) */
   char                  disaggrmergemethod; /**< method on merging blocks in disaggregation */
   int                   maxproprounds;      /**< limit on number of propagation rounds for a single constraint within one round of SCIP propagation during solve */
   int                   maxproproundspresolve; /**< limit on number of propagation rounds for a single constraint within one presolving round */
   SCIP_Real             sepanlpmincont;     /**< minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation */
   SCIP_Bool             enfocutsremovable;  /**< are cuts added during enforcement removable from the LP in the same node? */
   SCIP_Bool             gaugecuts;          /**< should convex quadratics generated strong cuts via gauge function? */
   SCIP_Bool             projectedcuts;      /**< should convex quadratics generated strong cuts via projections? */
   char                  interiorcomputation;/**< how the interior point should be computed: 'a'ny point per constraint, 'm'ost interior per constraint */
   char                  branchscoring;      /**< method to use to compute score of branching candidates */
   int                   enfolplimit;        /**< maximum number of enforcement round before declaring the LP relaxation
                                              * infeasible (-1: no limit); WARNING: if this parameter is not set to -1,
                                              * SCIP might declare sub-optimal solutions optimal or feasible instances
                                              * infeasible; thus, the result returned by SCIP might be incorrect!
                                              */
   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subnlp heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the trysol heuristic, if available */
   SCIP_EVENTHDLR*       eventhdlr;          /**< our handler for variable bound change events */
   int                   newsoleventfilterpos; /**< filter position of new solution event handler, if caught */
   SCIP_Bool             sepanlp;            /**< where linearization of the NLP relaxation solution added? */
   SCIP_NODE*            lastenfonode;       /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenforounds;        /**< counter on number of enforcement rounds for the current node */
   SCIP_QUADCONSUPGRADE** quadconsupgrades;  /**< quadratic constraint upgrade methods for specializing quadratic constraints */
   int                   quadconsupgradessize; /**< size of quadconsupgrade array */
   int                   nquadconsupgrades;  /**< number of quadratic constraint upgrade methods */

   BILINESTIMATOR*       bilinestimators;    /**< array containing all required information for using stronger estimators for each bilinear term in all quadratic constraints */
   int                   nbilinterms;        /**< number of bilinear terms in all quadratic constraints */

   SCIP_Bool             usebilinineqbranch; /**< should linear inequalities be considered when computing the branching scores for bilinear terms? */
   SCIP_Bool             storedbilinearterms; /**< did we already try to store all bilinear terms? */

   SCIP_Real             minscorebilinterms; /**< minimal required score in order to use linear inequalities for tighter bilinear relaxations */
   SCIP_Real             mincurvcollectbilinterms;/**< minimal curvature of constraints to be considered when returning bilinear terms to other plugins */
   int                   bilinineqmaxseparounds; /**< maximum number of separation rounds to use linear inequalities for the bilinear term relaxation in a local node */
};


/*
 * local methods for managing quadratic constraint update methods
 */


/** checks whether a quadratic constraint upgrade method has already be registered */
static
SCIP_Bool conshdlrdataHasUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd)),  /**< method to call for upgrading quadratic constraint */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(quadconsupgd != NULL);
   assert(conshdlrname != NULL);

   for( i = conshdlrdata->nquadconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->quadconsupgrades[i]->quadconsupgd == quadconsupgd )
      {
         SCIPwarningMessage(scip, "Try to add already known upgrade message for constraint handler <%s>.\n", conshdlrname);
         return TRUE;
      }
   }

   return FALSE;
}

/*
 * Local methods
 */

/** translate from one value of infinity to another
 * 
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/** catches variable bound change events on a linear variable in a quadratic constraint */
static
SCIP_RETCODE catchLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
   )
{
   SCIP_CONSDATA*  consdata;
   SCIP_QUADVAREVENTDATA* eventdata;
   SCIP_EVENTTYPE  eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);
   assert(consdata->lineventdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventdata) );

   eventdata->cons = cons;
   eventdata->varidx = linvarpos;

   eventtype = SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_GBDCHANGED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->linvars[linvarpos], eventtype, eventhdlr, (SCIP_EVENTDATA*)eventdata, &eventdata->filterpos) );

   consdata->lineventdata[linvarpos] = eventdata;

   /* invalidate activity information
    * NOTE: It could happen that a constraint gets temporary deactivated and some variable bounds change. In this case
    *       we do not recognize those bound changes with the variable events and thus we have to recompute the activities.
    */
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   return SCIP_OKAY;
}

/** drops variable bound change events on a linear variable in a quadratic constraint */
static
SCIP_RETCODE dropLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE  eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);
   assert(consdata->lineventdata != NULL);
   assert(consdata->lineventdata[linvarpos] != NULL);
   assert(consdata->lineventdata[linvarpos]->cons == cons);
   assert(consdata->lineventdata[linvarpos]->varidx == linvarpos);
   assert(consdata->lineventdata[linvarpos]->filterpos >= 0);

   eventtype = SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_GBDCHANGED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest
       * since we also want to keep activities in consdata up-to-date, we also need to know when the corresponding bound is relaxed */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->linvars[linvarpos], eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata->lineventdata[linvarpos], consdata->lineventdata[linvarpos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->lineventdata[linvarpos]);  /*lint !e866 */

   return SCIP_OKAY;
}

/** catches variable bound change events on a quadratic variable in a quadratic constraint */
static
SCIP_RETCODE catchQuadVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   quadvarpos          /**< position of variable in quadratic variables array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_QUADVAREVENTDATA* eventdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(quadvarpos >= 0);
   assert(quadvarpos < consdata->nquadvars);
   assert(consdata->quadvarterms[quadvarpos].eventdata == NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventdata) );

   eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_GBDCHANGED;
#ifdef CHECKIMPLINBILINEAR
   eventtype |= SCIP_EVENTTYPE_IMPLADDED;
#endif
   eventdata->cons = cons;
   eventdata->varidx   = -quadvarpos-1;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->quadvarterms[quadvarpos].var, eventtype, eventhdlr, (SCIP_EVENTDATA*)eventdata, &eventdata->filterpos) );

   consdata->quadvarterms[quadvarpos].eventdata = eventdata;

   /* invalidate activity information
    * NOTE: It could happen that a constraint gets temporary deactivated and some variable bounds change. In this case
    *       we do not recognize those bound changes with the variable events and thus we have to recompute the activities.
    */
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   return SCIP_OKAY;
}

/** catches variable bound change events on a quadratic variable in a quadratic constraint */
static
SCIP_RETCODE dropQuadVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   quadvarpos          /**< position of variable in quadratic variables array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(quadvarpos >= 0);
   assert(quadvarpos < consdata->nquadvars);
   assert(consdata->quadvarterms[quadvarpos].eventdata != NULL);
   assert(consdata->quadvarterms[quadvarpos].eventdata->cons == cons);
   assert(consdata->quadvarterms[quadvarpos].eventdata->varidx == -quadvarpos-1);
   assert(consdata->quadvarterms[quadvarpos].eventdata->filterpos >= 0);

   eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_GBDCHANGED;
#ifdef CHECKIMPLINBILINEAR
   eventtype |= SCIP_EVENTTYPE_IMPLADDED;
#endif

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->quadvarterms[quadvarpos].var, eventtype, eventhdlr, (SCIP_EVENTDATA*)consdata->quadvarterms[quadvarpos].eventdata, consdata->quadvarterms[quadvarpos].eventdata->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->quadvarterms[quadvarpos].eventdata);

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
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lineventdata == NULL);

   /* we will update isremovedfixings, so reset it to TRUE first */
   consdata->isremovedfixings = TRUE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize) );
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( catchLinearVarEvents(scip, eventhdlr, cons, i) );

      var = consdata->linvars[i];
      consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var)
         && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      assert(consdata->quadvarterms[i].eventdata == NULL);

      SCIP_CALL( catchQuadVarEvents(scip, eventhdlr, cons, i) );

      var = consdata->quadvarterms[i].var;
      consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var)
         && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   consdata->ispropagated = FALSE;

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

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->lineventdata != NULL )
   {
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         if( consdata->lineventdata[i] != NULL )
         {
            SCIP_CALL( dropLinearVarEvents(scip, eventhdlr, cons, i) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize);
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->quadvarterms[i].eventdata != NULL )
      {
         SCIP_CALL( dropQuadVarEvents(scip, eventhdlr, cons, i) );
      }
   }

   return SCIP_OKAY;
}

/** locks a linear variable in a constraint */
static
SCIP_RETCODE lockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to lock a variable */
   SCIP_VAR*             var,                /**< variable to lock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** unlocks a linear variable in a constraint */
static
SCIP_RETCODE unlockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to unlock a variable */
   SCIP_VAR*             var,                /**< variable to unlock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** locks a quadratic variable in a constraint */
static
SCIP_RETCODE lockQuadraticVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to lock a variable */
   SCIP_VAR*             var                 /**< variable to lock */
   )
{
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** unlocks a quadratic variable in a constraint */
static
SCIP_RETCODE unlockQuadraticVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to unlock a variable */
   SCIP_VAR*             var                 /**< variable to unlock */
   )
{
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** computes the minimal and maximal activity for the linear part in a constraint data
 *
 *  Only sums up terms that contribute finite values.
 *  Gives the number of terms that contribute infinite values.
 *  Only computes those activities where the corresponding side of the constraint is finite.
 */
static
void consdataUpdateLinearActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             intervalinfty       /**< infinity value used in interval operations */
   )
{  /*lint --e{666}*/
   SCIP_ROUNDMODE prevroundmode;
   int       i;
   SCIP_Real bnd;

   assert(scip != NULL);
   assert(consdata != NULL);

   /* if variable bounds are not strictly consistent, then the activity update methods may yield inconsistent activities
    * in this case, we also recompute the activities
    */
   if( consdata->minlinactivity != SCIP_INVALID && consdata->maxlinactivity != SCIP_INVALID &&  /*lint !e777 */
      (consdata->minlinactivityinf > 0 || consdata->maxlinactivityinf > 0 || consdata->minlinactivity <= consdata->maxlinactivity) )
   {
      /* activities should be up-to-date */
      assert(consdata->minlinactivityinf >= 0);
      assert(consdata->maxlinactivityinf >= 0);
      return;
   }

   consdata->minlinactivityinf = 0;
   consdata->maxlinactivityinf = 0;

   /* if lhs is -infinite, then we do not compute a maximal activity, so we set it to  infinity
    * if rhs is  infinite, then we do not compute a minimal activity, so we set it to -infinity
    */
   consdata->minlinactivity = SCIPisInfinity(scip,  consdata->rhs) ? -intervalinfty : 0.0;
   consdata->maxlinactivity = SCIPisInfinity(scip, -consdata->lhs) ?  intervalinfty : 0.0;

   if( consdata->nlinvars == 0 )
      return;

   /* if the activities computed here should be still up-to-date after bound changes,
    * variable events need to be caught */
   assert(consdata->lineventdata != NULL);

   prevroundmode = SCIPintervalGetRoundingMode();

   if( !SCIPisInfinity(scip,  consdata->rhs) )
   {
      /* compute minimal activity only if there is a finite right hand side */
      SCIPintervalSetRoundingModeDownwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(consdata->lineventdata[i] != NULL);
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = MIN(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         else
         {
            bnd = MAX(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         consdata->minlinactivity += consdata->lincoefs[i] * bnd;
      }
   }

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* compute maximal activity only if there is a finite left hand side */
      SCIPintervalSetRoundingModeUpwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(consdata->lineventdata[i] != NULL);
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = MAX(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         else
         {
            bnd = MIN(SCIPvarGetLbLocal(consdata->linvars[i]), SCIPvarGetUbLocal(consdata->linvars[i]));
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         consdata->maxlinactivity += consdata->lincoefs[i] * bnd;
      }
   }

   SCIPintervalSetRoundingMode(prevroundmode);

   assert(consdata->minlinactivityinf > 0 || consdata->maxlinactivityinf > 0 || consdata->minlinactivity <= consdata->maxlinactivity);
}

/** update the linear activities after a change in the lower bound of a variable */
static
void consdataUpdateLinearActivityLbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with lower bounds at infinity */
   assert(!SCIPisInfinity(scip, oldbnd));
   assert(!SCIPisInfinity(scip, newbnd));

   /* @todo since we check the linear activity for consistency later anyway, we may skip changing the rounding mode here */

   /* assume lhs <= a*x + y <= rhs, then the following bound changes can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */

   if( coef > 0.0 )
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->minlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->maxlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** update the linear activities after a change in the upper bound of a variable */
static
void consdataUpdateLinearActivityUbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with upper bounds at -infinity */
   assert(!SCIPisInfinity(scip, -oldbnd));
   assert(!SCIPisInfinity(scip, -newbnd));

   /* @todo since we check the linear activity for consistency later anyway, we may skip changing the rounding mode here */

   /* assume lhs <= a*x + y <= rhs, then the following bound changes can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */

   if( coef > 0.0 )
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->maxlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777 */
         return;

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->minlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** returns whether a quadratic variable domain can be reduced to its lower or upper bound; this is the case if the
 *  quadratic variable is in just one single quadratic constraint and (sqrcoef > 0 and LHS = -infinity), or
 *  (sqrcoef < 0 and RHS = +infinity) hold
 */
static
SCIP_Bool hasQuadvarHpProperty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   idx                 /**< index of quadratic variable */
   )
{
   SCIP_VAR* var;
   SCIP_Real quadcoef;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(idx >= 0 && idx < consdata->nquadvars);

   var = consdata->quadvarterms[idx].var;
   assert(var != NULL);

   quadcoef = consdata->quadvarterms[idx].sqrcoef;
   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   return SCIPvarGetNLocksDown(var) == 1 && SCIPvarGetNLocksUp(var) == 1 && SCIPisZero(scip, SCIPvarGetObj(var))
      && SCIPvarGetType(var) != SCIP_VARTYPE_BINARY && ((quadcoef < 0.0 && !haslhs) || (quadcoef > 0.0 && !hasrhs));
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;
   int varidx;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   cons = ((SCIP_QUADVAREVENTDATA*)eventdata)->cons;
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   varidx = ((SCIP_QUADVAREVENTDATA*)eventdata)->varidx;
   assert(varidx <  0 ||  varidx   < consdata->nlinvars);
   assert(varidx >= 0 || -varidx-1 < consdata->nquadvars);

   eventtype = SCIPeventGetType(event);

   /* process local bound changes */
   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      if( varidx < 0 )
      {
         /* mark activity bounds for quad term as not up to date anymore */
         SCIPintervalSetEmpty(&consdata->quadactivitybounds);
      }
      else
      {
         /* update activity bounds for linear terms */
         if( eventtype & SCIP_EVENTTYPE_LBCHANGED )
            consdataUpdateLinearActivityLbChange(scip, consdata, consdata->lincoefs[varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));
         else
            consdataUpdateLinearActivityUbChange(scip, consdata, consdata->lincoefs[varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));
      }

      if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
      {
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
         consdata->ispropagated = FALSE;
      }
   }

   /* process global bound changes */
   if( eventtype & SCIP_EVENTTYPE_GBDCHANGED )
   {
      SCIP_VAR* var;

      var = varidx < 0 ? consdata->quadvarterms[-varidx-1].var : consdata->linvars[varidx];
      assert(var != NULL);

      if( varidx < 0 )
      {
         SCIP_QUADVARTERM* quadvarterm;

         quadvarterm = &consdata->quadvarterms[-varidx-1];

         /* if an integer variable x with a x^2 is tightened to [0,1], then we can replace the x^2 by x, which is done in mergeAndCleanQuadVarTerms()
          * we currently do this only if the binary variable does not show up in any bilinear terms
          */
         if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING && SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) &&
            quadvarterm->sqrcoef != 0.0 && quadvarterm->nadjbilin == 0 )
         {
            consdata->quadvarsmerged = FALSE;
            consdata->initialmerge = FALSE;
         }
      }

      if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
         consdata->isremovedfixings = FALSE;
   }

   /* process variable fixing event */
   if( eventtype & SCIP_EVENTTYPE_VARFIXED )
   {
      consdata->isremovedfixings = FALSE;
   }

#ifdef CHECKIMPLINBILINEAR
   if( eventtype & SCIP_EVENTTYPE_IMPLADDED )
   {
      assert(varidx < 0); /* we catch impladded events only for quadratic variables */
      /* if variable is binary (quite likely if an implication has been added) and occurs in a bilinear term, then mark that we should check implications */
      if( SCIPvarIsBinary(SCIPeventGetVar(event)) && consdata->quadvarterms[-varidx-1].nadjbilin > 0 )
         consdata->isimpladded = TRUE;
   }
#endif

   return SCIP_OKAY;
}

/** ensures, that linear vars and coefs arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureLinearVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nlinvars <= consdata->linvarssize);

   if( num > consdata->linvarssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->linvars,  consdata->linvarssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lincoefs, consdata->linvarssize, newsize) );
      if( consdata->lineventdata != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize, newsize) );
      }
      consdata->linvarssize = newsize;
   }
   assert(num <= consdata->linvarssize);

   return SCIP_OKAY;
}

/** ensures, that quadratic variable terms array can store at least num entries */
static
SCIP_RETCODE consdataEnsureQuadVarTermsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nquadvars <= consdata->quadvarssize);

   if( num > consdata->quadvarssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->quadvarterms, consdata->quadvarssize, newsize) );
      consdata->quadvarssize = newsize;
   }
   assert(num <= consdata->quadvarssize);

   return SCIP_OKAY;
}

/** ensures, that adjacency array can store at least num entries */
static
SCIP_RETCODE consdataEnsureAdjBilinSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_QUADVARTERM*     quadvarterm,        /**< quadratic variable term */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(quadvarterm != NULL);
   assert(quadvarterm->nadjbilin <= quadvarterm->adjbilinsize);

   if( num > quadvarterm->adjbilinsize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &quadvarterm->adjbilin, quadvarterm->adjbilinsize, newsize) );
      quadvarterm->adjbilinsize = newsize;
   }
   assert(num <= quadvarterm->adjbilinsize);

   return SCIP_OKAY;
}

/** ensures, that bilinear term arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureBilinSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nbilinterms <= consdata->bilintermssize);

   if( num > consdata->bilintermssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->bilinterms, consdata->bilintermssize, newsize) );
      consdata->bilintermssize = newsize;
   }
   assert(num <= consdata->bilintermssize);

   return SCIP_OKAY;
}

/** creates empty constraint data structure */
static
SCIP_RETCODE consdataCreateEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< a buffer to store pointer to new constraint data */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->lhs = -SCIPinfinity(scip);
   (*consdata)->rhs =  SCIPinfinity(scip);

   (*consdata)->linvarssorted  = TRUE;
   (*consdata)->linvarsmerged  = TRUE;
   (*consdata)->quadvarssorted = TRUE;
   (*consdata)->quadvarsmerged = TRUE;
   (*consdata)->bilinsorted    = TRUE;
   (*consdata)->bilinmerged    = TRUE;

   (*consdata)->isremovedfixings = TRUE;
   (*consdata)->ispropagated     = TRUE;
   (*consdata)->initialmerge     = FALSE;

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   (*consdata)->isgaugeavailable = FALSE;
   (*consdata)->isedavailable = FALSE;

   return SCIP_OKAY;
}

/** creates constraint data structure */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< a buffer to store pointer to new constraint data */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            lincoefs,           /**< array of coefficients of linear variables */
   int                   nquadvars,          /**< number of quadratic variables */
   SCIP_QUADVARTERM*     quadvarterms,       /**< array of quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms */
   SCIP_BILINTERM*       bilinterms,         /**< array of bilinear terms */
   SCIP_Bool             capturevars         /**< whether we should capture variables */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);

   assert(nlinvars == 0 || linvars  != NULL);
   assert(nlinvars == 0 || lincoefs != NULL);
   assert(nquadvars == 0 || quadvarterms != NULL);
   assert(nbilinterms == 0 || bilinterms != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;

   if( nlinvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->linvars,  linvars,  nlinvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->lincoefs, lincoefs, nlinvars) );
      (*consdata)->nlinvars = nlinvars;
      (*consdata)->linvarssize = nlinvars;

      if( capturevars )
         for( i = 0; i < nlinvars; ++i )
         {
            SCIP_CALL( SCIPcaptureVar(scip, linvars[i]) );
         }
   }
   else
   {
      (*consdata)->linvarssorted = TRUE;
      (*consdata)->linvarsmerged = TRUE;
      (*consdata)->minlinactivity = 0.0;
      (*consdata)->maxlinactivity = 0.0;
      (*consdata)->minlinactivityinf = 0;
      (*consdata)->maxlinactivityinf = 0;
   }

   if( nquadvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->quadvarterms, quadvarterms, nquadvars) );

      for( i = 0; i < nquadvars; ++i )
      {
         (*consdata)->quadvarterms[i].eventdata = NULL;
         if( quadvarterms[i].nadjbilin )
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->quadvarterms[i].adjbilin, quadvarterms[i].adjbilin, quadvarterms[i].nadjbilin) );
            (*consdata)->quadvarterms[i].adjbilinsize = quadvarterms[i].nadjbilin;
         }
         else
         {
            assert((*consdata)->quadvarterms[i].nadjbilin == 0);
            (*consdata)->quadvarterms[i].adjbilin = NULL;
            (*consdata)->quadvarterms[i].adjbilinsize = 0;
         }
         if( capturevars )
         {
            SCIP_CALL( SCIPcaptureVar(scip, quadvarterms[i].var) );
         }
      }

      (*consdata)->nquadvars = nquadvars;
      (*consdata)->quadvarssize = nquadvars;
      SCIPintervalSetEmpty(&(*consdata)->quadactivitybounds);
   }
   else
   {
      (*consdata)->quadvarssorted = TRUE;
      (*consdata)->quadvarsmerged = TRUE;
      SCIPintervalSet(&(*consdata)->quadactivitybounds, 0.0);
   }

   if( nbilinterms > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->bilinterms, bilinterms, nbilinterms) );
      (*consdata)->nbilinterms = nbilinterms;
      (*consdata)->bilintermssize = nbilinterms;
   }
   else
   {
      (*consdata)->bilinsorted = TRUE;
      (*consdata)->bilinmerged = TRUE;
   }

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->activity = SCIP_INVALID;
   (*consdata)->lhsviol  = SCIPisInfinity(scip, -lhs) ? 0.0 : SCIP_INVALID;
   (*consdata)->rhsviol  = SCIPisInfinity(scip,  rhs) ? 0.0 : SCIP_INVALID;

   (*consdata)->isgaugeavailable = FALSE;

   return SCIP_OKAY;
}

/** frees constraint data structure */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data to free */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free sepa arrays, may exists if constraint is deleted in solving stage */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->sepaquadvars,     (*consdata)->nquadvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->sepabilinvar2pos, (*consdata)->nbilinterms);

   /* release linear variables and free linear part */
   if( (*consdata)->linvarssize > 0 )
   {
      for( i = 0; i < (*consdata)->nlinvars; ++i )
      {
         assert((*consdata)->lineventdata == NULL || (*consdata)->lineventdata[i] == NULL); /* variable events should have been dropped earlier */
         SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->linvars[i]) );
      }
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->linvars,  (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->lincoefs, (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lineventdata, (*consdata)->linvarssize);
   }
   assert((*consdata)->linvars == NULL);
   assert((*consdata)->lincoefs == NULL);
   assert((*consdata)->lineventdata == NULL);

   /* release quadratic variables and free quadratic variable term part */
   for( i = 0; i < (*consdata)->nquadvars; ++i )
   {
      assert((*consdata)->quadvarterms[i].eventdata == NULL); /* variable events should have been dropped earlier */
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->quadvarterms[i].adjbilin, (*consdata)->quadvarterms[i].adjbilinsize);
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->quadvarterms[i].var) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->quadvarterms, (*consdata)->quadvarssize);

   /* free bilinear terms */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->bilinterms, (*consdata)->bilintermssize);

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   /* free interior point information, may exists if constraint is deleted in solving stage */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->interiorpoint, (*consdata)->nquadvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->gaugecoefs, (*consdata)->nquadvars);

   /* free eigen decomposition information */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->eigenvalues, (*consdata)->nquadvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->eigenvectors, (int)((*consdata)->nquadvars*(*consdata)->nquadvars));
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->bp, (*consdata)->nquadvars);

   /* free unique indices of bilinear terms array */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->bilintermsidx, (*consdata)->nbilinterms);

   SCIPfreeBlockMemory(scip, consdata);
   *consdata = NULL;

   return SCIP_OKAY;
}

/** sorts linear part of constraint data */
static
void consdataSortLinearVars(
   SCIP_CONSDATA*        consdata            /**< quadratic constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->linvarssorted )
      return;

   if( consdata->nlinvars <= 1 )
   {
      consdata->linvarssorted = TRUE;
      return;
   }

   if( consdata->lineventdata == NULL )
   {
      SCIPsortPtrReal((void**)consdata->linvars, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);
   }
   else
   {
      int i;

      SCIPsortPtrPtrReal((void**)consdata->linvars, (void**)consdata->lineventdata, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);

      /* update variable indices in event data */
      for( i = 0; i < consdata->nlinvars; ++i )
         if( consdata->lineventdata[i] != NULL )
            consdata->lineventdata[i]->varidx = i;
   }

   consdata->linvarssorted = TRUE;
}

#ifdef SCIP_DISABLED_CODE /* no-one needs this routine currently */
/** returns the position of variable in the linear coefficients array of a constraint, or -1 if not found */
static
int consdataFindLinearVar(
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   SCIP_VAR*             var                 /**< variable to search for */
   )
{
   int pos;

   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->nlinvars == 0 )
      return -1;

   consdataSortLinearVars(consdata);

   if( !SCIPsortedvecFindPtr((void**)consdata->linvars, SCIPvarComp, (void*)var, consdata->nlinvars, &pos) )
      pos = -1;

   return pos;
}
#endif

/** index comparison method for quadratic variable terms: compares two indices of the quadratic variable set in the quadratic constraint */
static
SCIP_DECL_SORTINDCOMP(quadVarTermComp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nquadvars);
   assert(0 <= ind2 && ind2 < consdata->nquadvars);

   return SCIPvarCompare(consdata->quadvarterms[ind1].var, consdata->quadvarterms[ind2].var);
}

/** sorting of quadratic variable terms */
static
SCIP_RETCODE consdataSortQuadVarTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< quadratic constraint data */
   )
{
   int* perm;
   int  i;
   int  nexti;
   int  v;
   SCIP_QUADVARTERM quadterm;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->quadvarssorted )
      return SCIP_OKAY;

   if( consdata->nquadvars == 0 )
   {
      consdata->quadvarssorted = TRUE;
      return SCIP_OKAY;
   }

   /* get temporary memory to store the sorted permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nquadvars) );

   /* call bubble sort */
   SCIPsort(perm, quadVarTermComp, (void*)consdata, consdata->nquadvars);

   /* permute the quadratic variable terms according to the resulting permutation */
   for( v = 0; v < consdata->nquadvars; ++v )
   {
      if( perm[v] != v )
      {
         quadterm = consdata->quadvarterms[v];

         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < consdata->nquadvars);
            assert(perm[i] != i);
            consdata->quadvarterms[i] = consdata->quadvarterms[perm[i]];
            if( consdata->quadvarterms[i].eventdata != NULL )
            {
               consdata->quadvarterms[i].eventdata->varidx = -i-1;
            }
            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         consdata->quadvarterms[i] = quadterm;
         if( consdata->quadvarterms[i].eventdata != NULL )
         {
            consdata->quadvarterms[i].eventdata->varidx = -i-1;
         }
         perm[i] = i;
      }
   }
   consdata->quadvarssorted = TRUE;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}

/** returns the position of variable in the quadratic variable terms array of a constraint, or -1 if not found */
static
SCIP_RETCODE consdataFindQuadVarTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< quadratic constraint data */
   SCIP_VAR*             var,                /**< variable to search for */
   int*                  pos                 /**< buffer where to store position of var in quadvarterms array, or -1 if not found */
   )
{
   int left;
   int right;
   int cmpres;

   assert(consdata != NULL);
   assert(var != NULL);
   assert(pos != NULL);

   if( consdata->nquadvars == 0 )
   {
      *pos = -1;
      return SCIP_OKAY;
   }

   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   left = 0;
   right = consdata->nquadvars - 1;
   while( left <= right )
   {
      int middle;

      middle = (left+right)/2;
      assert(0 <= middle && middle < consdata->nquadvars);

      cmpres = SCIPvarCompare(var, consdata->quadvarterms[middle].var);

      if( cmpres < 0 )
         right = middle - 1;
      else if( cmpres > 0 )
         left  = middle + 1;
      else
      {
         *pos = middle;
         return SCIP_OKAY;
      }
   }
   assert(left == right+1);

   *pos = -1;

   return SCIP_OKAY;
}

/** index comparison method for bilinear terms: compares two index pairs of the bilinear term set in the quadratic constraint */
static
SCIP_DECL_SORTINDCOMP(bilinTermComp)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;
   int var1cmp;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nbilinterms);
   assert(0 <= ind2 && ind2 < consdata->nbilinterms);

   var1cmp = SCIPvarCompare(consdata->bilinterms[ind1].var1, consdata->bilinterms[ind2].var1);
   if( var1cmp != 0 )
      return var1cmp;

   return SCIPvarCompare(consdata->bilinterms[ind1].var2, consdata->bilinterms[ind2].var2);
}

#ifndef NDEBUG
/** checks if all bilinear terms are sorted correctly */
static
SCIP_Bool consdataCheckBilinTermsSort(
   SCIP_CONSDATA* consdata
   )
{
   int i;

   assert(consdata != NULL);

   /* nothing to check if the bilinear terms have not been sorted yet */
   if( !consdata->bilinsorted )
      return TRUE;

   for( i = 0; i < consdata->nbilinterms - 1; ++i )
   {
      if( bilinTermComp(consdata, i, i+1) > 0 )
         return FALSE;
   }
   return TRUE;
}
#endif

/** sorting of bilinear terms */
static
SCIP_RETCODE consdataSortBilinTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< quadratic constraint data */
   )
{
   int* perm;
   int* invperm;
   int  i;
   int  nexti;
   int  v;
   SCIP_BILINTERM bilinterm;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->bilinsorted )
      return SCIP_OKAY;

   if( consdata->nbilinterms == 0 )
   {
      consdata->bilinsorted = TRUE;
      return SCIP_OKAY;
   }

   /* get temporary memory to store the sorted permutation and the inverse permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm,    consdata->nbilinterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &invperm, consdata->nbilinterms) );

   /* call bubble sort */
   SCIPsort(perm, bilinTermComp, (void*)consdata, consdata->nbilinterms);

   /* compute inverted permutation */
   for( v = 0; v < consdata->nbilinterms; ++v )
   {
      assert(0 <= perm[v] && perm[v] < consdata->nbilinterms);
      invperm[perm[v]] = v;
   }

   /* permute the bilinear terms according to the resulting permutation */
   for( v = 0; v < consdata->nbilinterms; ++v )
   {
      if( perm[v] != v )
      {
         bilinterm = consdata->bilinterms[v];

         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < consdata->nbilinterms);
            assert(perm[i] != i);
            consdata->bilinterms[i] = consdata->bilinterms[perm[i]];
            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         consdata->bilinterms[i] = bilinterm;
         perm[i] = i;
      }
   }

   /* update the adjacency information in the quadratic variable terms */
   for( v = 0; v < consdata->nquadvars; ++v )
      for( i = 0; i < consdata->quadvarterms[v].nadjbilin; ++i )
         consdata->quadvarterms[v].adjbilin[i] = invperm[consdata->quadvarterms[v].adjbilin[i]];

   consdata->bilinsorted = TRUE;
   assert(consdataCheckBilinTermsSort(consdata));

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &invperm);
   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}

/** moves a linear variable from one position to another */
static
void consdataMoveLinearVar(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   oldpos,             /**< position of variable that shall be moved */
   int                   newpos              /**< new position of variable */
   )
{
   assert(consdata != NULL);
   assert(oldpos >= 0);
   assert(oldpos < consdata->nlinvars);
   assert(newpos >= 0);
   assert(newpos < consdata->linvarssize);

   if( newpos == oldpos )
      return;

   consdata->linvars [newpos] = consdata->linvars [oldpos];
   consdata->lincoefs[newpos] = consdata->lincoefs[oldpos];

   if( consdata->lineventdata != NULL )
   {
      assert(newpos >= consdata->nlinvars || consdata->lineventdata[newpos] == NULL);

      consdata->lineventdata[newpos] = consdata->lineventdata[oldpos];
      consdata->lineventdata[newpos]->varidx = newpos;

      consdata->lineventdata[oldpos] = NULL;
   }

   consdata->linvarssorted = FALSE;
}   

/** moves a quadratic variable from one position to another */
static
void consdataMoveQuadVarTerm(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   oldpos,             /**< position of variable that shall be moved */
   int                   newpos              /**< new position of variable */
   )
{
   assert(consdata != NULL);
   assert(oldpos >= 0);
   assert(oldpos < consdata->nquadvars);
   assert(newpos >= 0);
   assert(newpos < consdata->quadvarssize);

   if( newpos == oldpos )
      return;

   assert(newpos >= consdata->nquadvars || consdata->quadvarterms[newpos].eventdata == NULL);

   consdata->quadvarterms[newpos] = consdata->quadvarterms[oldpos];

   if( consdata->quadvarterms[newpos].eventdata != NULL )
   {
      consdata->quadvarterms[newpos].eventdata->varidx = -newpos-1;
      consdata->quadvarterms[oldpos].eventdata = NULL;
   }

   consdata->quadvarssorted = FALSE;
}   

/** adds linear coefficient in quadratic constraint */
static
SCIP_RETCODE addLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             coef                /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);

   /* ignore coefficient if it is nearly zero */
   if( SCIPisZero(scip, coef) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars+1) );
   consdata->linvars [consdata->nlinvars] = var;
   consdata->lincoefs[consdata->nlinvars] = coef;

   ++consdata->nlinvars;

   /* catch variable events */
   if( SCIPconsIsEnabled(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      assert(consdata->lineventdata != NULL);
      consdata->lineventdata[consdata->nlinvars-1] = NULL;

      /* catch bound change events of variable */
      SCIP_CALL( catchLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, consdata->nlinvars-1) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockLinearVariable(scip, cons, var, coef) );

   /* capture new variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;
   consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var)
      && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   if( consdata->nlinvars == 1 )
      consdata->linvarssorted = TRUE;
   else
      consdata->linvarssorted = consdata->linvarssorted && (SCIPvarCompare(consdata->linvars[consdata->nlinvars-2], consdata->linvars[consdata->nlinvars-1]) == -1);
   /* always set too FALSE since the new linear variable should be checked if already existing as quad var term */ 
   consdata->linvarsmerged = FALSE;

   return SCIP_OKAY;
}

/** deletes linear coefficient at given position from quadratic constraint data */
static
SCIP_RETCODE delLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nlinvars);

   var  = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );

   /* if we catch variable events, drop the events on the variable */
   if( consdata->lineventdata != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      SCIP_CALL( dropLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->linvars[pos]) );

   /* move the last variable to the free slot */
   consdataMoveLinearVar(consdata, consdata->nlinvars-1, pos);

   --consdata->nlinvars;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved  = FALSE;

   return SCIP_OKAY;
}

/** changes linear coefficient value at given position of quadratic constraint */
static
SCIP_RETCODE chgLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos,                /**< position of linear coefficient to change */
   SCIP_Real             newcoef             /**< new value of linear coefficient */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisZero(scip, newcoef));

   conshdlrdata = NULL;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos);
   assert(pos < consdata->nlinvars);
   assert(!SCIPisZero(scip, newcoef));

   var = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* if necessary, remove the rounding locks and event catching of the variable */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         assert(SCIPconsIsTransformed(cons));

         /* remove rounding locks for variable with old coefficient */
         SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );
      }

      if( consdata->lineventdata[pos] != NULL )
      {
         /* get event handler */
         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         assert(conshdlrdata->eventhdlr != NULL);

         /* drop bound change events of variable */
         SCIP_CALL( dropLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
      }
   }

   /* change the coefficient */
   consdata->lincoefs[pos] = newcoef;

   /* if necessary, install the rounding locks and event catching of the variable again */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         /* install rounding locks for variable with new coefficient */
         SCIP_CALL( lockLinearVariable(scip, cons, var, newcoef) );
      }

      if( conshdlrdata != NULL )
      {
         assert(SCIPconsIsEnabled(cons));

         /* catch bound change events of variable */
         SCIP_CALL( catchLinearVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
      }
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;

   return SCIP_OKAY;
}

/** adds quadratic variable term to quadratic constraint */
static
SCIP_RETCODE addQuadVarTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             lincoef,            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef             /**< square coefficient of variable */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;
   SCIP_QUADVARTERM* quadvarterm;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureQuadVarTermsSize(scip, consdata, consdata->nquadvars+1) );

   quadvarterm = &consdata->quadvarterms[consdata->nquadvars];
   quadvarterm->var       = var;
   quadvarterm->lincoef   = lincoef;
   quadvarterm->sqrcoef   = sqrcoef;
   quadvarterm->adjbilinsize = 0;
   quadvarterm->nadjbilin = 0;
   quadvarterm->adjbilin  = NULL;
   quadvarterm->eventdata = NULL;

   ++consdata->nquadvars;

   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* catch variable events, if we do so */
   if( SCIPconsIsEnabled(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* catch bound change events of variable */
      SCIP_CALL( catchQuadVarEvents(scip, conshdlrdata->eventhdlr, cons, consdata->nquadvars-1) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockQuadraticVariable(scip, cons, var) );

   consdata->ispropagated = FALSE;
   consdata->ispresolved  = FALSE;
   consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var)
      && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   if( consdata->nquadvars == 1 )
      consdata->quadvarssorted = TRUE;
   else
      consdata->quadvarssorted = consdata->quadvarssorted &&
         (SCIPvarCompare(consdata->quadvarterms[consdata->nquadvars-2].var, consdata->quadvarterms[consdata->nquadvars-1].var) == -1);
   /* also set to FALSE if nquadvars == 1, since the new variable should be checked for linearity and other stuff in mergeAndClean ... */ 
   consdata->quadvarsmerged = FALSE;

   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** deletes quadratic variable term at given position from quadratic constraint data */
static
SCIP_RETCODE delQuadVarTermPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos                 /**< position of term to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nquadvars);

   var = consdata->quadvarterms[pos].var;
   assert(var != NULL);
   assert(consdata->quadvarterms[pos].nadjbilin == 0);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockQuadraticVariable(scip, cons, var) );

   /* if we catch variable events, drop the events on the variable */
   if( consdata->quadvarterms[pos].eventdata != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      SCIP_CALL( dropQuadVarEvents(scip, conshdlrdata->eventhdlr, cons, pos) );
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->quadvarterms[pos].var) );

   /* free adjacency array */
   SCIPfreeBlockMemoryArrayNull(scip, &consdata->quadvarterms[pos].adjbilin, consdata->quadvarterms[pos].adjbilinsize);

   /* move the last variable term to the free slot */
   consdataMoveQuadVarTerm(consdata, consdata->nquadvars-1, pos);

   --consdata->nquadvars;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = FALSE;
   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** replace variable in quadratic variable term at given position of quadratic constraint data
 *
 * Allows to replace x by coef*y+offset, thereby maintaining linear and square coefficients and bilinear terms.
 */
static
SCIP_RETCODE replaceQuadVarTermPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   pos,                /**< position of term to replace */
   SCIP_VAR*             var,                /**< new variable */
   SCIP_Real             coef,               /**< linear coefficient of new variable */
   SCIP_Real             offset              /**< offset of new variable */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real constant;

   int i;
   SCIP_VAR* var2;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(pos >= 0);
   assert(pos <  consdata->nquadvars);

   quadvarterm = &consdata->quadvarterms[pos];

   /* remove rounding locks for old variable */
   SCIP_CALL( unlockQuadraticVariable(scip, cons, quadvarterm->var) );

   /* if we catch variable events, drop the events on the old variable */
   if( quadvarterm->eventdata != NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* get event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      eventhdlr = conshdlrdata->eventhdlr;

      /* drop bound change events of variable */
      SCIP_CALL( dropQuadVarEvents(scip, eventhdlr, cons, pos) );
   }
   else
   {
      eventhdlr = NULL;
   }

   /* compute constant and put into lhs/rhs */
   constant = quadvarterm->lincoef * offset + quadvarterm->sqrcoef * offset * offset;
   if( constant != 0.0 )
   {
      /* maintain constant part */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         consdata->lhs -= constant;
      if( !SCIPisInfinity(scip,  consdata->rhs) )
         consdata->rhs -= constant;
   }

   /* update linear and square coefficient */
   quadvarterm->lincoef *= coef;
   quadvarterm->lincoef += 2.0 * quadvarterm->sqrcoef * coef * offset;
   quadvarterm->sqrcoef *= coef * coef;

   /* update bilinear terms */
   for( i = 0; i < quadvarterm->nadjbilin; ++i )
   {
      bilinterm = &consdata->bilinterms[quadvarterm->adjbilin[i]];

      if( bilinterm->var1 == quadvarterm->var )
      {
         bilinterm->var1 = var;
         var2 = bilinterm->var2;
      }
      else
      {
         assert(bilinterm->var2 == quadvarterm->var);
         bilinterm->var2 = var;
         var2 = bilinterm->var1;
      }

      if( var == var2 )
      {
         /* looks like we actually have a square term here */
         quadvarterm->lincoef += bilinterm->coef * offset;
         quadvarterm->sqrcoef += bilinterm->coef * coef;
         /* deleting bilinear terms is expensive, since it requires updating adjacency information
          * thus, for now we just set the coefficient to 0.0 and clear in later when the bilinear terms are merged */
         bilinterm->coef = 0.0;
         continue;
      }

      /* swap var1 and var2 if they are in wrong order */
      if( SCIPvarCompare(bilinterm->var1, bilinterm->var2) > 0 )
      {
         SCIP_VAR* tmp;
         tmp = bilinterm->var1;
         bilinterm->var1 = bilinterm->var2;
         bilinterm->var2 = tmp;
      }
      assert(SCIPvarCompare(bilinterm->var1, bilinterm->var2) == -1);

      if( offset != 0.0 )
      {
         /* need to find var2 and add offset*bilinterm->coef to linear coefficient */
         int var2pos;

         var2pos = 0;
         while( consdata->quadvarterms[var2pos].var != var2 )
         {
            ++var2pos;
            assert(var2pos < consdata->nquadvars);
         }

         consdata->quadvarterms[var2pos].lincoef += bilinterm->coef * offset;
      }

      bilinterm->coef *= coef;
   }

   /* release old variable */
   SCIP_CALL( SCIPreleaseVar(scip, &quadvarterm->var) );

   /* set new variable */
   quadvarterm->var = var;

   /* capture new variable */
   SCIP_CALL( SCIPcaptureVar(scip, quadvarterm->var) );

   /* catch variable events, if we do so */
   if( eventhdlr != NULL )
   {
      assert(SCIPconsIsEnabled(cons));
      
      /* catch bound change events of variable */
      SCIP_CALL( catchQuadVarEvents(scip, eventhdlr, cons, pos) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockQuadraticVariable(scip, cons, var) );

   consdata->isremovedfixings = consdata->isremovedfixings && SCIPvarIsActive(var)
      && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   consdata->quadvarssorted = (consdata->nquadvars == 1);
   consdata->quadvarsmerged = FALSE;
   consdata->bilinsorted &= (quadvarterm->nadjbilin == 0);  /*lint !e514*/
   consdata->bilinmerged &= (quadvarterm->nadjbilin == 0);  /*lint !e514*/

   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = FALSE;
   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** adds a bilinear term to quadratic constraint */
static
SCIP_RETCODE addBilinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   var1pos,            /**< position of first variable in quadratic variables array */
   int                   var2pos,            /**< position of second variable in quadratic variables array */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;

   assert(scip != NULL);
   assert(cons != NULL);

   if( var1pos == var2pos )
   {
      SCIPerrorMessage("tried to add bilinear term where both variables are the same\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if the bilinear terms are sorted */
   assert(consdataCheckBilinTermsSort(consdata));

   assert(var1pos >= 0);
   assert(var1pos < consdata->nquadvars);
   assert(var2pos >= 0);
   assert(var2pos < consdata->nquadvars);

   SCIP_CALL( consdataEnsureBilinSize(scip, consdata, consdata->nbilinterms + 1) );

   bilinterm = &consdata->bilinterms[consdata->nbilinterms];
   if( SCIPvarCompare(consdata->quadvarterms[var1pos].var, consdata->quadvarterms[var2pos].var) < 0 )
   {
      bilinterm->var1 = consdata->quadvarterms[var1pos].var;
      bilinterm->var2 = consdata->quadvarterms[var2pos].var;
   }
   else
   {
      bilinterm->var1 = consdata->quadvarterms[var2pos].var;
      bilinterm->var2 = consdata->quadvarterms[var1pos].var;
   }
   bilinterm->coef = coef;

   if( bilinterm->var1 == bilinterm->var2 )
   {
      SCIPerrorMessage("tried to add bilinear term where both variables are the same, but appear at different positions in quadvarterms array\n");
      return SCIP_INVALIDDATA;
   }
   assert(SCIPvarCompare(bilinterm->var1, bilinterm->var2) == -1);

   SCIP_CALL( consdataEnsureAdjBilinSize(scip, &consdata->quadvarterms[var1pos], consdata->quadvarterms[var1pos].nadjbilin + 1) );
   SCIP_CALL( consdataEnsureAdjBilinSize(scip, &consdata->quadvarterms[var2pos], consdata->quadvarterms[var2pos].nadjbilin + 1) );

   consdata->quadvarterms[var1pos].adjbilin[consdata->quadvarterms[var1pos].nadjbilin] = consdata->nbilinterms;
   consdata->quadvarterms[var2pos].adjbilin[consdata->quadvarterms[var2pos].nadjbilin] = consdata->nbilinterms;
   ++consdata->quadvarterms[var1pos].nadjbilin;
   ++consdata->quadvarterms[var2pos].nadjbilin;

   ++consdata->nbilinterms;

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved  = FALSE;
   if( consdata->nbilinterms == 1 )
   {
      consdata->bilinsorted = TRUE;

      /* we have to take care of the bilinear term in mergeAndCleanBilinearTerms() if the coefficient is zero */
      consdata->bilinmerged = !SCIPisZero(scip, consdata->bilinterms[0].coef);
   }
   else
   {
      consdata->bilinsorted = consdata->bilinsorted
         && (bilinTermComp(consdata, consdata->nbilinterms-2, consdata->nbilinterms-1) <= 0);
      consdata->bilinmerged = FALSE;
   }

   consdata->iscurvchecked = FALSE;

   /* check if the bilinear terms are sorted */
   assert(consdataCheckBilinTermsSort(consdata));

   return SCIP_OKAY;
}

/** removes a set of bilinear terms and updates adjacency information in quad var terms
 *
 * Note: this function sorts the given array termposs.
 */
static
SCIP_RETCODE removeBilinearTermsPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int                   nterms,             /**< number of terms to delete */
   int*                  termposs            /**< indices of terms to delete */
   )
{
   SCIP_CONSDATA* consdata;
   int* newpos;
   int i;
   int j;
   int offset;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nterms == 0 || termposs != NULL);

   if( nterms == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPsortInt(termposs, nterms);

   SCIP_CALL( SCIPallocBufferArray(scip, &newpos, consdata->nbilinterms) );

   i = 0;
   offset = 0;
   for( j = 0; j < consdata->nbilinterms; ++j )
   {
      /* if j'th term is deleted, increase offset and continue */
      if( i < nterms && j == termposs[i] )
      {
         ++offset;
         ++i;
         newpos[j] = -1;
         continue;
      }

      /* otherwise, move it forward and remember new position */
      if( offset > 0 )
         consdata->bilinterms[j-offset] = consdata->bilinterms[j];
      newpos[j] = j - offset;
   }
   assert(offset == nterms);

   /* update adjacency and activity information in quad var terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      offset = 0;
      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         assert(consdata->quadvarterms[i].adjbilin[j] < consdata->nbilinterms);
         if( newpos[consdata->quadvarterms[i].adjbilin[j]] == -1 )
         {
            /* corresponding bilinear term was deleted, thus increase offset */
            ++offset;
         }
         else
         {
            /* update index of j'th bilinear term and store at position j-offset */
            consdata->quadvarterms[i].adjbilin[j-offset] = newpos[consdata->quadvarterms[i].adjbilin[j]];
         }
      }
      consdata->quadvarterms[i].nadjbilin -= offset;
      /* some bilinear term was removed, so invalidate activity bounds */
   }

   consdata->nbilinterms -= nterms;

   SCIPfreeBufferArray(scip, &newpos);

   /* some quad vars may be linear now */
   consdata->quadvarsmerged = FALSE;

   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = FALSE;
   consdata->iscurvchecked = FALSE;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   return SCIP_OKAY;
}

/** merges quad var terms that correspond to the same variable and does additional cleanup
 *
 *  If a quadratic variable terms is actually linear, makes a linear term out of it
 *  also replaces squares of binary variables by the binary variables, i.e., adds sqrcoef to lincoef.
 */
static
SCIP_RETCODE mergeAndCleanQuadVarTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_CONSDATA* consdata;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( consdata->quadvarsmerged )
      return SCIP_OKAY;

   if( consdata->nquadvars == 0 )
   {
      consdata->quadvarsmerged = TRUE;
      return SCIP_OKAY;
   }

   i = 0;
   while( i < consdata->nquadvars )
   {
      /* make sure quad var terms are sorted (do this in every round, since we may move variables around) */
      SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

      quadvarterm = &consdata->quadvarterms[i];

      for( j = i+1; j < consdata->nquadvars && consdata->quadvarterms[j].var == quadvarterm->var; ++j )
      {
         /* add quad var term j to current term i */
         quadvarterm->lincoef += consdata->quadvarterms[j].lincoef;
         quadvarterm->sqrcoef += consdata->quadvarterms[j].sqrcoef;
         if( consdata->quadvarterms[j].nadjbilin > 0 )
         {
            /* move adjacency information from j to i */
            SCIP_CALL( consdataEnsureAdjBilinSize(scip, quadvarterm, quadvarterm->nadjbilin + consdata->quadvarterms[j].nadjbilin) );
            BMScopyMemoryArray(&quadvarterm->adjbilin[quadvarterm->nadjbilin], consdata->quadvarterms[j].adjbilin, consdata->quadvarterms[j].nadjbilin);  /*lint !e866*/
            quadvarterm->nadjbilin += consdata->quadvarterms[j].nadjbilin;
            consdata->quadvarterms[j].nadjbilin = 0;
         }
         consdata->quadvarterms[j].lincoef = 0.0;
         consdata->quadvarterms[j].sqrcoef = 0.0;
         /* mark that activity information in quadvarterm is not up to date anymore */
      }

      /* remove quad var terms i+1..j-1 backwards */
      for( j = j-1; j > i; --j )
      {
         SCIP_CALL( delQuadVarTermPos(scip, cons, j) );
      }

      /* for binary variables, x^2 = x
       * however, we may destroy convexity of a quadratic term that involves also bilinear terms
       * thus, we do this step only if the variable does not appear in any bilinear term */
      if( quadvarterm->sqrcoef != 0.0 && SCIPvarIsBinary(quadvarterm->var) && quadvarterm->nadjbilin == 0 )
      {
         SCIPdebugMsg(scip, "replace square of binary variable by itself: <%s>^2 --> <%s>\n", SCIPvarGetName(quadvarterm->var), SCIPvarGetName(quadvarterm->var));
         quadvarterm->lincoef += quadvarterm->sqrcoef;
         quadvarterm->sqrcoef = 0.0;

         /* invalidate nonlinear row */
         if( consdata->nlrow != NULL )
         {
            SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
         }
      }

      /* if its 0.0 or linear, get rid of it */
      if( SCIPisZero(scip, quadvarterm->sqrcoef) && quadvarterm->nadjbilin == 0 )
      {
         if( !SCIPisZero(scip, quadvarterm->lincoef) )
         {
            /* seem to be a linear term now, thus add as linear term */
            SCIP_CALL( addLinearCoef(scip, cons, quadvarterm->var, quadvarterm->lincoef) );
         }
         /* remove term at pos i */
         SCIP_CALL( delQuadVarTermPos(scip, cons, i) );
      }
      else
      {
         ++i;
      }
   }

   consdata->quadvarsmerged = TRUE;
   SCIPintervalSetEmpty(&consdata->quadactivitybounds);

   return SCIP_OKAY;
}

/** merges entries with same linear variable into one entry and cleans up entries with coefficient 0.0 */
static
SCIP_RETCODE mergeAndCleanLinearVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real newcoef;
   int i;
   int j;
   int qvarpos;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( consdata->linvarsmerged )
      return SCIP_OKAY;

   if( consdata->nlinvars == 0 )
   {
      consdata->linvarsmerged = TRUE;
      return SCIP_OKAY;
   }

   i = 0;
   while( i < consdata->nlinvars )
   {
      /* make sure linear variables are sorted (do this in every round, since we may move variables around) */
      consdataSortLinearVars(consdata);

      /* sum up coefficients that correspond to variable i */
      newcoef = consdata->lincoefs[i];
      for( j = i+1; j < consdata->nlinvars && consdata->linvars[i] == consdata->linvars[j]; ++j )
         newcoef += consdata->lincoefs[j];
      /* delete the additional variables in backward order */ 
      for( j = j-1; j > i; --j )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, j) );
      }

      /* check if there is already a quadratic variable term with this variable */
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->linvars[i], &qvarpos) );
      if( qvarpos >= 0)
      {
         /* add newcoef to linear coefficient of quadratic variable and mark linear variable as to delete */
         assert(qvarpos < consdata->nquadvars);
         assert(consdata->quadvarterms[qvarpos].var == consdata->linvars[i]);
         consdata->quadvarterms[qvarpos].lincoef += newcoef;
         newcoef = 0.0;
         SCIPintervalSetEmpty(&consdata->quadactivitybounds);
      }

      /* delete also entry at position i, if it became zero (or was zero before) */
      if( SCIPisZero(scip, newcoef) )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, i) );
      }
      else
      {
         SCIP_CALL( chgLinearCoefPos(scip, cons, i, newcoef) );
         ++i;
      }
   }

   consdata->linvarsmerged = TRUE;

   return SCIP_OKAY;
}

/** merges bilinear terms with same variables into a single term, removes bilinear terms with coefficient 0.0 */
static
SCIP_RETCODE mergeAndCleanBilinearTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   int i;
   int j;
   int* todelete;
   int ntodelete;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   /* check if the bilinear terms are sorted */
   assert(consdataCheckBilinTermsSort(consdata));

   if( consdata->bilinmerged )
      return SCIP_OKAY;

   if( consdata->nbilinterms == 0 )
   {
      consdata->bilinmerged = TRUE;
      return SCIP_OKAY;
   }

   /* alloc memory for array of terms that need to be deleted finally */
   ntodelete = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &todelete, consdata->nbilinterms) );

   /* make sure bilinear terms are sorted */
   SCIP_CALL( consdataSortBilinTerms(scip, consdata) );

   i = 0;
   while( i < consdata->nbilinterms )
   {
      bilinterm = &consdata->bilinterms[i];

      /* sum up coefficients that correspond to same variables as term i */
      for( j = i+1; j < consdata->nbilinterms && bilinterm->var1 == consdata->bilinterms[j].var1 && bilinterm->var2 == consdata->bilinterms[j].var2; ++j )
      {
         bilinterm->coef += consdata->bilinterms[j].coef;
         todelete[ntodelete++] = j;
      }

      /* delete also entry at position i, if it became zero (or was zero before) */
      if( SCIPisZero(scip, bilinterm->coef) )
      {
         todelete[ntodelete++] = i;
      }

      /* continue with term after the current series */
      i = j;
   }

   /* delete bilinear terms */
   SCIP_CALL( removeBilinearTermsPos(scip, cons, ntodelete, todelete) );

   SCIPfreeBufferArray(scip, &todelete);

   consdata->bilinmerged = TRUE;

   /* check if the bilinear terms are sorted */
   assert(consdataCheckBilinTermsSort(consdata));

   return SCIP_OKAY;
}

/** removes fixes (or aggregated) variables from a quadratic constraint */
static
SCIP_RETCODE removeFixedVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real bilincoef;
   SCIP_Real coef;
   SCIP_Real offset;
   SCIP_VAR* var;
   SCIP_VAR* var2;
   int var2pos;
   int i;
   int j;
   int k;

   SCIP_Bool have_change;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   have_change = FALSE;
   i = 0;
   while( i < consdata->nlinvars )
   {
      var = consdata->linvars[i];

      if( SCIPvarIsActive(var) && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
      {
         ++i;
         continue;
      }

      have_change = TRUE;

      coef = consdata->lincoefs[i];
      offset = 0.0;

      if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
      {
         offset = coef * (SCIPvarGetLbGlobal(var) + SCIPvarGetUbGlobal(var)) / 2.0;
         coef = 0.0;
      }
      else
      {
         SCIP_CALL( SCIPgetProbvarSum(scip, &var, &coef, &offset) );
      }

      SCIPdebugMsg(scip, "  linear term %g*<%s> is replaced by %g * <%s> + %g\n", consdata->lincoefs[i], SCIPvarGetName(consdata->linvars[i]),
         coef, SCIPvarGetName(var), offset);

      /* delete previous variable (this will move another variable to position i) */
      SCIP_CALL( delLinearCoefPos(scip, cons, i) );

      /* put constant part into bounds */
      if( offset != 0.0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhs -= offset;
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhs -= offset;
      }

      /* nothing left to do if variable had been fixed */
      if( coef == 0.0 )
         continue;

      /* if GetProbvar gave a linear variable, just add it
       * if it's a multilinear variable, add it's disaggregated variables */
      if( SCIPvarIsActive(var) )
      {
         SCIP_CALL( addLinearCoef(scip, cons, var, coef) );
      }
      else
      {
         int        naggrs;
         SCIP_VAR** aggrvars;
         SCIP_Real* aggrscalars;
         SCIP_Real  aggrconstant;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

         naggrs = SCIPvarGetMultaggrNVars(var);
         aggrvars = SCIPvarGetMultaggrVars(var);
         aggrscalars = SCIPvarGetMultaggrScalars(var);
         aggrconstant = SCIPvarGetMultaggrConstant(var);

         SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars + naggrs) );

         for( j = 0; j < naggrs; ++j )
         {
            SCIP_CALL( addLinearCoef(scip, cons, aggrvars[j], coef * aggrscalars[j]) );
         }

         if( aggrconstant != 0.0 )
         {
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= coef * aggrconstant;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= coef * aggrconstant;
         }
      }
   }

   i = 0;
   while( i < consdata->nquadvars )
   {
      var = consdata->quadvarterms[i].var;

      if( SCIPvarIsActive(var) && !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
      {
         ++i;
         continue;
      }

      have_change = TRUE;

      coef   = 1.0;
      offset = 0.0;

      if( !SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
      {
         SCIP_CALL( SCIPgetProbvarSum(scip, &var, &coef, &offset) );
      }
      else
      {
         coef = 0.0;
         offset = (SCIPvarGetLbGlobal(var) + SCIPvarGetUbGlobal(var)) / 2.0;
      }

      SCIPdebugMsg(scip, "  quadratic variable <%s> with status %d is replaced by %g * <%s> + %g\n", SCIPvarGetName(consdata->quadvarterms[i].var),
         SCIPvarGetStatus(consdata->quadvarterms[i].var), coef, SCIPvarGetName(var), offset);

      /* handle fixed variable */
      if( coef == 0.0 )
      {
         /* if not fixed to 0.0, add to linear coefs of vars in bilinear terms, and deal with linear and square term as constant */ 
         if( offset != 0.0 )
         {
            for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
            {
               bilinterm = &consdata->bilinterms[consdata->quadvarterms[i].adjbilin[j]];

               var2 = bilinterm->var1 == consdata->quadvarterms[i].var ? bilinterm->var2 : bilinterm->var1;
               assert(var2 != consdata->quadvarterms[i].var);

               var2pos = 0;
               while( consdata->quadvarterms[var2pos].var != var2 )
               {
                  ++var2pos;
                  assert(var2pos < consdata->nquadvars);
               }
               consdata->quadvarterms[var2pos].lincoef += bilinterm->coef * offset;
            }

            offset = consdata->quadvarterms[i].lincoef * offset + consdata->quadvarterms[i].sqrcoef * offset * offset;
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= offset;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= offset;
         }

         /* remove bilinear terms */
         SCIP_CALL( removeBilinearTermsPos(scip, cons, consdata->quadvarterms[i].nadjbilin, consdata->quadvarterms[i].adjbilin) );

         /* delete quad. var term i */
         SCIP_CALL( delQuadVarTermPos(scip, cons, i) );

         continue;
      }

      assert(var != NULL);

      /* if GetProbvar gave an active variable, replace the quad var term so that it uses the new variable */
      if( SCIPvarIsActive(var) )
      {
         /* replace x by coef*y+offset */
         SCIP_CALL( replaceQuadVarTermPos(scip, cons, i, var, coef, offset) );

         continue;
      }
      else
      {
         /* if GetProbVar gave a multi-aggregated variable, add new quad var terms and new bilinear terms
          * x is replaced by coef * (sum_i a_ix_i + b) + offset
          * lcoef * x + scoef * x^2 + bcoef * x * y ->
          *   (b*coef + offset) * (lcoef + (b*coef + offset) * scoef)
          * + sum_i a_i*coef * (lcoef + 2 (b*coef + offset) * scoef) x_i
          * + sum_i (a_i*coef)^2 * scoef * x_i^2
          * + 2 sum_{i,j, i<j} (a_i a_j coef^2 scoef) x_i x_j
          * + bcoef * (b*coef + offset + coef * sum_i a_ix_i) y 
          */
         int        naggrs;
         SCIP_VAR** aggrvars;     /* x_i */
         SCIP_Real* aggrscalars;  /* a_i */
         SCIP_Real  aggrconstant; /* b */
         int nquadtermsold;

         SCIP_Real lcoef;
         SCIP_Real scoef;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

         naggrs = SCIPvarGetMultaggrNVars(var);
         aggrvars = SCIPvarGetMultaggrVars(var);
         aggrscalars = SCIPvarGetMultaggrScalars(var);
         aggrconstant = SCIPvarGetMultaggrConstant(var);

         lcoef = consdata->quadvarterms[i].lincoef;
         scoef = consdata->quadvarterms[i].sqrcoef;

         nquadtermsold = consdata->nquadvars;

         SCIP_CALL( consdataEnsureQuadVarTermsSize(scip, consdata, consdata->nquadvars + naggrs) );

         /* take care of constant part */
         if( aggrconstant != 0.0 || offset != 0.0 )
         {
            SCIP_Real constant;
            constant = (aggrconstant * coef + offset) * (lcoef + (aggrconstant * coef + offset) * scoef);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= constant;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= constant;
         }

         /* add x_i's with linear and square coefficients */
         for( j = 0; j < naggrs; ++j )
         {
            SCIP_CALL( addQuadVarTerm(scip, cons, aggrvars[j],
                  coef * aggrscalars[j] * (lcoef + 2.0 * scoef * (coef * aggrconstant + offset)),
                  coef * coef * aggrscalars[j] * aggrscalars[j] * scoef) );
         }

         /* ensure space for bilinear terms */
         SCIP_CALL( consdataEnsureBilinSize(scip, consdata, consdata->nquadvars + (scoef != 0.0 ? (naggrs * (naggrs-1))/2 : 0) + consdata->quadvarterms[j].nadjbilin * naggrs) );

         /* add x_j*x_k's */
         if( scoef != 0.0 )
         {
            for( j = 0; j < naggrs; ++j )
               for( k = 0; k < j; ++k )
               {
                  assert(aggrvars[j] != aggrvars[k]);
                  SCIP_CALL( addBilinearTerm(scip, cons, nquadtermsold + j, nquadtermsold + k, 
                        2.0 * aggrscalars[j] * aggrscalars[k] * coef * coef * scoef) );
               }
         }

         /* add x_i*y's */
         for( k = 0; k < consdata->quadvarterms[i].nadjbilin; ++k )
         {
            bilinterm = &consdata->bilinterms[consdata->quadvarterms[i].adjbilin[k]];
            bilincoef = bilinterm->coef;   /* copy coef, as bilinterm pointer may become invalid by realloc in addBilinearTerm() below */
            var2 = (bilinterm->var1 == consdata->quadvarterms[i].var) ? bilinterm->var2 : bilinterm->var1;
            assert(var2 != consdata->quadvarterms[i].var);

            /* this is not efficient, but we cannot sort the quadratic terms here, since we currently iterate over them */
            var2pos = 0;
            while( consdata->quadvarterms[var2pos].var != var2 )
            {
               ++var2pos;
               assert(var2pos < consdata->nquadvars);
            }

            for( j = 0; j < naggrs; ++j )
            {
               if( aggrvars[j] == var2 )
               { /* x_i == y, so we have a square term here */
                  consdata->quadvarterms[var2pos].sqrcoef += bilincoef * coef * aggrscalars[j];
               }
               else
               { /* x_i != y, so we need to add a bilinear term here */
                  SCIP_CALL( addBilinearTerm(scip, cons, nquadtermsold + j, var2pos, bilincoef * coef * aggrscalars[j]) );
               }
            }

            consdata->quadvarterms[var2pos].lincoef += bilincoef * (aggrconstant * coef + offset);
         }

         /* remove bilinear terms */
         SCIP_CALL( removeBilinearTermsPos(scip, cons, consdata->quadvarterms[i].nadjbilin, consdata->quadvarterms[i].adjbilin) );

         /* delete quad. var term i */
         SCIP_CALL( delQuadVarTermPos(scip, cons, i) );
      }
   }

   consdata->isremovedfixings = TRUE;

   SCIPdebugMsg(scip, "removed fixations from <%s>\n  -> ", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

#ifndef NDEBUG
   for( i = 0; i < consdata->nlinvars; ++i )
      assert(SCIPvarIsActive(consdata->linvars[i]));

   for( i = 0; i < consdata->nquadvars; ++i )
      assert(SCIPvarIsActive(consdata->quadvarterms[i].var));
#endif

   if( !have_change )
      return SCIP_OKAY;

   /* some quadratic variable may have been replaced by an already existing linear variable
    * in this case, we want the linear variable to be removed, which happens in mergeAndCleanLinearVars
    */ 
   consdata->linvarsmerged = FALSE;

   SCIP_CALL( mergeAndCleanBilinearTerms(scip, cons) );
   SCIP_CALL( mergeAndCleanQuadVarTerms(scip, cons) );
   SCIP_CALL( mergeAndCleanLinearVars(scip, cons) );

#ifndef NDEBUG
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      assert(consdata->bilinterms[i].var1 != consdata->bilinterms[i].var2);
      assert(consdata->bilinterms[i].coef != 0.0);
      assert(SCIPvarCompare(consdata->bilinterms[i].var1, consdata->bilinterms[i].var2) < 0);
   }
#endif

   return SCIP_OKAY;
}

/** create a nonlinear row representation of the constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< quadratic constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int        nquadvars;     /* number of variables in quadratic terms */
   SCIP_VAR** quadvars;      /* variables in quadratic terms */
   int        nquadelems;    /* number of quadratic elements (square and bilinear terms) */
   SCIP_QUADELEM* quadelems; /* quadratic elements (square and bilinear terms) */
   int        nquadlinterms; /* number of linear terms using variables that are in quadratic terms */
   SCIP_VAR** quadlinvars;   /* variables of linear terms using variables that are in quadratic terms */
   SCIP_Real* quadlincoefs;  /* coefficients of linear terms using variables that are in quadratic terms */
   int i;
   int idx1;
   int idx2;
   int lincnt;
   int elcnt;
   SCIP_VAR* lastvar;
   int lastvaridx;
   SCIP_EXPRCURV curvature;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   nquadvars = consdata->nquadvars;
   nquadelems = consdata->nbilinterms;
   nquadlinterms = 0;
   for( i = 0; i < nquadvars; ++i )
   {
      if( consdata->quadvarterms[i].sqrcoef != 0.0 )
         ++nquadelems;
      if( !SCIPisZero(scip, consdata->quadvarterms[i].lincoef) )
         ++nquadlinterms;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &quadvars,  nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nquadelems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadlinvars,  nquadlinterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadlincoefs, nquadlinterms) );

   lincnt = 0;
   elcnt = 0;
   for( i = 0; i < nquadvars; ++i )
   {
      quadvars[i] = consdata->quadvarterms[i].var;

      if( consdata->quadvarterms[i].sqrcoef != 0.0 )
      {
         assert(elcnt < nquadelems);
         quadelems[elcnt].idx1 = i;
         quadelems[elcnt].idx2 = i;
         quadelems[elcnt].coef = consdata->quadvarterms[i].sqrcoef;
         ++elcnt;
      }

      if( !SCIPisZero(scip, consdata->quadvarterms[i].lincoef) )
      {
         assert(lincnt < nquadlinterms);
         quadlinvars [lincnt] = consdata->quadvarterms[i].var;
         quadlincoefs[lincnt] = consdata->quadvarterms[i].lincoef;
         ++lincnt;
      }
   }
   assert(lincnt == nquadlinterms);

   /* bilinear terms are sorted first by first variable, then by second variable
    * thus, it makes sense to remember the index of the previous first variable for the case a series of bilinear terms with the same first var appears */
   lastvar = NULL;
   lastvaridx = -1;
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      if( lastvar == consdata->bilinterms[i].var1 )
      {
         assert(lastvaridx >= 0);
         assert(consdata->quadvarterms[lastvaridx].var == consdata->bilinterms[i].var1);
      }
      else
      {
         lastvar = consdata->bilinterms[i].var1;
         SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, lastvar, &lastvaridx) );
      }
      idx1 = lastvaridx;

      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->bilinterms[i].var2, &idx2) );

      assert(elcnt < nquadelems);
      quadelems[elcnt].idx1 = MIN(idx1, idx2);
      quadelems[elcnt].idx2 = MAX(idx1, idx2);
      quadelems[elcnt].coef = consdata->bilinterms[i].coef;
      ++elcnt;
   }
   assert(elcnt == nquadelems);

   /* set curvature for the nonlinear row */
   if( consdata->isconcave && consdata->isconvex )
   {
      assert(consdata->nbilinterms == 0 && consdata->nquadvars == 0);
      curvature = SCIP_EXPRCURV_LINEAR;
   }
   else if( consdata->isconcave )
      curvature = SCIP_EXPRCURV_CONCAVE;
   else if( consdata->isconvex )
      curvature = SCIP_EXPRCURV_CONVEX;
   else
      curvature = SCIP_EXPRCURV_UNKNOWN;

   SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
         consdata->nlinvars, consdata->linvars, consdata->lincoefs,
         nquadvars, quadvars, nquadelems, quadelems,
         NULL, consdata->lhs, consdata->rhs,
         curvature) );

   SCIP_CALL( SCIPaddLinearCoefsToNlRow(scip, consdata->nlrow, nquadlinterms, quadlinvars, quadlincoefs) );

   SCIPfreeBufferArray(scip, &quadlincoefs);
   SCIPfreeBufferArray(scip, &quadlinvars);
   SCIPfreeBufferArray(scip, &quadelems);
   SCIPfreeBufferArray(scip, &quadvars);

   return SCIP_OKAY;
}

/** solve constraint as presolving */
static
SCIP_RETCODE presolveSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_RESULT*          result,             /**< to store result of solve: cutoff, success, or do-not-find */
   SCIP_Bool*            redundant,          /**< to store whether constraint is redundant now (should be deleted) */
   int*                  naggrvars           /**< counter on number of variable aggregations */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);
   assert(redundant != NULL);

   *result = SCIP_DIDNOTFIND;
   *redundant = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* if constraint is an equality with two variables, at least one of them binary,
    * and linear after fixing the binary, then we can aggregate the variables */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) && consdata->nlinvars == 0 && consdata->nquadvars == 2 &&
      ((SCIPvarIsBinary(consdata->quadvarterms[0].var) && consdata->quadvarterms[1].sqrcoef == 0.0) ||
       (SCIPvarIsBinary(consdata->quadvarterms[1].var) && consdata->quadvarterms[0].sqrcoef == 0.0)) )
   {
      SCIP_Bool infeasible;
      SCIP_Bool aggregated;
      SCIP_Real a;
      SCIP_Real b;
      SCIP_Real c;
      SCIP_VAR* x;
      SCIP_VAR* y;
      int binvaridx;

      /* constraint is a*(x+x^2) + b*y + c*x*y = rhs, with x binary variable
       * x = 0 -> b*y == rhs
       * x = 1 -> (b+c)*y == rhs - a
       *
       * if b != 0 and b+c != 0, then y = (rhs-a)/(b+c) * x + rhs/b * (1-x) = ((rhs-a)/(b+c) - rhs/b) * x + rhs/b
       */

      binvaridx = (SCIPvarIsBinary(consdata->quadvarterms[0].var) && consdata->quadvarterms[1].sqrcoef == 0.0) ? 0 : 1;

      x = consdata->quadvarterms[binvaridx].var;
      a = consdata->quadvarterms[binvaridx].sqrcoef + consdata->quadvarterms[binvaridx].lincoef;

      y = consdata->quadvarterms[1-binvaridx].var;
      b = consdata->quadvarterms[1-binvaridx].lincoef;

      assert(consdata->nbilinterms <= 1);  /* should actually be 1, since constraint is otherwise linear */
      c = (consdata->nbilinterms == 1) ? consdata->bilinterms[0].coef : 0.0;

      if( !SCIPisZero(scip, b) && !SCIPisZero(scip, b+c) )
      {
         SCIPdebugMsg(scip, "<%s> = 0 -> %g*<%s> = %g  and  <%s> = 1 -> %g*<%s> = %g\n", SCIPvarGetName(x), b, SCIPvarGetName(y), consdata->rhs,
            SCIPvarGetName(x), b+c, SCIPvarGetName(y), consdata->rhs - a);
         SCIPdebugMsg(scip, "=> attempt aggregation <%s> = %g*<%s> + %g\n", SCIPvarGetName(y), (consdata->rhs-a)/(b+c) - consdata->rhs/b,
            SCIPvarGetName(x), consdata->rhs/b);

         SCIP_CALL( SCIPaggregateVars(scip, x, y, (consdata->rhs-a)/(b+c) - consdata->rhs/b, -1.0, -consdata->rhs/b, &infeasible, redundant, &aggregated) );
         if( infeasible )
            *result = SCIP_CUTOFF;
         else if( *redundant || aggregated )
         {
            /* aggregated (or were already aggregated), so constraint is now redundant */
            *result = SCIP_SUCCESS;
            *redundant = TRUE;

            if( aggregated )
               ++*naggrvars;
         }
      }

      /* @todo if b is 0 or b+c is 0, or lhs != rhs, then could replace by varbound constraint */
   }

   return SCIP_OKAY;
}


/** reformulates products of binary variables as AND constraint
 *
 *  For a product x*y, with x and y binary variables, the product is replaced by a new auxiliary variable z and the constraint z = {x and y} is added.
 */
static
SCIP_RETCODE presolveTryAddAND(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  naddconss           /**< buffer where to add the number of AND constraints added */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   char               name[SCIP_MAXSTRLEN];
   SCIP_VAR*          vars[2];
   SCIP_VAR*          auxvar;
   SCIP_CONS*         andcons;
   int                i;
   int                ntodelete;
   int*               todelete;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(naddconss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if no binary variables, then we will find nothing to reformulate here
    * (note that this does not count in integer variables with {0,1} bounds...)
    */
   if( SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;

   /* if user does not like AND very much, then return */
   if( conshdlrdata->empathy4and < 2 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nbilinterms == 0 )
      return SCIP_OKAY;

   /* get array to store indices of bilinear terms that shall be deleted */
   SCIP_CALL( SCIPallocBufferArray(scip, &todelete, consdata->nbilinterms) );
   ntodelete = 0;

   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      vars[0] = consdata->bilinterms[i].var1;
      if( !SCIPvarIsBinary(vars[0]) )
         continue;

      vars[1] = consdata->bilinterms[i].var2;
      if( !SCIPvarIsBinary(vars[1]) )
         continue;

      /* create auxiliary variable */ 
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s_%s_%s", SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]), SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, 
            SCIPvarIsInitial(vars[0]) || SCIPvarIsInitial(vars[1]), SCIPvarIsRemovable(vars[0]) && SCIPvarIsRemovable(vars[1]), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );
#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real var0val;
         SCIP_Real var1val;
         SCIP_CALL( SCIPdebugGetSolVal(scip, vars[0], &var0val) );
         SCIP_CALL( SCIPdebugGetSolVal(scip, vars[1], &var1val) );
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, var0val * var1val) );
      }
#endif

      /* create AND-constraint auxvar = x and y, need to be enforced as not redundant */
      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%sAND%s", SCIPvarGetName(vars[0]), SCIPvarGetName(vars[1]));
      SCIP_CALL( SCIPcreateConsAnd(scip, &andcons, name, auxvar, 2, vars,
            SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
            SCIPconsIsSeparated(cons), TRUE, TRUE,
            SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, andcons) );
      SCIPdebugMsg(scip, "added AND constraint: ");
      SCIPdebugPrintCons(scip, andcons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &andcons) );
      ++*naddconss;

      /* add bilincoef * auxvar to linear terms */
      SCIP_CALL( addLinearCoef(scip, cons, auxvar, consdata->bilinterms[i].coef) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

      /* remember that we have to delete this bilinear term */
      assert(ntodelete < consdata->nbilinterms);
      todelete[ntodelete++] = i;
   }

   /* remove bilinear terms that have been replaced */
   SCIP_CALL( removeBilinearTermsPos(scip, cons, ntodelete, todelete) );
   SCIPfreeBufferArray(scip, &todelete);

   return SCIP_OKAY;
}

/** gets bounds of variable y if x takes a certain value; checks whether x = xval has implications on y */
static
SCIP_RETCODE getImpliedBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< variable which implications to check */
   SCIP_Bool             xval,               /**< value of x to check for (TRUE for 1, FALSE for 0) */
   SCIP_VAR*             y,                  /**< variable to check if bounds can be reduced */
   SCIP_INTERVAL*        resultant           /**< buffer to store bounds on y */
   )
{
   SCIP_VAR** implvars;
   SCIP_BOUNDTYPE* impltypes;
   SCIP_Real* implbounds;
   int nimpls;
   int pos;

   assert(scip != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(resultant != NULL);

   SCIPintervalSetBounds(resultant, MIN(SCIPvarGetLbGlobal(y), SCIPvarGetUbGlobal(y)), MAX(SCIPvarGetLbGlobal(y), SCIPvarGetUbGlobal(y)));  /*lint !e666 */

   if( !SCIPvarIsBinary(x) || !SCIPvarIsActive(x) )
      return SCIP_OKAY;

   /* check in cliques for binary to binary implications */
   if( SCIPvarIsBinary(y) )
   {
      resultant->inf = MAX(resultant->inf, MIN(resultant->sup, 0.0));
      resultant->sup = MIN(resultant->sup, MAX(resultant->inf, 1.0));

      if( SCIPhaveVarsCommonClique(scip, x, xval, y, TRUE, FALSE) )
      {
         resultant->sup = MIN(resultant->sup, MAX(resultant->inf, 0.0));
      }
      else if( SCIPhaveVarsCommonClique(scip, x, xval, y, FALSE, FALSE) )
      {
         resultant->inf = MAX(resultant->inf, MIN(resultant->sup, 1.0));
      }

      return SCIP_OKAY;
   }

   /* analyze implications for x = xval */
   nimpls = SCIPvarGetNImpls(x, xval);
   if( nimpls == 0 )
      return SCIP_OKAY;

   implvars   = SCIPvarGetImplVars  (x, xval);
   impltypes  = SCIPvarGetImplTypes (x, xval);
   implbounds = SCIPvarGetImplBounds(x, xval);

   assert(implvars != NULL);
   assert(impltypes != NULL);
   assert(implbounds != NULL);

   /* find implications */
   if( !SCIPsortedvecFindPtr((void**)implvars, SCIPvarComp, (void*)y, nimpls, &pos) )
      return SCIP_OKAY;

   /* if there are several implications on y, go to the first one */
   while( pos > 0 && implvars[pos-1] == y )
      --pos;

   /* update implied lower and upper bounds on y
    * but make sure that resultant will not be empty, due to tolerances
    */
   while( pos < nimpls && implvars[pos] == y )
   {
      if( impltypes[pos] == SCIP_BOUNDTYPE_LOWER )
         resultant->inf = MAX(resultant->inf, MIN(resultant->sup, implbounds[pos]));
      else
         resultant->sup = MIN(resultant->sup, MAX(resultant->inf, implbounds[pos]));
      ++pos;
   }

   assert(resultant->sup >= resultant->inf);

   return SCIP_OKAY;
}

/** Reformulates products of binary times bounded continuous variables as system of linear inequalities (plus auxiliary variable).
 * 
 *  For a product x*y, with y a binary variable and x a continous variable with finite bounds,
 *  an auxiliary variable z and the inequalities \f$ x^L y \leq z \leq x^U y \f$ and \f$ x - (1-y) x^U \leq z \leq x - (1-y) x^L \f$ are added.
 * 
 *  If x is a linear term consisting of more than one variable, it is split up in groups of linear terms of length at most maxnrvar.
 *  For each product of linear term of length at most maxnrvar with y, an auxiliary z and linear inequalities are added.
 * 
 *  If y is a binary variable, the AND constraint \f$ z = x \wedge y \f$ may be added instead of linear constraints.
 */
static
SCIP_RETCODE presolveTryAddLinearReform(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  naddconss           /**< buffer where to add the number of auxiliary constraints added */
   )
{  /*lint --e{666} */
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR**         xvars;
   SCIP_Real*         xcoef;
   SCIP_INTERVAL      xbndszero;
   SCIP_INTERVAL      xbndsone;
   SCIP_INTERVAL      act0;
   SCIP_INTERVAL      act1;
   int                nxvars;
   SCIP_VAR*          y;
   SCIP_VAR*          bvar;
   char               name[SCIP_MAXSTRLEN];
   int                nbilinterms;
   SCIP_VAR*          auxvar;
   SCIP_CONS*         auxcons;
   int                i;
   int                j;
   int                k;
   int                bilinidx;
   SCIP_Real          bilincoef;
   SCIP_Real          mincoef;
   SCIP_Real          maxcoef;
   int*               todelete;
   int                ntodelete;
   int                maxnrvar;
   SCIP_Bool          integral;
   SCIP_Longint       gcd;
   SCIP_Bool          auxvarinitial;
   SCIP_Bool          auxvarremovable;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(naddconss != NULL);

   /* if no binary variables, then we will find nothing to reformulate here
    * (note that this does not count in integer variables with {0,1} bounds...)
    */
   if( SCIPgetNBinVars(scip) == 0 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   maxnrvar = conshdlrdata->replacebinaryprodlength;
   if( maxnrvar == 0 )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   xvars = NULL;
   xcoef = NULL;
   todelete = NULL;
   gcd = 0;

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      y = consdata->quadvarterms[i].var;
      if( !SCIPvarIsBinary(y) )
         continue;

      nbilinterms = consdata->quadvarterms[i].nadjbilin;
      if( nbilinterms == 0 )
         continue;

      SCIP_CALL( SCIPreallocBufferArray(scip, &xvars, MIN(maxnrvar, nbilinterms)+2) ); /* add 2 for later use when creating linear constraints */
      SCIP_CALL( SCIPreallocBufferArray(scip, &xcoef, MIN(maxnrvar, nbilinterms)+2) );

      /* alloc array to store indices of bilinear terms that shall be deleted */
      SCIP_CALL( SCIPreallocBufferArray(scip, &todelete, nbilinterms) );
      ntodelete = 0;

      auxvarinitial = SCIPvarIsInitial(y);
      auxvarremovable = SCIPvarIsRemovable(y);

      /* setup a list of bounded variables x_i with coefficients a_i that are multiplied with binary y: y*(sum_i a_i*x_i)
       * and compute range of sum_i a_i*x_i for the cases y = 0 and y = 1
       * we may need several rounds if maxnrvar < nbilinterms
       */
      j = 0;
      do
      {
         nxvars = 0;
         SCIPintervalSet(&xbndszero, 0.0);
         SCIPintervalSet(&xbndsone, 0.0);

         mincoef = SCIPinfinity(scip);
         maxcoef = 0.0;
         integral = TRUE;

         /* collect at most maxnrvar variables for x term */
         for( ; j < nbilinterms && nxvars < maxnrvar; ++j )
         {
            bilinidx = consdata->quadvarterms[i].adjbilin[j];
            assert(bilinidx >= 0);
            assert(bilinidx < consdata->nbilinterms);

            bvar = consdata->bilinterms[bilinidx].var1;
            if( bvar == y )
               bvar = consdata->bilinterms[bilinidx].var2;
            assert(bvar != y);

            /* skip products with unbounded variables */
            if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(bvar)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(bvar)) )
            {
               SCIPdebugMsg(scip, "skip reform of <%s><%s> due to unbounded second variable [%g,%g]\n",
                  SCIPvarGetName(y), SCIPvarGetName(bvar), SCIPvarGetLbGlobal(bvar), SCIPvarGetUbGlobal(bvar));
               continue;
            }

            /* skip products with non-binary variables if binreformbinaryonly is set */
            if( conshdlrdata->binreformbinaryonly && !SCIPvarIsBinary(bvar) )
            {
               SCIPdebugMsg(scip, "skip reform of <%s><%s> because second variable is not binary\n",
                  SCIPvarGetName(y), SCIPvarGetName(bvar));
               continue;
            }

            bilincoef = consdata->bilinterms[bilinidx].coef;
            assert(bilincoef != 0.0);

            /* get activity of bilincoef * x if y = 0 */
            SCIP_CALL( getImpliedBounds(scip, y, FALSE, bvar, &act0) );
            SCIPintervalMulScalar(SCIPinfinity(scip), &act0, act0, bilincoef);

            /* get activity of bilincoef * x if y = 1 */
            SCIP_CALL( getImpliedBounds(scip, y,  TRUE, bvar, &act1) );
            SCIPintervalMulScalar(SCIPinfinity(scip), &act1, act1, bilincoef);

            /* skip products that give rise to very large coefficients (big big-M's) */
            if( SCIPfeastol(scip) * REALABS(act0.inf) >= conshdlrdata->binreformmaxcoef || SCIPfeastol(scip) * REALABS(act0.sup) >= conshdlrdata->binreformmaxcoef )
            {
               SCIPdebugMsg(scip, "skip reform of %g<%s><%s> due to huge activity [%g,%g] for <%s> = 0.0\n",
                  bilincoef, SCIPvarGetName(y), SCIPvarGetName(bvar), SCIPintervalGetInf(act0), SCIPintervalGetSup(act0), SCIPvarGetName(y));
               continue;
            }
            if( SCIPfeastol(scip) * REALABS(act1.inf) >= conshdlrdata->binreformmaxcoef || SCIPfeastol(scip) * REALABS(act1.sup) >= conshdlrdata->binreformmaxcoef )
            {
               SCIPdebugMsg(scip, "skip reform of %g<%s><%s> due to huge activity [%g,%g] for <%s> = 1.0\n",
                  bilincoef, SCIPvarGetName(y), SCIPvarGetName(bvar), SCIPintervalGetInf(act1), SCIPintervalGetSup(act1), SCIPvarGetName(y));
               continue;
            }
            if( !SCIPisZero(scip, MIN(REALABS(act0.inf), REALABS(act0.sup))) &&
               SCIPfeastol(scip) * MAX(REALABS(act0.inf), REALABS(act0.sup)) / MIN(REALABS(act0.inf), REALABS(act0.sup)) >= conshdlrdata->binreformmaxcoef )
            {
               SCIPdebugMsg(scip, "skip reform of %g<%s><%s> due to huge activity ratio %g for <%s> = 0.0\n", bilincoef, SCIPvarGetName(y), SCIPvarGetName(bvar),
                  MAX(REALABS(act0.inf), REALABS(act0.sup)) / MIN(REALABS(act0.inf), REALABS(act0.sup)), SCIPvarGetName(y));
               continue;
            }
            if( !SCIPisZero(scip, MIN(REALABS(act1.inf), REALABS(act1.sup))) &&
               SCIPfeastol(scip) * MAX(REALABS(act1.inf), REALABS(act1.sup)) / MIN(REALABS(act1.inf), REALABS(act1.sup)) >= conshdlrdata->binreformmaxcoef )
            {
               SCIPdebugMsg(scip, "skip reform of %g<%s><%s> due to huge activity ratio %g for <%s> = 0.0\n", bilincoef, SCIPvarGetName(y), SCIPvarGetName(bvar),
                  MAX(REALABS(act1.inf), REALABS(act1.sup)) / MIN(REALABS(act1.inf), REALABS(act1.sup)), SCIPvarGetName(y));
               continue;
            }

            /* add bvar to x term */  
            xvars[nxvars] = bvar;
            xcoef[nxvars] = bilincoef;
            ++nxvars;

            /* update bounds on x term */
            SCIPintervalAdd(SCIPinfinity(scip), &xbndszero, xbndszero, act0);
            SCIPintervalAdd(SCIPinfinity(scip), &xbndsone,  xbndsone,  act1);

            if( REALABS(bilincoef) < mincoef )
               mincoef = ABS(bilincoef);
            if( REALABS(bilincoef) > maxcoef )
               maxcoef = ABS(bilincoef);

            /* update whether all coefficients will be integral and if so, compute their gcd */
            integral &= (SCIPvarGetType(bvar) < SCIP_VARTYPE_CONTINUOUS) && SCIPisIntegral(scip, bilincoef);  /*lint !e514 */
            if( integral )
            {
               if( nxvars == 1 )
                  gcd = (SCIP_Longint)SCIPround(scip, REALABS(bilincoef));
               else
                  gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)SCIPround(scip, REALABS(bilincoef)));
            }

            /* if bvar is initial, then also the auxiliary variable should be initial
             * if bvar is not removable, then also the auxiliary variable should not be removable
             */
            auxvarinitial |= SCIPvarIsInitial(bvar);
            auxvarremovable &= SCIPvarIsRemovable(bvar);

            /* remember that we have to remove this bilinear term later */
            assert(ntodelete < nbilinterms);
            todelete[ntodelete++] = bilinidx;
         }

         if( nxvars == 0 ) /* all (remaining) x_j seem to be unbounded */
            break;

         assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(xbndszero)));
         assert(!SCIPisInfinity(scip,  SCIPintervalGetSup(xbndszero)));
         assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(xbndsone)));
         assert(!SCIPisInfinity(scip,  SCIPintervalGetSup(xbndsone)));

#ifdef SCIP_DEBUG
         if( SCIPintervalGetInf(xbndszero) != SCIPintervalGetInf(xbndsone) || /*lint !e777*/
            +SCIPintervalGetSup(xbndszero) != SCIPintervalGetSup(xbndsone) ) /*lint !e777*/
         {
            SCIPdebugMsg(scip, "got different bounds for y = 0: [%g, %g] and y = 1: [%g, %g]\n", xbndszero.inf, xbndszero.sup, xbndsone.inf, xbndsone.sup);
         }
#endif

         if( nxvars == 1 && conshdlrdata->empathy4and >= 1 && SCIPvarIsBinary(xvars[0]) )
         {
            /* product of two binary variables, replace by auxvar and AND constraint */
            /* add auxiliary variable z */
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s_%s_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_IMPLINT,
                  auxvarinitial, auxvarremovable, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );

#ifdef WITH_DEBUG_SOLUTION
            if( SCIPdebugIsMainscip(scip) )
            {
               SCIP_Real var0val;
               SCIP_Real var1val;
               SCIP_CALL( SCIPdebugGetSolVal(scip, xvars[0], &var0val) );
               SCIP_CALL( SCIPdebugGetSolVal(scip, y, &var1val) );
               SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, var0val * var1val) );
            }
#endif

            /* add constraint z = x and y; need to be enforced, as it is not redundant w.r.t. existing constraints */
            xvars[1] = y;
            (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%sAND%s_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateConsAnd(scip, &auxcons, name, auxvar, 2, xvars,
                  SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
                  SCIPconsIsSeparated(cons), TRUE, TRUE,
                  SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMsg(scip, "added AND constraint: ");
            SCIPdebugPrintCons(scip, auxcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            ++*naddconss;

            /* add linear term coef*auxvar */
            SCIP_CALL( addLinearCoef(scip, cons, auxvar, xcoef[0]) );

            /* forget about auxvar */
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
         else
         {
            /* product of binary variable with more than one binary or with continuous variables or with binary and user
             * did not like AND -> replace by auxvar and linear constraints */
            SCIP_Real scale;

            /* scale auxiliary constraint by some nice value,
             * if all coefficients are integral, take a value that preserves integrality (-> gcd), so we can make the auxiliary variable impl. integer
             */
            if( integral )
            {
               scale = (SCIP_Real)gcd;
               assert(scale >= 1.0);
            }
            else if( nxvars == 1 )
            {
               /* scaling by the only coefficient gives auxiliary variable = x * y, which thus will be implicit integral provided y is not continuous */
               assert(mincoef == maxcoef);  /*lint !e777 */
               scale = mincoef;
               integral = SCIPvarGetType(xvars[0]) < SCIP_VARTYPE_CONTINUOUS;
            }
            else
            {
               scale = 1.0;
               if( maxcoef < 0.5 )
                  scale = maxcoef;
               if( mincoef > 2.0 )
                  scale = mincoef;
               if( scale != 1.0 )
                  scale = SCIPselectSimpleValue(scale / 2.0, 1.5 * scale, MAXDNOM);
            }
            assert(scale > 0.0);
            assert(!SCIPisInfinity(scip, scale));

            /* if x-term is always negative for y = 1, negate scale so we get a positive auxiliary variable; maybe this is better sometimes? */
            if( !SCIPisPositive(scip, SCIPintervalGetSup(xbndsone)) )
               scale = -scale;

            SCIPdebugMsg(scip, "binary reformulation using scale %g, nxvars = %d, integral = %u\n", scale, nxvars, integral);
            if( scale != 1.0 )
            {
               SCIPintervalDivScalar(SCIPinfinity(scip), &xbndszero, xbndszero, scale);
               SCIPintervalDivScalar(SCIPinfinity(scip), &xbndsone,  xbndsone, scale);
               for( k = 0; k < nxvars; ++k )
                  xcoef[k] /= scale;
            }

            /* add auxiliary variable z */
            if( nxvars == 1 )
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s_%s_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            else
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "prod%s_%s_more_%s", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, MIN(0., SCIPintervalGetInf(xbndsone)), MAX(0., SCIPintervalGetSup(xbndsone)),
                  0.0, integral ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_CONTINUOUS,
                  auxvarinitial, auxvarremovable, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );

            /* compute value of auxvar in debug solution */
#ifdef WITH_DEBUG_SOLUTION
            if( SCIPdebugIsMainscip(scip) )
            {
               SCIP_Real debugval;
               SCIP_Real varval;

               SCIP_CALL( SCIPdebugGetSolVal(scip, y, &varval) );
               if( SCIPisZero(scip, varval) )
               {
                  SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, 0.0) );
               }
               else
               {
                  assert(SCIPisEQ(scip, varval, 1.0));

                  debugval = 0.0;
                  for( k = 0; k < nxvars; ++k )
                  {
                     SCIP_CALL( SCIPdebugGetSolVal(scip, xvars[k], &varval) );
                     debugval += xcoef[k] * varval;
                  }
                  SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
               }
            }
#endif

            /* add auxiliary constraints
             * it seems to be advantageous to make the varbound constraints initial and the linear constraints not initial
             * maybe because it is more likely that a binary variable takes value 0 instead of 1, and thus the varbound constraints
             * are more often active, compared to the linear constraints added below
             * also, the varbound constraints are more sparse than the linear cons
             */
            if( SCIPisNegative(scip, SCIPintervalGetInf(xbndsone)) )
            {
               /* add 0 <= z - xbndsone.inf * y constraint (as varbound constraint), need to be enforced as not redundant */
               if( nxvars == 1 )
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_1", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               else
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_1", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetInf(xbndsone), 0.0, SCIPinfinity(scip),
                     SCIPconsIsInitial(cons) /*&& conshdlrdata->binreforminitial*/,
                     SCIPconsIsSeparated(cons), TRUE, TRUE,
                     SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMsg(scip, "added varbound constraint: ");
               SCIPdebugPrintCons(scip, auxcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
               ++*naddconss;
            }
            if( SCIPisPositive(scip, SCIPintervalGetSup(xbndsone)) )
            {
               /* add z - xbndsone.sup * y <= 0 constraint (as varbound constraint), need to be enforced as not redundant */
               if( nxvars == 1 )
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_2", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               else
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_2", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
               SCIP_CALL( SCIPcreateConsVarbound(scip, &auxcons, name, auxvar, y, -SCIPintervalGetSup(xbndsone), -SCIPinfinity(scip), 0.0,
                     SCIPconsIsInitial(cons) /*&& conshdlrdata->binreforminitial*/,
                     SCIPconsIsSeparated(cons), TRUE, TRUE,
                     SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                     SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebugMsg(scip, "added varbound constraint: ");
               SCIPdebugPrintCons(scip, auxcons, NULL);
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
               ++*naddconss;
            }

            /* add xbndszero.inf <= sum_i a_i*x_i + xbndszero.inf * y - z constraint, need to be enforced as not redundant */
            xvars[nxvars]   = y;
            xvars[nxvars+1] = auxvar;
            xcoef[nxvars]   = SCIPintervalGetInf(xbndszero);
            xcoef[nxvars+1] = -1;

            if( nxvars == 1 )
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_3", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            else
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_3", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, nxvars+2, xvars, xcoef, SCIPintervalGetInf(xbndszero), SCIPinfinity(scip),
                  SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
                  SCIPconsIsSeparated(cons), TRUE, TRUE,
                  SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMsg(scip, "added linear constraint: ");
            SCIPdebugPrintCons(scip, auxcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            ++*naddconss;

            /* add sum_i a_i*x_i + xbndszero.sup * y - z <= xbndszero.sup constraint, need to be enforced as not redundant */
            xcoef[nxvars] = SCIPintervalGetSup(xbndszero);

            if( nxvars == 1 )
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s_%s_4", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            else
               (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "linreform%s*%s*more_%s_4", SCIPvarGetName(y), SCIPvarGetName(xvars[0]), SCIPconsGetName(cons));
            SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, name, nxvars+2, xvars, xcoef, -SCIPinfinity(scip), SCIPintervalGetSup(xbndszero),
                  SCIPconsIsInitial(cons) && conshdlrdata->binreforminitial,
                  SCIPconsIsSeparated(cons), TRUE, TRUE,
                  SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
                  SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );
            SCIPdebugMsg(scip, "added linear constraint: ");
            SCIPdebugPrintCons(scip, auxcons, NULL);
            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            ++*naddconss;

            /* add linear term scale*auxvar to this constraint */
            SCIP_CALL( addLinearCoef(scip, cons, auxvar, scale) );

            /* forget about auxvar */
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
         }
      }
      while( j < nbilinterms );

      /* remove bilinear terms that have been replaced */
      SCIP_CALL( removeBilinearTermsPos(scip, cons, ntodelete, todelete) );
   }
   SCIPdebugMsg(scip, "resulting quadratic constraint: ");
   SCIPdebugPrintCons(scip, cons, NULL);

   SCIPfreeBufferArrayNull(scip, &xvars);
   SCIPfreeBufferArrayNull(scip, &xcoef);
   SCIPfreeBufferArrayNull(scip, &todelete);

   return SCIP_OKAY;
}

/** tries to automatically convert a quadratic constraint (or a part of it) into a more specific and more specialized constraint */
static
SCIP_RETCODE presolveUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_Bool*            upgraded,           /**< buffer to store whether constraint was upgraded */
   int*                  nupgdconss,         /**< buffer to increase if constraint was upgraded */
   int*                  naddconss,          /**< buffer to increase with number of additional constraints created during upgrade */
   SCIP_PRESOLTIMING     presoltiming        /**< current presolving timing */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real lincoef;
   SCIP_Real quadcoef;
   SCIP_Real lb;
   SCIP_Real ub;
   int nbinlin;
   int nbinquad;
   int nintlin;
   int nintquad;
   int nimpllin;
   int nimplquad;
   int ncontlin;
   int ncontquad;
   SCIP_Bool integral;
   int i;
   int j;
   SCIP_CONS** upgdconss;
   int upgdconsssize;
   int nupgdconss_;

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

   /* if there are no upgrade methods, we can also stop */
   if( conshdlrdata->nquadconsupgrades == 0 )
      return SCIP_OKAY;

   upgdconsssize = 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &upgdconss, upgdconsssize) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* calculate some statistics on quadratic constraint */
   nbinlin   = 0;
   nbinquad  = 0;
   nintlin   = 0;
   nintquad  = 0;
   nimpllin  = 0;
   nimplquad = 0;
   ncontlin  = 0;
   ncontquad = 0;
   integral  = TRUE;
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      var = consdata->linvars[i];
      lincoef = consdata->lincoefs[i];
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);
      assert(!SCIPisZero(scip, lincoef));

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nbinlin++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nintlin++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef);
         nimpllin++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisRelEQ(scip, lb, ub) && SCIPisIntegral(scip, lincoef * lb);
         ncontlin++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      var = consdata->quadvarterms[i].var;
      lincoef  = consdata->quadvarterms[i].lincoef;
      quadcoef = consdata->quadvarterms[i].sqrcoef;
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nbinquad++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nintquad++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, lincoef) && SCIPisIntegral(scip, quadcoef);
         nimplquad++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisRelEQ(scip, lb, ub) && SCIPisIntegral(scip, lincoef * lb + quadcoef * lb * lb);
         ncontquad++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
   }

   if( integral )
   {
      for( i = 0; i < consdata->nbilinterms && integral; ++i )
      {
         if( SCIPvarGetType(consdata->bilinterms[i].var1) < SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(consdata->bilinterms[i].var2) < SCIP_VARTYPE_CONTINUOUS )
            integral = integral && SCIPisIntegral(scip, consdata->bilinterms[i].coef);
         else
            integral = FALSE;
      }
   }

   /* call the upgrading methods */

   SCIPdebugMsg(scip, "upgrading quadratic constraint <%s> (%d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nquadconsupgrades);
   SCIPdebugMsg(scip, " binlin=%d binquad=%d intlin=%d intquad=%d impllin=%d implquad=%d contlin=%d contquad=%d integral=%u\n",
      nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral);
   SCIPdebugPrintCons(scip, cons, NULL);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nquadconsupgrades; ++i )
   {
      if( !conshdlrdata->quadconsupgrades[i]->active )
         continue;

      SCIP_CALL( conshdlrdata->quadconsupgrades[i]->quadconsupgd(scip, cons,
            nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral,
            &nupgdconss_, upgdconss, upgdconsssize, presoltiming) );

      while( nupgdconss_ < 0 )
      {
         /* upgrade function requires more memory: resize upgdconss and call again */
         assert(-nupgdconss_ > upgdconsssize);
         upgdconsssize = -nupgdconss_;
         SCIP_CALL( SCIPreallocBufferArray(scip, &upgdconss, -nupgdconss_) );

         SCIP_CALL( conshdlrdata->quadconsupgrades[i]->quadconsupgd(scip, cons,
               nbinlin, nbinquad, nintlin, nintquad, nimpllin, nimplquad, ncontlin, ncontquad, integral,
               &nupgdconss_, upgdconss, upgdconsssize, presoltiming) );

         assert(nupgdconss_ != 0);
      }

      if( nupgdconss_ > 0 )
      {
         /* got upgrade */
         SCIPdebugPrintCons(scip, cons, NULL);
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

/** helper function for presolveDisaggregate */
static
SCIP_RETCODE presolveDisaggregateMarkComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   quadvaridx,         /**< index of quadratic variable to mark */
   SCIP_HASHMAP*         var2component,      /**< variables to components mapping */
   int                   componentnr,        /**< the component number to mark to */
   int*                  componentsize       /**< buffer to store size of component (incremented by 1) */
   )
{
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_VAR* othervar;
   int othervaridx;
   int i;

   assert(consdata != NULL);
   assert(quadvaridx >= 0);
   assert(quadvaridx < consdata->nquadvars);
   assert(var2component != NULL);
   assert(componentnr >= 0);

   quadvarterm = &consdata->quadvarterms[quadvaridx];

   if( SCIPhashmapExists(var2component, quadvarterm->var) )
   {
      /* if we saw the variable before, then it should have the same component number */
      assert((int)(size_t)SCIPhashmapGetImage(var2component, quadvarterm->var) == componentnr);
      return SCIP_OKAY;
   }

   /* assign component number to variable */
   SCIP_CALL( SCIPhashmapInsert(var2component, quadvarterm->var, (void*)(size_t)componentnr) );
   ++*componentsize;

   /* assign same component number to all variables this variable is multiplied with */
   for( i = 0; i < quadvarterm->nadjbilin; ++i )
   {
      othervar = consdata->bilinterms[quadvarterm->adjbilin[i]].var1 == quadvarterm->var ?
         consdata->bilinterms[quadvarterm->adjbilin[i]].var2 : consdata->bilinterms[quadvarterm->adjbilin[i]].var1;
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, othervar, &othervaridx) );
      assert(othervaridx >= 0);
      SCIP_CALL( presolveDisaggregateMarkComponent(scip, consdata, othervaridx, var2component, componentnr, componentsize) );
   }

   return SCIP_OKAY;
}

/** merges components in variables connectivity graph */
static
SCIP_RETCODE presolveDisaggregateMergeComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_HASHMAP*         var2component,      /**< variables to component mapping */
   int                   nvars,              /**< number of variables */
   int*                  ncomponents,        /**< number of components */
   int*                  componentssize      /**< size of components */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_HASHMAPENTRY* entry;
   int maxncomponents;
   int* oldcompidx;
   int* newcompidx;
   int i;
   int oldcomponent;
   int newcomponent;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(var2component != NULL);
   assert(ncomponents != NULL);
   assert(componentssize != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   maxncomponents = conshdlrdata->maxdisaggrsize;
   assert(maxncomponents > 0);

   /* if already not too many components, then nothing to do */
   if( *ncomponents <= maxncomponents )
      return SCIP_OKAY;

   /*
   printf("component sizes before:");
   for( i = 0; i < *ncomponents; ++i )
      printf(" %d", componentssize[i]);
   printf("\n");
   */

   SCIP_CALL( SCIPallocBufferArray(scip, &oldcompidx, *ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newcompidx, *ncomponents) );

   for( i = 0; i < *ncomponents; ++i )
      oldcompidx[i] = i;

   switch( conshdlrdata->disaggrmergemethod )
   {
      case 's' :
         /* sort components by size, increasing order */
         SCIPsortIntInt(componentssize, oldcompidx, *ncomponents);
         break;
      case 'b' :
      case 'm' :
         /* sort components by size, decreasing order */
         SCIPsortDownIntInt(componentssize, oldcompidx, *ncomponents);
         break;
      default :
         SCIPerrorMessage("invalid value for constraints/quadratic/disaggrmergemethod parameter");
         return SCIP_PARAMETERWRONGVAL;
   }

   SCIPdebugMsg(scip, "%-30s: % 4d components of size % 4d to % 4d, median: % 4d\n", SCIPgetProbName(scip), *ncomponents, componentssize[0], componentssize[*ncomponents-1], componentssize[*ncomponents/2]);

   if( conshdlrdata->disaggrmergemethod == 'm' )
   {
      SCIP_Real targetsize;
      int count = 0;

      /* a minimal component size we should reach to have all components roughly the same size */
      targetsize = nvars / maxncomponents;  /*lint !e653*/
      for( i = 0; i < *ncomponents; ++i )
      {
         newcompidx[oldcompidx[i]] = i;
         count += componentssize[i];

         /* fill with small components until we reach targetsize
          * Since targetsize might be fractional, we also add another component if
          * the number of variables remaining (=nvars-count) is larger than
          * what we expect to put into the remaining components (=targetsize * (maxncomponents - i-1)).
          * Thus, from time to time, a component is made larger than the targetsize to avoid
          * having to add much into the last component.
          */
         while( i < *ncomponents-1 && (componentssize[i] + componentssize[*ncomponents-1] <= targetsize ||
            nvars - count > targetsize * (maxncomponents - i)) )
         {
            /* map last (=smallest) component to component i */
            newcompidx[oldcompidx[*ncomponents-1]] = i;

            /* increase size of component i accordingly */
            componentssize[i] += componentssize[*ncomponents-1];
            count += componentssize[*ncomponents-1];

            /* forget about last component */
            --*ncomponents;
         }
      }
      assert(count == nvars);
   }
   else
   {
      /* get inverse permutation */
      for( i = 0; i < *ncomponents; ++i )
         newcompidx[oldcompidx[i]] = i;
   }

   /* assign new component numbers to variables, cutting off at maxncomponents */
   for( i = 0; i < SCIPhashmapGetNEntries(var2component); ++i )
   {
      entry = SCIPhashmapGetEntry(var2component, i);
      if( entry == NULL )
         continue;

      oldcomponent = (int)(size_t)SCIPhashmapEntryGetImage(entry);

      newcomponent = newcompidx[oldcomponent];
      if( newcomponent >= maxncomponents )
      {
         newcomponent = maxncomponents-1;
         ++componentssize[maxncomponents-1];
      }

      SCIPhashmapEntrySetImage(entry, (void*)(size_t)newcomponent); /*lint !e571*/
   }
   if( *ncomponents > maxncomponents )
      *ncomponents = maxncomponents;

   /*
   printf("component sizes after :");
   for( i = 0; i < *ncomponents; ++i )
      printf(" %d", componentssize[i]);
   printf("\n");
   */

   SCIPfreeBufferArray(scip, &newcompidx);
   SCIPfreeBufferArray(scip, &oldcompidx);

   return SCIP_OKAY;
}

/** compute the next highest power of 2 for a 32-bit argument
 *
 * Source: https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
 *
 * @note Returns 0 for v=0.
 */
static
unsigned int nextPowerOf2(
   unsigned int          v                   /**< input */
   )
{
   v--;
   v |= v >> 1;
   v |= v >> 2;
   v |= v >> 4;
   v |= v >> 8;
   v |= v >> 16;
   v++;

   return v;
}


/** for quadratic constraints that consists of a sum of quadratic terms, disaggregates the sum into a set of constraints by introducing auxiliary variables
 *
 * Assume the quadratic constraint can be written in the form
 *   lhs <= b'x + sum_{k=1..p} q_k(x_k) <= rhs
 * where x_k denotes a subset of the variables in x and these subsets are pairwise disjunct
 * and q_k(.) is a quadratic form.
 * p is selected as large as possible, but to be <= conshdlrdata->maxdisaggrsize.
 *
 * Without additional scaling, the constraint is disaggregated into
 *   lhs <= b'x + sum_k c_k z_k <= rhs
 *   c_k z_k ~ q_k(x)
 * where "~" is either "<=", "==", or ">=", depending on whether lhs or rhs are infinite.
 * Further, c_k is chosen to be the maximal absolute value of the coefficients of the quadratic terms in q_k(x).
 * This is done to ensure that z_k takes values with a similar magnitute as the variables in x_k (better for separation).
 *
 * However, a solution of this disaggregated system can violate the original constraint by (p+1)*epsilon
 * (assuming unscaled violations are used, which is the default).
 * Therefore, all constraints are scaled by p+1:
 *   (p+1)*lhs <= (p+1)*b'x + (p+1) * sum_k c_k z_k <= (p+1) * rhs
 *   (p+1)*c_k z_k ~ (p+1)*q_k(x)
 */
static
SCIP_RETCODE presolveDisaggregate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   int*                  naddconss           /**< pointer to counter of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP* var2component;
   int* componentssize;
   int ncomponents;
   int i;
   int comp;
   SCIP_CONS** auxconss;
   SCIP_VAR** auxvars;
   SCIP_Real* auxcoefs;
#ifdef WITH_DEBUG_SOLUTION
   SCIP_Real* auxsolvals; /* value of auxiliary variable in debug solution */
#endif
   SCIP_Real scale;
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(naddconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* skip if constraint has been already disaggregated */
   if( consdata->isdisaggregated )
      return SCIP_OKAY;

   consdata->isdisaggregated = TRUE;

   /* make sure there are no quadratic variables without coefficients */
   SCIP_CALL( mergeAndCleanBilinearTerms(scip, cons) );
   SCIP_CALL( mergeAndCleanQuadVarTerms(scip, cons) );

   if( consdata->nquadvars <= 1 )
      return SCIP_OKAY;

   /* sort quadratic variable terms here, so we can later search in it without reordering the array */
   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   /* check how many quadratic terms with non-overlapping variables we have
    * in other words, the number of components in the sparsity graph of the quadratic term matrix
    */
   ncomponents = 0;
   SCIP_CALL( SCIPhashmapCreate(&var2component, SCIPblkmem(scip), consdata->nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &componentssize, consdata->nquadvars) );
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      /* if variable was marked already, skip it */
      if( SCIPhashmapExists(var2component, (void*)consdata->quadvarterms[i].var) )
         continue;

      /* start a new component with variable i */
      componentssize[ncomponents] = 0;
      SCIP_CALL( presolveDisaggregateMarkComponent(scip, consdata, i, var2component, ncomponents, componentssize + ncomponents) );
      ++ncomponents;
   }

   assert(ncomponents >= 1);

   /* if there is only one component, we cannot disaggregate
    * @todo we could still split the constraint into several while keeping the number of variables sharing several constraints as small as possible
    */
   if( ncomponents == 1 )
   {
      SCIPhashmapFree(&var2component);
      SCIPfreeBufferArray(scip, &componentssize);
      return SCIP_OKAY;
   }

   /* merge some components, if necessary */
   SCIP_CALL( presolveDisaggregateMergeComponents(scip, conshdlr, var2component, consdata->nquadvars, &ncomponents, componentssize) );

   SCIPfreeBufferArray(scip, &componentssize);

   /* scale all new constraints (ncomponents+1 many) by ncomponents+1 (or its next power of 2), so violations sum up to at most epsilon */
   scale = nextPowerOf2((unsigned int)ncomponents + 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &auxconss, ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxvars,  ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxcoefs, ncomponents) );
#ifdef WITH_DEBUG_SOLUTION
   SCIP_CALL( SCIPallocClearBufferArray(scip, &auxsolvals, ncomponents) );
#endif

   /* create auxiliary variables and empty constraints for each component */
   for( comp = 0; comp < ncomponents; ++comp )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_comp%d", SCIPconsGetName(cons), comp);

      SCIP_CALL( SCIPcreateVar(scip, &auxvars[comp], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), FALSE, NULL, NULL, NULL, NULL, NULL) );

      SCIP_CALL( SCIPcreateConsQuadratic2(scip, &auxconss[comp], name, 0, NULL, NULL, 0, NULL, 0, NULL,
            (SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : 0.0),
            (SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : 0.0),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), TRUE,
            TRUE, SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );

      auxcoefs[comp] = SCIPinfinity(scip);
   }

   /* add quadratic variables to each component constraint
    * delete adjacency information */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      assert(SCIPhashmapExists(var2component, consdata->quadvarterms[i].var));

      comp = (int)(size_t) SCIPhashmapGetImage(var2component, consdata->quadvarterms[i].var);
      assert(comp >= 0);
      assert(comp < ncomponents);

      /* add variable term to corresponding constraint */
      SCIP_CALL( SCIPaddQuadVarQuadratic(scip, auxconss[comp], consdata->quadvarterms[i].var, scale * consdata->quadvarterms[i].lincoef, scale * consdata->quadvarterms[i].sqrcoef) );

      /* reduce coefficient of aux variable */
      if( !SCIPisZero(scip, consdata->quadvarterms[i].lincoef) && ABS(consdata->quadvarterms[i].lincoef) < auxcoefs[comp] )
         auxcoefs[comp] = REALABS(consdata->quadvarterms[i].lincoef);
      if( !SCIPisZero(scip, consdata->quadvarterms[i].sqrcoef) && ABS(consdata->quadvarterms[i].sqrcoef) < auxcoefs[comp] )
         auxcoefs[comp] = REALABS(consdata->quadvarterms[i].sqrcoef);

      SCIPfreeBlockMemoryArray(scip, &consdata->quadvarterms[i].adjbilin, consdata->quadvarterms[i].adjbilinsize);
      consdata->quadvarterms[i].nadjbilin = 0;
      consdata->quadvarterms[i].adjbilinsize = 0;

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real debugvarval;

         SCIP_CALL( SCIPdebugGetSolVal(scip, consdata->quadvarterms[i].var, &debugvarval) );
         auxsolvals[comp] += consdata->quadvarterms[i].lincoef * debugvarval + consdata->quadvarterms[i].sqrcoef * debugvarval * debugvarval;
      }
#endif
   }

   /* add bilinear terms to each component constraint */
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      assert(SCIPhashmapExists(var2component, consdata->bilinterms[i].var1));
      assert(SCIPhashmapExists(var2component, consdata->bilinterms[i].var2));

      comp = (int)(size_t) SCIPhashmapGetImage(var2component, consdata->bilinterms[i].var1);
      assert(comp == (int)(size_t) SCIPhashmapGetImage(var2component, consdata->bilinterms[i].var2));
      assert(!SCIPisZero(scip, consdata->bilinterms[i].coef));

      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, auxconss[comp], 
            consdata->bilinterms[i].var1, consdata->bilinterms[i].var2, scale * consdata->bilinterms[i].coef) );

      if( ABS(consdata->bilinterms[i].coef) < auxcoefs[comp] )
         auxcoefs[comp] = ABS(consdata->bilinterms[i].coef);

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real debugvarval1;
         SCIP_Real debugvarval2;

         SCIP_CALL( SCIPdebugGetSolVal(scip, consdata->bilinterms[i].var1, &debugvarval1) );
         SCIP_CALL( SCIPdebugGetSolVal(scip, consdata->bilinterms[i].var2, &debugvarval2) );
         auxsolvals[comp] += consdata->bilinterms[i].coef * debugvarval1 * debugvarval2;
      }
#endif
   }

   /* forget about bilinear terms in cons */
   SCIPfreeBlockMemoryArray(scip, &consdata->bilinterms, consdata->bilintermssize);
   consdata->nbilinterms = 0;
   consdata->bilintermssize = 0;

   /* remove quadratic variable terms from cons */
   for( i = consdata->nquadvars - 1; i >= 0; --i )
   {
      SCIP_CALL( delQuadVarTermPos(scip, cons, i) );
   }
   assert(consdata->nquadvars == 0);

   /* scale remaining linear variables and sides by scale */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( chgLinearCoefPos(scip, cons, i, scale * consdata->lincoefs[i]) );
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      consdata->lhs *= scale;
      assert(!SCIPisInfinity(scip, -consdata->lhs) );
   }
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      consdata->rhs *= scale;
      assert(!SCIPisInfinity(scip, consdata->rhs) );
   }

   /* add auxiliary variables to auxiliary constraints
    * add aux vars and constraints to SCIP 
    * add aux vars to this constraint
    * set value of aux vars in debug solution, if any
    */
   SCIPdebugMsg(scip, "add %d constraints for disaggregation of quadratic constraint <%s>\n", ncomponents, SCIPconsGetName(cons));
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars + ncomponents) );
   for( comp = 0; comp < ncomponents; ++comp )
   {
      SCIP_CONSDATA* auxconsdata;

      SCIP_CALL( SCIPaddLinearVarQuadratic(scip, auxconss[comp], auxvars[comp], -scale * auxcoefs[comp]) );

      SCIP_CALL( SCIPaddVar(scip, auxvars[comp]) );

      SCIP_CALL( SCIPaddCons(scip, auxconss[comp]) );
      SCIPdebugPrintCons(scip, auxconss[comp], NULL);

      SCIP_CALL( addLinearCoef(scip, cons, auxvars[comp], scale * auxcoefs[comp]) );

      /* mark that the constraint should not further be disaggregated */
      auxconsdata = SCIPconsGetData(auxconss[comp]);
      assert(auxconsdata != NULL);
      auxconsdata->isdisaggregated = TRUE;

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         /* auxvar should take value from auxsolvals in debug solution, but we also scaled auxvar by auxcoefs[comp] */
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvars[comp], auxsolvals[comp] / auxcoefs[comp]) );
      }
#endif

      SCIP_CALL( SCIPreleaseCons(scip, &auxconss[comp]) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvars[comp]) );
   }
   *naddconss += ncomponents;

   SCIPdebugPrintCons(scip, cons, NULL);

   SCIPfreeBufferArray(scip, &auxconss);
   SCIPfreeBufferArray(scip, &auxvars);
   SCIPfreeBufferArray(scip, &auxcoefs);
#ifdef WITH_DEBUG_SOLUTION
   SCIPfreeBufferArray(scip, &auxsolvals);
#endif
   SCIPhashmapFree(&var2component);

   return SCIP_OKAY;
}

#ifdef CHECKIMPLINBILINEAR
/** checks if there are bilinear terms x*y with a binary variable x and an implication x = {0,1} -> y = 0
 *
 *  In this case, the bilinear term can be removed (x=0 case) or replaced by y (x=1 case).
 */
static
SCIP_RETCODE presolveApplyImplications(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   int*                  nbilinremoved       /**< buffer to store number of removed bilinear terms */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* x;
   SCIP_VAR* y;
   SCIP_INTERVAL implbnds;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nbilinremoved != NULL);

   *nbilinremoved = 0;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "apply implications in <%s>\n", SCIPconsGetName(cons));

   /* sort quadvarterms in case we need to search */
   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      x = consdata->quadvarterms[i].var;
      assert(x != NULL);

      if( consdata->quadvarterms[i].nadjbilin == 0 )
         continue;

      if( !SCIPvarIsBinary(x) )
         continue;

      if( !SCIPvarIsActive(x) )
         continue;

      if( SCIPvarGetNImpls(x, TRUE) == 0 && SCIPvarGetNImpls(x, FALSE) == 0 )
         continue;

      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         k = consdata->quadvarterms[i].adjbilin[j];
         assert(k >= 0);
         assert(k < consdata->nbilinterms);

         if( consdata->bilinterms[k].coef == 0.0 )
            continue;

         y = consdata->bilinterms[k].var1 == x ? consdata->bilinterms[k].var2 : consdata->bilinterms[k].var1;
         assert(x != y);

         SCIP_CALL( getImpliedBounds(scip, x, TRUE, y, &implbnds) );
         if( SCIPisZero(scip, implbnds.inf) && SCIPisZero(scip, implbnds.sup) )
         {
            /* if x = 1 implies y = 0, then we can remove the bilinear term x*y, since it is always 0
             * we only set the coefficient to 0.0 here and mark the bilinterms as not merged */
            SCIPdebugMsg(scip, "remove bilinear term %g<%s><%s> from <%s> due to implication\n", consdata->bilinterms[k].coef, SCIPvarGetName(x), SCIPvarGetName(y), SCIPconsGetName(cons));
            consdata->bilinterms[k].coef = 0.0;
            consdata->bilinmerged = FALSE;
            ++*nbilinremoved;
            continue;
         }

         SCIP_CALL( getImpliedBounds(scip, x, FALSE, y, &implbnds) );
         if( SCIPisZero(scip, implbnds.inf) && SCIPisZero(scip, implbnds.sup) )
         {
            /* if x = 0 implies y = 0, then we can replace the bilinear term x*y by y
             * we only move the coefficient to the linear coef of y here and mark the bilinterms as not merged */
            SCIPdebugMsg(scip, "replace bilinear term %g<%s><%s> by %g<%s> in <%s> due to implication\n", consdata->bilinterms[k].coef, SCIPvarGetName(x), SCIPvarGetName(y), consdata->bilinterms[k].coef, SCIPvarGetName(y), SCIPconsGetName(cons));
            assert(consdata->quadvarssorted);
            SCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, cons, y, consdata->bilinterms[k].coef) );
            consdata->bilinterms[k].coef = 0.0;
            consdata->bilinmerged = FALSE;
            ++*nbilinremoved;
         }
      }
   }

   if( *nbilinremoved > 0 )
   {
      SCIP_CALL( mergeAndCleanBilinearTerms(scip, cons) );

      /* invalidate nonlinear row */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }

      consdata->ispropagated  = FALSE;
      consdata->ispresolved   = FALSE;
      consdata->iscurvchecked = FALSE;
   }

   consdata->isimpladded = FALSE;

   return SCIP_OKAY;
}
#endif

/** checks a quadratic constraint for convexity and/or concavity without checking multivariate functions */
static
void checkCurvatureEasy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_Bool*            determined,         /**< pointer to store whether the curvature could be determined */
   SCIP_Bool             checkmultivariate   /**< whether curvature will be checked later on for multivariate functions */
   )
{
   SCIP_CONSDATA* consdata;
   int nquadvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(determined != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nquadvars = consdata->nquadvars;
   *determined = TRUE;

   if( consdata->iscurvchecked )
      return;

   SCIPdebugMsg(scip, "Checking curvature of constraint <%s> without multivariate functions\n", SCIPconsGetName(cons));

   consdata->maxnonconvexity = 0.0;
   if( nquadvars == 1 )
   {
      assert(consdata->nbilinterms == 0);
      consdata->isconvex      = !SCIPisNegative(scip, consdata->quadvarterms[0].sqrcoef);
      consdata->isconcave     = !SCIPisPositive(scip, consdata->quadvarterms[0].sqrcoef);
      consdata->iscurvchecked = TRUE;

      if( !SCIPisInfinity(scip, -consdata->lhs) && consdata->quadvarterms[0].sqrcoef > 0.0 )
         consdata->maxnonconvexity =  consdata->quadvarterms[0].sqrcoef;
      if( !SCIPisInfinity(scip,  consdata->rhs) && consdata->quadvarterms[0].sqrcoef < 0.0 )
         consdata->maxnonconvexity = -consdata->quadvarterms[0].sqrcoef;
   }
   else if( nquadvars == 0 )
   {
      consdata->isconvex = TRUE;
      consdata->isconcave = TRUE;
      consdata->iscurvchecked = TRUE;
   }
   else if( consdata->nbilinterms == 0 )
   {
      int v;

      consdata->isconvex = TRUE;
      consdata->isconcave = TRUE;

      for( v = nquadvars - 1; v >= 0; --v )
      {
         consdata->isconvex  = consdata->isconvex  && !SCIPisNegative(scip, consdata->quadvarterms[v].sqrcoef);
         consdata->isconcave = consdata->isconcave && !SCIPisPositive(scip, consdata->quadvarterms[v].sqrcoef);

         if( !SCIPisInfinity(scip, -consdata->lhs) &&  consdata->quadvarterms[v].sqrcoef > consdata->maxnonconvexity )
            consdata->maxnonconvexity =  consdata->quadvarterms[0].sqrcoef;
         if( !SCIPisInfinity(scip,  consdata->rhs) && -consdata->quadvarterms[v].sqrcoef > consdata->maxnonconvexity )
            consdata->maxnonconvexity = -consdata->quadvarterms[0].sqrcoef;
      }

      consdata->iscurvchecked = TRUE;
   }
   else if( !checkmultivariate )
   {
      consdata->isconvex  = FALSE;
      consdata->isconcave = FALSE;
      consdata->iscurvchecked = TRUE;
      consdata->maxnonconvexity = SCIPinfinity(scip);
   }
   else
      *determined = FALSE;
}

/** checks a quadratic constraint for convexity and/or concavity */
static
SCIP_RETCODE checkCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_Bool             checkmultivariate   /**< whether curvature should also be checked for multivariate functions */
   )
{
   SCIP_CONSDATA* consdata;
   double*        matrix;
   SCIP_HASHMAP*  var2index;
   int            i;
   int            n;
   int            nn;
   int            row;
   int            col;
   double*        alleigval;
   SCIP_Bool      determined;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   n = consdata->nquadvars;

   if( consdata->iscurvchecked )
      return SCIP_OKAY;

   /* easy checks for curvature detection */
   checkCurvatureEasy(scip, cons, &determined, checkmultivariate);

   /* if curvature was already detected stop */
   if( determined )
   {
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "Checking curvature of constraint <%s> with multivariate functions\n", SCIPconsGetName(cons));

   if( n == 2 )
   {
      SCIP_Real tracehalf;
      SCIP_Real discriminantroot;

      /* compute eigenvalues by hand */
      assert(consdata->nbilinterms == 1);
      consdata->isconvex =
         consdata->quadvarterms[0].sqrcoef >= 0 &&
         consdata->quadvarterms[1].sqrcoef >= 0 &&
         4 * consdata->quadvarterms[0].sqrcoef * consdata->quadvarterms[1].sqrcoef >= consdata->bilinterms[0].coef * consdata->bilinterms[0].coef;
      consdata->isconcave = 
         consdata->quadvarterms[0].sqrcoef <= 0 &&
         consdata->quadvarterms[1].sqrcoef <= 0 &&
         4 * consdata->quadvarterms[0].sqrcoef * consdata->quadvarterms[1].sqrcoef >= consdata->bilinterms[0].coef * consdata->bilinterms[0].coef;

      /* store largest eigenvalue causing nonconvexity according to sides */
      tracehalf = (consdata->quadvarterms[0].sqrcoef + consdata->quadvarterms[1].sqrcoef) / 2.0;
      discriminantroot = consdata->quadvarterms[0].sqrcoef * consdata->quadvarterms[1].sqrcoef - SQR(consdata->bilinterms[0].coef / 2.0);
      discriminantroot = SQR(tracehalf) - discriminantroot;
      assert(!SCIPisNegative(scip, discriminantroot));
      discriminantroot = SQRT(MAX(0.0, discriminantroot));

      consdata->maxnonconvexity = 0.0;
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         consdata->maxnonconvexity = MAX(consdata->maxnonconvexity, tracehalf + discriminantroot);
      if( !SCIPisInfinity(scip, consdata->rhs) )
         consdata->maxnonconvexity = MAX(consdata->maxnonconvexity, discriminantroot - tracehalf);

      consdata->iscurvchecked = TRUE;
      return SCIP_OKAY;
   }

   /* do not check curvature if n is too large */
   nn = n * n;
   if( nn < 0 || (unsigned) (int) nn > UINT_MAX / sizeof(SCIP_Real) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "cons_quadratic - n is too large to check the curvature\n");
      consdata->isconvex = FALSE;
      consdata->isconcave = FALSE;
      consdata->iscurvchecked = TRUE;
      consdata->maxnonconvexity = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   /* lower triangular of quadratic term matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, nn) );
   BMSclearMemoryArray(matrix, nn);

   consdata->isconvex  = TRUE;
   consdata->isconcave = TRUE;
   consdata->maxnonconvexity = 0.0;

   SCIP_CALL( SCIPhashmapCreate(&var2index, SCIPblkmem(scip), n) );
   for( i = 0; i < n; ++i )
   {
      if( consdata->quadvarterms[i].nadjbilin > 0 )
      {
         SCIP_CALL( SCIPhashmapInsert(var2index, consdata->quadvarterms[i].var, (void*)(size_t)i) );
         matrix[i*n + i] = consdata->quadvarterms[i].sqrcoef;
      }
      else
      {
         /* if pure square term, then update maximal nonconvex eigenvalue, as it will not be considered in lapack call below */
         if( !SCIPisInfinity(scip, -consdata->lhs) && consdata->quadvarterms[i].sqrcoef > consdata->maxnonconvexity )
            consdata->maxnonconvexity = consdata->quadvarterms[i].sqrcoef;
         if( !SCIPisInfinity(scip, consdata->rhs) && -consdata->quadvarterms[i].sqrcoef > consdata->maxnonconvexity )
            consdata->maxnonconvexity = -consdata->quadvarterms[i].sqrcoef;
      }
      /* nonzero elements on diagonal tell a lot about convexity/concavity */
      if( SCIPisNegative(scip, consdata->quadvarterms[i].sqrcoef) )
         consdata->isconvex  = FALSE;
      if( SCIPisPositive(scip, consdata->quadvarterms[i].sqrcoef) )
         consdata->isconcave = FALSE;
   }

   /* skip lapack call, if we know already that we are indefinite
    * NOTE: this will leave out updating consdata->maxnonconvexity, so that it only provides a lower bound in this case
    */
   if( !consdata->isconvex && !consdata->isconcave )
   {
      SCIPfreeBufferArray(scip, &matrix);
      SCIPhashmapFree(&var2index);
      consdata->iscurvchecked = TRUE;
      /* make sure that maxnonconvexity is strictly different from zero if nonconvex
       * TODO one could think about doing some eigenvalue estimation here (Gershgorin)
       */
      consdata->maxnonconvexity = MAX(1000.0, consdata->maxnonconvexity);
      return SCIP_OKAY;
   }

   if( SCIPisIpoptAvailableIpopt() )
   {
      for( i = 0; i < consdata->nbilinterms; ++i )
      {
         assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var1));
         assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var2));
         row = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var1);
         col = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var2);
         if( row < col )
            matrix[row * n + col] = consdata->bilinterms[i].coef/2;
         else
            matrix[col * n + row] = consdata->bilinterms[i].coef/2;
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &alleigval, n) );
      /* @todo Can we compute only min and max eigen value?
       * @todo Can we estimate the numerical error? 
       * @todo Trying a cholesky factorization may be much faster. 
       */
      if( LapackDsyev(FALSE, n, matrix, alleigval) != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Failed to compute eigenvalues of quadratic coefficient matrix of constraint %s. Assuming matrix is indefinite.\n", SCIPconsGetName(cons));
         consdata->isconvex = FALSE;
         consdata->isconcave = FALSE;
      }
      else
      {
         /* deconvexification reformulates a stricly convex quadratic function in binaries such that it becomes not-strictly convex
          * by adding the -lambda*(x^2-x) terms for lambda the smallest eigenvalue of the matrix
          * the result is still a convex form "but less so" (ref. papers by Guignard et.al.), but with hopefully tighter value for the continuous relaxation
          */
#ifdef DECONVEXIFY
         SCIP_Bool allbinary;
         printf("cons <%s>[%g,%g] spectrum = [%g,%g]\n", SCIPconsGetName(cons), consdata->lhs, consdata->rhs, alleigval[0], alleigval[n-1]);
#endif
         consdata->isconvex  &= !SCIPisNegative(scip, alleigval[0]);   /*lint !e514*/
         consdata->isconcave &= !SCIPisPositive(scip, alleigval[n-1]); /*lint !e514*/
         consdata->iscurvchecked = TRUE;
#ifdef DECONVEXIFY
         for( i = 0; i < consdata->nquadvars; ++i )
            if( !SCIPvarIsBinary(consdata->quadvarterms[i].var) )
               break;
         allbinary = i == consdata->nquadvars;

         if( !SCIPisInfinity(scip, consdata->rhs) && alleigval[0] > 0.1 && allbinary )
         {
            printf("deconvexify cons <%s> by shifting hessian by %g\n", SCIPconsGetName(cons), alleigval[0]);
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               consdata->quadvarterms[i].sqrcoef -= alleigval[0];
               consdata->quadvarterms[i].lincoef += alleigval[0];
            }
         }

         if( !SCIPisInfinity(scip, consdata->lhs) && alleigval[n-1] < -0.1 && allbinary )
         {
            printf("deconcavify cons <%s> by shifting hessian by %g\n", SCIPconsGetName(cons), alleigval[n-1]);
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               consdata->quadvarterms[i].sqrcoef -= alleigval[n-1];
               consdata->quadvarterms[i].lincoef += alleigval[n-1];
            }
         }
#endif
      }

      /* update largest eigenvalue causing nonconvexity according to sides */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         consdata->maxnonconvexity = MAX(consdata->maxnonconvexity, alleigval[n-1]);
      if( !SCIPisInfinity(scip, consdata->rhs) )
         consdata->maxnonconvexity = MAX(consdata->maxnonconvexity, -alleigval[0]);

      SCIPfreeBufferArray(scip, &alleigval);
   }
   else
   {
      consdata->isconvex = FALSE;
      consdata->isconcave = FALSE;
      consdata->iscurvchecked = TRUE; /* set to TRUE since it does not help to repeat this procedure again and again (that will not bring Ipopt in) */
      consdata->maxnonconvexity = SCIPinfinity(scip);
   }

   SCIPhashmapFree(&var2index);
   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}

/** check whether indefinite constraint function is factorable and store corresponding coefficients */
static
SCIP_RETCODE checkFactorable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_BILINTERM* bilinterm;
   SCIP_CONSDATA* consdata;
   SCIP_Real* a;
   SCIP_Real* eigvals;
   SCIP_Real sigma1;
   SCIP_Real sigma2;
   SCIP_Bool success;
   int n;
   int i;
   int idx1;
   int idx2;
   int posidx;
   int negidx;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->factorleft == NULL);
   assert(consdata->factorright == NULL);

   /* we don't need this if there are no bilinear terms */
   if( consdata->nbilinterms == 0 )
      return SCIP_OKAY;

   /* write constraint as lhs <= linear + x'^T A x' <= rhs where x' = (x,1) and
    * A = ( Q     b/2 )
    *     ( b^T/2  0  )
    * compute an eigenvalue factorization of A and check if there are one positive and one negative eigenvalue
    * if so, then let sigma1^2 and -sigma2^2 be these eigenvalues and v1 and v2 be the first two rows of the inverse eigenvector matrix
    * thus, x'^T A x' = sigma1^2 (v1^T x')^2 - sigma2^2 (v2^T x')^2
    *                 = (sigma1 (v1^T x') - sigma2 (v2^T x')) * (sigma1 (v1^T x') + sigma2 (v2^T x'))
    * we then store sigma1 v1^T - sigma2 v2^T as left factor coef, and sigma1 v1^T + sigma2 v2^T as right factor coef
    */

   /* if we already know that there are only positive or only negative eigenvalues, then don't try */
   if( consdata->iscurvchecked && (consdata->isconvex || consdata->isconcave) )
      return SCIP_OKAY;

   n = consdata->nquadvars + 1;

   /* @todo handle case n=3 explicitly */

   /* skip too large matrices */
   if( n > 50 )
      return SCIP_OKAY;

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

   SCIP_CALL( SCIPallocBufferArray(scip, &a, n*n) );
   BMSclearMemoryArray(a, n*n);

   /* set lower triangular entries of A corresponding to bilinear terms */
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      bilinterm = &consdata->bilinterms[i];

      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var1, &idx1) );
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var2, &idx2) );
      assert(idx1 >= 0);
      assert(idx2 >= 0);
      assert(idx1 != idx2);

      a[MIN(idx1,idx2) * n + MAX(idx1,idx2)] = bilinterm->coef / 2.0;
   }

   /* set lower triangular entries of A corresponding to square and linear terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      a[i*n + i]   = consdata->quadvarterms[i].sqrcoef;
      a[i*n + n-1] = consdata->quadvarterms[i].lincoef / 2.0;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &eigvals, n) );
   if( LapackDsyev(TRUE, n, a, eigvals) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors of augmented quadratic form matrix for constraint <%s>.\n", SCIPconsGetName(cons));
      goto CLEANUP;
   }

   /* check if there is exactly one positive and one negative eigenvalue */
   posidx = -1;
   negidx = -1;
   for( i = 0; i < n; ++i )
   {
      if( SCIPisPositive(scip, eigvals[i]) )
      {
         if( posidx == -1 )
            posidx = i;
         else
            break;
      }
      else if( SCIPisNegative(scip, eigvals[i]) )
      {
         if( negidx == -1 )
            negidx = i;
         else
            break;
      }
   }
   if( i < n || posidx == -1 || negidx == -1 )
   {
      SCIPdebugMsg(scip, "Augmented quadratic form of constraint <%s> is not factorable.\n", SCIPconsGetName(cons));
      goto CLEANUP;
   }
   assert(SCIPisPositive(scip, eigvals[posidx]));
   assert(SCIPisNegative(scip, eigvals[negidx]));

   /* compute factorleft and factorright */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->factorleft,  consdata->nquadvars + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->factorright, consdata->nquadvars + 1) );

   /* eigenvectors are stored in a, inverse eigenvector matrix is transposed of a
    * it seems that v1 and v2 are at &a[posidx*n] and &a[negidx*n]
    */
   sigma1 = sqrt( eigvals[posidx]);
   sigma2 = sqrt(-eigvals[negidx]);
   for( i = 0; i < n; ++i )
   {
      consdata->factorleft[i]  = sigma1 * a[posidx * n + i] - sigma2 * a[negidx * n + i];
      consdata->factorright[i] = sigma1 * a[posidx * n + i] + sigma2 * a[negidx * n + i];
      /* set almost-zero elements to zero */
      if( SCIPisZero(scip, consdata->factorleft[i]) )
         consdata->factorleft[i] = 0.0;
      if( SCIPisZero(scip, consdata->factorright[i]) )
         consdata->factorright[i] = 0.0;
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "constraint <%s> has factorable quadratic form: (%g", SCIPconsGetName(cons), consdata->factorleft[n-1]);
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->factorleft[i] != 0.0 )
         SCIPdebugMsgPrint(scip, " %+g<%s>", consdata->factorleft[i], SCIPvarGetName(consdata->quadvarterms[i].var));
   }
   SCIPdebugMsgPrint(scip, ") * (%g", consdata->factorright[n-1]);
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->factorright[i] != 0.0 )
         SCIPdebugMsgPrint(scip, " %+g<%s>", consdata->factorright[i], SCIPvarGetName(consdata->quadvarterms[i].var));
   }
   SCIPdebugMsgPrint(scip, ")\n");
#endif

   /* check whether factorleft * factorright^T is matrix of augmented quadratic form
    * we check here only the nonzero entries from the quadratic form
    */
   success = TRUE;

   /* check bilinear terms */
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      bilinterm = &consdata->bilinterms[i];

      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var1, &idx1) );
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, bilinterm->var2, &idx2) );

      if( !SCIPisRelEQ(scip, consdata->factorleft[idx1] * consdata->factorright[idx2] + consdata->factorleft[idx2] * consdata->factorright[idx1], bilinterm->coef) )
      {
         success = FALSE;
         break;
      }
   }

   /* set lower triangular entries of A corresponding to square and linear terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( !SCIPisRelEQ(scip, consdata->factorleft[i] * consdata->factorright[i], consdata->quadvarterms[i].sqrcoef) )
      {
         success = FALSE;
         break;
      }

      if( !SCIPisRelEQ(scip, consdata->factorleft[n-1] * consdata->factorright[i] + consdata->factorleft[i] * consdata->factorright[n-1], consdata->quadvarterms[i].lincoef) )
      {
         success = FALSE;
         break;
      }
   }

   if( !success )
   {
      SCIPdebugMsg(scip, "Factorization not accurate enough. Dropping it.\n");
      SCIPfreeBlockMemoryArray(scip, &consdata->factorleft,  consdata->nquadvars + 1);
      SCIPfreeBlockMemoryArray(scip, &consdata->factorright, consdata->nquadvars + 1);
   }

 CLEANUP:
   SCIPfreeBufferArray(scip, &a);
   SCIPfreeBufferArray(scip, &eigvals);

   return SCIP_OKAY;
}

/** computes activity and violation of a constraint
 *
 * If solution violates bounds by more than feastol, the violation is still computed, but *solviolbounds is set to TRUE
 */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool*            solviolbounds       /**< buffer to store whether quadratic variables in solution are outside their bounds by more than feastol */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real varval;
   SCIP_Real varval2;
   SCIP_Real absviol;
   SCIP_Real relviol;
   SCIP_VAR* var;
   SCIP_VAR* var2;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(solviolbounds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *solviolbounds = FALSE;
   consdata->activity = 0.0;
   consdata->lhsviol = 0.0;
   consdata->rhsviol = 0.0;

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_Real activity;

      var = consdata->linvars[i];
      varval = SCIPgetSolVal(scip, sol, var);
      activity = consdata->lincoefs[i] * varval;

      /* the contribution of a variable with |varval| = +inf is +inf when activity > 0.0, -inf when activity < 0.0, and
       * 0.0 otherwise
       */
      if( SCIPisInfinity(scip, REALABS(varval)) )
      {
         if( activity > 0.0 && !SCIPisInfinity(scip, consdata->rhs) )
         {
            consdata->activity = SCIPinfinity(scip);
            consdata->rhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }

         if( activity < 0.0 && !SCIPisInfinity(scip, -consdata->lhs) )
         {
            consdata->activity = -SCIPinfinity(scip);
            consdata->lhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }
      }

      consdata->activity += activity;
   }

   for( j = 0; j < consdata->nquadvars; ++j )
   {
      SCIP_Real activity;

      var = consdata->quadvarterms[j].var;
      varval = SCIPgetSolVal(scip, sol, var);
      activity = (consdata->quadvarterms[j].lincoef + consdata->quadvarterms[j].sqrcoef * varval) * varval;

      /* the contribution of a variable with |varval| = +inf is +inf when activity > 0.0, -inf when activity < 0.0, and
       * 0.0 otherwise
       */
      if( SCIPisInfinity(scip, REALABS(varval)) )
      {
         if( activity > 0.0 && !SCIPisInfinity(scip, consdata->rhs) )
         {
            consdata->activity = SCIPinfinity(scip);
            consdata->rhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }

         if( activity < 0.0 && !SCIPisInfinity(scip, -consdata->lhs) )
         {
            consdata->activity = -SCIPinfinity(scip);
            consdata->lhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }
      }

      /* project onto local box, in case the LP solution is slightly outside the bounds (which is not our job to enforce) */
      if( sol == NULL )
      {
         /* with non-initial columns, variables can shortly be a column variable before entering the LP and have value 0.0 in this case, which might violated the variable bounds */
         if( (!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) && !SCIPisFeasGE(scip, varval, SCIPvarGetLbLocal(var))) ||
             (!SCIPisInfinity(scip,  SCIPvarGetUbLocal(var)) && !SCIPisFeasLE(scip, varval, SCIPvarGetUbLocal(var))) )
            *solviolbounds = TRUE;
         else
         {
            varval = MAX(SCIPvarGetLbLocal(var), MIN(SCIPvarGetUbLocal(var), varval));
            activity = (consdata->quadvarterms[j].lincoef + consdata->quadvarterms[j].sqrcoef * varval) * varval;
         }
      }

      consdata->activity += activity;
   }

   for( j = 0; j < consdata->nbilinterms; ++j )
   {
      SCIP_Real activity;

      var = consdata->bilinterms[j].var1;
      var2 = consdata->bilinterms[j].var2;
      varval = SCIPgetSolVal(scip, sol, var);
      varval2 = SCIPgetSolVal(scip, sol, var2);

      /* project onto local box, in case the LP solution is slightly outside the bounds (which is not our job to enforce) */
      if( sol == NULL )
      {
         /* with non-initial columns, variables can shortly be a column variable before entering the LP and have value 0.0 in this case, which might violated the variable bounds */
         if( (!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) && !SCIPisFeasGE(scip, varval, SCIPvarGetLbLocal(var))) ||
             (!SCIPisInfinity(scip,  SCIPvarGetUbLocal(var)) && !SCIPisFeasLE(scip, varval, SCIPvarGetUbLocal(var))) )
            *solviolbounds = TRUE;
         else
            varval = MAX(SCIPvarGetLbLocal(var), MIN(SCIPvarGetUbLocal(var), varval));

         /* with non-initial columns, variables can shortly be a column variable before entering the LP and have value 0.0 in this case, which might violated the variable bounds */
         if( (!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var2)) && !SCIPisFeasGE(scip, varval2, SCIPvarGetLbLocal(var2))) ||
             (!SCIPisInfinity(scip,  SCIPvarGetUbLocal(var2)) && !SCIPisFeasLE(scip, varval2, SCIPvarGetUbLocal(var2))) )
            *solviolbounds = TRUE;
         else
            varval2 = MAX(SCIPvarGetLbLocal(var2), MIN(SCIPvarGetUbLocal(var2), varval2));
      }

      activity = consdata->bilinterms[j].coef * varval * varval2;

      /* consider var*var2 as a new variable and handle it as it would appear linearly */
      if( SCIPisInfinity(scip, REALABS(varval*varval2)) )
      {
         if( activity > 0.0 && !SCIPisInfinity(scip, consdata->rhs) )
         {
            consdata->activity = SCIPinfinity(scip);
            consdata->rhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }

         if( activity < 0.0 && !SCIPisInfinity(scip, -consdata->lhs) )
         {
            consdata->activity = -SCIPinfinity(scip);
            consdata->lhsviol = SCIPinfinity(scip);
            return SCIP_OKAY;
         }
      }

      consdata->activity += activity;
   }

   absviol = 0.0;
   relviol = 0.0;
   /* compute absolute violation left hand side */
   if( consdata->activity < consdata->lhs && !SCIPisInfinity(scip, -consdata->lhs) )
   {
      consdata->lhsviol = consdata->lhs - consdata->activity;
      absviol = consdata->lhsviol;
      relviol = SCIPrelDiff(consdata->lhs, consdata->activity);
   }
   else
      consdata->lhsviol = 0.0;

   /* compute absolute violation right hand side */
   if( consdata->activity > consdata->rhs && !SCIPisInfinity(scip,  consdata->rhs) )
   {
      consdata->rhsviol = consdata->activity - consdata->rhs;
      absviol = consdata->rhsviol;
      relviol = SCIPrelDiff(consdata->activity, consdata->rhs);
   }
   else
      consdata->rhsviol = 0.0;

   /* update absolute and relative violation of the solution */
   if( sol != NULL )
      SCIPupdateSolConsViolation(scip, sol, absviol, relviol);

   return SCIP_OKAY;
}

/** computes violation of a set of constraints */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool*            solviolbounds,      /**< buffer to store whether quadratic variables in solution are outside their bounds by more than feastol in some constraint */
   SCIP_CONS**           maxviolcon          /**< buffer to store constraint with largest violation, or NULL if solution is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      viol;
   SCIP_Real      maxviol;
   SCIP_Bool      solviolbounds1;
   int            c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(solviolbounds != NULL);
   assert(maxviolcon != NULL);

   *solviolbounds = FALSE;
   *maxviolcon = NULL;

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      SCIP_CALL( computeViolation(scip, conss[c], sol, &solviolbounds1) );
      *solviolbounds |= solviolbounds1;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if( viol > maxviol && SCIPisGT(scip, viol, SCIPfeastol(scip)) )
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }
   }

   return SCIP_OKAY;
}


/** index comparison method for bilinear terms */
static
SCIP_DECL_SORTINDCOMP(bilinTermComp2)
{  /*lint --e{715}*/
   SCIP_BILINTERM* bilinterms = (SCIP_BILINTERM*)dataptr;
   int var1cmp;

   assert(bilinterms != NULL);

   var1cmp = SCIPvarCompare(bilinterms[ind1].var1, bilinterms[ind2].var1);
   if( var1cmp != 0 )
      return var1cmp;

   return SCIPvarCompare(bilinterms[ind1].var2, bilinterms[ind2].var2);
}

/** volume comparison method for bilinear terms; prioritizes bilinear products with a larger volume */
static
SCIP_DECL_SORTINDCOMP(bilinTermCompVolume)
{  /*lint --e{715}*/
   SCIP_BILINTERM* bilinterms = (SCIP_BILINTERM*)dataptr;
   SCIP_Real vol1;
   SCIP_Real vol2;

   assert(bilinterms != NULL);

   vol1 = (SCIPvarGetUbLocal(bilinterms[ind1].var1) - SCIPvarGetLbLocal(bilinterms[ind1].var1))
      * (SCIPvarGetUbLocal(bilinterms[ind1].var2) - SCIPvarGetLbLocal(bilinterms[ind1].var2));
   vol2 = (SCIPvarGetUbLocal(bilinterms[ind2].var1) - SCIPvarGetLbLocal(bilinterms[ind2].var1))
      * (SCIPvarGetUbLocal(bilinterms[ind2].var2) - SCIPvarGetLbLocal(bilinterms[ind2].var2));

   if( vol1 > vol2 )
      return -1;
   else if( vol1 < vol2 )
      return 1;
   return bilinTermComp2(dataptr, ind1, ind2);
}

/** helper function to sort all bilinear terms in the constraint handler data */
static
SCIP_RETCODE sortAllBilinTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BILINTERM*       bilinterms,         /**< array containing all bilinear terms */
   int                   nbilinterms,        /**< total number of bilinear terms */
   SCIP_CONS**           bilinconss,         /**< array for mapping each term to its constraint */
   int*                  bilinposs           /**< array for mapping each term to its position in the corresponding
                                               *  bilinconss constraint */
   )
{
   int* perm;
   int  i;
   int  nexti;
   int  v;
   SCIP_BILINTERM bilinterm;
   SCIP_CONS* bilincons;
   int bilinpos;

   assert(scip != NULL);
   assert(bilinterms != NULL);
   assert(nbilinterms > 0);
   assert(bilinconss != NULL);
   assert(bilinposs != NULL);

   /* get temporary memory to store the sorted permutation and the inverse permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nbilinterms) );

   /* call quicksort */
   SCIPsort(perm, bilinTermCompVolume, (void*)bilinterms, nbilinterms);

   /* permute the bilinear terms according to the resulting permutation */
   for( v = 0; v < nbilinterms; ++v )
   {
      if( perm[v] != v )
      {
         bilinterm = bilinterms[v];
         bilincons = bilinconss[v];
         bilinpos = bilinposs[v];

         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < nbilinterms);
            assert(perm[i] != i);

            bilinterms[i] = bilinterms[perm[i]];
            bilinconss[i] = bilinconss[perm[i]];
            bilinposs[i] = bilinposs[perm[i]];

            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         bilinterms[i] = bilinterm;
         bilinconss[i] = bilincons;
         bilinposs[i] = bilinpos;
         perm[i] = i;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &perm);

   return SCIP_OKAY;
}

/** stores all bilinear terms in the quadratic constraint handler data; in addition, for each bilinear term we store
 *  the number of nonconvex constraints that require to over- or underestimate this term, which only depends on the
 *  lhs, rhs, and the bilinear coefficient
 */
static
SCIP_RETCODE storeAllBilinearTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_BILINTERM* bilinterms;
   SCIP_CONS** bilincons;
   int* bilinpos;
   int nbilinterms;
   int pos;
   int c;
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(conss != NULL);

   /* check for all cases for which we don't want to spend time for collecting all bilinear terms */
   if( nconss == 0 || conshdlrdata->storedbilinearterms || SCIPgetSubscipDepth(scip) != 0 || SCIPgetDepth(scip) >= 1
      || SCIPinProbing(scip) || SCIPinDive(scip) )
      return SCIP_OKAY;

   assert(conshdlrdata->bilinestimators == NULL);
   assert(conshdlrdata->nbilinterms == 0);

   conshdlrdata->storedbilinearterms = TRUE;
   nbilinterms = 0;

   /* count the number of bilinear terms (including duplicates) */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      nbilinterms += consdata->nbilinterms;
   }

   /* no bilinear terms available -> stop */
   if( nbilinterms == 0 )
      return SCIP_OKAY;

   /* allocate temporary memory for sorting all bilinear terms (including duplicates) */
   SCIP_CALL( SCIPallocBufferArray(scip, &bilinterms, nbilinterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bilincons, nbilinterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bilinpos, nbilinterms) );

   /* copy all bilinear terms; note that we need separate entries for x*y and y*x */
   pos = 0;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);

      /* allocate memory to store the later computed indices of each bilinear term in the bilinterms array of the
       * constraint handler data
       */
      if( consdata->nbilinterms > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->bilintermsidx, consdata->nbilinterms) );
      }

      for( i = 0; i < consdata->nbilinterms; ++i )
      {
         assert(consdata->bilinterms != NULL);
         assert(consdata->bilinterms[i].var1 != consdata->bilinterms[i].var2);

         /* add xy */
         bilinterms[pos] = consdata->bilinterms[i];
         bilincons[pos] = conss[c];
         bilinpos[pos] = i;
         ++pos;

         /* invalidate bilinear term index */
         assert(consdata->bilintermsidx != NULL);
         consdata->bilintermsidx[i] = -1;
      }
   }
   assert(pos == nbilinterms);

   /* sorts all bilinear terms (including duplicates) */
   SCIP_CALL( sortAllBilinTerms(scip, bilinterms, nbilinterms, bilincons, bilinpos) );

   /* count the number of bilinear terms without duplicates */
   conshdlrdata->nbilinterms = nbilinterms;
   for( i = 0; i < nbilinterms - 1; ++i )
   {
      assert(bilinTermCompVolume((void*)bilinterms, i, i+1) != 0 || bilinTermComp2((void*)bilinterms, i, i+1) <= 0);

      if( bilinTermComp2((void*)bilinterms, i, i+1) == 0 )
         --(conshdlrdata->nbilinterms);
   }
   assert(conshdlrdata->nbilinterms <= nbilinterms && conshdlrdata->nbilinterms > 0);

   /* store all information for each bilinear term into the constraint handler data */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &conshdlrdata->bilinestimators, conshdlrdata->nbilinterms) );

   /* filter duplicates and update entries in the corresponding constraint datas */
   pos = 0;
   for( i = 0; i < nbilinterms; ++i )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(bilincons[i]);
      SCIP_VAR* x;
      SCIP_Bool haslhs = !SCIPisInfinity(scip, -consdata->lhs);
      SCIP_Bool hasrhs = !SCIPisInfinity(scip, consdata->rhs);

      assert(consdata != NULL);
      assert(bilinpos[i] >= 0 && bilinpos[i] < consdata->nbilinterms);

      /* check for a new bilinear term */
      if( i == 0 || bilinTermComp2((void*)bilinterms, i-1, i) != 0 )
      {
         conshdlrdata->bilinestimators[pos].x = bilinterms[i].var1;
         conshdlrdata->bilinestimators[pos].y = bilinterms[i].var2;
         conshdlrdata->bilinestimators[pos].lastimprfac = 0.0;
         conshdlrdata->bilinestimators[pos].maxnonconvexity = 0.0;
         ++pos;
      }

      /* store whether under- or overestimation is needed for each bilinear term; note that we do not consider convex
       * constraints because they will not be used in separated generateCutNonConvex(), which is the only function that
       * uses a term-wise relaxation
       */
      if( SCIPisPositive(scip, bilinterms[i].coef) )
      {
         conshdlrdata->bilinestimators[pos-1].nunderest += (hasrhs && !consdata->isconvex) ? 1 : 0;
         conshdlrdata->bilinestimators[pos-1].noverest += (haslhs && !consdata->isconcave) ? 1 : 0;
         conshdlrdata->bilinestimators[pos-1].maxnonconvexity = MAX(conshdlrdata->bilinestimators[pos-1].maxnonconvexity, consdata->maxnonconvexity);
      }
      else
      {
         assert(SCIPisNegative(scip, bilinterms[i].coef));
         conshdlrdata->bilinestimators[pos-1].nunderest += (haslhs && !consdata->isconcave) ? 1 : 0;
         conshdlrdata->bilinestimators[pos-1].noverest += (hasrhs && !consdata->isconvex) ? 1 : 0;
         conshdlrdata->bilinestimators[pos-1].maxnonconvexity = MAX(conshdlrdata->bilinestimators[pos-1].maxnonconvexity, consdata->maxnonconvexity);
      }

      /* update index of bilinear term in the constraint data */
      x = consdata->bilinterms[bilinpos[i]].var1;

      assert(pos > 0);
      if( x == conshdlrdata->bilinestimators[pos-1].x )
      {
         assert(consdata->bilinterms[bilinpos[i]].var2 == conshdlrdata->bilinestimators[pos-1].y);
         consdata->bilintermsidx[bilinpos[i]] = pos-1;
      }
   }
   assert(pos == conshdlrdata->nbilinterms);

#ifndef NDEBUG
   /* check whether
    * - all bilintermsidx entries have been set
    * - variables in bilinear terms of each constraint data and the constraint handler data match
    */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      for( i = 0; i < consdata->nbilinterms; ++i )
      {
         SCIP_VAR* x = consdata->bilinterms[i].var1;
         SCIP_VAR* y = consdata->bilinterms[i].var2;
         int idx = consdata->bilintermsidx[i];

         assert(idx >= 0 && idx < conshdlrdata->nbilinterms);
         assert(x == conshdlrdata->bilinestimators[idx].x);
         assert(y == conshdlrdata->bilinestimators[idx].y);

         /* at least one direction is important if the constraint is not convex */
         if( !SCIPisInfinity(scip, consdata->rhs) && !consdata->isconvex )
            assert(conshdlrdata->bilinestimators[idx].nunderest + conshdlrdata->bilinestimators[idx].noverest > 0);
         if( !SCIPisInfinity(scip, -consdata->lhs) && !consdata->isconcave )
            assert(conshdlrdata->bilinestimators[idx].nunderest + conshdlrdata->bilinestimators[idx].noverest > 0);
      }
   }
#endif

   /* free memory */
   SCIPfreeBufferArray(scip, &bilinpos);
   SCIPfreeBufferArray(scip, &bilincons);
   SCIPfreeBufferArray(scip, &bilinterms);

   return SCIP_OKAY;
}

/** frees memory allocated in storeAllBilinearTerms() */
static
SCIP_RETCODE freeAllBilinearTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss              /**< number of constraints */

   )
{
   int c;

   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[c]); /*lint !e613*/
      assert(consdata != NULL);

      SCIPfreeBlockMemoryArrayNull(scip, &consdata->bilintermsidx, consdata->nbilinterms);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->bilinestimators, conshdlrdata->nbilinterms);

   conshdlrdata->nbilinterms = 0;
   conshdlrdata->storedbilinearterms = FALSE;

   return SCIP_OKAY;
}

/** tries to compute cut for multleft * <coefleft, x'> * multright <= rhs / (multright * <coefright, x'>) where x'=(x,1) */
static
SCIP_RETCODE generateCutFactorableDo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_Real             multleft,           /**< multiplicator on lhs */
   SCIP_Real*            coefleft,           /**< coefficient for factor on lhs */
   SCIP_Real             multright,          /**< multiplicator on both sides */
   SCIP_Real*            coefright,          /**< coefficient for factor that goes to rhs */
   SCIP_Real             rightminactivity,   /**< minimal activity of <coefright, x> */
   SCIP_Real             rightmaxactivity,   /**< maximal activity of <coefright, x> */
   SCIP_Real             rhs,                /**< denominator on rhs */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store cut coefs and constant */
   SCIP_Bool*            success             /**< buffer to indicate whether a cut was successfully computed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real constant;
   SCIP_Real coef;
   int i;

   assert(rowprep != NULL);
   assert(rightminactivity * multright > 0.0);
   assert(rightmaxactivity * multright > 0.0);
   assert(multright == 1.0 || multright == -1.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   rowprep->sidetype = SCIP_SIDETYPE_RIGHT;

   if( rhs > 0.0 )
   {
      /* if rhs > 0.0, then rhs / (multright * <coefright, x'>) is convex, thus need secant:
       *  1 / multright*<coefright, x'> <= 1/minact + 1/maxact - 1/(minact * maxact) multright*<coefright, x'>
       *  where [minact, maxact] = multright * [rightminactivity, rightmaxactivity]
       *
       * assuming multright is either -1 or 1, and substituting gives
       *  multright/rightminactivity + multright/rightmaxactivity  - multright/(rightminactivity * rightmaxactivity) *<coefright, x'>
       *
       * multiplying by rhs, gives the estimate
       *  rhs / (multright * <coefright, x'>) <= rhs * multright * (1/rightminactivity + 1/rightmaxactivity - 1/(rightminactivity * rightmaxactivity) * <coefright, x'>)
       */

      /* cannot do if unbounded */
      if( SCIPisInfinity(scip, rightmaxactivity * multright) )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      assert(SCIPisFeasLE(scip, rightminactivity, rightmaxactivity));

      constant = multleft * multright * coefleft[consdata->nquadvars];
      constant -= rhs * multright * (1.0 / rightminactivity + 1.0 / rightmaxactivity);
      constant += rhs * multright * coefright[consdata->nquadvars] / (rightminactivity * rightmaxactivity);

      SCIPaddRowprepConstant(rowprep, constant);

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         coef = multleft * multright * coefleft[i];
         coef += rhs * multright / (rightminactivity * rightmaxactivity) * coefright[i];
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, consdata->quadvarterms[i].var, coef) );
      }

      (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%s_factorablesecant_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
   }
   else
   {
      SCIP_Real refvalue;

      /* if rhs < 0.0, then rhs / (multright * <coefright, x'>) is convex, thus need linearization:
       *     rhs / (multright * <coefright, x'>)
       *  <= rhs / (multright * <coefright, ref'>) - rhs / (multright * <coefright, ref'>)^2 * (multright * <coefright, x'> - multright * <coefright, ref'>)
       *   = 2*rhs / (multright * <coefright, ref'>) - rhs / (multright * <coefright, ref'>)^2 * (multright * <coefright, x'>)
       *
       *  where ref' = (ref, 1)
       */

      /* compute <coefright, ref'> */
      refvalue = coefright[consdata->nquadvars];
      for( i = 0; i < consdata->nquadvars; ++i )
         refvalue += coefright[i] * ref[i];

      /* should not happen, since we checked activity of <coefright,x> before, and assume ref within bounds */
      assert(!SCIPisZero(scip, refvalue));

      constant  = multleft * multright * coefleft[consdata->nquadvars];
      constant -= 2.0 * rhs / (multright * refvalue);
      constant += rhs / (refvalue * refvalue) * multright * coefright[consdata->nquadvars];

      SCIPaddRowprepConstant(rowprep, constant);

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         coef = multleft * multright * coefleft[i];
         coef += rhs / (refvalue * refvalue) * multright * coefright[i];
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, consdata->quadvarterms[i].var, coef) );
      }

      (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%s_factorablelinearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));
   }

   /* @todo does not always need to be local */
   rowprep->local = TRUE;
   *success = TRUE;

   return SCIP_OKAY;
}

/** tries to generate a cut if constraint quadratic function is factorable and there are no linear variables
 * (ax+b)(cx+d) <= rhs and cx+d >= 0 -> (ax+b) <= rhs / (cx+d), where the right hand side is concave and can be linearized
 */
static
SCIP_RETCODE generateCutFactorable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_ROWPREP*         rowprep,            /**< data structure to store cut coefficients */
   SCIP_Bool*            success             /**< buffer to indicate whether a cut was successfully computed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real leftminactivity;
   SCIP_Real leftmaxactivity;
   SCIP_Real rightminactivity;
   SCIP_Real rightmaxactivity;
   SCIP_Real multleft;
   SCIP_Real multright;
   SCIP_Real rhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nlinvars == 0);
   assert(consdata->factorleft != NULL);
   assert(consdata->factorright != NULL);

   *success = FALSE;

   leftminactivity = consdata->factorleft[consdata->nquadvars];
   leftmaxactivity = consdata->factorleft[consdata->nquadvars];
   rightminactivity = consdata->factorright[consdata->nquadvars];
   rightmaxactivity = consdata->factorright[consdata->nquadvars];
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( !SCIPisInfinity(scip, -leftminactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, leftmaxactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }

      if( !SCIPisInfinity(scip, -rightminactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, rightmaxactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }
   }

   /* write violated constraints as multleft * factorleft * factorright <= rhs */
   if( violside == SCIP_SIDETYPE_RIGHT )
   {
      rhs = consdata->rhs;
      multleft = 1.0;
   }
   else
   {
      rhs = -consdata->lhs;
      multleft = -1.0;
   }

   if( SCIPisZero(scip, rhs) )
   {
      /* @todo do something for rhs == 0.0? */
      return SCIP_OKAY;
   }

   if( !SCIPisFeasPositive(scip, leftminactivity) && !SCIPisFeasNegative(scip, leftmaxactivity) )
   {
      /* left factor has 0 within activity bounds, or is very close, at least */
      if( !SCIPisFeasPositive(scip, rightminactivity) && !SCIPisFeasNegative(scip, rightmaxactivity) )
      {
         /* right factor also has 0 within activity bounds, or is very close, at least
          * -> cannot separate
          */
         return SCIP_OKAY;
      }

      /* write violated constraint as multleft * factorleft * multright * (multright * factorright) <= rhs
       *   such that multright * factorright > 0.0
       */
      if( rightminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorleft * multright <= rhs / (factorright * multright) */
      SCIP_CALL( generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorleft, multright, consdata->factorright, rightminactivity, rightmaxactivity, rhs, rowprep, success) );
   }
   else if( !SCIPisFeasPositive(scip, rightminactivity) && !SCIPisFeasNegative(scip, rightmaxactivity) )
   {
      /* left factor is bounded away from 0
       * right factor has 0 within activity bounds, or is very close, at least
       * -> so divide by left factor
       */

      /* write violated constraint as multleft * factorright * multright * (multright * factorleft) <= rhs
       *   such that multright * factorleft > 0.0
       */
      if( leftminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorright * multright <= rhs / (factorleft * multright) */
      SCIP_CALL( generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorright, multright, consdata->factorleft, leftminactivity, leftmaxactivity, rhs, rowprep, success) );
   }
   else if( SCIPisInfinity(scip, -leftminactivity) || SCIPisInfinity(scip, leftmaxactivity) ||
      (!SCIPisInfinity(scip, -rightminactivity) && !SCIPisInfinity(scip, rightmaxactivity) && rightmaxactivity - rightminactivity < leftmaxactivity - leftminactivity) )
   {
      /* both factors are bounded away from 0, but the right one has a smaller activity range, so divide by that one */

      /* write violated constraint as multleft * factorleft * multright * (multright * factorright) <= rhs
       *   such that multright * factorright > 0.0
       */
      if( rightminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorleft * multright <= rhs / (factorright * multright) */
      SCIP_CALL( generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorleft, multright, consdata->factorright, rightminactivity, rightmaxactivity, rhs, rowprep, success) );
   }
   else
   {
      /* both factors are bounded away from 0, but the left one has a smaller activity range, so divide by that one */

      /* write violated constraint as multleft * factorright * multright * (multright * factorleft) <= rhs
       *   such that multright * factorleft > 0.0
       */
      if( leftminactivity < 0.0 )
         multright = -1.0;
      else
         multright =  1.0;

      /* generate cut for multleft * factorright * multright <= rhs / (factorleft * multright) */
      SCIP_CALL( generateCutFactorableDo(scip, cons, ref, multleft, consdata->factorright, multright, consdata->factorleft, leftminactivity, leftmaxactivity, rhs, rowprep, success) );
   }

   return SCIP_OKAY;
}

/* finds intersections of a parametric line (x,y) = (x0,y0) + t [(x1,y1) - (x0,y0)] on curves x*y = wl and x*y = wu;
 * returns TRUE if unsuccessful and FALSE otherwise
 */
static
SCIP_Bool generateCutLTIfindIntersection(
   SCIP*                 scip,
   SCIP_Real             x0,
   SCIP_Real             y0_,
   SCIP_Real             x1,
   SCIP_Real             y1_,
   SCIP_Real             wl,
   SCIP_Real             wu,
   SCIP_Real*            xl,
   SCIP_Real*            yl,
   SCIP_Real*            xu,
   SCIP_Real*            yu
   )
{
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real c;
   SCIP_Real tl;
   SCIP_Real tu;

   assert(wl == SCIP_INVALID || (xl != NULL && yl != NULL));  /*lint !e777 */
   assert(wu == SCIP_INVALID || (xu != NULL && yu != NULL));  /*lint !e777 */

   /* The parametric line is of the form
    *
    *  x = x0 + t (x1-x0)
    *  y = y0 + t (y1-y0)
    *
    * and for that to satisfy xy = wl and xy = wu we must have
    *
    * x0 y0 + t [x0 (y1-y0) + y0 (x1-x0)] + t^2 (x1-x0) (y1-y0) = wl
    *                                                           = wu
    *
    * or a t^2 + b t + c - wl = 0 for proper values of a,b,c.
    *    a t^2 + b t + c - wu = 0
    *
    * Because of the way this procedure will be used, one of the two
    * solutions found we must always use the minimum nonnegative one
    */

   a = (x1 - x0) * (y1_ - y0_);
   c = x0 * y0_;
   b = x0 * y1_ + y0_ * x1 - 2.0 * c;

   tl = 0.0;
   tu = 0.0;

   if( !SCIPisZero(scip, (SCIP_Real)a) )
   {
      if( wl != SCIP_INVALID )  /*lint !e777 */
      {
         SCIP_Real tl1;
         SCIP_Real tl2;
         SCIP_Real denom;
         SCIP_Real q;

         if( b * b - 4.0 * a * (c - wl) < 0.0 )
         {
            SCIPdebugMsg(scip, "probable numerical difficulties, give up\n");
            return TRUE;
         }

         denom = sqrt(b * b - 4.0 * a * (c - wl));
         q = -0.5 * (b + COPYSIGN(denom, b));
         tl1 = q / a;
         tl2 = (c - wl) / q;

         /* choose the smallest non-negative root */
         tl = (tl1 >= 0.0 && (tl2 < 0.0 || tl1 < tl2)) ? tl1 : tl2;
      }

      if( wu != SCIP_INVALID )  /*lint !e777 */
      {
         SCIP_Real tu1;
         SCIP_Real tu2;
         SCIP_Real denom;
         SCIP_Real q;

         if( b * b - 4.0 * a * (c - wu) < 0.0 )
         {
            SCIPdebugMsg(scip, "probable numerical difficulties, give up\n");
            return TRUE;
         }

         denom = sqrt(b * b - 4.0 * a * (c - wu));
         q = -0.5 * (b + COPYSIGN(denom, b));
         tu1 = q / a;
         tu2 = (c - wu) / q;

         /* choose the smallest non-negative root */
         tu = (tu1 >= 0.0 && (tu2 < 0.0 || tu1 < tu2)) ? tu1 : tu2;
      }
   }
   else if( !SCIPisZero(scip, (SCIP_Real)b) )
   {
      if( wl != SCIP_INVALID )  /*lint !e777 */
         tl = (wl - c) / b;
      if( wu != SCIP_INVALID )  /*lint !e777 */
         tu = (wu - c) / b;
   }
   else
   {
      /* no or infinitely many solutions */
      return TRUE;
   }

   if( wl != SCIP_INVALID )  /*lint !e777 */
   {
      *xl = (SCIP_Real)(x0  + tl * (x1  - x0 ));
      *yl = (SCIP_Real)(y0_ + tl * (y1_ - y0_));

      if( !SCIPisRelEQ(scip, *xl * *yl, wl) )
      {
         SCIPdebugMsg(scip, "probable numerical difficulties, give up\n");
         return TRUE;
      }
   }

   if( wu != SCIP_INVALID )  /*lint !e777 */
   {
      *xu = (SCIP_Real)(x0  + tu * (x1 -  x0));
      *yu = (SCIP_Real)(y0_ + tu * (y1_ - y0_));

      if( !SCIPisRelEQ(scip, *xu * *yu, wu) )
      {
         SCIPdebugMsg(scip, "probable numerical difficulties, give up\n");
         return TRUE;
      }
   }

   /* do not use the computed points if one of the components is infinite */
   if( (xu != NULL && SCIPisInfinity(scip, *xu)) || (xl != NULL && SCIPisInfinity(scip, -*xl)) ||
      (yu != NULL && SCIPisInfinity(scip, *yu)) || (yl != NULL && SCIPisInfinity(scip, -*yl)) )
   {
      SCIPdebugMsg(scip, "probable numerical difficulties, give up\n");
      return TRUE;
   }

   return FALSE;
}

/** generate coefficients for a plane through points (x1, y1_, x1*y1) and (x2, y2, x2*y2)
 *  such that intersecting it with one of them (the first if whichuse is FALSE, the second otherwise)
 *  gives a tangent to the curve x*y = k
 *
 *  Returns TRUE on error and FALSE on success.
 */
static
SCIP_Bool generateCutLTIgenMulCoeff(
   SCIP*                 scip,
   SCIP_Real             x1,
   SCIP_Real             y1_,
   SCIP_Real             x2,
   SCIP_Real             y2,
   SCIP_Bool             whichuse,
   SCIP_Real*            cx,
   SCIP_Real*            cy,
   SCIP_Real*            cw
   )
{
   SCIP_Real xd;
   SCIP_Real yd;
   SCIP_Real xo;
   SCIP_Real yo;

   assert(cx != NULL);
   assert(cy != NULL);
   assert(cw != NULL);

   /* the x-y slope of this constraint must be tangent to a curve x*y = k at (xD,yD) */
   if( !whichuse )
   {
      xd = x1;
      xo = x2;
      yd = y1_;
      yo = y2;
   }
   else
   {
      xd = x2;
      xo = x1;
      yd = y2;
      yo = y1_;
   }

   *cx = yd;
   *cy = xd;

   /* lift it so that it touches the other curve */

   /* if the two points are on the same curve, then no cut */
   if( SCIPisZero(scip, xo * yo - xd * yd) )
      return TRUE;

   /* should ALWAYS be negative */
   *cw = (2.0 * xd * yd - (*cx * xo + *cy * yo)) / (xo * yo - xd * yd);

   return FALSE;
}

/** computes coefficients of a lifted-tangent inequality for x*y = w
 *
 *  The code is an adaptation of the methods in exprMul-upperHull.cpp in Couenne/stable/0.4 rev773,
 *  written by P. Belotti and licensed under Eclipse Public License.
 */
static
void generateCutLTIcomputeCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             xl,                 /**< lower bound on x */
   SCIP_Real             xu,                 /**< upper bound on x */
   SCIP_Real             x0,                 /**< reference point for x */
   SCIP_Real             yl,                 /**< lower bound on y */
   SCIP_Real             yu,                 /**< upper bound on y */
   SCIP_Real             y0_,                /**< reference point for y */
   SCIP_Real             wl,                 /**< lower bound on w */
   SCIP_Real             wu,                 /**< upper bound on w */
   SCIP_Real             w0,                 /**< reference point for w */
   SCIP_Real*            cx,                 /**< buffer where to store cut coefficient for x */
   SCIP_Real*            cy,                 /**< buffer where to store cut coefficient for y */
   SCIP_Real*            cw,                 /**< buffer where to store cut coefficient for w */
   SCIP_Real*            c0,                 /**< buffer where to store cut left-hand-side */
   SCIP_Bool*            success             /**< buffer where to indicate whether cut coefficients were computed */
   )
{
   SCIP_Bool flipx;
   SCIP_Bool flipy;
   SCIP_Bool flipw;
   SCIP_Real tmp;
   SCIP_Real xlow;
   SCIP_Real ylow;
   SCIP_Real xupp;
   SCIP_Real yupp;
   SCIP_Real c0x;
   SCIP_Real c0y;
   SCIP_Real c0w;

   assert(scip != NULL);
   assert(cx != NULL);
   assert(cy != NULL);
   assert(cw != NULL);
   assert(c0 != NULL);
   assert(success != NULL);

   *success = FALSE;
   *cx = 0.0;
   *cy = 0.0;
   *cw = 0.0;
   *c0 = 0.0;

   SCIPdebugMsg(scip, "entering points:\n");
   SCIPdebugMsg(scip, "x: %9g\t[%9g\t%9g]\n", x0, xl, xu);
   SCIPdebugMsg(scip, "y: %9g\t[%9g\t%9g]\n", y0_, yl, yu);
   SCIPdebugMsg(scip, "w: %9g\t[%9g\t%9g]\n", w0, wl, wu);

   /* generateCutLTI should have recognized these */
   assert(wl >= 0.0 || wu <= 0.0);
   assert(!SCIPisInfinity(scip, -wl));
   assert(!SCIPisInfinity(scip,  wu));

   assert(SCIPisFeasGE(scip, x0, xl));
   assert(SCIPisFeasLE(scip, x0, xu));
   assert(SCIPisFeasGE(scip, y0_, yl));
   assert(SCIPisFeasLE(scip, y0_, yu));

   /* preliminary bound tightening */
   if( wl >= 0.0 )
   {
      if( xl >= 0.0 || yl >= 0.0 || SCIPisLT(scip, xl * yl, wl) )
      {
         xl = MAX(xl, 0.0);
         yl = MAX(yl, 0.0);
      }
      else if( xu <= 0.0 || yu <= 0.0 || SCIPisLT(scip, xu * yu, wl) )
      {
         xu = MIN(xu, 0.0);
         yu = MIN(yu, 0.0);
      }
      else
      {
         /* both variables have mixed sign (xl < 0 && xu > 0 && yl < 0 && yu > 0) and both xl*yl and xu*yu are feasible
          * cannot generate cut for this
          */
         return;
      }
   }
   else
   {
      if( xl >= 0.0 || yu <= 0.0 || SCIPisGT(scip, xl * yu, wu) )
      {
         xl = MAX(xl, 0.0);
         yu = MIN(yu, 0.0);
      }
      else if( xu <= 0.0 || yl >= 0.0 || SCIPisGT(scip, xu * yl, wu))
      {
         xu = MIN(xu, 0.0);
         yl = MAX(yl, 0.0);
      }
      else
      {
         /* both variables have mixed sign (xl < 0 && xu > 0 && yl < 0 && yu > 0) and both xl*yu and xu*yl are feasible
          * cannot generate cut for this
          */
         return;
      }
   }

   /* if x or y is fixed now or even infeasible, then do not think about a cut */
   if( SCIPisGE(scip, xl, xu) || SCIPisGE(scip, yl, yu) )
      return;

   /* reduce to positive orthant by flipping variables */
   if( xl < 0.0 )
   {
      flipx = TRUE;
      tmp = xu;
      xu = -xl;
      xl = -tmp;
      x0 = -x0;
   }
   else
      flipx = FALSE;

   if( yl < 0.0 )
   {
      flipy = TRUE;
      tmp = yu;
      yu = -yl;
      yl = -tmp;
      y0_ = -y0_;
   }
   else
      flipy = FALSE;

   if( flipx ^ flipy )
   {
      flipw = TRUE;
      tmp = wu;
      wu = -wl;
      wl = -tmp;
      w0 = -w0;
   }
   else
      flipw = FALSE;

   /* project refpoint into box not only for numerical reasons, but also due to preliminary bound tightening above */
   x0 = MIN(xu, MAX(x0, xl));
   y0_ = MIN(yu, MAX(y0_, yl));
   w0 = MIN(wu, MAX(w0, wl));

   SCIPdebugMsg(scip, "reduced points:\n");
   SCIPdebugMsg(scip, "x: %9g\t[%9g\t%9g]\n", x0, xl, xu);
   SCIPdebugMsg(scip, "y: %9g\t[%9g\t%9g]\n", y0_, yl, yu);
   SCIPdebugMsg(scip, "w: %9g\t[%9g\t%9g]\n", w0, wl, wu);

   if( SCIPisGE(scip, xl * yl, wl) && SCIPisLE(scip, xu * yu, wu) )
   {
      SCIPdebugMsg(scip, "box for x and y inside feasible region -> nothing to separate\n");
      return;
   }
   if( SCIPisGE(scip, x0 * y0_, w0) )
   {
      SCIPdebugMsg(scip, "point to separate not below curve -> cannot separate\n");
      return;
   }

   /* find intersections of halfline from origin
    * return if no proper point could be found
    */
   if( generateCutLTIfindIntersection(scip, 0.0, 0.0, x0, y0_, wl, wu, &xlow, &ylow, &xupp, &yupp) )
      return;

   SCIPdebugMsg(scip, "intersections:\n");
   SCIPdebugMsg(scip, "lower: %9g\t%9g\tprod %9g\n", xlow, ylow, xlow*ylow);
   SCIPdebugMsg(scip, "upper: %9g\t%9g\tprod %9g\n", xupp, yupp, xupp*yupp);

   /* Case 1: If both are outside of bounding box, either NW or SE, then McCormick is sufficient, so return */
   if( (xlow <= xl && yupp >= yu) || (ylow <= yl && xupp >= xu) )
      return;

   /* There will be at least one cut. Define coefficients and rhs ---will have to change them back if (flipX || flipY) */
   if( xlow >= xl && xupp <= xu && ylow >= yl && yupp <= yu )
   {
      /* Case 2: both are inside. Easy lifting... */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, FALSE, cx, cy, cw) )
         return;

      c0x = *cx * xlow;
      c0y = *cy * ylow;
      c0w = *cw * wl;
   }
   else if( xlow >= xl && ylow >= yl && (xupp > xu || yupp > yu) )
   {
      /* Case 3a and 3b: through lower curve, but not upper. */
      if( yupp > yu )
      {
         /* upper intersect is North; place it within box */
         assert(!SCIPisInfinity(scip, yu));
         yupp = yu;
         xupp = wu / yu;
      }
      else
      {
         /* upper intersect is East; place it within box */
         assert(!SCIPisInfinity(scip, xu));
         xupp = xu;
         yupp = wu / xu;
      }

      /* find intersection on low curve on half line through new point and (x0,y0_) */
      if( generateCutLTIfindIntersection(scip, xupp, yupp, x0, y0_, wl, SCIP_INVALID, &xlow, &ylow, NULL, NULL) )
         return;

      /* check whether McCormick is sufficient */
      if( xlow < xl || ylow < yl )
         return;

      /* lift inequality on lower point */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, FALSE, cx, cy, cw) )
         return;

      c0x = *cx * xlow;
      c0y = *cy * ylow;
      c0w = *cw * wl;
   }
   else if( xupp <= xu && yupp <= yu && (xlow < xl || ylow < yl) )
   {
      /* Case 4a and 4b: viceversa (lift for validity) */
      if( ylow < yl )
      {
         /* upper intersect is South; place it within box */
         assert(!SCIPisZero(scip, yl));
         ylow = yl;
         xlow = wl / yl;
      }
      else
      {
         /* upper intersect is West; place it within box */
         assert(!SCIPisZero(scip, xl));
         xlow = xl;
         ylow = wl / xl;
      }

      /* find intersection on low curve on half line through new point and (x0,y0) */
      if( generateCutLTIfindIntersection(scip, xlow, ylow, x0, y0_, SCIP_INVALID, wu, NULL, NULL, &xupp, &yupp) )
         return;

      /* check whether McCormick is sufficient */
      if( xupp > xu || yupp > yu )
         return;

      /* lift inequality on UPPER point */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, TRUE, cx, cy, cw) )
         return;

      c0x = *cx * xupp;
      c0y = *cy * yupp;
      c0w = *cw * wu;
   }
   else if( (xlow < xl && xupp > xu) || (ylow < yl && yupp > yu) )
   {
      /* Case 5: both outside of bounding box, N and S or W and E. */
#if 0
      SCIP_Real xlow2;
      SCIP_Real ylow2;
      SCIP_Real xupp2;
      SCIP_Real yupp2;
#endif

      if( ylow < yl )
      {
         /* upper intersect is South; place it within box */
         assert(!SCIPisZero(scip, yl));
         assert(!SCIPisZero(scip, yu));
         ylow = yl;
         yupp = yu;
         xlow = wl / yl;
         xupp = wu / yu;
      }
      else
      {
         /* upper intersect is West; place it within box */
         assert(!SCIPisZero(scip, xl));
         assert(!SCIPisZero(scip, xu));
         xlow = xl;
         xupp = xu;
         ylow = wl / xl;
         yupp = wu / xu;
      }

      SCIPdebugMsg(scip, "New intersections:\n");
      SCIPdebugMsg(scip, "lower: %9g\t%9g\tprod %9g\n", xlow, ylow, xlow*ylow);
      SCIPdebugMsg(scip, "upper: %9g\t%9g\tprod %9g\n", xupp, yupp, xupp*yupp);

#if 1
      /* Nothing to find. Just separate two inequalities at the same point, just using different support */
      if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, FALSE, cx, cy, cw) )
      {
         if( generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp, yupp, TRUE, cx, cy, cw) )
            return;

         c0x = *cx * xupp;
         c0y = *cy * yupp;
         c0w = *cw * wu;
      }
      else
      {
         c0x = *cx * xlow;
         c0y = *cy * ylow;
         c0w = *cw * wl;
      }

#else
      /* find the intersection on the lower (upper) curve on the line through xLP and the upper (lower) point
       * this does not seem to work (cuts off solution at nous2), so it is disabled for now
       */
      if( generateCutLTIfindIntersection(scip, xlow, ylow, x0, y0_, SCIP_INVALID, wu, NULL, NULL, &xupp2, &yupp2) ||
         generateCutLTIgenMulCoeff(scip, xlow, ylow, xupp2, yupp2, FALSE, cx, cx, cw) )
      {
         if( generateCutLTIfindIntersection(scip, xupp, yupp, x0, y0_, wl, SCIP_INVALID, &xlow2, &ylow2, NULL, NULL) ||
            generateCutLTIgenMulCoeff(scip, xlow2, ylow2, xupp, yupp, TRUE, cx, cy, cw) )
            return;

         c0x = *cx * xupp;
         c0y = *cy * yupp;
         c0w = *cw * wu;
      }
      else
      {
         c0x = *cx * xlow;
         c0y = *cy * ylow;
         c0w = *cw * wl;
      }
#endif
   }
   else
   {
      SCIPdebugMsg(scip, "points are in a weird position:\n");
      SCIPdebugMsg(scip, "lower: %9g\t%9g\tprod %9g\n", xlow, ylow, xlow*ylow);
      SCIPdebugMsg(scip, "upper: %9g\t%9g\tprod %9g\n", xupp, yupp, xupp*yupp);

      return;
   }

   SCIPdebugMsg(scip, "cut w.r.t. reduced points: %gx-%g %+gy-%g %+gw-%g >= 0\n",
      *cx, c0x, *cy, c0y, *cw, c0w);

   /* re-transform back into original variables */
   if( flipx )
      *cx = -*cx;
   if( flipy )
      *cy = -*cy;
   if( flipw )
      *cw = -*cw;

   *c0 = c0x + c0y + c0w;

   *success = TRUE;
}

/** tries to generate a cut if constraint quadratic function is factorable and there are linear variables
 *
 *  Computes what is called a lifted tangent inequality described in@n
 *   Belotti, Miller, Namazifar, Lifted inequalities for bounded products of variables, SIAG/OPT Views-and-News 22:1, 2011
 */
static
SCIP_RETCODE generateCutLTI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_SOL*             sol,                /**< solution that shall be cutoff, NULL for LP solution */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store cut data */
   SCIP_Bool*            success             /**< buffer to indicate whether a cut was successfully computed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real leftminactivity;
   SCIP_Real leftmaxactivity;
   SCIP_Real leftrefactivity;
   SCIP_Real rightminactivity;
   SCIP_Real rightmaxactivity;
   SCIP_Real rightrefactivity;
   SCIP_Real rhsminactivity;
   SCIP_Real rhsmaxactivity;
   SCIP_Real rhsrefactivity;
   SCIP_Real coefleft;
   SCIP_Real coefright;
   SCIP_Real coefrhs;
   SCIP_Real cutlhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);
   /* currently only separate LP solution or solutions given as SCIP_SOL, i.e., no cutgeneration during initlp */
   assert(sol != NULL || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nlinvars > 0);
   assert(consdata->factorleft != NULL);
   assert(consdata->factorright != NULL);

   *success = FALSE;
   rowprep->sidetype = SCIP_SIDETYPE_LEFT;

   /* write violated constraints as factorleft * factorright '==' rhs
    * where rhs are constraint sides - activity bound of linear part
    */
   rhsminactivity = consdata->lhs;
   rhsmaxactivity = consdata->rhs;
   rhsrefactivity = (violside == SCIP_SIDETYPE_LEFT ? consdata->lhs : consdata->rhs);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( !SCIPisInfinity(scip, -rhsminactivity) )
      {
         if( consdata->lincoefs[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->linvars[i])) )
               rhsminactivity = -SCIPinfinity(scip);
            else
               rhsminactivity -= consdata->lincoefs[i] * SCIPvarGetLbLocal(consdata->linvars[i]);
         }
         else
         {
            assert(consdata->lincoefs[i] > 0.0);
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->linvars[i])) )
               rhsminactivity = -SCIPinfinity(scip);
            else
               rhsminactivity -= consdata->lincoefs[i] * SCIPvarGetUbLocal(consdata->linvars[i]);
         }
      }
      if( !SCIPisInfinity(scip, rhsmaxactivity) )
      {
         if( consdata->lincoefs[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->linvars[i])) )
               rhsmaxactivity = SCIPinfinity(scip);
            else
               rhsmaxactivity -= consdata->lincoefs[i] * SCIPvarGetUbLocal(consdata->linvars[i]);
         }
         else
         {
            assert(consdata->lincoefs[i] > 0.0);
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->linvars[i])) )
               rhsmaxactivity = SCIPinfinity(scip);
            else
               rhsmaxactivity -= consdata->lincoefs[i] * SCIPvarGetLbLocal(consdata->linvars[i]);
         }
      }
      rhsrefactivity -= consdata->lincoefs[i] * SCIPgetSolVal(scip, sol, consdata->linvars[i]);
   }

   if( SCIPisInfinity(scip, -rhsminactivity) || SCIPisInfinity(scip, rhsmaxactivity) )
   {
      /* if right hand side is unbounded, then cannot do LTI */
      return SCIP_OKAY;
   }

   if( !SCIPisFeasPositive(scip, rhsminactivity) && !SCIPisFeasNegative(scip, rhsmaxactivity) )
   {
      /* if right hand side has 0 inside activity, then cannot do anything
       * if it has 0.0 as min or max activity, then a usual McCormick should be sufficient, too
       */
      return SCIP_OKAY;
   }

   leftminactivity = consdata->factorleft[consdata->nquadvars];
   leftmaxactivity = consdata->factorleft[consdata->nquadvars];
   leftrefactivity = consdata->factorleft[consdata->nquadvars];
   rightminactivity = consdata->factorright[consdata->nquadvars];
   rightmaxactivity = consdata->factorright[consdata->nquadvars];
   rightrefactivity = consdata->factorright[consdata->nquadvars];
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( !SCIPisInfinity(scip, -leftminactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftminactivity = -SCIPinfinity(scip);
            else
               leftminactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, leftmaxactivity) )
      {
         if( consdata->factorleft[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorleft[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               leftmaxactivity = SCIPinfinity(scip);
            else
               leftmaxactivity += consdata->factorleft[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }
      leftrefactivity += consdata->factorleft[i] * ref[i];

      if( !SCIPisInfinity(scip, -rightminactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightminactivity = -SCIPinfinity(scip);
            else
               rightminactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
      }
      if( !SCIPisInfinity(scip, rightmaxactivity) )
      {
         if( consdata->factorright[i] > 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetUbLocal(consdata->quadvarterms[i].var);
         }
         else if( consdata->factorright[i] < 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->quadvarterms[i].var)) )
               rightmaxactivity = SCIPinfinity(scip);
            else
               rightmaxactivity += consdata->factorright[i] * SCIPvarGetLbLocal(consdata->quadvarterms[i].var);
         }
      }
      rightrefactivity += consdata->factorright[i] * ref[i];
   }

   /* if activities exceed "opposite" infinity, huge bounds seem to be involved, for which the below method is not prepared */
   if( SCIPisInfinity(scip, leftminactivity)  || SCIPisInfinity(scip, -leftmaxactivity) ||
       SCIPisInfinity(scip, rightminactivity) || SCIPisInfinity(scip, -rightmaxactivity) )
      return SCIP_OKAY;

   /* if any of the factors is essentially fixed, give up and do usual method (numerically less sensitive, I hope) */
   if( SCIPisRelEQ(scip, leftminactivity, leftmaxactivity) || SCIPisRelEQ(scip, rightminactivity, rightmaxactivity) )
      return SCIP_OKAY;

   /* success can only be expected for separation of violated x*y <= w, assuming x>=0, y>=0
    * @todo we should check this early? */

   /* call Couenne magic */
   generateCutLTIcomputeCoefs(scip,
      leftminactivity, leftmaxactivity, leftrefactivity,
      rightminactivity, rightmaxactivity, rightrefactivity,
      rhsminactivity, rhsmaxactivity, rhsrefactivity,
      &coefleft, &coefright, &coefrhs, &cutlhs,
      success);

   if( !*success )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "LTI for x[%g,%g] * y[%g,%g] = w[%g,%g]: %gx %+gy %+gw >= %g;  feas: %g\n",
      leftminactivity, leftmaxactivity, rightminactivity, rightmaxactivity, rhsminactivity, rhsmaxactivity,
      coefleft, coefright, coefrhs, cutlhs,
      coefleft * leftrefactivity + coefright * rightrefactivity + coefrhs * rhsrefactivity - cutlhs
      );

   if( coefleft * leftrefactivity + coefright * rightrefactivity + coefrhs * rhsrefactivity >= cutlhs )
   {
      SCIPdebugMsg(scip, "does not cutoff point? :-(\n");
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* setup cut coefs for
    *   coefleft * leftfactor + coefright * rightfactor + coefrhs * w >= cutlhs, where conslhs - lincoefs <= w <= consrhs - lincoefs
    */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, consdata->quadvarterms[i].var, coefleft * consdata->factorleft[i] + coefright * consdata->factorright[i]) );
   }
   SCIPaddRowprepConstant(rowprep, coefleft * consdata->factorleft[i] + coefright * consdata->factorright[i]);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, consdata->linvars[i], -coefrhs * consdata->lincoefs[i]) );
   }
   if( coefrhs > 0.0 )
   {
      /* use coefrhs * w <= coefrhs * (consrhs - lincoefs) */
      assert(!SCIPisInfinity(scip, consdata->rhs));
      SCIPaddRowprepConstant(rowprep, coefrhs * consdata->rhs);
   }
   else
   {
      /* use coefrhs * w <= coeflhs * (conslhs - lincoefs) */
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      SCIPaddRowprepConstant(rowprep, coefrhs * consdata->lhs);
   }
   SCIPaddRowprepSide(rowprep, cutlhs);

   rowprep->local = TRUE;

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%s_lti_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   *success = TRUE;

   return SCIP_OKAY;
}

/** computes cut coefficients by linearizing a quadratic function */
static
SCIP_RETCODE generateCutConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< side for which to generate cut */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store cut data */
   SCIP_Bool*            success             /**< buffer to indicate whether a cut was successfully computed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real constant;
   SCIP_Real coef;
   SCIP_Real coef2;
   SCIP_VAR* var;
   int var2pos;
   int j;
   int k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *success = TRUE;

   /* do first-order Taylor for each term */
   for( j = 0; j < consdata->nquadvars && *success; ++j )
   {
      var = consdata->quadvarterms[j].var;

      /* initialize coefficients to linear coefficients of quadratic variables */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, consdata->quadvarterms[j].lincoef) );

      /* add linearization of square term */
      coef = 0.0;
      constant = 0.0;
      SCIPaddSquareLinearization(scip, consdata->quadvarterms[j].sqrcoef, ref[j],
         consdata->quadvarterms[j].nadjbilin == 0 && SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS, &coef, &constant, success);
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef) );
      SCIPaddRowprepConstant(rowprep, constant);

      /* add linearization of bilinear terms that have var as first variable */
      for( k = 0; k < consdata->quadvarterms[j].nadjbilin && *success; ++k )
      {
         bilinterm = &consdata->bilinterms[consdata->quadvarterms[j].adjbilin[k]];
         if( bilinterm->var1 != var )
            continue;
         assert(bilinterm->var2 != var);
         assert(consdata->sepabilinvar2pos != NULL);

         var2pos = consdata->sepabilinvar2pos[consdata->quadvarterms[j].adjbilin[k]];
         assert(var2pos >= 0);
         assert(var2pos < consdata->nquadvars);
         assert(consdata->quadvarterms[var2pos].var == bilinterm->var2);

         coef = 0.0;
         coef2 = 0.0;
         constant = 0.0;
         SCIPaddBilinLinearization(scip, bilinterm->coef, ref[j], ref[var2pos], &coef, &coef2, &constant, success);
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef) );
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, bilinterm->var2, coef2) );
         SCIPaddRowprepConstant(rowprep, constant);
      }
   }

   if( !*success )
   {
      SCIPdebugMsg(scip, "no success in linearization of <%s> in reference point\n", SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   rowprep->sidetype = violside;
   SCIPaddRowprepSide(rowprep, violside == SCIP_SIDETYPE_LEFT ? consdata->lhs : consdata->rhs);

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%s_side%d_linearization_%d", SCIPconsGetName(cons), violside, SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

/** helper function to update the best relaxation for a bilinear term when using valid linear inequalities */
static
void updateBilinearRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR* RESTRICT    x,                  /**< first variable */
   SCIP_VAR* RESTRICT    y,                  /**< second variable */
   SCIP_Real             bilincoef,          /**< coefficient of the bilinear term */
   SCIP_SIDETYPE         violside,           /**< side of quadratic constraint that is violated */
   SCIP_Real             refx,               /**< reference point for the x variable */
   SCIP_Real             refy,               /**< reference point for the y variable */
   SCIP_Real* RESTRICT   ineqs,              /**< coefficients of each linear inequality; stored as triple (xcoef,ycoef,constant) */
   int                   nineqs,             /**< total number of inequalities */
   SCIP_Real             mccormickval,       /**< value of the McCormick relaxation at the reference point */
   SCIP_Real* RESTRICT   bestcoefx,          /**< pointer to update the x coefficient */
   SCIP_Real* RESTRICT   bestcoefy,          /**< pointer to update the y coefficient */
   SCIP_Real* RESTRICT   bestconst,          /**< pointer to update the constant */
   SCIP_Real* RESTRICT   bestval,            /**< value of the best relaxation that have been found so far */
   SCIP_Bool*            success             /**< buffer to store whether we found a better relaxation */
   )
{
   SCIP_Real constshift[2] = {0.0, 0.0};
   SCIP_Real constant;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real lbx;
   SCIP_Real ubx;
   SCIP_Real lby;
   SCIP_Real uby;
   SCIP_Bool update;
   SCIP_Bool overestimate;
   int i;

   assert(x != y);
   assert(!SCIPisZero(scip, bilincoef));
   assert(nineqs >= 0 && nineqs <= 2);
   assert(bestcoefx != NULL);
   assert(bestcoefy != NULL);
   assert(bestconst != NULL);
   assert(bestval != NULL);

   /* no inequalities available */
   if( nineqs == 0 )
      return;
   assert(ineqs != NULL);

   lbx = SCIPvarGetLbLocal(x);
   ubx = SCIPvarGetUbLocal(x);
   lby = SCIPvarGetLbLocal(y);
   uby = SCIPvarGetUbLocal(y);
   overestimate = (violside == SCIP_SIDETYPE_LEFT);

   /* check cases for which we can't compute a tighter relaxation */
   if( SCIPisFeasLE(scip, refx, lbx) || SCIPisFeasGE(scip, refx, ubx)
      || SCIPisFeasLE(scip, refy, lby) || SCIPisFeasGE(scip, refy, uby) )
      return;

   /* due to the feasibility tolerances of the LP and NLP solver, it might possible that the reference point is
    * violating the linear inequalities; to ensure that we compute a valid underestimate, we relax the linear
    * inequality by changing its constant part
    */
   for( i = 0; i < nineqs; ++i )
   {
      constshift[i] = MAX(0.0, ineqs[3*i] * refx - ineqs[3*i+1] * refy - ineqs[3*i+2]);
      SCIPdebugMsg(scip, "constant shift of inequality %d = %.16f\n", constshift[i]);
   }

   /* try to use both inequalities */
   if( nineqs == 2 )
   {
      SCIPcomputeBilinEnvelope2(scip, bilincoef, lbx, ubx, refx, lby, uby, refy, overestimate, ineqs[0], ineqs[1],
         ineqs[2] + constshift[0], ineqs[3], ineqs[4], ineqs[5] + constshift[1], &xcoef, &ycoef, &constant, &update);

      if( update )
      {
         SCIP_Real val = xcoef * refx + ycoef * refy + constant;
         SCIP_Real relimpr = 1.0 - (REALABS(val - bilincoef * refx * refy) + 1e-4) / (REALABS(*bestval - bilincoef * refx * refy) + 1e-4);
         SCIP_Real absimpr = REALABS(val - (*bestval));

         /* update relaxation if possible */
         if( relimpr > 0.05 && absimpr > 1e-3 && ((overestimate && SCIPisRelLT(scip, val, *bestval)) || (!overestimate && SCIPisRelGT(scip, val, *bestval))) )
         {
            *bestcoefx = xcoef;
            *bestcoefy = ycoef;
            *bestconst = constant;
            *bestval = val;
            *success = TRUE;
         }
      }
   }

   /* use inequalities individually */
   for( i = 0; i < nineqs; ++i )
   {
      SCIPcomputeBilinEnvelope1(scip, bilincoef, lbx, ubx, refx, lby, uby, refy, overestimate, ineqs[3*i], ineqs[3*i+1],
         ineqs[3*i+2] + constshift[i], &xcoef, &ycoef, &constant, &update);

      if( update )
      {
         SCIP_Real val = xcoef * refx + ycoef * refy + constant;
         SCIP_Real relimpr = 1.0 - (REALABS(val - bilincoef * refx * refy) + 1e-4) / (REALABS(mccormickval - bilincoef * refx * refy) + 1e-4);
         SCIP_Real absimpr = REALABS(val - (*bestval));

         /* update relaxation if possible */
         if( relimpr > 0.05 && absimpr > 1e-3 && ((overestimate && SCIPisRelLT(scip, val, *bestval)) || (!overestimate && SCIPisRelGT(scip, val, *bestval))) )
         {
            *bestcoefx = xcoef;
            *bestcoefy = ycoef;
            *bestconst = constant;
            *bestval = val;
            *success = TRUE;
         }
      }
   }
}

/* returns the interiority of a reference point w.r.t. given bounds */
static
SCIP_Real getInteriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lbx,                /**< lower bound of the first variable */
   SCIP_Real             ubx,                /**< upper bound of the first variable  */
   SCIP_Real             refx,               /**< reference point of the first variable */
   SCIP_Real             lby,                /**< lower bound of the second variable */
   SCIP_Real             uby,                /**< upper bound of the second variable  */
   SCIP_Real             refy                /**< reference point of the second variable */
   )
{
   SCIP_Real interiorityx;
   SCIP_Real interiorityy;

   interiorityx = MIN(refx-lbx, ubx-refx) / MAX(ubx-lbx, SCIPepsilon(scip)); /*lint !e666*/
   interiorityy = MIN(refy-lby, uby-refy) / MAX(uby-lby, SCIPepsilon(scip)); /*lint !e666*/

   return 2.0*MIN(interiorityx, interiorityy);
}

/** computes cut coefficients for a nonconvex quadratic function */
static
SCIP_RETCODE generateCutNonConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< side for which to generate cut */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to store cut data */
   SCIP_Bool*            success             /**< buffer to indicate whether a cut was successfully computed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real sqrcoef;
   SCIP_Real coef;
   SCIP_Real coef2;
   SCIP_Real constant;
   SCIP_VAR* var;
   int var2pos;
   int j;
   int k;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   rowprep->local = TRUE;
   *success = TRUE;

   /* underestimate (secant, McCormick) or linearize each term separately */
   for( j = 0; j < consdata->nquadvars && *success; ++j )
   {
      var = consdata->quadvarterms[j].var;

      /* initialize coefficients to linear coefficients of quadratic variables */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, consdata->quadvarterms[j].lincoef) );

      sqrcoef = consdata->quadvarterms[j].sqrcoef;
      if( sqrcoef != 0.0 )
      {
         coef = 0.0;
         constant = 0.0;
         if( (violside == SCIP_SIDETYPE_LEFT  && sqrcoef <= 0.0) || (violside == SCIP_SIDETYPE_RIGHT && sqrcoef > 0.0) )
         {
            /* convex -> linearize */
            SCIPaddSquareLinearization(scip, sqrcoef, ref[j], SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS, &coef,
               &constant, success);
         }
         else
         {
            /* not convex -> secant approximation */
            SCIPaddSquareSecant(scip, sqrcoef, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), ref[j], &coef,
               &constant, success);
         }
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef) );
         SCIPaddRowprepConstant(rowprep, constant);
      }

      /* relax each bilinear term */
      for( k = 0; k < consdata->quadvarterms[j].nadjbilin && (*success); ++k )
      {
         SCIP_VAR* x;
         SCIP_VAR* y;
         SCIP_Real refx;
         SCIP_Real refy;
         SCIP_Real lbx;
         SCIP_Real ubx;
         SCIP_Real lby;
         SCIP_Real uby;
         int idx;

         idx = consdata->quadvarterms[j].adjbilin[k];
         bilinterm = &consdata->bilinterms[idx];
         if( bilinterm->var1 != var )
            continue;
         assert(bilinterm->var2 != var);
         assert(consdata->sepabilinvar2pos != NULL);

         var2pos = consdata->sepabilinvar2pos[consdata->quadvarterms[j].adjbilin[k]];
         assert(var2pos >= 0);
         assert(var2pos < consdata->nquadvars);
         assert(consdata->quadvarterms[var2pos].var == bilinterm->var2);

         /* get data of the variables in the bilinear term */
         x = var;
         y = bilinterm->var2;
         refx = ref[j];
         refy = ref[var2pos];
         lbx = SCIPvarGetLbLocal(x);
         ubx = SCIPvarGetUbLocal(x);
         lby = SCIPvarGetLbLocal(y);
         uby = SCIPvarGetUbLocal(y);
         SCIPdebugMsg(scip, "bilinear term %g %s %s with (%g,%g) in [%g,%g]x[%g,%g] overestimate=%u\n", bilinterm->coef,
            SCIPvarGetName(x), SCIPvarGetName(y), refx, refy, lbx, ubx, lby, uby, violside == SCIP_SIDETYPE_LEFT);

         /* use the McCormick relaxation for under- or overestimating the bilinear term */
         coef = 0.0;
         coef2 = 0.0;
         constant = 0.0;
         SCIPaddBilinMcCormick(scip, bilinterm->coef, lbx, ubx, refx, lby, uby, refy,
            violside == SCIP_SIDETYPE_LEFT, &coef, &coef2, &constant, success);
         SCIPdebugMsg(scip, "McCormick = %g (%u)\n", refx * coef + refy * coef2 + constant, *success);

         /* tries to compute a tighter relaxation for xy by using valid linear inequalities */
         if( conshdlrdata->bilinestimators != NULL && ubx - lbx >= 0.1 && uby - lby >= 0.1
            && (SCIPgetNSepaRounds(scip) <= conshdlrdata->bilinineqmaxseparounds || SCIPgetDepth(scip) == 0) )
         {
            BILINESTIMATOR* bilinestimator;
            SCIP_Real mccormick;
            SCIP_Real score;
            int bilintermidx;

            mccormick = refx * coef + refy * coef2 + constant;
            score = getInteriority(scip, lbx, ubx, refx, lby, uby, refy);

            /* get data for bilinear term */
            bilintermidx = consdata->bilintermsidx[idx];
            assert(conshdlrdata->bilinestimators != NULL);
            bilinestimator = &(conshdlrdata->bilinestimators[bilintermidx]);
            assert(bilinestimator->x == x);
            assert(bilinestimator->y == y);

            /* reset the last improvement factor (used for getting better branching decisions) */
            bilinestimator->lastimprfac = 0.0;

            /* compute tighter relaxation for xy if the current score is large enough */
            if( SCIPisGE(scip, score, conshdlrdata->minscorebilinterms)
               && bilinestimator->nineqoverest + bilinestimator->ninequnderest > 0 )
            {
               SCIP_Real bestval = mccormick;
               SCIP_Bool updaterelax = FALSE;

               /*
                * note that we check the sign of the bilinear coefficient together with violside in
                * updateBilinearRelaxation in order to decide whether a valid under- or overestimate can be computed
                */

               /* use overestimates */
               updateBilinearRelaxation(scip, x, y, bilinterm->coef, violside, refx, refy, bilinestimator->ineqoverest,
                  bilinestimator->nineqoverest, mccormick, &coef, &coef2, &constant, &bestval, &updaterelax);

               /* use underestimates */
               updateBilinearRelaxation(scip, x, y, bilinterm->coef, violside, refx, refy, bilinestimator->inequnderest,
                  bilinestimator->ninequnderest, mccormick, &coef, &coef2, &constant, &bestval, &updaterelax);

               SCIPdebugMsg(scip, "found better relaxation value: %u (%g)\n", updaterelax, bestval);

               /* check whether the new relaxation is under- or overestimating xy properly */
               if( updaterelax )
               {
                  /* update improvement factor */
                  bilinestimator->lastimprfac = 1.0 - REALABS(bestval - bilinterm->coef * refx * refy) / REALABS(mccormick - bilinterm->coef * refx * refy);

#ifndef NDEBUG
                  assert(SCIPisEQ(scip, bestval, coef * refx + coef2 * refy + constant));
                  if( violside == SCIP_SIDETYPE_LEFT )
                  {
                     assert(SCIPisRelGE(scip, bestval, bilinterm->coef * refx * refy));
                     assert(SCIPisRelLE(scip, bestval, mccormick));
                  }
                  else
                  {
                     assert(SCIPisRelLE(scip, bestval, bilinterm->coef * refx * refy));
                     assert(SCIPisRelGE(scip, bestval, mccormick));
                  }
#endif
               }
            }
         }

         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef) );
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, bilinterm->var2, coef2) );
         SCIPaddRowprepConstant(rowprep, constant);
      }
   }

   if( !*success )
   {
      SCIPdebugMsg(scip, "no success to find estimator for nonconvex <%s>\n", SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   rowprep->sidetype = violside;
   SCIPaddRowprepSide(rowprep, violside == SCIP_SIDETYPE_LEFT ? consdata->lhs : consdata->rhs);

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%s_side%d_estimation_%d", SCIPconsGetName(cons), violside, SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

/** generates a cut based on linearization (if convex) or McCormick (if nonconvex) in a given reference point */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            ref,                /**< reference solution where to generate the cut */
   SCIP_SOL*             sol,                /**< point that we aim to separate, or NULL for LP solution */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real*            efficacy,           /**< buffer to store efficacy of row in reference solution, or NULL if not of interest */
   SCIP_Bool             checkcurvmultivar,  /**< are we allowed to check the curvature of a multivariate quadratic function, if not done yet */
   SCIP_Real             minefficacy         /**< minimal required efficacy */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_ROWPREP*  rowprep;
   SCIP_Bool      success;
   SCIP_Real      viol = 0.0;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(ref != NULL);
   assert(row != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(violside != SCIP_SIDETYPE_LEFT  || !SCIPisInfinity(scip, -consdata->lhs));
   assert(violside != SCIP_SIDETYPE_RIGHT || !SCIPisInfinity(scip,  consdata->rhs));

   *row = NULL;

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, SCIPconsIsLocal(cons)) );
   success = FALSE;

   /* if constraint function is factorable, then try to use factorable form to generate cut */
   if( consdata->factorleft != NULL )
   {
      if( consdata->nlinvars == 0 )
      {
         SCIP_CALL( generateCutFactorable(scip, cons, violside, ref, rowprep, &success) );
      }
      else if( sol != NULL || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* generateCutLTI needs reference values also for the linear variables, which we only have if sol is given or LP has been solved */
         SCIP_CALL( generateCutLTI(scip, cons, violside, ref, sol, rowprep, &success) );
      }
   }

   /* if constraint is not factorable or failed to generate cut, try default method */
   if( !success )
   {
      SCIP_CALL( checkCurvature(scip, cons, checkcurvmultivar) );

      if( (violside == SCIP_SIDETYPE_LEFT && consdata->isconcave) || (violside == SCIP_SIDETYPE_RIGHT && consdata->isconvex) )
      {
         SCIP_CALL( generateCutConvex(scip, cons, violside, ref, rowprep, &success) );
      }
      else
      {
         SCIP_CALL( generateCutNonConvex(scip, conshdlrdata, cons, violside, ref, rowprep, &success) );
      }

      SCIP_CALL( SCIPaddRowprepTerms(scip, rowprep, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
   }

   /* check if reference point violates cut at least a little bit */
   if( success && !SCIPisInfinity(scip, -minefficacy) )
   {
      viol = SCIPgetRowprepViolation(scip, rowprep, sol);
      if( viol <= 0.0 ) /*lint !e644*/
      {
         SCIPdebugMsg(scip, "skip cut for constraint <%s> because efficacy %g too low (< %g)\n", SCIPconsGetName(cons), viol, minefficacy);
         success = FALSE;
      }
   }

   /* cleanup and improve cut */
   if( success )
   {
      SCIP_Real coefrange;

      /* merge terms */
      SCIPmergeRowprepTerms(scip, rowprep);

      /* improve coefficients */
      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, conshdlrdata->cutmaxrange, minefficacy, &coefrange, &viol) );
      success = coefrange <= conshdlrdata->cutmaxrange;
   }

   /* check that side is finite */ /*lint --e{514} */
   success &= !SCIPisInfinity(scip, REALABS(rowprep->side)); /*lint !e514*/

   /* check whether maximal coef is finite, if any */ /*lint --e{514} */
   success &= (rowprep->nvars == 0) || !SCIPisInfinity(scip, REALABS(rowprep->coefs[0])); /*lint !e514*/

   /* check if reference point violates cut sufficiently */
   if( success && !SCIPisInfinity(scip, -minefficacy) && viol < minefficacy ) /*lint !e644*/
   {
      SCIPdebugMsg(scip, "skip cut for constraint <%s> because efficacy %g too low (< %g)\n", SCIPconsGetName(cons), viol, minefficacy);
      success = FALSE;
   }

   /* generate row */
   if( success )
   {
      SCIP_CALL( SCIPgetRowprepRowCons(scip, row, rowprep, SCIPconsGetHdlr(cons)) );

      SCIPdebugMsg(scip, "found cut <%s>, lhs=%g, rhs=%g, mincoef=%g, maxcoef=%g, range=%g, nnz=%d, efficacy=%g\n",
         SCIProwGetName(*row), SCIProwGetLhs(*row), SCIProwGetRhs(*row),
         rowprep->nvars > 0 ? rowprep->coefs[rowprep->nvars-1] : 0.0, rowprep->nvars > 0 ? rowprep->coefs[0] : 0.0,
         rowprep->nvars > 0 ? rowprep->coefs[0]/rowprep->coefs[rowprep->nvars-1] : 1.0,
         SCIProwGetNNonz(*row), viol);  /*lint !e414 */

      if( efficacy != NULL )
         *efficacy = viol;
   }

   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** computes eigen decomposition of A, where \f$ f(x) = x^T A x + b^T x \f$.
 *
 * The eigen decomposition is given by A = P D P^T, where D is diagonal formed by the eigenvalues and P is orthonormal
 * whose columns are the eigenvectors; we also compute b^T * P, in case one needs the change of variables P^T x = y <=>
 * x = P y We store P^T in an array, specifically, in consdata->eigenvectors we store P^T row-wise, i.e., the first row
 * of P^T is stored in eigenvector[0..n-1], the second row is stored in eigenvectors[n..2n-1], etc; equivalently, the
 * first eigenvector is eigenvector[0..n-1], the second one is eigenvectors[n..2n-1], etc.
 *
 * @todo: - at the moment of writing, checkCurvature computes the eigenvalues (and vectors) for determining curvature
 *          when it can't to it via other considerations. so one could try to merge both methods together.
 *        - it seems that if A is of the form [I 0; 0 A'], one only needs to compute the decomposition for A' so one
 *          could do better in terms of memory and speed. For instance, when the matrix is diagonal, the eigenvectors
 *          are the identity matrix and the eigenvalues are readily available from the constraint, so one could adapt
 *          the functions that uses the eigenvectors in this particular case. One could also think about storing the
 *          eigenvectors in a sparse fashion, though eigenvectors are seldom sparse.
 */
static
SCIP_RETCODE computeED(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int n;
   int nn;
   int row;
   int col;
   int i;
   int j;
   double*        matrix;
   SCIP_HASHMAP*  var2index;

   SCIPdebugMsg(scip, "computing ED for cons %s\n", SCIPconsGetName(cons));

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* function has to be convex with finite rhs or concave with finite lhs */
   assert((consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs)) ||
         (consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)));

   /* can't compute eigenvectors without IPOPT */
   if( !SCIPisIpoptAvailableIpopt() )
   {
      consdata->isedavailable = FALSE;
      return SCIP_OKAY;
   }

   /* @todo: - it seems that if A is of the form [I 0; 0 A'], one only needs to compute the decomposition for A'
    *          so one could do better in terms of memory and speed
    *        - if n too big don't compute SVD
    */
   n = consdata->nquadvars;

   /* do not compute eigendecomposition if n is too large */
   nn = n * n;
   if( nn < 0 || (unsigned) (int) nn > UINT_MAX / sizeof(SCIP_Real) )
   {
      SCIPdebugMsg(scip, "n is too large to compute eigendecomposition\n");
      consdata->isedavailable = FALSE;
      return SCIP_OKAY;
   }

   /* we just need to pass the upper triangle of A since it is symmetric; build it here */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eigenvectors, nn) );
   matrix = consdata->eigenvectors;
   BMSclearMemoryArray(matrix, nn);

   /* @todo if we are called in solving stage (or late from initsol), we can avoid the hashmap by using sepabilinvar2pos */
   SCIP_CALL( SCIPhashmapCreate(&var2index, SCIPblkmem(scip), n) );

   for( i = 0; i < n; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(var2index, consdata->quadvarterms[i].var, (void*)(size_t)i) );
      matrix[i*n + i] = consdata->quadvarterms[i].sqrcoef;
#ifdef DEBUG_PROJ
      printf("inserting in position %d, value %g\n", i*n + i, consdata->quadvarterms[i].sqrcoef);
#endif
   }

   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var1));
      assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var2));
      row = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var1);
      col = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var2);
      if( row < col )
      {
         matrix[row * n + col] = consdata->bilinterms[i].coef/2;
#ifdef DEBUG_PROJ
         printf("inserting in position %d, value %g\n", row*n + col, consdata->bilinterms[i].coef/2);
#endif
      }
      else
      {
         matrix[col * n + row] = consdata->bilinterms[i].coef/2;
#ifdef DEBUG_PROJ
         printf("inserting in position %d, value %g\n", col*n + row, consdata->bilinterms[i].coef/2);
#endif
      }
   }

#ifdef DEBUG_PROJ
   printf("matrix built:\n");
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < n; j++ )
         printf("%g ", matrix[i*n + j]);
      printf("\n");
   }
#endif

   /* compute eigenvalues and eigenvectors */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eigenvalues, n) );

   if( LapackDsyev(TRUE, n, matrix, consdata->eigenvalues) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "couldn't compute ED for cons %s\n", SCIPconsGetName(cons));
      consdata->isedavailable = FALSE;
   }
   else
   {
      consdata->isedavailable = TRUE;

      /* compute b^T*P */
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &consdata->bp, n) );
      for( i = 0; i < n; i++ )
         for( j = 0; j < n; j++ )
            consdata->bp[i] += consdata->quadvarterms[j].lincoef * matrix[i*n + j];

#ifdef DEBUG_PROJ
      printf("eigenvalues:\n");
      for( j = 0; j < n; j++ )
         printf("%g ", consdata->eigenvalues[j]);

      printf("\neigenvectors (P^T):\n");
      for( i = 0; i < n; i++ )
      {
         for( j = 0; j < n; j++ )
            printf("%g ", matrix[i*n + j]);
         printf("\n");
      }

      printf("b*P^T:\n");
      for( j = 0; j < n; j++ )
         printf("%g ", consdata->bp[j]);
      printf("svd computed successfully\n");
#endif
   }

   SCIPhashmapFree(&var2index);

   return SCIP_OKAY;
}

/** computes an interior point for the quadratic part of the convex constraint
 *
 *  There are different methods for computing the interior point
 *  - 'a'ny: solves min 0, f(x) <= rhs, x in bounds
 *  - 'm'ost interior: solves min f(x), x in bounds
 *
 * @todo: other methods for computing an interior point?
 */
static
SCIP_RETCODE computeInteriorPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   char                  method,             /**< method for computing interior point ('a' any point, 'm'ost interior) */
   SCIP_Bool*            success             /**< buffer to store if an interior point was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_QUADELEM* nlrowquadelems;
   SCIP_NLPIPROBLEM* prob;
   SCIP_NLPI* nlpi;
   SCIP_Real* interiorpoint;
   SCIP_Real* lbs;
   SCIP_Real* ubs;
   SCIP_Real* lincoefs;
   SCIP_Real nlpiside;
   char probname[SCIP_MAXSTRLEN];
   int* lininds;
   int nlrownquadelems;
   int nquadvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   assert(success != NULL);
   *success = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert((consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs)) ||
         (consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)));

   /* need an NLP solver */
   if( SCIPgetNNlpis(scip) == 0 )
      return SCIP_OKAY;

   nlpi = NULL;
   prob = NULL;
   lbs = NULL;
   ubs = NULL;
   lincoefs = NULL;
   lininds = NULL;

#ifdef SCIP_DEBUG_INT
   SCIPinfoMessage(scip, NULL, "Computing interior point for\n");
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
   SCIPinfoMessage(scip, NULL, ";\n");
#endif

   /* in the convex case, we try to find an interior point of x^T A x + b^T x <= rhs - maximum activity linear part
    * in the concave case: lhs - minimum activity linear part <= x^T A x + b^T x; we compute activities ourselves,
    * since consdata->max(min)linactivity are only computed when lhs (rhs) is finite and this not always holds
    */
   if( consdata->isconvex )
   {
      /* compute maximum activity */
      nlpiside = 0;
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         if( consdata->lincoefs[i] >= 0.0 )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->linvars[i]) ) )
               nlpiside = SCIPinfinity(scip);
            else
               nlpiside += consdata->lincoefs[i] * SCIPvarGetUbLocal(consdata->linvars[i]);
         }
         else
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->linvars[i]) ) )
               nlpiside = SCIPinfinity(scip);
            else
               nlpiside += consdata->lincoefs[i] * SCIPvarGetLbLocal(consdata->linvars[i]);
         }

         if( SCIPisInfinity(scip, nlpiside) )
         {
            SCIPdebugMsg(scip, "maximum activity is infinity: there is no interior point for fun <= rhs - maxlinactivity!\n");
            return SCIP_OKAY;
         }
      }

      if( consdata->nlinvars == 0 )
         nlpiside = INTERIOR_EPS;

      nlpiside = consdata->rhs - nlpiside;
   }
   else
   {
      /* compute minimum activity */
      nlpiside = 0;
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         if( consdata->lincoefs[i] >= 0.0 )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->linvars[i])) )
               nlpiside = -SCIPinfinity(scip);
            else
               nlpiside += consdata->lincoefs[i] * SCIPvarGetLbLocal(consdata->linvars[i]);
         }
         else
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->linvars[i])) )
               nlpiside = -SCIPinfinity(scip);
            else
               nlpiside += consdata->lincoefs[i] * SCIPvarGetUbLocal(consdata->linvars[i]);
         }

         if( SCIPisInfinity(scip,  -nlpiside) )
         {
            SCIPdebugMsg(scip, "minimum activity is -infinity: there is no interior point for fun >= lhs - minlinactivity!\n");
            return SCIP_OKAY;
         }
      }

      if( consdata->nlinvars == 0 )
         nlpiside = INTERIOR_EPS;

      nlpiside = consdata->lhs - nlpiside;
   }

   nquadvars = consdata->nquadvars;

   /* if we are looking for any interior point and the 0 is one, then use it */
   if( method == 'a' && ((consdata->isconvex && SCIPisGE(scip, nlpiside, 0.0))
            || (consdata->isconcave && SCIPisLE(scip, nlpiside, 0.0))) )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(consdata->interiorpoint), nquadvars) );

      *success = TRUE;
      goto TERMINATE;
   }

   /* build nlrow */
   if( consdata->nlrow == NULL )
   {
      SCIP_CALL( createNlRow(scip, cons) );
      assert(consdata->nlrow != NULL);
   }

   nlpi = SCIPgetNlpis(scip)[0];
   assert(nlpi != NULL);

   /* initializing the subproblem */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_subquad", SCIPgetProbName(scip));
   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &prob, probname) );
   assert(prob != NULL);

#ifdef SCIP_DEBUG_INT
   SCIP_CALL( SCIPnlpiSetIntPar(nlpi, prob, SCIP_NLPPAR_VERBLEVEL, 0) );
#endif
   /* TODO: maybe one should set some generous iteration limit and/or a timelimit (remaining scip solve time)? */

   /* ask for memory to store data needed to create vars and linear coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &lbs, nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nquadvars) );

   /* get bounds and linear coefficients */
   for( i = 0; i < nquadvars; i++ )
   {
      lbs[i] = SCIPvarGetLbGlobal(consdata->quadvarterms[i].var);
      ubs[i] = SCIPvarGetUbGlobal(consdata->quadvarterms[i].var);

      lincoefs[i] = consdata->quadvarterms[i].lincoef;
      lininds[i] = i;
   }

   /* add vars */
   SCIP_CALL( SCIPnlpiAddVars(nlpi, prob, nquadvars, lbs, ubs, NULL) );

   /* get nlrow info */
   nlrownquadelems = SCIPnlrowGetNQuadElems(consdata->nlrow);
   nlrowquadelems = SCIPnlrowGetQuadElems(consdata->nlrow);

#ifndef NDEBUG
   {
      SCIP_VAR** nlrowquadvars;

      nlrowquadvars = SCIPnlrowGetQuadVars(consdata->nlrow);
      for( i = 0; i < nlrownquadelems; i++ )
      {
         assert(nlrowquadvars[nlrowquadelems[i].idx1] == consdata->quadvarterms[nlrowquadelems[i].idx1].var);
         assert(nlrowquadvars[nlrowquadelems[i].idx2] == consdata->quadvarterms[nlrowquadelems[i].idx2].var);
      }
   }
#endif

   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s", SCIPconsGetName(cons));

   switch( method )
   {
      case 'a':
         /* add constraint */
         if( consdata->isconvex )
         {
            SCIP_CALL( SCIPnlpiAddConstraints(nlpi, prob, 1, NULL, &nlpiside, &nquadvars, &lininds, &lincoefs,
                     &nlrownquadelems, &nlrowquadelems, NULL, NULL, NULL) );
         }
         else
         {
            SCIP_CALL( SCIPnlpiAddConstraints(nlpi, prob, 1, &nlpiside, NULL, &nquadvars, &lininds, &lincoefs,
                     &nlrownquadelems, &nlrowquadelems, NULL, NULL, NULL) );
         }
         break;

      case 'm':
         /* add objective */
         if( consdata->isconvex )
         {
            SCIP_CALL( SCIPnlpiSetObjective(nlpi, prob, nquadvars, lininds, lincoefs,
                     nlrownquadelems, nlrowquadelems, NULL, NULL, 0.0) );
         }
         else
         {
            /* NLPI assumes minimization: change signs */
            for( i = 0; i < nquadvars; i++ )
               lincoefs[i] *= -1;

            /* WARNING: this pointer is not ours, information should be restored! */
            for( i = 0; i < nlrownquadelems; i++ )
               nlrowquadelems->coef *= -1;

            SCIP_CALL( SCIPnlpiSetObjective(nlpi, prob, nquadvars, lininds, lincoefs,
                     nlrownquadelems, nlrowquadelems, NULL, NULL, 0.0) );

            /* WARNING: restore information! */
            for( i = 0; i < nlrownquadelems; i++ )
               nlrowquadelems->coef *= -1;
         }
         break;

      default:
         SCIPerrorMessage("undefined method for computing interior point: %c\n", method);
         return SCIP_INVALIDDATA;
   }

   /* set NLP tolerances; we don't really need an optimal solution to this NLP */
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, prob, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip)) ); /*lint !e666*/
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, prob, SCIP_NLPPAR_RELOBJTOL, MAX(SCIPfeastol(scip), SCIPdualfeastol(scip))) ); /*lint !e666*/

   /* solve NLP problem */
   SCIP_CALL( SCIPnlpiSolve(nlpi, prob) );

   /* check termination status */
   if( SCIPnlpiGetTermstat(nlpi, prob) != SCIP_NLPTERMSTAT_OKAY )
   {
      SCIPdebugMsg(scip, "cons <%s>: NLP Solver termination status not okay: %d\n",
            SCIPconsGetName(cons), SCIPnlpiGetTermstat(nlpi, prob));
      *success = FALSE;
      goto TERMINATE;
   }

   /* check solution status */
   switch( SCIPnlpiGetSolstat(nlpi, prob) )
   {
      case SCIP_NLPSOLSTAT_GLOBOPT:
      case SCIP_NLPSOLSTAT_LOCOPT:
      case SCIP_NLPSOLSTAT_FEASIBLE:
         /* fallthrough */
         SCIPdebugMsg(scip, "cons <%s>: found an interior point.  solution status: %d, termination status: %d\n",
               SCIPconsGetName(cons), SCIPnlpiGetSolstat(nlpi, prob), SCIPnlpiGetTermstat(nlpi, prob));
         break;

      case SCIP_NLPSOLSTAT_LOCINFEASIBLE:
      case SCIP_NLPSOLSTAT_GLOBINFEASIBLE:
      case SCIP_NLPSOLSTAT_UNKNOWN:
         /* fallthrough */
         /* TODO: we could still use the point, and let evaluateGauge decide whether the point is interior or not */
         SCIPdebugMsg(scip, "cons <%s>: failed to find an interior point.  solution status: %d, termination status: %d\n",
               SCIPconsGetName(cons), SCIPnlpiGetSolstat(nlpi, prob), SCIPnlpiGetTermstat(nlpi, prob));
         goto TERMINATE;

      case SCIP_NLPSOLSTAT_UNBOUNDED:
      default:
         /* fallthrough */
         SCIPerrorMessage("cons <%s>: undefined behaviour of NLP Solver.  solution status: %d, termination status: %d\n",
               SCIPconsGetName(cons), SCIPnlpiGetSolstat(nlpi, prob), SCIPnlpiGetTermstat(nlpi, prob));
         SCIPABORT();
         goto TERMINATE; /*lint !e527*/
   }

   /* fetch solution
    * note: nlpiGetSolution (at least for IPOPT) makes interiorpoint point to the internal solution stored in the
    * nlpi problem data structure; we need to copy it here because it will be destroyed once the problem is free'd
    */
   SCIP_CALL( SCIPnlpiGetSolution(nlpi, prob, &interiorpoint, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->interiorpoint), nquadvars) );

   for( i = 0; i < nquadvars; i++ )
   {
      if( SCIPisFeasZero(scip, interiorpoint[i]) )
         consdata->interiorpoint[i] = 0.0;
      else
         consdata->interiorpoint[i] = interiorpoint[i];
   }

   *success = TRUE;

TERMINATE:

#ifdef SCIP_DEBUG_INT
   printf("Computation of interior point for cons <%s>:\n", SCIPconsGetName(cons));
   printf(" - has %d linear variables\n", consdata->nlinvars);
   if( consdata->isconvex )
   {
      printf(" - is convex. rhs: %g maximum activity of linear variables: %g\n", consdata->rhs, consdata->rhs - nlpiside);
      printf(" - searched for point whose quadratic part is <= %g\n", nlpiside);
   }
   else
   {
      printf(" - is concave. lhs: %g minimum activity of linear variables: %g\n", consdata->lhs, consdata->lhs - nlpiside);
      printf(" - searched for point whose quadratic part is >= %g\n", nlpiside);
   }

   if( *success )
   {
      if( prob == NULL )
      {
         printf("Computation successful, 0 is interior point.\n");
         for( i = 0; i < nquadvars; i++ )
         {
            assert(consdata->interiorpoint[i] == 0.0);
         }
      }
      else
      {
         printf("Computation successful, NLP soltat: %d, termstat: %d\nPoint found:\n",
               SCIPnlpiGetSolstat(nlpi, prob), SCIPnlpiGetTermstat(nlpi, prob));
         for( i = 0; i < nquadvars; i++ )
         {
            printf("%s = %g\n", SCIPvarGetName(consdata->quadvarterms[i].var), consdata->interiorpoint[i]);
         }
      }
   }
   else
   {
      printf("Computation failed. NLP soltat: %d, termstat: %d\n",
               SCIPnlpiGetSolstat(nlpi, prob), SCIPnlpiGetTermstat(nlpi, prob));
      printf("run with SCIP_DEBUG for more info\n");
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, ";\n");
      /* FIXME: instance camshape100 says that there is no interior point (interior empty)
       * is there something intelligent that can be said?
       */
   }
#endif

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &lbs);
   SCIPfreeBufferArrayNull(scip, &ubs);
   SCIPfreeBufferArrayNull(scip, &lininds);
   SCIPfreeBufferArrayNull(scip, &lincoefs);

   if( prob != NULL )
   {
      SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &prob) );
   }

   return SCIP_OKAY;
}

/** compute gauge function of the set \f$S - s_0\f$ where \f$ S = \{ x : f(x) \le c \}\f$ and \f$ s_0 \in \mathring S\f$.
 *
 * Here, \f$ f(x) \f$ is a purely quadratic (i.e, all \f$x\f$ variables appear in a bilinear or quadratic term).
 * Explicitly, \f$ f(x) = \pm x^T A x \pm b^T x \f$ depending whether \f$A\f$
 * is positive semidefinite (+) or negative semidefinite (-).
 * The constant \f$c\f$ is rhs - maximum activity of the purely linear part of the constraint
 * if \f$A \succeq 0\f$ and minimum activity - lhs if \f$A \preceq 0\f$.
 * This is computed only at INITSOL.
 *
 * The method does:
 * 1. compute interior point
 * 2. compute gauge function
 */
static
SCIP_RETCODE computeGauge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_BILINTERM* bilinterm;
   SCIP_Bool success;
   SCIP_Bool convex;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->gaugecuts);

   /* function has to be convex with finite rhs or concave with finite lhs */
   convex = consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs);
   assert(convex || (consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)));

   SCIPdebugMsg(scip, "cons %s: is %s\n", SCIPconsGetName(cons), convex ? "convex" : "concave");

   /* 1. */
   SCIP_CALL( computeInteriorPoint(scip, cons, conshdlrdata->interiorcomputation, &success) );

   /* if success, compute gaugecoefs (b_gauge) and gaugeconst (c_gauge) */
   if( !success )
   {
      SCIPdebugMsg(scip, "failed to compute gauge function\n");
      consdata->isgaugeavailable = FALSE;
      return SCIP_OKAY;
   }

   /* 2.
    * we are going to evaluate the function at interiorpoint; so, we need to compute interiorpoint^T A interiorpoint;
    * therefore, we need a mechanism that for a given variable, it returns its interior point value
    * fortunately, sepabilinvar2pos in consdata gives us all the information that we need
    */

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(consdata->gaugecoefs), consdata->nquadvars) );

   /* compute value of quadratic part at interior point, build map and compute gaugeconst (c_gauge) */
   consdata->interiorpointval = 0;
   consdata->gaugeconst = 0;
   for( i = 0; i < consdata->nquadvars; i++ )
   {
      SCIP_Real val;
      SCIP_Real val2;

      val = consdata->interiorpoint[i];
      quadvarterm = &consdata->quadvarterms[i];

      consdata->interiorpointval += (quadvarterm->lincoef + quadvarterm->sqrcoef * val) * val;
      consdata->gaugeconst += quadvarterm->sqrcoef * val * val;

      for( j = 0; j < quadvarterm->nadjbilin; ++j )
      {
         int bilintermidx;

         bilintermidx = quadvarterm->adjbilin[j];
         bilinterm = &consdata->bilinterms[bilintermidx];

         if( bilinterm->var1 != quadvarterm->var )
            continue;

         /* the index of the variable associated with var2 in bilinterm should be given by sepabilinvar2pos */
         assert(consdata->sepabilinvar2pos != NULL); /* this should have been computed in INITSOL */
         assert(consdata->quadvarterms[consdata->sepabilinvar2pos[bilintermidx]].var == bilinterm->var2);

         val2 = consdata->interiorpoint[consdata->sepabilinvar2pos[bilintermidx]];

         consdata->interiorpointval += bilinterm->coef * val * val2;
         consdata->gaugeconst += bilinterm->coef * val * val2;
      }
   }

   /* compute gaugecoefs (b_gauge = b + 2 * A * interiorpoint) */
   for( i = 0; i < consdata->nquadvars; i++ )
   {
      quadvarterm = &consdata->quadvarterms[i];
      consdata->gaugecoefs[i] += quadvarterm->lincoef + 2.0 * quadvarterm->sqrcoef * consdata->interiorpoint[i];

      for( j = 0; j < quadvarterm->nadjbilin; j++ )
      {
         int varpos;
         int bilintermidx;

         bilintermidx = quadvarterm->adjbilin[j];
         bilinterm = &consdata->bilinterms[bilintermidx];

         if( bilinterm->var1 == quadvarterm->var )
         {
            varpos = consdata->sepabilinvar2pos[bilintermidx];

            /* the index of the variable associated with var2 in bilinterm should be given by sepabilinvar2pos */
            assert(consdata->quadvarterms[varpos].var == bilinterm->var2);

            consdata->gaugecoefs[i] += bilinterm->coef * consdata->interiorpoint[varpos];
            consdata->gaugecoefs[varpos] += bilinterm->coef * consdata->interiorpoint[i];
         }
      }
   }

#ifdef SCIP_DEBUG_INT
   printf("quadratic part at interior point: %g\n", consdata->interiorpointval);

   for( j = 0; j < consdata->nquadvars; j++ )
   {
      printf("b_gauge[%s] = %g\n", SCIPvarGetName(consdata->quadvarterms[j].var), consdata->gaugecoefs[j]);
   }
   printf("c_gauge = %g\n", consdata->gaugeconst);
#endif

   SCIPdebugMsg(scip, "gauge function computed successfully\n");
   consdata->isgaugeavailable = TRUE;

   return SCIP_OKAY;
}

/** evaluates gauge function of the set \f$S - s_0\f$ where \f$ S = \{ x : f(x) \le c \}\f$ and \f$ s_0 \in \mathring S\f$.
 *
 * \f$ S = \{ x : f(x) \le c \}\f$ at \f$sol - s_0\f$;
 * see computeGauge() for more details
 *
 * @todo Think about if user should tell that function is convex or ...
 */
static
SCIP_RETCODE evaluateGauge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             refsol,             /**< reference point where to generate cut, or NULL if sol should be used */
   SCIP_Real*            gaugeval,           /**< buffer to store the value of the gauge function */
   SCIP_Bool*            success             /**< buffer to store if evaluation was successful */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real side;
   SCIP_Real aterm;
   SCIP_Real bterm;
   SCIP_Real cterm;
   SCIP_Bool convex;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->isgaugeavailable);

   *success = FALSE;

   convex = consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs);

   SCIPdebugMsg(scip, "cons %s: is %s\n", SCIPconsGetName(cons), convex ? "convex" : "concave");

   /* evaluate gauge function at x0 = (refsol - interior point)
    *
    * compute aterm = side - function(interior point)
    */
   if( convex )
   {
      side = consdata->rhs;
      for( i = 0; i < consdata->nlinvars; i++ )
         side -= SCIPgetSolVal(scip, refsol, consdata->linvars[i]) * consdata->lincoefs[i];

      aterm = side - consdata->interiorpointval;

      /* it can happen that the interior point is not really interior, since we are not so strict at the moment of
       * computing the interior point, which makes sense in the case that the constraint is quadratic <= linear expr,
       * since we compute a point in quadratic <= min linear expr and it might be that this set consists of a single
       * point which will not be interior. furthermore, if this set is empty, we could just take any point and it could
       * happen that for some value of linear expr, the point is actually interior, but for many it could not be.
       * also, if min linear expr = -infinity, we might have computed an interior point using some finite value.
       * the point will not be an interior point, if and only if aterm is negative.
       */
#ifdef SCIP_DEBUG_GAUGE
      if( SCIPisLE(scip, aterm, 0.0) )
      {
         printf("For current level, there is no interior point. ");
         printf("rhs: %g level: %15.20g interiorpointval: %15.20g\n", consdata->rhs, side, consdata->interiorpointval);
         if( consdata->nlinvars == 1 )
         {
            SCIP_VAR* var;

            var = consdata->linvars[0];
            printf("var <%s> = %g in [%15.20g, %15.20g] is linpart\n", SCIPvarGetName(var),
                  SCIPgetSolVal(scip, refsol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         }
      }
      else
      {
         printf("For current level, there is interior point. ");
         printf("rhs: %g level: %15.20g interiorpointval: %15.20g\n", consdata->rhs, side, consdata->interiorpointval);
      }
#endif
      if( !SCIPisPositive(scip, aterm) )
      {
         *gaugeval = -1.0;
         return SCIP_OKAY;
      }
   }
   else
   {
      side = consdata->lhs;
      for( i = 0; i < consdata->nlinvars; i++ )
         side -= SCIPgetSolVal(scip, refsol, consdata->linvars[i]) * consdata->lincoefs[i];

      aterm = side - consdata->interiorpointval;

#ifdef SCIP_DEBUG_GAUGE
      if( SCIPisGE(scip, aterm, 0.0) )
      {
         printf("For current level, there is no interior point. ");
         printf("lhs: %g level: %15.20g interiorpointval: %15.20g\n", consdata->lhs, side, consdata->interiorpointval);
         if( consdata->nlinvars == 1 )
         {
            SCIP_VAR* var;

            var = consdata->linvars[0];
            printf("var <%s> = %g in [%15.20g, %15.20g] is linpart\n", SCIPvarGetName(var),
                  SCIPgetSolVal(scip, refsol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         }
      }
      else
      {
         printf("For current level, there is interior point. ");
         printf("lhs: %g level: %15.20g interiorpointval: %15.20g\n", consdata->lhs, side, consdata->interiorpointval);
      }
#endif
      if( !SCIPisNegative(scip, aterm) )
      {
         *gaugeval = -1.0;
         return SCIP_OKAY;
      }
   }

   /* compute bterm = b_gauge^T * refsol - f(interiorpoint) - c_gauge
    * compute cterm = f(refsol) - b_gauge^T * refsol + c_gauge */
   bterm = -consdata->interiorpointval - consdata->gaugeconst;
   cterm = consdata->gaugeconst;
   for( i = 0; i < consdata->nquadvars; i++ )
   {
      SCIP_Real val;

      val = SCIPgetSolVal(scip, refsol, consdata->quadvarterms[i].var);
      bterm += consdata->gaugecoefs[i] * val;
      cterm -= consdata->gaugecoefs[i] * val;
      cterm += (consdata->quadvarterms[i].lincoef + consdata->quadvarterms[i].sqrcoef * val) * val;
   }

   for( i = 0; i < consdata->nbilinterms; i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      var1 = consdata->bilinterms[i].var1;
      var2 = consdata->bilinterms[i].var2;
      cterm += consdata->bilinterms[i].coef * SCIPgetSolVal(scip, refsol, var1) * SCIPgetSolVal(scip, refsol, var2);
   }

   /* now compute gauge */
   if( convex && cterm < 0.0 )
   {
      assert(SCIPisZero(scip, cterm));
      cterm = 0.0;
   }
   else if( !convex && cterm > 0.0 )
   {
      assert(SCIPisZero(scip, cterm));
      cterm = 0.0;
   }
   assert(bterm*bterm + 4*aterm*cterm >= 0);

   if( convex )
   {
      *gaugeval = bterm + sqrt(bterm*bterm + 4 * aterm * cterm);
      *gaugeval = *gaugeval / (2 * aterm);
   }
   else
   {
      *gaugeval = bterm - sqrt(bterm*bterm + 4 * aterm * cterm);
      *gaugeval = *gaugeval / (2 * aterm);
   }
   assert(!SCIPisNegative(scip, *gaugeval));
   *success = TRUE;

#ifdef SCIP_DEBUG_GAUGE
   printf("Gauge's aterm = %g, bterm = %g, cterm = %g\n", aterm, bterm, cterm);
#endif
   return SCIP_OKAY;
}

/** compute projection of refsol onto feasible region of cons; stores the projection in ref
 *
 * This method solves
 * \f[
 *      \min \{ ||x - \bar x||^2 : x^T A x + 2 b^T x \le c \}
 * \f]
 * where \f$ \bar x \f$ is refsol.
 * Note that \f$ \bar x \f$ is not feasible, so the optimal solution actually satisfies
 * \f[
 *      \min \{ ||x - \bar x||^2 : x^T A x + 2 b^T x = c \}
 * \f]
 * Using the eigendecomposition \f$ A = P D P^T \f$, the change of variables \f$ y = P^T x
 * \f$ and the optimality conditions, this reduces to finding \f$ \rho \f$ such that
 * \f[
 *      y(\rho) = (I + \rho D)^{-1} (\bar y - \rho \bar b)
 * \f]
 * makes the constraint active. In the previous formula, \f$ \bar y = P^T \bar x\f$ and \f$ \bar b = P^T b \f$.  If \f$
 * D \neq 0 \f$, the function
 * \f[
 *   \varphi(\rho) := y(\rho)^T D y(\rho) + 2 \bar b^T y(\rho) - c
 * \f]
 * is strictly convex. So this method actually computes the unique 0 of this function using Newton's method.
 */
static
SCIP_RETCODE computeReferencePointProjection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             refsol,             /**< the given point to project, or NULL if LP solution should be used */
   SCIP_Real*            ref                 /**< array to store reference point */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* pt; /* stores P^T */
   SCIP_Real* bp;
   SCIP_Real* D;
   SCIP_Real* y0_;
   SCIP_Real* yrho;
   SCIP_Real* yrhoprime;
   SCIP_Real c;
   SCIP_Real c1;
   SCIP_Real c2;
   SCIP_Real rho;
   SCIP_Real phirho;
   SCIP_Real phirhoprime;
   SCIP_Bool isconcave;
   int iter;
   int i;
   int j;
   int n;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->isedavailable);

   SCIPdebugMessage("computing projection\n");

   /* get the data we need */
   pt = consdata->eigenvectors;
   D  = consdata->eigenvalues;
   n  = consdata->nquadvars;
   bp = consdata->bp;
   c  = consdata->rhs;
   c1 = 0;
   c2 = 0;
   for( i = 0; i < consdata->nlinvars; i++ )
   {
      c1 += consdata->lincoefs[i] * SCIPgetSolVal(scip, refsol, consdata->linvars[i]);
      c2 -= consdata->lincoefs[i] * consdata->lincoefs[i];
   }
   c2 /= 2.0;

   /* determine if convex or concave */
   isconcave = consdata->isconcave;
   assert((isconcave && !SCIPisInfinity(scip, -consdata->lhs)) || !SCIPisInfinity(scip, consdata->rhs));

   SCIP_CALL( SCIPallocClearBufferArray(scip, &y0_, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &yrho, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &yrhoprime, n) );

   /* change data if function is concave */
   if( isconcave )
   {
      c  = -consdata->lhs;
      c1 = - c1;
      for( i = 0; i < n; i++ )
      {
         D[i]  = -D[i];
         bp[i] = -bp[i];
      }
   }

   /* change coordinates: compute y(0) = x_0' * P */
   for( i = 0; i < n; i++ )
      for( j = 0; j < n; j++ )
         y0_[i] += SCIPgetSolVal(scip, refsol, consdata->quadvarterms[j].var) * pt[i*n + j];

#ifdef DEBUG_PROJ
   /* debug output */
   printf("\nP^T:\n");
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < n; j++ )
         printf("%g ", pt[i*n + j]);
      printf("\n");
   }
   printf("x_0: ");
   for( i = 0; i < n; i++ )
      printf("%g ", SCIPgetSolVal(scip, refsol, consdata->quadvarterms[i].var));
   printf("\n");
   printf("P^T x_0: ");
   for( i = 0; i < n; i++ )
      printf("%g ", y0_[i]);
   printf("\n");
   printf("P^T b: ");
   for( i = 0; i < n; i++ )
      printf("%g ", bp[i]);
   printf("\n");
   printf("<d,linvars> = %g\n", c1);
   printf("-norm(d)^2/2 = %g\n", c2);
#endif

   /* perform newton's method:  rho^+ = rho - phi(rho)/phi'(rho) */
   rho = 0.0;
   phirho = c;
   phirhoprime = 1.0;
   for( iter = 0; iter < 9; iter++ )
   {
      assert(phirhoprime != 0.0);
      rho = rho - (phirho - c)/ phirhoprime;

      /* compute phi(rho) and phi'(rho):
       * note that formulas were deduced for constraints of the form x' A x + 2 b x, so we use b/2 in the formulas:
       * c1        = <lin_coefs, sol_lin_vars>
       * c2        = - norm(lin_coefs)^2/2
       * y(rho)    = (I + rho * D)^-1 * (y(0) - rho * bp/2)
       * y'(rho)   = -(I + rho * D)^-2 * (D y(0) + bp/2)
       * phi(rho)  = <y(rho), D * y(rho) + pb> + c1 + c2*rho
       * phi'(rho) = <y'(rho), 2 * D * y(rho) + pb> + c2
       */
      phirho = 0.0;
      phirhoprime = 0.0;
      for( i = 0; i < n; i++ )
      {
         assert(1.0 + rho * D[i] != 0.0);
         yrho[i]      = (y0_[i] - rho * bp[i]/2.0) / (1.0 + rho * D[i]);
         yrhoprime[i] = -(D[i] * y0_[i] + bp[i]/2.0) / ( (1.0 + rho * D[i])*(1.0 + rho * D[i]) );
         phirho      += yrho[i] * (yrho[i] * D[i] + bp[i]);
         phirhoprime += yrhoprime[i] * (2 * D[i] * yrho[i] + bp[i]);
      }
      phirho      += c2 * rho + c1;
      phirhoprime += c2;
#ifdef DEBUG_PROJ
      printf("iteration %d: rho = %g, phirho = %g, phirho' = %g\n", iter, rho, phirho, phirhoprime);
#endif
   }

   /* come back to the original coordinates: new ref point is P*yrho */
   for( i = 0; i < n; i++ )
   {
      ref[i] = 0.0;

      for( j = 0; j < n; j++ )
         ref[i] += pt[j*n + i] * yrho[j];
   }

   /* change data back if function is concave */
   if( isconcave )
   {
      for( i = 0; i < n; i++ )
      {
         D[i]  = -D[i];
         bp[i] = -bp[i];
      }
   }

#ifdef SCIP_DISABLED_CODE
   /* project onto bounds; this is important for some cut generation methods such as generateCutLTI */
   for( j = 0; j < consdata->nquadvars; ++j )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_VAR* var;

      var = consdata->quadvarterms[j].var;
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);
      /* do not like variables at infinity */
      assert(!SCIPisInfinity(scip,  lb));
      assert(!SCIPisInfinity(scip, -ub));

      ref[j] = MIN(ub, MAX(lb, ref[j])); /* project value into bounds */
   }
#endif

#ifdef DEBUG_PROJ
   printf("modified reference point by a projection:\n");
   for( j = 0; j < consdata->nquadvars; ++j )
   {
      printf("%s = %g\n", SCIPvarGetName(consdata->quadvarterms[j].var), ref[j]);
   }
#endif

   SCIPfreeBufferArray(scip, &y0_);
   SCIPfreeBufferArray(scip, &yrho);
   SCIPfreeBufferArray(scip, &yrhoprime);

   return SCIP_OKAY;
}

/** compute reference point suggested by gauge function */
static
SCIP_RETCODE computeReferencePointGauge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             refsol,             /**< reference point where to compute gauge, or NULL if LP solution should be used */
   SCIP_Real*            ref,                /**< array to store reference point */
   SCIP_Bool*            success             /**< buffer to store whether we succeeded computing reference point */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real gaugeval;
   SCIP_Real intpoint;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_VAR* var;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->isgaugeavailable);

   SCIPdebugMsg(scip, "evaluating gauge\n");
   SCIP_CALL( evaluateGauge(scip, conshdlr, cons, refsol, &gaugeval, success) );

   if( !(*success) )
   {
#ifdef SCIP_DEBUG_GAUGE
      printf("Couldn't evaluate gauge!\n");
#endif
      return SCIP_OKAY;
   }

#ifdef SCIP_DEBUG_GAUGE
   {
      SCIP_Real level;

      level = consdata->rhs;
      for( j = 0; j < consdata->nlinvars; j++ )
         level -= SCIPgetSolVal(scip, refsol, consdata->linvars[j]) * consdata->lincoefs[j];

      printf("Summary:\n");
      printf("For cons <%s>: gauge at level %g evaluated at (refsol - intpoint) is %.10f\n",
            SCIPconsGetName(cons), level, gaugeval);
      printf("refsol - intpoint:\n");

      for( j = 0; j < consdata->nquadvars; ++j )
      {
         SCIP_VAR* vvar;
         vvar = consdata->quadvarterms[j].var;
         printf("%s: % 20.15g  - %g = %g\n", SCIPvarGetName(vvar), SCIPgetSolVal(scip, refsol, vvar),
               consdata->interiorpoint[j], SCIPgetSolVal(scip, refsol, vvar) - consdata->interiorpoint[j]);
      }
      if( SCIPisFeasLE(scip, gaugeval, 1.0) )
         printf("refsol is in the closure of the region (gaugeval <= 1), don't modify reference point\n");
   }
#endif

   /* scale gauge value so that final point is close to the boundary, but not on the boundary (weakens the cut) */
   gaugeval *= GAUGESCALE;

   /* if the point is not sufficiently violated, we don't modify it */
   if( SCIPisFeasLE(scip, gaugeval, 1.0) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* set reference to (refsol - interior point)/gaugeval + interior point and project onto bounds this is important for
    * some cut generation methods such as generateCutLTI
    * @todo remove the projection onto the bounds; generateCutLTI shouldn't be called for convex constraints
    */
   for( j = 0; j < consdata->nquadvars; ++j )
   {
      var = consdata->quadvarterms[j].var;
      lb  = SCIPvarGetLbLocal(var);
      ub  = SCIPvarGetUbLocal(var);
      /* do not like variables at infinity */
      assert(!SCIPisInfinity(scip,  lb));
      assert(!SCIPisInfinity(scip, -ub));

      intpoint = consdata->interiorpoint[j];
      ref[j] = (SCIPgetSolVal(scip, refsol, var) - intpoint) / gaugeval + intpoint;
      ref[j] = MIN(ub, MAX(lb, ref[j])); /* project value into bounds */
   }

#ifdef SCIP_DEBUG_GAUGE
   printf("successful application of guage: %g\n", gaugeval);
   printf("modified reference point:\n");
   for( j = 0; j < consdata->nquadvars; ++j )
   {
      printf("%s = % 20.15g\n", SCIPvarGetName(consdata->quadvarterms[j].var), ref[j]);
   }
#endif

   return SCIP_OKAY;
}

/** generates a cut based on linearization (if convex) or McCormick (if nonconvex) in a solution
 * @note mode indicates whether we should modify the point we want to cutoff (sol) via gauge or projection,
 * or if just normal linearization should be use, or the default way (whatever is specified via settings)
 */
static
SCIP_RETCODE generateCutSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution where to generate cut, or NULL if LP solution should be used */
   SCIP_SOL*             refsol,             /**< reference point where to generate cut, or NULL if sol should be used */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real*            efficacy,           /**< buffer to store efficacy of row in reference solution, or NULL if not of interest */
   SCIP_Bool             checkcurvmultivar,  /**< are we allowed to check the curvature of a multivariate quadratic function, if not done yet */
   SCIP_Real             minefficacy,        /**< minimal required efficacy */
   char                  mode                /**< mode of execution 'g'auge, 'p'rojection, 'l'inearization gradient, 'd'efault */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR*  var;
   SCIP_Real  lb;
   SCIP_Real  ub;
   SCIP_Real* ref;
   SCIP_Bool success;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( refsol == NULL )
      refsol = sol;

   /* get reference point */
   SCIP_CALL( SCIPallocBufferArray(scip, &ref, consdata->nquadvars) );
   success = FALSE;

   if( mode == 'd')
   {
      if( (consdata->isconvex && violside == SCIP_SIDETYPE_RIGHT) ||
            (consdata->isconcave && violside == SCIP_SIDETYPE_LEFT) )
      {
         if( conshdlrdata->gaugecuts && consdata->isgaugeavailable )
         {
            SCIP_CALL( computeReferencePointGauge(scip, conshdlr, cons, refsol, ref, &success) );
         }
         else if( conshdlrdata->projectedcuts && consdata->isedavailable )
         {
            SCIPdebugMessage("use the projection of refsol onto the region defined by the constraint as reference point\n");
            SCIP_CALL( computeReferencePointProjection(scip, conshdlr, cons, refsol, ref) );
            success = TRUE;
         }
      }

      if( success )
      {
         SCIP_CALL( generateCut(scip, conshdlr, cons, ref, sol, violside, row, efficacy, checkcurvmultivar, minefficacy) );

         /* if cut fails, try again without modifying reference point */
         if( *row == NULL || (efficacy != NULL && !SCIPisGT(scip, *efficacy, minefficacy)) || !SCIPisCutApplicable(scip, *row) )  /*lint !e644 */
         {
            SCIPdebugMsg(scip, "%s cut fail, try without modifying\n", conshdlrdata->gaugecuts ? "gauge" : "projected");
            success = FALSE;
         }
      }

      /* note that this is not the same as calling this method with mode 'l', 'l' assume convex/concave function */
      if( !success )
      {
         for( j = 0; j < consdata->nquadvars; ++j )
         {
            var = consdata->quadvarterms[j].var;
            lb  = SCIPvarGetLbLocal(var);
            ub  = SCIPvarGetUbLocal(var);
            /* do not like variables at infinity */
            assert(!SCIPisInfinity(scip,  lb));
            assert(!SCIPisInfinity(scip, -ub));

            ref[j] = SCIPgetSolVal(scip, refsol, var);
            ref[j] = MIN(ub, MAX(lb, ref[j])); /* project value into bounds */
         }

         SCIP_CALL( generateCut(scip, conshdlr, cons, ref, sol, violside, row, efficacy, checkcurvmultivar, minefficacy) );
      }
   }
   /* gauge cut */
   if( mode == 'g' )
   {
      assert((consdata->isconvex && violside == SCIP_SIDETYPE_RIGHT) || (consdata->isconcave && violside == SCIP_SIDETYPE_LEFT));
      if( conshdlrdata->gaugecuts && consdata->isgaugeavailable )
      {
         SCIP_CALL( computeReferencePointGauge(scip, conshdlr, cons, refsol, ref, &success) );
      }
      if( success )
      {
         SCIP_CALL( generateCut(scip, conshdlr, cons, ref, sol, violside, row, efficacy, checkcurvmultivar, minefficacy) );
      }
   }
   /* projection cut */
   if( mode == 'p' )
   {
      assert((consdata->isconvex && violside == SCIP_SIDETYPE_RIGHT) || (consdata->isconcave && violside == SCIP_SIDETYPE_LEFT));
      if( conshdlrdata->projectedcuts && consdata->isedavailable )
      {
         SCIP_CALL( computeReferencePointProjection(scip, conshdlr, cons, refsol, ref) );
         SCIP_CALL( generateCut(scip, conshdlr, cons, ref, sol, violside, row, efficacy, checkcurvmultivar, minefficacy) );
      }
   }
   /* gradient linearization cut at refsol */
   if( mode == 'l' )
   {
      assert((consdata->isconvex && violside == SCIP_SIDETYPE_RIGHT) || (consdata->isconcave && violside == SCIP_SIDETYPE_LEFT));
      for( j = 0; j < consdata->nquadvars; ++j )
      {
         var = consdata->quadvarterms[j].var;
         lb  = SCIPvarGetLbLocal(var);
         ub  = SCIPvarGetUbLocal(var);
         /* do not like variables at infinity */
         assert(!SCIPisInfinity(scip,  lb));
         assert(!SCIPisInfinity(scip, -ub));

         ref[j] = SCIPgetSolVal(scip, refsol, var);
         ref[j] = MIN(ub, MAX(lb, ref[j])); /* project value into bounds */
      }
      SCIP_CALL( generateCut(scip, conshdlr, cons, ref, sol, violside, row, efficacy, checkcurvmultivar, minefficacy) );
   }

   SCIPfreeBufferArray(scip, &ref);

   return SCIP_OKAY;
}

/** tries to find a cut that intersects with an unbounded ray of the LP
 *
 *  For convex functions, we do this by linearizing in the feasible solution of the LPI.
 *  For nonconvex functions, we just call generateCutSol with the unbounded solution as reference point.
 */
static
SCIP_RETCODE generateCutUnboundedLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SIDETYPE         violside,           /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real*            rowrayprod,         /**< buffer to store product of ray with row coefficients, or NULL if not of interest */
   SCIP_Bool             checkcurvmultivar   /**< are we allowed to check the curvature of a multivariate quadratic function, if not done yet */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_BILINTERM* bilinterm;
   SCIP_VAR*  var;
   SCIP_Real* ref;
   SCIP_Real  matrixrayprod;
   SCIP_Real  linrayprod;
   SCIP_Real  quadrayprod;
   SCIP_Real  rayval;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(row  != NULL);
   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *row = NULL;

   if( !SCIPhasPrimalRay(scip) )
   {
      SCIPdebugMsg(scip, "do not have primal ray, thus cannot resolve unboundedness\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( checkCurvature(scip, cons, checkcurvmultivar) );
   if( (!consdata->isconvex && violside == SCIP_SIDETYPE_RIGHT) ||
      (!consdata->isconcave && violside == SCIP_SIDETYPE_LEFT) )
   {
      /* if not convex, just call generateCut and hope it's getting something useful */
      SCIP_CALL( generateCutSol(scip, conshdlr, cons, NULL, NULL, violside, row, NULL, FALSE, -SCIPinfinity(scip), 'd') );

      /* compute product of cut coefficients with ray, if required */
      if( *row != NULL && rowrayprod != NULL )
      {
         *rowrayprod = 0.0;
         for( i = 0; i < SCIProwGetNNonz(*row); ++i )
         {
            assert(SCIProwGetCols(*row)[i] != NULL);
            var = SCIPcolGetVar(SCIProwGetCols(*row)[i]);
            assert(var != NULL);

            *rowrayprod += SCIProwGetVals(*row)[i] * SCIPgetPrimalRayVal(scip, var);
         }
      }

      return SCIP_OKAY;
   }

   /* we seek for a linearization of the quadratic function such that it intersects with the unbounded ray
    * that is, we need a reference point ref such that for the gradient g of xAx+bx in ref, we have
    *   <g, ray> > 0.0 if rhs is finite and <g, ray> < 0.0 if lhs is finite
    * Since g = 2*A*ref + b, we have <g, ray> = <2*A*ref + b, ray> = <ref, 2*A*ray> + <b,ray>
    * initially, for finite rhs, we set ref_i = 1.0 if (A*ray)_i > 0.0 and ref_i = -1.0 if (A*ray)_i < 0.0 (for finite lhs analog)
    * <ref, 2*A*ray> + <b,ray> is sufficiently larger 0.0, we call generateCut for this point, otherwise, we scale up ref
    */

   quadrayprod = 0.0; /* <ref, 2*A*ray> */
   linrayprod = 0.0;  /* <b, ray> */
   SCIP_CALL( SCIPallocBufferArray(scip, &ref, consdata->nquadvars) );
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      var = consdata->quadvarterms[i].var;
      rayval = SCIPgetPrimalRayVal(scip, var);

      /* compute i-th entry of (2*A*ray) */
      matrixrayprod = 2.0 * consdata->quadvarterms[i].sqrcoef * rayval;
      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         bilinterm = &consdata->bilinterms[consdata->quadvarterms[i].adjbilin[j]];
         matrixrayprod += bilinterm->coef * SCIPgetPrimalRayVal(scip, bilinterm->var1 == var ? bilinterm->var2 : bilinterm->var1);
      }

      if( SCIPisPositive(scip, matrixrayprod) )
         ref[i] = (violside == SCIP_SIDETYPE_RIGHT ?  1.0 : -1.0);
      else if( SCIPisNegative(scip, matrixrayprod) )
         ref[i] = (violside == SCIP_SIDETYPE_RIGHT ? -1.0 :  1.0);
      else
         ref[i] = 0.0;

      quadrayprod += matrixrayprod * ref[i];
      linrayprod += consdata->quadvarterms[i].lincoef * rayval;
   }
   assert((violside == SCIP_SIDETYPE_RIGHT && quadrayprod >= 0.0) || (violside == SCIP_SIDETYPE_LEFT && quadrayprod <= 0.0));

   if( SCIPisZero(scip, quadrayprod) )
   {
      SCIPdebugMsg(scip, "ray is zero along cons <%s>\n", SCIPconsGetName(cons));
      SCIPfreeBufferArray(scip, &ref);
      return SCIP_OKAY;
   }

   /* add linear part to linrayprod */
   for( i = 0; i < consdata->nlinvars; ++i )
      linrayprod += consdata->lincoefs[i] * SCIPgetPrimalRayVal(scip, consdata->linvars[i]);

   SCIPdebugMsg(scip, "initially have <b,ray> = %g and <ref, 2*A*ref> = %g\n", linrayprod, quadrayprod);

   /* we scale the refpoint up, such that <ref, 2*A*ray> >= -2*<b, ray> (rhs finite) or <ref, 2*A*ray> <= -2*<b, ray> (lhs finite), if <b,ray> is not zero
    * if <b,ray> is zero, then we scale refpoint up if |<ref, 2*A*ray>| < 1.0
    */
   if( (!SCIPisZero(scip, linrayprod) && violside == SCIP_SIDETYPE_RIGHT && quadrayprod < -2*linrayprod) ||
      ( !SCIPisZero(scip, linrayprod) && violside == SCIP_SIDETYPE_LEFT  && quadrayprod > -2*linrayprod) ||
      (SCIPisZero(scip, linrayprod) && REALABS(quadrayprod) < 1.0) )
   {
      SCIP_Real scale;

      if( !SCIPisZero(scip, linrayprod) )
         scale = 2*REALABS(linrayprod/quadrayprod);  /*lint !e795 */
      else
         scale = 1.0/REALABS(quadrayprod);

      SCIPdebugMsg(scip, "scale refpoint by %g\n", scale);
      for( i = 0; i < consdata->nquadvars; ++i )
         ref[i] *= scale;
      quadrayprod *= scale;
   }

   if( rowrayprod != NULL )
      *rowrayprod = quadrayprod + linrayprod;

   SCIPdebugMsg(scip, "calling generateCut, expecting ray product %g\n", quadrayprod + linrayprod);
   SCIP_CALL( generateCut(scip, conshdlr, cons, ref, NULL, violside, row, NULL, FALSE, -SCIPinfinity(scip)) );

   SCIPfreeBufferArray(scip, &ref);

   return SCIP_OKAY;
}

/** processes a cut for constraint cons, i.e., checks numerics and possibly adds cut to sepastore */
static
SCIP_RETCODE processCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< cut to process */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             efficacy,           /**< efficacy of row in reference solution */
   SCIP_Real             minefficacy,        /**< minimal efficacy */
   SCIP_Bool             inenforcement,      /**< whether we are in constraint enforcement */
   SCIP_Real*            bestefficacy,       /**< buffer to store best efficacy of a cut that was added to the LP, if found; or NULL if not of interest */
   SCIP_RESULT*          result              /**< result of separation */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(row != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);
   assert(cons != NULL);

   /* no cut to process */
   if( *row == NULL )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPisGT(scip, efficacy, minefficacy) && SCIPisCutApplicable(scip, *row) )  /*lint !e644 */
   {
      SCIP_Bool infeasible;

      /* cut cuts off solution */
      SCIP_CALL( SCIPaddRow(scip, *row, FALSE /* forcecut */, &infeasible) );
      if( infeasible )
      {
         SCIPdebugMessage("cut for constraint <%s> is infeasible -> cutoff.\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
      }
      else
      {
         SCIPdebugMessage("add cut with efficacy %g for constraint <%s> violated by %g\n", efficacy,
               SCIPconsGetName(cons), consdata->lhsviol+consdata->rhsviol);
         *result = SCIP_SEPARATED;
      }
      SCIP_CALL( SCIPresetConsAge(scip, cons) );

      /* mark row as not removable from LP for current node, if in enforcement */
      if( inenforcement && !conshdlrdata->enfocutsremovable )
         SCIPmarkRowNotRemovableLocal(scip, *row);
   }
   if( bestefficacy != NULL && efficacy > *bestefficacy )
      *bestefficacy = efficacy;

   SCIP_CALL( SCIPreleaseRow (scip, row) );
   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 *
 *  assumes that constraint violations have been computed
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Real             minefficacy,        /**< minimal efficacy of a cut if it should be added to the LP */
   SCIP_Bool             inenforcement,      /**< whether we are in constraint enforcement */
   SCIP_RESULT*          result,             /**< result of separation */
   SCIP_Real*            bestefficacy        /**< buffer to store best efficacy of a cut that was added to the LP, if found; or NULL if not of interest */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          efficacy;
   SCIP_SIDETYPE      violside;
   int                c;
   SCIP_ROW*          row;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( bestefficacy != NULL )
      *bestefficacy = 0.0;

   row = NULL;
   /* loop over both sides of each constraint */
   for( c = 0, violside = SCIP_SIDETYPE_LEFT; c < nconss; c = (violside == SCIP_SIDETYPE_LEFT ? c : c+1), violside = (violside == SCIP_SIDETYPE_LEFT ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT) )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* if side not violated, then go on */
      if( !SCIPisGT(scip, violside == SCIP_SIDETYPE_LEFT ? consdata->lhsviol : consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      /* we are not feasible anymore */
      if( *result == SCIP_FEASIBLE )
         *result = SCIP_DIDNOTFIND;

      /* generate cut */
      if( sol == NULL && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
      {
         /* if the LP is unbounded, then we need a cut that cuts into the direction of a hopefully existing primal ray
          * that is, assume a ray r is given such that p + t*r is feasible for the LP for all t >= t_0 and some p
          * given a cut lhs <= <c,x> <= rhs, we check whether it imposes an upper bound on t and thus bounds the ray
          * this is given if rhs < infinity and <c,r> > 0, since we then enforce <c,p+t*r> = <c,p> + t<c,r> <= rhs, i.e., t <= (rhs - <c,p>)/<c,r>
          * similar, lhs > -infinity and <c,r> < 0 is good
          */
         SCIP_Real rayprod;

         rayprod = 0.0; /* for compiler */
         SCIP_CALL( generateCutUnboundedLP(scip, conshdlr, conss[c], violside, &row, &rayprod, conshdlrdata->checkcurvature) );

         if( row != NULL )
         {
            if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPisPositive(scip, rayprod) )
               efficacy =  rayprod;
            else if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) && SCIPisNegative(scip, rayprod) )
               efficacy = -rayprod;
            else
               efficacy = 0.0;

            SCIP_CALL( processCut(scip, &row, conshdlr, conss[c], efficacy, minefficacy, inenforcement, bestefficacy, result) );
         }
         continue;
      }
      else
      {
         SCIP_CALL( generateCutSol(scip, conshdlr, conss[c], sol, NULL, violside, &row, &efficacy,
            conshdlrdata->checkcurvature, minefficacy, 'd') );

         SCIP_CALL( processCut(scip, &row, conshdlr, conss[c], efficacy, minefficacy, inenforcement, bestefficacy, result) );
      }

      if( *result == SCIP_CUTOFF )
         break;

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */
      if( c >= nusefulconss && *result == SCIP_SEPARATED )
         break;
   }

   return SCIP_OKAY;
}

/** adds linearizations cuts for convex constraints w.r.t. a given reference point to cutpool and sepastore
 *
 *  - If separatedlpsol is not NULL, then a cut that separates the LP solution is added to the sepastore and is forced to enter the LP.
 *  - If separatedlpsol is not NULL, but cut does not separate the LP solution, then it is added to the cutpool only.
 *  - If separatedlpsol is NULL, then cut is added to cutpool only.
 */
static
SCIP_RETCODE addLinearizationCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             ref,                /**< reference point where to linearize, or NULL for LP solution */
   SCIP_Bool*            separatedlpsol,     /**< buffer to store whether a cut that separates the current LP solution was found and added to LP,
                                              *   or NULL if adding to cutpool only */
   SCIP_Real             minefficacy         /**< minimal efficacy of a cut when checking for separation of LP solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool addedtolp;
   SCIP_ROW* row;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( separatedlpsol != NULL )
      *separatedlpsol = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613 */

      if( SCIPconsIsLocal(conss[c]) || !SCIPconsIsEnabled(conss[c]) )  /*lint !e613 */
         continue;

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
      assert(consdata != NULL);

      if( consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIP_CALL( generateCutSol(scip, conshdlr, conss[c], NULL, ref, SCIP_SIDETYPE_RIGHT, &row, NULL,
               conshdlrdata->checkcurvature, -SCIPinfinity(scip), 'l') );  /*lint !e613 */
      }
      else if( consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         SCIP_CALL( generateCutSol(scip, conshdlr, conss[c], NULL, ref, SCIP_SIDETYPE_LEFT,  &row, NULL,
               conshdlrdata->checkcurvature, -SCIPinfinity(scip), 'l') );  /*lint !e613 */
      }
      else
         continue;

      if( row == NULL )
         continue;

      addedtolp = FALSE;

      /* if caller wants, then check if cut separates LP solution and add to sepastore if so */
      if( separatedlpsol != NULL )
      {
         SCIP_Real efficacy;

         efficacy = -SCIPgetRowLPFeasibility(scip, row);
         if( efficacy >= minefficacy )
         {
            SCIP_Bool infeasible;

            *separatedlpsol = TRUE;
            addedtolp = TRUE;
            SCIP_CALL( SCIPaddRow(scip, row, TRUE, &infeasible) );
            assert( ! infeasible );
            SCIPdebugMsg(scip, "added linearization cut <%s> to LP, efficacy = %g\n", SCIProwGetName(row), efficacy);
         }
      }

      if( !SCIProwIsLocal(row) && !addedtolp )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
         SCIPdebugMsg(scip, "added linearization cut <%s> to cutpool\n", SCIProwGetName(row));
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS**    conss;
   int            nconss;
   SCIP_SOL*      sol;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) != 0);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   nconss = SCIPconshdlrGetNConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* we are only interested in solution coming from some heuristic other than trysol, but not from the tree
    * the reason for ignoring trysol solutions is that they may come from an NLP solve in sepalp, where we already added linearizations,
    * or are from the tree, but postprocessed via proposeFeasibleSolution
    */
   if( SCIPsolGetHeur(sol) == NULL || SCIPsolGetHeur(sol) == conshdlrdata->trysolheur )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   SCIPdebugMsg(scip, "caught new sol event %" SCIP_EVENTTYPE_FORMAT " from heur <%s>; have %d conss\n", SCIPeventGetType(event), SCIPheurGetName(SCIPsolGetHeur(sol)), nconss);

   SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, sol, NULL, 0.0) );

   return SCIP_OKAY;
}

/** registers branching candidates according to convexification gap rule
 *
 * That is, computes for every nonconvex term the gap between the terms value in the LP solution and the value of the underestimator
 * as it would be (and maybe has been) constructed by the separation routines of this constraint handler. Then it registers all
 * variables occurring in each term with the computed gap. If variables appear in more than one term, they are registered several times.
 */
static
SCIP_RETCODE registerBranchingCandidatesGap(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
   int                j;
   SCIP_Bool          xbinary;
   SCIP_Bool          ybinary;
   SCIP_Bool          xunbounded;
   SCIP_Bool          yunbounded;
   SCIP_VAR*          x;
   SCIP_VAR*          y;
   SCIP_Real          xlb;
   SCIP_Real          xub;
   SCIP_Real          xval;
   SCIP_Real          ylb;
   SCIP_Real          yub;
   SCIP_Real          yval;
   SCIP_Real          gap;
   SCIP_Real          coef_;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   *nnotify = 0;
   yval = SCIP_INVALID;
   xval = SCIP_INVALID;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlr != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->nquadvars )
         continue;

      if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || consdata->isconcave) &&
         ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || consdata->isconvex ) )
         continue;
      SCIPdebugMsg(scip, "cons <%s> violation: %g %g  convex: %u %u\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->isconvex, consdata->isconcave);

      /* square terms */
      for( j = 0; j < consdata->nquadvars; ++j )
      {
         x = consdata->quadvarterms[j].var;
         if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && consdata->quadvarterms[j].sqrcoef < 0) ||
            ( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && consdata->quadvarterms[j].sqrcoef > 0) )
         {
            xlb = SCIPvarGetLbLocal(x);
            xub = SCIPvarGetUbLocal(x);
            if( SCIPisRelEQ(scip, xlb, xub) )
            {
               SCIPdebugMsg(scip, "ignore fixed variable <%s>[%g, %g], diff %g\n", SCIPvarGetName(x), xlb, xub, xub-xlb);
               continue;
            }

            xval = SCIPgetSolVal(scip, sol, x);

            /* if variable is at bounds, then no need to branch, since secant is exact there */
            if( SCIPisLE(scip, xval, xlb) || SCIPisGE(scip, xval, xub) )
               continue;

            if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
               gap = SCIPinfinity(scip);
            else
               gap = (xval-xlb)*(xub-xval)/(1+2*ABS(xval));
            assert(!SCIPisFeasNegative(scip, gap));
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, MAX(gap, 0.0), SCIP_INVALID) );
            ++*nnotify;
         }
      }

      /* bilinear terms */
      for( j = 0; j < consdata->nbilinterms; ++j )
      {
         /* if any of the variables if fixed, then it actually behaves like a linear term, so we don't need to branch on it */
         x = consdata->bilinterms[j].var1;
         xlb = SCIPvarGetLbLocal(x);
         xub = SCIPvarGetUbLocal(x);
         if( SCIPisRelEQ(scip, xlb, xub) )
            continue;

         y = consdata->bilinterms[j].var2;
         ylb = SCIPvarGetLbLocal(y);
         yub = SCIPvarGetUbLocal(y);
         if( SCIPisRelEQ(scip, ylb, yub) )
            continue;

         xunbounded = SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub);
         yunbounded = SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub);

         /* compute gap, if both variable are bounded */
         gap = SCIPinfinity(scip);
         if( !xunbounded && !yunbounded )
         {
            xval = SCIPgetSolVal(scip, sol, x);
            yval = SCIPgetSolVal(scip, sol, y);

            /* if both variables are at one of its bounds, then no need to branch, since McCormick is exact there */
            if( (SCIPisLE(scip, xval, xlb) || SCIPisGE(scip, xval, xub)) &&
               ( SCIPisLE(scip, yval, ylb) || SCIPisGE(scip, yval, yub)) )
               continue;

            xval = MAX(xlb, MIN(xval, xub));
            yval = MAX(ylb, MIN(yval, yub));

            coef_ = SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) ? -consdata->bilinterms[j].coef : consdata->bilinterms[j].coef;
            if( coef_ > 0.0 )
            {
               if( (xub-xlb)*yval + (yub-ylb)*xval <= xub*yub - xlb*ylb )
                  gap = (xval*yval - xlb*yval - ylb*xval + xlb*ylb) / (1+sqrt(xval*xval + yval*yval));
               else
                  gap = (xval*yval - xval*yub - yval*xub + xub*yub) / (1+sqrt(xval*xval + yval*yval));
            }
            else
            { /* coef_ < 0 */
               if( (xub-xlb)*yval - (yub-ylb)*xval <= xub*ylb - xlb*yub )
                  gap = -(xval*yval - xval*ylb - yval*xub + xub*ylb) / (1+sqrt(xval*xval + yval*yval));
               else
                  gap = -(xval*yval - xval*yub - yval*xlb + xlb*yub) / (1+sqrt(xval*xval + yval*yval));
            }

            assert(!SCIPisNegative(scip, gap / MAX3(MAX(REALABS(xlb), REALABS(xub)), MAX(REALABS(ylb), REALABS(yub)), 1.0)));  /*lint !e666*/
            if( gap < 0.0 )
               gap = 0.0;

            /* use tighter relaxation when using linear inequalities to adjust the branching scores for bilinear terms */
            if( consdata->bilintermsidx != NULL && conshdlrdata->usebilinineqbranch )
            {
               BILINESTIMATOR* bilinestimator;
               int bilinidx;

               assert(conshdlrdata->bilinestimators != NULL);

               bilinidx = consdata->bilintermsidx[j];
               assert(bilinidx >= 0 && bilinidx < conshdlrdata->nbilinterms);

               bilinestimator = &conshdlrdata->bilinestimators[bilinidx];
               assert(bilinestimator != NULL);
               assert(bilinestimator->x == x);
               assert(bilinestimator->y == y);

               if( SCIPisGT(scip, bilinestimator->lastimprfac, 0.0) )
                  gap *= MAX(0.0, 1.0 - bilinestimator->lastimprfac);
            }
         }

         /* if one of the variables is binary or integral with domain width 1, then branching on this makes the term linear, so prefer this */
         xbinary = SCIPvarIsBinary(x) || (SCIPvarIsIntegral(x) && xub - xlb < 1.5);
         ybinary = SCIPvarIsBinary(y) || (SCIPvarIsIntegral(y) && yub - ylb < 1.5);
         if( xbinary )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( ybinary )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( xbinary || ybinary )
            continue;

         /* if one of the variables is unbounded, then branch on it first */
         if( xunbounded )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( yunbounded )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( xunbounded || yunbounded )
            continue;

         /* if both variables are integral, prefer the one with the smaller domain, so variable gets fixed soon
          * does not seem to work well on tln instances, so disable for now and may look at it later again
          */
#ifdef BRANCHTOLINEARITY
         if( SCIPvarIsIntegral(x) && SCIPvarIsIntegral(y) )
         {
            if( SCIPisLT(scip, xub-xlb, yub-ylb) )
            {
               SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
               ++*nnotify;
               continue;
            }
            if( SCIPisGT(scip, xub-xlb, yub-ylb) )
            {
               SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
               ++*nnotify;
               continue;
            }
         }
#endif

         /* in the regular case, suggest those variables which are not at its bounds for branching
          * this is, because after branching both variables will be one the bounds, and McCormick will be exact then */
         if( !SCIPisLE(scip, xval, xlb) && !SCIPisGE(scip, xval, xub) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, x, gap, SCIP_INVALID) );
            ++*nnotify;
         }
         if( !SCIPisLE(scip, yval, ylb) && !SCIPisGE(scip, yval, yub) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, y, gap, SCIP_INVALID) );
            ++*nnotify;
         }
      }
   }

   SCIPdebugMsg(scip, "registered %d branching candidates\n", *nnotify);

   return SCIP_OKAY;
}

/** registers branching candidates according to constraint violation rule
 *
 * That is, registers all variables appearing in nonconvex terms^1 with a score that is the violation of the constraint.
 * This is the same rule as is applied in cons_nonlinear and other nonlinear constraint handlers.
 *
 * 1) We mean all quadratic variables that appear either in a nonconvex square term or in a bilinear term, if the constraint
 * itself is nonconvex. (and this under the assumption that the rhs is violated; for violated lhs, swap terms)
 */
static
SCIP_RETCODE registerBranchingCandidatesViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA*     consdata;
   SCIP_QUADVARTERM*  quadvarterm;
   int                c;
   int                j;
   SCIP_VAR*          x;
   SCIP_Real          xlb;
   SCIP_Real          xub;
   SCIP_Real          xval;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->nquadvars )
         continue;

      if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || consdata->isconcave) &&
         ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || consdata->isconvex ) )
         continue;
      SCIPdebugMsg(scip, "cons %s violation: %g %g  convex: %u %u\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->isconvex, consdata->isconcave);

      for( j = 0; j < consdata->nquadvars; ++j )
      {
         quadvarterm = &consdata->quadvarterms[j];
         if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && quadvarterm->sqrcoef < 0) ||
             (SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && quadvarterm->sqrcoef > 0) ||
             quadvarterm->nadjbilin > 0 )
         {
            x = quadvarterm->var;
            xlb = SCIPvarGetLbLocal(x);
            xub = SCIPvarGetUbLocal(x);

            if( quadvarterm->nadjbilin == 0 )
            {
               xval = SCIPgetSolVal(scip, sol, x);

               /* if variable is at bounds and only in a nonconvex square term, then no need to branch, since secant is exact there */
               if( SCIPisLE(scip, xval, xlb) || SCIPisGE(scip, xval, xub) )
                  continue;
            }

            if( SCIPisRelEQ(scip, xlb, xub) )
            {
               SCIPdebugMsg(scip, "ignore fixed variable <%s>[%g, %g], diff %g\n", SCIPvarGetName(x), xlb, xub, xub-xlb);
               continue;
            }

            SCIP_CALL( SCIPaddExternBranchCand(scip, x, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++*nnotify;
         }
      }
   }

   SCIPdebugMsg(scip, "registered %d branching candidates\n", *nnotify);

   return SCIP_OKAY;
}

/** registers branching candidates according to centrality rule
 *
 * That is, registers all variables appearing in nonconvex terms^1 with a score that is given by the distance of the
 * variable value from its bounds. This rule should not make sense, as the distance to the bounds is also (often) considered
 * by the branching rule later on.
 *
 * 1) We mean all quadratic variables that appear either in a nonconvex square term or in a bilinear term, if the constraint
 * itself is nonconvex. (and this under the assumption that the rhs is violated; for violated lhs, swap terms)
 */
static
SCIP_RETCODE registerBranchingCandidatesCentrality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA*     consdata;
   SCIP_QUADVARTERM*  quadvarterm;
   int                c;
   int                j;
   SCIP_VAR*          x;
   SCIP_Real          xlb;
   SCIP_Real          xub;
   SCIP_Real          xval;
   SCIP_Real          score;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->nquadvars )
         continue;

      if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || consdata->isconcave) &&
         ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || consdata->isconvex ) )
         continue;
      SCIPdebugMsg(scip, "cons %s violation: %g %g  convex: %u %u\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->isconvex, consdata->isconcave);

      for( j = 0; j < consdata->nquadvars; ++j )
      {
         quadvarterm = &consdata->quadvarterms[j];
         if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && quadvarterm->sqrcoef < 0) ||
             (SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && quadvarterm->sqrcoef > 0) ||
             quadvarterm->nadjbilin > 0 )
         {
            x = quadvarterm->var;
            xlb = SCIPvarGetLbLocal(x);
            xub = SCIPvarGetUbLocal(x);

            if( SCIPisRelEQ(scip, xlb, xub) )
            {
               SCIPdebugMsg(scip, "ignore fixed variable <%s>[%g, %g], diff %g\n", SCIPvarGetName(x), xlb, xub, xub-xlb);
               continue;
            }

            xval = SCIPgetSolVal(scip, sol, x);
            xval = MAX(xlb, MIN(xub, xval));

            /* compute relative difference of xval to each of its bounds
             * and scale such that if xval were in the middle, we get a score of 1
             * and if xval is on one its bounds, the score is 0
             */
            if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
            {
               if( (!SCIPisInfinity(scip, -xlb) && SCIPisEQ(scip, xval, xlb)) || (!SCIPisInfinity(scip, xub) && SCIPisEQ(scip, xval, xub)) )
                  score = 0.0;
               else
                  score = 1.0;
            }
            else
            {
               score = 4.0 * (xval - xlb) * (xub - xval) / ((xub - xlb) * (xub - xlb));
            }

            SCIP_CALL( SCIPaddExternBranchCand(scip, x, score, SCIP_INVALID) );
            ++*nnotify;
         }
      }
   }

   SCIPdebugMsg(scip, "registered %d branching candidates\n", *nnotify);

   return SCIP_OKAY;
}

/** registers branching candidates */
static
SCIP_RETCODE registerBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   switch( conshdlrdata->branchscoring )
   {
      case 'g' :
         SCIP_CALL( registerBranchingCandidatesGap(scip, conshdlr, conss, nconss, sol, nnotify) );
         break;

      case 'v' :
         SCIP_CALL( registerBranchingCandidatesViolation(scip, conshdlr, conss, nconss, sol, nnotify) );
         break;

      case 'c' :
         SCIP_CALL( registerBranchingCandidatesCentrality(scip, conshdlr, conss, nconss, sol, nnotify) );
         break;

      default :
         SCIPerrorMessage("invalid branchscoring selection");
         SCIPABORT();
         return SCIP_ERROR; /*lint !e527*/
   }

   return SCIP_OKAY;
}


/** registers a quadratic variable from a violated constraint as branching candidate that has a large absolute value in the (LP) relaxation */
static
SCIP_RETCODE registerLargeRelaxValueVariableForBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_VAR**            brvar               /**< buffer to store branching variable */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_Real           val;
   SCIP_Real           brvarval;
   int                 i;
   int                 c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   *brvar = NULL;
   brvarval = -1.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         /* do not propose fixed variables */
         if( SCIPisRelEQ(scip, SCIPvarGetLbLocal(consdata->quadvarterms[i].var), SCIPvarGetUbLocal(consdata->quadvarterms[i].var)) )
            continue;
         val = SCIPgetSolVal(scip, sol, consdata->quadvarterms[i].var);
         if( ABS(val) > brvarval )
         {
            brvarval = ABS(val);
            *brvar = consdata->quadvarterms[i].var;
         }
      }
   }

   if( *brvar != NULL )
   {
      SCIP_CALL( SCIPaddExternBranchCand(scip, *brvar, brvarval, SCIP_INVALID) );
   }

   return SCIP_OKAY;
}

/** replaces violated quadratic constraints where all quadratic variables are fixed by linear constraints */
static
SCIP_RETCODE replaceByLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            addedcons,          /**< buffer to store whether a linear constraint was added */
   SCIP_Bool*            reduceddom,         /**< whether a domain has been reduced */
   SCIP_Bool*            infeasible          /**< whether we detected infeasibility */
   )
{
   SCIP_CONS*          cons;
   SCIP_CONSDATA*      consdata;
   SCIP_RESULT         checkresult;
   SCIP_VAR*           var;
   SCIP_Bool           tightened;
   SCIP_Real           constant;
   SCIP_Real           val1;
   SCIP_Real           val2;
   int                 i;
   int                 c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);
   assert(addedcons != NULL);
   assert(reduceddom != NULL);
   assert(infeasible != NULL);

   *addedcons = FALSE;
   *reduceddom = FALSE;
   *infeasible = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      constant = 0.0;

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         var = consdata->quadvarterms[i].var;

         /* variables should be fixed if constraint is violated */
         assert(SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

         val1 = (SCIPvarGetUbLocal(var) + SCIPvarGetLbLocal(var)) / 2.0;
         constant += (consdata->quadvarterms[i].lincoef + consdata->quadvarterms[i].sqrcoef * val1) * val1;

         SCIPdebugMessage("<%s>: [%.20g, %.20g]\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

         /* if variable is not fixed w.r.t. absolute eps yet, then try to fix it
          * (SCIPfixVar() doesn't allow for small tightenings, so tighten lower and upper bound separately)
          */
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, val1, TRUE, infeasible, &tightened) );
            if( *infeasible )
            {
               SCIPdebugMsg(scip, "Fixing almost fixed variable <%s> lead to infeasibility.\n", SCIPvarGetName(var));
               return SCIP_OKAY;
            }
            if( tightened )
            {
               SCIPdebugMsg(scip, "Tightened lower bound of almost fixed variable <%s>.\n", SCIPvarGetName(var));
               *reduceddom = TRUE;
            }

            SCIP_CALL( SCIPtightenVarUb(scip, var, val1, TRUE, infeasible, &tightened) );
            if( *infeasible )
            {
               SCIPdebugMsg(scip, "Fixing almost fixed variable <%s> lead to infeasibility.\n", SCIPvarGetName(var));
               return SCIP_OKAY;
            }
            if( tightened )
            {
               SCIPdebugMsg(scip, "Tightened upper bound of almost fixed variable <%s>.\n", SCIPvarGetName(var));
               *reduceddom = TRUE;
            }
         }
      }

      /* if some quadratic variable was fixed now, then restart node (next enfo round) */
      if( *reduceddom )
         return SCIP_OKAY;

      for( i = 0; i < consdata->nbilinterms; ++i )
      {
         val1 = (SCIPvarGetUbLocal(consdata->bilinterms[i].var1) + SCIPvarGetLbLocal(consdata->bilinterms[i].var1)) / 2.0;
         val2 = (SCIPvarGetUbLocal(consdata->bilinterms[i].var2) + SCIPvarGetLbLocal(consdata->bilinterms[i].var2)) / 2.0;
         constant += consdata->bilinterms[i].coef * val1 * val2;
      }

      /* check if we have a bound change */
      if ( consdata->nlinvars == 1 )
      {
         SCIP_Real coef;
         SCIP_Real lhs;
         SCIP_Real rhs;

         coef = *consdata->lincoefs;
         var = *consdata->linvars;

         assert( ! SCIPisZero(scip, coef) );

         /* compute lhs/rhs, divide already by |coef| */
         if ( SCIPisInfinity(scip, -consdata->lhs) )
            lhs = -SCIPinfinity(scip);
         else
            lhs = (consdata->lhs - constant) / REALABS(coef);

         if ( SCIPisInfinity(scip, consdata->rhs) )
            rhs = SCIPinfinity(scip);
         else
            rhs = (consdata->rhs - constant) / REALABS(coef);

         SCIPdebugMsg(scip, "Linear constraint with one variable: %.20g <= %g <%s> <= %.20g\n", lhs, coef > 0.0 ? 1.0 : -1.0, SCIPvarGetName(var), rhs);

         SCIPdebugMessage("<%s>: [%.20g, %.20g]\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

         if ( coef < 0.0 )
         {
            /* swap lhs and rhs, with negated sign */
            SCIP_Real h;
            h = rhs;
            rhs = -lhs;
            lhs = -h;
         }
         SCIPdebugMsg(scip, "Linear constraint is a bound: %.20g <= <%s> <= %.20g\n", lhs, SCIPvarGetName(var), rhs);

         if( SCIPisInfinity(scip, -rhs) || SCIPisInfinity(scip, lhs) )
         {
            SCIPdebugMsg(scip, "node will marked as infeasible since lb/ub of %s is +/-infinity\n",
               SCIPvarGetName(var));

            *infeasible = TRUE;
            return SCIP_OKAY;
         }

         if ( ! SCIPisInfinity(scip, -lhs) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, lhs, TRUE, infeasible, &tightened) );
            if ( *infeasible )
            {
               SCIPdebugMsg(scip, "Lower bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               SCIPdebugMsg(scip, "Lower bound changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }

         if ( ! SCIPisInfinity(scip, rhs) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, rhs, TRUE, infeasible, &tightened) );
            if ( *infeasible )
            {
               SCIPdebugMsg(scip, "Upper bound leads to infeasibility.\n");
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               SCIPdebugMsg(scip, "Upper bound changed.\n");
               *reduceddom = TRUE;
               return SCIP_OKAY;
            }
         }
      }
      else
      {
         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, SCIPconsGetName(conss[c]),
               consdata->nlinvars, consdata->linvars, consdata->lincoefs,
               (SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : (consdata->lhs - constant)),
               (SCIPisInfinity(scip,  consdata->rhs) ?  SCIPinfinity(scip) : (consdata->rhs - constant)),
               SCIPconsIsInitial(conss[c]), SCIPconsIsSeparated(conss[c]), SCIPconsIsEnforced(conss[c]),
               SCIPconsIsChecked(conss[c]), SCIPconsIsPropagated(conss[c]),  TRUE,
               SCIPconsIsModifiable(conss[c]), SCIPconsIsDynamic(conss[c]), SCIPconsIsRemovable(conss[c]),
               SCIPconsIsStickingAtNode(conss[c])) );

         SCIPdebugMsg(scip, "replace quadratic constraint <%s> by linear constraint after all quadratic vars have been fixed\n", SCIPconsGetName(conss[c]) );
         SCIPdebugPrintCons(scip, cons, NULL);

         SCIP_CALL( SCIPcheckCons(scip, cons, NULL, FALSE, FALSE, FALSE, &checkresult) );

         if( checkresult != SCIP_INFEASIBLE && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIPdebugMsg(scip, "linear constraint is feasible and LP optimal, thus do not add\n");
         }
         else
         {
            SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
            *addedcons = TRUE;
         }
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
      SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
   }

   return SCIP_OKAY;
}

/** tightens a lower bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new lower bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds             /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(intervalinfty > 0.0);
   assert(bnd > -intervalinfty);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   /* new bound is no improvement */
   if( SCIPisHugeValue(scip, -bnd) || SCIPisLE(scip, bnd, SCIPvarGetLbLocal(var)) )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   /* new lower bound is very low (between -intervalinfty and -SCIPinfinity()) */
   if( SCIPisInfinity(scip, -bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarLb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarLb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMsg(scip, "%s found constraint <%s> infeasible due to tightened lower bound %g for variable <%s>\n",
         SCIPinProbing(scip) ? "in probing" : "", SCIPconsGetName(cons), bnd, SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMsg(scip, "%s tightened lower bound of variable <%s> in constraint <%s> to %g\n",
         SCIPinProbing(scip) ? "in probing" : "", SCIPvarGetName(var), SCIPconsGetName(cons), bnd);
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** tightens an upper bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new upper bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds             /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(intervalinfty > 0.0);
   assert(bnd < intervalinfty);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   /* new bound is no improvement */
   if( SCIPisHugeValue(scip, bnd) || SCIPisGE(scip, bnd, SCIPvarGetUbLocal(var)) )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, -bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   /* new upper bound is very high (between SCIPinfinity() and intervalinfty) */
   if( SCIPisInfinity(scip, bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarUb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarUb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMsg(scip, "%s found constraint <%s> infeasible due to tightened upper bound %g for variable <%s>\n",
         SCIPinProbing(scip) ? "in probing" : "", SCIPconsGetName(cons), bnd, SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMsg(scip, "%s tightened upper bound of variable <%s> in constraint <%s> to %g\n",
         SCIPinProbing(scip) ? "in probing" : "", SCIPvarGetName(var), SCIPconsGetName(cons), bnd);
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** solves a quadratic equation \f$ a x^2 + b x \in rhs \f$ (with b an interval) and reduces bounds on x or deduces infeasibility if possible */
static
SCIP_RETCODE propagateBoundsQuadVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where we currently propagate */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             var,                /**< variable which bounds with might tighten */
   SCIP_Real             a,                  /**< coefficient in square term */
   SCIP_INTERVAL         b,                  /**< coefficient in linear term */
   SCIP_INTERVAL         rhs,                /**< right hand side of quadratic equation */
   SCIP_RESULT*          result,             /**< result of propagation */
   int*                  nchgbds             /**< buffer where to add number of tightened bounds */
   )
{
   SCIP_INTERVAL newrange;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   /* compute solution of a*x^2 + b*x \in rhs */
   if( a == 0.0 && SCIPintervalGetInf(b) == 0.0 && SCIPintervalGetSup(b) == 0.0 )
   {
      /* relatively easy case: 0.0 \in rhs, thus check if infeasible or just redundant */
      if( SCIPintervalGetInf(rhs) > 0.0 || SCIPintervalGetSup(rhs) < 0.0 )
      {
         SCIPdebugMsg(scip, "found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *result = SCIP_CUTOFF;
      }
      return SCIP_OKAY;
   }
   else if( SCIPvarGetLbLocal(var) >= 0.0 )
   {
      SCIP_INTERVAL a_;

      /* need only positive solutions */
      SCIPintervalSet(&a_, a);
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &newrange, a_, b, rhs);
   }
   else if( SCIPvarGetUbLocal(var) <= 0.0 )
   {
      /* need only negative solutions */
      SCIP_INTERVAL a_;
      SCIP_INTERVAL tmp;
      SCIPintervalSet(&a_, a);
      SCIPintervalSetBounds(&tmp, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));
      SCIPintervalSolveUnivariateQuadExpressionPositive(intervalinfty, &tmp, a_, tmp, rhs);
      if( SCIPintervalIsEmpty(intervalinfty, tmp) )
      {
         SCIPdebugMsg(scip, "found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(var));
         *result = SCIP_CUTOFF;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         return SCIP_OKAY;
      }
      SCIPintervalSetBounds(&newrange, -SCIPintervalGetSup(tmp), -SCIPintervalGetInf(tmp));
   }
   else
   {
      /* need both positive and negative solution */
      SCIP_INTERVAL a_;
      SCIPintervalSet(&a_, a);
      SCIPintervalSolveUnivariateQuadExpression(intervalinfty, &newrange, a_, b, rhs);
   }

   /* SCIPdebugMsg(scip, "%g x^2 + [%g, %g] x in [%g, %g] -> [%g, %g]\n", a, b.inf, b.sup, rhs.inf, rhs.sup, newrange.inf, newrange.sup); */

   if( SCIPisInfinity(scip, SCIPintervalGetInf(newrange)) || SCIPisInfinity(scip, -SCIPintervalGetSup(newrange)) )
   {
      /* domain outside [-infty, +infty] -> declare node infeasible */
      SCIPdebugMsg(scip, "found <%s> infeasible because propagated domain of quadratic variable <%s> is outside of (-infty, +infty)\n",
         SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      return SCIP_OKAY;
   }

   if( SCIPintervalIsEmpty(intervalinfty, newrange) )
   {
      SCIPdebugMsg(scip, "found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(var));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -SCIPintervalGetInf(newrange)) )
   {
      SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, SCIPintervalGetInf(newrange), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip,  SCIPintervalGetSup(newrange)) )
   {
      SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, SCIPintervalGetSup(newrange), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/* The new version below computes potentially tighter bounds, but also always adds a small safety area since it is not implemented roundingsafe.
 * This may be a reason why it gives worse results on one of two instances.
 * Further, I have only very few instances where one can expect a difference.
 */
#ifndef PROPBILINNEW
/** tries to deduce domain reductions for x in xsqrcoef x^2 + xlincoef x + ysqrcoef y^2 + ylincoef y + bilincoef x y \\in rhs
 *
 *  @note Domain reductions for y are not deduced.
 */
static
SCIP_RETCODE propagateBoundsBilinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint, where the bilinear term belongs to */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_Real             xsqrcoef,           /**< square coefficient of x */
   SCIP_Real             xlincoef,           /**< linear coefficient of x */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             ysqrcoef,           /**< square coefficient of y */
   SCIP_Real             ylincoef,           /**< linear coefficient of y */
   SCIP_Real             bilincoef,          /**< bilinear coefficient of x*y */
   SCIP_INTERVAL         rhs,                /**< right hand side of quadratic equation */
   SCIP_RESULT*          result,             /**< pointer to store result of domain propagation */
   int*                  nchgbds             /**< counter to increment if domain reductions are found */
   )
{
   SCIP_INTERVAL myrhs;
   SCIP_INTERVAL varbnds;
   SCIP_INTERVAL lincoef;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);
   assert(bilincoef != 0.0);

   if( SCIPintervalIsEntire(intervalinfty, rhs) )
      return SCIP_OKAY;

   /* try to find domain reductions for x */
   SCIPintervalSetBounds(&varbnds, MIN(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y)), MAX(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y)));  /*lint !e666 */

   /* put ysqrcoef*y^2 + ylincoef * y into rhs */
   if( SCIPintervalGetSup(rhs) >= intervalinfty )
   {
      /* if rhs is unbounded by above, it is sufficient to get an upper bound on ysqrcoef*y^2 + ylincoef * y */
      SCIP_ROUNDMODE roundmode;
      SCIP_Real      tmp;

      SCIPintervalSet(&lincoef, ylincoef);
      tmp = SCIPintervalQuadUpperBound(intervalinfty, ysqrcoef, lincoef, varbnds);
      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();
      SCIPintervalSetBounds(&myrhs, SCIPintervalGetInf(rhs) - tmp, intervalinfty);
      SCIPintervalSetRoundingMode(roundmode);
   }
   else if( SCIPintervalGetInf(rhs) <= -intervalinfty )
   {
      /* if rhs is unbounded by below, it is sufficient to get a  lower bound on ysqrcoef*y^2 + ylincoef * y */
      SCIP_ROUNDMODE roundmode;
      SCIP_Real      tmp;

      SCIPintervalSet(&lincoef, -ylincoef);
      tmp = -SCIPintervalQuadUpperBound(intervalinfty, -ysqrcoef, lincoef, varbnds);
      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();
      SCIPintervalSetBounds(&myrhs, -intervalinfty, SCIPintervalGetSup(rhs) - tmp);
      SCIPintervalSetRoundingMode(roundmode);
   }
   else
   {
      /* if rhs is bounded, we need both bounds on ysqrcoef*y^2 + ylincoef * y */
      SCIP_INTERVAL tmp;

      SCIPintervalSet(&lincoef, ylincoef);
      SCIPintervalQuad(intervalinfty, &tmp, ysqrcoef, lincoef, varbnds);
      SCIPintervalSub(intervalinfty, &myrhs, rhs, tmp);
   }

   /* create equation xsqrcoef * x^2 + (xlincoef + bilincoef * [ylb, yub]) * x \in myrhs */
   SCIPintervalMulScalar(intervalinfty, &lincoef, varbnds, bilincoef);
   SCIPintervalAddScalar(intervalinfty, &lincoef, lincoef, xlincoef);

   /* propagate bounds on x */
   SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, x, xsqrcoef, lincoef, myrhs, result, nchgbds) );

   return SCIP_OKAY;
}
#else
/** tries to deduce domain reductions for x in xsqrcoef x^2 + xlincoef x + ysqrcoef y^2 + ylincoef y + bilincoef x y \\in rhs
 *
 *  @note Domain reductions for y are not deduced.
 */
static
SCIP_RETCODE propagateBoundsBilinearTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint, where the bilinear term belongs to */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_Real             xsqrcoef,           /**< square coefficient of x */
   SCIP_Real             xlincoef,           /**< linear coefficient of x */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             ysqrcoef,           /**< square coefficient of y */
   SCIP_Real             ylincoef,           /**< linear coefficient of y */
   SCIP_Real             bilincoef,          /**< bilinear coefficient of x*y */
   SCIP_INTERVAL         rhs,                /**< right hand side of quadratic equation */
   SCIP_RESULT*          result,             /**< pointer to store result of domain propagation */
   int*                  nchgbds             /**< counter to increment if domain reductions are found */
   )
{
   SCIP_INTERVAL xbnds;
   SCIP_INTERVAL ybnds;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(x != y);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);
   assert(bilincoef != 0.0);

   if( SCIPintervalIsEntire(intervalinfty, rhs) )
      return SCIP_OKAY;

   SCIPintervalSetBounds(&xbnds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x))),   /*lint !e666*/
      +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x))));  /*lint !e666*/
   SCIPintervalSetBounds(&ybnds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y))),   /*lint !e666*/
      +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y))));  /*lint !e666*/

   /* try to find domain reductions for x */
   SCIPintervalSolveBivariateQuadExpressionAllScalar(intervalinfty, &xbnds, xsqrcoef, ysqrcoef, bilincoef, xlincoef, ylincoef, rhs, xbnds, ybnds);

   if( SCIPintervalIsEmpty(intervalinfty, xbnds) )
   {
      SCIPdebugMsg(scip, "found <%s> infeasible due to domain propagation for quadratic variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(x));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip, -SCIPintervalGetInf(xbnds)) )
   {
      SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, x, SCIPintervalGetInf(xbnds), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   if( !SCIPisInfinity(scip,  SCIPintervalGetSup(xbnds)) )
   {
      SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, x, SCIPintervalGetSup(xbnds), result, nchgbds) );
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}
#endif

/** computes the minimal and maximal activity for the quadratic part in a constraint data
 *
 *  Only sums up terms that contribute finite values.
 *  Gives the number of terms that contribute infinite values.
 *  Only computes those activities where the corresponding side of the constraint is finite.
 */
static
void propagateBoundsGetQuadActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             intervalinfty,      /**< infinity value used in interval operations */
   SCIP_Real*            minquadactivity,    /**< minimal activity of quadratic variable terms where only terms with finite minimal activity contribute */
   SCIP_Real*            maxquadactivity,    /**< maximal activity of quadratic variable terms where only terms with finite maximal activity contribute */
   int*                  minactivityinf,     /**< number of quadratic variables that contribute -infinity to minimal activity */
   int*                  maxactivityinf,     /**< number of quadratic variables that contribute +infinity to maximal activity */
   SCIP_INTERVAL*        quadactcontr        /**< contribution of each quadratic variables to quadactivity */
   )
{  /*lint --e{666}*/
   SCIP_ROUNDMODE prevroundmode;
   int       i;
   int       j;
   int       k;
   SCIP_INTERVAL tmp;
   SCIP_Real bnd;
   SCIP_INTERVAL xrng;
   SCIP_INTERVAL lincoef;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(minquadactivity != NULL);
   assert(maxquadactivity != NULL);
   assert(minactivityinf != NULL);
   assert(maxactivityinf != NULL);
   assert(quadactcontr != NULL);

   /* if lhs is -infinite, then we do not compute a maximal activity, so we set it to  infinity
    * if rhs is  infinite, then we do not compute a minimal activity, so we set it to -infinity
    */
   *minquadactivity = SCIPisInfinity(scip,  consdata->rhs) ? -intervalinfty : 0.0;
   *maxquadactivity = SCIPisInfinity(scip, -consdata->lhs) ?  intervalinfty : 0.0;

   *minactivityinf = 0;
   *maxactivityinf = 0;

   if( consdata->nquadvars == 0 )
   {
      SCIPintervalSet(&consdata->quadactivitybounds, 0.0);
      return;
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      /* there should be no quadratic variables fixed at -/+ infinity due to our locks */
      assert(!SCIPisInfinity(scip,  SCIPvarGetLbLocal(consdata->quadvarterms[i].var)));
      assert(!SCIPisInfinity(scip, -SCIPvarGetUbLocal(consdata->quadvarterms[i].var)));

      SCIPintervalSetBounds(&quadactcontr[i], -intervalinfty, intervalinfty);

      SCIPintervalSetBounds(&xrng,
         -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->quadvarterms[i].var), SCIPvarGetUbLocal(consdata->quadvarterms[i].var))),
         +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->quadvarterms[i].var), SCIPvarGetUbLocal(consdata->quadvarterms[i].var))));

      SCIPintervalSet(&lincoef, consdata->quadvarterms[i].lincoef);
      for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
      {
         k = consdata->quadvarterms[i].adjbilin[j];
         if( consdata->bilinterms[k].var1 != consdata->quadvarterms[i].var )
            continue; /* handle this term later */

         SCIPintervalSetBounds(&tmp,
            -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
            +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));
         SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);
         SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);
      }

      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         /* compute maximal activity only if there is a finite left hand side */
         bnd = SCIPintervalQuadUpperBound(intervalinfty, consdata->quadvarterms[i].sqrcoef, lincoef, xrng);
         if( bnd >= intervalinfty )
         {
            ++*maxactivityinf;
         }
         else if( SCIPisInfinity(scip, -bnd) )
         {
            /* if maximal activity is below value for -infinity, let's take -1e10 as upper bound on maximal activity
             * @todo Something better?
             */
            bnd = -sqrt(SCIPinfinity(scip));
            *maxquadactivity += bnd;
            quadactcontr[i].sup = bnd;
         }
         else
         {
            prevroundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeUpwards();
            *maxquadactivity += bnd;
            SCIPintervalSetRoundingMode(prevroundmode);
            quadactcontr[i].sup = bnd;
         }
      }

      if( !SCIPisInfinity(scip,  consdata->rhs) )
      {
         /* compute minimal activity only if there is a finite right hand side */
         SCIPintervalSetBounds(&lincoef, -SCIPintervalGetSup(lincoef), -SCIPintervalGetInf(lincoef));
         bnd = -SCIPintervalQuadUpperBound(intervalinfty, -consdata->quadvarterms[i].sqrcoef, lincoef, xrng);

         if( bnd <= -intervalinfty )
         {
            ++*minactivityinf;
         }
         else if( SCIPisInfinity(scip, bnd) )
         {
            /* if minimal activity is above value for infinity, let's take 1e10 as lower bound on minimal activity
             * @todo Something better?
             */
            bnd = sqrt(SCIPinfinity(scip));
            *minquadactivity += bnd;
            quadactcontr[i].inf = bnd;
         }
         else
         {
            prevroundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeDownwards();
            *minquadactivity += bnd;
            SCIPintervalSetRoundingMode(prevroundmode);
            quadactcontr[i].inf = bnd;
         }
      }

   }

   SCIPintervalSetBounds(&consdata->quadactivitybounds,
      (*minactivityinf > 0 ? -intervalinfty : *minquadactivity),
      (*maxactivityinf > 0 ?  intervalinfty : *maxquadactivity));
   assert(!SCIPintervalIsEmpty(intervalinfty, consdata->quadactivitybounds));
}

/** propagates bounds on a quadratic constraint */
static
SCIP_RETCODE propagateBoundsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation call */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   SCIP_Bool*            redundant           /**< buffer where to store whether constraint has been found to be redundant */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA*     consdata;
   SCIP_INTERVAL      consbounds;    /* lower and upper bounds of constraint */
   SCIP_INTERVAL      consactivity;  /* activity of linear plus quadratic part */
   SCIP_Real          intervalinfty; /* infinity used for interval computation */
   SCIP_Real          minquadactivity; /* lower bound on finite activities of quadratic part */
   SCIP_Real          maxquadactivity; /* upper bound on finite activities of quadratic part */
   int                quadminactinf; /* number of quadratic variables that contribute -infinity to minimal activity of quadratic term */
   int                quadmaxactinf; /* number of quadratic variables that contribute +infinity to maximal activity of quadratic term */
   SCIP_INTERVAL*     quadactcontr;  /* contribution of each quadratic variable term to quadactivity */

   SCIP_VAR*          var;
   SCIP_INTERVAL      rhs;           /* right hand side of quadratic equation */
   SCIP_INTERVAL      tmp;
   SCIP_ROUNDMODE     roundmode;
   SCIP_Real          bnd;
   int                i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(redundant != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN;
   *redundant = FALSE;

   *result = SCIP_DIDNOTFIND;

   intervalinfty = 1000 * SCIPinfinity(scip) * SCIPinfinity(scip);

   quadactcontr = NULL;
   quadminactinf = -1;
   quadmaxactinf = -1;

   SCIPdebugMsg(scip, "start domain propagation for constraint <%s>\n", SCIPconsGetName(cons));

   /* make sure we have activity of linear term and that they are consistent */
   consdataUpdateLinearActivity(scip, consdata, intervalinfty);
   assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
   assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
   assert(consdata->minlinactivityinf >= 0);
   assert(consdata->maxlinactivityinf >= 0);

   /* sort quadratic variable terms, in case we need to search for terms occuring in bilinear terms later
    * we sort here already, since we rely on a constant variable order during this method
    */
   if( consdata->nbilinterms > 0 )
   {
      SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );
   }

   /* compute activity of quad term part, if not up to date
    * in that case, we also collect the contribution of each quad var term for later */
   if( SCIPintervalIsEmpty(intervalinfty, consdata->quadactivitybounds) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &quadactcontr, consdata->nquadvars) );
      propagateBoundsGetQuadActivity(scip, consdata, intervalinfty, &minquadactivity, &maxquadactivity, &quadminactinf, &quadmaxactinf, quadactcontr);
      assert(!SCIPintervalIsEmpty(intervalinfty, consdata->quadactivitybounds));
   }

   SCIPdebugMsg(scip, "linear activity: [%g, %g]   quadratic activity: [%g, %g]\n",
      (consdata->minlinactivityinf > 0 ? -intervalinfty : consdata->minlinactivity),
      (consdata->maxlinactivityinf > 0 ?  intervalinfty : consdata->maxlinactivity),
      consdata->quadactivitybounds.inf, consdata->quadactivitybounds.sup);

   /* extend constraint bounds by epsilon to avoid some numerical difficulties */
   SCIPintervalSetBounds(&consbounds,
      -infty2infty(SCIPinfinity(scip), intervalinfty, -consdata->lhs+SCIPepsilon(scip)),
      +infty2infty(SCIPinfinity(scip), intervalinfty,  consdata->rhs+SCIPepsilon(scip)));

   /* check redundancy and infeasibility */
   SCIPintervalSetBounds(&consactivity, consdata->minlinactivityinf > 0 ? -intervalinfty : consdata->minlinactivity,
      consdata->maxlinactivityinf > 0 ? intervalinfty : consdata->maxlinactivity);
   SCIPintervalAdd(intervalinfty, &consactivity, consactivity, consdata->quadactivitybounds);
   if( SCIPintervalIsSubsetEQ(intervalinfty, consactivity, consbounds) )
   {
      SCIPdebugMsg(scip, "found constraint <%s> to be redundant: sides: [%g, %g], activity: [%g, %g]\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity));
      *redundant = TRUE;
      goto CLEANUP;
   }

   /* was SCIPintervalAreDisjoint(consbounds, consactivity), but that would allow violations up to eps only
    * we need to decide feasibility w.r.t. feastol (but still want to propagate w.r.t. eps)
    */
   if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisGT(scip, consdata->lhs-SCIPfeastol(scip), SCIPintervalGetSup(consactivity))) ||
       (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisLT(scip, consdata->rhs+SCIPfeastol(scip), SCIPintervalGetInf(consactivity))) )
   {
      SCIPdebugMsg(scip, "found constraint <%s> to be infeasible; sides: [%g, %g], activity: [%g, %g], infeas: %g\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity),
         MAX(consdata->lhs - SCIPintervalGetSup(consactivity), SCIPintervalGetInf(consactivity) - consdata->rhs));
      *result = SCIP_CUTOFF;
      goto CLEANUP;
   }

   /* propagate linear part \in rhs = consbounds - quadactivity (use the one from consdata, since that includes infinities) */
   SCIPintervalSub(intervalinfty, &rhs, consbounds, consdata->quadactivitybounds);
   if( !SCIPintervalIsEntire(intervalinfty, rhs) )
   {
      SCIP_Real coef;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         coef = consdata->lincoefs[i];
         var  = consdata->linvars[i];

         /* skip fixed variables
          * @todo is that a good or a bad idea?
          *   we can't expect much more tightening, but may detect infeasiblity, but shouldn't the check on the constraints activity detect that?
          */
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            continue;

         /* due to large variable bounds and large coefficients, it might happen that the activity of the linear part
          * exceeds +/-SCIPinfinity() after updating the activities in consdataUpdateLinearActivity{Lb,Ub}Change; in
          * order to detect this case we need to check whether the value of consdata->{min,max}linactivity is infinite
          * (see #1433)
          */
         if( coef > 0.0 )
         {
            if( SCIPintervalGetSup(rhs) < intervalinfty )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the upper bound on var x */
               if( consdata->minlinactivityinf == 0 && !SCIPisInfinity(scip, -consdata->minlinactivity) )
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
                  /* tighten upper bound on x to (rhs.sup - (minlinactivity - coef * xlb)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd += coef * SCIPvarGetLbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the minimal linear activity equal -infinity, so
                   * we tighten upper bound on x to just (rhs.sup - minlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetInf(rhs) > -intervalinfty )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the lower bound on var x */
               if( consdata->maxlinactivityinf == 0 && !SCIPisInfinity(scip, consdata->maxlinactivity) )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* tighten lower bound on x to (rhs.inf - (maxlinactivity - coef * xub)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd += coef * SCIPvarGetUbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (rhs.inf - maxlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is +infinity and x is not solely responsible for this */
            }
         }
         else
         {
            assert(coef < 0.0 );
            if( SCIPintervalGetInf(rhs) > -intervalinfty )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the upper bound on var x */
               if( consdata->maxlinactivityinf == 0  && !SCIPisInfinity(scip, consdata->maxlinactivity) )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetLbLocal(var)));
                  /* compute upper bound on x to (maxlinactivity - coef * xlb) - rhs.inf / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd += (-coef) * SCIPvarGetLbLocal(var);
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (maxlinactivity - rhs.inf) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetSup(rhs) < intervalinfty )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
               /* try to tighten the lower bound on var x */
               if( consdata->minlinactivityinf == 0 && !SCIPisInfinity(scip, -consdata->minlinactivity) )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* compute lower bound on x to (minlinactivity - coef * xub) - rhs.sup / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd += (-coef) * SCIPvarGetUbLocal(var);
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal -infinity, so
                   * we tighten lower bound on x to just (minlinactivity - rhs.sup) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, intervalinfty, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }
         }
      }
      if( *result == SCIP_CUTOFF )
         goto CLEANUP;
   }

   /* propagate quadratic part \in rhs = consbounds - linactivity */
   assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777 */
   assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777 */
   consdataUpdateLinearActivity(scip, consdata, intervalinfty); /* make sure, activities of linear part did not become invalid by above bound changes, if any */
   assert(consdata->minlinactivityinf > 0 || consdata->maxlinactivityinf > 0 || consdata->minlinactivity <= consdata->maxlinactivity);
   SCIPintervalSetBounds(&tmp,
      (consdata->minlinactivityinf > 0 ? -intervalinfty : consdata->minlinactivity),
      (consdata->maxlinactivityinf > 0 ?  intervalinfty : consdata->maxlinactivity));
   SCIPintervalSub(intervalinfty, &rhs, consbounds, tmp);
   if( !SCIPintervalIsEntire(intervalinfty, rhs) )
   {
      if( consdata->nquadvars == 1 )
      {
         /* quadratic part is just a*x^2+b*x -> a common case that we treat directly */
         SCIP_INTERVAL lincoef;    /* linear coefficient of quadratic equation */

         assert(consdata->nbilinterms == 0);

         var = consdata->quadvarterms[0].var;
         SCIPintervalSet(&lincoef, consdata->quadvarterms[0].lincoef);

         /* propagate a*x^2 + b*x \in rhs */
         SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, consdata->quadvarterms[0].sqrcoef, lincoef, rhs, result, nchgbds) );
      }
      else if( consdata->nbilinterms == 1 && consdata->nquadvars == 2 )
      {
         /* quadratic part is just ax*x^2+bx*x + ay*y^2+by*y + c*xy -> a common case that we treat directly */
         assert(consdata->bilinterms[0].var1 == consdata->quadvarterms[0].var || consdata->bilinterms[0].var1 == consdata->quadvarterms[1].var);
         assert(consdata->bilinterms[0].var2 == consdata->quadvarterms[0].var || consdata->bilinterms[0].var2 == consdata->quadvarterms[1].var);

         /* find domain reductions for x from a_x x^2 + b_x x + a_y y^2 + b_y y + c x y \in rhs */
         SCIP_CALL( propagateBoundsBilinearTerm(scip, cons, intervalinfty,
               consdata->quadvarterms[0].var, consdata->quadvarterms[0].sqrcoef, consdata->quadvarterms[0].lincoef,
               consdata->quadvarterms[1].var, consdata->quadvarterms[1].sqrcoef, consdata->quadvarterms[1].lincoef,
               consdata->bilinterms[0].coef,
               rhs, result, nchgbds) );
         if( *result != SCIP_CUTOFF )
         {
            /* find domain reductions for y from a_x x^2 + b_x x + a_y y^2 + b_y y + c x y \in rhs */
            SCIP_CALL( propagateBoundsBilinearTerm(scip, cons, intervalinfty,
                  consdata->quadvarterms[1].var, consdata->quadvarterms[1].sqrcoef, consdata->quadvarterms[1].lincoef,
                  consdata->quadvarterms[0].var, consdata->quadvarterms[0].sqrcoef, consdata->quadvarterms[0].lincoef,
                  consdata->bilinterms[0].coef,
                  rhs, result, nchgbds) );
         }
      }
      else
      {
         /* general case */

         /* compute "advanced" information on quad var term activities, if not up-to-date */
         if( quadminactinf == -1  )
         {
            assert(quadactcontr == NULL);
            SCIP_CALL( SCIPallocBufferArray(scip, &quadactcontr, consdata->nquadvars) );
            propagateBoundsGetQuadActivity(scip, consdata, intervalinfty, &minquadactivity, &maxquadactivity, &quadminactinf, &quadmaxactinf, quadactcontr);
         }
         assert(quadactcontr != NULL);
         assert(quadminactinf >= 0);
         assert(quadmaxactinf >= 0);

         /* if the quad activities are not hopelessly unbounded on useful sides, try to deduce domain reductions on quad vars */
         if( (SCIPintervalGetSup(rhs) <  intervalinfty && quadminactinf <= 1) ||
            ( SCIPintervalGetInf(rhs) > -intervalinfty && quadmaxactinf <= 1) )
         {
            SCIP_INTERVAL lincoef;
            SCIP_INTERVAL rhs2;
            int j;
            int k;

            for( i = 0; i < consdata->nquadvars; ++i )
            {
               var = consdata->quadvarterms[i].var;

               /* skip fixed variables
                * @todo is that a good or a bad idea?
                *   we can't expect much more tightening, but may detect infeasiblity, but shouldn't the check on the constraints activity detect that?
                */
               if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                  continue;

               /* compute rhs2 such that we can propagate quadvarterm(x_i) \in rhs2 */

               /* setup rhs2.sup = rhs.sup - (quadactivity.inf - quadactcontr[i].inf), if everything were finite
                * if only quadactcontr[i].inf is infinite (i.e., the other i are all finite), we just get rhs2.sup = rhs.sup
                * otherwise we get rhs2.sup = infinity */
               if( SCIPintervalGetSup(rhs) < intervalinfty )
               {
                  if( quadminactinf == 0 || (quadminactinf == 1 && SCIPintervalGetInf(quadactcontr[i]) <= -intervalinfty) )
                  {
                     roundmode = SCIPintervalGetRoundingMode();
                     SCIPintervalSetRoundingModeUpwards();
                     rhs2.sup = rhs.sup - minquadactivity;  /*lint !e644*/
                     /* if the residual quad min activity w.r.t. quad var term i is finite and nonzero, so add it to right hand side */
                     if( quadminactinf == 0 && SCIPintervalGetInf(quadactcontr[i]) != 0.0 )
                        rhs2.sup += SCIPintervalGetInf(quadactcontr[i]);
                     SCIPintervalSetRoundingMode(roundmode);
                  }
                  else
                  {
                     /* there are either >= 2 quad var terms contributing -infinity, or there is one which is not i */
                     rhs2.sup = intervalinfty;
                  }
               }
               else
               {
                  rhs2.sup = intervalinfty;
               }

               /* setup rhs2.inf = rhs.inf - (quadactivity.sup - quadactcontr[i].sup), see also above */
               if( SCIPintervalGetInf(rhs) > -intervalinfty )
               {
                  if( quadmaxactinf == 0 || (quadmaxactinf == 1 && SCIPintervalGetSup(quadactcontr[i]) >= intervalinfty) )
                  {
                     roundmode = SCIPintervalGetRoundingMode();
                     SCIPintervalSetRoundingModeDownwards();
                     rhs2.inf = rhs.inf - maxquadactivity;  /*lint !e644*/
                     /* if the residual quad max activity w.r.t. quad var term i is finite and nonzero, so add it to right hand side */
                     if( quadmaxactinf == 0 && SCIPintervalGetSup(quadactcontr[i]) != 0.0 )
                        rhs2.inf += SCIPintervalGetSup(quadactcontr[i]);
                     SCIPintervalSetRoundingMode(roundmode);
                  }
                  else
                  {
                     /* there are either >= 2 quad var terms contributing infinity, or there is one which is not i */
                     rhs2.inf = -intervalinfty;
                  }
               }
               else
               {
                  rhs2.inf = -intervalinfty;
               }
               assert(!SCIPintervalIsEmpty(intervalinfty, rhs2));

               /* if rhs2 is entire, then there is nothing we could propagate */
               if( SCIPintervalIsEntire(intervalinfty, rhs2) )
                  continue;

               /* assemble linear coefficient for quad equation a*x^2 + b*x \in rhs2 */
               SCIPintervalSet(&lincoef, consdata->quadvarterms[i].lincoef);
               for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
               {
                  k = consdata->quadvarterms[i].adjbilin[j];
#if 1
                  if( consdata->bilinterms[k].var1 == var )
                  {
                     /* bilinear term k contributes to the activity of quad var term i, so just add bounds to linear coef */
                     SCIPintervalSetBounds(&tmp,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
                        +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));
                     SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);
                     SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);
                  }
                  else
                  {
                     /* bilinear term k does not contribute to the activity of quad var term i
                      * so bounds on term k are contained in rhs2
                      * if they are finite, we try to remove them from rhs2 and update lincoef instead
                      * if the bounds on bilinear term k as added to rhs2 are old due to recent bound tightening, we may not do best possible, but still correct
                      * HOWEVER: when computing rhs2, we may not just have added the bounds for the bilinear term, but for the associated quadratic term
                      *   for this complete term, we used SCIPintervalQuad to compute the bounds
                      *   since we do not want to repeat a call to SCIPintervalQuad for that quadratic term with bilinear term k removed,
                      *   we only remove the bounds for the bilinear term k from rhs2 if the associated quadratic term consists only of this bilinear term,
                      *   i.e., the quadratic term corresponding to var1 should be only var1*var2, but have no square or linear coefs or other bilinear terms
                      *   (for efficiency reasons, we check here only if there are any other bilinear terms than var1*var2 associated with var1, even if they are not associated with the quad var term for var1)
                      */
                     SCIP_INTERVAL me;
                     SCIP_INTERVAL bilinbounds;
                     int otherpos;

                     assert(consdata->bilinterms[k].var2 == var);

                     assert(consdata->quadvarssorted);
                     SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->bilinterms[k].var1, &otherpos) );
                     assert(otherpos >= 0);
                     assert(consdata->quadvarterms[otherpos].var == consdata->bilinterms[k].var1);

                     if( (consdata->quadvarterms[otherpos].sqrcoef != 0.0) || consdata->quadvarterms[otherpos].lincoef != 0.0 ||
                          consdata->quadvarterms[otherpos].nadjbilin > 1 )
                        continue;

                     /* set tmp to bounds of other variable and multiply with bilin coef */
                     SCIPintervalSetBounds(&tmp,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var1), SCIPvarGetUbLocal(consdata->bilinterms[k].var1))),
                        +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var1), SCIPvarGetUbLocal(consdata->bilinterms[k].var1))));
                     SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);

                     /* set me to bounds of i'th variable */
                     SCIPintervalSetBounds(&me,
                        -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
                        +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));

                     /* remove me*tmp from rhs2 */

                     roundmode = SCIPintervalGetRoundingMode();

                     if( rhs2.inf > -intervalinfty )
                     {
                        /* need upward rounding for SCIPintervalMulSup */
                        SCIPintervalSetRoundingModeUpwards();
                        SCIPintervalMulSup(intervalinfty, &bilinbounds, me, tmp);
                        /* rhs2.inf += bilinbounds.sup, but we are in upward rounding */
                        if( bilinbounds.sup < intervalinfty )
                           rhs2.inf = SCIPintervalNegateReal(SCIPintervalNegateReal(rhs2.inf) - bilinbounds.sup);
                     }

                     if( rhs2.sup <  intervalinfty )
                     {
                        /* need downward rounding for SCIPintervalMulInf */
                        SCIPintervalSetRoundingModeDownwards();
                        SCIPintervalMulInf(intervalinfty, &bilinbounds, me, tmp);
                        /* rhs2.sup += bilinbounds.inf, but we are in downward rounding */
                        if( bilinbounds.inf > -intervalinfty )
                           rhs2.sup = SCIPintervalNegateReal(SCIPintervalNegateReal(rhs2.sup) - bilinbounds.inf);
                     }

                     SCIPintervalSetRoundingMode(roundmode);

                     /* add tmp to lincoef */
                     SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);
                  }
#else
                  if( consdata->bilinterms[k].var1 != var )
                     continue; /* this term does not contribute to the activity of quad var term i */

                  SCIPintervalSetBounds(&tmp,
                     -infty2infty(SCIPinfinity(scip), intervalinfty, -MIN(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))),
                     +infty2infty(SCIPinfinity(scip), intervalinfty,  MAX(SCIPvarGetLbLocal(consdata->bilinterms[k].var2), SCIPvarGetUbLocal(consdata->bilinterms[k].var2))));
                  SCIPintervalMulScalar(intervalinfty, &tmp, tmp, consdata->bilinterms[k].coef);
                  SCIPintervalAdd(intervalinfty, &lincoef, lincoef, tmp);
#endif
               }

               /* deduce domain reductions for x_i */
               SCIP_CALL( propagateBoundsQuadVar(scip, cons, intervalinfty, var, consdata->quadvarterms[i].sqrcoef, lincoef, rhs2, result, nchgbds) );
               if( *result == SCIP_CUTOFF )
                  goto CLEANUP;
            }
         }
      }
   }

 CLEANUP:
   SCIPfreeBufferArrayNull(scip, &quadactcontr);

   return SCIP_OKAY;
}

/** calls domain propagation for a set of constraints */
static
SCIP_RETCODE propagateBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation calls */
   int*                  nchgbds             /**< buffer where to add the the number of changed bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT propresult;
   SCIP_Bool   redundant;
   int         c;
   int         roundnr;
   SCIP_Bool   success;
   int         maxproprounds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);
   assert(nchgbds != NULL);

   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;
   roundnr = 0;
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      maxproprounds = conshdlrdata->maxproproundspresolve;
   else
      maxproprounds = conshdlrdata->maxproprounds;

   do
   {
      success = FALSE;
      ++roundnr;

      SCIPdebugMsg(scip, "starting domain propagation round %d of %d for %d constraints\n", roundnr, maxproprounds, nconss);

      for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
      {
         assert(conss != NULL);
         if( !SCIPconsIsEnabled(conss[c]) )
            continue;

         if( SCIPconsIsMarkedPropagate(conss[c]) )
         {
            /* unmark constraint for propagation */
            SCIP_CALL( SCIPunmarkConsPropagate(scip, conss[c]) );

            SCIP_CALL( propagateBoundsCons(scip, conshdlr, conss[c], &propresult, nchgbds, &redundant) );
            if( propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN )
            {
               *result = propresult;
               success = TRUE;
            }
            if( redundant )
            {
               SCIPdebugMsg(scip, "deleting constraint <%s> locally\n", SCIPconsGetName(conss[c]));
               SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
            }
         }
      }

   }
   while( success && *result != SCIP_CUTOFF && roundnr < maxproprounds );

   return SCIP_OKAY;
}

/** checks for a linear variable that can be increase or decreased without harming feasibility */
static
void consdataFindUnlockedLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int i;
   int poslock;
   int neglock;

   consdata->linvar_maydecrease = -1;
   consdata->linvar_mayincrease = -1;

   /* check for a linear variable that can be increase or decreased without harming feasibility */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      /* compute locks of i'th linear variable */
      assert(consdata->lincoefs[i] != 0.0);
      if( consdata->lincoefs[i] > 0.0 )
      {
         poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
      }
      else
      {
         poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - neglock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can decrease x without harming other constraints */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_maydecrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_maydecrease]) / consdata->lincoefs[consdata->linvar_maydecrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_maydecrease = i;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - poslock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can increase x without harm */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_mayincrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_mayincrease]) / consdata->lincoefs[consdata->linvar_mayincrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_mayincrease = i;
      }
   }

#ifdef SCIP_DEBUG
   if( consdata->linvar_mayincrease >= 0 )
   {
      SCIPdebugMsg(scip, "may increase <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_mayincrease]));
   }
   if( consdata->linvar_maydecrease >= 0 )
   {
      SCIPdebugMsg(scip, "may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_maydecrease]));
   }
#endif
}

/** Given a solution where every quadratic constraint is either feasible or can be made feasible by
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
   SCIP_CONSDATA* consdata;
   SCIP_SOL* newsol;
   SCIP_VAR* var;
   int c;
   SCIP_Real viol;
   SCIP_Real delta;
   SCIP_Real gap;
   SCIP_Bool solviolbounds;

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
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* recompute violation of solution in case solution has changed
       * get absolution violation and sign
       * @todo do this only if solution has changed
       */
      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conss[c], newsol, &solviolbounds) );  /*lint !e613*/
         assert(!solviolbounds);
         viol = consdata->lhs - consdata->activity;
      }
      else if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conss[c], newsol, &solviolbounds) );  /*lint !e613*/
         assert(!solviolbounds);
         viol = consdata->rhs - consdata->activity;
      }
      else
         continue; /* constraint is satisfied */

      assert(viol != 0.0);
      if( consdata->linvar_mayincrease >= 0 &&
         ((viol > 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) || (viol < 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0)) )
      {
         /* have variable where increasing makes the constraint less violated */
         var = consdata->linvars[consdata->linvar_mayincrease];
         /* compute how much we would like to increase var */
         delta = viol / consdata->lincoefs[consdata->linvar_mayincrease];
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
            /*lint --e{613} */
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy lhs-violation %g of cons <%s>\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->lincoefs[consdata->linvar_mayincrease] * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->linvar_maydecrease >= 0 &&
         ((viol > 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) || (viol < 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0)) )
      {
         /* have variable where decreasing makes constraint less violated */
         var = consdata->linvars[consdata->linvar_maydecrease];
         /* compute how much we would like to decrease var */
         delta = viol / consdata->lincoefs[consdata->linvar_maydecrease];
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
            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            /*lint --e{613} */
            SCIPdebugMsg(scip, "increase <%s> by %g to %g to remedy rhs-violation %g of cons <%s>\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var), viol, SCIPconsGetName(conss[c]));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->lincoefs[consdata->linvar_maydecrease] * delta;
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

/** helper function to enforce constraints */
static
SCIP_RETCODE enforceConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*             sol,                /**< solution to enforce (NULL for the LP solution) */
   SCIP_Bool             solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Real          maxviol;
   SCIP_RESULT        propresult;
   SCIP_RESULT        separateresult;
   int                nchgbds;
   int                nnotify;
   SCIP_Real          sepaefficacy;
   SCIP_Bool          solviolbounds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(nusefulconss >= 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, sol, &solviolbounds, &maxviolcon) );

   if( maxviolcon == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   if( solviolbounds )
   {
      /* if LP solution violates variable bounds, then this should be because a row was added that
       * introduced this variable newly to the LP, in which case it gets value 0.0; the row should
       * have been added to resolve an infeasibility, so solinfeasible should be TRUE
       * see also issue #627
       */
      assert(solinfeasible);
      /* however, if solinfeasible is actually not TRUE, then better cut off the node to avoid that SCIP
       * stops because infeasible cannot be resolved */
      /*lint --e{774} */
      if( !solinfeasible )
         *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   consdata = SCIPconsGetData(maxviolcon);
   assert(consdata != NULL);
   maxviol = consdata->lhsviol + consdata->rhsviol;
   assert(SCIPisGT(scip, maxviol, SCIPfeastol(scip)));

   SCIPdebugMsg(scip, "enforcement with max violation %g in cons <%s> for %s solution\n", maxviol, SCIPconsGetName(maxviolcon),
         sol == NULL ? "LP" : "relaxation");

   /* if we are above the 100'th enforcement round for this node, something is strange
    * (maybe the LP / relaxator does not think that the cuts we add are violated, or we do ECP on a high-dimensional convex function)
    * in this case, check if some limit is hit or SCIP should stop for some other reason and terminate enforcement by creating a dummy node
    * (in optimized more, returning SCIP_INFEASIBLE in *result would be sufficient, but in debug mode this would give an assert in scip.c)
    * the reason to wait for 100 rounds is to avoid calls to SCIPisStopped in normal runs, which may be expensive
    * we only increment nenforounds until 101 to avoid an overflow
    */
   if( conshdlrdata->lastenfonode == SCIPgetCurrentNode(scip) )
   {
      if( conshdlrdata->nenforounds > 100 )
      {
         if( SCIPisStopped(scip) )
         {
            SCIP_NODE* child;

            SCIP_CALL( SCIPcreateChild(scip, &child, 1.0, SCIPnodeGetEstimate(SCIPgetCurrentNode(scip))) );
            *result = SCIP_BRANCHED;

            return SCIP_OKAY;
         }
      }

      ++conshdlrdata->nenforounds;

      /* cut off the current subtree, if a limit on the enforcement rounds should be applied. At this point, feasible
       * solutions might get cut off; the enfolplimit parameter should therefore only be set if SCIP is used as a
       * heuristic solver and when the returned result (infeasible, optimal, the gap) can be ignored
       */
      if( conshdlrdata->enfolplimit != -1 && conshdlrdata->nenforounds > conshdlrdata->enfolplimit )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "cut off subtree because enforcement limit was reached; this might lead to incorrect results\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   else
   {
      conshdlrdata->lastenfonode = SCIPgetCurrentNode(scip);
      conshdlrdata->nenforounds = 0;
   }

   /* run domain propagation */
   nchgbds = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &nchgbds) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      SCIPdebugMsg(scip, "propagation succeeded (%s)\n", propresult == SCIP_CUTOFF ? "cutoff" : "reduceddom");
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we would like a cut that is efficient enough that it is not redundant in the LP (>lpfeastol)
    * however, we also don't want very weak cuts, so try to reach at least feastol (=lpfeastol by default, though)
    */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, SCIPfeastol(scip), TRUE, &separateresult, &sepaefficacy) );
   if( separateresult == SCIP_CUTOFF )
   {
      SCIPdebugMsg(scip, "separation found cutoff\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   if( separateresult == SCIP_SEPARATED )
   {
      SCIPdebugMsg(scip, "separation succeeded (bestefficacy = %g, minefficacy = %g)\n", sepaefficacy, SCIPfeastol(scip));
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* we are not feasible, the whole node is not infeasible, and we cannot find a good cut
    * -> collect variables for branching
    */

   SCIPdebugMsg(scip, "separation failed (bestefficacy = %g < %g = minefficacy ); max viol: %g\n", sepaefficacy, SCIPfeastol(scip), maxviol);

   /* find branching candidates */
   SCIP_CALL( registerBranchingCandidates(scip, conshdlr, conss, nconss, sol, &nnotify) );

   if( nnotify == 0 && !solinfeasible && SCIPfeastol(scip) > SCIPlpfeastol(scip) )
   {
      /* fallback 1: we also have no branching candidates, so try to find a weak cut */
      SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, SCIPlpfeastol(scip), TRUE, &separateresult, &sepaefficacy) );
      if( separateresult == SCIP_CUTOFF )
      {
         SCIPdebugMsg(scip, "separation found cutoff\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( separateresult == SCIP_SEPARATED )
      {
         SCIPdebugMsg(scip, "separation fallback succeeded, efficacy = %g\n", sepaefficacy);
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   }

   if( nnotify == 0 && !solinfeasible )
   {
      /* fallback 2: separation probably failed because of numerical difficulties with a convex constraint;
       *  if noone declared solution infeasible yet and we had not even found a weak cut, try to resolve by branching
       */
      SCIP_VAR* brvar = NULL;
      SCIP_CALL( registerLargeRelaxValueVariableForBranching(scip, conss, nconss, sol, &brvar) );
      if( brvar == NULL )
      {
         /* fallback 3: all quadratic variables seem to be fixed -> replace by linear constraint */
         SCIP_Bool addedcons;
         SCIP_Bool reduceddom;
         SCIP_Bool infeasible;

         SCIP_CALL( replaceByLinearConstraints(scip, conss, nconss, &addedcons, &reduceddom, &infeasible) );
         /* if the linear constraints are actually feasible, then adding them and returning SCIP_CONSADDED confuses SCIP
          * when it enforces the new constraints again and nothing resolves the infeasibility that we declare here
          * thus, we only add them if considered violated, and otherwise claim the solution is feasible (but print a
          * warning) */
         if ( infeasible )
            *result = SCIP_CUTOFF;
         else if ( addedcons )
            *result = SCIP_CONSADDED;
         else if ( reduceddom )
            *result = SCIP_REDUCEDDOM;
         else
         {
            *result = SCIP_FEASIBLE;
            SCIPwarningMessage(scip, "could not enforce feasibility by separating or branching; declaring solution with viol %g as feasible\n", maxviol);
            assert(!SCIPisInfinity(scip, maxviol));
         }
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMsg(scip, "Could not find any usual branching variable candidate. Proposed variable <%s> with LP value %g for branching.\n",
            SCIPvarGetName(brvar), SCIPgetSolVal(scip, sol, brvar));
         nnotify = 1;
      }
   }

   assert(*result == SCIP_INFEASIBLE && (solinfeasible || nnotify > 0));
   return SCIP_OKAY;
}

/** tries to upgrade a nonlinear constraint into a quadratic constraint */
static
SCIP_DECL_NONLINCONSUPGD(nonlinconsUpgdQuadratic)
{
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   int i;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   node = SCIPgetExprgraphNodeNonlinear(scip, cons);

   /* no interest in linear constraints */
   if( node == NULL )
      return SCIP_OKAY;

   /* if a quadratic expression has been simplified, then all children of the node should be variables */
   if( !SCIPexprgraphAreAllNodeChildrenVars(node) )
      return SCIP_OKAY;

   switch( SCIPexprgraphGetNodeOperator(node) )
   {
   case SCIP_EXPR_VARIDX:
   case SCIP_EXPR_CONST:
   case SCIP_EXPR_PLUS:
   case SCIP_EXPR_MINUS:
   case SCIP_EXPR_SUM:
   case SCIP_EXPR_LINEAR:
      /* these should not appear as exprgraphnodes after constraint presolving */
      return SCIP_OKAY;

   case SCIP_EXPR_DIV:
   case SCIP_EXPR_SQRT:
   case SCIP_EXPR_REALPOWER:
   case SCIP_EXPR_INTPOWER:
   case SCIP_EXPR_SIGNPOWER:
   case SCIP_EXPR_EXP:
   case SCIP_EXPR_LOG:
   case SCIP_EXPR_SIN:
   case SCIP_EXPR_COS:
   case SCIP_EXPR_TAN:
      /* case SCIP_EXPR_ERF: */
      /* case SCIP_EXPR_ERFI: */
   case SCIP_EXPR_MIN:
   case SCIP_EXPR_MAX:
   case SCIP_EXPR_ABS:
   case SCIP_EXPR_SIGN:
   case SCIP_EXPR_PRODUCT:
   case SCIP_EXPR_POLYNOMIAL:
   case SCIP_EXPR_USER:
      /* these do not look like an quadratic expression (assuming the expression graph simplifier did run) */
      return SCIP_OKAY;

   case SCIP_EXPR_MUL:
   case SCIP_EXPR_SQUARE:
   case SCIP_EXPR_QUADRATIC:
      /* these mean that we have something quadratic */
      break;

   case SCIP_EXPR_PARAM:
   case SCIP_EXPR_LAST:
   default:
      SCIPwarningMessage(scip, "unexpected expression operator %d in nonlinear constraint <%s>\n", SCIPexprgraphGetNodeOperator(node), SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   /* setup a quadratic constraint */

   if( upgdconsssize < 1 )
   {
      /* request larger upgdconss array */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsQuadratic(scip, &upgdconss[0], SCIPconsGetName(cons),
         SCIPgetNLinearVarsNonlinear(scip, cons), SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
         0, NULL, 0, NULL,
         SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
   assert(!SCIPconsIsStickingAtNode(cons));

   exprgraph = SCIPgetExprgraphNonlinear(scip, SCIPconsGetHdlr(cons));

   /* add variables from expression tree as "quadratic" variables to quadratic constraint */
   for( i = 0; i < SCIPexprgraphGetNodeNChildren(node); ++i )
   {
      assert(SCIPexprgraphGetNodeChildren(node)[i] != NULL);
      SCIP_CALL( SCIPaddQuadVarQuadratic(scip, upgdconss[0], (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[i]), 0.0, 0.0) );
   }

   switch( SCIPexprgraphGetNodeOperator(node) )
   {
   case SCIP_EXPR_MUL:
      /* expression is product of two variables, so add bilinear term to constraint */
      assert(SCIPexprgraphGetNodeNChildren(node) == 2);

      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, upgdconss[0],
            (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[0]),
            (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[1]),
            1.0) );

      break;

   case SCIP_EXPR_SQUARE:
      /* expression is square of a variable, so change square coefficient of quadratic variable */
      assert(SCIPexprgraphGetNodeNChildren(node) == 1);

      SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, upgdconss[0],
            (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[0]),
            1.0) );

      break;

   case SCIP_EXPR_QUADRATIC:
   {
      /* expression is quadratic */
      SCIP_QUADELEM* quadelems;
      int nquadelems;
      SCIP_Real* lincoefs;

      quadelems  = SCIPexprgraphGetNodeQuadraticQuadElements(node);
      nquadelems = SCIPexprgraphGetNodeQuadraticNQuadElements(node);
      lincoefs   = SCIPexprgraphGetNodeQuadraticLinearCoefs(node);

      SCIPaddConstantQuadratic(scip, upgdconss[0], SCIPexprgraphGetNodeQuadraticConstant(node));

      if( lincoefs != NULL )
         for( i = 0; i < SCIPexprgraphGetNodeNChildren(node); ++i )
            if( lincoefs[i] != 0.0 )
            {
               /* linear term */
               SCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, upgdconss[0],
                     (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[i]),
                     lincoefs[i]) );
            }

      for( i = 0; i < nquadelems; ++i )
      {
         assert(quadelems[i].idx1 < SCIPexprgraphGetNodeNChildren(node));
         assert(quadelems[i].idx2 < SCIPexprgraphGetNodeNChildren(node));

         if( quadelems[i].idx1 == quadelems[i].idx2 )
         {
            /* square term */
            SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, upgdconss[0],
                  (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[quadelems[i].idx1]),
                  quadelems[i].coef) );
         }
         else
         {
            /* bilinear term */
            SCIP_CALL( SCIPaddBilinTermQuadratic(scip, upgdconss[0],
                  (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[quadelems[i].idx1]),
                  (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(node)[quadelems[i].idx2]),
                  quadelems[i].coef) );
         }
      }

      break;
   }

   default:
      SCIPerrorMessage("you should not be here\n");
      return SCIP_ERROR;
   }  /*lint !e788 */

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyQuadratic)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrQuadratic(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nquadconsupgrades; ++i )
   {
      assert(conshdlrdata->quadconsupgrades[i] != NULL);
      SCIPfreeBlockMemory(scip, &conshdlrdata->quadconsupgrades[i]); /*lint !e866*/
   }
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->quadconsupgrades, conshdlrdata->quadconsupgradessize);
   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitQuadratic)
{  /*lint --e{715} */
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitQuadratic)
{  /*lint --e{715} */
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   return SCIP_OKAY;
}

/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return SCIP_OKAY;
}
#endif

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA*     consdata;
   int                c;
#ifndef NDEBUG
   int                i;
#endif

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->isremovedfixings )
      {
         SCIP_CALL( removeFixedVariables(scip, conss[c]) );
      }

      /* make sure we do not have duplicate bilinear terms, quad var terms, or linear vars */
      SCIP_CALL( mergeAndCleanBilinearTerms(scip, conss[c]) );
      SCIP_CALL( mergeAndCleanQuadVarTerms(scip, conss[c]) );
      SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );

      assert(consdata->isremovedfixings);
      assert(consdata->linvarsmerged);
      assert(consdata->quadvarsmerged);
      assert(consdata->bilinmerged);

#ifndef NDEBUG
      for( i = 0; i < consdata->nlinvars; ++i )
         assert(SCIPvarIsActive(consdata->linvars[i]));

      for( i = 0; i < consdata->nquadvars; ++i )
         assert(SCIPvarIsActive(consdata->quadvarterms[i].var));
#endif

      /* tell SCIP that we have something nonlinear */
      if( SCIPconsIsAdded(conss[c]) && consdata->nquadvars > 0 )
         SCIPenableNLP(scip);
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin)
 *
 *  @note Also called from consEnableQuadratic during solving stage.
 */
static
SCIP_DECL_CONSINITSOL(consInitsolQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check for a linear variable that can be increase or decreased without harming feasibility */
      consdataFindUnlockedLinearVar(scip, consdata);

      /* setup lincoefsmin, lincoefsmax */
      consdata->lincoefsmin = SCIPinfinity(scip);
      consdata->lincoefsmax = 0.0;
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         consdata->lincoefsmin = MIN(consdata->lincoefsmin, REALABS(consdata->lincoefs[i]));  /*lint !e666 */
         consdata->lincoefsmax = MAX(consdata->lincoefsmax, REALABS(consdata->lincoefs[i]));  /*lint !e666 */
      }

      /* add nlrow representation to NLP, if NLP had been constructed */
      if( SCIPisNLPConstructed(scip) && SCIPconsIsEnabled(conss[c]) )
      {
         if( consdata->nlrow == NULL )
         {
            /* compute curvature for the quadratic constraint if not done yet */
            SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );

            SCIP_CALL( createNlRow(scip, conss[c]) );
            assert(consdata->nlrow != NULL);
         }
         SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
      }

      /* setup sepaquadvars and sepabilinvar2pos */
      assert(consdata->sepaquadvars == NULL);
      assert(consdata->sepabilinvar2pos == NULL);
      if( consdata->nquadvars > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->sepaquadvars,     consdata->nquadvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->sepabilinvar2pos, consdata->nbilinterms) );

         /* make sure, quadratic variable terms are sorted */
         SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );

         for( i = 0; i < consdata->nquadvars; ++i )
            consdata->sepaquadvars[i] = consdata->quadvarterms[i].var;

         for( i = 0; i < consdata->nbilinterms; ++i )
         {
            SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, consdata->bilinterms[i].var2, &consdata->sepabilinvar2pos[i]) );
         }
      }

      if( conshdlrdata->checkfactorable )
      {
         /* check if constraint function is factorable, i.e., can be written as product of two linear functions */
         SCIP_CALL( checkFactorable(scip, conss[c]) );
      }

      /* compute gauge function using interior points per constraint, only when there are quadratic variables */
      if( conshdlrdata->gaugecuts && SCIPgetSubscipDepth(scip) == 0 && consdata->nquadvars > 0 )
      {
         SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */
         if( (consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs)) ||
               (consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)) )
         {
            SCIP_CALL( computeGauge(scip, conshdlr, conss[c]) );
         }
      }

      /* compute eigendecomposition for convex quadratics */
      if( conshdlrdata->projectedcuts && SCIPgetSubscipDepth(scip) == 0 && consdata->nquadvars > 0 )
      {
         SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */
         if( (consdata->isconvex && !SCIPisInfinity(scip, consdata->rhs)) ||
               (consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)) )
         {
            SCIP_CALL( computeED(scip, conshdlr, conss[c]) );
         }
      }

      /* mark constraint for propagation */
      SCIP_CALL( SCIPmarkConsPropagate(scip, conss[c]) );
      consdata->ispropagated = FALSE;
   }

   if( SCIPgetStage(scip) != SCIP_STAGE_INITSOLVE )
   {
      /* if called from consEnableQuadratic, then don't do below */
      return SCIP_OKAY;
   }

   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 && conshdlrdata->linearizeheursol )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      /* @todo Should we catch every new solution or only new *best* solutions */
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   if( nconss != 0 && !SCIPisIpoptAvailableIpopt() && !SCIPisInRestart(scip) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Quadratic constraint handler does not have LAPACK for eigenvalue computation. Will assume that matrices (with size > 2x2) are indefinite.\n");
   }

   /* reset flags and counters */
   conshdlrdata->sepanlp = FALSE;
   conshdlrdata->lastenfonode = NULL;
   conshdlrdata->nenforounds = 0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed)
 *
 *  @note Also called from consDisableQuadratic during solving stage.
 */
static
SCIP_DECL_CONSEXITSOL(consExitsolQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }

      assert(!SCIPconsIsEnabled(conss[c]) || consdata->sepaquadvars     != NULL || consdata->nquadvars == 0);   /*lint !e613 */
      assert(!SCIPconsIsEnabled(conss[c]) || consdata->sepabilinvar2pos != NULL || consdata->nquadvars == 0);   /*lint !e613 */
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->sepaquadvars,     consdata->nquadvars);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->sepabilinvar2pos, consdata->nbilinterms);

      SCIPfreeBlockMemoryArrayNull(scip, &consdata->factorleft,  consdata->nquadvars + 1);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->factorright, consdata->nquadvars + 1);

      SCIPfreeBlockMemoryArrayNull(scip, &consdata->interiorpoint, consdata->nquadvars);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->gaugecoefs, consdata->nquadvars);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->eigenvalues, consdata->nquadvars);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->eigenvectors, (int)(consdata->nquadvars*consdata->nquadvars));
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->bp, consdata->nquadvars);
   }

   if( SCIPgetStage(scip) != SCIP_STAGE_EXITSOLVE )
   {
      /* if called from consDisableQuadratic, then don't do below */
      return SCIP_OKAY;
   }

   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   /* free all stored bilinear terms in the constraint handler and constraint data; note that we might not want to
    * recollect all bilinear terms and therefore keep them even if consDisableQuadratic is called
    */
   SCIP_CALL( freeAllBilinearTerms(scip, conshdlrdata, conss, nconss) );

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteQuadratic)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(SCIPconsGetData(cons) == *consdata);

   SCIP_CALL( consdataFree(scip, consdata) );

   assert(*consdata == NULL);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransQuadratic)
{  
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int            i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->lhs, sourcedata->rhs,
         sourcedata->nlinvars, sourcedata->linvars, sourcedata->lincoefs,
         sourcedata->nquadvars, sourcedata->quadvarterms,
         sourcedata->nbilinterms, sourcedata->bilinterms,
         FALSE) );

   for( i = 0; i < targetdata->nlinvars; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->linvars[i], &targetdata->linvars[i]) );
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->linvars[i]) );
   }

   for( i = 0; i < targetdata->nquadvars; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->quadvarterms[i].var, &targetdata->quadvarterms[i].var) );
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->quadvarterms[i].var) );
   }

   for( i = 0; i < targetdata->nbilinterms; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->bilinterms[i].var1, &targetdata->bilinterms[i].var1) );
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->bilinterms[i].var2, &targetdata->bilinterms[i].var2) );

      if( SCIPvarCompare(targetdata->bilinterms[i].var1, targetdata->bilinterms[i].var2) > 0 )
      {
         SCIP_VAR* tmp;
         tmp = targetdata->bilinterms[i].var2;
         targetdata->bilinterms[i].var2 = targetdata->bilinterms[i].var1;
         targetdata->bilinterms[i].var1 = tmp;
      }
   }

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   SCIPdebugMsg(scip, "created transformed quadratic constraint ");
   SCIPdebugPrintCons(scip, *targetcons, NULL);

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR*          var;
   SCIP_ROW*          row;
   SCIP_Real*         x;
   int                c;
   int                i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *infeasible = FALSE;

   for( c = 0; c < nconss && !(*infeasible); ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613 */

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
      assert(consdata != NULL);

      row = NULL;

      if( consdata->nquadvars == 0 )
      {
         /* if we are actually linear, add the constraint as row to the LP */
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(conss[c]), SCIPconsGetName(conss[c]), consdata->lhs, consdata->rhs,
               SCIPconsIsLocal(conss[c]), FALSE , TRUE) );  /*lint !e613 */
         SCIP_CALL( SCIPaddVarsToRow(scip, row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
         SCIP_CALL( SCIPreleaseRow (scip, &row) );
         continue;
      }

      /* alloc memory for reference point */
      SCIP_CALL( SCIPallocBufferArray(scip, &x, consdata->nquadvars) );

      /* for convex parts, add linearizations in 5 points */
      if( (consdata->isconvex && !SCIPisInfinity(scip,  consdata->rhs)) ||
         (consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs)) )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real lambda;
         int k;

         for( k = 0; k < 5; ++k )
         {
            lambda = 0.1 * (k+1); /* lambda = 0.1, 0.2, 0.3, 0.4, 0.5 */
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               var = consdata->quadvarterms[i].var;
               lb = SCIPvarGetLbGlobal(var);
               ub = SCIPvarGetUbGlobal(var);

               if( ub > -INITLPMAXVARVAL )
                  lb = MAX(lb, -INITLPMAXVARVAL);
               if( lb <  INITLPMAXVARVAL )
                  ub = MIN(ub,  INITLPMAXVARVAL);

               /* make bounds finite */
               if( SCIPisInfinity(scip, -lb) )
                  lb = MIN(-10.0, ub - 0.1*REALABS(ub));  /*lint !e666 */
               if( SCIPisInfinity(scip,  ub) )
                  ub = MAX( 10.0, lb + 0.1*REALABS(lb));  /*lint !e666 */

               if( SCIPvarGetBestBoundType(var) == SCIP_BOUNDTYPE_LOWER )
                  x[i] = lambda * ub + (1.0 - lambda) * lb;
               else
                  x[i] = lambda * lb + (1.0 - lambda) * ub;
            }

            SCIP_CALL( generateCut(scip, conshdlr, conss[c], x, NULL, consdata->isconvex ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT, &row, NULL,
                  FALSE, -SCIPinfinity(scip)) );  /*lint !e613 */
            if( row != NULL )
            {
               SCIPdebugMsg(scip, "initlp adds row <%s> for lambda = %g of conss <%s>\n", SCIProwGetName(row), lambda, SCIPconsGetName(conss[c]));  /*lint !e613 */
               SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow (scip, &row) );
            }
         }
      }

      /* for concave parts, add underestimator w.r.t. at most 2 reference points */
      if( !(*infeasible) && ((! consdata->isconvex && !SCIPisInfinity(scip,  consdata->rhs))
            || (! consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs))) )
      {
         SCIP_Bool unbounded;
         SCIP_Bool possquare;
         SCIP_Bool negsquare;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real lambda;
         int k;

         unbounded = FALSE; /* whether there are unbounded variables */
         possquare = FALSE; /* whether there is a positive square term */
         negsquare = FALSE; /* whether there is a negative square term */
         for( k = 0; k < 2; ++k )
         {
            /* Set reference point to 0 projected on bounds for unbounded variables or in between lower and upper bound
             * for bounded variables in the first round, we set it closer to the best bound for one part of the
             * variables, in the second closer to the best bound for the other part of the variables.
             * Additionally, we use slightly different weights for each variable.
             * The reason for the latter is, that for a bilinear term with bounded variables, there are always two linear underestimators
             * if the same weight is used for both variables of a product, then rounding and luck decides which underestimator is chosen
             * of course, the possible number of cuts is something in the order of 2^nquadvars, and we choose two of them here.
             */
            for( i = 0; i < consdata->nquadvars; ++i )
            {
               var = consdata->quadvarterms[i].var;
               lb = SCIPvarGetLbGlobal(var);
               ub = SCIPvarGetUbGlobal(var);

               if( SCIPisInfinity(scip, -lb) )
               {
                  if( SCIPisInfinity(scip, ub) )
                     x[i] = 0.0;
                  else
                     x[i] = MIN(0.0, ub);
                  unbounded = TRUE;
               }
               else
               {
                  if( SCIPisInfinity(scip, ub) )
                  {
                     x[i] = MAX(0.0, lb);
                     unbounded = TRUE;
                  }
                  else
                  {
                     lambda = 0.4 + 0.2 * ((i+k)%2) + 0.01 * i / (double)consdata->nquadvars;
                     x[i] = lambda * SCIPvarGetBestBoundLocal(var) + (1.0-lambda) * SCIPvarGetWorstBoundLocal(var);
                  }
               }

               possquare |= consdata->quadvarterms[i].sqrcoef > 0.0;  /*lint !e514 */
               negsquare |= consdata->quadvarterms[i].sqrcoef < 0.0;  /*lint !e514 */
            }

            if( !consdata->isconvex  && !SCIPisInfinity(scip,  consdata->rhs) )
            {
               SCIP_CALL( generateCut(scip, conshdlr, conss[c], x, NULL, SCIP_SIDETYPE_RIGHT, &row, NULL,
                     conshdlrdata->checkcurvature, -SCIPinfinity(scip)) );  /*lint !e613 */
               if( row != NULL )
               {
                  SCIPdebugMsg(scip, "initlp adds row <%s> for rhs of conss <%s>, round %d\n", SCIProwGetName(row), SCIPconsGetName(conss[c]), k);  /*lint !e613 */
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
                  SCIP_CALL( SCIPreleaseRow (scip, &row) );
               }
            }
            if( !(*infeasible) && !consdata->isconcave && !SCIPisInfinity(scip, -consdata->lhs) )
            {
               SCIP_CALL( generateCut(scip, conshdlr, conss[c], x, NULL, SCIP_SIDETYPE_LEFT, &row, NULL,
                     conshdlrdata->checkcurvature, -SCIPinfinity(scip)) );  /*lint !e613 */
               if( row != NULL )
               {
                  SCIPdebugMsg(scip, "initlp adds row <%s> for lhs of conss <%s>, round %d\n", SCIProwGetName(row), SCIPconsGetName(conss[c]), k);  /*lint !e613 */
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
                  SCIP_CALL( SCIPreleaseRow (scip, &row) );
               }
            }

            /* if there are unbounded variables, then there is typically only at most one possible underestimator, so don't try another round
             * similar, if there are no bilinear terms and no linearizations of square terms, then the reference point does not matter, so don't do another round */
            if( unbounded ||
               (consdata->nbilinterms == 0 && (!possquare || SCIPisInfinity(scip,  consdata->rhs))) ||
               (consdata->nbilinterms == 0 && (!negsquare || SCIPisInfinity(scip, -consdata->lhs))) )
               break;
         }
      }

      SCIPfreeBufferArray(scip, &x);
   }

   /* store all bilinear terms into constraint handler data; this code is not in initsolve because the sub-NLP
    * heuristic triggers this callback and should not collect all bilinear terms
    */
   SCIP_CALL( storeAllBilinearTerms(scip, conshdlrdata, conss, nconss) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpQuadratic)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool          solviolbounds;
   SCIP_CONS*         maxviolcon;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, &solviolbounds, &maxviolcon) );

   /* don't try to separate solutions that violate variable bounds */
   if( solviolbounds )
      return SCIP_OKAY;

   /* if nothing violated, then nothing to separate */
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* at root, check if we want to solve the NLP relaxation and use its solutions as reference point
    * if there is something convex, then linearizing in the solution of the NLP relaxation can be very useful
    */
   if( SCIPgetDepth(scip) == 0 && !conshdlrdata->sepanlp &&
      (SCIPgetNContVars(scip) >= conshdlrdata->sepanlpmincont * SCIPgetNVars(scip) ||
         (SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY && conshdlrdata->sepanlpmincont <= 1.0)) &&
      SCIPisNLPConstructed(scip) && SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_CONSDATA*  consdata;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Bool       solvednlp;
      int c;

      solstat = SCIPgetNLPSolstat(scip);
      solvednlp = FALSE;
      if( solstat == SCIP_NLPSOLSTAT_UNKNOWN )
      {
         /* NLP is not solved yet, so we might want to do this
          * but first check whether there is a violated constraint side which corresponds to a convex function
          */
         for( c = 0; c < nconss; ++c )
         {
            assert(conss[c] != NULL);  /*lint !e613 */

            consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
            assert(consdata != NULL);

            /* skip feasible constraints */
            if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
               continue;

            /* make sure curvature has been checked */
            SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkcurvature) );  /*lint !e613 */

            if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && consdata->isconvex) ||
               ( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && consdata->isconcave) )
               break;
         }

         if( c < nconss )
         {
            /* try to solve NLP and update solstat */

            /* ensure linear conss are in NLP */
            if( conshdlrdata->subnlpheur != NULL )
            {
               SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, conshdlrdata->subnlpheur, TRUE, TRUE) );
            }

            /* set LP solution as starting values, if available */
            if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );
            }

            /* SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, 1) ); */
            SCIP_CALL( SCIPsolveNLP(scip) );

            solstat = SCIPgetNLPSolstat(scip);
            SCIPdebugMsg(scip, "solved NLP relax, solution status: %d\n", solstat);

            solvednlp = TRUE;
         }
      }

      conshdlrdata->sepanlp = TRUE;

      if( solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      {
         SCIPdebugMsg(scip, "NLP relaxation is globally infeasible, thus can cutoff node\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         /* if we have feasible NLP solution, generate linearization cuts there */
         SCIP_Bool lpsolseparated;
         SCIP_SOL* nlpsol;

         SCIP_CALL( SCIPcreateNLPSol(scip, &nlpsol, NULL) );
         assert(nlpsol != NULL);

         /* if we solved the NLP and solution is integral, then pass it to trysol heuristic */
         if( solvednlp && conshdlrdata->trysolheur != NULL )
         {
            int nfracvars;

            nfracvars = 0;
            if( SCIPgetNBinVars(scip) > 0 || SCIPgetNIntVars(scip) > 0 )
            {
               SCIP_CALL( SCIPgetNLPFracVars(scip, NULL, NULL, NULL, &nfracvars, NULL) );
            }

            if( nfracvars == 0 )
            {
               SCIPdebugMsg(scip, "pass solution with obj. value %g to trysol\n", SCIPgetSolOrigObj(scip, nlpsol));
               SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, nlpsol) );
            }
         }

         SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, nlpsol, &lpsolseparated, SCIPgetSepaMinEfficacy(scip)) );

         SCIP_CALL( SCIPfreeSol(scip, &nlpsol) );

         /* if a cut that separated the LP solution was added, then return, otherwise continue with usual separation in LP solution */
         if( lpsolseparated )
         {
            SCIPdebugMsg(scip, "linearization cuts separate LP solution\n");
            *result = SCIP_SEPARATED;

            return SCIP_OKAY;
         }
      }
   }
   /* if we do not want to try solving the NLP, or have no NLP, or have no NLP solver, or solving the NLP failed,
    * or separating with NLP solution as reference point failed, then try (again) with LP solution as reference point
    */

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, SCIPgetSepaMinEfficacy(scip), FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolQuadratic)
{
   SCIP_Bool          solviolbounds;
   SCIP_CONS*         maxviolcon;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(sol != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conss, nconss, sol, &solviolbounds, &maxviolcon) );

   /* don't separate solution that are outside variable bounds */
   if( solviolbounds )
      return SCIP_OKAY;

   /* if nothing violated, then nothing to separate */
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, SCIPgetSepaMinEfficacy(scip), FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpQuadratic)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, solinfeasible, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxQuadratic)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, solinfeasible, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsQuadratic)
{  /*lint --e{715}*/
   SCIP_Bool          solviolbounds;
   SCIP_CONS*         maxviolcon;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        propresult;
   SCIP_VAR*          var;
   int                c;
   int                i;
   int                nchgbds;
   int                nnotify;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, &solviolbounds, &maxviolcon) );

   /* pseudo solutions should be within bounds by definition */
   assert(!solviolbounds);

   if( maxviolcon == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   SCIPdebugMsg(scip, "enfops with max violation in cons <%s>\n", SCIPconsGetName(maxviolcon));

   /* run domain propagation */
   nchgbds = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &nchgbds) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we are not feasible and we cannot proof that the whole node is infeasible
    * -> collect all variables in violated constraints for branching
    */
   nnotify = 0;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         var = consdata->linvars[i];
         if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++nnotify;
         }
      }

      for( i = 0; i < consdata->nquadvars; ++i )
      {
         var = consdata->quadvarterms[i].var;
         if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++nnotify;
         }
      }
   }

   if( nnotify == 0 )
   {
      SCIP_Bool addedcons;
      SCIP_Bool reduceddom;
      SCIP_Bool infeasible;

      /* if no branching candidate found, then all variables are almost fixed
       * calling replaceByLinearConstraints() should lead to fix all almost-fixed quadratic variables, and possibly replace some quad. conss by linear ones
       */
      SCIP_CALL( replaceByLinearConstraints(scip, conss, nconss, &addedcons, &reduceddom, &infeasible) );
      if( addedcons )
      {
         *result = SCIP_CONSADDED;
         return SCIP_OKAY;
      }
      if( reduceddom )
      {
         *result = SCIP_REDUCEDDOM;
         return SCIP_OKAY;
      }
      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      SCIPdebugMsg(scip, "All variables in violated constraints fixed (up to epsilon). Cannot find branching candidate. Forcing solution of LP.\n");
      *result = SCIP_SOLVELP;
   }

   assert(*result == SCIP_SOLVELP || (*result == SCIP_INFEASIBLE && nnotify > 0));
   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropQuadratic)
{
   int         nchgbds;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   nchgbds = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nmarkedconss, result, &nchgbds) );

   return SCIP_OKAY;
}  /*lint !e715 */

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolQuadratic)
{  /*lint --e{715,788}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        solveresult;
   SCIP_Bool          redundant;
   SCIP_Bool          havechange;
   SCIP_Bool          doreformulations;
   int                c;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   /* if other presolvers did not find enough changes for another presolving round and we are in exhaustive presolving,
    * then try the reformulations (replacing products with binaries, disaggregation, setting default variable bounds)
    * otherwise, we wait with these
    */
   doreformulations = ((presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0) && SCIPisPresolveFinished(scip);
   SCIPdebugMsg(scip, "presolving will %swait with reformulation\n", doreformulations ? "not " : "");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      SCIPdebugMsg(scip, "process constraint <%s>\n", SCIPconsGetName(conss[c]));
      SCIPdebugPrintCons(scip, conss[c], NULL);

      if( !consdata->initialmerge )
      {
         SCIP_CALL( mergeAndCleanBilinearTerms(scip, conss[c]) );
         SCIP_CALL( mergeAndCleanQuadVarTerms(scip, conss[c]) );
         SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );
         consdata->initialmerge = TRUE;
      }

      havechange = FALSE;
#ifdef CHECKIMPLINBILINEAR
      if( consdata->isimpladded && (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 )
      {
         int nbilinremoved;
         SCIP_CALL( presolveApplyImplications(scip, conss[c], &nbilinremoved) );
         if( nbilinremoved > 0 )
         {
            *nchgcoefs += nbilinremoved;
            havechange = TRUE;
            *result = SCIP_SUCCESS;
         }
         assert(!consdata->isimpladded);
      }
#endif
      /* call upgrade methods if the constraint has not been presolved yet or there has been a bound tightening or possibly be a change in variable type
       * we want to do this before (multi)aggregated variables are replaced, since that may change structure, e.g., introduce bilinear terms
       */
      if( !consdata->ispresolved || !consdata->ispropagated || nnewchgvartypes > 0 )
      {
         SCIP_Bool upgraded;

         SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss, presoltiming) );
         if( upgraded )
         {
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      if( !consdata->isremovedfixings )
      {
         SCIP_CALL( removeFixedVariables(scip, conss[c]) );
         assert(consdata->isremovedfixings);
         havechange = TRUE;
      }

      /* try to "solve" the constraint, e.g., reduce to a variable aggregation */
      SCIP_CALL( presolveSolve(scip, conss[c], &solveresult, &redundant, naggrvars) );
      if( solveresult == SCIP_CUTOFF )
      {
         SCIPdebugMsg(scip, "solving constraint <%s> says problem is infeasible in presolve\n", SCIPconsGetName(conss[c]));
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( redundant )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         ++*ndelconss;
         *result = SCIP_SUCCESS;
         break;
      }
      if( solveresult == SCIP_SUCCESS )
      {
         *result = SCIP_SUCCESS;
         havechange = TRUE;
      }

      /* @todo divide constraint by gcd of coefficients if all are integral */

      if( doreformulations )
      {
         int naddconss_old;

         naddconss_old = *naddconss;

         SCIP_CALL( presolveTryAddAND(scip, conshdlr, conss[c], naddconss) );
         assert(*naddconss >= naddconss_old);

         if( *naddconss == naddconss_old )
         {
            /* user not so empathic about AND, or we don't have products of two binaries, so try this more general reformulation */
            SCIP_CALL( presolveTryAddLinearReform(scip, conshdlr, conss[c], naddconss) );
            assert(*naddconss >= naddconss_old);
         }

         if( conshdlrdata->maxdisaggrsize > 1 )
         {
            /* try disaggregation, if enabled */
            SCIP_CALL( presolveDisaggregate(scip, conshdlr, conss[c], naddconss) );
         }

         if( *naddconss > naddconss_old )
         {
            /* if something happened, report success and cleanup constraint */
            *result = SCIP_SUCCESS;
            havechange = TRUE;
            SCIP_CALL( mergeAndCleanBilinearTerms(scip, conss[c]) );
            SCIP_CALL( mergeAndCleanQuadVarTerms(scip, conss[c]) );
            SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );
         }
      }

      if( consdata->nlinvars == 0 && consdata->nquadvars == 0 )
      {
         /* all variables fixed or removed, constraint function is 0.0 now */
         if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasPositive(scip, consdata->lhs)) ||
            ( !SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasNegative(scip, consdata->rhs)) )
         { /* left hand side positive or right hand side negative */
            SCIPdebugMsg(scip, "constraint <%s> is constant and infeasible\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            ++*ndelconss;
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         /* left and right hand side are consistent */
         SCIPdebugMsg(scip, "constraint <%s> is constant and feasible, deleting\n", SCIPconsGetName(conss[c]));
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         ++*ndelconss;
         *result = SCIP_SUCCESS;
         continue;
      }

      if( (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 && !consdata->ispropagated )
      {
         /* try domain propagation if there were bound changes or constraint has changed (in which case, processVarEvents may have set ispropagated to false) */
         SCIP_RESULT propresult;
         int roundnr;

         roundnr = 0;
         do
         {
            ++roundnr;

            SCIPdebugMsg(scip, "starting domain propagation round %d of %d\n", roundnr, conshdlrdata->maxproproundspresolve);

            if( !consdata->ispropagated )
            {
               consdata->ispropagated = TRUE;

               SCIP_CALL( propagateBoundsCons(scip, conshdlr, conss[c], &propresult, nchgbds, &redundant) );

               if( propresult == SCIP_CUTOFF )
               {
                  SCIPdebugMsg(scip, "propagation on constraint <%s> says problem is infeasible in presolve\n",
                     SCIPconsGetName(conss[c]));
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }

               /* delete constraint if found redundant by bound tightening */
               if( redundant )
               {
                  SCIP_CALL( SCIPdelCons(scip, conss[c]) );
                  ++*ndelconss;
                  *result = SCIP_SUCCESS;
                  break;
               }

               if( propresult == SCIP_REDUCEDDOM )
               {
                  *result = SCIP_SUCCESS;
                  havechange = TRUE;
               }
            }
         }
         while( !consdata->ispropagated && roundnr < conshdlrdata->maxproproundspresolve );

         if( redundant )
            continue;
      }

      /* check if we have a single linear continuous variable that we can make implicit integer */
      if( (nnewchgvartypes != 0 || havechange || !consdata->ispresolved)
         && (SCIPisEQ(scip, consdata->lhs, consdata->rhs) && SCIPisIntegral(scip, consdata->lhs)) )
      {
         int       ncontvar;
         SCIP_VAR* candidate;
         SCIP_Bool fail;

         fail = FALSE;
         candidate = NULL;
         ncontvar = 0;

         for( i = 0; !fail && i < consdata->nlinvars; ++i )
         {
            if( !SCIPisIntegral(scip, consdata->lincoefs[i]) )
            {
               fail = TRUE;
            }
            else if( SCIPvarGetType(consdata->linvars[i]) == SCIP_VARTYPE_CONTINUOUS )
            {
               if( ncontvar > 0 ) /* now at 2nd continuous variable */
                  fail = TRUE;
               else if( SCIPisEQ(scip, ABS(consdata->lincoefs[i]), 1.0) )
                  candidate = consdata->linvars[i];
               ++ncontvar;
            }
         }
         for( i = 0; !fail && i < consdata->nquadvars; ++i )
            fail = SCIPvarGetType(consdata->quadvarterms[i].var) == SCIP_VARTYPE_CONTINUOUS ||
               !SCIPisIntegral(scip, consdata->quadvarterms[i].lincoef) ||
               !SCIPisIntegral(scip, consdata->quadvarterms[i].sqrcoef);
         for( i = 0; !fail && i < consdata->nbilinterms; ++i )
            fail = !SCIPisIntegral(scip, consdata->bilinterms[i].coef);

         if( !fail && candidate != NULL )
         {
            SCIP_Bool infeasible;

            SCIPdebugMsg(scip, "make variable <%s> implicit integer due to constraint <%s>\n", SCIPvarGetName(candidate), SCIPconsGetName(conss[c]));

            SCIP_CALL( SCIPchgVarType(scip, candidate, SCIP_VARTYPE_IMPLINT, &infeasible) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, "infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(candidate));
               *result = SCIP_CUTOFF;

               return SCIP_OKAY;
            }

            ++(*nchgvartypes);
            *result = SCIP_SUCCESS;
            havechange = TRUE;
         }
      }

      /* call upgrade methods again if constraint has been changed */
      if( havechange )
      {
         SCIP_Bool upgraded;

         SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss, presoltiming) );
         if( upgraded )
         {
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      /* fix quadratic variables with proper square coefficients contained in a single quadratic constraint to their
       * upper or lower bounds
       */
      if( (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 && conshdlrdata->checkquadvarlocks != 'd'
         && SCIPisPresolveFinished(scip) )
      {
         SCIP_CONS* cons;
         SCIP_VAR* vars[2];
         SCIP_BOUNDTYPE boundtypes[2];
         SCIP_Real bounds[2];
         char name[SCIP_MAXSTRLEN];

         /* merge variables in order to get correct locks for quadratic variables */
         if( !consdata->initialmerge )
         {
            SCIP_CALL( mergeAndCleanBilinearTerms(scip, conss[c]) );
            SCIP_CALL( mergeAndCleanQuadVarTerms(scip, conss[c]) );
            SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );
            consdata->initialmerge = TRUE;
         }

         for( i = 0; i < consdata->nquadvars; ++i )
         {
            if( hasQuadvarHpProperty(scip, consdata, i) )
            {
               SCIP_VAR* var;

               var = consdata->quadvarterms[i].var;
               assert(var != NULL);

               /* try to change the variable type to binary */
               if( conshdlrdata->checkquadvarlocks == 't' && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0) && SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0) )
               {
                  SCIP_Bool infeasible;

                  assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
                  SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );

                  if( infeasible )
                  {
                     SCIPdebugMsg(scip, "detect infeasibility after changing variable <%s> to binary type\n", SCIPvarGetName(var));
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
               }
               /* add bound disjunction constraint if bounds of variable are finite */
               else if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) && !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
               {
                  vars[0] = var;
                  vars[1] = var;
                  boundtypes[0] = SCIP_BOUNDTYPE_LOWER;
                  boundtypes[1] = SCIP_BOUNDTYPE_UPPER;
                  bounds[0] = SCIPvarGetUbGlobal(var);
                  bounds[1] = SCIPvarGetLbGlobal(var);

                  SCIPdebugMsg(scip, "add bound disjunction constraint for %s\n", SCIPvarGetName(var));

                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quadvarbnddisj_%s", SCIPvarGetName(var));
                  SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, 2, vars, boundtypes, bounds, TRUE, TRUE,
                        TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

                  SCIP_CALL( SCIPaddCons(scip, cons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }

               *result = SCIP_SUCCESS;
            }
         }
      }

      consdata->ispresolved = TRUE;
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool      haslb;
   SCIP_Bool      hasub;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslb = !SCIPisInfinity(scip, -consdata->lhs);
   hasub = !SCIPisInfinity(scip, consdata->rhs);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( consdata->lincoefs[i] > 0 )
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
      }
   }

   for( i = 0; i < consdata->nquadvars; ++i )
   {
      /* @todo try to be more clever, but variable locks that depend on the bounds of other variables are not trival to maintain */
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->quadvarterms[i].var, nlockspos+nlocksneg, nlockspos+nlocksneg) );
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableQuadratic)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));
   assert(SCIPconsIsActive(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMsg(scip, "enable cons <%s>\n", SCIPconsGetName(cons));

   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITPRESOLVE )
   {
      /* merge duplicate bilinear terms, move quad terms that are linear to linear vars */
      SCIP_CALL( mergeAndCleanBilinearTerms(scip, cons) );
      SCIP_CALL( mergeAndCleanQuadVarTerms(scip, cons) );
      SCIP_CALL( mergeAndCleanLinearVars(scip, cons) );
   }

   /* catch variable events */
   SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );

   /* initialize solving data */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( consInitsolQuadratic(scip, conshdlr, &cons, 1) );
   }

   return SCIP_OKAY;
}

/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMsg(scip, "disable cons <%s>\n", SCIPconsGetName(cons));

   /* free solving data */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( consExitsolQuadratic(scip, conshdlr, &cons, 1, FALSE) );
   }

   /* drop variable events */
   SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->nlinvars == 0 && consdata->nquadvars == 0 )
   {
      SCIPinfoMessage(scip, file, "0 ");
   }
   else
   {
      SCIP_VAR*** monomialvars;
      SCIP_Real** monomialexps;
      SCIP_Real*  monomialcoefs;
      int*        monomialnvars;
      int         nmonomials;
      int         monomialssize;
      int         j;

      monomialssize = consdata->nlinvars + 2 * consdata->nquadvars + consdata->nbilinterms;
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars,  monomialssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialexps,  monomialssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialcoefs, monomialssize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomialnvars, monomialssize) );

      nmonomials = 0;
      for( j = 0; j < consdata->nlinvars; ++j )
      {
         assert(nmonomials < monomialssize);

         SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 1) );  /*lint !e866 */

         monomialvars[nmonomials][0] = consdata->linvars[j];
         monomialexps[nmonomials] = NULL;
         monomialcoefs[nmonomials] = consdata->lincoefs[j];
         monomialnvars[nmonomials] = 1;
         ++nmonomials;
      }

      for( j = 0; j < consdata->nquadvars; ++j )
      {
         if( consdata->quadvarterms[j].lincoef != 0.0 )
         {
            assert(nmonomials < monomialssize);

            SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 1) );  /*lint !e866 */

            monomialvars[nmonomials][0] = consdata->quadvarterms[j].var;
            monomialexps[nmonomials] = NULL;
            monomialcoefs[nmonomials] = consdata->quadvarterms[j].lincoef;
            monomialnvars[nmonomials] = 1;
            ++nmonomials;
         }

         if( consdata->quadvarterms[j].sqrcoef != 0.0 )
         {
            assert(nmonomials < monomialssize);

            SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 1) );  /*lint !e866 */
            SCIP_CALL( SCIPallocBufferArray(scip, &monomialexps[nmonomials], 1) );  /*lint !e866 */

            monomialvars[nmonomials][0] = consdata->quadvarterms[j].var;
            monomialexps[nmonomials][0] = 2.0;
            monomialcoefs[nmonomials] = consdata->quadvarterms[j].sqrcoef;
            monomialnvars[nmonomials] = 1;
            ++nmonomials;
         }
      }

      for( j = 0; j < consdata->nbilinterms; ++j )
      {
         assert(nmonomials < monomialssize);

         SCIP_CALL( SCIPallocBufferArray(scip, &monomialvars[nmonomials], 2) );  /*lint !e866 */

         monomialvars[nmonomials][0] = consdata->bilinterms[j].var1;
         monomialvars[nmonomials][1] = consdata->bilinterms[j].var2;
         monomialexps[nmonomials] = NULL;
         monomialcoefs[nmonomials] = consdata->bilinterms[j].coef;
         monomialnvars[nmonomials] = 2;
         ++nmonomials;
      }

      SCIP_CALL( SCIPwriteVarsPolynomial(scip, file, monomialvars, monomialexps, monomialcoefs, monomialnvars, nmonomials, TRUE) );

      for( j = 0; j < nmonomials; ++j )
      {
         SCIPfreeBufferArray(scip, &monomialvars[j]);
         SCIPfreeBufferArrayNull(scip, &monomialexps[j]);
      }

      SCIPfreeBufferArray(scip, &monomialvars);
      SCIPfreeBufferArray(scip, &monomialexps);
      SCIPfreeBufferArray(scip, &monomialcoefs);
      SCIPfreeBufferArray(scip, &monomialnvars);
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   }
   else
   {
      /* should be ignored by parser */
      SCIPinfoMessage(scip, file, " [free]");
   }

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   int                c;
   SCIP_Bool          maypropfeasible; /* whether we may be able to propose a feasible solution */
   SCIP_Bool          solviolbounds;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   maxviol = 0.0;
   maypropfeasible = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL) &&
      SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED && SCIPgetStage(scip) <= SCIP_STAGE_SOLVING;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, &solviolbounds) );
      assert(!solviolbounds);  /* see also issue #627 */

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");
            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g (scaled: %.15g)\n", consdata->lhs - consdata->activity, consdata->lhsviol);
            }
            if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g (scaled: %.15g)\n", consdata->activity - consdata->rhs, consdata->rhsviol);
            }
         }
         if( (conshdlrdata->subnlpheur == NULL || sol == NULL) && !maypropfeasible && !completely )
            return SCIP_OKAY;
         if( consdata->lhsviol > maxviol || consdata->rhsviol > maxviol )
            maxviol = consdata->lhsviol + consdata->rhsviol;

         /* do not try to shift linear variables if activity is at infinity (leads to setting variable to infinity in solution, which is not allowed) */
         if( maypropfeasible && SCIPisInfinity(scip, REALABS(consdata->activity)) )
            maypropfeasible = FALSE;

         if( maypropfeasible )
         {
            /* update information on linear variables that may be in- or decreased, if initsolve has not done so yet */
            if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED && SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE )
               consdataFindUnlockedLinearVar(scip, consdata);

            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               /* check if there is a variable which may help to get the left hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)));
               /* check if there is a variable which may help to get the right hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0) )
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

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyQuadratic)
{
   SCIP_CONSDATA*    consdata;
   SCIP_CONSDATA*    targetconsdata;
   SCIP_VAR**        linvars;
   SCIP_QUADVARTERM* quadvarterms;
   SCIP_BILINTERM*   bilinterms;
   int               i;
   int               j;
   int               k;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);

   linvars = NULL;
   quadvarterms = NULL;
   bilinterms = NULL;

   *valid = TRUE;

   if( consdata->nlinvars != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &linvars, consdata->nlinvars) );
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->linvars[i], &linvars[i], varmap, consmap, global, valid) );
         assert(!(*valid) || linvars[i] != NULL);

         /* we do not copy, if a variable is missing */
         if( !(*valid) )
            goto TERMINATE;
      }
   }

   if( consdata->nbilinterms != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &bilinterms, consdata->nbilinterms) );
   }

   if( consdata->nquadvars != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &quadvarterms, consdata->nquadvars) );
      for( i = 0; i < consdata->nquadvars; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->quadvarterms[i].var, &quadvarterms[i].var, varmap, consmap, global, valid) );
         assert(!(*valid) || quadvarterms[i].var != NULL);

         /* we do not copy, if a variable is missing */
         if( !(*valid) )
            goto TERMINATE;

         quadvarterms[i].lincoef   = consdata->quadvarterms[i].lincoef;
         quadvarterms[i].sqrcoef   = consdata->quadvarterms[i].sqrcoef;
         quadvarterms[i].eventdata = NULL;
         quadvarterms[i].nadjbilin = consdata->quadvarterms[i].nadjbilin;
         quadvarterms[i].adjbilin  = consdata->quadvarterms[i].adjbilin;

         assert(consdata->nbilinterms != 0 || consdata->quadvarterms[i].nadjbilin == 0);

         for( j = 0; j < consdata->quadvarterms[i].nadjbilin; ++j )
         {
            assert(bilinterms != NULL);

            k = consdata->quadvarterms[i].adjbilin[j];
            assert(consdata->bilinterms[k].var1 != NULL);
            assert(consdata->bilinterms[k].var2 != NULL);
            if( consdata->bilinterms[k].var1 == consdata->quadvarterms[i].var )
            {
               assert(consdata->bilinterms[k].var2 != consdata->quadvarterms[i].var);
               bilinterms[k].var1 = quadvarterms[i].var;
            }
            else
            {
               assert(consdata->bilinterms[k].var2 == consdata->quadvarterms[i].var);
               bilinterms[k].var2 = quadvarterms[i].var;
            }
            bilinterms[k].coef = consdata->bilinterms[k].coef;
         }
      }
   }

   assert(stickingatnode == FALSE);
   SCIP_CALL( SCIPcreateConsQuadratic2(scip, cons, name ? name : SCIPconsGetName(sourcecons),
         consdata->nlinvars, linvars, consdata->lincoefs,
         consdata->nquadvars, quadvarterms,
         consdata->nbilinterms, bilinterms,
         consdata->lhs, consdata->rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

   /* copy information on curvature */
   targetconsdata = SCIPconsGetData(*cons);
   targetconsdata->isconvex      = consdata->isconvex;
   targetconsdata->isconcave     = consdata->isconcave;
   targetconsdata->iscurvchecked = consdata->iscurvchecked;

 TERMINATE:
   SCIPfreeBufferArrayNull(sourcescip, &quadvarterms);
   SCIPfreeBufferArrayNull(sourcescip, &bilinterms);
   SCIPfreeBufferArrayNull(sourcescip, &linvars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseQuadratic)
{  /*lint --e{715}*/
   SCIP_VAR*** monomialvars;
   SCIP_Real** monomialexps;
   SCIP_Real*  monomialcoefs;
   char*       endptr;
   int*        monomialnvars;
   int         nmonomials;

   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   (*success) = FALSE;

   /* return of string empty */
   if( !*str )
      return SCIP_OKAY;

   /* ignore whitespace */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, &endptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_OKAY;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the first coefficient, but not a left-hand-side */
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

   SCIP_CALL( SCIPparseVarsPolynomial(scip, str, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, &nmonomials, &endptr, success) );

   if( *success )
   {
      /* check for right hand side */
      str = endptr;

      /* ignore whitespace */
      while( isspace((unsigned char)*str) )
         ++str;

      if( *str && str[0] == '<' && str[1] == '=' )
      {
         /* we seem to get a right-hand-side */
         str += 2;

         if( !SCIPstrToRealValue(str, &rhs, &endptr) )
         {
            SCIPerrorMessage("error parsing right-hand-side from %s\n", str);
            *success = FALSE;
         }
      }
      else if( *str && str[0] == '>' && str[1] == '=' )
      {
         /* we seem to get a left-hand-side */
         str += 2;

         /* we should not have a left-hand-side already */
         assert(SCIPisInfinity(scip, -lhs));

         if( !SCIPstrToRealValue(str, &lhs, &endptr) )
         {
            SCIPerrorMessage("error parsing left-hand-side from %s\n", str);
            *success = FALSE;
         }
      }
      else if( *str && str[0] == '=' && str[1] == '=' )
      {
         /* we seem to get a left- and right-hand-side */
         str += 2;

         /* we should not have a left-hand-side already */
         assert(SCIPisInfinity(scip, -lhs));

         if( !SCIPstrToRealValue(str, &lhs, &endptr) )
         {
            SCIPerrorMessage("error parsing left-hand-side from %s\n", str);
            *success = FALSE;
         }
         else
         {
            rhs = lhs;
         }
      }
   }

   if( *success )
   {
      int i;

      /* setup constraint */
      assert(stickingatnode == FALSE);
      SCIP_CALL( SCIPcreateConsQuadratic(scip, cons, name, 0, NULL, NULL,
            0, NULL, NULL, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );

      for( i = 0; i < nmonomials; ++i )
      {
         if( monomialnvars[i] == 0 )
         {
            /* constant monomial */
            SCIPaddConstantQuadratic(scip, *cons, monomialcoefs[i]);
         }
         else if( monomialnvars[i] == 1 && monomialexps[i][0] == 1.0 )
         {
            /* linear monomial */
            SCIP_CALL( SCIPaddLinearVarQuadratic(scip, *cons, monomialvars[i][0], monomialcoefs[i]) );
         }
         else if( monomialnvars[i] == 1 && monomialexps[i][0] == 2.0 )
         {
            /* square monomial */
            SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, monomialvars[i][0], 0.0, monomialcoefs[i]) );
         }
         else if( monomialnvars[i] == 2 && monomialexps[i][0] == 1.0 && monomialexps[i][1] == 1.0 )
         {
            /* bilinear term */
            SCIP_VAR* var1;
            SCIP_VAR* var2;
            int pos;

            var1 = monomialvars[i][0];
            var2 = monomialvars[i][1];
            if( var1 == var2 )
            {
               /* actually a square term */
               SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, var1, 0.0, monomialcoefs[i]) );
            }
            else
            {
               SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, *cons, var1, &pos) );
               if( pos == -1 )
               {
                  SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, var1, 0.0, 0.0) );
               }

               SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, *cons, var2, &pos) );
               if( pos == -1 )
               {
                  SCIP_CALL( SCIPaddQuadVarQuadratic(scip, *cons, var2, 0.0, 0.0) );
               }
            }

            SCIP_CALL( SCIPaddBilinTermQuadratic(scip, *cons, var1, var2, monomialcoefs[i]) );
         }
         else
         {
            SCIPerrorMessage("polynomial in quadratic constraint does not have degree at most 2\n");
            *success = FALSE;
            SCIP_CALL( SCIPreleaseCons(scip, cons) );
            break;
         }
      }
   }

   SCIPfreeParseVarsPolynomialData(scip, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, nmonomials);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nlinvars + consdata->nquadvars )
      (*success) = FALSE;
   else
   {
      int i;

      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->linvars, consdata->nlinvars);

      for( i = 0; i < consdata->nquadvars; ++i )
         vars[consdata->nlinvars+i] = consdata->quadvarterms[i].var;

      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nlinvars + consdata->nquadvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for quadratic constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create quadratic constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpQuadratic, consEnfopsQuadratic, consCheckQuadratic, consLockQuadratic,
         conshdlrdata) );
   assert(conshdlr != NULL);


   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyQuadratic, consCopyQuadratic) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteQuadratic) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableQuadratic) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableQuadratic) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitQuadratic) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreQuadratic) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolQuadratic) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeQuadratic) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsQuadratic) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsQuadratic) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitQuadratic) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolQuadratic) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpQuadratic) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseQuadratic) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolQuadratic, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintQuadratic) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropQuadratic, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpQuadratic, consSepasolQuadratic, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransQuadratic) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxQuadratic) );

   /* add quadratic constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/replacebinaryprod",
         "max. length of linear term which when multiplied with a binary variables is replaced by an auxiliary variable and a linear reformulation (0 to turn off)",
         &conshdlrdata->replacebinaryprodlength, FALSE, INT_MAX, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/empathy4and",
         "empathy level for using the AND constraint handler: 0 always avoid using AND; 1 use AND sometimes; 2 use AND as often as possible",
         &conshdlrdata->empathy4and, FALSE, 0, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/binreforminitial",
         "whether to make non-varbound linear constraints added due to replacing products with binary variables initial",
         &conshdlrdata->binreforminitial, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/binreformbinaryonly",
         "whether to consider only binary variables when replacing products with binary variables",
         &conshdlrdata->binreformbinaryonly, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/binreformmaxcoef",
         "limit (as factor on 1/feastol) on coefficients and coef. range in linear constraints created when replacing products with binary variables",
         &conshdlrdata->binreformmaxcoef, TRUE, 1e-4, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/cutmaxrange",
         "maximal coef range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",
         &conshdlrdata->cutmaxrange, TRUE, 1e+7, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/mincurvcollectbilinterms",
         "minimal curvature of constraints to be considered when returning bilinear terms to other plugins",
         &conshdlrdata->mincurvcollectbilinterms, TRUE, 0.8, -SCIPinfinity(scip), SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/linearizeheursol",
         "whether linearizations of convex quadratic constraints should be added to cutpool in a solution found by some heuristic",
         &conshdlrdata->linearizeheursol, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkcurvature",
         "whether multivariate quadratic functions should be checked for convexity/concavity",
         &conshdlrdata->checkcurvature, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkfactorable",
         "whether constraint functions should be checked to be factorable",
         &conshdlrdata->checkfactorable, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/checkquadvarlocks",
         "whether quadratic variables contained in a single constraint should be forced to be at their lower or upper bounds ('d'isable, change 't'ype, add 'b'ound disjunction)",
         &conshdlrdata->checkquadvarlocks, TRUE, 't', "bdt", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/linfeasshift",
         "whether to try to make solutions in check function feasible by shifting a linear variable (esp. useful if constraint was actually objective function)",
         &conshdlrdata->linfeasshift, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxdisaggrsize",
         "maximum number of created constraints when disaggregating a quadratic constraint (<= 1: off)",
         &conshdlrdata->maxdisaggrsize, FALSE, 1, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/disaggrmergemethod",
         "strategy how to merge independent blocks to reach maxdisaggrsize limit (keep 'b'iggest blocks and merge others; keep 's'mallest blocks and merge other; merge small blocks into bigger blocks to reach 'm'ean sizes)",
         &conshdlrdata->disaggrmergemethod, TRUE, 'm', "bms", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a single constraint within one round of SCIP propagation during solve",
         &conshdlrdata->maxproprounds, TRUE, 1, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproproundspresolve",
         "limit on number of propagation rounds for a single constraint within one round of SCIP presolve",
         &conshdlrdata->maxproproundspresolve, TRUE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/enfolplimit",
         "maximum number of enforcement rounds before declaring the LP relaxation infeasible (-1: no limit); WARNING: changing this parameter might lead to incorrect results!",
         &conshdlrdata->enfolplimit, TRUE, -1, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/sepanlpmincont",
         "minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation",
         &conshdlrdata->sepanlpmincont, FALSE, 1.0, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/enfocutsremovable",
         "are cuts added during enforcement removable from the LP in the same node?",
         &conshdlrdata->enfocutsremovable, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/gaugecuts",
         "should convex quadratics generated strong cuts via gauge function?",
         &conshdlrdata->gaugecuts, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/interiorcomputation",
         "how the interior point for gauge cuts should be computed: 'a'ny point per constraint, 'm'ost interior per constraint",
         &conshdlrdata->interiorcomputation, TRUE, 'a', "am", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/projectedcuts",
         "should convex quadratics generated strong cuts via projections?",
         &conshdlrdata->projectedcuts, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/branchscoring",
         "which score to give branching candidates: convexification 'g'ap, constraint 'v'iolation, 'c'entrality of variable value in domain",
         &conshdlrdata->branchscoring, TRUE, 'g', "cgv", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/usebilinineqbranch",
         "should linear inequalities be consindered when computing the branching scores for bilinear terms?",
         &conshdlrdata->usebilinineqbranch, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/minscorebilinterms",
         "minimal required score in order to use linear inequalities for tighter bilinear relaxations",
         &conshdlrdata->minscorebilinterms, FALSE, 0.01, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/bilinineqmaxseparounds",
         "maximum number of separation rounds to use linear inequalities for the bilinear term relaxation in a local node",
         &conshdlrdata->bilinineqmaxseparounds, TRUE, 3, 0, INT_MAX, NULL, NULL) );

   conshdlrdata->eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(conshdlrdata->eventhdlr),CONSHDLR_NAME"_boundchange", "signals a bound change to a quadratic constraint",
         processVarEvent, NULL) );
   assert(conshdlrdata->eventhdlr != NULL);

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, CONSHDLR_NAME"_newsolution", "handles the event that a new primal solution has been found",
         processNewSolutionEvent, NULL) );

   /* include the quadratic constraint upgrade in the nonlinear constraint handler */
   SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, nonlinconsUpgdQuadratic, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );

   return SCIP_OKAY;
}

/** includes a quadratic constraint update method into the quadratic constraint handler */
SCIP_RETCODE SCIPincludeQuadconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd)),  /**< method to call for upgrading quadratic constraint */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method be active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*        conshdlr;
   SCIP_CONSHDLRDATA*    conshdlrdata;
   SCIP_QUADCONSUPGRADE* quadconsupgrade;
   char                  paramname[SCIP_MAXSTRLEN];
   char                  paramdesc[SCIP_MAXSTRLEN];
   int                   i;

   assert(quadconsupgd != NULL);
   assert(conshdlrname != NULL );

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdataHasUpgrade(scip, conshdlrdata, quadconsupgd, conshdlrname) )
   {
      /* create a quadratic constraint upgrade data object */
      SCIP_CALL( SCIPallocBlockMemory(scip, &quadconsupgrade) );
      quadconsupgrade->quadconsupgd = quadconsupgd;
      quadconsupgrade->priority     = priority;
      quadconsupgrade->active       = active;

      /* insert quadratic constraint upgrade method into constraint handler data */
      assert(conshdlrdata->nquadconsupgrades <= conshdlrdata->quadconsupgradessize);
      if( conshdlrdata->nquadconsupgrades+1 > conshdlrdata->quadconsupgradessize )
      {
         int newsize;

         newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->nquadconsupgrades+1);
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->quadconsupgrades, conshdlrdata->quadconsupgradessize, newsize) );
         conshdlrdata->quadconsupgradessize = newsize;
      }
      assert(conshdlrdata->nquadconsupgrades+1 <= conshdlrdata->quadconsupgradessize);

      for( i = conshdlrdata->nquadconsupgrades; i > 0 && conshdlrdata->quadconsupgrades[i-1]->priority < quadconsupgrade->priority; --i )
         conshdlrdata->quadconsupgrades[i] = conshdlrdata->quadconsupgrades[i-1];
      assert(0 <= i && i <= conshdlrdata->nquadconsupgrades);
      conshdlrdata->quadconsupgrades[i] = quadconsupgrade;
      conshdlrdata->nquadconsupgrades++;

      /* adds parameter to turn on and off the upgrading step */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/" CONSHDLR_NAME "/upgrade/%s", conshdlrname);
      (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable quadratic upgrading for constraint handler <%s>", conshdlrname);
      SCIP_CALL( SCIPaddBoolParam(scip,
            paramname, paramdesc,
            &quadconsupgrade->active, FALSE, active, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** Creates and captures a quadratic constraint.
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_j z_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation (u) */
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
   SCIP_CONSDATA* consdata;
   SCIP_HASHMAP*  quadvaridxs;
   SCIP_Real      sqrcoef;
   int i;
   int var1pos;
   int var2pos;

   int nbilinterms;

   assert(linvars != NULL || nlinvars == 0);
   assert(lincoefs != NULL || nlinvars == 0);
   assert(quadvars1 != NULL || nquadterms == 0);
   assert(quadvars2 != NULL || nquadterms == 0);
   assert(quadcoefs != NULL || nquadterms == 0);

   assert(modifiable == FALSE); /* we do not support column generation */

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data and constraint */
   SCIP_CALL( consdataCreateEmpty(scip, &consdata) );

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   /* add quadratic variables and remember their indices */
   SCIP_CALL( SCIPhashmapCreate(&quadvaridxs, SCIPblkmem(scip), nquadterms) );
   nbilinterms = 0;
   for( i = 0; i < nquadterms; ++i )
   {
      if( SCIPisZero(scip, quadcoefs[i]) )  /*lint !e613*/
         continue;

      /* if it is actually a square term, remember it's coefficient */
      /* cppcheck-suppress nullPointer */
      if( quadvars1[i] == quadvars2[i] )   /*lint !e613*/
         sqrcoef = quadcoefs[i];   /*lint !e613 */
      else
         sqrcoef = 0.0;

      /* add quadvars1[i], if not in there already */
      if( !SCIPhashmapExists(quadvaridxs, quadvars1[i]) )  /*lint !e613*/
      {
         SCIP_CALL( addQuadVarTerm(scip, *cons, quadvars1[i], 0.0, sqrcoef) );   /*lint !e613*/
         assert(consdata->nquadvars >= 0);
         assert(consdata->quadvarterms[consdata->nquadvars-1].var == quadvars1[i]);  /*lint !e613*/

         SCIP_CALL( SCIPhashmapInsert(quadvaridxs, quadvars1[i], (void*)(size_t)(consdata->nquadvars-1)) );   /*lint !e613*/
      }
      else if( !SCIPisZero(scip, sqrcoef) )
      {
         /* if it's there already, but we got a square coefficient, add it to the previous one */
         var1pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, quadvars1[i]);   /*lint !e613*/
         assert(consdata->quadvarterms[var1pos].var == quadvars1[i]);   /*lint !e613*/
         consdata->quadvarterms[var1pos].sqrcoef += sqrcoef;
      }

      /* cppcheck-suppress nullPointer */
      if( quadvars1[i] == quadvars2[i] )  /*lint !e613*/
         continue;

      /* add quadvars2[i], if not in there already */
      if( !SCIPhashmapExists(quadvaridxs, quadvars2[i]) )   /*lint !e613*/
      {
         assert(sqrcoef == 0.0);
         SCIP_CALL( addQuadVarTerm(scip, *cons, quadvars2[i], 0.0, 0.0) );   /*lint !e613*/
         assert(consdata->nquadvars >= 0);
         assert(consdata->quadvarterms[consdata->nquadvars-1].var == quadvars2[i]);  /*lint !e613*/

         SCIP_CALL( SCIPhashmapInsert(quadvaridxs, quadvars2[i], (void*)(size_t)(consdata->nquadvars-1)) );  /*lint !e613*/
      }

      ++nbilinterms;
   }

   /* add bilinear terms, if we saw any */
   if( nbilinterms > 0 )
   {
      SCIP_CALL( consdataEnsureBilinSize(scip, consdata, nbilinterms) );
      for( i = 0; i < nquadterms; ++i )
      {
         if( SCIPisZero(scip, quadcoefs[i]) )  /*lint !e613*/
            continue;

         /* square terms have been taken care of already */
         if( quadvars1[i] == quadvars2[i] )  /*lint !e613 */
            continue;

         assert(SCIPhashmapExists(quadvaridxs, quadvars1[i]));  /*lint !e613*/
         assert(SCIPhashmapExists(quadvaridxs, quadvars2[i]));  /*lint !e613*/

         var1pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, quadvars1[i]);  /*lint !e613*/
         var2pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, quadvars2[i]);  /*lint !e613*/

         SCIP_CALL( addBilinearTerm(scip, *cons, var1pos, var2pos, quadcoefs[i]) );  /*lint !e613*/
      }
   }

   /* add linear variables */
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, nlinvars) );
   for( i = 0; i < nlinvars; ++i )
   {
      if( SCIPisZero(scip, lincoefs[i]) )  /*lint !e613*/
         continue;

      /* if it's a linear coefficient for a quadratic variable, add it there, otherwise add as linear variable */
      if( SCIPhashmapExists(quadvaridxs, linvars[i]) )  /*lint !e613*/
      {
         var1pos = (int) (size_t) SCIPhashmapGetImage(quadvaridxs, linvars[i]);  /*lint !e613*/
         assert(consdata->quadvarterms[var1pos].var == linvars[i]);  /*lint !e613*/
         consdata->quadvarterms[var1pos].lincoef += lincoefs[i];  /*lint !e613*/
      }
      else
      {
         SCIP_CALL( addLinearCoef(scip, *cons, linvars[i], lincoefs[i]) );  /*lint !e613*/
      }
   }

   SCIPhashmapFree(&quadvaridxs);

   SCIPdebugMsg(scip, "created quadratic constraint ");
   SCIPdebugPrintCons(scip, *cons, NULL);

   return SCIP_OKAY;
}

/** creates and captures a quadratic constraint with all its
 *  flags set to their default values.
 *
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_j z_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms (a_j) */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (ell) */
   SCIP_Real             rhs                 /**< right hand side of quadratic equation (u) */
   )
{
   SCIP_CALL( SCIPcreateConsQuadratic(scip, cons, name, nlinvars, linvars, lincoefs,
         nquadterms, quadvars1, quadvars2, quadcoefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** Creates and captures a quadratic constraint.
 *
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_k v_k w_k \leq u.
 * \f]
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadvarterms,      /**< number of quadratic terms (m) */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms (p) */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             lhs,                /**< constraint left hand side (ell) */
   SCIP_Real             rhs,                /**< constraint right hand side (u) */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(modifiable == FALSE); /* we do not support column generation */
   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadvarterms == 0 || quadvarterms != NULL);
   assert(nbilinterms == 0 || bilinterms != NULL);

   /* find the quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("quadratic constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, lhs, rhs,
         nlinvars, linvars, lincoefs, nquadvarterms, quadvarterms, nbilinterms, bilinterms,
         TRUE) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   return SCIP_OKAY;
}

/** creates and captures a quadratic constraint in its most basic version, i.e.,
 *  all constraint flags are set to their default values.
 *
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_k v_k w_k \leq u.
 * \f]
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) */
   int                   nquadvarterms,      /**< number of quadratic terms (m) */
   SCIP_QUADVARTERM*     quadvarterms,       /**< quadratic variable terms */
   int                   nbilinterms,        /**< number of bilinear terms (p) */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms */
   SCIP_Real             lhs,                /**< constraint left hand side (ell) */
   SCIP_Real             rhs                 /**< constraint right hand side (u) */
   )
{
   SCIP_CALL( SCIPcreateConsQuadratic2(scip, cons, name, nlinvars, linvars, lincoefs,
         nquadvarterms, quadvarterms, nbilinterms, bilinterms, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}


/** Adds a constant to the constraint function, that is, subtracts a constant from both sides */
void SCIPaddConstantQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             constant            /**< constant to subtract from both sides */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, REALABS(constant)));

   /* nlrow and solving data (see initsol) may become invalid when changing constraint */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsEnabled(cons) )
   {
      SCIPerrorMessage("Cannot modify enabled constraint in solving stage.\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lhs <= consdata->rhs);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
      consdata->lhs -= constant;
   if( !SCIPisInfinity(scip,  consdata->rhs) )
      consdata->rhs -= constant;

   if( consdata->lhs > consdata->rhs )
   {
      assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));
      consdata->lhs = consdata->rhs;
   }
}

/** Adds a linear variable with coefficient to a quadratic constraint. */
SCIP_RETCODE SCIPaddLinearVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient of variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   /* nlrow and solving data (see initsol) may become invalid when changing constraint */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsEnabled(cons) )
   {
      SCIPerrorMessage("Cannot modify enabled constraint in solving stage.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( addLinearCoef(scip, cons, var, coef) );

   return SCIP_OKAY;
}

/** Adds a quadratic variable with linear and square coefficient to a quadratic constraint. */
SCIP_RETCODE SCIPaddQuadVarQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             lincoef,            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef             /**< square coefficient of variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(lincoef)));
   assert(!SCIPisInfinity(scip, REALABS(sqrcoef)));

   /* nlrow and solving data (see initsol) may become invalid when changing constraint */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsEnabled(cons) )
   {
      SCIPerrorMessage("Cannot modify enabled constraint in solving stage.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( addQuadVarTerm(scip, cons, var, lincoef, sqrcoef) );

   return SCIP_OKAY;
}

/** Adds a linear coefficient for a quadratic variable.
 *
 *  Variable will be added with square coefficient 0.0 if not existing yet.
 */
SCIP_RETCODE SCIPaddQuadVarLinearCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to linear coefficient of variable */
   )
{
   SCIP_CONSDATA* consdata;
   int pos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   if( SCIPisZero(scip, coef) )
      return SCIP_OKAY;

   /* nlrow and solving data (see initsol) may become invalid when changing constraint */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsEnabled(cons) )
   {
      SCIPerrorMessage("Cannot modify enabled constraint in solving stage.\n");
      return SCIP_INVALIDCALL;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var, &pos) );
   if( pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var, coef, 0.0) );
      return SCIP_OKAY;
   }
   assert(pos < consdata->nquadvars);
   assert(consdata->quadvarterms[pos].var == var);

   consdata->quadvarterms[pos].lincoef += coef;

   /* update flags and invalid activities */
   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = consdata->ispresolved && !SCIPisZero(scip, consdata->quadvarterms[pos].lincoef);

   SCIPintervalSetEmpty(&consdata->quadactivitybounds);
   consdata->activity = SCIP_INVALID;

   return SCIP_OKAY;
}

/** Adds a square coefficient for a quadratic variable.
 *
 *  Variable will be added with linear coefficient 0.0 if not existing yet.
 */
SCIP_RETCODE SCIPaddSquareCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< value to add to square coefficient of variable */
   )
{
   SCIP_CONSDATA* consdata;
   int pos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   if( SCIPisZero(scip, coef) )
      return SCIP_OKAY;

   /* nlrow and solving data (see initsol) may become invalid when changing constraint */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsEnabled(cons) )
   {
      SCIPerrorMessage("Cannot modify enabled constraint in solving stage.\n");
      return SCIP_INVALIDCALL;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var, &pos) );
   if( pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var, 0.0, coef) );
      return SCIP_OKAY;
   }
   assert(pos < consdata->nquadvars);
   assert(consdata->quadvarterms[pos].var == var);

   consdata->quadvarterms[pos].sqrcoef += coef;

   /* update flags and invalid activities */
   consdata->isconvex      = FALSE;
   consdata->isconcave     = FALSE;
   consdata->iscurvchecked = FALSE;
   consdata->ispropagated  = FALSE;
   consdata->ispresolved   = consdata->ispresolved && !SCIPisZero(scip, consdata->quadvarterms[pos].sqrcoef);

   SCIPintervalSetEmpty(&consdata->quadactivitybounds);
   consdata->activity = SCIP_INVALID;

   return SCIP_OKAY;
}

/** Adds a bilinear term to a quadratic constraint.
 *
 *  Variables will be added with linear and square coefficient 0.0 if not existing yet.
 *  If variables are equal, only the square coefficient of the variable is updated.
 */
SCIP_RETCODE SCIPaddBilinTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   )
{
   SCIP_CONSDATA* consdata;
   int            var1pos;
   int            var2pos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   /* nlrow and solving data (see initsol) may become invalid when changing constraint */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPconsIsEnabled(cons) )
   {
      SCIPerrorMessage("Cannot modify enabled constraint in solving stage.\n");
      return SCIP_INVALIDCALL;
   }

   if( var1 == var2 )
   {
      SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, var1, coef) );
      return SCIP_OKAY;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var1, &var1pos) );
   if( var1pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var1, 0.0, 0.0) );
      var1pos = consdata->nquadvars-1;
   }

   if( !consdata->quadvarssorted )
   {
      SCIP_CALL( consdataSortQuadVarTerms(scip, consdata) );
      /* sorting may change the position of var1 */
      SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var1, &var1pos) );
      assert(var1pos >= 0);
   }

   assert(consdata->quadvarssorted);
   SCIP_CALL( consdataFindQuadVarTerm(scip, consdata, var2, &var2pos) );
   if( var2pos < 0 )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var2, 0.0, 0.0) );
      var2pos = consdata->nquadvars-1;
   }

   assert(consdata->quadvarterms[var1pos].var == var1);
   assert(consdata->quadvarterms[var2pos].var == var2);

   SCIP_CALL( addBilinearTerm(scip, cons, var1pos, var2pos, coef) );

   return SCIP_OKAY;
}

/** Gets the quadratic constraint as a nonlinear row representation. */
SCIP_RETCODE SCIPgetNlRowQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);

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

/** Gets the number of variables in the linear term of a quadratic constraint. */
int SCIPgetNLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nlinvars;
}

/** Gets the variables in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_VAR** SCIPgetLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->linvars;
}

/** Gets the coefficients in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
SCIP_Real* SCIPgetCoefsLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lincoefs;
}

/** Gets the number of quadratic variable terms of a quadratic constraint.
 */
int SCIPgetNQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nquadvars;
}

/** Gets the quadratic variable terms of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarTermsQuadratic.
 */
SCIP_QUADVARTERM* SCIPgetQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->quadvarterms;
}

/** Ensures that quadratic variable terms are sorted. */
SCIP_RETCODE SCIPsortQuadVarTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   SCIP_CALL( consdataSortQuadVarTerms(scip, SCIPconsGetData(cons)) );

   return SCIP_OKAY;
}

/** Finds the position of a quadratic variable term for a given variable.
 *
 *  @note If the quadratic variable terms have not been sorted before, then a search may reorder the current order of the terms.
 */
SCIP_RETCODE SCIPfindQuadVarTermQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable to search for */
   int*                  pos                 /**< buffer to store position of quadvarterm for var, or -1 if not found */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(var != NULL);
   assert(pos != NULL);

   SCIP_CALL( consdataFindQuadVarTerm(scip, SCIPconsGetData(cons), var, pos) );

   return SCIP_OKAY;
}

/** Gets the number of bilinear terms of a quadratic constraint. */
int SCIPgetNBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nbilinterms;
}

/** Gets the bilinear terms of a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
SCIP_BILINTERM* SCIPgetBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->bilinterms;
}

/** Gets the left hand side of a quadratic constraint. */
SCIP_Real SCIPgetLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lhs;
}

/** Gets the right hand side of a quadratic constraint. */
SCIP_Real SCIPgetRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhs;
}

/** get index of a variable in linvars that may be decreased without making any other constraint infeasible, or -1 if none */
int SCIPgetLinvarMayDecreaseQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA*     consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check for a linear variable that can be increase or decreased without harming feasibility */
   consdataFindUnlockedLinearVar(scip, consdata);

   return consdata->linvar_maydecrease;
}

/** get index of a variable in linvars that may be increased without making any other constraint infeasible, or -1 if none */
int SCIPgetLinvarMayIncreaseQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA*     consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check for a linear variable that can be increase or decreased without harming feasibility */
   consdataFindUnlockedLinearVar(scip, consdata);

   return consdata->linvar_mayincrease;
}

/** Check the quadratic function of a quadratic constraint for its semi-definiteness, if not done yet. */
SCIP_RETCODE SCIPcheckCurvatureQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);

   SCIP_CALL( checkCurvature(scip, cons, TRUE) );

   return SCIP_OKAY;
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) convex. */
SCIP_Bool SCIPisConvexQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_Bool determined;

   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   checkCurvatureEasy(scip, cons, &determined, FALSE);
   assert(determined);

   return (SCIPconsGetData(cons)->isconvex);
}

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) concave. */
SCIP_Bool SCIPisConcaveQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_Bool determined;

   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   checkCurvatureEasy(scip, cons, &determined, FALSE);
   assert(determined);

   return (SCIPconsGetData(cons)->isconcave);
}

/** Computes the violation of a constraint by a solution */
SCIP_RETCODE SCIPgetViolationQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< pointer to store violation of constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool solviolbounds;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violation != NULL);

   SCIP_CALL( computeViolation(scip, cons, sol, &solviolbounds) );
   /* we don't care here whether the solution violated variable bounds */

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violation = MAX(consdata->lhsviol, consdata->rhsviol);

   return SCIP_OKAY;
}

/** Indicates whether the quadratic constraint is local w.r.t. the current local bounds.
 *
 *  That is, checks whether each variable with a square term is fixed and for each bilinear term at least one variable is fixed.
 */
SCIP_Bool SCIPisLinearLocalQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check all square terms */
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( consdata->quadvarterms[i].sqrcoef == 0.0 )
         continue;

      var1 = consdata->quadvarterms[i].var;
      assert(var1 != NULL);

      if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) )
         return FALSE;
   }

   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      var1 = consdata->bilinterms[i].var1;
      var2 = consdata->bilinterms[i].var2;

      assert(var1 != NULL);
      assert(var2 != NULL);

      if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) &&
         ! SCIPisRelEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
         return FALSE;
   }

   return TRUE;
}

/** Adds the constraint to an NLPI problem. */
SCIP_RETCODE SCIPaddToNlpiProblemQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< NLPI problem where to add constraint */
   SCIP_HASHMAP*         scipvar2nlpivar,    /**< mapping from SCIP variables to variable indices in NLPI */
   SCIP_Bool             names               /**< whether to pass constraint names to NLPI */
   )
{
   SCIP_CONSDATA* consdata;
   int            nlininds;
   int*           lininds;
   SCIP_Real*     linvals;
   int            nquadelems;
   SCIP_QUADELEM* quadelems;
   SCIP_VAR*      othervar;
   const char*    name;
   int            j;
   int            l;
   int            lincnt;
   int            quadcnt;
   int            idx1;
   int            idx2;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nlpi != NULL);
   assert(nlpiprob != NULL);
   assert(scipvar2nlpivar != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* count nonzeros in quadratic part */
   nlininds = consdata->nlinvars;
   nquadelems = consdata->nbilinterms;
   for( j = 0; j < consdata->nquadvars; ++j )
   {
      if( consdata->quadvarterms[j].sqrcoef != 0.0 )
         ++nquadelems;
      if( consdata->quadvarterms[j].lincoef != 0.0 )
         ++nlininds;
   }

   /* setup linear part */
   lininds = NULL;
   linvals = NULL;
   lincnt  = 0;
   if( nlininds > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nlininds) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nlininds) );

      for( j = 0; j < consdata->nlinvars; ++j )
      {
         linvals[j] = consdata->lincoefs[j];
         assert(SCIPhashmapExists(scipvar2nlpivar, consdata->linvars[j]));
         lininds[j] = (int) (size_t) SCIPhashmapGetImage(scipvar2nlpivar, consdata->linvars[j]);
      }

      lincnt = consdata->nlinvars;
   }

   /* setup quadratic part */
   quadelems = NULL;
   if( nquadelems > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nquadelems) );
   }
   quadcnt = 0;

   for( j = 0; j < consdata->nquadvars; ++j )
   {
      assert(SCIPhashmapExists(scipvar2nlpivar, consdata->quadvarterms[j].var));
      idx1 = (int)(size_t)SCIPhashmapGetImage(scipvar2nlpivar, consdata->quadvarterms[j].var);
      if( consdata->quadvarterms[j].lincoef != 0.0 )
      {
         assert(lininds != NULL);
         assert(linvals != NULL);
         lininds[lincnt] = idx1;
         linvals[lincnt] = consdata->quadvarterms[j].lincoef;
         ++lincnt;
      }

      if( consdata->quadvarterms[j].sqrcoef != 0.0 )
      {
         assert(quadcnt < nquadelems);
         assert(quadelems != NULL);
         quadelems[quadcnt].idx1 = idx1;
         quadelems[quadcnt].idx2 = idx1;
         quadelems[quadcnt].coef = consdata->quadvarterms[j].sqrcoef;
         ++quadcnt;
      }

      for( l = 0; l < consdata->quadvarterms[j].nadjbilin; ++l )
      {
         othervar = consdata->bilinterms[consdata->quadvarterms[j].adjbilin[l]].var2;
         /* if othervar is on position 2, then we process this bilinear term later (or it was processed already) */
         if( othervar == consdata->quadvarterms[j].var )
            continue;

         assert(quadcnt < nquadelems);
         assert(quadelems != NULL);
         assert(SCIPhashmapExists(scipvar2nlpivar, othervar));
         idx2 = (int)(size_t)SCIPhashmapGetImage(scipvar2nlpivar, othervar);
         quadelems[quadcnt].idx1 = MIN(idx1, idx2);
         quadelems[quadcnt].idx2 = MAX(idx1, idx2);
         quadelems[quadcnt].coef = consdata->bilinterms[consdata->quadvarterms[j].adjbilin[l]].coef;
         ++quadcnt;
      }
   }

   assert(quadcnt == nquadelems);
   assert(lincnt  == nlininds);

   name = names ? SCIPconsGetName(cons) : NULL;

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1,
         &consdata->lhs, &consdata->rhs,
         &nlininds, &lininds, &linvals ,
         &nquadelems, &quadelems,
         NULL, NULL, &name) );

   SCIPfreeBufferArrayNull(scip, &quadelems);
   SCIPfreeBufferArrayNull(scip, &lininds);
   SCIPfreeBufferArrayNull(scip, &linvals);

   return SCIP_OKAY;
}


/** sets the left hand side of a quadratic constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint.
 */
SCIP_RETCODE SCIPchgLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, lhs));

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not quadratic\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM || !SCIPconsIsOriginal(cons) )
   {
      SCIPerrorMessage("method may only be called during problem creation stage for original constraints\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* adjust value to not be smaller than -inf */
   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);

   /* check for lhs <= rhs */
   if( !SCIPisLE(scip, lhs, consdata->rhs) )
      return SCIP_INVALIDDATA;

   consdata->lhs = lhs;

   return SCIP_OKAY;
}

/** sets the right hand side of a quadratic constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint.
 */
SCIP_RETCODE SCIPchgRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, -rhs));

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not quadratic\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM || !SCIPconsIsOriginal(cons) )
   {
      SCIPerrorMessage("method may only be called during problem creation stage for original constraints\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* adjust value to not be greater than inf */
   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);

   /* check for lhs <= rhs */
   if( !SCIPisLE(scip, consdata->lhs, rhs) )
      return SCIP_INVALIDDATA;

   consdata->rhs = rhs;

   return SCIP_OKAY;
}

/** gets the feasibility of the quadratic constraint in the given solution */
SCIP_RETCODE SCIPgetFeasibilityQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_Real*            feasibility         /**< pointer to store the feasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool solviolbounds;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(feasibility != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not quadratic\n");
      SCIPABORT();
   }

   SCIP_CALL( computeViolation(scip, cons, sol, &solviolbounds) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, -consdata->lhs) )
      *feasibility = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -consdata->lhs) )
      *feasibility = (consdata->rhs - consdata->activity);
   else if( SCIPisInfinity(scip, consdata->rhs) )
      *feasibility = (consdata->activity - consdata->lhs);
   else
   {
      assert(!SCIPisInfinity(scip, -consdata->rhs));
      assert(!SCIPisInfinity(scip, consdata->lhs));
      *feasibility = MIN( consdata->rhs - consdata->activity, consdata->activity - consdata->lhs );
   }

   return SCIP_OKAY;
}

/** gets the activity of the quadratic constraint in the given solution */
SCIP_RETCODE SCIPgetActivityQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol,                /**< solution, or NULL to use current node's solution */
   SCIP_Real*            activity            /**< pointer to store the activity */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool solviolbounds;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(activity != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not quadratic\n");
      SCIPABORT();
   }

   SCIP_CALL( computeViolation(scip, cons, sol, &solviolbounds) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *activity = consdata->activity;

   return SCIP_OKAY;
}

/** changes the linear coefficient value for a given quadratic variable in a quadratic constraint data; if not
 *  available, it adds it
 *
 *  @note this is only allowed for original constraints and variables in problem creation stage
 */
SCIP_RETCODE SCIPchgLinearCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< quadratic variable */
   SCIP_Real             coef                /**< new coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool found;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0  )
   {
      SCIPerrorMessage("constraint is not quadratic\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM || !SCIPconsIsOriginal(cons) || !SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("method may only be called during problem creation stage for original constraints and variables\n");
      return SCIP_INVALIDDATA;
   }

   consdata =  SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check all quadratic variables */
   found = FALSE;
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( var == consdata->quadvarterms[i].var )
      {
         if( found || SCIPisZero(scip, coef) )
         {
            consdata->quadvarterms[i].lincoef = 0.0;

            /* remember to merge quadratic variable terms */
            consdata->quadvarsmerged = FALSE;
         }
         else
            consdata->quadvarterms[i].lincoef = coef;

         found = TRUE;
      }
   }

   /* check all linear variables */
   i = 0;
   while( i < consdata->nlinvars )
   {
      if( var == consdata->linvars[i] )
      {
         if( found || SCIPisZero(scip, coef) )
         {
            SCIP_CALL( delLinearCoefPos(scip, cons, i) );

            /* decrease i by one since otherwise we would skip the coefficient which has been switched to position i */
            i--;
         }
         else
         {
            SCIP_CALL( chgLinearCoefPos(scip, cons, i, coef) );
         }

         found = TRUE;
      }
      i++;
   }

   /* add linear term if necessary */
   if( !found && !SCIPisZero(scip, coef) )
   {
      SCIP_CALL( addLinearCoef(scip, cons, var, coef) );
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;

   SCIPintervalSetEmpty(&consdata->quadactivitybounds);
   consdata->activity = SCIP_INVALID;

   return SCIP_OKAY;
}

/** changes the square coefficient value for a given quadratic variable in a quadratic constraint data; if not
 *  available, it adds it
 *
 *  @note this is only allowed for original constraints and variables in problem creation stage
 */
SCIP_RETCODE SCIPchgSquareCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< quadratic variable */
   SCIP_Real             coef                /**< new coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool found;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not quadratic\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM || !SCIPconsIsOriginal(cons) || !SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("method may only be called during problem creation stage for original constraints and variables\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* find the quadratic variable and change its quadratic coefficient */
   found = FALSE;
   for( i = 0; i < consdata->nquadvars; ++i )
   {
      if( var == consdata->quadvarterms[i].var )
      {
         consdata->quadvarterms[i].sqrcoef = (found || SCIPisZero(scip, coef)) ? 0.0 : coef;
         found = TRUE;
      }
   }

   /* add bilinear term if necessary */
   if( !found && !SCIPisZero(scip, coef) )
   {
      SCIP_CALL( addQuadVarTerm(scip, cons, var, 0.0, coef) );
   }

   /* update flags and invalidate activities */
   consdata->isconvex = FALSE;
   consdata->isconcave = FALSE;
   consdata->iscurvchecked = FALSE;
   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;

   SCIPintervalSetEmpty(&consdata->quadactivitybounds);
   consdata->activity = SCIP_INVALID;

   /* remember to merge quadratic variable terms */
   consdata->quadvarsmerged = FALSE;

   return SCIP_OKAY;
}

/** changes the bilinear coefficient value for a given quadratic variable in a quadratic constraint data; if not
 *  available, it adds it
 *
 *  @note this is only allowed for original constraints and variables in problem creation stage
 */
SCIP_RETCODE SCIPchgBilinCoefQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var1,               /**< first variable */
   SCIP_VAR*             var2,               /**< second variable */
   SCIP_Real             coef                /**< coefficient of bilinear term */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool found;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var1 != NULL);
   assert(var2 != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not quadratic\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM || !SCIPconsIsOriginal(cons) || !SCIPvarIsOriginal(var1) || !SCIPvarIsOriginal(var2) )
   {
      SCIPerrorMessage("method may only be called during problem creation stage for original constraints and variables\n");
      return SCIP_INVALIDDATA;
   }

   if( var1 == var2 )
   {
      SCIP_CALL( SCIPchgSquareCoefQuadratic(scip, cons, var1, coef) );
      return SCIP_OKAY;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* search array of bilinear terms */
   found = FALSE;
   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      if( (consdata->bilinterms[i].var1 == var1 && consdata->bilinterms[i].var2 == var2) ||
            (consdata->bilinterms[i].var1 == var2 && consdata->bilinterms[i].var2 == var1)  )
      {
         if( found || SCIPisZero(scip, coef) )
         {
            consdata->bilinterms[i].coef = 0.0;

            /* remember to merge bilinear terms */
            consdata->bilinmerged = FALSE;
         }
         else
            consdata->bilinterms[i].coef = coef;
         found = TRUE;
      }
   }

   /* add bilinear term if necessary */
   if( !found )
   {
      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, var1, var2, coef) );
   }

   /* update flags and invalidate activities */
   consdata->isconvex = FALSE;
   consdata->isconcave = FALSE;
   consdata->iscurvchecked = FALSE;
   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;

   SCIPintervalSetEmpty(&consdata->quadactivitybounds);
   consdata->activity = SCIP_INVALID;

   return SCIP_OKAY;
}

/** returns the total number of bilinear terms that are contained in all quadratic constraints */
int SCIPgetNAllBilinearTermsQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
      return 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nbilinterms;
}

/** returns all bilinear terms that are contained in all quadratic constraints */
SCIP_RETCODE SCIPgetAllBilinearTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR** RESTRICT   x,                  /**< array to store first variable of each bilinear term */
   SCIP_VAR** RESTRICT   y,                  /**< array to second variable of each bilinear term */
   int* RESTRICT         nbilinterms,        /**< buffer to store the total number of bilinear terms */
   int* RESTRICT         nunderests,         /**< array to store the total number of constraints that require to underestimate a bilinear term */
   int* RESTRICT         noverests,          /**< array to store the total number of constraints that require to overestimate a bilinear term */
   SCIP_Real*            maxnonconvexity     /**< largest absolute value of nonconvex eigenvalues of all quadratic constraints containing a bilinear term */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   int i;

   assert(scip != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(nbilinterms != NULL);
   assert(nunderests != NULL);
   assert(noverests!= NULL);
   assert(maxnonconvexity != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      *nbilinterms = 0;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nbilinterms; ++i )
   {
      x[i] = conshdlrdata->bilinestimators[i].x;
      y[i] = conshdlrdata->bilinestimators[i].y;
      nunderests[i] = conshdlrdata->bilinestimators[i].nunderest;
      noverests[i] = conshdlrdata->bilinestimators[i].noverest;
      maxnonconvexity[i] = conshdlrdata->bilinestimators[i].maxnonconvexity;
   }

   *nbilinterms = conshdlrdata->nbilinterms;

   return SCIP_OKAY;
}

/** helper function to compute the violation of an inequality of the form xcoef * x <= ycoef * y + constant for two
 *  corner points of the domain [lbx,ubx]x[lby,uby]
 */
static
void getIneqViol(
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             xcoef,              /**< x-coefficient */
   SCIP_Real             ycoef,              /**< y-coefficient */
   SCIP_Real             constant,           /**< constant */
   SCIP_Real*            viol1,              /**< buffer to store the violation of the first corner point */
   SCIP_Real*            viol2               /**< buffer to store the violation of the second corner point */
   )
{
   SCIP_Real norm;

   assert(viol1 != NULL);
   assert(viol2 != NULL);

   norm = SQRT(SQR(xcoef) + SQR(ycoef));

   /* inequality can be used for underestimating xy if and only if xcoef * ycoef > 0 */
   if( xcoef * ycoef >= 0 )
   {
      /* violation for top-left and bottom-right corner */
      *viol1 = MAX(0, (xcoef * SCIPvarGetLbLocal(x)  - ycoef * SCIPvarGetUbLocal(y) - constant) / norm); /*lint !e666*/
      *viol2 = MAX(0, (xcoef * SCIPvarGetUbLocal(x)  - ycoef * SCIPvarGetLbLocal(y) - constant) / norm); /*lint !e666*/
   }
   else
   {
      /* violation for top-right and bottom-left corner */
      *viol1 = MAX(0, (xcoef * SCIPvarGetUbLocal(x)  - ycoef * SCIPvarGetUbLocal(y) - constant) / norm); /*lint !e666*/
      *viol2 = MAX(0, (xcoef * SCIPvarGetLbLocal(x)  - ycoef * SCIPvarGetLbLocal(y) - constant) / norm); /*lint !e666*/
   }

   return;
}

/** adds a globally valid inequality of the form xcoef x <= ycoef y + constant for a bilinear term (x,y)
 *
 *  @note the indices of bilinear terms match with the entries of bilinear terms returned by SCIPgetAllBilinearTermsQuadratic
 */
SCIP_RETCODE SCIPaddBilinearIneqQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   int                   idx,                /**< index of the bilinear term */
   SCIP_Real             xcoef,              /**< x coefficient */
   SCIP_Real             ycoef,              /**< y coefficient */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool*            success             /**< buffer to store whether inequality has been accepted */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   BILINESTIMATOR* bilinest;
   SCIP_Real* ineqs;
   SCIP_Real viol1 = 0.0;
   SCIP_Real viol2 = 0.0;
   int* nineqs;
   int i;

   assert(scip != NULL);
   assert(x != NULL);
   assert(y != NULL);
   assert(idx >= 0);
   assert(xcoef != SCIP_INVALID); /*lint !e777 */
   assert(ycoef != SCIP_INVALID); /*lint !e777 */
   assert(constant != SCIP_INVALID); /*lint !e777 */
   assert(success != NULL);

   *success = FALSE;

   /* ignore inequalities that only yield to a (possible) bound tightening */
   if( SCIPisFeasZero(scip, xcoef) || SCIPisFeasZero(scip, ycoef) )
      return SCIP_OKAY;

   /* get constraint handler and its data */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(idx < conshdlrdata->nbilinterms);

   bilinest = &conshdlrdata->bilinestimators[idx];
   assert(bilinest != NULL);
   assert(bilinest->x == x);
   assert(bilinest->y == y);

   SCIPdebugMsg(scip, "add bilinear term inequality: %g %s <= %g %s + %g\n", xcoef, SCIPvarGetName(bilinest->x),
      ycoef, SCIPvarGetName(bilinest->y), constant);

   if( xcoef * ycoef > 0.0 )
   {
      ineqs = bilinest->inequnderest;
      nineqs = &bilinest->ninequnderest;
   }
   else
   {
      ineqs = bilinest->ineqoverest;
      nineqs = &bilinest->nineqoverest;
   }

   /* compute violation of the inequality of the important corner points */
   getIneqViol(x, y, xcoef, ycoef, constant, &viol1, &viol2);
   SCIPdebugMsg(scip, "violations of inequality = (%g,%g)\n", viol1, viol2);

   /* inequality does not cut off one of the important corner points */
   if( SCIPisFeasLE(scip, MAX(viol1, viol2), 0.0) )
      return SCIP_OKAY;

   /* check whether inequality exists already */
   for( i = 0; i < *nineqs; ++i )
   {
      if( SCIPisFeasEQ(scip, xcoef, ineqs[3*i]) && SCIPisFeasEQ(scip, ycoef, ineqs[3*i+1])
         && SCIPisFeasEQ(scip, constant, ineqs[3*i+2]) )
      {
         SCIPdebugMsg(scip, "inequality already found -> skip\n");
         return SCIP_OKAY;
      }
   }

   /* add inequality if we found less than two so far; otherwise compare the violations to decide which which
    * inequality might be replaced
    */
   if( *nineqs < 2 )
   {
      ineqs[3*(*nineqs)] = xcoef;
      ineqs[3*(*nineqs) + 1] = ycoef;
      ineqs[3*(*nineqs) + 2] = constant;
      ++(*nineqs);
      *success = TRUE;
   }
   else
   {
      SCIP_Real viols1[2] = {0.0, 0.0};
      SCIP_Real viols2[2] = {0.0, 0.0};
      SCIP_Real bestviol;
      int pos = -1;

      assert(*nineqs == 2);

      /* compute resulting violations of both corner points when replacing an existing inequality
       *
       * given the violations (v1,w1), (v2,w2), (v3,w3) we select two inequalities i and j that
       * maximize max{vi,vj} + max{wi,wj} this measurement guarantees that select inequalities that
       * separate both important corner points
       */
      getIneqViol(x, y, ineqs[0], ineqs[1], ineqs[2], &viols1[0], &viols2[0]);
      getIneqViol(x, y, ineqs[3], ineqs[4], ineqs[5], &viols1[1], &viols2[1]);
      bestviol = MAX(viols1[0], viols1[1]) + MAX(viols2[0], viols2[1]);

      for( i = 0; i < 2; ++i )
      {
         SCIP_Real viol = MAX(viol1, viols1[i]) + MAX(viol2, viols2[i]);
         if( SCIPisGT(scip, viol, bestviol) )
         {
            bestviol = viol;
            /* remember inequality that should be replaced */
            pos = 1 - i;
         }
      }

      /* replace inequality at pos when replacing an existing inequality improved the total violation */
      if( pos != -1 )
      {
         assert(pos >= 0 && pos < 2);
         ineqs[3*pos] = xcoef;
         ineqs[3*pos+1] = ycoef;
         ineqs[3*pos+2] = constant;
         *success = TRUE;
      }
   }
   SCIPdebugMsg(scip, "accepted inequality? %u\n", *success);

   return SCIP_OKAY;
}


/** creates a SCIP_ROWPREP datastructure
 *
 * Initial cut represents 0 <= 0.
 */
SCIP_RETCODE SCIPcreateRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep,            /**< buffer to store pointer to rowprep */
   SCIP_SIDETYPE         sidetype,           /**< whether cut will be or lower-equal or larger-equal type */
   SCIP_Bool             local               /**< whether cut will be valid only locally */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, rowprep) );
   BMSclearMemory(*rowprep);

   (*rowprep)->sidetype = sidetype;
   (*rowprep)->local = local;

   return SCIP_OKAY;
}

/** frees a SCIP_ROWPREP datastructure */
void SCIPfreeRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep             /**< pointer that stores pointer to rowprep */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(*rowprep != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*rowprep)->vars, (*rowprep)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*rowprep)->coefs, (*rowprep)->varssize);
   SCIPfreeBlockMemory(scip, rowprep);
}

/** creates a copy of a SCIP_ROWPREP datastructure */
SCIP_RETCODE SCIPcopyRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        target,             /**< buffer to store pointer of rowprep copy */
   SCIP_ROWPREP*         source              /**< rowprep to copy */
   )
{
   assert(scip != NULL);
   assert(target != NULL);
   assert(source != NULL);

   SCIP_CALL( SCIPduplicateBlockMemory(scip, target, source) );
   if( source->coefs != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*target)->coefs, source->coefs, source->varssize) );
   }
   if( source->vars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*target)->vars, source->vars, source->varssize) );
   }

   return SCIP_OKAY;
}

/** ensures that rowprep has space for at least given number of additional terms
 *
 * Useful when knowing in advance how many terms will be added.
 */
SCIP_RETCODE SCIPensureRowprepSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   size                /**< number of additional terms for which to alloc space in rowprep */
   )
{
   int oldsize;

   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(size >= 0);

   if( rowprep->varssize >= rowprep->nvars + size )
      return SCIP_OKAY;  /* already enough space left */

   /* realloc vars and coefs array */
   oldsize = rowprep->varssize;
   rowprep->varssize = SCIPcalcMemGrowSize(scip, rowprep->nvars + size);

   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &rowprep->vars,  oldsize, rowprep->varssize) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &rowprep->coefs, oldsize, rowprep->varssize) );

   return SCIP_OKAY;
}

/** prints a rowprep */
void SCIPprintRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be printed */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   )
{
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);

   if( *rowprep->name != '\0' )
   {
      SCIPinfoMessage(scip, file, "[%s](%c) ", rowprep->name, rowprep->local ? 'l' : 'g');
   }

   for( i = 0; i < rowprep->nvars; ++i )
   {
      SCIPinfoMessage(scip, file, "%+g*<%s> ", rowprep->coefs[i], SCIPvarGetName(rowprep->vars[i]));
   }

   SCIPinfoMessage(scip, file, rowprep->sidetype == SCIP_SIDETYPE_LEFT ? ">= %g\n" : "<= %g\n", rowprep->side);
}

/** adds a term coef*var to a rowprep */
SCIP_RETCODE SCIPaddRowprepTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             coef                /**< coefficient to add */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(var != NULL);

   if( coef == 0.0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, 1) );
   assert(rowprep->varssize > rowprep->nvars);

   rowprep->vars[rowprep->nvars] = var;
   rowprep->coefs[rowprep->nvars] = coef;
   ++rowprep->nvars;

   return SCIP_OKAY;
}

/** adds several terms coef*var to a rowprep */
SCIP_RETCODE SCIPaddRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   int                   nvars,              /**< number of terms to add */
   SCIP_VAR**            vars,               /**< variables to add */
   SCIP_Real*            coefs               /**< coefficients to add */
   )
{
   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(vars != NULL || nvars == 0);
   assert(coefs != NULL || nvars == 0);

   if( nvars == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nvars) );
   assert(rowprep->varssize >= rowprep->nvars + nvars);

   /*lint --e{866} */
   BMScopyMemoryArray(rowprep->vars + rowprep->nvars, vars, nvars);
   BMScopyMemoryArray(rowprep->coefs + rowprep->nvars, coefs, nvars);
   rowprep->nvars += nvars;

   return SCIP_OKAY;
}

#ifdef NDEBUG
#undef SCIPaddRowprepSide
#undef SCIPaddRowprepConstant
#endif

/** adds constant value to side of rowprep */
void SCIPaddRowprepSide(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             side                /**< constant value to be added to side */
   )
{
   assert(rowprep != NULL);

   rowprep->side += side;
}

/** adds constant term to rowprep
 *
 * Substracts constant from side.
 */
void SCIPaddRowprepConstant(
   SCIP_ROWPREP*         rowprep,            /**< rowprep */
   SCIP_Real             constant            /**< constant value to be added */
   )
{
   assert(rowprep != NULL);

   SCIPaddRowprepSide(rowprep, -constant);
}

/** computes violation of cut in a given solution */
SCIP_Real SCIPgetRowprepViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_SOL*             sol                 /**< solution or NULL for LP solution */
   )
{
   SCIP_Real activity;
   int i;

   activity = -rowprep->side;
   for( i = 0; i < rowprep->nvars; ++i )
   {
      /* Loose variable have the best bound as LP solution value.
       * HOWEVER, they become column variables when they are added to a row (via SCIPaddVarsToRow below).
       * When this happens, their LP solution value changes to 0.0!
       * So when calculating the row activity for an LP solution, we treat loose variable as if they were already column variables.
       */
      if( sol != NULL || SCIPvarGetStatus(rowprep->vars[i]) != SCIP_VARSTATUS_LOOSE )
         activity += rowprep->coefs[i] * SCIPgetSolVal(scip, sol, rowprep->vars[i]);
   }

   if( rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
      /* cut is activity <= 0.0 -> violation is activity, if positive */
      return MAX(activity, 0.0);
   else
      /* cut is activity >= 0.0 -> violation is -activity, if positive */
      return MAX(-activity, 0.0);
}

/** Merge terms that use same variable and eliminate zero coefficients.
 *
 * Terms are sorted by variable (@see SCIPvarComp) after return.
 */
void SCIPmergeRowprepTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep             /**< rowprep to be cleaned up */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(rowprep != NULL);

   if( rowprep->nvars <= 1 )
      return;

   /* sort terms by variable index */
   SCIPsortPtrReal((void**)rowprep->vars, rowprep->coefs, SCIPvarComp, rowprep->nvars);

   /* merge terms with same variable, drop 0 coefficients */
   i = 0;
   j = 1;
   while( j < rowprep->nvars )
   {
      if( rowprep->vars[i] == rowprep->vars[j] )
      {
         /* merge term j into term i */
         rowprep->coefs[i] += rowprep->coefs[j];
         ++j;
         continue;
      }

      if( rowprep->coefs[i] == 0.0 )
      {
         /* move term j to position i */
         rowprep->coefs[i] = rowprep->coefs[j];
         rowprep->vars[i] = rowprep->vars[j];
         ++j;
         continue;
      }

      /* move term j to position i+1 and move on */
      if( j != i+1 )
      {
         rowprep->vars[i+1] = rowprep->vars[j];
         rowprep->coefs[i+1] = rowprep->coefs[j];
      }
      ++i;
      ++j;
   }

   /* remaining term can have coef zero -> forget about it */
   if( rowprep->coefs[i] == 0.0 )
      --i;

   /* i points to last term */
   rowprep->nvars = i+1;
}

/** sort cut terms by absolute value of coefficients, from largest to smallest */
static
SCIP_RETCODE rowprepCleanupSortTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep             /**< rowprep to be sorted */
   )
{
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);

   /* special treatment for cuts with few variables */
   switch( rowprep->nvars )
   {
      case 0:
      case 1:
         break;

      case 2:
      {
         if( REALABS(rowprep->coefs[0]) < REALABS(rowprep->coefs[1]) )
         {
            SCIP_Real tmp1;
            SCIP_VAR* tmp2;

            tmp1 = rowprep->coefs[0];
            rowprep->coefs[0] = rowprep->coefs[1];
            rowprep->coefs[1] = tmp1;

            tmp2 = rowprep->vars[0];
            rowprep->vars[0] = rowprep->vars[1];
            rowprep->vars[1] = tmp2;
         }
         break;
      }

      default :
      {
         SCIP_Real* abscoefs;

         SCIP_CALL( SCIPallocBufferArray(scip, &abscoefs, rowprep->nvars) );
         for( i = 0; i < rowprep->nvars; ++i )
            abscoefs[i] = REALABS(rowprep->coefs[i]);
         SCIPsortDownRealRealPtr(abscoefs, rowprep->coefs, (void**)rowprep->vars, rowprep->nvars);
         SCIPfreeBufferArray(scip, &abscoefs);
      }
   }

   /* forget about coefs that are exactly zero (unlikely to have some) */
   while( rowprep->nvars > 0 && rowprep->coefs[rowprep->nvars-1] == 0.0 )
      --rowprep->nvars;

   return SCIP_OKAY;
}

/** try to improve coef range by aggregating cut with variable bounds
 *
 * Assumes terms have been sorted by rowprepCleanupSortTerms().
 */
static
void rowprepCleanupImproveCoefrange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             maxcoefrange        /**< maximal allowed coefficients range */
   )
{
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real ref;
   SCIP_Real coef;
   SCIP_Real mincoef;
   SCIP_Real maxcoef;
   SCIP_Real loss[2];
   int maxcoefidx;
   int pos;

   maxcoefidx = 0;
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }
   else
      mincoef = maxcoef = 1.0;

   /* eliminate minimal or maximal coefs as long as coef range is too large
    * this is likely going to eliminate coefs that are within eps of 0.0
    * if not, then we do so after scaling (or should we enforce this here?)
    */
   while( maxcoef / mincoef > maxcoefrange  )
   {
      SCIPdebugMsg(scip, "cut coefficients have very large range: mincoef = %g maxcoef = %g\n", mincoef, maxcoef);

      /* max/min can only be > 1 if there is more than one var
       * we need this below for updating the max/min coef after eliminating a term
       */
      assert(rowprep->nvars > 1);

      /* try to reduce coef range by aggregating with variable bounds
       * that is, eliminate a term like a*x from a*x + ... <= side by adding -a*x <= -a*lb(x)
       * with ref(x) the reference point we try to eliminate, this would weaken the cut by a*(lb(x)-ref(x))
       *
       * we consider eliminating either the term with maximal or the one with minimal coefficient,
       * taking the one that leads to the least weakening of the cut
       *
       * TODO (suggested by @bzfserra, see !496):
       * - Also one could think of not completely removing the coefficient but do an aggregation that makes the coefficient look better. For instance:
       *   say you have $`a x + 0.1 y \leq r`$ and $`y`$ has only an upper bound, $`y \leq b`$,
       *   then you can't really remove $`y`$. However, you could aggregate it with $`0.9 \cdot (y \leq b)`$ to get
       *   $`a x + y \leq r + 0.9 b`$, which has better numerics (and hopefully still cuts the point... actually, if for the point you want to separate, $`y^* = b`$, then the violation is the same)
       */

      for( pos = 0; pos < 2; ++pos )
      {
         var = rowprep->vars[pos ? rowprep->nvars-1 : maxcoefidx];
         coef = rowprep->coefs[pos ? rowprep->nvars-1 : maxcoefidx];
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         ref = SCIPgetSolVal(scip, sol, var);
         assert(coef != 0.0);

         /* make sure reference point is something reasonable within the bounds, preferable the value from the solution */
         if( SCIPisInfinity(scip, REALABS(ref)) )
            ref = 0.0;
         ref = MAX(lb, MIN(ub, ref));

         /* check whether we can eliminate coef*var from rowprep and how much we would loose w.r.t. ref(x) */
         if( ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) )
         {
            /* we would need to aggregate with -coef*var <= -coef*lb(x) */
            if( SCIPisInfinity(scip, -lb) )
               loss[pos] = SCIP_INVALID;
            else
               loss[pos] = REALABS(coef) * (ref - lb);
         }
         else
         {
            assert((coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT));
            /* we would need to aggregate with -coef*var >= -coef*ub(x) */
            if( SCIPisInfinity(scip, ub) )
               loss[pos] = SCIP_INVALID;
            else
               loss[pos] = REALABS(coef) * (ub - ref);
         }
         assert(loss[pos] >= 0.0);  /* assuming SCIP_INVALID >= 0 */

         SCIPdebugMsg(scip, "aggregating %g*<%s> %c= ... with <%s>[%g] %c= %g looses %g\n",
            coef, SCIPvarGetName(var), rowprep->sidetype == SCIP_SIDETYPE_RIGHT ? '<' : '>',
            SCIPvarGetName(var), ref,
            ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) ? '>' : '<',
            ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) ? lb : ub, loss[pos]);
      }

      /*lint --e{777} */
      if( loss[0] == SCIP_INVALID && loss[1] == SCIP_INVALID )
         break;  /* cannot eliminate coefficient */

      /* select position with smaller loss */
      pos = (loss[1] == SCIP_INVALID || loss[1] > loss[0]) ? 0 : 1;

      /* now do the actual elimination */
      var = rowprep->vars[pos ? rowprep->nvars-1 : maxcoefidx];
      coef = rowprep->coefs[pos ? rowprep->nvars-1 : maxcoefidx];

      /* eliminate coef*var from rowprep: increase side */
      if( ((coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT)) )
      {
         /* we aggregate with -coef*var <= -coef*lb(x) */
         assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
         SCIPaddRowprepConstant(rowprep, coef * SCIPvarGetLbLocal(var));
         rowprep->local |= SCIPisGT(scip, SCIPvarGetLbLocal(var), SCIPvarGetLbGlobal(var));
      }
      else
      {
         assert((coef < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT) || (coef > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT));
         /* we aggregate with -coef*var >= -coef*ub(x) */
         assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
         SCIPaddRowprepConstant(rowprep, coef * SCIPvarGetUbLocal(var));
         rowprep->local |= SCIPisLT(scip, SCIPvarGetUbLocal(var), SCIPvarGetUbGlobal(var));
      }

      /* eliminate coef*var from rowprep: remove coef */
      if( pos == 0 )
      {
         /* set first term to zero */
         rowprep->coefs[maxcoefidx] = 0.0;

         /* update index */
         ++maxcoefidx;

         /* update maxcoef */
         maxcoef = REALABS(rowprep->coefs[maxcoefidx]);
      }
      else
      {
         /* forget last term */
         --rowprep->nvars;

         /* update mincoef */
         mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
      }
   }

   /* if maximal coefs were removed, then there are now 0's in the beginning of the coefs array
    * -> move all remaining coefs and vars up front
    */
   if( maxcoefidx > 0 )
   {
      int i;
      for( i = maxcoefidx; i < rowprep->nvars; ++i )
      {
         rowprep->vars[i-maxcoefidx] = rowprep->vars[i];
         rowprep->coefs[i-maxcoefidx] = rowprep->coefs[i];
      }
      rowprep->nvars -= maxcoefidx;
   }
}


/** scales up rowprep if it seems useful */
static
void rowprepCleanupScaleup(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol,               /**< violation of cut in sol (input and output) */
   SCIP_Real             minviol             /**< minimal violation we try to achieve */
   )
{
   SCIP_Real scalefactor;
   SCIP_Real mincoef;
   SCIP_Real maxcoef;

   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(viol != NULL);

   /* if violation is very small than better don't scale up */
   if( *viol < ROWPREP_SCALEUP_VIOLNONZERO )
      return;

   /* if violation is already above minviol, then nothing to do */
   if( *viol >= minviol )
      return;

   /* if violation is sufficiently positive (>10*eps), but has not reached minviol,
    * then consider scaling up to reach approx MINVIOLFACTOR*minviol
    */
   scalefactor = ROWPREP_SCALEUP_MINVIOLFACTOR * minviol / *viol;

   /* scale by approx. scalefactor, if minimal coef is not so large yet and maximal coef and rhs don't get huge by doing so (or have been so before) */
   mincoef = rowprep->nvars > 0 ? REALABS(rowprep->coefs[rowprep->nvars-1]) : 1.0;
   maxcoef = rowprep->nvars > 0 ? REALABS(rowprep->coefs[0]) : 1.0;
   if( mincoef < ROWPREP_SCALEUP_MAXMINCOEF && scalefactor * maxcoef < ROWPREP_SCALEUP_MAXMAXCOEF && scalefactor * REALABS(rowprep->side) < ROWPREP_SCALEUP_MAXSIDE )
   {
      int scaleexp;

      /* SCIPinfoMessage(scip, NULL, "scale up by ~%g, viol=%g: ", scalefactor, myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */

      /* SCIPscaleRowprep returns the actually applied scale factor */
      scaleexp = SCIPscaleRowprep(rowprep, scalefactor);
      *viol = ldexp(*viol, scaleexp);

      /* SCIPinfoMessage(scip, NULL, "scaled up by %g, viol=%g: ", ldexp(1.0, scaleexp), myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */
   }
}

/** scales down rowprep if it improves coefs and keeps rowprep violated */
static
void rowprepCleanupScaledown(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol,               /**< violation of cut in sol (input and output) */
   SCIP_Real             minviol             /**< minimal violation we try to keep */
   )
{
   SCIP_Real scalefactor;

   /* if maxcoef < ROWPREP_SCALEDOWN_MINMAXCOEF (or no terms), then don't consider scaling down */
   if( rowprep->nvars == 0 || REALABS(rowprep->coefs[0]) < ROWPREP_SCALEDOWN_MINMAXCOEF )
      return;

   /* consider scaling down so that maxcoef ~ 10 */
   scalefactor = 10.0 / REALABS(rowprep->coefs[0]);

   /* if minimal violation would be lost by scaling down, then increase scalefactor such that minviol is still reached */
   if( *viol > minviol && scalefactor * *viol < minviol )
   {
      assert(minviol > 0.0);  /* since viol >= 0, the if-condition should ensure that minviol > 0 */
      assert(*viol > 0.0);    /* since minviol > 0, the if-condition ensures viol > 0 */
      scalefactor = ROWPREP_SCALEUP_MINVIOLFACTOR * minviol / *viol;
   }

   /* scale by approx. scalefactor if scaling down and minimal coef does not get too small
    * myviol < minviol (-> scalefactor > 1) or mincoef < feastol before scaling is possible, in which case we also don't scale down
    */
   if( scalefactor < 1.0 && scalefactor * REALABS(rowprep->coefs[rowprep->nvars-1]) > ROWPREP_SCALEDOWN_MINCOEF )
   {
      int scaleexp;

      /* SCIPinfoMessage(scip, NULL, "scale down by ~%g, viol=%g: ", scalefactor, myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */

      scaleexp = SCIPscaleRowprep(rowprep, scalefactor);
      *viol = ldexp(*viol, scaleexp);

      /* SCIPinfoMessage(scip, NULL, "scaled down by %g, viol=%g: ", ldexp(1.0, scaleexp), myviol);
         SCIPprintRowprep(scip, rowprep, NULL); */
   }
}

/** rounds almost integral coefs to integrals, thereby trying to relax the cut */
static
void rowprepCleanupIntegralCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol                /**< violation of cut in sol (input), set to SCIP_INVALID if some coef changed */
   )
{
   SCIP_Real coef;
   SCIP_Real roundcoef;
   int i;

   assert(scip != NULL);
   assert(rowprep != NULL);
   assert(viol != NULL);

   /* Coefficients smaller than epsilon are rounded to 0.0 when added to row and
    * coefficients very close to integral values are rounded to integers when added to LP.
    * Both cases can be problematic if variable value is very large (bad numerics).
    * Thus, we anticipate by rounding coef here, but also modify constant so that cut is still valid (if possible),
    * i.e., bound coef[i]*x by round(coef[i])*x + (coef[i]-round(coef[i])) * bound(x).
    * Or in other words, we aggregate with the variable bound.
    *
    * If the required bound of x is not finite, then only round coef (introduces an error).
    * @TODO If only the opposite bound is available, then one could move the coefficient
    *   away from the closest integer so that the SCIP_ROW won't try to round it.
    */
   for( i = 0; i < rowprep->nvars; ++i )
   {
      coef = rowprep->coefs[i];
      roundcoef = SCIPround(scip, coef);
      if( coef != roundcoef && SCIPisEQ(scip, coef, roundcoef) ) /*lint !e777*/
      {
         SCIP_Real xbnd;
         SCIP_VAR* var;

         var = rowprep->vars[i];
         if( rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
            if( rowprep->local )
               xbnd = coef > roundcoef ? SCIPvarGetLbLocal(var)  : SCIPvarGetUbLocal(var);
            else
               xbnd = coef > roundcoef ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
         else
            if( rowprep->local )
               xbnd = coef > roundcoef ? SCIPvarGetUbLocal(var)  : SCIPvarGetLbLocal(var);
            else
               xbnd = coef > roundcoef ? SCIPvarGetUbGlobal(var) : SCIPvarGetLbGlobal(var);

         if( !SCIPisInfinity(scip, REALABS(xbnd)) )
         {
            /* if there is a bound, then relax row side so rounding coef will not introduce an error */
            SCIPdebugMsg(scip, "var <%s> [%g,%g] has almost integral coef %.20g, round coefficient to %g and add constant %g\n",
               SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), coef, roundcoef, (coef-roundcoef) * xbnd);
            SCIPaddRowprepConstant(rowprep, (coef-roundcoef) * xbnd);
         }
         else
         {
            /* if there is no bound, then we make the coef integral, too, even though this will introduce an error
             * however, SCIP_ROW would do this anyway, but doing this here might eliminate some epsilon coefs (so they don't determine mincoef below)
             * and helps to get a more accurate row violation value
             */
            SCIPdebugMsg(scip, "var <%s> [%g,%g] has almost integral coef %.20g, round coefficient to %g without relaxing side (!)\n",
               SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), coef, roundcoef);
         }
         rowprep->coefs[i] = roundcoef;
         *viol = SCIP_INVALID;
      }
   }

   /* forget about coefs that became exactly zero by the above step */
   while( rowprep->nvars > 0 && rowprep->coefs[rowprep->nvars-1] == 0.0 )
      --rowprep->nvars;
}

/** relaxes almost zero side */
static
void rowprepCleanupSide(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be improve */
   SCIP_Real*            viol                /**< violation of cut in sol (input), set to SCIP_INVALID if some coef changed */
   )
{
   /* SCIP_ROW handling will replace a side close to 0 by 0.0, even if that makes the row more restrictive
    * we thus relax the side here so that it will either be 0 now or will not be rounded to 0 later
    */
   if( !SCIPisZero(scip, rowprep->side) )
      return;

   if( rowprep->side > 0.0 && rowprep->sidetype == SCIP_SIDETYPE_RIGHT )
      rowprep->side =  1.1*SCIPepsilon(scip);
   else if( rowprep->side < 0.0 && rowprep->sidetype == SCIP_SIDETYPE_LEFT )
      rowprep->side = -1.1*SCIPepsilon(scip);
   else
      rowprep->side = 0.0;

   *viol = SCIP_INVALID;
}

/* Cleans up and attempts to improve rowprep
 *
 * Drops small or large coefficients if coefrange is too large, if this can be done by relaxing the cut.
 * Scales coefficients and side up to reach minimal violation, if possible.
 * Scaling is omitted if violation is very small (ROWPREP_SCALEUP_VIOLNONZERO) or
 * maximal coefficient would become huge (ROWPREP_SCALEUP_MAXMAXCOEF).
 * Scales coefficients and side down if they are large and if the minimal violation is still reached.
 * Rounds coefficients close to integral values to integrals, if this can be done by relaxing the cut.
 * Rounds side within epsilon of 0 to 0.0 or +/-1.1*epsilon, whichever relaxes the cut least.
 *
 * After return, the terms in the rowprep will be sorted by absolute value of coefficient, in decreasing order.
 */
SCIP_RETCODE SCIPcleanupRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be cleaned */
   SCIP_SOL*             sol,                /**< solution that we try to cut off, or NULL for LP solution */
   SCIP_Real             maxcoefrange,       /**< maximal allowed coefficients range */
   SCIP_Real             minviol,            /**< minimal absolute violation the row should achieve (w.r.t. sol) */
   SCIP_Real*            coefrange,          /**< buffer to store coefrange of cleaned up cut, or NULL if not of interest */
   SCIP_Real*            viol                /**< buffer to store absolute violation of cleaned up cut in sol, or NULL if not of interest */
   )
{
   SCIP_Real myviol;
#ifdef SCIP_DEBUG
   SCIP_Real mincoef = 1.0;
   SCIP_Real maxcoef = 1.0;
#endif

   assert(maxcoefrange > 1.0);   /* not much interesting otherwise */

   /* sort term by absolute value of coef. */
   SCIP_CALL( rowprepCleanupSortTerms(scip, rowprep) );

#ifdef SCIP_DEBUG
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }

   SCIPinfoMessage(scip, NULL, "starting cleanup, coefrange %g: ", maxcoef/mincoef);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* improve coefficient range by aggregating out variables */
   rowprepCleanupImproveCoefrange(scip, rowprep, sol, maxcoefrange);

   /* get current violation in sol */
   myviol = SCIPgetRowprepViolation(scip, rowprep, sol);
   assert(myviol >= 0.0);

#ifdef SCIP_DEBUG
   if( rowprep->nvars > 0 )
   {
      maxcoef = REALABS(rowprep->coefs[0]);
      mincoef = REALABS(rowprep->coefs[rowprep->nvars-1]);
   }

   SCIPinfoMessage(scip, NULL, "improved coefrange to %g, viol %g: ", maxcoef / mincoef, myviol);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* if there is interest in achieving some minimal violation, then possibly scale up to increase violation, updates myviol */
   if( minviol > 0.0 )
   {
      /* first, try to achieve scip's minefficacy (typically 1e-4) */
      if( SCIPgetSepaMinEfficacy(scip) > minviol )
         rowprepCleanupScaleup(scip, rowprep, &myviol, SCIPgetSepaMinEfficacy(scip));
      /* in case scip minefficacy could not be reached or was smaller than minviol, try with the given minviol */
      rowprepCleanupScaleup(scip, rowprep, &myviol, minviol);
   }

   /* scale down to improve numerics, updates myviol */
   rowprepCleanupScaledown(scip, rowprep, &myviol, MAX(SCIPgetSepaMinEfficacy(scip), minviol)); /*lint !e666*/

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "applied scaling, viol %g: ", myviol);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* turn almost-integral coefs to integral values, may set myviol to SCIP_INVALID */
   rowprepCleanupIntegralCoefs(scip, rowprep, &myviol);

   /* relax almost-zero side, may set myviol to SCIP_INVALID */
   rowprepCleanupSide(scip, rowprep, &myviol);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "adjusted almost-integral coefs and sides, viol %g: ", myviol);
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

   /* compute final coefrange, if requested by caller */
   if( coefrange != NULL )
   {
      if( rowprep->nvars > 0 )
         *coefrange = REALABS(rowprep->coefs[0]) / REALABS(rowprep->coefs[rowprep->nvars-1]);
      else
         *coefrange = 1.0;
   }

   /* If we updated myviol correctly, then it should coincide with freshly computed violation.
    * I leave this assert off for now, since getting the tolerance in the EQ correctly is not trivial. We recompute viol below anyway.
    */
   /* assert(myviol == SCIP_INVALID || SCIPisEQ(scip, myviol, SCIPgetRowprepViolation(scip, rowprep, sol))); */

   /* compute final violation, if requested by caller */
   if( viol != NULL )  /*lint --e{777} */
      *viol = myviol == SCIP_INVALID ? SCIPgetRowprepViolation(scip, rowprep, sol) : myviol;

   return SCIP_OKAY;
}

/** scales a rowprep
 *
 * @return Exponent of actually applied scaling factor, if written as 2^x.
 */
int SCIPscaleRowprep(
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be scaled */
   SCIP_Real             factor              /**< suggested scale factor */
   )
{
   double v;
   int expon;
   int i;

   assert(rowprep != NULL);
   assert(factor > 0.0);

   /* write factor as v*2^expon with v in [0.5,1) */
   v = frexp(factor, &expon);
   /* adjust to v'*2^expon with v' in (0.5,1] by v'=v if v > 0.5, v'=1 if v=0.5 */
   if( v == 0.5 )
      --expon;

   /* multiply each coefficient by 2^expon */
   for( i = 0; i < rowprep->nvars; ++i )
      rowprep->coefs[i] = ldexp(rowprep->coefs[i], expon);

   /* multiply side by 2^expon */
   rowprep->side = ldexp(rowprep->side, expon);

   return expon;
}

/** generates a SCIP_ROW from a rowprep */
SCIP_RETCODE SCIPgetRowprepRowCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   assert(scip != NULL);
   assert(row != NULL);
   assert(rowprep != NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, row, conshdlr, rowprep->name,
      rowprep->sidetype == SCIP_SIDETYPE_LEFT  ? rowprep->side : -SCIPinfinity(scip),
      rowprep->sidetype == SCIP_SIDETYPE_RIGHT ? rowprep->side :  SCIPinfinity(scip),
      rowprep->local && (SCIPgetDepth(scip) > 0), FALSE, TRUE) );

   SCIP_CALL( SCIPaddVarsToRow(scip, *row, rowprep->nvars, rowprep->vars, rowprep->coefs) );

   return SCIP_OKAY;
}

/** generates a SCIP_ROW from a rowprep */
SCIP_RETCODE SCIPgetRowprepRowSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store pointer to new row */
   SCIP_ROWPREP*         rowprep,            /**< rowprep to be turned into a row */
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(scip != NULL);
   assert(row != NULL);
   assert(rowprep != NULL);

   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, row, sepa, rowprep->name,
      rowprep->sidetype == SCIP_SIDETYPE_LEFT  ? rowprep->side : -SCIPinfinity(scip),
      rowprep->sidetype == SCIP_SIDETYPE_RIGHT ? rowprep->side :  SCIPinfinity(scip),
      rowprep->local && (SCIPgetDepth(scip) > 0), FALSE, TRUE) );

   SCIP_CALL( SCIPaddVarsToRow(scip, *row, rowprep->nvars, rowprep->vars, rowprep->coefs) );

   return SCIP_OKAY;
}
