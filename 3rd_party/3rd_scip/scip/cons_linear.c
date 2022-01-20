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

/**@file   cons_linear.c
 * @brief Constraint handler for linear constraints in their most general form, \f$lhs <= a^T x <= rhs\f$.
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Michael Winkler
 * @author Gerald Gamrath
 * @author Domenico Salvagnin
 *
 *  Linear constraints are separated with a high priority, because they are easy
 *  to separate. Instead of using the global cut pool, the same effect can be
 *  implemented by adding linear constraints to the root node, such that they are
 *  separated each time, the linear constraints are separated. A constraint
 *  handler, which generates linear constraints in this way should have a lower
 *  separation priority than the linear constraint handler, and it should have a
 *  separation frequency that is a multiple of the frequency of the linear
 *  constraint handler. In this way, it can be avoided to separate the same cut
 *  twice, because if a separation run of the handler is always preceded by a
 *  separation of the linear constraints, the priorily added constraints are
 *  always satisfied.
 *
 *  Linear constraints are enforced and checked with a very low priority. Checking
 *  of (many) linear constraints is much more involved than checking the solution
 *  values for integrality. Because we are separating the linear constraints quite
 *  often, it is only necessary to enforce them for integral solutions. A constraint
 *  handler which generates pool cuts in its enforcing method should have an
 *  enforcing priority smaller than that of the linear constraint handler to avoid
 *  regenerating constraints which already exist.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc.h"
#include "scip/debug.h"

#define CONSHDLR_NAME          "linear"
#define CONSHDLR_DESC          "linear constraints of the form  lhs <= a^T x <= rhs"
#define CONSHDLR_SEPAPRIORITY   +100000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING    (SCIP_PRESOLTIMING_FAST | SCIP_PRESOLTIMING_EXHAUSTIVE) /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "linear"
#define EVENTHDLR_DESC         "bound change event handler for linear constraints"

#define CONFLICTHDLR_NAME      "linear"
#define CONFLICTHDLR_DESC      "conflict handler creating linear constraints"
#define CONFLICTHDLR_PRIORITY  -1000000

#define DEFAULT_TIGHTENBOUNDSFREQ       1 /**< multiplier on propagation frequency, how often the bounds are tightened */
#define DEFAULT_MAXROUNDS               5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT          -1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS            50 /**< maximal number of cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT       200 /**< maximal number of cuts separated per separation round in root node */
#define DEFAULT_PRESOLPAIRWISE       TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define DEFAULT_PRESOLUSEHASHING     TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define DEFAULT_NMINCOMPARISONS    200000 /**< number for minimal pairwise presolving comparisons */
#define DEFAULT_MINGAINPERNMINCOMP  1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise 
                                           *   comparison round */
#define DEFAULT_SORTVARS             TRUE /**< should variables be sorted after presolve w.r.t their coefficient absolute for faster
                                           *  propagation? */
#define DEFAULT_CHECKRELMAXABS      FALSE /**< should the violation for a constraint with side 0.0 be checked relative
                                           *   to 1.0 (FALSE) or to the maximum absolute value in the activity (TRUE)? */
#define DEFAULT_MAXAGGRNORMSCALE      0.0 /**< maximal allowed relative gain in maximum norm for constraint aggregation
                                           *   (0.0: disable constraint aggregation) */
#define DEFAULT_MAXEASYACTIVITYDELTA  1e6 /**< maximum activity delta to run easy propagation on linear constraint
                                           *   (faster, but numerically less stable) */
#define DEFAULT_MAXCARDBOUNDDIST      0.0 /**< maximal relative distance from current node's dual bound to primal bound compared
                                           *   to best node's dual bound for separating knapsack cardinality cuts */
#define DEFAULT_SEPARATEALL         FALSE /**< should all constraints be subject to cardinality cut generation instead of only
                                           *   the ones with non-zero dual value? */
#define DEFAULT_AGGREGATEVARIABLES   TRUE /**< should presolving search for redundant variables in equations */
#define DEFAULT_SIMPLIFYINEQUALITIES TRUE /**< should presolving try to simplify inequalities */
#define DEFAULT_DUALPRESOLVING       TRUE /**< should dual presolving steps be performed? */
#define DEFAULT_SINGLETONSTUFFING    TRUE /**< should stuffing of singleton continuous variables be performed? */
#define DEFAULT_SINGLEVARSTUFFING   FALSE /**< should single variable stuffing be performed, which tries to fulfill
                                           *   constraints using the cheapest variable? */
#define DEFAULT_DETECTCUTOFFBOUND    TRUE /**< should presolving try to detect constraints parallel to the objective
                                           *   function defining an upper bound and prevent these constraints from
                                           *   entering the LP */
#define DEFAULT_DETECTLOWERBOUND     TRUE /**< should presolving try to detect constraints parallel to the objective
                                           *   function defining a lower bound and prevent these constraints from
                                           *   entering the LP */
#define DEFAULT_DETECTPARTIALOBJECTIVE TRUE/**< should presolving try to detect subsets of constraints parallel to the
                                            *   objective function */
#define DEFAULT_RANGEDROWPROPAGATION TRUE /**< should we perform ranged row propagation */
#define DEFAULT_RANGEDROWARTCONS     TRUE /**< should presolving and propagation extract sub-constraints from ranged rows and equations? */
#define DEFAULT_RANGEDROWMAXDEPTH INT_MAX /**< maximum depth to apply ranged row propagation */
#define DEFAULT_RANGEDROWFREQ           1 /**< frequency for applying ranged row propagation */

#define DEFAULT_MULTAGGRREMOVE      FALSE /**< should multi-aggregations only be performed if the constraint can be
                                           *   removed afterwards? */

#define MAXDNOM                   10000LL /**< maximal denominator for simple rational fixed values */
#define MAXSCALEDCOEF               1e+03 /**< maximal coefficient value after scaling */
#define MAXSCALEDCOEFINTEGER        1e+05 /**< maximal coefficient value after scaling if all variables are of integral
                                           *   type
                                           */

#define MAXVALRECOMP                1e+06 /**< maximal abolsute value we trust without recomputing the activity */
#define MINVALRECOMP                1e-05 /**< minimal abolsute value we trust without recomputing the activity */


#define QUADCONSUPGD_PRIORITY     1000000 /**< priority of the constraint handler for upgrading of quadratic constraints */
#define NONLINCONSUPGD_PRIORITY   1000000 /**< priority of the constraint handler for upgrading of nonlinear constraints */

/* @todo add multi-aggregation of variables that are in exactly two equations (, if not numerically an issue),
 *       maybe in fullDualPresolve(), see convertLongEquality()
 */


/** constraint data for linear constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of row (for ranged rows) */
   SCIP_Real             rhs;                /**< right hand side of row */
   SCIP_Real             maxabsval;          /**< maximum absolute value of all coefficients */
   SCIP_Real             minabsval;          /**< minimal absolute value of all coefficients */
   SCIP_Real             minactivity;        /**< minimal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             maxactivity;        /**< maximal value w.r.t. the variable's local bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             lastminactivity;    /**< last minimal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             lastmaxactivity;    /**< last maximal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             glbminactivity;     /**< minimal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             glbmaxactivity;     /**< maximal value w.r.t. the variable's global bounds for the constraint's
                                              *   activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real             lastglbminactivity; /**< last global minimal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             lastglbmaxactivity; /**< last global maximal activity which was computed by complete summation
                                              *   over all contributing values */
   SCIP_Real             maxactdelta;        /**< maximal activity contribution of a single variable, or SCIP_INVALID if invalid */
   SCIP_VAR*             maxactdeltavar;     /**< variable with maximal activity contribution, or NULL if invalid */
   uint64_t              possignature;       /**< bit signature of coefficients that may take a positive value */
   uint64_t              negsignature;       /**< bit signature of coefficients that may take a negative value */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   SCIP_VAR**            vars;               /**< variables of constraint entries */
   SCIP_Real*            vals;               /**< coefficients of constraint entries */
   SCIP_EVENTDATA**      eventdata;          /**< event data for bound change events of the variables */
   int                   minactivityneginf;  /**< number of coefficients contributing with neg. infinite value to minactivity */
   int                   minactivityposinf;  /**< number of coefficients contributing with pos. infinite value to minactivity */
   int                   maxactivityneginf;  /**< number of coefficients contributing with neg. infinite value to maxactivity */
   int                   maxactivityposinf;  /**< number of coefficients contributing with pos. infinite value to maxactivity */
   int                   minactivityneghuge; /**< number of coefficients contributing with huge neg. value to minactivity */
   int                   minactivityposhuge; /**< number of coefficients contributing with huge pos. value to minactivity */
   int                   maxactivityneghuge; /**< number of coefficients contributing with huge neg. value to maxactivity */
   int                   maxactivityposhuge; /**< number of coefficients contributing with huge pos. value to maxactivity */
   int                   glbminactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbminactivity */
   int                   glbminactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbminactivity */
   int                   glbmaxactivityneginf;/**< number of coefficients contrib. with neg. infinite value to glbmaxactivity */
   int                   glbmaxactivityposinf;/**< number of coefficients contrib. with pos. infinite value to glbmaxactivity */
   int                   glbminactivityneghuge;/**< number of coefficients contrib. with huge neg. value to glbminactivity */
   int                   glbminactivityposhuge;/**< number of coefficients contrib. with huge pos. value to glbminactivity */
   int                   glbmaxactivityneghuge;/**< number of coefficients contrib. with huge neg. value to glbmaxactivity */
   int                   glbmaxactivityposhuge;/**< number of coefficients contrib. with huge pos. value to glbmaxactivity */
   int                   varssize;           /**< size of the vars- and vals-arrays */
   int                   nvars;              /**< number of nonzeros in constraint */
   int                   nbinvars;           /**< the number of binary variables in the constraint, only valid after
                                              *   sorting in stage >= SCIP_STAGE_INITSOLVE
                                              */
   unsigned int          boundstightened:2;  /**< is constraint already propagated with bound tightening? */
   unsigned int          rangedrowpropagated:2; /**< did we perform ranged row propagation on this constraint?
                                                 *   (0: no, 1: yes, 2: with potentially adding artificial constraint */
   unsigned int          validmaxabsval:1;   /**< is the maximum absolute value valid? */
   unsigned int          validminabsval:1;   /**< is the minimum absolute value valid? */
   unsigned int          validactivities:1;  /**< are the activity bounds (local and global) valid? */
   unsigned int          validminact:1;      /**< is the local minactivity valid? */
   unsigned int          validmaxact:1;      /**< is the local maxactivity valid? */
   unsigned int          validglbminact:1;   /**< is the global minactivity valid? */
   unsigned int          validglbmaxact:1;   /**< is the global maxactivity valid? */
   unsigned int          presolved:1;        /**< is constraint already presolved? */
   unsigned int          removedfixings:1;   /**< are all fixed variables removed from the constraint? */
   unsigned int          validsignature:1;   /**< is the bit signature valid? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          normalized:1;       /**< is the constraint in normalized form? */
   unsigned int          upgradetried:1;     /**< was the constraint already tried to be upgraded? */
   unsigned int          upgraded:1;         /**< is the constraint upgraded and will it be removed after preprocessing? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          merged:1;           /**< are the constraint's equal variables already merged? */
   unsigned int          cliquesadded:1;     /**< were the cliques of the constraint already extracted? */
   unsigned int          implsadded:1;       /**< were the implications of the constraint already extracted? */
   unsigned int          binvarssorted:1;    /**< are binary variables sorted w.r.t. the absolute of their coefficient? */
   unsigned int          varsdeleted:1;      /**< were variables deleted after last cleanup? */
   unsigned int          hascontvar:1;       /**< does the constraint contain at least one continuous variable? */
   unsigned int          hasnonbinvar:1;     /**< does the constraint contain at least one non-binary variable? */
   unsigned int          hasnonbinvalid:1;   /**< is the information stored in hasnonbinvar and hascontvar valid? */
   unsigned int          checkabsolute:1;    /**< should the constraint be checked w.r.t. an absolute feasibilty tolerance? */
};

/** event data for bound change event */
struct SCIP_EventData
{
   SCIP_CONS*            cons;               /**< linear constraint to process the bound change for */
   int                   varpos;             /**< position of variable in vars array */
   int                   filterpos;          /**< position of event in variable's event filter */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_LINCONSUPGRADE** linconsupgrades;    /**< linear constraint upgrade methods for specializing linear constraints */
   SCIP_Real             maxaggrnormscale;   /**< maximal allowed relative gain in maximum norm for constraint aggregation
                                              *   (0.0: disable constraint aggregation) */
   SCIP_Real             maxcardbounddist;   /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for separating knapsack cardinality cuts */
   SCIP_Real             mingainpernmincomp; /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */
   SCIP_Real             maxeasyactivitydelta;/**< maximum activity delta to run easy propagation on linear constraint
                                               *   (faster, but numerically less stable) */
   int                   linconsupgradessize;/**< size of linconsupgrade array */
   int                   nlinconsupgrades;   /**< number of linear constraint upgrade methods */
   int                   tightenboundsfreq;  /**< multiplier on propagation frequency, how often the bounds are tightened */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in root node */
   int                   nmincomparisons;    /**< number for minimal pairwise presolving comparisons */
   int                   naddconss;          /**< number of added constraints */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
   SCIP_Bool             separateall;        /**< should all constraints be subject to cardinality cut generation instead of only
                                              *   the ones with non-zero dual value? */
   SCIP_Bool             aggregatevariables; /**< should presolving search for redundant variables in equations */
   SCIP_Bool             simplifyinequalities;/**< should presolving try to cancel down or delete coefficients in inequalities */
   SCIP_Bool             dualpresolving;     /**< should dual presolving steps be performed? */
   SCIP_Bool             singletonstuffing;  /**< should stuffing of singleton continuous variables be performed? */
   SCIP_Bool             singlevarstuffing;  /**< should single variable stuffing be performed, which tries to fulfill
                                              *   constraints using the cheapest variable? */
   SCIP_Bool             sortvars;           /**< should binary variables be sorted for faster propagation? */
   SCIP_Bool             checkrelmaxabs;     /**< should the violation for a constraint with side 0.0 be checked relative
                                              *   to 1.0 (FALSE) or to the maximum absolute value in the activity (TRUE)? */
   SCIP_Bool             detectcutoffbound;  /**< should presolving try to detect constraints parallel to the objective
                                              *   function defining an upper bound and prevent these constraints from
                                              *   entering the LP */
   SCIP_Bool             detectlowerbound;   /**< should presolving try to detect constraints parallel to the objective
                                              *   function defining a lower bound and prevent these constraints from
                                              *   entering the LP */
   SCIP_Bool             detectpartialobjective;/**< should presolving try to detect subsets of constraints parallel to
                                                 *   the objective function */
   SCIP_Bool             rangedrowpropagation;/**< should presolving and propagation try to improve bounds, detect
                                               *   infeasibility, and extract sub-constraints from ranged rows and
                                               *   equations */
   SCIP_Bool             rangedrowartcons;   /**< should presolving and propagation extract sub-constraints from ranged rows and equations?*/
   int                   rangedrowmaxdepth;  /**< maximum depth to apply ranged row propagation */
   int                   rangedrowfreq;      /**< frequency for applying ranged row propagation */
   SCIP_Bool             multaggrremove;     /**< should multi-aggregations only be performed if the constraint can be
                                              *   removed afterwards? */
};

/** linear constraint update method */
struct SCIP_LinConsUpgrade
{
   SCIP_DECL_LINCONSUPGD((*linconsupgd));    /**< method to call for upgrading linear constraint */
   int                   priority;           /**< priority of upgrading method */
   SCIP_Bool             active;             /**< is upgrading enabled */
};


/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1_RHS        = 1,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the right hand side of the inequality */
   PROPRULE_1_LHS        = 2,                /**< activity residuals of all other variables tighten bounds of single
                                              *   variable due to the left hand side of the inequality */
   PROPRULE_1_RANGEDROW  = 3,                /**< fixed variables and gcd of all left variables tighten bounds of a
					      *   single variable in this reanged row */
   PROPRULE_INVALID      = 0                 /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;

/** inference information */
struct InferInfo
{
   union
   {
      struct
      {
         unsigned int    proprule:8;         /**< propagation rule that was applied */
         unsigned int    pos:24;             /**< variable position, the propagation rule was applied at */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};
typedef struct InferInfo INFERINFO;

/** converts an integer into an inference information */
static
INFERINFO intToInferInfo(
   int                   i                   /**< integer to convert */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asint = i;

   return inferinfo;
}

/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}

/** returns the propagation rule stored in the inference information */
static
int inferInfoGetProprule(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (int) inferinfo.val.asbits.proprule;
}

/** returns the position stored in the inference information */
static
int inferInfoGetPos(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (int) inferinfo.val.asbits.pos;
}

/** constructs an inference information out of a propagation rule and a position number */
static
INFERINFO getInferInfo(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   pos                 /**< variable position, the propagation rule was applied at */
   )
{
   INFERINFO inferinfo;

   assert(pos >= 0);
   /* in the inferinfo struct only 24 bits for 'pos' are reserved */
   assert(pos < (1<<24));

   inferinfo.val.asbits.proprule = (unsigned int) proprule; /*lint !e641*/
   inferinfo.val.asbits.pos = (unsigned int) pos; /*lint !e732*/

   return inferinfo;
}

/** constructs an inference information out of a propagation rule and a position number, returns info as int */
static
int getInferInt(
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   int                   pos                 /**< variable position, the propagation rule was applied at */
   )
{
   return inferInfoToInt(getInferInfo(proprule, pos));
}


/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that linconsupgrades array can store at least num entries */
static
SCIP_RETCODE conshdlrdataEnsureLinconsupgradesSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< linear constraint handler data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->nlinconsupgrades <= conshdlrdata->linconsupgradessize);

   if( num > conshdlrdata->linconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->linconsupgrades, conshdlrdata->linconsupgradessize, newsize) );
      conshdlrdata->linconsupgradessize = newsize;
   }
   assert(num <= conshdlrdata->linconsupgradessize);

   return SCIP_OKAY;
}

/** ensures, that vars and vals arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);

   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vals, consdata->varssize, newsize) );
      if( consdata->eventdata != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->eventdata, consdata->varssize, newsize) );
      }
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}


/*
 * local methods for managing linear constraint update methods
 */

/** creates a linear constraint upgrade data object */
static
SCIP_RETCODE linconsupgradeCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSUPGRADE** linconsupgrade,     /**< pointer to store the linear constraint upgrade */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int                   priority            /**< priority of upgrading method */
   )
{
   assert(scip != NULL);
   assert(linconsupgrade != NULL);
   assert(linconsupgd != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, linconsupgrade) );
   (*linconsupgrade)->linconsupgd = linconsupgd;
   (*linconsupgrade)->priority = priority;
   (*linconsupgrade)->active = TRUE;

   return SCIP_OKAY;
}

/** frees a linear constraint upgrade data object */
static
void linconsupgradeFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSUPGRADE** linconsupgrade      /**< pointer to the linear constraint upgrade */
   )
{
   assert(scip != NULL);
   assert(linconsupgrade != NULL);
   assert(*linconsupgrade != NULL);

   SCIPfreeBlockMemory(scip, linconsupgrade);
}

/** creates constraint handler data for linear constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata,       /**< pointer to store the constraint handler data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, conshdlrdata) );
   (*conshdlrdata)->linconsupgrades = NULL;
   (*conshdlrdata)->linconsupgradessize = 0;
   (*conshdlrdata)->nlinconsupgrades = 0;
   (*conshdlrdata)->naddconss = 0;

   /* set event handler for updating linear constraint activity bounds */
   (*conshdlrdata)->eventhdlr = eventhdlr;

   return SCIP_OKAY;
}

/** frees constraint handler data for linear constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   for( i = 0; i < (*conshdlrdata)->nlinconsupgrades; ++i )
   {
      linconsupgradeFree(scip, &(*conshdlrdata)->linconsupgrades[i]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*conshdlrdata)->linconsupgrades, (*conshdlrdata)->linconsupgradessize);

   SCIPfreeBlockMemory(scip, conshdlrdata);
}

/** creates a linear constraint upgrade data object */
static
SCIP_Bool conshdlrdataHasUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(linconsupgd != NULL);
   assert(conshdlrname != NULL);

   for( i = conshdlrdata->nlinconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->linconsupgrades[i]->linconsupgd == linconsupgd )
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage(scip, "Try to add already known upgrade message %p for constraint handler %s.\n", linconsupgd, conshdlrname);
#endif
         return TRUE;
      }
   }

   return FALSE;
}

/** adds a linear constraint update method to the constraint handler's data */
static
SCIP_RETCODE conshdlrdataIncludeUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_LINCONSUPGRADE*  linconsupgrade      /**< linear constraint upgrade method */
   )
{
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(linconsupgrade != NULL);

   SCIP_CALL( conshdlrdataEnsureLinconsupgradesSize(scip, conshdlrdata, conshdlrdata->nlinconsupgrades+1) );

   for( i = conshdlrdata->nlinconsupgrades;
        i > 0 && conshdlrdata->linconsupgrades[i-1]->priority < linconsupgrade->priority; --i )
   {
      conshdlrdata->linconsupgrades[i] = conshdlrdata->linconsupgrades[i-1];
   }
   assert(0 <= i && i <= conshdlrdata->nlinconsupgrades);
   conshdlrdata->linconsupgrades[i] = linconsupgrade;
   conshdlrdata->nlinconsupgrades++;

   return SCIP_OKAY;
}

/*
 * local methods
 */

/** installs rounding locks for the given variable associated to the given coefficient in the linear constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable associated to the given coefficient in the linear constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, val));

   if( SCIPisPositive(scip, val) )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons,
            !SCIPisInfinity(scip, consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** creates event data for variable at given position, and catches events */
/**! [SnippetDebugAssertions] */
static
SCIP_RETCODE consCatchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars != NULL);
   assert(consdata->vars[pos] != NULL);
   assert(SCIPvarIsTransformed(consdata->vars[pos]));
   assert(consdata->eventdata != NULL);
   assert(consdata->eventdata[pos] == NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &(consdata->eventdata[pos])) ); /*lint !e866*/
   consdata->eventdata[pos]->cons = cons;
   consdata->eventdata[pos]->varpos = pos;

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED
         | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_VARDELETED,
         eventhdlr, consdata->eventdata[pos], &consdata->eventdata[pos]->filterpos) );

   return SCIP_OKAY;
}
/**! [SnippetDebugAssertions] */

/** deletes event data for variable at given position, and drops events */
static
SCIP_RETCODE consDropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< array position of variable to catch bound change events for */
   )
{
   SCIP_CONSDATA* consdata;
   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(0 <= pos && pos < consdata->nvars);
   assert(consdata->vars[pos] != NULL);
   assert(consdata->eventdata != NULL);
   assert(consdata->eventdata[pos] != NULL);
   assert(consdata->eventdata[pos]->cons == cons);
   assert(consdata->eventdata[pos]->varpos == pos);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[pos],
         SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED
         | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_VARDELETED,
         eventhdlr, consdata->eventdata[pos], consdata->eventdata[pos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->eventdata[pos]); /*lint !e866*/

   return SCIP_OKAY;
}

/** catches bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consCatchAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->eventdata == NULL);

   /* allocate eventdata array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->eventdata, consdata->varssize) );
   assert(consdata->eventdata != NULL);
   BMSclearMemoryArray(consdata->eventdata, consdata->nvars);

   /* catch event for every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( consCatchEvent(scip, cons, eventhdlr, i) );
   }

   return SCIP_OKAY;
}

/** drops bound change events for all variables in transformed linear constraint */
static
SCIP_RETCODE consDropAllEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->eventdata != NULL);

   /* drop event of every single variable */
   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      SCIP_CALL( consDropEvent(scip, cons, eventhdlr, i) );
   }

   /* free eventdata array */
   SCIPfreeBlockMemoryArray(scip, &consdata->eventdata, consdata->varssize);
   assert(consdata->eventdata == NULL);

   return SCIP_OKAY;
}

/** creates a linear constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to linear constraint data */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs                 /**< right hand side of row */
   )
{
   int v;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, lhs) )
      lhs = SCIPinfinity(scip);

   if( SCIPisGT(scip, lhs, rhs) )
   {
      SCIPwarningMessage(scip, "left hand side of linear constraint greater than right hand side\n");
      SCIPwarningMessage(scip, " -> lhs=%g, rhs=%g\n", lhs, rhs);
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->varssize = 0;
   (*consdata)->nvars = nvars;
   (*consdata)->hascontvar = FALSE;
   (*consdata)->hasnonbinvar = FALSE;
   (*consdata)->hasnonbinvalid = TRUE;
   (*consdata)->vars = NULL;
   (*consdata)->vals = NULL;

   constant = 0.0;
   if( nvars > 0 )
   {
      int k;

      SCIP_VAR** varsbuffer;
      SCIP_Real* valsbuffer;

      /* copy variables into temporary buffer */
      SCIP_CALL( SCIPallocBufferArray(scip, &varsbuffer, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &valsbuffer, nvars) );
      k = 0;

      /* loop over variables and sort out fixed ones */
      for( v = 0; v < nvars; ++v )
      {
         SCIP_VAR* var;
         SCIP_Real val;

         var = vars[v];
         val = vals[v];

         assert(var != NULL);
         if( !SCIPisZero(scip, val) )
         {
            /* treat fixed variable as a constant if problem compression is enabled */
            if( SCIPisConsCompressionEnabled(scip) && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
            {
               constant += SCIPvarGetLbGlobal(var) * val;
            }
            else
            {
               varsbuffer[k] = var;
               valsbuffer[k] = val;
               k++;

               /* update hascontvar and hasnonbinvar flags */
               if( !(*consdata)->hascontvar )
               {
                  SCIP_VARTYPE vartype = SCIPvarGetType(var);

                  if( vartype != SCIP_VARTYPE_BINARY )
                  {
                     (*consdata)->hasnonbinvar = TRUE;

                     if( vartype == SCIP_VARTYPE_CONTINUOUS )
                        (*consdata)->hascontvar = TRUE;
                  }
               }
            }
         }
      }
      (*consdata)->nvars = k;

      if( k > 0 )
      {
         /* copy the possibly reduced buffer arrays into block */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, varsbuffer, k) );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vals, valsbuffer, k) );
         (*consdata)->varssize = k;
      }
      /* free temporary buffer */
      SCIPfreeBufferArray(scip, &valsbuffer);
      SCIPfreeBufferArray(scip, &varsbuffer);
   }

   (*consdata)->eventdata = NULL;

   /* due to compressed copying, we may have fixed variables contributing to the left and right hand side */
   if( !SCIPisZero(scip, constant) )
   {
      if( !SCIPisInfinity(scip, REALABS(lhs)) )
         lhs -= constant;

      if( !SCIPisInfinity(scip, REALABS(rhs)) )
         rhs -= constant;
   }

   (*consdata)->row = NULL;
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   (*consdata)->maxabsval = SCIP_INVALID;
   (*consdata)->minabsval = SCIP_INVALID;
   (*consdata)->minactivity = SCIP_INVALID;
   (*consdata)->maxactivity = SCIP_INVALID;
   (*consdata)->lastminactivity = SCIP_INVALID;
   (*consdata)->lastmaxactivity = SCIP_INVALID;
   (*consdata)->maxactdelta = SCIP_INVALID;
   (*consdata)->maxactdeltavar = NULL;
   (*consdata)->minactivityneginf = -1;
   (*consdata)->minactivityposinf = -1;
   (*consdata)->maxactivityneginf = -1;
   (*consdata)->maxactivityposinf = -1;
   (*consdata)->minactivityneghuge = -1;
   (*consdata)->minactivityposhuge = -1;
   (*consdata)->maxactivityneghuge = -1;
   (*consdata)->maxactivityposhuge = -1;
   (*consdata)->glbminactivity = SCIP_INVALID;
   (*consdata)->glbmaxactivity = SCIP_INVALID;
   (*consdata)->lastglbminactivity = SCIP_INVALID;
   (*consdata)->lastglbmaxactivity = SCIP_INVALID;
   (*consdata)->glbminactivityneginf = -1;
   (*consdata)->glbminactivityposinf = -1;
   (*consdata)->glbmaxactivityneginf = -1;
   (*consdata)->glbmaxactivityposinf = -1;
   (*consdata)->glbminactivityneghuge = -1;
   (*consdata)->glbminactivityposhuge = -1;
   (*consdata)->glbmaxactivityneghuge = -1;
   (*consdata)->glbmaxactivityposhuge = -1;
   (*consdata)->possignature = 0;
   (*consdata)->negsignature = 0;
   (*consdata)->validmaxabsval = FALSE;
   (*consdata)->validminabsval = FALSE;
   (*consdata)->validactivities = FALSE;
   (*consdata)->validminact = FALSE;
   (*consdata)->validmaxact = FALSE;
   (*consdata)->validglbminact = FALSE;
   (*consdata)->validglbmaxact = FALSE;
   (*consdata)->boundstightened = 0;
   (*consdata)->presolved = FALSE;
   (*consdata)->removedfixings = FALSE;
   (*consdata)->validsignature = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->normalized = FALSE;
   (*consdata)->upgradetried = FALSE;
   (*consdata)->upgraded = FALSE;
   (*consdata)->sorted = (nvars <= 1);
   (*consdata)->merged = (nvars <= 1);
   (*consdata)->cliquesadded = FALSE;
   (*consdata)->implsadded = FALSE;
   (*consdata)->binvarssorted = FALSE;
   (*consdata)->nbinvars = -1;
   (*consdata)->varsdeleted = FALSE;
   (*consdata)->rangedrowpropagated = 0;
   (*consdata)->checkabsolute = FALSE;

   if( SCIPisTransformed(scip) )
   {
      /* get transformed variables */
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   /* capture variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      assert(!SCIPisZero(scip, (*consdata)->vals[v]));
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[v]) );
   }

   return SCIP_OKAY;
}

/** frees a linear constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to linear constraint data */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->varssize >= 0);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* release variables */
   for( v = 0; v < (*consdata)->nvars; v++ )
   {
      assert((*consdata)->vars[v] != NULL);
      assert(!SCIPisZero(scip, (*consdata)->vals[v]));
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[v])) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vals, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** prints linear constraint in CIP format to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->nvars == 0 )
      SCIPinfoMessage(scip, file, "0");
   else
   {
      /* post linear sum of the linear constraint */
      SCIP_CALL( SCIPwriteVarsLinearsum(scip, file, consdata->vars, consdata->vals, consdata->nvars, TRUE) );
   }

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

/** prints linear constraint and contained solution values of variables to file stream */
static
SCIP_RETCODE consPrintConsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol,                /**< solution to print */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  [%s] <%s>: ", SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), SCIPconsGetName(cons));

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->nvars == 0 )
      SCIPinfoMessage(scip, file, "0");
   else
   {
      int v;

      /* post linear sum of the linear constraint */
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( consdata->vals != NULL )
         {
            if( consdata->vals[v] == 1.0 )
            {
               if( v > 0 )
                  SCIPinfoMessage(scip, file, " +");
            }
            else if( consdata->vals[v] == -1.0 )
               SCIPinfoMessage(scip, file, " -");
            else
               SCIPinfoMessage(scip, file, " %+.9g", consdata->vals[v]);
         }
         else if( consdata->nvars > 0 )
            SCIPinfoMessage(scip, file, " +");

         /* print variable name */
         SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[v], TRUE) );

         SCIPinfoMessage(scip, file, " (%+.9g)", SCIPgetSolVal(scip, sol, consdata->vars[v]));
      }
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   SCIPinfoMessage(scip, file, ";\n");

   return SCIP_OKAY;
}

/** invalidates activity bounds, such that they are recalculated in next get */
static
void consdataInvalidateActivities(
   SCIP_CONSDATA*        consdata            /**< linear constraint */
   )
{
   assert(consdata != NULL);

   consdata->validactivities = FALSE;
   consdata->validminact = FALSE;
   consdata->validmaxact = FALSE;
   consdata->validglbminact = FALSE;
   consdata->validglbmaxact = FALSE;
   consdata->validmaxabsval = FALSE;
   consdata->validminabsval = FALSE;
   consdata->hasnonbinvalid = FALSE;
   consdata->minactivity = SCIP_INVALID;
   consdata->maxactivity = SCIP_INVALID;
   consdata->lastminactivity = SCIP_INVALID;
   consdata->lastmaxactivity = SCIP_INVALID;
   consdata->maxabsval = SCIP_INVALID;
   consdata->minabsval = SCIP_INVALID;
   consdata->maxactdelta = SCIP_INVALID;
   consdata->maxactdeltavar = NULL;
   consdata->minactivityneginf = -1;
   consdata->minactivityposinf = -1;
   consdata->maxactivityneginf = -1;
   consdata->maxactivityposinf = -1;
   consdata->minactivityneghuge = -1;
   consdata->minactivityposhuge = -1;
   consdata->maxactivityneghuge = -1;
   consdata->maxactivityposhuge = -1;
   consdata->glbminactivity = SCIP_INVALID;
   consdata->glbmaxactivity = SCIP_INVALID;
   consdata->lastglbminactivity = SCIP_INVALID;
   consdata->lastglbmaxactivity = SCIP_INVALID;
   consdata->glbminactivityneginf = -1;
   consdata->glbminactivityposinf = -1;
   consdata->glbmaxactivityneginf = -1;
   consdata->glbmaxactivityposinf = -1;
   consdata->glbminactivityneghuge = -1;
   consdata->glbminactivityposhuge = -1;
   consdata->glbmaxactivityneghuge = -1;
   consdata->glbmaxactivityposhuge = -1;
}

/** compute the pseudo activity of a constraint */
static
SCIP_Real consdataComputePseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   int pseudoactivityposinf;
   int pseudoactivityneginf;
   SCIP_Real pseudoactivity;
   SCIP_Real bound;
   SCIP_Real val;

   pseudoactivity = 0;
   pseudoactivityposinf = 0;
   pseudoactivityneginf = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      val = consdata->vals[i];
      bound = (SCIPvarGetBestBoundType(consdata->vars[i]) == SCIP_BOUNDTYPE_LOWER) ? SCIPvarGetLbLocal(consdata->vars[i]) : SCIPvarGetUbLocal(consdata->vars[i]);
      if( SCIPisInfinity(scip, bound) )
      {
         if( val > 0.0 )
            pseudoactivityposinf++;
         else
            pseudoactivityneginf++;
      }
      else
      {
         if( SCIPisInfinity(scip, -bound) )
         {
            if( val > 0.0 )
               pseudoactivityneginf++;
            else
               pseudoactivityposinf++;
         }
         else
            pseudoactivity += val * bound;
      }
   }

   if( pseudoactivityneginf > 0 && pseudoactivityposinf > 0 )
      return SCIP_INVALID;
   else if( pseudoactivityneginf > 0 )
      return -SCIPinfinity(scip);
   else if( pseudoactivityposinf > 0 )
      return SCIPinfinity(scip);

   return pseudoactivity;
}

/** recompute the minactivity of a constraint */
static
void consdataRecomputeMinactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;

   consdata->minactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetLbLocal(consdata->vars[i]) : SCIPvarGetUbLocal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->vals[i] * bound) && !SCIPisHugeValue(scip, -consdata->vals[i] * bound) )
         consdata->minactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validminact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastminactivity = consdata->minactivity;
}

/** recompute the maxactivity of a constraint */
static
void consdataRecomputeMaxactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;

   consdata->maxactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetUbLocal(consdata->vars[i]) : SCIPvarGetLbLocal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->vals[i] * bound) && !SCIPisHugeValue(scip, -consdata->vals[i] * bound) )
         consdata->maxactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validmaxact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastmaxactivity = consdata->maxactivity;
}

/** recompute the global minactivity of a constraint */
static
void consdataRecomputeGlbMinactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;

   consdata->glbminactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetLbGlobal(consdata->vars[i]) : SCIPvarGetUbGlobal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->vals[i] * bound) && !SCIPisHugeValue(scip, -consdata->vals[i] * bound) )
         consdata->glbminactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validglbminact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastglbminactivity = consdata->glbminactivity;
}

/** recompute the global maxactivity of a constraint */
static
void consdataRecomputeGlbMaxactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;
   SCIP_Real bound;

   consdata->glbmaxactivity = 0;

   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      bound = (consdata->vals[i] > 0.0 ) ? SCIPvarGetUbGlobal(consdata->vars[i]) : SCIPvarGetLbGlobal(consdata->vars[i]);
      if( !SCIPisInfinity(scip, bound) && !SCIPisInfinity(scip, -bound)
         && !SCIPisHugeValue(scip, consdata->vals[i] * bound) && !SCIPisHugeValue(scip, -consdata->vals[i] * bound) )
         consdata->glbmaxactivity += consdata->vals[i] * bound;
   }

   /* the activity was just computed from scratch and is valid now */
   consdata->validglbmaxact = TRUE;

   /* the activity was just computed from scratch, mark it to be reliable */
   consdata->lastglbmaxactivity = consdata->glbmaxactivity;
}

/** calculates maximum absolute value of coefficients */
static
void consdataCalcMaxAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   SCIP_Real absval;
   int i;

   assert(consdata != NULL);
   assert(!consdata->validmaxabsval);
   assert(consdata->maxabsval >= SCIP_INVALID);

   consdata->validmaxabsval = TRUE;
   consdata->maxabsval = 0.0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      absval = consdata->vals[i];
      absval = REALABS(absval);
      if( absval > consdata->maxabsval )
         consdata->maxabsval = absval;
   }
}

/** calculates minimum absolute value of coefficients */
static
void consdataCalcMinAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   SCIP_Real absval;
   int i;

   assert(consdata != NULL);
   assert(!consdata->validminabsval);
   assert(consdata->minabsval >= SCIP_INVALID);

   consdata->validminabsval = TRUE;

   if( consdata->nvars > 0 )
      consdata->minabsval = REALABS(consdata->vals[0]);
   else
      consdata->minabsval = 0.0;

   for( i = 1; i < consdata->nvars; ++i )
   {
      absval = consdata->vals[i];
      absval = REALABS(absval);
      if( absval < consdata->minabsval )
         consdata->minabsval = absval;
   }
}

/** checks the type of all variables of the constraint and sets hasnonbinvar and hascontvar flags accordingly */
static
void consdataCheckNonbinvar(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int v;

   assert(!consdata->hasnonbinvalid);
   consdata->hasnonbinvar = FALSE;
   consdata->hascontvar = FALSE;

   for( v = consdata->nvars - 1; v >= 0; --v )
   {
      SCIP_VARTYPE vartype = SCIPvarGetType(consdata->vars[v]);

      if( vartype != SCIP_VARTYPE_BINARY )
      {
         consdata->hasnonbinvar = TRUE;

         if( vartype == SCIP_VARTYPE_CONTINUOUS )
         {
            consdata->hascontvar = TRUE;
            break;
         }
      }
   }
   assert(consdata->hascontvar || v < 0);

   consdata->hasnonbinvalid = TRUE;
}


#ifdef CHECKMAXACTDELTA
/* checks that the stored maximal activity delta (if not invalid) is correct */
static
void checkMaxActivityDelta(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   if( consdata->maxactdelta != SCIP_INVALID )
   {
      SCIP_Real maxactdelta = 0.0;
      SCIP_Real domain;
      SCIP_Real delta;
      SCIP_Real lb;
      SCIP_Real ub;
      int v;

      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         lb = SCIPvarGetLbLocal(consdata->vars[v]);
         ub = SCIPvarGetUbLocal(consdata->vars[v]);

         if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
         {
            maxactdelta = SCIPinfinity(scip);
            break;
         }

         domain = ub - lb;
         delta = REALABS(consdata->vals[v]) * domain;

         if( delta > maxactdelta )
         {
            maxactdelta = delta;
         }
      }
      assert(SCIPisFeasEQ(scip, maxactdelta, consdata->maxactdelta));
   }
}
#else
#define checkMaxActivityDelta(scip, consdata) /**/
#endif

/** recompute maximal activity contribution for a single variable */
static
void consdataRecomputeMaxActivityDelta(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   SCIP_Real delta;
   int v;

   consdata->maxactdelta = 0.0;

   if( !consdata->hasnonbinvalid )
      consdataCheckNonbinvar(consdata);

   /* easy case, the problem consists only of binary variables */
   if( !consdata->hasnonbinvar )
   {
      for( v = consdata->nvars - 1; v >= 0; --v )
      {
         if( SCIPvarGetLbLocal(consdata->vars[v]) < 0.5 && SCIPvarGetUbLocal(consdata->vars[v]) > 0.5 )
         {
            delta = REALABS(consdata->vals[v]);

            if( delta > consdata->maxactdelta )
            {
               consdata->maxactdelta = delta;
               consdata->maxactdeltavar = consdata->vars[v];
            }
         }
      }
      return;
   }

   for( v = consdata->nvars - 1; v >= 0; --v )
   {
      SCIP_Real domain;
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(consdata->vars[v]);
      ub = SCIPvarGetUbLocal(consdata->vars[v]);

      if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
      {
         consdata->maxactdelta = SCIPinfinity(scip);
         consdata->maxactdeltavar = consdata->vars[v];
         break;
      }

      domain = ub - lb;
      delta = REALABS(consdata->vals[v]) * domain;

      if( delta > consdata->maxactdelta )
      {
         consdata->maxactdelta = delta;
         consdata->maxactdeltavar = consdata->vars[v];
      }
   }
}


/** updates activities for a change in a bound */
static
void consdataUpdateActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed; can be NULL for global bound changes */
   SCIP_Real             oldbound,           /**< old bound of variable */
   SCIP_Real             newbound,           /**< new bound of variable */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_BOUNDTYPE        boundtype,          /**< type of the bound change */
   SCIP_Bool             global,             /**< is it a global or a local bound change? */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   SCIP_Real* activity;
   SCIP_Real* lastactivity;
   int* activityposinf;
   int* activityneginf;
   int* activityposhuge;
   int* activityneghuge;
   SCIP_Real oldcontribution;
   SCIP_Real newcontribution;
   SCIP_Real delta;
   SCIP_Bool validact;
   SCIP_Bool finitenewbound;
   SCIP_Bool hugevalnewcont;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(global || (var != NULL));
   assert(consdata->validactivities);
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->lastminactivity < SCIP_INVALID);
   assert(consdata->lastmaxactivity < SCIP_INVALID);
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);
   assert(consdata->minactivityneghuge >= 0);
   assert(consdata->minactivityposhuge >= 0);
   assert(consdata->maxactivityneghuge >= 0);
   assert(consdata->maxactivityposhuge >= 0);
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->lastglbminactivity < SCIP_INVALID);
   assert(consdata->lastglbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);
   assert(consdata->glbminactivityneghuge >= 0);
   assert(consdata->glbminactivityposhuge >= 0);
   assert(consdata->glbmaxactivityneghuge >= 0);
   assert(consdata->glbmaxactivityposhuge >= 0);

   delta = 0.0;

   /* we are updating global activities */
   if( global )
   {
      /* depending on the boundtype and the coefficient, we choose the activity to be updated:
       * lower bound + pos. coef: update minactivity
       * lower bound + neg. coef: update maxactivity, positive and negative infinity counters have to be switched
       * upper bound + pos. coef: update maxactivity
       * upper bound + neg. coef: update minactivity, positive and negative infinity counters have to be switched
       */
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         if( val > 0.0 )
         {
            activity = &(consdata->glbminactivity);
            lastactivity = &(consdata->lastglbminactivity);
            activityposinf = &(consdata->glbminactivityposinf);
            activityneginf = &(consdata->glbminactivityneginf);
            activityposhuge = &(consdata->glbminactivityposhuge);
            activityneghuge = &(consdata->glbminactivityneghuge);
            validact = consdata->validglbminact;
         }
         else
         {
            activity = &(consdata->glbmaxactivity);
            lastactivity = &(consdata->lastglbmaxactivity);
            activityposinf = &(consdata->glbmaxactivityneginf);
            activityneginf = &(consdata->glbmaxactivityposinf);
            activityposhuge = &(consdata->glbmaxactivityposhuge);
            activityneghuge = &(consdata->glbmaxactivityneghuge);
            validact = consdata->validglbmaxact;
         }
      }
      else
      {
         if( val > 0.0 )
         {
            activity = &(consdata->glbmaxactivity);
            lastactivity = &(consdata->lastglbmaxactivity);
            activityposinf = &(consdata->glbmaxactivityposinf);
            activityneginf = &(consdata->glbmaxactivityneginf);
            activityposhuge = &(consdata->glbmaxactivityposhuge);
            activityneghuge = &(consdata->glbmaxactivityneghuge);
            validact = consdata->validglbmaxact;
         }
         else
         {
            activity = &(consdata->glbminactivity);
            lastactivity = &(consdata->lastglbminactivity);
            activityposinf = &(consdata->glbminactivityneginf);
            activityneginf = &(consdata->glbminactivityposinf);
            activityposhuge = &(consdata->glbminactivityposhuge);
            activityneghuge = &(consdata->glbminactivityneghuge);
            validact = consdata->validglbminact;
         }
      }
   }
   /* we are updating local activities */
   else
   {
      /* depending on the boundtype and the coefficient, we choose the activity to be updated:
       * lower bound + pos. coef: update minactivity
       * lower bound + neg. coef: update maxactivity, positive and negative infinity counters have to be switched
       * upper bound + pos. coef: update maxactivity
       * upper bound + neg. coef: update minactivity, positive and negative infinity counters have to be switched
       */
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
      {
         if( val > 0.0 )
         {
            activity = &(consdata->minactivity);
            lastactivity = &(consdata->lastminactivity);
            activityposinf = &(consdata->minactivityposinf);
            activityneginf = &(consdata->minactivityneginf);
            activityposhuge = &(consdata->minactivityposhuge);
            activityneghuge = &(consdata->minactivityneghuge);
            validact = consdata->validminact;
         }
         else
         {
            activity = &(consdata->maxactivity);
            lastactivity = &(consdata->lastmaxactivity);
            activityposinf = &(consdata->maxactivityneginf);
            activityneginf = &(consdata->maxactivityposinf);
            activityposhuge = &(consdata->maxactivityposhuge);
            activityneghuge = &(consdata->maxactivityneghuge);
            validact = consdata->validmaxact;
         }
      }
      else
      {
         if( val > 0.0 )
         {
            activity = &(consdata->maxactivity);
            lastactivity = &(consdata->lastmaxactivity);
            activityposinf = &(consdata->maxactivityposinf);
            activityneginf = &(consdata->maxactivityneginf);
            activityposhuge = &(consdata->maxactivityposhuge);
            activityneghuge = &(consdata->maxactivityneghuge);
            validact = consdata->validmaxact;
         }
         else
         {
            activity = &(consdata->minactivity);
            lastactivity = &(consdata->lastminactivity);
            activityposinf = &(consdata->minactivityneginf);
            activityneginf = &(consdata->minactivityposinf);
            activityposhuge = &(consdata->minactivityposhuge);
            activityneghuge = &(consdata->minactivityneghuge);
            validact = consdata->validminact;
         }
      }
   }

   oldcontribution = val * oldbound;
   newcontribution = val * newbound;
   hugevalnewcont = SCIPisHugeValue(scip, REALABS(newcontribution));
   finitenewbound = !SCIPisInfinity(scip, REALABS(newbound));

   if( SCIPisInfinity(scip, REALABS(oldbound)) )
   {
      /* old bound was +infinity */
      if( oldbound > 0.0 )
      {
         assert((*activityposinf) >= 1);

         /* we only have to do something if the new bound is not again +infinity */
         if( finitenewbound || newbound < 0.0 )
         {
            /* decrease the counter for positive infinite contributions */
            (*activityposinf)--;

            /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
            if( !finitenewbound && newbound < 0.0 )
               (*activityneginf)++;
            else if( hugevalnewcont )
            {
               /* if the contribution of this variable is too large, increase the counter for huge values */
               if( newcontribution > 0.0 )
                  (*activityposhuge)++;
               else
                  (*activityneghuge)++;
            }
            /* "normal case": just add the contribution to the activity */
            else
               delta = newcontribution;
         }
      }
      /* old bound was -infinity */
      else
      {
         assert(oldbound < 0.0);
         assert((*activityneginf) >= 1);

         /* we only have to do something ig the new bound is not again -infinity */
         if( finitenewbound || newbound > 0.0 )
         {
            /* decrease the counter for negative infinite contributions */
            (*activityneginf)--;

            /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
            if( !finitenewbound && newbound > 0.0 )
               (*activityposinf)++;
            else if( hugevalnewcont )
            {
               /* if the contribution of this variable is too large, increase the counter for huge values */
               if( newcontribution > 0.0 )
                  (*activityposhuge)++;
               else
                  (*activityneghuge)++;
            }
            /* "normal case": just add the contribution to the activity */
            else
               delta = newcontribution;
         }
      }
   }
   else if( SCIPisHugeValue(scip, REALABS(oldcontribution)) )
   {
      /* old contribution was too large and positive */
      if( oldcontribution > 0.0 )
      {
         assert((*activityposhuge) >= 1);

         /* decrease the counter for huge positive contributions; it might be increased again later,
          * but checking here that the bound is not huge again would not handle a change from a huge to an infinite bound
          */
         (*activityposhuge)--;

         if( !finitenewbound )
         {
            /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
            if( newbound > 0.0 )
               (*activityposinf)++;
            /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
            else
               (*activityneginf)++;
         }
         else if( hugevalnewcont )
         {
            /* if the contribution of this variable is too large and positive, increase the corresponding counter */
            if( newcontribution > 0.0 )
               (*activityposhuge)++;
            /* if the contribution of this variable is too large and negative, increase the corresponding counter */
            else
               (*activityneghuge)++;
         }
         /* "normal case": just add the contribution to the activity */
         else
            delta = newcontribution;
      }
      /* old contribution was too large and negative */
      else
      {
         assert(oldcontribution < 0.0);
         assert((*activityneghuge) >= 1);

         /* decrease the counter for huge negative contributions; it might be increased again later,
          * but checking here that the bound is not huge again would not handle a change from a huge to an infinite bound
          */
         (*activityneghuge)--;

         if( !finitenewbound )
         {
            /* if the bound changed to +infinity, increase the counter for positive infinite contributions */
            if( newbound > 0.0 )
               (*activityposinf)++;
            /* if the bound changed to -infinity, increase the counter for negative infinite contributions */
            else
               (*activityneginf)++;
         }
         else if( hugevalnewcont )
         {
            /* if the contribution of this variable is too large and positive, increase the corresponding counter */
            if( newcontribution > 0.0 )
               (*activityposhuge)++;
            /* if the contribution of this variable is too large and negative, increase the corresponding counter */
            else
               (*activityneghuge)++;
         }
         /* "normal case": just add the contribution to the activity */
         else
            delta = newcontribution;
      }
   }
   /* old bound was finite and not too large */
   else
   {
      if( !finitenewbound )
      {
         /* if the new bound is +infinity, the old contribution has to be subtracted
          * and the counter for positive infinite contributions has to be increased
          */
         if( newbound > 0.0 )
         {
            (*activityposinf)++;
            delta = -oldcontribution;
         }
         /* if the new bound is -infinity, the old contribution has to be subtracted
          * and the counter for negative infinite contributions has to be increased
          */
         else
         {
            assert(newbound < 0.0 );

            (*activityneginf)++;
            delta = -oldcontribution;
         }
      }
      /* if the contribution of this variable is too large, increase the counter for huge values */
      else if( hugevalnewcont )
      {
         if( newcontribution > 0.0 )
         {
            (*activityposhuge)++;
            delta = -oldcontribution;
         }
         else
         {
            (*activityneghuge)++;
            delta = -oldcontribution;
         }
      }
      /* "normal case": just update the activity */
      else
         delta = newcontribution - oldcontribution;
   }

   /* update the activity, if the current value is valid and there was a change in the finite part */
   if( validact && (delta != 0.0) )
   {
      /* if the absolute value of the activity is increased, this is regarded as reliable,
       * otherwise, we check whether we can still trust the updated value
       */
      (*activity) = (*activity) + delta;
      assert(!SCIPisInfinity(scip, -(*activity)) && !SCIPisInfinity(scip, *activity));

      if( REALABS((*lastactivity)) < REALABS(*activity) )
      {
         (*lastactivity) = (*activity);
      }
      else
      {
         if( checkreliability && SCIPisUpdateUnreliable(scip, (*activity), (*lastactivity)) )
         {
            SCIPdebugMsg(scip, "%s activity of linear constraint unreliable after update: %16.9g\n",
               (global ? "global " : ""), (*activity));

            /* mark the activity that was just changed and is not reliable anymore to be invalid */
            if( global )
            {
               if( (boundtype == SCIP_BOUNDTYPE_LOWER) == (val > 0.0) )
                  consdata->validglbminact = FALSE;
               else
                  consdata->validglbmaxact = FALSE;
            }
            else
            {
               if( (boundtype == SCIP_BOUNDTYPE_LOWER) == (val > 0.0) )
                  consdata->validminact = FALSE;
               else
                  consdata->validmaxact = FALSE;
            }
         }
      }
   }
}

/** updates minimum and maximum activity for a change in lower bound */
static
void consdataUpdateActivitiesLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, var, oldlb, newlb, val, SCIP_BOUNDTYPE_LOWER, FALSE, checkreliability);

      assert(!SCIPisInfinity(scip, -consdata->minactivity) && !SCIPisInfinity(scip, consdata->minactivity));
      assert(!SCIPisInfinity(scip, -consdata->maxactivity) && !SCIPisInfinity(scip, consdata->maxactivity));
   }
}

/** updates minimum and maximum activity for a change in upper bound */
static
void consdataUpdateActivitiesUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable that has been changed */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub,              /**< new upper bound of variable */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, var, oldub, newub, val, SCIP_BOUNDTYPE_UPPER, FALSE, checkreliability);

      assert(!SCIPisInfinity(scip, -consdata->minactivity) && !SCIPisInfinity(scip, consdata->minactivity));
      assert(!SCIPisInfinity(scip, -consdata->maxactivity) && !SCIPisInfinity(scip, consdata->maxactivity));
   }
}

/** updates minimum and maximum global activity for a change in the global lower bound */
static
void consdataUpdateActivitiesGlbLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real             oldlb,              /**< old lower bound of variable */
   SCIP_Real             newlb,              /**< new lower bound of variable */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, NULL, oldlb, newlb, val, SCIP_BOUNDTYPE_LOWER, TRUE, checkreliability);

      assert(!SCIPisInfinity(scip, -consdata->glbminactivity) && !SCIPisInfinity(scip, consdata->glbminactivity));
      assert(!SCIPisInfinity(scip, -consdata->glbmaxactivity) && !SCIPisInfinity(scip, consdata->glbmaxactivity));
   }
}

/** updates minimum and maximum global activity for a change in global upper bound */
static
void consdataUpdateActivitiesGlbUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real             oldub,              /**< old upper bound of variable */
   SCIP_Real             newub,              /**< new upper bound of variable */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->validactivities )
   {
      consdataUpdateActivities(scip, consdata, NULL, oldub, newub, val, SCIP_BOUNDTYPE_UPPER, TRUE, checkreliability);

      assert(!SCIPisInfinity(scip, -consdata->glbminactivity) && !SCIPisInfinity(scip, consdata->glbminactivity));
      assert(!SCIPisInfinity(scip, -consdata->glbmaxactivity) && !SCIPisInfinity(scip, consdata->glbmaxactivity));
   }
}

/** updates minimum and maximum activity and maximum absolute value for coefficient addition */
static
void consdataUpdateAddCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      assert(consdata->maxabsval < SCIP_INVALID);

      absval = REALABS(val);
      consdata->maxabsval = MAX(consdata->maxabsval, absval);
   }

   if( consdata->validminabsval )
   {
      SCIP_Real absval;

      assert(consdata->minabsval < SCIP_INVALID);

      absval = REALABS(val);
      consdata->minabsval = MIN(consdata->minabsval, absval);
   }

   /* update minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateActivitiesLb(scip, consdata, var, 0.0, SCIPvarGetLbLocal(var), val, checkreliability);
      consdataUpdateActivitiesUb(scip, consdata, var, 0.0, SCIPvarGetUbLocal(var), val, checkreliability);
      consdataUpdateActivitiesGlbLb(scip, consdata, 0.0, SCIPvarGetLbGlobal(var), val, checkreliability);
      consdataUpdateActivitiesGlbUb(scip, consdata, 0.0, SCIPvarGetUbGlobal(var), val, checkreliability);
   }
}

/** updates minimum and maximum activity for coefficient deletion, invalidates maximum absolute value if necessary */
static
void consdataUpdateDelCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val,                /**< coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* invalidate maximum absolute value, if this coefficient was the maximum */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      absval = REALABS(val);

      if( SCIPisEQ(scip, absval, consdata->maxabsval) )
      {
         consdata->validmaxabsval = FALSE;
         consdata->maxabsval = SCIP_INVALID;
      }
   }

   /* invalidate minimum absolute value, if this coefficient was the minimum */
   if( consdata->validminabsval )
   {
      SCIP_Real absval;

      absval = REALABS(val);

      if( SCIPisEQ(scip, absval, consdata->minabsval) )
      {
         consdata->validminabsval = FALSE;
         consdata->minabsval = SCIP_INVALID;
      }
   }

   /* update minimal and maximal activity */
   if( consdata->validactivities )
   {
      assert(consdata->minactivity < SCIP_INVALID);
      assert(consdata->maxactivity < SCIP_INVALID);
      assert(consdata->glbminactivity < SCIP_INVALID);
      assert(consdata->glbmaxactivity < SCIP_INVALID);

      consdataUpdateActivitiesLb(scip, consdata, var, SCIPvarGetLbLocal(var), 0.0, val, checkreliability);
      consdataUpdateActivitiesUb(scip, consdata, var, SCIPvarGetUbLocal(var), 0.0, val, checkreliability);
      consdataUpdateActivitiesGlbLb(scip, consdata, SCIPvarGetLbGlobal(var), 0.0, val, checkreliability);
      consdataUpdateActivitiesGlbUb(scip, consdata, SCIPvarGetUbGlobal(var), 0.0, val, checkreliability);
   }
}

/** updates minimum and maximum activity for coefficient change, invalidates maximum absolute value if necessary */
static
void consdataUpdateChgCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             oldval,             /**< old coefficient of constraint entry */
   SCIP_Real             newval,             /**< new coefficient of constraint entry */
   SCIP_Bool             checkreliability    /**< should the reliability of the recalculated activity be checked? */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);

   /* old zero coefficients should be handled by consdataUpdateAddCoef() */
   assert(!SCIPisZero(scip, oldval));

   /* new zero coefficients should be handled by consdataUpdateDelCoef() */
   assert(!SCIPisZero(scip, newval));

   /* update maximum absolute value */
   if( consdata->validmaxabsval )
   {
      SCIP_Real absval;

      absval = REALABS(newval);

      if( SCIPisGE(scip, absval, consdata->maxabsval) )
      {
         consdata->maxabsval = absval;
      }
      else
      {
         absval = REALABS(oldval);

         /* invalidate maximum absolute value */
         if( SCIPisEQ(scip, absval, consdata->maxabsval) )
         {
            consdata->validmaxabsval = FALSE;
            consdata->maxabsval = SCIP_INVALID;
         }
      }
   }

   /* update minimum absolute value */
   if( consdata->validminabsval )
   {
      SCIP_Real absval;

      absval = REALABS(newval);

      if( SCIPisLE(scip, absval, consdata->minabsval) )
      {
         consdata->minabsval = absval;
      }
      else
      {
         absval = REALABS(oldval);

         /* invalidate minimum absolute value */
         if( SCIPisEQ(scip, absval, consdata->minabsval) )
         {
            consdata->validminabsval = FALSE;
            consdata->minabsval = SCIP_INVALID;
         }
      }
   }


   /* update maximum activity delta */
   if( !SCIPisInfinity(scip, consdata->maxactdelta ) )
   {
      SCIP_Real domain;
      SCIP_Real delta;

      assert(!SCIPisInfinity(scip, SCIPvarGetLbLocal(var)));
      assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));

      domain = SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var);
      delta = REALABS(newval) * domain;

      if( delta > consdata->maxactdelta )
      {
         consdata->maxactdelta = delta;
         consdata->maxactdeltavar = var;
      }
      else
      {
         /* reset maximal activity delta, so that it will be recalculated on the next real propagation */
         if( consdata->maxactdeltavar == var )
            consdata->maxactdelta = SCIP_INVALID;
      }
   }

   /* @todo do something more clever here, e.g. if oldval * newval >= 0, do the update directly */
   consdataUpdateDelCoef(scip, consdata, var, oldval, checkreliability);
   consdataUpdateAddCoef(scip, consdata, var, newval, checkreliability);
}

/** returns the maximum absolute value of all coefficients in the constraint */
static
SCIP_Real consdataGetMaxAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validmaxabsval )
      consdataCalcMaxAbsval(consdata);
   assert(consdata->validmaxabsval);
   assert(consdata->maxabsval < SCIP_INVALID);

   return consdata->maxabsval;
}

/** returns the minimum absolute value of all coefficients in the constraint */
static
SCIP_Real consdataGetMinAbsval(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validminabsval )
      consdataCalcMinAbsval(consdata);
   assert(consdata->validminabsval);
   assert(consdata->minabsval < SCIP_INVALID);

   return consdata->minabsval;
}

/** calculates minimum and maximum local and global activity for constraint from scratch;
 *  additionally recalculates maximum absolute value of coefficients
 */
static
void consdataCalcActivities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(!consdata->validactivities);
   assert(consdata->minactivity >= SCIP_INVALID || consdata->validminact);
   assert(consdata->maxactivity >= SCIP_INVALID || consdata->validmaxact);
   assert(consdata->glbminactivity >= SCIP_INVALID || consdata->validglbminact);
   assert(consdata->glbmaxactivity >= SCIP_INVALID || consdata->validglbmaxact);

   consdata->validmaxabsval = TRUE;
   consdata->validminabsval = TRUE;
   consdata->validactivities = TRUE;
   consdata->validminact = TRUE;
   consdata->validmaxact = TRUE;
   consdata->validglbminact = TRUE;
   consdata->validglbmaxact = TRUE;
   consdata->maxabsval = 0.0;
   consdata->minabsval = (consdata->nvars == 0 ? 0.0 : REALABS(consdata->vals[0]));
   consdata->minactivity = 0.0;
   consdata->maxactivity = 0.0;
   consdata->lastminactivity = 0.0;
   consdata->lastmaxactivity = 0.0;
   consdata->minactivityneginf = 0;
   consdata->minactivityposinf = 0;
   consdata->maxactivityneginf = 0;
   consdata->maxactivityposinf = 0;
   consdata->minactivityneghuge = 0;
   consdata->minactivityposhuge = 0;
   consdata->maxactivityneghuge = 0;
   consdata->maxactivityposhuge = 0;
   consdata->glbminactivity = 0.0;
   consdata->glbmaxactivity = 0.0;
   consdata->lastglbminactivity = 0.0;
   consdata->lastglbmaxactivity = 0.0;
   consdata->glbminactivityneginf = 0;
   consdata->glbminactivityposinf = 0;
   consdata->glbmaxactivityneginf = 0;
   consdata->glbmaxactivityposinf = 0;
   consdata->glbminactivityneghuge = 0;
   consdata->glbminactivityposhuge = 0;
   consdata->glbmaxactivityneghuge = 0;
   consdata->glbmaxactivityposhuge = 0;

   for( i = 0; i < consdata->nvars; ++i )
      consdataUpdateAddCoef(scip, consdata, consdata->vars[i], consdata->vals[i], FALSE);

   consdata->lastminactivity = consdata->minactivity;
   consdata->lastmaxactivity = consdata->maxactivity;
   consdata->lastglbminactivity = consdata->glbminactivity;
   consdata->lastglbmaxactivity = consdata->glbmaxactivity;
}

/** gets minimal activity for constraint and given values of counters for infinite and huge contributions
 *  and (if needed) delta to subtract from stored finite part of activity in case of a residual activity
 */
static
void getMinActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   int                   posinf,             /**< number of coefficients contributing pos. infinite value */
   int                   neginf,             /**< number of coefficients contributing neg. infinite value */
   int                   poshuge,            /**< number of coefficients contributing huge pos. value */
   int                   neghuge,            /**< number of coefficients contributing huge neg. value */
   SCIP_Real             delta,              /**< value to subtract from stored minactivity
                                              *   (contribution of the variable set to zero when getting residual activity) */
   SCIP_Bool             global,             /**< should the global or local minimal activity be returned? */
   SCIP_Bool             goodrelax,          /**< should a good relaxation be computed or are relaxed acticities ignored, anyway? */
   SCIP_Real*            minactivity,        /**< pointer to store the minimal activity */
   SCIP_Bool*            isrelax,            /**< pointer to store whether the activity is a relaxation,
                                              *   i.e. is <= the exact minactivity (in case of huge contributing values) */
   SCIP_Bool*            issettoinfinity     /**< pointer to store whether minactivity was set to infinity or calculated */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(posinf >= 0);
   assert(neginf >= 0);
   assert(poshuge >= 0);
   assert(neghuge >= 0);
   assert(minactivity != NULL);
   assert(isrelax != NULL);
   assert(issettoinfinity != NULL);

   /* if we have pos. infinite contributions, the minactivity is +infty */
   if( posinf > 0 )
   {
      *minactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have neg. (and no pos.) infinite contributions, the minactivity is -infty */
   else if( neginf > 0 )
   {
      *minactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have neg. huge contributions, we only know that -infty is a relaxation of the minactivity */
   else if( neghuge > 0 )
   {
      *minactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   /* we do not need a good relaxation and we have positve huge contributions, so we just return -infty as activity */
   else if( !goodrelax && poshuge > 0 )
   {
      *minactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   else
   {
      SCIP_Real tmpactivity;

      /* recompute minactivity if it is not valid */
      if( global )
      {
         if( !consdata->validglbminact )
            consdataRecomputeGlbMinactivity(scip, consdata);
         assert(consdata->validglbminact);

         tmpactivity = consdata->glbminactivity;
      }
      else
      {
         if( !consdata->validminact )
            consdataRecomputeMinactivity(scip, consdata);
         assert(consdata->validminact);

         tmpactivity = consdata->minactivity;
      }

      /* we have no infinite and no neg. huge contributions, but pos. huge contributions;
       * a feasible relaxation of the minactivity is the number of positive huge contributions
       * times the minimum value counting as "huge" plus finite (and non-huge) part of minactivity - delta
       */
      if( poshuge > 0 )
      {
         *minactivity = 1.0 * poshuge * SCIPgetHugeValue(scip) + (tmpactivity - delta);
         *issettoinfinity = FALSE;
         *isrelax = TRUE;
      }
      /* all counters are zero, so the minactivity is just stored and we subtract the delta */
      else
      {
         *minactivity = tmpactivity - delta;
         *issettoinfinity = FALSE;
         *isrelax = FALSE;
      }
   }
}

/** gets maximal activity for constraint and given values of counters for infinite and huge contributions
 *  and (if needed) delta to subtract from stored finite part of activity in case of a residual activity
 */
static
void getMaxActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   int                   posinf,             /**< number of coefficients contributing pos. infinite value */
   int                   neginf,             /**< number of coefficients contributing neg. infinite value */
   int                   poshuge,            /**< number of coefficients contributing huge pos. value */
   int                   neghuge,            /**< number of coefficients contributing huge neg. value */
   SCIP_Real             delta,              /**< value to subtract from stored maxactivity
                                              *   (contribution of the variable set to zero when getting residual activity) */
   SCIP_Bool             global,             /**< should the global or local maximal activity be returned? */
   SCIP_Bool             goodrelax,          /**< should a good relaxation be computed or are relaxed acticities ignored, anyway? */
   SCIP_Real*            maxactivity,        /**< pointer to store the maximal activity */
   SCIP_Bool*            isrelax,            /**< pointer to store whether the activity is a relaxation,
                                              *   i.e. is >= the exact maxactivity (in case of huge contributing values) */
   SCIP_Bool*            issettoinfinity     /**< pointer to store whether maxactivity was set to infinity or calculated */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(posinf >= 0);
   assert(neginf >= 0);
   assert(poshuge >= 0);
   assert(neghuge >= 0);
   assert(maxactivity != NULL);
   assert(isrelax != NULL);
   assert(issettoinfinity != NULL);

   /* if we have neg. infinite contributions, the maxactivity is -infty */
   if( neginf > 0 )
   {
      *maxactivity = -SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have pos. (and no neg.) infinite contributions, the maxactivity is +infty */
   else if( posinf > 0 )
   {
      *maxactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = FALSE;
   }
   /* if we have pos. huge contributions, we only know that +infty is a relaxation of the maxactivity */
   else if( poshuge > 0 )
   {
      *maxactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   /* we do not need a good relaxation and we have positve huge contributions, so we just return +infty as activity */
   else if( !goodrelax && neghuge > 0 )
   {
      *maxactivity = SCIPinfinity(scip);
      *issettoinfinity = TRUE;
      *isrelax = TRUE;
   }
   else
   {
      SCIP_Real tmpactivity;

      /* recompute maxactivity if it is not valid */
      if( global )
      {
         if( !consdata->validglbmaxact )
            consdataRecomputeGlbMaxactivity(scip, consdata);
         assert(consdata->validglbmaxact);

         tmpactivity = consdata->glbmaxactivity;
      }
      else
      {
         if( !consdata->validmaxact )
            consdataRecomputeMaxactivity(scip, consdata);
         assert(consdata->validmaxact);

         tmpactivity = consdata->maxactivity;
      }

      /* we have no infinite, and no pos. huge contributions, but neg. huge contributions;
       * a feasible relaxation of the maxactivity is minus the number of negative huge contributions
       * times the minimum value counting as "huge" plus the finite (and non-huge) part of maxactivity minus delta
       */
      if( neghuge > 0 )
      {
         *maxactivity = -1.0 * neghuge * SCIPgetHugeValue(scip) + tmpactivity - delta;
         *issettoinfinity = FALSE;
         *isrelax = TRUE;
      }
      /* all counters are zero, so the maxactivity is just stored and we subtract the delta */
      else
      {
         *maxactivity = tmpactivity - delta;
         *issettoinfinity = FALSE;
         *isrelax = FALSE;
      }
   }
}

/** gets activity bounds for constraint */
static
void consdataGetActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            minactivity,        /**< pointer to store the minimal activity */
   SCIP_Real*            maxactivity,        /**< pointer to store the maximal activity */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned minactivity is just a relaxation,
                                              *   i.e. <= the exact minactivity (in case of huge contributions),
                                              *   or equal to the exact minimal activity */
   SCIP_Bool*            maxisrelax          /**< pointer to store whether the returned maxactivity is just a relaxation,
                                              *   i.e. >= the exact maxactivity (in case of huge contributions),
                                              *   or equal to the exact maximal activity */
   )
{
   SCIP_Bool issettoinfinity;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(minactivity != NULL);
   assert(maxactivity != NULL);

   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validminact);
      assert(consdata->validmaxact);
   }
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);

   getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
      consdata->minactivityposhuge, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
      minactivity, minisrelax, &issettoinfinity);

   getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
      consdata->maxactivityposhuge, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
      maxactivity, maxisrelax, &issettoinfinity);
}

/** calculates activity bounds for constraint after setting variable to zero */
static
void consdataGetReliableResidualActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             cancelvar,          /**< variable to calculate activity residual for */
   SCIP_Real*            resactivity,        /**< pointer to store the residual activity */
   SCIP_Bool             isminresact,        /**< should minimal or maximal residual activity be calculated? */
   SCIP_Bool             useglobalbounds     /**< should global or local bounds be used? */
   )
{
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(cancelvar != NULL);
   assert(resactivity != NULL);

   *resactivity = 0.0;

   for( v = 0; v < consdata->nvars; ++v )
   {
      var = consdata->vars[v];
      assert(var != NULL);
      if( var == cancelvar )
         continue;

      val = consdata->vals[v];

      if( useglobalbounds )
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
      }
      else
      {
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
      }

      assert(!SCIPisZero(scip, val));
      assert(SCIPisLE(scip, lb, ub));

      if( val > 0.0 )
      {
         if( isminresact )
         {
            assert(!SCIPisInfinity(scip, -lb));
            assert(!SCIPisHugeValue(scip, REALABS(val*lb)));
            *resactivity += val*lb;
         }
         else
         {
            assert(!SCIPisInfinity(scip, ub));
            assert(!SCIPisHugeValue(scip, REALABS(val*ub)));
            *resactivity += val*ub;
         }
      }
      else
      {
         if( isminresact)
         {
            assert(!SCIPisInfinity(scip, ub));
            assert(!SCIPisHugeValue(scip, REALABS(val*ub)));
            *resactivity += val*ub;
         }
         else
         {
            assert(!SCIPisInfinity(scip, -lb));
            assert(!SCIPisHugeValue(scip, REALABS(val*lb)));
            *resactivity += val*lb;
         }
      }
   }
   assert(!SCIPisInfinity(scip, *resactivity) && !SCIPisInfinity(scip, -(*resactivity)));
}

/** gets activity bounds for constraint after setting variable to zero */
static
void consdataGetActivityResiduals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             var,                /**< variable to calculate activity residual for */
   SCIP_Real             val,                /**< coefficient value of variable in linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   SCIP_Real*            maxresactivity,     /**< pointer to store the maximal residual activity */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned residual minactivity is just a
                                              *   relaxation, i.e. <= the exact residual minactivity (in case of huge
                                              *   contributions), or equal to the exact residual minactivity */
   SCIP_Bool*            maxisrelax,         /**< pointer to store whether the returned residual maxactivity is just a
                                              *   relaxation, i.e. <= the exact residual maxactivity (in case of huge
                                              *   contributions), or equal to the exact residual minactivity */
   SCIP_Bool*            isminsettoinfinity, /**< pointer to store whether minresactivity was set to infinity or calculated */
   SCIP_Bool*            ismaxsettoinfinity  /**< pointer to store whether maxresactivity was set to infinity or calculated */
   )
{
   SCIP_Real minactbound;
   SCIP_Real maxactbound;
   SCIP_Real absval;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);
   assert(minisrelax != NULL);
   assert(maxisrelax != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   /* get activity bounds of linear constraint */
   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validminact);
      assert(consdata->validmaxact);
   }
   assert(consdata->minactivity < SCIP_INVALID);
   assert(consdata->maxactivity < SCIP_INVALID);
   assert(consdata->minactivityneginf >= 0);
   assert(consdata->minactivityposinf >= 0);
   assert(consdata->maxactivityneginf >= 0);
   assert(consdata->maxactivityposinf >= 0);
   assert(consdata->minactivityneghuge >= 0);
   assert(consdata->minactivityposhuge >= 0);
   assert(consdata->maxactivityneghuge >= 0);
   assert(consdata->maxactivityposhuge >= 0);

   if( val > 0.0 )
   {
      minactbound = SCIPvarGetLbLocal(var);
      maxactbound = SCIPvarGetUbLocal(var);
      absval = val;
   }
   else
   {
      minactbound = -SCIPvarGetUbLocal(var);
      maxactbound = -SCIPvarGetLbLocal(var);
      absval = -val;
   }

   /* get/compute minactivity by calling getMinActivity() with updated counters for infinite and huge values
    * and contribution of variable set to zero that has to be subtracted from finite part of activity
    */
   if( SCIPisInfinity(scip, minactbound) )
   {
      assert(consdata->minactivityposinf >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf - 1, consdata->minactivityneginf,
         consdata->minactivityposhuge, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else if( SCIPisInfinity(scip, -minactbound) )
   {
      assert(consdata->minactivityneginf >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf - 1,
         consdata->minactivityposhuge, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, minactbound * absval) )
   {
      assert(consdata->minactivityposhuge >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
         consdata->minactivityposhuge - 1, consdata->minactivityneghuge, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, -minactbound * absval) )
   {
      assert(consdata->minactivityneghuge >= 1);

      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
         consdata->minactivityposhuge, consdata->minactivityneghuge - 1, 0.0, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }
   else
   {
      getMinActivity(scip, consdata, consdata->minactivityposinf, consdata->minactivityneginf,
         consdata->minactivityposhuge, consdata->minactivityneghuge, absval * minactbound, FALSE, goodrelax,
         minresactivity, minisrelax, isminsettoinfinity);
   }

   /* get/compute maxactivity by calling getMaxActivity() with updated counters for infinite and huge values
    * and contribution of variable set to zero that has to be subtracted from finite part of activity
    */
   if( SCIPisInfinity(scip, -maxactbound) )
   {
      assert(consdata->maxactivityneginf >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf - 1,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else if( SCIPisInfinity(scip, maxactbound) )
   {
      assert(consdata->maxactivityposinf >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf - 1, consdata->maxactivityneginf,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, absval * maxactbound) )
   {
      assert(consdata->maxactivityposhuge >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
         consdata->maxactivityposhuge - 1, consdata->maxactivityneghuge, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else if( SCIPisHugeValue(scip, -absval * maxactbound) )
   {
      assert(consdata->maxactivityneghuge >= 1);

      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge - 1, 0.0, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
   else
   {
      getMaxActivity(scip, consdata, consdata->maxactivityposinf, consdata->maxactivityneginf,
         consdata->maxactivityposhuge, consdata->maxactivityneghuge, absval * maxactbound, FALSE, goodrelax,
         maxresactivity, maxisrelax, ismaxsettoinfinity);
   }
}

/** gets global activity bounds for constraint */
static
void consdataGetGlbActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            glbminactivity,     /**< pointer to store the minimal activity, or NULL, if not needed */
   SCIP_Real*            glbmaxactivity,     /**< pointer to store the maximal activity, or NULL, if not needed */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned minactivity is just a relaxation,
                                              *   i.e. <= the exact minactivity (in case of huge contributions),
                                              *   or equal to the exact minimal activity */
   SCIP_Bool*            maxisrelax,         /**< pointer to store whether the returned maxactivity is just a relaxation,
                                              *   i.e. >= the exact maxactivity (in case of huge contributions),
                                              *   or equal to the exact maximal activity */
   SCIP_Bool*            isminsettoinfinity, /**< pointer to store whether minresactivity was set to infinity or calculated */
   SCIP_Bool*            ismaxsettoinfinity  /**< pointer to store whether maxresactivity was set to infinity or calculated */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert((glbminactivity != NULL && minisrelax != NULL && isminsettoinfinity != NULL)
      || (glbmaxactivity != NULL && maxisrelax != NULL && ismaxsettoinfinity != NULL));

   if( !consdata->validactivities )
   {
      consdataCalcActivities(scip, consdata);
      assert(consdata->validglbminact);
      assert(consdata->validglbmaxact);
   }
   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);
   assert(consdata->glbminactivityneghuge >= 0);
   assert(consdata->glbminactivityposhuge >= 0);
   assert(consdata->glbmaxactivityneghuge >= 0);
   assert(consdata->glbmaxactivityposhuge >= 0);

   if( glbminactivity != NULL )
   {
      assert(isminsettoinfinity != NULL);
      assert(minisrelax != NULL);

      getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
         consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
         glbminactivity, minisrelax, isminsettoinfinity);
   }

   if( glbmaxactivity != NULL )
   {
      assert(ismaxsettoinfinity != NULL);
      assert(maxisrelax != NULL);

      getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
         consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
         glbmaxactivity, maxisrelax, ismaxsettoinfinity);
   }
}

/** gets global activity bounds for constraint after setting variable to zero */
static
void consdataGetGlbActivityResiduals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   SCIP_VAR*             var,                /**< variable to calculate activity residual for */
   SCIP_Real             val,                /**< coefficient value of variable in linear constraint */
   SCIP_Bool             goodrelax,          /**< if we have huge contributions, do we need a good relaxation or are
                                              *   relaxed acticities ignored, anyway? */
   SCIP_Real*            minresactivity,     /**< pointer to store the minimal residual activity, or NULL, if not needed */
   SCIP_Real*            maxresactivity,     /**< pointer to store the maximal residual activity, or NULL, if not needed */
   SCIP_Bool*            minisrelax,         /**< pointer to store whether the returned residual minactivity is just a
                                              *   relaxation, i.e. <= the exact residual minactivity (in case of huge
                                              *   contributions), or equal to the exact residual minactivity */
   SCIP_Bool*            maxisrelax,         /**< pointer to store whether the returned residual maxactivity is just a
                                              *   relaxation, i.e. <= the exact residual maxactivity (in case of huge
                                              *   contributions), or equal to the exact residual minactivity */
   SCIP_Bool*            isminsettoinfinity, /**< pointer to store whether minresactivity was set to infinity or calculated */
   SCIP_Bool*            ismaxsettoinfinity  /**< pointer to store whether maxresactivity was set to infinity or calculated */
   )
{
   SCIP_Real minactbound;
   SCIP_Real maxactbound;
   SCIP_Real absval;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);
   assert((minresactivity != NULL && minisrelax != NULL && isminsettoinfinity != NULL )
      || (maxresactivity != NULL && maxisrelax != NULL && ismaxsettoinfinity != NULL));

   /* get activity bounds of linear constraint */
   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);

   assert(consdata->glbminactivity < SCIP_INVALID);
   assert(consdata->glbmaxactivity < SCIP_INVALID);
   assert(consdata->glbminactivityneginf >= 0);
   assert(consdata->glbminactivityposinf >= 0);
   assert(consdata->glbmaxactivityneginf >= 0);
   assert(consdata->glbmaxactivityposinf >= 0);

   if( val > 0.0 )
   {
      minactbound = SCIPvarGetLbGlobal(var);
      maxactbound = SCIPvarGetUbGlobal(var);
      absval = val;
   }
   else
   {
      minactbound = -SCIPvarGetUbGlobal(var);
      maxactbound = -SCIPvarGetLbGlobal(var);
      absval = -val;
   }

   if( minresactivity != NULL )
   {
      assert(isminsettoinfinity != NULL);
      assert(minisrelax != NULL);

      /* get/compute minactivity by calling getMinActivity() with updated counters for infinite and huge values
       * and contribution of variable set to zero that has to be subtracted from finite part of activity
       */
      if( SCIPisInfinity(scip, minactbound) )
      {
         assert(consdata->glbminactivityposinf >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf - 1, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else if( SCIPisInfinity(scip, -minactbound) )
      {
         assert(consdata->glbminactivityneginf >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf - 1,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, minactbound * absval) )
      {
         assert(consdata->glbminactivityposhuge >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge - 1, consdata->glbminactivityneghuge, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, -minactbound * absval) )
      {
         assert(consdata->glbminactivityneghuge >= 1);

         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge - 1, 0.0, TRUE, goodrelax,
            minresactivity, minisrelax, isminsettoinfinity);
      }
      else
      {
         getMinActivity(scip, consdata, consdata->glbminactivityposinf, consdata->glbminactivityneginf,
            consdata->glbminactivityposhuge, consdata->glbminactivityneghuge, absval * minactbound, TRUE,
            goodrelax, minresactivity, minisrelax, isminsettoinfinity);
      }
   }

   if( maxresactivity != NULL )
   {
      assert(ismaxsettoinfinity != NULL);
      assert(maxisrelax != NULL);

      /* get/compute maxactivity by calling getMaxActivity() with updated counters for infinite and huge values
       * and contribution of variable set to zero that has to be subtracted from finite part of activity
       */
      if( SCIPisInfinity(scip, -maxactbound) )
      {
         assert(consdata->glbmaxactivityneginf >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf - 1,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else if( SCIPisInfinity(scip, maxactbound) )
      {
         assert(consdata->glbmaxactivityposinf >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf - 1, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, absval * maxactbound) )
      {
         assert(consdata->glbmaxactivityposhuge >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge - 1, consdata->glbmaxactivityneghuge, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else if( SCIPisHugeValue(scip, -absval * maxactbound) )
      {
         assert(consdata->glbmaxactivityneghuge >= 1);

         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge - 1, 0.0, TRUE, goodrelax,
            maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
      else
      {
         getMaxActivity(scip, consdata, consdata->glbmaxactivityposinf, consdata->glbmaxactivityneginf,
            consdata->glbmaxactivityposhuge, consdata->glbmaxactivityneghuge, absval * maxactbound, TRUE,
            goodrelax, maxresactivity, maxisrelax, ismaxsettoinfinity);
      }
   }
}

/** calculates the activity of the linear constraint for given solution */
static
SCIP_Real consdataGetActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol                 /**< solution to get activity for, NULL to current solution */
   )
{
   SCIP_Real activity;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( sol == NULL && !SCIPhasCurrentNodeLP(scip) )
      activity = consdataComputePseudoActivity(scip, consdata);
   else
   {
      SCIP_Real solval;
      int nposinf;
      int nneginf;
      SCIP_Bool negsign;
      int v;

      activity = 0.0;
      nposinf = 0;
      nneginf = 0;

      for( v = 0; v < consdata->nvars; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[v]);

         if( consdata->vals[v] < 0 )
            negsign = TRUE;
         else 
            negsign = FALSE;

         if( (SCIPisInfinity(scip, solval) && !negsign) || (SCIPisInfinity(scip, -solval) && negsign) )
            ++nposinf;
         else if( (SCIPisInfinity(scip, solval) && negsign) || (SCIPisInfinity(scip, -solval) && !negsign) )
            ++nneginf;
         else
            activity += consdata->vals[v] * solval;
      }
      assert(nneginf >= 0 && nposinf >= 0);

      SCIPdebugMsg(scip, "activity of linear constraint: %.15g, %d positive infinity values, %d negative infinity values \n", activity, nposinf, nneginf);

      /* check for amount of infinity values and correct the activity */
      if( nposinf > 0 && nneginf > 0 )
         activity = (consdata->rhs + consdata->lhs) / 2;
      else if( nposinf > 0 )
         activity = SCIPinfinity(scip);
      else if( nneginf > 0 )
         activity = -SCIPinfinity(scip);

      SCIPdebugMsg(scip, "corrected activity of linear constraint: %.15g\n", activity);
   }

   if( activity == SCIP_INVALID ) /*lint !e777*/
      return activity;
   else if( activity < 0 )
      activity = MAX(activity, -SCIPinfinity(scip)); /*lint !e666*/
   else
      activity = MIN(activity, SCIPinfinity(scip)); /*lint !e666*/

   return activity;
}

/** calculates the feasibility of the linear constraint for given solution */
static
SCIP_Real consdataGetFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_SOL*             sol                 /**< solution to get feasibility for, NULL to current solution */
   )
{
   SCIP_Real activity;

   assert(scip != NULL);
   assert(consdata != NULL);

   activity = consdataGetActivity(scip, consdata, sol);

   if( activity == SCIP_INVALID ) /*lint !e777*/
      return -SCIPinfinity(scip);

   return MIN(consdata->rhs - activity, activity - consdata->lhs);
}

/** updates bit signatures after adding a single coefficient */
static
void consdataUpdateSignatures(
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   pos                 /**< position of coefficient to update signatures for */
   )
{
   uint64_t varsignature;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real val;

   assert(consdata != NULL);
   assert(consdata->validsignature);

   varsignature = SCIPhashSignature64(SCIPvarGetIndex(consdata->vars[pos]));
   lb = SCIPvarGetLbGlobal(consdata->vars[pos]);
   ub = SCIPvarGetUbGlobal(consdata->vars[pos]);
   val = consdata->vals[pos];
   if( (val > 0.0 && ub > 0.0) || (val < 0.0 && lb < 0.0) )
      consdata->possignature |= varsignature;
   if( (val > 0.0 && lb < 0.0) || (val < 0.0 && ub > 0.0) )
      consdata->negsignature |= varsignature;
}

/** calculates the bit signatures of the given constraint data */
static
void consdataCalcSignatures(
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(consdata != NULL);

   if( !consdata->validsignature )
   {
      int i;

      consdata->validsignature = TRUE;
      consdata->possignature = 0;
      consdata->negsignature = 0;
      for( i = 0; i < consdata->nvars; ++i )
         consdataUpdateSignatures(consdata, i);
   }
}

/** index comparison method of linear constraints: compares two indices of the variable set in the linear constraint */
static
SCIP_DECL_SORTINDCOMP(consdataCompVar)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);

   return SCIPvarCompare(consdata->vars[ind1], consdata->vars[ind2]);
}

/** permutes the constraint's variables according to a given permutation. */
static
void permSortConsdata(
   SCIP_CONSDATA*        consdata,           /**< the constraint data */
   int*                  perm,               /**< the target permutation */
   int                   nvars               /**< the number of variables */
   )
{  /*lint --e{715}*/
   SCIP_VAR* varv;
   SCIP_EVENTDATA* eventdatav;
   SCIP_Real valv;
   int v;
   int i;
   int nexti;

   assert(perm != NULL);
   assert(consdata != NULL);

   /* permute the variables in the linear constraint according to the target permutation */
   eventdatav = NULL;
   for( v = 0; v < nvars; ++v )
   {
      if( perm[v] != v )
      {
         varv = consdata->vars[v];
         valv = consdata->vals[v];
         if( consdata->eventdata != NULL )
            eventdatav = consdata->eventdata[v];
         i = v;
         do
         {
            assert(0 <= perm[i] && perm[i] < nvars);
            assert(perm[i] != i);
            consdata->vars[i] = consdata->vars[perm[i]];
            consdata->vals[i] = consdata->vals[perm[i]];
            if( consdata->eventdata != NULL )
            {
               consdata->eventdata[i] = consdata->eventdata[perm[i]];
               consdata->eventdata[i]->varpos = i;
            }
            nexti = perm[i];
            perm[i] = i;
            i = nexti;
         }
         while( perm[i] != v );
         consdata->vars[i] = varv;
         consdata->vals[i] = valv;
         if( consdata->eventdata != NULL )
         {
            consdata->eventdata[i] = eventdatav;
            consdata->eventdata[i]->varpos = i;
         }
         perm[i] = i;
      }
   }
#ifdef SCIP_DEBUG
   /* check sorting */
   for( v = 0; v < nvars; ++v )
   {
      assert(perm[v] == v);
      assert(consdata->eventdata == NULL || consdata->eventdata[v]->varpos == v);
   }
#endif
}

/** sorts linear constraint's variables depending on the stage of the solving process:
 * - during PRESOLVING
 *       sorts variables by binaries, integers, implicit integers, and continuous variables,
 *       and the variables of the same type by non-decreasing variable index
 *
 * - during SOLVING
 *       sorts binary variables of the remaining problem w.r.t the absolute of their coefficient.
 *       This fastens the propagation time of the constraint handler.
 */
static
SCIP_RETCODE consdataSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< linear constraint data */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   /* check if there are variables for sorting */
   if( consdata->nvars <= 1 )
   {
      consdata->sorted = TRUE;
      consdata->binvarssorted = TRUE;
   }
   else if( SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE && !consdata->sorted )
   {
      int* perm;

      /* get temporary memory to store the sorted permutation */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nvars) );

      /* call sorting method  */
      SCIPsort(perm, consdataCompVar, (void*)consdata, consdata->nvars);

      permSortConsdata(consdata, perm, consdata->nvars);

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &perm);

      consdata->sorted = TRUE;
      consdata->binvarssorted = FALSE;
   }
   else if( SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && !consdata->binvarssorted )
   {
      SCIP_EVENTDATA** eventdata;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      int nvars;
      int v;
      int lastbin;

      nvars = consdata->nvars;
      vars = consdata->vars;
      vals = consdata->vals;
      eventdata = consdata->eventdata;
      assert(vars != NULL || nvars == 0);
      assert(vals != NULL || nvars == 0);

      lastbin = 0;
      /* count binary variables and permute variables such that binaries appear first in the sorted vars array */
      for( v = 0; v < nvars; ++v )
      {
         assert( vars != NULL); /* for flexelint */
         assert( vals != NULL); /* for flexelint */
         if( SCIPvarIsBinary(vars[v]) )
         {
            /* swap variable at the end of the binary variables, if necessary */
            if( lastbin < v )
            {
               SCIP_VAR* tmpvar;
               SCIP_Real tmpval;

               tmpvar = vars[lastbin];
               tmpval = vals[lastbin];

               vars[lastbin] = vars[v];
               vals[lastbin] = vals[v];

               vars[v] = tmpvar;
               vals[v] = tmpval;

               if( eventdata != NULL )
               {
                  SCIP_EVENTDATA* tmpeventdata;

                  tmpeventdata = eventdata[lastbin];
                  eventdata[lastbin] = eventdata[v];
                  eventdata[lastbin]->varpos = lastbin;
                  eventdata[v] = tmpeventdata;
                  eventdata[v]->varpos = v;
               }
               assert(SCIPvarIsBinary(vars[lastbin]));
            }
#ifndef NDEBUG
            else
               assert(lastbin == v);
#endif
            ++lastbin;
         }
      }
      consdata->nbinvars = lastbin;

#ifndef NDEBUG
      /* check sorting */
      for( v = 0; v < nvars; ++v )
      {
         assert(vars != NULL); /* for flexelint */
         assert(eventdata == NULL || eventdata[v]->varpos == v);
         assert((v >= consdata->nbinvars && !SCIPvarIsBinary(vars[v])) || (v < consdata->nbinvars && SCIPvarIsBinary(vars[v])));
      }
#endif

      if( consdata->nbinvars > 1 )
      {
         SCIP_Real* absvals;
         int*       perm;

         assert(lastbin == consdata->nbinvars);
         assert(lastbin <= nvars);
         assert(vals != NULL);

         /* initialize absolute coefficients and the target permutation for binary variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &absvals, lastbin) );
         SCIP_CALL( SCIPallocBufferArray(scip, &perm, lastbin) );

         for( v = 0; v < lastbin; ++v )
         {
            absvals[v] = ABS(vals[v]);
            perm[v] = v;
         }

         /* execute the sorting */
         SCIPsortDownRealInt(absvals, perm, lastbin);

         permSortConsdata(consdata, perm, lastbin);

         /* free temporary arrays */
         SCIPfreeBufferArray(scip, &perm);
         SCIPfreeBufferArray(scip, &absvals);
      }
      consdata->binvarssorted = TRUE;

      /* presolve sorting cannot be guaranteed after binary sorting */
      consdata->sorted = (consdata->sorted && consdata->nbinvars == 0);
   }
   assert(SCIPgetStage(scip) < SCIP_STAGE_INITSOLVE || consdata->binvarssorted);
   assert(SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE || consdata->sorted);

   return SCIP_OKAY;
}


/*
 * local linear constraint handler methods
 */

/** sets left hand side of linear constraint */
static
SCIP_RETCODE chgLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, lhs));

   /* adjust value to not be smaller than -inf */
   if ( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->lhs, lhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */   
   if( SCIPisEQ(scip, lhs, consdata->rhs) )
   {
      consdata->rhs = lhs;
      assert(consdata->row == NULL);
   }

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the left hand side switched from -infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }
      }
      else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -lhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the left hand side switched from a non-infinite value to -infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
         }
      }
   }

   /* check whether the left hand side is increased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( !SCIPisInfinity(scip, -lhs) && SCIPisGT(scip, lhs, consdata->lhs) )
   {
      consdata->boundstightened = 0;
      consdata->presolved = FALSE;
      consdata->cliquesadded = FALSE;
      consdata->implsadded = FALSE;

      /* mark the constraint for propagation */
      if( SCIPconsIsTransformed(cons) )
      {
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
      }
   }

   /* set new left hand side and update constraint data */
   consdata->lhs = lhs;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->rangedrowpropagated = 0;


   /* update the lhs of the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, consdata->row, lhs) );
   }

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, -rhs));

   /* adjust value to not be larger than inf */
   if ( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || (consdata->vars != NULL && consdata->vals != NULL));
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->rhs, rhs) )
      return SCIP_OKAY;

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */   
   if( SCIPisEQ(scip, rhs, consdata->lhs) )
   {
      consdata->lhs = rhs;
      assert(consdata->row == NULL);
   }

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      if( SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }
      }
      else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, rhs) )
      {
         SCIP_VAR** vars;
         SCIP_Real* vals;
         int v;

         /* the right hand side switched from a non-infinite value to infinity -> remove rounding locks */
         vars = consdata->vars;
         vals = consdata->vals;

         for( v = 0; v < consdata->nvars; ++v )
         {
            assert(vars[v] != NULL);
            assert(!SCIPisZero(scip, vals[v]));

            if( SCIPisPositive(scip, vals[v]) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, FALSE, TRUE) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, vars[v], cons, TRUE, FALSE) );
            }
         }
      }
   }

   /* check whether the right hand side is decreased, if and only if that's the case we maybe can propagate, tighten and add more cliques */
   if( !SCIPisInfinity(scip, rhs) && SCIPisLT(scip, rhs, consdata->rhs) )
   {
      consdata->boundstightened = 0;
      consdata->presolved = FALSE;
      consdata->cliquesadded = FALSE;
      consdata->implsadded = FALSE;

      /* mark the constraint for propagation */
      if( SCIPconsIsTransformed(cons) )
      {
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
      }
   }

   /* set new right hand side and update constraint data */
   consdata->rhs = rhs;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->rangedrowpropagated = 0;

   /* update the rhs of the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, consdata->row, rhs) );
   }

   return SCIP_OKAY;
}

/** adds coefficient in linear constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   /* ignore coefficient if it is nearly zero */
   if( SCIPisZero(scip, val) )
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

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   consdata->vals[consdata->nvars] = val;
   consdata->nvars++;

   /* capture variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   /* if we are in transformed problem, the variable needs an additional event data */
   if( transformed )
   {
      if( consdata->eventdata != NULL )
      {
         SCIP_CONSHDLR* conshdlr;
         SCIP_CONSHDLRDATA* conshdlrdata;

         /* check for event handler */
         conshdlr = SCIPconsGetHdlr(cons);
         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         assert(conshdlrdata->eventhdlr != NULL);

         /* initialize eventdata array */
         consdata->eventdata[consdata->nvars-1] = NULL;

         /* catch bound change events of variable */
         SCIP_CALL( consCatchEvent(scip, cons, conshdlrdata->eventhdlr, consdata->nvars-1) );
      }

      /* update minimum and maximum activities */
      consdataUpdateAddCoef(scip, consdata, var, val, FALSE);

      /* update maximum activity delta */
      if( !SCIPisInfinity(scip, consdata->maxactdelta ) )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);

         if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
         {
            consdata->maxactdelta = SCIPinfinity(scip);
            consdata->maxactdeltavar = var;
         }
         else
         {
            SCIP_Real domain = ub - lb;
            SCIP_Real delta = REALABS(val) * domain;

            if( delta > consdata->maxactdelta )
            {
               consdata->maxactdelta = delta;
               consdata->maxactdeltavar = var;
            }
         }
      }
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockRounding(scip, cons, var, val) );

   /* mark the constraint for propagation */
   if( transformed )
   {
      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->boundstightened = 0;
   consdata->presolved = FALSE;
   consdata->removedfixings = consdata->removedfixings && SCIPvarIsActive(var);

   if( consdata->validsignature )
      consdataUpdateSignatures(consdata, consdata->nvars-1);

   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   consdata->implsadded = FALSE;
   consdata->rangedrowpropagated = 0;

   if( consdata->nvars == 1 )
   {
      consdata->binvarssorted = TRUE;
      consdata->sorted = TRUE;
      consdata->merged = TRUE;
   }
   else
   {
      consdata->binvarssorted = consdata->binvarssorted && !SCIPvarIsBinary(var);
      consdata->sorted = consdata->sorted
         && (SCIPvarCompare(consdata->vars[consdata->nvars-2], consdata->vars[consdata->nvars-1]) <= 0);
      consdata->merged = FALSE;
   }

   /* update hascontvar and hasnonbinvar flags */
   if( consdata->hasnonbinvalid && !consdata->hascontvar )
   {
      SCIP_VARTYPE vartype = SCIPvarGetType(var);

      if( vartype != SCIP_VARTYPE_BINARY )
      {
         consdata->hasnonbinvar = TRUE;

         if( vartype == SCIP_VARTYPE_CONTINUOUS )
            consdata->hascontvar = TRUE;
      }
   }

   /* add the new coefficient to the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, val) );
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from linear constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, var, val) );

   /* if we are in transformed problem, delete the event data of the variable */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      /* check for event handler */
      conshdlr = SCIPconsGetHdlr(cons);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* drop bound change events of variable */
      if( consdata->eventdata != NULL )
      {
         SCIP_CALL( consDropEvent(scip, cons, conshdlrdata->eventhdlr, pos) );
         assert(consdata->eventdata[pos] == NULL);
      }
   }

   /* move the last variable to the free slot */
   if( pos != consdata->nvars-1 )
   {
      consdata->binvarssorted = consdata->binvarssorted && !SCIPvarIsBinary(consdata->vars[pos]);

      consdata->vars[pos] = consdata->vars[consdata->nvars-1];
      consdata->vals[pos] = consdata->vals[consdata->nvars-1];

      if( consdata->eventdata != NULL )
      {
         consdata->eventdata[pos] = consdata->eventdata[consdata->nvars-1];
         assert(consdata->eventdata[pos] != NULL);
         consdata->eventdata[pos]->varpos = pos;
      }
      consdata->sorted = consdata->sorted && (pos + 2 >= consdata->nvars || (SCIPvarCompare(consdata->vars[pos], consdata->vars[pos + 1]) <= 0));
   }
   consdata->nvars--;

   /* if at most one variable is left, the activities should be recalculated (to correspond exactly to the bounds
    * of the remaining variable, or give exactly 0.0)
    */
   if( consdata->nvars <= 1 )
      consdataInvalidateActivities(consdata);
   else
   {
      if( SCIPconsIsTransformed(cons) )
      {
         /* if we are in transformed problem, update minimum and maximum activities */
         consdataUpdateDelCoef(scip, consdata, var, val, TRUE);

         /* if the variable defining the maximal activity delta was removed from the constraint, the maximal activity
          * delta needs to be recalculated on the next real propagation
          */
         if( consdata->maxactdeltavar == var )
         {
            consdata->maxactdelta = SCIP_INVALID;
            consdata->maxactdeltavar = NULL;
         }
      }
   }

   /* mark the constraint for propagation */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->boundstightened = 0;
   consdata->presolved = FALSE;
   consdata->validsignature = FALSE;
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   consdata->implsadded = FALSE;
   consdata->rangedrowpropagated = 0;

   /* check if hasnonbinvar flag might be incorrect now */
   if( consdata->hasnonbinvar && SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
   {
      consdata->hasnonbinvalid = FALSE;
   }

   /* delete coefficient from the LP row */
   if( consdata->row != NULL )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, var, -val) );
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   return SCIP_OKAY;
}

/** changes coefficient value at given position of linear constraint data */
static
SCIP_RETCODE chgCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< position of coefficient to delete */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisZero(scip, newval));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(!SCIPisZero(scip, newval));

   var = consdata->vars[pos];
   val = consdata->vals[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* if necessary, update the rounding locks of the variable */
   if( SCIPconsIsLocked(cons) && newval * val < 0.0 )
   {
      assert(SCIPconsIsTransformed(cons));

      /* remove rounding locks for variable with old coefficient */
      SCIP_CALL( unlockRounding(scip, cons, var, val) );

      /* install rounding locks for variable with new coefficient */
      SCIP_CALL( lockRounding(scip, cons, var, newval) );
   }

   /* change the value */
   consdata->vals[pos] = newval;

   consdata->binvarssorted = consdata->binvarssorted && !SCIPvarIsBinary(var);

   /* update minimum and maximum activities */
   if( SCIPconsIsTransformed(cons) )
      consdataUpdateChgCoef(scip, consdata, var, val, newval, TRUE);

   /* mark the constraint for propagation */
   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   }

   consdata->boundstightened = 0;
   consdata->presolved = FALSE;
   consdata->validsignature = consdata->validsignature && (newval * val > 0.0);
   consdata->changed = TRUE;
   consdata->normalized = FALSE;
   consdata->upgradetried = FALSE;
   consdata->cliquesadded = FALSE;
   consdata->implsadded = FALSE;
   consdata->rangedrowpropagated = 0;

   return SCIP_OKAY;
}

/** scales a linear constraint with a constant scalar */
static
SCIP_RETCODE scaleCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint to scale */
   SCIP_Real             scalar              /**< value to scale constraint with */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real newval;
   SCIP_Real absscalar;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);
   assert(!SCIPisEQ(scip, scalar, 1.0));

   /* scale the coefficients */
   for( i = consdata->nvars - 1; i >= 0; --i )
   {
      newval = scalar * consdata->vals[i];

      /* because SCIPisScalingIntegral uses another integrality check as SCIPfeasFloor, we add an additional 0.5 before
       * flooring down our new value
       */
      if( SCIPisScalingIntegral(scip, consdata->vals[i], scalar) )
         newval = SCIPfeasFloor(scip, newval + 0.5);

      if( SCIPisZero(scip, newval) )
      {
         SCIPwarningMessage(scip, "coefficient %.15g of variable <%s> in linear constraint <%s> scaled to zero (scalar: %.15g)\n",
            consdata->vals[i], SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons), scalar);
         SCIP_CALL( delCoefPos(scip, cons, i) );
      }
      else
         consdata->vals[i] = newval;
   }

   /* scale the sides */
   if( scalar < 0.0 )
   {
      SCIP_Real lhs;

      lhs = consdata->lhs;
      consdata->lhs = -consdata->rhs;
      consdata->rhs = -lhs;
   }
   absscalar = REALABS(scalar);
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      newval = absscalar * consdata->lhs;

      /* because SCIPisScalingIntegral uses another integrality check as SCIPfeasFloor, we add an additional 0.5 before
       * flooring down our new value
       */
      if( SCIPisScalingIntegral(scip, consdata->lhs, absscalar) )
         consdata->lhs = SCIPfeasFloor(scip, newval + 0.5);
      else
         consdata->lhs = newval;
   }
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      newval = absscalar * consdata->rhs;

      /* because SCIPisScalingIntegral uses another integrality check as SCIPfeasCeil, we subtract 0.5 before ceiling up
       * our new value
       */
      if( SCIPisScalingIntegral(scip, consdata->rhs, absscalar) )
         consdata->rhs = SCIPfeasCeil(scip, newval - 0.5);
      else
         consdata->rhs = newval;
   }

   consdataInvalidateActivities(consdata);
   consdata->cliquesadded = FALSE;
   consdata->implsadded = FALSE;

   return SCIP_OKAY;
}

/* perform deletion of variables in all constraints of the constraint handler */
static
SCIP_RETCODE performVarDeletions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int v;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* iterate over all constraints */
   for( i = 0; i < nconss; i++ )
   {
      consdata = SCIPconsGetData(conss[i]);

      /* constraint is marked, that some of its variables were deleted */
      if( consdata->varsdeleted )
      {
         /* iterate over all variables of the constraint and delete them from the constraint */
         for( v = consdata->nvars - 1; v >= 0; --v )
         {
            if( SCIPvarIsDeleted(consdata->vars[v]) )
            {
               SCIP_CALL( delCoefPos(scip, conss[i], v) );
            }
         }
         consdata->varsdeleted = FALSE;
      }
   }

   return SCIP_OKAY;
}


/** normalizes a linear constraint with the following rules:
 *  - if all coefficients have them same absolute value, change them to (-)1.0
 *  - multiplication with +1 or -1:
 *      Apply the following rules in the given order, until the sign of the factor is determined. Later rules only apply,
 *      if the current rule doesn't determine the sign):
 *        1. the right hand side must not be negative
 *        2. the right hand side must not be infinite
 *        3. the absolute value of the right hand side must be greater than that of the left hand side
 *        4. the number of positive coefficients must not be smaller than the number of negative coefficients
 *        5. multiply with +1
 *  - rationals to integrals
 *      Try to identify a rational representation of the fractional coefficients, and multiply all coefficients
 *      by the smallest common multiple of all denominators to get integral coefficients.
 *      Forbid large denominators due to numerical stability.
 *  - division by greatest common divisor
 *      If all coefficients are integral, divide them by the greatest common divisor.
 */
static
SCIP_RETCODE normalizeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint to normalize */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Longint scm;
   SCIP_Longint nominator;
   SCIP_Longint denominator;
   SCIP_Longint gcd;
   SCIP_Longint maxmult;
   SCIP_Real epsilon;
   SCIP_Real feastol;
   SCIP_Real maxabsval;
   SCIP_Real minabsval;
   SCIP_Bool success;
   SCIP_Bool onlyintegral;
   int nvars;
   int mult;
   int nposcoeffs;
   int nnegcoeffs;
   int i;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   /* we must not change a modifiable constraint in any way */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* get constraint data */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint is already normalized */
   if( consdata->normalized )
      return SCIP_OKAY;

   /* get coefficient arrays */
   vals = consdata->vals;
   nvars = consdata->nvars;
   vars = consdata->vars;
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   if( nvars == 0 )
   {
      consdata->normalized = TRUE;
      return SCIP_OKAY;
   }

   assert(vars != NULL);
   assert(vals != NULL);

   /* get maximal and minimal absolute coefficient */
   maxabsval = consdataGetMaxAbsval(consdata);
   minabsval = consdataGetMinAbsval(consdata);

   /* return if scaling by maxval will eliminate coefficients */
   if( SCIPisZero(scip, minabsval/maxabsval) )
      return SCIP_OKAY;

   /* check if all coefficients are in absolute value equal, and not 1.0 */
   if( !SCIPisEQ(scip, maxabsval, 1.0) )
   {
      SCIP_Bool abscoefsequ;

      abscoefsequ = TRUE;

      for( v = nvars - 1; v >= 0; --v )
      {
         if( !SCIPisEQ(scip, REALABS(vals[v]), maxabsval) )
         {
            abscoefsequ = FALSE;
            break;
         }
      }

      /* all coefficients are in absolute value equal, so change them to (-)1.0 */
      if( abscoefsequ )
      {
         SCIPdebugMsg(scip, "divide linear constraint with %g, because all coefficents are in absolute value the same\n", maxabsval);
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIP_CALL( scaleCons(scip, cons, 1/maxabsval) );

         if( consdata->validmaxabsval )
         {
            if( !SCIPisEQ(scip, consdata->maxabsval, 1.0) )
               consdata->maxabsval = 1.0;
            if( !SCIPisEQ(scip, consdata->minabsval, 1.0) )
               consdata->minabsval = 1.0;

            maxabsval = 1.0;
         }
         else
         {
            /* get maximal absolute coefficient */
            maxabsval = consdataGetMaxAbsval(consdata);
         }

         /* get new consdata information, because scalecons() might have deleted variables */
         vals = consdata->vals;
         nvars = consdata->nvars;
         vars = consdata->vars;

         assert(nvars == 0 || vars != NULL);
         assert(nvars == 0 || vals != NULL);
      }
   }

   /* nvars might have changed */
   if( nvars == 0 )
   {
      consdata->normalized = TRUE;
      return SCIP_OKAY;
   }

   assert(vars != NULL);
   assert(vals != NULL);

   /* calculate the maximal multiplier for common divisor calculation:
    *   |p/q - val| < epsilon  and  q < feastol/epsilon  =>  |p - q*val| < feastol
    * which means, a value of feastol/epsilon should be used as maximal multiplier;
    * additionally, we don't want to scale the constraint if this would lead to too
    * large coefficients
    */
   epsilon = SCIPepsilon(scip) * 0.9;  /* slightly decrease epsilon to be safe in rational conversion below */
   feastol = SCIPfeastol(scip);
   maxmult = (SCIP_Longint)(feastol/epsilon + feastol);
   maxmult = MIN(maxmult, (SCIP_Longint)( MAXSCALEDCOEF/MAX(maxabsval, 1.0)));

   if( !consdata->hasnonbinvalid )
      consdataCheckNonbinvar(consdata);

   /* if all variables are of integral type we will allow a greater multiplier */
   if( !consdata->hascontvar )
   {
      if( SCIPvarGetType(vars[nvars - 1]) != SCIP_VARTYPE_CONTINUOUS )
      {
	 maxmult = (SCIP_Longint) (MAXSCALEDCOEFINTEGER/(MAX(maxabsval, 1.0)));
      }
   }
   else
   {
      SCIP_Bool foundcont;

      foundcont = FALSE;

      for( v = nvars - 1; v >= 0; --v )
      {
	 if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
	 {
	    foundcont = TRUE;
	    break;
	 }
      }

      if( !foundcont )
      {
	 maxmult = (SCIP_Longint) (MAXSCALEDCOEFINTEGER/(MAX(maxabsval, 1.0)));
      }
   }

   /*
    * multiplication with +1 or -1
    */
   mult = 0;

   /* 1. the right hand side must not be negative */
   if( SCIPisPositive(scip, consdata->lhs) )
      mult = +1;
   else if( SCIPisNegative(scip, consdata->rhs) )
      mult = -1;

   if( mult == 0 )
   {
      /* 2. the right hand side must not be infinite */
      if( SCIPisInfinity(scip, -consdata->lhs) )
         mult = +1;
      else if( SCIPisInfinity(scip, consdata->rhs) )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 3. the absolute value of the right hand side must be greater than that of the left hand side */
      if( SCIPisGT(scip, REALABS(consdata->rhs), REALABS(consdata->lhs)) )
         mult = +1;
      else if( SCIPisLT(scip, REALABS(consdata->rhs), REALABS(consdata->lhs)) )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 4. the number of positive coefficients must not be smaller than the number of negative coefficients */
      nposcoeffs = 0;
      nnegcoeffs = 0;
      for( i = 0; i < nvars; ++i )
      {
         if( vals[i] > 0.0 )
            nposcoeffs++;
         else
            nnegcoeffs++;
      }
      if( nposcoeffs > nnegcoeffs )
         mult = +1;
      else if( nposcoeffs < nnegcoeffs )
         mult = -1;
   }

   if( mult == 0 )
   {
      /* 5. multiply with +1 */
      mult = +1;
   }

   assert(mult == +1 || mult == -1);
   if( mult == -1 )
   {
      /* scale the constraint with -1 */
      SCIPdebugMsg(scip, "multiply linear constraint with -1.0\n");
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( scaleCons(scip, cons, -1.0) );

      /* scalecons() can delete variables, but scaling with -1 should not do that */
      assert(nvars == consdata->nvars);
   }

   /*
    * rationals to integrals
    *
    * @todo try scaling only on behalf of non-continuous variables
    */
   success = TRUE;
   scm = 1;
   for( i = 0; i < nvars && success && scm <= maxmult; ++i )
   {
      if( !SCIPisIntegral(scip, vals[i]) )
      {
         /* epsilon has been slightly decreased above - to be on the safe side */
         success = SCIPrealToRational(vals[i], -epsilon, epsilon , maxmult, &nominator, &denominator);
         if( success )
            scm = SCIPcalcSmaComMul(scm, denominator);
      }
   }
   assert(scm >= 1);

   /* it might be that we have really big coefficients, but all are integral, in that case we want to divide them by
    * their greatest common divisor
    */
   onlyintegral = TRUE;
   if( scm == 1 )
   {
      for( i = nvars - 1; i >= 0; --i )
      {
         if( !SCIPisIntegral(scip, vals[i]) )
         {
            onlyintegral = FALSE;
            break;
         }
      }
   }

   success = success && (scm <= maxmult || (scm == 1 && onlyintegral));
   if( success && scm != 1 )
   {
      /* scale the constraint with the smallest common multiple of all denominators */
      SCIPdebugMsg(scip, "scale linear constraint with %" SCIP_LONGINT_FORMAT " to make coefficients integral\n", scm);
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIP_CALL( scaleCons(scip, cons, (SCIP_Real)scm) );

      if( consdata->validmaxabsval )
      {
         consdata->maxabsval *= REALABS((SCIP_Real)scm);
         if( !SCIPisIntegral(scip, consdata->maxabsval) )
         {
            consdata->validmaxabsval = FALSE;
            consdata->maxabsval = SCIP_INVALID;
            consdataCalcMaxAbsval(consdata);
         }
      }

      if( consdata->validminabsval )
      {
         consdata->minabsval *= REALABS((SCIP_Real)scm);
         if( !SCIPisIntegral(scip, consdata->minabsval) )
         {
            consdata->validminabsval = FALSE;
            consdata->minabsval = SCIP_INVALID;
            consdataCalcMinAbsval(consdata);
         }
      }

      /* get new consdata information, because scalecons() might have deleted variables */
      vals = consdata->vals;
      nvars = consdata->nvars;
      assert(nvars == 0 || vals != NULL);
   }

   /*
    * division by greatest common divisor
    */
   if( success && nvars >= 1 )
   {
      /* all coefficients are integral: divide them by their greatest common divisor */
      assert(SCIPisIntegral(scip, vals[0]));

      gcd = (SCIP_Longint)(REALABS(vals[0]) + feastol);
      for( i = 1; i < nvars && gcd > 1; ++i )
      {
         assert(SCIPisIntegral(scip, vals[i]));
         gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[i]) + feastol));
      }

      if( gcd > 1 )
      {
         /* divide the constraint by the greatest common divisor of the coefficients */
         SCIPdebugMsg(scip, "divide linear constraint by greatest common divisor %" SCIP_LONGINT_FORMAT "\n", gcd);
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIP_CALL( scaleCons(scip, cons, 1.0/(SCIP_Real)gcd) );

         if( consdata->validmaxabsval )
         {
            consdata->maxabsval /= REALABS((SCIP_Real)gcd);
         }
         if( consdata->validminabsval )
         {
            consdata->minabsval /= REALABS((SCIP_Real)gcd);
         }
      }
   }

   /* mark constraint to be normalized */
   consdata->normalized = TRUE;

   SCIPdebugMsg(scip, "normalized constraint:\n");
   SCIPdebugPrintCons(scip, cons, NULL);

   return SCIP_OKAY;
}

/** replaces multiple occurrences of a variable by a single coefficient */
static
SCIP_RETCODE mergeMultiples(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real valsum;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->merged )
      return SCIP_OKAY;

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata) );

   /* go backwards through the constraint looking for multiple occurrences of the same variable;
    * backward direction is necessary, since delCoefPos() modifies the given position and
    * the subsequent ones
    */
   v = consdata->nvars-1;
   while( v >= 1 )
   {
      var = consdata->vars[v];
      if( consdata->vars[v-1] == var )
      {
         valsum = consdata->vals[v];
         do
         {
            SCIP_CALL( delCoefPos(scip, cons, v) );
            --v;
            valsum += consdata->vals[v];
         }
         while( v >= 1 && consdata->vars[v-1] == var );

         /* modify the last existing occurrence of the variable */
         assert(consdata->vars[v] == var);
         if( SCIPisZero(scip, valsum) )
         {
            SCIP_CALL( delCoefPos(scip, cons, v) );

            /* if the variable defining the maximal activity delta was removed from the constraint, the maximal activity
             * delta needs to be recalculated on the next real propagation
             */
            if( consdata->maxactdeltavar == var )
            {
               consdata->maxactdelta = SCIP_INVALID;
               consdata->maxactdeltavar = NULL;
            }
         }
         else
         {
            SCIP_CALL( chgCoefPos(scip, cons, v, valsum) );
         }
      }
      --v;
   }

   consdata->merged = TRUE;

   return SCIP_OKAY;
}

/** replaces all fixed and aggregated variables by their non-fixed counterparts */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility is detected; or NULL if this
                                              *   information is not needed; in this case, we apply all fixings
                                              *   instead of stopping after the first infeasible one */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_VAR** aggrvars;
   SCIP_Real val;
   SCIP_Real* aggrscalars;
   SCIP_Real fixedval;
   SCIP_Real aggrconst;
   int v;
   int naggrvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   if( infeasible != NULL )
      *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !consdata->removedfixings )
   {
      SCIP_Real lhssubtrahend;
      SCIP_Real rhssubtrahend;

      lhssubtrahend = 0.0;
      rhssubtrahend = 0.0;

      SCIPdebugMsg(scip, "applying fixings:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      v = 0;
      while( v < consdata->nvars )
      {
         var = consdata->vars[v];
         val = consdata->vals[v];
         assert(SCIPvarIsTransformed(var));

         switch( SCIPvarGetStatus(var) )
         {
         case SCIP_VARSTATUS_ORIGINAL:
            SCIPerrorMessage("original variable in transformed linear constraint\n");
            return SCIP_INVALIDDATA;

         case SCIP_VARSTATUS_LOOSE:
         case SCIP_VARSTATUS_COLUMN:
            ++v;
            break;

         case SCIP_VARSTATUS_FIXED:
            assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)));
            fixedval = SCIPvarGetLbGlobal(var);
            if( !SCIPisInfinity(scip, -consdata->lhs) )
            {
               if( SCIPisInfinity(scip, ABS(fixedval)) )
               {
                  if( val * fixedval > 0.0 )
                  {
                     SCIP_CALL( chgLhs(scip, cons, -SCIPinfinity(scip)) );
                  }
                  else
                  {
                     if( infeasible != NULL )
                     {
                        /* if lhs gets infinity it means that the problem is infeasible */
                        *infeasible = TRUE;
                        return SCIP_OKAY;
                     }
                     else
                     {
                        SCIP_CALL( chgLhs(scip, cons, SCIPinfinity(scip)) );
                     }
                  }
               }
               else
                  lhssubtrahend += val * fixedval;
            }
            if( !SCIPisInfinity(scip, consdata->rhs) )
            {
               if( SCIPisInfinity(scip, ABS(fixedval)) )
               {
                  if( val * fixedval > 0.0 )
                  {
                     if( infeasible != NULL )
                     {
                        /* if rhs gets -infinity it means that the problem is infeasible */
                        *infeasible = TRUE;
                        return SCIP_OKAY;
                     }
                     else
                     {
                        SCIP_CALL( chgRhs(scip, cons, -SCIPinfinity(scip)) );
                     }
                  }
                  else
                  {
                     SCIP_CALL( chgRhs(scip, cons, SCIPinfinity(scip)) );
                  }
               }
               else
                  rhssubtrahend += val * fixedval;
            }
            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_AGGREGATED:
	 {
	    SCIP_VAR* activevar = SCIPvarGetAggrVar(var);
	    SCIP_Real activescalar = val * SCIPvarGetAggrScalar(var);
	    SCIP_Real activeconstant = val * SCIPvarGetAggrConstant(var);

	    assert(activevar != NULL);
	    SCIP_CALL( SCIPgetProbvarSum(scip, &activevar, &activescalar, &activeconstant) );
	    assert(activevar != NULL);

	    if( !SCIPisZero(scip, activescalar) )
	    {
	       SCIP_CALL( addCoef(scip, cons, activevar, activescalar) );
	    }

	    if( !SCIPisZero(scip, activeconstant) )
	    {
	       if( !SCIPisInfinity(scip, -consdata->lhs) )
		  lhssubtrahend += activeconstant;
	       if( !SCIPisInfinity(scip, consdata->rhs) )
		  rhssubtrahend += activeconstant;
	    }

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;
	 }
         case SCIP_VARSTATUS_MULTAGGR:
            SCIP_CALL( SCIPflattenVarAggregationGraph(scip, var) );
            naggrvars = SCIPvarGetMultaggrNVars(var);
            aggrvars = SCIPvarGetMultaggrVars(var);
            aggrscalars = SCIPvarGetMultaggrScalars(var);
            for( i = 0; i < naggrvars; ++i )
            {
               SCIP_CALL( addCoef(scip, cons, aggrvars[i], val * aggrscalars[i]) );
            }
            aggrconst = SCIPvarGetMultaggrConstant(var);

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * aggrconst;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val * aggrconst;

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         case SCIP_VARSTATUS_NEGATED:
            SCIP_CALL( addCoef(scip, cons, SCIPvarGetNegationVar(var), -val) );
            aggrconst = SCIPvarGetNegationConstant(var);

            if( !SCIPisInfinity(scip, -consdata->lhs) )
               lhssubtrahend += val * aggrconst;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               rhssubtrahend += val * aggrconst;

            SCIP_CALL( delCoefPos(scip, cons, v) );
            break;

         default:
            SCIPerrorMessage("unknown variable status\n");
            SCIPABORT();
            return SCIP_INVALIDDATA;  /*lint !e527*/
         }
      }

      if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->lhs) )
      {
         /* for large numbers that are relatively equal, substraction can lead to cancellation,
          * causing wrong fixings of other variables --> better use a real zero here;
          * for small numbers, polishing the difference might lead to wrong results -->
          * better use the exact difference in this case
          */
         if( SCIPisEQ(scip, lhssubtrahend, consdata->lhs) && SCIPisFeasGE(scip, REALABS(lhssubtrahend), 1.0) )
         {
            SCIP_CALL( chgLhs(scip, cons, 0.0) );
         }
         else
         {
            SCIP_CALL( chgLhs(scip, cons, consdata->lhs - lhssubtrahend) );
         }
      }
      if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->rhs))
      {

         /* for large numbers that are relatively equal, substraction can lead to cancellation,
          * causing wrong fixings of other variables --> better use a real zero here;
          * for small numbers, polishing the difference might lead to wrong results -->
          * better use the exact difference in this case
          */
         if( SCIPisEQ(scip, rhssubtrahend, consdata->rhs ) && SCIPisFeasGE(scip, REALABS(rhssubtrahend), 1.0) )
         {
            SCIP_CALL( chgRhs(scip, cons, 0.0) );
         }
         else
         {
            SCIP_CALL( chgRhs(scip, cons, consdata->rhs - rhssubtrahend) );
         }
      }

      consdata->removedfixings = TRUE;

      SCIPdebugMsg(scip, "after fixings:\n");
      SCIPdebugPrintCons(scip, cons, NULL);

      /* if aggregated variables have been replaced, multiple entries of the same variable are possible and we have
       * to clean up the constraint
       */
      SCIP_CALL( mergeMultiples(scip, cons) );

      SCIPdebugMsg(scip, "after merging:\n");
      SCIPdebugPrintCons(scip, cons, NULL);
   }
   assert(consdata->removedfixings);

#ifndef NDEBUG
   /* check, if all fixings are applied */
   for( v = 0; v < consdata->nvars; ++v )
      assert(SCIPvarIsActive(consdata->vars[v]));
#endif

   return SCIP_OKAY;
}

/** for each variable in the linear constraint, except the inferred variable, adds one bound to the conflict analysis'
 *  candidate store (bound depends on sign of coefficient and whether the left or right hand side was the reason for the
 *  inference variable's bound change); the conflict analysis can be initialized with the linear constraint being the
 *  conflict detecting constraint by using NULL as inferred variable
 */
static
SCIP_RETCODE addConflictBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   int                   inferpos,           /**< position of the inferred variable in the vars array */
   SCIP_Bool             reasonisrhs         /**< is the right hand side responsible for the bound change? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   assert(consdata != NULL);

   vars = consdata->vars;
   vals = consdata->vals;
   nvars = consdata->nvars;

   assert(vars != NULL || nvars == 0);
   assert(vals != NULL || nvars == 0);

   assert(-1 <= inferpos && inferpos < nvars);
   assert((infervar == NULL) == (inferpos == -1));
   assert(inferpos == -1 || vars[inferpos] == infervar); /*lint !e613*/

   /* for each variable, add the bound to the conflict queue, that is responsible for the minimal or maximal
    * residual value, depending on whether the left or right hand side is responsible for the bound change:
    *  - if the right hand side is the reason, the minimal residual activity is responsible
    *  - if the left hand side is the reason, the maximal residual activity is responsible
    */

   /* if the variable is integral we only need to add reason bounds until the propagation could be applied */
   if( infervar == NULL || SCIPvarIsIntegral(infervar) )
   {
      SCIP_Real minresactivity;
      SCIP_Real maxresactivity;
      SCIP_Bool minisrelax;
      SCIP_Bool maxisrelax;
      SCIP_Bool isminsettoinfinity;
      SCIP_Bool ismaxsettoinfinity;

      minresactivity = -SCIPinfinity(scip);
      maxresactivity = SCIPinfinity(scip);

      /* calculate the minimal and maximal global activity of all other variables involved in the constraint */
      if( infervar != NULL )
      {
         assert(vals != NULL); /* for flexelint */
         if( reasonisrhs )
            consdataGetGlbActivityResiduals(scip, consdata, infervar, vals[inferpos], FALSE, &minresactivity, NULL,
               &minisrelax, NULL, &isminsettoinfinity, NULL);
         else
            consdataGetGlbActivityResiduals(scip, consdata, infervar, vals[inferpos], FALSE, NULL, &maxresactivity,
               NULL, &maxisrelax, NULL, &ismaxsettoinfinity);
      }
      else
      {
         if( reasonisrhs )
            consdataGetGlbActivityBounds(scip, consdata, FALSE, &minresactivity, NULL,
               &minisrelax, NULL, &isminsettoinfinity, NULL);
         else
            consdataGetGlbActivityBounds(scip, consdata, FALSE, NULL, &maxresactivity,
               NULL, &maxisrelax, NULL, &ismaxsettoinfinity);
      }

      /* we can only do something clever, if the residual activity is finite and not relaxed */
      if( (reasonisrhs && !isminsettoinfinity && !minisrelax) || (!reasonisrhs && !ismaxsettoinfinity && !maxisrelax) ) /*lint !e644*/
      {
         SCIP_Real rescap;
         SCIP_Bool resactisinf;

         resactisinf = FALSE;

         /* calculate the residual capacity that would be left, if the variable would be set to one more / one less
          * than its inferred bound
          */
         if( infervar != NULL )
         {
            assert(vals != NULL); /* for flexelint */

            if( reasonisrhs )
            {
               if( SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastglbminactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, infervar, &minresactivity, TRUE, TRUE);
                  if( SCIPisInfinity(scip, -minresactivity) )
                     resactisinf = TRUE;
               }
               rescap = consdata->rhs - minresactivity;
            }
            else
            {
               if( SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastglbmaxactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, infervar, &maxresactivity, FALSE, TRUE);
                  if( SCIPisInfinity(scip, maxresactivity) )
                     resactisinf = TRUE;
               }
               rescap = consdata->lhs - maxresactivity;
            }

            if( reasonisrhs == (vals[inferpos] > 0.0) )
               rescap -= vals[inferpos] * (SCIPgetVarUbAtIndex(scip, infervar, bdchgidx, TRUE) + 1.0);
            else
               rescap -= vals[inferpos] * (SCIPgetVarLbAtIndex(scip, infervar, bdchgidx, TRUE) - 1.0);
         }
         else
            rescap = (reasonisrhs ? consdata->rhs - minresactivity : consdata->lhs - maxresactivity);

         if( !resactisinf )
         {
            /* now add bounds as reasons until the residual capacity is exceeded */
            for( i = 0; i < nvars; ++i )
            {
               assert( vars != NULL && vals != NULL ); /* for lint */

               /* zero coefficients and the infered variable can be ignored */
               if( vars[i] == infervar || SCIPisZero(scip, vals[i]) )
                  continue;

               /* check if the residual capacity is exceeded */
               if( (reasonisrhs && SCIPisFeasNegative(scip, rescap))
                  || (!reasonisrhs && SCIPisFeasPositive(scip, rescap)) )
                  break;

               /* update the residual capacity due to the local bound of this variable */
               if( reasonisrhs == (vals[i] > 0.0) )
               {
                  /* rhs is reason and coeff is positive, or lhs is reason and coeff is negative -> lower bound */
                  SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
                  rescap -= vals[i] * (SCIPgetVarLbAtIndex(scip, vars[i], bdchgidx, FALSE) - SCIPvarGetLbGlobal(vars[i]));
               }
               else
               {
                  /* lhs is reason and coeff is positive, or rhs is reason and coeff is negative -> upper bound */
                  SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
                  rescap -= vals[i] * (SCIPgetVarUbAtIndex(scip, vars[i], bdchgidx, FALSE) - SCIPvarGetUbGlobal(vars[i]));
               }
            }
            return SCIP_OKAY;
         }
      }
   }

   /* for a bound change on a continuous variable, all locally changed bounds are responsible */
   for( i = 0; i < nvars; ++i )
   {
      assert(vars != NULL); /* for flexelint */
      assert(vals != NULL); /* for flexelint */

      /* zero coefficients and the infered variable can be ignored */
      if( vars[i] == infervar || SCIPisZero(scip, vals[i]) )
         continue;

      if( reasonisrhs == (vals[i] > 0.0) )
      {
         /* rhs is reason and coeff is positive, or lhs is reason and coeff is negative -> lower bound is responsible */
         SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
      }
      else
      {
         /* lhs is reason and coeff is positive, or rhs is reason and coeff is negative -> upper bound is responsible */
         SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
      }
   }

   return SCIP_OKAY;
}

/** for each variable in the linear ranged row constraint, except the inferred variable, adds the bounds of all fixed
 *  variables to the conflict analysis' candidate store; the conflict analysis can be initialized
 *  with the linear constraint being the conflict detecting constraint by using NULL as inferred variable
 */
static
SCIP_RETCODE addConflictFixedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   int                   inferpos            /**< position of the inferred variable in the vars array, or -1 */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   nvars = consdata->nvars;
   assert(vars != NULL || nvars == 0);
   assert(-1 <= inferpos && inferpos < nvars);
   assert((infervar == NULL) == (inferpos == -1));
   assert(inferpos == -1 || vars != NULL);
   assert(inferpos == -1 || vars[inferpos] == infervar); /*lint !e613*/

   /* collect all fixed variables */
   for( v = nvars - 1; v >= 0; --v )
   {
      assert(vars != NULL); /* for flexelint */

      /* need to add old bounds before propagation of inferrence variable */
      if( vars[v] == infervar )
      {
	 assert(vars[v] != NULL);

	 if( !SCIPisEQ(scip, SCIPgetVarLbAtIndex(scip, vars[v], bdchgidx, FALSE), SCIPvarGetLbGlobal(vars[v])) )
	 {
	    /* @todo get boundchange index before this last boundchange and correct the index */
	    SCIP_CALL( SCIPaddConflictLb(scip, vars[v], bdchgidx) );
	 }

	 if( !SCIPisEQ(scip, SCIPgetVarUbAtIndex(scip, vars[v], bdchgidx, FALSE), SCIPvarGetUbGlobal(vars[v])) )
	 {
	    /* @todo get boundchange index before this last boundchange and correct the index */
	    SCIP_CALL( SCIPaddConflictUb(scip, vars[v], bdchgidx) );
	 }

	 continue;
      }

      /* check for fixed variables */
      if( SCIPisEQ(scip, SCIPgetVarLbAtIndex(scip, vars[v], bdchgidx, FALSE), SCIPgetVarUbAtIndex(scip, vars[v], bdchgidx, FALSE)) )
      {
	 /* add all bounds of fixed variables which lead to the boundchange of the given inference variable */
         SCIP_CALL( SCIPaddConflictLb(scip, vars[v], bdchgidx) );
         SCIP_CALL( SCIPaddConflictUb(scip, vars[v], bdchgidx) );
      }
   }

   return SCIP_OKAY;
}

/** add reasoning variables to conflict candidate queue which led to the conflict */
static
SCIP_RETCODE addConflictReasonVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables reasoning the infeasibility */
   int                   nvars,              /**< number of variables reasoning the infeasibility */
   SCIP_VAR*             var,                /**< variable which was tried to fix/tighten, or NULL */
   SCIP_Real             bound               /**< bound of variable which was tried to apply, or SCIP_INVALID */
   )
{
   int v;

   assert(scip != NULL);

   /* collect all variables for which the local bounds differ from their global bounds */
   for( v = nvars - 1; v >= 0; --v )
   {
      assert(vars != NULL);

      /* check for local bound changes variables */
      if( !SCIPisEQ(scip, SCIPvarGetLbLocal(vars[v]), SCIPvarGetLbGlobal(vars[v])) )
      {
	 /* add conflict bound */
         SCIP_CALL( SCIPaddConflictLb(scip, vars[v], 0) );
      }

      if( !SCIPisEQ(scip, SCIPvarGetUbLocal(vars[v]), SCIPvarGetUbGlobal(vars[v])) )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, vars[v], 0) );
      }
   }

   if( var != NULL )
   {
      if( bound < SCIPvarGetLbLocal(var) )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, var, 0) );
      }

      if( bound > SCIPvarGetUbLocal(var) )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, var, 0) );
      }
   }

   return SCIP_OKAY;
}

/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) activity residuals of all other variables tighten bounds of single variable
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   INFERINFO             inferinfo,          /**< inference information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
#ifndef NDEBUG
   SCIP_Real* vals;
#endif
   int nvars;
   int inferpos;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   nvars = consdata->nvars;
#ifndef NDEBUG
   vals = consdata->vals;
   assert(vars != NULL);
   assert(vals != NULL);
#endif

   /* get the position of the inferred variable in the vars array */
   inferpos = inferInfoGetPos(inferinfo);
   if( inferpos >= nvars || vars[inferpos] != infervar )
   {
      /* find inference variable in constraint */
      /**@todo use a binary search here; the variables can be sorted by variable index */
      for( inferpos = 0; inferpos < nvars && vars[inferpos] != infervar; ++inferpos )
      {}
   }
   assert(inferpos < nvars);
   assert(vars[inferpos] == infervar);
   assert(!SCIPisZero(scip, vals[inferpos]));

   switch( inferInfoGetProprule(inferinfo) )
   {
   case PROPRULE_1_RHS:
      /* the bound of the variable was tightened, because the minimal or maximal residual activity of the linear
       * constraint (only taking the other variables into account) didn't leave enough space for a larger
       * domain in order to not exceed the right hand side of the inequality
       */
      assert((vals[inferpos] > 0.0) == (boundtype == SCIP_BOUNDTYPE_UPPER));
      SCIP_CALL( addConflictBounds(scip, cons, infervar, bdchgidx, inferpos, TRUE) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_1_LHS:
      /* the bound of the variable was tightened, because the minimal or maximal residual activity of the linear
       * constraint (only taking the other variables into account) didn't leave enough space for a larger
       * domain in order to not fall below the left hand side of the inequality
       */
      assert((vals[inferpos] > 0.0) == (boundtype == SCIP_BOUNDTYPE_LOWER));
      SCIP_CALL( addConflictBounds(scip, cons, infervar, bdchgidx, inferpos, FALSE) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_1_RANGEDROW:
      /* the bound of the variable was tightened, because some variables were already fixed and the leftover only allow
       * the given inference variable to their bounds in this given ranged row
       */

      /* check that we really have a ranged row here */
      assert(!SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs));
      SCIP_CALL( addConflictFixedVars(scip, cons, infervar, bdchgidx, inferpos) );
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in linear constraint <%s> at position %d for %s bound of variable <%s>\n",
         inferInfoGetProprule(inferinfo), SCIPconsGetName(cons), inferInfoGetPos(inferinfo),
         boundtype == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper", SCIPvarGetName(infervar));
      SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      SCIPinfoMessage(scip, NULL, ";\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bounds on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< conflict detecting constraint */
   SCIP_Bool             reasonisrhs         /**< is the right hand side responsible for the conflict? */
   )
{
   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   /* add the conflicting bound for each variable of infeasible constraint to conflict candidate queue */
   SCIP_CALL( addConflictBounds(scip, cons, NULL, NULL, -1, reasonisrhs) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** check if there is any hope of tightening some bounds */
static
SCIP_Bool canTightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   int infcountmin;
   int infcountmax;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   infcountmin = consdata->minactivityneginf
      + consdata->minactivityposinf
      + consdata->minactivityneghuge
      + consdata->minactivityposhuge;
   infcountmax = consdata->maxactivityneginf
      + consdata->maxactivityposinf
      + consdata->maxactivityneghuge
      + consdata->maxactivityposhuge;

   if( infcountmin > 1 && infcountmax > 1 )
      return FALSE;

   return TRUE;
}

/** tighten upper bound */
static
SCIP_RETCODE tightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< variable position */
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldub,              /**< old upper bound */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count the total number of tightened bounds */
   SCIP_Bool             force               /**< should a possible bound change be forced even if below bound strengthening tolerance */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, newub));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   var = consdata->vars[pos];
   assert(var != NULL);

   lb = SCIPvarGetLbLocal(var);
   newub = SCIPadjustedVarUb(scip, var, newub);

   if( force || SCIPisUbBetter(scip, newub, lb, oldub) )
   {
      SCIP_VARTYPE vartype;

      SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, activity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newub=%.15g\n",
         SCIPconsGetName(cons), SCIPvarGetName(var), lb, oldub, consdata->vals[pos], consdata->minactivity, consdata->maxactivity, consdata->lhs, consdata->rhs, newub);

      vartype = SCIPvarGetType(var);

      /* tighten upper bound */
      SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(proprule, pos), force, &infeasible, &tightened) );

      if( infeasible )
      {
         SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

         /* analyze conflict */
         SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

         *cutoff = TRUE;
      }
      else if( tightened )
      {
         assert(SCIPisFeasLE(scip, SCIPvarGetUbLocal(var), oldub));
         SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), lb, SCIPvarGetUbLocal(var));

         (*nchgbds)++;

         /* if variable type was changed we might be able to upgrade the constraint */
         if( vartype != SCIPvarGetType(var) )
            consdata->upgradetried = FALSE;
      }
   }
   return SCIP_OKAY;
}

/** tighten lower bound */
static
SCIP_RETCODE tightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< variable position */
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count the total number of tightened bounds */
   SCIP_Bool             force               /**< should a possible bound change be forced even if below bound strengthening tolerance */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real ub;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, newlb));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   var = consdata->vars[pos];
   assert(var != NULL);

   ub = SCIPvarGetUbLocal(var);
   newlb = SCIPadjustedVarLb(scip, var, newlb);

   if( force || SCIPisLbBetter(scip, newlb, oldlb, ub) )
   {
      SCIP_VARTYPE vartype;

      SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, activity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newlb=%.15g\n",
         SCIPconsGetName(cons), SCIPvarGetName(var), oldlb, ub, consdata->vals[pos], consdata->minactivity, consdata->maxactivity, consdata->lhs, consdata->rhs, newlb);

      vartype = SCIPvarGetType(var);

      /* tighten lower bound */
      SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(proprule, pos), force, &infeasible, &tightened) );

      if( infeasible )
      {
         SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

         /* analyze conflict */
         SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

         *cutoff = TRUE;
      }
      else if( tightened )
      {
         assert(SCIPisFeasGE(scip, SCIPvarGetLbLocal(var), oldlb));
         SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetLbLocal(var), ub);

         (*nchgbds)++;

         /* if variable type was changed we might be able to upgrade the constraint */
         if( vartype != SCIPvarGetType(var) )
            consdata->upgradetried = FALSE;
      }
   }
   return SCIP_OKAY;
}

/** tightens bounds of a single variable due to activity bounds (easy case) */
static
SCIP_RETCODE tightenVarBoundsEasy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< position of the variable in the vars array */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count the total number of tightened bounds */
   SCIP_Bool             force               /**< should a possible bound change be forced even if below bound strengthening tolerance */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   *cutoff = FALSE;

   var = consdata->vars[pos];
   assert(var != NULL);

   /* we cannot tighten bounds of multi-aggregated variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   val = consdata->vals[pos];
   lhs = consdata->lhs;
   rhs = consdata->rhs;
   assert(!SCIPisZero(scip, val));
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPisLE(scip, lb, ub));

   /* recompute activities if needed */
   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->validactivities);
   if( !consdata->validminact )
      consdataRecomputeMinactivity(scip, consdata);
   assert(consdata->validminact);
   if( !consdata->validmaxact )
      consdataRecomputeMaxactivity(scip, consdata);
   assert(consdata->validmaxact);

   if( val > 0.0 )
   {
      /* check, if we can tighten the variable's upper bound */
      if( !SCIPisInfinity(scip, rhs) )
      {
         SCIP_Real slack;
         SCIP_Real alpha;

         /* if the minactivity is larger than the right hand side by feasibility epsilon, the constraint is infeasible */
         if( SCIPisFeasLT(scip, rhs, consdata->minactivity) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, minactivity=%.15g > rhs=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), consdata->minactivity, rhs);

            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         slack = rhs - consdata->minactivity;

         /* if the slack is zero in tolerances (or negative, but not enough to make the constraint infeasible), we set
          * it to zero
          */
         if( !SCIPisPositive(scip, slack) )
            slack = 0.0;

         alpha = val * (ub - lb);
         assert(!SCIPisNegative(scip, alpha));

         if( SCIPisSumGT(scip, alpha, slack)  || (force && SCIPisGT(scip, alpha, slack)) )
         {
            SCIP_Real newub;

            /* compute new upper bound */
            newub = lb + (slack / val);

            SCIP_CALL( tightenVarUb(scip, cons, pos, PROPRULE_1_RHS, newub, ub, cutoff, nchgbds, force) );

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               return SCIP_OKAY;
            }

            /* collect the new upper bound which is needed for the lower bound computation */
            ub = SCIPvarGetUbLocal(var);
         }
      }

      /* check, if we can tighten the variable's lower bound */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_Real slack;
         SCIP_Real alpha;

         /* if the maxactivity is smaller than the left hand side by feasibility epsilon, the constraint is infeasible */
         if( SCIPisFeasLT(scip, consdata->maxactivity, lhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, maxactivity=%.15g < lhs=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), consdata->maxactivity, lhs);

            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         slack = consdata->maxactivity - lhs;

         /* if the slack is zero in tolerances (or negative, but not enough to make the constraint infeasible), we set
          * it to zero
          */
         if( !SCIPisPositive(scip, slack) )
            slack = 0.0;

         alpha = val * (ub - lb);
         assert(!SCIPisNegative(scip, alpha));

         if( SCIPisSumGT(scip, alpha, slack) || (force && SCIPisGT(scip, alpha, slack)) )
         {
            SCIP_Real newlb;

            /* compute new lower bound */
            newlb = ub - (slack / val);

            SCIP_CALL( tightenVarLb(scip, cons, pos, PROPRULE_1_LHS, newlb, lb, cutoff, nchgbds, force) );

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               return SCIP_OKAY;
            }
         }
      }
   }
   else
   {
      /* check, if we can tighten the variable's lower bound */
      if( !SCIPisInfinity(scip, rhs) )
      {
         SCIP_Real slack;
         SCIP_Real alpha;

         /* if the minactivity is larger than the right hand side by feasibility epsilon, the constraint is infeasible */
         if( SCIPisFeasLT(scip, rhs, consdata->minactivity) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, minactivity=%.15g > rhs=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), consdata->minactivity, rhs);

            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         slack = rhs - consdata->minactivity;

         /* if the slack is zero in tolerances (or negative, but not enough to make the constraint infeasible), we set
          * it to zero
          */
         if( !SCIPisPositive(scip, slack) )
            slack = 0.0;

         alpha = val * (lb - ub);
         assert(!SCIPisNegative(scip, alpha));

         if( SCIPisSumGT(scip, alpha, slack) || (force && SCIPisGT(scip, alpha, slack)) )
         {
            SCIP_Real newlb;

            /* compute new lower bound */
            newlb = ub + slack / val;

            SCIP_CALL( tightenVarLb(scip, cons, pos, PROPRULE_1_RHS, newlb, lb, cutoff, nchgbds, force) );

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               return SCIP_OKAY;
            }
            /* collect the new lower bound which is needed for the upper bound computation */
            lb = SCIPvarGetLbLocal(var);
         }
      }

      /* check, if we can tighten the variable's upper bound */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_Real slack;
         SCIP_Real alpha;

         /* if the maxactivity is smaller than the left hand side by feasibility epsilon, the constraint is infeasible */
         if( SCIPisFeasLT(scip, consdata->maxactivity, lhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, maxactivity=%.15g < lhs=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), consdata->maxactivity, lhs);

            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         slack = consdata->maxactivity - lhs;

         /* if the slack is zero in tolerances (or negative, but not enough to make the constraint infeasible), we set
          * it to zero
          */
         if( !SCIPisPositive(scip, slack) )
            slack = 0.0;

         alpha = val * (lb - ub);
         assert(!SCIPisNegative(scip, alpha));

         if( SCIPisSumGT(scip, alpha, slack) || (force && SCIPisGT(scip, alpha, slack)) )
         {
            SCIP_Real newub;

            /* compute new upper bound */
            newub = lb - (slack / val);

            SCIP_CALL( tightenVarUb(scip, cons, pos, PROPRULE_1_LHS, newub, ub, cutoff, nchgbds, force) );

            if( *cutoff )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               return SCIP_OKAY;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bounds on given ranged row constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflictRangedRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< conflict detecting constraint */
   SCIP_VAR**            vars,               /**< variables reasoning the infeasibility */
   int                   nvars,              /**< number of variables reasoning the infeasibility */
   SCIP_VAR*             var,                /**< variable which was tried to fix/tighten, or NULL */
   SCIP_Real             bound               /**< bound of variable which was tried to apply, or SCIP_INVALID */
   )
{
#ifndef NDEBUG
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs));
#endif

   /* conflict analysis can only be applied in solving stage and if it is turned on */
   if( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && !SCIPinProbing(scip)) || !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   /* add the conflicting fixed variables of this ranged row constraint to conflict candidate queue */
   SCIP_CALL( addConflictFixedVars(scip, cons, NULL, NULL, -1) );

   /* add reasoning variables to conflict candidate queue which led to the conflict */
   SCIP_CALL( addConflictReasonVars(scip, vars, nvars, var, bound) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagate ranged rows
 *
 *  Check ranged rows for possible solutions, possibly detect infeasibility, fix variables due to having only one possible
 *  solution, tighten bounds if having only two possible solutions or add constraints which propagate a subset of
 *  variables better.
 *
 *  Example:
 *  c1: 12 x1 + 9  x2 - x3 = 0  with x1, x2 free and 1 <= x3 <= 2
 *
 *  x3 needs to be a multiple of 3, so the instance is infeasible.
 *
 *  Example:
 *  c1: 12 x1 + 9  x2 - x3 = 1  with x1, x2 free and 1 <= x3 <= 2
 *
 *  The only posible value for x3 is 2, so the variable will be fixed.
 *
 *  @todo add holes if possible
 */
static
SCIP_RETCODE rangedRowPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  nchgbds,            /**< pointer to count the number of bound changes */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** infcheckvars;
   SCIP_Real* infcheckvals;
   SCIP_Real minactinfvars;
   SCIP_Real maxactinfvars;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real feastol;
   SCIP_Real fixedact;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real absminbincoef;
   SCIP_Longint gcd;
   SCIP_Longint gcdtmp;
   SCIP_Bool minactinfvarsinvalid;
   SCIP_Bool maxactinfvarsinvalid;
   SCIP_Bool possiblegcd;
   SCIP_Bool gcdisone;
   SCIP_Bool addartconss;
   int ninfcheckvars;
   int nunfixedvars;
   int nfixedconsvars;
   int ncontvars;
   int pos;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(nchgbds != NULL);
   assert(naddconss != NULL);

   /* modifiable constraint can be changed so we do not have all necessary information */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* we already did full ranged row propagation */
   if( consdata->rangedrowpropagated == 2 )
      return SCIP_OKAY;

   /* at least three variables are needed */
   if( consdata->nvars < 3 )
      return SCIP_OKAY;

   /* do nothing on normal inequalities */
   if( SCIPisInfinity(scip, -consdata->lhs) || SCIPisInfinity(scip, consdata->rhs) )
      return SCIP_OKAY;

   /* get constraint handler data */
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   addartconss = conshdlrdata->rangedrowartcons && SCIPgetDepth(scip) < 1 && !SCIPinProbing(scip) && !SCIPinRepropagation(scip);

   /* we may add artificial constraints */
   if( addartconss )
      consdata->rangedrowpropagated = 2;
   /* we are not allowed to add artificial constraints during propagation; if nothing changed on this constraint since
    * the last rangedrowpropagation, we can stop; otherwise, we mark this constraint to be rangedrowpropagated without
    * artificial constraints
    */
   else
   {
      if( consdata->rangedrowpropagated > 0 )
         return SCIP_OKAY;

      consdata->rangedrowpropagated = 1;
   }
   fixedact = 0;
   nfixedconsvars = 0;
   /* calculate fixed activity and number of fixed variables */
   for( v = consdata->nvars - 1; v >= 0; --v )
   {
      /* all zero coefficients should be eliminated */
      assert(!SCIPisZero(scip, consdata->vals[v]));

      if( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])) )
      {
         fixedact += SCIPvarGetLbLocal(consdata->vars[v]) * consdata->vals[v];
	 ++nfixedconsvars;
      }
   }

   /* do not work with huge fixed activities */
   if( SCIPisHugeValue(scip, REALABS(fixedact)) )
      return SCIP_OKAY;

   /* compute lhs and rhs for unfixed variables only and get number of unfixed variables */
   assert(!SCIPisInfinity(scip, -fixedact) && !SCIPisInfinity(scip, fixedact));
   lhs = consdata->lhs - fixedact;
   rhs = consdata->rhs - fixedact;
   nunfixedvars = consdata->nvars - nfixedconsvars;

   /* allocate temporary memory for variables and coefficients which may lead to infeasibility */
   SCIP_CALL( SCIPallocBufferArray(scip, &infcheckvars, nunfixedvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &infcheckvals, nunfixedvars) );

   absminbincoef = SCIP_REAL_MAX;
   ncontvars = 0;
   gcdisone = TRUE;
   possiblegcd = TRUE;

   /* we now partition all unfixed variables in two groups:
    *
    * the first one contains all integral variable with integral
    * coefficient so that all variables in this group will have a gcd greater than 1, this group will be implicitly
    * given
    *
    * the second group will contain all left unfixed variables and will be saved as infcheckvars with corresponding
    * coefficients as infcheckvals, the order of these variables should be the same as in the consdata object
    */

   /* find first integral variables with integral coefficient greater than 1, thereby collecting all other unfixed
    * variables
    */
   ninfcheckvars = 0;
   v = -1;
   pos = -1;
   do
   {
      ++v;

      /* partition the variables, do not change the order of collection, because it might be used later on */
      while( v < consdata->nvars && (SCIPvarGetType(consdata->vars[v]) == SCIP_VARTYPE_CONTINUOUS ||
            !SCIPisIntegral(scip, consdata->vals[v]) || SCIPisEQ(scip, REALABS(consdata->vals[v]), 1.0)) )
      {
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])) )
         {
            if( SCIPvarGetType(consdata->vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            {
               ++ncontvars;
            }
            else if( SCIPvarIsBinary(consdata->vars[v]) )
            {
               SCIP_Real absval;

               absval = REALABS(consdata->vals[v]);

               if( absminbincoef > absval )
                  absminbincoef = absval;
            }

            gcdisone = gcdisone && SCIPisEQ(scip, REALABS(consdata->vals[v]), 1.0);
            possiblegcd = FALSE;
            infcheckvars[ninfcheckvars] = consdata->vars[v];
            infcheckvals[ninfcheckvars] = consdata->vals[v];
            ++ninfcheckvars;

	    if( pos == -1 )
	       pos = v;
         }
         ++v;
      }
   }
   while( v < consdata->nvars && SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])) );

   /* if the first group of variables is empty, we stop */
   /* @todo try to propagate/split up a constraint of the form:
    *       x_1 + ... + x_m + a_1*y_1 + ... + a_n*y_n = k + c,
    *       with k \in Z, c \in (d,d + 1], d \in Z, (a_1*y_1 + ... + a_n*y_n) \in (c-1 + d,d + 1]
    */
   if( v == consdata->nvars )
      goto TERMINATE;

   /* we need at least two non-continuous variables */
   if( ncontvars + 2 > nunfixedvars )
      goto TERMINATE;

   assert(!SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])));
   assert(SCIPisIntegral(scip, consdata->vals[v]) && SCIPvarGetType(consdata->vars[v]) != SCIP_VARTYPE_CONTINUOUS && REALABS(consdata->vals[v]) > 1.5);

   feastol = SCIPfeastol(scip);

   gcd = (SCIP_Longint)(REALABS(consdata->vals[v]) + feastol);
   assert(gcd >= 2);

   /* go on to partition the variables, do not change the order of collection, because it might be used later on;
    * calculate gcd over the first part of variables */
   for( ; v < consdata->nvars; ++v )
   {
      if( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])) )
         continue;

      if( SCIPvarIsBinary(consdata->vars[v]) )
      {
         SCIP_Real absval;

         absval = REALABS(consdata->vals[v]);

         if( absminbincoef > absval )
            absminbincoef = absval;
      }

      if( !SCIPisIntegral(scip, consdata->vals[v]) || SCIPvarGetType(consdata->vars[v]) == SCIP_VARTYPE_CONTINUOUS ||
         SCIPisEQ(scip, REALABS(consdata->vals[v]), 1.0) )
      {
         if( SCIPvarGetType(consdata->vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            ++ncontvars;

         gcdisone = gcdisone && SCIPisEQ(scip, REALABS(consdata->vals[v]), 1.0);
         possiblegcd = FALSE;
         infcheckvars[ninfcheckvars] = consdata->vars[v];
         infcheckvals[ninfcheckvars] = consdata->vals[v];

         ++ninfcheckvars;

	 if( pos == -1 )
	    pos = v;
      }
      else
      {
         assert(REALABS(consdata->vals[v]) > 1.5);

         gcdtmp = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(consdata->vals[v]) + feastol));
         assert(gcdtmp >= 1);

         if( gcdtmp == 1 )
         {
            infcheckvars[ninfcheckvars] = consdata->vars[v];
            infcheckvals[ninfcheckvars] = consdata->vals[v];

            ++ninfcheckvars;

	    if( pos == -1 )
	       pos = v;
         }
         else
            gcd = gcdtmp;
      }
   }
   assert(gcd >= 2);

   /* it should not happen that all variables are of integral type and have a gcd >= 2, this should be done by
    * normalizeCons() */
   if( ninfcheckvars == 0 )
      goto TERMINATE;

   assert(pos >= 0);

   minactinfvarsinvalid = FALSE;
   maxactinfvarsinvalid = FALSE;
   maxactinfvars = 0.0;
   minactinfvars = 0.0;

   /* calculate activities over all infcheckvars */
   for( v = ninfcheckvars - 1; v >= 0; --v )
   {
      lb = SCIPvarGetLbLocal(infcheckvars[v]);
      ub = SCIPvarGetUbLocal(infcheckvars[v]);

      if( SCIPisInfinity(scip, -lb) )
      {
         if( infcheckvals[v] < 0.0 )
            maxactinfvarsinvalid = TRUE;
         else
            minactinfvarsinvalid = TRUE;
      }
      else
      {
         if( infcheckvals[v] < 0.0 )
            maxactinfvars += infcheckvals[v] * lb;
         else
            minactinfvars += infcheckvals[v] * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         if( infcheckvals[v] > 0.0 )
            maxactinfvarsinvalid = TRUE;
         else
            minactinfvarsinvalid = TRUE;
      }
      else
      {
         if( infcheckvals[v] > 0.0 )
            maxactinfvars += infcheckvals[v] * ub;
         else
            minactinfvars += infcheckvals[v] * ub;
      }

      /* better abort on to big values */
      if( SCIPisHugeValue(scip, -minactinfvars) )
	 minactinfvarsinvalid = TRUE;
      if( SCIPisHugeValue(scip, maxactinfvars) )
	 maxactinfvarsinvalid = TRUE;

      if( minactinfvarsinvalid || maxactinfvarsinvalid )
         goto TERMINATE;
   }
   assert(!minactinfvarsinvalid && !maxactinfvarsinvalid);

   SCIPdebugMsg(scip, "minactinfvarsinvalid = %u, minactinfvars = %g, maxactinfvarsinvalid = %u, maxactinfvars = %g, gcd = %lld, ninfcheckvars = %d, ncontvars = %d\n",
      minactinfvarsinvalid, minactinfvars, maxactinfvarsinvalid, maxactinfvars, gcd, ninfcheckvars, ncontvars);

   /* @todo maybe we took the wrong variables as infcheckvars we could try to exchange integer variables */
   /* @todo if minactinfvarsinvalid or maxactinfvarsinvalid are true, try to exchange both partitions to maybe get valid
    *       activities */
   /* @todo calculate minactivity and maxactivity for all non-intcheckvars, and use this for better bounding,
    *       !!!note!!!
    *       that therefore the conflict variables in addConflictFixedVars() need to be extended by all variables which
    *       are not at their global bound
    */

   /* check if between left hand side and right hand side exist a feasible point, if not the constraint leads to
    * infeasibility */
   if( !SCIPisIntegral(scip, (lhs - maxactinfvars) / gcd) &&
      SCIPisGT(scip, SCIPceil(scip, (lhs - maxactinfvars) / gcd) * gcd, rhs - minactinfvars) )
   {
      SCIPdebugMsg(scip, "no feasible value exist, constraint <%s> lead to infeasibility", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);

      /* start conflict analysis */
      /* @todo improve conflict analysis by adding relaxed bounds */
      SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars, NULL, SCIP_INVALID) );

      *cutoff = TRUE;
   }
   else if( ncontvars == 0 )
   {
      SCIP_Longint gcdinfvars = -1;

      /* check for gcd over all infcheckvars */
      if( possiblegcd )
      {
         v = ninfcheckvars - 1;
         gcdinfvars = (SCIP_Longint)(REALABS(infcheckvals[v]) + feastol);
         assert(gcdinfvars >= 2);

         for( ; v >= 0 && gcdinfvars >= 2; --v )
         {
            gcdinfvars = SCIPcalcGreComDiv(gcdinfvars, (SCIP_Longint)(REALABS(infcheckvals[v]) + feastol));
         }
      }
      else if( gcdisone )
         gcdinfvars = 1;

      SCIPdebugMsg(scip, "gcdinfvars =%lld, possiblegcd = %u\n", gcdinfvars, possiblegcd);

      /* compute solutions for this ranged row, if all variables are of integral type with integral coefficients */
      if( gcdinfvars >= 1 )
      {
         SCIP_Real value;
         SCIP_Real value2;
         SCIP_Real minvalue = SCIP_INVALID;
         SCIP_Real maxvalue = SCIP_INVALID;
         int nsols = 0;

         value = SCIPceil(scip, minactinfvars - SCIPfeastol(scip));

         /* check how many possible solutions exist */
         while( SCIPisLE(scip, value, maxactinfvars) )
         {
	    value2 = value + gcd * (SCIPceil(scip, (lhs - value) / gcd));

	    if( SCIPisGE(scip, value2, lhs) && SCIPisLE(scip, value2, rhs) )
	    {
	       ++nsols;

	       /* early termination if we found more than two solutions */
	       if( nsols == 3 )
		  break;

	       if( minvalue == SCIP_INVALID ) /*lint !e777*/
		  minvalue = value;

	       maxvalue = value;
	    }
            value += gcdinfvars;
         }
         assert(nsols < 2 || minvalue <= maxvalue);

         /* determine last possible solution for better bounding */
         if( nsols == 3 )
         {
#ifndef NDEBUG
            SCIP_Real secondsolval = maxvalue;
#endif
            value = SCIPfloor(scip, maxactinfvars + SCIPfeastol(scip));

            /* check how many possible solutions exist */
            while( SCIPisGE(scip, value, minactinfvars) )
            {
               value2 = value + gcd * (SCIPfloor(scip, (rhs - value) / gcd));

               if( SCIPisGE(scip, value2, lhs) && SCIPisLE(scip, value2, rhs) )
               {
                  maxvalue = value;
                  assert(maxvalue > minvalue);
                  break;
               }
               value -= gcdinfvars;
            }
            assert(maxvalue > secondsolval);
         }

         SCIPdebugMsg(scip, "here nsols %s %d, minsolvalue = %g, maxsolvalue = %g, ninfcheckvars = %d, nunfixedvars = %d\n",
            nsols > 2 ? ">=" : "=", nsols, minvalue, maxvalue, ninfcheckvars, nunfixedvars);

         /* no possible solution found */
         if( nsols == 0 )
         {
	    SCIPdebugMsg(scip, "gcdinfvars = %lld, gcd = %lld, correctedlhs = %g, correctedrhs = %g\n",
               gcdinfvars, gcd, lhs, rhs);
            SCIPdebugMsg(scip, "no solution found; constraint <%s> lead to infeasibility\n", SCIPconsGetName(cons));
            SCIPdebugPrintCons(scip, cons, NULL);

	    /* start conflict analysis */
            /* @todo improve conflict analysis by adding relaxed bounds */
	    SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars, NULL, SCIP_INVALID) );

            *cutoff = TRUE;
         }
         /* if only one solution exist we can extract a new constraint or fix variables */
         else if( nsols == 1 )
         {
	    assert(minvalue == maxvalue); /*lint !e777*/

	    /* we can fix the only variable in our second set of variables */
            if( ninfcheckvars == 1 )
            {
               SCIP_Bool fixed;

	       assert(SCIPisEQ(scip, (SCIP_Real)gcdinfvars, REALABS(infcheckvals[0])));

               SCIPdebugMsg(scip, "fixing single variable <%s> with bounds [%.15g,%.15g] to %.15g\n",
                  SCIPvarGetName(infcheckvars[0]), SCIPvarGetLbLocal(infcheckvars[0]),
                  SCIPvarGetUbLocal(infcheckvars[0]), maxvalue/infcheckvals[0]);

               /* fix variable to only possible value */
               SCIP_CALL( SCIPinferVarFixCons(scip, infcheckvars[0], maxvalue/infcheckvals[0], cons,
                     getInferInt(PROPRULE_1_RANGEDROW, pos), TRUE, cutoff, &fixed) );

	       if( *cutoff )
	       {
		  /* start conflict analysis */
                  /* @todo improve conflict analysis by adding relaxed bounds */
		  SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars, NULL, SCIP_INVALID) );
	       }

               if( fixed )
                  ++(*nfixedvars);
            }
            else
            {
               /* check for exactly one unfixed variable which is not part of the infcheckvars */
               if( ninfcheckvars == nunfixedvars - 1 )
               {
                  SCIP_Real bound;
                  SCIP_Bool foundvar = FALSE;
                  SCIP_Bool fixed;
		  int w = 0;

                  assert(ninfcheckvars > 0);

                  /* find variable which is not an infcheckvar and fix it */
                  for( v = 0; v < consdata->nvars - 1; ++v )
		  {
                     if( !SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])) )
		     {
			if( w >= ninfcheckvars || consdata->vars[v] != infcheckvars[w] )
			{
#ifndef NDEBUG
                           int v2 = v + 1;
                           int w2 = w;

                           assert((nfixedconsvars == 0) ? (consdata->nvars - v - 1 == ninfcheckvars - w) : TRUE);

                           for( ; v2 < consdata->nvars && w2 < ninfcheckvars; ++v2 )
                           {
                              if( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v2]), SCIPvarGetUbLocal(consdata->vars[v2])) )
                                 continue;

                              assert(consdata->vars[v2] == infcheckvars[w2]);
                              ++w2;
                           }
                           assert(w2 == ninfcheckvars);
#endif
			   assert(SCIPisEQ(scip, (SCIP_Real)gcd, REALABS(consdata->vals[v])));

			   foundvar = TRUE;

			   if( consdata->vals[v] < 0 )
			   {
			      bound = SCIPfloor(scip, (lhs - maxvalue) / consdata->vals[v]);
			   }
			   else
			   {
			      bound = SCIPceil(scip, (lhs - maxvalue) / consdata->vals[v]);
			   }

			   SCIPdebugMsg(scip, "fixing variable <%s> with bounds [%.15g,%.15g] to %.15g\n",
                              SCIPvarGetName(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v]),
                              SCIPvarGetUbLocal(consdata->vars[v]), bound);

			   /* fix variable to only possible value */
			   SCIP_CALL( SCIPinferVarFixCons(scip, consdata->vars[v], bound, cons,
                                 getInferInt(PROPRULE_1_RANGEDROW, v), TRUE, cutoff, &fixed) );

			   if( *cutoff )
			   {
			      /* start conflict analysis */
                              /* @todo improve conflict analysis by adding relaxed bounds */
			      SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars,
                                    consdata->vars[v], bound) );
			   }

			   if( fixed )
			      ++(*nfixedvars);

			   break;
			}

			++w;
		     }
		  }

                  /* maybe last variable was the not infcheckvar */
                  if( !foundvar )
                  {
                     assert(v == consdata->nvars - 1);
		     assert(SCIPisEQ(scip, (SCIP_Real)gcd, REALABS(consdata->vals[v])));

		     if( consdata->vals[v] < 0 )
		     {
			bound = SCIPfloor(scip, (lhs - maxvalue) / consdata->vals[v]);
		     }
		     else
		     {
			bound = SCIPceil(scip, (lhs - maxvalue) / consdata->vals[v]);
		     }

                     SCIPdebugMsg(scip, "fixing variable <%s> with bounds [%.15g,%.15g] to %.15g\n",
                        SCIPvarGetName(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v]),
                        SCIPvarGetUbLocal(consdata->vars[v]), bound);

                     /* fix variable to only possible value */
		     SCIP_CALL( SCIPinferVarFixCons(scip, consdata->vars[v], bound, cons,
                           getInferInt(PROPRULE_1_RANGEDROW, v), TRUE, cutoff, &fixed) );

		     if( *cutoff )
		     {
			/* start conflict analysis */
                        /* @todo improve conflict analysis by adding relaxed bounds */
                        SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars,
                              consdata->vars[v], bound) );
		     }

                     if( fixed )
                        ++(*nfixedvars);
                  }
               }
               else if( addartconss && (SCIPisGT(scip, minvalue, minactinfvars) || SCIPisLT(scip, maxvalue, maxactinfvars)) )
               {
                  /* aggregation possible if we have two variables, but this will be done later on */
                  SCIP_CONS* newcons;
                  char name[SCIP_MAXSTRLEN];

                  /* create, add, and release new artificial constraint */
                  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_artcons_%d", SCIPconsGetName(cons), conshdlrdata->naddconss);
                  ++conshdlrdata->naddconss;

		  SCIPdebugMsg(scip, "adding artificial constraint %s\n", name);

                  SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, ninfcheckvars, infcheckvars, infcheckvals,
                        maxvalue, maxvalue, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE) );
                  SCIP_CALL( SCIPaddConsLocal(scip, newcons, NULL) );

                  SCIPdebugPrintCons(scip, newcons, NULL);

                  SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

                  ++(*naddconss);
               }
            }
         }
         /* at least two solutions */
         else
         {
            /* @todo if we found more then one solution, we may reduced domains due to dualpresolving? */

	    /* only one variable in the second set, so we can bound this variables */
            if( ninfcheckvars == 1 )
            {
               SCIP_Bool tightened;
               SCIP_Real newlb;
               SCIP_Real newub;

	       assert(SCIPisEQ(scip, (SCIP_Real)gcdinfvars, REALABS(infcheckvals[0])));

               if( infcheckvals[0] < 0 )
               {
                  newlb = maxvalue/infcheckvals[0];
                  newub = minvalue/infcheckvals[0];
               }
               else
               {
                  newlb = minvalue/infcheckvals[0];
                  newub = maxvalue/infcheckvals[0];
               }
               assert(newlb < newub);

	       if( newlb > SCIPvarGetLbLocal(infcheckvars[0]) )
	       {
		  /* update lower bound of variable */
		  SCIPdebugMsg(scip, "tightening lower bound of variable <%s> from %g to %g\n",
                     SCIPvarGetName(infcheckvars[0]), SCIPvarGetLbLocal(infcheckvars[0]), newlb);

		  /* tighten variable lower bound to minimal possible value */
		  SCIP_CALL( SCIPinferVarLbCons(scip, infcheckvars[0], newlb, cons,
                        getInferInt(PROPRULE_1_RANGEDROW, pos), TRUE, cutoff, &tightened) );

		  if( *cutoff )
		  {
		     /* start conflict analysis */
                     /* @todo improve conflict analysis by adding relaxed bounds */
		     SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars, NULL, SCIP_INVALID) );
		  }

		  if( tightened )
		     ++(*nchgbds);
	       }

	       if( newub < SCIPvarGetUbLocal(infcheckvars[0]) )
	       {
		  /* update upper bound of variable */
		  SCIPdebugMsg(scip, "tightening upper bound of variable <%s> from %g to %g\n",
                     SCIPvarGetName(infcheckvars[0]), SCIPvarGetUbLocal(infcheckvars[0]), newub);

		  /* tighten variable upper bound to maximal possible value */
		  SCIP_CALL( SCIPinferVarUbCons(scip, infcheckvars[0], newub, cons,
                        getInferInt(PROPRULE_1_RANGEDROW, pos), TRUE, cutoff, &tightened) );

		  if( *cutoff )
		  {
		     /* start conflict analysis */
                     /* @todo improve conflict analysis by adding relaxed bounds */
		     SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars, NULL, SCIP_INVALID) );
		  }

		  if( tightened )
		     ++(*nchgbds);
	       }
            }
            /* check if we have only one not infcheckvars, if so we can tighten this variable */
            else if( ninfcheckvars == nunfixedvars - 1 )
            {
               SCIP_Bool foundvar = FALSE;
               SCIP_Bool tightened;
               SCIP_Real newlb;
               SCIP_Real newub;
	       int w = 0;

               assert(ninfcheckvars > 0);
               assert(minvalue < maxvalue);

               /* find variable which is not an infcheckvar and fix it */
               for( v = 0; v < consdata->nvars - 1; ++v )
	       {
                  if( !SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])) )
		  {
		     if( w >= ninfcheckvars || consdata->vars[v] != infcheckvars[w] )
		     {
#ifndef NDEBUG
                        int v2 = v + 1;
                        int w2 = w;

                        assert((nfixedconsvars == 0) ? (consdata->nvars - v - 1 == ninfcheckvars - w) : TRUE);

                        for( ; v2 < consdata->nvars && w2 < ninfcheckvars; ++v2 )
                        {
                           if( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v2]), SCIPvarGetUbLocal(consdata->vars[v2])) )
                              continue;

                           assert(consdata->vars[v2] == infcheckvars[w2]);
                           ++w2;
                        }
                        assert(w2 == ninfcheckvars);
#endif

			assert(SCIPisEQ(scip, (SCIP_Real)gcd, REALABS(consdata->vals[v])));
			foundvar = TRUE;

			if( consdata->vals[v] < 0 )
			{
			   newlb = SCIPfloor(scip, (rhs - minvalue) / consdata->vals[v]);
			   newub = SCIPfloor(scip, (lhs - maxvalue) / consdata->vals[v]);
			}
			else
			{
			   newlb = SCIPceil(scip, (lhs - maxvalue) / consdata->vals[v]);
			   newub = SCIPceil(scip, (rhs - minvalue) / consdata->vals[v]);
			}
			assert(SCIPisLE(scip, newlb, newub));

			if( newlb > SCIPvarGetLbLocal(consdata->vars[v]) )
			{
			   /* update lower bound of variable */
			   SCIPdebugMsg(scip, "tightening lower bound of variable <%s> from %g to %g\n",
                              SCIPvarGetName(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v]), newlb);

			   /* tighten variable lower bound to minimal possible value */
			   SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vars[v], newlb, cons,
                                 getInferInt(PROPRULE_1_RANGEDROW, v), TRUE, cutoff, &tightened) );

			   if( *cutoff )
			   {
			      /* start conflict analysis */
                              /* @todo improve conflict analysis by adding relaxed bounds */
			      SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars,
                                    consdata->vars[v], newlb) );
			   }

			   if( tightened )
			      ++(*nchgbds);
			}

			if( newub < SCIPvarGetUbLocal(consdata->vars[v]) )
			{
			   /* update upper bound of variable */
			   SCIPdebugMsg(scip, "tightening upper bound of variable <%s> from %g to %g\n",
                              SCIPvarGetName(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v]), newub);

			   /* tighten variable upper bound to maximal possible value */
			   SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vars[v], newub, cons,
                                 getInferInt(PROPRULE_1_RANGEDROW, v), TRUE, cutoff, &tightened) );

			   if( *cutoff )
			   {
			      /* start conflict analysis */
                              /* @todo improve conflict analysis by adding relaxed bounds */
			      SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars,
                                    consdata->vars[v], newub) );
			   }

			   if( tightened )
			      ++(*nchgbds);
			}

			break;
		     }

		     ++w;
		  }
	       }

               /* maybe last variable was the not infcheckvar */
               if( !foundvar )
               {
                  assert(v == consdata->nvars - 1);
		  assert(SCIPisEQ(scip, (SCIP_Real)gcd, REALABS(consdata->vals[v])));

		  if( consdata->vals[v] < 0 )
		  {
		     newlb = SCIPfloor(scip, (rhs - minvalue) / consdata->vals[v]);
		     newub = SCIPfloor(scip, (lhs - maxvalue) / consdata->vals[v]);
		  }
		  else
		  {
		     newlb = SCIPceil(scip, (lhs - maxvalue) / consdata->vals[v]);
		     newub = SCIPceil(scip, (rhs - minvalue) / consdata->vals[v]);
		  }
                  assert(SCIPisLE(scip, newlb, newub));

		  if( newlb > SCIPvarGetLbLocal(consdata->vars[v]) )
		  {
		     /* update lower bound of variable */
		     SCIPdebugMsg(scip, "tightening lower bound of variable <%s> from %g to %g\n",
                        SCIPvarGetName(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v]), newlb);

                     /* tighten variable lower bound to minimal possible value */
		     SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vars[v], newlb, cons,
                           getInferInt(PROPRULE_1_RANGEDROW, v), TRUE, cutoff, &tightened) );

		     if( *cutoff )
		     {
			/* start conflict analysis */
                        /* @todo improve conflict analysis by adding relaxed bounds */
                        SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars, consdata->vars[v], newlb) );
		     }

		     if( tightened )
			++(*nchgbds);
		  }

		  if( newub < SCIPvarGetUbLocal(consdata->vars[v]) )
		  {
		     /* update upper bound of variable */
		     SCIPdebugMsg(scip, "tightening upper bound of variable <%s> from %g to %g\n",
                        SCIPvarGetName(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v]), newub);

		     /* tighten variable upper bound to maximal possible value */
		     SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vars[v], newub, cons,
                           getInferInt(PROPRULE_1_RANGEDROW, v), TRUE, cutoff, &tightened) );

		     if( *cutoff )
		     {
			/* start conflict analysis */
                        /* @todo improve conflict analysis by adding relaxed bounds */
                        SCIP_CALL( analyzeConflictRangedRow(scip, cons, infcheckvars, ninfcheckvars, consdata->vars[v], newub) );
		     }

		     if( tightened )
			++(*nchgbds);
		  }
               }
            }
            /* at least two solutions and more than one variable, so we add a new constraint which bounds the feasible
             * region for our infcheckvars, if possible
             */
            else if( addartconss && (SCIPisGT(scip, minvalue, minactinfvars) || SCIPisLT(scip, maxvalue, maxactinfvars)) )
            {
               SCIP_CONS* newcons;
               char name[SCIP_MAXSTRLEN];
               SCIP_Real newlhs;
               SCIP_Real newrhs;

               assert(maxvalue > minvalue);

               if( SCIPisGT(scip, minvalue, minactinfvars) )
                  newlhs = minvalue;
               else
                  newlhs = -SCIPinfinity(scip);

               if( SCIPisLT(scip, maxvalue, maxactinfvars) )
                  newrhs = maxvalue;
               else
                  newrhs = SCIPinfinity(scip);

	       if( !SCIPisInfinity(scip, -newlhs) || !SCIPisInfinity(scip, newrhs) )
	       {
		  /* create, add, and release new artificial constraint */
		  (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_artcons1_%d", SCIPconsGetName(cons), conshdlrdata->naddconss);
		  ++conshdlrdata->naddconss;

		  SCIPdebugMsg(scip, "adding artificial constraint %s\n", name);

		  SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, ninfcheckvars, infcheckvars, infcheckvals, newlhs, newrhs,
			TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE) );
		  SCIP_CALL( SCIPaddConsLocal(scip, newcons, NULL) );

		  SCIPdebugPrintCons(scip, newcons, NULL);
		  SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

		  ++(*naddconss);
	       }
               /* @todo maybe add constraint for all variables which are not infcheckvars, lhs should be minvalue, rhs
                *       should be maxvalue */
            }
         }
      }
   }
   else if( addartconss && ncontvars < ninfcheckvars )
   {
      SCIP_Real maxact = 0.0;
      SCIP_Real minact = 0.0;
      int w = 0;

      /* compute activities of non-infcheckvars */
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( w < ninfcheckvars && consdata->vars[v] == infcheckvars[w] )
         {
            ++w;
            continue;
         }

         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->vars[v]), SCIPvarGetUbLocal(consdata->vars[v])) )
         {
            if( SCIPvarIsBinary(consdata->vars[v]) )
            {
               if( consdata->vals[v] > 0.0 )
                  maxact += consdata->vals[v];
               else
                  minact += consdata->vals[v];
            }
            else
            {
               SCIP_Real tmpval;

               assert(SCIPvarIsIntegral(consdata->vars[v]));

               if( consdata->vals[v] > 0.0 )
               {
                  tmpval = consdata->vals[v] * SCIPvarGetLbLocal(consdata->vars[v]);

                  if( SCIPisHugeValue(scip, -tmpval) )
                     break;

                  minact += tmpval;

                  tmpval = consdata->vals[v] * SCIPvarGetUbLocal(consdata->vars[v]);

                  if( SCIPisHugeValue(scip, tmpval) )
                     break;

                  maxact += tmpval;
               }
               else
               {
                  tmpval = consdata->vals[v] * SCIPvarGetUbLocal(consdata->vars[v]);

                  if( SCIPisHugeValue(scip, -tmpval) )
                     break;

                  minact += tmpval;

                  tmpval = consdata->vals[v] * SCIPvarGetLbLocal(consdata->vars[v]);

                  if( SCIPisHugeValue(scip, tmpval) )
                     break;

                  maxact += tmpval;
               }
            }
         }
      }
      if( v == consdata->nvars && !SCIPisHugeValue(scip, -minact) && !SCIPisHugeValue(scip, maxact) )
      {
         SCIP_CONS* newcons;
         char name[SCIP_MAXSTRLEN];
         SCIP_Real newlhs;
         SCIP_Real newrhs;

         assert(maxact > minact);
         assert(w == ninfcheckvars);

         newlhs = lhs - maxact;
         newrhs = rhs - minact;
         assert(newlhs < newrhs);

         /* create, add, and release new artificial constraint */
         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_artcons2_%d", SCIPconsGetName(cons), conshdlrdata->naddconss);
         ++conshdlrdata->naddconss;

         SCIPdebugMsg(scip, "adding artificial constraint %s\n", name);

         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, name, ninfcheckvars, infcheckvars, infcheckvals, newlhs, newrhs,
               TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddConsLocal(scip, newcons, NULL) );

         SCIPdebugPrintCons(scip, newcons, NULL);
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

         ++(*naddconss);
      }
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &infcheckvals);
   SCIPfreeBufferArray(scip, &infcheckvars);

   return SCIP_OKAY;
}

/** tightens bounds of a single variable due to activity bounds */
static
SCIP_RETCODE tightenVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int                   pos,                /**< position of the variable in the vars array */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count the total number of tightened bounds */
   SCIP_Bool             force               /**< should a possible bound change be forced even if below bound strengthening tolerance */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Bool minisrelax;
   SCIP_Bool maxisrelax;
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   *cutoff = FALSE;

   var = consdata->vars[pos];

   /* we cannot tighten bounds of multi-aggregated variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   val = consdata->vals[pos];
   lhs = consdata->lhs;
   rhs = consdata->rhs;
   consdataGetActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
      &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);
   assert(var != NULL);
   assert(!SCIPisZero(scip, val));
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPisLE(scip, lb, ub));

   if( val > 0.0 )
   {
      /* check, if we can tighten the variable's bounds */
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) && !minisrelax )
      {
         SCIP_Real newub;

         newub = (rhs - minresactivity)/val;

         if( !SCIPisInfinity(scip, newub) &&
            ((force && SCIPisLT(scip, newub, ub)) || (SCIPvarIsIntegral(var) && SCIPisFeasLT(scip, newub, ub)) || SCIPisUbBetter(scip, newub, lb, ub)) )
         {
            SCIP_Bool activityunreliable;
            activityunreliable = SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity);

            /* check minresactivities for reliability */
            if( activityunreliable )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
               newub = (rhs - minresactivity)/val;
               activityunreliable = SCIPisInfinity(scip, -minresactivity) ||
                  (!SCIPisUbBetter(scip, newub, lb, ub) && (!SCIPisFeasLT(scip, newub, ub) || !SCIPvarIsIntegral(var))
                     && (!force || !SCIPisLT(scip, newub, ub)));
            }

            if( !activityunreliable )
            {
               /* tighten upper bound */
               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newub=%.15g\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newub);
               SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_RHS, pos), force,
                     &infeasible, &tightened) );
               if( infeasible )
               {
                  SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

                  /* analyze conflict */
                  SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  ub = SCIPvarGetUbLocal(var); /* get bound again: it may be additionally modified due to integrality */
                  assert(SCIPisFeasLE(scip, ub, newub));
                  (*nchgbds)++;

                  SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
               }
            }
         }
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) && !maxisrelax )
      {
         SCIP_Real newlb;

         newlb = (lhs - maxresactivity)/val;
         if( !SCIPisInfinity(scip, -newlb) &&
            ((force && SCIPisGT(scip, newlb, lb)) || (SCIPvarIsIntegral(var) && SCIPisFeasGT(scip, newlb, lb)) || SCIPisLbBetter(scip, newlb, lb, ub)) )
         {
            /* check maxresactivities for reliability */
            if( SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
               newlb = (lhs - maxresactivity)/val;

               if( SCIPisInfinity(scip, maxresactivity) || (!SCIPisLbBetter(scip, newlb, lb, ub) 
                     && (!SCIPisFeasGT(scip, newlb, lb) || !SCIPvarIsIntegral(var)) 
                     && (!force || !SCIPisGT(scip, newlb, lb))) )
                  return SCIP_OKAY;
            }

            /* tighten lower bound */
            SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newlb=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newlb);
            SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_LHS, pos), force,
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               lb = SCIPvarGetLbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasGE(scip, lb, newlb));
               (*nchgbds)++;
               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
   }
   else
   {
      /* check, if we can tighten the variable's bounds */
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) && !minisrelax )
      {
         SCIP_Real newlb;

         newlb = (rhs - minresactivity)/val;
         if( !SCIPisInfinity(scip, -newlb) &&
            ((force && SCIPisGT(scip, newlb, lb)) || (SCIPvarIsIntegral(var) && SCIPisFeasGT(scip, newlb, lb)) || SCIPisLbBetter(scip, newlb, lb, ub)) )
         {
            SCIP_Bool activityunreliable;
            activityunreliable = SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity);
            /* check minresactivities for reliability */
            if( activityunreliable )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
               newlb = (rhs - minresactivity)/val;

               activityunreliable = SCIPisInfinity(scip, -minresactivity) 
                  || (!SCIPisLbBetter(scip, newlb, lb, ub) && (!SCIPisFeasGT(scip, newlb, lb) || !SCIPvarIsIntegral(var))
                     && (!force || !SCIPisGT(scip, newlb, lb)));
            }

            if( !activityunreliable )
            {
               /* tighten lower bound */
               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g] -> newlb=%.15g\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newlb);
               SCIP_CALL( SCIPinferVarLbCons(scip, var, newlb, cons, getInferInt(PROPRULE_1_RHS, pos), force,
                     &infeasible, &tightened) );
               if( infeasible )
               {
                  SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), newlb, ub);

                  /* analyze conflict */
                  SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

                  *cutoff = TRUE;
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  lb = SCIPvarGetLbLocal(var); /* get bound again: it may be additionally modified due to integrality */
                  assert(SCIPisFeasGE(scip, lb, newlb));
                  (*nchgbds)++;
                  SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
               }
            }
         }
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) && !maxisrelax )
      {
         SCIP_Real newub;

         newub = (lhs - maxresactivity)/val;
         if(  !SCIPisInfinity(scip, newub) &&
            ((force && SCIPisLT(scip, newub, ub)) || (SCIPvarIsIntegral(var) && SCIPisFeasLT(scip, newub, ub)) || SCIPisUbBetter(scip, newub, lb, ub)) )
         {
            /* check maxresactivities for reliability */
            if( SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
            {
               consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
               newub = (lhs - maxresactivity)/val;

               if( SCIPisInfinity(scip, maxresactivity) || (!SCIPisUbBetter(scip, newub, lb, ub) 
                     && (!SCIPisFeasLT(scip, newub, ub) && !SCIPvarIsIntegral(var))
                     && (!force || !SCIPisLT(scip, newub, ub))) )
                  return SCIP_OKAY;
            }

            /* tighten upper bound */
            SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, old bds=[%.15g,%.15g], val=%.15g, resactivity=[%.15g,%.15g], sides=[%.15g,%.15g], newub=%.15g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub, val, minresactivity, maxresactivity, lhs, rhs, newub);
            SCIP_CALL( SCIPinferVarUbCons(scip, var, newub, cons, getInferInt(PROPRULE_1_LHS, pos), force,
                  &infeasible, &tightened) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: cutoff  <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, newub);

               /* analyze conflict */
               SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               ub = SCIPvarGetUbLocal(var); /* get bound again: it may be additionally modified due to integrality */
               assert(SCIPisFeasLE(scip, ub, newub));
               (*nchgbds)++;
               SCIPdebugMsg(scip, "linear constraint <%s>: tighten <%s>, new bds=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), lb, ub);
            }
         }
      }
   }

   return SCIP_OKAY;
}

#define MAXTIGHTENROUNDS 10

/** tightens bounds of variables in constraint due to activity bounds */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             maxeasyactivitydelta,/**< maximum activity delta to run easy propagation on linear constraint */
   SCIP_Bool             sortvars,           /**< should variables be used in sorted order? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   unsigned int tightenmode;
   int nvars;
   int nrounds;
   int lastchange;
   int oldnchgbds;
#ifndef SCIP_DEBUG
   int oldnchgbdstotal;
#endif
   int v;
   SCIP_Bool force;
   SCIP_Bool easycase;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgbds != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   /* we cannot tighten variables' bounds, if the constraint may be not complete */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* if a constraint was created after presolve, then it may hold fixed variables
    * if there are even multi-aggregated variables, then we cannot do bound tightening on these
    * thus, ensure here again that variable fixings have been applied
    */
   SCIP_CALL( applyFixings(scip, cons, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;

   /* check if constraint has any chances of tightening bounds */
   if( !canTightenBounds(scip, cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;
   force = (nvars == 1) && !SCIPconsIsModifiable(cons);

   /* we are at the root node or during presolving */
   if( SCIPgetDepth(scip) < 1 )
      tightenmode = 2;
   else
      tightenmode = 1;

   /* stop if we already tightened the constraint and the tightening is not forced */
   if( !force && (consdata->boundstightened >= tightenmode) ) /*lint !e574*/
      return SCIP_OKAY;

   /* ensure that the variables are properly sorted */
   if( sortvars && SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE && !consdata->binvarssorted )
   {
      SCIP_CALL( consdataSort(scip, consdata) );
      assert(consdata->binvarssorted);
   }

   /* update maximal activity delta if necessary */
   if( consdata->maxactdelta == SCIP_INVALID ) /*lint !e777*/
      consdataRecomputeMaxActivityDelta(scip, consdata);

   assert(consdata->maxactdelta != SCIP_INVALID); /*lint !e777*/
   assert(!SCIPisFeasNegative(scip, consdata->maxactdelta));
   checkMaxActivityDelta(scip, consdata);

   /* this may happen if all variables are fixed */
   if( SCIPisFeasZero(scip, consdata->maxactdelta) )
      return SCIP_OKAY;

   if( !SCIPisInfinity(scip, consdata->maxactdelta) )
   {
      SCIP_Real slack;
      SCIP_Real surplus;
      SCIP_Real minactivity;
      SCIP_Real maxactivity;
      SCIP_Bool minisrelax;
      SCIP_Bool maxisrelax;

      /* use maximal activity delta to skip propagation (cannot deduce anything) */
      consdataGetActivityBounds(scip, consdata, FALSE, &minactivity, &maxactivity, &minisrelax, &maxisrelax);
      assert(!SCIPisInfinity(scip, minactivity));
      assert(!SCIPisInfinity(scip, -maxactivity));

      slack = (SCIPisInfinity(scip, consdata->rhs) || SCIPisInfinity(scip, -minactivity)) ? SCIPinfinity(scip) : (consdata->rhs - minactivity);
      surplus = (SCIPisInfinity(scip, -consdata->lhs) || SCIPisInfinity(scip, maxactivity)) ? SCIPinfinity(scip) : (maxactivity - consdata->lhs);

      /* check if the constraint will propagate */
      if( SCIPisLE(scip, consdata->maxactdelta, MIN(slack, surplus)) )
         return SCIP_OKAY;
   }

   /* check if we can use fast implementation for easy and numerically well behaved cases */
   easycase = SCIPisLT(scip, consdata->maxactdelta, maxeasyactivitydelta);

   /* as long as the bounds might be tightened again, try to tighten them; abort after a maximal number of rounds */
   lastchange = -1;
   oldnchgbds = 0;

#ifndef SCIP_DEBUG
   oldnchgbdstotal = *nchgbds;
#endif

   for( nrounds = 0; (force || consdata->boundstightened < tightenmode) && nrounds < MAXTIGHTENROUNDS; ++nrounds ) /*lint !e574*/
   {
      /* mark the constraint to have the variables' bounds tightened */
      consdata->boundstightened = (unsigned int)tightenmode;

      /* try to tighten the bounds of each variable in the constraint. During solving process, the binary variable
       * sorting enables skipping variables
       */
      v = 0;
      while( v < nvars && v != lastchange && !(*cutoff) )
      {
         oldnchgbds = *nchgbds;

         assert(!sortvars || SCIPgetStage(scip) < SCIP_STAGE_SOLVING || consdata->binvarssorted);

         if( easycase )
         {
            SCIP_CALL( tightenVarBoundsEasy(scip, cons, v, cutoff, nchgbds, force) );
         }
         else
         {
            SCIP_CALL( tightenVarBounds(scip, cons, v, cutoff, nchgbds, force) );
         }

         /* if there was no progress, skip the rest of the binary variables */
         if( *nchgbds > oldnchgbds )
         {
            lastchange = v;
            ++v;
         }
         else if( consdata->binvarssorted && v < consdata->nbinvars - 1
            && !SCIPisFeasEQ(scip, SCIPvarGetUbLocal(consdata->vars[v]), SCIPvarGetLbLocal(consdata->vars[v])) )
            v = consdata->nbinvars;
         else
            ++v;
      }

#ifndef SCIP_DEBUG
      SCIPdebugMessage("linear constraint <%s> found %d bound changes in round %d\n", SCIPconsGetName(cons),
         *nchgbds - oldnchgbdstotal, nrounds);
      oldnchgbdstotal += oldnchgbds;
#endif
   }

#ifndef NDEBUG
   if( force && SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      assert(*cutoff || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(consdata->vars[0]), SCIPvarGetUbLocal(consdata->vars[0])));
#endif

   return SCIP_OKAY;
}

/** checks linear constraint for feasibility of given solution or current solution */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_SOL*             sol,                /**< solution to be checked, or NULL for current solution */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             checkrelmaxabs,     /**< Should the violation for a constraint with side 0.0 be checked relative
                                              *   to 1.0 (FALSE) or to the maximum absolute value in the activity (TRUE)? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;
   SCIP_Real absviol;
   SCIP_Real relviol;
   SCIP_Real lhsviol;
   SCIP_Real rhsviol;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   SCIPdebugMsg(scip, "checking linear constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   if( consdata->row != NULL )
   {
      if( !checklprows && SCIProwIsInLP(consdata->row) )
         return SCIP_OKAY;
      else if( sol == NULL && !SCIPhasCurrentNodeLP(scip) )
         activity = consdataComputePseudoActivity(scip, consdata);
      else
         activity = SCIPgetRowSolActivity(scip, consdata->row, sol);
   }
   else
      activity = consdataGetActivity(scip, consdata, sol);

   SCIPdebugMsg(scip, "  consdata activity=%.15g (lhs=%.15g, rhs=%.15g, row=%p, checklprows=%u, rowinlp=%u, sol=%p, hascurrentnodelp=%u)\n",
      activity, consdata->lhs, consdata->rhs, (void*)consdata->row, checklprows,
      consdata->row == NULL ? 0 : SCIProwIsInLP(consdata->row), (void*)sol,
      consdata->row == NULL ? FALSE : SCIPhasCurrentNodeLP(scip));

   /* calculate absolute and relative bound violations */
   lhsviol = consdata->lhs - activity;
   rhsviol = activity - consdata->rhs;

   absviol = 0.0;
   relviol = 0.0;
   if( (lhsviol > 0) && (lhsviol > rhsviol) )
   {
      absviol = lhsviol;
      relviol = SCIPrelDiff(consdata->lhs, activity);
   }
   else if( rhsviol > 0 )
   {
      absviol = rhsviol;
      relviol = SCIPrelDiff(activity, consdata->rhs);
   }

   /* the activity of pseudo solutions may be invalid if it comprises positive and negative infinity contributions; we
    * return infeasible for safety
    */
   if( activity == SCIP_INVALID ) /*lint !e777*/
   {
      assert(sol == NULL);
      *violated = TRUE;

      /* set violation of invalid pseudo solutions */
      absviol = SCIP_INVALID;
      relviol = SCIP_INVALID;

      /* reset constraint age since we are in enforcement */
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }
   /* check with relative tolerances (the default) */
   else if( !consdata->checkabsolute && (SCIPisFeasLT(scip, activity, consdata->lhs) || SCIPisFeasGT(scip, activity, consdata->rhs)) )
   {
      /* the "normal" check: one of the two sides is violated */
      if( !checkrelmaxabs )
      {
         *violated = TRUE;

         /* only reset constraint age if we are in enforcement */
         if( sol == NULL )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
      }
      /* the (much) more complicated check: we try to disregard random noise and violations of a 0.0 side which are
       * small compared to the absolute values occurring in the activity
       */
      else
      {
         SCIP_Real maxabs;
         SCIP_Real coef;
         SCIP_Real absval;
         SCIP_Real solval;
         int v;

         maxabs = 1.0;

         /* compute maximum absolute value */
         for( v = 0; v < consdata->nvars; ++v )
         {
            if( consdata->vals != NULL )
            {
               coef = consdata->vals[v];
            }
            else
               coef = 1.0;

            solval = SCIPgetSolVal(scip, sol, consdata->vars[v]);
            absval = REALABS( coef * solval );
            maxabs = MAX( maxabs, absval );
         }

         /* regard left hand side, first */
         if( SCIPisFeasLT(scip, activity, consdata->lhs) )
         {
            /* check whether violation is random noise */
            if( (consdata->lhs - activity) <= (1e-15 * maxabs) )
            {
               SCIPdebugMsg(scip, "  lhs violated due to random noise: violation=%16.9g, maxabs=%16.9g\n",
                  consdata->lhs - activity, maxabs);
               SCIPdebug( SCIP_CALL( consPrintConsSol(scip, cons, sol, NULL) ) );

               /* only increase constraint age if we are in enforcement */
               if( sol == NULL )
               {
                  SCIP_CALL( SCIPincConsAge(scip, cons) );
               }
            }
            /* lhs is violated and lhs is 0.0: use relative tolerance w.r.t. largest absolute value */
            else if( SCIPisZero(scip, consdata->lhs) )
            {
               if( (consdata->lhs - activity) <= (SCIPfeastol(scip) * maxabs) )
               {
                  SCIPdebugMsg(scip, "  lhs violated absolutely (violation=%16.9g), but feasible when using relative tolerance w.r.t. maximum absolute value (%16.9g)\n",
                     consdata->lhs - activity, maxabs);
                  SCIPdebug( SCIP_CALL( consPrintConsSol(scip, cons, sol, NULL) ) );

                  /* only increase constraint age if we are in enforcement */
                  if( sol == NULL )
                  {
                     SCIP_CALL( SCIPincConsAge(scip, cons) );
                  }
               }
               else
               {
                  *violated = TRUE;

                  /* only reset constraint age if we are in enforcement */
                  if( sol == NULL )
                  {
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  }
               }
            }
            else
            {
               *violated = TRUE;

               /* only reset constraint age if we are in enforcement */
               if( sol == NULL )
               {
                  SCIP_CALL( SCIPresetConsAge(scip, cons) );
               }
            }
         }

         /* now regard right hand side */
         if( SCIPisFeasGT(scip, activity, consdata->rhs) )
         {
            /* check whether violation is random noise */
            if( (activity - consdata->rhs) <= (1e-15 * maxabs) )
            {
               SCIPdebugMsg(scip, "  rhs violated due to random noise: violation=%16.9g, maxabs=%16.9g\n",
                  activity - consdata->rhs, maxabs);
               SCIPdebug( SCIP_CALL( consPrintConsSol(scip, cons, sol, NULL) ) );

               /* only increase constraint age if we are in enforcement */
               if( sol == NULL )
               {
                  SCIP_CALL( SCIPincConsAge(scip, cons) );
               }
            }
            /* rhs is violated and rhs is 0.0, use relative tolerance w.r.t. largest absolute value */
            else if( SCIPisZero(scip, consdata->rhs) )
            {
               if( (activity - consdata->rhs) <= (SCIPfeastol(scip) * maxabs) )
               {
                  SCIPdebugMsg(scip, "  rhs violated absolutely (violation=%16.9g), but feasible when using relative tolerance w.r.t. maximum absolute value (%16.9g)\n",
                     activity - consdata->rhs, maxabs);
                  SCIPdebug( SCIP_CALL( consPrintConsSol(scip, cons, sol, NULL) ) );

                  /* only increase constraint age if we are in enforcement */
                  if( sol == NULL )
                  {
                     SCIP_CALL( SCIPincConsAge(scip, cons) );
                  }
               }
               else
               {
                  *violated = TRUE;

                  /* only reset constraint age if we are in enforcement */
                  if( sol == NULL )
                  {
                     SCIP_CALL( SCIPresetConsAge(scip, cons) );
                  }
               }
            }
            else
            {
               *violated = TRUE;

               /* only reset constraint age if we are in enforcement */
               if( sol == NULL )
               {
                  SCIP_CALL( SCIPresetConsAge(scip, cons) );
               }
            }
         }
      }
   }
   /* check with absolute tolerances */
   else if( consdata->checkabsolute &&
      ((!SCIPisInfinity(scip, -consdata->lhs) && SCIPisGT(scip, consdata->lhs-activity, SCIPfeastol(scip))) ||
       (!SCIPisInfinity(scip,  consdata->rhs) && SCIPisGT(scip, activity-consdata->rhs, SCIPfeastol(scip)))) )
   {
      *violated = TRUE;

      /* only reset constraint age if we are in enforcement */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }
   else
   {
      /* only increase constraint age if we are in enforcement */
      if( sol == NULL )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }
   }

   /* update absolute and relative violation of the solution */
   if( sol != NULL )
      SCIPupdateSolLPConsViolation(scip, sol, absviol, relviol);

   return SCIP_OKAY;
}

/** creates an LP row in a linear constraint data */
static
SCIP_RETCODE createRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &consdata->row, SCIPconsGetHdlr(cons), SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );

   SCIP_CALL( SCIPaddVarsToRow(scip, consdata->row, consdata->nvars, consdata->vars, consdata->vals) );

   return SCIP_OKAY;
}

/** adds linear constraint as cut to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff was found */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      /* convert consdata object into LP row */
      SCIP_CALL( createRow(scip, cons) );
   }
   assert(consdata->row != NULL);

   if( consdata->nvars == 0 )
   {
      SCIPdebugMsg(scip, "Empty linear constraint enters LP: <%s>\n", SCIPconsGetName(cons));
   }

   /* insert LP row as cut */
   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMsg(scip, "adding relaxation of linear constraint <%s>: ", SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintRow(scip, consdata->row, NULL)) );
      /* if presolving is turned off, the row might be trivial */
      if ( ! SCIPisInfinity(scip, -consdata->lhs) || ! SCIPisInfinity(scip, consdata->rhs) )
      {
         SCIP_CALL( SCIPaddRow(scip, consdata->row, FALSE, cutoff) );
      }
#ifndef NDEBUG
      else
      {
         int pr;
         int cr;
         SCIP_CALL( SCIPgetIntParam(scip, "presolving/maxrounds", &pr) );
         SCIP_CALL( SCIPgetIntParam(scip, "constraints/linear/maxprerounds", &cr) );
         assert( pr == 0 || cr == 0 );
      }
#endif
   }

   return SCIP_OKAY;
}

/** separates linear constraint: adds linear constraint as cut, if violated by given solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool             separatecards,      /**< should knapsack cardinality cuts be generated? */
   SCIP_Bool             separateall,        /**< should all constraints be subject to cardinality cut generation instead of only
                                              *   the ones with non-zero dual value? */
   int*                  ncuts,              /**< pointer to add up the number of found cuts */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool violated;
   int oldncuts;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);

   consdata = SCIPconsGetData(cons);
   assert(ncuts != NULL);
   assert(consdata != NULL);

   oldncuts = *ncuts;
   *cutoff = FALSE;

   SCIP_CALL( checkCons(scip, cons, sol, (sol != NULL), conshdlrdata->checkrelmaxabs, &violated) );

   if( violated )
   {
      /* insert LP row as cut */
      SCIP_CALL( addRelaxation(scip, cons, cutoff) );
      (*ncuts)++;
   }
   else if( !SCIPconsIsModifiable(cons) && separatecards )
   {
      /* relax linear constraint into knapsack constraint and separate lifted cardinality cuts */
      if( !separateall && sol == NULL )
      {
         /* we only want to call the knapsack cardinality cut separator for rows that have a non-zero dual solution */
         if( consdata->row != NULL && SCIProwIsInLP(consdata->row) )
         {
            SCIP_Real dualsol;

            dualsol = SCIProwGetDualsol(consdata->row);
            if( SCIPisFeasNegative(scip, dualsol) )
            {
               if( !SCIPisInfinity(scip, consdata->rhs) )
               {
                  SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                        consdata->vals, +1.0, consdata->rhs, sol, cutoff, ncuts) );
               }
            }
            else if( SCIPisFeasPositive(scip, dualsol) )
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) )
               {
                  SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                        consdata->vals, -1.0, -consdata->lhs, sol, cutoff, ncuts) );
               }
            }
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                  consdata->vals, +1.0, consdata->rhs, sol, cutoff, ncuts) );
         }
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPseparateRelaxedKnapsack(scip, cons, NULL, consdata->nvars, consdata->vars,
                  consdata->vals, -1.0, -consdata->lhs, sol, cutoff, ncuts) );
         }
      }
   }

   if( *ncuts > oldncuts )
   {
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** propagation method for linear constraints */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool             tightenbounds,      /**< should the variable's bounds be tightened? */
   SCIP_Bool             rangedrowpropagation,/**< should ranged row propagation be performed? */
   SCIP_Real             maxeasyactivitydelta,/**< maximum activity delta to run easy propagation on linear constraint */
   SCIP_Bool             sortvars,           /**< should variable sorting for faster propagation be used? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Bool minactisrelax;
   SCIP_Bool maxactisrelax;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   /*SCIPdebugMsg(scip, "propagating linear constraint <%s>\n", SCIPconsGetName(cons));*/

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->eventdata == NULL )
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPconsGetHdlr(cons);
      assert(conshdlr != NULL);

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* catch bound change events of variables */
      SCIP_CALL( consCatchAllEvents(scip, cons, conshdlrdata->eventhdlr) );
      assert(consdata->eventdata != NULL);
   }

   *cutoff = FALSE;

   /* we can only infer activity bounds of the linear constraint, if it is not modifiable */
   if( !SCIPconsIsModifiable(cons) )
   {
      /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
      if( !SCIPinRepropagation(scip) )
      {
         SCIP_CALL( SCIPincConsAge(scip, cons) );
      }

      /* tighten the variable's bounds */
      if( tightenbounds )
      {
         int oldnchgbds;

         oldnchgbds = *nchgbds;

         SCIP_CALL( tightenBounds(scip, cons, maxeasyactivitydelta, sortvars, cutoff, nchgbds) );

         if( *nchgbds > oldnchgbds )
         {
            SCIP_CALL( SCIPresetConsAge(scip, cons) );
         }
      }

      /* propagate ranged rows */
      if( rangedrowpropagation && tightenbounds && !(*cutoff) )
      {
         int nfixedvars;
         int naddconss;
         /* cppcheck-suppress unassignedVariable */
         int oldnchgbds;

         nfixedvars = 0;
         naddconss = 0;
         SCIPdebug( oldnchgbds = *nchgbds; )

         SCIP_CALL( rangedRowPropagation(scip, cons, cutoff, &nfixedvars, nchgbds, &naddconss) );

         if( *cutoff )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is infeasible\n", SCIPconsGetName(cons));
         }
         else
         {
            SCIPdebugMsg(scip, "linear constraint <%s> found %d bound changes and %d fixings\n", SCIPconsGetName(cons), *nchgbds - oldnchgbds, nfixedvars);
         }

         if( nfixedvars > 0 )
            *nchgbds += 2*nfixedvars;
      }

      /* check constraint for infeasibility and redundancy */
      if( !(*cutoff) )
      {
         consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

         if( SCIPisFeasGT(scip, minactivity, consdata->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is infeasible (rhs): activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);

            /* analyze conflict */
            SCIP_CALL( analyzeConflict(scip, cons, TRUE) );

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is infeasible (lhs): activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);

            /* analyze conflict */
            SCIP_CALL( analyzeConflict(scip, cons, FALSE) );

            SCIP_CALL( SCIPresetConsAge(scip, cons) );
            *cutoff = TRUE;
         }
         else if( SCIPisGE(scip, minactivity, consdata->lhs) && SCIPisLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);

            /* remove the constraint locally unless it has become empty, in which case it is removed globally */
            if( consdata->nvars > 0 )
               SCIP_CALL( SCIPdelConsLocal(scip, cons) );
            else
               SCIP_CALL( SCIPdelCons(scip, cons) );
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Presolving methods
 */

/** converts all variables with fixed domain into FIXED variables */
static
SCIP_RETCODE fixVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars          /**< pointer to count the total number of fixed variables */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_VARSTATUS varstatus;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool fixed;
   SCIP_Bool infeasible;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( v = 0; v < consdata->nvars; ++v )
   {
      assert(consdata->vars != NULL);
      var = consdata->vars[v];
      varstatus = SCIPvarGetStatus(var);

      if( varstatus != SCIP_VARSTATUS_FIXED )
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
         if( SCIPisEQ(scip, lb, ub) )
         {
            SCIP_Real fixval;

            fixval = SCIPselectSimpleValue(lb, ub, MAXDNOM);
            SCIPdebugMsg(scip, "converting variable <%s> with fixed bounds [%.15g,%.15g] into fixed variable fixed at %.15g\n",
               SCIPvarGetName(var), lb, ub, fixval);
            SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, " -> infeasible fixing\n");
               *cutoff = TRUE;
               return SCIP_OKAY;
            }
            if( fixed )
               (*nfixedvars)++;
         }
      }
   }

   SCIP_CALL( applyFixings(scip, cons, &infeasible) );

   if( infeasible )
   {
      SCIPdebugMsg(scip, " -> infeasible fixing\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   assert(consdata->removedfixings);

   return SCIP_OKAY;
}

#define MAX_CLIQUE_NONZEROS_PER_CONS 1000000

/** extracts cliques of the constraint and adds them to SCIP
 *
 *  The following clique extraction mechanism are implemeneted
 *
 *  1. collect binary variables and sort them in non increasing order, then
 *
 *     a) if the constraint has a finite right hand side and the negative infinity counters for the minactivity are zero
 *        then add the variables as a clique for which all successive pairs of coefficients fullfill the following
 *        condition
 *
 *            minactivity + vals[i] + vals[i+1] > rhs
 *
 *        and also add the binary to binary implication also for non-successive variables for which the same argument
 *        holds
 *
 *            minactivity + vals[i] + vals[j] > rhs
 *
 *        e.g. 5.3 x1 + 3.6 x2 + 3.3 x3 + 2.1 x4 <= 5.5 (all x are binary) would lead to the clique (x1, x2, x3) and the
 *             binary to binary implications x1 = 1 => x4 = 0 and x2 = 1 => x4 = 0
 *
 *     b) if the constraint has a finite left hand side and the positive infinity counters for the maxactivity are zero
 *        then add the variables as a clique for which all successive pairs of coefficients fullfill the follwoing
 *        condition
 *
 *            maxactivity + vals[i] + vals[i-1] < lhs
 *
 *        and also add the binary to binary implication also for non-successive variables for which the same argument
 *        holds
 *
 *            maxactivity + vals[i] + vals[j] < lhs
 *
 *        e.g. you could multiply the above example by -1
 *
 *     c) the constraint has a finite right hand side and a finite minactivity then add the variables as a negated
 *        clique(clique on the negated variables) for which all successive pairs of coefficients fullfill the following
 *        condition
 *
 *            minactivity - vals[i] - vals[i-1] > rhs
 *
 *        and also add the binary to binary implication also for non-successive variables for which the
 *        same argument holds
 *
 *            minactivity - vals[i] - vals[j] > rhs
 *
 *        e.g. -4 x1 -3 x2 - 2 x3 + 2 x4 <= -4 would lead to the (negated) clique (~x1, ~x2) and the binary to binary
 *             implication x1 = 0 => x3 = 1
 *
 *     d) the constraint has a finite left hand side and a finite maxactivity then add the variables as a negated
 *        clique(clique on the negated variables) for which all successive pairs of coefficients fullfill the following
 *        condition
 *
 *            maxactivity - vals[i] - vals[i+1] < lhs
 *
 *        and also add the binary to binary implication also for non-successive variables for which the same argument
 *        holds
 *
 *            maxactivity - vals[i] - vals[j] < lhs
 *
 *        e.g. you could multiply the above example by -1
 *
 *  2. if the linear constraint represents a set-packing or set-partitioning constraint, the whole constraint is added
 *     as clique, (this part is done at the end of the method)
 *
 */
static
SCIP_RETCODE extractCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             maxeasyactivitydelta,/**< maximum activity delta to run easy propagation on linear constraint */
   SCIP_Bool             sortvars,           /**< should variables be used in sorted order? */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  nchgbds,            /**< pointer to count the total number of tightened bounds */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if a cutoff was found */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_CONSDATA* consdata;
   SCIP_Bool lhsclique;
   SCIP_Bool rhsclique;
   SCIP_Bool finitelhs;
   SCIP_Bool finiterhs;
   SCIP_Bool finiteminact;
   SCIP_Bool finitemaxact;
   SCIP_Bool finitenegminact;
   SCIP_Bool finitenegmaxact;
   SCIP_Bool finiteposminact;
   SCIP_Bool finiteposmaxact;
   SCIP_Bool infeasible;
   SCIP_Bool stopped;
   int cliquenonzerosadded;
   int v;
   int i;
   int nposcoefs;
   int nnegcoefs;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);
   assert(nchgbds != NULL);
   assert(cutoff != NULL);
   assert(!SCIPconsIsDeleted(cons));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nvars < 2 )
      return SCIP_OKAY;

   /* add implications if posibble
    *
    * for now we only add binary to non-binary implications, and this is only done for the binary variable with the
    * maximal absolute contribution and also only if this variable would force all other variables to their bound
    * corresponding to the global minimal activity of the constraint
    */
   if( !consdata->implsadded )
   {
      /* sort variables by variable type */
      SCIP_CALL( consdataSort(scip, consdata) );

      /* @todo we might extract implications/cliques if SCIPvarIsBinary() variables exist and we have integer variables
       *       up front, might change sorting correspondingly
       */
      /* fast abort if no binaries exist */
      if( !SCIPvarIsBinary(consdata->vars[0]) )
         return SCIP_OKAY;

      nvars = consdata->nvars;
      vars = consdata->vars;
      vals = consdata->vals;

      /* recompute activities if needed */
      if( !consdata->validactivities )
         consdataCalcActivities(scip, consdata);
      assert(consdata->validactivities);

      finitelhs = !SCIPisInfinity(scip, -consdata->lhs);
      finiterhs = !SCIPisInfinity(scip, consdata->rhs);
      finitenegminact = (consdata->glbminactivityneginf == 0 && consdata->glbminactivityneghuge == 0);
      finitenegmaxact = (consdata->glbmaxactivityneginf == 0 && consdata->maxactivityneghuge == 0);
      finiteposminact = (consdata->glbminactivityposinf == 0 && consdata->glbminactivityposhuge == 0);
      finiteposmaxact = (consdata->glbmaxactivityposinf == 0 && consdata->glbmaxactivityposhuge == 0);
      finiteminact = (finitenegminact && finiteposminact);
      finitemaxact = (finitenegmaxact && finiteposmaxact);

      if( (finiterhs || finitelhs) && (finitenegminact || finiteposminact || finitenegmaxact || finiteposmaxact) )
      {
         SCIP_Real maxabscontrib = -1.0;
         SCIP_Bool posval = FALSE;
         SCIP_Bool allbinary = TRUE;
         int oldnchgbds = *nchgbds;
         int nbdchgs = 0;
         int nimpls = 0;
         int position = -1;

         /* we need a valid minimal/maximal activity to add cliques */
         if( (finitenegminact || finiteposminact) && !consdata->validglbminact )
         {
            consdataRecomputeGlbMinactivity(scip, consdata);
            assert(consdata->validglbminact);
         }

         if( (finitenegmaxact || finiteposmaxact) && !consdata->validglbmaxact )
         {
            consdataRecomputeGlbMaxactivity(scip, consdata);
            assert(consdata->validglbmaxact);
         }
         assert(consdata->validglbminact || consdata->validglbmaxact);

         /* @todo extend this to local/constraint probing */

         /* determine maximal contribution to the activity */
         for( v = nvars - 1; v >= 0; --v )
         {
            if( SCIPvarIsBinary(vars[v]) )
            {
               if( vals[v] > 0 )
               {
                  SCIP_Real value = vals[v] * SCIPvarGetUbGlobal(vars[v]);

                  if( value > maxabscontrib )
                  {
                     maxabscontrib = value;
                     position = v;
                     posval = TRUE;
                  }
               }
               else
               {
                  SCIP_Real value = vals[v] * SCIPvarGetLbGlobal(vars[v]);

                  value = REALABS(value);

                  if( value > maxabscontrib )
                  {
                     maxabscontrib = value;
                     position = v;
                     posval = FALSE;
                  }
               }
            }
            else
               allbinary = FALSE;
         }
         assert(0 <= position && position < nvars);

         if( !SCIPisEQ(scip, maxabscontrib, 1.0) && !allbinary )
         {
            /* if the right hand side and the minimal activity are finite and changing the variable with the biggest
             * influence to their bound forces all other variables to be at their minimal contribution, we can add these
             * implications
             */
            if( finiterhs && finiteminact && SCIPisEQ(scip, consdata->glbminactivity, consdata->rhs - maxabscontrib) )
            {
               for( v = nvars - 1; v >= 0; --v )
               {
                  /* binary to binary implications will be collected when extrating cliques */
                  if( !SCIPvarIsBinary(vars[v]) )
                  {
                     if( v != position )
                     {
                        if( vals[v] > 0 )
                        {
                           /* add implications */
                           SCIP_CALL( SCIPaddVarImplication(scip, vars[position], posval, vars[v], SCIP_BOUNDTYPE_UPPER, SCIPvarGetLbGlobal(vars[v]), &infeasible, &nbdchgs) );
                           ++nimpls;
                           *nchgbds += nbdchgs;
                        }
                        else
                        {
                           /* add implications */
                           SCIP_CALL( SCIPaddVarImplication(scip, vars[position], posval, vars[v], SCIP_BOUNDTYPE_LOWER, SCIPvarGetUbGlobal(vars[v]), &infeasible, &nbdchgs) );
                           ++nimpls;
                           *nchgbds += nbdchgs;
                        }

                        if( infeasible )
                        {
                           *cutoff = TRUE;
                           break;
                        }
                     }
                  }
                  /* stop when reaching a 'real' binary variable because the variables are sorted after their type */
                  else if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
                     break;
               }
            }

            /* if the left hand side and the maximal activity are finite and changing the variable with the biggest
             * influence to their bound forces all other variables to be at their minimal contribution, we can add these
             * implications
             */
            if( finitelhs && finitemaxact && SCIPisEQ(scip, consdata->glbmaxactivity, consdata->lhs - maxabscontrib) )
            {
               for( v = nvars - 1; v >= 0; --v )
               {
                  /* binary to binary implications will be collected when extrating cliques */
                  if( !SCIPvarIsBinary(vars[v]) )
                  {
                     if( v != position )
                     {
                        if( vals[v] > 0 )
                        {
                           /* add implications */
                           SCIP_CALL( SCIPaddVarImplication(scip, vars[position], posval, vars[v], SCIP_BOUNDTYPE_LOWER, SCIPvarGetUbGlobal(vars[v]), &infeasible, &nbdchgs) );
                           ++nimpls;
                           *nchgbds += nbdchgs;
                        }
                        else
                        {
                           /* add implications */
                           SCIP_CALL( SCIPaddVarImplication(scip, vars[position], posval, vars[v], SCIP_BOUNDTYPE_UPPER, SCIPvarGetLbGlobal(vars[v]), &infeasible, &nbdchgs) );
                           ++nimpls;
                           *nchgbds += nbdchgs;
                        }

                        if( infeasible )
                        {
                           *cutoff = TRUE;
                           break;
                        }
                     }
                  }
                  /* stop when reaching a 'real' binary variable because the variables are sorted after their type */
                  else if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY )
                     break;
               }
            }

            /* did we find some implications */
            if( nimpls > 0 )
            {
               SCIPdebugMsg(scip, "extracted %d implications from constraint %s which led to %d bound changes, %scutoff detetcted\n", nimpls, SCIPconsGetName(cons), *nchgbds - oldnchgbds, *cutoff ? "" : "no ");

               if( *cutoff )
                  return SCIP_OKAY;

               /* did we find some boundchanges, then we need to remove fixings and tighten the bounds further */
               if( *nchgbds - oldnchgbds > 0 )
               {
                  /* check for fixed variables */
                  SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );
                  if( *cutoff )
                     return SCIP_OKAY;

                  /* tighten variable's bounds */
                  SCIP_CALL( tightenBounds(scip, cons, maxeasyactivitydelta, sortvars, cutoff, nchgbds) );
                  if( *cutoff )
                     return SCIP_OKAY;

                  /* check for fixed variables */
                  SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );
                  if( *cutoff )
                     return SCIP_OKAY;
               }
            }
         }
      }

      consdata->implsadded = TRUE;
   }

   /* check if we already added the cliques of this constraint */
   if( consdata->cliquesadded )
      return SCIP_OKAY;

   consdata->cliquesadded = TRUE;
   cliquenonzerosadded = 0;
   stopped = FALSE;

   /* sort variables by variable type */
   SCIP_CALL( consdataSort(scip, consdata) );

   nvars = consdata->nvars;
   vars = consdata->vars;
   vals = consdata->vals;

   /**@todo extract more cliques, implications and variable bounds from linear constraints */

   /* recompute activities if needed */
   if( !consdata->validactivities )
      consdataCalcActivities(scip, consdata);
   assert(consdata->validactivities);

   finitelhs = !SCIPisInfinity(scip, -consdata->lhs);
   finiterhs = !SCIPisInfinity(scip, consdata->rhs);
   finitenegminact = (consdata->glbminactivityneginf == 0 && consdata->glbminactivityneghuge == 0);
   finitenegmaxact = (consdata->glbmaxactivityneginf == 0 && consdata->maxactivityneghuge == 0);
   finiteposminact = (consdata->glbminactivityposinf == 0 && consdata->glbminactivityposhuge == 0);
   finiteposmaxact = (consdata->glbmaxactivityposinf == 0 && consdata->glbmaxactivityposhuge == 0);
   finiteminact = (finitenegminact && finiteposminact);
   finitemaxact = (finitenegmaxact && finiteposmaxact);

   /* 1. we wheck whether some variables do not fit together into this constraint and add the corresponding clique
    *    information
    */
   if( (finiterhs || finitelhs) && (finitenegminact || finiteposminact || finitenegmaxact || finiteposmaxact) )
   {
      SCIP_VAR** binvars;
      SCIP_Real* binvarvals;
      int nposbinvars = 0;
      int nnegbinvars = 0;
      int allonebinary = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &binvarvals, nvars) );

      /* collect binary variables */
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarIsBinary(vars[i]) )
         {
            assert(!SCIPisZero(scip, vals[i]));

            if( SCIPisEQ(scip, REALABS(vals[i]), 1.0) )
               ++allonebinary;

            binvars[nposbinvars + nnegbinvars] = vars[i];
            binvarvals[nposbinvars + nnegbinvars] = vals[i];

            if( SCIPisPositive(scip, vals[i]) )
               ++nposbinvars;
            else
               ++nnegbinvars;

            assert(nposbinvars + nnegbinvars <= nvars);
         }
         /* stop searching for binary variables, because the constraint data is sorted */
         else if( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS )
            break;
      }
      assert(nposbinvars + nnegbinvars <= nvars);

      /* setppc constraints will be handled later; we need at least two binary variables with same sign to extract
       * cliques
       */
      if( allonebinary < nvars && (nposbinvars >= 2 || nnegbinvars >= 2) )
      {
         SCIP_Real threshold;
         int oldnchgbds = *nchgbds;
         int nbdchgs;
         int jstart;
         int j;

         /* we need a valid minimal/maximal activity to add cliques */
         if( (finitenegminact || finiteposminact) && !consdata->validglbminact )
         {
            consdataRecomputeGlbMinactivity(scip, consdata);
            assert(consdata->validglbminact);
         }

         if( (finitenegmaxact || finiteposmaxact) && !consdata->validglbmaxact )
         {
            consdataRecomputeGlbMaxactivity(scip, consdata);
            assert(consdata->validglbmaxact);
         }
         assert(consdata->validglbminact || consdata->validglbmaxact);

         /* sort coefficients non-increasing to be faster in the clique search */
         SCIPsortDownRealPtr(binvarvals, (void**) binvars, nposbinvars + nnegbinvars);

         /* case a) */
         if( finiterhs && finitenegminact && nposbinvars >= 2 )
         {
            /* compute value that needs to be exceeded */
            threshold = consdata->rhs - consdata->glbminactivity;

            i = 0;
            j = i + 1;
#if 0 /* assertion should only holds when constraints were fully propagated and boundstightened */
            /* check that it is possible to choose binvar[i], otherwise it should have been fixed to zero */
            assert(SCIPisFeasLE(scip, binvarvals[i], threshold));
#endif
            /* check if at least two variables are in a clique */
            if( SCIPisFeasGT(scip, binvarvals[i] + binvarvals[j], threshold) )
            {
               ++j;
               /* check for extending the clique */
               while( j < nposbinvars )
               {
                  if( !SCIPisFeasGT(scip, binvarvals[j-1] + binvarvals[j], threshold) )
                     break;
                  ++j;
               }
               assert(j >= 2);

               /* add clique with at least two variables */
               SCIP_CALL( SCIPaddClique(scip, &(binvars[i]), NULL, j - i, FALSE, &infeasible, &nbdchgs) );

               if( infeasible )
                  *cutoff = TRUE;

               *nchgbds += nbdchgs;

               cliquenonzerosadded += j;
               if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                  stopped = TRUE;

               /* exchange the last variable in the clique if possible and add all new ones */
               if( !stopped && !(*cutoff) && j < nposbinvars )
               {
                  SCIP_VAR** clqvars;
                  int lastfit = j - 2;
                  assert(lastfit >= i);

                  /* copy all 'main'-clique variables */
                  SCIP_CALL( SCIPduplicateBufferArray(scip, &clqvars, &(binvars[i]), j - i) );

                  /* iterate up to the end with j and up to the front with lastfit, and check for different cliques */
                  while( lastfit >= i && j < nposbinvars )
                  {
                     /* check if two variables are in a clique */
                     if( SCIPisFeasGT(scip, binvarvals[lastfit] + binvarvals[j], threshold) )
                     {
                        clqvars[lastfit + 1] = binvars[j];

                        /* add clique with at least two variables */
                        SCIP_CALL( SCIPaddClique(scip, clqvars, NULL, lastfit - i + 2, FALSE, &infeasible, &nbdchgs) );

                        if( infeasible )
                        {
                           *cutoff = TRUE;
                           break;
                        }

                        *nchgbds += nbdchgs;

                        cliquenonzerosadded += (lastfit - i + 2);
                        if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                        {
                           stopped = TRUE;
                           break;
                        }

                        ++j;
                     }
                     else
                        --lastfit;
                  }

                  SCIPfreeBufferArray(scip, &clqvars);
               }
            }
         }

         /* did we find some boundchanges, then we need to remove fixings and tighten the bounds further */
         if( !stopped && !*cutoff && *nchgbds - oldnchgbds > 0 )
         {
            /* check for fixed variables */
            SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );

            if( !*cutoff )
            {
               /* tighten variable's bounds */
               SCIP_CALL( tightenBounds(scip, cons, maxeasyactivitydelta, sortvars, cutoff, nchgbds) );

               if( !*cutoff )
               {
                  /* check for fixed variables */
                  SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );

                  if( !*cutoff )
                  {
                     /* sort variables by variable type */
                     SCIP_CALL( consdataSort(scip, consdata) );

                     /* recompute activities if needed */
                     if( !consdata->validactivities )
                        consdataCalcActivities(scip, consdata);
                     assert(consdata->validactivities);

                     nvars = consdata->nvars;
                     vars = consdata->vars;
                     vals = consdata->vals;
                     nposbinvars = 0;
                     nnegbinvars = 0;
                     allonebinary = 0;

                     /* update binary variables */
                     for( i = 0; i < nvars; ++i )
                     {
                        if( SCIPvarIsBinary(vars[i]) )
                        {
                           assert(!SCIPisZero(scip, vals[i]));

                           if( SCIPisEQ(scip, REALABS(vals[i]), 1.0) )
                              ++allonebinary;

                           binvars[nposbinvars + nnegbinvars] = vars[i];
                           binvarvals[nposbinvars + nnegbinvars] = vals[i];

                           if( SCIPisPositive(scip, vals[i]) )
                              ++nposbinvars;
                           else
                              ++nnegbinvars;

                           assert(nposbinvars + nnegbinvars <= nvars);
                        }
                        /* stop searching for binary variables, because the constraint data is sorted */
                        else if( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS )
                           break;
                     }
                     assert(nposbinvars + nnegbinvars <= nvars);
                  }
               }
            }

            oldnchgbds = *nchgbds;
         }

         /* case b) */
         if( !stopped && !(*cutoff) && finitelhs && finiteposmaxact && nnegbinvars >= 2 )
         {
            /* compute value that needs to be deceeded */
            threshold = consdata->lhs - consdata->glbmaxactivity;

            i = nposbinvars + nnegbinvars - 1;
            j = i - 1;
#if 0 /* assertion should only holds when constraints were fully propagated and boundstightened */
            /* check that it is possible to choose binvar[i], otherwise it should have been fixed to zero */
            assert(SCIPisFeasGE(scip, binvarvals[i], threshold));
#endif
            /* check if two variables are in a clique */
            if( SCIPisFeasLT(scip, binvarvals[i] + binvarvals[j], threshold) )
            {
               --j;
               /* check for extending the clique */
               while( j >= nposbinvars )
               {
                  if( !SCIPisFeasLT(scip, binvarvals[j+1] + binvarvals[j], threshold) )
                     break;
                  --j;
               }
               jstart = j;

               assert(i - j >= 2);
               /* add clique with at least two variables */
               SCIP_CALL( SCIPaddClique(scip, &(binvars[j+1]), NULL, i - j, FALSE, &infeasible, &nbdchgs) );

               if( infeasible )
                  *cutoff = TRUE;

               *nchgbds += nbdchgs;

               cliquenonzerosadded += (i - j);
               if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                  stopped = TRUE;

               /* exchange the last variable in the clique if possible and add all new ones */
               if( !stopped && !(*cutoff) && jstart >= nposbinvars )
               {
                  SCIP_VAR** clqvars;
                  int lastfit = jstart + 1;
                  assert(lastfit < i);

                  /* copy all 'main'-clique variables */
                  SCIP_CALL( SCIPduplicateBufferArray(scip, &clqvars, &(binvars[lastfit]), i - j) );
                  ++lastfit;

                  /* iterate up to the front with j and up to the end with lastfit, and check for different cliques */
                  while( lastfit <= i && j >= nposbinvars )
                  {
                     /* check if two variables are in a clique */
                     if( SCIPisFeasLT(scip, binvarvals[lastfit] + binvarvals[j], threshold) )
                     {
                        assert(lastfit - jstart - 2 >= 0 && lastfit - jstart - 2 < i);
                        clqvars[lastfit - jstart - 2] = binvars[j];

                        assert(i - lastfit + 2 >= 2);
                        /* add clique with at least two variables */
                        SCIP_CALL( SCIPaddClique(scip, &(clqvars[lastfit - jstart - 2]), NULL, i - lastfit + 2, FALSE, &infeasible, &nbdchgs) );

                        if( infeasible )
                        {
                           *cutoff = TRUE;
                           break;
                        }

                        *nchgbds += nbdchgs;

                        cliquenonzerosadded += (i - lastfit + 2);
                        if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                        {
                           stopped = TRUE;
                           break;
                        }

                        --j;
                     }
                     else
                        ++lastfit;
                  }

                  SCIPfreeBufferArray(scip, &clqvars);
               }
            }
         }

         /* did we find some boundchanges, then we need to remove fixings and tighten the bounds further */
         if( !stopped && !*cutoff && *nchgbds - oldnchgbds > 0 )
         {
            /* check for fixed variables */
            SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );

            if( !*cutoff )
            {
               /* tighten variable's bounds */
               SCIP_CALL( tightenBounds(scip, cons, maxeasyactivitydelta, sortvars, cutoff, nchgbds) );

               if( !*cutoff )
               {
                  /* check for fixed variables */
                  SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );

                  if( !*cutoff )
                  {
                     /* sort variables by variable type */
                     SCIP_CALL( consdataSort(scip, consdata) );

                     /* recompute activities if needed */
                     if( !consdata->validactivities )
                        consdataCalcActivities(scip, consdata);
                     assert(consdata->validactivities);

                     nvars = consdata->nvars;
                     vars = consdata->vars;
                     vals = consdata->vals;
                     nposbinvars = 0;
                     nnegbinvars = 0;
                     allonebinary = 0;

                     /* update binary variables */
                     for( i = 0; i < nvars; ++i )
                     {
                        if( SCIPvarIsBinary(vars[i]) )
                        {
                           assert(!SCIPisZero(scip, vals[i]));

                           if( SCIPisEQ(scip, REALABS(vals[i]), 1.0) )
                              ++allonebinary;

                           binvars[nposbinvars + nnegbinvars] = vars[i];
                           binvarvals[nposbinvars + nnegbinvars] = vals[i];

                           if( SCIPisPositive(scip, vals[i]) )
                              ++nposbinvars;
                           else
                              ++nnegbinvars;

                           assert(nposbinvars + nnegbinvars <= nvars);
                        }
                        /* stop searching for binary variables, because the constraint data is sorted */
                        else if( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS )
                           break;
                     }
                     assert(nposbinvars + nnegbinvars <= nvars);
                  }
               }
            }

            oldnchgbds = *nchgbds;
         }

         /* case c) */
         if( !(*cutoff) && finiterhs && finiteminact && nnegbinvars >= 2 )
         {
            SCIP_Bool* values;

            /* initialize clique values array for adding a negated clique */
            SCIP_CALL( SCIPallocBufferArray(scip, &values, nnegbinvars) );
            BMSclearMemoryArray(values, nnegbinvars);

            /* compute value that needs to be exceeded */
            threshold = consdata->rhs - consdata->glbminactivity;

            i = nposbinvars + nnegbinvars - 1;
            j = i - 1;

#if 0 /* assertion should only holds when constraints were fully propagated and boundstightened */
            /* check if the variable should not have already been fixed to one */
            assert(!SCIPisFeasGT(scip, binvarvals[i], threshold));
#endif

            if( SCIPisFeasGT(scip, -binvarvals[i] - binvarvals[j], threshold) )
            {
               --j;
               /* check for extending the clique */
               while( j >= nposbinvars )
               {
                  if( !SCIPisFeasGT(scip, -binvarvals[j+1] - binvarvals[j], threshold) )
                     break;
                  --j;
               }
               jstart = j;

               assert(i - j >= 2);
               /* add negated clique with at least two variables */
               SCIP_CALL( SCIPaddClique(scip, &(binvars[j+1]), values, i - j, FALSE, &infeasible, &nbdchgs) );

               if( infeasible )
                  *cutoff = TRUE;

               *nchgbds += nbdchgs;

               cliquenonzerosadded += (i - j);
               if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                  stopped = TRUE;

               /* exchange the last variable in the clique if possible and add all new ones */
               if( !stopped && !(*cutoff) && jstart >= nposbinvars )
               {
                  SCIP_VAR** clqvars;
                  int lastfit = j + 1;
                  assert(lastfit < i);

                  /* copy all 'main'-clique variables */
                  SCIP_CALL( SCIPduplicateBufferArray(scip, &clqvars, &(binvars[lastfit]), i - j) );
                  ++lastfit;

                  /* iterate up to the front with j and up to the end with lastfit, and check for different cliques */
                  while( lastfit <= i && j >= nposbinvars )
                  {
                     /* check if two variables are in a negated clique */
                     if( SCIPisFeasGT(scip, -binvarvals[lastfit] - binvarvals[j], threshold) )
                     {
                        assert(lastfit - jstart - 2 >= 0 && lastfit - jstart - 2 < i);
                        clqvars[lastfit - jstart - 2] = binvars[j];

                        assert(i - lastfit + 2 >= 2);
                        /* add clique with at least two variables */
                        SCIP_CALL( SCIPaddClique(scip, &(clqvars[lastfit - jstart - 2]), values, i - lastfit + 2, FALSE, &infeasible, &nbdchgs) );

                        if( infeasible )
                        {
                           *cutoff = TRUE;
                           break;
                        }

                        *nchgbds += nbdchgs;

                        cliquenonzerosadded += (i - lastfit + 2);
                        if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                        {
                           stopped = TRUE;
                           break;
                        }

                        --j;
                     }
                     else
                        ++lastfit;
                  }

                  SCIPfreeBufferArray(scip, &clqvars);
               }
            }

            SCIPfreeBufferArray(scip, &values);
         }

         /* did we find some boundchanges, then we need to remove fixings and tighten the bounds further */
         if( !stopped && !*cutoff && *nchgbds - oldnchgbds > 0 )
         {
            /* check for fixed variables */
            SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );

            if( !*cutoff )
            {
               /* tighten variable's bounds */
               SCIP_CALL( tightenBounds(scip, cons, maxeasyactivitydelta, sortvars, cutoff, nchgbds) );

               if( !*cutoff )
               {
                  /* check for fixed variables */
                  SCIP_CALL( fixVariables(scip, cons, cutoff, nfixedvars) );

                  if( !*cutoff )
                  {
                     /* sort variables by variable type */
                     SCIP_CALL( consdataSort(scip, consdata) );

                     /* recompute activities if needed */
                     if( !consdata->validactivities )
                        consdataCalcActivities(scip, consdata);
                     assert(consdata->validactivities);

                     nvars = consdata->nvars;
                     vars = consdata->vars;
                     vals = consdata->vals;
                     nposbinvars = 0;
                     nnegbinvars = 0;
                     allonebinary = 0;

                     /* update binary variables */
                     for( i = 0; i < nvars; ++i )
                     {
                        if( SCIPvarIsBinary(vars[i]) )
                        {
                           assert(!SCIPisZero(scip, vals[i]));

                           if( SCIPisEQ(scip, REALABS(vals[i]), 1.0) )
                              ++allonebinary;

                           binvars[nposbinvars + nnegbinvars] = vars[i];
                           binvarvals[nposbinvars + nnegbinvars] = vals[i];

                           if( SCIPisPositive(scip, vals[i]) )
                              ++nposbinvars;
                           else
                              ++nnegbinvars;

                           assert(nposbinvars + nnegbinvars <= nvars);
                        }
                        /* stop searching for binary variables, because the constraint data is sorted */
                        else if( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS )
                           break;
                     }
                     assert(nposbinvars + nnegbinvars <= nvars);
                  }
               }
            }
         }

         /* case d) */
         if( !stopped && !(*cutoff) && finitelhs && finitemaxact && nposbinvars >= 2 )
         {
            SCIP_Bool* values;

            /* initialize clique values array for adding a negated clique */
            SCIP_CALL( SCIPallocBufferArray(scip, &values, nposbinvars) );
            BMSclearMemoryArray(values, nposbinvars);

            /* compute value that needs to be exceeded */
            threshold = consdata->lhs - consdata->glbmaxactivity;

            i = 0;
            j = i + 1;

#if 0 /* assertion should only holds when constraints were fully propagated and boundstightened */
            /* check if the variable should not have already been fixed to one */
            assert(!SCIPisFeasLT(scip, -binvarvals[i], threshold));
#endif

            if( SCIPisFeasLT(scip, -binvarvals[i] - binvarvals[j], threshold) )
            {
               ++j;
               /* check for extending the clique */
               while( j < nposbinvars )
               {
                  if( !SCIPisFeasLT(scip, -binvarvals[j-1] - binvarvals[j], threshold) )
                     break;
                  ++j;
               }
               assert(j >= 2);

               /* add negated clique with at least two variables */
               SCIP_CALL( SCIPaddClique(scip, &(binvars[i]), values, j - i, FALSE, &infeasible, &nbdchgs) );

               if( infeasible )
                  *cutoff = TRUE;

               *nchgbds += nbdchgs;

               cliquenonzerosadded += j;
               if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                  stopped = TRUE;

               /* exchange the last variable in the clique if possible and add all new ones */
               if( !stopped && !(*cutoff) && j < nposbinvars )
               {
                  SCIP_VAR** clqvars;
                  int lastfit = j - 2;
                  assert(lastfit >= i);

                  /* copy all 'main'-clique variables */
                  SCIP_CALL( SCIPduplicateBufferArray(scip, &clqvars, &(binvars[i]), j - i) );

                  /* iterate up to the end with j and up to the front with lastfit, and check for different cliques */
                  while( lastfit >= i && j < nposbinvars )
                  {
                     /* check if two variables are in a negated clique */
                     if( SCIPisFeasLT(scip, -binvarvals[lastfit] - binvarvals[j], threshold) )
                     {
                        clqvars[lastfit + 1] = binvars[j];

                        /* add clique with at least two variables */
                        SCIP_CALL( SCIPaddClique(scip, clqvars, values, lastfit - i + 2, FALSE, &infeasible, &nbdchgs) );

                        if( infeasible )
                        {
                           *cutoff = TRUE;
                           break;
                        }

                        *nchgbds += nbdchgs;

                        cliquenonzerosadded += (lastfit - i + 2);
                        if( cliquenonzerosadded >= MAX_CLIQUE_NONZEROS_PER_CONS )
                           break;

                        ++j;
                     }
                     else
                        --lastfit;
                  }

                  SCIPfreeBufferArray(scip, &clqvars);
               }
            }

            SCIPfreeBufferArray(scip, &values);
         }
      }

      SCIPfreeBufferArray(scip, &binvarvals);
      SCIPfreeBufferArray(scip, &binvars);

      if( *cutoff )
         return SCIP_OKAY;
   }

   /* 2. we only check if the constraint is a set packing / partitioning constraint */

   /* check if all variables are binary, if the coefficients are +1 or -1, and if the right hand side is equal
    * to 1 - number of negative coefficients, or if the left hand side is equal to number of positive coefficients - 1
    */
   nposcoefs = 0;
   nnegcoefs = 0;
   for( i = 0; i < nvars; ++i )
   {
      if( !SCIPvarIsBinary(vars[i]) )
         return SCIP_OKAY;
      else if( SCIPisEQ(scip, vals[i], +1.0) )
         nposcoefs++;
      else if( SCIPisEQ(scip, vals[i], -1.0) )
         nnegcoefs++;
      else
         return SCIP_OKAY;
   }

   lhsclique = SCIPisEQ(scip, consdata->lhs, (SCIP_Real)nposcoefs - 1.0);
   rhsclique = SCIPisEQ(scip, consdata->rhs, 1.0 - (SCIP_Real)nnegcoefs);

   if( lhsclique || rhsclique )
   {
      SCIP_Bool* values;
      int nbdchgs;

      SCIPdebugMsg(scip, "linear constraint <%s>: adding clique with %d vars (%d pos, %d neg)\n",
         SCIPconsGetName(cons), nvars, nposcoefs, nnegcoefs);
      SCIP_CALL( SCIPallocBufferArray(scip, &values, nvars) );

      for( i = 0; i < nvars; ++i )
         values[i] = (rhsclique == (vals[i] > 0.0));

      SCIP_CALL( SCIPaddClique(scip, vars, values, nvars, SCIPisEQ(scip, consdata->lhs, consdata->rhs), &infeasible, &nbdchgs) );

      if( infeasible )
         *cutoff = TRUE;

      *nchgbds += nbdchgs;
      SCIPfreeBufferArray(scip, &values);
   }

   return SCIP_OKAY;
}

/** tightens left and right hand side of constraint due to integrality */
static
SCIP_RETCODE tightenSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool integral;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgsides != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisIntegral(scip, consdata->lhs) || !SCIPisIntegral(scip, consdata->rhs) )
   {
      integral = TRUE;
      for( i = 0; i < consdata->nvars && integral; ++i )
      {
         integral = SCIPisIntegral(scip, consdata->vals[i])
            && (SCIPvarGetType(consdata->vars[i]) != SCIP_VARTYPE_CONTINUOUS);
      }
      if( integral )
      {
         SCIPdebugMsg(scip, "linear constraint <%s>: make sides integral: sides=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
         if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisIntegral(scip, consdata->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons, SCIPfeasCeil(scip, consdata->lhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisIntegral(scip, consdata->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons, SCIPfeasFloor(scip, consdata->rhs)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         SCIPdebugMsg(scip, "linear constraint <%s>: new integral sides: sides=[%.15g,%.15g]\n",
            SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
      }
   }

   return SCIP_OKAY;
}

/** tightens coefficients of binary, integer, and implicit integer variables due to activity bounds in presolving:
 *  given an inequality  lhs <= a*x + ai*xi <= rhs, with a non-continuous variable  li <= xi <= ui
 *  let minact := min{a*x + ai*xi}, maxact := max{a*x + ai*xi}
 *  (i) ai >= 0:
 *      if  minact + ai >= lhs  and  maxact - ai <= rhs: (**)
 *       - a deviation from the lower/upper bound of xi would make the left/right hand side redundant
 *       - ai, lhs and rhs can be changed to have the same redundancy effect and the same results for
 *         xi fixed to its bounds, but with a reduced ai and tightened sides to tighten the LP relaxation
 *       - change coefficients:
 *           ai'  := max(lhs - minact, maxact - rhs)
 *           lhs' := lhs - (ai - ai')*li
 *           rhs' := rhs - (ai - ai')*ui
 * (ii) ai < 0:
 *      if  minact - ai >= lhs  and  maxact + ai <= rhs: (***)
 *       - a deviation from the upper/lower bound of xi would make the left/right hand side redundant
 *       - ai, lhs and rhs can be changed to have the same redundancy effect and the same results for
 *         xi fixed to its bounds, but with a reduced ai and tightened sides to tighten the LP relaxation
 *       - change coefficients:
 *           ai'  := min(rhs - maxact, minact - lhs)
 *           lhs' := lhs - (ai - ai')*ui
 *           rhs' := rhs - (ai - ai')*li
 *
 *  We further try to remove redundant variable from the constraint;
 *  Variables which fulfill conditions (**) or (***) are called surely non-redundant variables.
 *  A deviation of only one from their bound makes the lhs/rhs feasible (i.e., redundant), even if all other 
 *  variables are set to their "worst" bound. If all variables which are not surely non-redundant cannot make 
 *  the lhs/rhs redundant, even if they are set to their "best" bound, they can be removed from the constraint.
 *  E.g., for binary variables and an inequality x_1 +x_2 +10y_1 +10y_2 >= 5, setting either of the y_i to one 
 *  suffices to fulfill the inequality, whereas the x_i do not contribute to feasibility and can be removed.
 *
 *  @todo use also some tightening procedures for (knapsack) constraints with non-integer coefficients, see
 *        cons_knapsack.c the following methods detectRedundantVars() and tightenWeights()
 */
static
SCIP_RETCODE consdataTightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgcoefs,          /**< pointer to count total number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of side changes */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real minactivity; /* minimal value w.r.t. the variable's local bounds for the constraint's
                           * activity, ignoring the coefficients contributing with infinite value */
   SCIP_Real maxactivity; /* maximal value w.r.t. the variable's local bounds for the constraint's
                           * activity, ignoring the coefficients contributing with infinite value */
   SCIP_Bool minactisrelax; /* do huge finite values contribute to the minactivity? */
   SCIP_Bool maxactisrelax; /* do huge finite values contribute to the maxactivity? */
   SCIP_Real minleftactivity; /* minimal activity without surely non-redundant variables. */
   SCIP_Real maxleftactivity; /* maximal activity without surely non-redundant variables. */
   SCIP_Real aggrlhs; /* lhs without minimal activity of surely non-redundant variables. */
   SCIP_Real aggrrhs; /* rhs without maximal activity of surely non-redundant variables. */
   SCIP_Real lval; /* candidate for new value arising from considering the left hand side */
   SCIP_Real rval; /* candidate for new value arising from considering the left hand side */
   SCIP_Real val;
   SCIP_Real newval;
   SCIP_Real newlhs;
   SCIP_Real newrhs;
   SCIP_Real lb;
   SCIP_Real ub;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* @todo Is this still needed with automatic recomputation of activities? */
   /* if the maximal coefficient is too large, recompute the activities */
   if( (consdata->validmaxabsval && consdata->maxabsval > MAXVALRECOMP)
      || (consdata->validminabsval && consdata->minabsval < MINVALRECOMP) )
   {
      consdataRecomputeMinactivity(scip, consdata);
      consdataRecomputeMaxactivity(scip, consdata);
   }

   /* get the minimal and maximal activity of the constraint */
   consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

   minleftactivity = 0.0;
   maxleftactivity = 0.0;

   /* try to tighten each coefficient */
   i = 0;
   while( i < consdata->nvars )
   {
      var = consdata->vars[i];

      /* get coefficient and variable's bounds */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      val = consdata->vals[i];
      assert(!SCIPisZero(scip, val));

      /* check sign of coefficient */
      if( val >= 0.0 )
      {
         /* check, if a deviation from lower/upper bound would make lhs/rhs redundant */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS &&
            SCIPisGE(scip, minactivity + val, consdata->lhs) && SCIPisLE(scip, maxactivity - val, consdata->rhs) )
         {
            /* change coefficients:
             *   ai'  := max(lhs - minact, maxact - rhs)
             *   lhs' := lhs - (ai - ai')*li
             *   rhs' := rhs - (ai - ai')*ui
             */

            lval = consdata->lhs - minactivity;
            rval = maxactivity - consdata->rhs;

            /* Try to avoid cancellation, if there are only two variables */
            if( consdata->nvars == 2 )
            {
               SCIP_Real otherval;
               otherval = consdata->vals[1-i];

               if( !SCIPisInfinity(scip, -consdata->lhs) && consdata->minactivityneginf + consdata->minactivityneginf == 0 )
               {
                  lval = consdata->lhs - val*lb;
                  lval -= otherval > 0.0 ? otherval * SCIPvarGetLbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetUbLocal(consdata->vars[1-i]);
               }

               if( !SCIPisInfinity(scip,consdata->rhs) && consdata->maxactivityneginf + consdata->maxactivityneginf == 0 )
               {
                  rval = val*ub - consdata->rhs;
                  rval += otherval > 0.0 ? otherval * SCIPvarGetUbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetLbLocal(consdata->vars[1-i]);
               }
            }

            newval = MAX(lval, rval);
            assert(SCIPisSumRelLE(scip, newval, val));

            /* Try to avoid cancellation in computation of lhs/rhs */
            newlhs = consdata->lhs - val * lb;
            newlhs += newval * lb;
            newrhs = consdata->rhs - val * ub;
            newrhs += newval * ub;

            if( !SCIPisSumRelEQ(scip, newval, val) )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: change coefficient %+.15g<%s> to %+.15g<%s>, act=[%.15g,%.15g], side=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var),
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               if( SCIPisZero(scip, newval) )
               {
                  SCIP_CALL( delCoefPos(scip, cons, i) );
                  i--;
               }
               else
               {
                  SCIP_CALL( chgCoefPos(scip, cons, i, newval) );
               }
               (*nchgcoefs)++;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
               {
                  SCIPdebugMsg(scip, "linear constraint <%s>: change lhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

                  SCIP_CALL( chgLhs(scip, cons, newlhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->lhs, newlhs));
               }

               if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
               {
                  SCIPdebugMsg(scip, "linear constraint <%s>: change rhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

                  SCIP_CALL( chgRhs(scip, cons, newrhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->rhs, newrhs));
               }
            }
         }
         else
         {
            if( !SCIPisInfinity(scip, -minleftactivity) )
            {
               assert(!SCIPisInfinity(scip, val));
               assert(!SCIPisInfinity(scip, lb));
               if( SCIPisInfinity(scip, -lb) )
                  minleftactivity = -SCIPinfinity(scip);
               else
                  minleftactivity += val * lb;
            }

            if( !SCIPisInfinity(scip, maxleftactivity) )
            {
               assert(!SCIPisInfinity(scip, val));
               assert(!SCIPisInfinity(scip, -ub));
               if( SCIPisInfinity(scip,ub) )
                  maxleftactivity = SCIPinfinity(scip);
               else
                  maxleftactivity += val * ub;
            }
         }
      }
      else
      {
         /* check, if a deviation from lower/upper bound would make lhs/rhs redundant */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS &&
            SCIPisGE(scip, minactivity - val, consdata->lhs) && SCIPisLE(scip, maxactivity + val, consdata->rhs) )
         {
            /* change coefficients:
             *   ai'  := min(rhs - maxact, minact - lhs)
             *   lhs' := lhs - (ai - ai')*ui
             *   rhs' := rhs - (ai - ai')*li
             */

            lval = minactivity - consdata->lhs;
            rval = consdata->rhs - maxactivity;

            /* Try to avoid cancellation, if there are only two variables */
            if( consdata->nvars == 2 )
            {
               SCIP_Real otherval;
               otherval = consdata->vals[1-i];

               if( !SCIPisInfinity(scip,-consdata->lhs) && consdata->minactivityneginf + consdata->minactivityneginf == 0 )
               {
                  lval = val*ub - consdata->lhs;
                  lval += otherval > 0.0 ? otherval * SCIPvarGetLbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetUbLocal(consdata->vars[1-i]);
               }

               if( !SCIPisInfinity(scip,consdata->rhs) && consdata->maxactivityneginf + consdata->maxactivityneginf == 0 )
               {
                  rval = consdata->rhs - val*lb;
                  rval -= otherval > 0.0 ? otherval * SCIPvarGetUbLocal(consdata->vars[1-i]) : otherval * SCIPvarGetLbLocal(consdata->vars[1-i]);
               }
            }

            newval = MIN(lval, rval);
            assert(SCIPisSumRelGE(scip, newval, val));

            /* Try to avoid cancellation in computation of lhs/rhs */
            newlhs = consdata->lhs - val * ub;
            newlhs += newval * ub;
            newrhs = consdata->rhs - val * lb;
            newrhs += newval * lb;

            if( !SCIPisSumRelEQ(scip, newval, val) )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: change coefficient %+.15g<%s> to %+.15g<%s>, act=[%.15g,%.15g], side=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), val, SCIPvarGetName(var), newval, SCIPvarGetName(var),
                  minactivity, maxactivity, consdata->lhs, consdata->rhs);

               /* update the coefficient and the activity bounds */
               if( SCIPisZero(scip, newval) )
               {
                  SCIP_CALL( delCoefPos(scip, cons, i) );
                  i--;
               }
               else
               {
                  SCIP_CALL( chgCoefPos(scip, cons, i, newval) );
               }
               (*nchgcoefs)++;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisEQ(scip, newlhs, consdata->lhs) )
               {
                  SCIPdebugMsg(scip, "linear constraint <%s>: change lhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);

                  SCIP_CALL( chgLhs(scip, cons, newlhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->lhs, newlhs));
               }

               if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisEQ(scip, newrhs, consdata->rhs) )
               {
                  SCIPdebugMsg(scip, "linear constraint <%s>: change rhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);

                  SCIP_CALL( chgRhs(scip, cons, newrhs) );
                  (*nchgsides)++;
                  assert(SCIPisEQ(scip, consdata->rhs, newrhs));
               }
            }
         }
         else
         {
            if( !SCIPisInfinity(scip, -minleftactivity) )
            {
               assert(!SCIPisInfinity(scip, -val));
               assert(!SCIPisInfinity(scip, -ub));
               if( SCIPisInfinity(scip, ub) )
                  minleftactivity = -SCIPinfinity(scip);
               else
                  minleftactivity += val * ub;
            }

            if( !SCIPisInfinity(scip, maxleftactivity) )
            {
               assert(!SCIPisInfinity(scip, -val));
               assert(!SCIPisInfinity(scip, lb));
               if( SCIPisInfinity(scip, -lb) )
                  maxleftactivity = SCIPinfinity(scip);
               else
                  maxleftactivity += val * lb;
            }
         }
      }
      ++i;
   }

   SCIPdebugMsg(scip, "minleftactivity = %.15g, rhs = %.15g\n",
      minleftactivity, consdata->rhs);
   SCIPdebugMsg(scip, "maxleftactivity = %.15g, lhs = %.15g\n",
      maxleftactivity, consdata->lhs);

   /* minleft == \infty  ==>  minactivity == \infty */
   assert(!SCIPisInfinity(scip, -minleftactivity) || SCIPisInfinity(scip, -minactivity));
   assert(!SCIPisInfinity(scip, maxleftactivity) || SCIPisInfinity(scip, maxactivity));

   /* if the lhs is finite, we will check in the following whether the not non-redundant variables can make lhs feasible;
    * this is not valid, if the minactivity is -\infty (aggrlhs would be minus infinity in the following computation)
    * or if huge values contributed to the minactivity, because the minactivity is then just a relaxation
    * (<= the exact minactivity), and we might falsely claim variables to be redundant in the following
    */
   assert(!SCIPisInfinity(scip, minactivity));
   if( !SCIPisInfinity(scip, -consdata->lhs) && (SCIPisInfinity(scip, -minactivity) || minactisrelax) )
      return SCIP_OKAY;

   /* if the rhs is finite, we will check in the following whether the not non-redundant variables can make rhs feasible;
    * this is not valid, if the maxactivity is \infty (aggrrhs would be infinity in the following computation)
    * or if huge values contributed to the maxactivity, because the maxactivity is then just a relaxation
    * (>= the exact maxactivity), and we might falsely claim variables to be redundant in the following
    */
   assert(!SCIPisInfinity(scip, -maxactivity));
   if( !SCIPisInfinity(scip, consdata->rhs) && (SCIPisInfinity(scip, maxactivity) || maxactisrelax) )
      return SCIP_OKAY;

   /* correct lhs and rhs by min/max activity of surely non-redundant variables 
    * surely non-redundant variables are all those where a deviation from the bound makes the lhs/rhs redundant
    */
   aggrlhs = consdata->lhs - minactivity + minleftactivity;
   aggrrhs = consdata->rhs - maxactivity + maxleftactivity;

   /* check if the constraint contains variables which are redundant. The reasoning is the following:
    * Each non-redundant variable can make the lhs/rhs feasible with a deviation of only one in the bound.
    * If _all_ variables which are not non-redundant together cannot make lhs/rhs feasible, 
    * they can be removed from the constraint.
    * aggrrhs may contain some near-infinity value, but only if rhs is infinity.
    */
   if( (SCIPisInfinity(scip, -consdata->lhs) || SCIPisFeasLT(scip, maxleftactivity, aggrlhs))
      && (SCIPisInfinity(scip, consdata->rhs) || SCIPisFeasGT(scip, minleftactivity, aggrrhs)) )
   {
      SCIP_Real minleftactivitypart;
      SCIP_Real maxleftactivitypart;

      assert(!SCIPisInfinity(scip, -consdata->lhs) || !SCIPisInfinity(scip, consdata->rhs));

      /* try to remove redundant variables from constraint */
      i = 0;
      while( i < consdata->nvars )
      {
         var = consdata->vars[i];
         minleftactivitypart = 0.0;
         maxleftactivitypart = 0.0;
         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);

         /* get coefficient and variable's bounds */
         val = consdata->vals[i];
         assert(!SCIPisZero(scip, val));

         /* check sign of coefficient */
         if( val >= 0.0 )
         {     
            /* negation of condition above in case of positive val */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || 
               SCIPisLT(scip, minactivity + val, consdata->lhs) || SCIPisGT(scip, maxactivity - val, consdata->rhs) )
            {
               SCIPdebugMsg(scip, "minactivity = %g\tval = %g\tlhs = %g\n", minactivity, val, consdata->lhs);
               SCIPdebugMsg(scip, "maxactivity = %g\tval = %g\trhs = %g\n", maxactivity, val, consdata->rhs);
               SCIPdebugMsg(scip, "linear constraint <%s>: remove variable <%s> with coefficient <%g> from constraint since it is redundant\n",
                  SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), val);

               minleftactivitypart = val * lb;
               maxleftactivitypart = val * ub;

               SCIP_CALL( delCoefPos(scip, cons, i) );
               i--;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, FALSE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               /* we return above if the condition does not hold and deleting a variable cannot increase the number of
                * huge contributions
                */
               assert(!minactisrelax || SCIPisInfinity(scip, -consdata->lhs));
               assert(!maxactisrelax || SCIPisInfinity(scip, consdata->rhs));
            }
         }
         else 
         {
            /* negation of condition above in case of negative val */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || 
               SCIPisLT(scip, minactivity - val, consdata->lhs) || SCIPisGT(scip, maxactivity + val, consdata->rhs) )
            {
               SCIPdebugMsg(scip, "linear constraint <%s>: remove variable <%s> with coefficient <%g> from constraint since it is redundant\n",
                  SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[i]), val);

               minleftactivitypart = val * ub;
               maxleftactivitypart = val * lb;

               SCIP_CALL( delCoefPos(scip, cons, i) );
               i--;

               /* get the new minimal and maximal activity of the constraint */
               consdataGetActivityBounds(scip, consdata, FALSE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);

               /* we return above if the condition does not hold and deleting a variable cannot increase the number of
                * huge contributions
                */
               assert(!minactisrelax || SCIPisInfinity(scip, -consdata->lhs));
               assert(!maxactisrelax || SCIPisInfinity(scip, consdata->rhs));
            }
         }

         /* the following update step is needed in every iteration cause otherwise it is possible that the surely none-
          * redundant variables could get deleted, 
          * e.g. y_1 + 16y_2 >= 25, y1 with bounds [9,12], y2 with bounds [0,2], minactivity would be 9, it follows that
          * y_2 is surely not redundant and y_1 is redundant so we would first delete y1 and without updating the sides
          * we would also delete y2 and as a result we would have gotten infeasibility */
         /* adjust lhs and right hand side */
         newlhs = consdata->lhs - minleftactivitypart;
         newrhs = consdata->rhs - maxleftactivitypart;

         if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisFeasEQ(scip, newlhs, consdata->lhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s>: change lhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->lhs, newlhs);
            SCIP_CALL( chgLhs(scip, cons, newlhs) );
            ++(*nchgsides);
            assert(SCIPisEQ(scip, consdata->lhs, newlhs));
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisFeasEQ(scip, newrhs, consdata->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s>: change rhs %.15g to %.15g\n", SCIPconsGetName(cons), consdata->rhs, newrhs);
            SCIP_CALL( chgRhs(scip, cons, newrhs) );
            ++(*nchgsides);
            assert(SCIPisEQ(scip, consdata->rhs, newrhs));
         }
         ++i;
      }
   }

   return SCIP_OKAY;
}

/* processes equality with only one variable by fixing the variable and deleting the constraint */
static
SCIP_RETCODE convertUnaryEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real fixval;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 1);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   /* calculate the value to fix the variable to */
   var = consdata->vars[0];
   val = consdata->vals[0];
   assert(!SCIPisZero(scip, val));
   fixval = SCIPselectSimpleValue(consdata->lhs/val - 0.9 * SCIPepsilon(scip),
      consdata->rhs/val + 0.9 * SCIPepsilon(scip), MAXDNOM);
   SCIPdebugMsg(scip, "linear equality <%s>: fix <%s> == %.15g\n",
      SCIPconsGetName(cons), SCIPvarGetName(var), fixval);

   /* fix variable */
   SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeasible, &fixed) );
   if( infeasible )
   {
      SCIPdebugMsg(scip, " -> infeasible fixing\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }
   if( fixed )
      (*nfixedvars)++;

   /* disable constraint */
   SCIP_CALL( SCIPdelCons(scip, cons) );
   if( !consdata->upgraded )
      (*ndelconss)++;

   return SCIP_OKAY;
}

/* processes equality with exactly two variables by aggregating one of the variables and deleting the constraint */
static
SCIP_RETCODE convertBinaryEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   SCIPdebugMsg(scip, "linear constraint <%s>: aggregate %.15g<%s> + %.15g<%s> == %.15g\n",
      SCIPconsGetName(cons), consdata->vals[0], SCIPvarGetName(consdata->vars[0]),
      consdata->vals[1], SCIPvarGetName(consdata->vars[1]), consdata->rhs);

   /* aggregate the equality */
   SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], consdata->vals[0], consdata->vals[1],
         consdata->rhs, &infeasible, &redundant, &aggregated) );

   /* check for infeasibility of aggregation */
   if( infeasible )
   {
      SCIPdebugMsg(scip, " -> infeasible aggregation\n");
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* count the aggregation */
   if( aggregated )
      (*naggrvars)++;

   /* delete the constraint, if it is redundant */
   if( redundant )
   {
      SCIP_CALL( SCIPdelCons(scip, cons) );

      if( !consdata->upgraded )
         (*ndelconss)++;
   }

   return SCIP_OKAY;
}

/** calculates the new lhs and rhs of the constraint after the given variable is aggregated out */
static
void getNewSidesAfterAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_VAR*             slackvar,           /**< variable to be aggregated out */
   SCIP_Real             slackcoef,          /**< coefficient of variable in constraint */
   SCIP_Real*            newlhs,             /**< pointer to store new lhs of constraint */
   SCIP_Real*            newrhs              /**< pointer to store new rhs of constraint */
   )
{
   SCIP_Real slackvarlb;
   SCIP_Real slackvarub;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(newlhs != NULL);
   assert(newrhs != NULL);
   assert(!SCIPisInfinity(scip, -consdata->lhs));
   assert(!SCIPisInfinity(scip, consdata->rhs));

   slackvarlb = SCIPvarGetLbGlobal(slackvar);
   slackvarub = SCIPvarGetUbGlobal(slackvar);
   if( slackcoef > 0.0 )
   {
      if( SCIPisInfinity(scip, -slackvarlb) )
         *newrhs = SCIPinfinity(scip);
      else
         *newrhs = consdata->rhs - slackcoef * slackvarlb;
      if( SCIPisInfinity(scip, slackvarub) )
         *newlhs = -SCIPinfinity(scip);
      else
         *newlhs = consdata->lhs - slackcoef * slackvarub;
   }
   else
   {
      if( SCIPisInfinity(scip, -slackvarlb) )
         *newlhs = -SCIPinfinity(scip);
      else
         *newlhs = consdata->rhs - slackcoef * slackvarlb;
      if( SCIPisInfinity(scip, slackvarub) )
         *newrhs = SCIPinfinity(scip);
      else
         *newrhs = consdata->lhs - slackcoef * slackvarub;
   }
   assert(SCIPisLE(scip, *newlhs, *newrhs));
}

#define MAXMULTIAGGRQUOTIENT 1e+03

/* processes equality with more than two variables by multi-aggregating one of the variables and converting the equality
 * into an inequality; if multi-aggregation is not possible, tries to identify one continuous or integer variable that
 * is implicitly integral by this constraint
 *
 * @todo Check whether a more clever way of avoiding aggregation of variables containing implicitly integer variables
 *       can help.
 */
static
SCIP_RETCODE convertLongEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< linear constraint handler data */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_VARTYPE bestslacktype;
   SCIP_VARTYPE slacktype;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real bestslackdomrng;
   SCIP_Real minabsval;
   SCIP_Real maxabsval;
   SCIP_Bool bestremovescons;
   SCIP_Bool coefszeroone;
   SCIP_Bool coefsintegral;
   SCIP_Bool varsintegral;
   SCIP_Bool infeasible;
   SCIP_Bool samevar;
   int supinf;                               /* counter for infinite contributions to the supremum of a possible
                                              * multi-aggregation
                                              */
   int infinf;                               /* counter for infinite contributions to the infimum of a possible
                                              * multi-aggregation
                                              */
   int maxnlocksstay;
   int maxnlocksremove;
   int bestslackpos;
   int bestnlocks;
   int ncontvars;
   int contvarpos;
   int nintvars;
   int nimplvars;
   int intvarpos;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(naggrvars != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars > 2);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   SCIPdebugMsg(scip, "linear constraint <%s>: try to multi-aggregate equality\n", SCIPconsGetName(cons));

   /* We do not want to increase the total number of non-zeros due to the multi-aggregation.
    * Therefore, we have to restrict the number of locks of a variable that is aggregated out.
    *   maxnlocksstay:   maximal sum of lock numbers if the constraint does not become redundant after the aggregation
    *   maxnlocksremove: maximal sum of lock numbers if the constraint can be deleted after the aggregation
    */
   lhs = consdata->lhs;
   rhs = consdata->rhs;
   maxnlocksstay = 0;
   if( consdata->nvars == 3 )
   {
      /* If the constraint becomes redundant, 3 non-zeros are removed, and we get 1 additional non-zero for each
       * constraint the variable appears in. Thus, the variable must appear in at most 3 other constraints.
       */
      maxnlocksremove = 3;
   }
   else if( consdata->nvars == 4 )
   {
      /* If the constraint becomes redundant, 4 non-zeros are removed, and we get 2 additional non-zeros for each
       * constraint the variable appears in. Thus, the variable must appear in at most 2 other constraints.
       */
      maxnlocksremove = 2;
   }
   else
   {
      /* If the constraint is redundant but has more than 4 variables, we can only accept one other constraint. */
      maxnlocksremove = 1;
   }

   /* the locks on this constraint can be ignored */
   if( SCIPconsIsChecked(cons) )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         maxnlocksstay++;
         maxnlocksremove++;
      }
      if( !SCIPisInfinity(scip, rhs) )
      {
         maxnlocksstay++;
         maxnlocksremove++;
      }
   }

   /* look for a slack variable s to convert a*x + s == b into lhs <= a*x <= rhs */
   vars = consdata->vars;
   vals = consdata->vals;
   bestslackpos = -1;
   bestslacktype = SCIP_VARTYPE_BINARY;
   bestnlocks = INT_MAX;
   bestremovescons = FALSE;
   bestslackdomrng = 0.0;
   coefszeroone = TRUE;
   coefsintegral = TRUE;
   varsintegral = TRUE;
   ncontvars = 0;
   contvarpos = -1;
   nintvars = 0;
   nimplvars = 0;
   intvarpos = -1;
   minabsval = SCIPinfinity(scip);
   maxabsval = -1.0;
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real val;
      SCIP_Real absval;
      SCIP_Real varlb;
      SCIP_Real varub;
      SCIP_Bool iscont;
      int nlocks;

      assert(vars != NULL);
      assert(vals != NULL);

      var = vars[v];
      assert(!SCIPconsIsChecked(cons) || SCIPvarGetNLocksDown(var) >= 1); /* because variable is locked in this equality */
      assert(!SCIPconsIsChecked(cons) || SCIPvarGetNLocksUp(var) >= 1);
      varlb = SCIPvarGetLbGlobal(var);
      varub = SCIPvarGetUbGlobal(var);

      val = vals[v];
      absval = REALABS(val);
      assert(SCIPisPositive(scip, absval));

      /* calculate minimal and maximal absolute value */
      if( absval < minabsval )
         minabsval = absval;
      if( absval > maxabsval )
         maxabsval = absval;

      /* do not try to multi aggregate, when numerical bad */
      if( maxabsval / minabsval > MAXMULTIAGGRQUOTIENT )
         return SCIP_OKAY;

      slacktype = SCIPvarGetType(var);
      coefszeroone = coefszeroone && SCIPisEQ(scip, absval, 1.0);
      coefsintegral = coefsintegral && SCIPisIntegral(scip, val);
      varsintegral = varsintegral && (slacktype != SCIP_VARTYPE_CONTINUOUS);
      iscont = (slacktype == SCIP_VARTYPE_CONTINUOUS || slacktype == SCIP_VARTYPE_IMPLINT);

      /* update candidates for continuous -> implint and integer -> implint conversion */
      if( slacktype == SCIP_VARTYPE_CONTINUOUS )
      {
         ncontvars++;
         contvarpos = v;
      }
      else if( slacktype == SCIP_VARTYPE_IMPLINT )
      {
         ++nimplvars;
      }
      else if( slacktype == SCIP_VARTYPE_INTEGER )
      {
         nintvars++;
         intvarpos = v;
      }

      /* check, if variable is already fixed or aggregated */
      if( !SCIPvarIsActive(var) )
         continue;

      /* check, if variable is used in too many other constraints, even if this constraint could be deleted */
      nlocks = SCIPvarGetNLocksDown(var) + SCIPvarGetNLocksUp(var);

      if( nlocks > maxnlocksremove )
         continue;

      /* check, if variable can be used as a slack variable */
      if( (iscont || (coefsintegral && varsintegral && SCIPisEQ(scip, absval, 1.0))) &&
         !SCIPdoNotMultaggrVar(scip, var) )
      {
         SCIP_Bool better;
         SCIP_Bool equal;
         SCIP_Real slackdomrng;

         if( SCIPisInfinity(scip, varub) || SCIPisInfinity(scip, -varlb) )
            slackdomrng = SCIPinfinity(scip);
         /* we do not want to perform multi-aggregation due to numerics, if the bounds are huge */
         else if( SCIPisHugeValue(scip, varub) || SCIPisHugeValue(scip, -varlb) )
            return SCIP_OKAY;
         else
         {
            slackdomrng = (varub - varlb)*absval;
            assert(!SCIPisInfinity(scip, slackdomrng));
         }
         equal = FALSE;
         better = (slacktype > bestslacktype) || (bestslackpos == -1);
         if( !better && slacktype == bestslacktype )
         {
            better = (nlocks < bestnlocks);
            if( nlocks == bestnlocks && !bestremovescons )
            {
               better = SCIPisGT(scip, slackdomrng, bestslackdomrng);
               equal = !better && SCIPisGE(scip, slackdomrng, bestslackdomrng);
            }
         }

         if( better || equal )
         {
            SCIP_Real minresactivity;
            SCIP_Real maxresactivity;
            SCIP_Real newlhs;
            SCIP_Real newrhs;
            SCIP_Bool removescons;
            SCIP_Bool minisrelax;
            SCIP_Bool maxisrelax;
            SCIP_Bool isminsettoinfinity;
            SCIP_Bool ismaxsettoinfinity;

            /* check if the constraint becomes redundant after multi-aggregation */
            consdataGetActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
               &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);

            /* do not perform the multi-aggregation due to numerics, if we have huge contributions in the residual
             * activity
             */
            if( minisrelax || maxisrelax )
               continue;

            getNewSidesAfterAggregation(scip, consdata, var, val, &newlhs, &newrhs);
            removescons = (SCIPisFeasLE(scip, newlhs, minresactivity) && SCIPisFeasLE(scip, maxresactivity, newrhs));

            /* check resactivities for reliability */
            if( removescons )
            {
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity) )
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);

               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity)
                  && SCIPisFeasLE(scip, newlhs, minresactivity))
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);

               removescons = (SCIPisFeasLE(scip, newlhs, minresactivity) && SCIPisFeasLE(scip, maxresactivity, newrhs));
            }

            /* if parameter multaggrremove is set to TRUE, only aggregate when this removes constraint */
            if( conshdlrdata->multaggrremove && !removescons )
               continue;

            /* prefer variables that make the constraints redundant */
            if( bestremovescons && !removescons )
               continue;

            /* if the constraint does not become redundant, only accept the variable if it does not appear in
             * other constraints
             */
            if( !removescons && nlocks > maxnlocksstay )
               continue;

            better = better || (!bestremovescons && removescons);
            if( better )
            {
               bestslackpos = v;
               bestslacktype = slacktype;
               bestnlocks = nlocks;
               bestslackdomrng = slackdomrng;
               bestremovescons = removescons;
            }
         }
      }
   }

   /* if all coefficients and variables are integral, the right hand side must also be integral */
   if( coefsintegral && varsintegral && !SCIPisFeasIntegral(scip, consdata->rhs) )
   {
      SCIPdebugMsg(scip, "linear equality <%s> is integer infeasible\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   supinf = 0;
   infinf = 0;
   samevar = FALSE;

   /* check whether the the infimum and the supremum of the multi-aggregation can be get infinite */
   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v != bestslackpos )
      {
         if( SCIPisPositive(scip, consdata->vals[v]) )
         {
            if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[v])) )
            {
               ++supinf;
               if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[v])) )
               {
                  ++infinf;
                  samevar = TRUE;
               }
            }
            else if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[v])) )
               ++infinf;

         }
         else if( SCIPisNegative(scip, consdata->vals[v]) )
         {
            if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[v])) )
            {
               ++supinf;
               if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[v])) )
               {
                  ++infinf;
                  samevar = TRUE;
               }
            }
            else if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[v])) )
               ++infinf;
         }
      }
   }
   assert(!samevar || (supinf > 0 && infinf > 0));

   /* If the infimum and the supremum of a multi-aggregation are both infinite, then the multi-aggregation might not be resolvable.
    * E.g., consider the equality z = x-y. If x and y are both fixed to +infinity, the value for z is not determined */
   if( (samevar && (supinf > 1 || infinf > 1)) || (!samevar && supinf > 0 && infinf > 0) )
   {
      SCIPdebugMsg(scip, "do not perform multi-aggregation: infimum and supremum are both infinite\n");
      return SCIP_OKAY;
   }

   /* if the slack variable is of integer type, and the constraint itself may take fractional values,
    * we cannot aggregate the variable, because the integrality condition would get lost
    * Similarly, if there are implicitly integral variables we cannot aggregate, since we might
    * loose the integrality condition for this variable.
    */
   if( bestslackpos >= 0
      && (bestslacktype == SCIP_VARTYPE_CONTINUOUS || bestslacktype == SCIP_VARTYPE_IMPLINT
         || (coefsintegral && varsintegral && nimplvars == 0)) )
   {
      SCIP_VAR* slackvar;
      SCIP_Real* scalars;
      SCIP_Real slackcoef;
      SCIP_Real aggrconst;
      SCIP_Real newlhs;
      SCIP_Real newrhs;
      SCIP_Bool aggregated;

      /* we found a slack variable that only occurs in at most one other constraint:
       *   a_1*x_1 + ... + a_k*x_k + a'*s == rhs  ->  s == rhs - a_1/a'*x_1 - ... - a_k/a'*x_k
       */
      assert(bestslackpos < consdata->nvars);

      /* do not multi aggregate binary variables */
      if( SCIPvarIsBinary(vars[bestslackpos]) )
         return SCIP_OKAY;

      /* convert equality into inequality by deleting the slack variable:
       *  x + a*s == b, l <= s <= u   ->  b - a*u <= x <= b - a*l
       */
      slackvar = vars[bestslackpos];
      slackcoef = vals[bestslackpos];
      assert(!SCIPisZero(scip, slackcoef));
      aggrconst = consdata->rhs/slackcoef;

      getNewSidesAfterAggregation(scip, consdata, slackvar, slackcoef, &newlhs, &newrhs);
      assert(SCIPisLE(scip, newlhs, newrhs));
      SCIP_CALL( chgLhs(scip, cons, newlhs) );
      SCIP_CALL( chgRhs(scip, cons, newrhs) );
      SCIP_CALL( delCoefPos(scip, cons, bestslackpos) );

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &scalars, consdata->nvars) );

      /* set up the multi-aggregation */
      SCIPdebugMsg(scip, "linear constraint <%s>: multi-aggregate <%s> ==", SCIPconsGetName(cons), SCIPvarGetName(slackvar));
      for( v = 0; v < consdata->nvars; ++v )
      {
         scalars[v] = -consdata->vals[v]/slackcoef;
         SCIPdebugMsgPrint(scip, " %+.15g<%s>", scalars[v], SCIPvarGetName(vars[v]));
      }
      SCIPdebugMsgPrint(scip, " %+.15g, bounds of <%s>: [%.15g,%.15g], nlocks=%d, maxnlocks=%d, removescons=%u\n",
         aggrconst, SCIPvarGetName(slackvar), SCIPvarGetLbGlobal(slackvar), SCIPvarGetUbGlobal(slackvar),
         bestnlocks, bestremovescons ? maxnlocksremove : maxnlocksstay, bestremovescons);

      /* perform the multi-aggregation */
      SCIP_CALL( SCIPmultiaggregateVar(scip, slackvar, consdata->nvars, vars, scalars, aggrconst,
            &infeasible, &aggregated) );
      assert(aggregated);

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &scalars);

      /* check for infeasible aggregation */
      if( infeasible )
      {
         SCIPdebugMsg(scip, "linear constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(cons));
         *cutoff = TRUE;
         return SCIP_OKAY;
      }

      (*naggrvars)++;

      /* delete the constraint if it became redundant */
      if( bestremovescons )
      {
         SCIPdebugMsg(scip, "linear constraint <%s>: redundant after multi-aggregation\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );

         if( !consdata->upgraded )
            (*ndelconss)++;
      }
   }
   else if( ncontvars == 1 )
   {
      SCIP_VAR* var;

      assert(0 <= contvarpos && contvarpos < consdata->nvars);
      var = vars[contvarpos];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

      if( coefsintegral && SCIPisFeasIntegral(scip, consdata->rhs) )
      {
         /* upgrade continuous variable to an implicit one, if the absolute value of the coefficient is one */
         if( SCIPisEQ(scip, REALABS(vals[contvarpos]), 1.0) )
         {
            /* convert the continuous variable with coefficient 1.0 into an implicit integer variable */
            SCIPdebugMsg(scip, "linear constraint <%s>: converting continuous variable <%s> to implicit integer variable\n",
               SCIPconsGetName(cons), SCIPvarGetName(var));
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT, &infeasible) );
            if( infeasible )
            {
               SCIPdebugMsg(scip, "infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(var));
               *cutoff = TRUE;

               return SCIP_OKAY;
            }
         }
         /* aggregate continuous variable to an implicit one, if the absolute value of the coefficient is unequal to one */
         /* @todo check if the aggregation coefficient should be in some range(, which is not too big) */
         else if( !SCIPdoNotAggr(scip) )
         {
            SCIP_VAR* newvar;
            SCIP_Real absval;
            char newvarname[SCIP_MAXSTRLEN];
            SCIP_Bool redundant;
            SCIP_Bool aggregated;

            absval = REALABS(vals[contvarpos]);

            (void) SCIPsnprintf(newvarname, SCIP_MAXSTRLEN, "%s_impl", SCIPvarGetName(var));

            /* create new implicit variable for aggregation */
            SCIP_CALL( SCIPcreateVar(scip, &newvar, newvarname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
                  SCIP_VARTYPE_IMPLINT, SCIPvarIsInitial(var), SCIPvarIsRemovable(var), NULL, NULL, NULL, NULL, NULL) );

            /* add new variable to problem */
            SCIP_CALL( SCIPaddVar(scip, newvar) );

#ifdef WITH_DEBUG_SOLUTION
            if( SCIPdebugIsMainscip(scip) )
            {
               SCIP_Real varval;
               SCIP_CALL( SCIPdebugGetSolVal(scip, var, &varval) );
               SCIP_CALL( SCIPdebugAddSolVal(scip, newvar, absval * varval) );
            }
#endif

            /* convert the continuous variable with coefficient 1.0 into an implicit integer variable */
            SCIPdebugMsg(scip, "linear constraint <%s>: aggregating continuous variable <%s> to newly created implicit integer variable <%s>, aggregation factor = %g\n",
               SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetName(newvar), absval);

            /* aggregate continuous and implicit variable */
            SCIP_CALL( SCIPaggregateVars(scip, var, newvar, absval, -1.0, 0.0, &infeasible, &redundant, &aggregated) );

            if( infeasible )
            {
               SCIPdebugMsg(scip, "infeasible aggregation of variable <%s> to implicit variable <%s>, domain is empty\n",
                  SCIPvarGetName(var), SCIPvarGetName(newvar));
               *cutoff = TRUE;

               /* release implicit variable */
               SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

               return SCIP_OKAY;
            }

            /* release implicit variable */
            SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

            if( aggregated )
               (*naggrvars)++;
            else
               return SCIP_OKAY;
         }

         /* we do not have any event on vartype changes, so we need to manually force this constraint to be presolved
          * again
          */
         consdata->boundstightened = 0;
         consdata->rangedrowpropagated = 0;
         consdata->presolved = FALSE;
      }
   }
   else if( ncontvars == 0 && nimplvars == 0 && nintvars == 1 && !coefszeroone )
   {
      SCIP_VAR* var;

      /* this seems to help for rococo instances, but does not for rout (where all coefficients are +/- 1.0)
       *  -> we don't convert integers into implints if the row is a 0/1-row
       */
      assert(varsintegral);
      assert(0 <= intvarpos && intvarpos < consdata->nvars);
      var = vars[intvarpos];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      if( coefsintegral
         && SCIPisEQ(scip, REALABS(vals[intvarpos]), 1.0)
         && SCIPisFeasIntegral(scip, consdata->rhs) )
      {
         /* convert the integer variable with coefficient 1.0 into an implicit integer variable */
         SCIPdebugMsg(scip, "linear constraint <%s>: converting integer variable <%s> to implicit integer variable\n",
            SCIPconsGetName(cons), SCIPvarGetName(var));
         SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT, &infeasible) );
         if( infeasible )
         {
            SCIPdebugMsg(scip, "infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(var));
            *cutoff = TRUE;

            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** checks if the given variables and their coefficient are equal (w.r.t. scaling factor) to the objective function */
static
SCIP_Bool checkEqualObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   SCIP_Real*            scale,              /**< pointer to store the scaling factor between the constraint and the
					      *   objective function */
   SCIP_Real*            offset              /**< pointer to store the offset of the objective function resulting by
					      *   this constraint */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real objval;
   SCIP_Bool negated;
   int nvars;
   int v;

   vars = consdata->vars;
   nvars = consdata->nvars;

   assert(vars != NULL);

   for( v = 0; v < nvars; ++v )
   {
      negated = FALSE;
      var = vars[v];
      assert(var != NULL);

      if( SCIPvarIsNegated(var) )
      {
         negated = TRUE;
         var = SCIPvarGetNegatedVar(var);
         assert(var != NULL);
      }

      objval = SCIPvarGetObj(var);

      /* if a variable has a zero objective coefficient the linear constraint is not a subset of the objective
       * function
       */
      if( SCIPisZero(scip, objval) )
         return FALSE;
      else
      {
         SCIP_Real val;

         val = consdata->vals[v];

         if( negated )
         {
            if( v == 0 )
            {
               /* the first variable defines the scale */
               (*scale) = val / -objval;

               (*offset) += val;
            }
            else if( SCIPisEQ(scip, -objval * (*scale), val) )
               (*offset) += val;
            else
               return FALSE;
         }
         else if( v == 0 )
         {
            /* the first variable defines the scale */
            (*scale) = val / objval;
         }
         else if( !SCIPisEQ(scip, objval * (*scale), val) )
            return FALSE;
      }
   }

   return TRUE;
}

/** check if the linear equality constraint is equal to a subset of the objective function; if so we can remove the
 *  objective coefficients and add an objective offset
 */
static
SCIP_RETCODE checkPartialObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear equation constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< linear constraint handler data */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real offset;
   SCIP_Real scale;
   SCIP_Bool applicable;
   int nobjvars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(SCIPisEQ(scip, consdata->lhs, consdata->rhs));

   nvars = consdata->nvars;
   nobjvars = SCIPgetNObjVars(scip);

   /* check if the linear equality constraints does not have more variables than the objective function */
   if( nvars > nobjvars || nvars == 0 )
      return SCIP_OKAY;

   /* check for allowance of algorithm */
   if( (nvars < nobjvars && !conshdlrdata->detectpartialobjective) ||
      (nvars == nobjvars && (!conshdlrdata->detectcutoffbound || !conshdlrdata->detectlowerbound)) )
      return SCIP_OKAY;

   offset = consdata->rhs;
   scale = 1.0;

   /* checks if the variables and their coefficients are equal (w.r.t. scaling factor) to the objective function */
   applicable = checkEqualObjective(scip, consdata, &scale, &offset);

   if( applicable )
   {
      SCIP_VAR** vars;

      vars = consdata->vars;
      assert(vars != NULL);

      offset /= scale;

      SCIPdebugMsg(scip, "linear equality constraint <%s> == %g (offset %g) is a subset of the objective function\n",
         SCIPconsGetName(cons), consdata->rhs, offset);

      /* set all objective coefficient to zero */
      for( v = 0; v < nvars; ++v )
      {
         SCIP_CALL( SCIPchgVarObj(scip, vars[v], 0.0) );
      }

      /* add an objective offset */
      SCIP_CALL( SCIPaddObjoffset(scip, offset) );
   }

   return SCIP_OKAY;
}

/** updates the cutoff if the given primal bound (which is implied by the given constraint) is better */
static
SCIP_RETCODE updateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             primalbound         /**< feasible primal bound */
   )
{
   SCIP_Real cutoffbound;

   /* increase the cutoff bound value by an epsilon to ensue that solution with the value of the cutoff bound are still
    * accepted
    */
   cutoffbound = primalbound + SCIPcutoffbounddelta(scip);

   if( cutoffbound < SCIPgetCutoffbound(scip) )
   {
      SCIPdebugMsg(scip, "update cutoff bound <%g>\n", cutoffbound);

      SCIP_CALL( SCIPupdateCutoffbound(scip, cutoffbound) );
   }
   else
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* we cannot disable the enforcement and propagation on ranged rows, because the cutoffbound could only have
       * resulted from one side
       */
      if( SCIPisInfinity(scip, -consdata->lhs) || SCIPisInfinity(scip, consdata->rhs) )
      {
         /* in case the cutoff bound is worse then the currently known one, we additionally avoid enforcement and
          * propagation
          */
         SCIP_CALL( SCIPsetConsEnforced(scip, cons, FALSE) );
         SCIP_CALL( SCIPsetConsPropagated(scip, cons, FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** check if the linear constraint is parallel to objective function; if so update the cutoff bound and avoid that the
 *  constraint enters the LP by setting the initial and separated flag to FALSE
 */
static
SCIP_RETCODE checkParallelObjective(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< linear constraint handler data */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real offset;
   SCIP_Real scale;
   SCIP_Bool applicable;
   int nobjvars;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* ignore equalities since these are covered by the method checkPartialObjective() */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      return SCIP_OKAY;

   nvars = consdata->nvars;
   nobjvars = SCIPgetNObjVars(scip);

   /* check if the linear inequality constraints has the same number of variables as the objective function and if the
    * initial and/or separated flag is set to FALSE
    */
   if( nvars != nobjvars || (!SCIPconsIsInitial(cons) && !SCIPconsIsSeparated(cons)) )
      return SCIP_OKAY;

   offset = 0.0;
   scale = 1.0;

   /* There are no variables in the ojective function and in the constraint. Thus, the constraint is redundant or proves
    * infeasibility. Since we have a pure feasibility problem, we do not want to set a cutoff or lower bound.
    */
   if( nobjvars == 0 )
      return SCIP_OKAY;

   /* checks if the variables and their coefficients are equal (w.r.t. scaling factor) to the objective function */
   applicable = checkEqualObjective(scip, consdata, &scale, &offset);

   if( applicable )
   {
      SCIP_Bool rhsfinite = !SCIPisInfinity(scip, consdata->rhs);
      SCIP_Bool lhsfinite = !SCIPisInfinity(scip, -consdata->lhs);

      if( SCIPisPositive(scip, scale) )
      {
         if( conshdlrdata->detectcutoffbound && rhsfinite )
         {
            SCIP_Real primalbound;

            primalbound = (consdata->rhs - offset) / scale;

            SCIPdebugMsg(scip, "constraint <%s> is parallel to objective function and provides a cutoff bound <%g>\n",
               SCIPconsGetName(cons), primalbound);

            SCIP_CALL( updateCutoffbound(scip, cons, primalbound) );
         }

         if( conshdlrdata->detectlowerbound && lhsfinite )
         {
            SCIP_Real lowerbound;

            lowerbound = (consdata->lhs - offset) / scale;

            SCIPdebugMsg(scip, "constraint <%s> is parallel to objective function and provides a lower bound <%g>\n",
               SCIPconsGetName(cons), lowerbound);

            SCIP_CALL( SCIPupdateLocalLowerbound(scip, lowerbound) );
         }

         if( (conshdlrdata->detectcutoffbound && (conshdlrdata->detectlowerbound || !lhsfinite)) ||
            (conshdlrdata->detectlowerbound && !rhsfinite) )
         {
            /* avoid that the linear constraint enters the LP since it is parallel to the objective function */
            SCIP_CALL( SCIPsetConsInitial(scip, cons, FALSE) );
            SCIP_CALL( SCIPsetConsSeparated(scip, cons, FALSE) );
         }
      }
      else
      {
         if( conshdlrdata->detectlowerbound && rhsfinite )
         {
            SCIP_Real lowerbound;

            lowerbound = (consdata->rhs - offset) / scale;

            SCIPdebugMsg(scip, "constraint <%s> is parallel to objective function and provides a lower bound <%g>\n",
               SCIPconsGetName(cons), lowerbound);

            SCIP_CALL( SCIPupdateLocalLowerbound(scip, lowerbound) );
         }

         if( conshdlrdata->detectcutoffbound && lhsfinite )
         {
            SCIP_Real primalbound;

            primalbound = (consdata->lhs - offset) / scale;

            SCIPdebugMsg(scip, "constraint <%s> is parallel to objective function and provides a cutoff bound <%g>\n",
               SCIPconsGetName(cons), primalbound);

            SCIP_CALL( updateCutoffbound(scip, cons, primalbound) );
         }

         if( (conshdlrdata->detectcutoffbound && (conshdlrdata->detectlowerbound || !rhsfinite)) ||
            (conshdlrdata->detectlowerbound && !lhsfinite) )
         {
            /* avoid that the linear constraint enters the LP since it is parallel to the objective function */
            SCIP_CALL( SCIPsetConsInitial(scip, cons, FALSE) );
            SCIP_CALL( SCIPsetConsSeparated(scip, cons, FALSE) );
         }
      }
   }

   return SCIP_OKAY;
}

/** converts special equalities */
static
SCIP_RETCODE convertEquality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< linear constraint handler data */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlrdata != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->removedfixings);

   /* do nothing on inequalities */
   if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      return SCIP_OKAY;

   /* depending on the number of variables, call a special conversion method */
   if( consdata->nvars == 1 )
   {
      /* fix variable */
      SCIP_CALL( convertUnaryEquality(scip, cons, cutoff, nfixedvars, ndelconss) );
   }
   else if( consdata->nvars == 2 )
   {
      /* aggregate one of the variables */
      SCIP_CALL( convertBinaryEquality(scip, cons, cutoff, naggrvars, ndelconss) );
   }
   else
   {
      /* check if the equality is part of the objective function */
      SCIP_CALL( checkPartialObjective(scip, cons, conshdlrdata) );

      /* try to multi-aggregate one of the variables */
      SCIP_CALL( convertLongEquality(scip, conshdlrdata, cons, cutoff, naggrvars, ndelconss) );
   }

   return SCIP_OKAY;
}

/** returns whether the linear sum of all variables/coefficients except the given one divided by the given value is always
 *  integral
 */
static
SCIP_Bool consdataIsResidualIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint */
   int                   pos,                /**< position of variable to be left out */
   SCIP_Real             val                 /**< value to divide the coefficients by */
   )
{
   int v;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);

   for( v = 0; v < consdata->nvars; ++v )
   {
      if( v != pos && (!SCIPvarIsIntegral(consdata->vars[v]) || !SCIPisIntegral(scip, consdata->vals[v]/val)) )
         return FALSE;
   }

   return TRUE;
}

/* check if lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i 
 * check if rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i 
 */
static
void calculateMinvalAndMaxval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             side,               /**< lhs or rhs */
   SCIP_Real             val,                /**< coefficient */
   SCIP_Real             minresactivity,     /**< minimal residual activity */
   SCIP_Real             maxresactivity,     /**< maximal residual activity */
   SCIP_Real*            minval,             /**< pointer to store calculated minval */
   SCIP_Real*            maxval              /**< pointer to store calculated maxval */
   )
{
   assert(scip != NULL);
   assert(minval != NULL);
   assert(maxval != NULL);

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, ABS(maxresactivity)) )
         *minval = -maxresactivity;
      else
         *minval = (side - maxresactivity)/val;

      if( SCIPisInfinity(scip, ABS(minresactivity)) )
         *maxval = -minresactivity;
      else
         *maxval = (side - minresactivity)/val;
   }
   else
   {
      if( SCIPisInfinity(scip, ABS(minresactivity)) )
         *minval = minresactivity;
      else
         *minval = (side - minresactivity)/val;

      if( SCIPisInfinity(scip, ABS(maxresactivity)) )
         *maxval = maxresactivity;
      else
         *maxval = (side - maxresactivity)/val;
   }
}


/* applies dual presolving for variables that are locked only once in a direction, and this locking is due to a
 * linear inequality
 */
static
SCIP_RETCODE dualPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool lhsexists;
   SCIP_Bool rhsexists;
   SCIP_Bool bestisint;
   SCIP_Bool bestislhs;
   int bestpos;
   int i;
   int maxotherlocks;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);

   /* only process checked constraints (for which the locks are increased);
    * otherwise we would have to check for variables with nlocks == 0, and these are already processed by the
    * dualfix presolver
    */
   if( !SCIPconsIsChecked(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhsexists = !SCIPisInfinity(scip, -consdata->lhs);
   rhsexists = !SCIPisInfinity(scip, consdata->rhs);

   /* search for a single-locked variable which can be multi-aggregated; if a valid continuous variable was found, we
    * can use it safely for aggregation and break the search loop
    */
   bestpos = -1;
   bestisint = TRUE;
   bestislhs = FALSE;

   /* We only want to multi-aggregate variables, if they appear in maximal one additional constraint,
    * everything else would produce fill-in. Exceptions:
    * - If there are only two variables in the constraint from which the multi-aggregation arises, no fill-in will be
    *   produced.
    * - If there are three variables in the constraint, multi-aggregation in three additional constraints will remove
    *   six nonzeros (three from the constraint and the three entries of the multi-aggregated variable) and add
    *   six nonzeros (two variables per substitution).
    * - If there at most four variables in the constraint, multi-aggregation in two additional constraints will remove
    *   six nonzeros (four from the constraint and the two entries of the multi-aggregated variable) and add
    *   six nonzeros (three variables per substitution). God exists! 
    */
   if( consdata->nvars <= 2 )
      maxotherlocks = INT_MAX;
   else if( consdata->nvars == 3 )
      maxotherlocks = 3;
   else if( consdata->nvars == 4 )
      maxotherlocks = 2;
   else
      maxotherlocks = 1;

   /* if this constraint has both sides, it also provides a lock for the other side and thus we can allow one more lock */
   if( lhsexists && rhsexists && maxotherlocks < INT_MAX )
      maxotherlocks++;

   for( i = 0; i < consdata->nvars && bestisint; ++i )
   {
      SCIP_VAR* var;
      SCIP_Bool isint;
      SCIP_Real val;
      SCIP_Real obj;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool agglhs;
      SCIP_Bool aggrhs;

      var = consdata->vars[i];
      isint = (SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER);

      /* if we already found a candidate, skip integers */
      if( bestpos >= 0 && isint )
         continue;

      /* better do not multi-aggregate binary variables, since most plugins rely on their binary variables to be either
       * active, fixed, or single-aggregated with another binary variable
       */
      if( SCIPvarIsBinary(var) && consdata->nvars > 2 )
         continue;

      if ( SCIPdoNotMultaggrVar(scip, var) )
         continue;

      val = consdata->vals[i];
      obj = SCIPvarGetObj(var);
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      /* lhs <= a_0 * x_0 + a_1 * x_1 + ... + a_{n-1} * x_{n-1} <= rhs
       *
       * a_i >= 0, c_i >= 0, lhs exists, nlocksdown(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its lower bound
       *  - fix x_i to the smallest value for this constraint: x_i := lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i <= 0, c_i <= 0, lhs exists, nlocksup(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its upper bound
       *  - fix x_i to the largest value for this constraint: x_i := lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i >= 0, c_i <= 0, rhs exists, nlocksup(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its upper bound
       *  - fix x_i to the largest value for this constraint: x_i := rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *
       * a_i <= 0, c_i >= 0, rhs exists, nlocksdown(x_i) == 1:
       *  - constraint is the only one that forbids fixing the variable to its lower bound
       *  - fix x_i to the smallest value for this constraint: x_i := rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j
       *       
       * but: all this is only applicable, if the aggregated value is inside x_i's bounds for all possible values
       *      of all x_j
       * furthermore: we only want to apply this, if no fill-in will be produced
       */
      agglhs = lhsexists
         && ((val > 0.0 && !SCIPisNegative(scip, obj) && SCIPvarGetNLocksDown(var) == 1 
               && SCIPvarGetNLocksUp(var) <= maxotherlocks)
            || (val < 0.0 && !SCIPisPositive(scip, obj) && SCIPvarGetNLocksUp(var) == 1
               && SCIPvarGetNLocksDown(var) <= maxotherlocks));
      aggrhs = rhsexists
         && ((val > 0.0 && !SCIPisPositive(scip, obj) && SCIPvarGetNLocksUp(var) == 1 
               && SCIPvarGetNLocksDown(var) <= maxotherlocks)            
            || (val < 0.0 && !SCIPisNegative(scip, obj)  && SCIPvarGetNLocksDown(var) == 1 
               && SCIPvarGetNLocksUp(var) <= maxotherlocks));
      if( agglhs || aggrhs )
      {
         SCIP_Real minresactivity;
         SCIP_Real maxresactivity;
         SCIP_Real minval;
         SCIP_Real maxval;
         SCIP_Bool minisrelax;
         SCIP_Bool maxisrelax;
         SCIP_Bool isminsettoinfinity;
         SCIP_Bool ismaxsettoinfinity;

         /* calculate bounds for \sum_{j \neq i} a_j * x_j */
         consdataGetActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
            &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);
         assert(SCIPisLE(scip, minresactivity, maxresactivity));

         /* We called consdataGetActivityResiduals() saying that we do not need a good relaxation,
          * so whenever we have a relaxed activity, it should be relaxed to +/- infinity.
          * This is needed, because we do not want to rely on relaxed finite resactivities.
          */
         assert((!minisrelax || isminsettoinfinity) && (!maxisrelax || ismaxsettoinfinity));

         if( agglhs )
         {
            /* check if lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
            calculateMinvalAndMaxval(scip, consdata->lhs, val, minresactivity, maxresactivity, &minval, &maxval);

            assert(SCIPisLE(scip, minval, maxval));
            if( (!SCIPisInfinity(scip, -minval) && SCIPisFeasGE(scip, minval, lb)) &&
               (!SCIPisInfinity(scip, maxval) && SCIPisFeasLE(scip, maxval, ub)) )
            {
               SCIP_Real oldmaxresactivity;
               SCIP_Real oldminresactivity;
               SCIP_Bool recalculated;

               recalculated = FALSE;
               oldmaxresactivity = maxresactivity;
               oldminresactivity = minresactivity;

               /* check minresactivity for reliability */
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
                  recalculated = !SCIPisEQ(scip, oldminresactivity, minresactivity);
                  isminsettoinfinity = TRUE; /* here it means only that it was even calculated */
               }

               /* check maxresactivity for reliability */
               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
                  recalculated = recalculated || !SCIPisEQ(scip, oldmaxresactivity, maxresactivity);
                  ismaxsettoinfinity = TRUE; /* here it means only that it was even calculated */
               }

               /* minresactivity or maxresactivity wasn't reliable so recalculate min- and maxval*/
               if( recalculated )
               {
                  assert(SCIPisLE(scip, minresactivity, maxresactivity));

                  /* check again if lhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
                  calculateMinvalAndMaxval(scip, consdata->lhs, val, minresactivity, maxresactivity, &minval, &maxval);

                  assert(SCIPisLE(scip, minval, maxval));
               }

               if( !recalculated || (SCIPisFeasGE(scip, minval, lb) && SCIPisFeasLE(scip, maxval, ub)) )
               {
                  /* if the variable is integer, we have to check whether the integrality condition would always be satisfied
                   * in the multi-aggregation
                   */
                  if( !isint || (SCIPisIntegral(scip, consdata->lhs/val) && consdataIsResidualIntegral(scip, consdata, i, val)) )
                  {
                     bestpos = i;
                     bestisint = isint;
                     bestislhs = TRUE;
                     continue; /* no need to also look at the right hand side */
                  }
               }
            }
         }

         if( aggrhs )
         {
            /* check if rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
            calculateMinvalAndMaxval(scip, consdata->rhs, val, minresactivity, maxresactivity, &minval, &maxval);

            assert(SCIPisLE(scip,minval,maxval));
            if( (!SCIPisInfinity(scip, -minval) && SCIPisFeasGE(scip, minval, lb)) &&
               (!SCIPisInfinity(scip, maxval) && SCIPisFeasLE(scip, maxval, ub)) )
            {
               SCIP_Real oldmaxresactivity;
               SCIP_Real oldminresactivity;
               SCIP_Bool recalculated;

               recalculated = FALSE;
               oldmaxresactivity = maxresactivity;
               oldminresactivity = minresactivity;

               /* check minresactivity for reliability */
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastminactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, FALSE);
                  recalculated = !SCIPisEQ(scip, oldminresactivity, minresactivity);
               }

               /* check maxresactivity for reliability */
               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastmaxactivity) )
               {
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, FALSE);
                  recalculated = recalculated || !SCIPisEQ(scip, oldmaxresactivity, maxresactivity);
               }

               /* minresactivity or maxresactivity wasn't reliable so recalculate min- and maxval*/
               if( recalculated )
               {
                  /* check again if rhs/a_i - \sum_{j \neq i} a_j/a_i * x_j is always inside the bounds of x_i */
                  calculateMinvalAndMaxval(scip, consdata->rhs, val, minresactivity, maxresactivity, &minval, &maxval);
                  assert(SCIPisLE(scip,minval,maxval));
               }

               if( !recalculated || (SCIPisFeasGE(scip, minval, lb) && SCIPisFeasLE(scip, maxval, ub)) )
               {
                  /* if the variable is integer, we have to check whether the integrality condition would always be satisfied
                   * in the multi-aggregation
                   */
                  if( !isint || (SCIPisIntegral(scip, consdata->rhs/val) && consdataIsResidualIntegral(scip, consdata, i, val)) )
                  {
                     bestpos = i;
                     bestisint = isint;
                     bestislhs = FALSE;
                  }
               }
            }
         }
      }
   }

   if( bestpos >= 0 )
   {
      SCIP_VAR** aggrvars;
      SCIP_Real* aggrcoefs;
      SCIP_Real aggrconst;
      SCIP_VAR* bestvar;
      SCIP_Real bestval;
      SCIP_Real epsilon;
      int naggrs;
      int j;
      SCIP_Bool infeasible;
      SCIP_Bool aggregated;
      SCIP_Bool samevar;
      int supinf;                            /* counter for infinite contributions to the supremum of a possible
                                              * multi-aggregation
                                              */
      int infinf;                            /* counter for infinite contributions to the infimum of a possible
                                              * multi-aggregation
                                              */

      assert(!bestislhs || lhsexists);
      assert(bestislhs || rhsexists);

      bestvar = consdata->vars[bestpos];
      bestval = consdata->vals[bestpos];
      assert(bestisint ==
         (SCIPvarGetType(bestvar) == SCIP_VARTYPE_BINARY || SCIPvarGetType(bestvar) == SCIP_VARTYPE_INTEGER));

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &aggrvars, consdata->nvars-1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &aggrcoefs, consdata->nvars-1) );

      /* set up the multi-aggregation */
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMsg(scip, "linear constraint <%s> (dual): multi-aggregate <%s> ==", SCIPconsGetName(cons), SCIPvarGetName(bestvar));
      naggrs = 0;
      supinf = 0;
      infinf = 0;
      samevar = FALSE;
      epsilon = SCIPepsilon(scip);

      for( j = 0; j < consdata->nvars; ++j )
      {
         if( j != bestpos )
         {
            SCIP_Real absaggrcoef;

            aggrvars[naggrs] = consdata->vars[j];
            aggrcoefs[naggrs] = -consdata->vals[j]/consdata->vals[bestpos];
            SCIPdebugMsgPrint(scip, " %+.15g<%s>", aggrcoefs[naggrs], SCIPvarGetName(aggrvars[naggrs]));

            absaggrcoef = REALABS(aggrcoefs[naggrs]);

            /* do not try to multi aggregate, when numerical bad */
            if( absaggrcoef < epsilon )
            {
               SCIPdebugMsg(scip, "do not perform multi-aggregation: too large aggregation coefficients\n");

               /* free temporary memory */
               SCIPfreeBufferArray(scip, &aggrcoefs);
               SCIPfreeBufferArray(scip, &aggrvars);

               return SCIP_OKAY;
            }

            if( bestisint )
            {
               /* coefficient must be integral: round it to exact integral value */
               assert(SCIPisIntegral(scip, aggrcoefs[naggrs]));
               aggrcoefs[naggrs] = SCIPfloor(scip, aggrcoefs[naggrs]+0.5);
            }

            if( SCIPisPositive(scip, aggrcoefs[naggrs]) )
            {
               if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[j])) )
               {
                  ++supinf;
                  if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[j])) )
                  {
                     ++infinf;
                     samevar = TRUE;
                  }
               }
               else if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[j])) )
                  ++infinf;
            }
            else if( SCIPisNegative(scip, aggrcoefs[naggrs]) )
            {
               if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(consdata->vars[j])) )
               {
                  ++supinf;
                  if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[j])) )
                  {
                     ++infinf;
                     samevar = TRUE;
                  }
               }
               else if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(consdata->vars[j])) )
                  ++infinf;
            }

            naggrs++;
         }
      }
      assert(!samevar || (supinf > 0 && infinf > 0));

      aggrconst = (bestislhs ? consdata->lhs/bestval : consdata->rhs/bestval);
      SCIPdebugMsgPrint(scip, " %+.15g, bounds of <%s>: [%.15g,%.15g]\n", aggrconst, SCIPvarGetName(bestvar),
         SCIPvarGetLbGlobal(bestvar), SCIPvarGetUbGlobal(bestvar));
      assert(naggrs == consdata->nvars-1);

      /* right hand side must be integral: round it to exact integral value */
      if( bestisint )
      {
         assert(SCIPisIntegral(scip, aggrconst));
         aggrconst = SCIPfloor(scip, aggrconst+0.5);
      }

      aggregated = FALSE;
      infeasible = FALSE;

      /* perform the multi-aggregation */
      if( (samevar && supinf == 1 && infinf == 1) || (!samevar && (supinf == 0 || infinf == 0)) )
      {
         /* @todo if multi-aggregate makes them numerical trouble, avoid them if the coefficients differ to much, see
          * also convertLongEquality() early termination due to coefficients
          */
         SCIP_CALL( SCIPmultiaggregateVar(scip, bestvar, naggrs, aggrvars, aggrcoefs, aggrconst, &infeasible, &aggregated) );
      }
      else
      {
         /* If the infimum and the supremum of a multi-aggregation are both infinite, then the multi-aggregation might not be resolvable.
          * E.g., consider the equality z = x-y. If x and y are both fixed to +infinity, the value for z is not determined */
         SCIPdebugMsg(scip, "do not perform multi-aggregation: infimum and supremum are both infinite\n");
      }
      /* free temporary memory */
      SCIPfreeBufferArray(scip, &aggrcoefs);
      SCIPfreeBufferArray(scip, &aggrvars);

      /* check for infeasible aggregation */
      if( infeasible )
      {
         SCIPdebugMsg(scip, "linear constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(cons));
         *cutoff = TRUE;
         return SCIP_OKAY;
      }

      /* delete the constraint, if the aggregation was successful */
      if( aggregated )
      {
         SCIP_CALL( SCIPdelCons(scip, cons) );

         if( !consdata->upgraded )
            (*ndelconss)++;
         (*naggrvars)++;
      }
      else
      {
         SCIPdebugMsg(scip, "aggregation non successful!\n");
      }
   }

   return SCIP_OKAY;
}

#define BINWEIGHT  1
#define INTWEIGHT  4
#define CONTWEIGHT 8

/** gets weight for variable in a "weighted number of variables" sum */
static
int getVarWeight(
   SCIP_VAR*             var                 /**< variable to get weight for */
   )
{
   switch( SCIPvarGetType(var) )
   {
   case SCIP_VARTYPE_BINARY:
      return BINWEIGHT;
   case SCIP_VARTYPE_INTEGER:
   case SCIP_VARTYPE_IMPLINT:
      return INTWEIGHT;
   case SCIP_VARTYPE_CONTINUOUS:
      return CONTWEIGHT;
   default:
      SCIPerrorMessage("invalid variable type\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }
}

/** tries to aggregate variables in equations a^Tx = lhs
 *  in case there are at most two binary variables with an odd coefficient and all other
 *  variables are not continuous and have an even coefficient then:
 *  - exactly one odd binary variables 
 *    this binary variables y can be fixed to 0 if the lhs is even and to 1 if the lhs is odd    
 *     - lhs is odd ->  y = 1
 *     - lhs is even -> y = 0
 *  - exactly two odd binary variables 
 *    aggregate the two binary variables with odd coefficient 
 *     - lhs is odd -> exactly one of the variable has to be 1 -> var1 + var2 = 1
 *     - lhs is even -> both have to take the same value -> var1 - var2 = 0
 */
static
SCIP_RETCODE aggregateVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count number of fixed variables */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* check if the linear constraint is an equation with integral right hand side */
   if( !SCIPisEQ(scip, consdata->lhs, consdata->rhs) || !SCIPisIntegral(scip, consdata->lhs) )
      return SCIP_OKAY;

   /* try to fix and aggregated variables until nothing is possible anymore */
   do
   {
      int v;
      int nvars;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real lhs;
      SCIP_Bool lhsodd;

      SCIP_Bool infeasible;
      SCIP_Bool fixed;
      SCIP_Bool aggregated;
      SCIP_Bool redundant;

      SCIP_VAR* var1;
      SCIP_VAR* var2;
      int noddvars;

      success = FALSE;

      lhs = consdata->lhs;
      vars = consdata->vars;
      vals = consdata->vals;
      nvars = consdata->nvars;

      assert( !SCIPisInfinity(scip, ABS(lhs)) );

      var1 = NULL;
      var2 = NULL;
      noddvars = 0;

      /* search for binary variables with an odd coefficient */
      for( v = 0; v < nvars && noddvars < 3; ++v )
      {
         SCIP_Longint val;

         /* all coefficients and variables have to be integral */
         if( !SCIPisIntegral(scip, vals[v]) || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            return SCIP_OKAY;

         val = (SCIP_Longint)SCIPfeasFloor(scip, vals[v]);
         if( val % 2 != 0 )
         {
            /* the odd values have to belong to binary variables */
            if( !SCIPvarIsBinary(vars[v]) )
               return SCIP_OKAY;

            if( noddvars == 0 )
               var1 = vars[v];
            else
               var2 = vars[v];

            noddvars++;
         }
      }

      /* check lhs is odd or even */
      lhsodd = (((SCIP_Longint)SCIPfeasFloor(scip, lhs)) % 2 != 0);

      if( noddvars == 1 )
      {
         assert( var1 != NULL );

         SCIPdebugMsg(scip, "linear constraint <%s>: try fixing variable <%s> to <%g>\n",
            SCIPconsGetName(cons), SCIPvarGetName(var1), lhsodd ? 1.0 : 0.0);

         SCIP_CALL( SCIPfixVar(scip, var1, lhsodd? 1.0 : 0.0, &infeasible, &fixed) );

         /* check for infeasibility of fixing */
         if( infeasible )
         {
            SCIPdebugMsg(scip, " -> infeasible fixing\n");
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         if( fixed )
         {
            SCIPdebugMsg(scip, " -> feasible fixing\n");
            (*nfixedvars)++;
            success = TRUE;
         }
      }
      else if( noddvars == 2 )
      {
         assert( var1 != NULL );
         assert( var2 != NULL );

         /* aggregate the two variables with odd coefficient 
          * - lhs is odd -> exactly one of the variable has to be 1 -> var1 + var2 = 1
          * - lhs is even -> both have to take the same value -> var1 - var2 = 0
          */
         SCIPdebugMsg(scip, "linear constraint <%s>: try aggregation of variables <%s> and <%s>\n",
            SCIPconsGetName(cons), SCIPvarGetName(var1), SCIPvarGetName(var2));

         SCIP_CALL( SCIPaggregateVars(scip, var1, var2, 1.0, lhsodd ? 1.0 : -1.0,
               lhsodd ? 1.0 : 0.0, &infeasible, &redundant, &aggregated) );

         /* check for infeasibility of aggregation */
         if( infeasible )
         {
            SCIPdebugMsg(scip, " -> infeasible aggregation\n");
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         /* count the aggregation */
         if( aggregated )
         {
            SCIPdebugMsg(scip, " -> feasible aggregation\n");
            (*naggrvars)++;
            success = TRUE;
         }
      }

      if( success )
      {
         /* apply fixings and aggregation to successfully rerun this presolving step */
         SCIP_CALL( applyFixings(scip, cons, &infeasible) );

         if( infeasible )
         {
            SCIPdebugMsg(scip, " -> infeasible fixing\n");
            *cutoff = TRUE;
            return SCIP_OKAY;
         }

         /* normalize constraint */
         SCIP_CALL( normalizeCons(scip, cons) );
      }
   }
   while( success );

   return SCIP_OKAY;
}



/** sorting method for constraint data, compares two variables on given indices, continuous variables will be sorted to
 *  the end and for all other variables the sortation will be in non-increasing order of their absolute value of the
 *  coefficients
 */
static
SCIP_DECL_SORTINDCOMP(consdataCompSim)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;
   SCIP_VARTYPE vartype1;
   SCIP_VARTYPE vartype2;
   SCIP_Real value;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);

   vartype1 = SCIPvarGetType(consdata->vars[ind1]);
   vartype2 = SCIPvarGetType(consdata->vars[ind2]);

   if( vartype1 == SCIP_VARTYPE_CONTINUOUS )
   {
      /* continuous varibles will be sorted to the back */
      if( vartype2 != vartype1 )
         return +1;
      /* both variables are continuous */
      else
         return 0;
   }
   /* continuous variables will be sorted to the back */
   else if( vartype2 == SCIP_VARTYPE_CONTINUOUS )
      return -1;

   value = REALABS(consdata->vals[ind2]) - REALABS(consdata->vals[ind1]);

   /* for all non-continuous variables, the variables are sorted after decreasing absolute coefficients */
   return (value > 0 ? +1 : (value < 0 ? -1 : 0));
}

/** tries to simplify coefficients and delete variables in ranged row of the form lhs <= a^Tx <= rhs, e.g. using the greatest
 *  common divisor
 *
 *  1. lhs <= a^Tx <= rhs, forall a_i >= lhs, a_i <= rhs, and forall pairs a_i + a_j > rhs then we can change this
 *     constraint to 1^Tx = 1
 */
static
SCIP_RETCODE rangedRowSimplify(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides           /**< pointer to store the amount of changed sides */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real minval;
   SCIP_Real secondminval;
   SCIP_Real maxval;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nvars;
   int v;

   /* we must not change a modifiable constraint in any way */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   if( SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* do not check empty or bound-constraints */
   if( nvars < 2 )
      return SCIP_OKAY;

   vals = consdata->vals;
   vars = consdata->vars;
   assert(vars != NULL);
   assert(vals != NULL);

   lhs = consdata->lhs;
   rhs = consdata->rhs;
   assert(!SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs));
   assert(!SCIPisNegative(scip, rhs));

   minval = SCIP_INVALID;
   secondminval = SCIP_INVALID;
   maxval = -SCIP_INVALID;

   for( v = nvars - 1; v >= 0; --v )
   {
      if( SCIPvarIsBinary(vars[v]) )
      {
         if( minval > vals[v] || minval == SCIP_INVALID ) /*lint !e777*/
         {
            secondminval = minval;
            minval = vals[v];
         }
         else if( secondminval > vals[v] || secondminval == SCIP_INVALID ) /*lint !e777*/
            secondminval = vals[v];

         if( maxval < vals[v] || maxval == -SCIP_INVALID ) /*lint !e777*/
            maxval = vals[v];
      }
      else
         break;
   }

   /* check if all variables are binary */
   if( v == -1 )
   {
      if( SCIPisEQ(scip, minval, maxval) && SCIPisEQ(scip, lhs, rhs) )
         return SCIP_OKAY;

      /* check if we can and need to choose exactly one binary variable */
      if( SCIPisGE(scip, minval, lhs) && SCIPisLE(scip, maxval, rhs) && SCIPisGT(scip, minval + secondminval, rhs) )
      {
         /* change all coefficients to 1.0 */
         for( v = nvars - 1; v >= 0; --v )
         {
            SCIP_CALL( chgCoefPos(scip, cons, v, 1.0) );
         }
         (*nchgcoefs) += nvars;

         /* replace old right and left hand side with 1.0 */
         SCIP_CALL( chgRhs(scip, cons, 1.0) );
         SCIP_CALL( chgLhs(scip, cons, 1.0) );
         (*nchgsides) += 2;
      }
   }

   return SCIP_OKAY;
}

/** tries to simplify coefficients and delete variables in constraints of the form lhs <= a^Tx <= rhs
 *  for equations @see rangedRowSimplify() will be called
 *
 *  there are several different coefficient reduction steps which will be applied
 *
 *  1. We try to determine parts of the constraint which will not change anything on (in-)feasibility of the constraint
 *
 *     e.g. 5x1 + 5x2 + 3z1 <= 8 => 3z1 is redundant if all x are binary and -2 < 3z1 <= 3
 *
 *  2. We try to remove redundant fractional parts in a constraint
 *
 *     e.g. 5.2x1 + 5.1x2 + 3x3 <= 8.3  => will be changed to 5x1 + 5x2 + 3x3 <= 8 if all x are binary
 *
 *  3. We are using the greatest common divisor for further reductions
 *
 *     e.g. 10x1 + 5y2 + 5x3 + 3x4 <= 15  => will be changed to 2x1 + y2 + x3 + x4 <= 3 if all xi are binary and y2 is
 *          integral
 */
static
SCIP_RETCODE simplifyInequalities(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   int*                  nchgcoefs,          /**< pointer to store the amount of changed coefficients */
   int*                  nchgsides           /**< pointer to store the amount of changed sides */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int* perm;
   SCIP_Real minactsub;
   SCIP_Real maxactsub;
   SCIP_Real siderest;
   SCIP_Real feastol;
   SCIP_Real newcoef;
   SCIP_Real absval;
   SCIP_Real side;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Longint restcoef;
   SCIP_Longint oldgcd;
   SCIP_Longint rest;
   SCIP_Longint gcd;
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;
   SCIP_Bool isminrelax;
   SCIP_Bool ismaxrelax;
   SCIP_Bool allcoefintegral;
   SCIP_Bool onlybin;
   SCIP_Bool hasrhs;
   SCIP_Bool haslhs;
   int oldnchgcoefs;
   int oldnchgsides;
   int foundbin;
   int candpos;
   int candpos2;
   int offsetv;
   int nvars;
   int v;
   int w;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   /* we must not change a modifiable constraint in any way */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   if( SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   nvars = consdata->nvars;

   /* do not check empty or bound-constraints */
   if( nvars <= 2 )
      return SCIP_OKAY;

   /* update maximal activity delta if necessary */
   if( consdata->maxactdelta == SCIP_INVALID ) /*lint !e777*/
      consdataRecomputeMaxActivityDelta(scip, consdata);

   assert(consdata->maxactdelta != SCIP_INVALID); /*lint !e777*/
   assert(!SCIPisFeasNegative(scip, consdata->maxactdelta));
   checkMaxActivityDelta(scip, consdata);

   /* @todo the following might be too hard, check which steps can be applied and what code must be corrected
    *       accordingly
    */
   /* can only work with valid non-infinity activities per variable */
   if( SCIPisInfinity(scip, consdata->maxactdelta) )
      return SCIP_OKAY;

   /* @todo: change the following: due to vartype changes, the status of the normalization can be wrong, need an event
    *        but the eventsystem seems to be full
    */
   consdata->normalized = FALSE;

   /* normalize constraint */
   SCIP_CALL( normalizeCons(scip, cons) );
   assert(nvars == consdata->nvars);

   if( !consdata->normalized )
      return SCIP_OKAY;

   lhs = consdata->lhs;
   rhs = consdata->rhs;
   assert(!SCIPisInfinity(scip, -lhs) || !SCIPisInfinity(scip, rhs));
   assert(!SCIPisNegative(scip, rhs));

   if( !SCIPisInfinity(scip, -lhs) )
      haslhs = TRUE;
   else
      haslhs = FALSE;

   if( !SCIPisInfinity(scip, rhs) )
      hasrhs = TRUE;
   else
      hasrhs = FALSE;

   SCIPdebug( oldnchgcoefs = *nchgcoefs; )
   SCIPdebug( oldnchgsides = *nchgsides; )

   /* @todo also work on ranged rows */
   if( haslhs && hasrhs )
   {
      SCIP_CALL( rangedRowSimplify(scip, cons, nchgcoefs, nchgsides ) );

      return SCIP_OKAY;
   }
   assert(haslhs != hasrhs);

   /* if we have a normalized inequality (not ranged) the one side should be positive, @see normalizeCons() */
   assert(!hasrhs || !SCIPisNegative(scip, rhs));
   assert(!haslhs || !SCIPisNegative(scip, lhs));

   /* get temporary memory to store the sorted permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nvars) );

   /* call sorting method, order continuous variables to the end and all other variables after non-increasing absolute
    * value of their coefficients
    */
   SCIPsort(perm, consdataCompSim, (void*)consdata, nvars);

   /* perform sorting after permutation array */
   permSortConsdata(consdata, perm, nvars);
   consdata->sorted = FALSE;
   consdata->binvarssorted = FALSE;

   vars = consdata->vars;
   vals = consdata->vals;
   assert(vars != NULL);
   assert(vals != NULL);
   assert(consdata->validmaxabsval ? (SCIPisFeasEQ(scip, consdata->maxabsval, REALABS(vals[0])) || SCIPvarGetType(vars[nvars - 1]) == SCIP_VARTYPE_CONTINUOUS) : TRUE);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &perm);

   /* only check constraints with at least two non continuous variables */
   if( SCIPvarGetType(vars[1]) == SCIP_VARTYPE_CONTINUOUS )
      return SCIP_OKAY;

   /* do not process constraints when all coefficients are 1.0 */
   if( SCIPisEQ(scip, REALABS(vals[0]), 1.0) && ((hasrhs && SCIPisIntegral(scip, rhs)) || (haslhs && SCIPisIntegral(scip, lhs))) )
      return SCIP_OKAY;

   feastol = SCIPfeastol(scip);

   SCIPdebugMsg(scip, "starting simplification of coefficients\n");
   SCIPdebugPrintCons(scip, cons, NULL);

   /* get global activities */
   consdataGetGlbActivityBounds(scip, consdata, FALSE, &minactsub, &maxactsub,
      &isminrelax, &ismaxrelax, &isminsettoinfinity, &ismaxsettoinfinity);

   /* cannot work with infinite activities */
   if( isminsettoinfinity || ismaxsettoinfinity )
      return SCIP_OKAY;

   assert(!isminrelax);
   assert(!ismaxrelax);
   assert(maxactsub > minactsub);
   assert(!SCIPisInfinity(scip, -minactsub));
   assert(!SCIPisInfinity(scip, maxactsub));

   v = 0;
   offsetv = -1;
   side = haslhs ? lhs : rhs;

   /* we now determine coefficients as large as the side of the constraint to retrieve a better reduction where we
    * do not need to look at the large coefficients
    *
    * e.g.  all x are binary, z are positive integer
    *       c1: +5x1 + 5x2 + 3x3 + 3x4 + x5 >= 5   (x5 is redundant and does not change (in-)feasibility of this constraint)
    *       c2: +4x1 + 4x2 + 3x3 + 3x4 + x5 >= 4   (gcd (without the coefficient of x5) after the large coefficients is 3
    *       c3: +30x1 + 29x2 + 14x3 + 14z1 + 7x5 + 7x6 <= 30 (gcd (without the coefficient of x2) after the large coefficients is 7
    *
    *       can be changed to
    *
    *       c1: +6x1 + 6x2 + 3x3 + 3x4 >= 6        (will be changed to c1: +2x1 + 2x2 + x3 + x4 >= 2)
    *       c2: +6x1 + 6x2 + 3x3 + 3x4 + 3x5 >= 6  (will be changed to c2: +2x1 + 2x2 + x3 + x4 + x5 >= 2)
    *       c3: +28x1 + 28x2 + 14x3 + 14z1 + 7x5 + 7x6 <= 28 (will be changed to c3: +4x1 + 4x2 + 2x3 + 2z1 + x5 + x6 <= 4)
    */

   /* if the minimal activity is negative and we found more than one variable with a coefficient bigger than the left
    * hand side, we cannot apply the extra reduction step and need to reset v
    *
    * e.g. 7x1 + 7x2 - 4x3 - 4x4 >= 7 => xi = 1 for all i is not a solution, but if we would do a change on the
    *      coefficients due to the gcd on the "small" coefficients we would get 8x1 + 8x2 - 4x3 - 4x4 >= 8 were xi = 1
    *      for all i is a solution
    *
    * also redundancy of variables would not be correctly determined in such a case
    */
   if( nvars > 2 && SCIPisEQ(scip, vals[0], side) && !SCIPisNegative(scip, minactsub) )
   {
      v = 1;

      while( v < nvars && SCIPisEQ(scip, side, vals[v]) )
      {
         /* if we have integer variable with "side"-coefficients but also with a lower bound greater than 0 we stop this
          * extra step, which might have worked
          */
         if( SCIPvarGetLbGlobal(vars[v]) > 0.5 )
         {
            v = 0;
            break;
         }

         ++v;
      }

      /* easy and quick fix: if all coefficients were equal to the side, we cannot apply further simplifications */
      /* todo find numerically stable normalization conditions to scale this cons to have coefficients almost equal to 1 */
      if( v == nvars )
         return SCIP_OKAY;

      /* cannot work with continuous variables which have a big coefficient */
      if( v > 0 && SCIPvarGetType(vars[v - 1]) == SCIP_VARTYPE_CONTINUOUS )
         return SCIP_OKAY;

      /* big negative coefficient, do not try to use the extra coefficient reduction step */
      if( SCIPisEQ(scip, side, -vals[v]) )
         v = 0;

      /* all but one variable are processed or the next variable is continuous we cannot perform the extra coefficient
       * reduction
       */
      if( v == nvars - 1 || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
         v = 0;

      if( v > 0 )
      {
         assert(v < nvars);

         offsetv = v - 1;

         for( w = 0; w < v; ++w )
         {
            lb = SCIPvarGetLbGlobal(vars[w]);
            ub = SCIPvarGetUbGlobal(vars[w]);

            assert(vals[w] > 0);

            /* update residual activities */
            maxactsub -= ub * vals[w];
            minactsub -= lb * vals[w];
            assert(maxactsub > minactsub);
         }
      }
   }

   /* find and remove redundant variables which do not interact with the (in-)feasibility of this constraint
    *
    * e.g. let all x are binary and y1 is continuous with bounds [-3,1] then we can reduce
    *
    *        15x1 + 15x2 + 7x3 + 3x4 + y1 <= 26
    * to
    *        15x1 + 15x2 <= 26 <=> x1 + x2 <= 1
    */
   if( nvars > 2 && SCIPisIntegral(scip, vals[v]) )
   {
      SCIP_Bool redundant = FALSE;

      gcd = (SCIP_Longint)(REALABS(vals[v]) + feastol);
      assert(gcd >= 1);

      if( v == 0 )
      {
         lb = SCIPvarGetLbGlobal(vars[0]);
         ub = SCIPvarGetUbGlobal(vars[0]);

         /* update residual activities */
         if( vals[0] > 0 )
         {
            maxactsub -= ub * vals[0];
            minactsub -= lb * vals[0];
         }
         else
         {
            maxactsub -= lb * vals[0];
            minactsub -= ub * vals[0];
         }
         assert(maxactsub > minactsub);
         ++v;
      }

      siderest = -SCIP_INVALID;
      allcoefintegral = TRUE;

      /* check if some variables always fit into the given constraint */
      for( ; v < nvars - 1; ++v )
      {
         if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            break;

         if( !SCIPisIntegral(scip, vals[v]) )
         {
            allcoefintegral = FALSE;
            break;
         }

         /* calculate greatest common divisor for all general and binary variables */
         gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[v]) + feastol));

         if( gcd == 1 )
            break;

         lb = SCIPvarGetLbGlobal(vars[v]);
         ub = SCIPvarGetUbGlobal(vars[v]);

         assert(!SCIPisInfinity(scip, -lb));
         assert(!SCIPisInfinity(scip, ub));

         /* update residual activities */
         if( vals[v] > 0 )
         {
            maxactsub -= ub * vals[v];
            minactsub -= lb * vals[v];
         }
         else
         {
            maxactsub -= lb * vals[v];
            minactsub -= ub * vals[v];
         }
         assert(SCIPisGE(scip, maxactsub, minactsub));

         if( hasrhs )
         {
            /* determine the remainder of the right hand side and the gcd */
            siderest = rhs - SCIPfeasFloor(scip, rhs/gcd) * gcd;
         }
         else
         {
            /* determine the remainder of the left hand side and the gcd */
            siderest = lhs - SCIPfeasFloor(scip, lhs/gcd) * gcd;
            if( SCIPisZero(scip, siderest) )
               siderest = gcd;
         }

         /* early termination if the activities deceed the gcd */
         if( (offsetv == -1 && hasrhs && maxactsub <= siderest && SCIPisFeasGT(scip, minactsub, siderest - gcd)) || (haslhs && SCIPisFeasLT(scip, maxactsub, siderest) && minactsub >= siderest - gcd) )
         {
            redundant = TRUE;
            break;
         }
      }
      assert(v < nvars || (offsetv >= 0 && gcd > 1));

      if( !redundant )
      {
         if( hasrhs )
         {
            /* determine the remainder of the right hand side and the gcd */
            siderest = rhs - SCIPfeasFloor(scip, rhs/gcd) * gcd;
         }
         else
         {
            /* determine the remainder of the left hand side and the gcd */
            siderest = lhs - SCIPfeasFloor(scip, lhs/gcd) * gcd;
            if( SCIPisZero(scip, siderest) )
               siderest = gcd;
         }
      }
      else
         ++v;

      SCIPdebugMsg(scip, "stopped at pos %d (of %d), subactivities [%g, %g], redundant = %u, hasrhs = %u, siderest = %g, gcd = %" SCIP_LONGINT_FORMAT ", offset position for 'side' coefficients = %d\n",
            v, nvars, minactsub, maxactsub, redundant, hasrhs, siderest, gcd, offsetv);

      /* check if we can remove redundant variables */
      if( v < nvars && (redundant ||
            (offsetv == -1 && hasrhs && maxactsub <= siderest && SCIPisFeasGT(scip, minactsub, siderest - gcd)) ||
            (haslhs && SCIPisFeasLT(scip, maxactsub, siderest) && minactsub >= siderest - gcd)) )
      {
         SCIP_Real oldcoef;

         /* double check the redundancy */
#ifndef NDEBUG
         SCIP_Real tmpminactsub = 0.0;
         SCIP_Real tmpmaxactsub = 0.0;

         /* recompute residual activities */
         for( w = v; w < nvars; ++w )
         {
            lb = SCIPvarGetLbGlobal(vars[w]);
            ub = SCIPvarGetUbGlobal(vars[w]);

            assert(!SCIPisInfinity(scip, -lb));
            assert(!SCIPisInfinity(scip, ub));

            /* update residual activities */
            if( vals[w] > 0 )
            {
               tmpmaxactsub += ub * vals[w];
               tmpminactsub += lb * vals[w];
            }
            else
            {
               tmpmaxactsub += lb * vals[w];
               tmpminactsub += ub * vals[w];
            }
            assert(tmpmaxactsub >= tmpminactsub);
         }

         if( hasrhs )
         {
            assert(offsetv == -1);

            /* determine the remainder of the right hand side and the gcd */
            siderest = rhs - SCIPfeasFloor(scip, rhs/gcd) * gcd;
         }
         else
         {
            /* determine the remainder of the left hand side and the gcd */
            siderest = lhs - SCIPfeasFloor(scip, lhs/gcd) * gcd;
            if( SCIPisZero(scip, siderest) )
               siderest = gcd;
         }

         /* does the redundancy really is fulfilled */
         assert((hasrhs && SCIPisLE(scip, tmpmaxactsub, siderest) && tmpminactsub > siderest - gcd) || (haslhs && tmpmaxactsub < siderest && SCIPisGE(scip, tmpminactsub, siderest - gcd)));
#endif

         SCIPdebugMsg(scip, "removing %d last variables from constraint <%s>, because they never change anything on the feasibility of this constraint\n",
            nvars - v, SCIPconsGetName(cons));

         /* remove redundant variables */
         for( w = nvars - 1; w >= v; --w )
         {
            SCIP_CALL( delCoefPos(scip, cons, w) );
         }
         (*nchgcoefs) += (nvars - v);

         assert(w >= 0);

         oldcoef = vals[w];

         /* normalize constraint */
         SCIP_CALL( normalizeCons(scip, cons) );
         assert(vars == consdata->vars);
         assert(vals == consdata->vals);
         assert(w < consdata->nvars);

         /* compute new greatest common divisor due to normalization */
         gcd = (SCIP_Longint)(gcd / (oldcoef/vals[w]) + feastol);
         assert(gcd >= 1);

         /* update side */
         if( hasrhs )
         {
            /* replace old with new right hand side */
            SCIP_CALL( chgRhs(scip, cons, SCIPfeasFloor(scip, consdata->rhs)) );
            rhs = consdata->rhs;
         }
         else
         {
            if( SCIPisFeasGT(scip, oldcoef/vals[w], 1.0) )
            {
               SCIP_CALL( chgLhs(scip, cons, SCIPfeasCeil(scip, consdata->lhs)) );
               lhs = consdata->lhs;
            }
            else
               assert(offsetv == -1 || SCIPisEQ(scip, vals[offsetv], consdata->lhs));
         }
         ++(*nchgsides);

         assert(!hasrhs || !SCIPisNegative(scip, rhs));
         assert(!haslhs || !SCIPisNegative(scip, lhs));

         /* get new constraint data */
         nvars = consdata->nvars;
         assert(nvars >= 2);

         allcoefintegral = TRUE;

#ifndef NDEBUG
         /* check integrality */
         for( w = offsetv + 1; w < nvars; ++w )
         {
            assert(SCIPisIntegral(scip, vals[w]));
         }
#endif
         SCIPdebugPrintCons(scip, cons, NULL);
      }

      /* try to find a better gcd, when having large coefficients */
      if( offsetv >= 0 && gcd == 1 )
      {
         /* calculate greatest common divisor for all general variables */
         gcd = (SCIP_Longint)(REALABS(vals[nvars - 1]) + feastol);

         if( gcd > 1 )
         {
            gcd = -1;
            candpos = -1;

            for( v = nvars - 1; v > offsetv; --v )
            {
               assert(!SCIPisZero(scip, vals[v]));
               if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
                  break;

               if( !SCIPisIntegral(scip, vals[v]) )
               {
                  allcoefintegral = FALSE;
                  break;
               }

               oldgcd = gcd;

               if( gcd == -1 )
               {
                  gcd = (SCIP_Longint)(REALABS(vals[v]) + feastol);
                  assert(gcd >= 1);
               }
               else
               {
                  /* calculate greatest common divisor for all general and binary variables */
                  gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[v]) + feastol));
               }

               /* if the greatest commmon divisor has become 1, we might have found the possible coefficient to change or we
                * can stop searching
                */
               if( gcd == 1 )
               {
                  if( !SCIPvarIsBinary(vars[v]) )
                     break;

                  /* found candidate */
                  if( candpos == -1 )
                  {
                     gcd = oldgcd;
                     candpos = v;
                  }
                  /* two different binary variables lead to a gcd of one, so we cannot change a coefficient */
                  else
                     break;
               }
            }
            assert(v > offsetv || candpos > offsetv);
         }
         else
            candpos = -1;
      }
      else
         candpos = nvars - 1;

      /* check last coefficient for integrality */
      if( gcd > 1 && allcoefintegral && !redundant )
      {
         if( !SCIPisIntegral(scip, vals[nvars - 1]) )
            allcoefintegral = FALSE;
      }

      /* check for further necessary coefficient adjustments */
      if( offsetv >= 0 && gcd > 1 && allcoefintegral )
      {
         assert(offsetv + 1 < nvars);
         assert(0 <= candpos && candpos < nvars);

         if( SCIPvarGetType(vars[candpos]) != SCIP_VARTYPE_CONTINUOUS )
         {
            SCIP_Bool notchangable = FALSE;

#ifndef NDEBUG
            /* check integrality */
            for( w = offsetv + 1; w < nvars; ++w )
            {
               assert(SCIPisIntegral(scip, vals[w]));
            }
#endif

            if( vals[candpos] > 0 && SCIPvarIsBinary(vars[candpos]) &&
               SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[candpos]) + feastol)) < gcd )
            {
               /* determine the remainder of the side and the gcd */
               if( hasrhs )
                  rest = ((SCIP_Longint)(rhs + feastol)) % gcd;
               else
                  rest = ((SCIP_Longint)(lhs + feastol)) % gcd;
               assert(rest >= 0);
               assert(rest < gcd);

               /* determine the remainder of the coefficient candidate and the gcd */
               restcoef = ((SCIP_Longint)(vals[candpos] + feastol)) % gcd;
               assert(restcoef >= 1);
               assert(restcoef < gcd);

               if( hasrhs )
               {
                  /* calculate new coefficient */
                  if( restcoef > rest )
                     newcoef = vals[candpos] - restcoef + gcd;
                  else
                     newcoef = vals[candpos] - restcoef;
               }
               else
               {
                  /* calculate new coefficient */
                  if( rest == 0 || restcoef < rest )
                     newcoef = vals[candpos] - restcoef;
                  else
                     newcoef = vals[candpos] - restcoef + gcd;
               }


               /* done */

               /* new coeffcient must not be zero if we would loose the implication that a variable needs to be 0 if
                * another with the big coefficient was set to 1
                */
               if( hasrhs && SCIPisZero(scip, newcoef) )
               {
                  notchangable = TRUE;
               }
               else if( SCIPisZero(scip, newcoef) )
               {
                  /* delete old redundant coefficient */
                  SCIP_CALL( delCoefPos(scip, cons, candpos) );
                  ++(*nchgcoefs);
               }
               else
               {
                  /* replace old with new coefficient */
                  SCIP_CALL( chgCoefPos(scip, cons, candpos, newcoef) );
                  ++(*nchgcoefs);
               }
            }
            else if( vals[candpos] < 0 || !SCIPvarIsBinary(vars[candpos]) )
            {
               gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[candpos]) + feastol));
            }

            /* correct side and big coefficients */
            if( (!notchangable && hasrhs && ((!SCIPisFeasIntegral(scip, rhs) || SCIPcalcGreComDiv(gcd, (SCIP_Longint)(rhs + feastol)) < gcd) && (SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[candpos]) + feastol)) == gcd))) ||
               ( haslhs && (!SCIPisFeasIntegral(scip, lhs) || SCIPcalcGreComDiv(gcd, (SCIP_Longint)(lhs + feastol)) < gcd) && (SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[candpos]) + feastol)) == gcd)) )
            {
               if( haslhs )
               {
                  newcoef = (SCIP_Real)((SCIP_Longint)(SCIPfeasCeil(scip, lhs/gcd) * gcd + feastol));

                  SCIP_CALL( chgLhs(scip, cons, newcoef) );
                  ++(*nchgsides);
               }
               else
               {
                  assert(hasrhs);
                  newcoef = (SCIP_Real)((SCIP_Longint)(SCIPfeasFloor(scip, rhs/gcd) * gcd + feastol));

                  SCIP_CALL( chgRhs(scip, cons, newcoef) );
                  ++(*nchgsides);
               }

               /* correct coefficients up front */
               for( w = offsetv; w >= 0; --w )
               {
                  assert(vals[w] > 0);

                  SCIP_CALL( chgCoefPos(scip, cons, w, newcoef) );
               }
               (*nchgcoefs) += (offsetv + 1);
            }

            if( !notchangable )
            {
               /* normalize constraint */
               SCIP_CALL( normalizeCons(scip, cons) );
               assert(vars == consdata->vars);
               assert(vals == consdata->vals);

               /* get new constraint data */
               nvars = consdata->nvars;
               assert(nvars >= 2);

               SCIPdebugPrintCons(scip, cons, NULL);

               lhs = consdata->lhs;
               rhs = consdata->rhs;
               assert(!hasrhs || !SCIPisNegative(scip, rhs));
               assert(!haslhs || !SCIPisNegative(scip, lhs));
            }
         }
      }
   }

   /* @todo we still can remove continuous variables if they are redundant due to the non-integrality argument */
   /* no continuous variables are left over */
   if( SCIPvarGetType(vars[nvars - 1]) == SCIP_VARTYPE_CONTINUOUS )
      return SCIP_OKAY;

   onlybin = TRUE;
   allcoefintegral = TRUE;
   /* check if all variables are of binary type */
   for( v = nvars - 1; v >= 0; --v )
   {
      if( !SCIPvarIsBinary(vars[v]) )
         onlybin = FALSE;
      if( !SCIPisIntegral(scip, vals[v]) )
         allcoefintegral = FALSE;
   }

   /* check if the non-integrality part of all integral variables is smaller than the non-inegrality part of the right
    * hand side or bigger than the left hand side respectively, so we can make all of them integral
    *
    * @todo there are some steps missing ....
    */
   if( (hasrhs && !SCIPisFeasIntegral(scip, rhs)) || (haslhs && !SCIPisFeasIntegral(scip, lhs)) )
   {
      SCIP_Real val;
      SCIP_Real newval;
      SCIP_Real frac = 0.0;
      SCIP_Bool found = FALSE;

      if( hasrhs )
      {
         if( allcoefintegral )
         {
            /* replace old with new right hand side */
            SCIP_CALL( chgRhs(scip, cons, SCIPfloor(scip, rhs)) );
            ++(*nchgsides);
         }
         else
         {
            siderest = rhs - SCIPfloor(scip, rhs);

            /* try to round down all non-integral coefficients */
            for( v = nvars - 1; v >= 0; --v )
            {
               val = vals[v];

               /* add up all possible fractional parts */
               if( !SCIPisIntegral(scip, val) )
               {
                  lb = SCIPvarGetLbGlobal(vars[v]);
                  ub = SCIPvarGetUbGlobal(vars[v]);

                  /* at least one bound need to be at zero */
                  if( !onlybin && !SCIPisFeasZero(scip, lb) && !SCIPisFeasZero(scip, ub) )
                     return SCIP_OKAY;

                  /* swap bounds for 'standard' form */
                  if( !SCIPisFeasZero(scip, lb) )
                  {
                     ub = lb;
                     val *= -1;
                  }

                  found = TRUE;

                  frac += (val - SCIPfloor(scip, val)) * ub;

                  /* if we exceed the fractional part of the right hand side, we cannot tighten the coefficients
                   *
                   * e.g. 1.1x1 + 1.1x2 + 1.4x3 + 1.02x4 <= 2.4, here we cannot floor all fractionals because
                   *      x3, x4 set to 1 would be infeasible but feasible after flooring
                   */
                  if( SCIPisGT(scip, frac, siderest) )
                     return SCIP_OKAY;
               }
            }
            assert(v == -1);

            SCIPdebugMsg(scip, "rounding all non-integral coefficients and the right hand side down\n");

            /* round rhs and coefficients to integral values */
            if( found )
            {
               for( v = nvars - 1; v >= 0; --v )
               {
                  val = vals[v];

                  /* add the whole fractional part */
                  if( !SCIPisIntegral(scip, val) )
                  {
                     lb = SCIPvarGetLbGlobal(vars[v]);

                     if( SCIPisFeasZero(scip, lb) )
                        newval = SCIPfloor(scip, val);
                     else
                        newval = SCIPceil(scip, val);

                     if( SCIPisZero(scip, newval) )
                     {
                        /* delete old redundant coefficient */
                        SCIP_CALL( delCoefPos(scip, cons, v) );
                        ++(*nchgcoefs);
                     }
                     else
                     {
                        /* replace old with new coefficient */
                        SCIP_CALL( chgCoefPos(scip, cons, v, newval) );
                        ++(*nchgcoefs);
                     }
                  }
               }
            }

            /* replace old with new right hand side */
            SCIP_CALL( chgRhs(scip, cons, SCIPfloor(scip, rhs)) );
            ++(*nchgsides);
         }
      }
      else
      {
         if( allcoefintegral )
         {
            /* replace old with new left hand side */
            SCIP_CALL( chgLhs(scip, cons, SCIPceil(scip, lhs)) );
            ++(*nchgsides);
         }
         else
         {
            /* cannot floor left hand side to zero */
            if( SCIPisLT(scip, lhs, 1.0) )
               return SCIP_OKAY;

            siderest = lhs - SCIPfloor(scip, lhs);

            /* try to round down all non-integral coefficients */
            for( v = nvars - 1; v >= 0; --v )
            {
               val = vals[v];

               /* add up all possible fractional parts */
               if( !SCIPisIntegral(scip, val) )
               {
                  lb = SCIPvarGetLbGlobal(vars[v]);
                  ub = SCIPvarGetUbGlobal(vars[v]);

                  /* at least one bound need to be at zero */
                  if( !SCIPisFeasZero(scip, lb) && !SCIPisFeasZero(scip, ub) )
                     return SCIP_OKAY;

                  /* swap bounds for 'standard' form */
                  if( !SCIPisFeasZero(scip, lb) )
                  {
                     ub = lb;
                     val *= -1;
                  }

                  /* cannot floor to zero */
                  if( SCIPisLT(scip, val, 1.0) )
                     return SCIP_OKAY;

                  /* the fractional part on each variable need to exceed the fractional part on the left hand side */
                  if( SCIPisLT(scip, val - SCIPfloor(scip, val), siderest) )
                     return SCIP_OKAY;

                  found = TRUE;

                  frac += (val - SCIPfloor(scip, val)) * ub;

                  /* if we exceed the fractional part of the left hand side plus one by summing up all maximal
                   * fractional parts of the variables, we cannot tighten the coefficients
                   *
                   * e.g. 4.3x1 + 1.3x2 + 1.3x3 + 1.6x4 >= 4.2, here we cannot floor all fractionals because
                   *      x2-x4 set to 1 would be feasible but not after flooring
                   */
                  if( SCIPisGE(scip, frac, 1 + siderest) )
                     return SCIP_OKAY;
               }
               /* all coefficients need to be integral, otherwise we might do an invalid reduction */
               else
                  return SCIP_OKAY;
            }
            assert(v == -1);

            SCIPdebugMsg(scip, "rounding all non-integral coefficients and the left hand side down\n");

            /* round lhs and coefficients to integral values */
            if( found )
            {
               for( v = nvars - 1; v >= 0; --v )
               {
                  val = vals[v];

                  /* add the whole fractional part */
                  if( !SCIPisIntegral(scip, val) )
                  {
                     lb = SCIPvarGetLbGlobal(vars[v]);

                     if( SCIPisFeasZero(scip, lb) )
                        newval = SCIPfloor(scip, val);
                     else
                        newval = SCIPceil(scip, val);

                     if( SCIPisZero(scip, newval) )
                     {
                        /* delete old redundant coefficient */
                        SCIP_CALL( delCoefPos(scip, cons, v) );
                        ++(*nchgcoefs);
                     }
                     else
                     {
                        /* replace old with new coefficient */
                        SCIP_CALL( chgCoefPos(scip, cons, v, newval) );
                        ++(*nchgcoefs);
                     }
                  }
               }
            }

            /* replace old with new left hand side */
            SCIP_CALL( chgLhs(scip, cons, SCIPfloor(scip, lhs)) );
            ++(*nchgsides);
         }
      }

      /* normalize constraint */
      SCIP_CALL( normalizeCons(scip, cons) );
      assert(vars == consdata->vars);
      assert(vals == consdata->vals);

      rhs = consdata->rhs;
      lhs = consdata->lhs;

      assert(!hasrhs || !SCIPisNegative(scip, rhs));
      assert(!haslhs || !SCIPisNegative(scip, lhs));

      SCIPdebugPrintCons(scip, cons, NULL);

      nvars = consdata->nvars;
      if( nvars < 2 )
         return SCIP_OKAY;

      allcoefintegral = TRUE;
#ifndef NDEBUG
      /* debug check if all coefficients are really integral */
      for( v = nvars - 1; v >= 0; --v )
         assert(SCIPisIntegral(scip, vals[v]));
#endif
   }

   /* @todo following can also work on non integral coefficients, need more investigation */
   /* only check constraints with integral coefficients on all integral variables */
   if( !allcoefintegral )
      return SCIP_OKAY;

   /* we want to avoid numerical troubles, therefore we do not change non-integral sides */
   if( (hasrhs && !SCIPisIntegral(scip, rhs)) || (haslhs && !SCIPisIntegral(scip, lhs)) )
      return SCIP_OKAY;

   /* maximal absolute value of coefficients in constraint is one, so we cannot tighten it further */
   if( SCIPisEQ(scip, REALABS(vals[0]), 1.0) )
      return SCIP_OKAY;

   /* stop if the last coeffcients is one in absolute value and the variable is not binary */
   if( !SCIPvarIsBinary(vars[nvars - 1]) && SCIPisEQ(scip, REALABS(vals[nvars - 1]), 1.0) )
      return SCIP_OKAY;

   assert(nvars >= 2);

   /* start gcd procedure for all variables */
   do
   {
      SCIPdebug( oldnchgcoefs = *nchgcoefs; )
      SCIPdebug( oldnchgsides = *nchgsides; )

      /* stop if we have two coeffcients which are one in absolute value */
      if( SCIPisEQ(scip, REALABS(vals[nvars - 1]), 1.0) && SCIPisEQ(scip, REALABS(vals[nvars - 2]), 1.0) )
         return SCIP_OKAY;

      gcd = -1;

      /* calculate greatest common divisor over all integer variables */
      if( !onlybin )
      {
         foundbin = -1;

         for( v = nvars - 1; v >= 0; --v )
         {
            assert(!SCIPisZero(scip, vals[v]));
            assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_CONTINUOUS);

            if( SCIPvarIsBinary(vars[v]) )
            {
               if( foundbin == -1 )
                  foundbin = v;
               continue;
            }

            absval = REALABS(vals[v]);
            /* arithmetic precision can lead to the absolute value only being integral up to feasibility tolerance,
             * even though the value itself is feasible up to epsilon, but since we add feastol later, this is enough
             */
            assert(SCIPisFeasIntegral(scip, absval));

            if( gcd == -1 )
            {
               gcd = (SCIP_Longint)(absval + feastol);
               assert(gcd >= 1);
            }
            else
            {
               /* calculate greatest common divisor for all general variables */
               gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(absval + feastol));
            }
            if( gcd == 1 )
               break;
         }
      }
      else
         foundbin = nvars - 1;

      /* we need at least one binary variable and a gcd greater than 1 to try to perform further coefficient changes */
      if( gcd == 1 || foundbin == -1)
         return SCIP_OKAY;

      assert((onlybin && gcd == -1) || (!onlybin && gcd > 1));

      candpos = -1;
      candpos2 = -1;

      /* calculate greatest common divisor over all integer and binary variables and determine the candidate where we might
       * change the coefficient
       */
      for( v = foundbin; v >= 0; --v )
      {
         if( onlybin || SCIPvarIsBinary(vars[v]) )
         {
            absval = REALABS(vals[v]);
            /* arithmetic precision can lead to the absolute value only being integral up to feasibility tolerance,
             * even though the value itself is feasible up to epsilon, but since we add feastol later, this is enough
             */
            assert(SCIPisFeasIntegral(scip, absval));

            oldgcd = gcd;

            if( gcd == -1 )
            {
               gcd = (SCIP_Longint)(REALABS(vals[v]) + feastol);
               assert(gcd >= 1);
            }
            else
            {
               /* calculate greatest common divisor for all general and binary variables */
               gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[v]) + feastol));
            }

            /* if the greatest commmon divisor has become 1, we might have found the possible coefficient to change or we
             * can terminate
             */
            if( gcd == 1 )
            {
               /* found candidate */
               if( candpos == -1 )
               {
                  gcd = oldgcd;
                  candpos = v;

                  /* if we have only binary variables and both first coefficients have a gcd of 1, both are candidates for
                   * the coefficient change
                   */
                  if( onlybin && v == foundbin - 1 )
                     candpos2 = foundbin;
               }
               /* two different binary variables lead to a gcd of one, so we cannot change a coefficient */
               else
               {
                  if( onlybin && candpos == v + 1 && candpos2 == v + 2 )
                  {
                     assert(candpos2 == nvars - 1);

                     /* take new candidates */
                     candpos = candpos2;

                     /* recalculate gcd from scratch */
                     gcd = (SCIP_Longint)(REALABS(vals[v+1]) + feastol);
                     assert(gcd >= 1);

                     /* calculate greatest common divisor for all general and binary variables */
                     gcd = SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(vals[v]) + feastol));
                     if( gcd == 1 )
                        return SCIP_OKAY;
                  }
                  else
                     /* cannot determine a possible coefficient for reduction */
                     return SCIP_OKAY;
               }
            }
         }
      }
      assert(gcd >= 2);

      /* we should have found one coefficient, that led to a gcd of 1, otherwise we could normalize the constraint
       * further
       */
      assert(candpos >= 0 && candpos < nvars);

      /* all variables and all coefficients are integral, so the side should be too */
      assert((hasrhs && SCIPisIntegral(scip, rhs)) || (haslhs && SCIPisIntegral(scip, lhs)));

      /* check again, if we have a normalized inequality (not ranged) the one side should be positive,
       * @see normalizeCons()
       */
      assert(!hasrhs || !SCIPisNegative(scip, rhs));
      assert(!haslhs || !SCIPisNegative(scip, lhs));

      /* determine the remainder of the side and the gcd */
      if( hasrhs )
         rest = ((SCIP_Longint)(rhs + feastol)) % gcd;
      else
         rest = ((SCIP_Longint)(lhs + feastol)) % gcd;
      assert(rest >= 0);
      assert(rest < gcd);

      /* determine the remainder of the coefficient candidate and the gcd */
      if( vals[candpos] < 0 )
      {
         restcoef = ((SCIP_Longint)(vals[candpos] - feastol)) % gcd;
         assert(restcoef <= -1);
         restcoef += gcd;
      }
      else
         restcoef = ((SCIP_Longint)(vals[candpos] + feastol)) % gcd;
      assert(restcoef >= 1);
      assert(restcoef < gcd);

      if( hasrhs )
      {
         if( rest > 0 )
         {
            /* replace old with new right hand side */
            SCIP_CALL( chgRhs(scip, cons, rhs - rest) );
            ++(*nchgsides);
         }

         /* calculate new coefficient */
         if( restcoef > rest )
            newcoef = vals[candpos] - restcoef + gcd;
         else
            newcoef = vals[candpos] - restcoef;
      }
      else
      {
         if( rest > 0 )
         {
            /* replace old with new left hand side */
            SCIP_CALL( chgLhs(scip, cons, lhs - rest + gcd) );
            ++(*nchgsides);
         }

         /* calculate new coefficient */
         if( rest == 0 || restcoef < rest )
            newcoef = vals[candpos] - restcoef;
         else
            newcoef = vals[candpos] - restcoef + gcd;
      }
      assert(SCIPisZero(scip, newcoef) || SCIPcalcGreComDiv(gcd, (SCIP_Longint)(REALABS(newcoef) + feastol)) == gcd);

      SCIPdebugMsg(scip, "gcd = %" SCIP_LONGINT_FORMAT ", rest = %" SCIP_LONGINT_FORMAT ", restcoef = %" SCIP_LONGINT_FORMAT "; changing coef of variable <%s> to %g and %s by %" SCIP_LONGINT_FORMAT "\n", gcd, rest, restcoef, SCIPvarGetName(vars[candpos]), newcoef, hasrhs ? "reduced rhs" : "increased lhs", hasrhs ? rest : (rest > 0 ? gcd - rest : 0));

      if( SCIPisZero(scip, newcoef) )
      {
         /* delete redundant coefficient */
         SCIP_CALL( delCoefPos(scip, cons, candpos) );
      }
      else
      {
         /* replace old with new coefficient */
         SCIP_CALL( chgCoefPos(scip, cons, candpos, newcoef) );
      }
      ++(*nchgcoefs);

      /* now constraint can be normalized, might be directly done by dividing it by the gcd */
      SCIP_CALL( normalizeCons(scip, cons) );
      assert(vars == consdata->vars);
      assert(vals == consdata->vals);

      SCIPdebugPrintCons(scip, cons, NULL);

      rhs = consdata->rhs;
      lhs = consdata->lhs;
      assert(!hasrhs || !SCIPisNegative(scip, rhs));
      assert(!haslhs || !SCIPisNegative(scip, lhs));

      nvars = consdata->nvars;

      SCIPdebugMsg(scip, "we did %d coefficient changes and %d side changes on constraint %s when applying one round of the gcd algorithm\n", *nchgcoefs - oldnchgcoefs, *nchgsides - oldnchgsides, SCIPconsGetName(cons));
   }
   while( nvars >= 2 );

   return SCIP_OKAY;
}


/* tries to aggregate an (in)equality and an equality in order to decrease the number of variables in the (in)equality:
 *   cons0 := a * cons0 + b * cons1,
 * where a = val1[v] and b = -val0[v] for common variable v which removes most variable weight;
 * for numerical stability, we will only accept integral a and b;
 * the variable weight is a weighted sum over all included variables, where each binary variable weighs BINWEIGHT,
 * each integer or implicit integer variable weighs INTWEIGHT and each continuous variable weighs CONTWEIGHT
 */
static
SCIP_RETCODE aggregateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< (in)equality to modify */
   SCIP_CONS*            cons1,              /**< equality to use for aggregation of cons0 */
   int*                  commonidx0,         /**< array with indices of variables in cons0, that appear also in cons1 */
   int*                  commonidx1,         /**< array with indices of variables in cons1, that appear also in cons0 */
   int*                  diffidx0minus1,     /**< array with indices of variables in cons0, that don't appear in cons1 */
   int*                  diffidx1minus0,     /**< array with indices of variables in cons1, that don't appear in cons0 */
   int                   nvarscommon,        /**< number of variables, that appear in both constraints */
   int                   commonidxweight,    /**< variable weight sum of common variables */
   int                   diffidx0minus1weight, /**< variable weight sum of variables in cons0, that don't appear in cons1 */
   int                   diffidx1minus0weight, /**< variable weight sum of variables in cons1, that don't appear in cons0 */
   SCIP_Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   SCIP_Bool*            aggregated          /**< pointer to store whether an aggregation was made */
   )
{
   SCIP_CONSDATA* consdata0;
   SCIP_CONSDATA* consdata1;
   SCIP_Real a;
   SCIP_Real b;
   SCIP_Real aggrcoef;
   SCIP_Real scalarsum;
   SCIP_Real bestscalarsum;
   SCIP_Bool betterscalarsum;
   SCIP_Bool commonvarlindependent;  /* indicates whether coefficient vector of common variables in linearly dependent */
   int varweight;
   int nvars;
   int bestvarweight;
   int bestnvars;
   int bestv;
   int v;
   int i;

   assert(scip != NULL);
   assert(cons0 != NULL);
   assert(cons1 != NULL);
   assert(commonidx0 != NULL);
   assert(commonidx1 != NULL);
   assert(diffidx0minus1 != NULL);
   assert(diffidx1minus0 != NULL);
   assert(nvarscommon >= 1);
   assert(commonidxweight >= nvarscommon);
   assert(nchgcoefs != NULL);
   assert(aggregated != NULL);

   assert(SCIPconsIsActive(cons0));
   assert(SCIPconsIsActive(cons1));

   SCIPdebugMsg(scip, "try aggregation of <%s> and <%s>\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));

   /* cons0 is an (in)equality */
   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   assert(SCIPisLE(scip, consdata0->lhs, consdata0->rhs));
   assert(diffidx0minus1weight >= consdata0->nvars - nvarscommon);

   /* cons1 is an equality */
   consdata1 = SCIPconsGetData(cons1);
   assert(consdata1 != NULL);
   assert(consdata1->nvars >= 1);
   assert(SCIPisEQ(scip, consdata1->lhs, consdata1->rhs));
   assert(diffidx1minus0weight >= consdata1->nvars - nvarscommon);

   *aggregated = FALSE;

   /* search for the best common variable such that
    *   val1[var] * consdata0 - val0[var] * consdata1
    * has least weighted number of variables
    */
   bestvarweight = commonidxweight + diffidx0minus1weight;
   bestnvars = consdata0->nvars;
   bestv = -1;
   bestscalarsum = 0.0;
   commonvarlindependent = TRUE;
   for( v = 0; v < nvarscommon; ++v )
   {
      assert(consdata0->vars[commonidx0[v]] == consdata1->vars[commonidx1[v]]);
      a = consdata1->vals[commonidx1[v]];
      b = -consdata0->vals[commonidx0[v]];

      /* only try aggregation, if coefficients are integral (numerical stability) */
      if( SCIPisIntegral(scip, a) && SCIPisIntegral(scip, b) )
      {
         /* count the number of variables in the potential new constraint  a * consdata0 + b * consdata1 */
         varweight = diffidx0minus1weight + diffidx1minus0weight;
         nvars = consdata0->nvars + consdata1->nvars - 2*nvarscommon;
         scalarsum = REALABS(a) + REALABS(b);
         betterscalarsum = (scalarsum < bestscalarsum);
         for( i = 0; i < nvarscommon
                 && (varweight < bestvarweight || (varweight == bestvarweight && betterscalarsum)); ++i )
         {
            aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
            if( !SCIPisZero(scip, aggrcoef) )
            {
               varweight += getVarWeight(consdata0->vars[commonidx0[i]]);
               nvars++;
            }
         }
         if( varweight < bestvarweight || (varweight == bestvarweight && betterscalarsum) )
         {
            bestv = v;
            bestvarweight = varweight;
            bestnvars = nvars;
            bestscalarsum = scalarsum;
         }
      }

      /* update commonvarlindependent flag, if still TRUE:
       * v's common coefficient in cons1 / v's common coefficient in cons0 should be constant, i.e., equal 0's common coefficient in cons1 / 0's common coefficient in cons0
       */
      if( commonvarlindependent && v > 0 )
         commonvarlindependent = SCIPisEQ(scip,
            consdata1->vals[commonidx1[v]] * consdata0->vals[commonidx0[0]],
            consdata1->vals[commonidx1[0]] * consdata0->vals[commonidx0[v]]);
   }

   /* if better aggregation was found, create new constraint and delete old one */
   if( (bestv != -1 || commonvarlindependent) && SCIPconsGetNUpgradeLocks(cons0) == 0 )
   {
      SCIP_CONS* newcons;
      SCIP_CONSDATA* newconsdata;
      SCIP_VAR** newvars;
      SCIP_Real* newvals;
      SCIP_Real newlhs;
      SCIP_Real newrhs;
      int newnvars;

      if( bestv != -1 )
      {
         /* choose multipliers such that the multiplier for the (in)equality cons0 is positive */
         if( consdata1->vals[commonidx1[bestv]] > 0.0 )
         {
            a = consdata1->vals[commonidx1[bestv]];
            b = -consdata0->vals[commonidx0[bestv]];
         }
         else
         {
            a = -consdata1->vals[commonidx1[bestv]];
            b = consdata0->vals[commonidx0[bestv]];
         }
         assert(SCIPisIntegral(scip, a));
         assert(SCIPisPositive(scip, a));
         assert(SCIPisIntegral(scip, b));
         assert(!SCIPisZero(scip, b));
      }
      else
      {
         assert(commonvarlindependent);
         if( consdata1->vals[commonidx1[0]] > 0.0 )
         {
            a =  consdata1->vals[commonidx1[0]];
            b = -consdata0->vals[commonidx0[0]];
         }
         else
         {
            a = -consdata1->vals[commonidx1[0]];
            b =  consdata0->vals[commonidx0[0]];
         }
         assert(SCIPisPositive(scip, a));
         assert(!SCIPisZero(scip, b));

         /* if a/b is integral, then we can easily choose integer multipliers */
         if( SCIPisIntegral(scip, a/b) )
         {
            if( a/b > 0 )
            {
               a /= b;
               b = 1.0;
            }
            else
            {
               a /= -b;
               b = -1.0;
            }
         }

         /* setup best* variables that were not setup above because we are in the commonvarlindependent case */
         SCIPdebug( bestvarweight = diffidx0minus1weight + diffidx1minus0weight; )
         bestnvars = consdata0->nvars + consdata1->nvars - 2*nvarscommon;
      }

      SCIPdebugMsg(scip, "aggregate linear constraints <%s> := %.15g*<%s> + %.15g*<%s>  ->  nvars: %d -> %d, weight: %d -> %d\n",
         SCIPconsGetName(cons0), a, SCIPconsGetName(cons0), b, SCIPconsGetName(cons1),
         consdata0->nvars, bestnvars, commonidxweight + diffidx0minus1weight, bestvarweight);
      SCIPdebugPrintCons(scip, cons0, NULL);
      SCIPdebugPrintCons(scip, cons1, NULL);

      /* get temporary memory for creating the new linear constraint */
      SCIP_CALL( SCIPallocBufferArray(scip, &newvars, bestnvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newvals, bestnvars) );

      /* calculate the common coefficients, if we have not recognized linear dependency */
      newnvars = 0;
      if( !commonvarlindependent )
      {
         for( i = 0; i < nvarscommon; ++i )
         {
            assert(0 <= commonidx0[i] && commonidx0[i] < consdata0->nvars);
            assert(0 <= commonidx1[i] && commonidx1[i] < consdata1->nvars);

            aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
            if( !SCIPisZero(scip, aggrcoef) )
            {
               assert(newnvars < bestnvars);
               newvars[newnvars] = consdata0->vars[commonidx0[i]];
               newvals[newnvars] = aggrcoef;
               newnvars++;
            }
         }
      }
      else
      {
         /* if we recognized linear dependency of the common coefficients, then the aggregation coefficient should be 0.0 for every common variable */
#ifndef NDEBUG
         for( i = 0; i < nvarscommon; ++i )
         {
            assert(0 <= commonidx0[i] && commonidx0[i] < consdata0->nvars);
            assert(0 <= commonidx1[i] && commonidx1[i] < consdata1->nvars);

            aggrcoef = a * consdata0->vals[commonidx0[i]] + b * consdata1->vals[commonidx1[i]];
            assert(SCIPisZero(scip, aggrcoef));
         }
#endif
      }

      /* calculate the coefficients appearing in cons0 but not in cons1 */
      for( i = 0; i < consdata0->nvars - nvarscommon; ++i )
      {
         assert(0 <= diffidx0minus1[i] && diffidx0minus1[i] < consdata0->nvars);

         aggrcoef = a * consdata0->vals[diffidx0minus1[i]];
         assert(!SCIPisZero(scip, aggrcoef));
         assert(newnvars < bestnvars);
         newvars[newnvars] = consdata0->vars[diffidx0minus1[i]];
         newvals[newnvars] = aggrcoef;
         newnvars++;
      }

      /* calculate the coefficients appearing in cons1 but not in cons0 */
      for( i = 0; i < consdata1->nvars - nvarscommon; ++i )
      {
         assert(0 <= diffidx1minus0[i] && diffidx1minus0[i] < consdata1->nvars);

         aggrcoef = b * consdata1->vals[diffidx1minus0[i]];
         assert(!SCIPisZero(scip, aggrcoef));
         assert(newnvars < bestnvars);
         newvars[newnvars] = consdata1->vars[diffidx1minus0[i]];
         newvals[newnvars] = aggrcoef;
         newnvars++;
      }
      assert(newnvars == bestnvars);

      /* calculate the new left and right hand side of the (in)equality */
      assert(!SCIPisInfinity(scip, -consdata1->lhs));
      assert(!SCIPisInfinity(scip, consdata1->rhs));
      if( SCIPisInfinity(scip, -consdata0->lhs) )
         newlhs = -SCIPinfinity(scip);
      else
         newlhs = a * consdata0->lhs + b * consdata1->lhs;
      if( SCIPisInfinity(scip, consdata0->rhs) )
         newrhs = SCIPinfinity(scip);
      else
         newrhs = a * consdata0->rhs + b * consdata1->rhs;

      /* create the new linear constraint */
      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, SCIPconsGetName(cons0), newnvars, newvars, newvals, newlhs, newrhs,
            SCIPconsIsInitial(cons0), SCIPconsIsSeparated(cons0), SCIPconsIsEnforced(cons0),
            SCIPconsIsChecked(cons0), SCIPconsIsPropagated(cons0),
            SCIPconsIsLocal(cons0), SCIPconsIsModifiable(cons0),
            SCIPconsIsDynamic(cons0), SCIPconsIsRemovable(cons0), SCIPconsIsStickingAtNode(cons0)) );

      newconsdata = SCIPconsGetData(newcons);
      assert(newconsdata != NULL);

      /* copy the upgraded flag from the old cons0 to the new constraint */
      newconsdata->upgraded = consdata0->upgraded;

      /* normalize the new constraint */
      SCIP_CALL( normalizeCons(scip, newcons) );

      /* check, if we really want to use the new constraint instead of the old one:
       * use the new one, if the maximum norm doesn't grow too much
       */
      if( consdataGetMaxAbsval(SCIPconsGetData(newcons)) <= maxaggrnormscale * consdataGetMaxAbsval(consdata0) )
      {
         SCIPdebugMsg(scip, " -> aggregated to <%s>\n", SCIPconsGetName(newcons));
         SCIPdebugPrintCons(scip, newcons, NULL);

         /* update the statistics: we changed all coefficients */
         if( !consdata0->upgraded )
            (*nchgcoefs) += consdata0->nvars + consdata1->nvars - nvarscommon;
         *aggregated = TRUE;

         /* delete the old constraint, and add the new linear constraint to the problem */
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
      }

      /* release the new constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &newvals);
      SCIPfreeBufferArray(scip, &newvars);
   }

   return SCIP_OKAY;
}

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyLinearcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables and the 
 * coefficients are either equal or negated
 */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqLinearcons)
{
   SCIP* scip;
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   SCIP_Real cons1scale;
   SCIP_Real cons2scale;
   int i;

   assert(key1 != NULL);
   assert(key2 != NULL);
   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   assert(consdata1->sorted);
   assert(consdata2->sorted);

   scip = (SCIP*)userptr;
   assert(scip != NULL);

   /* if it is the same constraint we dont need to check anything */
   if( key1 == key2 )
      return TRUE;

   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   /* tests if variables are equal */
   for( i = 0; i < consdata1->nvars; ++i )
   {
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 ||
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 0);
   }

   /* compute scale before comparing coefficients of constraints */
   cons1scale = COPYSIGN(1.0/consdata1->maxabsval, consdata1->vals[0]);
   cons2scale = COPYSIGN(1.0/consdata2->maxabsval, consdata2->vals[0]);

   /* tests if coefficients are equal with the computed scale */
   for( i = 0; i < consdata1->nvars; ++i )
   {
      SCIP_Real val1;
      SCIP_Real val2;

      val1 = consdata1->vals[i] * cons1scale;
      val2 = consdata2->vals[i] * cons2scale;

      if( !SCIPisEQ(scip, val1, val2) )
         return FALSE;
   }

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValLinearcons)
{
   SCIP_CONSDATA* consdata;
   int minidx;
   int mididx;
   int maxidx;
   SCIP_Real scale;
#ifndef NDEBUG
   SCIP* scip;

   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   assert(key != NULL);
   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->nvars > 0);

   assert(consdata->sorted);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   assert(minidx >= 0 && minidx <= maxidx);
   scale = COPYSIGN(1.0/consdata->maxabsval, consdata->vals[0]);

   /* using only the variable indices as hash, since the values are compared by epsilon */
   return SCIPhashFour(consdata->nvars,
                       SCIPcombineTwoInt(minidx, SCIPrealHashCode(consdata->vals[0] * scale)),
                       SCIPcombineTwoInt(mididx, SCIPrealHashCode(consdata->vals[consdata->nvars / 2] * scale)),
                       SCIPcombineTwoInt(maxidx, SCIPrealHashCode(consdata->vals[consdata->nvars - 1] * scale))); /*lint !e571*/
}

/** returns the key for deciding which of two parallel constraints should be kept (smaller key should be kept);
 *  prefers non-upgraded constraints and as second criterion the constraint with the smallest position
 */
static
unsigned int getParallelConsKey(
   SCIP_CONS*            cons                /**< linear constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return (((unsigned int)consdata->upgraded)<<31) + (unsigned int)SCIPconsGetPos(cons); /*lint !e571*/
}

/** updates the hashtable such that out of all constraints in the hashtable that are detected
 *  to be parallel to *querycons, only one is kept in the hashtable and stored into *querycons,
 *  and all others are removed from the hashtable and stored in the given array
 */
static
SCIP_RETCODE retrieveParallelConstraints(
   SCIP_HASHTABLE*       hashtable,          /**< hashtable containing linear constraints */
   SCIP_CONS**           querycons,          /**< pointer to linear constraint used to look for duplicates in the hash table;
                                              *   upon return will contain the constraint that should be kept */
   SCIP_CONS**           parallelconss,      /**< array to return constraints that are parallel to the given;
                                              *   these constraints where removed from the hashtable */
   int*                  nparallelconss      /**< pointer to return number of parallel constraints */
   )
{
   SCIP_CONS* parallelcons;
   unsigned int querykey;

   *nparallelconss = 0;
   querykey = getParallelConsKey(*querycons);

   while( (parallelcons = (SCIP_CONS*)SCIPhashtableRetrieve(hashtable, (void*)(*querycons))) != NULL )
   {
      unsigned int conskey = getParallelConsKey(parallelcons);

      if( conskey < querykey )
      {
         parallelconss[(*nparallelconss)++] = *querycons;
         *querycons = parallelcons;
         querykey = conskey;
      }
      else
      {
         parallelconss[(*nparallelconss)++] = parallelcons;
      }

      /* if the constraint that just came out of the hash table is the one that is kept,
       * we do not need to look into the hashtable again, since the invariant is that
       * in the hashtable only pair-wise non-parallel constraints are contained.
       * For the original querycons, however, multiple constraints that compare equal (=parallel)
       * could be contained due to non-transitivity of the equality comparison.
       * Also we can return immediately, since parallelcons is already contained in the
       * hashtable and we do not need to remove and reinsert it.
       */
      if( *querycons == parallelcons )
         return SCIP_OKAY;

      /* remove parallelcons from the hashtable, since it will be replaced by querycons */
      SCIP_CALL( SCIPhashtableRemove(hashtable, (void*) parallelcons) );
   }

   /* in debug mode we make sure, that the hashtable cannot contain a constraint that
    * comnpares equal to querycons at this point
    */
#ifndef NDEBUG
   SCIP_CALL_ABORT( SCIPhashtableSafeInsert(hashtable, *querycons) );
#else
   SCIP_CALL( SCIPhashtableInsert(hashtable, *querycons) );
#endif

   return SCIP_OKAY;
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint 
 *  accordingly; in contrast to preprocessConstraintPairs(), it uses a hash table 
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  firstchange,        /**< pointer to store first changed constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_HASHTABLE* hashtable;
   SCIP_CONS** parallelconss;
   int nparallelconss;
   int hashtablesize;
   int c;

   assert(scip != NULL);
   assert(blkmem != NULL);
   assert(conss != NULL);
   assert(firstchange != NULL);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = nconss;
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyLinearcons, hashKeyEqLinearcons, hashKeyValLinearcons, (void*) scip) );

   SCIP_CALL( SCIPallocBufferArray(scip, &parallelconss, nconss) );

   /* check all constraints in the given set for redundancy */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONSDATA* consdata0;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      /* check for interuption */
      if( c % 1000 == 0 && SCIPisStopped(scip) )
         break;

      /* sorts the constraint */
      consdata0 = SCIPconsGetData(cons0);
      assert(consdata0 != NULL);
      SCIP_CALL( consdataSort(scip, consdata0) );
      assert(consdata0->sorted);

      /* get constraints from current hash table with same variables as cons0 and with coefficients equal
       * to the ones of cons0 when both are scaled such that maxabsval is 1.0 and the coefficient of the
       * first variable is positive
       * Also inserts cons0 into the hashtable.
       */
      SCIP_CALL( retrieveParallelConstraints(hashtable, &cons0, parallelconss, &nparallelconss) );

      if( nparallelconss != 0 )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

         int i;

         /* cons0 may have been changed in retrieveParallelConstraints() */
         consdata0 = SCIPconsGetData(cons0);

         lhs = consdata0->lhs;
         rhs = consdata0->rhs;

         for( i = 0; i < nparallelconss; ++i )
         {
            SCIP_CONS* consdel;
            SCIP_CONSDATA* consdatadel;
            SCIP_Real scale;

            consdel = parallelconss[i];
            consdatadel = SCIPconsGetData(consdel);

            assert(SCIPconsIsActive(consdel));
            assert(!SCIPconsIsModifiable(consdel));

            /* constraint found: create a new constraint with same coefficients and best left and right hand side;
             * delete old constraints afterwards
             */
            assert(consdatadel != NULL);
            assert(consdata0->nvars >= 1 && consdata0->nvars == consdatadel->nvars);

            assert(consdatadel->sorted);
            assert(consdata0->vars[0] == consdatadel->vars[0]);

            scale = consdata0->vals[0] / consdatadel->vals[0];
            assert(scale != 0.0);

            /* in debug mode, check that all coefficients are equal with respect to epsilon
             * if the constraints are in equilibrium scale
             */
#ifndef NDEBUG
            {
               int k;
               SCIP_Real scale0 = 1.0 / consdata0->maxabsval;
               SCIP_Real scaledel = COPYSIGN(1.0 / consdatadel->maxabsval, scale);

               for( k = 0; k < consdata0->nvars; ++k )
               {
                  assert(SCIPisEQ(scip, scale0 * consdata0->vals[k], scaledel * consdatadel->vals[k]));
               }
            }
#endif

            if( scale > 0.0 )
            {
               /* the coefficients of both constraints are parallel with a positive scale */
               SCIPdebugMsg(scip, "aggregate linear constraints <%s> and <%s> with equal coefficients into single ranged row\n",
                            SCIPconsGetName(cons0), SCIPconsGetName(consdel));
               SCIPdebugPrintCons(scip, cons0, NULL);
               SCIPdebugPrintCons(scip, consdel, NULL);

               if( ! SCIPisInfinity(scip, -consdatadel->lhs) )
                  lhs = MAX(scale * consdatadel->lhs, lhs);

               if( ! SCIPisInfinity(scip, consdatadel->rhs) )
                  rhs = MIN(scale * consdatadel->rhs, rhs);
            }
            else
            {
               /* the coefficients of both rows are negations */
               SCIPdebugMsg(scip, "aggregate linear constraints <%s> and <%s> with negated coefficients into single ranged row\n",
                            SCIPconsGetName(cons0), SCIPconsGetName(consdel));
               SCIPdebugPrintCons(scip, cons0, NULL);
               SCIPdebugPrintCons(scip, consdel, NULL);

               if( ! SCIPisInfinity(scip, consdatadel->rhs) )
                  lhs = MAX(scale * consdatadel->rhs, lhs);

               if( ! SCIPisInfinity(scip, -consdatadel->lhs) )
                  rhs = MIN(scale * consdatadel->lhs, rhs);
            }

            /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons0, consdel) );

            /* delete consdel */
            assert( ! consdata0->upgraded || consdatadel->upgraded );
            SCIP_CALL( SCIPdelCons(scip, consdel) );
            if( !consdatadel->upgraded )
               (*ndelconss)++;
         }

         if( SCIPisFeasLT(scip, rhs, lhs) )
         {
            SCIPdebugMsg(scip, "aggregated linear constraint <%s> is infeasible\n", SCIPconsGetName(cons0));
            *cutoff = TRUE;
            break;
         }

         /* ensure that lhs <= rhs holds without tolerances as we only allow such rows to enter the LP */
         if( lhs > rhs )
         {
            rhs = (lhs + rhs)/2;
            lhs = rhs;
         }

         /* update lhs and rhs of cons0 */
         SCIP_CALL( chgLhs(scip, cons0, lhs) );
         SCIP_CALL( chgRhs(scip, cons0, rhs) );

         /* update the first changed constraint to begin the next aggregation round with */
         if( consdata0->changed && SCIPconsGetPos(cons0) < *firstchange )
            *firstchange = SCIPconsGetPos(cons0);

         assert(SCIPconsIsActive(cons0));
      }
   }
#ifdef  SCIP_MORE_DEBUG
   SCIPinfoMessage(scip, NULL, "linear pairwise comparison hashtable statistics:\n");
   SCIPhashtablePrintStatistics(hashtable, SCIPgetMessagehdlr(scip));
#endif

   SCIPfreeBufferArray(scip, &parallelconss);

   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}

/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   SCIP_Real             maxaggrnormscale,   /**< maximal allowed relative gain in maximum norm for constraint aggregation */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgsides,          /**< pointer to count number of changed left/right hand sides */
   int*                  nchgcoefs           /**< pointer to count number of changed coefficients */
   )
{
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   int* commonidx0;
   int* commonidx1;
   int* diffidx0minus1;
   int* diffidx1minus0;
   uint64_t possignature0;
   uint64_t negsignature0;
   SCIP_Bool cons0changed;
   SCIP_Bool cons0isequality;
   int diffidx1minus0size;
   int c;
   SCIP_Real cons0lhs;
   SCIP_Real cons0rhs;
   SCIP_Bool cons0upgraded;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(cutoff != NULL);
   assert(ndelconss != NULL);
   assert(nchgsides != NULL);
   assert(nchgcoefs != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(cons0 != NULL);
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);
   cons0isequality = SCIPisEQ(scip, consdata0->lhs, consdata0->rhs);

   /* sort the constraint */
   SCIP_CALL( consdataSort(scip, consdata0) );

   /* calculate bit signatures of cons0 for potentially positive and negative coefficients */
   consdataCalcSignatures(consdata0);
   possignature0 = consdata0->possignature;
   negsignature0 = consdata0->negsignature;

   /* get temporary memory for indices of common variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &commonidx0, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &commonidx1, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &diffidx0minus1, consdata0->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &diffidx1minus0, consdata0->nvars) );
   diffidx1minus0size = consdata0->nvars;

   cons0lhs = consdata0->lhs;
   cons0rhs = consdata0->rhs;
   cons0upgraded = consdata0->upgraded;

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   consdata0->changed = FALSE;
   for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff) && conss[chkind] != NULL; ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      uint64_t possignature1;
      uint64_t negsignature1;
      SCIP_Bool cons0dominateslhs;
      SCIP_Bool cons1dominateslhs;
      SCIP_Bool cons0dominatesrhs;
      SCIP_Bool cons1dominatesrhs;
      SCIP_Bool cons1isequality;
      SCIP_Bool coefsequal;
      SCIP_Bool coefsnegated;
      SCIP_Bool tryaggregation;
      int nvarscommon;
      int nvars0minus1;
      int nvars1minus0;
      int commonidxweight;
      int diffidx0minus1weight;
      int diffidx1minus0weight;
      int v0;
      int v1;

      assert(cons0lhs == consdata0->lhs);  /*lint !e777*/
      assert(cons0rhs == consdata0->rhs);  /*lint !e777*/
      assert(cons0upgraded == consdata0->upgraded);

      cons1 = conss[c];

      /* cons1 has become inactive during presolving of constraint pairs */
      if( cons1 == NULL )
         continue;

      assert(SCIPconsIsActive(cons0) && !SCIPconsIsModifiable(cons0));
      assert(SCIPconsIsActive(cons1) && !SCIPconsIsModifiable(cons1));

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* SCIPdebugMsg(scip, "preprocess linear constraint pair <%s>[chgd:%d, upgd:%d] and <%s>[chgd:%d, upgd:%d]\n",
         SCIPconsGetName(cons0), cons0changed, cons0upgraded,
         SCIPconsGetName(cons1), consdata1->changed, consdata1->upgraded); */

      /* if both constraints didn't change since last pair processing, we can ignore the pair */
      if( !cons0changed && !consdata1->changed )
         continue;

      /* if both constraints are already upgraded, skip the pair; 
       * because changes on these constraints cannot be applied to the instance anymore */
      if( cons0upgraded && consdata1->upgraded )
         continue;

      assert(consdata1->nvars >= 1);

      /* sort the constraint */
      SCIP_CALL( consdataSort(scip, consdata1) );

      /* calculate bit signatures of cons1 for potentially positive and negative coefficients */
      consdataCalcSignatures(consdata1);
      possignature1 = consdata1->possignature;
      negsignature1 = consdata1->negsignature;

      /* the signatures give a quick test to check for domination and equality of coefficients */
      coefsequal = (possignature0 == possignature1) && (negsignature0 == negsignature1);
      coefsnegated = (possignature0 == negsignature1) && (negsignature0 == possignature1);
      cons0dominateslhs = SCIPisGE(scip, cons0lhs, consdata1->lhs)
         && ((possignature0 | possignature1) == possignature1)  /* possignature0 <= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature0); /* negsignature0 >= negsignature1 (as bit vector) */
      cons1dominateslhs = SCIPisGE(scip, consdata1->lhs, cons0lhs)
         && ((possignature0 | possignature1) == possignature0)  /* possignature0 >= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature1); /* negsignature0 <= negsignature1 (as bit vector) */
      cons0dominatesrhs = SCIPisLE(scip, cons0rhs, consdata1->rhs)
         && ((possignature0 | possignature1) == possignature0)  /* possignature0 >= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature1); /* negsignature0 <= negsignature1 (as bit vector) */
      cons1dominatesrhs = SCIPisLE(scip, consdata1->rhs, cons0rhs)
         && ((possignature0 | possignature1) == possignature1)  /* possignature0 <= possignature1 (as bit vector) */
         && ((negsignature0 | negsignature1) == negsignature0); /* negsignature0 >= negsignature1 (as bit vector) */
      cons1isequality = SCIPisEQ(scip, consdata1->lhs, consdata1->rhs);
      tryaggregation = (cons0isequality || cons1isequality) && (maxaggrnormscale > 0.0);
      if( !cons0dominateslhs && !cons1dominateslhs && !cons0dominatesrhs && !cons1dominatesrhs
         && !coefsequal && !coefsnegated && !tryaggregation )
         continue;

      /* make sure, we have enough memory for the index set of V_1 \ V_0 */
      if( tryaggregation && consdata1->nvars > diffidx1minus0size )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &diffidx1minus0, consdata1->nvars) );
         diffidx1minus0size = consdata1->nvars;
      }

      /* check consdata0 against consdata1:
       * - if lhs0 >= lhs1 and for each variable v and each solution value x_v val0[v]*x_v <= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. left hand side
       * - if rhs0 <= rhs1 and for each variable v and each solution value x_v val0[v]*x_v >= val1[v]*x_v,
       *   consdata0 dominates consdata1 w.r.t. right hand side
       * - if val0[v] == -val1[v] for all variables v, the two inequalities can be replaced by a single
       *   ranged row (or equality)
       * - if at least one constraint is an equality, count the weighted number of common variables W_c
       *   and the weighted number of variable in the difference sets W_0 = w(V_0 \ V_1), W_1 = w(V_1 \ V_0),
       *   where the weight of each variable depends on its type, such that aggregations in order to remove the
       *   number of continuous and integer variables are preferred:
       *   - if W_c > W_1, try to aggregate  consdata0 := a * consdata0 + b * consdata1  in order to decrease the
       *     variable weight in consdata0, where a = +/- val1[v] and b = -/+ val0[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons0
       *   - if W_c > W_0, try to aggregate  consdata1 := a * consdata1 + b * consdata0  in order to decrease the
       *     variable weight in consdata1, where a = +/- val0[v] and b = -/+ val1[v] for common v which leads to
       *     the smallest weight; for numerical stability, we will only accept integral a and b; the sign of a has
       *     to be positive to not switch the sense of the (in)equality cons1
       */

      /* check consdata0 against consdata1 for redundancy, or ranged row accumulation */
      nvarscommon = 0;
      commonidxweight = 0;
      nvars0minus1 = 0;
      diffidx0minus1weight = 0;
      nvars1minus0 = 0;
      diffidx1minus0weight = 0;
      v0 = 0;
      v1 = 0;
      while( (v0 < consdata0->nvars || v1 < consdata1->nvars)
         && (cons0dominateslhs || cons1dominateslhs || cons0dominatesrhs || cons1dominatesrhs
            || coefsequal || coefsnegated || tryaggregation) )
      {
         SCIP_VAR* var;
         SCIP_Real val0;
         SCIP_Real val1;
         int varcmp;

         /* test, if variable appears in only one or in both constraints */
         if( v0 < consdata0->nvars && v1 < consdata1->nvars )
            varcmp = SCIPvarCompare(consdata0->vars[v0], consdata1->vars[v1]);
         else if( v0 < consdata0->nvars )
            varcmp = -1;
         else
            varcmp = +1;

         switch( varcmp )
         {
         case -1:
            /* variable doesn't appear in consdata1 */
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = 0.0;
            if( tryaggregation )
            {
               diffidx0minus1[nvars0minus1] = v0;
               nvars0minus1++;
               diffidx0minus1weight += getVarWeight(var);
            }
            v0++;
            coefsequal = FALSE;
            coefsnegated = FALSE;
            break;

         case +1:
            /* variable doesn't appear in consdata0 */
            var = consdata1->vars[v1];
            val0 = 0.0;
            val1 = consdata1->vals[v1];
            if( tryaggregation )
            {
               diffidx1minus0[nvars1minus0] = v1;
               nvars1minus0++;
               diffidx1minus0weight += getVarWeight(var);
            }
            v1++;
            coefsequal = FALSE;
            coefsnegated = FALSE;
            break;

         case 0:
            /* variable appears in both constraints */
            assert(consdata0->vars[v0] == consdata1->vars[v1]);
            var = consdata0->vars[v0];
            val0 = consdata0->vals[v0];
            val1 = consdata1->vals[v1];
            if( tryaggregation )
            {
               commonidx0[nvarscommon] = v0;
               commonidx1[nvarscommon] = v1;
               nvarscommon++;
               commonidxweight += getVarWeight(var);
            }
            v0++;
            v1++;
            coefsequal = coefsequal && (SCIPisEQ(scip, val0, val1));
            coefsnegated = coefsnegated && (SCIPisEQ(scip, val0, -val1));
            break;

         default:
            SCIPerrorMessage("invalid comparison result\n");
            SCIPABORT();
            var = NULL;
            val0 = 0.0;
            val1 = 0.0;
         }
         assert(var != NULL);

         /* update domination criteria w.r.t. the coefficient and the variable's bounds */
         if( SCIPisGT(scip, val0, val1) )
         {
            if( SCIPisNegative(scip, SCIPvarGetLbGlobal(var)) )
            {
               cons0dominatesrhs = FALSE;
               cons1dominateslhs = FALSE;
            }
            if( SCIPisPositive(scip, SCIPvarGetUbGlobal(var)) )
            {
               cons0dominateslhs = FALSE;
               cons1dominatesrhs = FALSE;
            }
         }
         else if( SCIPisLT(scip, val0, val1) )
         {
            if( SCIPisNegative(scip, SCIPvarGetLbGlobal(var)) )
            {
               cons0dominateslhs = FALSE;
               cons1dominatesrhs = FALSE;
            }
            if( SCIPisPositive(scip, SCIPvarGetUbGlobal(var)) )
            {
               cons0dominatesrhs = FALSE;
               cons1dominateslhs = FALSE;
            }
         }
      }

      /* check for disaggregated ranged rows */
      if( coefsequal || coefsnegated )
      {
         SCIP_CONS* consstay;
         SCIP_CONS* consdel;
#ifndef NDEBUG
         SCIP_CONSDATA* consdatastay;
#endif
         SCIP_CONSDATA* consdatadel;
         SCIP_Real lhs;
         SCIP_Real rhs;
         int consinddel;

         /* the coefficients in both rows are either equal or negated: create a new constraint with same coefficients and
          * best left and right hand sides; delete the old constraints afterwards
          */
         SCIPdebugMsg(scip, "aggregate linear constraints <%s> and <%s> with %s coefficients into single ranged row\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1), coefsequal ? "equal" : "negated");
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         if( coefsequal )
         {
            /* the coefficients of both rows are equal */
            lhs = MAX(consdata0->lhs, consdata1->lhs);
            rhs = MIN(consdata0->rhs, consdata1->rhs);
         }
         else
         {
            /* the coefficients of both rows are negations */
            lhs = MAX(consdata0->lhs, -consdata1->rhs);
            rhs = MIN(consdata0->rhs, -consdata1->lhs);
         }
         if( SCIPisFeasLT(scip, rhs, lhs) )
         {
            SCIPdebugMsg(scip, "aggregated linear constraint <%s> is infeasible\n", SCIPconsGetName(cons0));
            *cutoff = TRUE;
            break;
         }

         /* check which constraint has to stay; 
          * changes applied to an upgraded constraint will not be considered in the instance */
         if( consdata0->upgraded )
         {
            assert(!consdata1->upgraded);
            consstay = cons1;
#ifndef NDEBUG
            consdatastay = consdata1;
#endif

            consdel = cons0;
            consdatadel = consdata0;
            consinddel = chkind;
         }
         else
         {
            consstay = cons0;
#ifndef NDEBUG
            consdatastay = consdata0;
#endif

            consdel = cons1;
            consdatadel = consdata1;
            consinddel = c;
         }

         /* update the sides of consstay */
         SCIP_CALL( chgLhs(scip, consstay, lhs) );
         SCIP_CALL( chgRhs(scip, consstay, rhs) );
         if( !consdata0->upgraded )
         {
            assert(consstay == cons0);
            cons0lhs = consdata0->lhs;
            cons0rhs = consdata0->rhs;
         }

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( SCIPupdateConsFlags(scip, consstay, consdel) );

         assert( !consdatastay->upgraded );
         /* delete consdel */
         SCIP_CALL( SCIPdelCons(scip, consdel) );
         conss[consinddel] = NULL;
         if( !consdatadel->upgraded )
            (*ndelconss)++;
         continue;
      }

      /* check for domination: remove dominated sides, but don't touch equalities as long as they are not totally
       * redundant
       */
      if( cons1dominateslhs && (!cons0isequality || cons1dominatesrhs || SCIPisInfinity(scip, consdata0->rhs) ) )
      {
         /* left hand side is dominated by consdata1: delete left hand side of consdata0 */
         SCIPdebugMsg(scip, "left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata1->lhs, consdata0->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant left hand side */
         if( !SCIPisInfinity(scip, -consdata0->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons0, -SCIPinfinity(scip)) );
            cons0lhs = consdata0->lhs;
            cons0isequality = FALSE;
            if( !consdata0->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( SCIPupdateConsFlags(scip, cons1, cons0) );

               (*nchgsides)++;
            }
         }
      }
      else if( cons0dominateslhs && (!cons1isequality || cons0dominatesrhs || SCIPisInfinity(scip, consdata1->rhs)) )
      {
         /* left hand side is dominated by consdata0: delete left hand side of consdata1 */
         SCIPdebugMsg(scip, "left hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebugPrintCons(scip, cons1, NULL);
         SCIPdebugPrintCons(scip, cons0, NULL);

         /* check for infeasibility */
         if( SCIPisFeasGT(scip, consdata0->lhs, consdata1->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant left hand side */
         if( !SCIPisInfinity(scip, -consdata1->lhs) )
         {
            SCIP_CALL( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
            cons1isequality = FALSE;
            if( !consdata1->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );

               (*nchgsides)++;
            }
         }
      }
      if( cons1dominatesrhs && (!cons0isequality || cons1dominateslhs || SCIPisInfinity(scip, -consdata0->lhs)) )
      {
         /* right hand side is dominated by consdata1: delete right hand side of consdata0 */
         SCIPdebugMsg(scip, "right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1));
         SCIPdebugPrintCons(scip, cons0, NULL);
         SCIPdebugPrintCons(scip, cons1, NULL);

         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata1->rhs, consdata0->lhs) )
         {
            SCIPdebugMsg(scip, "linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant right hand side */
         if( !SCIPisInfinity(scip, consdata0->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons0, SCIPinfinity(scip)) );
            cons0rhs = consdata0->rhs;
            cons0isequality = FALSE;
            if( !consdata0->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( SCIPupdateConsFlags(scip, cons1, cons0) );

               (*nchgsides)++;
            }
         }
      }
      else if( cons0dominatesrhs && (!cons1isequality || cons0dominateslhs || SCIPisInfinity(scip, -consdata1->lhs)) )
      {
         /* right hand side is dominated by consdata0: delete right hand side of consdata1 */
         SCIPdebugMsg(scip, "right hand side of linear constraint <%s> is dominated by <%s>:\n",
            SCIPconsGetName(cons1), SCIPconsGetName(cons0));
         SCIPdebugPrintCons(scip, cons1, NULL);
         SCIPdebugPrintCons(scip, cons0, NULL);

         /* check for infeasibility */
         if( SCIPisFeasLT(scip, consdata0->rhs, consdata1->lhs) )
         {
            SCIPdebugMsg(scip, "linear constraints <%s> and <%s> are infeasible\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* remove redundant right hand side */
         if( !SCIPisInfinity(scip, consdata1->rhs) )
         {
            SCIP_CALL( chgRhs(scip, cons1, SCIPinfinity(scip)) );
            cons1isequality = FALSE;
            if( !consdata1->upgraded )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );

               (*nchgsides)++;
            }
         }
      }

      /* check for now redundant constraints */
      if( SCIPisInfinity(scip, -consdata0->lhs) && SCIPisInfinity(scip, consdata0->rhs) )
      {
         /* consdata0 became redundant */
         SCIPdebugMsg(scip, "linear constraint <%s> is redundant\n", SCIPconsGetName(cons0));
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         conss[chkind] = NULL;
         if( !consdata0->upgraded )
         {
            /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons1, cons0) );

            (*ndelconss)++;
         }
         continue;
      }
      if( SCIPisInfinity(scip, -consdata1->lhs) && SCIPisInfinity(scip, consdata1->rhs) )
      {
         /* consdata1 became redundant */
         SCIPdebugMsg(scip, "linear constraint <%s> is redundant\n", SCIPconsGetName(cons1));
         SCIP_CALL( SCIPdelCons(scip, cons1) );
         conss[c] = NULL;
         if( !consdata1->upgraded )
         {
            /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
            SCIP_CALL( SCIPupdateConsFlags(scip, cons0, cons1) );

            (*ndelconss)++;
         }
         continue;
      }

      /* check, if we want to aggregate an (in)equality with an equality:
       *   consdata0 := a * consdata0 + b * consdata1  or  consdata1 := a * consdata1 + b * consdata0
       */
      if( tryaggregation )
      {
         SCIP_Bool aggregated;

         assert(consdata0->nvars == nvarscommon + nvars0minus1);
         assert(consdata1->nvars == nvarscommon + nvars1minus0);

         aggregated = FALSE;
         if( cons1isequality && !consdata0->upgraded && commonidxweight > diffidx1minus0weight )
         {
            /* W_c > W_1: try to aggregate  consdata0 := a * consdata0 + b * consdata1 */
            SCIP_CALL( aggregateConstraints(scip, cons0, cons1, commonidx0, commonidx1, diffidx0minus1, diffidx1minus0,
                  nvarscommon, commonidxweight, diffidx0minus1weight, diffidx1minus0weight, maxaggrnormscale,
                  nchgcoefs, &aggregated) );

            /* update array of active constraints */
            if( aggregated )
            {
               assert(!SCIPconsIsActive(cons0));
               assert(SCIPconsIsActive(cons1));
               conss[chkind] = NULL;
            }
         }
         if( !aggregated && cons0isequality && !consdata1->upgraded && commonidxweight > diffidx0minus1weight )
         {
            /* W_c > W_0: try to aggregate  consdata1 := a * consdata1 + b * consdata0 */
            SCIP_CALL( aggregateConstraints(scip, cons1, cons0, commonidx1, commonidx0, diffidx1minus0, diffidx0minus1,
                  nvarscommon, commonidxweight, diffidx1minus0weight, diffidx0minus1weight, maxaggrnormscale,
                  nchgcoefs, &aggregated) );

            /* update array of active constraints */
            if( aggregated )
            {
               assert(!SCIPconsIsActive(cons1));
               assert(SCIPconsIsActive(cons0));
               conss[c] = NULL;
            }
         }
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &diffidx1minus0);
   SCIPfreeBufferArray(scip, &diffidx0minus1);
   SCIPfreeBufferArray(scip, &commonidx1);
   SCIPfreeBufferArray(scip, &commonidx0);

   return SCIP_OKAY;
}

/** do stuffing presolving on a single constraint */
static
SCIP_RETCODE presolStuffing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Bool             singletonstuffing,  /**< should stuffing of singleton continuous variables be performed? */
   SCIP_Bool             singlevarstuffing,  /**< should single variable stuffing be performed, which tries to fulfill
                                              *   constraints using the cheapest variable? */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to count the total number of fixed variables */
   int*                  nchgbds             /**< pointer to count the total number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* ratios;
   int* varpos;
   SCIP_Bool* swapped;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_VAR* var;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real maxcondactivity;
   SCIP_Real mincondactivity;
   SCIP_Real rhs;
   SCIP_Real val;
   SCIP_Real obj;
   SCIP_Real factor;
   SCIP_Bool minactisrelax;
   SCIP_Bool maxactisrelax;
   SCIP_Bool tryfixing;
   int nsingletons;
   int idx;
   int v;
   int nvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nfixedvars != NULL);

   consdata = SCIPconsGetData(cons);

   /* we only want to run for inequalities */
   if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->lhs) )
      return SCIP_OKAY;

   if( singlevarstuffing )
   {
      consdataGetActivityBounds(scip, consdata, FALSE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);
   }
   else
   {
      minactivity = SCIP_INVALID;
      maxactivity = SCIP_INVALID;
   }

   /* we want to have a <= constraint, if the rhs is infinite, we implicitly multiply the constraint by -1,
    * the new maxactivity is minus the old minactivity then
    */
   if( SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIP_Real tmp;
      rhs = -consdata->lhs;
      factor = -1.0;
      tmp = maxactivity;
      maxactivity = -minactivity;
      minactivity = -tmp;
   }
   else
   {
      assert(SCIPisInfinity(scip, -consdata->lhs));
      rhs = consdata->rhs;
      factor = 1.0;
   }

   nvars = consdata->nvars;
   vars = consdata->vars;
   vals = consdata->vals;

   /* check for continuous singletons */
   if( singletonstuffing )
   {
      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];

         if( (SCIPvarGetNLocksUp(var) + SCIPvarGetNLocksDown(var)) == 1 &&
            SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            break;
      }
   }
   else
      /* we don't want to go into the next block */
      v = nvars;

   /* a singleton was found -> perform singleton variable stuffing */
   if( v < nvars )
   {
      assert(singletonstuffing);

      SCIP_CALL( SCIPallocBufferArray(scip, &varpos, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ratios, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &swapped, nvars) );

      tryfixing = TRUE;
      nsingletons = 0;
      mincondactivity = 0.0;
      maxcondactivity = 0.0;

      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
         obj = SCIPvarGetObj(var);
         val = factor * vals[v];

         assert(!SCIPisZero(scip, val));

         /* the variable is a singleton and continuous */
         if( (SCIPvarGetNLocksUp(var) + SCIPvarGetNLocksDown(var)) == 1 &&
            SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            if( SCIPisNegative(scip, obj) && val > 0 )
            {
               /* case 1: obj < 0 and coef > 0 */
               if( SCIPisInfinity(scip, -lb) )
               {
                  tryfixing = FALSE;
                  break;
               }

               maxcondactivity += val * lb;
               mincondactivity += val * lb;
               swapped[v] = FALSE;
               ratios[nsingletons] = obj / val;
               varpos[nsingletons] = v;
               nsingletons++;
            }
            else if( SCIPisPositive(scip, obj) && val < 0 )
            {
               /* case 2: obj > 0 and coef < 0 */
               if( SCIPisInfinity(scip, ub) )
               {
                  tryfixing = FALSE;
                  break;
               }
               /* multiply column by (-1) to become case 1.
                * now bounds are swapped: ub := -lb, lb := -ub
                */

               maxcondactivity += val * ub;
               mincondactivity += val * ub;
               swapped[v] = TRUE;
               ratios[nsingletons] = obj / val;
               varpos[nsingletons] = v;
               nsingletons++;
            }
            else if( val > 0 )
            {
               /* case 3: obj >= 0 and coef >= 0 is handled by duality fixing.
                *  we only consider the lower bound for the constants
                */
               assert(!SCIPisNegative(scip, obj));

               if( SCIPisInfinity(scip, -lb) )
               {
                  /* maybe unbounded */
                  tryfixing = FALSE;
                  break;
               }

               maxcondactivity += val * lb;
               mincondactivity += val * lb;
            }
            else
            {
               /* case 4: obj <= 0 and coef <= 0 is also handled by duality fixing.
                * we only consider the upper bound for the constants
                */
               assert(!SCIPisPositive(scip, obj));
               assert(val < 0);

               if( SCIPisInfinity(scip, ub) )
               {
                  /* maybe unbounded */
                  tryfixing = FALSE;
                  break;
               }

               maxcondactivity += val * ub;
               mincondactivity += val * ub;
            }
         }
         else
         {
            /* consider contribution of discrete variables, non-singleton
             * continuous variables and variables with more than one lock
             */
            if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
            {
               tryfixing = FALSE;
               break;
            }

            if( val > 0 )
            {
               maxcondactivity += val * ub;
               mincondactivity += val * lb;
            }
            else
            {
               maxcondactivity += val * lb;
               mincondactivity += val * ub;
            }
         }
      }
      if( tryfixing && nsingletons > 0 && (SCIPisGT(scip, rhs, maxcondactivity) || SCIPisLE(scip, rhs, mincondactivity)) )
      {
         SCIP_Real delta;
         SCIP_Bool tightened;
#ifdef SCIP_DEBUG
         int oldnfixedvars = *nfixedvars;
         int oldnchgbds = *nchgbds;
#endif

         SCIPsortRealInt(ratios, varpos, nsingletons);

         /* verify which singleton continuous variables can be fixed */
         for( v = 0; v < nsingletons; ++v )
         {
            idx = varpos[v];
            var = vars[idx];
            val = factor * vals[idx];
            lb = SCIPvarGetLbGlobal(var);
            ub = SCIPvarGetUbGlobal(var);

            assert(val > 0 || SCIPisPositive(scip, SCIPvarGetObj(var)));
            assert((val < 0) == swapped[idx]);
            val = REALABS(val);

            /* stop fixing if variable bounds are not finite */
            if( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) )
               break;

            assert((SCIPvarGetNLocksUp(var) + SCIPvarGetNLocksDown(var)) == 1);
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

            /* calculate the change in the row activities if this variable changes
             * its value from its worst to its best bound
             */
            if( swapped[idx] )
               delta = -(lb - ub) * val;
            else
               delta =  (ub - lb) * val;

            assert(!SCIPisNegative(scip, delta));

            if( SCIPisLE(scip, delta, rhs - maxcondactivity) )
            {
               if( swapped[idx] )
               {
                  SCIPdebugMsg(scip, "fix <%s> to its lower bound %g\n", SCIPvarGetName(var), lb);
                  SCIP_CALL( SCIPfixVar(scip, var, lb, cutoff, &tightened) );

               }
               else
               {
                  SCIPdebugMsg(scip, "fix <%s> to its upper bound %g\n", SCIPvarGetName(var), ub);
                  SCIP_CALL( SCIPfixVar(scip, var, ub, cutoff, &tightened) );
               }

               if( *cutoff )
                  break;
               if( tightened )
               {
                  (*nfixedvars)++;
               }
            }
            /* @note: we could in theory tighten the bound of the first singleton variable which does not fall into the above case,
             *        since it cannot be fully fixed. However, this is not needed and should be done by activity-based bound tightening
             *        anyway after all other continuous singleton columns were fixed; doing it here may introduce numerical
             *        troubles in case of large bounds.
             */
            else if( SCIPisLE(scip, rhs, mincondactivity) )
            {
               if( swapped[idx] )
               {
                  SCIPdebugMsg(scip, "fix <%s> to its upper bound %g\n", SCIPvarGetName(var), ub);
                  SCIP_CALL( SCIPfixVar(scip, var, ub, cutoff, &tightened) );
               }
               else
               {
                  SCIPdebugMsg(scip, "fix <%s> to its lower bound %g\n", SCIPvarGetName(var), lb);
                  SCIP_CALL( SCIPfixVar(scip, var, lb, cutoff, &tightened) );
               }

               if( *cutoff )
                  break;
               if( tightened )
               {
                  (*nfixedvars)++;
               }
            }

            maxcondactivity += delta;
            mincondactivity += delta;
         }

#ifdef SCIP_DEBUG
         if( *nfixedvars - oldnfixedvars > 0 || *nchgbds - oldnchgbds > 0 )
         {
            SCIPdebugMsg(scip, "### stuffing fixed %d variables and changed %d bounds\n", *nfixedvars - oldnfixedvars, *nchgbds - oldnchgbds);
         }
#endif
      }

      SCIPfreeBufferArray(scip, &swapped);
      SCIPfreeBufferArray(scip, &ratios);
      SCIPfreeBufferArray(scip, &varpos);
   }

   /* perform single-variable stuffing:
    * for a linear inequality
    *  a_1 x_1 + a_2 x_2 + ... + a_n x_n <= b
    * with a_i > 0 and objective coefficients c_i < 0,
    * setting all variables to their upper bound (giving us the maximal activity of the constraint) is worst w.r.t.
    * feasibility of the constraint. On the other hand, this gives the best objective function contribution of the
    * variables contained in the constraint. The maximal activity should be larger than the rhs, otherwise the constraint
    * is redundant.
    * Now we are searching for a variable x_k with maximal ratio c_k / a_k (note that all these ratios are negative), so
    * that by reducing the value of this variable we reduce the activity of the constraint while having the smallest
    * objective deterioration per activity unit. If x_k has no downlocks, is continuous, and can be reduced enough to
    * render the constraint feasible, and ALL other variables have only the one uplock installed by the current constraint,
    * we can reduce the upper bound of x_k such that the maxactivity equals the rhs and fix all other variables to their
    * upper bound.
    * Note that the others variables may have downlocks from other constraints, which we do not need to care
    * about since we are setting them to the highest possible value. Also, they may be integer or binary, because the
    * computed ratio is still a lower bound on the change in the objective caused by reducing those variable to reach
    * constraint feasibility. On the other hand, uplocks on x_k from other constraint do no interfer with the method.
    * With a slight adjustment, the procedure even works for integral x_k. If (maxactivity - rhs)/val is integral,
    * the variable gets an integral value in order to fulfill the constraint tightly, and we can just apply the procedure.
    * If (maxactivity - rhs)/val is fractional, we need to check, if overfulfilling the constraint by setting x_k to
    * ceil((maxactivity - rhs)/val) is still better than setting x_k to ceil((maxactivity - rhs)/val) - 1 and
    * filling the remaining gap in the constraint with the next-best variable. For this, we check that
    * c_k * ceil((maxactivity - rhs)/val) is still better than
    * c_k * floor((maxactivity - rhs)/val) + c_j * ((maxactivity - rhs) - (floor((maxactivity - rhs)/val) * val))/a_j.
    * In this case, the upper bound of x_k is decreased to ub_k - ceil(maxactivity - rhs).
    * If there are variables with a_i < 0 and c_i > 0, they are negated to obtain the above form, variables with same
    * sign of coefficients in constraint and objective prevent the use of this method.
    */
   if( singlevarstuffing && !SCIPisInfinity(scip, -minactivity) )
   {
      SCIP_Real bestratio = -SCIPinfinity(scip);
      SCIP_Real secondbestratio = -SCIPinfinity(scip);
      SCIP_Real ratio;
      int bestindex = -1;
      int bestuplocks = 0;
      int bestdownlocks = 1;
      int downlocks;
      int uplocks;
      int oldnfixedvars;
      int oldnchgbds;

      SCIPdebug( oldnfixedvars = *nfixedvars; )
      SCIPdebug( oldnchgbds = *nchgbds; )
      /* loop over all variables to identify the best and second-best ratio */
      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         obj = SCIPvarGetObj(var);
         val = factor * vals[v];

         assert(!SCIPisZero(scip, val));

         ratio = obj / val;

         /* if both objective and constraint push the variable to the same direction, we can do nothing here */
         if( !SCIPisNegative(scip, ratio) )
         {
            bestindex = -1;
            break;
         }

         if( val > 0 )
         {
            downlocks = SCIPvarGetNLocksDown(var);
            uplocks = SCIPvarGetNLocksUp(var);
         }
         else
         {
            downlocks = SCIPvarGetNLocksUp(var);
            uplocks = SCIPvarGetNLocksDown(var);
         }

         /* better ratio, update best candidate
          * @todo use some tolerance
          * @todo check size of domain and updated ratio for integer variables already?
          */
         if( ratio > bestratio || ((ratio == bestratio) && downlocks == 0 && (bestdownlocks > 0 /*lint !e777*/
                  || (SCIPvarGetType(vars[bestindex]) != SCIP_VARTYPE_CONTINUOUS
                     && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS))) )
         {
            /* best index becomes second-best*/
            if( bestindex != -1 )
            {
               /* second-best index must not have more than 1 uplock */
               if( bestuplocks > 1 )
               {
                  bestindex = -1;
                  break;
               }
               else
               {
                  secondbestratio = bestratio;
               }
            }
            bestdownlocks = downlocks;
            bestuplocks = uplocks;
            bestratio = ratio;
            bestindex = v;

            /* if this variable is the best in the end, we cannot do reductions since it has a downlocks,
             * if it is not the best, it has too many uplocks -> not applicable
             */
            if( bestdownlocks > 0 && bestuplocks > 1 )
            {
               bestindex = -1;
               break;
            }
         }
         else
         {
            /* non-best index must not have more than 1 uplock */
            if( uplocks > 1 )
            {
               bestindex = -1;
               break;
            }
            /* update second-best ratio */
            if( ratio > secondbestratio )
            {
               secondbestratio = ratio;
            }
         }
      }

      /* check if we can apply single variable stuffing */
      if( bestindex != -1 && bestdownlocks == 0 )
      {
         SCIP_Bool tightened = FALSE;
         SCIP_Real bounddelta;

         var = vars[bestindex];
         obj = SCIPvarGetObj(var);
         val = factor * vals[bestindex];
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
         tryfixing = TRUE;

         if( val < 0 )
         {
            assert(!SCIPisNegative(scip, obj));

            /* the best variable is integer, and we need to overfulfill the constraint when using just the variable */
            if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPisIntegral(scip, (maxactivity - rhs)/-val) )
            {
               SCIP_Real bestvarfloor = SCIPfloor(scip, (maxactivity - rhs)/-val);
               SCIP_Real activitydelta = (maxactivity - rhs) - (bestvarfloor * -val);
               assert(SCIPisPositive(scip, activitydelta));

               tryfixing = SCIPisLE(scip, obj, -activitydelta * secondbestratio);

               bounddelta = SCIPceil(scip, (maxactivity - rhs)/-val);
               assert(SCIPisPositive(scip, bounddelta));
            }
            else
               bounddelta = (maxactivity - rhs)/-val;

            tryfixing = tryfixing && SCIPisLE(scip, bounddelta, ub - lb);

            if( tryfixing )
            {
               assert(SCIPvarGetNLocksUp(var) == 0);

               if( SCIPisEQ(scip, lb + bounddelta, ub) )
               {
                  SCIPdebugMsg(scip, "fix var <%s> to %g\n", SCIPvarGetName(var), lb + bounddelta);
                  SCIP_CALL( SCIPfixVar(scip, var, lb + bounddelta, cutoff, &tightened) );
               }
               else
               {
                  SCIPdebugMsg(scip, "tighten the lower bound of <%s> from %g to %g (ub=%g)\n", SCIPvarGetName(var), lb, lb + bounddelta, ub);
                  SCIP_CALL( SCIPtightenVarLb(scip, var, lb + bounddelta, FALSE, cutoff, &tightened) );
               }
            }
         }
         else
         {
            assert(!SCIPisPositive(scip, obj));

            /* the best variable is integer, and we need to overfulfill the constraint when using just the variable */
            if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS  && !SCIPisIntegral(scip, (maxactivity - rhs)/val))
            {
               SCIP_Real bestvarfloor = SCIPfloor(scip, (maxactivity - rhs)/val);
               SCIP_Real activitydelta = (maxactivity - rhs) - (bestvarfloor * val);
               assert(SCIPisPositive(scip, activitydelta));

               tryfixing = SCIPisLE(scip, -obj, activitydelta * secondbestratio);

               bounddelta = SCIPceil(scip, (maxactivity - rhs)/val);
               assert(SCIPisPositive(scip, bounddelta));
            }
            else
               bounddelta = (maxactivity - rhs)/val;

            tryfixing = tryfixing && SCIPisLE(scip, bounddelta, ub - lb);

            if( tryfixing )
            {
               assert(SCIPvarGetNLocksDown(var) == 0);

               if( SCIPisEQ(scip, ub - bounddelta, lb) )
               {
                  SCIPdebugMsg(scip, "fix var <%s> to %g\n", SCIPvarGetName(var), ub - bounddelta);
                  SCIP_CALL( SCIPfixVar(scip, var, ub - bounddelta, cutoff, &tightened) );
               }
               else
               {
                  SCIPdebugMsg(scip, "tighten the upper bound of <%s> from %g to %g (lb=%g)\n", SCIPvarGetName(var), ub, ub - bounddelta, lb);
                  SCIP_CALL( SCIPtightenVarUb(scip, var, ub - bounddelta, FALSE, cutoff, &tightened) );
               }
            }
         }

         if( *cutoff )
            return SCIP_OKAY;
         if( tightened )
         {
            if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
               ++(*nfixedvars);
            else
               ++(*nchgbds);

            SCIPdebugMsg(scip, "cons <%s>: %g <=\n", SCIPconsGetName(cons), factor > 0 ? consdata->lhs : -consdata->rhs);
            for( v = 0; v < nvars; ++v )
            {
               SCIPdebugMsg(scip, "%+g <%s>([%g,%g],%g,[%d,%d],%s)\n", factor * vals[v], SCIPvarGetName(vars[v]), SCIPvarGetLbGlobal(vars[v]),
                  SCIPvarGetUbGlobal(vars[v]), SCIPvarGetObj(vars[v]), SCIPvarGetNLocksDown(vars[v]), SCIPvarGetNLocksUp(vars[v]),
                  SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS ? "C" : "I");
            }
            SCIPdebugMsg(scip, "<= %g\n", factor > 0 ? consdata->rhs : -consdata->lhs);

            for( v = 0; v < nvars; ++v )
            {
               if( v == bestindex )
                  continue;

               if( factor * vals[v] < 0 )
               {
                  assert(SCIPvarGetNLocksDown(vars[v]) == 1);
                  SCIPdebugMsg(scip, "fix <%s> to its lower bound (%g)\n", SCIPvarGetName(vars[v]), SCIPvarGetLbGlobal(vars[v]));
                  SCIP_CALL( SCIPfixVar(scip, vars[v], SCIPvarGetLbGlobal(vars[v]), cutoff, &tightened) );
               }
               else
               {
                  assert(SCIPvarGetNLocksUp(vars[v]) == 1);
                  SCIPdebugMsg(scip, "fix <%s> to its upper bound (%g)\n", SCIPvarGetName(vars[v]), SCIPvarGetUbGlobal(vars[v]));
                  SCIP_CALL( SCIPfixVar(scip, vars[v], SCIPvarGetUbGlobal(vars[v]), cutoff, &tightened) );
               }

               if( *cutoff )
                  return SCIP_OKAY;
               if( tightened )
                  ++(*nfixedvars);
            }
            SCIPdebugMsg(scip, "### new stuffing fixed %d vars, tightened %d bounds\n", *nfixedvars - oldnfixedvars, *nchgbds - oldnchgbds);
         }
      }
   }

   return SCIP_OKAY;
}

/** applies full dual presolving on variables that only appear in linear constraints */
static
SCIP_RETCODE fullDualPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nchgbds             /**< pointer to count the number of bound changes */
   )
{
   SCIP_Real* redlb;
   SCIP_Real* redub;
   int* nlocksdown;
   int* nlocksup;
   SCIP_Bool* isimplint;
   SCIP_VAR** origvars;
   SCIP_VAR** vars;
   SCIP_VAR** conscontvars;
   int nvars;
   int nbinvars;
   int nintvars;
   int ncontvars;
   int v;
   int c;

   /* we calculate redundancy bounds with the following meaning:
    *   redlb[v] == k : if x_v >= k, we can always round x_v down to x_v == k without violating any constraint
    *   redub[v] == k : if x_v <= k, we can always round x_v up to x_v == k without violating any constraint
    * then:
    *   c_v >= 0 : x_v <= redlb[v] is feasible due to optimality
    *   c_v <= 0 : x_v >= redub[v] is feasible due to optimality
    */

   /* Additionally, we detect continuous variables that are implicitly integral.
    * A continuous variable j is implicit integral if it only has only +/-1 coefficients,
    * and all constraints (including the bounds as trivial constraints) in which:
    *   c_j > 0: the variable is down-locked,
    *   c_j < 0: the variable is up-locked,
    *   c_j = 0: the variable appears
    * have, apart from j, only integer variables with integral coefficients and integral sides.
    * This is because then, the value of the variable is either determined by one of its bounds or
    * by one of these constraints, and in all cases, the value of the variable is integral.
    */

   assert(scip != NULL);
   assert(nconss == 0 || conss != NULL);
   assert(nchgbds != NULL);
   assert(!SCIPinProbing(scip));

   /* get active variables */
   nvars = SCIPgetNVars(scip);
   origvars = SCIPgetVars(scip);

   /* if the problem is a pure binary program, nothing can be achieved by full dual presolve */
   nbinvars = SCIPgetNBinVars(scip);
   if( nbinvars == nvars )
      return SCIP_OKAY;

   /* get number of continuous variables */
   ncontvars = SCIPgetNContVars(scip);
   nintvars = nvars - ncontvars;

   /* copy the variable array since this array might change during the curse of this algorithm */
   nvars = nvars - nbinvars;
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, &(origvars[nbinvars]), nvars) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &redlb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redub, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksdown, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nlocksup, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isimplint, ncontvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conscontvars, ncontvars) );

   /* initialize redundancy bounds */
   for( v = 0; v < nvars; ++v )
   {
      assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_BINARY);
      redlb[v] = SCIPvarGetLbGlobal(vars[v]);
      redub[v] = SCIPvarGetUbGlobal(vars[v]);
   }
   BMSclearMemoryArray(nlocksdown, nvars);
   BMSclearMemoryArray(nlocksup, nvars);

   /* Initialize isimplint array: variable may be implied integer if rounded to their best bound they are integral.
    * We better not use SCIPisFeasIntegral() in these checks.
    */
   for( v = 0; v < ncontvars; v++ )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_Real lb;
      SCIP_Real ub;

      var = vars[v + nintvars - nbinvars];
      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      obj = SCIPvarGetObj(var);
      if( SCIPisZero(scip, obj) )
         isimplint[v] = (SCIPisInfinity(scip, -lb) || SCIPisIntegral(scip, lb)) && (SCIPisInfinity(scip, ub) || SCIPisIntegral(scip, ub));
      else
      {
         if( SCIPisPositive(scip, obj) )
            isimplint[v] = (SCIPisInfinity(scip, -lb) || SCIPisIntegral(scip, lb));
         else
         {
            assert(SCIPisNegative(scip, obj));
            isimplint[v] = (SCIPisInfinity(scip, ub) || SCIPisIntegral(scip, ub));
         }
      }
   }

   /* scan all constraints */
   for( c = 0; c < nconss; ++c )
   {
      /* we only need to consider constraints that have been locked (i.e., checked constraints or constraints that are
       * part of checked disjunctions)
       */
      if( SCIPconsIsLocked(conss[c]) )
      {
         SCIP_CONSDATA* consdata;
         SCIP_Bool lhsexists;
         SCIP_Bool rhsexists;
         SCIP_Bool hasimpliedpotential;
         SCIP_Bool integralcoefs;
         int nlockspos;
         int contvarpos;
         int nconscontvars;
         int i;

         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /* get number of times the constraint was locked */
         nlockspos = SCIPconsGetNLocksPos(conss[c]);

         /* we do not want to include constraints with locked negation (this would be too weird) */
         if( SCIPconsGetNLocksNeg(conss[c]) > 0 )
         {
            /* mark all continuous variables as not being implicit integral */
            for( i = 0; i < consdata->nvars; ++i )
            {
               SCIP_VAR* var;

               var = consdata->vars[i];
               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  int contv;
                  contv = SCIPvarGetProbindex(var) - nintvars;
                  assert(0 <= contv && contv < ncontvars); /* variable should be active due to applyFixings() */
                  isimplint[contv] = FALSE;
               }
            }
            continue;
         }

         /* check for existing sides */
         lhsexists = !SCIPisInfinity(scip, -consdata->lhs);
         rhsexists = !SCIPisInfinity(scip, consdata->rhs);

         /* count locks and update redundancy bounds */
         contvarpos = -1;
         nconscontvars = 0;
         hasimpliedpotential = FALSE;
         integralcoefs = !SCIPconsIsModifiable(conss[c]);

         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_VAR* var;
            SCIP_Real val;
            SCIP_Real minresactivity;
            SCIP_Real maxresactivity;
            SCIP_Real newredlb;
            SCIP_Real newredub;
            SCIP_Bool minisrelax;
            SCIP_Bool maxisrelax;
            SCIP_Bool isminsettoinfinity;
            SCIP_Bool ismaxsettoinfinity;
            int arrayindex;

            var = consdata->vars[i];
            val = consdata->vals[i];

            /* check if still all integer variables have integral coefficients */
            if( SCIPvarIsIntegral(var) )
               integralcoefs = integralcoefs && SCIPisIntegral(scip, val);

            /* we do not need to process binary variables */
            if( SCIPvarIsBinary(var) )
               continue;

            if( SCIPconsIsModifiable(conss[c]) )
            {
               minresactivity = -SCIPinfinity(scip);
               maxresactivity =  SCIPinfinity(scip);
               isminsettoinfinity = TRUE;
               ismaxsettoinfinity = TRUE;
            }
            else
            {
               /* calculate residual activity bounds if variable would be fixed to zero */
               consdataGetGlbActivityResiduals(scip, consdata, var, val, FALSE, &minresactivity, &maxresactivity,
                  &minisrelax, &maxisrelax, &isminsettoinfinity, &ismaxsettoinfinity);

               /* We called consdataGetGlbActivityResiduals() saying that we do not need a good relaxation,
                * so whenever we have a relaxed activity, it should be relaxed to +/- infinity.
                * This is needed, because we do not want to rely on relaxed finite resactivities.
                */
               assert((!minisrelax || isminsettoinfinity) && (!maxisrelax || ismaxsettoinfinity));

               /* check minresactivity for reliability */
               if( !isminsettoinfinity && SCIPisUpdateUnreliable(scip, minresactivity, consdata->lastglbminactivity) )
                  consdataGetReliableResidualActivity(scip, consdata, var, &minresactivity, TRUE, TRUE);

               /* check maxresactivity for reliability */
               if( !ismaxsettoinfinity && SCIPisUpdateUnreliable(scip, maxresactivity, consdata->lastglbmaxactivity) )
                  consdataGetReliableResidualActivity(scip, consdata, var, &maxresactivity, FALSE, TRUE);
            }

            arrayindex = SCIPvarGetProbindex(var) - nbinvars;

            assert(0 <= arrayindex && arrayindex < nvars); /* variable should be active due to applyFixings() */

            newredlb = redlb[arrayindex];
            newredub = redub[arrayindex];
            if( val > 0.0 )
            {
               if( lhsexists )
               {
                  /* lhs <= d*x + a*y, d > 0  ->  redundant in y if  x >= (lhs - min{a*y})/d */
                  nlocksdown[arrayindex] += nlockspos;
                  newredlb = (isminsettoinfinity ? SCIPinfinity(scip) : (consdata->lhs - minresactivity)/val);
               }
               if( rhsexists )
               {
                  /* d*x + a*y <= rhs, d > 0  ->  redundant in y if  x <= (rhs - max{a*y})/d */
                  nlocksup[arrayindex] += nlockspos;
                  newredub = (ismaxsettoinfinity ? -SCIPinfinity(scip) : (consdata->rhs - maxresactivity)/val);
               }
            }
            else
            {
               if( lhsexists )
               {
                  /* lhs <= d*x + a*y, d < 0  ->  redundant in y if  x <= (lhs - min{a*y})/d */
                  nlocksup[arrayindex] += nlockspos;
                  newredub = (isminsettoinfinity ? -SCIPinfinity(scip) : (consdata->lhs - minresactivity)/val);
               }
               if( rhsexists )
               {
                  /* d*x + a*y <= rhs, d < 0  ->  redundant in y if  x >= (rhs - max{a*y})/d */
                  nlocksdown[arrayindex] += nlockspos;
                  newredlb = (ismaxsettoinfinity ? SCIPinfinity(scip) : (consdata->rhs - maxresactivity)/val);
               }
            }

            /* if the variable is integer, we have to round the value to the next integral value */
            if( SCIPvarIsIntegral(var) )
            {
               if( !SCIPisInfinity(scip, newredlb) )
                  newredlb = SCIPceil(scip, newredlb);
               if( !SCIPisInfinity(scip, -newredub) )
                  newredub = SCIPfloor(scip, newredub);
            }

            /* update redundancy bounds */
            redlb[arrayindex] = MAX(redlb[arrayindex], newredlb);
            redub[arrayindex] = MIN(redub[arrayindex], newredub);

            /* collect the continuous variables of the constraint */
            if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
            {
               int contv;

               assert(nconscontvars < ncontvars);
               contvarpos = i;
               conscontvars[nconscontvars] = var;
               nconscontvars++;

               contv = SCIPvarGetProbindex(var) - nintvars;
               assert(0 <= contv && contv < ncontvars);
               hasimpliedpotential = hasimpliedpotential || isimplint[contv];
            }
         }

         /* update implied integer status of continuous variables */
         if( hasimpliedpotential )
         {
            if( nconscontvars > 1 || !integralcoefs )
            {
               /* there is more than one continuous variable or the integer variables have fractional coefficients:
                * none of the continuous variables is implied integer
                */
               for( i = 0; i < nconscontvars; i++ )
               {
                  int contv;
                  contv = SCIPvarGetProbindex(conscontvars[i]) - nintvars;
                  assert(0 <= contv && contv < ncontvars);
                  isimplint[contv] = FALSE;
               }
            }
            else
            {
               SCIP_VAR* var;
               SCIP_Real val;
               SCIP_Real absval;
               int contv;

               /* there is exactly one continuous variable and the integer variables have integral coefficients:
                * this is the interesting case, and we have to check whether the coefficient is +/-1 and the corresponding
                * side(s) of the constraint is integral
                */
               assert(nconscontvars == 1);
               assert(0 <= contvarpos && contvarpos < consdata->nvars);
               var = consdata->vars[contvarpos];
               val = consdata->vals[contvarpos];
               contv = SCIPvarGetProbindex(var) - nintvars;
               assert(0 <= contv && contv < ncontvars);
               assert(isimplint[contv]);

               absval = REALABS(val);
               if( !SCIPisEQ(scip, absval, 1.0) )
                  isimplint[contv] =  FALSE;
               else
               {
                  SCIP_Real obj;

                  obj = SCIPvarGetObj(var);
                  if( obj * val >= 0.0 && lhsexists )
                  {
                     /* the variable may be blocked by the constraint's left hand side */
                     isimplint[contv] = isimplint[contv] && SCIPisIntegral(scip, consdata->lhs);
                  }
                  if( obj * val <= 0.0 && rhsexists )
                  {
                     /* the variable may be blocked by the constraint's left hand side */
                     isimplint[contv] = isimplint[contv] && SCIPisIntegral(scip, consdata->rhs);
                  }
               }
            }
         }
      }
   }

   /* check if any bounds can be tightened due to optimality */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      assert(SCIPvarGetType(vars[v]) != SCIP_VARTYPE_BINARY);
      assert(SCIPvarGetNLocksDown(vars[v]) >= nlocksdown[v]);
      assert(SCIPvarGetNLocksUp(vars[v]) >= nlocksup[v]);

      var = vars[v];
      obj = SCIPvarGetObj(var);
      if( obj >= 0.0 )
      {
         /* making the variable as small as possible does not increase the objective:
          * check if all down locks of the variables are due to linear constraints;
          * if largest bound to make constraints redundant is -infinity, we better do nothing for numerical reasons
          */
         if( SCIPvarGetNLocksDown(var) == nlocksdown[v]
            && !SCIPisInfinity(scip, -redlb[v])
            && redlb[v] < SCIPvarGetUbGlobal(var) )
         {
            SCIP_Real ub;

            /* if x_v >= redlb[v], we can always round x_v down to x_v == redlb[v] without violating any constraint
             *  -> tighten upper bound to x_v <= redlb[v]
             */
            SCIPdebugMsg(scip, "variable <%s> only locked down in linear constraints: dual presolve <%s>[%.15g,%.15g] <= %.15g\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
               redlb[v]);
            SCIP_CALL( SCIPtightenVarUb(scip, var, redlb[v], FALSE, &infeasible, &tightened) );
            assert(!infeasible);

            ub = SCIPvarGetUbGlobal(var);
            redub[v] = MIN(redub[v], ub);
            if( tightened )
               (*nchgbds)++;
         }
      }
      if( obj <= 0.0 )
      {
         /* making the variable as large as possible does not increase the objective:
          * check if all up locks of the variables are due to linear constraints;
          * if smallest bound to make constraints redundant is +infinity, we better do nothing for numerical reasons
          */
         if( SCIPvarGetNLocksUp(var) == nlocksup[v]
            && !SCIPisInfinity(scip, redub[v])
            && redub[v] > SCIPvarGetLbGlobal(var) )
         {
            SCIP_Real lb;

            /* if x_v <= redub[v], we can always round x_v up to x_v == redub[v] without violating any constraint
             *  -> tighten lower bound to x_v >= redub[v]
             */
            SCIPdebugMsg(scip, "variable <%s> only locked up in linear constraints: dual presolve <%s>[%.15g,%.15g] >= %.15g\n",
               SCIPvarGetName(var), SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var),
               redub[v]);
            SCIP_CALL( SCIPtightenVarLb(scip, var, redub[v], FALSE, &infeasible, &tightened) );
            assert(!infeasible);

            lb = SCIPvarGetLbGlobal(var);
            redlb[v] = MAX(redlb[v], lb);
            if( tightened )
               (*nchgbds)++;
         }
      }
   }

   /* upgrade continuous variables to implied integers */
   for( v = nintvars - nbinvars; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool infeasible;

      var = vars[v];
      assert(var != NULL);

      assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
      assert(SCIPvarGetNLocksDown(var) >= nlocksdown[v]);
      assert(SCIPvarGetNLocksUp(var) >= nlocksup[v]);
      assert(0 <= v - nintvars + nbinvars && v - nintvars + nbinvars < ncontvars);

      /* we can only conclude implied integrality if the variable appears in no other constraint */
      if( isimplint[v - nintvars + nbinvars]
         && SCIPvarGetNLocksDown(var) == nlocksdown[v]
         && SCIPvarGetNLocksUp(var) == nlocksup[v] )
      {

         /* since we locally copied the variable array we can change the variable type immediately */
         SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_IMPLINT, &infeasible) );

         if( infeasible )
         {
            SCIPdebugMsg(scip, "infeasible upgrade of variable <%s> to integral type, domain is empty\n", SCIPvarGetName(var));
            *cutoff = TRUE;

            break;
         }

         SCIPdebugMsg(scip, "dual presolve: converting continuous variable <%s>[%g,%g] to implicit integer\n",
            SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &conscontvars);
   SCIPfreeBufferArray(scip, &isimplint);
   SCIPfreeBufferArray(scip, &nlocksup);
   SCIPfreeBufferArray(scip, &nlocksdown);
   SCIPfreeBufferArray(scip, &redub);
   SCIPfreeBufferArray(scip, &redlb);

   SCIPfreeBufferArray(scip, &vars);

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
   SCIP_RESULT*          result              /**< pointer to store the result of the enforcing call */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool checkrelmaxabs;
   SCIP_Bool violated;
   SCIP_Bool cutoff = FALSE;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   checkrelmaxabs = conshdlrdata->checkrelmaxabs;

   SCIPdebugMsg(scip, "Enforcement method of linear constraints for %s solution\n", sol == NULL ? "LP" : "relaxation");

   /* check for violated constraints
    * LP is processed at current node -> we can add violated linear constraints to the SCIP_LP
    */
   *result = SCIP_FEASIBLE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], sol, FALSE, checkrelmaxabs, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], &cutoff) );
         if ( cutoff )
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;
      }
   }

   /* check all obsolete linear constraints for feasibility */
   for( c = nusefulconss; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], sol, FALSE, checkrelmaxabs, &violated) );

      if( violated )
      {
         /* insert LP row as cut */
         SCIP_CALL( addRelaxation(scip, conss[c], &cutoff) );
         if ( cutoff )
            *result = SCIP_CUTOFF;
         else
            *result = SCIP_SEPARATED;
      }
   }
   
   SCIPdebugMsg(scip, "-> constraints checked, %s\n", *result == SCIP_FEASIBLE ? "all constraints feasible" : "infeasibility detected");

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);

   /* check for event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);
   assert(nconss == 0 || conss != NULL);

   /* catch events for the constraints */
   for( c = 0; c < nconss; ++c )
   {
      /* catch all events */
      SCIP_CALL( consCatchAllEvents(scip, conss[c], conshdlrdata->eventhdlr) );
   }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitLinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);

   /* check for event handler */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* drop events for the constraints */
   for( c = nconss - 1; c >= 0; --c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->eventdata != NULL )
      {
         /* drop all events */
         SCIP_CALL( consDropAllEvents(scip, conss[c], conshdlrdata->eventhdlr) );
         assert(consdata->eventdata == NULL);
      }
   }

   return SCIP_OKAY;

}

/** is constraint ranged row, i.e., -inf < lhs < rhs < inf? */
static
SCIP_Bool isRangedRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   assert(scip != NULL);

   return !(SCIPisEQ(scip, lhs, rhs)
      || SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs) );
}

/** is constraint ranged row, i.e., -inf < lhs < rhs < inf? */
static
SCIP_Bool isFiniteNonnegativeIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x                   /**< value */
   )
{
   assert(scip != NULL);

   return (!SCIPisInfinity(scip, x) && !SCIPisNegative(scip, x) && SCIPisIntegral(scip, x));
}

/** performs linear constraint type classification as used for MIPLIB
 *
 *  iterates through all linear constraints and stores relevant statistics in the linear constraint statistics \p linconsstats.
 *
 *  @note only constraints are iterated that belong to the linear constraint handler. If the problem has been presolved already,
 *  constraints that were upgraded to more special types such as, e.g., varbound constraints, will not be shown correctly anymore.
 *  Similarly, if specialized constraints were created through the API, these are currently not present.
 */
SCIP_RETCODE SCIPclassifyConstraintTypesLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LINCONSSTATS*    linconsstats        /**< linear constraint type classification */
   )
{
   int c;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;

   assert(scip != NULL);
   assert(linconsstats != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
   {
      conss = SCIPgetConss(scip);
      nconss = SCIPgetNConss(scip);
   }
   else
   {
      conss = SCIPconshdlrGetConss(conshdlr);
      nconss = SCIPconshdlrGetNConss(conshdlr);
   }


   /* reset linear constraint type classification */
   SCIPlinConsStatsReset(linconsstats);

   /* loop through all constraints */
   for( c = 0; c < nconss; c++ )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int i;

      /* get constraint */
      cons = conss[c];
      assert(cons != NULL);

      /* skip constraints that are not handled by the constraint handler */
      if( SCIPconsGetHdlr(cons) != conshdlr )
         continue;

      /* get constraint data */
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      rhs = consdata->rhs;
      lhs = consdata->lhs;

      /* merge multiples and delete variables with zero coefficient */
      SCIP_CALL( mergeMultiples(scip, cons) );
      for( i = 0; i < consdata->nvars; i++ )
      {
         assert(!SCIPisZero(scip, consdata->vals[i]));
      }

      /* is constraint of type SCIP_CONSTYPE_EMPTY? */
      if( consdata->nvars == 0 )
      {
         SCIPdebugMsg(scip, "classified as EMPTY: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_EMPTY, 1);

         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_FREE? */
      if( SCIPisInfinity(scip, rhs) && SCIPisInfinity(scip, -lhs) )
      {
         SCIPdebugMsg(scip, "classified as FREE: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_FREE, 1);

         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_SINGLETON? */
      if( consdata->nvars == 1 )
      {
         SCIPdebugMsg(scip, "classified as SINGLETON: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_SINGLETON, isRangedRow(scip, lhs, rhs) ? 2 : 1);

         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_AGGREGATION? */
      if( consdata->nvars == 2 && SCIPisEQ(scip, lhs, rhs) )
      {
         SCIPdebugMsg(scip, "classified as AGGREGATION: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_AGGREGATION, 1);

         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{VARBOUND,PRECEDENCE}? */
      if( consdata->nvars == 2 )
      {
         SCIP_LINCONSTYPE constype;

         /* precedence constraints have the same coefficient, but with opposite sign for the same variable type */
         if( SCIPisEQ(scip, consdata->vals[0], -consdata->vals[1])
               && SCIPvarGetType(consdata->vars[0]) == SCIPvarGetType(consdata->vars[1]))
         {
            constype = SCIP_LINCONSTYPE_PRECEDENCE;
            SCIPdebugMsg(scip, "classified as PRECEDENCE: ");
         }
         else
         {
            constype = SCIP_LINCONSTYPE_VARBOUND;
            SCIPdebugMsg(scip, "classified as VARBOUND: ");
         }
         SCIPdebugPrintCons(scip, cons, NULL);

         SCIPlinConsStatsIncTypeCount(linconsstats, constype, isRangedRow(scip, lhs, rhs) ? 2 : 1);

         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{SETPARTITION, SETPACKING, SETCOVERING, CARDINALITY, INVKNAPSACK}? */
      {
         SCIP_Real scale;
         SCIP_Real b;
         SCIP_Bool unmatched;
         int nnegbinvars;

         unmatched = FALSE;
         nnegbinvars = 0;

         scale = REALABS(consdata->vals[0]);

         /* scan through variables and detect if all variables are binary and have a coefficient +/-1 */
         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(consdata->vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(scip, SCIPvarGetLbGlobal(consdata->vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(scip, SCIPvarGetUbGlobal(consdata->vars[i]), 2.0);
            unmatched = unmatched || !SCIPisEQ(scip, REALABS(consdata->vals[i]), scale);

            if( consdata->vals[i] < 0.0 )
               nnegbinvars++;
         }

         if( !unmatched )
         {
            if( SCIPisEQ(scip, lhs, rhs) )
            {
               b = rhs/scale + nnegbinvars;
               if( SCIPisEQ(scip, 1.0, b) )
               {
                  SCIPdebugMsg(scip, "classified as SETPARTITION: ");
                  SCIPdebugPrintCons(scip, cons, NULL);
                  SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_SETPARTITION, 1);

                  continue;
               }
               else if( SCIPisIntegral(scip, b) && !SCIPisNegative(scip, b) )
               {
                  SCIPdebugMsg(scip, "classified as CARDINALITY: ");
                  SCIPdebugPrintCons(scip, cons, NULL);
                  SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_CARDINALITY, 1);

                  continue;
               }
            }

            /* compute right hand side divided by scale */
            if( !SCIPisInfinity(scip, rhs) )
               b = rhs/scale + nnegbinvars;
            else
               b = SCIPinfinity(scip);

            if( SCIPisEQ(scip, 1.0, b) )
            {
               SCIPdebugMsg(scip, "classified as SETPACKING: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_SETPACKING, 1);

               /* relax right hand side to prevent further classifications */
               rhs = SCIPinfinity(scip);
            }
            else if( !SCIPisInfinity(scip, b) && SCIPisIntegral(scip, b) && !SCIPisNegative(scip, b) )
            {
               SCIPdebugMsg(scip, "classified as INVKNAPSACK: ");
               SCIPdebugPrintCons(scip, cons, NULL);

               SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_INVKNAPSACK, 1);;

               /* relax right hand side to prevent further classifications */
               rhs = SCIPinfinity(scip);
            }

            if( !SCIPisInfinity(scip, lhs) )
               b = lhs/scale + nnegbinvars;
            else
               b = SCIPinfinity(scip);

            if( SCIPisEQ(scip, 1.0, b) )
            {
               SCIPdebugMsg(scip, "classified as SETCOVERING: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_SETCOVERING, 1);

               /* relax left hand side to prevent further classifications */
               lhs = -SCIPinfinity(scip);
            }

            /* if both sides are infinite at this point, no further classification is necessary for this constraint */
            if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
               continue;
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{EQKNAPSACK, BINPACKING, KNAPSACK}? */
      /* @todo If coefficients or rhs are not integral, we currently do not check
       * if the constraint could be scaled (finitely), such that they are.
       */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         b = rhs;
         unmatched = FALSE;
         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(consdata->vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(scip, SCIPvarGetLbGlobal(consdata->vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(scip, SCIPvarGetUbGlobal(consdata->vars[i]), 2.0);
            unmatched = unmatched || !SCIPisIntegral(scip, consdata->vals[i]);

            if( SCIPisNegative(scip, consdata->vals[i]) )
               b -= consdata->vals[i];
         }
         unmatched = unmatched || !isFiniteNonnegativeIntegral(scip, b);

         if( !unmatched )
         {
            if( SCIPisEQ(scip, lhs, rhs) )
            {
               SCIPdebugMsg(scip, "classified as EQKNAPSACK: ");
               SCIPdebugPrintCons(scip, cons, NULL);

               SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_EQKNAPSACK, 1);

               continue;
            }
            else
            {
               SCIP_Bool matched;

               matched = FALSE;
               for( i = 0; i < consdata->nvars && !matched; i++ )
               {
                  matched = matched || SCIPisEQ(scip, b, REALABS(consdata->vals[i]));
               }

               SCIPdebugMsg(scip, "classified as %s: ", matched ? "BINPACKING" : "KNAPSACK");
               SCIPdebugPrintCons(scip, cons, NULL);
               SCIPlinConsStatsIncTypeCount(linconsstats, matched ? SCIP_LINCONSTYPE_BINPACKING : SCIP_LINCONSTYPE_KNAPSACK, 1);
            }

            /* check if finite left hand side allows for a second classification, relax already used right hand side */
            if( SCIPisInfinity(scip, -lhs) )
               continue;
            else
               rhs = SCIPinfinity(scip);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{INTKNAPSACK}? */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         unmatched = FALSE;

         b = rhs;
         unmatched = unmatched || !isFiniteNonnegativeIntegral(scip, b);

         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(consdata->vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisNegative(scip, SCIPvarGetLbGlobal(consdata->vars[i]));
            unmatched = unmatched || !SCIPisIntegral(scip, consdata->vals[i]);
            unmatched = unmatched || SCIPisNegative(scip, consdata->vals[i]);
         }

         if( !unmatched )
         {
            SCIPdebugMsg(scip, "classified as INTKNAPSACK: ");
            SCIPdebugPrintCons(scip, cons, NULL);
            SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_INTKNAPSACK, 1);

            /* check if finite left hand side allows for a second classification, relax already used right hand side */
            if( SCIPisInfinity(scip, -lhs) )
               continue;
            else
               rhs = SCIPinfinity(scip);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{MIXEDBINARY}? */
      {
         SCIP_Bool unmatched;

         unmatched = FALSE;
         for( i = 0; i < consdata->nvars && !unmatched; i++ )
         {
            if( SCIPvarGetType(consdata->vars[i]) != SCIP_VARTYPE_CONTINUOUS
               && (SCIPisLE(scip, SCIPvarGetLbGlobal(consdata->vars[i]), -1.0)
                  || SCIPisGE(scip, SCIPvarGetUbGlobal(consdata->vars[i]), 2.0)) )
               unmatched = TRUE;
         }

         if( !unmatched )
         {
            SCIPdebugMsg(scip, "classified as MIXEDBINARY (%d): ", isRangedRow(scip, lhs, rhs) ? 2 : 1);
            SCIPdebugPrintCons(scip, cons, NULL);
            SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_MIXEDBINARY, isRangedRow(scip, lhs, rhs) ? 2 : 1);

            continue;
         }
      }

      /* no special structure detected */
      SCIPdebugMsg(scip, "classified as GENERAL: ");
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPlinConsStatsIncTypeCount(linconsstats, SCIP_LINCONSTYPE_GENERAL, isRangedRow(scip, lhs, rhs) ? 2 : 1);
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreLinear)
{  /*lint --e{715}*/
   int c;
#ifdef SCIP_STATISTIC
   SCIP_CONSHDLRDATA* conshdlrdata;
   int ngoodconss;
   int nallconss;
#endif

   /* delete all linear constraints that were upgraded to a more specific constraint type;
    * make sure, only active variables remain in the remaining constraints
    */
   assert(scip != NULL);

#ifdef SCIP_STATISTIC
   /* count number of well behaved linear constraints */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ngoodconss = 0;
   nallconss = 0;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->upgraded )
         continue;

      nallconss++;

      consdataRecomputeMaxActivityDelta(scip, consdata);

      if( SCIPisLT(scip, consdata->maxactdelta, conshdlrdata->maxeasyactivitydelta) )
         ngoodconss++;
   }
   if( nallconss )
   {
      SCIPstatisticMessage("below threshold: %d / %d ratio= %g\n", ngoodconss, nallconss, (100.0 * ngoodconss / nallconss));
   }
#endif

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      if( SCIPconsIsDeleted(conss[c]) )
         continue;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->upgraded )
      {
         /* this is no problem reduction, because the upgraded constraint was added to the problem before, and the
          * (redundant) linear constraint was only kept in order to support presolving the the linear constraint handler
          */
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      }
      else
      {
         /* since we are not allowed to detect infeasibility in the exitpre stage, we dont give an infeasible pointer */
         SCIP_CALL( applyFixings(scip, conss[c], NULL) );
      }
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolLinear)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   /* if this is a restart, convert cutpool rows into linear constraints */
   if( restart )
   {
      int ncutsadded;

      ncutsadded = 0;

      /* create out of all active cuts in cutpool linear constraints */
      SCIP_CALL( SCIPconvertCutsToConss(scip, NULL, NULL, TRUE, &ncutsadded) );

      if( ncutsadded > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "(restart) converted %d cuts from the global cut pool into linear constraints\n", ncutsadded);
         /* an extra blank line should be printed separately since the buffer message handler only handles up to one
          * line correctly
          */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n");
      }
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveLinear)
{  /*lint --e{715}*/

   if( SCIPconsIsDeleted(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_CONSDATA* consdata;

      assert(scip != NULL);
      assert(conshdlr != NULL);
      assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

      /* get constraint data */
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* check for event handler */
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->eventhdlr != NULL);

      /* free event data */
      if( consdata->eventdata != NULL )
      {
         /* drop bound change events of variables */
         SCIP_CALL( consDropAllEvents(scip, cons, conshdlrdata->eventhdlr) );
      }
      assert(consdata->eventdata == NULL);
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   if( (*consdata)->eventdata != NULL )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* drop all events */
      SCIP_CALL( consDropAllEvents(scip, cons, conshdlrdata->eventhdlr) );
      assert((*consdata)->eventdata == NULL);
   }

   /* free linear constraint */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   /*debugMsg(scip, "Trans method of linear constraints\n");*/

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->row == NULL);  /* in original problem, there cannot be LP rows */

   /* create linear constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars, sourcedata->vals, sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpLinear)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   *infeasible = FALSE;

   for( c = 0; c < nconss && !(*infeasible); ++c )
   {
      assert(SCIPconsIsInitial(conss[c]));
      SCIP_CALL( addRelaxation(scip, conss[c], infeasible) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real loclowerbound;
   SCIP_Real glblowerbound;
   SCIP_Real cutoffbound;
   SCIP_Real maxbound;
   SCIP_Bool separatecards;
   SCIP_Bool cutoff;
   int c;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   /*debugMsg(scip, "Sepa method of linear constraints\n");*/

   *result = SCIP_DIDNOTRUN;

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   /* check if we want to produce knapsack cardinality cuts at this node */
   loclowerbound = SCIPgetLocalLowerbound(scip);
   glblowerbound = SCIPgetLowerbound(scip);
   cutoffbound = SCIPgetCutoffbound(scip);
   maxbound = glblowerbound + conshdlrdata->maxcardbounddist * (cutoffbound - glblowerbound);
   separatecards = SCIPisLE(scip, loclowerbound, maxbound);
   separatecards = separatecards && (SCIPgetNLPBranchCands(scip) > 0);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;
   cutoff = FALSE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss && ncuts < maxsepacuts && !cutoff; ++c )
   {
      /*debugMsg(scip, "separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      SCIP_CALL( separateCons(scip, conss[c], conshdlrdata, NULL, separatecards, conshdlrdata->separateall, &ncuts, &cutoff) );
   }

   /* adjust return value */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* combine linear constraints to get more cuts */
   /**@todo further cuts of linear constraints */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   int depth;
   int nrounds;
   int maxsepacuts;
   int ncuts;
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   depth = SCIPgetDepth(scip);
   nrounds = SCIPgetNSepaRounds(scip);

   /*debugMsg(scip, "Sepa method of linear constraints\n");*/

   *result = SCIP_DIDNOTRUN;

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && conshdlrdata->maxroundsroot >= 0 && nrounds >= conshdlrdata->maxroundsroot)
      || (depth > 0 && conshdlrdata->maxrounds >= 0 && nrounds >= conshdlrdata->maxrounds) )
      return SCIP_OKAY;

   /* get the maximal number of cuts allowed in a separation round */
   maxsepacuts = (depth == 0 ? conshdlrdata->maxsepacutsroot : conshdlrdata->maxsepacuts);

   *result = SCIP_DIDNOTFIND;
   ncuts = 0;
   cutoff = FALSE;

   /* check all useful linear constraints for feasibility */
   for( c = 0; c < nusefulconss && ncuts < maxsepacuts && !cutoff; ++c )
   {
      /*debugMsg(scip, "separating linear constraint <%s>\n", SCIPconsGetName(conss[c]));*/
      SCIP_CALL( separateCons(scip, conss[c], conshdlrdata, sol, TRUE, conshdlrdata->separateall, &ncuts, &cutoff) );
   }

   /* adjust return value */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   /* combine linear constraints to get more cuts */
   /**@todo further cuts of linear constraints */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpLinear)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxLinear)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, result) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool checkrelmaxabs;
   SCIP_Bool violated;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   checkrelmaxabs = conshdlrdata->checkrelmaxabs;

   SCIPdebugMsg(scip, "Enfops method of linear constraints\n");

   /* if the solution is infeasible anyway due to objective value, skip the enforcement */
   if( objinfeasible )
   {
      SCIPdebugMsg(scip, "-> pseudo solution is objective infeasible, return.\n");

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* check all linear constraints for feasibility */
   violated = FALSE;
   for( c = 0; c < nconss && !violated; ++c )
   {
      SCIP_CALL( checkCons(scip, conss[c], NULL, TRUE, checkrelmaxabs, &violated) );
   }

   if( violated )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   SCIPdebugMsg(scip, "-> constraints checked, %s\n", *result == SCIP_FEASIBLE ? "all constraints feasible" : "infeasibility detected");

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool checkrelmaxabs;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   checkrelmaxabs = conshdlrdata->checkrelmaxabs;

   /*debugMsg(scip, "Check method of linear constraints\n");*/

   /* check all linear constraints for feasibility */
   for( c = 0; c < nconss && (*result == SCIP_FEASIBLE || completely); ++c )
   {
      SCIP_Bool violated = FALSE;
      SCIP_CALL( checkCons(scip, conss[c], sol, checklprows, checkrelmaxabs, &violated) );

      if( violated )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIP_CONSDATA* consdata;
            SCIP_Real activity;

            consdata = SCIPconsGetData(conss[c]);
            assert( consdata != NULL);

            activity = consdataGetActivity(scip, consdata, sol);

            SCIP_CALL( consPrintConsSol(scip, conss[c], sol, NULL ) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( activity == SCIP_INVALID ) /*lint !e777*/
               SCIPinfoMessage(scip, NULL, "activity invalid due to positive and negative infinity contributions\n");
            else if( SCIPisFeasLT(scip, activity, consdata->lhs) )
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhs - activity);
            else if( SCIPisFeasGT(scip, activity, consdata->rhs) )
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", activity - consdata->rhs);
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool rangedrowpropagation = FALSE;
   SCIP_Bool tightenbounds;
   SCIP_Bool cutoff;

   int nchgbds;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /*debugMsg(scip, "Prop method of linear constraints\n");*/

   /* check, if we want to tighten variable's bounds (in probing, we always want to tighten the bounds) */
   if( SCIPinProbing(scip) )
      tightenbounds = TRUE;
   else
   {
      int depth;
      int propfreq;
      int tightenboundsfreq;
      int rangedrowfreq;

      depth = SCIPgetDepth(scip);
      propfreq = SCIPconshdlrGetPropFreq(conshdlr);
      tightenboundsfreq = propfreq * conshdlrdata->tightenboundsfreq;
      tightenbounds = (conshdlrdata->tightenboundsfreq >= 0)
         && ((tightenboundsfreq == 0 && depth == 0) || (tightenboundsfreq >= 1 && (depth % tightenboundsfreq == 0)));

      /* check if we want to do ranged row propagation */
      rangedrowpropagation = conshdlrdata->rangedrowpropagation;
      rangedrowpropagation = rangedrowpropagation && !SCIPinRepropagation(scip);
      rangedrowpropagation = rangedrowpropagation && (depth <= conshdlrdata->rangedrowmaxdepth);
      rangedrowfreq = propfreq * conshdlrdata->rangedrowfreq;
      rangedrowpropagation = rangedrowpropagation && (depth % rangedrowfreq == 0);
   }

   cutoff = FALSE;
   nchgbds = 0;

   /* process constraints marked for propagation */
   for( i = 0; i < nmarkedconss && !cutoff; i++ )
   {
      SCIP_CALL( SCIPunmarkConsPropagate(scip, conss[i]) );
      SCIP_CALL( propagateCons(scip, conss[i], tightenbounds, rangedrowpropagation,
            conshdlrdata->maxeasyactivitydelta, conshdlrdata->sortvars, &cutoff, &nchgbds) );
   }

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


#define MAXCONSPRESOLROUNDS 10
/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolLinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Bool minactisrelax;
   SCIP_Bool maxactisrelax;
   SCIP_Bool cutoff;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldndelconss;
   int oldnupgdconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int firstchange;
   int firstupgradetry;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /*debugMsg(scip, "Presol method of linear constraints\n");*/

   /* remember old preprocessing counters */
   cutoff = FALSE;
   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnupgdconss = *nupgdconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process single constraints */
   firstchange = INT_MAX;
   firstupgradetry = INT_MAX;
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      int npresolrounds;
      SCIP_Bool infeasible;

      infeasible = FALSE;

      cons = conss[c];
      assert(SCIPconsIsActive(cons));
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( consdata->eventdata == NULL )
      {
         /* catch bound change events of variables */
         SCIP_CALL( consCatchAllEvents(scip, cons, conshdlrdata->eventhdlr) );
         assert(consdata->eventdata != NULL);
      }

      /* constraint should not be already presolved in the initial round */
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || SCIPconsIsMarkedPropagate(cons));
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || consdata->boundstightened == 0);
      assert(SCIPgetNRuns(scip) > 0 || nrounds > 0 || !consdata->presolved);
      assert(!SCIPconsIsMarkedPropagate(cons) || !consdata->presolved);

      /* incorporate fixings and aggregations in constraint */
      SCIP_CALL( applyFixings(scip, cons, &infeasible) );

      if( infeasible )
      {
         SCIPdebugMsg(scip, " -> infeasible fixing\n");
         cutoff = TRUE;
         break;
      }

      assert(consdata->removedfixings);

      /* we can only presolve linear constraints, that are not modifiable */
      if( SCIPconsIsModifiable(cons) )
         continue;

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* remember the first constraint that was not yet tried to be upgraded, to begin the next upgrading round with */
      if( firstupgradetry == INT_MAX && !consdata->upgradetried )
         firstupgradetry = c;

      /* check, if constraint is already preprocessed */
      if( consdata->presolved )
         continue;

      assert(SCIPconsIsActive(cons));

      SCIPdebugMsg(scip, "presolving linear constraint <%s>\n", SCIPconsGetName(cons));
      SCIPdebugPrintCons(scip, cons, NULL);

      /* apply presolving as long as possible on the single constraint (however, abort after a certain number of rounds
       * to avoid nearly infinite cycling due to very small bound changes)
       */
      npresolrounds = 0;
      while( !consdata->presolved && npresolrounds < MAXCONSPRESOLROUNDS && !SCIPisStopped(scip) )
      {
         assert(!cutoff);
         npresolrounds++;

         /* mark constraint being presolved and propagated */
         consdata->presolved = TRUE;
         SCIP_CALL( SCIPunmarkConsPropagate(scip, cons) );

         /* normalize constraint */
         SCIP_CALL( normalizeCons(scip, cons) );

         /* tighten left and right hand side due to integrality */
         SCIP_CALL( tightenSides(scip, cons, nchgsides) );

         /* check bounds */
         if( SCIPisFeasGT(scip, consdata->lhs, consdata->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is infeasible: sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
            cutoff = TRUE;
            break;
         }

         /* tighten variable's bounds */
         SCIP_CALL( tightenBounds(scip, cons, conshdlrdata->maxeasyactivitydelta, conshdlrdata->sortvars, &cutoff, nchgbds) );
         if( cutoff )
            break;

         /* check for fixed variables */
         SCIP_CALL( fixVariables(scip, cons, &cutoff, nfixedvars) );
         if( cutoff )
            break;

         /* check constraint for infeasibility and redundancy */
         consdataGetActivityBounds(scip, consdata, TRUE, &minactivity, &maxactivity, &minactisrelax, &maxactisrelax);
         if( SCIPisFeasGT(scip, minactivity, consdata->rhs) || SCIPisFeasLT(scip, maxactivity, consdata->lhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is infeasible: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            cutoff = TRUE;
            break;
         }
         else if( SCIPisFeasGE(scip, minactivity, consdata->lhs) && SCIPisFeasLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( SCIPdelCons(scip, cons) );
            assert(!SCIPconsIsActive(cons));

            if( !consdata->upgraded )
               (*ndelconss)++;
            break;
         }
         else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasGE(scip, minactivity, consdata->lhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> left hand side is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( chgLhs(scip, cons, -SCIPinfinity(scip)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisFeasLE(scip, maxactivity, consdata->rhs) )
         {
            SCIPdebugMsg(scip, "linear constraint <%s> right hand side is redundant: activitybounds=[%.15g,%.15g], sides=[%.15g,%.15g]\n",
               SCIPconsGetName(cons), minactivity, maxactivity, consdata->lhs, consdata->rhs);
            SCIP_CALL( chgRhs(scip, cons, SCIPinfinity(scip)) );
            if( !consdata->upgraded )
               (*nchgsides)++;
         }
         assert(consdata->nvars >= 1); /* otherwise, it should be redundant or infeasible */

         /* handle empty constraint */
         if( consdata->nvars == 0 )
         {
            if( SCIPisFeasGT(scip, consdata->lhs, consdata->rhs) )
            {
               SCIPdebugMsg(scip, "empty linear constraint <%s> is infeasible: sides=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
               cutoff = TRUE;
            }
            else
            {
               SCIPdebugMsg(scip, "empty linear constraint <%s> is redundant: sides=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
               SCIP_CALL( SCIPdelCons(scip, cons) );
               assert(!SCIPconsIsActive(cons));

               if( !consdata->upgraded )
                  (*ndelconss)++;
            }
            break;
         }

         /* reduce big-M coefficients, that make the constraint redundant if the variable is on a bound */
         SCIP_CALL( consdataTightenCoefs(scip, cons, nchgcoefs, nchgsides) );

         /* try to simplify inequalities */
         if( conshdlrdata->simplifyinequalities )
         {
            SCIP_CALL( simplifyInequalities(scip, cons, nchgcoefs, nchgsides) );
         }

         /* aggregation variable in equations */
         if( conshdlrdata->aggregatevariables )
         {
            SCIP_CALL( aggregateVariables(scip, cons, &cutoff, nfixedvars, naggrvars, ndelconss) );
            if( cutoff )
               break;
         }
      }

      if( !cutoff && !SCIPisStopped(scip) )
      {
         /* perform ranged row propagation */
         if( conshdlrdata->rangedrowpropagation )
         {
            int lastnfixedvars;

            lastnfixedvars = *nfixedvars;

            SCIP_CALL( rangedRowPropagation(scip, cons, &cutoff, nfixedvars, nchgbds, naddconss) );
            if( !cutoff )
            {
               if( lastnfixedvars < *nfixedvars )
               {
                  SCIP_CALL( applyFixings(scip, cons, &cutoff) );
               }
            }
         }

         /* extract cliques from constraint */
	 if( !cutoff && SCIPconsIsActive(cons) )
	 {
            SCIP_CALL( extractCliques(scip, cons, conshdlrdata->maxeasyactivitydelta, conshdlrdata->sortvars,
                  nfixedvars, nchgbds, &cutoff) );

            /* check if the constraint got redundant or infeasible */
            if( !cutoff && SCIPconsIsActive(cons) && consdata->nvars == 0 )
            {
               if( SCIPisFeasGT(scip, consdata->lhs, consdata->rhs) )
               {
                  SCIPdebugMsg(scip, "empty linear constraint <%s> is infeasible: sides=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
                  cutoff = TRUE;
               }
               else
               {
                  SCIPdebugMsg(scip, "empty linear constraint <%s> is redundant: sides=[%.15g,%.15g]\n",
                     SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
                  SCIP_CALL( SCIPdelCons(scip, cons) );
                  assert(!SCIPconsIsActive(cons));

                  if( !consdata->upgraded )
                     (*ndelconss)++;
               }
            }
         }

         /* convert special equalities */
         if( !cutoff && SCIPconsIsActive(cons) )
         {
            SCIP_CALL( convertEquality(scip, cons, conshdlrdata, &cutoff, nfixedvars, naggrvars, ndelconss) );
         }

         /* apply dual presolving for variables that appear in only one constraint */
         if( !cutoff && SCIPconsIsActive(cons) && conshdlrdata->dualpresolving && SCIPallowDualReds(scip) )
         {
            SCIP_CALL( dualPresolve(scip, cons, &cutoff, nfixedvars, naggrvars, ndelconss) );
         }

         /* check if an inequality is parallel to the objective function */
         if( !cutoff && SCIPconsIsActive(cons) )
         {
            SCIP_CALL( checkParallelObjective(scip, cons, conshdlrdata) );
         }

         /* remember the first changed constraint to begin the next aggregation round with */
         if( firstchange == INT_MAX && consdata->changed )
            firstchange = c;

         /* remember the first constraint that was not yet tried to be upgraded, to begin the next upgrading round with */
         if( firstupgradetry == INT_MAX && !consdata->upgradetried )
            firstupgradetry = c;
      }

      /* singleton column stuffing */
      if( !cutoff && SCIPconsIsActive(cons) && SCIPconsIsChecked(cons) &&
         (conshdlrdata->singletonstuffing || conshdlrdata->singlevarstuffing) && SCIPallowDualReds(scip) )
      {
         SCIP_CALL( presolStuffing(scip, cons, conshdlrdata->singletonstuffing,
               conshdlrdata->singlevarstuffing, &cutoff, nfixedvars, nchgbds) );

         /* handle empty constraint */
         if( consdata->nvars == 0 )
         {
            if( SCIPisFeasGT(scip, consdata->lhs, consdata->rhs) )
            {
               SCIPdebugMsg(scip, "empty linear constraint <%s> is infeasible: sides=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
               cutoff = TRUE;
            }
            else
            {
               SCIPdebugMsg(scip, "empty linear constraint <%s> is redundant: sides=[%.15g,%.15g]\n",
                  SCIPconsGetName(cons), consdata->lhs, consdata->rhs);
               SCIP_CALL( SCIPdelCons(scip, cons) );
               assert(!SCIPconsIsActive(cons));

               if( !consdata->upgraded )
                  (*ndelconss)++;
            }
            break;
         }
      }
   }

   /* process pairs of constraints: check them for redundancy and try to aggregate them;
    * only apply this expensive procedure in exhaustive presolving timing
    */
   if( !cutoff && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 && (conshdlrdata->presolusehashing || conshdlrdata->presolpairwise) && !SCIPisStopped(scip) )
   {
      assert(firstchange >= 0);

      if( firstchange < nconss && conshdlrdata->presolusehashing )
      {
         /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
         SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, &cutoff,
               ndelconss, nchgsides) );
      }

      if( firstchange < nconss && conshdlrdata->presolpairwise )
      {
         SCIP_CONS** usefulconss;
         int nusefulconss;
         int firstchangenew;
         SCIP_Longint npaircomparisons;

         npaircomparisons = 0;
         oldndelconss = *ndelconss;
         oldnchgsides = *nchgsides;
         oldnchgcoefs = *nchgcoefs;

         /* allocate temporary memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &usefulconss, nconss) );

         nusefulconss = 0;
         firstchangenew = -1;
         for( c = 0; c < nconss; ++c )
         {
            /* update firstchange */
            if( c == firstchange )
               firstchangenew = nusefulconss;

            /* ignore inactive and modifiable constraints */
            if( !SCIPconsIsActive(conss[c]) || SCIPconsIsModifiable(conss[c]) )
               continue;

            usefulconss[nusefulconss] = conss[c];
            ++nusefulconss;
         }
         firstchange = firstchangenew;
         assert(firstchangenew >= 0 && firstchangenew <= nusefulconss);

         for( c = firstchange; c < nusefulconss && !cutoff && !SCIPisStopped(scip); ++c )
         {
            /* constraint has become inactive or modifiable during pairwise presolving */
            if( usefulconss[c] == NULL )
               continue;

            npaircomparisons += (SCIPconsGetData(conss[c])->changed) ? c : (c - firstchange); /*lint !e776*/

            assert(SCIPconsIsActive(usefulconss[c]) && !SCIPconsIsModifiable(usefulconss[c]));
            SCIP_CALL( preprocessConstraintPairs(scip, usefulconss, firstchange, c, conshdlrdata->maxaggrnormscale,
                  &cutoff, ndelconss, nchgsides, nchgcoefs) );

            if( npaircomparisons > conshdlrdata->nmincomparisons )
            {
               assert(npaircomparisons > 0);
               if( ((*ndelconss - oldndelconss) + (*nchgsides - oldnchgsides)/2.0 + (*nchgcoefs - oldnchgcoefs)/10.0) / ((SCIP_Real) npaircomparisons) < conshdlrdata->mingainpernmincomp )
                  break;
               oldndelconss = *ndelconss;
               oldnchgsides = *nchgsides;
               oldnchgcoefs = *nchgcoefs;
               npaircomparisons = 0;
            }
         }
         /* free temporary memory */
         SCIPfreeBufferArray(scip, &usefulconss);
      }
   }

   /* before upgrading, check whether we can apply some additional dual presolving, because a variable only appears
    * in linear constraints and we therefore have full information about it
    */
   if( !cutoff && firstupgradetry < nconss
      && *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars && *nchgbds == oldnchgbds && *ndelconss == oldndelconss
      && *nupgdconss == oldnupgdconss && *nchgcoefs == oldnchgcoefs && *nchgsides == oldnchgsides
      )
   {
      if( conshdlrdata->dualpresolving && SCIPallowDualReds(scip) && !SCIPisStopped(scip) )
      {
         SCIP_CALL( fullDualPresolve(scip, conss, nconss, &cutoff, nchgbds) );
      }
   }

   /* try to upgrade constraints into a more specific constraint type;
    * only upgrade constraints, if no reductions were found in this round (otherwise, the linear constraint handler
    * may find additional reductions before giving control away to other (less intelligent?) constraint handlers)
    */
   if( !cutoff && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 && SCIPisPresolveFinished(scip) )
   {
      for( c = firstupgradetry; c < nconss && !SCIPisStopped(scip); ++c )
      {
         cons = conss[c];

         /* don't upgrade modifiable constraints */
         if( SCIPconsIsModifiable(cons) )
            continue;

         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);

         /* only upgrade completely presolved constraints, that changed since the last upgrading call */
         if( consdata->upgradetried )
            continue;
         /* @todo force that upgrade will be performed later? */
         if( !consdata->presolved )
            continue;

         consdata->upgradetried = TRUE;
         if( SCIPconsIsActive(cons) )
         {
            SCIP_CONS* upgdcons;

            SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
            if( upgdcons != NULL )
            {
               /* add the upgraded constraint to the problem */
               SCIP_CALL( SCIPaddCons(scip, upgdcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &upgdcons) );
               (*nupgdconss)++;

               /* mark the linear constraint being upgraded and to be removed after presolving;
                * don't delete it directly, because it may help to preprocess other linear constraints
                */
               assert(!consdata->upgraded);
               consdata->upgraded = TRUE;

               /* delete upgraded inequalities immediately;
                * delete upgraded equalities, if we don't need it anymore for aggregation and redundancy checking
                */
               if( SCIPisLT(scip, consdata->lhs, consdata->rhs)
                  || !conshdlrdata->presolpairwise
                  || (conshdlrdata->maxaggrnormscale == 0.0) )
               {
                  SCIP_CALL( SCIPdelCons(scip, cons) );
               }
            }
         }
      }
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds || *ndelconss > oldndelconss
      || *nupgdconss > oldnupgdconss || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropLinear)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(cons != NULL);
   assert(result != NULL);

   SCIP_CALL( resolvePropagation(scip, cons, infervar, intToInferInfo(inferinfo), boundtype, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool haslhs;
   SCIP_Bool hasrhs;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslhs = !SCIPisInfinity(scip, -consdata->lhs);
   hasrhs = !SCIPisInfinity(scip, consdata->rhs);

   /* update rounding locks of every single variable */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPisPositive(scip, consdata->vals[i]) )
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
         }
         if( hasrhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos, nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}


/** variable deletion method of constraint handler */
static
SCIP_DECL_CONSDELVARS(consDelvarsLinear)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   if( nconss > 0 )
   {
      SCIP_CALL( performVarDeletions(scip, conshdlr, conss, nconss) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintLinear)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file) );

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyLinear)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_Real* sourcecoefs;
   const char* consname;
   int nvars;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsLinear(sourcescip, sourcecons);
   sourcecoefs = SCIPgetValsLinear(sourcescip, sourcecons); 
   nvars = SCIPgetNVarsLinear(sourcescip, sourcecons);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, nvars, sourcevars, sourcecoefs,
         SCIPgetLhsLinear(sourcescip, sourcecons), SCIPgetRhsLinear(sourcescip, sourcecons), varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   assert(cons != NULL || *valid == FALSE);

   /* @todo should also the checkabsolute flag of the constraint be copied? */

   return SCIP_OKAY;
}

/* find operators '<=', '==', '>=', [free] in input string and return those places. There should only be one operator,
 * except for ranged rows for which exactly two operators '<=' must be present
 */
static
SCIP_RETCODE findOperators(
   const char*           str,                /**< null terminated input string */
   char**                firstoperator,      /**< pointer to store the string starting at the first operator */
   char**                secondoperator,     /**< pointer to store the string starting at the second operator */
   SCIP_Bool*            success             /**< pointer to store if the line contains a valid operator order */
   )
{
   char* curr;

   assert(str != NULL);
   assert(firstoperator != NULL);
   assert(secondoperator != NULL);

   *firstoperator = NULL;
   *secondoperator = NULL;

   curr = (char*)str;
   *success = TRUE;

   /* loop over the input string to find all operators */
   while( *curr && *success )
   {
      SCIP_Bool found = FALSE;
      int increment = 1;

      /* try if we found a possible operator */
      switch( *curr )
      {
         case '<':
         case '=':
         case '>':

            /* check if the two characters curr[0,1] form an operator together */
            if( curr[1] == '=' )
            {
               found = TRUE;

               /* update increment to continue after this operator */
               increment = 2;
            }
            break;
         case '[':
            if( strncmp(curr, "[free]", 6) == 0 )
            {
               found = TRUE;

               /* update increment to continue after this operator */
               increment = 6;
            }
            break;
         default:
            break;

      }

      /* assign the found operator to the first or second pointer and check for violations of the linear constraint grammar */
      if( found )
      {
         if( *firstoperator == NULL )
         {
            *firstoperator = curr;
         }
         else
         {
            if( *secondoperator != NULL )
            {
               SCIPerrorMessage("Found more than two operators in line %s\n", str);
               *success = FALSE;
            }
            else if( strncmp(*firstoperator, "<=", 2) != 0 )
            {
               SCIPerrorMessage("Two operators in line that is not a ranged row: %s", str);
               *success = FALSE;
            }
            else if( strncmp(curr, "<=", 2) != 0 )
            {
               SCIPerrorMessage("Bad second operator, expected ranged row specification: %s", str);
               *success = FALSE;
            }

            *secondoperator = curr;
         }
      }

      curr += increment;
   }

   /* check if we did find at least one operator */
   if( *success )
   {
      if( *firstoperator == NULL )
      {
         SCIPerrorMessage("Could not find any operator in line %s\n", str);
         *success = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseLinear)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int        nvars;
   int        coefssize;
   int        requsize;
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   char*      endptr;
   char*      firstop;
   char*      secondop;
   SCIP_Bool  operatorsuccess;
   char*      lhsstrptr;
   char*      rhsstrptr;
   char*      varstrptr;

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

   /* find operators in the line first, all other remaining parsing depends on occurence of the operators '<=', '>=', '==',
    * and the special word [free]
    */
   SCIP_CALL( findOperators(str, &firstop, &secondop, &operatorsuccess) );

   /* if the grammar is not valid for parsing a linear constraint, return */
   if( ! operatorsuccess )
      return SCIP_OKAY;

   varstrptr = (char *)str;
   lhsstrptr = rhsstrptr = NULL;
   assert(firstop != NULL);

   /* assign the strings for parsing the left hand side, right hand side, and the linear variable sum */
   switch( *firstop )
   {
      case '<':
         assert(firstop[1] == '=');
         /* we have ranged row lhs <= a_1 x_1 + ... + a_n x_n <= rhs */
         if( secondop != NULL )
         {
            assert(secondop[0] == '<' && secondop[1] == '=');
            lhsstrptr = (char *)str;
            varstrptr = firstop + 2;
            rhsstrptr = secondop + 2;
         }
         else
         {
            /* we have an inequality with infinite left hand side a_1 x_1 + ... + a_n x_n <= rhs */
            lhsstrptr = NULL;
            varstrptr = (char *)str;
            rhsstrptr = firstop + 2;
         }
         break;
      case '>':
         assert(firstop[1] == '=');
         assert(secondop == NULL);
         /* we have a_1 x_1 + ... + a_n x_n >= lhs */
         lhsstrptr = firstop + 2;
         break;
      case '=':
         assert(firstop[1] == '=');
         assert(secondop == NULL);
         /* we have a_1 x_1 + ... + a_n x_n == lhs (rhs) */
         rhsstrptr = firstop + 2;
         lhsstrptr = firstop + 2;
         break;
      case '[':
         assert(strncmp(firstop, "[free]", 6) == 0);
         assert(secondop == NULL);
         /* nothing to assign in case of a free a_1 x_1 + ... + a_n x_n [free] */
         break;
      default:
         /* it should not be possible that a different character appears in that position */
         SCIPerrorMessage("Parsing has wrong operator character '%c', should be one of <=>[", *firstop);
         return SCIP_READERROR;
   }

   /* parse left hand side, if necessary */
   if( lhsstrptr != NULL )
   {
      if( ! SCIPparseReal(scip, lhsstrptr, &lhs, &endptr) )
      {
         SCIPerrorMessage("error parsing left hand side number from <%s>\n", lhsstrptr);
         return SCIP_OKAY;
      }

      /* in case of an equation, assign the left also to the right hand side */
      if( rhsstrptr == lhsstrptr )
         rhs = lhs;
   }

   /* parse right hand side, if different from left hand side */
   if( rhsstrptr != NULL && rhsstrptr != lhsstrptr )
   {
      if( ! SCIPparseReal(scip, rhsstrptr, &rhs, &endptr) )
      {
         SCIPerrorMessage("error parsing right hand side number from <%s>\n", lhsstrptr);
         return SCIP_OKAY;
      }
   }

   /* initialize buffers for storing the variables and coefficients */
   coefssize = 100;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  coefssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, coefssize) );

   assert(varstrptr != NULL);

   /* parse linear sum to get variables and coefficients */
   SCIP_CALL( SCIPparseVarsLinearsum(scip, varstrptr, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );

   if( *success && requsize > coefssize )
   {
      /* realloc buffers and try again */
      coefssize = requsize;
      SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  coefssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, coefssize) );

      SCIP_CALL( SCIPparseVarsLinearsum(scip, varstrptr, vars, coefs, &nvars, coefssize, &requsize, &endptr, success) );
      assert(!*success || requsize <= coefssize); /* if successful, then should have had enough space now */
   }

   if( !*success )
   {
      SCIPerrorMessage("no luck in parsing linear sum '%s'\n", varstrptr);
   }
   else
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, nvars, vars, coefs, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      assert(vars != NULL);

      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/**! [Callback for the number of variables]*/
/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsLinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}
/**! [Callback for the number of variables]*/

/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecLinear)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);

   cons = eventdata->cons;
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* we can skip events droped for deleted constraints */
   if( SCIPconsIsDeleted(cons) )
      return SCIP_OKAY;

   eventtype = SCIPeventGetType(event);
   var = SCIPeventGetVar(event);

   if( (eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != 0 )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_Real val;
      int varpos;

      varpos = eventdata->varpos;
      assert(0 <= varpos && varpos < consdata->nvars);
      oldbound = SCIPeventGetOldbound(event);
      newbound = SCIPeventGetNewbound(event);
      assert(var != NULL);
      assert(consdata->vars[varpos] == var);
      val = consdata->vals[varpos];

      /* we only need to update the activities if the constraint is active,
       * otherwise we mark them to be invalid
       */
      if( SCIPconsIsActive(cons) )
      {
         /* update the activity values */
         if( (eventtype & SCIP_EVENTTYPE_LBCHANGED) != 0 )
            consdataUpdateActivitiesLb(scip, consdata, var, oldbound, newbound, val, TRUE);
         else
         {
            assert((eventtype & SCIP_EVENTTYPE_UBCHANGED) != 0);
            consdataUpdateActivitiesUb(scip, consdata, var, oldbound, newbound, val, TRUE);
         }
      }
      else
         consdataInvalidateActivities(consdata);

      consdata->presolved = FALSE;
      consdata->rangedrowpropagated = 0;

      /* bound change can turn the constraint infeasible or redundant only if it was a tightening */
      if( (eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED) != 0 )
      {
         SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );

         /* reset maximal activity delta, so that it will be recalculated on the next real propagation */
         if( consdata->maxactdeltavar == var )
         {
            consdata->maxactdelta = SCIP_INVALID;
            consdata->maxactdeltavar = NULL;
         }

         /* check whether bound tightening might now be successful */
         if( consdata->boundstightened > 0)
         {
            switch( eventtype )
            {
            case SCIP_EVENTTYPE_LBTIGHTENED:
               if( (val > 0.0 ? !SCIPisInfinity(scip, consdata->rhs) : !SCIPisInfinity(scip, -consdata->lhs)) )
                  consdata->boundstightened = 0;
               break;
            case SCIP_EVENTTYPE_UBTIGHTENED:
               if( (val > 0.0 ? !SCIPisInfinity(scip, -consdata->lhs) : !SCIPisInfinity(scip, consdata->rhs)) )
                  consdata->boundstightened = 0;
               break;
            default:
               SCIPerrorMessage("invalid event type %d\n", eventtype);
               return SCIP_INVALIDDATA;
            }
         }
      }
      /* update maximal activity delta if a bound was relaxed */
      else if( !SCIPisInfinity(scip, consdata->maxactdelta) )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real domain;
         SCIP_Real delta;

         assert((eventtype & SCIP_EVENTTYPE_BOUNDRELAXED) != 0);

         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);

         domain = ub - lb;
         delta = REALABS(val) * domain;

         if( delta > consdata->maxactdelta )
         {
            consdata->maxactdelta = delta;
            consdata->maxactdeltavar = var;
         }
      }
   }
   else if( (eventtype & SCIP_EVENTTYPE_VARFIXED) != 0 )
   {
      /* we want to remove the fixed variable */
      consdata->presolved = FALSE;
      consdata->removedfixings = FALSE;
      consdata->rangedrowpropagated = 0;

      /* reset maximal activity delta, so that it will be recalculated on the next real propagation */
      if( consdata->maxactdeltavar == var )
      {
         consdata->maxactdelta = SCIP_INVALID;
         consdata->maxactdeltavar = NULL;
      }
   }

   else if( (eventtype & SCIP_EVENTTYPE_VARUNLOCKED) != 0 )
   {
      /* there is only one lock left: we may multi-aggregate the variable as slack of an equation */
      assert(SCIPvarGetNLocksDown(var) <= 1);
      assert(SCIPvarGetNLocksUp(var) <= 1);
      consdata->presolved = FALSE;
   }
   else if( (eventtype & SCIP_EVENTTYPE_GBDCHANGED) != 0 )
   {
      SCIP_Real oldbound;
      SCIP_Real newbound;
      SCIP_Real val;
      int varpos;

      varpos = eventdata->varpos;
      assert(0 <= varpos && varpos < consdata->nvars);
      oldbound = SCIPeventGetOldbound(event);
      newbound = SCIPeventGetNewbound(event);
      assert(var != NULL);
      assert(consdata->vars[varpos] == var);
      val = consdata->vals[varpos];

      consdata->rangedrowpropagated = 0;

      /* update the activity values */
      if( (eventtype & SCIP_EVENTTYPE_GLBCHANGED) != 0 )
         consdataUpdateActivitiesGlbLb(scip, consdata, oldbound, newbound, val, TRUE);
      else
      {
         assert((eventtype & SCIP_EVENTTYPE_GUBCHANGED) != 0);
         consdataUpdateActivitiesGlbUb(scip, consdata, oldbound, newbound, val, TRUE);
      }
   }
   else
   {
      assert((eventtype & SCIP_EVENTTYPE_VARDELETED) != 0);
      consdata->varsdeleted = TRUE;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of conflict handler
 */

static
SCIP_DECL_CONFLICTEXEC(conflictExecLinear)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real lhs;
   int i;

   assert(scip != NULL);
   assert(conflicthdlr != NULL);
   assert(strcmp(SCIPconflicthdlrGetName(conflicthdlr), CONFLICTHDLR_NAME) == 0);
   assert(bdchginfos != NULL || nbdchginfos == 0);
   assert(result != NULL);

   /* don't process already resolved conflicts */
   if( resolved )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* create array of variables and coefficients: sum_{i \in P} x_i - sum_{i \in N} x_i >= 1 - |N| */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbdchginfos) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nbdchginfos) );
   lhs = 1.0;
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(bdchginfos != NULL);

      vars[i] = SCIPbdchginfoGetVar(bdchginfos[i]);

      /* we can only treat binary variables */
      /**@todo extend linear conflict constraints to some non-binary cases */
      if( !SCIPvarIsBinary(vars[i]) )
         break;

      /* check whether the variable is fixed to zero (P) or one (N) in the conflict set */
      if( SCIPbdchginfoGetNewbound(bdchginfos[i]) < 0.5 )
         vals[i] = 1.0;
      else
      {
         vals[i] = -1.0;
         lhs -= 1.0;
      }
   }

   if( i == nbdchginfos )
   {
      SCIP_CONS* cons;
      SCIP_CONS* upgdcons;
      char consname[SCIP_MAXSTRLEN];

      /* create a constraint out of the conflict set */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cf%" SCIP_LONGINT_FORMAT, SCIPgetNConflictConssApplied(scip));
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, nbdchginfos, vars, vals, lhs, SCIPinfinity(scip),
            FALSE, separate, FALSE, FALSE, TRUE, local, FALSE, dynamic, removable, FALSE) );

      /* try to automatically convert a linear constraint into a more specific and more specialized constraint */
      SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
      if( upgdcons != NULL )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         cons = upgdcons;
      }

      /* add conflict to SCIP */
      SCIP_CALL( SCIPaddConflict(scip, node, cons, validnode, conftype, cutoffinvolved) );

      *result = SCIP_CONSADDED;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/*
 * Quadratic constraint upgrading
 */


/** upgrades quadratic constraints with only and at least one linear variables into a linear constraint
 */
static
SCIP_DECL_QUADCONSUPGD(upgradeConsQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSDATA* upgdconsdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nupgdconss != NULL);
   assert(upgdconss  != NULL);

   *nupgdconss = 0;

   SCIPdebugMsg(scip, "upgradeConsQuadratic called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   if( SCIPgetNQuadVarTermsQuadratic(scip, cons) > 0 )
      return SCIP_OKAY;
   if( SCIPgetNLinearVarsQuadratic(scip, cons) == 0 )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* signal that we need more memory */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsLinear(scip, &upgdconss[0], SCIPconsGetName(cons),
         SCIPgetNLinearVarsQuadratic(scip, cons),
         SCIPgetLinearVarsQuadratic(scip, cons),
         SCIPgetCoefsLinearVarsQuadratic(scip, cons),
         SCIPgetLhsQuadratic(scip, cons), SCIPgetRhsQuadratic(scip, cons),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

   upgdconsdata = SCIPconsGetData(upgdconss[0]);
   assert(upgdconsdata != NULL);

   /* check violation of this linear constraint with absolute tolerances, to be consistent with the original quadratic constraint */
   upgdconsdata->checkabsolute = TRUE;

   SCIPdebugMsg(scip, "created linear constraint:\n");
   SCIPdebugPrintCons(scip, upgdconss[0], NULL);

   return SCIP_OKAY;
}

/** tries to upgrade a nonlinear constraint into a linear constraint */
static
SCIP_DECL_NONLINCONSUPGD(upgradeConsNonlinear)
{
   SCIP_CONSDATA* upgdconsdata;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   /* no interest in nonlinear constraints */
   if( SCIPgetExprgraphNodeNonlinear(scip, cons) != NULL )
      return SCIP_OKAY;

   /* no interest in constant constraints */
   if( SCIPgetNLinearVarsNonlinear(scip, cons) == 0 )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* request larger upgdconss array */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsLinear(scip, &upgdconss[0], SCIPconsGetName(cons),
         SCIPgetNLinearVarsNonlinear(scip, cons), SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
         SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );

   upgdconsdata = SCIPconsGetData(upgdconss[0]);
   assert(upgdconsdata != NULL);

   /* check violation of this linear constraint with absolute tolerances, to be consistent with the original nonlinear constraint */
   upgdconsdata->checkabsolute = TRUE;

   SCIPdebugMsg(scip, "created linear constraint:\n");
   SCIPdebugPrintCons(scip, upgdconss[0], NULL);

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrLinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_CONFLICTHDLR* conflicthdlr;

   assert(scip != NULL);

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecLinear, NULL) );

   /* create conflict handler for linear constraints */
   SCIP_CALL( SCIPincludeConflicthdlrBasic(scip, &conflicthdlr, CONFLICTHDLR_NAME, CONFLICTHDLR_DESC, CONFLICTHDLR_PRIORITY,
         conflictExecLinear, NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata, eventhdlr) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpLinear, consEnfopsLinear, consCheckLinear, consLockLinear,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyLinear, consCopyLinear) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveLinear) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteLinear) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsLinear) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitLinear) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreLinear) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolLinear) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeLinear) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsLinear) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsLinear) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitLinear) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpLinear) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseLinear) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolLinear, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintLinear) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropLinear, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropLinear) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpLinear, consSepasolLinear, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransLinear) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxLinear) );

   if( SCIPfindConshdlr(scip, "quadratic") != NULL )
   {
      /* include function that upgrades quadratic constraint to linear constraints */
      SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, upgradeConsQuadratic, QUADCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   if( SCIPfindConshdlr(scip, "nonlinear") != NULL )
   {
      /* include the linear constraint upgrade in the nonlinear constraint handler */
      SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, upgradeConsNonlinear, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   /* add linear constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/tightenboundsfreq",
         "multiplier on propagation frequency, how often the bounds are tightened (-1: never, 0: only at root)",
         &conshdlrdata->tightenboundsfreq, TRUE, DEFAULT_TIGHTENBOUNDSFREQ, -1, SCIP_MAXTREEDEPTH, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxrounds",
         "maximal number of separation rounds per node (-1: unlimited)",
         &conshdlrdata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxroundsroot",
         "maximal number of separation rounds per node in the root node (-1: unlimited)",
         &conshdlrdata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxsepacuts",
         "maximal number of cuts separated per separation round",
         &conshdlrdata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxsepacutsroot",
         "maximal number of cuts separated per separation round in the root node",
         &conshdlrdata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance", 
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/nmincomparisons",
         "number for minimal pairwise presolve comparisons",
         &conshdlrdata->nmincomparisons, TRUE, DEFAULT_NMINCOMPARISONS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/mingainpernmincomparisons",
         "minimal gain per minimal pairwise presolve comparisons to repeat pairwise comparison round",
         &conshdlrdata->mingainpernmincomp, TRUE, DEFAULT_MINGAINPERNMINCOMP, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/maxaggrnormscale",
         "maximal allowed relative gain in maximum norm for constraint aggregation (0.0: disable constraint aggregation)",
         &conshdlrdata->maxaggrnormscale, TRUE, DEFAULT_MAXAGGRNORMSCALE, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/maxeasyactivitydelta",
         "maximum activity delta to run easy propagation on linear constraint (faster, but numerically less stable)",
         &conshdlrdata->maxeasyactivitydelta, TRUE, DEFAULT_MAXEASYACTIVITYDELTA, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/maxcardbounddist",
         "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for separating knapsack cardinality cuts",
         &conshdlrdata->maxcardbounddist, TRUE, DEFAULT_MAXCARDBOUNDDIST, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/separateall",
         "should all constraints be subject to cardinality cut generation instead of only the ones with non-zero dual value?",
         &conshdlrdata->separateall, FALSE, DEFAULT_SEPARATEALL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/aggregatevariables",
         "should presolving search for aggregations in equations",
         &conshdlrdata->aggregatevariables, TRUE, DEFAULT_AGGREGATEVARIABLES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/simplifyinequalities",
         "should presolving try to simplify inequalities",
         &conshdlrdata->simplifyinequalities, TRUE, DEFAULT_SIMPLIFYINEQUALITIES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/dualpresolving",
         "should dual presolving steps be performed?",
         &conshdlrdata->dualpresolving, TRUE, DEFAULT_DUALPRESOLVING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/singletonstuffing",
         "should stuffing of singleton continuous variables be performed?",
         &conshdlrdata->singletonstuffing, TRUE, DEFAULT_SINGLETONSTUFFING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/singlevarstuffing",
         "should single variable stuffing be performed, which tries to fulfill constraints using the cheapest variable?",
         &conshdlrdata->singlevarstuffing, TRUE, DEFAULT_SINGLEVARSTUFFING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/sortvars", "apply binaries sorting in decr. order of coeff abs value?",
         &conshdlrdata->sortvars, TRUE, DEFAULT_SORTVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/checkrelmaxabs",
         "should the violation for a constraint with side 0.0 be checked relative to 1.0 (FALSE) or to the maximum absolute value in the activity (TRUE)?",
         &conshdlrdata->checkrelmaxabs, TRUE, DEFAULT_CHECKRELMAXABS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/detectcutoffbound",
         "should presolving try to detect constraints parallel to the objective function defining an upper bound and prevent these constraints from entering the LP?",
         &conshdlrdata->detectcutoffbound, TRUE, DEFAULT_DETECTCUTOFFBOUND, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/detectlowerbound",
         "should presolving try to detect constraints parallel to the objective function defining a lower bound and prevent these constraints from entering the LP?",
         &conshdlrdata->detectlowerbound, TRUE, DEFAULT_DETECTLOWERBOUND, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/detectpartialobjective",
         "should presolving try to detect subsets of constraints parallel to the objective function?",
         &conshdlrdata->detectpartialobjective, TRUE, DEFAULT_DETECTPARTIALOBJECTIVE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/rangedrowpropagation",
         "should presolving and propagation try to improve bounds, detect infeasibility, and extract sub-constraints from ranged rows and equations?",
         &conshdlrdata->rangedrowpropagation, TRUE, DEFAULT_RANGEDROWPROPAGATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/rangedrowartcons",
         "should presolving and propagation extract sub-constraints from ranged rows and equations?",
         &conshdlrdata->rangedrowartcons, TRUE, DEFAULT_RANGEDROWARTCONS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/rangedrowmaxdepth",
         "maximum depth to apply ranged row propagation",
         &conshdlrdata->rangedrowmaxdepth, TRUE, DEFAULT_RANGEDROWMAXDEPTH, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/rangedrowfreq",
         "frequency for applying ranged row propagation",
         &conshdlrdata->rangedrowfreq, TRUE, DEFAULT_RANGEDROWFREQ, 1, SCIP_MAXTREEDEPTH, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/multaggrremove",
         "should multi-aggregations only be performed if the constraint can be removed afterwards?",
         &conshdlrdata->multaggrremove, TRUE, DEFAULT_MULTAGGRREMOVE, NULL, NULL) );

   return SCIP_OKAY;
}

/** includes a linear constraint update method into the linear constraint handler */
SCIP_RETCODE SCIPincludeLinconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int                   priority,           /**< priority of upgrading method */
   const char*           conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_LINCONSUPGRADE* linconsupgrade;
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(linconsupgd != NULL);
   assert(conshdlrname != NULL );

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check if linear constraint update method already exists in constraint handler data */
   if( !conshdlrdataHasUpgrade(scip, conshdlrdata, linconsupgd, conshdlrname) )
   {
      /* create a linear constraint upgrade data object */
      SCIP_CALL( linconsupgradeCreate(scip, &linconsupgrade, linconsupgd, priority) );

      /* insert linear constraint update method into constraint handler data */
      SCIP_CALL( conshdlrdataIncludeUpgrade(scip, conshdlrdata, linconsupgrade) );

      /* adds parameter to turn on and off the upgrading step */
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/linear/upgrade/%s", conshdlrname);
      (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable linear upgrading for constraint handler <%s>", conshdlrname);
      SCIP_CALL( SCIPaddBoolParam(scip,
            paramname, paramdesc,
            &linconsupgrade->active, FALSE, TRUE, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** creates and captures a linear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
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
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   /* find the linear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("linear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* for the solving process we need linear rows, containing only active variables; therefore when creating a linear
    * constraint after presolving we have to ensure that it holds active variables
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITPRESOLVE && nvars > 0 )
   {
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_Real constant = 0.0;
      int nconsvars;
      int requiredsize;

      nconsvars = nvars;
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consvars, vars, nconsvars) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, &consvals, vals, nconsvars) );

      /* get active variables for new constraint */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, nconsvars, &constant, &requiredsize, TRUE) );

      /* if space was not enough we need to resize the buffers */
      if( requiredsize > nconsvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, requiredsize, &constant, &requiredsize, TRUE) );
         assert(requiredsize <= nconsvars);
      }

      /* adjust sides and check that we do not subtract infinity values */
      if( SCIPisInfinity(scip, REALABS(constant)) )
      {
         if( constant < 0.0 )
         {
            if( SCIPisInfinity(scip, lhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
            if( SCIPisInfinity(scip, rhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }

            lhs = -SCIPinfinity(scip);
            rhs = -SCIPinfinity(scip);
         }
         else
         {
            if( SCIPisInfinity(scip, -lhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }
            if( SCIPisInfinity(scip, -rhs) )
            {
               SCIPfreeBufferArray(scip, &consvals);
               SCIPfreeBufferArray(scip, &consvars);

               SCIPerrorMessage("try to generate inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", name);

               SCIPABORT();
               return SCIP_INVALIDDATA; /*lint !e527*/
            }

            lhs = SCIPinfinity(scip);
            rhs = SCIPinfinity(scip);
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, REALABS(lhs)) )
            lhs -= constant;
         if( !SCIPisInfinity(scip, REALABS(rhs)) )
            rhs -= constant;

         if( SCIPisInfinity(scip, -lhs) )
            lhs = -SCIPinfinity(scip);
         else if( SCIPisInfinity(scip, lhs) )
            lhs = SCIPinfinity(scip);

         if( SCIPisInfinity(scip, rhs) )
            rhs = SCIPinfinity(scip);
         else if( SCIPisInfinity(scip, -rhs) )
            rhs = -SCIPinfinity(scip);
      }

      /* create constraint data */
      SCIP_CALL( consdataCreate(scip, &consdata, nconsvars, consvars, consvals, lhs, rhs) );
      assert(consdata != NULL);

      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }
   else
   {
      /* create constraint data */
      SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars, vals, lhs, rhs) );
      assert(consdata != NULL);
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a linear constraint
 *  in its most basic version, i. e., all constraint flags are set to their basic value as explained for the
 *  method SCIPcreateConsLinear(); all flags can be set via SCIPsetConsFLAGNAME-methods in scip.h
 *
 *  @see SCIPcreateConsLinear() for information about the basic constraint flag configuration
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, nvars, vars, vals, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates by copying and captures a linear constraint */
SCIP_RETCODE SCIPcopyConsLinear(
   SCIP*                 scip,               /**< target SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to store the created target constraint */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in source variable array */
   SCIP_VAR**            sourcevars,         /**< source variables of the linear constraints */
   SCIP_Real*            sourcecoefs,        /**< coefficient array of the linear constraint, or NULL if all coefficients are one */
   SCIP_Real             lhs,                /**< left hand side of the linear constraint */
   SCIP_Real             rhs,                /**< right hand side of the linear constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store if the copying was valid */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* coefs;

   SCIP_Real constant;
   int requiredsize;
   int v;
   SCIP_Bool success;

   if( SCIPisGT(scip, lhs, rhs) )
   {
      *valid = FALSE;
      return SCIP_OKAY;
   }

   (*valid) = TRUE;

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, 0, NULL, NULL, lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      return SCIP_OKAY;
   }

   /* duplicate variable array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, sourcevars, nvars) );

   /* duplicate coefficient array */
   if( sourcecoefs != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &coefs, sourcecoefs, nvars) );
   }
   else 
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
      for( v = 0; v < nvars; ++v )
         coefs[v] = 1.0;
   }

   constant = 0.0;

   /* transform source variable to active variables of the source SCIP since only these can be mapped to variables of
    * the target SCIP
    */
   if( !SCIPvarIsOriginal(vars[0]) )
   {
      SCIP_CALL( SCIPgetProbvarLinearSum(sourcescip, vars, coefs, &nvars, nvars, &constant, &requiredsize, TRUE) );

      if( requiredsize > nvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(sourcescip, vars, coefs, &nvars, requiredsize, &constant, &requiredsize, TRUE) );
         assert(requiredsize <= nvars);
      }
   }
   else
   {
      for( v = 0; v < nvars; ++v )
      {
         assert(SCIPvarIsOriginal(vars[v]));
         SCIP_CALL( SCIPvarGetOrigvarSum(&vars[v], &coefs[v], &constant) );
         assert(vars[v] != NULL);
      }
   }


   success = TRUE;
   /* map variables of the source constraint to variables of the target SCIP */
   for( v = 0; v < nvars && success; ++v )
   {
      SCIP_VAR* var;
      var = vars[v];

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &vars[v], varmap, consmap, global, &success) );
      assert(!(success) || vars[v] != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( success )
   {
      if( !SCIPisInfinity(scip, -lhs) )
         lhs -= constant;

      if( !SCIPisInfinity(scip, rhs) )
         rhs -= constant;

      SCIP_CALL( SCIPcreateConsLinear(scip, cons, name, nvars, vars, coefs, lhs, rhs, 
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }
   else
      *valid = FALSE;

   /* free buffer array */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** adds coefficient to linear constraint (if it is not zero) */
SCIP_RETCODE SCIPaddCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< coefficient of constraint entry */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   /* for the solving process we need linear rows, containing only active variables; therefore when creating a linear
    * constraint after presolving we have to ensure that it holds active variables
    */
   if( SCIPgetStage(scip) >= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_Real constant = 0.0;
      SCIP_Real rhs;
      SCIP_Real lhs;
      int nconsvars;
      int requiredsize;
      int v;

      nconsvars = 1;
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
      consvars[0] = var;
      consvals[0] = val;

      /* get active variables for new constraint */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, nconsvars, &constant, &requiredsize, TRUE) );

      /* if space was not enough we need to resize the buffers */
      if( requiredsize > nconsvars )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvars, requiredsize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, requiredsize) );

         SCIP_CALL( SCIPgetProbvarLinearSum(scip, consvars, consvals, &nconsvars, requiredsize, &constant, &requiredsize, TRUE) );
         assert(requiredsize <= nconsvars);
      }

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      lhs = consdata->lhs;
      rhs = consdata->rhs;

      /* adjust sides and check that we do not subtract infinity values */
      /* constant is infinite */
      if( SCIPisInfinity(scip, REALABS(constant)) )
      {
	 if( constant < 0.0 )
	 {
	    if( SCIPisInfinity(scip, lhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }
	    if( SCIPisInfinity(scip, rhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }

	    lhs = -SCIPinfinity(scip);
	    rhs = -SCIPinfinity(scip);
	 }
	 else
	 {
	    if( SCIPisInfinity(scip, -lhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite left hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }
	    if( SCIPisInfinity(scip, -rhs) )
	    {
	       SCIPfreeBufferArray(scip, &consvals);
	       SCIPfreeBufferArray(scip, &consvars);

	       SCIPerrorMessage("adding variable <%s> leads to inconsistent constraint <%s>, active variables leads to a infinite constant constradict the infinite right hand side of the constraint\n", SCIPvarGetName(var), SCIPconsGetName(cons));

	       SCIPABORT();
	       return SCIP_INVALIDDATA; /*lint !e527*/
	    }

	    lhs = SCIPinfinity(scip);
	    rhs = SCIPinfinity(scip);
	 }
      }
      /* constant is not infinite */
      else
      {
	 if( !SCIPisInfinity(scip, REALABS(lhs)) )
	    lhs -= constant;
	 if( !SCIPisInfinity(scip, REALABS(rhs)) )
	    rhs -= constant;

	 if( SCIPisInfinity(scip, -lhs) )
	    lhs = -SCIPinfinity(scip);
	 else if( SCIPisInfinity(scip, lhs) )
	    lhs = SCIPinfinity(scip);

	 if( SCIPisInfinity(scip, rhs) )
	    rhs = SCIPinfinity(scip);
	 else if( SCIPisInfinity(scip, -rhs) )
	    rhs = -SCIPinfinity(scip);
      }

      /* add all active variables to constraint */
      for( v = nconsvars - 1; v >= 0; --v )
      {
	 SCIP_CALL( addCoef(scip, cons, consvars[v], consvals[v]) );
      }

      /* update left and right hand sides */
      SCIP_CALL( chgLhs(scip, cons, lhs));
      SCIP_CALL( chgRhs(scip, cons, rhs));

      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }
   else
   {
      SCIP_CALL( addCoef(scip, cons, var, val) );
   }

   return SCIP_OKAY;
}

/** changes coefficient of variable in linear constraint; deletes the variable if coefficient is zero; adds variable if
 *  not yet contained in the constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
SCIP_RETCODE SCIPchgCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             val                 /**< new coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Bool found;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPgetStage(scip) > SCIP_STAGE_PROBLEM || !SCIPconsIsOriginal(cons) || !SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("method may only be called during problem creation stage for original constraints and variables\n");
      return SCIP_INVALIDDATA;
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   found = FALSE;
   i = 0;
   while( i < consdata->nvars )
   {
      if( vars[i] == var )
      {
         if( found || SCIPisZero(scip, val) )
         {
            SCIP_CALL( delCoefPos(scip, cons, i) );

            /* decrease i by one since otherwise we would skip the coefficient which has been switched to position i */
            i--;
         }
         else
         {
            SCIP_CALL( chgCoefPos(scip, cons, i, val) );
         }
         found = TRUE;
      }
      i++;
   }

   if( !found && !SCIPisZero(scip, val) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, val) );
   }

   return SCIP_OKAY;
}

/** deletes variable from linear constraint
 *
 *  @note This method may only be called during problem creation stage for an original constraint and variable.
 *
 *  @note This method requires linear time to search for occurences of the variable in the constraint data.
 */
SCIP_RETCODE SCIPdelCoefLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPchgCoefLinear(scip, cons, var, 0.0) );

   return SCIP_OKAY;
}

/** gets left hand side of linear constraint */
SCIP_Real SCIPgetLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of linear constraint */
SCIP_Real SCIPgetRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** changes left hand side of linear constraint */
SCIP_RETCODE SCIPchgLhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgLhs(scip, cons, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of linear constraint */
SCIP_RETCODE SCIPchgRhsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( chgRhs(scip, cons, rhs) );

   return SCIP_OKAY;
}

/** gets the number of variables in the linear constraint */
int SCIPgetNVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return -1;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets the array of variables in the linear constraint; the user must not modify this array! */
SCIP_VAR** SCIPgetVarsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the array of coefficient values in the linear constraint; the user must not modify this array! */
SCIP_Real* SCIPgetValsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vals;
}

/** gets the activity of the linear constraint in the given solution
 *
 *  @note if the solution contains values at infinity, this method will return SCIP_INVALID in case the activity
 *        comprises positive and negative infinity contributions
 */
SCIP_Real SCIPgetActivityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIPgetRowSolActivity(scip, consdata->row, sol);
   else
      return consdataGetActivity(scip, consdata, sol);
}

/** gets the feasibility of the linear constraint in the given solution */
SCIP_Real SCIPgetFeasibilityLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_SOL*             sol                 /**< solution, or NULL to use current node's solution */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIPgetRowSolFeasibility(scip, consdata->row, sol);
   else
      return consdataGetFeasibility(scip, consdata, sol);
}

/** gets the dual solution of the linear constraint in the current LP */
SCIP_Real SCIPgetDualsolLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(!SCIPconsIsOriginal(cons)); /* original constraints would always return 0 */

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the linear constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(!SCIPconsIsOriginal(cons)); /* original constraints would always return 0 */

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return SCIP_INVALID;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given linear constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
SCIP_RETCODE SCIPupgradeConsLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_CONS**           upgdcons            /**< pointer to store upgraded constraint, or NULL if not successful */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real poscoeffsum;
   SCIP_Real negcoeffsum;
   SCIP_Bool integral;
   int nposbin;
   int nnegbin;
   int nposint;
   int nnegint;
   int nposimpl;
   int nnegimpl;
   int nposimplbin;
   int nnegimplbin;
   int nposcont;
   int nnegcont;
   int ncoeffspone;
   int ncoeffsnone;
   int ncoeffspint;
   int ncoeffsnint;
   int ncoeffspfrac;
   int ncoeffsnfrac;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(upgdcons != NULL);

   *upgdcons = NULL;

   /* we cannot upgrade a modifiable linear constraint, since we don't know what additional coefficients to expect */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   /* check for upgradability */
   if( SCIPconsGetNUpgradeLocks(cons) > 0 )
      return SCIP_OKAY;

   /* get the constraint handler and check, if it's really a linear constraint */
   conshdlr = SCIPconsGetHdlr(cons);
   if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not linear\n");
      return SCIP_INVALIDDATA;
   }

   /* get constraint handler data and constraint data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check, if the constraint was already upgraded and will be deleted anyway after preprocessing */
   if( consdata->upgraded )
      return SCIP_OKAY;

   /* check, if the constraint is already stored as LP row */
   if( consdata->row != NULL )
   {
      if( SCIProwIsInLP(consdata->row) )
      {
         SCIPerrorMessage("cannot upgrade linear constraint that is already stored as row in the LP\n");
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   /* normalize constraint */
   SCIP_CALL( normalizeCons(scip, cons) );


   /*
    * calculate some statistics on linear constraint
    */

   nposbin = 0;
   nnegbin = 0;
   nposint = 0;
   nnegint = 0;
   nposimpl = 0;
   nnegimpl = 0;
   nposimplbin = 0;
   nnegimplbin = 0;
   nposcont = 0;
   nnegcont = 0;
   ncoeffspone = 0;
   ncoeffsnone = 0;
   ncoeffspint = 0;
   ncoeffsnint = 0;
   ncoeffspfrac = 0;
   ncoeffsnfrac = 0;
   integral = TRUE;
   poscoeffsum = 0.0;
   negcoeffsum = 0.0;

   for( i = 0; i < consdata->nvars; ++i )
   {
      var = consdata->vars[i];
      val = consdata->vals[i];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      assert(!SCIPisZero(scip, val));

      switch( SCIPvarGetType(var) )
      {
      case SCIP_VARTYPE_BINARY:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposbin++;
         else
            nnegbin++;
         break;
      case SCIP_VARTYPE_INTEGER:
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposint++;
         else
            nnegint++;
         break;
      case SCIP_VARTYPE_IMPLINT:
         if( SCIPvarIsBinary(var) )
         {
            if( val >= 0.0 )
               nposimplbin++;
            else
               nnegimplbin++;
         }
         if( !SCIPisZero(scip, lb) || !SCIPisZero(scip, ub) )
            integral = integral && SCIPisIntegral(scip, val);
         if( val >= 0.0 )
            nposimpl++;
         else
            nnegimpl++;
         break;
      case SCIP_VARTYPE_CONTINUOUS:
         integral = integral && SCIPisEQ(scip, lb, ub) && SCIPisIntegral(scip, val * lb);
         if( val >= 0.0 )
            nposcont++;
         else
            nnegcont++;
         break;
      default:
         SCIPerrorMessage("unknown variable type\n");
         return SCIP_INVALIDDATA;
      }
      if( SCIPisEQ(scip, val, 1.0) )
         ncoeffspone++;
      else if( SCIPisEQ(scip, val, -1.0) )
         ncoeffsnone++;
      else if( SCIPisIntegral(scip, val) )
      {
         if( SCIPisPositive(scip, val) )
            ncoeffspint++;
         else
            ncoeffsnint++;
      }
      else
      {
         if( SCIPisPositive(scip, val) )
            ncoeffspfrac++;
         else
            ncoeffsnfrac++;
      }
      if( SCIPisPositive(scip, val) )
         poscoeffsum += val;
      else
         negcoeffsum += val;
   }


   /*
    * call the upgrading methods
    */

   SCIPdebugMsg(scip, "upgrading linear constraint <%s> (%d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nlinconsupgrades);
   SCIPdebugMsg(scip, " +bin=%d -bin=%d +int=%d -int=%d +impl=%d -impl=%d +cont=%d -cont=%d +1=%d -1=%d +I=%d -I=%d +F=%d -F=%d possum=%.15g negsum=%.15g integral=%u\n",
      nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposcont, nnegcont,
      ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
      poscoeffsum, negcoeffsum, integral);

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nlinconsupgrades && *upgdcons == NULL; ++i )
   {
      if( conshdlrdata->linconsupgrades[i]->active )
      {
         SCIP_CALL( conshdlrdata->linconsupgrades[i]->linconsupgd(scip, cons, consdata->nvars,
               consdata->vars, consdata->vals, consdata->lhs, consdata->rhs,
               nposbin, nnegbin, nposint, nnegint, nposimpl, nnegimpl, nposimplbin, nnegimplbin, nposcont, nnegcont,
               ncoeffspone, ncoeffsnone, ncoeffspint, ncoeffsnint, ncoeffspfrac, ncoeffsnfrac,
               poscoeffsum, negcoeffsum, integral,
               upgdcons) );
      }
   }

#ifdef SCIP_DEBUG
   if( *upgdcons != NULL )
   {
      SCIPdebugPrintCons(scip, cons, NULL);
      SCIPdebugMsg(scip, " -> upgraded to constraint type <%s>\n", SCIPconshdlrGetName(SCIPconsGetHdlr(*upgdcons)));
      SCIPdebugPrintCons(scip, *upgdcons, NULL);
   }
#endif

   return SCIP_OKAY;
}
