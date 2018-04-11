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

/**@file    prop_obbt.c
 * @ingroup PROPAGATORS
 * @brief   optimization-based bound tightening propagator
 * @author  Stefan Weltge
 * @author  Benjamin Mueller
 */

/**@todo if bound tightenings of other propagators are the reason for lpsolstat != SCIP_LPSOLSTAT_OPTIMAL, resolve LP */
/**@todo only run more than once in root node if primal bound improved or many cuts were added to the LP */
/**@todo filter bounds of a variable already if SCIPisLbBetter()/SCIPisUbBetter() would return FALSE */
/**@todo improve warmstarting of LP solving */
/**@todo include bound value (finite/infinite) into getScore() function */
/**@todo use unbounded ray in filtering */
/**@todo do we want to run if the LP is unbounded, maybe for infinite variable bounds? */
/**@todo add first filter round in direction of objective function */
/**@todo implement conflict resolving callback by calling public method of genvbounds propagator, since the reason are
 *       exactly the variable bounds with nonnegative reduced costs stored in the right-hand side of the generated
 *       generalized variable bound (however, this only makes sense if we run locally)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_obbt.h"
#include "scip/prop_genvbounds.h"
#include "scip/debug.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_abspower.h"
#include "scip/cons_bivariate.h"

#define PROP_NAME                       "obbt"
#define PROP_DESC                       "optimization-based bound tightening propagator"
#define PROP_TIMING                     SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY                -1000000      /**< propagator priority */
#define PROP_FREQ                           0      /**< propagator frequency */
#define PROP_DELAY                       TRUE      /**< should propagation method be delayed, if other propagators
                                                    *   found reductions? */

#define DEFAULT_CREATE_GENVBOUNDS        TRUE      /**< should obbt try to provide genvbounds if possible? */
#define DEFAULT_FILTERING_NORM           TRUE      /**< should coefficients in filtering be normalized w.r.t. the
                                                    *   domains sizes? */
#define DEFAULT_APPLY_FILTERROUNDS      FALSE      /**< try to filter bounds in so-called filter rounds by solving
                                                    *   auxiliary LPs? */
#define DEFAULT_APPLY_TRIVIALFITLERING   TRUE      /**< should obbt try to use the LP solution to filter some bounds? */
#define DEFAULT_GENVBDSDURINGFILTER      TRUE      /**< try to genrate genvbounds during trivial and aggressive filtering? */
#define DEFAULT_DUALFEASTOL              1e-9      /**< feasibility tolerance for reduced costs used in obbt; this value
                                                    *   is used if SCIP's dual feastol is greater */
#define DEFAULT_CONDITIONLIMIT           -1.0      /**< maximum condition limit used in LP solver (-1.0: no limit) */
#define DEFAULT_BOUNDSTREPS             0.001      /**< minimal relative improve for strengthening bounds */
#define DEFAULT_FILTERING_MIN               2      /**< minimal number of filtered bounds to apply another filter
                                                    *   round */
#define DEFAULT_ITLIMITFACTOR            10.0      /**< multiple of root node LP iterations used as total LP iteration
                                                    *   limit for obbt (<= 0: no limit ) */
#define DEFAULT_MINITLIMIT              5000L      /**< minimum LP iteration limit */
#define DEFAULT_ONLYNONCONVEXVARS       FALSE      /**< only apply obbt on non-convex variables */
#define DEFAULT_TIGHTINTBOUNDSPROBING    TRUE      /**< should bounds of integral variables be tightened during
                                                    *   the probing mode? */
#define DEFAULT_TIGHTCONTBOUNDSPROBING  FALSE      /**< should bounds of continuous variables be tightened during
                                                    *   the probing mode? */
#define DEFAULT_ORDERINGALGO                1      /**< which type of ordering algorithm should we use?
                                                    *   (0: no, 1: greedy, 2: greedy reverse) */
#define OBBT_SCOREBASE                      5      /**< base that is used to calculate a bounds score value */
#define GENVBOUND_PROP_NAME             "genvbounds"
#define INTERVALINFTY                   1E+43      /**< value for infinity in interval operations */

#define DEFAULT_SEPARATESOL             FALSE      /**< should the obbt LP solution be separated? note that that by
                                                    *   separating solution OBBT will apply all bound tightenings
                                                    *   immediatly */
#define DEFAULT_SEPAMINITER                 0      /**< minimum number of iteration spend to separate an obbt LP solution */
#define DEFAULT_SEPAMAXITER                10      /**< maximum number of iteration spend to separate an obbt LP solution */
#define DEFAULT_GENVBDSDURINGSEPA        TRUE      /**< try to create genvbounds during separation process? */
#define DEFAULT_PROPAGATEFREQ               0      /**< trigger a propagation round after that many bound tightenings
                                                    *   (0: no propagation) */
#define DEFAULT_CREATE_BILININEQS        TRUE      /**< solve auxiliary LPs in order to find valid inequalities for bilinear terms? */
#define DEFAULT_ITLIMITFAC_BILININEQS     3.0      /**< multiple of OBBT LP limit used as total LP iteration limit for solving bilinear inequality LPs (< 0 for no limit) */
#define DEFAULT_MINNONCONVEXITY          1e-1      /**< minimum nonconvexity for choosing a bilinear term */
#define DEFAULT_RANDSEED                  149      /**< initial random seed */


/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/*
 * Data structures
 */

/** bound data */
struct Bound
{
   SCIP_VAR*             var;                /**< variable */
   SCIP_Real             newval;             /**< stores a probably tighter value for this bound */
   SCIP_BOUNDTYPE        boundtype;          /**< type of bound */
   unsigned int          score;              /**< score value that is used to group bounds */
   unsigned int          filtered:1;         /**< thrown out during pre-filtering step */
   unsigned int          found:1;            /**< stores whether a probably tighter value for this bound was found */
   unsigned int          done:1;             /**< has this bound been processed already? */
   unsigned int          nonconvex:1;        /**< is this bound affecting a nonconvex term? */
   int                   index;              /**< unique index */
};
typedef struct Bound BOUND;

/* all possible corners of a rectangular domain */
enum Corner
{
   LEFTBOTTOM  = 1,
   RIGHTBOTTOM = 2,
   RIGHTTOP    = 4,
   LEFTTOP     = 8,
   FILTERED    = 15
};
typedef enum Corner CORNER;

/** bilinear bound data */
struct BilinBound
{
   SCIP_VAR*             x;                  /**< first variable */
   SCIP_VAR*             y;                  /**< second variable */
   int                   filtered;           /**< corners that could be thrown out during pre-filtering step */
   unsigned int          done:1;             /**< has this bilinear term been processed already? */
   int                   nunderest;          /**< number of constraints that require to underestimate the bilinear term */
   int                   noverest;           /**< number of constraints that require to overestimate the bilinear term */
   int                   index;              /**< index of the bilinear term in the quadratic constraint handler */
   SCIP_Real             score;              /**< score value that is used to group bilinear term bounds */
};
typedef struct BilinBound BILINBOUND;

/** propagator data */
struct SCIP_PropData
{
   BOUND**               bounds;             /**< array of interesting bounds */
   BILINBOUND**          bilinbounds;        /**< array of interesting bilinear bounds */
   SCIP_ROW*             cutoffrow;          /**< pointer to current objective cutoff row */
   SCIP_PROP*            genvboundprop;      /**< pointer to genvbound propagator */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Longint          lastnode;           /**< number of last node where obbt was performed */
   SCIP_Longint          npropagatedomreds;  /**< number of domain reductions found during propagation */
   SCIP_Longint          nprobingiterations; /**< number of LP iterations during the probing mode */
   SCIP_Longint          nfilterlpiters;     /**< number of LP iterations spend for filtering */
   SCIP_Longint          minitlimit;         /**< minimum LP iteration limit */
   SCIP_Longint          itlimitbilin;       /**< total LP iterations limit for solving bilinear inequality LPs */
   SCIP_Longint          itusedbilin;        /**< total LP iterations used for solving bilinear inequality LPs */
   SCIP_Real             dualfeastol;        /**< feasibility tolerance for reduced costs used in obbt; this value is
                                              *   used if SCIP's dual feastol is greater */
   SCIP_Real             conditionlimit;     /**< maximum condition limit used in LP solver (-1.0: no limit) */
   SCIP_Real             boundstreps;        /**< minimal relative improve for strengthening bounds */
   SCIP_Real             itlimitfactor;      /**< LP iteration limit for obbt will be this factor times total LP
                                              *   iterations in root node */
   SCIP_Real             itlimitfactorbilin; /**< multiple of OBBT LP limit used as total LP iteration limit for solving bilinear inequality LPs (< 0 for no limit) */
   SCIP_Real             minnonconvexity;    /**< lower bound on minimum absolute value of nonconvex eigenvalues for a bilinear term */
   SCIP_Bool             applyfilterrounds;  /**< apply filter rounds? */
   SCIP_Bool             applytrivialfilter; /**< should obbt try to use the LP solution to filter some bounds? */
   SCIP_Bool             genvbdsduringfilter;/**< should we try to generate genvbounds during trivial and aggressive
                                              *   filtering? */
   SCIP_Bool             genvbdsduringsepa;  /**< try to create genvbounds during separation process? */
   SCIP_Bool             creategenvbounds;   /**< should obbt try to provide genvbounds if possible? */
   SCIP_Bool             normalize;          /**< should coefficients in filtering be normalized w.r.t. the domains
                                              *   sizes? */
   SCIP_Bool             onlynonconvexvars;  /**< only apply obbt on non-convex variables */
   SCIP_Bool             tightintboundsprobing; /**< should bounds of integral variables be tightened during
                                              *   the probing mode? */
   SCIP_Bool             tightcontboundsprobing;/**< should bounds of continuous variables be tightened during
                                              *   the probing mode? */
   SCIP_Bool             separatesol;        /**< should the obbt LP solution be separated? note that that by
                                              *   separating solution OBBT will apply all bound tightenings
                                              *   immediatly */
   SCIP_Bool             createbilinineqs;   /**< solve auxiliary LPs in order to find valid inequalities for bilinear terms? */
   int                   orderingalgo;       /**< which type of ordering algorithm should we use?
                                              *   (0: no, 1: greedy, 2: greedy reverse) */
   int                   nbounds;            /**< length of interesting bounds array */
   int                   nbilinbounds;       /**< length of interesting bilinear bounds array */
   int                   boundssize;         /**< size of bounds array */
   int                   nminfilter;         /**< minimal number of filtered bounds to apply another filter round */
   int                   nfiltered;          /**< number of filtered bounds by solving auxiliary variables */
   int                   ntrivialfiltered;   /**< number of filtered bounds because the LP value was equal to the bound */
   int                   nsolvedbounds;      /**< number of solved bounds during the loop in applyObbt() */
   int                   ngenvboundsprobing; /**< number of non-trivial genvbounds generated and added during obbt */
   int                   ngenvboundsaggrfil; /**< number of non-trivial genvbounds found during aggressive filtering */
   int                   ngenvboundstrivfil; /**< number of non-trivial genvbounds found during trivial filtering */
   int                   lastidx;            /**< index to store the last undone and unfiltered bound */
   int                   lastbilinidx;       /**< index to store the last undone and unfiltered bilinear bound */
   int                   sepaminiter;        /**< minimum number of iteration spend to separate an obbt LP solution */
   int                   sepamaxiter;        /**< maximum number of iteration spend to separate an obbt LP solution */
   int                   propagatefreq;      /**< trigger a propagation round after that many bound tightenings
                                              *   (0: no propagation) */
   int                   propagatecounter;   /**< number of bound tightenings since the last propagation round */
};


/*
 * Local methods
 */

/** solves the LP and handles errors */
static
SCIP_RETCODE solveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlimit,            /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            error,              /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            optimal             /**< was the LP solved to optimalilty? */
   )
{
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(itlimit == -1 || itlimit >= 0);
   assert(error != NULL);
   assert(optimal != NULL);

   *optimal = FALSE;
   *error = FALSE;

   retcode = SCIPsolveProbingLP(scip, itlimit, error, NULL);

   lpsolstat = SCIPgetLPSolstat(scip);

   /* an error should not kill the overall solving process */
   if( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "   error while solving LP in obbt propagator; LP solve terminated with code <%d>\n", retcode);
      SCIPwarningMessage(scip, "   this does not affect the remaining solution procedure --> continue\n");

      *error = TRUE;

      return SCIP_OKAY;
   }

   if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      assert(!*error);
      *optimal = TRUE;
   }
#ifdef SCIP_DEBUG
   else
   {
      switch( lpsolstat )
      {
      case SCIP_LPSOLSTAT_ITERLIMIT:
         SCIPdebugMsg(scip, "   reached lp iteration limit\n");
         break;
      case SCIP_LPSOLSTAT_TIMELIMIT:
         SCIPdebugMsg(scip, "   reached time limit while solving lp\n");
         break;
      case SCIP_LPSOLSTAT_UNBOUNDEDRAY:
         SCIPdebugMsg(scip, "   lp was unbounded\n");
         break;
      case SCIP_LPSOLSTAT_NOTSOLVED:
         SCIPdebugMsg(scip, "   lp was not solved\n");
         break;
      case SCIP_LPSOLSTAT_ERROR:
         SCIPdebugMsg(scip, "   an error occured during solving lp\n");
         break;
      case SCIP_LPSOLSTAT_INFEASIBLE:
      case SCIP_LPSOLSTAT_OBJLIMIT:
      case SCIP_LPSOLSTAT_OPTIMAL: /* should not appear because it is handled earlier */
      default:
         SCIPdebugMsg(scip, "   received an unexpected solstat during solving lp: %d\n", lpsolstat);
      }
   }
#endif

   return SCIP_OKAY;
}

/** adds the objective cutoff to the LP; must be in probing mode */
static
SCIP_RETCODE addObjCutoff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   SCIP_ROW* row;
   SCIP_VAR** vars;
   char rowname[SCIP_MAXSTRLEN];

   int nvars;
   int i;

   assert(scip != NULL);
   assert(SCIPinProbing(scip));
   assert(propdata != NULL);
   assert(propdata->cutoffrow == NULL);

   if( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
   {
      SCIPdebugMsg(scip, "no objective cutoff since there is no cutoff bound\n");
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "create objective cutoff and add it to the LP\n");

   /* get variables data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create objective cutoff row; set local flag to FALSE since primal cutoff is globally valid */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "obbt_objcutoff");
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, rowname, -SCIPinfinity(scip), SCIPgetCutoffbound(scip), FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], SCIPvarGetObj(vars[i])) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   /* add row to the LP */
   SCIP_CALL( SCIPaddRowProbing(scip, row) );

   propdata->cutoffrow = row;
   assert(SCIProwIsInLP(propdata->cutoffrow));

   return SCIP_OKAY;
}

/** determines, whether a variable is already locally fixed */
static
SCIP_Bool varIsFixedLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to check */
   )
{
   return SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
}

/** sets objective to minimize or maximize a single variable */
static
SCIP_RETCODE setObjProbing(
   SCIP*                 scip,
   SCIP_PROPDATA*        propdata,
   BOUND*                bound,
   SCIP_Real             coef
   )
{
#ifdef SCIP_DEBUG
   SCIP_VAR** vars;
   int nvars;
   int counter;
   int i;
#endif

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( bound != NULL );

   /* set the objective for bound->var */
   if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, coef) );
   }
   else
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, -coef) );
   }

#ifdef SCIP_DEBUG
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   counter = 0;

   for( i = 0; i < nvars; ++i )
   {
      if( SCIPgetVarObjProbing(scip, vars[i]) != 0.0 )
         ++counter;
   }

   assert((counter == 0 && coef == 0.0) || (counter == 1 && coef != 0.0));
#endif

   return SCIP_OKAY;
}

/** determines whether variable should be included in the right-hand side of the generalized variable bound */
static
SCIP_Bool includeVarGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to check */
   )
{
   SCIP_Real redcost;

   assert(scip != NULL);
   assert(var != NULL);

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return FALSE;

   redcost = SCIPgetVarRedcost(scip, var);
   assert(redcost != SCIP_INVALID); /*lint !e777 */

   if( redcost == SCIP_INVALID ) /*lint !e777 */
      return FALSE;

   if( redcost < SCIPdualfeastol(scip) && redcost > -SCIPdualfeastol(scip) )
      return FALSE;

   return TRUE;
}

/** returns number of LP iterations left (-1: no limit ) */
static
int getIterationsLeft(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint          nolditerations,     /**< iterations count at the beginning of the corresponding function */
   SCIP_Longint          itlimit             /**< LP iteration limit (-1: no limit) */
   )
{
   SCIP_Longint itsleft;

   assert(scip != NULL);
   assert(nolditerations >= 0);
   assert(itlimit == -1 || itlimit >= 0);

   if( itlimit == -1 )
   {
      SCIPdebugMsg(scip, "iterations left: unlimited\n");
      return -1;
   }
   else
   {
      itsleft = itlimit - ( SCIPgetNLPIterations(scip) - nolditerations );
      itsleft = MAX(itsleft, 0);
      itsleft = MIN(itsleft, INT_MAX);

      SCIPdebugMsg(scip, "iterations left: %d\n", (int) itsleft);
      return (int) itsleft;
   }
}

/** returns the objective coefficient for a variable's bound that will be chosen during filtering */
static
SCIP_Real getFilterCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_VAR*             var,                /**< variable */
   SCIP_BOUNDTYPE        boundtype           /**< boundtype to be filtered? */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(var != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   /* this function should not be called for fixed variables */
   assert(!varIsFixedLocal(scip, var));

   /* infinite bounds will not be reached */
   if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisInfinity(scip, -lb) )
      return 0.0;
   if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisInfinity(scip, ub) )
      return 0.0;

   if( propdata->normalize )
   {
      /* if the length of the domain is too large then the coefficient should be set to +/- 1.0 */
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisInfinity(scip, ub) )
         return 1.0;
      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisInfinity(scip, -lb) )
         return -1.0;

      /* otherwise the coefficient is +/- 1.0 / ( ub - lb ) */
      return boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 / (ub - lb) : -1.0 / (ub - lb);
   }
   else
   {
      return boundtype == SCIP_BOUNDTYPE_LOWER ? 1.0 : -1.0;
   }
}

/** creates a genvbound if the dual LP solution provides such information
 *
 *  Consider the problem
 *
 *     min { +/- x_i : obj * x <= z, lb <= Ax <= ub, l <= x <= u },
 *
 *  where z is the current cutoff bound. Let (mu, nu, gamma, alpha, beta) >= 0 be the optimal solution of the dual of
 *  problem (P), where the variables correspond to the primal inequalities in the following way:
 *
 *           Ax >=  lb    <->   mu
 *          -Ax >= -ub    <->   nu
 *     -obj * x >=  -z    <->   gamma
 *            x >=   l    <->   alpha
 *           -x >=  -u    <->   beta
 *
 *  Fixing these multipliers, by weak duality, we obtain the inequality
 *
 *     +/- x_i >= lb*mu - ub*nu - z*gamma + l*alpha - u*beta
 *
 *  that holds for all primal feasible points x with objective value at least z. Setting
 *
 *     c = lb*mu - ub*nu, redcost_k = alpha_k - beta_k
 *
 *  we obtain the inequality
 *
 *     +/- x_i >= sum ( redcost_k * x_k ) + (-gamma) * cutoff_bound + c,
 *
 *  that holds for all primal feasible points with objective value at least cutoff_bound. Therefore, the latter
 *  inequality can be added as a generalized variable bound.
 */
static
SCIP_RETCODE createGenVBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   BOUND*                bound,              /**< bound of x_i */
   SCIP_Bool*            found               /**< pointer to store if we have found a non-trivial genvbound */
   )
{
   assert(scip != NULL);
   assert(bound != NULL);
   assert(propdata != NULL);
   assert(propdata->genvboundprop != NULL);
   assert(found != NULL);

   *found = FALSE;

   /* make sure we are in probing mode having an optimal LP solution */
   assert(SCIPinProbing(scip));

   assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

   /* only genvbounds created in the root node are globally valid
    *
    * note: depth changes to one if we use the probing mode to solve the obbt LPs
    */
   assert(SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1));

   SCIPdebugMsg(scip, "      try to create a genvbound for <%s>...\n", SCIPvarGetName(bound->var));

   /* a genvbound with a multiplier for x_i would not help us */
   if( SCIPisZero(scip, SCIPgetVarRedcost(scip, bound->var)) )
   {
      SCIP_VAR** vars;                          /* global variables array */
      SCIP_VAR** genvboundvars;                 /* genvbound variables array */

      SCIP_VAR* xi;                             /* variable x_i */

      SCIP_Real* genvboundcoefs;                /* genvbound coefficients array */

      SCIP_Real gamma_dual;                     /* dual multiplier of objective cutoff */

      int k;                                    /* variable for indexing global variables array */
      int ncoefs;                               /* number of nonzero coefficients in genvbound */
      int nvars;                                /* number of global variables */

      /* set x_i */
      xi = bound->var;

      /* get variable data */
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

      /* count nonzero coefficients in genvbound */
      ncoefs = 0;
      for( k = 0; k < nvars; k++ )
      {
         if( includeVarGenVBound(scip, vars[k]) )
         {
            assert(vars[k] != xi);
            ncoefs++;
         }
      }

      /* get dual multiplier for the objective cutoff (set to zero if there is no) */
      if( propdata->cutoffrow == NULL )
      {
         gamma_dual = 0.0;
      }
      else
      {
         assert(!SCIPisInfinity(scip, SCIPgetCutoffbound(scip)));

         /* note that the objective cutoff is of the form
          *    -inf <= obj * x <= cutoff_bound
          * but we want the positive dual multiplier!
          */
         gamma_dual = -SCIProwGetDualsol(propdata->cutoffrow);
      }

      /* we need at least one nonzero coefficient or a nonzero dual multiplier for the objective cutoff */
      if( ncoefs > 0 || !SCIPisZero(scip, gamma_dual) )
      {
         SCIP_Bool addgenvbound;                /* if everything is fine with the redcosts and the bounds, add the genvbound */
         SCIP_Real c;                           /* helper variable to calculate constant term in genvbound */
         int idx;                               /* variable for indexing genvbound's coefficients array */

         /* add the bound if the bool is still TRUE after the loop */
         addgenvbound = TRUE;

         /* there should be no coefficient for x_i */
         assert(SCIPisZero(scip, SCIPgetVarRedcost(scip, xi)));

         /* allocate memory for storing the genvbounds right-hand side variables and coefficients */
         SCIP_CALL( SCIPallocBufferArray(scip, &(genvboundvars), ncoefs) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(genvboundcoefs), ncoefs) );

         /* set c = lb*mu - ub*nu - z*gamma + l*alpha - u*beta */
         c = SCIPgetLPObjval(scip);

         /* subtract ( - z * gamma ) from c */
         c += SCIPgetCutoffbound(scip) * gamma_dual;

         /* subtract ( l*alpha - u*beta ) from c and set the coefficients of the variables */
         idx = 0;
         for( k = 0; k < nvars; k++ )
         {
            SCIP_VAR* xk;

            xk = vars[k];

            if( includeVarGenVBound(scip, xk) )
            {
               SCIP_Real redcost;

               redcost = SCIPgetVarRedcost(scip, xk);

               assert(redcost != SCIP_INVALID); /*lint !e777 */
               assert(xk != xi);

               /* in this case dont add a genvbound */
               if( ( (redcost > SCIPdualfeastol(scip))  && SCIPisInfinity(scip, -SCIPvarGetLbLocal(xk)) ) ||
                  ( (redcost < -SCIPdualfeastol(scip))  && SCIPisInfinity(scip, SCIPvarGetUbLocal(xk)) ) )
               {
                  addgenvbound = FALSE;
                  break;
               }

               /* store coefficients */
               assert(idx < ncoefs);
               genvboundvars[idx] = xk;
               genvboundcoefs[idx] = redcost;
               idx++;

               /* if redcost > 0, then redcost = alpha_k, otherwise redcost = - beta_k */
               assert(redcost <= 0 || !SCIPisInfinity(scip, -SCIPvarGetLbLocal(xk)));
               assert(redcost >= 0 || !SCIPisInfinity(scip, SCIPvarGetUbLocal(xk)));
               c -= redcost > 0 ? redcost * SCIPvarGetLbLocal(xk) : redcost * SCIPvarGetUbLocal(xk);
            }
         }

         assert(!addgenvbound || idx == ncoefs);

         /* add genvbound */
         if( addgenvbound && !SCIPisInfinity(scip, -c) )
         {
            SCIPdebugMsg(scip, "         adding genvbound\n");
            SCIP_CALL( SCIPgenVBoundAdd(scip, propdata->genvboundprop, genvboundvars, xi, genvboundcoefs, ncoefs,
                  !SCIPisPositive(scip, gamma_dual) ? 0.0 : -gamma_dual, c, bound->boundtype) );

            *found = TRUE;
         }

         /* free arrays */
         SCIPfreeBufferArray(scip, &genvboundcoefs);
         SCIPfreeBufferArray(scip, &genvboundvars);
      }
      else
      {
         SCIPdebugMsg(scip, "         trivial genvbound, skipping\n");
      }
   }
   else
   {
      SCIPdebugMsg(scip, "         found multiplier for <%s>: %g, skipping\n",
         SCIPvarGetName(bound->var), SCIPgetVarRedcost(scip, bound->var));
   }

   return SCIP_OKAY;
}

/** exchange a bound which has been processed and updates the last undone and unfiltered bound index
 *  NOTE: this method has to be called after filtering or processing a bound
 */
static
void exchangeBounds(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   i                   /**< bound that was filtered or processed */
   )
{
   assert(i >= 0 && i < propdata->nbounds);
   assert(propdata->lastidx >= 0 && propdata->lastidx < propdata->nbounds);

   /* exchange the bounds */
   if( propdata->lastidx != i )
   {
      BOUND* tmp;

      tmp = propdata->bounds[i];
      propdata->bounds[i] = propdata->bounds[propdata->lastidx];
      propdata->bounds[propdata->lastidx] = tmp;
   }

   propdata->lastidx -= 1;
}

/** helper function to return a corner of the domain of two variables */
static
void getCorner(
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   CORNER                corner,             /**< corner */
   SCIP_Real*            px,                 /**< buffer to store point for x */
   SCIP_Real*            py                  /**< buffer to store point for y */
   )
{
   assert(x != NULL);
   assert(y != NULL);
   assert(px != NULL);
   assert(py != NULL);

   switch( corner )
   {
      case LEFTBOTTOM:
         *px = SCIPvarGetLbGlobal(x);
         *py = SCIPvarGetLbGlobal(y);
         break;
      case RIGHTBOTTOM:
         *px = SCIPvarGetUbGlobal(x);
         *py = SCIPvarGetLbGlobal(y);
         break;
      case LEFTTOP:
         *px = SCIPvarGetLbGlobal(x);
         *py = SCIPvarGetUbGlobal(y);
         break;
      case RIGHTTOP:
         *px = SCIPvarGetUbGlobal(x);
         *py = SCIPvarGetUbGlobal(y);
         break;
      case FILTERED:
         SCIPABORT();
   }
}

/** helper function to return the two end points of a diagonal */
static
void getCorners(
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   CORNER                corner,             /**< corner */
   SCIP_Real*            xs,                 /**< buffer to store start point for x */
   SCIP_Real*            ys,                 /**< buffer to store start point for y */
   SCIP_Real*            xt,                 /**< buffer to store end point for x */
   SCIP_Real*            yt                  /**< buffer to store end point for y */
   )
{
   assert(x != NULL);
   assert(y != NULL);
   assert(xs != NULL);
   assert(ys != NULL);
   assert(xt != NULL);
   assert(yt != NULL);

   /* get end point */
   getCorner(x,y, corner, xt, yt);

   /* get start point */
   switch( corner )
   {
      case LEFTBOTTOM:
         getCorner(x,y, RIGHTTOP, xs, ys);
         break;
      case RIGHTBOTTOM:
         getCorner(x,y, LEFTTOP, xs, ys);
         break;
      case LEFTTOP:
         getCorner(x,y, RIGHTBOTTOM, xs, ys);
         break;
      case RIGHTTOP:
         getCorner(x,y, LEFTBOTTOM, xs, ys);
         break;
      case FILTERED:
         SCIPABORT();
   }
}

/** trying to filter some bounds using the existing LP solution */
static
SCIP_RETCODE filterExistingLP(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   int*                  nfiltered,          /**< how many bounds were filtered this round? */
   BOUND*                currbound           /**< bound for which OBBT LP was solved (Note: might be NULL) */
   )
{
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(nfiltered != NULL);

   *nfiltered = 0;

   /* only apply filtering if an LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIPdebugMsg(scip, "can't filter using existing lp solution since it was not solved to optimality\n");
      return SCIP_OKAY;
   }

   /* check if a bound is tight */
   for( i = propdata->nbounds - 1; i >= 0; --i )
   {
      BOUND* bound;                          /* shortcut for current bound */

      SCIP_Real solval;                      /* the variables value in the current solution */
      SCIP_Real boundval;                    /* current local bound for the variable */

      bound = propdata->bounds[i];
      if( bound->filtered || bound->done )
         continue;

      boundval = bound->boundtype == SCIP_BOUNDTYPE_UPPER ?
         SCIPvarGetUbLocal(bound->var) : SCIPvarGetLbLocal(bound->var);
      solval = SCIPvarGetLPSol(bound->var);

      /* bound is tight; since this holds for all fixed variables, those are filtered here automatically; if the lp solution
       * is infinity, then also the bound is tight */
      if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER &&
               (SCIPisInfinity(scip, solval) || SCIPisFeasGE(scip, solval, boundval)))
            || (bound->boundtype == SCIP_BOUNDTYPE_LOWER &&
               (SCIPisInfinity(scip, -solval) || SCIPisFeasLE(scip, solval, boundval))) )
      {
         SCIP_BASESTAT basestat;

         /* mark bound as filtered */
         bound->filtered = TRUE;
         SCIPdebugMsg(scip, "trivial filtered var: %s boundval=%e solval=%e\n", SCIPvarGetName(bound->var), boundval, solval);

         /* get the basis status of the variable */
         basestat = SCIPcolGetBasisStatus(SCIPvarGetCol(bound->var));

         /* solve corresponding OBBT LP and try to generate a nontrivial genvbound */
         if( propdata->genvbdsduringfilter && currbound != NULL && basestat == SCIP_BASESTAT_BASIC )
         {
#ifndef NDEBUG
            int j;
#endif
            SCIP_Bool optimal;
            SCIP_Bool error;

            /* set objective coefficient of the bound */
            SCIP_CALL( SCIPchgVarObjProbing(scip, currbound->var, 0.0) );
            SCIP_CALL( setObjProbing(scip, propdata, bound, 1.0) );

#ifndef NDEBUG
            for( j = 0; j < SCIPgetNVars(scip); ++j )
            {
               SCIP_VAR* var;

               var = SCIPgetVars(scip)[j];
               assert(var != NULL);
               assert(SCIPisZero(scip, SCIPgetVarObjProbing(scip, var)) || var == bound->var);
            }
#endif

            /* solve the OBBT LP */
            propdata->nprobingiterations -= SCIPgetNLPIterations(scip);
            SCIP_CALL( solveLP(scip, -1, &error, &optimal) );
            propdata->nprobingiterations += SCIPgetNLPIterations(scip);
            assert(propdata->nprobingiterations >= 0);

            /* try to generate a genvbound if we have solved the OBBT LP */
            if( optimal && propdata->genvboundprop != NULL
                  && (SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1)) )
            {
               SCIP_Bool found;

               assert(!error);
               SCIP_CALL( createGenVBound(scip, propdata, bound, &found) );

               if( found )
               {
                  propdata->ngenvboundstrivfil += 1;
                  SCIPdebugMsg(scip, "found genvbound during trivial filtering\n");
               }
            }

            /* restore objective function */
            SCIP_CALL( setObjProbing(scip, propdata, bound, 0.0) );
            SCIP_CALL( setObjProbing(scip, propdata, currbound, 1.0) );
         }

         /* exchange bound i with propdata->bounds[propdata->lastidx] */
         if( propdata->lastidx >= 0 )
            exchangeBounds(propdata, i);

         /* increase number of filtered variables */
         (*nfiltered)++;
      }
   }

   /* try to filter bilinear bounds */
   for( i = propdata->lastbilinidx; i < propdata->nbilinbounds; ++i )
   {
      CORNER corners[4] = {LEFTTOP, LEFTBOTTOM, RIGHTTOP, RIGHTBOTTOM};
      BILINBOUND* bilinbound = propdata->bilinbounds[i];
      SCIP_Real solx;
      SCIP_Real soly;
      SCIPdebug(int oldfiltered;)
      int j;

      /* skip processed and filtered bounds */
      if( bilinbound->done || bilinbound->filtered == FILTERED ) /*lint !e641*/
         continue;

      SCIPdebug(oldfiltered = bilinbound->filtered;)
      solx = SCIPvarGetLPSol(bilinbound->x);
      soly = SCIPvarGetLPSol(bilinbound->y);

      /* check cases of unbounded solution values */
      if( SCIPisInfinity(scip, solx) )
         bilinbound->filtered = bilinbound->filtered | RIGHTTOP | RIGHTBOTTOM; /*lint !e641*/
      else if( SCIPisInfinity(scip, -solx) )
         bilinbound->filtered = bilinbound->filtered | LEFTTOP | LEFTBOTTOM; /*lint !e641*/

      if( SCIPisInfinity(scip, soly) )
         bilinbound->filtered = bilinbound->filtered | RIGHTTOP | LEFTTOP; /*lint !e641*/
      else if( SCIPisInfinity(scip, -soly) )
         bilinbound->filtered = bilinbound->filtered | RIGHTBOTTOM | LEFTBOTTOM; /*lint !e641*/

      /* check all corners */
      for( j = 0; j < 4; ++j )
      {
         SCIP_Real xt = SCIP_INVALID;
         SCIP_Real yt = SCIP_INVALID;

         getCorner(bilinbound->x, bilinbound->y, corners[j], &xt, &yt);

         if( (SCIPisInfinity(scip, REALABS(solx)) || SCIPisFeasEQ(scip, xt, solx))
            && (SCIPisInfinity(scip, REALABS(soly)) || SCIPisFeasEQ(scip, yt, soly)) )
            bilinbound->filtered = bilinbound->filtered | corners[j]; /*lint !e641*/
      }

#ifdef SCIP_DEBUG
      if( oldfiltered != bilinbound->filtered )
      {
         SCIP_VAR* x = bilinbound->x;
         SCIP_VAR* y = bilinbound->y;
         SCIPdebugMessage("filtered corners %d for (%s,%s) = (%g,%g) in [%g,%g]x[%g,%g]\n",
            bilinbound->filtered - oldfiltered, SCIPvarGetName(x), SCIPvarGetName(y), solx, soly,
            SCIPvarGetLbGlobal(x), SCIPvarGetUbGlobal(x), SCIPvarGetLbGlobal(y), SCIPvarGetUbGlobal(y));
      }
#endif
   }

   return SCIP_OKAY;
}

/** enforces one round of filtering */
static
SCIP_RETCODE filterRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   int                   itlimit,            /**< LP iteration limit (-1: no limit) */
   int*                  nfiltered,          /**< how many bounds were filtered this round */
   SCIP_Real*            objcoefs,           /**< array to store the nontrivial objective coefficients */
   int*                  objcoefsinds,       /**< array to store bound indices for which their corresponding variables
                                               *  has a nontrivial objective coefficient */
   int                   nobjcoefs           /**< number of nontrivial objective coefficients */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   SCIP_Bool error;
   SCIP_Bool optimal;

   int nvars;                                /* number of the problems variables */
   int i;

   assert(scip != NULL);
   assert(SCIPinProbing(scip));
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);
   assert(nfiltered != NULL);
   assert(objcoefs != NULL);
   assert(objcoefsinds != NULL);
   assert(nobjcoefs >= 0);

   *nfiltered = 0;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* solve LP */
   propdata->nfilterlpiters -= (int) SCIPgetNLPIterations(scip);
   SCIP_CALL( solveLP(scip, itlimit, &error, &optimal) );
   propdata->nfilterlpiters += (int) SCIPgetNLPIterations(scip);
   assert(propdata->nfilterlpiters >= 0);

   if( !optimal )
   {
      SCIPdebugMsg(scip, "skipping filter round since the LP was not solved to optimality\n");
      return SCIP_OKAY;
   }

   assert(!error);

   /* check if a bound is tight */
   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound;                          /* shortcut for current bound */

      SCIP_Real solval;                      /* the variables value in the current solution */
      SCIP_Real boundval;                    /* current local bound for the variable */

      bound = propdata->bounds[i];

      /* if bound is filtered it was handled already before */
      if( bound->filtered )
         continue;

      boundval = bound->boundtype == SCIP_BOUNDTYPE_UPPER ?
         SCIPvarGetUbLocal(bound->var) : SCIPvarGetLbLocal(bound->var);
      solval = SCIPvarGetLPSol(bound->var);

      /* bound is tight */
      if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisFeasGE(scip, solval, boundval))
         || (bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisFeasLE(scip, solval, boundval)) )
      {
         SCIP_Real objcoef;
         SCIP_BASESTAT basestat;

         /* mark bound as filtered */
         bound->filtered = TRUE;

         /* get the basis status of the variable */
         basestat = SCIPcolGetBasisStatus(SCIPvarGetCol(bound->var));

         /* increase number of filtered variables */
         (*nfiltered)++;

         /* solve corresponding OBBT LP and try to generate a nontrivial genvbound */
         if( propdata->genvbdsduringfilter && basestat == SCIP_BASESTAT_BASIC )
         {
            int j;

            /* set all objective coefficients to zero */
            for( j = 0; j < nobjcoefs; ++j )
            {
               BOUND* filterbound;

               filterbound = propdata->bounds[ objcoefsinds[j] ];
               assert(filterbound != NULL);

               SCIP_CALL( SCIPchgVarObjProbing(scip, filterbound->var, 0.0) );
            }

#ifndef NDEBUG
            for( j = 0; j < nvars; ++j )
               assert(SCIPisZero(scip, SCIPgetVarObjProbing(scip, vars[j])));
#endif

            /* set objective coefficient of the bound */
            SCIP_CALL( setObjProbing(scip, propdata, bound, 1.0) );

            /* solve the OBBT LP */
            propdata->nfilterlpiters -= (int) SCIPgetNLPIterations(scip);
            SCIP_CALL( solveLP(scip, -1, &error, &optimal) );
            propdata->nfilterlpiters += (int) SCIPgetNLPIterations(scip);
            assert(propdata->nfilterlpiters >= 0);

            /* try to generate a genvbound if we have solved the OBBT LP */
            if( optimal && propdata->genvboundprop != NULL
                  && (SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1)) )
            {
               SCIP_Bool found;

               assert(!error);
               SCIP_CALL( createGenVBound(scip, propdata, bound, &found) );

               if( found )
               {
                  propdata->ngenvboundsaggrfil += 1;
                  SCIPdebugMsg(scip, "found genvbound during aggressive filtering\n");
               }

            }

            /* restore objective function */
            for( j = 0; j < nobjcoefs; ++j )
            {
               BOUND* filterbound;

               filterbound = propdata->bounds[ objcoefsinds[j] ];
               assert(filterbound != NULL);

               /* NOTE: only restore coefficients of nonfiltered bounds */
               if( !filterbound->filtered )
               {
                  assert(!SCIPisZero(scip, objcoefs[j]));
                  SCIP_CALL( SCIPchgVarObjProbing(scip, propdata->bounds[ objcoefsinds[j] ]->var, objcoefs[j]) );
               }
            }
         }

         /* get the corresponding variable's objective coefficient */
         objcoef = SCIPgetVarObjProbing(scip, bound->var);

         /* change objective coefficient if it was set up for this bound */
          if( (bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisNegative(scip, objcoef))
             || (bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisPositive(scip, objcoef)) )
          {
             SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, 0.0) );
          }
      }
   }

   return SCIP_OKAY;
}

/** filter some bounds that are not improvable by solving auxiliary LPs */
static
SCIP_RETCODE filterBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint          itlimit             /**< LP iteration limit (-1: no limit) */
   )
{
   SCIP_VAR** vars;
   SCIP_Longint nolditerations;
   SCIP_Real* objcoefs;               /* array to store the nontrivial objective coefficients */
   int* objcoefsinds;                 /* array to store bound indices for which the corresponding variable
                                       * has a nontrivial objective coefficient */
   int nobjcoefs;                     /* number of nontrivial objective coefficients */
   int nleftiterations;
   int i;
   int nfiltered;
   int ntotalfiltered;
   int nvars;

   assert(scip != NULL);
   assert(SCIPinProbing(scip));
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   ntotalfiltered = 0;
   nolditerations = SCIPgetNLPIterations(scip);
   nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIPdebugMsg(scip, "start filter rounds\n");

   SCIP_CALL( SCIPallocBufferArray(scip, &objcoefs, propdata->nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objcoefsinds, propdata->nbounds) );
   nobjcoefs = 0;

   /*
    * 1.) Try first to filter lower bounds of interesting variables, whose bounds are not already filtered
    */

   for( i = 0; i < nvars; i++ )
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, vars[i], 0.0) );
   }

   for( i = 0; i < propdata->nbounds; i++ )
   {
      if( propdata->bounds[i]->boundtype == SCIP_BOUNDTYPE_LOWER && !propdata->bounds[i]->filtered
            && !propdata->bounds[i]->done )
      {
         SCIP_Real objcoef;

         objcoef = getFilterCoef(scip, propdata, propdata->bounds[i]->var, SCIP_BOUNDTYPE_LOWER);

         if( !SCIPisZero(scip, objcoef) )
         {
            SCIP_CALL( SCIPchgVarObjProbing(scip, propdata->bounds[i]->var, objcoef) );

            /* store nontrivial objective coefficients */
            objcoefs[nobjcoefs] = objcoef;
            objcoefsinds[nobjcoefs] = i;
            ++nobjcoefs;
         }
      }
   }

   do
   {
      SCIPdebugMsg(scip, "doing a lower bounds round\n");
      SCIP_CALL( filterRound(scip, propdata, nleftiterations, &nfiltered, objcoefs, objcoefsinds, nobjcoefs) );
      ntotalfiltered += nfiltered;
      SCIPdebugMsg(scip, "filtered %d more bounds in lower bounds round\n", nfiltered);

      /* update iterations left */
      nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
   }
   while( nfiltered >= propdata->nminfilter && ( nleftiterations == -1 ||  nleftiterations > 0 ) );

   /*
    * 2.) Now try to filter the remaining upper bounds of interesting variables, whose bounds are not already filtered
    */

   /* set all objective coefficients to zero */
   for( i = 0; i < nobjcoefs; i++ )
   {
      BOUND* bound;

      assert(objcoefsinds[i] >= 0 && objcoefsinds[i] < propdata->nbounds);
      bound = propdata->bounds[ objcoefsinds[i] ];
      assert(bound != NULL);
      SCIP_CALL( SCIPchgVarObjProbing(scip, bound->var, 0.0) );
   }

   /* reset number of nontrivial objective coefficients */
   nobjcoefs = 0;

#ifndef NDEBUG
   for( i = 0; i < nvars; ++i )
      assert(SCIPisZero(scip, SCIPgetVarObjProbing(scip, vars[i])));
#endif

   for( i = 0; i < propdata->nbounds; i++ )
   {
      if( propdata->bounds[i]->boundtype == SCIP_BOUNDTYPE_UPPER && !propdata->bounds[i]->filtered )
      {
         SCIP_Real objcoef;

         objcoef = getFilterCoef(scip, propdata, propdata->bounds[i]->var, SCIP_BOUNDTYPE_UPPER);

         if( !SCIPisZero(scip, objcoef) )
         {
            SCIP_CALL( SCIPchgVarObjProbing(scip, propdata->bounds[i]->var, objcoef) );

            /* store nontrivial objective coefficients */
            objcoefs[nobjcoefs] = objcoef;
            objcoefsinds[nobjcoefs] = i;
            ++nobjcoefs;
         }
      }
   }

   do
   {
      SCIPdebugMsg(scip, "doing an upper bounds round\n");
      SCIP_CALL( filterRound(scip, propdata, nleftiterations, &nfiltered, objcoefs, objcoefsinds, nobjcoefs) );
      SCIPdebugMsg(scip, "filtered %d more bounds in upper bounds round\n", nfiltered);
      ntotalfiltered += nfiltered;
      /* update iterations left */
      nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
   }
   while( nfiltered >= propdata->nminfilter && ( nleftiterations == -1 ||  nleftiterations > 0 ) );

   SCIPdebugMsg(scip, "filtered %d this round\n", ntotalfiltered);
   propdata->nfiltered += ntotalfiltered;

   /* free array */
   SCIPfreeBufferArray(scip, &objcoefsinds);
   SCIPfreeBufferArray(scip, &objcoefs);

   return SCIP_OKAY;
}

/** applies possible bound changes that were found */
static
SCIP_RETCODE applyBoundChgs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
#ifdef SCIP_DEBUG
   int ntightened;                           /* stores the number of successful bound changes */
#endif
   int i;

   assert(scip != NULL);
   assert(!SCIPinProbing(scip));
   assert(propdata != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND);

   SCIPdebug( ntightened = 0 );

   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound;                          /* shortcut to the current bound */
      SCIP_Bool infeas;                      /* stores wether a tightening approach forced an infeasibilty */
      SCIP_Bool tightened;                   /* stores wether a tightening approach was successful */

      bound = propdata->bounds[i];

      if( bound->found )
      {
         SCIPdebug( double oldbound = (bound->boundtype == SCIP_BOUNDTYPE_LOWER)
            ? SCIPvarGetLbLocal(bound->var)
            : SCIPvarGetUbLocal(bound->var) );

         if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, bound->var, bound->newval, FALSE, &infeas, &tightened) );
         }
         else
         {
            SCIP_CALL( SCIPtightenVarUb(scip, bound->var, bound->newval, FALSE, &infeas, &tightened) );
         }

         /* handle information about the success */
         if( infeas )
         {
            *result = SCIP_CUTOFF;
            SCIPdebugMsg(scip, "cut off\n");
            break;
         }

         if( tightened )
         {
            SCIPdebug( SCIPdebugMsg(scip, "tightended: %s old: %e new: %e\n" , SCIPvarGetName(bound->var), oldbound,
                  bound->newval) );
            *result = SCIP_REDUCEDDOM;
            SCIPdebug( ntightened++ );
         }
      }
   }

   SCIPdebug( SCIPdebugMsg(scip, "tightened bounds: %d\n", ntightened) );

   return SCIP_OKAY;
}

/** tries to tighten a bound in probing mode  */
static
SCIP_RETCODE tightenBoundProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   BOUND*                bound,              /**< bound that could be tightened */
   SCIP_Real             newval,             /**< new bound value */
   SCIP_Bool*            tightened           /**< was tightening successful? */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(SCIPinProbing(scip));
   assert(bound != NULL);
   assert(tightened != NULL);

   *tightened = FALSE;

   /* get old bounds */
   lb = SCIPvarGetLbLocal(bound->var);
   ub = SCIPvarGetUbLocal(bound->var);

   if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      /* round bounds new value if variable is integral */
      if( SCIPvarIsIntegral(bound->var) )
         newval = SCIPceil(scip, newval);

      /* ensure that we give consistent bounds to the LP solver */
      if( newval > ub )
         newval = ub;

      /* tighten if really better */
      if( SCIPisLbBetter(scip, newval, lb, ub) )
      {
         SCIP_CALL( SCIPchgVarLbProbing(scip, bound->var, newval) );
         *tightened = TRUE;
      }
   }
   else
   {
      /* round bounds new value if variable is integral */
      if( SCIPvarIsIntegral(bound->var) )
         newval = SCIPfloor(scip, newval);

      /* ensure that we give consistent bounds to the LP solver */
      if( newval < lb )
         newval = lb;

      /* tighten if really better */
      if( SCIPisUbBetter(scip, newval, lb, ub) )
      {
         SCIP_CALL( SCIPchgVarUbProbing(scip, bound->var, newval) );
         *tightened = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** comparison method for two bounds w.r.t. their scores */
static
SCIP_DECL_SORTPTRCOMP(compBoundsScore)
{
   BOUND* bound1 = (BOUND*) elem1;
   BOUND* bound2 = (BOUND*) elem2;

   return bound1->score == bound2->score ? 0 : ( bound1->score > bound2->score ? 1 : -1 );
}

/** comparison method for two bilinear term bounds w.r.t. their scores */
static
SCIP_DECL_SORTPTRCOMP(compBilinboundsScore)
{
   BILINBOUND* bound1 = (BILINBOUND*) elem1;
   BILINBOUND* bound2 = (BILINBOUND*) elem2;

   return bound1->score == bound2->score ? 0 : ( bound1->score > bound2->score ? 1 : -1 ); /*lint !e777*/
}

/** comparison method for two bounds w.r.t. their boundtype */
static
SCIP_DECL_SORTPTRCOMP(compBoundsBoundtype)
{
   int diff;
   BOUND* bound1 = (BOUND*) elem1;
   BOUND* bound2 = (BOUND*) elem2;

   /* prioritize undone bounds */
   diff = (!bound1->done ? 1 : 0) - (!bound2->done ? 1 : 0);
   if( diff != 0 )
      return diff;

   /* prioritize unfiltered bounds */
   diff = (!bound1->filtered ? 1 : 0) - (!bound2->filtered ? 1 : 0);
   if( diff != 0 )
      return diff;

   diff = (bound1->boundtype == SCIP_BOUNDTYPE_LOWER ? 1 : 0) - (bound2->boundtype == SCIP_BOUNDTYPE_LOWER ? 1 : 0);

   if( diff == 0 )
      return (bound1->score == bound2->score) ? 0 : (bound1->score > bound2->score ? 1 : -1);
   else
      return diff;
}

/** sort the propdata->bounds array with their distance or their boundtype key */
static
SCIP_RETCODE sortBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);

   SCIPdebugMsg(scip, "sort bounds\n");
   SCIPsortDownPtr((void**) propdata->bounds, compBoundsBoundtype, propdata->nbounds);

   return SCIP_OKAY;
}

/** evaluates a bound for the current LP solution */
static
SCIP_Real evalBound(
   SCIP*                 scip,
   BOUND*                bound
   )
{
   assert(scip != NULL);
   assert(bound != NULL);

   if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
      return REALABS( SCIPvarGetLPSol(bound->var) - SCIPvarGetLbLocal(bound->var) );
   else
      return REALABS( SCIPvarGetUbLocal(bound->var) - SCIPvarGetLPSol(bound->var) );
}

/** returns the index of the next undone and unfiltered bound with the smallest distance */
static
int nextBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Bool             convexphase         /**< consider only convex variables? */
   )
{
   SCIP_Real bestval;
   int bestidx;
   int k;

   assert(scip != NULL);
   assert(propdata != NULL);

   bestidx = -1;
   bestval = SCIPinfinity(scip);

   for( k = 0; k <= propdata->lastidx; ++k )
   {
      BOUND* tmpbound;
      tmpbound = propdata->bounds[k];

      assert(tmpbound != NULL);

      if( !tmpbound->filtered && !tmpbound->done && (tmpbound->nonconvex == !convexphase) )
      {
         SCIP_Real boundval;

         /* return the next bound which is not done or unfiltered yet */
         if( propdata->orderingalgo == 0 )
            return k;

         boundval = evalBound(scip, tmpbound);

         /* negate boundval if we use the reverse greedy algorithm */
         boundval = (propdata->orderingalgo == 2) ? -1.0 * boundval : boundval;

         if( bestidx == -1 || boundval < bestval )
         {
            bestidx = k;
            bestval = boundval;
         }
      }
   }

   return bestidx;
}

/** try to separate the solution of the last OBBT LP in order to learn better variable bounds; we apply additional
 *  separation rounds as long as the routine finds better bounds; because of dual degeneracy we apply a minimum number of
 *  separation rounds
 */
static
SCIP_RETCODE applySeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   BOUND*                currbound,          /**< current bound */
   SCIP_Longint*         nleftiterations,    /**< number of left iterations (-1 for no limit) */
   SCIP_Bool*            success             /**< pointer to store if we have found a better bound */
   )
{
   SCIP_Bool inroot;
   int i;

   assert(nleftiterations != NULL);
   assert(success != NULL);
   assert(SCIPinProbing(scip));

   *success = FALSE;

   /* check if we are originally in the root node */
   inroot = SCIPgetDepth(scip) == 1;

   for( i = 0; i <= propdata->sepamaxiter; ++i )
   {
      SCIP_Longint nlpiter;
      SCIP_Real oldval;
      SCIP_Bool cutoff;
      SCIP_Bool delayed;
      SCIP_Bool error;
      SCIP_Bool optimal;
      SCIP_Bool tightened;

      oldval = SCIPvarGetLPSol(currbound->var);

      /* find and store cuts to separate the current LP solution */
      SCIP_CALL( SCIPseparateSol(scip, NULL, inroot, TRUE, FALSE, &delayed, &cutoff) );
      SCIPdebugMsg(scip, "applySeparation() - ncuts = %d\n", SCIPgetNCuts(scip));

      /* leave if we did not found any cut */
      if( SCIPgetNCuts(scip) == 0 )
         break;

      /* apply cuts and resolve LP */
      SCIP_CALL( SCIPapplyCutsProbing(scip, &cutoff) );
      assert(SCIPgetNCuts(scip) == 0);
      SCIPdebug( nlpiter = SCIPgetNLPIterations(scip); )
      SCIP_CALL( solveLP(scip, (int) *nleftiterations, &error, &optimal) );
      SCIPdebug( nlpiter = SCIPgetNLPIterations(scip) - nlpiter; )
      SCIPdebugMsg(scip, "applySeparation() - optimal=%u error=%u lpiter=%" SCIP_LONGINT_FORMAT "\n", optimal, error, nlpiter);
      SCIPdebugMsg(scip, "oldval = %e newval = %e\n", oldval, SCIPvarGetLPSol(currbound->var));

      /* leave if we did not solve the LP to optimality or an error occured */
      if( error || !optimal )
         break;

      /* try to generate a genvbound */
      if( inroot && propdata->genvboundprop != NULL && propdata->genvbdsduringsepa )
      {
         SCIP_Bool found;
         SCIP_CALL( createGenVBound(scip, propdata, currbound, &found) );
         propdata->ngenvboundsprobing += found ? 1 : 0;
      }

      /* try to tight the variable bound */
      tightened = FALSE;
      if( !SCIPisEQ(scip, oldval, SCIPvarGetLPSol(currbound->var)) )
      {
         SCIP_CALL( tightenBoundProbing(scip, currbound, SCIPvarGetLPSol(currbound->var), &tightened) );
         SCIPdebugMsg(scip, "apply separation - tightened=%u oldval=%e newval=%e\n", tightened, oldval,
            SCIPvarGetLPSol(currbound->var));

         *success |= tightened;
      }

      /* leave the separation if we did not tighten the bound and proceed at least propdata->sepaminiter iterations */
      if( !tightened && i >= propdata->sepaminiter )
         break;
   }

   return SCIP_OKAY;
}

/** finds new variable bounds until no iterations left or all bounds have been checked */
static
SCIP_RETCODE findNewBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint*         nleftiterations,    /**< pointer to store the number of left iterations */
   SCIP_Bool             convexphase         /**< consider only convex variables? */
   )
{
   SCIP_Longint nolditerations;
   SCIP_Bool iterationsleft;
   BOUND* currbound;
   SCIP_Longint itlimit;
   int nextboundidx;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(nleftiterations != NULL);

   /* update the number of left iterations */
   nolditerations = SCIPgetNLPIterations(scip);
   itlimit = *nleftiterations;
   assert(*nleftiterations == getIterationsLeft(scip, nolditerations, itlimit));
   iterationsleft = (*nleftiterations == -1) || (*nleftiterations > 0);

   /* To improve the performance we sort the bound in such a way that the undone and
    * unfiltered bounds are at the end of propdata->bounds. We calculate and update
    * the position of the last unfiltered and undone bound in propdata->lastidx
    */
   if( !convexphase )
   {
      /* sort bounds */
      SCIP_CALL( sortBounds(scip, propdata) );

      /* if the first bound is filtered or done then there is no bound left */
      if( propdata->bounds[0]->done || propdata->bounds[0]->filtered )
      {
         SCIPdebugMsg(scip, "no unprocessed/unfiltered bound left\n");
         return SCIP_OKAY;
      }

      /* compute the last undone and unfiltered node */
      propdata->lastidx = 0;
      while( propdata->lastidx < propdata->nbounds - 1 && !propdata->bounds[propdata->lastidx]->done &&
            !propdata->bounds[propdata->lastidx]->filtered )
         ++propdata->lastidx;

      SCIPdebugMsg(scip, "lastidx = %d\n", propdata->lastidx);
   }

   /* find the first unprocessed bound */
   nextboundidx = nextBound(scip, propdata, convexphase);

   /* skip if there is no bound left */
   if( nextboundidx == -1 )
   {
      SCIPdebugMsg(scip, "no unprocessed/unfiltered bound left\n");
      return SCIP_OKAY;
   }

   currbound = propdata->bounds[nextboundidx];
   assert(!currbound->done && !currbound->filtered);

   /* main loop */
   while( iterationsleft &&  !SCIPisStopped(scip) )
   {
      SCIP_Bool optimal;
      SCIP_Bool error;
      int nfiltered;

      assert(currbound != NULL);
      assert(currbound->done == FALSE);
      assert(currbound->filtered == FALSE);

      /* do not visit currbound more than once */
      currbound->done = TRUE;
      exchangeBounds(propdata, nextboundidx);

      /* set objective for curr */
      SCIP_CALL( setObjProbing(scip, propdata, currbound, 1.0) );

      SCIPdebugMsg(scip, "before solving      Boundtype: %d , LB: %e , UB: %e\n",
         currbound->boundtype == SCIP_BOUNDTYPE_LOWER, SCIPvarGetLbLocal(currbound->var),
         SCIPvarGetUbLocal(currbound->var) );
      SCIPdebugMsg(scip, "before solving      var <%s>, LP value: %f\n",
         SCIPvarGetName(currbound->var), SCIPvarGetLPSol(currbound->var));

      SCIPdebugMsg(scip, "probing iterations before solve: %lld \n", SCIPgetNLPIterations(scip));

      propdata->nprobingiterations -= SCIPgetNLPIterations(scip);

      /* now solve the LP */
      SCIP_CALL( solveLP(scip, (int) *nleftiterations, &error, &optimal) );

      propdata->nprobingiterations += SCIPgetNLPIterations(scip);
      propdata->nsolvedbounds++;

      SCIPdebugMsg(scip, "probing iterations after solve: %lld \n", SCIPgetNLPIterations(scip));
      SCIPdebugMsg(scip, "OPT: %u ERROR: %u\n" , optimal, error);
      SCIPdebugMsg(scip, "after solving      Boundtype: %d , LB: %e , UB: %e\n",
         currbound->boundtype == SCIP_BOUNDTYPE_LOWER, SCIPvarGetLbLocal(currbound->var),
         SCIPvarGetUbLocal(currbound->var) );
      SCIPdebugMsg(scip, "after solving      var <%s>, LP value: %f\n",
         SCIPvarGetName(currbound->var), SCIPvarGetLPSol(currbound->var));

      /* update nleftiterations */
      *nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
      iterationsleft = (*nleftiterations == -1) || (*nleftiterations > 0);

      if( error )
      {
         SCIPdebugMsg(scip, "ERROR during LP solving\n");

         /* set the objective of currbound to zero to null the whole objective; otherwise the objective is wrong when
          * we call findNewBounds() for the convex phase
          */
         SCIP_CALL( SCIPchgVarObjProbing(scip, currbound->var, 0.0) );

         return SCIP_OKAY;
      }

      if( optimal )
      {
         SCIP_Bool success;

         currbound->newval = SCIPvarGetLPSol(currbound->var);
         currbound->found = TRUE;

         /* in root node we may want to create a genvbound (independent of tightening success) */
         if( (SCIPgetDepth(scip) == 0 || (SCIPinProbing(scip) && SCIPgetDepth(scip) == 1))
               && propdata->genvboundprop != NULL )
         {
            SCIP_Bool found;

            SCIP_CALL( createGenVBound(scip, propdata, currbound, &found) );

            if( found )
               propdata->ngenvboundsprobing += 1;
         }

         /* try to tighten bound in probing mode */
         success = FALSE;
         if( propdata->tightintboundsprobing && SCIPvarIsIntegral(currbound->var) )
         {
            SCIPdebugMsg(scip, "tightening bound %s = %e bounds: [%e, %e]\n", SCIPvarGetName(currbound->var),
                currbound->newval, SCIPvarGetLbLocal(currbound->var), SCIPvarGetUbLocal(currbound->var) );
            SCIP_CALL( tightenBoundProbing(scip, currbound, currbound->newval, &success) );
            SCIPdebugMsg(scip, "tightening bound %s\n", success ? "successful" : "not successful");
         }
         else if( propdata->tightcontboundsprobing && !SCIPvarIsIntegral(currbound->var) )
         {
            SCIPdebugMsg(scip, "tightening bound %s = %e bounds: [%e, %e]\n", SCIPvarGetName(currbound->var),
               currbound->newval, SCIPvarGetLbLocal(currbound->var), SCIPvarGetUbLocal(currbound->var) );
            SCIP_CALL( tightenBoundProbing(scip, currbound, currbound->newval, &success) );
            SCIPdebugMsg(scip, "tightening bound %s\n", success ? "successful" : "not successful");
         }

         /* separate current OBBT LP solution */
         if( iterationsleft && propdata->separatesol )
         {
            propdata->nprobingiterations -= SCIPgetNLPIterations(scip);
            SCIP_CALL( applySeparation(scip, propdata, currbound, nleftiterations, &success) );
            propdata->nprobingiterations += SCIPgetNLPIterations(scip);

            /* remember best solution value after solving additional separations LPs */
            if( success )
            {
#ifndef NDEBUG
               SCIP_Real newval = SCIPvarGetLPSol(currbound->var);

               /* round new bound if the variable is integral */
               if( SCIPvarIsIntegral(currbound->var) )
                  newval = currbound->boundtype == SCIP_BOUNDTYPE_LOWER ?
                     SCIPceil(scip, newval) : SCIPfloor(scip, newval);

               assert((currbound->boundtype == SCIP_BOUNDTYPE_LOWER &&
                     SCIPisGT(scip, newval, currbound->newval))
                  || (currbound->boundtype == SCIP_BOUNDTYPE_UPPER &&
                     SCIPisLT(scip, newval, currbound->newval)));
#endif

               currbound->newval = SCIPvarGetLPSol(currbound->var);
            }
         }

         /* filter bound candidates by using the current LP solution */
         if( propdata->applytrivialfilter )
         {
            SCIP_CALL( filterExistingLP(scip, propdata, &nfiltered, currbound) );
            SCIPdebugMsg(scip, "filtered %d bounds via inspecting present LP solution\n", nfiltered);
            propdata->ntrivialfiltered += nfiltered;
         }

         propdata->propagatecounter += success ? 1 : 0;

         /* propagate if we have found enough bound tightenings */
         if( propdata->propagatefreq != 0 && propdata->propagatecounter >= propdata->propagatefreq )
         {
            SCIP_Longint ndomredsfound;
            SCIP_Bool cutoff;

            SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, &ndomredsfound) );
            SCIPdebugMsg(scip, "propagation - cutoff %u  ndomreds %" SCIP_LONGINT_FORMAT "\n", cutoff, ndomredsfound);

            propdata->npropagatedomreds += ndomredsfound;
            propdata->propagatecounter = 0;
         }
      }

      /* set objective to zero */
      SCIP_CALL( setObjProbing(scip, propdata, currbound, 0.0) );

      /* find the first unprocessed bound */
      nextboundidx = nextBound(scip, propdata, convexphase);

      /* check if there is no unprocessed and unfiltered node left */
      if( nextboundidx == -1 )
      {
         SCIPdebugMsg(scip, "NO unvisited/unfiltered bound left!\n");
         break;
      }

      currbound = propdata->bounds[nextboundidx];
      assert(!currbound->done && !currbound->filtered);
   }

   if( iterationsleft )
   {
      SCIPdebugMsg(scip, "still iterations left: %" SCIP_LONGINT_FORMAT "\n", *nleftiterations);
   }
   else
   {
      SCIPdebugMsg(scip, "no iterations left\n");
   }

   return SCIP_OKAY;
}


/** main function of obbt */
static
SCIP_RETCODE applyObbt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint          itlimit,            /**< LP iteration limit (-1: no limit) */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* oldlbs;
   SCIP_Real* oldubs;
   SCIP_Longint lastnpropagatedomreds;
   SCIP_Longint nleftiterations;
   SCIP_Real oldconditionlimit;
   SCIP_Real oldboundstreps;
   SCIP_Real olddualfeastol;
   SCIP_Bool hasconditionlimit;
   SCIP_Bool continuenode;
   SCIP_Bool boundleft;
   int oldpolishing;
   int nfiltered;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);

   SCIPdebugMsg(scip, "apply obbt\n");

   oldlbs = NULL;
   oldubs = NULL;
   lastnpropagatedomreds = propdata->npropagatedomreds;
   nleftiterations = itlimit;
   continuenode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == propdata->lastnode;
   propdata->lastidx = -1;
   boundleft = FALSE;
   *result = SCIP_DIDNOTFIND;

   /* store old variable bounds if we use propagation during obbt */
   if( propdata->propagatefreq > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &oldlbs, propdata->nbounds) );
      SCIP_CALL( SCIPallocBufferArray(scip, &oldubs, propdata->nbounds) );
   }

   /* reset bound data structure flags; fixed variables are marked as filtered */
   for( i = 0; i < propdata->nbounds; i++ )
   {
      BOUND* bound = propdata->bounds[i];
      bound->found = FALSE;

      /* store old variable bounds */
      if( oldlbs != NULL && oldubs != NULL )
      {
         oldlbs[bound->index] = SCIPvarGetLbLocal(bound->var);
         oldubs[bound->index] = SCIPvarGetUbLocal(bound->var);
      }

      /* reset 'done' and 'filtered' flag in a new B&B node */
      if( !continuenode )
      {
         bound->done = FALSE;
         bound->filtered = FALSE;
      }

      /* mark fixed variables as filtered */
      bound->filtered |= varIsFixedLocal(scip, bound->var);

      /* check for an unprocessed bound */
      if( !bound->filtered && !bound->done )
         boundleft = TRUE;
   }

   /* no bound left to check */
   if( !boundleft )
      goto TERMINATE;

   /* filter variables via inspecting present LP solution */
   if( propdata->applytrivialfilter && !continuenode )
   {
      SCIP_CALL( filterExistingLP(scip, propdata, &nfiltered, NULL) );
      SCIPdebugMsg(scip, "filtered %d bounds via inspecting present LP solution\n", nfiltered);
      propdata->ntrivialfiltered += nfiltered;
   }

   /* store old dualfeasibletol */
   olddualfeastol = SCIPdualfeastol(scip);

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );
   SCIPdebugMsg(scip, "start probing\n");

   /* tighten dual feastol */
   if( propdata->dualfeastol < olddualfeastol )
   {
      SCIP_CALL( SCIPchgDualfeastol(scip, propdata->dualfeastol) );
   }

   /* tighten condition limit */
   hasconditionlimit = (SCIPgetRealParam(scip, "lp/conditionlimit", &oldconditionlimit) == SCIP_OKAY);
   if( !hasconditionlimit )
   {
      SCIPwarningMessage(scip, "obbt propagator could not set condition limit in LP solver - running without\n");
   }
   else if( propdata->conditionlimit > 0.0 && (oldconditionlimit < 0.0 || propdata->conditionlimit < oldconditionlimit) )
   {
      SCIP_CALL( SCIPsetRealParam(scip, "lp/conditionlimit", propdata->conditionlimit) );
   }

   /* tighten relative bound improvement limit */
   SCIP_CALL( SCIPgetRealParam(scip, "numerics/boundstreps", &oldboundstreps) );
   if( !SCIPisEQ(scip, oldboundstreps, propdata->boundstreps) )
   {
     SCIP_CALL( SCIPsetRealParam(scip, "numerics/boundstreps", propdata->boundstreps) );
   }

   /* add objective cutoff */
   SCIP_CALL( addObjCutoff(scip, propdata) );

   /* deactivate LP polishing */
   SCIP_CALL( SCIPgetIntParam(scip, "lp/solutionpolishing", &oldpolishing) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solutionpolishing", 0) );

   /* apply filtering */
   if( propdata->applyfilterrounds )
   {
      SCIP_CALL( filterBounds(scip, propdata, nleftiterations) );
   }

   /* set objective coefficients to zero */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( i = 0; i < nvars; ++i )
   {
      /* note that it is not possible to change the objective of non-column variables during probing; we have to take
       * care of the objective contribution of loose variables in createGenVBound()
       */
      if( SCIPvarGetObj(vars[i]) != 0.0 && SCIPvarGetStatus(vars[i]) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPchgVarObjProbing(scip, vars[i], 0.0) );
      }
   }

   /* find new bounds for the variables */
   SCIP_CALL( findNewBounds(scip, propdata, &nleftiterations, FALSE) );

   if( nleftiterations > 0 || itlimit < 0 )
   {
      SCIP_CALL( findNewBounds(scip, propdata, &nleftiterations, TRUE) );
   }

   /* reset dual feastol and condition limit */
   SCIP_CALL( SCIPchgDualfeastol(scip, olddualfeastol) );
   if( hasconditionlimit )
   {
      SCIP_CALL( SCIPsetRealParam(scip, "lp/conditionlimit", oldconditionlimit) );
   }

   /* update bound->newval if we have learned additional bound tightenings during SCIPpropagateProbing() */
   if( oldlbs != NULL && oldubs != NULL && propdata->npropagatedomreds - lastnpropagatedomreds > 0 )
   {
      assert(propdata->propagatefreq > 0);
      for( i = 0; i < propdata->nbounds; ++i )
      {
         BOUND* bound = propdata->bounds[i];

         /* it might be the case that a bound found by the additional propagation is better than the bound found after solving an OBBT
          * LP
          */
         if( bound->found )
         {
            if( bound->boundtype == SCIP_BOUNDTYPE_LOWER )
               bound->newval = MAX(bound->newval, SCIPvarGetLbLocal(bound->var)); /*lint !e666*/
            else
               bound->newval = MIN(bound->newval, SCIPvarGetUbLocal(bound->var)); /*lint !e666*/
         }
         else
         {
            SCIP_Real oldlb;
            SCIP_Real oldub;

            oldlb = oldlbs[bound->index];
            oldub = oldubs[bound->index];

            if( bound->boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisLbBetter(scip, SCIPvarGetLbLocal(bound->var), oldlb, oldub) )
            {
               SCIPdebugMsg(scip, "tighter lower bound due to propagation: %d - %e -> %e\n", i, oldlb, SCIPvarGetLbLocal(bound->var));
               bound->newval = SCIPvarGetLbLocal(bound->var);
               bound->found = TRUE;
            }

            if( bound->boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisUbBetter(scip, SCIPvarGetUbLocal(bound->var), oldlb, oldub) )
            {
               SCIPdebugMsg(scip, "tighter upper bound due to propagation: %d - %e -> %e\n", i, oldub, SCIPvarGetUbLocal(bound->var));
               bound->newval = SCIPvarGetUbLocal(bound->var);
               bound->found = TRUE;
            }
         }
      }
   }

   /* reset relative bound improvement limit */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/boundstreps", oldboundstreps) );

   /* reset original LP polishing setting */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solutionpolishing", oldpolishing) );

   /* end probing */
   SCIP_CALL( SCIPendProbing(scip) );
   SCIPdebugMsg(scip, "end probing!\n");

   /* release cutoff row if there is one */
   if( propdata->cutoffrow != NULL )
   {
      assert(!SCIProwIsInLP(propdata->cutoffrow));
      SCIP_CALL( SCIPreleaseRow(scip, &(propdata->cutoffrow)) );
   }

   /* apply buffered bound changes */
   SCIP_CALL( applyBoundChgs(scip, propdata, result) );

TERMINATE:
   SCIPfreeBufferArrayNull(scip, &oldubs);
   SCIPfreeBufferArrayNull(scip, &oldlbs);

   return SCIP_OKAY;
}

/** computes a valid inequality from the current LP relaxation for a bilinear term xy only involving x and y; the
 *  inequality is found by optimizing along the line connecting the points (xs,ys) and (xt,yt) over the currently given
 *  linear relaxation of the problem; this optimization problem is an LP
 *
 *  max lambda
 *  s.t. Ax <= b
 *       (x,y) = (xs,ys) + lambda ((xt,yt) - (xs,ys))
 *       lambda in [0,1]
 *
 *  which is equivalent to
 *
 *  max x
 *  s.t. (1) Ax <= b
 *       (2) (x - xs) / (xt - xs) = (y - ys) / (yt - ys)
 *
 *  Let x* be the optimal primal and (mu,theta) be the optimal dual solution of this LP. The KKT conditions imply that
 *  the aggregation of the linear constraints mu*Ax <= mu*b can be written as
 *
 *  x * (1 - theta) / (xt - xs) + y * theta / (yt - ys) = mu * Ax <= mu * b
 *
 *  <=> alpha * x + beta * y <= mu * b = alpha * (x*) + beta * (y*)
 *
 *  which is a valid inequality in the (x,y)-space; in order to avoid numerical difficulties when (xs,ys) is too close
 *  to (xt,yt), we scale constraint (1) by max{1,|xt-xs|,|yt-ys|} beforehand
 */
static
SCIP_RETCODE solveBilinearLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_Real             xs,                 /**< x-coordinate of the first point */
   SCIP_Real             ys,                 /**< y-coordinate of the first point */
   SCIP_Real             xt,                 /**< x-coordinate of the second point */
   SCIP_Real             yt,                 /**< y-coordinate of the second point */
   SCIP_Real*            xcoef,              /**< pointer to store the coefficient of x */
   SCIP_Real*            ycoef,              /**< pointer to store the coefficient of y */
   SCIP_Real*            constant,           /**< pointer to store the constant */
   SCIP_Longint          iterlim             /**< iteration limit (-1: for no limit) */
   )
{
   SCIP_ROW* row;
   SCIP_Real signx;
   SCIP_Real scale;
   SCIP_Real side;
   SCIP_Bool lperror;

   assert(xcoef != NULL);
   assert(ycoef != NULL);
   assert(constant != NULL);
   assert(SCIPinProbing(scip));

   *xcoef = SCIP_INVALID;
   *ycoef = SCIP_INVALID;
   *constant= SCIP_INVALID;

   SCIPdebugMsg(scip, "   solve bilinear LP for (%s,%s) from (%g,%g) to (%g,%g)\n", SCIPvarGetName(x), SCIPvarGetName(y), xs,
      ys, xt, yt);

   /* skip computations if (xs,ys) and (xt,yt) are too close to each other or contain too large values */
   if( SCIPisFeasEQ(scip, xs, xt) || SCIPisFeasEQ(scip, ys, yt)
      || SCIPisHugeValue(scip, REALABS(xs)) || SCIPisHugeValue(scip, REALABS(xt))
      || SCIPisHugeValue(scip, REALABS(ys)) || SCIPisHugeValue(scip, REALABS(yt)) )
   {
      SCIPdebugMsg(scip, "   -> skip: bounds are too close/large\n");
      return SCIP_OKAY;
   }

   /* compute scaler for the additional linear constraint */
   scale = MIN(MAX3(1.0, REALABS(xt-xs), REALABS(yt-ys)), 100.0); /*lint !e666*/

   /* set objective function */
   signx = (xs > xt) ? 1.0 : -1.0;
   SCIP_CALL( SCIPchgVarObjProbing(scip, x, signx) );

   /* create new probing node to remove the added LP row afterwards */
   SCIP_CALL( SCIPnewProbingNode(scip) );

   /* create row */
   side = scale * (xs/(xt-xs) - ys/(yt-ys));
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, "bilinrow", side, side, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, x, scale/(xt-xs)) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, y, -scale/(yt-ys)) );
   SCIP_CALL( SCIPaddRowProbing(scip, row) );

   /* solve probing LP */
#ifdef NDEBUG
   {
      SCIP_RETCODE retstat;
      retstat = SCIPsolveProbingLP(scip, iterlim, &lperror, NULL);
      if( retstat != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving LP in quadratic constraint handler; LP solve terminated with" \
            "code <%d>\n", retstat);
      }
   }
#else
   SCIP_CALL( SCIPsolveProbingLP(scip, (int)iterlim, &lperror, NULL) ); /*lint !e712*/
#endif

   SCIPdebugMsg(scip, "   solved probing LP -> lperror=%u lpstat=%d\n", lperror, SCIPgetLPSolstat(scip));

   /* collect dual and primal solution entries */
   if( !lperror  && SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Real xval = SCIPvarGetLPSol(x);
      SCIP_Real yval = SCIPvarGetLPSol(y);
      SCIP_Real mu = -SCIProwGetDualsol(row);

      SCIPdebugMsg(scip, "   primal=(%g,%g) dual=%g\n", xval, yval, mu);

      /* xcoef x + ycoef y <= constant */
      *xcoef  = -signx - (mu * scale) / (xt - xs);
      *ycoef = (mu * scale) / (yt - ys);
      *constant = (*xcoef) * xval + (*ycoef) * yval;

      /* xcoef x <= -ycoef y + constant */
      *ycoef = -(*ycoef);

      /* inequality is only useful when both coefficients are different from zero; normalize inequality if possible */
      if( !SCIPisFeasZero(scip, *xcoef) && !SCIPisFeasZero(scip, *ycoef) )
      {
         SCIP_Real val = REALABS(*xcoef);
         *xcoef /= val;
         *ycoef /= val;
         *constant /= val;
      }
      else
      {
         *xcoef = SCIP_INVALID;
         *ycoef = SCIP_INVALID;
         *constant = SCIP_INVALID;
      }
   }

   /* release row and backtrack probing node */
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIP_CALL( SCIPbacktrackProbing(scip, 0) );

   /* reset objective function */
   SCIP_CALL( SCIPchgVarObjProbing(scip, x, 0.0) );

   return SCIP_OKAY;
}

/* applies obbt for finding valid inequalities for bilinear terms; function works as follows:
 *
 *  1. start probing mode
 *  2. add objective cutoff (if possible)
 *  3. set objective function to zero
 *  4. set feasibility, optimality, and relative bound improvement tolerances of SCIP
 *  5. main loop
 *  6. restore old tolerances
 *
 */
static
SCIP_RETCODE applyObbtBilinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of the obbt propagator */
   SCIP_Longint          itlimit,            /**< LP iteration limit (-1: no limit) */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_VAR** vars;
   SCIP_Real oldfeastol;
   SCIP_Bool lperror;
   SCIP_Longint nolditerations;
   SCIP_Longint nleftiterations;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(itlimit == -1 || itlimit >= 0);
   assert(result != NULL);

   if( propdata->nbilinbounds <= 0 || SCIPgetDepth(scip) != 0 || propdata->lastbilinidx >= propdata->nbilinbounds )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "call applyObbtBilinear starting from %d\n", propdata->lastbilinidx);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   nolditerations = SCIPgetNLPIterations(scip);
   nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
   SCIPdebugMsg(scip, "iteration limit: %lld\n", nleftiterations);

   /* 1. start probing */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* 2. add objective cutoff */
   SCIP_CALL( addObjCutoff(scip, propdata) );

   /* 3. set objective function to zero */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, vars[i], 0.0) );
   }

   /* we need to solve the probing LP before creating new probing nodes in solveBilinearLP() */
   SCIP_CALL( SCIPsolveProbingLP(scip, (int)nleftiterations, &lperror, NULL) );

   /* 4. tighten LP feasibility tolerance to be at most feastol/10.0 */
   oldfeastol = SCIPchgRelaxfeastol(scip, SCIPfeastol(scip) / 10.0);

   /* 5. main loop */
   for( i = propdata->lastbilinidx; i < propdata->nbilinbounds
      && (nleftiterations > 0 || nleftiterations == -1)
      && (propdata->itlimitbilin < 0 || propdata->itlimitbilin > propdata->itusedbilin )
      && !SCIPisStopped(scip); ++i )
   {
      CORNER corners[4] = {LEFTBOTTOM, LEFTTOP, RIGHTTOP, RIGHTBOTTOM};
      BILINBOUND* bilinbound;
      int k;

      bilinbound = propdata->bilinbounds[i];
      assert(bilinbound != NULL);

      SCIPdebugMsg(scip, "process %d: %s %s done=%u filtered=%d nunderest=%d noverest=%d\n", i,
         SCIPvarGetName(bilinbound->x), SCIPvarGetName(bilinbound->y), bilinbound->done, bilinbound->filtered,
         bilinbound->nunderest, bilinbound->noverest);

      /* we already solved LPs for this bilinear term */
      if( bilinbound->done || bilinbound->filtered == (int)FILTERED )
         continue;

      /* iterate through all corners
       *
       * 0: (xs,ys)=(ubx,lby) (xt,yt)=(lbx,uby) -> underestimate
       * 1: (xs,ys)=(ubx,uby) (xt,yt)=(lbx,lby) -> overestimate
       * 2: (xs,ys)=(lbx,uby) (xt,yt)=(ubx,lby) -> underestimate
       * 3: (xs,ys)=(lbx,lby) (xt,yt)=(ubx,uby) -> overestimate
       */
      for( k = 0; k < 4; ++k )
      {
         CORNER corner = corners[k];
         SCIP_Real xcoef;
         SCIP_Real ycoef;
         SCIP_Real constant;
         SCIP_Real xs = SCIP_INVALID;
         SCIP_Real ys = SCIP_INVALID;
         SCIP_Real xt = SCIP_INVALID;
         SCIP_Real yt = SCIP_INVALID;

         /* skip corners that lead to an under- or overestimate that is not needed */
         if( ((corner == LEFTTOP || corner == RIGHTBOTTOM) && bilinbound->nunderest == 0)
            || ((corner == LEFTBOTTOM || corner == RIGHTTOP) && bilinbound->noverest == 0) )
            continue;

         /* check whether corner has been filtered already */
         if( (bilinbound->filtered & corner) != 0 ) /*lint !e641*/
            continue;

         /* get corners (xs,ys) and (xt,yt) */
         getCorners(bilinbound->x, bilinbound->y, corner, &xs, &ys, &xt, &yt);

         /* skip target corner points with too large values */
         if( SCIPisHugeValue(scip, REALABS(xt)) || SCIPisHugeValue(scip, REALABS(yt)) )
            continue;

         /* compute inequality */
         propdata->itusedbilin -= SCIPgetNLPIterations(scip);
         SCIP_CALL( solveBilinearLP(scip, bilinbound->x, bilinbound->y, xs, ys, xt, yt, &xcoef, &ycoef, &constant, -1L) );
         propdata->itusedbilin += SCIPgetNLPIterations(scip);

         /* update number of LP iterations */
         nleftiterations = getIterationsLeft(scip, nolditerations, itlimit);
         SCIPdebugMsg(scip, "LP iterations left: %lld\n", nleftiterations);

         /* add inequality to quadratic constraint handler if it separates (xt,yt) */
         if( !SCIPisHugeValue(scip, xcoef)  && !SCIPisFeasZero(scip, xcoef)
            && REALABS(ycoef) < 1e+3 && REALABS(ycoef) > 1e-3 /* TODO add a parameter for this */
            && SCIPisFeasGT(scip, (xcoef*xt - ycoef*yt - constant) / SQRT(SQR(xcoef) + SQR(ycoef) + SQR(constant)), 1e-2) )
         {
            SCIP_Bool success;

            SCIP_CALL( SCIPaddBilinearIneqQuadratic(scip, bilinbound->x, bilinbound->y, bilinbound->index, xcoef,
               ycoef, constant, &success) );

            /* check whether the inequality has been accepted by the quadratic constraint handler */
            if( success )
            {
               *result = SCIP_REDUCEDDOM;
               SCIPdebugMsg(scip, "   found %g x <= %g y + %g with violation %g\n", xcoef, ycoef, constant,
                  (xcoef*xt - ycoef*yt - constant) / SQRT(SQR(xcoef) + SQR(ycoef) + SQR(constant)));
            }
         }
      }

      /* mark the bound as processed */
      bilinbound->done = TRUE;
   }

   /* remember last unprocessed bilinear term */
   propdata->lastbilinidx = i;

   /* end probing */
   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) ); /* TODO necessary to solve LP here again? */
   SCIP_CALL( SCIPendProbing(scip) );

   /* release cutoff row if there is one */
   if( propdata->cutoffrow != NULL )
   {
      assert(!SCIProwIsInLP(propdata->cutoffrow));
      SCIP_CALL( SCIPreleaseRow(scip, &(propdata->cutoffrow)) );
   }

   /* 6. restore old tolerance */
   (void) SCIPchgRelaxfeastol(scip, oldfeastol);

   return SCIP_OKAY;
}

/** computes the score of a bound */
static
unsigned int getScore(
   SCIP*                 scip,               /**< SCIP data structure */
   BOUND*                bound,              /**< pointer of bound */
   int                   nlcount,            /**< number of nonlinear constraints containing the bounds variable */
   int                   maxnlcount          /**< maximal number of nonlinear constraints a variable appears in */
   )
{
   unsigned int score;                       /* score to be computed */

   assert(scip != NULL);
   assert(bound != NULL);
   assert(nlcount >= 0);
   assert(maxnlcount >= nlcount);

   /* score = ( nlcount * ( BASE - 1 ) / maxnlcount ) * BASE^2 + vartype * BASE + boundtype */
   score = (unsigned int) ( nlcount > 0 ? (OBBT_SCOREBASE * nlcount * ( OBBT_SCOREBASE - 1 )) / maxnlcount : 0 );
   switch( SCIPvarGetType(bound->var) )
   {
   case SCIP_VARTYPE_INTEGER:
      score += 1;
      break;
   case SCIP_VARTYPE_IMPLINT:
      score += 2;
      break;
   case SCIP_VARTYPE_CONTINUOUS:
      score += 3;
      break;
   case SCIP_VARTYPE_BINARY:
      score += 4;
      break;
   default:
      break;
   }

   score *= OBBT_SCOREBASE;
   if( bound->boundtype == SCIP_BOUNDTYPE_UPPER )
      score += 1;

   return score;
}

/** computes the score of a bilinear term bound */
static
SCIP_Real getScoreBilinBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   BILINBOUND*           bilinbound,         /**< bilinear term bound */
   int                   nbilinterms         /**< maximal number of bilinear terms in all quadratic constraints */
   )
{
   SCIP_Real lbx = SCIPvarGetLbLocal(bilinbound->x);
   SCIP_Real ubx = SCIPvarGetUbLocal(bilinbound->x);
   SCIP_Real lby = SCIPvarGetLbLocal(bilinbound->y);
   SCIP_Real uby = SCIPvarGetUbLocal(bilinbound->y);
   SCIP_Real score;

   assert(scip != NULL);
   assert(randnumgen != NULL);
   assert(bilinbound != NULL);

   /* consider how often a bilinear term is present in the problem */
   score = (bilinbound->noverest + bilinbound->nunderest) / (SCIP_Real)nbilinterms;

   /* penalize small variable domains; TODO tune the factor in the logarithm, maybe add a parameter for it */
   if( ubx - lbx < 0.5 )
      score += log(2.0*(ubx-lbx) + SCIPepsilon(scip));
   if( uby - lby < 0.5 )
      score += log(2.0*(uby-lby) + SCIPepsilon(scip));

   /* consider interiority of variables in the LP solution */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Real solx = SCIPvarGetLPSol(bilinbound->x);
      SCIP_Real soly = SCIPvarGetLPSol(bilinbound->y);
      SCIP_Real interiorityx = MIN(solx-lbx, ubx-solx) / MAX(ubx-lbx, SCIPepsilon(scip)); /*lint !e666*/
      SCIP_Real interiorityy = MIN(soly-lby, uby-soly) / MAX(uby-lby, SCIPepsilon(scip)); /*lint !e666*/

      score += interiorityx + interiorityy;
   }

   /* randomize score */
   score *= 1.0 + SCIPrandomGetReal(randnumgen, -SCIPepsilon(scip), SCIPepsilon(scip));

   return score;
}

/** count the variables which appear in non-convex term of nlrow  */
static
SCIP_RETCODE countNLRowVarsNonConvexity(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcounts,           /**< store the number each variable appears in a
                                              *   non-convex term */
   SCIP_NLROW*           nlrow               /**< nonlinear row */
   )
{
   int t;
   int nexprtreevars;
   SCIP_VAR** exprtreevars;
   SCIP_EXPRTREE* exprtree;

   assert(scip != NULL);
   assert(nlcounts != NULL);
   assert(nlrow != NULL);

   /* go through all quadratic terms */
   for( t = SCIPnlrowGetNQuadElems(nlrow) - 1; t >= 0; --t )
   {
      SCIP_QUADELEM* quadelem;
      SCIP_VAR* bilinvar1;
      SCIP_VAR* bilinvar2;

      /* get quadratic term */
      quadelem = &SCIPnlrowGetQuadElems(nlrow)[t];

      /* get involved variables */
      bilinvar1 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx1];
      bilinvar2 = SCIPnlrowGetQuadVars(nlrow)[quadelem->idx2];

      assert(bilinvar1 != NULL);
      assert(bilinvar2 != NULL);

      /* we have a non-convex square term */
      if( bilinvar1 == bilinvar2 && !(quadelem->coef >= 0 ? SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)) : SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow))) )
      {
         ++nlcounts[SCIPvarGetProbindex(bilinvar1)];
         ++nlcounts[SCIPvarGetProbindex(bilinvar2)];
      }

      /* bilinear terms are in general non-convex */
      if( bilinvar1 != bilinvar2 )
      {
         ++nlcounts[SCIPvarGetProbindex(bilinvar1)];
         ++nlcounts[SCIPvarGetProbindex(bilinvar2)];
      }
   }

   exprtree = SCIPnlrowGetExprtree(nlrow);
   if( exprtree != NULL )
   {
      nexprtreevars = SCIPexprtreeGetNVars(exprtree);
      exprtreevars = SCIPexprtreeGetVars(exprtree);

      /* assume that the expression tree represents a non-convex constraint */
      for( t = 0; t < nexprtreevars; ++t)
      {
         SCIP_VAR* var;
         var = exprtreevars[t];
         assert(var != NULL);

         ++nlcounts[SCIPvarGetProbindex(var)];
      }
   }

   return SCIP_OKAY;
}

/** count how often each variable appears in a non-convex term */
static
SCIP_RETCODE getNLPVarsNonConvexity(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcounts            /**< store the number each variable appears in a
                                              *   non-convex term */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nvars;
   int nconss;
   int i;

   assert(scip != NULL);
   assert(nlcounts != NULL);

   nvars = SCIPgetNVars(scip);
   BMSclearMemoryArray(nlcounts, nvars);

   /* quadratic constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "quadratic");
   if( conshdlr != NULL )
   {

      /*SCIPdebugMsg(scip, "cons_quadratic is there!\n");*/
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMsg(scip, "nconss(quadratic) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_Bool isnonconvex;

         isnonconvex = (!SCIPisConvexQuadratic(scip, conss[i]) && !SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, conss[i])))
            || (!SCIPisConcaveQuadratic(scip, conss[i]) && !SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, conss[i])));

         /* only check the nlrow if the constraint is not convex */
         if( isnonconvex )
         {
            SCIP_NLROW* nlrow;
            SCIP_CALL( SCIPgetNlRowQuadratic(scip, conss[i], &nlrow) );
            assert(nlrow != NULL);

            SCIP_CALL( countNLRowVarsNonConvexity(scip, nlcounts, nlrow) );
         }
      }
   }

   /* nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   if( conshdlr != NULL )
   {
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMsg(scip, "nconss(nonlinear) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_EXPRCURV curvature;
         SCIP_Bool isnonconvex;

         SCIP_CALL( SCIPgetCurvatureNonlinear(scip, conss[i], TRUE, &curvature) );

         isnonconvex = (curvature != SCIP_EXPRCURV_CONVEX && !SCIPisInfinity(scip, SCIPgetRhsNonlinear(scip, conss[i])))
            || (curvature != SCIP_EXPRCURV_CONCAVE && !SCIPisInfinity(scip, -SCIPgetLhsNonlinear(scip, conss[i])));

         /* only check the nlrow if the constraint is not convex */
         if( isnonconvex )
         {
            SCIP_NLROW* nlrow;
            SCIP_CALL( SCIPgetNlRowNonlinear(scip, conss[i], &nlrow) );
            assert(nlrow != NULL);

            SCIP_CALL( countNLRowVarsNonConvexity(scip, nlcounts, nlrow) );
         }
      }
   }

   /* bivariate constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "bivariate");
   if( conshdlr != NULL )
   {
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMsg(scip, "nconss(bivariate) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         SCIP_EXPRCURV curvature;
         SCIP_INTERVAL* varbounds;
         SCIP_EXPRTREE* exprtree;
         int j;

         exprtree = SCIPgetExprtreeBivariate(scip, conss[i]);
         if( exprtree != NULL )
         {
            SCIP_Bool isnonconvex;

            SCIP_CALL( SCIPallocBufferArray(scip, &varbounds, SCIPexprtreeGetNVars(exprtree)) );
            for( j = 0; j < SCIPexprtreeGetNVars(exprtree); ++j )
            {
               SCIP_VAR* var;
               var = SCIPexprtreeGetVars(exprtree)[j];

               SCIPintervalSetBounds(&varbounds[j],
                  -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))),    /*lint !e666*/
                  +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))) );  /*lint !e666*/
            }

            SCIP_CALL( SCIPexprtreeCheckCurvature(exprtree, SCIPinfinity(scip), varbounds, &curvature, NULL) );

            isnonconvex = (curvature != SCIP_EXPRCURV_CONVEX && !SCIPisInfinity(scip, SCIPgetRhsBivariate(scip, conss[i])))
               || (curvature != SCIP_EXPRCURV_CONCAVE && !SCIPisInfinity(scip, -SCIPgetLhsBivariate(scip, conss[i])));

            /* increase counter for all variables in the expression tree if the constraint is non-convex */
            if( isnonconvex )
            {
               for( j = 0; j < SCIPexprtreeGetNVars(exprtree); ++j )
               {
                  SCIP_VAR* var;
                  var = SCIPexprtreeGetVars(exprtree)[j];

                  ++nlcounts[SCIPvarGetProbindex(var)];
               }
            }
            SCIPfreeBufferArray(scip, &varbounds);
         }
      }
   }

   /* abspower constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "abspower");
   if( conshdlr != NULL )
   {
      nconss = SCIPconshdlrGetNActiveConss(conshdlr);
      conss = SCIPconshdlrGetConss(conshdlr);

      SCIPdebugMsg(scip, "nconss(abspower) = %d\n", nconss);

      for( i = 0; i < nconss; ++i )
      {
         /* constraint is non-convex in general */
         SCIP_NLROW* nlrow;
         SCIP_CALL( SCIPgetNlRowAbspower(scip, conss[i], &nlrow) );
         assert(nlrow != NULL);

         SCIP_CALL( countNLRowVarsNonConvexity(scip, nlcounts, nlrow) );
      }
   }

   return SCIP_OKAY;
}


/** determines whether a variable is interesting */
static
SCIP_Bool varIsInteresting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to check */
   int                   nlcount             /**< number of nonlinear constraints containing the variable
                                               *  or number of non-convex terms containing the variable
                                               * (depends on propdata->onlynonconvexvars)  */
   )
{
   assert(SCIPgetDepth(scip) == 0);

   return !SCIPvarIsBinary(var) && SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && nlcount > 0
      && !varIsFixedLocal(scip, var);
}

/** initializes interesting bounds */
static
SCIP_RETCODE initBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< data of the obbt propagator */
   )
{
   SCIP_VAR** vars;                          /* array of the problems variables */
   int* nlcount;                             /* array that stores in how many nonlinearities each variable appears */
   int* nccount;                             /* array that stores in how many nonconvexities each variable appears */

   int bdidx;                                /* bound index inside propdata->bounds */
   int maxnlcount;                           /* maximal number of nonlinear constraints a variable appears in */
   int nvars;                                /* number of the problems variables */
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(SCIPisNLPConstructed(scip));

   SCIPdebugMsg(scip, "initialize bounds\n");

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* count nonlinearities */
   assert(SCIPgetNNLPVars(scip) == nvars);

   SCIP_CALL( SCIPallocBufferArray(scip, &nlcount, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nccount, nvars) );

   SCIP_CALL( SCIPgetNLPVarsNonlinearity(scip, nlcount) );
   SCIP_CALL( getNLPVarsNonConvexity(scip, nccount) );

   maxnlcount = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( maxnlcount < nlcount[i] )
         maxnlcount = nlcount[i];
   }

   /* allocate interesting bounds array */
   propdata->boundssize = 2 * nvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->bounds), 2 * nvars) );

   /* get all interesting variables and their bounds */
   bdidx = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( varIsInteresting(scip, vars[i], (propdata->onlynonconvexvars ? nccount[i] : nlcount[i])) )
      {
         BOUND** bdaddress;

         /* create lower bound */
         bdaddress = &(propdata->bounds[bdidx]);
         SCIP_CALL( SCIPallocBlockMemory(scip, bdaddress) );
         propdata->bounds[bdidx]->boundtype = SCIP_BOUNDTYPE_LOWER;
         propdata->bounds[bdidx]->var = vars[i];
         propdata->bounds[bdidx]->found = FALSE;
         propdata->bounds[bdidx]->filtered = FALSE;
         propdata->bounds[bdidx]->newval = 0.0;
         propdata->bounds[bdidx]->score = getScore(scip, propdata->bounds[bdidx], nlcount[i], maxnlcount);
         propdata->bounds[bdidx]->done = FALSE;
         propdata->bounds[bdidx]->nonconvex = (nccount[i] > 0);
         propdata->bounds[bdidx]->index = bdidx;
         bdidx++;

         /* create upper bound */
         bdaddress = &(propdata->bounds[bdidx]);
         SCIP_CALL( SCIPallocBlockMemory(scip, bdaddress) );
         propdata->bounds[bdidx]->boundtype = SCIP_BOUNDTYPE_UPPER;
         propdata->bounds[bdidx]->var = vars[i];
         propdata->bounds[bdidx]->found = FALSE;
         propdata->bounds[bdidx]->filtered = FALSE;
         propdata->bounds[bdidx]->newval = 0.0;
         propdata->bounds[bdidx]->score = getScore(scip, propdata->bounds[bdidx], nlcount[i], maxnlcount);
         propdata->bounds[bdidx]->done = FALSE;
         propdata->bounds[bdidx]->nonconvex = (nccount[i] > 0);
         propdata->bounds[bdidx]->index = bdidx;
         bdidx++;
      }
   }

   /* set number of interesting bounds */
   propdata->nbounds = bdidx;

   /* collect all bilinear terms from quadratic constraint handler */
   if( propdata->nbounds > 0 && SCIPgetNAllBilinearTermsQuadratic(scip) > 0 && propdata->createbilinineqs )
   {
      SCIP_VAR** x;
      SCIP_VAR** y;
      SCIP_Real* maxnonconvexity;
      int* nunderest;
      int* noverest;
      int nbilins;
      int bilinidx;
      int nbilinterms;

      nbilins = SCIPgetNAllBilinearTermsQuadratic(scip);
      bilinidx = 0;
      nbilinterms = 0;


      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &x, nbilins) );
      SCIP_CALL( SCIPallocBufferArray(scip, &y, nbilins) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nunderest, nbilins) );
      SCIP_CALL( SCIPallocBufferArray(scip, &noverest, nbilins) );
      SCIP_CALL( SCIPallocBufferArray(scip, &maxnonconvexity, nbilins) );

      /* get data for bilinear terms */
      SCIP_CALL( SCIPgetAllBilinearTermsQuadratic(scip, x, y, &nbilins, nunderest, noverest, maxnonconvexity) );

      /* count the number of interesting bilinear terms */
      propdata->nbilinbounds = 0;
      for( i = 0; i < nbilins; ++i )
         if( nunderest[i] + noverest[i] > 0 && propdata->minnonconvexity <= maxnonconvexity[i]
            && varIsInteresting(scip, x[i], 1) && varIsInteresting(scip, y[i], 1) )
            ++(propdata->nbilinbounds);

      if( propdata->nbilinbounds == 0 )
         goto TERMINATE;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(propdata->bilinbounds), propdata->nbilinbounds) );
      BMSclearMemoryArray(propdata->bilinbounds, propdata->nbilinbounds);

      for( i = 0; i < nbilins; ++i )
      {
         if( nunderest[i] + noverest[i] > 0 && propdata->minnonconvexity <= maxnonconvexity[i]
            && varIsInteresting(scip, x[i], 1) && varIsInteresting(scip, y[i], 1) )
         {
            SCIP_CALL( SCIPallocBlockMemory(scip, &propdata->bilinbounds[bilinidx]) ); /*lint !e866*/
            BMSclearMemory(propdata->bilinbounds[bilinidx]); /*lint !e866*/

            propdata->bilinbounds[bilinidx]->x = x[i];
            propdata->bilinbounds[bilinidx]->y = y[i];
            propdata->bilinbounds[bilinidx]->nunderest = nunderest[i];
            propdata->bilinbounds[bilinidx]->noverest = noverest[i];
            propdata->bilinbounds[bilinidx]->index = i;
            ++bilinidx;

            /* count how often bilinear terms appear in quadratic constraints */
            nbilinterms += nunderest[i] + noverest[i];
         }
      }
      assert(propdata->nbilinbounds == bilinidx);

      /* compute scores for each term */
      for( i = 0; i < propdata->nbilinbounds; ++i )
      {
         propdata->bilinbounds[i]->score = getScoreBilinBound(scip, propdata->randnumgen, propdata->bilinbounds[i],
            nbilinterms);
         SCIPdebugMsg(scip, "score of %i = %g\n", i, propdata->bilinbounds[i]->score);
      }

      /* sort bounds according to decreasing score */
      if( propdata->nbilinbounds > 1 )
      {
         SCIPsortDownPtr((void**) propdata->bilinbounds, compBilinboundsScore, propdata->nbilinbounds);
      }

TERMINATE:
      /* free memory */
      SCIPfreeBufferArray(scip, &maxnonconvexity);
      SCIPfreeBufferArray(scip, &noverest);
      SCIPfreeBufferArray(scip, &nunderest);
      SCIPfreeBufferArray(scip, &y);
      SCIPfreeBufferArray(scip, &x);
   }

   /* free memory for buffering nonlinearities */
   assert(nlcount != NULL);
   assert(nccount != NULL);
   SCIPfreeBufferArray(scip, &nccount);
   SCIPfreeBufferArray(scip, &nlcount);

   /*  propdata->bounds array if empty */
   if( propdata->nbounds <= 0 )
   {
      assert(propdata->nbounds == 0);
      assert(propdata->boundssize >= 0 );
      SCIPfreeBlockMemoryArray(scip, &(propdata->bounds), propdata->boundssize);
   }

   SCIPdebugMsg(scip, "problem has %d/%d interesting bounds\n", propdata->nbounds, 2 * nvars);

   if( propdata->nbounds > 0 )
   {
      /* sort bounds according to decreasing score; although this initial order will be overruled by the distance
       * criterion later, gives a more well-defined starting situation for OBBT and might help to reduce solver
       * variability
       */
      SCIPsortDownPtr((void**) propdata->bounds, compBoundsScore, propdata->nbounds);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of propagator
 */

/** copy method for propagator plugins (called when SCIP copies plugins)
 *
 *  @note The UG framework assumes that all default plug-ins of SCIP implement a copy callback. We check
 *  SCIPgetSubscipDepth() in PROPEXEC to prevent the propagator to run in a sub-SCIP.
 */
static
SCIP_DECL_PROPCOPY(propCopyObbt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludePropObbt(scip) );

   return SCIP_OKAY;
}

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->bounds = NULL;
   propdata->nbounds = -1;
   propdata->boundssize = 0;
   propdata->cutoffrow = NULL;
   propdata->lastnode = -1;

   /* if genvbounds propagator is not available, we cannot create genvbounds */
   propdata->genvboundprop = propdata->creategenvbounds ? SCIPfindProp(scip, GENVBOUND_PROP_NAME) : NULL;

   SCIPdebugMsg(scip, "creating genvbounds: %s\n", propdata->genvboundprop != NULL ? "true" : "false");

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &propdata->randnumgen, DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Longint itlimit;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   *result = SCIP_DIDNOTRUN;

   /* do not run in: presolving, repropagation, probing mode, if no objective propagation is allowed  */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPinRepropagation(scip) || SCIPinProbing(scip) || !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* do not run propagator in a sub-SCIP */
   if( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   /* only run for nonlinear problems, i.e., if NLP is constructed */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "NLP not constructed, skipping obbt\n");
      return SCIP_OKAY;
   }

   /* only run if LP all columns are in the LP, i.e., the LP is a relaxation; e.g., do not run if pricers are active
    * since pricing is not performed in probing mode
    */
   if( !SCIPallColsInLP(scip) )
   {
      SCIPdebugMsg(scip, "not all columns in LP, skipping obbt\n");
      return SCIP_OKAY;
   }

   if( !SCIPallowObjProp(scip) )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* ensure that bounds are initialized */
   if( propdata->nbounds == -1 )
   {
      /* bounds must be initialized at root node */
      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( initBounds(scip, propdata) );
      }
      else
      {
         assert(!SCIPinProbing(scip));
         return SCIP_OKAY;
      }
   }
   assert(propdata->nbounds >= 0);

   /* do not run if there are no interesting bounds */
   /**@todo disable */
   if( propdata->nbounds <= 0 )
   {
      SCIPdebugMsg(scip, "there are no interesting bounds\n");
      return SCIP_OKAY;
   }

   /* only run once in a node != root */
   if( SCIPgetDepth(scip) > 0 && SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == propdata->lastnode )
   {
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "applying obbt for problem <%s> at depth %d\n", SCIPgetProbName(scip), SCIPgetDepth(scip));

   /* without an optimal LP solution we don't want to run; this may be because propagators with higher priority have
    * already found reductions or numerical troubles occured during LP solving
    */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPdebugMsg(scip, "aborting since no optimal LP solution is at hand\n");
      return SCIP_OKAY;
   }

   /* compute iteration limit */
   if( propdata->itlimitfactor > 0.0 )
      itlimit = (SCIP_Longint) MAX(propdata->itlimitfactor * SCIPgetNRootLPIterations(scip),
         propdata->minitlimit); /*lint !e666*/
   else
      itlimit = -1;

   /* apply obbt */
   SCIP_CALL( applyObbt(scip, propdata, itlimit, result) );
   assert(*result != SCIP_DIDNOTRUN);

   /* compute globally inequalities for bilinear terms */
   if( propdata->createbilinineqs )
   {
      /* set LP iteration limit */
      if( propdata->itlimitbilin == 0L )
      {
         /* no iteration limit if itlimit < 0 or itlimitfactorbilin < 0 */
         propdata->itlimitbilin = (itlimit < 0 || propdata->itlimitfactorbilin < 0)
            ? -1L : (SCIP_Longint)(itlimit * propdata->itlimitfactorbilin);
      }

      SCIP_CALL( applyObbtBilinear(scip, propdata, itlimit, result) );
   }

   /* set current node as last node */
   propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   return SCIP_OKAY;
}

/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropObbt)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int i;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* free random number generator */
   SCIPfreeRandom(scip, &propdata->randnumgen);
   propdata->randnumgen = NULL;

   /* note that because we reset filtered flags to false at each call to obbt, the same bound may be filtered multiple
    * times
    */
   SCIPstatisticMessage("DIVE-LP: %" SCIP_LONGINT_FORMAT "  NFILTERED: %d NTRIVIALFILTERED: %d NSOLVED: %d "
      "FILTER-LP: %" SCIP_LONGINT_FORMAT " NGENVB(dive): %d NGENVB(aggr.): %d NGENVB(triv.) %d\n",
      propdata->nprobingiterations, propdata->nfiltered, propdata->ntrivialfiltered, propdata->nsolvedbounds,
      propdata->nfilterlpiters, propdata->ngenvboundsprobing, propdata->ngenvboundsaggrfil, propdata->ngenvboundstrivfil);

   /* free bilinear bounds */
   if( propdata->nbilinbounds > 0 )
   {
      for( i = propdata->nbilinbounds - 1; i >= 0; --i )
      {
         SCIPfreeBlockMemory(scip, &propdata->bilinbounds[i]); /*lint !e866*/
      }
      SCIPfreeBlockMemoryArray(scip, &propdata->bilinbounds, propdata->nbilinbounds);
      propdata->nbilinbounds = 0;
   }

   /* free memory allocated for the bounds */
   if( propdata->nbounds > 0 )
   {
      /* free bounds */
      for( i = propdata->nbounds - 1; i >= 0; i-- )
      {
         SCIPfreeBlockMemory(scip, &(propdata->bounds[i])); /*lint !e866*/
      }
      SCIPfreeBlockMemoryArray(scip, &(propdata->bounds), propdata->boundssize);
   }

   propdata->nbounds = -1;
   propdata->itlimitbilin = 0;
   propdata->itusedbilin = 0;

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeBlockMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the obbt propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropObbt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create obbt propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   BMSclearMemory(propdata);

   /* initialize statistic variables */
   propdata->nprobingiterations = 0;
   propdata->nfiltered = 0;
   propdata->ntrivialfiltered = 0;
   propdata->nsolvedbounds = 0;
   propdata->ngenvboundsprobing = 0;
   propdata->ngenvboundsaggrfil = 0;
   propdata->ngenvboundstrivfil = 0;
   propdata->nfilterlpiters = 0;
   propdata->lastidx = -1;
   propdata->propagatecounter = 0;
   propdata->npropagatedomreds = 0;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecObbt, propdata) );

   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyObbt) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeObbt) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolObbt) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolObbt) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropObbt) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/creategenvbounds",
         "should obbt try to provide genvbounds if possible?",
         &propdata->creategenvbounds, TRUE, DEFAULT_CREATE_GENVBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/normalize",
         "should coefficients in filtering be normalized w.r.t. the domains sizes?",
         &propdata->normalize, TRUE, DEFAULT_FILTERING_NORM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/applyfilterrounds",
         "try to filter bounds in so-called filter rounds by solving auxiliary LPs?",
         &propdata->applyfilterrounds, TRUE, DEFAULT_APPLY_FILTERROUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/applytrivialfilter",
         "try to filter bounds with the LP solution after each solve?",
         &propdata->applytrivialfilter, TRUE, DEFAULT_APPLY_TRIVIALFITLERING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/genvbdsduringfilter",
         "should we try to generate genvbounds during trivial and aggressive filtering?",
         &propdata->genvbdsduringfilter, TRUE, DEFAULT_GENVBDSDURINGFILTER, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/genvbdsduringsepa",
         "try to create genvbounds during separation process?",
         &propdata->genvbdsduringsepa, TRUE, DEFAULT_GENVBDSDURINGSEPA, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/minfilter",
         "minimal number of filtered bounds to apply another filter round",
         &propdata->nminfilter, TRUE, DEFAULT_FILTERING_MIN, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/itlimitfactor",
         "multiple of root node LP iterations used as total LP iteration limit for obbt (<= 0: no limit )",
         &propdata->itlimitfactor, FALSE, DEFAULT_ITLIMITFACTOR, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/itlimitfactorbilin",
         "multiple of OBBT LP limit used as total LP iteration limit for solving bilinear inequality LPs (< 0 for no limit)",
         &propdata->itlimitfactorbilin, FALSE, DEFAULT_ITLIMITFAC_BILININEQS, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/minnonconvexity",
         "minimum absolute value of nonconvex eigenvalues for a bilinear term",
         &propdata->minnonconvexity, FALSE, DEFAULT_MINNONCONVEXITY, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "propagating/" PROP_NAME "/minitlimit",
         "minimum LP iteration limit",
         &propdata->minitlimit, FALSE, DEFAULT_MINITLIMIT, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/dualfeastol",
         "feasibility tolerance for reduced costs used in obbt; this value is used if SCIP's dual feastol is greater",
         &propdata->dualfeastol, FALSE, DEFAULT_DUALFEASTOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/conditionlimit",
         "maximum condition limit used in LP solver (-1.0: no limit)",
         &propdata->conditionlimit, FALSE, DEFAULT_CONDITIONLIMIT, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "propagating/" PROP_NAME "/boundstreps",
         "minimal relative improve for strengthening bounds",
         &propdata->boundstreps, FALSE, DEFAULT_BOUNDSTREPS, 0.0, 1.0, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/onlynonconvexvars",
         "only apply obbt on non-convex variables",
         &propdata->onlynonconvexvars, TRUE, DEFAULT_ONLYNONCONVEXVARS, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/tightintboundsprobing",
         "should integral bounds be tightened during the probing mode?",
         &propdata->tightintboundsprobing, TRUE, DEFAULT_TIGHTINTBOUNDSPROBING, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/tightcontboundsprobing",
         "should continuous bounds be tightened during the probing mode?",
         &propdata->tightcontboundsprobing, TRUE, DEFAULT_TIGHTCONTBOUNDSPROBING, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/createbilinineqs",
         "solve auxiliary LPs in order to find valid inequalities for bilinear terms?",
         &propdata->createbilinineqs, TRUE, DEFAULT_CREATE_BILININEQS, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/orderingalgo",
        "select the type of ordering algorithm which should be used (0: no special ordering, 1: greedy, 2: greedy reverse)",
        &propdata->orderingalgo, TRUE, DEFAULT_ORDERINGALGO, 0, 2, NULL, NULL) );

  SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/separatesol",
         "should the obbt LP solution be separated?",
         &propdata->separatesol, TRUE, DEFAULT_SEPARATESOL, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/sepaminiter",
        "minimum number of iteration spend to separate an obbt LP solution",
        &propdata->sepaminiter, TRUE, DEFAULT_SEPAMINITER, 0, INT_MAX, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/sepamaxiter",
        "maximum number of iteration spend to separate an obbt LP solution",
        &propdata->sepamaxiter, TRUE, DEFAULT_SEPAMAXITER, 0, INT_MAX, NULL, NULL) );

  SCIP_CALL( SCIPaddIntParam(scip, "propagating/" PROP_NAME "/propagatefreq",
        "trigger a propagation round after that many bound tightenings (0: no propagation)",
        &propdata->propagatefreq, TRUE, DEFAULT_PROPAGATEFREQ, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
