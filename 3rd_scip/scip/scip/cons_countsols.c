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

/**@file   cons_countsols.c
 * @brief  constraint handler for counting feasible solutions
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * If this constraint handler is activated than it counts or collects all feasible solutions. We refer to \ref COUNTER for
 * more details about using SCIP for counting feasible solutions.
 *
 * @todo In the last round of presolving we should check if variables exist, which have up and down lock one. In that
 *       case we know that these locks are coming from this constraint handler. Therefore, they are totally free and can
 *       be ignored in the branch and bound process. To get this result we have to store these variables in the
 *       constraint handler data structure (to remember this free dimensions) and fix them to any feasible value.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_and.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_countsols.h"
#include "scip/dialog_default.h"

/* depending if the GMP library is available we use a GMP data type or a SCIP_Longint */
#ifdef WITH_GMP
#include <gmp.h>
typedef mpz_t                Int;
#else
typedef SCIP_Longint         Int;
#endif

/* constraint handler properties */
#define CONSHDLR_NAME          "countsols"
#define CONSHDLR_DESC          "constraint to count feasible solutions"
#define CONSHDLR_ENFOPRIORITY  -9999999 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* default parameter settings */
#define DEFAULT_SPARSETEST         TRUE /**< sparse test on or off */
#define DEFAULT_DISCARDSOLS        TRUE /**< is it allowed to discard solutions */
#define DEFAULT_ACTIVE            FALSE /**< is the constraint handler active */
#define DEFAULT_COLLECT           FALSE /**< should the solutions be collected */
#define DEFAULT_SOLLIMIT           -1LL /**< counting stops, if the given number of solutions were found (-1: no limit) */

/* default column settings */
#define DISP_SOLS_NAME             "sols"
#define DISP_SOLS_DESC             "number of detected feasible solutions"
#define DISP_SOLS_HEADER           " sols "
#define DISP_SOLS_WIDTH            7
#define DISP_SOLS_PRIORITY         110000
#define DISP_SOLS_POSITION         100000
#define DISP_SOLS_STRIPLINE        TRUE

#define DISP_CUTS_NAME             "feasST"
#define DISP_CUTS_DESC             "number of detected non trivial feasible subtrees"
#define DISP_CUTS_HEADER           "feasST"
#define DISP_CUTS_WIDTH            6
#define DISP_CUTS_PRIORITY         110000
#define DISP_CUTS_POSITION         110000
#define DISP_CUTS_STRIPLINE        TRUE

/** creates and adds a constraint which cuts off the solution from the feasibility
 *  region
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sol             : solution to cut off
 *  - conshdlrdata    : constraint handler data
 */
#define CUTOFF_CONSTRAINT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SOL* sol, SCIP_CONSHDLRDATA* conshdlrdata)


/** constraint handler data */
struct SCIP_ConshdlrData
{
   /* solution data and statistic variables */
   SCIP_SPARSESOL**      solutions;          /**< array to store all solutions */
   int                   nsolutions;         /**< number of solution stored */
   int                   ssolutions;         /**< size of the solution array */
   int                   feasST;             /**< number of non trivial feasible subtrees */
   int                   nDiscardSols;       /**< number of discard solutions */
   int                   nNonSparseSols;     /**< number of non sparse solutions */
   Int                   nsols;              /**< number of solutions */
   CUTOFF_CONSTRAINT((*cutoffSolution));     /**< method for cutting of a solution */

   /* constraint handler parameters */
   SCIP_Longint          sollimit;           /**< counting stops, if the given number of solutions were found (-1: no limit) */
   SCIP_Bool             active;             /**< constraint handler active */
   SCIP_Bool             discardsols;        /**< allow to discard solutions */
   SCIP_Bool             sparsetest;         /**< allow to check for sparse solutions */
   SCIP_Bool             collect;            /**< should the solutions be collected */

   SCIP_Bool             warning;            /**< was the warning messages already posted? */

   /* specific problem data */
   SCIP_HASHMAP*         hashmap;            /**< hashmap to store position of active transformed problem variable in our vars array */
   SCIP_VAR**            allvars;            /**< array containing a copy of all variables before presolving */
   SCIP_VAR**            vars;               /**< array containing a copy of all active variables (after presolving) */
   int                   nallvars;           /**< number of all variables in the problem */
   int                   nvars;              /**< number of all active variables in the problem */
   SCIP_Bool             continuous;         /**< are there continuous variables */
};


/*
 * Local methods for handling the <Int> data structure
 */

/** allocates memory for the value pointer */
static
void allocInt(
   Int*                  value               /**< pointer to the value to allocate memory */
   )
{  /*lint --e{715}*/
#ifdef WITH_GMP
   mpz_init(*value);
#endif
}


/** sets the value pointer to the new value */
static
void setInt(
   Int*                  value,              /**< pointer to the value to initialize */
   SCIP_Longint          newvalue            /**< new value */
   )
{
   assert(newvalue < LONG_MAX);

#ifdef WITH_GMP
   mpz_set_si(*value, (long) newvalue);
#else
   (*value) = newvalue;
#endif
}


/** sets a power of 2 to the given value */
static
void setPowerOfTwo(
   Int*                  value,              /**< pointer to the value to increase */
   SCIP_Longint          exponent            /**< exponent for the base 2 */
   )
{
   assert(0 <= exponent && exponent < LONG_MAX);

#ifdef WITH_GMP
   mpz_ui_pow_ui(*value, 2UL, (unsigned long) exponent);
#else
   assert(exponent < 64);
   (*value) = (SCIP_Longint)1 << exponent;
#endif
}


/** free memory */
static
void freeInt(
   Int*                  value               /**< pointer to the value to free */
   )
{  /*lint --e{715}*/
#ifdef WITH_GMP
   mpz_clear(*value);
#endif
}


/** adds one to the given value */
static
void addOne(
   Int*                  value               /**< pointer to the value to increase */
   )
{
#ifdef WITH_GMP
   mpz_add_ui(*value, *value, 1UL);
#else
   (*value)++;
#endif
}


/** adds the summand to the given value */
static
void addInt(
   Int*                  value,              /**< pointer to the value to increase */
   Int*                  summand             /**< summand to add on */
   )
{
#ifdef WITH_GMP
   mpz_add(*value, *value, *summand);
#else
   (*value) += (*summand);
#endif
}


/** multiplies the factor by the given value */
static
void multInt(
   Int*                  value,              /**< pointer to the value to increase */
   SCIP_Longint          factor              /**< factor to multiply with */
   )
{
   assert(0 <= factor && factor < LONG_MAX);

#ifdef WITH_GMP
   mpz_mul_ui(*value, *value, (unsigned long) factor);
#else
   (*value) *= factor;
#endif
}


/** method for creating a string out of an Int which is a mpz_t or SCIP_Longint */
static
void toString(
   Int                   value,              /**< number */
   char**                buffer,             /**< pointer to buffer for storing the string */
   int                   buffersize          /**< length of the buffer */
   )
{  /*lint --e{715}*/
#ifdef WITH_GMP
   (void) mpz_get_str(*buffer, 10, value);
#else
   (void) SCIPsnprintf (*buffer, buffersize, "%" SCIP_LONGINT_FORMAT "", value);
#endif
}


/** method for creating a SCIP_Longing out of an Int */
static
SCIP_Longint getNCountedSols(
   Int                   value,              /**< number to convert */
   SCIP_Bool*            valid               /**< pointer to store if the return value is valid */
   )
{
#ifdef WITH_GMP
   *valid = FALSE;
   if( 0 != mpz_fits_sint_p(value) )
      (*valid) = TRUE;

   return mpz_get_si(value);
#else
   *valid = TRUE;
   return value;
#endif
}


/*
 * Local methods
 */


/** returns whether a given integer variable is unfixed in the local domain */
static
SCIP_Bool varIsUnfixedLocal(
   SCIP_VAR*             var                 /**< integer variable */
   )
{
   assert( var != NULL );
   assert( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS );
   assert( SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var) >= 0.0 );

   return ( SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var) > 0.5 );
}


/** creates the constraint handler data */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store constraint handler data */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, conshdlrdata) );

   (*conshdlrdata)->feasST = 0;
   (*conshdlrdata)->nDiscardSols = 0;
   (*conshdlrdata)->nNonSparseSols = 0;
   (*conshdlrdata)->solutions = NULL;
   (*conshdlrdata)->nsolutions = 0;
   (*conshdlrdata)->ssolutions = 0;

   allocInt(&(*conshdlrdata)->nsols); /*lint !e545*/

   (*conshdlrdata)->cutoffSolution = NULL;
   (*conshdlrdata)->warning = FALSE;
   (*conshdlrdata)->hashmap = NULL;
   (*conshdlrdata)->allvars = NULL;
   (*conshdlrdata)->vars = NULL;
   (*conshdlrdata)->nallvars = 0;
   (*conshdlrdata)->nvars = 0;
   (*conshdlrdata)->continuous = FALSE;

   return SCIP_OKAY;
}


#ifndef NDEBUG
/** check solution in original space */
static
void checkSolutionOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   SCIP_Bool feasible;
   SCIP_RETCODE retcode;

   /* turn off solution counting to be able to check the solution */
   conshdlrdata->active = FALSE;

   SCIPdebugMsg(scip, "check solution in original space before counting\n");

   feasible = FALSE;

   /* check solution in original space */
   retcode = SCIPcheckSolOrig(scip, sol, &feasible, TRUE, TRUE);
   assert(feasible);

   /* check return code manually */
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      SCIPABORT();
   }

   /* turn on solution counting to continue */
   conshdlrdata->active = TRUE;
}
#else
#define checkSolutionOrig(scip, sol, conshdlrdata) /**/
#endif

/** check if the current parameter setting is correct for a save counting process */
static
SCIP_RETCODE checkParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEUR** heuristics;
   int nheuristics;

   int h;
   int intvalue;

   SCIP_Bool valid;

   assert( scip != NULL );

   valid = TRUE;

   /* check if all heuristics are turned off */
   heuristics = SCIPgetHeurs(scip);
   nheuristics = SCIPgetNHeurs(scip);

   for( h = 0; h < nheuristics && valid; ++h )
   {
      if( SCIPheurGetFreq(heuristics[h]) != -1 )
         valid = FALSE;
   }

   if( valid )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "At least one heuristic is not turned off! Heuristic solutions are currently not accepted.\n");
   }

   /* check if restart is turned off */
   SCIP_CALL( SCIPgetIntParam(scip,  "presolving/maxrestarts", &intvalue) );
   if( intvalue != 0 )
   {
      /* need to disabled restarts, since collecting solutions won't work but also the captures for variables are not
       * correctly handled
       */
      SCIPwarningMessage(scip, "counting forces parameter <presolving/maxrestarts> to 0\n");
      if( SCIPisParamFixed(scip, "presolving/maxrestarts") )
      {
         SCIP_CALL( SCIPunfixParam(scip, "presolving/maxrestarts") );
      }

      SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   }

   return SCIP_OKAY;
}

/** creates and adds a constraints which cuts off the current solution from the feasibility region in the case there are
 *  only binary variables
 */
static
CUTOFF_CONSTRAINT(addBinaryCons)
{
   int v;

   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   int nvars;

   SCIP_Real value;
   SCIP_VAR* var;
   SCIP_CONS* cons;

   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );

   vars = conshdlrdata->vars;
   nvars = conshdlrdata->nvars;

   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );

   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];

      assert( var != NULL );
      assert( SCIPvarIsBinary(var) );

      value = SCIPgetSolVal(scip, sol, var);
      assert( SCIPisFeasIntegral(scip, value) );

      if( value > 0.5 )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, var, &consvars[v]) );
      }
      else
         consvars[v] = var;
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, "Setcovering created by countsols", nvars, consvars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* add and release constraint */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* free buffer array */
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


/** creates and adds a bound disjunction constraints which cuts off the current solution from the feasibility region; if
 *  only binary variables are involved, then a set covering constraint is created which is a special case of a bound
 *  disjunction constraint
 */
static
CUTOFF_CONSTRAINT(addIntegerCons)
{
   int v;

   SCIP_VAR** consvars;
   SCIP_VAR** vars;
   SCIP_Real* bounds;
   SCIP_BOUNDTYPE* boundtypes;
   int nvars;
   int nbinvars = 0;
   int nconsvars;
   SCIP_VAR* var;
   SCIP_Real value;

   SCIP_CONS* cons;

   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );

   vars = conshdlrdata->vars;
   nvars = conshdlrdata->nvars;

   nconsvars = nvars * 2;
   assert( nvars > 0 );

   /* allocate buffer memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, nconsvars) );

   nconsvars = 0;

   for( v = nvars - 1; v >= 0; --v )
   {
      var = vars[v];

      assert( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS );

      if( SCIPvarIsBinary(var) )
      {
         ++nbinvars;
         value = SCIPgetSolVal(scip, sol, var);
         assert( SCIPisFeasIntegral(scip, value) );

         if( value < 0.5 )
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars] = 1;
         }
         else
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars] = 0;
         }
      }
      else
      {
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real valueInt;

         assert( SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(var)) );
         assert( SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(var)) );
         assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, var)) );

         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         valueInt = SCIPgetSolVal(scip, sol, var);

         if( SCIPisFeasEQ(scip, valueInt, lb) )
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars] = lb + 1.0;
         }
         else if( SCIPisFeasEQ(scip, valueInt, ub) )
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars] = ub - 1.0;
         }
         else
         {
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars] = valueInt + 1.0;
            consvars[nconsvars] = var;
            ++nconsvars;
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars] = valueInt - 1.0;
         }
      }

      consvars[nconsvars] = var;
      ++nconsvars;
   }

   /* check if only binary variables appear in the constraint; if this is the case we
    * create a set covering constraint instead of a bound disjunction constraint
    */
   if( nvars == nbinvars )
   {
      for( v = nbinvars - 1; v >= 0; --v )
      {
         /* in the case the bound is zero we have use the negated variable */
         if( bounds[v] == 0)
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, consvars[v], &consvars[v]));
         }
      }

      SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, "Setcovering created by countsols", nbinvars, consvars,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, "Bounddisjunction created by countsols",
            nconsvars, consvars, boundtypes, bounds,
            FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   }

   /* add and release constraint locally */
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* free buffer memory */
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &boundtypes);

   return SCIP_OKAY;
}

/** collect given solution or local domains as sparse solution */
static
SCIP_RETCODE collectSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_SOL*             sol                 /**< solution, or NULL if local domains */
   )
{
   SCIP_SPARSESOL* solution;
   SCIP_Longint* lbvalues;
   SCIP_Longint* ubvalues;
   int nvars;
   int v;

   /* ensure size of solution array
    *
    * we use normal memory instead of block memory because this plugin is rarely used and the size of 'solutions'
    * can be arbitrary large and the change that the other blocks can be used is quite small
    */
   if( conshdlrdata->nsolutions == conshdlrdata->ssolutions )
   {
      if( conshdlrdata->ssolutions == 0 )
      {
         conshdlrdata->ssolutions = 100;
         SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrdata->solutions,  conshdlrdata->ssolutions) );
      }
      else
      {
         assert( conshdlrdata->ssolutions < INT_MAX / 2);
         conshdlrdata->ssolutions *= 2;
         SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->solutions,  conshdlrdata->ssolutions) );
      }
   }
   assert( conshdlrdata->nsolutions < conshdlrdata->ssolutions );

   /* get number of active variables */
   nvars = conshdlrdata->nvars;

   SCIPdebugMsg(scip, "creating solution number %d\n", conshdlrdata->nsolutions);

   /* create a solution */
   SCIP_CALL_FINALLY( SCIPsparseSolCreate(&solution, conshdlrdata->vars, nvars, FALSE), SCIPsparseSolFree(&solution) );
   assert(solution != NULL);

   lbvalues = SCIPsparseSolGetLbs(solution);
   ubvalues = SCIPsparseSolGetUbs(solution);
   assert(ubvalues != NULL);
   assert(lbvalues != NULL);

   for( v = nvars - 1; v >= 0; --v )
   {
      SCIP_VAR* var;

      var = conshdlrdata->vars[v];
      assert(var != NULL);

      if( sol == NULL )
      {
         lbvalues[v] = SCIPconvertRealToLongint(scip, SCIPvarGetLbLocal(var));
         ubvalues[v] = SCIPconvertRealToLongint(scip, SCIPvarGetUbLocal(var));
      }
      else
      {
         lbvalues[v] = SCIPconvertRealToLongint(scip, SCIPgetSolVal(scip, sol, var));
         ubvalues[v] = lbvalues[v];
      }

      SCIPdebugMsg(scip, "variable <%s> [%" SCIP_LONGINT_FORMAT ",%" SCIP_LONGINT_FORMAT "]\n",
         SCIPvarGetName(var), lbvalues[v], ubvalues[v]);
   }

   conshdlrdata->solutions[conshdlrdata->nsolutions] = solution;
   conshdlrdata->nsolutions++;

   return SCIP_OKAY;
}


/** counts the number of solutions represented by sol */
static
SCIP_RETCODE countSparseSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Bool             feasible,           /**< is solution feasible? */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_RESULT*          result              /**< pointer to store the result of the checking process */
   )
{
   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );
   assert( result != NULL );

   /* the result should be infeasible since we reject any solution; however, if the solution passes the sparse test the
    * result is set to SCIP_CUTOFF which cuts off the subtree initialized through the current node
    */
   assert(*result == SCIP_INFEASIBLE);

   if( feasible )
   {
      int v;

      Int newsols;

      SCIP_VAR** vars;
      int nvars;

      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;

      SCIPdebugMsg(scip, "counts number of solutions represented through the given one\n");

      /**@note aggregations and multi aggregations: we do not have to care about these things
       *       since we count solution from the transformed problem and therefore, SCIP does
       *       it for us
       */
      assert( SCIPgetNPseudoBranchCands(scip) != 0 );

      allocInt(&newsols); /*lint !e545*/

      /* set newsols to one */
      setInt(&newsols, 1LL); /*lint !e545*/

      if( SCIPgetNBinVars(scip) == SCIPgetNVars(scip) )
      {
         int npseudocands;

         npseudocands = SCIPgetNPseudoBranchCands(scip);

         /* sets a power of 2 to the number of solutions */
         setPowerOfTwo(&newsols, (SCIP_Longint) npseudocands); /*lint !e545*/
      }
      else
      {
         SCIP_VAR* origvar;
         SCIP_Real scalar = 1.0;
         SCIP_Real constant = 0.0;

         SCIP_CALL( SCIPgetPseudoBranchCands(scip, &vars, &nvars, NULL) );

         for( v = 0; v < nvars; ++v )
         {
            var = vars[v];
            origvar = var;
            /* get original variable to decide if we will count the domain; continuous variables aren't counted */
            SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

            if( origvar != NULL && SCIPvarGetType(origvar) != SCIP_VARTYPE_CONTINUOUS )
            {
               lb = SCIPvarGetLbLocal(var);
               ub = SCIPvarGetUbLocal(var);

               SCIPdebugMsg(scip, "variable <%s> Local Bounds are [%g,%g]\n", SCIPvarGetName(var), lb, ub);

               assert( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS );
               assert( SCIPisFeasIntegral(scip, lb) );
               assert( SCIPisFeasIntegral(scip, ub) );
               assert( SCIPisFeasIntegral(scip, ub - lb) );
               assert( SCIPisFeasLT(scip, lb, ub) );

               /* the number of integers laying in the interval [lb,ub] is (ub - lb + 1); to make everything integral we
                * add another 0.5 and cut the fractional part off
                */
               multInt(&newsols, (SCIP_Longint)(ub - lb + 1.5) );  /*lint !e545*/
            }
         }
      }

      *result = SCIP_CUTOFF;
      conshdlrdata->feasST++;

      if( conshdlrdata->collect )
      {
         SCIP_CALL( collectSolution(scip, conshdlrdata, NULL) );
      }

      addInt(&conshdlrdata->nsols, &newsols); /*lint !e545*/
      freeInt(&newsols); /*lint !e545*/
   }
   else if(!conshdlrdata->discardsols)
   {
      SCIP_CALL( conshdlrdata->cutoffSolution(scip, sol, conshdlrdata) );
      addOne(&conshdlrdata->nsols); /*lint !e545*/
      conshdlrdata->nNonSparseSols++;
      if( conshdlrdata->collect )
      {
         SCIP_CALL( collectSolution(scip, conshdlrdata, sol) );
      }
   }
   else
      conshdlrdata->nDiscardSols++;

   return SCIP_OKAY;
}


/** checks if the new solution is feasible for the logicor constraints */
static
SCIP_RETCODE checkLogicor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of enabled constraints */
   SCIP_Bool*            satisfied           /**< pointer to store if the logicor constraints a satisfied */
   )
{
   /**@note the logicor constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree.
    */

   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_Bool fixedone;
   int nvars;
   int c;
   int v;

   SCIPdebugMsg(scip, "check logicor %d constraints\n", nconss);

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"logicor") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );

   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );

   (*satisfied) = TRUE;
   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;

   for( ; c >= 0 && nconss > 0 && (*satisfied); --c )
   {
      SCIPdebugMsg(scip, "logicor constraint %d\n", c);

      if( !SCIPconsIsEnabled(conss[c]) )
         continue;

      nconss--;

      nvars = SCIPgetNVarsLogicor(scip, conss[c]);
      vars = SCIPgetVarsLogicor(scip, conss[c]);

      /* calculate the constraint's activity */
      fixedone = FALSE;
      for( v = 0; v < nvars && !fixedone; ++v )
      {
         assert(SCIPvarIsBinary(vars[v]));

         if( !varIsUnfixedLocal(vars[v] ) )
            fixedone = SCIPvarGetLbLocal(vars[v]) > 0.5;
      }

      if( !fixedone )
      {
         SCIPdebugMsg(scip, "constraint <%s> cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebugPrintCons(scip, conss[c], NULL);
         (*satisfied) = FALSE;
      }
      else
      {
         /* delete constraint from the problem locally since it is satisfied */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** checks if the new solution is feasible for the knapsack constraints */
static
SCIP_RETCODE checkKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of enabled constraints */
   SCIP_Bool*            satisfied           /**< pointer to store if the logicor constraints a satisfied */
   )
{
   /**@note the knapsack constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree.
    */

   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_Longint* weights;
   SCIP_Longint capacity;
   SCIP_Real capa;
   int nvars;
   int c;
   int v;

   SCIPdebugMsg(scip, "check knapsack %d constraints\n", nconss);

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"knapsack") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );

   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );

   (*satisfied) = TRUE;
   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;

   for( ; c >= 0 && nconss > 0 && (*satisfied); --c )
   {
      SCIPdebugMsg(scip, "knapsack constraint %d\n", c);

      if( !SCIPconsIsEnabled(conss[c]) )
         continue;

      nconss--;

      nvars = SCIPgetNVarsKnapsack(scip, conss[c]);
      vars = SCIPgetVarsKnapsack(scip, conss[c]);
      capacity = SCIPgetCapacityKnapsack(scip, conss[c]);
      weights = SCIPgetWeightsKnapsack(scip,conss[c]);

      SCIPdebugMsg(scip, "knapsack capacity = %" SCIP_LONGINT_FORMAT "\n", capacity);

      capa = capacity + 0.1;

      for( v = nvars - 1; v >= 0 && capa >= 0 ; --v )
      {
         SCIPdebug( SCIP_CALL( SCIPprintVar( scip, vars[v], NULL) ) );
         SCIPdebugMsg(scip, "weight = %" SCIP_LONGINT_FORMAT " :\n", weights[v]);
         assert( SCIPvarIsIntegral(vars[v]) );

         /* the weights should be greater or equal to zero */
         assert( weights[v] >= 0);

         if( !varIsUnfixedLocal(vars[v]) )
         {
            /* variables is fixed locally; therefore, subtract fixed variable value multiplied by
             * the weight;
             */
            capa -= weights[v] * SCIPvarGetLbLocal(vars[v]);
         }
         else if( weights[v] >= 1 )
         {
            /*  variable is unfixed and weight is greater than 0; therefore, subtract upper bound
             *  value multiplied by the weight
             */
            capa -= weights[v] * SCIPvarGetUbLocal(vars[v]);
         }
      }

      if( SCIPisFeasLT(scip, capa, 0.0) )
      {
         SCIPdebugMsg(scip, "constraint %s cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebugPrintCons(scip, conss[c], NULL);
         (*satisfied) = FALSE;
      }
      else
      {
         /* delete constraint from the problem locally since it is satisfied */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
      }
   }
   return SCIP_OKAY;
}


/** checks if the new solution is feasible for the bounddisjunction constraints */
static
SCIP_RETCODE checkBounddisjunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of enabled constraints */
   SCIP_Bool*            satisfied           /**< pointer to store if the logicor constraints a satisfied */
   )
{
   /**@note the bounddisjunction constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree
    */

   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_BOUNDTYPE* boundtypes;
   SCIP_Real* bounds;
   SCIP_Bool satisfiedbound;
   int nvars;
   int c;
   int v;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"bounddisjunction") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );

   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );

   (*satisfied) = TRUE;
   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;

   for( ; c >= 0 && nconss > 0 && (*satisfied); --c )
   {
      if( !SCIPconsIsEnabled(conss[c]) )
         continue;

      nconss--;
      satisfiedbound = FALSE;

      nvars = SCIPgetNVarsBounddisjunction(scip, conss[c]);
      vars = SCIPgetVarsBounddisjunction(scip, conss[c]);

      boundtypes = SCIPgetBoundtypesBounddisjunction(scip, conss[c]);
      bounds = SCIPgetBoundsBounddisjunction(scip, conss[c]);

      for( v = nvars-1; v >= 0 && !satisfiedbound; --v )
      {
         SCIPdebug( SCIPprintVar(scip, vars[v], NULL) );

         /* variable should be in right bounds to delete constraint */
         if( boundtypes[v] == SCIP_BOUNDTYPE_LOWER )
            satisfiedbound = SCIPisFeasGE(scip, SCIPvarGetLbLocal(vars[v]), bounds[v]);
         else
         {
            assert( boundtypes[v] == SCIP_BOUNDTYPE_UPPER );
            satisfiedbound = SCIPisFeasLE(scip, SCIPvarGetUbLocal(vars[v]), bounds[v]);
         }
      }

      if( !satisfiedbound )
      {
         SCIPdebugMsg(scip, "constraint %s cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebugPrintCons(scip, conss[c], NULL);
         (*satisfied) = FALSE;
      }
      else
      {
         /* delete constraint from the problem locally since it is satisfied */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
      }
   }
   return SCIP_OKAY;
}


/** checks if the new solution is feasible for the varbound constraints */
static
SCIP_RETCODE checkVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nconss,             /**< number of enabled constraints */
   SCIP_Bool*            satisfied           /**< pointer to store if the logicor constraints a satisfied */
   )
{
   /**@note the varbound constraints are not fully propagated; therefore, we have to check
    *       them by hand if they are satisfied or not; if a constraint is satisfied we
    *       delete it locally from the branch and bound tree.
    */

   SCIP_CONS** conss;
   SCIP_VAR* var;
   SCIP_VAR* vbdvar;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real coef;
   int c;

   SCIPdebugMsg(scip, "check varbound %d constraints\n", nconss);

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr),"varbound") == 0 );
   assert( nconss == SCIPconshdlrGetNEnabledConss(conshdlr) );

   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );

   (*satisfied) = TRUE;
   c = SCIPconshdlrGetNActiveConss(conshdlr) - 1;

   for( ; c >= 0 && nconss > 0 && (*satisfied); --c )
   {
      SCIPdebugMsg(scip, "varbound constraint %d\n", c);

      if( !SCIPconsIsEnabled(conss[c]) )
         continue;

      nconss--;

      var = SCIPgetVarVarbound(scip, conss[c]);
      vbdvar = SCIPgetVbdvarVarbound(scip, conss[c]);

      assert (SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);

      coef = SCIPgetVbdcoefVarbound(scip, conss[c]);
      lhs = SCIPgetLhsVarbound(scip, conss[c]);
      rhs = SCIPgetRhsVarbound(scip, conss[c]);

      /* variables y is fixed locally; therefore, subtract fixed variable value multiplied by
       * the coefficient;
       */
      if(SCIPisGT(scip, SCIPvarGetUbLocal(var), rhs - SCIPvarGetUbLocal(vbdvar) * coef )
         || !SCIPisGE(scip, SCIPvarGetLbLocal(var), lhs - SCIPvarGetLbLocal(vbdvar) * coef ) )
      {
         SCIPdebugMsg(scip, "constraint %s cannot be disabled\n", SCIPconsGetName(conss[c]));
         SCIPdebugPrintCons(scip, conss[c], NULL);
         SCIPdebugMsg(scip, "<%s>  lb: %.15g\t ub: %.15g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         SCIPdebugMsg(scip, "<%s>  lb: %.15g\t ub: %.15g\n", SCIPvarGetName(vbdvar), SCIPvarGetLbLocal(vbdvar), SCIPvarGetUbLocal(vbdvar));
         (*satisfied) = FALSE;
      }
      else
      {
         /* delete constraint from the problem locally since it is satisfied */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
      }
   }

   return SCIP_OKAY;
}


/** check if the current node initializes a non trivial unrestricted subtree */
static
SCIP_RETCODE checkFeasSubtree(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Bool*            feasible            /**< pointer to store the result of the check */
   )
{
   int h;

   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;

   SCIP_CONSHDLR* conshdlr;
   int nconss;

   SCIPdebugMsg(scip, "check if the sparse solution is feasible\n");

   assert( scip != NULL );
   assert( sol != NULL );
   assert( feasible != NULL );

   assert( SCIPgetNPseudoBranchCands(scip) != 0 );

   *feasible = FALSE;

   nconshdlrs = SCIPgetNConshdlrs(scip) - 1;
   conshdlrs = SCIPgetConshdlrs(scip);
   assert(conshdlrs != NULL);

   /* check each constraint handler if there are constraints which are not enabled */
   for( h = nconshdlrs ;  h >= 0 ; --h )
   {
      conshdlr = conshdlrs[h];
      assert( conshdlr != NULL );

      /* skip this constraints handler */
      if( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 )
         continue;

      nconss = SCIPconshdlrGetNEnabledConss(conshdlr);

      if( nconss > 0 )
      {
         SCIP_Bool satisfied;

         SCIPdebugMsg(scip, "constraint handler %s has %d active constraint(s)\n",
            SCIPconshdlrGetName(conshdlr), nconss );

         if( strcmp(SCIPconshdlrGetName(conshdlr), "logicor") == 0 )
         {
            SCIP_CALL( checkLogicor(scip, conshdlr, nconss, &satisfied) );
            if( !satisfied )
            {
               SCIPdebugMsg(scip, "a <logicor> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else if( strcmp(SCIPconshdlrGetName(conshdlr), "knapsack") == 0 )
         {
            SCIP_CALL( checkKnapsack(scip, conshdlr, nconss, &satisfied) );
            if( !satisfied )
            {
               SCIPdebugMsg(scip, "a <knapsack> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else if( strcmp(SCIPconshdlrGetName(conshdlr), "bounddisjunction") == 0 )
         {
            SCIP_CALL( checkBounddisjunction(scip, conshdlr, nconss, &satisfied) );
            if( !satisfied )
            {
               SCIPdebugMsg(scip, "a <bounddisjunction> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else if( strcmp(SCIPconshdlrGetName(conshdlr), "varbound") == 0 )
         {
            SCIP_CALL( checkVarbound(scip, conshdlr, nconss, &satisfied) );
            if( !satisfied )
            {
               SCIPdebugMsg(scip, "a <varbound> constraint cannot be disabled\n");
               return SCIP_OKAY;
            }
         }
         else
         {
            SCIPdebugMsg(scip, "sparse solution is infeasible since the following constraint (and maybe more) is(/are) enabled\n");
            SCIPdebugPrintCons(scip, SCIPconshdlrGetConss(conshdlr)[0], NULL);
            return SCIP_OKAY;
         }
      }
   }

   *feasible = TRUE;
   SCIPdebugMsg(scip, "sparse solution is feasible\n");

   return SCIP_OKAY;
}


/** check the given solution */
static
SCIP_RETCODE checkSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_RESULT*          result              /**< pointer to store the result of the checking process */
   )
{
   SCIP_Longint nsols;
   SCIP_Bool feasible;
   SCIP_Bool valid;

   SCIPdebugMsg(scip, "start to add sparse solution\n");

   assert( scip != NULL );
   assert( sol != NULL );
   assert( conshdlrdata != NULL );
   assert( result != NULL );

   /* the solution should not be found through a heuristic since in this case the informations of SCIP are not valid for
    * this solution
    */

   /**@todo it might be not necessary to check this assert since we can check in general all solutions of feasibility
    *       independently of the origin; however, the locally fixed technique does only work if the solution comes from
    *       the branch and bound tree; in case the solution comes from a heuristic we should try to sequentially fix the
    *       variables in the branch and bound tree and check after every fixing if all constraints are disabled; at the
    *       point where all constraints are disabled the unfixed variables are "stars" (arbitrary);
    */
   assert( SCIPgetNOrigVars(scip) != 0);
   assert( SCIPsolGetHeur(sol) == NULL);

   /* setting result to infeasible since we reject any solution; however, if the solution passes the sparse test or is
    * completely fixed the result is set to SCIP_CUTOFF which cuts off the subtree initialized through the current node
    */
   *result = SCIP_INFEASIBLE;

#ifdef SCIP_DEBUG
   {
      SCIP_VAR* var;
      SCIP_VAR** vars;
      int v;
      int nvars;

      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);

      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         SCIPdebugMsg(scip, "variables <%s> Local Bounds are [%g,%g] Global Bounds are [%g,%g]\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
      }
   }
#endif

   /* check if integer variables are completely fixed */
   if( SCIPgetNPseudoBranchCands(scip) == 0 )
   {
      /* check solution original space */
      checkSolutionOrig(scip, sol, conshdlrdata);

      addOne(&conshdlrdata->nsols); /*lint !e545*/
      conshdlrdata->nNonSparseSols++;

      SCIPdebugMsg(scip, "-> add one to number of solutions\n");

      if( conshdlrdata->collect )
      {
         SCIP_CALL( collectSolution(scip, conshdlrdata, sol) );
      }

      /* in case of continuous variables are present we explicitly cutoff the integer assignment since in case of
       * nonlinear constraint we want to avoid the count that integer assignment again
       */
      if( conshdlrdata->continuous )
      {
         SCIP_CALL( conshdlrdata->cutoffSolution(scip, sol, conshdlrdata) );
      }

      /* since all integer are fixed we cut off the subtree */
      *result = SCIP_CUTOFF;
   }
   else if( conshdlrdata->sparsetest )
   {
      SCIP_CALL( checkFeasSubtree(scip, sol, &feasible) ) ;
      SCIP_CALL( countSparseSol(scip, sol, feasible, conshdlrdata, result) );
   }

   /* transform the current number of solutions into a SCIP_Longint */
   nsols = getNCountedSols(conshdlrdata->nsols, &valid);

   /* check if the solution limit is hit and stop SCIP if this is the case */
   if( conshdlrdata->sollimit > -1 && (!valid || conshdlrdata->sollimit <= nsols) )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   assert( *result == SCIP_INFEASIBLE || *result == SCIP_CUTOFF );
   SCIPdebugMsg(scip, "result is %s\n", *result == SCIP_INFEASIBLE ? "SCIP_INFEASIBLE" : "SCIP_CUTOFF" );

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** creates the handler for countsols constraints and includes it in SCIP */
static
SCIP_RETCODE includeConshdlrCountsols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             dialogs             /**< sould count dialogs be added */
   );

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* in case the count constraint handler is active we avoid copying to ensure a safe count */
   if( conshdlrdata->active )
      *valid = FALSE;
   else
   {
      assert(scip != NULL);
      assert(conshdlr != NULL);
      assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

      /* call inclusion method of constraint handler and do not add counting dialogs */
      SCIP_CALL( includeConshdlrCountsols(scip, FALSE) );

      *valid = TRUE;
   }

   return SCIP_OKAY;
}

#define consCopyCountsols NULL

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free conshdlrdata */
   freeInt(&conshdlrdata->nsols); /*lint !e545*/

   assert( conshdlrdata->solutions == NULL );
   assert( conshdlrdata->nsolutions == 0 );
   assert( conshdlrdata->ssolutions == 0 );

   SCIPfreeBlockMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );

   /* reset counting variables */
   conshdlrdata->feasST = 0;             /* number of non trivial unrestricted subtrees */
   conshdlrdata->nDiscardSols = 0;       /* number of discard solutions */
   conshdlrdata->nNonSparseSols = 0;     /* number of non sparse solutions */
   setInt(&conshdlrdata->nsols, 0LL);    /* number of solutions */  /*lint !e545*/

   conshdlrdata->solutions = NULL;
   conshdlrdata->nsolutions = 0;
   conshdlrdata->ssolutions = 0;

   if( conshdlrdata->active )
   {
      SCIP_VAR** origvars;
      int norigvars;
      int nallvars;
      int v;

      origvars = SCIPgetOrigVars(scip);
      norigvars = SCIPgetNOrigVars(scip);

      /* get number of integral variables */
      conshdlrdata->nallvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->allvars, conshdlrdata->nallvars) );

      nallvars = 0;

      /* capture and lock all variables */
      for( v = 0; v < norigvars; ++v )
      {
         if( SCIPvarGetType(origvars[v]) != SCIP_VARTYPE_CONTINUOUS )
         {
            assert(nallvars < conshdlrdata->nallvars);

            SCIP_CALL( SCIPgetTransformedVar(scip, origvars[v], &conshdlrdata->allvars[nallvars]) );
            assert(conshdlrdata->allvars[nallvars] != NULL);

            /* capture variable to ensure that the variable will not be deleted */
            SCIP_CALL( SCIPcaptureVar(scip, conshdlrdata->allvars[nallvars]) );

            if( strncmp(SCIPvarGetName(conshdlrdata->allvars[nallvars]), "t_andresultant_", strlen("t_andresultant_")) != 0 )
            {
               /* lock variable to avoid dual reductions */
               SCIP_CALL( SCIPaddVarLocks(scip, conshdlrdata->allvars[nallvars], 1, 1) );
            }

            nallvars++;
         }
      }
      assert(nallvars == conshdlrdata->nallvars);

      /* check if continuous variables are present */
      conshdlrdata->continuous = SCIPgetNContVars(scip) > 0;
   }

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int s;
   int v;

   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );

   /* release variables to hashmap */
   for( v = conshdlrdata->nvars - 1; v >= 0; --v )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(conshdlrdata->vars[v])) );
   }

   if( conshdlrdata->hashmap != NULL)
   {
      /* free hashmap of active variables to pistions */
      SCIPhashmapFree(&(conshdlrdata->hashmap));
   }

   /* free active variables */
   SCIPfreeBlockMemoryArrayNull(scip, &(conshdlrdata->vars), conshdlrdata->nvars);
   conshdlrdata->nvars = 0;

   if( conshdlrdata->allvars != NULL )
   {
      /* release and unlock all variables */
      for( v = 0; v < conshdlrdata->nallvars; ++v )
      {
         if( strncmp(SCIPvarGetName(conshdlrdata->allvars[v]), "t_andresultant_", strlen("t_andresultant_")) != 0 )
         {
            /* remove the previously added variable locks */
            SCIP_CALL( SCIPaddVarLocks(scip, conshdlrdata->allvars[v], -1, -1) );
         }

         SCIP_CALL( SCIPreleaseVar(scip, &conshdlrdata->allvars[v]) );
      }

      SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->allvars, conshdlrdata->nallvars);
      conshdlrdata->nallvars = 0;
   }

   if( conshdlrdata->nsolutions > 0 )
   {
      for( s = conshdlrdata->nsolutions - 1; s >= 0 ; --s )
      {
         SCIPsparseSolFree(&(conshdlrdata->solutions[s]));
      }

      SCIPfreeMemoryArrayNull(scip, &conshdlrdata->solutions);
      conshdlrdata->nsolutions = 0;
      conshdlrdata->ssolutions = 0;

      assert( conshdlrdata->solutions == NULL );
   }
   conshdlrdata->continuous = FALSE;

   assert( conshdlrdata->solutions == NULL );
   assert( conshdlrdata->nsolutions == 0 );
   assert( conshdlrdata->ssolutions == 0 );

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The constraint handler may use this call to initialize its branch and bound specific data.
 */
static
SCIP_DECL_CONSINITSOL(consInitsolCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );

   if( conshdlrdata->active )
   {
      SCIP_VAR** vars;
      int v;

      assert(conshdlrdata->nsolutions == 0);
      assert(conshdlrdata->solutions == NULL);

      conshdlrdata->nvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
      vars = SCIPgetVars(scip);

      /* exclude upgrade continuous original variables */
      for( v = conshdlrdata->nvars - 1; v >= 0; --v )
      {
         SCIP_VAR* origvar;
         SCIP_Real scalar = 1.0;
         SCIP_Real constant = 0.0;

         origvar = vars[v];

         /* get original variable to decide if we will count the domain; continuous variables aren't counted */
         SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

         if( origvar != NULL && SCIPvarGetType(origvar) != SCIP_VARTYPE_CONTINUOUS )
            break;
      }
      conshdlrdata->nvars = v + 1;

      /* @todo we need to forbid variable downgrading, from integer type to implicit integer type, e.g. done in
       *       cons_linear
       */
#ifndef NDEBUG
      for( v = conshdlrdata->nvars - 1; v >= 0; --v )
      {
         SCIP_VAR* origvar;
         SCIP_Real scalar = 1.0;
         SCIP_Real constant = 0.0;

         origvar = vars[v];

         /* get original variable to decide if we will count the domain; continuous variables aren't counted */
         SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

         assert(origvar != NULL && SCIPvarGetType(origvar) != SCIP_VARTYPE_CONTINUOUS);
      }
#endif

      /* copy array of active variables */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(conshdlrdata->vars), vars, conshdlrdata->nvars) );

      /* store mapping from all active variables to their position afetr presolving because during solving new variables
       * might be added and therefore could destroy writing collected solutions
       */
      SCIP_CALL( SCIPhashmapCreate(&(conshdlrdata->hashmap), SCIPblkmem(scip), conshdlrdata->nvars + 1) );

      /* add variables to hashmap */
      for( v = conshdlrdata->nvars - 1; v >= 0; --v )
      {
         assert(SCIPvarGetProbindex(conshdlrdata->vars[v]) == v);
         SCIP_CALL( SCIPhashmapInsert(conshdlrdata->hashmap, conshdlrdata->vars[v], (void*) (size_t)(v+1)) );
         SCIP_CALL( SCIPcaptureVar(scip, conshdlrdata->vars[v]) );
      }

      /* check if the problem is binary (ignoring continuous variables) */
      if( SCIPgetNBinVars(scip) == (SCIPgetNVars(scip) - SCIPgetNContVars(scip)) )
         conshdlrdata->cutoffSolution = addBinaryCons;
      else
         conshdlrdata->cutoffSolution = addIntegerCons;
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(nconss == 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL );

   if( conshdlrdata->active && restart )
   {
      SCIPerrorMessage("When collecting and counting solutions restarts need to be disabled (presolving/maxrestarts = 0).\n");
      SCIPABORT();
      return SCIP_INVALIDCALL; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMsg(scip, "method SCIP_DECL_CONSENFOLP(consEnfolpCountsols)\n");

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( nconss == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if( conshdlrdata->active )
   {
      if( !solinfeasible )
      {
         SCIP_SOL* sol;

         SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL ) );

         SCIP_CALL( checkSolution(scip, sol, conshdlrdata, result) );
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }
      else
         *result = SCIP_INFEASIBLE;
   }
   else
      *result = SCIP_FEASIBLE;

   assert( !conshdlrdata->active || *result == SCIP_INFEASIBLE || *result == SCIP_CUTOFF );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxCountsols)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMsg(scip, "method SCIP_DECL_CONSENFORELAX(consEnfolpCountsols)\n");

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( nconss == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if( conshdlrdata->active )
   {
      if( !solinfeasible )
      {
         SCIP_CALL( checkSolution(scip, sol, conshdlrdata, result) );
      }
      else
         *result = SCIP_INFEASIBLE;
   }
   else
      *result = SCIP_FEASIBLE;

   assert( !conshdlrdata->active || *result == SCIP_INFEASIBLE || *result == SCIP_CUTOFF );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCountsols)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMsg(scip, "method SCIP_DECL_CONSENFOPS(consEnfopsCountsols)\n");

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( nconss == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );


   if( conshdlrdata->active )
   {
      if( !solinfeasible )
      {
         SCIP_SOL* sol;

         SCIP_CALL( SCIPcreatePseudoSol(scip, &sol, NULL ) );

         SCIP_CALL( checkSolution(scip, sol, conshdlrdata, result) );
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      }
      else
         *result = SCIP_INFEASIBLE;
   }
   else
      *result = SCIP_FEASIBLE;

   assert( !conshdlrdata->active || *result == SCIP_INFEASIBLE || *result == SCIP_CUTOFF );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckCountsols)
{  /*lint --e{715}*/
   /**@todo solutions which come from scip_check should be ignored since it is not clear who
    *       generated these solution; later we should analyze this problem */
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMsg(scip, "method SCIP_DECL_CONSCHECK(consCheckCountsols)\n");

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if( conshdlrdata->active )
   {
      if( !conshdlrdata->warning )
      {
         SCIPwarningMessage(scip, "a solution comes in over <SCIP_DECL_CONSCHECK(consCheckCountsols)>; currently these solutions are ignored\n");
         conshdlrdata->warning = TRUE;
      }

      *result = SCIP_INFEASIBLE;
   }
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockCountsols)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/*
 * Callback methods and local method for dialogs
 */

/** dialog execution method for the count command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCountPresolve)
{  /*lint --e{715}*/
   SCIP_Bool active;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );
   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", &active) );

   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
      /* activate constraint handler cons_countsols */
      if( !active )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", TRUE) );
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* presolve problem */
      SCIP_CALL( SCIPpresolve(scip) );

      /* reset cons_countsols activation */
      if( !active )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", FALSE) );
      }
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPdialogMessage(scip, NULL, "problem is already presolved\n");
      break;

   case SCIP_STAGE_SOLVED:
      SCIPdialogMessage(scip, NULL, "problem is already (pre)solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   } /*lint --e{616}*/

   SCIPdialogMessage(scip, NULL, "\n");
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the count command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCount)
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;
   SCIP_Bool active;

   SCIP_Bool valid;
   SCIP_Longint nsols;
   int displayprimalbound;
   int displaygap;
   int displaysols;
   int displayfeasST;
   int nrestarts;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );
   SCIPdialogMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", &active) );
   SCIP_CALL( SCIPgetIntParam(scip, "presolving/maxrestarts", &nrestarts) );

   if( nrestarts != 0 )
   {
      /* need to disabled restarts, since collecting solutions won't work but also the captures for variables are not
       * correctly handled
       */
      SCIPwarningMessage(scip, "counting forces parameter <presolving/maxrestarts> to 0\n");
      if( SCIPisParamFixed(scip, "presolving/maxrestarts") )
      {
         SCIP_CALL( SCIPunfixParam(scip, "presolving/maxrestarts") );
      }
      SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   }

   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
      /* activate constraint handler cons_countsols */
      if( !active )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", TRUE) );
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* presolve problem */
      SCIP_CALL( SCIPpresolve(scip) );
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
      /* reset activity status of constraint handler cons_countsols */
      if( !active )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", FALSE) );
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      /* check if the problem contains continuous variables */
      if( SCIPgetNContVars(scip) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "Problem contains continuous variables (after presolving). Counting projection to integral variables!\n");
      }

      /* turn off primal bound and gap column */
      SCIP_CALL( SCIPgetIntParam(scip, "display/primalbound/active", &displayprimalbound) );
      if( displayprimalbound != 0 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/primalbound/active", 0) );
      }
      SCIP_CALL( SCIPgetIntParam(scip, "display/gap/active", &displaygap) );
      if( displaygap != 0 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/gap/active", 0) );
      }

      /* turn on sols and feasST column */
      SCIP_CALL( SCIPgetIntParam(scip, "display/sols/active", &displaysols) );
      if( displayprimalbound != 2 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/sols/active", 2) );
      }
      SCIP_CALL( SCIPgetIntParam(scip, "display/feasST/active", &displayfeasST) );
      if( displayprimalbound != 2 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/feasST/active", 2) );
      }

      /* find the countsols constraint handler */
      assert( SCIPfindConshdlr(scip, CONSHDLR_NAME) != NULL );

      retcode =  SCIPcount(scip);

      valid = FALSE;
      nsols = SCIPgetNCountedSols(scip, &valid);

      if( valid )
         SCIPdialogMessage(scip, NULL, "Feasible Solutions : %" SCIP_LONGINT_FORMAT "", nsols);
      else
      {
         char* buffer;
         int buffersize = SCIP_MAXSTRLEN;
         int requiredsize;

         SCIP_CALL( SCIPallocBufferArray(scip, &buffer, buffersize) );
         SCIPgetNCountedSolsstr(scip, &buffer, buffersize, &requiredsize);

         if( requiredsize > buffersize )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &buffer, requiredsize) );
            SCIPgetNCountedSolsstr(scip, &buffer, buffersize, &requiredsize);
         }

         assert( buffersize >= requiredsize );
         SCIPdialogMessage(scip, NULL, "Feasible Solutions : %s", buffer);

         SCIPfreeBufferArray(scip, &buffer);
      }

      SCIPdialogMessage(scip, NULL, " (%d non-trivial feasible subtrees)\n", SCIPgetNCountedFeasSubtrees(scip));

      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

      /* reset display columns */
      if( displayprimalbound != 0 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/primalbound/active", displayprimalbound) );
      }
      if( displaygap != 0 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/gap/active", displaygap) );
      }

      /* reset sols and feasST column */
      if( displaysols != 2 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/sols/active", displaysols) );
      }
      if( displayfeasST != 2 )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "display/feasST/active", displayfeasST) );
      }

      /* reset cons_countsols activation */
      if( !active )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", FALSE) );
      }

      /* evaluate retcode */
      SCIP_CALL( retcode );
      break;

   case SCIP_STAGE_SOLVED:
      SCIPdialogMessage(scip, NULL, "problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   } /*lint --e{616}*/

   SCIPdialogMessage(scip, NULL, "\n");
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** comparison method for sorting variables by non-decreasing w.r.t. problem index */
static
SCIP_DECL_SORTPTRCOMP(varCompProbindex)
{
   SCIP_VAR* var1;
   SCIP_VAR* var2;

   var1 = (SCIP_VAR*)elem1;
   var2 = (SCIP_VAR*)elem2;

   assert(var1 != NULL);
   assert(var2 != NULL);

   if( SCIPvarGetProbindex(var1) < SCIPvarGetProbindex(var2) )
      return -1;
   else if( SCIPvarGetProbindex(var1) > SCIPvarGetProbindex(var2) )
      return +1;
   else
   {
      assert(var1 == var2 || (SCIPvarGetProbindex(var1) == -1 && SCIPvarGetProbindex(var2) == -1));
      return 0;
   }
}

/** expands the sparse solutions and writes them to the file */
static
SCIP_RETCODE writeExpandedSolutions(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file handler */
   SCIP_VAR**            allvars,            /**< SCIP variables */
   int                   nallvars,           /**< number of all variables */
   SCIP_VAR**            activevars,         /**< SCIP variables */
   int                   nactivevars,        /**< number of active variables */
   SCIP_HASHMAP*         hashmap,            /**< hashmap from active solution variable to the position in the active
                                              *   variables array
                                              */
   SCIP_SPARSESOL**      sols,               /**< sparse solutions to expands and write */
   int                   nsols               /**< number of sparse solutions */
   )
{
   SCIP_SPARSESOL* sparsesol;
   SCIP_VAR** vars;
   SCIP_Real* scalars;
   SCIP_Longint* sol;
   SCIP_Longint solcnt;
   int s;
   int v;

   assert(scip != NULL);
   assert(file != NULL);
   assert(hashmap != NULL);
   assert(allvars != NULL || nallvars == 0);
   assert(activevars != NULL || nactivevars == 0);
   assert(sols != NULL || nsols == 0);

   solcnt = 0;

   /* get memory to store active solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &sol, nactivevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nactivevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scalars, nactivevars) );

   /* loop over all sparse solutions */
   for( s = 0; s < nsols; ++s )
   {
      sparsesol = sols[s]; /*lint !e613*/
      assert(sparsesol != NULL);
      assert(SCIPsparseSolGetNVars(sparsesol) == nactivevars);

      /* get first solution of the sparse solution */
      SCIPsparseSolGetFirstSol(sparsesol, sol, nactivevars);

      do
      {
         SCIP_Real objval;

         solcnt++;

         /* print solution number */
         SCIPinfoMessage(scip, file, "%d(%" SCIP_LONGINT_FORMAT "), ", s+1, solcnt);

         objval = 0.0;

         /* write none active variables */
         for( v = 0; v < nallvars; ++v )
         {
            SCIP_Real constant;
            SCIP_Real realvalue;
            int requiredsize;
            int nvars;
            int idx;
            int i;

            vars[0] = allvars[v]; /*lint !e613*/
            scalars[0] = 1.0;
            nvars = 1;
            constant = 0.0;

            SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, &nvars, nallvars, &constant, &requiredsize, TRUE) );
            assert(requiredsize <= nallvars);
            assert(nvars <= nactivevars);

            realvalue = constant;

            for( i = 0; i < nvars; ++i )
            {
               assert(SCIPhashmapExists(hashmap, vars[i]));
               idx = ((int) (size_t)SCIPhashmapGetImage(hashmap, vars[i])) - 1;
               assert(0 <= idx && idx < nactivevars);
               assert(activevars[idx] == vars[i]); /*lint !e613*/

               objval += SCIPvarGetObj(vars[i]) * sol[idx];
               realvalue += scalars[i] * sol[idx];
            }
            assert(SCIPisIntegral(scip, realvalue));

            SCIPinfoMessage(scip, file, "%g, ", realvalue);
         }

         /* transform objective value into original problem space */
         objval = SCIPretransformObj(scip, objval);

         /* output the objective value of the solution */
         SCIPinfoMessage(scip, file, "%g\n", objval);
      }
      while( SCIPsparseSolGetNextSol(sparsesol, sol, nactivevars) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &scalars);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &sol);

   return SCIP_OKAY;
}

/** execution method of dialog for writing all solutions */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteAllsolutions)
{  /*lint --e{715}*/
   FILE* file;
   SCIP_Longint nsols;
   char* filename;
   char* word;
   SCIP_Bool endoffile;
   SCIP_Bool valid;

   assert( scip != NULL );

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "no problem available\n");
      break;
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_FREETRANS:
      SCIPdialogMessage(scip, NULL, "the counting process was not started yet\n");
      break;
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_EXITSOLVE:
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;
      int nsparsesols;

      valid = FALSE;
      nsols = SCIPgetNCountedSols(scip, &valid);

      /* find the countsols constraint handler */
      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
      assert( conshdlr != NULL );

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      nsparsesols = conshdlrdata->nsolutions;

      if( !valid )
      {
         /* too many solutions, output not "possible" */
         char* buffer;
         int buffersize;
         int requiredsize;

         buffersize = SCIP_MAXSTRLEN;

         SCIP_CALL( SCIPallocBufferArray(scip, &buffer, buffersize) );
         SCIPgetNCountedSolsstr(scip, &buffer, buffersize, &requiredsize);

         if( requiredsize > buffersize )
         {
            buffersize = requiredsize;
            SCIP_CALL( SCIPreallocBufferArray(scip, &buffer, requiredsize) );
            SCIPgetNCountedSolsstr(scip, &buffer, buffersize, &requiredsize);
         }

         assert( buffersize >= requiredsize );
         SCIPdialogMessage(scip, NULL, "no output, because of too many feasible solutions : %s\n", buffer);

         SCIPfreeBufferArray(scip, &buffer);
      }
      else if( nsols == 0 )
      {
         SCIPdialogMessage(scip, NULL, "there are no counted solutions\n");
      }
      else if( nsparsesols == 0 )
      {
         SCIPdialogMessage(scip, NULL, "there is no solution collect (set parameter <constraints/countsols/collect> to TRUE)\n");
      }
      else
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &word, &endoffile) );

         /* copy the filename for later use */
         SCIP_CALL( SCIPduplicateBufferArray(scip, &filename, word, (int)strlen(word)+1) );

         if( endoffile )
         {
            *nextdialog = NULL;
            return SCIP_OKAY;
         }

         SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

         if( filename[0] != '\0' )
         {
            file = fopen(filename, "w");

            if( file == NULL )
            {
               SCIPdialogMessage(scip, NULL, "error creating file <%s>\n", filename);
               SCIPdialoghdlrClearBuffer(dialoghdlr);
            }
            else
            {
               SCIP_SPARSESOL** sparsesols;
               SCIP_VAR** origvars;
               SCIP_VAR** allvars;
               int norigvars;
               int nvars;
               int v;

               SCIP_RETCODE retcode;

               /* get sparse solutions defined over the active variables */
               nvars = conshdlrdata->nvars;
               sparsesols = conshdlrdata->solutions;

               /* get original problem variables */
               retcode = SCIPallocBufferArray(scip, &origvars, SCIPgetNOrigVars(scip));
               if( retcode != SCIP_OKAY )
               {
                  fclose(file);
                  SCIP_CALL( retcode );
               }

               norigvars = 0;

               for( v = 0; v < SCIPgetNOrigVars(scip); ++v )
               {
                  if( SCIPvarGetType(SCIPgetOrigVars(scip)[v]) != SCIP_VARTYPE_CONTINUOUS )
                  {
                     origvars[norigvars] = SCIPgetOrigVars(scip)[v];
                     norigvars++;
                  }
               }
               assert(norigvars == conshdlrdata->nallvars);

               retcode = SCIPduplicateBufferArray(scip, &allvars, conshdlrdata->allvars, norigvars);
               if( retcode != SCIP_OKAY )
               {
                  fclose(file); /*lint !e449*/
                  SCIP_CALL( retcode );
               }

               /* sort original variables array and the corresponding transformed variables w.r.t. the problem index */
               SCIPsortDownPtrPtr((void**)allvars, (void**)origvars, varCompProbindex, norigvars);

               SCIPdialogMessage(scip, NULL, "saving %" SCIP_LONGINT_FORMAT " (%d) feasible solutions\n", nsols, nsparsesols);

               /* first row: output the names of the variables in the given ordering */
               SCIPinfoMessage(scip, file, "#, ");

               for( v = 0; v < norigvars; ++v )
               {
#ifndef NDEBUG
                  {
                     /* check if the original variable fits to the transform variable the constraint handler has */
                     SCIP_VAR* transvar;
                     SCIP_CALL( SCIPgetTransformedVar(scip, origvars[v], &transvar) );
                     assert(transvar != NULL);
                     assert(transvar == allvars[v]);
                  }
#endif
                  SCIPinfoMessage(scip, file, "%s, ", SCIPvarGetName(origvars[v]));
               }

               SCIPinfoMessage(scip, file, "objval\n");

               /* expand and write solution */
               retcode = writeExpandedSolutions(scip, file, allvars, conshdlrdata->nallvars, conshdlrdata->vars, nvars, conshdlrdata->hashmap, sparsesols, nsparsesols);
               if( retcode != SCIP_OKAY )
               {
                  fclose(file);
                  SCIP_CALL( retcode );
               }
               SCIPdialogMessage(scip, NULL, "written solutions information to file <%s>\n", filename);

               SCIPfreeBufferArray(scip, &allvars);
               SCIPfreeBufferArray(scip, &origvars);

               fclose(file);
            }

            /* free buffer array */
            SCIPfreeBufferArray(scip, &filename);
         }
      }
      break;
   }
   case SCIP_STAGE_FREE:
      SCIPerrorMessage("invalid call during SCIP_STAGE_FREE\n");
      return SCIP_ERROR;
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** create the interactive shell dialogs for the counting process */
static
SCIP_RETCODE createCountDialog(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DIALOG* root;
   SCIP_DIALOG* dialog;
   SCIP_DIALOG* setmenu;
   SCIP_DIALOG* submenu;

   /* includes or updates the default dialog menus in SCIP */
   SCIP_CALL( SCIPincludeDialogDefault(scip) );

   root = SCIPgetRootDialog(scip);
   assert( root != NULL );

   /* add dialog entry for counting */
   if( !SCIPdialogHasEntry(root, "count") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecCount, NULL, NULL,
            "count", "count number of feasible solutions", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* add dialog entry for counting */
   if( !SCIPdialogHasEntry(root, "countpresolve") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecCountPresolve, NULL, NULL,
            "countpresolve", "presolve instance before counting number of feasible solutions", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* search for the "write" sub menu to add "allsolutions" dialog */
   if( SCIPdialogFindEntry(root, "write", &submenu) != 1 )
   {
      SCIPerrorMessage("write sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert(submenu != NULL);

   /* add dialog "allsolutions" to sub menu "write" */
   if( !SCIPdialogHasEntry(submenu, "allsolutions") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, SCIPdialogExecWriteAllsolutions, NULL, NULL,
            "allsolutions", "write all counted primal solutions to file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* search for the "set" sub menu to find the "emphasis" sub menu */
   if( SCIPdialogFindEntry(root, "set", &setmenu) != 1 )
   {
      SCIPerrorMessage("set sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert(setmenu != NULL);

   return SCIP_OKAY;
}

/*
 * Callback methods for columns
 */

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputSols)
{  /*lint --e{715}*/
#ifndef NDEBUG
   SCIP_CONSHDLR* conshdlr;
#endif
   SCIP_Longint sols;
   SCIP_Bool valid;

   assert(disp != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_SOLS_NAME) == 0);
   assert(scip != NULL);

#ifndef NDEBUG
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   assert( SCIPconshdlrGetNConss(conshdlr) == 0 );
#endif

   sols = SCIPgetNCountedSols(scip, &valid);

   if( !valid )
   {
      SCIPinfoMessage(scip, file, "TooMany");
   }
   else
   {
      SCIPdispLongint(SCIPgetMessagehdlr(scip), file, sols, DISP_SOLS_WIDTH);
   }

   return SCIP_OKAY;
}


/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputFeasSubtrees)
{ /*lint --e{715}*/
#ifndef NDEBUG
   SCIP_CONSHDLR* conshdlr;
#endif

   assert(disp != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPdispGetName(disp), DISP_CUTS_NAME) == 0);

#ifndef NDEBUG
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   assert( SCIPconshdlrGetNConss(conshdlr) == 0 );
#endif

   SCIPdispLongint(SCIPgetMessagehdlr(scip), file, SCIPgetNCountedFeasSubtrees(scip), DISP_CUTS_WIDTH);

   return SCIP_OKAY;
}


/*
 * Interface methods of constraint handler
 */

/** creates the handler for countsols constraints and includes it in SCIP */
static
SCIP_RETCODE includeConshdlrCountsols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             dialogs             /**< sould count dialogs be added */
   )
{
   /* create countsol constraint handler data */
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

#ifdef WITH_GMP
   char gmpversion[20];
#endif

   /* create constraint handler specific data here */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpCountsols, consEnfopsCountsols, consCheckCountsols, consLockCountsols,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyCountsols, consCopyCountsols) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitCountsols) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolCountsols) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeCountsols) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitCountsols) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolCountsols) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxCountsols) );

   /* add countsols constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/active",
         "is the constraint handler active?",
         &conshdlrdata->active, FALSE, DEFAULT_ACTIVE, NULL, NULL));
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/sparsetest",
         "should the sparse solution test be turned on?",
         &conshdlrdata->sparsetest, FALSE, DEFAULT_SPARSETEST, NULL, NULL));
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/discardsols",
         "is it allowed to discard solutions?",
         &conshdlrdata->discardsols, FALSE, DEFAULT_DISCARDSOLS, NULL, NULL));
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/collect",
         "should the solutions be collected?",
         &conshdlrdata->collect, FALSE, DEFAULT_COLLECT, NULL, NULL));
   SCIP_CALL( SCIPaddLongintParam(scip,
         "constraints/" CONSHDLR_NAME "/sollimit",
         "counting stops, if the given number of solutions were found (-1: no limit)",
         &conshdlrdata->sollimit, FALSE, DEFAULT_SOLLIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL));

   /* create the interactive shell dialogs for the counting process  */
   if( dialogs )
   {
      SCIP_CALL( createCountDialog(scip) );
   }

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_SOLS_NAME, DISP_SOLS_DESC, DISP_SOLS_HEADER, SCIP_DISPSTATUS_OFF,
         NULL, NULL, NULL, NULL, NULL, NULL, dispOutputSols,
         NULL, DISP_SOLS_WIDTH, DISP_SOLS_PRIORITY, DISP_SOLS_POSITION, DISP_SOLS_STRIPLINE) );
   SCIP_CALL( SCIPincludeDisp(scip, DISP_CUTS_NAME, DISP_CUTS_DESC, DISP_CUTS_HEADER, SCIP_DISPSTATUS_OFF,
         NULL, NULL, NULL, NULL, NULL, NULL, dispOutputFeasSubtrees,
         NULL, DISP_CUTS_WIDTH, DISP_CUTS_PRIORITY, DISP_CUTS_POSITION, DISP_CUTS_STRIPLINE) );

#ifdef WITH_GMP
#ifdef mpir_version
   /* add info about using MPIR to external codes information */
   (void) SCIPsnprintf(gmpversion, (int) sizeof(gmpversion), "MPIR %s", mpir_version);
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, gmpversion, "Multiple Precision Integers and Rationals Library developed by W. Hart (mpir.org)") );
#else
   /* add info about using GMP to external codes information */
   (void) SCIPsnprintf(gmpversion, (int) sizeof(gmpversion), "GMP %s", gmp_version);
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, gmpversion, "GNU Multiple Precision Arithmetic Library developed by T. Granlund (gmplib.org)") );
#endif
#endif

   return SCIP_OKAY;
}

/** creates the handler for countsols constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCountsols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include constraint handler including the count dialog */
   SCIP_CALL( includeConshdlrCountsols(scip, TRUE) );

   return SCIP_OKAY;
}


/** execute counting */
SCIP_RETCODE SCIPcount(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool active;

   /* activate constraint handler cons_countsols */
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", &active) );
   if( !active )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", TRUE) );
   }

   /* check if the parameter setting allows a valid counting process */
   SCIP_CALL( checkParameters(scip) );

   /* start the solving process */
   SCIP_CALL( SCIPsolve(scip) );

   /* reset activity status of constraint handler cons_countsols */
   if( !active )
   {
      SCIP_CALL( SCIPsetBoolParam(scip, "constraints/" CONSHDLR_NAME "/active", FALSE) );
   }

   return SCIP_OKAY;
}


/** returns number of feasible solutions found as SCIP_Longint; if the number does not fit into
 *  a SCIP_Longint the valid flag is set to FALSE
 */
SCIP_Longint SCIPgetNCountedSols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            valid               /**< pointer to store if the return value is valid */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* find the countsols constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   return getNCountedSols(conshdlrdata->nsols, valid);
}


/** puts the number of counted solutions in the given char* buffer */
void SCIPgetNCountedSolsstr(
   SCIP*                 scip,               /**< SCIP data structure */
   char**                buffer,             /**< buffer to store the number for counted solutions */
   int                   buffersize,         /**< buffer size */
   int*                  requiredsize        /**< pointer to store the required size */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* find the countsols constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

#ifdef WITH_GMP
   /* size must be by two larger than the length of the string, since there need to be storage for a sign and a
    * null-termination
    */
   assert(0 <= (int) (mpz_sizeinbase( conshdlrdata->nsols, 10 ) + 2));
   *requiredsize = (int) (mpz_sizeinbase( conshdlrdata->nsols, 10 ) + 2);
   if( *requiredsize <= buffersize)
      toString(conshdlrdata->nsols, buffer, buffersize);
#else
   if( conshdlrdata->nsols < pow(10.0, (double)buffersize) )
   {
      toString(conshdlrdata->nsols, buffer, buffersize);
      *requiredsize = (int)strlen(*buffer);
   }
   else
      *requiredsize = 21;
#endif
}


/** returns number of counted non trivial feasible subtrees */
SCIP_Longint SCIPgetNCountedFeasSubtrees(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );

   /* find the countsols constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   return conshdlrdata->feasST;
}


/** Method to get the sparse solution.
 *
 *  @note You get the pointer to the sparse solutions stored in the constraint handler (not a copy).
 *
 *  @note The sparse solutions are stored w.r.t. the active variables. This are the variables which got not removed
 *        during presolving. For none active variables the value has to be computed depending on their aggregation
 *        type. See for more details about that \ref COLLECTALLFEASEBLES.
 */
void SCIPgetCountedSparseSols(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to active variable array defining to variable order */
   int*                  nvars,              /**< number of active variables */
   SCIP_SPARSESOL***     sols,               /**< pointer to the solutions */
   int*                  nsols               /**< pointer to number of solutions */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );

   /* find the countsols constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *vars = conshdlrdata->vars;
   *nvars = conshdlrdata->nvars;
   *sols = conshdlrdata->solutions;
   *nsols = conshdlrdata->nsolutions;
}

/** setting SCIP parameters for such that a valid counting process is possible */
SCIP_RETCODE SCIPsetParamsCountsols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPsetEmphasis(scip, SCIP_PARAMEMPHASIS_COUNTER, TRUE) );
   return SCIP_OKAY;
}
