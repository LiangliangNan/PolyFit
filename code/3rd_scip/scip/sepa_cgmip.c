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
/* #define SCIP_WRITEPROB */
/* #define SCIP_OUTPUT */
/**@file   sepa_cgmip.c
 * @brief  Chvatal-Gomory cuts computed via a sub-MIP
 * @author Marc Pfetsch
 *
 * Separate Chv&aacute;tal-Gomory cuts using a sub-MIP. The approach is based on the following papers.
 *
 * M. Fischetti and A. Lodi@n
 * Optimizing over the first Chv&aacute;tal closure,@n
 * in: M. J&uuml;nger and V. Kaibel (eds.) Integer Programming and Combinatorial Optimization IPCO 2005,@n
 * LNCS 3509, pp. 12-22. Springer, Berlin Heidelberg New York (2005)
 *
 * M. Fischetti and A. Lodi@n
 * Optimizing over the first Chv&aacute;tal closure,@n
 * Mathematical Programming 110, 3-20 (2007)
 *
 * P. Bonami, G. Cornu&eacute;jols, S. Dash, M. Fischetti, and A. Lodi@n
 * Projected Chv&aacute;tal-Gomory cuts for mixed integer linear programs,@n
 * Mathematical Programming 113, No. 2 (2008)
 *
 *
 * There are several versions to generate the final cut:
 *
 * - The CMIR-routines of SCIP can be used (if @p usecmir is true). One can determine which bound is
 *   used in the rounding operation (if cmirownbounds is true) or let SCIP choose the best. This
 *   version is generally numerically the most stable.
 * - If @p usestrongcg is true, we try to generate Strong-CG cuts (as done in sepa_strongcg.c).
 * - One can directly generate the CG-cut as computed (if @p usecmir and @p usestrongcg are
 *   false). The cut is not take from the solution of the MIP, but is recomputed, and some care (but
 *   not as much as in the first version) has been taken to create a valid cut.
 *
 * The computation time of the separation MIP is limited as follows:
 * - There is a node limit (parameters @a minnodelimit and @a maxnodelimit).
 * - There is a time limit (parameter @a timelimit).
 * - If paramter @a earlyterm is true, the separation is run until the first cut that is violated is
 *   found. (Note that these cuts are not necessarily added to the LP, because here also the norm of
 *   the cuts are taken into account - which cannot easily be included into the separation subscip.)
 *   Then the solution is continued for a certain number of nodes.
 *
 * @todo Check whether one can weaken the conditions on the continuous variables.
 * @todo Use pointers to originating separators to sort out cuts that should not be used.
 *
 * @warning This separator should be used carefully - it may require a long separation time.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_cgmip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/pub_misc.h"
#include "scip/pub_lp.h"


#define SEPA_NAME              "cgmip"
#define SEPA_DESC              "Chvatal-Gomory cuts via MIPs separator"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP           TRUE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        50 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXDEPTH             -1 /**< maximal depth at which the separator is applied */
#define DEFAULT_DECISIONTREE      FALSE /**< Use decision tree to turn separation on/off? */
#define DEFAULT_TIMELIMIT          1e20 /**< time limit for sub-MIP (set to infinity in order to be deterministic) */
#define DEFAULT_MEMORYLIMIT        1e20 /**< memory limit for sub-MIP */
#define DEFAULT_CUTCOEFBND       1000.0 /**< bounds on the values of the coefficients in the CG-cut */
#define DEFAULT_MINNODELIMIT      500LL /**< minimum number of nodes considered for sub-MIP (-1: unlimited) */
#define DEFAULT_MAXNODELIMIT     5000LL /**< maximum number of nodes considered for sub-MIP (-1: unlimited) */
#define DEFAULT_ONLYACTIVEROWS    FALSE /**< Use only active rows to generate cuts? */
#define DEFAULT_MAXROWAGE            -1 /**< maximal age of rows to consider if onlyactiverows is false */
#define DEFAULT_ONLYRANKONE       FALSE /**< Separate rank 1 inequalities w.r.t. CG-MIP separator? */
#define DEFAULT_ONLYINTVARS       FALSE /**< Generate cuts for problems with only integer variables? */
#define DEFAULT_CONTCONVERT       FALSE /**< Convert some integral variables to be continuous to reduce the size of the sub-MIP? */
#define DEFAULT_CONTCONVFRAC        0.1 /**< fraction of integral variables converted to be continuous (if contconvert) */
#define DEFAULT_CONTCONVMIN         100 /**< minimum number of integral variables before some are converted to be continuous */
#define DEFAULT_INTCONVERT        FALSE /**< Convert some integral variables attaining fractional values to have integral value? */
#define DEFAULT_INTCONVFRAC         0.1 /**< fraction of fractional integral variables converted to have integral value (if intconvert) */
#define DEFAULT_INTCONVMIN          100 /**< minimum number of integral variables before some are converted to have integral value */
#define DEFAULT_SKIPMULTBOUNDS     TRUE /**< Skip the upper bounds on the multipliers in the sub-MIP? */
#define DEFAULT_OBJLONE           FALSE /**< Should the objective of the sub-MIP only minimize the l1-norm of the multipliers? */
#define DEFAULT_OBJWEIGHT         1e-03 /**< objective weight for artificial variables */
#define DEFAULT_OBJWEIGHTSIZE      TRUE /**< Weight each row by its size? */
#define DEFAULT_DYNAMICCUTS        TRUE /**< Should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_USECMIR            TRUE /**< Use CMIR-generator (otherwise add cut directly)? */
#define DEFAULT_USESTRONGCG       FALSE /**< Use strong CG-function to strengthen cut? */
#define DEFAULT_CMIROWNBOUNDS     FALSE /**< Tell CMIR-generator which bounds to used in rounding? */
#define DEFAULT_USECUTPOOL         TRUE /**< Use cutpool to store CG-cuts even if the are not efficient? */
#define DEFAULT_PRIMALSEPARATION   TRUE /**< Only separate cuts that are tight for the best feasible solution? */
#define DEFAULT_EARLYTERM          TRUE /**< Terminate separation if a violated (but possibly sub-optimal) cut has been found? */
#define DEFAULT_ADDVIOLATIONCONS  FALSE /**< Add constraint to subscip that only allows violated cuts (otherwise add obj. limit)?*/
#define DEFAULT_ADDVIOLCONSHDLR   FALSE /**< Add constraint handler to filter out violated cuts? */
#define DEFAULT_CONSHDLRUSENORM    TRUE /**< Should the violation constraint handler use the norm of a cut to check for feasibility? */
#define DEFAULT_USEOBJUB          FALSE /**< Use upper bound on objective function (via primal solution)? */
#define DEFAULT_USEOBJLB          FALSE /**< Use lower bound on objective function (via lower bound)? */
#define DEFAULT_SUBSCIPFAST        TRUE /**< Should the settings for the sub-MIP be optimized for speed? */
#define DEFAULT_OUTPUT            FALSE /**< Should information about the sub-MIP and cuts be displayed? */
#define DEFAULT_RANDSEED            101 /**< start random seed for random number generation */

#define NROWSTOOSMALL                 5 /**< only separate if the number of rows is larger than this number */
#define NCOLSTOOSMALL                 5 /**< only separate if the number of columns is larger than this number */

#define EPSILONVALUE              1e-03 /**< epsilon value needed to model strict-inequalities */
#define BETAEPSILONVALUE          1e-02 /**< epsilon value for fracbeta - is larger than EPSILONVALUE for numerical stability */
#define STALLNODELIMIT           1000LL /**< number of stalling nodes if earlyterm is true */
#define CONSHDLRFULLNORM          FALSE /**< compute real cut and compute norm for this (if addviolconshdlr and conshdlrusenorm are true) */
#define MINEFFICACY                0.05 /**< minimum efficacy of a cut - compare set.c */
#define MAXNSOLS                   1000 /**< maximal number of solutions stored in sub-SCIP */
#define OBJWEIGHTRANGE             0.01 /**< maximal range of scaling of objective w.r.t. size of rows */

/* parameters used for CMIR-generation (taken from sepa_gomory) */
#define BOUNDSWITCH              0.9999
#define USEVBDS                    TRUE
#define POSTPROCESS                TRUE
#define MINFRAC                  0.0009 /* to allow a deviation of the same size as EPSILONVALUE */
#define MAXFRAC                  0.9991 /* to allow a deviation of the same size as EPSILONVALUE */
#define FIXINTEGRALRHS            FALSE
#define MAKECONTINTEGRAL          FALSE
#define MAXWEIGHTRANGE            1e+05 /**< maximal valid range max(|weights|)/min(|weights|) of row weights */

#define MAXAGGRLEN(nvars)         nvars      /**< currently very large to allow any generation; an alternative would be (0.1*(nvars)+1000) */

/** separator data */
struct SCIP_SepaData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxdepth;           /**< maximal depth at which the separator is applied */
   SCIP_Bool             decisiontree;       /**< Use decision tree to turn separation on/off? */
   SCIP_Real             timelimit;          /**< time limit for subscip */
   SCIP_Real             memorylimit;        /**< memory limit for subscip */
   SCIP_Longint          minnodelimit;       /**< minimum number of nodes considered for sub-MIP (-1: unlimited) */
   SCIP_Longint          maxnodelimit;       /**< maximum number of nodes considered for sub-MIP (-1: unlimited) */
   SCIP_Real             cutcoefbnd;         /**< bounds on the values of the coefficients in the CG-cut */
   SCIP_Bool             onlyactiverows;     /**< Use only active rows to generate cuts? */
   int                   maxrowage;          /**< maximal age of rows to consider if onlyactiverows is false */
   SCIP_Bool             onlyrankone;        /**< Separate only rank 1 inequalities w.r.t. CG-MIP separator? */
   SCIP_Bool             onlyintvars;        /**< Generate cuts for problems with only integer variables? */
   SCIP_Bool             allowlocal;         /**< Allow local cuts? */
   SCIP_Bool             contconvert;        /**< Convert some integral variables to be continuous to reduce the size of the sub-MIP? */
   SCIP_Real             contconvfrac;       /**< fraction of integral variables converted to be continuous (if contconvert) */
   int                   contconvmin;        /**< minimum number of integral variables before some are converted to be continuous */
   SCIP_Bool             intconvert;         /**< Convert some integral variables attaining fractional values to have integral value? */
   SCIP_Real             intconvfrac;        /**< fraction of frac. integral variables converted to have integral value (if intconvert) */
   int                   intconvmin;         /**< minimum number of integral variables before some are converted to have integral value */
   SCIP_Bool             skipmultbounds;     /**< Skip the upper bounds on the multipliers in the sub-MIP? */
   SCIP_Bool             objlone;            /**< Should the objective of the sub-MIP only minimize the l1-norm of the multipliers? */
   SCIP_Real             objweight;          /**< objective weight for artificial variables */
   SCIP_Bool             objweightsize;      /**< Weight each row by its size? */
   SCIP_Bool             dynamiccuts;        /**< Should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Bool             usecmir;            /**< Use CMIR-generator (otherwise add cut directly)? */
   SCIP_Bool             usestrongcg;        /**< Use strong CG-function to strengthen cut? */
   SCIP_Bool             cmirownbounds;      /**< Tell CMIR-generator which bounds to used in rounding? */
   SCIP_Bool             usecutpool;         /**< Use cutpool to store CG-cuts even if the are not efficient? */
   SCIP_Bool             primalseparation;   /**< Only separate cuts that are tight for the best feasible solution? */
   SCIP_Bool             earlyterm;          /**< Terminate separation if a violated (but possibly sub-optimal) cut has been found? */
   SCIP_Bool             addviolationcons;   /**< Add constraint to subscip that only allows violated cuts? */
   SCIP_Bool             addviolconshdlr;    /**< Add constraint handler to filter out violated cuts? */
   SCIP_Bool             conshdlrusenorm;    /**< Should the violation constraint handler use the cut-norm to check for feasibility? */
   SCIP_Bool             useobjub;           /**< Use upper bound on objective function (via primal solution)? */
   SCIP_Bool             useobjlb;           /**< Use lower bound on objective function (via lower bound)? */
   SCIP_Bool             subscipfast;        /**< Should the settings for the sub-MIP be optimized for speed? */
   SCIP_Bool             output;             /**< Should information about the sub-MIP and cuts be displayed? */
};


/** what happens for columns in the LP */
enum CGMIP_ColType
{
   colPresent    = 0,    /**< column is present in the separating MIP */
   colContinuous = 1,    /**< column corresponds to a continuous variable */
   colConverted  = 2,    /**< column is converted to be continuous */
   colAtUb       = 3,    /**< variable corresponding to column was at it's upper bound and was complemented */
   colAtLb       = 4     /**< variable corresponding to column was at it's lower bound (possibly complemented) */
};
typedef enum CGMIP_ColType CGMIP_COLTYPE;


/** data for the sub-MIP */
struct CGMIP_MIPData
{
   SCIP*                 subscip;            /**< pointer to (sub)SCIP data structure containing the auxiliary IP */
   unsigned int          m;                  /**< number of constraints of subscip */
   unsigned int          n;                  /**< number of variables of subscip */
   unsigned int          nrows;              /**< number of rows of original LP */
   unsigned int          ncols;              /**< number of columns of original LP */
   unsigned int          ntotalrows;         /**< number of total rows used (possibly including objective rows) */

   SCIP_VAR**            alpha;              /**< cut coefficient variable (NULL if not in separating MIP) */
   SCIP_VAR*             beta;               /**< rhs of cut */
   SCIP_VAR**            fracalpha;          /**< fractional part of lhs of cut (NULL if not present) */
   SCIP_VAR*             fracbeta;           /**< fractional part of rhs of cut */
   CGMIP_COLTYPE*        coltype;            /**< type for the columns */
   SCIP_Bool*            iscomplemented;     /**< whether the variable was complemented */
   SCIP_Bool*            isshifted;          /**< whether the variable was shifted to have 0 lower bound */

   SCIP_VAR**            ylhs;               /**< auxiliary row variables for lhs (NULL if not present) */
   SCIP_VAR**            yrhs;               /**< auxiliary row variables for rhs (NULL if not present) */

   SCIP_VAR**            z;                  /**< auxiliary variables for upper bounds (NULL if not present) */

   char                  normtype;           /**< type of norm to use for efficacy norm calculation */

   /* additional redundant data */
   SCIP_Bool             conshdlrusenorm;    /**< copy from sepadata */
   SCIP_Bool             conshdlrfullnorm;   /**< compute real cut and compute norm for this (if addviolconshdlr and conshdlrusenorm are true) */
   SCIP*                 scip;               /**< original SCIP */
   SCIP_SEPA*            sepa;               /**< CG-cut separator */
   SCIP_SEPADATA*        sepadata;           /**< CG-cut separator data */
};
typedef struct CGMIP_MIPData CGMIP_MIPDATA;


/*
 * constraint handler to filter out violated cuts
 */

/* constraint handler properties */
#define CONSHDLR_NAME          "violatedCuts"
#define CONSHDLR_DESC          "only allow solutions corresponding to violated cuts"

/** constraint handler data */
struct SCIP_ConshdlrData
{
   CGMIP_MIPDATA*        mipdata;            /**< data of separating sub-MIP */
};

/* temporary forward declaration */
static
SCIP_RETCODE computeCut(
   SCIP*                 scip,               /**< original scip */
   SCIP_SEPA*            sepa,               /**< separator */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< current solution for sub-MIP */
   SCIP_Real*            cutcoefs,           /**< coefficients of the cut */
   SCIP_Real*            cutrhs,             /**< rhs of the cut */
   SCIP_Bool*            localrowsused,      /**< pointer to store whether local rows were used in summation */
   SCIP_Bool*            localboundsused,    /**< pointer to store whether local bounds were used in summation */
   int *                 cutrank,            /**< pointer to store the cut rank */
   SCIP_Bool*            success             /**< whether we produced a valid cut */
   );

/** check whether cut corresponding to solution is violated */
static
SCIP_RETCODE solCutIsViolated(
   SCIP*                 scip,               /**< SCIP data structure */
   CGMIP_MIPDATA*        mipdata,            /**< data of separating sub-MIP */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_Bool*            violated            /**< pointer to store if the cut is violated */
   )
{
   SCIP_Real cutsqrnorm = 0.0;
   SCIP* subscip;
   SCIP_Real act;
   SCIP_Real norm;
   SCIP_Real val;
   SCIP_VAR* var;
   SCIP_Real rhs;
   unsigned int j;
   int len = 0;

   assert( mipdata != NULL );
   subscip = mipdata->subscip;
   assert( subscip != NULL );
   assert( violated != NULL );

   /* initialize activity and norm */
   act = 0.0;
   norm = 1.0;
   *violated = FALSE;

   /* compute activity and norm  */
   if ( mipdata->conshdlrusenorm )
   {
      /* check whether we should compute the full cut and then compute the norm */
      if ( mipdata->conshdlrfullnorm )
      {
         SCIP_Real* cutcoefs;
         SCIP_Bool localrowsused;
         SCIP_Bool localboundsused;
         SCIP_Bool success;
         SCIP_VAR** vars;
         int cutrank = 0;
         int nvars;

         /* get data */
         SCIP_CALL( SCIPgetVarsData(mipdata->scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
         assert(nvars >= 0);
         SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );

         /* compute coefficients */
         SCIP_CALL( computeCut(mipdata->scip, mipdata->sepa, mipdata, mipdata->sepadata, sol, cutcoefs, &rhs, &localrowsused, &localboundsused, &cutrank, &success) );

#ifdef SCIP_MORE_DEBUG
         for (j = 0; j < (unsigned int) nvars; ++j)
         {
            if ( ! SCIPisZero(scip, cutcoefs[j]) )
               SCIPinfoMessage(scip, NULL, "+ %f x%d", cutcoefs[j], j);
         }
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         /* ignore solution if cut was not valid */
         if ( ! success )
            return SCIP_OKAY;

         /* compute activity and Euclidean norm (todo: use arbitrary norm) */
         cutsqrnorm = 0.0;
         for (j = 0; j < (unsigned int) nvars; ++j)
         {
            if ( ! SCIPisZero(scip, cutcoefs[j]) )
            {
               act += cutcoefs[j] * SCIPvarGetLPSol(vars[j]);
               cutsqrnorm += SQR(cutcoefs[j]);
            }
         }
         norm = SQRT(cutsqrnorm);

         SCIPfreeBufferArray(scip, &cutcoefs);
      }
      else
      {
         switch ( mipdata->normtype )
         {
         case 'e':
            cutsqrnorm = 0.0;
            for (j = 0; j < mipdata->ncols; ++j)
            {
               var = mipdata->alpha[j];
               if ( var == NULL )
                  continue;

               val = SCIPgetSolVal(subscip, sol, var);
               if ( !SCIPisZero(scip, val) )
               {
                  act += val * SCIPvarGetObj(var);
                  cutsqrnorm += SQR(val);
               }
            }
            norm = SQRT(cutsqrnorm);
            break;
         case 'm':
            for (j = 0; j < mipdata->ncols; ++j)
            {
               var = mipdata->alpha[j];
               if ( var == NULL )
                  continue;

               val = SCIPgetSolVal(subscip, sol, var);
               if ( !SCIPisZero(scip, val) )
               {
                  act += val * SCIPvarGetObj(var);
                  if ( REALABS(val) > norm )
                     norm = REALABS(val);
               }
            }
            break;
         case 's':
            for (j = 0; j < mipdata->ncols; ++j)
            {
               var = mipdata->alpha[j];
               if ( var == NULL )
                  continue;

               val = SCIPgetSolVal(subscip, sol, var);
               if ( !SCIPisZero(scip, val) )
               {
                  act += val * SCIPvarGetObj(var);
                  norm += REALABS(val);
               }
            }
            break;
         case 'd':
            for (j = 0; j < mipdata->ncols; ++j)
            {
               var = mipdata->alpha[j];
               if ( var == NULL )
                  continue;

               val = SCIPgetSolVal(subscip, sol, var);
               if ( !SCIPisZero(scip, val) )
               {
                  act += val * SCIPvarGetObj(var);
                  ++len;
               }
            }
            if ( len > 0 )
               norm = 1.0;
            break;
         default:
            SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", mipdata->normtype);
            return SCIP_INVALIDDATA;
         }
         /* get rhs */
         rhs = SCIPgetSolVal(subscip, sol, mipdata->beta);
      }

      /* if norm is 0, the cut is trivial */
      if ( SCIPisZero(subscip, norm) )
         return SCIP_OKAY;
   }
   else
   {
      for (j = 0; j < mipdata->ncols; ++j)
      {
         var = mipdata->alpha[j];
         if ( var == NULL )
            continue;

         val = SCIPgetSolVal(subscip, sol, var);
         if ( !SCIPisZero(subscip, val) )
            act += SCIPvarGetObj(var) * val;
      }

      /* get rhs */
      rhs = SCIPgetSolVal(subscip, sol, mipdata->beta);
   }

#ifdef SCIP_DEBUG
   if ( SCIPisEfficacious(subscip, (act - rhs)/norm) )
   {
      SCIPdebugMsg(scip, "Violated cut from solution - act: %f, rhs: %f, norm: %f, eff.: %f\n", act, rhs, norm, (act-rhs)/norm);
   }
   else
   {
      SCIPdebugMsg(scip, "Rejected cut from solution - act: %f, rhs: %f, norm: %f, eff.: %f\n", act, rhs, norm, (act-rhs)/norm);
   }
#endif

   *violated = SCIPisEfficacious(subscip, (act - rhs)/norm);

   return SCIP_OKAY;
}


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeViolatedCuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpViolatedCuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( result != NULL );

   assert( SCIPgetNLPBranchCands(scip) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIP_CALL( solCutIsViolated(scip, conshdlrdata->mipdata, NULL, &violated) );

   if ( violated )
      *result = SCIP_FEASIBLE;
   else
      *result = SCIP_CUTOFF;  /* cutoff, since all integer variables are integer, but the solution is not feasible */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsViolatedCuts)
{  /*lint --e{715}*/
   assert( result != NULL );

   /* this function should better not be called, since we need an LP solution for the sub-MIP to
    * make sense, because of the multiplier variables. We therefore return SCIP_FEASIBLE. */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckViolatedCuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( sol != NULL );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIP_CALL( solCutIsViolated(scip, conshdlrdata->mipdata, sol, &violated) );

   if ( violated )
      *result = SCIP_FEASIBLE;
   else
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockViolatedCuts)
{  /*lint --e{715}*/
   /* do not lock variables */
   return SCIP_OKAY;
}


/** creates the violated CG-cut constraint handler and includes it in SCIP */
static
SCIP_RETCODE SCIPincludeConshdlrViolatedCut(
   SCIP*                 scip,               /**< SCIP data structure */
   CGMIP_MIPDATA*        mipdata             /**< data of separating sub-MIP */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->mipdata = mipdata;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         -1000000, -1000000, 100, FALSE,
         consEnfolpViolatedCuts, consEnfopsViolatedCuts, consCheckViolatedCuts, consLockViolatedCuts,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeViolatedCuts) );

   return SCIP_OKAY;
}


/*
 * local methods
 */


/** stores nonzero elements of dense coefficient vector as sparse vector and calculates activity and norm
 *
 *  copied from sepa_gomory.c
 */
static
SCIP_RETCODE storeCutInArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of problem variables */
   SCIP_Real*            cutcoefs,           /**< dense coefficient vector */
   SCIP_Real*            varsolvals,         /**< dense variable LP solution vector */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   int*                  cutinds,            /**< array to store variables of sparse cut vector */
   SCIP_Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*                  cutlen,             /**< pointer to store number of nonzero entries in cut */
   SCIP_Real*            cutact,             /**< pointer to store activity of cut */
   SCIP_Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   SCIP_Real val;
   SCIP_Real cutsqrnorm;
   SCIP_Real act;
   SCIP_Real norm;
   int len;
   int v;

   assert( nvars == 0 || cutcoefs != NULL );
   assert( nvars == 0 || varsolvals != NULL );
   assert( cutinds != NULL );
   assert( cutvals != NULL );
   assert( cutlen != NULL );
   assert( cutact != NULL );
   assert( cutnorm != NULL );

   len = 0;
   act = 0.0;
   norm = 0.0;
   switch ( normtype )
   {
   case 'e':
      cutsqrnorm = 0.0;
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutinds[len] = v;
            cutvals[len++] = val;
         }
      }
      norm = SQRT(cutsqrnorm);
      break;
   case 'm':
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            if ( REALABS(val) > norm )
               norm = REALABS(val);
            cutinds[len] = v;
            cutvals[len++] = val;
         }
      }
      break;
   case 's':
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm += REALABS(val);
            cutinds[len] = v;
            cutvals[len++] = val;
         }
      }
      break;
   case 'd':
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutinds[len] = v;
            cutvals[len++] = val;
         }
      }
      if ( len > 0 )
         norm = 1.0;
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", normtype);
      return SCIP_INVALIDDATA;
   }

   *cutlen = len;
   *cutact = act;
   *cutnorm = norm;

   return SCIP_OKAY;
}


/** Compute lhs/rhs for transformed column
 *
 *  Consider a variable \f$x_j\f$ and some row of the original system:
 *  \f[
 *       \gamma \leq a^T x \leq \delta, \quad \ell_j \leq x_j \leq u_j.
 *  \f]
 *  We perform the transformation
 *  \f[
 *       x_i' = \left\{
 *       \begin{array}{ll}
 *         s + \frac{1}{\sigma}\, x_j & \mbox{if }i = j\\
 *         x_i              & \mbox{otherwise},
 *       \end{array}
 *       \right.
 *  \f]
 *  where \f$s\f$ is the offset value and \f$\sigma\f$ is a scaling factor. The new system is
 *  \f[
 *     \gamma + \sigma\, a_j\,s \leq \sum_{i \neq j} a_i\, x_i' + \sigma a_j\, x_j' \leq \delta + \sigma\, a_j\, s
 *  \f]
 *  with bounds
 *  \f[
 *     \frac{1}{\sigma} \ell_j + s \leq x_j' \leq \frac{1}{\sigma} u_j + s, \qquad \mbox{ if }\sigma > 0
 *  \f]
 *  and
 *  \f[
 *     \frac{1}{\sigma} u_j + s \leq x_j' \leq \frac{1}{\sigma} \ell_j + s, \qquad \mbox{ if }\sigma < 0.
 *  \f]
 *
 *  This can be used as follows:
 *
 *  - If \f$x_j \geq \ell_j\f$ has a (nonzero) lower bound, one can use \f$s = -\ell_j\f$, \f$\sigma = 1\f$,
 *    and obtain \f$\gamma - a_j\,\ell_j \leq a^T x' \leq \delta - a_j\,\ell_j\f$, \f$0 \leq x_j' \leq u_j - \ell_j\f$.
 *
 *  - If \f$x_j \leq u_j\f$ has a (nonzero) upper bound, one can use \f$s = u_j\f$, \f$\sigma = -1\f$,
 *    and obtain \f$\gamma - a_j\,u_j \leq \sum_{i \neq j} a_i\, x_i' - a_j\, x_j' \leq \delta - a_j\, u_j\f$,
 *    \f$0 \leq x_j' \leq u_j - \ell_j\f$.
 */
static
SCIP_RETCODE transformColumn(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_COL*             col,                /**< column that should be complemented */
   SCIP_Real             offset,             /**< offset by which column should be shifted */
   SCIP_Real             sigma,              /**< scaling factor */
   SCIP_Real*            lhs,                /**< array of lhs of rows */
   SCIP_Real*            rhs,                /**< array rhs of rows */
   SCIP_Real*            lb,                 /**< pointer to lb of column */
   SCIP_Real*            ub,                 /**< pointer to ub of column */
   SCIP_Real*            primsol             /**< pointer to solution value */
   )
{
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   int pos, i;

   assert( scip != NULL );
   assert( lhs != NULL );
   assert( rhs != NULL );
   assert( col != NULL );

   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   assert( SCIPcolGetNLPNonz(col) == 0 || colrows != NULL );
   assert( SCIPcolGetNLPNonz(col) == 0 || colvals != NULL );
   assert( ! SCIPisZero(scip, sigma) );

   /* loop through rows that contain column */
   for (i = 0; i < SCIPcolGetNLPNonz(col); ++i)
   {
      SCIP_ROW* row;

      row = colrows[i];
      assert( row != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
         continue;

      pos = SCIProwGetLPPos(row);
      assert( 0 <= pos && pos < (int) mipdata->nrows );

      assert( ! SCIPisInfinity(scip, lhs[pos]) );
      if ( ! SCIPisInfinity(scip, -lhs[pos]) )
         lhs[pos] += sigma * colvals[i] * offset;

      assert( ! SCIPisInfinity(scip, -rhs[pos]) );
      if ( ! SCIPisInfinity(scip, rhs[pos]) )
         rhs[pos] += sigma * colvals[i] * offset;
   }

   /* check objective function */
   if ( sepadata->useobjub || sepadata->useobjlb )
   {
      assert( SCIPisEQ(scip, SCIPcolGetObj(col), SCIPvarGetObj(SCIPcolGetVar(col))) );
      assert( mipdata->ntotalrows == mipdata->nrows + 1 );

      if ( ! SCIPisInfinity(scip, -lhs[mipdata->nrows]) )
         lhs[mipdata->nrows] += sigma * SCIPcolGetObj(col) * offset;

      if ( ! SCIPisInfinity(scip, rhs[mipdata->nrows]) )
         rhs[mipdata->nrows] += sigma * SCIPcolGetObj(col) * offset;
   }

   /* correct lower and upper bounds and solution */
   if ( SCIPisNegative(scip, sigma) )
   {
      SCIP_Real l;

      assert( ! SCIPisInfinity(scip, -*ub) );
      if ( ! SCIPisInfinity(scip, *ub) )
         l = *ub/sigma + offset;
      else
         l = -SCIPinfinity(scip);

      assert( ! SCIPisInfinity(scip, *lb) );
      if ( ! SCIPisInfinity(scip, -*lb) )
         *ub = *lb/sigma + offset;
      else
         *ub = SCIPinfinity(scip);
      *lb = l;
   }
   else
   {
      assert( ! SCIPisInfinity(scip, *lb) );
      if ( ! SCIPisInfinity(scip, -*lb) )
         *lb = *lb/sigma + offset;
      assert( ! SCIPisInfinity(scip, -*ub) );
      if ( ! SCIPisInfinity(scip, *ub) )
         *ub = *ub/sigma + offset;
   }
   *primsol = *primsol/sigma + offset;

   return SCIP_OKAY;
}


/** compute objective coefficient for rows that are weighted by size
 *
 *  The objective is computed by multiplying a default value by
 *  \f[
 *  1 - (r_{\mbox{max}} - r) \frac{1 - a}{r_{\mbox{max}} - r_{\mbox{min}}},
 *  \f]
 *  where \f$r\f$ is the size of the current row, \f$a \in [0,1]\f$ is a parameter, and \f$r_{\mbox{max}}\f$ and
 *  \f$r_{\mbox{min}}\f$ are the maximal and minimal size of a row, respectively.
 *
 *  Thus, if \f$r = r_{\mbox{max}}\f$, we get 1 and if \f$r = r_{\mbox{min}}\f$, we get \f$a\f$.
 */
static
SCIP_Real computeObjWeightSize(
   int                   rowsize,            /**< size of current row */
   int                   minrowsize,         /**< maximal size of rows */
   int                   maxrowsize          /**< minimal size of rows */
   )
{
   SCIP_Real a;

   assert( maxrowsize > 0 );
   assert( minrowsize < INT_MAX );
   assert( minrowsize <= maxrowsize );
   assert( minrowsize <= rowsize && rowsize <= maxrowsize );

   if ( minrowsize == maxrowsize )
      return 1.0;

   a = (1.0 - OBJWEIGHTRANGE)/((SCIP_Real) (maxrowsize - minrowsize));

   return 1.0 - a * ((SCIP_Real) (maxrowsize - rowsize));
}


/** Creates a subscip representing the separating MIP.
 *
 *  Let the constraints of the original MIP be of the following form:
 *  \f[
 *    \begin{array}{l@{\;}ll}
 *      a \leq A x + & C r & \leq b\\
 *      \ell \leq x & & \leq u\\
 *      c \leq & r & \leq d\\
 *      x \in Z^n.
 *    \end{array}
 *  \f]
 *  Here, some of the bounds may have value \f$\infty\f$ or \f$-\infty\f$.  Written in
 *  \f$\leq\f$-form this becomes:
 *  \f[
 *    \begin{array}{r@{\;}l}
 *      \tilde{A} x + \tilde{C} r & \leq \tilde{b}\\
 *      -x & \leq -\ell\\
 *      x & \leq u\\
 *      -r & \leq -c\\
 *      r & \leq d\\
 *      x \in Z^n,
 *    \end{array}
 *  \f]
 *  where we use
 *  \f[
 *    \tilde{A} =
 *    \left[
 *    \begin{array}{r}
 *      -A \\
 *      A
 *    \end{array}
 *    \right],
 *    \quad
 *    \tilde{C} =
 *    \left[
 *    \begin{array}{r}
 *      - C\\
 *      C
 *    \end{array}
 *    \right]
 *    \qquad\mbox{ and }\qquad
 *    \tilde{b} =
 *    \left[
 *    \begin{array}{r}
 *      -a\\
 *      b
 *    \end{array}
 *    \right].
 *  \f]
 *  For the moment we assume that \f$c = 0\f$, i.e., the lower bounds on the continuous variables
 *  are 0.  To obtain a Chv&aacute;tal-Gomory cut we have to find nonnegative multipliers \f$y\f$,
 *  \f$\underline{z}\f$, and \f$\overline{z}\f$ such that
 *  \f[
 *      y^T \tilde{A} - \underline{z}^T + \overline{z}^T  \in Z \qquad\mbox{ and }\qquad
 *      y^T \tilde{C} \geq 0.
 *  \f]
 *  Note that we use zero multipliers for the bounds on the continuous variables \f$r\f$. Moreover,
 *  if some bounds are infinity, the corresponding multipliers are assumed to be 0. From these
 *  conditions, we obtain
 *  \f[
 *      (y^T \tilde{A} - \underline{z}^T + \overline{z}^T)\, x +
 *      y^T \tilde{C} \, r \leq
 *      y^T \tilde{b} - \underline{z}^T \ell + \overline{z}^T u.
 *  \f]
 *  Because \f$r \geq 0\f$, we can ignore the term \f$y^T \tilde{C} \, r \geq 0\f$ and obtain the
 *  following cut:
 *  \f[
 *      (y^T \tilde{A} - \underline{z}^T + \overline{z}^T )\, x \leq
 *      \lfloor y^T \tilde{b} - \underline{z}^T \ell + \overline{z}^T u \rfloor.
 *  \f]
 *  Assume that \f$\ell = 0\f$ for the meantime. Then the cut can be written as:
 *  \f[
 *      \lfloor y^T \tilde{A} + \overline{z}^T \rfloor \, x \leq
 *      \lfloor y^T \tilde{b} + \overline{z}^T u \rfloor.
 *  \f]
 *
 *  Following Fischetti and Lodi [2005], let \f$(x^*,r^*)\f$ be a fractional solution of the above
 *  original system.  The separating MIP created below is
 *  \f[
 *    \begin{array}{rlr@{\;}l}
 *       \max & \multicolumn{2}{@{}l}{(x^*)^T \alpha - \beta - w^T y} &\\
 *            & f = & \tilde{A}^T y + \overline{z} - \alpha & \\
 *            & \tilde{f} = & \tilde{b}^T y + u^T \overline{z} - \beta &\\
 *            & & \tilde{C}^T y & \geq 0\\
 *            & & 0 \leq f & \leq 1 - \epsilon \\
 *            & & 0 \leq \tilde{f} & \leq 1 - \epsilon\\
 *            & & 0 \leq y, \overline{z} & \leq 1 - \epsilon.\\
 *            & & \alpha \in Z^m, \beta & \in Z.
 *    \end{array}
 *  \f]
 *  Here, \f$w\f$ is a weight vector; it's idea is to make the sum over all components of \f$y\f$ as
 *  small as possible, in order to generate sparse cuts.
 *
 *  We perform the following additional computations:
 *
 *  - If the lower bounds on \f$x_i\f$ or \f$r_j\f$ are finite, we shift the variable to have a zero
 *    lower bound, i.e., we replace it by \f$x_i - \ell_i\f$ (or \f$r_j - u_j\f$). This is helpful in
 *    several ways: As seen above, the resulting inequalities/formulations simplify. Moreover, it
 *    allows to drop a variable if \f$x^*_i = 0\f$, see the next comment. If the lower bounds are not
 *    finite, but the upper bounds are finite, we can complement the variable. If the variables are
 *    free, the above formulation changes as follows: For free continuous variables, we require
 *    \f$\tilde{C}^T y = 0\f$. For a free integer variable \f$x_j\f$ (which rarely occurs in
 *    practice), we require \f$f_j = 0\f$, i.e., we force that \f$(\tilde{A}^T y + \overline{z})_j =
 *    \alpha_j\f$.
 *
 *  - If \f$x^*_j = 0 = \ell_j\f$ (after the above preprocessing), we drop variable \f$\alpha_j\f$
 *    from the formulation. Let \f$(\alpha^*, \beta^*, y^*, \overline{z}^*)\f$ be an
 *    optimal solution to the separating MIP. Then we can compute \f$\alpha_j =
 *    \lfloor(\tilde{A}_j^T y^* + \overline{z}^*)\rfloor\f$.
 *
 *  - If \f$x^*_i = u_i\f$, we complement the variable and drop it from the formulation, since the
 *    lower bound is 0 afterwards.
 *
 *  - If a variable has been shifted or complemented, we have to recompute \f$\beta\f$ with the
 *    original lhs/rhs.
 *
 *  - If a continuous variable \f$r_j\f$ is free, we have to force equality for the corresponding components in
 *    \f$y^T \tilde{C} \, r \geq 0\f$.
 *
 *  - If an integer variable \f$x_i\f$ is free, we are not allowed to round the cut down. In this
 *    case, the combintation of rows and bounds has to be integral. We force this by requiring that
 *    \f$f_i = 0\f$.
 *
 *  - If @p contconvert is true some integral variables are randomly treated as if they were
 *    continuous. This has the effect that in the resulting cut the corresponding coefficient has
 *    value 0. This makes the cuts more sparse. Moreover, the separation problems should become
 *    easier.
 *
 *  - If required, i.e., parameter @p primalseparation is true, we force a primal separation step. For
 *    this we require that the cut is tight at the currently best solution. To get reliable solutions
 *    we relax equality by EPSILONVALUE.
 *
 *  - If required (via parameters @p useobjub or @p useobjlb), we add a row corresponding to the objective function with
 *    respect to the current lower and upper bounds.
 */
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata             /**< data for sub-MIP */
   )
{
   SCIP* subscip;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* primsol;
   SCIP_Real multvarub;

   unsigned int cnt;
   unsigned int ucnt;
   unsigned int nshifted;
   unsigned int ncomplemented;
   unsigned int ncontconverted;
   unsigned int nintconverted;
   unsigned int nlbounds;
   unsigned int nubounds;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_CONS* cons;
   int nconsvars;
   char name[SCIP_MAXSTRLEN];

   int ncols;
   int nrows;
   int ntotalrows;
   int maxrowsize = 0;
   int minrowsize = INT_MAX;
   int i, j;

   assert( scip != NULL );
   assert( sepadata != NULL );

   assert( mipdata->subscip == NULL );

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert( ncols > 0 && nrows > 0 );

   mipdata->m = 0;
   mipdata->n = 0;
   mipdata->nrows = (unsigned int) nrows;
   mipdata->ncols = (unsigned int) ncols;
   mipdata->ntotalrows = mipdata->nrows;

   if ( sepadata->useobjub || sepadata->useobjlb )
      mipdata->ntotalrows = mipdata->nrows + 1;

   assert(mipdata->ntotalrows <= INT_MAX);
   ntotalrows = (int) mipdata->ntotalrows;

   /* copy value */
   mipdata->conshdlrusenorm = sepadata->conshdlrusenorm;

   /* create subscip */
   SCIP_CALL( SCIPcreate( &(mipdata->subscip) ) );
   subscip = mipdata->subscip;
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* add violation constraint handler if requested */
   if ( sepadata->addviolconshdlr )
   {
      SCIP_CALL( SCIPincludeConshdlrViolatedCut(subscip, mipdata) );
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sepa_cgmip separating MIP (%s)", SCIPgetProbName(scip));
   SCIP_CALL( SCIPcreateProb(subscip, name, NULL, NULL , NULL , NULL , NULL , NULL , NULL) );
   SCIP_CALL( SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE) );

   /* alloc memory for subscipdata elements */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->alpha), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->fracalpha), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->coltype), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->iscomplemented), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->isshifted), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->ylhs), ntotalrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->yrhs), ntotalrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->z), 2*ncols) );

   /* get temporary storage */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, ntotalrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, ntotalrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &primsol, ncols) );

   /* store lhs/rhs for complementing (see below) and compute maximal nonzeros of candidate rows */
   for (i = 0; i < nrows; ++i)
   {
      SCIP_Real val;
      SCIP_ROW* row;

      row = rows[i];
      assert( row != NULL );

      val = SCIProwGetLhs(row) - SCIProwGetConstant(row);
      if ( SCIProwIsIntegral(row) )
         val = SCIPfeasCeil(scip, val); /* row is integral: round left hand side up */
      lhs[i] = val;

      val = SCIProwGetRhs(row) - SCIProwGetConstant(row);
      if ( SCIProwIsIntegral(row) )
         val = SCIPfeasFloor(scip, val); /* row is integral: round right hand side down */
      rhs[i] = val;

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
         continue;

      /* skip rows that not have been active for a longer time */
      if ( ! sepadata->onlyactiverows && sepadata->maxrowage > 0 && SCIProwGetAge(row) > sepadata->maxrowage )
         continue;

      /* check whether we want to skip cut produced by the CGMIP separator */
      if ( sepadata->onlyrankone )
      {
         if ( SCIProwGetOriginSepa(row) == sepa )
            continue;
      }

      /* determine maximal row size: */
      val = SCIPgetRowLPActivity(scip, row);
      if ( ! SCIPisInfinity(scip, REALABS(lhs[i])) )
      {
         if ( ! sepadata->onlyactiverows || SCIPisFeasEQ(scip, val, SCIProwGetLhs(row)) )
         {
            if ( SCIProwGetNLPNonz(row) > maxrowsize )
               maxrowsize = SCIProwGetNLPNonz(row);
            if ( SCIProwGetNLPNonz(row) < minrowsize )
               minrowsize = SCIProwGetNLPNonz(row);
         }
      }
      else
      {
         if ( ! SCIPisInfinity(scip, rhs[i]) )
         {
            if ( ! sepadata->onlyactiverows || SCIPisFeasEQ(scip, val, SCIProwGetRhs(row)) )
            {
               if ( SCIProwGetNLPNonz(row) > maxrowsize )
                  maxrowsize = SCIProwGetNLPNonz(row);
               if ( SCIProwGetNLPNonz(row) < minrowsize )
                  minrowsize = SCIProwGetNLPNonz(row);
            }
         }
      }
   }
   assert( maxrowsize > 0 );
   assert( minrowsize < INT_MAX );

   /* add cuts for objective function if required */
   if ( sepadata->useobjub )
   {
      assert( mipdata->ntotalrows == mipdata->nrows + 1 );
      rhs[mipdata->nrows] = SCIPgetUpperbound(scip);
      assert( ! SCIPisObjIntegral(scip) || SCIPisFeasIntegral(scip, SCIPgetUpperbound(scip)) );

      if ( ! SCIPisInfinity(scip, SCIPgetUpperbound(scip)) && SCIPgetNObjVars(scip) > maxrowsize )
         maxrowsize = SCIPgetNObjVars(scip);
      if ( ! SCIPisInfinity(scip, SCIPgetUpperbound(scip)) && SCIPgetNObjVars(scip) < minrowsize )
         minrowsize = SCIPgetNObjVars(scip);
   }
   if ( sepadata->useobjlb )
   {
      assert( mipdata->ntotalrows == mipdata->nrows + 1 );

      if ( SCIPisObjIntegral(scip) )
         lhs[mipdata->nrows] = SCIPfeasCeil(scip, SCIPgetLowerbound(scip));
      else
         lhs[mipdata->nrows] = SCIPgetLowerbound(scip);

      if ( ! SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) && SCIPgetNObjVars(scip) > maxrowsize )
         maxrowsize = SCIPgetNObjVars(scip);
      if ( ! SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) && SCIPgetNObjVars(scip) < minrowsize )
         minrowsize = SCIPgetNObjVars(scip);
   }

   /* store lb/ub for complementing and perform preprocessing */
   nshifted = 0;
   ncomplemented = 0;
   ncontconverted = 0;
   nintconverted = 0;
   nlbounds = 0;
   nubounds = 0;
   for (j = 0; j < ncols; ++j)
   {
      SCIP_COL* col;
      SCIP_VAR* var;

      col = cols[j];
      assert( col != NULL );
      var = SCIPcolGetVar(col);
      assert( var != NULL );

      primsol[j] = SCIPcolGetPrimsol(col);
      assert( SCIPisEQ(scip, SCIPgetVarSol(scip, var), primsol[j]) );

      lb[j] = SCIPvarGetLbGlobal(var);
      assert( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPcolGetLb(col)) );

      /* if allowed, try to use stronger local bound */
      if ( sepadata->allowlocal && SCIPisGT(scip, SCIPvarGetLbLocal(var), lb[j]) )
         lb[j] = SCIPvarGetLbLocal(var);

      ub[j] = SCIPvarGetUbGlobal(var);
      assert( SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPcolGetUb(col)) );

      /* if allowed, try to use stronger local bound */
      if ( sepadata->allowlocal && SCIPisLT(scip, SCIPvarGetUbLocal(var), ub[j]) )
         ub[j] = SCIPvarGetUbLocal(var);

      mipdata->coltype[j] = colPresent;
      mipdata->iscomplemented[j] = FALSE;
      mipdata->isshifted[j] = FALSE;

      /* check status of column/variable */
      if ( SCIPcolIsIntegral(col) )
      {
         /* integral variables taking integral values are not interesting - will be substituted out below */
         if ( ! SCIPisFeasIntegral(scip, primsol[j]) )
         {
            /* possibly convert fractional integral variables to take integral values */
            if ( sepadata->intconvert && ncols >= sepadata->intconvmin )
            {
               /* randomly convert variables */
               if ( SCIPrandomGetReal(sepadata->randnumgen, 0.0, 1.0) <= sepadata->intconvfrac )
               {
                  assert( ! SCIPisInfinity(scip, ub[j]) || ! SCIPisInfinity(scip, -lb[j]) );

                  /* if both bounds are finite, take the closer one */
                  if ( ! SCIPisInfinity(scip, ub[j]) && ! SCIPisInfinity(scip, -lb[j]) )
                  {
                     assert( SCIPisFeasIntegral(scip, ub[j]) );
                     assert( SCIPisFeasIntegral(scip, lb[j]) );
                     assert( SCIPisFeasLT(scip, primsol[j], ub[j]) );
                     assert( SCIPisFeasGT(scip, primsol[j], lb[j]) );
                     if ( ub[j] - primsol[j] < primsol[j] - lb[j] )
                        primsol[j] = ub[j];
                     else
                        primsol[j] = lb[j];
                     ++nintconverted;
                  }
                  else
                  {
                     /* if only lower bound is finite */
                     if ( ! SCIPisInfinity(scip, -lb[j]) )
                     {
                        assert( SCIPisFeasIntegral(scip, lb[j]) );
                        primsol[j] = lb[j];
                        ++nintconverted;
                     }
                     else
                     {
                        assert( ! SCIPisInfinity(scip, ub[j]) );
                        assert( SCIPisFeasIntegral(scip, ub[j]) );
                        primsol[j] = ub[j];
                        ++nintconverted;
                     }
                  }
               }
            }
         }

         /* integral variables taking integral values are not interesting - will be substituted out below */
         if ( ! SCIPisFeasIntegral(scip, primsol[j]) )
         {
            /* possibly convert integral variables to be continuous */
            if ( sepadata->contconvert && ncols >= sepadata->contconvmin )
            {
               /* randomly convert variables */
               if ( SCIPrandomGetReal(sepadata->randnumgen, 0.0, 1.0) <= sepadata->contconvfrac )
               {
                  /* preprocessing is also performed for converted columns */
                  mipdata->coltype[j] = colConverted;
                  ++ncontconverted;
               }
            }
         }
      }
      else
      {
         /* detect continuous variables, but perform preprocessing for them */
         mipdata->coltype[j] = colContinuous;
      }

      /* if integer variable is at its upper bound -> complementing (this also generates a 0 lower bound) */
      if ( mipdata->coltype[j] == colPresent && SCIPisFeasEQ(scip, primsol[j], ub[j]) )
      {
         assert( ! SCIPisInfinity(scip, ub[j]) );
         SCIP_CALL( transformColumn(scip, sepadata, mipdata, col, ub[j], -1.0, lhs, rhs, &(lb[j]), &(ub[j]), &(primsol[j])) );
         mipdata->iscomplemented[j] = TRUE;
         mipdata->coltype[j] = colAtUb;
         ++nubounds;
      }
      else
      {
         /* if a variable has a finite nonzero lower bound -> shift */
         if ( ! SCIPisInfinity(scip, -lb[j]) )
         {
            if ( ! SCIPisZero(scip, lb[j]) )
            {
               SCIP_CALL( transformColumn(scip, sepadata, mipdata, col, -lb[j], 1.0, lhs, rhs, &(lb[j]), &(ub[j]), &(primsol[j])) );
               assert( SCIPisZero(scip, lb[j]) );
               mipdata->isshifted[j] = TRUE;
               ++nshifted;
            }

            /* if integer variable is at its lower bound */
            if ( mipdata->coltype[j] == colPresent && SCIPisZero(scip, primsol[j]) )
            {
               mipdata->coltype[j] = colAtLb;
               ++nlbounds;
            }
         }
         else
         {
            /* lower bound is minus-infinity -> check whether upper bound is finite */
            if ( ! SCIPisInfinity(scip, ub[j]) )
            {
               /* complement variable */
               SCIP_CALL( transformColumn(scip, sepadata, mipdata, col, ub[j], -1.0, lhs, rhs, &(lb[j]), &(ub[j]), &(primsol[j])) );
               assert( SCIPisZero(scip, lb[j]) );
               mipdata->iscomplemented[j] = TRUE;
               ++ncomplemented;

               /* if integer variable is at its lower bound */
               if ( mipdata->coltype[j] == colPresent && SCIPisZero(scip, primsol[j]) )
               {
                  mipdata->coltype[j] = colAtLb;
                  ++nlbounds;
               }
            }
         }
      }

      assert( SCIPisFeasLE(scip, lb[j], primsol[j]) );
      assert( SCIPisFeasLE(scip, primsol[j], ub[j]) );
   }

#ifndef NDEBUG
   if ( sepadata->intconvert && ncols >= sepadata->intconvmin )
   {
      SCIPdebugMsg(scip, "Converted %u fractional integral variables to have integral value.\n", nintconverted);
   }
   if ( sepadata->contconvert && ncols >= sepadata->contconvmin )
   {
      SCIPdebugMsg(scip, "Converted %u integral variables to be continuous.\n", ncontconverted);
   }
#endif
   SCIPdebugMsg(scip, "original variables: %d integral, %d continuous, %u shifted, %u complemented, %u at lb, %u at ub\n",
      SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip), SCIPgetNContVars(scip),
      nshifted, ncomplemented, nlbounds, nubounds);

   /* prepare upper bound on y-variables */
   if ( sepadata->skipmultbounds )
      multvarub = SCIPinfinity(scip);
   else
      multvarub =  1.0-EPSILONVALUE;

   /* create artificial variables for row combinations (y-variables) */
   cnt = 0;
   for (i = 0; i < nrows; ++i)
   {
      SCIP_ROW* row;

      row = rows[i];
      assert( row != NULL );

      mipdata->ylhs[i] = NULL;
      mipdata->yrhs[i] = NULL;

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
         continue;

      /* skip rows that not have been active for a longer time */
      if ( ! sepadata->onlyactiverows && sepadata->maxrowage > 0 && SCIProwGetAge(row) > sepadata->maxrowage )
         continue;

      /* check whether we want to skip cut produced by the CGMIP separator */
      if ( sepadata->onlyrankone )
      {
         if ( SCIProwGetOriginSepa(row) == sepa )
            continue;
      }

      /* if we have an equation */
      if ( SCIPisEQ(scip, lhs[i], rhs[i]) )
      {
         SCIP_Real weight = -sepadata->objweight;

         assert( ! SCIPisInfinity(scip, rhs[i]) );
         assert( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row)) ); /* equations should always be active */
         assert( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetRhs(row)) );

         if ( sepadata->objweightsize )
            weight = - sepadata->objweight * computeObjWeightSize(SCIProwGetNLPNonz(row), minrowsize, maxrowsize);

         /* create two variables for each equation */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yeq1_%d", i);
         SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->ylhs[i]), name, 0.0, multvarub,
               weight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(subscip, mipdata->ylhs[i]) );
         ++cnt;

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "Created variable <%s> for equation <%s>.\n", name, SCIProwGetName(row));
#endif

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yeq2_%d", i);
         SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->yrhs[i]), name, 0.0, multvarub,
               weight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(subscip, mipdata->yrhs[i]) );
         ++cnt;

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "Created variable <%s> for equation <%s>.\n", name, SCIProwGetName(row));
#endif
      }
      else
      {
         /* create variable for lhs of row if necessary */
         if ( ! SCIPisInfinity(scip, -lhs[i]) )
         {
            SCIP_Bool isactive = FALSE;
            SCIP_Real weight = 0.0;

            /* if the row is active, use objective weight equal to -sepadata->objweight */
            if ( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row)) )
            {
               isactive = TRUE;
               if ( sepadata->objweightsize )
                  weight = -sepadata->objweight * computeObjWeightSize(SCIProwGetNLPNonz(row), minrowsize, maxrowsize);
               else
                  weight = -sepadata->objweight;
            }

            if ( ! sepadata->onlyactiverows || isactive )
            {
               /* add variable */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ylhs_%d", i);
               SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->ylhs[i]), name, 0.0, multvarub,
                     weight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(subscip, mipdata->ylhs[i]) );
               ++cnt;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "Created variable <%s> for >= inequality <%s> (weight: %f).\n", name, SCIProwGetName(row), weight);
#endif
            }
         }

         /* create variable for rhs of row if necessary */
         if ( ! SCIPisInfinity(scip, rhs[i]) )
         {
            SCIP_Bool isactive = FALSE;
            SCIP_Real weight = 0.0;

            /* if the row is active, use objective weight equal to -sepadata->objweight */
            if ( SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetRhs(row)) )
            {
               isactive = TRUE;
               if ( sepadata->objweightsize )
                  weight = -sepadata->objweight * computeObjWeightSize(SCIProwGetNLPNonz(row), minrowsize, maxrowsize);
               else
                  weight = -sepadata->objweight;
            }

            if ( ! sepadata->onlyactiverows || isactive )
            {
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yrhs_%d", i);
               SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->yrhs[i]), name, 0.0, multvarub,
                     weight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(subscip, mipdata->yrhs[i]) );
               ++cnt;

#ifdef SCIP_MORE_DEBUG
               SCIPdebugMsg(scip, "Created variable <%s> for <= inequality <%s> (weight: %f).\n", name, SCIProwGetName(row), weight);
#endif
            }
         }
      }
   }
   assert( (int) cnt <= 2 * nrows );
   mipdata->n += cnt;

   /* create artificial variables for objective function (if required) (y-variables) */
   if ( sepadata->useobjub || sepadata->useobjlb )
   {
      SCIP_Real weight = 0.0;

      assert( mipdata->ntotalrows == mipdata->nrows + 1 );
      mipdata->ylhs[mipdata->nrows] = NULL;
      mipdata->yrhs[mipdata->nrows] = NULL;
      cnt = 0;

      if ( sepadata->objweightsize )
         weight = -sepadata->objweight * computeObjWeightSize(SCIPgetNObjVars(scip), minrowsize, maxrowsize);
      else
         weight = -sepadata->objweight;

      /* create variable for upper objective bound if necessary */
      if ( sepadata->useobjub && ! SCIPisInfinity(scip, rhs[mipdata->nrows]) )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yobjub");
         SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->yrhs[mipdata->nrows]), name, 0.0, multvarub,
               weight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(subscip, mipdata->yrhs[mipdata->nrows]) );
         ++cnt;

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "Created variable <%s> for upper bound on objective (weight: %f).\n", name, weight);
#endif
      }

      /* create variable for lower bound objective if necessary */
      if ( sepadata->useobjlb && ! SCIPisInfinity(scip, -lhs[mipdata->nrows]) )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yobjlb");
         SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->ylhs[mipdata->nrows]), name, 0.0, multvarub,
               weight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(subscip, mipdata->ylhs[mipdata->nrows]) );
         ++cnt;

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMsg(scip, "Created variable <%s> for lower bound on objective (weight: %f).\n", name, weight);
#endif
      }

      assert( (int) cnt <= 2 * ntotalrows );
      mipdata->n += cnt;
   }

   /* create alpha, bound, and fractional variables */
   cnt = 0;
   ucnt = 0;
   for (j = 0; j < ncols; ++j)
   {
      mipdata->z[j] = NULL;
      mipdata->alpha[j] = NULL;
      mipdata->fracalpha[j] = NULL;

      if ( mipdata->coltype[j] == colPresent )
      {
         SCIP_Real obj;

         if ( sepadata->objlone )
            obj = 0.0;
         else
            obj = primsol[j];

         /* create alpha variables */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "alpha_%d", j);
         SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->alpha[j]), name, -sepadata->cutcoefbnd, sepadata->cutcoefbnd, obj,
               SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(subscip, mipdata->alpha[j]) );
         ++cnt;

         /* create fractional variables */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "f_%d", j);
         if ( SCIPisInfinity(scip, -lb[j]) && SCIPisInfinity(scip, ub[j]) )
         {
            /* fix fractional value to be zero for free original variables */
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->fracalpha[j]), name, 0.0, 0.0, 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         }
         else
         {
            /* fractional value in [0, 1) for variables with finite bounds */
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->fracalpha[j]), name, 0.0, 1.0-EPSILONVALUE, 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         }
         SCIP_CALL( SCIPaddVar(subscip, mipdata->fracalpha[j]) );
         ++cnt;

         /* create variables for upper bounds */
         if ( ! SCIPisInfinity(scip, ub[j]) )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "zub_%d", j);
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->z[j]), name, 0.0, multvarub,
                  0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, mipdata->z[j]) );
            ++ucnt;
         }
      }
   }
   assert( (int) cnt <= 2 * ncols );
   assert( (int) ucnt <= ncols );

   /* create variable for the rhs of the cut */
   if ( sepadata->objlone )
   {
      SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->beta), "beta", -sepadata->cutcoefbnd, sepadata->cutcoefbnd, 0.0,
            SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   }
   else
   {
      SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->beta), "beta", -sepadata->cutcoefbnd, sepadata->cutcoefbnd, -1.0,
            SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   }
   SCIP_CALL( SCIPaddVar(subscip, mipdata->beta) );

   /* create fractional variable for the rhs */
   SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->fracbeta), "fracbeta", 0.0, 1.0-BETAEPSILONVALUE, 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(subscip, mipdata->fracbeta) );
   mipdata->n += cnt + ucnt + 2;

   /* get temporary storage */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, (int) mipdata->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, (int) mipdata->n) );

   /* create constraints for alpha variables of CG-cut */
   cnt = 0;
   for (j = 0; j < ncols; ++j)
   {
      SCIP_ROW** colrows;
      SCIP_Real* colvals;

      /* create ordinary part for all selected variables */
      if ( mipdata->coltype[j] == colPresent )
      {
         SCIP_Real sigma;

         assert( cols[j] != NULL );
         colrows = SCIPcolGetRows(cols[j]);
         colvals = SCIPcolGetVals(cols[j]);
         nconsvars = 0;

         if ( mipdata->iscomplemented[j] )
            sigma = -1.0;
         else
            sigma = 1.0;

         /* add part for columns */
         for (i = 0; i < SCIPcolGetNLPNonz(cols[j]); ++i)
         {
            SCIP_ROW* row;
            int pos;

            row = colrows[i];
            assert( row != NULL );

            /* skip modifiable rows and local rows, unless allowed */
            if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
               continue;

            pos = SCIProwGetLPPos(row);
            assert( 0 <= pos && pos < nrows );

            if ( mipdata->ylhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->ylhs[pos];
               consvals[nconsvars] = -sigma * colvals[i];
               ++nconsvars;
            }
            if ( mipdata->yrhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->yrhs[pos];
               consvals[nconsvars] = sigma * colvals[i];
               ++nconsvars;
            }
            assert( nconsvars <= (int) mipdata->n );
         }
         /* add part for upper bounds */
         if ( mipdata->z[j] != NULL )
         {
            assert( ! SCIPisInfinity(scip, ub[j]) );
            consvars[nconsvars] = mipdata->z[j];
            consvals[nconsvars] = 1.0;
            ++nconsvars;
         }
         assert( nconsvars <= (int) mipdata->n );

         /* add alpha variable */
         consvars[nconsvars] = mipdata->alpha[j];
         consvals[nconsvars] = -1.0;
         ++nconsvars;
         assert( nconsvars <= (int) mipdata->n );

         /* add fractional-alpha variable */
         consvars[nconsvars] = mipdata->fracalpha[j];
         consvals[nconsvars] = -1.0;
         ++nconsvars;
         assert( nconsvars <= (int) mipdata->n );

         /* check for lower and upper objective bounds */
         if ( (sepadata->useobjub || sepadata->useobjlb) && ! SCIPisZero(scip, SCIPcolGetObj(cols[j])) )
         {
            /* add lower objective bound */
            if ( mipdata->ylhs[mipdata->nrows] != NULL )
            {
               assert( sepadata->useobjlb );
               consvars[nconsvars] = mipdata->ylhs[mipdata->nrows];
               consvals[nconsvars] = -sigma * SCIPcolGetObj(cols[j]);
               ++nconsvars;
            }

            /* add upper objective bound */
            if ( mipdata->yrhs[mipdata->nrows] != NULL )
            {
               assert( sepadata->useobjub );
               consvars[nconsvars] = mipdata->yrhs[mipdata->nrows];
               consvals[nconsvars] = sigma * SCIPcolGetObj(cols[j]);
               ++nconsvars;
            }
         }

         /* add linear constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "alpha_%d", j);
         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, nconsvars, consvars, consvals, 0.0, 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         ++cnt;
      }
      /* generate part that makes sure that cut is valid for continuous variables */
      else if ( mipdata->coltype[j] == colContinuous || mipdata->coltype[j] == colConverted )
      {
         SCIP_Real sigma;
         SCIP_Real r;

         assert( cols[j] != NULL );
         colrows = SCIPcolGetRows(cols[j]);
         colvals = SCIPcolGetVals(cols[j]);
         nconsvars = 0;

         if ( mipdata->iscomplemented[j] )
            sigma = -1.0;
         else
            sigma = 1.0;

         /* add part for columns */
         for (i = 0; i < SCIPcolGetNLPNonz(cols[j]); ++i)
         {
            SCIP_ROW* row;
            int pos;

            row = colrows[i];
            assert( row != NULL );

            /* skip modifiable rows and local rows, unless allowed */
            if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
               continue;

            pos = SCIProwGetLPPos(row);
            assert( 0 <= pos && pos < nrows );

            if ( mipdata->ylhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->ylhs[pos];
               consvals[nconsvars] = -sigma * colvals[i];
               ++nconsvars;
            }
            if ( mipdata->yrhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->yrhs[pos];
               consvals[nconsvars] = sigma * colvals[i];
               ++nconsvars;
            }
            assert( nconsvars <= (int) mipdata->n );
         }

         /* check for lower and upper objective bounds */
         if ( (sepadata->useobjub || sepadata->useobjlb) && ! SCIPisZero(scip, SCIPcolGetObj(cols[j])) )
         {
            /* add lower objective bound */
            if ( mipdata->ylhs[mipdata->nrows] )
            {
               assert( sepadata->useobjlb );
               consvars[nconsvars] = mipdata->ylhs[mipdata->nrows];
               consvals[nconsvars] = -sigma * SCIPcolGetObj(cols[j]);
               ++nconsvars;
            }

            /* add upper objective bound */
            if ( mipdata->yrhs[mipdata->nrows] )
            {
               assert( sepadata->useobjub );
               consvars[nconsvars] = mipdata->yrhs[mipdata->nrows];
               consvals[nconsvars] = sigma * SCIPcolGetObj(cols[j]);
               ++nconsvars;
            }
         }

         /* add linear constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cont_%d", j);

         /* for free continuous variables require equality */
         r = SCIPinfinity(subscip);
         if ( SCIPisInfinity(scip, -lb[j]) && SCIPisInfinity(scip, ub[j]) )
            r = 0.0;
         else
            assert( SCIPisZero(scip, lb[j]) );

         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, nconsvars, consvars, consvals, 0.0, r,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         ++cnt;
      }
   }
   assert( (int) cnt <= ncols );
   mipdata->m += cnt;

   /* create constraints for rhs of cut */
   nconsvars = 0;

   /* first for the rows */
   for (i = 0; i < nrows; ++i)
   {
      assert( rows[i] != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(rows[i]) || (SCIProwIsLocal(rows[i]) && !sepadata->allowlocal) )
         continue;

      /* if lhs is there */
      if ( mipdata->ylhs[i] != NULL && ! SCIPisZero(scip, lhs[i]) )
      {
         assert( ! SCIPisInfinity(scip, -lhs[i]) );
         consvars[nconsvars] = mipdata->ylhs[i];
         consvals[nconsvars] = -lhs[i];
         ++nconsvars;
      }
      /* if rhs is there */
      if ( mipdata->yrhs[i] != NULL && ! SCIPisZero(scip, rhs[i]) )
      {
         assert( ! SCIPisInfinity(scip, rhs[i]) );
         consvars[nconsvars] = mipdata->yrhs[i];
         consvals[nconsvars] = rhs[i];
         ++nconsvars;
      }
      assert( nconsvars <= (int) mipdata->n );
   }

   if ( sepadata->useobjub || sepadata->useobjlb )
   {
      /* add lower objective bound */
      if ( mipdata->ylhs[mipdata->nrows] != NULL && ! SCIPisZero(scip, lhs[mipdata->nrows]) )
      {
         assert( sepadata->useobjlb );
         assert( ! SCIPisInfinity(scip, -lhs[mipdata->nrows]) );
         consvars[nconsvars] = mipdata->ylhs[mipdata->nrows];
         consvals[nconsvars] = -lhs[mipdata->nrows];
         ++nconsvars;
      }

      /* add upper objective bound */
      if ( mipdata->yrhs[mipdata->nrows] != NULL && ! SCIPisZero(scip, rhs[mipdata->nrows]) )
      {
         assert( sepadata->useobjub );
         assert( ! SCIPisInfinity(scip, rhs[mipdata->nrows]) );
         consvars[nconsvars] = mipdata->yrhs[mipdata->nrows];
         consvals[nconsvars] = rhs[mipdata->nrows];
         ++nconsvars;
      }
      assert( nconsvars <= (int) mipdata->n );
   }

   /* next for the columns */
   for (j = 0; j < ncols; ++j)
   {
      /* if ub is there */
      if ( mipdata->z[j] != NULL && ! SCIPisZero(scip, ub[j]) )
      {
         assert( mipdata->coltype[j] == colPresent );
         assert( ! SCIPisInfinity(scip, ub[j]) );
         consvars[nconsvars] = mipdata->z[j];
         consvals[nconsvars] = ub[j];
         ++nconsvars;
         assert( nconsvars <= (int) mipdata->n );
      }
   }
   /* add beta variable */
   consvars[nconsvars] = mipdata->beta;
   consvals[nconsvars] = -1.0;
   ++nconsvars;

   /* add fractional-beta variable */
   consvars[nconsvars] = mipdata->fracbeta;
   consvals[nconsvars] = -1.0;
   ++nconsvars;
   assert( nconsvars <= (int) mipdata->n );

   /* add linear constraint */
   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "beta", nconsvars, consvars, consvals, 0.0, 0.0,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   ++mipdata->m;

   /* add primal separation constraint if required */
   if ( sepadata->primalseparation )
   {
      SCIP_SOL* bestsol;
      bestsol = SCIPgetBestSol(scip);
      if ( bestsol != NULL )
      {
         nconsvars = 0;
         for (j = 0; j < ncols; ++j)
         {
            if ( mipdata->alpha[j] != NULL )
            {
               SCIP_Real val;
               assert( mipdata->coltype[j] == colPresent );

               val = SCIPgetSolVal(scip, bestsol, SCIPcolGetVar(cols[j]));
               consvars[nconsvars] = mipdata->alpha[j];
               consvals[nconsvars] = val;
               ++nconsvars;
               assert( nconsvars <= (int) mipdata->n );
            }
         }
         consvars[nconsvars] = mipdata->beta;
         consvals[nconsvars] = -1.0;
         ++nconsvars;

         /* add linear constraint - allow slight deviation from equality */
         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "primalseparation", nconsvars, consvars, consvals, -EPSILONVALUE, EPSILONVALUE,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         ++mipdata->m;
      }
   }

   /* add constraint to force violated cuts if required */
   if ( sepadata->addviolationcons )
   {
      nconsvars = 0;
      for (j = 0; j < ncols; ++j)
      {
         if ( mipdata->alpha[j] != NULL )
         {
            consvars[nconsvars] = mipdata->alpha[j];
            consvals[nconsvars] = primsol[j];
            ++nconsvars;
            assert( nconsvars <= (int) mipdata->n );
         }
      }
      consvars[nconsvars] = mipdata->beta;
      consvals[nconsvars] = -1.0;
      ++nconsvars;

      /* add linear constraint - allow slight deviation from equality */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "violationConstraint", nconsvars, consvars, consvals, MINEFFICACY, SCIPinfinity(subscip),
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      ++mipdata->m;
   }

   SCIPdebugMsg(scip, "Subscip has %u vars (%d integral, %d continuous), %u conss.\n",
      mipdata->n, SCIPgetNIntVars(subscip), SCIPgetNContVars(subscip), mipdata->m);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &consvals);

   SCIPfreeBufferArray(scip, &primsol);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);

   /* SCIPdebug( SCIP_CALL( SCIPprintOrigProblem(subscip, NULL, NULL, FALSE) ) ); */

#ifdef SCIP_WRITEPROB
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cgsepa%s%s%s%s_%s.lp",
         sepadata->objlone ? "_l1" : "",
         sepadata->addviolationcons ? "_vc" : "",
         sepadata->skipmultbounds ? "_ub" : "",
         sepadata->primalseparation ? "_ps" : "",
         SCIPgetProbName(scip));
      SCIP_CALL( SCIPwriteOrigProblem(subscip, name, "lp", FALSE) );
      SCIPinfoMessage(scip, NULL, "Wrote subscip to file <%s>.\n", name);
   }
#endif

   return SCIP_OKAY;
}


/** sets parameters for subscip */
static
SCIP_RETCODE subscipSetParams(
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_Bool*            success             /**< if setting was successful -> stop solution otherwise */
   )
{
   SCIP* subscip;

   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( success != NULL );

   *success = TRUE;

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   /* set objective limit, if no corresponding constraint has been added */
   if ( ! sepadata->addviolationcons && ! sepadata->addviolconshdlr )
   {
      SCIP_CALL( SCIPsetObjlimit(subscip, MINEFFICACY) );
   }

   /* do not abort subscip on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable memory saving mode: this is likely to result in the maximal depth being reached. This is because DFS
    * results in a repeated branching on the alpha-variables, which often have large bounds resulting in deep levels of
    * the tree. */
   SCIP_CALL( SCIPsetRealParam(subscip, "memory/savefac", 1.0) );

   /* set number of solutions stored */
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/maxsol", MAXNSOLS) );

   /* determine output to console */
#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1000) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/nsols/active", 2) );
#else
   if ( sepadata->output )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
      SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1000) );
      SCIP_CALL( SCIPsetIntParam(subscip, "display/nsols/active", 2) );
   }
   else
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   }
#endif

   if ( sepadata->subscipfast )
   {
      /* forbid recursive call of plugins solving subMIPs (also disables CG-separation) */
#ifdef SCIP_OUTPUT
      SCIP_CALL( SCIPsetSubscipsOff(subscip, FALSE) );
#else
      SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) ); /* quiet */
#endif
   }
   else
   {
      /* avoid recursive call */
      if ( ! SCIPisParamFixed(subscip, "separating/cgmip/freq") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "separating/cgmip/freq", -1) );
      }
   }

#ifdef SCIP_DISABLED_CODE
   /* the following possibly helps to improve performance (untested) */
   SCIP_CALL( SCIPsetEmphasis(subscip, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );
#else

   /* zirounding is often successful, so allow it some more calls */
   if( !SCIPisParamFixed(subscip, "heuristics/zirounding/minstopncalls") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/zirounding/minstopncalls", 10000) );
   }

   if ( sepadata->subscipfast )
   {
      /* set other heuristics */
      if( !SCIPisParamFixed(subscip, "heuristics/shifting/freq") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/shifting/freq", 3) );
      }
      if( !SCIPisParamFixed(subscip, "heuristics/simplerounding/freq") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/simplerounding/freq", 1) );
      }
      if( !SCIPisParamFixed(subscip, "heuristics/rounding/freq") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rounding/freq", 1) );
      }
      if( !SCIPisParamFixed(subscip, "heuristics/oneopt/freq") )
      {
         SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/oneopt/freq", 1) );
      }

      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/pscostdiving/freq", 1) ); */
      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/feaspump/freq", 3) ); */

      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/coefdiving/freq", -1) ); */
      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/fracdiving/freq", -1) ); */
      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/guideddiving/freq", -1) ); */
      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/linesearchdiving/freq", -1) ); */
      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/objpscostdiving/freq", -1) ); */
      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rootsoldiving/freq", -1) ); */
      /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/veclendiving/freq", -1) ); */

      /* use fast presolving */
      SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

      /* disable conflict analysis */
      /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) ); */
      /*     SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useinflp", 'o') ); */
      /*     SCIP_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') ); */
      /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) ); */
      /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) ); */

      /* use fast separation */
      SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_FAST, TRUE) );
   }
#endif

   return SCIP_OKAY;
}


/** solve subscip */
static
SCIP_RETCODE solveSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_Bool*            success             /**< if setting was successful -> stop solution otherwise */
   )
{
   SCIP* subscip;
   SCIP_RETCODE retcode;
   SCIP_STATUS status;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;
   SCIP_Longint nodelimit;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( success != NULL );

   *success = TRUE;

   subscip = mipdata->subscip;

   SCIP_CALL( SCIPcheckCopyLimits(scip, success) );

   if ( *success )
   {
      SCIP_CALL( SCIPcopyLimits(scip, subscip) );

      SCIP_CALL( SCIPgetRealParam(subscip, "limits/time", &timelimit) );
      SCIP_CALL( SCIPgetRealParam(subscip, "limits/memory", &memorylimit) );

      /* reduce time and memory limit if a smaller limit is stored in the separator data */
      if ( sepadata->timelimit < timelimit )
      {
         SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", sepadata->timelimit) );
      }
      if ( sepadata->memorylimit < memorylimit )
      {
         SCIP_CALL( SCIPsetRealParam(subscip, "limits/memorylimit", sepadata->memorylimit) );
      }
   }
   else
      return SCIP_OKAY;

   /* set nodelimit for subproblem */
   if ( sepadata->minnodelimit < 0 || sepadata->maxnodelimit < 0 )
      nodelimit = SCIP_LONGINT_MAX;
   else
   {
      assert( sepadata->minnodelimit >= 0 && sepadata->maxnodelimit >= 0 );
      nodelimit = SCIPgetNLPIterations(scip);
      nodelimit = MAX(sepadata->minnodelimit, nodelimit);
      nodelimit = MIN(sepadata->maxnodelimit, nodelimit);
   }
   assert( nodelimit >= 0 );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

   SCIPdebugMsg(scip, "Solving sub-SCIP (time limit: %f  mem limit: %f  node limit: %" SCIP_LONGINT_FORMAT ") ...\n", timelimit, memorylimit, nodelimit);

   /* disable statistic timing inside sub SCIP */
   if( !sepadata->output )
   {
      SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
   }

   /* check whether we want a complete solve */
   if ( ! sepadata->earlyterm )
   {
      retcode = SCIPsolve(subscip);

      /* errors in solving the subproblem should not kill the overall solving process;
       * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop. */
      if ( retcode != SCIP_OKAY )
      {
#ifndef NDEBUG
         SCIP_CALL( retcode );
#endif
         SCIPwarningMessage(scip, "Error while solving subproblem in CGMIP separator; sub-SCIP terminated with code <%d>\n", retcode);
         *success = FALSE;
         return SCIP_OKAY;
      }

      status = SCIPgetStatus(subscip);

#ifdef SCIP_OUTPUT
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#else
      if ( sepadata->output )
      {
         SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
      }
#endif

      /* if the solution process was terminated or the problem is infeasible (can happen because of violation constraint) */
      if ( status == SCIP_STATUS_TIMELIMIT || status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_INFORUNBD )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      /* all other statuses except optimal are invalid */
      if ( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_NODELIMIT )
      {
         SCIPerrorMessage("Solution of subscip for CG-separation returned with invalid status %d.\n", status);
         return SCIP_ERROR;
      }
   }
   else
   {
      /* otherwise we want a heuristic solve */

      /* -> solve until first solution is found */
      SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", 1) );
      retcode = SCIPsolve(subscip);
      SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", -1) );

      /* errors in solving the subproblem should not kill the overall solving process;
       * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop. */
      if ( retcode != SCIP_OKAY )
      {
#ifndef NDEBUG
         SCIP_CALL( retcode );
#endif
         SCIPwarningMessage(scip, "Error while solving subproblem in CGMIP separator; sub-SCIP terminated with code <%d>\n", retcode);
         *success = FALSE;
         return SCIP_OKAY;
      }

      status = SCIPgetStatus(subscip);

      /* if the solution process was terminated or the problem is infeasible (can happen because of violation constraint) */
      if ( status == SCIP_STATUS_TIMELIMIT || status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_NODELIMIT ||
           status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_INFORUNBD || status == SCIP_STATUS_MEMLIMIT )
      {
         /* output statistics before stopping */
#ifdef SCIP_OUTPUT
         SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#else
         if ( sepadata->output )
         {
            SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
         }
#endif
         *success = FALSE;
         return SCIP_OKAY;
      }

      /* all other statuses except optimal or bestsollimit are invalid - (problem cannot be infeasible) */
      if ( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_BESTSOLLIMIT )
      {
         SCIPerrorMessage("Solution of subscip for CG-separation returned with invalid status %d.\n", status);
         return SCIP_ERROR;
      }

      /* solve some more, if a feasible solution was found */
      if ( status == SCIP_STATUS_BESTSOLLIMIT )
      {
         SCIPdebugMsg(scip, "Continue solving separation problem ...\n");

         SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", STALLNODELIMIT) );
         retcode = SCIPsolve(subscip);
         SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", -1LL) );

         /* errors in solving the subproblem should not kill the overall solving process;
          * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop. */
         if ( retcode != SCIP_OKAY )
         {
#ifndef NDEBUG
            SCIP_CALL( retcode );
#endif
            SCIPwarningMessage(scip, "Error while solving subproblem in CGMIP separator; sub-SCIP terminated with code <%d>\n", retcode);
            *success = FALSE;
            return SCIP_OKAY;
         }

         status = SCIPgetStatus(subscip);
         assert( status != SCIP_STATUS_BESTSOLLIMIT );

#ifdef SCIP_OUTPUT
         SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#else
         if ( sepadata->output )
         {
            SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
         }
#endif

         /* if the solution process was terminated */
         if ( status == SCIP_STATUS_TIMELIMIT || status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_MEMLIMIT )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }

         /* all other statuses except optimal or bestsollimit are invalid */
         if ( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_STALLNODELIMIT && status != SCIP_STATUS_NODELIMIT )
         {
            SCIPerrorMessage("Solution of subscip for CG-separation returned with invalid status %d.\n", status);
            return SCIP_ERROR;
         }
      }
   }

   return SCIP_OKAY;
}

/** Computes cut from the given multipliers
 *
 *  Note that the cut computed here in general will not be the same as the one computed with the
 *  sub-MIP, because of numerical differences. Here, we only combine rows whose corresponding
 *  multiplier is positive w.r.t. the feasibility tolerance. In the sub-MIP, however, the rows are
 *  combined in any case. This makes a difference, if the coefficients in the matrix are large and
 *  hence yield a value that is larger than the tolerance.
 *
 *  Because of the transformations we have the following:
 *
 *  If variable \f$x_j\f$ was complemented, we have \f$x'_j = u_j - x_j\f$. If in the transformed
 *  system the lower bound is used, its corresponding multiplier is \f$y^T A'_j - \lfloor y^T A'_j
 *  \rfloor\f$, which corresponds to
 *  \f[
 *      y^T A'_j - \lfloor y^T A'_j \rfloor = - y^T A_j - \lfloor - y^T A_j \rfloor = - y^T A_j + \lceil y^T A_j \rceil
 *  \f]
 *  in the original system.
 *
 *  If such a variable was at its upper bound before the transformation, it is at its lower bound
 *  afterwards. Hence, its contribution to the cut is 0.
 *
 *  Note that if the original LP-solution does not satisfy some of the rows with equality the
 *  violation of the cut might be smaller than what is computed with the reduced sub-MIP.
 *
 *  Furthermore, note that if continuous variables have been shifted, the computed violated may be
 *  different as well, because the necessary changes in the lhs/rhs are not used here anymore.
 *
 *  @todo check if cut is correct if continuous variables have been shifted.
 */
static
SCIP_RETCODE computeCut(
   SCIP*                 scip,               /**< original scip */
   SCIP_SEPA*            sepa,               /**< separator */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< current solution for sub-MIP */
   SCIP_Real*            cutcoefs,           /**< coefficients of the cut */
   SCIP_Real*            cutrhs,             /**< rhs of the cut */
   SCIP_Bool*            localrowsused,      /**< pointer to store whether local rows were used in summation */
   SCIP_Bool*            localboundsused,    /**< pointer to store whether local bounds were used in summation */
   int*                  cutrank,            /**< pointer to store the cut rank */
   SCIP_Bool*            success             /**< whether we produced a valid cut */
   )
{
   SCIP* subscip;
   SCIP_VAR** vars;
   SCIP_ROW** rows;
   SCIP_COL** cols;
   SCIP_Real val;
   SCIP_Real maxabsweight;
   int nvars;
   int ncols;
   int nrows;
   int i;
   int j;

   assert( scip != NULL );
   assert( mipdata != NULL );
   assert( sepadata != NULL );
   assert( cutcoefs != NULL );
   assert( cutrhs != NULL );
   assert( localrowsused != NULL );
   assert( cutrank != NULL );
   assert( success != NULL );

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   /* get data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   assert( nrows == (int) mipdata->nrows );
   assert( ncols == (int) mipdata->ncols );

   /* initialize */
   *success = TRUE;
   *localrowsused = FALSE;
   *cutrank = 0;
   *localboundsused = FALSE;
   BMSclearMemoryArray(cutcoefs, nvars);
   *cutrhs = 0.0;

   /* find maximal absolute weight */
   maxabsweight = 0.0;
   for (i = 0; i < nrows; ++i)
   {
      SCIP_ROW* row;
      SCIP_Real absweight = 0.0;

      row = rows[i];
      assert( row != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
      {
         assert( mipdata->ylhs[i] == NULL && mipdata->yrhs[i] == NULL );
         continue;
      }

      /* get weight from solution (take larger of the values of lhs/rhs) */
      if ( mipdata->ylhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[i]);

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         if ( SCIPisFeasPositive(scip, val) )
            absweight = val;

         assert( ! sepadata->onlyrankone || SCIProwGetOriginSepa(row) != sepa );
      }
      if ( mipdata->yrhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[i]);

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, absweight) )
            absweight = val;

         assert( ! sepadata->onlyrankone || SCIProwGetOriginSepa(row) != sepa );
      }
      assert( ! SCIPisNegative(scip, absweight) );

      if ( absweight > maxabsweight )
         maxabsweight = absweight;
   }

   /* get weight from objective cuts */
   if ( sepadata->useobjub || sepadata->useobjlb )
   {
      SCIP_Real absweight = 0.0;

      assert( mipdata->ntotalrows == mipdata->nrows + 1 );

      if ( mipdata->ylhs[mipdata->nrows] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[mipdata->nrows]);
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         if ( SCIPisFeasPositive(scip, val) )
            absweight = val;
      }
      if ( mipdata->yrhs[mipdata->nrows] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[mipdata->nrows]);
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, absweight) )
            absweight = val;
      }

      if ( absweight > maxabsweight )
         maxabsweight = absweight;
   }

   /* calculate the row summation */
   for (i = 0; i < nrows; ++i)
   {
      SCIP_ROW* row;
      SCIP_Real weight;
      SCIP_Real absweight;
      SCIP_Bool uselhs;

      row = rows[i];
      assert( row != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
      {
         assert( mipdata->ylhs[i] == NULL && mipdata->yrhs[i] == NULL );
         continue;
      }

      /* get weight from solution */
      weight = 0.0;
      uselhs = FALSE;
      if ( mipdata->ylhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[i]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         if ( SCIPisFeasPositive(scip, val) )
         {
            uselhs = TRUE;
            weight = -val;
         }
      }
      if ( mipdata->yrhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[i]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, REALABS(weight)) )
            weight = val;
      }

      /* add row if weight is nonzero and lies within range */
      absweight = REALABS(weight);
      if ( ! SCIPisSumZero(scip, weight) && absweight * MAXWEIGHTRANGE >= maxabsweight )
      {
         SCIP_COL** rowcols;
         SCIP_Real* rowvals;

         rowcols = SCIProwGetCols(row);
         rowvals = SCIProwGetVals(row);

         /* add the row coefficients to the sum */
         for (j = 0; j < SCIProwGetNLPNonz(row); ++j)
         {
            int idx;
            SCIP_VAR* var;

            assert( rowcols[j] != NULL );
            var = SCIPcolGetVar(rowcols[j]);

            assert( var != NULL );
            assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );
            assert( SCIPvarGetCol(var) == rowcols[j] );

            idx = SCIPvarGetProbindex(var);
            assert( 0 <= idx && idx < nvars );

            cutcoefs[idx] += weight * rowvals[j];
         }

         /* compute rhs */
         if ( uselhs )
         {
            assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );
            val = SCIProwGetLhs(row) - SCIProwGetConstant(row);
            if ( SCIProwIsIntegral(row) )
               val = SCIPfeasCeil(scip, val); /* row is integral: round left hand side up */
         }
         else
         {
            assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );
            val = SCIProwGetRhs(row) - SCIProwGetConstant(row);
            if ( SCIProwIsIntegral(row) )
               val = SCIPfeasFloor(scip, val); /* row is integral: round right hand side down */
         }
         (*cutrhs) += weight * val;

         *localrowsused = *localrowsused || SCIProwIsLocal(row);

         if ( SCIProwGetRank(row) > *cutrank )
            *cutrank = SCIProwGetRank(row);
      }
   }
   /* add 1 to cutrank */
   ++(*cutrank);

   /* get weight from objective bounds */
   if ( sepadata->useobjub || sepadata->useobjlb )
   {
      SCIP_Real weight = 0.0;
      SCIP_Bool uselhs = FALSE;
      SCIP_Real absweight;

      assert( mipdata->ntotalrows == mipdata->nrows + 1 );

      if ( mipdata->ylhs[mipdata->nrows] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[mipdata->nrows]);
         assert( ! SCIPisFeasNegative(subscip, val) );
         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         if ( SCIPisFeasPositive(scip, val) )
         {
            uselhs = TRUE;
            weight = -val;
         }
      }
      if ( mipdata->yrhs[mipdata->nrows] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[mipdata->nrows]);
         assert( ! SCIPisFeasNegative(subscip, val) );
         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, REALABS(weight)) )
            weight = val;
      }

      /* add objective row if weight is nonzero and lies within range */
      absweight = REALABS(weight);
      if ( ! SCIPisSumZero(scip, weight) && absweight * MAXWEIGHTRANGE >= maxabsweight )
      {
         SCIP_Real obj = 0.0;

         /* add the objective row coefficients to the sum */
         for (j = 0; j < ncols; ++j)
         {
            obj = SCIPcolGetObj(cols[j]);
            if ( ! SCIPisZero(scip, obj) )
               cutcoefs[j] += weight * obj;
         }

         /* compute rhs */
         if ( uselhs )
         {
            val = SCIPgetLowerbound(scip);
            assert( ! SCIPisInfinity(scip, -val) );
            if ( SCIPisObjIntegral(scip) )
               val = SCIPfeasCeil(scip, val); /* objective is integral: round left hand side up */
         }
         else
         {
            val = SCIPgetUpperbound(scip);
            assert( ! SCIPisInfinity(scip, val) );
            if ( SCIPisObjIntegral(scip) )
               val = SCIPfeasFloor(scip, val); /* objective is integral: round right hand side down */
         }
         (*cutrhs) += weight * val;
      }
   }

   /* add upper bounds */
   for (j = 0; j < ncols; ++j)
   {
      assert( cols[j] != NULL );
      if ( mipdata->z[j] != NULL )
      {
         assert( mipdata->coltype[j] == colPresent );

         val = SCIPgetSolVal(subscip, sol, mipdata->z[j]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         /* if a bound has been used */
         if ( SCIPisSumPositive(subscip, val) )
         {
            SCIP_VAR* var;
            int idx;

            var = SCIPcolGetVar(cols[j]);

            assert( var != NULL );
            assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );
            assert( SCIPvarIsIntegral(var) );
            assert( SCIPvarGetCol(var) == cols[j] );

            idx = SCIPvarGetProbindex(var);
            assert( 0 <= idx && idx < nvars );

            /* check whether variable is complemented */
            if ( mipdata->iscomplemented[j] )
            {
               SCIP_Real lbnd;
               lbnd = SCIPvarGetLbGlobal(var);
               assert( ! SCIPisInfinity(scip, -lbnd) );
               assert( SCIPisIntegral(scip, lbnd) );
               assert( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPcolGetLb(cols[j])) );

               /* variable should not be free */
               assert( ! SCIPisInfinity(scip, -lbnd) || ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetLbLocal(var) - 0.5 > lbnd )
               {
                  lbnd = SCIPvarGetLbLocal(var);
                  assert( SCIPisIntegral(scip, lbnd) );
                  *localboundsused = TRUE;
               }

               cutcoefs[idx] -= val;
               *cutrhs -= lbnd * val;
            }
            else
            {
               SCIP_Real ubnd;
               ubnd = SCIPvarGetUbGlobal(var);
               assert( ! SCIPisInfinity(scip, ubnd) );
               assert( SCIPisIntegral(scip, ubnd) );
               assert( SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPcolGetUb(cols[j])) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetUbLocal(var) + 0.5 < ubnd )
               {
                  ubnd = SCIPvarGetUbLocal(var);
                  assert( SCIPisIntegral(scip, ubnd) );
                  *localboundsused = TRUE;
               }

               cutcoefs[idx] += val;
               *cutrhs += ubnd * val;
            }
         }
      }
   }

   /* check lower bounds for integral variables */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;
      int pos;

      var = vars[j];
      assert( var != NULL );
      pos = SCIPcolGetLPPos(SCIPvarGetCol(var));

      /* a variable may have status COLUMN, but the corresponding column may not (yet) be in the LP */
      if ( pos >= 0 && mipdata->coltype[pos] != colContinuous && mipdata->coltype[pos] != colConverted )
      {
         assert( pos < ncols );
         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );
         assert( SCIPvarIsIntegral(var) );

         /* check whether variable is complemented */
         if ( mipdata->iscomplemented[pos] )
         {
            assert( ! mipdata->isshifted[pos] );
            /* if the variable is complemented, the multiplier for the upper bound arises from the
               lower bound multiplier for the transformed problem - because of the minus-sign in the
               transformation this yields a round-up operation. */
            val = SCIPfeasCeil(scip, cutcoefs[j]) - cutcoefs[j];
            assert( ! SCIPisFeasNegative(scip, val) );

            /* only if variable needs to be rounded */
            if ( SCIPisSumPositive(scip, val) )
            {
               SCIP_Real ubnd;
               ubnd = SCIPvarGetUbGlobal(var);
               assert( ! SCIPisInfinity(scip, ubnd) );
               assert( SCIPisIntegral(scip, ubnd) );

               /* variable should not be free */
               assert( ! SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) || ! SCIPisInfinity(scip, ubnd) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetUbLocal(var) + 0.5 < ubnd )
               {
                  ubnd = SCIPvarGetUbLocal(var);
                  assert( SCIPisIntegral(scip, ubnd) );
                  *localboundsused = TRUE;
               }

               /* round cut coefficients, i.e., add val to cutcoefs[j] */
               cutcoefs[j] = SCIPfeasCeil(scip, cutcoefs[j]);

               /* correct rhs */
               if ( ! SCIPisSumZero(scip, ubnd) )
                  *cutrhs += ubnd * val;
            }
         }
         else
         {
            /* compute multiplier for lower bound: */
            val = cutcoefs[j] - SCIPfeasFloor(scip, cutcoefs[j]);
            assert( ! SCIPisFeasNegative(scip, val) );

            /* only if variable needs to be rounded */
            if ( SCIPisSumPositive(scip, val) )
            {
               SCIP_Real lbnd;
               lbnd = SCIPvarGetLbGlobal(var);
               assert( ! SCIPisInfinity(scip, -lbnd) );
               assert( SCIPisIntegral(scip, lbnd) );

               /* variable should not be free */
               assert( ! SCIPisInfinity(scip, -lbnd) || ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetLbLocal(var) - 0.5 > lbnd )
               {
                  lbnd = SCIPvarGetLbLocal(var);
                  assert( SCIPisIntegral(scip, lbnd) );
                  *localboundsused = TRUE;
               }

               /* round cut coefficients, i.e., subtract val from cutcoefs[j] */
               cutcoefs[j] = SCIPfeasFloor(scip, cutcoefs[j]);

               /* correct rhs */
               if ( ! SCIPisSumZero(scip, lbnd) )
                  *cutrhs -= lbnd * val;
            }
         }
      }
      else
      {
         /* force coefficients of all continuous variables or of variables not in the lp to zero */
         assert( pos == -1 || mipdata->coltype[pos] == colContinuous || mipdata->coltype[pos] == colConverted );

         /* check whether all coefficients for continuous or converted variables are nonnegative */
         if ( pos >= 0 )
         {
            if ( SCIPisNegative(scip, cutcoefs[j]) )
            {
               *success = FALSE;
               break;
            }
         }

         cutcoefs[j] = 0.0;
      }
   }

   /* round rhs */
   *cutrhs = SCIPfeasFloor(scip, *cutrhs);

   return SCIP_OKAY;
}

/** Create CG-cut directly from solution of sub-MIP */
static
SCIP_RETCODE createCGCutDirect(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_SOL*             sol,                /**< solution of sub-MIP */
   SCIP_Real*            cutcoefs,           /**< cut coefficients */
   int*                  cutinds,            /**< problem indices of variables appearing in cut */
   SCIP_Real*            cutvals,            /**< values of variables in cut */
   SCIP_Real*            varsolvals,         /**< solution value of variables */
   SCIP_Real*            weights,            /**< weights to compute cmir cut */
   int*                  nprevrows,          /**< number of previously generated rows */
   SCIP_ROW**            prevrows,           /**< previously generated rows */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   unsigned int*         ngen                /**< number of generated cuts */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Bool cutislocal;
   SCIP_Bool localrowsused;
   SCIP_Bool localboundsused;
   SCIP_Real cutrhs;
   SCIP_Real cutact;
   SCIP_Bool success;
   SCIP_VAR** vars;
   int cutrank = 0;
   int nvars;
   int k;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( sol != NULL );
   assert( cutcoefs != NULL );
   assert( cutinds != NULL );
   assert( cutvals != NULL );
   assert( varsolvals != NULL );
   assert( weights != NULL );
   assert( nprevrows != NULL );
   assert( prevrows != NULL );
   assert( cutoff != NULL );
   assert( ngen != NULL );

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   cutrhs = 0.0;
   localrowsused = FALSE;
   localboundsused = FALSE;
   *cutoff = FALSE;
   success = TRUE;

   /* compute coefficients */
   SCIP_CALL( computeCut(scip, sepa, mipdata, sepadata, sol, cutcoefs, &cutrhs, &localrowsused, &localboundsused, &cutrank, &success) );
   cutislocal = localrowsused || localboundsused;

   /* take next solution if cut was not valid */
   if ( ! success )
   {
      SCIPdebugMsg(scip, "cut not valid - skipping ...\n");
      return SCIP_OKAY;
   }

   /* compute activity */
   cutact = 0.0;
   for (k = 0; k < nvars; ++k)
      cutact += cutcoefs[k] * varsolvals[k];

#ifdef SCIP_DISABLED_CODE
   /* the following test should be treated with care because of numerical differences - see computeCut() */
   {
      /* check for correctness of computed values */
      SCIP* subscip;
      SCIP_Real obj = 0.0;
      SCIP_Real val;
      SCIP_Bool contVarShifted = FALSE;
      unsigned int j;
      SCIP_COL** cols;
      int ncols;

      subscip = mipdata->subscip;
      assert( subscip != NULL );

      SCIP_CALL( SCIPprintSol(subscip, sol, NULL, FALSE) );

      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      for (j = 0; j < mipdata->ncols; ++j)
      {
         if ( mipdata->coltype[j] == colPresent )
         {
            int idx;
            assert( mipdata->alpha[j] != NULL );
            val = SCIPgetSolVal(subscip, sol, mipdata->alpha[j]);
            assert( SCIPisFeasIntegral(subscip, val) );
            idx = SCIPvarGetProbindex(SCIPcolGetVar(cols[j]));
            assert( SCIPisFeasEQ(scip, val, cutcoefs[idx]) );
            obj += val * SCIPvarGetObj(mipdata->alpha[j]);
         }
         else
         {
            if ( (mipdata->coltype[j] == colContinuous || mipdata->coltype[j] == colConverted) && mipdata->isshifted[j] )
               contVarShifted = TRUE;
         }
      }
      assert( mipdata->beta != NULL );
      val = SCIPgetSolVal(subscip, sol, mipdata->beta);
      assert( SCIPisFeasIntegral(subscip, val) );
      obj += val * SCIPvarGetObj(mipdata->beta);
      assert( contVarShifted || SCIPisFeasEQ(scip, obj, cutact - cutrhs) );
   }
#endif

   /* if successful, convert dense cut into sparse row, and add the row as a cut */
   if ( SCIPisFeasGT(scip, cutact, cutrhs) )
   {
      SCIP_Real cutnorm;
      int cutlen;

      /* store the cut as sparse row, calculate activity and norm of cut */
      SCIP_CALL( storeCutInArrays(scip, nvars, cutcoefs, varsolvals, mipdata->normtype,
            cutinds, cutvals, &cutlen, &cutact, &cutnorm) );

      SCIPdebugMsg(scip, "act=%f, rhs=%f, norm=%f, eff=%f\n", cutact, cutrhs, cutnorm, (cutact - cutrhs)/cutnorm);

      /* if norm is 0, the cut is trivial */
      if ( SCIPisPositive(scip, cutnorm) )
      {
         SCIP_Bool violated = SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm);

         if ( violated || (sepadata->usecutpool && ! cutislocal ) )
         {
            SCIP_ROW* cut;

            /* create the cut */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cgcut%d_%u", SCIPgetNLPs(scip), *ngen);
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, name, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

            for( k = 0; k < cutlen; ++k )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[k]], cutvals[k]) );
            }

            /* set cut rank */
            SCIProwChgRank(cut, cutrank);

            SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

            /*SCIPdebug( SCIP_CALL( SCIPprintRow(scip, cut, NULL) ) );*/

            /* add cut to pool */
            if ( ! cutislocal )
            {
               assert( violated || sepadata->usecutpool );
               SCIP_CALL( SCIPaddPoolCut(scip, cut) );
            }

            /* add cut if it is violated */
            if ( violated )
            {
               /* check whether cut has been found before - may happend due to projection */
               for (k = 0; k < *nprevrows; ++k)
               {
                  SCIP_Real parval;

                  assert( prevrows[k] != NULL );
                  parval = SCIProwGetParallelism(cut, prevrows[k], 'e');
                  /* exit if row is parallel to existing cut and rhs is not better */
                  if ( SCIPisEQ(scip, parval, 1.0) && SCIPisGE(scip, cutrhs, SCIProwGetRhs(prevrows[k])) )
                     break;
               }

               /* if cut is new */
               if ( k >= *nprevrows )
               {
                  prevrows[*nprevrows] = cut;
                  ++(*nprevrows);

                  SCIPdebugMsg(scip, " -> CG-cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                     name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                     SCIPgetCutEfficacy(scip, NULL, cut),
                     SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                     SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
#ifdef SCIP_DEBUG
                  SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#else
                  if ( sepadata->output )
                  {
                     SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
                  }
#endif
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
                  ++(*ngen);
               }
               else
               {
                  SCIPdebugMsg(scip, "Cut already exists.\n");
                  /* release the row */
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }
            }
            else
            {
               /* release the row */
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** create CG-cut via CMIR-function */
static
SCIP_RETCODE createCGCutCMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_SOL*             sol,                /**< solution of sub-MIP */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row to use for creating MIR cut */
   SCIP_Real*            cutcoefs,           /**< cut coefficients */
   int*                  cutinds,            /**< problem indices of variables appearing in cut */
   SCIP_Real*            cutvals,            /**< values of variables in cut */
   SCIP_Real*            varsolvals,         /**< solution value of variables */
   SCIP_Real*            weights,            /**< weights to compute cmir cut */
   int*                  boundsfortrans,     /**< bounds for cmir function of NULL */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds for cmir function or NULL */
   int*                  nprevrows,          /**< number of previously generated rows */
   SCIP_ROW**            prevrows,           /**< previously generated rows */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   unsigned int*         ngen                /**< number of generated cuts */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Longint maxdnom;
   SCIP_Bool cutislocal;
   SCIP_Real maxscale;
   SCIP_Real cutrhs;
   SCIP_Real cutefficacy;
   SCIP_Bool success;
   SCIP_ROW** rows;
   SCIP_VAR** vars;
   SCIP* subscip;
   int nrows;
   int nvars;
   int k;
   int cutrank;
   int cutnnz;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( sol != NULL );
   assert( cutcoefs != NULL );
   assert( cutinds != NULL );
   assert( cutvals != NULL );
   assert( varsolvals != NULL );
   assert( weights != NULL );
   assert( nprevrows != NULL );
   assert( prevrows != NULL );
   assert( cutoff != NULL );
   assert( ngen != NULL );

   *cutoff = FALSE;
   subscip = mipdata->subscip;
   assert( subscip != NULL );

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert( nrows > 0 );
   assert( (int) mipdata->nrows == nrows );

   /* @todo more advanced settings - compare sepa_gomory.c */
   maxdnom = (SCIP_Longint) sepadata->cutcoefbnd+1;
   maxscale = 10000.0;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* generate weights */
   for (k = 0; k < nrows; ++k)
   {
      SCIP_Real val;

      weights[k] = 0;
      if ( mipdata->ylhs[k] != NULL )
      {
         assert( !SCIProwIsModifiable(rows[k]) && (!SCIProwIsLocal(rows[k]) || sepadata->allowlocal) );

         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[k]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         if ( SCIPisFeasPositive(subscip, val) )
            weights[k] = -val;
      }
      if ( mipdata->yrhs[k] != NULL )
      {
         assert( !SCIProwIsModifiable(rows[k]) && (!SCIProwIsLocal(rows[k]) || sepadata->allowlocal) );

         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[k]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, ABS(weights[k])) )
            weights[k] = val;
      }
   }

   /* set up data for bounds to use */
   if ( sepadata->cmirownbounds )
   {
      int typefortrans;

      assert( boundsfortrans != NULL );
      assert( boundtypesfortrans != NULL );

      if ( sepadata->allowlocal )
         typefortrans = -2;
      else
         typefortrans = -1;

      /* check all variables */
      for (k = 0; k < nvars; ++k)
      {
         int pos;
         SCIP_VAR* var;

         var = vars[k];
         assert( var != NULL );
         pos = SCIPcolGetLPPos(SCIPvarGetCol(var));

         if ( pos < 0 )
            continue;

         assert( pos < (int) mipdata->ncols );
         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );

         boundsfortrans[k] = typefortrans;
         boundtypesfortrans[k] = SCIP_BOUNDTYPE_LOWER;

         if ( mipdata->coltype[pos] == colContinuous || mipdata->coltype[pos] == colConverted )
         {
            assert( SCIPvarIsIntegral(var) || mipdata->coltype[pos] != colContinuous );
            continue;
         }

         /* check upper bound */
         if ( mipdata->z[pos] != NULL && SCIPisSumPositive(subscip, SCIPgetSolVal(subscip, sol, mipdata->z[pos])) )
         {
            /* check whether variable is complemented */
            if ( ! mipdata->iscomplemented[pos] )
               boundtypesfortrans[k] = SCIP_BOUNDTYPE_UPPER;
            /* otherwise use lower bound */
         }
         else
         {
            /* check whether variable is complemented */
            if ( mipdata->iscomplemented[pos] )
               boundtypesfortrans[k] = SCIP_BOUNDTYPE_UPPER;
            /* otherwise use lower bound */
         }
      }
   }

   /* create a MIR cut using the above calculated weights */
   cutefficacy = -1.0;
   cutrhs = -1.0;
   SCIP_CALL( SCIPaggrRowSumRows(scip, aggrrow, weights, NULL, -1, FALSE,
         sepadata->allowlocal, 2, (int) MAXAGGRLEN(nvars), &success) );

   if( !success )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcalcMIR(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, sepadata->allowlocal, FIXINTEGRALRHS, boundsfortrans,
         boundtypesfortrans, MINFRAC, MAXFRAC, 1.0, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy,
         &cutrank, &cutislocal, &success) );

   assert( sepadata->allowlocal || !cutislocal );
   SCIPdebugMsg(scip, "CMIR: success = %u, cut is%sefficacious (cutefficacy: %g, cutrhs: %g)\n", success,
      SCIPisEfficacious(scip, cutefficacy) ? " " : " not ", cutefficacy, cutrhs);

   /* if successful, convert dense cut into sparse row, and add the row as a cut
    * only if the cut if violated - if it is not violated we might store non-local cuts in the pool
    */
   if( success && (SCIPisEfficacious(scip, cutefficacy) || (sepadata->usecutpool && ! cutislocal)) )
   {
      SCIP_ROW* cut;

      /* create the cut */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cgcut%d_%u", SCIPgetNLPs(scip), *ngen);
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, name, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );

      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

      for( k = 0; k < cutnnz; ++k )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[k]], cutcoefs[k]) );
      }

      assert( success );

      /* set cut rank */
      SCIProwChgRank(cut, cutrank);

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#else
      if ( sepadata->output )
      {
         SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
      }
#endif

      /* try to scale the cut to integral values */
      SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                                     maxdnom, maxscale, MAKECONTINTEGRAL, &success) );

      /* if the cut could be made integral */
      if ( success )
      {
         SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

         /* add cut to pool */
         if ( ! cutislocal )
         {
            assert( SCIPisEfficacious(scip, cutefficacy) || sepadata->usecutpool );
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         }

         if ( ! SCIPisCutEfficacious(scip, NULL, cut) )
         {
            SCIPdebugMsg(scip, " -> CG-cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
                         name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                         SCIPgetCutEfficacy(scip, NULL, cut));

            /* release the row */
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
         else
         {
            /* check whether cut has been found before - may happend due to projection */
            for (k = 0; k < *nprevrows; ++k)
            {
               SCIP_Real parval;

               assert( prevrows[k] != NULL );
               parval = SCIProwGetParallelism(cut, prevrows[k], 'e');
               /* exit if row is parallel to existing cut and rhs is not better */
               if ( SCIPisEQ(scip, parval, 1.0) && SCIPisGE(scip, cutrhs, SCIProwGetRhs(prevrows[k])) )
                  break;
            }

            /* if cut is new */
            if ( k >= *nprevrows )
            {
               prevrows[*nprevrows] = cut;
               ++(*nprevrows);

               SCIPdebugMsg(scip, " -> CG-cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, rank=%d, min=%f, max=%f (range=%f)\n",
                            name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                            SCIPgetCutEfficacy(scip, NULL, cut), SCIProwGetRank(cut),
                            SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                            SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
               #ifdef SCIP_OUTPUT
               SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
               #else
               if ( sepadata->output )
               {
                  SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
               }
               #endif
               SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
               ++(*ngen);
            }
            else
            {
               SCIPdebugMsg(scip, "Cut already exists.\n");
               /* release the row */
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            }
         }
      }
      else
      {
         SCIPdebugMsg(scip, " -> CG-cut <%s> could not be scaled to integral coefficients: rhs=%f, eff=%f\n",
                      name, cutefficacy, cutrhs);

         /* release the row */
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      }
   }

   return SCIP_OKAY;
}


/** create CG-cut via strong-CG-function */
static
SCIP_RETCODE createCGCutStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_SOL*             sol,                /**< solution of sub-MIP */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row to use for creating MIR cut */
   SCIP_Real*            cutcoefs,           /**< cut coefficients */
   int*                  cutinds,            /**< problem indices of variables appearing in cut */
   SCIP_Real*            cutvals,            /**< values of variables in cut */
   SCIP_Real*            varsolvals,         /**< solution value of variables */
   SCIP_Real*            weights,            /**< weights to compute cmir cut */
   int*                  nprevrows,          /**< number of previously generated rows */
   SCIP_ROW**            prevrows,           /**< previously generated rows */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   unsigned int*         ngen                /**< number of generated cuts */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Longint maxdnom;
   SCIP_Bool cutislocal;
   SCIP_Real maxscale;
   SCIP_Real cutrhs;
   SCIP_Real cutefficacy;
   SCIP_Bool success;
   SCIP_ROW** rows;
   SCIP_VAR** vars;
   SCIP* subscip;
   int nrows;
   int nvars;
   int k;
   int cutrank;
   int cutnnz;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( sol != NULL );
   assert( cutcoefs != NULL );
   assert( cutinds != NULL );
   assert( cutvals != NULL );
   assert( varsolvals != NULL );
   assert( weights != NULL );
   assert( nprevrows != NULL );
   assert( prevrows != NULL );
   assert( cutoff != NULL );
   assert( ngen != NULL );

   *cutoff = FALSE;
   subscip = mipdata->subscip;
   assert( subscip != NULL );

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert( nrows > 0 );
   assert( (int) mipdata->nrows == nrows );

   /* @todo more advanced settings - compare sepa_gomory.c */
   maxdnom = (SCIP_Longint) sepadata->cutcoefbnd + 1;
   maxscale = 10000.0;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* generate weights */
   for (k = 0; k < nrows; ++k)
   {
      SCIP_Real val;

      weights[k] = 0;
      if ( mipdata->ylhs[k] != NULL )
      {
         assert( !SCIProwIsModifiable(rows[k]) && (!SCIProwIsLocal(rows[k]) || sepadata->allowlocal) );

         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[k]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         if ( SCIPisFeasPositive(subscip, val) )
            weights[k] = -val;
      }
      if ( mipdata->yrhs[k] != NULL )
      {
         assert( !SCIProwIsModifiable(rows[k]) && (!SCIProwIsLocal(rows[k]) || sepadata->allowlocal) );

         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[k]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         assert( sepadata->skipmultbounds || SCIPisFeasLT(subscip, val, 1.0) );
         val = SCIPfrac(scip, val);  /* take fractional value if variable has no upper bounds */

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, ABS(weights[k])) )
            weights[k] = val;
      }
   }

   /* create a strong CG cut out of the weighted LP rows using the B^-1 row as weights */
   cutefficacy = -1.0;
   cutrhs = -1.0;
   SCIP_CALL( SCIPaggrRowSumRows(scip, aggrrow, weights, NULL, -1, FALSE,
         sepadata->allowlocal, 1, (int) MAXAGGRLEN(nvars), &success) );

   if( !success )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcalcStrongCG(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, sepadata->allowlocal, MINFRAC, MAXFRAC,
         1.0, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank, &cutislocal, &success) );

   assert( sepadata->allowlocal || !cutislocal );
   SCIPdebugMsg(scip, "Strong-CG: success = %u, cut is%sefficacious (cutefficacy: %g, cutrhs: %g)\n", success,
      SCIPisEfficacious(scip, cutefficacy) ? " " : " not ", cutefficacy, cutrhs);

   /* if successful, convert dense cut into sparse row, and add the row as a cut
    * only if the cut if violated - if it is not violated we might store non-local cuts in the pool
    */
   if( success && (SCIPisEfficacious(scip, cutefficacy) || (sepadata->usecutpool && ! cutislocal)) )
   {
      SCIP_ROW* cut;

      /* create the cut */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cgcut%d_%u", SCIPgetNLPs(scip), *ngen);
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, name, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );

      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

      for( k = 0; k < cutnnz; ++k )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[k]], cutcoefs[k]) );
      }

      assert( success );

      /* set cut rank */
      SCIProwChgRank(cut, cutrank);

#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#else
      if ( sepadata->output )
      {
         SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
      }
#endif

      /* try to scale the cut to integral values */
      SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                                     maxdnom, maxscale, MAKECONTINTEGRAL, &success) );

      /* if the cut could be made integral */
      if ( success )
      {
         SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

         /* add cut to pool */
         if ( ! cutislocal )
         {
            assert( SCIPisEfficacious(scip, cutefficacy) || sepadata->usecutpool );
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
         }

         if ( ! SCIPisCutEfficacious(scip, NULL, cut) )
         {
            SCIPdebugMsg(scip, " -> CG-cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
                         name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                         SCIPgetCutEfficacy(scip, NULL, cut));

            /* release the row */
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
         else
         {
            /* check whether cut has been found before - may happend due to projection */
            for (k = 0; k < *nprevrows; ++k)
            {
               SCIP_Real parval;

               assert( prevrows[k] != NULL );
               parval = SCIProwGetParallelism(cut, prevrows[k], 'e');
               /* exit if row is parallel to existing cut and rhs is not better */
               if ( SCIPisEQ(scip, parval, 1.0) && SCIPisGE(scip, cutrhs, SCIProwGetRhs(prevrows[k])) )
                  break;
            }

            /* if cut is new */
            if ( k >= *nprevrows )
            {
               prevrows[*nprevrows] = cut;
               ++(*nprevrows);

               SCIPdebugMsg(scip, " -> CG-cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, rank=%d, min=%f, max=%f (range=%f)\n",
                            name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                            SCIPgetCutEfficacy(scip, NULL, cut), SCIProwGetRank(cut),
                            SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                            SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
               #ifdef SCIP_OUTPUT
               SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
               #else
               if ( sepadata->output )
               {
                  SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
               }
               #endif
               SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
               ++(*ngen);
            }
            else
            {
               SCIPdebugMsg(scip, "Cut already exists.\n");
               /* release the row */
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            }
         }
      }
      else
      {
         SCIPdebugMsg(scip, " -> CG-cut <%s> could not be scaled to integral coefficients: rhs=%f, eff=%f\n",
                      name, cutefficacy, cutrhs);

         /* release the row */
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      }
   }

   return SCIP_OKAY;
}


/** Create CG-cuts from solutions of sub-MIP */
static
SCIP_RETCODE createCGCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   unsigned int*         ngen                /**< number of generated cuts */
   )
{
   SCIP_BOUNDTYPE* boundtypesfortrans;
   SCIP_STAGE stage;
   SCIP_AGGRROW* aggrrow = NULL;
   SCIP_Real* varsolvals;
   SCIP_Real* weights;
   int* cutinds;
   SCIP_Real* cutvals;
   SCIP_Real* cutcoefs;
   SCIP_ROW** prevrows;
   SCIP_SOL** sols;
   SCIP_VAR** vars;
   SCIP* subscip;
   int* boundsfortrans;
   int nprevrows;
   int ntotalrows;
   int nsols;
   int nvars;
   int k;
   int s;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( cutoff != NULL );
   assert( ngen != NULL );

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   *cutoff = FALSE;
   *ngen = 0;

   /* check if solving was successful and get solutions */
   stage = SCIPgetStage(subscip);
   if ( stage == SCIP_STAGE_SOLVING || stage == SCIP_STAGE_SOLVED )
      nsols = SCIPgetNSols(subscip);
   else
      nsols = 0;

   /* only if solutions have been found */
   if ( nsols == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Generating CG-cuts from %d sols (cmir: %u, strong-cg: %u) ...\n", nsols, sepadata->usecmir, sepadata->usestrongcg);

   sols = SCIPgetSols(subscip);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* allocate temporary memory */
   assert(mipdata->ntotalrows <= INT_MAX);
   ntotalrows = (int)mipdata->ntotalrows;
   assert( ntotalrows >= SCIPgetNLPRows(scip) && ntotalrows <= SCIPgetNLPRows(scip) + 1 );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, ntotalrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prevrows, 2 * nsols) );

   if ( sepadata->usecmir || sepadata->usestrongcg )
   {
      SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );
   }

   /* prepare arrays for bound information, if requested */
   if ( sepadata->usecmir && sepadata->cmirownbounds )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &boundsfortrans, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundtypesfortrans, nvars) );
   }
   else
   {
      boundsfortrans = NULL;
      boundtypesfortrans = NULL;
   }

   /* get solution values */
   for (k = 0; k < nvars; ++k)
   {
      if ( SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_COLUMN )
         varsolvals[k] = SCIPvarGetLPSol(vars[k]);
      else
         varsolvals[k] = 0.0;
   }

   /* loop through solutions found */
   nprevrows = 0;
   for (s = 0; s < nsols; ++s)
   {
      SCIP_SOL* sol;
      sol = sols[s];

      /* generate cuts by the C-MIR and/or Strong-CG functions */
      if ( sepadata->usecmir )
      {
         SCIP_CALL( createCGCutCMIR(scip, sepa, sepadata, mipdata, sol, aggrrow, cutcoefs, cutinds, cutvals, varsolvals, weights,
               boundsfortrans, boundtypesfortrans, &nprevrows, prevrows, cutoff, ngen) );
      }

      if ( sepadata->usestrongcg )
      {
         SCIP_CALL( createCGCutStrongCG(scip, sepa, sepadata, mipdata, sol, aggrrow, cutcoefs, cutinds, cutvals, varsolvals, weights,
               &nprevrows, prevrows, cutoff, ngen) );
      }

      if ( ! sepadata->usecmir && ! sepadata->usestrongcg )
      {
         SCIP_CALL( createCGCutDirect(scip, sepa, sepadata, mipdata, sol, cutcoefs, cutinds, cutvals, varsolvals, weights,
               &nprevrows, prevrows, cutoff, ngen) );

         assert(! sepadata->usecmir && ! sepadata->usestrongcg);
      }
   }
   assert( nprevrows <= 2 * nsols );
   assert( sepadata->usecmir || nprevrows <= nsols );
   assert( sepadata->usestrongcg || nprevrows <= nsols );

   /* release rows */
   for (k = 0; k < nprevrows; ++k)
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(prevrows[k])) );
   }

   if ( sepadata->usecmir || sepadata->usestrongcg )
      SCIPaggrRowFree(scip, &aggrrow);

   /* free temporary memory */
   SCIPfreeBufferArrayNull(scip, &boundsfortrans);
   SCIPfreeBufferArrayNull(scip, &boundtypesfortrans);

   SCIPfreeBufferArray(scip, &prevrows);
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &cutvals);
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &cutcoefs);

   return SCIP_OKAY;
}


/** frees "subscip" data */
static
SCIP_RETCODE freeSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator data */
   CGMIP_MIPDATA*        mipdata             /**< data for sub-MIP */
   )
{
   SCIP_SEPADATA* sepadata;
   unsigned int i, j;
   SCIP* subscip;

   assert( scip != NULL );
   assert( sepa != NULL );
   assert( mipdata != NULL );

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIPdebugMsg(scip, "Freeing subscip ...\n");

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   for (j = 0; j < mipdata->ncols; ++j)
   {
      if ( mipdata->coltype[j] == colPresent )
      {
         assert( mipdata->alpha[j] != NULL );
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->alpha[j])) );
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->fracalpha[j])) );
      }
   }
   SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->beta)) );
   SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->fracbeta)) );

   for (i = 0; i < mipdata->nrows; ++i)
   {
      if ( mipdata->ylhs[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->ylhs[i])) );
      }
      if ( mipdata->yrhs[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->yrhs[i])) );
      }
   }

   if ( sepadata->useobjub || sepadata->useobjlb )
   {
      if ( mipdata->yrhs[mipdata->nrows] )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->yrhs[mipdata->nrows])) );
      }
      if ( mipdata->ylhs[mipdata->nrows] )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->ylhs[mipdata->nrows])) );
      }
   }

   for (j = 0; j < mipdata->ncols; ++j)
   {
      if ( mipdata->z[j] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->z[j])) );
      }
   }

   SCIP_CALL( SCIPfree(&(mipdata->subscip)) );

   SCIPfreeBlockMemoryArray(scip, &(mipdata->z), 2*mipdata->ncols); /*lint !e647*/
   SCIPfreeBlockMemoryArray(scip, &(mipdata->yrhs), mipdata->ntotalrows);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->ylhs), mipdata->ntotalrows);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->isshifted), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->iscomplemented), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->coltype), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->fracalpha), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->alpha), mipdata->ncols);

   return SCIP_OKAY;
}


/*
 * Callback methods
 */


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitCGMIP)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* create and initialize random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, DEFAULT_RANDSEED) );

   return SCIP_OKAY;
}

/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitCGMIP)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeRandom(scip, &sepadata->randnumgen);

   return SCIP_OKAY;
}

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyCGMIP)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaCGMIP(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeCGMIP)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIPfreeBlockMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpCGMIP)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   CGMIP_MIPDATA* mipdata;

   int depth;
   int ncalls;
   int ncols;
   int nrows;
   unsigned int ngen;
   SCIP_Bool success;
   SCIP_Bool cutoff = FALSE;

   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;
   ngen = 0;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);

   /* only call separator, if we are not close to terminating */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call separator up to a maximum depth */
   if ( sepadata->maxdepth >= 0 && depth > sepadata->maxdepth )
      return SCIP_OKAY;

   /* only call separator a given number of times at each node */
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   if ( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* skip separation if there are continuous variables, but only integers required */
   if ( SCIPgetNContVars(scip) > 0 && sepadata->onlyintvars )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* check for parameters */
   if ( ( sepadata->useobjub || sepadata->useobjlb ) && ( sepadata->usecmir || sepadata->usestrongcg ) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
         "WARNING - sepa_cgmip: Using objective function bounds and CMIR or Strong-CG functions is useless. Turning off usage of objective function bounds.\n");
      SCIP_CALL( SCIPsetBoolParam(scip, "separating/cgmip/useobjub", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "separating/cgmip/useobjlb", FALSE) );
   }
   sepadata->allowlocal = allowlocal;

   /* get LP data */
   ncols = SCIPgetNLPCols(scip);
   nrows = SCIPgetNLPRows(scip);
   if ( ncols <= NCOLSTOOSMALL || nrows <= NROWSTOOSMALL )
      return SCIP_OKAY;

   /* determine whether we should run the separation based on a decision tree */
   if ( sepadata->decisiontree )
   {
      SCIP_Bool separate;
      SCIP_Real firstlptime;

      separate = FALSE;
      firstlptime = SCIPgetFirstLPTime(scip);

      if ( nrows <= 136 && firstlptime <= 0.05 && ncols <= 143 )
         separate = TRUE;
      else if ( nrows <= 136 && 0.05 < firstlptime && firstlptime <= 0.15 && ncols <= 143 )
         separate = TRUE;
      else if ( 136 < nrows && nrows <= 332 && ncols <= 143 )
         separate = TRUE;
      else if ( 136 < nrows && nrows <= 332 && 655 < ncols && ncols <= 1290 )
         separate = TRUE;
      else if ( 333 < nrows && nrows <= 874 && 0.15 < firstlptime && firstlptime <= 0.25 && 2614 < ncols && ncols <= 5141 )
         separate = TRUE;
      else if ( 875 < nrows && nrows <= 1676 && firstlptime <= 0.05 && 143 < ncols && ncols <= 265 )
         separate = TRUE;
      else if ( 875 < nrows && nrows <= 1676 && firstlptime <= 0.05 && 265 < ncols && ncols <= 654 )
         separate = TRUE;
      else if ( 875 < nrows && nrows <= 1676 && 0.05 < firstlptime && firstlptime <= 0.15 )
         separate = TRUE;
      else if ( 875 < nrows && nrows <= 1676 &&  0.15 < firstlptime && firstlptime <= 0.25 && 1291 < ncols && ncols <= 2613 )
         separate = TRUE;
      else if ( nrows > 8146 && 0.75 < firstlptime && firstlptime <= 6.25 && 655 < ncols && ncols <= 1290 )
         separate = TRUE;
      else if ( nrows > 8146 && 0.75 < firstlptime && firstlptime <= 6.25 && 1291 < ncols && ncols <= 2613 )
         separate = TRUE;
      else if ( nrows > 8146 && firstlptime > 6.25 )
         separate = TRUE;

      if ( ! separate )
      {
         return SCIP_OKAY;
      }
   }

   /* preceed with separation */
   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "separating CG-cuts via sub-MIPs: %d cols, %d rows\n", ncols, nrows);

   /* prepare data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &mipdata) );
   mipdata->subscip = NULL;
   mipdata->alpha = NULL;
   mipdata->fracalpha = NULL;
   mipdata->beta = NULL;
   mipdata->fracbeta = NULL;
   mipdata->coltype = NULL;
   mipdata->iscomplemented = NULL;
   mipdata->isshifted = NULL;
   mipdata->ylhs = NULL;
   mipdata->yrhs = NULL;
   mipdata->z = NULL;
   mipdata->normtype = ' ';

   mipdata->conshdlrfullnorm = CONSHDLRFULLNORM;
   mipdata->scip = scip;
   mipdata->sepa = sepa;
   mipdata->sepadata = sepadata;

   /* get the type of norm to use for efficacy calculations */
   SCIP_CALL( SCIPgetCharParam(scip, "separating/efficacynorm", &mipdata->normtype) );

   /* create subscip */
   SCIP_CALL( createSubscip(scip, sepa, sepadata, mipdata) );

   /* set parameters */
   SCIP_CALL( subscipSetParams(sepadata, mipdata, &success) );

   if ( success && !SCIPisStopped(scip) )
   {
      /* solve subscip */
      SCIP_CALL( solveSubscip(scip, sepadata, mipdata, &success) );

      /* preceed if solution was successful */
      if ( success && ! SCIPisStopped(scip) )
      {
         SCIP_CALL( createCGCuts(scip, sepa, sepadata, mipdata, &cutoff, &ngen) );
      }
   }

   SCIP_CALL( freeSubscip(scip, sepa, mipdata) );
   SCIPfreeBlockMemory(scip, &mipdata);

   SCIPdebugMsg(scip, "Found %u CG-cuts.\n", ngen);

   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ngen > 0 )
      *result = SCIP_SEPARATED;

#ifdef SCIP_OUTPUT
   /* SCIP_CALL( SCIPwriteLP(scip, "cuts.lp") ); */
   /* SCIP_CALL( SCIPwriteMIP(scip, "cuts.lp", FALSE, TRUE) ); */
#endif

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the CGMIP MIR cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaCGMIP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );

   sepa = NULL;
   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpCGMIP, NULL, sepadata) );
   assert(sepa != NULL);

   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyCGMIP) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeCGMIP) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitCGMIP) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitCGMIP) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxrounds",
         "maximal number of cgmip separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxroundsroot",
         "maximal number of cgmip separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxdepth",
         "maximal depth at which the separator is applied (-1: unlimited)",
         &sepadata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/decisiontree",
         "Use decision tree to turn separation on/off?",
         &sepadata->decisiontree, FALSE, DEFAULT_DECISIONTREE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/timelimit",
         "time limit for sub-MIP",
         &sepadata->timelimit, TRUE, DEFAULT_TIMELIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/memorylimit",
         "memory limit for sub-MIP",
         &sepadata->memorylimit, TRUE, DEFAULT_MEMORYLIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "separating/" SEPA_NAME "/minnodelimit",
         "minimum number of nodes considered for sub-MIP (-1: unlimited)",
         &sepadata->minnodelimit, FALSE, DEFAULT_MINNODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "separating/" SEPA_NAME "/maxnodelimit",
         "maximum number of nodes considered for sub-MIP (-1: unlimited)",
         &sepadata->maxnodelimit, FALSE, DEFAULT_MAXNODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/cutcoefbnd",
         "bounds on the values of the coefficients in the CG-cut",
         &sepadata->cutcoefbnd, TRUE, DEFAULT_CUTCOEFBND, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/onlyactiverows",
         "Use only active rows to generate cuts?",
         &sepadata->onlyactiverows, FALSE, DEFAULT_ONLYACTIVEROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxrowage",
         "maximal age of rows to consider if onlyactiverows is false",
         &sepadata->maxrowage, FALSE, DEFAULT_MAXROWAGE, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/onlyrankone",
         "Separate only rank 1 inequalities w.r.t. CG-MIP separator?",
         &sepadata->onlyrankone, FALSE, DEFAULT_ONLYRANKONE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/onlyintvars",
         "Generate cuts for problems with only integer variables?",
         &sepadata->onlyintvars, FALSE, DEFAULT_ONLYINTVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/contconvert",
         "Convert some integral variables to be continuous to reduce the size of the sub-MIP?",
         &sepadata->contconvert, FALSE, DEFAULT_CONTCONVERT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/contconvfrac",
         "fraction of integral variables converted to be continuous (if contconvert)",
         &sepadata->contconvfrac, FALSE, DEFAULT_CONTCONVFRAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/contconvmin",
         "minimum number of integral variables before some are converted to be continuous",
         &sepadata->contconvmin, FALSE, DEFAULT_CONTCONVMIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/intconvert",
         "Convert some integral variables attaining fractional values to have integral value?",
         &sepadata->intconvert, FALSE, DEFAULT_INTCONVERT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/intconvfrac",
         "fraction of frac. integral variables converted to have integral value (if intconvert)",
         &sepadata->intconvfrac, FALSE, DEFAULT_INTCONVFRAC, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/intconvmin",
         "minimum number of integral variables before some are converted to have integral value",
         &sepadata->intconvmin, FALSE, DEFAULT_INTCONVMIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/skipmultbounds",
         "Skip the upper bounds on the multipliers in the sub-MIP?",
         &sepadata->skipmultbounds, FALSE, DEFAULT_SKIPMULTBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/objlone",
         "Should the objective of the sub-MIP minimize the l1-norm of the multipliers?",
         &sepadata->objlone, FALSE, DEFAULT_OBJLONE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/objweight",
         "weight used for the row combination coefficient in the sub-MIP objective",
         &sepadata->objweight, TRUE, DEFAULT_OBJWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/objweightsize",
         "Weight each row by its size?",
         &sepadata->objweightsize, FALSE, DEFAULT_OBJWEIGHTSIZE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/usecmir",
         "use CMIR-generator (otherwise add cut directly)?",
         &sepadata->usecmir, FALSE, DEFAULT_USECMIR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/usestrongcg",
         "use strong CG-function to strengthen cut?",
         &sepadata->usestrongcg, FALSE, DEFAULT_USESTRONGCG, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/cmirownbounds",
         "tell CMIR-generator which bounds to used in rounding?",
         &sepadata->cmirownbounds, FALSE, DEFAULT_CMIROWNBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/usecutpool",
         "use cutpool to store CG-cuts even if the are not efficient?",
         &sepadata->usecutpool, FALSE, DEFAULT_USECUTPOOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/primalseparation",
         "only separate cuts that are tight for the best feasible solution?",
         &sepadata->primalseparation, FALSE, DEFAULT_PRIMALSEPARATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/earlyterm",
         "terminate separation if a violated (but possibly sub-optimal) cut has been found?",
         &sepadata->earlyterm, FALSE, DEFAULT_EARLYTERM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/addviolationcons",
         "add constraint to subscip that only allows violated cuts (otherwise add obj. limit)?",
         &sepadata->addviolationcons, FALSE, DEFAULT_ADDVIOLATIONCONS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/addviolconshdlr",
         "add constraint handler to filter out violated cuts?",
         &sepadata->addviolconshdlr, FALSE, DEFAULT_ADDVIOLCONSHDLR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/conshdlrusenorm",
         "should the violation constraint handler use the norm of a cut to check for feasibility?",
         &sepadata->conshdlrusenorm, FALSE, DEFAULT_CONSHDLRUSENORM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/useobjub",
         "Use upper bound on objective function (via primal solution)?",
         &sepadata->useobjub, FALSE, DEFAULT_USEOBJUB, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/useobjlb",
         "Use lower bound on objective function (via primal solution)?",
         &sepadata->useobjlb, FALSE, DEFAULT_USEOBJLB, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/subscipfast",
         "Should the settings for the sub-MIP be optimized for speed?",
         &sepadata->subscipfast, FALSE, DEFAULT_SUBSCIPFAST, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/output",
         "Should information about the sub-MIP and cuts be displayed?",
         &sepadata->output, FALSE, DEFAULT_OUTPUT, NULL, NULL) );

   return SCIP_OKAY;
}
