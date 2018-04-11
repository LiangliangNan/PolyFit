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

/**@file   cons_orbisack.c
 * @brief  constraint handler for orbisack constraints
 * @author Christopher Hojny
 *
 *
 * The type of constraints of this constraint handler is described in cons_orbisack.h.
 *
 * The details of the method implemented here are described in the following papers:
 *
 * Describing Orbitopes by Linear Inequalities and Projection Based Tools@n
 * Andreas Loos,@n
 * PhD thesis, Otto-von-Guericke-Universitaet Magdeburg (2010).
 *
 * This thesis provides a complete linear description of orbisacks and a separation
 * routine for its inequalities.
 *
 * Polytopes Associated with Symmetry Handling@n
 * Christopher Hojny and Marc E. Pfetsch,@n
 * (2017), preprint available at http://www.optimization-online.org/DB_HTML/2017/01/5835.html
 *
 * This paper describes a linear time separation routine for so-called cover inequalities of
 * orbisacks.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_orbisack.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_setppc.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "orbisack"
#define CONSHDLR_DESC          "symmetry breaking constraint handler for orbisacks"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -1005200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -1005200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             5 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             5 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_EXHAUSTIVE

/* default parameters for separation routines: */
#define DEFAULT_ORBISEPARATION               FALSE     /**< whether orbisack inequalities should be separated */
#define DEFAULT_COVERSEPARATION               TRUE     /**< whether cover inequalities should be separated */

/* default parameters for constraints */
#define DEFAULT_COEFFBOUND               1000000.0     /**< maximum size of coefficients in orbisack inequalities */

#define DEFAULT_PPORBISACK        TRUE /**< whether we allow upgrading to packing/partitioning orbisacks */
#define DEFAULT_CHECKALWAYSFEAS    TRUE /**< whether check routine returns always SCIP_FEASIBLE */


/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             coverseparation;    /**< whether only cover inequalities should be separated */
   SCIP_Bool             orbiseparation;     /**< whether orbisack as well as cover inequalities should be separated */
   SCIP_Real             coeffbound;         /**< maximum size of coefficients in orbisack inequalities */
   SCIP_Bool             checkpporbisack;    /**< whether we allow upgrading to packing/partitioning orbisacks */
   SCIP_Bool             checkalwaysfeas;    /**< whether check routine returns always SCIP_FEASIBLE */
};

/** constraint data for orbisack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars1;              /**< first column of variable matrix */
   SCIP_VAR**            vars2;              /**< second column of variable matrix */
   int                   nrows;              /**< number of rows of variable matrix */
};


/*
 * Local methods
 */

/** frees orbisack constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to orbisack constraint data */
   )
{
   int nrows;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   nrows = (*consdata)->nrows;
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars2), nrows);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars1), nrows);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates orbisack constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_VAR*const*       vars1,              /**< first column of variable matrix */
   SCIP_VAR*const*       vars2,              /**< second column of variable matrix */
   int                   nrows               /**< number of rows in variable matrix */
   )
{
   int i;

   assert( consdata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars1, vars1, nrows) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars2, vars2, nrows) );

#ifndef NDEBUG
   {
      for (i = 0; i < nrows; ++i)
      {
         assert( SCIPvarIsBinary(vars1[i]) );
         assert( SCIPvarIsBinary(vars2[i]) );
      }
   }
#endif

   (*consdata)->nrows = nrows;

   /* get transformed variables, if we are in the transformed problem */
   if ( SCIPisTransformed(scip) )
   {
      /* Make sure that all variables cannot be multiaggregated (cannot be handled by cons_orbisack, since one cannot
       * easily eliminate single variables from an orbisack constraint. */
      for (i = 0; i < nrows; ++i)
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vars1[i], &(*consdata)->vars1[i]) );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars1[i]) );

         SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vars2[i], &(*consdata)->vars2[i]) );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars2[i]) );
      }
   }

   return SCIP_OKAY;
}


/** check wether an orbisack is even a packing/partitioning orbisack */
static
SCIP_RETCODE packingUpgrade(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< variables of second column */
   int                   nrows,              /**< number of rows of orbisack */
   SCIP_Bool*            success,            /**< memory address to store whether constraint can be upgraded */
   SCIP_Bool*            isparttype          /**< memory address to store whether upgraded orbisack is partitioning orbisack */
   )
{
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   SCIP_Bool* rowcovered;
   int nsetppcconss;
   int i;
   int c;

   assert( scip != NULL );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( success != NULL );
   assert( isparttype != NULL );

   *success = FALSE;

   /* get data of setppc conshdlr */
   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if ( setppcconshdlr == NULL )
   {
      SCIPwarningMessage(scip, "Check for upgrading orbisacks to packing/partitioning orbisacks not possible - setppc constraint handler not found.\n");
      return SCIP_OKAY;
   }
   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);

   /* upgrade cannot be successful */
   if ( nsetppcconss == 0 )
      return SCIP_OKAY;
   assert( setppcconss != NULL );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &rowcovered, nrows) );

   /* iterate over orbisack rows and check whether rows are contained in partitioning constraints */
   *success = TRUE;
   for (i = 0; i < nrows; ++i)
   {
      /* iterate over constraints */
      for (c = 0; c < nsetppcconss; ++c)
      {
         int nsetppcvars;
         SCIP_VAR** setppcvars;
         SCIP_VAR* var;
         int varidx1;
         int varidx2;

         /* check type */
         if ( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_COVERING ||
              SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PACKING )
            continue;
         assert( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PARTITIONING );

         /* get set packing/partitioning variables */
         nsetppcvars = SCIPgetNVarsSetppc(scip, setppcconss[c]);
         assert( nsetppcvars > 0 );

         /* partitioning constraint contains too much variables */
         if ( nsetppcvars != 2 )
            continue;
         assert( nsetppcvars == 2 );

         setppcvars = SCIPgetVarsSetppc(scip, setppcconss[c]);
         assert( setppcvars != NULL );

         /* check whether i-th row is contained in partitioning constraint */
         var = setppcvars[0];
         if ( SCIPvarIsNegated(var) )
            continue;

         varidx1 = SCIPvarGetProbindex(var);

         var = setppcvars[1];
         if ( SCIPvarIsNegated(var) )
            continue;

         varidx2 = SCIPvarGetProbindex(var);

         if ( (varidx1 == SCIPvarGetProbindex(vars1[i]) && varidx2 == SCIPvarGetProbindex(vars2[i])) ||
              (varidx2 == SCIPvarGetProbindex(vars1[i]) && varidx1 == SCIPvarGetProbindex(vars2[i])))
         {
            rowcovered[i] = TRUE;
            break;
         }
      }

      /* no partitioning constraint corresponds to row i */
      if ( c >= nsetppcconss )
         *success = FALSE;
   }

   if ( *success )
   {
      *isparttype = TRUE;
      SCIPfreeBufferArray(scip, &rowcovered);

      return SCIP_OKAY;
   }

   /* Iterate over orbisack rows and check whether rows are contained in packing constraints.
    * In particular, it is possible that variables in row i form a subset of variables in the
    * a packing/partitioning constraint. */
   *success = TRUE;
   for (i = 0; i < nrows; ++i)
   {
      /* skip already covered rows */
      if ( rowcovered[i] )
         continue;

      /* iterate over constraints */
      for (c = 0; c < nsetppcconss; ++c)
      {
         SCIP_VAR** setppcvars;
         SCIP_VAR* var;
         int nsetppcvars;
         int varidx;
         int nfound = 0;
         int j;

         /* check type */
         if ( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_COVERING )
            continue;
         assert( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PARTITIONING ||
            SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PACKING );

         /* get set packing/partitioning variables */
         nsetppcvars = SCIPgetNVarsSetppc(scip, setppcconss[c]);
         assert( nsetppcvars > 0 );

         setppcvars = SCIPgetVarsSetppc(scip, setppcconss[c]);
         assert( setppcvars != NULL );

         /* check whether all variables of the cycle are contained in setppc constraint */
         for (j = 0; j < nsetppcvars && nfound < 2; ++j)
         {
            var = setppcvars[j];

            if ( SCIPvarIsNegated(var) )
               continue;

            varidx = SCIPvarGetProbindex(var);

            if ( varidx == SCIPvarGetProbindex(vars1[i]) || varidx == SCIPvarGetProbindex(vars2[i]) )
               ++nfound;
         }

         if ( nfound == 2 )
            break;
      }

      /* row i is not contained in an setppc constraint */
      if ( c >= nsetppcconss )
         *success = FALSE;
   }

   /* we have found a packing orbisack */
   if ( *success )
      *isparttype = FALSE;

   SCIPfreeBufferArray(scip, &rowcovered);

   return SCIP_OKAY;
}


/** generate initial LP cut
 *
 *  We generate the inequality of the orbisack on the elements of the first row, i.e.,
 *  the inequality \f$-x_{1,1} + x_{1,2} \leq 0\f$.
 */
static
SCIP_RETCODE initLP(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_VAR* tmpvars[2];
   SCIP_ROW* row;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != 0 );
   assert( consdata->nrows > 0 );
   assert( consdata->vars1 != NULL );
   assert( consdata->vars2 != NULL );

   vars1 = consdata->vars1;
   vars2 = consdata->vars2;

   tmpvars[0] = vars1[0];
   tmpvars[1] = vars2[0];

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "orbisack0#0", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, tmpvars[0], -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, tmpvars[1], 1.0) );

   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** add orbisack cover inequality */
static
SCIP_RETCODE addOrbisackCover(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nrows,              /**< number of rows of orbisack */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< variables of second column */
   SCIP_Real*            coeffs1,            /**< coefficients of the variables of the first column of the inequality to be added */
   SCIP_Real*            coeffs2,            /**< coefficients of the variables of the second column of the inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( coeffs1 != NULL );
   assert( coeffs2 != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "orbisackcover", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
   for (i = 0; i < nrows; ++i)
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars1[i], coeffs1[i]) );
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars2[i], coeffs2[i]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** Separate lifted orbisack cover inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 *
 *  We iterate over the nrows-many cover inequalities which are potentially
 *  maximal w.r.t. their violation.
 */
static
SCIP_RETCODE separateOrbisackCovers(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nrows,              /**< number of rows of orbisack */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< variables of second column */
   SCIP_Real*            vals1,              /**< LP-solution for those variables in first column */
   SCIP_Real*            vals2,              /**< LP-solution for those variables in second column */
   int*                  ngen,               /**< number of separated covers */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_Real rhs = 0.0;
   SCIP_Real lhs = 0.0;
   SCIP_Real* coeff1;
   SCIP_Real* coeff2;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nrows > 0 );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   *infeasible = FALSE;
   *ngen = 0;

   /* allocate memory for inequality coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff2, nrows) );

   /* initialize coefficient matrix */
   for (i = 0; i < nrows; ++i)
   {
      coeff1[i] = 0.0;
      coeff2[i] = 0.0;
   }

   /* detect violated covers */
   for (i = 0; i < nrows; ++i)
   {
      /* cover inequality is violated */
      if ( SCIPisEfficacious(scip, -vals1[i] + vals2[i] + lhs - rhs) )
      {
         /* set coefficients for inequality */
         coeff1[i] = -1.0;
         coeff2[i] = 1.0;

         SCIP_CALL( addOrbisackCover(scip, cons, nrows, vars1, vars2, coeff1, coeff2, rhs, infeasible) );
         ++(*ngen);
         if ( *infeasible )
            break;

         /* reset coefficients for next inequality */
         coeff1[i] = 0.0;
         coeff2[i] = 0.0;
      }

      /* add argmax( 1 - vals[i][0], vals[i][1] ) as coefficient and ensure that both vars1[0] and vars2[0] are
       * contained in the LIFTED cover inequality */
      if ( SCIPisEfficacious(scip, 1.0 - vals1[i] - vals2[i]) )
      {
         coeff1[i] = -1.0;
         lhs = lhs - vals1[i];

         /* lifting */
         if ( i == 0 )
         {
            coeff2[0] = 1.0;
            lhs += vals2[i];
         }
      }
      else
      {
         coeff2[i] = 1.0;
         rhs += 1.0;
         lhs = lhs + vals2[i];

         /* lifting */
         if ( i == 0 )
         {
            coeff1[0] = -1.0;
            lhs -= vals1[i];
            rhs -= 1.0;
         }
      }
   }

   /* free coefficient matrix */
   SCIPfreeBufferArray(scip, &coeff2);
   SCIPfreeBufferArray(scip, &coeff1);

   return SCIP_OKAY;
}


/** add orbisack inequality */
static
SCIP_RETCODE addOrbisackInequality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nrows,              /**< number of rows of orbisack */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< variables of second column */
   SCIP_Real*            coeffs1,            /**< first column of coefficient matrix of inequality to be added */
   SCIP_Real*            coeffs2,            /**< second column of coefficient matrix of inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( coeffs1 != NULL );
   assert( coeffs2 != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "orbisack", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for (i = 0; i < nrows; ++i)
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars1[i], coeffs1[i]) );
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars2[i], coeffs2[i]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** separate orbisack inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 *
 *  We stop if we checked for each possible basement row, whether a cut could be added. If the coefficients grow too
 *  large, we start separating cover inequalities.
 *
 *  We implement the separation algorithm for orbisacks described in@n
 *  A. Loos. Describing Orbitopes by Linear Inequalities and Projection Based Tools.
 *  PhD thesis, Otto-von-Guericke-Universitaet Magdeburg, 2010.
 */
static
SCIP_RETCODE separateOrbisack(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nrows,              /**< number of rows of orbisack */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< variables of second column */
   SCIP_Real*            vals1,              /**< LP-solution for those variables in first column */
   SCIP_Real*            vals2,              /**< LP-solution for those variables in second column */
   SCIP_Bool             coverseparation,    /**< whether we separate cover inequalities */
   SCIP_Real             coeffbound,         /**< maximum size of coefficients in orbisack inequalities */
   int*                  ngen,               /**< pointer to store the number of generated cuts */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_Real* coeff1;
   SCIP_Real* coeff2;
   SCIP_Real rhs;
   SCIP_Real lhs;
   SCIP_Real valueA;
   SCIP_Real valueB;
   SCIP_Real valueC;
   int basement;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nrows > 0 );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( coeffbound >= 0.0 );
   assert( ngen != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;
   *ngen = 0;

   /* if there is only one row, all cuts are added by initLP */
   if ( nrows < 2 )
      return SCIP_OKAY;

   /* allocate memory for inequality coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff2, nrows) );

   /* initialize coefficient matrix row 0 */
   coeff1[0] = -1.0;
   coeff2[0] = 1.0;
   for (i = 2; i < nrows; ++i)
   {
      coeff1[i] = 0.0;
      coeff2[i] = 0.0;
   }

   /* initialize right-hand side and left-hand side (lhs for row 0) */
   rhs = 0.0;
   lhs = - vals1[0] + vals2[0];

   /* basement row of orbisack */
   basement = 1;

   /* update value of left-hand side and coefficients for basement row = 1 */
   lhs += - vals1[1] + vals2[1];
   coeff1[1] = -1.0;
   coeff2[1] = 1.0;

   /* check whether cut for basement row = 1 is violated */
   if ( SCIPisEfficacious(scip, lhs - rhs) )
   {
      SCIP_CALL( addOrbisackInequality(scip, cons, nrows, vars1, vars2, coeff1, coeff2, rhs, infeasible) );
      ++(*ngen);
   }

   /* check whether there exists a cut with basement rows > 1 that is violated */
   while ( basement < nrows - 1 && ! *infeasible )
   {
      valueA = lhs + vals1[basement] - vals1[basement + 1] + vals2[basement + 1] - rhs - 1.0; /*lint !e679, !e834*/
      valueB = lhs - vals2[basement] - vals1[basement + 1] + vals2[basement + 1] - rhs; /*lint !e679, !e834*/
      valueC = 2.0 * lhs + vals1[basement] - vals2[basement] - vals1[basement + 1] + vals2[basement + 1] - 2.0 * rhs; /*lint !e679, !e834*/

      /* update inequality */
      if ( valueA >= valueB && valueA >= valueC )
      {
         ++rhs;
         coeff1[basement] = 0.0;
         lhs += vals1[basement++];
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs += - vals1[basement] + vals2[basement];
      }
      else if ( valueB >= valueA && valueB >= valueC )
      {
         coeff2[basement] = 0.0;
         lhs -= vals2[basement++];
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs += - vals1[basement] + vals2[basement];
      }
      else
      {
         rhs *= 2.0;
         lhs = 0.0;
         for (i = 0; i < basement; ++i)
         {
            coeff1[i] = 2.0 * coeff1[i];
            coeff2[i] = 2.0 * coeff2[i];
            lhs += coeff1[i] * vals1[i] + coeff2[i] * vals2[i];
         }
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs -= vals1[basement];
         lhs += vals2[basement++];
         coeff1[basement] = -1.0;
         coeff2[basement] = 1.0;
         lhs -= vals1[basement];
         lhs += vals2[basement];
      }

      /* to avoid numerical troubles, we bound the size of coefficients and rhs */
      if ( rhs > coeffbound || -coeff1[0] > coeffbound || coeff2[0] > coeffbound )
      {
         /* avoid separating cover inequalities twice */
         if ( ! coverseparation )
         {
            int ncuts;
            SCIP_CALL( separateOrbisackCovers(scip, cons, nrows, vars1, vars2, vals1, vals2, &ncuts, infeasible) );
            *ngen += ncuts;
         }
         break;
      }

      /* if current inequality is violated */
      if ( SCIPisEfficacious(scip, lhs - rhs) )
      {
         SCIP_CALL( addOrbisackInequality(scip, cons, nrows, vars1, vars2, coeff1, coeff2, rhs, infeasible) );
         ++(*ngen);
      }
   }

   /* free allocated memory */
   SCIPfreeBufferArray(scip, &coeff2);
   SCIPfreeBufferArray(scip, &coeff1);

   return SCIP_OKAY;
}


/** propagation */
static
SCIP_RETCODE propVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   SCIP_Bool*            found,              /**< pointer to store whether a new propagation could be found */
   int*                  ngen                /**< pointer to store the number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   SCIP_Bool tightened = FALSE;
   SCIP_VAR* var;
   SCIP_Real ub;
   SCIP_Real lb;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   int* solu1;
   int* solu2;
   int nrows;
   int i;
   int r;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );
   assert( found != NULL );

   SCIPdebugMsg(scip, "Propagating variables of constraint <%s>.\n", SCIPconsGetName(cons));

   *ngen = 0;
   *infeasible = FALSE;
   *found = FALSE;

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars1 != NULL );
   assert( consdata->vars2 != NULL );
   assert( consdata->nrows > 0 );

   nrows = consdata->nrows;
   vars1 = consdata->vars1;
   vars2 = consdata->vars2;

   /* determine current solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &solu1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solu2, nrows) );

   for (i = 0; i < nrows; ++i)
   {
      /* determine value in first column */
      var = vars1[i];
      assert( SCIPvarIsBinary(var) );

      /* get local upper and lower bound on variable */
      ub = SCIPvarGetUbLocal(var);
      lb = SCIPvarGetLbLocal(var);

      /* if variable is fixed to 1 -> solu[i][j] = 1,
       * if it is fixed to 0 -> solu[i][j] = 0,
       * else, -> solu[i][j] = 2 */
      if ( lb > 0.5 )
         solu1[i] = 1;
      else if (ub < 0.5)
         solu1[i] = 0;
      else
         solu1[i] = 2;

      /* determine value in second column */
      var = vars2[i];
      assert( SCIPvarIsBinary(var) );

      /* get local upper and lower bound on variable */
      ub = SCIPvarGetUbLocal(var);
      lb = SCIPvarGetLbLocal(var);

      /* if variable is fixed to 1 -> solu[i][j] = 1,
       * if it is fixed to 0 -> solu[i][j] = 0,
       * else, -> solu[i][j] = 2 */
      if ( lb > 0.5 )
         solu2[i] = 1;
      else if (ub < 0.5)
         solu2[i] = 0;
      else
         solu2[i] = 2;
   }

   /* loop through all variables */
   for (i = 0; i < nrows; ++i)
   {
      /* get variables of first and second column */
      var1 = vars1[i];
      var2 = vars2[i];
      assert( var1 != NULL );
      assert( var2 != NULL );

      /* if variable in first column is fixed to 0 and variable in second column is fixed to 1 */
      if ( SCIPvarGetUbLocal(var1) < 0.5 && SCIPvarGetLbLocal(var2) > 0.5 )
      {
         SCIP_Bool nocritical = TRUE;

         SCIPdebugMsg(scip, "Check variable pair (%d,0) and (%d,1).\n", i, i);

         /* check whether there is a critical row above row i, otherwise the solution is infeasible
          * if there is a row without fixing (2) above the current row, we cannot obtain a result */
         for (r = 0; r < i; ++r)
         {
            if ( (solu1[r] == 1 && solu2[r] == 0) || solu1[r] == 2 || solu2[r] == 2 )
            {
               nocritical = FALSE;
               break;
            }
         }

         if ( nocritical )
         {
            SCIPdebugMsg(scip, " -> node infeasible (row was fixed to 0,1 but there was no critical row above).\n");

            /* perform conflict analysis */
            if ( SCIPisConflictAnalysisApplicable(scip) )
            {
               SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

               for (r = 0; r <= i; ++r)
               {
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars1[r]) );
                  SCIP_CALL( SCIPaddConflictBinvar(scip, vars2[r]) );
               }

               SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

               *infeasible = TRUE;

               /* free current solution */
               SCIPfreeBufferArray(scip, &solu2);
               SCIPfreeBufferArray(scip, &solu1);

               return SCIP_OKAY;
            }
         }
      }
      /* if variable in the first column is fixed to 0 and the variable in the second column is free */
      else if ( SCIPvarGetUbLocal(var1) < 0.5 && SCIPvarGetUbLocal(var2) > 0.5 )
      {
         SCIP_Bool allconstant = TRUE;

         SCIPdebugMsg(scip, "Check variable pair (%d,0) and (%d,1).\n", i, i);

         /* Check whether all rows above row i are constant. In this case, the variable in the second */
         /* column can be fixed to 0. If an entry above row i is unfixed or a row is not constant, we cannot */
         /* fix the second entry in row i. */
         for (r = 0; r < i; ++r)
         {
            if ( solu1[r] == 2 || solu2[r] == 2 || solu1[r] != solu2[r] )
            {
               allconstant = FALSE;
               break;
            }
         }

         /* fix variable in the second column to 0 */
         if ( allconstant )
         {
            assert( SCIPvarGetLbLocal(var2) < 0.5 );
            SCIP_CALL( SCIPinferVarUbCons(scip, var2, 0.0, cons, i, FALSE, infeasible, &tightened) ); /*lint !e713*/
            assert( ! *infeasible );

            *found = *found || tightened;
            if ( tightened )
               ++(*ngen);
         }
      }
   }

   /* free current solution */
   SCIPfreeBufferArray(scip, &solu2);
   SCIPfreeBufferArray(scip, &solu1);

   return SCIP_OKAY;
}


/** separate orbisack and cover inequalities */
static
SCIP_RETCODE separateInequalities(
   SCIP*                 scip,               /**< pointer to scip */
   SCIP_RESULT*          result,             /**< pointer to store the result of separation */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nrows,              /**< number of rows of orbisack */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< variables of second column */
   SCIP_Real*            vals1,              /**< LP-solution for those variables in first column */
   SCIP_Real*            vals2               /**< LP-solution for those variables in second column */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible = FALSE;
   int ngen1 = 0;
   int ngen2 = 0;

   assert( scip != NULL );
   assert( result != NULL );
   assert( cons != NULL );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( vals1 != NULL );
   assert( vals2 != NULL );

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->orbiseparation )
   {
      SCIP_CALL( separateOrbisack(scip, cons, nrows, vars1, vars2, vals1, vals2, FALSE, conshdlrdata->coeffbound, &ngen1, &infeasible) );
   }

   if ( ! infeasible && conshdlrdata->coverseparation )
   {
      SCIP_CALL( separateOrbisackCovers(scip, cons, nrows, vars1, vars2, vals1, vals2, &ngen2, &infeasible) );
   }

   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if ( ngen1 + ngen2 > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/*--------------------------------------------------------------------------------------------
 *--------------------------------- SCIP functions -------------------------------------------
 *--------------------------------------------------------------------------------------------*/

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteOrbisack)
{  /*lint --e{715}*/
   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( consdata != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeOrbisack)
{   /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransOrbisack)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* consdata = NULL;
   int nrows;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMsg(scip, "Transforming constraint.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->nrows > 0 );
   assert( sourcedata->vars1 != NULL );
   assert( sourcedata->vars2 != NULL );

   /* create transformed constraint data (copy data where necessary) */
   nrows = sourcedata->nrows;

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nrows = nrows;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars1, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars2, nrows) );

   SCIP_CALL( SCIPgetTransformedVars(scip, nrows, sourcedata->vars1, consdata->vars1) );
   SCIP_CALL( SCIPgetTransformedVars(scip, nrows, sourcedata->vars2, consdata->vars2) );

   /* create transformed constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpOrbisack)
{
   int c;

   assert( infeasible != NULL );
   *infeasible = FALSE;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != 0 );

      SCIPdebugMsg(scip, "Generating initial orbisack cut for constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( initLP(scip, conss[c], infeasible) );
      if ( *infeasible )
         break;

      SCIPdebugMsg(scip, "Generated initial orbisack cut.\n");
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for orbisack constraints.\n");

   *result = SCIP_DIDNOTRUN;

   /* if solution is not integer */
   if ( SCIPgetNLPBranchCands(scip) > 0 )
   {
      *result = SCIP_DIDNOTFIND;

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);

         /* get solution */
         SCIP_CALL( SCIPallocBufferArray(scip, &vals1, consdata->nrows) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals2, consdata->nrows) );

         SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nrows, consdata->vars1, vals1) );
         SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nrows, consdata->vars2, vals2) );

         SCIPdebugMsg(scip, "Separating orbisack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         SCIP_CALL( separateInequalities(scip, result, conss[c], consdata->nrows, consdata->vars1, consdata->vars2, vals1, vals2) );

         SCIPfreeBufferArray(scip, &vals2);
         SCIPfreeBufferArray(scip, &vals1);

         if ( *result == SCIP_CUTOFF )
            break;
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution */
static
SCIP_DECL_CONSSEPASOL(consSepasolOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for orbisack constraints\n");

   *result = SCIP_DIDNOTFIND;

   if ( nconss > 0 )
   {
      int nvals;

      nvals = SCIPgetNVars(scip)/2;
      SCIP_CALL( SCIPallocBufferArray(scip, &vals1, nvals) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals2, nvals) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);

         /* get solution */
         assert( consdata->nrows <= nvals );
         SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nrows, consdata->vars1, vals1) );
         SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nrows, consdata->vars2, vals2) );

         SCIPdebugMsg(scip, "Separating orbisack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         SCIP_CALL( separateInequalities(scip, result, conss[c], consdata->nrows, consdata->vars1, consdata->vars2, vals1, vals2) );
         if ( *result == SCIP_CUTOFF )
            break;
      }

      SCIPfreeBufferArray(scip, &vals2);
      SCIPfreeBufferArray(scip, &vals1);
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions
 *
 *  @pre It is assumed that the solution is integral (this can be ensured by appropriate priorities).
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   int ngen = 0;
   int c;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != 0 );

   SCIPdebugMsg(scip, "Enfolp method for orbisack constraints\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   if ( nconss > 0 )
   {
      int nvals;

      nvals = SCIPgetNVars(scip)/2;
      SCIP_CALL( SCIPallocBufferArray(scip, &vals1, nvals) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals2, nvals) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         /* get data of constraint */
         assert( conss[c] != 0 );
         consdata = SCIPconsGetData(conss[c]);

         /* get solution */
         assert( consdata->nrows <= nvals );
         SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nrows, consdata->vars1, vals1) );
         SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nrows, consdata->vars2, vals2) );

         SCIPdebugMsg(scip, "Enforcing orbisack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* Separate only cover inequalities to ensure that enforcing works correctly. */
         /* Otherwise, it may happen that infeasible solutions cannot be detected, since */
         /* we bound the size of the coefficients for the orbisack inequalities. */
         SCIP_CALL( separateOrbisackCovers(scip, conss[c], consdata->nrows, consdata->vars1, consdata->vars2, vals1, vals2, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            break;
         }

         SCIPdebugMsg(scip, "Generated orbisack inequalities for <%s>: %d\n", SCIPconsGetName(conss[c]), ngen);

         if ( ngen > 0 )
            *result = SCIP_SEPARATED;
      }

      SCIPfreeBufferArray(scip, &vals2);
      SCIPfreeBufferArray(scip, &vals1);
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsOrbisack)
{  /*lint --e{715}*/
   SCIP_Bool feasible = TRUE;
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for orbisack constraints (pseudo solutions) ...\n");

   *result = SCIP_FEASIBLE;

   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL);
      assert( consdata->nrows > 0 );
      assert( consdata->vars1 != NULL );
      assert( consdata->vars2 != NULL );

      SCIP_CALL( SCIPcheckSolutionOrbisack(scip, NULL, consdata->vars1, consdata->vars2, consdata->nrows, FALSE, &feasible) );

      if ( ! feasible )
      {
         *result = SCIP_INFEASIBLE;
         break;
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   int ngen = 0;
   int c;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != 0 );

   SCIPdebugMsg(scip, "Enforelax method for orbisack constraints.\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   if ( nconss > 0 )
   {
      int nvals;

      nvals = SCIPgetNVars(scip)/2;
      SCIP_CALL( SCIPallocBufferArray(scip, &vals1, nvals) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals2, nvals) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         /* get data of constraint */
         assert( conss[c] != 0 );
         consdata = SCIPconsGetData(conss[c]);

         /* get solution */
         assert( consdata->nrows <= nvals );
         SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nrows, consdata->vars1, vals1) );
         SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nrows, consdata->vars2, vals2) );

         SCIPdebugMsg(scip, "Enforcing orbisack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* Separate only cover inequalities to ensure that enforcing works correctly. */
         /* Otherwise, it may happen that infeasible solutions cannot be detected, since */
         /* we bound the size of the coefficients for the orbisack inequalities. */
         SCIP_CALL( separateOrbisackCovers(scip, conss[c], consdata->nrows, consdata->vars1, consdata->vars2, vals1, vals2, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            break;
         }

         SCIPdebugMsg(scip, "Generated orbisack inequalities for <%s>: %d\n", SCIPconsGetName(conss[c]), ngen);

         if ( ngen > 0 )
            *result = SCIP_SEPARATED;
      }

      SCIPfreeBufferArray(scip, &vals2);
      SCIPfreeBufferArray(scip, &vals1);
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckOrbisack)
{  /*lint --e{715}*/
   SCIP_Bool feasible = TRUE;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->checkalwaysfeas )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL);
      assert( consdata->nrows > 0 );
      assert( consdata->vars1 != NULL );
      assert( consdata->vars2 != NULL );

      SCIPdebugMsg(scip, "Check method for orbisack constraint <%s> (%d rows) ...\n", SCIPconsGetName(conss[c]), consdata->nrows);

      SCIP_CALL( SCIPcheckSolutionOrbisack(scip, sol, consdata->vars1, consdata->vars2, consdata->nrows, printreason, &feasible) );

      if ( ! feasible )
      {
         *result = SCIP_INFEASIBLE;
         SCIPdebugMsg(scip, "Solution is feasible.\n");
         break;
      }
   }

   if ( feasible )
      SCIPdebugMsg(scip, "Solution is feasible.\n");

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropOrbisack)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Propagation method of orbisack constraint handler.\n");

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool found = FALSE;
      int ngen = 0;

      assert( conss[c] != NULL );

      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &found, &ngen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if ( found )
         *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolOrbisack)
{  /*lint --e{715}*/
   int c;
   int ngen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Presolving method of orbisack constraint handler. Propagating orbisack inequalities.\n");

   *result = SCIP_DIDNOTFIND;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_Bool found = FALSE;
      int curngen = 0;

      assert( conss[c] != NULL );
      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &found, &curngen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      ngen += curngen;
   }

   if ( ngen > 0 )
   {
      *nfixedvars += ngen;
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis */
static
SCIP_DECL_CONSRESPROP(consRespropOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Propagation resolution method of orbisack constraint handler.\n");

   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->nrows > 0 );
   assert( consdata->vars1 != NULL );
   assert( consdata->vars2 != NULL );

   vars1 = consdata->vars1;
   vars2 = consdata->vars2;

   assert( 0 <= inferinfo && inferinfo < consdata->nrows );

   assert( vars2[inferinfo] == infervar );
   assert( SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, TRUE) < 0.5 );

   if ( SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, FALSE) > 0.5 && SCIPvarGetUbAtIndex(vars2[inferinfo], bdchgidx, TRUE) < 0.5 )
   {
      SCIPdebugMsg(scip, " -> reason for setting x[%d][1] = 0 was fixing x[%d][0] to 0 and each row above is fixed to the same value.\n",
         inferinfo, inferinfo);

      for (i = 0; i < inferinfo; ++i)
      {
         SCIP_CALL( SCIPaddConflictUb(scip, vars1[i], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars1[i], bdchgidx) );
         SCIP_CALL( SCIPaddConflictUb(scip, vars2[i], bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, vars2[i], bdchgidx) );
      }
      SCIP_CALL( SCIPaddConflictUb(scip, vars1[inferinfo], bdchgidx) );

      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/** Lock variables
 *
 *  We assume we have only one global (void) constraint and lock all variables.
 *
 * - Orbisack constraints may get violated if the variables of the first column
 *   are rounded down, we therefor call SCIPaddVarLocks(..., nlockspos, nlocksneg).
 * - Orbisack constraints may get violated if the variables of the second column
 *   are rounded up , we therefor call SCIPaddVarLocks(..., nlocksneg, nlockspo ).
 */
static
SCIP_DECL_CONSLOCK(consLockOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nrows;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "Locking method for orbisack constraint handler.\n");

   /* get data of original constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->nrows > 0 );
   assert( consdata->vars1 != NULL );
   assert( consdata->vars2 != NULL );

   nrows = consdata->nrows;
   vars1 = consdata->vars1;
   vars2 = consdata->vars2;

   for (i = 0; i < nrows; ++i)
   {
      SCIP_CALL( SCIPaddVarLocks(scip, vars1[i], nlockspos, nlocksneg) );
      SCIP_CALL( SCIPaddVarLocks(scip, vars2[i], nlocksneg, nlockspos) );
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler
 *
 *  The constraint handler should output a representation of the constraint into the given text file.
 */
static
SCIP_DECL_CONSPRINT(consPrintOrbisack)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nrows;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars1 != NULL );
   assert( consdata->vars2 != NULL );
   assert( consdata->nrows > 0 );

   vars1 = consdata->vars1;
   vars2 = consdata->vars2;
   nrows = consdata->nrows;

   SCIPdebugMsg(scip, "Printing method for orbisack constraint handler\n");

   SCIPinfoMessage(scip, file, "orbisack(");

   for (i = 0; i < nrows; ++i)
   {
      SCIPinfoMessage(scip, file, "%s,%s", SCIPvarGetName(vars1[i]), SCIPvarGetName(vars2[i]));
      if ( i < nrows-1 )
         SCIPinfoMessage(scip, file, ".");
   }

   return SCIP_OKAY;
}


/** checks given solution for feasibility */
SCIP_RETCODE SCIPcheckSolutionOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to check for feasibility */
   SCIP_VAR**            vars1,              /**< variables of first column */
   SCIP_VAR**            vars2,              /**< variables of second column */
   int                   nrows,              /**< number of rows */
   SCIP_Bool             printreason,        /**< whether reason for infeasibility should be printed */
   SCIP_Bool*            feasible            /**< memory address to store whether sol is feasible */
   )
{
   int i;
   int val1;
   int val2;

   assert( scip != NULL );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( nrows > 0 );
   assert( feasible != NULL );

   *feasible = TRUE;

   /* find first non-constant row and check for feasibility */
   for (i = 0; i < nrows; ++i)
   {
      assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars1[i])) );
      assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars2[i])) );

      /* get values of i-th row */
      val1 = SCIPgetSolVal(scip, sol, vars1[i]) > 0.5 ? 1 : 0;
      val2 = SCIPgetSolVal(scip, sol, vars2[i]) > 0.5 ? 1 : 0;

      /* if row i is constrant */
      if ( val1 == val2 )
         continue;
      /* row i has type (1,0) -> feasible */
      else if ( val1 == 1 )
      {
         assert( val2 == 0 );
         break;
      }
      else /* infeasible */
      {
         if ( printreason )
            SCIPinfoMessage(scip, NULL, "First non-constant row %d is fixed to (0,1).\n", i);
         *feasible = FALSE;
         break;
      }
   }

   return SCIP_OKAY;
}


/** separate orbisack cover inequalities for a given solution */
SCIP_RETCODE SCIPseparateCoversOrbisack(
   SCIP*                 scip,               /**< pointer to scip */
   SCIP_CONS*            cons,               /**< pointer to constraint for which cover inequality should be added */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_VAR**            vars1,              /**< variables of first columns */
   SCIP_VAR**            vars2,              /**< variables of second columns */
   int                   nrows,              /**< number of rows */
   SCIP_Bool*            infeasible,         /**< memory address to store whether we detected infeasibility */
   int*                  ngen                /**< memory address to store number of generated cuts */
   )
{
   SCIP_ROW* row;
   SCIP_Real rhs = 0.0;
   SCIP_Real lhs = 0.0;
   SCIP_Real* solvals1;
   SCIP_Real* solvals2;
   SCIP_Real* coeff1;
   SCIP_Real* coeff2;
   int i;
   int j;

   assert( scip != NULL );
   assert( vars1 != NULL );
   assert( vars2 != NULL );
   assert( nrows > 0 );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   *infeasible = FALSE;
   *ngen = 0;

   /* allocate memory for inequality coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff2, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals1, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals2, nrows) );

   /* initialize coefficient matrix and solution values */
   for (i = 0; i < nrows; ++i)
   {
      coeff1[i] = 0.0;
      coeff2[i] = 0.0;
   }
   SCIP_CALL( SCIPgetSolVals(scip, sol, nrows, vars1, solvals1) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nrows, vars2, solvals2) );

   /* detect violated covers */
   for (i = 0; i < nrows; ++i)
   {
      /* cover inequality is violated */
      if ( SCIPisEfficacious(scip, -solvals1[i] + solvals2[i] + lhs - rhs) )
      {
         /* set coefficients for inequality */
         coeff1[i] = -1.0;
         coeff2[i] = 1.0;

         /* add inequality */
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "orbisackcover", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
         for (j = 0; j < nrows; ++j)
         {
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars1[j], coeff1[j]) );
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars2[j], coeff2[j]) );
         }
         SCIP_CALL( SCIPflushRowExtensions(scip, row) );

         SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         ++(*ngen);
         if ( *infeasible )
            break;

         /* reset coefficients for next inequality */
         coeff1[i] = 0.0;
         coeff2[i] = 0.0;
      }

      /* add argmax( 1 - solvals[i][0], solvals[i][1] ) as coefficient */
      if ( SCIPisEfficacious(scip, 1.0 - solvals1[i] - solvals2[i]) )
      {
         coeff1[i] = -1.0;
         lhs = lhs - solvals1[i];
      }
      else
      {
         coeff2[i] = 1.0;
         rhs += 1.0;
         lhs = lhs + solvals2[i];
      }
   }

   /* free coefficient matrix */
   SCIPfreeBufferArray(scip, &solvals2);
   SCIPfreeBufferArray(scip, &solvals1);
   SCIPfreeBufferArray(scip, &coeff2);
   SCIPfreeBufferArray(scip, &coeff1);

   return SCIP_OKAY;

}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( vars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( varssize < 2 * consdata->nrows )
      (*success) = FALSE;
   else
   {
      int cnt = 0;
      int i;

      for (i = 0; i < consdata->nrows; ++i)
      {
         vars[cnt++] = consdata->vars1[i];
         vars[cnt++] = consdata->vars2[i];
      }
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsOrbisack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   (*nvars) = 2 * consdata->nrows;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/** creates the handler for orbisack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrOrbisack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpOrbisack, consEnfopsOrbisack, consCheckOrbisack, consLockOrbisack,
         conshdlrdata) );
   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxOrbisack) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeOrbisack) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteOrbisack) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsOrbisack) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsOrbisack) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolOrbisack, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintOrbisack) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropOrbisack, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropOrbisack) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpOrbisack, consSepasolOrbisack, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransOrbisack) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpOrbisack) );

   /* separation methods */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/coverseparation",
         "Separate cover inequalities for orbisacks?",
         &conshdlrdata->coverseparation, TRUE, DEFAULT_COVERSEPARATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/orbiSeparation",
         "Separate orbisack inequalities?",
         &conshdlrdata->orbiseparation, TRUE, DEFAULT_ORBISEPARATION, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/coeffbound",
         "Maximum size of coefficients for orbisack inequalities",
         &conshdlrdata->coeffbound, TRUE, DEFAULT_COEFFBOUND, 0.0, DBL_MAX, NULL, NULL) );

   /* whether we allow upgrading to packing/partioning orbisack constraints*/
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkpporbisack",
         "Upgrade orbisack constraints to packing/partioning orbisacks?",
         &conshdlrdata->checkpporbisack, TRUE, DEFAULT_PPORBISACK, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkalwaysfeas",
         "Whether check routine returns always SCIP_FEASIBLE.",
         &conshdlrdata->checkalwaysfeas, TRUE, DEFAULT_CHECKALWAYSFEAS, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates and captures a orbisack constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*const*       vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR*const*       vars2,              /**< second column of matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of rows in variable matrix */
   SCIP_Bool             ispporbisack,       /**< whether the orbisack is a packing/partitioning orbisack */
   SCIP_Bool             isparttype,         /**< whether the orbisack is a partitioning orbisack */
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   SCIP_Bool success;
   SCIP_ORBITOPETYPE orbitopetype;
   int i;

   /* find the orbisack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("orbisack constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert( nrows > 0 );

   /* check for upgrade to packing/partitioning orbisacks*/
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   if ( ! ispporbisack && conshdlrdata->checkpporbisack )
   {
      SCIP_CALL( packingUpgrade(scip, vars1, vars2, nrows, &success, &isparttype) );

      if ( success )
         ispporbisack = TRUE;
   }

   /* create constraint, if it is a packing/partitioning orbisack, add orbitope constraint
    * instead of orbitsack constraint */
   if (  ispporbisack )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nrows) );
      for (i = 0; i < nrows; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &vars[i], 2) ); /*lint !e866*/
         vars[i][0] = vars1[i];
         vars[i][1] = vars2[i];
      }

      if ( isparttype )
         orbitopetype = SCIP_ORBITOPETYPE_PARTITIONING;
      else
         orbitopetype = SCIP_ORBITOPETYPE_PACKING;

      SCIP_CALL( SCIPcreateConsOrbitope(scip, cons, "pporbisack", vars, orbitopetype, nrows, 2, TRUE, initial, separate, enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );

      for (i = 0; i < nrows; ++i)
         SCIPfreeBufferArray(scip, &vars[i]);
      SCIPfreeBufferArray(scip, &vars);
   }
   else
   {
      /* create constraint data */
      SCIP_CALL( consdataCreate(scip, &consdata, vars1, vars2, nrows) );

      SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
            local, modifiable, dynamic, removable, stickingatnode) );
   }

   return SCIP_OKAY;
}


/** creates and captures an orbisack constraint in its most basic variant
 *
 *  All constraint flags set to their default values, which can be set afterwards using SCIPsetConsFLAGNAME() in scip.h.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicOrbisack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR**            vars1,              /**< first column of matrix of variables on which the symmetry acts */
   SCIP_VAR**            vars2,              /**< second column of matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of rows in constraint matrix */
   SCIP_Bool             ispporbisack,       /**< whether the orbisack is a packing/partitioning orbisack */
   SCIP_Bool             isparttype          /**< whether the orbisack is a partitioning orbisack */
   )
{
   SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, name, vars1, vars2, nrows, ispporbisack, isparttype,
         TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
