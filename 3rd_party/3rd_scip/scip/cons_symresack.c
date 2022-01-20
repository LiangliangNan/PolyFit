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

/**@file   cons_symresack.c
 * @brief  constraint handler for symresack constraints
 * @author Christopher Hojny
 *
 * The type of constraints of this constraint handler is described in cons_symresack.h.
 *
 * The details of the method implemented here are described in the following papers:
 *
 * Fundamental Domains for Integer Programs with Symmetries@n
 * Eric J. Friedman,@n
 * Combinatorial Optimization, volume 4616 of LNCS, 146-153 (2007)
 *
 * This paper describes an inequality to handle symmetries of a single permutation. This
 * so-called FD-inequality is the basic for the propagation routine of our implementation.
 *
 * Polytopes Associated with Symmetry Handling@n
 * Christopher Hojny and Marc E. Pfetsch,@n
 * (2017), preprint available at http://www.optimization-online.org/DB_HTML/2017/01/5835.html
 *
 * This paper describes an almost linear time separation routine for so-called cove
 * inequalities of symresacks. In our implementation, however, we use a separation routine with
 * quadratic worst case running time.
 *
 * Packing, Partitioning, and Covering Symresacks@n
 * Christopher Hojny,@n
 * (2017), preprint available at http://www.optimization-online.org/DB_HTML/2017/05/5990.html
 *
 * This paper introduces linearly many inequalities with ternary coefficients that suffice to
 * characterize the binary points contained in a packing and partitioning symresack completely.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_orbisack.h"
#include "scip/cons_setppc.h"
#include "scip/cons_symresack.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "symresack"
#define CONSHDLR_DESC          "symmetry breaking constraint handler relying on symresacks"
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

#define DEFAULT_PPSYMRESACK       FALSE /**< whether we allow upgrading to packing/partitioning symresacks */
#define DEFAULT_CHECKALWAYSFEAS    TRUE /**< whether check routine returns always SCIP_FEASIBLE */

/* macros for getting bounds of pseudo solutions in propagation */
#define ISFIXED0(x)   (SCIPvarGetUbLocal(x) < 0.5 ? TRUE : FALSE)
#define ISFIXED1(x)   (SCIPvarGetLbLocal(x) > 0.5 ? TRUE : FALSE)


/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             checkppsymresack;   /**< whether we allow upgrading to packing/partitioning symresacks */
   SCIP_Bool             checkalwaysfeas;    /**< whether check routine returns always SCIP_FEASIBLE */
};


/** constraint data for symresack constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables */
   int                   nvars;              /**< number of variables */
   int*                  perm;               /**< permutation associated to the symresack */
   int*                  invperm;            /**< inverse permutation */
   SCIP_Bool             ppupgrade;          /**< whether constraint is upgraded to packing/partitioning symresack */
#ifdef SCIP_DEBUG
   int                   debugcnt;           /**< counter to store number of added cover inequalities */
#endif

   /* data for upgraded symresack constraints */
   int                   ncycles;            /**< number of cycles in permutation */
   int**                 cycledecomposition; /**< cycle decomposition */
};


/*
 * Local methods
 */

/** frees a symresack constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to symresack constraint data */
   )
{
   int nvars;
   int i;

   assert( consdata != NULL );
   assert( *consdata != NULL );

   nvars = (*consdata)->nvars;

   if ( nvars == 0 )
   {
      SCIPfreeBlockMemory(scip, consdata);

      return SCIP_OKAY;
   }

   if ( (*consdata)->ppupgrade )
   {
      for (i = 0; i < (*consdata)->ncycles; ++i)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycledecomposition[i]), nvars + 1);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->cycledecomposition), (*consdata)->ncycles);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->invperm), nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->perm), nvars);

   for (i = 0; i < nvars; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->vars), nvars);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** check whether constraint can be upgraded to packing/partitioning symresack */
static
SCIP_RETCODE packingUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables affected by permutation */
   int                   nvars,              /**< length of permutation */
   SCIP_Bool*            upgrade             /**< pointer to store whether upgrade was successful */
   )
{
   SCIP_Bool* covered;
   SCIP_Bool descent;
   SCIP_CONSHDLR* setppcconshdlr;
   SCIP_CONS** setppcconss;
   SCIP_VAR* var;
   SCIP_Bool terminated = FALSE;
   int** cycledecomposition;
   int* indicesincycle;
   int nsetppcconss;
   int curcycle;
   int maxcyclelength;
   int ncycles = 0;
   int c;
   int i;
   int j;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );
   assert( upgrade != NULL );

   *upgrade = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &covered, nvars) );

   for (i = 0; i < nvars; ++i)
      covered[i] = FALSE;

   /* check wether permutation is monotone */
   for (i = 0; i < nvars; ++i)
   {
      /* skip checked indices */
      if ( covered[i] )
         continue;

      ++ncycles;
      j = i;
      descent = FALSE;

      do
      {
         covered[j] = TRUE;

         if ( perm[j] < j )
         {
            if ( ! descent )
               descent = TRUE;
            else
               break;
         }

         j = perm[j];
      }
      while ( j != i );

      /* if cycle is not monotone */
      if ( j != i )
      {
         SCIPfreeBufferArray(scip, &covered);

         return SCIP_OKAY;
      }
   }
   assert( ncycles <= nvars / 2 );

   /* each cycle is monotone; check for packing/partitioning type */
   for (i = 0; i < nvars; ++i)
      covered[i] = FALSE;

   /* compute cycle decomposition: row i stores in entry 0 the length of the cycle,
    * the remaining entries are the coordinates in the cycle */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycledecomposition, ncycles) );
   for (i = 0; i < ncycles; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cycledecomposition[i], nvars + 1) );
   }

   curcycle = 0;
   maxcyclelength = 0;
   for (i = 0; i < nvars; ++i)
   {
      int cyclelength = 0;

      /* skip checked indices */
      if ( covered[i] )
         continue;

      j = i;
      do
      {
         covered[j] = TRUE;
         cycledecomposition[curcycle][++cyclelength] = j;
         j = perm[j];
      }
      while ( j != i );

      cycledecomposition[curcycle][0] = cyclelength;
      ++curcycle;

      if ( maxcyclelength < cyclelength )
         maxcyclelength = cyclelength;
   }

   /* permutation can be upgraded -> check whether the symresack is of packing/partitioning type */
   setppcconshdlr = SCIPfindConshdlr(scip, "setppc");
   if ( setppcconshdlr == NULL )
   {
      SCIPerrorMessage("Setppc constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }
   setppcconss = SCIPconshdlrGetConss(setppcconshdlr);
   nsetppcconss = SCIPconshdlrGetNConss(setppcconshdlr);

   /* Check whether each cycle of the symresack is contained in a set packing/partitioning constraint.
    * To this end, we have to guarantee that all affected variables are not negated since permutations
    * are given w.r.t. original variables. */
   *upgrade = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &indicesincycle, maxcyclelength) );

   for (i = 0; i < ncycles && *upgrade && ! terminated; ++i)
   {
      int cyclelength;

      /* get indices of variables in current cycle */
      for (j = 0; j < cycledecomposition[i][0]; ++ j)
      {
         var = vars[cycledecomposition[i][j + 1]];

         if ( SCIPvarIsNegated(var) )
         {
            terminated = TRUE;
            break;
         }

         indicesincycle[j] = SCIPvarGetProbindex(var);
      }

      cyclelength = cycledecomposition[i][0];

      /* iterate over constraints
       *
       * @todo Improve the check by sorting the constraints in the setppcconss array
       * by type and number of contained variables. */
      for (c = 0; c < nsetppcconss; ++c)
      {
         int nsetppcvars;
         SCIP_VAR** setppcvars;
         int varidx;
         int nfound = 0;
         int k;

         /* check type */
         if ( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_COVERING )
            continue;
         assert( SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PARTITIONING || SCIPgetTypeSetppc(scip, setppcconss[c]) == SCIP_SETPPCTYPE_PACKING );

         /* get set packing/partitioning variables */
         nsetppcvars = SCIPgetNVarsSetppc(scip, setppcconss[c]);
         assert( nsetppcvars > 0 );

         setppcvars = SCIPgetVarsSetppc(scip, setppcconss[c]);
         assert( setppcvars != NULL );

         /* check whether all variables of the cycle are contained in setppc constraint */
         for (j = 0; j < nsetppcvars && nfound < cyclelength; ++j)
         {
            var = setppcvars[j];

            if ( SCIPvarIsNegated(var) )
               continue;

            varidx = SCIPvarGetProbindex(var);

            for (k = 0; k < cyclelength; ++k)
            {
               if ( varidx == indicesincycle[k] )
               {
                  ++nfound;
                  break;
               }
            }
         }
         assert( nfound <= cyclelength );

         if ( nfound == cyclelength )
            break;
      }

      /* row is not contained in a set packing/partitioning constraint */
      if ( c >= nsetppcconss )
         *upgrade = FALSE;
   }

   if ( *upgrade )
   {
      (*consdata)->ncycles = ncycles;
      (*consdata)->cycledecomposition = cycledecomposition;

      SCIPfreeBufferArray(scip, &indicesincycle);
      SCIPfreeBufferArray(scip, &covered);
   }
   else
   {
      SCIPfreeBufferArray(scip, &indicesincycle);
      SCIPfreeBufferArray(scip, &covered);
      for (i = 0; i < ncycles; ++i)
      {
         SCIPfreeBlockMemoryArray(scip, &cycledecomposition[i], nvars + 1);
      }
      SCIPfreeBlockMemoryArray(scip, &cycledecomposition, ncycles);
   }

   return SCIP_OKAY;
}


/** creates symresack constraint data
 *
 *  If the input data contain non-binary variables of fixed
 *  points, we delete these variables in a preprocessing step.
 */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store constraint data */
   SCIP_VAR*const*       inputvars,          /**< input variables of the constraint handler */
   int                   inputnvars,         /**< input number of variables of the constraint handler*/
   int*                  inputperm           /**< input permutation of the constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** vars;
   SCIP_Bool upgrade;
   int* indexcorrection;
   int* invperm;
   int* perm;
   int naffectedvariables;
   int i;
   int j = 0;

   assert( consdata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

#ifdef SCI_DEBUG
   consdata->debugcnt = 0;
#endif

   /* count the number of binary variables which are affected by the permutation */
   SCIP_CALL( SCIPallocBufferArray(scip, &indexcorrection, inputnvars) );
   indexcorrection[0] = -1;
   for (i = 0; i < inputnvars; ++i)
   {
      if ( inputperm[i] != i && SCIPvarIsBinary(inputvars[i]) )
      {
         if ( i == 0 )
            indexcorrection[i] = 0;
         else
            indexcorrection[i] = indexcorrection[i - 1] + 1;
      }
      else
      {
         if ( i > 0 )
            indexcorrection[i] = indexcorrection[i - 1];
      }
   }
   naffectedvariables = indexcorrection[inputnvars - 1] + 1;

   (*consdata)->nvars = naffectedvariables;

   /* Stop if we detect that the permutation fixes each binary point. */
   if ( naffectedvariables == 0 )
   {
      SCIPfreeBufferArrayNull(scip, &indexcorrection);
      return SCIP_OKAY;
   }

   /* remove fixed points from permutation representation */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, naffectedvariables) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &perm, naffectedvariables) );
   for (i = 0; i < inputnvars; ++i)
   {
      if ( i == 0 )
      {
         if ( indexcorrection[i] > -1 )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
      else
      {
         if ( indexcorrection[i] > indexcorrection[i - 1] )
         {
            vars[j] = inputvars[i];
            perm[j++] = indexcorrection[inputperm[i]];
         }
      }
   }
   SCIPfreeBufferArrayNull(scip, &indexcorrection);

   (*consdata)->vars = vars;
   (*consdata)->perm = perm;

   for (i = 0; i < naffectedvariables; ++i)
   {
      SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vars[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &invperm, naffectedvariables) );
   for (i = 0; i < naffectedvariables; ++i)
      invperm[perm[i]] = i;
   (*consdata)->invperm = invperm;

   /* check for upgrade to packing/partitioning symresacks*/
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("symresack constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   upgrade = FALSE;
   if ( conshdlrdata->checkppsymresack )
   {
      SCIP_CALL( packingUpgrade(scip, consdata, perm, vars, naffectedvariables, &upgrade) );
   }

   (*consdata)->ppupgrade = upgrade;

   /* get transformed variables, if we are in the transformed problem */
   if ( SCIPisTransformed(scip) )
   {
      /* Make sure that all variables cannot be multiaggregated (cannot be handled by cons_symresack, since one cannot
       * easily eliminate single variables from a symresack constraint. */
      for (i = 0; i < naffectedvariables; ++i)
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vars[i], &(*consdata)->vars[i]) );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, (*consdata)->vars[i]) );
      }
   }

   return SCIP_OKAY;
}


/** generate initial LP cut
 *
 *  We generate the ordering inequality for the pair \f$(1, \gamma^{-1}(1))\f$, i.e.,
 *  the inequality \f$-x_{1} + x_{\gamma^{-1}(1)} \leq 0\f$. This inequality is valid,
 *  because we guaranteed in a preprocessing step that all variables are binary.
 *
 *  Furthermore, we add facet inequalities of packing/partitioning symresacks if
 *  we deal with packing/partitioning symresacks.
 */
static
SCIP_RETCODE initLP(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   int nvars;
#ifdef SCIP_DEBUG
   char name[SCIP_MAXSTRLEN];
#endif
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;

   /* avoid stupid problems */
   if ( nvars <= 1 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   vars = consdata->vars;

   /* there are no fixed points */
   assert( consdata->invperm[0] != 0 );

   /* add ordering inequality */
#ifdef SCIP_DEBUG
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_init_%s", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#else
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#endif
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[0], -1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[consdata->invperm[0]], 1.0) );

   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   /* check whether we have a packing/partioning symresack */
   if ( consdata->ppupgrade && ! *infeasible )
   {
      SCIP_VAR** varsincons;
      SCIP_Real* coeffs;
      int** cycledecomposition;
      int ncycles;
      int nvarsincons;
      int nvarsincycle;
      int firstelemincycle;

      ncycles = consdata->ncycles;
      cycledecomposition = consdata->cycledecomposition;

      SCIP_CALL( SCIPallocBufferArray(scip, &varsincons, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nvars) );

      coeffs[0] = 1.0;

      /* add packing/partitioning symresack constraints */
      for (i = 0; i < ncycles; ++i)
      {
         assert( cycledecomposition[i][0] > 0 );

         nvarsincycle = cycledecomposition[i][0];
         varsincons[0] = vars[cycledecomposition[i][nvarsincycle]];
         firstelemincycle = cycledecomposition[i][1];

         assert( firstelemincycle == consdata->perm[cycledecomposition[i][nvarsincycle]] );

         nvarsincons = 1;

         /* add variables of other cycles to the constraint */
         for (j = 0; j < i; ++j)
         {
            nvarsincycle = cycledecomposition[j][0];
            for (k = 1; k <= nvarsincycle; ++k)
            {
               if ( cycledecomposition[j][k] < firstelemincycle )
               {
                  varsincons[nvarsincons] = vars[cycledecomposition[j][k]];
                  coeffs[nvarsincons++] = -1.0;
               }
               else
                  continue;
            }
         }

#ifdef SCIP_DEBUG
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ppSymresack_%d_%s", i, SCIPconsGetName(cons));
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#else
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "", -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
#endif
         SCIP_CALL( SCIPaddVarsToRow(scip, row, nvarsincons, varsincons, coeffs) );

         SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         if ( *infeasible )
            break;
      }

      SCIPfreeBufferArray(scip, &coeffs);
      SCIPfreeBufferArray(scip, &varsincons);
   }

   return SCIP_OKAY;
}


/** perform propagation of symresack constraint */
static
SCIP_RETCODE propVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to be propagated */
   SCIP_Bool*            infeasible,         /**< pointer to store whether it was detected that the node is infeasible */
   int*                  ngen                /**< pointer to store number of generated bound strengthenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* invperm;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   SCIPdebugMsg(scip, "Propagating variables of constraint <%s>.\n", SCIPconsGetName(cons));

   *ngen = 0;
   *infeasible = FALSE;

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->nvars > 0 );
   nvars = consdata->nvars;

   /* avoid trivial problems */
   if ( nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->invperm != NULL );
   vars = consdata->vars;
   invperm = consdata->invperm;

   /* loop through all variables */
   for (i = 0; i < nvars; ++i)
   {
      SCIP_VAR* var2;
      SCIP_VAR* var;
      int r;
      SCIP_Bool tightened;

      /* there are no fixed points */
      assert( invperm[i] != i );

      /* get variables of first and second column */
      var = vars[i];
      var2 = vars[invperm[i]];
      assert( var != NULL );
      assert( var2 != NULL );

      /* if first part of variable pair fixed to 0 and second part is fixed to 1 */
      if ( ISFIXED0(var) && ISFIXED1(var2) )
      {
         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);

         SCIPdebugMsg(scip, " -> node infeasible (pair was fixed to (0,1) but there was no pair of type (1,0) before) ---> lexicographical order violated, infeasible.\n");

         /* perform conflict analysis */
         if ( SCIPisConflictAnalysisApplicable(scip) )
         {
            SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

            for (r = 0; r <= i; ++r)
            {
               /* there are no fixed points */
               assert( invperm[r] != r );

               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[r]) );
               SCIP_CALL( SCIPaddConflictBinvar(scip, vars[invperm[r]]) );
            }

            SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );
         }

         *infeasible = TRUE;
         break;
      }
      /* if first part of the variable pair is fixed to 0 and the second part is free --> fix second part to 0 */
      else if ( ISFIXED0(var) && ( ! ISFIXED0(var2) ) )
      {
         assert( SCIPvarGetUbLocal(var) < 0.5 );
         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         assert( SCIPvarGetUbLocal(var2) > 0.5 );

         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);

         assert( SCIPvarGetLbLocal(var2) < 0.5 );
         SCIP_CALL( SCIPinferVarUbCons(scip, var2, 0.0, cons, i, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         if ( tightened )
            ++(*ngen);
      }
      /* if second part of the variable pair is fixed to 1 and the first part is free --> fix first part to 1 */
      else if ( ( ! ISFIXED1(var) ) && ISFIXED1(var2) )
      {
         assert( SCIPvarGetLbLocal(var) < 0.5 );
         assert( SCIPvarGetUbLocal(var) > 0.5 );
         assert( SCIPvarGetLbLocal(var2) > 0.5 );

         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);

         assert( SCIPvarGetUbLocal(var) > 0.5 );
         SCIP_CALL( SCIPinferVarLbCons(scip, var, 1.0, cons, i + nvars, FALSE, infeasible, &tightened) ); /*lint !e713*/
         assert( ! *infeasible );

         if ( tightened )
            ++(*ngen);
      }
      /* if solution is lexicographically maximal */
      else if ( ISFIXED1(var) && ISFIXED0(var2) )
      {
         assert( SCIPvarGetLbLocal(var) > 0.5 );
         assert( SCIPvarGetUbLocal(var2) < 0.5 );

         SCIPdebugMsg(scip, "Check variable pair (%d,%d).\n", i, invperm[i]);
         SCIPdebugMsg(scip, " -> node is feasible (pair was fixed to (1,0) and every earlier pair is constant).\n");

         break;
      }
      /* cannot apply propagation */
      else
         break;
   }

   return SCIP_OKAY;
}


/** add symresack cover inequality */
static
SCIP_RETCODE addSymresackInequality(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< variables */
   int*                  coeffs,             /**< coefficient vector of inequality to be added */
   SCIP_Real             rhs,                /**< right-hand side of inequality to be added */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_ROW* row;
   int i;
#ifdef SCIP_DEBUG
   SCIP_CONSDATA* consdata;
   char name[SCIP_MAXSTRLEN];
#endif

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars > 0 );
   assert( vars != NULL );
   assert( coeffs != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

#ifdef SCIP_DEBUG
   consdata = SCIPconsGetData(cons);
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_cover_%s_%d", SCIPconsGetName(cons), consdata->debugcnt);
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), name, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
   consdata->debugcnt += 1;
#else
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, SCIPconsGetHdlr(cons), "", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
#endif
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for (i = 0; i < nvars; ++i)
   {
      if ( coeffs[i] == 1 || coeffs[i] == -1 )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i], (SCIP_Real) coeffs[i]) );
      }
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );
   SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}


/** separate symresack cover inequalities
 *
 *  We currently do NOT enter cuts into the pool.
 */
static
SCIP_RETCODE separateSymresackCovers(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   const SCIP_CONSDATA*  consdata,           /**< constraint data */
   SCIP_Real*            vals,               /**< solution values of variables */
   int*                  ngen,               /**< pointer to store the number of separated covers */
   SCIP_Bool*            infeasible          /**< pointer to store whether we detected infeasibility */
   )
{
   SCIP_Real constobjective;
   SCIP_Real* sepaobjective;
   SCIP_Real tmpsoluobj = 0.0;
   SCIP_Real maxsoluobj = 0.0;
   int* tmpsolu;
   int* maxsolu;
   int* invperm;
   int* perm;
   int nvars;
   int crit;
   int i;

   *infeasible = FALSE;
   *ngen = 0;

   assert( scip != NULL );
   assert( consdata != NULL );

   /* we don't have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );
   assert( consdata->invperm != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   nvars = consdata->nvars;
   perm = consdata->perm;
   invperm = consdata->invperm;

   /* initialize objective */
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaobjective, nvars) );

   constobjective = 1.0; /* constant part of separation objective */
   for (i = 0; i < nvars; ++i)
   {
      if ( i < perm[i] )
      {
         sepaobjective[i] = vals[i];
         constobjective -= vals[i];
      }
      else
         sepaobjective[i] = vals[i] - 1.0;
   }

   /* allocate memory for temporary and global solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpsolu, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxsolu, nvars) );

   /* start separation procedure by iterating over critical rows */
   for (crit = 0; crit < nvars; ++crit)
   {
      /* there are no fixed points */
      assert( perm[crit] != crit );

      /* initialize temporary solution */
      for (i = 0; i < nvars; ++i)
         tmpsolu[i] = 2;
      tmpsoluobj = 0.0;

      /* perform fixings implied by the critical row */
      tmpsolu[crit] = 0;
      assert( invperm[crit] < nvars );

      tmpsolu[invperm[crit]] = 1;
      tmpsoluobj += sepaobjective[invperm[crit]];

      /* perform 1-fixings */
      i = invperm[crit];
      while ( i < crit )
      {
         i = invperm[i];
         tmpsolu[i] = 1;
         tmpsoluobj += sepaobjective[i];
      }

      /* row c cannot be critical */
      if ( i == crit )
         continue;

      assert( tmpsolu[crit] == 0 );

      /* perform 0-fixing */
      i = perm[crit];
      while ( i < crit )
      {
         tmpsolu[i] = 0;
         i = perm[i];
      }

      /* iterate over rows above the critical row */
      for (i = 0; i < crit; ++i)
      {
         SCIP_Real objimpact = 0.0;
         int j;

         /* skip already fixed entries */
         if ( tmpsolu[i] != 2 )
            continue;

         /* Check effect of fixing entry i to 1 and apply all implied fixing to other entries.
          *
          * Observe: Experiments indicate that entries are more often fixed to 1 than to 0.
          * For this reason, we apply the 1-fixings directly. If it turns out that the 1-fixings
          * have a negative impact on the objective, we undo these fixings afterwards and apply
          * 0-fixings instead. */

         /* check fixings in invperm direction */
         j = i;
         do
         {
            assert( tmpsolu[j] == 2 );
            tmpsolu[j] = 1;
            objimpact += sepaobjective[j];
            j = invperm[j];
         }
         while ( j < crit && j != i );

         /* if we do not detect a cycle */
         if ( j != i )
         {
            /* fix entry j since this is not done in the above do-while loop */
            assert( tmpsolu[j] == 2 );
            tmpsolu[j] = 1;
            objimpact += sepaobjective[j];

            /* check fixings in perm direction */
            j = perm[i];
            while ( j < crit )
            {
               assert( j != i );
               assert( tmpsolu[j] == 2 );
               tmpsolu[j] = 1;
               objimpact += sepaobjective[j];
               j = perm[j];
            }

            assert( j != crit );
         }

         /* if fixing entry i has a positive impact -> keep above fixings of entries to 1 */
         /* otherwise -> reset entries to 0 */
         if ( SCIPisEfficacious(scip, objimpact) )
            tmpsoluobj += objimpact;
         else
         {
            j = i;
            do
            {
               assert( tmpsolu[j] == 1 );
               tmpsolu[j] = 0;
               j = invperm[j];
            }
            while ( j < crit && j != i );

            /* if we do not detect a cycle */
            if ( j != i )
            {
               /* fix entry j since this is not done in the above do-while loop */
               assert( tmpsolu[j] == 1 );
               tmpsolu[j] = 0;

               /* check fixings in perm direction */
               j = perm[i];
               while ( j < crit )
               {
                  assert( j != i );
                  assert( tmpsolu[j] == 1 );
                  tmpsolu[j] = 0;
                  j = perm[j];
               }

               assert( j != crit );
            }
         }
      }

      /* iterate over unfixed entries below the critical row */
      for (i = crit + 1; i < nvars; ++i)
      {
         /* skip already fixed entries */
         if ( tmpsolu[i] != 2 )
            continue;

         if ( SCIPisEfficacious(scip, sepaobjective[i]) )
         {
            assert( tmpsolu[i] == 2 );
            tmpsolu[i] = 1;
            tmpsoluobj += sepaobjective[i];
         }
         else
         {
            assert( tmpsolu[i] == 2 );
            tmpsolu[i] = 0;
         }
      }

      /* check whether we have found a better solution which has positive separation objective*/
      if ( SCIPisEfficacious(scip, tmpsoluobj + constobjective - maxsoluobj) )
      {
         assert( SCIPisEfficacious(scip, tmpsoluobj + constobjective) );
         for (i = 0; i < nvars; ++i)
            maxsolu[i] = tmpsolu[i];
         maxsoluobj = tmpsoluobj + constobjective;
      }
   }

   /* Check whether the separation objective is positive, i.e., a violated cover was found. */
   if ( SCIPisEfficacious(scip, maxsoluobj) )
   {
      SCIP_Real rhs = -1.0;
      SCIP_Real lhs = 0.0;

      for (i = 0; i < nvars; ++i)
      {
         if ( i < perm[i] )
         {
            maxsolu[i] = maxsolu[i] - 1;
            lhs += vals[i] * maxsolu[i];
         }
         else
         {
            lhs += vals[i] * maxsolu[i];
            rhs += maxsolu[i];
         }
      }

      assert( SCIPisGT(scip, lhs, rhs) );

      /* add cover inequality */
      SCIP_CALL( addSymresackInequality(scip, cons, nvars, consdata->vars, maxsolu, rhs, infeasible) );

      if ( ! *infeasible )
         ++(*ngen);
   }

   SCIPfreeBufferArrayNull(scip, &maxsolu);
   SCIPfreeBufferArrayNull(scip, &tmpsolu);
   SCIPfreeBufferArrayNull(scip, &sepaobjective);

   return SCIP_OKAY;
}


static
SCIP_RETCODE checkSymresackSolution(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constrained for which we check the solution */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_RESULT*          result,             /**< pointer to store whether we detected infeasibility */
   SCIP_Bool             printreason         /**< whether reason for infeasibility should be printed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* invperm;
   int nvars;
   int i;

   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);

   /* we don't have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->invperm != NULL );

   SCIPdebugMsg(scip, "Check method for symresack constraint <%s> (%d rows) ...\n", SCIPconsGetName(cons), consdata->nvars);

   nvars = consdata->nvars;
   vars = consdata->vars;
   invperm = consdata->invperm;

   /* detect first non-constant pair of variables */
   for (i = 0; i < nvars; ++i)
   {
      SCIP_Real solval;
      int val1;
      int val2;

      /* there are no fixed points */
      assert( invperm[i] != i );

      /* get value of variable i and its inverse */
      solval = SCIPgetSolVal(scip, sol, vars[i]);
      assert( SCIPisFeasIntegral(scip, solval) );
      if ( solval > 0.5 )
         val1 = 1;
      else
         val1 = 0;

      solval = SCIPgetSolVal(scip, sol, vars[invperm[i]]);
      assert( SCIPisFeasIntegral(scip, solval) );
      if ( solval > 0.5 )
         val2 = 1;
      else
         val2 = 0;

      /* if we detected a constant pair */
      if ( val1 == val2 )
         continue;
      /* pair is (1,0) --> lexicographically maximal */
      else if ( val1 > val2 )
         break;

      /* pair is (0,1) --> solution is infeasible */
      assert( val2 > val1 );
      SCIPdebugMsg(scip, "Solution is infeasible.\n");
      *result = SCIP_INFEASIBLE;

      if ( printreason )
         SCIPinfoMessage(scip, NULL, "First non-constant pair (%d, %d) of variables has pattern (0,1).\n", i, invperm[i]);

      break;
   }

   return SCIP_OKAY;
}


/** Upgrade symresack constraints to orbisacks */
static
SCIP_RETCODE orbisackUpgrade(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            inputvars,          /**< permuted variables array */
   int                   nvars,              /**< size of perm array */
   SCIP_Bool*            upgrade,            /**< whether constraint was upgraded */
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
   SCIP_VAR** vars1;
   SCIP_VAR** vars2;
   int nrows = 0;
   int i;

   assert( scip != NULL );
   assert( perm != NULL );
   assert( nvars > 0 );
   assert( inputvars != NULL );
   assert( upgrade != NULL );

   *upgrade = TRUE;

   /* check whether orbisack conshdlr is available */
   conshdlr = SCIPfindConshdlr(scip, "orbisack");
   if ( conshdlr == NULL )
   {
      *upgrade = FALSE;
      SCIPdebugMsg(scip, "Cannot check whether symresack constraint can be upgraded to orbisack constraint. ");
      SCIPdebugMsg(scip, "---> Orbisack constraint handler not found.\n");

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vars1, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars2, nvars) );

   /* check whether permutation is a composition of 2-cycles */
   for (i = 0; i < nvars; ++i)
   {
      /* ignore non-binary variables */
      if ( ! SCIPvarIsBinary(inputvars[i]) )
         continue;

      if ( perm[perm[i]] != i )
      {
         *upgrade = FALSE;
         break;
      }

      if ( perm[i] > i )
      {
         vars1[nrows] = inputvars[i];
         vars2[nrows++] = inputvars[perm[i]];

         assert( nrows <= nvars );
      }
   }

   /* if permutation can be upgraded to an orbisack */
   if ( nrows == 0 )
      *upgrade = FALSE;
   else if ( *upgrade )
   {
      SCIP_CALL( SCIPcreateConsOrbisack(scip, cons, name, vars1, vars2, nrows, FALSE, FALSE,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   SCIPfreeBufferArray(scip, &vars2);
   SCIPfreeBufferArray(scip, &vars1);

   return SCIP_OKAY;
}


/** creates a symmetry breaking constraint
 *
 * Depending on the given permutation, either an orbisack or symresack constraint
 * is created.
 */
SCIP_RETCODE SCIPcreateSymbreakCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
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
   SCIP_Bool upgrade = FALSE;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( perm != NULL );
   assert( vars != NULL );
   assert( nvars > 0 );

   SCIP_CALL( orbisackUpgrade(scip, cons, name, perm, vars, nvars, &upgrade,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   if ( ! upgrade )
   {
      SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, perm, vars, nvars,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }


   return SCIP_OKAY;
}


/*--------------------------------------------------------------------------------------------
 *--------------------------------- SCIP functions -------------------------------------------
 *--------------------------------------------------------------------------------------------*/

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSymresack)
{   /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( consdata != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}


/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeSymresack)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSymresack)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* consdata = NULL;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMsg(scip, "Transforming constraint.\n");

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);
   assert( sourcedata->nvars != 0 );
   assert( sourcedata->vars != NULL );
   assert( sourcedata->perm != NULL );
   assert( sourcedata->invperm != NULL );
#ifndef NDEBUG
   if ( sourcedata->ppupgrade )
   {
      assert( sourcedata->ncycles != 0 );
      assert( sourcedata->cycledecomposition != NULL );
      for (i = 0; i < sourcedata->ncycles; ++i)
      {
         assert( sourcedata->cycledecomposition[i] != NULL );
         assert( sourcedata->cycledecomposition[i][0] != 0 );
      }
   }
#endif

   /* create transformed constraint data (copy data where necessary) */
   nvars = sourcedata->nvars;

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

#ifdef SCIP_DEBUG
   consdata->debugcnt = sourcedata->debugcnt;
#endif
   consdata->nvars = nvars;

   if ( nvars > 0 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars) );
      SCIP_CALL( SCIPgetTransformedVars(scip, nvars, sourcedata->vars, consdata->vars) );
      for (i = 0; i < nvars; ++i)
      {
         SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
      }

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->perm, sourcedata->perm, nvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->invperm, sourcedata->invperm, nvars) );

      consdata->ppupgrade = sourcedata->ppupgrade;

      if ( sourcedata->ppupgrade )
      {
         consdata->ncycles = sourcedata->ncycles;
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycledecomposition, sourcedata->cycledecomposition, sourcedata->ncycles) );
         for (i = 0; i < sourcedata->ncycles; ++i)
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->cycledecomposition[i], sourcedata->cycledecomposition[i], nvars + 1) ); /*lint !e866*/
         }
      }
   }

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
SCIP_DECL_CONSINITLP(consInitlpSymresack)
{
   int c;

   assert( infeasible != NULL );
   *infeasible = FALSE;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      assert( conss[c] != NULL );

      SCIPdebugMsg(scip, "Generating initial symresack cut for constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( initLP(scip, conss[c], infeasible) );
      if ( *infeasible )
         break;
   }
   SCIPdebugMsg(scip, "Generated initial symresack cuts.\n");

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int ntotalvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   /* if solution is not integer */
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   if ( nconss == 0 )
      return SCIP_OKAY;

   ntotalvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, ntotalvars) );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      SCIPdebugMsg(scip, "Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      assert( consdata->nvars <= ntotalvars );
      SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, vals) );
      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         SCIPfreeBufferArray(scip, &vals);

         return SCIP_OKAY;
      }

      if ( ngen > 0 )
         *result = SCIP_SEPARATED;

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution */
static
SCIP_DECL_CONSSEPASOL(consSepasolSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int ntotalvars;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation method for symresack constraints\n");

   *result = SCIP_DIDNOTRUN;

   if ( nconss == 0 )
      return SCIP_OKAY;

   ntotalvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, ntotalvars) );

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      SCIPdebugMsg(scip, "Separating symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      /* get data of constraint */
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);

      /* get solution */
      assert( consdata->nvars <= ntotalvars );
      SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, vals) );
      SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );


      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         SCIPfreeBufferArray(scip, &vals);

         return SCIP_OKAY;
      }

      if ( ngen > 0 )
         *result = SCIP_SEPARATED;

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions.
 *
 *  To check feasibility, we separate cover inequalities.
 *
 *  @pre It is assumed that the solution is integral (this can be ensured by appropriate priorities).
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (lp solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   if ( nconss > 0 )
   {
      SCIP_Real* vals;
      int ntotalvars;

      ntotalvars = SCIPgetNVars(scip);
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, ntotalvars) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_Bool infeasible = FALSE;
         int ngen = 0;

         SCIPdebugMsg(scip, "Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);

         /* get solution */
         assert( consdata->nvars <= ntotalvars );
         SCIP_CALL( SCIPgetSolVals(scip, NULL, consdata->nvars, consdata->vars, vals) );
         SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIPfreeBufferArray(scip, &vals);

            return SCIP_OKAY;
         }

         /* SCIPdebugMsg(scip, "Generated symresack inequalities for <%s>: %d\n", SCIPconsGetName(conss[c]), ngen); */

         if ( ngen > 0 )
            *result = SCIP_SEPARATED;
      }
      SCIPfreeBufferArray(scip, &vals);
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSymresack)
{  /*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (pseudo solutions) ...\n");

   *result = SCIP_FEASIBLE;

   if ( objinfeasible || solinfeasible )
      return SCIP_OKAY;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CALL( checkSymresackSolution(scip, conss[c], NULL, result, FALSE) );

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions
 *
 *  To check feasibility, we separate cover inequalities.
 *
 */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSymresack)
{   /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Enforcing method for symresack constraints (relaxation solutions) ...\n");

   /* we have a negative priority, so we should come after the integrality conshdlr. */
   assert( SCIPgetNLPBranchCands(scip) == 0 );

   *result = SCIP_FEASIBLE;

   if ( nconss > 0 )
   {
      SCIP_Real* vals;
      int ntotalvars;

      ntotalvars = SCIPgetNVars(scip);
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, ntotalvars) );

      /* loop through constraints */
      for (c = 0; c < nconss; ++c)
      {
         SCIP_Bool infeasible = FALSE;
         int ngen = 0;

         SCIPdebugMsg(scip, "Enforcing symresack constraint <%s> ...\n", SCIPconsGetName(conss[c]));

         /* get data of constraint */
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);

          /* get solution */
         assert( consdata->nvars <= ntotalvars );
         SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, vals) );
         SCIP_CALL( separateSymresackCovers(scip, conss[c], consdata, vals, &ngen, &infeasible) );

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            SCIPfreeBufferArray(scip, &vals);

            return SCIP_OKAY;
         }

         if ( ngen > 0 )
            *result = SCIP_SEPARATED;
      }
      SCIPfreeBufferArray(scip, &vals);
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSymresack)
{   /*lint --e{715}*/
   int c;
   SCIP_CONSHDLRDATA* conshdlrdata;

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
      SCIP_CALL( checkSymresackSolution(scip, conss[c], sol, result, printreason) );

      if ( *result == SCIP_INFEASIBLE )
         break;
   }

   if ( *result == SCIP_FEASIBLE )
      SCIPdebugMsg(scip, "Solution is feasible.\n");

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSymresack)
{  /*lint --e{715}*/
   int c;
   SCIP_Bool success = FALSE;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Propagation method of symresack constraint handler.\n");

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      int ngen = 0;

      assert( conss[c] != NULL );

      SCIP_CALL( propVariables(scip, conss[c], &infeasible, &ngen) );

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      success = success || ( ngen > 0 );

      *result = SCIP_DIDNOTFIND;
   }

   if ( success )
   {
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSymresack)
{  /*lint --e{715}*/
   int c;
   SCIP_Bool success = FALSE;
   int oldndelconss;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   oldndelconss = *ndelconss;

   SCIPdebugMsg(scip, "Presolving method of symresack constraint handler. Propagating symresack inequalities.\n");
   *result = SCIP_DIDNOTRUN;

   /* loop through constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_Bool infeasible = FALSE;
      SCIP_CONSDATA* consdata;
      int ngen = 0;

      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* avoid trivial problems */
      if ( consdata->nvars == 0 )
      {
         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         (*ndelconss)++;
      }
      else
      {
         SCIP_CALL( propVariables(scip, conss[c], &infeasible, &ngen) );
      }

      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if ( ngen > 0 )
      {
         *nfixedvars += ngen;
         success = TRUE;
      }

      *result = SCIP_DIDNOTFIND;
   }

   if ( *ndelconss > oldndelconss ||  success )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** Propagation resolution for conflict analysis */
static
SCIP_DECL_CONSRESPROP(consRespropSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* invperm;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Propagation resolution method of symresack constraint handler.\n");

   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* we don't have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->invperm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   invperm = consdata->invperm;

   assert( 0 <= inferinfo && inferinfo < (2 * nvars - 1) );

   /* if first part of variable pair was fixed to 0 */
   if ( inferinfo < nvars )
   {
      assert( vars[invperm[inferinfo]] == infervar );
      assert( SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, FALSE) > 0.5
         && SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, TRUE) < 0.5 );

      if ( SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, FALSE) > 0.5
         && SCIPvarGetUbAtIndex(vars[invperm[inferinfo]], bdchgidx, TRUE) < 0.5 )
      {
         SCIPdebugMsg(scip, " -> reason for setting x[%d] = 0 was fixing x[%d] to 0 ", invperm[inferinfo], inferinfo);
         SCIPdebugMsg(scip, "and each pair of binary variables before (%d,%d) which are not fixed points is constant.\n",
            inferinfo, invperm[inferinfo]);

         SCIP_CALL( SCIPaddConflictUb(scip, vars[inferinfo], bdchgidx) );

         for (i = 0; i < inferinfo; ++i)
         {
            /* there are no fixed points */
            assert( invperm[i] != i );

            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invperm[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[i]], bdchgidx) );
         }
      }
   }
   /* if second part of variable pair was fixed to 1 */
   else
   {
      int inferinfo2;

      inferinfo2 = inferinfo - nvars;
      assert( vars[inferinfo2] == infervar );
      assert( SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, FALSE) < 0.5
         && SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, TRUE) > 0.5 );

      if ( SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, FALSE) < 0.5
         && SCIPvarGetLbAtIndex(vars[inferinfo2], bdchgidx, TRUE) > 0.5 )
      {
         SCIPdebugMsg(scip, " -> reason for setting x[%d] = 1 was fixing x[%d] to 1 ", inferinfo2, invperm[inferinfo2]);
         SCIPdebugMsg(scip, "and each pair of binary variables before (%d,%d) which are not fixed points is constant.\n",
            inferinfo2, invperm[inferinfo2]);

         SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[inferinfo2]], bdchgidx) );

         for (i = 0; i < inferinfo2; ++i)
         {
            /* there are no fixed points */
            assert( invperm[i] != i );

            SCIP_CALL( SCIPaddConflictUb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[i], bdchgidx) );
            SCIP_CALL( SCIPaddConflictUb(scip, vars[invperm[i]], bdchgidx) );
            SCIP_CALL( SCIPaddConflictLb(scip, vars[invperm[i]], bdchgidx) );
         }
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** lock variables
 *
 *  We assume we have only one global (void) constraint and lock all binary variables
 *  which do not correspond to fixed points of the permutation.
 *
 * - Symresack constraints may get violated if the variables with a negative coefficient
 *   in the FD inequality are rounded down, we therefor call
 *   SCIPaddVarLocks(..., nlockspos, nlocksneg).
 * - Symresack constraints may get violated if the variables with a positive coefficient
 *   in the FD inequality are rounded up, we therefor call
 *   SCIPaddVarLocks(..., nlocksneg, nlockspo ).
 */
static
SCIP_DECL_CONSLOCK(consLockSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int* perm;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "Locking method for symresack constraint handler.\n");

   /* get data of original constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* we don't have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
      return SCIP_OKAY;

   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );

   nvars = consdata->nvars;
   vars = consdata->vars;
   perm = consdata->perm;

   for (i = 0; i < nvars; ++i)
   {
      /* due to clean-up in consdataCreate, there are no fixed points */
      assert( perm[i] != i );

      if ( perm[i] > i )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlockspos, nlocksneg) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, vars[i], nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** constraint display method of constraint handler
 *
 *  The constraint handler should output a representation of the constraint into the given text file.
 */
static
SCIP_DECL_CONSPRINT(consPrintSymresack)
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   SCIP_Bool* covered;
   int* perm;
   int nvars;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMsg(scip, "Printing method for symresack constraint handler\n");

   /* we don't have to take care of trivial constraints */
   if ( consdata->nvars < 2 )
   {
      SCIPinfoMessage(scip, file, "symresack()");
      return SCIP_OKAY;
   }

   assert( consdata->vars != NULL );
   assert( consdata->perm != NULL );

   vars = consdata->vars;
   nvars = consdata->nvars;
   perm = consdata->perm;

   SCIP_CALL( SCIPallocBufferArray(scip, &covered, nvars) );
   for (i = 0; i < nvars; ++i)
      covered[i] = FALSE;

   if ( consdata->ppupgrade )
      SCIPinfoMessage(scip, file, "ppSymresack(");
   else
      SCIPinfoMessage(scip, file, "symresack(");

   for (i = 0; i < nvars; ++i)
   {
      if ( covered[i] )
         continue;

      /* print cycle of perm containing i */
      SCIPinfoMessage(scip, file, "[%s", SCIPvarGetName(vars[i]));
      covered[i] = TRUE;
      j = perm[i];
      while ( j != i )
      {
         SCIPinfoMessage(scip, file, ",%s", SCIPvarGetName(vars[j]));
         covered[j] = TRUE;
         j = perm[j];
      }
      SCIPinfoMessage(scip, file, "]");
   }
   SCIPinfoMessage(scip, file, ")");

   SCIPfreeBufferArray(scip, &covered);

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( vars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( varssize < consdata->nvars )
      (*success) = FALSE;
   else
   {
      int cnt = 0;
      int i;

      for (i = 0; i < consdata->nvars; ++i)
         vars[cnt++] = consdata->vars[i];
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSymresack)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   assert( success != NULL );
   assert( nvars != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   (*nvars) = consdata->nvars;
   (*success) = TRUE;

   return SCIP_OKAY;
}


/** creates the handler for symresack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSymresack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSymresack, consEnfopsSymresack, consCheckSymresack, consLockSymresack,
         conshdlrdata) );
   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSymresack) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSymresack) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSymresack) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSymresack) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSymresack) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSymresack, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSymresack) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSymresack, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSymresack) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSymresack, consSepasolSymresack, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSymresack) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpSymresack) );

   /* whether we allow upgrading to packing/partioning symresack constraints*/
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/ppsymresack",
         "Upgrade symresack constraints to packing/partioning symresacks?",
         &conshdlrdata->checkppsymresack, TRUE, DEFAULT_PPSYMRESACK, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkalwaysfeas",
         "Whether check routine returns always SCIP_FEASIBLE.",
         &conshdlrdata->checkalwaysfeas, TRUE, DEFAULT_CHECKALWAYSFEAS, NULL, NULL) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates and captures a symresack constraint
 *
 *  In a presolving step, we check whether the permutation acts only on binary points. Otherwise, we eliminate
 *  the non-binary variables from the permutation.
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars,              /**< number of variables in vars array */
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

   assert( cons != NULL );
   assert( nvars > 0 );

   /* find the symresack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Symresack constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, vars, nvars, perm) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate && (! consdata->ppupgrade), enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}


/** creates and captures a symresack constraint
 *  in its most basic variant, i.e., with all constraint flags set to their default values
 *
 *  In a presolving step, we remove all fixed points and cycles that act on non-binary variables of the permutation
 *
 *  @note The constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons().
 */
SCIP_RETCODE SCIPcreateConsBasicSymresack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int*                  perm,               /**< permutation */
   SCIP_VAR**            vars,               /**< variables */
   int                   nvars               /**< number of variables in vars array */
   )
{
   SCIP_CALL( SCIPcreateConsSymresack(scip, cons, name, perm, vars, nvars,
         TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
