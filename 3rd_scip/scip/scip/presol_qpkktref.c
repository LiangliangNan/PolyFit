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

/**@file   presol_qpkktref.c
 * @brief  qpkktref presolver
 * @author Tobias Fischer
 *
 * This presolver tries to add the KKT conditions as additional (redundant) constraints to the (mixed-binary) quadratic
 * program
 *  \f[
 *  \begin{array}{ll}
 *  \min         & x^T Q x + c^T x + d \\
 *               & A x \leq b, \\
 *               & x \in \{0, 1\}^{p} \times R^{n-p}.
 * \end{array}
 * \f]
 *
 * We first check if the structure of the program is like (QP), see the documentation of the function
 * checkConsQuadraticProblem().
 *
 * If the problem is known to be bounded (all variables have finite lower and upper bounds), then we add the KKT
 * conditions. For a continuous QPs the KKT conditions have the form
 * \f[
 *  \begin{array}{ll}
 *   Q x + c + A^T \mu = 0,\\
 *   Ax \leq b,\\
 *   \mu_i \cdot (Ax - b)_i = 0,    & i \in \{1, \dots, m\},\\
 *   \mu \geq 0.
 * \end{array}
 * \f]
 * where \f$\mu\f$ are the Lagrangian variables. Each of the complementarity constraints \f$\mu_i \cdot (Ax - b)_i = 0\f$
 * is enforced via an SOS1 constraint for \f$\mu_i\f$ and an additional slack variable \f$s_i = (Ax - b)_i\f$.
 *
 * For mixed-binary QPs, the KKT-like conditions are
 * \f[
 *  \begin{array}{ll}
 *   Q x + c + A^T \mu + I_J \lambda = 0,\\
 *   Ax \leq b,\\
 *   x_j \in \{0,1\}                    & j \in J,\\
 *   (1 - x_j) \cdot z_j = 0            & j \in J,\\
 *   x_j \cdot (z_j - \lambda_j) = 0    & j \in J,\\
 *   \mu_i \cdot (Ax - b)_i = 0         & i \in \{1, \dots, m\},\\
 *   \mu \geq 0,
 * \end{array}
 * \f]
 * where \f$J = \{1,\dots, p\}\f$, \f$\mu\f$ and \f$\lambda\f$ are the Lagrangian variables, and \f$I_J\f$ is the
 * submatrix of the \f$n\times n\f$ identity matrix with columns indexed by \f$J\f$. For the derivation of the KKT-like
 * conditions, see
 *
 *  Branch-And-Cut for Complementarity and Cardinality Constrained Linear Programs,@n
 *  Tobias Fischer, PhD Thesis (2016)
 *
 * Algorithmically:
 *
 * - we handle the quadratic term variables of the quadratic constraint like in the method
 *   presolveAddKKTQuadQuadraticTerms()
 * - we handle the bilinear term variables of the quadratic constraint like in the method presolveAddKKTQuadBilinearTerms()
 * - we handle the linear term variables of the quadratic constraint like in the method presolveAddKKTQuadLinearTerms()
 * - we handle linear constraints in the method presolveAddKKTLinearConss()
 * - we handle aggregated variables in the method presolveAddKKTAggregatedVars()
 *
 * we have a hashmap from each variable to the index of the dual constraint in the KKT conditions.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/presol_qpkktref.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_sos1.h"


#define PRESOL_NAME            "qpkktref"
#define PRESOL_DESC            "adds KKT conditions to (mixed-binary) quadratic programs"
#define PRESOL_PRIORITY              -1 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers);
                                         *   combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_MEDIUM /* timing of the presolver (fast, medium, or exhaustive) */


/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Bool             addkktbinary;       /**< if TRUE then allow binary variables for KKT update */
   SCIP_Bool             updatequadbounded;  /**< if TRUE then only apply the update to QPs with bounded variables; if
                                              *   the variables are not bounded then a finite optimal solution might not
                                              *   exist and the KKT conditions would then be invalid */
   SCIP_Bool             updatequadindef;    /**< if TRUE then apply quadratic constraint update even if the quadratic 
                                              *   constraint matrix is known to be indefinite */
};


/*
 * Local methods
 */

/** for a linear constraint \f$a^T x \leq b\f$, create the complementarity constraint \f$\mu \cdot s = 0\f$, where
 *  \f$s = b - a^T x\f$ and \f$\mu\f$ is the dual variable associated to the constraint \f$a^T x \leq b\f$
 */
static
SCIP_RETCODE createKKTComplementarityLinear(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           namepart,           /**< name of linear constraint */
   SCIP_VAR**            vars,               /**< variables of linear constraint */
   SCIP_Real*            vals,               /**< coefficients of variables in linear constraint */
   SCIP_Real             lhs,                /**< left hand side of linear constraint */
   SCIP_Real             rhs,                /**< right hand side of linear constraint */
   int                   nvars,              /**< number of variables of linear constraint */
   SCIP_VAR*             dualvar,            /**< dual variable associated to linear constraint */
   SCIP_Bool             takelhs,            /**< whether to consider the lhs or the rhs of the constraint */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* KKTlincons;
   SCIP_CONS* sos1cons;
   SCIP_VAR* slack;
   SCIP_Real slackcoef;
   SCIP_Real eqval;

   assert( scip != NULL );
   assert( namepart != NULL );
   assert( vars != NULL );
   assert( vals != NULL );
   assert( dualvar != NULL );
   assert( ! takelhs || ! SCIPisInfinity(scip, -lhs) );
   assert( takelhs || ! SCIPisInfinity(scip, rhs) );
   assert( naddconss != NULL );

   if( takelhs )
   {
      eqval = lhs;
      slackcoef = -1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_lhs_%s", namepart);
   }
   else
   {
      eqval = rhs;
      slackcoef = 1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_rhs_%s", namepart);
   }

   /* create slack variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &slack, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

   /* add skack variable */
   SCIP_CALL( SCIPaddVar(scip, slack) );

   /* create a new linear constraint */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTlin_%s_%d", namepart, takelhs);
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &KKTlincons, name, nvars, vars, vals, eqval, eqval) );

   /* add slack variable to linear constraint */
   SCIP_CALL( SCIPaddCoefLinear(scip, KKTlincons, slack, slackcoef) );

   /* create SOS1 (complementarity) constraint involving dual variable of linear constraint and slack variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_lin_%s_%d", namepart, takelhs);
   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons, name, 0, NULL, NULL) );

   /* add slack and dual variable to SOS1 constraint */
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, slack, 1.0) );
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, dualvar, 2.0) );

   /* add/release constraints */
   SCIP_CALL( SCIPaddCons(scip, sos1cons) );
   SCIP_CALL( SCIPaddCons(scip, KKTlincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &sos1cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &KKTlincons) );
   *naddconss = *naddconss + 2;

   /* release slack variable */
   SCIP_CALL( SCIPreleaseVar(scip, &slack) );

   return SCIP_OKAY;
}

/** create complementarity constraints of KKT conditions associated to bounds of variables
 * - for an upper bound constraint \f$x_i \leq u_i\f$, create the complementarity constraint \f$\mu_i \cdot s_i = 0\f$,
 *   where \f$s_i = u_i - x_i\f$ and \f$\mu_i\f$ is the dual variable of the upper bound constraint
 * - for a lower bound constraint \f$x_i \geq l_i\f$, create the complementarity constraint \f$\lambda_i \cdot w_i = 0\f$,
 *   where \f$w_i = x_i - l_i\f$
 *   and \f$\lambda_i\f$ is the dual variable of the lower bound constraint
 */
static
SCIP_RETCODE createKKTComplementarityBounds(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             dualvar,            /**< dual variable associated to bound of variable */
   SCIP_Bool             takelb,             /**< whether to consider the lower or upper bound of variable */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* KKTlincons;
   SCIP_CONS* sos1cons;
   SCIP_VAR* slack;
   SCIP_Real slackcoef;
   SCIP_Real eqval;

   assert( scip != NULL );
   assert( var != NULL );
   assert( dualvar != NULL );
   assert( ! takelb || ! SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) );
   assert( takelb || ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) );
   assert( naddconss != NULL );

   if( takelb )
   {
      eqval = SCIPvarGetLbGlobal(var);
      slackcoef = -1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_lb_%s", SCIPvarGetName(var));
   }
   else
   {
      eqval = SCIPvarGetUbGlobal(var);
      slackcoef = 1.0;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "slack_ub_%s", SCIPvarGetName(var));
   }

   /* create complementarity constraint; if bound is nonzero, we additionally need to introduce a slack variable */
   if( SCIPisFeasZero(scip, eqval) && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      /* create SOS1 (complementarity) constraint involving dual variable of linear constraint and slack variable */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bound%s_%d", SCIPvarGetName(var), takelb);
      SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons, name, 0, NULL, NULL) );

      /* add slack and dual variable to SOS1 constraint */
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, var, 1.0) );
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, dualvar, 2.0) );

      /* add/release constraint */
      SCIP_CALL( SCIPaddCons(scip, sos1cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &sos1cons) );
      ++(*naddconss);
   }
   else
   {
      /* create slack variable */
      SCIP_CALL( SCIPcreateVarBasic(scip, &slack, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

      /* add skack variable */
      SCIP_CALL( SCIPaddVar(scip, slack) );

      /* create a new linear constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKT_bound%s_%d", SCIPvarGetName(var), takelb);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &KKTlincons, name, 0, NULL, NULL, eqval, eqval) );

      /* add slack variable to linear constraint */
      SCIP_CALL( SCIPaddCoefLinear(scip, KKTlincons, var, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, KKTlincons, slack, slackcoef) );

      /* create SOS1 (complementarity) constraint involving dual variable of linear constraint and slack variable */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bound%s_%d", SCIPvarGetName(var), takelb);
      SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons, name, 0, NULL, NULL) );

      /* add slack and dual variable to SOS1 constraint */
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, slack, 1.0) );
      SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons, dualvar, 2.0) );

      /* add/release constraints */
      SCIP_CALL( SCIPaddCons(scip, sos1cons) );
      SCIP_CALL( SCIPaddCons(scip, KKTlincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &sos1cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &KKTlincons) );
      *naddconss = *naddconss + 2;

      /* release slack variable */
      SCIP_CALL( SCIPreleaseVar(scip, &slack) );
   }

   return SCIP_OKAY;
}

/** create the complementarity constraints of the KKT-like conditions associated to a binary variable \f$x_i\f$;
 *  these are \f$(1 - x_i) \cdot z_i = 0\f$ and \f$x_i \cdot (z_i - \lambda_i) = 0\f$, where \f$z_i\f$ and
 *  \f$\lambda_i\f$ are dual variables
 */
static
SCIP_RETCODE createKKTComplementarityBinary(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR*             var,                /**< variable */
   SCIP_VAR*             dualbin1,           /**< first dual variable associated to binary variable */
   SCIP_VAR*             dualbin2,           /**< second dual variable associated to binary variable */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* conslinbin1;
   SCIP_CONS* conslinbin2;
   SCIP_CONS* sos1cons1;
   SCIP_CONS* sos1cons2;
   SCIP_VAR* slackbin1;
   SCIP_VAR* slackbin2;

   assert( scip != NULL );
   assert( var != NULL );
   assert( dualbin1 != NULL );
   assert( dualbin2 != NULL );
   assert( naddconss != NULL );

   /* create first slack variable associated to binary constraint; domain [-inf, inf] */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_slackbin1", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateVarBasic(scip, &slackbin1, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
        SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, slackbin1) );
   assert( slackbin1 != NULL );

   /* create a new linear constraint: dualbin1 - dualbin2 = slackbin */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTBinary1_%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conslinbin1, name, 0, NULL, NULL, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin1, dualbin1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin1, dualbin2, -1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin1, slackbin1, -1.0) );
   SCIP_CALL( SCIPaddCons(scip, conslinbin1) );
   SCIP_CALL( SCIPreleaseCons(scip, &conslinbin1) );
   ++(*naddconss);

   /* create SOS1 (complementarity) constraint involving binary variable and slack variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bin1%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons1, name, 0, NULL, NULL) );

   /* add slack and dual variable to SOS1 constraint */
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons1, var, 1.0) );
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons1, slackbin1, 2.0) );

   /* add/release constraint */
   SCIP_CALL( SCIPaddCons(scip, sos1cons1) );
   SCIP_CALL( SCIPreleaseCons(scip, &sos1cons1) );
   ++(*naddconss);

   /* release slack variable */
   SCIP_CALL( SCIPreleaseVar(scip, &slackbin1) );


   /* create second slack variable associated to binary constraint; domain [0, inf] */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_slackbin2", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateVarBasic(scip, &slackbin2, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, slackbin2) );
   assert( slackbin2 != NULL );

   /* create a new linear constraint: 1.0 - var = slackbin2 */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTBinary2_%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &conslinbin2, name, 0, NULL, NULL, 1.0, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin2, var, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, conslinbin2, slackbin2, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, conslinbin2) );
   SCIP_CALL( SCIPreleaseCons(scip, &conslinbin2) );
   ++(*naddconss);

   /* create SOS1 (complementarity) constraint involving first dual variable and slack variable */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTsos1_bin2%s", SCIPvarGetName(var));
   SCIP_CALL( SCIPcreateConsBasicSOS1(scip, &sos1cons2, name, 0, NULL, NULL) );

   /* add slack and dual variable to SOS1 constraint */
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons2, dualbin1, 1.0) );
   SCIP_CALL( SCIPaddVarSOS1(scip, sos1cons2, slackbin2, 2.0) );

   /* add/release constraint */
   SCIP_CALL( SCIPaddCons(scip, sos1cons2) );
   SCIP_CALL( SCIPreleaseCons(scip, &sos1cons2) );
   ++(*naddconss);

   /* release slack variable */
   SCIP_CALL( SCIPreleaseVar(scip, &slackbin2) );

   return SCIP_OKAY;
}

/** create/get dual constraint of KKT conditions associated to primal variable @n@n
 * if variable does not already exist in hashmap then
 * 1. create dual constraint for variable
 * 2. create a dual variable \f$\mu_i\f$ for the upper bound constraint \f$x_i \leq u_i\f$
 * 3. create a dual variable \f$\lambda_i\f$ for the lower bound constraint \f$x_i \geq l_i\f$
 * 4. create the complementarity constraint \f$\mu_i \cdot s_i = 0\f$, where \f$s_i = u_i - x_i\f$
 * 5. create the complementarity constraint \f$\lambda_i \cdot w_i = 0\f$, where \f$w_i = x_i - l_i\f$
 * 6. add objective coefficients of dual variables
 * 7. the treatment of binary variables needs special care see the documentation of createKKTComplementarityBinary()
 *
 * if variable exists in hasmap then the dual constraint associated to the variable has already been created and is returned
 */
static
SCIP_RETCODE createKKTDualCons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   SCIP_CONS**           dualcons,           /**< dual constraint associated to variable */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   SCIP_VAR* dualub = NULL;     /* dual variable associated to upper bound constraint */
   SCIP_VAR* duallb = NULL;     /* dual variable associated to lower bound constraint */
   SCIP_VAR* dualbin1 = NULL;   /* first dual variable associated to binary variable */
   SCIP_VAR* dualbin2 = NULL;   /* second dual variable associated to binary variable */

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( var != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );

   /* if variable exists in hashmap */
   if( SCIPhashmapExists(varhash, var) )
   {
      int ind;
      ind = (int) (size_t) SCIPhashmapGetImage(varhash, var);
      *dualcons = dualconss[ind];
   }
   else
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      /* create dual variables corresponding to the bounds of the variables; binary variables have to be treated in a
       * different way */
      if( SCIPvarIsBinary(var) )
      {
         /* create first dual variable associated to binary constraint; the domain of dualbin is [-inf,inf]; the objective
          * coefficient is -0.5 */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_bin1", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVarBasic(scip, &dualbin1, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
              SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, dualbin1) );
         assert( dualbin1 != NULL );
         SCIP_CALL( SCIPaddCoefLinear(scip, objcons, dualbin1, -0.5) );

         /* create second variable associated to binary constraint; the domain of dualbin2 is [-inf,inf]; the objective
          * coefficient is zero */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_bin2", SCIPvarGetName(var));
         SCIP_CALL( SCIPcreateVarBasic(scip, &dualbin2, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
              SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, dualbin2) );
         assert( dualbin2 != NULL );
      }
      else
      {
         if( ! SCIPisInfinity(scip, -lb) )
         {
            /* create dual variable associated to lower bound; the domain of duallb is [0,inf] */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_lb", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateVarBasic(scip, &duallb, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, duallb) );
            assert( duallb != NULL );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallb, 0.5 * lb) );
         }

         if( ! SCIPisInfinity(scip, ub) )
         {
            /* create dual variable associated to upper bound; the domain of dualub is [0,inf] */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_ub", SCIPvarGetName(var));
            SCIP_CALL( SCIPcreateVarBasic(scip, &dualub, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, dualub) );
            assert( dualub != NULL );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, dualub, -0.5 * ub) );
         }
      }

      /* add variable in map  */
      SCIP_CALL( SCIPhashmapInsert(varhash, var, (void*) (size_t) *ndualconss) );/*lint !e571*/
      assert( *ndualconss == (int) (size_t) SCIPhashmapGetImage(varhash, var) );

      /* create a new linear constraint */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "KKTref_%s", SCIPvarGetName(var));
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, dualcons, name, 0, NULL, NULL, 0.0, 0.0) );

      /* add dual constraint to array for later use */
      dualconss[(*ndualconss)++] = *dualcons;


      /* add dual variables to dual constraints and create complementarity constraints; binary variables have to be
       * treated in a different way */
      if( SCIPvarIsBinary(var) )
      {
         /* add coefficient of second dual variable corresponding to binary variable */
         SCIP_CALL( SCIPaddCoefLinear(scip, *dualcons, dualbin2, 1.0) );

         /* create complementarity constraints */
         SCIP_CALL( createKKTComplementarityBinary(scip, var, dualbin1, dualbin2, naddconss) );

         SCIP_CALL( SCIPreleaseVar(scip, &dualbin1) );
         SCIP_CALL( SCIPreleaseVar(scip, &dualbin2) );
      }
      else
      {
         if( duallb != NULL )
         {
            /* add dual variable corresponding to lower bound of variable */
            SCIP_CALL( SCIPaddCoefLinear(scip, *dualcons, duallb, -1.0) );

            /* create complementarity constraint between slack variable of lower bound constraint and dual variable of
             * lower bound */
            SCIP_CALL( createKKTComplementarityBounds(scip, var, duallb, TRUE, naddconss) );

            SCIP_CALL( SCIPreleaseVar(scip, &duallb) );
         }

         if( dualub != NULL )
         {
            /* add dual variable corresponding to upper bound of variable */
            SCIP_CALL( SCIPaddCoefLinear(scip, *dualcons, dualub, 1.0) );

            /* create complementarity constraint between slack variable of upper bound constraint and dual variable of
             * upper bound */
            SCIP_CALL( createKKTComplementarityBounds(scip, var, dualub, FALSE, naddconss) );

            SCIP_CALL( SCIPreleaseVar(scip, &dualub) );
         }
      }
   }
   assert( *dualcons != NULL );

   return SCIP_OKAY;
}

/** handle (a single) linear constraint for quadratic constraint update
 * 1. create the dual constraints (i.e., the two rows of \f$Q x + c + A^T \mu = 0\f$) associated to the variables of the
 *    linear constraint, if not done already
 * 2. create the dual variables and the complementarity constraints for the lower and upper bound constraints of the
 *    variables of the linear constraint, if not done already
 * 3. create the dual variable \f$\mu_i\f$ associated to this linear constraint
 * 4. create the complementarity constraint \f$\mu_i \cdot (Ax - b)_i = 0\f$ associated to this linear constraint
 * 5. add objective coefficients of dual variables
 *
 * for steps 1 and 2 see the documentation of createKKTDualCons() for further information.@n
 * for step 4 see the documentation of the function createKKTComplementarityLinear() for further information.
 */
static
SCIP_RETCODE presolveAddKKTLinearCons(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   const char*           namepart,           /**< name of linear constraint */
   SCIP_VAR**            vars,               /**< variables of linear constraint */
   SCIP_Real*            vals,               /**< coefficients of variables in linear constraint */
   SCIP_Real             lhs,                /**< left hand side of linear constraint */
   SCIP_Real             rhs,                /**< right hand side of linear constraint */
   int                   nvars,              /**< number of variables of linear constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   int i;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( namepart != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( vars != NULL );
   assert( vals != NULL );
   assert( namepart != NULL );
   assert( naddconss != NULL );

   /* differ between left hand side and right hand side case (i=0 -> lhs; i=1 -> rhs) */
   for( i = 0; i < 2; ++i )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR* duallin = NULL;
      int j;

      /* skip one iteration if lhs equals rhs */
      if( i == 0 && SCIPisFeasEQ(scip, lhs, rhs) )
         continue;

      /* create dual variable corresponding to linear constraint */
      if( i == 0 )
      {
         assert( ! SCIPisFeasEQ(scip, lhs, rhs) );

         if( SCIPisInfinity(scip, -lhs) )
            continue;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_lhs", namepart);
         SCIP_CALL( SCIPcreateVarBasic(scip, &duallin, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, duallin) );
         SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallin, 0.5 * lhs) );

         /* create complementarity constraint between dual variable and slack variable of linear constraint */
         SCIP_CALL( createKKTComplementarityLinear(scip, namepart, vars, vals, lhs, rhs, nvars, duallin, TRUE,
              naddconss) );
      }
      else
      {
         if( SCIPisInfinity(scip, rhs) )
            continue;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dual_%s_rhs", namepart);
         if( SCIPisFeasEQ(scip, lhs, rhs) )
         {
            SCIP_CALL( SCIPcreateVarBasic(scip, &duallin, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
                 SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, duallin) );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallin, -0.5 * rhs) );
         }
         else
         {
            SCIP_CALL( SCIPcreateVarBasic(scip, &duallin, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, duallin) );
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, duallin, -0.5 * rhs) );

            /* create complementarity constraint between dual variable and slack variable of linear constraint */
            SCIP_CALL( createKKTComplementarityLinear(scip, namepart, vars, vals, lhs, rhs, nvars, duallin, FALSE,
                 naddconss) );
         }
      }
      assert( duallin != NULL );


      /* loop through variables of linear constraint */
      for( j = 0; j < nvars; ++j )
      {
         SCIP_CONS* dualcons = NULL;  /* dual constraint associated to variable */
         SCIP_VAR* var;

         var = vars[j];

         /* create/get dual constraint associated to variable;
          * if variable does not already exist in hashmap then create dual variables for its bounds */
         SCIP_CALL( createKKTDualCons(scip, objcons, var, varhash, dualconss, ndualconss, &dualcons, naddconss) );
         assert( dualcons != NULL );

         /* add dual variable corresponding to linear constraint */
         if( i == 0 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, dualcons, duallin, -vals[j]) );
         }
         else
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, dualcons, duallin, vals[j]) );
         }
      }

      /* release dual variable */
      SCIP_CALL( SCIPreleaseVar(scip, &duallin) );
   }

   return SCIP_OKAY;
}

/** handle linear constraints for quadratic constraint update, see the documentation of the function
 *  presolveAddKKTLinearCons() for an explanation
 */
static
SCIP_RETCODE presolveAddKKTLinearConss(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_CONS**           savelinconss,       /**< copy of array with linear constraints */
   int                   nlinconss,          /**< number of linear constraints */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss,          /**< buffer to increase with number of created additional constraints */
   int*                  ndelconss           /**< buffer to increase with number of deleted constraints */
   )
{
   int c;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );
   assert( ndelconss != NULL );

   /* loop through linear constraints */
   for( c = 0; c < nlinconss; ++c )
   {
      SCIP_CONS* lincons;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nvars;

      /* get data of constraint */
      lincons = savelinconss[c];
      assert( lincons != NULL );
      lhs = SCIPgetLhsLinear(scip, lincons);
      rhs = SCIPgetRhsLinear(scip, lincons);
      nvars = SCIPgetNVarsLinear(scip, lincons);
      vars = SCIPgetVarsLinear(scip, lincons);
      vals = SCIPgetValsLinear(scip, lincons);

      /* handle linear constraint for quadratic constraint update */
      SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(lincons),
            vars, vals, lhs, rhs, nvars, varhash, dualconss, ndualconss, naddconss) );
   }

   /* remove linear constraints if lhs != rhs, since they are now redundant; their feasibility is already expressed
    * by s >= 0, where s is the new slack variable that we introduced for these linear constraints */
   for( c = nlinconss-1; c >= 0; --c )
   {
      SCIP_CONS* lincons;

      lincons = savelinconss[c];
      assert( savelinconss[c] != NULL );

      if( ! SCIPisFeasEQ(scip, SCIPgetLhsLinear(scip, lincons), SCIPgetRhsLinear(scip, lincons)) )
      {
         SCIP_CALL( SCIPdelCons(scip, savelinconss[c]) );
         ++(*ndelconss);
      }
   }

   return SCIP_OKAY;
}

/** handle knapsack constraints for quadratic constraint update, see the documentation of the function
 *  presolveAddKKTLinearCons() for an explanation
 */
static
SCIP_RETCODE presolveAddKKTKnapsackConss(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss,          /**< buffer to increase with number of created additional constraints */
   int*                  ndelconss           /**< buffer to increase with number of deleted constraints */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );
   assert( ndelconss != NULL );

   conshdlr = SCIPfindConshdlr(scip, "knapsack");
   if( conshdlr == NULL )
      return SCIP_OKAY;

   nconss = SCIPconshdlrGetNConss(conshdlr);
   conss = SCIPconshdlrGetConss(conshdlr);

   /* loop through knapsack constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Longint* weights;
      SCIP_Real* vals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nvars;
      int v;

      /* get data of constraint */
      cons = conss[c];
      assert( cons != NULL );
      lhs = -SCIPinfinity(scip);
      rhs = (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons);
      nvars = SCIPgetNVarsKnapsack(scip, cons);
      vars = SCIPgetVarsKnapsack(scip, cons);
      weights = SCIPgetWeightsKnapsack(scip, cons);

      /* set coefficients of variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      for( v = 0; v < nvars; ++v )
         vals[v] = (SCIP_Real) weights[v];

      /* handle linear constraint for quadratic constraint update */
      SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
            vars, vals, lhs, rhs, nvars, varhash, dualconss, ndualconss, naddconss) );

      /* free buffer array */
      SCIPfreeBufferArray(scip, &vals);
   }

   /* remove knapsack constraints, since they are now redundant; their feasibility is already expressed
    * by s >= 0, where s is the new slack variable that we introduced for these linear constraints */
   for( c = nconss-1; c >= 0; --c )
   {
      assert( conss[c] != NULL );
      SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      ++(*ndelconss);
   }

   return SCIP_OKAY;
}

/** handle set packing constraints for quadratic constraint update, see the documentation of the function
 *  presolveAddKKTLinearCons() for an explanation
 */
static
SCIP_RETCODE presolveAddKKTSetppcConss(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss,          /**< buffer to increase with number of created additional constraints */
   int*                  ndelconss           /**< buffer to increase with number of deleted constraints */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );
   assert( ndelconss != NULL );

   conshdlr = SCIPfindConshdlr(scip, "setppc");
   if( conshdlr == NULL )
      return SCIP_OKAY;

   nconss = SCIPconshdlrGetNConss(conshdlr);
   conss = SCIPconshdlrGetConss(conshdlr);

   /* loop through linear constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_SETPPCTYPE type;
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nvars;
      int v;

      /* get data of constraint */
      cons = conss[c];
      assert( cons != NULL );

      /* get setppc type */
      type = SCIPgetTypeSetppc(scip, cons);
      lhs = -SCIPinfinity(scip);
      rhs = SCIPinfinity(scip);
      switch( type )
      {
         case SCIP_SETPPCTYPE_PARTITIONING:
            lhs = 1.0;
            rhs = 1.0;
            break;
         case SCIP_SETPPCTYPE_PACKING:
            rhs = 1.0;
            break;
         case SCIP_SETPPCTYPE_COVERING:
            lhs = 1.0;
            break;
         default:
            SCIPerrorMessage("unknown setppc type\n");
            return SCIP_INVALIDDATA;
      }

      nvars = SCIPgetNVarsSetppc(scip, cons);
      vars = SCIPgetVarsSetppc(scip, cons);

      /* set coefficients of variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      for( v = 0; v < nvars; ++v )
         vals[v] = 1.0;

      /* handle linear constraint for quadratic constraint update */
      SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
            vars, vals, lhs, rhs, nvars, varhash, dualconss, ndualconss, naddconss) );

      /* free buffer array */
      SCIPfreeBufferArray(scip, &vals);
   }

   /* remove set packing constraints if lhs != rhs, since they are now redundant; their feasibility is already expressed
    * by s >= 0, where s is the new slack variable that we introduced for these linear constraints */
   for( c = nconss-1; c >= 0; --c )
   {
      assert( conss[c] != NULL );

      if( SCIPgetTypeSetppc(scip, conss[c]) != SCIP_SETPPCTYPE_PARTITIONING )
      {
         assert( SCIPgetTypeSetppc(scip, conss[c]) == SCIP_SETPPCTYPE_PACKING
              || SCIPgetTypeSetppc(scip, conss[c]) == SCIP_SETPPCTYPE_COVERING );

         SCIP_CALL( SCIPdelCons(scip, conss[c]) );
         ++(*ndelconss);
      }
   }

   return SCIP_OKAY;
}

/** handle varbound constraints for quadratic constraint update, see the documentation of the function
 *  presolveAddKKTLinearCons() for an explanation
 */
static
SCIP_RETCODE presolveAddKKTVarboundConss(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss,          /**< buffer to increase with number of created additional constraints */
   int*                  ndelconss           /**< buffer to increase with number of deleted constraints */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );
   assert( ndelconss != NULL );

   conshdlr = SCIPfindConshdlr(scip, "varbound");
   if( conshdlr == NULL )
      return SCIP_OKAY;

   nconss = SCIPconshdlrGetNConss(conshdlr);
   conss = SCIPconshdlrGetConss(conshdlr);

   /* loop through linear constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nvars;

      /* allocate buffer arrays */
     SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
     SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

      /* get data of constraint */
      cons = conss[c];
      assert( cons != NULL );

      lhs = SCIPgetLhsVarbound(scip, cons);
      rhs = SCIPgetRhsVarbound(scip, cons);
      vars[0] = SCIPgetVarVarbound(scip, cons);
      vars[1] = SCIPgetVbdvarVarbound(scip, cons);
      vals[0] = 1.0;
      vals[1] = SCIPgetVbdcoefVarbound(scip, cons);
      nvars = 2;

      /* handle linear constraint for quadratic constraint update */
      SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
            vars, vals, lhs, rhs, nvars, varhash, dualconss, ndualconss, naddconss) );

      /* free buffer array */
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* remove varbound constraints if lhs != rhs, since they are now redundant; their feasibility is already expressed
    * by s >= 0, where s is the new slack variable that we introduced for these linear constraints */
   for( c = nconss-1; c >= 0; --c )
   {
      SCIP_CONS* cons;

      cons = conss[c];
      assert( cons != NULL );

      if( ! SCIPisFeasEQ(scip, SCIPgetLhsVarbound(scip, cons), SCIPgetRhsVarbound(scip, cons)) )
      {
         SCIP_CALL( SCIPdelCons(scip, cons) );
         ++(*ndelconss);
      }
   }

   return SCIP_OKAY;
}

/** handle logicor constraints for quadratic constraint update, see the documentation of the function
 *  presolveAddKKTLinearCons() for an explanation
 */
static
SCIP_RETCODE presolveAddKKTLogicorConss(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss,          /**< buffer to increase with number of created additional constraints */
   int*                  ndelconss           /**< buffer to increase with number of deleted constraints */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   int c;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );
   assert( ndelconss != NULL );

   conshdlr = SCIPfindConshdlr(scip, "logicor");
   if( conshdlr == NULL )
      return SCIP_OKAY;

   nconss = SCIPconshdlrGetNConss(conshdlr);
   conss = SCIPconshdlrGetConss(conshdlr);

   /* loop through linear constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nvars;
      int v;

      /* get data of constraint */
      cons = conss[c];
      assert( cons != NULL );

      /* get setppc type */
      lhs = 1.0;
      rhs = SCIPinfinity(scip);

      nvars = SCIPgetNVarsLogicor(scip, cons);
      vars = SCIPgetVarsLogicor(scip, cons);

      /* set coefficients of variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      for( v = 0; v < nvars; ++v )
         vals[v] = 1.0;

      /* handle linear constraint for quadratic constraint update */
      SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPconsGetName(cons),
            vars, vals, lhs, rhs, nvars, varhash, dualconss, ndualconss, naddconss) );

      /* free buffer array */
      SCIPfreeBufferArray(scip, &vals);
   }

   /* remove logicor constraints, since they are now redundant; their feasibility is already expressed
    * by s >= 0, where s is the new slack variable that we introduced for these linear constraints */
   for( c = nconss-1; c >= 0; --c )
   {
      assert( conss[c] != NULL );

      SCIP_CALL( SCIPdelCons(scip, conss[c]) );
      ++(*ndelconss);
   }

   return SCIP_OKAY;
}

/** handle aggregated variables for quadratic constraint update @n
 *  we apply the function presolveAddKKTLinearCons() to the aggregation constraint, see the documentation of this
 *  function for further information
 */
static
SCIP_RETCODE presolveAddKKTAggregatedVars(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_VAR**            agrvars,            /**< aggregated variables */
   int                   nagrvars,           /**< number of aggregated variables */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   int v;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( agrvars != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );

   /* loop through variables */
   for( v = 0; v < nagrvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_VAR** vars = NULL;
      SCIP_Real* vals = NULL;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nvars;

      var = agrvars[v];

      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED )
      {
         SCIP_Real constant;

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

         /* get aggregation variable */
         constant = SCIPvarGetAggrConstant(var);
         vars[0] = SCIPvarGetAggrVar(var);
         vals[0] = SCIPvarGetAggrScalar(var);
         vars[1] = var;
         vals[1] = -1.0;
         lhs = -constant;
         rhs = -constant;
         nvars = 2;
      }
      else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      {
         SCIP_Real* scalars;
         SCIP_VAR** multvars;
         SCIP_Real constant;
         int nmultvars;
         int nbuffer;
         int j;

         nmultvars = SCIPvarGetMultaggrNVars(var);
         nbuffer = nmultvars+1;

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nbuffer) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, nbuffer) );

         /* get aggregation variables */
         multvars = SCIPvarGetMultaggrVars(var);
         scalars = SCIPvarGetMultaggrScalars(var);
         constant = SCIPvarGetMultaggrConstant(var);

         /* add multi-aggregated variables to array */
         for( j = 0; j < nmultvars; ++j )
         {
            vars[j] = multvars[j];
            vals[j] = scalars[j];
         }

         /* add new variable to array */
         vars[nmultvars] = var;
         vals[nmultvars] = -1.0;
         lhs = -constant;
         rhs = -constant;
         nvars = nmultvars + 1;
      }
      else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
      {
         SCIP_VAR* negvar;
         SCIP_Real negconst;

         /* get negation variable and negation offset */
         negvar = SCIPvarGetNegationVar(var);
         negconst = SCIPvarGetNegationConstant(var);

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

         vars[0] = negvar;
         vars[1] = var;
         vals[0] = 1.0;
         vals[1] = 1.0;
         lhs = negconst;
         rhs = negconst;
         nvars = 2;
      }
      else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);
         assert( SCIPisFeasEQ(scip, lb, ub) );

         if( SCIPisFeasZero(scip, lb) && SCIPisFeasZero(scip, ub) )
            continue;
         else
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &vars, 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vals, 1) );

            vars[0] = var;
            vals[0] = 1.0;
            lhs = lb;
            rhs = lb;
            nvars = 1;
         }
      }
      else
      {
         SCIPerrorMessage("unexpected variable status\n");
         return SCIP_ERROR;
      }

      if( nvars > 0 )
      {
         /* handle aggregation constraint for quadratic constraint update */
         SCIP_CALL( presolveAddKKTLinearCons(scip, objcons, SCIPvarGetName(var),
               vars, vals, lhs, rhs, nvars, varhash, dualconss, ndualconss, naddconss) );
      }

      SCIPfreeBufferArrayNull(scip, &vars);
      SCIPfreeBufferArrayNull(scip, &vals);
   }

   return SCIP_OKAY;
}

/** handle bilinear terms of quadratic constraint for quadratic constraint update
 *
 * For the two variables of each bilinear term
 * 1. create the dual constraints (i.e., the two rows of \f$Q x + c + A^T \mu = 0\f$) associated to these variables, if not
 *    done already
 * 2. create the dual variables and the complementarity constraints for the lower and upper bound constraints of the two
 *    variables of the bilinear term, if not done already
 * 3. add the coefficient \f$Q_{ij}\f$ of the bilinear term to the dual constraint
 *
 * for steps 1 and 2 see the documentation of createKKTDualCons() for further information.
 **/
static
SCIP_RETCODE presolveAddKKTQuadBilinearTerms(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_BILINTERM*       bilinterms,         /**< bilinear terms of quadratic constraint */
   int                   nbilinterms,        /**< number of bilinear terms of quadratic constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_Real             scale,              /**< scale factor of quadratic constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   int j;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );

   /* return if there are no bilinear terms */
   if( bilinterms == NULL )
      return SCIP_OKAY;

   /* loop through bilinear terms of quadratic constraint */
   for( j = 0; j < nbilinterms; ++j )
   {
      int i;

      /* quadratic matrix has to be symmetric; therefore, split bilinear terms into two parts */
      for( i = 0; i < 2; ++i )
      {
         SCIP_CONS* dualcons = NULL;  /* dual constraint associated to variable */
         SCIP_VAR* bilvar1;
         SCIP_VAR* bilvar2;

         if( i == 0 )
         {
            bilvar1 = bilinterms[j].var1;
            bilvar2 = bilinterms[j].var2;
         }
         else
         {
            bilvar1 = bilinterms[j].var2;
            bilvar2 = bilinterms[j].var1;
         }

         /* create/get dual constraint associated to variable 'bilvar1';
          * if variable does not already exist in hashmap then create dual variables for its bounds */
         SCIP_CALL( createKKTDualCons(scip, objcons, bilvar1, varhash, dualconss, ndualconss, &dualcons, naddconss) );
         assert( dualcons != NULL );

         /* add variable to dual constraint */
         assert( ! SCIPisFeasZero(scip, scale) );
         SCIP_CALL( SCIPaddCoefLinear(scip, dualcons, bilvar2, bilinterms[j].coef / scale) );
      }
   }

   return SCIP_OKAY;
}

/** handle quadratic terms of quadratic constraint for quadratic constraint update
 *
 * For each quadratic term variable
 * 1. create the dual constraint (i.e., a row of \f$Q x + c + A^T \mu = 0\f$) associated to this variable, if not done
 *    already
 * 2. create the dual variables and the complementarity constraints for the lower and upper bound constraints of this
 *    variable, if not done already
 * 3. add the coefficient \f$Q_{ii}\f$ of this variable to the dual constraint
 *
 * for steps 1 and 2 see the documentation of createKKTDualCons() for further information.
 **/
static
SCIP_RETCODE presolveAddKKTQuadQuadraticTerms(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_QUADVARTERM*     quadterms,          /**< quadratic terms of quadratic constraint */
   int                   nquadterms,         /**< number of quadratic terms of quadratic constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_Real             scale,              /**< scale factor of quadratic constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   int j;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( varhash != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );

   /* return if there are no quadratic terms */
   if( quadterms == NULL )
      return SCIP_OKAY;

   /* loop through quadratic terms */
   for( j = 0; j < nquadterms; ++j )
   {
      SCIP_CONS* dualcons = NULL;  /* dual constraint associated to variable */
      SCIP_VAR* quadvar;

      quadvar = quadterms[j].var;

      /* create/get dual constraint associated to variable 'bilvar1';
       * if variable does not already exist in hashmap then create dual variables for its bounds */
      SCIP_CALL( createKKTDualCons(scip, objcons, quadvar, varhash, dualconss, ndualconss, &dualcons, naddconss) );
      assert( dualcons != NULL );

      /* add variable to dual constraint */
      assert( ! SCIPisFeasZero(scip, scale) );
      SCIP_CALL( SCIPaddCoefLinear(scip, dualcons, quadvar, (SCIP_Real) quadterms[j].sqrcoef * 2 / scale) );
   }

   return SCIP_OKAY;
}

/** handle linear terms of quadratic constraint for quadratic constraint update
 *
 * For each linear term variable
 * 1. create the dual constraint (i.e., a row of \f$Q x + c + A^T \mu = 0\f$) associated to this variable, if not done
 *    already
 * 2. create the dual variables and the complementarity constraints for the lower and upper bound constraints of this
 *    variable, if not done already
 * 3. add the right hand side \f$-c_i\f$ to the dual constraint
 * 4. add \f$c_i\f$ to the objective constraint \f$1/2 ( c^T x + b^T \mu) = t\f$, where t is the objective variable
 *
 * for steps 1 and 2 see the documentation of createKKTDualCons() for further information.
 **/
static
SCIP_RETCODE presolveAddKKTQuadLinearTerms(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            objcons,            /**< objective constraint */
   SCIP_VAR**            lintermvars,        /**< linear terms of quadratic constraint */
   SCIP_Real*            lintermcoefs,       /**< coefficients of linear terms of quadratic constraint */
   int                   nlintermvars,       /**< number of linear terms of quadratic constraints */
   SCIP_QUADVARTERM*     quadterms,          /**< quadratic terms of quadratic constraint */
   int                   nquadterms,         /**< number of quadratic terms of quadratic constraint */
   SCIP_HASHMAP*         varhash,            /**< hash map from variable to index of linear constraint */
   SCIP_VAR*             objvar,             /**< variable of objective function */
   SCIP_Real             scale,              /**< scale factor of quadratic constraint */
   SCIP_CONS**           dualconss,          /**< array with dual constraints */
   int*                  ndualconss,         /**< pointer to store number of dual constraints */
   int*                  naddconss           /**< buffer to increase with number of created additional constraints */
   )
{
   int j;

   assert( scip != NULL );
   assert( objcons != NULL );
   assert( lintermcoefs != NULL );
   assert( varhash != NULL );
   assert( objvar != NULL );
   assert( dualconss != NULL );
   assert( ndualconss != NULL );
   assert( naddconss != NULL );

   /* loop through linear terms of quadratic constraint */
   if( lintermvars != NULL )
   {
      assert( lintermcoefs != NULL );
      for( j = 0; j < nlintermvars; ++j )
      {
         SCIP_VAR* var;

         var = lintermvars[j];

         if( var != objvar )
         {
            SCIP_CONS* dualcons = NULL;  /* dual constraint associated to variable */
            SCIP_Real coef;

            /* create/get dual constraint associated to variable;
             * if variable does not already exist in hashmap then create dual variables for its bounds */
            SCIP_CALL( createKKTDualCons(scip, objcons, var, varhash, dualconss, ndualconss, &dualcons, naddconss) );
            assert( dualcons != NULL );

            /* get coefficient of variable in quadratic constraint */
            coef = lintermcoefs[j];

            /* change lhs and rhs of dual constraint */
            assert( ! SCIPisFeasZero(scip, scale) );
            SCIP_CALL( SCIPchgLhsLinear(scip, dualcons, SCIPgetLhsLinear(scip, dualcons) - coef / scale) );
            SCIP_CALL( SCIPchgRhsLinear(scip, dualcons, SCIPgetRhsLinear(scip, dualcons) - coef / scale) );

            /* add variable to objective constraint */
            SCIP_CALL( SCIPaddCoefLinear(scip, objcons, var, coef / (scale * 2)) );
         }
      }
   }

   /* loop through linear terms that are part of a quadratic term of quadratic constraint */
   if( quadterms != NULL )
   {
      for( j = 0; j < nquadterms; ++j )
      {
         SCIP_CONS* dualcons;
         SCIP_Real coef;
         SCIP_VAR* var;
         int ind;

         var = quadterms[j].var;
         coef = quadterms[j].lincoef;
         assert( var != objvar );

         /* get dual constraint associated to variable (has already been created in function
          * presolveAddKKTQuadQuadraticTerms() */
         assert( SCIPhashmapExists(varhash, var) );
         ind = (int) (size_t) SCIPhashmapGetImage(varhash, var);
         dualcons = dualconss[ind];
         assert( dualcons != NULL );

         /* change lhs and rhs of dual constraint */
         assert( ! SCIPisFeasZero(scip, scale) );
         SCIP_CALL( SCIPchgLhsLinear(scip, dualcons, SCIPgetLhsLinear(scip, dualcons) -coef / scale) );
         SCIP_CALL( SCIPchgRhsLinear(scip, dualcons, SCIPgetRhsLinear(scip, dualcons) -coef / scale) );

         /* add variable to objective constraint */
         SCIP_CALL( SCIPaddCoefLinear(scip, objcons, var, coef / (scale * 2)) );
      }
   }

   return SCIP_OKAY;
}

/** checks for a given constraint whether it is the objective function of a (mixed-binary) quadratic program
 * \f[
 *  \begin{array}{ll}
 *  \min         & z \\
 *  s.t.         & x^T Q x + c^T x + d <= z \\
 *               & A x \leq b, \\
 *               & x \in \{0, 1\}^{p} \times R^{n-p},
 * \end{array}
 * \f]
 * which is equivalent to
 * \f[
 *  \begin{array}{ll}
 *  \min         & x^T Q x + c^T x + d \\
 *  s.t.         & A x \leq b, \\
 *               & x \in \{0, 1\}^{p} \times R^{n-p}.
 * \end{array}
 * \f]
 *
 *
 * We check whether
 * 1. there is a single quadratic constraint that can be written as \f$x^T Q x + c^T x + d \leq t\f$
 * 2. all other constraints are linear
 * 3. all integer variables are binary if allowbinary = TRUE, or all variables are continuous if allowbinary = FALSE
 * 4. t is the only variable in the objective and doesn't appear in any other constraint
 */
static
SCIP_RETCODE checkConsQuadraticProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        quadconshdlr,       /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< quadratic constraint */
   SCIP_Bool             allowbinary,        /**< if TRUE then allow binary variables in the problem, if FALSE then all
                                              *   variables have to be continuous */
   SCIP_VAR**            objvar,             /**< pointer to store the objective variable @p z */
   SCIP_Real*            scale,              /**< pointer to store the value by which we have to scale the quadratic
                                              *   constraint such that the objective variable @p z has coefficient -1 */
   SCIP_Real*            objrhs,             /**< pointer to store the right hand side @p -d of the objective constraint */
   SCIP_Bool*            isqp                /**< pointer to store whether the problem is a (mixed-binary) QP */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_VAR** lintermvars;
   SCIP_Real* lintermcoefs;
   int nconss = 0;
   SCIP_Real quadlhs;
   SCIP_Real quadrhs;
   SCIP_Real coef;
   SCIP_Real obj;
   int mayincrease;
   int maydecrease;
   int objind = -1;

   *objrhs = 0.0;
   *scale = 0.0;
   *isqp = FALSE;
   *objvar = NULL;

   /* desired structure: there exists only one variable with nonzero objective value; this is the objective variable 'z' */
   if( SCIPgetNObjVars(scip) != 1 )
      return SCIP_OKAY;

   /* desired structure: all integer variables are binary; if the parameter 'allowbinary' is set to FALSE, then all
    * variables have to be continuous */
   if( SCIPgetNIntVars(scip) > 0 || ( ! allowbinary && SCIPgetNBinVars(scip) > 0 ) )
      return SCIP_OKAY;

   /* desired structure: there exists only one quadratic constraint */
   if( SCIPconshdlrGetNConss(quadconshdlr) != 1 )
      return SCIP_OKAY;

   /* desired structure: the constraint has to take one of the three forms
    * i)   x^T Q x + c^T x <= d
    * ii)  x^T Q x + c^T x >= d
    * iii) x^T Q x + c^T x == d
    * the case a <= x^T Q x + c^T x <= d with 'a' and 'd' finite and a != d is not allowed.
    */
   quadlhs = SCIPgetLhsQuadratic(scip, cons);
   quadrhs = SCIPgetRhsQuadratic(scip, cons);
   if( ! SCIPisFeasEQ(scip, quadlhs, quadrhs) && ! SCIPisInfinity(scip, -quadlhs) && ! SCIPisInfinity(scip, quadrhs) )
      return SCIP_OKAY;

   /* get number of linear constraints (including special cases of linear constraints) */
   conshdlr = SCIPfindConshdlr(scip, "linear");
   if( conshdlr != NULL )
      nconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "setppc");
   if( conshdlr != NULL )
      nconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "knapsack");
   if( conshdlr != NULL )
      nconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "varbound");
   if( conshdlr != NULL )
      nconss += SCIPconshdlrGetNConss(conshdlr);

   conshdlr = SCIPfindConshdlr(scip, "logicor");
   if( conshdlr != NULL )
      nconss += SCIPconshdlrGetNConss(conshdlr);

   /* desired structure: all the nonquadratic constraints are linear constraints */
   if( nconss != SCIPgetNConss(scip) - 1 )
      return SCIP_OKAY;

   /* get variables that are in the linear term of the quadratic constraint */
   lintermvars = SCIPgetLinearVarsQuadratic(scip, cons);
   lintermcoefs = SCIPgetCoefsLinearVarsQuadratic(scip, cons);

   /* compute the objective shift of the QP. Note that
    *
    * min z     s.t.    x^T Q x + c^T x <= d + z
    *                   Ax <= b
    *
    * is equivalent to
    *
    * min x^T Q x + c^T x - d    s.t.    Ax <= b
    *
    * Here, -d is the objective shift. We define b to be the right hand side of the objective constraint.
    */
   if( ! SCIPisInfinity(scip, -quadlhs) )
      *objrhs = quadlhs;
   else
      *objrhs = quadrhs;
   assert( ! SCIPisInfinity(scip, REALABS(*objrhs)) );

   /* search for the objective variable 'objvar' in the linear term of quadratic constraint (it is already known that
    * at most one variable has a nonzero objective value); additionally, check the sign of the objective variable */
   maydecrease = SCIPgetLinvarMayDecreaseQuadratic(scip, cons);
   mayincrease = SCIPgetLinvarMayIncreaseQuadratic(scip, cons);
   if( maydecrease < 0 && mayincrease < 0 )
      return SCIP_OKAY;
   else if( maydecrease >= 0 )
   {
      objind = maydecrease;

      /* if both mayincrease and maydecrese are nonnegative, then check objective coefficient */
      if( mayincrease >= 0 && SCIPisFeasZero(scip, SCIPvarGetObj(lintermvars[maydecrease])) )
         objind = mayincrease;
   }
   else
      objind = mayincrease;
   assert( objind < SCIPgetNLinearVarsQuadratic(scip, cons) );

   *objvar = lintermvars[objind];
   coef = lintermcoefs[objind];
   obj = SCIPvarGetObj(*objvar);

   /* check sign of coefficient */
   if( SCIPisFeasPositive(scip, obj)
          && ( ( SCIPisFeasNegative(scip, coef) && SCIPisFeasEQ(scip, quadrhs, *objrhs) )
               || ( SCIPisFeasPositive(scip, coef) && SCIPisFeasEQ(scip, quadlhs, *objrhs) )
             )
      )
      *scale = -1.0/coef; /* value by which we have to scale the quadratic constraint such that the objective variable
                           * has coefficient -1 */
   else if( SCIPisFeasNegative(scip, obj)
             && ( ( SCIPisFeasNegative(scip, coef) && SCIPisFeasEQ(scip, quadlhs, *objrhs) )
                  || ( SCIPisFeasPositive(scip, coef) && SCIPisFeasEQ(scip, quadrhs, *objrhs) )
                )
           )
      *scale = 1.0/coef; /* value by which we have to scale the quadratic constraint such that the objective variable
                          * has coefficient 1 */
   else
      return SCIP_OKAY;
   assert( *objvar != NULL && ! SCIPisFeasZero(scip, SCIPvarGetObj(*objvar)) );
   assert( ! SCIPisFeasZero(scip, *scale) );

   /* scale the right hand side of the objective constraint */
   *objrhs = (*objrhs)/(*scale); /*lint !e414*/

   /* check whether 'objvar' is part of a linear constraint; if this is true then return
    * whether 'objvar' is part of a linear constraint can be deduced from the variable locks */
   if( SCIPisFeasEQ(scip, quadlhs, quadrhs) )
   {
      if( SCIPvarGetNLocksDown(*objvar) != 1 || SCIPvarGetNLocksUp(*objvar) != 1 )
         return SCIP_OKAY;
   }
   else
   {
      assert( SCIPisInfinity(scip, -quadlhs) || SCIPisInfinity(scip, quadrhs) );

      if( ( SCIPvarGetNLocksDown(*objvar) != 1 || SCIPvarGetNLocksUp(*objvar) != 0 )
           && ( SCIPvarGetNLocksDown(*objvar) != 0 || SCIPvarGetNLocksUp(*objvar) != 1 ) )
         return SCIP_OKAY;
   }

   *isqp = TRUE;

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyQPKKTref)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolQPKKTref(scip) );

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeQPKKTref)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecQPKKTref)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_CONSHDLR* linconshdlr;
   SCIP_CONSHDLR* quadconshdlr;
   SCIP_CONS** conss;
   SCIP_CONS* cons;

   SCIP_QUADVARTERM* quadterms;
   SCIP_BILINTERM* bilinterms;
   int nquadterms;
   int nbilinterms;

   SCIP_VAR** lintermvars;
   SCIP_Real* lintermcoefs;
   int nlintermvars;

   SCIP_CONS** savelinconss = NULL;
   SCIP_CONS** linconss = NULL;
   int nlinconss = 0;

   SCIP_HASHMAP* varhash; /* hash map from variable to index of dual constraint */
   SCIP_CONS** dualconss; /* constraints associated to the Lagrangean function */
   int ndualconss = 0;

   SCIP_CONS* objcons;
   SCIP_VAR* objvar;
   SCIP_Real scale;
   SCIP_Real objrhs;
   SCIP_Bool isqp;
   int j;

   assert( scip != NULL );
   assert( naddconss != NULL );
   assert( ndelconss != NULL );

   /* desired structure: there exists only one quadratic constraint */
   quadconshdlr = SCIPfindConshdlr(scip,"quadratic");
   if( quadconshdlr == NULL || SCIPconshdlrGetNConss(quadconshdlr) != 1 )
      return SCIP_OKAY;

   /* get quadratic constraint */
   conss = SCIPconshdlrGetConss(quadconshdlr);
   cons = conss[0];
   assert( cons != NULL );

   SCIPdebugMsg(scip, "tries to add the KKT conditions for constraint <%s>\n", SCIPconsGetName(cons));

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* desired structure: matrix associated to quadratic constraint is indefinite;
    * otherwise, the problem usually can be solved faster by standard methods. */
   SCIP_CALL( SCIPcheckCurvatureQuadratic(scip, cons) );
   if( ! presoldata->updatequadindef && ( SCIPisConvexQuadratic(scip, cons) || SCIPisConcaveQuadratic(scip, cons) ) )
   {
      SCIPdebugMsg(scip, "quadratic constraint update failed, since matrix associated to quadratic constraint <%s> is not \
           indefinite.\n", SCIPconsGetName(cons) );
      return SCIP_OKAY;
   }

   /* first, check whether the problem is equivalent to
    *
    * min   z
    * s.t.  x^T Q x + c^T x <= b + z
    *       x \in \{0, 1\}^{p} \times R^{n-p}.
    *
    */
   SCIP_CALL( checkConsQuadraticProblem(scip, quadconshdlr, cons, presoldata->addkktbinary, &objvar, &scale, &objrhs, &isqp) );
   if( ! isqp )
      return SCIP_OKAY;
   assert( objvar != NULL );

   /* get constraint handler data of linear constraints */
   linconshdlr = SCIPfindConshdlr(scip, "linear");

   /* get linear constraints and number of linear constraints */
   if( linconshdlr != NULL )
   {
      nlinconss = SCIPconshdlrGetNConss(linconshdlr);
      linconss = SCIPconshdlrGetConss(linconshdlr);
   }

   /* get variables that are in the linear term of the quadratic constraint */
   nlintermvars = SCIPgetNLinearVarsQuadratic(scip, cons);
   lintermvars = SCIPgetLinearVarsQuadratic(scip, cons);
   lintermcoefs = SCIPgetCoefsLinearVarsQuadratic(scip, cons);

   /* get bilinear terms */
   bilinterms = SCIPgetBilinTermsQuadratic(scip, cons);
   nbilinterms = SCIPgetNBilinTermsQuadratic(scip, cons);

   /* get quadratic terms */
   quadterms = SCIPgetQuadVarTermsQuadratic(scip, cons);
   nquadterms = SCIPgetNQuadVarTermsQuadratic(scip, cons);

   /* the update is only valid if a finite optimal solution of the problem exists,
    * since only finite optimal solutions satisfy the KKT conditions;
    * we check whether all variables have finite bounds, otherwise we return */
   if( presoldata->updatequadbounded )
   {
      /* check linear term variables */
      for( j = 0; j < nlintermvars; ++j )
      {
         SCIP_VAR* var;

         var = lintermvars[j];
         if( var != objvar &&
              ( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
            )
         {
            SCIPdebugMsg(scip, "failed adding the KKT conditions, since not all variables to quadratic constraint <%s> are \
                 bounded.\n", SCIPconsGetName(cons) );
            return SCIP_OKAY;
         }
      }

      /* check linear term variables */
      for( j = 0; j < nbilinterms; ++j )
      {
         SCIP_VAR* bilvar1;
         SCIP_VAR* bilvar2;

         bilvar1 = bilinterms[j].var1;
         bilvar2 = bilinterms[j].var2;
         if( ( bilvar1 != objvar && ( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(bilvar1)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(bilvar1)) ) )
            || ( bilvar2 != objvar && ( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(bilvar2)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(bilvar2)) ) ) )
         {
            SCIPdebugMsg(scip, "failed adding the KKT conditions, since not all variables to quadratic constraint <%s> \
                 are bounded.\n", SCIPconsGetName(cons) );
            return SCIP_OKAY;
         }
      }

      /* check quadratic term variables */
      for( j = 0; j < nquadterms; ++j )
      {
         SCIP_VAR* var;

         var = quadterms[j].var;
         if( var != objvar
              && ( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
            )
         {
            SCIPdebugMsg(scip, "failed adding the KKT conditions, since not all variables to quadratic constraint <%s> \
                 are bounded.\n", SCIPconsGetName(cons) );
            return SCIP_OKAY;
         }
      }
   }


   /* add KKT constraints */

   /* set up hash map */
   SCIP_CALL( SCIPhashmapCreate(&varhash, SCIPblkmem(scip), SCIPgetNVars(scip) + SCIPgetNFixedVars(scip)) );

   /* allocate buffer array */
   SCIP_CALL( SCIPallocBufferArray(scip, &dualconss, 2 * SCIPgetNVars(scip) + 2 * SCIPgetNFixedVars(scip)) ); /*lint !e647*/

   /* duplicate linconss for later use, since in the following, we create new linear constraints */
   if( linconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &savelinconss, linconss, nlinconss) );
   }

   /* create new objective constraint */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &objcons, "objcons", 0, NULL, NULL, objrhs, objrhs) );
   if( SCIPisFeasNegative(scip, SCIPvarGetObj(objvar)) )
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, objcons, objvar, 1.0) );
   }
   else
   {
      SCIP_CALL( SCIPaddCoefLinear(scip, objcons, objvar, -1.0) );
   }

   /* handle linear constraints */
   if( savelinconss != NULL )
   {
      SCIP_CALL( presolveAddKKTLinearConss(scip, objcons, savelinconss, nlinconss, varhash, dualconss, &ndualconss,
            naddconss, ndelconss) );
   }

   /* handle set packing constraints */
   SCIP_CALL( presolveAddKKTSetppcConss(scip, objcons, varhash, dualconss, &ndualconss, naddconss, ndelconss) );

   /* handle knapsack constraints */
   SCIP_CALL( presolveAddKKTKnapsackConss(scip, objcons, varhash, dualconss, &ndualconss, naddconss, ndelconss) );

   /* handle varbound constraints */
   SCIP_CALL( presolveAddKKTVarboundConss(scip, objcons, varhash, dualconss, &ndualconss, naddconss, ndelconss) );

   /* handle logicor constraints */
   SCIP_CALL( presolveAddKKTLogicorConss(scip, objcons, varhash, dualconss, &ndualconss, naddconss, ndelconss) );

   /* handle linear constraints associated to aggregations of variables */
   if( SCIPgetNFixedVars(scip) > 0 )
   {
      SCIP_CALL( presolveAddKKTAggregatedVars(scip, objcons, SCIPgetFixedVars(scip), SCIPgetNFixedVars(scip),
            varhash, dualconss, &ndualconss, naddconss) );
   }

   /* handle bilinear terms of quadratic constraint */
   SCIP_CALL( presolveAddKKTQuadBilinearTerms(scip, objcons, bilinterms, nbilinterms, varhash, scale, dualconss,
        &ndualconss, naddconss) );

   /* handle quadratic terms of quadratic constraint */
   SCIP_CALL( presolveAddKKTQuadQuadraticTerms(scip, objcons, quadterms, nquadterms, varhash, scale, dualconss,
        &ndualconss, naddconss) );

   /* handle linear terms of quadratic constraint */
   SCIP_CALL( presolveAddKKTQuadLinearTerms(scip, objcons, lintermvars, lintermcoefs, nlintermvars, quadterms, nquadterms,
        varhash, objvar, scale, dualconss, &ndualconss, naddconss) );


   /* add/release objective constraint */
   SCIP_CALL( SCIPaddCons(scip, objcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &objcons) );
   ++(*naddconss);

   /* add/release dual constraints associated to the KKT conditions */
   for( j = 0; j < ndualconss; ++j )
   {
      SCIP_CALL( SCIPaddCons(scip, dualconss[j]) );
      SCIP_CALL( SCIPreleaseCons(scip, &dualconss[j]) );
   }
   *naddconss = *naddconss + ndualconss;

   /* free buffer array */
   SCIPfreeBufferArrayNull(scip, &savelinconss);
   SCIPfreeBufferArray(scip, &dualconss);

   /* free hash map */
   SCIPhashmapFree(&varhash);

   if( SCIPgetNBinVars(scip) > 0 )
      SCIPdebugMsg(scip, "added the KKT conditions to the mixed-binary quadratic program\n");
   else
      SCIPdebugMsg(scip, "added the KKT conditions to the quadratic program\n");

   /*SCIP_CALL( SCIPwriteTransProblem(scip, "trafoQP.lp", NULL, FALSE ) );*/

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the QP KKT reformulation presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolQPKKTref(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol= NULL;

   /* alloc presolve data object */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecQPKKTref, presoldata) );
   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyQPKKTref) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeQPKKTref) );

   /* add qpkktref presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "presolving/" PRESOL_NAME "/addkktbinary",
         "if TRUE then allow binary variables for KKT update",
         &presoldata->addkktbinary, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "presolving/" PRESOL_NAME "/updatequadbounded",
         "if TRUE then only apply the update to QPs with bounded variables; if the variables are not bounded then a \
         finite optimal solution might not exist and the KKT conditions would then be invalid",
         &presoldata->updatequadbounded, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "presolving/" PRESOL_NAME "/updatequadindef",
         "if TRUE then apply quadratic constraint update even if the quadratic constraint matrix is known to be indefinite",
         &presoldata->updatequadindef, TRUE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
