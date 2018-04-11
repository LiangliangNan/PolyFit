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

/**@file   cons_soc.c
 * @brief  constraint handler for second order cone constraints \f$\sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} \leq \alpha_{n+1}\, (x_{n+1}+\beta_{n+1})\f$
 * @author Stefan Vigerske
 * @author Marc Pfetsch
 * @author Felipe Serrano
 * @author Benjamin Mueller
 *
 * @todo rhsvar == NULL is supported in some routines, but not everywhere
 * @todo merge square terms with same variables in presol/exitpre
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <assert.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define SCIP_PRIVATE_ROWPREP

#include "scip/cons_soc.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/intervalarith.h"
#include "nlpi/nlpi.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi_ipopt.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "soc"
#define CONSHDLR_DESC          "constraint handler for second order cone constraints"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -40 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING      SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */
#define CONSHDLR_PRESOLTIMING     SCIP_PRESOLTIMING_ALWAYS /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */

#define QUADCONSUPGD_PRIORITY     10000 /**< priority of the constraint handler for upgrading of quadratic constraints */

#define UPGSCALE 10 /* scale factor used in general upgrades of quadratic cons to soc */

/*
 * Data structures
 */

/** Eventdata for variable bound change events. */
struct VarEventData
{
   SCIP_CONS*            cons;               /**< the constraint */
   int                   varidx;             /**< the index of a variable on the left hand side which bound change is caught, or -1 for variable on right hand side */
   int                   filterpos;          /**< position of corresponding event in event filter */
};
typedef struct VarEventData VAREVENTDATA;

/** constraint data for soc constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables on left hand side (n) */
   SCIP_VAR**            vars;               /**< variables on left hand side (x_i) */
   SCIP_Real*            coefs;              /**< coefficients for variables on left hand side (alpha_i) */
   SCIP_Real*            offsets;            /**< offsets for variables on left hand side (beta_i) */
   SCIP_Real             constant;           /**< constant on left hand side (gamma) */

   SCIP_VAR*             rhsvar;             /**< variable on right hand side (x_{n+1}) */
   SCIP_Real             rhscoeff;           /**< coefficient of square term on right hand side (alpha_{n+1}) */
   SCIP_Real             rhsoffset;          /**< offset for variable on right hand side (beta_{n+1}) */

   SCIP_NLROW*           nlrow;              /**< nonlinear row representation of constraint */

   SCIP_Real             lhsval;             /**< value of left hand side in current point */
   SCIP_Real             violation;          /**< violation of constraint in current point */

   VAREVENTDATA*         lhsbndchgeventdata; /**< eventdata for bound change events on left  hand side variables */
   VAREVENTDATA          rhsbndchgeventdata; /**< eventdata for bound change event  on right hand side variable  */
   SCIP_Bool             isapproxadded;      /**< has a linear outer approximation be added? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subNLP heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the trysol heuristic, if available */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   int                   newsoleventfilterpos;/**< filter position of new solution event handler, if caught */
   SCIP_Bool             haveexprint;        /**< indicates whether an expression interpreter is available */
   SCIP_Bool             sepanlp;            /**< where linearization of the NLP relaxation solution added? */

   SCIP_Bool             glineur;            /**< is the Glineur outer approx preferred to Ben-Tal Nemirovski? */
   SCIP_Bool             projectpoint;       /**< is the point in which a cut is generated projected onto the feasible set? */
   int                   nauxvars;           /**< number of auxiliary variables to use when creating a linear outer approx. of a SOC3 constraint */
   SCIP_Bool             sparsify;           /**< whether to sparsify cuts */
   SCIP_Real             sparsifymaxloss;    /**< maximal loss in cut efficacy by sparsification */
   SCIP_Real             sparsifynzgrowth;   /**< growth rate of maximal allowed nonzeros in cuts in sparsification */
   SCIP_Bool             linfeasshift;       /**< whether to try to make solutions feasible in check by shifting the variable on the right hand side */
   char                  nlpform;            /**< formulation of SOC constraint in NLP */
   SCIP_Real             sepanlpmincont;     /**< minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation */
   SCIP_Bool             enfocutsremovable;  /**< are cuts added during enforcement removable from the LP in the same node? */
   SCIP_Bool             generalsocupg;      /**< try to upgrade more general quadratics to soc? */
   SCIP_Bool             disaggregate;       /**< try to completely disaggregate soc? */

   SCIP_NODE*            lastenfonode;       /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenforounds;        /**< counter on number of enforcement rounds for the current node */
};


/*
 * Local methods
 */

/** catch left hand side variable events */
static
SCIP_RETCODE catchLhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   varidx              /**< index of the variable which events to catch */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(varidx >= 0);
   assert(varidx < consdata->nvars);
   assert(consdata->lhsbndchgeventdata != NULL);

   consdata->lhsbndchgeventdata[varidx].cons = cons;
   consdata->lhsbndchgeventdata[varidx].varidx   = varidx;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[varidx], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
         (SCIP_EVENTDATA*)&consdata->lhsbndchgeventdata[varidx], &consdata->lhsbndchgeventdata[varidx].filterpos) );

   /* since bound changes were not catched before, a possibly stored activity may have become outdated */
   SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );

   return SCIP_OKAY;
}

/** catch right hand side variable events */
static
SCIP_RETCODE catchRhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);

   consdata->rhsbndchgeventdata.cons = cons;
   consdata->rhsbndchgeventdata.varidx   = -1;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->rhsvar, SCIP_EVENTTYPE_UBTIGHTENED, eventhdlr,
         (SCIP_EVENTDATA*)&consdata->rhsbndchgeventdata, &consdata->rhsbndchgeventdata.filterpos) );

   /* since bound changes were not catched before, a possibly stored activity may have become outdated */
   SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );

   return SCIP_OKAY;
}

/** catch variables events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int            i;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(consdata->lhsbndchgeventdata == NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->lhsbndchgeventdata, consdata->nvars) );

   for( i = 0; i < consdata->nvars; ++i )
   {
      if( consdata->vars[i] != NULL )
      {
         SCIP_CALL( catchLhsVarEvents(scip, eventhdlr, cons, i) );
      }
   }

   if( consdata->rhsvar != NULL )
   {
      SCIP_CALL( catchRhsVarEvents(scip, eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** drop left hand side variable events */
static
SCIP_RETCODE dropLhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   varidx              /**< index of the variable which events to catch */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(varidx >= 0);
   assert(varidx < consdata->nvars);
   assert(consdata->lhsbndchgeventdata != NULL);
   assert(consdata->lhsbndchgeventdata[varidx].varidx == varidx);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[varidx], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr,
         (SCIP_EVENTDATA*)&consdata->lhsbndchgeventdata[varidx], consdata->lhsbndchgeventdata[varidx].filterpos) );

   return SCIP_OKAY;
}

/** drop right hand side variable events */
static
SCIP_RETCODE dropRhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(consdata->rhsbndchgeventdata.varidx == -1);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->rhsvar, SCIP_EVENTTYPE_UBTIGHTENED, eventhdlr,
         (SCIP_EVENTDATA*)&consdata->rhsbndchgeventdata, consdata->rhsbndchgeventdata.filterpos) );

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip      != NULL);
   assert(eventhdlr != NULL);
   assert(cons      != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);

   for( i = 0; i < consdata->nvars; ++i )
   {
      if( consdata->vars[i] != NULL )
      {
         SCIP_CALL( dropLhsVarEvents(scip, eventhdlr, cons, i) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &consdata->lhsbndchgeventdata, consdata->nvars);

   if( consdata->rhsvar != NULL )
   {
      SCIP_CALL( dropRhsVarEvents(scip, eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** process variable bound tightening event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_CONS* cons;

   assert(scip      != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   cons = ((VAREVENTDATA*)eventdata)->cons;
   assert(cons != NULL);

   SCIP_CALL( SCIPmarkConsPropagate(scip, cons) );
   /* @todo look at bounds on x_i to decide whether propagation makes sense */

   return SCIP_OKAY;
}

/** create a nonlinear row representation of the constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< SOC constraint handler */
   SCIP_CONS*            cons                /**< SOC constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   char nlpform;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   nlpform = conshdlrdata->nlpform;
   if( nlpform == 'a' )
   {
      /* if the user let us choose, then we take 's' for "small" SOC constraints, but 'q' for large ones,
       * since the 's' form leads to nvars^2 elements in Hessian, while the 'q' form yields only n elements
       * however, if there is no expression interpreter, then the NLPI may have trouble, so we always use 'q' in this case
       */
      if( consdata->nvars < 100 && conshdlrdata->haveexprint )
         nlpform = 's';
      else
         nlpform = 'q';
   }

   switch( nlpform )
   {
   case 'e':
   {
      /* construct expression exp(\sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} - alpha_{n+1}(x_{n+1} + beta_{n+1})) */

      if( consdata->nvars > 0 )
      {
         SCIP_EXPR* expr;
         SCIP_EXPR* exprterm;
         SCIP_EXPR* expr2;
         SCIP_EXPRTREE* exprtree;

         if( consdata->constant != 0.0 )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_CONST, consdata->constant) );  /* gamma */
         }
         else
         {
            exprterm = NULL;
         }

         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, i) );  /* x_i */
            if( consdata->offsets[i] != 0.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->offsets[i]) );  /* beta_i */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_PLUS, expr, expr2) );  /* x_i + beta_i */
            }
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_SQUARE, expr) );  /* (x_i + beta_i)^2 */
            if( consdata->coefs[i] != 1.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, SQR(consdata->coefs[i])) );  /* (alpha_i)^2 */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL, expr, expr2) );  /* (alpha_i)^2 * (x_i + beta_i)^2 */
            }
            if( exprterm != NULL )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_PLUS, exprterm, expr) );
            }
            else
            {
               exprterm = expr;
            }
         }

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_SQRT, exprterm) );  /* sqrt(gamma + sum_i (...)^2) */

         if( consdata->rhsvar != NULL )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, consdata->nvars) );  /* x_{n+1} */
            if( consdata->rhsoffset != 0.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->rhsoffset) );  /* beta_{n+1} */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_PLUS, expr, expr2) );  /* x_{n+1} + beta_{n+1} */
            }
            if( consdata->rhscoeff != 1.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->rhscoeff) );  /* alpha_{n+1} */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL, expr, expr2) );  /* alpha_{n+1} * (x_{n+1} + beta_{n+1}) */
            }
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_CONST, consdata->rhscoeff * consdata->rhsoffset) );
         }
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_MINUS, exprterm, expr) ); /* sqrt(gamma + sum_i (...)^2) - alpha_{n+1} * (x_{n+1} + beta_{n+1}) */

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_EXP, exprterm) ); /* exp(sqrt(gamma + sum_i (...)^2) - alpha_{n+1} * (x_{n+1} + beta_{n+1})) */

         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exprterm, consdata->nvars+1, 0, NULL) );

         SCIP_CALL( SCIPexprtreeSetVars(exprtree, consdata->nvars, consdata->vars) );
         SCIP_CALL( SCIPexprtreeAddVars(exprtree, 1, &consdata->rhsvar) );

         SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons),
               0.0,
               0, NULL, NULL,
               0, NULL, 0, NULL,
               exprtree, -SCIPinfinity(scip), 1.0,
               SCIP_EXPRCURV_CONVEX) );

         SCIP_CALL( SCIPexprtreeFree(&exprtree) );

         break;
      }
      /* if there are no left-hand-side variables, then we let the 's' case handle it */
   } /*lint -fallthrough */

   case 's':
   {
      /* construct expression \sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} */

      SCIP_EXPR* expr;
      SCIP_EXPR* exprterm;
      SCIP_EXPR* expr2;
      SCIP_EXPRTREE* exprtree;
      SCIP_Real lincoef;

      if( consdata->constant != 0.0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_CONST, consdata->constant) );  /* gamma */
      }
      else
      {
         exprterm = NULL;
      }

      for( i = 0; i < consdata->nvars; ++i )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, i) );  /* x_i */
         if( consdata->offsets[i] != 0.0 )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->offsets[i]) );  /* beta_i */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_PLUS, expr, expr2) );  /* x_i + beta_i */
         }
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_SQUARE, expr) );  /* (x_i + beta_i)^2 */
         if( consdata->coefs[i] != 1.0 )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, SQR(consdata->coefs[i])) );  /* (alpha_i)^2 */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL, expr, expr2) );  /* (alpha_i)^2 * (x_i + beta_i)^2 */
         }
         if( exprterm != NULL )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_PLUS, exprterm, expr) );
         }
         else
         {
            exprterm = expr;
         }
      }

      if( exprterm != NULL )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_SQRT, exprterm) );  /* sqrt(gamma + sum_i (...)^2) */
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exprterm, consdata->nvars, 0, NULL) );
         SCIP_CALL( SCIPexprtreeSetVars(exprtree, consdata->nvars, consdata->vars) );
      }
      else
      {
         assert(consdata->nvars == 0);
         assert(consdata->constant == 0.0);
         exprtree = NULL;
      }

      /* linear and constant part is -\alpha_{n+1} (x_{n+1}+\beta_{n+1}) */
      lincoef = -consdata->rhscoeff;
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons),
            -consdata->rhscoeff * consdata->rhsoffset,
            1, &consdata->rhsvar, &lincoef,
            0, NULL, 0, NULL,
            exprtree, -SCIPinfinity(scip), 0.0,
            SCIP_EXPRCURV_CONVEX) );

      SCIP_CALL( SCIPexprtreeFree(&exprtree) );

      break;
   }

   case 'q':
   {
      /* construct quadratic form gamma + sum_{i=1}^{n} (alpha_i (x_i + beta_i))^2 <= (alpha_{n+1} (x_{n+1} + beta_{n+1})^2 */
      SCIP_QUADELEM sqrterm;
      SCIP_EXPRCURV curvature;
      SCIP_Real rhs;
      int rhsvarpos;

      /* the expression is convex if alpha_{n+1} is zero */
      curvature = consdata->rhscoeff == 0.0 ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_UNKNOWN;

      /* create initial empty row with left hand side variables */
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            0, NULL, NULL,
            consdata->nvars, consdata->vars, 0, NULL,
            NULL, -SCIPinfinity(scip), 0.0,
            curvature) );

      /* add gamma + sum_{i=1}^{n} (alpha_i x_i)^2 + 2 alpha_i beta_i x_i + beta_i^2 */
      rhs = -consdata->constant;
      for( i = 0; i < consdata->nvars; ++i )
      {
         sqrterm.idx1 = i;
         sqrterm.idx2 = i;
         sqrterm.coef = consdata->coefs[i] * consdata->coefs[i];
         SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, sqrterm) );

         if( ! SCIPisZero(scip, consdata->offsets[i]) )
         {
            rhs -= consdata->offsets[i] * consdata->offsets[i];
            SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, consdata->vars[i], 2.0 * consdata->coefs[i] * consdata->offsets[i]) );
         }
      }

      /* add rhsvar to quadvars of nlrow, if not there yet */
      rhsvarpos = SCIPnlrowSearchQuadVar(consdata->nlrow, consdata->rhsvar);
      if( rhsvarpos == -1 )
      {
         SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, consdata->rhsvar) );
         rhsvarpos = SCIPnlrowSearchQuadVar(consdata->nlrow, consdata->rhsvar);
         assert(rhsvarpos >= 0);
      }

      /* add -(alpha_{n+1} x_{n+1))^2 - 2 alpha_{n+1} beta_{n+1} x_{n+1} - beta_{n+1}^2 */
      sqrterm.idx1 = rhsvarpos;
      sqrterm.idx2 = rhsvarpos;
      sqrterm.coef = -consdata->rhscoeff * consdata->rhscoeff;
      SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, sqrterm) );

      if( consdata->rhsoffset != 0.0 )
      {
         rhs += consdata->rhsoffset * consdata->rhsoffset;
         SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, consdata->rhsvar, -2.0 * consdata->rhscoeff * consdata->rhsoffset) );
      }

      SCIP_CALL( SCIPchgNlRowRhs(scip, consdata->nlrow, rhs) );

      break;
   }

   case 'd':
   {
      /* construct division form (gamma + sum_{i=1}^n (alpha_i(x_i+beta_i))^2)/(alpha_{n+1}(x_{n+1}+beta_{n+1})) <= alpha_{n+1}(x_{n+1}+beta_{n+1}) */
      SCIP_EXPRTREE* exprtree;
      SCIP_EXPR* expr;
      SCIP_EXPR* nominator;
      SCIP_EXPR* denominator;
      SCIP_EXPR** exprs;
      SCIP_EXPRDATA_MONOMIAL** monomials;
      SCIP_Real lincoef;
      SCIP_Real one;
      SCIP_Real two;

      SCIP_CALL( SCIPallocBufferArray(scip, &exprs,     consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomials, consdata->nvars) );
      one = 1.0;
      two = 2.0;

      for( i = 0; i < consdata->nvars; ++i )
      {
         /* put x_i + beta_i into exprs[i] */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprs[i], SCIP_EXPR_VARIDX, i) );
         if( consdata->offsets[i] != 0.0 )
         {
            SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &exprs[i], 1, &exprs[i], &one, consdata->offsets[i]) );
         }

         /* create monomial alpha_i^2 y_i^2, where y_i will be x_i + beta_i */
         SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomials[i], consdata->coefs[i] * consdata->coefs[i], 1, &i, &two) );
      }

      /* setup polynomial expression for gamma + sum_{i=1}^n alpha_i^2 (x_i+beta_i)^2 */
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), &nominator, consdata->nvars, exprs, consdata->nvars, monomials, consdata->constant, FALSE) );  /*lint !e850 */

      SCIPfreeBufferArray(scip, &monomials);
      SCIPfreeBufferArray(scip, &exprs);

      /* setup alpha_{n+1}(x_{n+1}+beta_{n+1})
       * assert that this term is >= 0.0 (otherwise constraint is infeasible anyway) */
      assert(consdata->rhsvar != NULL);
      assert((consdata->rhscoeff >= 0.0 && !SCIPisNegative(scip, SCIPvarGetLbGlobal(consdata->rhsvar) + consdata->rhsoffset)) ||
         (consdata->rhscoeff <= 0.0 && !SCIPisPositive(scip, SCIPvarGetUbGlobal(consdata->rhsvar) + consdata->rhsoffset)));
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &denominator, SCIP_EXPR_VARIDX, consdata->nvars) );
      if( consdata->rhscoeff != 1.0 || consdata->rhsoffset != 0.0 )
      {
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &denominator, 1, &denominator, &consdata->rhscoeff, consdata->rhscoeff * consdata->rhsoffset) );
      }

      /* setup nominator/denominator */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, nominator, denominator) );

      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 0, 0, NULL) );
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, consdata->nvars, consdata->vars) );
      SCIP_CALL( SCIPexprtreeAddVars(exprtree, 1, &consdata->rhsvar) );

      /* linear and constant part is -\alpha_{n+1} (x_{n+1}+\beta_{n+1}) */
      lincoef = -consdata->rhscoeff;
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons),
            -consdata->rhscoeff * consdata->rhsoffset,
            1, &consdata->rhsvar, &lincoef,
            0, NULL, 0, NULL,
            exprtree, -SCIPinfinity(scip), 0.0,
            SCIP_EXPRCURV_UNKNOWN) );

      SCIP_CALL( SCIPexprtreeFree(&exprtree) );

      break;
   }

   default:
      SCIPerrorMessage("unknown value for nlp formulation parameter\n");
      return SCIP_ERROR;
   }

   SCIPdebugMsg(scip, "created nonlinear row representation of SOC constraint\n");
   SCIPdebugPrintCons(scip, cons, NULL);
   SCIPdebug( SCIP_CALL( SCIPprintNlRow(scip, consdata->nlrow, NULL) ) );

   return SCIP_OKAY;
}

/** evaluates the left hand side of a SOC constraint */ 
static
SCIP_RETCODE evalLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to evaluate */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*      var;
   SCIP_Real      val;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->lhsval = consdata->constant;

   for( i = 0; i < consdata->nvars; ++i )
   {
      var = consdata->vars[i];
      val = SCIPgetSolVal(scip, sol, var);

      if( SCIPisInfinity(scip, val) || SCIPisInfinity(scip, -val) )
      {
         consdata->lhsval = SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      val = consdata->coefs[i] * (val + consdata->offsets[i]);
      consdata->lhsval += val * val;      
   }
   consdata->lhsval = sqrt(consdata->lhsval);

   return SCIP_OKAY;
}

/** computes violation of a SOC constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to evaluate */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real rhsval;
   SCIP_Real rhs;
   SCIP_Real absviol;
   SCIP_Real relviol;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( evalLhs(scip, cons, sol) );

   if( SCIPisInfinity(scip, consdata->lhsval) )
   {
      /* infinity <= infinity is feasible
       * infinity <= finite value is not feasible and has violation infinity
       */
      if( (consdata->rhscoeff > 0.0 && SCIPisInfinity(scip,  SCIPgetSolVal(scip, sol, consdata->rhsvar))) ||
         ( consdata->rhscoeff < 0.0 && SCIPisInfinity(scip, -SCIPgetSolVal(scip, sol, consdata->rhsvar))) )
         consdata->violation = 0.0;
      else
         consdata->violation = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   rhsval = SCIPgetSolVal(scip, sol, consdata->rhsvar);
   if( SCIPisInfinity(scip,  rhsval) )
   {
      consdata->violation = consdata->rhscoeff > 0.0 ? 0.0 : SCIPinfinity(scip);
      return SCIP_OKAY;
   }
   if( SCIPisInfinity(scip, -rhsval) )
   {
      consdata->violation = consdata->rhscoeff < 0.0 ? 0.0 : SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   rhs = consdata->rhscoeff * (rhsval + consdata->rhsoffset);
   consdata->violation = consdata->lhsval - rhs;
   absviol = consdata->violation;
   relviol = SCIPrelDiff(consdata->lhsval, rhs);
   if( consdata->violation <= 0.0 )
   {
      /* constraint is not violated for sure */
      consdata->violation = 0.0;
      return SCIP_OKAY;
   }

   if( sol != NULL )
      SCIPupdateSolConsViolation(scip, sol, absviol, relviol);

   return SCIP_OKAY;
}

/** computes violations for a set of constraints */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to evaluate */
   int                   nconss,             /**< number of constraints to evaluate */
   SCIP_SOL*             sol,                /**< solution to evaluate, or NULL if LP solution should be used */
   SCIP_CONS**           maxviolcons         /**< a buffer to store pointer to maximal violated constraint, or NULL if of no interest */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      maxviol = 0.0;
   int            c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   if( maxviolcons != NULL )
      *maxviolcons = NULL;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conshdlr, conss[c], sol) );  /*lint !e613*/
      if( maxviolcons != NULL )
      {
         consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
         assert(consdata != NULL);
         if( consdata->violation > maxviol && SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) )
         {
            maxviol      = consdata->violation;
            *maxviolcons = conss[c];  /*lint !e613*/
         }
      }
   }

   return SCIP_OKAY;
}

/** generate supporting hyperplane in a given solution */
static
SCIP_RETCODE generateCutSol(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_ROWPREP**        rowprep             /**< place to store cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      val;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(rowprep != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(SCIPisPositive(scip, consdata->lhsval)); /* do not like to linearize in 0 */
   assert(!SCIPisInfinity(scip, consdata->lhsval));

   SCIP_CALL( SCIPcreateRowprep(scip, rowprep, SCIP_SIDETYPE_RIGHT, SCIPconsIsLocal(cons)) );
   SCIP_CALL( SCIPensureRowprepSize(scip, *rowprep, consdata->nvars+1) );
   (void) SCIPsnprintf((*rowprep)->name, SCIP_MAXSTRLEN, "%s_linearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   for( i = 0; i < consdata->nvars; ++i )
   {
      val  = SCIPgetSolVal(scip, sol, consdata->vars[i]) + consdata->offsets[i];
      val *= consdata->coefs[i] * consdata->coefs[i];

      SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, consdata->vars[i], val / consdata->lhsval) );

      val *= SCIPgetSolVal(scip, sol, consdata->vars[i]);
      SCIPaddRowprepSide(*rowprep, val);
   }
   (*rowprep)->side /= consdata->lhsval;
   (*rowprep)->side -= consdata->lhsval - consdata->rhscoeff * consdata->rhsoffset;

   SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, consdata->rhsvar, -consdata->rhscoeff) );

   return SCIP_OKAY;
}

/** generate supporting hyperplane in a given point */
static
SCIP_RETCODE generateCutPoint(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            x,                  /**< point (lhs-vars) where to generate cut */
   SCIP_ROWPREP**        rowprep             /**< place to store cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      lhsval;
   SCIP_Real      val;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(rowprep != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhsval = consdata->constant;
   for( i = 0; i < consdata->nvars; ++i )
   {
      assert(!SCIPisInfinity(scip, ABS(x[i])));

      val = consdata->coefs[i] * (x[i] + consdata->offsets[i]);
      lhsval += val * val;      
   }
   lhsval = sqrt(lhsval);

   if( SCIPisZero(scip, lhsval) )
   { /* do not like to linearize in 0 */
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateRowprep(scip, rowprep, SCIP_SIDETYPE_RIGHT, SCIPconsIsLocal(cons)) );
   SCIP_CALL( SCIPensureRowprepSize(scip, *rowprep, consdata->nvars+1) );
   (void) SCIPsnprintf((*rowprep)->name, SCIP_MAXSTRLEN, "%s_linearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   for( i = 0; i < consdata->nvars; ++i )
   {
      val  = x[i] + consdata->offsets[i];
      if( SCIPisZero(scip, val) )
         continue;

      val *= consdata->coefs[i] * consdata->coefs[i];

      SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, consdata->vars[i], val / lhsval) );

      val *= x[i];
      SCIPaddRowprepSide(*rowprep, val);
   }
   (*rowprep)->side /= lhsval;
   (*rowprep)->side -= lhsval - consdata->rhscoeff * consdata->rhsoffset;

   SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, consdata->rhsvar, -consdata->rhscoeff) );

   return SCIP_OKAY;
}

/** generate supporting hyperplane w.r.t. solution projected on feasible set 
 * 
 * Instead of linearizing the SOC constraint in the given solution point, this function projects the point
 * first onto the feasible set of the SOC constraint (w.r.t. euclidean norm (scaled by alpha))
 * and linearizes the SOC constraint there.
 * The hope is that this produces somewhat tighter cuts.
 * 
 * The projection has only be computed for the case gamma = 0.
 * If gamma > 0, generateCut is called. 
 * 
 * Let \f$\hat x\f$ be sol or the LP solution if sol == NULL.
 * Then the projected point \f$\tilde x\f$ is given by
 * \f[
 *   \tilde x_i = \frac{\hat x_i + \lambda \beta_i}{1-\lambda},  \quad i=1,\ldots, n; \quad
 *   \tilde x_{n+1} = \frac{\hat x_{n+1} - \lambda \beta_{n+1}}{1+\lambda}
 * \f]
 * where
 * \f[
 *   \lambda = \frac{1-A}{1+A}, \qquad 
 *   A = \frac{\alpha_{n+1}(\hat x_{n+1}+\beta_{n+1})}{\sqrt{\sum_{i=1}^n (\alpha_i(\hat x_i+\beta_i))^2}}
 * \f]
 * 
 * If lambda is very close to 1, generateCut is called.
 * 
 * The generated cut is very similar to the unprojected form.
 * The only difference is in the right hand side, which is (in the case beta = 0) multiplied by 1/(1-lambda).
 */
static
SCIP_RETCODE generateCutProjectedPoint(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_ROWPREP**        rowprep             /**< place to store cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      val;
   SCIP_Real      A, lambda;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(rowprep != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(SCIPisPositive(scip, consdata->lhsval)); /* do not like to linearize in 0 */
   assert(!SCIPisInfinity(scip, consdata->lhsval));

   if( !SCIPisZero(scip, consdata->constant) )
   {  /* have not thought about this case yet */
      SCIP_CALL( generateCutSol(scip, cons, sol, rowprep) );
      return SCIP_OKAY;
   }

   A  = consdata->rhscoeff * (SCIPgetSolVal(scip, sol, consdata->rhsvar) + consdata->rhsoffset);
   A /= consdata->lhsval;

   lambda = (1.0 - A) / (1.0 + A);

   assert(!SCIPisNegative(scip, lambda)); /* otherwise A > 1, so constraint is not violated */

   SCIPdebugMsg(scip, "A = %g \t lambda = %g\n", A, lambda);

   if( SCIPisFeasEQ(scip, lambda, 1.0) )
   {  /* avoid numerical difficulties when dividing by (1-lambda) below */ 
      SCIP_CALL( generateCutSol(scip, cons, sol, rowprep) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateRowprep(scip, rowprep, SCIP_SIDETYPE_RIGHT, SCIPconsIsLocal(cons)) );
   SCIP_CALL( SCIPensureRowprepSize(scip, *rowprep, consdata->nvars+1) );
   (void) SCIPsnprintf((*rowprep)->name, SCIP_MAXSTRLEN, "%s_linearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   for( i = 0; i < consdata->nvars; ++i )
   {
      val  = SCIPgetSolVal(scip, sol, consdata->vars[i]) + consdata->offsets[i];
      val *= consdata->coefs[i] * consdata->coefs[i];

      SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, consdata->vars[i], val / consdata->lhsval) );

      val *= SCIPgetSolVal(scip, sol, consdata->vars[i]) + lambda * consdata->offsets[i];
      SCIPaddRowprepSide(*rowprep, val);
   }
   (*rowprep)->side /= consdata->lhsval;
   (*rowprep)->side-= consdata->lhsval;
   (*rowprep)->side /= 1.0 - lambda;
   (*rowprep)->side -= consdata->rhscoeff * consdata->rhsoffset;

   SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, consdata->rhsvar, -consdata->rhscoeff) );

   return SCIP_OKAY;
}

/** generates sparsified supporting hyperplane */
static
SCIP_RETCODE generateSparseCut(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_ROWPREP**        rowprep,            /**< place to store cut */
   SCIP_Real             minefficacy         /**< minimal efficacy for a cut to be accepted */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real*     x;
   SCIP_Real*     dist;  /* distance to 0 */
   int*           ind;   /* indices */
   int            i;
   int            maxnz, nextmaxnz;
   SCIP_Real      efficacy;
   SCIP_Real      goodefficacy;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(rowprep != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(SCIPisPositive(scip, consdata->lhsval)); /* do not like to linearize in 0 */
   assert(!SCIPisInfinity(scip, consdata->lhsval));

   if( consdata->nvars <= 3 )
   {
      SCIP_CALL( generateCutSol(scip, cons, sol, rowprep) );
      return SCIP_OKAY;
   }

   goodefficacy = MAX((1.0-conshdlrdata->sparsifymaxloss) * consdata->violation, minefficacy);

   SCIP_CALL( SCIPallocBufferArray(scip, &x,    consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dist, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind,  consdata->nvars) );

   SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, x) );
   /* distance to "-offset" * alpha_i^2 should indicate loss when moving refpoint to x[i] = -offset[i] */
   for( i = 0; i < consdata->nvars; ++i )
   {
      ind[i] = i;
      dist[i]  = ABS(x[i] + consdata->offsets[i]);
      dist[i] *= consdata->coefs[i] * consdata->coefs[i];
   }

   /* sort variables according to dist */
   SCIPsortRealInt(dist, ind, consdata->nvars);

   maxnz = 2;
   /* set all except biggest maxnz entries in x to -offset */
   for( i = 0; i < consdata->nvars - maxnz; ++i )
      x[ind[i]] = -consdata->offsets[i];

   do
   {
      /* @todo speed up a bit by computing efficacy of new cut from efficacy of old cut
       * generate row only if efficient enough
       */
      SCIP_CALL( generateCutPoint(scip, cons, x, rowprep) );

      if( *rowprep != NULL )
      {
         efficacy = SCIPgetRowprepViolation(scip, *rowprep, sol);
         if( SCIPisGT(scip, efficacy, goodefficacy) ||
            (maxnz >= consdata->nvars && SCIPisGT(scip, efficacy, minefficacy)) )
         {
            /* cut cuts off solution and is efficient enough */
            SCIPdebugMsg(scip, "accepted cut with %d of %d nonzeros, efficacy = %g\n", maxnz, consdata->nvars, efficacy);
            break;
         }

         SCIPfreeRowprep(scip, rowprep);
      }

      /* cut also not efficient enough if generated in original refpoint (that's bad) */
      if( maxnz >= consdata->nvars )
         break;

      nextmaxnz = (int)(conshdlrdata->sparsifynzgrowth * maxnz);
      if( nextmaxnz == consdata->nvars - 1)
         nextmaxnz = consdata->nvars;
      else if( nextmaxnz == maxnz )
         ++nextmaxnz;

      /* restore entries of x that are nonzero in next attempt */
      for( i = MAX(0, consdata->nvars - nextmaxnz); i < consdata->nvars - maxnz; ++i )
         x[ind[i]] = SCIPgetSolVal(scip, sol, consdata->vars[ind[i]]);

      maxnz = nextmaxnz;
   }
   while( TRUE );  /*lint !e506*/

   SCIPfreeBufferArray(scip, &x);
   SCIPfreeBufferArray(scip, &dist);
   SCIPfreeBufferArray(scip, &ind);

   return SCIP_OKAY;
}

/** separates a point, if possible */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_Bool             inenforcement,      /**< whether we are in constraint enforcement */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a fixing leads to a cutoff */
   SCIP_Bool*            success             /**< buffer to store whether the point was separated */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          minefficacy;
   int                c;
   SCIP_ROW*          row;
   SCIP_ROWPREP*      rowprep;

   assert(scip    != NULL);
   assert(conss   != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(cutoff != NULL);
   assert(success != NULL);

   *cutoff = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *success = FALSE;

   minefficacy = inenforcement ? SCIPlpfeastol(scip) : SCIPgetSepaMinEfficacy(scip);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) && !SCIPisInfinity(scip, consdata->violation) )
      {
         SCIP_Real efficacy;

         rowprep = NULL;

         /* generate cut */
         if( conshdlrdata->sparsify )
         {
            SCIP_CALL( generateSparseCut(scip, conshdlr, conss[c], sol, &rowprep, minefficacy) );  /*lint !e613*/
         }
         else if( conshdlrdata->projectpoint )
         {
            SCIP_CALL( generateCutProjectedPoint(scip, conss[c], sol, &rowprep) );  /*lint !e613*/
         }
         else
         {
            SCIP_CALL( generateCutSol(scip, conss[c], sol, &rowprep) );  /*lint !e613*/
         }

         if( rowprep == NULL )
            continue;

         /* NOTE: The way that rowprep was constructed, there should be no need to call SCIPmergeRowprep,
          * since no variable gets added twice. However, if rowprep were replacing multiaggregated variables
          * (as there can exist for soc cons), then SCIPmergeRowprep would be necessary.
          */
         /* cleanup rowprep (there is no limit on coefrange for cons_soc) TODO add a coefrange limit? */
         SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIPinfinity(scip), minefficacy, NULL, &efficacy) );

         if( SCIPisLE(scip, efficacy, minefficacy) )
         {
            SCIPfreeRowprep(scip, &rowprep);
            continue;
         }

         /* cut cuts off solution and efficient enough */
         SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, conshdlr) );
         if( SCIPisCutApplicable(scip, row) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );  /*lint !e613*/

            *success = TRUE;

            SCIPdebugMsg(scip, "added cut with efficacy %g\n", SCIPgetCutEfficacy(scip, sol, row));

            /* mark row as not removable from LP for current node, if in enforcement */
            if( inenforcement && !conshdlrdata->enfocutsremovable )
               SCIPmarkRowNotRemovableLocal(scip, row);
         }

         SCIP_CALL( SCIPreleaseRow (scip, &row) );
         SCIPfreeRowprep(scip, &rowprep);
      }

      if( *cutoff )
         break;

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */ 
      if( c >= nusefulconss && *success )
         break;
   }

   return SCIP_OKAY;
}

/** adds linearizations cuts for convex constraints w.r.t. a given reference point to cutpool and sepastore
 *
 * If separatedlpsol is not NULL, then a cut that separates the LP solution is added to the sepastore and is forced to enter the LP.
 * If separatedlpsol is not NULL, but cut does not separate the LP solution, then it is added to the cutpool only.
 * If separatedlpsol is NULL, then cut is added to cutpool only.
 */
static
SCIP_RETCODE addLinearizationCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             ref,                /**< reference point where to linearize, or NULL for LP solution */
   SCIP_Bool*            separatedlpsol,     /**< buffer to store whether a cut that separates the current LP solution was found and added to LP, or NULL if adding to cutpool only */
   SCIP_Real             minefficacy,        /**< minimal efficacy of a cut when checking for separation of LP solution */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff was detected */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool addedtolp;
   SCIP_ROWPREP* rowprep;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(cutoff != NULL);
   *cutoff = FALSE;

   if( separatedlpsol != NULL )
      *separatedlpsol = FALSE;

   for( c = 0; c < nconss && !(*cutoff); ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613 */

      if( SCIPconsIsLocal(conss[c]) )  /*lint !e613 */
         continue;

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
      assert(consdata != NULL);

      SCIP_CALL( evalLhs(scip, conss[c], ref) );  /*lint !e613 */
      if( !SCIPisPositive(scip, consdata->lhsval) || SCIPisInfinity(scip, consdata->lhsval) )
      {
         SCIPdebugMsg(scip, "skip adding linearization for <%s> since lhs is %g\n", SCIPconsGetName(conss[c]), consdata->lhsval);  /*lint !e613 */
         continue;
      }

      SCIP_CALL( generateCutSol(scip, conss[c], ref, &rowprep) );  /*lint !e613 */

      if( rowprep == NULL )
         continue;

      addedtolp = FALSE;

      /* if caller wants, then check if cut separates LP solution and add to sepastore if so */
      if( separatedlpsol != NULL )
      {
         if( SCIPgetRowprepViolation(scip, rowprep, NULL) >= minefficacy )
         {
            SCIP_ROW* row;

            SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, conshdlr) );
            SCIP_CALL( SCIPaddRow(scip, row, TRUE, cutoff) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );

            *separatedlpsol = TRUE;
            addedtolp = TRUE;
         }
      }

      if( !addedtolp && !rowprep->local )
      {
         SCIP_ROW* row;

         SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, conshdlr) );
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }

      SCIPfreeRowprep(scip, &rowprep);
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
   SCIP_Bool      cutoff;

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

   SCIPdebugMsg(scip, "caught new sol event %" SCIP_EVENTTYPE_FORMAT " from heur <%s>; have %d conss\n", SCIPeventGetType(event),
      SCIPheurGetName(SCIPsolGetHeur(sol)), nconss);

   SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, sol, NULL, 0.0, &cutoff) );
   /* ignore cutoff, cannot return status */

   return SCIP_OKAY;
}

/** removes fixed variables, replace aggregated and negated variables
 *
 * repeats replacements until no further change is found;
 * takes care of capture/release and locks, but not of variable events (assumes that var events are not caught yet) 
 */
static
SCIP_RETCODE presolveRemoveFixedVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for signpower constraints */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  ndelconss,          /**< counter for number of deleted constraints */
   int*                  nupgdconss,         /**< counter for number of upgraded constraints */
   int*                  nchgbds,            /**< counter for number of bound changes */
   int*                  nfixedvars,         /**< counter for number of fixed variables */
   SCIP_Bool*            iscutoff,           /**< to store whether constraint cannot be satisfied */
   SCIP_Bool*            isdeleted           /**< to store whether constraint is redundant and can be deleted */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool      havechange;
   SCIP_Bool      haveremovedvar;
   int            i;
   SCIP_VAR*      x;
   SCIP_Real      coef;
   SCIP_Real      offset;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(iscutoff != NULL);
   assert(isdeleted != NULL);

   *iscutoff  = FALSE;
   *isdeleted = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMsg(scip, "remove fixed variables from constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   havechange     = FALSE;
   haveremovedvar = FALSE;

   /* process variables on left hand side */
   for( i = 0; i < consdata->nvars; ++i )
   {
      x = consdata->vars[i];
      assert(x != NULL);
      assert(SCIPvarGetStatus(x) != SCIP_VARSTATUS_ORIGINAL);

      if( SCIPvarIsActive(x) || SCIPvarGetStatus(x) == SCIP_VARSTATUS_MULTAGGR )
         continue;

      havechange = TRUE;

      /* drop variable event and unlock and release variable */
      SCIP_CALL( dropLhsVarEvents(scip, conshdlrdata->eventhdlr, cons, i) );
      SCIP_CALL( SCIPunlockVarCons(scip, x, cons, TRUE, TRUE) );
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->vars[i]) );

      coef = 1.0;
      offset = consdata->offsets[i];
      SCIP_CALL( SCIPgetProbvarSum(scip, &x, &coef, &offset) );

      SCIPdebugMsg(scip, "  lhs term at position %d is replaced by %g * <%s> + %g\n",
         i, coef, SCIPvarGetName(x), offset);

      /* if variable has been fixed, add (alpha*offset)^2 to gamma and continue */
      if( coef == 0.0 || x == NULL )
      {
         consdata->constant  += consdata->coefs[i] * consdata->coefs[i] * offset * offset;
         consdata->offsets[i] = 0.0;
         haveremovedvar = TRUE;
         continue;
      }

      assert(SCIPvarIsActive(x) || SCIPvarGetStatus(x) == SCIP_VARSTATUS_MULTAGGR);

      /* replace coefs[i] * (vars[i] + offsets[i]) by coefs[i]*coef * (x + offsets[i]/coef) */
      consdata->offsets[i] = offset;
      if( coef != 1.0 )
      {
         consdata->coefs[i]    = REALABS(coef * consdata->coefs[i]);
         consdata->offsets[i] /= coef;
      }
      consdata->vars[i] = x;

      /* capture and lock new variable, catch variable events */
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
      SCIP_CALL( SCIPlockVarCons(scip, consdata->vars[i], cons, TRUE, TRUE) );
      SCIP_CALL( catchLhsVarEvents(scip, conshdlrdata->eventhdlr, cons, i) );
   }

   /* process variable on right hand side */
   x = consdata->rhsvar;
   assert(x != NULL);
   if( !SCIPvarIsActive(x) && SCIPvarGetStatus(x) != SCIP_VARSTATUS_MULTAGGR )
   {
      havechange = TRUE;

      /* drop variable event and unlock and release variable */
      SCIP_CALL( dropRhsVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->rhsvar) );
      SCIP_CALL( SCIPunlockVarCons(scip, x, cons, consdata->rhscoeff > 0.0, consdata->rhscoeff < 0.0) );

      coef = 1.0;
      offset = 0.0;
      SCIP_CALL( SCIPgetProbvarSum(scip, &x, &coef, &offset) );

      SCIPdebugMsg(scip, "  rhs variable is replaced by %g * <%s> + %g\n", coef, SCIPvarGetName(x), offset);

      if( coef == 0.0 || x == NULL )
      {
         /* if variable has been fixed, add offset to rhsoffset */
         consdata->rhsoffset += offset;
      }
      else
      {
         /* replace rhscoef * (rhsvar + rhsoffset) by rhscoef*coef * (x + offset/coef + rhsoffset/coef) */
         assert(SCIPvarIsActive(x) || SCIPvarGetStatus(x) == SCIP_VARSTATUS_MULTAGGR);

         consdata->rhsoffset = (consdata->rhsoffset + offset) / coef;
         consdata->rhscoeff *= coef;
         consdata->rhsvar = x;

         /* capture and lock new variable, catch variable events */
         SCIP_CALL( SCIPcaptureVar(scip, consdata->rhsvar) );
         SCIP_CALL( SCIPlockVarCons(scip, consdata->rhsvar, cons, consdata->rhscoeff > 0.0, consdata->rhscoeff < 0.0) );
         SCIP_CALL( catchRhsVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      }
   }

   if( !havechange )
      return SCIP_OKAY;

   /* free nonlinear row representation */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* if a variable has been removed, close gaps in vars array */
   if( haveremovedvar )
   {
      int oldnvars;

      /* due to the realloc of the block memory below and the way we store the eventdata in consdata, we best drop all events here and catch them again below */
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );

      oldnvars = consdata->nvars;
      for( i = 0; i < consdata->nvars; ++i )
      {
         /* forget about empty places at end of vars array */
         while( consdata->nvars && consdata->vars[consdata->nvars-1] == NULL )
            --consdata->nvars;

         /* all variables at index >= i have been removed */
         if( i == consdata->nvars )
            break;

         if( consdata->vars[i] != NULL )
            continue;

         /* move variable from position nvars-1 to position i */

         assert(consdata->nvars >= 1);
         assert(consdata->vars[consdata->nvars-1] != NULL);

         consdata->vars[i]    = consdata->vars[consdata->nvars-1];
         consdata->offsets[i] = consdata->offsets[consdata->nvars-1];
         consdata->coefs[i]   = consdata->coefs[consdata->nvars-1];

         --consdata->nvars;
      }

      assert(consdata->nvars < oldnvars);

      /* shrink arrays in consdata */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars,    oldnvars, consdata->nvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->offsets, oldnvars, consdata->nvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->coefs,   oldnvars, consdata->nvars) );

      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   SCIPdebugMsg(scip, "\t-> ");
   SCIPdebugPrintCons(scip, cons, NULL);

   if( consdata->nvars == 0 )
   { /* all variables on left hand size have been removed, remaining constraint is sqrt(gamma) <= ... */
      assert(!SCIPisNegative(scip, consdata->constant));
      if( consdata->rhsvar == NULL )
      { /* also rhsvar has been removed, remaining constraint is sqrt(gamma) <= rhscoeff * rhsoffset */
         if( SCIPisFeasLE(scip, sqrt(consdata->constant), consdata->rhscoeff*consdata->rhsoffset) )
         {
            SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing all variables\n", SCIPconsGetName(cons));
         }
         else
         {
            SCIPdebugMsg(scip, "found problem infeasible after fixing all variables in <%s>\n", SCIPconsGetName(cons));
            *iscutoff = TRUE;
         }
         ++*ndelconss;
      }
      else if( !SCIPvarIsActive(consdata->rhsvar) )
      { /* remaining constraint is sqrt(gamma) - rhscoeff * rhsoffset <= rhscoeff * rhsvar, and rhsvar is probably multi-aggregated */
         SCIP_CONS* lincons;

         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 1, &consdata->rhsvar, &consdata->rhscoeff,
               sqrt(consdata->constant) - consdata->rhscoeff * consdata->rhsoffset, SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
         ++*nupgdconss;
      }
      else if( consdata->rhscoeff > 0.0 )
      { /* remaining constraint is sqrt(gamma) / rhscoeff - rhsoffset <= rhsvar */
         SCIP_Bool tightened;
         SCIP_CALL( SCIPtightenVarLb(scip, consdata->rhsvar, sqrt(consdata->constant) / consdata->rhscoeff - consdata->rhsoffset, TRUE, iscutoff, &tightened) );
         if( *iscutoff )
         {
            SCIPdebugMsg(scip, "found problem infeasible after fixing all lhs variables in <%s> and tightening lower bound of rhs var\n", SCIPconsGetName(cons));
         }
         else if( tightened )
         {
            SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing all lhs variables and tightening lower bound of rhs var\n", SCIPconsGetName(cons));
            ++*nchgbds;
         }
         else
         {
            SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing all lhs variables\n", SCIPconsGetName(cons));
         }
         ++*ndelconss;
      }
      else
      { /* remaining constraint is sqrt(gamma) / rhscoeff - rhsoffset >= rhsvar */
         SCIP_Bool tightened;
         SCIP_CALL( SCIPtightenVarUb(scip, consdata->rhsvar, sqrt(consdata->constant) / consdata->rhscoeff - consdata->rhsoffset, TRUE, iscutoff, &tightened) );
         if( *iscutoff )
         {
            SCIPdebugMsg(scip, "found problem infeasible after fixing all lhs variables in <%s> and tightening upper bound of rhs var\n", SCIPconsGetName(cons));
         }
         else if( tightened )
         {
            SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing all lhs variables and tightening upper bound of rhs var\n", SCIPconsGetName(cons));
            ++*nchgbds;
         }
         else
         {
            SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing all lhs variables\n", SCIPconsGetName(cons));
         }
         ++*ndelconss;
      }
      SCIP_CALL( SCIPdelCons(scip, cons) );
      *isdeleted = TRUE;
      return SCIP_OKAY;
   }

   if( consdata->rhsvar == NULL )
   { /* constraint becomes sum_i (alpha_i*(x_i+beta_i))^2 <= (rhscoeff*rhsoffset)^2 - gamma */
      if( consdata->nvars > 1 )
      { /* upgrade to quadratic constraint */
         SCIP_CONS* quadcons;
         SCIP_QUADVARTERM* quadvarterms;
         SCIP_Real  rhs;

         SCIP_CALL( SCIPallocBufferArray(scip, &quadvarterms, consdata->nvars) );
         BMSclearMemoryArray(quadvarterms, consdata->nvars);
         rhs = consdata->rhscoeff * consdata->rhsoffset;
         rhs = rhs * rhs - consdata->constant;

         for( i = 0; i < consdata->nvars; ++i )
         {
            quadvarterms[i].var = consdata->vars[i];
            quadvarterms[i].sqrcoef = consdata->coefs[i] * consdata->coefs[i];
            if( consdata->offsets[i] != 0.0 )
            {
               quadvarterms[i].lincoef = 2 * consdata->offsets[i] * quadvarterms[i].sqrcoef;
               rhs -= quadvarterms[i].sqrcoef * consdata->offsets[i]*consdata->offsets[i];
            }
         }

         assert(!SCIPconsIsStickingAtNode(cons));
         SCIP_CALL( SCIPcreateConsQuadratic2(scip, &quadcons, SCIPconsGetName(cons), 0, NULL, NULL,
               consdata->nvars, quadvarterms, 0, NULL, -SCIPinfinity(scip), rhs,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddCons(scip, quadcons) );
         SCIPdebugMsg(scip, "upgraded <%s> to quadratic constraint: ", SCIPconsGetName(cons));
         SCIPdebugPrintCons(scip, quadcons, NULL);

         SCIP_CALL( SCIPreleaseCons(scip, &quadcons) );

         SCIPfreeBufferArray(scip, &quadvarterms);

         ++*nupgdconss;
      }
      else if( !SCIPvarIsActive(consdata->vars[0]) )
      { /* constraint is |alpha*(x+beta)| <= sqrt((rhscoeff*rhsoffset)^2 - gamma), but x is probably multaggr. -> turn into ranged linear constraint */
         SCIP_CONS* lincons;

         /* create constraint alpha*x <=  sqrt((rhscoeff*rhsoffset)^2 - gamma) - alpha*beta
          *                   alpha*x >= -sqrt((rhscoeff*rhsoffset)^2 - gamma) - alpha*beta */
         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 1, &consdata->vars[0], &consdata->coefs[0],
               -sqrt(consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset * consdata->rhsoffset - consdata->constant) - consdata->coefs[0] * consdata->offsets[0],
               +sqrt(consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset * consdata->rhsoffset - consdata->constant) - consdata->coefs[0] * consdata->offsets[0],
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

         ++*nupgdconss;
      }
      else
      { /* constraint is |alpha*(x+beta)| <= sqrt((rhscoeff*rhsoffset)^2 - gamma) -> propagate bounds */
         SCIP_Bool tightened;
         SCIP_Real rhs;
         assert(consdata->nvars == 1); /* case == 0 handled before */
         rhs = consdata->rhscoeff * consdata->rhsoffset;
         rhs = rhs * rhs;
         if( SCIPisNegative(scip, rhs - consdata->constant) )
         { /* take this as infeasible */
            SCIPdebugMsg(scip, "found problem infeasible after fixing rhs and all except one lhs variables in <%s>\n", SCIPconsGetName(cons));
            *iscutoff = TRUE;
         }
         else
         {
            rhs -= consdata->constant;
            rhs  = rhs < 0.0 ? 0.0 : sqrt(rhs);

            if( SCIPisZero(scip, rhs) )
            { /* constraint is x = -beta */
               SCIP_CALL( SCIPfixVar(scip, consdata->vars[0], -consdata->offsets[0], iscutoff, &tightened) );
               if( *iscutoff )
               {
                  SCIPdebugMsg(scip, "found problem infeasible after fixing rhs and all except one lhs variables and fixing remaining lhs var in <%s>\n", SCIPconsGetName(cons));
               }
               else if( tightened )
               {
                  SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing rhs and all except one lhs variables and fixing remaining lhs var\n", SCIPconsGetName(cons));
                  ++*nfixedvars;
               }
               else
               {
                  SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing rhs and all except one lhs variables and fixing remaining lhs var\n", SCIPconsGetName(cons));
               }
            }
            else
            { /* constraint is -rhs/|alpha| - beta <= x <= rhs/|alpha| - beta */
               rhs /= ABS(consdata->coefs[0]);
               SCIP_CALL( SCIPtightenVarLb(scip, consdata->vars[0], -rhs - consdata->offsets[0], TRUE, iscutoff, &tightened) );
               if( *iscutoff )
               {
                  SCIPdebugMsg(scip, "found problem infeasible after fixing rhs and all except one lhs variables and tightening lower bound of remaining lhs var in <%s>\n", SCIPconsGetName(cons));
               }
               else
               {
                  if( tightened )
                     ++*nchgbds;
                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vars[0], rhs - consdata->offsets[0], TRUE, iscutoff, &tightened) );
                  if( *iscutoff )
                  {
                     SCIPdebugMsg(scip, "found problem infeasible after fixing rhs and all except one lhs variables and tightening upper bound of remaining lhs var in <%s>\n", SCIPconsGetName(cons));
                  }
                  else if( tightened )
                     ++*nchgbds;
               }
               if( !*iscutoff )
               {
                  SCIPdebugMsg(scip, "remove redundant constraint <%s> after fixing rhs and all except one lhs variables and tightening bounds on remaining lhs var\n", SCIPconsGetName(cons));
               }
            }
         }
         ++*ndelconss;
      }
      *isdeleted = TRUE;
      SCIP_CALL( SCIPdelCons(scip, cons) );
      return SCIP_OKAY;
   }

   if( consdata->nvars == 1 && SCIPisZero(scip, consdata->constant) )
   { /* one variable on lhs left and no constant, constraint becomes |alpha*(x+beta)| <= rhscoef*(rhsvar+rhsoffset) -> upgrade to two linear constraints */
      SCIP_CONS* lincons;
      SCIP_VAR*  vars[2];
      SCIP_Real  coefs[2];
      SCIP_Real  rhs;
      assert(consdata->rhsvar != NULL); /* case == NULL has been handled before */

      vars[0] = consdata->vars[0];
      vars[1] = consdata->rhsvar;
      coefs[0] = consdata->coefs[0];
      coefs[1] = -consdata->rhscoeff;
      rhs = consdata->rhscoeff * consdata->rhsoffset - coefs[0] * consdata->offsets[0];

      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 2, vars, coefs, -SCIPinfinity(scip), rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

      coefs[0] = -coefs[0];
      rhs = consdata->rhscoeff * consdata->rhsoffset - coefs[0] * consdata->offsets[0];

      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 2, vars, coefs, -SCIPinfinity(scip), rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

      SCIPdebugMsg(scip, "upgraded <%s> to two linear constraint\n", SCIPconsGetName(cons));

      ++*nupgdconss;
      SCIP_CALL( SCIPdelCons(scip, cons) );
      *isdeleted = TRUE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** adds the linear outer-approximation of Glineur et.al. for a SOC constraint of dimension 3
 *
 * Input is the data for a constraint \f$\sqrt{(\alpha_1(x_1+offset1))^2 + (\alpha_2(x_2+offset2))^2)} \leq \alpha_3(x_3+offset3)\f$.
 * Here \f$\alpha3 > 0\f$, and the lower bound of \f$x_3 \geq -offset3\f$.
 * Also x2 = NULL is allowed, in which case the second term is assumed to be constant, and \f$offset2 \neq 0\f$ is needed.
 */
static
SCIP_RETCODE presolveCreateGlineurApproxDim3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< original constraint */
   SCIP_VAR*             x1,                 /**< variable x1 */
   SCIP_VAR*             x2,                 /**< variable x2 */
   SCIP_VAR*             x3,                 /**< variable x3 */
   SCIP_Real             alpha1,             /**< coefficient of x1 */
   SCIP_Real             alpha2,             /**< coefficient of x2 */
   SCIP_Real             alpha3,             /**< coefficient of x3 */
   SCIP_Real             offset1,            /**< offset of x1 */
   SCIP_Real             offset2,            /**< offset of x2 */
   SCIP_Real             offset3,            /**< offset of x3 */
   int                   N,                  /**< size of linear approximation, need to be >= 1 */
   const char*           basename,           /**< string to use for building variable and constraint names */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   SCIP_CONS*     lincons;
   SCIP_VAR*      vars[3];
   SCIP_Real      vals[3];
   char           varname[255];
   char           linname[255];
   int            i;
   SCIP_VAR**     avars;
   SCIP_VAR**     bvars;
   SCIP_Real      val;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x1   != NULL);
   assert(x2   != NULL || !SCIPisZero(scip, offset2));
   assert(x3   != NULL);
   assert(SCIPisPositive(scip, alpha3));
   assert(SCIPisGE(scip, SCIPconsIsLocal(cons) ? SCIPvarGetLbLocal(x3) : SCIPvarGetLbGlobal(x3), -offset3));
   assert(basename != NULL);
   assert(N >= 1);
   assert(naddconss != NULL);

   SCIPdebugMsg(scip, "Creating linear Glineur outer-approximation for <%s>.\n", basename);
   SCIPdebugMsg(scip, "sqr(%g(%s+%g)) + sqr(%g(%s+%g)) <= sqr(%g(%s+%g)).\n",
      alpha1, SCIPvarGetName(x1), offset1, alpha2, x2 ? SCIPvarGetName(x2) : "0", offset2, alpha3, SCIPvarGetName(x3), offset3);

   SCIP_CALL( SCIPallocBufferArray(scip, &avars, N+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bvars, N+1) );

   /* create additional variables; we do not use avars[0] and bvars[0] */
   for( i = 1; i <= N; ++i )
   {
      (void) SCIPsnprintf(varname, 255, "soc#%s_a%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &avars[i], varname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, 
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, avars[i]) );

      (void) SCIPsnprintf(varname, 255, "soc#%s_b%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &bvars[i], varname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, bvars[i]) );
   }

   /* create linear constraints for the first case
    * cos(pi) = -1, sin(pi) = 0
    * -> a_1  = - alpha1 (x1 + offset1)    ->  -alpha1*x1 - a_1  =  alpha1*offset1 
    * -> b_1 >= | alpha2 (x2 + offset2) |  ->   alpha2*x2 - b_1 <= -alpha2*offset2
    *                                           alpha2*x2 + b_1 >= -alpha2*offset2
    */

   vars[0] = x1;
   vals[0] = -alpha1;
   vars[1] = avars[1];
   vals[1] = -1.0;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, 0);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, alpha1*offset1, alpha1*offset1,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIPdebugPrintCons(scip, lincons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++(*naddconss);

   if( x2 != NULL )
   {
      vars[0] = x2;
      vals[0] = alpha2;
      vars[1] = bvars[1];
      vals[1] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -SCIPinfinity(scip), -alpha2*offset2,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);

      vars[0] = x2;
      vals[0] = alpha2;
      vars[1] = bvars[1];
      vals[1] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha2*offset2, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);
   }
   else
   { /* x2 == NULL ->  b_1 >= |alpha2*offset2| */
      SCIP_Bool infeas;
      SCIP_Bool tightened;
      SCIP_CALL( SCIPtightenVarLb(scip, bvars[1], ABS(alpha2 * offset2), TRUE, &infeas, &tightened) );
      if( infeas == TRUE )
      {
         SCIPwarningMessage(scip, "creating glineur outer approximation of SOC3 constraint found problem infeasible.\n");
      }
   }

   /* create intermediate linear constraints */
   val = M_PI;
   for( i = 1; i < N; ++i )
   {
      val /= 2.0;

      vars[0] = avars[i];
      vals[0] = cos(val);
      vars[1] = bvars[i];
      vals[1] = sin(val);
      vars[2] = avars[i+1]; /*lint !e679*/
      vals[2] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);

      vars[0] = avars[i];
      vals[0] = -sin(val);
      vars[1] = bvars[i];
      vals[1] = cos(val);
      vars[2] = bvars[i+1]; /*lint !e679*/
      vals[2] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, -SCIPinfinity(scip), 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);

      vars[0] = avars[i];
      vals[0] = -sin(val);
      vars[1] = bvars[i];
      vals[1] = cos(val);
      vars[2] = bvars[i+1]; /*lint !e679*/
      vals[2] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebugPrintCons(scip, lincons, NULL);
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);
   }

   /* create last linear constraint */
   val /= 2.0;
   vars[0] = avars[N];
   vals[0] = -cos(val);
   vars[1] = bvars[N];
   vals[1] = -sin(val);
   vars[2] = x3;
   vals[2] = alpha3;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, N);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, -alpha3*offset3, -alpha3*offset3,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIPdebugPrintCons(scip, lincons, NULL);
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++(*naddconss);

   for( i = 1; i <= N; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &avars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &bvars[i]) );
   }
   SCIPfreeBufferArray(scip, &avars);
   SCIPfreeBufferArray(scip, &bvars);

   return SCIP_OKAY;
}

/** adds the linear outer-approximation of Ben-Tal and Nemirovski for a SOC constraint of dimension 3
 *
 * Input is the data for a constraint \f$\sqrt{constant + (\alpha_1(x_1+offset1))^2 + (\alpha_2(x_2+offset2))^2)} \leq \alpha_3(x_3+offset3)\f$.
 * Here \f$\alpha3 > 0.0\f$, and the lower bound of \f$x_3 \geq -offset3\f$.
 * Also x2 = NULL is allowed, in which case the second term is assumed to be constant, and \f$offset2 \neq 0\f$ is needed.
 * */
static
SCIP_RETCODE presolveCreateBenTalNemirovskiApproxDim3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< original constraint */
   SCIP_VAR*             x1,                 /**< variable x1 */
   SCIP_VAR*             x2,                 /**< variable x2 */
   SCIP_VAR*             x3,                 /**< variable x3 */
   SCIP_Real             alpha1,             /**< coefficient of x1 */
   SCIP_Real             alpha2,             /**< coefficient of x2 */
   SCIP_Real             alpha3,             /**< coefficient of x3 */
   SCIP_Real             offset1,            /**< offset of x1 */
   SCIP_Real             offset2,            /**< offset of x2 */
   SCIP_Real             offset3,            /**< offset of x3 */
   int                   N,                  /**< size of linear approximation, need to be >= 1 */
   const char*           basename,           /**< string to use for building variable and constraint names */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   SCIP_CONS*     lincons;
   SCIP_VAR*      vars[3];
   SCIP_Real      vals[3];
   char           varname[255];
   char           linname[255];
   int            i;
   SCIP_VAR**     avars;
   SCIP_VAR**     bvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x1   != NULL);
   assert(x2   != NULL || !SCIPisZero(scip, offset2));
   assert(x3   != NULL);
   assert(SCIPisPositive(scip, alpha3));
   assert(SCIPisGE(scip, SCIPconsIsLocal(cons) ? SCIPvarGetLbLocal(x3) : SCIPvarGetLbGlobal(x3), -offset3));
   assert(basename != NULL);
   assert(N >= 1);
   assert(naddconss != NULL);

   SCIPdebugMsg(scip, "Creating linear Ben-Tal Nemirovski outer-approximation for <%s>.\n", basename);

   SCIP_CALL( SCIPallocBufferArray(scip, &avars, N+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bvars, N+1) );

   /* create additional variables */
   for( i = 0; i <= N; ++i )
   {
      (void) SCIPsnprintf(varname, 255, "soc#%s_a%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &avars[i], varname, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsLocal(cons), TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, avars[i]) );

      (void) SCIPsnprintf(varname, 255, "soc#%s_b%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &bvars[i], varname, 0.0, SCIPinfinity(scip), 0.0, 
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsLocal(cons), TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, bvars[i]) );
   }

   /* create first linear constraints - split into two because of the absolute value */
   vars[0] = avars[0];
   vals[0] = 1.0;
   vars[1] = x1;
   vals[1] = -alpha1;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, 0);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, alpha1 * offset1, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++(*naddconss);

   vars[0] = avars[0];
   vals[0] = 1.0;
   vars[1] = x1;
   vals[1] = alpha1;

   (void) SCIPsnprintf(linname, 255, "soc#%s#A%d", basename, 0);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha1 * offset1, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++(*naddconss);

   if( x2 != NULL )
   {
      vars[0] = bvars[0];
      vals[0] = 1.0;
      vars[1] = x2;
      vals[1] = -alpha2;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, alpha2 * offset2, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);

      vars[0] = bvars[0];
      vals[0] = 1.0;
      vars[1] = x2;
      vals[1] = alpha2;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha2 * offset2, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);
   }
   else
   { /* second summand is just a constant */
      if( SCIPconsIsLocal(cons) )
      {
         SCIP_CALL( SCIPchgVarLbNode(scip, NULL, bvars[0], ABS(alpha2 * offset2)) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarLbGlobal(scip, bvars[0], ABS(alpha2 * offset2)) );
      }
   }

   /* create intermediate linear constraints */
   for( i = 1; i <= N; ++i )
   {
      SCIP_Real val;

      val = M_PI / pow(2.0, (double) (i+1));

      vars[0] = avars[i-1];
      vals[0] = cos(val);
      vars[1] = bvars[i-1];
      vals[1] = sin(val);
      vars[2] = avars[i];
      vals[2] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);

      vars[0] = avars[i-1];
      vals[0] = sin(val);
      vars[1] = bvars[i-1];
      vals[1] = -cos(val);
      vars[2] = bvars[i];
      vals[2] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);

      vars[0] = avars[i-1];
      vals[0] = -sin(val);
      vars[1] = bvars[i-1];
      vals[1] = cos(val);
      vars[2] = bvars[i];
      vals[2] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++(*naddconss);
   }

   /* create last linear constraints */
   vars[0] = x3;
   vals[0] = alpha3;
   vars[1] = avars[N];
   vals[1] = -1.0;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, N);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha3 * offset3, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++(*naddconss);

   vars[0] = avars[N];
   vals[0] = tan( M_PI / pow(2.0, (double) (N+1)) );
   vars[1] = bvars[N];
   vals[1] = -1.0;

   (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, i);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, 0.0, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++(*naddconss);

   for( i = 0; i <= N; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &avars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &bvars[i]) );
   }
   SCIPfreeBufferArray(scip, &avars);
   SCIPfreeBufferArray(scip, &bvars);

   return SCIP_OKAY;
}

/** adds a linear outer approx for a three dimensional SOC constraint
 * 
 * chooses between Ben-Tan/Nemirovski and Glineur and calls the corresponding function
 */
static
SCIP_RETCODE presolveCreateOuterApproxDim3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< original constraint */
   SCIP_VAR*             x1,                 /**< variable x1 */
   SCIP_VAR*             x2,                 /**< variable x2 */
   SCIP_VAR*             x3,                 /**< variable x3 */
   SCIP_Real             alpha1,             /**< coefficient of x1 */
   SCIP_Real             alpha2,             /**< coefficient of x2 */
   SCIP_Real             alpha3,             /**< coefficient of x3 */
   SCIP_Real             offset1,            /**< offset of x1 */
   SCIP_Real             offset2,            /**< offset of x2 */
   SCIP_Real             offset3,            /**< offset of x3 */
   int                   N,                  /**< size of linear approximation, need to be >= 1 */
   SCIP_Bool             glineur,            /**< whether to prefer Glineur to Ben-Tal Nemirovski */
   const char*           basename,           /**< string to use for building variable and constraint names */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   if( glineur )
   {
      SCIP_CALL( presolveCreateGlineurApproxDim3(scip, cons, x1, x2, x3, alpha1, alpha2, alpha3, offset1, offset2, offset3, N, basename, naddconss) );
   }
   else
   {
      SCIP_CALL( presolveCreateBenTalNemirovskiApproxDim3(scip, cons, x1, x2, x3, alpha1, alpha2, alpha3, offset1, offset2, offset3, N, basename, naddconss) );
   }

   return SCIP_OKAY;
}

/** adds linear outer approximation of Ben-Tal and Nemirovski for a constraint \f$\gamma + \sum_{i=1}^n (\alpha_i (x_i + \beta_i))^2 \leq (\alpha_{n+1} (x_{n+1} + \beta_{n+1}))^2\f$ to the LP
 *
 * if n > 2, calls same function recursively;
 * if n = 2, calls presolveCreateBenTalNemirovskiApproxDim3
 */
static
SCIP_RETCODE presolveCreateOuterApprox(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nlhsvars,           /**< number of variables on left hand side (n) */
   SCIP_VAR**            lhsvars,            /**< variables on left hand side (x_i) */
   SCIP_Real*            lhscoefs,           /**< coefficients of variables on left hand side (alpha_i) */
   SCIP_Real*            lhsoffsets,         /**< offsets of variable on left hand side (beta_i) */
   SCIP_VAR*             rhsvar,             /**< variable on right hand side (y) */
   SCIP_Real             rhscoeff,           /**< coefficient of variable on right hand side (alpha_{n+1}) */
   SCIP_Real             rhsoffset,          /**< offset of variable on right hand side (beta_{n+1}) */
   SCIP_Real             constant,           /**< constant term (gamma) */
   const char*           basename,           /**< prefix for variable and constraint name */
   SCIP_CONS*            origcons,           /**< original constraint for which this SOC3 set is added */
   int                   soc3_nr_auxvars,    /**< number of auxiliary variables to use for a SOC3 constraint, or 0 if automatic */
   SCIP_Bool             glineur,            /**< whether Glineur should be preferred to Ben-Tal Nemirovski */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   char       name[255];
   SCIP_VAR*  auxvar1;
   SCIP_VAR*  auxvar2;

   assert(scip     != NULL);
   assert(lhsvars  != NULL);
   assert(nlhsvars >= 1);
   assert(lhscoefs != NULL);
   assert(lhsoffsets != NULL);
   assert(rhsvar   != NULL);
   assert(basename != NULL);
   assert(!SCIPisNegative(scip, constant));
   assert(naddconss != NULL);

   if( nlhsvars == 1 )
   { /* end of recursion */
      assert(SCIPisPositive(scip, constant));
      SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
            lhsvars[0],    NULL,           rhsvar,
            lhscoefs[0],   1.0,            rhscoeff,
            lhsoffsets[0], sqrt(constant), rhsoffset,
            soc3_nr_auxvars, glineur, basename, naddconss) );

      return SCIP_OKAY;
   }

   if( nlhsvars == 2 && SCIPisZero(scip, constant) )
   { /* end of recursion */
      assert(lhsvars[0] != NULL);
      assert(lhsvars[1] != NULL);
      assert(rhsvar     != NULL);
      SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
            lhsvars[0],    lhsvars[1],    rhsvar,
            lhscoefs[0],   lhscoefs[1],   rhscoeff,
            lhsoffsets[0], lhsoffsets[1], rhsoffset,
            soc3_nr_auxvars, glineur, basename, naddconss) );

      return SCIP_OKAY;
   }

   if( nlhsvars == 3 || (nlhsvars == 2 && !SCIPisZero(scip, constant)) )
   { 
      /* a bit special case too */
      /* for first two variables on lhs, create a new aux.var and a new SOC3 */
      (void) SCIPsnprintf(name, 255, "%s#z1", basename);
      SCIP_CALL( SCIPcreateVar(scip, &auxvar1, name, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar1) );

      /* constraint alpha_0 (x_0+beta0)^2 + alpha_1 (x_1+beta1)^2 <= auxvar^2 */
      SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
            lhsvars[0],    lhsvars[1],    auxvar1,
            lhscoefs[0],   lhscoefs[1],   1.0,
            lhsoffsets[0], lhsoffsets[1], 0.0,
            soc3_nr_auxvars, glineur, name, naddconss) );

      (void) SCIPsnprintf(name, 255, "%s_soc3", basename);
      if( nlhsvars == 3 )
      { /* create new constraint alpha_2 (x_2+beta2)^2 + auxvar^2 <= (rhscoeff * (rhsvar+rhsoffset))^2 */
         SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
               lhsvars[2],    auxvar1, rhsvar,
               lhscoefs[2],   1.0,     rhscoeff,
               lhsoffsets[2], 0.0,     rhsoffset,
               soc3_nr_auxvars, glineur, name, naddconss) );
      }
      else
      { /* create new constraint auxvar^2 + sqrt(constant)^2 <= (rhscoeff * (rhsvar+rhsoffset))^2 */
         SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
               auxvar1, NULL,           rhsvar,
               1.0,     1.0,            rhscoeff,
               0.0,     sqrt(constant), rhsoffset,
               soc3_nr_auxvars, glineur, name, naddconss) );
      }

      SCIP_CALL( SCIPreleaseVar(scip, &auxvar1) );

      return SCIP_OKAY;
   }

   /* nlhsvars >= 4 */

   (void) SCIPsnprintf(name, 255, "%s#z1", basename);
   SCIP_CALL( SCIPcreateVar(scip, &auxvar1, name, 0.0, SCIPinfinity(scip), 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar1) );

   /* approx for left half of lhs */
   SCIP_CALL( presolveCreateOuterApprox(scip,
         nlhsvars/2, lhsvars, lhscoefs, lhsoffsets,
         auxvar1, 1.0, 0.0,
         constant, name, origcons, soc3_nr_auxvars, glineur, naddconss) );

   (void) SCIPsnprintf(name, 255, "%s#z2", basename);
   SCIP_CALL( SCIPcreateVar(scip, &auxvar2, name, 0., SCIPinfinity(scip), 0.0, 
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar2) );

   /* approx for right half of lhs */
   SCIP_CALL( presolveCreateOuterApprox(scip,
         nlhsvars-nlhsvars/2, &lhsvars[nlhsvars/2], &lhscoefs[nlhsvars/2], &lhsoffsets[nlhsvars/2],
         auxvar2, 1.0, 0.0,
         0.0, name, origcons, soc3_nr_auxvars, glineur, naddconss) );

   /* SOC constraint binding both auxvar's */
   (void)SCIPsnprintf(name, 255, "%s_soc3", basename);
   SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
         auxvar1, auxvar2, rhsvar,
         1.0,     1.0,     rhscoeff,
         0.0,     0.0,     rhsoffset,
         soc3_nr_auxvars, glineur, name, naddconss) );

   SCIP_CALL( SCIPreleaseVar(scip, &auxvar1) );
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar2) );

   return SCIP_OKAY;
}

/** propagates variable bounds */
static
SCIP_RETCODE propagateBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_RESULT*          result,             /**< buffer to store result of propagation */
   int*                  nchgbds,            /**< buffer where to add number of tightened bounds */
   SCIP_Bool*            redundant           /**< buffer to indicate whether constraint was marked for deletion because of redundancy */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL  lhsrange;
   SCIP_INTERVAL* lhsranges;
   SCIP_INTERVAL  rhsrange;
   SCIP_INTERVAL  a, b, c;
   SCIP_ROUNDMODE roundmode;
   SCIP_Bool      infeas, tightened;
   int            i;
   SCIP_Real      lb, ub;

   assert(scip   != NULL);
   assert(cons   != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(redundant != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *redundant = FALSE;

   if( !SCIPconsIsMarkedPropagate(cons) )
   {
      SCIPdebugMsg(scip, "skip propagation for constraint %s\n", SCIPconsGetName(cons));
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   else
   {
      SCIPdebugMsg(scip, "try propagation for constraint %s\n", SCIPconsGetName(cons));
   }

   *result = SCIP_DIDNOTFIND;
   SCIP_CALL( SCIPunmarkConsPropagate(scip, cons) );

   /* @todo do something clever to decide whether propagation should be tried */

   /* compute activity on lhs */
   SCIPintervalSet(&lhsrange, consdata->constant);
   SCIP_CALL( SCIPallocBufferArray(scip, &lhsranges, consdata->nvars) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      lb = SCIPcomputeVarLbLocal(scip, consdata->vars[i]) - SCIPepsilon(scip);
      ub = SCIPcomputeVarUbLocal(scip, consdata->vars[i]) + SCIPepsilon(scip);
      SCIPintervalSetBounds(&lhsranges[i], MIN(lb, ub), MAX(lb, ub));
      if( consdata->offsets[i] != 0.0 )
         SCIPintervalAddScalar(SCIPinfinity(scip), &lhsranges[i], lhsranges[i], consdata->offsets[i]);
      if( consdata->coefs[i]   != 1.0 )
         SCIPintervalMulScalar(SCIPinfinity(scip), &lhsranges[i], lhsranges[i], consdata->coefs[i]);
      SCIPintervalSquare(SCIPinfinity(scip), &lhsranges[i], lhsranges[i]);

      SCIPintervalAdd(SCIPinfinity(scip), &lhsrange, lhsrange, lhsranges[i]);
   }

   /* compute activity on rhs: rhsrange = sqr(rhscoeff * (rhsvar + rhsoffset) ) */
   lb = SCIPcomputeVarLbLocal(scip, consdata->rhsvar) - SCIPepsilon(scip);
   ub = SCIPcomputeVarUbLocal(scip, consdata->rhsvar) + SCIPepsilon(scip);
   SCIPintervalSetBounds(&rhsrange, MIN(lb, ub), MAX(lb, ub));

   if( consdata->rhsoffset != 0.0 )
      SCIPintervalAddScalar(SCIPinfinity(scip), &rhsrange, rhsrange, consdata->rhsoffset);
   if( consdata->rhscoeff  != 1.0 )
      SCIPintervalMulScalar(SCIPinfinity(scip), &rhsrange, rhsrange, consdata->rhscoeff);
   SCIPintervalSquare(SCIPinfinity(scip), &rhsrange, rhsrange);

   /* check for infeasibility */
   if( SCIPisGT(scip, lhsrange.inf-SCIPfeastol(scip), rhsrange.sup) )
   {
      SCIPdebugMsg(scip, "propagation found constraint <%s> infeasible: lhs = [%.15g,%.15g]-feastol-eps > rhs = [%.15g,%.15g]\n",
         SCIPconsGetName(cons), lhsrange.inf, lhsrange.sup, rhsrange.inf, rhsrange.sup);
      *result = SCIP_CUTOFF;
      goto TERMINATE;
   }

   /* check for redundancy: max(lhsrange) <= min(rhsrange) */
   if( SCIPisLE(scip, lhsrange.sup, rhsrange.inf) )
   {
      SCIPdebugMsg(scip, "propagation found constraint <%s> redundant: lhs = [%.15g,%.15g] <= rhs = [%.15g,%.15g]\n",
         SCIPconsGetName(cons), lhsrange.inf, lhsrange.sup, rhsrange.inf, rhsrange.sup);
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      goto TERMINATE;
   }

   /* try to tighten variable on rhs */
   if( SCIPvarGetStatus(consdata->rhsvar) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPintervalSquareRoot(SCIPinfinity(scip), &a, lhsrange);
      if( consdata->rhscoeff != 1.0 )
         SCIPintervalDivScalar(SCIPinfinity(scip), &a, a, consdata->rhscoeff);
      if( consdata->rhsoffset != 0.0 )
         SCIPintervalSubScalar(SCIPinfinity(scip), &a, a, consdata->rhsoffset);
      SCIP_CALL( SCIPtightenVarLb(scip, consdata->rhsvar, SCIPintervalGetInf(a), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMsg(scip, "propagation found constraint <%s> infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         goto TERMINATE;
      }
      if( tightened )
      {
         SCIPdebugMsg(scip, "propagation tightened bounds of rhs variable <%s> in constraint <%s>\n", SCIPvarGetName(consdata->rhsvar), SCIPconsGetName(cons));
         *result = SCIP_REDUCEDDOM;
         ++*nchgbds;
      }
   }

   /* try to tighten variables on lhs */
   SCIPintervalSub(SCIPinfinity(scip), &b, rhsrange, lhsrange);  /*lint !e644 */
   for( i = 0; i < consdata->nvars; ++i )
   {
      if( SCIPvarGetStatus(consdata->vars[i]) == SCIP_VARSTATUS_MULTAGGR )
         continue;

      roundmode = SCIPintervalGetRoundingMode();
      if( !SCIPisInfinity(scip, b.sup) )
      {
         SCIPintervalSetRoundingModeUpwards();
         a.sup = b.sup + lhsranges[i].inf;
      }
      else
      {
         a.sup = SCIPinfinity(scip);
      }
      if( !SCIPisInfinity(scip, -b.inf) )
      {
         SCIPintervalSetRoundingModeDownwards();
         a.inf = b.inf + lhsranges[i].sup;
      }
      else
      {
         a.inf = -SCIPinfinity(scip);
      }
      SCIPintervalSetRoundingMode(roundmode);
      SCIPintervalSquareRoot(SCIPinfinity(scip), &a, a);

      assert(consdata->coefs[i] >= 0.0); /* should be ensured in create and presolveRemoveFixed */

      c = a;
      if( consdata->coefs[i]   != 1.0 )
         SCIPintervalDivScalar(SCIPinfinity(scip), &c, c, consdata->coefs[i]);
      if( consdata->offsets[i] != 0.0 )
         SCIPintervalSubScalar(SCIPinfinity(scip), &c, c, consdata->offsets[i]);

      SCIP_CALL( SCIPtightenVarUb(scip, consdata->vars[i], SCIPintervalGetSup(c), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMsg(scip, "propagation found constraint <%s> infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         goto TERMINATE;
      }
      if( tightened )
      {
         SCIPdebugMsg(scip, "propagation tightened bounds of lhs variable <%s> in constraint <%s>\n", SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons));
         *result = SCIP_REDUCEDDOM;
         ++*nchgbds;
      }

      c = a;
      SCIPintervalDivScalar(SCIPinfinity(scip), &c, c, -consdata->coefs[i]);
      if( consdata->offsets[i] != 0.0 )
         SCIPintervalSubScalar(SCIPinfinity(scip), &c, c, consdata->offsets[i]);

      SCIP_CALL( SCIPtightenVarLb(scip, consdata->vars[i], SCIPintervalGetInf(c), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMsg(scip, "propagation found constraint <%s> infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
         goto TERMINATE;
      }
      if( tightened )
      {
         SCIPdebugMsg(scip, "propagation tightened bounds of lhs variable <%s> in constraint <%s>\n", SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons));
         *result = SCIP_REDUCEDDOM;
         ++*nchgbds;
      }
   }

TERMINATE:
   SCIPfreeBufferArray(scip, &lhsranges);

   if( *result != SCIP_DIDNOTFIND )
   {
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** tries to adjust a solution such that it satisfies a given constraint by increasing the value for the constraints right hand side variable */
static
SCIP_RETCODE polishSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to polish */
   SCIP_Bool*            success             /**< buffer to store whether polishing was successful */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real rhsval;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sol  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, consdata->rhscoeff));

   /* compute minimal rhs variable value so that constraint is satisfied */
   if( !SCIPisInfinity(scip, consdata->lhsval) )
      rhsval = consdata->lhsval / consdata->rhscoeff - consdata->rhsoffset;
   else
      rhsval = consdata->rhscoeff > 0.0 ? SCIPinfinity(scip) : -SCIPinfinity(scip);

   if( consdata->rhscoeff > 0.0 )
   {
      assert(SCIPvarMayRoundUp(consdata->rhsvar));

      /* round rhsval up, if variable is integral */
      if( SCIPvarIsIntegral(consdata->rhsvar) && !SCIPisInfinity(scip, rhsval) )
         rhsval = SCIPceil(scip, rhsval);

      /* if new value is above upper bound, we are lost */
      if( SCIPisGT(scip, rhsval, SCIPvarGetUbGlobal(consdata->rhsvar)) )
      {
         *success = FALSE;
      }
      else
      {
         /* if new value is larger then current one, increase to new value */
         if( rhsval > SCIPgetSolVal(scip, sol, consdata->rhsvar) )
         {
            SCIPdebugMsg(scip, "increase <%s> to %g\n", SCIPvarGetName(consdata->rhsvar), rhsval);
            SCIP_CALL( SCIPsetSolVal(scip, sol, consdata->rhsvar, rhsval) );
         }

         *success = TRUE;
      }
   }
   else
   {
      assert(SCIPvarMayRoundDown(consdata->rhsvar));

      /* round rhsval down, if variable is integral */
      if( SCIPvarIsIntegral(consdata->rhsvar) )
         rhsval = SCIPfloor(scip, rhsval);

      /* if new value is below lower bound, we are lost */
      if( SCIPisLT(scip, rhsval, SCIPvarGetLbGlobal(consdata->rhsvar)) )
      {
         *success = FALSE;
      }
      else
      {
         /* if new value is below current one, decrease to new value */
         if( rhsval < SCIPgetSolVal(scip, sol, consdata->rhsvar) )
         {
            SCIPdebugMsg(scip, "decrease <%s> to %g\n", SCIPvarGetName(consdata->rhsvar), rhsval);
            SCIP_CALL( SCIPsetSolVal(scip, sol, consdata->rhsvar, rhsval) );
         }

         *success = TRUE;
      }
   }

   SCIPdebugMsg(scip, "polishing solution for constraint <%s> was %ssuccessful\n", SCIPconsGetName(cons), *success ? "" : "not ");

   return SCIP_OKAY;
}

/** disaggregates a (sufficiently large) SOC constraint into smaller ones; for each term on the lhs we add a quadratic
 *  constraint \f$(\alpha_i * (x_i + \beta_i))^2 \leq \alpha_{n+1} (x_{n+1} + \beta_{n+1})\, z_i\f$ and a single linear constraint
 *  \f$\sum_i z_i \leq \alpha_{n+1}\, (x_{n+1} + \beta_{n+1})\f$; each quadratic constraint might be upgraded to a SOC; since the
 *  violations of all quadratic constraints sum up we scale each constraint by the number of lhs terms + 1
 *
 *  @todo if rhsvar is NULL, then the disaggregation does not produce further cones. Should it then be upgraded
 *  to a quadratic and let the quadratic desaggregate it?
 *  The code assumes now that the rhsvar is not NULL in order build the direct SOC -> SOC disaggregation
 */
static
SCIP_RETCODE disaggregate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int*                  naddconss,          /**< pointer to count total number of added constraints */
   int*                  ndelconss,          /**< pointer to count total number of deleted constraints */
   SCIP_Bool*            success             /**< pointer to store whether disaggregation was successful */
   )
{
   SCIP_CONS* discons;
   SCIP_VAR** disvars;
   SCIP_VAR** sumvars;
   SCIP_VAR** difvars;
   SCIP_Real* discoefs;
   SCIP_VAR*  lhsvars[2];
   SCIP_VAR*  aggvars[2];
   SCIP_Real  coefs[2];
   SCIP_Real  offsets[2];
   SCIP_Real  scalars[2];
   char name[SCIP_MAXSTRLEN];
   SCIP_Real constant;
   SCIP_Real scale;
   SCIP_Bool infeas;
   int ndisvars;
   int i;

   assert(naddconss != NULL);
   assert(ndelconss != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* disaggregation does not make much sense if there are too few variables */
   if( consdata->nvars < 3 )
   {
      SCIPdebugMsg(scip, "can not disaggregate too small soc constraint %s\n", SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   if( consdata->rhsvar == NULL )
   {
      SCIPdebugMsg(scip, "can not disaggregate directly into a soc without rhs var %s\n", SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   /* there are at most n + 2 many linear varibles */
   SCIP_CALL( SCIPallocBufferArray(scip, &disvars, consdata->nvars + 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sumvars, consdata->nvars + 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &difvars, consdata->nvars + 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &discoefs, consdata->nvars + 2) );
   ndisvars = 0;

   scale = 1.0 * (consdata->nvars + 1)/4.0;

   /* add (*) (alpha_i * (x_i + beta_i))^2 <= alpha_{n+1} * (x_{n+1} + beta_{n+1}) * z_i:
    * create sumvar = alpha_{n+1} * (x_{n+1} + beta_{n+1}) + z_i (multiagg)
    * create difvar = alpha_{n+1} * (x_{n+1} + beta_{n+1}) - z_i (multiagg)
    * note that (*) is equiv to sqrt( (2 * alpha_i * (x_i + beta_i))^2 + difvar^2) <= sumvar
    * scaling give us: sqrt( (2 * scale * alpha_i * (x_i + beta_i))^2 + (scale * difvar)^2) <= scale * sumvar
    */
   aggvars[0] = consdata->rhsvar;
   scalars[0] = consdata->rhscoeff;
   for( i = 0; i < consdata->nvars; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%s_%d", SCIPvarGetName(consdata->vars[i]), i);
      SCIP_CALL( SCIPcreateVar(scip, &disvars[i], name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE,
            NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, disvars[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedisS_%s_%d", SCIPvarGetName(consdata->vars[i]), i);
      SCIP_CALL( SCIPcreateVar(scip, &sumvars[i], name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE,
            NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, sumvars[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedisD_%s_%d", SCIPvarGetName(consdata->vars[i]), i);
      SCIP_CALL( SCIPcreateVar(scip, &difvars[i], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, difvars[i]) );

      aggvars[1] = disvars[i];
      scalars[1] = 1.0;
      constant    = consdata->rhscoeff * consdata->rhsoffset;
      SCIP_CALL( SCIPmultiaggregateVar(scip, sumvars[i], 2, aggvars, scalars, constant, &infeas, success) );
      /* @todo what shall we do if multiagg fails? */
      assert(!infeas && *success);

      scalars[1] = -1.0;
      SCIP_CALL( SCIPmultiaggregateVar(scip, difvars[i], 2, aggvars, scalars, constant, &infeas, success) );
      assert(!infeas && *success);

      /* create soc */
      lhsvars[0] = difvars[i];
      coefs[0]   = scale;
      offsets[0] = 0.0;
      lhsvars[1] = consdata->vars[i];
      coefs[1]   = scale * 2 * consdata->coefs[i];
      offsets[1] = consdata->offsets[i];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "consdis_%s_%d", SCIPconsGetName(cons), i);
      SCIP_CALL( SCIPcreateConsBasicSOC(scip, &discons, name, 2, lhsvars, coefs, offsets, 0.0, sumvars[i], scale, 0.0) );
      SCIP_CALL( SCIPaddCons(scip, discons) );
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintCons(scip, discons, NULL) );
#endif
      SCIP_CALL( SCIPreleaseCons(scip, &discons) );
      ++(*naddconss);

      /* linear coefficient in the linear constraint */
      discoefs[ndisvars] = 1.0;
      ++ndisvars;
   }
   assert(ndisvars == consdata->nvars);

   /* add gamma <= alpha_{n+1} * (x_{n+1} + beta_{n+1}) * z_i
    * sumvar and difvar are the same as before, but the equivalent soc now is
    * sqrt(4 * gamma + difvar^2) <= sumvar
    * scaling give us: sqrt( (4 * scale^2 * gamma + (scale * difvar)^2) <= scale * sumvar
    */
   if( !SCIPisZero(scip, consdata->constant) )
   {

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_const_%s", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &disvars[ndisvars], name, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, disvars[ndisvars]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedisS_const_%s", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &sumvars[ndisvars], name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE,
            NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, sumvars[ndisvars]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedisD_const_%s", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &difvars[ndisvars], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, difvars[ndisvars]) );

      aggvars[1] = disvars[i];
      scalars[1] = 1.0;
      constant   = consdata->rhscoeff * consdata->rhsoffset;
      SCIP_CALL( SCIPmultiaggregateVar(scip, sumvars[i], 2, aggvars, scalars, constant, &infeas, success) );
      assert(!infeas && *success);

      scalars[1] = -1.0;
      SCIP_CALL( SCIPmultiaggregateVar(scip, difvars[i], 2, aggvars, scalars, constant, &infeas, success) );
      assert(!infeas && *success);

      /* create soc */
      lhsvars[0] = difvars[ndisvars];
      coefs[0]   = scale;
      offsets[0] = 0.0;
      constant   = 4.0 * SQR(scale) * consdata->constant;
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "consdis_%s_constant", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateConsBasicSOC(scip, &discons, name, 1, lhsvars, coefs, offsets, constant,
               sumvars[ndisvars], scale, 0.0) );
      SCIP_CALL( SCIPaddCons(scip, discons) );
      SCIP_CALL( SCIPreleaseCons(scip, &discons) );
      ++(*naddconss);

      /* linear coefficient in the linear constraint */
      discoefs[ndisvars] = 1.0;
      ++ndisvars;
   }

   /* create linear constraint sum z_i <= alpha_{n+1} * (x_{n+1} + beta_{n+1}); first add extra coefficient for the rhs */
   discoefs[ndisvars] = -1.0 * consdata->rhscoeff;
   disvars[ndisvars] = consdata->rhsvar;

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "consdis_linear_%s", SCIPconsGetName(cons));
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &discons, name, ndisvars + 1, disvars, discoefs, -SCIPinfinity(scip),
            consdata->rhscoeff * consdata->rhsoffset) );

   SCIP_CALL( SCIPaddCons(scip, discons) );
   SCIP_CALL( SCIPreleaseCons(scip, &discons) );
   ++(*naddconss);

   /* release all variables */
   for( i = ndisvars - 1; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &disvars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &sumvars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &difvars[i]) );
   }
   SCIPfreeBufferArray(scip, &discoefs);
   SCIPfreeBufferArray(scip, &difvars);
   SCIPfreeBufferArray(scip, &sumvars);
   SCIPfreeBufferArray(scip, &disvars);

   /* delete constraint */
   SCIP_CALL( SCIPdelCons(scip, cons) );
   ++(*ndelconss);

   *success = TRUE;

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
   SCIP_CONSDATA*     consdata;
   SCIP_CONS*         maxviolcons;
   SCIP_Bool          success;
   SCIP_Bool          cutoff;
   SCIP_Bool          redundant;
   int                nbndchg;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result   != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, sol, &maxviolcons) );

   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   /* if we are above the 100'th enforcement round for this node, something is strange
    * (maybe the LP does not think that the cuts we add are violated, or we do ECP on a high-dimensional convex function)
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
      else
         ++conshdlrdata->nenforounds;
   }
   else
   {
      conshdlrdata->lastenfonode = SCIPgetCurrentNode(scip);
      conshdlrdata->nenforounds = 0;
   }

   /* try separation, this should usually work */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, TRUE, &cutoff, &success) );
   if( cutoff )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }
   if( success )
   {
      SCIPdebugMsg(scip, "enforced by separation\n");
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* try propagation */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      if( !SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) )
         continue;

      nbndchg = 0;
      SCIP_CALL( propagateBounds(scip, conss[c], result, &nbndchg, &redundant) );  /*lint !e613*/
      assert(!redundant); /* constraint should not be violated and redundant simultaneously (unless solution is far out of bounds) */
      if( *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM )
      {
         SCIPdebugMsg(scip, "enforced by %s\n", *result == SCIP_CUTOFF ? "cutting off node" : "reducing domain");
         return SCIP_OKAY;
      }
   }

   SCIPwarningMessage(scip, "could not enforce feasibility by separating or branching; declaring solution with viol %g feasible\n", SCIPconsGetData(maxviolcons)->violation);
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/*
 * Quadratic constraint upgrading
 */


#ifdef QUADCONSUPGD_PRIORITY
/** tries to upgrade a quadratic constraint to a SOC constraint
 *
 *  We treat as special cases, quadratic with no bilinear terms and hyperbolic quadratic
 *  constraints with exactly on bilinear component containing nonnegative variables. For this we use the formula:
 *  \f[
 *    x^T x \leq yz,\; y \geq 0,\; z \geq 0
 *    \qquad\Leftrightarrow\qquad
 *    \left\| \left(\begin{array}{c} x \\ \frac{1}{2}(y - z)\end{array}\right) \right\| \leq \frac{1}{2}(y + z).
 *  \f]
 *
 * @todo implement more general hyperbolic upgrade, e.g., for -x^T x + yz >= 0 or x^T x <= ax + by + cyz
 */
static
SCIP_DECL_QUADCONSUPGD(upgradeConsQuadratic)
{
   int            nquadvars;
   int            nbilinterms;
   SCIP_QUADVARTERM* term;
   SCIP_QUADVARTERM* quadterms;
   SCIP_BILINTERM* bilinterm;
   SCIP_VAR**     quadvars;
   SCIP_VAR**     lhsvars;
   SCIP_Real*     lhscoefs;
   SCIP_Real*     lhsoffsets;
   SCIP_Real      lhsconstant;
   int            lhscount;
   int            lhsnvars;
   SCIP_VAR*      rhsvar; 
   SCIP_Real      rhscoef;
   SCIP_Real      rhsoffset;
   SCIP_VAR*      bilinvar1;
   SCIP_VAR*      bilinvar2;
   SCIP_Real      bilincoef;
   int            i;
   int            j;
   int            negeigpos;
   SCIP_Real*     a;
   SCIP_Real*     bp;
   SCIP_Real*     eigvals;
   SCIP_Bool      infeas;
   SCIP_Bool      success;
   SCIP_Bool      trygeneralupg;
   int            nneg;
   int            npos;
   SCIP_Bool      rhsvarfound;
   SCIP_Bool      rhsissoc;
   SCIP_Bool      lhsissoc;
   SCIP_Real      rhsvarlb;
   SCIP_Real      rhsvarub;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nupgdconss != NULL);
   assert(upgdconss  != NULL);

   *nupgdconss = 0;

   SCIPdebugMsg(scip, "upgradeConsQuadratic called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   /* currently do not support linear parts in upgrading of SOC constraints */
   if( SCIPgetNLinearVarsQuadratic(scip, cons) > 0 )
      return SCIP_OKAY;

   nbilinterms = SCIPgetNBilinTermsQuadratic(scip, cons);
   nquadvars = SCIPgetNQuadVarTermsQuadratic(scip, cons);

   /* currently, a proper SOC constraint needs at least 3 variables */
   if( nquadvars < 3 )
      return SCIP_OKAY;

   /* reserve space */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhsvars,    nquadvars - 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhscoefs,   nquadvars - 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhsoffsets, nquadvars - 1) );

   /* initialize data */
   a = NULL;
   bp = NULL;
   eigvals = NULL;
   quadvars = NULL;
   bilinvar1 = NULL;
   bilinvar2 = NULL;
   lhsconstant = 0.0;
   lhscount = 0;
   rhsvar = NULL;
   rhscoef = 0.0;
   rhsoffset = 0.0;

   trygeneralupg = FALSE;

   /* if more than one bilinear term is present, we are in the general case */
   if( nbilinterms > 1 )
   {
      trygeneralupg = TRUE;
      goto GENERALUPG;
   }

   /* check hyperbolic part */
   if ( nbilinterms == 1 )
   {
      bilinterm = SCIPgetBilinTermsQuadratic(scip, cons);
      bilinvar1 = bilinterm->var1;
      bilinvar2 = bilinterm->var2;
      bilincoef = bilinterm->coef;

      /* the variables in the bilinear term need to be nonnegative */
      if ( SCIPisNegative(scip, SCIPvarGetLbGlobal(bilinvar1)) || SCIPisNegative(scip, SCIPvarGetLbGlobal(bilinvar2)) )
      {
         trygeneralupg = TRUE;
         goto GENERALUPG;
      }

      /* we need a rhs */
      if ( ! SCIPisZero(scip, SCIPgetRhsQuadratic(scip, cons)) )
      {
         trygeneralupg = TRUE;
         goto GENERALUPG;
      }

      /* we only allow for -1.0 bilincoef */
      if ( ! SCIPisEQ(scip, bilincoef, -1.0) )
      {
         trygeneralupg = TRUE;
         goto GENERALUPG;
      }

      /* check that bilinear terms do not appear in the rest and quadratic terms have postive sqrcoef have no lincoef */
      quadterms = SCIPgetQuadVarTermsQuadratic(scip, cons);
      for (i = 0; i < nquadvars; ++i)
      {
         term = &quadterms[i];

         if( ! SCIPisZero(scip, term->lincoef) || SCIPisNegative(scip, term->sqrcoef) )
         {
            trygeneralupg = TRUE;
            goto GENERALUPG;
         }

         if ( (term->var == bilinvar1 || term->var == bilinvar2) && ! SCIPisZero(scip, term->sqrcoef) )
         {
            trygeneralupg = TRUE;
            goto GENERALUPG;
         }
      }
   }


   quadterms = SCIPgetQuadVarTermsQuadratic(scip, cons);
   assert( quadterms != NULL );

   if( ! SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons)) )
   {
      /* try whether constraint on right hand side is SOC */
      lhsconstant = -SCIPgetRhsQuadratic(scip, cons);

      for( i = 0; i < nquadvars; ++i )
      {
         term = &quadterms[i];

         /* skip terms with 0 coefficients */
         if ( SCIPisZero(scip, term->sqrcoef) )
            continue;

         if( term->sqrcoef > 0.0 )
         {
            if( lhscount >= nquadvars - 1 )
            { /* too many variables on lhs, i.e., all variables seem to have positive coefficient */
               rhsvar = NULL;
               break;
            }

            lhsvars[lhscount]  = term->var;
            lhscoefs[lhscount] = sqrt(term->sqrcoef);

            if( term->lincoef != 0.0 )
            {
               lhsoffsets[lhscount] = term->lincoef / (2 * term->sqrcoef);
               lhsconstant -= term->lincoef * term->lincoef / (4 * term->sqrcoef);
            }
            else
            {
               lhsoffsets[lhscount] = 0.0;
            }

            ++lhscount;
         }
         else if( rhsvar != NULL ||
               (SCIPisLT(scip, SCIPcomputeVarLbGlobal(scip, term->var), -term->lincoef / (2 * term->sqrcoef))
                && SCIPisGT(scip, SCIPcomputeVarUbGlobal(scip, term->var), -term->lincoef / (2 * term->sqrcoef))) )
         { /* second variable with negative coefficient -> cannot be SOC */
            /* or rhs side changes sign */
            rhsvar = NULL;
            break;
         }
         else
         {
            rhsvar       = term->var;
            rhsoffset    = term->lincoef / (2 * term->sqrcoef);
            lhsconstant -= term->lincoef * term->lincoef / (4 * term->sqrcoef);

            if( SCIPisGE(scip, SCIPcomputeVarLbGlobal(scip, term->var), -term->lincoef / (2*term->sqrcoef)) )
               rhscoef = sqrt(-term->sqrcoef);
            else
               rhscoef = -sqrt(-term->sqrcoef);
         }
      }
   }

   /* treat hyberbolic case */
   if ( nbilinterms == 1 )
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR* auxvarsum;
      SCIP_VAR* auxvardiff;
      SCIP_CONS* couplingcons;
      SCIP_VAR* consvars[3];
      SCIP_Real consvals[3];

      /* can only upgrade if rhs is 0 */
      if ( rhsvar != NULL )
         goto cleanup;

      SCIPdebugMsg(scip, "found hyberbolic quadratic constraint <%s> to be SOC\n", SCIPconsGetName(cons));

      /* check if upgdconss is long enough to store upgrade constraints: we need two if we will have a quadratic
       * constraint for the left hand side left */
      *nupgdconss = SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) ? 1 : 2;
      if ( *nupgdconss > upgdconsssize )
      {
         /* signal that we need more memory and return */
         *nupgdconss = -*nupgdconss;
         goto cleanup;
      }

      /* create auxiliary variable for sum (nonnegative) */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soc#%s_s", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &auxvarsum, name, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvarsum) );

      /* create auxiliary variable for difference (free) */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soc#%s_d", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &auxvardiff, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvardiff) );

      /* add coupling constraint for sum */
      consvars[0] = bilinvar1;
      consvars[1] = bilinvar2;
      consvars[2] = auxvarsum;
      consvals[0] = 1.0;
      consvals[1] = 1.0;
      consvals[2] = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soc#%s_cs", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateConsLinear(scip, &couplingcons, name, 3, consvars, consvals, 0.0, 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), TRUE,
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, couplingcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &couplingcons) );

      /* add coupling constraint for difference */
      consvars[0] = bilinvar1;
      consvars[1] = bilinvar2;
      consvars[2] = auxvardiff;
      consvals[0] = 1.0;
      consvals[1] = -1.0;
      consvals[2] = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "soc#%s_cd", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateConsLinear(scip, &couplingcons, name, 3, consvars, consvals, 0.0, 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), TRUE,
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, couplingcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &couplingcons) );

      /* add difference variable to constraint */
      lhsvars[lhscount] = auxvardiff;
      lhscoefs[lhscount] = 0.5;
      lhsoffsets[lhscount] = 0.0;

      SCIP_CALL( SCIPcreateConsSOC(scip, &upgdconss[0], SCIPconsGetName(cons),
            lhscount + 1, lhsvars, lhscoefs, lhsoffsets, MAX(lhsconstant, 0.0),
            auxvarsum, 0.5, 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      SCIPdebugPrintCons(scip, upgdconss[0], NULL);

      /* create constraint that is equal to cons except that rhs is now infinity */
      if( ! SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) )
      {
         SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[1], SCIPconsGetName(cons),
               SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
               SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
               SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
               SCIPgetLhsQuadratic(scip, cons), SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      }
      SCIP_CALL( SCIPreleaseVar(scip, &auxvardiff) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvarsum) );

      for (i = 0; i < lhscount; ++i)
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, lhsvars[i]) );
      }

      goto cleanup;
   }

   if( rhsvar != NULL && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant) )
   { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax right hand side */
      SCIPdebugMsg(scip, "found right hand side of constraint <%s> to be SOC\n", SCIPconsGetName(cons));

      /* check if upgdconss is long enough to store upgrade constraints:
       * we need two if we will have a quadratic constraint for the left hand side left */
      *nupgdconss = SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) ? 1 : 2;
      if( *nupgdconss > upgdconsssize )
      {
         /* signal that we need more memory and return */
         *nupgdconss = -*nupgdconss;
         goto cleanup;
      }

      SCIP_CALL( SCIPcreateConsSOC(scip, &upgdconss[0], SCIPconsGetName(cons),
            lhscount, lhsvars, lhscoefs, lhsoffsets, MAX(lhsconstant, 0.0),
            rhsvar, rhscoef, rhsoffset,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      SCIPdebugPrintCons(scip, upgdconss[0], NULL);

      /* create constraint that is equal to cons except that rhs is now infinity */
      if( !SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) )
      {
         SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[1], SCIPconsGetName(cons),
               SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
               SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
               SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
               SCIPgetLhsQuadratic(scip, cons), SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      }
   }
   else if( !SCIPisInfinity(scip, - SCIPgetLhsQuadratic(scip, cons)) )
   { /* if the first failed, try if constraint on left hand side is SOC (using negated coefficients) */
      lhscount = 0;
      rhsvar = NULL;

      lhsconstant = SCIPgetLhsQuadratic(scip, cons);

      for( i = 0; i < nquadvars; ++i )
      {
         term = &SCIPgetQuadVarTermsQuadratic(scip, cons)[i];

         /* if there is a linear variable that is still considered as quadratic (constraint probably not presolved yet),
          * then give up
          */
         if( term->sqrcoef == 0.0 )
            goto cleanup;

         if( term->sqrcoef < 0.0 )
         {
            if( lhscount >= nquadvars - 1 )
            { /* too many variables on lhs, i.e., all variables seem to have negative coefficient */
               rhsvar = NULL;
               break;
            }

            lhsvars[lhscount]  = term->var;
            lhscoefs[lhscount] = sqrt(-term->sqrcoef);

            if( term->lincoef != 0.0 )
            {
               lhsoffsets[lhscount] = term->lincoef / (2 * term->sqrcoef);
               lhsconstant += term->lincoef * term->lincoef / (4 * term->sqrcoef);
            }
            else
            {
               lhsoffsets[lhscount] = 0.0;
            }

            ++lhscount;
         }
         else if( rhsvar != NULL ||
               (SCIPisLT(scip, SCIPcomputeVarLbGlobal(scip, term->var), -term->lincoef / (2 * term->sqrcoef))
                && SCIPisGT(scip, SCIPcomputeVarUbGlobal(scip, term->var), -term->lincoef / (2 * term->sqrcoef))) )
         { /* second variable with positive coefficient -> cannot be SOC */
            /* or rhs side changes sign */
            rhsvar = NULL;
            break;
         }
         else
         {
            rhsvar       = term->var;
            rhsoffset    = term->lincoef / (2 * term->sqrcoef);
            lhsconstant += term->lincoef * term->lincoef / (4 * term->sqrcoef);

            if( SCIPisGE(scip, SCIPcomputeVarLbGlobal(scip, term->var), -term->lincoef / (2*term->sqrcoef)) )
               rhscoef = sqrt(term->sqrcoef);
            else
               rhscoef = -sqrt(term->sqrcoef);
         }
      }

      if( rhsvar && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant) )
      { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax left hand side */
         SCIPdebugMsg(scip, "found left hand side of constraint <%s> to be SOC\n", SCIPconsGetName(cons));

         /* check if upgdconss is long enough to store upgrade constraints:
          * we need two if we will have a quadratic constraint for the right hand side left */
         *nupgdconss = SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons)) ? 1 : 2;
         if( *nupgdconss > upgdconsssize )
         {
            /* signal that we need more memory and return */
            *nupgdconss = -*nupgdconss;
            goto cleanup;
         }

         SCIP_CALL( SCIPcreateConsSOC(scip, &upgdconss[0], SCIPconsGetName(cons),
               lhscount, lhsvars, lhscoefs, lhsoffsets, MAX(lhsconstant, 0.0),
               rhsvar, rhscoef, rhsoffset,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
         SCIPdebugPrintCons(scip, upgdconss[0], NULL);

         /* create constraint that is equal to cons except that lhs is now -infinity */
         if( !SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons)) )
         {
            SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[1], SCIPconsGetName(cons),
                  SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
                  SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
                  SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
                  -SCIPinfinity(scip), SCIPgetRhsQuadratic(scip, cons),
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
         }
      }
   }

GENERALUPG:
   if( !trygeneralupg )
      goto cleanup;

   /* find the soc constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdata->generalsocupg )
      goto cleanup;

   SCIPdebugMsg(scip, "Trying general method of upgrade to a soc const\n");

   rhsvarlb = 1.0;
   rhsvarub = 0.0;

#ifndef SCIP_STATISTIC
   /* skip large matrices (needs to compute eigenvalues/vectors) according to presolve timing */
   if( (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 && nquadvars > 10 )
      goto cleanup;
   if( (presoltiming & SCIP_PRESOLTIMING_MEDIUM) != 0 && nquadvars > 50 )
      goto cleanup;
#endif

   /* need routine to compute eigenvalues/eigenvectors */
   if( !SCIPisIpoptAvailableIpopt() )
      goto cleanup;

   /* build lower triangular part the A matrix */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &a, nquadvars*nquadvars) ); /*lint !e647*/
   SCIP_CALL( SCIPallocClearBufferArray(scip, &bp, nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadvars, nquadvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eigvals, nquadvars) );

   /* make sure quadratic variables terms are sorted */
   SCIP_CALL( SCIPsortQuadVarTermsQuadratic(scip, cons) );

   /* set lower triangular entries of A corresponding to square terms */
   for( i = 0; i < nquadvars; ++i )
   {
      term = &SCIPgetQuadVarTermsQuadratic(scip, cons)[i];
      a[i*nquadvars + i] = term->sqrcoef;
      quadvars[i] = term->var;
   }

   /* set lower triangular entries of A corresponding to bilinear terms */
   for( i = 0; i < nbilinterms; ++i )
   {
      int idx1;
      int idx2;

      bilinterm = &SCIPgetBilinTermsQuadratic(scip, cons)[i];

      SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, bilinterm->var1, &idx1) );
      SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, bilinterm->var2, &idx2) );

      assert(idx1 >= 0);
      assert(idx2 >= 0);
      assert(idx1 != idx2);

      a[MIN(idx1,idx2) * nquadvars + MAX(idx1,idx2)] = bilinterm->coef / 2.0;
   }

   /* compute eigenvalues and vectors, A = PDP^t
    * note: a stores P^t
    */
   if( LapackDsyev(TRUE, nquadvars, a, eigvals) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors for constraint <%s>.\n", SCIPconsGetName(cons));
      goto cleanup;
   }

   /* set small eigenvalues to 0 and compute b*P */
   nneg = 0;
   npos = 0;
   for( i = 0; i < nquadvars; ++i )
   {
      for( j = 0; j < nquadvars; ++j )
      {
         term = &SCIPgetQuadVarTermsQuadratic(scip, cons)[j];
         bp[i] += term->lincoef * a[i * nquadvars + j];
      }

      if( SCIPisZero(scip, eigvals[i]) )
      {
         /* if there is a purely linear variable, the constraint can't be written as a SOC */
         if( !SCIPisZero(scip, bp[i]) )
         {
            goto cleanup;
         }

         bp[i] = 0.0;
         eigvals[i] = 0.0;
      }
      else if( eigvals[i] > 0.0 )
         npos++;
      else
         nneg++;
   }

   /* currently, a proper SOC constraint needs at least 3 variables */
   if( npos + nneg < 3 )
      goto cleanup;

   /* determine whether rhs or lhs of cons is potentially SOC, if any */
   rhsissoc = nneg == 1 && !SCIPisInfinity(scip,  SCIPgetRhsQuadratic(scip, cons));
   lhsissoc = npos == 1 && !SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons));

   /* if non is potentially SOC, stop */
   if( !rhsissoc && !lhsissoc )
      goto cleanup;

   assert(rhsissoc != lhsissoc);


   if( lhsissoc )
   {
      /* lhs is potentially SOC, change signs */
      lhsconstant = SCIPgetLhsQuadratic(scip, cons);

      for( i = 0; i < nquadvars; ++i )
      {
         eigvals[i] = -eigvals[i];
         bp[i] = -bp[i];
      }
   }
   else
   {
      lhsconstant = -SCIPgetRhsQuadratic(scip, cons);
   }

   /* we have lhsconstant + x^t A x + b x <= 0 and A has a single negative eigenvalue; try to build soc */
   negeigpos = -1;
   rhsvarfound = FALSE;
   for( i = 0; i < nquadvars; ++i )
   {
      if( SCIPisZero(scip, eigvals[i]) )
         continue;

      if( eigvals[i] > 0.0 )
      {
         lhscoefs[lhscount] = sqrt(eigvals[i]) * UPGSCALE;
         lhsoffsets[lhscount] = bp[i] / (2 * eigvals[i]);
         lhsconstant -= bp[i] * bp[i] / (4 * eigvals[i]);

         ++lhscount;
      }
      else
      {
         /* should enter here only once */
         assert(rhsvarfound == FALSE);

         rhsoffset = bp[i] / (2 * eigvals[i]);

         /* the constraint can only be a soc if the resulting rhs var does not change var; the rhs var is going to be a
          * multiaggregated variable, so estimate its bounds
          */
         rhsvarlb = 0.0;
         rhsvarub = 0.0;
         for( j = 0; j < nquadvars; ++j )
         {
            SCIP_Real aux;

            if( a[i * nquadvars + j] > 0.0 )
            {
               aux = SCIPcomputeVarLbGlobal(scip, quadvars[j]);
               assert(!SCIPisInfinity(scip, aux));
            }
            else
            {
               aux = SCIPcomputeVarUbGlobal(scip, quadvars[j]);
               assert(!SCIPisInfinity(scip, -aux));
            }

            if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
            {
               rhsvarlb = -SCIPinfinity(scip);
               break;
            }
            else
               rhsvarlb += a[i * nquadvars + j] * aux;
         }
         rhsvarlb += rhsoffset;

         for( j = 0; j < nquadvars; ++j )
         {
            SCIP_Real aux;

            if( a[i * nquadvars + j] > 0.0 )
            {
               aux = SCIPcomputeVarUbGlobal(scip, quadvars[j]);
               assert(!SCIPisInfinity(scip, -aux));
            }
            else
            {
               aux = SCIPcomputeVarLbGlobal(scip, quadvars[j]);
               assert(!SCIPisInfinity(scip, aux));
            }

            if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
            {
               rhsvarub = SCIPinfinity(scip);
               break;
            }
            else
               rhsvarub += a[i * nquadvars + j] * aux;
         }
         rhsvarub += rhsoffset;

         /* since we are just interested in obtaining an interval that contains the real bounds and is tight enough so
          * that we can identify that the rhsvar does not change sign, we swap the bounds in case of numerical troubles
          */
         if( rhsvarub < rhsvarlb )
         {
            assert(SCIPisEQ(scip, rhsvarub, rhsvarlb));
            SCIPswapReals(&rhsvarub, &rhsvarlb);
         }

         /* check whether rhsvar changes sign */
         if( SCIPisGE(scip, rhsvarlb, 0.0) || SCIPisLE(scip, rhsvarub, 0.0) )
         {
            rhsvarfound  = TRUE;
            negeigpos    = i;
            lhsconstant -= rhsoffset * rhsoffset * eigvals[i];
            rhscoef      = SCIPisLE(scip, rhsvarub, 0.0) ? -sqrt(-eigvals[i]) : sqrt(-eigvals[i]);
         }
         else
         {
            SCIPdebugMsg(scip, "Failed because rhsvar [%g, %g] changes sign.\n", rhsvarlb, rhsvarub);
            rhsvarfound = FALSE;
            break;
         }
      }
   }

   if( rhsvarfound && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant) )
   {
      /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax right hand side */
      SCIPdebugMsg(scip, "found right hand side of constraint <%s> to be SOC\n", SCIPconsGetName(cons));

      /* check if upgdconss is long enough to store upgrade constraints:
       * we need two if we will have a quadratic constraint for the left hand side left */
      if( rhsissoc )
         *nupgdconss = SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) ? 1 : 2;
      else
         *nupgdconss = SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons)) ? 1 : 2;
      if( *nupgdconss > upgdconsssize )
      {
         /* signal that we need more memory and return */
         *nupgdconss = -*nupgdconss;
         goto cleanup;
      }

      /* create lhs and rhs vars */
      lhsnvars = 0;
      for( i = 0; i < nquadvars; ++i )
      {
         if( eigvals[i] <= 0.0 )
            continue;

         SCIP_CALL( SCIPcreateVarBasic(scip, &lhsvars[lhsnvars], NULL,
                  -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, lhsvars[lhsnvars]) );

         SCIP_CALL( SCIPmultiaggregateVar(scip, lhsvars[lhsnvars], nquadvars, quadvars, &(a[i * nquadvars]),
                  lhsoffsets[lhsnvars], &infeas, &success) );

         if( infeas || !success )
         {
            SCIPdebugMsg(scip, "Problem with aggregation while trying to upgrade <%s>.\n", SCIPconsGetName(cons) );

            /* release all created vars so far */
            for( j = 0; j <= lhsnvars; ++j )
               SCIP_CALL( SCIPreleaseVar(scip, &lhsvars[j]) );

            goto cleanup;
         }
         lhsnvars++;
      }
      assert(lhsnvars == lhscount);
      assert(negeigpos >= 0);

      SCIP_CALL( SCIPcreateVarBasic(scip, &rhsvar, NULL,
               rhsvarlb, rhsvarub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, rhsvar) );
      SCIP_CALL( SCIPmultiaggregateVar(scip, rhsvar, nquadvars, quadvars, &(a[negeigpos * nquadvars]),
               rhsoffset, &infeas, &success) );

      if( infeas || !success )
      {
         SCIPdebugMsg(scip, "Problem with aggregation while trying to upgrade <%s>.\n", SCIPconsGetName(cons) );

         /* release all created vars */
         SCIP_CALL( SCIPreleaseVar(scip, &rhsvar) );
         for( j = 0; j < lhsnvars; ++j )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &lhsvars[j]) );
         }

         goto cleanup;
      }

      SCIP_CALL( SCIPcreateConsSOC(scip, &upgdconss[0], SCIPconsGetName(cons),
               lhscount, lhsvars, lhscoefs, NULL, MAX(lhsconstant, 0.0) * UPGSCALE * UPGSCALE,
               rhsvar, rhscoef * UPGSCALE, 0.0,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      SCIPdebugPrintCons(scip, upgdconss[0], NULL);

      /* release created variables */
      SCIP_CALL( SCIPreleaseVar(scip, &rhsvar) );
      for( i = 0; i < lhsnvars; ++i )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &lhsvars[i]) );
      }

      /* create constraint that is equal to cons except that rhs/lhs is now +/-infinity */
      if( rhsissoc && !SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) )
      {
         SCIPdebugMsg(scip, "rhs is soc, keep quadratic\n");
         SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[1], SCIPconsGetName(cons),
                  SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
                  SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
                  SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
                  SCIPgetLhsQuadratic(scip, cons), SCIPinfinity(scip),
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      }
      else if( lhsissoc && !SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip,cons)) )
      {
         SCIPdebugMsg(scip, "lhs is soc, keep quadratic\n");
         SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[1], SCIPconsGetName(cons),
                  SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
                  SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
                  SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
                  -SCIPinfinity(scip), SCIPgetRhsQuadratic(scip, cons),
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      }
   }
#ifdef SCIP_DEBUG
   else
   {
      if( lhscount < 2 )
      {
         SCIPdebugMsg(scip, "Failed because there are not enough lhsvars (%d)\n", lhscount);
      }
      if( SCIPisNegative(scip, lhsconstant) )
      {
         SCIPdebugMsg(scip, "Failed because lhsconstant is negative (%g)\n", lhsconstant);
      }
   }
#endif

 cleanup:
   SCIPfreeBufferArray(scip, &lhsoffsets);
   SCIPfreeBufferArray(scip, &lhscoefs);
   SCIPfreeBufferArray(scip, &lhsvars);
   SCIPfreeBufferArrayNull(scip, &a);
   SCIPfreeBufferArrayNull(scip, &bp);
   SCIPfreeBufferArrayNull(scip, &quadvars);
   SCIPfreeBufferArrayNull(scip, &eigvals);

   return SCIP_OKAY;
} /*lint !e715*/
#endif

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySOC)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSOC(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSOC)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitSOC)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur  = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur  = SCIPfindHeur(scip, "trysol");
   conshdlrdata->haveexprint = (strcmp(SCIPexprintGetName(), "NONE") != 0);

   /* mark constraints for propagation */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( SCIPmarkConsPropagate(scip, conss[c]) );  /*lint !e613*/
   }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitSOC)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur  = NULL;
   conshdlrdata->trysolheur  = NULL;
   conshdlrdata->haveexprint = FALSE;

   return SCIP_OKAY;
}

/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreSOC)
{  /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   /* tell SCIP that we have something nonlinear */
   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsAdded(conss[c]) ) /*lint !e613*/
      {
         SCIPenableNLP(scip);
         break;
      }
   }

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolSOC)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   /* add nlrow representation to NLP, if NLP has been enabled */
   if( SCIPisNLPConstructed(scip) )
   {
      for( c = 0; c < nconss; ++c )
      {
         if( SCIPconsIsEnabled(conss[c]) )
         {
            consdata = SCIPconsGetData(conss[c]);
            assert(consdata != NULL);

            if( consdata->nlrow == NULL )
            {
               SCIP_CALL( createNlRow(scip, conshdlr, conss[c]) );
               assert(consdata->nlrow != NULL);
            }
            SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
         }
      }
   }

   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   /* reset flags and counters */
   conshdlrdata->sepanlp = FALSE;
   conshdlrdata->lastenfonode = NULL;
   conshdlrdata->nenforounds = 0;

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolSOC)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSOC)
{
   int i;

   assert(scip      != NULL);
   assert(conshdlr  != NULL);
   assert(cons      != NULL);
   assert(consdata  != NULL);
   assert(*consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   SCIPdebugMsg(scip, "Deleting SOC constraint <%s>.\n", SCIPconsGetName(cons) );

   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   for( i = 0; i < (*consdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars,    (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->coefs,   (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->offsets, (*consdata)->nvars);
   assert((*consdata)->lhsbndchgeventdata == NULL);

   if( (*consdata)->rhsvar != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->rhsvar) );
   }

   /* free nonlinear row representation
    * normally released in exitsol, but constraint may be deleted early (e.g., if found redundant)
    */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransSOC)
{  
   SCIP_CONSDATA*     consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     sourcedata;
   char               s[SCIP_MAXSTRLEN];
   int                i;

   assert(scip       != NULL);
   assert(conshdlr   != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMsg(scip, "Transforming SOC constraint: <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata       != NULL);
   assert(sourcedata->vars != NULL);
   assert(sourcedata->coefs != NULL);
   assert(sourcedata->offsets != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = sourcedata->nvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, consdata->nvars) );
   SCIP_CALL( SCIPgetTransformedVars(scip, consdata->nvars, sourcedata->vars, consdata->vars) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
   }

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->coefs,   sourcedata->coefs,   consdata->nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->offsets, sourcedata->offsets, consdata->nvars) );
   consdata->constant = sourcedata->constant;

   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->rhsvar, &consdata->rhsvar) );
   consdata->rhsoffset = sourcedata->rhsoffset;
   consdata->rhscoeff  = sourcedata->rhscoeff;

   SCIP_CALL( SCIPcaptureVar(scip, consdata->rhsvar) );

   consdata->nlrow = NULL;
   consdata->lhsbndchgeventdata = NULL;
   consdata->isapproxadded = FALSE;

   /* create transformed constraint with the same flags */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *targetcons) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSOC)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          sepasuccess;
   SCIP_Bool cutoff;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, NULL, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* at root, check if we want to solve the NLP relaxation and use its solutions as reference point
    * if there is something convex, then linearizing in the solution of the NLP relaxation can be very useful
    */
   if( SCIPgetDepth(scip) == 0 && !conshdlrdata->sepanlp &&
      (SCIPgetNContVars(scip) >= conshdlrdata->sepanlpmincont * SCIPgetNVars(scip) || (SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY && conshdlrdata->sepanlpmincont <= 1.0)) &&
      SCIPisNLPConstructed(scip) && SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_NLPSOLSTAT solstat;
      SCIP_Bool       solvednlp;

      solstat = SCIPgetNLPSolstat(scip);
      solvednlp = FALSE;
      if( solstat == SCIP_NLPSOLSTAT_UNKNOWN )
      {
         /* NLP is not solved yet, so we do it now and update solstat */

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
            int nfracvars = 0;

            if( SCIPgetNBinVars(scip) > 0 || SCIPgetNIntVars(scip) > 0 )
            {
               SCIP_CALL( SCIPgetNLPFracVars(scip, NULL, NULL, NULL, &nfracvars, NULL) );
            }

            if( nfracvars == 0 )
            {
               SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, nlpsol) );
            }
         }

         SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, nlpsol, &lpsolseparated, SCIPgetSepaMinEfficacy(scip), &cutoff) );

         SCIP_CALL( SCIPfreeSol(scip, &nlpsol) );

         if( cutoff )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

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

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, FALSE, &cutoff, &sepasuccess) );
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if ( sepasuccess )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolSOC)
{  
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          sepasuccess;
   SCIP_Bool          cutoff;

   assert(scip     != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);
   assert(sol      != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, sol, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, FALSE, &cutoff, &sepasuccess) );
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if ( sepasuccess )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSOC)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, NULL, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSOC)
{  /*lint --e{715}*/
   SCIP_CALL( enforceConstraint(scip, conshdlr, conss, nconss, nusefulconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSOC)
{
   SCIP_CONS*         maxviolcons;

   assert(scip     != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   SCIP_CALL( computeViolations(scip, conshdlr, conss, nconss, NULL, &maxviolcons) );

   if( maxviolcons == NULL )
      *result = SCIP_FEASIBLE;

   *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
} /*lint !e715*/


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSOC)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   SCIP_Bool          dolinfeasshift;
   SCIP_SOL*          polishedsol;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result     = SCIP_FEASIBLE;
   maxviol     = 0.0;

   dolinfeasshift = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL);
   polishedsol = NULL;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conshdlr, conss[c], sol) );  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* if feasible, just continue */
      if( !SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) )
         continue;

      *result = SCIP_INFEASIBLE;

      if( consdata->violation > maxviol )
         maxviol = consdata->violation;

      if( printreason )
      {
         SCIP_Real unscaledviol;

         unscaledviol  = consdata->lhsval;
         if( !SCIPisInfinity(scip, unscaledviol) )
            unscaledviol -= consdata->rhscoeff * (SCIPgetSolVal(scip, sol, consdata->rhsvar) + consdata->rhsoffset);

         SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );  /*lint !e613*/            
         SCIPinfoMessage(scip, NULL, ";\n\tviolation: %g (scaled: %g)\n", unscaledviol, consdata->violation);
      }

      /* if we do linear feasibility shifting, then try to adjust solution */
      if( dolinfeasshift )
      {
         if( SCIPvarGetStatus(consdata->rhsvar) != SCIP_VARSTATUS_MULTAGGR &&
            !SCIPisInfinity(scip, REALABS(consdata->lhsval)) &&
            (  (consdata->rhscoeff > 0.0 && SCIPvarMayRoundUp  (consdata->rhsvar)) ||
               (consdata->rhscoeff < 0.0 && SCIPvarMayRoundDown(consdata->rhsvar)) ) )
         {
            SCIP_Bool success;

            if( polishedsol == NULL )
            {
               if( sol != NULL )
               {
                  SCIP_CALL( SCIPcreateSolCopy(scip, &polishedsol, sol) );
               }
               else
               {
                  SCIP_CALL( SCIPcreateLPSol(scip, &polishedsol, NULL) );
               }
               SCIP_CALL( SCIPunlinkSol(scip, polishedsol) );
            }
            SCIP_CALL( polishSolution(scip, conss[c], polishedsol, &success) );  /*lint !e613*/

            /* disable solution polishing if we failed for this constraint */
            dolinfeasshift = success;
         }
         else /* if locks of the variable are bad or rhs is multi-aggregated, disable solution polishing */
         {
            dolinfeasshift = FALSE;
         }
      }

      /* if solution polishing is off and there is no NLP heuristic or we just check the LP solution,
       * then there is no need to check remaining constraints (NLP heuristic will pick up LP solution anyway) */
      if( !dolinfeasshift && (conshdlrdata->subnlpheur == NULL || sol == NULL) && !completely )
         break;
   }

   /* if we failed to polish solution, clear the solution */ 
   if( !dolinfeasshift && polishedsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &polishedsol) );
   }

   if( polishedsol != NULL )
   {
      assert(*result == SCIP_INFEASIBLE);
      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, polishedsol) );
      SCIP_CALL( SCIPfreeSol(scip, &polishedsol) );
   }
   else if( conshdlrdata->subnlpheur != NULL && sol != NULL && *result == SCIP_INFEASIBLE && !SCIPisInfinity(scip, maxviol) )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSOC)
{
   SCIP_RESULT propresult;
   int         c;
   int         nchgbds;
   SCIP_Bool   redundant;

   assert(scip     != NULL);
   assert(conss    != NULL || ((nconss == 0) && (nmarkedconss == 0)));
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;
   nchgbds = 0;

   for( c = 0; c < nmarkedconss && *result != SCIP_CUTOFF; ++c )
   {
      SCIP_CALL( propagateBounds(scip, conss[c], &propresult, &nchgbds, &redundant) );  /*lint !e613*/
      if( propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN )
         *result = propresult;
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSOC)
{
   SCIP_CONSHDLRDATA*  conshdlrdata;
   SCIP_CONSDATA*      consdata;
   int                 c;
   SCIP_RESULT         propresult;
   SCIP_Bool           iscutoff;
   SCIP_Bool           isdeleted;

   assert(scip     != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(conshdlr != NULL);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      SCIP_CALL( presolveRemoveFixedVariables(scip, conshdlr, conss[c], ndelconss, nupgdconss, nchgbds, nfixedvars, &iscutoff, &isdeleted) );  /*lint !e613*/
      if( iscutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( isdeleted )
      {
         /* conss[c] has been deleted */
         *result = SCIP_SUCCESS;
         continue;
      }

      if( conshdlrdata->nauxvars > 0 && !consdata->isapproxadded && consdata->nvars > 1 )
      {
         SCIP_CALL( presolveCreateOuterApprox(scip, consdata->nvars, consdata->vars, consdata->coefs, consdata->offsets,
               consdata->rhsvar, consdata->rhscoeff, consdata->rhscoeff, consdata->constant, SCIPconsGetName(conss[c]), conss[c],
               conshdlrdata->nauxvars, conshdlrdata->glineur, naddconss) );  /*lint !e613*/
         consdata->isapproxadded = TRUE;
      }

      if( (presoltiming & SCIP_PRESOLTIMING_FAST) != 0 )
      {
         SCIP_Bool redundant;

         SCIP_CALL( propagateBounds(scip, conss[c], &propresult, nchgbds, &redundant) );  /*lint !e613*/
         switch( propresult )
         {
            case SCIP_DIDNOTRUN:
            case SCIP_DIDNOTFIND:
               break;
            case SCIP_REDUCEDDOM:
               *result = SCIP_SUCCESS;
               break;
            case SCIP_CUTOFF:
               *result = SCIP_CUTOFF;
               SCIPdebugMsg(scip, "infeasible in presolve due to propagation for constraint %s\n", SCIPconsGetName(conss[c]));  /*lint !e613*/
               return SCIP_OKAY;
            default:
               SCIPerrorMessage("unexpected result from propagation: %d\n", propresult);
               return SCIP_ERROR;
         } /*lint !e788*/
         if( redundant )
            ++*ndelconss;
      }

      /* disaggregate each lhs term to a quadratic constraint by using auxiliary variables */
      if( conshdlrdata->disaggregate && (presoltiming & SCIP_PRESOLTIMING_EXHAUSTIVE) != 0 )
      {
         SCIP_Bool success;

         SCIP_CALL( disaggregate(scip, conss[c], consdata, naddconss, ndelconss, &success) ); /*lint !e613*/

         if( success )
         {
            SCIPdebugMsg(scip, "disaggregated SOC constraint\n");

            /* conss[c] has been deleted */
            *result = SCIP_SUCCESS;
            continue;
         }
      }
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSOC)
{
   SCIP_CONSDATA* consdata;
   int            i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(cons     != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMsg(scip, "Locking constraint <%s>.\n", SCIPconsGetName(cons));

   /* Changing variables x_i, i <= n, in both directions can lead to an infeasible solution. */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   /* Rounding x_{n+1} up will not violate a solution. */
   if( consdata->rhsvar != NULL )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->rhsvar, consdata->rhscoeff > 0.0 ? nlockspos : nlocksneg, consdata->rhscoeff > 0.0 ? nlocksneg : nlockspos) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSOC)
{  
   SCIP_CONSDATA* consdata;
   int            i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(cons     != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPinfoMessage(scip, file, "sqrt( ");
   if( consdata->constant != 0.0 )
   {
      SCIPinfoMessage(scip, file, "%.15g", consdata->constant);
   }

   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIPinfoMessage(scip, file, "+ (%.15g*(", consdata->coefs[i]);
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[i], TRUE) );
      SCIPinfoMessage(scip, file, "%+.15g))^2 ", consdata->offsets[i]);
   }

   SCIPinfoMessage(scip, file, ") <= ");
   if( consdata->rhsvar != NULL )
   {
      SCIPinfoMessage(scip, file, "%.15g*(", consdata->rhscoeff);
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->rhsvar, TRUE) );
      SCIPinfoMessage(scip, file, "%+.15g)", consdata->rhsoffset);
   }
   else
   {
      SCIPinfoMessage(scip, file, "%.15g", consdata->rhscoeff*consdata->rhsoffset);
   }

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySOC)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR**     vars;
   SCIP_VAR*      rhsvar;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);
   assert(stickingatnode == FALSE);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);

   *valid = TRUE; 

   SCIP_CALL( SCIPallocBufferArray(sourcescip, &vars, consdata->nvars) );

   /* map variables to active variables of the target SCIP */   
   for( i = 0; i < consdata->nvars && *valid; ++i )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->vars[i], &vars[i], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[i] != NULL);
   }

   /* map rhs variable to active variable of the target SCIP */   
   if( *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->rhsvar, &rhsvar, varmap, consmap, global, valid) );
      assert(!(*valid) || rhsvar != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsSOC(scip, cons, name ? name : SCIPconsGetName(sourcecons),
            consdata->nvars, vars, consdata->coefs, consdata->offsets, consdata->constant,
            rhsvar, consdata->rhscoeff, consdata->rhsoffset,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );  /*lint !e644 */
   }

   SCIPfreeBufferArray(sourcescip, &vars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseSOC)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real* offsets;
   int nvars;
   int varssize;
   SCIP_VAR* rhsvar;
   SCIP_Real rhscoef;
   SCIP_Real rhsoffset;
   SCIP_Real constant;
   SCIP_Real coef;
   SCIP_Real offset;
   char* endptr;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* check that string starts with "sqrt( " */
   if( strncmp(str, "sqrt( ", 6) != 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected 'sqrt( ' at begin of soc constraint string '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str += 6;

   *success = TRUE;

   /* check if we have a constant in the beginning */
   if( SCIPstrToRealValue(str, &constant, &endptr) )
      str = endptr;
   else
      constant = 0.0;

   nvars = 0;
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,    varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs,   varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, varssize) );

   /* read '+ (coef*(var+offset))^2' on lhs, as long as possible */
   while( *str != '\0' )
   {
      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;

      /* stop if no more coefficients */
      if( strncmp(str, "+ (", 3) != 0 )
         break;

      str += 3;

      /* parse coef */
      if( !SCIPstrToRealValue(str, &coef, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected coefficient at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str = endptr;

      if( strncmp(str, "*(", 2) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '*(' at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str += 2;

      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, str, &var, &endptr) );
      if( var == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
         *success = FALSE;
         break;
      }
      str = endptr;

      /* parse offset */
      if( !SCIPstrToRealValue(str, &offset, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected offset at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str = endptr;

      if( strncmp(str, "))^2", 4) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '))^2' at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str += 4;

      if( varssize <= nvars )
      {
         varssize = SCIPcalcMemGrowSize(scip, varssize+1);
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars,    varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs,   varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &offsets, varssize) );
      }
      vars[nvars]    = var;
      coefs[nvars]   = coef;
      offsets[nvars] = offset;
      ++nvars;
   }

   if( strncmp(str, ") <=", 4) != 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected ') <=' at begin of '%s'\n", str);
      *success = FALSE;
   }
   str += 4;

   /* read rhs coef*(var+offset) or just a constant */

   /* parse coef */
   if( *success )
   {
      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;

      if( !SCIPstrToRealValue(str, &rhscoef, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected coefficient at begin of '%s'\n", str);
         *success = FALSE;
      }
      str = endptr;

      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;
   }

   /* parse *(var+offset) */
   if( *str != '\0' )
   {
      if( *success )
      {
         if( strncmp(str, "*(", 2) != 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '*(' at begin of '%s'\n", str);
            *success = FALSE;
         }
         else
         {
            str += 2;
         }
      }

      /* parse variable name */
      if( *success )
      {
         SCIP_CALL( SCIPparseVarName(scip, str, &rhsvar, &endptr) );
         if( rhsvar == NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
            *success = FALSE;
         }
         else
         {
            str = endptr;
         }
      }

      /* parse offset */
      if( *success )
      {
         if( !SCIPstrToRealValue(str, &rhsoffset, &endptr) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected offset at begin of '%s'\n", str);
            *success = FALSE;
         }
         else
         {
            str = endptr;
         }
      }

      if( *success )
      {
         if( *str != ')' )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected ')' at begin of '%s'\n", str);
            *success = FALSE;
         }
      }
   }
   else if( *success )
   {
      /* only a constant at right hand side */
      rhsoffset = rhscoef;  /*lint !e644*/
      rhscoef = 1.0;
      rhsvar = NULL;
   }

   if( *success )
   {
      assert(!stickingatnode);
      SCIP_CALL( SCIPcreateConsSOC(scip, cons, name, nvars, vars, coefs, offsets, constant, rhsvar, rhscoef, rhsoffset,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );  /*lint !e644 */
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &offsets);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSOC)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( varssize < consdata->nvars + 1 )
      (*success) = FALSE;
   else
   {
      BMScopyMemoryArray(vars, consdata->vars, consdata->nvars);
      vars[consdata->nvars] = consdata->rhsvar;
      (*success) = TRUE;
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variable (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSOC)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   (*nvars) = consdata->nvars + 1;
   (*success) = TRUE;

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for second order cone constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSOC(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_EVENTHDLR* eventhdlr;

   /* create constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, CONSHDLR_NAME"_boundchange",
         "signals a bound change to a second order cone constraint",
         processVarEvent, NULL) );
   conshdlrdata->eventhdlr = eventhdlr;

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, CONSHDLR_NAME"_newsolution",
         "handles the event that a new primal solution has been found",
         processNewSolutionEvent, NULL) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSOC, consEnfopsSOC, consCheckSOC, consLockSOC,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySOC, consCopySOC) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSOC) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitSOC) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSOC) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolSOC) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSOC) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSOC) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSOC) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitSOC) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolSOC) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSOC) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSOC, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSOC) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSOC, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSOC, consSepasolSOC, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSOC) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSOC) );

   if( SCIPfindConshdlr(scip,"quadratic") != NULL )
   {
      /* notify function that upgrades quadratic constraint to SOC's */
      SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, upgradeConsQuadratic, QUADCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   /* add soc constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/projectpoint",
         "whether the reference point of a cut should be projected onto the feasible set of the SOC constraint",
         &conshdlrdata->projectpoint,     TRUE,  FALSE,         NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "constraints/" CONSHDLR_NAME "/nauxvars",
         "number of auxiliary variables to use when creating a linear outer approx. of a SOC3 constraint; 0 to turn off",
         &conshdlrdata->nauxvars,         FALSE, 0, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/glineur",
         "whether the Glineur Outer Approximation should be used instead of Ben-Tal Nemirovski",
         &conshdlrdata->glineur,          FALSE, TRUE,          NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/sparsify",
         "whether to sparsify cuts",
         &conshdlrdata->sparsify,         TRUE,  FALSE,         NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/sparsifymaxloss",
         "maximal loss in cut efficacy by sparsification",
         &conshdlrdata->sparsifymaxloss,  TRUE,  0.2, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/sparsifynzgrowth",
         "growth rate of maximal allowed nonzeros in cuts in sparsification",
         &conshdlrdata->sparsifynzgrowth, TRUE,  1.3, 1.000001, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/linfeasshift",
         "whether to try to make solutions feasible in check by shifting the variable on the right hand side",
         &conshdlrdata->linfeasshift,     FALSE, TRUE,          NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/" CONSHDLR_NAME "/nlpform",
         "which formulation to use when adding a SOC constraint to the NLP (a: automatic, q: nonconvex quadratic form, s: convex sqrt form, e: convex exponential-sqrt form, d: convex division form)",
         &conshdlrdata->nlpform,          FALSE, 'a', "aqsed", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/sepanlpmincont",
         "minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation",
         &conshdlrdata->sepanlpmincont, FALSE, 1.0, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/enfocutsremovable",
         "are cuts added during enforcement removable from the LP in the same node?",
         &conshdlrdata->enfocutsremovable, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/generalsocupgrade",
         "try to upgrade more general quadratics to soc?",
         &conshdlrdata->generalsocupg, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/" CONSHDLR_NAME "/disaggregate",
         "try to completely disaggregate soc?",
         &conshdlrdata->disaggregate, TRUE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a second order cone constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSOC(
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
   SCIP_Real             rhsoffset,          /**< offset of variable on right hand side (beta_{n+1}) */
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
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(modifiable == FALSE); /* we do not support column generation */

   /* find the soc constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("SOC constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert(vars     != NULL);
   assert(nvars    >= 1);
   assert(constant >= 0.0);
   assert(!SCIPisInfinity(scip, ABS(rhsoffset)));
   assert(!SCIPisInfinity(scip, constant));
   assert(rhsvar == NULL || rhscoeff <= 0.0 || SCIPisGE(scip, local ? SCIPcomputeVarLbLocal(scip, rhsvar) : SCIPcomputeVarLbGlobal(scip, rhsvar), -rhsoffset));
   assert(rhsvar == NULL || rhscoeff >= 0.0 || SCIPisLE(scip, local ? SCIPcomputeVarUbLocal(scip, rhsvar) : SCIPcomputeVarUbGlobal(scip, rhsvar), -rhsoffset));

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = nvars;
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->vars, vars, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPcaptureVar(scip, vars[i]) );
   }

   if( coefs != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->coefs, coefs, nvars) );
      for( i = 0; i < nvars; ++i )
      {
         if( consdata->coefs[i] < 0.0 )
            consdata->coefs[i] = -consdata->coefs[i];
      }
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->coefs, nvars) );
      for( i = 0; i < nvars; ++i )
         consdata->coefs[i] = 1.0;
   }

   if( offsets != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->offsets, offsets, nvars) );
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->offsets, nvars) );
      BMSclearMemoryArray(consdata->offsets, nvars);
   }

   consdata->constant  = constant;
   consdata->rhsvar    = rhsvar;
   consdata->rhscoeff  = rhscoeff;
   consdata->rhsoffset = rhsoffset;

   if( rhsvar != NULL )
   {
      SCIP_CALL( SCIPcaptureVar(scip, rhsvar) );
   }

   consdata->nlrow = NULL;

   consdata->lhsbndchgeventdata = NULL;
   consdata->isapproxadded       = FALSE;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   if( SCIPisTransformed(scip) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );
   }

   return SCIP_OKAY;
}

/** creates and captures a second order cone constraint with all its constraint flags
 *  set to their default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicSOC(
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
   SCIP_CALL( SCIPcreateConsSOC(scip, cons, name, nvars, vars, coefs, offsets, constant,
         rhsvar, rhscoeff, rhsoffset,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** Gets the SOC constraint as a nonlinear row representation. */
SCIP_RETCODE SCIPgetNlRowSOC(
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
      SCIP_CALL( createNlRow(scip, SCIPconsGetHdlr(cons), cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** Gets the number of variables on the left hand side of a SOC constraint. */
int SCIPgetNLhsVarsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nvars;
}

/** Gets the variables on the left hand side of a SOC constraint. */
SCIP_VAR** SCIPgetLhsVarsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->vars;
}

/** Gets the coefficients of the variables on the left hand side of a SOC constraint, or NULL if all are equal to 1.0. */
SCIP_Real* SCIPgetLhsCoefsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->coefs;
}

/** Gets the offsets of the variables on the left hand side of a SOC constraint, or NULL if all are equal to 0.0. */
SCIP_Real* SCIPgetLhsOffsetsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->offsets;
}

/** Gets the constant on the left hand side of a SOC constraint. */
SCIP_Real SCIPgetLhsConstantSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->constant;
}

/** Gets the variable on the right hand side of a SOC constraint. */
SCIP_VAR* SCIPgetRhsVarSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhsvar;
}

/** Gets the coefficient of the variable on the right hand side of a SOC constraint. */
SCIP_Real SCIPgetRhsCoefSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhscoeff;
}

/** Gets the offset of the variables on the right hand side of a SOC constraint. */
SCIP_Real SCIPgetRhsOffsetSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhsoffset;
}

/** Adds the constraint to an NLPI problem.
 *  Uses nonconvex formulation as quadratic function.
 */
SCIP_RETCODE SCIPaddToNlpiProblemSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SOC constraint */
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
   int            j;
   int            lincnt;
   SCIP_Real      lhs;
   SCIP_Real      rhs;
   const char*    name;

   assert(scip     != NULL);
   assert(cons     != NULL);
   assert(nlpi     != NULL);
   assert(nlpiprob != NULL);
   assert(scipvar2nlpivar != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhs = -SCIPinfinity(scip);
   rhs = -consdata->constant;

   /* count how length is the linear part, i.e., how many offsets we have */
   nlininds = consdata->rhsoffset != 0.0 ? 1 : 0;
   for( j = 0; j < consdata->nvars; ++j )
      if( consdata->offsets[j] != 0.0 )
         ++nlininds;

   lininds = NULL;
   linvals = NULL;
   if( nlininds )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nlininds) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nlininds) );
   }
   lincnt = 0;

   nquadelems = consdata->nvars + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nquadelems) );

   for( j = 0; j < consdata->nvars; ++j )
   {
      quadelems[j].idx1 = (int) (size_t) SCIPhashmapGetImage(scipvar2nlpivar, consdata->vars[j]);
      quadelems[j].idx2 = quadelems[j].idx1;
      quadelems[j].coef = consdata->coefs[j] * consdata->coefs[j];

      if( consdata->offsets[j] != 0.0 )
      {
         assert(lininds != NULL);
         assert(linvals != NULL);
         lininds[lincnt] = quadelems[j].idx1;
         linvals[lincnt] = 2 * quadelems[j].coef * consdata->offsets[j];
         ++lincnt;

         rhs -= quadelems[j].coef * consdata->offsets[j] * consdata->offsets[j];
      }
   }
   quadelems[consdata->nvars].idx1 = (int) (size_t) SCIPhashmapGetImage(scipvar2nlpivar, consdata->rhsvar);
   quadelems[consdata->nvars].idx2 = quadelems[consdata->nvars].idx1;
   quadelems[consdata->nvars].coef = - consdata->rhscoeff * consdata->rhscoeff;

   if( consdata->rhsoffset != 0.0 )
   {
      assert(lininds != NULL);
      assert(linvals != NULL);
      lininds[lincnt] = quadelems[consdata->nvars].idx1;
      linvals[lincnt] = -2.0 * consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset;
      ++lincnt;

      rhs += consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset * consdata->rhsoffset;
   }

   assert(lincnt == nlininds);

   name = names ? SCIPconsGetName(cons) : NULL;

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1,
         &lhs, &rhs,
         &nlininds, &lininds, &linvals,
         &nquadelems, &quadelems,
         NULL, NULL, &name) );

   SCIPfreeBufferArrayNull(scip, &lininds);
   SCIPfreeBufferArrayNull(scip, &linvals);
   SCIPfreeBufferArray(scip, &quadelems);

   return SCIP_OKAY;
}
