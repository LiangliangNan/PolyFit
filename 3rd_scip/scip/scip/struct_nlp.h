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

/**@file   struct_nlp.h
 * @ingroup INTERNALAPI
 * @brief  datastructures for NLP management
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 *
 *  In SCIP, the NLP is defined as follows:
 *
 *   min         const + obj * x + <x, Qx> + f(x)
 *        lhs <= const + A   * x                  <= rhs
 *        lhs <= const + A   * x + <x, Qx> + f(x) <= rhs
 *        lb  <=               x                  <= ub
 *
 *  where the linear rows and variable bounds are managed by the LP
 *  and the nonlinear rows are managed by the NLP.
 *
 *  The row activities are defined as
 *     activity = A * x + const
 *  for a linear row and as
 *     activity = f(x) + <x, Qx> + A * x + const
 *  for a nonlinear row.
 *  The activities must therefore be in the range of [lhs,rhs].
 *
 *  The main datastructures for storing an NLP are the nonlinear rows.
 *  A nonlinear row can live on its own (if it was created by a separator),
 *  or as relaxation of a constraint. Thus, it has a nuses-counter, and is
 *  deleted, if not needed any more.
 *  In difference to columns of an LP, nonlinear rows are defined
 *  with respect SCIP variables.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_NLP_H__
#define __SCIP_STRUCT_NLP_H__

#include "scip/def.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "nlpi/type_nlpi.h"
#include "nlpi/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** NLP row */
struct SCIP_NlRow
{
   /* sides */
   SCIP_Real             lhs;                /**< left hand side */
   SCIP_Real             rhs;                /**< right hand side */

   /* constant part */
   SCIP_Real             constant;           /**< constant value */

   /* linear part */
   int                   nlinvars;           /**< number of linear variables */
   int                   linvarssize;        /**< size of arrays storing linear part of row */
   SCIP_VAR**            linvars;            /**< linear variables */
   double*               lincoefs;           /**< coefficients of linear variables */
   SCIP_Bool             linvarssorted;      /**< are the linear coefficients sorted (by variable indices?) */

   /* quadratic part */
   int                   nquadvars;          /**< number of variables in quadratic terms */
   int                   quadvarssize;       /**< size of array storing quadratic variables of row */
   SCIP_VAR**            quadvars;           /**< variables in quadratic term */
   SCIP_HASHMAP*         quadvarshash;       /**< hash map from variable to indices in quadvars */
   int                   nquadelems;         /**< number of entries in quadratic matrix */
   int                   quadelemssize;      /**< size of quadratic elements array */
   SCIP_QUADELEM*        quadelems;          /**< entries in quadratic matrix */
   SCIP_Bool             quadelemssorted;    /**< are the quadratic elements sorted? */

   /* nonquadratic part */
   SCIP_EXPRTREE*        exprtree;           /**< expression tree representing nonquadratic part */

   /* miscellaneous */
   char*                 name;               /**< name */
   int                   nuses;              /**< number of times, this row is referenced */
   SCIP_Real             activity;           /**< row activity value in NLP, or SCIP_INVALID if not yet calculated */
   SCIP_Longint          validactivitynlp;   /**< NLP number for which activity value is valid */
   SCIP_Real             pseudoactivity;     /**< row activity value in pseudo solution, or SCIP_INVALID if not yet calculated */
   SCIP_Longint          validpsactivitydomchg; /**< domain change number for which pseudo activity value is valid */
   SCIP_Real             minactivity;        /**< minimal activity value w.r.t. the variables' bounds, or SCIP_INVALID */
   SCIP_Real             maxactivity;        /**< maximal activity value w.r.t. the variables' bounds, or SCIP_INVALID */
   SCIP_Longint          validactivitybdsdomchg; /**< domain change number for which activity bound values are valid */
   int                   nlpindex;           /**< index of this row in NLP, or -1 if not added */
   int                   nlpiindex;          /**< index of this row in NLPI problem, or -1 if not in there */
   SCIP_Real             dualsol;            /**< dual value associated with row in last NLP solve */
   SCIP_EXPRCURV         curvature;          /**< curvature of the nonlinear row */
};

/** current NLP data */
struct SCIP_Nlp
{
   /* NLP solver */
   SCIP_NLPI*            solver;             /**< interface to NLP solver, or NULL if no NLP solvers are available */
   SCIP_NLPIPROBLEM*     problem;            /**< problem in NLP solver */

   /* status */
   int                   nunflushedvaradd;   /**< number of variable additions not flushed to NLPI problem yet */
   int                   nunflushedvardel;   /**< number of variable deletions not flushed to NLPI problem yet */
   int                   nunflushednlrowadd; /**< number of nonlinear row additions not flushed to NLPI problem yet */
   int                   nunflushednlrowdel; /**< number of nonlinear row deletions not flushed to NLPI problem yet */
   SCIP_Bool             isrelax;            /**< is the current NLP a relaxation of a SCIP problem? */
   SCIP_Bool             indiving;           /**< are we currently in diving mode? */

   /* variables in problem */
   int                   nvars;              /**< number of variables */
   int                   sizevars;           /**< allocated space for variables */
   SCIP_VAR**            vars;               /**< variables */
   SCIP_HASHMAP*         varhash;            /**< variable hash: map SCIP_VAR* to index of variable in NLP */
   /* variables in NLPI problem */
   int                   nvars_solver;       /**< number of variables in NLPI problem */
   int                   sizevars_solver;    /**< allocated space for variables in NLPI problem */
   int*                  varmap_nlp2nlpi;    /**< index of variables in NLPI problem, or -1 if variable has not been added to NLPI problem yet */
   int*                  varmap_nlpi2nlp;    /**< index of a NLPI problem variable in NLP (varmap_nlp2nlpi[varmap_nlpi2nlp[i]] == i for i = 0..nvarssolver-1), or -1 if variable has been deleted from NLP */

   /* nonlinear rows in problem */
   int                   nnlrows;            /**< number of nonlinear rows */
   int                   sizenlrows;         /**< allocated space for nonlinear rows */
   SCIP_NLROW**          nlrows;             /**< nonlinear rows */
   /* nonlinear rows in NLPI problem */
   int                   nnlrows_solver;     /**< number of nonlinear rows in solver */
   int                   sizenlrows_solver;  /**< allocated space for nonlinear rows in solver */
   int*                  nlrowmap_nlpi2nlp;  /**< index of a NLPI row in NLP (nlrows[nlrowmap_nlpi2nlp[i]]->nlpiidx == i for i = 0..nnlrows_solver-1), or -1 if row has been deleted from NLP */

   /* objective function */
   SCIP_Bool             objflushed;         /**< is the objective in the NLPI up to date? */
   SCIP_NLROW*           divingobj;          /**< objective function during diving */

   /* initial guess */
   SCIP_Bool             haveinitguess;      /**< is an initial guess available? */
   SCIP_Real*            initialguess;       /**< initial guess of primal values to use in next NLP solve, if available */

   /* solution of NLP */
   SCIP_Real             primalsolobjval;    /**< objective function value of primal solution */
   SCIP_NLPSOLSTAT       solstat;            /**< status of NLP solution (feasible, optimal, unknown...) */
   SCIP_NLPTERMSTAT      termstat;           /**< termination status of NLP (normal, some limit reached, ...) */
   SCIP_Real*            varlbdualvals;      /**< dual values associated with variable lower bounds */
   SCIP_Real*            varubdualvals;      /**< dual values associated with variable upper bounds */

   /* event handling */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   int                   globalfilterpos;    /**< position of event handler in event handler filter */

   /* fractional variables in last NLP solution */
   SCIP_VAR**            fracvars;           /**< fractional variables */
   SCIP_Real*            fracvarssol;        /**< values of the fractional variables */
   SCIP_Real*            fracvarsfrac;       /**< fractionality of the fractional variables  */
   int                   nfracvars;          /**< number of fractional variables */
   int                   npriofracvars;      /**< number of fractional variables with highest branching priority */
   int                   fracvarssize;       /**< size of fracvars* arrays */
   SCIP_Longint          validfracvars;      /**< the NLP solve for which the fractional variables are valid, or -1 if never setup */

   /* miscellaneous */
   char*                 name;               /**< problem name */
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_NLP_H__ */
